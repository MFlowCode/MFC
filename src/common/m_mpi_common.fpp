!>
!! @file
!! @brief Contains module m_mpi_common

#:include 'case.fpp'
#:include 'macros.fpp'

#if defined(MFC_SIMULATION) && defined(MFC_MPI)
#define MFC_IGR_BATCHED_EXCHANGE
#endif
#if defined(MFC_SIMULATION) && defined(MFC_MPI)
#if defined(MFC_HALO_SINGLE_ROUND) || (defined(MFC_NVSHMEM) && defined(__NVCOMPILER))
#define MFC_HALO_DEFERRED_COMPLETION
#endif
#endif

!> @brief MPI communication layer: domain decomposition, halo exchange, reductions, and parallel I/O setup
module m_mpi_common

#ifdef MFC_MPI
    use mpi  !< Message passing interface (MPI) module
#endif
#if defined(MFC_NVSHMEM) && defined(__NVCOMPILER)
    use, intrinsic :: iso_c_binding, only: c_double, c_float, c_f_pointer, c_int, c_long_long, c_ptr, c_size_t, c_sizeof
#endif

    use m_derived_types
    use m_global_parameters
    use m_helper
    use ieee_arithmetic
    use m_nvtx
    use m_constants, only: recon_type_weno
#if defined(MFC_NVSHMEM) && defined(__NVCOMPILER)
    use m_nvshmem_interfaces
#endif

    implicit none

    private :: s_apply_decomposition_policies

    integer, private :: v_size
    $:GPU_DECLARE(create='[v_size]')

    real(wp), private, allocatable, dimension(:) :: buff_send  !< Primitive variable send buffer for halo exchange
    !> Primitive variable receive buffer for halo exchange Variables for EL bubbles communication
    real(wp), private, allocatable, dimension(:) :: buff_recv
#ifdef MFC_IGR_BATCHED_EXCHANGE
    real(wp), private, allocatable, dimension(:) :: igr_buff_send, igr_buff_recv
    $:GPU_DECLARE(create='[igr_buff_send, igr_buff_recv]')
    integer(kind=8), private :: igr_face_off(0:5), igr_face_cnt(0:5)
    $:GPU_DECLARE(create='[igr_face_off]')
#endif
    integer, parameter :: halo_kind = wp
#define MFC_HALO_SEND buff_send
#define MFC_HALO_RECV buff_recv
#define MFC_MPI_HALO_P mpi_p
#if defined(MFC_NVSHMEM) && defined(__NVCOMPILER)
#define MFC_ABSOLUTE_HALO_VALUE(R) real(nvshmem_inbox(nvshmem_recv_offset + (R)), kind=stp)
#else
#define MFC_ABSOLUTE_HALO_VALUE(R) real(MFC_HALO_RECV(R), kind=stp)
#endif
#ifdef MFC_HALO_SINGLE_ROUND
#if defined(MFC_NVSHMEM) && defined(__NVCOMPILER) && defined(MFC_SIMULATION)
#define MFC_HALO_26_RECV(R) nvshmem_26_inbox(nvshmem_26_recv_base + (R))
#else
#define MFC_HALO_26_RECV(R) MFC_HALO_RECV(R)
#endif
#endif
#if defined(MFC_NVSHMEM) && defined(__NVCOMPILER) && defined(MFC_SIMULATION) && defined(MFC_MPI)
#define MFC_IGR_RECV(R) nvshmem_igr_inbox(nvshmem_igr_recv_base + (R))
#else
#define MFC_IGR_RECV(R) igr_buff_recv(R)
#endif
    type(int_bounds_info) :: comm_coords(3)
    integer               :: comm_size(3)
    !> q_beta indices to communicate: 1=void fraction, 2=d(beta)/dt, 5=energy source
    integer :: beta_vars(1:3) = [1, 2, 5]
    $:GPU_DECLARE(create='[comm_coords, comm_size, beta_vars]')

#ifdef MFC_HALO_SINGLE_ROUND
    integer, parameter :: halo_26_max_neighbors = 26
    integer            :: halo_26_beg(3,-1:1), halo_26_extent(3,-1:1)
    integer            :: halo_26_dir(3, halo_26_max_neighbors)
    integer            :: halo_26_point_offset(halo_26_max_neighbors), halo_26_point_count(halo_26_max_neighbors)
    integer            :: halo_26_total_points, halo_26_neighbor_count
#if defined(MFC_NVSHMEM) && defined(__NVCOMPILER)
    !> Active-box trimmed exchange (NVSHMEM only): per-neighbor sub-box, in the region-local coordinates of the shared neighbor
    !! enumeration, of the points the receiving rank's window plus stencil margin actually consumes. Side 1 enumerates what this
    !! rank packs for each partner (the partner's consumption), side 2 what this rank unpacks (its own consumption). Both sides
    !! default to the full region, so an untrimmed run enumerates the full layout; only s_mpi_apply_26_trim ever shrinks them.
    integer :: halo_26_trim_lo(3, 2, halo_26_max_neighbors)
    integer :: halo_26_trim_hi(3, 2, halo_26_max_neighbors)
    integer :: halo_26_trim_cnt(2, halo_26_max_neighbors)
    integer :: halo_26_trim_off(2,0:halo_26_max_neighbors)
    integer :: halo_26_pack_points, halo_26_recv_points
    logical :: halo_26_trim_active = .false.
    logical :: halo_26_trim_requested = .false.
    $:GPU_DECLARE(create='[halo_26_trim_lo, halo_26_trim_hi, halo_26_trim_cnt, halo_26_trim_off, halo_26_pack_points, &
                  & halo_26_recv_points, halo_26_trim_active]')
#define MFC_26_PACK_POINTS halo_26_pack_points
#define MFC_26_RECV_POINTS halo_26_recv_points
#else
#define MFC_26_PACK_POINTS halo_26_total_points
#define MFC_26_RECV_POINTS halo_26_total_points
#endif
#ifdef MFC_HALO_INDEX_CACHE
    integer, allocatable :: halo_26_indices(:,:)
    $:GPU_DECLARE(create='[halo_26_indices]')
#define MFC_26_POINT_INDICES s_get_26_cached_indices
#else
#define MFC_26_POINT_INDICES s_get_26_neighbor_indices
#endif
    $:GPU_DECLARE(create='[halo_26_beg, halo_26_extent, halo_26_dir, halo_26_point_offset, halo_26_total_points, halo_26_neighbor_count]')

    integer :: halo_26_send_req(halo_26_max_neighbors), halo_26_recv_req(halo_26_max_neighbors)
#if defined(MFC_MPI)
    ! Persistent-request state for the rdma_mpi GPU-direct arm. Requests are bound per
    ! distinct v_size (one set per exchange width), started at each exchange's post phase,
    ! and completed lazily at the next exchange's post phase (recvs are waited explicitly
    ! in the wait phase first; only the sends are left outstanding).
    integer, parameter :: halo_26_rdma_max_vsizes = 4
    integer            :: halo_26_rdma_vsize(halo_26_rdma_max_vsizes) = 0
    integer            :: halo_26_rdma_nvsizes = 0
    integer            :: halo_26_rdma_recv(halo_26_max_neighbors*halo_26_rdma_max_vsizes) = MPI_REQUEST_NULL
    integer            :: halo_26_rdma_send(halo_26_max_neighbors*halo_26_rdma_max_vsizes) = MPI_REQUEST_NULL
    integer            :: halo_26_rdma_pending_slot = 0
    logical            :: halo_26_rdma_pending = .false.
#endif
#endif

    !> Weighted partition, enabled when MFC_SFC_WEIGHTS names a readable file. Per-block cost weights drive a weight-balanced
    !! Cartesian decomposition computed once at decomposition time.
    integer, private               :: sfc_gcells(1:3) = 0  !< global cell counts per direction
    integer, private               :: sfc_nblk(1:3) = 0    !< weight blocks per direction
    integer, private               :: sfc_blksz = 1        !< block edge length in cells
    integer, private               :: sfc_snap(1:3) = 1    !< cut snapping stride in cells (1 = none)
    real(wp), private, allocatable :: sfc_wblk(:)          !< per-block weight sums, x fastest
    logical, private               :: sfc_active = .false.

#ifndef __NVCOMPILER_GPU_UNIFIED_MEM
    $:GPU_DECLARE(create='[buff_send, buff_recv]')
#endif

    integer(kind=8) :: halo_size
    $:GPU_DECLARE(create='[halo_size]')

    logical, private :: exchange_all_chemistry_temperatures = .false.
    logical, private :: use_rdma_transport = .false.

#if defined(MFC_NVSHMEM) && defined(__NVCOMPILER)
#ifdef MFC_SINGLE_PRECISION
    real(c_float), pointer, contiguous, private :: nvshmem_inbox(:)
#else
    real(c_double), pointer, contiguous, private :: nvshmem_inbox(:)
#endif
#ifdef MFC_HALO_SINGLE_ROUND
#ifdef MFC_SINGLE_PRECISION
    real(c_float), pointer, contiguous, private :: nvshmem_26_inbox(:)
#else
    real(c_double), pointer, contiguous, private :: nvshmem_26_inbox(:)
#endif
    integer(kind=8), private :: nvshmem_26_recv_base
    $:GPU_DECLARE(create='[nvshmem_26_recv_base]')
    !> Bootstrap 26-neighbor transport configuration, retained so the active-box trim can re-run mfc_nvshmem_configure_26 with
    !! trimmed point counts on the same destinations, channels, and full-size region offsets.
    integer, private              :: nvshmem_26_partner_offsets(halo_26_max_neighbors)
    integer(c_int), private       :: nvshmem_26_dst_pe(halo_26_max_neighbors), nvshmem_26_channel(halo_26_max_neighbors)
    integer(c_long_long), private :: nvshmem_26_src_points(halo_26_max_neighbors), nvshmem_26_dst_points(halo_26_max_neighbors)
    !> Device-resident epoch exchange (MFC_HALO_EPOCH_DEVICE=1). The 26-neighbor exchange's epoch, parity, and receive base live in
    !! device memory, every exchange launch is parameter-free, and the puts use the k_put_26_dev block copy, so the exchange is
    !! CUDA-graph capturable without per-exchange host updates or put retirement. The mode targets node-local P2P fabrics, where the
    !! block put completes with its kernel and stream order alone makes repacking the source buffer safe; multinode proxy puts use
    !! the retire path. Public so m_rhs can extend its CUDA-graph capture over the communication.
    logical :: halo_26_epoch_device = .false.
    !> Set once the first exchange has bound nvshmem_26_recv_base's device address into the backend. The bind happens lazily at the
    !! first exchange, inside the main data region and before any CUDA-graph capture, so a captured region never records the
    !! cudaMemcpyToSymbol.
    logical, private :: halo_26_epoch_bound = .false.
#endif
#ifdef MFC_SIMULATION
#ifdef MFC_SINGLE_PRECISION
    real(c_float), pointer, contiguous, private :: nvshmem_igr_inbox(:)
#else
    real(c_double), pointer, contiguous, private :: nvshmem_igr_inbox(:)
#endif
    integer(kind=8), private :: nvshmem_igr_recv_base
    $:GPU_DECLARE(create='[nvshmem_igr_recv_base]')
#endif
#endif

contains

    !> Initialize the module.
    impure subroutine s_initialize_mpi_common_module(exchange_all_chemistry_temperatures_in, use_rdma_transport_in)

        logical, intent(in) :: exchange_all_chemistry_temperatures_in
        logical, intent(in) :: use_rdma_transport_in

#ifdef MFC_HALO_SINGLE_ROUND
        integer :: axis, dx, dy, dz, point_count
#ifdef MFC_HALO_INDEX_CACHE
        integer :: cache_point, cache_j, cache_k, cache_l, cache_sj, cache_sk, cache_sl
#endif
        integer :: grid_dims(3)
#endif
#if defined(MFC_NVSHMEM) && defined(__NVCOMPILER)
        integer(kind=8)   :: nvshmem_local_halo_elements
        integer(kind=8)   :: nvshmem_halo_elements
        integer(c_size_t) :: nvshmem_inbox_elements
        integer(c_int)    :: nvshmem_status
        type(c_ptr)       :: nvshmem_inbox_cptr
#ifdef MFC_HALO_SINGLE_ROUND
        integer(c_size_t)    :: nvshmem_26_inbox_elements
        type(c_ptr)          :: nvshmem_26_inbox_cptr
        integer              :: neighbor_idx
        integer              :: nvshmem_26_reqs(2*halo_26_max_neighbors), nvshmem_26_nreq
        integer(c_long_long) :: nvshmem_26_point_count(halo_26_max_neighbors)
        character(len=8)     :: epoch_env
        integer              :: epoch_env_status
#endif
        integer :: ierr
#endif

        exchange_all_chemistry_temperatures = exchange_all_chemistry_temperatures_in
        use_rdma_transport = use_rdma_transport_in

#ifdef MFC_MPI
        ! Allocating buff_send/recv and. Please note that for the sake of simplicity, both variables are provided sufficient storage
        ! to hold the largest buffer in the computational domain.

        if (qbmm .and. .not. polytropic) then
            v_size = sys_size + 2*nb*nnode
        else if (heat_conduction .or. (chemistry .and. (chem_params%diffusion .or. exchange_all_chemistry_temperatures))) then
            v_size = sys_size + 1
        else
            v_size = sys_size
        end if

        if (n > 0) then
            if (p > 0) then
                halo_size = nint(-1._wp + 1._wp*buff_size*(v_size)*(m + 2*buff_size + 1)*(n + 2*buff_size + 1)*(p + 2*buff_size &
                                 & + 1)/(cells_bounds%mnp_min + 2*buff_size + 1))
            else
                halo_size = -1 + buff_size*(v_size)*(cells_bounds%mn_max + 2*buff_size + 1)
            end if
        else
            halo_size = -1 + buff_size*(v_size)
        end if

#ifdef MFC_HALO_SINGLE_ROUND
        grid_dims = (/m, n, p/)
        halo_26_beg = 0
        halo_26_extent = 0
        halo_26_dir = 0
        halo_26_point_offset = 0
        halo_26_point_count = 0

        do axis = 1, 3
            halo_26_beg(axis, 0) = 0
            halo_26_extent(axis, 0) = grid_dims(axis) + 1
            if (axis <= num_dims) then
                if (nidx(axis)%beg == -1) then
                    halo_26_beg(axis, -1) = -buff_size
                    halo_26_extent(axis, -1) = buff_size
                else
                    halo_26_beg(axis, 0) = -buff_size
                    halo_26_extent(axis, 0) = halo_26_extent(axis, 0) + buff_size
                end if
                if (nidx(axis)%end == 1) then
                    halo_26_beg(axis, 1) = grid_dims(axis) + 1
                    halo_26_extent(axis, 1) = buff_size
                else
                    halo_26_extent(axis, 0) = halo_26_extent(axis, 0) + buff_size
                end if
            end if
        end do

        halo_26_total_points = 0
        halo_26_neighbor_count = 0
        do dz = nidx(3)%beg, nidx(3)%end
            do dy = nidx(2)%beg, nidx(2)%end
                do dx = nidx(1)%beg, nidx(1)%end
                    if (abs(dx) + abs(dy) + abs(dz) == 0) cycle
                    point_count = halo_26_extent(1, dx)*halo_26_extent(2, dy)*halo_26_extent(3, dz)
                    halo_26_neighbor_count = halo_26_neighbor_count + 1
                    halo_26_dir(:,halo_26_neighbor_count) = (/dx, dy, dz/)
                    halo_26_point_offset(halo_26_neighbor_count) = halo_26_total_points
                    halo_26_point_count(halo_26_neighbor_count) = point_count
                    halo_26_total_points = halo_26_total_points + point_count
                end do
            end do
        end do
        halo_size = max(halo_size, int(v_size, kind=8)*int(halo_26_total_points, kind=8) - 1_8)

#ifdef MFC_HALO_INDEX_CACHE
        @:ALLOCATE(halo_26_indices(0:max(0, halo_26_total_points - 1), 1:6))
        halo_26_indices = 0
        do cache_point = 0, halo_26_total_points - 1
            call s_get_26_neighbor_indices(cache_point, cache_j, cache_k, cache_l, cache_sj, cache_sk, cache_sl)
            halo_26_indices(cache_point,:) = [cache_j, cache_k, cache_l, cache_sj, cache_sk, cache_sl]
        end do
        $:GPU_UPDATE(device='[halo_26_indices]')
#endif
#endif

        $:GPU_UPDATE(device='[halo_size, v_size]')
#ifdef MFC_HALO_SINGLE_ROUND
        $:GPU_UPDATE(device='[halo_26_beg, halo_26_extent, halo_26_dir, halo_26_point_offset, halo_26_total_points, halo_26_neighbor_count]')
#if defined(MFC_NVSHMEM) && defined(__NVCOMPILER)
        ! Trimmed-exchange state starts as an exact copy of the full layout. MFC_HALO_26_TRIM=1
        ! opts the run into shrinking it from the active-box window (s_mpi_apply_26_trim).
        block
            character(len=8) :: trim_env
            integer          :: trim_env_status, trim_axis, trim_dir, trim_n, trim_cnt
            halo_26_trim_active = .false.
            halo_26_trim_requested = .false.
            call get_environment_variable('MFC_HALO_26_TRIM', value=trim_env, status=trim_env_status)
            if (trim_env_status == 0 .and. trim(trim_env) == '1') halo_26_trim_requested = .true.
            do trim_n = 1, halo_26_neighbor_count
                trim_cnt = 1
                do trim_axis = 1, 3
                    trim_dir = halo_26_dir(trim_axis, trim_n)
                    halo_26_trim_lo(trim_axis, 1, trim_n) = 0
                    halo_26_trim_hi(trim_axis, 1, trim_n) = halo_26_extent(trim_axis, trim_dir) - 1
                    halo_26_trim_lo(trim_axis, 2, trim_n) = 0
                    halo_26_trim_hi(trim_axis, 2, trim_n) = halo_26_extent(trim_axis, trim_dir) - 1
                    trim_cnt = trim_cnt*halo_26_extent(trim_axis, trim_dir)
                end do
                halo_26_trim_cnt(1, trim_n) = trim_cnt
                halo_26_trim_cnt(2, trim_n) = trim_cnt
            end do
            do trim_n = 1, halo_26_neighbor_count
                halo_26_trim_off(1, trim_n) = halo_26_point_offset(trim_n)
                halo_26_trim_off(2, trim_n) = halo_26_point_offset(trim_n)
            end do
            halo_26_trim_off(1, 0) = 0
            halo_26_trim_off(2, 0) = 0
            halo_26_pack_points = halo_26_total_points
            halo_26_recv_points = halo_26_total_points
            $:GPU_UPDATE(device='[halo_26_trim_lo, halo_26_trim_hi, halo_26_trim_cnt, halo_26_trim_off, halo_26_pack_points, &
                         & halo_26_recv_points, halo_26_trim_active]')
        end block
#endif
#endif

#ifdef MFC_IGR_BATCHED_EXCHANGE
        block
            integer(kind=8) :: face_elems(0:5), fi
            face_elems(0:1) = int(buff_size, 8)*int(n + 1, 8)*int(p + 1, 8)
            face_elems(2:3) = int(buff_size, 8)*int(m + 1, 8)*int(p + 1, 8)
            face_elems(4:5) = int(buff_size, 8)*int(m + 1, 8)*int(n + 1, 8)
            igr_face_off(0) = 0_8
            do fi = 0, 5
                igr_face_cnt(fi) = face_elems(fi)
                if (fi > 0) igr_face_off(fi) = igr_face_off(fi - 1) + face_elems(fi - 1)
            end do
            @:ALLOCATE(igr_buff_send(0:igr_face_off(5) + face_elems(5) - 1), igr_buff_recv(0:igr_face_off(5) + face_elems(5) - 1))
        end block
        $:GPU_UPDATE(device='[igr_face_off]')
#endif

#ifndef __NVCOMPILER_GPU_UNIFIED_MEM
        @:ALLOCATE(buff_send(0:halo_size), buff_recv(0:halo_size))
#else
        allocate (buff_send(0:halo_size), buff_recv(0:halo_size))
        $:GPU_ENTER_DATA(create='[capture:buff_send]')
        $:GPU_ENTER_DATA(create='[capture:buff_recv]')
#endif

#if defined(MFC_NVSHMEM) && defined(__NVCOMPILER)
        nvshmem_local_halo_elements = halo_size + 1
        call MPI_ALLREDUCE(nvshmem_local_halo_elements, nvshmem_halo_elements, 1, MPI_INTEGER8, MPI_MAX, MPI_COMM_WORLD, ierr)
#ifdef MFC_SINGLE_PRECISION
        if (mfc_nvshmem_halo_value_bytes() /= int(c_sizeof(1.0_c_float), &
            & c_int)) call s_mpi_abort('The NVSHMEM backend was compiled with a different halo precision than the solver.')
#else
        if (mfc_nvshmem_halo_value_bytes() /= int(c_sizeof(1.0_c_double), &
            & c_int)) call s_mpi_abort('The NVSHMEM backend was compiled with a different halo precision than the solver.')
#endif
        nvshmem_status = mfc_nvshmem_initialize(int(nvshmem_halo_elements, c_size_t), int(MPI_COMM_WORLD, c_int), &
                                                & nvshmem_inbox_cptr, nvshmem_inbox_elements)
        if (nvshmem_status /= 0) call s_mpi_abort('Unable to initialize the NVSHMEM halo backend.')
        call c_f_pointer(nvshmem_inbox_cptr, nvshmem_inbox, [nvshmem_inbox_elements])
        nvshmem_inbox(0:) => nvshmem_inbox
#ifdef MFC_HALO_SINGLE_ROUND
        nvshmem_status = mfc_nvshmem_initialize_26(int(nvshmem_halo_elements, c_size_t), nvshmem_26_inbox_cptr, &
            & nvshmem_26_inbox_elements)
        if (nvshmem_status /= 0) call s_mpi_abort('Unable to initialize the NVSHMEM 26-neighbor halo backend.')
        call c_f_pointer(nvshmem_26_inbox_cptr, nvshmem_26_inbox, [nvshmem_26_inbox_elements])
        nvshmem_26_inbox(0:) => nvshmem_26_inbox
        ! MFC_COMM_BEGIN name=nvshmem_bootstrap mode=mpi_only pattern=setup
        nvshmem_26_nreq = 0
        do neighbor_idx = 1, halo_26_neighbor_count
            dx = halo_26_dir(1, neighbor_idx)
            dy = halo_26_dir(2, neighbor_idx)
            dz = halo_26_dir(3, neighbor_idx)
            call MPI_IRECV(nvshmem_26_partner_offsets(neighbor_idx), 1, MPI_INTEGER, neighbor_ranks(dx, dy, dz), &
                           & f_26_neighbor_tag(dx, dy, dz), MPI_COMM_WORLD, nvshmem_26_reqs(nvshmem_26_nreq + 1), ierr)
            call MPI_ISEND(halo_26_point_offset(neighbor_idx), 1, MPI_INTEGER, neighbor_ranks(dx, dy, dz), f_26_neighbor_tag(-dx, &
                           & -dy, -dz), MPI_COMM_WORLD, nvshmem_26_reqs(nvshmem_26_nreq + 2), ierr)
            nvshmem_26_nreq = nvshmem_26_nreq + 2
        end do
        if (nvshmem_26_nreq > 0) then
            call MPI_WAITALL(nvshmem_26_nreq, nvshmem_26_reqs(1:nvshmem_26_nreq), MPI_STATUSES_IGNORE, ierr)
        end if

        do neighbor_idx = 1, halo_26_neighbor_count
            dx = halo_26_dir(1, neighbor_idx)
            dy = halo_26_dir(2, neighbor_idx)
            dz = halo_26_dir(3, neighbor_idx)
            nvshmem_26_dst_pe(neighbor_idx) = int(neighbor_ranks(dx, dy, dz), c_int)
            #:if defined('MFC_SIMULATION')
#include "nvshmem_halo_26_channel.inc"
            #:endif
            nvshmem_26_src_points(neighbor_idx) = int(halo_26_point_offset(neighbor_idx), c_long_long)
            nvshmem_26_dst_points(neighbor_idx) = int(nvshmem_26_partner_offsets(neighbor_idx), c_long_long)
            nvshmem_26_point_count(neighbor_idx) = int(halo_26_point_count(neighbor_idx), c_long_long)
        end do
        nvshmem_status = mfc_nvshmem_configure_26(int(halo_26_neighbor_count, c_int), nvshmem_26_dst_pe, nvshmem_26_src_points, &
            & nvshmem_26_dst_points, nvshmem_26_point_count, nvshmem_26_channel)
        if (nvshmem_status /= 0) call s_mpi_abort('Unable to configure the NVSHMEM 26-neighbor topology.')
        call get_environment_variable('MFC_HALO_EPOCH_DEVICE', epoch_env, status=epoch_env_status)
        halo_26_epoch_device = epoch_env_status == 0 .and. trim(epoch_env) == '1'
#endif
#ifdef MFC_SIMULATION
        block
            integer(kind=8)       :: igr_local_elements, igr_max_elements, partner_off(0:5)
            integer               :: fi, fdir, floc, fpartner, igr_nreq, n_act
            integer               :: igr_reqs(0:11)
            integer(c_size_t)     :: nvshmem_igr_inbox_elements
            type(c_ptr)           :: nvshmem_igr_inbox_cptr
            integer(c_int)        :: igr_dst_pe(0:5), igr_face_id(0:5)
            integer(c_long_long)  :: igr_src_off(0:5), igr_dst_off(0:5), igr_cnt(0:5)
            type(int_bounds_info) :: fb(1:3)

            igr_local_elements = igr_face_off(5) + igr_face_cnt(5)
            call MPI_ALLREDUCE(igr_local_elements, igr_max_elements, 1, MPI_INTEGER8, MPI_MAX, MPI_COMM_WORLD, ierr)
            nvshmem_status = mfc_nvshmem_initialize_igr(int(igr_max_elements, c_size_t), nvshmem_igr_inbox_cptr, &
                & nvshmem_igr_inbox_elements)
            if (nvshmem_status /= 0) call s_mpi_abort('Unable to initialize the NVSHMEM IGR backend.')
            call c_f_pointer(nvshmem_igr_inbox_cptr, nvshmem_igr_inbox, [nvshmem_igr_inbox_elements])
            nvshmem_igr_inbox(0:) => nvshmem_igr_inbox

            fb(1) = bc_x; fb(2) = bc_y; fb(3) = bc_z
            partner_off = -1_8
            igr_nreq = 0
            do fdir = 1, num_dims
                do floc = -1, 1, 2
                    fi = 2*(fdir - 1) + f_logical_to_int(floc == 1)
                    fpartner = fb(fdir)%beg
                    if (floc == 1) fpartner = fb(fdir)%end
                    if (fpartner < 0) cycle
                    call MPI_IRECV(partner_off(fi), 1, MPI_INTEGER8, fpartner, 6100 + ieor(fi, 1), MPI_COMM_WORLD, &
                                   & igr_reqs(igr_nreq), ierr)
                    igr_nreq = igr_nreq + 1
                    call MPI_ISEND(igr_face_off(fi), 1, MPI_INTEGER8, fpartner, 6100 + fi, MPI_COMM_WORLD, igr_reqs(igr_nreq), ierr)
                    igr_nreq = igr_nreq + 1
                end do
            end do
            if (igr_nreq > 0) call MPI_WAITALL(igr_nreq, igr_reqs(0:igr_nreq - 1), MPI_STATUSES_IGNORE, ierr)

            n_act = 0
            do fdir = 1, num_dims
                do floc = -1, 1, 2
                    fi = 2*(fdir - 1) + f_logical_to_int(floc == 1)
                    fpartner = fb(fdir)%beg
                    if (floc == 1) fpartner = fb(fdir)%end
                    if (fpartner < 0) cycle
                    igr_dst_pe(n_act) = int(fpartner, c_int)
                    igr_face_id(n_act) = int(fi, c_int)
                    igr_src_off(n_act) = int(igr_face_off(fi), c_long_long)
                    igr_dst_off(n_act) = int(partner_off(fi), c_long_long)
                    igr_cnt(n_act) = int(igr_face_cnt(fi), c_long_long)
                    n_act = n_act + 1
                end do
            end do
            nvshmem_status = mfc_nvshmem_configure_igr(int(n_act, c_int), igr_dst_pe, igr_face_id, igr_src_off, igr_dst_off, &
                & igr_cnt)
            if (nvshmem_status /= 0) call s_mpi_abort('Unable to configure the NVSHMEM IGR topology.')
        end block
        ! MFC_COMM_END
#endif
#endif
#endif

        $:GPU_UPDATE(device='[beta_vars]')

    end subroutine s_initialize_mpi_common_module

    !> Initialize the MPI execution environment and query the number of processors and local rank.
    impure subroutine s_mpi_initialize

#ifdef MFC_MPI
        integer :: ierr  !< Generic flag used to identify and report MPI errors

        call MPI_INIT(ierr)

        if (ierr /= MPI_SUCCESS) then
            print '(A)', 'Unable to initialize MPI environment. Exiting.'
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        end if

        call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)

        call MPI_COMM_RANK(MPI_COMM_WORLD, proc_rank, ierr)
#else
        num_procs = 1
        proc_rank = 0
#endif

        $:GPU_UPDATE(device='[num_procs, proc_rank]')

    end subroutine s_mpi_initialize

    !> Set up MPI I/O data views and variable pointers for parallel file output.
    impure subroutine s_initialize_mpi_data(q_cons_vf, ib_markers, ib_mpi_data, beta, qbmm_pb, qbmm_mv)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        type(integer_field), optional, intent(in)           :: ib_markers
        type(mpi_io_ib_var), optional, intent(inout)        :: ib_mpi_data
        type(scalar_field), intent(in), optional            :: beta
        type(pres_field), intent(in), optional              :: qbmm_pb, qbmm_mv
        integer, dimension(num_dims)                        :: sizes_glb, sizes_loc

#ifdef MFC_MPI
        integer :: i, j
        integer :: ierr  !< Generic flag used to identify and report MPI errors
        integer :: alt_sys
        logical :: bind_qbmm_fields

        if (present(qbmm_pb) .neqv. present(qbmm_mv)) then
            call s_mpi_abort('QBMM MPI I/O requires both pressure and moment fields.')
        end if
        bind_qbmm_fields = qbmm .and. .not. polytropic .and. present(qbmm_pb) .and. present(qbmm_mv)

        if (present(beta)) then
            alt_sys = sys_size + 1
        else
            alt_sys = sys_size
        end if

        do i = 1, sys_size
            MPI_IO_DATA%var(i)%sf => q_cons_vf(i)%sf(0:m,0:n,0:p)
        end do

        if (present(beta)) then
            MPI_IO_DATA%var(alt_sys)%sf => beta%sf(0:m,0:n,0:p)
        end if

        ! Additional variables pb and mv for non-polytropic qbmm
        if (bind_qbmm_fields) then
            do i = 1, nb
                do j = 1, nnode
                    MPI_IO_DATA%var(sys_size + (i - 1)*nnode + j)%sf => qbmm_pb%sf(0:m,0:n,0:p,j, i)
                    MPI_IO_DATA%var(sys_size + (i - 1)*nnode + j + nb*nnode)%sf => qbmm_mv%sf(0:m,0:n,0:p,j, i)
                end do
            end do
        end if

        ! Define global(g) and local(l) sizes for flow variables
        sizes_glb(1) = m_glb + 1; sizes_loc(1) = m + 1
        if (n > 0) then
            sizes_glb(2) = n_glb + 1; sizes_loc(2) = n + 1
            if (p > 0) then
                sizes_glb(num_dims) = p_glb + 1; sizes_loc(num_dims) = p + 1
            end if
        end if

        ! Define the view for each variable
        do i = 1, alt_sys
            call MPI_TYPE_CREATE_SUBARRAY(num_dims, sizes_glb, sizes_loc, start_idx, MPI_ORDER_FORTRAN, mpi_p, &
                                          & MPI_IO_DATA%view(i), ierr)
            call MPI_TYPE_COMMIT(MPI_IO_DATA%view(i), ierr)
        end do

        if (bind_qbmm_fields) then
            do i = sys_size + 1, sys_size + 2*nb*nnode
                call MPI_TYPE_CREATE_SUBARRAY(num_dims, sizes_glb, sizes_loc, start_idx, MPI_ORDER_FORTRAN, mpi_p, &
                                              & MPI_IO_DATA%view(i), ierr)
                call MPI_TYPE_COMMIT(MPI_IO_DATA%view(i), ierr)
            end do
        end if

        if (present(ib_markers) .neqv. present(ib_mpi_data)) then
            call s_mpi_abort('Immersed-boundary MPI I/O requires both marker and descriptor fields.')
        end if

        if (present(ib_markers)) then
            ib_mpi_data%var%sf => ib_markers%sf(0:m,0:n,0:p)
            call MPI_TYPE_CREATE_SUBARRAY(num_dims, sizes_glb, sizes_loc, start_idx, MPI_ORDER_FORTRAN, MPI_INTEGER, &
                                          & ib_mpi_data%view, ierr)
            call MPI_TYPE_COMMIT(ib_mpi_data%view, ierr)
        end if
#endif

    end subroutine s_initialize_mpi_data

    !> Set up MPI I/O data views for downsampled (coarsened) parallel file output.
    subroutine s_initialize_mpi_data_ds(m_ds, n_ds, p_ds, q_cons_vf)

        integer, intent(in)                                           :: m_ds, n_ds, p_ds
        type(scalar_field), dimension(sys_size), intent(in), optional :: q_cons_vf
        integer, dimension(num_dims)                                  :: sizes_loc
        integer, dimension(3)                                         :: sf_start_idx

#ifdef MFC_MPI
        integer :: i, ierr

        sf_start_idx = (/0, 0, 0/)

        if (present(q_cons_vf)) then
            do i = 1, sys_size
                MPI_IO_DATA%var(i)%sf => q_cons_vf(i)%sf(-1:m_ds + 1,-1:n_ds + 1,-1:p_ds + 1)
            end do
        end if
        ! Define global(g) and local(l) sizes for flow variables
        sizes_loc(1) = m_ds + 3
        if (n > 0) then
            sizes_loc(2) = n_ds + 3
            if (p > 0) then
                sizes_loc(num_dims) = p_ds + 3
            end if
        end if

        ! Define the view for each variable
        do i = 1, sys_size
            call MPI_TYPE_CREATE_SUBARRAY(num_dims, sizes_loc, sizes_loc, sf_start_idx, MPI_ORDER_FORTRAN, mpi_p, &
                                          & MPI_IO_DATA%view(i), ierr)
            call MPI_TYPE_COMMIT(MPI_IO_DATA%view(i), ierr)
        end do
#endif

    end subroutine s_initialize_mpi_data_ds

    !> Gather variable-length real vectors from all MPI ranks onto the root process.
    impure subroutine s_mpi_gather_data(my_vector, counts, gathered_vector, root)

        integer, intent(in)                     :: counts              !< Array of vector lengths for each process
        real(wp), intent(in), dimension(counts) :: my_vector           !< Input vector on each process
        integer, intent(in)                     :: root                !< Rank of the root process
        real(wp), allocatable, intent(out)      :: gathered_vector(:)  !< Gathered vector on the root process
        integer                                 :: i
        integer                                 :: ierr                !< Generic flag used to identify and report MPI errors
        integer, allocatable                    :: recounts(:), displs(:)

#ifdef MFC_MPI
        allocate (recounts(num_procs))

        call MPI_GATHER(counts, 1, MPI_INTEGER, recounts, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)

        allocate (displs(size(recounts)))

        displs(1) = 0

        do i = 2, size(recounts)
            displs(i) = displs(i - 1) + recounts(i - 1)
        end do

        allocate (gathered_vector(sum(recounts)))
        call MPI_GATHERV(my_vector, counts, mpi_p, gathered_vector, recounts, displs, mpi_p, root, MPI_COMM_WORLD, ierr)
#endif

    end subroutine s_mpi_gather_data

    !> Gather per-rank time step wall-clock times onto rank 0 for performance reporting.
    impure subroutine mpi_bcast_time_step_values(proc_time, time_avg)

        real(wp), dimension(0:num_procs - 1), intent(inout) :: proc_time
        real(wp), intent(inout)                             :: time_avg

#ifdef MFC_MPI
        integer :: ierr  !< Generic flag used to identify and report MPI errors

        call MPI_GATHER(time_avg, 1, mpi_p, proc_time(0), 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
#endif

    end subroutine mpi_bcast_time_step_values

    !> Print a case file error with the prohibited condition and message, then abort execution.
    impure subroutine s_prohibit_abort(condition, message)

        character(len=*), intent(in) :: condition, message

        print *, ""
        print *, "CASE FILE ERROR"
        print *, "  - Prohibited condition: ", trim(condition)
        if (len_trim(message) > 0) then
            print *, "  - Note: ", trim(message)
        end if
        print *, ""
        call s_mpi_abort(code=CASE_FILE_ERROR_CODE)

    end subroutine s_prohibit_abort

    !> The goal of this subroutine is to determine the global extrema of the stability criteria in the computational domain. This is
    !! performed by sifting through the local extrema of each stability criterion. Note that each of the local extrema is from a
    !! single process, within its assigned section of the computational domain. Finally, note that the global extrema values are
    !! only bookkept on the rank 0 processor. `max_loc`/`max_glb` hold the max-reduced criteria (ICFL, VCFL, CCFL, TCFL) and
    !! `min_loc`/`min_glb` hold the min-reduced criteria (Rc).
    impure subroutine s_mpi_reduce_stability_criteria_extrema(max_loc, min_loc, bubs_loc, max_glb, min_glb, bubs_glb)

        real(wp), dimension(:), intent(in)  :: max_loc
        real(wp), dimension(:), intent(in)  :: min_loc
        integer, intent(in)                 :: bubs_loc
        real(wp), dimension(:), intent(out) :: max_glb
        real(wp), dimension(:), intent(out) :: min_glb
        integer, intent(out)                :: bubs_glb

        max_glb = max_loc
        min_glb = min_loc

#ifdef MFC_MPI
        block
            integer :: ierr

            bubs_glb = 0
            call MPI_REDUCE(max_loc(1), max_glb(1), 1, mpi_p, MPI_MAX, 0, MPI_COMM_WORLD, ierr)  ! ICFL

            if (viscous) then
                call MPI_REDUCE(max_loc(2), max_glb(2), 1, mpi_p, MPI_MAX, 0, MPI_COMM_WORLD, ierr)  ! VCFL
                call MPI_REDUCE(min_loc(1), min_glb(1), 1, mpi_p, MPI_MIN, 0, MPI_COMM_WORLD, ierr)  ! Rc
            end if

            if (surface_tension) then
                call MPI_REDUCE(max_loc(3), max_glb(3), 1, mpi_p, MPI_MAX, 0, MPI_COMM_WORLD, ierr)  ! CCFL
            end if

            if (heat_conduction) then
                call MPI_REDUCE(max_loc(4), max_glb(4), 1, mpi_p, MPI_MAX, 0, MPI_COMM_WORLD, ierr)  ! TCFL
            end if

            if (bubbles_lagrange) then
                call MPI_REDUCE(bubs_loc, bubs_glb, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            end if
        end block
#else
        bubs_glb = 0
        if (bubbles_lagrange) bubs_glb = bubs_loc
#endif

    end subroutine s_mpi_reduce_stability_criteria_extrema

    !> Reduce a local integer value to its global sum across all MPI ranks.
    subroutine s_mpi_reduce_int_sum(var_loc, sum)

        integer, intent(in)  :: var_loc
        integer, intent(out) :: sum

#ifdef MFC_MPI
        integer :: ierr  !< Generic flag used to identify and report MPI errors

        call MPI_REDUCE(var_loc, sum, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
#else
        sum = var_loc
#endif

    end subroutine s_mpi_reduce_int_sum

    !> Reduce a local real value to its global sum across all MPI ranks.
    impure subroutine s_mpi_allreduce_sum(var_loc, var_glb)

        real(wp), intent(in)  :: var_loc
        real(wp), intent(out) :: var_glb

#ifdef MFC_MPI
        integer :: ierr  !< Generic flag used to identify and report MPI errors

        call MPI_ALLREDUCE(var_loc, var_glb, 1, mpi_p, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

    end subroutine s_mpi_allreduce_sum

    !> Reduce an array of vectors to their global sums across all MPI ranks.
    impure subroutine s_mpi_allreduce_vectors_sum(var_loc, var_glb, num_vectors, vector_length)

        integer, intent(in)                     :: num_vectors, vector_length
        real(wp), dimension(:,:), intent(in)    :: var_loc
        real(wp), dimension(:,:), intent(inout) :: var_glb

#ifdef MFC_MPI
        integer :: ierr  !< Generic flag used to identify and report MPI errors

        if (loc(var_loc) == loc(var_glb)) then
            call MPI_Allreduce(MPI_IN_PLACE, var_glb, num_vectors*vector_length, mpi_p, MPI_SUM, MPI_COMM_WORLD, ierr)
        else
            call MPI_Allreduce(var_loc, var_glb, num_vectors*vector_length, mpi_p, MPI_SUM, MPI_COMM_WORLD, ierr)
        end if
#else
        var_glb(1:num_vectors,1:vector_length) = var_loc(1:num_vectors,1:vector_length)
#endif

    end subroutine s_mpi_allreduce_vectors_sum

    !> Reduce a local integer value to its global sum across all MPI ranks.
    impure subroutine s_mpi_allreduce_integer_sum(var_loc, var_glb)

        integer(kind=8), intent(in)  :: var_loc
        integer(kind=8), intent(out) :: var_glb

#ifdef MFC_MPI
        integer :: ierr  !< Generic flag used to identify and report MPI errors

        call MPI_ALLREDUCE(var_loc, var_glb, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
        var_glb = var_loc
#endif

    end subroutine s_mpi_allreduce_integer_sum

    !> Reduce an integer(8) array to its element-wise global sum across all MPI ranks.
    impure subroutine s_mpi_allreduce_integer_sum_vec(var_loc, var_glb)

        integer(kind=8), dimension(:), intent(in)  :: var_loc
        integer(kind=8), dimension(:), intent(out) :: var_glb

#ifdef MFC_MPI
        integer :: ierr  !< Generic flag used to identify and report MPI errors

        call MPI_ALLREDUCE(var_loc, var_glb, size(var_loc), MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
        var_glb = var_loc
#endif

    end subroutine s_mpi_allreduce_integer_sum_vec

    !> Reduce a local real value to its global minimum across all MPI ranks.
    impure subroutine s_mpi_allreduce_min(var_loc, var_glb)

        real(wp), intent(in)  :: var_loc
        real(wp), intent(out) :: var_glb

#ifdef MFC_MPI
        integer :: ierr  !< Generic flag used to identify and report MPI errors

        call MPI_ALLREDUCE(var_loc, var_glb, 1, mpi_p, MPI_MIN, MPI_COMM_WORLD, ierr)
#endif

    end subroutine s_mpi_allreduce_min

    !> Reduce a local real vector to its elementwise global minimum across all MPI ranks.
    impure subroutine s_mpi_allreduce_min_vec(var_loc, var_glb)

        real(wp), dimension(:), intent(in)  :: var_loc
        real(wp), dimension(:), intent(out) :: var_glb

#ifdef MFC_MPI
        integer :: ierr  !< Generic flag used to identify and report MPI errors

        call MPI_ALLREDUCE(var_loc, var_glb, size(var_loc), mpi_p, MPI_MIN, MPI_COMM_WORLD, ierr)
#else
        var_glb = var_loc
#endif

    end subroutine s_mpi_allreduce_min_vec

    !> Reduce a local real value to its global maximum across all MPI ranks.
    impure subroutine s_mpi_allreduce_max(var_loc, var_glb)

        real(wp), intent(in)  :: var_loc
        real(wp), intent(out) :: var_glb

#ifdef MFC_MPI
        integer :: ierr  !< Generic flag used to identify and report MPI errors

        call MPI_ALLREDUCE(var_loc, var_glb, 1, mpi_p, MPI_MAX, MPI_COMM_WORLD, ierr)
#endif

    end subroutine s_mpi_allreduce_max

    !> Reduce a local real value to its global minimum across all ranks
    impure subroutine s_mpi_reduce_min(var_loc)

        real(wp), intent(inout) :: var_loc

#ifdef MFC_MPI
        integer  :: ierr  !< Generic flag used to identify and report MPI errors
        real(wp) :: var_glb

        call MPI_REDUCE(var_loc, var_glb, 1, mpi_p, MPI_MIN, 0, MPI_COMM_WORLD, ierr)

        call MPI_BCAST(var_glb, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)

        var_loc = var_glb
#endif

    end subroutine s_mpi_reduce_min

    !> Reduce a 2-element variable to its global maximum value with the owning processor rank (MPI_MAXLOC).
    !> Reduce a local value to its global maximum with location (rank) across all ranks
    impure subroutine s_mpi_reduce_maxloc(var_loc)

        real(wp), dimension(2), intent(inout) :: var_loc

#ifdef MFC_MPI
        integer                :: ierr     !< Generic flag used to identify and report MPI errors
        real(wp), dimension(2) :: var_glb  !< Reduced (max value, rank) pair
        call MPI_REDUCE(var_loc, var_glb, 1, mpi_2p, MPI_MAXLOC, 0, MPI_COMM_WORLD, ierr)

        call MPI_BCAST(var_glb, 1, mpi_2p, 0, MPI_COMM_WORLD, ierr)

        var_loc = var_glb
#endif

    end subroutine s_mpi_reduce_maxloc

    !> The subroutine terminates the MPI execution environment.
    impure subroutine s_mpi_abort(prnt, code)

        character(len=*), intent(in), optional :: prnt
        integer, intent(in), optional          :: code

#ifdef MFC_MPI
        integer :: ierr  !< Generic flag used to identify and report MPI errors
#endif

        if (present(prnt)) then
            print *, prnt
            call flush (6)
        end if

#ifndef MFC_MPI
        if (present(code)) then
            stop code
        else
            stop 1
        end if
#else
        if (present(code)) then
            call MPI_ABORT(MPI_COMM_WORLD, code, ierr)
        else
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        end if
#endif

    end subroutine s_mpi_abort

    !> Abort naming the file if MPI_FILE_OPEN failed. File errors default to MPI_ERRORS_RETURN, so an unchecked failure carries on
    !! with an invalid handle and crashes somewhere unrelated.
    impure subroutine s_check_mpi_file_open(ierr, file_loc)

        integer, intent(in)          :: ierr
        character(len=*), intent(in) :: file_loc

#ifdef MFC_MPI
        character(len=MPI_MAX_ERROR_STRING) :: msg
        integer                             :: msg_len, ierr2

        if (ierr == MPI_SUCCESS) return
        call MPI_ERROR_STRING(ierr, msg, msg_len, ierr2)
        call s_mpi_abort('MPI_FILE_OPEN failed on ' // trim(file_loc) // ': ' // msg(1:msg_len))
#endif

    end subroutine s_check_mpi_file_open

    !> Halts all processes until all have reached barrier.
    impure subroutine s_mpi_barrier

#ifdef MFC_MPI
        integer :: ierr  !< Generic flag used to identify and report MPI errors

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif

    end subroutine s_mpi_barrier

    !> The subroutine finalizes the MPI execution environment.
    impure subroutine s_mpi_finalize

#ifdef MFC_MPI
        integer :: ierr  !< Generic flag used to identify and report MPI errors

        call MPI_FINALIZE(ierr)
#endif

    end subroutine s_mpi_finalize

#ifdef MFC_HALO_SINGLE_ROUND
#ifdef MFC_HALO_INDEX_CACHE
    subroutine s_get_26_cached_indices(point_idx, j, k, l, sj, sk, sl)

        $:GPU_ROUTINE(function_name='s_get_26_cached_indices', parallelism='[seq]', cray_inline=True)
        integer, intent(in)  :: point_idx
        integer, intent(out) :: j, k, l, sj, sk, sl

        j = halo_26_indices(point_idx, 1)
        k = halo_26_indices(point_idx, 2)
        l = halo_26_indices(point_idx, 3)
        sj = halo_26_indices(point_idx, 4)
        sk = halo_26_indices(point_idx, 5)
        sl = halo_26_indices(point_idx, 6)

    end subroutine s_get_26_cached_indices
#endif

    subroutine s_get_26_neighbor_indices(point_idx, j, k, l, sj, sk, sl)

        $:GPU_ROUTINE(function_name='s_get_26_neighbor_indices', parallelism='[seq]', cray_inline=True)
        integer, intent(in)  :: point_idx
        integer, intent(out) :: j, k, l, sj, sk, sl
        integer              :: neighbor_idx, local_idx, dx, dy, dz, lo, hi, mid

        lo = 1
        hi = halo_26_neighbor_count
        do while (lo < hi)
            mid = (lo + hi + 1)/2
            if (point_idx >= halo_26_point_offset(mid)) then
                lo = mid
            else
                hi = mid - 1
            end if
        end do
        neighbor_idx = lo

        dx = halo_26_dir(1, neighbor_idx)
        dy = halo_26_dir(2, neighbor_idx)
        dz = halo_26_dir(3, neighbor_idx)
        local_idx = point_idx - halo_26_point_offset(neighbor_idx)

        j = halo_26_beg(1, dx) + mod(local_idx, halo_26_extent(1, dx))
        local_idx = local_idx/halo_26_extent(1, dx)
        k = halo_26_beg(2, dy) + mod(local_idx, halo_26_extent(2, dy))
        local_idx = local_idx/halo_26_extent(2, dy)
        l = halo_26_beg(3, dz) + local_idx

        sj = j
        if (dx == -1) sj = j + buff_size
        if (dx == 1) sj = j - buff_size
        sk = k
        if (dy == -1) sk = k + buff_size
        if (dy == 1) sk = k - buff_size
        sl = l
        if (dz == -1) sl = l + buff_size
        if (dz == 1) sl = l - buff_size

    end subroutine s_get_26_neighbor_indices

#if defined(MFC_NVSHMEM) && defined(__NVCOMPILER)
    !> Decode one point of the trimmed 26-neighbor enumeration into the region cell indices, the source-cell shift, and the buffer
    !! position. The trimmed points of each neighbor are packed densely from that neighbor's full region offset, so the buffer index
    !! is the region offset plus the point's local trimmed index; the receiving rank decodes the identical sequence for its mirrored
    !! region, which is what keeps the one-sided contiguous put and the unpack aligned.
    subroutine s_get_26_trimmed_indices(point_idx, receiving, j, k, l, sj, sk, sl, buf_idx)

        $:GPU_ROUTINE(function_name='s_get_26_trimmed_indices', parallelism='[seq]', cray_inline=True)
        integer, intent(in)  :: point_idx
        logical, intent(in)  :: receiving
        integer, intent(out) :: j, k, l, sj, sk, sl, buf_idx
        integer              :: side, nb, local, dx, dy, dz, e1, e2, lo_i, hi_i, mid

        side = 1
        if (receiving) side = 2

        lo_i = 1
        hi_i = halo_26_neighbor_count
        do while (lo_i < hi_i)
            mid = (lo_i + hi_i + 1)/2
            if (point_idx >= halo_26_trim_off(side, mid)) then
                lo_i = mid
            else
                hi_i = mid - 1
            end if
        end do
        nb = lo_i
        local = point_idx - halo_26_trim_off(side, nb)

        dx = halo_26_dir(1, nb)
        dy = halo_26_dir(2, nb)
        dz = halo_26_dir(3, nb)
        e1 = halo_26_trim_hi(1, side, nb) - halo_26_trim_lo(1, side, nb) + 1
        e2 = halo_26_trim_hi(2, side, nb) - halo_26_trim_lo(2, side, nb) + 1
        j = halo_26_beg(1, dx) + halo_26_trim_lo(1, side, nb) + mod(local, e1)
        local = local/e1
        k = halo_26_beg(2, dy) + halo_26_trim_lo(2, side, nb) + mod(local, e2)
        l = halo_26_beg(3, dz) + halo_26_trim_lo(3, side, nb) + local/e2

        sj = j
        if (dx == -1) sj = j + buff_size
        if (dx == 1) sj = j - buff_size
        sk = k
        if (dy == -1) sk = k + buff_size
        if (dy == 1) sk = k - buff_size
        sl = l
        if (dz == -1) sl = l + buff_size
        if (dz == 1) sl = l - buff_size

        buf_idx = halo_26_point_offset(nb) + (point_idx - halo_26_trim_off(side, nb))

    end subroutine s_get_26_trimmed_indices
#endif

    integer function f_26_neighbor_tag(dx, dy, dz) result(tag)

        integer, intent(in) :: dx, dy, dz

        tag = 32 + (dz + 1)*9 + (dy + 1)*3 + (dx + 1)

    end function f_26_neighbor_tag

#if defined(MFC_MPI) && defined(MFC_SIMULATION) && !(defined(MFC_NVSHMEM) && defined(__NVCOMPILER))
    !> Completes the rdma_mpi 26-neighbor exchange's outstanding persistent requests. Invoked
    !! lazily at the next exchange's post phase (and before any other use of MFC_HALO_SEND /
    !! MFC_HALO_RECV, and at shutdown) so the host never blocks for remote data mid-exchange.
    subroutine s_halo_26_rdma_complete()

    integer :: ierr, complete_slot

    if (.not. halo_26_rdma_pending) return
    complete_slot = halo_26_rdma_pending_slot
    ! MFC_COMM_BEGIN name=halo_26_rdma_complete mode=mpi_only pattern=batch
    call MPI_WAITALL(halo_26_neighbor_count, halo_26_rdma_recv((complete_slot - 1)*halo_26_max_neighbors + 1), &
                     & MPI_STATUSES_IGNORE, ierr)
    call MPI_WAITALL(halo_26_neighbor_count, halo_26_rdma_send((complete_slot - 1)*halo_26_max_neighbors + 1), &
                     & MPI_STATUSES_IGNORE, ierr)
    ! MFC_COMM_END
    halo_26_rdma_pending = .false.

end subroutine s_halo_26_rdma_complete
#endif

subroutine s_mpi_sendrecv_variables_buffers_26(q_comm, nVar, pb_in, mv_in, q_T_sf, phase)

    type(scalar_field), dimension(1:), intent(inout) :: q_comm
    real(stp), optional, dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_in, mv_in
    type(scalar_field), optional, intent(inout) :: q_T_sf
    integer, intent(in) :: nVar
    integer, optional, intent(in) :: phase
    integer :: i, j, k, l, q, r, r0, point_idx, sj, sk, sl, node_idx
    integer :: halo_phase
    logical :: qbmm_comm, chem_diff_comm, do_post, do_wait

#ifdef MFC_MPI
    integer :: neighbor_idx, dx, dy, dz, partner, buffer_offset, buffer_count
    integer :: send_tag, recv_tag, ierr
    integer :: rdma_slot, rdma_scan
#if defined(MFC_NVSHMEM) && defined(__NVCOMPILER)
    integer(c_int)    :: nvshmem_status
    integer(c_size_t) :: nvshmem_26_recv_offset
#endif
#endif

    if (halo_26_neighbor_count == 0) return

    halo_phase = 0
    if (present(phase)) halo_phase = phase
#if defined(MFC_NVSHMEM) && defined(__NVCOMPILER)
#else
    if (.not. use_rdma_transport) then
        if (halo_phase == 2) return
        halo_phase = 0
    end if
#endif
    do_post = halo_phase /= 2
    do_wait = halo_phase /= 1

    qbmm_comm = present(pb_in) .and. present(mv_in) .and. qbmm .and. .not. polytropic
    chem_diff_comm = present(q_T_sf) .and. (heat_conduction .or. (chemistry .and. (chem_params%diffusion &
                             & .or. exchange_all_chemistry_temperatures)))

    if (qbmm_comm) then
        v_size = nVar + 2*nb*nnode
    else if (chem_diff_comm) then
        v_size = nVar + 1
    else
        v_size = nVar
    end if
#if defined(MFC_MPI) && defined(MFC_SIMULATION) && !(defined(MFC_NVSHMEM) && defined(__NVCOMPILER))
        ! rdma_mpi GPU-direct arm: lazily complete the previous exchange's requests here,
        ! before this exchange's pack overwrites MFC_HALO_SEND.
        if (use_rdma_transport .and. do_post) call s_halo_26_rdma_complete()
#endif
    if (do_post) then
        $:GPU_UPDATE(device='[v_size]')
#if defined(MFC_NVSHMEM) && defined(__NVCOMPILER) && defined(MFC_HALO_SINGLE_ROUND)
        ! Deferred completion of the previous 26-neighbor exchange: retire its puts
        ! here, before this pack overwrites the source buffer. The device-epoch mode's
        ! k_put_26_dev block puts complete with their kernel, so stream order alone
        ! protects the source buffer and the retire is skipped.
        if (.not. halo_26_epoch_device) then
            nvshmem_status = mfc_nvshmem_retire_26()
            if (nvshmem_status /= 0) call s_mpi_abort('NVSHMEM 26-neighbor put retirement failed.')
        end if
#endif
        call nvtxStartRange("RHS-COMM-26-PACKBUF")
        $:GPU_PARALLEL_LOOP(private='[i, j, k, l, r, r0, sj, sk, sl]', firstprivate='[nVar]')
        do point_idx = 0, MFC_26_PACK_POINTS - 1
#if defined(MFC_NVSHMEM) && defined(__NVCOMPILER)
            if (halo_26_trim_active) then
                call s_get_26_trimmed_indices(point_idx, .false., j, k, l, sj, sk, sl, r0)
            else
#endif
                call MFC_26_POINT_INDICES(point_idx, j, k, l, sj, sk, sl)
                r0 = point_idx
#if defined(MFC_NVSHMEM) && defined(__NVCOMPILER)
            end if
#endif
            $:GPU_LOOP(parallelism='[seq]')
            do i = 1, nVar
                r = (i - 1) + v_size*r0
                MFC_HALO_SEND(r) = real(q_comm(i)%sf(sj, sk, sl), kind=halo_kind)
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        if (chem_diff_comm) then
            $:GPU_PARALLEL_LOOP(private='[j, k, l, r, sj, sk, sl]', firstprivate='[nVar]')
            do point_idx = 0, halo_26_total_points - 1
                call MFC_26_POINT_INDICES(point_idx, j, k, l, sj, sk, sl)
                r = nVar + v_size*point_idx
                MFC_HALO_SEND(r) = real(q_T_sf%sf(sj, sk, sl), kind=halo_kind)
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        if (qbmm_comm) then
            $:GPU_PARALLEL_LOOP(private='[i, j, k, l, q, r, sj, sk, sl, node_idx]', firstprivate='[nVar]')
            do point_idx = 0, halo_26_total_points - 1
                call MFC_26_POINT_INDICES(point_idx, j, k, l, sj, sk, sl)
                $:GPU_LOOP(parallelism='[seq]')
                do i = 0, nb*nnode - 1
                    node_idx = mod(i, nnode) + 1
                    q = i/nnode + 1
                    r = nVar + i + v_size*point_idx
                    MFC_HALO_SEND(r) = real(pb_in(sj, sk, sl, node_idx, q), kind=halo_kind)
                    r = nVar + nb*nnode + i + v_size*point_idx
                    MFC_HALO_SEND(r) = real(mv_in(sj, sk, sl, node_idx, q), kind=halo_kind)
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if
        call nvtxEndRange
    end if

#ifdef MFC_MPI
#ifdef MFC_SIMULATION
#if defined(MFC_NVSHMEM) && defined(__NVCOMPILER)
    #:if defined('MFC_SIMULATION')
#include "nvshmem_halo_26_exchange.inc"
    #:endif
#else
    ! MFC_COMM_BEGIN name=halo_26 mode=device pattern=batch
    #:for rdma_mpi_enabled in [False, True]
        if (use_rdma_transport .eqv. ${'.true.' if rdma_mpi_enabled else '.false.'}$) then
            #:if rdma_mpi_enabled
                if (do_post) then
                    $:GPU_WAIT()
                end if
                #:call GPU_HOST_DATA(use_device_addr='[MFC_HALO_SEND, MFC_HALO_RECV]', no_async_wait=True)
                    call nvtxStartRange("RHS-COMM-26-SENDRECV-RDMA")
                    if (do_post) then
                        ! Persistent requests, one set per distinct v_size; they are started
                        ! here and completed lazily at the next exchange's post phase.
                        rdma_slot = 0
                        do rdma_scan = 1, halo_26_rdma_nvsizes
                            if (halo_26_rdma_vsize(rdma_scan) == v_size) rdma_slot = rdma_scan
                        end do
                        if (rdma_slot == 0) then
                            halo_26_rdma_nvsizes = halo_26_rdma_nvsizes + 1
                            rdma_slot = halo_26_rdma_nvsizes
                            halo_26_rdma_vsize(rdma_slot) = v_size
                            do neighbor_idx = 1, halo_26_neighbor_count
                                dx = halo_26_dir(1, neighbor_idx)
                                dy = halo_26_dir(2, neighbor_idx)
                                dz = halo_26_dir(3, neighbor_idx)
                                partner = neighbor_ranks(dx, dy, dz)
                                buffer_offset = v_size*halo_26_point_offset(neighbor_idx)
                                buffer_count = v_size*halo_26_point_count(neighbor_idx)
                                recv_tag = f_26_neighbor_tag(dx, dy, dz)
                                call MPI_RECV_INIT(MFC_HALO_RECV(buffer_offset), buffer_count, MFC_MPI_HALO_P, partner, recv_tag, &
                                                   & MPI_COMM_WORLD, &
                                                   & halo_26_rdma_recv((rdma_slot - 1)*halo_26_max_neighbors + neighbor_idx), ierr)
                                send_tag = f_26_neighbor_tag(-dx, -dy, -dz)
                                call MPI_SEND_INIT(MFC_HALO_SEND(buffer_offset), buffer_count, MFC_MPI_HALO_P, partner, send_tag, &
                                                   & MPI_COMM_WORLD, &
                                                   & halo_26_rdma_send((rdma_slot - 1)*halo_26_max_neighbors + neighbor_idx), ierr)
                            end do
                        end if
                        call MPI_STARTALL(halo_26_neighbor_count, halo_26_rdma_recv((rdma_slot - 1)*halo_26_max_neighbors + 1), &
                                          & ierr)
                        call MPI_STARTALL(halo_26_neighbor_count, halo_26_rdma_send((rdma_slot - 1)*halo_26_max_neighbors + 1), &
                                          & ierr)
                        halo_26_rdma_pending = .true.
                        halo_26_rdma_pending_slot = rdma_slot
                    end if
                    if (do_wait) then
                        ! The host waits only for the receives: UCX's GPU-direct protocols
                        ! (rendezvous handshake and host-staged eager alike) need receiver
                        ! progress to land the data, so the wait cannot move to the device.
                        ! The sends stay outstanding and complete lazily at the next
                        ! exchange's post phase, and no queue drain follows the unpack.
                        if (halo_26_rdma_pending) then
                            call MPI_WAITALL(halo_26_neighbor_count, &
                                             & halo_26_rdma_recv((halo_26_rdma_pending_slot - 1)*halo_26_max_neighbors + 1), &
                                             & MPI_STATUSES_IGNORE, ierr)
                        end if
                    end if
                    call nvtxEndRange
                #:endcall GPU_HOST_DATA
            #:else
                call nvtxStartRange("RHS-COMM-26-DEV2HOST")
                $:GPU_UPDATE(host='[MFC_HALO_SEND]')
                call nvtxEndRange

                call nvtxStartRange("RHS-COMM-26-SENDRECV")
                do neighbor_idx = 1, halo_26_neighbor_count
                    dx = halo_26_dir(1, neighbor_idx)
                    dy = halo_26_dir(2, neighbor_idx)
                    dz = halo_26_dir(3, neighbor_idx)
                    partner = neighbor_ranks(dx, dy, dz)
                    buffer_offset = v_size*halo_26_point_offset(neighbor_idx)
                    buffer_count = v_size*halo_26_point_count(neighbor_idx)
                    recv_tag = f_26_neighbor_tag(dx, dy, dz)
                    call MPI_IRECV(MFC_HALO_RECV(buffer_offset), buffer_count, MFC_MPI_HALO_P, partner, recv_tag, MPI_COMM_WORLD, &
                                   & halo_26_recv_req(neighbor_idx), ierr)
                end do
                do neighbor_idx = 1, halo_26_neighbor_count
                    dx = halo_26_dir(1, neighbor_idx)
                    dy = halo_26_dir(2, neighbor_idx)
                    dz = halo_26_dir(3, neighbor_idx)
                    partner = neighbor_ranks(dx, dy, dz)
                    buffer_offset = v_size*halo_26_point_offset(neighbor_idx)
                    buffer_count = v_size*halo_26_point_count(neighbor_idx)
                    send_tag = f_26_neighbor_tag(-dx, -dy, -dz)
                    call MPI_ISEND(MFC_HALO_SEND(buffer_offset), buffer_count, MFC_MPI_HALO_P, partner, send_tag, MPI_COMM_WORLD, &
                                   & halo_26_send_req(neighbor_idx), ierr)
                end do
                call MPI_WAITALL(halo_26_neighbor_count, halo_26_recv_req, MPI_STATUSES_IGNORE, ierr)
                call MPI_WAITALL(halo_26_neighbor_count, halo_26_send_req, MPI_STATUSES_IGNORE, ierr)
                call nvtxEndRange

                call nvtxStartRange("RHS-COMM-26-HOST2DEV")
                $:GPU_UPDATE(device='[MFC_HALO_RECV]')
                call nvtxEndRange
            #:endif
        end if
    #:endfor
    ! MFC_COMM_END
#endif
#endif
#endif

    if (.not. do_wait) return

    call nvtxStartRange("RHS-COMM-26-UNPACKBUF")
    $:GPU_PARALLEL_LOOP(private='[i, j, k, l, r, r0, sj, sk, sl]', firstprivate='[nVar]')
    do point_idx = 0, MFC_26_RECV_POINTS - 1
#if defined(MFC_NVSHMEM) && defined(__NVCOMPILER)
        if (halo_26_trim_active) then
            call s_get_26_trimmed_indices(point_idx, .true., j, k, l, sj, sk, sl, r0)
        else
#endif
            call MFC_26_POINT_INDICES(point_idx, j, k, l, sj, sk, sl)
            r0 = point_idx
#if defined(MFC_NVSHMEM) && defined(__NVCOMPILER)
        end if
#endif
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, nVar
            r = (i - 1) + v_size*r0
            q_comm(i)%sf(j, k, l) = real(MFC_HALO_26_RECV(r), kind=stp)
#if defined(__INTEL_COMPILER)
            if (ieee_is_nan(q_comm(i)%sf(j, k, l))) then
                print *, "Error", j, k, l, i
                call s_mpi_abort("NaN(s) in recv")
            end if
#endif
        end do
    end do
    $:END_GPU_PARALLEL_LOOP()

    if (chem_diff_comm) then
        $:GPU_PARALLEL_LOOP(private='[j, k, l, r, sj, sk, sl]', firstprivate='[nVar]')
        do point_idx = 0, halo_26_total_points - 1
            call MFC_26_POINT_INDICES(point_idx, j, k, l, sj, sk, sl)
            r = nVar + v_size*point_idx
            q_T_sf%sf(j, k, l) = real(MFC_HALO_26_RECV(r), kind=stp)
        end do
        $:END_GPU_PARALLEL_LOOP()
    end if

    if (qbmm_comm) then
        $:GPU_PARALLEL_LOOP(private='[i, j, k, l, q, r, sj, sk, sl, node_idx]', firstprivate='[nVar]')
        do point_idx = 0, halo_26_total_points - 1
            call MFC_26_POINT_INDICES(point_idx, j, k, l, sj, sk, sl)
            $:GPU_LOOP(parallelism='[seq]')
            do i = 0, nb*nnode - 1
                node_idx = mod(i, nnode) + 1
                q = i/nnode + 1
                r = nVar + i + v_size*point_idx
                pb_in(j, k, l, node_idx, q) = real(MFC_HALO_26_RECV(r), kind=stp)
                r = nVar + nb*nnode + i + v_size*point_idx
                mv_in(j, k, l, node_idx, q) = real(MFC_HALO_26_RECV(r), kind=stp)
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
    end if
    call nvtxEndRange

#ifdef MFC_MPI
#ifdef MFC_SIMULATION
#if defined(MFC_NVSHMEM) && defined(__NVCOMPILER)
    #:if defined('MFC_SIMULATION')
#include "nvshmem_halo_26_complete.inc"
    #:endif
#endif
#endif
#endif

end subroutine s_mpi_sendrecv_variables_buffers_26

#if defined(MFC_NVSHMEM) && defined(__NVCOMPILER) && defined(MFC_MPI)
!> Shrink the 26-neighbor exchange to the sub-boxes the receiving ranks' active-box windows consume. <share> carries, per rank, each
!! axis' window bounds followed by its edge-proximity flags (window within buff_size of the low and high domain edges), gathered by
!! m_active_box at every window update; every rank derives every partner's sub-boxes from the same integers, so the packing
!! enumeration a sender uses and the unpacking enumeration its receiver uses stay identical without negotiating message sizes. The
!! inboxes, region offsets, channels, and flow-control epochs are unchanged: a trimmed put simply writes a shorter contiguous run at
!! the same region offset, and regions no consumer reads are never written or read. Eligibility excludes every feature whose kernels
!! read exchanged ghost cells beyond the fused columns' stencil reach; those keep the full-size exchange.
subroutine s_mpi_apply_26_trim(share)

    integer, intent(in)  :: share(1:,1:)
    integer              :: n, a, side, rnk, d, lo, hi, cnt, off, flag
    integer(c_int)       :: nvshmem_status
    integer(c_long_long) :: trim_counts(halo_26_max_neighbors)

    if (.not. halo_26_trim_requested) return
    if (.not. (num_procs > 1 .and. .not. qbmm .and. .not. chemistry .and. .not. igr .and. .not. ib .and. .not. surface_tension &
        & .and. .not. bubbles_lagrange .and. .not. bubbles_euler .and. .not. hypoelasticity .and. .not. mhd &
        & .and. .not. relativity .and. .not. acoustic_source .and. .not. cont_damage .and. .not. cyl_coord .and. .not. viscous)) &
        & return

    do side = 1, 2
        off = 0
        do n = 1, halo_26_neighbor_count
            if (side == 1) then
                rnk = neighbor_ranks(halo_26_dir(1, n), halo_26_dir(2, n), halo_26_dir(3, n)) + 1
            else
                rnk = proc_rank + 1
            end if
            cnt = 1
            do a = 1, 3
                d = halo_26_dir(a, n)
                if (share(2*a - 1, rnk) > share(2*a, rnk)) then
                    ! Inverted window bounds mark an empty consumer window. Until some rank
                    ! engages (share slot 21), empty windows compute the full domain and the
                    ! halos stay untrimmed; afterwards an empty window computes nothing and
                    ! trims to zero.
                    if (any(share(21,:) == 1)) then
                        lo = 0
                        hi = -1
                    else
                        lo = 0
                        hi = halo_26_extent(a, d) - 1
                    end if
                else if (d == 0) then
                    ! Transverse axis clamps the consumer's window plus the stencil margin
                    ! into the shared region extent (local coordinates).
                    lo = max(0, share(2*a - 1, rnk) - buff_size - halo_26_beg(a, 0))
                    hi = min(halo_26_extent(a, 0) - 1, share(2*a, rnk) + buff_size - halo_26_beg(a, 0))
                else
                    ! Normal axis sends the whole ghost slab when the consumer's window
                    ! reaches within the stencil margin of the matching domain edge. The
                    ! pack side feeds the partner's ghost slab opposite to the direction.
                    ! The share payload interleaves the edge flags per axis as lo, hi.
                    flag = merge(2*a + 5, 2*a + 6, (side == 1) .eqv. (d == 1))
                    if (share(flag, rnk) == 0) cnt = 0
                    lo = 0
                    hi = halo_26_extent(a, d) - 1
                end if
                halo_26_trim_lo(a, side, n) = lo
                halo_26_trim_hi(a, side, n) = hi
                if (lo > hi) cnt = 0
                cnt = cnt*(hi - lo + 1)
            end do
            halo_26_trim_cnt(side, n) = cnt
            halo_26_trim_off(side, n) = off
            off = off + cnt
        end do
        halo_26_trim_off(side, 0) = 0
        if (side == 1) then
            halo_26_pack_points = off
        else
            halo_26_recv_points = off
        end if
    end do
    halo_26_trim_active = .true.
    $:GPU_UPDATE(device='[halo_26_trim_lo, halo_26_trim_hi, halo_26_trim_cnt, halo_26_trim_off, halo_26_pack_points, &
                 & halo_26_recv_points, halo_26_trim_active]')

    do n = 1, halo_26_neighbor_count
        trim_counts(n) = int(halo_26_trim_cnt(1, n), c_long_long)
    end do
    nvshmem_status = mfc_nvshmem_configure_26(int(halo_26_neighbor_count, c_int), nvshmem_26_dst_pe, nvshmem_26_src_points, &
        & nvshmem_26_dst_points, trim_counts, nvshmem_26_channel)
    if (nvshmem_status /= 0) call s_mpi_abort('Unable to reconfigure the NVSHMEM 26-neighbor trim.')

end subroutine s_mpi_apply_26_trim
#endif
#endif

#ifdef MFC_HALO_DEFERRED_COMPLETION
!> Complete the previous exchange's outstanding sends or puts before buff_send and buff_recv are reused.
subroutine s_complete_deferred_halo_exchanges()

#if defined(MFC_NVSHMEM) && defined(__NVCOMPILER)
    integer(c_int) :: nvshmem_status

#ifdef MFC_HALO_SINGLE_ROUND
    ! The device-epoch mode's block puts complete with their kernel, so only the other modes retire.
    if (.not. halo_26_epoch_device) then
        nvshmem_status = mfc_nvshmem_retire_26()
        if (nvshmem_status /= 0) call s_mpi_abort('NVSHMEM 26-neighbor put retirement failed.')
    end if
#endif
    nvshmem_status = mfc_nvshmem_retire()
    if (nvshmem_status /= 0) call s_mpi_abort('NVSHMEM halo put retirement failed.')
#else
    call s_halo_26_rdma_complete()
#endif

end subroutine s_complete_deferred_halo_exchanges
#endif

!> The goal of this procedure is to populate the buffers of the cell-average conservative variables by communicating with the
!! neighboring processors.
subroutine s_mpi_sendrecv_variables_buffers(q_comm, mpi_dir, pbc_loc, nVar, pb_in, mv_in, q_T_sf)

    type(scalar_field), dimension(1:), intent(inout) :: q_comm
    real(stp), optional, dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_in, mv_in
    integer, intent(in) :: mpi_dir, pbc_loc, nVar
    integer :: depth                !< Ghost layers exchanged along mpi_dir
    integer :: i, j, k, l, r, q     !< Generic loop iterators
    integer :: buffer_counts(1:3), buffer_count
    type(int_bounds_info) :: boundary_conditions(1:3)
    integer :: beg_end(1:2), grid_dims(1:3)
    integer :: dst_proc, src_proc, recv_tag, send_tag
    logical :: beg_end_geq_0, qbmm_comm, chem_diff_comm
    integer :: pack_offset, unpack_offset
    integer :: t, ix, iy, iz, s, v  !< Flat-loop decode temporaries
    type(scalar_field), optional, intent(inout) :: q_T_sf

#if defined(MFC_NVSHMEM) && defined(__NVCOMPILER)
    integer(c_int)    :: nvshmem_status
    integer(c_size_t) :: nvshmem_recv_offset
    integer           :: send_channel, recv_channel
#endif

#ifdef MFC_MPI
    integer :: ierr  !< Generic flag used to identify and report MPI errors

#ifdef MFC_HALO_DEFERRED_COMPLETION
    call s_complete_deferred_halo_exchanges()
#endif

    call nvtxStartRange("RHS-COMM-PACKBUF")

    depth = buff_size
    qbmm_comm = .false.
    chem_diff_comm = .false.

    if (present(pb_in) .and. present(mv_in) .and. qbmm .and. .not. polytropic) then
        qbmm_comm = .true.
        v_size = nVar + 2*nb*nnode
        buffer_counts = (/buff_size*v_size*(n + 1)*(p + 1), buff_size*v_size*(m + 2*buff_size + 1)*(p + 1), &
                         & buff_size*v_size*(m + 2*buff_size + 1)*(n + 2*buff_size + 1)/)
    else if (present(q_T_sf) .and. (heat_conduction .or. (chemistry .and. (chem_params%diffusion &
             & .or. exchange_all_chemistry_temperatures)))) then
        ! Consumers that convert over ghost-inclusive bounds request temperature exchange for every chemistry run.
        ! The temperature Newton guess must be valid at rank seams even when diffusion is disabled:
        ! an unexchanged seam ghost is an uninitialized guess -> NaN T/pres/c in the output
        chem_diff_comm = .true.
        v_size = nVar + 1
        buffer_counts = (/buff_size*v_size*(n + 1)*(p + 1), buff_size*v_size*(m + 2*buff_size + 1)*(p + 1), &
                         & buff_size*v_size*(m + 2*buff_size + 1)*(n + 2*buff_size + 1)/)
    else
        v_size = nVar
        buffer_counts = (/depth*v_size*(n + 1)*(p + 1), depth*v_size*(m + 2*buff_size + 1)*(p + 1), &
                         & depth*v_size*(m + 2*buff_size + 1)*(n + 2*buff_size + 1)/)
    end if

    $:GPU_UPDATE(device='[v_size]')

    buffer_count = buffer_counts(mpi_dir)
    boundary_conditions = (/bc_x, bc_y, bc_z/)
    beg_end = (/boundary_conditions(mpi_dir)%beg, boundary_conditions(mpi_dir)%end/)
    beg_end_geq_0 = beg_end(max(pbc_loc, 0) - pbc_loc + 1) >= 0

    ! Implements: pbc_loc bc_x >= 0 -> [send/recv]_tag [dst/src]_proc -1 (=0) 0 -> [1,0] [0,0] | 0 0 [1,0] [beg,beg] -1 (=0) 1
    ! -> [0,0] [1,0] | 0 1 [0,0] [end,beg] +1 (=1) 0 -> [0,1] [1,1] | 1 0 [0,1] [end,end] +1 (=1) 1 -> [1,1] [0,1] | 1 1 [1,1]
    ! [beg,end]

    send_tag = f_logical_to_int(.not. f_xor(beg_end_geq_0, pbc_loc == 1))
    recv_tag = f_logical_to_int(pbc_loc == 1)

    dst_proc = beg_end(1 + f_logical_to_int(f_xor(pbc_loc == 1, beg_end_geq_0)))
    src_proc = beg_end(1 + f_logical_to_int(pbc_loc == 1))

    grid_dims = (/m, n, p/)

    pack_offset = 0
    if (f_xor(pbc_loc == 1, beg_end_geq_0)) then
        pack_offset = grid_dims(mpi_dir) - depth + 1
    end if

    unpack_offset = 0
    if (pbc_loc == 1) then
        unpack_offset = grid_dims(mpi_dir) + depth + 1
    end if

    ! Pack Buffer to Send
    #:for mpi_dir in [1, 2, 3]
        if (mpi_dir == ${mpi_dir}$) then
            #:if mpi_dir == 1
                ! q enumerates q_comm's nVar variables; qbmm and chemistry fill the slots beyond nVar below.
                $:GPU_PARALLEL_LOOP(private='[i, q, t, ix, iy, iz]')
                do q = 0, (buffer_count/v_size)*nVar - 1
                    i = mod(q, nVar) + 1
                    t = q/nVar
                    r = (i - 1) + v_size*t
                    ix = mod(t, (depth)); t = t/(depth)
                    iy = mod(t, (n + 1)); iz = t/(n + 1)
                    MFC_HALO_SEND(r) = real(q_comm(i)%sf(ix + (pack_offset), iy + (0), iz + (0)), kind=halo_kind)
                end do
                $:END_GPU_PARALLEL_LOOP()

                if (chem_diff_comm) then
                    $:GPU_PARALLEL_LOOP(private='[t, ix, iy, iz]')
                    do s = 0, buffer_count/v_size - 1
                        ix = mod(s, (buff_size)); t = s/(buff_size)
                        iy = mod(t, (n + 1)); iz = t/(n + 1)
                        MFC_HALO_SEND(nVar + v_size*s) = real(q_T_sf%sf(ix + (pack_offset), iy + (0), iz + (0)), kind=halo_kind)
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if

                if (qbmm_comm) then
                    $:GPU_PARALLEL_LOOP(private='[t, ix, iy, iz, v, s]')
                    do r = 0, nb*nnode*(buffer_count/v_size) - 1
                        v = mod(r, nb*nnode)
                        s = r/(nb*nnode)
                        ix = mod(s, (buff_size)); t = s/(buff_size)
                        iy = mod(t, (n + 1)); iz = t/(n + 1)
                        MFC_HALO_SEND(nVar + (0) + v + v_size*s) = real(pb_in(ix + (pack_offset), iy + (0), iz + (0), mod(v, &
                                      & nnode) + 1, v/nnode + 1), kind=halo_kind)
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                    $:GPU_PARALLEL_LOOP(private='[t, ix, iy, iz, v, s]')
                    do r = 0, nb*nnode*(buffer_count/v_size) - 1
                        v = mod(r, nb*nnode)
                        s = r/(nb*nnode)
                        ix = mod(s, (buff_size)); t = s/(buff_size)
                        iy = mod(t, (n + 1)); iz = t/(n + 1)
                        MFC_HALO_SEND(nVar + (nb*nnode) + v + v_size*s) = real(mv_in(ix + (pack_offset), iy + (0), iz + (0), &
                                      & mod(v, nnode) + 1, v/nnode + 1), kind=halo_kind)
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
            #:elif mpi_dir == 2
                $:GPU_PARALLEL_LOOP(private='[i, q, t, ix, iy, iz]')
                do q = 0, (buffer_count/v_size)*nVar - 1
                    i = mod(q, nVar) + 1
                    t = q/nVar
                    r = (i - 1) + v_size*t
                    ix = mod(t, (m + 2*buff_size + 1)); t = t/(m + 2*buff_size + 1)
                    iy = mod(t, (depth)); iz = t/(depth)
                    MFC_HALO_SEND(r) = real(q_comm(i)%sf(ix + (-buff_size), iy + (pack_offset), iz + (0)), kind=halo_kind)
                end do
                $:END_GPU_PARALLEL_LOOP()

                if (chem_diff_comm) then
                    $:GPU_PARALLEL_LOOP(private='[t, ix, iy, iz]')
                    do s = 0, buffer_count/v_size - 1
                        ix = mod(s, (m + 2*buff_size + 1)); t = s/(m + 2*buff_size + 1)
                        iy = mod(t, (buff_size)); iz = t/(buff_size)
                        MFC_HALO_SEND(nVar + v_size*s) = real(q_T_sf%sf(ix + (-buff_size), iy + (pack_offset), iz + (0)), &
                                      & kind=halo_kind)
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if

                if (qbmm_comm) then
                    $:GPU_PARALLEL_LOOP(private='[t, ix, iy, iz, v, s]')
                    do r = 0, nb*nnode*(buffer_count/v_size) - 1
                        v = mod(r, nb*nnode)
                        s = r/(nb*nnode)
                        ix = mod(s, (m + 2*buff_size + 1)); t = s/(m + 2*buff_size + 1)
                        iy = mod(t, (buff_size)); iz = t/(buff_size)
                        MFC_HALO_SEND(nVar + (0) + v + v_size*s) = real(pb_in(ix + (-buff_size), iy + (pack_offset), iz + (0), &
                                      & mod(v, nnode) + 1, v/nnode + 1), kind=halo_kind)
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                    $:GPU_PARALLEL_LOOP(private='[t, ix, iy, iz, v, s]')
                    do r = 0, nb*nnode*(buffer_count/v_size) - 1
                        v = mod(r, nb*nnode)
                        s = r/(nb*nnode)
                        ix = mod(s, (m + 2*buff_size + 1)); t = s/(m + 2*buff_size + 1)
                        iy = mod(t, (buff_size)); iz = t/(buff_size)
                        MFC_HALO_SEND(nVar + (nb*nnode) + v + v_size*s) = real(mv_in(ix + (-buff_size), iy + (pack_offset), &
                                      & iz + (0), mod(v, nnode) + 1, v/nnode + 1), kind=halo_kind)
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
            #:else
                $:GPU_PARALLEL_LOOP(private='[i, q, t, ix, iy, iz]')
                do q = 0, (buffer_count/v_size)*nVar - 1
                    i = mod(q, nVar) + 1
                    t = q/nVar
                    r = (i - 1) + v_size*t
                    ix = mod(t, (m + 2*buff_size + 1)); t = t/(m + 2*buff_size + 1)
                    iy = mod(t, (n + 2*buff_size + 1)); iz = t/(n + 2*buff_size + 1)
                    MFC_HALO_SEND(r) = real(q_comm(i)%sf(ix + (-buff_size), iy + (-buff_size), iz + (pack_offset)), kind=halo_kind)
                end do
                $:END_GPU_PARALLEL_LOOP()

                if (chem_diff_comm) then
                    $:GPU_PARALLEL_LOOP(private='[t, ix, iy, iz]')
                    do s = 0, buffer_count/v_size - 1
                        ix = mod(s, (m + 2*buff_size + 1)); t = s/(m + 2*buff_size + 1)
                        iy = mod(t, (n + 2*buff_size + 1)); iz = t/(n + 2*buff_size + 1)
                        MFC_HALO_SEND(nVar + v_size*s) = real(q_T_sf%sf(ix + (-buff_size), iy + (-buff_size), &
                                      & iz + (pack_offset)), kind=halo_kind)
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if

                if (qbmm_comm) then
                    $:GPU_PARALLEL_LOOP(private='[t, ix, iy, iz, v, s]')
                    do r = 0, nb*nnode*(buffer_count/v_size) - 1
                        v = mod(r, nb*nnode)
                        s = r/(nb*nnode)
                        ix = mod(s, (m + 2*buff_size + 1)); t = s/(m + 2*buff_size + 1)
                        iy = mod(t, (n + 2*buff_size + 1)); iz = t/(n + 2*buff_size + 1)
                        MFC_HALO_SEND(nVar + (0) + v + v_size*s) = real(pb_in(ix + (-buff_size), iy + (-buff_size), &
                                      & iz + (pack_offset), mod(v, nnode) + 1, v/nnode + 1), kind=halo_kind)
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                    $:GPU_PARALLEL_LOOP(private='[t, ix, iy, iz, v, s]')
                    do r = 0, nb*nnode*(buffer_count/v_size) - 1
                        v = mod(r, nb*nnode)
                        s = r/(nb*nnode)
                        ix = mod(s, (m + 2*buff_size + 1)); t = s/(m + 2*buff_size + 1)
                        iy = mod(t, (n + 2*buff_size + 1)); iz = t/(n + 2*buff_size + 1)
                        MFC_HALO_SEND(nVar + (nb*nnode) + v + v_size*s) = real(mv_in(ix + (-buff_size), iy + (-buff_size), &
                                      & iz + (pack_offset), mod(v, nnode) + 1, v/nnode + 1), kind=halo_kind)
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
            #:endif
        end if
    #:endfor
    call nvtxEndRange  ! Packbuf

    ! Send/Recv
#if defined(MFC_NVSHMEM) && defined(__NVCOMPILER)
    #:if defined('MFC_SIMULATION')
#include "nvshmem_halo_axis_exchange.inc"
    #:endif
#else
    ! MFC_COMM_BEGIN name=halo_axis mode=device pattern=sendrecv
    #:for rdma_mpi in [False, True]
        if (use_rdma_transport .eqv. ${'.true.' if rdma_mpi else '.false.'}$) then
            #:if rdma_mpi
                #:call GPU_HOST_DATA(use_device_addr='[buff_send, buff_recv]')
                    call nvtxStartRange("RHS-COMM-SENDRECV-RDMA")

                    call MPI_SENDRECV(buff_send, buffer_count, mpi_p, dst_proc, send_tag, buff_recv, buffer_count, mpi_p, &
                                      & src_proc, recv_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                    call nvtxEndRange  ! RHS-MPI-SENDRECV-(NO)-RDMA
                #:endcall GPU_HOST_DATA
                $:GPU_WAIT()
            #:else
                call nvtxStartRange("RHS-COMM-DEV2HOST")
                $:GPU_UPDATE(host='[buff_send]')
                call nvtxEndRange
                call nvtxStartRange("RHS-COMM-SENDRECV-NO-RMDA")

                call MPI_SENDRECV(buff_send, buffer_count, mpi_p, dst_proc, send_tag, buff_recv, buffer_count, mpi_p, src_proc, &
                                  & recv_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                call nvtxEndRange  ! RHS-MPI-SENDRECV-(NO)-RDMA

                call nvtxStartRange("RHS-COMM-HOST2DEV")
                $:GPU_UPDATE(device='[buff_recv]')
                call nvtxEndRange
            #:endif
        end if
    #:endfor
    ! MFC_COMM_END
#endif

    ! Unpack Received Buffer
    call nvtxStartRange("RHS-COMM-UNPACKBUF")
    #:for mpi_dir in [1, 2, 3]
        if (mpi_dir == ${mpi_dir}$) then
            #:if mpi_dir == 1
                ! Same enumeration as the pack; qbmm and chemistry unpack the slots beyond nVar below.
                $:GPU_PARALLEL_LOOP(private='[i, q, t, ix, iy, iz]')
                do q = 0, (buffer_count/v_size)*nVar - 1
                    i = mod(q, nVar) + 1
                    t = q/nVar
                    r = (i - 1) + v_size*t
                    ix = mod(t, (depth)); t = t/(depth)
                    iy = mod(t, (n + 1)); iz = t/(n + 1)
                    q_comm(i)%sf(ix + (unpack_offset - depth), iy + (0), iz + (0)) = MFC_ABSOLUTE_HALO_VALUE(r)
#if defined(__INTEL_COMPILER)
                    if (ieee_is_nan(q_comm(i)%sf(ix + (unpack_offset - depth), iy + (0), iz + (0)))) then
                        print *, "Error", ix, iy, iz, i
                        call s_mpi_abort("NaN(s) in recv")
                    end if
#endif
                end do
                $:END_GPU_PARALLEL_LOOP()

                if (chem_diff_comm) then
                    $:GPU_PARALLEL_LOOP(private='[t, ix, iy, iz]')
                    do s = 0, buffer_count/v_size - 1
                        ix = mod(s, (buff_size)); t = s/(buff_size)
                        iy = mod(t, (n + 1)); iz = t/(n + 1)
                        q_T_sf%sf(ix + (unpack_offset - buff_size), iy + (0), iz + (0)) = MFC_ABSOLUTE_HALO_VALUE(nVar + v_size*s)
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if

                if (qbmm_comm) then
                    $:GPU_PARALLEL_LOOP(private='[t, ix, iy, iz, v, s]')
                    do r = 0, nb*nnode*(buffer_count/v_size) - 1
                        v = mod(r, nb*nnode)
                        s = r/(nb*nnode)
                        ix = mod(s, (buff_size)); t = s/(buff_size)
                        iy = mod(t, (n + 1)); iz = t/(n + 1)
                        pb_in(ix + (unpack_offset - buff_size), iy + (0), iz + (0), mod(v, nnode) + 1, &
                              & v/nnode + 1) = MFC_ABSOLUTE_HALO_VALUE(nVar + (0) + v + v_size*s)
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                    $:GPU_PARALLEL_LOOP(private='[t, ix, iy, iz, v, s]')
                    do r = 0, nb*nnode*(buffer_count/v_size) - 1
                        v = mod(r, nb*nnode)
                        s = r/(nb*nnode)
                        ix = mod(s, (buff_size)); t = s/(buff_size)
                        iy = mod(t, (n + 1)); iz = t/(n + 1)
                        mv_in(ix + (unpack_offset - buff_size), iy + (0), iz + (0), mod(v, nnode) + 1, &
                              & v/nnode + 1) = MFC_ABSOLUTE_HALO_VALUE(nVar + (nb*nnode) + v + v_size*s)
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
            #:elif mpi_dir == 2
                $:GPU_PARALLEL_LOOP(private='[i, q, t, ix, iy, iz]')
                do q = 0, (buffer_count/v_size)*nVar - 1
                    i = mod(q, nVar) + 1
                    t = q/nVar
                    r = (i - 1) + v_size*t
                    ix = mod(t, (m + 2*buff_size + 1)); t = t/(m + 2*buff_size + 1)
                    iy = mod(t, (depth)); iz = t/(depth)
                    q_comm(i)%sf(ix + (-buff_size), iy + (unpack_offset - depth), iz + (0)) = MFC_ABSOLUTE_HALO_VALUE(r)
#if defined(__INTEL_COMPILER)
                    if (ieee_is_nan(q_comm(i)%sf(ix + (-buff_size), iy + (unpack_offset - depth), iz + (0)))) then
                        print *, "Error", ix, iy, iz, i
                        call s_mpi_abort("NaN(s) in recv")
                    end if
#endif
                end do
                $:END_GPU_PARALLEL_LOOP()

                if (chem_diff_comm) then
                    $:GPU_PARALLEL_LOOP(private='[t, ix, iy, iz]')
                    do s = 0, buffer_count/v_size - 1
                        ix = mod(s, (m + 2*buff_size + 1)); t = s/(m + 2*buff_size + 1)
                        iy = mod(t, (buff_size)); iz = t/(buff_size)
                        q_T_sf%sf(ix + (-buff_size), iy + (unpack_offset - buff_size), &
                                  & iz + (0)) = MFC_ABSOLUTE_HALO_VALUE(nVar + v_size*s)
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if

                if (qbmm_comm) then
                    $:GPU_PARALLEL_LOOP(private='[t, ix, iy, iz, v, s]')
                    do r = 0, nb*nnode*(buffer_count/v_size) - 1
                        v = mod(r, nb*nnode)
                        s = r/(nb*nnode)
                        ix = mod(s, (m + 2*buff_size + 1)); t = s/(m + 2*buff_size + 1)
                        iy = mod(t, (buff_size)); iz = t/(buff_size)
                        pb_in(ix + (-buff_size), iy + (unpack_offset - buff_size), iz + (0), mod(v, nnode) + 1, &
                              & v/nnode + 1) = MFC_ABSOLUTE_HALO_VALUE(nVar + (0) + v + v_size*s)
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                    $:GPU_PARALLEL_LOOP(private='[t, ix, iy, iz, v, s]')
                    do r = 0, nb*nnode*(buffer_count/v_size) - 1
                        v = mod(r, nb*nnode)
                        s = r/(nb*nnode)
                        ix = mod(s, (m + 2*buff_size + 1)); t = s/(m + 2*buff_size + 1)
                        iy = mod(t, (buff_size)); iz = t/(buff_size)
                        mv_in(ix + (-buff_size), iy + (unpack_offset - buff_size), iz + (0), mod(v, nnode) + 1, &
                              & v/nnode + 1) = MFC_ABSOLUTE_HALO_VALUE(nVar + (nb*nnode) + v + v_size*s)
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
            #:else
                $:GPU_PARALLEL_LOOP(private='[i, q, t, ix, iy, iz]')
                do q = 0, (buffer_count/v_size)*nVar - 1
                    i = mod(q, nVar) + 1
                    t = q/nVar
                    r = (i - 1) + v_size*t
                    ix = mod(t, (m + 2*buff_size + 1)); t = t/(m + 2*buff_size + 1)
                    iy = mod(t, (n + 2*buff_size + 1)); iz = t/(n + 2*buff_size + 1)
                    q_comm(i)%sf(ix + (-buff_size), iy + (-buff_size), iz + (unpack_offset - depth)) = MFC_ABSOLUTE_HALO_VALUE(r)
#if defined(__INTEL_COMPILER)
                    if (ieee_is_nan(q_comm(i)%sf(ix + (-buff_size), iy + (-buff_size), iz + (unpack_offset - depth)))) then
                        print *, "Error", ix, iy, iz, i
                        call s_mpi_abort("NaN(s) in recv")
                    end if
#endif
                end do
                $:END_GPU_PARALLEL_LOOP()

                if (chem_diff_comm) then
                    $:GPU_PARALLEL_LOOP(private='[t, ix, iy, iz]')
                    do s = 0, buffer_count/v_size - 1
                        ix = mod(s, (m + 2*buff_size + 1)); t = s/(m + 2*buff_size + 1)
                        iy = mod(t, (n + 2*buff_size + 1)); iz = t/(n + 2*buff_size + 1)
                        q_T_sf%sf(ix + (-buff_size), iy + (-buff_size), &
                                  & iz + (unpack_offset - buff_size)) = MFC_ABSOLUTE_HALO_VALUE(nVar + v_size*s)
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if

                if (qbmm_comm) then
                    $:GPU_PARALLEL_LOOP(private='[t, ix, iy, iz, v, s]')
                    do r = 0, nb*nnode*(buffer_count/v_size) - 1
                        v = mod(r, nb*nnode)
                        s = r/(nb*nnode)
                        ix = mod(s, (m + 2*buff_size + 1)); t = s/(m + 2*buff_size + 1)
                        iy = mod(t, (n + 2*buff_size + 1)); iz = t/(n + 2*buff_size + 1)
                        pb_in(ix + (-buff_size), iy + (-buff_size), iz + (unpack_offset - buff_size), mod(v, nnode) + 1, &
                              & v/nnode + 1) = MFC_ABSOLUTE_HALO_VALUE(nVar + (0) + v + v_size*s)
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                    $:GPU_PARALLEL_LOOP(private='[t, ix, iy, iz, v, s]')
                    do r = 0, nb*nnode*(buffer_count/v_size) - 1
                        v = mod(r, nb*nnode)
                        s = r/(nb*nnode)
                        ix = mod(s, (m + 2*buff_size + 1)); t = s/(m + 2*buff_size + 1)
                        iy = mod(t, (n + 2*buff_size + 1)); iz = t/(n + 2*buff_size + 1)
                        mv_in(ix + (-buff_size), iy + (-buff_size), iz + (unpack_offset - buff_size), mod(v, nnode) + 1, &
                              & v/nnode + 1) = MFC_ABSOLUTE_HALO_VALUE(nVar + (nb*nnode) + v + v_size*s)
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
            #:endif
        end if
    #:endfor
    call nvtxEndRange
#if defined(MFC_NVSHMEM) && defined(__NVCOMPILER)
    #:if defined('MFC_SIMULATION')
#include "nvshmem_halo_axis_complete.inc"
    #:endif
#endif
#endif

end subroutine s_mpi_sendrecv_variables_buffers

#ifdef MFC_IGR_BATCHED_EXCHANGE
subroutine s_pack_F_igr_face(jac_sf, face, t1, t2, depth)

    $:GPU_ROUTINE(function_name='s_pack_F_igr_face', parallelism='[seq]', cray_inline=True)
    type(scalar_field), dimension(1:), intent(in) :: jac_sf
    integer, intent(in)                           :: face, t1, t2, depth
    integer                                       :: dir, j, pack_off
    integer(kind=8)                               :: base

    dir = face/2 + 1

    pack_off = 0
    if (mod(face, 2) == 1) then
        if (dir == 1) then
            pack_off = m - depth + 1
        else if (dir == 2) then
            pack_off = n - depth + 1
        else
            pack_off = p - depth + 1
        end if
    end if

    base = igr_face_off(face) + int(depth, 8)*(int(t1, 8) + int(merge(n + 1, m + 1, dir == 1), 8)*int(t2, 8))

    if (dir == 1) then
        do j = 0, depth - 1
            igr_buff_send(base + int(j, 8)) = real(jac_sf(1)%sf(j + pack_off, t1, t2), kind=wp)
        end do
    else if (dir == 2) then
        do j = 0, depth - 1
            igr_buff_send(base + int(j, 8)) = real(jac_sf(1)%sf(t1, j + pack_off, t2), kind=wp)
        end do
    else
        do j = 0, depth - 1
            igr_buff_send(base + int(j, 8)) = real(jac_sf(1)%sf(t1, t2, j + pack_off), kind=wp)
        end do
    end if

end subroutine s_pack_F_igr_face

subroutine s_mpi_sendrecv_F_igr_batched(jac_sf, depth, active_faces)

    type(scalar_field), dimension(1:), intent(inout) :: jac_sf
    integer, intent(in)                              :: depth, active_faces
    integer                                          :: dir, loc, face, partner, nreq, ierr
    integer                                          :: requests(0:11), unpack_off
    integer(kind=8)                                  :: r, base
    integer                                          :: e, s, t1, t2
    type(int_bounds_info)                            :: bounds(1:3)
    integer                                          :: dims(1:3)

#if defined(MFC_NVSHMEM) && defined(__NVCOMPILER)
    integer(c_int)    :: nvshmem_status
    integer(c_size_t) :: nvshmem_igr_recv_offset
#endif
    if (active_faces == 0) return

    bounds(1) = bc_x; bounds(2) = bc_y; bounds(3) = bc_z
    dims = (/m, n, p/)
    nreq = 0

#if defined(MFC_NVSHMEM) && defined(__NVCOMPILER)
    #:if defined('MFC_SIMULATION')
#include "nvshmem_igr_exchange.inc"
    #:endif
#else
    ! MFC_COMM_BEGIN name=igr mode=device pattern=batch
    #:for rdma_mpi_enabled in [False, True]
        if (use_rdma_transport .eqv. ${'.true.' if rdma_mpi_enabled else '.false.'}$) then
            #:if rdma_mpi_enabled
                #:call GPU_HOST_DATA(use_device_addr='[igr_buff_send, igr_buff_recv]')
                    do dir = 1, num_dims
                        do loc = -1, 1, 2
                            face = 2*(dir - 1) + f_logical_to_int(loc == 1)
                            if (.not. btest(active_faces, face)) cycle
                            partner = bounds(dir)%beg
                            if (loc == 1) partner = bounds(dir)%end
                            call MPI_IRECV(igr_buff_recv(igr_face_off(face)), int(igr_face_cnt(face)), mpi_p, partner, &
                                           & 6000 + 2*(dir - 1) + f_logical_to_int(loc == -1), MPI_COMM_WORLD, requests(nreq), ierr)
                            nreq = nreq + 1
                            call MPI_ISEND(igr_buff_send(igr_face_off(face)), int(igr_face_cnt(face)), mpi_p, partner, &
                                           & 6000 + 2*(dir - 1) + f_logical_to_int(loc == 1), MPI_COMM_WORLD, requests(nreq), ierr)
                            nreq = nreq + 1
                        end do
                    end do
                    if (nreq > 0) call MPI_WAITALL(nreq, requests(0:nreq - 1), MPI_STATUSES_IGNORE, ierr)
                #:endcall GPU_HOST_DATA
            #:else
                $:GPU_UPDATE(host='[igr_buff_send]')
                do dir = 1, num_dims
                    do loc = -1, 1, 2
                        face = 2*(dir - 1) + f_logical_to_int(loc == 1)
                        if (.not. btest(active_faces, face)) cycle
                        partner = bounds(dir)%beg
                        if (loc == 1) partner = bounds(dir)%end
                        call MPI_IRECV(igr_buff_recv(igr_face_off(face)), int(igr_face_cnt(face)), mpi_p, partner, &
                                       & 6000 + 2*(dir - 1) + f_logical_to_int(loc == -1), MPI_COMM_WORLD, requests(nreq), ierr)
                        nreq = nreq + 1
                        call MPI_ISEND(igr_buff_send(igr_face_off(face)), int(igr_face_cnt(face)), mpi_p, partner, &
                                       & 6000 + 2*(dir - 1) + f_logical_to_int(loc == 1), MPI_COMM_WORLD, requests(nreq), ierr)
                        nreq = nreq + 1
                    end do
                end do
                if (nreq > 0) call MPI_WAITALL(nreq, requests(0:nreq - 1), MPI_STATUSES_IGNORE, ierr)
                $:GPU_UPDATE(device='[igr_buff_recv]')
            #:endif
        end if
    #:endfor
    ! MFC_COMM_END
#endif

    do dir = 1, num_dims
        do loc = -1, 1, 2
            face = 2*(dir - 1) + f_logical_to_int(loc == 1)
            if (btest(active_faces, face)) then
                base = igr_face_off(face)
                unpack_off = 0
                if (loc == 1) unpack_off = dims(dir) + depth + 1
                #:for d in [1, 2, 3]
                    if (dir == ${d}$) then
                        #:if d == 1
                            $:GPU_PARALLEL_LOOP(private='[s, t1, t2, r]')
                            do e = 0, depth*((n + 1))*(p + 1) - 1
                                s = mod(e, depth)
                                t1 = mod(e/depth, (n + 1))
                                t2 = e/(depth*(n + 1))
                                r = base + int(e, 8)
                                jac_sf(1)%sf(s - depth + unpack_off, t1, t2) = real(MFC_IGR_RECV(r), kind=stp)
                            end do
                            $:END_GPU_PARALLEL_LOOP()
                        #:elif d == 2
                            $:GPU_PARALLEL_LOOP(private='[s, t1, t2, r]')
                            do e = 0, depth*((m + 1))*(p + 1) - 1
                                s = mod(e, depth)
                                t1 = mod(e/depth, (m + 1))
                                t2 = e/(depth*(m + 1))
                                r = base + int(e, 8)
                                jac_sf(1)%sf(t1, s - depth + unpack_off, t2) = real(MFC_IGR_RECV(r), kind=stp)
                            end do
                            $:END_GPU_PARALLEL_LOOP()
                        #:else
                            $:GPU_PARALLEL_LOOP(private='[s, t1, t2, r]')
                            do e = 0, depth*((m + 1))*(n + 1) - 1
                                s = mod(e, depth)
                                t1 = mod(e/depth, (m + 1))
                                t2 = e/(depth*(m + 1))
                                r = base + int(e, 8)
                                jac_sf(1)%sf(t1, t2, s - depth + unpack_off) = real(MFC_IGR_RECV(r), kind=stp)
                            end do
                            $:END_GPU_PARALLEL_LOOP()
                        #:endif
                    end if
                #:endfor
            end if
        end do
    end do

#if defined(MFC_NVSHMEM) && defined(__NVCOMPILER)
    #:if defined('MFC_SIMULATION')
#include "nvshmem_igr_complete.inc"
    #:endif
#endif

end subroutine s_mpi_sendrecv_F_igr_batched
#endif

!> The goal of this procedure is to populate the buffers of the cell-average conservative variables by communicating with the
!! neighboring processors.
!! @param q_cons_vf Cell-average conservative variables
!! @param mpi_dir MPI communication coordinate direction
!! @param pbc_loc Processor boundary condition (PBC) location
subroutine s_mpi_reduce_beta_variables_buffers(q_comm, kahan_comp, mpi_dir, pbc_loc, nVar)

    type(scalar_field), dimension(1:), intent(inout) :: q_comm
    type(scalar_field), dimension(1:), intent(inout) :: kahan_comp
    integer, intent(in)                              :: mpi_dir, pbc_loc, nVar
    integer                                          :: i, t, ix, iy, iz, r, q  !< Generic loop iterators
    integer                                          :: lb_size
    integer                                          :: buffer_counts(1:3), buffer_count
    type(int_bounds_info)                            :: boundary_conditions(1:3)
    integer                                          :: beg_end(1:2), grid_dims(1:3)
    integer                                          :: dst_proc, src_proc, recv_tag, send_tag
    logical                                          :: replace_buff
    integer                                          :: pack_offset, unpack_offset
    real(wp)                                         :: y_kahan, t_kahan

#ifdef MFC_MPI
    integer :: ierr  !< Generic flag used to identify and report MPI errors

#ifdef MFC_HALO_DEFERRED_COMPLETION
    call s_complete_deferred_halo_exchanges()
#endif

    call nvtxStartRange("BETA-COMM-PACKBUF")

    ! Set bounds for each dimension Always include the full buffer range for each existing dimension. The Gaussian smearing
    ! kernel writes to buffer cells even at physical boundaries, and these contributions must be communicated to neighbors in
    ! other directions via ADD operations.
    comm_coords(1)%beg = -mapcells - 1
    comm_coords(1)%end = m + mapcells + 1
    comm_coords(2)%beg = merge(-mapcells - 1, 0, n > 0)
    comm_coords(2)%end = merge(n + mapcells + 1, n, n > 0)
    comm_coords(3)%beg = merge(-mapcells - 1, 0, p > 0)
    comm_coords(3)%end = merge(p + mapcells + 1, p, p > 0)

    ! Compute sizes
    comm_size(1) = comm_coords(1)%end - comm_coords(1)%beg + 1
    comm_size(2) = comm_coords(2)%end - comm_coords(2)%beg + 1
    comm_size(3) = comm_coords(3)%end - comm_coords(3)%beg + 1

    ! Buffer counts using the conditional sizes
    v_size = nVar
    lb_size = 2*(mapcells + 1)  ! Size of the buffer region for beta variables (-mapcells - 1, mapcells)
    buffer_counts = (/lb_size*v_size*comm_size(2)*comm_size(3), lb_size*v_size*comm_size(1)*comm_size(3), &
                     & lb_size*v_size*comm_size(1)*comm_size(2)/)

    $:GPU_UPDATE(device='[v_size, comm_coords, comm_size]')

    buffer_count = buffer_counts(mpi_dir)
    boundary_conditions = (/bc_x, bc_y, bc_z/)
    beg_end = (/boundary_conditions(mpi_dir)%beg, boundary_conditions(mpi_dir)%end/)
    grid_dims = (/m, n, p/)

    if (pbc_loc == -1) then  ! PBC at the beginning
        ! Phase 1: Rightward accumulation Send END buffer to right neighbor, recv from left into BEG, ADD
        pack_offset = grid_dims(mpi_dir) + 1
        unpack_offset = 0
        dst_proc = merge(beg_end(2), MPI_PROC_NULL, beg_end(2) >= 0)
        src_proc = merge(beg_end(1), MPI_PROC_NULL, beg_end(1) >= 0)
        send_tag = 0
        recv_tag = 0
        replace_buff = .false.
    else
        ! Phase 2: Leftward distribution Send BEG buffer to left neighbor, recv from right into END, REPLACE
        pack_offset = 0
        unpack_offset = grid_dims(mpi_dir) + 1
        dst_proc = merge(beg_end(1), MPI_PROC_NULL, beg_end(1) >= 0)
        src_proc = merge(beg_end(2), MPI_PROC_NULL, beg_end(2) >= 0)
        send_tag = 1
        recv_tag = 1
        replace_buff = .true.
    end if

    ! Pack Buffer to Send
    #:for mpi_dir in [1, 2, 3]
        if (mpi_dir == ${mpi_dir}$) then
            #:if mpi_dir == 1
                $:GPU_PARALLEL_LOOP(private='[i, t, ix, iy, iz]')
                do r = 0, buffer_count - 1
                    i = mod(r, v_size) + 1
                    t = r/v_size
                    ix = mod(t, (lb_size)); t = t/(lb_size)
                    iy = mod(t, (comm_size(2))); iz = t/(comm_size(2))
                    buff_send(r) = real(q_comm(beta_vars(i))%sf(ix + (pack_offset - mapcells - 1), iy + (comm_coords(2)%beg), &
                              & iz + (comm_coords(3)%beg)), &
                              & kind=wp) - real(kahan_comp(beta_vars(i))%sf(ix + (pack_offset - mapcells - 1), &
                              & iy + (comm_coords(2)%beg), iz + (comm_coords(3)%beg)), kind=wp)
                end do
                $:END_GPU_PARALLEL_LOOP()
            #:elif mpi_dir == 2
                $:GPU_PARALLEL_LOOP(private='[i, t, ix, iy, iz]')
                do r = 0, buffer_count - 1
                    i = mod(r, v_size) + 1
                    t = r/v_size
                    ix = mod(t, (comm_size(1))); t = t/(comm_size(1))
                    iy = mod(t, (lb_size)); iz = t/(lb_size)
                    buff_send(r) = real(q_comm(beta_vars(i))%sf(ix + (comm_coords(1)%beg), iy + (pack_offset - mapcells - 1), &
                              & iz + (comm_coords(3)%beg)), &
                              & kind=wp) - real(kahan_comp(beta_vars(i))%sf(ix + (comm_coords(1)%beg), &
                              & iy + (pack_offset - mapcells - 1), iz + (comm_coords(3)%beg)), kind=wp)
                end do
                $:END_GPU_PARALLEL_LOOP()
            #:else
                $:GPU_PARALLEL_LOOP(private='[i, t, ix, iy, iz]')
                do r = 0, buffer_count - 1
                    i = mod(r, v_size) + 1
                    t = r/v_size
                    ix = mod(t, (comm_size(1))); t = t/(comm_size(1))
                    iy = mod(t, (comm_size(2))); iz = t/(comm_size(2))
                    buff_send(r) = real(q_comm(beta_vars(i))%sf(ix + (comm_coords(1)%beg), iy + (comm_coords(2)%beg), &
                              & iz + (pack_offset - mapcells - 1)), &
                              & kind=wp) - real(kahan_comp(beta_vars(i))%sf(ix + (comm_coords(1)%beg), iy + (comm_coords(2)%beg), &
                              & iz + (pack_offset - mapcells - 1)), kind=wp)
                end do
                $:END_GPU_PARALLEL_LOOP()
            #:endif
        end if
    #:endfor
    call nvtxEndRange  ! Packbuf

    ! MFC_COMM_BEGIN name=beta_axis mode=mpi_only pattern=sendrecv
    #:for rdma_mpi in [False, True]
        if (use_rdma_transport .eqv. ${'.true.' if rdma_mpi else '.false.'}$) then
            #:if rdma_mpi
                #:call GPU_HOST_DATA(use_device_addr='[buff_send, buff_recv]')
                    call nvtxStartRange("BETA-COMM-SENDRECV-RDMA")

                    call MPI_SENDRECV(buff_send, buffer_count, mpi_p, dst_proc, send_tag, buff_recv, buffer_count, mpi_p, &
                                      & src_proc, recv_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                    call nvtxEndRange  ! BETA-MPI-SENDRECV-(NO)-RDMA
                #:endcall GPU_HOST_DATA
                $:GPU_WAIT()
            #:else
                call nvtxStartRange("BETA-COMM-DEV2HOST")
                $:GPU_UPDATE(host='[buff_send]')
                call nvtxEndRange
                call nvtxStartRange("BETA-COMM-SENDRECV-NO-RMDA")

                call MPI_SENDRECV(buff_send, buffer_count, mpi_p, dst_proc, send_tag, buff_recv, buffer_count, mpi_p, src_proc, &
                                  & recv_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                call nvtxEndRange  ! BETA-MPI-SENDRECV-(NO)-RDMA

                call nvtxStartRange("BETA-COMM-HOST2DEV")
                $:GPU_UPDATE(device='[buff_recv]')
                call nvtxEndRange
            #:endif
        end if
    #:endfor
    ! MFC_COMM_END

    ! Unpack Received Buffer (skip if no source rank)
    call nvtxStartRange("BETA-COMM-UNPACKBUF")
    if (src_proc /= MPI_PROC_NULL) then
        #:for mpi_dir in [1, 2, 3]
            if (mpi_dir == ${mpi_dir}$) then
                #:if mpi_dir == 1
                    $:GPU_PARALLEL_LOOP(private='[i, t, ix, iy, iz, y_kahan, t_kahan]', copyin='[replace_buff]')
                    do r = 0, buffer_count - 1
                        i = mod(r, v_size) + 1
                        t = r/v_size
                        ix = mod(t, (lb_size)); t = t/(lb_size)
                        iy = mod(t, (comm_size(2))); iz = t/(comm_size(2))
                        if (replace_buff) then
                            q_comm(beta_vars(i))%sf(ix + (unpack_offset - mapcells - 1), iy + (comm_coords(2)%beg), &
                                   & iz + (comm_coords(3)%beg)) = real(buff_recv(r), kind=stp)
                            kahan_comp(beta_vars(i))%sf(ix + (unpack_offset - mapcells - 1), iy + (comm_coords(2)%beg), &
                                       & iz + (comm_coords(3)%beg)) = real(q_comm(beta_vars(i))%sf(ix + (unpack_offset - mapcells &
                                       & - 1), iy + (comm_coords(2)%beg), iz + (comm_coords(3)%beg)), kind=wp) - buff_recv(r)
                        else
                            y_kahan = buff_recv(r) - real(kahan_comp(beta_vars(i))%sf(ix + (unpack_offset - mapcells - 1), &
                                                & iy + (comm_coords(2)%beg), iz + (comm_coords(3)%beg)), kind=wp)
                            t_kahan = real(q_comm(beta_vars(i))%sf(ix + (unpack_offset - mapcells - 1), &
                                           & iy + (comm_coords(2)%beg), iz + (comm_coords(3)%beg)), kind=wp) + y_kahan
                            kahan_comp(beta_vars(i))%sf(ix + (unpack_offset - mapcells - 1), iy + (comm_coords(2)%beg), &
                                       & iz + (comm_coords(3)%beg)) = (t_kahan - q_comm(beta_vars(i))%sf(ix + (unpack_offset &
                                       & - mapcells - 1), iy + (comm_coords(2)%beg), iz + (comm_coords(3)%beg))) - y_kahan
                            q_comm(beta_vars(i))%sf(ix + (unpack_offset - mapcells - 1), iy + (comm_coords(2)%beg), &
                                   & iz + (comm_coords(3)%beg)) = t_kahan
                        end if
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                #:elif mpi_dir == 2
                    $:GPU_PARALLEL_LOOP(private='[i, t, ix, iy, iz, y_kahan, t_kahan]', copyin='[replace_buff]')
                    do r = 0, buffer_count - 1
                        i = mod(r, v_size) + 1
                        t = r/v_size
                        ix = mod(t, (comm_size(1))); t = t/(comm_size(1))
                        iy = mod(t, (lb_size)); iz = t/(lb_size)
                        if (replace_buff) then
                            q_comm(beta_vars(i))%sf(ix + (comm_coords(1)%beg), iy + (unpack_offset - mapcells - 1), &
                                   & iz + (comm_coords(3)%beg)) = real(buff_recv(r), kind=stp)
                            kahan_comp(beta_vars(i))%sf(ix + (comm_coords(1)%beg), iy + (unpack_offset - mapcells - 1), &
                                       & iz + (comm_coords(3)%beg)) = real(q_comm(beta_vars(i))%sf(ix + (comm_coords(1)%beg), &
                                       & iy + (unpack_offset - mapcells - 1), iz + (comm_coords(3)%beg)), kind=wp) - buff_recv(r)
                        else
                            y_kahan = buff_recv(r) - real(kahan_comp(beta_vars(i))%sf(ix + (comm_coords(1)%beg), &
                                                & iy + (unpack_offset - mapcells - 1), iz + (comm_coords(3)%beg)), kind=wp)
                            t_kahan = real(q_comm(beta_vars(i))%sf(ix + (comm_coords(1)%beg), &
                                           & iy + (unpack_offset - mapcells - 1), iz + (comm_coords(3)%beg)), kind=wp) + y_kahan
                            kahan_comp(beta_vars(i))%sf(ix + (comm_coords(1)%beg), iy + (unpack_offset - mapcells - 1), &
                                       & iz + (comm_coords(3)%beg)) = (t_kahan - q_comm(beta_vars(i))%sf(ix + (comm_coords(1)%beg) &
                                       & , iy + (unpack_offset - mapcells - 1), iz + (comm_coords(3)%beg))) - y_kahan
                            q_comm(beta_vars(i))%sf(ix + (comm_coords(1)%beg), iy + (unpack_offset - mapcells - 1), &
                                   & iz + (comm_coords(3)%beg)) = t_kahan
                        end if
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                #:else
                    $:GPU_PARALLEL_LOOP(private='[i, t, ix, iy, iz, y_kahan, t_kahan]', copyin='[replace_buff]')
                    do r = 0, buffer_count - 1
                        i = mod(r, v_size) + 1
                        t = r/v_size
                        ix = mod(t, (comm_size(1))); t = t/(comm_size(1))
                        iy = mod(t, (comm_size(2))); iz = t/(comm_size(2))
                        if (replace_buff) then
                            q_comm(beta_vars(i))%sf(ix + (comm_coords(1)%beg), iy + (comm_coords(2)%beg), &
                                   & iz + (unpack_offset - mapcells - 1)) = real(buff_recv(r), kind=stp)
                            kahan_comp(beta_vars(i))%sf(ix + (comm_coords(1)%beg), iy + (comm_coords(2)%beg), &
                                       & iz + (unpack_offset - mapcells - 1)) = real(q_comm(beta_vars(i))%sf(ix &
                                       & + (comm_coords(1)%beg), iy + (comm_coords(2)%beg), iz + (unpack_offset - mapcells - 1)), &
                                       & kind=wp) - buff_recv(r)
                        else
                            y_kahan = buff_recv(r) - real(kahan_comp(beta_vars(i))%sf(ix + (comm_coords(1)%beg), &
                                                & iy + (comm_coords(2)%beg), iz + (unpack_offset - mapcells - 1)), kind=wp)
                            t_kahan = real(q_comm(beta_vars(i))%sf(ix + (comm_coords(1)%beg), iy + (comm_coords(2)%beg), &
                                           & iz + (unpack_offset - mapcells - 1)), kind=wp) + y_kahan
                            kahan_comp(beta_vars(i))%sf(ix + (comm_coords(1)%beg), iy + (comm_coords(2)%beg), &
                                       & iz + (unpack_offset - mapcells - 1)) = (t_kahan - q_comm(beta_vars(i))%sf(ix &
                                       & + (comm_coords(1)%beg), iy + (comm_coords(2)%beg), &
                                       & iz + (unpack_offset - mapcells - 1))) - y_kahan
                            q_comm(beta_vars(i))%sf(ix + (comm_coords(1)%beg), iy + (comm_coords(2)%beg), &
                                   & iz + (unpack_offset - mapcells - 1)) = t_kahan
                        end if
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                #:endif
            end if
        #:endfor
    end if
    call nvtxEndRange
#endif

end subroutine s_mpi_reduce_beta_variables_buffers

!> The purpose of this procedure is to optimally decompose the computational domain among the available processors. This is
!! performed by attempting to award each processor, in each of the coordinate directions, approximately the same number of cells,
!! and then recomputing the affected global parameters.
subroutine s_mpi_decompose_computational_domain(write_silo_ghost_offsets, adjust_local_domains, output_offsets, local_domains)

    logical, intent(in)                                          :: write_silo_ghost_offsets
    logical, intent(in)                                          :: adjust_local_domains
    type(int_bounds_info), dimension(3), intent(inout), optional :: output_offsets
    type(bounds_info), dimension(3), intent(inout), optional     :: local_domains

#ifdef MFC_MPI
    !> Non-optimal number of processors in the x-, y- and z-directions
    real(wp)              :: tmp_num_procs_x, tmp_num_procs_y, tmp_num_procs_z
    real(wp)              :: fct_min        !< Processor factorization (fct) minimization parameter
    integer               :: MPI_COMM_CART  !< Cartesian processor topology communicator
    integer               :: rem_cells      !< Remaining cells after distribution among processors
    integer               :: rem_cells_by_dim(3)
    integer               :: recon_order    !< WENO or MUSCL reconstruction order
    integer               :: i, j, k, l     !< Generic loop iterators
    integer               :: ierr           !< Generic flag used to identify and report MPI errors
    integer, dimension(3) :: runtime_process_grid, cart_axes, physical_process_grid
    integer, dimension(3) :: cart_process_grid, cart_coords, cart_neighbor_coords
    logical               :: have_runtime_grid, have_runtime_axis_order
#ifdef MFC_SIMULATION
    character(len=128) :: runtime_grid_text, runtime_axis_order_text
    integer            :: runtime_grid_length, runtime_grid_status, runtime_axis_order_length, runtime_axis_order_status
    integer            :: runtime_read_status
#endif

    ! temp array to store neighbor rank coordinates
    integer, dimension(1:num_dims) :: neighbor_coords

    rem_cells_by_dim = 0

    ! Zeroing out communication needs for moving EL bubbles/particles
    nidx(1)%beg = 0; nidx(1)%end = 0
    nidx(2)%beg = 0; nidx(2)%end = 0
    nidx(3)%beg = 0; nidx(3)%end = 0

    if (recon_type == recon_type_weno) then
        recon_order = weno_order
    else
        recon_order = muscl_order
    end if

    if (num_procs == 1 .and. parallel_io) then
        do i = 1, num_dims
            start_idx(i) = 0
        end do
        return
    end if

    if (num_procs > 1) call s_sfc_read_weights()

    if (igr) then
        recon_order = igr_order
    end if

    runtime_process_grid = 1
    cart_axes = (/1, 2, 3/)
    have_runtime_grid = .false.
    have_runtime_axis_order = .false.
#ifdef MFC_SIMULATION
    call get_environment_variable('MFC_DECOMP_PROCESS_GRID', runtime_grid_text, length=runtime_grid_length, &
                                  & status=runtime_grid_status)
    if (runtime_grid_status == 0 .and. runtime_grid_length > 0) then
        read (runtime_grid_text(1:runtime_grid_length), *, iostat=runtime_read_status) runtime_process_grid
        if (runtime_read_status /= 0) call s_mpi_abort('Invalid runtime Cartesian process grid.')
        have_runtime_grid = .true.
    end if

    call get_environment_variable('MFC_DECOMP_AXIS_ORDER', runtime_axis_order_text, length=runtime_axis_order_length, &
                                  & status=runtime_axis_order_status)
    if (runtime_axis_order_status == 0 .and. runtime_axis_order_length > 0) then
        read (runtime_axis_order_text(1:runtime_axis_order_length), *, iostat=runtime_read_status) cart_axes
        if (runtime_read_status /= 0) call s_mpi_abort('Invalid runtime Cartesian axis order.')
        have_runtime_axis_order = .true.
    end if
#endif

    if (have_runtime_grid) then
        if (any(runtime_process_grid < 1) .or. product(runtime_process_grid) /= num_procs) then
            call s_mpi_abort('Runtime Cartesian process grid does not match the MPI rank count.')
        end if
        if ((num_dims < 3 .and. runtime_process_grid(3) /= 1) .or. (num_dims < 2 .and. runtime_process_grid(2) /= 1)) then
            call s_mpi_abort('Runtime Cartesian process grid decomposes an inactive direction.')
        end if
        if (fft_wrt .and. runtime_process_grid(1) /= 1) then
            call s_mpi_abort('Runtime Cartesian process grid is incompatible with fft_wrt.')
        end if
        if (cyl_coord .and. num_dims == 3 .and. runtime_process_grid(3) /= 1) then
            call s_mpi_abort('Runtime Cartesian process grid is incompatible with cylindrical pencil blocking.')
        end if
        if ((m + 1)/runtime_process_grid(1) < num_stcls_min*recon_order .or. (num_dims > 1 .and. (n + 1)/runtime_process_grid(2) &
            & < num_stcls_min*recon_order) .or. (num_dims > 2 .and. (p + 1)/runtime_process_grid(3) < num_stcls_min*recon_order)) &
            & then
            call s_mpi_abort('Runtime Cartesian process grid leaves too few cells for reconstruction.')
        end if
    end if

    if (have_runtime_axis_order) then
        do i = 1, 3
            if (cart_axes(i) < 1 .or. cart_axes(i) > 3) then
                call s_mpi_abort('Runtime Cartesian axis order is not a permutation.')
            end if
            do j = i + 1, 3
                if (cart_axes(i) == cart_axes(j)) call s_mpi_abort('Runtime Cartesian axis order is not a permutation.')
            end do
        end do
        if ((num_dims < 3 .and. cart_axes(3) /= 3) .or. (num_dims < 2 .and. cart_axes(2) /= 2)) then
            call s_mpi_abort('Runtime Cartesian axis order includes an inactive direction.')
        end if
    end if

    ! The other executables decompose with the default grid, so per-rank files would not match.
    if ((have_runtime_grid .or. have_runtime_axis_order) .and. .not. parallel_io) then
        call s_mpi_abort('MFC_DECOMP_PROCESS_GRID and MFC_DECOMP_AXIS_ORDER require parallel_io = T.')
    end if

    ! 3D Cartesian Processor Topology
    if (n > 0) then
        if (p > 0) then
            if (have_runtime_grid) then
                num_procs_x = runtime_process_grid(1)
                num_procs_y = runtime_process_grid(2)
                num_procs_z = runtime_process_grid(3)
                ierr = 0
            else if (fft_wrt) then
                ! Initial estimate of optimal processor topology
                num_procs_x = 1
                num_procs_y = 1
                num_procs_z = num_procs
                ierr = -1

                ! Benchmarking the quality of this initial guess
                tmp_num_procs_y = num_procs_y
                tmp_num_procs_z = num_procs_z
                fct_min = 10._wp*abs((n + 1)/tmp_num_procs_y - (p + 1)/tmp_num_procs_z)

                ! Optimization of the initial processor topology
                do i = 1, num_procs
                    if (mod(num_procs, i) == 0 .and. (n + 1)/i >= num_stcls_min*recon_order) then
                        tmp_num_procs_y = i
                        tmp_num_procs_z = num_procs/i

                        if (fct_min >= abs((n + 1)/tmp_num_procs_y - (p + 1)/tmp_num_procs_z) .and. (p + 1) &
                            & /tmp_num_procs_z >= num_stcls_min*recon_order) then
                            num_procs_y = i
                            num_procs_z = num_procs/i
                            fct_min = abs((n + 1)/tmp_num_procs_y - (p + 1)/tmp_num_procs_z)
                            ierr = 0
                        end if
                    end if
                end do
            else
                if (cyl_coord .and. p > 0) then
                    ! Pencil blocking for cylindrical coordinates (Fourier filter near axis)

                    ! Initial values of the processor factorization optimization
                    num_procs_x = 1
                    num_procs_y = num_procs
                    num_procs_z = 1
                    ierr = -1

                    ! Computing minimization variable for these initial values
                    tmp_num_procs_x = num_procs_x
                    tmp_num_procs_y = num_procs_y
                    tmp_num_procs_z = num_procs_z
                    fct_min = 10._wp*abs((m + 1)/tmp_num_procs_x - (n + 1)/tmp_num_procs_y)

                    ! Searching for optimal computational domain distribution
                    do i = 1, num_procs
                        if (mod(num_procs, i) == 0 .and. (m + 1)/i >= num_stcls_min*recon_order) then
                            tmp_num_procs_x = i
                            tmp_num_procs_y = num_procs/i

                            if (fct_min >= abs((m + 1)/tmp_num_procs_x - (n + 1)/tmp_num_procs_y) .and. (n + 1) &
                                & /tmp_num_procs_y >= num_stcls_min*recon_order) then
                                num_procs_x = i
                                num_procs_y = num_procs/i
                                fct_min = abs((m + 1)/tmp_num_procs_x - (n + 1)/tmp_num_procs_y)
                                ierr = 0
                            end if
                        end if
                    end do
                else
                    ! Initial estimate of optimal processor topology
                    num_procs_x = 1
                    num_procs_y = 1
                    num_procs_z = num_procs
                    ierr = -1

                    ! Benchmarking the quality of this initial guess
                    tmp_num_procs_x = num_procs_x
                    tmp_num_procs_y = num_procs_y
                    tmp_num_procs_z = num_procs_z
                    fct_min = 10._wp*abs((m + 1)/tmp_num_procs_x - (n + 1)/tmp_num_procs_y) + 10._wp*abs((n + 1)/tmp_num_procs_y &
                                         & - (p + 1)/tmp_num_procs_z)

                    ! Optimization of the initial processor topology
                    do i = 1, num_procs
                        if (mod(num_procs, i) == 0 .and. (m + 1)/i >= num_stcls_min*recon_order) then
                            do j = 1, num_procs/i
                                if (mod(num_procs/i, j) == 0 .and. (n + 1)/j >= num_stcls_min*recon_order) then
                                    tmp_num_procs_x = i
                                    tmp_num_procs_y = j
                                    tmp_num_procs_z = num_procs/(i*j)

                                    if (fct_min >= abs((m + 1)/tmp_num_procs_x - (n + 1)/tmp_num_procs_y) + abs((n + 1) &
                                        & /tmp_num_procs_y - (p + 1)/tmp_num_procs_z) .and. (p + 1) &
                                        & /tmp_num_procs_z >= num_stcls_min*recon_order) then
                                        num_procs_x = i
                                        num_procs_y = j
                                        num_procs_z = num_procs/(i*j)
                                        fct_min = abs((m + 1)/tmp_num_procs_x - (n + 1)/tmp_num_procs_y) + abs((n + 1) &
                                                      & /tmp_num_procs_y - (p + 1)/tmp_num_procs_z)
                                        ierr = 0
                                    end if
                                end if
                            end do
                        end if
                    end do
                    if (sfc_active) call s_sfc_select_topology(recon_order)
                end if
            end if

            ! Verifying that a valid decomposition of the computational domain has been established. If not, the simulation exits.
            if (proc_rank == 0 .and. ierr == -1) then
                call s_mpi_abort('Unsupported combination of values ' // 'of num_procs, m, n, p and ' &
                                 & // 'weno/muscl/igr_order. Exiting.')
            end if

            ! Creating new communicator using the Cartesian topology, its axes in cart_axes order
            physical_process_grid = (/num_procs_x, num_procs_y, num_procs_z/)
            do i = 1, 3
                cart_process_grid(i) = physical_process_grid(cart_axes(i))
            end do
            call MPI_CART_CREATE(MPI_COMM_WORLD, 3, cart_process_grid, (/.true., .true., .true./), .false., MPI_COMM_CART, ierr)

            ! Finding the Cartesian coordinates of the local process
            call MPI_CART_COORDS(MPI_COMM_CART, proc_rank, 3, cart_coords, ierr)
            do i = 1, 3
                proc_coords(cart_axes(i)) = cart_coords(i)
            end do

            ! Global Parameters for z-direction

            ! Number of remaining cells
            rem_cells = mod(p + 1, num_procs_z)
            rem_cells_by_dim(3) = rem_cells

            ! Optimal number of cells per processor
            p = (p + 1)/num_procs_z - 1

            ! Distributing the remaining cells
            do i = 1, rem_cells
                if (proc_coords(3) == i - 1) then
                    p = p + 1; exit
                end if
            end do

            ! Boundary condition at the beginning
            if (proc_coords(3) > 0 .or. (bc_z%beg == BC_PERIODIC .and. num_procs_z > 1)) then
                proc_coords(3) = proc_coords(3) - 1
                do l = 1, num_dims
                    cart_coords(l) = proc_coords(cart_axes(l))
                end do
                call MPI_CART_RANK(MPI_COMM_CART, cart_coords, bc_z%beg, ierr)
                proc_coords(3) = proc_coords(3) + 1
                nidx(3)%beg = -1
            end if

            ! Boundary condition at the end
            if (proc_coords(3) < num_procs_z - 1 .or. (bc_z%end == BC_PERIODIC .and. num_procs_z > 1)) then
                proc_coords(3) = proc_coords(3) + 1
                do l = 1, num_dims
                    cart_coords(l) = proc_coords(cart_axes(l))
                end do
                call MPI_CART_RANK(MPI_COMM_CART, cart_coords, bc_z%end, ierr)
                proc_coords(3) = proc_coords(3) - 1
                nidx(3)%end = 1
            end if

            ! Beginning and end sub-domain boundary locations
            if (parallel_io) then
                if (proc_coords(3) < rem_cells) then
                    start_idx(3) = (p + 1)*proc_coords(3)
                else
                    start_idx(3) = (p + 1)*proc_coords(3) + rem_cells
                end if
            end if

            ! 2D Cartesian Processor Topology
        else
            if (have_runtime_grid) then
                num_procs_x = runtime_process_grid(1)
                num_procs_y = runtime_process_grid(2)
                ierr = 0
            else
                ! Initial estimate of optimal processor topology
                num_procs_x = 1
                num_procs_y = num_procs
                ierr = -1

                ! Benchmarking the quality of this initial guess
                tmp_num_procs_x = num_procs_x
                tmp_num_procs_y = num_procs_y
                fct_min = 10._wp*abs((m + 1)/tmp_num_procs_x - (n + 1)/tmp_num_procs_y)

                ! Optimization of the initial processor topology
                do i = 1, num_procs
                    if (mod(num_procs, i) == 0 .and. (m + 1)/i >= num_stcls_min*recon_order) then
                        tmp_num_procs_x = i
                        tmp_num_procs_y = num_procs/i

                        if (fct_min >= abs((m + 1)/tmp_num_procs_x - (n + 1)/tmp_num_procs_y) .and. (n + 1) &
                            & /tmp_num_procs_y >= num_stcls_min*recon_order) then
                            num_procs_x = i
                            num_procs_y = num_procs/i
                            fct_min = abs((m + 1)/tmp_num_procs_x - (n + 1)/tmp_num_procs_y)
                            ierr = 0
                        end if
                    end if
                end do
                if (sfc_active) call s_sfc_select_topology(recon_order)
            end if

            ! Verifying that a valid decomposition of the computational domain has been established. If not, the simulation exits.
            if (proc_rank == 0 .and. ierr == -1) then
                call s_mpi_abort('Unsupported combination of values ' // 'of num_procs, m, n and ' &
                                 & // 'weno/muscl/igr_order. Exiting.')
            end if

            ! Creating new communicator using the Cartesian topology
            physical_process_grid = (/num_procs_x, num_procs_y, 1/)
            do i = 1, 2
                cart_process_grid(i) = physical_process_grid(cart_axes(i))
            end do
            call MPI_CART_CREATE(MPI_COMM_WORLD, 2, cart_process_grid(1:2), (/.true., .true./), .false., MPI_COMM_CART, ierr)

            ! Finding the Cartesian coordinates of the local process
            call MPI_CART_COORDS(MPI_COMM_CART, proc_rank, 2, cart_coords, ierr)
            do i = 1, 2
                proc_coords(cart_axes(i)) = cart_coords(i)
            end do
        end if

        ! Global Parameters for y-direction

        ! Number of remaining cells
        rem_cells = mod(n + 1, num_procs_y)
        rem_cells_by_dim(2) = rem_cells

        ! Optimal number of cells per processor
        n = (n + 1)/num_procs_y - 1

        ! Distributing the remaining cells
        do i = 1, rem_cells
            if (proc_coords(2) == i - 1) then
                n = n + 1; exit
            end if
        end do

        ! Boundary condition at the beginning
        if (proc_coords(2) > 0 .or. (bc_y%beg == BC_PERIODIC .and. num_procs_y > 1)) then
            proc_coords(2) = proc_coords(2) - 1
            do l = 1, num_dims
                cart_coords(l) = proc_coords(cart_axes(l))
            end do
            call MPI_CART_RANK(MPI_COMM_CART, cart_coords, bc_y%beg, ierr)
            proc_coords(2) = proc_coords(2) + 1
            nidx(2)%beg = -1
        end if

        ! Boundary condition at the end
        if (proc_coords(2) < num_procs_y - 1 .or. (bc_y%end == BC_PERIODIC .and. num_procs_y > 1)) then
            proc_coords(2) = proc_coords(2) + 1
            do l = 1, num_dims
                cart_coords(l) = proc_coords(cart_axes(l))
            end do
            call MPI_CART_RANK(MPI_COMM_CART, cart_coords, bc_y%end, ierr)
            proc_coords(2) = proc_coords(2) - 1
            nidx(2)%end = 1
        end if

        ! Beginning and end sub-domain boundary locations
        if (parallel_io) then
            if (proc_coords(2) < rem_cells) then
                start_idx(2) = (n + 1)*proc_coords(2)
            else
                start_idx(2) = (n + 1)*proc_coords(2) + rem_cells
            end if
        end if

        ! 1D Cartesian Processor Topology
    else
        ! Optimal processor topology
        num_procs_x = num_procs
        if (have_runtime_grid) num_procs_x = runtime_process_grid(1)

        ! Creating new communicator using the Cartesian topology
        call MPI_CART_CREATE(MPI_COMM_WORLD, 1, (/num_procs_x/), (/.true./), .false., MPI_COMM_CART, ierr)

        ! Finding the Cartesian coordinates of the local process
        call MPI_CART_COORDS(MPI_COMM_CART, proc_rank, 1, proc_coords, ierr)
    end if

    ! Global Parameters for x-direction

    ! Number of remaining cells
    rem_cells = mod(m + 1, num_procs_x)
    rem_cells_by_dim(1) = rem_cells

    ! Optimal number of cells per processor
    m = (m + 1)/num_procs_x - 1

    ! Distributing the remaining cells
    do i = 1, rem_cells
        if (proc_coords(1) == i - 1) then
            m = m + 1; exit
        end if
    end do

    call s_update_cell_bounds(cells_bounds, m, n, p)

    ! Boundary condition at the beginning
    if (proc_coords(1) > 0 .or. (bc_x%beg == BC_PERIODIC .and. num_procs_x > 1)) then
        proc_coords(1) = proc_coords(1) - 1
        do l = 1, num_dims
            cart_coords(l) = proc_coords(cart_axes(l))
        end do
        call MPI_CART_RANK(MPI_COMM_CART, cart_coords, bc_x%beg, ierr)
        proc_coords(1) = proc_coords(1) + 1
        nidx(1)%beg = -1
    end if

    ! Boundary condition at the end
    if (proc_coords(1) < num_procs_x - 1 .or. (bc_x%end == BC_PERIODIC .and. num_procs_x > 1)) then
        proc_coords(1) = proc_coords(1) + 1
        do l = 1, num_dims
            cart_coords(l) = proc_coords(cart_axes(l))
        end do
        call MPI_CART_RANK(MPI_COMM_CART, cart_coords, bc_x%end, ierr)
        proc_coords(1) = proc_coords(1) - 1
        nidx(1)%end = 1
    end if

    ! Beginning and end sub-domain boundary locations
    if (parallel_io) then
        if (proc_coords(1) < rem_cells) then
            start_idx(1) = (m + 1)*proc_coords(1)
        else
            start_idx(1) = (m + 1)*proc_coords(1) + rem_cells
        end if
    end if

    if (sfc_active) then
        call s_sfc_apply_weighted_partition(recon_order, adjust_local_domains .and. (.not. parallel_io), local_domains)
        call s_apply_decomposition_policies((/num_procs_x, num_procs_y, num_procs_z/), rem_cells_by_dim, (/m, n, p/), (/m_glb, &
                                            & n_glb, p_glb/), write_silo_ghost_offsets, .false., output_offsets, local_domains)
    else
        call s_apply_decomposition_policies((/num_procs_x, num_procs_y, num_procs_z/), rem_cells_by_dim, (/m, n, p/), (/m_glb, &
                                            & n_glb, p_glb/), write_silo_ghost_offsets, &
                                            & adjust_local_domains .and. (.not. parallel_io), output_offsets, local_domains)
    end if

    @:ALLOCATE(neighbor_ranks(nidx(1)%beg:nidx(1)%end, nidx(2)%beg:nidx(2)%end, nidx(3)%beg:nidx(3)%end))
    do k = nidx(3)%beg, nidx(3)%end
        do j = nidx(2)%beg, nidx(2)%end
            do i = nidx(1)%beg, nidx(1)%end
                if (abs(i) + abs(j) + abs(k) > 0) then
                    neighbor_coords(1) = proc_coords(1) + i
                    if (num_dims > 1) neighbor_coords(2) = proc_coords(2) + j
                    if (num_dims > 2) neighbor_coords(3) = proc_coords(3) + k
                    do l = 1, num_dims
                        cart_neighbor_coords(l) = neighbor_coords(cart_axes(l))
                    end do
                    call MPI_CART_RANK(MPI_COMM_CART, cart_neighbor_coords, neighbor_ranks(i, j, k), ierr)
                end if
            end do
        end do
    end do
#endif

end subroutine s_mpi_decompose_computational_domain

#ifdef MFC_MPI
!> Read the env-gated per-block cost-weights file (MFC_SFC_WEIGHTS). The header line is "g1 g2 g3 block snap1 snap2 snap3" (global
!! cells per direction, block edge length in cells, per-direction cut snapping stride with 1 = none), followed by
!! sfc_nblk(1)*sfc_nblk(2)*sfc_nblk(3) block weight sums, x fastest. Weights are relative cost units; a uniform file reproduces the
!! equal-cell decomposition exactly.
impure subroutine s_sfc_read_weights()

    character(len=path_len) :: wpath
    integer                 :: iunit, ierr, i, nblocks

    call get_environment_variable('MFC_SFC_WEIGHTS', value=wpath, status=ierr)
    if (ierr /= 0) return
    sfc_active = .true.

    open (newunit=iunit, file=trim(wpath), status='old', action='read', iostat=ierr)
    if (ierr /= 0) call s_mpi_abort('MFC_SFC_WEIGHTS: cannot open ' // trim(wpath))
    read (iunit, *, iostat=ierr) sfc_gcells(1), sfc_gcells(2), sfc_gcells(3), sfc_blksz, sfc_snap(1), sfc_snap(2), sfc_snap(3)
    if (ierr /= 0) call s_mpi_abort('MFC_SFC_WEIGHTS: malformed header.')
    if (sfc_blksz < 1 .or. any(sfc_snap < 1)) call s_mpi_abort('MFC_SFC_WEIGHTS: invalid block or snap size.')
    if (any(sfc_gcells /= (/m_glb + 1, n_glb + 1, p_glb + 1/))) then
        call s_mpi_abort('MFC_SFC_WEIGHTS: global cell counts do not match the case.')
    end if
    sfc_nblk = (sfc_gcells + sfc_blksz - 1)/sfc_blksz
    nblocks = product(sfc_nblk)

    @:ALLOCATE(sfc_wblk(0:nblocks - 1))
    read (iunit, *, iostat=ierr) (sfc_wblk(i), i=0, nblocks - 1)
    close (iunit)
    if (ierr /= 0) then
        @:DEALLOCATE(sfc_wblk)
        call s_mpi_abort('MFC_SFC_WEIGHTS: truncated weights list.')
    end if

end subroutine s_sfc_read_weights

!> Cumulative offsets splitting the weight marginal w into n_parts contiguous chunks of near-equal weight, each at least l_min
!! entries. off(0) = 0, off(n_parts) = size(w). A degenerate marginal (total <= 0) falls back to the equal split.
pure function f_sfc_weighted_offsets(w, n_parts, l_min) result(off)

    real(wp), dimension(0:), intent(in) :: w
    integer, intent(in)                 :: n_parts, l_min
    integer, dimension(0:n_parts)       :: off
    real(wp)                            :: csum, total
    integer                             :: g, i, r, q, rem

    g = size(w)
    q = g/n_parts
    rem = mod(g, n_parts)
    do r = 0, n_parts
        off(r) = r*q + min(r, rem)
    end do
    if (n_parts == 1) return
    total = sum(w)
    if (total <= 0._wp) return
    r = 1
    csum = 0._wp
    do i = 0, g - 1
        csum = csum + w(i)
        do while (r < n_parts .and. csum >= real(r, wp)*total/real(n_parts, wp))
            off(r) = i + 1
            r = r + 1
        end do
    end do
    do while (r < n_parts)
        off(r) = g
        r = r + 1
    end do
    do r = 1, n_parts - 1
        if (off(r) < off(r - 1) + l_min) off(r) = off(r - 1) + l_min
        if (off(r) > g - (n_parts - r)*l_min) off(r) = g - (n_parts - r)*l_min
    end do

end function f_sfc_weighted_offsets

!> Block-marginal weights along direction d (sum of sfc_wblk over the two other directions), spread uniformly over the cells of each
!! block.
pure subroutine s_sfc_block_marginal(d, wm)

    integer, intent(in)                  :: d
    real(wp), dimension(0:), intent(out) :: wm
    real(wp)                             :: wb(0:sfc_nblk(d) - 1)
    integer                              :: bx, by, bz, i, nb, cells_in_blk

    wb = 0._wp
    do bz = 0, sfc_nblk(3) - 1
        do by = 0, sfc_nblk(2) - 1
            do bx = 0, sfc_nblk(1) - 1
                if (d == 1) then
                    wb(bx) = wb(bx) + sfc_wblk(bx + sfc_nblk(1)*(by + sfc_nblk(2)*bz))
                else if (d == 2) then
                    wb(by) = wb(by) + sfc_wblk(bx + sfc_nblk(1)*(by + sfc_nblk(2)*bz))
                else
                    wb(bz) = wb(bz) + sfc_wblk(bx + sfc_nblk(1)*(by + sfc_nblk(2)*bz))
                end if
            end do
        end do
    end do

    do i = 0, sfc_gcells(d) - 1
        nb = i/sfc_blksz
        cells_in_blk = min(sfc_blksz, sfc_gcells(d) - nb*sfc_blksz)
        wm(i) = wb(nb)/real(cells_in_blk, wp)
    end do

end subroutine s_sfc_block_marginal

!> Cumulative cell cut offsets along direction d for num_parts ranks, weight-balanced from the block marginal, optionally snapped to
!! sfc_snap(d)-cell boundaries, with every piece at least l_min cells. Every rank computes the same offsets.
pure function f_sfc_axis_offsets(d, num_parts, l_min) result(off)

    integer, intent(in)             :: d, num_parts, l_min
    integer, dimension(0:num_parts) :: off
    real(wp)                        :: wm(0:sfc_gcells(d) - 1)
    integer                         :: r

    call s_sfc_block_marginal(d, wm)
    off = f_sfc_weighted_offsets(wm, num_parts, l_min)
    if (sfc_snap(d) > 1) then
        do r = 1, num_parts - 1
            off(r) = sfc_snap(d)*(off(r)/sfc_snap(d))
        end do
        do r = 1, num_parts - 1
            if (off(r) < off(r - 1) + l_min) off(r) = off(r - 1) + l_min
            if (off(r) > sfc_gcells(d) - (num_parts - r)*l_min) off(r) = sfc_gcells(d) - (num_parts - r)*l_min
        end do
    end if

end function f_sfc_axis_offsets

!> Index of the partition (0-based) containing cell position pos, from cumulative cut offsets, by linear scan over the few parts.
pure function f_sfc_owner(off, pos) result(r)

    integer, dimension(0:), intent(in) :: off
    integer, intent(in)                :: pos
    integer                            :: r

    r = 0
    do while (r < ubound(off, 1) .and. off(r + 1) <= pos)
        r = r + 1
    end do

end function f_sfc_owner

!> Predicted per-rank weight shares of the weighted partition for the (complete) process grid pd; the return value is the maximum
!! share over the mean.
pure function f_sfc_partition_imbalance(pd, l_min) result(imb)

    integer, dimension(3), intent(in) :: pd
    integer, intent(in)               :: l_min
    real(wp)                          :: imb
    integer                           :: offx(0:pd(1)), offy(0:pd(2)), offz(0:pd(3))
    real(wp)                          :: wsum(0:product(pd) - 1)
    real(wp)                          :: wtot
    integer                           :: bx, by, bz, r

    offx = f_sfc_axis_offsets(1, pd(1), l_min)
    offy = f_sfc_axis_offsets(2, pd(2), l_min)
    offz = f_sfc_axis_offsets(3, pd(3), l_min)
    wsum = 0._wp
    do bz = 0, sfc_nblk(3) - 1
        do by = 0, sfc_nblk(2) - 1
            do bx = 0, sfc_nblk(1) - 1
                r = f_sfc_owner(offx, bx*sfc_blksz) + pd(1)*(f_sfc_owner(offy, by*sfc_blksz) + pd(2)*f_sfc_owner(offz, &
                                & bz*sfc_blksz))
                wsum(r) = wsum(r) + sfc_wblk(bx + sfc_nblk(1)*(by + sfc_nblk(2)*bz))
            end do
        end do
    end do
    wtot = sum(wsum)
    if (wtot <= 0._wp) then
        imb = 1._wp
    else
        imb = maxval(wsum)/wtot*real(product(pd), wp)
    end if

end function f_sfc_partition_imbalance

!> Replace the default process grid with the factorization of num_procs whose weighted partition minimizes the predicted per-rank
!! weight imbalance, over the grids leaving every direction at least num_stcls_min*recon_order cells per rank. Ties keep the grid
!! seen first, so a uniform weights file keeps the default decomposition.
impure subroutine s_sfc_select_topology(recon_order)

    integer, intent(in) :: recon_order
    integer             :: l_min, pd(3), cand(3), i, j
    real(wp)            :: imb, imb_best

    l_min = num_stcls_min*recon_order
    pd = 1
    pd(1) = num_procs_x
    if (num_dims >= 2) pd(2) = num_procs_y
    if (num_dims >= 3) pd(3) = num_procs_z
    imb_best = f_sfc_partition_imbalance(pd, l_min)
    cand = pd

    if (num_dims == 3) then
        do i = 1, num_procs
            if (mod(num_procs, i) /= 0 .or. sfc_gcells(1)/i < l_min) cycle
            do j = 1, num_procs/i
                if (mod(num_procs/i, j) /= 0 .or. sfc_gcells(2)/j < l_min) cycle
                if (sfc_gcells(3)/(num_procs/(i*j)) < l_min) cycle
                pd = (/i, j, num_procs/(i*j)/)
                imb = f_sfc_partition_imbalance(pd, l_min)
                if (imb < imb_best) then
                    imb_best = imb
                    cand = pd
                end if
            end do
        end do
    else
        do i = 1, num_procs
            if (mod(num_procs, i) /= 0 .or. sfc_gcells(1)/i < l_min) cycle
            if (sfc_gcells(2)/(num_procs/i) < l_min) cycle
            pd = (/i, num_procs/i, 1/)
            imb = f_sfc_partition_imbalance(pd, l_min)
            if (imb < imb_best) then
                imb_best = imb
                cand = pd
            end if
        end do
    end if

    num_procs_x = cand(1)
    num_procs_y = cand(2)
    num_procs_z = cand(3)

end subroutine s_sfc_select_topology

!> Override the equal-cell per-axis splits with weight-balanced splits (identically computed on every rank), optionally moving the
!! local domain bounds to the weighted cuts.
impure subroutine s_sfc_apply_weighted_partition(recon_order, adjust_local_domains, local_domains)

    integer, intent(in)                                      :: recon_order
    logical, intent(in)                                      :: adjust_local_domains
    type(bounds_info), dimension(3), intent(inout), optional :: local_domains
    real(wp)                                                 :: spacing
    integer                                                  :: l_min, offx(0:num_procs_x), offy(0:num_procs_y), offz(0:num_procs_z)

    l_min = num_stcls_min*recon_order
    offx = f_sfc_axis_offsets(1, num_procs_x, l_min)
    offy = f_sfc_axis_offsets(2, num_procs_y, l_min)
    offz = f_sfc_axis_offsets(3, num_procs_z, l_min)

    m = offx(proc_coords(1) + 1) - offx(proc_coords(1)) - 1
    if (num_dims >= 2) n = offy(proc_coords(2) + 1) - offy(proc_coords(2)) - 1
    if (num_dims >= 3) p = offz(proc_coords(3) + 1) - offz(proc_coords(3)) - 1
    if (parallel_io) then
        start_idx(1) = offx(proc_coords(1))
        if (num_dims >= 2) start_idx(2) = offy(proc_coords(2))
        if (num_dims >= 3) start_idx(3) = offz(proc_coords(3))
    end if
    call s_update_cell_bounds(cells_bounds, m, n, p)
    if (adjust_local_domains) then
        spacing = (local_domains(1)%end - local_domains(1)%beg)/real(sfc_gcells(1), wp)
        dx_min = spacing
        local_domains(1)%end = local_domains(1)%end - spacing*real(sfc_gcells(1) - offx(proc_coords(1) + 1), wp)
        local_domains(1)%beg = local_domains(1)%beg + spacing*real(offx(proc_coords(1)), wp)
        if (num_dims >= 2) then
            spacing = (local_domains(2)%end - local_domains(2)%beg)/real(sfc_gcells(2), wp)
            dy_min = spacing
            local_domains(2)%end = local_domains(2)%end - spacing*real(sfc_gcells(2) - offy(proc_coords(2) + 1), wp)
            local_domains(2)%beg = local_domains(2)%beg + spacing*real(offy(proc_coords(2)), wp)
        end if
        if (num_dims >= 3) then
            spacing = (local_domains(3)%end - local_domains(3)%beg)/real(sfc_gcells(3), wp)
            dz_min = spacing
            local_domains(3)%end = local_domains(3)%end - spacing*real(sfc_gcells(3) - offz(proc_coords(3) + 1), wp)
            local_domains(3)%beg = local_domains(3)%beg + spacing*real(offz(proc_coords(3)), wp)
        end if
    end if

end subroutine s_sfc_apply_weighted_partition
#endif

!> Apply executable-configured output and local-domain policies after the shared Cartesian decomposition.
subroutine s_apply_decomposition_policies(proc_counts, remainders, local_cells, global_cells, write_silo_ghost_offsets, &
    & adjust_local_domains, output_offsets, local_domains)

    integer, dimension(3), intent(in)                            :: proc_counts, remainders, local_cells, global_cells
    logical, intent(in)                                          :: write_silo_ghost_offsets, adjust_local_domains
    type(int_bounds_info), dimension(3), intent(inout), optional :: output_offsets
    type(bounds_info), dimension(3), intent(inout), optional     :: local_domains
    integer                                                      :: dim
    real(wp)                                                     :: domain_beg, domain_end, spacing

    if (write_silo_ghost_offsets .and. .not. present(output_offsets)) then
        call s_mpi_abort('Silo ghost-offset policy requires output offset storage.')
    end if
    if (adjust_local_domains .and. .not. present(local_domains)) then
        call s_mpi_abort('Local-domain adjustment policy requires domain storage.')
    end if

    if (present(output_offsets)) then
        do dim = 1, 3
            output_offsets(dim)%beg = 0
            output_offsets(dim)%end = 0
        end do

        if (write_silo_ghost_offsets) then
            do dim = 1, num_dims
                if (proc_coords(dim) > 0) output_offsets(dim)%beg = 2
                if (proc_coords(dim) < proc_counts(dim) - 1) output_offsets(dim)%end = 2
            end do
        end if
    end if

    if (adjust_local_domains) then
        do dim = 1, num_dims
            domain_beg = local_domains(dim)%beg
            domain_end = local_domains(dim)%end
            spacing = (domain_end - domain_beg)/real(global_cells(dim) + 1, wp)

            select case (dim)
            case (1)
                dx_min = spacing
            case (2)
                dy_min = spacing
            case (3)
                dz_min = spacing
            end select

            if (proc_coords(dim) < remainders(dim)) then
                local_domains(dim)%beg = domain_beg + spacing*real((local_cells(dim) + 1)*proc_coords(dim), wp)
                local_domains(dim)%end = domain_end - spacing*real((local_cells(dim) + 1)*(proc_counts(dim) - proc_coords(dim) &
                              & - 1) - (proc_counts(dim) - remainders(dim)), wp)
            else
                local_domains(dim)%beg = domain_beg + spacing*real((local_cells(dim) + 1)*proc_coords(dim) + remainders(dim), wp)
                local_domains(dim)%end = domain_end - spacing*real((local_cells(dim) + 1)*(proc_counts(dim) - proc_coords(dim) &
                              & - 1), wp)
            end if
        end do
    end if

end subroutine s_apply_decomposition_policies

!> The goal of this procedure is to populate the buffers of the grid variables by communicating with the neighboring processors.
!! Note that only the buffers of the cell-width distributions are handled in such a way. This is because the buffers of
!! cell-boundary locations may be calculated directly from those of the cell-width distributions.
subroutine s_mpi_sendrecv_grid_variable_buffer(cell_boundaries, cell_centers, cell_widths, num_cells, bc_bounds, pbc_loc, offset)

    integer, intent(in)               :: num_cells, pbc_loc
    type(int_bounds_info), intent(in) :: bc_bounds, offset
    ! Contiguous so that passing an element to MPI is a plain address, with no descriptor
    ! or copy-in/copy-out. Every actual argument is a whole module array. The attribute must
    ! be on ALL THREE, and on every frame that forwards them: CCE 19 IPA drops the stores to
    ! the non-contiguous dummies when one call mixes contiguous and non-contiguous arrays.
    real(wp), contiguous, intent(inout) :: cell_boundaries(-1 - offset%beg:)
    real(wp), contiguous, intent(inout) :: cell_centers(-buff_size:)
    real(wp), contiguous, intent(inout) :: cell_widths(-buff_size:)

#ifdef MFC_MPI
    integer :: ierr
    integer :: i

    ! MFC_COMM_BEGIN name=grid_setup mode=mpi_only pattern=sendrecv
    if (pbc_loc == -1) then
        if (bc_bounds%end >= 0) then
            call MPI_SENDRECV(cell_widths(num_cells - buff_size + 1), buff_size, mpi_p, bc_bounds%end, 0, &
                              & cell_widths(-buff_size), buff_size, mpi_p, bc_bounds%beg, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                              & ierr)
        else
            call MPI_SENDRECV(cell_widths(0), buff_size, mpi_p, bc_bounds%beg, 1, cell_widths(-buff_size), buff_size, mpi_p, &
                              & bc_bounds%beg, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        end if
        do i = 1, offset%beg
            cell_boundaries(-1 - i) = cell_boundaries(-i) - cell_widths(-i)
        end do
        do i = 1, buff_size
            cell_centers(-i) = cell_centers(1 - i) - (cell_widths(1 - i) + cell_widths(-i))/2._wp
        end do
    else
        if (bc_bounds%beg >= 0) then
            call MPI_SENDRECV(cell_widths(0), buff_size, mpi_p, bc_bounds%beg, 1, cell_widths(num_cells + 1), buff_size, mpi_p, &
                              & bc_bounds%end, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        else
            call MPI_SENDRECV(cell_widths(num_cells - buff_size + 1), buff_size, mpi_p, bc_bounds%end, 0, &
                              & cell_widths(num_cells + 1), buff_size, mpi_p, bc_bounds%end, 1, MPI_COMM_WORLD, &
                              & MPI_STATUS_IGNORE, ierr)
        end if
        do i = 1, offset%end
            cell_boundaries(num_cells + i) = cell_boundaries(num_cells + i - 1) + cell_widths(num_cells + i)
        end do
        do i = 1, buff_size
            cell_centers(num_cells + i) = cell_centers(num_cells + i - 1) + (cell_widths(num_cells + i - 1) &
                         & + cell_widths(num_cells + i))/2._wp
        end do
    end if
    ! MFC_COMM_END
#endif

end subroutine s_mpi_sendrecv_grid_variable_buffer

!> Populate the local cell-boundary, cell-center, and cell-width arrays in one direction directly from the global cell-boundary
!! array. This guarantees that every rank sees bitwise-identical values at any shared physical cell or boundary
subroutine s_apply_grid_from_global_dim(x_cb_glb, m_dim_glb, m_dim, sidx, bc_beg, bc_end, cb_lo, cb_hi, cw_lo, cw_hi, x_cb_loc, &
                                        & x_cc_loc, dx_loc)

    integer, intent(in)     :: m_dim_glb, m_dim, sidx, bc_beg, bc_end
    integer, intent(in)     :: cb_lo, cb_hi, cw_lo, cw_hi
    real(wp), intent(in)    :: x_cb_glb(-1:m_dim_glb)
    real(wp), intent(inout) :: x_cb_loc(-1 - cb_lo:m_dim + cb_hi)
    real(wp), intent(inout) :: x_cc_loc(-cw_lo:m_dim + cw_hi)
    real(wp), intent(inout) :: dx_loc(-cw_lo:m_dim + cw_hi)
    real(wp)                :: domain_len
    integer                 :: i, gidx, lo, hi

    domain_len = x_cb_glb(m_dim_glb) - x_cb_glb(-1)

    ! Interior cell boundaries sliced directly from the global list
    do i = -1, m_dim
        x_cb_loc(i) = x_cb_glb(sidx + i)
    end do

    ! Left ghost cell boundaries
    if (bc_beg >= 0) then
        if (sidx == 0) then
            ! Leftmost rank with a neighbor -> periodic+multirank, so wrap from the global right end
            do i = 1, cb_lo
                x_cb_loc(-1 - i) = x_cb_glb(m_dim_glb - i) - domain_len
            end do
        else
            do i = 1, cb_lo
                gidx = sidx - 1 - i
                if (gidx >= -1) then
                    x_cb_loc(-1 - i) = x_cb_glb(gidx)
                else
                    x_cb_loc(-1 - i) = x_cb_glb(m_dim_glb + 1 + gidx) - domain_len
                end if
            end do
        end if
    end if

    ! Right ghost cell boundaries
    if (bc_end >= 0) then
        if (sidx + m_dim == m_dim_glb) then
            ! Rightmost rank with a neighbor -> periodic+multirank, wrap from the global left end
            do i = 1, cb_hi
                x_cb_loc(m_dim + i) = x_cb_glb(i - 1) + domain_len
            end do
        else
            do i = 1, cb_hi
                gidx = sidx + m_dim + i
                if (gidx <= m_dim_glb) then
                    x_cb_loc(m_dim + i) = x_cb_glb(gidx)
                else
                    x_cb_loc(m_dim + i) = x_cb_glb(gidx - m_dim_glb - 1) + domain_len
                end if
            end do
        end if
    end if

    ! Recompute dx and x_cc over the range where x_cb is now valid using one formula so values are bitwise-identical
    if (bc_beg >= 0) then
        lo = -min(cw_lo, cb_lo)
    else
        lo = 0
    end if

    if (bc_end >= 0) then
        hi = m_dim + min(cw_hi, cb_hi)
    else
        hi = m_dim
    end if

    do i = lo, hi
        dx_loc(i) = x_cb_loc(i) - x_cb_loc(i - 1)
        x_cc_loc(i) = (x_cb_loc(i) + x_cb_loc(i - 1))/2._wp
    end do

end subroutine s_apply_grid_from_global_dim

!> Module deallocation and/or disassociation procedures
impure subroutine s_finalize_mpi_common_module

#ifdef MFC_MPI
    if (allocated(sfc_wblk)) then
        @:DEALLOCATE(sfc_wblk)
    end if
#ifdef MFC_HALO_SINGLE_ROUND
#if defined(MFC_SIMULATION) && !(defined(MFC_NVSHMEM) && defined(__NVCOMPILER))
        ! Retire the last exchange's outstanding requests and free the persistent sets.
        call s_halo_26_rdma_complete()
    ! MFC_COMM_BEGIN name=halo_26_rdma_free mode=mpi_only pattern=teardown
    block
        integer :: free_slot, free_idx, free_ierr
        do free_slot = 1, halo_26_rdma_nvsizes
            do free_idx = 1, halo_26_neighbor_count
                if (halo_26_rdma_recv((free_slot - 1)*halo_26_max_neighbors + free_idx) /= MPI_REQUEST_NULL) then
                    call MPI_REQUEST_FREE(halo_26_rdma_recv((free_slot - 1)*halo_26_max_neighbors + free_idx), free_ierr)
                end if
                if (halo_26_rdma_send((free_slot - 1)*halo_26_max_neighbors + free_idx) /= MPI_REQUEST_NULL) then
                    call MPI_REQUEST_FREE(halo_26_rdma_send((free_slot - 1)*halo_26_max_neighbors + free_idx), free_ierr)
                end if
            end do
        end do
    end block
    ! MFC_COMM_END
#endif
#ifdef MFC_HALO_INDEX_CACHE
    @:DEALLOCATE(halo_26_indices)
#endif
#endif
#if defined(MFC_NVSHMEM) && defined(__NVCOMPILER)
    nullify (nvshmem_inbox)
    call mfc_nvshmem_finalize()
#endif
#ifdef MFC_IGR_BATCHED_EXCHANGE
    @:DEALLOCATE(igr_buff_send, igr_buff_recv)
#endif
#ifndef __NVCOMPILER_GPU_UNIFIED_MEM
    @:DEALLOCATE(buff_send, buff_recv)
#else
    $:GPU_EXIT_DATA(delete='[buff_send, buff_recv]')
    deallocate (buff_send, buff_recv)
#endif
#endif

end subroutine s_finalize_mpi_common_module

end module m_mpi_common
