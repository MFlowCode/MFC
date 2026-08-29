!>
!!@file
!!@brief Contains module m_sfc_partition

#:include 'macros.fpp'

!> @brief Analysis-only weighted space-filling-curve partitioner.
module m_sfc_partition

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_mpi_common
    use m_load_weight, only: load_weight
    use m_box, only: f_morton

    implicit none

    private
    public :: s_initialize_sfc_partition_module, s_finalize_sfc_partition_module, s_compute_sfc_partition, s_report_sfc_partition

    integer               :: n_tiles_x, n_tiles_y, n_tiles_z, n_tiles
    real(wp), allocatable :: tile_weight(:)  !< global per-tile aggregated cost (linear index)
    integer, allocatable  :: tile_rank(:)    !< proposed owning rank per tile
    integer, allocatable  :: sfc_order(:)    !< tile linear indices in Morton order
    real(wp)              :: cur_w_max       !< current static per-rank max weight (existing decomposition)

contains

    impure subroutine s_initialize_sfc_partition_module

        if (.not. sfc_partition_wrt) return
        n_tiles_x = (m_glb + 1 + partition_tile_size - 1)/partition_tile_size
        n_tiles_y = (n_glb + 1 + partition_tile_size - 1)/partition_tile_size
        n_tiles_z = (p_glb + 1 + partition_tile_size - 1)/partition_tile_size
        n_tiles = n_tiles_x*n_tiles_y*n_tiles_z
        @:ALLOCATE(tile_weight(0:n_tiles - 1))
        @:ALLOCATE(tile_rank(0:n_tiles - 1))
        @:ALLOCATE(sfc_order(1:n_tiles))

    end subroutine s_initialize_sfc_partition_module

    impure subroutine s_finalize_sfc_partition_module

        if (.not. sfc_partition_wrt) return
        @:DEALLOCATE(tile_weight)
        @:DEALLOCATE(tile_rank)
        @:DEALLOCATE(sfc_order)

    end subroutine s_finalize_sfc_partition_module

    !> Aggregate per-cell load weights (computed by the caller) into global per-tile weights via MPI_ALLREDUCE.
    impure subroutine s_compute_sfc_partition()

        real(wp), allocatable :: tile_local(:)
        integer               :: j, k, l, gx, gy, gz, tx, ty, tz, t, ierr
        real(wp)              :: my_w
        integer               :: sfc_start(3)  !< bounds-safe copy of start_idx (0 for inactive dims)

        if (.not. sfc_partition_wrt) return

        ! 3-element start offset safe to access regardless of num_dims (start_idx is allocated
        ! (1:num_dims) only; higher dims are always 0).
        sfc_start = 0
        sfc_start(1) = start_idx(1)
        if (num_dims >= 2) sfc_start(2) = start_idx(2)
        if (num_dims >= 3) sfc_start(3) = start_idx(3)

        allocate (tile_local(0:n_tiles - 1)); tile_local = 0._wp
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    gx = sfc_start(1) + j
                    gy = sfc_start(2) + k
                    gz = sfc_start(3) + l
                    tx = gx/partition_tile_size; ty = gy/partition_tile_size; tz = gz/partition_tile_size
                    t = (tz*n_tiles_y + ty)*n_tiles_x + tx
                    tile_local(t) = tile_local(t) + real(load_weight%sf(j, k, l), wp)
                end do
            end do
        end do
#ifdef MFC_MPI
        call MPI_ALLREDUCE(tile_local, tile_weight, n_tiles, mpi_p, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
        tile_weight = tile_local
#endif
        my_w = sum(tile_local)
#ifdef MFC_MPI
        call MPI_ALLREDUCE(my_w, cur_w_max, 1, mpi_p, MPI_MAX, MPI_COMM_WORLD, ierr)
#else
        cur_w_max = my_w
#endif
        deallocate (tile_local)

        call s_build_sfc_order(sfc_order)

        block
            real(wp), allocatable :: wsfc(:)
            real(wp)              :: lo, hi, mid, acc
            integer               :: i, r, it
            allocate (wsfc(n_tiles))
            do i = 1, n_tiles; wsfc(i) = tile_weight(sfc_order(i)); end do
            ! binary search the smallest feasible max-load bound
            lo = maxval(wsfc); hi = sum(wsfc)
            do it = 1, 200
                if (hi - lo <= 1.e-12_wp*max(hi, 1._wp)) exit
                mid = 0.5_wp*(lo + hi)
                if (f_segments_needed(wsfc, mid) <= num_procs) then; hi = mid; else; lo = mid; end if
            end do
            ! assign ranks greedily with bound=hi, capping at num_procs-1 for the tail
            r = 0; acc = 0._wp
            do i = 1, n_tiles
                if (acc + wsfc(i) > hi .and. acc > 0._wp .and. r < num_procs - 1) then
                    r = r + 1; acc = wsfc(i)
                else
                    acc = acc + wsfc(i)
                end if
                tile_rank(sfc_order(i)) = r
            end do
            deallocate (wsfc)
        end block

    end subroutine s_compute_sfc_partition

    !> Greedy count of contiguous segments (each <= bound) over SFC-ordered weights.
    pure integer function f_segments_needed(wsfc, bound) result(nseg)

        real(wp), intent(in) :: wsfc(:)
        real(wp), intent(in) :: bound
        real(wp)             :: acc; integer :: i

        nseg = 1; acc = 0._wp
        do i = 1, size(wsfc)
            if (acc + wsfc(i) > bound .and. acc > 0._wp) then
                nseg = nseg + 1; acc = wsfc(i)
            else
                acc = acc + wsfc(i)
            end if
        end do

    end function f_segments_needed

    !> Fills order(1:n_tiles) with tile linear indices sorted by Morton code.
    impure subroutine s_build_sfc_order(order)

        integer, intent(out)         :: order(:)
        integer(kind=8), allocatable :: code(:)
        integer                      :: tx, ty, tz, t, i
        integer                      :: width, lo_r, mid_r, hi_r, ia, ib_m, iw
        integer, allocatable         :: work(:)

        allocate (code(0:n_tiles - 1))
        do tz = 0, n_tiles_z - 1
            do ty = 0, n_tiles_y - 1
                do tx = 0, n_tiles_x - 1
                    t = (tz*n_tiles_y + ty)*n_tiles_x + tx
                    code(t) = f_morton(tx, ty, tz)
                end do
            end do
        end do
        ! Index sort by Morton code: bottom-up mergesort, O(n_tiles log n_tiles). A selection loop is
        ! prohibitive at fine tile sizes (n_tiles reaches 1e6+ on large 3D grids). Codes are unique
        ! (one per tile), so the order is deterministic on every rank.
        do i = 1, n_tiles
            order(i) = i - 1
        end do
        allocate (work(1:n_tiles))
        width = 1
        do while (width < n_tiles)
            do lo_r = 1, n_tiles, 2*width
                mid_r = min(lo_r + width - 1, n_tiles)
                hi_r = min(lo_r + 2*width - 1, n_tiles)
                if (mid_r >= hi_r) cycle
                ia = lo_r; ib_m = mid_r + 1; iw = lo_r
                do while (ia <= mid_r .and. ib_m <= hi_r)
                    if (code(order(ia)) <= code(order(ib_m))) then
                        work(iw) = order(ia); ia = ia + 1
                    else
                        work(iw) = order(ib_m); ib_m = ib_m + 1
                    end if
                    iw = iw + 1
                end do
                do while (ia <= mid_r); work(iw) = order(ia); ia = ia + 1; iw = iw + 1; end do
                do while (ib_m <= hi_r); work(iw) = order(ib_m); ib_m = ib_m + 1; iw = iw + 1; end do
                order(lo_r:hi_r) = work(lo_r:hi_r)
            end do
            width = 2*width
        end do
        deallocate (work)
        deallocate (code)

    end subroutine s_build_sfc_order

    !> Print current static imbalance, predicted post-balance imbalance, and gain ratio.
    impure subroutine s_report_sfc_partition

        real(wp), allocatable :: rank_w(:)
        real(wp)              :: w_sum, w_max, w_mean, imb_new, imb_cur, gain
        integer               :: t

        if (.not. sfc_partition_wrt) return
        allocate (rank_w(0:num_procs - 1)); rank_w = 0._wp
        do t = 0, n_tiles - 1
            rank_w(tile_rank(t)) = rank_w(tile_rank(t)) + tile_weight(t)
        end do
        w_sum = sum(rank_w); w_max = maxval(rank_w); w_mean = w_sum/real(num_procs, wp)
        imb_new = w_max/max(w_mean, tiny(1._wp))
        imb_cur = cur_w_max/max(w_mean, tiny(1._wp))
        gain = imb_cur/max(imb_new, tiny(1._wp))
        if (proc_rank == 0) then
            print '(A,F8.3,A,F8.3,A,F8.3,A,I0,A,I0,A)', '[sfc_partition] imbalance current=', imb_cur, ' predicted=', imb_new, &
                & ' gain=', gain, '  (', n_tiles, ' tiles over ', num_procs, ' ranks)'
        end if
        deallocate (rank_w)

    end subroutine s_report_sfc_partition

end module m_sfc_partition
