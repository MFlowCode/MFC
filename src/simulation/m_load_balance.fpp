!>
!!@file
!!@brief Contains module m_load_balance

#:include 'macros.fpp'

!> @brief Weighted static Cartesian decomposition (Tier-2 sub-project C).
module m_load_balance

    use m_derived_types
    use m_global_parameters
    use m_mpi_common

    implicit none

    private
    public :: f_weighted_splits, s_load_balance_rebalance

contains

    !> Cumulative offsets splitting marginal w into n_parts contiguous chunks of near-equal weight, each at least l_min cells.
    !! off(0)=0, off(n_parts)=size(w). Feasibility (size(w) >= n_parts*l_min) must be checked by the caller; this function is pure
    !! and cannot call s_mpi_abort.
    pure function f_weighted_splits(w, n_parts, l_min) result(off)

        real(wp), dimension(0:), intent(in) :: w
        integer, intent(in)                 :: n_parts, l_min
        integer, dimension(0:n_parts)       :: off
        real(wp)                            :: csum, total, target_w
        integer                             :: g, i, r

        g = size(w)
        off(0) = 0
        off(n_parts) = g
        if (n_parts == 1) return

        total = sum(w)
        ! ideal weighted boundaries: smallest i with cumulative >= r*total/n_parts
        r = 1
        csum = 0._wp
        do i = 0, g - 1
            csum = csum + w(i)
            do while (r < n_parts .and. csum >= real(r, wp)*total/real(n_parts, wp))
                off(r) = i + 1
                r = r + 1
            end do
        end do
        do while (r < n_parts)  ! zero-weight tail: park remaining boundaries at g
            off(r) = g; r = r + 1
        end do
        ! enforce the l_min floor, left to right, keeping strictly increasing
        do r = 1, n_parts - 1
            if (off(r) < r*l_min) off(r) = r*l_min
            if (off(r) > g - (n_parts - r)*l_min) off(r) = g - (n_parts - r)*l_min
            if (off(r) <= off(r - 1)) off(r) = off(r - 1) + l_min
        end do

    end function f_weighted_splits

    !> Returns true if the cumulative offsets off(0:parts) differ from the equal-split boundaries for a global extent of g cells
    !! distributed over parts ranks (matching the remainder distribution used by s_mpi_decompose_computational_domain).
    pure function f_offsets_differ_from_equal(off, g, parts) result(differs)

        integer, dimension(0:), intent(in) :: off
        integer, intent(in)                :: g, parts
        logical                            :: differs
        integer                            :: r

        differs = .false.
        do r = 0, parts
            if (off(r) /= r*(g/parts) + min(r, mod(g, parts))) then
                differs = .true.
                return
            end if
        end do

    end function f_offsets_differ_from_equal

    !> Read the first advection variable at the equal layout from the restart file and build global per-axis marginals. wx(0:m_glb),
    !! wy(0:n_glb), wz(0:p_glb) are allocated here; caller deallocates. Requires eqn_idx%adv%beg to be valid (call
    !! s_initialize_eqn_idx before invoking this).
    impure subroutine s_probe_field_marginals(wx, wy, wz)

        real(wp), allocatable, dimension(:), intent(out) :: wx, wy, wz

#ifdef MFC_MPI
        real(stp), allocatable, dimension(:,:,:) :: probe
        real(wp), allocatable, dimension(:)      :: lx, ly, lz
        integer(MPI_OFFSET_KIND)                 :: m_MOK, n_MOK, p_MOK, WP_MOK, MOK, var_MOK, disp
        integer, dimension(3)                    :: sizes_glb, sizes_loc
        integer                                  :: view, ifile, ierr, v, data_size
        integer                                  :: lb_start(3), j, k, l
        character(LEN=path_len + 2*name_len)     :: file_loc
        character(len=10)                        :: step_str
        logical                                  :: file_exist
#endif

        allocate (wx(0:m_glb), wy(0:n_glb), wz(0:p_glb))
        wx = 1._wp; wy = 1._wp; wz = 1._wp

#ifdef MFC_MPI
        ! Build the restart file path (mirrors s_read_parallel_data_files non-file_per_process path)
        if (cfl_dt) then
            write (step_str, '(I0)') n_start
        else
            write (step_str, '(I0)') t_step_start
        end if
        write (file_loc, '(A)') trim(step_str) // '.dat'
        file_loc = trim(case_dir) // '/restart_data' // trim(mpiiofs) // trim(file_loc)

        inquire (FILE=trim(file_loc), EXIST=file_exist)
        if (.not. file_exist) then
            if (proc_rank == 0) print *, '[load_balance] probe: restart file missing, using equal decomposition: ' // trim(file_loc)
            return
        end if

        ! Pick first advection variable (void fraction / bubble-concentration signal)
        v = eqn_idx%adv%beg

        m_MOK = int(m_glb + 1, MPI_OFFSET_KIND)
        n_MOK = int(n_glb + 1, MPI_OFFSET_KIND)
        p_MOK = int(p_glb + 1, MPI_OFFSET_KIND)
        WP_MOK = int(storage_size(0._stp)/8, MPI_OFFSET_KIND)
        MOK = int(1._wp, MPI_OFFSET_KIND)
        var_MOK = int(v, MPI_OFFSET_KIND)
        disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

        ! Full 3-element assignment (conforming for all num_dims); MPI calls slice (1:num_dims).
        sizes_glb = [m_glb + 1, n_glb + 1, p_glb + 1]
        sizes_loc = [m + 1, n + 1, p + 1]

        call MPI_TYPE_CREATE_SUBARRAY(num_dims, sizes_glb(1:num_dims), sizes_loc(1:num_dims), start_idx(1:num_dims), &
                                      & MPI_ORDER_FORTRAN, mpi_io_p, view, ierr)
        call MPI_TYPE_COMMIT(view, ierr)

        call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(file_loc), MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
        call MPI_FILE_SET_VIEW(ifile, disp, mpi_io_p, view, 'native', mpi_info_int, ierr)

        data_size = (m + 1)*(n + 1)*(p + 1)
        allocate (probe(0:m,0:n,0:p))
        call MPI_FILE_READ_ALL(ifile, probe, data_size*mpi_io_type, mpi_io_p, MPI_STATUS_IGNORE, ierr)
        call MPI_FILE_CLOSE(ifile, ierr)
        call MPI_TYPE_FREE(view, ierr)

        ! Bin into local per-axis marginals (lb_start guards inactive dimensions)
        lb_start = 0
        lb_start(1) = start_idx(1)
        if (num_dims >= 2) lb_start(2) = start_idx(2)
        if (num_dims >= 3) lb_start(3) = start_idx(3)

        allocate (lx(0:m_glb), ly(0:n_glb), lz(0:p_glb))
        lx = 0._wp; ly = 0._wp; lz = 0._wp

        do l = 0, p
            do k = 0, n
                do j = 0, m
                    lx(lb_start(1) + j) = lx(lb_start(1) + j) + real(probe(j, k, l), wp)
                    ly(lb_start(2) + k) = ly(lb_start(2) + k) + real(probe(j, k, l), wp)
                    lz(lb_start(3) + l) = lz(lb_start(3) + l) + real(probe(j, k, l), wp)
                end do
            end do
        end do

        call MPI_ALLREDUCE(lx, wx, m_glb + 1, mpi_p, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(ly, wy, n_glb + 1, mpi_p, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(lz, wz, p_glb + 1, mpi_p, MPI_SUM, MPI_COMM_WORLD, ierr)

        deallocate (probe, lx, ly, lz)
#endif

    end subroutine s_probe_field_marginals

    !> One-shot weighted re-decomposition at init (no-op if load_balance is off or the weighted splits match the equal splits on
    !! every axis). Reads one advection variable from the restart file at the equal layout to build per-axis weight marginals.
    impure subroutine s_load_balance_rebalance()

        real(wp), allocatable, dimension(:) :: wx, wy, wz
        integer, allocatable, dimension(:)  :: ox, oy, oz
        integer                             :: lmin, recon_order
        logical                             :: changed

        if (.not. load_balance) return

        ! Populate eqn_idx so that s_probe_field_marginals can pick eqn_idx%adv%beg.
        ! s_initialize_eqn_idx is pure index arithmetic (no allocations); it is safe to call
        ! here and again at its normal site in s_initialize_global_parameters_module.
        call s_initialize_eqn_idx(nmom, nb)

        if (recon_type == recon_type_weno) then
            recon_order = weno_order
        else
            recon_order = muscl_order
        end if
        if (igr) recon_order = igr_order

        lmin = num_stcls_min*recon_order
        call s_probe_field_marginals(wx, wy, wz)

        ! Only axes actually split across >1 ranks must satisfy the min-cells floor;
        ! a single-rank (incl. collapsed 1D/2D) axis owns all its cells and is always feasible.
        @:PROHIBIT(num_procs_x > 1 .and. (m_glb + 1) < num_procs_x*lmin, "load_balance: x-axis too small for min cells per rank")
        @:PROHIBIT(num_procs_y > 1 .and. (n_glb + 1) < num_procs_y*lmin, "load_balance: y-axis too small for min cells per rank")
        @:PROHIBIT(num_procs_z > 1 .and. (p_glb + 1) < num_procs_z*lmin, "load_balance: z-axis too small for min cells per rank")

        allocate (ox(0:num_procs_x), oy(0:num_procs_y), oz(0:num_procs_z))
        ox = f_weighted_splits(wx, num_procs_x, lmin)
        oy = f_weighted_splits(wy, num_procs_y, lmin)
        oz = f_weighted_splits(wz, num_procs_z, lmin)

        changed = f_offsets_differ_from_equal(ox, m_glb + 1, num_procs_x) .or. f_offsets_differ_from_equal(oy, n_glb + 1, &
                                              & num_procs_y) .or. f_offsets_differ_from_equal(oz, p_glb + 1, num_procs_z)

        if (proc_rank == 0) then
            print *, '[load_balance] x-offsets:', ox
            if (num_dims >= 2) print *, '[load_balance] y-offsets:', oy
            if (num_dims >= 3) print *, '[load_balance] z-offsets:', oz
        end if

        if (changed) call s_apply_weighted_offsets(ox, oy, oz)

        deallocate (wx, wy, wz, ox, oy, oz)

    end subroutine s_load_balance_rebalance

end module m_load_balance
