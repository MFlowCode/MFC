!>
!!@file
!!@brief Contains module m_load_balance

#:include 'macros.fpp'

!> @brief Weighted static Cartesian decomposition.
module m_load_balance

    use m_derived_types
    use m_global_parameters
    use m_mpi_common
    use m_box, only: f_equal_splits, f_weighted_splits

    implicit none

    private
    public :: s_load_balance_rebalance

contains

    !> True if cumulative offsets off(0:parts) differ from the equal-split boundaries for g cells over parts ranks (same remainder
    !! distribution as s_mpi_decompose_computational_domain).
    pure function f_offsets_differ_from_equal(off, g, parts) result(differs)

        integer, dimension(0:), intent(in) :: off
        integer, intent(in)                :: g, parts
        logical                            :: differs

        differs = any(off /= f_equal_splits(g, parts))

    end function f_offsets_differ_from_equal

    !> True iff, for one axis's cumulative offsets, every part's AMR-block intersection satisfies the mirror-decomposition scratch
    !! cap 2*(isect cells) - 1 <= part cells - 1 (the per-dim bound s_initialize_amr_module aborts on). Per-axis-part checks equal
    !! per-rank-box checks: a rank owns the block iff its intersection is nonempty in every active dim, and each per-dim
    !! intersection depends only on that axis's part index.
    pure function f_amr_isect_fits(off, n_parts, pbeg, pend) result(ok)

        integer, dimension(0:), intent(in) :: off
        integer, intent(in)                :: n_parts, pbeg, pend
        logical                            :: ok
        integer                            :: r, ilo, ihi

        ok = .true.
        do r = 0, n_parts - 1
            ilo = max(pbeg, off(r)); ihi = min(pend, off(r + 1) - 1)
            if (ihi >= ilo .and. 2*(ihi - ilo + 1) - 1 > off(r + 1) - off(r) - 1) ok = .false.
        end do

    end function f_amr_isect_fits

    !> Read the first advection variable at the equal layout from the restart file, build global per-axis marginals. Allocates
    !! wx(0:m_glb), wy(0:n_glb), wz(0:p_glb); caller deallocates. Requires a valid eqn_idx%adv%beg (call s_initialize_eqn_idx
    !! first).
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
        ! Restart file path (mirrors s_read_parallel_data_files non-file_per_process path)
        if (cfl_dt) then
            write (step_str, '(I0)') n_start
        else
            write (step_str, '(I0)') t_step_start
        end if
        write (file_loc, '(A)') trim(step_str) // '.dat'
        file_loc = trim(case_dir) // '/restart_data' // trim(mpiiofs) // trim(file_loc)

        inquire (FILE=trim(file_loc), EXIST=file_exist)
        if (.not. file_exist) then
#ifdef MFC_DEBUG
            if (proc_rank == 0) print *, '[load_balance] probe: restart file missing, using equal decomposition: ' // trim(file_loc)
#endif
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

        real(wp), allocatable, dimension(:) :: wx, wy, wz, vx, vy, vz
        integer, allocatable, dimension(:)  :: ox, oy, oz
        real(wp)                            :: amr_w(3), ff, scale
        integer                             :: lmin, recon_order, pc(3), try
        logical                             :: changed

        if (.not. load_balance) return

        ! Populate eqn_idx so s_probe_field_marginals can pick eqn_idx%adv%beg. s_initialize_eqn_idx
        ! is pure index arithmetic (no allocations); safe to call here and again at its normal site
        ! in s_initialize_global_parameters_module.
        call s_initialize_eqn_idx(nmom, nb)

        if (recon_type == recon_type_weno) then
            recon_order = weno_order
        else
            recon_order = muscl_order
        end if
        if (igr) recon_order = igr_order

        lmin = num_stcls_min*recon_order
        if (bubbles_euler) then
            ! EE-bubble source cost is flat across void magnitude (calibration note in m_load_weight),
            ! so the first-advection-alpha marginal is not a load signal here: use uniform marginals
            ! and let only the AMR fine-work injection move split planes.
            allocate (wx(0:m_glb), wy(0:n_glb), wz(0:p_glb))
            wx = 1._wp; wy = 1._wp; wz = 1._wp
        else
            call s_probe_field_marginals(wx, wy, wz)
        end if

        ! Only axes split across >1 ranks must satisfy the min-cells floor; a single-rank axis
        ! (incl. collapsed 1D/2D) owns all its cells and is always feasible.
        @:PROHIBIT(num_procs_x > 1 .and. (m_glb + 1) < num_procs_x*lmin, "load_balance: x-axis too small for min cells per rank")
        @:PROHIBIT(num_procs_y > 1 .and. (n_glb + 1) < num_procs_y*lmin, "load_balance: y-axis too small for min cells per rank")
        @:PROHIBIT(num_procs_z > 1 .and. (p_glb + 1) < num_procs_z*lmin, "load_balance: z-axis too small for min cells per rank")

        allocate (ox(0:num_procs_x), oy(0:num_procs_y), oz(0:num_procs_z))
        allocate (vx(0:m_glb), vy(0:n_glb), vz(0:p_glb))

        ! AMR fine-work injection: each block-covered coarse cell costs an extra 2**num_dims (interleaved;
        ! x2 subcycled) fine-cell RHS evals per coarse step, so the block's axis projections gain
        ! fine_factor * (block transverse cells) in each marginal's own units (per-cell scale =
        ! sum(w_axis)/total cells per axis: the probe's marginal sums are all equal, the missing-file
        ! fallback's are per-index and differ across axes).
        amr_w = 0._wp
        pc = 1
        if (amr) then
            pc(1) = amr_block_end(1) - amr_block_beg(1) + 1
            if (n_glb > 0) pc(2) = amr_block_end(2) - amr_block_beg(2) + 1
            if (p_glb > 0) pc(3) = amr_block_end(3) - amr_block_beg(3) + 1
            ff = real(2**(num_dims + merge(1, 0, amr_subcycle)), wp)/(real(m_glb + 1, wp)*real(n_glb + 1, wp)*real(p_glb + 1, wp))
            amr_w(1) = ff*sum(wx)*real(pc(2), wp)*real(pc(3), wp)
            amr_w(2) = ff*sum(wy)*real(pc(1), wp)*real(pc(3), wp)
            amr_w(3) = ff*sum(wz)*real(pc(1), wp)*real(pc(2), wp)
        end if

        ! Deterministic feasibility clamp: if the weighted boxes would violate the AMR scratch cap on any
        ! rank, halve the amr contribution and re-split (<= 3 retries; the final fallback drops it, governed
        ! by s_initialize_amr_module's own init check). Pure arithmetic on rank-identical arrays, so every
        ! rank takes the same branches.
        scale = 1._wp
        do try = 1, 4
            if (try == 4) scale = 0._wp
            vx = wx; vy = wy; vz = wz
            if (amr .and. scale > 0._wp) then
                vx(amr_block_beg(1):amr_block_end(1)) = vx(amr_block_beg(1):amr_block_end(1)) + scale*amr_w(1)
                if (n_glb > 0) vy(amr_block_beg(2):amr_block_end(2)) = vy(amr_block_beg(2):amr_block_end(2)) + scale*amr_w(2)
                if (p_glb > 0) vz(amr_block_beg(3):amr_block_end(3)) = vz(amr_block_beg(3):amr_block_end(3)) + scale*amr_w(3)
            end if
            ox = f_weighted_splits(vx, num_procs_x, lmin)
            oy = f_weighted_splits(vy, num_procs_y, lmin)
            oz = f_weighted_splits(vz, num_procs_z, lmin)
            if (.not. amr) exit
            if (f_amr_isect_fits(ox, num_procs_x, amr_block_beg(1), amr_block_end(1)) .and. (n_glb == 0 .or. f_amr_isect_fits(oy, &
                & num_procs_y, amr_block_beg(2), amr_block_end(2))) .and. (p_glb == 0 .or. f_amr_isect_fits(oz, num_procs_z, &
                & amr_block_beg(3), amr_block_end(3)))) exit
            scale = 0.5_wp*scale
        end do
#ifdef MFC_DEBUG
        if (amr .and. scale < 1._wp .and. proc_rank == 0) then
            print *, '[load_balance] amr weight softened to fit the fine block per rank; scale =', scale
        end if
#endif

        changed = f_offsets_differ_from_equal(ox, m_glb + 1, num_procs_x) .or. f_offsets_differ_from_equal(oy, n_glb + 1, &
                                              & num_procs_y) .or. f_offsets_differ_from_equal(oz, p_glb + 1, num_procs_z)

        if (changed) call s_apply_weighted_offsets(ox, oy, oz)

        deallocate (wx, wy, wz, vx, vy, vz, ox, oy, oz)

    end subroutine s_load_balance_rebalance

end module m_load_balance
