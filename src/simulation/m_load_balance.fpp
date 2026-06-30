!>
!!@file
!!@brief Contains module m_load_balance

#:include 'macros.fpp'

!> @brief Weighted static Cartesian decomposition (Tier-2 sub-project C).
module m_load_balance

    use m_derived_types
    use m_global_parameters
    use m_mpi_common
    use m_load_weight, only: load_weight, s_compute_load_weight

    implicit none

    private
    public :: f_weighted_splits, s_compute_load_marginals, s_load_balance_rebalance

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

    !> Compute global per-axis marginals of the per-cell load weight. Allocates wx(0:m_glb), wy(0:n_glb), wz(0:p_glb); caller
    !! deallocates.
    impure subroutine s_compute_load_marginals(q_cons_vf, wx, wy, wz)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        real(wp), allocatable, dimension(:), intent(out)    :: wx, wy, wz
        real(wp), allocatable, dimension(:)                 :: lx, ly, lz
        integer                                             :: lb_start(3)  !< bounds-safe copy of start_idx (0 for inactive dims)
        integer                                             :: j, k, l, ierr

        call s_compute_load_weight(q_cons_vf)
        $:GPU_UPDATE(host='[load_weight%sf]')

        ! Build a 3-element start offset safe to access regardless of num_dims.
        ! start_idx is allocated (1:num_dims) only; higher dims are always 0.
        lb_start = 0
        lb_start(1) = start_idx(1)
        if (num_dims >= 2) lb_start(2) = start_idx(2)
        if (num_dims >= 3) lb_start(3) = start_idx(3)

        allocate (wx(0:m_glb), lx(0:m_glb)); lx = 0._wp
        allocate (wy(0:n_glb), ly(0:n_glb)); ly = 0._wp
        allocate (wz(0:p_glb), lz(0:p_glb)); lz = 0._wp

        do l = 0, p
            do k = 0, n
                do j = 0, m
                    lx(lb_start(1) + j) = lx(lb_start(1) + j) + real(load_weight%sf(j, k, l), wp)
                    ly(lb_start(2) + k) = ly(lb_start(2) + k) + real(load_weight%sf(j, k, l), wp)
                    lz(lb_start(3) + l) = lz(lb_start(3) + l) + real(load_weight%sf(j, k, l), wp)
                end do
            end do
        end do
#ifdef MFC_MPI
        call MPI_ALLREDUCE(lx, wx, m_glb + 1, mpi_p, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(ly, wy, n_glb + 1, mpi_p, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(lz, wz, p_glb + 1, mpi_p, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
        wx = lx; wy = ly; wz = lz
#endif
        deallocate (lx, ly, lz)

    end subroutine s_compute_load_marginals

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

    !> One-shot weighted re-decomposition at init (no-op if load_balance is off or the weighted splits match the equal splits on
    !! every axis).
    impure subroutine s_load_balance_rebalance(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        real(wp), allocatable, dimension(:)                 :: wx, wy, wz
        integer, allocatable, dimension(:)                  :: ox, oy, oz
        integer                                             :: lmin, recon_order
        logical                                             :: changed

        if (.not. load_balance) return

        if (recon_type == recon_type_weno) then
            recon_order = weno_order
        else
            recon_order = muscl_order
        end if
        if (igr) recon_order = igr_order

        lmin = num_stcls_min*recon_order
        call s_compute_load_marginals(q_cons_vf, wx, wy, wz)

        @:PROHIBIT((m_glb + 1) < num_procs_x*lmin, "load_balance: x-axis too small for min cells per rank")
        @:PROHIBIT((n_glb + 1) < num_procs_y*lmin, "load_balance: y-axis too small for min cells per rank")
        @:PROHIBIT((p_glb + 1) < num_procs_z*lmin, "load_balance: z-axis too small for min cells per rank")

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
