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
    public :: f_weighted_splits, s_compute_load_marginals

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

end module m_load_balance
