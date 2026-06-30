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
    public :: f_weighted_splits

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

end module m_load_balance
