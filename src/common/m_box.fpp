!>
!!@file
!!@brief Contains module m_box

#:include 'macros.fpp'

!> @brief Owned domain-decomposition Box abstraction and partition arithmetic (v1: one box per rank).
module m_box

    use m_derived_types, only: t_box
    use m_global_parameters, only: wp

    implicit none

    private
    public :: t_box, f_equal_splits, f_weighted_splits, f_box_from_splits

contains

    !> Cumulative equal-cell offsets for g cells over n_parts ranks: off(r) = r*(g/n_parts) + min(r, mod(g,n_parts)). Reproduces
    !! MFC's block distribution (remainder to the first ranks) exactly. Pure integer path.
    pure function f_equal_splits(g, n_parts) result(off)

        integer, intent(in)           :: g, n_parts
        integer, dimension(0:n_parts) :: off
        integer                       :: q, rem, r

        q = g/n_parts
        rem = mod(g, n_parts)
        do r = 0, n_parts
            off(r) = r*q + min(r, rem)
        end do

    end function f_equal_splits

    !> Cumulative offsets splitting marginal w into n_parts contiguous chunks of near-equal weight, each >= l_min cells. off(0)=0,
    !! off(n_parts)=size(w). Feasibility (size(w) >= n_parts*l_min) is the caller's responsibility (pure; no abort).
    pure function f_weighted_splits(w, n_parts, l_min) result(off)

        real(wp), dimension(0:), intent(in) :: w
        integer, intent(in)                 :: n_parts, l_min
        integer, dimension(0:n_parts)       :: off
        real(wp)                            :: csum, total
        integer                             :: g, i, r

        g = size(w)
        off(0) = 0
        off(n_parts) = g
        if (n_parts == 1) return
        total = sum(w)
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
            off(r) = g; r = r + 1
        end do
        do r = 1, n_parts - 1
            if (off(r) < r*l_min) off(r) = r*l_min
            if (off(r) > g - (n_parts - r)*l_min) off(r) = g - (n_parts - r)*l_min
            if (off(r) <= off(r - 1)) off(r) = off(r - 1) + l_min
        end do

    end function f_weighted_splits

    !> Assemble this rank's box from per-axis cumulative offsets and the rank's Cartesian coords (0-based). lo(d) =
    !! off_d(coords(d)); hi(d) = off_d(coords(d)+1) - 1. Works for collapsed axes (off_d = [0,1] -> lo=hi=0).
    pure function f_box_from_splits(off_x, off_y, off_z, coords) result(box)

        integer, dimension(0:), intent(in) :: off_x, off_y, off_z
        integer, intent(in)                :: coords(3)
        type(t_box)                        :: box

        box%lo(1) = off_x(coords(1)); box%hi(1) = off_x(coords(1) + 1) - 1
        box%lo(2) = off_y(coords(2)); box%hi(2) = off_y(coords(2) + 1) - 1
        box%lo(3) = off_z(coords(3)); box%hi(3) = off_z(coords(3) + 1) - 1

    end function f_box_from_splits

end module m_box
