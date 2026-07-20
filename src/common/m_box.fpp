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
    public :: t_box, f_equal_splits, f_weighted_splits, f_box_from_splits, f_morton

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
    !! off(n_parts)=size(w). Feasibility (size(w) >= n_parts*l_min) is the caller's responsibility (pure; no abort). A degenerate
    !! marginal (sum(w) <= 0) falls back to the equal split.
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
        if (total <= 0._wp) then
            off = f_equal_splits(g, n_parts)
            return
        end if
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
        ! Enforce the l_min floor: gap push first (off(0)=0 makes it imply off(r) >= r*l_min inductively), then the
        ! upper clamp - which cannot re-break the gap, since off(r-1) <= g - (n_parts-r+1)*l_min after its own pass.
        do r = 1, n_parts - 1
            if (off(r) < off(r - 1) + l_min) off(r) = off(r - 1) + l_min
            if (off(r) > g - (n_parts - r)*l_min) off(r) = g - (n_parts - r)*l_min
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

    !> 3D Morton (Z-order) key interleaving the bits of (ix, iy, iz); collapsed dims contribute 0. 21 bits/dim (fits a 64-bit key
    !! for grids up to 2^21 cells/dim). Shared by the AMR block partition (m_amr) and the SFC-partition diagnostic
    !! (m_sfc_partition); negative coordinates clamp to 0.
    pure integer(kind=8) function f_morton(ix, iy, iz) result(key)
        integer, intent(in) :: ix, iy, iz
        integer             :: b
        integer(kind=8)     :: xx, yy, zz

        xx = int(max(ix, 0), 8); yy = int(max(iy, 0), 8); zz = int(max(iz, 0), 8)
        key = 0_8
        do b = 0, 20
            key = ior(key, ishft(iand(ishft(xx, -b), 1_8), 3*b))
            key = ior(key, ishft(iand(ishft(yy, -b), 1_8), 3*b + 1))
            key = ior(key, ishft(iand(ishft(zz, -b), 1_8), 3*b + 2))
        end do

    end function f_morton

end module m_box
