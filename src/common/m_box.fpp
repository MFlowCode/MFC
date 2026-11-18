!>
!!@file
!!@brief Contains module m_box

#:include 'macros.fpp'

!> @brief Box type, per-rank partition arithmetic (one box per rank) and the Morton key.
module m_box

    use m_derived_types, only: t_box
    use m_global_parameters, only: wp

    implicit none

    private
    public :: t_box, f_equal_splits, f_box_from_splits, f_morton

contains

    !> Cumulative equal-cell offsets for g cells over n_parts ranks, remainder to the first ranks (MFC's block distribution).
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

    !> This rank's box from per-axis cumulative offsets and its 0-based Cartesian coords (collapsed axes: off_d = [0,1]).
    pure function f_box_from_splits(off_x, off_y, off_z, coords) result(box)

        integer, dimension(0:), intent(in) :: off_x, off_y, off_z
        integer, intent(in)                :: coords(3)
        type(t_box)                        :: box

        box%lo(1) = off_x(coords(1)); box%hi(1) = off_x(coords(1) + 1) - 1
        box%lo(2) = off_y(coords(2)); box%hi(2) = off_y(coords(2) + 1) - 1
        box%lo(3) = off_z(coords(3)); box%hi(3) = off_z(coords(3) + 1) - 1

    end function f_box_from_splits

    !> 3D Morton (Z-order) key of (ix, iy, iz), 21 bits/dim; negative coords clamp to 0.
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
