!>
!! @file m_derived_types.f90
!! @brief Contains module m_derived_types

!> @brief This file contains the definitions of all of the custom-defined
!!              types used in the pre-process code.
module m_helper

    ! Dependencies =============================================================
    
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    ! ==========================================================================

    implicit none

    real(kind(0d0)) :: cart_y, cart_z
    real(kind(0d0)) :: sph_phi !<
    !! Variables to be used to hold cell locations in Cartesian coordinates if
    !! 3D simulation is using cylindrical coordinates

    private; public :: s_convert_cylindrical_to_cartesian_coord, &
        s_convert_cylindrical_to_spherical_coord, &
        f_r

contains

    subroutine s_convert_cylindrical_to_cartesian_coord(cyl_y, cyl_z)

        real(kind(0d0)), intent(IN) :: cyl_y, cyl_z

        cart_y = cyl_y*sin(cyl_z)
        cart_z = cyl_y*cos(cyl_z)

    end subroutine s_convert_cylindrical_to_cartesian_coord ! --------------

    subroutine s_convert_cylindrical_to_spherical_coord(cyl_x, cyl_y)

        real(kind(0d0)), intent(IN) :: cyl_x, cyl_y

        sph_phi = atan(cyl_y/cyl_x)

    end subroutine s_convert_cylindrical_to_spherical_coord ! --------------

    !> Archimedes spiral function
        !! @param myth Angle
        !! @param offset Thickness
        !! @param a Starting position
    function f_r(myth, offset, a)
        real(kind(0d0)), intent(IN) :: myth, offset, a
        real(kind(0d0)) :: b
        real(kind(0d0)) :: f_r

        !r(th) = a + b*th

        b = 2.d0*a/(2.d0*pi)
        f_r = a + b*myth + offset
    end function f_r

end module m_helper