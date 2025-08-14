#:include 'macros.fpp'

!> @brief This module contains subroutines that read, and check consistency
!!              of, the user provided inputs and data.

#:include 'macros.fpp'

module m_check_patches

    ! Dependencies
    use m_derived_types          !< Definitions of the derived types

    use m_global_parameters      !< Global parameters for the code

    use m_mpi_proxy              !< Message passing interface (MPI) module proxy

    use m_data_output            !< Procedures to write the grid data and the
                                 !! conservative variables to files

#ifdef MFC_MPI
    use mpi                      !< Message passing interface (MPI) module
#endif

    use m_compile_specific

    use m_helper_basic           !< Functions to compare floating point numbers

    use m_helper

    implicit none

    private; public :: s_check_patches

    character(len=10) :: iStr

contains

    impure subroutine s_check_patches

        integer :: i
        character(len=10) :: num_patches_str

        call s_int_to_str(num_patches, num_patches_str)

        do i = 1, num_patches_max
            if (i <= num_patches) then
                ! call s_check_patch_geometry(i)
                call s_int_to_str(i, iStr)
                @:PROHIBIT(patch_icpp(i)%geometry == 6, "Invalid patch geometry number. "// &
                    "patch_icpp("//trim(iStr)//")%geometry is deprecated.")
                @:PROHIBIT(patch_icpp(i)%geometry == 7, "Invalid patch geometry number. "// &
                    "patch_icpp("//trim(iStr)//")%geometry is deprecated.")
                @:PROHIBIT(patch_icpp(i)%geometry == 13, "Invalid patch geometry number. "// &
                    "patch_icpp("//trim(iStr)//")%geometry is deprecated.")
                @:PROHIBIT(patch_icpp(i)%geometry == 15, "Invalid patch geometry number. "// &
                    "patch_icpp("//trim(iStr)//")%geometry is deprecated.")
                @:PROHIBIT(patch_icpp(i)%geometry == dflt_int, "Invalid patch geometry number. "// &
                    "patch_icpp("//trim(iStr)//")%geometry must be set.")

                ! Constraints on the geometric initial condition patch parameters
                if (patch_icpp(i)%geometry == 1) then
                    call s_check_line_segment_patch_geometry(i)
                elseif (patch_icpp(i)%geometry == 2) then
                    call s_check_circle_patch_geometry(i)
                elseif (patch_icpp(i)%geometry == 3) then
                    call s_check_rectangle_patch_geometry(i)
                elseif (patch_icpp(i)%geometry == 4) then
                    call s_check_line_sweep_patch_geometry(i)
                elseif (patch_icpp(i)%geometry == 5) then
                    call s_check_ellipse_patch_geometry(i)
                elseif (patch_icpp(i)%geometry == 8) then
                    call s_check_sphere_patch_geometry(i)
                elseif (patch_icpp(i)%geometry == 9) then
                    call s_check_cuboid_patch_geometry(i)
                elseif (patch_icpp(i)%geometry == 10) then
                    call s_check_cylinder_patch_geometry(i)
                elseif (patch_icpp(i)%geometry == 11) then
                    call s_check_plane_sweep_patch_geometry(i)
                elseif (patch_icpp(i)%geometry == 12) then
                    call s_check_ellipsoid_patch_geometry(i)
                elseif (patch_icpp(i)%geometry == 14) then
                    call s_check_spherical_harmonic_patch_geometry(i)
                elseif (patch_icpp(i)%geometry == 20) then
                    call s_check_2D_TaylorGreen_vortex_patch_geometry(i)
                elseif (patch_icpp(i)%geometry == 21) then
                    call s_check_model_geometry(i)
                else
                    call s_prohibit_abort("Invalid patch geometry number", "patch_icpp("//trim(iStr)//")%geometry "// &
                                          "must be between 1 and 21")
                end if
            else
                @:PROHIBIT(patch_icpp(i)%geometry /= dflt_int, "Inactive patch defined. "// &
                    "patch_icpp("//trim(iStr)//")%geometry not be set for inactive patches. "// &
                    "Patch "//trim(iStr)//" is inactive as the number of patches is "//trim(num_patches_str))
                call s_check_inactive_patch_geometry(i)
            end if
        end do

        ! Constraints on overwrite rights initial condition patch parameters
        do i = 1, num_patches
            if (i <= num_patches) then
                call s_check_active_patch_alteration_rights(i)
            else
                call s_check_inactive_patch_alteration_rights(i)
            end if
        end do

        ! Constraints on smoothing initial condition patch parameters
        do i = 1, num_patches
            if (i > 1 .and. (patch_icpp(i)%geometry == 2 .or. &
                             patch_icpp(i)%geometry == 3 .or. &
                             patch_icpp(i)%geometry == 4 .or. &
                             patch_icpp(i)%geometry == 5 .or. &
                             patch_icpp(i)%geometry == 8 .or. &
                             patch_icpp(i)%geometry == 0 .or. &
                             patch_icpp(i)%geometry == 0 .or. &
                             patch_icpp(i)%geometry == 9 .or. &
                             patch_icpp(i)%geometry == 10 .or. &
                             patch_icpp(i)%geometry == 11 .or. &
                             patch_icpp(i)%geometry == 12 .or. &
                             patch_icpp(i)%geometry == 14)) then
                call s_check_supported_patch_smoothing(i)
            else
                call s_check_unsupported_patch_smoothing(i)
            end if
        end do

        ! Constraints on flow variables initial condition patch parameters
        do i = 1, num_patches
            if (i <= num_patches) then
                call s_check_active_patch_primitive_variables(i)
            else
                call s_check_inactive_patch_primitive_variables(i)
            end if
        end do

    end subroutine s_check_patches

    !> This subroutine checks the line segment patch input
        !!  @param patch_id Patch identifier
    impure subroutine s_check_line_segment_patch_geometry(patch_id)

        integer, intent(in) :: patch_id
        call s_int_to_str(patch_id, iStr)

        @:PROHIBIT(n > 0, "Line segment patch "//trim(iStr)//": n must be zero")
        @:PROHIBIT(patch_icpp(patch_id)%length_x <= 0._wp, "Line segment patch "//trim(iStr)//": length_x must be greater than zero")
        @:PROHIBIT(f_is_default(patch_icpp(patch_id)%x_centroid), "Line segment patch "//trim(iStr)//": x_centroid must be set")
        @:PROHIBIT(cyl_coord, "Line segment patch "//trim(iStr)//": cyl_coord is not supported")

    end subroutine s_check_line_segment_patch_geometry

    !>  This subroutine checks the circle patch input
        !!  @param patch_id Patch identifier
    impure subroutine s_check_circle_patch_geometry(patch_id)

        integer, intent(in) :: patch_id
        call s_int_to_str(patch_id, iStr)

        @:PROHIBIT(n == 0, "Circle patch "//trim(iStr)//": n must be zero")
        @:PROHIBIT(p > 0, "Circle patch "//trim(iStr)//": p must be greater than zero")
        @:PROHIBIT(patch_icpp(patch_id)%radius <= 0._wp, "Circle patch "//trim(iStr)//": radius must be greater than zero")
        @:PROHIBIT(f_is_default(patch_icpp(patch_id)%x_centroid), "Circle patch "//trim(iStr)//": x_centroid must be set")
        @:PROHIBIT(f_is_default(patch_icpp(patch_id)%y_centroid), "Circle patch "//trim(iStr)//": y_centroid must be set")

    end subroutine s_check_circle_patch_geometry

    !>  This subroutine checks the rectangle patch input
        !!  @param patch_id Patch identifier
    impure subroutine s_check_rectangle_patch_geometry(patch_id)

        integer, intent(in) :: patch_id
        call s_int_to_str(patch_id, iStr)

        @:PROHIBIT(n == 0, "Rectangle patch "//trim(iStr)//": n must be greater than zero")
        @:PROHIBIT(p > 0, "Rectangle patch "//trim(iStr)//": p must be zero")
        @:PROHIBIT(f_is_default(patch_icpp(patch_id)%x_centroid), "Rectangle patch "//trim(iStr)//": x_centroid must be set")
        @:PROHIBIT(f_is_default(patch_icpp(patch_id)%y_centroid), "Rectangle patch "//trim(iStr)//": y_centroid must be set")
        @:PROHIBIT(patch_icpp(patch_id)%length_x <= 0._wp, "Rectangle patch "//trim(iStr)//": length_x must be greater than zero")
        @:PROHIBIT(patch_icpp(patch_id)%length_y <= 0._wp, "Rectangle patch "//trim(iStr)//": length_y must be greater than zero")

    end subroutine s_check_rectangle_patch_geometry

    !> This subroutine checks the line sweep patch input
        !!  @param patch_id Patch identifier
    impure subroutine s_check_line_sweep_patch_geometry(patch_id)

        integer, intent(in) :: patch_id
        call s_int_to_str(patch_id, iStr)

        @:PROHIBIT(n == 0, "Line sweep patch "//trim(iStr)//": n must be greater than zero")
        @:PROHIBIT(p > 0, "Line sweep patch "//trim(iStr)//": p must be zero")
        @:PROHIBIT(f_is_default(patch_icpp(patch_id)%x_centroid), "Line sweep patch "//trim(iStr)//": x_centroid must be set")
        @:PROHIBIT(f_is_default(patch_icpp(patch_id)%y_centroid), "Line sweep patch "//trim(iStr)//": y_centroid must be set")
        @:PROHIBIT(f_is_default(patch_icpp(patch_id)%normal(1)), "Line sweep patch "//trim(iStr)//": normal(1) must be set")
        @:PROHIBIT(f_is_default(patch_icpp(patch_id)%normal(2)), "Line sweep patch "//trim(iStr)//": normal(2) must be set")
        @:PROHIBIT(.not. f_is_default(patch_icpp(patch_id)%normal(3)), "Line sweep patch "//trim(iStr)//": normal(3) must not be set")

    end subroutine s_check_line_sweep_patch_geometry

    !>  This subroutine checks the ellipse patch input
        !!  @param patch_id Patch identifier
    impure subroutine s_check_ellipse_patch_geometry(patch_id)

        integer, intent(in) :: patch_id
        call s_int_to_str(patch_id, iStr)

        @:PROHIBIT(n == 0, "Ellipse patch "//trim(iStr)//": n must be greater than zero")
        @:PROHIBIT(p > 0, "Ellipse patch "//trim(iStr)//": p must be zero")
        @:PROHIBIT(f_is_default(patch_icpp(patch_id)%x_centroid), "Ellipse patch "//trim(iStr)//": x_centroid must be set")
        @:PROHIBIT(f_is_default(patch_icpp(patch_id)%y_centroid), "Ellipse patch "//trim(iStr)//": y_centroid must be set")
        @:PROHIBIT(patch_icpp(patch_id)%radii(1) <= 0._wp, "Ellipse patch "//trim(iStr)//": radii(1) must be greater than zero")
        @:PROHIBIT(patch_icpp(patch_id)%radii(2) <= 0._wp, "Ellipse patch "//trim(iStr)//": radii(2) must be greater than zero")
        @:PROHIBIT(.not. f_is_default(patch_icpp(patch_id)%radii(3)), "Ellipse patch "//trim(iStr)//": radii(3) must not be set")

    end subroutine s_check_ellipse_patch_geometry

    !>  This subroutine checks the model patch input
        !!  @param patch_id Patch identifier
    impure subroutine s_check_2D_TaylorGreen_vortex_patch_geometry(patch_id)

        integer, intent(in) :: patch_id
        call s_int_to_str(patch_id, iStr)

        @:PROHIBIT(n == 0, "Taylor Green vortex patch "//trim(iStr)//": n must be greater than zero")
        @:PROHIBIT(p > 0, "Taylor Green vortex patch "//trim(iStr)//": p must be zero")
        @:PROHIBIT(f_is_default(patch_icpp(patch_id)%x_centroid), "Taylor Green vortex patch "//trim(iStr)//": x_centroid must be set")
        @:PROHIBIT(f_is_default(patch_icpp(patch_id)%y_centroid), "Taylor Green vortex patch "//trim(iStr)//": y_centroid must be set")
        @:PROHIBIT(patch_icpp(patch_id)%length_x <= 0._wp, "Taylor Green vortex patch "//trim(iStr)//": length_x must be greater than zero")
        @:PROHIBIT(patch_icpp(patch_id)%length_y <= 0._wp, "Taylor Green vortex patch "//trim(iStr)//": length_y must be greater than zero")
        @:PROHIBIT(patch_icpp(patch_id)%vel(2) <= 0._wp, "Taylor Green vortex patch "//trim(iStr)//": vel(2) must be greater than zero")

    end subroutine s_check_2D_TaylorGreen_vortex_patch_geometry

    !> This subroutine checks the model patch input
        !!  @param patch_id Patch identifier
    impure subroutine s_check_sphere_patch_geometry(patch_id)

        integer, intent(in) :: patch_id
        call s_int_to_str(patch_id, iStr)

        @:PROHIBIT(p == 0, "Sphere patch "//trim(iStr)//": p must be greater than zero")
        @:PROHIBIT(patch_icpp(patch_id)%radius <= 0._wp, "Sphere patch "//trim(iStr)//": radius must be greater than zero")
        @:PROHIBIT(f_is_default(patch_icpp(patch_id)%x_centroid), "Sphere patch "//trim(iStr)//": x_centroid must be set")
        @:PROHIBIT(f_is_default(patch_icpp(patch_id)%y_centroid), "Sphere patch "//trim(iStr)//": y_centroid must be set")
        @:PROHIBIT(f_is_default(patch_icpp(patch_id)%z_centroid), "Sphere patch "//trim(iStr)//": z_centroid must be set")

    end subroutine s_check_sphere_patch_geometry

    !>  This subroutine checks the model patch input
        !!  @param patch_id Patch identifier
    impure subroutine s_check_spherical_harmonic_patch_geometry(patch_id)
        integer, intent(in) :: patch_id

        call s_int_to_str(patch_id, iStr)

        @:PROHIBIT(p == 0, "Spherical harmonic patch "//trim(iStr)//": p must be greater than zero")
        @:PROHIBIT(patch_icpp(patch_id)%radius <= 0._wp, "Spherical harmonic patch "//trim(iStr)//": radius must be greater than zero")
        @:PROHIBIT(f_is_default(patch_icpp(patch_id)%x_centroid), "Spherical harmonic patch "//trim(iStr)//": x_centroid must be set")
        @:PROHIBIT(f_is_default(patch_icpp(patch_id)%y_centroid), "Spherical harmonic patch "//trim(iStr)//": y_centroid must be set")
        @:PROHIBIT(f_is_default(patch_icpp(patch_id)%z_centroid), "Spherical harmonic patch "//trim(iStr)//": z_centroid must be set")
        @:PROHIBIT(.not. f_approx_in_array(patch_icpp(patch_id)%epsilon, (/1._wp, 2._wp, 3._wp, 4._wp, 5._wp/)), &
            "Spherical harmonic patch "//trim(iStr)//": epsilon must be one of 1, 2, 3, 4, 5")
        @:PROHIBIT(patch_icpp(patch_id)%beta < 0._wp, &
            "Spherical harmonic patch "//trim(iStr)//": beta must be greater than or equal to zero")
        @:PROHIBIT(patch_icpp(patch_id)%beta > patch_icpp(patch_id)%epsilon, &
            "Spherical harmonic patch "//trim(iStr)//": beta must be less than or equal to epsilon")

    end subroutine s_check_spherical_harmonic_patch_geometry

    !>  This subroutine checks the model patch input
        !!  @param patch_id Patch identifier
    impure subroutine s_check_cuboid_patch_geometry(patch_id)

        ! Patch identifier
        integer, intent(in) :: patch_id
        call s_int_to_str(patch_id, iStr)

        @:PROHIBIT(p == 0, "Cuboid patch "//trim(iStr)//": p must be greater than zero")
        @:PROHIBIT(f_is_default(patch_icpp(patch_id)%x_centroid), "Cuboid patch "//trim(iStr)//": x_centroid must be set")
        @:PROHIBIT(f_is_default(patch_icpp(patch_id)%y_centroid), "Cuboid patch "//trim(iStr)//": y_centroid must be set")
        @:PROHIBIT(f_is_default(patch_icpp(patch_id)%z_centroid), "Cuboid patch "//trim(iStr)//": z_centroid must be set")
        @:PROHIBIT(patch_icpp(patch_id)%length_x <= 0._wp, "Cuboid patch "//trim(iStr)//": length_x must be greater than zero")
        @:PROHIBIT(patch_icpp(patch_id)%length_y <= 0._wp, "Cuboid patch "//trim(iStr)//": length_y must be greater than zero")
        @:PROHIBIT(patch_icpp(patch_id)%length_z <= 0._wp, "Cuboid patch "//trim(iStr)//": length_z must be greater than zero")

    end subroutine s_check_cuboid_patch_geometry

    !>  This subroutine checks the model patch input
        !!  @param patch_id Patch identifier
    impure subroutine s_check_cylinder_patch_geometry(patch_id)

        ! Patch identifier
        integer, intent(in) :: patch_id
        call s_int_to_str(patch_id, iStr)

        @:PROHIBIT(p == 0, "Cylinder patch "//trim(iStr)//": p must be greater than zero")
        @:PROHIBIT(f_is_default(patch_icpp(patch_id)%x_centroid), "Cylinder patch "//trim(iStr)//": x_centroid must be set")
        @:PROHIBIT(f_is_default(patch_icpp(patch_id)%y_centroid), "Cylinder patch "//trim(iStr)//": y_centroid must be set")
        @:PROHIBIT(f_is_default(patch_icpp(patch_id)%z_centroid), "Cylinder patch "//trim(iStr)//": z_centroid must be set")
        @:PROHIBIT(patch_icpp(patch_id)%radius <= 0._wp, "Cylinder patch "//trim(iStr)//": radius must be greater than zero")

        ! Check if exactly one length is defined
        @:PROHIBIT(count([ &
            patch_icpp(patch_id)%length_x > 0._wp, &
            patch_icpp(patch_id)%length_y > 0._wp, &
            patch_icpp(patch_id)%length_z > 0._wp &
            ]) /= 1, "Cylinder patch "//trim(iStr)//": Exactly one of length_x, length_y, or length_z must be defined and positive")

        ! Ensure the defined length is positive
        @:PROHIBIT( &
            (.not. f_is_default(patch_icpp(patch_id)%length_x) .and. patch_icpp(patch_id)%length_x <= 0._wp) .or. &
            (.not. f_is_default(patch_icpp(patch_id)%length_y) .and. patch_icpp(patch_id)%length_y <= 0._wp) .or. &
            (.not. f_is_default(patch_icpp(patch_id)%length_z) .and. patch_icpp(patch_id)%length_z <= 0._wp), &
            "Cylinder patch "//trim(iStr)//": The defined length_{} must be greater than zero")

    end subroutine s_check_cylinder_patch_geometry

    !>  This subroutine checks the model patch input
        !!  @param patch_id Patch identifier
    impure subroutine s_check_plane_sweep_patch_geometry(patch_id)

        ! Patch identifier
        integer, intent(in) :: patch_id
        call s_int_to_str(patch_id, iStr)

        @:PROHIBIT(p == 0, "Plane sweep patch "//trim(iStr)//": p must be greater than zero")
        @:PROHIBIT(f_is_default(patch_icpp(patch_id)%x_centroid), "Plane sweep patch "//trim(iStr)//": x_centroid must be set")
        @:PROHIBIT(f_is_default(patch_icpp(patch_id)%y_centroid), "Plane sweep patch "//trim(iStr)//": y_centroid must be set")
        @:PROHIBIT(f_is_default(patch_icpp(patch_id)%z_centroid), "Plane sweep patch "//trim(iStr)//": z_centroid must be set")
        @:PROHIBIT(f_is_default(patch_icpp(patch_id)%normal(1)), "Plane sweep patch "//trim(iStr)//": normal(1) must be set")
        @:PROHIBIT(f_is_default(patch_icpp(patch_id)%normal(2)), "Plane sweep patch "//trim(iStr)//": normal(2) must be set")
        @:PROHIBIT(f_is_default(patch_icpp(patch_id)%normal(3)), "Plane sweep patch "//trim(iStr)//": normal(3) must be set")

    end subroutine s_check_plane_sweep_patch_geometry

    !> This subroutine checks the model patch input
        !!  @param patch_id Patch identifier
    impure subroutine s_check_ellipsoid_patch_geometry(patch_id)

        integer, intent(in) :: patch_id
        call s_int_to_str(patch_id, iStr)

        @:PROHIBIT(p == 0, "Ellipsoid patch "//trim(iStr)//": p must be greater than zero")
        @:PROHIBIT(f_is_default(patch_icpp(patch_id)%x_centroid), "Ellipsoid patch "//trim(iStr)//": x_centroid must be set")
        @:PROHIBIT(f_is_default(patch_icpp(patch_id)%y_centroid), "Ellipsoid patch "//trim(iStr)//": y_centroid must be set")
        @:PROHIBIT(f_is_default(patch_icpp(patch_id)%z_centroid), "Ellipsoid patch "//trim(iStr)//": z_centroid must be set")
        @:PROHIBIT(patch_icpp(patch_id)%radii(1) <= 0._wp, "Ellipsoid patch "//trim(iStr)//": radii(1) must be greater than zero")
        @:PROHIBIT(patch_icpp(patch_id)%radii(2) <= 0._wp, "Ellipsoid patch "//trim(iStr)//": radii(2) must be greater than zero")
        @:PROHIBIT(patch_icpp(patch_id)%radii(3) <= 0._wp, "Ellipsoid patch "//trim(iStr)//": radii(3) must be greater than zero")

    end subroutine s_check_ellipsoid_patch_geometry

    !!>  This subroutine verifies that the geometric parameters of
        !!      the inactive patch remain unaltered by the user inputs.
        !!  @param patch_id Patch identifier
    impure subroutine s_check_inactive_patch_geometry(patch_id)

        integer, intent(in) :: patch_id
        call s_int_to_str(patch_id, iStr)

        @:PROHIBIT(.not. f_is_default(patch_icpp(patch_id)%x_centroid), "Inactive patch "//trim(iStr)//": x_centroid must not be set")
        @:PROHIBIT(.not. f_is_default(patch_icpp(patch_id)%y_centroid), "Inactive patch "//trim(iStr)//": y_centroid must not be set")
        @:PROHIBIT(.not. f_is_default(patch_icpp(patch_id)%z_centroid), "Inactive patch "//trim(iStr)//": z_centroid must not be set")
        @:PROHIBIT(.not. f_is_default(patch_icpp(patch_id)%length_x), "Inactive patch "//trim(iStr)//": length_x must not be set")
        @:PROHIBIT(.not. f_is_default(patch_icpp(patch_id)%length_y), "Inactive patch "//trim(iStr)//": length_y must not be set")
        @:PROHIBIT(.not. f_is_default(patch_icpp(patch_id)%length_z), "Inactive patch "//trim(iStr)//": length_z must not be set")
        @:PROHIBIT(.not. f_is_default(patch_icpp(patch_id)%radius), "Inactive patch "//trim(iStr)//": radius must not be set")
        @:PROHIBIT(.not. f_is_default(patch_icpp(patch_id)%epsilon), "Inactive patch "//trim(iStr)//": epsilon must not be set")
        @:PROHIBIT(.not. f_is_default(patch_icpp(patch_id)%beta), "Inactive patch "//trim(iStr)//": beta must not be set")
        @:PROHIBIT(.not. f_is_default(patch_icpp(patch_id)%normal(1)), "Inactive patch "//trim(iStr)//": normal(1) must not be set")
        @:PROHIBIT(.not. f_is_default(patch_icpp(patch_id)%normal(2)), "Inactive patch "//trim(iStr)//": normal(2) must not be set")
        @:PROHIBIT(.not. f_is_default(patch_icpp(patch_id)%normal(3)), "Inactive patch "//trim(iStr)//": normal(3) must not be set")
        @:PROHIBIT(.not. f_is_default(patch_icpp(patch_id)%radii(1)), "Inactive patch "//trim(iStr)//": radii(1) must not be set")
        @:PROHIBIT(.not. f_is_default(patch_icpp(patch_id)%radii(2)), "Inactive patch "//trim(iStr)//": radii(2) must not be set")
        @:PROHIBIT(.not. f_is_default(patch_icpp(patch_id)%radii(3)), "Inactive patch "//trim(iStr)//": radii(3) must not be set")

    end subroutine s_check_inactive_patch_geometry

    !>  This subroutine verifies the active patch's right to overwrite the preceding patches
        !!  @param patch_id Patch identifier
    impure subroutine s_check_active_patch_alteration_rights(patch_id)

        integer, intent(in) :: patch_id
        call s_int_to_str(patch_id, iStr)

        @:PROHIBIT(.not. patch_icpp(patch_id)%alter_patch(0), "Patch "//trim(iStr)//": alter_patch(0) must be true")
        @:PROHIBIT(any(patch_icpp(patch_id)%alter_patch(patch_id:)), "Patch "//trim(iStr)// &
            ":alter_patch(i) must be false for i >= "//trim(iStr)//". Only preceding patches can be altered")

    end subroutine s_check_active_patch_alteration_rights

    !>  This subroutine verifies that inactive patches cannot overwrite other patches
        !!  @param patch_id Patch identifier
    impure subroutine s_check_inactive_patch_alteration_rights(patch_id)

        ! Patch identifier
        integer, intent(in) :: patch_id
        call s_int_to_str(patch_id, iStr)

        @:PROHIBIT(.not. patch_icpp(patch_id)%alter_patch(0), "Inactive patch "//trim(iStr)//": cannot have alter_patch(0) altered")
        @:PROHIBIT(any(patch_icpp(patch_id)%alter_patch(1:)), "Inactive patch "//trim(iStr)//": cannot have any alter_patch(i) enabled")

    end subroutine s_check_inactive_patch_alteration_rights

    !> This subroutine checks the smoothing parameters
        !!  @param patch_id Patch identifier
    impure subroutine s_check_supported_patch_smoothing(patch_id)

        integer, intent(in) :: patch_id
        call s_int_to_str(patch_id, iStr)

        if (patch_icpp(patch_id)%smoothen) then
            @:PROHIBIT(patch_icpp(patch_id)%smooth_patch_id >= patch_id, &
                "Smoothen enabled. Patch "//trim(iStr)//": smooth_patch_id must be less than patch_id")
            @:PROHIBIT(patch_icpp(patch_id)%smooth_patch_id == 0, &
                "Smoothen enabled. Patch "//trim(iStr)//": smooth_patch_id must be greater than zero")
            @:PROHIBIT(patch_icpp(patch_id)%smooth_coeff <= 0._wp, &
                "Smoothen enabled. Patch "//trim(iStr)//": smooth_coeff must be greater than zero")
        else
            @:PROHIBIT(patch_icpp(patch_id)%smooth_patch_id /= patch_id, &
                "Smoothen disabled. Patch "//trim(iStr)//": smooth_patch_id must be equal to patch_id")
            @:PROHIBIT(.not. f_is_default(patch_icpp(patch_id)%smooth_coeff), &
                "Smoothen disabled. Patch "//trim(iStr)//": smooth_coeff must not be set")
        end if

    end subroutine s_check_supported_patch_smoothing

    !> This subroutine verifies that inactive patches cannot be smoothed
        !!  @param patch_id Patch identifier
    impure subroutine s_check_unsupported_patch_smoothing(patch_id)

        ! Patch identifier
        integer, intent(in) :: patch_id
        call s_int_to_str(patch_id, iStr)

        @:PROHIBIT(patch_icpp(patch_id)%smoothen, &
            "Inactive patch "//trim(iStr)//": cannot have smoothen enabled")
        @:PROHIBIT(patch_icpp(patch_id)%smooth_patch_id /= patch_id, &
            "Inactive patch "//trim(iStr)//": smooth_patch_id must be equal to patch_id")
        @:PROHIBIT(.not. f_is_default(patch_icpp(patch_id)%smooth_coeff), &
            "Inactive patch "//trim(iStr)//": smooth_coeff must not be set")

    end subroutine s_check_unsupported_patch_smoothing

    !>  This subroutine checks the primitive variables
        !!  @param patch_id Patch identifier
    impure subroutine s_check_active_patch_primitive_variables(patch_id)

        integer, intent(in) :: patch_id

        logical, dimension(3) :: is_set_B

        call s_int_to_str(patch_id, iStr)

        @:PROHIBIT(f_is_default(patch_icpp(patch_id)%vel(1)), &
            "Patch "//trim(iStr)//": vel(1) must be set")
        @:PROHIBIT(n == 0 .and. (.not. f_is_default(patch_icpp(patch_id)%vel(2))) .and. (.not. f_approx_equal(patch_icpp(patch_id)%vel(2) , 0._wp)) .and. (.not. mhd), &
            "Patch "//trim(iStr)//": vel(2) must not be set when n = 0")
        @:PROHIBIT(n > 0 .and. f_is_default(patch_icpp(patch_id)%vel(2)), &
            "Patch "//trim(iStr)//": vel(2) must be set when n > 0")
        @:PROHIBIT(p == 0 .and. (.not. f_is_default(patch_icpp(patch_id)%vel(3))) .and. (.not. f_approx_equal(patch_icpp(patch_id)%vel(3) , 0._wp)) .and. (.not. mhd), &
            "Patch "//trim(iStr)//": vel(3) must not be set when p = 0")
        @:PROHIBIT(p > 0 .and. f_is_default(patch_icpp(patch_id)%vel(3)), &
            "Patch "//trim(iStr)//": vel(3) must be set when p > 0")
        @:PROHIBIT(mhd .and. (f_is_default(patch_icpp(patch_id)%vel(2)) .or. f_is_default(patch_icpp(patch_id)%vel(3))), &
            "Patch "//trim(iStr)//": All velocities (vel(1:3)) must be set when mhd = true")
        @:PROHIBIT(model_eqns == 1 .and. patch_icpp(patch_id)%rho <= 0._wp, &
            "Patch "//trim(iStr)//": rho must be greater than zero when model_eqns = 1")
        @:PROHIBIT(model_eqns == 1 .and. patch_icpp(patch_id)%gamma <= 0._wp, &
            "Patch "//trim(iStr)//": gamma must be greater than zero when model_eqns = 1")
        @:PROHIBIT(model_eqns == 1 .and. patch_icpp(patch_id)%pi_inf < 0._wp, &
            "Patch "//trim(iStr)//": pi_inf must be greater than or equal to zero when model_eqns = 1")
        @:PROHIBIT(patch_icpp(patch_id)%geometry == 5 .and. patch_icpp(patch_id)%pi_inf > 0, &
            "Patch "//trim(iStr)//": pi_inf must be less than or equal to zero when geometry = 5")
        @:PROHIBIT(model_eqns == 2 .and. any(patch_icpp(patch_id)%alpha_rho(1:num_fluids) < 0._wp), &
            "Patch "//trim(iStr)//": alpha_rho(1:num_fluids) must be greater than or equal to zero when model_eqns = 2")

        is_set_B(1) = .not. f_is_default(patch_icpp(patch_id)%Bx)
        is_set_B(2) = .not. f_is_default(patch_icpp(patch_id)%By)
        is_set_B(3) = .not. f_is_default(patch_icpp(patch_id)%Bz)

        @:PROHIBIT(.not. mhd .and. any(is_set_B), &
            "Bx, By, and Bz must not be set if MHD is not enabled")
        @:PROHIBIT(mhd .and. n == 0 .and. is_set_B(1), "Bx must not be set in 1D MHD simulations")
        @:PROHIBIT(mhd .and. n > 0 .and. .not. is_set_B(1), "Bx must be set in 2D/3D MHD simulations")
        @:PROHIBIT(mhd .and. .not. (is_set_B(2) .and. is_set_B(3)), "By and Bz must be set in all MHD simulations")

        if (model_eqns == 2 .and. num_fluids < num_fluids_max) then
            @:PROHIBIT(.not. f_all_default(patch_icpp(patch_id)%alpha_rho(num_fluids + 1:)), &
                "Patch "//trim(iStr)//": alpha_rho(i) must not be set for i > num_fluids")
            @:PROHIBIT(.not. f_all_default(patch_icpp(patch_id)%alpha(num_fluids + 1:)), &
                "Patch "//trim(iStr)//": alpha(i) must not be set for i > num_fluids")
            @:PROHIBIT(f_is_default(patch_icpp(patch_id)%alpha(num_fluids)), &
                "Patch "//trim(iStr)//": alpha(num_fluids) must be set")
        end if

        if (chemistry) then
            !@:ASSERT(all(patch_icpp(patch_id)%Y(1:num_species) >=       0._wp), "Patch " // trim(iStr) // ".")
            !@:ASSERT(any(patch_icpp(patch_id)%Y(1:num_species) >  verysmall), "Patch " // trim(iStr) // ".")
        end if

    end subroutine s_check_active_patch_primitive_variables

    !>  This subroutine verifies that the primitive variables
        !!      associated with the given inactive patch remain unaltered
        !!      by the user inputs.
        !!  @param patch_id Patch identifier
    impure subroutine s_check_inactive_patch_primitive_variables(patch_id)

        integer, intent(in) :: patch_id
        call s_int_to_str(patch_id, iStr)

        @:PROHIBIT(.not. f_all_default(patch_icpp(patch_id)%alpha_rho), &
            "Inactive patch "//trim(iStr)//": alpha_rho must not be set")
        @:PROHIBIT(.not. f_is_default(patch_icpp(patch_id)%rho), &
            "Inactive patch "//trim(iStr)//": rho must not be set")
        @:PROHIBIT(.not. f_all_default(patch_icpp(patch_id)%vel), &
            "Inactive patch "//trim(iStr)//": vel must not be set")
        @:PROHIBIT(.not. f_is_default(patch_icpp(patch_id)%pres), &
            "Inactive patch "//trim(iStr)//": pres must not be set")
        @:PROHIBIT(.not. f_all_default(patch_icpp(patch_id)%alpha), &
            "Inactive patch "//trim(iStr)//": alpha must not be set")
        @:PROHIBIT(.not. f_is_default(patch_icpp(patch_id)%gamma), &
            "Inactive patch "//trim(iStr)//": gamma must not be set")
        @:PROHIBIT(.not. f_is_default(patch_icpp(patch_id)%pi_inf), &
            "Inactive patch "//trim(iStr)//": pi_inf must not be set")

    end subroutine s_check_inactive_patch_primitive_variables

    impure subroutine s_check_model_geometry(patch_id)

        integer, intent(in) :: patch_id

        logical :: file_exists

        inquire (file=patch_icpp(patch_id)%model_filepath, exist=file_exists)

        @:PROHIBIT(.not. file_exists, "Model file "//trim(patch_icpp(patch_id)%model_filepath)// &
            " requested by patch "//trim(iStr)//" does not exist")

    end subroutine s_check_model_geometry

end module m_check_patches
