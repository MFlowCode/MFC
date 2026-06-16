!>
!! @file
!! @brief Contains module m_check_ib_patches

!> @brief Validates geometry parameters and constraints for immersed boundary patches

#:include 'macros.fpp'

module m_check_ib_patches

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_data_output
#ifdef MFC_MPI
    use mpi  !< Message passing interface (MPI) module
#endif

    use m_compile_specific
    use m_helper_basic
    use m_helper

    implicit none

    private
    public :: s_check_ib_patches

    character(len=10) :: iStr

contains

    !> Validate the geometry parameters of all active and inactive immersed boundary patches.
    impure subroutine s_check_ib_patches

        integer :: i

        do i = 1, num_patches_max
            if (i <= num_ibs) then
                call s_int_to_str(i, iStr)
                @:PROHIBIT(patch_ib(i)%geometry == dflt_int, "IB patch undefined. patch_ib("//trim(iStr)//")%geometry must be set.")
                call s_check_ib_patch_thermal(i)

                ! Constraints on the geometric initial condition patch parameters
                if (patch_ib(i)%geometry == 2) then
                    call s_check_circle_ib_patch_geometry(i)
                else if (patch_ib(i)%geometry == 3) then
                    call s_check_rectangle_ib_patch_geometry(i)
                else if (patch_ib(i)%geometry == 8) then
                    call s_check_sphere_ib_patch_geometry(i)
                else if (patch_ib(i)%geometry == 9) then
                    call s_check_cuboid_ib_patch_geometry(i)
                else if (patch_ib(i)%geometry == 4) then
                    call s_check_airfoil_ib_patch_geometry(i)
                else if (patch_ib(i)%geometry == 11) then
                    call s_check_3d_airfoil_ib_patch_geometry(i)
                else if (patch_ib(i)%geometry == 10) then
                    call s_check_cylinder_ib_patch_geometry(i)
                else if (patch_ib(i)%geometry == 5 .or. patch_ib(i)%geometry == 12) then
                    call s_check_model_ib_patch_geometry(i)
                else if (patch_ib(i)%geometry == 6) then
                    call s_check_ellipse_ib_patch_geometry(i)
                else
                    call s_prohibit_abort("Invalid IB patch", &
                                          & "patch_ib(" // trim(iStr) // ")%geometry must be " // "2-4, 8-10, 11 or 12.")
                end if
            else
                @:PROHIBIT(patch_ib(i)%geometry /= dflt_int, &
                           & "Inactive IB patch defined. " // "patch_ib(" // trim(iStr) &
                           & // ")%geometry must not be set for inactive patches.")
                call s_check_inactive_ib_patch_geometry(i)
            end if
        end do

    end subroutine s_check_ib_patches

    impure subroutine s_check_ib_patch_thermal(patch_id)

        integer, intent(in) :: patch_id

        call s_int_to_str(patch_id, iStr)

        if (.not. patch_ib(patch_id)%isothermal) return

        @:PROHIBIT(f_is_default(patch_ib(patch_id)%Twall) .or. patch_ib(patch_id)%Twall <= 0._wp, &
                   & 'in immersed boundary patch ' // trim(iStr) // ': patch_ib('//trim(iStr)//')%Twall must be set and > 0')

    end subroutine s_check_ib_patch_thermal

    !> Verify that the geometric parameters of the circle patch have been consistently inputted.

    impure subroutine s_check_circle_ib_patch_geometry(patch_id)

        integer, intent(in) :: patch_id

        call s_int_to_str(patch_id, iStr)

        @:PROHIBIT(n == 0 .or. p > 0 .or. patch_ib(patch_id)%radius <= 0._wp .or. f_is_default(patch_ib(patch_id)%x_centroid) &
                   & .or. f_is_default(patch_ib(patch_id)%y_centroid), 'in circle IB patch ' // trim(iStr))

    end subroutine s_check_circle_ib_patch_geometry

    !> Verify that the geometric parameters of the ellipse patch have been consistently inputted.

    impure subroutine s_check_ellipse_ib_patch_geometry(patch_id)

        integer, intent(in) :: patch_id

        call s_int_to_str(patch_id, iStr)

        @:PROHIBIT(n == 0 .or. p > 0 .or. patch_ib(patch_id)%length_x <= 0._wp .or. patch_ib(patch_id)%length_y <= 0._wp &
                   & .or. f_is_default(patch_ib(patch_id)%x_centroid) .or. f_is_default(patch_ib(patch_id)%y_centroid), &
                   & 'in ellipse IB patch ' // trim(iStr))

    end subroutine s_check_ellipse_ib_patch_geometry

    !> Verify that the geometric parameters of the airfoil patch have been consistently inputted.

    impure subroutine s_check_airfoil_ib_patch_geometry(patch_id)

        integer, intent(in) :: patch_id

        call s_int_to_str(patch_id, iStr)

        @:PROHIBIT(n == 0 .or. p > 0 .or. patch_ib(patch_id)%airfoil_id <= 0 &
                   & .or. ib_airfoil(patch_ib(patch_id)%airfoil_id)%c <= 0._wp &
                   & .or. ib_airfoil(patch_ib(patch_id)%airfoil_id)%p <= 0._wp &
                   & .or. ib_airfoil(patch_ib(patch_id)%airfoil_id)%t <= 0._wp &
                   & .or. ib_airfoil(patch_ib(patch_id)%airfoil_id)%m <= 0._wp .or. f_is_default(patch_ib(patch_id)%x_centroid) &
                   & .or. f_is_default(patch_ib(patch_id)%y_centroid), 'in airfoil IB patch ' // trim(iStr))

    end subroutine s_check_airfoil_ib_patch_geometry

    !> Verify that the geometric parameters of the 3D airfoil patch have been consistently inputted.

    impure subroutine s_check_3d_airfoil_ib_patch_geometry(patch_id)

        integer, intent(in) :: patch_id

        call s_int_to_str(patch_id, iStr)

        @:PROHIBIT(n == 0 .or. p == 0 .or. patch_ib(patch_id)%airfoil_id <= 0 &
                   & .or. ib_airfoil(patch_ib(patch_id)%airfoil_id)%c <= 0._wp &
                   & .or. ib_airfoil(patch_ib(patch_id)%airfoil_id)%p <= 0._wp &
                   & .or. ib_airfoil(patch_ib(patch_id)%airfoil_id)%t <= 0._wp &
                   & .or. ib_airfoil(patch_ib(patch_id)%airfoil_id)%m <= 0._wp .or. f_is_default(patch_ib(patch_id)%x_centroid) &
                   & .or. f_is_default(patch_ib(patch_id)%y_centroid) .or. f_is_default(patch_ib(patch_id)%z_centroid) &
                   & .or. f_is_default(patch_ib(patch_id)%length_z), 'in 3d airfoil IB patch ' // trim(iStr))

    end subroutine s_check_3d_airfoil_ib_patch_geometry

    !> Verify that the geometric parameters of the rectangle patch have been consistently inputted.

    impure subroutine s_check_rectangle_ib_patch_geometry(patch_id)

        integer, intent(in) :: patch_id

        call s_int_to_str(patch_id, iStr)

        @:PROHIBIT(n == 0 .or. p > 0 .or. f_is_default(patch_ib(patch_id)%x_centroid) &
                   & .or. f_is_default(patch_ib(patch_id)%y_centroid) .or. patch_ib(patch_id)%length_x <= 0._wp &
                   & .or. patch_ib(patch_id)%length_y <= 0._wp, 'in rectangle IB patch ' // trim(iStr))

    end subroutine s_check_rectangle_ib_patch_geometry

    !> Verify that the geometric parameters of the sphere patch have been consistently inputted.

    impure subroutine s_check_sphere_ib_patch_geometry(patch_id)

        integer, intent(in) :: patch_id

        call s_int_to_str(patch_id, iStr)

        @:PROHIBIT(n == 0 .or. p == 0 .or. f_is_default(patch_ib(patch_id)%x_centroid) &
                   & .or. f_is_default(patch_ib(patch_id)%y_centroid) .or. f_is_default(patch_ib(patch_id)%z_centroid) &
                   & .or. patch_ib(patch_id)%radius <= 0._wp, 'in sphere IB patch ' // trim(iStr))

    end subroutine s_check_sphere_ib_patch_geometry

    !> Verify that the geometric parameters of the cuboid patch have been consistently inputted.

    impure subroutine s_check_cuboid_ib_patch_geometry(patch_id)

        integer, intent(in) :: patch_id

        call s_int_to_str(patch_id, iStr)

        @:PROHIBIT(n == 0 .or. p == 0 .or. f_is_default(patch_ib(patch_id)%x_centroid) &
                   & .or. f_is_default(patch_ib(patch_id)%y_centroid) .or. f_is_default(patch_ib(patch_id)%z_centroid) &
                   & .or. patch_ib(patch_id)%length_x <= 0._wp .or. patch_ib(patch_id)%length_y <= 0._wp &
                   & .or. patch_ib(patch_id)%length_z <= 0._wp, 'in cuboid IB patch ' // trim(iStr))

    end subroutine s_check_cuboid_ib_patch_geometry

    !> Verify that the geometric parameters of the cylinder patch have been consistently inputted.

    impure subroutine s_check_cylinder_ib_patch_geometry(patch_id)

        integer, intent(in) :: patch_id

        call s_int_to_str(patch_id, iStr)

        @:PROHIBIT(p == 0 .or. f_is_default(patch_ib(patch_id)%x_centroid) .or. f_is_default(patch_ib(patch_id)%y_centroid) &
                   & .or. f_is_default(patch_ib(patch_id)%z_centroid) .or. (patch_ib(patch_id)%length_x <= 0._wp &
                   & .and. patch_ib(patch_id)%length_y <= 0._wp .and. patch_ib(patch_id)%length_z <= 0._wp) &
                   & .or. patch_ib(patch_id)%radius <= 0._wp, 'in cylinder IB patch ' // trim(iStr))

        @:PROHIBIT((patch_ib(patch_id)%length_x > 0._wp .and. ((.not. f_is_default(patch_ib(patch_id)%length_y)) &
                   & .or. (.not. f_is_default(patch_ib(patch_id)%length_z)))) .or. (patch_ib(patch_id)%length_y > 0._wp &
                   & .and. ((.not. f_is_default(patch_ib(patch_id)%length_x)) &
                   & .or. (.not. f_is_default(patch_ib(patch_id)%length_z)))) .or. (patch_ib(patch_id)%length_z > 0._wp &
                   & .and. ((.not. f_is_default(patch_ib(patch_id)%length_x)) &
                   & .or. (.not. f_is_default(patch_ib(patch_id)%length_y)))), 'in cylinder IB patch ' // trim(iStr))

    end subroutine s_check_cylinder_ib_patch_geometry

    !> Verify that the geometric parameters of the model patch have been consistently inputted.

    impure subroutine s_check_model_ib_patch_geometry(patch_id)

        integer, intent(in) :: patch_id
        integer             :: mid
        character(len=10)   :: midStr

        call s_int_to_str(patch_id, iStr)

        mid = patch_ib(patch_id)%model_id
        call s_int_to_str(mid, midStr)

        @:PROHIBIT(mid <= 0 .or. mid > num_stl_models, &
                   & 'patch_ib('//trim(iStr)//')%model_id='//trim(midStr)//' must be in [1, num_stl_models]')

        @:PROHIBIT(stl_models(mid)%model_filepath == dflt_char, 'Empty model file path for stl_models('//trim(midStr)//')')

        @:PROHIBIT(stl_models(mid)%model_scale(1) <= 0._wp .or. stl_models(mid)%model_scale(2) <= 0._wp &
                   & .or. stl_models(mid)%model_scale(3) <= 0._wp, 'Negative scale in stl_models(' // trim(midStr) // ')')

    end subroutine s_check_model_ib_patch_geometry

    !> Verify that inactive IB patch geometry parameters remain at defaults
    impure subroutine s_check_inactive_ib_patch_geometry(patch_id)

        integer, intent(in) :: patch_id

        call s_int_to_str(patch_id, iStr)

        @:PROHIBIT((.not. f_is_default(patch_ib(patch_id)%x_centroid)) .or. (.not. f_is_default(patch_ib(patch_id)%y_centroid)) &
                   & .or. (.not. f_is_default(patch_ib(patch_id)%z_centroid)) &
                   & .or. (.not. f_is_default(patch_ib(patch_id)%length_x)) .or. (.not. f_is_default(patch_ib(patch_id)%length_y)) &
                   & .or. (.not. f_is_default(patch_ib(patch_id)%length_z)) .or. (.not. f_is_default(patch_ib(patch_id)%radius)) &
                   & .or. patch_ib(patch_id)%isothermal .or. (.not. f_is_default(patch_ib(patch_id)%Twall)), &
                   & 'in inactive IB patch ' // trim(iStr))

    end subroutine s_check_inactive_ib_patch_geometry

end module m_check_ib_patches
