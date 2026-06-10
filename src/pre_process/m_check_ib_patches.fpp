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

                ! Constraints on the geometric initial condition patch parameters
                if (patch_ib(i)%geometry == 2 .or. patch_ib(i)%geometry == 8 .or. patch_ib(i)%geometry == 10) then
                    call s_check_circle_ib_patch_geometry(i)
                else if (patch_ib(i)%geometry == 3 .or. patch_ib(i)%geometry == 9) then
                    call s_check_rectangle_ib_patch_geometry(i)
                else if (patch_ib(i)%geometry == 4 .or. patch_ib(i)%geometry == 11) then
                    call s_check_airfoil_ib_patch_geometry(i)
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

    !> Verify that the geometric parameters of the circle patch have been consistently inputted.

    impure subroutine s_check_circle_ib_patch_geometry(patch_id)

        integer, intent(in) :: patch_id

        call s_int_to_str(patch_id, iStr)

        @:PROHIBIT(n == 0 .or. patch_ib(patch_id)%radius <= 0._wp .or. f_is_default(patch_ib(patch_id)%x_centroid) &
                   & .or. f_is_default(patch_ib(patch_id)%y_centroid), 'in circle/sphere/cylinder IB patch ' // trim(iStr))

        ! Spheres and cylinders (3D) additionally require a z-centroid
        if (p > 0) then
            @:PROHIBIT(f_is_default(patch_ib(patch_id)%z_centroid), 'in 3D sphere/cylinder IB patch ' // trim(iStr))
        end if

        ! Cylinders are extruded along exactly one positive length
        if (patch_ib(patch_id)%geometry == 10) then
            @:PROHIBIT(count([patch_ib(patch_id)%length_x > 0._wp, patch_ib(patch_id)%length_y > 0._wp, &
                       & patch_ib(patch_id)%length_z > 0._wp]) /= 1, &
                       & 'in cylinder IB patch ' // trim(iStr) &
                       & // ': exactly one of length_x, length_y, or length_z must be defined and positive')
            @:PROHIBIT((.not. f_is_default(patch_ib(patch_id)%length_x) .and. patch_ib(patch_id)%length_x <= 0._wp) &
                       & .or. (.not. f_is_default(patch_ib(patch_id)%length_y) .and. patch_ib(patch_id)%length_y <= 0._wp) &
                       & .or. (.not. f_is_default(patch_ib(patch_id)%length_z) .and. patch_ib(patch_id)%length_z <= 0._wp), &
                       & 'in cylinder IB patch ' // trim(iStr) // ': the defined lengths must be greater than zero')
        end if

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

        @:PROHIBIT(n == 0 .or. patch_ib(patch_id)%airfoil_id <= 0 .or. ib_airfoil(patch_ib(patch_id)%airfoil_id)%c <= 0._wp &
                   & .or. ib_airfoil(patch_ib(patch_id)%airfoil_id)%p <= 0._wp &
                   & .or. ib_airfoil(patch_ib(patch_id)%airfoil_id)%t <= 0._wp &
                   & .or. ib_airfoil(patch_ib(patch_id)%airfoil_id)%m <= 0._wp .or. f_is_default(patch_ib(patch_id)%x_centroid) &
                   & .or. f_is_default(patch_ib(patch_id)%y_centroid), 'in airfoil IB patch ' // trim(iStr))

        ! Additional checks strictly for 3D airfoils
        if (p > 0) then
            @:PROHIBIT(f_is_default(patch_ib(patch_id)%z_centroid) .or. f_is_default(patch_ib(patch_id)%length_z), &
                       & 'in 3D airfoil IB patch ' // trim(iStr))
        end if

    end subroutine s_check_airfoil_ib_patch_geometry

    !> Verify that the geometric parameters of the rectangle patch have been consistently inputted.

    impure subroutine s_check_rectangle_ib_patch_geometry(patch_id)

        integer, intent(in) :: patch_id

        call s_int_to_str(patch_id, iStr)

        @:PROHIBIT(n == 0 .or. f_is_default(patch_ib(patch_id)%x_centroid) .or. f_is_default(patch_ib(patch_id)%y_centroid) &
                   & .or. patch_ib(patch_id)%length_x <= 0._wp .or. patch_ib(patch_id)%length_y <= 0._wp, &
                   & 'in rectangle/cuboid IB patch ' // trim(iStr))

        ! If the simulation is 3D, also mandate Z lengths and centroids
        if (p > 0) then
            @:PROHIBIT(f_is_default(patch_ib(patch_id)%z_centroid) .or. patch_ib(patch_id)%length_z <= 0._wp, &
                       & 'in 3D cuboid IB patch ' // trim(iStr))
        end if

    end subroutine s_check_rectangle_ib_patch_geometry

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
                   & .or. (.not. f_is_default(patch_ib(patch_id)%length_z)) .or. (.not. f_is_default(patch_ib(patch_id)%radius)), &
                   & 'in inactive IB patch ' // trim(iStr))

    end subroutine s_check_inactive_ib_patch_geometry

end module m_check_ib_patches
