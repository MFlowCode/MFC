!> @brief This module contains subroutines that read, and check consistency
!!              of, the user provided inputs and data.

#:include 'macros.fpp'

module m_check_ib_patches

    ! Dependencies =============================================================
    use m_derived_types          !< Definitions of the derived types

    use m_global_parameters      !< Global parameters

    use m_mpi_proxy              !< Message passing interface (MPI) module proxy

    use m_data_output            !< Procedures to write the grid data and the
                                 !! conservative variables to files

#ifdef MFC_MPI
    use mpi                      !< Message passing interface (MPI) module
#endif

    use m_compile_specific

    use m_helper_basic           !< Functions to compare floating point numbers

    use m_helper
    ! ==========================================================================

    implicit none

    private; 
    public :: s_check_ib_patches

    character(len=10) :: iStr

contains

    subroutine s_check_ib_patches

        integer :: i

        do i = 1, num_patches_max
            if (i <= num_ibs) then
                ! call s_check_patch_geometry(i)
                call s_int_to_str(i, iStr)

                ! Constraints on the geometric initial condition patch parameters
                if (patch_ib(i)%geometry == 2) then
                    call s_check_circle_ib_patch_geometry(i)
                else if (patch_ib(i)%geometry == 3) then
                    call s_check_rectangle_ib_patch_geometry(i)
                else if (patch_ib(i)%geometry == 8) then
                    call s_check_sphere_ib_patch_geometry(i)
                else if (patch_ib(i)%geometry == 4) then
                    call s_check_airfoil_ib_patch_geometry(i)
                else if (patch_ib(i)%geometry == 11) then
                    call s_check_3d_airfoil_ib_patch_geometry(i)
                else if (patch_ib(i)%geometry == 10) then
                    call s_check_cylinder_ib_patch_geometry(i)
                else if (patch_ib(i)%geometry == dflt_int) then
                    call s_prohibit_abort("IB patch undefined", &
                                          "patch_ib("//trim(iStr)//")%geometry must be set.")
                else
                    call s_prohibit_abort("Invalid IB patch", &
                                          "patch_ib("//trim(iStr)//")%geometry must be "// &
                                          "2-4, 8, 10, or 11.")
                end if
            else
                if (patch_ib(i)%geometry == dflt_int) then
                    call s_check_inactive_ib_patch_geometry(i)
                else
                    call s_prohibit_abort("Inactive IB patch defined", &
                                          "patch_ib("//trim(iStr)//")%geometry "// &
                                          "must not be set for inactive patches.")
                end if
            end if
        end do

    end subroutine s_check_ib_patches

    !>  This subroutine verifies that the geometric parameters of
        !!      the circle patch have consistently been inputted by the
        !!      user.
        !!  @param patch_id Patch identifier
    subroutine s_check_circle_ib_patch_geometry(patch_id)

        integer, intent(in) :: patch_id

        call s_int_to_str(patch_id, iStr)

        @:PROHIBIT(n == 0 .or. p > 0 &
            .or. patch_ib(patch_id)%radius <= 0d0 &
            .or. f_is_default(patch_ib(patch_id)%x_centroid) &
            .or. f_is_default(patch_ib(patch_id)%y_centroid), &
            'in circle IB patch '//trim(iStr))

    end subroutine s_check_circle_ib_patch_geometry

    !>  This subroutine verifies that the geometric parameters of
        !!      the airfoil patch have consistently been inputted by the
        !!      user.
        !!  @param patch_id Patch identifier
    subroutine s_check_airfoil_ib_patch_geometry(patch_id)

        integer, intent(in) :: patch_id

        call s_int_to_str(patch_id, iStr)

        @:PROHIBIT(n == 0 .or. p > 0 &
            .or. patch_ib(patch_id)%c <= 0d0 &
            .or. patch_ib(patch_id)%p <= 0d0 &
            .or. patch_ib(patch_id)%t <= 0d0 &
            .or. patch_ib(patch_id)%m <= 0d0 &
            .or. f_is_default(patch_ib(patch_id)%x_centroid) &
            .or. f_is_default(patch_ib(patch_id)%y_centroid), &
            'in airfoil IB patch '//trim(iStr))

    end subroutine s_check_airfoil_ib_patch_geometry

    !>  This subroutine verifies that the geometric parameters of
        !!      the 3d airfoil patch have consistently been inputted by the
        !!      user.
        !!  @param patch_id Patch identifier
    subroutine s_check_3d_airfoil_ib_patch_geometry(patch_id)

        integer, intent(in) :: patch_id

        call s_int_to_str(patch_id, iStr)

        @:PROHIBIT(n == 0 .or. p == 0 &
            .or. patch_ib(patch_id)%c <= 0d0 &
            .or. patch_ib(patch_id)%p <= 0d0 &
            .or. patch_ib(patch_id)%t <= 0d0 &
            .or. patch_ib(patch_id)%m <= 0d0 &
            .or. f_is_default(patch_ib(patch_id)%x_centroid) &
            .or. f_is_default(patch_ib(patch_id)%y_centroid) &
            .or. f_is_default(patch_ib(patch_id)%z_centroid) &
            .or. f_is_default(patch_ib(patch_id)%length_z), &
            'in 3d airfoil IB patch '//trim(iStr))

    end subroutine s_check_3d_airfoil_ib_patch_geometry

    !>  This subroutine verifies that the geometric parameters of
        !!      the rectangle patch have consistently been inputted by
        !!      the user.
        !!  @param patch_id Patch identifier
    subroutine s_check_rectangle_ib_patch_geometry(patch_id)

        integer, intent(in) :: patch_id

        call s_int_to_str(patch_id, iStr)

        @:PROHIBIT(n == 0 .or. p > 0 &
            .or. &
            f_is_default(patch_ib(patch_id)%x_centroid) &
            .or. &
            f_is_default(patch_ib(patch_id)%y_centroid) &
            .or. &
            patch_ib(patch_id)%length_x <= 0d0 &
            .or. &
            patch_ib(patch_id)%length_y <= 0d0, &
            'in rectangle IB patch '//trim(iStr))

    end subroutine s_check_rectangle_ib_patch_geometry

    !>  This subroutine verifies that the geometric parameters of
        !!      the sphere patch have consistently been inputted by
        !!      the user.
        !!  @param patch_id Patch identifier
    subroutine s_check_sphere_ib_patch_geometry(patch_id)

        integer, intent(in) :: patch_id

        call s_int_to_str(patch_id, iStr)

        @:PROHIBIT(n == 0 .or. p == 0 &
            .or. &
            f_is_default(patch_ib(patch_id)%x_centroid) &
            .or. &
            f_is_default(patch_ib(patch_id)%y_centroid) &
            .or. &
            f_is_default(patch_ib(patch_id)%z_centroid) &
            .or. &
            patch_ib(patch_id)%radius <= 0d0, &
            'in sphere IB patch '//trim(iStr))

    end subroutine s_check_sphere_ib_patch_geometry

    !>  This subroutine verifies that the geometric parameters of
        !!      the cylinder patch have consistently been inputted by
        !!      the user.
        !!  @param patch_id Patch identifier
    subroutine s_check_cylinder_ib_patch_geometry(patch_id)

        integer, intent(in) :: patch_id

        call s_int_to_str(patch_id, iStr)

        @:PROHIBIT(p == 0 &
            .or. &
            f_is_default(patch_ib(patch_id)%x_centroid) &
            .or. &
            f_is_default(patch_ib(patch_id)%y_centroid) &
            .or. &
            f_is_default(patch_ib(patch_id)%z_centroid) &
            .or. &
            (patch_ib(patch_id)%length_x <= 0d0 .and. &
             patch_ib(patch_id)%length_y <= 0d0 .and. &
             patch_ib(patch_id)%length_z <= 0d0) &
            .or. &
            patch_ib(patch_id)%radius <= 0d0, &
            'in cylinder IB patch '//trim(iStr))

        @:PROHIBIT( &
            (patch_ib(patch_id)%length_x > 0d0 .and. &
             ((.not. f_is_default(patch_ib(patch_id)%length_y)) .or. &
              (.not. f_is_default(patch_ib(patch_id)%length_z)))) &
            .or. &
            (patch_ib(patch_id)%length_y > 0d0 .and. &
             ((.not. f_is_default(patch_ib(patch_id)%length_x)) .or. &
              (.not. f_is_default(patch_ib(patch_id)%length_z)))) &
            .or. &
            (patch_ib(patch_id)%length_z > 0d0 .and. &
             ((.not. f_is_default(patch_ib(patch_id)%length_x)) .or. &
              (.not. f_is_default(patch_ib(patch_id)%length_y)))), &
            'in cylinder IB patch '//trim(iStr))

    end subroutine s_check_cylinder_ib_patch_geometry

    !!>  This subroutine verifies that the geometric parameters of
        !!      the inactive patch remain unaltered by the user inputs.
        !!  @param patch_id Patch identifier
    subroutine s_check_inactive_ib_patch_geometry(patch_id)

        integer, intent(in) :: patch_id

        call s_int_to_str(patch_id, iStr)

        @:PROHIBIT((.not. f_is_default(patch_ib(patch_id)%x_centroid)) &
            .or. &
            (.not. f_is_default(patch_ib(patch_id)%y_centroid)) &
            .or. &
            (.not. f_is_default(patch_ib(patch_id)%z_centroid)) &
            .or. &
            (.not. f_is_default(patch_ib(patch_id)%length_x)) &
            .or. &
            (.not. f_is_default(patch_ib(patch_id)%length_y)) &
            .or. &
            (.not. f_is_default(patch_ib(patch_id)%length_z)) &
            .or. &
            (.not. f_is_default(patch_ib(patch_id)%radius)), &
            'in inactive IB patch '//trim(iStr))

    end subroutine s_check_inactive_ib_patch_geometry

end module m_check_ib_patches
