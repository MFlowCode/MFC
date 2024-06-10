!> @brief This module contains subroutines that read, and check consistency
!!              of, the user provided inputs and data.
module m_check_ib_patches

    ! Dependencies =============================================================
    use m_derived_types          !< Definitions of the derived types

    use m_global_parameters      !< Global parameters for the code

    use m_mpi_proxy              !< Message passing interface (MPI) module proxy

    use m_data_output            !< Procedures to write the grid data and the
                                 !! conservative variables to files

#ifdef MFC_MPI
    use mpi                      !< Message passing interface (MPI) module
#endif

    use m_compile_specific

    use m_helper
    ! ==========================================================================

    implicit none

    private; public :: s_check_ib_patches

    character(len=10) :: iStr

contains

    subroutine s_check_ib_patches()

        ! integer, intent(in) :: i

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
                    call s_check_3D_airfoil_ib_patch_geometry(i)
                else if (patch_ib(i)%geometry == 10) then
                    call s_check_cylinder_ib_patch_geometry(i)
                else
                    call s_mpi_abort('Unsupported choice of the '// &
                                     'geometry of active patch '//trim(iStr)// &
                                     ' detected. Exiting ...')
                end if
            else
                if (patch_ib(i)%geometry == dflt_int) then
                    call s_check_inactive_ib_patch_geometry(i)
                else
                    call s_mpi_abort('Unsupported choice of the '// &
                                     'geometry of inactive patch '//trim(iStr)// &
                                     ' detected. Exiting ...')
                end if
            end if
        end do

    end subroutine s_check_ib_patches

    !>  This subroutine verifies that the geometric parameters of
        !!      the circle patch have consistently been inputted by the
        !!      user.
        !!  @param patch_id Patch identifier
    subroutine s_check_circle_ib_patch_geometry(patch_id) ! -------------------

        integer, intent(IN) :: patch_id
        call s_int_to_str(patch_id, iStr)

        ! Constraints on the geometric parameters of the circle patch
        if (n == 0 .or. p > 0 .or. patch_ib(patch_id)%radius <= 0d0 &
            .or. &
            patch_ib(patch_id)%x_centroid == dflt_real &
            .or. &
            patch_ib(patch_id)%y_centroid == dflt_real) then

            call s_mpi_abort('Inconsistency(ies) detected in '// &
                             'geometric parameters of circle '// &
                             'patch '//trim(iStr)//'. Exiting ...')

        end if

    end subroutine s_check_circle_ib_patch_geometry ! -------------------------

    subroutine s_check_airfoil_ib_patch_geometry(patch_id) ! -------------------

        integer, intent(IN) :: patch_id
        call s_int_to_str(patch_id, iStr)

        ! Constraints on the geometric parameters of the circle patch
        if (n == 0 .or. p > 0 .or. patch_ib(patch_id)%c <= 0d0 &
            .or. patch_ib(patch_id)%p <= 0d0 .or. patch_ib(patch_id)%t <= 0d0 &
            .or. patch_ib(patch_id)%m <= 0d0 .or. patch_ib(patch_id)%x_centroid == dflt_real &
            .or. patch_ib(patch_id)%y_centroid == dflt_real) then

            call s_mpi_abort('Inconsistency(ies) detected in '// &
                             'geometric parameters of airfoil '// &
                             'patch '//trim(iStr)//'. Exiting ...')

        end if

    end subroutine s_check_airfoil_ib_patch_geometry ! -------------------------

    subroutine s_check_3d_airfoil_ib_patch_geometry(patch_id) ! -------------------

        integer, intent(IN) :: patch_id
        call s_int_to_str(patch_id, iStr)

        ! Constraints on the geometric parameters of the circle patch
        if (n == 0 .or. p == 0 .or. patch_ib(patch_id)%c <= 0d0 &
            .or. patch_ib(patch_id)%p <= 0d0 .or. patch_ib(patch_id)%t <= 0d0 &
            .or. patch_ib(patch_id)%m <= 0d0 .or. patch_ib(patch_id)%x_centroid == dflt_real &
            .or. patch_ib(patch_id)%y_centroid == dflt_real .or. patch_ib(patch_id)%z_centroid == dflt_real &
            .or. patch_ib(patch_id)%length_z == dflt_real) then

            call s_mpi_abort('Inconsistency(ies) detected in '// &
                             'geometric parameters of airfoil '// &
                             'patch '//trim(iStr)//'. Exiting ...')

        end if

    end subroutine s_check_3d_airfoil_ib_patch_geometry ! -------------------------

    !>  This subroutine verifies that the geometric parameters of
        !!      the rectangle patch have consistently been inputted by
        !!      the user.
        !!  @param patch_id Patch identifier
    subroutine s_check_rectangle_ib_patch_geometry(patch_id) ! ----------------

        integer, intent(IN) :: patch_id
        call s_int_to_str(patch_id, iStr)

        ! Constraints on the geometric parameters of the rectangle patch
        if (n == 0 .or. p > 0 &
            .or. &
            patch_ib(patch_id)%x_centroid == dflt_real &
            .or. &
            patch_ib(patch_id)%y_centroid == dflt_real &
            .or. &
            patch_ib(patch_id)%length_x <= 0d0 &
            .or. &
            patch_ib(patch_id)%length_y <= 0d0) then

            call s_mpi_abort('Inconsistency(ies) detected in '// &
                             'geometric parameters of rectangle '// &
                             'patch '//trim(iStr)//'. Exiting ...')

        end if

    end subroutine s_check_rectangle_ib_patch_geometry ! ----------------------

    !>  This subroutine verifies that the geometric parameters of
        !!      the sphere patch have consistently been inputted by
        !!      the user.
        !!  @param patch_id Patch identifier
    subroutine s_check_sphere_ib_patch_geometry(patch_id) ! ----------------

        integer, intent(IN) :: patch_id
        call s_int_to_str(patch_id, iStr)

        ! Constraints on the geometric parameters of the rectangle patch
        if (n == 0 .or. p == 0 &
            .or. &
            patch_ib(patch_id)%x_centroid == dflt_real &
            .or. &
            patch_ib(patch_id)%y_centroid == dflt_real &
            .or. &
            patch_ib(patch_id)%z_centroid == dflt_real &
            .or. &
            patch_ib(patch_id)%radius <= 0d0) then

            call s_mpi_abort('Inconsistency(ies) detected in '// &
                             'geometric parameters of rectangle '// &
                             'patch '//trim(iStr)//'. Exiting ...')

        end if

    end subroutine s_check_sphere_ib_patch_geometry ! ----------------------

    !>  This subroutine verifies that the geometric parameters of
        !!      the cylinder patch have consistently been inputted by
        !!      the user.
        !!  @param patch_id Patch identifier
    subroutine s_check_cylinder_ib_patch_geometry(patch_id) ! -----------------

        ! Patch identifier
        integer, intent(IN) :: patch_id
        call s_int_to_str(patch_id, iStr)

        ! Constraints on the geometric parameters of the cylinder patch
        if (p == 0 &
            .or. &
            patch_ib(patch_id)%x_centroid == dflt_real &
            .or. &
            patch_ib(patch_id)%y_centroid == dflt_real &
            .or. &
            patch_ib(patch_id)%z_centroid == dflt_real &
            .or. &
            (patch_ib(patch_id)%length_x <= 0d0 .and. &
             patch_ib(patch_id)%length_y <= 0d0 .and. &
             patch_ib(patch_id)%length_z <= 0d0) &
            .or. &
            (patch_ib(patch_id)%length_x > 0d0 .and. &
             (patch_ib(patch_id)%length_y /= dflt_real .or. &
              patch_ib(patch_id)%length_z /= dflt_real)) &
            .or. &
            (patch_ib(patch_id)%length_y > 0d0 .and. &
             (patch_ib(patch_id)%length_x /= dflt_real .or. &
              patch_ib(patch_id)%length_z /= dflt_real)) &
            .or. &
            (patch_ib(patch_id)%length_z > 0d0 .and. &
             (patch_ib(patch_id)%length_x /= dflt_real .or. &
              patch_ib(patch_id)%length_y /= dflt_real)) &
            .or. &
            patch_ib(patch_id)%radius <= 0d0) then

            call s_mpi_abort('Inconsistency(ies) detected in '// &
                             'geometric parameters of cylinder '// &
                             'patch '//trim(iStr)//'. Exiting ...')

        end if

    end subroutine s_check_cylinder_ib_patch_geometry ! -----------------------

    !!>  This subroutine verifies that the geometric parameters of
        !!      the inactive patch remain unaltered by the user inputs.
        !!  @param patch_id Patch identifier
    subroutine s_check_inactive_ib_patch_geometry(patch_id) ! -----------------

        integer, intent(IN) :: patch_id
        call s_int_to_str(patch_id, iStr)

        ! Constraints on the geometric parameters of the inactive patch
        if (patch_ib(patch_id)%x_centroid /= dflt_real &
            .or. &
            patch_ib(patch_id)%y_centroid /= dflt_real &
            .or. &
            patch_ib(patch_id)%z_centroid /= dflt_real &
            .or. &
            patch_ib(patch_id)%length_x /= dflt_real &
            .or. &
            patch_ib(patch_id)%length_y /= dflt_real &
            .or. &
            patch_ib(patch_id)%length_z /= dflt_real &
            .or. &
            patch_ib(patch_id)%radius /= dflt_real) then

            call s_mpi_abort('Inconsistency(ies) detected in '// &
                             'geometric parameters of inactive '// &
                             'patch '//trim(iStr)//'. Exiting ...')

        end if

    end subroutine s_check_inactive_ib_patch_geometry ! -----------------------

end module m_check_ib_patches
