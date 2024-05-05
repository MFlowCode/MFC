!> @brief This module contains subroutines that read, and check consistency
!!              of, the user provided inputs and data.
module m_check_patches

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

    private; public :: s_check_patches

    character(len=10) :: iStr

contains

    subroutine s_check_patches()

        ! integer, intent(in) :: i

        integer :: i

        do i = 1, num_patches_max
            if (i <= num_patches) then
                ! call s_check_patch_geometry(i)
                call s_int_to_str(i, iStr)
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
                elseif (patch_icpp(i)%geometry == 6) then
                    call s_mpi_abort('Unimplemented choice of geometry 6'// &
                                     ' (formerly "Vortex") of active patch '//trim(iStr)// &
                                     ' detected. Exiting ...')
                elseif (patch_icpp(i)%geometry == 7) then
                    call s_check_2D_analytical_patch_geometry(i)
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
                elseif (patch_icpp(i)%geometry == 13) then
                    call s_check_3D_analytical_patch_geometry(i)
                elseif (patch_icpp(i)%geometry == 14) then
                    call s_check_spherical_harmonic_patch_geometry(i)
                elseif (patch_icpp(i)%geometry == 15) then
                    call s_check_1d_analytical_patch_geometry(i)
                elseif (patch_icpp(i)%geometry == 16) then
                    print *, '1d pressure sinusoid'
                elseif (patch_icpp(i)%geometry == 17) then
                    print *, '2d spiral'
                elseif (patch_icpp(i)%geometry == 18) then
                    print *, '2d var circle'
                elseif (patch_icpp(i)%geometry == 19) then
                    print *, '3d var circle'
                elseif (patch_icpp(i)%geometry == 20) then
                    call s_check_2D_TaylorGreen_vortex_patch_geometry(i)
                elseif (patch_icpp(i)%geometry == 21) then
                    call s_check_model_geometry(i)
                else
                    call s_mpi_abort('Unsupported choice of the '// &
                                     'geometry of active patch '//trim(iStr)// &
                                     ' detected. Exiting ...')
                end if
            else
                if (patch_icpp(i)%geometry == dflt_int) then
                    call s_check_inactive_patch_geometry(i)
                else
                    call s_mpi_abort('Unsupported choice of the '// &
                                     'geometry of inactive patch '//trim(iStr)// &
                                     ' detected. Exiting ...')
                end if
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
                             patch_icpp(i)%geometry == 12)) then
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

    !> This subroutine verifies that the geometric parameters of
        !!      the line segment patch have consistently been inputted by
        !!      the user.
        !!  @param patch_id Patch identifier
    subroutine s_check_line_segment_patch_geometry(patch_id) ! -------------

        integer, intent(IN) :: patch_id
        call s_int_to_str(patch_id, iStr)

        ! Constraints on the geometric parameters of the line segment patch
        if (n > 0 .or. patch_icpp(patch_id)%length_x <= 0d0 &
            .or. &
            patch_icpp(patch_id)%x_centroid == dflt_real &
            .or. &
            cyl_coord) then

            call s_mpi_abort('Inconsistency(ies) detected in '// &
                             'geometric parameters of line segment '// &
                             'patch '//trim(iStr)//'. Exiting ...')

        end if

    end subroutine s_check_line_segment_patch_geometry ! -------------------

    !>  This subroutine verifies that the geometric parameters of
        !!      the circle patch have consistently been inputted by the
        !!      user.
        !!  @param patch_id Patch identifier
    subroutine s_check_circle_patch_geometry(patch_id) ! -------------------

        integer, intent(IN) :: patch_id
        call s_int_to_str(patch_id, iStr)

        ! Constraints on the geometric parameters of the circle patch
        if (n == 0 .or. p > 0 .or. patch_icpp(patch_id)%radius <= 0d0 &
            .or. &
            patch_icpp(patch_id)%x_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%y_centroid == dflt_real) then

            call s_mpi_abort('Inconsistency(ies) detected in '// &
                             'geometric parameters of circle '// &
                             'patch '//trim(iStr)//'. Exiting ...')

        end if

    end subroutine s_check_circle_patch_geometry ! -------------------------

    !>  This subroutine verifies that the geometric parameters of
        !!      the rectangle patch have consistently been inputted by
        !!      the user.
        !!  @param patch_id Patch identifier
    subroutine s_check_rectangle_patch_geometry(patch_id) ! ----------------

        integer, intent(IN) :: patch_id
        call s_int_to_str(patch_id, iStr)

        ! Constraints on the geometric parameters of the rectangle patch
        if (n == 0 .or. p > 0 &
            .or. &
            patch_icpp(patch_id)%x_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%y_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%length_x <= 0d0 &
            .or. &
            patch_icpp(patch_id)%length_y <= 0d0) then

            call s_mpi_abort('Inconsistency(ies) detected in '// &
                             'geometric parameters of rectangle '// &
                             'patch '//trim(iStr)//'. Exiting ...')

        end if

    end subroutine s_check_rectangle_patch_geometry ! ----------------------

    !> This subroutine verifies that the geometric parameters of
        !!      the line sweep patch have consistently been inputted by
        !!      the user.
        !!  @param patch_id Patch identifier
    subroutine s_check_line_sweep_patch_geometry(patch_id) ! ---------------

        integer, intent(IN) :: patch_id
        call s_int_to_str(patch_id, iStr)

        ! Constraints on the geometric parameters of the line sweep patch
        if (n == 0 .or. p > 0 &
            .or. &
            patch_icpp(patch_id)%x_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%y_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%normal(1) == dflt_real &
            .or. &
            patch_icpp(patch_id)%normal(2) == dflt_real &
            .or. &
            patch_icpp(patch_id)%normal(3) /= dflt_real) then

            call s_mpi_abort('Inconsistency(ies) detected in '// &
                             'geometric parameters of line sweep '// &
                             'patch '//trim(iStr)//'. Exiting ...')

        end if

    end subroutine s_check_line_sweep_patch_geometry ! ---------------------

    !>  This subroutine verifies that the geometric parameters of
        !!      the ellipse patch have consistently been inputted by
        !!      the user.
        !!  @param patch_id Patch identifier
    subroutine s_check_ellipse_patch_geometry(patch_id) ! ------------------

        integer, intent(IN) :: patch_id
        call s_int_to_str(patch_id, iStr)

        ! Constraints on the geometric parameters of the ellipse patch
        if (n == 0 .or. p > 0 &
            .or. &
            patch_icpp(patch_id)%x_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%y_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%radii(1) == dflt_real &
            .or. &
            patch_icpp(patch_id)%radii(2) == dflt_real &
            .or. &
            patch_icpp(patch_id)%radii(3) /= dflt_real) then

            call s_mpi_abort('Inconsistency(ies) detected in '// &
                             'geometric parameters of ellipse '// &
                             'patch '//trim(iStr)//'. Exiting ...')

        end if

    end subroutine s_check_ellipse_patch_geometry ! ------------------------

    !>  This subroutine verifies that the geometric parameters of
        !!      the Taylor Green vortex patch have been entered by the user
        !!      consistently.
        !!  @param patch_id Patch identifier
    subroutine s_check_2D_TaylorGreen_vortex_patch_geometry(patch_id) ! --------

        integer, intent(IN) :: patch_id
        call s_int_to_str(patch_id, iStr)

        ! Constraints on the TaylorGreen vortex patch geometric parameters
        if (n == 0 .or. p > 0 &
            .or. &
            patch_icpp(patch_id)%x_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%y_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%length_x <= 0d0 &
            .or. &
            patch_icpp(patch_id)%length_y <= 0d0 &
            .or. &
            patch_icpp(patch_id)%vel(2) <= 0d0) then

            call s_mpi_abort('Inconsistency(ies) detected in '// &
                             'geometric parameters of Taylor Green '// &
                             'vortex patch '//trim(iStr)//'. Exiting ...')

        end if

    end subroutine s_check_2D_TaylorGreen_vortex_patch_geometry! --------------

    !>  This subroutine verifies that the geometric parameters of
        !!      the analytical patch have consistently been inputted by
        !!      the user.
        !!  @param patch_id Patch identifier
    subroutine s_check_1D_analytical_patch_geometry(patch_id) ! ---------------

        integer, intent(IN) :: patch_id
        call s_int_to_str(patch_id, iStr)

        ! Constraints on the geometric parameters of the analytical patch
        if (n > 0 .or. p > 0 &
            .or. &
            (model_eqns /= 4 .and. model_eqns /= 2) &
            .or. &
            patch_icpp(patch_id)%x_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%length_x <= 0d0) then

            call s_mpi_abort('Inconsistency(ies) detected in '// &
                             'geometric parameters of 1D analytical '// &
                             'patch '//trim(iStr)//'. Exiting...')
        end if
    end subroutine s_check_1D_analytical_patch_geometry ! ---------------------

    !>  This subroutine verifies that the geometric parameters of
        !!      the analytical patch have consistently been inputted by
        !!      the user.
        !!  @param patch_id Patch identifier
    subroutine s_check_2D_analytical_patch_geometry(patch_id) ! ---------------

        integer, intent(IN) :: patch_id
        call s_int_to_str(patch_id, iStr)

        ! Constraints on the geometric parameters of the analytical patch
        if (n == 0 .or. p > 0 &
            .or. &
            patch_icpp(patch_id)%x_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%y_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%length_x <= 0d0 &
            .or. &
            patch_icpp(patch_id)%length_y <= 0d0) then

            call s_mpi_abort('Inconsistency(ies) detected in '// &
                             'geometric parameters of 2D analytical '// &
                             'patch '//trim(iStr)//'. Exiting...')
        end if
    end subroutine s_check_2D_analytical_patch_geometry ! ---------------------

    !>  This subroutine verifies that the geometric parameters of
        !!      the analytical patch have consistently been inputted by
        !!      the user.
        !!  @param patch_id Patch identifier
    subroutine s_check_3D_analytical_patch_geometry(patch_id) ! ---------------

        integer, intent(IN) :: patch_id
        call s_int_to_str(patch_id, iStr)

        ! Constraints on the geometric parameters of the analytical patch
        if (p == 0 &
            .or. &
            patch_icpp(patch_id)%x_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%y_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%z_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%length_x <= 0d0 &
            .or. &
            patch_icpp(patch_id)%length_y <= 0d0 &
            .or. &
            patch_icpp(patch_id)%length_z <= 0d0) then

            call s_mpi_abort('Inconsistency(ies) detected in '// &
                             'geometric parameters of 3D analytical '// &
                             'patch '//trim(iStr)//'. Exiting...')
        end if
    end subroutine s_check_3D_analytical_patch_geometry ! ---------------------

    !> This subroutine verifies that the geometric parameters of
        !!      the sphere patch have consistently been inputted by the
        !!      user.
        !!  @param patch_id Patch identifier
    subroutine s_check_sphere_patch_geometry(patch_id) ! -------------------

        integer, intent(IN) :: patch_id
        call s_int_to_str(patch_id, iStr)

        ! Constraints on the geometric parameters of the sphere patch
        if (p == 0 &
            .or. &
            patch_icpp(patch_id)%radius <= 0d0 &
            .or. &
            patch_icpp(patch_id)%x_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%y_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%z_centroid == dflt_real) then

            call s_mpi_abort('Inconsistency(ies) detected in '// &
                             'geometric parameters of sphere '// &
                             'patch '//trim(iStr)//'. Exiting ...')

        end if

    end subroutine s_check_sphere_patch_geometry ! -------------------------

    !>  This subroutine verifies that the geometric parameters of
        !!      the spherical harmonic  patch have consistently been
        !!      inputted by the user.
        !!  @param patch_id Patch identifier
    subroutine s_check_spherical_harmonic_patch_geometry(patch_id) ! -------

        integer, intent(IN) :: patch_id
        call s_int_to_str(patch_id, iStr)

        ! Constraints on the geometric parameters of the spherical harmonic patch
        if (p == 0 &
            .or. &
            patch_icpp(patch_id)%radius <= 0d0 &
            .or. &
            patch_icpp(patch_id)%x_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%y_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%z_centroid == dflt_real &
            .or. &
            all(patch_icpp(patch_id)%epsilon /= (/1d0, 2d0, 3d0, 4d0, 5d0/)) &
            .or. &
            patch_icpp(patch_id)%beta < 0d0 &
            .or. &
            patch_icpp(patch_id)%beta > patch_icpp(patch_id)%epsilon) then

            call s_mpi_abort('Inconsistency(ies) detected in '// &
                             'geometric parameters of spherical '// &
                             'harmonic patch '//trim(iStr)//'. Exiting ...')

        end if

    end subroutine s_check_spherical_harmonic_patch_geometry ! -------------

    !>  This subroutine verifies that the geometric parameters of
        !!      the cuboid patch have consistently been inputted by the
        !!      user.
        !!  @param patch_id Patch identifier
    subroutine s_check_cuboid_patch_geometry(patch_id) ! -------------------

        ! Patch identifier
        integer, intent(IN) :: patch_id
        call s_int_to_str(patch_id, iStr)

        ! Constraints on the geometric parameters of the cuboid patch
        if (p == 0 &
            .or. &
            patch_icpp(patch_id)%x_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%y_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%z_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%length_x <= 0d0 &
            .or. &
            patch_icpp(patch_id)%length_y <= 0d0 &
            .or. &
            patch_icpp(patch_id)%length_z <= 0d0) then

            call s_mpi_abort('Inconsistency(ies) detected in '// &
                             'geometric parameters of cuboid '// &
                             'patch '//trim(iStr)//'. Exiting ...')

        end if

    end subroutine s_check_cuboid_patch_geometry ! -------------------------

    !>  This subroutine verifies that the geometric parameters of
        !!      the cylinder patch have consistently been inputted by the
        !!      user.
        !!  @param patch_id Patch identifier
    subroutine s_check_cylinder_patch_geometry(patch_id) ! -----------------

        ! Patch identifier
        integer, intent(IN) :: patch_id
        call s_int_to_str(patch_id, iStr)

        ! Constraints on the geometric parameters of the cylinder patch
        if (p == 0 &
            .or. &
            patch_icpp(patch_id)%x_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%y_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%z_centroid == dflt_real &
            .or. &
            (patch_icpp(patch_id)%length_x <= 0d0 .and. &
             patch_icpp(patch_id)%length_y <= 0d0 .and. &
             patch_icpp(patch_id)%length_z <= 0d0) &
            .or. &
            (patch_icpp(patch_id)%length_x > 0d0 .and. &
             (patch_icpp(patch_id)%length_y /= dflt_real .or. &
              patch_icpp(patch_id)%length_z /= dflt_real)) &
            .or. &
            (patch_icpp(patch_id)%length_y > 0d0 .and. &
             (patch_icpp(patch_id)%length_x /= dflt_real .or. &
              patch_icpp(patch_id)%length_z /= dflt_real)) &
            .or. &
            (patch_icpp(patch_id)%length_z > 0d0 .and. &
             (patch_icpp(patch_id)%length_x /= dflt_real .or. &
              patch_icpp(patch_id)%length_y /= dflt_real)) &
            .or. &
            patch_icpp(patch_id)%radius <= 0d0) then

            call s_mpi_abort('Inconsistency(ies) detected in '// &
                             'geometric parameters of cylinder '// &
                             'patch '//trim(iStr)//'. Exiting ...')

        end if

    end subroutine s_check_cylinder_patch_geometry ! -----------------------

    !>  This subroutine verifies that the geometric parameters of
        !!      the plane sweep patch have consistently been inputted by
        !!      the user.
        !!  @param patch_id Patch identifier
    subroutine s_check_plane_sweep_patch_geometry(patch_id) ! --------------

        ! Patch identifier
        integer, intent(IN) :: patch_id
        call s_int_to_str(patch_id, iStr)

        ! Constraints on the geometric parameters of the plane sweep patch
        if (p == 0 &
            .or. &
            patch_icpp(patch_id)%x_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%y_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%z_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%normal(1) == dflt_real &
            .or. &
            patch_icpp(patch_id)%normal(2) == dflt_real &
            .or. &
            patch_icpp(patch_id)%normal(3) == dflt_real) then

            call s_mpi_abort('Inconsistency(ies) detected in '// &
                             'geometric parameters of plane sweep '// &
                             'patch '//trim(iStr)//'. Exiting ...')

        end if

    end subroutine s_check_plane_sweep_patch_geometry ! --------------------

    !> This subroutine verifies that the geometric parameters of
        !!      the ellipsoid patch have consistently been inputted by
        !!      the user.
        !!  @param patch_id Patch identifier
    subroutine s_check_ellipsoid_patch_geometry(patch_id) ! ----------------

        integer, intent(IN) :: patch_id
        call s_int_to_str(patch_id, iStr)

        ! Constraints on the geometric parameters of the ellipsoid patch
        if (p == 0 &
            .or. &
            patch_icpp(patch_id)%x_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%y_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%z_centroid == dflt_real &
            .or. &
            patch_icpp(patch_id)%radii(1) == dflt_real &
            .or. &
            patch_icpp(patch_id)%radii(2) == dflt_real &
            .or. &
            patch_icpp(patch_id)%radii(3) == dflt_real) then

            call s_mpi_abort('Inconsistency(ies) detected in '// &
                             'geometric parameters of ellipsoid '// &
                             'patch '//trim(iStr)//'. Exiting ...')

        end if

    end subroutine s_check_ellipsoid_patch_geometry ! ----------------------

    !!>  This subroutine verifies that the geometric parameters of
        !!      the inactive patch remain unaltered by the user inputs.
        !!  @param patch_id Patch identifier
    subroutine s_check_inactive_patch_geometry(patch_id) ! -----------------

        integer, intent(IN) :: patch_id
        call s_int_to_str(patch_id, iStr)

        ! Constraints on the geometric parameters of the inactive patch
        if (patch_icpp(patch_id)%x_centroid /= dflt_real &
            .or. &
            patch_icpp(patch_id)%y_centroid /= dflt_real &
            .or. &
            patch_icpp(patch_id)%z_centroid /= dflt_real &
            .or. &
            patch_icpp(patch_id)%length_x /= dflt_real &
            .or. &
            patch_icpp(patch_id)%length_y /= dflt_real &
            .or. &
            patch_icpp(patch_id)%length_z /= dflt_real &
            .or. &
            patch_icpp(patch_id)%radius /= dflt_real &
            .or. &
            patch_icpp(patch_id)%epsilon /= dflt_real &
            .or. &
            patch_icpp(patch_id)%beta /= dflt_real &
            .or. &
            patch_icpp(patch_id)%normal(1) /= dflt_real &
            .or. &
            patch_icpp(patch_id)%normal(2) /= dflt_real &
            .or. &
            patch_icpp(patch_id)%normal(3) /= dflt_real &
            .or. &
            patch_icpp(patch_id)%radii(1) /= dflt_real &
            .or. &
            patch_icpp(patch_id)%radii(2) /= dflt_real &
            .or. &
            patch_icpp(patch_id)%radii(3) /= dflt_real) then

            call s_mpi_abort('Inconsistency(ies) detected in '// &
                             'geometric parameters of inactive '// &
                             'patch '//trim(iStr)//'. Exiting ...')

        end if

    end subroutine s_check_inactive_patch_geometry ! -----------------------

    !>  This subroutine verifies that any rights granted to the
        !!      given active patch, to overwrite the preceding active
        !!      patches, were consistently inputted by the user.
        !!  @param patch_id Patch identifier
    subroutine s_check_active_patch_alteration_rights(patch_id) ! ----------

        integer, intent(IN) :: patch_id
        call s_int_to_str(patch_id, iStr)

        ! Constraints on the alteration rights of an active patch
        if (patch_icpp(patch_id)%alter_patch(0) .eqv. .false. &
            .or. &
            any(patch_icpp(patch_id)%alter_patch(patch_id:))) then

            call s_mpi_abort('Inconsistency(ies) detected in '// &
                             'alteration rights of active '// &
                             'patch '//trim(iStr)//'. Exiting ...')

        end if

    end subroutine s_check_active_patch_alteration_rights ! ----------------

    !>  This subroutine verifies that the rights of the given
        !!      inactive patch, to overwrite the preceding patches,
        !!      remain unaltered by the user inputs.
        !!  @param patch_id Patch identifier
    subroutine s_check_inactive_patch_alteration_rights(patch_id) ! --------

        ! Patch identifier
        integer, intent(IN) :: patch_id
        call s_int_to_str(patch_id, iStr)

        ! Constraints on the alteration rights of an inactive patch
        if (patch_icpp(patch_id)%alter_patch(0) .eqv. .false. &
            .or. &
            any(patch_icpp(patch_id)%alter_patch(1:))) then

            call s_mpi_abort('Inconsistency(ies) detected in '// &
                             'alteration rights of inactive '// &
                             'patch '//trim(iStr)//'. Exiting ...')

        end if

    end subroutine s_check_inactive_patch_alteration_rights ! --------------

    !> This subroutine verifies that the smoothing parameters of
        !!      the given patch, which supports the smoothing out of its
        !!      boundaries, have consistently been inputted by the user.
        !!  @param patch_id Patch identifier
    subroutine s_check_supported_patch_smoothing(patch_id) ! ---------------

        integer, intent(IN) :: patch_id
        call s_int_to_str(patch_id, iStr)

        ! Constraints on the smoothing parameters of a supported patch
        if ((patch_icpp(patch_id)%smoothen &
             .and. &
             (patch_icpp(patch_id)%smooth_patch_id >= patch_id &
              .or. &
              patch_icpp(patch_id)%smooth_patch_id == 0 &
              .or. &
              patch_icpp(patch_id)%smooth_coeff <= 0d0)) &
            .or. &
            ((patch_icpp(patch_id)%smoothen .neqv. .true.) &
             .and. &
             (patch_icpp(patch_id)%smooth_patch_id /= patch_id &
              .or. &
              patch_icpp(patch_id)%smooth_coeff /= dflt_real))) then

            call s_mpi_abort('Inconsistency(ies) detected in '// &
                             'smoothing parameters of supported '// &
                             'patch '//trim(iStr)//'. Exiting ...')

        end if

    end subroutine s_check_supported_patch_smoothing ! ---------------------

    !> This subroutine verifies that the smoothing parameters of
        !!      the given patch, which does not support the smoothing out
        !!          of its boundaries, remain unaltered by the user inputs.
        !!  @param patch_id Patch identifier
    subroutine s_check_unsupported_patch_smoothing(patch_id) ! -------------

        ! Patch identifier
        integer, intent(IN) :: patch_id
        ! call s_int_to_str(patch_id, iStr)

        ! Constraints on the smoothing parameters of an unsupported patch
        if (patch_icpp(patch_id)%smoothen &
            .or. &
            patch_icpp(patch_id)%smooth_patch_id /= patch_id &
            .or. &
            patch_icpp(patch_id)%smooth_coeff /= dflt_real) then

            call s_mpi_abort('Inconsistency(ies) detected in '// &
                             'smoothing parameters of unsupported '// &
                             'patch '//trim(iStr)//'. Exiting ...')

        end if

    end subroutine s_check_unsupported_patch_smoothing ! -------------------

    !>  This subroutine verifies that the primitive variables
        !!      associated with the given active patch are physically
        !!      consistent.
        !!  @param patch_id Patch identifier
    subroutine s_check_active_patch_primitive_variables(patch_id) ! --------

        integer, intent(IN) :: patch_id
        call s_int_to_str(patch_id, iStr)

        ! Constraints on the primitive variables of an active patch
        if (patch_icpp(patch_id)%vel(1) == dflt_real &
            .or. &
            (n == 0 .and. patch_icpp(patch_id)%vel(2) /= dflt_real .and. patch_icpp(patch_id)%vel(2) /= 0) &
            .or. &
            (n > 0 .and. patch_icpp(patch_id)%vel(2) == dflt_real) &
            .or. &
            (p == 0 .and. patch_icpp(patch_id)%vel(3) /= dflt_real .and. patch_icpp(patch_id)%vel(3) /= 0) &
            .or. &
            (p > 0 .and. patch_icpp(patch_id)%vel(3) == dflt_real) &
            !                               .OR.                           &
            !                 patch_icpp(patch_id)%pres <= 0d0             &
            .or. &
            (model_eqns == 1 .and. &
             (patch_icpp(patch_id)%rho <= 0d0 .or. &
              patch_icpp(patch_id)%gamma <= 0d0 .or. &
              patch_icpp(patch_id)%pi_inf < 0d0)) &
            .or. &
            (patch_icpp(patch_id)%geometry == 5 &
             .and. &
             patch_icpp(patch_id)%pi_inf > 0) &
            .or. &
            (model_eqns == 2 &
             .and. &
             (any(patch_icpp(patch_id)%alpha_rho(1:num_fluids) < 0d0)))) then

            call s_mpi_abort('Inconsistency(ies) detected in '// &
                             'primitive variables of active '// &
                             'patch '//trim(iStr)//'. Exiting ...')

        end if

        if (model_eqns == 2 .and. num_fluids < num_fluids) then

            if (any(patch_icpp(patch_id)%alpha_rho(num_fluids + 1:) &
                    /= dflt_real) &
                .or. &
                any(patch_icpp(patch_id)%alpha(num_fluids + 1:) &
                    /= dflt_real) &
                .or. &
                (patch_icpp(patch_id)%alpha(num_fluids) == dflt_real)) then

                call s_mpi_abort('Inconsistency(ies) detected in '// &
                                 'primitive variables of active '// &
                                 'patch '//trim(iStr)//'. Exiting ...')

            end if

        end if

    end subroutine s_check_active_patch_primitive_variables ! --------------

    !>  This subroutine verifies that the primitive variables
        !!      associated with the given inactive patch remain unaltered
        !!      by the user inputs.
        !!  @param patch_id Patch identifier
    subroutine s_check_inactive_patch_primitive_variables(patch_id) ! ------

        integer, intent(IN) :: patch_id
        call s_int_to_str(patch_id, iStr)

        ! Constraints on the primitive variables of an inactive patch
        if (any(patch_icpp(patch_id)%alpha_rho /= dflt_real) &
            .or. &
            patch_icpp(patch_id)%rho /= dflt_real &
            .or. &
            any(patch_icpp(patch_id)%vel /= dflt_real) &
            .or. &
            patch_icpp(patch_id)%pres /= dflt_real &
            .or. &
            any(patch_icpp(patch_id)%alpha /= dflt_real) &
            .or. &
            patch_icpp(patch_id)%gamma /= dflt_real &
            .or. &
            patch_icpp(patch_id)%pi_inf /= dflt_real) then

            call s_mpi_abort('Inconsistency(ies) detected in '// &
                             'primitive variables of inactive '// &
                             'patch '//trim(iStr)//'. Exiting ...')

        end if

    end subroutine s_check_inactive_patch_primitive_variables ! ------------

    subroutine s_check_model_geometry(patch_id) ! ------------------------------

        integer, intent(IN) :: patch_id

        logical :: file_exists

        inquire (file=patch_icpp(patch_id)%model%filepath, exist=file_exists)

        if (.not. file_exists) then

            print '(A,I0,A)', 'Model file '//trim(patch_icpp(patch_id)%model%filepath)// &
                ' requested by patch ', patch_id, ' does not exist. Exiting ...'

            call s_mpi_abort()

        end if

    end subroutine s_check_model_geometry ! -----------------------------------

end module m_check_patches
