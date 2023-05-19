# 1 "/Users/anand/MFC/src/pre_process/m_initial_condition.fpp"
!>
!! @file m_initial_condition.f90
!! @brief Contains module m_initial_condition

!> @brief This module provides a platform that is analagous to constructive
!!              solid geometry techniques and in this way allows for the creation
!!              of a wide variety of initial conditions. Several 1D, 2D and 3D
!!              fundamental geometries are included that may further be combined
!!              into more complex shapes. This is achieved by carefully setting
!!              up the order in which the patches are laid out in the domain and
!!              specifying the priority that each patch has over the preceeding
!!              ones. The resulting shapes may be identified both by the values
!!              of their primitive variables and the associated patch identities.
!!              Note that the user may choose to read in and modify a preexisting
!!              initial condition. The module m_start_up.f90 is responsible for
!!             reading in the relevant data files.
module m_initial_condition

    ! Dependencies =============================================================
    use m_derived_types         ! Definitions of the derived types

    use m_global_parameters     ! Global parameters for the code

    use m_variables_conversion  ! Subroutines to change the state variables from
    ! one form to another

    use m_patches

    use m_assign_variables
    ! ==========================================================================
    ! ==========================================================================

    implicit none

    ! NOTE: The abstract interface allows for the declaration of a pointer to
    ! a procedure such that the choice of the model equations does not have to
    ! be queried every time the patch primitive variables are to be assigned in
    ! a cell in the computational domain.
    type(scalar_field), allocatable, dimension(:) :: q_prim_vf !< primitive variables

    type(scalar_field), allocatable, dimension(:) :: q_cons_vf !< conservative variables

    integer, allocatable, dimension(:, :, :) :: patch_id_fp !<
    !! Bookkepping variable used to track the patch identities (id) associated
    !! with each of the cells in the computational domain. Note that only one
    !! patch identity may be associated with any one cell.

contains

    !> Computation of parameters, allocation procedures, and/or
        !!              any other tasks needed to properly setup the module
    subroutine s_initialize_initial_condition_module() ! -------------------

        integer :: i !< generic loop iterator

        ! Allocating the primitive and conservative variables
        allocate (q_prim_vf(1:sys_size))
        allocate (q_cons_vf(1:sys_size))

        do i = 1, sys_size
            allocate (q_prim_vf(i)%sf(0:m, 0:n, 0:p))
            allocate (q_cons_vf(i)%sf(0:m, 0:n, 0:p))
        end do

        ! Allocating the patch identities bookkeeping variable
        allocate (patch_id_fp(0:m, 0:n, 0:p))

        ! Setting default values for conservative and primitive variables so
        ! that in the case that the initial condition is wrongly laid out on
        ! the grid the simulation component will catch the problem on start-
        ! up. The conservative variables do not need to be similarly treated
        ! since they are computed directly from the primitive variables.
        do i = 1, sys_size
            q_cons_vf(i)%sf = dflt_real
            q_prim_vf(i)%sf = dflt_real
        end do

        ! Setting default values for patch identities bookkeeping variable.
        ! This is necessary to avoid any confusion in the assessment of the
        ! extent of application that the overwrite permissions give a patch
        ! when it is being applied in the domain.
        patch_id_fp = 0

    end subroutine s_initialize_initial_condition_module ! -----------------

    !>  This subroutine peruses the patches and depending on the
        !!              type of geometry associated with a particular patch, it
        !!              calls the related subroutine to setup the said geometry
        !!              on the grid using the primitive variables included with
        !!              the patch parameters. The subroutine is complete once the
        !!              primitive variables are converted to conservative ones.
    subroutine s_generate_initial_condition() ! ----------------------------

        integer :: i  !< Generic loop operator

        ! Converting the conservative variables to the primitive ones given
        ! preexisting initial condition data files were read in on start-up
        if (old_ic) then
            call s_convert_conservative_to_primitive_variables(q_cons_vf, &
                                                               q_prim_vf)
        end if

        !  3D Patch Geometries =============================================
        if (p > 0) then

            do i = 1, num_patches

                ! Spherical patch
                if (patch_icpp(i)%geometry == 8) then
                    call s_sphere(i, patch_id_fp, q_prim_vf)

                    ! Cuboidal patch
                elseif (patch_icpp(i)%geometry == 9) then
                    call s_cuboid(i, patch_id_fp, q_prim_vf)

                    ! Cylindrical patch
                elseif (patch_icpp(i)%geometry == 10) then
                    call s_cylinder(i, patch_id_fp, q_prim_vf)

                    ! Swept plane patch
                elseif (patch_icpp(i)%geometry == 11) then
                    call s_sweep_plane(i, patch_id_fp, q_prim_vf)

                    ! Ellipsoidal patch
                elseif (patch_icpp(i)%geometry == 12) then
                    call s_ellipsoid(i, patch_id_fp, q_prim_vf)

                    ! Analytical function patch for testing purposes
                elseif (patch_icpp(i)%geometry == 13) then
                    call s_3D_analytical(i, patch_id_fp, q_prim_vf)

                    ! Spherical harmonic patch
                elseif (patch_icpp(i)%geometry == 14) then
                    call s_spherical_harmonic(i, patch_id_fp, q_prim_vf)

                    ! 3D Modified circular patch
                elseif (patch_icpp(i)%geometry == 19) then
                    call s_3dvarcircle(i, patch_id_fp, q_prim_vf)

                end if

            end do

            ! ==================================================================

            ! 2D Patch Geometries ==============================================
        elseif (n > 0) then

            do i = 1, num_patches

                ! Circular patch
                if (patch_icpp(i)%geometry == 2) then
                    call s_circle(i, patch_id_fp, q_prim_vf)

                    ! Rectangular patch
                elseif (patch_icpp(i)%geometry == 3) then
                    call s_rectangle(i, patch_id_fp, q_prim_vf)

                    ! Swept line patch
                elseif (patch_icpp(i)%geometry == 4) then
                    call s_sweep_line(i, patch_id_fp, q_prim_vf)

                    ! Elliptical patch
                elseif (patch_icpp(i)%geometry == 5) then
                    call s_ellipse(i, patch_id_fp, q_prim_vf)

                    ! Isentropic vortex patch
                elseif (patch_icpp(i)%geometry == 6) then
                    call s_isentropic_vortex(i, patch_id_fp, q_prim_vf)

                    ! Analytical function patch for testing purposes
                elseif (patch_icpp(i)%geometry == 7) then
                    call s_2D_analytical(i, patch_id_fp, q_prim_vf)

                    ! Spiral patch
                elseif (patch_icpp(i)%geometry == 17) then
                    call s_spiral(i, patch_id_fp, q_prim_vf)

                    ! Modified circular patch
                elseif (patch_icpp(i)%geometry == 18) then
                    call s_varcircle(i, patch_id_fp, q_prim_vf)

                end if

            end do

            ! ==================================================================

            ! 1D Patch Geometries ==============================================
        else

            do i = 1, num_patches

                ! Line segment patch
                if (patch_icpp(i)%geometry == 1) then
                    call s_line_segment(i, patch_id_fp, q_prim_vf)

                    ! 1d analytical
                elseif (patch_icpp(i)%geometry == 15) then
                    call s_1d_analytical(i, patch_id_fp, q_prim_vf)

                    ! 1d bubble screen with sinusoidal pressure pulse
                elseif (patch_icpp(i)%geometry == 16) then
                    call s_1d_bubble_pulse(i, patch_id_fp, q_prim_vf)
                end if

            end do

        end if
        ! ==================================================================

        if (perturb_flow) call s_perturb_surrounding_flow()
        if (perturb_sph) call s_perturb_sphere()

        ! Converting the primitive variables to the conservative ones
        call s_convert_primitive_to_conservative_variables(q_prim_vf, &
                                                           q_cons_vf)

    end subroutine s_generate_initial_condition ! --------------------------

    subroutine s_perturb_sphere() ! ----------------------------------------

        integer :: i, j, k, l !< generic loop operators

        real(kind(0d0)) :: perturb_alpha
        real(kind(0d0)) :: alpha_unadv
        real(kind(0d0)) :: rand_real
        call random_seed()

        do k = 0, p
            do j = 0, n
                do i = 0, m
                    call random_number(rand_real)

                    perturb_alpha = q_prim_vf(E_idx + perturb_sph_fluid)%sf(i, j, k)

                    ! Perturb partial density fields to match perturbed volume fraction fields
!                        IF ((perturb_alpha >= 25d-2) .AND. (perturb_alpha <= 75d-2)) THEN
                    if ((perturb_alpha /= 0d0) .and. (perturb_alpha /= 1d0)) then

                        ! Derive new partial densities
                        do l = 1, num_fluids
                            q_prim_vf(l)%sf(i, j, k) = q_prim_vf(E_idx + l)%sf(i, j, k)*fluid_rho(l)
                        end do

                    end if
                end do
            end do
        end do

    end subroutine s_perturb_sphere ! --------------------------------------

    subroutine s_perturb_surrounding_flow() ! ------------------------------

        integer :: i, j, k, l !<  generic loop iterators

        real(kind(0d0)) :: perturb_alpha
        real(kind(0d0)) :: rand_real
        call random_seed()

        ! Perturb partial density or velocity of surrounding flow by some random small amount of noise
        do k = 0, p
            do j = 0, n
                do i = 0, m

                    perturb_alpha = q_prim_vf(E_idx + perturb_flow_fluid)%sf(i, j, k)
                    ! IF (perturb_alpha == 1d0) THEN
                    ! Perturb partial density
!                            CALL RANDOM_NUMBER(rand_real)
!                            rand_real = rand_real / 1d2 / 1d3
!                            q_prim_vf(perturb_flow_fluid)%sf(i,j,k) = q_prim_vf(perturb_flow_fluid)%sf(i,j,k) + rand_real
                    ! Perturb velocity
                    call random_number(rand_real)
                    rand_real = rand_real*1.d-2
                    q_prim_vf(mom_idx%beg)%sf(i, j, k) = (1.d0 + rand_real)*q_prim_vf(mom_idx%beg)%sf(i, j, k)
                    q_prim_vf(mom_idx%end)%sf(i, j, k) = rand_real*q_prim_vf(mom_idx%beg)%sf(i, j, k)
                    if (bubbles) then
                        q_prim_vf(alf_idx)%sf(i, j, k) = (1.d0 + rand_real)*q_prim_vf(alf_idx)%sf(i, j, k)
                    end if
                    ! END IF
                end do
            end do
        end do

    end subroutine s_perturb_surrounding_flow ! ----------------------------

    !>  Deallocation procedures for the module
    subroutine s_finalize_initial_condition_module() ! ---------------------

        integer :: i !< Generic loop iterator

        ! Dellocating the primitive and conservative variables
        do i = 1, sys_size
            deallocate (q_prim_vf(i)%sf)
            deallocate (q_cons_vf(i)%sf)
        end do

        deallocate (q_prim_vf)
        deallocate (q_cons_vf)

        ! Deallocating the patch identities bookkeeping variable
        deallocate (patch_id_fp)

    end subroutine s_finalize_initial_condition_module ! -------------------

end module m_initial_condition
