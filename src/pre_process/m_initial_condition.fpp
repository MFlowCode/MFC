!>
!! @file m_initial_condition.f90
!! @brief Contains module m_initial_condition

!> @brief This module provides a platform that is analogous to constructive
!!              solid geometry techniques and in this way allows for the creation
!!              of a wide variety of initial conditions. Several 1D, 2D and 3D
!!              fundamental geometries are included that may further be combined
!!              into more complex shapes. This is achieved by carefully setting
!!              up the order in which the patches are laid out in the domain and
!!              specifying the priority that each patch has over the preceding
!!              ones. The resulting shapes may be identified both by the values
!!              of their primitive variables and the associated patch identities.
!!              Note that the user may choose to read in and modify a preexisting
!!              initial condition. The module m_start_up.f90 is responsible for
!!             reading in the relevant data files.
module m_initial_condition

    use m_derived_types         ! Definitions of the derived types

    use m_global_parameters     ! Global parameters for the code

    use m_mpi_proxy              !< Message passing interface (MPI) module proxy

    use m_helper

    use m_variables_conversion  ! Subroutines to change the state variables from
    ! one form to another

    use m_patches

    use m_compute_levelset      ! Subroutines to calculate levelsets for IBs

    use m_assign_variables

    use m_perturbation          ! Subroutines to perturb initial flow fields

    use m_chemistry

    implicit none

    ! NOTE: The abstract interface allows for the declaration of a pointer to
    ! a procedure such that the choice of the model equations does not have to
    ! be queried every time the patch primitive variables are to be assigned in
    ! a cell in the computational domain.
    type(scalar_field), allocatable, dimension(:) :: q_prim_vf !< primitive variables

    type(scalar_field), allocatable, dimension(:) :: q_cons_vf !< conservative variables

    type(scalar_field) :: q_T_sf !< Temperature field

    integer, allocatable, dimension(:, :, :) :: patch_id_fp !<
    !! Bookkepping variable used to track the patch identities (id) associated
    !! with each of the cells in the computational domain. Note that only one
    !! patch identity may be associated with any one cell.

    type(integer_field) :: ib_markers !<
    !! Bookkepping variable used to track whether a given cell is within an
    !! immersed boundary. The default is 0, otherwise the value is assigned
    !! to the patch ID of the immersed boundary.

    type(levelset_field) :: levelset
    type(levelset_norm_field) :: levelset_norm

contains

    !> Computation of parameters, allocation procedures, and/or
        !!              any other tasks needed to properly setup the module
    subroutine s_initialize_initial_condition_module

        integer :: i !< generic loop iterator

        ! Allocating the primitive and conservative variables
        allocate (q_prim_vf(1:sys_size))
        allocate (q_cons_vf(1:sys_size))

        do i = 1, sys_size
            allocate (q_prim_vf(i)%sf(0:m, 0:n, 0:p))
            allocate (q_cons_vf(i)%sf(0:m, 0:n, 0:p))
        end do

        if (chemistry) then
            allocate (q_T_sf%sf(0:m, 0:n, 0:p))
        end if

        ! Allocating the patch identities bookkeeping variable
        allocate (patch_id_fp(0:m, 0:n, 0:p))

        allocate (ib_markers%sf(0:m, 0:n, 0:p))

        allocate (levelset%sf(0:m, 0:n, 0:p, 1:num_ibs))
        allocate (levelset_norm%sf(0:m, 0:n, 0:p, 1:num_ibs, 1:3))

        if (qbmm .and. .not. polytropic) then
            !Allocate bubble pressure pb and vapor mass mv for non-polytropic qbmm at all quad nodes and R0 bins
            allocate (pb%sf(0:m, &
                            0:n, &
                            0:p, 1:nnode, 1:nb))
            allocate (mv%sf(0:m, &
                            0:n, &
                            0:p, 1:nnode, 1:nb))
        end if

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
        ib_markers%sf = 0

    end subroutine s_initialize_initial_condition_module

    !>  This subroutine peruses the patches and depending on the
        !!              type of geometry associated with a particular patch, it
        !!              calls the related subroutine to setup the said geometry
        !!              on the grid using the primitive variables included with
        !!              the patch parameters. The subroutine is complete once the
        !!              primitive variables are converted to conservative ones.
    subroutine s_generate_initial_condition

        integer :: i  !< Generic loop operator

        character(len=10) :: iStr

        ! First, compute the temperature field from the conservative variables.
        if (chemistry) call s_compute_q_T_sf(q_T_sf, q_cons_vf, idwbuff)

        ! Converting the conservative variables to the primitive ones given
        ! preexisting initial condition data files were read in on start-up
        if (old_ic) then
            call s_convert_conservative_to_primitive_variables(q_cons_vf, &
                                                               q_T_sf, &
                                                               q_prim_vf, &
                                                               idwbuff)
        end if

        !  3D Patch Geometries
        if (p > 0) then

            do i = 1, num_patches

                if (proc_rank == 0) then
                    print *, 'Processing patch', i
                end if

                !> ICPP Patches
                !> @{
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

                    ! 3D STL patch
                elseif (patch_icpp(i)%geometry == 21) then
                    call s_model(i, patch_id_fp, q_prim_vf)

                end if

            end do
            !> @}

            !> IB Patches
            !> @{
            ! Spherical patch
            do i = 1, num_ibs
                if (proc_rank == 0) then
                    print *, 'Processing 3D ib patch ', i
                end if

                if (patch_ib(i)%geometry == 8) then
                    call s_sphere(i, ib_markers%sf, q_prim_vf, ib)
                    call s_sphere_levelset(levelset, levelset_norm, i)
                    ! Cylindrical patch
                elseif (patch_ib(i)%geometry == 9) then
                    call s_cuboid(i, ib_markers%sf, q_prim_vf, ib)
                    call s_cuboid_levelset(levelset, levelset_norm, i)
                elseif (patch_ib(i)%geometry == 10) then
                    call s_cylinder(i, ib_markers%sf, q_prim_vf, ib)
                    call s_cylinder_levelset(levelset, levelset_norm, i)
                elseif (patch_ib(i)%geometry == 11) then
                    call s_3D_airfoil(i, ib_markers%sf, q_prim_vf, ib)
                    call s_3D_airfoil_levelset(levelset, levelset_norm, i)

                    ! STL+IBM patch
                elseif (patch_ib(i)%geometry == 12) then
                    call s_model(i, ib_markers%sf, q_prim_vf, ib, levelset, levelset_norm)
                end if
            end do
            !> @}

            ! 2D Patch Geometries
        elseif (n > 0) then

            do i = 1, num_patches

                if (proc_rank == 0) then
                    print *, 'Processing patch', i
                end if

                !> ICPP Patches
                !> @{
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

                    ! Unimplemented patch (formerly isentropic vortex)
                elseif (patch_icpp(i)%geometry == 6) then
                    call s_mpi_abort('This used to be the isentropic vortex patch, '// &
                                     'which no longer exists. See Examples. Exiting.')

                    ! Analytical function patch for testing purposes
                elseif (patch_icpp(i)%geometry == 7) then
                    call s_2D_analytical(i, patch_id_fp, q_prim_vf)

                    ! Spherical Harmonic Patch
                elseif (patch_icpp(i)%geometry == 14) then
                    call s_spherical_harmonic(i, patch_id_fp, q_prim_vf)

                    ! Spiral patch
                elseif (patch_icpp(i)%geometry == 17) then
                    call s_spiral(i, patch_id_fp, q_prim_vf)

                    ! Modified circular patch
                elseif (patch_icpp(i)%geometry == 18) then
                    call s_varcircle(i, patch_id_fp, q_prim_vf)

                    ! TaylorGreen vortex patch
                elseif (patch_icpp(i)%geometry == 20) then
                    call s_2D_TaylorGreen_vortex(i, patch_id_fp, q_prim_vf)

                    ! STL patch
                elseif (patch_icpp(i)%geometry == 21) then
                    call s_model(i, patch_id_fp, q_prim_vf)

                end if
                !> @}
            end do

            !> IB Patches
            !> @{
            do i = 1, num_ibs
                if (proc_rank == 0) then
                    print *, 'Processing 2D ib patch ', i
                end if
                if (patch_ib(i)%geometry == 2) then
                    call s_circle(i, ib_markers%sf, q_prim_vf, ib)
                    call s_circle_levelset(levelset, levelset_norm, i)
                    ! Rectangular patch
                elseif (patch_ib(i)%geometry == 3) then
                    call s_rectangle(i, ib_markers%sf, q_prim_vf, ib)
                    call s_rectangle_levelset(levelset, levelset_norm, i)
                elseif (patch_ib(i)%geometry == 4) then
                    call s_airfoil(i, ib_markers%sf, q_prim_vf, ib)
                    call s_airfoil_levelset(levelset, levelset_norm, i)
                    ! STL+IBM patch
                elseif (patch_ib(i)%geometry == 5) then
                    call s_model(i, ib_markers%sf, q_prim_vf, ib, levelset, levelset_norm)
                end if
            end do
            !> @}

            ! 1D Patch Geometries
        else

            do i = 1, num_patches

                if (proc_rank == 0) then
                    print *, 'Processing patch', i
                end if

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

        if (perturb_flow) call s_perturb_surrounding_flow(q_prim_vf)
        if (perturb_sph) call s_perturb_sphere(q_prim_vf)
        if (mixlayer_perturb) call s_superposition_instability_wave(q_prim_vf)

        ! Converting the primitive variables to the conservative ones
        call s_convert_primitive_to_conservative_variables(q_prim_vf, q_cons_vf)

        if (chemistry) call s_compute_q_T_sf(q_T_sf, q_cons_vf, idwint)

        if (qbmm .and. .not. polytropic) then
            !Initialize pb and mv
            call s_initialize_mv(q_cons_vf, mv%sf)
            call s_initialize_pb(q_cons_vf, mv%sf, pb%sf)
        end if

    end subroutine s_generate_initial_condition

    !>  Deallocation procedures for the module
    subroutine s_finalize_initial_condition_module

        integer :: i !< Generic loop iterator

        ! Dellocating the primitive and conservative variables
        do i = 1, sys_size
            deallocate (q_prim_vf(i)%sf)
            deallocate (q_cons_vf(i)%sf)
        end do

        deallocate (q_prim_vf)
        deallocate (q_cons_vf)

        if (chemistry) then
            deallocate (q_T_sf%sf)
        end if

        ! Deallocating the patch identities bookkeeping variable
        deallocate (patch_id_fp)
        deallocate (ib_markers%sf)

    end subroutine s_finalize_initial_condition_module

end module m_initial_condition
