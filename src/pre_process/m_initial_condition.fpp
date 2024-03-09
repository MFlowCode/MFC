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

    ! Dependencies =============================================================
    use m_derived_types         ! Definitions of the derived types

    use m_global_parameters     ! Global parameters for the code

    use m_mpi_proxy              !< Message passing interface (MPI) module proxy

    use m_helper

    use m_variables_conversion  ! Subroutines to change the state variables from
    ! one form to another

    use m_patches

    use m_assign_variables

    use m_eigen_solver          ! Subroutines to solve eigenvalue problem for
    ! complex general matrix

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

    type(integer_field) :: ib_markers !<
    !! Bookkepping variable used to track whether a given cell is within an
    !! immersed boundary. The default is 0, otherwise the value is assigned
    !! to the patch ID of the immersed boundary.

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

        allocate (ib_markers%sf(0:m, 0:n, 0:p))

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

    end subroutine s_initialize_initial_condition_module ! -----------------

    !>  This subroutine peruses the patches and depending on the
        !!              type of geometry associated with a particular patch, it
        !!              calls the related subroutine to setup the said geometry
        !!              on the grid using the primitive variables included with
        !!              the patch parameters. The subroutine is complete once the
        !!              primitive variables are converted to conservative ones.
    subroutine s_generate_initial_condition() ! ----------------------------

        integer :: i  !< Generic loop operator

        character(len=10) :: iStr

        ! Converting the conservative variables to the primitive ones given
        ! preexisting initial condition data files were read in on start-up
        if (old_ic) then
            call s_convert_conservative_to_primitive_variables(q_cons_vf, &
                                                               q_prim_vf)
        end if

        !  3D Patch Geometries =============================================
        if (p > 0) then

            do i = 1, num_patches

                if (proc_rank == 0) then
                    print *, 'Processing patch', i
                end if

                !> IB Patches
                !> @{
                ! Spherical patch
                if (patch_ib(i)%geometry == 8) then
                    call s_sphere(i, ib_markers%sf, q_prim_vf, .true.)
                    ! Cylindrical patch
                elseif (patch_ib(i)%geometry == 10) then
                    call s_cylinder(i, ib_markers%sf, q_prim_vf, .true.)

                elseif (patch_ib(i)%geometry == 11) then
                    call s_3D_airfoil(i, ib_markers%sf, q_prim_vf, .true.)
                end if
                !> @}

                !> ICPP Patches
                !> @{
                ! Spherical patch
                if (patch_icpp(i)%geometry == 8) then
                    call s_sphere(i, patch_id_fp, q_prim_vf, .false.)

                    ! Cuboidal patch
                elseif (patch_icpp(i)%geometry == 9) then
                    call s_cuboid(i, patch_id_fp, q_prim_vf)

                    ! Cylindrical patch
                elseif (patch_icpp(i)%geometry == 10) then
                    call s_cylinder(i, patch_id_fp, q_prim_vf, .false.)

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

            ! ==================================================================

            ! 2D Patch Geometries ==============================================
        elseif (n > 0) then

            do i = 1, num_patches

                if (proc_rank == 0) then
                    print *, 'Processing patch', i
                end if

                !> IB Patches
                !> @{
                ! Circular patch
                if (patch_ib(i)%geometry == 2) then
                    call s_circle(i, ib_markers%sf, q_prim_vf, .true.)

                    ! Rectangular patch
                elseif (patch_ib(i)%geometry == 3) then
                    call s_rectangle(i, ib_markers%sf, q_prim_vf, .true.)

                elseif (patch_ib(i)%geometry == 4) then
                    call s_airfoil(i, ib_markers%sf, q_prim_vf, .true.)
                end if
                !> @}

                !> ICPP Patches
                !> @{
                ! Circular patch
                if (patch_icpp(i)%geometry == 2) then
                    call s_circle(i, patch_id_fp, q_prim_vf, .false.)

                    ! Rectangular patch
                elseif (patch_icpp(i)%geometry == 3) then
                    call s_rectangle(i, patch_id_fp, q_prim_vf, .false.)

                    ! Swept line patch
                elseif (patch_icpp(i)%geometry == 4) then
                    call s_sweep_line(i, patch_id_fp, q_prim_vf)

                    ! Elliptical patch
                elseif (patch_icpp(i)%geometry == 5) then
                    call s_ellipse(i, patch_id_fp, q_prim_vf)

                    ! Unimplemented patch (formerly isentropic vortex)
                elseif (patch_icpp(i)%geometry == 6) then
                    call s_mpi_abort('This used to be the isentropic vortex patch, '// &
                                     'which no longer exists. See Examples. Exiting ...')

                    ! Analytical function patch for testing purposes
                elseif (patch_icpp(i)%geometry == 7) then
                    call s_2D_analytical(i, patch_id_fp, q_prim_vf)

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

            ! ==================================================================

            ! 1D Patch Geometries ==============================================
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
        ! ==================================================================

        if (perturb_flow) call s_perturb_surrounding_flow()
        if (perturb_sph) call s_perturb_sphere()
        if (instability_wave) call s_superposition_instability_wave()

        ! Converting the primitive variables to the conservative ones
        call s_convert_primitive_to_conservative_variables(q_prim_vf, &
                                                           q_cons_vf)

        if (qbmm .and. .not. polytropic) then
            !Initialize pb and mv
            call s_initialize_mv(q_cons_vf, mv%sf)
            call s_initialize_pb(q_cons_vf, mv%sf, pb%sf)
        end if

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
                    rand_real = rand_real*perturb_flow_mag
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

    !>  This subroutine computes velocity perturbations for a temporal mixing
        !!              layer with hypertangent mean streamwise velocity profile
        !!              obtained from linear stability analysis. For a 2D case,
        !!              instability waves with spatial wavenumbers, (4,0), (2,0),
        !!              and (1,0) are superposed. For a 3D waves, (4,4), (4,-4),
        !!              (2,2), (2,-2), (1,1), (1,-1) areadded on top of 2D waves.
    subroutine s_superposition_instability_wave() ! ------------------------
        real(kind(0d0)), dimension(5, 0:m, 0:n, 0:p) :: wave, wave1, wave2, wave_tmp
        real(kind(0d0)) :: tr, ti
        real(kind(0d0)) :: Lx, Lz
        integer :: i, j, k

        Lx = x_domain%end - x_domain%beg
        if (p > 0) then
            Lz = z_domain%end - z_domain%beg
        end if

        wave = 0d0
        wave1 = 0d0
        wave2 = 0d0

        ! Compute 2D waves
        call s_instability_wave(2*pi*4.0/Lx, 0d0, tr, ti, wave_tmp, 0d0)
        wave1 = wave1 + wave_tmp
        call s_instability_wave(2*pi*2.0/Lx, 0d0, tr, ti, wave_tmp, 0d0)
        wave1 = wave1 + wave_tmp
        call s_instability_wave(2*pi*1.0/Lx, 0d0, tr, ti, wave_tmp, 0d0)
        wave1 = wave1 + wave_tmp
        wave = wave1*0.05

        if (p > 0) then
            ! Compute 3D waves with phase shifts.
            call s_instability_wave(2*pi*4.0/Lx, 2*pi*4.0/Lz, tr, ti, wave_tmp, 2*pi*11d0/31d0)
            wave2 = wave2 + wave_tmp
            call s_instability_wave(2*pi*2.0/Lx, 2*pi*2.0/Lz, tr, ti, wave_tmp, 2*pi*13d0/31d0)
            wave2 = wave2 + wave_tmp
            call s_instability_wave(2*pi*1.0/Lx, 2*pi*1.0/Lz, tr, ti, wave_tmp, 2*pi*17d0/31d0)
            wave2 = wave2 + wave_tmp
            call s_instability_wave(2*pi*4.0/Lx, -2*pi*4.0/Lz, tr, ti, wave_tmp, 2*pi*19d0/31d0)
            wave2 = wave2 + wave_tmp
            call s_instability_wave(2*pi*2.0/Lx, -2*pi*2.0/Lz, tr, ti, wave_tmp, 2*pi*23d0/31d0)
            wave2 = wave2 + wave_tmp
            call s_instability_wave(2*pi*1.0/Lx, -2*pi*1.0/Lz, tr, ti, wave_tmp, 2*pi*29d0/31d0)
            wave2 = wave2 + wave_tmp
            wave = wave + 0.15*wave2
        end if

        ! Superpose velocity perturbuations (instability waves) to the velocity field
        do k = 0, p
            do j = 0, n
                do i = 0, m
                    q_prim_vf(mom_idx%beg)%sf(i, j, k) = q_prim_vf(mom_idx%beg)%sf(i, j, k) + wave(2, i, j, k)
                    q_prim_vf(mom_idx%beg + 1)%sf(i, j, k) = q_prim_vf(mom_idx%beg + 1)%sf(i, j, k) + wave(3, i, j, k)
                    if (p > 0) then
                        q_prim_vf(mom_idx%beg + 2)%sf(i, j, k) = q_prim_vf(mom_idx%beg + 2)%sf(i, j, k) + wave(4, i, j, k)
                    end if
                end do
            end do
        end do

    end subroutine s_superposition_instability_wave ! ----------------------

    !>  This subroutine computes instability waves for a given set of spatial
        !!              wavenumbers (alpha, beta) in x and z directions.
        !!              The eigenvalue problem is derived from the linearized
        !!              Euler equations with parallel mean flow assumption
        !!              (See Sandham 1989 PhD thesis for details).
    subroutine s_instability_wave(alpha, beta, tr, ti, wave, shift)
        real(kind(0d0)), intent(in) :: alpha, beta !<  spatial wavenumbers
        real(kind(0d0)), dimension(0:n) :: rho_mean, u_mean, t_mean !<  mean profiles
        real(kind(0d0)), dimension(0:n) :: drho_mean, du_mean, dt_mean !< y-derivatives of mean profiles
        real(kind(0d0)), dimension(0:n, 0:n) :: d !< differential operator in y dir
        real(kind(0d0)), dimension(0:5*(n + 1) - 1, 0:5*(n + 1) - 1) :: ar, ai, br, bi, ci !< matrices for eigenvalue problem
        real(kind(0d0)), dimension(0:5*(n + 1) - 1, 0:5*(n + 1) - 1) :: zr, zi !< eigenvectors
        real(kind(0d0)), dimension(0:5*(n + 1) - 1) :: wr, wi !< eigenvalues
        real(kind(0d0)), dimension(0:5*(n + 1) - 1) :: fv1, fv2, fv3 !< temporary memory
        real(kind(0d0)) :: tr, ti !< most unstable eigenvalue
        real(kind(0d0)), dimension(0:5*(n + 1) - 1) :: vr, vi, vnr, vni !< most unstable eigenvector and normalized one
        real(kind(0d0)), dimension(5, 0:m, 0:n, 0:p) :: wave !< instability wave
        real(kind(0d0)) :: shift !< phase shift
        real(kind(0d0)) :: gam, pi_inf, rho1, mach, c1
        integer :: ierr
        integer :: j, k, l !<  generic loop iterators
        integer :: ii, jj !< block matrix indices

        ! Set fluid flow properties
        gam = 1.+1./fluid_pp(1)%gamma
        pi_inf = fluid_pp(1)%pi_inf*(gam - 1.)/gam
        if (bubbles .and. num_fluids == 1) then
            rho1 = patch_icpp(1)%alpha_rho(1)/(1d0 - patch_icpp(1)%alpha(1))
        else
            rho1 = patch_icpp(1)%alpha_rho(1)/patch_icpp(1)%alpha(1)
        end if
        c1 = sqrt((gam*(patch_icpp(1)%pres + pi_inf))/rho1)
        mach = 1./c1

        ! Assign mean profiles
        do j = 0, n
            u_mean(j) = tanh(y_cc(j))
            t_mean(j) = 1 + 0.5*(gam - 1)*mach**2*(1 - u_mean(j)**2)
            rho_mean(j) = 1/T_mean(j)
        end do

        ! Compute differential operator in y-dir
        ! based on 4th order central difference (inner)
        ! and 2nd order central difference (near boundaries)
        dy = y_cc(1) - y_cc(0)
        d = 0d0
        d(1, 0) = -1/(2*dy)
        d(1, 2) = 1/(2*dy)
        do j = 2, n - 2
            d(j, j - 2) = 1/(12*dy)
            d(j, j - 1) = -8/(12*dy)
            d(j, j + 1) = 8/(12*dy)
            d(j, j + 2) = -1/(12*dy)
        end do
        d(n - 1, n - 2) = -1/(2*dy)
        d(n - 1, n) = 1/(2*dy)

        ! Compute y-derivatives of rho, u, T
        do j = 0, n
            drho_mean(j) = 0
            du_mean(j) = 0
            dt_mean(j) = 0
            do k = 0, n
                drho_mean(j) = drho_mean(j) + d(j, k)*rho_mean(k)
                du_mean(j) = du_mean(j) + d(j, k)*u_mean(k)
                dt_mean(j) = dt_mean(j) + d(j, k)*t_mean(k)
            end do
        end do

        ! Compute B and C, then A = B + C
        ! B includes terms without differential operator, and
        ! C includes terms with differential operator
        br = 0d0
        bi = 0d0
        ci = 0d0
        do j = 0, n
            ii = 1; jj = 1; br((ii - 1)*(n + 1) + j, (jj - 1)*(n + 1) + j) = alpha*u_mean(j); 
            ii = 1; jj = 2; br((ii - 1)*(n + 1) + j, (jj - 1)*(n + 1) + j) = alpha*rho_mean(j); 
            ii = 1; jj = 3; bi((ii - 1)*(n + 1) + j, (jj - 1)*(n + 1) + j) = -drho_mean(j); 
            ii = 1; jj = 4; br((ii - 1)*(n + 1) + j, (jj - 1)*(n + 1) + j) = beta*rho_mean(j); 
            ii = 2; jj = 1; br((ii - 1)*(n + 1) + j, (jj - 1)*(n + 1) + j) = alpha*t_mean(j)/(rho_mean(j)*gam*mach**2); 
            ii = 2; jj = 2; br((ii - 1)*(n + 1) + j, (jj - 1)*(n + 1) + j) = alpha*u_mean(j); 
            ii = 2; jj = 3; bi((ii - 1)*(n + 1) + j, (jj - 1)*(n + 1) + j) = -du_mean(j); 
            ii = 2; jj = 5; br((ii - 1)*(n + 1) + j, (jj - 1)*(n + 1) + j) = alpha/(gam*mach**2); 
            ii = 3; jj = 1; bi((ii - 1)*(n + 1) + j, (jj - 1)*(n + 1) + j) = -dt_mean(j)/(rho_mean(j)*gam*mach**2); 
            ii = 3; jj = 3; br((ii - 1)*(n + 1) + j, (jj - 1)*(n + 1) + j) = alpha*u_mean(j); 
            ii = 3; jj = 5; bi((ii - 1)*(n + 1) + j, (jj - 1)*(n + 1) + j) = -drho_mean(j)/(rho_mean(j)*gam*mach**2); 
            ii = 4; jj = 1; br((ii - 1)*(n + 1) + j, (jj - 1)*(n + 1) + j) = beta*t_mean(j)/(rho_mean(j)*gam*mach**2); 
            ii = 4; jj = 4; br((ii - 1)*(n + 1) + j, (jj - 1)*(n + 1) + j) = alpha*u_mean(j); 
            ii = 4; jj = 5; br((ii - 1)*(n + 1) + j, (jj - 1)*(n + 1) + j) = beta/(gam*mach**2); 
            ii = 5; jj = 2; br((ii - 1)*(n + 1) + j, (jj - 1)*(n + 1) + j) = (gam - 1)*alpha/rho_mean(j); 
            ii = 5; jj = 3; bi((ii - 1)*(n + 1) + j, (jj - 1)*(n + 1) + j) = -dt_mean(j); 
            ii = 5; jj = 4; br((ii - 1)*(n + 1) + j, (jj - 1)*(n + 1) + j) = (gam - 1)*beta/rho_mean(j); 
            ii = 5; jj = 5; br((ii - 1)*(n + 1) + j, (jj - 1)*(n + 1) + j) = alpha*u_mean(j); 
            do k = 0, n
                ii = 1; jj = 3; ci((ii - 1)*(n + 1) + j, (jj - 1)*(n + 1) + k) = -rho_mean(j)*d(j, k); 
                ii = 3; jj = 1; ci((ii - 1)*(n + 1) + j, (jj - 1)*(n + 1) + k) = -t_mean(j)*d(j, k)/(rho_mean(j)*gam*mach**2); 
                ii = 3; jj = 5; ci((ii - 1)*(n + 1) + j, (jj - 1)*(n + 1) + k) = -d(j, k)/(gam*mach**2); 
                ii = 5; jj = 3; ci((ii - 1)*(n + 1) + j, (jj - 1)*(n + 1) + k) = -(gam - 1)*d(j, k)/rho_mean(j); 
            end do
        end do
        ar = br
        ai = bi + ci

        ! Compute eigenvalues and eigenvectors
        call cg(5*(n + 1), 5*(n + 1), ar, ai, wr, wi, zr, zi, fv1, fv2, fv3, ierr)

        ! Generate instability wave
        call s_generate_wave(5*(n + 1), wr, wi, zr, zi, alpha, beta, wave, shift)

    end subroutine s_instability_wave

    !>  This subroutine generates an instability wave using the most unstable
        !!              eigenvalue and corresponding eigenvector among the
        !!              given set of eigenvalues and eigenvectors.
    subroutine s_generate_wave(nl, wr, wi, zr, zi, alpha, beta, wave, shift)
        integer nl
        real(kind(0d0)), dimension(0:nl - 1) :: wr, wi !< eigenvalues
        real(kind(0d0)), dimension(0:nl - 1, 0:nl - 1) :: zr, zi !< eigenvectors
        real(kind(0d0)), dimension(0:nl - 1) :: vr, vi, vnr, vni !< most unstable eigenvector
        real(kind(0d0)), dimension(5, 0:m, 0:n, 0:p) :: wave
        real(kind(0d0)) :: alpha, beta, ang, shift
        real(kind(0d0)) :: norm
        real(kind(0d0)) :: tr, ti, cr, ci !< temporary memory
        integer idx
        integer i, j, k

        ! Find the most unstable eigenvalue and corresponding eigenvector
        k = 0
        do i = 1, nl - 1
            if (wi(i) > wi(k)) then
                k = i
            end if
        end do
        vr = zr(:, k)
        vi = zi(:, k)

        ! Normalize the eigenvector by its component with the largest modulus.
        norm = 0d0
        do i = 0, nl - 1
            if (dsqrt(vr(i)**2 + vi(i)**2) > norm) then
                idx = i
                norm = dsqrt(vr(i)**2 + vi(i)**2)
            end if
        end do

        tr = vr(idx)
        ti = vi(idx)
        do i = 0, nl - 1
            call cdiv(vr(i), vi(i), tr, ti, cr, ci)
            vnr(i) = cr
            vni(i) = ci
        end do

        ! Generate an instability wave
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    if (beta == 0) then
                        ang = alpha*x_cc(i)
                    else
                        ang = alpha*x_cc(i) + beta*z_cc(k) + shift
                    end if
                    wave(1, i, j, k) = vnr(j)*cos(ang) - vni(j)*sin(ang)                 ! rho
                    wave(2, i, j, k) = vnr((n + 1) + j)*cos(ang) - vni((n + 1) + j)*sin(ang)     ! u
                    wave(3, i, j, k) = vnr(2*(n + 1) + j)*cos(ang) - vni(2*(n + 1) + j)*sin(ang) ! v
                    wave(4, i, j, k) = vnr(3*(n + 1) + j)*cos(ang) - vni(3*(n + 1) + j)*sin(ang) ! w
                    wave(5, i, j, k) = vnr(4*(n + 1) + j)*cos(ang) - vni(4*(n + 1) + j)*sin(ang) ! T
                end do
            end do
        end do

    end subroutine s_generate_wave

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
        deallocate (ib_markers%sf)

    end subroutine s_finalize_initial_condition_module ! -------------------

end module m_initial_condition
