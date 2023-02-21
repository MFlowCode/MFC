!>
!! @file m_global_parameters.f90
!! @brief Contains module m_global_parameters

!> @brief This module contains all of the parameters characterizing the
!!              computational domain, simulation algorithm, initial condition
!!              and the stiffened equation of state.
module m_global_parameters

    ! Dependencies =============================================================
#ifdef MFC_MPI
    use mpi                     ! Message passing interface (MPI) module
#endif

    use m_derived_types         ! Definitions of the derived types

    ! ==========================================================================

    implicit none

    ! Logistics ================================================================
    integer :: num_procs            !< Number of processors
    character(LEN=path_len) :: case_dir             !< Case folder location
    logical :: old_grid             !< Use existing grid data
    logical :: old_ic               !< Use existing IC data
    integer :: t_step_old           !< Existing IC/grid folder
    ! ==========================================================================

    ! Computational Domain Parameters ==========================================

    integer :: proc_rank !< Rank of the local processor

    integer :: m
    integer :: n
    integer :: p !<
    !! Number of cells in the x-, y- and z-coordinate directions

    integer :: m_glb, n_glb, p_glb !<
    !! Global number of cells in each direction

    integer :: num_dims !< Number of spatial dimensions

    logical :: cyl_coord
    integer :: grid_geometry !<
    !! Cylindrical coordinates (either axisymmetric or full 3D)

    real(kind(0d0)), allocatable, dimension(:) :: x_cc, y_cc, z_cc !<
    !! Locations of cell-centers (cc) in x-, y- and z-directions, respectively

    real(kind(0d0)), allocatable, dimension(:) :: x_cb, y_cb, z_cb !<
    !! Locations of cell-boundaries (cb) in x-, y- and z-directions, respectively

    real(kind(0d0)) :: dx, dy, dz !<
    !! Minimum cell-widths in the x-, y- and z-coordinate directions

    type(bounds_info) :: x_domain, y_domain, z_domain !<
    !! Locations of the domain bounds in the x-, y- and z-coordinate directions

    logical :: stretch_x, stretch_y, stretch_z !<
    !! Grid stretching flags for the x-, y- and z-coordinate directions

    ! Parameters of the grid stretching function for the x-, y- and z-coordinate
    ! directions. The "a" parameters are a measure of the rate at which the grid
    ! is stretched while the remaining parameters are indicative of the location
    ! on the grid at which the stretching begins.
    real(kind(0d0)) :: a_x, a_y, a_z
    integer :: loops_x, loops_y, loops_z
    real(kind(0d0)) :: x_a, y_a, z_a
    real(kind(0d0)) :: x_b, y_b, z_b

    ! ==========================================================================

    ! Simulation Algorithm Parameters ==========================================
    integer :: model_eqns      !< Multicomponent flow model
    integer :: num_fluids      !< Number of different fluids present in the flow
    logical :: adv_alphan      !< Advection of the last volume fraction
    logical :: mpp_lim         !< Alpha limiter
    integer :: sys_size        !< Number of unknowns in the system of equations
    integer :: weno_order      !< Order of accuracy for the WENO reconstruction
    logical :: hypoelasticity  !< activate hypoelasticity

    ! Annotations of the structure, i.e. the organization, of the state vectors
    type(int_bounds_info) :: cont_idx                   !< Indexes of first & last continuity eqns.
    type(int_bounds_info) :: mom_idx                    !< Indexes of first & last momentum eqns.
    integer :: E_idx                      !< Index of total energy equation
    integer :: alf_idx                    !< Index of void fraction
    type(int_bounds_info) :: adv_idx                    !< Indexes of first & last advection eqns.
    type(int_bounds_info) :: internalEnergies_idx       !< Indexes of first & last internal energy eqns.
    type(bub_bounds_info) :: bub_idx                    !< Indexes of first & last bubble variable eqns.
    integer :: gamma_idx                  !< Index of specific heat ratio func. eqn.
    integer :: pi_inf_idx                 !< Index of liquid stiffness func. eqn.
    type(int_bounds_info) :: stress_idx                 !< Indexes of elastic shear stress eqns.

    type(int_bounds_info) :: bc_x, bc_y, bc_z !<
    !! Boundary conditions in the x-, y- and z-coordinate directions

    logical :: parallel_io !< Format of the data files
    integer :: precision !< Precision of output files

    ! Perturb density of surrounding air so as to break symmetry of grid
    logical :: perturb_flow
    integer :: perturb_flow_fluid   !< Fluid to be perturbed with perturb_flow flag
    logical :: perturb_sph
    integer :: perturb_sph_fluid    !< Fluid to be perturbed with perturb_sph flag
    real(kind(0d0)), dimension(num_fluids_max) :: fluid_rho

    integer, allocatable, dimension(:) :: proc_coords !<
    !! Processor coordinates in MPI_CART_COMM

    integer, allocatable, dimension(:) :: start_idx !<
    !! Starting cell-center index of local processor in global grid

#ifdef MFC_MPI

    type(mpi_io_var), public :: MPI_IO_DATA

    character(LEN=name_len) :: mpiiofs
    integer :: mpi_info_int !<
    !! MPI info for parallel IO with Lustre file systems

#endif

    integer, private :: ierr
    ! ==========================================================================

    ! Initial Condition Parameters =============================================
    integer :: num_patches     !< Number of patches composing initial condition

    type(ic_patch_parameters), dimension(num_patches_max) :: patch_icpp !<
    !! Database of the initial condition patch parameters (icpp) for each of the
    !! patches employed in the configuration of the initial condition. Note that
    !! the maximum allowable number of patches, num_patches_max, may be changed
    !! in the module m_derived_types.f90.
    ! ==========================================================================

    ! Fluids Physical Parameters ===============================================
    type(physical_parameters), dimension(num_fluids_max) :: fluid_pp !<
    !! Database of the physical parameters of each of the fluids that is present
    !! in the flow. These include the stiffened gas equation of state parameters,
    !! the Reynolds numbers and the Weber numbers.

    ! ==========================================================================

    real(kind(0d0)) :: rhoref, pref !< Reference parameters for Tait EOS

    !> @name Bubble modeling
    !> @{
    integer :: nb
    real(kind(0d0)) :: R0ref
    real(kind(0d0)) :: Ca, Web, Re_inv
    real(kind(0d0)), dimension(:), allocatable :: weight, R0, V0
    logical :: bubbles
    logical :: qbmm      !< Quadrature moment method
    integer :: nmom  !< Number of carried moments
    real(kind(0d0)) :: sigR, sigV, rhoRV !< standard deviations in R/V
    !> @}

    !> @name Non-polytropic bubble gas compression
    !> @{
    logical :: polytropic
    logical :: polydisperse
    integer :: thermal  !1 = adiabatic, 2 = isotherm, 3 = transfer
    real(kind(0d0)) :: R_n, R_v, phi_vn, phi_nv, Pe_c, Tw
    real(kind(0d0)), dimension(:), allocatable :: k_n, k_v, pb0, mass_n0, mass_v0, Pe_T
    real(kind(0d0)), dimension(:), allocatable :: Re_trans_T, Re_trans_c, Im_trans_T, Im_trans_c, omegaN
    real(kind(0d0)) :: poly_sigma
    integer :: dist_type !1 = binormal, 2 = lognormal-normal
    integer :: R0_type   !1 = simpson
    !> @}

    !> @name Index variables used for m_variables_conversion
    !> @{
    integer :: momxb, momxe
    integer :: advxb, advxe
    integer :: contxb, contxe
    integer :: intxb, intxe
    integer :: bubxb, bubxe
    integer :: strxb, strxe
    !> @}

    integer, allocatable, dimension(:, :, :) :: logic_grid


contains

    !>  Assigns default values to user inputs prior to reading
        !!              them in. This allows for an easier consistency check of
        !!              these parameters once they are read from the input file.
    subroutine s_assign_default_values_to_user_inputs() ! ------------------

        integer :: i !< Generic loop operator

        ! Logistics
        case_dir = '.'
        old_grid = .false.
        old_ic = .false.
        t_step_old = dflt_int

        ! Computational domain parameters
        m = dflt_int; n = 0; p = 0

        cyl_coord = .false.

        x_domain%beg = dflt_real
        x_domain%end = dflt_real
        y_domain%beg = dflt_real
        y_domain%end = dflt_real
        z_domain%beg = dflt_real
        z_domain%end = dflt_real

        stretch_x = .false.
        stretch_y = .false.
        stretch_z = .false.

        a_x = dflt_real
        a_y = dflt_real
        a_z = dflt_real
        loops_x = 1
        loops_y = 1
        loops_z = 1
        x_a = dflt_real
        x_b = dflt_real
        y_a = dflt_real
        y_b = dflt_real
        z_a = dflt_real
        z_b = dflt_real

        ! Simulation algorithm parameters
        model_eqns = dflt_int
        num_fluids = dflt_int
        adv_alphan = .false.
        weno_order = dflt_int

        hypoelasticity = .false.

        bc_x%beg = dflt_int
        bc_x%end = dflt_int
        bc_y%beg = dflt_int
        bc_y%end = dflt_int
        bc_z%beg = dflt_int
        bc_z%end = dflt_int

        parallel_io = .false.
        precision = 2
        perturb_flow = .false.
        perturb_flow_fluid = dflt_int
        perturb_sph = .false.
        perturb_sph_fluid = dflt_int
        fluid_rho = dflt_real

        ! Initial condition parameters
        num_patches = dflt_int

        do i = 1, num_patches_max
            patch_icpp(i)%geometry = dflt_int
            patch_icpp(i)%x_centroid = dflt_real
            patch_icpp(i)%y_centroid = dflt_real
            patch_icpp(i)%z_centroid = dflt_real
            patch_icpp(i)%length_x = dflt_real
            patch_icpp(i)%length_y = dflt_real
            patch_icpp(i)%length_z = dflt_real
            patch_icpp(i)%radius = dflt_real
            patch_icpp(i)%epsilon = dflt_real
            patch_icpp(i)%beta = dflt_real
            patch_icpp(i)%normal = dflt_real
            patch_icpp(i)%radii = dflt_real
            patch_icpp(i)%alter_patch = .false.
            patch_icpp(i)%alter_patch(0) = .true.
            patch_icpp(i)%smoothen = .false.
            patch_icpp(i)%smooth_patch_id = i
            patch_icpp(i)%smooth_coeff = dflt_real
            patch_icpp(i)%alpha_rho = dflt_real
            patch_icpp(i)%rho = dflt_real
            patch_icpp(i)%vel = dflt_real
            patch_icpp(i)%pres = dflt_real
            patch_icpp(i)%alpha = dflt_real
            patch_icpp(i)%gamma = dflt_real
            patch_icpp(i)%pi_inf = dflt_real
            patch_icpp(i)%tau_e = 0d0
            !should get all of r0's and v0's
            patch_icpp(i)%r0 = dflt_real
            patch_icpp(i)%v0 = dflt_real

            patch_icpp(i)%p0 = dflt_real
            patch_icpp(i)%m0 = dflt_real
        end do

        ! Tait EOS
        rhoref = dflt_real
        pref = dflt_real

        ! Bubble modeling
        bubbles = .false.
        polytropic = .true.
        polydisperse = .false.

        thermal = dflt_int
        R0ref = dflt_real
        nb = dflt_int

        Ca = dflt_real
        Re_inv = dflt_real
        Web = dflt_real
        poly_sigma = dflt_real

        qbmm = .false.
        nmom = 1
        sigR = dflt_real
        sigV = dflt_real
        rhoRV = 0d0
        dist_type = dflt_int
        R0_type = dflt_int

        R_n = dflt_real
        R_v = dflt_real
        phi_vn = dflt_real
        phi_nv = dflt_real
        Pe_c = dflt_real
        Tw = dflt_real

        ! Fluids physical parameters
        do i = 1, num_fluids_max
            fluid_pp(i)%gamma = dflt_real
            fluid_pp(i)%pi_inf = dflt_real
            fluid_pp(i)%mul0 = dflt_real
            fluid_pp(i)%ss = dflt_real
            fluid_pp(i)%pv = dflt_real
            fluid_pp(i)%gamma_v = dflt_real
            fluid_pp(i)%M_v = dflt_real
            fluid_pp(i)%mu_v = dflt_real
            fluid_pp(i)%k_v = dflt_real
            fluid_pp(i)%G = 0d0
        end do

    end subroutine s_assign_default_values_to_user_inputs ! ----------------

    !> Computation of parameters, allocation procedures, and/or
        !! any other tasks needed to properly setup the module
    subroutine s_initialize_global_parameters_module() ! ----------------------

        integer :: i, j, fac

        ! Determining the layout of the state vectors and overall size of
        ! the system of equations, given the dimensionality and choice of
        ! the equations of motion

        ! Gamma/Pi_inf Model ===============================================
        if (model_eqns == 1) then

            ! Setting number of fluids
            num_fluids = 1

            ! Annotating structure of the state and flux vectors belonging
            ! to the system of equations defined by the selected number of
            ! spatial dimensions and the gamma/pi_inf model
            cont_idx%beg = 1
            cont_idx%end = cont_idx%beg
            mom_idx%beg = cont_idx%end + 1
            mom_idx%end = cont_idx%end + num_dims
            E_idx = mom_idx%end + 1
            adv_idx%beg = E_idx + 1
            adv_idx%end = adv_idx%beg + 1
            gamma_idx = adv_idx%beg
            pi_inf_idx = adv_idx%end
            sys_size = adv_idx%end

            ! ==================================================================

            ! Volume Fraction Model (5-equation model) =========================
        else if (model_eqns == 2) then

            ! Annotating structure of the state and flux vectors belonging
            ! to the system of equations defined by the selected number of
            ! spatial dimensions and the volume fraction model
            cont_idx%beg = 1
            cont_idx%end = num_fluids
            mom_idx%beg = cont_idx%end + 1
            mom_idx%end = cont_idx%end + num_dims
            E_idx = mom_idx%end + 1
            adv_idx%beg = E_idx + 1
            adv_idx%end = E_idx + num_fluids

            sys_size = adv_idx%end

            if (bubbles) then
                alf_idx = adv_idx%end
            else
                alf_idx = 1
            end if

            if (bubbles) then
                bub_idx%beg = sys_size + 1
                if (qbmm) then
                    if (nnode == 4) then
                        nmom = 6 !! Already set as a parameter
                    end if
                    bub_idx%end = adv_idx%end + nb*nmom
                else
                    if (.not. polytropic) then
                        bub_idx%end = sys_size + 4*nb
                    else
                        bub_idx%end = sys_size + 2*nb
                    end if
                end if
                sys_size = bub_idx%end

                allocate (weight(nb), R0(nb), V0(nb))
                allocate (bub_idx%rs(nb), bub_idx%vs(nb))
                allocate (bub_idx%ps(nb), bub_idx%ms(nb))

                if (qbmm) then
                    allocate (bub_idx%moms(nb, nmom))
                    allocate (bub_idx%fullmom(nb, 0:nmom, 0:nmom))

                    do i = 1, nb
                        do j = 1, nmom
                            bub_idx%moms(i, j) = bub_idx%beg + (j - 1) + (i - 1)*nmom
                        end do
                        bub_idx%fullmom(i, 0, 0) = bub_idx%moms(i, 1)
                        bub_idx%fullmom(i, 1, 0) = bub_idx%moms(i, 2)
                        bub_idx%fullmom(i, 0, 1) = bub_idx%moms(i, 3)
                        bub_idx%fullmom(i, 2, 0) = bub_idx%moms(i, 4)
                        bub_idx%fullmom(i, 1, 1) = bub_idx%moms(i, 5)
                        bub_idx%fullmom(i, 0, 2) = bub_idx%moms(i, 6)
                        bub_idx%rs(i) = bub_idx%fullmom(i, 1, 0)
                    end do
                else
                    do i = 1, nb
                        if (.not. polytropic) then
                            fac = 4
                        else
                            fac = 2
                        end if

                        bub_idx%rs(i) = bub_idx%beg + (i - 1)*fac
                        bub_idx%vs(i) = bub_idx%rs(i) + 1

                        if (.not. polytropic) then
                            bub_idx%ps(i) = bub_idx%vs(i) + 1
                            bub_idx%ms(i) = bub_idx%ps(i) + 1
                        end if
                    end do
                end if

                if (nb == 1) then
                    weight(:) = 1d0
                    R0(:) = 1d0
                    V0(:) = 1d0
                else if (nb > 1) then
                    if (R0_type == 1) then
                        call s_simpson
                    else
                        print *, 'Invalid R0 type - abort'
                        stop
                    end if
                    V0(:) = 1d0
                else
                    stop 'Invalid value of nb'
                end if

                print *, 'R0 weights: ', weight(:)
                print *, 'R0 abscissas: ', R0(:)

                if (.not. polytropic) then
                    call s_initialize_nonpoly
                else
                    rhoref = 1.d0
                    pref = 1.d0
                end if
            end if

            if (hypoelasticity) then
                stress_idx%beg = sys_size + 1
                stress_idx%end = sys_size + (num_dims*(num_dims + 1))/2
                ! number of stresses is 1 in 1D, 3 in 2D, 6 in 3D
                sys_size = stress_idx%end
            end if

            ! ==================================================================

            ! Volume Fraction Model (6-equation model) =========================
        else if (model_eqns == 3) then

            ! Annotating structure of the state and flux vectors belonging
            ! to the system of equations defined by the selected number of
            ! spatial dimensions and the volume fraction model
            cont_idx%beg = 1
            cont_idx%end = num_fluids
            mom_idx%beg = cont_idx%end + 1
            mom_idx%end = cont_idx%end + num_dims
            E_idx = mom_idx%end + 1
            adv_idx%beg = E_idx + 1
            adv_idx%end = E_idx + num_fluids
            internalEnergies_idx%beg = adv_idx%end + 1
            internalEnergies_idx%end = adv_idx%end + num_fluids
            sys_size = internalEnergies_idx%end
            !========================
        else if (model_eqns == 4) then
            ! 4 equation model with subgrid bubbles
            cont_idx%beg = 1 ! one continuity equation
            cont_idx%end = 1 ! num_fluids
            mom_idx%beg = cont_idx%end + 1 ! one momentum equation in each direction
            mom_idx%end = cont_idx%end + num_dims
            E_idx = mom_idx%end + 1 ! one energy equation
            adv_idx%beg = E_idx + 1
            adv_idx%end = adv_idx%beg !one volume advection equation
            alf_idx = adv_idx%end
            sys_size = alf_idx !adv_idx%end

            if (bubbles) then
                bub_idx%beg = sys_size + 1
                bub_idx%end = sys_size + 2*nb
                if (.not. polytropic) then
                    bub_idx%end = sys_size + 4*nb
                end if
                sys_size = bub_idx%end

                allocate (bub_idx%rs(nb), bub_idx%vs(nb))
                allocate (bub_idx%ps(nb), bub_idx%ms(nb))
                allocate (weight(nb), R0(nb), V0(nb))

                do i = 1, nb
                    if (.not. polytropic) then
                        fac = 4
                    else
                        fac = 2
                    end if

                    bub_idx%rs(i) = bub_idx%beg + (i - 1)*fac
                    bub_idx%vs(i) = bub_idx%rs(i) + 1

                    if (.not. polytropic) then
                        bub_idx%ps(i) = bub_idx%vs(i) + 1
                        bub_idx%ms(i) = bub_idx%ps(i) + 1
                    end if
                end do

                if (nb == 1) then
                    weight(:) = 1d0
                    R0(:) = 1d0
                    V0(:) = 0d0
                else if (nb > 1) then
                    if (R0_type == 1) then
                        call s_simpson
                    else
                        print *, 'Invalid R0 type - abort'
                        stop
                    end if
                    V0(:) = 1d0
                else
                    stop 'Invalid value of nb'
                end if

                if (.not. polytropic) then
                    call s_initialize_nonpoly
                else
                    rhoref = 1.d0
                    pref = 1.d0
                end if
            end if
        end if

        momxb = mom_idx%beg
        momxe = mom_idx%end
        advxb = adv_idx%beg
        advxe = adv_idx%end
        contxb = cont_idx%beg
        contxe = cont_idx%end
        bubxb = bub_idx%beg
        bubxe = bub_idx%end
        strxb = stress_idx%beg
        strxe = stress_idx%end
        intxb = internalEnergies_idx%beg
        intxe = internalEnergies_idx%end

        ! ==================================================================

#ifdef MFC_MPI

        allocate (MPI_IO_DATA%view(1:sys_size))
        allocate (MPI_IO_DATA%var(1:sys_size))

        do i = 1, sys_size
            allocate (MPI_IO_DATA%var(i)%sf(0:m, 0:n, 0:p))
            MPI_IO_DATA%var(i)%sf => null()
        end do

#endif

        ! Allocating grid variables for the x-direction
        allocate (x_cc(0:m), x_cb(-1:m))
        ! Allocating grid variables for the y- and z-directions
        if (n > 0) then
            allocate (y_cc(0:n), y_cb(-1:n))
            if (p > 0) then
                allocate (z_cc(0:p), z_cb(-1:p))
            end if
        end if

        if (cyl_coord .neqv. .true.) then ! Cartesian grid
            grid_geometry = 1
        elseif (cyl_coord .and. p == 0) then ! Axisymmetric cylindrical grid
            grid_geometry = 2
        else ! Fully 3D cylindrical grid
            grid_geometry = 3
        end if

        allocate (logic_grid(0:m, 0:n, 0:p))

    end subroutine s_initialize_global_parameters_module ! --------------------

    !> Initializes and computes bubble properties
        !! for non-polytropic processes
    subroutine s_initialize_nonpoly
        integer :: ir
        real(kind(0.d0)) :: rhol0
        real(kind(0.d0)) :: pl0
        real(kind(0.d0)) :: uu
        real(kind(0.d0)) :: D_m
        real(kind(0.d0)) :: temp
        real(kind(0.d0)) :: omega_ref
        real(kind(0.d0)), dimension(Nb) :: chi_vw0
        real(kind(0.d0)), dimension(Nb) :: cp_m0
        real(kind(0.d0)), dimension(Nb) :: k_m0
        real(kind(0.d0)), dimension(Nb) :: rho_m0
        real(kind(0.d0)), dimension(Nb) :: x_vw
        ! polytropic index used to compute isothermal natural frequency
        real(kind(0.d0)), parameter :: k_poly = 1.d0
        ! universal gas constant
        real(kind(0.d0)), parameter :: Ru = 8314.d0

        ! liquid physical properties
        real(kind(0.d0)) :: mul0, ss, pv, gamma_v, M_v, mu_v

        ! gas physical properties
        real(kind(0.d0)) :: gamma_m, gamma_n, M_n, mu_n

        rhol0 = rhoref
        pl0 = pref

        allocate (pb0(nb), mass_n0(nb), mass_v0(nb), Pe_T(nb))
        allocate (k_n(nb), k_v(nb), omegaN(nb))
        allocate (Re_trans_T(nb), Re_trans_c(nb), Im_trans_T(nb), Im_trans_c(nb))

        pb0(:) = dflt_real
        mass_n0(:) = dflt_real
        mass_v0(:) = dflt_real
        Pe_T(:) = dflt_real
        omegaN(:) = dflt_real

        mul0 = fluid_pp(1)%mul0
        ss = fluid_pp(1)%ss
        pv = fluid_pp(1)%pv
        gamma_v = fluid_pp(1)%gamma_v
        M_v = fluid_pp(1)%M_v
        mu_v = fluid_pp(1)%mu_v
        k_v(:) = fluid_pp(1)%k_v

        gamma_n = fluid_pp(2)%gamma_v
        M_n = fluid_pp(2)%M_v
        mu_n = fluid_pp(2)%mu_v
        k_n(:) = fluid_pp(2)%k_v

        gamma_m = gamma_n
        if (thermal == 2) gamma_m = 1.d0 !isothermal

        temp = 293.15d0
        D_m = 0.242d-4
        uu = DSQRT(pl0/rhol0)

        omega_ref = 3.d0*k_poly*Ca + 2.d0*(3.d0*k_poly - 1.d0)/Web

            !!! thermal properties !!!
        ! gas constants
        R_n = Ru/M_n
        R_v = Ru/M_v
        ! phi_vn & phi_nv (phi_nn = phi_vv = 1)
        phi_vn = (1.d0 + DSQRT(mu_v/mu_n)*(M_n/M_v)**(0.25d0))**2 &
                 /(DSQRT(8.d0)*DSQRT(1.d0 + M_v/M_n))
        phi_nv = (1.d0 + DSQRT(mu_n/mu_v)*(M_v/M_n)**(0.25d0))**2 &
                 /(DSQRT(8.d0)*DSQRT(1.d0 + M_n/M_v))
        ! internal bubble pressure
        pb0 = pl0 + 2.d0*ss/(R0ref*R0)

        ! mass fraction of vapor
        chi_vw0 = 1.d0/(1.d0 + R_v/R_n*(pb0/pv - 1.d0))
        ! specific heat for gas/vapor mixture
        cp_m0 = chi_vw0*R_v*gamma_v/(gamma_v - 1.d0) &
                + (1.d0 - chi_vw0)*R_n*gamma_n/(gamma_n - 1.d0)
        ! mole fraction of vapor
        x_vw = M_n*chi_vw0/(M_v + (M_n - M_v)*chi_vw0)
        ! thermal conductivity for gas/vapor mixture
        k_m0 = x_vw*k_v/(x_vw + (1.d0 - x_vw)*phi_vn) &
               + (1.d0 - x_vw)*k_n/(x_vw*phi_nv + 1.d0 - x_vw)
        ! mixture density
        rho_m0 = pv/(chi_vw0*R_v*temp)

        ! mass of gas/vapor computed using dimensional quantities
        mass_n0 = 4.d0*(pb0 - pv)*pi/(3.d0*R_n*temp*rhol0)*R0**3
        mass_v0 = 4.d0*pv*pi/(3.d0*R_v*temp*rhol0)*R0**3
        ! Peclet numbers
        Pe_T = rho_m0*cp_m0*uu*R0ref/k_m0
        Pe_c = uu*R0ref/D_m
        ! nondimensional properties
        R_n = rhol0*R_n*temp/pl0
        R_v = rhol0*R_v*temp/pl0
        k_n = k_n/k_m0
        k_v = k_v/k_m0
        pb0 = pb0/pl0
        pv = pv/pl0

        print *, 'pb0 nondim/final', pb0

        ! bubble wall temperature, normalized by T0, in the liquid
        ! keeps a constant (cold liquid assumption)
        Tw = 1.d0
        ! natural frequencies
        omegaN = DSQRT(3.d0*k_poly*Ca + 2.d0*(3.d0*k_poly - 1.d0)/(Web*R0))/R0

        pl0 = 1.d0
        do ir = 1, Nb
            call s_transcoeff(omegaN(ir)*R0(ir), Pe_T(ir)*R0(ir), &
                              Re_trans_T(ir), Im_trans_T(ir))
            call s_transcoeff(omegaN(ir)*R0(ir), Pe_c*R0(ir), &
                              Re_trans_c(ir), Im_trans_c(ir))
        end do
        Im_trans_T = 0d0
        Im_trans_c = 0d0

        rhoref = 1.d0
        pref = 1.d0
    end subroutine s_initialize_nonpoly

    !> Computes the transfer coefficient for the non-polytropic bubble compression process
        !! @param omega natural frqeuencies
        !! @param peclet Peclet number
        !! @param Re_trans Real part of the transport coefficients
        !! @param Im_trans Imaginary part of the transport coefficients
    subroutine s_transcoeff(omega, peclet, Re_trans, Im_trans)

        real(kind(0.d0)), intent(IN) :: omega
        real(kind(0.d0)), intent(IN) :: peclet
        real(kind(0.d0)), intent(OUT) :: Re_trans
        real(kind(0.d0)), intent(OUT) :: Im_trans
        complex :: trans, c1, c2, c3
        complex :: imag = (0., 1.)
        real(kind(0.d0)) :: f_transcoeff

        c1 = imag*omega*peclet
        c2 = CSQRT(c1)
        c3 = (CEXP(c2) - CEXP(-c2))/(CEXP(c2) + CEXP(-c2)) ! TANH(c2)
        trans = ((c2/c3 - 1.d0)**(-1) - 3.d0/c1)**(-1) ! transfer function

        Re_trans = dble(trans)
        Im_trans = aimag(trans)

    end subroutine s_transcoeff

    subroutine s_initialize_parallel_io() ! --------------------------------

        num_dims = 1 + min(1, n) + min(1, p)

        allocate (proc_coords(1:num_dims))

        if (parallel_io .neqv. .true.) return

#ifdef MFC_MPI

        ! Option for Lustre file system (Darter/Comet/Stampede)
        write (mpiiofs, '(A)') '/lustre_'
        mpiiofs = trim(mpiiofs)
        call MPI_INFO_CREATE(mpi_info_int, ierr)
        call MPI_INFO_SET(mpi_info_int, 'romio_ds_write', 'disable', ierr)

        ! Option for UNIX file system (Hooke/Thomson)
        ! WRITE(mpiiofs, '(A)') '/ufs_'
        ! mpiiofs = TRIM(mpiiofs)
        ! mpi_info_int = MPI_INFO_NULL

        allocate (start_idx(1:num_dims))

#endif

    end subroutine s_initialize_parallel_io ! ------------------------------

    subroutine s_finalize_global_parameters_module() ! ------------------------

        integer :: i

        ! Deallocating grid variables for the x-direction
        deallocate (x_cc, x_cb)
        ! Deallocating grid variables for the y- and z-directions
        if (n > 0) then
            deallocate (y_cc, y_cb)
            if (p > 0) then
                deallocate (z_cc, z_cb)
            end if
        end if

        deallocate (proc_coords)

#ifdef MFC_MPI

        if (parallel_io) then
            deallocate (start_idx)
            do i = 1, sys_size
                MPI_IO_DATA%var(i)%sf => null()
            end do

            deallocate (MPI_IO_DATA%var)
            deallocate (MPI_IO_DATA%view)
        end if

#endif

    end subroutine s_finalize_global_parameters_module ! ----------------------

    !> Computes the bubble number density n from the conservative variables
        !! @param vftmp is the void fraction
        !! @param nRtmp is the bubble number  density times the bubble radii
        !! @param ntmp is the output number bubble density
    subroutine s_comp_n_from_cons(vftmp, nRtmp, ntmp)
        real(kind(0.d0)), intent(IN) :: vftmp
        real(kind(0.d0)), dimension(nb), intent(IN) :: nRtmp
        real(kind(0.d0)), intent(OUT) :: ntmp
        real(kind(0.d0)) :: nR3

        call s_quad(nRtmp**3.d0, nR3)  !returns itself if NR0 = 1
        ntmp = DSQRT((4.d0*pi/3.d0)*nR3/vftmp)
        ! ntmp = 1d0

    end subroutine s_comp_n_from_cons

    !> Computes the bubble number density n from the primitive variables
        !! @param vftmp is the void fraction
        !! @param Rtmp is the  bubble radii
        !! @param ntmp is the output number bubble density
    subroutine s_comp_n_from_prim(vftmp, Rtmp, ntmp)
        real(kind(0.d0)), intent(IN) :: vftmp
        real(kind(0.d0)), dimension(nb), intent(IN) :: Rtmp
        real(kind(0.d0)), intent(OUT) :: ntmp
        real(kind(0.d0)) :: R3

        call s_quad(Rtmp**3.d0, R3)  !returns itself if NR0 = 1
        ntmp = (3.d0/(4.d0*pi))*vftmp/R3
        ! ntmp = 1d0

    end subroutine s_comp_n_from_prim

    !> Computes the quadrature for polydisperse bubble populations
        !! @param func is the bubble dynamic variables for each bin
        !! @param mom is the computed moment
    subroutine s_quad(func, mom)

        real(kind(0.d0)), dimension(nb), intent(IN) :: func 
        real(kind(0.d0)), intent(OUT) :: mom

        mom = dot_product(weight, func)

    end subroutine s_quad

    !> Computes the Simpson weights for quadrature
    subroutine s_simpson

        integer :: ir
        real(kind(0.d0)) :: R0mn
        real(kind(0.d0)) :: R0mx
        real(kind(0.d0)) :: dphi
        real(kind(0.d0)) :: tmp
        real(kind(0.d0)) :: sd
        real(kind(0.d0)), dimension(nb) :: phi

        ! nondiml. min. & max. initial radii for numerical quadrature
        !sd   = 0.05D0
        !R0mn = 0.75D0
        !R0mx = 1.3D0

        !sd   = 0.3D0
        !R0mn = 0.3D0
        !R0mx = 6.D0

        !sd   = 0.7D0
        !R0mn = 0.12D0
        !R0mx = 150.D0

        sd = poly_sigma
        R0mn = 0.8d0*DEXP(-2.8d0*sd)
        R0mx = 0.2d0*DEXP(9.5d0*sd) + 1.d0

        ! phi = ln( R0 ) & return R0
        do ir = 1, nb
            phi(ir) = DLOG(R0mn) &
                      + dble(ir - 1)*DLOG(R0mx/R0mn)/dble(nb - 1)
            R0(ir) = DEXP(phi(ir))
        end do
        dphi = phi(2) - phi(1)

        ! weights for quadrature using Simpson's rule
        do ir = 2, nb - 1
            ! Gaussian
            tmp = DEXP(-0.5d0*(phi(ir)/sd)**2)/DSQRT(2.d0*pi)/sd
            if (mod(ir, 2) == 0) then
                weight(ir) = tmp*4.d0*dphi/3.d0
            else
                weight(ir) = tmp*2.d0*dphi/3.d0
            end if
        end do

        tmp = DEXP(-0.5d0*(phi(1)/sd)**2)/DSQRT(2.d0*pi)/sd
        weight(1) = tmp*dphi/3.d0
        tmp = DEXP(-0.5d0*(phi(nb)/sd)**2)/DSQRT(2.d0*pi)/sd
        weight(nb) = tmp*dphi/3.d0

    end subroutine s_simpson

end module m_global_parameters