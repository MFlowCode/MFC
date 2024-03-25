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
    integer :: t_step_old, t_step_start           !< Existing IC/grid folder
    ! ==========================================================================

    ! Computational Domain Parameters ==========================================

    integer :: proc_rank !< Rank of the local processor

    integer :: m
    integer :: n
    integer :: p !<
    !! Number of cells in the x-, y- and z-coordinate directions

    integer(8) :: nGlobal ! Global number of cells in the domain

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
    logical :: relax           !< activate phase change
    integer :: relax_model     !< Relax Model
    real(kind(0d0)) :: palpha_eps     !< trigger parameter for the p relaxation procedure, phase change model
    real(kind(0d0)) :: ptgalpha_eps   !< trigger parameter for the pTg relaxation procedure, phase change model
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
    integer :: n_idx                      !< Index of number density
    type(int_bounds_info) :: adv_idx                    !< Indexes of first & last advection eqns.
    type(int_bounds_info) :: internalEnergies_idx       !< Indexes of first & last internal energy eqns.
    type(bub_bounds_info) :: bub_idx                    !< Indexes of first & last bubble variable eqns.
    integer :: gamma_idx                  !< Index of specific heat ratio func. eqn.
    integer :: pi_inf_idx                 !< Index of liquid stiffness func. eqn.
    type(int_bounds_info) :: stress_idx                 !< Indexes of elastic shear stress eqns.

    type(int_bounds_info) :: bc_x, bc_y, bc_z !<
    !! Boundary conditions in the x-, y- and z-coordinate directions

    logical :: parallel_io !< Format of the data files
    logical :: file_per_process !< type of data output
    integer :: precision !< Precision of output files

    logical :: vel_profile !< Set hyperbolic tangent streamwise velocity profile
    logical :: instability_wave !< Superimpose instability waves to surrounding fluid flow

    real(kind(0d0)) :: pi_fac !< Factor for artificial pi_inf

    ! Perturb density of surrounding air so as to break symmetry of grid
    logical :: perturb_flow
    integer :: perturb_flow_fluid   !< Fluid to be perturbed with perturb_flow flag
    real(kind(0d0)) :: perturb_flow_mag   !< Magnitude of perturbation with perturb_flow flag
    logical :: perturb_sph
    integer :: perturb_sph_fluid    !< Fluid to be perturbed with perturb_sph flag
    real(kind(0d0)), dimension(num_fluids_max) :: fluid_rho

    integer, allocatable, dimension(:) :: proc_coords !<
    !! Processor coordinates in MPI_CART_COMM

    integer, allocatable, dimension(:) :: start_idx !<
    !! Starting cell-center index of local processor in global grid

#ifdef MFC_MPI

    type(mpi_io_var), public :: MPI_IO_DATA
    type(mpi_io_ib_var), public :: MPI_IO_IB_DATA
    type(mpi_io_airfoil_ib_var), public :: MPI_IO_airfoil_IB_DATA

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
    logical :: adv_n !< Solve the number density equation and compute alpha from number density
    !> @}

    !> @name Immersed Boundaries
    !> @{
    logical :: ib           !< Turn immersed boundaries on
    integer :: num_ibs      !< Number of immersed boundaries
    integer :: Np

    type(ib_patch_parameters), dimension(num_patches_max) :: patch_ib

    type(probe_parameters), allocatable, dimension(:) :: airfoil_grid_u, airfoil_grid_l
    !! Database of the immersed boundary patch parameters for each of the
    !! patches employed in the configuration of the initial condition. Note that
    !! the maximum allowable number of patches, num_patches_max, may be changed
    !! in the module m_derived_types.f90.
    ! ==========================================================================

    !> @}

    !> @name Non-polytropic bubble gas compression
    !> @{
    logical :: polytropic
    logical :: polydisperse
    integer :: thermal  !1 = adiabatic, 2 = isotherm, 3 = transfer
    real(kind(0d0)) :: R_n, R_v, phi_vn, phi_nv, Pe_c, Tw, pv, M_n, M_v
    real(kind(0d0)), dimension(:), allocatable :: k_n, k_v, pb0, mass_n0, mass_v0, Pe_T
    real(kind(0d0)), dimension(:), allocatable :: Re_trans_T, Re_trans_c, Im_trans_T, Im_trans_c, omegaN
    real(kind(0d0)) :: mul0, ss, gamma_v, mu_v
    real(kind(0d0)) :: gamma_m, gamma_n, mu_n
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

    type(pres_field) :: pb
    type(pres_field) :: mv

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
        t_step_start = dflt_int

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
        relax = .false.
        relax_model = dflt_int
        palpha_eps = dflt_real
        ptgalpha_eps = dflt_real
        num_fluids = dflt_int
        adv_alphan = .false.
        weno_order = dflt_int

        hypoelasticity = .false.

        bc_x%beg = dflt_int; bc_x%end = dflt_int
        bc_y%beg = dflt_int; bc_y%end = dflt_int
        bc_z%beg = dflt_int; bc_z%end = dflt_int

        #:for DIM in ['x', 'y', 'z']
            #:for DIR in [1, 2, 3]
                bc_${DIM}$%vb${DIR}$ = 0d0
                bc_${DIM}$%ve${DIR}$ = 0d0
            #:endfor
        #:endfor

        parallel_io = .false.
        file_per_process = .false.
        precision = 2
        vel_profile = .false.
        instability_wave = .false.
        perturb_flow = .false.
        perturb_flow_fluid = dflt_int
        perturb_flow_mag = dflt_real
        perturb_sph = .false.
        perturb_sph_fluid = dflt_int
        fluid_rho = dflt_real

        ! Initial condition parameters
        num_patches = dflt_int

        do i = 1, num_patches_max
            patch_icpp(i)%geometry = dflt_int
            patch_icpp(i)%model%scale(:) = 1d0
            patch_icpp(i)%model%translate(:) = 0d0
            patch_icpp(i)%model%filepath(:) = ' '
            patch_icpp(i)%model%spc = 10
            patch_icpp(i)%model%threshold = 0.9d0
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
            patch_icpp(i)%cv = 0d0
            patch_icpp(i)%qv = 0d0
            patch_icpp(i)%qvp = 0d0
            patch_icpp(i)%tau_e = 0d0
            !should get all of r0's and v0's
            patch_icpp(i)%r0 = dflt_real
            patch_icpp(i)%v0 = dflt_real

            patch_icpp(i)%p0 = dflt_real
            patch_icpp(i)%m0 = dflt_real

            patch_icpp(i)%hcid = dflt_int
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

        adv_n = .false.

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

        pi_fac = 1d0

        ! Immersed Boundaries
        ib = .false.
        num_ibs = dflt_int

        do i = 1, num_patches_max
            patch_ib(i)%geometry = dflt_int
            patch_ib(i)%x_centroid = dflt_real
            patch_ib(i)%y_centroid = dflt_real
            patch_ib(i)%z_centroid = dflt_real
            patch_ib(i)%length_x = dflt_real
            patch_ib(i)%length_y = dflt_real
            patch_ib(i)%length_z = dflt_real
            patch_ib(i)%radius = dflt_real
            patch_ib(i)%theta = dflt_real
            patch_ib(i)%c = dflt_real
            patch_ib(i)%t = dflt_real
            patch_ib(i)%m = dflt_real
            patch_ib(i)%p = dflt_real
            patch_ib(i)%slip = .false.
        end do

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
            fluid_pp(i)%cv = 0d0
            fluid_pp(i)%qv = 0d0
            fluid_pp(i)%qvp = 0d0
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

                if (adv_n) then
                    n_idx = bub_idx%end + 1
                    sys_size = n_idx
                end if

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
                    V0(:) = 1d0
                    !R0 and weight initialized in s_simpson
                else
                    stop 'Invalid value of nb'
                end if

                !Initialize pref,rhoref for polytropic qbmm (done in s_initialize_nonpoly for non-polytropic)
                if (.not. qbmm) then
                    if (polytropic) then
                        rhoref = 1.d0
                        pref = 1.d0
                    end if
                end if

                !Initialize pb0,pv,pref,rhoref for polytropic qbmm (done in s_initialize_nonpoly for non-polytropic)
                if (qbmm) then
                    if (polytropic) then
                        allocate (pb0(nb))
                        if (Web == dflt_real) then
                            pb0 = pref
                            pb0 = pb0/pref
                            pref = 1d0
                        end if
                        rhoref = 1d0
                    end if
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
                    V0(:) = 1d0
                else
                    stop 'Invalid value of nb'
                end if

                if (polytropic) then
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

        if (qbmm .and. .not. polytropic) then
            allocate (MPI_IO_DATA%view(1:sys_size + 2*nb*4))
            allocate (MPI_IO_DATA%var(1:sys_size + 2*nb*4))
        else
            allocate (MPI_IO_DATA%view(1:sys_size))
            allocate (MPI_IO_DATA%var(1:sys_size))
        end if

        do i = 1, sys_size
            allocate (MPI_IO_DATA%var(i)%sf(0:m, 0:n, 0:p))
            MPI_IO_DATA%var(i)%sf => null()
        end do
        if (qbmm .and. .not. polytropic) then
            do i = sys_size + 1, sys_size + 2*nb*4
                allocate (MPI_IO_DATA%var(i)%sf(0:m, 0:n, 0:p))
                MPI_IO_DATA%var(i)%sf => null()
            end do
        end if

        if (ib) allocate (MPI_IO_IB_DATA%var%sf(0:m, 0:n, 0:p))

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

        if (ib) deallocate (MPI_IO_IB_DATA%var%sf)

#endif

    end subroutine s_finalize_global_parameters_module ! ----------------------

end module m_global_parameters
