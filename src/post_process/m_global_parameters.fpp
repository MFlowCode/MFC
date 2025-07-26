!>
!! @file m_global_parameters.f90
!! @brief Contains module m_global_parameters

#:include 'case.fpp'

!> @brief This module contains all of the parameters characterizing the
!!      computational domain, simulation algorithm, stiffened equation of
!!      state and finally, the formatted database file(s) structure.
module m_global_parameters

#ifdef MFC_MPI
    use mpi                     !< Message passing interface (MPI) module
#endif

    use m_derived_types         !< Definitions of the derived types

    use m_helper_basic          !< Functions to compare floating point numbers

    use m_thermochem, only: num_species, species_names

    implicit none

    !> @name Logistics
    !> @{
    integer :: num_procs            !< Number of processors
    character(LEN=path_len) :: case_dir             !< Case folder location
    !> @}

    ! Computational Domain Parameters

    integer :: proc_rank !< Rank of the local processor

    !> @name Number of cells in the x-, y- and z-coordinate directions
    !> @{
    integer :: m, m_root
    integer :: n
    integer :: p
    !> @}

    !> @name Max and min number of cells in a direction of each combination of x-,y-, and z-
    type(cell_num_bounds) :: cells_bounds

    integer(8) :: nGlobal ! Total number of cells in global domain

    !> @name Cylindrical coordinates (either axisymmetric or full 3D)
    !> @{
    logical :: cyl_coord
    integer :: grid_geometry
    !> @}

    !> @name Global number of cells in each direction
    !> @{
    integer :: m_glb, n_glb, p_glb
    !> @}

    integer :: num_dims !< Number of spatial dimensions
    integer :: num_vels !< Number of velocity components (different from num_dims for mhd)

    !> @name Cell-boundary locations in the x-, y- and z-coordinate directions
    !> @{
    real(wp), allocatable, dimension(:) :: x_cb, x_root_cb, y_cb, z_cb
    real(wp), allocatable, dimension(:) :: x_cb_s, y_cb_s, z_cb_s
    !> @}

    !> @name Cell-center locations in the x-, y- and z-coordinate directions
    !> @{
    real(wp), allocatable, dimension(:) :: x_cc, x_root_cc, y_cc, z_cc
    real(sp), allocatable, dimension(:) :: x_root_cc_s, x_cc_s
    !> @}

    !> Cell-width distributions in the x-, y- and z-coordinate directions
    !> @{
    real(wp), allocatable, dimension(:) :: dx, dy, dz
    !> @}

    integer :: buff_size !<
    !! Number of cells in buffer region. For the variables which feature a buffer
    !! region, this region is used to store information outside the computational
    !! domain based on the boundary conditions.

    integer :: t_step_start  !< First time-step directory
    integer :: t_step_stop   !< Last time-step directory
    integer :: t_step_save   !< Interval between consecutive time-step directory

    !> @name IO options for adaptive time-stepping
    !> @{
    logical :: cfl_adap_dt, cfl_const_dt, cfl_dt
    real(wp) :: t_save
    real(wp) :: t_stop
    real(wp) :: cfl_target
    integer :: n_save
    integer :: n_start
    !> @}

    ! NOTE: The variables m_root, x_root_cb and x_root_cc contain the grid data
    ! of the defragmented computational domain. They are only used in 1D. For
    ! serial simulations, they are equal to m, x_cb and x_cc, respectively.

    !> @name Simulation Algorithm Parameters
    !> @{
    integer :: model_eqns      !< Multicomponent flow model
    integer :: num_fluids      !< Number of different fluids present in the flow
    logical :: relax           !< phase change
    integer :: relax_model     !< Phase change relaxation model
    logical :: mpp_lim         !< Maximum volume fraction limiter
    integer :: sys_size        !< Number of unknowns in the system of equations
    integer :: recon_type      !< Which type of reconstruction to use
    integer :: weno_order      !< Order of accuracy for the WENO reconstruction
    integer :: muscl_order     !< Order of accuracy for the MUSCL reconstruction
    logical :: mixture_err     !< Mixture error limiter
    logical :: alt_soundspeed  !< Alternate sound speed
    logical :: mhd             !< Magnetohydrodynamics
    logical :: relativity      !< Relativity for RMHD
    logical :: hypoelasticity  !< Turn hypoelasticity on
    logical :: hyperelasticity !< Turn hyperelasticity on
    logical :: elasticity      !< elasticity modeling, true for hyper or hypo
    integer :: b_size          !< Number of components in the b tensor
    integer :: tensor_size     !< Number of components in the nonsymmetric tensor
    logical :: cont_damage     !< Continuum damage modeling
    logical :: igr             !< enable IGR
    integer :: igr_order       !< IGR reconstruction order
    logical, parameter :: chemistry = .${chemistry}$. !< Chemistry modeling
    !> @}

    integer :: avg_state       !< Average state evaluation method

    !> @name Annotations of the structure, i.e. the organization, of the state vectors
    !> @{
    type(int_bounds_info) :: cont_idx              !< Indexes of first & last continuity eqns.
    type(int_bounds_info) :: mom_idx               !< Indexes of first & last momentum eqns.
    integer :: E_idx                               !< Index of energy equation
    integer :: n_idx                               !< Index of number density
    integer :: beta_idx                            !< Index of lagrange bubbles beta
    type(int_bounds_info) :: adv_idx               !< Indexes of first & last advection eqns.
    type(int_bounds_info) :: internalEnergies_idx  !< Indexes of first & last internal energy eqns.
    type(bub_bounds_info) :: bub_idx               !< Indexes of first & last bubble variable eqns.
    integer :: gamma_idx                           !< Index of specific heat ratio func. eqn.
    integer :: alf_idx                             !< Index of specific heat ratio func. eqn.
    integer :: pi_inf_idx                          !< Index of liquid stiffness func. eqn.
    type(int_bounds_info) :: B_idx                 !< Indexes of first and last magnetic field eqns.
    type(int_bounds_info) :: stress_idx            !< Indices of elastic stresses
    type(int_bounds_info) :: xi_idx                !< Indexes of first and last reference map eqns.
    integer :: c_idx                               !< Index of color function
    type(int_bounds_info) :: species_idx           !< Indexes of first & last concentration eqns.
    integer :: damage_idx                          !< Index of damage state variable (D) for continuum damage model
    !> @}

    ! Cell Indices for the (local) interior points (O-m, O-n, 0-p).
    ! Stands for "InDices With BUFFer".
    type(int_bounds_info) :: idwint(1:3)

    ! Cell Indices for the entire (local) domain. In simulation, this includes
    ! the buffer region. idwbuff and idwint are the same otherwise.
    ! Stands for "InDices With BUFFer".
    type(int_bounds_info) :: idwbuff(1:3)

    integer :: num_bc_patches
    logical :: bc_io
    !> @name Boundary conditions in the x-, y- and z-coordinate directions
    !> @{
    type(int_bounds_info) :: bc_x, bc_y, bc_z
    !> @}

    integer :: shear_num !! Number of shear stress components
    integer, dimension(3) :: shear_indices !<
    !! Indices of the stress components that represent shear stress
    integer :: shear_BC_flip_num !<
    !! Number of shear stress components to reflect for boundary conditions
    integer, dimension(3, 2) :: shear_BC_flip_indices !<
    !! Indices of shear stress components to reflect for boundary conditions.
    !! Size: (1:3, 1:shear_BC_flip_num) for (x/y/z, [indices])

    logical :: parallel_io    !< Format of the data files
    logical :: sim_data
    logical :: file_per_process !< output format

    integer, allocatable, dimension(:) :: proc_coords !<
    !! Processor coordinates in MPI_CART_COMM

    integer, allocatable, dimension(:) :: start_idx !<
    !! Starting cell-center index of local processor in global grid

    integer :: num_ibs  !< Number of immersed boundaries

#ifdef MFC_MPI

    type(mpi_io_var), public :: MPI_IO_DATA
    type(mpi_io_ib_var), public :: MPI_IO_IB_DATA
    type(mpi_io_levelset_var), public :: MPI_IO_levelset_DATA
    type(mpi_io_levelset_norm_var), public :: MPI_IO_levelsetnorm_DATA
    real(wp), allocatable, dimension(:, :), public :: MPI_IO_DATA_lg_bubbles

#endif

    !> @name MPI info for parallel IO with Lustre file systems
    !> @{
    character(LEN=name_len) :: mpiiofs
    integer :: mpi_info_int
    !> @}

    type(physical_parameters), dimension(num_fluids_max) :: fluid_pp !<
    !! Database of the physical parameters of each of the fluids that is present
    !! in the flow. These include the stiffened gas equation of state parameters,
    !! the Reynolds numbers and the Weber numbers.

    real(wp), allocatable, dimension(:) :: adv !< Advection variables

    ! Formatted Database File(s) Structure Parameters

    integer :: format !< Format of the database file(s)

    integer :: precision !< Floating point precision of the database file(s)

    logical :: output_partial_domain !< Specify portion of domain to output for post-processing

    type(bounds_info) :: x_output, y_output, z_output !< Portion of domain to output for post-processing
    type(int_bounds_info) :: x_output_idx, y_output_idx, z_output_idx !< Indices of domain to output for post-processing

    !> @name Size of the ghost zone layer in the x-, y- and z-coordinate directions.
    !! The definition of the ghost zone layers is only necessary when using the
    !! Silo database file format in multidimensions. These zones provide VisIt
    !! with the subdomain connectivity information that it requires in order to
    !! produce smooth plots.
    !> @{
    type(int_bounds_info) :: offset_x, offset_y, offset_z
    !> @}

    !> @name The list of all possible flow variables that may be written to a database
    !! file. It includes partial densities, density, momentum, velocity, energy,
    !! pressure, volume fraction(s), specific heat ratio function, specific heat
    !! ratio, liquid stiffness function, liquid stiffness, primitive variables,
    !! conservative variables, speed of sound, the vorticity,
    !! and the numerical Schlieren function.
    !> @{
    logical, dimension(num_fluids_max) :: alpha_rho_wrt
    logical :: rho_wrt
    logical, dimension(3) :: mom_wrt
    logical, dimension(3) :: vel_wrt
    integer :: flux_lim
    logical, dimension(3) :: flux_wrt
    logical :: E_wrt
    logical :: pres_wrt
    logical, dimension(num_fluids_max) :: alpha_wrt
    logical :: gamma_wrt
    logical :: heat_ratio_wrt
    logical :: pi_inf_wrt
    logical :: pres_inf_wrt
    logical :: prim_vars_wrt
    logical :: cons_vars_wrt
    logical :: c_wrt
    logical, dimension(3) :: omega_wrt
    logical :: qm_wrt
    logical :: schlieren_wrt
    logical :: cf_wrt
    logical :: ib
    logical :: chem_wrt_Y(1:num_species)
    logical :: chem_wrt_T
    !> @}

    real(wp), dimension(num_fluids_max) :: schlieren_alpha    !<
    !! Amplitude coefficients of the numerical Schlieren function that are used
    !! to adjust the intensity of numerical Schlieren renderings for individual
    !! fluids. This enables waves and interfaces of varying strengths and in all
    !! of the fluids to be made simultaneously visible on a single plot.

    integer :: fd_order !<
    !! The order of the finite-difference (fd) approximations of the first-order
    !! derivatives that need to be evaluated when vorticity and/or the numerical
    !! Schlieren function are to be outputted to the formatted database file(s).

    integer :: fd_number !<
    !! The finite-difference number is given by MAX(1, fd_order/2). Essentially,
    !! it is a measure of the half-size of the finite-difference stencil for the
    !! selected order of accuracy.

    !> @name Reference parameters for Tait EOS
    !> @{
    real(wp) :: rhoref, pref
    !> @}

    !> @name Bubble modeling variables and parameters
    !> @{
    integer :: nb
    real(wp) :: R0ref
    real(wp) :: Ca, Web, Re_inv
    real(wp), dimension(:), allocatable :: weight, R0
    logical :: bubbles_euler
    logical :: qbmm
    logical :: polytropic
    logical :: polydisperse
    logical :: adv_n
    integer :: thermal  !< 1 = adiabatic, 2 = isotherm, 3 = transfer
    real(wp) :: R_n, R_v, phi_vn, phi_nv, Pe_c, Tw, G, pv, M_n, M_v
    real(wp), dimension(:), allocatable :: k_n, k_v, pb0, mass_n0, mass_v0, Pe_T
    real(wp), dimension(:), allocatable :: Re_trans_T, Re_trans_c, Im_trans_T, Im_trans_c, omegaN
    real(wp) :: mul0, ss, gamma_v, mu_v
    real(wp) :: gamma_m, gamma_n, mu_n
    real(wp) :: poly_sigma
    real(wp) :: sigR
    integer :: nmom
    !> @}

    !> @name surface tension coefficient
    !> @{

    real(wp) :: sigma
    logical :: surface_tension
    !> #}

    !> @name Index variables used for m_variables_conversion
    !> @{
    integer :: momxb, momxe
    integer :: advxb, advxe
    integer :: contxb, contxe
    integer :: intxb, intxe
    integer :: bubxb, bubxe
    integer :: strxb, strxe
    integer :: xibeg, xiend
    integer :: chemxb, chemxe
    !> @}

    !> @name Lagrangian bubbles
    !> @{
    logical :: bubbles_lagrange
    !> @}

    real(wp) :: Bx0 !< Constant magnetic field in the x-direction (1D)

contains

    !> Assigns default values to user inputs prior to reading
        !!      them in. This allows for an easier consistency check of
        !!      these parameters once they are read from the input file.
    impure subroutine s_assign_default_values_to_user_inputs

        integer :: i !< Generic loop iterator

        ! Logistics
        case_dir = '.'

        ! Computational domain parameters
        m = dflt_int; n = 0; p = 0
        call s_update_cell_bounds(cells_bounds, m, n, p)

        m_root = dflt_int
        cyl_coord = .false.

        t_step_start = dflt_int
        t_step_stop = dflt_int
        t_step_save = dflt_int

        cfl_adap_dt = .false.
        cfl_const_dt = .false.
        cfl_dt = .false.
        cfl_target = dflt_real
        t_save = dflt_real
        n_start = dflt_int
        t_stop = dflt_real

        ! Simulation algorithm parameters
        model_eqns = dflt_int
        num_fluids = dflt_int
        recon_type = WENO_TYPE
        weno_order = dflt_int
        muscl_order = dflt_int
        mixture_err = .false.
        alt_soundspeed = .false.
        relax = .false.
        relax_model = dflt_int

        mhd = .false.
        relativity = .false.

        hypoelasticity = .false.
        hyperelasticity = .false.
        elasticity = .false.
        b_size = dflt_int
        tensor_size = dflt_int
        cont_damage = .false.
        igr = .false.

        bc_x%beg = dflt_int; bc_x%end = dflt_int
        bc_y%beg = dflt_int; bc_y%end = dflt_int
        bc_z%beg = dflt_int; bc_z%end = dflt_int
        bc_io = .false.
        num_bc_patches = dflt_int

        #:for DIM in ['x', 'y', 'z']
            #:for DIR in [1, 2, 3]
                bc_${DIM}$%vb${DIR}$ = 0._wp
                bc_${DIM}$%ve${DIR}$ = 0._wp
            #:endfor
        #:endfor

        ! Fluids physical parameters
        do i = 1, num_fluids_max
            fluid_pp(i)%gamma = dflt_real
            fluid_pp(i)%pi_inf = dflt_real
            fluid_pp(i)%cv = 0._wp
            fluid_pp(i)%qv = 0._wp
            fluid_pp(i)%qvp = 0._wp
            fluid_pp(i)%G = dflt_real
        end do

        ! Formatted database file(s) structure parameters
        format = dflt_int

        precision = dflt_int

        alpha_rho_wrt = .false.
        rho_wrt = .false.
        mom_wrt = .false.
        vel_wrt = .false.
        chem_wrt_Y = .false.
        chem_wrt_T = .false.
        flux_lim = dflt_int
        flux_wrt = .false.
        parallel_io = .false.
        file_per_process = .false.
        E_wrt = .false.
        pres_wrt = .false.
        alpha_wrt = .false.
        gamma_wrt = .false.
        heat_ratio_wrt = .false.
        pi_inf_wrt = .false.
        pres_inf_wrt = .false.
        prim_vars_wrt = .false.
        cons_vars_wrt = .false.
        c_wrt = .false.
        omega_wrt = .false.
        qm_wrt = .false.
        schlieren_wrt = .false.
        sim_data = .false.
        cf_wrt = .false.
        ib = .false.

        schlieren_alpha = dflt_real

        fd_order = dflt_int
        avg_state = dflt_int

        ! Tait EOS
        rhoref = dflt_real
        pref = dflt_real

        ! Bubble modeling
        bubbles_euler = .false.
        qbmm = .false.
        R0ref = dflt_real
        nb = dflt_int
        polydisperse = .false.
        poly_sigma = dflt_real
        sigR = dflt_real
        sigma = dflt_real
        surface_tension = .false.
        adv_n = .false.

        ! Lagrangian bubbles modeling
        bubbles_lagrange = .false.

        ! IBM
        num_ibs = dflt_int

        ! Output partial domain
        output_partial_domain = .false.
        x_output%beg = dflt_real
        x_output%end = dflt_real
        y_output%beg = dflt_real
        y_output%end = dflt_real
        z_output%beg = dflt_real
        z_output%end = dflt_real

        ! MHD
        Bx0 = dflt_real

    end subroutine s_assign_default_values_to_user_inputs

    !>  Computation of parameters, allocation procedures, and/or
        !!      any other tasks needed to properly setup the module
    impure subroutine s_initialize_global_parameters_module

        integer :: i, j, fac

        ! Setting m_root equal to m in the case of a 1D serial simulation
        if (n == 0) m_root = m_glb

        ! Gamma/Pi_inf Model
        if (model_eqns == 1) then

            ! Setting number of fluids
            num_fluids = 1

            ! Annotating structure of the state and flux vectors belonging
            ! to the system of equations defined by the selected number of
            ! spatial dimensions and the gamma/pi_inf model
            cont_idx%beg = 1
            cont_idx%end = cont_idx%beg
            mom_idx%beg = cont_idx%end + 1
            mom_idx%end = cont_idx%end + num_vels
            E_idx = mom_idx%end + 1
            adv_idx%beg = E_idx + 1
            adv_idx%end = adv_idx%beg + 1
            gamma_idx = adv_idx%beg
            pi_inf_idx = adv_idx%end
            sys_size = adv_idx%end

            ! Volume Fraction Model (5-equation model)
        else if (model_eqns == 2) then

            ! Annotating structure of the state and flux vectors belonging
            ! to the system of equations defined by the selected number of
            ! spatial dimensions and the volume fraction model
            cont_idx%beg = 1
            cont_idx%end = num_fluids
            mom_idx%beg = cont_idx%end + 1
            mom_idx%end = cont_idx%end + num_vels
            E_idx = mom_idx%end + 1
            adv_idx%beg = E_idx + 1
            if (igr) then
                if (num_fluids == 1) then
                    adv_idx%end = adv_idx%beg
                else
                    adv_idx%end = E_idx + num_fluids - 1
                end if
            else
                adv_idx%end = E_idx + num_fluids
            end if

            sys_size = adv_idx%end

            if (bubbles_euler) then
                alf_idx = adv_idx%end
            else
                alf_idx = 1
            end if

            if (qbmm) then
                nmom = 6
            end if

            if (bubbles_euler) then

                bub_idx%beg = sys_size + 1
                if (qbmm) then
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

                allocate (bub_idx%rs(nb), bub_idx%vs(nb))
                allocate (bub_idx%ps(nb), bub_idx%ms(nb))
                allocate (weight(nb), R0(nb))

                if (qbmm) then
                    allocate (bub_idx%moms(nb, nmom))
                    do i = 1, nb
                        do j = 1, nmom
                            bub_idx%moms(i, j) = bub_idx%beg + (j - 1) + (i - 1)*nmom
                        end do
                        bub_idx%rs(i) = bub_idx%moms(i, 2)
                        bub_idx%vs(i) = bub_idx%moms(i, 3)
                    end do
                else
                    do i = 1, nb
                        if (polytropic .neqv. .true.) then
                            fac = 4
                        else
                            fac = 2
                        end if

                        bub_idx%rs(i) = bub_idx%beg + (i - 1)*fac
                        bub_idx%vs(i) = bub_idx%rs(i) + 1

                        if (polytropic .neqv. .true.) then
                            bub_idx%ps(i) = bub_idx%vs(i) + 1
                            bub_idx%ms(i) = bub_idx%ps(i) + 1
                        end if
                    end do
                end if

                if (nb == 1) then
                    weight(:) = 1._wp
                    R0(:) = 1._wp
                else if (nb < 1) then
                    stop 'Invalid value of nb'
                end if

                if (polytropic .neqv. .true.) then
                    !call s_initialize_nonpoly
                else
                    rhoref = 1._wp
                    pref = 1._wp
                end if

            end if

            if (bubbles_lagrange) then
                beta_idx = sys_size + 1
                sys_size = beta_idx
            end if

            if (mhd) then
                B_idx%beg = sys_size + 1
                if (n == 0) then
                    B_idx%end = sys_size + 2 ! 1D: By, Bz
                else
                    B_idx%end = sys_size + 3 ! 2D/3D: Bx, By, Bz
                end if
                sys_size = B_idx%end
            end if

            ! Volume Fraction Model (6-equation model)
        else if (model_eqns == 3) then

            ! Annotating structure of the state and flux vectors belonging
            ! to the system of equations defined by the selected number of
            ! spatial dimensions and the volume fraction model
            cont_idx%beg = 1
            cont_idx%end = num_fluids
            mom_idx%beg = cont_idx%end + 1
            mom_idx%end = cont_idx%end + num_vels
            E_idx = mom_idx%end + 1
            adv_idx%beg = E_idx + 1
            adv_idx%end = E_idx + num_fluids
            internalEnergies_idx%beg = adv_idx%end + 1
            internalEnergies_idx%end = adv_idx%end + num_fluids
            sys_size = internalEnergies_idx%end
            alf_idx = 1 ! dummy, cannot actually have a void fraction

        else if (model_eqns == 4) then
            cont_idx%beg = 1 ! one continuity equation
            cont_idx%end = 1 !num_fluids
            mom_idx%beg = cont_idx%end + 1 ! one momentum equation in each
            mom_idx%end = cont_idx%end + num_vels
            E_idx = mom_idx%end + 1 ! one energy equation
            adv_idx%beg = E_idx + 1
            adv_idx%end = adv_idx%beg !one volume advection equation
            alf_idx = adv_idx%end
            sys_size = alf_idx !adv_idx%end

            if (bubbles_euler) then
                bub_idx%beg = sys_size + 1
                bub_idx%end = sys_size + 2*nb
                if (polytropic .neqv. .true.) then
                    bub_idx%end = sys_size + 4*nb
                end if
                sys_size = bub_idx%end

                allocate (bub_idx%rs(nb), bub_idx%vs(nb))
                allocate (bub_idx%ps(nb), bub_idx%ms(nb))
                allocate (weight(nb), R0(nb))

                do i = 1, nb
                    if (polytropic .neqv. .true.) then
                        fac = 4
                    else
                        fac = 2
                    end if

                    bub_idx%rs(i) = bub_idx%beg + (i - 1)*fac
                    bub_idx%vs(i) = bub_idx%rs(i) + 1

                    if (polytropic .neqv. .true.) then
                        bub_idx%ps(i) = bub_idx%vs(i) + 1
                        bub_idx%ms(i) = bub_idx%ps(i) + 1
                    end if
                end do

                if (nb == 1) then
                    weight(:) = 1._wp
                    R0(:) = 1._wp
                else if (nb < 1) then
                    stop 'Invalid value of nb'
                end if

                if (polytropic) then
                    rhoref = 1._wp
                    pref = 1._wp
                end if
            end if
        end if

        if (model_eqns == 2 .or. model_eqns == 3) then

            if (hypoelasticity .or. hyperelasticity) then
                elasticity = .true.
                stress_idx%beg = sys_size + 1
                stress_idx%end = sys_size + (num_dims*(num_dims + 1))/2
                if (cyl_coord) stress_idx%end = stress_idx%end + 1
                ! number of stresses is 1 in 1D, 3 in 2D, 4 in 2D-Axisym, 6 in 3D
                sys_size = stress_idx%end

                ! shear stress index is 2 for 2D and 2,4,5 for 3D
                if (num_dims == 1) then
                    shear_num = 0
                else if (num_dims == 2) then
                    shear_num = 1
                    shear_indices(1) = stress_idx%beg - 1 + 2
                    shear_BC_flip_num = 1
                    shear_BC_flip_indices(1:2, 1) = shear_indices(1)
                    ! Both x-dir and y-dir: flip tau_xy only
                else if (num_dims == 3) then
                    shear_num = 3
                    shear_indices(1:3) = stress_idx%beg - 1 + (/2, 4, 5/)
                    shear_BC_flip_num = 2
                    shear_BC_flip_indices(1, 1:2) = shear_indices((/1, 2/))
                    shear_BC_flip_indices(2, 1:2) = shear_indices((/1, 3/))
                    shear_BC_flip_indices(3, 1:2) = shear_indices((/2, 3/))
                    ! x-dir: flip tau_xy and tau_xz
                    ! y-dir: flip tau_xy and tau_yz
                    ! z-dir: flip tau_xz and tau_yz
                end if
            end if

            if (hyperelasticity) then
                xi_idx%beg = sys_size + 1
                xi_idx%end = sys_size + num_dims
                ! adding three more equations for the \xi field and the elastic energy
                sys_size = xi_idx%end + 1
                ! number of entries in the symmetric btensor plus the jacobian
                b_size = (num_dims*(num_dims + 1))/2 + 1
                tensor_size = num_dims**2 + 1
            end if

            if (surface_tension) then
                c_idx = sys_size + 1
                sys_size = c_idx
            end if

            if (cont_damage) then
                damage_idx = sys_size + 1
                sys_size = damage_idx
            else
                damage_idx = dflt_int
            end if

        end if

        if (chemistry) then
            species_idx%beg = sys_size + 1
            species_idx%end = sys_size + num_species
            sys_size = species_idx%end
        else
            species_idx%beg = 1
            species_idx%end = 1
        end if

        if (output_partial_domain) then
            x_output_idx%beg = 0
            x_output_idx%end = 0
            y_output_idx%beg = 0
            y_output_idx%end = 0
            z_output_idx%beg = 0
            z_output_idx%end = 0
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
        xibeg = xi_idx%beg
        xiend = xi_idx%end
        chemxb = species_idx%beg
        chemxe = species_idx%end

#ifdef MFC_MPI
        allocate (MPI_IO_DATA%view(1:sys_size))
        allocate (MPI_IO_DATA%var(1:sys_size))
        do i = 1, sys_size
            allocate (MPI_IO_DATA%var(i)%sf(0:m, 0:n, 0:p))
            MPI_IO_DATA%var(i)%sf => null()
        end do

        if (ib) allocate (MPI_IO_IB_DATA%var%sf(0:m, 0:n, 0:p))
#endif

        ! Size of the ghost zone layer is non-zero only when post-processing
        ! the raw simulation data of a parallel multidimensional computation
        ! in the Silo-HDF5 format. If this is the case, one must also verify
        ! whether the raw simulation data is 2D or 3D. In the 2D case, size
        ! of the z-coordinate direction ghost zone layer must be zeroed out.
        if (num_procs == 1 .or. format /= 1 .or. n == 0) then

            offset_x%beg = 0
            offset_x%end = 0
            offset_y%beg = 0
            offset_y%end = 0
            offset_z%beg = 0
            offset_z%end = 0

        elseif (p == 0) then

            offset_z%beg = 0
            offset_z%end = 0

        end if

        ! Determining the finite-difference number and the buffer size. Note
        ! that the size of the buffer is unrelated to the order of the WENO
        ! scheme. Rather, it is directly dependent on maximum size of ghost
        ! zone layers and possibly the order of the finite difference scheme
        ! used for the computation of vorticity and/or numerical Schlieren
        ! function.
        buff_size = max(offset_x%beg, offset_x%end, offset_y%beg, &
                        offset_y%end, offset_z%beg, offset_z%end)

        if (any(omega_wrt) .or. schlieren_wrt .or. qm_wrt) then
            fd_number = max(1, fd_order/2)
            buff_size = buff_size + fd_number
        end if

        ! Configuring Coordinate Direction Indexes
        idwint(1)%beg = 0; idwint(2)%beg = 0; idwint(3)%beg = 0
        idwint(1)%end = m; idwint(2)%end = n; idwint(3)%end = p

        idwbuff(1)%beg = -buff_size
        if (num_dims > 1) then; idwbuff(2)%beg = -buff_size; else; idwbuff(2)%beg = 0; end if
        if (num_dims > 2) then; idwbuff(3)%beg = -buff_size; else; idwbuff(3)%beg = 0; end if

        idwbuff(1)%end = idwint(1)%end - idwbuff(1)%beg
        idwbuff(2)%end = idwint(2)%end - idwbuff(2)%beg
        idwbuff(3)%end = idwint(3)%end - idwbuff(3)%beg

        ! Allocating single precision grid variables if needed
        if (precision == 1) then
            allocate (x_cb_s(-1 - offset_x%beg:m + offset_x%end))
            if (n > 0) then
                allocate (y_cb_s(-1 - offset_y%beg:n + offset_y%end))
                if (p > 0) then
                    allocate (z_cb_s(-1 - offset_z%beg:p + offset_z%end))
                end if
            end if
        else
            allocate (x_cc_s(-buff_size:m + buff_size))
        end if

        ! Allocating the grid variables in the x-coordinate direction
        allocate (x_cb(-1 - offset_x%beg:m + offset_x%end))
        allocate (x_cc(-buff_size:m + buff_size))
        allocate (dx(-buff_size:m + buff_size))

        ! Allocating grid variables in the y- and z-coordinate directions
        if (n > 0) then

            allocate (y_cb(-1 - offset_y%beg:n + offset_y%end))
            allocate (y_cc(-buff_size:n + buff_size))
            allocate (dy(-buff_size:n + buff_size))

            if (p > 0) then
                allocate (z_cb(-1 - offset_z%beg:p + offset_z%end))
                allocate (z_cc(-buff_size:p + buff_size))
                allocate (dz(-buff_size:p + buff_size))
            end if

            ! Allocating the grid variables, only used for the 1D simulations,
            ! and containing the defragmented computational domain grid data
        else

            allocate (x_root_cb(-1:m_root))
            allocate (x_root_cc(0:m_root))

            if (precision == 1) then
                allocate (x_root_cc_s(0:m_root))
            end if

        end if

        allocate (adv(num_fluids))

        if (cyl_coord .neqv. .true.) then ! Cartesian grid
            grid_geometry = 1
        elseif (cyl_coord .and. p == 0) then ! Axisymmetric cylindrical grid
            grid_geometry = 2
        else ! Fully 3D cylindrical grid
            grid_geometry = 3
        end if

    end subroutine s_initialize_global_parameters_module

    !> Subroutine to initialize parallel infrastructure
    impure subroutine s_initialize_parallel_io

#ifdef MFC_MPI
        integer :: ierr !< Generic flag used to identify and report MPI errors
#endif

        num_dims = 1 + min(1, n) + min(1, p)

        if (mhd) then
            num_vels = 3
        else
            num_vels = num_dims
        end if

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

    end subroutine s_initialize_parallel_io

    !> Deallocation procedures for the module
    impure subroutine s_finalize_global_parameters_module

        integer :: i

        ! Deallocating the grid variables for the x-coordinate direction
        deallocate (x_cc, x_cb, dx)

        ! Deallocating grid variables for the y- and z-coordinate directions
        if (n > 0) then
            deallocate (y_cc, y_cb, dy)
            if (p > 0) then
                deallocate (z_cc, z_cb, dz)
            end if
        else
            ! Deallocating the grid variables, only used for the 1D simulations,
            ! and containing the defragmented computational domain grid data
            deallocate (x_root_cb, x_root_cc)
        end if

        deallocate (proc_coords)

        deallocate (adv)

#ifdef MFC_MPI

        if (parallel_io) then
            deallocate (start_idx)
            do i = 1, sys_size
                MPI_IO_DATA%var(i)%sf => null()
            end do

            deallocate (MPI_IO_DATA%var)
            deallocate (MPI_IO_DATA%view)
        end if

        if (ib) MPI_IO_IB_DATA%var%sf => null()
#endif

    end subroutine s_finalize_global_parameters_module

end module m_global_parameters
