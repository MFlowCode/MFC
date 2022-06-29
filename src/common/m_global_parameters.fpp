!>
!! @file m_global_parameters.f90
!! @brief Contains module m_global_parameters
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

#:include 'case.fpp'

!> @brief This module contains all of the parameters characterizing the
!!              computational domain, simulation algorithm, initial condition
!!              and the stiffened equation of state.
module m_global_parameters

    ! Dependencies =============================================================
    use m_derived_types         ! Definitions of the derived types

    ! ==========================================================================

    implicit none

    ! Logistics ================================================================
    integer                    :: num_procs                  !< Number of processors
    integer,         parameter :: num_stcls_min = 5          !< Mininum # of stencils
    integer,         parameter :: path_len      = 400        !< Maximum path length
    integer,         parameter :: name_len      = 50         !< Maximum name length
    CHARACTER ,      parameter :: dflt_char     = ' '        !< Default string value
    real(kind(0d0)), parameter :: dflt_real     = -1d6       !< Default real value
    integer,         parameter :: dflt_int      = -100       !< Default integer value
    integer,         parameter :: fourier_rings = 5          !< Fourier filter ring limit
    REAL(KIND(0d0)), parameter :: sgm_eps       = 1d-16      !< Segmentation tolerance
    REAL(KIND(0d0)), parameter :: small_alf     = 1d-7       !< Small alf tolerance
    character(LEN=path_len)    :: case_dir      = "${CASE['logistics']['case_dir']}$"       !< Case folder location
    logical,         parameter :: run_time_info = .${CASE['logistics']['run_time_info']}$.  !< Run-time output flag
    logical,         parameter :: debug         = .${CASE['logistics']['debug']}$.          !< Debug mode print statements
    logical,         parameter :: old_grid      = .${CASE['logistics']['old_grid']}$.       !< Use existing grid data
    logical,         parameter :: old_ic        = .${CASE['logistics']['old_ic']}$.         !< Use existing IC data
    integer,         parameter :: t_step_old    = ${CASE['logistics']['t_step_old']}$       !< Existing IC/grid folder
!$acc declare create(small_alf, dflt_real, dflt_int, sgm_eps, fourier_rings)
    ! ==========================================================================

    ! Computational Domain Parameters ==========================================

    integer :: proc_rank !< Rank of the local processor

    integer, parameter :: m_glb = ${CASE['domain']['cells']['x']}$
    integer, parameter :: n_glb = ${CASE['domain']['cells']['y']}$
    integer, parameter :: p_glb = ${CASE['domain']['cells']['z']}$
    !! Global number of cells in each direction

    integer :: m_root = dflt_int
    integer :: m = m_glb
    integer :: n = n_glb
    integer :: p = p_glb
    !! Number of cells in the x-, y- and z-coordinate directions

    logical, parameter :: cyl_coord     = .${CASE['domain']['cyl_coord']}$.
    integer, parameter :: grid_geometry = ${CASE['autogen']['grid_geometry']}$!<
!$acc declare create(cyl_coord, grid_geometry)
    !! Cylindrical coordinates (either axisymmetric or full 3D)

    real(kind(0d0)), target, allocatable, dimension(:) :: x_cc, x_root_cc, y_cc, z_cc !<
    !! Locations of cell-centers (cc) in x-, y- and z-directions, respectively

    real(kind(0d0)), target, allocatable, dimension(:) :: x_cb, x_root_cb, y_cb, z_cb !<
    real(kind(0d0)), allocatable, dimension(:) :: coarse_x_cb, coarse_y_cb, coarse_z_cb
    !! Locations of cell-boundaries (cb) in x-, y- and z-directions, respectively

    real(kind(0d0)) :: dx_min, dy_min, dz_min !<
    !! Minimum cell-widths in the x-, y- and z-coordinate directions

    !> @name Cell-width distributions in the x-, y- and z-directions, respectively
    !> @{
    REAL(KIND(0d0)), TARGET, ALLOCATABLE, DIMENSION(:) :: dx, dy, dz
    !> @}

    REAL(KIND(0d0)) :: dt = ${CASE['domain']['time']['dt']}$ !< Size of the time-step
!$acc declare create(x_cb, y_cb, z_cb, x_cc, y_cc, z_cc, dx, dy, dz, dt, m, n, p)

    !> @name Starting time-step iteration, stopping time-step iteration and the number
    !! of time-step iterations between successive solution backups, respectively
    !> @{
    INTEGER, parameter :: t_step_start = ${CASE['domain']['time']['begin']}$
    INTEGER, parameter :: t_step_stop  = ${CASE['domain']['time']['end']}$
    INTEGER, parameter :: t_step_save  = ${CASE['domain']['time']['save']}$
    !> @}

    type(bounds_info) :: x_domain = ${CASE['autogen']['x_domain']}$
    type(bounds_info) :: y_domain = ${CASE['autogen']['y_domain']}$
    type(bounds_info) :: z_domain = ${CASE['autogen']['z_domain']}$
    !! Locations of the domain bounds in the x-, y- and z-coordinate directions

    logical, parameter :: stretch_x = .${CASE['domain']['stretch']['x']['stretch']}$.
    logical, parameter :: stretch_y = .${CASE['domain']['stretch']['y']['stretch']}$.
    logical, parameter :: stretch_z = .${CASE['domain']['stretch']['z']['stretch']}$.
    !<
    !! Grid stretching flags for the x-, y- and z-coordinate directions

    ! Parameters of the grid stretching function for the x-, y- and z-coordinate
    ! directions. The "a" parameters are a measure of the rate at which the grid
    ! is stretched while the remaining parameters are indicative of the location
    ! on the grid at which the stretching begins.
    real(kind(0d0)), parameter :: a_x     = ${CASE['domain']['stretch']['x']['rate']}$
    real(kind(0d0)), parameter :: a_y     = ${CASE['domain']['stretch']['y']['rate']}$
    real(kind(0d0)), parameter :: a_z     = ${CASE['domain']['stretch']['z']['rate']}$
    integer,         parameter :: loops_x = ${CASE['domain']['stretch']['x']['loops']}$
    integer,         parameter :: loops_y = ${CASE['domain']['stretch']['y']['loops']}$
    integer,         parameter :: loops_z = ${CASE['domain']['stretch']['z']['loops']}$
    real(kind(0d0))            :: x_a     = ${CASE['domain']['stretch']['x']['begin_neg']}$
    real(kind(0d0))            :: y_a     = ${CASE['domain']['stretch']['y']['begin_neg']}$
    real(kind(0d0))            :: z_a     = ${CASE['domain']['stretch']['z']['begin_neg']}$
    real(kind(0d0))            :: x_b     = ${CASE['domain']['stretch']['x']['begin_pos']}$
    real(kind(0d0))            :: y_b     = ${CASE['domain']['stretch']['y']['begin_pos']}$
    real(kind(0d0))            :: z_b     = ${CASE['domain']['stretch']['z']['begin_pos']}$

    ! ==========================================================================

    ! Note: We assume 1 + min(1, n) + min(1, p) == 1 + min(1, n_glb) + min(1, p_glb)

    ! Simulation Algorithm Parameters ==========================================
    INTEGER, parameter :: model_eqns       = ${CASE['algorithm']['model']}$             !< Multicomponent flow model
    INTEGER, parameter :: num_dims         = ${1 + min(1, CASE['domain']['cells']['y']) + min(1, CASE['domain']['cells']['z'])}$!< Number of spatial dimensions
    LOGICAL, parameter :: adv_alphan       = .${CASE['algorithm']['adv_alphan']}$.      !< Advection of the last volume fraction
    LOGICAL, parameter :: mpp_lim          = .${CASE['algorithm']['mpp_lim']}$.         !< Mixture physical parameters (MPP) limits
    INTEGER, parameter :: time_stepper     = ${CASE['algorithm']['time_stepper']}$      !< Time-stepper algorithm
    INTEGER, parameter :: weno_vars        = ${CASE['algorithm']['weno']['variables']}$ !< WENO-reconstructed state variables type
    INTEGER, parameter :: weno_order       = ${CASE['algorithm']['weno']['order']}$     !< Order of the WENO reconstruction
    INTEGER, parameter :: weno_polyn       = (weno_order - 1) / 2  !< Degree of the WENO polynomials (polyn)
    REAL(KIND(0d0)), parameter :: weno_eps = ${CASE['algorithm']['weno']['epsilon']}$      !< Binding for the WENO nonlinear weights
    LOGICAL, parameter :: char_decomp      = .${CASE['algorithm']['char_decomp']}$. !< Characteristic decomposition
    LOGICAL, parameter :: mapped_weno      = .${CASE['algorithm']['weno']['mapped']}$.    !< WENO with mapping of nonlinear weights
    LOGICAL, parameter :: mp_weno          = .${CASE['algorithm']['weno']['mp']}$.      !< Monotonicity preserving (MP) WENO
    LOGICAL, parameter :: weno_avg         = .${CASE['algorithm']['weno']['average']}$.  !< Average left/right cell-boundary states
    LOGICAL, parameter :: weno_Re_flux     = .${CASE['algorithm']['weno']['Re_flux']}$.  !< WENO reconstruct velocity gradients for viscous stress tensor
    logical, parameter :: weno_flat        = .${CASE['algorithm']['weno']['flat']}$.
    INTEGER, parameter :: riemann_solver   = ${CASE['algorithm']['riemann_solver']}$   !< Riemann solver algorithm
    INTEGER, parameter :: wave_speeds      = ${CASE['algorithm']['wave_speeds']}$      !< Wave speeds estimation method
    INTEGER, parameter :: avg_state        = ${CASE['algorithm']['avg_state']}$        !< Average state evaluation method
    LOGICAL, parameter :: commute_err      = .${CASE['algorithm']['commute_err']}$. !< Commutative error correction
    LOGICAL, parameter :: split_err        = .${CASE['algorithm']['split_err']}$.!< Dimensional splitting error correction
    LOGICAL, parameter :: alt_crv          = .${CASE['algorithm']['alt_crv']}$. !< Alternate curvature definition
    LOGICAL, parameter :: alt_soundspeed   = .${CASE['algorithm']['alt_soundspeed']}$. !< Alternate mixture sound speed
    LOGICAL, parameter :: regularization   = .${CASE['algorithm']['regularization']}$. !< Regularization terms of Tiwari (2013)
    REAL(KIND(0d0)), parameter :: reg_eps  = ${CASE['algorithm']['reg_eps']}$     !< User-defined interface thickness parameter for regularization terms
    LOGICAL, parameter :: null_weights     = .${CASE['algorithm']['null_weights']}$.!< Null undesired WENO weights
    LOGICAL, parameter :: mixture_err      = .${CASE['algorithm']['mixture_err']}$.!< Mixture properties correction
    LOGICAL, parameter :: tvd_riemann_flux = .${CASE['algorithm']['tvd_riemann_flux']}$.!< Apply TVD flux limiter to left and right states inside Riemann solver
    LOGICAL, parameter :: tvd_rhs_flux     = .${CASE['algorithm']['tvd_rhs_flux']}$.!< Apply TVD flux limiter to to intercell fluxes outside Riemann solver
    LOGICAL, parameter :: tvd_wave_speeds  = .${CASE['algorithm']['tvd_wave_speeds']}$. !< Use TVD wavespeeds when computing fluxes inside Riemann solver
    INTEGER, parameter :: flux_lim         = ${CASE['algorithm']['flux_lim']}$ !< Choice of flux limiter
    LOGICAL, parameter :: We_riemann_flux  = .${CASE['algorithm']['We_riemann_flux']}$. !< Account for capillary effects in the Riemann solver
    LOGICAL, parameter :: We_rhs_flux      = .${CASE['algorithm']['We_rhs_flux']}$. !< Account for capillary effects using the conservative formulation in RHS
    LOGICAL, parameter :: We_src           = .${CASE['algorithm']['We_src']}$. !< Account for capillary effects in non-conservative formulation in RHS
    LOGICAL, parameter :: We_wave_speeds   = .${CASE['algorithm']['We_wave_speeds']}$. !< Account for capillary effects when computing the contact wave speed
    LOGICAL, parameter :: lsq_deriv        = .${CASE['algorithm']['lsq_deriv']}$. !< Use linear least squares to calculate normals and curvatures
    LOGICAL, parameter :: hypoelasticity   = .${CASE['algorithm']['hypoelasticity']}$. !< Hypoelastic modeling
!$acc declare create(weno_polyn, mpp_lim, num_fluids, model_eqns, num_dims, mixture_err, alt_soundspeed, avg_state, mapped_weno, mp_weno, weno_eps)

    INTEGER         :: cpu_start, cpu_end, cpu_rate

    !> @name Boundary conditions (BC) in the x-, y- and z-directions, respectively
    !> @{
    TYPE(int_bounds_info), parameter :: bc_x_glb = int_bounds_info(${CASE['algorithm']['boundary']['x']['begin']}$, ${CASE['algorithm']['boundary']['x']['end']}$)
    TYPE(int_bounds_info), parameter :: bc_y_glb = int_bounds_info(${CASE['algorithm']['boundary']['y']['begin']}$, ${CASE['algorithm']['boundary']['y']['end']}$)
    TYPE(int_bounds_info), parameter :: bc_z_glb = int_bounds_info(${CASE['algorithm']['boundary']['z']['begin']}$, ${CASE['algorithm']['boundary']['z']['end']}$)

    TYPE(int_bounds_info) :: bc_x = bc_x_glb
    TYPE(int_bounds_info) :: bc_y = bc_y_glb
    TYPE(int_bounds_info) :: bc_z = bc_z_glb

    !> @}

    !> @name Annotations of the structure of the state and flux vectors in terms of the
    !! size and the configuration of the system of equations to which they belong
    !> @{
    INTEGER               :: sys_size              !< Number of unknowns in system of eqns.
    TYPE(int_bounds_info) :: cont_idx              !< Indexes of first & last continuity eqns.
    TYPE(int_bounds_info) :: mom_idx               !< Indexes of first & last momentum eqns.
    INTEGER               :: E_idx                 !< Index of energy equation
    TYPE(int_bounds_info) :: adv_idx               !< Indexes of first & last advection eqns.
    TYPE(int_bounds_info) :: internalEnergies_idx  !< Indexes of first & last internal energy eqns.
    TYPE(bub_bounds_info) :: bub_idx               !< Indexes of first & last bubble variable eqns.
    INTEGER               :: alf_idx               !< Index of void fraction
    INTEGER               :: gamma_idx             !< Index of specific heat ratio func. eqn.
    INTEGER               :: pi_inf_idx            !< Index of liquid stiffness func. eqn.
    TYPE(int_bounds_info) :: stress_idx            !< Indexes of first and last shear stress eqns.
    !> @}
!$acc declare create(bub_idx)


    !> @name The number of fluids, along with their identifying indexes, respectively,
    !! for which viscous effects, e.g. the shear and/or the volume Reynolds (Re)
    !! numbers, will be non-negligible.
    !> @{
    INTEGER,              DIMENSION(2)   :: Re_size
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: Re_idx
    !> @}
!$acc declare create(Re_size, Re_idx)


    !> @name The number of idiosyncratic material interfaces, along with the indexes of
    !! the fluids comprising them, respectively, for which the effects of surface
    !! tension, e.g. the Weber (We) numbers, will be non-negligible.
    !> @{
    INTEGER                              :: We_size
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: We_idx
    !> @}

    !> @name The number of fluids, along with their identifying indexes, respectively,
    !! for which the physical and the geometric curvatures (crv) of the material
    !! interfaces will be calculated.
    !> @{
    INTEGER                            :: crv_size
    INTEGER, ALLOCATABLE, DIMENSION(:) :: crv_idx
    !> @}

    !> @name The coordinate direction indexes and flags (flg), respectively, for which
    !! the configurations will be determined with respect to a working direction
    !! and that will be used to isolate the contributions, in that direction, in
    !! the dimensionally split system of equations.
    !> @{
    INTEGER        , DIMENSION(3) :: dir_idx
    REAL(KIND(0d0)), DIMENSION(3) :: dir_flg
    !> @}
!$acc declare create(dir_idx, dir_flg)

    !! extra coordinate direction index used if hypoelasticity = true
    INTEGER        , DIMENSION(3) :: dir_idx_tau

    REAL(KIND(0d0)) :: wa_flg !<
    !! The WENO average (WA) flag regulates whether the calculation of any cell-
    !! average spatial derivatives is carried out in each cell by utilizing the
    !! arithmetic mean of the left and right, WENO-reconstructed, cell-boundary
    !! values or simply, the unaltered left and right, WENO-reconstructed, cell-
    !! boundary values.

    INTEGER :: buff_size !<
    !! The number of cells that are necessary to be able to store enough boundary
    !! conditions data to march the solution in the physical computational domain
    !! to the next time-step.


    integer :: startx, starty, startz

!$acc declare create(sys_size, buff_size, startx, starty, startz, E_idx, gamma_idx, pi_inf_idx, alf_idx)


    ! END: Simulation Algorithm Parameters =====================================

    ! Formatted Database File(s) Structure Parameters ==========================
    
    integer, parameter :: format       = ${CASE['database']['format']}$ !< Format of the database file(s)
    logical, parameter :: coarsen_silo = .${CASE['database']['coarsen_silo']}$.
    logical, parameter :: parallel_io  = .${CASE['database']['parallel_io']}$. !< Format of the data files
    integer            :: precision    = ${CASE['database']['precision']}$ !< Precision of output files

    !> @name Size of the ghost zone layer in the x-, y- and z-coordinate directions.
    !! The definition of the ghost zone layers is only necessary when using the
    !! Silo database file format in multidimensions. These zones provide VisIt
    !! with the subdomain connectivity information that it requires in order to
    !! produce smooth plots.
    !> @{
    type(int_bounds_info) :: offset_x, offset_y, offset_z
    !> @}

    ! Perturb density of surrounding air so as to break symmetry of grid
    logical, parameter :: perturb_flow       = .${CASE['algorithm']['perturb_flow']}$.
    integer, parameter :: perturb_flow_fluid = ${CASE['algorithm']['perturb_flow_fluid']}$   !< Fluid to be perturbed with perturb_flow flag
    logical, parameter :: perturb_sph        = .${CASE['algorithm']['perturb_sph']}$.
    integer, parameter :: perturb_sph_fluid  = ${CASE['algorithm']['perturb_sph_fluid']}$  !< Fluid to be perturbed with perturb_sph flag
    real(kind(0d0)), dimension(num_fluids), parameter :: fluid_rho = ${CASE['algorithm']['fluid_rho']}$ 

    integer, allocatable, dimension(:) :: proc_coords !<
    !! Processor coordinates in MPI_CART_COMM

    integer, allocatable, dimension(:) :: start_idx !<
    !! Starting cell-center index of local processor in global grid

    type(mpi_io_var), public :: MPI_IO_DATA

    character(LEN=name_len) :: mpiiofs
    integer :: mpi_info_int !<
    !! MPI info for parallel IO with Lustre file systems

    integer, private :: ierr
    ! ==========================================================================

    ! Initial Condition Parameters =============================================
    type(ic_patch_parameters), dimension(num_patches) :: patch_icpp = ${CASE['autogen']['patch_icpp']}$

    !<
    !! Database of the initial condition patch parameters (icpp) for each of the
    !! patches employed in the configuration of the initial condition.
    ! ==========================================================================

    ! Fluids Physical Parameters ===============================================
    type(physical_parameters), dimension(num_fluids_alloc), parameter :: fluid_pp = ${CASE['autogen']['fluid_pp']}$

    !<
    !! Database of the physical parameters of each of the fluids that is present
    !! in the flow. These include the stiffened gas equation of state parameters,
    !! the Reynolds numbers and the Weber numbers.

    ! ==========================================================================

    INTEGER, parameter :: fd_order = ${CASE['database']['fd_order']}$ !<
    !! The order of the finite-difference (fd) approximations of the first-order
    !! derivatives that need to be evaluated when the CoM or flow probe data
    !! files are to be written at each time step

    INTEGER, parameter :: fd_number = ${max(1, CASE['database']['fd_order'] // 2)}$ !<
    !! The finite-difference number is given by MAX(1, fd_order/2). Essentially,
    !! it is a measure of the half-size of the finite-difference stencil for the
    !! selected order of accuracy.

    LOGICAL, DIMENSION(num_fluids_alloc) :: com_wrt, cb_wrt
    LOGICAL, parameter :: probe_wrt     = .${CASE['database']['write']['probe']}$.
    LOGICAL, parameter :: integral_wrt  = .${CASE['database']['write']['integral']}$.
    INTEGER, parameter :: num_probes    = ${len(CASE['database']['probes'])}$
    INTEGER, parameter :: num_integrals = ${len(CASE['database']['integrals'])}$
    TYPE(probe_parameters),    DIMENSION(num_probes),    parameter :: probe    = ${CASE['autogen']['probe']}$
    TYPE(integral_parameters), DIMENSION(num_integrals), parameter :: integral = ${CASE['autogen']['integral']}$
    REAL(KIND(0d0)), DIMENSION(5) :: threshold_mf = dflt_real
    INTEGER, DIMENSION(5) :: moment_order = (/ dflt_int, dflt_int, dflt_int, dflt_int, dflt_int /)

    !> @name The list of all possible flow variables that may be written to a database
    !! file. It includes partial densities, density, momentum, velocity, energy,
    !! pressure, volume fraction(s), specific heat ratio function, specific heat
    !! ratio, liquid stiffness function, liquid stiffness, primitive variables,
    !! conservative variables, speed of sound, the vorticity,
    !! and the numerical Schlieren function.
    !> @{
    logical, dimension(num_fluids_alloc) :: alpha_rho_wrt  = .${CASE['database']['write']['alpha_rho']}$.
    logical, parameter                   :: rho_wrt        = .${CASE['database']['write']['rho']}$.
    logical, dimension(3)                :: mom_wrt        = .${CASE['database']['write']['mom']}$.
    logical, dimension(3)                :: vel_wrt        = .${CASE['database']['write']['velocity']}$.
    logical, dimension(3)                :: flux_wrt       = .${CASE['database']['write']['flux']}$.
    logical, parameter                   :: E_wrt          = .${CASE['database']['write']['E']}$.
    logical, parameter                   :: pres_wrt       = .${CASE['database']['write']['pressure']}$.
    logical, dimension(num_fluids_alloc) :: alpha_wrt      = .${CASE['database']['write']['alpha']}$.
    logical, parameter                   :: gamma_wrt      = .${CASE['database']['write']['gamma']}$.
    logical, parameter                   :: heat_ratio_wrt = .${CASE['database']['write']['heat_ratio']}$.
    logical, parameter                   :: pi_inf_wrt     = .${CASE['database']['write']['pi_inf']}$.
    logical, parameter                   :: pres_inf_wrt   = .${CASE['database']['write']['pressure_inf']}$.
    logical, parameter                   :: prim_vars_wrt  = .${CASE['database']['write']['prim_vars']}$.
    logical, parameter                   :: cons_vars_wrt  = .${CASE['database']['write']['cons_vars']}$.
    logical, parameter                   :: c_wrt          = .${CASE['database']['write']['c']}$.
    logical, dimension(3)                :: omega_wrt      = .${CASE['database']['write']['omega']}$.
    logical, parameter                   :: schlieren_wrt  = .${CASE['database']['write']['schlieren']}$.
    !> @}


    real(kind(0d0)), dimension(num_fluids_alloc) :: schlieren_alpha = dflt_real   !<
    !! Amplitude coefficients of the numerical Schlieren function that are used
    !! to adjust the intensity of numerical Schlieren renderings for individual
    !! fluids. This enables waves and interfaces of varying strenghts and in all
    !! of the fluids to be made simulatenously visible on a single plot.

    real(kind(0d0)) :: rhoref = ${CASE['algorithm']['rhoref']}$
    real(kind(0d0)) :: pref   = ${CASE['algorithm']['pref']}$
    !< Reference parameters for Tait EOS
!$acc declare create(rhoref, pref)


    !> @name Bubble modeling
    !> @{
    INTEGER,         parameter :: nb     = ${CASE['bubbles']['number']}$  !< Number of eq. bubble sizes
    REAL(KIND(0d0)), parameter :: R0ref  = ${CASE['bubbles']['R0ref']}$  !< Reference bubble size
    REAL(KIND(0d0)), parameter :: Ca     = ${CASE['bubbles']['cavitation']}$       !< Cavitation number
    REAL(KIND(0d0)), parameter :: Web    = ${CASE['bubbles']['weber']}$     !< Weber number
    REAL(KIND(0d0)), parameter :: Re_inv = ${CASE['bubbles']['Re_inv']}$  !< Inverse Reynolds number
    REAL(KIND(0d0)), DIMENSION(:), ALLOCATABLE :: weight !< Simpson quadrature weights
    REAL(KIND(0d0)), DIMENSION(:), ALLOCATABLE :: R0     !< Bubble sizes
    REAL(KIND(0d0)), DIMENSION(:), ALLOCATABLE :: V0     !< Bubble velocities
    LOGICAL,         parameter :: bubbles      = .${CASE['bubbles']['bubbles']}$. !< Bubbles on/off
    LOGICAL,         parameter :: polytropic   = .${CASE['bubbles']['polytropic']}$.  !< Polytropic  switch
    LOGICAL,         parameter :: polydisperse = .${CASE['bubbles']['polydisperse']}$. !< Polydisperse bubbles

    INTEGER, parameter :: bubble_model = ${CASE['bubbles']['model']}$ !< Gilmore or Keller--Miksis bubble model
    INTEGER, parameter :: thermal      = ${CASE['bubbles']['thermal']}$ !< Thermal behavior. 1 = adiabatic, 2 = isotherm, 3 = transfer
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:,:) :: ptil  !< Pressure modification
    REAL(KIND(0d0)), parameter :: poly_sigma = ${CASE['bubbles']['poly_sigma']}$ !< log normal sigma for polydisperse PDF

    LOGICAL,         parameter :: qbmm      = .${CASE['bubbles']['qbmm']}$. !< Quadrature moment method
    INTEGER,         parameter :: nmom      = 6 !< Number of carried moments per R0 location
    INTEGER,         parameter :: nnode     = ${CASE['bubbles']['nnode']}$ !< Number of QBMM nodes
    INTEGER,         parameter :: nmomsp    = 4 !< Number of moments required by ensemble-averaging
    INTEGER,         parameter :: nmomtot   = nmom*nb !< Total number of carried moments moments/transport equations
    integer,         parameter :: dist_type = ${CASE['bubbles']['distribution']}$ !1 = binormal, 2 = lognormal-normal
    INTEGER,         parameter :: R0_type   = ${CASE['bubbles']['R0_type']}$ !< R0 distribution type
    real(kind(0d0)), parameter :: sigR      = ${CASE['bubbles']['sigR']}$
    real(kind(0d0)), parameter :: sigV      = ${CASE['bubbles']['sigV']}$
    real(kind(0d0)), parameter :: rhoRV     = ${CASE['bubbles']['rhoRV']}$

    TYPE(scalar_field), ALLOCATABLE, DIMENSION(:)     :: mom_sp
    TYPE(scalar_field), ALLOCATABLE, DIMENSION(:,:,:) :: mom_3d
    !> @}
!$acc declare create(nb, R0ref, Ca, Web, Re_inv, weight, R0, V0, bubbles, polytropic, polydisperse, qbmm, nmom, nnode, nmomsp, nmomtot, R0_type, ptil, bubble_model, thermal, poly_sigma)
!$acc declare create(mom_sp, mom_3d)


    !> @name Physical bubble parameters (see Ando 2010, Preston 2007)
    !> @{
    REAL(KIND(0d0)) :: R_n    = dflt_real, R_v    = dflt_real
    REAL(KIND(0d0)) :: phi_vn = dflt_real, phi_nv = dflt_real
    REAL(KIND(0d0)) :: Pe_c   = dflt_real, Tw     = dflt_real
    REAL(KIND(0d0)) :: pv, M_n, M_v
    REAL(KIND(0d0)), DIMENSION(:), ALLOCATABLE :: k_n, k_v, pb0, mass_n0, mass_v0, Pe_T
    REAL(KIND(0d0)), DIMENSION(:), ALLOCATABLE :: Re_trans_T, Re_trans_c, Im_trans_T, Im_trans_c, omegaN
    REAL(KIND(0d0)) :: mul0, ss, gamma_v, mu_v, G
    REAL(KIND(0d0)) :: gamma_m, gamma_n, mu_n
    REAL(KIND(0d0)) :: gam
    !> @}
!$acc declare create(R_n, R_v, phi_vn, phi_nv, Pe_c, Tw, pv, M_n, M_v, k_n, k_v, pb0, mass_n0, mass_v0, Pe_T, Re_trans_T, Re_trans_c, Im_trans_T, Im_trans_c, omegaN , mul0, ss, gamma_v, mu_v, gamma_m, gamma_n, mu_n, gam)


    !> @name Acoustic monopole parameters
    !> @{
    !< Monopole switch
    LOGICAL, parameter :: monopole = .${CASE['acoustic']['monopole']}$.
    !< Number of monopoles
    INTEGER, parameter :: num_mono = ${len(CASE['acoustic']['monopoles'])}$
    !< Monopole parameters
    TYPE(mono_parameters), DIMENSION(num_mono), parameter :: mono = ${CASE['autogen']['mono']}$
    !> @}
!$acc declare create(monopole, mono, num_mono)


    integer, allocatable, dimension(:, :, :) :: logic_grid

    REAL(KIND(0d0))            :: mytime       !< Current simulation time
    REAL(KIND(0d0))            :: dt0          !< Initial time step size
    REAL(KIND(0d0)), parameter :: finaltime = t_step_stop*${CASE['domain']['time']['dt']}$ !< Final simulation time

    logical, parameter :: riemann_flat = .${CASE['algorithm']['riemann_flat']}$.
    logical, parameter :: cu_mpi       = .${CASE['logistics']['cu_mpi']}$.
    logical, parameter :: cu_tensor    = .${CASE['logistics']['cu_tensor']}$.

    ! Mathematical and Physical Constants ======================================
    real(kind(0d0)), parameter :: pi = 3.141592653589793d0
!$acc declare create(pi)
    ! ==========================================================================

contains
    !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_global_parameters_module() ! -------------------
!$acc routine seq

        integer :: tmp_idx !< Temporary indexes storage
        integer :: i, j !< Generic loop iterators
        integer :: k !< Generic counter
        integer :: fac
        integer :: i1, i2, i3

        type(int_bounds_info) :: ix, iy, iz

!$acc update device(weno_polyn)
!$acc update device(nb)


        ! Setting m_root equal to m in the case of a 1D serial simulation
        if (num_procs == 1 .and. n == 0) then
            m_root = m
        end if

        ! Initializing the number of fluids for which viscous effects will
        ! be non-negligible, the number of distinctive material interfaces
        ! for which surface tension will be important and also, the number
        ! of fluids for which the physical and geometric curvatures of the
        ! interfaces will be computed
        Re_size = 0

        ! Gamma/Pi_inf Model ===============================================
        if (model_eqns == 1) then

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

            ! Volume Fraction Model ============================================
        else

            ! Annotating structure of the state and flux vectors belonging
            ! to the system of equations defined by the selected number of
            ! spatial dimensions and the volume fraction model
            if (model_eqns == 2) then
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
                    alf_idx = 0
                end if

                if (bubbles) then
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
                    ! print*, 'alf idx', alf_idx
                    ! print*, 'bub -idx beg end', bub_idx%beg, bub_idx%end

                    allocate (weight(nb), R0(nb), V0(nb))
                    allocate (bub_idx%rs(nb), bub_idx%vs(nb))
                    allocate (bub_idx%ps(nb), bub_idx%ms(nb))

                    if (num_fluids == 1) then
                        gam = 1.d0/fluid_pp(num_fluids + 1)%gamma + 1.d0
                    else
                        gam = 1.d0/fluid_pp(num_fluids)%gamma + 1.d0
                    end if

                    if (qbmm) then
                        allocate (bub_idx%moms(nb, nmom))
                        allocate ( bub_idx%fullmom(nb,0:nmom,0:nmom) )

                        do i = 1, nb
                            do j = 1, nmom
                                bub_idx%moms(i, j) = bub_idx%beg + (j - 1) + (i - 1)*nmom
                            end do
                            bub_idx%rs(i) = bub_idx%moms(i, 2)
                            bub_idx%vs(i) = bub_idx%moms(i, 3)

                            bub_idx%fullmom(i,0,0) = bub_idx%moms(i,1)
                            bub_idx%fullmom(i,1,0) = bub_idx%moms(i,2)
                            bub_idx%fullmom(i,0,1) = bub_idx%moms(i,3)
                            bub_idx%fullmom(i,2,0) = bub_idx%moms(i,4)
                            bub_idx%fullmom(i,1,1) = bub_idx%moms(i,5)
                            bub_idx%fullmom(i,0,2) = bub_idx%moms(i,6)

                            bub_idx%rs(i) = bub_idx%fullmom(i,1,0)
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
                            print*, 'Invalid R0 type - abort'
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


            else if (model_eqns == 3) then
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
            else if (model_eqns == 4) then
                cont_idx%beg = 1 ! one continuity equation
                cont_idx%end = 1 !num_fluids
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
                        if (polytropic) then
                            fac = 2
                        else
                            fac = 4
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
                            print*, 'Invalid R0 type - abort'
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

            ! Determining the number of fluids for which the shear and the
            ! volume Reynolds numbers, e.g. viscous effects, are important
            do i = 1, num_fluids
                if (fluid_pp(i)%Re(1) > 0) Re_size(1) = Re_size(1) + 1
                if (fluid_pp(i)%Re(2) > 0) Re_size(2) = Re_size(2) + 1
            end do

            ! Bookkeeping the indexes of any viscous fluids and any pairs of
            ! fluids whose interface will support effects of surface tension
            if (any(Re_size > 0)) then

                allocate (Re_idx(1:2, 1:maxval(Re_size)))

                k = 0
                do i = 1, num_fluids
                    if (fluid_pp(i)%Re(1) > 0) then
                        k = k + 1; Re_idx(1, k) = i
                    end if
                end do

                k = 0
                do i = 1, num_fluids
                    if (fluid_pp(i)%Re(2) > 0) then
                        k = k + 1; Re_idx(2, k) = i
                    end if
                end do

            end if

        end if
        ! END: Volume Fraction Model =======================================

        allocate (MPI_IO_DATA%view(1:sys_size))
        allocate (MPI_IO_DATA%var(1:sys_size))

        do i = 1, sys_size
            allocate (MPI_IO_DATA%var(i)%sf(0:m, 0:n, 0:p))
            MPI_IO_DATA%var(i)%sf => null()
        end do

!$acc update device(Re_size, Re_idx)
        ! Determining the number of cells that are needed in order to store
        ! sufficient boundary conditions data as to iterate the solution in
        ! the physical computational domain from one time-step iteration to
        ! the next one
        if (any(Re_size > 0)) then
            buff_size = 2*weno_polyn + 2
        else
            buff_size = weno_polyn + 2
        end if





        ! Configuring Coordinate Direction Indexes =========================
        if (bubbles) then
            ix%beg = -buff_size; iy%beg = 0; iz%beg = 0
            if (n > 0) iy%beg = -buff_size; if (p > 0) iz%beg = -buff_size
            ix%end = m - ix%beg; iy%end = n - iy%beg; iz%end = p - iz%beg
            allocate (ptil(ix%beg:ix%end, &
                           iy%beg:iy%end, &
                           iz%beg:iz%end))
        end if

        if (probe_wrt) then
            buff_size = buff_size + fd_number
        end if


        startx = -buff_size
        starty = 0
        startz = 0
        if(n > 0) then
            starty = -buff_size
        end if
        if(p > 0) then
            startz = -buff_size
        end if

!$acc update device(startx, starty, startz)

        ! Allocating grid variables for the x-, y- and z-directions
        allocate (x_cb(-1 - buff_size:m + buff_size))
        allocate (x_cc(-buff_size:m + buff_size))
        allocate (dx(-buff_size:m + buff_size))

        if (n > 0 .and. n /= dflt_int) then
            allocate (y_cb(-1 - buff_size:n + buff_size))
            allocate (y_cc(-buff_size:n + buff_size))
            allocate (dy(-buff_size:n + buff_size))
            if (p > 0 .and. p /= dflt_int) then
                allocate (z_cb(-1 - buff_size:p + buff_size))
                allocate (z_cc(-buff_size:p + buff_size))
                allocate (dz(-buff_size:p + buff_size))
            end if

            ! Allocating the grid variables, only used for the 1D simulations,
            ! and containing the defragmented computational domain grid data
        else
            allocate (x_root_cb(-1:m_root))
            allocate (x_root_cc(0:m_root))
        end if

        if (coarsen_silo) then
            allocate (coarse_x_cb(-1 - offset_x%beg:(m/2) + offset_x%end))

            if (n > 0 .and. n /= dflt_int) then
                allocate (coarse_y_cb(-1 - offset_y%beg:(n/2) + offset_y%end))
                
                if (p > 0 .and. p /= dflt_int) then
                    allocate (coarse_z_cb(-1 - offset_z%beg:(p/2) + offset_z%end))
                end if
            end if
        end if

        allocate (logic_grid(0:m, 0:n, 0:p))

    end subroutine s_initialize_global_parameters_module ! -----------------

    !> Initializes non-polydisperse bubble modeling
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

        real(kind(0.d0)), parameter :: k_poly = 1.d0 !<
            !! polytropic index used to compute isothermal natural frequency

        real(kind(0.d0)), parameter :: Ru = 8314.d0 !<
            !! universal gas constant

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
        if (thermal == 2) gamma_m = 1.d0

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

    !>  Computes transfer coefficient for non-polydisperse bubble modeling (Preston 2007)
        !!  @param omega Frequency
        !!  @param peclet Peclet number
        !!  @param Re_trans Real part of transfer coefficient
        !!  @param Im_trans Imaginary part of transfer coefficient
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

    !> Initializes parallel infrastructure
    subroutine s_initialize_parallel_io() ! --------------------------------

        allocate (proc_coords(1:num_dims))

        if (parallel_io .neqv. .true.) return

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

    end subroutine s_initialize_parallel_io ! ------------------------------

    !> Module deallocation and/or disassociation procedures
    subroutine s_finalize_global_parameters_module() ! ---------------------
!$acc routine seq

        integer :: i

        ! Deallocating the variables bookkeeping the indexes of any viscous
        ! fluids and any pairs of fluids whose interfaces supported effects
        ! of surface tension
        if (any(Re_size > 0)) deallocate (Re_idx)

        ! Deallocating grid variables for the x-, y- and z-directions
        deallocate(x_cb, x_cc, dx)

        if (n > 0 .and. n /= dflt_int) then
            deallocate(y_cb, y_cc, dy)

            if (p > 0 .and. p /= dflt_int) then
                deallocate(z_cb, z_cc, dz)
            end if
        else
            deallocate(x_root_cb, x_root_cc)
        end if

        if (coarsen_silo) then
            deallocate (coarse_x_cb)

            if (n > 0 .and. n /= dflt_int) then
                deallocate (coarse_y_cb)

                if (p > 0 .and. p /= dflt_int) then
                    deallocate (coarse_z_cb)
                end if
            end if
        end if

        deallocate (proc_coords)
        if (parallel_io) then
            deallocate (start_idx)
            do i = 1, sys_size
                MPI_IO_DATA%var(i)%sf => null()
            end do

            deallocate (MPI_IO_DATA%var)
            deallocate (MPI_IO_DATA%view)
        end if

    end subroutine s_finalize_global_parameters_module ! -------------------

    !>  Computes the bubble number density n from the conservative variables
        !!  \f$ n = \sqrt{ \frac{4 \pi}{3} } \frac{ nR^3}{\alpha} \f$
        !! @param vftmp is the void fraction
        !! @param nRtmp is the bubble number  density times the bubble radii
        !! @param ntmp is the output number bubble density
    subroutine s_comp_n_from_cons(vftmp, nRtmp, ntmp)
!$acc routine seq

        real(kind(0d0)), intent(IN) :: vftmp
        real(kind(0d0)), dimension(:), intent(IN) :: nRtmp
        real(kind(0d0)), intent(OUT) :: ntmp
        real(kind(0d0)) :: nR3
        integer :: i

        nR3 = 0d0
        do i = 1, nb
            nR3 = nR3 + weight(i)*(nRtmp(i)**3d0)
        end do


        !if (nR3 < 0d0) then
            ! DO i = 1,nb
            ! IF (nRtmp(i) < small_alf) THEN
            ! nRtmp(i) = small_alf
            ! END IF
            ! END DO
            ! nR3 = 1.d-12
            !print *, vftmp, nR3, nRtmp(:)
         !   stop 'nR3 is negative'
        !end if
        !if (vftmp < 0d0) then
            ! vftmp = small_alf
            ! ntmp = DSQRT( (4.d0*pi/3.d0)*nR3/1.d-12 )
            !print *, vftmp, nR3, nRtmp(:)
         !   stop 'vf negative'
        !end if

        ntmp = DSQRT((4.d0*pi/3.d0)*nR3/vftmp)

    end subroutine s_comp_n_from_cons

    !> Computes the bubble number density n from the primitive variables
        !!  \f$ n = \sqrt{ \frac{3}{4 \pi} } \frac{ \alpha }{ R^3} \f$
        !! @param vftmp is the void fraction
        !! @param Rtmp is the  bubble radii
        !! @param ntmp is the output number bubble density
    subroutine s_comp_n_from_prim(vftmp, Rtmp, ntmp)
!$acc routine seq

        real(kind(0.d0)), intent(IN) :: vftmp
        real(kind(0.d0)), dimension(:), intent(IN) :: Rtmp
        real(kind(0.d0)), intent(OUT) :: ntmp
        real(kind(0.d0)) :: R3
        integer :: i

        R3 = 0d0
        do i = 1, nb
            R3 = R3 + weight(i)*(Rtmp(i)**3d0)
        end do

        IF ( R3 < 0d0 ) THEN
            !PRINT*, vftmp, R3, Rtmp(:)
            STOP 'R3 is negative'
        END IF
        IF (vftmp < 0d0) THEN
            !PRINT*, vftmp, R3, Rtmp(:)
            STOP 'vf negative'
        END IF

        ntmp = (3.d0/(4.d0*pi))*vftmp/R3

    end subroutine s_comp_n_from_prim

    !> Computes the quadrature for polydisperse bubble populations
        !! @param func is the bubble dynamic variables for each bin
        !! @param mom is the computed moment
    subroutine s_quad(func, mom)
!$acc routine seq

        real(kind(0.d0)), dimension(:), intent(IN) :: func
        real(kind(0.d0)), intent(OUT) :: mom
        integer :: i

        mom = 0d0
        do i = 1, nb
            mom = mom + weight(i)*func(i)
        end do


    end subroutine s_quad

    !> Computes the Simpson weights for quadrature
    subroutine s_simpson
!$acc routine seq

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
