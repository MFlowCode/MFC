!>
!! @file m_global_parameters.f90
!! @brief Contains module m_global_parameters

!> @brief This module contains all of the parameters characterizing the
!!      computational domain, simulation algorithm, stiffened equation of
!!      state and finally, the formatted database file(s) structure.
module m_global_parameters

    ! Dependencies =============================================================
#ifdef MFC_MPI
    use mpi                     !< Message passing interface (MPI) module
#endif

    use m_derived_types         !< Definitions of the derived types
    ! ==========================================================================

    implicit none

    !> @name Logistics
    !> @{
    integer :: num_procs            !< Number of processors
    integer, parameter :: num_stcls_min = 5    !< Mininum # of stencils
    integer, parameter :: path_len = 400  !< Maximum path length
    integer, parameter :: name_len = 50   !< Maximum name length
    real(kind(0d0)), parameter :: dflt_real = -1d6 !< Default real value
    integer, parameter :: dflt_int = -100 !< Default integer value
    real(kind(0d0)), parameter :: sgm_eps = 1d-16 !< Segmentation tolerance
    character(LEN=path_len) :: case_dir             !< Case folder location
    !> @}

    ! Computational Domain Parameters ==========================================

    integer :: proc_rank !< Rank of the local processor

    !> @name Number of cells in the x-, y- and z-coordinate directions
    !> @{
    integer :: m, m_root
    integer :: n
    integer :: p
    !> @}

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

    !> @name Cell-boundary locations in the x-, y- and z-coordinate directions
    !> @{
    real(kind(0d0)), allocatable, dimension(:) :: x_cb, x_root_cb, y_cb, z_cb
    real(kind(0d0)), allocatable, dimension(:) :: coarse_x_cb, coarse_y_cb, coarse_z_cb
    !> @}

    !> @name Cell-center locations in the x-, y- and z-coordinate directions
    !> @{
    real(kind(0d0)), allocatable, dimension(:) :: x_cc, x_root_cc, y_cc, z_cc
    !> @}

    !> Cell-width distributions in the x-, y- and z-coordinate directions
    !> @{
    real(kind(0d0)), allocatable, dimension(:) :: dx, dy, dz
    !> @}

    integer :: buff_size !<
    !! Number of cells in buffer region. For the variables which feature a buffer
    !! region, this region is used to store information outside the computational
    !! domain based on the boundary conditions.

    integer :: t_step_start  !< First time-step directory
    integer :: t_step_stop   !< Last time-step directory
    integer :: t_step_save   !< Interval between consecutive time-step directory

    ! NOTE: The variables m_root, x_root_cb and x_root_cc contain the grid data
    ! of the defragmented computational domain. They are only used in 1D. For
    ! serial simulations, they are equal to m, x_cb and x_cc, respectively.

    ! ==========================================================================

    !> @name Simulation Algorithm Parameters
    !> @{
    integer :: model_eqns      !< Multicomponent flow model
    integer :: num_fluids      !< Number of different fluids present in the flow
    logical :: adv_alphan      !< Advection of the last volume fraction
    logical :: mpp_lim         !< Maximum volume fraction limiter
    integer :: sys_size        !< Number of unknowns in the system of equations
    integer :: weno_order      !< Order of accuracy for the WENO reconstruction
    logical :: mixture_err     !< Mixture error limiter
    logical :: alt_soundspeed  !< Alternate sound speed
    logical :: hypoelasticity  !< Turn hypoelasticity on
    !> @}

    !> @name Annotations of the structure, i.e. the organization, of the state vectors
    !> @{
    type(int_bounds_info) :: cont_idx              !< Indexes of first & last continuity eqns.
    type(int_bounds_info) :: mom_idx               !< Indexes of first & last momentum eqns.
    integer :: E_idx                               !< Index of energy equation
    type(int_bounds_info) :: adv_idx               !< Indexes of first & last advection eqns.
    type(int_bounds_info) :: internalEnergies_idx  !< Indexes of first & last internal energy eqns.
    type(bub_bounds_info) :: bub_idx               !< Indexes of first & last bubble variable eqns.
    integer :: gamma_idx                           !< Index of specific heat ratio func. eqn.
    integer :: alf_idx                             !< Index of specific heat ratio func. eqn.
    integer :: pi_inf_idx                          !< Index of liquid stiffness func. eqn.
    type(int_bounds_info) :: stress_idx            !< Indices of elastic stresses
    !> @}

    !> @name Boundary conditions in the x-, y- and z-coordinate directions
    !> @{
    type(int_bounds_info) :: bc_x, bc_y, bc_z
    !> @}

    logical :: parallel_io    !< Format of the data files

    integer, allocatable, dimension(:) :: proc_coords !<
    !! Processor coordinates in MPI_CART_COMM

    integer, allocatable, dimension(:) :: start_idx !<
    !! Starting cell-center index of local processor in global grid

#ifdef MFC_MPI

    type(mpi_io_var), public :: MPI_IO_DATA

#endif

    !> @name MPI info for parallel IO with Lustre file systems
    !> @{
    character(LEN=name_len) :: mpiiofs
    integer :: mpi_info_int
    !> @}

    integer, private :: ierr
    ! ==========================================================================

    type(physical_parameters), dimension(num_fluids_max) :: fluid_pp !<
    !! Database of the physical parameters of each of the fluids that is present
    !! in the flow. These include the stiffened gas equation of state parameters,
    !! the Reynolds numbers and the Weber numbers.

    ! ==========================================================================

    ! Formatted Database File(s) Structure Parameters ==========================

    integer :: format !< Format of the database file(s)

    logical :: coarsen_silo

    integer :: precision !< Floating point precision of the database file(s)

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
    logical :: schlieren_wrt
    !> @}

    !> @name Options for Fourier decomposition in the azimuthal direction if 3D
    !! cylindrical coordinates are used
    !> @{
    logical :: fourier_decomp
    !> @}

    real(kind(0d0)), dimension(num_fluids_max) :: schlieren_alpha    !<
    !! Amplitude coefficients of the numerical Schlieren function that are used
    !! to adjust the intensity of numerical Schlieren renderings for individual
    !! fluids. This enables waves and interfaces of varying strenghts and in all
    !! of the fluids to be made simulatenously visible on a single plot.

    integer :: fd_order !<
    !! The order of the finite-difference (fd) approximations of the first-order
    !! derivatives that need to be evaluated when vorticity and/or the numerical
    !! Schlieren function are to be outputted to the formatted database file(s).

    integer :: fd_number !<
    !! The finite-difference number is given by MAX(1, fd_order/2). Essentially,
    !! it is a measure of the half-size of the finite-difference stencil for the
    !! selected order of accuracy.

    ! ==========================================================================

    !> @name Reference parameters for Tait EOS
    !> @{
    real(kind(0d0)) :: rhoref, pref
    !> @}

    !> @name Bubble modeling variables and parameters
    !> @{
    integer :: nb
    real(kind(0d0)) :: R0ref
    real(kind(0d0)) :: Ca, Web, Re_inv
    real(kind(0d0)), dimension(:), allocatable :: weight, R0, V0
    logical :: bubbles
    logical :: polytropic
    logical :: polydisperse
    integer :: thermal  !< 1 = adiabatic, 2 = isotherm, 3 = transfer
    real(kind(0d0)) :: R_n, R_v, phi_vn, phi_nv, Pe_c, Tw, G
    real(kind(0d0)), dimension(:), allocatable :: k_n, k_v, pb0, mass_n0, mass_v0, Pe_T
    real(kind(0d0)), dimension(:), allocatable :: Re_trans_T, Re_trans_c, Im_trans_T, Im_trans_c, omegaN
    real(kind(0d0)) :: poly_sigma
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

    ! Mathematical and Physical Constants ======================================
    real(kind(0d0)), parameter :: pi = 3.141592653589793d0
    ! ==========================================================================

contains

    !> Assigns default values to user inputs prior to reading
        !!      them in. This allows for an easier consistency check of
        !!      these parameters once they are read from the input file.
    subroutine s_assign_default_values_to_user_inputs() ! ------------------

        integer :: i !< Generic loop iterator

        ! Logistics
        case_dir = ' '

        ! Computational domain parameters
        m = dflt_int; n = 0; p = 0
        m_root = dflt_int
        cyl_coord = .false.

        t_step_start = dflt_int
        t_step_stop = dflt_int
        t_step_save = dflt_int

        ! Simulation algorithm parameters
        model_eqns = dflt_int
        num_fluids = dflt_int
        adv_alphan = .false.
        weno_order = dflt_int
        mixture_err = .false.
        alt_soundspeed = .false.
        hypoelasticity = .false.

        bc_x%beg = dflt_int
        bc_x%end = dflt_int
        bc_y%beg = dflt_int
        bc_y%end = dflt_int
        bc_z%beg = dflt_int
        bc_z%end = dflt_int

        ! Fluids physical parameters
        do i = 1, num_fluids_max
            fluid_pp(i)%gamma = dflt_real
            fluid_pp(i)%pi_inf = dflt_real
            fluid_pp(i)%G = dflt_real
        end do

        ! Formatted database file(s) structure parameters
        format = dflt_int

        precision = dflt_int

        coarsen_silo = .false.

        alpha_rho_wrt = .false.
        rho_wrt = .false.
        mom_wrt = .false.
        vel_wrt = .false.
        flux_lim = dflt_int
        flux_wrt = .false.
        parallel_io = .false.
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
        schlieren_wrt = .false.

        schlieren_alpha = dflt_real

        fourier_decomp = .false.

        fd_order = dflt_int

        ! Tait EOS
        rhoref = dflt_real
        pref = dflt_real

        ! Bubble modeling
        bubbles = .false.
        R0ref = dflt_real
        nb = dflt_int
        polydisperse = .false.
        poly_sigma = dflt_real

    end subroutine s_assign_default_values_to_user_inputs ! ----------------

    !>  Computation of parameters, allocation procedures, and/or
        !!      any other tasks needed to properly setup the module
    subroutine s_initialize_global_parameters_module() ! ----------------------

        integer :: i, fac

        ! Setting m_root equal to m in the case of a 1D serial simulation
        if (num_procs == 1 .and. n == 0) m_root = m

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
                alf_idx = 0
            end if

            if (bubbles) then
                bub_idx%beg = sys_size + 1
                bub_idx%end = sys_size + 2*nb
                if (polytropic .neqv. .true.) then
                    bub_idx%end = sys_size + 4*nb
                end if
                sys_size = bub_idx%end

                allocate (bub_idx%rs(nb), bub_idx%vs(nb))
                allocate (bub_idx%ps(nb), bub_idx%ms(nb))
                allocate (weight(nb), R0(nb), V0(nb))

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
                    weight(:) = 1d0
                    R0(:) = 1d0
                    V0(:) = 0d0
                else if (nb > 1) then
                    call s_simpson(nb)
                    V0(:) = 0d0
                else
                    stop 'Invalid value of nb'
                end if

                if (polytropic .neqv. .true.) then
                    call s_initialize_nonpoly
                else
                    rhoref = 1.d0
                    pref = 1.d0
                end if
            end if

            if (hypoelasticity) then
                stress_idx%beg = sys_size + 1
                stress_idx%end = sys_size + (num_dims*(num_dims + 1))/2
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
            if (adv_alphan .neqv. .true.) adv_idx%end = adv_idx%end - 1
            internalEnergies_idx%beg = adv_idx%end + 1
            internalEnergies_idx%end = adv_idx%end + num_fluids
            sys_size = internalEnergies_idx%end

        else if (model_eqns == 4) then
            cont_idx%beg = 1 ! one continuity equation
            cont_idx%end = 1 !num_fluids
            mom_idx%beg = cont_idx%end + 1 ! one momentum equation in each
            mom_idx%end = cont_idx%end + num_dims
            E_idx = mom_idx%end + 1 ! one energy equation
            adv_idx%beg = E_idx + 1
            adv_idx%end = adv_idx%beg !one volume advection equation
            alf_idx = adv_idx%end
            sys_size = alf_idx !adv_idx%end

            if (bubbles) then
                bub_idx%beg = sys_size + 1
                bub_idx%end = sys_size + 2*nb
                if (polytropic .neqv. .true.) then
                    bub_idx%end = sys_size + 4*nb
                end if
                sys_size = bub_idx%end

                allocate (bub_idx%rs(nb), bub_idx%vs(nb))
                allocate (bub_idx%ps(nb), bub_idx%ms(nb))
                allocate (weight(nb), R0(nb), V0(nb))

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
                    weight(:) = 1d0
                    R0(:) = 1d0
                    V0(:) = 0d0
                else if (nb > 1) then
                    call s_simpson(nb)
                    V0(:) = 0d0
                else
                    stop 'Invalid value of nb'
                end if

                if (polytropic .neqv. .true.) then
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

        if (any(omega_wrt) .or. schlieren_wrt) then
            fd_number = max(1, fd_order/2)
            buff_size = buff_size + fd_number
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

        end if

        if (coarsen_silo) then
            allocate (coarse_x_cb(-1 - offset_x%beg:(m/2) + offset_x%end))
            if (n > 0) then
                allocate (coarse_y_cb(-1 - offset_y%beg:(n/2) + offset_y%end))
                if (p > 0) allocate (coarse_z_cb(-1 - offset_z%beg:(p/2) + offset_z%end))
            end if
        end if

        if (cyl_coord .neqv. .true.) then ! Cartesian grid
            grid_geometry = 1
        elseif (cyl_coord .and. p == 0) then ! Axisymmetric cylindrical grid
            grid_geometry = 2
        else ! Fully 3D cylindrical grid
            grid_geometry = 3
        end if

    end subroutine s_initialize_global_parameters_module ! --------------------

    !> Subroutine to initialize variable for non-polytropic gas modeling processes
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

    !> Subroutine to compute the transfer coefficient for non-polytropic gas modeling
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

    !> Subroutine to initialize parallel infrastructure
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

    !> Deallocation procedures for the module
    subroutine s_finalize_global_parameters_module() ! -------------------

        integer :: i

        ! Deallocating the grid variables for the x-coordinate direction
        deallocate (x_cb, x_cc, dx)

        ! Deallocating grid variables for the y- and z-coordinate directions
        if (n > 0) then

            deallocate (y_cb, y_cc, dy)

            if (p > 0) deallocate (z_cb, z_cc, dz)

            ! Deallocating the grid variables, only used for the 1D simulations,
            ! and containing the defragmented computational domain grid data
        else

            deallocate (x_root_cb, x_root_cc)

        end if

        if (coarsen_silo) then
            deallocate (coarse_x_cb)
            if (n > 0) then
                deallocate (coarse_y_cb)
                if (p > 0) deallocate (coarse_z_cb)
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

    end subroutine s_finalize_global_parameters_module ! -----------------

    !> Computes the bubble number density n from the conservative variables
        !! @param vftmp is the void fraction
        !! @param nRtmp is the bubble number  density times the bubble radii
        !! @param ntmp is the output number bubble density
    subroutine s_comp_n_from_cons(vftmp, nRtmp, ntmp)

        real(kind(0.d0)), intent(IN) :: vftmp
        real(kind(0.d0)), dimension(nb), intent(IN) :: nRtmp
        real(kind(0.d0)), intent(OUT) :: ntmp
        real(kind(0.d0)) :: nR3

        call s_quad(nRtmp**3, nR3)  !returns itself if NR0 = 1
        ntmp = DSQRT((4.d0*pi/3.d0)*nR3/vftmp)

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

        call s_quad(Rtmp**3, R3)  !returns itself if NR0 = 1
        ntmp = (3.d0/(4.d0*pi))*vftmp/R3

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
        !! @param Npt is the number of bins that represent the polydisperse bubble population
    subroutine s_simpson(Npt)

        integer, intent(IN) :: Npt
        integer :: ir
        real(kind(0.d0)) :: R0mn
        real(kind(0.d0)) :: R0mx
        real(kind(0.d0)) :: dphi
        real(kind(0.d0)) :: tmp
        real(kind(0.d0)) :: sd
        real(kind(0.d0)), dimension(Npt) :: phi

        ! nondiml. min. & max. initial radii for numerical quadrature
        !sd   = 0.05D0
        !R0mn = 0.75D0
        !R0mx = 1.3D0

        sd = poly_sigma
        R0mn = 0.8d0*DEXP(-2.8d0*sd)
        R0mx = 0.2d0*DEXP(9.5d0*sd) + 1.d0

        ! phi = ln( R0 ) & return R0
        do ir = 1, Npt
            phi(ir) = DLOG(R0mn) &
                      + dble(ir - 1)*DLOG(R0mx/R0mn)/dble(Npt - 1)
            R0(ir) = DEXP(phi(ir))
        end do
        dphi = phi(2) - phi(1)

        ! weights for quadrature using Simpson's rule
        do ir = 2, Npt - 1
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
        tmp = DEXP(-0.5d0*(phi(Npt)/sd)**2)/DSQRT(2.d0*pi)/sd
        weight(Npt) = tmp*dphi/3.d0

    end subroutine s_simpson

end module m_global_parameters
