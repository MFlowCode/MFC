!!       __  _______________
!!      /  |/  / ____/ ____/
!!     / /|_/ / /_  / /
!!    / /  / / __/ / /___
!!   /_/  /_/_/    \____/
!!
!!  This file is part of MFC.
!!
!!  MFC is the legal property of its developers, whose names
!!  are listed in the copyright file included with this source
!!  distribution.
!!
!!  MFC is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published
!!  by the Free Software Foundation, either version 3 of the license
!!  or any later version.
!!
!!  MFC is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with MFC (LICENSE).
!!  If not, see <http://www.gnu.org/licenses/>.

!>
!! @file m_global_parameters.f90
!! @brief Contains module m_global_parameters
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief The module contains all of the parameters describing the program
!!              logistics, the computational domain and the simulation algorithm.
!!              Additionally, for the volume fraction model, physical parameters
!!              of each of the fluids present in the flow are located here. They
!!              include stiffened gas equation of state parameters, the Reynolds
!!              numbers and the Weber numbers.
module m_global_parameters

    ! Dependencies =============================================================
    use mpi                    !< Message passing interface (MPI) module

    use m_derived_types        !< Definitions of the derived types

    ! ==========================================================================

    implicit none

    ! Logistics ================================================================
    integer                    :: num_procs             !< Number of processors
    integer, parameter :: num_stcls_min = 5     !< Mininum # of stencils
    integer, parameter :: path_len = 400   !< Maximum path length
    integer, parameter :: name_len = 50    !< Maximum name length
    character, parameter :: dflt_char = ' '   !< Default string value
    real(kind(0d0)), parameter :: dflt_real = -1d6  !< Default real value
    integer, parameter :: dflt_int = -100  !< Default integer value
    real(kind(0d0)), parameter :: sgm_eps = 1d-16 !< Segmentation tolerance
    integer, parameter :: fourier_rings = 5     !< Fourier filter ring limit
    character(LEN=path_len)  :: case_dir              !< Case folder location
    logical                    :: run_time_info         !< Run-time output flag
    logical                    :: debug                 !< Debug mode print statements
    integer                    :: t_step_old            !< Existing IC/grid folder
    real(kind(0d0)), parameter :: small_alf = 1d-7 !< Small alf tolerance
    ! ==========================================================================

    ! Computational Domain Parameters ==========================================

    integer :: proc_rank !< Rank of the local processor

    !> @name Number of cells in the x-, y- and z-directions, respectively
    !> @{
    integer :: m, n, p
    !> @}

    !> @name Global number of cells in each direction
    !> @{
    integer :: m_glb, n_glb, p_glb
    !> @}

    !> @name Cylindrical coordinates (either axisymmetric or full 3D)
    !> @{
    logical :: cyl_coord
    integer :: grid_geometry
    !> @}

    !> @name Cell-boundary (CB) locations in the x-, y- and z-directions, respectively
    !> @{
    real(kind(0d0)), target, allocatable, dimension(:) :: x_cb, y_cb, z_cb
    !> @}

    !> @name Cell-center (CC) locations in the x-, y- and z-directions, respectively
    !> @{
    real(kind(0d0)), target, allocatable, dimension(:) :: x_cc, y_cc, z_cc
    !> @}

    !> @name Cell-width distributions in the x-, y- and z-directions, respectively
    !> @{
    real(kind(0d0)), target, allocatable, dimension(:) :: dx, dy, dz
    !> @}

    real(kind(0d0)) :: dt !< Size of the time-step

    !> @name Starting time-step iteration, stopping time-step iteration and the number
    !! of time-step iterations between successive solution backups, respectively
    !> @{
    integer :: t_step_start, t_step_stop, t_step_save
    !> @}

    ! ==========================================================================

    ! Simulation Algorithm Parameters ==========================================
    integer         :: model_eqns     !< Multicomponent flow model
    integer         :: num_dims       !< Number of spatial dimensions
    integer         :: num_fluids     !< Number of fluids in the flow
    logical         :: adv_alphan     !< Advection of the last volume fraction
    logical         :: mpp_lim        !< Mixture physical parameters (MPP) limits
    integer         :: time_stepper   !< Time-stepper algorithm
    integer         :: weno_vars      !< WENO-reconstructed state variables type
    integer         :: weno_order     !< Order of the WENO reconstruction
    integer         :: weno_polyn     !< Degree of the WENO polynomials (polyn)
    real(kind(0d0)) :: weno_eps       !< Binding for the WENO nonlinear weights
    logical         :: mapped_weno    !< WENO with mapping of nonlinear weights
    logical         :: mp_weno        !< Monotonicity preserving (MP) WENO
    logical         :: weno_Re_flux   !< WENO reconstruct velocity gradients for viscous stress tensor
    integer         :: riemann_solver !< Riemann solver algorithm
    integer         :: wave_speeds    !< Wave speeds estimation method
    integer         :: avg_state      !< Average state evaluation method
    logical         :: alt_soundspeed !< Alternate mixture sound speed
    logical         :: regularization !< Regularization terms of Tiwari (2013)
    real(kind(0d0)) :: reg_eps        !< User-defined interface thickness parameter for regularization terms
    logical         :: null_weights   !< Null undesired WENO weights
    logical         :: mixture_err    !< Mixture properties correction

    integer         :: cpu_start, cpu_end, cpu_rate

    !> @name Boundary conditions (BC) in the x-, y- and z-directions, respectively
    !> @{
    type(bounds_info) :: bc_x, bc_y, bc_z
    type(bounds_info) :: bc_x_glb, bc_y_glb, bc_z_glb
    !> @}

    logical :: parallel_io !< Format of the data files
    integer :: precision !< Precision of output files

    integer, allocatable, dimension(:) :: proc_coords !<
    !! Processor coordinates in MPI_CART_COMM

    integer, allocatable, dimension(:) :: start_idx !<
    !! Starting cell-center index of local processor in global grid

    type(mpi_io_var), public :: MPI_IO_DATA

    !> @name MPI info for parallel IO with Lustre file systems
    !> @{
    character(LEN=name_len) :: mpiiofs
    integer :: mpi_info_int
    !> @}

    integer, private :: ierr

    !> @name Annotations of the structure of the state and flux vectors in terms of the
    !! size and the configuration of the system of equations to which they belong
    !> @{
    integer           :: sys_size                  !< Number of unknowns in system of eqns.
    type(bounds_info) :: cont_idx                  !< Indexes of first & last continuity eqns.
    type(bounds_info) :: mom_idx                   !< Indexes of first & last momentum eqns.
    integer           :: E_idx                     !< Index of energy equation
    type(bounds_info) :: adv_idx                   !< Indexes of first & last advection eqns.
    type(bounds_info) :: internalEnergies_idx      !< Indexes of first & last internal energy eqns.
    type(bub_bounds_info) :: bub_idx               !< Indexes of first & last bubble variable eqns.
    integer               :: alf_idx               !< Index of void fraction
    integer           :: gamma_idx                 !< Index of specific heat ratio func. eqn.
    integer           :: pi_inf_idx                !< Index of liquid stiffness func. eqn.
    !> @}

    !> @name The number of fluids, along with their identifying indexes, respectively,
    !! for which viscous effects, e.g. the shear and/or the volume Reynolds (Re)
    !! numbers, will be non-negligible.
    !> @{
    integer, dimension(2)   :: Re_size
    integer, allocatable, dimension(:, :) :: Re_idx
    !> @}

    !> @name The coordinate direction indexes and flags (flg), respectively, for which
    !! the configurations will be determined with respect to a working direction
    !! and that will be used to isolate the contributions, in that direction, in
    !! the dimensionally split system of equations.
    !> @{
    integer, dimension(3) :: dir_idx
    real(kind(0d0)), dimension(3) :: dir_flg
    !> @}


    integer :: buff_size !<
    !! The number of cells that are necessary to be able to store enough boundary
    !! conditions data to march the solution in the physical computational domain
    !! to the next time-step.

    ! END: Simulation Algorithm Parameters =====================================

    ! Fluids Physical Parameters ===============================================

    type(physical_parameters), dimension(num_fluids_max) :: fluid_pp !<
    !! Database of the physical parameters of each of the fluids that is present
    !! in the flow. These include the stiffened gas equation of state parameters,
    !! the Reynolds numbers and the Weber numbers.

    ! ==========================================================================

    integer :: fd_order !<
    !! The order of the finite-difference (fd) approximations of the first-order
    !! derivatives that need to be evaluated when the CoM or flow probe data
    !! files are to be written at each time step

    integer :: fd_number !<
    !! The finite-difference number is given by MAX(1, fd_order/2). Essentially,
    !! it is a measure of the half-size of the finite-difference stencil for the
    !! selected order of accuracy.

    logical, dimension(num_fluids_max) :: com_wrt, cb_wrt
    logical :: probe_wrt
    logical :: integral_wrt
    integer :: num_probes
    integer :: num_integrals
    type(probe_parameters), dimension(num_probes_max) :: probe
    type(integral_parameters), dimension(num_probes_max) :: integral
    real(kind(0d0)), dimension(5) :: threshold_mf
    integer, dimension(5) :: moment_order

    !> @name Reference density and pressure for Tait EOS
    !> @{
    real(kind(0d0)) :: rhoref, pref
    !> @}

    !> @name Bubble modeling
    !> @{
    integer         :: nb       !< Number of eq. bubble sizes
    real(kind(0d0)) :: R0ref    !< Reference bubble size
    real(kind(0d0)) :: Ca       !< Cavitation number
    real(kind(0d0)) :: Web      !< Weber number
    real(kind(0d0)) :: Re_inv   !< Inverse Reynolds number
    real(kind(0d0)), dimension(:), allocatable :: weight !< Simpson quadrature weights
    real(kind(0d0)), dimension(:), allocatable :: R0     !< Bubble sizes
    real(kind(0d0)), dimension(:), allocatable :: V0     !< Bubble velocities
    logical         :: bubbles      !< Bubbles on/off
    logical         :: polytropic   !< Polytropic  switch
    logical         :: polydisperse !< Polydisperse bubbles

    integer         :: bubble_model !< Gilmore or Keller--Miksis bubble model
    integer         :: thermal      !< Thermal behavior. 1 = adiabatic, 2 = isotherm, 3 = transfer
    real(kind(0d0)), allocatable, dimension(:, :, :) :: ptil  !< Pressure modification
    real(kind(0d0)) :: poly_sigma  !< log normal sigma for polydisperse PDF

    logical         :: qbmm !< Quadrature moment method
    integer         :: nmom !< Number of carried moments per R0 location
    integer         :: nnode !< Number of QBMM nodes
    integer         :: nmomsp !< Number of moments required by ensemble-averaging
    integer         :: nmomtot !< Total number of carried moments moments/transport equations
    integer         :: R0_type

    type(scalar_field), allocatable, dimension(:) :: mom_sp
    type(scalar_field), allocatable, dimension(:, :, :) :: mom_3d
    !> @}

    !> @name Physical bubble parameters (see Ando 2010, Preston 2007)
    !> @{
    real(kind(0d0)) :: R_n, R_v, phi_vn, phi_nv, Pe_c, Tw, pv, M_n, M_v
    real(kind(0d0)), dimension(:), allocatable :: k_n, k_v, pb0, mass_n0, mass_v0, Pe_T
    real(kind(0d0)), dimension(:), allocatable :: Re_trans_T, Re_trans_c, Im_trans_T, Im_trans_c, omegaN
    real(kind(0d0)) :: mul0, ss, gamma_v, mu_v
    real(kind(0d0)) :: gamma_m, gamma_n, mu_n
    real(kind(0d0)) :: gam
    !> @}

    !> @name Acoustic monopole parameters
    !> @{
    logical         :: monopole !< Monopole switch
    type(mono_parameters), dimension(num_probes_max) :: mono !< Monopole parameters
    integer         :: num_mono !< Number of monopoles
    !> @}

    real(kind(0d0)) :: mytime       !< Current simulation time
    real(kind(0d0)) :: finaltime    !< Final simulation time
    real(kind(0d0)) :: t_tol        !< Tolerance
    real(kind(0d0)) :: dt0          !< Initial time step size

    ! ======================================================================

    ! Mathematical and Physical Constants ======================================
    ! REAL(KIND(0d0)), PARAMETER :: pi = 3.141592653589793d0 !< Pi
    real(kind(0d0)), parameter :: pi = 3.14159265358979311599796 !< Pi
    ! ==========================================================================

contains

    !> Assigns default values to the user inputs before reading
        !!  them in. This enables for an easier consistency check of
        !!  these parameters once they are read from the input file.
    subroutine s_assign_default_values_to_user_inputs() ! ------------------

        integer :: i, j !< Generic loop iterator

        ! Logistics
        case_dir = dflt_char
        run_time_info = .false.
        t_step_old = dflt_int

        debug = .false.

        ! Computational domain parameters
        m = dflt_int; n = dflt_int; p = dflt_int

        cyl_coord = .false.

        dt = dflt_real
        t_tol = dflt_real

        t_step_start = dflt_int
        t_step_stop = dflt_int
        t_step_save = dflt_int

        ! Simulation algorithm parameters
        model_eqns = dflt_int
        num_fluids = dflt_int
        adv_alphan = .false.
        mpp_lim = .false.
        time_stepper = dflt_int
        weno_vars = dflt_int
        weno_order = dflt_int
        weno_eps = dflt_real
        mapped_weno = .false.
        mp_weno = .false.
        weno_Re_flux = .false.
        riemann_solver = dflt_int
        wave_speeds = dflt_int
        avg_state = dflt_int
        alt_soundspeed = .false.
        regularization = .false.
        reg_eps = dflt_real
        null_weights = .false.
        mixture_err = .false.
        parallel_io = .false.
        precision = 2

        bc_x%beg = dflt_int; bc_x%end = dflt_int
        bc_y%beg = dflt_int; bc_y%end = dflt_int
        bc_z%beg = dflt_int; bc_z%end = dflt_int

        ! Fluids physical parameters
        do i = 1, num_fluids_max
            fluid_pp(i)%gamma = dflt_real
            fluid_pp(i)%pi_inf = dflt_real
            fluid_pp(i)%Re(:) = dflt_real
            fluid_pp(i)%mul0 = dflt_real
            fluid_pp(i)%ss = dflt_real
            fluid_pp(i)%pv = dflt_real
            fluid_pp(i)%gamma_v = dflt_real
            fluid_pp(i)%M_v = dflt_real
            fluid_pp(i)%mu_v = dflt_real
            fluid_pp(i)%k_v = dflt_real
        end do

        ! Tait EOS
        rhoref = dflt_real
        pref = dflt_real

        ! Bubble modeling
        bubbles = .false.
        bubble_model = 1
        polytropic = .true.
        polydisperse = .false.
        thermal = dflt_int
        R0ref = dflt_real
        nb = dflt_int
        R0_type = dflt_int

        ! User inputs for qbmm for simulation code
        qbmm = .false.
        nnode = 1

        Ca = dflt_real
        Re_inv = dflt_real
        Web = dflt_real
        poly_sigma = dflt_real

        ! Monopole source
        monopole = .false.
        num_mono = 1

        do j = 1, num_probes_max
            do i = 1, 3
                mono(j)%loc(i) = dflt_real
            end do
            mono(j)%mag = dflt_real
            mono(j)%length = dflt_real
            mono(j)%delay = dflt_real
            mono(j)%dir = 1.d0
            mono(j)%npulse = 1.d0
            mono(j)%pulse = 1
            mono(j)%support = 1
        end do

        fd_order = dflt_int
        com_wrt = .false.
        cb_wrt = .false.
        probe_wrt = .false.
        integral_wrt = .false.
        num_probes = dflt_int
        num_integrals = dflt_int

        do i = 1, num_probes_max
            probe(i)%x = dflt_real
            probe(i)%y = dflt_real
            probe(i)%z = dflt_real
        end do

        do i = 1, num_probes_max
            integral(i)%xmin = dflt_real
            integral(i)%xmax = dflt_real
            integral(i)%ymin = dflt_real
            integral(i)%ymax = dflt_real
            integral(i)%ymin = dflt_real
            integral(i)%ymax = dflt_real
        end do

        do i = 1, 5
            threshold_mf(i) = dflt_real
            moment_order(i) = dflt_int
        end do

    end subroutine s_assign_default_values_to_user_inputs ! ----------------

    !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_global_parameters_module() ! -------------------

        integer :: tmp_idx !< Temporary indexes storage
        integer :: i, j !< Generic loop iterators
        integer :: k !< Generic counter
        integer :: fac
        integer :: i1, i2, i3

        type(bounds_info) :: ix, iy, iz

        ! Determining the degree of the WENO polynomials
        weno_polyn = (weno_order - 1)/2

        ! Initializing the number of fluids for which viscous effects will
        ! be non-negligible, the number of distinctive material interfaces
        ! for which surface tension will be important and also, the number
        ! of fluids for which the physical and geometric curvatures of the
        ! interfaces will be computed
        Re_size = 0

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

                if ((adv_alphan .neqv. .true.) .and. &
                    (num_fluids > 1)) adv_idx%end = adv_idx%end - 1

                sys_size = adv_idx%end

                if (bubbles) then
                    alf_idx = adv_idx%end
                else
                    alf_idx = 0
                end if

                if (bubbles) then
                    bub_idx%beg = sys_size + 1
                    if (qbmm) then
                        nmomsp = 4 !number of special moments
                        if (nnode == 4) then
                            nmom = 6
                            nmomtot = nmom*nb
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

                        do i = 1, nb
                            do j = 1, nmom
                                bub_idx%moms(i, j) = bub_idx%beg + (j - 1) + (i - 1)*nmom
                            end do
                            bub_idx%rs(i) = bub_idx%moms(i, 1)
                            bub_idx%vs(i) = bub_idx%moms(i, 2)
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
            fd_number = max(1, fd_order/2)
            buff_size = buff_size + fd_number
        end if

        if (cyl_coord .neqv. .true.) then ! Cartesian grid
            grid_geometry = 1
        elseif (cyl_coord .and. p == 0) then ! Axisymmetric cylindrical grid
            grid_geometry = 2
        else ! Fully 3D cylindrical grid
            grid_geometry = 3
        end if

        ! Allocating grid variables for the x-, y- and z-directions
        allocate (x_cb(-1 - buff_size:m + buff_size))
        allocate (x_cc(-buff_size:m + buff_size))
        allocate (dx(-buff_size:m + buff_size))

        if (n == 0) return; allocate (y_cb(-1 - buff_size:n + buff_size))
        allocate (y_cc(-buff_size:n + buff_size))
        allocate (dy(-buff_size:n + buff_size))

        if (p == 0) return; allocate (z_cb(-1 - buff_size:p + buff_size))
        allocate (z_cc(-buff_size:p + buff_size))
        allocate (dz(-buff_size:p + buff_size))

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

        num_dims = 1 + min(1, n) + min(1, p)

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

        integer :: i

        ! Deallocating the variables bookkeeping the indexes of any viscous
        ! fluids and any pairs of fluids whose interfaces supported effects
        ! of surface tension
        if (any(Re_size > 0)) deallocate (Re_idx)

        ! Deallocating grid variables for the x-, y- and z-directions
        deallocate (x_cb, x_cc, dx)

        if (n == 0) return; deallocate (y_cb, y_cc, dy)

        if (p == 0) return; deallocate (z_cb, z_cc, dz)

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

        real(kind(0.d0)), intent(IN) :: vftmp
        real(kind(0.d0)), dimension(nb), intent(IN) :: nRtmp
        real(kind(0.d0)), intent(OUT) :: ntmp
        real(kind(0.d0)) :: nR3
        integer :: i

        call s_quad(nRtmp**3d0, nR3)

        if (nR3 < 0d0) then
            ! DO i = 1,nb
            ! IF (nRtmp(i) < small_alf) THEN
            ! nRtmp(i) = small_alf
            ! END IF
            ! END DO
            ! nR3 = 1.d-12
            print *, vftmp, nR3, nRtmp(:)
            stop 'nR3 is negative'
        end if
        if (vftmp < 0d0) then
            ! vftmp = small_alf
            ! ntmp = DSQRT( (4.d0*pi/3.d0)*nR3/1.d-12 )
            print *, vftmp, nR3, nRtmp(:)
            stop 'vf negative'
        end if

        ntmp = DSQRT((4.d0*pi/3.d0)*nR3/vftmp)

    end subroutine s_comp_n_from_cons

    !> Computes the bubble number density n from the primitive variables
        !!  \f$ n = \sqrt{ \frac{3}{4 \pi} } \frac{ \alpha }{ R^3} \f$
        !! @param vftmp is the void fraction
        !! @param Rtmp is the  bubble radii
        !! @param ntmp is the output number bubble density
    subroutine s_comp_n_from_prim(vftmp, Rtmp, ntmp)

        real(kind(0.d0)), intent(IN) :: vftmp
        real(kind(0.d0)), dimension(nb), intent(IN) :: Rtmp
        real(kind(0.d0)), intent(OUT) :: ntmp
        real(kind(0.d0)) :: R3

        call s_quad(Rtmp**3d0, R3)

        if (R3 < 0d0) then
            print *, vftmp, R3, Rtmp(:)
            stop 'R3 is negative'
        end if
        if (vftmp < 0d0) then
            print *, vftmp, R3, Rtmp(:)
            stop 'vf negative'
        end if

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
