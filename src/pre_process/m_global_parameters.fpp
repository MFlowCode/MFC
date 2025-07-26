!>
!! @file m_global_parameters.f90
!! @brief Contains module m_global_parameters

#:include 'case.fpp'

!> @brief This module contains all of the parameters characterizing the
!!              computational domain, simulation algorithm, initial condition
!!              and the stiffened equation of state.
module m_global_parameters

#ifdef MFC_MPI
    use mpi                     ! Message passing interface (MPI) module
#endif

    use m_derived_types         ! Definitions of the derived types

    use m_helper_basic          ! Functions to compare floating point numbers

    use m_thermochem, only: num_species

    implicit none

    ! Logistics
    integer :: num_procs            !< Number of processors
    character(LEN=path_len) :: case_dir             !< Case folder location
    logical :: old_grid             !< Use existing grid data
    logical :: old_ic, non_axis_sym               !< Use existing IC data
    integer :: t_step_old, t_step_start           !< Existing IC/grid folder

    logical :: cfl_adap_dt, cfl_const_dt, cfl_dt
    integer :: n_start, n_start_old

    ! Computational Domain Parameters

    integer :: proc_rank !< Rank of the local processor

    !! Number of cells in the x-, y- and z-coordinate directions
    integer :: m
    integer :: n
    integer :: p

    !> @name Max and min number of cells in a direction of each combination of x-,y-, and z-
    type(cell_num_bounds) :: cells_bounds

    integer(8) :: nGlobal !< Global number of cells in the domain

    integer :: m_glb, n_glb, p_glb !< Global number of cells in each direction

    integer :: num_dims !< Number of spatial dimensions
    integer :: num_vels !< Number of velocity components (different from num_dims for mhd)

    logical :: cyl_coord
    integer :: grid_geometry !< Cylindrical coordinates (either axisymmetric or full 3D)

    real(wp), allocatable, dimension(:) :: x_cc, y_cc, z_cc !<
    !! Locations of cell-centers (cc) in x-, y- and z-directions, respectively

    real(wp), allocatable, dimension(:) :: x_cb, y_cb, z_cb !<
    !! Locations of cell-boundaries (cb) in x-, y- and z-directions, respectively

    real(wp) :: dx, dy, dz !<
    !! Minimum cell-widths in the x-, y- and z-coordinate directions

    type(bounds_info) :: x_domain, y_domain, z_domain !<
    !! Locations of the domain bounds in the x-, y- and z-coordinate directions

    logical :: stretch_x, stretch_y, stretch_z !<
    !! Grid stretching flags for the x-, y- and z-coordinate directions

    ! Parameters of the grid stretching function for the x-, y- and z-coordinate
    ! directions. The "a" parameters are a measure of the rate at which the grid
    ! is stretched while the remaining parameters are indicative of the location
    ! on the grid at which the stretching begins.
    real(wp) :: a_x, a_y, a_z
    integer :: loops_x, loops_y, loops_z
    real(wp) :: x_a, y_a, z_a
    real(wp) :: x_b, y_b, z_b

    ! Simulation Algorithm Parameters
    integer :: model_eqns            !< Multicomponent flow model
    logical :: relax                 !< activate phase change
    integer :: relax_model           !< Relax Model
    real(wp) :: palpha_eps           !< trigger parameter for the p relaxation procedure, phase change model
    real(wp) :: ptgalpha_eps         !< trigger parameter for the pTg relaxation procedure, phase change model
    integer :: num_fluids            !< Number of different fluids present in the flow
    logical :: mpp_lim               !< Alpha limiter
    integer :: sys_size              !< Number of unknowns in the system of equations
    integer :: recon_type            !< Reconstruction Type
    integer :: weno_polyn            !< Degree of the WENO polynomials (polyn)
    integer :: muscl_polyn           !< Degree of the MUSCL polynomials (polyn)
    integer :: weno_order            !< Order of accuracy for the WENO reconstruction
    integer :: muscl_order           !< Order of accuracy for the MUSCL reconstruction
    logical :: hypoelasticity        !< activate hypoelasticity
    logical :: hyperelasticity       !< activate hyperelasticity
    logical :: elasticity            !< elasticity modeling, true for hyper or hypo
    logical :: mhd                   !< Magnetohydrodynamics
    logical :: relativity            !< Relativity for RMHD
    integer :: b_size                !< Number of components in the b tensor
    integer :: tensor_size           !< Number of components in the nonsymmetric tensor
    logical :: pre_stress            !< activate pre_stressed domain
    logical :: cont_damage           !< continuum damage modeling
    logical :: igr                   !< Use information geometric regularization
    integer :: igr_order             !< IGR reconstruction order
    logical, parameter :: chemistry = .${chemistry}$. !< Chemistry modeling

    ! Annotations of the structure, i.e. the organization, of the state vectors
    type(int_bounds_info) :: cont_idx              !< Indexes of first & last continuity eqns.
    type(int_bounds_info) :: mom_idx               !< Indexes of first & last momentum eqns.
    integer :: E_idx                               !< Index of total energy equation
    integer :: alf_idx                             !< Index of void fraction
    integer :: n_idx                               !< Index of number density
    type(int_bounds_info) :: adv_idx               !< Indexes of first & last advection eqns.
    type(int_bounds_info) :: internalEnergies_idx  !< Indexes of first & last internal energy eqns.
    type(bub_bounds_info) :: bub_idx               !< Indexes of first & last bubble variable eqns.
    integer :: gamma_idx                           !< Index of specific heat ratio func. eqn.
    integer :: pi_inf_idx                          !< Index of liquid stiffness func. eqn.
    type(int_bounds_info) :: B_idx                 !< Indexes of first and last magnetic field eqns.
    type(int_bounds_info) :: stress_idx            !< Indexes of elastic shear stress eqns.
    type(int_bounds_info) :: xi_idx                !< Indexes of first and last reference map eqns.
    integer :: c_idx                               !< Index of the color function
    type(int_bounds_info) :: species_idx           !< Indexes of first & last concentration eqns.
    integer :: damage_idx                          !< Index of damage state variable (D) for continuum damage model

    ! Cell Indices for the (local) interior points (O-m, O-n, 0-p).
    ! Stands for "InDices With BUFFer".
    type(int_bounds_info) :: idwint(1:3)

    ! Cell Indices for the entire (local) domain. In simulation and post_process,
    ! this includes the buffer region. idwbuff and idwint are the same otherwise.
    ! Stands for "InDices With BUFFer".
    type(int_bounds_info) :: idwbuff(1:3)

    type(int_bounds_info) :: bc_x, bc_y, bc_z !<
    !! Boundary conditions in the x-, y- and z-coordinate directions

    integer :: shear_num !! Number of shear stress components
    integer, dimension(3) :: shear_indices !<
    !! Indices of the stress components that represent shear stress
    integer :: shear_BC_flip_num !<
    !! Number of shear stress components to reflect for boundary conditions
    integer, dimension(3, 2) :: shear_BC_flip_indices !<
    !! Indices of shear stress components to reflect for boundary conditions.
    !! Size: (1:3, 1:shear_BC_flip_num) for (x/y/z, [indices])

    logical :: parallel_io !< Format of the data files
    logical :: file_per_process !< type of data output
    integer :: precision !< Precision of output files

    logical :: mixlayer_vel_profile !< Set hyperbolic tangent streamwise velocity profile
    real(wp) :: mixlayer_vel_coef !< Coefficient for the hyperbolic tangent streamwise velocity profile
    logical :: mixlayer_perturb !< Superimpose instability waves to surrounding fluid flow
    integer :: mixlayer_perturb_nk  !< Number of Fourier modes for perturbation with mixlayer_perturb flag
    real(wp) :: mixlayer_perturb_k0  !< Peak wavenumber of prescribed energy spectra with mixlayer_perturb flag
                                     !! Default value (k0 = 0.4446) is most unstable mode obtained from linear stability analysis
                                     !! See Michalke (1964, JFM) for details

    real(wp) :: pi_fac !< Factor for artificial pi_inf

    logical :: viscous
    logical :: bubbles_lagrange

    ! Perturb density of surrounding air so as to break symmetry of grid
    logical :: perturb_flow
    integer :: perturb_flow_fluid   !< Fluid to be perturbed with perturb_flow flag
    real(wp) :: perturb_flow_mag   !< Magnitude of perturbation with perturb_flow flag
    logical :: perturb_sph
    integer :: perturb_sph_fluid    !< Fluid to be perturbed with perturb_sph flag
    real(wp), dimension(num_fluids_max) :: fluid_rho

    logical :: elliptic_smoothing
    integer :: elliptic_smoothing_iters

    integer, allocatable, dimension(:) :: proc_coords !<
    !! Processor coordinates in MPI_CART_COMM

    integer, allocatable, dimension(:) :: start_idx !<
    !! Starting cell-center index of local processor in global grid

#ifdef MFC_MPI

    type(mpi_io_var), public :: MPI_IO_DATA
    type(mpi_io_ib_var), public :: MPI_IO_IB_DATA
    type(mpi_io_airfoil_ib_var), public :: MPI_IO_airfoil_IB_DATA
    type(mpi_io_levelset_var), public :: MPI_IO_levelset_DATA
    type(mpi_io_levelset_norm_var), public :: MPI_IO_levelsetnorm_DATA

    character(LEN=name_len) :: mpiiofs
    integer :: mpi_info_int !<
    !! MPI info for parallel IO with Lustre file systems

#endif

    ! Initial Condition Parameters
    integer :: num_patches     !< Number of patches composing initial condition

    type(ic_patch_parameters), dimension(num_patches_max) :: patch_icpp !<
    !! Database of the initial condition patch parameters (icpp) for each of the
    !! patches employed in the configuration of the initial condition. Note that
    !! the maximum allowable number of patches, num_patches_max, may be changed
    !! in the module m_derived_types.f90.

    integer :: num_bc_patches  !< Number of boundary condition patches
    logical :: bc_io !< whether or not to save BC data
    type(bc_patch_parameters), dimension(num_bc_patches_max) :: patch_bc
    !! Database of the boundary condition patch parameters for each of the patches
    !! employed in the configuration of the boundary conditions

    ! Fluids Physical Parameters
    type(physical_parameters), dimension(num_fluids_max) :: fluid_pp !<
    !! Database of the physical parameters of each of the fluids that is present
    !! in the flow. These include the stiffened gas equation of state parameters,
    !! the Reynolds numbers and the Weber numbers.

    real(wp) :: rhoref, pref !< Reference parameters for Tait EOS

    !> @name Bubble modeling
    !> @{
    integer :: nb
    real(wp) :: R0ref
    real(wp) :: Ca, Web, Re_inv
    real(wp), dimension(:), allocatable :: weight, R0
    logical :: bubbles_euler
    logical :: qbmm      !< Quadrature moment method
    integer :: nmom  !< Number of carried moments
    real(wp) :: sigR, sigV, rhoRV !< standard deviations in R/V
    logical :: adv_n !< Solve the number density equation and compute alpha from number density
    !> @}

    !> @name Immersed Boundaries
    !> @{
    logical :: ib           !< Turn immersed boundaries on
    integer :: num_ibs      !< Number of immersed boundaries
    integer :: Np

    type(ib_patch_parameters), dimension(num_patches_max) :: patch_ib

    type(vec3_dt), allocatable, dimension(:) :: airfoil_grid_u, airfoil_grid_l
    !! Database of the immersed boundary patch parameters for each of the
    !! patches employed in the configuration of the initial condition. Note that
    !! the maximum allowable number of patches, num_patches_max, may be changed
    !! in the module m_derived_types.f90.

    !> @}

    !> @name Non-polytropic bubble gas compression
    !> @{
    logical :: polytropic
    logical :: polydisperse
    integer :: thermal  !1 = adiabatic, 2 = isotherm, 3 = transfer
    real(wp) :: R_n, R_v, phi_vn, phi_nv, Pe_c, Tw, pv, M_n, M_v
    real(wp), dimension(:), allocatable :: k_n, k_v, pb0, mass_n0, mass_v0, Pe_T
    real(wp), dimension(:), allocatable :: Re_trans_T, Re_trans_c, Im_trans_T, Im_trans_c, omegaN
    real(wp) :: mul0, ss, gamma_v, mu_v
    real(wp) :: gamma_m, gamma_n, mu_n
    real(wp) :: poly_sigma
    integer :: dist_type !1 = binormal, 2 = lognormal-normal
    !> @}

    !> @name Surface Tension Modeling
    !> @{
    real(wp) :: sigma
    logical :: surface_tension
    !> @}

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

    integer, allocatable, dimension(:, :, :) :: logic_grid

    type(pres_field) :: pb
    type(pres_field) :: mv

    real(wp) :: Bx0 !< Constant magnetic field in the x-direction (1D)

    integer :: buff_size !<
    !! The number of cells that are necessary to be able to store enough boundary
    !! conditions data to march the solution in the physical computational domain
    !! to the next time-step.

contains

    !>  Assigns default values to user inputs prior to reading
        !!              them in. This allows for an easier consistency check of
        !!              these parameters once they are read from the input file.
    impure subroutine s_assign_default_values_to_user_inputs

        integer :: i !< Generic loop operator

        ! Logistics
        case_dir = '.'
        old_grid = .false.
        old_ic = .false.
        t_step_old = dflt_int
        t_step_start = dflt_int

        cfl_adap_dt = .false.
        cfl_const_dt = .false.
        cfl_dt = .false.
        n_start = dflt_int

        ! Computational domain parameters
        m = dflt_int; n = 0; p = 0

        call s_update_cell_bounds(cells_bounds, m, n, p)

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
        recon_type = WENO_TYPE
        weno_order = dflt_int
        igr = .false.
        igr_order = dflt_int
        muscl_order = dflt_int

        hypoelasticity = .false.
        hyperelasticity = .false.
        elasticity = .false.
        pre_stress = .false.
        b_size = dflt_int
        tensor_size = dflt_int
        cont_damage = .false.

        mhd = .false.
        relativity = .false.

        bc_x%beg = dflt_int; bc_x%end = dflt_int
        bc_y%beg = dflt_int; bc_y%end = dflt_int
        bc_z%beg = dflt_int; bc_z%end = dflt_int

        #:for DIM in ['x', 'y', 'z']
            #:for DIR in [1, 2, 3]
                bc_${DIM}$%vb${DIR}$ = 0._wp
                bc_${DIM}$%ve${DIR}$ = 0._wp
            #:endfor
        #:endfor

        parallel_io = .false.
        file_per_process = .false.
        precision = 2
        viscous = .false.
        bubbles_lagrange = .false.
        mixlayer_vel_profile = .false.
        mixlayer_vel_coef = 1._wp
        mixlayer_perturb = .false.
        mixlayer_perturb_nk = 100
        mixlayer_perturb_k0 = 0.4446_wp
        perturb_flow = .false.
        perturb_flow_fluid = dflt_int
        perturb_flow_mag = dflt_real
        perturb_sph = .false.
        perturb_sph_fluid = dflt_int
        fluid_rho = dflt_real
        elliptic_smoothing_iters = dflt_int
        elliptic_smoothing = .false.

        ! Initial condition parameters
        num_patches = dflt_int

        do i = 1, num_patches_max
            patch_icpp(i)%geometry = dflt_int
            patch_icpp(i)%model_scale(:) = 1._wp
            patch_icpp(i)%model_translate(:) = 0._wp
            patch_icpp(i)%model_filepath(:) = dflt_char
            patch_icpp(i)%model_spc = num_ray
            patch_icpp(i)%model_threshold = ray_tracing_threshold
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
            patch_icpp(i)%cv = 0._wp
            patch_icpp(i)%qv = 0._wp
            patch_icpp(i)%qvp = 0._wp
            patch_icpp(i)%tau_e = 0._wp
            patch_icpp(i)%Bx = dflt_real
            patch_icpp(i)%By = dflt_real
            patch_icpp(i)%Bz = dflt_real
            patch_icpp(i)%a(2) = dflt_real
            patch_icpp(i)%a(3) = dflt_real
            patch_icpp(i)%a(4) = dflt_real
            patch_icpp(i)%a(5) = dflt_real
            patch_icpp(i)%a(6) = dflt_real
            patch_icpp(i)%a(7) = dflt_real
            patch_icpp(i)%a(8) = dflt_real
            patch_icpp(i)%a(9) = dflt_real
            patch_icpp(i)%non_axis_sym = .false.

            !should get all of r0's and v0's
            patch_icpp(i)%r0 = dflt_real
            patch_icpp(i)%v0 = dflt_real

            patch_icpp(i)%p0 = dflt_real
            patch_icpp(i)%m0 = dflt_real

            patch_icpp(i)%hcid = dflt_int

            if (chemistry) then
                patch_icpp(i)%Y(:) = 0._wp
            end if
        end do

        num_bc_patches = 0
        bc_io = .false.

        do i = 1, num_bc_patches_max
            patch_bc(i)%geometry = dflt_int
            patch_bc(i)%type = dflt_int
            patch_bc(i)%dir = dflt_int
            patch_bc(i)%loc = dflt_int
            patch_bc(i)%centroid(:) = dflt_real
            patch_bc(i)%length(:) = dflt_real
            patch_bc(i)%radius = dflt_real
        end do

        ! Tait EOS
        rhoref = dflt_real
        pref = dflt_real

        ! Bubble modeling
        bubbles_euler = .false.
        polytropic = .true.
        polydisperse = .false.

        thermal = dflt_int
        R0ref = dflt_real
        nb = dflt_int

        Ca = dflt_real
        Re_inv = dflt_real
        Web = dflt_real
        poly_sigma = dflt_real
        surface_tension = .false.

        adv_n = .false.

        qbmm = .false.
        nmom = 1
        sigR = dflt_real
        sigV = dflt_real
        rhoRV = 0._wp
        dist_type = dflt_int

        R_n = dflt_real
        R_v = dflt_real
        phi_vn = dflt_real
        phi_nv = dflt_real
        Pe_c = dflt_real
        Tw = dflt_real

        ! surface tension modeling
        sigma = dflt_real
        pi_fac = 1._wp

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

            ! Proper default values for translating STL models
            patch_ib(i)%model_scale(:) = 1._wp
            patch_ib(i)%model_translate(:) = 0._wp
            patch_ib(i)%model_rotate(:) = 0._wp
            patch_ib(i)%model_filepath(:) = dflt_char
            patch_ib(i)%model_spc = num_ray
            patch_ib(i)%model_threshold = ray_tracing_threshold
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
            fluid_pp(i)%cv = 0._wp
            fluid_pp(i)%qv = 0._wp
            fluid_pp(i)%qvp = 0._wp
            fluid_pp(i)%G = 0._wp
        end do

        Bx0 = dflt_real

    end subroutine s_assign_default_values_to_user_inputs

    !> Computation of parameters, allocation procedures, and/or
        !! any other tasks needed to properly setup the module
    impure subroutine s_initialize_global_parameters_module

        integer :: i, j, fac

        if (recon_type == WENO_TYPE) then
            weno_polyn = (weno_order - 1)/2
        elseif (recon_type == MUSCL_TYPE) then
            muscl_polyn = muscl_order
        end if

        ! Determining the layout of the state vectors and overall size of
        ! the system of equations, given the dimensionality and choice of
        ! the equations of motion

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

            if (bubbles_euler) then
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

                allocate (weight(nb), R0(nb))
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
                    weight(:) = 1._wp
                    R0(:) = 1._wp
                else if (nb < 1) then
                    stop 'Invalid value of nb'
                end if

                !Initialize pref,rhoref for polytropic qbmm (done in s_initialize_nonpoly for non-polytropic)
                if (.not. qbmm) then
                    if (polytropic) then
                        rhoref = 1._wp
                        pref = 1._wp
                    end if
                end if

                !Initialize pb0,pv,pref,rhoref for polytropic qbmm (done in s_initialize_nonpoly for non-polytropic)
                if (qbmm) then
                    if (polytropic) then
                        allocate (pb0(nb))
                        if ((f_is_default(Web))) then
                            pb0 = pref
                            pb0 = pb0/pref
                            pref = 1._wp
                        end if
                        rhoref = 1._wp
                    end if
                end if
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

        else if (model_eqns == 4) then
            ! 4 equation model with subgrid bubbles_euler
            cont_idx%beg = 1 ! one continuity equation
            cont_idx%end = 1 ! num_fluids
            mom_idx%beg = cont_idx%end + 1 ! one momentum equation in each direction
            mom_idx%end = cont_idx%end + num_vels
            E_idx = mom_idx%end + 1 ! one energy equation
            adv_idx%beg = E_idx + 1
            adv_idx%end = adv_idx%beg !one volume advection equation
            alf_idx = adv_idx%end
            sys_size = alf_idx !adv_idx%end

            if (bubbles_euler) then
                bub_idx%beg = sys_size + 1
                bub_idx%end = sys_size + 2*nb
                if (.not. polytropic) then
                    bub_idx%end = sys_size + 4*nb
                end if
                sys_size = bub_idx%end

                allocate (bub_idx%rs(nb), bub_idx%vs(nb))
                allocate (bub_idx%ps(nb), bub_idx%ms(nb))
                allocate (weight(nb), R0(nb))

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
                ! number of entries in the symmetric btensor plus the jacobian
                b_size = (num_dims*(num_dims + 1))/2 + 1
                tensor_size = num_dims**2 + 1
                xi_idx%beg = sys_size + 1
                xi_idx%end = sys_size + num_dims
                ! adding three more equations for the \xi field and the elastic energy
                sys_size = xi_idx%end + 1
            end if

            if (surface_tension) then
                c_idx = sys_size + 1
                sys_size = c_idx
            end if

            if (cont_damage) then
                damage_idx = sys_size + 1
                sys_size = damage_idx
            end if

        end if

        if (chemistry) then
            species_idx%beg = sys_size + 1
            species_idx%end = sys_size + num_species
            sys_size = species_idx%end
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

        call s_configure_coordinate_bounds(recon_type, weno_polyn, muscl_polyn, buff_size, &
                                           idwint, idwbuff, viscous, &
                                           bubbles_lagrange, m, n, p, &
                                           num_dims, igr)

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

        if (ib) then
            allocate (MPI_IO_IB_DATA%var%sf(0:m, 0:n, 0:p))
            allocate (MPI_IO_levelset_DATA%var%sf(0:m, 0:n, 0:p, 1:num_ibs))
            allocate (MPI_IO_levelsetnorm_DATA%var%sf(0:m, 0:n, 0:p, 1:num_ibs, 1:3))
        end if
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

    end subroutine s_initialize_global_parameters_module

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

    impure subroutine s_finalize_global_parameters_module

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

    end subroutine s_finalize_global_parameters_module

end module m_global_parameters
