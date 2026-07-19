!>
!! @file
!! @brief Contains module m_global_parameters

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief Global parameters for the computational domain, fluid properties, and simulation algorithm configuration
module m_global_parameters

#ifdef MFC_MPI
    use mpi  !< Message passing interface (MPI) module
#endif

    use m_derived_types
    use m_helper_basic
    ! Shared state: generated_decls, generated_case_opt_decls, sys_size, eqn_idx, b_size, tensor_size, chemistry, elasticity,
    ! shear_*
    use m_global_parameters_common
    ! $:USE_GPU_MODULE()

    implicit none

    real(wp) :: wall_time = 0
    real(wp) :: wall_time_avg = 0

    ! Logistics
    integer :: num_procs  !< Number of processors
    ! Computational Domain Parameters
    integer :: proc_rank  !< Rank of the local processor
    $:GPU_DECLARE(create='[num_procs, proc_rank]')

    !> @name Max and min number of cells in a direction of each combination of x-,y-, and z-
    type(cell_num_bounds) :: cells_bounds

    !> @name Global number of cells in each direction
    !> @{
    integer :: m_glb, n_glb, p_glb
    !> @}

    !> @name Cylindrical coordinates (either axisymmetric or full 3D)
    !> @{
    integer :: grid_geometry
    !> @}
    $:GPU_DECLARE(create='[grid_geometry]')

    !> @name Cell-boundary (CB) locations in the x-, y- and z-directions, respectively
    !> @{
    real(wp), target, allocatable, dimension(:) :: x_cb, y_cb, z_cb
    type(bounds_info), dimension(3)             :: glb_bounds
    !> @}

    !> @name Cell-center (CC) locations in the x-, y- and z-directions, respectively
    !> @{
    real(wp), target, allocatable, dimension(:) :: x_cc, y_cc, z_cc
    !> @}
    !> @name Cell-width distributions in the x-, y- and z-directions, respectively
    !> @{
    real(wp), target, allocatable, dimension(:) :: dx, dy, dz
    !> @}

    $:GPU_DECLARE(create='[x_cb, y_cb, z_cb, x_cc, y_cc, z_cc, dx, dy, dz]')

    ! dt, m, n, p, cfl_target: GPU-declared via generated_decls.fpp (registered params)
    $:GPU_DECLARE(create='[glb_bounds]')

    logical :: cfl_dt
    ! Simulation Algorithm Parameters generated_case_opt_decls.fpp: now in m_global_parameters_common

    integer :: hyper_model  !< hyperelasticity solver algorithm
    ! elasticity, chemistry: in m_global_parameters_common
    logical                :: shear_stress  !< Shear stresses
    logical                :: bulk_stress   !< Bulk stresses
    logical                :: bodyForces
    real(wp), dimension(3) :: accel_bf
    $:GPU_DECLARE(create='[accel_bf]')
    ! $:GPU_DECLARE(create='[k_x,w_x,p_x,g_x,k_y,w_y,p_y,g_y,k_z,w_z,p_z,g_z]')

    ! Synthetic turbulence (scalars auto-generated in generated_decls.fpp; their
    ! GPU_DECLARE lines live in m_global_parameters_common)
    integer, dimension(num_synth_shells_max)     :: synth_n_waves_per_shell
    real(wp), dimension(num_synth_shells_max)    :: synth_k_shell, synth_amp_shell
    real(wp), dimension(num_turb_sources_max, 3) :: turb_pos, synth_L
    $:GPU_DECLARE(create='[synth_n_waves_per_shell, synth_k_shell, synth_amp_shell]')
    $:GPU_DECLARE(create='[turb_pos, synth_L]')

    integer :: cpu_start, cpu_end, cpu_rate

    $:GPU_DECLARE(create='[hyper_model]')
    $:GPU_DECLARE(create='[shear_stress, bulk_stress]')

    logical               :: bc_io
    logical, dimension(3) :: periodic_bc
    !> @name Boundary conditions (BC) in the x-, y- and z-directions, respectively
    !> @{
    type(int_bounds_info) :: bc_x, bc_y, bc_z
    type(bc_xyz_info)     :: bc
    !> @}
    !> @name Original boundary conditions preserved for immersed boundary code
    !> (bc_x/y/z get overwritten with MPI neighbor ranks during decomposition)
    !> @{
    type(int_bounds_info) :: ib_bc_x, ib_bc_y, ib_bc_z
    !> @}
#if defined(MFC_OpenACC)
    $:GPU_DECLARE(create='[bc_x%vb1, bc_x%vb2, bc_x%vb3, bc_x%ve1, bc_x%ve2, bc_x%ve3]')
    $:GPU_DECLARE(create='[bc_y%vb1, bc_y%vb2, bc_y%vb3, bc_y%ve1, bc_y%ve2, bc_y%ve3]')
    $:GPU_DECLARE(create='[bc_z%vb1, bc_z%vb2, bc_z%vb3, bc_z%ve1, bc_z%ve2, bc_z%ve3]')
    $:GPU_DECLARE(create='[ib_bc_x%beg, ib_bc_y%beg, ib_bc_z%beg]')
#elif defined(MFC_OpenMP)
    $:GPU_DECLARE(create='[bc_x, bc_y, bc_z]')
    $:GPU_DECLARE(create='[ib_bc_x, ib_bc_y, ib_bc_z]')
#endif
    $:GPU_DECLARE(create='[bc]')
    type(bounds_info) :: neighbor_domain_x, neighbor_domain_y, neighbor_domain_z
    integer           :: num_gbl_ibs, num_local_ibs
    $:GPU_DECLARE(create='[neighbor_domain_x, neighbor_domain_y, neighbor_domain_z, num_gbl_ibs]')

    ! proc_coords, start_idx, mpiiofs, mpi_info_int: in m_global_parameters_common
    ! down_sample: GPU-declared via generated_decls.fpp (registered param)

    !> @name MPI domain-decomposition state for Lagrangian-bubble exchange (#1290)
    !> @{
    type(bounds_info), allocatable, dimension(:) :: pcomm_coords    !< Local rank physical domain bounds
    type(int_bounds_info), dimension(3)          :: nidx            !< Neighbor index offsets per direction
    integer, allocatable, dimension(:,:,:)       :: neighbor_ranks  !< MPI ranks of neighbors
    $:GPU_DECLARE(create='[pcomm_coords]')
    !> @}
    type(mpi_io_var), public                      :: MPI_IO_DATA
    type(mpi_io_ib_var), public                   :: MPI_IO_IB_DATA
    type(mpi_io_airfoil_ib_var), public           :: MPI_IO_airfoil_IB_DATA
    type(mpi_io_levelset_var), public             :: MPI_IO_levelset_DATA
    type(mpi_io_levelset_norm_var), public        :: MPI_IO_levelsetnorm_DATA
    real(wp), allocatable, dimension(:,:), public :: MPI_IO_DATA_lag_bubbles

    ! sys_size, eqn_idx, b_size, tensor_size: in m_global_parameters_common (GPU_DECLARE there too)
    type(qbmm_idx_info) :: qbmm_idx  !< QBMM moment index mappings (allocatable; GPU-managed separately).

    ! Cell Indices for the (local) interior points (O-m, O-n, 0-p). Stands for "InDices With INTerior".
    type(int_bounds_info) :: idwint(1:3)
    $:GPU_DECLARE(create='[idwint]')

    ! Cell Indices for the entire (local) domain. In simulation and post_process, this includes the buffer region. idwbuff and
    ! idwint are the same otherwise. Stands for "InDices With BUFFer".
    type(int_bounds_info) :: idwbuff(1:3)
    $:GPU_DECLARE(create='[idwbuff]')

    !> @name The number of fluids, along with their identifying indexes, respectively, for which viscous effects, e.g. the shear
    !! and/or the volume Reynolds (Re) numbers, will be non-negligible.
    !> @{
    integer, dimension(2)                :: Re_size
    integer                              :: Re_size_max
    integer, allocatable, dimension(:,:) :: Re_idx
    !> @}

    $:GPU_DECLARE(create='[Re_size, Re_size_max, Re_idx]')

    !> @name Herschel-Bulkley non-Newtonian viscosity: per-fluid flags and parameter arrays.
    !> @{
    logical                             :: any_non_newtonian  !< .true. if any fluid is non-Newtonian
    logical, allocatable, dimension(:)  :: is_non_newtonian   !< per-fluid NN flag
    real(wp), allocatable, dimension(:) :: hb_tau0, hb_K, hb_nn, hb_m_arr
    real(wp), allocatable, dimension(:) :: hb_mu_min, hb_mu_max
    real(wp), allocatable, dimension(:) :: fluid_inv_re       !< per-fluid Newtonian inverse-Re
    !> @}

    $:GPU_DECLARE(create='[any_non_newtonian, is_non_newtonian, hb_tau0, hb_K, hb_nn, hb_m_arr, hb_mu_min, hb_mu_max, fluid_inv_re]')

    ! WENO averaging flag: use arithmetic mean or unaltered WENO-reconstructed cell-boundary values
    !> @{
    real(wp) :: wa_flg
    !> @}

    $:GPU_DECLARE(create='[wa_flg]')

    !> @name The coordinate direction indexes and flags (flg), respectively, for which the configurations will be determined with
    !! respect to a working direction and that will be used to isolate the contributions, in that direction, in the dimensionally
    !! split system of equations.
    !> @{
    integer, dimension(3)  :: dir_idx
    real(wp), dimension(3) :: dir_flg
    integer, dimension(3)  :: dir_idx_tau  !< used for hypoelasticity=true
    !> @}

    $:GPU_DECLARE(create='[dir_idx, dir_flg, dir_idx_tau]')

    integer :: buff_size  !< Number of ghost cells for boundary condition storage
    $:GPU_DECLARE(create='[buff_size]')

    ! shear_num, shear_indices, shear_BC_flip_num, shear_BC_flip_indices: in m_global_parameters_common

    ! END: Simulation Algorithm Parameters

    ! Fluids Physical Parameters fluid_pp, bub_pp: auto-generated in generated_decls.fpp

    integer :: fd_number  !< Finite-difference half-stencil size: MAX(1, fd_order/2)
    $:GPU_DECLARE(create='[fd_number]')

    !> @name Centered finite-difference coefficients in x-, y- and z-coordinate directions
    !> @{
    real(wp), allocatable, dimension(:,:) :: fd_coeff_x
    real(wp), allocatable, dimension(:,:) :: fd_coeff_y
    real(wp), allocatable, dimension(:,:) :: fd_coeff_z
    !> @}
    $:GPU_DECLARE(create='[fd_coeff_x, fd_coeff_y, fd_coeff_z]')

    ! probe, integral: auto-generated in generated_decls.fpp

    !> @name Reference density and pressure for Tait EOS
    !> @{
    !> @name Immersed Boundaries
    !> patch_ib, ib_airfoil, stl_models, particle_cloud: auto-generated in generated_decls.fpp
    !> @{
    integer, dimension(num_local_ibs_max) :: local_ib_patch_ids  !< lookup table of IBs in the local compute domain
    integer, allocatable, dimension(:,:,:) :: ib_neighbor_ranks  !< MPI ranks of neighborhood domains, indexed (-N:N,-N:N,-N:N)
    type(ib_airfoil_grid), dimension(num_ib_airfoils_max) :: ib_airfoil_grids  !< Per-airfoil computed surface grids

    $:GPU_DECLARE(create='[ib_airfoil_grids]')
    !> @}

    !> @name Bubble modeling
    !> @{
    #:if MFC_CASE_OPTIMIZATION
        integer, parameter :: nb = ${nb}$  !< Number of eq. bubble sizes
    #:else
        integer :: nb
    #:endif

    real(wp) :: Eu  !< Euler number
    $:GPU_DECLARE(create='[Eu]')

    real(wp), dimension(:), allocatable :: weight  !< Simpson quadrature weights
    real(wp), dimension(:), allocatable :: R0      !< Bubble sizes
    $:GPU_DECLARE(create='[weight, R0]')

    real(wp), allocatable, dimension(:,:,:) :: ptil  !< Pressure modification
    $:GPU_DECLARE(create='[ptil]')

    integer, parameter :: nmom = 6  !< Number of carried moments per R0 location
    integer            :: nmomsp    !< Number of moments required by ensemble-averaging
    integer            :: nmomtot   !< Total number of carried moments moments/transport equations
    $:GPU_DECLARE(create='[nmomsp, nmomtot]')

    #:if not MFC_CASE_OPTIMIZATION
        $:GPU_DECLARE(create='[nb]')
    #:endif

    type(scalar_field), allocatable, dimension(:)     :: mom_sp
    type(scalar_field), allocatable, dimension(:,:,:) :: mom_3d
    $:GPU_DECLARE(create='[mom_sp, mom_3d]')
    !> @}

    ! chem_params: auto-generated in generated_decls.fpp

    !> @name Physical bubble parameters (see Ando 2010, Preston 2007)
    !> @{
    real(wp) :: phi_vg, phi_gv, Pe_c, Tw, k_vl, k_gl
    $:GPU_DECLARE(create='[phi_vg, phi_gv, Pe_c, Tw, k_vl, k_gl]')

    real(wp), dimension(:), allocatable :: pb0, mass_g0, mass_v0, Pe_T, k_v, k_g
    real(wp), dimension(:), allocatable :: Re_trans_T, Re_trans_c, Im_trans_T, Im_trans_c, omegaN
    $:GPU_DECLARE(create='[pb0, mass_g0, mass_v0, Pe_T, k_v, k_g]')
    $:GPU_DECLARE(create='[Re_trans_T, Re_trans_c, Im_trans_T, Im_trans_c, omegaN]')

    real(wp) :: gam, gam_m
    $:GPU_DECLARE(create='[gam, gam_m]')

    real(wp) :: p0ref, rho0ref, T0ref, ss, pv, vd, mu_l, mu_v, mu_g, gam_v, gam_g, M_v, M_g, cp_v, cp_g, R_v, R_g
    $:GPU_DECLARE(create='[p0ref, rho0ref, T0ref, ss, pv, vd, mu_l, mu_v, mu_g, gam_v, gam_g, M_v, M_g, cp_v, cp_g, R_v, R_g]')
    !> @}

    ! acoustic: auto-generated in generated_decls.fpp

    !> @name Surface tension parameters
    !> @{
    !> @}

    real(wp), allocatable, dimension(:) :: gammas, gs_min, pi_infs, ps_inf, cvs, qvs, qvps
    $:GPU_DECLARE(create='[gammas, gs_min, pi_infs, ps_inf, cvs, qvs, qvps]')

    real(wp), allocatable, dimension(:) :: jwl_As, jwl_Bs, jwl_R1s, jwl_R2s, jwl_omegas, jwl_rho0s, jwl_E0s
    real(wp), allocatable, dimension(:) :: jwl_air_e0s, jwl_air_rho0s, jwl_air_gammas, jwl_air_pi_infs
    $:GPU_DECLARE(create='[jwl_As, jwl_Bs, jwl_R1s, jwl_R2s, jwl_omegas, jwl_rho0s, jwl_E0s]')
    $:GPU_DECLARE(create='[jwl_air_e0s, jwl_air_rho0s, jwl_air_gammas, jwl_air_pi_infs]')

    real(wp)                                    :: mytime     !< Current simulation time
    real(wp)                                    :: finaltime  !< Final simulation time
    type(pres_field), allocatable, dimension(:) :: pb_ts
    type(pres_field), allocatable, dimension(:) :: mv_ts

    $:GPU_DECLARE(create='[mytime, pb_ts, mv_ts]')

    !> @name lagrangian subgrid bubble parameters
    !> lag_params: auto-generated in generated_decls.fpp
    !> @{!
    ! lag_params (decl + GPU_DECLARE) auto-generated in generated_decls.fpp; bubbles_lagrange GPU-declared in
    ! m_global_parameters_common
    integer :: n_el_bubs_loc, n_el_bubs_glb  !< Number of Lagrangian bubbles (local and global)
    logical :: moving_lag_bubbles
    logical :: lag_pressure_force
    logical :: lag_gravity_force
    integer :: lag_vel_model, lag_drag_model
    $:GPU_DECLARE(create='[n_el_bubs_loc, n_el_bubs_glb]')
    $:GPU_DECLARE(create='[moving_lag_bubbles, lag_vel_model, lag_drag_model]')
    $:GPU_DECLARE(create='[lag_pressure_force, lag_gravity_force]')
    !> @}

    !> @name Continuum damage model parameters
    !> @{!
    !> @}

    !> @name MHD Hyperbolic cleaning parameters
    !> @{!
    !> @}

contains

    !> Assigns default values to the user inputs before reading them in. This enables for an easier consistency check of these
    !! parameters once they are read from the input file.
    impure subroutine s_assign_default_values_to_user_inputs

        integer :: i, j  !< Generic loop iterator

        ! Shared defaults (case_dir, m/n/p, cyl_coord, cfl flags, model_eqns, elasticity, BC blocks,
        ! recon/weno/muscl/num_fluids/igr/mhd/relativity under case-opt guard, Tait EOS, bubble flags,
        ! IB flags, parallel I/O flags, fft_wrt)

        call s_assign_common_defaults

        ! Boundary conditions (bc_x/y/z are per-target declarations, not visible in common)
        bc_x%beg = dflt_int; bc_x%end = dflt_int
        bc_y%beg = dflt_int; bc_y%end = dflt_int
        bc_z%beg = dflt_int; bc_z%end = dflt_int

        #:for DIM in ['x', 'y', 'z']
            #:for DIR in [1, 2, 3]
                bc_${DIM}$%vb${DIR}$ = 0._wp
                bc_${DIM}$%ve${DIR}$ = 0._wp
            #:endfor
        #:endfor

        #:for dir in ['x', 'y', 'z']
            bc_${dir}$%isothermal_in = .false.
            bc_${dir}$%isothermal_out = .false.
            bc_${dir}$%Twall_in = dflt_real
            bc_${dir}$%Twall_out = dflt_real
        #:endfor

        call s_update_cell_bounds(cells_bounds, m, n, p)

        ! Logistics (sim-specific)
        run_time_info = .false.
        t_step_old = dflt_int

        ! Computational domain parameters (sim-specific)
        dt = dflt_real
        cfl_dt = .false.
        cfl_target = dflt_real

        t_step_stop = dflt_int
        t_step_save = dflt_int
        t_step_print = 1

        t_stop = dflt_real
        t_save = dflt_real

        ! NVIDIA UVM options
        nv_uvm_out_of_core = .false.
        nv_uvm_igr_temps_on_gpu = 3  ! => jac, jac_rhs, and jac_old on GPU (default)
        nv_uvm_pref_gpu = .false.

        ! Simulation algorithm parameters (sim-specific)
        mpp_lim = .false.
        time_stepper = dflt_int
        muscl_eps = dflt_real
        weno_eps = dflt_real
        teno_CT = dflt_real
        mp_weno = .false.
        weno_avg = .false.
        weno_Re_flux = .false.
        riemann_solver = dflt_int
        low_Mach = 0
        wave_speeds = dflt_int
        avg_state = dflt_int
        alt_soundspeed = .false.
        null_weights = .false.
        mixture_err = .false.
        precision = 2
        palpha_eps = dflt_real
        ptgalpha_eps = dflt_real
        int_comp = 0
        ic_eps = dflt_ic_eps
        ic_beta = dflt_ic_beta
        hyper_model = dflt_int
        rdma_mpi = .false.
        shear_stress = .false.
        bulk_stress = .false.
        any_non_newtonian = .false.
        num_igr_iters = dflt_num_igr_iters
        num_igr_warm_start_iters = dflt_num_igr_warm_start_iters
        alf_factor = dflt_alf_factor

        #:if not MFC_CASE_OPTIMIZATION
            mapped_weno = .false.
            wenoz = .false.
            teno = .false.
            wenoz_q = dflt_real
            igr_order = dflt_int
            igr_pres_lim = .false.
            viscous = .false.
            igr_iter_solver = 1
        #:endif

        chem_params%diffusion = .false.
        chem_params%reactions = .false.
        chem_params%gamma_method = 1
        chem_params%transport_model = 1

        chem_params%reaction_substeps = 0
        chem_params%adap_substeps = .false.
        chem_params%reaction_substeps_max = 0

        num_bc_patches = 0
        bc_io = .false.
        periodic_bc = .false.

        ! bc_x/y/z (incl. vb/ve loop) already defaulted above; glb_bounds is #1290's grid-derived global extent
        glb_bounds(1)%beg = dflt_real; glb_bounds(1)%end = dflt_real
        glb_bounds(2)%beg = dflt_real; glb_bounds(2)%end = dflt_real
        glb_bounds(3)%beg = dflt_real; glb_bounds(3)%end = dflt_real

        ! Fluids physical parameters (sim-specific; Re(:) and G=0._wp differ from post)
        do i = 1, num_fluids_max
            fluid_pp(i)%gamma = dflt_real
            fluid_pp(i)%pi_inf = dflt_real
            fluid_pp(i)%cv = 0._wp
            fluid_pp(i)%qv = 0._wp
            fluid_pp(i)%qvp = 0._wp
            fluid_pp(i)%Re(:) = dflt_real
            fluid_pp(i)%G = 0._wp
            fluid_pp(i)%non_newtonian = .false.
            fluid_pp(i)%K = dflt_real
            fluid_pp(i)%nn = dflt_real
            fluid_pp(i)%tau0 = 0._wp
            fluid_pp(i)%hb_m = dflt_real
            fluid_pp(i)%mu_min = dflt_real
            fluid_pp(i)%mu_max = dflt_real
            fluid_pp(i)%mu_bulk = dflt_real
            fluid_pp(i)%eos = eos_stiffened_gas
            call s_assign_jwl_fluid_defaults(fluid_pp(i))
        end do

        ! Subgrid bubble parameters (bub_pp struct + scalar companions; scalar companions are
        ! per-target manual declarations not in m_global_parameters_common scope)
        bub_pp%R0ref = dflt_real; R0ref = dflt_real
        bub_pp%p0ref = dflt_real; p0ref = dflt_real
        bub_pp%rho0ref = dflt_real; rho0ref = dflt_real
        bub_pp%T0ref = dflt_real; T0ref = dflt_real
        bub_pp%ss = dflt_real; ss = dflt_real
        bub_pp%pv = dflt_real; pv = dflt_real
        bub_pp%vd = dflt_real; vd = dflt_real
        bub_pp%mu_l = dflt_real; mu_l = dflt_real
        bub_pp%mu_v = dflt_real; mu_v = dflt_real
        bub_pp%mu_g = dflt_real; mu_g = dflt_real
        bub_pp%gam_v = dflt_real; gam_v = dflt_real
        bub_pp%gam_g = dflt_real; gam_g = dflt_real
        bub_pp%M_v = dflt_real; M_v = dflt_real
        bub_pp%M_g = dflt_real; M_g = dflt_real
        bub_pp%k_v = dflt_real
        bub_pp%k_g = dflt_real
        bub_pp%cp_v = dflt_real; cp_v = dflt_real
        bub_pp%cp_g = dflt_real; cp_g = dflt_real
        bub_pp%R_v = dflt_real; R_v = dflt_real
        bub_pp%R_g = dflt_real; R_g = dflt_real

        ! Immersed Boundaries (sim-specific extras)
        ib_neighborhood_radius = 1
        collision_model = 0
        coefficient_of_restitution = dflt_real
        collision_time = dflt_real
        ib_coefficient_of_friction = dflt_real
        ib_state_wrt = .false.
        many_ib_patch_parallelism = .false.

        ! Bubble modeling (sim-specific)
        bubble_model = 1
        polytropic = .true.
        thermal = dflt_int

        #:if not MFC_CASE_OPTIMIZATION
            nb = 1
            muscl_lim = dflt_int
        #:endif

        adv_n = .false.
        adap_dt = .false.
        adap_dt_tol = dflt_adap_dt_tol
        adap_dt_max_iters = dflt_adap_dt_max_iters

        pi_fac = 1._wp

        Eu = dflt_real
        Ca = dflt_real
        Re_inv = dflt_real
        Web = dflt_real

        ! Acoustic source
        acoustic_source = .false.
        num_source = dflt_int

        bodyForces = .false.
        bf_x = .false.; bf_y = .false.; bf_z = .false.
        !> amplitude, frequency, and phase shift sinusoid in each direction
        #:for dir in ['x', 'y', 'z']
            #:for param in ['k', 'w', 'p', 'g']
                ${param}$_${dir}$ = dflt_real
            #:endfor
        #:endfor

        synthetic_turbulence = .false.
        synth_seed = 1234
        synth_n_shells = dflt_int
        num_turbulent_sources = 0
        synth_U_inf = dflt_real
        synth_n_waves_per_shell = 0
        synth_k_shell = dflt_real
        synth_amp_shell = dflt_real
        turb_pos = dflt_real
        synth_L = dflt_real

        do j = 1, num_probes_max
            acoustic(j)%pulse = dflt_int
            acoustic(j)%support = dflt_int
            acoustic(j)%dipole = .false.
            do i = 1, 3
                acoustic(j)%loc(i) = dflt_real
            end do
            acoustic(j)%mag = dflt_real
            acoustic(j)%length = dflt_real
            acoustic(j)%height = dflt_real
            acoustic(j)%wavelength = dflt_real
            acoustic(j)%frequency = dflt_real
            acoustic(j)%gauss_sigma_dist = dflt_real
            acoustic(j)%gauss_sigma_time = dflt_real
            acoustic(j)%npulse = dflt_real
            acoustic(j)%dir = dflt_real
            acoustic(j)%delay = dflt_real
            acoustic(j)%foc_length = dflt_real
            acoustic(j)%aperture = dflt_real
            acoustic(j)%element_spacing_angle = dflt_real
            acoustic(j)%element_polygon_ratio = dflt_real
            acoustic(j)%rotate_angle = dflt_real
            acoustic(j)%num_elements = dflt_int
            acoustic(j)%element_on = dflt_int
            acoustic(j)%bb_num_freq = dflt_int
            acoustic(j)%bb_lowest_freq = dflt_real
            acoustic(j)%bb_bandwidth = dflt_real
        end do

        fd_order = dflt_int
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
            integral(i)%zmin = dflt_real
            integral(i)%zmax = dflt_real
        end do

        ! GRCBC flags
        #:for dir in ['x', 'y', 'z']
            bc_${dir}$%grcbc_in = .false.
            bc_${dir}$%grcbc_out = .false.
            bc_${dir}$%grcbc_vel_out = .false.
        #:endfor

        ! Lagrangian subgrid bubble model
        lag_params%solver_approach = dflt_int
        lag_params%cluster_type = dflt_int
        lag_params%pressure_corrector = .false.
        lag_params%smooth_type = dflt_int
        lag_params%heatTransfer_model = .false.
        lag_params%massTransfer_model = .false.
        lag_params%write_bubbles = .false.
        lag_params%write_bubbles_stats = .false.
        lag_params%write_void_evol = .false.
        lag_params%nBubs_glb = dflt_int
        lag_params%vel_model = dflt_int
        lag_params%drag_model = dflt_int
        lag_params%pressure_force = .true.
        lag_params%gravity_force = .false.
        lag_params%kahan_summation = .true.
        lag_params%epsilonb = 1._wp
        lag_params%charwidth = dflt_real
        lag_params%charNz = dflt_int
        lag_params%valmaxvoid = dflt_real
        lag_params%input_path = 'input/lag_bubbles.dat'
        moving_lag_bubbles = .false.
        lag_vel_model = dflt_int

        ! Continuum damage model
        tau_star = dflt_real
        cont_damage_s = dflt_real
        alpha_bar = dflt_real

        ! MHD (sim-specific extras beyond common Bx0)
        hyper_cleaning_speed = dflt_real
        hyper_cleaning_tau = dflt_real

        do i = 1, num_ib_airfoils_max
            ib_airfoil(i)%c = dflt_real
            ib_airfoil(i)%p = dflt_real
            ib_airfoil(i)%t = dflt_real
            ib_airfoil(i)%m = dflt_real
            ib_airfoil_grids(i)%Np = 0
        end do

        num_particle_clouds = 0
        do i = 1, num_particle_clouds_max
            particle_cloud(i)%x_centroid = 0._wp
            particle_cloud(i)%y_centroid = 0._wp
            particle_cloud(i)%z_centroid = 0._wp
            particle_cloud(i)%length_x = dflt_real
            particle_cloud(i)%length_y = dflt_real
            particle_cloud(i)%length_z = dflt_real
            particle_cloud(i)%num_particles = 0
            particle_cloud(i)%radius = dflt_real
            particle_cloud(i)%mass = dflt_real
            particle_cloud(i)%min_spacing = 0._wp
            particle_cloud(i)%moving_ibm = 0
            particle_cloud(i)%seed = 0
            particle_cloud(i)%packing_method = dflt_int
        end do

        do i = 1, num_ib_patches_max_namelist
            patch_ib(i)%gbl_patch_id = i
            patch_ib(i)%geometry = dflt_int
            patch_ib(i)%x_centroid = 0._wp
            patch_ib(i)%y_centroid = 0._wp
            patch_ib(i)%z_centroid = 0._wp
            patch_ib(i)%length_x = dflt_real
            patch_ib(i)%length_y = dflt_real
            patch_ib(i)%length_z = dflt_real
            patch_ib(i)%radius = dflt_real
            patch_ib(i)%airfoil_id = 0
            patch_ib(i)%model_id = 0
            patch_ib(i)%slip = .false.

            ! Variables to handle moving immersed boundaries, defaulting to no movement
            patch_ib(i)%moving_ibm = 0
            patch_ib(i)%vel(:) = 0._wp
            patch_ib(i)%angles(:) = 0._wp
            patch_ib(i)%angular_vel(:) = 0._wp
            patch_ib(i)%mass = dflt_real
            patch_ib(i)%moment = dflt_real
            patch_ib(i)%centroid_offset(:) = 0._wp

            ! sets values of a rotation matrix which can be used when calculating rotations
            patch_ib(i)%rotation_matrix = 0._wp
            patch_ib(i)%rotation_matrix(1, 1) = 1._wp
            patch_ib(i)%rotation_matrix(2, 2) = 1._wp
            patch_ib(i)%rotation_matrix(3, 3) = 1._wp
            patch_ib(i)%rotation_matrix_inverse = patch_ib(i)%rotation_matrix
        end do

        num_stl_models = 0

        do i = 1, num_stl_models_max
            stl_models(i)%model_filepath(:) = dflt_char
            stl_models(i)%model_translate(:) = 0._wp
            stl_models(i)%model_scale(:) = 1._wp
            stl_models(i)%model_threshold = ray_tracing_threshold
        end do

    end subroutine s_assign_default_values_to_user_inputs

    !> Initialize the global parameters module
    impure subroutine s_initialize_global_parameters_module

        integer :: i, j, k
        integer :: fac

        #:if not MFC_CASE_OPTIMIZATION
            ! Determining the degree of the WENO polynomials
            if (recon_type == recon_type_weno) then
                weno_polyn = (weno_order - 1)/2
                if (teno) then
                    weno_num_stencils = weno_order - 3
                else
                    weno_num_stencils = weno_polyn
                end if
            else if (recon_type == recon_type_muscl) then
                muscl_polyn = muscl_order
            end if
            $:GPU_UPDATE(device='[weno_polyn, muscl_polyn]')
            $:GPU_UPDATE(device='[weno_num_stencils]')
            $:GPU_UPDATE(device='[nb]')
            $:GPU_UPDATE(device='[num_dims, num_vels, num_fluids]')
            $:GPU_UPDATE(device='[igr, igr_order, igr_iter_solver]')
        #:endif

        ! muscl_eps: use per-limiter defaults when user did not set it
        if (f_is_default(muscl_eps)) then
            if (muscl_lim <= 2) then
                muscl_eps = 1e-9_wp  ! minmod, MC
            else
                muscl_eps = 1e-6_wp  ! Van Albada, Van Leer, SUPERBEE
            end if
        end if

        ! Initialize counts: viscous fluids, surface-tension interfaces, curvature interfaces
        Re_size = 0
        Re_size_max = 0

        ! Populate eqn_idx, sys_size, b_size, tensor_size, elasticity, shear_* (shared logic)
        call s_initialize_eqn_idx(nmom, nb)

        ! sim-only: GPU update for shear state after s_initialize_eqn_idx populated it
        if (model_eqns == model_eqns_5eq .or. model_eqns == model_eqns_6eq) then
            if (hypoelasticity .or. hyperelasticity) then
                $:GPU_UPDATE(device='[shear_num, shear_indices, shear_BC_flip_num, shear_BC_flip_indices]')
            end if
        end if

        ! Per-target (sim): nmomsp/nmomtot for qbmm, qbmm_idx alloc/fill, gam, Re_idx
        if (model_eqns == model_eqns_5eq .and. bubbles_euler) then
            if (qbmm) then
                nmomsp = 4  ! number of special moments
                if (nnode == 4) nmomtot = nmom*nb
            end if

            @:ALLOCATE(qbmm_idx%rs(nb), qbmm_idx%vs(nb))
            @:ALLOCATE(qbmm_idx%ps(nb), qbmm_idx%ms(nb))

            gam = bub_pp%gam_g

            if (qbmm) then
                @:ALLOCATE(qbmm_idx%moms(nb, nmom))
                do i = 1, nb
                    do j = 1, nmom
                        qbmm_idx%moms(i, j) = eqn_idx%bub%beg + (j - 1) + (i - 1)*nmom
                    end do
                    qbmm_idx%rs(i) = qbmm_idx%moms(i, 2)
                    qbmm_idx%vs(i) = qbmm_idx%moms(i, 3)
                end do
            else
                do i = 1, nb
                    if (.not. polytropic) then
                        fac = 4
                    else
                        fac = 2
                    end if

                    qbmm_idx%rs(i) = eqn_idx%bub%beg + (i - 1)*fac
                    qbmm_idx%vs(i) = qbmm_idx%rs(i) + 1

                    if (.not. polytropic) then
                        qbmm_idx%ps(i) = qbmm_idx%vs(i) + 1
                        qbmm_idx%ms(i) = qbmm_idx%ps(i) + 1
                    end if
                end do
            end if
        end if

        if (model_eqns == model_eqns_4eq .and. bubbles_euler) then
            @:ALLOCATE(qbmm_idx%rs(nb), qbmm_idx%vs(nb))
            @:ALLOCATE(qbmm_idx%ps(nb), qbmm_idx%ms(nb))

            do i = 1, nb
                if (polytropic) then
                    fac = 2
                else
                    fac = 4
                end if

                qbmm_idx%rs(i) = eqn_idx%bub%beg + (i - 1)*fac
                qbmm_idx%vs(i) = qbmm_idx%rs(i) + 1

                if (.not. polytropic) then
                    qbmm_idx%ps(i) = qbmm_idx%vs(i) + 1
                    qbmm_idx%ms(i) = qbmm_idx%ps(i) + 1
                end if
            end do
        end if

        ! sim-only: Re_idx (non-gamma-law models only)
        if (model_eqns /= model_eqns_gamma_law) then
            ! Count fluids with non-negligible viscous effects (Re > 0)
            do i = 1, num_fluids
                if (fluid_pp(i)%Re(1) > 0) Re_size(1) = Re_size(1) + 1
                if (fluid_pp(i)%Re(2) > 0) Re_size(2) = Re_size(2) + 1
            end do

            if (Re_size(1) > 0._wp) shear_stress = .true.
            if (Re_size(2) > 0._wp) bulk_stress = .true.

            Re_size_max = maxval(Re_size)

            $:GPU_UPDATE(device='[Re_size, Re_size_max, shear_stress, bulk_stress]')

            ! Bookkeeping the indexes of any viscous fluids
            if (viscous) then
                @:ALLOCATE(Re_idx(1:2, 1:Re_size_max))

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

        ! Herschel-Bulkley non-Newtonian viscosity: gather per-fluid parameters into device arrays
        @:ALLOCATE(is_non_newtonian(1:num_fluids))
        @:ALLOCATE(hb_tau0(1:num_fluids), hb_K(1:num_fluids), hb_nn(1:num_fluids), hb_m_arr(1:num_fluids))
        @:ALLOCATE(hb_mu_min(1:num_fluids), hb_mu_max(1:num_fluids))
        @:ALLOCATE(fluid_inv_re(1:num_fluids))

        any_non_newtonian = .false.
        do i = 1, num_fluids
            is_non_newtonian(i) = fluid_pp(i)%non_newtonian
            if (is_non_newtonian(i)) any_non_newtonian = .true.
            hb_tau0(i) = fluid_pp(i)%tau0
            hb_K(i) = fluid_pp(i)%K
            hb_nn(i) = fluid_pp(i)%nn
            hb_m_arr(i) = fluid_pp(i)%hb_m
            hb_mu_min(i) = fluid_pp(i)%mu_min
            hb_mu_max(i) = fluid_pp(i)%mu_max
            if (fluid_pp(i)%Re(1) > 0._wp) then
                fluid_inv_re(i) = 1._wp/fluid_pp(i)%Re(1)
            else
                fluid_inv_re(i) = 0._wp
            end if
        end do
        $:GPU_UPDATE(device='[any_non_newtonian, is_non_newtonian, hb_tau0, hb_K, hb_nn, hb_m_arr, hb_mu_min, hb_mu_max, fluid_inv_re]')

        if (bubbles_euler .and. qbmm .and. .not. polytropic) then
            allocate (MPI_IO_DATA%view(1:sys_size + 2*nb*nnode))
            allocate (MPI_IO_DATA%var(1:sys_size + 2*nb*nnode))
        else if (bubbles_lagrange) then
            allocate (MPI_IO_DATA%view(1:sys_size + 1))
            allocate (MPI_IO_DATA%var(1:sys_size + 1))
        else
            allocate (MPI_IO_DATA%view(1:sys_size))
            allocate (MPI_IO_DATA%var(1:sys_size))
        end if

        if (.not. down_sample) then
            do i = 1, sys_size
                allocate (MPI_IO_DATA%var(i)%sf(0:m,0:n,0:p))
                MPI_IO_DATA%var(i)%sf => null()
            end do
        end if
        if (bubbles_euler .and. qbmm .and. .not. polytropic) then
            do i = sys_size + 1, sys_size + 2*nb*nnode
                allocate (MPI_IO_DATA%var(i)%sf(0:m,0:n,0:p))
                MPI_IO_DATA%var(i)%sf => null()
            end do
        else if (bubbles_lagrange) then
            do i = 1, sys_size + 1
                allocate (MPI_IO_DATA%var(i)%sf(0:m,0:n,0:p))
                MPI_IO_DATA%var(i)%sf => null()
            end do
        end if

        ! Configure WENO averaging flag (arithmetic mean vs. unaltered values)
        wa_flg = 0._wp; if (weno_avg) wa_flg = 1._wp
        $:GPU_UPDATE(device='[wa_flg]')

        ! Resort to default WENO-JS if no other WENO scheme is selected
        #:if not MFC_CASE_OPTIMIZATION
            wenojs = .not. (mapped_weno .or. wenoz .or. teno)
        #:endif

        if (ib) allocate (MPI_IO_IB_DATA%var%sf(0:m,0:n,0:p))

        if (elasticity .or. mhd .or. probe_wrt .or. ib .or. bubbles_lagrange) then
            fd_number = max(1, fd_order/2)
        end if

        call s_configure_coordinate_bounds(recon_type, weno_polyn, muscl_polyn, igr_order, buff_size, idwint, idwbuff, viscous, &
                                           & bubbles_lagrange, m, n, p, num_dims, igr, ib, fd_number)
        $:GPU_UPDATE(device='[idwint, idwbuff]')

        ! Configuring Coordinate Direction Indexes
        if (bubbles_euler) then
            @:ALLOCATE(ptil( idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end))
        end if

        $:GPU_UPDATE(device='[fd_order, fd_number]')

        if (cyl_coord .neqv. .true.) then  ! Cartesian grid
            grid_geometry = 1
        else if (cyl_coord .and. p == 0) then  ! Axisymmetric cylindrical grid
            grid_geometry = 2
        else
            grid_geometry = 3
        end if

        $:GPU_UPDATE(device='[sys_size, buff_size, eqn_idx, adv_n, adap_dt, pi_fac, adap_dt_tol, adap_dt_max_iters]')
        $:GPU_UPDATE(device='[b_size, tensor_size]')

        $:GPU_UPDATE(device='[cfl_target, m, n, p]')

        $:GPU_UPDATE(device='[alt_soundspeed, acoustic_source, num_source]')
        $:GPU_UPDATE(device='[dt, sys_size, buff_size, pref, rhoref, eqn_idx, mpp_lim, bubbles_euler, hypoelasticity, &
                     & alt_soundspeed, avg_state, model_eqns, mixture_err, grid_geometry, cyl_coord, mp_weno, weno_eps, teno_CT, &
                     & hyperelasticity, hyper_model, elasticity, low_Mach]')

        $:GPU_UPDATE(device='[Bx0]')

        $:GPU_UPDATE(device='[chem_params]')

        $:GPU_UPDATE(device='[cont_damage, tau_star, cont_damage_s, alpha_bar]')

        $:GPU_UPDATE(device='[hyper_cleaning, hyper_cleaning_speed, hyper_cleaning_tau]')

        #:if not MFC_CASE_OPTIMIZATION
            $:GPU_UPDATE(device='[wenojs, mapped_weno, wenoz, teno]')
            $:GPU_UPDATE(device='[wenoz_q]')
            $:GPU_UPDATE(device='[mhd, relativity]')
            $:GPU_UPDATE(device='[muscl_order, muscl_lim]')
            $:GPU_UPDATE(device='[igr, igr_order]')
            $:GPU_UPDATE(device='[num_fluids, num_dims, viscous, num_vels, nb, muscl_lim]')
        #:endif

        $:GPU_UPDATE(device='[int_comp, ic_eps, ic_beta]')
        $:GPU_UPDATE(device='[muscl_eps]')
        $:GPU_UPDATE(device='[dir_idx, dir_flg, dir_idx_tau]')

        $:GPU_UPDATE(device='[relax, relax_model, palpha_eps, ptgalpha_eps]')

        if (synthetic_turbulence) then
            $:GPU_UPDATE(device='[synthetic_turbulence, num_turbulent_sources]')
            $:GPU_UPDATE(device='[synth_U_inf, synth_n_waves_per_shell, synth_k_shell, synth_amp_shell]')
            $:GPU_UPDATE(device='[turb_pos, synth_L]')
        end if

        ! Allocating grid variables for the x-, y- and z-directions
        @:ALLOCATE(x_cb(-1 - buff_size:m + buff_size))
        @:ALLOCATE(x_cc(-buff_size:m + buff_size))
        @:ALLOCATE(dx(-buff_size:m + buff_size))
        @:PREFER_GPU(x_cb)
        @:PREFER_GPU(x_cc)
        @:PREFER_GPU(dx)

        if (n == 0) return
        @:ALLOCATE(y_cb(-1 - buff_size:n + buff_size))
        @:ALLOCATE(y_cc(-buff_size:n + buff_size))
        @:ALLOCATE(dy(-buff_size:n + buff_size))
        @:PREFER_GPU(y_cb)
        @:PREFER_GPU(y_cc)
        @:PREFER_GPU(dy)

        if (p == 0) return
        @:ALLOCATE(z_cb(-1 - buff_size:p + buff_size))
        @:ALLOCATE(z_cc(-buff_size:p + buff_size))
        @:ALLOCATE(dz(-buff_size:p + buff_size))
        @:PREFER_GPU(z_cb)
        @:PREFER_GPU(z_cc)
        @:PREFER_GPU(dz)

    end subroutine s_initialize_global_parameters_module

    !> Initializes parallel infrastructure
    impure subroutine s_initialize_parallel_io

        ! proc_coords/start_idx/mpiiofs/mpi_info_int setup moved into the shared routine
        call s_initialize_parallel_io_common

        ! #1290: per-rank physical comm-domain bounds for Lagrangian-bubble exchange
        @:ALLOCATE(pcomm_coords(1:num_dims))

    end subroutine s_initialize_parallel_io

    !> Module deallocation and/or disassociation procedures
    impure subroutine s_finalize_global_parameters_module

        integer :: i

        ! Deallocating the variables bookkeeping the indexes of any viscous fluids and any pairs of fluids whose interfaces
        ! supported effects of surface tension

        if (viscous) then
            @:DEALLOCATE(Re_idx)
        end if

        ! Herschel-Bulkley non-Newtonian viscosity arrays (always allocated)
        @:DEALLOCATE(is_non_newtonian)
        @:DEALLOCATE(hb_tau0, hb_K, hb_nn, hb_m_arr)
        @:DEALLOCATE(hb_mu_min, hb_mu_max)
        @:DEALLOCATE(fluid_inv_re)

        if (bubbles_euler) then
            @:DEALLOCATE(ptil)
            @:DEALLOCATE(qbmm_idx%rs, qbmm_idx%vs, qbmm_idx%ps, qbmm_idx%ms)
            if (qbmm) then
                @:DEALLOCATE(qbmm_idx%moms)
            end if
        end if

        @:DEALLOCATE(pcomm_coords)

        ! Shared: deallocate proc_coords and start_idx
        call s_finalize_global_parameters_common

        if (parallel_io) then
            if (bubbles_lagrange) then
                do i = 1, sys_size + 1
                    MPI_IO_DATA%var(i)%sf => null()
                end do
            else
                do i = 1, sys_size
                    MPI_IO_DATA%var(i)%sf => null()
                end do
            end if

            deallocate (MPI_IO_DATA%var)
            deallocate (MPI_IO_DATA%view)
        end if

        if (ib) MPI_IO_IB_DATA%var%sf => null()

        ! Deallocating grid variables for the x-, y- and z-directions
        @:DEALLOCATE(x_cb, x_cc, dx)

        if (n == 0) return
        @:DEALLOCATE(y_cb, y_cc, dy)

        if (p == 0) return
        @:DEALLOCATE(z_cb, z_cc, dz)

        if (allocated(neighbor_ranks)) then
            @:DEALLOCATE(neighbor_ranks)
        end if

    end subroutine s_finalize_global_parameters_module

end module m_global_parameters
