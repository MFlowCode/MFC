!>
!! @file
!! @brief Contains module m_global_parameters

#:include 'case.fpp'

!> @brief Defines global parameters for the computational domain, simulation algorithm, and initial conditions
module m_global_parameters

#ifdef MFC_MPI
    use mpi  ! Message passing interface (MPI) module
#endif

    use m_derived_types  ! Definitions of the derived types
    use m_helper_basic  ! Functions to compare floating point numbers
    ! Shared state: generated_decls, sys_size, eqn_idx, b_size, tensor_size, chemistry, elasticity, shear_*
    use m_global_parameters_common

    implicit none

    ! Logistics
    integer :: num_procs     !< Number of processors
    logical :: non_axis_sym  !< Use existing IC data
    logical :: cfl_dt

    ! Computational Domain Parameters

    integer :: proc_rank  !< Rank of the local processor Number of cells in the x-, y- and z-coordinate directions

    !> @name Max and min number of cells in a direction of each combination of x-,y-, and z-
    type(cell_num_bounds) :: cells_bounds
    integer(kind=8)       :: nGlobal              !< Global number of cells in the domain
    integer               :: m_glb, n_glb, p_glb  !< Global number of cells in each direction
    integer               :: grid_geometry        !< Cylindrical coordinates (either axisymmetric or full 3D)
    !> Locations of cell-centers (cc) in x-, y- and z-directions, respectively
    real(wp), allocatable, dimension(:) :: x_cc, y_cc, z_cc
    !> Locations of cell-boundaries (cb) in x-, y- and z-directions, respectively
    real(wp), allocatable, dimension(:) :: x_cb, y_cb, z_cb
    real(wp) :: dx, dy, dz                             !< Minimum cell-widths in the x-, y- and z-coordinate directions
    type(bounds_info) :: x_domain, y_domain, z_domain  !< Locations of the domain bounds in the x-, y- and z-coordinate directions
    !> Global (pre-decomposition) domain bounds, needed by s_generate_serial_grid to stretch the grid using the full domain length
    !! rather than a local processor's sub-domain length
    type(bounds_info) :: x_domain_glb, y_domain_glb, z_domain_glb

    ! Simulation Algorithm Parameters
    ! sys_size, eqn_idx, b_size, tensor_size, chemistry, elasticity, shear_*: in m_global_parameters_common
    ! weno_polyn, muscl_polyn, num_dims, num_vels: in m_global_parameters_common
    ! Annotations of the structure, i.e. the organization, of the state vectors
    type(qbmm_idx_info) :: qbmm_idx  !< QBMM moment index mappings.
    ! Cell Indices for the (local) interior points (O-m, O-n, 0-p). Stands for "InDices With BUFFer".
    type(int_bounds_info) :: idwint(1:3)

    ! Cell indices (InDices With BUFFer): includes buffer except in pre_process
    type(int_bounds_info) :: idwbuff(1:3)
    type(int_bounds_info) :: bc_x, bc_y, bc_z  !< Boundary conditions in the x-, y- and z-coordinate directions
    ! simplex_params: auto-generated in generated_decls.fpp

    ! fluid_rho (perturbs surrounding-air density to break grid symmetry): auto-generated in generated_decls.fpp
    ! proc_coords, start_idx, mpiiofs, mpi_info_int: in m_global_parameters_common
#ifdef MFC_MPI
    type(mpi_io_var), public :: MPI_IO_DATA
#endif

    ! Initial Condition Parameters patch_icpp, patch_bc: auto-generated in generated_decls.fpp
    logical :: bc_io  !< whether or not to save BC data

    ! Fluids Physical Parameters fluid_pp, bub_pp: auto-generated in generated_decls.fpp
    type(chemistry_parameters) :: chem_params
    !> @name Bubble modeling
    !> @{
    real(wp)                            :: Eu
    real(wp), dimension(:), allocatable :: weight, R0
    integer                             :: nmom  !< Number of carried moments
    !> @}

    !> @name Immersed Boundaries
    !> @{
    ! patch_ib, ib_airfoil, stl_models: auto-generated in generated_decls.fpp
    !> Per-airfoil computed surface grids (unused in pre_process)
    type(ib_airfoil_grid), allocatable, dimension(:) :: ib_airfoil_grids
    !> @}

    !> @name Non-polytropic bubble gas compression
    !> @{
    real(wp) :: phi_vg, phi_gv, Pe_c, Tw, k_vl, k_gl
    real(wp) :: gam_m
    real(wp), dimension(:), allocatable :: pb0, mass_g0, mass_v0, Pe_T, k_v, k_g
    real(wp), dimension(:), allocatable :: Re_trans_T, Re_trans_c, Im_trans_T, Im_trans_c, omegaN
    real(wp) :: p0ref, rho0ref, T0ref, ss, pv, vd, mu_l, mu_v, mu_g, gam_v, gam_g, M_v, M_g, cp_v, cp_g, R_v, R_g
    !> @}

    integer, allocatable, dimension(:,:,:) :: logic_grid
    type(pres_field)                       :: pb
    type(pres_field)                       :: mv
    integer                                :: buff_size  !< Number of ghost cells for boundary condition storage

contains

    !> Assigns default values to user inputs prior to reading them in. This allows for an easier consistency check of these
    !! parameters once they are read from the input file.
    impure subroutine s_assign_default_values_to_user_inputs

        integer :: i  !< Generic loop operator

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

        ! Logistics (pre-specific)
        file_extension = '000000'
        files_dir = './'
        old_grid = .false.
        old_ic = .false.
        t_step_old = dflt_int
        cfl_dt = .false.

        ! Computational domain parameters (pre-specific)
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

        ! Simulation algorithm parameters (pre-specific)
        palpha_eps = dflt_real
        ptgalpha_eps = dflt_real
        igr_order = dflt_int
        pre_stress = .false.

        precision = 2
        viscous = .false.
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

        simplex_perturb = .false.
        simplex_params%perturb_vel(:) = .false.
        simplex_params%perturb_vel_freq(:) = dflt_real
        simplex_params%perturb_vel_scale(:) = dflt_real
        simplex_params%perturb_vel_offset(:,:) = dflt_real
        simplex_params%perturb_dens(:) = .false.
        simplex_params%perturb_dens_freq(:) = dflt_real
        simplex_params%perturb_dens_scale(:) = dflt_real
        simplex_params%perturb_dens_offset(:,:) = dflt_real

        ! Initial condition parameters
        num_patches = dflt_int

        do i = 1, num_patches_max
            patch_icpp(i)%geometry = dflt_int
            patch_icpp(i)%model_id = 0
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
            patch_icpp(i)%fourier_cos(:) = 0._wp
            patch_icpp(i)%fourier_sin(:) = 0._wp
            patch_icpp(i)%modal_clip_r_to_min = .false.
            patch_icpp(i)%modal_r_min = 1.e-12_wp
            patch_icpp(i)%modal_use_exp_form = .false.
            patch_icpp(i)%sph_har_coeff(:,:) = 0._wp

            ! should get all of r0's and v0's
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

        ! Bubble modeling (pre-specific)
        polytropic = .true.
        thermal = dflt_int
        nb = dflt_int

        Eu = dflt_real
        Ca = dflt_real
        Re_inv = dflt_real
        Web = dflt_real

        nmom = 1
        sigR = dflt_real
        sigV = dflt_real
        rhoRV = 0._wp
        dist_type = dflt_int

        R_g = dflt_real
        R_v = dflt_real
        phi_vg = dflt_real
        phi_gv = dflt_real
        Pe_c = dflt_real
        Tw = dflt_real

        pi_fac = 1._wp

        do i = 1, num_ib_patches_max_namelist
            patch_ib(i)%geometry = dflt_int
            patch_ib(i)%x_centroid = dflt_real
            patch_ib(i)%y_centroid = dflt_real
            patch_ib(i)%z_centroid = dflt_real
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

        do i = 1, num_ib_airfoils_max
            ib_airfoil(i)%c = dflt_real
            ib_airfoil(i)%p = dflt_real
            ib_airfoil(i)%t = dflt_real
            ib_airfoil(i)%m = dflt_real
        end do

        num_stl_models = 0

        do i = 1, num_stl_models_max
            stl_models(i)%model_filepath(:) = dflt_char
            stl_models(i)%model_translate(:) = 0._wp
            stl_models(i)%model_scale(:) = 1._wp
            stl_models(i)%model_threshold = ray_tracing_threshold
        end do

        chem_params%gamma_method = 1
        chem_params%transport_model = 1

        chem_params%reaction_substeps = 0
        chem_params%adap_substeps = .false.
        chem_params%reaction_substeps_max = 0

        ! Fluids physical parameters
        do i = 1, num_fluids_max
            fluid_pp(i)%gamma = dflt_real
            fluid_pp(i)%pi_inf = dflt_real
            fluid_pp(i)%cv = 0._wp
            fluid_pp(i)%qv = 0._wp
            fluid_pp(i)%qvp = 0._wp
            fluid_pp(i)%G = 0._wp
            fluid_pp(i)%non_newtonian = .false.
            fluid_pp(i)%K = dflt_real
            fluid_pp(i)%nn = dflt_real
            fluid_pp(i)%tau0 = 0._wp
            fluid_pp(i)%hb_m = dflt_real
            fluid_pp(i)%mu_min = dflt_real
            fluid_pp(i)%mu_max = dflt_real
            fluid_pp(i)%mu_bulk = dflt_real
        end do

        ! Subgrid bubble parameters
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

    end subroutine s_assign_default_values_to_user_inputs

    !> Computation of parameters, allocation procedures, and/or any other tasks needed to properly setup the module
    impure subroutine s_initialize_global_parameters_module

        integer :: i, j, fac

        if (recon_type == recon_type_weno) then
            weno_polyn = (weno_order - 1)/2
        else if (recon_type == recon_type_muscl) then
            muscl_polyn = muscl_order
        end if

        ! Gamma/Pi_inf: force num_fluids=1 (pre_process-specific side effect of the gamma-law model)
        if (model_eqns == model_eqns_gamma_law) num_fluids = 1

        ! Pre-process sets nmom to 6 for qbmm before the shared eqn_idx setup
        ! (guards match the original site: 5-equation bubbles with 4-node qbmm)
        if (model_eqns == model_eqns_5eq .and. bubbles_euler .and. qbmm .and. nnode == 4) nmom = 6

        ! Populate eqn_idx, sys_size, b_size, tensor_size, elasticity, shear_* (shared logic)
        call s_initialize_eqn_idx(nmom, nb)

        ! Per-target (pre_process): qbmm_idx allocations and fills
        if (model_eqns == model_eqns_5eq .and. bubbles_euler) then
            allocate (qbmm_idx%rs(nb), qbmm_idx%vs(nb))
            allocate (qbmm_idx%ps(nb), qbmm_idx%ms(nb))

            if (qbmm) then
                allocate (qbmm_idx%moms(nb, nmom))
                allocate (qbmm_idx%fullmom(nb,0:nmom,0:nmom))

                do i = 1, nb
                    do j = 1, nmom
                        qbmm_idx%moms(i, j) = eqn_idx%bub%beg + (j - 1) + (i - 1)*nmom
                    end do
                    qbmm_idx%fullmom(i, 0, 0) = qbmm_idx%moms(i, 1)
                    qbmm_idx%fullmom(i, 1, 0) = qbmm_idx%moms(i, 2)
                    qbmm_idx%fullmom(i, 0, 1) = qbmm_idx%moms(i, 3)
                    qbmm_idx%fullmom(i, 2, 0) = qbmm_idx%moms(i, 4)
                    qbmm_idx%fullmom(i, 1, 1) = qbmm_idx%moms(i, 5)
                    qbmm_idx%fullmom(i, 0, 2) = qbmm_idx%moms(i, 6)
                    qbmm_idx%rs(i) = qbmm_idx%fullmom(i, 1, 0)
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
            allocate (qbmm_idx%rs(nb), qbmm_idx%vs(nb))
            allocate (qbmm_idx%ps(nb), qbmm_idx%ms(nb))
            allocate (weight(nb), R0(nb))

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

        call s_configure_coordinate_bounds(recon_type, weno_polyn, muscl_polyn, igr_order, buff_size, idwint, idwbuff, viscous, &
                                           & bubbles_lagrange, m, n, p, num_dims, igr, ib)

#ifdef MFC_MPI
        if (qbmm .and. .not. polytropic) then
            allocate (MPI_IO_DATA%view(1:sys_size + 2*nb*nnode))
            allocate (MPI_IO_DATA%var(1:sys_size + 2*nb*nnode))
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
        if (qbmm .and. .not. polytropic) then
            do i = sys_size + 1, sys_size + 2*nb*nnode
                allocate (MPI_IO_DATA%var(i)%sf(0:m,0:n,0:p))
                MPI_IO_DATA%var(i)%sf => null()
            end do
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

        if (cyl_coord .neqv. .true.) then  ! Cartesian grid
            grid_geometry = 1
        else if (cyl_coord .and. p == 0) then  ! Axisymmetric cylindrical grid
            grid_geometry = 2
        else  ! Fully 3D cylindrical grid
            grid_geometry = 3
        end if

        if (.not. igr) then
            allocate (logic_grid(0:m,0:n,0:p))
        end if

    end subroutine s_initialize_global_parameters_module

    !> Configure MPI parallel I/O settings and allocate processor coordinate arrays.
    impure subroutine s_initialize_parallel_io

        call s_initialize_parallel_io_common

    end subroutine s_initialize_parallel_io

    !> Deallocate all global grid, index, and equation-of-state parameter arrays.
    impure subroutine s_finalize_global_parameters_module

        integer :: i

        if (bubbles_euler) then
            deallocate (qbmm_idx%rs, qbmm_idx%vs, qbmm_idx%ps, qbmm_idx%ms)
            if (qbmm) deallocate (qbmm_idx%moms, qbmm_idx%fullmom)
        end if

        ! Deallocating grid variables for the x-direction
        deallocate (x_cc, x_cb)
        ! Deallocating grid variables for the y- and z-directions
        if (n > 0) then
            deallocate (y_cc, y_cb)
            if (p > 0) then
                deallocate (z_cc, z_cb)
            end if
        end if

        ! Shared: deallocate proc_coords and start_idx
        call s_finalize_global_parameters_common

#ifdef MFC_MPI
        if (parallel_io) then
            do i = 1, sys_size
                MPI_IO_DATA%var(i)%sf => null()
            end do

            deallocate (MPI_IO_DATA%var)
            deallocate (MPI_IO_DATA%view)
        end if
#endif

    end subroutine s_finalize_global_parameters_module

end module m_global_parameters
