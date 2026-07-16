!>
!! @file
!! @brief Contains module m_global_parameters

#:include 'case.fpp'

!> @brief Global parameters for the post-process: domain geometry, equation of state, and output database settings
module m_global_parameters

#ifdef MFC_MPI
    use mpi  !< Message passing interface (MPI) module
#endif

    use m_derived_types
    use m_helper_basic
    use m_thermochem, only: species_names
    use m_constants, only: format_silo, precision_single
    ! Shared state: generated_decls, num_dims, num_vels, sys_size, eqn_idx, b_size, tensor_size, chemistry, elasticity, shear_*
    use m_global_parameters_common

    implicit none

    !> @name Logistics
    !> @{
    integer :: num_procs  !< Number of processors
    !> @}

    ! Computational Domain Parameters

    integer :: proc_rank  !< Rank of the local processor
    !> @name Number of cells in the x-, y- and z-coordinate directions
    !> @{
    integer :: m_root
    !> @}

    !> @name Max and min number of cells in a direction of each combination of x-,y-, and z-
    type(cell_num_bounds) :: cells_bounds
    integer(kind=8)       :: nGlobal  !< Total number of cells in global domain

    !> @name Cylindrical coordinates (either axisymmetric or full 3D)
    !> @{
    integer :: grid_geometry
    !> @}

    !> @name Global number of cells in each direction
    !> @{
    integer :: m_glb, n_glb, p_glb
    !> @}

    ! num_dims, num_vels: in m_global_parameters_common
    !> @name Cell-boundary locations in the x-, y- and z-coordinate directions
    !> @{
    real(wp), allocatable, dimension(:) :: x_cb, x_root_cb, y_cb, z_cb
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

    integer :: buff_size  !< Number of ghost cells for boundary condition storage
    !> @name IO options for adaptive time-stepping
    !> @{
    logical :: cfl_dt
    integer :: n_save
    !> @}

    ! NOTE: m_root, x_root_cb, x_root_cc = defragmented grid (1D only; equals m, x_cb, x_cc in serial)

    !> @name Simulation Algorithm Parameters
    !> @{
    ! sys_size, elasticity, b_size, tensor_size, chemistry, eqn_idx: in m_global_parameters_common
    !> @}

    !> @name Annotations of the structure, i.e. the organization, of the state vectors
    !> @{
    type(qbmm_idx_info) :: qbmm_idx  !< QBMM moment index mappings.
    integer             :: beta_idx  !< Index of lagrange bubbles beta
    !> @}

    ! Cell Indices for the (local) interior points (O-m, O-n, 0-p). Stands for "InDices With BUFFer".
    type(int_bounds_info) :: idwint(1:3)

    ! Cell indices (InDices With BUFFer): includes buffer in simulation only
    type(int_bounds_info) :: idwbuff(1:3)
    logical               :: bc_io
    !> @name Boundary conditions in the x-, y- and z-coordinate directions
    !> @{
    type(int_bounds_info) :: bc_x, bc_y, bc_z
    type(bc_xyz_info)     :: bc
    !> @}

    ! shear_num, shear_indices, shear_BC_flip_num, shear_BC_flip_indices: in m_global_parameters_common
    ! proc_coords, start_idx, mpiiofs, mpi_info_int: in m_global_parameters_common
    type(int_bounds_info), dimension(3)                    :: nidx  !< Neighbor index offsets per direction (#1290 decomposition)
    integer, allocatable, dimension(:,:,:)                 :: neighbor_ranks  !< MPI ranks of neighbors (#1290 decomposition)
    type(ib_airfoil_parameters), allocatable, dimension(:) :: ib_airfoil  !< Per-airfoil NACA parameters (unused in post_process)
    !> Per-airfoil computed surface grids (unused in post_process)
    type(ib_airfoil_grid), allocatable, dimension(:) :: ib_airfoil_grids

#ifdef MFC_MPI
    type(mpi_io_var), public                      :: MPI_IO_DATA
    type(mpi_io_ib_var), public                   :: MPI_IO_IB_DATA
    type(mpi_io_levelset_var), public             :: MPI_IO_levelset_DATA
    type(mpi_io_levelset_norm_var), public        :: MPI_IO_levelsetnorm_DATA
    real(wp), allocatable, dimension(:,:), public :: MPI_IO_DATA_lg_bubbles
#endif

    ! fluid_pp, bub_pp: auto-generated in generated_decls.fpp
    real(wp), allocatable, dimension(:) :: adv  !< Advection variables
    ! Formatted Database File(s) Structure Parameters

    type(bounds_info)     :: x_output, y_output, z_output              !< Portion of domain to output for post-processing
    type(int_bounds_info) :: x_output_idx, y_output_idx, z_output_idx  !< Indices of domain to output for post-processing
    !> @name Size of the ghost zone layer in the x-, y- and z-coordinate directions. The definition of the ghost zone layers is only
    !! necessary when using the Silo database file format in multidimensions. These zones provide VisIt with the subdomain
    !! connectivity information that it requires in order to produce smooth plots.
    !> @{
    type(int_bounds_info) :: offset_x, offset_y, offset_z
    !> @}

    ! alpha_rho_wrt, mom_wrt, vel_wrt, flux_wrt, alpha_rho_e_wrt, alpha_wrt,
    ! omega_wrt, chem_wrt_Y, schlieren_alpha: auto-generated in generated_decls.fpp
    integer                    :: fd_number  !< Finite-difference half-stencil size: MAX(1, fd_order/2)
    type(chemistry_parameters) :: chem_params
    !> @name Bubble modeling variables and parameters
    !> @{
    real(wp) :: Eu
    real(wp), dimension(:), allocatable :: weight, R0
    real(wp) :: phi_vg, phi_gv, Pe_c, Tw, k_vl, k_gl
    real(wp) :: gam_m
    real(wp), dimension(:), allocatable :: pb0, mass_g0, mass_v0, Pe_T, k_v, k_g
    real(wp), dimension(:), allocatable :: Re_trans_T, Re_trans_c, Im_trans_T, Im_trans_c, omegaN
    real(wp) :: p0ref, rho0ref, T0ref, ss, pv, vd, mu_l, mu_v, mu_g, gam_v, gam_g, M_v, M_g, cp_v, cp_g, R_v, R_g
    real(wp) :: G
    integer :: nmom
    !> @}

    real(wp) :: wall_time, wall_time_avg  !< Wall time measurements

contains

    !> Assigns default values to user inputs prior to reading them in. This allows for an easier consistency check of these
    !! parameters once they are read from the input file.
    impure subroutine s_assign_default_values_to_user_inputs

        integer :: i  !< Generic loop iterator

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

        ! Computational domain parameters (post-specific)
        m_root = dflt_int

        t_step_stop = dflt_int
        t_step_save = dflt_int

        cfl_dt = .false.
        cfl_target = dflt_real
        t_save = dflt_real
        t_stop = dflt_real

        ! AMR: post_process overlays the refined fine blocks when this is on (default off)
        amr = .false.

        ! Simulation algorithm parameters (post-specific)
        mixture_err = .false.
        alt_soundspeed = .false.

        bc_io = .false.
        num_bc_patches = dflt_int

        chem_params%gamma_method = 1
        chem_params%transport_model = 1

        ! Fluids physical parameters (post-specific; G = dflt_real differs from pre/sim)
        do i = 1, num_fluids_max
            fluid_pp(i)%gamma = dflt_real
            fluid_pp(i)%pi_inf = dflt_real
            fluid_pp(i)%cv = 0._wp
            fluid_pp(i)%qv = 0._wp
            fluid_pp(i)%qvp = 0._wp
            fluid_pp(i)%G = dflt_real
            fluid_pp(i)%non_newtonian = .false.
            fluid_pp(i)%K = dflt_real
            fluid_pp(i)%nn = dflt_real
            fluid_pp(i)%tau0 = 0._wp
            fluid_pp(i)%hb_m = dflt_real
            fluid_pp(i)%mu_min = dflt_real
            fluid_pp(i)%mu_max = dflt_real
            fluid_pp(i)%mu_bulk = dflt_real
        end do

        ! Subgrid bubble parameters (bub_pp struct + scalar companions; bub_pp%R0ref is set in common
        ! via R0ref; the scalar companions are per-target manual declarations)
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

        ! Formatted database file(s) structure parameters (post-specific)
        format = dflt_int

        precision = dflt_int

        alpha_rho_wrt = .false.
        alpha_rho_e_wrt = .false.
        rho_wrt = .false.
        mom_wrt = .false.
        vel_wrt = .false.
        chem_wrt_Y = .false.
        chem_wrt_T = .false.
        flux_lim = dflt_int
        flux_wrt = .false.
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
        liutex_wrt = .false.
        schlieren_wrt = .false.
        sim_data = .false.
        cf_wrt = .false.
        ib_state_wrt = .false.
        lag_txt_wrt = .false.
        lag_header = .true.
        lag_db_wrt = .false.
        lag_id_wrt = .true.
        lag_pos_wrt = .true.
        lag_pos_prev_wrt = .false.
        lag_vel_wrt = .true.
        lag_rad_wrt = .true.
        lag_rvel_wrt = .false.
        lag_r0_wrt = .false.
        lag_rmax_wrt = .false.
        lag_rmin_wrt = .false.
        lag_dphidt_wrt = .false.
        lag_pres_wrt = .false.
        lag_mv_wrt = .false.
        lag_mg_wrt = .false.
        lag_betaT_wrt = .false.
        lag_betaC_wrt = .false.

        schlieren_alpha = dflt_real

        fd_order = dflt_int
        avg_state = dflt_int

        ! Bubble modeling (post-specific)
        nb = dflt_int
        sigR = dflt_real

        ! Output partial domain (post-specific)
        output_partial_domain = .false.
        x_output%beg = dflt_real
        x_output%end = dflt_real
        y_output%beg = dflt_real
        y_output%end = dflt_real
        z_output%beg = dflt_real
        z_output%end = dflt_real

    end subroutine s_assign_default_values_to_user_inputs

    !> Computation of parameters, allocation procedures, and/or any other tasks needed to properly setup the module
    impure subroutine s_initialize_global_parameters_module

        integer :: i, j, fac

        ! Setting m_root equal to m in the case of a 1D serial simulation

        if (n == 0) m_root = m_glb

        ! Gamma/Pi_inf: force num_fluids=1 (post_process-specific side effect of the gamma-law model)
        if (model_eqns == model_eqns_gamma_law) num_fluids = 1

        ! post_process sets nmom to 6 for qbmm before the shared eqn_idx setup
        ! (guard matches the original site: inside the 5-equation branch)
        if (model_eqns == model_eqns_5eq .and. qbmm) nmom = 6

        ! Populate eqn_idx, sys_size, b_size, tensor_size, elasticity, shear_* (shared logic)
        call s_initialize_eqn_idx(nmom, nb)

        ! post-only: 6eq alf is a dummy (no void fraction in 6eq)
        if (model_eqns == model_eqns_6eq) eqn_idx%alf = 1

        ! post-only: set default indices for disabled fields (used by post-processing consumers)
        if (model_eqns == model_eqns_5eq .or. model_eqns == model_eqns_6eq) then
            if (.not. cont_damage) eqn_idx%damage = dflt_int
            if (.not. hyper_cleaning) eqn_idx%psi = dflt_int
        end if

        ! post-only: species defaults when chemistry is off
        if (.not. chemistry) then
            eqn_idx%species%beg = 1
            eqn_idx%species%end = 1
        end if

        ! Per-target (post_process): beta_idx for bubbles_lagrange (5eq only, after main eqn_idx setup)
        if (model_eqns == model_eqns_5eq .and. bubbles_lagrange) then
            beta_idx = sys_size + 1
            sys_size = beta_idx
        end if

        ! Per-target (post_process): qbmm_idx allocations and fills
        if (model_eqns == model_eqns_5eq .and. bubbles_euler) then
            allocate (qbmm_idx%rs(nb), qbmm_idx%vs(nb))
            allocate (qbmm_idx%ps(nb), qbmm_idx%ms(nb))

            if (qbmm) then
                allocate (qbmm_idx%moms(nb, nmom))
                do i = 1, nb
                    do j = 1, nmom
                        qbmm_idx%moms(i, j) = eqn_idx%bub%beg + (j - 1) + (i - 1)*nmom
                    end do
                    qbmm_idx%rs(i) = qbmm_idx%moms(i, 2)
                    qbmm_idx%vs(i) = qbmm_idx%moms(i, 3)
                end do
            else
                do i = 1, nb
                    if (polytropic .neqv. .true.) then
                        fac = 4
                    else
                        fac = 2
                    end if

                    qbmm_idx%rs(i) = eqn_idx%bub%beg + (i - 1)*fac
                    qbmm_idx%vs(i) = qbmm_idx%rs(i) + 1

                    if (polytropic .neqv. .true.) then
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
                if (polytropic .neqv. .true.) then
                    fac = 4
                else
                    fac = 2
                end if

                qbmm_idx%rs(i) = eqn_idx%bub%beg + (i - 1)*fac
                qbmm_idx%vs(i) = qbmm_idx%rs(i) + 1

                if (polytropic .neqv. .true.) then
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

        if (output_partial_domain) then
            x_output_idx%beg = 0
            x_output_idx%end = 0
            y_output_idx%beg = 0
            y_output_idx%end = 0
            z_output_idx%beg = 0
            z_output_idx%end = 0
        end if

#ifdef MFC_MPI
        if (qbmm .and. .not. polytropic) then
            allocate (MPI_IO_DATA%view(1:sys_size + 2*nb*nnode))
            allocate (MPI_IO_DATA%var(1:sys_size + 2*nb*nnode))
        else
            allocate (MPI_IO_DATA%view(1:sys_size))
            allocate (MPI_IO_DATA%var(1:sys_size))
        end if

        do i = 1, sys_size
            if (down_sample) then
                allocate (MPI_IO_DATA%var(i)%sf(-1:m + 1,-1:n + 1,-1:p + 1))
            else
                allocate (MPI_IO_DATA%var(i)%sf(0:m,0:n,0:p))
            end if
            MPI_IO_DATA%var(i)%sf => null()
        end do
        if (qbmm .and. .not. polytropic) then
            do i = sys_size + 1, sys_size + 2*nb*nnode
                allocate (MPI_IO_DATA%var(i)%sf(0:m,0:n,0:p))
                MPI_IO_DATA%var(i)%sf => null()
            end do
        end if

        if (ib) allocate (MPI_IO_IB_DATA%var%sf(0:m,0:n,0:p))
#endif

        ! Size of the ghost zone layer is non-zero only when post-processing the raw simulation data of a parallel multidimensional
        ! computation in the Silo-HDF5 format. If this is the case, one must also verify whether the raw simulation data is 2D or
        ! 3D. In the 2D case, size of the z-coordinate direction ghost zone layer must be zeroed out.
        if (num_procs == 1 .or. format /= format_silo) then
            offset_x%beg = 0
            offset_x%end = 0
            offset_y%beg = 0
            offset_y%end = 0
            offset_z%beg = 0
            offset_z%end = 0
        else if (n == 0) then
            offset_y%beg = 0
            offset_y%end = 0
            offset_z%beg = 0
            offset_z%end = 0
        else if (p == 0) then
            offset_z%beg = 0
            offset_z%end = 0
        end if

        ! Determining the finite-difference number and the buffer size. Note that the size of the buffer is unrelated to the order
        ! of the WENO scheme. Rather, it is directly dependent on maximum size of ghost zone layers and possibly the order of the
        ! finite difference scheme used for the computation of vorticity and/or numerical Schlieren function.
        buff_size = max(offset_x%beg, offset_x%end, offset_y%beg, offset_y%end, offset_z%beg, offset_z%end)

        if (any(omega_wrt) .or. schlieren_wrt .or. qm_wrt .or. liutex_wrt) then
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
        allocate (x_cc_s(-buff_size:m + buff_size))

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

            ! Allocating the grid variables, only used for the 1D simulations, and containing the defragmented computational domain
            ! grid data
        else
            allocate (x_root_cb(-1:m_root))
            allocate (x_root_cc(0:m_root))

            if (precision == precision_single) then
                allocate (x_root_cc_s(0:m_root))
            end if
        end if

        allocate (adv(num_fluids))

        if (cyl_coord .neqv. .true.) then  ! Cartesian grid
            grid_geometry = 1
        else if (cyl_coord .and. p == 0) then  ! Axisymmetric cylindrical grid
            grid_geometry = 2
        else  ! Fully 3D cylindrical grid
            grid_geometry = 3
        end if

    end subroutine s_initialize_global_parameters_module

    !> Subroutine to initialize parallel infrastructure
    impure subroutine s_initialize_parallel_io

        call s_initialize_parallel_io_common

    end subroutine s_initialize_parallel_io

    !> Deallocation procedures for the module
    impure subroutine s_finalize_global_parameters_module

        integer :: i

        if (bubbles_euler) then
            deallocate (qbmm_idx%rs, qbmm_idx%vs, qbmm_idx%ps, qbmm_idx%ms)
            if (qbmm) deallocate (qbmm_idx%moms)
        end if

        ! Deallocating the grid variables for the x-coordinate direction
        deallocate (x_cc, x_cb, dx)

        ! Deallocating grid variables for the y- and z-coordinate directions
        if (n > 0) then
            deallocate (y_cc, y_cb, dy)
            if (p > 0) then
                deallocate (z_cc, z_cb, dz)
            end if
        else
            ! Deallocating the grid variables, only used for the 1D simulations, and containing the defragmented computational
            ! domain grid data
            deallocate (x_root_cb, x_root_cc)
        end if

        ! Shared: deallocate proc_coords and start_idx
        call s_finalize_global_parameters_common

        deallocate (adv)

#ifdef MFC_MPI
        if (parallel_io) then
            do i = 1, sys_size
                MPI_IO_DATA%var(i)%sf => null()
            end do

            deallocate (MPI_IO_DATA%var)
            deallocate (MPI_IO_DATA%view)
        end if

        if (ib) MPI_IO_IB_DATA%var%sf => null()
#endif

        if (allocated(neighbor_ranks)) deallocate (neighbor_ranks)

    end subroutine s_finalize_global_parameters_module

end module m_global_parameters
