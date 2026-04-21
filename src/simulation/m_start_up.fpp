!>
!! @file
!! @brief Contains module m_start_up

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief Reads input files, loads initial conditions and grid data, and orchestrates solver initialization and finalization
module m_start_up

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_mpi_common
    use m_variables_conversion
    use m_weno
    use m_muscl
    use m_riemann_solvers
    use m_cbc
    use m_boundary_common
    use m_acoustic_src
    use m_rhs
    use m_chemistry
    use m_data_output
    use m_time_steppers
    use m_qbmm
    use m_derived_variables
    use m_hypoelastic
    use m_hyperelastic
    use m_phase_change
    use m_viscous
    use m_bubbles_EE
    use m_bubbles_EL
    use ieee_arithmetic
    use m_helper_basic
    use m_helper

    $:USE_GPU_MODULE()

    use m_nvtx
    use m_ibm
    use m_compile_specific
    use m_checker_common
    use m_checker
    use m_surface_tension
    use m_body_forces
    use m_sim_helpers
    use m_igr

    implicit none

    private; public :: s_read_input_file, s_check_input_file, s_read_data_files, s_read_serial_data_files, &
        & s_read_parallel_data_files, s_initialize_internal_energy_equations, s_initialize_modules, s_initialize_gpu_vars, &
        & s_initialize_mpi_domain, s_finalize_modules, s_perform_time_step, s_save_data, s_save_performance_metrics

    type(scalar_field), allocatable, dimension(:) :: q_cons_temp
    real(wp)                                      :: dt_init

contains

    !> Read data files. Dispatch subroutine that replaces procedure pointer.
    impure subroutine s_read_data_files(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf

        if (.not. parallel_io) then
            call s_read_serial_data_files(q_cons_vf)
        else
            call s_read_parallel_data_files(q_cons_vf)
        end if

    end subroutine s_read_data_files

    !> Verify the input file exists and read it
    impure subroutine s_read_input_file

        character(LEN=name_len), parameter :: file_path = './simulation.inp'
        logical                            :: file_exist  !< Logical used to check the existence of the input file
        integer                            :: iostatus
        ! Integer to check iostat of file read

        character(len=1000) :: line

        namelist /user_inputs/ case_dir, run_time_info, m, n, p, dt, &
            t_step_start, t_step_stop, t_step_save, t_step_print, &
            model_eqns, mpp_lim, time_stepper, weno_eps, &
            rdma_mpi, teno_CT, mp_weno, weno_avg, &
            riemann_solver, low_Mach, wave_speeds, avg_state, &
            bc_x, bc_y, bc_z, &
            x_a, y_a, z_a, x_b, y_b, z_b, &
            x_domain, y_domain, z_domain, &
            hypoelasticity, &
            ib, num_ibs, patch_ib, &
            collision_model, coefficient_of_restitution, collision_time, &
            ib_coefficient_of_friction, ib_state_wrt, &
            fluid_pp, bub_pp, probe_wrt, prim_vars_wrt, &
            fd_order, probe, num_probes, t_step_old, &
            alt_soundspeed, mixture_err, weno_Re_flux, &
            null_weights, precision, parallel_io, cyl_coord, &
            rhoref, pref, bubbles_euler, bubble_model, &
            R0ref, chem_params, &
        #:if not MFC_CASE_OPTIMIZATION
            nb, mapped_weno, wenoz, teno, wenoz_q, weno_order, &
            num_fluids, mhd, relativity, igr_order, viscous, &
            igr_iter_solver, igr, igr_pres_lim, &
            recon_type, muscl_order, muscl_lim, &
        #:endif
        Ca, Web, Re_inv, acoustic_source, acoustic, num_source, polytropic, thermal, integral, integral_wrt, num_integrals, &
            & polydisperse, poly_sigma, qbmm, relax, relax_model, palpha_eps, ptgalpha_eps, file_per_process, sigma, pi_fac, &
            & adv_n, adap_dt, adap_dt_tol, adap_dt_max_iters, bf_x, bf_y, bf_z, k_x, k_y, k_z, w_x, w_y, w_z, p_x, p_y, p_z, g_x, &
            & g_y, g_z, n_start, t_save, t_stop, cfl_adap_dt, cfl_const_dt, cfl_target, surface_tension, bubbles_lagrange, &
            & lag_params, hyperelasticity, R0ref, num_bc_patches, Bx0, cont_damage, tau_star, cont_damage_s, alpha_bar, &
            & hyper_cleaning, hyper_cleaning_speed, hyper_cleaning_tau, alf_factor, num_igr_iters, num_igr_warm_start_iters, &
            & int_comp, ic_eps, ic_beta, nv_uvm_out_of_core, nv_uvm_igr_temps_on_gpu, nv_uvm_pref_gpu, down_sample, fft_wrt

        inquire (FILE=trim(file_path), EXIST=file_exist)

        if (file_exist) then
            open (1, FILE=trim(file_path), form='formatted', ACTION='read', STATUS='old')
            read (1, NML=user_inputs, iostat=iostatus)

            if (iostatus /= 0) then
                backspace (1)
                read (1, fmt='(A)') line
                print *, 'Invalid line in namelist: ' // trim(line)
                call s_mpi_abort('Invalid line in simulation.inp. It is ' // 'likely due to a datatype mismatch. Exiting.')
            end if

            close (1)

            if ((bf_x) .or. (bf_y) .or. (bf_z)) then
                bodyForces = .true.
            end if

            m_glb = m
            n_glb = n
            p_glb = p

            call s_update_cell_bounds(cells_bounds, m, n, p)

            if (cfl_adap_dt .or. cfl_const_dt) cfl_dt = .true.

            if (any((/bc_x%beg, bc_x%end, bc_y%beg, bc_y%end, bc_z%beg, bc_z%end/) == -17) .or. num_bc_patches > 0) then
                bc_io = .true.
            end if
        else
            call s_mpi_abort(trim(file_path) // ' is missing. Exiting.')
        end if

    end subroutine s_read_input_file

    !> Validate that all user-provided inputs form a consistent simulation configuration
    impure subroutine s_check_input_file

        character(LEN=path_len) :: file_path
        logical                 :: file_exist

        file_path = trim(case_dir) // '/.'

        call my_inquire(file_path, file_exist)

        if (file_exist .neqv. .true.) then
            call s_mpi_abort(trim(file_path) // ' is missing. Exiting.')
        end if

        call s_check_inputs_common()
        call s_check_inputs()

    end subroutine s_check_input_file

    !> Read serial initial condition and grid data files and compute cell-width distributions
    impure subroutine s_read_serial_data_files(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        character(LEN=path_len + 2*name_len) :: t_step_dir  !< Relative path to the starting time-step directory
        character(LEN=path_len + 3*name_len) :: file_path   !< Relative path to the grid and conservative variables data files
        logical :: file_exist
        integer :: i, r

        if (cfl_dt) then
            write (t_step_dir, '(A,I0,A,I0)') trim(case_dir) // '/p_all/p', proc_rank, '/', n_start
        else
            write (t_step_dir, '(A,I0,A,I0)') trim(case_dir) // '/p_all/p', proc_rank, '/', t_step_start
        end if

        file_path = trim(t_step_dir) // '/.'
        call my_inquire(file_path, file_exist)

        if (file_exist .neqv. .true.) then
            call s_mpi_abort(trim(file_path) // ' is missing. Exiting.')
        end if

        if (bc_io) then
            call s_read_serial_boundary_condition_files(t_step_dir, bc_type)
        else
            call s_assign_default_bc_type(bc_type)
        end if

        file_path = trim(t_step_dir) // '/x_cb.dat'

        inquire (FILE=trim(file_path), EXIST=file_exist)

        if (file_exist) then
            open (2, FILE=trim(file_path), form='unformatted', ACTION='read', STATUS='old')
            read (2) x_cb(-1:m); close (2)
        else
            call s_mpi_abort(trim(file_path) // ' is missing. Exiting.')
        end if

        dx(0:m) = x_cb(0:m) - x_cb(-1:m - 1)
        x_cc(0:m) = x_cb(-1:m - 1) + dx(0:m)/2._wp

        if (ib) then
            do i = 1, num_ibs
                if (patch_ib(i)%c > 0) then
                    Np = int((patch_ib(i)%p*patch_ib(i)%c/dx(0))*20) + int(((patch_ib(i)%c - patch_ib(i)%p*patch_ib(i)%c)/dx(0)) &
                             & *20) + 1
                end if
            end do
        end if

        if (n > 0) then
            file_path = trim(t_step_dir) // '/y_cb.dat'

            inquire (FILE=trim(file_path), EXIST=file_exist)

            if (file_exist) then
                open (2, FILE=trim(file_path), form='unformatted', ACTION='read', STATUS='old')
                read (2) y_cb(-1:n); close (2)
            else
                call s_mpi_abort(trim(file_path) // ' is missing. Exiting.')
            end if

            dy(0:n) = y_cb(0:n) - y_cb(-1:n - 1)
            y_cc(0:n) = y_cb(-1:n - 1) + dy(0:n)/2._wp
        end if

        if (p > 0) then
            file_path = trim(t_step_dir) // '/z_cb.dat'

            inquire (FILE=trim(file_path), EXIST=file_exist)

            if (file_exist) then
                open (2, FILE=trim(file_path), form='unformatted', ACTION='read', STATUS='old')
                read (2) z_cb(-1:p); close (2)
            else
                call s_mpi_abort(trim(file_path) // ' is missing. Exiting.')
            end if

            dz(0:p) = z_cb(0:p) - z_cb(-1:p - 1)
            z_cc(0:p) = z_cb(-1:p - 1) + dz(0:p)/2._wp
        end if

        do i = 1, sys_size
            write (file_path, '(A,I0,A)') trim(t_step_dir) // '/q_cons_vf', i, '.dat'
            inquire (FILE=trim(file_path), EXIST=file_exist)
            if (file_exist) then
                open (2, FILE=trim(file_path), form='unformatted', ACTION='read', STATUS='old')
                read (2) q_cons_vf(i)%sf(0:m,0:n,0:p); close (2)
            else
                call s_mpi_abort(trim(file_path) // ' is missing. Exiting.')
            end if
        end do

        if (bubbles_euler .or. elasticity) then
            ! Read pb and mv for non-polytropic qbmm
            if (qbmm .and. .not. polytropic) then
                do i = 1, nb
                    do r = 1, nnode
                        write (file_path, '(A,I0,A)') trim(t_step_dir) // '/pb', sys_size + (i - 1)*nnode + r, '.dat'
                        inquire (FILE=trim(file_path), EXIST=file_exist)
                        if (file_exist) then
                            open (2, FILE=trim(file_path), form='unformatted', ACTION='read', STATUS='old')
                            read (2) pb_ts(1)%sf(0:m,0:n,0:p,r, i); close (2)
                        else
                            call s_mpi_abort(trim(file_path) // ' is missing. Exiting.')
                        end if
                    end do
                end do
                do i = 1, nb
                    do r = 1, nnode
                        write (file_path, '(A,I0,A)') trim(t_step_dir) // '/mv', sys_size + (i - 1)*nnode + r, '.dat'
                        inquire (FILE=trim(file_path), EXIST=file_exist)
                        if (file_exist) then
                            open (2, FILE=trim(file_path), form='unformatted', ACTION='read', STATUS='old')
                            read (2) mv_ts(1)%sf(0:m,0:n,0:p,r, i); close (2)
                        else
                            call s_mpi_abort(trim(file_path) // ' is missing. Exiting.')
                        end if
                    end do
                end do
            end if
        end if

    end subroutine s_read_serial_data_files

    !> Read parallel initial condition and grid data files via MPI I/O
    impure subroutine s_read_parallel_data_files(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf

#ifdef MFC_MPI
        real(wp), allocatable, dimension(:)  :: x_cb_glb, y_cb_glb, z_cb_glb
        integer                              :: ifile, ierr, data_size
        integer, dimension(MPI_STATUS_SIZE)  :: status
        integer(KIND=MPI_OFFSET_KIND)        :: disp
        integer(KIND=MPI_OFFSET_KIND)        :: m_MOK, n_MOK, p_MOK
        integer(KIND=MPI_OFFSET_KIND)        :: WP_MOK, var_MOK
        integer(KIND=MPI_OFFSET_KIND)        :: MOK
        character(LEN=path_len + 2*name_len) :: file_loc
        logical                              :: file_exist
        character(len=10)                    :: t_step_start_string
        integer                              :: i, j

        ! Downsampled data variables
        integer :: m_ds, n_ds, p_ds
        integer :: m_glb_ds, n_glb_ds, p_glb_ds
        integer :: m_glb_read, n_glb_read, p_glb_read  !< data size of read

        allocate (x_cb_glb(-1:m_glb))
        allocate (y_cb_glb(-1:n_glb))
        allocate (z_cb_glb(-1:p_glb))

        file_loc = trim(case_dir) // '/restart_data' // trim(mpiiofs) // 'x_cb.dat'
        inquire (FILE=trim(file_loc), EXIST=file_exist)

        if (down_sample) then
            m_ds = int((m + 1)/3) - 1
            n_ds = int((n + 1)/3) - 1
            p_ds = int((p + 1)/3) - 1

            m_glb_ds = int((m_glb + 1)/3) - 1
            n_glb_ds = int((n_glb + 1)/3) - 1
            p_glb_ds = int((p_glb + 1)/3) - 1
        end if

        if (file_exist) then
            data_size = m_glb + 2
            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
            call MPI_FILE_READ(ifile, x_cb_glb, data_size, mpi_p, status, ierr)
            call MPI_FILE_CLOSE(ifile, ierr)
        else
            call s_mpi_abort('File ' // trim(file_loc) // ' is missing. Exiting.')
        end if

        x_cb(-1:m) = x_cb_glb((start_idx(1) - 1):(start_idx(1) + m))
        dx(0:m) = x_cb(0:m) - x_cb(-1:m - 1)
        x_cc(0:m) = x_cb(-1:m - 1) + dx(0:m)/2._wp

        if (ib) then
            do i = 1, num_ibs
                if (patch_ib(i)%c > 0) then
                    Np = int((patch_ib(i)%p*patch_ib(i)%c/dx(0))*20) + int(((patch_ib(i)%c - patch_ib(i)%p*patch_ib(i)%c)/dx(0)) &
                             & *20) + 1
                    allocate (MPI_IO_airfoil_IB_DATA%var(1:2*Np))
                end if
            end do
        end if

        if (n > 0) then
            file_loc = trim(case_dir) // '/restart_data' // trim(mpiiofs) // 'y_cb.dat'
            inquire (FILE=trim(file_loc), EXIST=file_exist)

            if (file_exist) then
                data_size = n_glb + 2
                call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
                call MPI_FILE_READ(ifile, y_cb_glb, data_size, mpi_p, status, ierr)
                call MPI_FILE_CLOSE(ifile, ierr)
            else
                call s_mpi_abort('File ' // trim(file_loc) // ' is missing. Exiting.')
            end if

            y_cb(-1:n) = y_cb_glb((start_idx(2) - 1):(start_idx(2) + n))
            dy(0:n) = y_cb(0:n) - y_cb(-1:n - 1)
            y_cc(0:n) = y_cb(-1:n - 1) + dy(0:n)/2._wp

            if (p > 0) then
                file_loc = trim(case_dir) // '/restart_data' // trim(mpiiofs) // 'z_cb.dat'
                inquire (FILE=trim(file_loc), EXIST=file_exist)

                if (file_exist) then
                    data_size = p_glb + 2
                    call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
                    call MPI_FILE_READ(ifile, z_cb_glb, data_size, mpi_p, status, ierr)
                    call MPI_FILE_CLOSE(ifile, ierr)
                else
                    call s_mpi_abort('File ' // trim(file_loc) // 'is missing. Exiting.')
                end if

                z_cb(-1:p) = z_cb_glb((start_idx(3) - 1):(start_idx(3) + p))
                dz(0:p) = z_cb(0:p) - z_cb(-1:p - 1)
                z_cc(0:p) = z_cb(-1:p - 1) + dz(0:p)/2._wp
            end if
        end if

        if (file_per_process) then
            if (cfl_dt) then
                call s_int_to_str(n_start, t_step_start_string)
                write (file_loc, '(I0,A1,I7.7,A)') n_start, '_', proc_rank, '.dat'
            else
                call s_int_to_str(t_step_start, t_step_start_string)
                write (file_loc, '(I0,A1,I7.7,A)') t_step_start, '_', proc_rank, '.dat'
            end if
            file_loc = trim(case_dir) // '/restart_data/lustre_' // trim(t_step_start_string) // trim(mpiiofs) // trim(file_loc)
            inquire (FILE=trim(file_loc), EXIST=file_exist)

            if (file_exist) then
                call MPI_FILE_OPEN(MPI_COMM_SELF, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

                if (down_sample) then
                    call s_initialize_mpi_data_ds(q_cons_vf)
                else
                    if (ib) then
                        call s_initialize_mpi_data(q_cons_vf, ib_markers)
                    else
                        call s_initialize_mpi_data(q_cons_vf)
                    end if
                end if

                if (down_sample) then
                    data_size = (m_ds + 3)*(n_ds + 3)*(p_ds + 3)
                    m_glb_read = m_glb_ds + 1
                    n_glb_read = n_glb_ds + 1
                    p_glb_read = p_glb_ds + 1
                else
                    data_size = (m + 1)*(n + 1)*(p + 1)
                    m_glb_read = m_glb + 1
                    n_glb_read = n_glb + 1
                    p_glb_read = p_glb + 1
                end if

                m_MOK = int(m_glb_read + 1, MPI_OFFSET_KIND)
                n_MOK = int(m_glb_read + 1, MPI_OFFSET_KIND)
                p_MOK = int(m_glb_read + 1, MPI_OFFSET_KIND)
                WP_MOK = int(storage_size(0._stp)/8, MPI_OFFSET_KIND)
                MOK = int(1._wp, MPI_OFFSET_KIND)

                if (bubbles_euler .or. elasticity) then
                    do i = 1, sys_size
                        var_MOK = int(i, MPI_OFFSET_KIND)

                        call MPI_FILE_READ(ifile, MPI_IO_DATA%var(i)%sf, data_size*mpi_io_type, mpi_io_p, status, ierr)
                    end do
                    ! Read pb and mv for non-polytropic qbmm
                    if (qbmm .and. .not. polytropic) then
                        do i = sys_size + 1, sys_size + 2*nb*nnode
                            var_MOK = int(i, MPI_OFFSET_KIND)

                            call MPI_FILE_READ(ifile, MPI_IO_DATA%var(i)%sf, data_size*mpi_io_type, mpi_io_p, status, ierr)
                        end do
                    end if
                else
                    if (down_sample) then
                        do i = 1, sys_size
                            var_MOK = int(i, MPI_OFFSET_KIND)

                            call MPI_FILE_READ(ifile, q_cons_temp(i)%sf, data_size*mpi_io_type, mpi_io_p, status, ierr)
                        end do
                    else
                        do i = 1, sys_size
                            var_MOK = int(i, MPI_OFFSET_KIND)

                            call MPI_FILE_READ(ifile, MPI_IO_DATA%var(i)%sf, data_size*mpi_io_type, mpi_io_p, status, ierr)
                        end do
                    end if
                end if

                call s_mpi_barrier()

                call MPI_FILE_CLOSE(ifile, ierr)
            else
                call s_mpi_abort('File ' // trim(file_loc) // ' is missing. Exiting.')
            end if
        else
            if (cfl_dt) then
                write (file_loc, '(I0,A)') n_start, '.dat'
            else
                write (file_loc, '(I0,A)') t_step_start, '.dat'
            end if
            file_loc = trim(case_dir) // '/restart_data' // trim(mpiiofs) // trim(file_loc)
            inquire (FILE=trim(file_loc), EXIST=file_exist)

            if (file_exist) then
                call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

                if (ib) then
                    call s_initialize_mpi_data(q_cons_vf, ib_markers)
                else
                    call s_initialize_mpi_data(q_cons_vf)
                end if

                data_size = (m + 1)*(n + 1)*(p + 1)

                m_MOK = int(m_glb + 1, MPI_OFFSET_KIND)
                n_MOK = int(n_glb + 1, MPI_OFFSET_KIND)
                p_MOK = int(p_glb + 1, MPI_OFFSET_KIND)
                WP_MOK = int(storage_size(0._stp)/8, MPI_OFFSET_KIND)
                MOK = int(1._wp, MPI_OFFSET_KIND)

                if (bubbles_euler .or. elasticity) then
                    do i = 1, sys_size
                        var_MOK = int(i, MPI_OFFSET_KIND)
                        disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                        call MPI_FILE_SET_VIEW(ifile, disp, mpi_io_p, MPI_IO_DATA%view(i), 'native', mpi_info_int, ierr)
                        call MPI_FILE_READ(ifile, MPI_IO_DATA%var(i)%sf, data_size*mpi_io_type, mpi_io_p, status, ierr)
                    end do
                    ! Read pb and mv for non-polytropic qbmm
                    if (qbmm .and. .not. polytropic) then
                        do i = sys_size + 1, sys_size + 2*nb*nnode
                            var_MOK = int(i, MPI_OFFSET_KIND)
                            disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                            call MPI_FILE_SET_VIEW(ifile, disp, mpi_io_p, MPI_IO_DATA%view(i), 'native', mpi_info_int, ierr)
                            call MPI_FILE_READ(ifile, MPI_IO_DATA%var(i)%sf, data_size*mpi_io_type, mpi_io_p, status, ierr)
                        end do
                    end if
                else
                    do i = 1, sys_size
                        var_MOK = int(i, MPI_OFFSET_KIND)

                        disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                        call MPI_FILE_SET_VIEW(ifile, disp, mpi_io_p, MPI_IO_DATA%view(i), 'native', mpi_info_int, ierr)
                        call MPI_FILE_READ_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size*mpi_io_type, mpi_io_p, status, ierr)
                    end do
                end if

                call s_mpi_barrier()

                call MPI_FILE_CLOSE(ifile, ierr)
            else
                call s_mpi_abort('File ' // trim(file_loc) // ' is missing. Exiting.')
            end if
        end if

        deallocate (x_cb_glb, y_cb_glb, z_cb_glb)

        if (bc_io) then
            call s_read_parallel_boundary_condition_files(bc_type)
        else
            call s_assign_default_bc_type(bc_type)
        end if
#endif

    end subroutine s_read_parallel_data_files

    !> Initialize internal-energy equations from phase mass, mixture momentum, and total energy
    subroutine s_initialize_internal_energy_equations(v_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: v_vf
        real(wp)                                               :: rho
        real(wp)                                               :: dyn_pres
        real(wp)                                               :: gamma
        real(wp)                                               :: pi_inf
        real(wp)                                               :: qv
        real(wp), dimension(2)                                 :: Re
        real(wp)                                               :: pres, T
        integer                                                :: i, j, k, l, c
        real(wp), dimension(num_species)                       :: rhoYks
        real(wp)                                               :: pres_mag

        pres_mag = 0._wp

        T = dflt_T_guess

        do j = 0, m
            do k = 0, n
                do l = 0, p
                    call s_convert_to_mixture_variables(v_vf, j, k, l, rho, gamma, pi_inf, qv, Re)

                    dyn_pres = 0._wp
                    do i = eqn_idx%mom%beg, eqn_idx%mom%end
                        dyn_pres = dyn_pres + 5.e-1_wp*v_vf(i)%sf(j, k, l)*v_vf(i)%sf(j, k, l)/max(rho, sgm_eps)
                    end do

                    if (chemistry) then
                        do c = 1, num_species
                            rhoYks(c) = v_vf(eqn_idx%species%beg + c - 1)%sf(j, k, l)
                        end do
                    end if

                    if (mhd) then
                        if (n == 0) then
                            pres_mag = 0.5_wp*(Bx0**2 + v_vf(eqn_idx%B%beg)%sf(j, k, l)**2 + v_vf(eqn_idx%B%beg + 1)%sf(j, k, l)**2)
                        else
                            pres_mag = 0.5_wp*(v_vf(eqn_idx%B%beg)%sf(j, k, l)**2 + v_vf(eqn_idx%B%beg + 1)%sf(j, k, &
                                               & l)**2 + v_vf(eqn_idx%B%beg + 2)%sf(j, k, l)**2)
                        end if
                    end if

                    call s_compute_pressure(v_vf(eqn_idx%E)%sf(j, k, l), 0._stp, dyn_pres, pi_inf, gamma, rho, qv, rhoYks, pres, &
                                            & T, pres_mag=pres_mag)

                    do i = 1, num_fluids
                        v_vf(i + eqn_idx%int_en%beg - 1)%sf(j, k, l) = v_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, &
                             & l)*(gammas(i)*pres + pi_infs(i)) + v_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l)*qvs(i)
                    end do
                end do
            end do
        end do

    end subroutine s_initialize_internal_energy_equations

    !> Advance the simulation by one time step, handling CFL-based dt and time-stepper dispatch
    impure subroutine s_perform_time_step(t_step, time_avg)

        integer, intent(inout)  :: t_step
        real(wp), intent(inout) :: time_avg
        integer                 :: i, eta_hh, eta_mm, eta_ss
        real(wp)                :: eta_sec

        if (cfl_dt) then
            if (cfl_const_dt .and. t_step == 0) call s_compute_dt()

            if (cfl_adap_dt) call s_compute_dt()

            if (t_step == 0) dt_init = dt

            if (dt < 1.e-3_wp*dt_init .and. cfl_adap_dt .and. proc_rank == 0) then
                print *, "Delta t = ", dt
                call s_mpi_abort("Delta t has become too small")
            end if
        end if

        if (cfl_dt) then
            if ((mytime + dt) >= t_stop) then
                dt = t_stop - mytime
                $:GPU_UPDATE(device='[dt]')
            end if
        else
            if ((mytime + dt) >= finaltime) then
                dt = finaltime - mytime
                $:GPU_UPDATE(device='[dt]')
            end if
        end if

        if (cfl_dt) then
            if (proc_rank == 0 .and. mod(t_step - t_step_start, t_step_print) == 0) then
                eta_sec = wall_time_avg*(t_stop - mytime)/max(dt, tiny(dt))
                eta_hh = int(eta_sec)/3600
                eta_mm = mod(int(eta_sec), 3600)/60
                eta_ss = mod(int(eta_sec), 60)
                print '(" [", I3, "%] Time ", ES16.6, " dt = ", ES16.6, " @ Time Step = ", I8,  " Time Avg = ", ES16.6,  " Time/step = ", ES12.6, " ETA (HH:MM:SS) = ", I0, ":", I2.2, ":", I2.2)', &
                    & int(ceiling(100._wp*(mytime/t_stop))), mytime, dt, t_step, wall_time_avg, wall_time, eta_hh, eta_mm, eta_ss
            end if
        else
            if (proc_rank == 0 .and. mod(t_step - t_step_start, t_step_print) == 0) then
                eta_sec = wall_time_avg*real(t_step_stop - t_step, wp)
                eta_hh = int(eta_sec)/3600
                eta_mm = mod(int(eta_sec), 3600)/60
                eta_ss = mod(int(eta_sec), 60)
                print '(" [", I3, "%]  Time step ", I8, " of ", I0, " @ t_step = ", I8,  " Time Avg = ", ES12.6,  " Time/step= ", ES12.6, " ETA (HH:MM:SS) = ", I0, ":", I2.2, ":", I2.2)', &
                    & int(ceiling(100._wp*(real(t_step - t_step_start)/(t_step_stop - t_step_start + 1)))), &
                    & t_step - t_step_start + 1, t_step_stop - t_step_start + 1, t_step, wall_time_avg, wall_time, eta_hh, &
                    & eta_mm, eta_ss
            end if
        end if

        if (probe_wrt) then
            do i = 1, sys_size
                $:GPU_UPDATE(host='[q_cons_ts(1)%vf(i)%sf]')
            end do
        end if

        ! Total-variation-diminishing (TVD) Runge-Kutta (RK) time-steppers
        if (any(time_stepper == (/1, 2, 3/))) then
            call s_tvd_rk(t_step, time_avg, time_stepper)
        end if

        ! Advance time after RK so source terms see current-step time
        mytime = mytime + dt

        if (relax) call s_infinite_relaxation_k(q_cons_ts(1)%vf)

        ! Time-stepping loop controls
        t_step = t_step + 1

    end subroutine s_perform_time_step

    !> Collect per-process wall-clock times and write aggregate performance metrics to file
    impure subroutine s_save_performance_metrics(time_avg, time_final, io_time_avg, io_time_final, proc_time, io_proc_time, &

        & file_exists)

        real(wp), intent(inout)               :: time_avg, time_final
        real(wp), intent(inout)               :: io_time_avg, io_time_final
        real(wp), dimension(:), intent(inout) :: proc_time
        real(wp), dimension(:), intent(inout) :: io_proc_time
        logical, intent(inout)                :: file_exists
        real(wp)                              :: grind_time

        call s_mpi_barrier()

        if (num_procs > 1) then
            call mpi_bcast_time_step_values(proc_time, time_avg)

            call mpi_bcast_time_step_values(io_proc_time, io_time_avg)
        end if

        if (proc_rank == 0) then
            time_final = 0._wp
            io_time_final = 0._wp
            if (num_procs == 1) then
                time_final = time_avg
                io_time_final = io_time_avg
            else
                time_final = maxval(proc_time)
                io_time_final = maxval(io_proc_time)
            end if

            grind_time = time_final*1.0e9_wp/(real(sys_size, wp)*real(maxval((/1, m_glb/)), wp)*real(maxval((/1, n_glb/)), &
                                              & wp)*real(maxval((/1, p_glb/)), wp))

            print *, "Performance:", grind_time, "ns/gp/eq/rhs"
            inquire (FILE='time_data.dat', EXIST=file_exists)
            if (file_exists) then
                open (1, file='time_data.dat', position='append', status='old')
            else
                open (1, file='time_data.dat', status='new')
                write (1, '(A10, A15, A15)') "Ranks", "s/step", "ns/gp/eq/rhs"
            end if

            write (1, '(I10, 2(F15.8))') num_procs, time_final, grind_time

            close (1)

            inquire (FILE='io_time_data.dat', EXIST=file_exists)
            if (file_exists) then
                open (1, file='io_time_data.dat', position='append', status='old')
            else
                open (1, file='io_time_data.dat', status='new')
                write (1, '(A10, A15)') "Ranks", "s/step"
            end if

            write (1, '(I10, F15.8)') num_procs, io_time_final
            close (1)
        end if

    end subroutine s_save_performance_metrics

    !> Save conservative variable data to disk at the current time step
    impure subroutine s_save_data(t_step, start, finish, io_time_avg, nt)

        integer, intent(inout)  :: t_step
        real(wp), intent(inout) :: start, finish, io_time_avg
        integer, intent(inout)  :: nt
        integer(kind=8)         :: i, j, k, l
        integer                 :: stor
        integer                 :: save_count

        if (down_sample) then
            call s_populate_variables_buffers(bc_type, q_cons_ts(1)%vf)
        end if

        stor = 1

        if (time_stepper /= 1) then
            $:GPU_PARALLEL_LOOP(collapse=4, copyin='[idwbuff]')
            do i = 1, sys_size
                do l = idwbuff(3)%beg, idwbuff(3)%end
                    do k = idwbuff(2)%beg, idwbuff(2)%end
                        do j = idwbuff(1)%beg, idwbuff(1)%end
                            q_cons_ts(2)%vf(i)%sf(j, k, l) = q_cons_ts(1)%vf(i)%sf(j, k, l)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            stor = 2
        end if

        call cpu_time(start)
        call nvtxStartRange("SAVE-DATA")
        do i = 1, sys_size
#ifndef FRONTIER_UNIFIED
            $:GPU_UPDATE(host='[q_cons_ts(stor)%vf(i)%sf]')
#endif
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        if (ieee_is_nan(real(q_cons_ts(stor)%vf(i)%sf(j, k, l), kind=wp))) then
                            print *, "NaN(s) in timestep output.", j, k, l, i, proc_rank, t_step, m, n, p
                            call s_mpi_abort("NaN(s) in timestep output.")
                        end if
                    end do
                end do
            end do
        end do

        if (qbmm .and. .not. polytropic) then
            $:GPU_UPDATE(host='[pb_ts(1)%sf]')
            $:GPU_UPDATE(host='[mv_ts(1)%sf]')
        end if

        if (cfl_dt) then
            save_count = int(mytime/t_save)
        else
            save_count = t_step
        end if

        if (bubbles_lagrange) then
            $:GPU_UPDATE(host='[lag_id, mtn_pos, mtn_posPrev, mtn_vel, intfc_rad, intfc_vel, bub_R0, Rmax_stats, Rmin_stats, &
                         & bub_dphidt, gas_p, gas_mv, gas_mg, gas_betaT, gas_betaC]')
            do i = 1, nBubs
                if (ieee_is_nan(intfc_rad(i, 1)) .or. intfc_rad(i, 1) <= 0._wp) then
                    call s_mpi_abort("Bubble radius is negative or NaN, please reduce dt.")
                end if
            end do

            $:GPU_UPDATE(host='[q_beta(1)%sf]')
            call s_write_data_files(q_cons_ts(stor)%vf, q_T_sf, q_prim_vf, save_count, bc_type, q_beta(1))
            $:GPU_UPDATE(host='[Rmax_stats, Rmin_stats, gas_p, gas_mv, intfc_vel]')
            call s_write_restart_lag_bubbles(save_count)  ! parallel
            if (lag_params%write_bubbles_stats) call s_write_lag_bubble_stats()
        else
            call s_write_data_files(q_cons_ts(stor)%vf, q_T_sf, q_prim_vf, save_count, bc_type)
        end if

        ! Write IB kinematic state for restart
        if (ib) call s_write_ib_state_file(t_step)

        call nvtxEndRange
        call cpu_time(finish)
        if (cfl_dt) then
            nt = mytime/t_save
        else
            nt = int((t_step - t_step_start)/(t_step_save))
        end if

        if (nt == 1) then
            io_time_avg = abs(finish - start)
        else
            io_time_avg = (abs(finish - start) + io_time_avg*(nt - 1))/nt
        end if

    end subroutine s_save_data

    !> Initialize all simulation sub-modules in the required dependency order
    impure subroutine s_initialize_modules

        integer :: m_ds, n_ds, p_ds
        integer :: i

        call s_initialize_global_parameters_module()
        #:if USING_AMD
            #:for BC in {-5, -6, -7, -8, -9, -10, -11, -12, -13}
                @:PROHIBIT(any((/bc_x%beg, bc_x%end, bc_y%beg, bc_y%end, bc_z%beg, &
                           & bc_z%end/) == ${BC}$) .and. eqn_idx%adv%end > 20 .and. (.not. chemistry), &
                           & "CBC module with AMD compiler requires eqn_idx%adv%end <= 20 when case optimization is turned off")
                @:PROHIBIT(any((/bc_x%beg, bc_x%end, bc_y%beg, bc_y%end, bc_z%beg, &
                           & bc_z%end/) == ${BC}$) .and. sys_size > 20 .and. (chemistry), &
                           & "CBC module with AMD compiler and chemistry requires sys_size <= 20 when case optimization is turned off")
            #:endfor
        #:endif
        if (bubbles_euler .or. bubbles_lagrange) then
            call s_initialize_bubbles_model()
        end if
        call s_initialize_mpi_common_module()
        call s_initialize_mpi_proxy_module()
        call s_initialize_variables_conversion_module()
        if (grid_geometry == 3) call s_initialize_fftw_module()

        if (bubbles_euler) call s_initialize_bubbles_EE_module()
        if (ib) call s_initialize_ibm_module()
        if (qbmm) call s_initialize_qbmm_module()

        if (acoustic_source) then
            call s_initialize_acoustic_src()
        end if

        if (viscous .and. (.not. igr)) then
            call s_initialize_viscous_module()
        end if

        call s_initialize_rhs_module()

        if (surface_tension) call s_initialize_surface_tension_module()

        if (relax) call s_initialize_phasechange_module()

        call s_initialize_data_output_module()
        call s_initialize_derived_variables_module()
        call s_initialize_time_steppers_module()

        call s_initialize_boundary_common_module()

        if (down_sample) then
            m_ds = int((m + 1)/3) - 1
            n_ds = int((n + 1)/3) - 1
            p_ds = int((p + 1)/3) - 1

            allocate (q_cons_temp(1:sys_size))
            do i = 1, sys_size
                allocate (q_cons_temp(i)%sf(-1:m_ds + 1,-1:n_ds + 1,-1:p_ds + 1))
            end do
        end if

        if (down_sample) then
            call s_read_data_files(q_cons_temp)
            call s_upsample_data(q_cons_ts(1)%vf, q_cons_temp)
            do i = 1, sys_size
                $:GPU_UPDATE(device='[q_cons_ts(1)%vf(i)%sf]')
            end do
            do i = 1, sys_size
                deallocate (q_cons_temp(i)%sf)
            end do
            deallocate (q_cons_temp)
        else
            call s_read_data_files(q_cons_ts(1)%vf)
        end if

        call s_populate_grid_variables_buffers()

        if (model_eqns == 3) call s_initialize_internal_energy_equations(q_cons_ts(1)%vf)
        if (ib) then
            if (t_step_start > 0) call s_read_ib_restart_data(t_step_start)
            call s_ibm_setup()
            if (t_step_start == 0) then
                call s_write_ib_data_file(0)
                call s_write_ib_state_file(0)
            end if
        end if
        if (bodyForces) call s_initialize_body_forces_module()
        if (acoustic_source) call s_precalculate_acoustic_spatial_sources()

        ! Initialize the Temperature cache.
        if (chemistry) call s_compute_q_T_sf(q_T_sf, q_cons_ts(1)%vf, idwint)

        ! Computation of parameters, allocation of memory, association of pointers, and/or execution of any other tasks that are
        ! needed to properly configure the modules. The preparations below DO DEPEND on the grid being complete.
        if (igr .or. dummy) then
            call s_initialize_igr_module()
        end if
        if (.not. igr .or. dummy) then
            if (recon_type == WENO_TYPE) then
                call s_initialize_weno_module()
            else if (recon_type == MUSCL_TYPE) then
                call s_initialize_muscl_module()
            end if
            call s_initialize_cbc_module()
            call s_initialize_riemann_solvers_module()
        end if

        call s_initialize_derived_variables()
        if (bubbles_lagrange) call s_initialize_bubbles_EL_module(q_cons_ts(1)%vf)

        if (hypoelasticity) call s_initialize_hypoelastic_module()
        if (hyperelasticity) call s_initialize_hyperelastic_module()

    end subroutine s_initialize_modules

    !> Set up the MPI execution environment, bind GPUs, and decompose the computational domain
    impure subroutine s_initialize_mpi_domain

        integer :: ierr

#ifdef MFC_GPU
        real(wp) :: starttime, endtime
        integer  :: num_devices, local_size, num_nodes, ppn, my_device_num
        integer  :: dev, devNum, local_rank
#ifdef MFC_MPI
        integer :: local_comm
#endif
#if defined(MFC_OpenACC)
        integer(acc_device_kind) :: devtype
#endif
#endif

        call s_mpi_initialize()

#ifdef MFC_GPU
#ifndef MFC_MPI
        local_size = 1
        local_rank = 0
#else
        call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, local_comm, ierr)
        call MPI_Comm_size(local_comm, local_size, ierr)
        call MPI_Comm_rank(local_comm, local_rank, ierr)
#endif
#if defined(MFC_OpenACC)
        devtype = acc_get_device_type()
        devNum = acc_get_num_devices(devtype)
        dev = mod(local_rank, devNum)

        call acc_set_device_num(dev, devtype)
#elif defined(MFC_OpenMP)
        devNum = omp_get_num_devices()
        dev = mod(local_rank, devNum)
        call omp_set_default_device(dev)
#endif
#endif

        if (proc_rank == 0) then
            call s_assign_default_values_to_user_inputs()
            call s_read_input_file()
            call s_check_input_file()

            print '(" Simulating a ", A, " ", I0, "x", I0, "x", I0, " case on ", I0, " rank(s) ", A, ".")', &
            #:if not MFC_CASE_OPTIMIZATION
                "regular", &
            #:else
                "case-optimized", &
            #:endif
            m, n, p, num_procs, &
#if defined(MFC_OpenACC)
            "with OpenACC offloading"
#elif defined(MFC_OpenMP)
            "with OpenMP offloading"
#else
            "on CPUs"
#endif
        end if

        call s_mpi_bcast_user_inputs()

        ! Save original BCs before decomposition overwrites them with MPI neighbor ranks
        ib_bc_x = bc_x
        ib_bc_y = bc_y
        ib_bc_z = bc_z

        call s_initialize_parallel_io()

        call s_mpi_decompose_computational_domain()

    end subroutine s_initialize_mpi_domain

    !> Transfer initial conservative variable and model parameter data to the GPU device
    subroutine s_initialize_gpu_vars

        integer :: i

        if (.not. down_sample) then
            do i = 1, sys_size
                $:GPU_UPDATE(device='[q_cons_ts(1)%vf(i)%sf]')
            end do
        end if

        if (qbmm .and. .not. polytropic) then
            $:GPU_UPDATE(device='[pb_ts(1)%sf, mv_ts(1)%sf]')
        end if
        if (chemistry) then
            $:GPU_UPDATE(device='[q_T_sf%sf]')
        end if

        $:GPU_UPDATE(device='[chem_params]')

        $:GPU_UPDATE(device='[R0ref, p0ref, rho0ref, ss, pv, vd, mu_l, mu_v, mu_g, gam_v, gam_g, M_v, M_g, R_v, R_g, Tw, cp_v, &
                     & cp_g, k_vl, k_gl, gam, gam_m, Eu, Ca, Web, Re_inv, Pe_c, phi_vg, phi_gv, omegaN, bubbles_euler, &
                     & polytropic, polydisperse, qbmm, ptil, bubble_model, thermal, poly_sigma, adv_n, adap_dt, adap_dt_tol, &
                     & adap_dt_max_iters, eqn_idx%n, pi_fac, low_Mach]')

        if (bubbles_euler) then
            $:GPU_UPDATE(device='[weight, R0]')
            if (.not. polytropic) then
                $:GPU_UPDATE(device='[pb0, Pe_T, k_g, k_v, mass_g0, mass_v0, Re_trans_T, Re_trans_c, Im_trans_T, Im_trans_c]')
            else if (qbmm) then
                $:GPU_UPDATE(device='[pb0]')
            end if
        end if

        $:GPU_UPDATE(device='[adv_n, adap_dt, adap_dt_tol, adap_dt_max_iters, pi_fac, low_Mach]')

        $:GPU_UPDATE(device='[acoustic_source, num_source]')
        $:GPU_UPDATE(device='[sigma, surface_tension]')

        $:GPU_UPDATE(device='[dx, dy, dz, x_cb, x_cc, y_cb, y_cc, z_cb, z_cc]')
        $:GPU_UPDATE(device='[bc_x%beg, bc_x%end, bc_y%beg, bc_y%end, bc_z%beg, bc_z%end]')
        $:GPU_UPDATE(device='[bc_x%vb1, bc_x%vb2, bc_x%vb3, bc_x%ve1, bc_x%ve2, bc_x%ve3]')
        $:GPU_UPDATE(device='[bc_y%vb1, bc_y%vb2, bc_y%vb3, bc_y%ve1, bc_y%ve2, bc_y%ve3]')
        $:GPU_UPDATE(device='[bc_z%vb1, bc_z%vb2, bc_z%vb3, bc_z%ve1, bc_z%ve2, bc_z%ve3]')

        $:GPU_UPDATE(device='[bc_x%grcbc_in, bc_x%grcbc_out, bc_x%grcbc_vel_out]')
        $:GPU_UPDATE(device='[bc_y%grcbc_in, bc_y%grcbc_out, bc_y%grcbc_vel_out]')
        $:GPU_UPDATE(device='[bc_z%grcbc_in, bc_z%grcbc_out, bc_z%grcbc_vel_out]')

        $:GPU_UPDATE(device='[bc_x%isothermal_in, bc_x%isothermal_out]')
        $:GPU_UPDATE(device='[bc_y%isothermal_in, bc_y%isothermal_out]')
        $:GPU_UPDATE(device='[bc_z%isothermal_in, bc_z%isothermal_out]')
        $:GPU_UPDATE(device='[bc_x%Twall_in, bc_x%Twall_out, bc_y%Twall_in, bc_y%Twall_out, bc_z%Twall_in, bc_z%Twall_out]')

        $:GPU_UPDATE(device='[relax, relax_model]')
        if (relax) then
            $:GPU_UPDATE(device='[palpha_eps, ptgalpha_eps]')
        end if

        if (ib) then
            $:GPU_UPDATE(device='[ib_markers%sf]')
        end if
        #:if not MFC_CASE_OPTIMIZATION
            $:GPU_UPDATE(device='[igr, nb, igr_order]')
        #:endif
        #:if USING_AMD
            block
                use m_thermochem, only: molecular_weights
                use m_chemistry, only: molecular_weights_nonparameter
                molecular_weights_nonparameter(:) = molecular_weights(:)
                $:GPU_UPDATE(device='[molecular_weights_nonparameter]')
            end block
        #:endif

    end subroutine s_initialize_gpu_vars

    !> Finalize and deallocate all simulation sub-modules in reverse initialization order
    impure subroutine s_finalize_modules

        call s_finalize_time_steppers_module()
        if (hypoelasticity) call s_finalize_hypoelastic_module()
        if (hyperelasticity) call s_finalize_hyperelastic_module()
        call s_finalize_derived_variables_module()
        call s_finalize_data_output_module()
        call s_finalize_rhs_module()
        if (igr) then
            call s_finalize_igr_module()
        else
            call s_finalize_cbc_module()
            call s_finalize_riemann_solvers_module()
            if (recon_type == WENO_TYPE) then
                call s_finalize_weno_module()
            else if (recon_type == MUSCL_TYPE) then
                call s_finalize_muscl_module()
            end if
        end if
        call s_finalize_variables_conversion_module()
        if (grid_geometry == 3) call s_finalize_fftw_module
        call s_finalize_mpi_common_module()
        call s_finalize_global_parameters_module()
        call s_finalize_boundary_common_module()
        if (relax) call s_finalize_relaxation_solver_module()
        if (bubbles_lagrange) call s_finalize_lagrangian_solver()
        if (viscous .and. (.not. igr)) then
            call s_finalize_viscous_module()
        end if
        call s_finalize_mpi_proxy_module()

        if (surface_tension) call s_finalize_surface_tension_module()
        if (bodyForces) call s_finalize_body_forces_module()
        if (ib) call s_finalize_ibm_module()

        call s_mpi_finalize()

    end subroutine s_finalize_modules

    !> @brief Reads IB kinematic state from restart_data/ib_state.dat on restart. Rank 0 reads the last num_ibs records and
    !! broadcasts to all ranks. Overwrites patch_ib vel, angular_vel, angles, and centroid.
    impure subroutine s_read_ib_restart_data(t_step)

        integer, intent(in)                  :: t_step
        character(len=path_len + 2*name_len) :: file_loc
        integer                              :: i, ios, file_unit, ierr
        integer, parameter                   :: NFIELDS_PER_IB = 20
        real(wp)                             :: ib_buf(NFIELDS_PER_IB)
        logical                              :: file_exist

        write (file_loc, '(A,I0,A)') '/restart_data/ib_state_', t_step, '.dat'
        file_loc = trim(case_dir) // trim(file_loc)

        if (proc_rank == 0) then
            inquire (FILE=trim(file_loc), EXIST=file_exist)
            if (.not. file_exist) then
                call s_mpi_abort('Cannot open IB state file for restart: ' // trim(file_loc))
            end if

            open (newunit=file_unit, file=trim(file_loc), form='unformatted', access='stream', status='old', iostat=ios)
            if (ios /= 0) call s_mpi_abort('Error opening IB state restart file: ' // trim(file_loc))

            do i = 1, num_ibs
                read (file_unit, iostat=ios) ib_buf
                if (ios /= 0) call s_mpi_abort('Error reading IB state restart file')

                patch_ib(i)%vel = ib_buf(8:10)
                patch_ib(i)%angular_vel = ib_buf(11:13)
                patch_ib(i)%angles = ib_buf(14:16)
                patch_ib(i)%x_centroid = ib_buf(17)
                patch_ib(i)%y_centroid = ib_buf(18)
                patch_ib(i)%z_centroid = ib_buf(19)
            end do

            close (file_unit)
        end if

#ifdef MFC_MPI
        do i = 1, num_ibs
            call MPI_BCAST(patch_ib(i)%vel, 3, mpi_p, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_ib(i)%angular_vel, 3, mpi_p, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_ib(i)%angles, 3, mpi_p, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_ib(i)%x_centroid, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_ib(i)%y_centroid, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_ib(i)%z_centroid, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
        end do
#endif

    end subroutine s_read_ib_restart_data

end module m_start_up
