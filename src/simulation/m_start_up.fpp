!>
!! @file m_start_up.f90
!! @brief Contains module m_start_up

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief The purpose of the module is primarily to read in the files that
!!              contain the inputs, the initial condition data and the grid data
!!              that are provided by the user. The module is additionally tasked
!!              with verifying the consistency of the user inputs and completing
!!              the grid variablesThe purpose of the module is primarily to read
!!              in the files that
!!              contain the inputs, the initial condition data and the grid data
!!              that are provided by the user. The module is additionally tasked
!!              with verifying the consistency of the user inputs and completing
!!              the grid variables. This module also also allocating, initializing
!!              I/O, and deallocating the relevant variables on both cpus and gpus as well as
!!              setting up the time stepping, domain decomposition and I/O procedures.
module m_start_up

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_mpi_common

    use m_variables_conversion !< State variables type conversion procedures

    use m_weno                 !< Weighted and essentially non-oscillatory (WENO)
                               !! schemes for spatial reconstruction of variables

    use m_muscl                !< Monotonic Upstream-centered (MUSCL)
                               !! schemes for convservation laws

    use m_riemann_solvers      !< Exact and approximate Riemann problem solvers

    use m_cbc                  !< Characteristic boundary conditions (CBC)

    use m_boundary_common

    use m_acoustic_src      !< Acoustic source calculations

    use m_rhs                  !< Right-hane-side (RHS) evaluation procedures

    use m_chemistry            !< Chemistry module

    use m_data_output          !< Run-time info & solution data output procedures

    use m_time_steppers        !< Time-stepping algorithms

    use m_qbmm                 !< Quadrature MOM

    use m_derived_variables     !< Procedures used to compute quantities derived
                                !! from the conservative and primitive variables
    use m_hypoelastic

    use m_hyperelastic

    use m_phase_change          !< Phase-change module

    use m_viscous

    use m_bubbles_EE            !< Ensemble-averaged bubble dynamics routines

    use m_bubbles_EL            !< Lagrange bubble dynamics routines

    use ieee_arithmetic

    use m_helper_basic          !< Functions to compare floating point numbers

#ifdef MFC_OpenACC
    use openacc
#endif

    use m_nvtx

    use m_ibm

    use m_compile_specific

    use m_checker_common

    use m_checker

    use m_surface_tension

    use m_body_forces

    use m_sim_helpers

    use m_mhd

    use m_igr

    implicit none

    private; public :: s_read_input_file, &
 s_check_input_file, &
 s_read_data_files, &
 s_read_serial_data_files, &
 s_read_parallel_data_files, &
 s_initialize_internal_energy_equations, &
 s_initialize_modules, s_initialize_gpu_vars, &
 s_initialize_mpi_domain, s_finalize_modules, &
 s_perform_time_step, s_save_data, &
 s_save_performance_metrics

    type(scalar_field), allocatable, dimension(:) :: q_cons_temp
    $:GPU_DECLARE(create='[q_cons_temp]')

    real(wp) :: dt_init

contains

   !> Read data files. Dispatch subroutine that replaces procedure pointer.
        !! @param q_cons_vf Conservative variables
    impure subroutine s_read_data_files(q_cons_vf)

        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: q_cons_vf


        if (.not. parallel_io) then
            call s_read_serial_data_files(q_cons_vf)
        else
            call s_read_parallel_data_files(q_cons_vf)
        end if

    end subroutine s_read_data_files

    !>  The purpose of this procedure is to first verify that an
        !!      input file has been made available by the user. Provided
        !!      that this is so, the input file is then read in.
    impure subroutine s_read_input_file

        ! Relative path to the input file provided by the user
        character(LEN=name_len), parameter :: file_path = './simulation.inp'

        logical :: file_exist !<
            !! Logical used to check the existence of the input file

        integer :: iostatus
            !! Integer to check iostat of file read

        character(len=1000) :: line

        ! Namelist of the global parameters which may be specified by user
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
            fluid_pp, probe_wrt, prim_vars_wrt, &
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
            Ca, Web, Re_inv, &
            acoustic_source, acoustic, num_source, &
            polytropic, thermal, &
            integral, integral_wrt, num_integrals, &
            polydisperse, poly_sigma, qbmm, &
            relax, relax_model, &
            palpha_eps, ptgalpha_eps, &
            file_per_process, sigma, &
            pi_fac, adv_n, adap_dt, adap_dt_tol, adap_dt_max_iters, &
            bf_x, bf_y, bf_z, &
            k_x, k_y, k_z, w_x, w_y, w_z, p_x, p_y, p_z, &
            g_x, g_y, g_z, n_start, t_save, t_stop, &
            cfl_adap_dt, cfl_const_dt, cfl_target, &
            surface_tension, bubbles_lagrange, lag_params, &
            hyperelasticity, R0ref, num_bc_patches, Bx0, powell, &
            cont_damage, tau_star, cont_damage_s, alpha_bar, &
            alf_factor, num_igr_iters, num_igr_warm_start_iters, &
            int_comp, ic_eps, ic_beta, nv_uvm_out_of_core, &
            nv_uvm_igr_temps_on_gpu, nv_uvm_pref_gpu, down_sample

        ! Checking that an input file has been provided by the user. If it
        ! has, then the input file is read in, otherwise, simulation exits.
        inquire (FILE=trim(file_path), EXIST=file_exist)

        if (file_exist) then
            open (1, FILE=trim(file_path), &
                  FORM='formatted', &
                  ACTION='read', &
                  STATUS='old')
            read (1, NML=user_inputs, iostat=iostatus)

            if (iostatus /= 0) then
                backspace (1)
                read (1, fmt='(A)') line
                print *, 'Invalid line in namelist: '//trim(line)
                call s_mpi_abort('Invalid line in simulation.inp. It is '// &
                                 'likely due to a datatype mismatch. Exiting.')
            end if

            close (1)

            if ((bf_x) .or. (bf_y) .or. (bf_z)) then
                bodyForces = .true.
            endif

            ! Store m,n,p into global m,n,p
            m_glb = m
            n_glb = n
            p_glb = p

            call s_update_cell_bounds(cells_bounds, m, n, p)

            if (cfl_adap_dt .or. cfl_const_dt) cfl_dt = .true.

            if (any((/bc_x%beg, bc_x%end, bc_y%beg, bc_y%end, bc_z%beg, bc_z%end/) == -17) .or. &
                num_bc_patches > 0) then
                bc_io = .true.
            endif

        else
            call s_mpi_abort(trim(file_path)//' is missing. Exiting.')
        end if

    end subroutine s_read_input_file

    !> The goal of this procedure is to verify that each of the
    !!      user provided inputs is valid and that their combination
    !!      constitutes a meaningful configuration for the simulation.
    impure subroutine s_check_input_file

        ! Relative path to the current directory file in the case directory
        character(LEN=path_len) :: file_path

        ! Logical used to check the existence of the current directory file
        logical :: file_exist

        ! Logistics
        file_path = trim(case_dir)//'/.'

        call my_inquire(file_path, file_exist)

        if (file_exist .neqv. .true.) then
            call s_mpi_abort(trim(file_path)//' is missing. Exiting.')
        end if

        call s_check_inputs_common()
        call s_check_inputs()

    end subroutine s_check_input_file

        !!              initial condition and grid data files. The cell-average
        !!              conservative variables constitute the former, while the
        !!              cell-boundary locations in x-, y- and z-directions make
        !!              up the latter. This procedure also calculates the cell-
        !!              width distributions from the cell-boundary locations.
        !! @param q_cons_vf Cell-averaged conservative variables
    impure subroutine s_read_serial_data_files(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf


        character(LEN=path_len + 2*name_len) :: t_step_dir !<
            !! Relative path to the starting time-step directory

        character(LEN=path_len + 3*name_len) :: file_path !<
            !! Relative path to the grid and conservative variables data files

        logical :: file_exist !<
        ! Logical used to check the existence of the data files

        integer :: i, r !< Generic loop iterator

        ! Confirming that the directory from which the initial condition and
        ! the grid data files are to be read in exists and exiting otherwise
        if (cfl_dt) then
            write (t_step_dir, '(A,I0,A,I0)') &
                trim(case_dir)//'/p_all/p', proc_rank, '/', n_start
        else
            write (t_step_dir, '(A,I0,A,I0)') &
                trim(case_dir)//'/p_all/p', proc_rank, '/', t_step_start
        end if

        file_path = trim(t_step_dir)//'/.'
        call my_inquire(file_path, file_exist)

        if (file_exist .neqv. .true.) then
            call s_mpi_abort(trim(file_path)//' is missing. Exiting.')
        end if

        if (bc_io) then
            call s_read_serial_boundary_condition_files(t_step_dir, bc_type)
        else
            call s_assign_default_bc_type(bc_type)
        end if

        ! Cell-boundary Locations in x-direction
        file_path = trim(t_step_dir)//'/x_cb.dat'

        inquire (FILE=trim(file_path), EXIST=file_exist)

        if (file_exist) then
            open (2, FILE=trim(file_path), &
                  FORM='unformatted', &
                  ACTION='read', &
                  STATUS='old')
            read (2) x_cb(-1:m); close (2)
        else
            call s_mpi_abort(trim(file_path)//' is missing. Exiting.')
        end if

        dx(0:m) = x_cb(0:m) - x_cb(-1:m - 1)
        x_cc(0:m) = x_cb(-1:m - 1) + dx(0:m)/2._wp

        if (ib) then
            do i = 1, num_ibs
                if (patch_ib(i)%c > 0) then
                    Np = int((patch_ib(i)%p*patch_ib(i)%c/dx(0))*20) + int(((patch_ib(i)%c - patch_ib(i)%p*patch_ib(i)%c)/dx(0))*20) + 1
                end if
            end do
        end if

        ! Cell-boundary Locations in y-direction
        if (n > 0) then

            file_path = trim(t_step_dir)//'/y_cb.dat'

            inquire (FILE=trim(file_path), EXIST=file_exist)

            if (file_exist) then
                open (2, FILE=trim(file_path), &
                      FORM='unformatted', &
                      ACTION='read', &
                      STATUS='old')
                read (2) y_cb(-1:n); close (2)
            else
                call s_mpi_abort(trim(file_path)//' is missing. Exiting.')
            end if

            dy(0:n) = y_cb(0:n) - y_cb(-1:n - 1)
            y_cc(0:n) = y_cb(-1:n - 1) + dy(0:n)/2._wp

        end if

        ! Cell-boundary Locations in z-direction
        if (p > 0) then

            file_path = trim(t_step_dir)//'/z_cb.dat'

            inquire (FILE=trim(file_path), EXIST=file_exist)

            if (file_exist) then
                open (2, FILE=trim(file_path), &
                      FORM='unformatted', &
                      ACTION='read', &
                      STATUS='old')
                read (2) z_cb(-1:p); close (2)
            else
                call s_mpi_abort(trim(file_path)//' is missing. Exiting.')
            end if

            dz(0:p) = z_cb(0:p) - z_cb(-1:p - 1)
            z_cc(0:p) = z_cb(-1:p - 1) + dz(0:p)/2._wp

        end if

        do i = 1, sys_size
            write (file_path, '(A,I0,A)') &
                trim(t_step_dir)//'/q_cons_vf', i, '.dat'
            inquire (FILE=trim(file_path), EXIST=file_exist)
            if (file_exist) then
                open (2, FILE=trim(file_path), &
                      FORM='unformatted', &
                      ACTION='read', &
                      STATUS='old')
                read (2) q_cons_vf(i)%sf(0:m, 0:n, 0:p); close (2)
            else
                call s_mpi_abort(trim(file_path)//' is missing. Exiting.')
            end if
        end do

        if (bubbles_euler .or. elasticity) then
            ! Read pb and mv for non-polytropic qbmm
            if (qbmm .and. .not. polytropic) then
                do i = 1, nb
                    do r = 1, nnode
                        write (file_path, '(A,I0,A)') &
                            trim(t_step_dir)//'/pb', sys_size + (i - 1)*nnode + r, '.dat'
                        inquire (FILE=trim(file_path), EXIST=file_exist)
                        if (file_exist) then
                            open (2, FILE=trim(file_path), &
                                  FORM='unformatted', &
                                  ACTION='read', &
                                  STATUS='old')
                            read (2) pb_ts(1)%sf(0:m, 0:n, 0:p, r, i); close (2)
                        else
                            call s_mpi_abort(trim(file_path)//' is missing. Exiting.')
                        end if
                    end do
                end do
                do i = 1, nb
                    do r = 1, nnode
                        write (file_path, '(A,I0,A)') &
                            trim(t_step_dir)//'/mv', sys_size + (i - 1)*nnode + r, '.dat'
                        inquire (FILE=trim(file_path), EXIST=file_exist)
                        if (file_exist) then
                            open (2, FILE=trim(file_path), &
                                  FORM='unformatted', &
                                  ACTION='read', &
                                  STATUS='old')
                            read (2) mv_ts(1)%sf(0:m, 0:n, 0:p, r, i); close (2)
                        else
                            call s_mpi_abort(trim(file_path)//' is missing. Exiting.')
                        end if
                    end do
                end do
            end if
        end if

        ! Read IBM Data
        if (ib) then
            ! Read IB markers
            write (file_path, '(A,I0,A)') &
                trim(t_step_dir)//'/ib.dat'
            inquire (FILE=trim(file_path), EXIST=file_exist)
            if (file_exist) then
                open (2, FILE=trim(file_path), &
                        FORM='unformatted', &
                        ACTION='read', &
                        STATUS='old')
                read (2) ib_markers%sf(0:m, 0:n, 0:p); close (2)
            else
                call s_mpi_abort(trim(file_path)//' is missing. Exiting.')
            end if

            ! Read Levelset
            write (file_path, '(A)') &
                trim(t_step_dir)//'/levelset.dat'
            inquire (FILE=trim(file_path), EXIST=file_exist)
            if (file_exist) then
                open (2, FILE=trim(file_path), &
                        FORM='unformatted', &
                        ACTION='read', &
                        STATUS='old')
                read (2) levelset%sf(0:m, 0:n, 0:p, 1:num_ibs); close (2)
                ! print*, 'check', STL_levelset(106, 50, 0, 1)
            else
                call s_mpi_abort(trim(file_path)//' is missing. Exiting.')
            end if

            ! Read Levelset Norm
            write (file_path, '(A)') &
                trim(t_step_dir)//'/levelset_norm.dat'
            inquire (FILE=trim(file_path), EXIST=file_exist)
            if (file_exist) then
                open (2, FILE=trim(file_path), &
                        FORM='unformatted', &
                        ACTION='read', &
                        STATUS='old')
                read (2) levelset_norm%sf(0:m, 0:n, 0:p, 1:num_ibs, 1:3); close (2)
            else
                call s_mpi_abort(trim(file_path)//' is missing. Exiting.')
            end if

            do i = 1, num_ibs
                if (patch_ib(i)%c > 0) then
                    allocate (airfoil_grid_u(1:Np))
                    allocate (airfoil_grid_l(1:Np))

                    write (file_path, '(A)') &
                        trim(t_step_dir)//'/airfoil_u.dat'
                    inquire (FILE=trim(file_path), EXIST=file_exist)
                    if (file_exist) then
                        open (2, FILE=trim(file_path), &
                              FORM='unformatted', &
                              ACTION='read', &
                              STATUS='old')
                        read (2) airfoil_grid_u; close (2)
                    else
                        call s_mpi_abort(trim(file_path)//' is missing. Exiting.')
                    end if

                    write (file_path, '(A)') &
                        trim(t_step_dir)//'/airfoil_l.dat'
                    inquire (FILE=trim(file_path), EXIST=file_exist)
                    if (file_exist) then
                        open (2, FILE=trim(file_path), &
                              FORM='unformatted', &
                              ACTION='read', &
                              STATUS='old')
                        read (2) airfoil_grid_l; close (2)
                    else
                        call s_mpi_abort(trim(file_path)//' is missing. Exiting.')
                    end if
                end if
            end do

        end if

    end subroutine s_read_serial_data_files

        !! @param q_cons_vf Conservative variables
    impure subroutine s_read_parallel_data_files(q_cons_vf)

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: q_cons_vf

#ifdef MFC_MPI

        real(wp), allocatable, dimension(:) :: x_cb_glb, y_cb_glb, z_cb_glb

        integer :: ifile, ierr, data_size
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(KIND=MPI_OFFSET_KIND) :: disp
        integer(KIND=MPI_OFFSET_KIND) :: m_MOK, n_MOK, p_MOK
        integer(KIND=MPI_OFFSET_KIND) :: WP_MOK, var_MOK, str_MOK
        integer(KIND=MPI_OFFSET_KIND) :: NVARS_MOK
        integer(KIND=MPI_OFFSET_KIND) :: MOK

        character(LEN=path_len + 2*name_len) :: file_loc
        logical :: file_exist

        character(len=10) :: t_step_start_string

        integer :: i, j

        ! Downsampled data variables
        integer :: m_ds, n_ds, p_ds
        integer :: m_glb_ds, n_glb_ds, p_glb_ds
        integer :: m_glb_read, n_glb_read, p_glb_read ! data size of read

        allocate (x_cb_glb(-1:m_glb))
        allocate (y_cb_glb(-1:n_glb))
        allocate (z_cb_glb(-1:p_glb))

        ! Read in cell boundary locations in x-direction
        file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//'x_cb.dat'
        inquire (FILE=trim(file_loc), EXIST=file_exist)

        if(down_sample) then
            m_ds = INT((m+1)/3) - 1
            n_ds = INT((n+1)/3) - 1
            p_ds = INT((p+1)/3) - 1

            m_glb_ds = INT((m_glb+1)/3) - 1
            n_glb_ds = INT((n_glb+1)/3) - 1
            p_glb_ds = INT((p_glb+1)/3) - 1
        end if

        if (file_exist) then
            data_size = m_glb + 2
            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
            call MPI_FILE_READ(ifile, x_cb_glb, data_size, mpi_p, status, ierr)
            call MPI_FILE_CLOSE(ifile, ierr)
        else
            call s_mpi_abort('File '//trim(file_loc)//' is missing. Exiting.')
        end if

        ! Assigning local cell boundary locations
        x_cb(-1:m) = x_cb_glb((start_idx(1) - 1):(start_idx(1) + m))
        ! Computing the cell width distribution
        dx(0:m) = x_cb(0:m) - x_cb(-1:m - 1)
        ! Computing the cell center locations
        x_cc(0:m) = x_cb(-1:m - 1) + dx(0:m)/2._wp

        if (ib) then
            do i = 1, num_ibs
                if (patch_ib(i)%c > 0) then
                    Np = int((patch_ib(i)%p*patch_ib(i)%c/dx(0))*20) + int(((patch_ib(i)%c - patch_ib(i)%p*patch_ib(i)%c)/dx(0))*20) + 1
                    allocate (MPI_IO_airfoil_IB_DATA%var(1:2*Np))
                    print *, "HERE Np", Np
                end if
            end do
        end if

        if (n > 0) then
            ! Read in cell boundary locations in y-direction
            file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//'y_cb.dat'
            inquire (FILE=trim(file_loc), EXIST=file_exist)

            if (file_exist) then
                data_size = n_glb + 2
                call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
                call MPI_FILE_READ(ifile, y_cb_glb, data_size, mpi_p, status, ierr)
                call MPI_FILE_CLOSE(ifile, ierr)
            else
                call s_mpi_abort('File '//trim(file_loc)//' is missing. Exiting.')
            end if

            ! Assigning local cell boundary locations
            y_cb(-1:n) = y_cb_glb((start_idx(2) - 1):(start_idx(2) + n))
            ! Computing the cell width distribution
            dy(0:n) = y_cb(0:n) - y_cb(-1:n - 1)
            ! Computing the cell center locations
            y_cc(0:n) = y_cb(-1:n - 1) + dy(0:n)/2._wp

            if (p > 0) then
                ! Read in cell boundary locations in z-direction
                file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//'z_cb.dat'
                inquire (FILE=trim(file_loc), EXIST=file_exist)

                if (file_exist) then
                    data_size = p_glb + 2
                    call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
                    call MPI_FILE_READ(ifile, z_cb_glb, data_size, mpi_p, status, ierr)
                    call MPI_FILE_CLOSE(ifile, ierr)
                else
                    call s_mpi_abort('File '//trim(file_loc)//'is missing. Exiting.')
                end if

                ! Assigning local cell boundary locations
                z_cb(-1:p) = z_cb_glb((start_idx(3) - 1):(start_idx(3) + p))
                ! Computing the cell width distribution
                dz(0:p) = z_cb(0:p) - z_cb(-1:p - 1)
                ! Computing the cell center locations
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
            file_loc = trim(case_dir)//'/restart_data/lustre_'//trim(t_step_start_string)//trim(mpiiofs)//trim(file_loc)
            inquire (FILE=trim(file_loc), EXIST=file_exist)

            if (file_exist) then
                call MPI_FILE_OPEN(MPI_COMM_SELF, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

                ! Initialize MPI data I/O
                if(down_sample) then
                    call s_initialize_mpi_data_ds(q_cons_vf)
                else
                    if (ib) then
                        call s_initialize_mpi_data(q_cons_vf, ib_markers, &
                            levelset, levelset_norm)
                    else
                        call s_initialize_mpi_data(q_cons_vf)
                    end if
                end if

                if(down_sample) then
                    ! Size of local arrays
                    data_size = (m_ds + 3)*(n_ds + 3)*(p_ds + 3)
                    m_glb_read = m_glb_ds + 1
                    n_glb_read = n_glb_ds + 1
                    p_glb_read = p_glb_ds + 1
                else
                    ! Size of local arrays
                    data_size = (m + 1)*(n + 1)*(p + 1)
                    m_glb_read = m_glb + 1
                    n_glb_read = n_glb + 1
                    p_glb_read = p_glb + 1
                end if

                ! Resize some integers so MPI can read even the biggest file
                m_MOK = int(m_glb_read + 1, MPI_OFFSET_KIND)
                n_MOK = int(m_glb_read + 1, MPI_OFFSET_KIND)
                p_MOK = int(m_glb_read + 1, MPI_OFFSET_KIND)
                WP_MOK = int(8._wp, MPI_OFFSET_KIND)
                MOK = int(1._wp, MPI_OFFSET_KIND)
                str_MOK = int(name_len, MPI_OFFSET_KIND)
                NVARS_MOK = int(sys_size, MPI_OFFSET_KIND)

                ! Read the data for each variable
                if (bubbles_euler .or. elasticity) then

                    do i = 1, sys_size!adv_idx%end
                        var_MOK = int(i, MPI_OFFSET_KIND)

                        call MPI_FILE_READ(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                           mpi_p, status, ierr)
                    end do
                    !Read pb and mv for non-polytropic qbmm
                    if (qbmm .and. .not. polytropic) then
                        do i = sys_size + 1, sys_size + 2*nb*nnode
                            var_MOK = int(i, MPI_OFFSET_KIND)

                            call MPI_FILE_READ(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                               mpi_p, status, ierr)
                        end do
                    end if
                else
                    if(down_sample) then
                        do i = 1, sys_size
                            var_MOK = int(i, MPI_OFFSET_KIND)

                            call MPI_FILE_READ(ifile, q_cons_temp(i)%sf, data_size, &
                                               mpi_p, status, ierr)
                        end do
                    else
                        do i = 1, sys_size
                            var_MOK = int(i, MPI_OFFSET_KIND)

                            call MPI_FILE_READ(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                               mpi_p, status, ierr)
                        end do
                    end if
                end if

                call s_mpi_barrier()

                call MPI_FILE_CLOSE(ifile, ierr)

                if (ib) then
                    ! Read IB Markers
                    write (file_loc, '(A)') 'ib.dat'
                    file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
                    inquire (FILE=trim(file_loc), EXIST=file_exist)

                    if (file_exist) then

                        call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

                        disp = 0

                        call MPI_FILE_SET_VIEW(ifile, disp, MPI_INTEGER, MPI_IO_IB_DATA%view, &
                                               'native', mpi_info_int, ierr)
                        call MPI_FILE_READ(ifile, MPI_IO_IB_DATA%var%sf, data_size, &
                                           MPI_INTEGER, status, ierr)

                    else
                        call s_mpi_abort('File '//trim(file_loc)//' is missing. Exiting.')
                    end if

                    ! Read Levelset
                    write (file_loc, '(A)') 'levelset.dat'
                    file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
                    inquire (FILE=trim(file_loc), EXIST=file_exist)

                    if (file_exist) then

                        call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

                        disp = 0

                        call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, MPI_IO_levelset_DATA%view, &
                                               'native', mpi_info_int, ierr)
                        call MPI_FILE_READ(ifile, MPI_IO_levelset_DATA%var%sf, data_size * num_ibs, &
                                           mpi_p, status, ierr)

                    else
                        call s_mpi_abort('File '//trim(file_loc)//' is missing. Exiting.')
                    end if

                    ! Read Levelset Norm
                    write (file_loc, '(A)') 'levelset_norm.dat'
                    file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
                    inquire (FILE=trim(file_loc), EXIST=file_exist)

                    if (file_exist) then

                        call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

                        disp = 0

                        call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, MPI_IO_levelsetnorm_DATA%view, &
                                               'native', mpi_info_int, ierr)
                        call MPI_FILE_READ(ifile, MPI_IO_levelsetnorm_DATA%var%sf, data_size * num_ibs * 3, &
                                           mpi_p, status, ierr)

                    else
                        call s_mpi_abort('File '//trim(file_loc)//' is missing. Exiting.')
                    end if

                end if

            else
                call s_mpi_abort('File '//trim(file_loc)//' is missing. Exiting.')
            end if
        else

            ! Open the file to read conservative variables
            if (cfl_dt) then
                write (file_loc, '(I0,A)') n_start, '.dat'
            else
                write (file_loc, '(I0,A)') t_step_start, '.dat'
            end if
            file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
            inquire (FILE=trim(file_loc), EXIST=file_exist)

            if (file_exist) then
                call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

                ! Initialize MPI data I/O

                if (ib) then
                    call s_initialize_mpi_data(q_cons_vf, ib_markers, &
                        levelset, levelset_norm)
                else

                    call s_initialize_mpi_data(q_cons_vf)

                end if


                ! Size of local arrays
                data_size = (m + 1)*(n + 1)*(p + 1)

                ! Resize some integers so MPI can read even the biggest file
                m_MOK = int(m_glb + 1, MPI_OFFSET_KIND)
                n_MOK = int(n_glb + 1, MPI_OFFSET_KIND)
                p_MOK = int(p_glb + 1, MPI_OFFSET_KIND)
                WP_MOK = int(8._wp, MPI_OFFSET_KIND)
                MOK = int(1._wp, MPI_OFFSET_KIND)
                str_MOK = int(name_len, MPI_OFFSET_KIND)
                NVARS_MOK = int(sys_size, MPI_OFFSET_KIND)

                ! Read the data for each variable
                if (bubbles_euler .or. elasticity) then
                    do i = 1, sys_size !adv_idx%end
                        var_MOK = int(i, MPI_OFFSET_KIND)
                        ! Initial displacement to skip at beginning of file
                        disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                        call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, MPI_IO_DATA%view(i), &
                                               'native', mpi_info_int, ierr)
                        call MPI_FILE_READ(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                           mpi_p, status, ierr)
                    end do
                    !Read pb and mv for non-polytropic qbmm
                    if (qbmm .and. .not. polytropic) then
                        do i = sys_size + 1, sys_size + 2*nb*nnode
                            var_MOK = int(i, MPI_OFFSET_KIND)
                            ! Initial displacement to skip at beginning of file
                            disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                            call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, MPI_IO_DATA%view(i), &
                                                   'native', mpi_info_int, ierr)
                            call MPI_FILE_READ(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                               mpi_p, status, ierr)
                        end do
                    end if
                else
                    do i = 1, sys_size
                        var_MOK = int(i, MPI_OFFSET_KIND)

                        ! Initial displacement to skip at beginning of file
                        disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                        call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, MPI_IO_DATA%view(i), &
                                               'native', mpi_info_int, ierr)
                        call MPI_FILE_READ(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                           mpi_p, status, ierr)

                    end do
                end if

                call s_mpi_barrier()

                call MPI_FILE_CLOSE(ifile, ierr)

                if (ib) then

                    ! Read IB Markers
                    write (file_loc, '(A)') 'ib.dat'
                    file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
                    inquire (FILE=trim(file_loc), EXIST=file_exist)

                    if (file_exist) then

                        call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

                        disp = 0

                        call MPI_FILE_SET_VIEW(ifile, disp, MPI_INTEGER, MPI_IO_IB_DATA%view, &
                                               'native', mpi_info_int, ierr)
                        call MPI_FILE_READ(ifile, MPI_IO_IB_DATA%var%sf, data_size, &
                                           MPI_INTEGER, status, ierr)

                    else
                        call s_mpi_abort('File '//trim(file_loc)//' is missing. Exiting.')
                    end if

                    ! Read Levelset
                    write (file_loc, '(A)') 'levelset.dat'
                    file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
                    inquire (FILE=trim(file_loc), EXIST=file_exist)

                    if (file_exist) then

                        call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

                        disp = 0

                        call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, MPI_IO_levelset_DATA%view, &
                                               'native', mpi_info_int, ierr)
                        call MPI_FILE_READ(ifile, MPI_IO_levelset_DATA%var%sf, data_size * num_ibs, &
                                           mpi_p, status, ierr)

                    else
                        call s_mpi_abort('File '//trim(file_loc)//' is missing. Exiting.')
                    end if

                    ! Read Levelset Norm
                    write (file_loc, '(A)') 'levelset_norm.dat'
                    file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
                    inquire (FILE=trim(file_loc), EXIST=file_exist)

                    if (file_exist) then

                        call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

                        disp = 0

                        call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, MPI_IO_levelsetnorm_DATA%view, &
                                               'native', mpi_info_int, ierr)
                        call MPI_FILE_READ(ifile, MPI_IO_levelsetnorm_DATA%var%sf, data_size * num_ibs * 3, &
                                           mpi_p, status, ierr)

                    else
                        call s_mpi_abort('File '//trim(file_loc)//' is missing. Exiting.')
                    end if

                end if

            else
                call s_mpi_abort('File '//trim(file_loc)//' is missing. Exiting.')
            end if

        end if

        if (ib) then

            do j = 1, num_ibs
                if (patch_ib(j)%c > 0) then

                    print *, "HERE Np", Np

                    allocate (airfoil_grid_u(1:Np))
                    allocate (airfoil_grid_l(1:Np))

                    write (file_loc, '(A)') 'airfoil_l.dat'
                    file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
                    inquire (FILE=trim(file_loc), EXIST=file_exist)
                    if (file_exist) then

                        call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

                        ! Initial displacement to skip at beginning of file
                        disp = 0

                        call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, MPI_IO_airfoil_IB_DATA%view(1), &
                                               'native', mpi_info_int, ierr)
                        call MPI_FILE_READ(ifile, MPI_IO_airfoil_IB_DATA%var(1:Np), 3*Np, &
                                           mpi_p, status, ierr)

                    end if

                    write (file_loc, '(A)') 'airfoil_u.dat'
                    file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
                    inquire (FILE=trim(file_loc), EXIST=file_exist)
                    if (file_exist) then

                        call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

                        ! Initial displacement to skip at beginning of file
                        disp = 0

                        call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, MPI_IO_airfoil_IB_DATA%view(2), &
                                               'native', mpi_info_int, ierr)
                        call MPI_FILE_READ(ifile, MPI_IO_airfoil_IB_DATA%var(Np + 1:2*Np), 3*Np, &
                                           mpi_p, status, ierr)
                    end if

                    do i = 1, Np
                        airfoil_grid_l(i)%x = MPI_IO_airfoil_IB_DATA%var(i)%x
                        airfoil_grid_l(i)%y = MPI_IO_airfoil_IB_DATA%var(i)%y
                    end do

                    do i = 1, Np
                        airfoil_grid_u(i)%x = MPI_IO_airfoil_IB_DATA%var(Np + i)%x
                        airfoil_grid_u(i)%y = MPI_IO_airfoil_IB_DATA%var(Np + i)%y
                    end do

                end if
            end do
        end if

        deallocate (x_cb_glb, y_cb_glb, z_cb_glb)

        if (bc_io) then
            call s_read_parallel_boundary_condition_files(bc_type)
        else
            call s_assign_default_bc_type(bc_type)
        end if

#endif

    end subroutine s_read_parallel_data_files

        !> The purpose of this procedure is to initialize the
        !!      values of the internal-energy equations of each phase
        !!      from the mass of each phase, the mixture momentum and
        !!      mixture-total-energy equations.
        !! @param v_vf conservative variables
    subroutine s_initialize_internal_energy_equations(v_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: v_vf

        real(wp) :: rho
        real(wp) :: dyn_pres
        real(wp) :: gamma
        real(wp) :: pi_inf
        real(wp) :: qv
        real(wp), dimension(2) :: Re
        real(wp) :: pres, T

        integer :: i, j, k, l, c

        real(wp), dimension(num_species) :: rhoYks

        real(wp) :: pres_mag

        pres_mag = 0._wp

        T = dflt_T_guess

        do j = 0, m
            do k = 0, n
                do l = 0, p

                    call s_convert_to_mixture_variables(v_vf, j, k, l, rho, gamma, pi_inf, qv, Re)

                    dyn_pres = 0._wp
                    do i = mom_idx%beg, mom_idx%end
                        dyn_pres = dyn_pres + 5.e-1_wp*v_vf(i)%sf(j, k, l)*v_vf(i)%sf(j, k, l) &
                                   /max(rho, sgm_eps)
                    end do

                    if (chemistry) then
                        do c = 1, num_species
                            rhoYks(c) = v_vf(chemxb + c - 1)%sf(j, k, l)
                        end do
                    end if

                    if (mhd) then
                        if (n == 0) then
                            pres_mag = 0.5_wp*(Bx0**2 + v_vf(B_idx%beg)%sf(j, k, l)**2 + v_vf(B_idx%beg+1)%sf(j, k, l)**2)
                        else
                            pres_mag = 0.5_wp*(v_vf(B_idx%beg)%sf(j, k, l)**2 + v_vf(B_idx%beg+1)%sf(j, k, l)**2 + v_vf(B_idx%beg+2)%sf(j, k, l)**2)
                        end if
                    end if

                    call s_compute_pressure(v_vf(E_idx)%sf(j, k, l), 0._wp, &
                                            dyn_pres, pi_inf, gamma, rho, qv, rhoYks, pres, T, pres_mag = pres_mag)

                    do i = 1, num_fluids
                        v_vf(i + internalEnergies_idx%beg - 1)%sf(j, k, l) = v_vf(i + adv_idx%beg - 1)%sf(j, k, l)* &
                                                                             (fluid_pp(i)%gamma*pres + fluid_pp(i)%pi_inf) &
                                                                             + v_vf(i + cont_idx%beg - 1)%sf(j, k, l)*fluid_pp(i)%qv
                    end do

                end do
            end do
        end do

    end subroutine s_initialize_internal_energy_equations

    impure subroutine s_perform_time_step(t_step, time_avg)
        integer, intent(inout) :: t_step
        real(wp), intent(inout) :: time_avg


        integer :: i

        if (cfl_dt) then
            if (cfl_const_dt .and. t_step == 0) call s_compute_dt()

            if (cfl_adap_dt) call s_compute_dt()

            if (t_step == 0) dt_init = dt

            if (dt < 1.e-3_wp*dt_init .and. cfl_adap_dt .and. proc_rank == 0) then
                print*, "Delta t = ", dt
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
                print '(" [", I3, "%] Time ", ES16.6, " dt = ", ES16.6, " @ Time Step = ", I8,  " Time Avg = ", ES16.6,  " Time/step = ", ES12.6, "")', &
                    int(ceiling(100._wp*(mytime/t_stop))), &
                    mytime, &
                    dt, &
                    t_step, &
                    wall_time_avg, &
                    wall_time
            end if
        else
            if (proc_rank == 0 .and. mod(t_step - t_step_start, t_step_print) == 0) then
                print '(" [", I3, "%]  Time step ", I8, " of ", I0, " @ t_step = ", I8,  " Time Avg = ", ES12.6,  " Time/step= ", ES12.6, "")', &
                   int(ceiling(100._wp*(real(t_step - t_step_start)/(t_step_stop - t_step_start + 1)))), &
                    t_step - t_step_start + 1, &
                    t_step_stop - t_step_start + 1, &
                    t_step, &
                    wall_time_avg, &
                    wall_time
            end if
        end if

        if (probe_wrt) then
            do i = 1, sys_size
                $:GPU_UPDATE(host='[q_cons_ts(1)%vf(i)%sf]')
            end do
        end if

        call s_compute_derived_variables(t_step)


#ifdef DEBUG
        print *, 'Computed derived vars'
#endif

        mytime = mytime + dt

        ! Total-variation-diminishing (TVD) Runge-Kutta (RK) time-steppers
        if (time_stepper == 1) then
            call s_1st_order_tvd_rk(t_step, time_avg)
        elseif (time_stepper == 2) then
            call s_2nd_order_tvd_rk(t_step, time_avg)
        elseif (time_stepper == 3 .and. (.not. adap_dt)) then
            call s_3rd_order_tvd_rk(t_step, time_avg)
        elseif (time_stepper == 3 .and. adap_dt) then
            call s_strang_splitting(t_step, time_avg)
        end if

        if (relax) call s_infinite_relaxation_k(q_cons_ts(1)%vf)

        ! Time-stepping loop controls

        t_step = t_step + 1

    end subroutine s_perform_time_step

    impure subroutine s_save_performance_metrics(time_avg, time_final, io_time_avg, io_time_final, proc_time, io_proc_time, file_exists)

        real(wp), intent(inout) :: time_avg, time_final
        real(wp), intent(inout) :: io_time_avg, io_time_final
        real(wp), dimension(:), intent(inout) :: proc_time
        real(wp), dimension(:), intent(inout) :: io_proc_time
        logical, intent(inout) :: file_exists

        real(wp) :: grind_time

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

            grind_time = time_final * 1.0e9_wp / &
                (real(sys_size, wp) * real(maxval((/1, m_glb/)), wp) * &
                 real(maxval((/1, n_glb/)), wp) * real(maxval((/1, p_glb/)), wp))

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

    impure subroutine s_save_data(t_step, start, finish, io_time_avg, nt)
        integer, intent(inout) :: t_step
        real(wp), intent(inout) :: start, finish, io_time_avg
        integer, intent(inout) :: nt

        integer :: i, j, k, l

        integer :: save_count

        call cpu_time(start)
        call nvtxStartRange("SAVE-DATA")
        do i = 1, sys_size
            $:GPU_UPDATE(host='[q_cons_ts(1)%vf(i)%sf]')
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        if (ieee_is_nan(q_cons_ts(1)%vf(i)%sf(j, k, l))) then
                            print *, "NaN(s) in timestep output.", j, k, l, i, proc_rank, t_step, m, n, p
                            error stop "NaN(s) in timestep output."
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
            $:GPU_UPDATE(host='[intfc_rad]')
            do i = 1, nBubs
                if (ieee_is_nan(intfc_rad(i, 1)) .or. intfc_rad(i, 1) <= 0._wp) then
                    call s_mpi_abort("Bubble radius is negative or NaN, please reduce dt.")
                end if
            end do

            $:GPU_UPDATE(host='[q_beta%vf(1)%sf]')
            call s_write_data_files(q_cons_ts(1)%vf, q_T_sf, q_prim_vf, save_count, bc_type, q_beta%vf(1))
            $:GPU_UPDATE(host='[Rmax_stats,Rmin_stats,gas_p,gas_mv,intfc_vel]')
            call s_write_restart_lag_bubbles(save_count) !parallel
            if (lag_params%write_bubbles_stats) call s_write_lag_bubble_stats()
        else
            call s_write_data_files(q_cons_ts(1)%vf, q_T_sf, q_prim_vf, save_count, bc_type)
        end if

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

    impure subroutine s_initialize_modules

        integer :: m_ds, n_ds, p_ds
        integer :: i, j, k, l, x_id, y_id, z_id, ix, iy, iz
        real(wp) :: temp1, temp2, temp3, temp4

        call s_initialize_global_parameters_module()
        !Quadrature weights and nodes for polydisperse simulations
        if (bubbles_euler .and. nb > 1) then
            call s_simpson(weight, R0)
        end if
        !Initialize variables for non-polytropic (Preston) model
        if (bubbles_euler .and. .not. polytropic) then
            call s_initialize_nonpoly()
        end if
        !Initialize pb based on surface tension for qbmm (polytropic)
        if (qbmm .and. polytropic .and. (.not. f_is_default(Web))) then
            pb0(:) = pref + 2._wp*fluid_pp(1)%ss/(R0(:)*R0ref)
            pb0 = pb0/pref
            pref = 1._wp
        end if

        call s_initialize_mpi_common_module()
        call s_initialize_mpi_proxy_module()
        call s_initialize_variables_conversion_module()
        if (grid_geometry == 3) call s_initialize_fftw_module()

        if(bubbles_euler) call s_initialize_bubbles_EE_module()
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

        if(down_sample) then
            m_ds = INT((m+1)/3) - 1
            n_ds = INT((n+1)/3) - 1
            p_ds = INT((p+1)/3) - 1

            allocate(q_cons_temp(1:sys_size))
            do i = 1, sys_size
                allocate(q_cons_temp(i)%sf(-1:m_ds+1,-1:n_ds+1,-1:p_ds+1))
            end do
        end if

        ! Reading in the user provided initial condition and grid data
        if(down_sample) then
            call s_read_data_files(q_cons_temp)
            call s_upsample_data(q_cons_ts(1)%vf, q_cons_temp)
            do i = 1, sys_size
                $:GPU_UPDATE(device='[q_cons_ts(1)%vf(i)%sf]')
            end do
        else
            call s_read_data_files(q_cons_ts(1)%vf)
        end if

        ! Populating the buffers of the grid variables using the boundary conditions
        call s_populate_grid_variables_buffers()

        if (model_eqns == 3) call s_initialize_internal_energy_equations(q_cons_ts(1)%vf)
        if (ib) call s_ibm_setup()
        if (bodyForces) call s_initialize_body_forces_module()
        if (acoustic_source) call s_precalculate_acoustic_spatial_sources()

        ! Initialize the Temperature cache.
        if (chemistry) call s_compute_q_T_sf(q_T_sf, q_cons_ts(1)%vf, idwint)

        ! Computation of parameters, allocation of memory, association of pointers,
        ! and/or execution of any other tasks that are needed to properly configure
        ! the modules. The preparations below DO DEPEND on the grid being complete.
        if (igr) then
            call s_initialize_igr_module()
        else
            if (recon_type == WENO_TYPE) then
                call s_initialize_weno_module()
            elseif (recon_type == MUSCL_TYPE) then
                call s_initialize_muscl_module()
            end if
            call s_initialize_cbc_module()
            call s_initialize_riemann_solvers_module()
        end if

        call s_initialize_derived_variables()
        if (bubbles_lagrange) call s_initialize_bubbles_EL_module(q_cons_ts(1)%vf)

        if (hypoelasticity) call s_initialize_hypoelastic_module()
        if (hyperelasticity) call s_initialize_hyperelastic_module()

        if (mhd .and. powell) call s_initialize_mhd_powell_module

    end subroutine s_initialize_modules

    impure subroutine s_initialize_mpi_domain
        integer :: ierr
#ifdef MFC_OpenACC
        real(wp) :: starttime, endtime
        integer :: num_devices, local_size, num_nodes, ppn, my_device_num
        integer :: dev, devNum, local_rank
#ifdef MFC_MPI
        integer :: local_comm
#endif
        integer(acc_device_kind) :: devtype
#endif

        ! Initializing MPI execution environment

        call s_mpi_initialize()

        ! Bind GPUs if OpenACC is enabled
#ifdef MFC_OpenACC
#ifndef MFC_MPI
        local_size = 1
        local_rank = 0
#else
        call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, &
                                 MPI_INFO_NULL, local_comm, ierr)
        call MPI_Comm_size(local_comm, local_size, ierr)
        call MPI_Comm_rank(local_comm, local_rank, ierr)
#endif

        devtype = acc_get_device_type()
        devNum = acc_get_num_devices(devtype)
        dev = mod(local_rank, devNum)

        call acc_set_device_num(dev, devtype)
#endif

        ! The rank 0 processor assigns default values to the user inputs prior to
        ! reading them in from the input file. Next, the user inputs are read and
        ! their consistency is checked. The identification of any inconsistencies
        ! will result in the termination of the simulation.
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
#ifdef MFC_OpenACC
!&<
                "with OpenACC offloading"
!&>
#else
                "on CPUs"
#endif
        end if

        ! Broadcasting the user inputs to all of the processors and performing the
        ! parallel computational domain decomposition. Neither procedure has to be
        ! carried out if the simulation is in fact not truly executed in parallel.

        call s_mpi_bcast_user_inputs()

        call s_initialize_parallel_io()

        call s_mpi_decompose_computational_domain()

    end subroutine s_initialize_mpi_domain

    subroutine s_initialize_gpu_vars
        integer :: i
        !Update GPU DATA
        if (.not. down_sample) then
            do i = 1, sys_size
                $:GPU_UPDATE(device='[q_cons_ts(1)%vf(i)%sf]')
            end do
        end if

        if (qbmm .and. .not. polytropic) then
            $:GPU_UPDATE(device='[pb_ts(1)%sf,mv_ts(1)%sf]')
        end if
        if (chemistry) then
            $:GPU_UPDATE(device='[q_T_sf%sf]')
        end if
        $:GPU_UPDATE(device='[nb,R0ref,Ca,Web,Re_inv,weight,R0, &
            & bubbles_euler,polytropic,polydisperse,qbmm, &
            & ptil,bubble_model,thermal,poly_sigma,adv_n,adap_dt, &
            & adap_dt_tol,adap_dt_max_iters,n_idx,pi_fac,low_Mach]')
        $:GPU_UPDATE(device='[R_n,R_v,phi_vn,phi_nv,Pe_c,Tw,pv,M_n, &
            & M_v,k_n,k_v,pb0,mass_n0,mass_v0,Pe_T,Re_trans_T, &
            & Re_trans_c,Im_trans_T,Im_trans_c,omegaN,mul0,ss, &
            & gamma_v,mu_v,gamma_m,gamma_n,mu_n,gam]')

        $:GPU_UPDATE(device='[acoustic_source, num_source]')
        $:GPU_UPDATE(device='[sigma, surface_tension]')

        $:GPU_UPDATE(device='[dx,dy,dz,x_cb,x_cc,y_cb,y_cc,z_cb,z_cc]')

        $:GPU_UPDATE(device='[bc_x%vb1,bc_x%vb2,bc_x%vb3,bc_x%ve1,bc_x%ve2,bc_x%ve3]')
        $:GPU_UPDATE(device='[bc_y%vb1,bc_y%vb2,bc_y%vb3,bc_y%ve1,bc_y%ve2,bc_y%ve3]')
        $:GPU_UPDATE(device='[bc_z%vb1,bc_z%vb2,bc_z%vb3,bc_z%ve1,bc_z%ve2,bc_z%ve3]')

        $:GPU_UPDATE(device='[bc_x%grcbc_in,bc_x%grcbc_out,bc_x%grcbc_vel_out]')
        $:GPU_UPDATE(device='[bc_y%grcbc_in,bc_y%grcbc_out,bc_y%grcbc_vel_out]')
        $:GPU_UPDATE(device='[bc_z%grcbc_in,bc_z%grcbc_out,bc_z%grcbc_vel_out]')

        $:GPU_UPDATE(device='[relax, relax_model]')
        if (relax) then
            $:GPU_UPDATE(device='[palpha_eps, ptgalpha_eps]')
        end if

        if (ib) then
            $:GPU_UPDATE(device='[ib_markers%sf]')
        end if

        $:GPU_UPDATE(device='[igr, igr_order]')

    end subroutine s_initialize_gpu_vars

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
            elseif (recon_type == MUSCL_TYPE) then
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

        if (surface_tension)  call s_finalize_surface_tension_module()
        if (bodyForces) call s_finalize_body_forces_module()
        if (mhd .and. powell) call s_finalize_mhd_powell_module

        ! Terminating MPI execution environment
        call s_mpi_finalize()
    end subroutine s_finalize_modules

end module m_start_up
