!>
!! @file
!! @brief Contains module m_start_up

#:include 'macros.fpp'

!> @brief Reads and validates user inputs, loads existing grid/IC data, and initializes pre-process modules
module m_start_up

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_mpi_common
    use m_variables_conversion
    use m_grid
    use m_initial_condition
    use m_data_output
    use m_compile_specific
    use m_icpp_patches
    use m_assign_variables
    use m_phase_change
    use m_helper_basic
    use m_helper

#ifdef MFC_MPI
    use mpi
#endif

    use m_check_patches
    use m_check_ib_patches
    use m_helper
    use m_checker_common
    use m_checker
    use m_boundary_common
    use m_boundary_conditions

    implicit none

    private
    public :: s_read_input_file, s_check_input_file, s_read_grid_data_files, s_read_ic_data_files, s_read_serial_grid_data_files, &
        & s_read_serial_ic_data_files, s_read_parallel_grid_data_files, s_read_parallel_ic_data_files, s_check_grid_data_files, &
        & s_initialize_modules, s_initialize_mpi_domain, s_finalize_modules, s_apply_initial_condition, s_save_data, s_read_grid

    abstract interface

        !> Abstract interface for reading grid data files in serial or parallel.
        impure subroutine s_read_abstract_grid_data_files

        end subroutine s_read_abstract_grid_data_files

        !> Abstract interface for reading initial condition data files in serial or parallel.
        impure subroutine s_read_abstract_ic_data_files(q_cons_vf_in)

            import :: scalar_field, integer_field, sys_size, pres_field

            type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf_in

        end subroutine s_read_abstract_ic_data_files
    end interface

    character(LEN=path_len + name_len) :: proc_rank_dir  !< Location of the folder associated with the rank of the local processor
    character(LEN=path_len + 2*name_len), private :: t_step_dir  !< Path to preexisting time-step folder for restart
    procedure(s_read_abstract_grid_data_files), pointer :: s_read_grid_data_files => null()
    procedure(s_read_abstract_ic_data_files), pointer :: s_read_ic_data_files => null()

contains

    !> Reads the configuration file pre_process.inp, in order to populate the parameters in module m_global_parameters.f90 with the
    !! user provided inputs
    impure subroutine s_read_input_file

        character(LEN=name_len) :: file_loc
        logical                 :: file_check
        integer                 :: iostatus
        character(len=1000)     :: line

        namelist /user_inputs/ case_dir, old_grid, old_ic, t_step_old, t_step_start, m, n, p, x_domain, y_domain, z_domain, &
            & stretch_x, stretch_y, stretch_z, a_x, a_y, a_z, x_a, y_a, z_a, x_b, y_b, z_b, model_eqns, num_fluids, mpp_lim, &
            & weno_order, bc_x, bc_y, bc_z, num_patches, hypoelasticity, mhd, patch_icpp, fluid_pp, bub_pp, precision, &
            & parallel_io, mixlayer_vel_profile, mixlayer_vel_coef, mixlayer_perturb, mixlayer_perturb_nk, mixlayer_perturb_k0, &
            & pi_fac, perturb_flow, perturb_flow_fluid, perturb_flow_mag, perturb_sph, perturb_sph_fluid, fluid_rho, cyl_coord, &
            & loops_x, loops_y, loops_z, rhoref, pref, bubbles_euler, R0ref, nb, polytropic, thermal, Ca, Web, Re_inv, &
            & polydisperse, poly_sigma, qbmm, sigR, sigV, dist_type, rhoRV, file_per_process, relax, relax_model, palpha_eps, &
            & ptgalpha_eps, ib, num_ibs, patch_ib, sigma, adv_n, cfl_adap_dt, cfl_const_dt, n_start, n_start_old, &
            & surface_tension, hyperelasticity, pre_stress, elliptic_smoothing, elliptic_smoothing_iters, viscous, &
            & bubbles_lagrange, num_bc_patches, patch_bc, Bx0, relativity, cont_damage, igr, igr_order, down_sample, recon_type, &
            & muscl_order, hyper_cleaning, simplex_perturb, simplex_params, fft_wrt

        file_loc = 'pre_process.inp'
        inquire (FILE=trim(file_loc), EXIST=file_check)

        if (file_check) then
            open (1, FILE=trim(file_loc), form='formatted', STATUS='old', ACTION='read')
            read (1, NML=user_inputs, iostat=iostatus)
            if (iostatus /= 0) then
                backspace (1)
                read (1, fmt='(A)') line
                print *, 'Invalid line in namelist: ' // trim(line)
                call s_mpi_abort('Invalid line in pre_process.inp. It is ' // 'likely due to a datatype mismatch. Exiting.')
            end if
            close (1)

            call s_update_cell_bounds(cells_bounds, m, n, p)

            m_glb = m
            n_glb = n
            p_glb = p

            nGlobal = int(m_glb + 1, kind=8)*int(n_glb + 1, kind=8)*int(p_glb + 1, kind=8)

            if (cfl_adap_dt .or. cfl_const_dt) cfl_dt = .true.

            if (any((/bc_x%beg, bc_x%end, bc_y%beg, bc_y%end, bc_z%beg, bc_z%end/) == BC_DIRICHLET) .or. num_bc_patches > 0) then
                bc_io = .true.
            end if
        else
            call s_mpi_abort('File pre_process.inp is missing. Exiting.')
        end if

    end subroutine s_read_input_file

    !> Checking that the user inputs make sense, i.e. that the individual choices are compatible with the code's options and that
    !! the combination of these choices results into a valid configuration for the pre-process
    impure subroutine s_check_input_file

        character(LEN=len_trim(case_dir)) :: file_loc
        logical                           :: dir_check

        case_dir = adjustl(case_dir)

        file_loc = trim(case_dir) // '/.'

        call my_inquire(file_loc, dir_check)

        if (dir_check .neqv. .true.) then
            print '(A)', 'WARNING: Ensure that compiler flags/choices in Makefiles match your compiler! '
            print '(A)', 'WARNING: Ensure that preprocessor flags are enabled! '
            call s_mpi_abort('Unsupported choice for the value of case_dir.' // 'Exiting.')
        end if

        call s_check_inputs_common()
        call s_check_inputs()

        call s_check_patches()

        if (ib) call s_check_ib_patches()

    end subroutine s_check_input_file

    !> The goal of this subroutine is to read in any preexisting grid data as well as based on the imported grid, complete the
    !! necessary global computational domain parameters.
    impure subroutine s_read_serial_grid_data_files

        character(LEN=len_trim(case_dir) + 3*name_len) :: file_loc
        logical                                        :: dir_check
        logical                                        :: file_check

        write (proc_rank_dir, '(A,I0)') '/p_all/p', proc_rank
        proc_rank_dir = trim(case_dir) // trim(proc_rank_dir)

        write (t_step_dir, '(A,I0)') '/', t_step_start
        t_step_dir = trim(proc_rank_dir) // trim(t_step_dir)

        file_loc = trim(t_step_dir) // '/.'
        call my_inquire(file_loc, dir_check)

        if (dir_check .neqv. .true.) then
            call s_mpi_abort('Time-step folder ' // trim(t_step_dir) // ' is missing. Exiting.')
        end if

        file_loc = trim(t_step_dir) // '/x_cb.dat'
        inquire (FILE=trim(file_loc), EXIST=file_check)

        if (file_check) then
            open (1, FILE=trim(file_loc), form='unformatted', STATUS='old', ACTION='read')
            read (1) x_cb(-1:m)
            close (1)
        else
            call s_mpi_abort('File x_cb.dat is missing in ' // trim(t_step_dir) // '. Exiting.')
        end if

        x_cc(0:m) = (x_cb(0:m) + x_cb(-1:(m - 1)))/2._wp

        dx = minval(x_cb(0:m) - x_cb(-1:m - 1))
        if (num_procs > 1) call s_mpi_reduce_min(dx)

        x_domain%beg = x_cb(-1)
        x_domain%end = x_cb(m)

        if (n > 0) then
            file_loc = trim(t_step_dir) // '/y_cb.dat'
            inquire (FILE=trim(file_loc), EXIST=file_check)

            if (file_check) then
                open (1, FILE=trim(file_loc), form='unformatted', STATUS='old', ACTION='read')
                read (1) y_cb(-1:n)
                close (1)
            else
                call s_mpi_abort('File y_cb.dat is missing in ' // trim(t_step_dir) // '. Exiting.')
            end if

            y_cc(0:n) = (y_cb(0:n) + y_cb(-1:(n - 1)))/2._wp

            dy = minval(y_cb(0:n) - y_cb(-1:n - 1))
            if (num_procs > 1) call s_mpi_reduce_min(dy)

            y_domain%beg = y_cb(-1)
            y_domain%end = y_cb(n)

            if (p > 0) then
                file_loc = trim(t_step_dir) // '/z_cb.dat'
                inquire (FILE=trim(file_loc), EXIST=file_check)

                if (file_check) then
                    open (1, FILE=trim(file_loc), form='unformatted', STATUS='old', ACTION='read')
                    read (1) z_cb(-1:p)
                    close (1)
                else
                    call s_mpi_abort('File z_cb.dat is missing in ' // trim(t_step_dir) // '. Exiting.')
                end if

                z_cc(0:p) = (z_cb(0:p) + z_cb(-1:(p - 1)))/2._wp

                dz = minval(z_cb(0:p) - z_cb(-1:p - 1))
                if (num_procs > 1) call s_mpi_reduce_min(dz)

                z_domain%beg = z_cb(-1)
                z_domain%end = z_cb(p)
            end if
        end if

        ! Clean processor dir and create time-step dir (unless reading preexisting IC)
        if (old_ic .neqv. .true.) then
            call s_delete_directory(trim(proc_rank_dir))
            call s_create_directory(trim(proc_rank_dir) // '/0')
        end if

    end subroutine s_read_serial_grid_data_files

    !> Cell-boundary data are checked for consistency by looking at the (non-)uniform cell-width distributions for all the active
    !! coordinate directions and making sure that all of the cell-widths are positively valued
    impure subroutine s_check_grid_data_files

        if (any(x_cb(0:m) - x_cb(-1:m - 1) <= 0._wp)) then
            call s_mpi_abort('x_cb.dat in ' // trim(t_step_dir) // ' contains non-positive cell-spacings. Exiting.')
        end if

        if (n > 0) then
            if (any(y_cb(0:n) - y_cb(-1:n - 1) <= 0._wp)) then
                call s_mpi_abort('y_cb.dat in ' // trim(t_step_dir) // ' contains non-positive cell-spacings. ' // 'Exiting.')
            end if

            if (p > 0) then
                if (any(z_cb(0:p) - z_cb(-1:p - 1) <= 0._wp)) then
                    call s_mpi_abort('z_cb.dat in ' // trim(t_step_dir) // ' contains non-positive cell-spacings' // ' .Exiting.')
                end if
            end if
        end if

    end subroutine s_check_grid_data_files

    !> The goal of this subroutine is to read in any preexisting initial condition data files so that they may be used by the
    !! pre-process as a starting point in the creation of an all new initial condition.
    impure subroutine s_read_serial_ic_data_files(q_cons_vf_in)

        type(scalar_field), dimension(sys_size), intent(inout)   :: q_cons_vf_in
        character(LEN=len_trim(case_dir) + 3*name_len)           :: file_loc
        character(LEN=int(floor(log10(real(sys_size, wp)))) + 1) :: file_num
        logical                                                  :: file_check
        integer                                                  :: i, r

        do i = 1, sys_size
            write (file_num, '(I0)') i
            file_loc = trim(t_step_dir) // '/q_cons_vf' // trim(file_num) // '.dat'
            inquire (FILE=trim(file_loc), EXIST=file_check)

            if (file_check) then
                open (1, FILE=trim(file_loc), form='unformatted', STATUS='old', ACTION='read')
                read (1) q_cons_vf_in(i)%sf
                close (1)
            else
                call s_mpi_abort('File q_cons_vf' // trim(file_num) // '.dat is missing in ' // trim(t_step_dir) // '. Exiting.')
            end if
        end do

        if (qbmm .and. .not. polytropic) then
            do i = 1, nb
                do r = 1, nnode
                    write (file_num, '(I0)') sys_size + r + (i - 1)*nnode
                    file_loc = trim(t_step_dir) // '/pb' // trim(file_num) // '.dat'
                    inquire (FILE=trim(file_loc), EXIST=file_check)

                    if (file_check) then
                        open (1, FILE=trim(file_loc), form='unformatted', STATUS='old', ACTION='read')
                        read (1) pb%sf(:,:,:,r, i)
                        close (1)
                    else
                        call s_mpi_abort('File pb' // trim(file_num) // '.dat is missing in ' // trim(t_step_dir) // '. Exiting.')
                    end if
                end do
            end do

            do i = 1, nb
                do r = 1, nnode
                    write (file_num, '(I0)') sys_size + r + (i - 1)*nnode
                    file_loc = trim(t_step_dir) // '/mv' // trim(file_num) // '.dat'
                    inquire (FILE=trim(file_loc), EXIST=file_check)

                    if (file_check) then
                        open (1, FILE=trim(file_loc), form='unformatted', STATUS='old', ACTION='read')
                        read (1) mv%sf(:,:,:,r, i)
                        close (1)
                    else
                        call s_mpi_abort('File mv' // trim(file_num) // '.dat is missing in ' // trim(t_step_dir) // '. Exiting.')
                    end if
                end do
            end do
        end if

        ! Since the preexisting grid and initial condition data files have been read in, the directory associated with the rank of
        ! the local process may be cleaned out to make room for new pre-process data. In addition, the time-step folder that will
        ! contain the new grid and initial condition data are also generated.
        call s_delete_directory(trim(proc_rank_dir))
        call s_create_directory(trim(proc_rank_dir) // '/0')

    end subroutine s_read_serial_ic_data_files

    !> Cell-boundary data are checked for consistency by looking at the (non-)uniform cell-width distributions for all the active
    !! coordinate directions and making sure that all of the cell-widths are positively valued
    impure subroutine s_read_parallel_grid_data_files

#ifdef MFC_MPI
        real(wp), allocatable, dimension(:)  :: x_cb_glb, y_cb_glb, z_cb_glb
        integer                              :: ifile, ierr, data_size
        integer, dimension(MPI_STATUS_SIZE)  :: status
        character(LEN=path_len + 2*name_len) :: file_loc
        logical                              :: file_exist

        allocate (x_cb_glb(-1:m_glb))
        allocate (y_cb_glb(-1:n_glb))
        allocate (z_cb_glb(-1:p_glb))

        file_loc = trim(case_dir) // '/restart_data' // trim(mpiiofs) // 'x_cb.dat'
        inquire (FILE=trim(file_loc), EXIST=file_exist)

        if (file_exist) then
            data_size = m_glb + 2
            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
            call MPI_FILE_READ_ALL(ifile, x_cb_glb, data_size, mpi_p, status, ierr)
            call MPI_FILE_CLOSE(ifile, ierr)
        else
            call s_mpi_abort('File ' // trim(file_loc) // ' is missing. Exiting. ')
        end if

        x_cb(-1:m) = x_cb_glb((start_idx(1) - 1):(start_idx(1) + m))
        x_cc(0:m) = (x_cb(0:m) + x_cb(-1:(m - 1)))/2._wp
        dx = minval(x_cb(0:m) - x_cb(-1:(m - 1)))
        if (num_procs > 1) call s_mpi_reduce_min(dx)
        x_domain%beg = x_cb(-1)
        x_domain%end = x_cb(m)

        if (n > 0) then
            file_loc = trim(case_dir) // '/restart_data' // trim(mpiiofs) // 'y_cb.dat'
            inquire (FILE=trim(file_loc), EXIST=file_exist)

            if (file_exist) then
                data_size = n_glb + 2
                call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
                call MPI_FILE_READ_ALL(ifile, y_cb_glb, data_size, mpi_p, status, ierr)
                call MPI_FILE_CLOSE(ifile, ierr)
            else
                call s_mpi_abort('File ' // trim(file_loc) // ' is missing. Exiting. ')
            end if

            y_cb(-1:n) = y_cb_glb((start_idx(2) - 1):(start_idx(2) + n))
            y_cc(0:n) = (y_cb(0:n) + y_cb(-1:(n - 1)))/2._wp
            dy = minval(y_cb(0:n) - y_cb(-1:(n - 1)))
            if (num_procs > 1) call s_mpi_reduce_min(dy)
            y_domain%beg = y_cb(-1)
            y_domain%end = y_cb(n)

            if (p > 0) then
                file_loc = trim(case_dir) // '/restart_data' // trim(mpiiofs) // 'z_cb.dat'
                inquire (FILE=trim(file_loc), EXIST=file_exist)

                if (file_exist) then
                    data_size = p_glb + 2
                    call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
                    call MPI_FILE_READ_ALL(ifile, z_cb_glb, data_size, mpi_p, status, ierr)
                    call MPI_FILE_CLOSE(ifile, ierr)
                else
                    call s_mpi_abort('File ' // trim(file_loc) // ' is missing. Exiting. ')
                end if

                z_cb(-1:p) = z_cb_glb((start_idx(3) - 1):(start_idx(3) + p))
                z_cc(0:p) = (z_cb(0:p) + z_cb(-1:(p - 1)))/2._wp
                dz = minval(z_cb(0:p) - z_cb(-1:(p - 1)))
                if (num_procs > 1) call s_mpi_reduce_min(dz)
                z_domain%beg = z_cb(-1)
                z_domain%end = z_cb(p)
            end if
        end if

        deallocate (x_cb_glb, y_cb_glb, z_cb_glb)
#endif

    end subroutine s_read_parallel_grid_data_files

    !> The goal of this subroutine is to read in any preexisting initial condition data files so that they may be used by the
    !! pre-process as a starting point in the creation of an all new initial condition.
    impure subroutine s_read_parallel_ic_data_files(q_cons_vf_in)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf_in

#ifdef MFC_MPI
        integer                              :: ifile, ierr, data_size
        integer, dimension(MPI_STATUS_SIZE)  :: status
        integer(KIND=MPI_OFFSET_KIND)        :: disp
        integer(KIND=MPI_OFFSET_KIND)        :: m_MOK, n_MOK, p_MOK
        integer(KIND=MPI_OFFSET_KIND)        :: WP_MOK, var_MOK
        integer(KIND=MPI_OFFSET_KIND)        :: MOK
        character(LEN=path_len + 2*name_len) :: file_loc
        logical                              :: file_exist
        integer                              :: i

        if (cfl_adap_dt) then
            write (file_loc, '(I0,A)') n_start, '.dat'
        else
            write (file_loc, '(I0,A)') t_step_start, '.dat'
        end if
        file_loc = trim(restart_dir) // trim(mpiiofs) // trim(file_loc)
        inquire (FILE=trim(file_loc), EXIST=file_exist)

        if (file_exist) then
            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

            call s_initialize_mpi_data(q_cons_vf_in)

            data_size = (m + 1)*(n + 1)*(p + 1)

            ! Resize some integers so MPI can read even the biggest files
            m_MOK = int(m_glb + 1, MPI_OFFSET_KIND)
            n_MOK = int(n_glb + 1, MPI_OFFSET_KIND)
            p_MOK = int(p_glb + 1, MPI_OFFSET_KIND)
            WP_MOK = int(storage_size(0._stp)/8, MPI_OFFSET_KIND)
            MOK = int(1._wp, MPI_OFFSET_KIND)

            do i = 1, sys_size
                var_MOK = int(i, MPI_OFFSET_KIND)

                ! Initial displacement to skip at beginning of file
                disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, MPI_IO_DATA%view(i), 'native', mpi_info_int, ierr)
                call MPI_FILE_READ(ifile, MPI_IO_DATA%var(i)%sf, data_size, mpi_p, status, ierr)
            end do

            if (qbmm .and. .not. polytropic) then
                do i = sys_size + 1, sys_size + 2*nb*nnode
                    var_MOK = int(i, MPI_OFFSET_KIND)

                    ! Initial displacement to skip at beginning of file
                    disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                    call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, MPI_IO_DATA%view(i), 'native', mpi_info_int, ierr)
                    call MPI_FILE_READ(ifile, MPI_IO_DATA%var(i)%sf, data_size, mpi_p, status, ierr)
                end do
            end if

            call s_mpi_barrier()

            call MPI_FILE_CLOSE(ifile, ierr)
        else
            call s_mpi_abort('File ' // trim(file_loc) // ' is missing. Exiting. ')
        end if

        call s_mpi_barrier()
#endif

    end subroutine s_read_parallel_ic_data_files

    !> Initialize all pre-process modules, allocate data structures, and set I/O procedure pointers.
    impure subroutine s_initialize_modules

        call s_initialize_global_parameters_module()
        if (bubbles_euler .or. bubbles_lagrange) then
            call s_initialize_bubbles_model()
        end if
        call s_initialize_mpi_common_module()
        call s_initialize_data_output_module()
        call s_initialize_variables_conversion_module()
        call s_initialize_grid_module()
        call s_initialize_initial_condition_module()
        call s_initialize_perturbation_module()
        call s_initialize_assign_variables_module()
        call s_initialize_boundary_common_module()
        if (relax) call s_initialize_phasechange_module()

        ! Create the D directory if it doesn't exit, to store the serial data files
        call s_create_directory('D')

        if (parallel_io .neqv. .true.) then
            s_generate_grid => s_generate_serial_grid
            s_read_grid_data_files => s_read_serial_grid_data_files
            s_read_ic_data_files => s_read_serial_ic_data_files
            s_write_data_files => s_write_serial_data_files
        else
            s_generate_grid => s_generate_parallel_grid
            s_read_grid_data_files => s_read_parallel_grid_data_files
            s_read_ic_data_files => s_read_parallel_ic_data_files
            s_write_data_files => s_write_parallel_data_files
        end if

    end subroutine s_initialize_modules

    !> Read an existing grid from data files or generate a new grid from user inputs.
    impure subroutine s_read_grid()

        if (old_grid) then
            call s_read_grid_data_files()
            call s_check_grid_data_files()
        else
            if (parallel_io .neqv. .true.) then
                call s_generate_grid()
            else
                if (proc_rank == 0) call s_generate_grid()
                call s_mpi_barrier()
                call s_read_grid_data_files()
                call s_check_grid_data_files()
            end if
        end if

    end subroutine s_read_grid

    !> Generate or read the initial condition, apply relaxation if needed, and write output data files.
    impure subroutine s_apply_initial_condition(start, finish)

        real(wp), intent(inout) :: start, finish
        integer                 :: j, k, l
        real(wp)                :: r2

        call cpu_time(start)

        if (old_ic) call s_read_ic_data_files(q_cons_vf)

        call s_generate_initial_condition()

        ! hard-coded psi
        if (hyper_cleaning) then
            @:ASSERT(eqn_idx%psi > 0, "hyper_cleaning requires eqn_idx%psi to be set")
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        r2 = x_cc(j)**2
                        if (n > 0) r2 = r2 + y_cc(k)**2
                        if (p > 0) r2 = r2 + z_cc(l)**2
                        q_cons_vf(eqn_idx%psi)%sf(j, k, l) = 1.0e-2_wp*exp(-r2/(2.0_wp*0.05_wp**2))
                        q_prim_vf(eqn_idx%psi)%sf(j, k, l) = q_cons_vf(eqn_idx%psi)%sf(j, k, l)
                    end do
                end do
            end do
        end if

        if (relax) then
            if (proc_rank == 0) then
                print *, 'initial condition might have been altered due to enforcement of pTg-equilibrium (relax = "T" activated)'
            end if

            call s_infinite_relaxation_k(q_cons_vf)
        end if

        if (chemistry) then
            call s_write_data_files(q_cons_vf, q_prim_vf, bc_type, q_T_sf)
        else
            call s_write_data_files(q_cons_vf, q_prim_vf, bc_type)
        end if

        call cpu_time(finish)

    end subroutine s_apply_initial_condition

    !> Gather processor timing data and write elapsed wall-clock time to a summary file.
    impure subroutine s_save_data(proc_time, time_avg, time_final, file_exists)

        real(wp), dimension(:), intent(inout) :: proc_time
        real(wp), intent(inout)               :: time_avg, time_final
        logical, intent(inout)                :: file_exists

        call s_mpi_barrier()

        if (num_procs > 1) then
            call mpi_bcast_time_step_values(proc_time, time_avg)
        end if

        if (proc_rank == 0) then
            time_final = 0._wp
            if (num_procs == 1) then
                time_final = time_avg
                print *, "Elapsed Time", time_final
            else
                time_final = maxval(proc_time)
                print *, "Elapsed Time", time_final
            end if
            inquire (FILE='pre_time_data.dat', EXIST=file_exists)
            if (file_exists) then
                open (1, file='pre_time_data.dat', position='append', status='old')
                write (1, *) num_procs, time_final
                close (1)
            else
                open (1, file='pre_time_data.dat', status='new')
                write (1, *) num_procs, time_final
                close (1)
            end if
        end if

    end subroutine s_save_data

    !> Initialize MPI, read and validate user inputs on rank 0, and decompose the computational domain.
    impure subroutine s_initialize_mpi_domain

        call s_mpi_initialize()

        if (proc_rank == 0) then
            call s_assign_default_values_to_user_inputs()
            call s_read_input_file()
            call s_check_input_file()

            print '(" Pre-processing a ", I0, "x", I0, "x", I0, " case on ", I0, " rank(s)")', m, n, p, num_procs
        end if

        ! Broadcasting the user inputs to all of the processors and performing the parallel computational domain decomposition.
        ! Neither procedure has to be carried out if pre-process is in fact not truly executed in parallel.
        call s_mpi_bcast_user_inputs()
        call s_initialize_parallel_io()
        call s_mpi_decompose_computational_domain()

    end subroutine s_initialize_mpi_domain

    !> Finalize all pre-process modules, deallocate resources, and shut down MPI.
    impure subroutine s_finalize_modules

        s_generate_grid => null()
        s_read_grid_data_files => null()
        s_read_ic_data_files => null()
        s_write_data_files => null()

        call s_finalize_mpi_common_module()
        call s_finalize_grid_module()
        call s_finalize_variables_conversion_module()
        call s_finalize_data_output_module()
        call s_finalize_global_parameters_module()
        call s_finalize_assign_variables_module()
        call s_finalize_perturbation_module()
        call s_finalize_boundary_common_module()
        if (relax) call s_finalize_relaxation_solver_module()
        call s_finalize_initial_condition_module()
        call s_mpi_finalize()

    end subroutine s_finalize_modules

end module m_start_up
