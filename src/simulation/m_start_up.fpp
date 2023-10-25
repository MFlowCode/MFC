!>
!! @file m_start_up.f90
!! @brief Contains module m_start_up

#:include 'case.fpp'

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

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    use m_weno                 !< Weighted and essentially non-oscillatory (WENO)
                               !! schemes for spatial reconstruction of variables

    use m_riemann_solvers      !< Exact and approximate Riemann problem solvers

    use m_cbc                  !< Characteristic boundary conditions (CBC)

    use m_monopole             !< Monopole calculations

    use m_rhs                  !< Right-hand-side (RHS) evaluation procedures

    use m_data_output          !< Run-time info & solution data output procedures

    use m_time_steppers        !< Time-stepping algorithms

    use m_qbmm                 !< Quadrature MOM

    use m_derived_variables     !< Procedures used to compute quantites derived
                                !! from the conservative and primitive variables

    use m_hypoelastic

    use m_viscous

    use m_bubbles

    use ieee_arithmetic

#ifdef MFC_OpenACC
    use openacc
#endif

    use m_nvtx

    use m_helper

    use m_compile_specific

    use m_checker
    ! ==========================================================================

    implicit none

    private; public ::   s_read_input_file, &
                         s_check_input_file, &
                         s_read_data_files, &
                         s_read_serial_data_files, &
                         s_read_parallel_data_files, &
                         s_populate_grid_variables_buffers, &
                         s_initialize_internal_energy_equations, &
                         s_initialize_modules, s_initialize_gpu_vars, &
                         s_initialize_mpi_domain, s_finalize_modules, &
                         s_perform_time_step, s_save_data, &
                         s_save_performance_metrics

    abstract interface ! ===================================================

        !! @param q_cons_vf  Conservative variables
        subroutine s_read_abstract_data_files(q_cons_vf) ! -----------

            import :: scalar_field, sys_size, pres_field

            type(scalar_field), &
                dimension(sys_size), &
                intent(INOUT) :: q_cons_vf

        end subroutine s_read_abstract_data_files ! -----------------

    end interface ! ========================================================

    type(scalar_field), allocatable, dimension(:) :: grad_x_vf, grad_y_vf, grad_z_vf, norm_vf

    procedure(s_read_abstract_data_files), pointer :: s_read_data_files => null()

contains

    !>  The purpose of this procedure is to first verify that an
        !!      input file has been made available by the user. Provided
        !!      that this is so, the input file is then read in.
    subroutine s_read_input_file() ! ---------------------------------------

        ! Relative path to the input file provided by the user
        character(LEN=name_len) :: file_path = './simulation.inp'

        logical :: file_exist !<
            !! Logical used to check the existence of the input file

        integer :: iostatus
            !! Integer to check iostat of file read

        ! Namelist of the global parameters which may be specified by user
        namelist /user_inputs/ case_dir, run_time_info, m, n, p, dt, &
            t_step_start, t_step_stop, t_step_save, &
            model_eqns, num_fluids, adv_alphan, &
            mpp_lim, time_stepper,  weno_eps, weno_flat, &
            riemann_flat, cu_mpi, cu_tensor, &
            mapped_weno, mp_weno, weno_avg, &
            riemann_solver, wave_speeds, avg_state, &
            bc_x, bc_y, bc_z, &
            hypoelasticity, &
            fluid_pp, probe_wrt, prim_vars_wrt, &
            fd_order, probe, num_probes, t_step_old, &
            alt_soundspeed, mixture_err, weno_Re_flux, &
            null_weights, precision, parallel_io, cyl_coord, &
            rhoref, pref, bubbles, bubble_model, &
            R0ref, &
#:if not MFC_CASE_OPTIMIZATION
            nb, weno_order, &
#:endif
            Ca, Web, Re_inv, &
            monopole, mono, num_mono, &
            polytropic, thermal, &
            integral, integral_wrt, num_integrals, &
            polydisperse, poly_sigma, qbmm, &
            R0_type

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
                call s_mpi_abort('Invalid line in simulation.inp. It is '// &
                    'likely due to a datatype mismatch. Exiting ...')
            end if

            close (1)

            ! Store m,n,p into global m,n,p
            m_glb = m
            n_glb = n
            p_glb = p

        else
            call s_mpi_abort(trim(file_path)//' is missing. Exiting ...')
        end if

    end subroutine s_read_input_file ! -------------------------------------

    !> The goal of this procedure is to verify that each of the
    !!      user provided inputs is valid and that their combination
    !!      consitutes a meaningful configuration for the simulation.
    subroutine s_check_input_file() ! --------------------------------------

        ! Relative path to the current directory file in the case directory
        character(LEN=path_len) :: file_path

        ! Logical used to check the existence of the current directory file
        logical :: file_exist

        ! Generic loop iterators
        integer :: i, j

        ! Logistics ========================================================
        file_path = trim(case_dir)//'/.'

        call my_inquire(file_path, file_exist)

        if (file_exist .neqv. .true.) then
            call s_mpi_abort(trim(file_path)//' is missing. Exiting ...')
        end if
        ! ==================================================================

        call s_check_inputs()

    end subroutine s_check_input_file ! ------------------------------------

        !!              initial condition and grid data files. The cell-average
        !!              conservative variables constitute the former, while the
        !!              cell-boundary locations in x-, y- and z-directions make
        !!              up the latter. This procedure also calculates the cell-
        !!              width distributions from the cell-boundary locations.
        !! @param q_cons_vf Cell-averaged conservative variables
    subroutine s_read_serial_data_files(q_cons_vf) ! ------------------------------

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
        write (t_step_dir, '(A,I0,A,I0)') &
            trim(case_dir)//'/p_all/p', proc_rank, '/', t_step_start

        file_path = trim(t_step_dir)//'/.'
        call my_inquire(file_path, file_exist)

        if (file_exist .neqv. .true.) then 
            call s_mpi_abort(trim(file_path)//' is missing. Exiting ...')
        end if

        ! Cell-boundary Locations in x-direction ===========================
        file_path = trim(t_step_dir)//'/x_cb.dat'

        inquire (FILE=trim(file_path), EXIST=file_exist)

        if (file_exist) then
            open (2, FILE=trim(file_path), &
                  FORM='unformatted', &
                  ACTION='read', &
                  STATUS='old')
            read (2) x_cb(-1:m); close (2)
        else
            call s_mpi_abort(trim(file_path)//' is missing. Exiting ...')
        end if

        dx(0:m) = x_cb(0:m) - x_cb(-1:m - 1)
        x_cc(0:m) = x_cb(-1:m - 1) + dx(0:m)/2d0
        ! ==================================================================

        ! Cell-boundary Locations in y-direction ===========================
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
                call s_mpi_abort(trim(file_path)//' is missing. Exiting ...')
            end if

            dy(0:n) = y_cb(0:n) - y_cb(-1:n - 1)
            y_cc(0:n) = y_cb(-1:n - 1) + dy(0:n)/2d0

        end if
        ! ==================================================================

        ! Cell-boundary Locations in z-direction ===========================
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
                call s_mpi_abort(trim(file_path)//' is missing. Exiting ...')
            end if

            dz(0:p) = z_cb(0:p) - z_cb(-1:p - 1)
            z_cc(0:p) = z_cb(-1:p - 1) + dz(0:p)/2d0

        end if
        ! ==================================================================

        ! Cell-average Conservative Variables ==============================
        if ((bubbles .neqv. .true.) .and. (hypoelasticity .neqv. .true.)) then
            do i = 1, adv_idx%end
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
                    call s_mpi_abort(trim(file_path)//' is missing. Exiting ...')
                end if
            end do
        else
            !make sure to read bubble variables
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
                    call s_mpi_abort(trim(file_path)//' is missing. Exiting ...')
                end if
            end do
            !Read pb and mv for non-polytropic qbmm
            if(qbmm .and. .not. polytropic) then
                do i = 1, nb
                    do r = 1, nnode
                        write (file_path, '(A,I0,A)') &
                            trim(t_step_dir)//'/pb', sys_size + (i-1)*nnode + r, '.dat'
                        inquire (FILE=trim(file_path), EXIST=file_exist)
                        if (file_exist) then
                            open (2, FILE=trim(file_path), &
                                  FORM='unformatted', &
                                  ACTION='read', &
                                  STATUS='old')
                            read (2) pb_ts(1)%sf(0:m, 0:n, 0:p, r, i); close (2)
                        else
                            call s_mpi_abort(trim(file_path)//' is missing. Exiting ...')
                        end if
                    end do
                end do                
                do i = 1, nb
                    do r = 1, nnode
                        write (file_path, '(A,I0,A)') &
                            trim(t_step_dir)//'/mv', sys_size + (i-1)*nnode + r , '.dat'
                        inquire (FILE=trim(file_path), EXIST=file_exist)
                        if (file_exist) then
                            open (2, FILE=trim(file_path), &
                                  FORM='unformatted', &
                                  ACTION='read', &
                                  STATUS='old')
                            read (2) mv_ts(1)%sf(0:m, 0:n, 0:p, r, i); close (2)
                        else
                            call s_mpi_abort(trim(file_path)//' is missing. Exiting ...')
                        end if
                    end do
                end do 
            end if
        end if
        ! ==================================================================

    end subroutine s_read_serial_data_files ! -------------------------------------

        !! @param q_cons_vf Conservative variables
    subroutine s_read_parallel_data_files(q_cons_vf) ! ---------------------------

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: q_cons_vf

#ifdef MFC_MPI

        real(kind(0d0)), allocatable, dimension(:) :: x_cb_glb, y_cb_glb, z_cb_glb

        integer :: ifile, ierr, data_size
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(KIND=MPI_OFFSET_KIND) :: disp
        integer(KIND=MPI_OFFSET_KIND) :: m_MOK, n_MOK, p_MOK
        integer(KIND=MPI_OFFSET_KIND) :: WP_MOK, var_MOK, str_MOK
        integer(KIND=MPI_OFFSET_KIND) :: NVARS_MOK
        integer(KIND=MPI_OFFSET_KIND) :: MOK

        character(LEN=path_len + 2*name_len) :: file_loc
        logical :: file_exist

        integer :: i

        allocate (x_cb_glb(-1:m_glb))
        allocate (y_cb_glb(-1:n_glb))
        allocate (z_cb_glb(-1:p_glb))

        ! Read in cell boundary locations in x-direction
        file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//'x_cb.dat'
        inquire (FILE=trim(file_loc), EXIST=file_exist)

        if (file_exist) then
            data_size = m_glb + 2
            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
            call MPI_FILE_READ(ifile, x_cb_glb, data_size, MPI_DOUBLE_PRECISION, status, ierr)
            call MPI_FILE_CLOSE(ifile, ierr)
        else
            call s_mpi_abort('File '//trim(file_loc)//' is missing. Exiting...')
        end if

        ! Assigning local cell boundary locations
        x_cb(-1:m) = x_cb_glb((start_idx(1) - 1):(start_idx(1) + m))
        ! Computing the cell width distribution
        dx(0:m) = x_cb(0:m) - x_cb(-1:m - 1)
        ! Computing the cell center locations
        x_cc(0:m) = x_cb(-1:m - 1) + dx(0:m)/2d0

        if (n > 0) then
            ! Read in cell boundary locations in y-direction
            file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//'y_cb.dat'
            inquire (FILE=trim(file_loc), EXIST=file_exist)

            if (file_exist) then
                data_size = n_glb + 2
                call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
                call MPI_FILE_READ(ifile, y_cb_glb, data_size, MPI_DOUBLE_PRECISION, status, ierr)
                call MPI_FILE_CLOSE(ifile, ierr)
            else
                call s_mpi_abort('File '//trim(file_loc)//' is missing. Exiting...')
            end if

            ! Assigning local cell boundary locations
            y_cb(-1:n) = y_cb_glb((start_idx(2) - 1):(start_idx(2) + n))
            ! Computing the cell width distribution
            dy(0:n) = y_cb(0:n) - y_cb(-1:n - 1)
            ! Computing the cell center locations
            y_cc(0:n) = y_cb(-1:n - 1) + dy(0:n)/2d0

            if (p > 0) then
                ! Read in cell boundary locations in z-direction
                file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//'z_cb.dat'
                inquire (FILE=trim(file_loc), EXIST=file_exist)

                if (file_exist) then
                    data_size = p_glb + 2
                    call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
                    call MPI_FILE_READ(ifile, z_cb_glb, data_size, MPI_DOUBLE_PRECISION, status, ierr)
                    call MPI_FILE_CLOSE(ifile, ierr)
                else
                    call s_mpi_abort( 'File '//trim(file_loc)//'is missing. Exiting...')
                end if

                ! Assigning local cell boundary locations
                z_cb(-1:p) = z_cb_glb((start_idx(3) - 1):(start_idx(3) + p))
                ! Computing the cell width distribution
                dz(0:p) = z_cb(0:p) - z_cb(-1:p - 1)
                ! Computing the cell center locations
                z_cc(0:p) = z_cb(-1:p - 1) + dz(0:p)/2d0

            end if
        end if

        ! Open the file to read conservative variables
        write (file_loc, '(I0,A)') t_step_start, '.dat'
        file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
        inquire (FILE=trim(file_loc), EXIST=file_exist)

        if (file_exist) then
            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

            ! Initialize MPI data I/O


            call s_initialize_mpi_data(q_cons_vf)

            ! Size of local arrays
            data_size = (m + 1)*(n + 1)*(p + 1)

            ! Resize some integers so MPI can read even the biggest file
            m_MOK = int(m_glb + 1, MPI_OFFSET_KIND)
            n_MOK = int(n_glb + 1, MPI_OFFSET_KIND)
            p_MOK = int(p_glb + 1, MPI_OFFSET_KIND)
            WP_MOK = int(8d0, MPI_OFFSET_KIND)
            MOK = int(1d0, MPI_OFFSET_KIND)
            str_MOK = int(name_len, MPI_OFFSET_KIND)
            NVARS_MOK = int(sys_size, MPI_OFFSET_KIND)

            ! Read the data for each variable
            if (bubbles .or. hypoelasticity) then

                do i = 1, sys_size!adv_idx%end
                    var_MOK = int(i, MPI_OFFSET_KIND)
                    ! Initial displacement to skip at beginning of file
                    disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                    call MPI_FILE_SET_VIEW(ifile, disp, MPI_DOUBLE_PRECISION, MPI_IO_DATA%view(i), &
                                           'native', mpi_info_int, ierr)
                    call MPI_FILE_READ(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                       MPI_DOUBLE_PRECISION, status, ierr)
                end do
                !Read pb and mv for non-polytropic qbmm
                if(qbmm .and. .not. polytropic) then
                    do i = sys_size + 1, sys_size + 2*nb*nnode
                        var_MOK = int(i, MPI_OFFSET_KIND)
                        ! Initial displacement to skip at beginning of file
                        disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                        call MPI_FILE_SET_VIEW(ifile, disp, MPI_DOUBLE_PRECISION, MPI_IO_DATA%view(i), &
                                               'native', mpi_info_int, ierr)
                        call MPI_FILE_READ(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                           MPI_DOUBLE_PRECISION, status, ierr)
                    end do
                end if
            else
                do i = 1, adv_idx%end
                    var_MOK = int(i, MPI_OFFSET_KIND)

                    ! Initial displacement to skip at beginning of file
                    disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                    call MPI_FILE_SET_VIEW(ifile, disp, MPI_DOUBLE_PRECISION, MPI_IO_DATA%view(i), &
                                           'native', mpi_info_int, ierr)
                    call MPI_FILE_READ(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                       MPI_DOUBLE_PRECISION, status, ierr)
                end do
            end if

            call s_mpi_barrier()

            call MPI_FILE_CLOSE(ifile, ierr)
        else
            call s_mpi_abort('File '//trim(file_loc)//' is missing. Exiting...')
        end if

        deallocate (x_cb_glb, y_cb_glb, z_cb_glb)

#endif

    end subroutine s_read_parallel_data_files ! -------------------------------

    !> The purpose of this subroutine is to populate the buffers
        !!          of the grid variables, which are constituted of the cell-
        !!          boundary locations and cell-width distributions, based on
        !!          the boundary conditions.
    subroutine s_populate_grid_variables_buffers() ! -----------------------

        integer :: i !< Generic loop iterator

        ! Population of Buffers in x-direction =============================

        ! Populating cell-width distribution buffer, at the beginning of the
        ! coordinate direction, based on the selected boundary condition. In
        ! order, these are the ghost-cell extrapolation, symmetry, periodic,
        ! and processor boundary conditions.
        if (bc_x%beg <= -3) then
            do i = 1, buff_size
                dx(-i) = dx(0)
            end do
        elseif (bc_x%beg == -2) then
            do i = 1, buff_size
                dx(-i) = dx(i - 1)
            end do
        elseif (bc_x%beg == -1) then
            do i = 1, buff_size
                dx(-i) = dx(m - (i - 1))
            end do
        else
            call s_mpi_sendrecv_grid_variables_buffers(1, -1)
        end if

        ! Computing the cell-boundary locations buffer, at the beginning of
        ! the coordinate direction, from the cell-width distribution buffer
        do i = 1, buff_size
            x_cb(-1 - i) = x_cb(-i) - dx(-i)
        end do
        ! Computing the cell-center locations buffer, at the beginning of
        ! the coordinate direction, from the cell-width distribution buffer
        do i = 1, buff_size
            x_cc(-i) = x_cc(1 - i) - (dx(1 - i) + dx(-i))/2d0
        end do

        ! Populating the cell-width distribution buffer, at the end of the
        ! coordinate direction, based on desired boundary condition. These
        ! include, in order, ghost-cell extrapolation, symmetry, periodic,
        ! and processor boundary conditions.
        if (bc_x%end <= -3) then
            do i = 1, buff_size
                dx(m + i) = dx(m)
            end do
        elseif (bc_x%end == -2) then
            do i = 1, buff_size
                dx(m + i) = dx(m - (i - 1))
            end do
        elseif (bc_x%end == -1) then
            do i = 1, buff_size
                dx(m + i) = dx(i - 1)
            end do
        else
            call s_mpi_sendrecv_grid_variables_buffers(1, 1)
        end if

        ! Populating the cell-boundary locations buffer, at the end of the
        ! coordinate direction, from buffer of the cell-width distribution
        do i = 1, buff_size
            x_cb(m + i) = x_cb(m + (i - 1)) + dx(m + i)
        end do
        ! Populating the cell-center locations buffer, at the end of the
        ! coordinate direction, from buffer of the cell-width distribution
        do i = 1, buff_size
            x_cc(m + i) = x_cc(m + (i - 1)) + (dx(m + (i - 1)) + dx(m + i))/2d0
        end do

        ! END: Population of Buffers in x-direction ========================

        ! Population of Buffers in y-direction =============================

        ! Populating cell-width distribution buffer, at the beginning of the
        ! coordinate direction, based on the selected boundary condition. In
        ! order, these are the ghost-cell extrapolation, symmetry, periodic,
        ! and processor boundary conditions.
        if (n == 0) then
            return
        elseif (bc_y%beg <= -3 .and. bc_y%beg /= -14) then
            do i = 1, buff_size
                dy(-i) = dy(0)
            end do
        elseif (bc_y%beg == -2 .or. bc_y%beg == -14) then
            do i = 1, buff_size
                dy(-i) = dy(i - 1)
            end do
        elseif (bc_y%beg == -1) then
            do i = 1, buff_size
                dy(-i) = dy(n - (i - 1))
            end do
        else
            call s_mpi_sendrecv_grid_variables_buffers(2, -1)
        end if

        ! Computing the cell-boundary locations buffer, at the beginning of
        ! the coordinate direction, from the cell-width distribution buffer
        do i = 1, buff_size
            y_cb(-1 - i) = y_cb(-i) - dy(-i)
        end do
        ! Computing the cell-center locations buffer, at the beginning of
        ! the coordinate direction, from the cell-width distribution buffer
        do i = 1, buff_size
            y_cc(-i) = y_cc(1 - i) - (dy(1 - i) + dy(-i))/2d0
        end do

        ! Populating the cell-width distribution buffer, at the end of the
        ! coordinate direction, based on desired boundary condition. These
        ! include, in order, ghost-cell extrapolation, symmetry, periodic,
        ! and processor boundary conditions.
        if (bc_y%end <= -3) then
            do i = 1, buff_size
                dy(n + i) = dy(n)
            end do
        elseif (bc_y%end == -2) then
            do i = 1, buff_size
                dy(n + i) = dy(n - (i - 1))
            end do
        elseif (bc_y%end == -1) then
            do i = 1, buff_size
                dy(n + i) = dy(i - 1)
            end do
        else
            call s_mpi_sendrecv_grid_variables_buffers(2, 1)
        end if

        ! Populating the cell-boundary locations buffer, at the end of the
        ! coordinate direction, from buffer of the cell-width distribution
        do i = 1, buff_size
            y_cb(n + i) = y_cb(n + (i - 1)) + dy(n + i)
        end do
        ! Populating the cell-center locations buffer, at the end of the
        ! coordinate direction, from buffer of the cell-width distribution
        do i = 1, buff_size
            y_cc(n + i) = y_cc(n + (i - 1)) + (dy(n + (i - 1)) + dy(n + i))/2d0
        end do

        ! END: Population of Buffers in y-direction ========================

        ! Population of Buffers in z-direction =============================

        ! Populating cell-width distribution buffer, at the beginning of the
        ! coordinate direction, based on the selected boundary condition. In
        ! order, these are the ghost-cell extrapolation, symmetry, periodic,
        ! and processor boundary conditions.
        if (p == 0) then
            return
        elseif (bc_z%beg <= -3) then
            do i = 1, buff_size
                dz(-i) = dz(0)
            end do
        elseif (bc_z%beg == -2) then
            do i = 1, buff_size
                dz(-i) = dz(i - 1)
            end do
        elseif (bc_z%beg == -1) then
            do i = 1, buff_size
                dz(-i) = dz(p - (i - 1))
            end do
        else
            call s_mpi_sendrecv_grid_variables_buffers(3, -1)
        end if

        ! Computing the cell-boundary locations buffer, at the beginning of
        ! the coordinate direction, from the cell-width distribution buffer
        do i = 1, buff_size
            z_cb(-1 - i) = z_cb(-i) - dz(-i)
        end do
        ! Computing the cell-center locations buffer, at the beginning of
        ! the coordinate direction, from the cell-width distribution buffer
        do i = 1, buff_size
            z_cc(-i) = z_cc(1 - i) - (dz(1 - i) + dz(-i))/2d0
        end do

        ! Populating the cell-width distribution buffer, at the end of the
        ! coordinate direction, based on desired boundary condition. These
        ! include, in order, ghost-cell extrapolation, symmetry, periodic,
        ! and processor boundary conditions.
        if (bc_z%end <= -3) then
            do i = 1, buff_size
                dz(p + i) = dz(p)
            end do
        elseif (bc_z%end == -2) then
            do i = 1, buff_size
                dz(p + i) = dz(p - (i - 1))
            end do
        elseif (bc_z%end == -1) then
            do i = 1, buff_size
                dz(p + i) = dz(i - 1)
            end do
        else
            call s_mpi_sendrecv_grid_variables_buffers(3, 1)
        end if

        ! Populating the cell-boundary locations buffer, at the end of the
        ! coordinate direction, from buffer of the cell-width distribution
        do i = 1, buff_size
            z_cb(p + i) = z_cb(p + (i - 1)) + dz(p + i)
        end do
        ! Populating the cell-center locations buffer, at the end of the
        ! coordinate direction, from buffer of the cell-width distribution
        do i = 1, buff_size
            z_cc(p + i) = z_cc(p + (i - 1)) + (dz(p + (i - 1)) + dz(p + i))/2d0
        end do

        ! END: Population of Buffers in z-direction ========================

    end subroutine s_populate_grid_variables_buffers ! ---------------------

    !> The purpose of this procedure is to initialize the
        !!      values of the internal-energy equations of each phase
        !!      from the mass of each phase, the mixture momentum and
        !!      mixture-total-energy equations.
        !! @param v_vf conservative variables
    subroutine s_initialize_internal_energy_equations(v_vf) !---------------

        type(scalar_field), dimension(sys_size), intent(INOUT) :: v_vf
        real(kind(0d0)) :: rho
        real(kind(0d0)) :: dyn_pres
        real(kind(0d0)) :: gamma
        real(kind(0d0)) :: pi_inf
        real(kind(0d0)), dimension(2) :: Re
        real(kind(0d0)) :: pres

        integer :: i, j, k, l

        do j = 0, m
            do k = 0, n
                do l = 0, p

                    call s_convert_to_mixture_variables(v_vf, j, k, l, rho, gamma, pi_inf, Re)

                    dyn_pres = 0d0
                    do i = mom_idx%beg, mom_idx%end
                        dyn_pres = dyn_pres + 5d-1*v_vf(i)%sf(j, k, l)*v_vf(i)%sf(j, k, l) &
                                   /max(rho, sgm_eps)
                    end do

                    call s_compute_pressure(v_vf(E_idx)%sf(j, k, l), 0d0, &
                        dyn_pres, pi_inf, gamma, rho, pres)

                    do i = 1, num_fluids
                        v_vf(i + internalEnergies_idx%beg - 1)%sf(j, k, l) = v_vf(i + adv_idx%beg - 1)%sf(j, k, l)* &
                                                                             (fluid_pp(i)%gamma*pres + fluid_pp(i)%pi_inf)
                    end do

                end do
            end do
        end do

    end subroutine s_initialize_internal_energy_equations !-----------------

    subroutine s_perform_time_step(t_step, time_avg, time_final, io_time_avg, io_time_final, proc_time, io_proc_time, file_exists, start, finish, nt)
        integer, intent(INOUT) :: t_step
        real(kind(0d0)), intent(INOUT) :: time_avg, time_final
        real(kind(0d0)), intent(INOUT) :: io_time_avg, io_time_final
        real(kind(0d0)),  dimension(:), intent(INOUT) :: proc_time
        real(kind(0d0)),  dimension(:), intent(INOUT) :: io_proc_time
        logical, intent(INOUT) :: file_exists
        real(kind(0d0)), intent(INOUT) :: start, finish
        integer, intent(INOUT) :: nt

        integer :: i, j, k, l

        if (proc_rank == 0) then
            print '(" ["I3"%]  Time step "I8" of "I0" @ t_step = "I0"")',                             &
                  int(ceiling(100d0*(real(t_step - t_step_start)/(t_step_stop - t_step_start + 1)))), &
                  t_step      - t_step_start + 1,                                                     &
                  t_step_stop - t_step_start + 1,                                                     &
                  t_step
        end if
        mytime = mytime + dt

        if (probe_wrt) then
            do i = 1, sys_size
                !$acc update host(q_cons_ts(1)%vf(i)%sf)
            end do
        end if

        call s_compute_derived_variables(t_step)

#ifdef DEBUG
        print *, 'Computed derived vars'
#endif

        ! Total-variation-diminishing (TVD) Runge-Kutta (RK) time-steppers
        if (time_stepper == 1) then
            call s_1st_order_tvd_rk(t_step, time_avg)
        elseif (time_stepper == 2) then
            call s_2nd_order_tvd_rk(t_step, time_avg)
        elseif (time_stepper == 3) then
            call s_3rd_order_tvd_rk(t_step, time_avg)
        end if

        ! Time-stepping loop controls
        if ((mytime + dt) >= finaltime) dt = finaltime - mytime 
        t_step = t_step + 1

    end subroutine s_perform_time_step

    subroutine s_save_performance_metrics(t_step, time_avg, time_final, io_time_avg, io_time_final, proc_time, io_proc_time, file_exists, start, finish, nt)

        integer, intent(INOUT) :: t_step
        real(kind(0d0)), intent(INOUT) :: time_avg, time_final
        real(kind(0d0)), intent(INOUT) :: io_time_avg, io_time_final
        real(kind(0d0)),  dimension(:), intent(INOUT) :: proc_time
        real(kind(0d0)),  dimension(:), intent(INOUT) :: io_proc_time
        logical, intent(INOUT) :: file_exists
        real(kind(0d0)), intent(INOUT) :: start, finish
        integer, intent(INOUT) :: nt

        call s_mpi_barrier()

        if (num_procs > 1) then
            call mpi_bcast_time_step_values(proc_time, time_avg)

            call mpi_bcast_time_step_values(io_proc_time, io_time_avg)
        end if

        if (proc_rank == 0) then
            time_final = 0d0
            io_time_final = 0d0
            if (num_procs == 1) then
                time_final = time_avg
                io_time_final = io_time_avg
                print *, "Final Time", time_final
            else
                time_final = maxval(proc_time)
                io_time_final = maxval(io_proc_time)
                print *, "Final Time", time_final
            end if
            inquire (FILE='time_data.dat', EXIST=file_exists)
            if (file_exists) then
                open (1, file='time_data.dat', position='append', status='old')
                write (1, *) num_procs, time_final
                close (1)
            else
                open (1, file='time_data.dat', status='new')
                write (1, *) num_procs, time_final
                close (1)
            end if

            inquire (FILE='io_time_data.dat', EXIST=file_exists)
            if (file_exists) then
                open (1, file='io_time_data.dat', position='append', status='old')
                write (1, *) num_procs, io_time_final
                close (1)
            else
                open (1, file='io_time_data.dat', status='new')
                write (1, *) num_procs, io_time_final
                close (1)
            end if

        end if

    end subroutine s_save_performance_metrics

    subroutine s_save_data(t_step, start, finish, io_time_avg, nt)
        real(kind(0d0)), intent(INOUT) ::  start, finish, io_time_avg
        integer, intent(INOUT) :: t_step, nt
        integer :: i, j, k, l

        if (mod(t_step - t_step_start, t_step_save) == 0 .or. t_step == t_step_stop) then

            call cpu_time(start)
            !  call nvtxStartRange("I/O")
            do i = 1, sys_size
                !$acc update host(q_cons_ts(1)%vf(i)%sf)
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            if(ieee_is_nan(q_cons_ts(1)%vf(i)%sf(j, k, l))) then
                                print *, "NaN(s) in timestep output.", j, k, l, i,  proc_rank, t_step, m, n, p                                
                                error stop "NaN(s) in timestep output."
                            end if
                        end do
                    end do
                end do
            end do

            if(qbmm .and. .not. polytropic) then
                !$acc update host(pb_ts(1)%sf)
                !$acc update host(mv_ts(1)%sf)
            end if

            call s_write_data_files(q_cons_ts(1)%vf, q_prim_vf, t_step)
            !  call nvtxEndRange
            call cpu_time(finish)
            nt = int((t_step - t_step_start)/(t_step_save))
            if (nt == 1) then
                io_time_avg = abs(finish - start)
            else
                io_time_avg = (abs(finish - start) + io_time_avg*(nt - 1))/nt
            end if
        end if

    end subroutine s_save_data    

    subroutine s_initialize_modules()
        call s_initialize_global_parameters_module()
        !Quadrature weights and nodes for polydisperse simulations
        if(bubbles .and. nb > 1 .and. R0_type == 1) then
            call s_simpson
        end if
        !Initialize variables for non-polytropic (Preston) model
        if(bubbles .and. .not. polytropic) then
            call s_initialize_nonpoly()
        end if
        !Initialize pb based on surface tension for qbmm (polytropic)
        if(qbmm .and. polytropic .and. Web /= dflt_real) then
            pb0 = pref + 2d0 * fluid_pp(1)%ss / (R0*R0ref)
            pb0 = pb0 / pref
            pref = 1d0                              
        end if

#if defined(MFC_OpenACC) && defined(MFC_MEMORY_DUMP)
        call acc_present_dump()
#endif

        call s_initialize_mpi_proxy_module()
        call s_initialize_variables_conversion_module()
        if (grid_geometry == 3) call s_initialize_fftw_module()
        call s_initialize_riemann_solvers_module()

        if(bubbles) call s_initialize_bubbles_module()

        if (qbmm) call s_initialize_qbmm_module()

#if defined(MFC_OpenACC) && defined(MFC_MEMORY_DUMP)
        call acc_present_dump()
#endif

        if (monopole) then
            call s_initialize_monopole_module()
        end if
        if (any(Re_size > 0)) then
            call s_initialize_viscous_module()
        end if
        call s_initialize_rhs_module()

#if defined(MFC_OpenACC) && defined(MFC_MEMORY_DUMP)
        call acc_present_dump()
#endif

        if (hypoelasticity) call s_initialize_hypoelastic_module()
        call s_initialize_data_output_module()
        call s_initialize_derived_variables_module()
        call s_initialize_time_steppers_module()

#if defined(MFC_OpenACC) && defined(MFC_MEMORY_DUMP)
        call acc_present_dump()
#endif

        ! Associate pointers for serial or parallel I/O
        if (parallel_io .neqv. .true.) then
            s_read_data_files => s_read_serial_data_files
            s_write_data_files => s_write_serial_data_files
        else
            s_read_data_files => s_read_parallel_data_files
            s_write_data_files => s_write_parallel_data_files
        end if

        ! Reading in the user provided initial condition and grid data
        call s_read_data_files(q_cons_ts(1)%vf)
        if (model_eqns == 3) call s_initialize_internal_energy_equations(q_cons_ts(1)%vf)

        ! Populating the buffers of the grid variables using the boundary conditions
        call s_populate_grid_variables_buffers()

        ! Computation of parameters, allocation of memory, association of pointers,
        ! and/or execution of any other tasks that are needed to properly configure
        ! the modules. The preparations below DO DEPEND on the grid being complete.
        call s_initialize_weno_module()

#if defined(MFC_OpenACC) && defined(MFC_MEMORY_DUMP)
        print *, "[MEM-INST] After: s_initialize_weno_module"
        call acc_present_dump()
#endif

        call s_initialize_cbc_module()

        call s_initialize_derived_variables()

    end subroutine s_initialize_modules

    subroutine s_initialize_mpi_domain()
        integer :: ierr
#ifdef MFC_OpenACC
        real(kind(0d0)) :: starttime, endtime
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
        devNum  = acc_get_num_devices(devtype)
        dev     = mod(local_rank, devNum)

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
            print '(" Simulating a "I0"x"I0"x"I0" case on "I0" rank(s)")', m, n, p, num_procs
        end if

        ! Broadcasting the user inputs to all of the processors and performing the
        ! parallel computational domain decomposition. Neither procedure has to be
        ! carried out if the simulation is in fact not truly executed in parallel.

        call s_mpi_bcast_user_inputs()
        call s_initialize_parallel_io()
        call s_mpi_decompose_computational_domain()

    end subroutine s_initialize_mpi_domain

    subroutine s_initialize_gpu_vars()
        integer :: i
        !Update GPU DATA
        !$acc update device(dt, dx, dy, dz, x_cc, y_cc, z_cc, x_cb, y_cb, z_cb)
        !$acc update device(sys_size, buff_size)
        !$acc update device(m, n, p)
        !$acc update device(momxb, momxe, bubxb, bubxe, advxb, advxe, contxb, contxe, strxb, strxe)
        do i = 1, sys_size
        !$acc update device(q_cons_ts(1)%vf(i)%sf)
        end do
        if(qbmm .and. .not. polytropic) then
        !$acc update device(pb_ts(1)%sf, mv_ts(1)%sf)
        end if
        !$acc update device(dt, sys_size, pref, rhoref, gamma_idx, pi_inf_idx, E_idx, alf_idx, stress_idx, mpp_lim, bubbles, hypoelasticity, alt_soundspeed, avg_state, num_fluids, model_eqns, num_dims, mixture_err, nb, weight, grid_geometry, cyl_coord, mapped_weno, mp_weno, weno_eps)
        !$acc update device(nb, R0ref, Ca, Web, Re_inv, weight, R0, V0, bubbles, polytropic, polydisperse, qbmm, R0_type, ptil, bubble_model, thermal, poly_sigma)
        !$acc update device(R_n, R_v, phi_vn, phi_nv, Pe_c, Tw, pv, M_n, M_v, k_n, k_v, pb0, mass_n0, mass_v0, Pe_T, Re_trans_T, Re_trans_c, Im_trans_T, Im_trans_c, omegaN , mul0, ss, gamma_v, mu_v, gamma_m, gamma_n, mu_n, gam)
        !$acc update device(monopole, num_mono)
    end subroutine s_initialize_gpu_vars


    subroutine s_finalize_modules()
        ! Disassociate pointers for serial and parallel I/O
        s_read_data_files => null()
        s_write_data_files => null()

        call s_finalize_time_steppers_module()
        call s_finalize_derived_variables_module()
        call s_finalize_data_output_module()
        call s_finalize_rhs_module()
        call s_finalize_cbc_module()
        call s_finalize_riemann_solvers_module()
        call s_finalize_weno_module()
        call s_finalize_variables_conversion_module()
        if (grid_geometry == 3) call s_finalize_fftw_module
        call s_finalize_mpi_proxy_module()
        call s_finalize_global_parameters_module()

        if (any(Re_size > 0)) then
            call s_finalize_viscous_module()
        end if

        ! Terminating MPI execution environment
        call s_mpi_finalize()
    end subroutine s_finalize_modules

end module m_start_up
