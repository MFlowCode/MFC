!! @file m_data_output.f90
!! @brief Contains module m_data_output

#:include 'macros.fpp'

!> @brief The primary purpose of this module is to output the grid and the
!!              conservative variables data at the chosen time-step interval. In
!!              addition, this module is also in charge of outputting a run-time
!!              information file which summarizes the time-dependent behavior !of
!!              the stability criteria. The latter include the inviscid Courant–
!!              Friedrichs–Lewy (ICFL), viscous CFL (VCFL), capillary CFL (CCFL)
!!              and cell Reynolds (Rc) numbers.
module m_data_output

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    use m_compile_specific

    use m_helper

    use m_helper_basic         !< Functions to compare floating point numbers

    use m_sim_helpers

    use m_delay_file_access

    use m_ibm

    use m_boundary_common

    implicit none

    private; 
    public :: s_initialize_data_output_module, &
              s_open_run_time_information_file, &
              s_open_com_files, &
              s_open_probe_files, &
              s_write_run_time_information, &
              s_write_data_files, &
              s_write_serial_data_files, &
              s_write_parallel_data_files, &
              s_write_com_files, &
              s_write_probe_files, &
              s_close_run_time_information_file, &
              s_close_com_files, &
              s_close_probe_files, &
              s_finalize_data_output_module

    real(wp), allocatable, dimension(:, :, :) :: icfl_sf  !< ICFL stability criterion
    real(wp), allocatable, dimension(:, :, :) :: vcfl_sf  !< VCFL stability criterion
    real(wp), allocatable, dimension(:, :, :) :: ccfl_sf  !< CCFL stability criterion
    real(wp), allocatable, dimension(:, :, :) :: Rc_sf  !< Rc stability criterion
    real(wp), public, allocatable, dimension(:, :) :: c_mass
    $:GPU_DECLARE(create='[icfl_sf,vcfl_sf,ccfl_sf,Rc_sf,c_mass]')

    real(wp) :: icfl_max_loc, icfl_max_glb !< ICFL stability extrema on local and global grids
    real(wp) :: vcfl_max_loc, vcfl_max_glb !< VCFL stability extrema on local and global grids
    real(wp) :: ccfl_max_loc, ccfl_max_glb !< CCFL stability extrema on local and global grids
    real(wp) :: Rc_min_loc, Rc_min_glb !< Rc   stability extrema on local and global grids
    $:GPU_DECLARE(create='[icfl_max_loc,icfl_max_glb,vcfl_max_loc,vcfl_max_glb]')
    $:GPU_DECLARE(create='[ccfl_max_loc,ccfl_max_glb,Rc_min_loc,Rc_min_glb]')

    !> @name ICFL, VCFL, CCFL and Rc stability criteria extrema over all the time-steps
    !> @{
    real(wp) :: icfl_max !< ICFL criterion maximum
    real(wp) :: vcfl_max !< VCFL criterion maximum
    real(wp) :: ccfl_max !< CCFL criterion maximum
    real(wp) :: Rc_min !< Rc criterion maximum
    !> @}

    type(scalar_field), allocatable, dimension(:) :: q_cons_temp

contains

    !> Write data files. Dispatch subroutine that replaces procedure pointer.
        !! @param q_cons_vf Conservative variables
        !! @param q_prim_vf Primitive variables
        !! @param t_step Current time step
    impure subroutine s_write_data_files(q_cons_vf, q_T_sf, q_prim_vf, t_step, bc_type, beta)

        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: q_cons_vf

        type(scalar_field), &
            intent(inout) :: q_T_sf

        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: q_prim_vf

        integer, intent(in) :: t_step

        type(scalar_field), &
            intent(inout), optional :: beta

        type(integer_field), &
            dimension(1:num_dims, -1:1), &
            intent(in) :: bc_type

        if (.not. parallel_io) then
            call s_write_serial_data_files(q_cons_vf, q_T_sf, q_prim_vf, t_step, bc_type, beta)
        else
            call s_write_parallel_data_files(q_cons_vf, t_step, bc_type, beta)
        end if

    end subroutine s_write_data_files

    !>  The purpose of this subroutine is to open a new or pre-
        !!          existing run-time information file and append to it the
        !!      basic header information relevant to current simulation.
        !!      In general, this requires generating a table header for
        !!      those stability criteria which will be written at every
        !!      time-step.
    impure subroutine s_open_run_time_information_file

        character(LEN=name_len), parameter :: file_name = 'run_time.inf' !<
            !! Name of the run-time information file

        character(LEN=path_len + name_len) :: file_path !<
            !! Relative path to a file in the case directory

        character(LEN=8) :: file_date !<
            !! Creation date of the run-time information file

        ! Opening the run-time information file
        file_path = trim(case_dir)//'/'//trim(file_name)

        open (3, FILE=trim(file_path), &
              FORM='formatted', &
              STATUS='replace')

        write (3, '(A)') 'Description: Stability information at '// &
            'each time-step of the simulation. This'
        write (3, '(13X,A)') 'data is composed of the inviscid '// &
            'Courant–Friedrichs–Lewy (ICFL)'
        write (3, '(13X,A)') 'number, the viscous CFL (VCFL) number, '// &
            'the capillary CFL (CCFL)'
        write (3, '(13X,A)') 'number and the cell Reynolds (Rc) '// &
            'number. Please note that only'
        write (3, '(13X,A)') 'those stability conditions pertinent '// &
            'to the physics included in'
        write (3, '(13X,A)') 'the current computation are displayed.'

        call date_and_time(DATE=file_date)

        write (3, '(A)') 'Date: '//file_date(5:6)//'/'// &
            file_date(7:8)//'/'// &
            file_date(3:4)

        write (3, '(A)') ''; write (3, '(A)') ''

        ! Generating table header for the stability criteria to be outputted
        write (3, '(13X,A9,13X,A10,13X,A10,13X,A10)', advance="no") &
            trim('Time-step'), trim('dt'), trim('Time'), trim('ICFL Max')

        if (viscous) then
            write (3, '(13X,A10,13X,A16)', advance="no") &
                trim('VCFL Max'), trim('Rc Min')
        end if

        write (3, *) ! new line

    end subroutine s_open_run_time_information_file

    !>  This opens a formatted data file where the root processor
        !!      can write out the CoM information
    impure subroutine s_open_com_files()

        character(len=path_len + 3*name_len) :: file_path !<
            !! Relative path to the CoM file in the case directory
        integer :: i !< Generic loop iterator

        do i = 1, num_fluids
            ! Generating the relative path to the CoM data file
            write (file_path, '(A,I0,A)') '/fluid', i, '_com.dat'
            file_path = trim(case_dir)//trim(file_path)
            ! Creating the formatted data file and setting up its
            ! structure
            open (i + 120, file=trim(file_path), &
                  form='formatted', &
                  position='append', &
                  status='unknown')
            if (n == 0) then
                write (i + 120, '(A)') '    Non-Dimensional Time '// &
                    '    Total Mass '// &
                    '    x-loc '// &
                    '    Total Volume    '
            elseif (p == 0) then
                write (i + 120, '(A)') '    Non-Dimensional Time '// &
                    '    Total Mass '// &
                    '    x-loc '// &
                    '    y-loc '// &
                    '    Total Volume    '
            else
                write (i + 120, '(A)') '    Non-Dimensional Time '// &
                    '    Total Mass '// &
                    '    x-loc '// &
                    '    y-loc '// &
                    '    z-loc '// &
                    '    Total Volume    '
            end if
        end do
    end subroutine s_open_com_files

    !>  This opens a formatted data file where the root processor
        !!      can write out flow probe information
    impure subroutine s_open_probe_files

        character(LEN=path_len + 3*name_len) :: file_path !<
            !! Relative path to the probe data file in the case directory

        integer :: i !< Generic loop iterator
        logical :: file_exist

        do i = 1, num_probes
            ! Generating the relative path to the data file
            write (file_path, '(A,I0,A)') '/D/probe', i, '_prim.dat'
            file_path = trim(case_dir)//trim(file_path)

            ! Creating the formatted data file and setting up its
            ! structure
            inquire (file=trim(file_path), exist=file_exist)

            if (file_exist) then
                open (i + 30, FILE=trim(file_path), &
                      FORM='formatted', &
                      STATUS='old', &
                      POSITION='append')
            else
                open (i + 30, FILE=trim(file_path), &
                      FORM='formatted', &
                      STATUS='unknown')
            end if
        end do

        if (integral_wrt) then
            do i = 1, num_integrals
                write (file_path, '(A,I0,A)') '/D/integral', i, '_prim.dat'
                file_path = trim(case_dir)//trim(file_path)

                open (i + 70, FILE=trim(file_path), &
                      FORM='formatted', &
                      POSITION='append', &
                      STATUS='unknown')
            end do
        end if

    end subroutine s_open_probe_files

    !>  The goal of the procedure is to output to the run-time
        !!      information file the stability criteria extrema in the
        !!      entire computational domain and at the given time-step.
        !!      Moreover, the subroutine is also in charge of tracking
        !!      these stability criteria extrema over all time-steps.
        !!  @param q_prim_vf Cell-average primitive variables
        !!  @param t_step Current time step
    impure subroutine s_write_run_time_information(q_prim_vf, t_step)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        integer, intent(in) :: t_step

        real(wp) :: rho        !< Cell-avg. density
        real(wp), dimension(num_vels) :: vel        !< Cell-avg. velocity
        real(wp) :: vel_sum    !< Cell-avg. velocity sum
        real(wp) :: pres       !< Cell-avg. pressure
        real(wp), dimension(num_fluids) :: alpha      !< Cell-avg. volume fraction
        real(wp) :: gamma      !< Cell-avg. sp. heat ratio
        real(wp) :: pi_inf     !< Cell-avg. liquid stiffness function
        real(wp) :: c          !< Cell-avg. sound speed
        real(wp) :: H          !< Cell-avg. enthalpy
        real(wp), dimension(2) :: Re         !< Cell-avg. Reynolds numbers
        integer :: j, k, l

        ! Computing Stability Criteria at Current Time-step
        $:GPU_PARALLEL_LOOP(collapse=3, private='[vel, alpha, Re]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    call s_compute_enthalpy(q_prim_vf, pres, rho, gamma, pi_inf, Re, H, alpha, vel, vel_sum, j, k, l)

                    call s_compute_speed_of_sound(pres, rho, gamma, pi_inf, H, alpha, vel_sum, 0._wp, c)

                    if (viscous) then
                        call s_compute_stability_from_dt(vel, c, rho, Re, j, k, l, icfl_sf, vcfl_sf, Rc_sf)
                    else
                        call s_compute_stability_from_dt(vel, c, rho, Re, j, k, l, icfl_sf)
                    end if

                end do
            end do
        end do

        ! end: Computing Stability Criteria at Current Time-step

        ! Determining local stability criteria extrema at current time-step

#ifdef _CRAYFTN
        $:GPU_UPDATE(host='[icfl_sf]')

        if (viscous) then
            $:GPU_UPDATE(host='[vcfl_sf,Rc_sf]')
        end if

        icfl_max_loc = maxval(icfl_sf)

        if (viscous) then
            vcfl_max_loc = maxval(vcfl_sf)
            Rc_min_loc = minval(Rc_sf)
        end if
#else
        #:call GPU_PARALLEL(copyout='[icfl_max_loc]', copyin='[icfl_sf]')
            icfl_max_loc = maxval(icfl_sf)
        #:endcall GPU_PARALLEL
        if (viscous) then
            #:call GPU_PARALLEL(copyout='[vcfl_max_loc, Rc_min_loc]', copyin='[vcfl_sf,Rc_sf]')
                vcfl_max_loc = maxval(vcfl_sf)
                Rc_min_loc = minval(Rc_sf)
            #:endcall GPU_PARALLEL
        end if
#endif

        ! Determining global stability criteria extrema at current time-step
        if (num_procs > 1) then
            call s_mpi_reduce_stability_criteria_extrema(icfl_max_loc, &
                                                         vcfl_max_loc, &
                                                         Rc_min_loc, &
                                                         icfl_max_glb, &
                                                         vcfl_max_glb, &
                                                         Rc_min_glb)
        else
            icfl_max_glb = icfl_max_loc
            if (viscous) vcfl_max_glb = vcfl_max_loc
            if (viscous) Rc_min_glb = Rc_min_loc
        end if

        ! Determining the stability criteria extrema over all the time-steps
        if (icfl_max_glb > icfl_max) icfl_max = icfl_max_glb

        if (viscous) then
            if (vcfl_max_glb > vcfl_max) vcfl_max = vcfl_max_glb
            if (Rc_min_glb < Rc_min) Rc_min = Rc_min_glb
        end if

        ! Outputting global stability criteria extrema at current time-step
        if (proc_rank == 0) then
            write (3, '(13X,I9,13X,F10.6,13X,F10.6,13X,F10.6)', advance="no") &
                t_step, dt, mytime, icfl_max_glb

            if (viscous) then
                write (3, '(13X,F10.6,13X,ES16.6)', advance="no") &
                    vcfl_max_glb, &
                    Rc_min_glb
            end if

            write (3, *) ! new line

            if (.not. f_approx_equal(icfl_max_glb, icfl_max_glb)) then
                call s_mpi_abort('ICFL is NaN. Exiting.')
            elseif (icfl_max_glb > 1._wp) then
                print *, 'icfl', icfl_max_glb
                call s_mpi_abort('ICFL is greater than 1.0. Exiting.')
            end if

            if (viscous) then
                if (.not. f_approx_equal(vcfl_max_glb, vcfl_max_glb)) then
                    call s_mpi_abort('VCFL is NaN. Exiting.')
                elseif (vcfl_max_glb > 1._wp) then
                    print *, 'vcfl', vcfl_max_glb
                    call s_mpi_abort('VCFL is greater than 1.0. Exiting.')
                end if
            end if
        end if

        call s_mpi_barrier()

    end subroutine s_write_run_time_information

    !>  The goal of this subroutine is to output the grid and
        !!      conservative variables data files for given time-step.
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param q_prim_vf Cell-average primitive variables
        !!  @param t_step Current time-step
    impure subroutine s_write_serial_data_files(q_cons_vf, q_T_sf, q_prim_vf, t_step, bc_type, beta)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), intent(inout) :: q_T_sf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        integer, intent(in) :: t_step
        type(scalar_field), intent(inout), optional :: beta
        type(integer_field), dimension(1:num_dims, -1:1), intent(in) :: bc_type

        character(LEN=path_len + 2*name_len) :: t_step_dir !<
            !! Relative path to the current time-step directory

        character(LEN=path_len + 3*name_len) :: file_path !<
            !! Relative path to the grid and conservative variables data files

        logical :: file_exist !<
            !! Logical used to check existence of current time-step directory

        character(LEN=15) :: FMT

        integer :: i, j, k, l, r

        real(wp) :: gamma, lit_gamma, pi_inf, qv !< Temporary EOS params

        ! Creating or overwriting the time-step root directory
        write (t_step_dir, '(A,I0,A,I0)') trim(case_dir)//'/p_all'

        ! Creating or overwriting the current time-step directory
        write (t_step_dir, '(a,i0,a,i0)') trim(case_dir)//'/p_all/p', &
            proc_rank, '/', t_step

        file_path = trim(t_step_dir)//'/.'
        call my_inquire(file_path, file_exist)
        if (file_exist) call s_delete_directory(trim(t_step_dir))
        call s_create_directory(trim(t_step_dir))

        ! Writing the grid data file in the x-direction
        file_path = trim(t_step_dir)//'/x_cb.dat'

        open (2, FILE=trim(file_path), &
              FORM='unformatted', &
              STATUS='new')
        write (2) x_cb(-1:m); close (2)

        ! Writing the grid data files in the y- and z-directions
        if (n > 0) then

            file_path = trim(t_step_dir)//'/y_cb.dat'

            open (2, FILE=trim(file_path), &
                  FORM='unformatted', &
                  STATUS='new')
            write (2) y_cb(-1:n); close (2)

            if (p > 0) then

                file_path = trim(t_step_dir)//'/z_cb.dat'

                open (2, FILE=trim(file_path), &
                      FORM='unformatted', &
                      STATUS='new')
                write (2) z_cb(-1:p); close (2)

            end if

        end if

        ! Writing the conservative variables data files
        do i = 1, sys_size
            write (file_path, '(A,I0,A)') trim(t_step_dir)//'/q_cons_vf', &
                i, '.dat'

            open (2, FILE=trim(file_path), &
                  FORM='unformatted', &
                  STATUS='new')

            write (2) q_cons_vf(i)%sf(0:m, 0:n, 0:p); close (2)
        end do

        if (qbmm .and. .not. polytropic) then
            do i = 1, nb
                do r = 1, nnode
                    write (file_path, '(A,I0,A)') trim(t_step_dir)//'/pb', &
                        sys_size + (i - 1)*nnode + r, '.dat'

                    open (2, FILE=trim(file_path), &
                          FORM='unformatted', &
                          STATUS='new')

                    write (2) pb_ts(1)%sf(0:m, 0:n, 0:p, r, i); close (2)
                end do
            end do

            do i = 1, nb
                do r = 1, nnode
                    write (file_path, '(A,I0,A)') trim(t_step_dir)//'/mv', &
                        sys_size + (i - 1)*nnode + r, '.dat'

                    open (2, FILE=trim(file_path), &
                          FORM='unformatted', &
                          STATUS='new')

                    write (2) mv_ts(1)%sf(0:m, 0:n, 0:p, r, i); close (2)
                end do
            end do
        end if

        ! Writing the IB markers
        if (ib) then
            write (file_path, '(A,I0,A)') trim(t_step_dir)//'/ib.dat'

            open (2, FILE=trim(file_path), &
                  FORM='unformatted', &
                  STATUS='new')

            write (2) ib_markers%sf; close (2)
        end if

        gamma = fluid_pp(1)%gamma
        lit_gamma = 1._wp/fluid_pp(1)%gamma + 1._wp
        pi_inf = fluid_pp(1)%pi_inf
        qv = fluid_pp(1)%qv

        if (precision == 1) then
            FMT = "(2F30.3)"
        else
            FMT = "(2F40.14)"
        end if

        ! writing an output directory
        write (t_step_dir, '(A,I0,A,I0)') trim(case_dir)//'/D'
        file_path = trim(t_step_dir)//'/.'

        inquire (FILE=trim(file_path), EXIST=file_exist)

        if (.not. file_exist) call s_create_directory(trim(t_step_dir))

        if (prim_vars_wrt .or. (n == 0 .and. p == 0)) then
            call s_convert_conservative_to_primitive_variables(q_cons_vf, q_T_sf, q_prim_vf, idwint)
            do i = 1, sys_size
                $:GPU_UPDATE(host='[q_prim_vf(i)%sf(:,:,:)]')
            end do
            ! q_prim_vf(bubxb) stores the value of nb needed in riemann solvers, so replace with true primitive value (=1._wp)
            if (qbmm) then
                q_prim_vf(bubxb)%sf = 1._wp
            end if
        end if

        !1D
        if (n == 0 .and. p == 0) then

            if (model_eqns == 2 .and. (.not. igr)) then
                do i = 1, sys_size
                    write (file_path, '(A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir)//'/prim.', i, '.', proc_rank, '.', t_step, '.dat'

                    open (2, FILE=trim(file_path))
                    do j = 0, m
                        ! todo: revisit change here
                        if (((i >= adv_idx%beg) .and. (i <= adv_idx%end))) then
                            write (2, FMT) x_cb(j), q_cons_vf(i)%sf(j, 0, 0)
                        else
                            write (2, FMT) x_cb(j), q_prim_vf(i)%sf(j, 0, 0)
                        end if
                    end do
                    close (2)
                end do
            end if

            do i = 1, sys_size
                write (file_path, '(A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir)//'/cons.', i, '.', proc_rank, '.', t_step, '.dat'

                open (2, FILE=trim(file_path))
                do j = 0, m
                    write (2, FMT) x_cb(j), q_cons_vf(i)%sf(j, 0, 0)
                end do
                close (2)
            end do

            if (qbmm .and. .not. polytropic) then
                do i = 1, nb
                    do r = 1, nnode
                        write (file_path, '(A,I0,A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir)//'/pres.', i, '.', r, '.', proc_rank, '.', t_step, '.dat'

                        open (2, FILE=trim(file_path))
                        do j = 0, m
                            write (2, FMT) x_cb(j), pb_ts(1)%sf(j, 0, 0, r, i)
                        end do
                        close (2)
                    end do
                end do
                do i = 1, nb
                    do r = 1, nnode
                        write (file_path, '(A,I0,A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir)//'/mv.', i, '.', r, '.', proc_rank, '.', t_step, '.dat'

                        open (2, FILE=trim(file_path))
                        do j = 0, m
                            write (2, FMT) x_cb(j), mv_ts(1)%sf(j, 0, 0, r, i)
                        end do
                        close (2)
                    end do
                end do
            end if
        end if

        if (precision == 1) then
            FMT = "(3F30.7)"
        else
            FMT = "(3F40.14)"
        end if

        ! 2D
        if ((n > 0) .and. (p == 0)) then
            do i = 1, sys_size
                write (file_path, '(A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir)//'/cons.', i, '.', proc_rank, '.', t_step, '.dat'
                open (2, FILE=trim(file_path))
                do j = 0, m
                    do k = 0, n
                        write (2, FMT) x_cb(j), y_cb(k), q_cons_vf(i)%sf(j, k, 0)
                    end do
                    write (2, *)
                end do
                close (2)
            end do

            if (present(beta)) then
                write (file_path, '(A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir)//'/beta.', i, '.', proc_rank, '.', t_step, '.dat'
                open (2, FILE=trim(file_path))
                do j = 0, m
                    do k = 0, n
                        write (2, FMT) x_cb(j), y_cb(k), beta%sf(j, k, 0)
                    end do
                    write (2, *)
                end do
                close (2)
            end if

            if (qbmm .and. .not. polytropic) then
                do i = 1, nb
                    do r = 1, nnode
                        write (file_path, '(A,I0,A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir)//'/pres.', i, '.', r, '.', proc_rank, '.', t_step, '.dat'

                        open (2, FILE=trim(file_path))
                        do j = 0, m
                            do k = 0, n
                                write (2, FMT) x_cb(j), y_cb(k), pb_ts(1)%sf(j, k, 0, r, i)
                            end do
                        end do
                        close (2)
                    end do
                end do
                do i = 1, nb
                    do r = 1, nnode
                        write (file_path, '(A,I0,A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir)//'/mv.', i, '.', r, '.', proc_rank, '.', t_step, '.dat'

                        open (2, FILE=trim(file_path))
                        do j = 0, m
                            do k = 0, n
                                write (2, FMT) x_cb(j), y_cb(k), mv_ts(1)%sf(j, k, 0, r, i)
                            end do
                        end do
                        close (2)
                    end do
                end do
            end if

            if (prim_vars_wrt .and. (.not. igr)) then
                do i = 1, sys_size
                    write (file_path, '(A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir)//'/prim.', i, '.', proc_rank, '.', t_step, '.dat'

                    open (2, FILE=trim(file_path))

                    do j = 0, m
                        do k = 0, n
                            if (((i >= cont_idx%beg) .and. (i <= cont_idx%end)) &
                                .or. &
                                ((i >= adv_idx%beg) .and. (i <= adv_idx%end)) &
                                ) then
                                write (2, FMT) x_cb(j), y_cb(k), q_cons_vf(i)%sf(j, k, 0)
                            else
                                write (2, FMT) x_cb(j), y_cb(k), q_prim_vf(i)%sf(j, k, 0)
                            end if
                        end do
                        write (2, *)
                    end do
                    close (2)
                end do
            end if
        end if

        if (precision == 1) then
            FMT = "(4F30.7)"
        else
            FMT = "(4F40.14)"
        end if

        ! 3D
        if (p > 0) then
            do i = 1, sys_size
                write (file_path, '(A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir)//'/cons.', i, '.', proc_rank, '.', t_step, '.dat'
                open (2, FILE=trim(file_path))
                do j = 0, m
                    do k = 0, n
                        do l = 0, p
                            write (2, FMT) x_cb(j), y_cb(k), z_cb(l), q_cons_vf(i)%sf(j, k, l)
                        end do
                        write (2, *)
                    end do
                    write (2, *)
                end do
                close (2)
            end do

            if (present(beta)) then
                write (file_path, '(A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir)//'/beta.', i, '.', proc_rank, '.', t_step, '.dat'
                open (2, FILE=trim(file_path))
                do j = 0, m
                    do k = 0, n
                        do l = 0, p
                            write (2, FMT) x_cb(j), y_cb(k), z_cb(l), beta%sf(j, k, l)
                        end do
                        write (2, *)
                    end do
                    write (2, *)
                end do
                close (2)
            end if

            if (qbmm .and. .not. polytropic) then
                do i = 1, nb
                    do r = 1, nnode
                        write (file_path, '(A,I0,A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir)//'/pres.', i, '.', r, '.', proc_rank, '.', t_step, '.dat'

                        open (2, FILE=trim(file_path))
                        do j = 0, m
                            do k = 0, n
                                do l = 0, p
                                    write (2, FMT) x_cb(j), y_cb(k), z_cb(l), pb_ts(1)%sf(j, k, l, r, i)
                                end do
                            end do
                        end do
                        close (2)
                    end do
                end do
                do i = 1, nb
                    do r = 1, nnode
                        write (file_path, '(A,I0,A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir)//'/mv.', i, '.', r, '.', proc_rank, '.', t_step, '.dat'

                        open (2, FILE=trim(file_path))
                        do j = 0, m
                            do k = 0, n
                                do l = 0, p
                                    write (2, FMT) x_cb(j), y_cb(k), z_cb(l), mv_ts(1)%sf(j, k, l, r, i)
                                end do
                            end do
                        end do
                        close (2)
                    end do
                end do
            end if

            if (prim_vars_wrt .and. (.not. igr)) then
                do i = 1, sys_size
                    write (file_path, '(A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir)//'/prim.', i, '.', proc_rank, '.', t_step, '.dat'

                    open (2, FILE=trim(file_path))

                    do j = 0, m
                        do k = 0, n
                            do l = 0, p
                                if (((i >= cont_idx%beg) .and. (i <= cont_idx%end)) &
                                    .or. &
                                    ((i >= adv_idx%beg) .and. (i <= adv_idx%end)) &
                                    .or. &
                                    ((i >= chemxb) .and. (i <= chemxe)) &
                                    ) then
                                    write (2, FMT) x_cb(j), y_cb(k), z_cb(l), q_cons_vf(i)%sf(j, k, l)
                                else
                                    write (2, FMT) x_cb(j), y_cb(k), z_cb(l), q_prim_vf(i)%sf(j, k, l)
                                end if
                            end do
                            write (2, *)
                        end do
                        write (2, *)
                    end do
                    close (2)
                end do
            end if
        end if

    end subroutine s_write_serial_data_files

    !>  The goal of this subroutine is to output the grid and
        !!      conservative variables data files for given time-step.
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param t_step Current time-step
        !!  @param beta Eulerian void fraction from lagrangian bubbles
    impure subroutine s_write_parallel_data_files(q_cons_vf, t_step, bc_type, beta)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        integer, intent(in) :: t_step
        type(scalar_field), intent(inout), optional :: beta
        type(integer_field), &
            dimension(1:num_dims, -1:1), &
            intent(in) :: bc_type

#ifdef MFC_MPI

        integer :: ifile, ierr, data_size
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(kind=MPI_OFFSET_kind) :: disp
        integer(kind=MPI_OFFSET_kind) :: m_MOK, n_MOK, p_MOK
        integer(kind=MPI_OFFSET_kind) :: WP_MOK, var_MOK, str_MOK
        integer(kind=MPI_OFFSET_kind) :: NVARS_MOK
        integer(kind=MPI_OFFSET_kind) :: MOK

        character(LEN=path_len + 2*name_len) :: file_loc
        logical :: file_exist, dir_check
        character(len=10) :: t_step_string

        integer :: i !< Generic loop iterator

        integer :: alt_sys !< Altered system size for the lagrangian subgrid bubble model

        ! Down sampling variables
        integer :: m_ds, n_ds, p_ds
        integer :: m_glb_ds, n_glb_ds, p_glb_ds
        integer :: m_glb_save, n_glb_save, p_glb_save ! Global save size

        if (down_sample) then
            call s_populate_variables_buffers(bc_type, q_cons_vf)
            call s_downsample_data(q_cons_vf, q_cons_temp, &
                                   m_ds, n_ds, p_ds, m_glb_ds, n_glb_ds, p_glb_ds)
        end if

        if (present(beta)) then
            alt_sys = sys_size + 1
        else
            alt_sys = sys_size
        end if

        if (file_per_process) then

            call s_int_to_str(t_step, t_step_string)

            ! Initialize MPI data I/O
            if (down_sample) then
                call s_initialize_mpi_data_ds(q_cons_temp)
            else
                if (ib) then
                    call s_initialize_mpi_data(q_cons_vf, ib_markers, levelset, levelset_norm)
                else
                    call s_initialize_mpi_data(q_cons_vf)
                end if
            end if

            if (proc_rank == 0) then
                file_loc = trim(case_dir)//'/restart_data/lustre_'//trim(t_step_string)
                call my_inquire(file_loc, dir_check)
                if (dir_check .neqv. .true.) then
                    call s_create_directory(trim(file_loc))
                end if
                call s_create_directory(trim(file_loc))
            end if
            call s_mpi_barrier()
            call DelayFileAccess(proc_rank)

            ! Initialize MPI data I/O
            call s_initialize_mpi_data(q_cons_vf)

            ! Open the file to write all flow variables
            write (file_loc, '(I0,A,i7.7,A)') t_step, '_', proc_rank, '.dat'
            file_loc = trim(case_dir)//'/restart_data/lustre_'//trim(t_step_string)//trim(mpiiofs)//trim(file_loc)
            inquire (FILE=trim(file_loc), EXIST=file_exist)
            if (file_exist .and. proc_rank == 0) then
                call MPI_FILE_DELETE(file_loc, mpi_info_int, ierr)
            end if
            call MPI_FILE_OPEN(MPI_COMM_SELF, file_loc, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
                               mpi_info_int, ifile, ierr)

            if (down_sample) then
                ! Size of local arrays
                data_size = (m_ds + 3)*(n_ds + 3)*(p_ds + 3)
                m_glb_save = m_glb_ds + 1
                n_glb_save = n_glb_ds + 1
                p_glb_save = p_glb_ds + 1
            else
                ! Size of local arrays
                data_size = (m + 1)*(n + 1)*(p + 1)
                m_glb_save = m_glb + 1
                n_glb_save = n_glb + 1
                p_glb_save = p_glb + 1
            end if

            ! Resize some integers so MPI can write even the biggest files
            m_MOK = int(m_glb_save + 1, MPI_OFFSET_KIND)
            n_MOK = int(n_glb_save + 1, MPI_OFFSET_KIND)
            p_MOK = int(p_glb_save + 1, MPI_OFFSET_KIND)
            WP_MOK = int(8._wp, MPI_OFFSET_KIND)
            MOK = int(1._wp, MPI_OFFSET_KIND)
            str_MOK = int(name_len, MPI_OFFSET_KIND)
            NVARS_MOK = int(sys_size, MPI_OFFSET_KIND)

            if (bubbles_euler) then
                ! Write the data for each variable
                do i = 1, sys_size
                    var_MOK = int(i, MPI_OFFSET_KIND)

                    call MPI_FILE_WRITE_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                            mpi_p, status, ierr)
                end do
                !Write pb and mv for non-polytropic qbmm
                if (qbmm .and. .not. polytropic) then
                    do i = sys_size + 1, sys_size + 2*nb*nnode
                        var_MOK = int(i, MPI_OFFSET_KIND)

                        call MPI_FILE_WRITE_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                                mpi_p, status, ierr)
                    end do
                end if
            else
                if (down_sample) then
                    do i = 1, sys_size !TODO: check if correct (sys_size
                        var_MOK = int(i, MPI_OFFSET_KIND)

                        call MPI_FILE_WRITE_ALL(ifile, q_cons_temp(i)%sf, data_size, &
                                                mpi_p, status, ierr)
                    end do
                else
                    do i = 1, sys_size !TODO: check if correct (sys_size
                        var_MOK = int(i, MPI_OFFSET_KIND)

                        call MPI_FILE_WRITE_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                                mpi_p, status, ierr)
                    end do
                end if
            end if

            call MPI_FILE_CLOSE(ifile, ierr)
        else
            ! Initialize MPI data I/O

            if (ib) then
                call s_initialize_mpi_data(q_cons_vf, ib_markers, levelset, levelset_norm)
            elseif (present(beta)) then
                call s_initialize_mpi_data(q_cons_vf, beta=beta)
            else
                call s_initialize_mpi_data(q_cons_vf)
            end if

            write (file_loc, '(I0,A)') t_step, '.dat'
            file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
            inquire (FILE=trim(file_loc), EXIST=file_exist)
            if (file_exist .and. proc_rank == 0) then
                call MPI_FILE_DELETE(file_loc, mpi_info_int, ierr)
            end if
            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
                               mpi_info_int, ifile, ierr)

            ! Size of local arrays
            data_size = (m + 1)*(n + 1)*(p + 1)

            ! Resize some integers so MPI can write even the biggest files
            m_MOK = int(m_glb + 1, MPI_OFFSET_KIND)
            n_MOK = int(n_glb + 1, MPI_OFFSET_KIND)
            p_MOK = int(p_glb + 1, MPI_OFFSET_KIND)
            WP_MOK = int(8._wp, MPI_OFFSET_KIND)
            MOK = int(1._wp, MPI_OFFSET_KIND)
            str_MOK = int(name_len, MPI_OFFSET_KIND)
            NVARS_MOK = int(alt_sys, MPI_OFFSET_KIND)

            if (bubbles_euler) then
                ! Write the data for each variable
                do i = 1, sys_size
                    var_MOK = int(i, MPI_OFFSET_KIND)

                    ! Initial displacement to skip at beginning of file
                    disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                    call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, MPI_IO_DATA%view(i), &
                                           'native', mpi_info_int, ierr)
                    call MPI_FILE_WRITE_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                            mpi_p, status, ierr)
                end do
                !Write pb and mv for non-polytropic qbmm
                if (qbmm .and. .not. polytropic) then
                    do i = sys_size + 1, sys_size + 2*nb*nnode
                        var_MOK = int(i, MPI_OFFSET_KIND)

                        ! Initial displacement to skip at beginning of file
                        disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                        call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, MPI_IO_DATA%view(i), &
                                               'native', mpi_info_int, ierr)
                        call MPI_FILE_WRITE_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                                mpi_p, status, ierr)
                    end do
                end if
            else
                do i = 1, sys_size !TODO: check if correct (sys_size
                    var_MOK = int(i, MPI_OFFSET_KIND)

                    ! Initial displacement to skip at beginning of file
                    disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                    call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, MPI_IO_DATA%view(i), &
                                           'native', mpi_info_int, ierr)
                    call MPI_FILE_WRITE_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                            mpi_p, status, ierr)
                end do
            end if

            ! Correction for the lagrangian subgrid bubble model
            if (present(beta)) then
                var_MOK = int(sys_size + 1, MPI_OFFSET_KIND)

                ! Initial displacement to skip at beginning of file
                disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, MPI_IO_DATA%view(sys_size + 1), &
                                       'native', mpi_info_int, ierr)
                call MPI_FILE_WRITE_ALL(ifile, MPI_IO_DATA%var(sys_size + 1)%sf, data_size, &
                                        mpi_p, status, ierr)
            end if

            call MPI_FILE_CLOSE(ifile, ierr)
        end if
#endif

    end subroutine s_write_parallel_data_files

    !>  This writes a formatted data file where the root processor
    !!      can write out the CoM information
    !!  @param t_step Current time-step
    !!  @param q_com Center of mass information
    !!  @param moments Higher moment information
    impure subroutine s_write_com_files(t_step, c_mass_in)

        integer, intent(in) :: t_step
        real(wp), dimension(num_fluids, 5), intent(in) :: c_mass_in
        integer :: i !< Generic loop iterator
        real(wp) :: nondim_time !< Non-dimensional time

        ! Non-dimensional time calculation
        if (t_step_old /= dflt_int) then
            nondim_time = real(t_step + t_step_old, wp)*dt
        else
            nondim_time = real(t_step, wp)*dt
        end if

        if (proc_rank == 0) then
            if (n == 0) then ! 1D simulation
                do i = 1, num_fluids ! Loop through fluids
                    write (i + 120, '(6X,4F24.12)') &
                        nondim_time, &
                        c_mass_in(i, 1), &
                        c_mass_in(i, 2), &
                        c_mass_in(i, 5)
                end do
            elseif (p == 0) then ! 2D simulation
                do i = 1, num_fluids ! Loop through fluids
                    write (i + 120, '(6X,5F24.12)') &
                        nondim_time, &
                        c_mass_in(i, 1), &
                        c_mass_in(i, 2), &
                        c_mass_in(i, 3), &
                        c_mass_in(i, 5)
                end do
            else ! 3D simulation
                do i = 1, num_fluids ! Loop through fluids
                    write (i + 120, '(6X,6F24.12)') &
                        nondim_time, &
                        c_mass_in(i, 1), &
                        c_mass_in(i, 2), &
                        c_mass_in(i, 3), &
                        c_mass_in(i, 4), &
                        c_mass_in(i, 5)
                end do
            end if
        end if

    end subroutine s_write_com_files

    !>  This writes a formatted data file for the flow probe information
        !!  @param t_step Current time-step
        !!  @param q_cons_vf Conservative variables
        !!  @param accel_mag Acceleration magnitude information
    impure subroutine s_write_probe_files(t_step, q_cons_vf, accel_mag)

        integer, intent(in) :: t_step
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        real(wp), dimension(0:m, 0:n, 0:p), intent(in) :: accel_mag

        real(wp), dimension(-1:m) :: distx
        real(wp), dimension(-1:n) :: disty
        real(wp), dimension(-1:p) :: distz

        ! The cell-averaged partial densities, density, velocity, pressure,
        ! volume fractions, specific heat ratio function, liquid stiffness
        ! function, and sound speed.
        real(wp) :: lit_gamma, nbub
        real(wp) :: rho
        real(wp), dimension(num_vels) :: vel
        real(wp) :: pres
        real(wp) :: ptilde
        real(wp) :: ptot
        real(wp) :: alf
        real(wp) :: alfgr
        real(wp), dimension(num_fluids) :: alpha
        real(wp) :: gamma
        real(wp) :: pi_inf
        real(wp) :: qv
        real(wp) :: c
        real(wp) :: M00, M10, M01, M20, M11, M02
        real(wp) :: varR, varV
        real(wp), dimension(Nb) :: nR, R, nRdot, Rdot
        real(wp) :: nR3
        real(wp) :: accel
        real(wp) :: int_pres
        real(wp) :: max_pres
        real(wp), dimension(2) :: Re
        real(wp), dimension(6) :: tau_e
        real(wp) :: G_local
        real(wp) :: dyn_p, T
        real(wp) :: damage_state

        integer :: i, j, k, l, s, d !< Generic loop iterator

        real(wp) :: nondim_time !< Non-dimensional time

        real(wp) :: tmp !<
            !! Temporary variable to store quantity for mpi_allreduce

        integer :: npts !< Number of included integral points
        real(wp) :: rad, thickness !< For integral quantities
        logical :: trigger !< For integral quantities

        real(wp) :: rhoYks(1:num_species)

        T = dflt_T_guess

        ! Non-dimensional time calculation
        if (time_stepper == 23) then
            nondim_time = mytime
        else
            if (t_step_old /= dflt_int) then
                nondim_time = real(t_step + t_step_old, wp)*dt
            else
                nondim_time = real(t_step, wp)*dt
            end if
        end if

        do i = 1, num_probes
            ! Zeroing out flow variables for all processors
            rho = 0._wp
            do s = 1, num_vels
                vel(s) = 0._wp
            end do
            pres = 0._wp
            gamma = 0._wp
            pi_inf = 0._wp
            qv = 0._wp
            c = 0._wp
            accel = 0._wp
            nR = 0._wp; R = 0._wp
            nRdot = 0._wp; Rdot = 0._wp
            nbub = 0._wp
            M00 = 0._wp
            M10 = 0._wp
            M01 = 0._wp
            M20 = 0._wp
            M11 = 0._wp
            M02 = 0._wp
            varR = 0._wp; varV = 0._wp
            alf = 0._wp
            do s = 1, (num_dims*(num_dims + 1))/2
                tau_e(s) = 0._wp
            end do
            damage_state = 0._wp

            ! Find probe location in terms of indices on a
            ! specific processor
            if (n == 0) then ! 1D simulation
                if ((probe(i)%x >= x_cb(-1)) .and. (probe(i)%x <= x_cb(m))) then
                    do s = -1, m
                        distx(s) = x_cb(s) - probe(i)%x
                        if (distx(s) < 0._wp) distx(s) = 1000._wp
                    end do
                    j = minloc(distx, 1)
                    if (j == 1) j = 2 ! Pick first point if probe is at edge
                    k = 0
                    l = 0

                    if (chemistry) then
                        do d = 1, num_species
                            rhoYks(d) = q_cons_vf(chemxb + d - 1)%sf(j - 2, k, l)
                        end do
                    end if

                    ! Computing/Sharing necessary state variables
                    if (elasticity) then
                        call s_convert_to_mixture_variables(q_cons_vf, j - 2, k, l, &
                                                            rho, gamma, pi_inf, qv, &
                                                            Re, G_local, fluid_pp(:)%G)
                    else
                        call s_convert_to_mixture_variables(q_cons_vf, j - 2, k, l, &
                                                            rho, gamma, pi_inf, qv)
                    end if
                    do s = 1, num_vels
                        vel(s) = q_cons_vf(cont_idx%end + s)%sf(j - 2, k, l)/rho
                    end do

                    dyn_p = 0.5_wp*rho*dot_product(vel, vel)

                    if (elasticity) then
                        if (cont_damage) then
                            damage_state = q_cons_vf(damage_idx)%sf(j - 2, k, l)
                            G_local = G_local*max((1._wp - damage_state), 0._wp)
                        end if

                        call s_compute_pressure( &
                            q_cons_vf(1)%sf(j - 2, k, l), &
                            q_cons_vf(alf_idx)%sf(j - 2, k, l), &
                            dyn_p, pi_inf, gamma, rho, qv, rhoYks(:), pres, T, &
                            q_cons_vf(stress_idx%beg)%sf(j - 2, k, l), &
                            q_cons_vf(mom_idx%beg)%sf(j - 2, k, l), G_local)
                    else
                        call s_compute_pressure( &
                            q_cons_vf(1)%sf(j - 2, k, l), &
                            q_cons_vf(alf_idx)%sf(j - 2, k, l), &
                            dyn_p, pi_inf, gamma, rho, qv, rhoYks(:), pres, T)
                    end if

                    if (model_eqns == 4) then
                        lit_gamma = 1._wp/fluid_pp(1)%gamma + 1._wp
                    else if (elasticity) then
                        tau_e(1) = q_cons_vf(stress_idx%end)%sf(j - 2, k, l)/rho
                    end if

                    if (bubbles_euler) then
                        alf = q_cons_vf(alf_idx)%sf(j - 2, k, l)
                        if (num_fluids == 3) then
                            alfgr = q_cons_vf(alf_idx - 1)%sf(j - 2, k, l)
                        end if
                        do s = 1, nb
                            nR(s) = q_cons_vf(bub_idx%rs(s))%sf(j - 2, k, l)
                            nRdot(s) = q_cons_vf(bub_idx%vs(s))%sf(j - 2, k, l)
                        end do

                        if (adv_n) then
                            nbub = q_cons_vf(n_idx)%sf(j - 2, k, l)
                        else
                            nR3 = 0._wp
                            do s = 1, nb
                                nR3 = nR3 + weight(s)*(nR(s)**3._wp)
                            end do

                            nbub = sqrt((4._wp*pi/3._wp)*nR3/alf)
                        end if
#ifdef DEBUG
                        print *, 'In probe, nbub: ', nbub
#endif
                        if (qbmm) then
                            M00 = q_cons_vf(bub_idx%moms(1, 1))%sf(j - 2, k, l)/nbub
                            M10 = q_cons_vf(bub_idx%moms(1, 2))%sf(j - 2, k, l)/nbub
                            M01 = q_cons_vf(bub_idx%moms(1, 3))%sf(j - 2, k, l)/nbub
                            M20 = q_cons_vf(bub_idx%moms(1, 4))%sf(j - 2, k, l)/nbub
                            M11 = q_cons_vf(bub_idx%moms(1, 5))%sf(j - 2, k, l)/nbub
                            M02 = q_cons_vf(bub_idx%moms(1, 6))%sf(j - 2, k, l)/nbub

                            M10 = M10/M00
                            M01 = M01/M00
                            M20 = M20/M00
                            M11 = M11/M00
                            M02 = M02/M00

                            varR = M20 - M10**2._wp
                            varV = M02 - M01**2._wp
                        end if
                        R(:) = nR(:)/nbub
                        Rdot(:) = nRdot(:)/nbub

                        ptilde = ptil(j - 2, k, l)
                        ptot = pres - ptilde
                    end if

                    ! Compute mixture sound Speed
                    call s_compute_speed_of_sound(pres, rho, gamma, pi_inf, &
                                                  ((gamma + 1._wp)*pres + pi_inf)/rho, alpha, 0._wp, 0._wp, c)

                    accel = accel_mag(j - 2, k, l)
                end if
            elseif (p == 0) then ! 2D simulation
                if (chemistry) then
                    do d = 1, num_species
                        rhoYks(d) = q_cons_vf(chemxb + d - 1)%sf(j - 2, k - 2, l)
                    end do
                end if

                if ((probe(i)%x >= x_cb(-1)) .and. (probe(i)%x <= x_cb(m))) then
                    if ((probe(i)%y >= y_cb(-1)) .and. (probe(i)%y <= y_cb(n))) then
                        do s = -1, m
                            distx(s) = x_cb(s) - probe(i)%x
                            if (distx(s) < 0._wp) distx(s) = 1000._wp
                        end do
                        do s = -1, n
                            disty(s) = y_cb(s) - probe(i)%y
                            if (disty(s) < 0._wp) disty(s) = 1000._wp
                        end do
                        j = minloc(distx, 1)
                        k = minloc(disty, 1)
                        if (j == 1) j = 2 ! Pick first point if probe is at edge
                        if (k == 1) k = 2 ! Pick first point if probe is at edge
                        l = 0

                        ! Computing/Sharing necessary state variables
                        call s_convert_to_mixture_variables(q_cons_vf, j - 2, k - 2, l, &
                                                            rho, gamma, pi_inf, qv, &
                                                            Re, G_local, fluid_pp(:)%G)
                        do s = 1, num_vels
                            vel(s) = q_cons_vf(cont_idx%end + s)%sf(j - 2, k - 2, l)/rho
                        end do

                        dyn_p = 0.5_wp*rho*dot_product(vel, vel)

                        if (elasticity) then
                            if (cont_damage) then
                                damage_state = q_cons_vf(damage_idx)%sf(j - 2, k - 2, l)
                                G_local = G_local*max((1._wp - damage_state), 0._wp)
                            end if

                            call s_compute_pressure( &
                                q_cons_vf(1)%sf(j - 2, k - 2, l), &
                                q_cons_vf(alf_idx)%sf(j - 2, k - 2, l), &
                                dyn_p, pi_inf, gamma, rho, qv, &
                                rhoYks, &
                                pres, &
                                T, &
                                q_cons_vf(stress_idx%beg)%sf(j - 2, k - 2, l), &
                                q_cons_vf(mom_idx%beg)%sf(j - 2, k - 2, l), G_local)
                        else
                            call s_compute_pressure(q_cons_vf(E_idx)%sf(j - 2, k - 2, l), &
                                                    q_cons_vf(alf_idx)%sf(j - 2, k - 2, l), &
                                                    dyn_p, pi_inf, gamma, rho, qv, &
                                                    rhoYks, pres, T)
                        end if

                        if (model_eqns == 4) then
                            lit_gamma = 1._wp/fluid_pp(1)%gamma + 1._wp
                        else if (elasticity) then
                            do s = 1, 3
                                tau_e(s) = q_cons_vf(s)%sf(j - 2, k - 2, l)/rho
                            end do
                        end if

                        if (bubbles_euler) then
                            alf = q_cons_vf(alf_idx)%sf(j - 2, k - 2, l)
                            do s = 1, nb
                                nR(s) = q_cons_vf(bub_idx%rs(s))%sf(j - 2, k - 2, l)
                                nRdot(s) = q_cons_vf(bub_idx%vs(s))%sf(j - 2, k - 2, l)
                            end do

                            if (adv_n) then
                                nbub = q_cons_vf(n_idx)%sf(j - 2, k - 2, l)
                            else
                                nR3 = 0._wp
                                do s = 1, nb
                                    nR3 = nR3 + weight(s)*(nR(s)**3._wp)
                                end do

                                nbub = sqrt((4._wp*pi/3._wp)*nR3/alf)
                            end if

                            R(:) = nR(:)/nbub
                            Rdot(:) = nRdot(:)/nbub
                        end if
                        ! Compute mixture sound speed
                        call s_compute_speed_of_sound(pres, rho, gamma, pi_inf, &
                                                      ((gamma + 1._wp)*pres + pi_inf)/rho, alpha, 0._wp, 0._wp, c)

                    end if
                end if
            else ! 3D
                if ((probe(i)%x >= x_cb(-1)) .and. (probe(i)%x <= x_cb(m))) then
                    if ((probe(i)%y >= y_cb(-1)) .and. (probe(i)%y <= y_cb(n))) then
                        if ((probe(i)%z >= z_cb(-1)) .and. (probe(i)%z <= z_cb(p))) then
                            do s = -1, m
                                distx(s) = x_cb(s) - probe(i)%x
                                if (distx(s) < 0._wp) distx(s) = 1000._wp
                            end do
                            do s = -1, n
                                disty(s) = y_cb(s) - probe(i)%y
                                if (disty(s) < 0._wp) disty(s) = 1000._wp
                            end do
                            do s = -1, p
                                distz(s) = z_cb(s) - probe(i)%z
                                if (distz(s) < 0._wp) distz(s) = 1000._wp
                            end do
                            j = minloc(distx, 1)
                            k = minloc(disty, 1)
                            l = minloc(distz, 1)
                            if (j == 1) j = 2 ! Pick first point if probe is at edge
                            if (k == 1) k = 2 ! Pick first point if probe is at edge
                            if (l == 1) l = 2 ! Pick first point if probe is at edge

                            ! Computing/Sharing necessary state variables
                            call s_convert_to_mixture_variables(q_cons_vf, j - 2, k - 2, l - 2, &
                                                                rho, gamma, pi_inf, qv, &
                                                                Re, G_local, fluid_pp(:)%G)
                            do s = 1, num_vels
                                vel(s) = q_cons_vf(cont_idx%end + s)%sf(j - 2, k - 2, l - 2)/rho
                            end do

                            dyn_p = 0.5_wp*rho*dot_product(vel, vel)

                            if (chemistry) then
                                do d = 1, num_species
                                    rhoYks(d) = q_cons_vf(chemxb + d - 1)%sf(j - 2, k - 2, l - 2)
                                end do
                            end if

                            if (elasticity) then
                                if (cont_damage) then
                                    damage_state = q_cons_vf(damage_idx)%sf(j - 2, k - 2, l - 2)
                                    G_local = G_local*max((1._wp - damage_state), 0._wp)
                                end if

                                call s_compute_pressure( &
                                    q_cons_vf(1)%sf(j - 2, k - 2, l - 2), &
                                    q_cons_vf(alf_idx)%sf(j - 2, k - 2, l - 2), &
                                    dyn_p, pi_inf, gamma, rho, qv, &
                                    rhoYks, pres, T, &
                                    q_cons_vf(stress_idx%beg)%sf(j - 2, k - 2, l - 2), &
                                    q_cons_vf(mom_idx%beg)%sf(j - 2, k - 2, l - 2), G_local)
                            else
                                call s_compute_pressure(q_cons_vf(E_idx)%sf(j - 2, k - 2, l - 2), &
                                                        q_cons_vf(alf_idx)%sf(j - 2, k - 2, l - 2), &
                                                        dyn_p, pi_inf, gamma, rho, qv, &
                                                        rhoYks, pres, T)
                            end if

                            ! Compute mixture sound speed
                            call s_compute_speed_of_sound(pres, rho, gamma, pi_inf, &
                                                          ((gamma + 1._wp)*pres + pi_inf)/rho, alpha, 0._wp, 0._wp, c)

                            accel = accel_mag(j - 2, k - 2, l - 2)
                        end if
                    end if
                end if
            end if
            if (num_procs > 1) then
                #:for VAR in ['rho','pres','gamma','pi_inf','qv','c','accel']
                    tmp = ${VAR}$
                    call s_mpi_allreduce_sum(tmp, ${VAR}$)
                #:endfor

                do s = 1, num_vels
                    tmp = vel(s)
                    call s_mpi_allreduce_sum(tmp, vel(s))
                end do

                if (bubbles_euler) then
                    #:for VAR in ['alf','alfgr','nbub','nR(1)','nRdot(1)','M00','R(1)','Rdot(1)','ptilde','ptot']
                        tmp = ${VAR}$
                        call s_mpi_allreduce_sum(tmp, ${VAR}$)
                    #:endfor

                    if (qbmm) then
                        #:for VAR in ['varR','varV','M10','M01','M20','M02']
                            tmp = ${VAR}$
                            call s_mpi_allreduce_sum(tmp, ${VAR}$)
                        #:endfor
                    end if
                end if

                if (elasticity) then
                    do s = 1, (num_dims*(num_dims + 1))/2
                        tmp = tau_e(s)
                        call s_mpi_allreduce_sum(tmp, tau_e(s))
                    end do
                end if

                if (cont_damage) then
                    tmp = damage_state
                    call s_mpi_allreduce_sum(tmp, damage_state)
                end if
            end if
            if (proc_rank == 0) then
                if (n == 0) then
                    if (bubbles_euler .and. (num_fluids <= 2)) then
                        if (qbmm) then
                            write (i + 30, '(6x,f12.6,14f28.16)') &
                                nondim_time, &
                                rho, &
                                vel(1), &
                                pres, &
                                alf, &
                                R(1), &
                                Rdot(1), &
                                nR(1), &
                                nRdot(1), &
                                varR, &
                                varV, &
                                M10, &
                                M01, &
                                M20, &
                                M02
                        else
                            write (i + 30, '(6x,f12.6,8f24.8)') &
                                nondim_time, &
                                rho, &
                                vel(1), &
                                pres, &
                                alf, &
                                R(1), &
                                Rdot(1), &
                                nR(1), &
                                nRdot(1)
                            ! ptilde, &
                            ! ptot
                        end if
                    else if (bubbles_euler .and. (num_fluids == 3)) then
                        write (i + 30, '(6x,f12.6,f24.8,f24.8,f24.8,f24.8,f24.8,'// &
                               'f24.8,f24.8,f24.8,f24.8,f24.8, f24.8)') &
                            nondim_time, &
                            rho, &
                            vel(1), &
                            pres, &
                            alf, &
                            alfgr, &
                            nR(1), &
                            nRdot(1), &
                            R(1), &
                            Rdot(1), &
                            ptilde, &
                            ptot
                    else if (bubbles_euler .and. num_fluids == 4) then
                        write (i + 30, '(6x,f12.6,f24.8,f24.8,f24.8,f24.8,'// &
                               'f24.8,f24.8,f24.8,f24.8,f24.8,f24.8,f24.8,f24.8,f24.8)') &
                            nondim_time, &
                            q_cons_vf(1)%sf(j - 2, 0, 0), &
                            q_cons_vf(2)%sf(j - 2, 0, 0), &
                            q_cons_vf(3)%sf(j - 2, 0, 0), &
                            q_cons_vf(4)%sf(j - 2, 0, 0), &
                            q_cons_vf(5)%sf(j - 2, 0, 0), &
                            q_cons_vf(6)%sf(j - 2, 0, 0), &
                            q_cons_vf(7)%sf(j - 2, 0, 0), &
                            q_cons_vf(8)%sf(j - 2, 0, 0), &
                            q_cons_vf(9)%sf(j - 2, 0, 0), &
                            q_cons_vf(10)%sf(j - 2, 0, 0), &
                            nbub, &
                            R(1), &
                            Rdot(1)
                    else
                        write (i + 30, '(6X,F12.6,F24.8,F24.8,F24.8)') &
                            nondim_time, &
                            rho, &
                            vel(1), &
                            pres
                    end if
                elseif (p == 0) then
                    if (bubbles_euler) then
                        write (i + 30, '(6X,10F24.8)') &
                            nondim_time, &
                            rho, &
                            vel(1), &
                            vel(2), &
                            pres, &
                            alf, &
                            nR(1), &
                            nRdot(1), &
                            R(1), &
                            Rdot(1)
                    else if (elasticity) then
                        write (i + 30, '(6X,F12.6,F24.8,F24.8,F24.8,F24.8,'// &
                               'F24.8,F24.8,F24.8)') &
                            nondim_time, &
                            rho, &
                            vel(1), &
                            vel(2), &
                            pres, &
                            tau_e(1), &
                            tau_e(2), &
                            tau_e(3)
                    else
                        write (i + 30, '(6X,F12.6,F24.8,F24.8,F24.8)') &
                            nondim_time, &
                            rho, &
                            vel(1), &
                            pres
                        print *, 'time =', nondim_time, 'rho =', rho, 'pres =', pres
                    end if
                else
                    write (i + 30, '(6X,F12.6,F24.8,F24.8,F24.8,F24.8,'// &
                           'F24.8,F24.8,F24.8,F24.8,F24.8,'// &
                           'F24.8)') &
                        nondim_time, &
                        rho, &
                        vel(1), &
                        vel(2), &
                        vel(3), &
                        pres, &
                        gamma, &
                        pi_inf, &
                        qv, &
                        c, &
                        accel
                end if
            end if
        end do

        if (integral_wrt .and. bubbles_euler) then
            if (n == 0) then ! 1D simulation
                do i = 1, num_integrals
                    int_pres = 0._wp
                    max_pres = 0._wp
                    k = 0; l = 0
                    npts = 0
                    do j = 1, m
                        pres = 0._wp
                        do s = 1, num_vels
                            vel(s) = 0._wp
                        end do
                        rho = 0._wp
                        pres = 0._wp
                        gamma = 0._wp
                        pi_inf = 0._wp
                        qv = 0._wp

                        if ((integral(i)%xmin <= x_cb(j)) .and. (integral(i)%xmax >= x_cb(j))) then
                            npts = npts + 1
                            call s_convert_to_mixture_variables(q_cons_vf, j, k, l, &
                                                                rho, gamma, pi_inf, qv, Re)
                            do s = 1, num_vels
                                vel(s) = q_cons_vf(cont_idx%end + s)%sf(j, k, l)/rho
                            end do

                            pres = ( &
                                   (q_cons_vf(E_idx)%sf(j, k, l) - &
                                    0.5_wp*(q_cons_vf(mom_idx%beg)%sf(j, k, l)**2._wp)/rho)/ &
                                   (1._wp - q_cons_vf(alf_idx)%sf(j, k, l)) - &
                                   pi_inf - qv &
                                   )/gamma
                            int_pres = int_pres + (pres - 1._wp)**2._wp
                        end if
                    end do
                    int_pres = sqrt(int_pres/(1._wp*npts))

                    if (num_procs > 1) then
                        tmp = int_pres
                        call s_mpi_allreduce_sum(tmp, int_pres)
                    end if

                    if (proc_rank == 0) then
                        if (bubbles_euler .and. (num_fluids <= 2)) then
                            write (i + 70, '(6x,f12.6,f24.8)') &
                                nondim_time, int_pres
                        end if
                    end if
                end do
            elseif (p == 0) then
                if (num_integrals /= 3) then
                    call s_mpi_abort('Incorrect number of integrals')
                end if

                rad = integral(1)%xmax
                thickness = integral(1)%xmin

                do i = 1, num_integrals
                    int_pres = 0._wp
                    max_pres = 0._wp
                    l = 0
                    npts = 0
                    do j = 1, m
                        do k = 1, n
                            trigger = .false.
                            if (i == 1) then
                                !inner portion
                                if (sqrt(x_cb(j)**2._wp + y_cb(k)**2._wp) < (rad - 0.5_wp*thickness)) &
                                    trigger = .true.
                            elseif (i == 2) then
                                !net region
                                if (sqrt(x_cb(j)**2._wp + y_cb(k)**2._wp) > (rad - 0.5_wp*thickness) .and. &
                                    sqrt(x_cb(j)**2._wp + y_cb(k)**2._wp) < (rad + 0.5_wp*thickness)) &
                                    trigger = .true.
                            elseif (i == 3) then
                                !everything else
                                if (sqrt(x_cb(j)**2._wp + y_cb(k)**2._wp) > (rad + 0.5_wp*thickness)) &
                                    trigger = .true.
                            end if

                            pres = 0._wp
                            do s = 1, num_vels
                                vel(s) = 0._wp
                            end do
                            rho = 0._wp
                            pres = 0._wp
                            gamma = 0._wp
                            pi_inf = 0._wp
                            qv = 0._wp

                            if (trigger) then
                                npts = npts + 1
                                call s_convert_to_mixture_variables(q_cons_vf, j, k, l, &
                                                                    rho, gamma, pi_inf, qv, Re)
                                do s = 1, num_vels
                                    vel(s) = q_cons_vf(cont_idx%end + s)%sf(j, k, l)/rho
                                end do

                                pres = ( &
                                       (q_cons_vf(E_idx)%sf(j, k, l) - &
                                        0.5_wp*(q_cons_vf(mom_idx%beg)%sf(j, k, l)**2._wp)/rho)/ &
                                       (1._wp - q_cons_vf(alf_idx)%sf(j, k, l)) - &
                                       pi_inf - qv &
                                       )/gamma
                                int_pres = int_pres + abs(pres - 1._wp)
                                max_pres = max(max_pres, abs(pres - 1._wp))
                            end if

                        end do
                    end do

                    if (npts > 0) then
                        int_pres = int_pres/(1._wp*npts)
                    else
                        int_pres = 0._wp
                    end if

                    if (num_procs > 1) then
                        tmp = int_pres
                        call s_mpi_allreduce_sum(tmp, int_pres)

                        tmp = max_pres
                        call s_mpi_allreduce_max(tmp, max_pres)
                    end if

                    if (proc_rank == 0) then
                        if (bubbles_euler .and. (num_fluids <= 2)) then
                            write (i + 70, '(6x,f12.6,f24.8,f24.8)') &
                                nondim_time, int_pres, max_pres
                        end if
                    end if
                end do
            end if
        end if

    end subroutine s_write_probe_files

    !>  The goal of this subroutine is to write to the run-time
        !!      information file basic footer information applicable to
        !!      the current computation and to close the file when done.
        !!      The footer contains the stability criteria extrema over
        !!      all of the time-steps and the simulation run-time.
    impure subroutine s_close_run_time_information_file

        real(wp) :: run_time !< Run-time of the simulation

        ! Writing the footer of and closing the run-time information file
        write (3, '(A)') '    '
        write (3, '(A)') ''

        write (3, '(A,F9.6)') 'ICFL Max: ', icfl_max
        if (viscous) write (3, '(A,F9.6)') 'VCFL Max: ', vcfl_max
        if (viscous) write (3, '(A,F10.6)') 'Rc Min: ', Rc_min

        call cpu_time(run_time)

        write (3, '(A)') ''
        write (3, '(A,I0,A)') 'Run-time: ', int(anint(run_time)), 's'
        write (3, '(A)') '    '
        close (3)

    end subroutine s_close_run_time_information_file

    !> Closes communication files
    impure subroutine s_close_com_files()

        integer :: i !< Generic loop iterator
        do i = 1, num_fluids
            close (i + 120)
        end do

    end subroutine s_close_com_files

    !> Closes probe files
    impure subroutine s_close_probe_files

        integer :: i !< Generic loop iterator

        do i = 1, num_probes
            close (i + 30)
        end do

    end subroutine s_close_probe_files

    !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    impure subroutine s_initialize_data_output_module

        integer :: i, m_ds, n_ds, p_ds

        ! Allocating/initializing ICFL, VCFL, CCFL and Rc stability criteria
        if (run_time_info) then
            @:ALLOCATE(icfl_sf(0:m, 0:n, 0:p))
            icfl_max = 0._wp

            if (viscous) then
                @:ALLOCATE(vcfl_sf(0:m, 0:n, 0:p))
                @:ALLOCATE(Rc_sf  (0:m, 0:n, 0:p))

                vcfl_max = 0._wp
                Rc_min = 1.e3_wp
            end if
        end if

        if (probe_wrt) then
            @:ALLOCATE(c_mass(num_fluids,5))
        end if

        if (down_sample) then
            m_ds = int((m + 1)/3) - 1
            n_ds = int((n + 1)/3) - 1
            p_ds = int((p + 1)/3) - 1

            allocate (q_cons_temp(1:sys_size))
            do i = 1, sys_size
                allocate (q_cons_temp(i)%sf(-1:m_ds + 1, -1:n_ds + 1, -1:p_ds + 1))
            end do
        end if

    end subroutine s_initialize_data_output_module

    !> Module deallocation and/or disassociation procedures
    impure subroutine s_finalize_data_output_module

        integer :: i

        if (probe_wrt) then
            @:DEALLOCATE(c_mass)
        end if

        if (run_time_info) then
            ! Deallocating the ICFL, VCFL, CCFL, and Rc stability criteria
            @:DEALLOCATE(icfl_sf)
            if (viscous) then
                @:DEALLOCATE(vcfl_sf, Rc_sf)
            end if
        end if

        if (down_sample) then
            do i = 1, sys_size
                deallocate (q_cons_temp(i)%sf)
            end do
            deallocate (q_cons_temp)
        end if

    end subroutine s_finalize_data_output_module

end module m_data_output
