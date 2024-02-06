!>
!! @file m_data_output.f90
!! @brief Contains module m_data_output

#:include 'macros.fpp'
#:include 'inline_conversions.fpp'

!> @brief The primary purpose of this module is to output the grid and the
!!              conservative variables data at the chosen time-step interval. In
!!              addition, this module is also in charge of outputting a run-time
!!              information file which summarizes the time-dependent behavior !of
!!              the stability criteria. The latter include the inviscid Courant–
!!              Friedrichs–Lewy (ICFL), viscous CFL (VCFL), capillary CFL (CCFL)
!!              and cell Reynolds (Rc) numbers.
module m_data_output

    !  Dependencies ============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    use m_compile_specific

    use m_helper

    use m_delay_file_access

    use m_ibm
    ! ==========================================================================

    implicit none

    private; public :: s_initialize_data_output_module, &
 s_open_run_time_information_file, &
 s_open_probe_files, &
 s_write_run_time_information, &
 s_write_data_files, &
 s_write_serial_data_files, &
 s_write_parallel_data_files, &
 s_write_probe_files, &
 s_close_run_time_information_file, &
 s_close_probe_files, &
 s_finalize_data_output_module

    abstract interface ! ===================================================

        !> Write data files
        !! @param q_cons_vf Conservative variables
        !! @param t_step Current time step
        subroutine s_write_abstract_data_files(q_cons_vf, q_prim_vf, t_step)

            import :: scalar_field, sys_size, pres_field

            type(scalar_field), &
                dimension(sys_size), &
                intent(IN) :: q_cons_vf

            type(scalar_field), &
                dimension(sys_size), &
                intent(INOUT) :: q_prim_vf

            integer, intent(IN) :: t_step

        end subroutine s_write_abstract_data_files ! -------------------
    end interface ! ========================================================

    real(kind(0d0)), allocatable, dimension(:, :, :) :: icfl_sf  !< ICFL stability criterion
    real(kind(0d0)), allocatable, dimension(:, :, :) :: vcfl_sf  !< VCFL stability criterion
    real(kind(0d0)), allocatable, dimension(:, :, :) :: ccfl_sf  !< CCFL stability criterion
    real(kind(0d0)), allocatable, dimension(:, :, :) :: Rc_sf  !< Rc stability criterion

    !$acc declare create(icfl_sf, vcfl_sf, ccfl_sf, Rc_sf)

    real(kind(0d0)) :: icfl_max_loc, icfl_max_glb !< ICFL stability extrema on local and global grids
    real(kind(0d0)) :: vcfl_max_loc, vcfl_max_glb !< VCFL stability extrema on local and global grids
    real(kind(0d0)) :: ccfl_max_loc, ccfl_max_glb !< CCFL stability extrema on local and global grids
    real(kind(0d0)) :: Rc_min_loc, Rc_min_glb !< Rc   stability extrema on local and global grids

    !$acc declare create(icfl_max_loc, icfl_max_glb, vcfl_max_loc, vcfl_max_glb, ccfl_max_loc, ccfl_max_glb, Rc_min_loc, Rc_min_glb)

    !> @name ICFL, VCFL, CCFL and Rc stability criteria extrema over all the time-steps
    !> @{
    real(kind(0d0)) :: icfl_max !< ICFL criterion maximum
    real(kind(0d0)) :: vcfl_max !< VCFL criterion maximum
    real(kind(0d0)) :: ccfl_max !< CCFL criterion maximum
    real(kind(0d0)) :: Rc_min !< Rc criterion maximum
    !> @}

    procedure(s_write_abstract_data_files), pointer :: s_write_data_files => null()

contains

    !>  The purpose of this subroutine is to open a new or pre-
        !!          existing run-time information file and append to it the
        !!      basic header information relevant to current simulation.
        !!      In general, this requires generating a table header for
        !!      those stability criteria which will be written at every
        !!      time-step.
    subroutine s_open_run_time_information_file() ! ------------------------

        character(LEN=name_len) :: file_name = 'run_time.inf' !<
            !! Name of the run-time information file

        character(LEN=path_len + name_len) :: file_path !<
            !! Relative path to a file in the case directory

        character(LEN=8) :: file_date !<
            !! Creation date of the run-time information file

        logical :: file_exist !<
            !! Logical used to check existence of run-time information file

        ! Opening the run-time information file
        file_path = trim(case_dir)//'/'//trim(file_name)

        inquire (FILE=trim(file_path), EXIST=file_exist)

        open (1, FILE=trim(file_path), &
              FORM='formatted', &
              POSITION='append', &
              STATUS='unknown')

        ! Generating file header for a new run-time information file
        if (file_exist .neqv. .true.) then

            write (1, '(A)') 'Description: Stability information at '// &
                'each time-step of the simulation. This'
            write (1, '(13X,A)') 'data is composed of the inviscid '// &
                'Courant–Friedrichs–Lewy (ICFL)'
            write (1, '(13X,A)') 'number, the viscous CFL (VCFL) number, '// &
                'the capillary CFL (CCFL)'
            write (1, '(13X,A)') 'number and the cell Reynolds (Rc) '// &
                'number. Please note that only'
            write (1, '(13X,A)') 'those stability conditions pertinent '// &
                'to the physics included in'
            write (1, '(13X,A)') 'the current computation are displayed.'

            call date_and_time(DATE=file_date)

            write (1, '(A)') 'Date: '//file_date(5:6)//'/'// &
                file_date(7:8)//'/'// &
                file_date(3:4)

        end if

        write (1, '(A)') ''; write (1, '(A)') ''

        ! Generating table header for the stability criteria to be outputted
        if (any(Re_size > 0)) then
            write (1, '(A)') '==== Time-steps ====== Time ======= ICFL '// &
                'Max ==== VCFL Max ====== Rc Min ======='
        else
            write (1, '(A)') '=========== Time-steps ============== Time '// &
                '============== ICFL Max ============='
        end if

    end subroutine s_open_run_time_information_file ! ----------------------

    !>  This opens a formatted data file where the root processor
        !!      can write out flow probe information
    subroutine s_open_probe_files() ! --------------------------------------

        character(LEN=path_len + 3*name_len) :: file_path !<
            !! Relative path to the probe data file in the case directory

        integer :: i !< Generic loop iterator

        do i = 1, num_probes
            ! Generating the relative path to the data file
            write (file_path, '(A,I0,A)') '/D/probe', i, '_prim.dat'
            file_path = trim(case_dir)//trim(file_path)

            ! Creating the formatted data file and setting up its
            ! structure
            open (i + 30, FILE=trim(file_path), &
                  FORM='formatted', &
                  STATUS='unknown')
            ! POSITION = 'append', &
            !WRITE(i+30,'(A,I0,A)') 'Probe ',i, ' located at:'
            !WRITE(i+30,'(A,F10.6)') 'x = ',probe(i)%x
            !WRITE(i+30,'(A,F10.6)') 'y = ',probe(i)%y
            !WRITE(i+30,'(A,F10.6)') 'z = ',probe(i)%z
            !WRITE(i+30, *)
            !WRITE(i+30,'(A)') '=== Non-Dimensional Time ' // &
            !                '=== Density ' // &
            !                '=== Velocity ' // &
            !                '=== Pressure ' // &
            !                '=== Gamma ' // &
            !                '=== Stiffness ' // &
            !                '=== Sound Speed ' // &
            !                '=== Acceleration ==='
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

    end subroutine s_open_probe_files ! ------------------------------------

    !>  The goal of the procedure is to output to the run-time
        !!      information file the stability criteria extrema in the
        !!      entire computational domain and at the given time-step.
        !!      Moreover, the subroutine is also in charge of tracking
        !!      these stability criteria extrema over all time-steps.
        !!  @param q_prim_vf Cell-average primitive variables
        !!  @param t_step Current time step
    subroutine s_write_run_time_information(q_prim_vf, t_step) ! -----------

        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf
        integer, intent(IN) :: t_step

        real(kind(0d0)), dimension(num_fluids) :: alpha_rho  !< Cell-avg. partial density
        real(kind(0d0)) :: rho        !< Cell-avg. density
        real(kind(0d0)), dimension(num_dims) :: vel        !< Cell-avg. velocity
        real(kind(0d0)) :: vel_sum    !< Cell-avg. velocity sum
        real(kind(0d0)) :: pres       !< Cell-avg. pressure
        real(kind(0d0)), dimension(num_fluids) :: alpha      !< Cell-avg. volume fraction
        real(kind(0d0)) :: gamma      !< Cell-avg. sp. heat ratio
        real(kind(0d0)) :: pi_inf     !< Cell-avg. liquid stiffness function
        real(kind(0d0)) :: qv         !< Cell-avg. fluid reference energy
        real(kind(0d0)) :: c          !< Cell-avg. sound speed
        real(kind(0d0)) :: E          !< Cell-avg. energy
        real(kind(0d0)) :: H          !< Cell-avg. enthalpy
        real(kind(0d0)), dimension(2) :: Re         !< Cell-avg. Reynolds numbers

        ! ICFL, VCFL, CCFL and Rc stability criteria extrema for the current
        ! time-step and located on both the local (loc) and the global (glb)
        ! computational domains

        real(kind(0d0)) :: blkmod1, blkmod2 !<
            !! Fluid bulk modulus for Woods mixture sound speed

        integer :: i, j, k, l, q !< Generic loop iterators

        integer :: Nfq
        real(kind(0d0)) :: fltr_dtheta   !<
            !! Modified dtheta accounting for Fourier filtering in azimuthal direction.

        ! Computing Stability Criteria at Current Time-step ================
        !$acc parallel loop collapse(3) gang vector default(present) private(alpha_rho, vel, alpha, Re)
        do l = 0, p
            do k = 0, n
                do j = 0, m

                    do i = 1, num_fluids
                        alpha_rho(i) = q_prim_vf(i)%sf(j, k, l)
                        alpha(i) = q_prim_vf(E_idx + i)%sf(j, k, l)
                    end do

                    if (bubbles) then
                        call s_convert_species_to_mixture_variables_bubbles_acc(rho, gamma, pi_inf, qv, alpha, alpha_rho, Re, j, k, l)
                    else
                        call s_convert_species_to_mixture_variables_acc(rho, gamma, pi_inf, qv, alpha, alpha_rho, Re, j, k, l)
                    end if

                    do i = 1, num_dims
                        vel(i) = q_prim_vf(contxe + i)%sf(j, k, l)
                    end do

                    vel_sum = 0d0
                    do i = 1, num_dims
                        vel_sum = vel_sum + vel(i)**2d0
                    end do

                    pres = q_prim_vf(E_idx)%sf(j, k, l)

                    E = gamma*pres + pi_inf + 5d-1*rho*vel_sum + qv

                    H = (E + pres)/rho

                    ! Compute mixture sound speed
                    call s_compute_speed_of_sound(pres, rho, gamma, pi_inf, H, alpha, vel_sum, c)

                    if (grid_geometry == 3) then
                        if (k == 0) then
                            fltr_dtheta = 2d0*pi*y_cb(0)/3d0
                        elseif (k <= fourier_rings) then
                            Nfq = min(floor(2d0*real(k, kind(0d0))*pi), (p + 1)/2 + 1)
                            fltr_dtheta = 2d0*pi*y_cb(k - 1)/real(Nfq, kind(0d0))
                        else
                            fltr_dtheta = y_cb(k - 1)*dz(l)
                        end if
                    end if

                    if (p > 0) then
                        !3D
                        if (grid_geometry == 3) then
                            icfl_sf(j, k, l) = dt/min(dx(j)/(abs(vel(1)) + c), &
                                                      dy(k)/(abs(vel(2)) + c), &
                                                      fltr_dtheta/(abs(vel(3)) + c))
                        else
                            icfl_sf(j, k, l) = dt/min(dx(j)/(abs(vel(1)) + c), &
                                                      dy(k)/(abs(vel(2)) + c), &
                                                      dz(l)/(abs(vel(3)) + c))
                        end if

                        if (any(Re_size > 0)) then

                            if (grid_geometry == 3) then
                                vcfl_sf(j, k, l) = maxval(dt/Re/rho) &
                                                   /min(dx(j), dy(k), fltr_dtheta)**2d0

                                Rc_sf(j, k, l) = min(dx(j)*(abs(vel(1)) + c), &
                                                     dy(k)*(abs(vel(2)) + c), &
                                                     fltr_dtheta*(abs(vel(3)) + c)) &
                                                 /maxval(1d0/Re)
                            else
                                vcfl_sf(j, k, l) = maxval(dt/Re/rho) &
                                                   /min(dx(j), dy(k), dz(l))**2d0

                                Rc_sf(j, k, l) = min(dx(j)*(abs(vel(1)) + c), &
                                                     dy(k)*(abs(vel(2)) + c), &
                                                     dz(l)*(abs(vel(3)) + c)) &
                                                 /maxval(1d0/Re)
                            end if

                        end if

                    elseif (n > 0) then
                        !2D
                        icfl_sf(j, k, l) = dt/min(dx(j)/(abs(vel(1)) + c), &
                                                  dy(k)/(abs(vel(2)) + c))

                        if (any(Re_size > 0)) then

                            vcfl_sf(j, k, l) = maxval(dt/Re/rho)/min(dx(j), dy(k))**2d0

                            Rc_sf(j, k, l) = min(dx(j)*(abs(vel(1)) + c), &
                                                 dy(k)*(abs(vel(2)) + c)) &
                                             /maxval(1d0/Re)

                        end if

                    else
                        !1D
                        icfl_sf(j, k, l) = (dt/dx(j))*(abs(vel(1)) + c)

                        if (any(Re_size > 0)) then

                            vcfl_sf(j, k, l) = maxval(dt/Re/rho)/dx(j)**2d0

                            Rc_sf(j, k, l) = dx(j)*(abs(vel(1)) + c)/maxval(1d0/Re)

                        end if

                    end if

                end do
            end do
        end do
        ! END: Computing Stability Criteria at Current Time-step ===========

        ! Determining local stability criteria extrema at current time-step

        !$acc kernels
        icfl_max_loc = maxval(icfl_sf)
        !$acc end kernels

        if (any(Re_size > 0)) then
            !$acc kernels
            vcfl_max_loc = maxval(vcfl_sf)
            Rc_min_loc = minval(Rc_sf)
            !$acc end kernels
        end if

        ! Determining global stability criteria extrema at current time-step
        if (num_procs > 1) then
            call s_mpi_reduce_stability_criteria_extrema(icfl_max_loc, &
                                                         vcfl_max_loc, &
                                                         ccfl_max_loc, &
                                                         Rc_min_loc, &
                                                         icfl_max_glb, &
                                                         vcfl_max_glb, &
                                                         ccfl_max_glb, &
                                                         Rc_min_glb)
        else
            icfl_max_glb = icfl_max_loc
            if (any(Re_size > 0)) vcfl_max_glb = vcfl_max_loc
            if (any(Re_size > 0)) Rc_min_glb = Rc_min_loc
        end if

        ! Determining the stability criteria extrema over all the time-steps
        if (icfl_max_glb > icfl_max) icfl_max = icfl_max_glb

        if (any(Re_size > 0)) then
            if (vcfl_max_glb > vcfl_max) vcfl_max = vcfl_max_glb
            if (Rc_min_glb < Rc_min) Rc_min = Rc_min_glb
        end if

        ! Outputting global stability criteria extrema at current time-step
        if (proc_rank == 0) then
            if (any(Re_size > 0)) then
                write (1, '(6X,I8,6X,F10.6,6X,F9.6,6X,F9.6,6X,F10.6)') &
                    t_step, t_step*dt, icfl_max_glb, &
                    vcfl_max_glb, &
                    Rc_min_glb
            else
                write (1, '(13X,I8,14X,F10.6,13X,F9.6)') &
                    t_step, t_step*dt, icfl_max_glb
            end if

            if (icfl_max_glb /= icfl_max_glb) then
                call s_mpi_abort('ICFL is NaN. Exiting ...')
            elseif (icfl_max_glb > 1d0) then
                print *, 'icfl', icfl_max_glb
                call s_mpi_abort('ICFL is greater than 1.0. Exiting ...')
            end if

            if (vcfl_max_glb /= vcfl_max_glb) then
                call s_mpi_abort('VCFL is NaN. Exiting ...')
            elseif (vcfl_max_glb > 1d0) then
                print *, 'vcfl', vcfl_max_glb
                call s_mpi_abort('VCFL is greater than 1.0. Exiting ...')
            end if
        end if

        call s_mpi_barrier()

    end subroutine s_write_run_time_information ! --------------------------

    !>  The goal of this subroutine is to output the grid and
        !!      conservative variables data files for given time-step.
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param t_step Current time-step
    subroutine s_write_serial_data_files(q_cons_vf, q_prim_vf, t_step) ! ---------------------

        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_prim_vf

        integer, intent(IN) :: t_step

        character(LEN=path_len + 2*name_len) :: t_step_dir !<
            !! Relative path to the current time-step directory

        character(LEN=path_len + 3*name_len) :: file_path !<
            !! Relative path to the grid and conservative variables data files

        logical :: file_exist !<
            !! Logical used to check existence of current time-step directory

        character(LEN=15) :: FMT

        integer :: i, j, k, l, ii, r!< Generic loop iterators

        real(kind(0d0)), dimension(nb) :: nRtmp         !< Temporary bubble concentration
        real(kind(0d0)) :: nbub, nR3, vftmp                         !< Temporary bubble number density
        real(kind(0d0)) :: gamma, lit_gamma, pi_inf, qv !< Temporary EOS params
        real(kind(0d0)) :: rho                          !< Temporary density
        real(kind(0d0)), dimension(2) :: Re !< Temporary Reynolds number
        real(kind(0d0)) :: E_e                          !< Temp. elastic energy contribution

        ! Creating or overwriting the time-step root directory
        write (t_step_dir, '(A,I0,A,I0)') trim(case_dir)//'/p_all'

        ! Creating or overwriting the current time-step directory
        write (t_step_dir, '(A,I0,A,I0)') trim(case_dir)//'/p_all/p', &
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
        lit_gamma = 1d0/fluid_pp(1)%gamma + 1d0
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
            call s_convert_conservative_to_primitive_variables(q_cons_vf, q_prim_vf)
            do i = 1, sys_size
                !$acc update host(q_prim_vf(i)%sf(:,:,:))
            end do
            ! q_prim_vf(bubxb) stores the value of nb needed in riemann solvers, so replace with true primitive value (=1d0)
            if (qbmm) then
                q_prim_vf(bubxb)%sf = 1d0
            end if
        end if

        !1D
        if (n == 0 .and. p == 0) then

            if (model_eqns == 2) then
                do i = 1, sys_size
                    write (file_path, '(A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir)//'/prim.', i, '.', proc_rank, '.', t_step, '.dat'

                    open (2, FILE=trim(file_path))
                    do j = 0, m
                        if (((i >= cont_idx%beg) .and. (i <= cont_idx%end)) &
                            .or. &
                            ((i >= adv_idx%beg) .and. (i <= adv_idx%end)) &
                            ) then
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

            if (prim_vars_wrt) then
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

            if (prim_vars_wrt) then
                do i = 1, sys_size
                    write (file_path, '(A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir)//'/prim.', i, '.', proc_rank, '.', t_step, '.dat'

                    open (2, FILE=trim(file_path))

                    do j = 0, m
                        do k = 0, n
                            do l = 0, p
                                if (((i >= cont_idx%beg) .and. (i <= cont_idx%end)) &
                                    .or. &
                                    ((i >= adv_idx%beg) .and. (i <= adv_idx%end)) &
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

    end subroutine s_write_serial_data_files ! ------------------------------------

    !>  The goal of this subroutine is to output the grid and
        !!      conservative variables data files for given time-step.
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param t_step Current time-step
    subroutine s_write_parallel_data_files(q_cons_vf, q_prim_vf, t_step) ! --

        type(scalar_field), &
            dimension(sys_size), &
            intent(IN) :: q_cons_vf

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: q_prim_vf

        integer, intent(IN) :: t_step

#ifdef MFC_MPI

        integer :: ifile, ierr, data_size
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(KIND=MPI_OFFSET_KIND) :: disp
        integer(KIND=MPI_OFFSET_KIND) :: m_MOK, n_MOK, p_MOK
        integer(KIND=MPI_OFFSET_KIND) :: WP_MOK, var_MOK, str_MOK
        integer(KIND=MPI_OFFSET_KIND) :: NVARS_MOK
        integer(KIND=MPI_OFFSET_KIND) :: MOK

        character(LEN=path_len + 2*name_len) :: file_loc
        logical :: file_exist, dir_check
        character(len=10) :: t_step_string

        integer :: i !< Generic loop iterator

        if (file_per_process) then

            call s_int_to_str(t_step, t_step_string)

            ! Initialize MPI data I/O

            if (ib) then
                call s_initialize_mpi_data(q_cons_vf, ib_markers)
            else
                call s_initialize_mpi_data(q_cons_vf)
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

            ! Size of local arrays
            data_size = (m + 1)*(n + 1)*(p + 1)

            ! Resize some integers so MPI can write even the biggest files
            m_MOK = int(m_glb + 1, MPI_OFFSET_KIND)
            n_MOK = int(n_glb + 1, MPI_OFFSET_KIND)
            p_MOK = int(p_glb + 1, MPI_OFFSET_KIND)
            WP_MOK = int(8d0, MPI_OFFSET_KIND)
            MOK = int(1d0, MPI_OFFSET_KIND)
            str_MOK = int(name_len, MPI_OFFSET_KIND)
            NVARS_MOK = int(sys_size, MPI_OFFSET_KIND)

            if (bubbles) then
                ! Write the data for each variable
                do i = 1, sys_size
                    var_MOK = int(i, MPI_OFFSET_KIND)

                    call MPI_FILE_WRITE_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                            MPI_DOUBLE_PRECISION, status, ierr)
                end do
                !Write pb and mv for non-polytropic qbmm
                if (qbmm .and. .not. polytropic) then
                    do i = sys_size + 1, sys_size + 2*nb*nnode
                        var_MOK = int(i, MPI_OFFSET_KIND)

                        call MPI_FILE_WRITE_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                                MPI_DOUBLE_PRECISION, status, ierr)
                    end do
                end if
            else
                do i = 1, sys_size !TODO: check if correct (sys_size
                    var_MOK = int(i, MPI_OFFSET_KIND)

                    call MPI_FILE_WRITE_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                            MPI_DOUBLE_PRECISION, status, ierr)
                end do
            end if

            call MPI_FILE_CLOSE(ifile, ierr)
        else
            ! Initialize MPI data I/O

            call s_initialize_mpi_data(q_cons_vf)

            ! Open the file to write all flow variables
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
            WP_MOK = int(8d0, MPI_OFFSET_KIND)
            MOK = int(1d0, MPI_OFFSET_KIND)
            str_MOK = int(name_len, MPI_OFFSET_KIND)
            NVARS_MOK = int(sys_size, MPI_OFFSET_KIND)

            if (bubbles) then
                ! Write the data for each variable
                do i = 1, sys_size
                    var_MOK = int(i, MPI_OFFSET_KIND)

                    ! Initial displacement to skip at beginning of file
                    disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                    call MPI_FILE_SET_VIEW(ifile, disp, MPI_DOUBLE_PRECISION, MPI_IO_DATA%view(i), &
                                           'native', mpi_info_int, ierr)
                    call MPI_FILE_WRITE_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                            MPI_DOUBLE_PRECISION, status, ierr)
                end do
                !Write pb and mv for non-polytropic qbmm
                if (qbmm .and. .not. polytropic) then
                    do i = sys_size + 1, sys_size + 2*nb*nnode
                        var_MOK = int(i, MPI_OFFSET_KIND)

                        ! Initial displacement to skip at beginning of file
                        disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                        call MPI_FILE_SET_VIEW(ifile, disp, MPI_DOUBLE_PRECISION, MPI_IO_DATA%view(i), &
                                               'native', mpi_info_int, ierr)
                        call MPI_FILE_WRITE_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                                MPI_DOUBLE_PRECISION, status, ierr)
                    end do
                end if
            else
                do i = 1, sys_size !TODO: check if correct (sys_size
                    var_MOK = int(i, MPI_OFFSET_KIND)

                    ! Initial displacement to skip at beginning of file
                    disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                    call MPI_FILE_SET_VIEW(ifile, disp, MPI_DOUBLE_PRECISION, MPI_IO_DATA%view(i), &
                                           'native', mpi_info_int, ierr)
                    call MPI_FILE_WRITE_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                            MPI_DOUBLE_PRECISION, status, ierr)
                end do
            end if

            call MPI_FILE_CLOSE(ifile, ierr)
        end if

        if (ib) then
            var_MOK = int(sys_size + 1, MPI_OFFSET_KIND)

            ! Initial displacement to skip at beginning of file
            disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

            call MPI_FILE_SET_VIEW(ifile, disp, MPI_INTEGER, MPI_IO_IB_DATA%view, &
                                   'native', mpi_info_int, ierr)
            call MPI_FILE_WRITE_ALL(ifile, MPI_IO_IB_DATA%var%sf, data_size, &
                                    MPI_DOUBLE_PRECISION, status, ierr)
        end if

        call MPI_FILE_CLOSE(ifile, ierr)

#endif

    end subroutine s_write_parallel_data_files ! ---------------------------

    !>  This writes a formatted data file for the flow probe information
        !!  @param t_step Current time-step
        !!  @param q_cons_vf Conservative variables
        !!  @param accel_mag Acceleration magnitude information
    subroutine s_write_probe_files(t_step, q_cons_vf, accel_mag) ! -----------

        integer, intent(IN) :: t_step
        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        real(kind(0d0)), dimension(0:m, 0:n, 0:p), intent(IN) :: accel_mag

        real(kind(0d0)), dimension(-1:m) :: distx
        real(kind(0d0)), dimension(-1:n) :: disty
        real(kind(0d0)), dimension(-1:p) :: distz

        ! The cell-averaged partial densities, density, velocity, pressure,
        ! volume fractions, specific heat ratio function, liquid stiffness
        ! function, and sound speed.
        real(kind(0d0)) :: lit_gamma, nbub
        real(kind(0d0)) :: rho
        real(kind(0d0)), dimension(num_dims) :: vel
        real(kind(0d0)) :: pres
        real(kind(0d0)) :: ptilde
        real(kind(0d0)) :: ptot
        real(kind(0d0)) :: alf
        real(kind(0d0)) :: alfgr
        real(kind(0d0)), dimension(num_fluids) :: alpha
        real(kind(0d0)) :: gamma
        real(kind(0d0)) :: pi_inf
        real(kind(0d0)) :: qv
        real(kind(0d0)) :: c
        real(kind(0d0)) :: M00, M10, M01, M20, M11, M02
        real(kind(0d0)) :: varR, varV
        real(kind(0d0)), dimension(Nb) :: nR, R, nRdot, Rdot
        real(kind(0d0)) :: nR3
        real(kind(0d0)) :: accel
        real(kind(0d0)) :: int_pres
        real(kind(0d0)) :: max_pres
        real(kind(0d0)), dimension(2) :: Re
        real(kind(0d0)) :: E_e
        real(kind(0d0)), dimension(6) :: tau_e
        real(kind(0d0)) :: G

        integer :: i, j, k, l, s, q !< Generic loop iterator

        real(kind(0d0)) :: nondim_time !< Non-dimensional time

        real(kind(0d0)) :: tmp !<
            !! Temporary variable to store quantity for mpi_allreduce

        real(kind(0d0)) :: blkmod1, blkmod2 !<
            !! Fluid bulk modulus for Woods mixture sound speed

        integer :: npts !< Number of included integral points
        real(kind(0d0)) :: rad, thickness !< For integral quantities
        logical :: trigger !< For integral quantities

        ! Non-dimensional time calculation
        if (time_stepper == 23) then
            nondim_time = mytime
        else
            if (t_step_old /= dflt_int) then
                nondim_time = real(t_step + t_step_old, kind(0d0))*dt
            else
                nondim_time = real(t_step, kind(0d0))*dt !*1.d-5/10.0761131451d0
            end if
        end if

        do i = 1, num_probes
            ! Zeroing out flow variables for all processors
            rho = 0d0
            do s = 1, num_dims
                vel(s) = 0d0
            end do
            pres = 0d0
            gamma = 0d0
            pi_inf = 0d0
            qv = 0d0
            c = 0d0
            accel = 0d0
            nR = 0d0; R = 0d0
            nRdot = 0d0; Rdot = 0d0
            nbub = 0d0
            M00 = 0d0
            M10 = 0d0
            M01 = 0d0
            M20 = 0d0
            M11 = 0d0
            M02 = 0d0
            varR = 0d0; varV = 0d0
            alf = 0d0
            do s = 1, (num_dims*(num_dims + 1))/2
                tau_e(s) = 0d0
            end do

            ! Find probe location in terms of indices on a
            ! specific processor
            if (n == 0) then ! 1D simulation
                if ((probe(i)%x >= x_cb(-1)) .and. (probe(i)%x <= x_cb(m))) then
                    do s = -1, m
                        distx(s) = x_cb(s) - probe(i)%x
                        if (distx(s) < 0d0) distx(s) = 1000d0
                    end do
                    j = minloc(distx, 1)
                    if (j == 1) j = 2 ! Pick first point if probe is at edge
                    k = 0
                    l = 0

                    ! Computing/Sharing necessary state variables
                    if (hypoelasticity) then
                        call s_convert_to_mixture_variables(q_cons_vf, j - 2, k, l, &
                                                            rho, gamma, pi_inf, qv, &
                                                            Re, G, fluid_pp(:)%G)
                    else
                        call s_convert_to_mixture_variables(q_cons_vf, j - 2, k, l, &
                                                            rho, gamma, pi_inf, qv)
                    end if
                    do s = 1, num_dims
                        vel(s) = q_cons_vf(cont_idx%end + s)%sf(j - 2, k, l)/rho
                    end do

                    if (hypoelasticity) then
                        call s_compute_pressure( &
                            q_cons_vf(1)%sf(j - 2, k, l), &
                            q_cons_vf(alf_idx)%sf(j - 2, k, l), &
                            0.5d0*(q_cons_vf(2)%sf(j - 2, k, l)**2.d0)/ &
                            q_cons_vf(1)%sf(j - 2, k, l), &
                            pi_inf, gamma, rho, qv, pres, &
                            q_cons_vf(stress_idx%beg)%sf(j - 2, k, l), &
                            q_cons_vf(mom_idx%beg)%sf(j - 2, k, l), G)
                    else
                        call s_compute_pressure( &
                            q_cons_vf(1)%sf(j - 2, k, l), &
                            q_cons_vf(alf_idx)%sf(j - 2, k, l), &
                            0.5d0*(q_cons_vf(2)%sf(j - 2, k, l)**2.d0)/ &
                            q_cons_vf(1)%sf(j - 2, k, l), &
                            pi_inf, gamma, rho, qv, pres)
                    end if

                    if (model_eqns == 4) then
                        lit_gamma = 1d0/fluid_pp(1)%gamma + 1d0
                    else if (hypoelasticity) then
                        tau_e(1) = q_cons_vf(stress_idx%end)%sf(j - 2, k, l)/rho
                    end if

                    if (bubbles) then
                        alf = q_cons_vf(alf_idx)%sf(j - 2, k, l)
                        if (num_fluids == 3) then
                            alfgr = q_cons_vf(alf_idx - 1)%sf(j - 2, k, l)
                        end if
                        do s = 1, nb
                            nR(s) = q_cons_vf(bub_idx%rs(s))%sf(j - 2, k, l)
                            nRdot(s) = q_cons_vf(bub_idx%vs(s))%sf(j - 2, k, l)
                        end do
                        !call comp_n_from_cons(alf, nR, nbub)

                        nR3 = 0d0
                        do s = 1, nb
                            nR3 = nR3 + weight(s)*(nR(s)**3d0)
                        end do

                        nbub = DSQRT((4.d0*pi/3.d0)*nR3/alf)

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

                            varR = M20 - M10**2d0
                            varV = M02 - M01**2d0
                        end if
                        R(:) = nR(:)/nbub
                        Rdot(:) = nRdot(:)/nbub

                        ptilde = ptil(j - 2, k, l)
                        ptot = pres - ptilde
                    end if

                    ! Compute mixture sound Speed
                    call s_compute_speed_of_sound(pres, rho, gamma, pi_inf, &
                                                  ((gamma + 1d0)*pres + pi_inf)/rho, alpha, 0d0, c)

                    accel = accel_mag(j - 2, k, l)
                end if
            elseif (p == 0) then ! 2D simulation
                if ((probe(i)%x >= x_cb(-1)) .and. (probe(i)%x <= x_cb(m))) then
                    if ((probe(i)%y >= y_cb(-1)) .and. (probe(i)%y <= y_cb(n))) then
                        do s = -1, m
                            distx(s) = x_cb(s) - probe(i)%x
                            if (distx(s) < 0d0) distx(s) = 1000d0
                        end do
                        do s = -1, n
                            disty(s) = y_cb(s) - probe(i)%y
                            if (disty(s) < 0d0) disty(s) = 1000d0
                        end do
                        j = minloc(distx, 1)
                        k = minloc(disty, 1)
                        if (j == 1) j = 2 ! Pick first point if probe is at edge
                        if (k == 1) k = 2 ! Pick first point if probe is at edge
                        l = 0

                        ! Computing/Sharing necessary state variables
                        call s_convert_to_mixture_variables(q_cons_vf, j - 2, k - 2, l, &
                                                            rho, gamma, pi_inf, qv, &
                                                            Re, G, fluid_pp(:)%G)
                        do s = 1, num_dims
                            vel(s) = q_cons_vf(cont_idx%end + s)%sf(j - 2, k - 2, l)/rho
                        end do

                        call s_compute_pressure( &
                            q_cons_vf(1)%sf(j - 2, k - 2, l), &
                            q_cons_vf(alf_idx)%sf(j - 2, k - 2, l), &
                            0.5d0*(q_cons_vf(2)%sf(j - 2, k - 2, l)**2.d0)/ &
                            q_cons_vf(1)%sf(j - 2, k - 2, l), &
                            pi_inf, gamma, rho, qv, pres, &
                            q_cons_vf(stress_idx%beg)%sf(j - 2, k - 2, l), &
                            q_cons_vf(mom_idx%beg)%sf(j - 2, k - 2, l), G)

                        if (model_eqns == 4) then
                            lit_gamma = 1d0/fluid_pp(1)%gamma + 1d0
                        else if (hypoelasticity) then
                            do s = 1, 3
                                tau_e(s) = q_cons_vf(s)%sf(j - 2, k - 2, l)/rho
                            end do
                        end if

                        if (bubbles) then
                            alf = q_cons_vf(alf_idx)%sf(j - 2, k - 2, l)
                            do s = 1, nb
                                nR(s) = q_cons_vf(bub_idx%rs(s))%sf(j - 2, k - 2, l)
                                nRdot(s) = q_cons_vf(bub_idx%vs(s))%sf(j - 2, k - 2, l)
                            end do
                            !call comp_n_from_cons(alf, nR, nbub)

                            nR3 = 0d0
                            do s = 1, nb
                                nR3 = nR3 + weight(s)*(nR(s)**3d0)
                            end do

                            nbub = DSQRT((4.d0*pi/3.d0)*nR3/alf)

                            R(:) = nR(:)/nbub
                            Rdot(:) = nRdot(:)/nbub
                        end if

                        ! Compute mixture sound speed
                        call s_compute_speed_of_sound(pres, rho, gamma, pi_inf, &
                                                      ((gamma + 1d0)*pres + pi_inf)/rho, alpha, 0d0, c)

                        accel = accel_mag(j - 2, k - 2, l)
                    end if
                end if
            else ! 3D simulation
                if ((probe(i)%x >= x_cb(-1)) .and. (probe(i)%x <= x_cb(m))) then
                    if ((probe(i)%y >= y_cb(-1)) .and. (probe(i)%y <= y_cb(n))) then
                        if ((probe(i)%z >= z_cb(-1)) .and. (probe(i)%z <= z_cb(p))) then
                            do s = -1, m
                                distx(s) = x_cb(s) - probe(i)%x
                                if (distx(s) < 0d0) distx(s) = 1000d0
                            end do
                            do s = -1, n
                                disty(s) = y_cb(s) - probe(i)%y
                                if (disty(s) < 0d0) disty(s) = 1000d0
                            end do
                            do s = -1, p
                                distz(s) = z_cb(s) - probe(i)%z
                                if (distz(s) < 0d0) distz(s) = 1000d0
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
                                                                Re, G, fluid_pp(:)%G)
                            do s = 1, num_dims
                                vel(s) = q_cons_vf(cont_idx%end + s)%sf(j - 2, k - 2, l - 2)/rho
                            end do

                            call s_compute_pressure(q_cons_vf(E_idx)%sf(j - 2, k - 2, l - 2), &
                                                    0d0, 0.5d0*rho*dot_product(vel, vel), pi_inf, gamma, rho, qv, pres)

                            ! Compute mixture sound speed
                            call s_compute_speed_of_sound(pres, rho, gamma, pi_inf, &
                                                          ((gamma + 1d0)*pres + pi_inf)/rho, alpha, 0d0, c)

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

                do s = 1, num_dims
                    tmp = vel(s)
                    call s_mpi_allreduce_sum(tmp, vel(s))
                end do

                if (bubbles) then
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

                if (hypoelasticity) then
                    do s = 1, (num_dims*(num_dims + 1))/2
                        tmp = tau_e(s)
                        call s_mpi_allreduce_sum(tmp, tau_e(s))
                    end do
                end if
            end if

            if (proc_rank == 0) then
                if (n == 0) then
                    if (bubbles .and. (num_fluids <= 2)) then
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
                    else if (bubbles .and. (num_fluids == 3)) then
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
                    else if (bubbles .and. num_fluids == 4) then
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
                    if (bubbles) then
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
                    else if (hypoelasticity) then
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

        if (integral_wrt .and. bubbles) then
            if (n == 0) then ! 1D simulation
                do i = 1, num_integrals
                    int_pres = 0d0
                    max_pres = 0d0
                    k = 0; l = 0
                    npts = 0
                    do j = 1, m
                        pres = 0d0
                        do s = 1, num_dims
                            vel(s) = 0d0
                        end do
                        rho = 0d0
                        pres = 0d0
                        gamma = 0d0
                        pi_inf = 0d0
                        qv = 0d0

                        if ((integral(i)%xmin <= x_cb(j)) .and. (integral(i)%xmax >= x_cb(j))) then
                            npts = npts + 1
                            call s_convert_to_mixture_variables(q_cons_vf, j, k, l, &
                                                                rho, gamma, pi_inf, qv, Re)
                            do s = 1, num_dims
                                vel(s) = q_cons_vf(cont_idx%end + s)%sf(j, k, l)/rho
                            end do

                            pres = ( &
                                   (q_cons_vf(E_idx)%sf(j, k, l) - &
                                    0.5d0*(q_cons_vf(mom_idx%beg)%sf(j, k, l)**2.d0)/rho)/ &
                                   (1.d0 - q_cons_vf(alf_idx)%sf(j, k, l)) - &
                                   pi_inf - qv &
                                   )/gamma
                            int_pres = int_pres + (pres - 1.d0)**2.d0
                        end if
                    end do
                    int_pres = dsqrt(int_pres/(1.d0*npts))

                    if (num_procs > 1) then
                        tmp = int_pres
                        call s_mpi_allreduce_sum(tmp, int_pres)
                    end if

                    if (proc_rank == 0) then
                        if (bubbles .and. (num_fluids <= 2)) then
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
                    int_pres = 0d0
                    max_pres = 0d0
                    l = 0
                    npts = 0
                    do j = 1, m
                        do k = 1, n
                            trigger = .false.
                            if (i == 1) then
                                !inner portion
                                if (dsqrt(x_cb(j)**2.d0 + y_cb(k)**2.d0) < (rad - 0.5d0*thickness)) &
                                    trigger = .true.
                            elseif (i == 2) then
                                !net region
                                if (dsqrt(x_cb(j)**2.d0 + y_cb(k)**2.d0) > (rad - 0.5d0*thickness) .and. &
                                    dsqrt(x_cb(j)**2.d0 + y_cb(k)**2.d0) < (rad + 0.5d0*thickness)) &
                                    trigger = .true.
                            elseif (i == 3) then
                                !everything else
                                if (dsqrt(x_cb(j)**2.d0 + y_cb(k)**2.d0) > (rad + 0.5d0*thickness)) &
                                    trigger = .true.
                            end if

                            pres = 0d0
                            do s = 1, num_dims
                                vel(s) = 0d0
                            end do
                            rho = 0d0
                            pres = 0d0
                            gamma = 0d0
                            pi_inf = 0d0
                            qv = 0d0

                            if (trigger) then
                                npts = npts + 1
                                call s_convert_to_mixture_variables(q_cons_vf, j, k, l, &
                                                                    rho, gamma, pi_inf, qv, Re)
                                do s = 1, num_dims
                                    vel(s) = q_cons_vf(cont_idx%end + s)%sf(j, k, l)/rho
                                end do

                                pres = ( &
                                       (q_cons_vf(E_idx)%sf(j, k, l) - &
                                        0.5d0*(q_cons_vf(mom_idx%beg)%sf(j, k, l)**2.d0)/rho)/ &
                                       (1.d0 - q_cons_vf(alf_idx)%sf(j, k, l)) - &
                                       pi_inf - qv &
                                       )/gamma
                                int_pres = int_pres + abs(pres - 1.d0)
                                max_pres = max(max_pres, abs(pres - 1.d0))
                            end if

                        end do
                    end do

                    if (npts > 0) then
                        int_pres = int_pres/(1.d0*npts)
                    else
                        int_pres = 0.d0
                    end if

                    if (num_procs > 1) then
                        tmp = int_pres
                        call s_mpi_allreduce_sum(tmp, int_pres)

                        tmp = max_pres
                        call s_mpi_allreduce_max(tmp, max_pres)
                    end if

                    if (proc_rank == 0) then
                        if (bubbles .and. (num_fluids <= 2)) then
                            write (i + 70, '(6x,f12.6,f24.8,f24.8)') &
                                nondim_time, int_pres, max_pres
                        end if
                    end if
                end do
            end if
        end if

    end subroutine s_write_probe_files ! -----------------------------------

    @:s_compute_speed_of_sound()

    !>  The goal of this subroutine is to write to the run-time
        !!      information file basic footer information applicable to
        !!      the current computation and to close the file when done.
        !!      The footer contains the stability criteria extrema over
        !!      all of the time-steps and the simulation run-time.
    subroutine s_close_run_time_information_file() ! -----------------------

        real(kind(0d0)) :: run_time !< Run-time of the simulation

        ! Writing the footer of and closing the run-time information file
        write (1, '(A)') '----------------------------------------'// &
            '----------------------------------------'
        write (1, '(A)') ''

        write (1, '(A,F9.6)') 'ICFL Max: ', icfl_max
        if (any(Re_size > 0)) write (1, '(A,F9.6)') 'VCFL Max: ', vcfl_max
        if (any(Re_size > 0)) write (1, '(A,F10.6)') 'Rc Min: ', Rc_min

        call cpu_time(run_time)

        write (1, '(A)') ''
        write (1, '(A,I0,A)') 'Run-time: ', int(anint(run_time)), 's'
        write (1, '(A)') '========================================'// &
            '========================================'
        close (1)

    end subroutine s_close_run_time_information_file ! ---------------------

    !> Closes probe files
    subroutine s_close_probe_files() ! -------------------------------------

        integer :: i !< Generic loop iterator

        do i = 1, num_probes
            close (i + 30)
        end do

    end subroutine s_close_probe_files ! -----------------------------------

    !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_data_output_module() ! -------------------------

        type(int_bounds_info) :: ix, iy, iz

        integer :: i !< Generic loop iterator

        ! Allocating/initializing ICFL, VCFL, CCFL and Rc stability criteria
        @:ALLOCATE(icfl_sf(0:m, 0:n, 0:p))
        icfl_max = 0d0

        if (any(Re_size > 0)) then
            @:ALLOCATE(vcfl_sf(0:m, 0:n, 0:p))
            @:ALLOCATE(Rc_sf  (0:m, 0:n, 0:p))

            vcfl_max = 0d0
            Rc_min = 1d3
        end if

        ! Associating the procedural pointer to the appropriate subroutine
        ! that will be utilized in the conversion to the mixture variables

        if (model_eqns == 1) then        ! Gamma/pi_inf model
            s_convert_to_mixture_variables => &
                s_convert_mixture_to_mixture_variables
        elseif (bubbles) then           ! Volume fraction for bubbles
            s_convert_to_mixture_variables => &
                s_convert_species_to_mixture_variables_bubbles
        else                            ! Volume fraction model
            s_convert_to_mixture_variables => &
                s_convert_species_to_mixture_variables
        end if

        if (parallel_io .neqv. .true.) then
            s_write_data_files => s_write_serial_data_files
        else
            s_write_data_files => s_write_parallel_data_files
        end if

    end subroutine s_initialize_data_output_module ! -----------------------

    !> Module deallocation and/or disassociation procedures
    subroutine s_finalize_data_output_module() ! ---------------------------

        integer :: i !< Generic loop iterator

        ! Deallocating the ICFL, VCFL, CCFL, and Rc stability criteria
        @:DEALLOCATE(icfl_sf)
        if (any(Re_size > 0)) then
            @:DEALLOCATE(vcfl_sf, Rc_sf)
        end if

        ! Disassociating the pointer to the procedure that was utilized to
        ! to convert mixture or species variables to the mixture variables
        s_convert_to_mixture_variables => null()
        s_write_data_files => null()

    end subroutine s_finalize_data_output_module ! -------------------------

end module m_data_output
