!>
!! @file p_main.f90
!! @brief Contains program p_main

!> @brief  Quasi-conservative, shock- and interface- capturing finite-volume
!!              scheme for the multicomponent Navier-Stokes equations. The system
!!              is augmented with the relevant advection equations to capture the
!!              material interfaces and closed by the stiffened equation of state
!!              as well as any required mixture relations. The effects of surface
!!              tension are included and modeled through a volume force that acts
!!              across the diffuse material interface regions. The implementation
!!              specifics of surface tension may be found in the work by Perigaud
!!              and Saurel (2005). Note that both viscous and capillarity effects
!!              are only available in the volume fraction model.
program p_main

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_start_up             !< Reading and checking procedures for the input
    !< and the initial condition and grid data files

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

#ifdef _OPENACC
    use openacc
#endif

    use m_nvtx

    ! ==========================================================================

    implicit none

    integer :: err_code, ierr

    integer :: t_step, i !< Iterator for the time-stepping loop
    real(kind(0d0)) :: time_avg, time_final
    real(kind(0d0)) :: io_time_avg, io_time_final
    real(kind(0d0)), allocatable, dimension(:) :: proc_time
    real(kind(0d0)), allocatable, dimension(:) :: io_proc_time
    logical :: file_exists
    real(kind(0d0)) :: start, finish
    integer :: nt

#ifdef _OPENACC
    real(kind(0d0)) :: starttime, endtime
    integer :: num_devices, local_size, num_nodes, ppn, my_device_num
    integer :: dev, devNum, local_rank
#ifdef MFC_MPI
    integer :: local_comm
#endif
    integer(acc_device_kind) :: devtype
#endif

    call system_clock(COUNT=cpu_start, COUNT_RATE=cpu_rate)

    ! Initializing MPI execution environment

    call s_mpi_initialize()

! Bind GPUs if OpenACC is enabled
#ifdef _OPENACC
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

    ! Computation of parameters, allocation of memory, association of pointers,
    ! and/or the execution of any other tasks needed to properly configure the
    ! modules. The preparations below DO NOT DEPEND on the grid being complete.
    call s_initialize_global_parameters_module()

#if defined(_OPENACC) && defined(MFC_MEMORY_DUMP)
    call acc_present_dump()
#endif // defined(_OPENACC) && defined(MFC_MEMORY_DUMP)

    call s_initialize_mpi_proxy_module()
    call s_initialize_variables_conversion_module()
    if (grid_geometry == 3) call s_initialize_fftw_module()
    call s_initialize_start_up_module()
    call s_initialize_riemann_solvers_module()

    if(bubbles) call s_initialize_bubbles_module()

    if (qbmm) call s_initialize_qbmm_module()

#if defined(_OPENACC) && defined(MFC_MEMORY_DUMP)
    call acc_present_dump()
#endif // defined(_OPENACC) && defined(MFC_MEMORY_DUMP)

    if (monopole) then
        call s_initialize_monopole_module()
    end if
    if (any(Re_size > 1)) then
        call s_initialize_viscous_module()
    end if
    call s_initialize_rhs_module()

#if defined(_OPENACC) && defined(MFC_MEMORY_DUMP)
    call acc_present_dump()
#endif // defined(_OPENACC) && defined(MFC_MEMORY_DUMP)

    if (hypoelasticity) call s_initialize_hypoelastic_module()
    call s_initialize_data_output_module()
    call s_initialize_derived_variables_module()
    call s_initialize_time_steppers_module()

#if defined(_OPENACC) && defined(MFC_MEMORY_DUMP)
    call acc_present_dump()
#endif // defined(_OPENACC) && defined(MFC_MEMORY_DUMP)

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

#if defined(_OPENACC) && defined(MFC_MEMORY_DUMP)
    print *, "[MEM-INST] After: s_initialize_weno_module"
    call acc_present_dump()
#endif // defined(_OPENACC) && defined(MFC_MEMORY_DUMP)

    call s_initialize_cbc_module()

    call s_initialize_derived_variables()

    allocate (proc_time(0:num_procs - 1))
    allocate (io_proc_time(0:num_procs - 1))

!$acc update device(dt, dx, dy, dz, x_cc, y_cc, z_cc, x_cb, y_cb, z_cb)
!$acc update device(sys_size, buff_size)
!$acc update device(m, n, p)
    do i = 1, sys_size
!$acc update device(q_cons_ts(1)%vf(i)%sf)
    end do

    ! Setting the time-step iterator to the first time-step
    t_step = t_step_start
    if (t_step == 0) then
        mytime = 0d0
    else
        mytime = t_step*dt
    end if
    finaltime = t_step_stop*dt

    ! Time-stepping Loop =======================================================
    do
        if (proc_rank == 0) then
            print '(" ["I3"%]  Time step "I8" of "I0" @ t_step = "I0"")',       &
                  int(100*(real(t_step + 1)/(t_step_stop - t_step_start + 1))), &
                  t_step      - t_step_start + 1,                               &
                  t_step_stop - t_step_start + 1,                               &
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

        if (t_step == t_step_stop) then

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

            exit
        else
            if ((mytime + dt) >= finaltime) dt = finaltime - mytime
            t_step = t_step + 1
        end if

        !if (num_procs > 1) then
        !    print*, "m_compress timings:"
        !    print*, " - nCalls_time  ", nCalls_time
        !    print*, " - s_compress   ", (compress_time   / nCalls_time), "s"
        !    print*, " - mpi_sendrecv ", (mpi_time        / nCalls_time), "s"
        !    print*, " - s_decompress ", (decompress_time / nCalls_time), "s"
        !end if

        ! print*, 'Write data files'
        ! Backing up the grid and conservative variables data
        if (mod(t_step - t_step_start, t_step_save) == 0 .or. t_step == t_step_stop) then

            call cpu_time(start)
            !  call nvtxStartRange("I/O")
            do i = 1, sys_size
!$acc update host(q_cons_ts(1)%vf(i)%sf)
            end do
            call s_write_data_files(q_cons_ts(1)%vf, t_step)
            !  call nvtxEndRange
            call cpu_time(finish)
            nt = int((t_step - t_step_start)/(t_step_save))
            if (nt == 1) then
                io_time_avg = abs(finish - start)
            else
                io_time_avg = (abs(finish - start) + io_time_avg*(nt - 1))/nt
            end if
        end if

        call system_clock(cpu_end)
    end do
    ! ==========================================================================

    ! Disassociate pointers for serial and parallel I/O
    s_read_data_files => null()
    s_write_data_files => null()

    ! Deallocation and/or disassociation procedures for the modules

    deallocate (proc_time)

    call s_finalize_time_steppers_module()
    call s_finalize_derived_variables_module()
    call s_finalize_data_output_module()
    call s_finalize_rhs_module()
    call s_finalize_cbc_module()
    call s_finalize_riemann_solvers_module()
    call s_finalize_weno_module()
    call s_finalize_start_up_module()
    call s_finalize_variables_conversion_module()
    if (grid_geometry == 3) call s_finalize_fftw_module
    call s_finalize_mpi_proxy_module()
    call s_finalize_global_parameters_module()

    if (any(Re_size > 0)) then
        call s_finalize_viscous_module()
    end if

    ! Terminating MPI execution environment
    call s_mpi_finalize()

end program p_main
