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

    use m_global_parameters    !< Definitions of the global parameters

    use m_start_up

    use m_time_steppers

    use m_nvtx

    implicit none

    integer :: t_step !< Iterator for the time-stepping loop
    real(wp) :: time_avg, time_final
    real(wp) :: io_time_avg, io_time_final
    real(wp), allocatable, dimension(:) :: proc_time
    real(wp), allocatable, dimension(:) :: io_proc_time
    logical :: file_exists
    real(wp) :: start, finish
    integer :: nt

    call system_clock(COUNT=cpu_start, COUNT_RATE=cpu_rate)

    call nvtxStartRange("INIT")

    !Initialize MPI
    call nvtxStartRange("INIT-MPI")
    call s_initialize_mpi_domain()
    call nvtxEndRange

    !Initialize Modules
    call nvtxStartRange("INIT-MODULES")
    call s_initialize_modules()
    call nvtxEndRange

    allocate (proc_time(0:num_procs - 1))
    allocate (io_proc_time(0:num_procs - 1))

    call nvtxStartRange("INIT-GPU-VARS")
    call s_initialize_gpu_vars()
    call nvtxEndRange

    ! Setting the time-step iterator to the first time-step
    if (cfl_dt) then
        t_step = 0
        mytime = t_save*n_start
    else
        t_step = t_step_start
        if (t_step == 0) then
            mytime = 0._wp
        else
            mytime = t_step*dt
        end if
        finaltime = t_step_stop*dt
    end if

    call nvtxEndRange ! INIT

    call nvtxStartRange("SIMULATION-TIME-MARCH")
    ! Time-stepping Loop
    do

        if (cfl_dt) then
            if (mytime >= t_stop) then
                call s_save_performance_metrics(time_avg, time_final, io_time_avg, &
                                                io_time_final, proc_time, io_proc_time, file_exists)
                exit
            end if
        else
            if (t_step == t_step_stop) then
                call s_save_performance_metrics(time_avg, time_final, io_time_avg, &
                                                io_time_final, proc_time, io_proc_time, file_exists)
                exit
            end if
        end if

        call s_perform_time_step(t_step, time_avg)

        if (cfl_dt) then
            if (abs(mod(mytime, t_save)) < dt .or. mytime >= t_stop) then
                call s_save_data(t_step, start, finish, io_time_avg, nt)
            end if
        else
            if (mod(t_step - t_step_start, t_step_save) == 0 .or. t_step == t_step_stop) then
                call s_save_data(t_step, start, finish, io_time_avg, nt)
            end if
        end if

        call system_clock(cpu_end)
    end do

    call nvtxEndRange ! Simulation

    deallocate (proc_time, io_proc_time)

    call nvtxStartRange("FINALIZE-MODULES")
    call s_finalize_modules()
    call nvtxEndRange

end program p_main
