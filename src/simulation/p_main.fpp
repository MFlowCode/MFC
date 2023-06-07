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

    use m_global_parameters    !< Definitions of the global parameters

    use m_start_up

    ! ==========================================================================

    implicit none

    integer :: t_step !< Iterator for the time-stepping loop
    real(kind(0d0)) :: time_avg, time_final
    real(kind(0d0)) :: io_time_avg, io_time_final
    real(kind(0d0)), allocatable, dimension(:) :: proc_time
    real(kind(0d0)), allocatable, dimension(:) :: io_proc_time
    logical :: file_exists
    real(kind(0d0)) :: start, finish
    integer :: nt

    call system_clock(COUNT=cpu_start, COUNT_RATE=cpu_rate)

    !Initialize MPI
    call s_initialize_mpi_domain()

    !Initialize Modules
    call s_initialize_modules()

    allocate (proc_time(0:num_procs - 1))
    allocate (io_proc_time(0:num_procs - 1))

    call s_initialize_gpu_vars()

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
        if (t_step == t_step_stop) then
            call s_save_performance_metrics(t_step, time_avg, time_final, io_time_avg, &
                io_time_final, proc_time, io_proc_time, file_exists, start, finish, nt)
            exit 
        end if

        call s_perform_time_step(t_step, time_avg, time_final, io_time_avg, io_time_final, &
         proc_time, io_proc_time, file_exists, start, finish, nt)

        if (mod(t_step - t_step_start, t_step_save) == 0 .or. t_step == t_step_stop) then
            call s_save_data(t_step, start, finish, io_time_avg, nt)
        end if

        call system_clock(cpu_end)
    end do
    ! ==========================================================================

    deallocate(proc_time, io_proc_time)

    call s_finalize_modules()

end program p_main
