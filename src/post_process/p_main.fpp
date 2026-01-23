!>
!! @file p_main.f90
!! @brief Contains program p_main

!> @brief The post-process restructures raw unformatted data, outputted by
!!              the simulation, into a formatted database, Silo-HDF5 or Binary,
!!              chosen by the user. The user may also specify which variables to
!!              include in the database. The choices range from any one of the
!!              primitive and conservative variables, as well as quantities that
!!              can be derived from those such as the unadvected volume fraction,
!!              specific heat ratio, liquid stiffness, speed of sound, vorticity
!!              and the numerical Schlieren function.
program p_main

    use m_global_parameters     !< Global parameters for the code
    use m_start_up

    implicit none

    integer :: t_step !< Iterator for the main time-stepping loop

    character(LEN=name_len) :: varname !<
    !! Generic storage for the name(s) of the flow variable(s) that will be added
    !! to the formatted database file(s)

    real(wp) :: pres
    real(wp) :: c
    real(wp) :: H
    real(wp) :: start, finish

    call s_initialize_mpi_domain()

    call s_initialize_modules()

    if (cfl_dt) then
        t_step = n_start
        n_save = int(t_stop/t_save) + 1
    else
        ! Setting the time-step iterator to the first time step to be post-processed
        t_step = t_step_start
    end if

    ! Time-Marching Loop
    do

        ! If all time-steps are not ready to be post-processed and one rank is
        ! faster than another, the slower rank processing the last available
        ! step might be killed when the faster rank attempts to process the
        ! first missing step, before the slower rank finishes writing the last
        ! available step. To avoid this, we force synchronization here.
        call s_mpi_barrier()

        call cpu_time(start)

        call s_perform_time_step(t_step)

        call s_save_data(t_step, varname, pres, c, H)

        call cpu_time(finish)

        wall_time = abs(finish - start)

        if (t_step >= 2) then
            wall_time_avg = (wall_time + (t_step - 2)*wall_time_avg)/(t_step - 1)
        else
            wall_time_avg = 0._wp
        end if

        if (cfl_dt) then
            if (t_step == n_save - 1) then
                exit
            end if
        else
            ! Modifies the time-step iterator so that it may reach the final time-
            ! step to be post-processed, in the case that this one is not originally
            ! attainable through constant incrementation from the first time-step.
            ! This modification is performed upon reaching the final time-step. In
            ! case that it is not needed, the post-processor is done and may exit.
            if ((t_step_stop - t_step) < t_step_save .and. t_step_stop /= t_step) then
                t_step = t_step_stop - t_step_save
            elseif (t_step == t_step_stop) then
                exit
            end if
        end if

        if (cfl_dt) then
            t_step = t_step + 1
        else
            ! Incrementing time-step iterator to next time-step to be post-processed
            t_step = t_step + t_step_save
        end if

    end do
    ! END: Time-Marching Loop

    close (11)

    call s_finalize_modules()

end program p_main
