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

    ! Dependencies =============================================================
    use m_global_parameters     !< Global parameters for the code

    use m_start_up
    ! ==========================================================================

    implicit none

    integer :: t_step !< Iterator for the main time-stepping loop

    character(LEN=name_len) :: varname !<
    !! Generic storage for the name(s) of the flow variable(s) that will be added
    !! to the formatted database file(s)

    real(kind(0d0)) :: pres
    real(kind(0d0)) :: c
    real(kind(0d0)) :: H 

    call s_initialize_mpi_domain()

    call s_initialize_modules()

    ! Setting the time-step iterator to the first time step to be post-processed
    t_step = t_step_start

    ! Time-Marching Loop =======================================================
    do
        call s_perform_time_step(t_step)

        call s_save_data(t_step, varname, pres, c, H)

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

        ! Incrementing time-step iterator to next time-step to be post-processed
        t_step = t_step + t_step_save

    end do
    ! END: Time-Marching Loop ==================================================

    close (11)

    call s_finalize_modules()

end program p_main
