!>
!! @file p_main.f90
!! @brief Contains program p_main

!> @brief This program takes care of setting up the initial condition and
!!              grid data for the multicomponent flow code.
program p_main

    use m_global_parameters     !< Global parameters for the code

    use m_start_up

    implicit none

    logical :: file_exists
    real(wp) :: start, finish, time_avg, time_final
    real(wp), allocatable, dimension(:) :: proc_time

    call random_seed()

    call s_initialize_mpi_domain()

    ! Initialization of the MPI environment

    call s_initialize_modules()

    call s_read_grid()

    allocate (proc_time(0:num_procs - 1))

    call s_apply_initial_condition(start, finish)

    time_avg = abs(finish - start)

    call s_save_data(proc_time, time_avg, time_final, file_exists)

    deallocate (proc_time)

    call s_finalize_modules()

end program p_main
