!>
!! @file p_main.f90
!! @brief Contains program p_main

!> @brief This program takes care of setting up the initial condition and
!!              grid data for the multicomponent flow code.
program p_main

    ! Dependencies =============================================================
    use m_derived_types         !< Definitions of the derived types

    use m_global_parameters     !< Global parameters for the code

    use m_mpi_proxy             !< Message passing interface (MPI) module proxy

    use m_variables_conversion  !< Subroutines to change the state variables from
                                !! one form to another

    use m_start_up              !< Procedures to read in and check consistency of
                                !! the user provided inputs and data

    use m_grid                  !< Procedures to generate (non-)uniform grids

    use m_initial_condition     !< Procedures to generate initial condition

    use m_data_output           !< Procedures to write the grid data and the
                                !! conservative variables to files

    use m_compile_specific      !< Compile-specific procedures

    use m_patches

    use m_assign_variables
    ! ==========================================================================

    implicit none

    integer :: i
    logical :: file_exists
    real(wp) :: start, finish, time_avg, time_final
    real(wp), allocatable, dimension(:) :: proc_time

    ! Initialization of the MPI environment

    call s_mpi_initialize()

    ! Rank 0 processor assigns default values to user inputs prior to reading
    ! those in from the input file. Next, the user inputs are read in and their
    ! consistency is checked. The detection of any inconsistencies automatically
    ! leads to the termination of the pre-process.

    if (proc_rank == 0) then
        call s_assign_default_values_to_user_inputs()
        call s_read_input_file()
        call s_check_input_file()
    
        print '(" Pre-processing a "I0"x"I0"x"I0" case on "I0" rank(s)")', m, n, p, num_procs
    end if

    ! Broadcasting the user inputs to all of the processors and performing the
    ! parallel computational domain decomposition. Neither procedure has to be
    ! carried out if pre-process is in fact not truly executed in parallel.
    call s_mpi_bcast_user_inputs()
    call s_initialize_parallel_io()
    call s_mpi_decompose_computational_domain()

    ! Computation of parameters, allocation procedures, and/or any other tasks
    ! needed to properly setup the modules
    call s_initialize_global_parameters_module()
    call s_initialize_data_output_module()
    call s_initialize_variables_conversion_module()
    call s_initialize_start_up_module()
    call s_initialize_grid_module()
    call s_initialize_initial_condition_module()
    call s_initialize_assign_variables_module()

    ! Associate pointers for serial or parallel I/O
    if (parallel_io .neqv. .true.) then
        s_generate_grid => s_generate_serial_grid
        s_read_grid_data_files => s_read_serial_grid_data_files
        s_read_ic_data_files => s_read_serial_ic_data_files
        s_write_data_files => s_write_serial_data_files

        ! Create the D directory if it doesn't exit, to store
        ! the serial data files
        call s_create_directory('D')
    else
        s_generate_grid => s_generate_parallel_grid
        s_read_grid_data_files => s_read_parallel_grid_data_files
        s_read_ic_data_files => s_read_parallel_ic_data_files
        s_write_data_files => s_write_parallel_data_files
    end if

    ! Setting up the grid and the initial condition. If the grid is read in from
    ! preexisting grid data files, it is checked for consistency. If the grid is
    ! not read in, it is generated from scratch according to the inputs provided
    ! by the user. The initial condition may also be read in. It in turn is not
    ! checked for consistency since it WILL further be edited by the pre-process
    ! and also because it may be incomplete at the time it is read in. Finally,
    ! when the grid and initial condition are completely setup, they are written
    ! to their respective data files.

    ! Setting up grid and initial condition

    allocate (proc_time(0:num_procs - 1))

    call cpu_time(start)

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

    if (old_ic) call s_read_ic_data_files(q_cons_vf)

    call s_generate_initial_condition()

    call s_write_data_files(q_cons_vf)

    call cpu_time(finish)

    time_avg = abs(finish - start)

    call s_mpi_barrier()

    if (num_procs > 1) then
        call mpi_bcast_time_step_values(proc_time, time_avg)
    end if

    if (proc_rank == 0) then
        time_final = 0._wp
        if (num_procs == 1) then
            time_final = time_avg
            print *, "Final Time", time_final
        else
            time_final = maxval(proc_time)
            print *, "Final Time", time_final
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

    ! Disassociate pointers for serial and parallel I/O
    s_generate_grid => null()
    s_read_grid_data_files => null()
    s_read_ic_data_files => null()
    s_write_data_files => null()

    ! Deallocation procedures for the modules
    call s_finalize_grid_module()
    call s_finalize_start_up_module()
    call s_finalize_variables_conversion_module()
    call s_finalize_data_output_module()
    call s_finalize_global_parameters_module()
    call s_finialize_assign_variables_module()

    ! Finalization of the MPI environment
    call s_mpi_finalize()

end program p_main
