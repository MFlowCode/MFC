!>
!! @file m_data_output.f90
!! @brief Contains module m_data_output

#:include 'inline_conversions.fpp'

!> @brief This module takes care of writing the grid and initial condition
!!              data files into the "0" time-step directory located in the folder
!!              associated with the rank of the local processor, which is a sub-
!!              directory of the case folder specified by the user in the input
!!              file pre_process.inp.
module m_data_output

    ! Dependencies =============================================================
    use m_derived_types         !< Definitions of the derived types

    use m_global_parameters     !< Global parameters for the code

    use m_helper

    use m_mpi_proxy             !< Message passing interface (MPI) module proxy

#ifdef MFC_MPI
    use mpi                     !< Message passing interface (MPI) module
#endif

    use m_compile_specific

    use m_variables_conversion

    use m_helper
    ! ==========================================================================

    implicit none

    private; public :: s_write_serial_data_files, &
 s_write_parallel_data_files, &
 s_write_data_files, &
 s_initialize_data_output_module, &
 s_finalize_data_output_module

    abstract interface ! ===================================================

        !>  Interface for the conservative data
        !! @param q_cons_vf The conservative variables
        subroutine s_write_abstract_data_files(q_cons_vf)

            import :: scalar_field, sys_size

            ! Conservative variables
            type(scalar_field), &
                dimension(sys_size), &
                intent(IN) :: q_cons_vf

        end subroutine s_write_abstract_data_files ! -------------------
    end interface ! ========================================================

    character(LEN=path_len + 2*name_len), private :: t_step_dir !<
    !! Time-step folder into which grid and initial condition data will be placed

    character(LEN=path_len + 2*name_len), public :: restart_dir !<
    !! Restart data folder

    procedure(s_write_abstract_data_files), pointer :: s_write_data_files => null()

contains

    !>  Writes grid and initial condition data files to the "0"
        !!  time-step directory in the local processor rank folder
        !! @param q_cons_vf The conservative variables
    subroutine s_write_serial_data_files(q_cons_vf) ! -----------
        type(scalar_field), &
            dimension(sys_size), &
            intent(IN) :: q_cons_vf

        logical :: file_exist !< checks if file exists

        character(LEN=15) :: FMT
        character(LEN=3) :: status

        character(LEN= &
                  int(floor(log10(real(sys_size, wp)))) + 1) :: file_num !< Used to store
            !! the number, in character form, of the currently
            !! manipulated conservative variable data file

        character(LEN=len_trim(t_step_dir) + name_len) :: file_loc !<
            !! Generic string used to store the address of a particular file

        integer :: i, j, k, l !< Generic loop iterator
        integer :: t_step

        real(wp), dimension(nb) :: nRtmp         !< Temporary bubble concentration
        real(wp) :: nbub                         !< Temporary bubble number density
        real(wp) :: gamma, lit_gamma, pi_inf     !< Temporary EOS params
        real(wp) :: rho                          !< Temporary density
        real(wp) :: pres                         !< Temporary pressure

        real(wp) :: nR3
        real(wp) :: ntmp

        t_step = 0

        ! Outputting the Locations of the Cell-boundaries ==================

        if (old_grid) then
            status = 'old'
        else
            status = 'new'
        end if

        ! x-coordinate direction
        file_loc = trim(t_step_dir)//'/x_cb.dat'
        open (1, FILE=trim(file_loc), FORM='unformatted', STATUS=status)
        write (1) x_cb(-1:m)
        close (1)

        ! y- and z-coordinate directions
        if (n > 0) then
            ! y-coordinate direction
            file_loc = trim(t_step_dir)//'/y_cb.dat'
            open (1, FILE=trim(file_loc), FORM='unformatted', &
                  STATUS=status)
            write (1) y_cb(-1:n)
            close (1)

            ! z-coordinate direction
            if (p > 0) then
                file_loc = trim(t_step_dir)//'/z_cb.dat'
                open (1, FILE=trim(file_loc), FORM='unformatted', &
                      STATUS=status)
                write (1) z_cb(-1:p)
                close (1)
            end if
        end if
        ! ==================================================================

        ! Outputting Conservative Variables ================================
        do i = 1, sys_size
            write (file_num, '(I0)') i
            file_loc = trim(t_step_dir)//'/q_cons_vf'//trim(file_num) &
                       //'.dat'
            open (1, FILE=trim(file_loc), FORM='unformatted', &
                  STATUS=status)
            write (1) q_cons_vf(i)%sf
            close (1)
        end do
        ! ==================================================================

        gamma = fluid_pp(1)%gamma
        lit_gamma = 1._wp/fluid_pp(1)%gamma + 1._wp
        pi_inf = fluid_pp(1)%pi_inf

        if (precision == 1) then
            FMT = "(2F30.3)"
        else
            FMT = "(2F40.14)"
        end if

        write (t_step_dir, '(A,I0,A,I0)') trim(case_dir)//'/D'
        file_loc = trim(t_step_dir)//'/.'

        inquire (FILE=trim(file_loc), EXIST=file_exist)

        if (.not. file_exist) call s_create_directory(trim(t_step_dir))

        !1D
        if (n == 0 .and. p == 0) then
            if (model_eqns == 2) then
                do i = 1, sys_size
                    write (file_loc, '(A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir)//'/prim.', i, '.', proc_rank, '.', t_step, '.dat'

                    open (2, FILE=trim(file_loc))
                    do j = 0, m
                        call s_convert_to_mixture_variables(q_cons_vf, j, 0, 0, rho, gamma, pi_inf)

                        lit_gamma = 1._wp/gamma + 1._wp

                        if (((i >= cont_idx%beg) .and. (i <= cont_idx%end)) &
                            .or. &
                            ((i >= adv_idx%beg) .and. (i <= adv_idx%end)) &
                            ) then
                            write (2, FMT) x_cb(j), q_cons_vf(i)%sf(j, 0, 0)
                        else if (i == mom_idx%beg) then !u
                            write (2, FMT) x_cb(j), q_cons_vf(mom_idx%beg)%sf(j, 0, 0)/rho
                        else if (i == stress_idx%beg) then !tau_e
                            write (2, FMT) x_cb(j), q_cons_vf(stress_idx%beg)%sf(j, 0, 0)/rho
                        else if (i == E_idx) then !p
                            call s_compute_pressure( &
                                q_cons_vf(E_idx)%sf(j, 0, 0), &
                                q_cons_vf(alf_idx)%sf(j, 0, 0), &
                                0.5_wp*(q_cons_vf(mom_idx%beg)%sf(j, 0, 0)**2._wp)/rho, &
                                pi_inf, gamma, rho, pres)
                            write (2, FMT) x_cb(j), pres
                        else if ((i >= bub_idx%beg) .and. (i <= bub_idx%end) .and. bubbles) then
                            do k = 1, nb
                                nRtmp(k) = q_cons_vf(bub_idx%rs(k))%sf(j, 0, 0)
                            end do
                            
                            call s_comp_n_from_cons(q_cons_vf(alf_idx)%sf(j, 0, 0), nRtmp, nbub, weight)

                            write (2, FMT) x_cb(j), q_cons_vf(i)%sf(j, 0, 0)/nbub
                        end if
                    end do
                    close (2)
                end do
            end if

            do i = 1, sys_size
                write (file_loc, '(A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir)//'/cons.', i, '.', proc_rank, '.', t_step, '.dat'

                open (2, FILE=trim(file_loc))
                do j = 0, m
                    write (2, FMT) x_cb(j), q_cons_vf(i)%sf(j, 0, 0)
                end do
                close (2)
            end do
        end if

        if (precision == 1) then
            FMT = "(3F30.7)"
        else
            FMT = "(3F40.14)"
        end if

        ! 2D
        if ((n > 0) .and. (p == 0)) then
            do i = 1, sys_size
                write (file_loc, '(A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir)//'/cons.', i, '.', proc_rank, '.', t_step, '.dat'
                open (2, FILE=trim(file_loc))
                do j = 0, m
                    do k = 0, n
                        write (2, FMT) x_cb(j), y_cb(k), q_cons_vf(i)%sf(j, k, 0)
                    end do
                    write (2, *)
                end do
                close (2)
            end do
        end if

        ! 3D
        if (p > 0) then
            do i = 1, sys_size
                write (file_loc, '(A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir)//'/cons.', i, '.', proc_rank, '.', t_step, '.dat'
                open (2, FILE=trim(file_loc))
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
        end if

    end subroutine s_write_serial_data_files ! ------------------------------------

    !> Writes grid and initial condition data files in parallel to the "0"
        !!  time-step directory in the local processor rank folder
        !! @param q_cons_vf The conservative variables
    subroutine s_write_parallel_data_files(q_cons_vf) ! --

        ! Conservative variables
        type(scalar_field), &
            dimension(sys_size), &
            intent(IN) :: q_cons_vf

#ifdef MFC_MPI

        integer :: ifile, ierr, data_size
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(KIND=MPI_OFFSET_KIND) :: disp
        integer(KIND=MPI_OFFSET_KIND) :: m_MOK, n_MOK, p_MOK
        integer(KIND=MPI_OFFSET_KIND) :: WP_MOK, var_MOK, str_MOK
        integer(KIND=MPI_OFFSET_KIND) :: NVARS_MOK
        integer(KIND=MPI_OFFSET_KIND) :: MOK

        character(LEN=path_len + 2*name_len) :: file_loc
        logical :: file_exist

        ! Generic loop iterator
        integer :: i

        ! Initialize MPI data I/O
        call s_initialize_mpi_data(q_cons_vf)

        ! Open the file to write all flow variables
        write (file_loc, '(I0,A)') t_step_start, '.dat'
        file_loc = trim(restart_dir)//trim(mpiiofs)//trim(file_loc)
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
        NVARS_MOK = int(sys_size, MPI_OFFSET_KIND)

        ! Write the data for each variable
        if (bubbles) then
            do i = 1, sys_size! adv_idx%end
                var_MOK = int(i, MPI_OFFSET_KIND)

                ! Initial displacement to skip at beginning of file
                disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, MPI_IO_DATA%view(i), &
                                       'native', mpi_info_int, ierr)
                call MPI_FILE_WRITE_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                        mpi_p, status, ierr)
            end do
        else
            do i = 1, sys_size !TODO: check if this is right
!            do i = 1, adv_idx%end
                var_MOK = int(i, MPI_OFFSET_KIND)

                ! Initial displacement to skip at beginning of file
                disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, MPI_IO_DATA%view(i), &
                                       'native', mpi_info_int, ierr)
                call MPI_FILE_WRITE_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                        mpi_p, status, ierr)
            end do
        end if

        call MPI_FILE_CLOSE(ifile, ierr)

#endif

    end subroutine s_write_parallel_data_files ! ---------------------------

    !> Computation of parameters, allocation procedures, and/or
        !!              any other tasks needed to properly setup the module
    subroutine s_initialize_data_output_module() ! ----------------------------
        ! Generic string used to store the address of a particular file
        character(LEN=len_trim(case_dir) + 2*name_len) :: file_loc

        ! Generic logical used to check the existence of directories
        logical :: dir_check

        if (parallel_io .neqv. .true.) then
            ! Setting the address of the time-step directory
            write (t_step_dir, '(A,I0,A)') '/p_all/p', proc_rank, '/0'
            t_step_dir = trim(case_dir)//trim(t_step_dir)

            ! Checking the existence of the time-step directory, removing it, if
            ! it exists, and creating a new copy. Note that if preexisting grid
            ! and/or initial condition data are to be read in from the very same
            ! location, then the above described steps are not executed here but
            ! rather in the module m_start_up.f90.
            if (old_grid .neqv. .true.) then

                file_loc = trim(t_step_dir)//'/'

                call my_inquire(file_loc, dir_check)

                if (dir_check) call s_delete_directory(trim(t_step_dir))

                call s_create_directory(trim(t_step_dir))

            end if

            s_write_data_files => s_write_serial_data_files
        else
            write (restart_dir, '(A)') '/restart_data'
            restart_dir = trim(case_dir)//trim(restart_dir)

            if ((old_grid .neqv. .true.) .and. (proc_rank == 0)) then

                file_loc = trim(restart_dir)//'/'
                call my_inquire(file_loc, dir_check)

                if (dir_check) call s_delete_directory(trim(restart_dir))
                call s_create_directory(trim(restart_dir))
            end if

            call s_mpi_barrier()

            s_write_data_files => s_write_parallel_data_files

        end if

    end subroutine s_initialize_data_output_module ! --------------------------

    !> Resets s_write_data_files pointer
    subroutine s_finalize_data_output_module() ! ---------------------------

        s_write_data_files => null()

    end subroutine s_finalize_data_output_module ! -------------------------

end module m_data_output
