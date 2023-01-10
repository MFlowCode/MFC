!>
!! @file m_data_input.f90
!> @brief Contains module m_data_input

!> @brief This module features procedures, which for a specific time-step,
!!             read in the raw simulation data for the grid and the conservative
!!             variables and fill out their buffer regions.
module m_data_input

    ! Dependencies =============================================================
#ifdef MFC_MPI
    use mpi                     !< Message passing interface (MPI) module
#endif

    use m_derived_types         !< Definitions of the derived types

    use m_global_parameters     !< Global parameters for the code

    use m_mpi_proxy             !< Message passing interface (MPI) module proxy

    use m_compile_specific
    ! ==========================================================================

    implicit none

    private; public :: s_initialize_data_input_module, &
 s_read_data_files, &
 s_read_serial_data_files, &
 s_read_parallel_data_files, &
 s_populate_grid_variables_buffer_regions, &
 s_populate_conservative_variables_buffer_regions, &
 s_finalize_data_input_module

    abstract interface ! ===================================================

        !> Subroutine for reading data files
        !!  @param t_step Current time-step to input
        subroutine s_read_abstract_data_files(t_step) ! ------------

            integer, intent(IN) :: t_step

        end subroutine s_read_abstract_data_files ! ----------------

    end interface ! ========================================================

    type(scalar_field), allocatable, dimension(:), public :: q_cons_vf !<
    !! Conservative variables

    type(scalar_field), allocatable, dimension(:), public :: q_prim_vf !<
    !! Primitive variables

    procedure(s_read_abstract_data_files), pointer :: s_read_data_files => null()

contains

    !>  This subroutine is called at each time-step that has to
        !!      be post-processed in order to read the raw data files
        !!      present in the corresponding time-step directory and to
        !!      populate the associated grid and conservative variables.
        !!  @param t_step Current time-step
    subroutine s_read_serial_data_files(t_step) ! -----------------------------

        integer, intent(IN) :: t_step

        character(LEN=len_trim(case_dir) + 2*name_len) :: t_step_dir !<
            !! Location of the time-step directory associated with t_step

        character(LEN=len_trim(case_dir) + 3*name_len) :: file_loc !<
            !! Generic string used to store the location of a particular file

        character(LEN= &
                  int(floor(log10(real(sys_size, kind(0d0))))) + 1) :: file_num !<
            !! Used to store the variable position, in character form, of the
            !! currently manipulated conservative variable file

        logical :: dir_check !<
            !! Generic logical used to test the existence of a particular folder

        logical :: file_check  !<
            !! Generic logical used to test the existence of a particular file

        integer :: i !< Generic loop iterator

        ! Setting location of time-step folder based on current time-step
        write (t_step_dir, '(A,I0,A,I0)') '/p_all/p', proc_rank, '/', t_step
        t_step_dir = trim(case_dir)//trim(t_step_dir)

        ! Inquiring as to the existence of the time-step directory
        file_loc = trim(t_step_dir)//'/.'

        call my_inquire(file_loc, dir_check)

        ! If the time-step directory is missing, the post-process exits.
        if (dir_check .neqv. .true.) then
            print '(A)', 'Time-step folder '//trim(t_step_dir)// &
                ' is missing. Exiting ...'
            call s_mpi_abort()
        end if

        ! Reading the Grid Data File for the x-direction ===================

        ! Checking whether x_cb.dat exists
        file_loc = trim(t_step_dir)//'/x_cb.dat'
        inquire (FILE=trim(file_loc), EXIST=file_check)

        ! Reading x_cb.dat if it exists, exiting otherwise
        if (file_check) then
            open (1, FILE=trim(file_loc), FORM='unformatted', &
                  STATUS='old', ACTION='read')
            read (1) x_cb(-1:m)
            close (1)
        else
            print '(A)', 'File x_cb.dat is missing in '// &
                trim(t_step_dir)//'. Exiting ...'
            call s_mpi_abort()
        end if

        ! Computing the cell-width distribution
        dx(0:m) = x_cb(0:m) - x_cb(-1:m - 1)

        ! Computing the cell-center locations
        x_cc(0:m) = x_cb(-1:m - 1) + dx(0:m)/2d0

        ! ==================================================================

        ! Reading the Grid Data File for the y-direction ===================

        if (n > 0) then

            ! Checking whether y_cb.dat exists
            file_loc = trim(t_step_dir)//'/y_cb.dat'
            inquire (FILE=trim(file_loc), EXIST=file_check)

            ! Reading y_cb.dat if it exists, exiting otherwise
            if (file_check) then
                open (1, FILE=trim(file_loc), FORM='unformatted', &
                      STATUS='old', ACTION='read')
                read (1) y_cb(-1:n)
                close (1)
            else
                print '(A)', 'File y_cb.dat is missing in '// &
                    trim(t_step_dir)//'. Exiting ...'
                call s_mpi_abort()
            end if

            ! Computing the cell-width distribution
            dy(0:n) = y_cb(0:n) - y_cb(-1:n - 1)

            ! Computing the cell-center locations
            y_cc(0:n) = y_cb(-1:n - 1) + dy(0:n)/2d0

            ! ==================================================================

            ! Reading the Grid Data File for the z-direction ===================

            if (p > 0) then

                ! Checking whether z_cb.dat exists
                file_loc = trim(t_step_dir)//'/z_cb.dat'
                inquire (FILE=trim(file_loc), EXIST=file_check)

                ! Reading z_cb.dat if it exists, exiting otherwise
                if (file_check) then
                    open (1, FILE=trim(file_loc), FORM='unformatted', &
                          STATUS='old', ACTION='read')
                    read (1) z_cb(-1:p)
                    close (1)
                else
                    print '(A)', 'File z_cb.dat is missing in '// &
                        trim(t_step_dir)//'. Exiting ...'
                    call s_mpi_abort()
                end if

                ! Computing the cell-width distribution
                dz(0:p) = z_cb(0:p) - z_cb(-1:p - 1)

                ! Computing the cell-center locations
                z_cc(0:p) = z_cb(-1:p - 1) + dz(0:p)/2d0

            end if

        end if

        ! ==================================================================

        ! Reading the Conservative Variables Data Files ====================
        do i = 1, sys_size

            ! Checking whether the data file associated with the variable
            ! position of currently manipulated conservative variable exists
            write (file_num, '(I0)') i
            file_loc = trim(t_step_dir)//'/q_cons_vf'// &
                       trim(file_num)//'.dat'
            inquire (FILE=trim(file_loc), EXIST=file_check)

            ! Reading the data file if it exists, exiting otherwise
            if (file_check) then
                open (1, FILE=trim(file_loc), FORM='unformatted', &
                      STATUS='old', ACTION='read')
                read (1) q_cons_vf(i)%sf(0:m, 0:n, 0:p)
                close (1)
            else
                print '(A)', 'File q_cons_vf'//trim(file_num)// &
                    '.dat is missing in '//trim(t_step_dir)// &
                    '. Exiting ...'
                call s_mpi_abort()
            end if

        end do

        ! ==================================================================

    end subroutine s_read_serial_data_files ! ---------------------------------

    !>  This subroutine is called at each time-step that has to
        !!      be post-processed in order to parallel-read the raw data files
        !!      present in the corresponding time-step directory and to
        !!      populate the associated grid and conservative variables.
        !!  @param t_step Current time-step
    subroutine s_read_parallel_data_files(t_step) ! ---------------------------

        integer, intent(IN) :: t_step

#ifdef MFC_MPI

        real(kind(0d0)), allocatable, dimension(:) :: x_cb_glb, y_cb_glb, z_cb_glb

        integer :: ifile, ierr, data_size
        integer, dimension(MPI_STATUS_SIZE) :: status
        real(kind(0d0)) :: start, finish
        integer(KIND=MPI_OFFSET_KIND) :: disp
        integer(KIND=MPI_OFFSET_KIND) :: m_MOK, n_MOK, p_MOK
        integer(KIND=MPI_OFFSET_KIND) :: WP_MOK, var_MOK, str_MOK
        integer(KIND=MPI_OFFSET_KIND) :: NVARS_MOK
        integer(KIND=MPI_OFFSET_KIND) :: MOK

        character(LEN=path_len + 2*name_len) :: file_loc
        logical :: file_exist

        integer :: i

        allocate (x_cb_glb(-1:m_glb))
        allocate (y_cb_glb(-1:n_glb))
        allocate (z_cb_glb(-1:p_glb))

        ! Read in cell boundary locations in x-direction
        file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//'x_cb.dat'
        inquire (FILE=trim(file_loc), EXIST=file_exist)

        if (file_exist) then
            data_size = m_glb + 2
            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
            call MPI_FILE_READ(ifile, x_cb_glb, data_size, MPI_DOUBLE_PRECISION, status, ierr)
            call MPI_FILE_CLOSE(ifile, ierr)
        else
            print '(A)', 'File ', trim(file_loc), ' is missing. Exiting...'
            call s_mpi_abort()
        end if

        ! Assigning local cell boundary locations
        x_cb(-1:m) = x_cb_glb((start_idx(1) - 1):(start_idx(1) + m))
        ! Computing the cell width distribution
        dx(0:m) = x_cb(0:m) - x_cb(-1:m - 1)
        ! Computing the cell center location
        x_cc(0:m) = x_cb(-1:m - 1) + dx(0:m)/2d0

        if (n > 0) then
            ! Read in cell boundary locations in y-direction
            file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//'y_cb.dat'
            inquire (FILE=trim(file_loc), EXIST=file_exist)

            if (file_exist) then
                data_size = n_glb + 2
                call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
                call MPI_FILE_READ(ifile, y_cb_glb, data_size, MPI_DOUBLE_PRECISION, status, ierr)
                call MPI_FILE_CLOSE(ifile, ierr)
            else
                print '(A)', 'File ', trim(file_loc), ' is missing. Exiting...'
                call s_mpi_abort()
            end if

            ! Assigning local cell boundary locations
            y_cb(-1:n) = y_cb_glb((start_idx(2) - 1):(start_idx(2) + n))
            ! Computing the cell width distribution
            dy(0:n) = y_cb(0:n) - y_cb(-1:n - 1)
            ! Computing the cell center location
            y_cc(0:n) = y_cb(-1:n - 1) + dy(0:n)/2d0

            if (p > 0) then
                ! Read in cell boundary locations in z-direction
                file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//'z_cb.dat'
                inquire (FILE=trim(file_loc), EXIST=file_exist)

                if (file_exist) then
                    data_size = p_glb + 2
                    call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
                    call MPI_FILE_READ(ifile, z_cb_glb, data_size, MPI_DOUBLE_PRECISION, status, ierr)
                    call MPI_FILE_CLOSE(ifile, ierr)
                else
                    print '(A)', 'File ', trim(file_loc), ' is missing. Exiting...'
                    call s_mpi_abort()
                end if

                ! Assigning local cell boundary locations
                z_cb(-1:p) = z_cb_glb((start_idx(3) - 1):(start_idx(3) + p))
                ! Computing the cell width distribution
                dz(0:p) = z_cb(0:p) - z_cb(-1:p - 1)
                ! Computing the cell center location
                z_cc(0:p) = z_cb(-1:p - 1) + dz(0:p)/2d0
            end if
        end if

        ! Open the file to read conservative variables
        write (file_loc, '(I0,A)') t_step, '.dat'
        file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
        inquire (FILE=trim(file_loc), EXIST=file_exist)

        if (file_exist) then
            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

            ! Initialize MPI data I/O
            call s_initialize_mpi_data(q_cons_vf)

            ! Size of local arrays
            data_size = (m + 1)*(n + 1)*(p + 1)

            ! Resize some integers so MPI can read even the biggest file
            m_MOK = int(m_glb + 1, MPI_OFFSET_KIND)
            n_MOK = int(n_glb + 1, MPI_OFFSET_KIND)
            p_MOK = int(p_glb + 1, MPI_OFFSET_KIND)
            WP_MOK = int(8d0, MPI_OFFSET_KIND)
            MOK = int(1d0, MPI_OFFSET_KIND)
            str_MOK = int(name_len, MPI_OFFSET_KIND)
            NVARS_MOK = int(sys_size, MPI_OFFSET_KIND)

            ! Read the data for each variable
            if (bubbles .or. hypoelasticity) then
                do i = 1, sys_size
                    var_MOK = int(i, MPI_OFFSET_KIND)

                    ! Initial displacement to skip at beginning of file
                    disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                    call MPI_FILE_SET_VIEW(ifile, disp, MPI_DOUBLE_PRECISION, MPI_IO_DATA%view(i), &
                                           'native', mpi_info_int, ierr)
                    call MPI_FILE_READ_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                           MPI_DOUBLE_PRECISION, status, ierr)
                end do
            else
                do i = 1, adv_idx%end
                    var_MOK = int(i, MPI_OFFSET_KIND)

                    ! Initial displacement to skip at beginning of file
                    disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                    call MPI_FILE_SET_VIEW(ifile, disp, MPI_DOUBLE_PRECISION, MPI_IO_DATA%view(i), &
                                           'native', mpi_info_int, ierr)
                    call MPI_FILE_READ_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                           MPI_DOUBLE_PRECISION, status, ierr)
                end do
            end if

            call s_mpi_barrier()

            call MPI_FILE_CLOSE(ifile, ierr)
        else
            print '(A)', 'File ', trim(file_loc), ' is missing. Exiting...'
            call s_mpi_abort()
        end if

        deallocate (x_cb_glb, y_cb_glb, z_cb_glb)

#endif

    end subroutine s_read_parallel_data_files ! -------------------------------

    !>  The following subroutine populates the buffer regions of
        !!      the cell-width spacings, the cell-boundary locations and
        !!      the cell-center locations. Note that the buffer regions
        !!      of the last two variables should be interpreted slightly
        !!      differently than usual. They are really ghost zones that
        !!      are used in aiding the multidimensional visualization of
        !!      Silo database files, in VisIt, when processor boundary
        !!      conditions are present.
    subroutine s_populate_grid_variables_buffer_regions() ! ----------------

        integer :: i !< Generic loop iterator

        ! Populating Buffer Regions in the x-direction =====================

        ! Ghost-cell extrapolation BC at the beginning
        if (bc_x%beg <= -3) then

            do i = 1, buff_size
                dx(-i) = dx(0)
            end do

            ! Symmetry BC at the beginning
        elseif (bc_x%beg == -2) then

            do i = 1, buff_size
                dx(-i) = dx(i - 1)
            end do

            ! Periodic BC at the beginning
        elseif (bc_x%beg == -1) then

            do i = 1, buff_size
                dx(-i) = dx((m + 1) - i)
            end do

            ! Processor BC at the beginning
        else

            call s_mpi_sendrecv_grid_vars_buffer_regions('beg', 'x')

            do i = 1, offset_x%beg
                x_cb(-1 - i) = x_cb(-i) - dx(-i)
            end do

            do i = 1, buff_size
                x_cc(-i) = x_cc(1 - i) - (dx(1 - i) + dx(-i))/2d0
            end do

        end if

        ! Ghost-cell extrapolation BC at the end
        if (bc_x%end <= -3) then

            do i = 1, buff_size
                dx(m + i) = dx(m)
            end do

            ! Symmetry BC at the end
        elseif (bc_x%end == -2) then

            do i = 1, buff_size
                dx(m + i) = dx((m + 1) - i)
            end do

            ! Periodic BC at the end
        elseif (bc_x%end == -1) then

            do i = 1, buff_size
                dx(m + i) = dx(i - 1)
            end do

            ! Processor BC at the end
        else

            call s_mpi_sendrecv_grid_vars_buffer_regions('end', 'x')

            do i = 1, offset_x%end
                x_cb(m + i) = x_cb(m + (i - 1)) + dx(m + i)
            end do

            do i = 1, buff_size
                x_cc(m + i) = x_cc(m + (i - 1)) + (dx(m + (i - 1)) + dx(m + i))/2d0
            end do

        end if

        ! END: Populating Buffer Regions in the x-direction ================

        ! Populating Buffer Regions in the y-direction =====================

        if (n > 0) then

            ! Ghost-cell extrapolation BC at the beginning
            if (bc_y%beg <= -3 .and. bc_y%beg /= -13) then

                do i = 1, buff_size
                    dy(-i) = dy(0)
                end do

                ! Symmetry BC at the beginning
            elseif (bc_y%beg == -2 .or. bc_y%beg == -13) then

                do i = 1, buff_size
                    dy(-i) = dy(i - 1)
                end do

                ! Periodic BC at the beginning
            elseif (bc_y%beg == -1) then

                do i = 1, buff_size
                    dy(-i) = dy((n + 1) - i)
                end do

                ! Processor BC at the beginning
            else

                call s_mpi_sendrecv_grid_vars_buffer_regions('beg', 'y')

                do i = 1, offset_y%beg
                    y_cb(-1 - i) = y_cb(-i) - dy(-i)
                end do

                do i = 1, buff_size
                    y_cc(-i) = y_cc(1 - i) - (dy(1 - i) + dy(-i))/2d0
                end do

            end if

            ! Ghost-cell extrapolation BC at the end
            if (bc_y%end <= -3) then

                do i = 1, buff_size
                    dy(n + i) = dy(n)
                end do

                ! Symmetry BC at the end
            elseif (bc_y%end == -2) then

                do i = 1, buff_size
                    dy(n + i) = dy((n + 1) - i)
                end do

                ! Periodic BC at the end
            elseif (bc_y%end == -1) then

                do i = 1, buff_size
                    dy(n + i) = dy(i - 1)
                end do

                ! Processor BC at the end
            else

                call s_mpi_sendrecv_grid_vars_buffer_regions('end', 'y')

                do i = 1, offset_y%end
                    y_cb(n + i) = y_cb(n + (i - 1)) + dy(n + i)
                end do

                do i = 1, buff_size
                    y_cc(n + i) = y_cc(n + (i - 1)) + (dy(n + (i - 1)) + dy(n + i))/2d0
                end do

            end if

            ! END: Populating Buffer Regions in the y-direction ================

            ! Populating Buffer Regions in the z-direction =====================

            if (p > 0) then

                ! Ghost-cell extrapolation BC at the beginning
                if (bc_z%beg <= -3) then

                    do i = 1, buff_size
                        dz(-i) = dz(0)
                    end do

                    ! Symmetry BC at the beginning
                elseif (bc_z%beg == -2) then

                    do i = 1, buff_size
                        dz(-i) = dz(i - 1)
                    end do

                    ! Periodic BC at the beginning
                elseif (bc_z%beg == -1) then

                    do i = 1, buff_size
                        dz(-i) = dz((p + 1) - i)
                    end do

                    ! Processor BC at the beginning
                else

                    call s_mpi_sendrecv_grid_vars_buffer_regions('beg', 'z')

                    do i = 1, offset_z%beg
                        z_cb(-1 - i) = z_cb(-i) - dz(-i)
                    end do

                    do i = 1, buff_size
                        z_cc(-i) = z_cc(1 - i) - (dz(1 - i) + dz(-i))/2d0
                    end do

                end if

                ! Ghost-cell extrapolation BC at the end
                if (bc_z%end <= -3) then

                    do i = 1, buff_size
                        dz(p + i) = dz(p)
                    end do

                    ! Symmetry BC at the end
                elseif (bc_z%end == -2) then

                    do i = 1, buff_size
                        dz(p + i) = dz((p + 1) - i)
                    end do

                    ! Periodic BC at the end
                elseif (bc_z%end == -1) then

                    do i = 1, buff_size
                        dz(p + i) = dz(i - 1)
                    end do

                    ! Processor BC at the end
                else

                    call s_mpi_sendrecv_grid_vars_buffer_regions('end', 'z')

                    do i = 1, offset_z%end
                        z_cb(p + i) = z_cb(p + (i - 1)) + dz(p + i)
                    end do

                    do i = 1, buff_size
                        z_cc(p + i) = z_cc(p + (i - 1)) + (dz(p + (i - 1)) + dz(p + i))/2d0
                    end do

                end if

            end if

        end if

        ! END: Populating Buffer Regions in the z-direction ================

    end subroutine s_populate_grid_variables_buffer_regions ! --------------

    !>  The purpose of this procedure is to populate the buffers
        !!      of the cell-average conservative variables, depending on
        !!      the boundary conditions.
    subroutine s_populate_conservative_variables_buffer_regions() ! --------

        integer :: i, j, k !< Generic loop iterators

        ! Populating Buffer Regions in the x-direction =====================

        ! Ghost-cell extrapolation BC at the beginning
        if (bc_x%beg <= -3) then

            do j = 1, buff_size
                do i = 1, sys_size
                    q_cons_vf(i)%sf(-j, 0:n, 0:p) = q_cons_vf(i)%sf(0, 0:n, 0:p)
                end do
            end do

            ! Symmetry BC at the beginning
        elseif (bc_x%beg == -2) then

            do j = 1, buff_size

                ! Density or partial densities
                do i = 1, cont_idx%end
                    q_cons_vf(i)%sf(-j, 0:n, 0:p) = &
                        q_cons_vf(i)%sf(j - 1, 0:n, 0:p)
                end do

                ! x-component of momentum
                q_cons_vf(mom_idx%beg)%sf(-j, 0:n, 0:p) = &
                    -q_cons_vf(mom_idx%beg)%sf(j - 1, 0:n, 0:p)

                ! Remaining momentum component(s), if any, as well as the
                ! energy and the variable(s) from advection equation(s)
                do i = mom_idx%beg + 1, sys_size
                    q_cons_vf(i)%sf(-j, 0:n, 0:p) = &
                        q_cons_vf(i)%sf(j - 1, 0:n, 0:p)
                end do

            end do

            ! Periodic BC at the beginning
        elseif (bc_x%beg == -1) then

            do j = 1, buff_size
                do i = 1, sys_size
                    q_cons_vf(i)%sf(-j, 0:n, 0:p) = &
                        q_cons_vf(i)%sf((m + 1) - j, 0:n, 0:p)
                end do
            end do

            ! Processor BC at the beginning
        else

            call s_mpi_sendrecv_cons_vars_buffer_regions(q_cons_vf, &
                                                         'beg', 'x')

        end if

        ! Ghost-cell extrapolation BC at the end
        if (bc_x%end <= -3) then

            do j = 1, buff_size
                do i = 1, sys_size
                    q_cons_vf(i)%sf(m + j, 0:n, 0:p) = &
                        q_cons_vf(i)%sf(m, 0:n, 0:p)
                end do
            end do

            ! Symmetry BC at the end
        elseif (bc_x%end == -2) then

            do j = 1, buff_size

                ! Density or partial densities
                do i = 1, cont_idx%end
                    q_cons_vf(i)%sf(m + j, 0:n, 0:p) = &
                        q_cons_vf(i)%sf((m + 1) - j, 0:n, 0:p)
                end do

                ! x-component of momentum
                q_cons_vf(mom_idx%beg)%sf(m + j, 0:n, 0:p) = &
                    -q_cons_vf(mom_idx%beg)%sf((m + 1) - j, 0:n, 0:p)

                ! Remaining momentum component(s), if any, as well as the
                ! energy and the variable(s) from advection equation(s)
                do i = mom_idx%beg + 1, sys_size
                    q_cons_vf(i)%sf(m + j, 0:n, 0:p) = &
                        q_cons_vf(i)%sf((m + 1) - j, 0:n, 0:p)
                end do

            end do

            ! Perodic BC at the end
        elseif (bc_x%end == -1) then

            do j = 1, buff_size
                do i = 1, sys_size
                    q_cons_vf(i)%sf(m + j, 0:n, 0:p) = &
                        q_cons_vf(i)%sf(j - 1, 0:n, 0:p)
                end do
            end do

            ! Processor BC at the end
        else

            call s_mpi_sendrecv_cons_vars_buffer_regions(q_cons_vf, &
                                                         'end', 'x')

        end if

        ! END: Populating Buffer Regions in the x-direction ================

        ! Populating Buffer Regions in the y-direction =====================

        if (n > 0) then

            ! Ghost-cell extrapolation BC at the beginning
            if (bc_y%beg <= -3 .and. bc_y%beg /= -13) then

                do j = 1, buff_size
                    do i = 1, sys_size
                        q_cons_vf(i)%sf(:, -j, 0:p) = q_cons_vf(i)%sf(:, 0, 0:p)
                    end do
                end do

                ! Axis BC at the beginning
            elseif (bc_y%beg == -13) then

                do j = 1, buff_size
                    do k = 0, p
                        if (z_cc(k) < pi) then
                            do i = 1, mom_idx%beg
                                q_cons_vf(i)%sf(:, -j, k) = &
                                    q_cons_vf(i)%sf(:, j - 1, k + ((p + 1)/2))
                            end do

                            q_cons_vf(mom_idx%beg + 1)%sf(:, -j, k) = &
                                -q_cons_vf(mom_idx%beg + 1)%sf(:, j - 1, k + ((p + 1)/2))

                            q_cons_vf(mom_idx%end)%sf(:, -j, k) = &
                                -q_cons_vf(mom_idx%end)%sf(:, j - 1, k + ((p + 1)/2))

                            do i = E_idx, sys_size
                                q_cons_vf(i)%sf(:, -j, k) = &
                                    q_cons_vf(i)%sf(:, j - 1, k + ((p + 1)/2))
                            end do
                        else
                            do i = 1, mom_idx%beg
                                q_cons_vf(i)%sf(:, -j, k) = &
                                    q_cons_vf(i)%sf(:, j - 1, k - ((p + 1)/2))
                            end do

                            q_cons_vf(mom_idx%beg + 1)%sf(:, -j, k) = &
                                -q_cons_vf(mom_idx%beg + 1)%sf(:, j - 1, k - ((p + 1)/2))

                            q_cons_vf(mom_idx%end)%sf(:, -j, k) = &
                                -q_cons_vf(mom_idx%end)%sf(:, j - 1, k - ((p + 1)/2))

                            do i = E_idx, sys_size
                                q_cons_vf(i)%sf(:, -j, k) = &
                                    q_cons_vf(i)%sf(:, j - 1, k - ((p + 1)/2))
                            end do
                        end if
                    end do
                end do

                ! Symmetry BC at the beginning
            elseif (bc_y%beg == -2) then

                do j = 1, buff_size

                    ! Density or partial densities and x-momentum component
                    do i = 1, mom_idx%beg
                        q_cons_vf(i)%sf(:, -j, 0:p) = &
                            q_cons_vf(i)%sf(:, j - 1, 0:p)
                    end do

                    ! y-component of momentum
                    q_cons_vf(mom_idx%beg + 1)%sf(:, -j, 0:p) = &
                        -q_cons_vf(mom_idx%beg + 1)%sf(:, j - 1, 0:p)

                    ! Remaining z-momentum component, if any, as well as the
                    ! energy and variable(s) from advection equation(s)
                    do i = mom_idx%beg + 2, sys_size
                        q_cons_vf(i)%sf(:, -j, 0:p) = &
                            q_cons_vf(i)%sf(:, j - 1, 0:p)
                    end do

                end do

                ! Periodic BC at the beginning
            elseif (bc_y%beg == -1) then

                do j = 1, buff_size
                    do i = 1, sys_size
                        q_cons_vf(i)%sf(:, -j, 0:p) = &
                            q_cons_vf(i)%sf(:, (n + 1) - j, 0:p)
                    end do
                end do

                ! Processor BC at the beginning
            else

                call s_mpi_sendrecv_cons_vars_buffer_regions(q_cons_vf, &
                                                             'beg', 'y')

            end if

            ! Ghost-cell extrapolation BC at the end
            if (bc_y%end <= -3) then

                do j = 1, buff_size
                    do i = 1, sys_size
                        q_cons_vf(i)%sf(:, n + j, 0:p) = &
                            q_cons_vf(i)%sf(:, n, 0:p)
                    end do
                end do

                ! Symmetry BC at the end
            elseif (bc_y%end == -2) then

                do j = 1, buff_size

                    ! Density or partial densities and x-momentum component
                    do i = 1, mom_idx%beg
                        q_cons_vf(i)%sf(:, n + j, 0:p) = &
                            q_cons_vf(i)%sf(:, (n + 1) - j, 0:p)
                    end do

                    ! y-component of momentum
                    q_cons_vf(mom_idx%beg + 1)%sf(:, n + j, 0:p) = &
                        -q_cons_vf(mom_idx%beg + 1)%sf(:, (n + 1) - j, 0:p)

                    ! Remaining z-momentum component, if any, as well as the
                    ! energy and variable(s) from advection equation(s)
                    do i = mom_idx%beg + 2, sys_size
                        q_cons_vf(i)%sf(:, n + j, 0:p) = &
                            q_cons_vf(i)%sf(:, (n + 1) - j, 0:p)
                    end do

                end do

                ! Perodic BC at the end
            elseif (bc_y%end == -1) then

                do j = 1, buff_size
                    do i = 1, sys_size
                        q_cons_vf(i)%sf(:, n + j, 0:p) = &
                            q_cons_vf(i)%sf(:, j - 1, 0:p)
                    end do
                end do

                ! Processor BC at the end
            else

                call s_mpi_sendrecv_cons_vars_buffer_regions(q_cons_vf, &
                                                             'end', 'y')

            end if

            ! END: Populating Buffer Regions in the y-direction ================

            ! Populating Buffer Regions in the z-direction =====================

            if (p > 0) then

                ! Ghost-cell extrapolation BC at the beginning
                if (bc_z%beg <= -3) then

                    do j = 1, buff_size
                        do i = 1, sys_size
                            q_cons_vf(i)%sf(:, :, -j) = q_cons_vf(i)%sf(:, :, 0)
                        end do
                    end do

                    ! Symmetry BC at the beginning
                elseif (bc_z%beg == -2) then

                    do j = 1, buff_size

                        ! Density or the partial densities and the momentum
                        ! components in x- and y-directions
                        do i = 1, mom_idx%beg + 1
                            q_cons_vf(i)%sf(:, :, -j) = &
                                q_cons_vf(i)%sf(:, :, j - 1)
                        end do

                        ! z-component of momentum
                        q_cons_vf(mom_idx%end)%sf(:, :, -j) = &
                            -q_cons_vf(mom_idx%end)%sf(:, :, j - 1)

                        ! Energy and advection equation(s) variable(s)
                        do i = E_idx, sys_size
                            q_cons_vf(i)%sf(:, :, -j) = &
                                q_cons_vf(i)%sf(:, :, j - 1)
                        end do

                    end do

                    ! Periodic BC at the beginning
                elseif (bc_z%beg == -1) then

                    do j = 1, buff_size
                        do i = 1, sys_size
                            q_cons_vf(i)%sf(:, :, -j) = &
                                q_cons_vf(i)%sf(:, :, (p + 1) - j)
                        end do
                    end do

                    ! Processor BC at the beginning
                else

                    call s_mpi_sendrecv_cons_vars_buffer_regions(q_cons_vf, &
                                                                 'beg', 'z')

                end if

                ! Ghost-cell extrapolation BC at the end
                if (bc_z%end <= -3) then

                    do j = 1, buff_size
                        do i = 1, sys_size
                            q_cons_vf(i)%sf(:, :, p + j) = &
                                q_cons_vf(i)%sf(:, :, p)
                        end do
                    end do

                    ! Symmetry BC at the end
                elseif (bc_z%end == -2) then

                    do j = 1, buff_size

                        ! Density or the partial densities and the momentum
                        ! components in x- and y-directions
                        do i = 1, mom_idx%beg + 1
                            q_cons_vf(i)%sf(:, :, p + j) = &
                                q_cons_vf(i)%sf(:, :, (p + 1) - j)
                        end do

                        ! z-component of momentum
                        q_cons_vf(mom_idx%end)%sf(:, :, p + j) = &
                            -q_cons_vf(mom_idx%end)%sf(:, :, (p + 1) - j)

                        ! Energy and advection equation(s) variable(s)
                        do i = E_idx, sys_size
                            q_cons_vf(i)%sf(:, :, p + j) = &
                                q_cons_vf(i)%sf(:, :, (p + 1) - j)
                        end do

                    end do

                    ! Perodic BC at the end
                elseif (bc_z%end == -1) then

                    do j = 1, buff_size
                        do i = 1, sys_size
                            q_cons_vf(i)%sf(:, :, p + j) = &
                                q_cons_vf(i)%sf(:, :, j - 1)
                        end do
                    end do

                    ! Processor BC at the end
                else

                    call s_mpi_sendrecv_cons_vars_buffer_regions(q_cons_vf, &
                                                                 'end', 'z')

                end if

            end if

        end if

        ! END: Populating Buffer Regions in the z-direction ================

    end subroutine s_populate_conservative_variables_buffer_regions ! ------

    !>  Computation of parameters, allocation procedures, and/or
        !!      any other tasks needed to properly setup the module
    subroutine s_initialize_data_input_module() ! -----------------------------

        integer :: i !< Generic loop iterator

        ! Allocating the parts of the conservative and primitive variables
        ! that do not require the direct knowledge of the dimensionality of
        ! the simulation
        allocate (q_cons_vf(1:sys_size))
        allocate (q_prim_vf(1:sys_size))

        ! Allocating the parts of the conservative and primitive variables
        ! that do require the direct knowledge of the dimensionality of the
        ! simulation

        ! Simulation is at least 2D
        if (n > 0) then

            ! Simulation is 3D
            if (p > 0) then

                do i = 1, sys_size
                    allocate (q_cons_vf(i)%sf(-buff_size:m + buff_size, &
                                              -buff_size:n + buff_size, &
                                              -buff_size:p + buff_size))
                    allocate (q_prim_vf(i)%sf(-buff_size:m + buff_size, &
                                              -buff_size:n + buff_size, &
                                              -buff_size:p + buff_size))
                end do

                ! Simulation is 2D
            else

                do i = 1, sys_size
                    allocate (q_cons_vf(i)%sf(-buff_size:m + buff_size, &
                                              -buff_size:n + buff_size, &
                                              0:0))
                    allocate (q_prim_vf(i)%sf(-buff_size:m + buff_size, &
                                              -buff_size:n + buff_size, &
                                              0:0))
                end do

            end if

            ! Simulation is 1D
        else

            do i = 1, sys_size
                allocate (q_cons_vf(i)%sf(-buff_size:m + buff_size, &
                                          0:0, &
                                          0:0))
                allocate (q_prim_vf(i)%sf(-buff_size:m + buff_size, &
                                          0:0, &
                                          0:0))
            end do

        end if

        if (parallel_io .neqv. .true.) then
            s_read_data_files => s_read_serial_data_files
        else
            s_read_data_files => s_read_parallel_data_files
        end if

    end subroutine s_initialize_data_input_module ! ---------------------------

    !> Deallocation procedures for the module
    subroutine s_finalize_data_input_module() ! --------------------------

        integer :: i !< Generic loop iterator

        ! Deallocating the conservative and primitive variables
        do i = 1, sys_size
            deallocate (q_cons_vf(i)%sf)
            deallocate (q_prim_vf(i)%sf)
        end do

        deallocate (q_cons_vf)
        deallocate (q_prim_vf)

        s_read_data_files => null()

    end subroutine s_finalize_data_input_module ! ------------------------

end module m_data_input
