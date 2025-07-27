!>
!! @file m_data_input.f90
!> @brief Contains module m_data_input

!> @brief This module features procedures, which for a specific time-step,
!!             read in the raw simulation data for the grid and the conservative
!!             variables and fill out their buffer regions.
module m_data_input

#ifdef MFC_MPI
    use mpi                     !< Message passing interface (MPI) module
#endif

    use m_derived_types         !< Definitions of the derived types

    use m_global_parameters     !< Global parameters for the code

    use m_mpi_proxy             !< Message passing interface (MPI) module proxy

    use m_mpi_common

    use m_compile_specific

    use m_boundary_common

    use m_helper

    implicit none

    private; public :: s_initialize_data_input_module, &
 s_read_data_files, &
 s_read_serial_data_files, &
 s_read_parallel_data_files, &
 s_finalize_data_input_module

    abstract interface

        !> Subroutine for reading data files
        !!  @param t_step Current time-step to input
        impure subroutine s_read_abstract_data_files(t_step)

            implicit none

            integer, intent(in) :: t_step

        end subroutine s_read_abstract_data_files

    end interface

    type(scalar_field), allocatable, dimension(:), public :: q_cons_vf !<
    !! Conservative variables

    type(scalar_field), allocatable, dimension(:), public :: q_prim_vf !<
    !! Primitive variables

    type(integer_field), allocatable, dimension(:, :), public :: bc_type !<
    !! Boundary condition identifiers

    type(scalar_field), public :: q_T_sf !<
    !! Temperature field

    ! type(scalar_field), public :: ib_markers !<
    type(integer_field), public :: ib_markers

    procedure(s_read_abstract_data_files), pointer :: s_read_data_files => null()

contains

    !> Helper subroutine to read grid data files for a given direction
    !!  @param t_step_dir Directory containing the time-step data
    !!  @param direction Direction name ('x', 'y', 'z')
    !!  @param cb_array Cell boundary array to populate
    !!  @param d_array Cell width array to populate
    !!  @param cc_array Cell center array to populate
    !!  @param size_dim Size of the dimension
    impure subroutine s_read_grid_data_direction(t_step_dir, direction, cb_array, d_array, cc_array, size_dim)

        character(len=*), intent(in) :: t_step_dir
        character(len=1), intent(in) :: direction
        real(wp), dimension(-1:), intent(out) :: cb_array
        real(wp), dimension(0:), intent(out) :: d_array
        real(wp), dimension(0:), intent(out) :: cc_array
        integer, intent(in) :: size_dim

        character(LEN=len_trim(t_step_dir) + 10) :: file_loc
        logical :: file_check

        ! Checking whether direction_cb.dat exists
        file_loc = trim(t_step_dir)//'/'//direction//'_cb.dat'
        inquire (FILE=trim(file_loc), EXIST=file_check)

        ! Reading direction_cb.dat if it exists, exiting otherwise
        if (file_check) then
            open (1, FILE=trim(file_loc), FORM='unformatted', &
                  STATUS='old', ACTION='read')
            read (1) cb_array(-1:size_dim)
            close (1)
        else
            call s_mpi_abort('File '//direction//'_cb.dat is missing in '// &
                             trim(t_step_dir)//'. Exiting.')
        end if

        ! Computing the cell-width distribution
        d_array(0:size_dim) = cb_array(0:size_dim) - cb_array(-1:size_dim - 1)

        ! Computing the cell-center locations
        cc_array(0:size_dim) = cb_array(-1:size_dim - 1) + d_array(0:size_dim)/2._wp

    end subroutine s_read_grid_data_direction

#ifdef MFC_MPI
    !> Helper subroutine to setup MPI data I/O parameters
    !!  @param data_size Local array size (output)
    !!  @param m_MOK, n_MOK, p_MOK MPI offset kinds for dimensions (output)
    !!  @param WP_MOK, MOK, str_MOK, NVARS_MOK Other MPI offset kinds (output)
    impure subroutine s_setup_mpi_io_params(data_size, m_MOK, n_MOK, p_MOK, WP_MOK, MOK, str_MOK, NVARS_MOK)

        integer, intent(out) :: data_size
        integer(KIND=MPI_OFFSET_KIND), intent(out) :: m_MOK, n_MOK, p_MOK
        integer(KIND=MPI_OFFSET_KIND), intent(out) :: WP_MOK, MOK, str_MOK, NVARS_MOK

        ! Initialize MPI data I/O
        if (ib) then
            call s_initialize_mpi_data(q_cons_vf, ib_markers)
        else
            call s_initialize_mpi_data(q_cons_vf)
        end if

        ! Size of local arrays
        data_size = (m + 1)*(n + 1)*(p + 1)

        ! Resize some integers so MPI can read even the biggest file
        m_MOK = int(m_glb + 1, MPI_OFFSET_KIND)
        n_MOK = int(n_glb + 1, MPI_OFFSET_KIND)
        p_MOK = int(p_glb + 1, MPI_OFFSET_KIND)
        WP_MOK = int(8._wp, MPI_OFFSET_KIND)
        MOK = int(1._wp, MPI_OFFSET_KIND)
        str_MOK = int(name_len, MPI_OFFSET_KIND)
        NVARS_MOK = int(sys_size, MPI_OFFSET_KIND)

    end subroutine s_setup_mpi_io_params
#endif

    !> Helper subroutine to read IB data files
    !!  @param file_loc_base Base file location for IB data
    impure subroutine s_read_ib_data_files(file_loc_base)

        character(len=*), intent(in) :: file_loc_base

        character(LEN=len_trim(file_loc_base) + 20) :: file_loc
        logical :: file_exist
        integer :: ifile, ierr, data_size

#ifdef MFC_MPI
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(KIND=MPI_OFFSET_KIND) :: disp
#endif
        if (.not. ib) return

        if (parallel_io) then
            write (file_loc, '(A)') trim(file_loc_base)//'ib.dat'
        else
            write (file_loc, '(A)') trim(file_loc_base)//'/ib.dat'
        end if
        inquire (FILE=trim(file_loc), EXIST=file_exist)

        if (file_exist) then
            if (parallel_io) then
#ifdef MFC_MPI
                call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

                data_size = (m + 1)*(n + 1)*(p + 1)
                disp = 0

                call MPI_FILE_SET_VIEW(ifile, disp, MPI_INTEGER, MPI_IO_IB_DATA%view, &
                                       'native', mpi_info_int, ierr)
                call MPI_FILE_READ(ifile, MPI_IO_IB_DATA%var%sf, data_size, &
                                   MPI_INTEGER, status, ierr)

                call MPI_FILE_CLOSE(ifile, ierr)
#endif
            else
                open (2, FILE=trim(file_loc), &
                      FORM='unformatted', &
                      ACTION='read', &
                      STATUS='old')
                read (2) ib_markers%sf(0:m, 0:n, 0:p)
                close (2)
            end if
        else
            call s_mpi_abort('File '//trim(file_loc)//' is missing. Exiting.')
        end if

    end subroutine s_read_ib_data_files

    !> Helper subroutine to allocate field arrays for given dimensionality
    !!  @param start_idx Starting index for allocation
    !!  @param end_x, end_y, end_z End indices for each dimension
    impure subroutine s_allocate_field_arrays(local_start_idx, end_x, end_y, end_z)

        integer, intent(in) :: local_start_idx, end_x, end_y, end_z
        integer :: i

        do i = 1, sys_size
            allocate (q_cons_vf(i)%sf(local_start_idx:end_x, local_start_idx:end_y, local_start_idx:end_z))
            allocate (q_prim_vf(i)%sf(local_start_idx:end_x, local_start_idx:end_y, local_start_idx:end_z))
        end do

        if (ib) then
            allocate (ib_markers%sf(local_start_idx:end_x, local_start_idx:end_y, local_start_idx:end_z))
        end if

        if (chemistry) then
            allocate (q_T_sf%sf(local_start_idx:end_x, local_start_idx:end_y, local_start_idx:end_z))
        end if

    end subroutine s_allocate_field_arrays

    !>  This subroutine is called at each time-step that has to
        !!      be post-processed in order to read the raw data files
        !!      present in the corresponding time-step directory and to
        !!      populate the associated grid and conservative variables.
        !!  @param t_step Current time-step
    impure subroutine s_read_serial_data_files(t_step)

        integer, intent(in) :: t_step

        character(LEN=len_trim(case_dir) + 2*name_len) :: t_step_dir !<
            !! Location of the time-step directory associated with t_step

        character(LEN=len_trim(case_dir) + 3*name_len) :: file_loc !<
            !! Generic string used to store the location of a particular file

        character(LEN= &
                  int(floor(log10(real(sys_size, wp)))) + 1) :: file_num !<
            !! Used to store the variable position, in character form, of the
            !! currently manipulated conservative variable file

        character(LEN=len_trim(case_dir) + 2*name_len) :: t_step_ib_dir !<
        !! Location of the time-step directory associated with t_step

        logical :: dir_check !<
            !! Generic logical used to test the existence of a particular folder

        logical :: file_check  !<
            !! Generic logical used to test the existence of a particular file

        integer :: i !< Generic loop iterator

        ! Setting location of time-step folder based on current time-step
        write (t_step_dir, '(A,I0,A,I0)') '/p_all/p', proc_rank, '/', t_step
        t_step_dir = trim(case_dir)//trim(t_step_dir)

        write (t_step_ib_dir, '(A,I0,A,I0)') '/p_all/p', proc_rank, '/', 0
        t_step_ib_dir = trim(case_dir)//trim(t_step_ib_dir)

        ! Inquiring as to the existence of the time-step directory
        file_loc = trim(t_step_dir)//'/.'
        call my_inquire(file_loc, dir_check)

        ! If the time-step directory is missing, the post-process exits.
        if (dir_check .neqv. .true.) then
            call s_mpi_abort('Time-step folder '//trim(t_step_dir)// &
                             ' is missing. Exiting.')
        end if

        file_loc = trim(t_step_ib_dir)//'/.'
        call my_inquire(file_loc, dir_check)

        ! If the time-step directory is missing, the post-process exits.
        if (dir_check .neqv. .true.) then
            call s_mpi_abort('Time-step folder '//trim(t_step_ib_dir)// &
                             ' is missing. Exiting.')
        end if

        if (bc_io) then
            call s_read_serial_boundary_condition_files(t_step_dir, bc_type)
        else
            call s_assign_default_bc_type(bc_type)
        end if

        ! Reading the Grid Data Files using helper subroutine
        call s_read_grid_data_direction(t_step_dir, 'x', x_cb, dx, x_cc, m)

        if (n > 0) then
            call s_read_grid_data_direction(t_step_dir, 'y', y_cb, dy, y_cc, n)

            if (p > 0) then
                call s_read_grid_data_direction(t_step_dir, 'z', z_cb, dz, z_cc, p)
            end if
        end if

        ! Reading the Conservative Variables Data Files
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
                call s_mpi_abort('File q_cons_vf'//trim(file_num)// &
                                 '.dat is missing in '//trim(t_step_dir)// &
                                 '. Exiting.')
            end if

        end do

        ! Reading IB data using helper subroutine
        call s_read_ib_data_files(t_step_ib_dir)

    end subroutine s_read_serial_data_files

    !>  This subroutine is called at each time-step that has to
        !!      be post-processed in order to parallel-read the raw data files
        !!      present in the corresponding time-step directory and to
        !!      populate the associated grid and conservative variables.
        !!  @param t_step Current time-step
    impure subroutine s_read_parallel_data_files(t_step)

        integer, intent(in) :: t_step

#ifdef MFC_MPI

        real(wp), allocatable, dimension(:) :: x_cb_glb, y_cb_glb, z_cb_glb

        integer :: ifile, ierr, data_size
        integer, dimension(MPI_STATUS_SIZE) :: status

        integer(KIND=MPI_OFFSET_KIND) :: disp
        integer(KIND=MPI_OFFSET_KIND) :: m_MOK, n_MOK, p_MOK
        integer(KIND=MPI_OFFSET_KIND) :: WP_MOK, var_MOK, str_MOK
        integer(KIND=MPI_OFFSET_KIND) :: NVARS_MOK
        integer(KIND=MPI_OFFSET_KIND) :: MOK

        character(LEN=path_len + 2*name_len) :: file_loc
        logical :: file_exist

        character(len=10) :: t_step_string

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
            call MPI_FILE_READ(ifile, x_cb_glb, data_size, mpi_p, status, ierr)
            call MPI_FILE_CLOSE(ifile, ierr)
        else
            call s_mpi_abort('File '//trim(file_loc)//' is missing. Exiting.')
        end if

        ! Assigning local cell boundary locations
        x_cb(-1:m) = x_cb_glb((start_idx(1) - 1):(start_idx(1) + m))
        ! Computing the cell width distribution
        dx(0:m) = x_cb(0:m) - x_cb(-1:m - 1)
        ! Computing the cell center location
        x_cc(0:m) = x_cb(-1:m - 1) + dx(0:m)/2._wp

        if (n > 0) then
            ! Read in cell boundary locations in y-direction
            file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//'y_cb.dat'
            inquire (FILE=trim(file_loc), EXIST=file_exist)

            if (file_exist) then
                data_size = n_glb + 2
                call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
                call MPI_FILE_READ(ifile, y_cb_glb, data_size, mpi_p, status, ierr)
                call MPI_FILE_CLOSE(ifile, ierr)
            else
                call s_mpi_abort('File '//trim(file_loc)//' is missing. Exiting.')
            end if

            ! Assigning local cell boundary locations
            y_cb(-1:n) = y_cb_glb((start_idx(2) - 1):(start_idx(2) + n))
            ! Computing the cell width distribution
            dy(0:n) = y_cb(0:n) - y_cb(-1:n - 1)
            ! Computing the cell center location
            y_cc(0:n) = y_cb(-1:n - 1) + dy(0:n)/2._wp

            if (p > 0) then
                ! Read in cell boundary locations in z-direction
                file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//'z_cb.dat'
                inquire (FILE=trim(file_loc), EXIST=file_exist)

                if (file_exist) then
                    data_size = p_glb + 2
                    call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
                    call MPI_FILE_READ(ifile, z_cb_glb, data_size, mpi_p, status, ierr)
                    call MPI_FILE_CLOSE(ifile, ierr)
                else
                    call s_mpi_abort('File '//trim(file_loc)//' is missing. Exiting.')
                end if

                ! Assigning local cell boundary locations
                z_cb(-1:p) = z_cb_glb((start_idx(3) - 1):(start_idx(3) + p))
                ! Computing the cell width distribution
                dz(0:p) = z_cb(0:p) - z_cb(-1:p - 1)
                ! Computing the cell center location
                z_cc(0:p) = z_cb(-1:p - 1) + dz(0:p)/2._wp
            end if
        end if

        call s_read_parallel_conservative_data(t_step, m_MOK, n_MOK, p_MOK, WP_MOK, MOK, str_MOK, NVARS_MOK)

        deallocate (x_cb_glb, y_cb_glb, z_cb_glb)

        if (bc_io) then
            call s_read_parallel_boundary_condition_files(bc_type)
        else
            call s_assign_default_bc_type(bc_type)
        end if

#endif

    end subroutine s_read_parallel_data_files

#ifdef MFC_MPI
    !> Helper subroutine to read parallel conservative variable data
    !!  @param t_step Current time-step
    !!  @param m_MOK, n_MOK, p_MOK MPI offset kinds for dimensions
    !!  @param WP_MOK, MOK, str_MOK, NVARS_MOK Other MPI offset kinds
    impure subroutine s_read_parallel_conservative_data(t_step, m_MOK, n_MOK, p_MOK, WP_MOK, MOK, str_MOK, NVARS_MOK)

        integer, intent(in) :: t_step
        integer(KIND=MPI_OFFSET_KIND), intent(inout) :: m_MOK, n_MOK, p_MOK
        integer(KIND=MPI_OFFSET_KIND), intent(inout) :: WP_MOK, MOK, str_MOK, NVARS_MOK

        integer :: ifile, ierr, data_size
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(KIND=MPI_OFFSET_KIND) :: disp, var_MOK
        character(LEN=path_len + 2*name_len) :: file_loc
        logical :: file_exist
        character(len=10) :: t_step_string
        integer :: i

        if (file_per_process) then
            call s_int_to_str(t_step, t_step_string)
            ! Open the file to read conservative variables
            write (file_loc, '(I0,A1,I7.7,A)') t_step, '_', proc_rank, '.dat'
            file_loc = trim(case_dir)//'/restart_data/lustre_'//trim(t_step_string)//trim(mpiiofs)//trim(file_loc)
            inquire (FILE=trim(file_loc), EXIST=file_exist)

            if (file_exist) then
                call MPI_FILE_OPEN(MPI_COMM_SELF, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

                call s_setup_mpi_io_params(data_size, m_MOK, n_MOK, p_MOK, WP_MOK, MOK, str_MOK, NVARS_MOK)

                ! Read the data for each variable
                if (bubbles_euler .or. elasticity .or. mhd) then
                    do i = 1, sys_size
                        var_MOK = int(i, MPI_OFFSET_KIND)
                        call MPI_FILE_READ_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                               mpi_p, status, ierr)
                    end do
                else
                    do i = 1, adv_idx%end
                        var_MOK = int(i, MPI_OFFSET_KIND)
                        call MPI_FILE_READ_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                               mpi_p, status, ierr)
                    end do
                end if

                call s_mpi_barrier()
                call MPI_FILE_CLOSE(ifile, ierr)

                call s_read_ib_data_files(trim(case_dir)//'/restart_data'//trim(mpiiofs))
            else
                call s_mpi_abort('File '//trim(file_loc)//' is missing. Exiting.')
            end if
        else
            ! Open the file to read conservative variables
            write (file_loc, '(I0,A)') t_step, '.dat'
            file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//trim(file_loc)
            inquire (FILE=trim(file_loc), EXIST=file_exist)

            if (file_exist) then
                call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

                call s_setup_mpi_io_params(data_size, m_MOK, n_MOK, p_MOK, WP_MOK, MOK, str_MOK, NVARS_MOK)

                ! Read the data for each variable
                do i = 1, sys_size
                    var_MOK = int(i, MPI_OFFSET_KIND)

                    ! Initial displacement to skip at beginning of file
                    disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                    call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, MPI_IO_DATA%view(i), &
                                           'native', mpi_info_int, ierr)
                    call MPI_FILE_READ_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size, &
                                           mpi_p, status, ierr)
                end do

                call s_mpi_barrier()
                call MPI_FILE_CLOSE(ifile, ierr)

                call s_read_ib_data_files(trim(case_dir)//'/restart_data'//trim(mpiiofs))
            else
                call s_mpi_abort('File '//trim(file_loc)//' is missing. Exiting.')
            end if
        end if
    end subroutine s_read_parallel_conservative_data
#endif

    !>  Computation of parameters, allocation procedures, and/or
        !!      any other tasks needed to properly setup the module
    impure subroutine s_initialize_data_input_module

        integer :: i !< Generic loop iterator

        ! Allocating the parts of the conservative and primitive variables
        ! that do not require the direct knowledge of the dimensionality of
        ! the simulation
        allocate (q_cons_vf(1:sys_size))
        allocate (q_prim_vf(1:sys_size))

        ! Allocating the parts of the conservative and primitive variables
        ! that do require the direct knowledge of the dimensionality of
        ! the simulation using helper subroutine

        ! Simulation is at least 2D
        if (n > 0) then
            ! Simulation is 3D
            if (p > 0) then
                call s_allocate_field_arrays(-buff_size, m + buff_size, n + buff_size, p + buff_size)
            else
                ! Simulation is 2D
                call s_allocate_field_arrays(-buff_size, m + buff_size, n + buff_size, 0)
            end if
        else
            ! Simulation is 1D
            call s_allocate_field_arrays(-buff_size, m + buff_size, 0, 0)
        end if

        ! Allocating arrays to store the bc types
        allocate (bc_type(1:num_dims, -1:1))

        allocate (bc_type(1, -1)%sf(0:0, 0:n, 0:p))
        allocate (bc_type(1, 1)%sf(0:0, 0:n, 0:p))
        if (n > 0) then
            allocate (bc_type(2, -1)%sf(-buff_size:m + buff_size, 0:0, 0:p))
            allocate (bc_type(2, 1)%sf(-buff_size:m + buff_size, 0:0, 0:p))
            if (p > 0) then
                allocate (bc_type(3, -1)%sf(-buff_size:m + buff_size, -buff_size:n + buff_size, 0:0))
                allocate (bc_type(3, 1)%sf(-buff_size:m + buff_size, -buff_size:n + buff_size, 0:0))
            end if
        end if

        if (parallel_io .neqv. .true.) then
            s_read_data_files => s_read_serial_data_files
        else
            s_read_data_files => s_read_parallel_data_files
        end if

    end subroutine s_initialize_data_input_module

    !> Deallocation procedures for the module
    impure subroutine s_finalize_data_input_module

        integer :: i !< Generic loop iterator

        ! Deallocating the conservative and primitive variables
        do i = 1, sys_size
            deallocate (q_cons_vf(i)%sf)
            deallocate (q_prim_vf(i)%sf)
        end do

        deallocate (q_cons_vf)
        deallocate (q_prim_vf)

        if (ib) then
            deallocate (ib_markers%sf)
        end if

        if (chemistry) then
            deallocate (q_T_sf%sf)
        end if

        deallocate (bc_type(1, -1)%sf, bc_type(1, 1)%sf)
        if (n > 0) then
            deallocate (bc_type(2, -1)%sf, bc_type(2, 1)%sf)
            if (p > 0) then
                deallocate (bc_type(3, -1)%sf, bc_type(3, 1)%sf)
            end if
        end if

        deallocate (bc_type)

        s_read_data_files => null()

    end subroutine s_finalize_data_input_module

end module m_data_input
