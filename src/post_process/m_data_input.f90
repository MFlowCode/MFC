!>
!! @file
!> @brief Contains module m_data_input

!> @brief Reads raw simulation grid and conservative-variable data for a given time-step and fills buffer regions
module m_data_input

#ifdef MFC_MPI
    use mpi
#endif

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_mpi_common
    use m_compile_specific
    use m_boundary_common
    use m_helper

    implicit none

    private; public :: s_initialize_data_input_module, s_read_data_files, s_read_serial_data_files, s_read_parallel_data_files, &
        & s_read_amr_data, s_free_amr_data, s_finalize_data_input_module

    abstract interface

        !> Subroutine for reading data files
        impure subroutine s_read_abstract_data_files(t_step)

            implicit none

            integer, intent(in) :: t_step

        end subroutine s_read_abstract_data_files
    end interface

    type(scalar_field), allocatable, dimension(:), public    :: q_cons_vf  !< Conservative variables
    type(scalar_field), allocatable, dimension(:), public    :: q_cons_temp
    type(scalar_field), allocatable, dimension(:), public    :: q_prim_vf  !< Primitive variables
    type(integer_field), allocatable, dimension(:,:), public :: bc_type    !< Boundary condition identifiers
    type(scalar_field), public                               :: q_T_sf     !< Temperature field
    ! type(scalar_field), public :: ib_markers !<
    type(integer_field), public :: ib_markers

    !> One AMR fine-block piece owned by this rank, held for visualization overlay of the refined solution.
    type, public :: amr_fine_block
        integer                                       :: lo(3), hi(3)      !< global coarse-index region bounds of the parent block
        integer                                       :: m, n, p           !< local fine extents (interior 0:m, 0:n, 0:p)
        real(wp), allocatable, dimension(:)           :: x_cb, y_cb, z_cb  !< reconstructed fine cell boundaries
        type(scalar_field), allocatable, dimension(:) :: q_cons            !< fine conservative state
    end type amr_fine_block

    type(amr_fine_block), allocatable, dimension(:), public :: amr_fine  !< this rank's owned block pieces
    integer, public :: amr_num_fine  !< number of block pieces this rank owns (<= file's block count)

    procedure(s_read_abstract_data_files), pointer :: s_read_data_files => null()

contains

    !> Helper subroutine to read grid data files for a given direction
    impure subroutine s_read_grid_data_direction(t_step_dir, direction, cb_array, d_array, cc_array, size_dim)

        character(len=*), intent(in)             :: t_step_dir
        character(len=1), intent(in)             :: direction
        real(wp), dimension(-1:), intent(out)    :: cb_array
        real(wp), dimension(0:), intent(out)     :: d_array
        real(wp), dimension(0:), intent(out)     :: cc_array
        integer, intent(in)                      :: size_dim
        character(LEN=len_trim(t_step_dir) + 10) :: file_loc
        logical                                  :: file_check

        file_loc = trim(t_step_dir) // '/' // direction // '_cb.dat'
        inquire (FILE=trim(file_loc), EXIST=file_check)

        if (file_check) then
            open (1, FILE=trim(file_loc), form='unformatted', STATUS='old', ACTION='read')
            read (1) cb_array(-1:size_dim)
            close (1)
        else
            call s_mpi_abort('File ' // direction // '_cb.dat is missing in ' // trim(t_step_dir) // '. Exiting.')
        end if

        d_array(0:size_dim) = cb_array(0:size_dim) - cb_array(-1:size_dim - 1)
        cc_array(0:size_dim) = cb_array(-1:size_dim - 1) + d_array(0:size_dim)/2._wp

    end subroutine s_read_grid_data_direction

#ifdef MFC_MPI
    !> Helper subroutine to setup MPI data I/O parameters
    impure subroutine s_setup_mpi_io_params(data_size, m_MOK, n_MOK, p_MOK, WP_MOK, MOK, str_MOK, NVARS_MOK)

        integer, intent(out)                       :: data_size
        integer(KIND=MPI_OFFSET_KIND), intent(out) :: m_MOK, n_MOK, p_MOK
        integer(KIND=MPI_OFFSET_KIND), intent(out) :: WP_MOK, MOK, str_MOK, NVARS_MOK

        if (ib) then
            call s_initialize_mpi_data(q_cons_vf, ib_markers)
        else
            call s_initialize_mpi_data(q_cons_vf)
        end if

        data_size = (m + 1)*(n + 1)*(p + 1)

        m_MOK = int(m_glb + 1, MPI_OFFSET_KIND)
        n_MOK = int(n_glb + 1, MPI_OFFSET_KIND)
        p_MOK = int(p_glb + 1, MPI_OFFSET_KIND)
        WP_MOK = int(storage_size(0._stp)/8, MPI_OFFSET_KIND)
        MOK = int(1._wp, MPI_OFFSET_KIND)
        str_MOK = int(name_len, MPI_OFFSET_KIND)
        NVARS_MOK = int(sys_size, MPI_OFFSET_KIND)

    end subroutine s_setup_mpi_io_params
#endif

    !> Helper subroutine to read IB data files
    impure subroutine s_read_ib_data_files(file_loc_base, t_step)

        character(len=*), intent(in)                :: file_loc_base
        integer, intent(in), optional               :: t_step
        character(LEN=len_trim(file_loc_base) + 20) :: file_loc
        logical                                     :: file_exist
        integer                                     :: ifile, ierr, data_size

#ifdef MFC_MPI
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(KIND=MPI_OFFSET_KIND)       :: disp
        integer(KIND=MPI_OFFSET_KIND)       :: m_MOK, n_MOK, p_MOK, MOK, WP_MOK, var_MOK
        integer                             :: save_index
#endif

        if (.not. ib) return

        if (parallel_io) then
            write (file_loc, '(A)') trim(file_loc_base) // 'ib.dat'
        else
            write (file_loc, '(A)') trim(file_loc_base) // '/ib_data.dat'
        end if
        inquire (FILE=trim(file_loc), EXIST=file_exist)

        if (file_exist) then
            if (parallel_io) then
#ifdef MFC_MPI
                call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

                m_MOK = int(m_glb + 1, MPI_OFFSET_KIND)
                n_MOK = int(n_glb + 1, MPI_OFFSET_KIND)
                p_MOK = int(p_glb + 1, MPI_OFFSET_KIND)
                MOK = int(1._wp, MPI_OFFSET_KIND)
                WP_MOK = int(storage_size(0._stp)/8, MPI_OFFSET_KIND)
                save_index = t_step/t_step_save  ! get the number of saves done to this point

                data_size = (m + 1)*(n + 1)*(p + 1)
                var_MOK = int(sys_size + 1, MPI_OFFSET_KIND)
                if (t_step == 0) then
                    disp = 0
                else
                    disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1 + int(save_index, MPI_OFFSET_KIND))
                end if

                call MPI_FILE_SET_VIEW(ifile, disp, MPI_INTEGER, MPI_IO_IB_DATA%view, 'native', mpi_info_int, ierr)
                call MPI_FILE_READ(ifile, MPI_IO_IB_DATA%var%sf, data_size, MPI_INTEGER, status, ierr)

                call MPI_FILE_CLOSE(ifile, ierr)
#endif
            else
                open (2, FILE=trim(file_loc), form='unformatted', ACTION='read', STATUS='old')
                read (2) ib_markers%sf(0:m,0:n,0:p)
                close (2)
            end if
        else
            call s_mpi_abort('File ' // trim(file_loc) // ' is missing. Exiting.')
        end if

    end subroutine s_read_ib_data_files

    !> Helper subroutine to allocate field arrays for given dimensionality
    impure subroutine s_allocate_field_arrays(local_start_idx, end_x, end_y, end_z)

        integer, intent(in) :: local_start_idx, end_x, end_y, end_z
        integer             :: i

        do i = 1, sys_size
            allocate (q_cons_vf(i)%sf(local_start_idx:end_x,local_start_idx:end_y,local_start_idx:end_z))
            allocate (q_prim_vf(i)%sf(local_start_idx:end_x,local_start_idx:end_y,local_start_idx:end_z))
        end do

        if (ib) then
            allocate (ib_markers%sf(local_start_idx:end_x,local_start_idx:end_y,local_start_idx:end_z))
        end if

        if (chemistry) then
            allocate (q_T_sf%sf(local_start_idx:end_x,local_start_idx:end_y,local_start_idx:end_z))
        end if

    end subroutine s_allocate_field_arrays

    !> Read the raw data files present in the corresponding time-step directory and to populate the associated grid and conservative
    !! variables.
    impure subroutine s_read_serial_data_files(t_step)

        integer, intent(in)                                      :: t_step
        character(LEN=len_trim(case_dir) + 2*name_len)           :: t_step_dir
        character(LEN=len_trim(case_dir) + 3*name_len)           :: file_loc
        character(LEN=int(floor(log10(real(sys_size, wp)))) + 1) :: file_num
        logical                                                  :: dir_check
        logical                                                  :: file_check
        integer                                                  :: i

        write (t_step_dir, '(A,I0,A,I0)') '/p_all/p', proc_rank, '/', t_step
        t_step_dir = trim(case_dir) // trim(t_step_dir)

        file_loc = trim(t_step_dir) // '/.'
        call my_inquire(file_loc, dir_check)

        if (dir_check .neqv. .true.) then
            call s_mpi_abort('Time-step folder ' // trim(t_step_dir) // ' is missing. Exiting.')
        end if

        if (bc_io) then
            call s_read_serial_boundary_condition_files(t_step_dir, bc_type)
        else
            call s_assign_default_bc_type(bc_type)
        end if

        ! Pass explicit slices so the dummy `dimension(-1:)` / `dimension(0:)` arguments map to the correct interior indices of the
        ! actual arrays. Without slicing, when offset_x%beg or buff_size > 0 (i.e. format=1 parallel 3D ranks), Fortran's
        ! assumed-shape re-mapping shifts the read by that many slots and leaves the last interior cells uninitialized - corrupting
        ! downstream ghost-cell extrapolation.
        call s_read_grid_data_direction(t_step_dir, 'x', x_cb(-1:m), dx(0:m), x_cc(0:m), m)

        if (n > 0) then
            call s_read_grid_data_direction(t_step_dir, 'y', y_cb(-1:n), dy(0:n), y_cc(0:n), n)

            if (p > 0) then
                call s_read_grid_data_direction(t_step_dir, 'z', z_cb(-1:p), dz(0:p), z_cc(0:p), p)
            end if
        end if

        do i = 1, sys_size
            write (file_num, '(I0)') i
            file_loc = trim(t_step_dir) // '/q_cons_vf' // trim(file_num) // '.dat'
            inquire (FILE=trim(file_loc), EXIST=file_check)

            if (file_check) then
                open (1, FILE=trim(file_loc), form='unformatted', STATUS='old', ACTION='read')
                read (1) q_cons_vf(i)%sf(0:m,0:n,0:p)
                close (1)
            else if (bubbles_lagrange .and. i == beta_idx) then
                ! beta (Lagrangian void fraction) is not written by pre_process for t_step_start; initialize to zero.
                q_cons_vf(i)%sf(0:m,0:n,0:p) = 0._wp
            else
                call s_mpi_abort('File q_cons_vf' // trim(file_num) // '.dat is missing in ' // trim(t_step_dir) // '. Exiting.')
            end if
        end do

        call s_read_ib_data_files(t_step_dir)

    end subroutine s_read_serial_data_files

    !> Parallel-read the raw data files present in the corresponding time-step directory and to populate the associated grid and
    !! conservative variables.
    impure subroutine s_read_parallel_data_files(t_step)

        integer, intent(in) :: t_step

#ifdef MFC_MPI
        real(wp), allocatable, dimension(:)  :: x_cb_glb, y_cb_glb, z_cb_glb
        integer                              :: ifile, ierr, data_size, filetype, stride
        integer, dimension(MPI_STATUS_SIZE)  :: status
        integer(KIND=MPI_OFFSET_KIND)        :: disp
        integer(KIND=MPI_OFFSET_KIND)        :: m_MOK, n_MOK, p_MOK
        integer(KIND=MPI_OFFSET_KIND)        :: WP_MOK, var_MOK, str_MOK
        integer(KIND=MPI_OFFSET_KIND)        :: NVARS_MOK
        integer(KIND=MPI_OFFSET_KIND)        :: MOK
        integer(kind=MPI_OFFSET_KIND)        :: offset
        character(LEN=path_len + 2*name_len) :: file_loc
        logical                              :: file_exist
        character(len=10)                    :: t_step_string
        integer                              :: i

        allocate (x_cb_glb(-1:m_glb))
        allocate (y_cb_glb(-1:n_glb))
        allocate (z_cb_glb(-1:p_glb))

        if (down_sample) then
            stride = 3
        else
            stride = 1
        end if

        file_loc = trim(case_dir) // '/restart_data' // trim(mpiiofs) // 'x_cb.dat'
        inquire (FILE=trim(file_loc), EXIST=file_exist)

        if (file_exist) then
            data_size = m_glb + 2
            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

            call MPI_TYPE_VECTOR(data_size, 1, stride, mpi_p, filetype, ierr)
            call MPI_TYPE_COMMIT(filetype, ierr)

            offset = 0
            call MPI_FILE_SET_VIEW(ifile, offset, mpi_p, filetype, 'native', mpi_info_int, ierr)

            call MPI_FILE_READ(ifile, x_cb_glb, data_size, mpi_p, status, ierr)
            call MPI_FILE_CLOSE(ifile, ierr)
        else
            call s_mpi_abort('File ' // trim(file_loc) // ' is missing. Exiting.')
        end if

        x_cb(-1:m) = x_cb_glb((start_idx(1) - 1):(start_idx(1) + m))
        dx(0:m) = x_cb(0:m) - x_cb(-1:m - 1)
        x_cc(0:m) = x_cb(-1:m - 1) + dx(0:m)/2._wp

        if (n > 0) then
            file_loc = trim(case_dir) // '/restart_data' // trim(mpiiofs) // 'y_cb.dat'
            inquire (FILE=trim(file_loc), EXIST=file_exist)

            if (file_exist) then
                data_size = n_glb + 2
                call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

                call MPI_TYPE_VECTOR(data_size, 1, stride, mpi_p, filetype, ierr)
                call MPI_TYPE_COMMIT(filetype, ierr)

                offset = 0
                call MPI_FILE_SET_VIEW(ifile, offset, mpi_p, filetype, 'native', mpi_info_int, ierr)

                call MPI_FILE_READ(ifile, y_cb_glb, data_size, mpi_p, status, ierr)
                call MPI_FILE_CLOSE(ifile, ierr)
            else
                call s_mpi_abort('File ' // trim(file_loc) // ' is missing. Exiting.')
            end if

            y_cb(-1:n) = y_cb_glb((start_idx(2) - 1):(start_idx(2) + n))
            dy(0:n) = y_cb(0:n) - y_cb(-1:n - 1)
            y_cc(0:n) = y_cb(-1:n - 1) + dy(0:n)/2._wp

            if (p > 0) then
                file_loc = trim(case_dir) // '/restart_data' // trim(mpiiofs) // 'z_cb.dat'
                inquire (FILE=trim(file_loc), EXIST=file_exist)

                if (file_exist) then
                    data_size = p_glb + 2
                    call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

                    call MPI_TYPE_VECTOR(data_size, 1, stride, mpi_p, filetype, ierr)
                    call MPI_TYPE_COMMIT(filetype, ierr)

                    offset = 0
                    call MPI_FILE_SET_VIEW(ifile, offset, mpi_p, filetype, 'native', mpi_info_int, ierr)

                    call MPI_FILE_READ(ifile, z_cb_glb, data_size, mpi_p, status, ierr)
                    call MPI_FILE_CLOSE(ifile, ierr)
                else
                    call s_mpi_abort('File ' // trim(file_loc) // ' is missing. Exiting.')
                end if

                z_cb(-1:p) = z_cb_glb((start_idx(3) - 1):(start_idx(3) + p))
                dz(0:p) = z_cb(0:p) - z_cb(-1:p - 1)
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
    impure subroutine s_read_parallel_conservative_data(t_step, m_MOK, n_MOK, p_MOK, WP_MOK, MOK, str_MOK, NVARS_MOK)

        integer, intent(in)                          :: t_step
        integer(KIND=MPI_OFFSET_KIND), intent(inout) :: m_MOK, n_MOK, p_MOK
        integer(KIND=MPI_OFFSET_KIND), intent(inout) :: WP_MOK, MOK, str_MOK, NVARS_MOK
        integer                                      :: ifile, ierr, data_size
        integer, dimension(MPI_STATUS_SIZE)          :: status
        integer(KIND=MPI_OFFSET_KIND)                :: disp, var_MOK
        character(LEN=path_len + 2*name_len)         :: file_loc
        logical                                      :: file_exist
        character(len=10)                            :: t_step_string
        integer                                      :: i

        if (file_per_process) then
            call s_int_to_str(t_step, t_step_string)
            write (file_loc, '(I0,A1,I7.7,A)') t_step, '_', proc_rank, '.dat'
            file_loc = trim(case_dir) // '/restart_data/lustre_' // trim(t_step_string) // trim(mpiiofs) // trim(file_loc)
            inquire (FILE=trim(file_loc), EXIST=file_exist)

            if (file_exist) then
                call MPI_FILE_OPEN(MPI_COMM_SELF, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

                if (down_sample) then
                    call s_initialize_mpi_data_ds(q_cons_temp)
                else
                    if (ib) then
                        call s_initialize_mpi_data(q_cons_vf, ib_markers)
                    else
                        call s_initialize_mpi_data(q_cons_vf)
                    end if
                end if

                if (down_sample) then
                    data_size = (m + 3)*(n + 3)*(p + 3)
                else
                    data_size = (m + 1)*(n + 1)*(p + 1)
                end if

                m_MOK = int(m_glb + 1, MPI_OFFSET_KIND)
                n_MOK = int(n_glb + 1, MPI_OFFSET_KIND)
                p_MOK = int(p_glb + 1, MPI_OFFSET_KIND)
                WP_MOK = int(storage_size(0._stp)/8, MPI_OFFSET_KIND)
                MOK = int(1._wp, MPI_OFFSET_KIND)
                str_MOK = int(name_len, MPI_OFFSET_KIND)
                NVARS_MOK = int(sys_size, MPI_OFFSET_KIND)

                if (bubbles_euler .or. elasticity .or. mhd) then
                    do i = 1, sys_size
                        var_MOK = int(i, MPI_OFFSET_KIND)
                        call MPI_FILE_READ_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size*mpi_io_type, mpi_io_p, status, ierr)
                    end do
                else
                    do i = 1, sys_size
                        var_MOK = int(i, MPI_OFFSET_KIND)
                        call MPI_FILE_READ_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size*mpi_io_type, mpi_io_p, status, ierr)
                    end do
                end if

                call s_mpi_barrier()
                call MPI_FILE_CLOSE(ifile, ierr)

                if (down_sample) then
                    do i = 1, sys_size
                        q_cons_vf(i)%sf(0:m,0:n,0:p) = q_cons_temp(i)%sf(0:m,0:n,0:p)
                    end do
                end if

                call s_read_ib_data_files(trim(case_dir) // '/restart_data' // trim(mpiiofs), t_step)
            else
                call s_mpi_abort('File ' // trim(file_loc) // ' is missing. Exiting.')
            end if
        else
            write (file_loc, '(I0,A)') t_step, '.dat'
            file_loc = trim(case_dir) // '/restart_data' // trim(mpiiofs) // trim(file_loc)
            inquire (FILE=trim(file_loc), EXIST=file_exist)

            if (file_exist) then
                call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)

                call s_setup_mpi_io_params(data_size, m_MOK, n_MOK, p_MOK, WP_MOK, MOK, str_MOK, NVARS_MOK)

                do i = 1, sys_size
                    var_MOK = int(i, MPI_OFFSET_KIND)

                    disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                    call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, MPI_IO_DATA%view(i), 'native', mpi_info_int, ierr)
                    call MPI_FILE_READ_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size*mpi_io_type, mpi_io_p, status, ierr)
                end do

                call s_mpi_barrier()
                call MPI_FILE_CLOSE(ifile, ierr)

                call s_read_ib_data_files(trim(case_dir) // '/restart_data' // trim(mpiiofs), t_step)
            else
                call s_mpi_abort('File ' // trim(file_loc) // ' is missing. Exiting.')
            end if
        end if

    end subroutine s_read_parallel_conservative_data
#endif

    !> Computation of parameters, allocation procedures, and/or any other tasks needed to properly setup the module
    impure subroutine s_initialize_data_input_module

        integer :: i

        allocate (q_cons_vf(1:sys_size))
        allocate (q_prim_vf(1:sys_size))
        allocate (q_cons_temp(1:sys_size))

        if (n > 0) then
            if (p > 0) then
                call s_allocate_field_arrays(-buff_size, m + buff_size, n + buff_size, p + buff_size)
                if (down_sample) then
                    do i = 1, sys_size
                        allocate (q_cons_temp(i)%sf(-1:m + 1,-1:n + 1,-1:p + 1))
                    end do
                end if
            else
                call s_allocate_field_arrays(-buff_size, m + buff_size, n + buff_size, 0)
            end if
        else
            call s_allocate_field_arrays(-buff_size, m + buff_size, 0, 0)
        end if

        allocate (bc_type(1:num_dims,1:2))

        allocate (bc_type(1, 1)%sf(0:0,0:n,0:p))
        allocate (bc_type(1, 2)%sf(0:0,0:n,0:p))
        if (n > 0) then
            allocate (bc_type(2, 1)%sf(-buff_size:m + buff_size,0:0,0:p))
            allocate (bc_type(2, 2)%sf(-buff_size:m + buff_size,0:0,0:p))
            if (p > 0) then
                allocate (bc_type(3, 1)%sf(-buff_size:m + buff_size,-buff_size:n + buff_size,0:0))
                allocate (bc_type(3, 2)%sf(-buff_size:m + buff_size,-buff_size:n + buff_size,0:0))
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

        integer :: i

        do i = 1, sys_size
            deallocate (q_cons_vf(i)%sf)
            deallocate (q_prim_vf(i)%sf)
            if (down_sample) then
                deallocate (q_cons_temp(i)%sf)
            end if
        end do

        deallocate (q_cons_vf)
        deallocate (q_prim_vf)
        deallocate (q_cons_temp)

        if (ib) then
            deallocate (ib_markers%sf)
        end if

        if (chemistry) then
            deallocate (q_T_sf%sf)
        end if

        deallocate (bc_type(1, 1)%sf, bc_type(1, 2)%sf)
        if (n > 0) then
            deallocate (bc_type(2, 1)%sf, bc_type(2, 2)%sf)
            if (p > 0) then
                deallocate (bc_type(3, 1)%sf, bc_type(3, 2)%sf)
            end if
        end if

        deallocate (bc_type)

        s_read_data_files => null()

        call s_free_amr_data()

    end subroutine s_finalize_data_input_module

    !> Reconstruct level-1 fine cell boundaries fcb(-1:nfine) by bisecting the coarse cells of pcb, starting at this rank's LOCAL
    !! coarse index lo_local. Mirrors s_build_level_coords (m_amr) but keeps only the boundaries needed by the mesh.
    pure subroutine s_amr_reconstruct_fine_cb(pcb, pcb_lb, lo_local, nfine, fcb)

        real(wp), intent(in)                 :: pcb(:)
        integer, intent(in)                  :: pcb_lb, lo_local, nfine
        real(wp), allocatable, intent(inout) :: fcb(:)
        integer                              :: fi, c, off
        real(wp)                             :: xl, xr, xm

        off = 1 - pcb_lb  ! pcb(j) = coarse_cb(j + pcb_lb - 1); coarse_cb(c) = pcb(c + off)
        do fi = 0, nfine
            c = lo_local + fi/2
            xl = pcb(c - 1 + off)
            xr = pcb(c + off)
            xm = 0.5_wp*(xl + xr)
            if (mod(fi, 2) == 0) then
                fcb(fi - 1) = xl
                fcb(fi) = xm
            else
                fcb(fi) = xr
            end if
        end do

    end subroutine s_amr_reconstruct_fine_cb

    !> Populate block slot `k` metadata, allocate its conservative fields, and reconstruct its fine coordinates from the coarse cell
    !! boundaries (x_cb/y_cb/z_cb, already read for this t_step). isect_lo is the block's global coarse origin; sidx is this rank's
    !! global coarse origin (0 for a single-rank/no-MPI run).
    impure subroutine s_setup_amr_block(k, reg, isect_lo, sidx, fm, fn, fp)

        integer, intent(in) :: k, reg(6), isect_lo(3), sidx(3), fm, fn, fp
        integer             :: i

        amr_fine(k)%lo = reg(1:3)
        amr_fine(k)%hi = reg(4:6)
        amr_fine(k)%m = fm; amr_fine(k)%n = fn; amr_fine(k)%p = fp

        allocate (amr_fine(k)%q_cons(1:sys_size))
        do i = 1, sys_size
            allocate (amr_fine(k)%q_cons(i)%sf(0:fm,0:fn,0:fp))
        end do

        ! isect_lo is GLOBAL; sidx is this rank's global origin (0 for a single-rank/no-MPI run), so
        ! isect_lo - sidx is the LOCAL coarse index whose x_cb slice the bisection reads.
        allocate (amr_fine(k)%x_cb(-1:fm))
        call s_amr_reconstruct_fine_cb(x_cb, lbound(x_cb, 1), isect_lo(1) - sidx(1), fm, amr_fine(k)%x_cb)
        if (n > 0) then
            allocate (amr_fine(k)%y_cb(-1:fn))
            call s_amr_reconstruct_fine_cb(y_cb, lbound(y_cb, 1), isect_lo(2) - sidx(2), fn, amr_fine(k)%y_cb)
        end if
        if (p > 0) then
            allocate (amr_fine(k)%z_cb(-1:fp))
            call s_amr_reconstruct_fine_cb(z_cb, lbound(z_cb, 1), isect_lo(3) - sidx(3), fp, amr_fine(k)%z_cb)
        end if

    end subroutine s_setup_amr_block

    !> Read the AMR fine-level restart file for t_step (mirrors s_read_amr_restart in m_amr, serial and parallel branches) and store
    !! this rank's owned block pieces for the post-process overlay. No-op when amr is off or the file is absent.
    impure subroutine s_read_amr_data(t_step)

        integer, intent(in)                  :: t_step
        character(LEN=path_len + 3*name_len) :: file_loc
        logical                              :: file_exist
        integer                              :: k, i, nblk, ghdr(3), reg(6), rm, rn, rp
        integer                              :: sidx(3), ext(3), isect_lo(3), isect_hi(3), fm, fn, fp, d
        integer                              :: have_loc, have_glb
        logical                              :: owns

#ifdef MFC_MPI
        integer                             :: ifile, ierr, cnt, idx, fi, fj, fk, ibytes, sbytes
        integer                             :: myext(3)
        integer, allocatable                :: wext(:), rext(:)
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(kind=MPI_OFFSET_KIND)       :: my_cnt, my_off, tot_cnt, disp0, ddisp
        real(stp), allocatable              :: buf(:)
#endif

        call s_free_amr_data()
        if (.not. amr) return

        ! this rank's subdomain in global coarse indices (mirror s_amr_compute_isect's sidx/ext). start_idx is
        ! allocated only under MPI; for a single-rank/no-MPI run the global origin is 0.
        sidx = 0; ext = 0
        ext(1) = m
        if (n > 0) ext(2) = n
        if (p > 0) ext(3) = p
#ifdef MFC_MPI
        sidx(1) = start_idx(1)
        if (n > 0) sidx(2) = start_idx(2)
        if (p > 0) sidx(3) = start_idx(3)
#endif

        if (.not. parallel_io) then
            write (file_loc, '(A,I0,A,I0,A)') trim(case_dir) // '/p_all/p', proc_rank, '/', t_step, '/amr_fine.dat'
        else
            write (file_loc, '(A,I0,A)') 'amr_', t_step, '.dat'
            file_loc = trim(case_dir) // '/restart_data' // trim(mpiiofs) // trim(file_loc)
        end if
        inquire (FILE=trim(file_loc), EXIST=file_exist)
        ! all ranks must agree: in serial (per-rank-file) mode a partially present p_all tree would
        ! otherwise mix fine-overlay and coarse-only ranks with no message unless rank 0 was the
        ! missing one (mirrors the sim reader's allreduce-min agreement)
        have_loc = merge(1, 0, file_exist)
        call s_mpi_allreduce_integer_min(have_loc, have_glb)
        if (have_glb == 0) then
            if (proc_rank == 0 .and. file_exist) print '(A,I0,A)', ' [amr] post: AMR fine-block file(s) at t_step ', t_step, &
                & ' are missing on some ranks; writing the coarse mesh only'
            if (proc_rank == 0 .and. .not. file_exist) print '(A,I0,A)', ' [amr] post: no AMR fine-block file at t_step ', &
                & t_step, '; writing the coarse mesh only'
            return
        end if

        if (.not. parallel_io) then
            open (2, FILE=trim(file_loc), form='unformatted', ACTION='read', STATUS='old')
            read (2) ghdr
            ! the layout offsets depend on both: a stale file from a different run configuration
            ! would otherwise misalign every record (mirrors the sim reader's header validation)
            if (ghdr(1) /= num_procs) then
                call s_mpi_abort('amr post: the AMR fine-block file was written with a different rank count; ' &
                                 & // 'run post_process with the same number of ranks as the simulation')
            end if
            if (ghdr(3) /= sys_size) then
                call s_mpi_abort('amr post: the AMR fine-block file was written with a different number of ' &
                                 & // 'conserved variables; the physics configuration must match the simulation')
            end if
            nblk = ghdr(2)
            allocate (amr_fine(nblk))
            amr_num_fine = 0
            do k = 1, nblk
                read (2) reg, rm, rn, rp
                do d = 1, 3
                    isect_lo(d) = max(reg(d), sidx(d))
                    isect_hi(d) = min(reg(3 + d), sidx(d) + ext(d))
                end do
                owns = isect_lo(1) <= isect_hi(1)
                if (n > 0) owns = owns .and. isect_lo(2) <= isect_hi(2)
                if (p > 0) owns = owns .and. isect_lo(3) <= isect_hi(3)
                if (.not. owns) cycle  ! writer emitted no data record for a block this rank does not own
                fm = 2*(isect_hi(1) - isect_lo(1) + 1) - 1
                fn = 0; fp = 0
                if (n > 0) fn = 2*(isect_hi(2) - isect_lo(2) + 1) - 1
                if (p > 0) fp = 2*(isect_hi(3) - isect_lo(3) + 1) - 1
                if (fm /= rm .or. fn /= rn .or. fp /= rp) call s_mpi_abort('amr post: fine-extent mismatch vs the ' &
                    & // 'AMR restart file (identical rank count and decomposition required). Exiting.')
                amr_num_fine = amr_num_fine + 1
                call s_setup_amr_block(amr_num_fine, reg, isect_lo, sidx, fm, fn, fp)
                do i = 1, sys_size
                    read (2) amr_fine(amr_num_fine)%q_cons(i)%sf(0:fm,0:fn,0:fp)
                end do
            end do
            close (2)
        else
#ifdef MFC_MPI
            ibytes = storage_size(0)/8; sbytes = storage_size(0._stp)/8
            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
            call MPI_FILE_READ_AT_ALL(ifile, int(0, MPI_OFFSET_KIND), ghdr, 3, MPI_INTEGER, status, ierr)
            ! the parallel layout is per-rank slices concatenated in WRITER rank order: a different
            ! rank count or decomposition would silently misalign every block - fail closed instead
            if (ghdr(1) /= num_procs) then
                call s_mpi_abort('amr post: the AMR fine-block file was written with a different rank count; ' &
                                 & // 'run post_process with the same number of ranks as the simulation')
            end if
            if (ghdr(3) /= sys_size) then
                call s_mpi_abort('amr post: the AMR fine-block file was written with a different number of ' &
                                 & // 'conserved variables; the physics configuration must match the simulation')
            end if
            nblk = ghdr(2)
            allocate (amr_fine(nblk))
            allocate (wext(3*num_procs), rext(3*num_procs))
            amr_num_fine = 0
            disp0 = int(3*ibytes, MPI_OFFSET_KIND)
            do k = 1, nblk
                call MPI_FILE_READ_AT_ALL(ifile, disp0, reg, 6, MPI_INTEGER, status, ierr)
                call MPI_FILE_READ_AT_ALL(ifile, disp0 + int(6*ibytes, MPI_OFFSET_KIND), wext, 3*num_procs, MPI_INTEGER, status, &
                                          & ierr)
                do d = 1, 3
                    isect_lo(d) = max(reg(d), sidx(d))
                    isect_hi(d) = min(reg(3 + d), sidx(d) + ext(d))
                end do
                owns = isect_lo(1) <= isect_hi(1)
                if (n > 0) owns = owns .and. isect_lo(2) <= isect_hi(2)
                if (p > 0) owns = owns .and. isect_lo(3) <= isect_hi(3)
                fm = 2*max(isect_hi(1) - isect_lo(1) + 1, 0) - 1
                fn = 0; fp = 0
                if (n > 0) fn = 2*max(isect_hi(2) - isect_lo(2) + 1, 0) - 1
                if (p > 0) fp = 2*max(isect_hi(3) - isect_lo(3) + 1, 0) - 1
                ! validate the writer's per-rank layout against this run's decomposition (catches a
                ! load_balance simulation, whose weighted splits post never reproduces)
                myext = 0
                if (owns) myext = [fm, fn, fp]
                call MPI_ALLGATHER(myext, 3, MPI_INTEGER, rext, 3, MPI_INTEGER, MPI_COMM_WORLD, ierr)
                if (any(rext /= wext)) then
                    call s_mpi_abort('amr post: the per-rank fine-block layout in the file does not match this ' &
                                     & // 'decomposition (e.g. the simulation used load_balance); the AMR overlay ' &
                                     & // 'requires the writing decomposition')
                end if
                cnt = sys_size*(fm + 1)*(fn + 1)*(fp + 1)
                if (.not. owns) cnt = 0
                my_cnt = int(cnt, MPI_OFFSET_KIND)
                my_off = int(0, MPI_OFFSET_KIND)
                call MPI_EXSCAN(my_cnt, my_off, 1, MPI_OFFSET, MPI_SUM, MPI_COMM_WORLD, ierr)
                if (proc_rank == 0) my_off = int(0, MPI_OFFSET_KIND)
                call MPI_ALLREDUCE(my_cnt, tot_cnt, 1, MPI_OFFSET, MPI_SUM, MPI_COMM_WORLD, ierr)
                ddisp = disp0 + int((6 + 3*num_procs)*ibytes, MPI_OFFSET_KIND)
                allocate (buf(max(cnt, 1)))
                call MPI_FILE_READ_AT_ALL(ifile, ddisp + my_off*int(sbytes, MPI_OFFSET_KIND), buf, cnt*mpi_io_type, mpi_io_p, &
                                          & status, ierr)
                if (owns) then
                    amr_num_fine = amr_num_fine + 1
                    call s_setup_amr_block(amr_num_fine, reg, isect_lo, sidx, fm, fn, fp)
                    idx = 0
                    do i = 1, sys_size
                        do fk = 0, fp
                            do fj = 0, fn
                                do fi = 0, fm
                                    idx = idx + 1
                                    amr_fine(amr_num_fine)%q_cons(i)%sf(fi, fj, fk) = buf(idx)
                                end do
                            end do
                        end do
                    end do
                end if
                deallocate (buf)
                disp0 = ddisp + tot_cnt*int(sbytes, MPI_OFFSET_KIND)
            end do
            call MPI_FILE_CLOSE(ifile, ierr)
#endif
        end if

        if (proc_rank == 0) print '(A,I0,A,I0,A)', ' [amr] post: read ', ghdr(2), ' fine block(s) (', amr_num_fine, &
            & ' owned by rank 0)'

    end subroutine s_read_amr_data

    !> Release the stored AMR fine-block pieces.
    impure subroutine s_free_amr_data()

        integer :: k, i

        if (allocated(amr_fine)) then
            do k = 1, amr_num_fine
                if (allocated(amr_fine(k)%q_cons)) then
                    do i = 1, sys_size
                        if (associated(amr_fine(k)%q_cons(i)%sf)) deallocate (amr_fine(k)%q_cons(i)%sf)
                    end do
                    deallocate (amr_fine(k)%q_cons)
                end if
                if (allocated(amr_fine(k)%x_cb)) deallocate (amr_fine(k)%x_cb)
                if (allocated(amr_fine(k)%y_cb)) deallocate (amr_fine(k)%y_cb)
                if (allocated(amr_fine(k)%z_cb)) deallocate (amr_fine(k)%z_cb)
            end do
            deallocate (amr_fine)
        end if
        amr_num_fine = 0

    end subroutine s_free_amr_data

end module m_data_input
