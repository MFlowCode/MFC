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
    use m_amr_restart_io
    use m_mpi_proxy
    use m_mpi_common
    use m_compile_specific
    use m_boundary_common
    use m_boundary_io
    use m_helper
    use m_helper_basic, only: f_save_index

    implicit none

    private; public :: s_initialize_data_input_module, s_read_data_files, s_read_serial_data_files, s_read_parallel_data_files, &
        & s_read_amr_data, s_free_amr_data, s_finalize_data_input_module, f_save_exists

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
    type(integer_field), public                              :: ib_markers

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
            call s_initialize_mpi_data(q_cons_vf, ib_markers=ib_markers, ib_mpi_data=MPI_IO_IB_DATA)
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

    !> Does a saved restart exist for this index? Under cfl_dt a step can cross two t_save boundaries, so save indices may be
    !! skipped legitimately and the reader must skip them. Only the shared-file layout is checked: with file_per_process the answer
    !! would not be rank-uniform, so that path stays fail-closed.
    impure function f_save_exists(t_step) result(present_)

        integer, intent(in)                  :: t_step
        logical                              :: present_
        character(LEN=path_len + 2*name_len) :: floc
        character(LEN=name_len)              :: fnum

        present_ = .true.
        if (parallel_io) then
            if (file_per_process) return
            write (fnum, '(I0,A)') t_step, '.dat'
            floc = trim(case_dir) // '/restart_data' // trim(mpiiofs) // trim(fnum)
            inquire (FILE=trim(floc), EXIST=present_)
        else
            ! the serial layout skips the same save indices; every rank saves the same steps, so this stays rank-uniform
            write (floc, '(A,I0,A,I0,A)') trim(case_dir) // '/p_all/p', proc_rank, '/', t_step, '/.'
            call my_inquire(floc, present_)
        end if

    end function f_save_exists

    !> Helper subroutine to read IB data files
    impure subroutine s_read_ib_data_files(file_loc_base, t_step)

        character(len=*), intent(in)         :: file_loc_base
        integer, intent(in), optional        :: t_step
        character(LEN=path_len + 2*name_len) :: file_loc
        logical                              :: file_exist
        integer                              :: ifile, ierr, data_size
        character(len=10)                    :: t_step_string

#ifdef MFC_MPI
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(KIND=MPI_OFFSET_KIND)       :: disp
        integer(KIND=MPI_OFFSET_KIND)       :: m_MOK, n_MOK, p_MOK, MOK, WP_MOK, var_MOK
        integer                             :: save_index
#endif

        if (.not. ib) return

        if (parallel_io .and. file_per_process) then
            call s_int_to_str(t_step, t_step_string)
            write (file_loc, '(A,I0,A,i7.7,A)') 'ib_markers_', t_step, '_', proc_rank, '.dat'
            file_loc = trim(case_dir) // '/restart_data/lustre_' // trim(t_step_string) // '/' // trim(file_loc)
        else if (parallel_io) then
            write (file_loc, '(A)') trim(file_loc_base) // 'ib.dat'
        else
            write (file_loc, '(A)') trim(file_loc_base) // '/ib_data.dat'
        end if
        inquire (FILE=trim(file_loc), EXIST=file_exist)

        if (file_exist) then
            if (parallel_io .and. file_per_process) then
#ifdef MFC_MPI
                call MPI_FILE_OPEN(MPI_COMM_SELF, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
                call s_check_mpi_file_open(ierr, file_loc)

                data_size = (m + 1)*(n + 1)*(p + 1)

                call MPI_FILE_READ(ifile, MPI_IO_IB_DATA%var%sf, data_size, MPI_INTEGER, status, ierr)

                call MPI_FILE_CLOSE(ifile, ierr)
#endif
            else if (parallel_io) then
#ifdef MFC_MPI
                call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
                call s_check_mpi_file_open(ierr, file_loc)

                m_MOK = int(m_glb + 1, MPI_OFFSET_KIND)
                n_MOK = int(n_glb + 1, MPI_OFFSET_KIND)
                p_MOK = int(p_glb + 1, MPI_OFFSET_KIND)
                MOK = int(1._wp, MPI_OFFSET_KIND)
                WP_MOK = int(storage_size(0._stp)/8, MPI_OFFSET_KIND)
                save_index = f_save_index(t_step, cfl_dt, t_step_save)

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

        if (chemistry .or. heat_conduction) then
            allocate (q_T_sf%sf(local_start_idx:end_x,local_start_idx:end_y,local_start_idx:end_z))
            q_T_sf%sf = 0._wp  ! Buffer population reads the interior before anything writes it
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
        integer(kind=8)                      :: file_bytes, bytes_needed
        character(len=10)                    :: case_m_str, file_m_str
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
        inquire (FILE=trim(file_loc), EXIST=file_exist, SIZE=file_bytes)

        ! The grid file holds one cell boundary per value, so its size says which grid wrote the restart. Without
        ! this check a case file whose resolution does not match the run reads past the end of every restart
        ! file and post-processes silently, exiting 0 with NaN-filled output that is indistinguishable from
        ! success until someone plots it. The strided read down_sample performs touches stride*(m_glb + 1) + 1
        ! boundaries of a full-resolution file, so it needs more of the file, not less; only the un-strided read
        ! pins the size exactly, since down-sampling three grids of different size can land on the same m_glb.
        if (file_exist) then
            bytes_needed = (int(stride, 8)*int(m_glb + 1, 8) + 1_8)*int(storage_size(0._wp)/8, 8)
            if (file_bytes < bytes_needed .or. (.not. down_sample .and. file_bytes /= bytes_needed)) then
                call s_int_to_str(m_glb, case_m_str)
                call s_int_to_str(int(file_bytes/int(storage_size(0._wp)/8, 8)) - 2, file_m_str)
                call s_mpi_abort('Restart grid mismatch: this case has m = ' // trim(case_m_str) // ' but ' // trim(file_loc) &
                                 & // ' was written with m = ' // trim(file_m_str) &
                                 & // '. Post-processing must use the same grid as the run that wrote the ' &
                                 & // 'restart files, or it reads past the end of every file and writes NaN.')
            end if
            data_size = m_glb + 2
            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
            call s_check_mpi_file_open(ierr, file_loc)

            call MPI_TYPE_VECTOR(data_size, 1, stride, mpi_p, filetype, ierr)
            call MPI_TYPE_COMMIT(filetype, ierr)

            offset = 0
            call MPI_FILE_SET_VIEW(ifile, offset, mpi_p, filetype, 'native', mpi_info_int, ierr)

            call MPI_FILE_READ(ifile, x_cb_glb, data_size, mpi_p, status, ierr)
            call MPI_FILE_CLOSE(ifile, ierr)
        else
            call s_mpi_abort('File ' // trim(file_loc) // ' is missing. Exiting.')
        end if

        ! Bitwise-consistent grid distribution from the global file
        call s_apply_grid_from_global_dim(x_cb_glb, m_glb, m, start_idx(1), bc_x%beg, bc_x%end, offset_x%beg, offset_x%end, &
                                          & buff_size, buff_size, x_cb, x_cc, dx)

        if (n > 0) then
            file_loc = trim(case_dir) // '/restart_data' // trim(mpiiofs) // 'y_cb.dat'
            inquire (FILE=trim(file_loc), EXIST=file_exist)

            if (file_exist) then
                data_size = n_glb + 2
                call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
                call s_check_mpi_file_open(ierr, file_loc)

                call MPI_TYPE_VECTOR(data_size, 1, stride, mpi_p, filetype, ierr)
                call MPI_TYPE_COMMIT(filetype, ierr)

                offset = 0
                call MPI_FILE_SET_VIEW(ifile, offset, mpi_p, filetype, 'native', mpi_info_int, ierr)

                call MPI_FILE_READ(ifile, y_cb_glb, data_size, mpi_p, status, ierr)
                call MPI_FILE_CLOSE(ifile, ierr)
            else
                call s_mpi_abort('File ' // trim(file_loc) // ' is missing. Exiting.')
            end if

            call s_apply_grid_from_global_dim(y_cb_glb, n_glb, n, start_idx(2), bc_y%beg, bc_y%end, offset_y%beg, offset_y%end, &
                                              & buff_size, buff_size, y_cb, y_cc, dy)

            if (p > 0) then
                file_loc = trim(case_dir) // '/restart_data' // trim(mpiiofs) // 'z_cb.dat'
                inquire (FILE=trim(file_loc), EXIST=file_exist)

                if (file_exist) then
                    data_size = p_glb + 2
                    call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
                    call s_check_mpi_file_open(ierr, file_loc)

                    call MPI_TYPE_VECTOR(data_size, 1, stride, mpi_p, filetype, ierr)
                    call MPI_TYPE_COMMIT(filetype, ierr)

                    offset = 0
                    call MPI_FILE_SET_VIEW(ifile, offset, mpi_p, filetype, 'native', mpi_info_int, ierr)

                    call MPI_FILE_READ(ifile, z_cb_glb, data_size, mpi_p, status, ierr)
                    call MPI_FILE_CLOSE(ifile, ierr)
                else
                    call s_mpi_abort('File ' // trim(file_loc) // ' is missing. Exiting.')
                end if

                call s_apply_grid_from_global_dim(z_cb_glb, p_glb, p, start_idx(3), bc_z%beg, bc_z%end, offset_z%beg, &
                                                  & offset_z%end, buff_size, buff_size, z_cb, z_cc, dz)
            end if
        end if

        call s_read_parallel_conservative_data(t_step, m_MOK, n_MOK, p_MOK, WP_MOK, MOK, str_MOK, NVARS_MOK)

        deallocate (x_cb_glb, y_cb_glb, z_cb_glb)

        if (bc_io) then
            call s_read_parallel_boundary_condition_files(bc_type, strict=.false.)
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
                call s_check_mpi_file_open(ierr, file_loc)

                if (down_sample) then
                    call s_initialize_mpi_data_ds(m, n, p, q_cons_temp)
                else
                    if (ib) then
                        call s_initialize_mpi_data(q_cons_vf, ib_markers=ib_markers, ib_mpi_data=MPI_IO_IB_DATA)
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

                if (bubbles_euler .or. hypoelasticity .or. mhd) then
                    do i = 1, sys_size
                        var_MOK = int(i, MPI_OFFSET_KIND)
                        call MPI_FILE_READ_ALL(ifile, MPI_IO_DATA%var(i)%sf(:,:,:), data_size*mpi_io_type, mpi_io_p, status, ierr)
                    end do
                else
                    do i = 1, sys_size
                        var_MOK = int(i, MPI_OFFSET_KIND)
                        call MPI_FILE_READ_ALL(ifile, MPI_IO_DATA%var(i)%sf(:,:,:), data_size*mpi_io_type, mpi_io_p, status, ierr)
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
                call s_check_mpi_file_open(ierr, file_loc)

                call s_setup_mpi_io_params(data_size, m_MOK, n_MOK, p_MOK, WP_MOK, MOK, str_MOK, NVARS_MOK)

                do i = 1, sys_size
                    var_MOK = int(i, MPI_OFFSET_KIND)

                    disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                    call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, MPI_IO_DATA%view(i), 'native', mpi_info_int, ierr)
                    call MPI_FILE_READ_ALL(ifile, MPI_IO_DATA%var(i)%sf(:,:,:), data_size*mpi_io_type, mpi_io_p, status, ierr)
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

        if (chemistry .or. heat_conduction) then
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

    !> Populate block slot `k` metadata, allocate its conservative fields, and reconstruct its fine coordinates from the coarse cell
    !! boundaries (x_cb/y_cb/z_cb, already read for this t_step). isect_lo is the block's global coarse origin; sidx is this rank's
    !! global coarse origin (0 for a single-rank/no-MPI run). rr is the refinement factor for this block (amr_ref_ratio**level).
    impure subroutine s_setup_amr_block(k, reg, isect_lo, sidx, fm, fn, fp, rr)

        integer, intent(in) :: k, reg(6), isect_lo(3), sidx(3), fm, fn, fp, rr
        integer             :: i

        amr_fine(k)%lo = reg(1:3)
        amr_fine(k)%hi = reg(4:6)
        amr_fine(k)%m = fm; amr_fine(k)%n = fn; amr_fine(k)%p = fp

        allocate (amr_fine(k)%q_cons(1:sys_size))
        do i = 1, sys_size
            allocate (amr_fine(k)%q_cons(i)%sf(0:fm,0:fn,0:fp))
            ! zeroed: the reader fills 1..sys_size_io but the conversion reads the whole band (post-only beta_idx slot, void
            ! fraction 0 in fine blocks)
            amr_fine(k)%q_cons(i)%sf = 0._wp
        end do

        ! isect_lo - sidx is the local coarse index whose x_cb slice the rr-way subdivision reads
        allocate (amr_fine(k)%x_cb(-1:fm))
        call s_amr_subdivide_cb(x_cb, lbound(x_cb, 1), isect_lo(1) - sidx(1), fm, rr, amr_fine(k)%x_cb)
        if (n > 0) then
            allocate (amr_fine(k)%y_cb(-1:fn))
            call s_amr_subdivide_cb(y_cb, lbound(y_cb, 1), isect_lo(2) - sidx(2), fn, rr, amr_fine(k)%y_cb)
        end if
        if (p > 0) then
            allocate (amr_fine(k)%z_cb(-1:fp))
            call s_amr_subdivide_cb(z_cb, lbound(z_cb, 1), isect_lo(3) - sidx(3), fp, rr, amr_fine(k)%z_cb)
        end if

    end subroutine s_setup_amr_block

    !> This rank's subdomain in global coarse cells: origin sidx (0 for a single rank or no MPI) and extent ext.
    impure subroutine s_amr_local_domain(sidx, ext)

        integer, intent(out) :: sidx(3), ext(3)

        sidx = 0; ext = 0
        ext(1) = m
        if (n > 0) ext(2) = n
        if (p > 0) ext(3) = p
#ifdef MFC_MPI
        sidx(1) = start_idx(1)
        if (n > 0) sidx(2) = start_idx(2)
        if (p > 0) sidx(3) = start_idx(3)
#endif

    end subroutine s_amr_local_domain

    !> Intersection of the block region reg with this rank's subdomain (sidx, ext); owns = non-empty in every active dimension.
    pure subroutine s_amr_block_isect(reg, sidx, ext, isect_lo, isect_hi, owns)

        integer, intent(in)  :: reg(6), sidx(3), ext(3)
        integer, intent(out) :: isect_lo(3), isect_hi(3)
        logical, intent(out) :: owns
        integer              :: d

        do d = 1, 3
            isect_lo(d) = max(reg(d), sidx(d))
            isect_hi(d) = min(reg(3 + d), sidx(d) + ext(d))
        end do
        owns = isect_lo(1) <= isect_hi(1)
        if (n > 0) owns = owns .and. isect_lo(2) <= isect_hi(2)
        if (p > 0) owns = owns .and. isect_lo(3) <= isect_hi(3)

    end subroutine s_amr_block_isect

    !> Refinement factor of a block from its level and fine x-extent fmf (cells - 1); aborts on an inconsistent header.
    impure integer function f_amr_block_ratio(reg, lvl, fmf) result(rr)

        integer, intent(in) :: reg(6), lvl, fmf
        integer             :: cw

        cw = max(reg(4) - reg(1) + 1, 1)
        rr = (fmf + 1)/cw
        if (lvl < 1 .or. rr < 2 .or. mod(fmf + 1, &
            & cw) /= 0) call s_mpi_abort('amr post: malformed fine-block header (level/' &
            & // 'extent inconsistent); the AMR restart writer and reader header layouts differ')

    end function f_amr_block_ratio

    !> The kept sub-box of a block: fine-cell offset foff from the block origin and count fcnt of the intersection.
    pure subroutine s_amr_block_piece(reg, isect_lo, isect_hi, rr, foff, fcnt)

        integer, intent(in)  :: reg(6), isect_lo(3), isect_hi(3), rr
        integer, intent(out) :: foff(3), fcnt(3)

        foff = (isect_lo - reg(1:3))*rr
        fcnt = (isect_hi - isect_lo + 1)*rr
        if (n == 0) then; foff(2) = 0; fcnt(2) = 1; end if
        if (p == 0) then; foff(3) = 0; fcnt(3) = 1; end if

    end subroutine s_amr_block_piece

    !> Read the AMR fine-level restart file for t_step (layout: m_amr_restart_io) and keep this rank's pieces of each block for the
    !! post-process overlay. No-op when amr is off or the file is absent.
    impure subroutine s_read_amr_data(t_step)

        integer, intent(in)                  :: t_step
        character(LEN=path_len + 3*name_len) :: file_loc
        logical                              :: file_exist, owns
        integer                              :: k, i, nblk, ghdr(amr_restart_ghdr_ints), reg(6), lvl, rm, rn, rp
        integer                              :: sidx(3), ext(3), isect_lo(3), isect_hi(3), foff(3), fcnt(3), rr
        integer                              :: have_loc, have_glb, ncoarse, ncoarse_glb
        real(stp), allocatable               :: sbuf(:,:,:,:)

#ifdef MFC_MPI
        integer                             :: ifile, ierr, cnt, fi, fj, fk
        integer(8)                          :: idx, cnt8, plane
        type(t_amr_restart_catalog)         :: cat
        integer, dimension(MPI_STATUS_SIZE) :: status
        real(stp), allocatable              :: buf(:)
#endif

        call s_free_amr_data()
        if (.not. amr) return

        call s_amr_local_domain(sidx, ext)

        file_loc = f_amr_restart_path(t_step)
        inquire (FILE=trim(file_loc), EXIST=file_exist)
        ! all ranks must agree: with per-rank files a partially present tree would mix fine-overlay and coarse-only ranks
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
            call s_amr_restart_check_header(ghdr, sys_size_io, 'amr post', nblk)
            allocate (amr_fine(nblk))
            amr_num_fine = 0
            ncoarse = 0
            do k = 1, nblk
                read (2) reg, lvl, rm, rn, rp  ! region(6), level, fine extents (mirrors s_write_amr_restart)
                call s_amr_block_isect(reg, sidx, ext, isect_lo, isect_hi, owns)
                ! data presence comes from the file: the writer emits a record only for this rank's own blocks (rm < 0
                ! otherwise), and an owner need not be the rank whose subdomain holds the block
                if (rm < 0) then
                    if (owns .and. lvl >= 1) ncoarse = ncoarse + 1  ! refined here, but its data lives in another rank's file
                    cycle
                end if
                rr = f_amr_block_ratio(reg, lvl, rm)
                allocate (sbuf(0:rm,0:rn,0:rp,1:sys_size_io))
                do i = 1, sys_size_io
                    read (2) sbuf(0:rm,0:rn,0:rp,i)
                end do
                if (owns) then
                    call s_amr_block_piece(reg, isect_lo, isect_hi, rr, foff, fcnt)
                    amr_num_fine = amr_num_fine + 1
                    call s_setup_amr_block(amr_num_fine, reg, isect_lo, sidx, fcnt(1) - 1, fcnt(2) - 1, fcnt(3) - 1, rr)
                    do i = 1, sys_size_io
                        amr_fine(amr_num_fine)%q_cons(i)%sf(:,:,:) = sbuf(foff(1):foff(1) + fcnt(1) - 1, &
                                 & foff(2):foff(2) + fcnt(2) - 1,foff(3):foff(3) + fcnt(3) - 1,i)
                    end do
                end if
                deallocate (sbuf)
            end do
            close (2)
            ! with per-rank files a block owned by another rank cannot be rendered here
            call s_mpi_allreduce_integer_max(ncoarse, ncoarse_glb)
            if (proc_rank == 0 .and. ncoarse_glb > 0) print '(A)', &
                & ' [amr] post: WARNING: with parallel_io = F each rank ' &
                & // 'reads only its own restart file, and some refined regions belong to a block another rank owns; those ' &
                & // 'regions are written at coarse resolution. Use parallel_io = T for a complete AMR overlay at np > 1.'
        else
#ifdef MFC_MPI
            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
            call MPI_FILE_READ_AT_ALL(ifile, int(0, MPI_OFFSET_KIND), ghdr, amr_restart_ghdr_ints, MPI_INTEGER, status, ierr)
            call s_amr_restart_check_header(ghdr, sys_size_io, 'amr post', nblk)
            allocate (amr_fine(nblk))
            amr_num_fine = 0
            call s_amr_restart_read_catalog(ifile, nblk, sys_size_io, cat)
            do k = 1, nblk
                reg = cat%reg(:,k); lvl = cat%lvl(k)
                call s_amr_block_isect(reg, sidx, ext, isect_lo, isect_hi, owns)
                owns = owns .and. cat%own(1, k) > 0  ! a block with no owner has no chunk
                rr = 1
                if (owns) rr = f_amr_block_ratio(reg, lvl, cat%own(2, k))
                cnt8 = 0
                if (owns) cnt8 = int(sys_size_io, 8)*(cat%own(2, k) + 1)*(cat%own(3, k) + 1)*(cat%own(4, k) + 1)
                if (cnt8 > huge(cnt)) call s_mpi_abort('amr post: a fine block exceeds the MPI read count limit')
                cnt = int(cnt8)
                allocate (buf(max(cnt, 1)))
                ! collective: every rank calls it, non-participants with count 0
                call MPI_FILE_READ_AT_ALL(ifile, cat%data_disp(k), buf, cnt*mpi_io_type, mpi_io_p, status, ierr)
                if (owns) then
                    call s_amr_block_piece(reg, isect_lo, isect_hi, rr, foff, fcnt)
                    amr_num_fine = amr_num_fine + 1
                    call s_setup_amr_block(amr_num_fine, reg, isect_lo, sidx, fcnt(1) - 1, fcnt(2) - 1, fcnt(3) - 1, rr)
                    ! file order: i, fk, fj, fi (fastest) over the full block extent
                    plane = int(cat%own(3, k) + 1, 8)*(cat%own(2, k) + 1)
                    do i = 1, sys_size_io
                        do fk = 0, fcnt(3) - 1
                            do fj = 0, fcnt(2) - 1
                                do fi = 0, fcnt(1) - 1
                                    idx = (i - 1)*plane*(cat%own(4, &
                                           & k) + 1) + (fk + foff(3))*plane + (fj + foff(2))*int(cat%own(2, k) + 1, &
                                           & 8) + (fi + foff(1)) + 1
                                    amr_fine(amr_num_fine)%q_cons(i)%sf(fi, fj, fk) = buf(idx)
                                end do
                            end do
                        end do
                    end do
                end if
                deallocate (buf)
            end do
            call MPI_FILE_CLOSE(ifile, ierr)
#endif
        end if

        if (proc_rank == 0) print '(A,I0,A,I0,A)', ' [amr] post: read ', nblk, ' fine block(s); rank 0 keeps ', amr_num_fine, &
            & ' piece(s)'

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
