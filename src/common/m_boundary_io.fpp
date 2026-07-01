!>
!! @file
!! @brief Contains module m_boundary_io

!> @brief Boundary condition restart I/O, capillary/IGR buffer population, and grid-variable buffers
#:include 'case.fpp'
#:include 'macros.fpp'

module m_boundary_io

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_constants
    use m_delay_file_access
    use m_compile_specific
    use m_boundary_primitives

    implicit none

#ifdef MFC_MPI
    integer, dimension(1:3,1:2) :: MPI_BC_TYPE_TYPE
    integer, dimension(1:3,1:2) :: MPI_BC_BUFFER_TYPE
#endif

contains

    !> Create MPI derived datatypes for boundary condition type arrays and buffer arrays used in parallel I/O.
    impure subroutine s_create_mpi_types(bc_type)

        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type

#ifdef MFC_MPI
        integer               :: dir, loc
        integer, dimension(3) :: sf_start_idx, sf_extents_loc
        integer               :: ierr

        do dir = 1, num_dims
            do loc = 1, 2
                sf_start_idx = (/0, 0, 0/)
                sf_extents_loc = shape(bc_type(dir, loc)%sf)

                call MPI_TYPE_CREATE_SUBARRAY(num_dims, sf_extents_loc, sf_extents_loc, sf_start_idx, MPI_ORDER_FORTRAN, &
                                              & MPI_INTEGER, MPI_BC_TYPE_TYPE(dir, loc), ierr)
                call MPI_TYPE_COMMIT(MPI_BC_TYPE_TYPE(dir, loc), ierr)
            end do
        end do

        do dir = 1, num_dims
            do loc = 1, 2
                sf_start_idx = (/0, 0, 0/)
                sf_extents_loc = shape(bc_buffers(dir, loc)%sf)

                call MPI_TYPE_CREATE_SUBARRAY(num_dims, sf_extents_loc*mpi_io_type, sf_extents_loc*mpi_io_type, sf_start_idx, &
                                              & MPI_ORDER_FORTRAN, mpi_io_p, MPI_BC_BUFFER_TYPE(dir, loc), ierr)
                call MPI_TYPE_COMMIT(MPI_BC_BUFFER_TYPE(dir, loc), ierr)
            end do
        end do
#endif

    end subroutine s_create_mpi_types

    !> Write boundary condition type and buffer data to serial (unformatted) restart files.
    subroutine s_write_serial_boundary_condition_files(q_prim_vf, bc_type, step_dirpath, old_grid_in, q_T_sf)

        type(scalar_field), dimension(sys_size), intent(in)        :: q_prim_vf
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        logical, intent(in)                                        :: old_grid_in
        character(LEN=*), intent(in)                               :: step_dirpath
        integer                                                    :: dir, loc
        character(len=path_len)                                    :: file_path
        character(len=10)                                          :: status
        type(scalar_field), optional, intent(in)                   :: q_T_sf

        if (old_grid_in) then
            status = 'old'
        else
            status = 'new'
        end if

        call s_pack_boundary_condition_buffers(q_prim_vf, q_T_sf)

        file_path = trim(step_dirpath) // '/bc_type.dat'
        open (1, FILE=trim(file_path), form='unformatted', STATUS=status)
        do dir = 1, num_dims
            do loc = 1, 2
                write (1) bc_type(dir, loc)%sf
            end do
        end do
        close (1)

        file_path = trim(step_dirpath) // '/bc_buffers.dat'
        open (1, FILE=trim(file_path), form='unformatted', STATUS=status)
        do dir = 1, num_dims
            do loc = 1, 2
                write (1) bc_buffers(dir, loc)%sf
            end do
        end do
        close (1)

    end subroutine s_write_serial_boundary_condition_files

    !> Write boundary condition type and buffer data to per-rank parallel files using MPI I/O.
    subroutine s_write_parallel_boundary_condition_files(q_prim_vf, bc_type, q_T_sf)

        type(scalar_field), dimension(sys_size), intent(in)        :: q_prim_vf
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        integer                                                    :: dir, loc
        character(len=path_len)                                    :: file_loc, file_path
        type(scalar_field), intent(in), optional                   :: q_T_sf

#ifdef MFC_MPI
        integer          :: ierr
        integer          :: file_id
        character(len=7) :: proc_rank_str
        logical          :: dir_check
        integer          :: nelements

        call s_pack_boundary_condition_buffers(q_prim_vf, q_T_sf)

        file_loc = trim(case_dir) // '/restart_data/boundary_conditions'
        if (proc_rank == 0) then
            call my_inquire(file_loc, dir_check)
            if (dir_check .neqv. .true.) then
                call s_create_directory(trim(file_loc))
            end if
        end if

        call s_create_mpi_types(bc_type)

        call s_mpi_barrier()

        call DelayFileAccess(proc_rank)

        write (proc_rank_str, '(I7.7)') proc_rank
        file_path = trim(file_loc) // '/bc_' // trim(proc_rank_str) // '.dat'
        call MPI_File_open(MPI_COMM_SELF, trim(file_path), MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, file_id, ierr)

        ! Write bc_types
        do dir = 1, num_dims
            do loc = 1, 2
#ifdef MFC_MIXED_PRECISION
                nelements = sizeof(bc_type(dir, loc)%sf)
                call MPI_File_write_all(file_id, bc_type(dir, loc)%sf, nelements, MPI_BYTE, MPI_STATUS_IGNORE, ierr)
#else
                nelements = sizeof(bc_type(dir, loc)%sf)/4
                call MPI_File_write_all(file_id, bc_type(dir, loc)%sf, nelements, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
#endif
            end do
        end do

        ! Write bc_buffers
        do dir = 1, num_dims
            do loc = 1, 2
                nelements = sizeof(bc_buffers(dir, loc)%sf)*mpi_io_type/stp
                call MPI_File_write_all(file_id, bc_buffers(dir, loc)%sf, nelements, mpi_io_p, MPI_STATUS_IGNORE, ierr)
            end do
        end do

        call MPI_File_close(file_id, ierr)
#endif

    end subroutine s_write_parallel_boundary_condition_files

    !> Read boundary condition type and buffer data from serial (unformatted) restart files.
    subroutine s_read_serial_boundary_condition_files(step_dirpath, bc_type)

        character(LEN=*), intent(in)                                  :: step_dirpath
        type(integer_field), dimension(1:num_dims,1:2), intent(inout) :: bc_type
        integer                                                       :: dir, loc
        logical                                                       :: file_exist
        character(len=path_len)                                       :: file_path

        ! Read bc_types

        file_path = trim(step_dirpath) // '/bc_type.dat'
        inquire (FILE=trim(file_path), EXIST=file_exist)
        if (.not. file_exist) then
            call s_mpi_abort(trim(file_path) // ' is missing. Exiting.')
        end if

        open (1, FILE=trim(file_path), form='unformatted', STATUS='unknown')
        do dir = 1, num_dims
            do loc = 1, 2
                read (1) bc_type(dir, loc)%sf
                $:GPU_UPDATE(device='[bc_type(dir, loc)%sf]')
            end do
        end do
        close (1)

        ! Read bc_buffers
        file_path = trim(step_dirpath) // '/bc_buffers.dat'
        inquire (FILE=trim(file_path), EXIST=file_exist)
        if (.not. file_exist) then
            call s_mpi_abort(trim(file_path) // ' is missing. Exiting.')
        end if

        open (1, FILE=trim(file_path), form='unformatted', STATUS='unknown')
        do dir = 1, num_dims
            do loc = 1, 2
                read (1) bc_buffers(dir, loc)%sf
                $:GPU_UPDATE(device='[bc_buffers(dir, loc)%sf]')
            end do
        end do
        close (1)

    end subroutine s_read_serial_boundary_condition_files

    !> Read boundary condition type and buffer data from per-rank parallel files using MPI I/O.
    subroutine s_read_parallel_boundary_condition_files(bc_type)

        type(integer_field), dimension(1:num_dims,1:2), intent(inout) :: bc_type
        integer                                                       :: dir, loc
        character(len=path_len)                                       :: file_loc, file_path

#ifdef MFC_MPI
        integer          :: ierr
        integer          :: file_id
        character(len=7) :: proc_rank_str
        logical          :: dir_check
        integer          :: nelements

        file_loc = trim(case_dir) // '/restart_data/boundary_conditions'

        if (proc_rank == 0) then
            call my_inquire(file_loc, dir_check)
            if (dir_check .neqv. .true.) then
                call s_mpi_abort(trim(file_loc) // ' is missing. Exiting.')
            end if
        end if

        call s_create_mpi_types(bc_type)

        call s_mpi_barrier()

        call DelayFileAccess(proc_rank)

        write (proc_rank_str, '(I7.7)') proc_rank
        file_path = trim(file_loc) // '/bc_' // trim(proc_rank_str) // '.dat'
        call MPI_File_open(MPI_COMM_SELF, trim(file_path), MPI_MODE_RDONLY, MPI_INFO_NULL, file_id, ierr)

        ! Read bc_types
        do dir = 1, num_dims
            do loc = 1, 2
#ifdef MFC_MIXED_PRECISION
                nelements = sizeof(bc_type(dir, loc)%sf)
                call MPI_File_read_all(file_id, bc_type(dir, loc)%sf, nelements, MPI_BYTE, MPI_STATUS_IGNORE, ierr)
#else
                nelements = sizeof(bc_type(dir, loc)%sf)/4
                call MPI_File_read_all(file_id, bc_type(dir, loc)%sf, nelements, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
#endif
                $:GPU_UPDATE(device='[bc_type(dir, loc)%sf]')
            end do
        end do

        ! Read bc_buffers
        do dir = 1, num_dims
            do loc = 1, 2
                nelements = sizeof(bc_buffers(dir, loc)%sf)*mpi_io_type/stp
                call MPI_File_read_all(file_id, bc_buffers(dir, loc)%sf, nelements, mpi_io_p, MPI_STATUS_IGNORE, ierr)
                $:GPU_UPDATE(device='[bc_buffers(dir, loc)%sf]')
            end do
        end do

        call MPI_File_close(file_id, ierr)
#endif

    end subroutine s_read_parallel_boundary_condition_files

    !> Pack primitive variable boundary slices into bc_buffers arrays for serialization.
    subroutine s_pack_boundary_condition_buffers(q_prim_vf, q_T_sf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        integer                                             :: i, j, k
        type(scalar_field), intent(in), optional            :: q_T_sf

        do k = 0, p
            do j = 0, n
                do i = 1, sys_size
                    bc_buffers(1, 1)%sf(i, j, k) = q_prim_vf(i)%sf(0, j, k)
                    bc_buffers(1, 2)%sf(i, j, k) = q_prim_vf(i)%sf(m, j, k)
                end do
                if (chemistry .and. present(q_T_sf)) then
                    bc_buffers(1, 1)%sf(sys_size + 1, j, k) = q_T_sf%sf(0, j, k)
                    bc_buffers(1, 2)%sf(sys_size + 1, j, k) = q_T_sf%sf(m, j, k)
                end if
            end do
        end do

        #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
            if (n > 0) then
                do k = 0, p
                    do j = 1, sys_size
                        do i = 0, m
                            bc_buffers(2, 1)%sf(i, j, k) = q_prim_vf(j)%sf(i, 0, k)
                            bc_buffers(2, 2)%sf(i, j, k) = q_prim_vf(j)%sf(i, n, k)
                        end do
                    end do
                    if (chemistry .and. present(q_T_sf)) then
                        do i = 0, m
                            bc_buffers(2, 1)%sf(i, sys_size + 1, k) = q_T_sf%sf(i, 0, k)
                            bc_buffers(2, 2)%sf(i, sys_size + 1, k) = q_T_sf%sf(i, n, k)
                        end do
                    end if
                end do

                #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                    if (p > 0) then
                        do k = 1, sys_size
                            do j = 0, n
                                do i = 0, m
                                    bc_buffers(3, 1)%sf(i, j, k) = q_prim_vf(k)%sf(i, j, 0)
                                    bc_buffers(3, 2)%sf(i, j, k) = q_prim_vf(k)%sf(i, j, p)
                                end do
                            end do
                        end do
                        if (chemistry .and. present(q_T_sf)) then
                            do j = 0, n
                                do i = 0, m
                                    bc_buffers(3, 1)%sf(i, j, sys_size + 1) = q_T_sf%sf(i, j, 0)
                                    bc_buffers(3, 2)%sf(i, j, sys_size + 1) = q_T_sf%sf(i, j, p)
                                end do
                            end do
                        end if
                    end if
                #:endif
            end if
        #:endif

    end subroutine s_pack_boundary_condition_buffers

    !> Initialize the per-cell boundary condition type arrays with the global default BC values.
    subroutine s_assign_default_bc_type(bc_type)

        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type

        bc_type(1, 1)%sf(:,:,:) = int(min(bc_x%beg, 0), kind=1)
        bc_type(1, 2)%sf(:,:,:) = int(min(bc_x%end, 0), kind=1)
        $:GPU_UPDATE(device='[bc_type(1, 1)%sf, bc_type(1, 2)%sf]')

        #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
            if (n > 0) then
                bc_type(2, 1)%sf(:,:,:) = int(min(bc_y%beg, 0), kind=1)
                bc_type(2, 2)%sf(:,:,:) = int(min(bc_y%end, 0), kind=1)
                $:GPU_UPDATE(device='[bc_type(2, 1)%sf, bc_type(2, 2)%sf]')
                #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                    if (p > 0) then
                        bc_type(3, 1)%sf(:,:,:) = int(min(bc_z%beg, 0), kind=1)
                        bc_type(3, 2)%sf(:,:,:) = int(min(bc_z%end, 0), kind=1)
                        $:GPU_UPDATE(device='[bc_type(3, 1)%sf, bc_type(3, 2)%sf]')
                    end if
                #:endif
            end if
        #:endif

    end subroutine s_assign_default_bc_type

end module m_boundary_io
