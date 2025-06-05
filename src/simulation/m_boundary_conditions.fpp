! @file m_boundary_conditions.fpp
! @brief Contains module m_boundary_conditions

!> @brief This module contains
module m_boundary_conditions

    use m_derived_types

    use m_global_parameters
#ifdef MFC_MPI
    use mpi
#endif
    use m_delay_file_access

    use m_compile_specific

    use m_boundary_common

contains

    impure subroutine s_read_serial_boundary_condition_files(step_dirpath, bc_type)

        character(LEN=*), intent(in) :: step_dirpath

        type(integer_field), dimension(1:num_dims, -1:1), intent(inout) :: bc_type

        integer :: dir, loc
        logical :: file_exist
        character(len=path_len) :: file_path

        ! Read bc_types
        file_path = trim(step_dirpath)//'/bc_type.dat'
        inquire (FILE=trim(file_path), EXIST=file_exist)
        if (.not. file_exist) then
            call s_mpi_abort(trim(file_path)//' is missing. Exiting.')
        end if

        open (1, FILE=trim(file_path), FORM='unformatted', STATUS='unknown')
        do dir = 1, num_dims
            do loc = -1, 1, 2
                read (1) bc_type(dir, loc)%sf
                !$acc update device(bc_type(dir, loc)%sf)
            end do
        end do
        close (1)

        ! Read bc_buffers
        file_path = trim(step_dirpath)//'/bc_buffers.dat'
        inquire (FILE=trim(file_path), EXIST=file_exist)
        if (.not. file_exist) then
            call s_mpi_abort(trim(file_path)//' is missing. Exiting.')
        end if

        open (1, FILE=trim(file_path), FORM='unformatted', STATUS='unknown')
        do dir = 1, num_dims
            do loc = -1, 1, 2
                read (1) bc_buffers(dir, loc)%sf
                !$acc update device(bc_buffers(dir, loc)%sf)
            end do
        end do
        close (1)

    end subroutine s_read_serial_boundary_condition_files

    impure subroutine s_read_parallel_boundary_condition_files(bc_type)

        type(integer_field), dimension(1:num_dims, -1:1), intent(inout) :: bc_type

        integer :: dir, loc
        character(len=path_len) :: file_loc, file_path

#ifdef MFC_MPI
        integer :: ierr
        integer :: file_id
        integer :: offset
        character(len=7) :: proc_rank_str
        logical :: dir_check

        file_loc = trim(case_dir)//'/restart_data/boundary_conditions'

        if (proc_rank == 0) then
            call my_inquire(file_loc, dir_check)
            if (dir_check .neqv. .true.) then
                call s_mpi_abort(trim(file_loc)//' is missing. Exiting.')
            end if
        end if

        call s_create_mpi_types(bc_type)

        call s_mpi_barrier()

        call DelayFileAccess(proc_rank)

        write (proc_rank_str, '(I7.7)') proc_rank
        file_path = trim(file_loc)//'/bc_'//trim(proc_rank_str)//'.dat'
        call MPI_File_open(MPI_COMM_SELF, trim(file_path), MPI_MODE_RDONLY, MPI_INFO_NULL, file_id, ierr)

        offset = 0

        ! Read bc_types
        do dir = 1, num_dims
            do loc = -1, 1, 2
                call MPI_File_set_view(file_id, int(offset, KIND=MPI_ADDRESS_KIND), MPI_INTEGER, MPI_BC_TYPE_TYPE(dir, loc), 'native', MPI_INFO_NULL, ierr)
                call MPI_File_read_all(file_id, bc_type(dir, loc)%sf, 1, MPI_BC_TYPE_TYPE(dir, loc), MPI_STATUS_IGNORE, ierr)
                offset = offset + sizeof(bc_type(dir, loc)%sf)
                !$acc update device(bc_type(dir, loc)%sf)
            end do
        end do

        ! Read bc_buffers
        do dir = 1, num_dims
            do loc = -1, 1, 2
                call MPI_File_set_view(file_id, int(offset, KIND=MPI_ADDRESS_KIND), mpi_p, MPI_BC_BUFFER_TYPE(dir, loc), 'native', MPI_INFO_NULL, ierr)
                call MPI_File_read_all(file_id, bc_buffers(dir, loc)%sf, 1, MPI_BC_BUFFER_TYPE(dir, loc), MPI_STATUS_IGNORE, ierr)
                offset = offset + sizeof(bc_buffers(dir, loc)%sf)
                !$acc update device(bc_buffers(dir, loc)%sf)
            end do
        end do

        call MPI_File_close(file_id, ierr)
#endif

    end subroutine s_read_parallel_boundary_condition_files

end module m_boundary_conditions
