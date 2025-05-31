!>
!! @file m_boundary_conditions.fpp
!! @brief Contains module m_boundary_conditions

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

    implicit none

    real(wp) :: x_centroid, y_centroid, z_centroid
    real(wp) :: length_x, length_y, length_z
    real(wp) :: radius
    type(bounds_info) :: x_boundary, y_boundary, z_boundary  !<

    private; public :: s_apply_boundary_patches, &
 s_write_serial_boundary_condition_files, &
 s_write_parallel_boundary_condition_files

contains
    impure subroutine s_line_segment_bc(patch_id, q_prim_vf, bc_type)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        type(integer_field), dimension(1:num_dims, -1:1), intent(inout) :: bc_type
        integer, intent(in) :: patch_id

        integer :: j

        ! Patch is a vertical line at x_beg or x_end
        if (patch_bc(patch_id)%dir == 1) then
            y_centroid = patch_bc(patch_id)%centroid(2)
            length_y = patch_bc(patch_id)%length(2)

            y_boundary%beg = y_centroid - 0.5_wp*length_y
            y_boundary%end = y_centroid + 0.5_wp*length_y

            ! Patch is a vertical line at x_beg and x_beg is a domain boundary
            #:for BOUND, X, LOC in [('beg', '-i', -1), ('end', 'm+i', 1)]
                if (patch_bc(patch_id)%loc == ${LOC}$ .and. bc_x%${BOUND}$ < 0) then
                    do j = 0, n
                        if (y_cc(j) > y_boundary%beg .and. y_cc(j) < y_boundary%end) then
                            bc_type(1, ${LOC}$)%sf(0, j, 0) = patch_bc(patch_id)%type
                        end if
                    end do
                end if
            #:endfor
        end if

        ! Patch is a vertical line at y_beg or y_end
        if (patch_bc(patch_id)%dir == 2) then
            x_centroid = patch_bc(patch_id)%centroid(1)
            length_x = patch_bc(patch_id)%length(1)

            x_boundary%beg = x_centroid - 0.5_wp*length_x
            x_boundary%end = x_centroid + 0.5_wp*length_x

            ! Patch is a vertical line at x_beg and x_beg is a domain boundary
            #:for BOUND, Y, LOC in [('beg', '-i', -1), ('end', 'n+i', 1)]
                if (patch_bc(patch_id)%loc == ${LOC}$ .and. bc_y%${BOUND}$ < 0) then
                    do j = 0, m
                        if (x_cc(j) > x_boundary%beg .and. x_cc(j) < x_boundary%end) then
                            bc_type(2, ${LOC}$)%sf(j, 0, 0) = patch_bc(patch_id)%type
                        end if
                    end do
                end if
            #:endfor
        end if

    end subroutine s_line_segment_bc

    impure subroutine s_circle_bc(patch_id, q_prim_vf, bc_type)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        type(integer_field), dimension(1:num_dims, -1:1), intent(inout) :: bc_type

        integer, intent(in) :: patch_id

        integer :: j, k
        if (patch_bc(patch_id)%dir == 1) then
            y_centroid = patch_bc(patch_id)%centroid(2)
            z_centroid = patch_bc(patch_id)%centroid(3)
            radius = patch_bc(patch_id)%radius
            ! Patch is a circle at x_beg and x_beg is a domain boundary
            #:for BOUND, X, LOC in [('beg', '-i', -1), ('end', 'm+i', 1)]
                if (patch_bc(patch_id)%loc == ${LOC}$ .and. bc_x%${BOUND}$ < 0) then
                    do k = 0, p
                        do j = 0, n
                            if ((z_cc(k) - z_centroid)**2._wp + &
                                (y_cc(j) - y_centroid)**2._wp <= radius**2._wp) then
                                bc_type(1, -1)%sf(0, j, k) = patch_bc(patch_id)%type
                            end if
                        end do
                    end do
                end if
            #:endfor
        end if
        if (patch_bc(patch_id)%dir == 2) then
            x_centroid = patch_bc(patch_id)%centroid(1)
            z_centroid = patch_bc(patch_id)%centroid(3)
            radius = patch_bc(patch_id)%radius
            ! Patch is a circle at y_beg and y_beg is a domain boundary
            #:for BOUND, Y, LOC in [('beg', '-i', -1), ('end', 'n+i', 1)]
                if (patch_bc(patch_id)%loc == ${LOC}$ .and. bc_y%${BOUND}$ < 0) then
                    do k = 0, p
                        do j = 0, m
                            if ((z_cc(k) - z_centroid)**2._wp + &
                                (x_cc(j) - x_centroid)**2._wp <= radius**2._wp) then
                                bc_type(2, -1)%sf(j, 0, k) = patch_bc(patch_id)%type
                            end if
                        end do
                    end do
                end if
            #:endfor
        end if
        if (patch_bc(patch_id)%dir == 3) then
            x_centroid = patch_bc(patch_id)%centroid(1)
            y_centroid = patch_bc(patch_id)%centroid(2)
            radius = patch_bc(patch_id)%radius
            #:for BOUND, Z, LOC in [('beg', '-i', -1), ('end', 'p+i', 1)]
                if (patch_bc(patch_id)%loc == ${LOC}$ .and. bc_z%${BOUND}$ < 0) then
                    do k = 0, n
                        do j = 0, m
                            if ((y_cc(k) - y_centroid)**2._wp + &
                                (x_cc(j) - x_centroid)**2._wp <= radius**2._wp) then
                                bc_type(3, -1)%sf(j, k, 0) = patch_bc(patch_id)%type
                            end if
                        end do
                    end do
                end if
            #:endfor
        end if

    end subroutine s_circle_bc

    impure subroutine s_rectangle_bc(patch_id, q_prim_vf, bc_type)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        type(integer_field), dimension(1:num_dims, -1:1), intent(inout) :: bc_type

        integer, intent(in) :: patch_id
        integer :: j, k
        if (patch_bc(patch_id)%dir == 1) then
            y_centroid = patch_bc(patch_id)%centroid(2)
            z_centroid = patch_bc(patch_id)%centroid(3)
            length_y = patch_bc(patch_id)%length(2)
            length_z = patch_bc(patch_id)%length(3)

            y_boundary%beg = y_centroid - 0.5_wp*length_y
            y_boundary%end = y_centroid + 0.5_wp*length_y

            z_boundary%beg = z_centroid - 0.5_wp*length_z
            z_boundary%end = z_centroid + 0.5_wp*length_z
            ! Patch is a circle at x_beg and x_beg is a domain boundary
            #:for BOUND, X, LOC in [('beg', '-i', -1), ('end', 'm+i', 1)]
                if (patch_bc(patch_id)%loc == ${LOC}$ .and. bc_x%${BOUND}$ < 0) then
                    do k = 0, p
                        do j = 0, n
                            if (y_boundary%beg <= y_cc(j) .and. &
                                y_boundary%end >= y_cc(j) .and. &
                                z_boundary%beg <= z_cc(k) .and. &
                                z_boundary%end >= z_cc(k)) then
                                bc_type(1, -1)%sf(0, j, k) = patch_bc(patch_id)%type
                            end if
                        end do
                    end do
                end if
            #:endfor
        end if
        if (patch_bc(patch_id)%dir == 2) then
            x_centroid = patch_bc(patch_id)%centroid(1)
            z_centroid = patch_bc(patch_id)%centroid(3)
            length_x = patch_bc(patch_id)%length(1)
            length_z = patch_bc(patch_id)%length(3)

            x_boundary%beg = x_centroid - 0.5_wp*length_x
            x_boundary%end = x_centroid + 0.5_wp*length_x

            z_boundary%beg = z_centroid - 0.5_wp*length_z
            z_boundary%end = z_centroid + 0.5_wp*length_z
            ! Patch is a circle at y_beg and y_beg is a domain boundary
            #:for BOUND, Y, LOC in [('beg', '-i', -1), ('end', 'n+i', 1)]
                if (patch_bc(patch_id)%loc == ${LOC}$ .and. bc_y%${BOUND}$ < 0) then
                    do k = 0, p
                        do j = 0, m
                            if (x_boundary%beg <= x_cc(j) .and. &
                                x_boundary%end >= x_cc(j) .and. &
                                z_boundary%beg <= z_cc(k) .and. &
                                z_boundary%end >= z_cc(k)) then
                                bc_type(2, -1)%sf(j, 0, k) = patch_bc(patch_id)%type
                            end if
                        end do
                    end do
                end if
            #:endfor
        end if
        if (patch_bc(patch_id)%dir == 3) then
            x_centroid = patch_bc(patch_id)%centroid(1)
            y_centroid = patch_bc(patch_id)%centroid(2)
            length_x = patch_bc(patch_id)%length(1)
            length_y = patch_bc(patch_id)%length(2)

            x_boundary%beg = x_centroid - 0.5_wp*length_x
            x_boundary%end = x_centroid + 0.5_wp*length_x

            y_boundary%beg = y_centroid - 0.5_wp*length_y
            y_boundary%end = y_centroid + 0.5_wp*length_y
            #:for BOUND, Z, LOC in [('beg', '-i', -1), ('end', 'p+i', 1)]
                if (patch_bc(patch_id)%loc == ${LOC}$ .and. bc_z%${BOUND}$ < 0) then
                    do k = 0, n
                        do j = 0, m
                            if (x_boundary%beg <= x_cc(j) .and. &
                                x_boundary%end >= x_cc(j) .and. &
                                y_boundary%beg <= y_cc(k) .and. &
                                y_boundary%end >= y_cc(k)) then
                                bc_type(3, -1)%sf(j, k, 0) = patch_bc(patch_id)%type
                            end if
                        end do
                    end do
                end if
            #:endfor
        end if

    end subroutine s_rectangle_bc

    impure subroutine s_apply_boundary_patches(q_prim_vf, bc_type)

        type(scalar_field), dimension(sys_size) :: q_prim_vf
        type(integer_field), dimension(1:num_dims, -1:1) :: bc_type
        integer :: i

        !< Apply 2D patches to 3D domain
        if (p > 0) then
            do i = 1, num_bc_patches
                if (proc_rank == 0) then
                    print *, 'Processing boundary condition patch', i
                end if

                if (patch_bc(i)%geometry == 2) then
                    call s_circle_bc(i, q_prim_vf, bc_type)
                elseif (patch_bc(i)%geometry == 3) then
                    call s_rectangle_bc(i, q_prim_vf, bc_type)
                end if
            end do
            !< Apply 1D patches to 2D domain
        elseif (n > 0) then
            do i = 1, num_bc_patches
                if (proc_rank == 0) then
                    print *, 'Processing boundary condition patch', i
                end if

                if (patch_bc(i)%geometry == 1) then
                    call s_line_segment_bc(i, q_prim_vf, bc_type)
                end if
            end do
        end if

    end subroutine s_apply_boundary_patches

    impure subroutine s_write_serial_boundary_condition_files(q_prim_vf, bc_type, step_dirpath)

        type(scalar_field), dimension(sys_size) :: q_prim_vf
        type(integer_field), dimension(1:num_dims, -1:1) :: bc_type

        character(LEN=*), intent(in) :: step_dirpath

        integer :: dir, loc
        character(len=path_len) :: file_path

        character(len=10) :: status

        if (old_grid) then
            status = 'old'
        else
            status = 'new'
        end if

        call s_pack_boundary_condition_buffers(q_prim_vf)

        file_path = trim(step_dirpath)//'/bc_type.dat'
        open (1, FILE=trim(file_path), FORM='unformatted', STATUS=status)
        do dir = 1, num_dims
            do loc = -1, 1, 2
                write (1) bc_type(dir, loc)%sf
            end do
        end do
        close (1)

        file_path = trim(step_dirpath)//'/bc_buffers.dat'
        open (1, FILE=trim(file_path), FORM='unformatted', STATUS=status)
        do dir = 1, num_dims
            do loc = -1, 1, 2
                write (1) bc_buffers(dir, loc)%sf
            end do
        end do
        close (1)

    end subroutine s_write_serial_boundary_condition_files

    impure subroutine s_write_parallel_boundary_condition_files(q_prim_vf, bc_type)

        type(scalar_field), dimension(sys_size) :: q_prim_vf
        type(integer_field), dimension(1:num_dims, -1:1) :: bc_type

        integer :: dir, loc
        character(len=path_len) :: file_loc, file_path

#ifdef MFC_MPI
        integer :: ierr
        integer :: file_id
        integer :: offset
        character(len=7) :: proc_rank_str
        logical :: dir_check

        call s_pack_boundary_condition_buffers(q_prim_vf)

        file_loc = trim(case_dir)//'/restart_data/boundary_conditions'
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
        file_path = trim(file_loc)//'/bc_'//trim(proc_rank_str)//'.dat'
        call MPI_File_open(MPI_COMM_SELF, trim(file_path), MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, file_id, ierr)

        offset = 0

        ! Write bc_types
        do dir = 1, num_dims
            do loc = -1, 1, 2
                call MPI_File_set_view(file_id, int(offset, KIND=MPI_ADDRESS_KIND), MPI_INTEGER, MPI_BC_TYPE_TYPE(dir, loc), 'native', MPI_INFO_NULL, ierr)
                call MPI_File_write_all(file_id, bc_type(dir, loc)%sf, 1, MPI_BC_TYPE_TYPE(dir, loc), MPI_STATUS_IGNORE, ierr)
                offset = offset + sizeof(bc_type(dir, loc)%sf)
            end do
        end do

        ! Write bc_buffers
        do dir = 1, num_dims
            do loc = -1, 1, 2
                call MPI_File_set_view(file_id, int(offset, KIND=MPI_ADDRESS_KIND), mpi_p, MPI_BC_BUFFER_TYPE(dir, loc), 'native', MPI_INFO_NULL, ierr)
                call MPI_File_write_all(file_id, bc_buffers(dir, loc)%sf, 1, MPI_BC_BUFFER_TYPE(dir, loc), MPI_STATUS_IGNORE, ierr)
                offset = offset + sizeof(bc_buffers(dir, loc)%sf)
            end do
        end do

        call MPI_File_close(file_id, ierr)
#endif

    end subroutine s_write_parallel_boundary_condition_files

    impure subroutine s_pack_boundary_condition_buffers(q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        integer :: i, j, k

        do k = 0, p
            do j = 0, n
                do i = 1, sys_size
                    bc_buffers(1, -1)%sf(i, j, k) = q_prim_vf(i)%sf(0, j, k)
                    bc_buffers(1, 1)%sf(i, j, k) = q_prim_vf(i)%sf(m, j, k)
                end do
            end do
        end do

        if (n > 0) then
            do k = 0, p
                do j = 1, sys_size
                    do i = 0, m
                        bc_buffers(2, -1)%sf(i, j, k) = q_prim_vf(j)%sf(i, 0, k)
                        bc_buffers(2, 1)%sf(i, j, k) = q_prim_vf(j)%sf(i, n, k)
                    end do
                end do
            end do

            if (p > 0) then
                do k = 1, sys_size
                    do j = 0, n
                        do i = 0, m
                            bc_buffers(3, -1)%sf(i, j, k) = q_prim_vf(k)%sf(i, j, 0)
                            bc_buffers(3, 1)%sf(i, j, k) = q_prim_vf(k)%sf(i, j, p)
                        end do
                    end do
                end do
            end if
        end if

    end subroutine s_pack_boundary_condition_buffers

end module m_boundary_conditions
