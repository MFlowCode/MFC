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
    use m_mpi_common, only: s_mpi_allreduce_min, s_mpi_allreduce_max
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

    !> Populate ghost cell buffers for the color function and its divergence used in capillary surface tension.
    impure subroutine s_populate_capillary_buffers(c_divs, bc_type, bc)

        type(scalar_field), dimension(num_dims + 1), intent(inout) :: c_divs
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        type(bc_xyz_info), intent(in)                              :: bc
        integer                                                    :: k, l

        !> x-direction

        if (bc%x%beg >= 0) then
            call s_mpi_sendrecv_variables_buffers(c_divs, 1, -1, num_dims + 1)
        else
            $:GPU_PARALLEL_LOOP(private='[l, k]', collapse=2)
            do l = 0, p
                do k = 0, n
                    select case (bc_type(1, 1)%sf(0, k, l))
                    case (BC_PERIODIC)
                        call s_color_function_periodic(c_divs, 1, -1, k, l)
                    case (BC_REFLECTIVE)
                        call s_color_function_reflective(c_divs, 1, -1, k, l)
                    case default
                        call s_color_function_ghost_cell_extrapolation(c_divs, 1, -1, k, l)
                    end select
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        if (bc%x%end >= 0) then
            call s_mpi_sendrecv_variables_buffers(c_divs, 1, 1, num_dims + 1)
        else
            $:GPU_PARALLEL_LOOP(private='[l, k]', collapse=2)
            do l = 0, p
                do k = 0, n
                    select case (bc_type(1, 2)%sf(0, k, l))
                    case (BC_PERIODIC)
                        call s_color_function_periodic(c_divs, 1, 1, k, l)
                    case (BC_REFLECTIVE)
                        call s_color_function_reflective(c_divs, 1, 1, k, l)
                    case default
                        call s_color_function_ghost_cell_extrapolation(c_divs, 1, 1, k, l)
                    end select
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        if (n == 0) return

        #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
            !> y-direction
            if (bc%y%beg >= 0) then
                call s_mpi_sendrecv_variables_buffers(c_divs, 2, -1, num_dims + 1)
            else
                $:GPU_PARALLEL_LOOP(private='[l, k]', collapse=2)
                do l = 0, p
                    do k = -buff_size, m + buff_size
                        select case (bc_type(2, 1)%sf(k, 0, l))
                        case (BC_PERIODIC)
                            call s_color_function_periodic(c_divs, 2, -1, k, l)
                        case (BC_REFLECTIVE)
                            call s_color_function_reflective(c_divs, 2, -1, k, l)
                        case default
                            call s_color_function_ghost_cell_extrapolation(c_divs, 2, -1, k, l)
                        end select
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            if (bc%y%end >= 0) then
                call s_mpi_sendrecv_variables_buffers(c_divs, 2, 1, num_dims + 1)
            else
                $:GPU_PARALLEL_LOOP(private='[l, k]', collapse=2)
                do l = 0, p
                    do k = -buff_size, m + buff_size
                        select case (bc_type(2, 2)%sf(k, 0, l))
                        case (BC_PERIODIC)
                            call s_color_function_periodic(c_divs, 2, 1, k, l)
                        case (BC_REFLECTIVE)
                            call s_color_function_reflective(c_divs, 2, 1, k, l)
                        case default
                            call s_color_function_ghost_cell_extrapolation(c_divs, 2, 1, k, l)
                        end select
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endif

        if (p == 0) return

        #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
            !> z-direction
            if (bc%z%beg >= 0) then
                call s_mpi_sendrecv_variables_buffers(c_divs, 3, -1, num_dims + 1)
            else
                $:GPU_PARALLEL_LOOP(private='[l, k]', collapse=2)
                do l = -buff_size, n + buff_size
                    do k = -buff_size, m + buff_size
                        select case (bc_type(3, 1)%sf(k, l, 0))
                        case (BC_PERIODIC)
                            call s_color_function_periodic(c_divs, 3, -1, k, l)
                        case (BC_REFLECTIVE)
                            call s_color_function_reflective(c_divs, 3, -1, k, l)
                        case default
                            call s_color_function_ghost_cell_extrapolation(c_divs, 3, -1, k, l)
                        end select
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            if (bc%z%end >= 0) then
                call s_mpi_sendrecv_variables_buffers(c_divs, 3, 1, num_dims + 1)
            else
                $:GPU_PARALLEL_LOOP(private='[l, k]', collapse=2)
                do l = -buff_size, n + buff_size
                    do k = -buff_size, m + buff_size
                        select case (bc_type(3, 2)%sf(k, l, 0))
                        case (BC_PERIODIC)
                            call s_color_function_periodic(c_divs, 3, 1, k, l)
                        case (BC_REFLECTIVE)
                            call s_color_function_reflective(c_divs, 3, 1, k, l)
                        case default
                            call s_color_function_ghost_cell_extrapolation(c_divs, 3, 1, k, l)
                        end select
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endif

    end subroutine s_populate_capillary_buffers

    !> Apply periodic boundary conditions to the color function and its divergence fields.
    subroutine s_color_function_periodic(c_divs, bc_dir, bc_loc, k, l)

        $:GPU_ROUTINE(function_name='s_color_function_periodic', parallelism='[seq]', cray_inline=True)
        type(scalar_field), dimension(num_dims + 1), intent(inout) :: c_divs
        integer, intent(in)                                        :: bc_dir, bc_loc
        integer, intent(in)                                        :: k, l
        integer                                                    :: j, i

        if (bc_dir == 1) then  !< x-direction
            if (bc_loc == -1) then  ! bc_x%beg
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        c_divs(i)%sf(-j, k, l) = c_divs(i)%sf(m - (j - 1), k, l)
                    end do
                end do
            else  !< bc_x%end
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        c_divs(i)%sf(m + j, k, l) = c_divs(i)%sf(j - 1, k, l)
                    end do
                end do
            end if
        else if (bc_dir == 2) then  !< y-direction
            if (bc_loc == -1) then  !< bc_y%beg
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        c_divs(i)%sf(k, -j, l) = c_divs(i)%sf(k, n - (j - 1), l)
                    end do
                end do
            else  !< bc_y%end
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        c_divs(i)%sf(k, n + j, l) = c_divs(i)%sf(k, j - 1, l)
                    end do
                end do
            end if
        else if (bc_dir == 3) then  !< z-direction
            if (bc_loc == -1) then  !< bc_z%beg
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        c_divs(i)%sf(k, l, -j) = c_divs(i)%sf(k, l, p - (j - 1))
                    end do
                end do
            else  !< bc_z%end
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        c_divs(i)%sf(k, l, p + j) = c_divs(i)%sf(k, l, j - 1)
                    end do
                end do
            end if
        end if

    end subroutine s_color_function_periodic

    !> Apply reflective boundary conditions to the color function and its divergence fields.
    subroutine s_color_function_reflective(c_divs, bc_dir, bc_loc, k, l)

        $:GPU_ROUTINE(function_name='s_color_function_reflective', parallelism='[seq]', cray_inline=True)
        type(scalar_field), dimension(num_dims + 1), intent(inout) :: c_divs
        integer, intent(in)                                        :: bc_dir, bc_loc
        integer, intent(in)                                        :: k, l
        integer                                                    :: j, i

        if (bc_dir == 1) then  !< x-direction
            if (bc_loc == -1) then  ! bc_x%beg
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        if (i == bc_dir) then
                            c_divs(i)%sf(-j, k, l) = -c_divs(i)%sf(j - 1, k, l)
                        else
                            c_divs(i)%sf(-j, k, l) = c_divs(i)%sf(j - 1, k, l)
                        end if
                    end do
                end do
            else  !< bc_x%end
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        if (i == bc_dir) then
                            c_divs(i)%sf(m + j, k, l) = -c_divs(i)%sf(m - (j - 1), k, l)
                        else
                            c_divs(i)%sf(m + j, k, l) = c_divs(i)%sf(m - (j - 1), k, l)
                        end if
                    end do
                end do
            end if
        else if (bc_dir == 2) then  !< y-direction
            if (bc_loc == -1) then  !< bc_y%beg
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        if (i == bc_dir) then
                            c_divs(i)%sf(k, -j, l) = -c_divs(i)%sf(k, j - 1, l)
                        else
                            c_divs(i)%sf(k, -j, l) = c_divs(i)%sf(k, j - 1, l)
                        end if
                    end do
                end do
            else  !< bc_y%end
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        if (i == bc_dir) then
                            c_divs(i)%sf(k, n + j, l) = -c_divs(i)%sf(k, n - (j - 1), l)
                        else
                            c_divs(i)%sf(k, n + j, l) = c_divs(i)%sf(k, n - (j - 1), l)
                        end if
                    end do
                end do
            end if
        else if (bc_dir == 3) then  !< z-direction
            if (bc_loc == -1) then  !< bc_z%beg
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        if (i == bc_dir) then
                            c_divs(i)%sf(k, l, -j) = -c_divs(i)%sf(k, l, j - 1)
                        else
                            c_divs(i)%sf(k, l, -j) = c_divs(i)%sf(k, l, j - 1)
                        end if
                    end do
                end do
            else  !< bc_z%end
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        if (i == bc_dir) then
                            c_divs(i)%sf(k, l, p + j) = -c_divs(i)%sf(k, l, p - (j - 1))
                        else
                            c_divs(i)%sf(k, l, p + j) = c_divs(i)%sf(k, l, p - (j - 1))
                        end if
                    end do
                end do
            end if
        end if

    end subroutine s_color_function_reflective

    !> Extrapolate the color function and its divergence into ghost cells by copying boundary values.
    subroutine s_color_function_ghost_cell_extrapolation(c_divs, bc_dir, bc_loc, k, l)

        $:GPU_ROUTINE(function_name='s_color_function_ghost_cell_extrapolation', parallelism='[seq]', cray_inline=True)
        type(scalar_field), dimension(num_dims + 1), intent(inout) :: c_divs
        integer, intent(in)                                        :: bc_dir, bc_loc
        integer, intent(in)                                        :: k, l
        integer                                                    :: j, i

        if (bc_dir == 1) then  !< x-direction
            if (bc_loc == -1) then  ! bc_x%beg
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        c_divs(i)%sf(-j, k, l) = c_divs(i)%sf(0, k, l)
                    end do
                end do
            else  !< bc_x%end
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        c_divs(i)%sf(m + j, k, l) = c_divs(i)%sf(m, k, l)
                    end do
                end do
            end if
        else if (bc_dir == 2) then  !< y-direction
            if (bc_loc == -1) then  !< bc_y%beg
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        c_divs(i)%sf(k, -j, l) = c_divs(i)%sf(k, 0, l)
                    end do
                end do
            else  !< bc_y%end
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        c_divs(i)%sf(k, n + j, l) = c_divs(i)%sf(k, n, l)
                    end do
                end do
            end if
        else if (bc_dir == 3) then  !< z-direction
            if (bc_loc == -1) then  !< bc_z%beg
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        c_divs(i)%sf(k, l, -j) = c_divs(i)%sf(k, l, 0)
                    end do
                end do
            else  !< bc_z%end
                do i = 1, num_dims + 1
                    do j = 1, buff_size
                        c_divs(i)%sf(k, l, p + j) = c_divs(i)%sf(k, l, p)
                    end do
                end do
            end if
        end if

    end subroutine s_color_function_ghost_cell_extrapolation

    !> Populate ghost cell buffers for the Jacobian scalar field used in the IGR elliptic solver.
    impure subroutine s_populate_F_igr_buffers(bc_type, jac_sf)

        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        type(scalar_field), dimension(1:), intent(inout)           :: jac_sf
        integer                                                    :: j, k, l

        if (bc_x%beg >= 0) then
            call s_mpi_sendrecv_variables_buffers(jac_sf, 1, -1, 1)
        else
            $:GPU_PARALLEL_LOOP(private='[l, k]', collapse=2)
            do l = 0, p
                do k = 0, n
                    select case (bc_type(1, 1)%sf(0, k, l))
                    case (BC_PERIODIC)
                        do j = 1, buff_size
                            jac_sf(1)%sf(-j, k, l) = jac_sf(1)%sf(m - j + 1, k, l)
                        end do
                    case (BC_REFLECTIVE)
                        do j = 1, buff_size
                            jac_sf(1)%sf(-j, k, l) = jac_sf(1)%sf(j - 1, k, l)
                        end do
                    case default
                        do j = 1, buff_size
                            jac_sf(1)%sf(-j, k, l) = jac_sf(1)%sf(0, k, l)
                        end do
                    end select
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        if (bc_x%end >= 0) then
            call s_mpi_sendrecv_variables_buffers(jac_sf, 1, 1, 1)
        else
            $:GPU_PARALLEL_LOOP(private='[l, k]', collapse=2)
            do l = 0, p
                do k = 0, n
                    select case (bc_type(1, 2)%sf(0, k, l))
                    case (BC_PERIODIC)
                        do j = 1, buff_size
                            jac_sf(1)%sf(m + j, k, l) = jac_sf(1)%sf(j - 1, k, l)
                        end do
                    case (BC_REFLECTIVE)
                        do j = 1, buff_size
                            jac_sf(1)%sf(m + j, k, l) = jac_sf(1)%sf(m - (j - 1), k, l)
                        end do
                    case default
                        do j = 1, buff_size
                            jac_sf(1)%sf(m + j, k, l) = jac_sf(1)%sf(m, k, l)
                        end do
                    end select
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
            if (n == 0) then
                return
            else if (bc_y%beg >= 0) then
                call s_mpi_sendrecv_variables_buffers(jac_sf, 2, -1, 1)
            else
                $:GPU_PARALLEL_LOOP(private='[l, k]', collapse=2)
                do l = 0, p
                    do k = idwbuff(1)%beg, idwbuff(1)%end
                        select case (bc_type(2, 1)%sf(k, 0, l))
                        case (BC_PERIODIC)
                            do j = 1, buff_size
                                jac_sf(1)%sf(k, -j, l) = jac_sf(1)%sf(k, n - j + 1, l)
                            end do
                        case (BC_REFLECTIVE)
                            do j = 1, buff_size
                                jac_sf(1)%sf(k, -j, l) = jac_sf(1)%sf(k, j - 1, l)
                            end do
                        case default
                            do j = 1, buff_size
                                jac_sf(1)%sf(k, -j, l) = jac_sf(1)%sf(k, 0, l)
                            end do
                        end select
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            if (bc_y%end >= 0) then
                call s_mpi_sendrecv_variables_buffers(jac_sf, 2, 1, 1)
            else
                $:GPU_PARALLEL_LOOP(private='[l, k]', collapse=2)
                do l = 0, p
                    do k = idwbuff(1)%beg, idwbuff(1)%end
                        select case (bc_type(2, 2)%sf(k, 0, l))
                        case (BC_PERIODIC)
                            do j = 1, buff_size
                                jac_sf(1)%sf(k, n + j, l) = jac_sf(1)%sf(k, j - 1, l)
                            end do
                        case (BC_REFLECTIVE)
                            do j = 1, buff_size
                                jac_sf(1)%sf(k, n + j, l) = jac_sf(1)%sf(k, n - (j - 1), l)
                            end do
                        case default
                            do j = 1, buff_size
                                jac_sf(1)%sf(k, n + j, l) = jac_sf(1)%sf(k, n, l)
                            end do
                        end select
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endif

        #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
            if (p == 0) then
                return
            else if (bc_z%beg >= 0) then
                call s_mpi_sendrecv_variables_buffers(jac_sf, 3, -1, 1)
            else
                $:GPU_PARALLEL_LOOP(private='[l, k]', collapse=2)
                do l = idwbuff(2)%beg, idwbuff(2)%end
                    do k = idwbuff(1)%beg, idwbuff(1)%end
                        select case (bc_type(3, 1)%sf(k, l, 0))
                        case (BC_PERIODIC)
                            do j = 1, buff_size
                                jac_sf(1)%sf(k, l, -j) = jac_sf(1)%sf(k, l, p - j + 1)
                            end do
                        case (BC_REFLECTIVE)
                            do j = 1, buff_size
                                jac_sf(1)%sf(k, l, -j) = jac_sf(1)%sf(k, l, j - 1)
                            end do
                        case default
                            do j = 1, buff_size
                                jac_sf(1)%sf(k, l, -j) = jac_sf(1)%sf(k, l, 0)
                            end do
                        end select
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            if (bc_z%end >= 0) then
                call s_mpi_sendrecv_variables_buffers(jac_sf, 3, 1, 1)
            else
                $:GPU_PARALLEL_LOOP(private='[l, k]', collapse=2)
                do l = idwbuff(2)%beg, idwbuff(2)%end
                    do k = idwbuff(1)%beg, idwbuff(1)%end
                        select case (bc_type(3, 2)%sf(k, l, 0))
                        case (BC_PERIODIC)
                            do j = 1, buff_size
                                jac_sf(1)%sf(k, l, p + j) = jac_sf(1)%sf(k, l, j - 1)
                            end do
                        case (BC_REFLECTIVE)
                            do j = 1, buff_size
                                jac_sf(1)%sf(k, l, p + j) = jac_sf(1)%sf(k, l, p - (j - 1))
                            end do
                        case default
                            do j = 1, buff_size
                                jac_sf(1)%sf(k, l, p + j) = jac_sf(1)%sf(k, l, p)
                            end do
                        end select
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endif

    end subroutine s_populate_F_igr_buffers

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

    !> Populate the buffers of the grid variables, which are constituted of the cell-boundary locations and cell-width
    !! distributions, based on the boundary conditions.
    subroutine s_populate_grid_variables_buffers

        integer :: i

#ifdef MFC_SIMULATION
        ! Required for compatibility between codes
        type(int_bounds_info) :: offset_x, offset_y, offset_z

        offset_x%beg = buff_size; offset_x%end = buff_size
        offset_y%beg = buff_size; offset_y%end = buff_size
        offset_z%beg = buff_size; offset_z%end = buff_size

        ! Global domain bounds
#ifdef MFC_MPI
        call s_mpi_allreduce_min(x_cb(-1), glb_bounds(1)%beg)
        call s_mpi_allreduce_max(x_cb(m), glb_bounds(1)%end)
        if (n > 0) then
            call s_mpi_allreduce_min(y_cb(-1), glb_bounds(2)%beg)
            call s_mpi_allreduce_max(y_cb(n), glb_bounds(2)%end)
            if (p > 0) then
                call s_mpi_allreduce_min(z_cb(-1), glb_bounds(3)%beg)
                call s_mpi_allreduce_max(z_cb(p), glb_bounds(3)%end)
            end if
        end if
#else
        glb_bounds(1)%beg = x_cb(-1); glb_bounds(1)%end = x_cb(m)
        if (n > 0) then
            glb_bounds(2)%beg = y_cb(-1); glb_bounds(2)%end = y_cb(n)
            if (p > 0) then
                glb_bounds(3)%beg = z_cb(-1); glb_bounds(3)%end = z_cb(p)
            end if
        end if
#endif
        $:GPU_UPDATE(device='[glb_bounds]')
#endif

#ifndef MFC_PRE_PROCESS
        ! Population of Buffers in x-direction

        ! Populating cell-width distribution buffer at bc_x%beg
        if (bc_x%beg >= 0) then
            call s_mpi_sendrecv_grid_variables_buffers(1, -1)
        else if (bc_x%beg <= BC_GHOST_EXTRAP) then
            do i = 1, buff_size
                dx(-i) = dx(0)
            end do
        else if (bc_x%beg == BC_REFLECTIVE) then
            do i = 1, buff_size
                dx(-i) = dx(i - 1)
            end do
        else if (bc_x%beg == BC_PERIODIC) then
            do i = 1, buff_size
                dx(-i) = dx(m - (i - 1))
            end do
        end if

        ! Computing the cell-boundary and center locations buffer at bc_x%beg
        do i = 1, offset_x%beg
            x_cb(-1 - i) = x_cb(-i) - dx(-i)
        end do

        do i = 1, buff_size
            x_cc(-i) = x_cc(1 - i) - (dx(1 - i) + dx(-i))/2._wp
        end do

        ! Populating the cell-width distribution buffer at bc_x%end
        if (bc_x%end >= 0) then
            call s_mpi_sendrecv_grid_variables_buffers(1, 1)
        else if (bc_x%end <= BC_GHOST_EXTRAP) then
            do i = 1, buff_size
                dx(m + i) = dx(m)
            end do
        else if (bc_x%end == BC_REFLECTIVE) then
            do i = 1, buff_size
                dx(m + i) = dx(m - (i - 1))
            end do
        else if (bc_x%end == BC_PERIODIC) then
            do i = 1, buff_size
                dx(m + i) = dx(i - 1)
            end do
        end if

        ! Populating the cell-boundary and center locations buffer at bc_x%end
        do i = 1, offset_x%end
            x_cb(m + i) = x_cb(m + (i - 1)) + dx(m + i)
        end do

        do i = 1, buff_size
            x_cc(m + i) = x_cc(m + (i - 1)) + (dx(m + (i - 1)) + dx(m + i))/2._wp
        end do

        ! Population of Buffers in y-direction

        ! Populating cell-width distribution buffer at bc_y%beg
        if (n == 0) then
            return
        else if (bc_y%beg >= 0) then
            call s_mpi_sendrecv_grid_variables_buffers(2, -1)
        else if (bc_y%beg <= BC_GHOST_EXTRAP .and. bc_y%beg /= BC_AXIS) then
            do i = 1, buff_size
                dy(-i) = dy(0)
            end do
        else if (bc_y%beg == BC_REFLECTIVE .or. bc_y%beg == BC_AXIS) then
            do i = 1, buff_size
                dy(-i) = dy(i - 1)
            end do
        else if (bc_y%beg == BC_PERIODIC) then
            do i = 1, buff_size
                dy(-i) = dy(n - (i - 1))
            end do
        end if

        ! Computing the cell-boundary and center locations buffer at bc_y%beg
        do i = 1, offset_y%beg
            y_cb(-1 - i) = y_cb(-i) - dy(-i)
        end do

        do i = 1, buff_size
            y_cc(-i) = y_cc(1 - i) - (dy(1 - i) + dy(-i))/2._wp
        end do

        ! Populating the cell-width distribution buffer at bc_y%end
        if (bc_y%end >= 0) then
            call s_mpi_sendrecv_grid_variables_buffers(2, 1)
        else if (bc_y%end <= BC_GHOST_EXTRAP) then
            do i = 1, buff_size
                dy(n + i) = dy(n)
            end do
        else if (bc_y%end == BC_REFLECTIVE) then
            do i = 1, buff_size
                dy(n + i) = dy(n - (i - 1))
            end do
        else if (bc_y%end == BC_PERIODIC) then
            do i = 1, buff_size
                dy(n + i) = dy(i - 1)
            end do
        end if

        ! Populating the cell-boundary and center locations buffer at bc_y%end
        do i = 1, offset_y%end
            y_cb(n + i) = y_cb(n + (i - 1)) + dy(n + i)
        end do

        do i = 1, buff_size
            y_cc(n + i) = y_cc(n + (i - 1)) + (dy(n + (i - 1)) + dy(n + i))/2._wp
        end do

        ! Population of Buffers in z-direction

        ! Populating cell-width distribution buffer at bc_z%beg
        if (p == 0) then
            return
        else if (Bc_z%beg >= 0) then
            call s_mpi_sendrecv_grid_variables_buffers(3, -1)
        else if (bc_z%beg <= BC_GHOST_EXTRAP) then
            do i = 1, buff_size
                dz(-i) = dz(0)
            end do
        else if (bc_z%beg == BC_REFLECTIVE) then
            do i = 1, buff_size
                dz(-i) = dz(i - 1)
            end do
        else if (bc_z%beg == BC_PERIODIC) then
            do i = 1, buff_size
                dz(-i) = dz(p - (i - 1))
            end do
        end if

        ! Computing the cell-boundary and center locations buffer at bc_z%beg
        do i = 1, offset_z%beg
            z_cb(-1 - i) = z_cb(-i) - dz(-i)
        end do

        do i = 1, buff_size
            z_cc(-i) = z_cc(1 - i) - (dz(1 - i) + dz(-i))/2._wp
        end do

        ! Populating the cell-width distribution buffer at bc_z%end
        if (bc_z%end >= 0) then
            call s_mpi_sendrecv_grid_variables_buffers(3, 1)
        else if (bc_z%end <= BC_GHOST_EXTRAP) then
            do i = 1, buff_size
                dz(p + i) = dz(p)
            end do
        else if (bc_z%end == BC_REFLECTIVE) then
            do i = 1, buff_size
                dz(p + i) = dz(p - (i - 1))
            end do
        else if (bc_z%end == BC_PERIODIC) then
            do i = 1, buff_size
                dz(p + i) = dz(i - 1)
            end do
        end if

        ! Populating the cell-boundary and center locations buffer at bc_z%end
        do i = 1, offset_z%end
            z_cb(p + i) = z_cb(p + (i - 1)) + dz(p + i)
        end do

        do i = 1, buff_size
            z_cc(p + i) = z_cc(p + (i - 1)) + (dz(p + (i - 1)) + dz(p + i))/2._wp
        end do
#endif

    end subroutine s_populate_grid_variables_buffers

end module m_boundary_io
