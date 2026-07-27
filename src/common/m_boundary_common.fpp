!>
!! @file
!! @brief Contains module m_boundary_common

!> @brief Noncharacteristic and processor boundary condition application for ghost cells and buffer regions
#:include 'case.fpp'
#:include 'macros.fpp'

module m_boundary_common

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_mpi_common
    use m_constants
    use m_boundary_primitives
    use m_boundary_io

    implicit none

    private; public :: s_initialize_boundary_common_module, s_populate_variables_buffers, s_populate_capillary_buffers, &
        & s_populate_F_igr_buffers, s_populate_grid_variables_buffers, s_finalize_boundary_common_module, s_populate_beta_buffers

    public :: bc_buffers

#ifdef MFC_MPI
    public :: MPI_BC_TYPE_TYPE, MPI_BC_BUFFER_TYPE
#endif

    !> Lagrangian-bubble beta (void-fraction) buffer bounds (#1290)
    type(int_bounds_info), dimension(3) :: beta_bc_bounds
    $:GPU_DECLARE(create='[beta_bc_bounds]')

contains

    !> Allocate and set up boundary condition buffer arrays for all coordinate directions.
    impure subroutine s_initialize_boundary_common_module(use_dirichlet_buffers)

        integer                       :: i, j, sys_size_alloc
        logical, intent(in), optional :: use_dirichlet_buffers

        dirichlet_from_buffers = .false.
        if (present(use_dirichlet_buffers)) dirichlet_from_buffers = use_dirichlet_buffers
        $:GPU_UPDATE(device='[dirichlet_from_buffers]')

        @:ALLOCATE(bc_buffers(1:3, 1:2))

        if (bc_io) then
            sys_size_alloc = sys_size
            if (chemistry) sys_size_alloc = sys_size + 1

            @:ALLOCATE(bc_buffers(1, 1)%sf(1:sys_size_alloc, 0:n, 0:p))
            @:ALLOCATE(bc_buffers(1, 2)%sf(1:sys_size_alloc, 0:n, 0:p))
            #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
                if (n > 0) then
                    @:ALLOCATE(bc_buffers(2,1)%sf(-buff_size:m+buff_size,1:sys_size_alloc,0:p))
                    @:ALLOCATE(bc_buffers(2,2)%sf(-buff_size:m+buff_size,1:sys_size_alloc,0:p))
                    #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                        if (p > 0) then
                            @:ALLOCATE(bc_buffers(3,1)%sf(-buff_size:m+buff_size,-buff_size:n+buff_size,1:sys_size_alloc))
                            @:ALLOCATE(bc_buffers(3,2)%sf(-buff_size:m+buff_size,-buff_size:n+buff_size,1:sys_size_alloc))
                        end if
                    #:endif
                end if
            #:endif
            do i = 1, num_dims
                do j = 1, 2
                    @:ACC_SETUP_SFs(bc_buffers(i,j))
                end do
            end do
        end if

        if (bubbles_lagrange) then
            beta_bc_bounds(1)%beg = -mapcells - 1
            beta_bc_bounds(1)%end = m + mapcells + 1
            ! n > 0 always for bubbles_lagrange
            beta_bc_bounds(2)%beg = -mapcells - 1
            beta_bc_bounds(2)%end = n + mapcells + 1
            if (p == 0) then
                beta_bc_bounds(3)%beg = 0
                beta_bc_bounds(3)%end = 0
            else
                beta_bc_bounds(3)%beg = -mapcells - 1
                beta_bc_bounds(3)%end = p + mapcells + 1
            end if
        end if
        $:GPU_UPDATE(device='[beta_bc_bounds]')

    end subroutine s_initialize_boundary_common_module

    !> Populate the buffers of the primitive variables based on the selected boundary conditions.
    impure subroutine s_populate_variables_buffers(bc_type, q_prim_vf, pb_in, mv_in, q_T_sf)

        type(scalar_field), dimension(sys_size), intent(inout)                                               :: q_prim_vf
        real(stp), optional, dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_in, mv_in
        type(integer_field), dimension(1:num_dims,1:2), intent(in)                                           :: bc_type
        type(scalar_field), optional, intent(inout)                                                          :: q_T_sf

        call s_populate_bc_direction(1, -1, bc_x, bc_type(1, 1), q_prim_vf, pb_in, mv_in, q_T_sf)
        call s_populate_bc_direction(1, 1, bc_x, bc_type(1, 2), q_prim_vf, pb_in, mv_in, q_T_sf)

        if (n == 0) return

        #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
            call s_populate_bc_direction(2, -1, bc_y, bc_type(2, 1), q_prim_vf, pb_in, mv_in, q_T_sf)
            call s_populate_bc_direction(2, 1, bc_y, bc_type(2, 2), q_prim_vf, pb_in, mv_in, q_T_sf)
        #:endif

        if (p == 0) return

        #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
            call s_populate_bc_direction(3, -1, bc_z, bc_type(3, 1), q_prim_vf, pb_in, mv_in, q_T_sf)
            call s_populate_bc_direction(3, 1, bc_z, bc_type(3, 2), q_prim_vf, pb_in, mv_in, q_T_sf)
        #:endif

    end subroutine s_populate_variables_buffers

    !> Populate the variable buffers along one direction and location, via MPI exchange for processor boundaries or by dispatching
    !! the per-cell BC routines over the boundary face.
    impure subroutine s_populate_bc_direction(bc_dir, bc_loc, bc_bounds, bc_type_edge, q_prim_vf, pb_in, mv_in, q_T_sf)

        integer, intent(in) :: bc_dir, bc_loc
        type(int_bounds_info), intent(in) :: bc_bounds
        type(integer_field), intent(in) :: bc_type_edge
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(stp), optional, dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_in, mv_in
        type(scalar_field), optional, intent(inout) :: q_T_sf
        integer :: bc_edge, k_beg, k_end, l_beg, l_end
        integer :: bc_code, k, l

        if (bc_loc == -1) then
            bc_edge = bc_bounds%beg
        else
            bc_edge = bc_bounds%end
        end if

        ! BC type codes defined in m_constants.fpp; non-negative values are MPI boundaries
        if (bc_edge >= 0) then
            call s_mpi_sendrecv_variables_buffers(q_prim_vf, bc_dir, bc_loc, sys_size, pb_in, mv_in, q_T_sf)
            return
        end if

        if (bc_dir == 1) then
            k_beg = 0; k_end = n; l_beg = 0; l_end = p
        else if (bc_dir == 2) then
            k_beg = -buff_size; k_end = m + buff_size; l_beg = 0; l_end = p
        else
            k_beg = -buff_size; k_end = m + buff_size; l_beg = -buff_size; l_end = n + buff_size
        end if

        $:GPU_PARALLEL_LOOP(private='[l, k, bc_code]', collapse=2)
        do l = l_beg, l_end
            do k = k_beg, k_end
                if (bc_dir == 1) then
                    bc_code = int(bc_type_edge%sf(0, k, l))
                else if (bc_dir == 2) then
                    bc_code = int(bc_type_edge%sf(k, 0, l))
                else
                    bc_code = int(bc_type_edge%sf(k, l, 0))
                end if

                select case (bc_code)
                case (BC_CHAR_SUP_OUTFLOW:BC_GHOST_EXTRAP)
                    call s_ghost_cell_extrapolation(q_prim_vf, bc_dir, bc_loc, k, l, q_T_sf)
                case (BC_AXIS)
                    if (bc_dir == 2 .and. bc_loc == -1) call s_axis(q_prim_vf, pb_in, mv_in, k, l)
                case (BC_REFLECTIVE)
                    call s_symmetry(q_prim_vf, bc_dir, bc_loc, k, l, pb_in, mv_in, q_T_sf)
                case (BC_PERIODIC)
                    call s_periodic(q_prim_vf, bc_dir, bc_loc, k, l, pb_in, mv_in, q_T_sf)
                case (BC_SLIP_WALL)
                    call s_slip_wall(q_prim_vf, bc_dir, bc_loc, k, l, q_T_sf)
                case (BC_NO_SLIP_WALL)
                    call s_no_slip_wall(q_prim_vf, bc_dir, bc_loc, k, l, q_T_sf)
                case (BC_DIRICHLET)
                    call s_dirichlet(q_prim_vf, bc_dir, bc_loc, k, l, q_T_sf)
                end select

                if (qbmm .and. (.not. polytropic) .and. present(pb_in) .and. present(mv_in) .and. (bc_code <= BC_GHOST_EXTRAP) &
                    & .and. .not. (bc_dir == 2 .and. bc_loc == -1 .and. bc_code == BC_AXIS)) then
                    call s_qbmm_extrapolation(bc_dir, bc_loc, k, l, pb_in, mv_in)
                end if
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_populate_bc_direction

    !> Populate ghost cell buffers for the color function and its divergence used in capillary surface tension.
    impure subroutine s_populate_capillary_buffers(c_divs, bc_type, bc)

        type(scalar_field), dimension(num_dims + 1), intent(inout) :: c_divs
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        type(bc_xyz_info), intent(in)                              :: bc

        call s_populate_capillary_bc_direction(1, -1, bc%x, bc_type(1, 1), c_divs)
        call s_populate_capillary_bc_direction(1, 1, bc%x, bc_type(1, 2), c_divs)

        if (n == 0) return

        #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
            call s_populate_capillary_bc_direction(2, -1, bc%y, bc_type(2, 1), c_divs)
            call s_populate_capillary_bc_direction(2, 1, bc%y, bc_type(2, 2), c_divs)
        #:endif

        if (p == 0) return

        #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
            call s_populate_capillary_bc_direction(3, -1, bc%z, bc_type(3, 1), c_divs)
            call s_populate_capillary_bc_direction(3, 1, bc%z, bc_type(3, 2), c_divs)
        #:endif

    end subroutine s_populate_capillary_buffers

    !> Populate ghost cell buffers for one capillary BC direction and location, via MPI exchange for processor boundaries or by
    !! dispatching the per-cell capillary BC routines over the boundary face.
    impure subroutine s_populate_capillary_bc_direction(bc_dir, bc_loc, bc_bounds, bc_type_edge, c_divs)

        integer, intent(in)                                        :: bc_dir, bc_loc
        type(int_bounds_info), intent(in)                          :: bc_bounds
        type(scalar_field), dimension(num_dims + 1), intent(inout) :: c_divs
        type(integer_field), intent(in)                            :: bc_type_edge
        integer                                                    :: bc_edge, k_beg, k_end, l_beg, l_end, k, l, bc_code

        if (bc_loc == -1) then
            bc_edge = bc_bounds%beg
        else
            bc_edge = bc_bounds%end
        end if

        if (bc_edge >= 0) then
            call s_mpi_sendrecv_variables_buffers(c_divs, bc_dir, bc_loc, num_dims + 1)
            return
        end if

        if (bc_dir == 1) then
            k_beg = 0; k_end = n; l_beg = 0; l_end = p
        else if (bc_dir == 2) then
            k_beg = -buff_size; k_end = m + buff_size; l_beg = 0; l_end = p
        else
            k_beg = -buff_size; k_end = m + buff_size; l_beg = -buff_size; l_end = n + buff_size
        end if

        $:GPU_PARALLEL_LOOP(private='[l, k, bc_code]', collapse=2)
        do l = l_beg, l_end
            do k = k_beg, k_end
                if (bc_dir == 1) then
                    bc_code = int(bc_type_edge%sf(0, k, l))
                else if (bc_dir == 2) then
                    bc_code = int(bc_type_edge%sf(k, 0, l))
                else
                    bc_code = int(bc_type_edge%sf(k, l, 0))
                end if

                select case (bc_code)
                case (BC_PERIODIC)
                    call s_color_function_periodic(c_divs, bc_dir, bc_loc, k, l)
                case (BC_REFLECTIVE)
                    call s_color_function_reflective(c_divs, bc_dir, bc_loc, k, l)
                case default
                    call s_color_function_ghost_cell_extrapolation(c_divs, bc_dir, bc_loc, k, l)
                end select
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_populate_capillary_bc_direction

    !> Populate ghost cell buffers for the Jacobian scalar field used in the IGR elliptic solver.
    impure subroutine s_populate_F_igr_buffers(bc_type, jac_sf)

        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        type(scalar_field), dimension(1:), intent(inout)           :: jac_sf

        call s_populate_F_igr_bc_direction(1, -1, bc_x, bc_type(1, 1), jac_sf)
        call s_populate_F_igr_bc_direction(1, 1, bc_x, bc_type(1, 2), jac_sf)

        if (n == 0) return

        #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
            call s_populate_F_igr_bc_direction(2, -1, bc_y, bc_type(2, 1), jac_sf)
            call s_populate_F_igr_bc_direction(2, 1, bc_y, bc_type(2, 2), jac_sf)
        #:endif

        if (p == 0) return

        #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
            call s_populate_F_igr_bc_direction(3, -1, bc_z, bc_type(3, 1), jac_sf)
            call s_populate_F_igr_bc_direction(3, 1, bc_z, bc_type(3, 2), jac_sf)
        #:endif

    end subroutine s_populate_F_igr_buffers

    !> Populate ghost cell buffers for one IGR Jacobian BC direction and location, via MPI exchange for processor boundaries or by
    !! dispatching the per-cell IGR Jacobian BC routines over the boundary face.
    impure subroutine s_populate_F_igr_bc_direction(bc_dir, bc_loc, bc_bounds, bc_type_edge, jac_sf)

        integer, intent(in)                              :: bc_dir, bc_loc
        type(int_bounds_info), intent(in)                :: bc_bounds
        type(integer_field), intent(in)                  :: bc_type_edge
        type(scalar_field), dimension(1:), intent(inout) :: jac_sf
        integer                                          :: bc_edge, k_beg, k_end, l_beg, l_end, k, l, j, bc_code

        if (bc_loc == -1) then
            bc_edge = bc_bounds%beg
        else
            bc_edge = bc_bounds%end
        end if

        if (bc_edge >= 0) then
            call s_mpi_sendrecv_variables_buffers(jac_sf, bc_dir, bc_loc, 1)
            return
        end if

        if (bc_dir == 1) then
            k_beg = 0; k_end = n; l_beg = 0; l_end = p
        else if (bc_dir == 2) then
            k_beg = idwbuff(1)%beg; k_end = idwbuff(1)%end; l_beg = 0; l_end = p
        else
            k_beg = idwbuff(1)%beg; k_end = idwbuff(1)%end; l_beg = idwbuff(2)%beg; l_end = idwbuff(2)%end
        end if

        $:GPU_PARALLEL_LOOP(private='[l, k, bc_code]', collapse=2)
        do l = l_beg, l_end
            do k = k_beg, k_end
                if (bc_dir == 1) then
                    bc_code = int(bc_type_edge%sf(0, k, l))
                else if (bc_dir == 2) then
                    bc_code = int(bc_type_edge%sf(k, 0, l))
                else
                    bc_code = int(bc_type_edge%sf(k, l, 0))
                end if

                select case (bc_code)
                case (BC_PERIODIC)
                    call s_F_igr_periodic(jac_sf, bc_dir, bc_loc, k, l)
                case (BC_REFLECTIVE)
                    call s_F_igr_reflective(jac_sf, bc_dir, bc_loc, k, l)
                case default
                    call s_F_igr_ghost_cell_extrapolation(jac_sf, bc_dir, bc_loc, k, l)
                end select
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_populate_F_igr_bc_direction

    !> Populate the buffers of the grid variables, which are constituted of the cell-boundary locations and cell-width
    !! distributions, based on the boundary conditions.
    subroutine s_populate_grid_variables_buffers(x_cb_in, x_cc_in, dx_in, x_offset, y_offset, z_offset, y_cb_in, y_cc_in, dy_in, &
        & z_cb_in, z_cc_in, dz_in, global_bounds)

        type(int_bounds_info), intent(in)                        :: x_offset, y_offset, z_offset
        real(wp), contiguous, intent(inout)                      :: x_cb_in(-1 - x_offset%beg:)
        real(wp), contiguous, intent(inout)                      :: x_cc_in(-buff_size:), dx_in(-buff_size:)
        real(wp), optional, contiguous, intent(inout)            :: y_cb_in(-1 - y_offset%beg:), z_cb_in(-1 - z_offset%beg:)
        real(wp), optional, contiguous, intent(inout)            :: y_cc_in(-buff_size:), dy_in(-buff_size:)
        real(wp), optional, contiguous, intent(inout)            :: z_cc_in(-buff_size:), dz_in(-buff_size:)
        type(bounds_info), optional, dimension(3), intent(inout) :: global_bounds

        if (present(global_bounds)) then
#ifdef MFC_MPI
            call s_mpi_allreduce_min(x_cb_in(-1), global_bounds(1)%beg)
            call s_mpi_allreduce_max(x_cb_in(m), global_bounds(1)%end)
            if (n > 0) then
                call s_mpi_allreduce_min(y_cb_in(-1), global_bounds(2)%beg)
                call s_mpi_allreduce_max(y_cb_in(n), global_bounds(2)%end)
                if (p > 0) then
                    call s_mpi_allreduce_min(z_cb_in(-1), global_bounds(3)%beg)
                    call s_mpi_allreduce_max(z_cb_in(p), global_bounds(3)%end)
                end if
            end if
#else
            global_bounds(1)%beg = x_cb_in(-1); global_bounds(1)%end = x_cb_in(m)
            if (n > 0) then
                global_bounds(2)%beg = y_cb_in(-1); global_bounds(2)%end = y_cb_in(n)
                if (p > 0) then
                    global_bounds(3)%beg = z_cb_in(-1); global_bounds(3)%end = z_cb_in(p)
                end if
            end if
#endif
        end if

        call s_populate_grid_bc_direction(x_cb_in, x_cc_in, dx_in, m, 1, -1, bc_x, x_offset)
        call s_populate_grid_bc_direction(x_cb_in, x_cc_in, dx_in, m, 1, 1, bc_x, x_offset)

        if (n == 0) return

        #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
            call s_populate_grid_bc_direction(y_cb_in, y_cc_in, dy_in, n, 2, -1, bc_y, y_offset)
            call s_populate_grid_bc_direction(y_cb_in, y_cc_in, dy_in, n, 2, 1, bc_y, y_offset)
        #:endif

        if (p == 0) return

        #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
            call s_populate_grid_bc_direction(z_cb_in, z_cc_in, dz_in, p, 3, -1, bc_z, z_offset)
            call s_populate_grid_bc_direction(z_cb_in, z_cc_in, dz_in, p, 3, 1, bc_z, z_offset)
        #:endif

    end subroutine s_populate_grid_variables_buffers

    !> Populate cell-boundary, cell-center, and cell-width buffers for one coordinate direction.
    subroutine s_populate_grid_bc_direction(cell_boundaries, cell_centers, cell_widths, num_cells, bc_dir, bc_loc, bc_bounds, &
                                            & offset)

        integer, intent(in)                 :: num_cells, bc_dir, bc_loc
        type(int_bounds_info), intent(in)   :: bc_bounds, offset
        real(wp), contiguous, intent(inout) :: cell_boundaries(-1 - offset%beg:)
        real(wp), contiguous, intent(inout) :: cell_centers(-buff_size:), cell_widths(-buff_size:)
        integer                             :: bc_edge, i, source_index

        if (bc_loc == -1) then
            bc_edge = bc_bounds%beg
        else
            bc_edge = bc_bounds%end
        end if

        if (bc_edge >= 0) then
            call s_mpi_sendrecv_grid_variable_buffer(cell_boundaries, cell_centers, cell_widths, num_cells, bc_bounds, bc_loc, &
                & offset)
            return
        end if

        if (bc_edge == BC_AXIS .and. (bc_dir /= 2 .or. bc_loc == 1)) return

        do i = 1, buff_size
            if (bc_loc == -1) then
                select case (bc_edge)
                case (BC_PERIODIC)
                    source_index = num_cells - i + 1
                case (BC_REFLECTIVE, BC_AXIS)
                    source_index = i - 1
                case default
                    source_index = 0
                end select
                cell_widths(-i) = cell_widths(source_index)
            else
                select case (bc_edge)
                case (BC_PERIODIC)
                    source_index = i - 1
                case (BC_REFLECTIVE)
                    source_index = num_cells - i + 1
                case default
                    source_index = num_cells
                end select
                cell_widths(num_cells + i) = cell_widths(source_index)
            end if
        end do

        if (bc_loc == -1) then
            do i = 1, offset%beg
                cell_boundaries(-1 - i) = cell_boundaries(-i) - cell_widths(-i)
            end do
            do i = 1, buff_size
                cell_centers(-i) = cell_centers(1 - i) - (cell_widths(1 - i) + cell_widths(-i))/2._wp
            end do
        else
            do i = 1, offset%end
                cell_boundaries(num_cells + i) = cell_boundaries(num_cells + i - 1) + cell_widths(num_cells + i)
            end do
            do i = 1, buff_size
                cell_centers(num_cells + i) = cell_centers(num_cells + i - 1) + (cell_widths(num_cells + i - 1) &
                             & + cell_widths(num_cells + i))/2._wp
            end do
        end if

    end subroutine s_populate_grid_bc_direction

    !> Deallocate boundary condition buffer arrays allocated during module initialization.
    subroutine s_finalize_boundary_common_module()

        if (bc_io) then
            @:DEALLOCATE(bc_buffers(1, 1)%sf)
            @:DEALLOCATE(bc_buffers(1, 2)%sf)
            #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
                if (n > 0) then
                    @:DEALLOCATE(bc_buffers(2, 1)%sf)
                    @:DEALLOCATE(bc_buffers(2, 2)%sf)
                    #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                        if (p > 0) then
                            @:DEALLOCATE(bc_buffers(3, 1)%sf)
                            @:DEALLOCATE(bc_buffers(3, 2)%sf)
                        end if
                    #:endif
                end if
            #:endif
        end if
        @:DEALLOCATE(bc_buffers)

    end subroutine s_finalize_boundary_common_module

    !> Populate ghost cell buffers of the Lagrangian-bubble beta (void fraction) variables based on the boundary conditions.
    impure subroutine s_populate_beta_buffers(q_beta, kahan_comp, bc_type, nvar)

        type(scalar_field), dimension(:), intent(inout)            :: q_beta
        type(scalar_field), dimension(:), intent(inout)            :: kahan_comp
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        integer, intent(in)                                        :: nvar

        call s_populate_beta_bc_direction(1, -1, bc%x, bc_type(1, 1), q_beta, kahan_comp, nvar)
        call s_populate_beta_bc_direction(1, 1, bc%x, bc_type(1, 2), q_beta, kahan_comp, nvar)

        ! n > 0 always for bubbles_lagrange
        #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
            call s_populate_beta_bc_direction(2, -1, bc%y, bc_type(2, 1), q_beta, kahan_comp, nvar)
            call s_populate_beta_bc_direction(2, 1, bc%y, bc_type(2, 2), q_beta, kahan_comp, nvar)
        #:endif

        if (p == 0) return

        #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
            call s_populate_beta_bc_direction(3, -1, bc%z, bc_type(3, 1), q_beta, kahan_comp, nvar)
            call s_populate_beta_bc_direction(3, 1, bc%z, bc_type(3, 2), q_beta, kahan_comp, nvar)
        #:endif

    end subroutine s_populate_beta_buffers

    !> Populate beta variable buffers for one direction and location, by dispatching the per-cell beta BC routines over the boundary
    !! face and performing the paired MPI reduction for processor boundaries.
    impure subroutine s_populate_beta_bc_direction(bc_dir, bc_loc, bc_bounds, bc_type_edge, q_beta, kahan_comp, nvar)

        integer, intent(in)                             :: bc_dir, bc_loc
        type(int_bounds_info), intent(in)               :: bc_bounds
        type(integer_field), intent(in)                 :: bc_type_edge
        type(scalar_field), dimension(:), intent(inout) :: q_beta
        type(scalar_field), dimension(:), intent(inout) :: kahan_comp
        integer, intent(in)                             :: nvar
        integer                                         :: bc_edge, k_beg, k_end, l_beg, l_end, k, l, bc_code

        if (bc_loc == -1) then
            bc_edge = bc_bounds%beg
        else
            bc_edge = bc_bounds%end
        end if

        if (bc_edge < 0) then
            if (bc_dir == 1) then
                k_beg = beta_bc_bounds(2)%beg; k_end = beta_bc_bounds(2)%end
                l_beg = beta_bc_bounds(3)%beg; l_end = beta_bc_bounds(3)%end
            else if (bc_dir == 2) then
                k_beg = beta_bc_bounds(1)%beg; k_end = beta_bc_bounds(1)%end
                l_beg = beta_bc_bounds(3)%beg; l_end = beta_bc_bounds(3)%end
            else
                k_beg = beta_bc_bounds(1)%beg; k_end = beta_bc_bounds(1)%end
                l_beg = beta_bc_bounds(2)%beg; l_end = beta_bc_bounds(2)%end
            end if

            $:GPU_PARALLEL_LOOP(private='[l, k, bc_code]', collapse=2)
            do l = l_beg, l_end
                do k = k_beg, k_end
                    ! bc_type is not allocated over the beta ghost extents in x and y, so those directions dispatch on the
                    ! domain-edge BC; in z it is allocated with buff_size (>= mapcells + 1) ghost layers and dispatches per cell.
                    if (bc_dir == 3) then
                        bc_code = int(bc_type_edge%sf(k, l, 0))
                    else
                        bc_code = bc_edge
                    end if

                    select case (bc_code)
                    case (BC_PERIODIC)
                        call s_beta_periodic(q_beta, kahan_comp, bc_dir, bc_loc, k, l, nvar)
                    case (BC_REFLECTIVE)
                        call s_beta_reflective(q_beta, kahan_comp, bc_dir, bc_loc, k, l, nvar)
                    end select
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        ! The beta reduction is a paired exchange (rightward accumulate at bc_loc = -1, leftward distribute at bc_loc = 1), so it
        ! must run at both locations whenever either edge of the direction is a processor boundary.
        if (bc_bounds%beg >= 0 .or. bc_bounds%end >= 0) then
            call s_mpi_reduce_beta_variables_buffers(q_beta, kahan_comp, bc_dir, bc_loc, nvar)
        end if

    end subroutine s_populate_beta_bc_direction

end module m_boundary_common
