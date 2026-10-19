!>
!! @file
!! @brief Contains module m_boundary_common

!> @brief Noncharacteristic and processor boundary condition application for ghost cells and buffer regions
#:include 'case.fpp'
#:include 'macros.fpp'

#if defined(MFC_SIMULATION) && defined(MFC_MPI)
#define MFC_IGR_BATCHED_EXCHANGE
#endif

#:def APPLY_VARIABLE_BC()
    select case (bc_code)
    case (BC_CHAR_SUP_OUTFLOW:BC_GHOST_EXTRAP)
        call s_ghost_cell_extrapolation(q_prim_vf, bc_dir, bc_loc, k, l, i, j, q_T_sf)
    case (BC_AXIS)
        if (bc_dir == 2 .and. bc_loc == -1) call s_axis(q_prim_vf, pb_in, mv_in, k, l, i, j, q_T_sf)
    case (BC_REFLECTIVE)
        call s_symmetry(q_prim_vf, bc_dir, bc_loc, k, l, i, j, pb_in, mv_in, q_T_sf)
    case (BC_PERIODIC)
        call s_periodic(q_prim_vf, bc_dir, bc_loc, k, l, i, j, pb_in, mv_in, q_T_sf)
    case (BC_SLIP_WALL)
        call s_slip_wall(q_prim_vf, bc_dir, bc_loc, k, l, i, j, q_T_sf)
    case (BC_NO_SLIP_WALL)
        call s_no_slip_wall(q_prim_vf, bc_dir, bc_loc, k, l, i, j, q_T_sf)
    case (BC_DIRICHLET)
        call s_dirichlet(q_prim_vf, bc_dir, bc_loc, k, l, i, j, q_T_sf)
    end select

    if (qbmm .and. (.not. polytropic) .and. present(pb_in) .and. present(mv_in) .and. (bc_code <= BC_GHOST_EXTRAP) &
        & .and. .not. (bc_dir == 2 .and. bc_loc == -1 .and. bc_code == BC_AXIS)) then
        call s_qbmm_extrapolation(bc_dir, bc_loc, k, l, i, j, pb_in, mv_in)
    end if
#:enddef

#:def PB_SHELL_LOOP()
    $:GPU_PARALLEL_LOOP(private='[slab_cell, face_offset, ix, iy, iz, sx, sy, sz, i, valid_cell, reflect_x, reflect_y, reflect_z, &
                        & value]', firstprivate='[z_beg, x_shell_cells, y_shell_cells]')
    do cell = 0, shell_cells - 1
        if (cell < x_shell_cells) then
            face_offset = mod(cell, 2*buff_size)
            slab_cell = cell/(2*buff_size)
            if (face_offset < buff_size) then
                ix = -face_offset - 1
            else
                ix = m + face_offset - buff_size + 1
            end if
            iy = mod(slab_cell, y_full) - buff_size
            iz = slab_cell/y_full + z_beg
        else if (cell < x_shell_cells + y_shell_cells) then
            slab_cell = cell - x_shell_cells
            ix = mod(slab_cell, m + 1)
            slab_cell = slab_cell/(m + 1)
            face_offset = mod(slab_cell, 2*buff_size)
            if (face_offset < buff_size) then
                iy = -face_offset - 1
            else
                iy = n + face_offset - buff_size + 1
            end if
            iz = slab_cell/(2*buff_size) + z_beg
        else
            slab_cell = cell - x_shell_cells - y_shell_cells
            ix = mod(slab_cell, m + 1)
            slab_cell = slab_cell/(m + 1)
            iy = mod(slab_cell, n + 1)
            face_offset = slab_cell/(n + 1)
            if (face_offset < buff_size) then
                iz = -face_offset - 1
            else
                iz = p + face_offset - buff_size + 1
            end if
        end if

        sx = ix; sy = iy; sz = iz
        valid_cell = .true.
        reflect_x = .false.; reflect_y = .false.; reflect_z = .false.
        if (ix < 0) then
            if (bc_x%beg >= 0) then
                valid_cell = .false.
            else if (bc_x%beg == BC_REFLECTIVE) then
                sx = -ix - 1; reflect_x = .true.
            else
                sx = 0
            end if
        else if (ix > m) then
            if (bc_x%end >= 0) then
                valid_cell = .false.
            else if (bc_x%end == BC_REFLECTIVE) then
                sx = m - (ix - m - 1); reflect_x = .true.
            else
                sx = m
            end if
        end if
        #:if num_dims > 1
            if (iy < 0) then
                if (bc_y%beg >= 0) then
                    valid_cell = .false.
                else if (bc_y%beg == BC_REFLECTIVE) then
                    sy = -iy - 1; reflect_y = .true.
                else
                    sy = 0
                end if
            else if (iy > n) then
                if (bc_y%end >= 0) then
                    valid_cell = .false.
                else if (bc_y%end == BC_REFLECTIVE) then
                    sy = n - (iy - n - 1); reflect_y = .true.
                else
                    sy = n
                end if
            end if
        #:endif
        #:if num_dims > 2
            if (iz < 0) then
                if (bc_z%beg >= 0) then
                    valid_cell = .false.
                else if (bc_z%beg == BC_REFLECTIVE) then
                    sz = -iz - 1; reflect_z = .true.
                else
                    sz = 0
                end if
            else if (iz > p) then
                if (bc_z%end >= 0) then
                    valid_cell = .false.
                else if (bc_z%end == BC_REFLECTIVE) then
                    sz = p - (iz - p - 1); reflect_z = .true.
                else
                    sz = p
                end if
            end if
        #:endif
        if (valid_cell) then
            do i = 1, sys_size
                value = q_prim_vf(i)%sf(sx, sy, sz)
                if (reflect_x .and. i == eqn_idx%mom%beg) value = -value
                #:if num_dims > 1
                    if (reflect_y .and. i == eqn_idx%mom%beg + 1) value = -value
                #:endif
                #:if num_dims > 2
                    if (reflect_z .and. i == eqn_idx%mom%beg + 2) value = -value
                #:endif
                q_prim_vf(i)%sf(ix, iy, iz) = value
            end do
            if ((chemistry .or. heat_conduction) .and. present(q_T_sf)) q_T_sf%sf(ix, iy, iz) = q_T_sf%sf(sx, sy, sz)
        end if
    end do
    $:END_GPU_PARALLEL_LOOP()
#:enddef

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

#ifdef MFC_HALO_SINGLE_ROUND
    !> The single exchange runs only when no rank has a Dirichlet face (see s_check_single_round_bcs).
    logical :: single_round_checked = .false.
    logical :: single_round_ok = .true.
#endif

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
            if (chemistry .or. heat_conduction) sys_size_alloc = sys_size + 1

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
    impure subroutine s_populate_variables_buffers(bc_type, q_prim_vf, pb_in, mv_in, q_T_sf, phase)

        type(scalar_field), dimension(sys_size), intent(inout)                                               :: q_prim_vf
        real(stp), optional, dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_in, mv_in
        type(integer_field), dimension(1:num_dims,1:2), intent(in)                                           :: bc_type
        type(scalar_field), optional, intent(inout)                                                          :: q_T_sf
        integer, optional, intent(in)                                                                        :: phase

#ifdef MFC_HALO_SINGLE_ROUND
        if (.not. single_round_checked) call s_check_single_round_bcs(bc_type)

        ! Phase 2 completes an exchange whose physical boundaries were populated before it was posted.
        if (present(phase)) then
            if (phase == 2) then
                if (single_round_ok) call s_mpi_sendrecv_variables_buffers_26(q_prim_vf, sys_size, pb_in, mv_in, q_T_sf, phase)
                return
            end if
        end if
#endif

#if defined(MFC_PERSISTENT_BC) && defined(MFC_HALO_SINGLE_ROUND)
        if (single_round_ok .and. f_use_persistent_bc(bc_type)) then
            call s_populate_variables_buffers_openacc(bc_type, q_prim_vf, pb_in, mv_in, q_T_sf)
            call s_mpi_sendrecv_variables_buffers_26(q_prim_vf, sys_size, pb_in, mv_in, q_T_sf, phase)
            return
        end if
#endif

        call s_populate_bc_direction(1, bc_x, bc_type(1, 1), bc_type(1, 2), q_prim_vf, pb_in, mv_in, q_T_sf)

        if (n > 0) then
            #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
                call s_populate_bc_direction(2, bc_y, bc_type(2, 1), bc_type(2, 2), q_prim_vf, pb_in, mv_in, q_T_sf)
            #:endif

            if (p > 0) then
                #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                    call s_populate_bc_direction(3, bc_z, bc_type(3, 1), bc_type(3, 2), q_prim_vf, pb_in, mv_in, q_T_sf)
                #:endif
            end if
        end if

#ifdef MFC_HALO_SINGLE_ROUND
        ! Physical boundaries are populated first so their edge and corner values travel with the single exchange.
        if (single_round_ok) call s_mpi_sendrecv_variables_buffers_26(q_prim_vf, sys_size, pb_in, mv_in, q_T_sf, phase)
#endif

    end subroutine s_populate_variables_buffers

#ifdef MFC_HALO_SINGLE_ROUND
    !> Decide once, identically on every rank, whether the single exchange applies. A Dirichlet ghost corner next to a processor
    !! face holds this rank's boundary buffer in the staged exchange but the neighbor's in the single exchange, so any Dirichlet
    !! face selects the staged exchange.
    impure subroutine s_check_single_round_bcs(bc_type)

        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        integer(kind=8)                                            :: dirichlet_loc, dirichlet_glb
        integer                                                    :: i, j

        dirichlet_loc = 0_8
        do i = 1, num_dims
            do j = 1, 2
                dirichlet_loc = dirichlet_loc + int(count(bc_type(i, j)%sf == BC_DIRICHLET), kind=8)
            end do
        end do
        call s_mpi_allreduce_integer_sum(dirichlet_loc, dirichlet_glb)
        single_round_ok = dirichlet_glb == 0_8
        single_round_checked = .true.

    end subroutine s_check_single_round_bcs
#endif

#ifdef MFC_PERSISTENT_BC
    !> Whether a physical boundary code is one the persistent kernel implements exactly. Processor boundaries (non-negative) are
    !! handled by the MPI exchange, not the fill, so they are always admissible here.
    pure logical function f_pb_bc_supported(code)

        integer, intent(in) :: code

        f_pb_bc_supported = code >= 0 .or. code == BC_PERIODIC .or. code == BC_REFLECTIVE .or. (code <= BC_GHOST_EXTRAP &
            & .and. code >= BC_CHAR_SUP_OUTFLOW) .or. code == BC_AXIS .or. code == BC_SLIP_WALL .or. code == BC_NO_SLIP_WALL &
            & .or. code == BC_DIRICHLET

    end function f_pb_bc_supported

    pure logical function f_pb_shell_bc_supported(code)

        integer, intent(in) :: code

        f_pb_shell_bc_supported = code >= 0 .or. code == BC_REFLECTIVE .or. (code >= BC_CHAR_SUP_OUTFLOW &
            & .and. code <= BC_GHOST_EXTRAP)

    end function f_pb_shell_bc_supported

    logical function f_use_persistent_bc(bc_type)

        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type

        f_use_persistent_bc = .false.
        if (num_dims < 2) return
        if (qbmm) return
        if (num_dims < 3 .and. hypoelasticity) return
        if (.not. f_pb_bc_supported(bc_x%beg)) return
        if (.not. f_pb_bc_supported(bc_x%end)) return
        if (.not. f_pb_bc_supported(bc_y%beg)) return
        if (.not. f_pb_bc_supported(bc_y%end)) return
        if (num_dims == 3) then
            if (.not. f_pb_bc_supported(bc_z%beg)) return
            if (.not. f_pb_bc_supported(bc_z%end)) return
        end if
        f_use_persistent_bc = .true.

    end function f_use_persistent_bc
#endif

#ifdef MFC_PERSISTENT_BC
    !> Fill physical boundary faces. Case-optimized uniform reflective/extrapolation boundaries map the complete ghost shell
    !! directly to interior source cells in one kernel; all other boundaries retain direction-ordered kernel barriers and primitive
    !! dispatch.
    impure subroutine s_populate_variables_buffers_openacc(bc_type, q_prim_vf, pb_in, mv_in, q_T_sf)

        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(stp), optional, dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_in, mv_in
        type(scalar_field), optional, intent(inout) :: q_T_sf
        integer :: bc_dir, iloc, i, j, k, l, bc_loc, bc_code
        integer :: loc_beg, loc_end, var_end, k_beg, k_end, l_beg, l_end
        integer :: nloc, extent_k, extent_l, kk, ll, offset
        integer :: gx, gy, gz, sx, sy, sz, mx, my, mz, ex, ey, ez
        integer :: dim_end, nrm, lsrc, bc_iloc, icomp, q
        integer :: m_b, n_b, p_b, buff_size_b, sys_size_b, mom_beg_b, mom_end_b
        integer(kind=8) :: off, work
        real(wp) :: wall_vel, Twall
        real(stp) :: bc_val
        logical :: flip, isothermal

        #:if MFC_CASE_OPTIMIZATION
            integer   :: cell, slab_cell, shell_cells, x_shell_cells, y_shell_cells
            integer   :: y_full, z_full, face_offset
            integer   :: ix, iy, iz, sx, sy, sz, z_beg
            logical   :: simple_uniform, has_physical, valid_cell, reflect_x, reflect_y, reflect_z
            real(stp) :: value
        #:endif

        var_end = sys_size
        if ((chemistry .or. heat_conduction) .and. present(q_T_sf)) var_end = sys_size + 1
        m_b = m; n_b = n; p_b = p; buff_size_b = buff_size; sys_size_b = sys_size
        mom_beg_b = eqn_idx%mom%beg; mom_end_b = eqn_idx%mom%end
        #:if MFC_CASE_OPTIMIZATION
            ! Without boundary patches, each physical bc_type face is uniform and equal to its scalar boundary code.
            simple_uniform = .not. hypoelasticity .and. num_bc_patches == 0 .and. f_pb_shell_bc_supported(bc_x%beg) &
                & .and. f_pb_shell_bc_supported(bc_x%end)
            has_physical = bc_x%beg < 0 .or. bc_x%end < 0
            #:if num_dims > 1
                simple_uniform = simple_uniform .and. f_pb_shell_bc_supported(bc_y%beg) .and. f_pb_shell_bc_supported(bc_y%end)
                has_physical = has_physical .or. bc_y%beg < 0 .or. bc_y%end < 0
            #:endif
            #:if num_dims > 2
                simple_uniform = simple_uniform .and. f_pb_shell_bc_supported(bc_z%beg) .and. f_pb_shell_bc_supported(bc_z%end)
                has_physical = has_physical .or. bc_z%beg < 0 .or. bc_z%end < 0
            #:endif

            if (simple_uniform .and. has_physical) then
                z_beg = 0; z_full = 1
                #:if num_dims > 2
                    z_beg = -buff_size; z_full = p + 1 + 2*buff_size
                #:endif
                y_full = 1
                #:if num_dims > 1
                    y_full = n + 1 + 2*buff_size
                #:endif
                x_shell_cells = 2*buff_size*y_full*z_full
                y_shell_cells = 0
                #:if num_dims > 1
                    y_shell_cells = (m + 1)*2*buff_size*z_full
                #:endif
                shell_cells = x_shell_cells + y_shell_cells
                #:if num_dims > 2
                    shell_cells = shell_cells + (m + 1)*(n + 1)*2*buff_size
                #:endif
                $:PB_SHELL_LOOP()
                return
            end if
        #:endif

        do bc_dir = 1, num_dims
            ! Processor boundaries are exchanged by MPI after this routine returns, so their faces are skipped here.
            loc_beg = 1; loc_end = 2
            if (bc_dir == 1) then
                if (bc_x%beg >= 0) loc_beg = 2
                if (bc_x%end >= 0) loc_end = 1
            else if (bc_dir == 2) then
                if (bc_y%beg >= 0) loc_beg = 2
                if (bc_y%end >= 0) loc_end = 1
            else
                if (bc_z%beg >= 0) loc_beg = 2
                if (bc_z%end >= 0) loc_end = 1
            end if

            ! Per-direction transverse extent: dir 1 stays in the interior, dirs 2 and 3 spread into the ghost region.
            if (bc_dir == 1) then
                k_beg = 0; k_end = n; l_beg = 0; l_end = p
            else if (bc_dir == 2) then
                k_beg = -buff_size; k_end = m + buff_size; l_beg = 0; l_end = p
            else
                k_beg = -buff_size; k_end = m + buff_size; l_beg = -buff_size; l_end = n + buff_size
            end if

            if (loc_beg <= loc_end) then
                nloc = loc_end - loc_beg + 1
                extent_k = k_end - k_beg + 1
                extent_l = l_end - l_beg + 1
                work = int(nloc, 8)*int(var_end, 8)*int(extent_l, 8)*int(extent_k, 8)
                $:GPU_PARALLEL_LOOP(private='[offset, kk, ll, i, iloc, k, l, bc_loc, bc_code, j, q]', firstprivate='[extent_k, &
                                    & extent_l, var_end, loc_beg, k_beg, l_beg, m_b, n_b, p_b, buff_size_b, sys_size_b, &
                                    & mom_beg_b, mom_end_b]')
                do off = 0_8, work - 1_8
                    offset = int(off)
                    kk = mod(offset, extent_k)
                    offset = offset/extent_k
                    ll = mod(offset, extent_l)
                    offset = offset/extent_l
                    i = mod(offset, var_end) + 1
                    offset = offset/var_end
                    iloc = loc_beg + offset

                    k = k_beg + kk
                    l = l_beg + ll
                    bc_loc = 2*iloc - 3

                    if (bc_dir == 1) then
                        if (iloc == 1) then
                            bc_code = int(bc_type(1, 1)%sf(0, k, l))
                        else
                            bc_code = int(bc_type(1, 2)%sf(0, k, l))
                        end if
                    else if (bc_dir == 2) then
                        if (iloc == 1) then
                            bc_code = int(bc_type(2, 1)%sf(k, 0, l))
                        else
                            bc_code = int(bc_type(2, 2)%sf(k, 0, l))
                        end if
                    else
                        if (iloc == 1) then
                            bc_code = int(bc_type(3, 1)%sf(k, l, 0))
                        else
                            bc_code = int(bc_type(3, 2)%sf(k, l, 0))
                        end if
                    end if

                    if (bc_dir == 1) then
                        dim_end = m_b
                    else if (bc_dir == 2) then
                        dim_end = n_b
                    else
                        dim_end = p_b
                    end if

                    do j = 1, buff_size_b
                        if (bc_loc == -1) then
                            nrm = -j
                        else
                            nrm = dim_end + j
                        end if
                        if (bc_dir == 1) then
                            gx = nrm; gy = k; gz = l
                        else if (bc_dir == 2) then
                            gx = k; gy = nrm; gz = l
                        else
                            gx = k; gy = l; gz = nrm
                        end if

                        select case (bc_code)
                        case (BC_CHAR_SUP_OUTFLOW:BC_GHOST_EXTRAP)
                            nrm = merge(0, dim_end, bc_loc == -1)
                            if (bc_dir == 1) then
                                sx = nrm; sy = k; sz = l
                            else if (bc_dir == 2) then
                                sx = k; sy = nrm; sz = l
                            else
                                sx = k; sy = l; sz = nrm
                            end if
                            if (i <= sys_size_b) then
                                q_prim_vf(i)%sf(gx, gy, gz) = q_prim_vf(i)%sf(sx, sy, sz)
                            else if (i == sys_size_b + 1 .and. present(q_T_sf)) then
                                q_T_sf%sf(gx, gy, gz) = q_T_sf%sf(sx, sy, sz)
                            end if
                        case (BC_AXIS)
                            if (bc_dir == 2 .and. bc_loc == -1) then
                                if (z_cc(l) < pi) then
                                    lsrc = l + (p_b + 1)/2
                                else
                                    lsrc = l - (p_b + 1)/2
                                end if
                                if (i <= sys_size_b) then
                                    if (i == mom_beg_b + 1 .or. i == mom_end_b) then
                                        q_prim_vf(i)%sf(k, -j, l) = -q_prim_vf(i)%sf(k, j - 1, lsrc)
                                    else
                                        q_prim_vf(i)%sf(k, -j, l) = q_prim_vf(i)%sf(k, j - 1, lsrc)
                                    end if
                                else if (i == sys_size_b + 1 .and. present(q_T_sf)) then
                                    q_T_sf%sf(k, -j, l) = q_T_sf%sf(k, j - 1, lsrc)
                                end if
                            end if
                        case (BC_REFLECTIVE)
                            nrm = merge(j - 1, dim_end - (j - 1), bc_loc == -1)
                            if (bc_dir == 1) then
                                sx = nrm; sy = k; sz = l
                            else if (bc_dir == 2) then
                                sx = k; sy = nrm; sz = l
                            else
                                sx = k; sy = l; sz = nrm
                            end if
                            if (i <= sys_size_b) then
                                flip = i == mom_beg_b + bc_dir - 1
                                if (hypoelasticity) then
                                    do q = 1, shear_BC_flip_num
                                        if (i == shear_BC_flip_indices(bc_dir, q)) flip = .true.
                                    end do
                                end if
                                if (flip) then
                                    q_prim_vf(i)%sf(gx, gy, gz) = -q_prim_vf(i)%sf(sx, sy, sz)
                                else
                                    q_prim_vf(i)%sf(gx, gy, gz) = q_prim_vf(i)%sf(sx, sy, sz)
                                end if
                            end if
                            if (i == sys_size_b + 1 .and. present(q_T_sf)) q_T_sf%sf(gx, gy, gz) = q_T_sf%sf(sx, sy, sz)
                        case (BC_PERIODIC)
                            nrm = merge(dim_end - (j - 1), j - 1, bc_loc == -1)
                            if (bc_dir == 1) then
                                sx = nrm; sy = k; sz = l
                            else if (bc_dir == 2) then
                                sx = k; sy = nrm; sz = l
                            else
                                sx = k; sy = l; sz = nrm
                            end if
                            if (i <= sys_size_b) then
                                q_prim_vf(i)%sf(gx, gy, gz) = q_prim_vf(i)%sf(sx, sy, sz)
                            end if
                            if (i == sys_size_b + 1 .and. present(q_T_sf)) q_T_sf%sf(gx, gy, gz) = q_T_sf%sf(sx, sy, sz)
                        case (BC_SLIP_WALL, BC_NO_SLIP_WALL)
                            if (bc_loc == -1) then
                                nrm = j - 1
                            else
                                nrm = dim_end - (j - 1)
                            end if
                            if (bc_dir == 1) then
                                mx = nrm; my = k; mz = l
                                ex = merge(0, dim_end, bc_loc == -1); ey = k; ez = l
                            else if (bc_dir == 2) then
                                mx = k; my = nrm; mz = l
                                ex = k; ey = merge(0, dim_end, bc_loc == -1); ez = l
                            else
                                mx = k; my = l; mz = nrm
                                ex = k; ey = l; ez = merge(0, dim_end, bc_loc == -1)
                            end if

                            wall_vel = 0._wp
                            isothermal = .false.
                            Twall = 0._wp
                            if (bc_dir == 1) then
                                if (bc_loc == -1) then
                                    isothermal = bc_x%isothermal_in; Twall = bc_x%Twall_in
                                    if (bc_code == BC_SLIP_WALL) then
                                        wall_vel = bc_x%vb1
                                    else
                                        icomp = i - mom_beg_b + 1
                                        if (icomp == 1) wall_vel = bc_x%vb1
                                        if (icomp == 2) wall_vel = bc_x%vb2
                                        if (icomp == 3) wall_vel = bc_x%vb3
                                    end if
                                else
                                    isothermal = bc_x%isothermal_out; Twall = bc_x%Twall_out
                                    if (bc_code == BC_SLIP_WALL) then
                                        wall_vel = bc_x%ve1
                                    else
                                        icomp = i - mom_beg_b + 1
                                        if (icomp == 1) wall_vel = bc_x%ve1
                                        if (icomp == 2) wall_vel = bc_x%ve2
                                        if (icomp == 3) wall_vel = bc_x%ve3
                                    end if
                                end if
                            else if (bc_dir == 2) then
                                if (bc_loc == -1) then
                                    isothermal = bc_y%isothermal_in; Twall = bc_y%Twall_in
                                    if (bc_code == BC_SLIP_WALL) then
                                        wall_vel = bc_y%vb2
                                    else
                                        icomp = i - mom_beg_b + 1
                                        if (icomp == 1) wall_vel = bc_y%vb1
                                        if (icomp == 2) wall_vel = bc_y%vb2
                                        if (icomp == 3) wall_vel = bc_y%vb3
                                    end if
                                else
                                    isothermal = bc_y%isothermal_out; Twall = bc_y%Twall_out
                                    if (bc_code == BC_SLIP_WALL) then
                                        wall_vel = bc_y%ve2
                                    else
                                        icomp = i - mom_beg_b + 1
                                        if (icomp == 1) wall_vel = bc_y%ve1
                                        if (icomp == 2) wall_vel = bc_y%ve2
                                        if (icomp == 3) wall_vel = bc_y%ve3
                                    end if
                                end if
                            else
                                if (bc_loc == -1) then
                                    isothermal = bc_z%isothermal_in; Twall = bc_z%Twall_in
                                    if (bc_code == BC_SLIP_WALL) then
                                        wall_vel = bc_z%vb3
                                    else
                                        icomp = i - mom_beg_b + 1
                                        if (icomp == 1) wall_vel = bc_z%vb1
                                        if (icomp == 2) wall_vel = bc_z%vb2
                                        if (icomp == 3) wall_vel = bc_z%vb3
                                    end if
                                else
                                    isothermal = bc_z%isothermal_out; Twall = bc_z%Twall_out
                                    if (bc_code == BC_SLIP_WALL) then
                                        wall_vel = bc_z%ve3
                                    else
                                        icomp = i - mom_beg_b + 1
                                        if (icomp == 1) wall_vel = bc_z%ve1
                                        if (icomp == 2) wall_vel = bc_z%ve2
                                        if (icomp == 3) wall_vel = bc_z%ve3
                                    end if
                                end if
                            end if

                            if (i <= sys_size_b) then
                                icomp = i - mom_beg_b + 1
                                flip = icomp == bc_dir
                                if (bc_code == BC_NO_SLIP_WALL) flip = icomp >= 1 .and. icomp <= num_dims
                                if (flip) then
                                    q_prim_vf(i)%sf(gx, gy, gz) = -q_prim_vf(i)%sf(mx, my, mz) + 2._wp*wall_vel
                                else
                                    q_prim_vf(i)%sf(gx, gy, gz) = q_prim_vf(i)%sf(ex, ey, ez)
                                end if
                            end if
                            if (i == sys_size_b + 1 .and. present(q_T_sf)) then
                                if (isothermal) then
                                    q_T_sf%sf(gx, gy, gz) = 2._wp*Twall - q_T_sf%sf(mx, my, mz)
                                else
                                    q_T_sf%sf(gx, gy, gz) = q_T_sf%sf(ex, ey, ez)
                                end if
                            end if
                        case (BC_DIRICHLET)
                            bc_iloc = (bc_loc + 3)/2
                            if (bc_dir == 1) then
                                bc_val = bc_buffers(1, bc_iloc)%sf(i, k, l)
                            else if (bc_dir == 2) then
                                bc_val = bc_buffers(2, bc_iloc)%sf(k, i, l)
                            else
                                bc_val = bc_buffers(3, bc_iloc)%sf(k, l, i)
                            end if
                            if (i <= sys_size_b) then
                                q_prim_vf(i)%sf(gx, gy, gz) = bc_val*f_dir_vel(i, bc_dir)
                            else if (i == sys_size_b + 1 .and. present(q_T_sf)) then
                                q_T_sf%sf(gx, gy, gz) = bc_val
                            end if
                        end select
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        end do

    end subroutine s_populate_variables_buffers_openacc
#endif

    !> Populate the variable buffers along one direction, both locations in one device loop, via MPI exchange for processor
    !! boundaries and the per-cell BC routines over the physical boundary faces.
    impure subroutine s_populate_bc_direction(bc_dir, bc_bounds, bc_type_beg, bc_type_end, q_prim_vf, pb_in, mv_in, q_T_sf)

        integer, intent(in) :: bc_dir
        type(int_bounds_info), intent(in) :: bc_bounds
        type(integer_field), intent(in) :: bc_type_beg, bc_type_end
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(stp), optional, dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_in, mv_in
        type(scalar_field), optional, intent(inout) :: q_T_sf
        integer :: k_beg, k_end, l_beg, l_end, loc_beg, loc_end, var_end
        integer :: bc_code, k, l, i, j, iloc, bc_loc

        loc_beg = 1
        loc_end = 2

        if (bc_bounds%beg >= 0) then
            loc_beg = 2
#ifdef MFC_HALO_SINGLE_ROUND
            if (.not. single_round_ok) call s_mpi_sendrecv_variables_buffers(q_prim_vf, bc_dir, -1, sys_size, pb_in, mv_in, q_T_sf)
#else
            call s_mpi_sendrecv_variables_buffers(q_prim_vf, bc_dir, -1, sys_size, pb_in, mv_in, q_T_sf)
#endif
        end if

        if (bc_bounds%end >= 0) then
            loc_end = 1
#ifdef MFC_HALO_SINGLE_ROUND
            if (.not. single_round_ok) call s_mpi_sendrecv_variables_buffers(q_prim_vf, bc_dir, 1, sys_size, pb_in, mv_in, q_T_sf)
#else
            call s_mpi_sendrecv_variables_buffers(q_prim_vf, bc_dir, 1, sys_size, pb_in, mv_in, q_T_sf)
#endif
        end if

        if (loc_beg > loc_end) return

        if (bc_dir == 1) then
            k_beg = 0; k_end = n; l_beg = 0; l_end = p
        else if (bc_dir == 2) then
            k_beg = -buff_size; k_end = m + buff_size; l_beg = 0; l_end = p
        else
            k_beg = -buff_size; k_end = m + buff_size; l_beg = -buff_size; l_end = n + buff_size
        end if
        var_end = sys_size
        if ((chemistry .or. heat_conduction) .and. present(q_T_sf)) var_end = sys_size + 1
        if (qbmm .and. (.not. polytropic) .and. present(pb_in) .and. present(mv_in)) var_end = max(var_end, nb*nnode)
        if (n > 0) then
            $:GPU_PARALLEL_LOOP(private='[l, k, i, j, bc_code, iloc, bc_loc]', collapse=4)
            do iloc = loc_beg, loc_end
                do i = 1, var_end
                    do l = l_beg, l_end
                        do k = k_beg, k_end
                            bc_loc = 2*iloc - 3
                            if (bc_dir == 1) then
                                if (iloc == 1) then
                                    bc_code = int(bc_type_beg%sf(0, k, l))
                                else
                                    bc_code = int(bc_type_end%sf(0, k, l))
                                end if
                            else if (bc_dir == 2) then
                                if (iloc == 1) then
                                    bc_code = int(bc_type_beg%sf(k, 0, l))
                                else
                                    bc_code = int(bc_type_end%sf(k, 0, l))
                                end if
                            else
                                if (iloc == 1) then
                                    bc_code = int(bc_type_beg%sf(k, l, 0))
                                else
                                    bc_code = int(bc_type_end%sf(k, l, 0))
                                end if
                            end if

                            do j = 1, buff_size
                                @:APPLY_VARIABLE_BC()
                            end do
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            return
        end if

        $:GPU_PARALLEL_LOOP(private='[l, k, i, j, bc_code, iloc, bc_loc]', collapse=5)
        do iloc = loc_beg, loc_end
            do i = 1, var_end
                do j = 1, buff_size
                    do l = l_beg, l_end
                        do k = k_beg, k_end
                            bc_loc = 2*iloc - 3
                            if (bc_dir == 1) then
                                if (iloc == 1) then
                                    bc_code = int(bc_type_beg%sf(0, k, l))
                                else
                                    bc_code = int(bc_type_end%sf(0, k, l))
                                end if
                            else if (bc_dir == 2) then
                                if (iloc == 1) then
                                    bc_code = int(bc_type_beg%sf(k, 0, l))
                                else
                                    bc_code = int(bc_type_end%sf(k, 0, l))
                                end if
                            else
                                if (iloc == 1) then
                                    bc_code = int(bc_type_beg%sf(k, l, 0))
                                else
                                    bc_code = int(bc_type_end%sf(k, l, 0))
                                end if
                            end if

                            @:APPLY_VARIABLE_BC()
                        end do
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_populate_bc_direction

    !> Populate ghost cell buffers for the color function and its divergence used in capillary surface tension.
    impure subroutine s_populate_capillary_buffers(c_divs, bc_type, bc)

        type(scalar_field), dimension(num_dims + 1), intent(inout) :: c_divs
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        type(bc_xyz_info), intent(in)                              :: bc

        call s_populate_capillary_bc_direction(1, bc%x, bc_type(1, 1), bc_type(1, 2), c_divs)

        if (n == 0) return

        #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
            call s_populate_capillary_bc_direction(2, bc%y, bc_type(2, 1), bc_type(2, 2), c_divs)
        #:endif

        if (p == 0) return

        #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
            call s_populate_capillary_bc_direction(3, bc%z, bc_type(3, 1), bc_type(3, 2), c_divs)
        #:endif

    end subroutine s_populate_capillary_buffers

    !> Populate ghost cell buffers for one capillary BC direction, both locations at once, via MPI exchange for processor boundaries
    !! or by dispatching the per-cell capillary BC routines over both boundary faces in one device loop.
    impure subroutine s_populate_capillary_bc_direction(bc_dir, bc_bounds, bc_type_beg, bc_type_end, c_divs)

        integer, intent(in)                                        :: bc_dir
        type(int_bounds_info), intent(in)                          :: bc_bounds
        type(scalar_field), dimension(num_dims + 1), intent(inout) :: c_divs
        type(integer_field), intent(in)                            :: bc_type_beg, bc_type_end
        integer                                                    :: k_beg, k_end, l_beg, l_end, loc_beg, loc_end
        integer                                                    :: k, l, iloc, bc_loc, bc_code

        loc_beg = 1
        loc_end = 2

        if (bc_bounds%beg >= 0) then
            loc_beg = 2
            call s_mpi_sendrecv_variables_buffers(c_divs, bc_dir, -1, num_dims + 1)
        end if

        if (bc_bounds%end >= 0) then
            loc_end = 1
            call s_mpi_sendrecv_variables_buffers(c_divs, bc_dir, 1, num_dims + 1)
        end if

        if (loc_beg > loc_end) return

        if (bc_dir == 1) then
            k_beg = 0; k_end = n; l_beg = 0; l_end = p
        else if (bc_dir == 2) then
            k_beg = -buff_size; k_end = m + buff_size; l_beg = 0; l_end = p
        else
            k_beg = -buff_size; k_end = m + buff_size; l_beg = -buff_size; l_end = n + buff_size
        end if

        $:GPU_PARALLEL_LOOP(private='[l, k, bc_code, iloc, bc_loc]', collapse=3)
        do iloc = loc_beg, loc_end
            do l = l_beg, l_end
                do k = k_beg, k_end
                    bc_loc = 2*iloc - 3
                    if (bc_dir == 1) then
                        if (iloc == 1) then
                            bc_code = int(bc_type_beg%sf(0, k, l))
                        else
                            bc_code = int(bc_type_end%sf(0, k, l))
                        end if
                    else if (bc_dir == 2) then
                        if (iloc == 1) then
                            bc_code = int(bc_type_beg%sf(k, 0, l))
                        else
                            bc_code = int(bc_type_end%sf(k, 0, l))
                        end if
                    else
                        if (iloc == 1) then
                            bc_code = int(bc_type_beg%sf(k, l, 0))
                        else
                            bc_code = int(bc_type_end%sf(k, l, 0))
                        end if
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
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_populate_capillary_bc_direction

    !> Populate ghost cell buffers for the Jacobian scalar field used in the IGR elliptic solver.
    impure subroutine s_populate_F_igr_buffers(bc_type, jac_sf, halo_depth)

        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        type(scalar_field), dimension(1:), intent(inout)           :: jac_sf
        integer, optional, intent(in)                              :: halo_depth

#ifdef MFC_IGR_BATCHED_EXCHANGE
        call s_populate_F_igr_buffers_batched(bc_type, jac_sf, halo_depth)
        return
#endif

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

#ifdef MFC_IGR_BATCHED_EXCHANGE
    !> Populate the IGR Jacobian ghost cells with all processor faces packed in the boundary kernel and exchanged in one batch.
    impure subroutine s_populate_F_igr_buffers_batched(bc_type, jac_sf, halo_depth)

        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        type(scalar_field), dimension(1:), intent(inout)           :: jac_sf
        integer, optional, intent(in)                              :: halo_depth
        type(int_bounds_info), dimension(1:3)                      :: bc_bounds
        integer                                                    :: dir, iloc, face, bc_edge, depth
        integer                                                    :: fill_faces, pack_faces
        integer                                                    :: k, l, k_max, l_max, bc_dir, bc_loc, bc_code
        logical                                                    :: on_face

        depth = buff_size
        if (present(halo_depth)) depth = max(1, min(halo_depth, buff_size))

        bc_bounds(1) = bc_x; bc_bounds(2) = bc_y; bc_bounds(3) = bc_z
        fill_faces = 0
        pack_faces = 0
        do dir = 1, num_dims
            do iloc = 1, 2
                bc_edge = bc_bounds(dir)%beg
                if (iloc == 2) bc_edge = bc_bounds(dir)%end

                if (bc_edge < 0) then
                    fill_faces = ibset(fill_faces, 2*(dir - 1) + iloc - 1)
                else
                    pack_faces = ibset(pack_faces, 2*(dir - 1) + iloc - 1)
                end if
            end do
        end do
        k_max = n
        l_max = p
        if (num_dims > 1) k_max = max(m, n)
        if (num_dims > 2) l_max = max(n, p)

        $:GPU_PARALLEL_LOOP(private='[face, l, k, iloc, bc_dir, bc_loc, bc_code, on_face]', collapse=3)
        do face = 0, 2*num_dims - 1
            do l = 0, l_max
                do k = 0, k_max
                    iloc = mod(face, 2) + 1
                    bc_dir = face/2 + 1
                    bc_loc = 2*iloc - 3
                    bc_code = 0

                    if (bc_dir == 1) then
                        on_face = (k <= n .and. l <= p)
                    else if (bc_dir == 2) then
                        on_face = (k <= m .and. l <= p)
                    else
                        on_face = (k <= m .and. l <= n)
                    end if

                    if (on_face .and. btest(fill_faces, face)) then
                        if (bc_dir == 1) then
                            bc_code = int(bc_type(bc_dir, iloc)%sf(0, k, l))
                        else if (bc_dir == 2) then
                            bc_code = int(bc_type(bc_dir, iloc)%sf(k, 0, l))
                        else
                            bc_code = int(bc_type(bc_dir, iloc)%sf(k, l, 0))
                        end if

                        select case (bc_code)
                        case (BC_PERIODIC)
                            call s_F_igr_periodic(jac_sf, bc_dir, bc_loc, k, l)
                        case (BC_REFLECTIVE)
                            call s_F_igr_reflective(jac_sf, bc_dir, bc_loc, k, l)
                        case default
                            call s_F_igr_ghost_cell_extrapolation(jac_sf, bc_dir, bc_loc, k, l)
                        end select
                    else if (on_face .and. btest(pack_faces, face)) then
                        call s_pack_F_igr_face(jac_sf, face, k, l, depth)
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        call s_mpi_sendrecv_F_igr_batched(jac_sf, depth, pack_faces)

    end subroutine s_populate_F_igr_buffers_batched
#endif

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
