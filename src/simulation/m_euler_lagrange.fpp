!>
!! @file
!! @brief Contains module m_euler_lagrange

#:include 'macros.fpp'

!> @brief Utilities shared by the Euler-Lagrange bubble and solid-particle solvers: cell location, domain and cell-volume helpers,
!! boundary handling, RK stages and Fornberg gradients. Their file I/O is in m_data_output (writes) and m_start_up (reads).
module m_euler_lagrange

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_mpi_common
    use m_compile_specific

    implicit none

    private
    public :: s_locate_cell, s_get_cell, f_lag_entity_in_physical_domain, s_get_char_vol, s_get_lag_restart_point, &
        & s_set_lag_comm_coords, s_count_lag_entities_glb, s_apply_lag_entity_domain_bc, s_wrap_lag_entities_periodic, &
        & f_lag_rk_stage, f_lag_entity_in_wall, s_initialize_lag_gradient, s_finalize_lag_gradient, s_lag_gradient_dir

    !> Fornberg first-derivative weights at each cell center along x, y and z, one field per stencil point
    type(scalar_field), dimension(:), allocatable :: weights_x_grad, weights_y_grad, weights_z_grad
    integer                                       :: nWeights_grad  !< Number of stencil points
    $:GPU_DECLARE(create='[weights_x_grad, weights_y_grad, weights_z_grad, nWeights_grad]')
    !> Lagrangian volume fraction above which a cell counts in the cloud statistics
    real(wp), parameter :: void_stat_min = 5.e-11_wp

contains

    !> Save index and time the run starts from: zero on a fresh start, the restart point otherwise.
    subroutine s_get_lag_restart_point(save_count, qtime)

        integer, intent(out)  :: save_count
        real(wp), intent(out) :: qtime

        if (cfl_dt) then
            save_count = n_start
            qtime = n_start*t_save
        else
            save_count = t_step_start
            qtime = t_step_start*dt
        end if

    end subroutine s_get_lag_restart_point

    !> Physical coordinates for triggering MPI communication of a Lagrangian entity.
    impure subroutine s_set_lag_comm_coords()

        pcomm_coords(1)%beg = x_cb(-1)
        pcomm_coords(1)%end = x_cb(m)
        $:GPU_UPDATE(device='[pcomm_coords(1)]')
        if (n > 0) then
            pcomm_coords(2)%beg = y_cb(-1)
            pcomm_coords(2)%end = y_cb(n)
            $:GPU_UPDATE(device='[pcomm_coords(2)]')
            if (num_dims == 3) then
                pcomm_coords(3)%beg = z_cb(-1)
                pcomm_coords(3)%end = z_cb(p)
                $:GPU_UPDATE(device='[pcomm_coords(3)]')
            end if
        end if

    end subroutine s_set_lag_comm_coords

    !> Sum the local count n_loc over the ranks into n_glb and abort if no entity lies in the domain.
    impure subroutine s_count_lag_entities_glb(n_loc, n_glb, path)

        integer, intent(in)          :: n_loc
        integer, intent(out)         :: n_glb
        character(len=*), intent(in) :: path

        if (num_procs > 1) then
            call s_mpi_reduce_int_sum(n_loc, n_glb)
        else
            n_glb = n_loc
        end if
        if (proc_rank == 0 .and. n_glb == 0) call s_mpi_abort('No Lagrangian entities in the domain. Check ' // trim(path))

    end subroutine s_count_lag_entities_glb

    !> Locate the cell holding pos by bisection on the cell boundaries, and its computational coordinates scoord (cell index plus
    !! the fraction of the cell width). Positions beyond the buffer region clamp to its edge cell.
    subroutine s_locate_cell(pos, cell, scoord)

        $:GPU_ROUTINE(function_name='s_locate_cell',parallelism='[seq]', cray_inline=True)

        real(wp), dimension(3), intent(in)  :: pos
        real(wp), dimension(3), intent(out) :: scoord
        integer, dimension(3), intent(out)  :: cell
        integer                             :: lo, hi, mid

        cell(3) = 0
        #:for D, CB, N in [(1, 'x_cb', 'm'), (2, 'y_cb', 'n'), (3, 'z_cb', 'p')]
            if (${D}$ < 3 .or. num_dims == 3) then
                ! Smallest cell in the buffered range whose right boundary exceeds pos
                lo = -buff_size; hi = ${N}$ + buff_size
                do while (lo < hi)
                    mid = lo + (hi - lo)/2
                    if (pos(${D}$) >= ${CB}$(mid)) then
                        lo = mid + 1
                    else
                        hi = mid
                    end if
                end do
                cell(${D}$) = lo
            end if
        #:endfor

        ! Cell 0 starts at the domain boundary, so the center of cell i is at scoord = i + 1/2
        scoord(1) = cell(1) + (pos(1) - x_cb(cell(1) - 1))/dx(cell(1))
        scoord(2) = cell(2) + (pos(2) - y_cb(cell(2) - 1))/dy(cell(2))
        scoord(3) = 0._wp
        if (num_dims == 3) scoord(3) = cell(3) + (pos(3) - z_cb(cell(3) - 1))/dz(cell(3))

    end subroutine s_locate_cell

    !> Convert computational (real) cell coordinates to integer cell indices
    subroutine s_get_cell(s_cell, get_cell)

        $:GPU_ROUTINE(function_name='s_get_cell',parallelism='[seq]', cray_inline=True)

        real(wp), dimension(3), intent(in) :: s_cell
        integer, dimension(3), intent(out) :: get_cell
        integer                            :: i

        get_cell(:) = int(s_cell(:))
        do i = 1, num_dims
            if (s_cell(i) < 0._wp) get_cell(i) = get_cell(i) - 1
        end do

    end subroutine s_get_cell

    !> True if pos lies in this rank's physical domain (ghost cells excluded).
    logical function f_lag_entity_in_physical_domain(pos_part) result(in_domain)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), dimension(3), intent(in) :: pos_part

        in_domain = pos_part(1) < x_cb(m) .and. pos_part(1) >= x_cb(-1) .and. pos_part(2) < y_cb(n) .and. pos_part(2) >= y_cb(-1)

        if (num_dims == 3) then
            in_domain = in_domain .and. pos_part(3) < z_cb(p) .and. pos_part(3) >= z_cb(-1)
        end if

    end function f_lag_entity_in_physical_domain

    !> True if cell is a ghost cell beyond one of this rank's wall or symmetry boundaries.
    logical function f_lag_entity_in_wall(cell) result(in_wall)

        integer, dimension(3), intent(in) :: cell

        in_wall = (f_is_lag_entity_wall(bc_x%beg) .and. cell(1) < 0) .or. (f_is_lag_entity_wall(bc_x%end) .and. cell(1) > m) &
                   & .or. (f_is_lag_entity_wall(bc_y%beg) .and. cell(2) < 0) .or. (f_is_lag_entity_wall(bc_y%end) .and. cell(2) > n)
        if (num_dims == 3) in_wall = in_wall .or. (f_is_lag_entity_wall(bc_z%beg) .and. cell(3) < 0) &
            & .or. (f_is_lag_entity_wall(bc_z%end) .and. cell(3) > p)

    end function f_lag_entity_in_wall

    !> Allocate the Fornberg first-derivative weights for an n_weights-point stencil and compute them on this rank's grid.
    impure subroutine s_initialize_lag_gradient(n_weights)

        integer, intent(in) :: n_weights
        integer             :: i

        nWeights_grad = n_weights
        $:GPU_UPDATE(device='[nWeights_grad]')
        @:ALLOCATE(weights_x_grad(1:nWeights_grad), weights_y_grad(1:nWeights_grad), weights_z_grad(1:nWeights_grad))
        do i = 1, nWeights_grad
            @:ALLOCATE(weights_x_grad(i)%sf(idwbuff(1)%beg:idwbuff(1)%end,1:1,1:1))
            @:ALLOCATE(weights_y_grad(i)%sf(idwbuff(2)%beg:idwbuff(2)%end,1:1,1:1))
            @:ALLOCATE(weights_z_grad(i)%sf(idwbuff(3)%beg:idwbuff(3)%end,1:1,1:1))
            @:ACC_SETUP_SFs(weights_x_grad(i), weights_y_grad(i), weights_z_grad(i))
        end do

        call s_compute_fornberg_fd_weights((nWeights_grad - 1)/2)

    end subroutine s_initialize_lag_gradient

    !> Deallocate the Fornberg first-derivative weights.
    impure subroutine s_finalize_lag_gradient()

        integer :: i

        do i = 1, nWeights_grad
            @:DEALLOCATE(weights_x_grad(i)%sf, weights_y_grad(i)%sf, weights_z_grad(i)%sf)
        end do
        @:DEALLOCATE(weights_x_grad, weights_y_grad, weights_z_grad)

    end subroutine s_finalize_lag_gradient

    !> Derivative of a cell-centered field along dir with the precomputed Fornberg finite-difference weights.
    subroutine s_lag_gradient_dir(q, dq, dir)

        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:), intent(in)  :: q
        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:), intent(out) :: dq
        integer, intent(in)                                                                :: dir
        integer                                                                            :: i, j, k, a, npts, s_idx
        integer, dimension(3)                                                              :: lo, hi

        npts = (nWeights_grad - 1)/2
        ! The stencil needs npts cells on each side along dir only
        lo = [idwbuff(1)%beg, idwbuff(2)%beg, idwbuff(3)%beg]
        hi = [idwbuff(1)%end, idwbuff(2)%end, idwbuff(3)%end]
        lo(dir) = lo(dir) + npts
        hi(dir) = hi(dir) - npts

        #:for DIR, X, IDX, SHIFT in [(1, 'x', 'i', 'i + a, j, k'), (2, 'y', 'j', 'i, j + a, k'), (3, 'z', 'k', 'i, j, k + a')]
            if (dir == ${DIR}$) then
                $:GPU_PARALLEL_LOOP(private='[i, j, k, s_idx, a]', collapse=3, copyin='[npts, lo, hi]')
                do k = lo(3), hi(3)
                    do j = lo(2), hi(2)
                        do i = lo(1), hi(1)
                            dq(i, j, k) = 0._wp
                            do a = -npts, npts
                                s_idx = a + npts + 1
                                dq(i, j, k) = dq(i, j, k) + weights_${X}$_grad(s_idx)%sf(${IDX}$, 1, 1)*q(${SHIFT}$)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

    end subroutine s_lag_gradient_dir

    !> Fornberg first-derivative weights at every cell center in each direction, computed once on the host and copied to the device.
    impure subroutine s_compute_fornberg_fd_weights(npts)

        integer, intent(in) :: npts
        integer             :: i, j, k, a, m_order
        integer             :: s_idx
        real(wp)            :: x0, y0, z0
        real(wp)            :: x_stencil(nWeights_grad)
        real(wp)            :: c(nWeights_grad,0:1)

        m_order = 1  ! first derivative

        ! Computed once on the host: the run-time-sized stencil arrays must not be GPU-private
        do i = idwbuff(1)%beg + npts, idwbuff(1)%end - npts
            do a = -npts, npts
                s_idx = a + npts + 1
                x_stencil(s_idx) = x_cc(i + a)
            end do
            x0 = x_cc(i)

            call s_fornberg_weights(x0, x_stencil, nWeights_grad, m_order, c)

            do a = -npts, npts
                s_idx = a + npts + 1
                weights_x_grad(s_idx)%sf(i, 1, 1) = c(s_idx, 1)
            end do
        end do

        do j = idwbuff(2)%beg + npts, idwbuff(2)%end - npts
            do a = -npts, npts
                s_idx = a + npts + 1
                x_stencil(s_idx) = y_cc(j + a)
            end do
            y0 = y_cc(j)

            call s_fornberg_weights(y0, x_stencil, nWeights_grad, m_order, c)

            do a = -npts, npts
                s_idx = a + npts + 1
                weights_y_grad(s_idx)%sf(j, 1, 1) = c(s_idx, 1)
            end do
        end do

        if (num_dims == 3) then
            do k = idwbuff(3)%beg + npts, idwbuff(3)%end - npts
                do a = -npts, npts
                    s_idx = a + npts + 1
                    x_stencil(s_idx) = z_cc(k + a)
                end do
                z0 = z_cc(k)

                call s_fornberg_weights(z0, x_stencil, nWeights_grad, m_order, c)

                do a = -npts, npts
                    s_idx = a + npts + 1
                    weights_z_grad(s_idx)%sf(k, 1, 1) = c(s_idx, 1)
                end do
            end do
        end if

        do i = 1, nWeights_grad
            $:GPU_UPDATE(device='[weights_x_grad(i)%sf, weights_y_grad(i)%sf, weights_z_grad(i)%sf]')
        end do

    end subroutine s_compute_fornberg_fd_weights

    !> Fornberg (1988) finite-difference weights at x0 on the given stencil, for derivatives up to m_order.
    subroutine s_fornberg_weights(x0, stencil, npts, m_order, coeffs)

        $:GPU_ROUTINE(parallelism='[seq]')

        integer, intent(in)   :: npts  ! number of stencil points
        integer, intent(in)   :: m_order  ! highest derivative order
        real(wp), intent(in)  :: x0  ! evaluation point
        real(wp), intent(in)  :: stencil(npts)  ! stencil coordinates
        real(wp), intent(out) :: coeffs(npts,0:m_order)
        integer               :: i, j, k, mn
        real(wp)              :: c(5)  ! Running products and stencil offsets of the recurrence

        coeffs = 0.0_wp
        c(1) = 1.0_wp
        c(4) = stencil(1) - x0
        coeffs(1, 0) = 1.0_wp

        do i = 2, npts
            mn = min(i - 1, m_order)
            c(2) = 1.0_wp
            c(5) = c(4)
            c(4) = stencil(i) - x0

            do j = 1, i - 1
                c(3) = stencil(i) - stencil(j)
                c(2) = c(2)*c(3)

                if (j == i - 1) then
                    do k = mn, 1, -1
                        coeffs(i, k) = c(1)*(k*coeffs(i - 1, k - 1) - c(5)*coeffs(i - 1, k))/c(2)
                    end do
                    coeffs(i, 0) = -c(1)*c(5)*coeffs(i - 1, 0)/c(2)
                end if

                do k = mn, 1, -1
                    coeffs(j, k) = (c(4)*coeffs(j, k) - k*coeffs(j, k - 1))/c(3)
                end do
                coeffs(j, 0) = c(4)*coeffs(j, 0)/c(3)
            end do

            c(1) = c(2)
        end do

    end subroutine s_fornberg_weights

    !> One TVD RK stage in the Shu-Osher form of m_time_steppers, from stage value u_cur, step start u_n and stage RHS rhs.
    function f_lag_rk_stage(u_cur, u_n, rhs, coef) result(u_new)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)               :: u_cur, u_n, rhs
        real(wp), dimension(4), intent(in) :: coef
        real(wp)                           :: u_new

        u_new = (coef(1)*u_cur + coef(2)*u_n + coef(3)*dt*rhs)/coef(4)

    end function f_lag_rk_stage

    !> Apply this rank's domain boundary in direction dir to one Lagrangian entity at stage position pos (pos_prev before the
    !! stage): clamp it offset inside a wall, flag a periodic crossing in wrap_loc (-1 beg, 1 end, else 0), or clear keep once it
    !! leaves the domain.
    subroutine s_apply_lag_entity_domain_bc(dir, pos, pos_prev, offset, keep, wrap_loc)

        $:GPU_ROUTINE(function_name='s_apply_lag_entity_domain_bc',parallelism='[seq]', cray_inline=True)

        integer, intent(in)     :: dir
        real(wp), intent(inout) :: pos
        real(wp), intent(in)    :: pos_prev, offset
        integer, intent(inout)  :: keep
        integer, intent(out)    :: wrap_loc
        integer                 :: bc_beg, bc_end
        real(wp)                :: lo, hi

        if (dir == 1) then
            bc_beg = bc_x%beg; bc_end = bc_x%end
        else if (dir == 2) then
            bc_beg = bc_y%beg; bc_end = bc_y%end
        else
            bc_beg = bc_z%beg; bc_end = bc_z%end
        end if
        lo = pcomm_coords(dir)%beg
        hi = pcomm_coords(dir)%end

        wrap_loc = 0
        if (f_is_lag_entity_wall(bc_beg) .and. pos < lo + offset) then
            pos = lo + offset
        else if (f_is_lag_entity_wall(bc_end) .and. pos > hi - offset) then
            pos = hi - offset
        else if (bc_beg == BC_PERIODIC .and. pos < lo .and. pos_prev >= lo) then
            wrap_loc = -1
        else if (bc_end == BC_PERIODIC .and. pos > hi .and. pos_prev <= hi) then
            wrap_loc = 1
        else if (.not. (pos >= lo .and. pos < hi)) then
            keep = 0
        end if

    end subroutine s_apply_lag_entity_domain_bc

    !> True for the boundary types a Lagrangian entity cannot cross.
    logical function f_is_lag_entity_wall(bc) result(is_wall)

        $:GPU_ROUTINE(parallelism='[seq]')

        integer, intent(in) :: bc

        is_wall = bc == BC_REFLECTIVE .or. bc == BC_CHAR_SLIP_WALL .or. bc == BC_SLIP_WALL .or. bc == BC_NO_SLIP_WALL

    end function f_is_lag_entity_wall

    !> Shift the entities flagged by s_apply_lag_entity_domain_bc across the periodic domain, both time levels (host data).
    impure subroutine s_wrap_lag_entities_periodic(n_loc, wrap_loc, pos, pos_prev)

        integer, intent(in)                       :: n_loc
        integer, dimension(:,:), intent(in)       :: wrap_loc
        real(wp), dimension(:,:,:), intent(inout) :: pos, pos_prev
        real(wp)                                  :: offset
        integer                                   :: k, i

        do k = 1, n_loc
            do i = 1, num_dims
                if (wrap_loc(k, i) /= 0) then
                    offset = real(wrap_loc(k, i), wp)*(glb_bounds(i)%end - glb_bounds(i)%beg)
                    pos(k, i,1:2) = pos(k, i,1:2) - offset
                    pos_prev(k, i,1:2) = pos_prev(k, i,1:2) - offset
                end if
            end do
        end do

    end subroutine s_wrap_lag_entities_periodic

    !> Volume of cell (cellx, celly, cellz); in planar 2D the cell has depth charwidth.
    subroutine s_get_char_vol(cellx, celly, cellz, charwidth, Charvol)

        $:GPU_ROUTINE(function_name='s_get_char_vol',parallelism='[seq]', cray_inline=True)

        integer, intent(in)   :: cellx, celly, cellz
        real(wp), intent(in)  :: charwidth
        real(wp), intent(out) :: Charvol

        if (num_dims == 3) then
            Charvol = dx(cellx)*dy(celly)*dz(cellz)
        else
            if (cyl_coord) then
                Charvol = dx(cellx)*dy(celly)*y_cc(celly)*2._wp*pi
            else
                Charvol = dx(cellx)*dy(celly)*charwidth
            end if
        end if

    end subroutine s_get_char_vol

end module m_euler_lagrange
