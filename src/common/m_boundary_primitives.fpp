!>
!! @file
!! @brief Contains module m_boundary_primitives

!> @brief Per-cell noncharacteristic boundary condition primitives applied in the ghost cells
#:include 'case.fpp'
#:include 'macros.fpp'

module m_boundary_primitives

    use m_derived_types
    use m_global_parameters
    use m_constants
    use m_mpi_common

    implicit none

    type(scalar_field), dimension(:,:), allocatable :: bc_buffers
    $:GPU_DECLARE(create='[bc_buffers]')

    logical :: dirichlet_from_buffers = .false.
    $:GPU_DECLARE(create='[dirichlet_from_buffers]')

    !> f_vel_ramp per direction, applied to Dirichlet ghost velocities; set by the time stepper when a ramp is active
    real(wp), dimension(3) :: bc_vel_ramp = 1._wp
    $:GPU_DECLARE(create='[bc_vel_ramp]')

    public :: f_vel_ramp, bc_vel_ramp

contains

    !> Number of interior cells along the boundary-normal direction of a face: m, n or p.
    function f_bc_normal_extent(bc_dir) result(dim_end)

        $:GPU_ROUTINE(function_name='f_bc_normal_extent', parallelism='[seq]', cray_inline=True)
        integer, intent(in) :: bc_dir
        integer             :: dim_end

        if (bc_dir == 1) then
            dim_end = m
        else if (bc_dir == 2) then
            dim_end = n
        else
            dim_end = p
        end if

    end function f_bc_normal_extent

    subroutine s_bc_cell_index(bc_dir, k, l, nrm, ix, iy, iz)

        $:GPU_ROUTINE(function_name='s_bc_cell_index', parallelism='[seq]', cray_inline=True)
        integer, intent(in)  :: bc_dir, k, l, nrm
        integer, intent(out) :: ix, iy, iz

        if (bc_dir == 1) then
            ix = nrm; iy = k; iz = l
        else if (bc_dir == 2) then
            ix = k; iy = nrm; iz = l
        else
            ix = k; iy = l; iz = nrm
        end if

    end subroutine s_bc_cell_index

    !> Velocity scaling for a GRCBC inflow that is ramping up: unity unless `bc_[x,y,z]%vel_in_ramp` is set, so an unramped case is
    !! untouched. Takes the time as an argument so the caller can evaluate it on the device from `mytime`, which the time stepper
    !! already places there, rather than computing it on the host and copying the result every Runge-Kutta stage.
    pure function f_vel_ramp(tau, t0, frac0, t) result(f)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: tau, t0, frac0, t
        real(wp)             :: f

        if (tau > 0._wp) then
            f = frac0 + (1._wp - frac0)*0.5_wp*(1._wp + tanh(6._wp*(t - t0)/tau - 3._wp))
        else
            f = 1._wp
        end if

    end function f_vel_ramp

    !> Fill ghost cells by copying the nearest boundary cell value along the specified direction.
    subroutine s_ghost_cell_extrapolation(q_prim_vf, bc_dir, bc_loc, k, l, i, j, q_T_sf)

        $:GPU_ROUTINE(function_name='s_ghost_cell_extrapolation', parallelism='[seq]', cray_inline=True)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        integer, intent(in)                                    :: bc_dir, bc_loc
        integer, intent(in)                                    :: k, l, i, j
        type(scalar_field), optional, intent(inout)            :: q_T_sf
        integer                                                :: gx, gy, gz, sx, sy, sz, dim_end

        dim_end = f_bc_normal_extent(bc_dir)
        if (bc_loc == -1) then
            call s_bc_cell_index(bc_dir, k, l, -j, gx, gy, gz)
            call s_bc_cell_index(bc_dir, k, l, 0, sx, sy, sz)
        else
            call s_bc_cell_index(bc_dir, k, l, dim_end + j, gx, gy, gz)
            call s_bc_cell_index(bc_dir, k, l, dim_end, sx, sy, sz)
        end if

        if (i <= sys_size) then
            q_prim_vf(i)%sf(gx, gy, gz) = q_prim_vf(i)%sf(sx, sy, sz)
        else if (i == sys_size + 1 .and. (chemistry .or. heat_conduction) .and. present(q_T_sf)) then
            q_T_sf%sf(gx, gy, gz) = q_T_sf%sf(sx, sy, sz)
        end if

    end subroutine s_ghost_cell_extrapolation

    !> Apply reflective (symmetry) boundary conditions by mirroring primitive variables and flipping the normal velocity component.
    subroutine s_symmetry(q_prim_vf, bc_dir, bc_loc, k, l, i, j, pb_in, mv_in, q_T_sf)

        $:GPU_ROUTINE(parallelism='[seq]')
        type(scalar_field), dimension(sys_size), intent(inout)                                               :: q_prim_vf
        real(stp), optional, dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_in, mv_in
        integer, intent(in)                                                                                  :: bc_dir, bc_loc
        integer, intent(in)                                                                                  :: k, l, i, j
        type(scalar_field), optional, intent(inout)                                                          :: q_T_sf
        integer                                                                                              :: gx, gy, gz
        integer                                                                                              :: sx, sy, sz
        integer                                                                                              :: dim_end, nrm
        integer                                                                                              :: q, ib, qn
        logical                                                                                              :: flip

        dim_end = f_bc_normal_extent(bc_dir)
        if (bc_loc == -1) then
            call s_bc_cell_index(bc_dir, k, l, -j, gx, gy, gz)
            nrm = j - 1
        else
            call s_bc_cell_index(bc_dir, k, l, dim_end + j, gx, gy, gz)
            nrm = dim_end - (j - 1)
        end if
        call s_bc_cell_index(bc_dir, k, l, nrm, sx, sy, sz)

        if (i <= sys_size) then
            flip = (i == eqn_idx%mom%beg + bc_dir - 1)
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
        else if (i == sys_size + 1 .and. (chemistry .or. heat_conduction) .and. present(q_T_sf)) then
            q_T_sf%sf(gx, gy, gz) = q_T_sf%sf(sx, sy, sz)
        end if

        if (qbmm .and. .not. polytropic .and. present(pb_in) .and. present(mv_in)) then
            if (i <= nb*nnode) then
                ib = (i - 1)/nnode + 1
                qn = mod(i - 1, nnode) + 1
                pb_in(gx, gy, gz, qn, ib) = pb_in(sx, sy, sz, qn, ib)
                mv_in(gx, gy, gz, qn, ib) = mv_in(sx, sy, sz, qn, ib)
            end if
        end if

    end subroutine s_symmetry

    !> Apply periodic boundary conditions by copying values from the opposite domain boundary.
    subroutine s_periodic(q_prim_vf, bc_dir, bc_loc, k, l, i, j, pb_in, mv_in, q_T_sf)

        $:GPU_ROUTINE(parallelism='[seq]')
        type(scalar_field), dimension(sys_size), intent(inout)                                               :: q_prim_vf
        real(stp), optional, dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_in, mv_in
        integer, intent(in)                                                                                  :: bc_dir, bc_loc
        integer, intent(in)                                                                                  :: k, l, i, j
        type(scalar_field), optional, intent(inout)                                                          :: q_T_sf
        integer                                                                                              :: gx, gy, gz
        integer                                                                                              :: sx, sy, sz
        integer                                                                                              :: dim_end, nrm
        integer                                                                                              :: ib, qn

        dim_end = f_bc_normal_extent(bc_dir)
        if (bc_loc == -1) then
            call s_bc_cell_index(bc_dir, k, l, -j, gx, gy, gz)
            nrm = dim_end - (j - 1)
        else
            call s_bc_cell_index(bc_dir, k, l, dim_end + j, gx, gy, gz)
            nrm = j - 1
        end if
        call s_bc_cell_index(bc_dir, k, l, nrm, sx, sy, sz)

        if (i <= sys_size) then
            q_prim_vf(i)%sf(gx, gy, gz) = q_prim_vf(i)%sf(sx, sy, sz)
        else if (i == sys_size + 1 .and. (chemistry .or. heat_conduction) .and. present(q_T_sf)) then
            q_T_sf%sf(gx, gy, gz) = q_T_sf%sf(sx, sy, sz)
        end if

        if (qbmm .and. .not. polytropic .and. present(pb_in) .and. present(mv_in)) then
            if (i <= nb*nnode) then
                ib = (i - 1)/nnode + 1
                qn = mod(i - 1, nnode) + 1
                pb_in(gx, gy, gz, qn, ib) = pb_in(sx, sy, sz, qn, ib)
                mv_in(gx, gy, gz, qn, ib) = mv_in(sx, sy, sz, qn, ib)
            end if
        end if

    end subroutine s_periodic

    !> Apply axis boundary conditions for cylindrical coordinates by reflecting values across the axis with azimuthal phase shift.
    subroutine s_axis(q_prim_vf, pb_in, mv_in, k, l, i, j, q_T_sf)

        $:GPU_ROUTINE(parallelism='[seq]')
        type(scalar_field), dimension(sys_size), intent(inout)                                               :: q_prim_vf
        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), optional, intent(inout) :: pb_in, mv_in
        integer, intent(in)                                                                                  :: k, l, i, j
        integer                                                                                              :: lsrc, ib, qn
        type(scalar_field), optional, intent(inout)                                                          :: q_T_sf

        ! Reflection across the axis pairs each azimuthal cell with the one half a revolution away.
        if (z_cc(l) < pi) then
            lsrc = l + ((p + 1)/2)
        else
            lsrc = l - ((p + 1)/2)
        end if

        if (i <= sys_size) then
            ! Both velocity components transverse to the axis change sign; the axial one and every scalar do not.
            if (i == eqn_idx%mom%beg + 1 .or. i == eqn_idx%mom%end) then
                q_prim_vf(i)%sf(k, -j, l) = -q_prim_vf(i)%sf(k, j - 1, lsrc)
            else
                q_prim_vf(i)%sf(k, -j, l) = q_prim_vf(i)%sf(k, j - 1, lsrc)
            end if
        else if (i == sys_size + 1 .and. (chemistry .or. heat_conduction) .and. present(q_T_sf)) then
            q_T_sf%sf(k, -j, l) = q_T_sf%sf(k, j - 1, lsrc)
        end if

        if (qbmm .and. .not. polytropic .and. present(pb_in) .and. present(mv_in)) then
            if (i <= nb*nnode) then
                ib = (i - 1)/nnode + 1
                qn = mod(i - 1, nnode) + 1
                pb_in(k, -j, l, qn, ib) = pb_in(k, j - 1, lsrc, qn, ib)
                mv_in(k, -j, l, qn, ib) = mv_in(k, j - 1, lsrc, qn, ib)
            end if
        end if

    end subroutine s_axis

    !> Apply slip wall boundary conditions by extrapolating scalars and reflecting the wall-normal velocity component.
    subroutine s_slip_wall(q_prim_vf, bc_dir, bc_loc, k, l, i, j, q_T_sf)

        $:GPU_ROUTINE(function_name='s_slip_wall',parallelism='[seq]', cray_inline=True)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        integer, intent(in)                                    :: bc_dir, bc_loc
        integer, intent(in)                                    :: k, l, i, j
        type(scalar_field), optional, intent(inout)            :: q_T_sf
        integer                                                :: gx, gy, gz, mx, my, mz, ex, ey, ez, dim_end
        real(wp)                                               :: wall_vel, Twall
        logical                                                :: isothermal

        ! Ghost cell, its mirror image across the wall, and the edge cell the scalars extrapolate from.
        dim_end = f_bc_normal_extent(bc_dir)
        if (bc_loc == -1) then
            call s_bc_cell_index(bc_dir, k, l, -j, gx, gy, gz)
            call s_bc_cell_index(bc_dir, k, l, j - 1, mx, my, mz)
            call s_bc_cell_index(bc_dir, k, l, 0, ex, ey, ez)
        else
            call s_bc_cell_index(bc_dir, k, l, dim_end + j, gx, gy, gz)
            call s_bc_cell_index(bc_dir, k, l, dim_end - (j - 1), mx, my, mz)
            call s_bc_cell_index(bc_dir, k, l, dim_end, ex, ey, ez)
        end if

        ! A slip wall prescribes only its own normal velocity component, so one scalar per face.
        if (bc_dir == 1) then
            if (bc_loc == -1) then
                wall_vel = bc_x%vb1; isothermal = bc_x%isothermal_in; Twall = bc_x%Twall_in
            else
                wall_vel = bc_x%ve1; isothermal = bc_x%isothermal_out; Twall = bc_x%Twall_out
            end if
        else if (bc_dir == 2) then
            if (bc_loc == -1) then
                wall_vel = bc_y%vb2; isothermal = bc_y%isothermal_in; Twall = bc_y%Twall_in
            else
                wall_vel = bc_y%ve2; isothermal = bc_y%isothermal_out; Twall = bc_y%Twall_out
            end if
        else
            if (bc_loc == -1) then
                wall_vel = bc_z%vb3; isothermal = bc_z%isothermal_in; Twall = bc_z%Twall_in
            else
                wall_vel = bc_z%ve3; isothermal = bc_z%isothermal_out; Twall = bc_z%Twall_out
            end if
        end if

        if (i <= sys_size) then
            if (i == eqn_idx%mom%beg + bc_dir - 1) then
                q_prim_vf(i)%sf(gx, gy, gz) = -q_prim_vf(i)%sf(mx, my, mz) + 2._wp*wall_vel
            else
                q_prim_vf(i)%sf(gx, gy, gz) = q_prim_vf(i)%sf(ex, ey, ez)
            end if
        else if (i == sys_size + 1 .and. (chemistry .or. heat_conduction) .and. present(q_T_sf)) then
            if (isothermal) then
                q_T_sf%sf(gx, gy, gz) = 2._wp*Twall - q_T_sf%sf(mx, my, mz)
            else
                q_T_sf%sf(gx, gy, gz) = q_T_sf%sf(ex, ey, ez)
            end if
        end if

    end subroutine s_slip_wall

    !> Apply no-slip wall boundary conditions by reflecting and negating all velocity components at the wall.
    subroutine s_no_slip_wall(q_prim_vf, bc_dir, bc_loc, k, l, i, j, q_T_sf)

        $:GPU_ROUTINE(function_name='s_no_slip_wall',parallelism='[seq]', cray_inline=True)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        integer, intent(in)                                    :: bc_dir, bc_loc
        integer, intent(in)                                    :: k, l, i, j
        type(scalar_field), optional, intent(inout)            :: q_T_sf
        integer                                                :: gx, gy, gz, mx, my, mz, ex, ey, ez, dim_end, icomp
        real(wp), dimension(3)                                 :: wall_vel
        real(wp)                                               :: Twall
        logical                                                :: isothermal

        ! Ghost cell, its mirror image across the wall, and the edge cell the scalars extrapolate from.
        dim_end = f_bc_normal_extent(bc_dir)
        if (bc_loc == -1) then
            call s_bc_cell_index(bc_dir, k, l, -j, gx, gy, gz)
            call s_bc_cell_index(bc_dir, k, l, j - 1, mx, my, mz)
            call s_bc_cell_index(bc_dir, k, l, 0, ex, ey, ez)
        else
            call s_bc_cell_index(bc_dir, k, l, dim_end + j, gx, gy, gz)
            call s_bc_cell_index(bc_dir, k, l, dim_end - (j - 1), mx, my, mz)
            call s_bc_cell_index(bc_dir, k, l, dim_end, ex, ey, ez)
        end if

        ! A no-slip wall prescribes all three velocity components, each from its own face.
        if (bc_dir == 1) then
            if (bc_loc == -1) then
                wall_vel(1) = bc_x%vb1; wall_vel(2) = bc_x%vb2; wall_vel(3) = bc_x%vb3
                isothermal = bc_x%isothermal_in; Twall = bc_x%Twall_in
            else
                wall_vel(1) = bc_x%ve1; wall_vel(2) = bc_x%ve2; wall_vel(3) = bc_x%ve3
                isothermal = bc_x%isothermal_out; Twall = bc_x%Twall_out
            end if
        else if (bc_dir == 2) then
            if (bc_loc == -1) then
                wall_vel(1) = bc_y%vb1; wall_vel(2) = bc_y%vb2; wall_vel(3) = bc_y%vb3
                isothermal = bc_y%isothermal_in; Twall = bc_y%Twall_in
            else
                wall_vel(1) = bc_y%ve1; wall_vel(2) = bc_y%ve2; wall_vel(3) = bc_y%ve3
                isothermal = bc_y%isothermal_out; Twall = bc_y%Twall_out
            end if
        else
            if (bc_loc == -1) then
                wall_vel(1) = bc_z%vb1; wall_vel(2) = bc_z%vb2; wall_vel(3) = bc_z%vb3
                isothermal = bc_z%isothermal_in; Twall = bc_z%Twall_in
            else
                wall_vel(1) = bc_z%ve1; wall_vel(2) = bc_z%ve2; wall_vel(3) = bc_z%ve3
                isothermal = bc_z%isothermal_out; Twall = bc_z%Twall_out
            end if
        end if

        if (i <= sys_size) then
            icomp = i - eqn_idx%mom%beg + 1
            if (icomp >= 1 .and. icomp <= num_dims) then
                q_prim_vf(i)%sf(gx, gy, gz) = -q_prim_vf(i)%sf(mx, my, mz) + 2._wp*wall_vel(icomp)
            else
                q_prim_vf(i)%sf(gx, gy, gz) = q_prim_vf(i)%sf(ex, ey, ez)
            end if
        else if (i == sys_size + 1 .and. (chemistry .or. heat_conduction) .and. present(q_T_sf)) then
            if (isothermal) then
                q_T_sf%sf(gx, gy, gz) = 2._wp*Twall - q_T_sf%sf(mx, my, mz)
            else
                q_T_sf%sf(gx, gy, gz) = q_T_sf%sf(ex, ey, ez)
            end if
        end if

    end subroutine s_no_slip_wall

    !> Apply Dirichlet boundary conditions by prescribing ghost cell values from stored boundary buffers.
    subroutine s_dirichlet(q_prim_vf, bc_dir, bc_loc, k, l, i, j, q_T_sf)

        $:GPU_ROUTINE(function_name='s_dirichlet',parallelism='[seq]', cray_inline=True)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        integer, intent(in)                                    :: bc_dir, bc_loc
        integer, intent(in)                                    :: k, l, i, j
        type(scalar_field), optional, intent(inout)            :: q_T_sf
        integer                                                :: gx, gy, gz, dim_end, bc_iloc
        real(stp)                                              :: bc_val

        if (.not. dirichlet_from_buffers) then
            call s_ghost_cell_extrapolation(q_prim_vf, bc_dir, bc_loc, k, l, i, j, q_T_sf)
            return
        end if

        ! bc_buffers holds sys_size variables, plus the temperature under chemistry or heat conduction.
        if (i > sys_size) then
            if (.not. (i == sys_size + 1 .and. (chemistry .or. heat_conduction) .and. present(q_T_sf))) return
        end if

        dim_end = f_bc_normal_extent(bc_dir)
        bc_iloc = (bc_loc + 3)/2
        if (bc_loc == -1) then
            call s_bc_cell_index(bc_dir, k, l, -j, gx, gy, gz)
        else
            call s_bc_cell_index(bc_dir, k, l, dim_end + j, gx, gy, gz)
        end if

        ! The buffer carries the variable index in the slot of the face normal.
        bc_val = 0._stp
        if (bc_dir == 1) then
            bc_val = bc_buffers(1, bc_iloc)%sf(i, k, l)
        else if (bc_dir == 2) then
            #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
                bc_val = bc_buffers(2, bc_iloc)%sf(k, i, l)
            #:endif
        else
            #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                bc_val = bc_buffers(3, bc_iloc)%sf(k, l, i)
            #:endif
        end if

        ! Ramp the inflow ghost velocities; f_dir_vel is 1 for every other variable
        if (i <= sys_size) then
            q_prim_vf(i)%sf(gx, gy, gz) = bc_val*f_dir_vel(i, bc_dir)
        else
            q_T_sf%sf(gx, gy, gz) = bc_val
        end if

    end subroutine s_dirichlet

    !> Ramp factor for variable i of a Dirichlet ghost cell on a dir-normal face: bc_vel_ramp(dir) for velocities, else 1
    pure function f_dir_vel(i, dir) result(f)

        $:GPU_ROUTINE(parallelism='[seq]')
        integer, intent(in) :: i, dir
        real(wp)            :: f

        f = merge(bc_vel_ramp(dir), 1._wp, i >= eqn_idx%mom%beg .and. i <= eqn_idx%mom%end)

    end function f_dir_vel

    !> Extrapolate QBMM bubble pressure and mass-vapor variables into ghost cells by copying boundary values.
    subroutine s_qbmm_extrapolation(bc_dir, bc_loc, k, l, i, j, pb_in, mv_in)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(stp), optional, dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_in, mv_in
        integer, intent(in)                                                                                  :: bc_dir, bc_loc
        integer, intent(in)                                                                                  :: k, l, i, j
        integer                                                                                              :: gx, gy, gz
        integer                                                                                              :: sx, sy, sz
        integer                                                                                              :: dim_end, ib, qn

        ! `i` indexes the (bubble, quadrature node) pairs here, not the primitive variables.
        if (i > nb*nnode) return
        ib = (i - 1)/nnode + 1
        qn = mod(i - 1, nnode) + 1

        dim_end = f_bc_normal_extent(bc_dir)
        if (bc_loc == -1) then
            call s_bc_cell_index(bc_dir, k, l, -j, gx, gy, gz)
            call s_bc_cell_index(bc_dir, k, l, 0, sx, sy, sz)
        else
            call s_bc_cell_index(bc_dir, k, l, dim_end + j, gx, gy, gz)
            call s_bc_cell_index(bc_dir, k, l, dim_end, sx, sy, sz)
        end if

        pb_in(gx, gy, gz, qn, ib) = pb_in(sx, sy, sz, qn, ib)
        mv_in(gx, gy, gz, qn, ib) = mv_in(sx, sy, sz, qn, ib)

    end subroutine s_qbmm_extrapolation

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

    !> Apply periodic boundary conditions to the IGR Jacobian field by copying values from the opposite domain boundary.
    subroutine s_F_igr_periodic(jac_sf, bc_dir, bc_loc, k, l)

        $:GPU_ROUTINE(function_name='s_F_igr_periodic', parallelism='[seq]', cray_inline=True)
        type(scalar_field), dimension(1:), intent(inout) :: jac_sf
        integer, intent(in)                              :: bc_dir, bc_loc
        integer, intent(in)                              :: k, l
        integer                                          :: j

        if (bc_dir == 1) then
            if (bc_loc == -1) then
                do j = 1, buff_size
                    jac_sf(1)%sf(-j, k, l) = jac_sf(1)%sf(m - j + 1, k, l)
                end do
            else
                do j = 1, buff_size
                    jac_sf(1)%sf(m + j, k, l) = jac_sf(1)%sf(j - 1, k, l)
                end do
            end if
        else if (bc_dir == 2) then
            if (bc_loc == -1) then
                do j = 1, buff_size
                    jac_sf(1)%sf(k, -j, l) = jac_sf(1)%sf(k, n - j + 1, l)
                end do
            else
                do j = 1, buff_size
                    jac_sf(1)%sf(k, n + j, l) = jac_sf(1)%sf(k, j - 1, l)
                end do
            end if
        else
            if (bc_loc == -1) then
                do j = 1, buff_size
                    jac_sf(1)%sf(k, l, -j) = jac_sf(1)%sf(k, l, p - j + 1)
                end do
            else
                do j = 1, buff_size
                    jac_sf(1)%sf(k, l, p + j) = jac_sf(1)%sf(k, l, j - 1)
                end do
            end if
        end if

    end subroutine s_F_igr_periodic

    !> Apply reflective boundary conditions to the IGR Jacobian field by mirroring values across the boundary.
    subroutine s_F_igr_reflective(jac_sf, bc_dir, bc_loc, k, l)

        $:GPU_ROUTINE(function_name='s_F_igr_reflective', parallelism='[seq]', cray_inline=True)
        type(scalar_field), dimension(1:), intent(inout) :: jac_sf
        integer, intent(in)                              :: bc_dir, bc_loc
        integer, intent(in)                              :: k, l
        integer                                          :: j

        if (bc_dir == 1) then
            if (bc_loc == -1) then
                do j = 1, buff_size
                    jac_sf(1)%sf(-j, k, l) = jac_sf(1)%sf(j - 1, k, l)
                end do
            else
                do j = 1, buff_size
                    jac_sf(1)%sf(m + j, k, l) = jac_sf(1)%sf(m - (j - 1), k, l)
                end do
            end if
        else if (bc_dir == 2) then
            if (bc_loc == -1) then
                do j = 1, buff_size
                    jac_sf(1)%sf(k, -j, l) = jac_sf(1)%sf(k, j - 1, l)
                end do
            else
                do j = 1, buff_size
                    jac_sf(1)%sf(k, n + j, l) = jac_sf(1)%sf(k, n - (j - 1), l)
                end do
            end if
        else
            if (bc_loc == -1) then
                do j = 1, buff_size
                    jac_sf(1)%sf(k, l, -j) = jac_sf(1)%sf(k, l, j - 1)
                end do
            else
                do j = 1, buff_size
                    jac_sf(1)%sf(k, l, p + j) = jac_sf(1)%sf(k, l, p - (j - 1))
                end do
            end if
        end if

    end subroutine s_F_igr_reflective

    !> Extrapolate the IGR Jacobian field into ghost cells by copying boundary values.
    subroutine s_F_igr_ghost_cell_extrapolation(jac_sf, bc_dir, bc_loc, k, l)

        $:GPU_ROUTINE(function_name='s_F_igr_ghost_cell_extrapolation', parallelism='[seq]', cray_inline=True)
        type(scalar_field), dimension(1:), intent(inout) :: jac_sf
        integer, intent(in)                              :: bc_dir, bc_loc
        integer, intent(in)                              :: k, l
        integer                                          :: j

        if (bc_dir == 1) then
            if (bc_loc == -1) then
                do j = 1, buff_size
                    jac_sf(1)%sf(-j, k, l) = jac_sf(1)%sf(0, k, l)
                end do
            else
                do j = 1, buff_size
                    jac_sf(1)%sf(m + j, k, l) = jac_sf(1)%sf(m, k, l)
                end do
            end if
        else if (bc_dir == 2) then
            if (bc_loc == -1) then
                do j = 1, buff_size
                    jac_sf(1)%sf(k, -j, l) = jac_sf(1)%sf(k, 0, l)
                end do
            else
                do j = 1, buff_size
                    jac_sf(1)%sf(k, n + j, l) = jac_sf(1)%sf(k, n, l)
                end do
            end if
        else
            if (bc_loc == -1) then
                do j = 1, buff_size
                    jac_sf(1)%sf(k, l, -j) = jac_sf(1)%sf(k, l, 0)
                end do
            else
                do j = 1, buff_size
                    jac_sf(1)%sf(k, l, p + j) = jac_sf(1)%sf(k, l, p)
                end do
            end if
        end if

    end subroutine s_F_igr_ghost_cell_extrapolation

    subroutine s_beta_periodic(q_beta, kahan_comp, bc_dir, bc_loc, k, l, nvar)

        $:GPU_ROUTINE(function_name='s_beta_periodic', parallelism='[seq]', cray_inline=True)
        type(scalar_field), dimension(num_dims + 1), intent(inout) :: q_beta
        type(scalar_field), dimension(num_dims + 1), intent(inout) :: kahan_comp
        integer, intent(in)                                        :: bc_dir, bc_loc
        integer, intent(in)                                        :: k, l
        integer, intent(in)                                        :: nvar
        integer                                                    :: j, i
        real(wp)                                                   :: y_kahan, t_kahan

        if (bc_dir == 1) then  !< x-direction
            if (bc_loc == -1) then  ! bc%x%beg
                do i = 1, nvar
                    do j = -mapCells - 1, mapCells
                        ! Kahan-compensated addition of ghost to interior
                        y_kahan = real(q_beta(beta_vars(i))%sf(m + j + 1, k, l), &
                                       & kind=wp) + kahan_comp(beta_vars(i))%sf(m + j + 1, k, l) - kahan_comp(beta_vars(i))%sf(j, &
                                       & k, l)
                        t_kahan = real(q_beta(beta_vars(i))%sf(j, k, l), kind=wp) + y_kahan
                        kahan_comp(beta_vars(i))%sf(j, k, l) = (t_kahan - q_beta(beta_vars(i))%sf(j, k, l)) - y_kahan
                        q_beta(beta_vars(i))%sf(j, k, l) = t_kahan
                    end do
                end do
            else
                do i = 1, nvar
                    do j = -mapcells, mapcells + 1
                        q_beta(beta_vars(i))%sf(m + j, k, l) = q_beta(beta_vars(i))%sf(j - 1, k, l)
                        kahan_comp(beta_vars(i))%sf(m + j, k, l) = kahan_comp(beta_vars(i))%sf(j - 1, k, l)
                    end do
                end do
            end if
        else if (bc_dir == 2) then  !< y-direction
            if (bc_loc == -1) then  !< bc%y%beg
                do i = 1, nvar
                    do j = -mapcells - 1, mapcells
                        y_kahan = real(q_beta(beta_vars(i))%sf(k, n + j + 1, l), kind=wp) + kahan_comp(beta_vars(i))%sf(k, &
                                       & n + j + 1, l) - kahan_comp(beta_vars(i))%sf(k, j, l)
                        t_kahan = real(q_beta(beta_vars(i))%sf(k, j, l), kind=wp) + y_kahan
                        kahan_comp(beta_vars(i))%sf(k, j, l) = (t_kahan - q_beta(beta_vars(i))%sf(k, j, l)) - y_kahan
                        q_beta(beta_vars(i))%sf(k, j, l) = t_kahan
                    end do
                end do
            else
                do i = 1, nvar
                    do j = -mapcells, mapcells + 1
                        q_beta(beta_vars(i))%sf(k, n + j, l) = q_beta(beta_vars(i))%sf(k, j - 1, l)
                        kahan_comp(beta_vars(i))%sf(k, n + j, l) = kahan_comp(beta_vars(i))%sf(k, j - 1, l)
                    end do
                end do
            end if
        else if (bc_dir == 3) then  !< z-direction
            if (bc_loc == -1) then  !< bc%z%beg
                do i = 1, nvar
                    do j = -mapcells - 1, mapcells
                        y_kahan = real(q_beta(beta_vars(i))%sf(k, l, p + j + 1), kind=wp) + kahan_comp(beta_vars(i))%sf(k, l, &
                                       & p + j + 1) - kahan_comp(beta_vars(i))%sf(k, l, j)
                        t_kahan = real(q_beta(beta_vars(i))%sf(k, l, j), kind=wp) + y_kahan
                        kahan_comp(beta_vars(i))%sf(k, l, j) = (t_kahan - q_beta(beta_vars(i))%sf(k, l, j)) - y_kahan
                        q_beta(beta_vars(i))%sf(k, l, j) = t_kahan
                    end do
                end do
            else
                do i = 1, nvar
                    do j = -mapcells, mapcells + 1
                        q_beta(beta_vars(i))%sf(k, l, p + j) = q_beta(beta_vars(i))%sf(k, l, j - 1)
                        kahan_comp(beta_vars(i))%sf(k, l, p + j) = kahan_comp(beta_vars(i))%sf(k, l, j - 1)
                    end do
                end do
            end if
        end if

    end subroutine s_beta_periodic

    subroutine s_beta_extrapolation(q_beta, bc_dir, bc_loc, k, l, nvar)

        $:GPU_ROUTINE(function_name='s_beta_extrapolation', parallelism='[seq]', cray_inline=True)
        type(scalar_field), dimension(num_dims + 1), intent(inout) :: q_beta
        integer, intent(in)                                        :: bc_dir, bc_loc
        integer, intent(in)                                        :: k, l
        integer, intent(in)                                        :: nvar
        integer                                                    :: j, i

        if (bc_dir == 1) then  !< x-direction
            if (bc_loc == -1) then  ! bc%x%beg
                do i = 1, nvar
                    do j = 1, buff_size
                        q_beta(beta_vars(i))%sf(-j, k, l) = 0._wp
                    end do
                end do
            else
                do i = 1, nvar
                    do j = 1, buff_size
                        q_beta(beta_vars(i))%sf(m + j, k, l) = 0._wp
                    end do
                end do
            end if
        else if (bc_dir == 2) then  !< y-direction
            if (bc_loc == -1) then  !< bc%y%beg
                do i = 1, nvar
                    do j = 1, buff_size
                        q_beta(beta_vars(i))%sf(k, -j, l) = 0._wp
                    end do
                end do
            else
                do i = 1, nvar
                    do j = 1, buff_size
                        q_beta(beta_vars(i))%sf(k, n + j, l) = 0._wp
                    end do
                end do
            end if
        else if (bc_dir == 3) then  !< z-direction
            if (bc_loc == -1) then  !< bc%z%beg
                do i = 1, nvar
                    do j = 1, buff_size
                        q_beta(beta_vars(i))%sf(k, l, -j) = 0._wp
                    end do
                end do
            else  !< bc%z%end
                do i = 1, nvar
                    do j = 1, buff_size
                        q_beta(beta_vars(i))%sf(k, l, p + j) = 0._wp
                    end do
                end do
            end if
        end if

    end subroutine s_beta_extrapolation

    subroutine s_beta_reflective(q_beta, kahan_comp, bc_dir, bc_loc, k, l, nvar)

        $:GPU_ROUTINE(function_name='s_beta_reflective', parallelism='[seq]', cray_inline=True)
        type(scalar_field), dimension(num_dims + 1), intent(inout) :: q_beta
        type(scalar_field), dimension(num_dims + 1), intent(inout) :: kahan_comp
        integer, intent(in)                                        :: bc_dir, bc_loc
        integer, intent(in)                                        :: k, l
        integer, intent(in)                                        :: nvar
        integer                                                    :: j, i
        real(wp)                                                   :: y_kahan, t_kahan

        ! Reflective BC for void fraction: 1) Fold ghost-cell contributions back onto their mirror interior cells (Kahan) 2) Set
        ! ghost cells = mirror of (now-folded) interior values

        if (bc_dir == 1) then  !< x-direction
            if (bc_loc == -1) then  ! bc%x%beg
                do i = 1, nvar
                    do j = 1, mapCells + 1
                        y_kahan = real(q_beta(beta_vars(i))%sf(-j, k, l), kind=wp) + kahan_comp(beta_vars(i))%sf(-j, k, &
                                       & l) - kahan_comp(beta_vars(i))%sf(j - 1, k, l)
                        t_kahan = real(q_beta(beta_vars(i))%sf(j - 1, k, l), kind=wp) + y_kahan
                        kahan_comp(beta_vars(i))%sf(j - 1, k, l) = (t_kahan - q_beta(beta_vars(i))%sf(j - 1, k, l)) - y_kahan
                        q_beta(beta_vars(i))%sf(j - 1, k, l) = t_kahan
                    end do
                    do j = 1, mapCells + 1
                        q_beta(beta_vars(i))%sf(-j, k, l) = q_beta(beta_vars(i))%sf(j - 1, k, l)
                        kahan_comp(beta_vars(i))%sf(-j, k, l) = kahan_comp(beta_vars(i))%sf(j - 1, k, l)
                    end do
                end do
            else  !< bc%x%end
                do i = 1, nvar
                    do j = 1, mapCells + 1
                        y_kahan = real(q_beta(beta_vars(i))%sf(m + j, k, l), kind=wp) + kahan_comp(beta_vars(i))%sf(m + j, k, &
                                       & l) - kahan_comp(beta_vars(i))%sf(m - (j - 1), k, l)
                        t_kahan = real(q_beta(beta_vars(i))%sf(m - (j - 1), k, l), kind=wp) + y_kahan
                        kahan_comp(beta_vars(i))%sf(m - (j - 1), k, l) = (t_kahan - q_beta(beta_vars(i))%sf(m - (j - 1), k, &
                                   & l)) - y_kahan
                        q_beta(beta_vars(i))%sf(m - (j - 1), k, l) = t_kahan
                    end do
                    do j = 1, mapCells + 1
                        q_beta(beta_vars(i))%sf(m + j, k, l) = q_beta(beta_vars(i))%sf(m - (j - 1), k, l)
                        kahan_comp(beta_vars(i))%sf(m + j, k, l) = kahan_comp(beta_vars(i))%sf(m - (j - 1), k, l)
                    end do
                end do
            end if
        else if (bc_dir == 2) then  !< y-direction
            if (bc_loc == -1) then  !< bc%y%beg
                do i = 1, nvar
                    do j = 1, mapCells + 1
                        y_kahan = real(q_beta(beta_vars(i))%sf(k, -j, l), kind=wp) + kahan_comp(beta_vars(i))%sf(k, -j, &
                                       & l) - kahan_comp(beta_vars(i))%sf(k, j - 1, l)
                        t_kahan = real(q_beta(beta_vars(i))%sf(k, j - 1, l), kind=wp) + y_kahan
                        kahan_comp(beta_vars(i))%sf(k, j - 1, l) = (t_kahan - q_beta(beta_vars(i))%sf(k, j - 1, l)) - y_kahan
                        q_beta(beta_vars(i))%sf(k, j - 1, l) = t_kahan
                    end do
                    do j = 1, mapCells + 1
                        q_beta(beta_vars(i))%sf(k, -j, l) = q_beta(beta_vars(i))%sf(k, j - 1, l)
                        kahan_comp(beta_vars(i))%sf(k, -j, l) = kahan_comp(beta_vars(i))%sf(k, j - 1, l)
                    end do
                end do
            else  !< bc%y%end
                do i = 1, nvar
                    do j = 1, mapCells + 1
                        y_kahan = real(q_beta(beta_vars(i))%sf(k, n + j, l), kind=wp) + kahan_comp(beta_vars(i))%sf(k, n + j, &
                                       & l) - kahan_comp(beta_vars(i))%sf(k, n - (j - 1), l)
                        t_kahan = real(q_beta(beta_vars(i))%sf(k, n - (j - 1), l), kind=wp) + y_kahan
                        kahan_comp(beta_vars(i))%sf(k, n - (j - 1), l) = (t_kahan - q_beta(beta_vars(i))%sf(k, n - (j - 1), &
                                   & l)) - y_kahan
                        q_beta(beta_vars(i))%sf(k, n - (j - 1), l) = t_kahan
                    end do
                    do j = 1, mapCells + 1
                        q_beta(beta_vars(i))%sf(k, n + j, l) = q_beta(beta_vars(i))%sf(k, n - (j - 1), l)
                        kahan_comp(beta_vars(i))%sf(k, n + j, l) = kahan_comp(beta_vars(i))%sf(k, n - (j - 1), l)
                    end do
                end do
            end if
        else if (bc_dir == 3) then  !< z-direction
            if (bc_loc == -1) then  !< bc%z%beg
                do i = 1, nvar
                    do j = 1, mapCells + 1
                        y_kahan = real(q_beta(beta_vars(i))%sf(k, l, -j), kind=wp) + kahan_comp(beta_vars(i))%sf(k, l, &
                                       & -j) - kahan_comp(beta_vars(i))%sf(k, l, j - 1)
                        t_kahan = real(q_beta(beta_vars(i))%sf(k, l, j - 1), kind=wp) + y_kahan
                        kahan_comp(beta_vars(i))%sf(k, l, j - 1) = (t_kahan - q_beta(beta_vars(i))%sf(k, l, j - 1)) - y_kahan
                        q_beta(beta_vars(i))%sf(k, l, j - 1) = t_kahan
                    end do
                    do j = 1, mapCells + 1
                        q_beta(beta_vars(i))%sf(k, l, -j) = q_beta(beta_vars(i))%sf(k, l, j - 1)
                        kahan_comp(beta_vars(i))%sf(k, l, -j) = kahan_comp(beta_vars(i))%sf(k, l, j - 1)
                    end do
                end do
            else  !< bc%z%end
                do i = 1, nvar
                    do j = 1, mapCells + 1
                        y_kahan = real(q_beta(beta_vars(i))%sf(k, l, p + j), kind=wp) + kahan_comp(beta_vars(i))%sf(k, l, &
                                       & p + j) - kahan_comp(beta_vars(i))%sf(k, l, p - (j - 1))
                        t_kahan = real(q_beta(beta_vars(i))%sf(k, l, p - (j - 1)), kind=wp) + y_kahan
                        kahan_comp(beta_vars(i))%sf(k, l, p - (j - 1)) = (t_kahan - q_beta(beta_vars(i))%sf(k, l, &
                                   & p - (j - 1))) - y_kahan
                        q_beta(beta_vars(i))%sf(k, l, p - (j - 1)) = t_kahan
                    end do
                    do j = 1, mapCells + 1
                        q_beta(beta_vars(i))%sf(k, l, p + j) = q_beta(beta_vars(i))%sf(k, l, p - (j - 1))
                        kahan_comp(beta_vars(i))%sf(k, l, p + j) = kahan_comp(beta_vars(i))%sf(k, l, p - (j - 1))
                    end do
                end do
            end if
        end if

    end subroutine s_beta_reflective

end module m_boundary_primitives
