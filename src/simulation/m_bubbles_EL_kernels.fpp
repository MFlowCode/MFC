!>
!! @file
!! @brief Contains module @ref m_bubbles_el_kernels "m_bubbles_EL_kernels"

#:include 'macros.fpp'

!> @brief Kernel functions (Gaussian, delta) that smear Lagrangian bubble effects onto the Eulerian grid
module m_bubbles_EL_kernels

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    implicit none

    ! Cell-centered pressure gradients (precomputed for translational motion)
    real(wp), allocatable, dimension(:, :, :) :: grad_p_x, grad_p_y, grad_p_z
    $:GPU_DECLARE(create='[grad_p_x, grad_p_y, grad_p_z]')

    ! Finite-difference coefficients for pressure gradient computation
    real(wp), allocatable, dimension(:, :) :: fd_coeff_x_pgrad
    real(wp), allocatable, dimension(:, :) :: fd_coeff_y_pgrad
    real(wp), allocatable, dimension(:, :) :: fd_coeff_z_pgrad
    $:GPU_DECLARE(create='[fd_coeff_x_pgrad, fd_coeff_y_pgrad, fd_coeff_z_pgrad]')

    ! Cell list for bubble-to-cell mapping (rebuilt each RK stage before smearing)
    integer, allocatable, dimension(:, :, :) :: cell_list_start  ! (0:m, 0:n, 0:p)
    integer, allocatable, dimension(:, :, :) :: cell_list_count  ! (0:m, 0:n, 0:p)
    integer, allocatable, dimension(:) :: cell_list_idx          ! (1:nBubs_glb) sorted bubble indices
    $:GPU_DECLARE(create='[cell_list_start, cell_list_count, cell_list_idx]')

contains

    !> The purpose of this subroutine is to smear the strength of the lagrangian
            !!      bubbles into the Eulerian framework using different approaches.
            !! @param nBubs Number of lagrangian bubbles in the current domain
            !! @param lbk_rad Radius of the bubbles
            !! @param lbk_vel Interface velocity of the bubbles
            !! @param lbk_s Computational coordinates of the bubbles
            !! @param lbk_pos Spatial coordinates of the bubbles
            !! @param updatedvar Eulerian variable to be updated
    subroutine s_smoothfunction(nBubs, lbk_rad, lbk_vel, lbk_s, lbk_pos, updatedvar, kcomp)

        integer, intent(in) :: nBubs
        real(wp), dimension(1:lag_params%nBubs_glb, 1:3, 1:2), intent(in) :: lbk_s, lbk_pos
        real(wp), dimension(1:lag_params%nBubs_glb, 1:2), intent(in) :: lbk_rad, lbk_vel
        type(scalar_field), dimension(:), intent(inout) :: updatedvar
        type(scalar_field), dimension(:), intent(inout) :: kcomp

        smoothfunc:select case(lag_params%smooth_type)
        case (1)
        call s_gaussian(nBubs, lbk_rad, lbk_vel, lbk_s, lbk_pos, updatedvar, kcomp)
        case (2)
        call s_deltafunc(nBubs, lbk_rad, lbk_vel, lbk_s, updatedvar, kcomp)
        end select smoothfunc

    end subroutine s_smoothfunction

    !> Builds a sorted cell list mapping each interior cell (0:m,0:n,0:p) to its
    !!      resident bubbles. Uses a counting-sort on the host (O(nBubs + N_cells)).
    !!      Must be called before s_gaussian each RK stage.
    !! @param nBubs Number of lagrangian bubbles in the current domain
    !! @param lbk_s Computational coordinates of the bubbles
    subroutine s_build_cell_list(nBubs, lbk_s)

        integer, intent(in) :: nBubs
        real(wp), dimension(1:lag_params%nBubs_glb, 1:3, 1:2), intent(in) :: lbk_s

        integer :: l, ci, cj, ck, idx
        real(wp), dimension(3) :: s_coord

        ! Bring current bubble positions to host
        $:GPU_UPDATE(host='[lbk_s]')

        ! Pass 1: zero counts and count bubbles per cell
        cell_list_count = 0
        do l = 1, nBubs
            s_coord(1:3) = lbk_s(l, 1:3, 2)
            ci = int(s_coord(1))
            cj = int(s_coord(2))
            ck = int(s_coord(3))
            ! Clamp to interior (bubbles should already be in [0:m,0:n,0:p])
            ci = max(0, min(ci, m))
            cj = max(0, min(cj, n))
            ck = max(0, min(ck, p))
            cell_list_count(ci, cj, ck) = cell_list_count(ci, cj, ck) + 1
        end do

        ! Prefix sum to compute start indices (1-based into cell_list_idx)
        idx = 1
        do ck = 0, p
            do cj = 0, n
                do ci = 0, m
                    cell_list_start(ci, cj, ck) = idx
                    idx = idx + cell_list_count(ci, cj, ck)
                end do
            end do
        end do

        ! Pass 2: place bubble indices into cell_list_idx
        ! Temporarily reuse cell_list_count as a running offset
        cell_list_count = 0
        do l = 1, nBubs
            s_coord(1:3) = lbk_s(l, 1:3, 2)
            ci = int(s_coord(1))
            cj = int(s_coord(2))
            ck = int(s_coord(3))
            ci = max(0, min(ci, m))
            cj = max(0, min(cj, n))
            ck = max(0, min(ck, p))
            cell_list_idx(cell_list_start(ci, cj, ck) + cell_list_count(ci, cj, ck)) = l
            cell_list_count(ci, cj, ck) = cell_list_count(ci, cj, ck) + 1
        end do

        ! Send cell list arrays to GPU
        $:GPU_UPDATE(device='[cell_list_start, cell_list_count, cell_list_idx]')

    end subroutine s_build_cell_list

    !> Cell-centric delta-function smearing using the cell list (no GPU atomics).
    !!      Each bubble only affects the cell it resides in. The outer GPU loop
    !!      iterates over interior cells and sums contributions from resident bubbles.
    subroutine s_deltafunc(nBubs, lbk_rad, lbk_vel, lbk_s, updatedvar, kcomp)

        integer, intent(in) :: nBubs
        real(wp), dimension(1:lag_params%nBubs_glb, 1:3, 1:2), intent(in) :: lbk_s
        real(wp), dimension(1:lag_params%nBubs_glb, 1:2), intent(in) :: lbk_rad, lbk_vel
        type(scalar_field), dimension(:), intent(inout) :: updatedvar
        type(scalar_field), dimension(:), intent(inout) :: kcomp

        real(wp) :: strength_vel, strength_vol
        real(wp) :: volpart, Vol
        real(wp) :: y_kahan, t_kahan
        integer :: i, j, k, lb, bub_idx

        $:GPU_PARALLEL_LOOP(collapse=3, private='[i,j,k,lb,bub_idx,volpart,Vol,strength_vel,strength_vol,y_kahan,t_kahan]')
        do k = 0, p
            do j = 0, n
                do i = 0, m

                    ! Cell volume
                    if (num_dims == 2) then
                        Vol = dx(i)*dy(j)*lag_params%charwidth
                        if (cyl_coord) Vol = dx(i)*dy(j)*y_cc(j)*2._wp*pi
                    else
                        Vol = dx(i)*dy(j)*dz(k)
                    end if

                    ! Loop over bubbles in this cell
                    $:GPU_LOOP(parallelism='[seq]')
                    do lb = cell_list_start(i, j, k), &
                        cell_list_start(i, j, k) + cell_list_count(i, j, k) - 1

                        bub_idx = cell_list_idx(lb)

                        volpart = 4._wp/3._wp*pi*lbk_rad(bub_idx, 2)**3._wp
                        strength_vol = volpart
                        strength_vel = 4._wp*pi*lbk_rad(bub_idx, 2)**2._wp*lbk_vel(bub_idx, 2)

                        ! Kahan summation for void fraction
                        y_kahan = real(strength_vol/Vol, kind=wp) - kcomp(1)%sf(i, j, k)
                        t_kahan = updatedvar(1)%sf(i, j, k) + y_kahan
                        kcomp(1)%sf(i, j, k) = (t_kahan - updatedvar(1)%sf(i, j, k)) - y_kahan
                        updatedvar(1)%sf(i, j, k) = t_kahan

                        ! Kahan summation for time derivative of void fraction
                        y_kahan = real(strength_vel/Vol, kind=wp) - kcomp(2)%sf(i, j, k)
                        t_kahan = updatedvar(2)%sf(i, j, k) + y_kahan
                        kcomp(2)%sf(i, j, k) = (t_kahan - updatedvar(2)%sf(i, j, k)) - y_kahan
                        updatedvar(2)%sf(i, j, k) = t_kahan

                        ! Product of two smeared functions
                        if (lag_params%cluster_type >= 4) then
                            y_kahan = real((strength_vol*strength_vel)/Vol, kind=wp) - kcomp(5)%sf(i, j, k)
                            t_kahan = updatedvar(5)%sf(i, j, k) + y_kahan
                            kcomp(5)%sf(i, j, k) = (t_kahan - updatedvar(5)%sf(i, j, k)) - y_kahan
                            updatedvar(5)%sf(i, j, k) = t_kahan
                        end if
                    end do

                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_deltafunc

    !> Cell-centric gaussian smearing using the cell list (no GPU atomics).
    !!      Each grid cell accumulates contributions from nearby bubbles looked up
    !!      via cell_list_start/count/idx.
    subroutine s_gaussian(nBubs, lbk_rad, lbk_vel, lbk_s, lbk_pos, updatedvar, kcomp)

        integer, intent(in) :: nBubs
        real(wp), dimension(1:lag_params%nBubs_glb, 1:3, 1:2), intent(in) :: lbk_s, lbk_pos
        real(wp), dimension(1:lag_params%nBubs_glb, 1:2), intent(in) :: lbk_rad, lbk_vel
        type(scalar_field), dimension(:), intent(inout) :: updatedvar
        type(scalar_field), dimension(:), intent(inout) :: kcomp

        real(wp), dimension(3) :: center, nodecoord, s_coord
        integer, dimension(3) :: cell, cellijk
        real(wp) :: stddsv, volpart
        real(wp) :: strength_vel, strength_vol
        real(wp) :: func, func2
        real(wp) :: y_kahan, t_kahan
        integer :: i, j, k, di, dj, dk, lb, bub_idx
        integer :: di_beg, di_end, dj_beg, dj_end, dk_beg, dk_end
        integer :: smear_x_beg, smear_x_end
        integer :: smear_y_beg, smear_y_end
        integer :: smear_z_beg, smear_z_end

        ! Extended grid range for smearing (includes buffer cells for MPI communication)
        smear_x_beg = -mapCells - 1
        smear_x_end = m + mapCells + 1
        smear_y_beg = merge(-mapCells - 1, 0, n > 0)
        smear_y_end = merge(n + mapCells + 1, n, n > 0)
        smear_z_beg = merge(-mapCells - 1, 0, p > 0)
        smear_z_end = merge(p + mapCells + 1, p, p > 0)

        $:GPU_PARALLEL_LOOP(collapse=3, &
            & private='[i,j,k,di,dj,dk,lb,bub_idx,center,nodecoord,s_coord,cell,cellijk,stddsv,volpart,strength_vel,strength_vol,func,func2,y_kahan,t_kahan,di_beg,di_end,dj_beg,dj_end,dk_beg,dk_end]', &
            & copyin='[smear_x_beg,smear_x_end,smear_y_beg,smear_y_end,smear_z_beg,smear_z_end]')
        do k = smear_z_beg, smear_z_end
            do j = smear_y_beg, smear_y_end
                do i = smear_x_beg, smear_x_end

                    cellijk(1) = i
                    cellijk(2) = j
                    cellijk(3) = k

                    nodecoord(1) = x_cc(i)
                    nodecoord(2) = y_cc(j)
                    nodecoord(3) = 0._wp
                    if (p > 0) nodecoord(3) = z_cc(k)

                    ! Neighbor cell range clamped to interior [0:m, 0:n, 0:p]
                    di_beg = max(i - mapCells, 0)
                    di_end = min(i + mapCells, m)
                    dj_beg = max(j - mapCells, 0)
                    dj_end = min(j + mapCells, n)
                    dk_beg = max(k - mapCells, 0)
                    dk_end = min(k + mapCells, p)

                    $:GPU_LOOP(parallelism='[seq]')
                    do dk = dk_beg, dk_end
                        $:GPU_LOOP(parallelism='[seq]')
                        do dj = dj_beg, dj_end
                            $:GPU_LOOP(parallelism='[seq]')
                            do di = di_beg, di_end
                                $:GPU_LOOP(parallelism='[seq]')
                                do lb = cell_list_start(di, dj, dk), &
                                    cell_list_start(di, dj, dk) + cell_list_count(di, dj, dk) - 1

                                    bub_idx = cell_list_idx(lb)

                                    ! Bubble properties
                                    volpart = 4._wp/3._wp*pi*lbk_rad(bub_idx, 2)**3._wp
                                    s_coord(1:3) = lbk_s(bub_idx, 1:3, 2)
                                    call s_get_cell(s_coord, cell)
                                    call s_compute_stddsv(cell, volpart, stddsv)

                                    strength_vol = volpart
                                    strength_vel = 4._wp*pi*lbk_rad(bub_idx, 2)**2._wp*lbk_vel(bub_idx, 2)

                                    center(1:2) = lbk_pos(bub_idx, 1:2, 2)
                                    center(3) = 0._wp
                                    if (p > 0) center(3) = lbk_pos(bub_idx, 3, 2)

                                    call s_applygaussian(center, cellijk, nodecoord, stddsv, 0._wp, func)

                                    ! Kahan summation for void fraction
                                    y_kahan = real(func*strength_vol, kind=wp) - kcomp(1)%sf(i, j, k)
                                    t_kahan = updatedvar(1)%sf(i, j, k) + y_kahan
                                    kcomp(1)%sf(i, j, k) = (t_kahan - updatedvar(1)%sf(i, j, k)) - y_kahan
                                    updatedvar(1)%sf(i, j, k) = t_kahan

                                    ! Kahan summation for time derivative of void fraction
                                    y_kahan = real(func*strength_vel, kind=wp) - kcomp(2)%sf(i, j, k)
                                    t_kahan = updatedvar(2)%sf(i, j, k) + y_kahan
                                    kcomp(2)%sf(i, j, k) = (t_kahan - updatedvar(2)%sf(i, j, k)) - y_kahan
                                    updatedvar(2)%sf(i, j, k) = t_kahan

                                    if (lag_params%cluster_type >= 4) then
                                        call s_applygaussian(center, cellijk, nodecoord, stddsv, 1._wp, func2)
                                        y_kahan = real(func2*strength_vol*strength_vel, kind=wp) - kcomp(5)%sf(i, j, k)
                                        t_kahan = updatedvar(5)%sf(i, j, k) + y_kahan
                                        kcomp(5)%sf(i, j, k) = (t_kahan - updatedvar(5)%sf(i, j, k)) - y_kahan
                                        updatedvar(5)%sf(i, j, k) = t_kahan
                                    end if

                                end do
                            end do
                        end do
                    end do

                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_gaussian

    !> The purpose of this subroutine is to apply the gaussian kernel function for each bubble (Maeda and Colonius, 2018)).
    subroutine s_applygaussian(center, cellaux, nodecoord, stddsv, strength_idx, func)
        $:GPU_ROUTINE(function_name='s_applygaussian',parallelism='[seq]', &
            & cray_inline=True)

        real(wp), dimension(3), intent(in) :: center
        integer, dimension(3), intent(in) :: cellaux
        real(wp), dimension(3), intent(in) :: nodecoord
        real(wp), intent(in) :: stddsv
        real(wp), intent(in) :: strength_idx
        real(wp), intent(out) :: func
        integer :: i

        real(wp) :: distance
        real(wp) :: theta, dtheta, L2, dzp, Lz2, zc
        real(wp) :: Nr, Nr_count

        distance = sqrt((center(1) - nodecoord(1))**2._wp + (center(2) - nodecoord(2))**2._wp + (center(3) - nodecoord(3))**2._wp)

        if (num_dims == 3) then
            !< 3D gaussian function
            func = exp(-0.5_wp*(distance/stddsv)**2._wp)/(sqrt(2._wp*pi)*stddsv)**3._wp
        else
            if (cyl_coord) then
                !< 2D cylindrical function:
                ! We smear particles in the azimuthal direction for given r
                theta = 0._wp
                Nr = ceiling(2._wp*pi*nodecoord(2)/(y_cb(cellaux(2)) - y_cb(cellaux(2) - 1)))
                dtheta = 2._wp*pi/Nr
                L2 = center(2)**2._wp + nodecoord(2)**2._wp - 2._wp*center(2)*nodecoord(2)*cos(theta)
                distance = sqrt((center(1) - nodecoord(1))**2._wp + L2)
                ! Factor 2._wp is for symmetry (upper half of the 2D field (+r) is considered)
                func = dtheta/2._wp/pi*exp(-0.5_wp*(distance/stddsv)**2._wp)/(sqrt(2._wp*pi)*stddsv)**3._wp
                Nr_count = 0._wp
                do while (Nr_count < Nr - 1._wp)
                    Nr_count = Nr_count + 1._wp
                    theta = Nr_count*dtheta
                    ! trigonometric relation
                    L2 = center(2)**2._wp + nodecoord(2)**2._wp - 2._wp*center(2)*nodecoord(2)*cos(theta)
                    distance = sqrt((center(1) - nodecoord(1))**2._wp + L2)
                    ! nodecoord(2)*dtheta is the azimuthal width of the cell
                    func = func + &
                           dtheta/2._wp/pi*exp(-0.5_wp*(distance/stddsv)**2._wp)/(sqrt(2._wp*pi)*stddsv)**(3._wp*(strength_idx + 1._wp))
                end do
            else
                !< 2D cartesian function: Equation (48) from Madea and Colonius 2018
                ! We smear particles considering a virtual depth (lag_params%charwidth) with lag_params%charNz cells
                dzp = (lag_params%charwidth/(lag_params%charNz + 1._wp))

                func = 0._wp
                do i = 0, lag_params%charNz
                    zc = (-lag_params%charwidth/2._wp + dzp*(0.5_wp + i)) ! Center of virtual cell i in z-direction
                    Lz2 = (center(3) - zc)**2._wp
                    distance = sqrt((center(1) - nodecoord(1))**2._wp + (center(2) - nodecoord(2))**2._wp + Lz2)
                    func = func + dzp/lag_params%charwidth*exp(-0.5_wp*(distance/stddsv)**2._wp)/(sqrt(2._wp*pi)*stddsv)**3._wp
                end do
            end if
        end if

    end subroutine s_applygaussian

    !> Calculates the standard deviation of the bubble being smeared in the Eulerian framework.
            !! @param cell Cell where the bubble is located
            !! @param volpart Volume of the bubble
            !! @param stddsv Standard deviaton
    subroutine s_compute_stddsv(cell, volpart, stddsv)
        $:GPU_ROUTINE(function_name='s_compute_stddsv',parallelism='[seq]', &
            & cray_inline=True)

        integer, dimension(3), intent(in) :: cell
        real(wp), intent(in) :: volpart
        real(wp), intent(out) :: stddsv

        real(wp) :: chardist, charvol
        real(wp) :: rad

        !< Compute characteristic distance
        chardist = sqrt(dx(cell(1))*dy(cell(2)))
        if (p > 0) chardist = (dx(cell(1))*dy(cell(2))*dz(cell(3)))**(1._wp/3._wp)

        !< Compute characteristic volume
        if (p > 0) then
            charvol = dx(cell(1))*dy(cell(2))*dz(cell(3))
        else
            if (cyl_coord) then
                charvol = dx(cell(1))*dy(cell(2))*y_cc(cell(2))*2._wp*pi
            else
                charvol = dx(cell(1))*dy(cell(2))*lag_params%charwidth
            end if
        end if

        !< Compute Standard deviaton
        if ((volpart/charvol) > 0.5_wp*lag_params%valmaxvoid .or. (lag_params%smooth_type == 1)) then
            rad = (3._wp*volpart/(4._wp*pi))**(1._wp/3._wp)
            stddsv = 1._wp*lag_params%epsilonb*max(chardist, rad)
        else
            stddsv = 0._wp
        end if

    end subroutine s_compute_stddsv

    !> The purpose of this procedure is to calculate the characteristic cell volume
            !! @param cellx x-direction cell index
            !! @param celly y-direction cell index
            !! @param cellz z-direction cell index
            !! @param Charvol Characteristic volume
    subroutine s_get_char_vol(cellx, celly, cellz, Charvol)
        $:GPU_ROUTINE(function_name='s_get_char_vol',parallelism='[seq]', &
            & cray_inline=True)

        integer, intent(in) :: cellx, celly, cellz
        real(wp), intent(out) :: Charvol

        if (p > 0) then
            Charvol = dx(cellx)*dy(celly)*dz(cellz)
        else
            if (cyl_coord) then
                Charvol = dx(cellx)*dy(celly)*y_cc(celly)*2._wp*pi
            else
                Charvol = dx(cellx)*dy(celly)*lag_params%charwidth
            end if
        end if

    end subroutine s_get_char_vol

    !> This subroutine transforms the computational coordinates of the bubble from
            !!      real type into integer.
            !! @param s_cell Computational coordinates of the bubble, real type
            !! @param get_cell Computational coordinates of the bubble, integer type
    subroutine s_get_cell(s_cell, get_cell)
        $:GPU_ROUTINE(function_name='s_get_cell',parallelism='[seq]', &
            & cray_inline=True)

        real(wp), dimension(3), intent(in) :: s_cell
        integer, dimension(3), intent(out) :: get_cell
        integer :: i

        get_cell(:) = int(s_cell(:))
        do i = 1, num_dims
            if (s_cell(i) < 0._wp) get_cell(i) = get_cell(i) - 1
        end do

    end subroutine s_get_cell

    !> Precomputes cell-centered pressure gradients (dp/dx, dp/dy, dp/dz) at all cell centers
        !!      using finite-difference coefficients of the specified order. This avoids
        !!      scattered memory accesses to the pressure field when computing translational
        !!      bubble forces.
        !! @param q_prim_vf Primitive variables (pressure is at index E_idx)
    subroutine s_compute_pressure_gradients(q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf

        integer :: i, j, k, r

        ! dp/dx at all cell centers
        $:GPU_PARALLEL_LOOP(private='[i,j,k,r]', collapse=3)
        do k = 0, p
            do j = 0, n
                do i = 0, m
                    grad_p_x(i, j, k) = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do r = -fd_number, fd_number
                        grad_p_x(i, j, k) = grad_p_x(i, j, k) + &
                                            q_prim_vf(E_idx)%sf(i + r, j, k)*fd_coeff_x_pgrad(r, i)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! dp/dy at all cell centers
        if (n > 0) then
            $:GPU_PARALLEL_LOOP(private='[i,j,k,r]', collapse=3)
            do k = 0, p
                do j = 0, n
                    do i = 0, m
                        grad_p_y(i, j, k) = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do r = -fd_number, fd_number
                            grad_p_y(i, j, k) = grad_p_y(i, j, k) + &
                                                q_prim_vf(E_idx)%sf(i, j + r, k)*fd_coeff_y_pgrad(r, j)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        ! dp/dz at all cell centers
        if (p > 0) then
            $:GPU_PARALLEL_LOOP(private='[i,j,k,r]', collapse=3)
            do k = 0, p
                do j = 0, n
                    do i = 0, m
                        grad_p_z(i, j, k) = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do r = -fd_number, fd_number
                            grad_p_z(i, j, k) = grad_p_z(i, j, k) + &
                                                q_prim_vf(E_idx)%sf(i, j, k + r)*fd_coeff_z_pgrad(r, k)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_compute_pressure_gradients

    !! This function interpolates the velocity of Eulerian field at the position
            !! of the bubble.
            !! @param pos Position of the bubble in directiion i
            !! @param cell Computational coordinates of the bubble
            !! @param i Direction of the velocity (1: x, 2: y, 3: z)
            !! @param q_prim_vf Eulerian field with primitive variables
            !! @return v Interpolated velocity at the position of the bubble
    function f_interpolate_velocity(pos, cell, i, q_prim_vf) result(v)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: pos
        integer, dimension(3), intent(in) :: cell
        integer, intent(in) :: i
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf

        real(wp) :: v
        real(wp), dimension(fd_order + 1) :: xi, eta, L

        if (fd_order == 2) then
            if (i == 1) then
                xi(1) = x_cc(cell(1) - 1)
                eta(1) = q_prim_vf(momxb)%sf(cell(1) - 1, cell(2), cell(3))
                xi(2) = x_cc(cell(1))
                eta(2) = q_prim_vf(momxb)%sf(cell(1), cell(2), cell(3))
                xi(3) = x_cc(cell(1) + 1)
                eta(3) = q_prim_vf(momxb)%sf(cell(1) + 1, cell(2), cell(3))
            elseif (i == 2) then
                xi(1) = y_cc(cell(2) - 1)
                eta(1) = q_prim_vf(momxb + 1)%sf(cell(1), cell(2) - 1, cell(3))
                xi(2) = y_cc(cell(2))
                eta(2) = q_prim_vf(momxb + 1)%sf(cell(1), cell(2), cell(3))
                xi(3) = y_cc(cell(2) + 1)
                eta(3) = q_prim_vf(momxb + 1)%sf(cell(1), cell(2) + 1, cell(3))
            elseif (i == 3) then
                xi(1) = z_cc(cell(3) - 1)
                eta(1) = q_prim_vf(momxe)%sf(cell(1), cell(2), cell(3) - 1)
                xi(2) = z_cc(cell(3))
                eta(2) = q_prim_vf(momxe)%sf(cell(1), cell(2), cell(3))
                xi(3) = z_cc(cell(3) + 1)
                eta(3) = q_prim_vf(momxe)%sf(cell(1), cell(2), cell(3) + 1)
            end if

            L(1) = ((pos - xi(2))*(pos - xi(3)))/((xi(1) - xi(2))*(xi(1) - xi(3)))
            L(2) = ((pos - xi(1))*(pos - xi(3)))/((xi(2) - xi(1))*(xi(2) - xi(3)))
            L(3) = ((pos - xi(1))*(pos - xi(2)))/((xi(3) - xi(1))*(xi(3) - xi(2)))

            v = L(1)*eta(1) + L(2)*eta(2) + L(3)*eta(3)
        elseif (fd_order == 4) then
            if (i == 1) then
                xi(1) = x_cc(cell(1) - 2)
                eta(1) = q_prim_vf(momxb)%sf(cell(1) - 2, cell(2), cell(3))
                xi(2) = x_cc(cell(1) - 1)
                eta(2) = q_prim_vf(momxb)%sf(cell(1) - 1, cell(2), cell(3))
                xi(3) = x_cc(cell(1))
                eta(3) = q_prim_vf(momxb)%sf(cell(1), cell(2), cell(3))
                xi(4) = x_cc(cell(1) + 1)
                eta(4) = q_prim_vf(momxb)%sf(cell(1) + 1, cell(2), cell(3))
                xi(5) = x_cc(cell(1) + 2)
                eta(5) = q_prim_vf(momxb)%sf(cell(1) + 2, cell(2), cell(3))
            elseif (i == 2) then
                xi(1) = y_cc(cell(2) - 2)
                eta(1) = q_prim_vf(momxb + 1)%sf(cell(1), cell(2) - 2, cell(3))
                xi(2) = y_cc(cell(2) - 1)
                eta(2) = q_prim_vf(momxb + 1)%sf(cell(1), cell(2) - 1, cell(3))
                xi(3) = y_cc(cell(2))
                eta(3) = q_prim_vf(momxb + 1)%sf(cell(1), cell(2), cell(3))
                xi(4) = y_cc(cell(2) + 1)
                eta(4) = q_prim_vf(momxb + 1)%sf(cell(1), cell(2) + 1, cell(3))
                xi(5) = y_cc(cell(2) + 2)
                eta(5) = q_prim_vf(momxb + 1)%sf(cell(1), cell(2) + 2, cell(3))
            elseif (i == 3) then
                xi(1) = z_cc(cell(3) - 2)
                eta(1) = q_prim_vf(momxe)%sf(cell(1), cell(2), cell(3) - 2)
                xi(2) = z_cc(cell(3) - 1)
                eta(2) = q_prim_vf(momxe)%sf(cell(1), cell(2), cell(3) - 1)
                xi(3) = z_cc(cell(3))
                eta(3) = q_prim_vf(momxe)%sf(cell(1), cell(2), cell(3))
                xi(4) = z_cc(cell(3) + 1)
                eta(4) = q_prim_vf(momxe)%sf(cell(1), cell(2), cell(3) + 1)
                xi(5) = z_cc(cell(3) + 2)
                eta(5) = q_prim_vf(momxe)%sf(cell(1), cell(2), cell(3) + 2)
            end if

            L(1) = ((pos - xi(2))*(pos - xi(3))*(pos - xi(4))*(pos - xi(5)))/ &
                   ((xi(1) - xi(2))*(xi(1) - xi(3))*(xi(1) - xi(3))*(xi(2) - xi(5)))
            L(2) = ((pos - xi(1))*(pos - xi(3))*(pos - xi(4))*(pos - xi(5)))/ &
                   ((xi(2) - xi(1))*(xi(2) - xi(3))*(xi(2) - xi(3))*(xi(2) - xi(5)))
            L(3) = ((pos - xi(1))*(pos - xi(2))*(pos - xi(4))*(pos - xi(5)))/ &
                   ((xi(3) - xi(1))*(xi(3) - xi(2))*(xi(3) - xi(4))*(xi(3) - xi(5)))
            L(4) = ((pos - xi(1))*(pos - xi(2))*(pos - xi(3))*(pos - xi(4)))/ &
                   ((xi(4) - xi(1))*(xi(4) - xi(2))*(xi(4) - xi(3))*(xi(4) - xi(5)))
            L(5) = ((pos - xi(1))*(pos - xi(2))*(pos - xi(3))*(pos - xi(4)))/ &
                   ((xi(5) - xi(1))*(xi(5) - xi(2))*(xi(5) - xi(3))*(xi(5) - xi(4)))

            v = L(1)*eta(1) + L(2)*eta(2) + L(3)*eta(3) + L(4)*eta(4) + L(5)*eta(5)
        end if

    end function f_interpolate_velocity

    !! This function calculates the force on a bubble
            !!      based on the pressure gradient, velocity, and drag model.
            !! @param pos Position of the bubble in direction i
            !! @param rad Radius of the bubble
            !! @param rdot Radial velocity of the bubble
            !! @param vel Velocity of the bubble
            !! @param mg Mass of the gas in the bubble
            !! @param mv Mass of the liquid in the bubble
            !! @param Re Reynolds number
            !! @param rho Density of the fluid
            !! @param cell Computational coordinates of the bubble
            !! @param i Direction of the velocity (1: x, 2: y, 3: z)
            !! @param q_prim_vf Eulerian field with primitive variables
            !! @return a Acceleration of the bubble in direction i
    function f_get_bubble_force(pos, rad, rdot, vel, mg, mv, Re, rho, cell, i, q_prim_vf) result(force)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: pos, rad, rdot, mg, mv, Re, rho, vel
        integer, dimension(3), intent(in) :: cell
        integer, intent(in) :: i
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf

        real(wp) :: dp, vol, force
        real(wp) :: v_rel

        if (fd_order > 1) then
            v_rel = vel - f_interpolate_velocity(pos, cell, i, q_prim_vf)
        else
            v_rel = vel - q_prim_vf(momxb + i - 1)%sf(cell(1), cell(2), cell(3))
        end if

        force = 0._wp

        if (lag_params%drag_model == 1) then ! Free slip Stokes drag
            force = force - (4._wp*pi*rad*v_rel)/Re
        else if (lag_params%drag_model == 2) then ! No slip Stokes drag
            force = force - (6._wp*pi*rad*v_rel)/Re
        else if (lag_params%drag_model == 3) then ! Levich drag
            force = force - (12._wp*pi*rad*v_rel)/Re
        end if

        if (lag_pressure_force) then
            ! Use precomputed cell-centered pressure gradients
            if (i == 1) then
                dp = grad_p_x(cell(1), cell(2), cell(3))
            elseif (i == 2) then
                dp = grad_p_y(cell(1), cell(2), cell(3))
            elseif (i == 3) then
                dp = grad_p_z(cell(1), cell(2), cell(3))
            end if

            vol = (4._wp/3._wp)*pi*(rad**3._wp)
            force = force - vol*dp
        end if

        if (lag_params%gravity_force) then
            force = force + (mg + mv)*accel_bf(i)
        end if

    end function f_get_bubble_force

end module m_bubbles_EL_kernels
