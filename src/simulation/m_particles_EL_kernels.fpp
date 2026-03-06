!>
!! @file m_particles_EL_kernels.f90
!! @brief Contains module m_particles_EL_kernels

#:include 'macros.fpp'

!> @brief This module contains kernel functions used to map the effect of the lagrangian particles
!!        in the Eulerian framework.
module m_particles_EL_kernels

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use ieee_arithmetic        !< For checking NaN

    implicit none

    ! Cell list for particle-to-cell mapping (rebuilt each RK stage before smearing)
    integer, allocatable, dimension(:, :, :) :: cell_list_start  ! (0:m, 0:n, 0:p)
    integer, allocatable, dimension(:, :, :) :: cell_list_count  ! (0:m, 0:n, 0:p)
    integer, allocatable, dimension(:) :: cell_list_idx          ! (1:nParticles_glb) sorted particle indices
    $:GPU_DECLARE(create='[cell_list_start, cell_list_count, cell_list_idx]')

contains

    !> The purpose of this subroutine is to smear the strength of the lagrangian
            !!      particles into the Eulerian framework using a gaussian approach.
            !! @param nParticles Number of lagrangian particles in the current domain
            !! @param lbk_rad Radius of the particles
            !! @param lbk_s Computational coordinates of the particles
            !! @param lbk_pos Spatial coordinates of the particles
            !! @param updatedvar Eulerian variable to be updated
            !! @param lbk_f_p Forces on the particles
    subroutine s_smoothfunction(nParticles, lbk_rad, lbk_vel, lbk_s, lbk_pos, lbk_f_p, lbk_g_sum, updatedvar, kcomp, lbk_linked_list, lbk_particle_head)

        integer, intent(in) :: nParticles
        real(wp), dimension(1:lag_params%nParticles_glb, 1:3, 1:2), intent(in) :: lbk_s, lbk_pos, lbk_vel
        real(wp), dimension(1:lag_params%nParticles_glb, 1:2), intent(in) :: lbk_rad
        real(wp), dimension(1:lag_params%nParticles_glb, 1:3), intent(in) :: lbk_f_p
        real(wp), dimension(1:lag_params%nParticles_glb), intent(inout) :: lbk_g_sum

        integer, dimension(1:lag_params%nParticles_glb), intent(in) :: lbk_linked_list
        integer, dimension(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end), intent(in) :: lbk_particle_head

        type(scalar_field), dimension(:), intent(inout) :: updatedvar
        type(scalar_field), dimension(:), intent(inout) :: kcomp

        ! call s_gaussian(nParticles, lbk_rad, lbk_s, lbk_pos, lbk_vel, updatedvar, lbk_f_p, ngh)
        call s_gaussian(nParticles, lbk_rad, lbk_vel, lbk_s, lbk_pos, lbk_f_p, lbk_g_sum, updatedvar, kcomp, lbk_linked_list, lbk_particle_head)

    end subroutine s_smoothfunction

    !> Cell-centric gaussian smearing using the cell list (no GPU atomics).
    !!      Each grid cell accumulates contributions from nearby particles looked up
    !!      via cell_list_start/count/idx.
    ! nParticles, lbk_rad, lbk_s, lbk_pos, lbk_vel, updatedvar, lbk_f_p)
    subroutine s_gaussian(nParticles, lbk_rad, lbk_vel, lbk_s, lbk_pos, lbk_f_p, lbk_g_sum, updatedvar, kcomp, lbk_linked_list, lbk_particle_head)

        integer, intent(in) :: nParticles
        real(wp), dimension(1:lag_params%nParticles_glb, 1:3, 1:2), intent(in) :: lbk_s, lbk_pos, lbk_vel
        real(wp), dimension(1:lag_params%nParticles_glb, 1:2), intent(in) :: lbk_rad
        real(wp), dimension(1:lag_params%nParticles_glb, 1:3), intent(in) :: lbk_f_p
        real(wp), dimension(1:lag_params%nParticles_glb), intent(inout) :: lbk_g_sum

        integer, dimension(1:lag_params%nParticles_glb), intent(in) :: lbk_linked_list
        integer, dimension(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end), intent(in) :: lbk_particle_head

        type(scalar_field), dimension(:), intent(inout) :: updatedvar
        type(scalar_field), dimension(:), intent(inout) :: kcomp

        real(wp), dimension(3) :: center, nodecoord, s_coord
        integer, dimension(3) :: cell, cellijk
        real(wp) :: stddsv, volpart, Vol_loc
        real(wp) :: strength_vel, strength_vol
        real(wp) :: func
        real(wp) :: y_kahan, t_kahan
        real(wp) :: fp_x, fp_y, fp_z, vp_x, vp_y, vp_z
        integer :: l, i, j, k, di, dj, dk, lb, part_idx
        integer :: di_beg, di_end, dj_beg, dj_end, dk_beg, dk_end
        integer :: smear_x_beg, smear_x_end
        integer :: smear_y_beg, smear_y_end
        integer :: smear_z_beg, smear_z_end
        integer :: mapCells_loc

        mapCells_loc = 1

        smear_x_beg = -mapCells_loc
        smear_x_end = m + mapCells_loc
        smear_y_beg = merge(-mapCells_loc, 0, n > 0)
        smear_y_end = merge(n + mapCells_loc, n, n > 0)
        smear_z_beg = merge(-mapCells_loc, 0, p > 0)
        smear_z_end = merge(p + mapCells_loc, p, p > 0)

        $:GPU_PARALLEL_LOOP(private='[l,volpart,s_coord,cell,stddsv,i,j,k,di,dj,dk,di_beg,di_end,dj_beg,dj_end,dk_beg,dk_end,nodecoord,center,cellijk,Vol_loc,func]')
        do l = 1, nParticles

            lbk_g_sum(l) = 0._wp

            volpart = 4._wp/3._wp*pi*lbk_rad(l, 2)**3._wp
            s_coord(1:3) = lbk_s(l, 1:3, 2)

            call s_get_cell(s_coord, cell)

            call s_compute_stddsv(cell, volpart, stddsv)

            i = cell(1)
            j = cell(2)
            k = cell(3)

            di_beg = i - mapCells_loc
            di_end = i + mapCells_loc
            dj_beg = j - mapCells_loc
            dj_end = j + mapCells_loc
            dk_beg = k - mapCells_loc
            dk_end = k + mapCells_loc

            $:GPU_LOOP(parallelism='[seq]')
            do di = di_beg, di_end
                $:GPU_LOOP(parallelism='[seq]')
                do dj = dj_beg, dj_end
                    $:GPU_LOOP(parallelism='[seq]')
                    do dk = dk_beg, dk_end

                        nodecoord(1) = x_cc(di)
                        nodecoord(2) = y_cc(dj)
                        nodecoord(3) = 0._wp
                        if (p > 0) nodecoord(3) = z_cc(dk)

                        cellijk(1) = di
                        cellijk(2) = dj
                        cellijk(3) = dk

                        center(1:2) = lbk_pos(l, 1:2, 2)
                        center(3) = 0._wp
                        if (p > 0) center(3) = lbk_pos(l, 3, 2)

                        Vol_loc = dx(cellijk(1))*dy(cellijk(2))
                        if (num_dims == 3) Vol_loc = Vol_loc*dz(cellijk(3))

                        call s_applygaussian(center, cellijk, nodecoord, stddsv, 0._wp, func)

                        lbk_g_sum(l) = lbk_g_sum(l) + func*Vol_loc
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        $:GPU_PARALLEL_LOOP(collapse=3, &
            & private='[i,j,k,di,dj,dk,lb,part_idx,center,nodecoord,s_coord,cell,cellijk,stddsv,volpart,strength_vol,fp_x,fp_y,fp_z,vp_x,vp_y,vp_z,func,y_kahan,t_kahan,di_beg,di_end,dj_beg,dj_end,dk_beg,dk_end]', &
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
                    di_beg = max(i - mapCells_loc, 0)
                    di_end = min(i + mapCells_loc, m)
                    dj_beg = max(j - mapCells_loc, 0)
                    dj_end = min(j + mapCells_loc, n)
                    dk_beg = max(k - mapCells_loc, 0)
                    dk_end = min(k + mapCells_loc, p)

                    $:GPU_LOOP(parallelism='[seq]')
                    do dk = dk_beg, dk_end
                        $:GPU_LOOP(parallelism='[seq]')
                        do dj = dj_beg, dj_end
                            $:GPU_LOOP(parallelism='[seq]')
                            do di = di_beg, di_end

                                part_idx = lbk_particle_head(di, dj, dk)

                                do while (part_idx /= -1)

                                    ! do lb = cell_list_start(di, dj, dk), &
                                    ! cell_list_start(di, dj, dk) + cell_list_count(di, dj, dk) - 1

                                    ! Particle properties
                                    volpart = 4._wp/3._wp*pi*lbk_rad(part_idx, 2)**3._wp
                                    s_coord(1:3) = lbk_s(part_idx, 1:3, 2)
                                    call s_get_cell(s_coord, cell)

                                    fp_x = -lbk_f_p(part_idx, 1)
                                    fp_y = -lbk_f_p(part_idx, 2)
                                    fp_z = -lbk_f_p(part_idx, 3)

                                    vp_x = lbk_vel(part_idx, 1, 2)
                                    vp_y = lbk_vel(part_idx, 2, 2)
                                    vp_z = lbk_vel(part_idx, 3, 2)

                                    call s_compute_stddsv(cell, volpart, stddsv)

                                    strength_vol = volpart
                                    ! strength_vel = 4._wp*pi*lbk_rad(part_idx, 2)**2._wp*lbk_vel(part_idx, 2)

                                    center(1:2) = lbk_pos(part_idx, 1:2, 2)
                                    center(3) = 0._wp
                                    if (p > 0) center(3) = lbk_pos(part_idx, 3, 2)

                                    call s_applygaussian(center, cellijk, nodecoord, stddsv, 0._wp, func)

                                    func = func/lbk_g_sum(part_idx)

                                    ! Kahan summation for void fraction
                                    y_kahan = real(func*strength_vol, kind=wp) - kcomp(1)%sf(i, j, k)
                                    t_kahan = updatedvar(1)%sf(i, j, k) + y_kahan
                                    kcomp(1)%sf(i, j, k) = (t_kahan - updatedvar(1)%sf(i, j, k)) - y_kahan
                                    updatedvar(1)%sf(i, j, k) = t_kahan

                                    if (lag_params%solver_approach == 2) then

                                        ! Kahan summation for alpha_p vp_x
                                        y_kahan = real(func*strength_vol*vp_x, kind=wp) - kcomp(2)%sf(i, j, k)
                                        t_kahan = updatedvar(2)%sf(i, j, k) + y_kahan
                                        kcomp(2)%sf(i, j, k) = (t_kahan - updatedvar(2)%sf(i, j, k)) - y_kahan
                                        updatedvar(2)%sf(i, j, k) = t_kahan

                                        ! Kahan summation for alpha_p vp_y
                                        y_kahan = real(func*strength_vol*vp_y, kind=wp) - kcomp(3)%sf(i, j, k)
                                        t_kahan = updatedvar(3)%sf(i, j, k) + y_kahan
                                        kcomp(3)%sf(i, j, k) = (t_kahan - updatedvar(3)%sf(i, j, k)) - y_kahan
                                        updatedvar(3)%sf(i, j, k) = t_kahan

                                        ! Kahan summation for alpha_p vp_z
                                        y_kahan = real(func*strength_vol*vp_z, kind=wp) - kcomp(4)%sf(i, j, k)
                                        t_kahan = updatedvar(4)%sf(i, j, k) + y_kahan
                                        kcomp(4)%sf(i, j, k) = (t_kahan - updatedvar(4)%sf(i, j, k)) - y_kahan
                                        updatedvar(4)%sf(i, j, k) = t_kahan

                                        ! Kahan summation for fp_x
                                        y_kahan = real(func*fp_x, kind=wp) - kcomp(5)%sf(i, j, k)
                                        t_kahan = updatedvar(5)%sf(i, j, k) + y_kahan
                                        kcomp(5)%sf(i, j, k) = (t_kahan - updatedvar(5)%sf(i, j, k)) - y_kahan
                                        updatedvar(5)%sf(i, j, k) = t_kahan

                                        ! Kahan summation for fp_y
                                        y_kahan = real(func*fp_y, kind=wp) - kcomp(6)%sf(i, j, k)
                                        t_kahan = updatedvar(6)%sf(i, j, k) + y_kahan
                                        kcomp(6)%sf(i, j, k) = (t_kahan - updatedvar(6)%sf(i, j, k)) - y_kahan
                                        updatedvar(6)%sf(i, j, k) = t_kahan

                                        ! Kahan summation for fp_z
                                        y_kahan = real(func*fp_z, kind=wp) - kcomp(7)%sf(i, j, k)
                                        t_kahan = updatedvar(7)%sf(i, j, k) + y_kahan
                                        kcomp(7)%sf(i, j, k) = (t_kahan - updatedvar(7)%sf(i, j, k)) - y_kahan
                                        updatedvar(7)%sf(i, j, k) = t_kahan

                                        ! Kahan summation for E_source
                                        y_kahan = real(func*0._wp, kind=wp) - kcomp(8)%sf(i, j, k)
                                        t_kahan = updatedvar(8)%sf(i, j, k) + y_kahan
                                        kcomp(8)%sf(i, j, k) = (t_kahan - updatedvar(8)%sf(i, j, k)) - y_kahan
                                        updatedvar(8)%sf(i, j, k) = t_kahan

                                    end if

                                    part_idx = lbk_linked_list(part_idx)
                                end do
                            end do
                        end do
                    end do

                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_gaussian

    !> The purpose of this subroutine is to apply the gaussian kernel function for each particle (Maeda and Colonius, 2018)).
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

    !> The purpose of this subroutine is to check if the current cell is outside the computational domain or not (including ghost cells).
            !! @param cellaux Tested cell to smear the particle effect in.
            !! @param celloutside If true, then cellaux is outside the computational domain.
    subroutine s_check_celloutside(cellaux, celloutside)
        $:GPU_ROUTINE(function_name='s_check_celloutside',parallelism='[seq]', &
            & cray_inline=True)

        integer, dimension(3), intent(inout) :: cellaux
        logical, intent(out) :: celloutside

        celloutside = .false.

        if (num_dims == 2) then
            if ((cellaux(1) < -buff_size) .or. (cellaux(2) < -buff_size)) then
                celloutside = .true.
            end if

            if ((cellaux(1) > m + buff_size) .or. (cellaux(2) > n + buff_size)) then
                celloutside = .true.
            end if
        else
            if ((cellaux(1) < -buff_size) .or. (cellaux(2) < -buff_size) .or. (cellaux(3) < -buff_size)) then
                celloutside = .true.
            end if

            if ((cellaux(1) > m + buff_size) .or. (cellaux(2) > n + buff_size) .or. (cellaux(3) > p + buff_size)) then
                celloutside = .true.
            end if
        end if

    end subroutine s_check_celloutside

    !> This subroutine transforms the computational coordinates of the particle from
            !!      real type into integer.
            !! @param s Computational coordinates of the particle, real type
            !! @param get_cell Computational coordinates of the particle, integer type
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

    !> The purpose of this procedure is to calculate the characteristic cell volume
            !! @param cell Computational coordinates (x, y, z)
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

    !! This function calculates the force on a particle
            !!      based on the pressure gradient, velocity, and drag model.
            !! @param pos Position of the particle
            !! @param rad Radius of the particle
            !! @param vel_p Velocity of the particle
            !! @param mass_p Particle mass
            !! @param Re Viscosity!
            !! @param rho Density of the fluid
            !! @param vol_frac Particle Volume Fraction
            !! @param cell Computational coordinates of the particle
            !! @param q_prim_vf Eulerian field with primitive variables
            !! @return a Acceleration of the particle in direction i
    subroutine s_get_particle_force(pos, rad, vel_p, mass_p, Re, gamm, vol_frac, cell, &
                                    q_prim_vf, fieldvars, wx, wy, wz, force, rmass_add)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: rad, mass_p, Re, gamm, vol_frac
        real(wp), dimension(3), intent(in) :: pos
        integer, dimension(3), intent(in) :: cell
        real(wp), dimension(3), intent(in) :: vel_p
        type(scalar_field), dimension(:), intent(in) :: fieldvars
        type(scalar_field), dimension(:), intent(in) :: wx, wy, wz
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf

        real(wp), dimension(3), intent(out) :: force
        real(wp), intent(out) :: rmass_add

        real(wp) :: a, vol, rho_fluid, pressure_fluid
        real(wp), dimension(3) :: v_rel, dp
        real(wp), dimension(fd_order) :: xi, eta, L
        real(wp) :: particle_diam, gas_mu, vmag, cson
        real(wp) :: slip_velocity_x, slip_velocity_y, slip_velocity_z, beta
        real(wp), dimension(3) :: fluid_vel
        integer :: dir

        !Added pass params
        real(wp) :: mach, Cam, udot_grad_rho, flux_f, flux_b, div_rhou, SDrho, vgradrho, drhodt
        real(wp), dimension(3) :: rhoDuDt, grad_rho, fam
        integer, dimension(3) :: p1

        force = 0._wp
        dp = 0._wp
        grad_rho = 0._wp
        drhodt = 0._wp
        fam = 0._wp
        fluid_vel = 0._wp
        v_rel = 0._wp
        rhoDuDt = 0._wp

        !!Interpolation - either even ordered barycentric or 0th order
        if (lag_params%interpolation_order > 1) then
            rho_fluid = f_interp_barycentric(pos, cell, q_prim_vf, 1, wx, wy, wz)
            pressure_fluid = f_interp_barycentric(pos, cell, q_prim_vf, E_idx, wx, wy, wz)
            do dir = 1, num_dims
                if (lag_params%pressure_force .or. lag_params%added_mass_model > 0) then
                    dp(dir) = f_interp_barycentric(pos, cell, fieldvars, dir, wx, wy, wz)
                end if
                if (lag_params%added_mass_model > 0) then
                    grad_rho(dir) = f_interp_barycentric(pos, cell, fieldvars, 3 + dir, wx, wy, wz)
                    drhodt = drhodt + f_interp_barycentric(pos, cell, fieldvars, 6 + dir, wx, wy, wz)
                end if
                fluid_vel(dir) = f_interp_barycentric(pos, cell, q_prim_vf, momxb + dir - 1, wx, wy, wz)
            end do
        else
            rho_fluid = q_prim_vf(1)%sf(cell(1), cell(2), cell(3))
            pressure_fluid = q_prim_vf(E_idx)%sf(cell(1), cell(2), cell(3))
            do dir = 1, num_dims
                if (lag_params%pressure_force .or. lag_params%added_mass_model > 0) then
                    dp(dir) = fieldvars(dir)%sf(cell(1), cell(2), cell(3))
                end if
                if (lag_params%added_mass_model > 0) then
                    grad_rho(dir) = fieldvars(3 + dir)%sf(cell(1), cell(2), cell(3))
                    drhodt = drhodt + fieldvars(6 + dir)%sf(cell(1), cell(2), cell(3))
                end if
                fluid_vel(dir) = q_prim_vf(momxb + dir - 1)%sf(cell(1), cell(2), cell(3))
            end do
        end if

        drhodt = -drhodt

        v_rel = vel_p - fluid_vel

        if (lag_params%qs_drag_model > 0 .or. lag_params%added_mass_model > 0) then
            ! Quasi-steady Drag Force Parameters
            slip_velocity_x = fluid_vel(1) - vel_p(1)
            slip_velocity_y = fluid_vel(2) - vel_p(2)
            if (num_dims == 3) then
                slip_velocity_z = fluid_vel(3) - vel_p(3)
                vmag = sqrt(slip_velocity_x*slip_velocity_x + slip_velocity_y*slip_velocity_y + &
                            slip_velocity_z*slip_velocity_z)
            elseif (num_dims == 2) then
                vmag = sqrt(slip_velocity_x*slip_velocity_x + slip_velocity_y*slip_velocity_y)
            end if
            particle_diam = rad*2._wp
            cson = sqrt((gamm*pressure_fluid)/rho_fluid) !gamma*P/rho
            gas_mu = Re
        end if

        if (lag_params%added_mass_model > 0) then
            rhoDuDt = -dp
            udot_grad_rho = dot_product(fluid_vel, grad_rho)
            vgradrho = dot_product(vel_p, grad_rho)
            SDrho = (drhodt + udot_grad_rho) !/(1._wp - vol_frac)
            mach = vmag/cson
        end if

        ! Step 1: Force component quasi-steady
        if (lag_params%qs_drag_model == 1) then
            beta = QS_Parmar(rho_fluid, cson, gas_mu, gamm, vmag, particle_diam, vol_frac)
            force = force - beta*v_rel
        else if (lag_params%qs_drag_model == 2) then
            beta = QS_Osnes(rho_fluid, cson, gas_mu, gamm, vmag, particle_diam, vol_frac)
            force = force - beta*v_rel
        else if (lag_params%qs_drag_model == 3) then
            beta = QS_ModifiedParmar(rho_fluid, cson, gas_mu, gamm, vmag, particle_diam, vol_frac)
            force = force - beta*v_rel
        else if (lag_params%qs_drag_model == 4) then
            beta = QS_Gidaspow(rho_fluid, cson, gas_mu, gamm, vmag, particle_diam, vol_frac)
            force = force - beta*v_rel
        else
            !No Quasi-Steady drag
        end if

        ! Step 1.1: Stokes drag
        if (lag_params%stokes_drag == 1) then ! Free slip Stokes drag
            force = force - 4._wp*pi*gas_mu*rad*v_rel
        elseif (lag_params%stokes_drag == 2) then ! No slip Stokes drag
            force = force - 6._wp*pi*gas_mu*rad*v_rel
        else
            !No stokes drag
        end if

        ! Step 2: Pressure Gradient Force
        if (lag_params%pressure_force) then
            vol = (4._wp/3._wp)*pi*(rad**3._wp)
            force = force - vol*dp
        end if

        ! Step 3: Gravitational Force
        if (lag_params%gravity_force) then
            force = force + (mass_p)*accel_bf
        end if

        ! Step 4: Added Mass Force
        if (lag_params%added_mass_model == 1) then
            vol = (4._wp/3._wp)*pi*(rad**3._wp)
            if (mach > 0.6_wp) then
                Cam = 1._wp + 1.8_wp*(0.6_wp**2) + 7.6_wp* &
                      (0.6_wp**4)
            else
                Cam = 1._wp + 1.8_wp*mach**2 + 7.6_wp*mach**4
            end if

            Cam = 0.5_wp*Cam*(1._wp + 0.68_wp*vol_frac**2)
            rmass_add = rho_fluid*vol*Cam !(1._wp-vol_frac)*rho_fluid*vol*Cam

            fam = Cam*vol*(vel_p*SDrho + rhoDuDt + &
                           fluid_vel*(vgradrho))

            do dir = 1, num_dims
                if (.not. ieee_is_finite(fam(dir))) then
                    fam(dir) = 0._wp
                    rmass_add = 0._wp
                end if
            end do
            force = force + fam
        else
            rmass_add = 0._wp
        end if

        ! do dir = 1, num_dims
        !     if (.not. ieee_is_finite(force(dir))) then
        !         force(dir) = 0._wp
        !     end if
        ! end do

    end subroutine s_get_particle_force

    function f_interp_barycentric(pos, cell, field_vf, field_index, wx, wy, wz) result(val)
        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), dimension(3), intent(in) :: pos
        integer, dimension(3), intent(in) :: cell
        type(scalar_field), dimension(:), intent(in) :: field_vf
        type(scalar_field), dimension(:), intent(in) :: wx, wy, wz
        integer, intent(in) :: field_index

        integer :: i, j, k, ix, jy, kz, npts, npts_z, N, a, b
        integer :: ix_count, jy_count, kz_count
        real(wp) :: weight, numerator, denominator, xBar, eps
        real(wp) :: val, local_min, local_max, prod_x, prod_y, prod_z

        i = cell(1)
        j = cell(2)
        k = cell(3)

        N = lag_params%interpolation_order
        npts = N/2
        npts_z = npts
        if (num_dims == 2) npts_z = 0
        eps = 1.e-12_wp
        numerator = 0._wp
        denominator = 0._wp

        ! if (abs(pos(1) - x_cc(i)) <= eps .and. &
        !     abs(pos(2) - y_cc(j)) <= eps .and. &
        !     abs(pos(3) - z_cc(k)) <= eps) then
        !     val = field_vf(field_index)%sf(i, j, k)
        !     return
        ! end if

        ix_count = 0
        do ix = i - npts, i + npts
            ix_count = ix_count + 1
            jy_count = 0
            do jy = j - npts, j + npts
                jy_count = jy_count + 1
                kz_count = 0
                do kz = k - npts_z, k + npts_z
                    kz_count = kz_count + 1
                    if (num_dims == 3) then
                        xBar = (pos(1) - x_cc(ix))*(pos(2) - y_cc(jy))*(pos(3) - z_cc(kz))
                        weight = wx(ix_count)%sf(i, 1, 1)*wy(jy_count)%sf(j, 1, 1)*wz(kz_count)%sf(k, 1, 1)
                    else
                        xBar = (pos(1) - x_cc(ix))*(pos(2) - y_cc(jy))
                        weight = wx(ix_count)%sf(i, 1, 1)*wy(jy_count)%sf(j, 1, 1)
                    end if
                    weight = weight/xBar
                    numerator = numerator + weight*field_vf(field_index)%sf(ix, jy, kz)
                    denominator = denominator + weight
                end do
            end do
        end do

        val = numerator/denominator

        if (.not. ieee_is_finite(val)) then
            val = field_vf(field_index)%sf(i, j, k)

        elseif (abs(val) <= eps) then
            val = 0._wp

        end if

    end function

    ! Quasi-steady force (Re_p and Ma_p corrections):
    !   Improved Drag Correlation for Spheres and Application
    !   to Shock-Tube Experiments
    !   - Parmar et al. (2010)
    !   - AIAA Journal
    !
    ! Quasi-steady force (phi corrections):
    !   The Added Mass, Basset, and Viscous Drag Coefficients
    !   in Nondilute Bubbly Liquids Undergoing Small-Amplitude
    !   Oscillatory Motion
    !   - Sangani et al. (1991)
    !   - Phys. Fluids A
    function QS_Parmar(rho, cson, mu_fluid, gamma, vmag, dp, volume_fraction) result(beta)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: rho, cson, mu_fluid, gamma, vmag, dp, volume_fraction
        real(wp) :: rcd1, rmacr, rcd_mcr, rcd_std, rmach_rat, rcd_M1
        real(wp) :: rcd_M2, C1, C2, C3, f1M, f2M, f3M, lrep, factor, cd, phi_corr
        real(wp) :: beta
        real(wp) :: rmachp, mp, phi, rep, re

        rmachp = vmag/cson
        mp = max(rmachp, 0.01_wp)
        phi = max(volume_fraction, 0.0001_wp)
        rep = vmag*dp*rho/mu_fluid
        re = max(rep, 0.1_wp)

        if (re < 1.e-14_wp) then
            rcd1 = 1.0_wp
        else
            rmacr = 0.6_wp ! Critical rmachp no
            rcd_mcr = (1.+0.15*re**(0.684)) + &
                      (re/24.0)*(0.513/(1.+483./re**(0.669)))
            if (mp <= rmacr) then
                rcd_std = (1.+0.15*re**(0.687)) + &
                          (re/24.0)*(0.42/(1.+42500./re**(1.16)))
                rmach_rat = mp/rmacr
                rcd1 = rcd_std + (rcd_mcr - rcd_std)*rmach_rat
            else if (mp <= 1.0) then
                rcd_M1 = (1.0 + 0.118*re**0.813) + &
                         (re/24.0)*0.69/(1.0 + 3550.0/re**.793)
                C1 = 6.48_wp
                C2 = 9.28_wp
                C3 = 12.21_wp
                f1M = -1.884_wp + 8.422_wp*mp - 13.70_wp*mp**2 + 8.162_wp*mp**3
                f2M = -2.228_wp + 10.35_wp*mp - 16.96_wp*mp**2 + 9.840_wp*mp**3
                f3M = 4.362_wp - 16.91_wp*mp + 19.84_wp*mp**2 - 6.296_wp*mp**3
                lrep = log(re)
                factor = f1M*(lrep - C2)*(lrep - C3)/((C1 - C2)*(C1 - C3)) &
                         + f2M*(lrep - C1)*(lrep - C3)/((C2 - C1)*(C2 - C3)) &
                         + f3M*(lrep - C1)*(lrep - C2)/((C3 - C1)*(C3 - C2))
                rcd1 = rcd_mcr + (rcd_M1 - rcd_mcr)*factor
            else if (mp < 1.75) then
                rcd_M1 = (1.0 + 0.118*re**0.813) + &
                         (re/24.0)*0.69/(1.0 + 3550.0/re**.793)
                rcd_M2 = (1.0 + 0.107*re**0.867) + &
                         (re/24.0)*0.646/(1.0 + 861.0/re**.634)
                C1 = 6.48_wp
                C2 = 8.93_wp
                C3 = 12.21_wp
                f1M = -2.963 + 4.392*mp - 1.169*mp**2 - 0.027*mp**3 &
                      - 0.233*exp((1.0 - mp)/0.011)
                f2M = -6.617 + 12.11*mp - 6.501*mp**2 + 1.182*mp**3 &
                      - 0.174*exp((1.0 - mp)/0.010)
                f3M = -5.866 + 11.57*mp - 6.665*mp**2 + 1.312*mp**3 &
                      - 0.350*exp((1.0 - mp)/0.012)
                lrep = log(re)
                factor = f1M*(lrep - C2)*(lrep - C3)/((C1 - C2)*(C1 - C3)) &
                         + f2M*(lrep - C1)*(lrep - C3)/((C2 - C1)*(C2 - C3)) &
                         + f3M*(lrep - C1)*(lrep - C2)/((C3 - C1)*(C3 - C2))
                rcd1 = rcd_M1 + (rcd_M2 - rcd_M1)*factor
            else
                rcd1 = (1.0 + 0.107*re**0.867) + &
                       (re/24.0)*0.646/(1.0 + 861.0/re**.634)
            end if ! mp
        end if    ! re

        ! Sangani's volume fraction correction for dilute random arrays
        ! Capping volume fraction at 0.5
        phi_corr = (1.0 + 5.94*min(phi, 0.5))

        cd = (24.0/re)*rcd1*phi_corr

        beta = rcd1*3.0*pi*mu_fluid*dp

        beta = beta*phi_corr

    end function QS_Parmar

    ! Quasi-steady force (Re_p and Ma_p corrections):
    !   Improved Drag Correlation for Spheres and Application
    !   to Shock-Tube Experiments
    !   - Parmar et al. (2010)
    !   - AIAA Journal
    !
    ! Quasi-steady force (phi corrections):
    !   Sangani et al. (1991) volume fraction correction overshoots
    !   the drag coefficient.
    !
    !   We adopt instead Osnes et al. (2023) volume fraction correction
    !   based on Tenneti et al. with one extra term.
    !
    !   At Mach=0, the drag coefficient from this subroutine matches very
    !   well with the one calculated using the Osnes subroutine, for various
    !   Reynolds numbers and volume fractions.
    function QS_ModifiedParmar(rho, cson, mu_fluid, gamma, vmag, dp, volume_fraction) result(beta)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: rho, cson, mu_fluid, gamma, vmag, dp, volume_fraction
        real(wp) :: rcd1, rmacr, rcd_mcr, rcd_std, rmach_rat, rcd_M1
        real(wp) :: rcd_M2, C1, C2, C3, f1M, f2M, f3M, lrep, factor, cd, phi_corr
        real(wp) :: b1, b2, b3
        real(wp) :: beta
        real(wp) :: rmachp, mp, phi, rep, re

        rmachp = vmag/cson
        mp = max(rmachp, 0.01_wp)
        phi = max(volume_fraction, 0.0001_wp)
        rep = vmag*dp*rho/mu_fluid
        re = max(rep, 0.1_wp)

        if (re < 1e-14) then
            rcd1 = 1.0
        else
            rmacr = 0.6 ! Critical rmachp no.
            rcd_mcr = (1.+0.15*re**(0.684)) + &
                      (re/24.0)*(0.513/(1.+483./re**(0.669)))
            if (mp <= rmacr) then
                rcd_std = (1.+0.15*re**(0.687)) + &
                          (re/24.0)*(0.42/(1.+42500./re**(1.16)))
                rmach_rat = mp/rmacr
                rcd1 = rcd_std + (rcd_mcr - rcd_std)*rmach_rat
            else if (mp <= 1.0) then
                rcd_M1 = (1.0 + 0.118*re**0.813) + &
                         (re/24.0)*0.69/(1.0 + 3550.0/re**.793)
                C1 = 6.48
                C2 = 9.28
                C3 = 12.21
                f1M = -1.884 + 8.422*mp - 13.70*mp**2 + 8.162*mp**3
                f2M = -2.228 + 10.35*mp - 16.96*mp**2 + 9.840*mp**3
                f3M = 4.362 - 16.91*mp + 19.84*mp**2 - 6.296*mp**3
                lrep = log(re)
                factor = f1M*(lrep - C2)*(lrep - C3)/((C1 - C2)*(C1 - C3)) &
                         + f2M*(lrep - C1)*(lrep - C3)/((C2 - C1)*(C2 - C3)) &
                         + f3M*(lrep - C1)*(lrep - C2)/((C3 - C1)*(C3 - C2))
                rcd1 = rcd_mcr + (rcd_M1 - rcd_mcr)*factor
            else if (mp < 1.75) then
                rcd_M1 = (1.0 + 0.118*re**0.813) + &
                         (re/24.0)*0.69/(1.0 + 3550.0/re**.793)
                rcd_M2 = (1.0 + 0.107*re**0.867) + &
                         (re/24.0)*0.646/(1.0 + 861.0/re**.634)
                C1 = 6.48
                C2 = 8.93
                C3 = 12.21
                f1M = -2.963 + 4.392*mp - 1.169*mp**2 - 0.027*mp**3 &
                      - 0.233*exp((1.0 - mp)/0.011)
                f2M = -6.617 + 12.11*mp - 6.501*mp**2 + 1.182*mp**3 &
                      - 0.174*exp((1.0 - mp)/0.010)
                f3M = -5.866 + 11.57*mp - 6.665*mp**2 + 1.312*mp**3 &
                      - 0.350*exp((1.0 - mp)/0.012)
                lrep = log(re)
                factor = f1M*(lrep - C2)*(lrep - C3)/((C1 - C2)*(C1 - C3)) &
                         + f2M*(lrep - C1)*(lrep - C3)/((C2 - C1)*(C2 - C3)) &
                         + f3M*(lrep - C1)*(lrep - C2)/((C3 - C1)*(C3 - C2))
                rcd1 = rcd_M1 + (rcd_M2 - rcd_M1)*factor
            else
                rcd1 = (1.0 + 0.107*re**0.867) + &
                       (re/24.0)*0.646/(1.0 + 861.0/re**.634)
            end if ! mp
        end if    ! re

        ! Osnes's volume fraction correction
        b1 = 5.81*phi/((1.0 - phi)**2) + &
             0.48*(phi**(1._wp/3._wp))/((1.0 - phi)**3)

        b2 = ((1.0 - phi)**2)*(phi**3)* &
             re*(0.95 + 0.61*(phi**3)/((1.0 - phi)*2))

        b3 = min(sqrt(20.0_wp*mp), 1.0_wp)* &
             (5.65*phi - 22.0*(phi**2) + 23.4*(phi**3))* &
             (1 + tanh((mp - (0.65 - 0.24*phi))/0.35))

        cd = (24.0/re)*rcd1

        cd = cd/(1.0 - phi) + b3 + (24.0/re)*(1.0 - phi)*(b1 + b2)

        beta = 3.0*pi*mu_fluid*dp*(re/24.0)*cd

    end function QS_ModifiedParmar

    ! QS Force calculated as a function of Re, Ma and phi
    !
    ! Use Osnes etal (2023) correlations
    ! A.N. Osnes, M. Vartdal, M. Khalloufi,
    !    J. Capecelatro, and S. Balachandar.
    ! Comprehensive quasi-steady force correlations for compressible flow
    !    through random particle suspensions.
    ! International Journal of Multiphase Flow, Vol. 165, 104485, (2023).
    ! doi: https://doi.org/10.1016/j.imultiphaseflow.2023.104485.
    !
    ! E. Loth, J.T. Daspit, M. Jeong, T. Nagata, and T. Nonomura.
    ! Supersonic and hypersonic drag coefficients for a sphere.
    ! AIAA Journal, Vol. 59(8), pp. 3261-3274, (2021).
    ! doi: https://doi.org/10.2514/1.J060153.
    !
    ! NOTE: Re<45 Rarefied fomu_fluidla of Loth et al has been redefined by Balachandar
    ! to avoid singularity as Ma -> 0.
    function QS_Osnes(rho, cson, mu_fluid, gamma, vmag, dp, volume_fraction) result(beta)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: rho, cson, mu_fluid, gamma, vmag, dp, volume_fraction
        real(wp) :: rmachp, mp, phi, rep, re
        real(wp) :: Knp, fKn, CD1, s, JM, CD2, cd_loth, CM, GM, HM, b1, b2, b3, cd, sgby2, JMt
        real(wp) :: beta

        rmachp = vmag/cson
        mp = max(rmachp, 0.01_wp)
        phi = max(volume_fraction, 0.0001_wp)
        rep = vmag*dp*rho/mu_fluid
        re = max(rep, 0.1_wp)

        ! Loth's correlation
        if (re <= 45.0_wp) then
            ! Rarefied-dominated regime
            Knp = sqrt(0.5_wp*pi*gamma)*mp/re
            if (Knp > 0.01_wp) then
                fKn = 1.0_wp/(1.0_wp + Knp*(2.514_wp + 0.8_wp*exp(-0.55_wp/Knp)))
            else
                fKn = 1.0_wp/(1.0_wp + Knp*(2.514_wp + 0.8_wp*exp(-0.55_wp/0.01_wp)))

            end if
            CD1 = (24.0_wp/re)*(1.0_wp + 0.15_wp*re**(0.687_wp))*fKn
            s = mp*sqrt(0.5_wp*gamma)
            sgby2 = sqrt(0.5_wp*gamma)
            if (mp <= 1) then
                !JMt = 2.26_wp*(mp**4) - 0.1_wp*(mp**3) + 0.14_wp*mp
                JMt = 2.26_wp*(mp**4) + 0.14_wp*mp
            else
                JMt = 1.6_wp*(mp**4) + 0.25_wp*(mp**3) + 0.11_wp*(mp**2) + 0.44_wp*mp
            end if
            !
            ! Refomu_fluidlated version of Loth et al. to avoid singularity at mp = 0
            !
            CD2 = (1.0_wp + 2.0_wp*(s**2))*exp(-s**2)*mp &
                  /((sgby2**3)*sqrt(pi)) &
                  + (4.0_wp*(s**4) + 4.0_wp*(s**2) - 1.0_wp) &
                  *erf(s)/(2.0_wp*(sgby2**4)) &
                  + (2.0_wp*(mp**3)/(3.0_wp*sgby2))*sqrt(pi)

            CD2 = CD2/(1.0_wp + (((CD2/JMt) - 1.0_wp)*sqrt(re/45.0_wp)))
            cd_loth = CD1/(1.0_wp + (mp**4)) &
                      + CD2/(1.0_wp + (mp**4))
        else
            ! Compression-dominated regime
            ! TLJ: coefficients tweaked to get continuous values
            !      on the two branches at the critical points
            if (mp < 1.5_wp) then
                CM = 1.65_wp + 0.65_wp*tanh(4._wp*mp - 3.4_wp)
            else
                !CM = 2.18_wp - 0.13_wp*tanh(0.9_wp*mp - 2.7_wp)
                CM = 2.18_wp - 0.12913149918318745_wp*tanh(0.9_wp*mp - 2.7_wp)
            end if
            if (mp < 0.8) then
                GM = 166.0_wp*(mp**3) + 3.29_wp*(mp**2) - 10.9_wp*mp + 20._wp
            else
                !GM = 5.0_wp + 40._wp*(mp**(-3))
                GM = 5.0_wp + 47.809331200000017_wp*(mp**(-3))
            end if
            if (mp < 1) then
                HM = 0.0239_wp*(mp**3) + 0.212_wp*(mp**2) &
                     - 0.074_wp*mp + 1._wp
            else
                !HM =   0.93_wp + 1.0_wp / (3.5_wp + (mp**5))
                HM = 0.93967777777777772_wp + 1.0_wp/(3.5_wp + (mp**5))
            end if

            cd_loth = (24.0_wp/re)*(1._wp + 0.15_wp*(re**(0.687_wp)))*HM + &
                      0.42_wp*CM/(1._wp + 42500._wp/re**(1.16_wp*CM) + GM/sqrt(re))

        end if

        b1 = 5.81_wp*phi/((1.0_wp - phi)**2) + &
             0.48_wp*(phi**(1._wp/3._wp))/((1.0_wp - phi)**3)

        b2 = ((1.0_wp - phi)**2)*(phi**3)* &
             re*(0.95_wp + 0.61_wp*(phi**3)/((1.0_wp - phi)*2))

        b3 = min(sqrt(20.0_wp*mp), 1.0_wp)* &
             (5.65_wp*phi - 22.0_wp*(phi**2) + 23.4_wp*(phi**3))* &
             (1._wp + tanh((mp - (0.65_wp - 0.24_wp*phi))/0.35_wp))

        cd = cd_loth/(1.0_wp - phi) + b3 + (24.0_wp/re)*(1.0_wp - phi)*(b1 + b2)

        beta = 3.0_wp*pi*mu_fluid*dp*(re/24.0_wp)*cd

    end function QS_Osnes

! Subroutine for Quasi-Steady Drag Model of Gidaspow
!
! D. Gidaspow, Multiphase Flow and Fluidization (Academic Press, 1994)
!
! Note: Model is provided per cell volume. We convert that to per
! particle using the particle volume fraction and volume
    function QS_Gidaspow(rho, cson, mu_fluid, gamma, vmag, dp, volume_fraction) result(beta)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: rho, cson, mu_fluid, gamma, vmag, dp, volume_fraction
        real(wp) :: cd, phifRep, phif
        real(wp) :: phi, rep, re
        real(wp) :: beta

        rep = vmag*dp*rho/mu_fluid
        phi = max(volume_fraction, 0.0001_wp)
        phif = max(1._wp - volume_fraction, 0.0001_wp)
        re = max(rep, 0.1_wp)

        phifRep = phif*re

        if (phifRep < 1000.0) then
            cd = 24.0/phifRep*(1.0 + 0.15*(phifRep)**0.687)
        else
            cd = 0.44
        end if

        if (phif < 0.8) then
            beta = 150.0*((phi**2)*mu_fluid)/(phif*dp**2) &
                   + 1.75*(rho*phi*vmag/dp)
        else
            beta = 0.75*cd*phi*rho*vmag/(dp*phif**1.65)
        end if

        beta = beta*(pi*dp**3)/(6.0*phi)

    end function QS_Gidaspow

end module m_particles_EL_kernels
