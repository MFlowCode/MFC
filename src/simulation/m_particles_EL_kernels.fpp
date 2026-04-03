!>
!! @file m_particles_EL_kernels.fpp
!! @brief Contains module m_particles_EL_kernels

#:include 'macros.fpp'

!> @brief This module contains kernel functions used to map the effect of the lagrangian particles in the Eulerian framework.
module m_particles_EL_kernels

    use m_mpi_proxy      !< Message passing interface (MPI) module proxy
    use ieee_arithmetic  !< For checking NaN
    use m_model          !< For getting f_model_random_number

    implicit none

    integer, parameter :: alphaf_id_loc = 1
    integer, parameter :: alphaupx_id_loc = 2   !< x particle momentum index
    integer, parameter :: alphaupy_id_loc = 3   !< y particle momentum index
    integer, parameter :: alphaupz_id_loc = 4   !< z particle momentum index
    integer, parameter :: alphaup2x_id_loc = 5  !< x particle velocity squared
    integer, parameter :: alphaup2y_id_loc = 6  !< y particle velocity squared
    integer, parameter :: alphaup2z_id_loc = 7  !< z particle velocity squared
    integer, parameter :: Smx_id_loc = 8
    integer, parameter :: Smy_id_loc = 9
    integer, parameter :: Smz_id_loc = 10
    integer, parameter :: SE_id_loc = 11
    integer, parameter :: dPx_id_loc = 1        !< Spatial pressure gradient in x, y, and z
    integer, parameter :: drhox_id_loc = 4      !< Spatial density gradient in x, y, and z
    integer, parameter :: dufxdx_id_loc = 7  ! du_x/dx
    integer, parameter :: dufxdy_id_loc = 8  ! du_x/dy
    integer, parameter :: dufxdz_id_loc = 9  ! du_x/dz
    integer, parameter :: dufydx_id_loc = 10  ! du_y/dx
    integer, parameter :: dufydy_id_loc = 11  ! du_y/dy
    integer, parameter :: dufydz_id_loc = 12  ! du_y/dz
    integer, parameter :: dufzdx_id_loc = 13  ! du_z/dx
    integer, parameter :: dufzdy_id_loc = 14  ! du_z/dy
    integer, parameter :: dufzdz_id_loc = 15  ! du_z/dz

contains

    ! !> The purpose of this subroutine is to compute each particles total contribution to the gaussian for proper normalization
    subroutine s_compute_gaussian_contribution(rad, pos, cell, func_s, func_s_sources, updatedvar_old)

        $:GPU_ROUTINE(function_name='s_compute_gaussian_contribution',parallelism='[seq]', cray_inline=True)

        type(scalar_field), dimension(:), intent(in) :: updatedvar_old
        real(wp), intent(in) :: rad
        real(wp), intent(in), dimension(3) :: pos
        integer, intent(in), dimension(3) :: cell
        real(wp), intent(out) :: func_s, func_s_sources
        real(wp) :: volpart, stddsv, Vol_loc, func, alpha_f
        real(wp), dimension(3) :: nodecoord, center
        integer :: ip, jp, kp, di, dj, dk, di_beg, di_end, dj_beg, dj_end, dk_beg, dk_end, mapCells_loc
        integer, dimension(3) :: cellijk

        mapCells_loc = 1

        volpart = (4._wp/3._wp)*pi*rad**3

        call s_compute_stddsv(cell, volpart, stddsv)

        ip = cell(1)
        jp = cell(2)
        kp = cell(3)

        di_beg = ip - mapCells_loc
        di_end = ip + mapCells_loc
        dj_beg = jp - mapCells_loc
        dj_end = jp + mapCells_loc
        dk_beg = kp
        dk_end = kp

        if (num_dims == 3) then
            dk_beg = kp - mapCells_loc
            dk_end = kp + mapCells_loc
        end if

        func_s = 0._wp
        func_s_sources = 0._wp
        do dk = dk_beg, dk_end
            do dj = dj_beg, dj_end
                do di = di_beg, di_end
                    nodecoord(1) = x_cc(di)
                    nodecoord(2) = y_cc(dj)
                    nodecoord(3) = 0._wp
                    if (p > 0) nodecoord(3) = z_cc(dk)

                    cellijk(1) = di
                    cellijk(2) = dj
                    cellijk(3) = dk

                    center(1:2) = pos(1:2)
                    center(3) = 0._wp
                    if (p > 0) center(3) = pos(3)

                    Vol_loc = dx(cellijk(1))*dy(cellijk(2))
                    if (num_dims == 3) Vol_loc = dx(cellijk(1))*dy(cellijk(2))*dz(cellijk(3))

                    call s_applygaussian(center, cellijk, nodecoord, stddsv, 0._wp, func)

                    alpha_f = updatedvar_old(alphaf_id_loc)%sf(cellijk(1), cellijk(2), cellijk(3))
                    func_s = func_s + (func*Vol_loc)
                    func_s_sources = func_s_sources + (func*Vol_loc*alpha_f)
                end do
            end do
        end do

    end subroutine s_compute_gaussian_contribution

    !> The purpose of this subroutine is to compute the gaussian smearing of particle volume fraction and source terms with atomic
    !! cell updates
    subroutine s_gaussian_atomic(rad, vel, pos, force_p, gauSum, gauSum_sources, cell, updatedvar, kcomp, ind_start, ind_end)

        $:GPU_ROUTINE(function_name='s_gaussian_atomic',parallelism='[seq]', cray_inline=True)

        real(wp), intent(in) :: rad, gauSum, gauSum_sources
        real(wp), intent(in), dimension(3) :: pos, vel, force_p
        integer, intent(in), dimension(3) :: cell
        type(scalar_field), dimension(:), intent(inout) :: updatedvar
        type(scalar_field), dimension(:), intent(inout) :: kcomp
        integer, intent(in) :: ind_start, ind_end
        real(wp) :: volpart, stddsv, Vol_loc, func, weight
        real(wp) :: fp_x, fp_y, fp_z, vp_x, vp_y, vp_z
        real(wp) :: addFun
        real(wp), dimension(3) :: nodecoord, center
        integer :: ip, jp, kp, di, dj, dk, di_beg, di_end, dj_beg, dj_end, dk_beg, dk_end, mapCells_loc, field_ind
        integer, dimension(3) :: cellijk

        mapCells_loc = 1

        volpart = (4._wp/3._wp)*pi*rad**3._wp

        call s_compute_stddsv(cell, volpart, stddsv)

        ip = cell(1)
        jp = cell(2)
        kp = cell(3)

        di_beg = ip - mapCells_loc
        di_end = ip + mapCells_loc
        dj_beg = jp - mapCells_loc
        dj_end = jp + mapCells_loc
        dk_beg = kp
        dk_end = kp

        if (num_dims == 3) then
            dk_beg = kp - mapCells_loc
            dk_end = kp + mapCells_loc
        end if

        fp_x = -force_p(1)
        fp_y = -force_p(2)
        fp_z = -force_p(3)

        vp_x = vel(1)
        vp_y = vel(2)
        vp_z = vel(3)

        center(1:2) = pos(1:2)
        center(3) = 0._wp
        if (p > 0) center(3) = pos(3)

        do dk = dk_beg, dk_end
            do dj = dj_beg, dj_end
                do di = di_beg, di_end
                    nodecoord(1) = x_cc(di)
                    nodecoord(2) = y_cc(dj)
                    nodecoord(3) = 0._wp
                    if (p > 0) nodecoord(3) = z_cc(dk)

                    cellijk(1) = di
                    cellijk(2) = dj
                    cellijk(3) = dk

                    Vol_loc = dx(cellijk(1))*dy(cellijk(2))
                    if (num_dims == 3) Vol_loc = Vol_loc*dz(cellijk(3))

                    call s_applygaussian(center, cellijk, nodecoord, stddsv, 0._wp, func)

                    weight = func/gauSum

                    do field_ind = ind_start, ind_end
                        if (field_ind == alphaf_id_loc) then
                            addFun = weight*volpart
                        else if (field_ind == alphaupx_id_loc) then
                            addFun = weight*volpart*vp_x
                        else if (field_ind == alphaupy_id_loc) then
                            addFun = weight*volpart*vp_y
                        else if (field_ind == alphaupz_id_loc) then
                            addFun = weight*volpart*vp_z
                        else if (field_ind == alphaup2x_id_loc) then
                            addFun = weight*volpart*vp_x**2
                        else if (field_ind == alphaup2y_id_loc) then
                            addFun = weight*volpart*vp_y**2
                        else if (field_ind == alphaup2z_id_loc) then
                            addFun = weight*volpart*vp_z**2
                        else if (field_ind == Smx_id_loc) then
                            weight = func/gauSum_sources
                            addFun = weight*fp_x
                        else if (field_ind == Smy_id_loc) then
                            weight = func/gauSum_sources
                            addFun = weight*fp_y
                        else if (field_ind == Smz_id_loc) then
                            weight = func/gauSum_sources
                            addFun = weight*fp_z
                        else if (field_ind == SE_id_loc) then
                            weight = func/gauSum_sources
                            addFun = 0._wp
                        end if

                        $:GPU_ATOMIC(atomic='update')
                        updatedvar(field_ind)%sf(cellijk(1), cellijk(2), cellijk(3)) = updatedvar(field_ind)%sf(cellijk(1), &
                                   & cellijk(2), cellijk(3)) + real(addFun, kind=stp)
                    end do
                end do
            end do
        end do

    end subroutine s_gaussian_atomic

    !> The purpose of this subroutine is to apply the gaussian kernel function for each particle (Maeda and Colonius, 2018)).
    subroutine s_applygaussian(center, cellaux, nodecoord, stddsv, strength_idx, func)

        $:GPU_ROUTINE(function_name='s_applygaussian',parallelism='[seq]', cray_inline=True)

        real(wp), dimension(3), intent(in) :: center
        integer, dimension(3), intent(in)  :: cellaux
        real(wp), dimension(3), intent(in) :: nodecoord
        real(wp), intent(in)               :: stddsv
        real(wp), intent(in)               :: strength_idx
        real(wp), intent(out)              :: func
        integer                            :: i
        real(wp)                           :: distance
        real(wp)                           :: theta, dtheta, L2, dzp, Lz2, zc
        real(wp)                           :: Nr, Nr_count

        distance = sqrt((center(1) - nodecoord(1))**2._wp + (center(2) - nodecoord(2))**2._wp + (center(3) - nodecoord(3))**2._wp)

        if (num_dims == 3) then
            !> 3D gaussian function
            func = exp(-0.5_wp*(distance/stddsv)**2._wp)/(sqrt(2._wp*pi)*stddsv)**3._wp
        else
            if (cyl_coord) then
                !> 2D cylindrical function:
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
                    func = func + dtheta/2._wp/pi*exp(-0.5_wp*(distance/stddsv)**2._wp)/(sqrt(2._wp*pi)*stddsv) &
                                                      & **(3._wp*(strength_idx + 1._wp))
                end do
            else
                !> 2D cartesian function: Equation (48) from Maeda and Colonius 2018
                ! We smear particles considering a virtual depth (lag_params%charwidth) with lag_params%charNz cells
                dzp = (lag_params%charwidth/(lag_params%charNz + 1._wp))

                func = 0._wp
                do i = 0, lag_params%charNz
                    zc = (-lag_params%charwidth/2._wp + dzp*(0.5_wp + i))  ! Center of virtual cell i in z-direction
                    Lz2 = (center(3) - zc)**2._wp
                    distance = sqrt((center(1) - nodecoord(1))**2._wp + (center(2) - nodecoord(2))**2._wp + Lz2)
                    func = func + dzp/lag_params%charwidth*exp(-0.5_wp*(distance/stddsv)**2._wp)/(sqrt(2._wp*pi)*stddsv)**3._wp
                end do
            end if
        end if

    end subroutine s_applygaussian

    !> Calculates the standard deviation of the particle being smeared in the Eulerian framework.
    !! @param cell Cell where the particle is located
    !! @param volpart Volume of the particle
    !! @param stddsv Standard deviaton
    subroutine s_compute_stddsv(cell, volpart, stddsv)

        $:GPU_ROUTINE(function_name='s_compute_stddsv',parallelism='[seq]', cray_inline=True)

        integer, dimension(3), intent(in) :: cell
        real(wp), intent(in)              :: volpart
        real(wp), intent(out)             :: stddsv
        real(wp)                          :: chardist, charvol
        real(wp)                          :: rad

        !> Compute characteristic distance
        chardist = sqrt(dx(cell(1))*dy(cell(2)))
        if (p > 0) chardist = (dx(cell(1))*dy(cell(2))*dz(cell(3)))**(1._wp/3._wp)

        !> Compute characteristic volume
        if (p > 0) then
            charvol = dx(cell(1))*dy(cell(2))*dz(cell(3))
        else
            if (cyl_coord) then
                charvol = dx(cell(1))*dy(cell(2))*y_cc(cell(2))*2._wp*pi
            else
                charvol = dx(cell(1))*dy(cell(2))*lag_params%charwidth
            end if
        end if

        rad = (3._wp*volpart/(4._wp*pi))**(1._wp/3._wp)
        stddsv = lag_params%epsilonb*max(chardist, rad)

    end subroutine s_compute_stddsv

    !> The purpose of this subroutine is to check if the current cell is outside the computational domain or not (including ghost
    !! cells).
    !! @param cellaux Tested cell to smear the particle effect in.
    !! @param celloutside If true, then cellaux is outside the computational domain.
    subroutine s_check_celloutside(cellaux, celloutside)

        $:GPU_ROUTINE(function_name='s_check_celloutside',parallelism='[seq]', cray_inline=True)

        integer, dimension(3), intent(inout) :: cellaux
        logical, intent(out)                 :: celloutside

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

    !> This subroutine transforms the computational coordinates of the particle from real type into integer.
    !! @param s Computational coordinates of the particle, real type
    !! @param get_cell Computational coordinates of the particle, integer type
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

    !> The purpose of this procedure is to calculate the characteristic cell volume
    !! @param cell Computational coordinates (x, y, z)
    !! @param Charvol Characteristic volume
    subroutine s_get_char_vol(cellx, celly, cellz, Charvol)

        $:GPU_ROUTINE(function_name='s_get_char_vol',parallelism='[seq]', cray_inline=True)

        integer, intent(in)   :: cellx, celly, cellz
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

    !! This function calculates the force on a particle based on the pressure gradient, velocity, and drag model.
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
    subroutine s_get_particle_force(pos, rad, vel_p, mass_p, Re, gamm, seed, fqsfluct, cell, q_prim_vf, q_cons_vf, q_particles, &
                                    & fieldvars, rhs_old, duidxj_id_loc, wx, wy, wz, force, rmass_add, new_seed, new_fqsfluct)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in)                                :: rad, mass_p, Re, gamm
        real(wp), dimension(3), intent(in)                  :: pos
        integer, dimension(3), intent(in)                   :: cell
        integer, intent(in)                                 :: seed
        real(wp), dimension(3), intent(in)                  :: vel_p, fqsfluct
        integer, dimension(3, 3), intent(in)                :: duidxj_id_loc
        type(scalar_field), dimension(:), intent(in)        :: q_particles
        type(scalar_field), dimension(:), intent(in)        :: fieldvars
        type(scalar_field), dimension(:), intent(in)        :: wx, wy, wz
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(in) :: rhs_old
        real(wp), dimension(3), intent(out)                 :: force, new_fqsfluct
        real(wp), intent(out)                               :: rmass_add
        integer, intent(out)                                :: new_seed
        integer                                             :: seed_loc
        real(wp)                                            :: a, vol, rho_fluid, pressure_fluid, alpha_f
        real(wp), dimension(3)                              :: v_rel, dp
        real(wp)                                            :: particle_diam, gas_mu, vmag, cson
        real(wp)                                            :: slip_velocity_x, slip_velocity_y, slip_velocity_z, beta
        real(wp)                                            :: vol_frac
        real(wp), dimension(3)                              :: fluid_vel
        integer                                             :: dir, l

        ! Added pass params
        real(wp)               :: mach, Cam, flux_f, flux_b, SDrho, vrel_gradrho, drhodt
        real(wp), dimension(3) :: rhoDuDt, grad_rho, fam, udot_gradu
        integer, dimension(3)  :: p1

        ! QS Fluct
        real(wp), dimension(3) :: vel_p_mean, vel2_p_mean

        ! Suth
        real(wp) :: fluid_temp, R_fluid, suth, tref

        force = 0._wp
        dp = 0._wp
        grad_rho = 0._wp
        fam = 0._wp
        fluid_vel = 0._wp
        v_rel = 0._wp
        rhoDuDt = 0._wp
        SDrho = 0._wp

        udot_gradu = 0._wp

        vel_p_mean = 0._wp
        vel2_p_mean = 0._wp

        !!Interpolation - either even ordered barycentric or 0th order
        if (lag_params%interpolation_order > 1) then
            rho_fluid = f_interp_barycentric(pos, cell, q_prim_vf, 1, wx, wy, wz)
            pressure_fluid = f_interp_barycentric(pos, cell, q_prim_vf, E_idx, wx, wy, wz)
            alpha_f = f_interp_barycentric(pos, cell, q_particles, alphaf_id_loc, wx, wy, wz)
            vol_frac = 1._wp - alpha_f

            do dir = 1, num_dims
                vel_p_mean(dir) = f_interp_barycentric(pos, cell, q_particles, alphaupx_id_loc + dir - 1, wx, wy, &
                           & wz)/max(vol_frac, 1.e-12_wp)
                vel2_p_mean(dir) = f_interp_barycentric(pos, cell, q_particles, alphaup2x_id_loc + dir - 1, wx, wy, &
                            & wz)/max(vol_frac, 1.e-12_wp)
                fluid_vel(dir) = f_interp_barycentric(pos, cell, q_prim_vf, momxb + dir - 1, wx, wy, wz)
            end do

            if (lag_params%added_mass_model > 0) then
                drhodt = rhs_old(1)%sf(cell(1), cell(2), cell(3))
            end if

            do dir = 1, num_dims
                if (lag_params%pressure_force .or. lag_params%added_mass_model > 0) then
                    dp(dir) = f_interp_barycentric(pos, cell, fieldvars, dPx_id_loc + dir - 1, wx, wy, wz)
                end if
                if (lag_params%added_mass_model > 0) then
                    grad_rho(dir) = f_interp_barycentric(pos, cell, fieldvars, drhox_id_loc + dir - 1, wx, wy, wz)
                    rhoDuDt(dir) = (rhs_old(momxb + dir - 1)%sf(cell(1), cell(2), cell(3)) - fluid_vel(dir)*drhodt)/rho_fluid
                    do l = 1, num_dims
                        udot_gradu(dir) = udot_gradu(dir) + fluid_vel(l)*f_interp_barycentric(pos, cell, fieldvars, &
                                   & duidxj_id_loc(dir, l), wx, wy, wz)
                    end do
                end if
            end do
        else
            rho_fluid = q_prim_vf(1)%sf(cell(1), cell(2), cell(3))
            pressure_fluid = q_prim_vf(E_idx)%sf(cell(1), cell(2), cell(3))
            alpha_f = q_particles(alphaf_id_loc)%sf(cell(1), cell(2), cell(3))
            vol_frac = 1._wp - alpha_f

            do dir = 1, num_dims
                vel_p_mean(dir) = q_particles(alphaupx_id_loc + dir - 1)%sf(cell(1), cell(2), cell(3))/max(vol_frac, 1.e-12_wp)
                vel2_p_mean(dir) = q_particles(alphaup2x_id_loc + dir - 1)%sf(cell(1), cell(2), cell(3))/max(vol_frac, 1.e-12_wp)
                fluid_vel(dir) = q_prim_vf(momxb + dir - 1)%sf(cell(1), cell(2), cell(3))
            end do

            if (lag_params%added_mass_model > 0) then
                drhodt = rhs_old(1)%sf(cell(1), cell(2), cell(3))
            end if

            do dir = 1, num_dims
                if (lag_params%pressure_force .or. lag_params%added_mass_model > 0) then
                    dp(dir) = fieldvars(dPx_id_loc + dir - 1)%sf(cell(1), cell(2), cell(3))
                end if
                if (lag_params%added_mass_model > 0) then
                    grad_rho(dir) = fieldvars(drhox_id_loc + dir - 1)%sf(cell(1), cell(2), cell(3))
                    rhoDuDt(dir) = (rhs_old(momxb + dir - 1)%sf(cell(1), cell(2), cell(3)) - fluid_vel(dir)*drhodt)/rho_fluid
                    do l = 1, num_dims
                        udot_gradu(dir) = udot_gradu(dir) + fluid_vel(l)*fieldvars(duidxj_id_loc(dir, l))%sf(cell(1), cell(2), &
                                   & cell(3))
                    end do
                end if
            end do
        end if

        v_rel = vel_p - fluid_vel

        if (lag_params%qs_drag_model > 0 .or. lag_params%added_mass_model > 0) then
            ! Quasi-steady Drag Force Parameters
            slip_velocity_x = fluid_vel(1) - vel_p(1)
            slip_velocity_y = fluid_vel(2) - vel_p(2)
            if (num_dims == 3) then
                slip_velocity_z = fluid_vel(3) - vel_p(3)
                vmag = sqrt(slip_velocity_x*slip_velocity_x + slip_velocity_y*slip_velocity_y + slip_velocity_z*slip_velocity_z)
            else if (num_dims == 2) then
                vmag = sqrt(slip_velocity_x*slip_velocity_x + slip_velocity_y*slip_velocity_y)
            end if
            particle_diam = rad*2._wp
            if (rho_fluid > 0._wp) then
                cson = sqrt((gamm*pressure_fluid)/rho_fluid)
            else
                cson = 1._wp
            end if
            gas_mu = Re

            if (.not. viscous) then
                tref = 273.15_wp
                R_fluid = 287.05_wp  ! J/kg-K for Air
                suth = 110.4_wp
                fluid_temp = pressure_fluid/(rho_fluid*R_fluid)

                gas_mu = gas_mu*sqrt(fluid_temp/tref)*(1.0_wp + suth/tref)/(1.0_wp + suth/fluid_temp)
            end if
        end if

        if (lag_params%added_mass_model > 0) then
            rhoDuDt = rho_fluid*(rhoDuDt + udot_gradu)
            vrel_gradrho = dot_product(-v_rel, grad_rho)
            SDrho = drhodt + vel_p(1)*grad_rho(1) + vel_p(2)*grad_rho(2) + vel_p(3)*grad_rho(3)
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
            ! No Quasi-Steady drag
        end if

        ! Step 1.1: Stokes drag
        if (lag_params%stokes_drag == 1) then  ! Free slip Stokes drag
            force = force - 4._wp*pi*gas_mu*rad*v_rel
        else if (lag_params%stokes_drag == 2) then  ! No slip Stokes drag
            force = force - 6._wp*pi*gas_mu*rad*v_rel
        else
            ! No stokes drag
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
                Cam = 1._wp + 1.8_wp*(0.6_wp**2) + 7.6_wp*(0.6_wp**4)
            else
                Cam = 1._wp + 1.8_wp*mach**2 + 7.6_wp*mach**4
            end if

            Cam = 0.5_wp*Cam*(1._wp + 0.68_wp*vol_frac**2)
            rmass_add = rho_fluid*vol*Cam

            fam = Cam*vol*(-v_rel*SDrho + rhoDuDt + fluid_vel*(vrel_gradrho))

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

        if (lag_params%qs_fluct_force) then
            ! Step 5: QS Fluctuations
            !> New addition: QS Fluctuations
            seed_loc = seed
            call s_compute_qs_fluctuations(vel_p, fluid_vel, rho_fluid, cson, gas_mu, particle_diam, vol_frac, vmag, vel_p_mean, &
                                           & vel2_p_mean, seed_loc, fqsfluct, dt, new_fqsfluct)
            force = force + new_fqsfluct
            new_seed = seed_loc
            !>
        end if

        do dir = 1, num_dims
            if (.not. ieee_is_finite(force(dir))) then
                force(dir) = 0._wp
            end if
        end do

    end subroutine s_get_particle_force

    ! Quasi-steady force fluctuations Osnes, Vartdal, Khalloufi, Capecelatro (2023) Comprehensive quasi-steady force correlations
    ! for compressible flow through random particle suspensions. International Journal of Multiphase Flows, Vo. 165, 104485.
    ! Lattanzi, Tavanashad, Subramaniam, Capecelatro (2022) Stochastic model for the hydrodynamic force in Euler-Lagrange
    ! silumations of particle-laden flows. Physical Review Fluids, Vol. 7, 014301. Note: To compute the granular temperature, we
    ! assume the velocity fluctuations are uncorrelated. Note: The means are computed using a box filter with an adaptive filter
    ! width. Compute mean using box filter for langevin model - not for feedback
    !
    ! The mean is calculated according to Lattanzi etal, Physical Review Fluids, 2022.
    subroutine s_compute_qs_fluctuations(vel_p, fluid_vel, rho_fluid, cson, gas_mu, particle_diam, vol_frac, vmag, vel_p_mean, &
                                         & vel2_p_mean, seed, fqs_fluct_old, dt_loc, fqs_fluct_new)
        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), dimension(3), intent(in)  :: vel_p, fluid_vel
        real(wp), intent(in)                :: rho_fluid, cson, gas_mu, particle_diam, vol_frac, vmag
        real(wp), dimension(3), intent(in)  :: vel_p_mean, vel2_p_mean
        integer, intent(inout)              :: seed
        real(wp), dimension(3), intent(in)  :: fqs_fluct_old
        real(wp), intent(in)                :: dt_loc
        real(wp), dimension(3), intent(out) :: fqs_fluct_new
        real(wp)                            :: upmean, vpmean, wpmean
        real(wp)                            :: u2pmean, v2pmean, w2pmean
        real(wp)                            :: rphip, rep, rmachp, theta, chi, tF_inv
        real(wp)                            :: fq, bq, Fs, sigD, sigT
        real(wp)                            :: aSDE, bSDE_CD, bSDE_CL
        real(wp)                            :: CD_prime, CD_frac, sigmoid_cf, f_CF
        real(wp)                            :: Z1, Z2, dW1, dW2, cosrand, sinrand
        real(wp), dimension(3)              :: avec, bvec, cvec, dvec, eunit
        real(wp), dimension(3)              :: slip_vel
        real(wp)                            :: denum, TwoPi
        real(wp), dimension(5)              :: UnifRnd

        fqs_fluct_new = 0._wp

        upmean = vel_p_mean(1)
        vpmean = vel_p_mean(2)
        wpmean = vel_p_mean(3)

        u2pmean = vel2_p_mean(1)
        v2pmean = vel2_p_mean(2)
        w2pmean = vel2_p_mean(3)

        TwoPi = 2._wp*pi

        ! Particle phase properties
        rphip = vol_frac  ! particle volume fraction
        slip_vel = fluid_vel - vel_p

        ! Particle Reynolds number
        rep = rho_fluid*vmag*particle_diam/gas_mu

        ! Particle Mach number
        rmachp = vmag/cson

        ! Granular temperature from Eulerian fields \theta = (<u_p^2> - <u_p>^2) / 3
        theta = ((u2pmean + v2pmean + w2pmean) - (upmean**2 + vpmean**2 + wpmean**2))/3._wp

        if (theta <= 1.e-12_wp) theta = 0._wp

        ! Fluctuating drag magnitude (Osnes Eqs 9-12)
        fq = 6.52_wp*rphip - 22.56_wp*(rphip**2) + 49.90_wp*(rphip**3)
        Fs = 3._wp*pi*gas_mu*particle_diam*(1._wp + 0.15_wp*((rep*(1._wp - rphip))**0.687_wp))*(1._wp - rphip)*vmag
        bq = min(sqrt(20._wp*rmachp), 1._wp)*0.55_wp*(rphip**0.7_wp)*(1._wp + tanh((rmachp - 0.5_wp)/0.2_wp))
        sigD = (fq + bq)*Fs

        ! Granular kinetic theory relaxation rate
        chi = (1._wp + 2.50_wp*rphip + 4.51_wp*(rphip**2) + 4.52_wp*(rphip**3))/((1._wp - (rphip/0.64_wp)**3)**0.68_wp)
        tF_inv = (24._wp*rphip*chi/particle_diam)*sqrt(theta/pi)

        ! SDE coefficients
        aSDE = tF_inv
        bSDE_CD = sigD*sqrt(2._wp*tF_inv)

        ! Unit slip velocity direction
        if (vmag > 1.e-8_wp) then
            avec = slip_vel/vmag

            ! Project old fluctuating force onto slip direction
            CD_prime = fqs_fluct_old(1)*avec(1) + fqs_fluct_old(2)*avec(2) + fqs_fluct_old(3)*avec(3)
            CD_frac = CD_prime/max(sigD, 1.e-30_wp)
        else
            avec = [1._wp, 0._wp, 0._wp]
            CD_prime = 0._wp
            sigD = 0._wp
            CD_frac = 0._wp
        end if

        ! Perpendicular fluctuation magnitude
        sigmoid_cf = 1._wp/(1._wp + exp(-CD_frac))
        f_CF = 0.39356905_wp*sigmoid_cf + 0.43758848_wp
        sigT = f_CF*sigD
        bSDE_CL = sigT*sqrt(2._wp*tF_inv)

        ! Build orthogonal basis
        eunit = [1._wp, 0._wp, 0._wp]
        if (abs(avec(2)) + abs(avec(3)) <= 1.e-8_wp) then
            eunit = [0._wp, 1._wp, 0._wp]
        else if (abs(avec(1)) + abs(avec(3)) <= 1.e-8_wp) then
            eunit = [0._wp, 0._wp, 1._wp]
        end if

        ! bvec = avec x eunit
        bvec(1) = avec(2)*eunit(3) - avec(3)*eunit(2)
        bvec(2) = avec(3)*eunit(1) - avec(1)*eunit(3)
        bvec(3) = avec(1)*eunit(2) - avec(2)*eunit(1)
        denum = max(1.e-8_wp, sqrt(bvec(1)**2 + bvec(2)**2 + bvec(3)**2))
        bvec = bvec/denum

        ! cvec = avec x bvec
        cvec(1) = avec(2)*bvec(3) - avec(3)*bvec(2)
        cvec(2) = avec(3)*bvec(1) - avec(1)*bvec(3)
        cvec(3) = avec(1)*bvec(2) - avec(2)*bvec(1)
        denum = max(1.e-8_wp, sqrt(cvec(1)**2 + cvec(2)**2 + cvec(3)**2))
        cvec = cvec/denum

        ! Generate random numbers
        UnifRnd(1) = f_model_random_number(seed)
        UnifRnd(2) = f_model_random_number(seed)
        UnifRnd(3) = f_model_random_number(seed)
        UnifRnd(4) = f_model_random_number(seed)
        UnifRnd(5) = f_model_random_number(seed)

        UnifRnd(1) = max(UnifRnd(1), 1.e-30_wp)
        UnifRnd(3) = max(UnifRnd(3), 1.e-30_wp)

        ! Box-Muller transform
        Z1 = sqrt(-2._wp*log(UnifRnd(1)))*cos(TwoPi*UnifRnd(2))
        Z2 = sqrt(-2._wp*log(UnifRnd(3)))*cos(TwoPi*UnifRnd(4))

        ! Scaled stochastic increments
        dW1 = sqrt(dt_loc)*Z1
        dW2 = sqrt(dt_loc)*Z2

        ! Random perpendicular direction
        cosrand = cos(TwoPi*UnifRnd(5))
        sinrand = sin(TwoPi*UnifRnd(5))
        dvec = bvec*cosrand + cvec*sinrand
        denum = max(1.e-8_wp, sqrt(dvec(1)**2 + dvec(2)**2 + dvec(3)**2))
        dvec = dvec/denum

        ! Ornstein-Uhlenbeck update
        fqs_fluct_new(1) = (1._wp - aSDE*dt_loc)*fqs_fluct_old(1) + bSDE_CD*dW1*avec(1) + bSDE_CL*dW2*dvec(1)
        fqs_fluct_new(2) = (1._wp - aSDE*dt_loc)*fqs_fluct_old(2) + bSDE_CD*dW1*avec(2) + bSDE_CL*dW2*dvec(2)
        if (num_dims == 3) then
            fqs_fluct_new(3) = (1._wp - aSDE*dt_loc)*fqs_fluct_old(3) + bSDE_CD*dW1*avec(3) + bSDE_CL*dW2*dvec(3)
        end if

    end subroutine s_compute_qs_fluctuations

    !> Interpolate an Eulerian field to a particle position using barycentric Lagrange interpolation with precomputed weights. Falls
    !! back to nearest-cell value if the interpolant is non-finite.
    !! @param pos Particle position
    !! @param cell Grid cell containing the particle
    !! @param field_vf Eulerian field to interpolate
    !! @param field_index Component index in field_vf
    function f_interp_barycentric(pos, cell, field_vf, field_index, wx, wy, wz) result(val)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), dimension(3), intent(in)           :: pos
        integer, dimension(3), intent(in)            :: cell
        type(scalar_field), dimension(:), intent(in) :: field_vf
        type(scalar_field), dimension(:), intent(in) :: wx, wy, wz
        integer, intent(in)                          :: field_index
        integer                                      :: i, j, k, ix, jy, kz, npts, npts_z, N, a, b
        integer                                      :: ix_count, jy_count, kz_count
        real(wp)                                     :: weight, numerator, denominator, xBar, eps
        real(wp)                                     :: val, local_min, local_max, prod_x, prod_y, prod_z

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

        ! if (abs(pos(1) - x_cc(i)) <= eps .and. & abs(pos(2) - y_cc(j)) <= eps .and. & abs(pos(3) - z_cc(k)) <= eps) then val =
        ! field_vf(field_index)%sf(i, j, k) return end if

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
        else if (abs(val) <= eps) then
            val = 0._wp
        end if

    end function f_interp_barycentric

    ! Quasi-steady force (Re_p and Ma_p corrections): Improved Drag Correlation for Spheres and Application to Shock-Tube
    ! Experiments - Parmar et al. (2010) - AIAA Journal
    !
    ! Quasi-steady force (phi corrections): The Added Mass, Basset, and Viscous Drag Coefficients in Nondilute Bubbly Liquids
    ! Undergoing Small-Amplitude Oscillatory Motion - Sangani et al. (1991) - Phys. Fluids A
    function QS_Parmar(rho, cson, mu_fluid, gamma, vmag, dp, volume_fraction) result(beta)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: rho, cson, mu_fluid, gamma, vmag, dp, volume_fraction
        real(wp)             :: rcd1, rmacr, rcd_mcr, rcd_std, rmach_rat, rcd_M1
        real(wp)             :: rcd_M2, C1, C2, C3, f1M, f2M, f3M, lrep, factor, cd, phi_corr
        real(wp)             :: beta
        real(wp)             :: rmachp, mp, phi, rep, re

        rmachp = vmag/cson
        mp = max(rmachp, 0.01_wp)
        phi = max(volume_fraction, 0.0001_wp)
        rep = vmag*dp*rho/mu_fluid
        re = max(rep, 0.1_wp)

        if (re < 1.e-14_wp) then
            rcd1 = 1.0_wp
        else
            rmacr = 0.6_wp  ! Critical rmachp no
            rcd_mcr = (1._wp + 0.15_wp*re**(0.684_wp)) + (re/24.0_wp)*(0.513_wp/(1._wp + 483._wp/re**(0.669_wp)))
            if (mp <= rmacr) then
                rcd_std = (1._wp + 0.15_wp*re**(0.687_wp)) + (re/24.0_wp)*(0.42_wp/(1._wp + 42500._wp/re**(1.16_wp)))
                rmach_rat = mp/rmacr
                rcd1 = rcd_std + (rcd_mcr - rcd_std)*rmach_rat
            else if (mp <= 1.0_wp) then
                rcd_M1 = (1.0_wp + 0.118_wp*re**0.813_wp) + (re/24.0_wp)*0.69_wp/(1.0_wp + 3550.0_wp/re**0.793_wp)
                C1 = 6.48_wp
                C2 = 9.28_wp
                C3 = 12.21_wp
                f1M = -1.884_wp + 8.422_wp*mp - 13.70_wp*mp**2 + 8.162_wp*mp**3
                f2M = -2.228_wp + 10.35_wp*mp - 16.96_wp*mp**2 + 9.840_wp*mp**3
                f3M = 4.362_wp - 16.91_wp*mp + 19.84_wp*mp**2 - 6.296_wp*mp**3
                lrep = log(re)
                factor = f1M*(lrep - C2)*(lrep - C3)/((C1 - C2)*(C1 - C3)) + f2M*(lrep - C1)*(lrep - C3)/((C2 - C1)*(C2 - C3)) &
                              & + f3M*(lrep - C1)*(lrep - C2)/((C3 - C1)*(C3 - C2))
                rcd1 = rcd_mcr + (rcd_M1 - rcd_mcr)*factor
            else if (mp < 1.75_wp) then
                rcd_M1 = (1.0_wp + 0.118_wp*re**0.813_wp) + (re/24.0_wp)*0.69_wp/(1.0_wp + 3550.0_wp/re**0.793_wp)
                rcd_M2 = (1.0_wp + 0.107_wp*re**0.867_wp) + (re/24.0_wp)*0.646_wp/(1.0_wp + 861.0_wp/re**0.634_wp)
                C1 = 6.48_wp
                C2 = 8.93_wp
                C3 = 12.21_wp
                f1M = -2.963_wp + 4.392_wp*mp - 1.169_wp*mp**2 - 0.027_wp*mp**3 - 0.233_wp*exp((1.0_wp - mp)/0.011_wp)
                f2M = -6.617_wp + 12.11_wp*mp - 6.501_wp*mp**2 + 1.182_wp*mp**3 - 0.174_wp*exp((1.0_wp - mp)/0.010_wp)
                f3M = -5.866_wp + 11.57_wp*mp - 6.665_wp*mp**2 + 1.312_wp*mp**3 - 0.350_wp*exp((1.0_wp - mp)/0.012_wp)
                lrep = log(re)
                factor = f1M*(lrep - C2)*(lrep - C3)/((C1 - C2)*(C1 - C3)) + f2M*(lrep - C1)*(lrep - C3)/((C2 - C1)*(C2 - C3)) &
                              & + f3M*(lrep - C1)*(lrep - C2)/((C3 - C1)*(C3 - C2))
                rcd1 = rcd_M1 + (rcd_M2 - rcd_M1)*factor
            else
                rcd1 = (1.0_wp + 0.107_wp*re**0.867_wp) + (re/24.0_wp)*0.646_wp/(1.0_wp + 861.0_wp/re**0.634_wp)
            end if  ! mp
        end if  ! re

        ! Sangani's volume fraction correction for dilute random arrays Capping volume fraction at 0.5
        phi_corr = (1.0_wp + 5.94_wp*min(phi, 0.5_wp))

        cd = (24.0_wp/re)*rcd1*phi_corr

        beta = rcd1*3.0_wp*pi*mu_fluid*dp

        beta = beta*phi_corr

    end function QS_Parmar

    ! Quasi-steady force (Re_p and Ma_p corrections): Improved Drag Correlation for Spheres and Application to Shock-Tube
    ! Experiments - Parmar et al. (2010) - AIAA Journal
    !
    ! Quasi-steady force (phi corrections): Sangani et al. (1991) volume fraction correction overshoots the drag coefficient.
    !
    ! We adopt instead Osnes et al. (2023) volume fraction correction based on Tenneti et al. with one extra term.
    !
    ! At Mach=0, the drag coefficient from this subroutine matches very well with the one calculated using the Osnes subroutine, for
    ! various Reynolds numbers and volume fractions.
    function QS_ModifiedParmar(rho, cson, mu_fluid, gamma, vmag, dp, volume_fraction) result(beta)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: rho, cson, mu_fluid, gamma, vmag, dp, volume_fraction
        real(wp)             :: rcd1, rmacr, rcd_mcr, rcd_std, rmach_rat, rcd_M1
        real(wp)             :: rcd_M2, C1, C2, C3, f1M, f2M, f3M, lrep, factor, cd, phi_corr
        real(wp)             :: b1, b2, b3
        real(wp)             :: beta
        real(wp)             :: rmachp, mp, phi, rep, re

        rmachp = vmag/cson
        mp = max(rmachp, 0.01_wp)
        phi = max(volume_fraction, 0.0001_wp)
        rep = vmag*dp*rho/mu_fluid
        re = max(rep, 0.1_wp)

        if (re < 1e-14_wp) then
            rcd1 = 1.0_wp
        else
            rmacr = 0.6_wp  ! Critical rmachp no.
            rcd_mcr = (1._wp + 0.15_wp*re**(0.684_wp)) + (re/24.0_wp)*(0.513_wp/(1._wp + 483._wp/re**(0.669_wp)))
            if (mp <= rmacr) then
                rcd_std = (1._wp + 0.15_wp*re**(0.687_wp)) + (re/24.0_wp)*(0.42_wp/(1._wp + 42500._wp/re**(1.16_wp)))
                rmach_rat = mp/rmacr
                rcd1 = rcd_std + (rcd_mcr - rcd_std)*rmach_rat
            else if (mp <= 1.0_wp) then
                rcd_M1 = (1.0_wp + 0.118_wp*re**0.813_wp) + (re/24.0_wp)*0.69_wp/(1.0_wp + 3550.0_wp/re**0.793_wp)
                C1 = 6.48_wp
                C2 = 9.28_wp
                C3 = 12.21_wp
                f1M = -1.884_wp + 8.422_wp*mp - 13.70_wp*mp**2 + 8.162_wp*mp**3
                f2M = -2.228_wp + 10.35_wp*mp - 16.96_wp*mp**2 + 9.840_wp*mp**3
                f3M = 4.362_wp - 16.91_wp*mp + 19.84_wp*mp**2 - 6.296_wp*mp**3
                lrep = log(re)
                factor = f1M*(lrep - C2)*(lrep - C3)/((C1 - C2)*(C1 - C3)) + f2M*(lrep - C1)*(lrep - C3)/((C2 - C1)*(C2 - C3)) &
                              & + f3M*(lrep - C1)*(lrep - C2)/((C3 - C1)*(C3 - C2))
                rcd1 = rcd_mcr + (rcd_M1 - rcd_mcr)*factor
            else if (mp < 1.75_wp) then
                rcd_M1 = (1.0_wp + 0.118_wp*re**0.813_wp) + (re/24.0_wp)*0.69_wp/(1.0_wp + 3550.0_wp/re**0.793_wp)
                rcd_M2 = (1.0_wp + 0.107_wp*re**0.867_wp) + (re/24.0_wp)*0.646_wp/(1.0_wp + 861.0_wp/re**0.634_wp)
                C1 = 6.48_wp
                C2 = 8.93_wp
                C3 = 12.21_wp
                f1M = -2.963_wp + 4.392_wp*mp - 1.169_wp*mp**2 - 0.027_wp*mp**3 - 0.233_wp*exp((1.0_wp - mp)/0.011_wp)
                f2M = -6.617_wp + 12.11_wp*mp - 6.501_wp*mp**2 + 1.182_wp*mp**3 - 0.174_wp*exp((1.0_wp - mp)/0.010_wp)
                f3M = -5.866_wp + 11.57_wp*mp - 6.665_wp*mp**2 + 1.312_wp*mp**3 - 0.350_wp*exp((1.0_wp - mp)/0.012_wp)
                lrep = log(re)
                factor = f1M*(lrep - C2)*(lrep - C3)/((C1 - C2)*(C1 - C3)) + f2M*(lrep - C1)*(lrep - C3)/((C2 - C1)*(C2 - C3)) &
                              & + f3M*(lrep - C1)*(lrep - C2)/((C3 - C1)*(C3 - C2))
                rcd1 = rcd_M1 + (rcd_M2 - rcd_M1)*factor
            else
                rcd1 = (1.0_wp + 0.107_wp*re**0.867_wp) + (re/24.0_wp)*0.646_wp/(1.0_wp + 861.0_wp/re**0.634_wp)
            end if  ! mp
        end if  ! re

        ! Osnes's volume fraction correction
        b1 = 5.81_wp*phi/((1.0_wp - phi)**2) + 0.48_wp*(phi**(1._wp/3._wp))/((1.0_wp - phi)**3)

        b2 = ((1.0_wp - phi)**2)*(phi**3)*re*(0.95_wp + 0.61_wp*(phi**3)/((1.0_wp - phi)**2))

        b3 = min(sqrt(20.0_wp*mp), &
                 & 1.0_wp)*(5.65_wp*phi - 22.0_wp*(phi**2) + 23.4_wp*(phi**3))*(1._wp + tanh((mp - (0.65_wp - 0.24_wp*phi)) &
                 & /0.35_wp))

        cd = (24.0_wp/re)*rcd1

        cd = cd/(1.0_wp - phi) + b3 + (24.0_wp/re)*(1.0_wp - phi)*(b1 + b2)

        beta = 3.0_wp*pi*mu_fluid*dp*(re/24.0_wp)*cd

    end function QS_ModifiedParmar

    ! QS Force calculated as a function of Re, Ma and phi
    !
    ! Use Osnes etal (2023) correlations A.N. Osnes, M. Vartdal, M. Khalloufi, J. Capecelatro, and S. Balachandar. Comprehensive
    ! quasi-steady force correlations for compressible flow through random particle suspensions. International Journal of Multiphase
    ! Flow, Vol. 165, 104485, (2023). doi: https://doi.org/10.1016/j.imultiphaseflow.2023.104485.
    !
    ! E. Loth, J.T. Daspit, M. Jeong, T. Nagata, and T. Nonomura. Supersonic and hypersonic drag coefficients for a sphere. AIAA
    ! Journal, Vol. 59(8), pp. 3261-3274, (2021). doi: https://doi.org/10.2514/1.J060153.
    !
    ! NOTE: Re<45 Rarefied formula of Loth et al has been redefined by Balachandar to avoid singularity as Ma -> 0.
    function QS_Osnes(rho, cson, mu_fluid, gamma, vmag, dp, volume_fraction) result(beta)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: rho, cson, mu_fluid, gamma, vmag, dp, volume_fraction
        real(wp)             :: rmachp, mp, phi, rep, re
        real(wp)             :: Knp, fKn, CD1, s, JM, CD2, cd_loth, CM, GM, HM, b1, b2, b3, cd, sgby2, JMt
        real(wp)             :: beta

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
            if (mp <= 1._wp) then
                ! JMt = 2.26_wp*(mp**4) - 0.1_wp*(mp**3) + 0.14_wp*mp
                JMt = 2.26_wp*(mp**4) + 0.14_wp*mp
            else
                JMt = 1.6_wp*(mp**4) + 0.25_wp*(mp**3) + 0.11_wp*(mp**2) + 0.44_wp*mp
            end if
            !
            ! Reformulated version of Loth et al. to avoid singularity at mp = 0
            !
            CD2 = (1.0_wp + 2.0_wp*(s**2))*exp(-s**2)*mp/((sgby2**3)*sqrt(pi)) + (4.0_wp*(s**4) + 4.0_wp*(s**2) - 1.0_wp)*erf(s) &
                   & /(2.0_wp*(sgby2**4)) + (2.0_wp*(mp**3)/(3.0_wp*sgby2))*sqrt(pi)

            CD2 = CD2/(1.0_wp + (((CD2/JMt) - 1.0_wp)*sqrt(re/45.0_wp)))
            cd_loth = CD1/(1.0_wp + (mp**4)) + CD2/(1.0_wp + (mp**4))
        else
            ! Compression-dominated regime TLJ: coefficients tweaked to get continuous values on the two branches at the critical
            ! points
            if (mp < 1.5_wp) then
                CM = 1.65_wp + 0.65_wp*tanh(4._wp*mp - 3.4_wp)
            else
                ! CM = 2.18_wp - 0.13_wp*tanh(0.9_wp*mp - 2.7_wp)
                CM = 2.18_wp - 0.12913149918318745_wp*tanh(0.9_wp*mp - 2.7_wp)
            end if
            if (mp < 0.8_wp) then
                GM = 166.0_wp*(mp**3) + 3.29_wp*(mp**2) - 10.9_wp*mp + 20._wp
            else
                ! GM = 5.0_wp + 40._wp*(mp**(-3))
                GM = 5.0_wp + 47.809331200000017_wp*(mp**(-3))
            end if
            if (mp < 1._wp) then
                HM = 0.0239_wp*(mp**3) + 0.212_wp*(mp**2) - 0.074_wp*mp + 1._wp
            else
                ! HM =   0.93_wp + 1.0_wp / (3.5_wp + (mp**5))
                HM = 0.93967777777777772_wp + 1.0_wp/(3.5_wp + (mp**5))
            end if

            cd_loth = (24.0_wp/re)*(1._wp + 0.15_wp*(re**(0.687_wp)))*HM + 0.42_wp*CM/(1._wp + 42500._wp/re**(1.16_wp*CM) &
                       & + GM/sqrt(re))
        end if

        b1 = 5.81_wp*phi/((1.0_wp - phi)**2) + 0.48_wp*(phi**(1._wp/3._wp))/((1.0_wp - phi)**3)

        b2 = ((1.0_wp - phi)**2)*(phi**3)*re*(0.95_wp + 0.61_wp*(phi**3)/((1.0_wp - phi)**2))

        b3 = min(sqrt(20.0_wp*mp), &
                 & 1.0_wp)*(5.65_wp*phi - 22.0_wp*(phi**2) + 23.4_wp*(phi**3))*(1._wp + tanh((mp - (0.65_wp - 0.24_wp*phi)) &
                 & /0.35_wp))

        cd = cd_loth/(1.0_wp - phi) + b3 + (24.0_wp/re)*(1.0_wp - phi)*(b1 + b2)

        beta = 3.0_wp*pi*mu_fluid*dp*(re/24.0_wp)*cd

    end function QS_Osnes

    ! Subroutine for Quasi-Steady Drag Model of Gidaspow
    !
    ! D. Gidaspow, Multiphase Flow and Fluidization (Academic Press, 1994)
    !
    ! Note: Model is provided per cell volume. We convert that to per particle using the particle volume fraction and volume
    function QS_Gidaspow(rho, cson, mu_fluid, gamma, vmag, dp, volume_fraction) result(beta)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: rho, cson, mu_fluid, gamma, vmag, dp, volume_fraction
        real(wp)             :: cd, phifRep, phif
        real(wp)             :: phi, rep, re
        real(wp)             :: beta

        rep = vmag*dp*rho/mu_fluid
        phi = max(volume_fraction, 0.0001_wp)
        phif = max(1._wp - volume_fraction, 0.0001_wp)
        re = max(rep, 0.1_wp)

        phifRep = phif*re

        if (phifRep < 1000.0_wp) then
            cd = 24.0_wp/phifRep*(1.0_wp + 0.15_wp*(phifRep)**0.687_wp)
        else
            cd = 0.44_wp
        end if

        if (phif < 0.8_wp) then
            beta = 150.0_wp*((phi**2)*mu_fluid)/(phif*dp**2) + 1.75_wp*(rho*phi*vmag/dp)
        else
            beta = 0.75_wp*cd*phi*rho*vmag/(dp*phif**1.65_wp)
        end if

        beta = beta*(pi*dp**3)/(6.0_wp*phi)

    end function QS_Gidaspow

end module m_particles_EL_kernels
