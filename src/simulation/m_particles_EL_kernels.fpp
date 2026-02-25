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

contains

    !> The purpose of this subroutine is to smear the strength of the lagrangian
            !!      particles into the Eulerian framework using a gaussian approach.
            !! @param nParticles Number of lagrangian particles in the current domain
            !! @param lbk_rad Radius of the particles
            !! @param lbk_s Computational coordinates of the particles
            !! @param lbk_pos Spatial coordinates of the particles
            !! @param updatedvar Eulerian variable to be updated
            !! @param lbk_f_p Forces on the particles
    subroutine s_smoothfunction(nParticles, lbk_rad, lbk_s, lbk_pos, lbk_vel, updatedvar, lbk_f_p)

        integer, intent(in) :: nParticles
        real(wp), dimension(1:lag_params%nParticles_glb, 1:3, 1:2), intent(in) :: lbk_s, lbk_pos, lbk_vel
        real(wp), dimension(1:lag_params%nParticles_glb, 1:2), intent(in) :: lbk_rad
        type(scalar_field), dimension(:), intent(inout) :: updatedvar
        real(wp), dimension(1:lag_params%nParticles_glb, 1:3), intent(in) :: lbk_f_p

        call s_gaussian(nParticles, lbk_rad, lbk_s, lbk_pos, lbk_vel, updatedvar, lbk_f_p)

    end subroutine s_smoothfunction

    !> The purpose of this procedure contains the algorithm to use the gaussian kernel function to map the effect of the particles.
    subroutine s_gaussian(nParticles, lbk_rad, lbk_s, lbk_pos, lbk_vel, updatedvar, lbk_f_p)

        integer, intent(in) :: nParticles
        real(wp), dimension(1:lag_params%nParticles_glb, 1:3, 1:2), intent(in) :: lbk_s, lbk_pos, lbk_vel
        real(wp), dimension(1:lag_params%nParticles_glb, 1:2), intent(in) :: lbk_rad
        real(wp), dimension(1:lag_params%nParticles_glb, 1:3), intent(in) :: lbk_f_p
        type(scalar_field), dimension(:), intent(inout) :: updatedvar

        real(wp), dimension(3) :: center
        integer, dimension(3) :: cell
        integer, dimension(3) :: cellaux
        real(wp) :: stddsv

        real(wp), dimension(3) :: nodecoord
        real(wp) :: addFun1, addFun2_x, addFun2_y, addFun2_z, addFun_E, func_sum
        real(wp) :: addFun_alphap_vp_x, addFun_alphap_vp_y, addFun_alphap_vp_z
        real(wp) :: func, volpart, Vol, Vol_loc, rad
        real(wp), dimension(3) :: s_coord
        integer :: l, i, j, k
        logical :: celloutside
        real(wp) :: fp_x, fp_y, fp_z, vp_x, vp_y, vp_z

        real(wp) :: rc, r2, weight, alpha_cut, vol_frac, g, chardist
        integer :: d
        integer :: ncx, ncy, ncz, ix, jy, kz

        d = 2
        if (num_dims == 3) d = 3

        $:GPU_PARALLEL_LOOP(private='[nodecoord,l,s_coord,cell,center,volpart,rad,cellaux]', copyin='[d]')
        do l = 1, nParticles
            nodecoord(1:3) = 0
            center(1:3) = 0._wp
            volpart = 4._wp/3._wp*pi*lbk_rad(l, 2)**3._wp
            s_coord(1:3) = lbk_s(l, 1:3, 2)
            center(1:2) = lbk_pos(l, 1:2, 2)
            rad = lbk_rad(l, 2)

            if (p > 0) center(3) = lbk_pos(l, 3, 2)
            call s_get_cell(s_coord, cell)

            if (num_dims == 2) then
                Vol = dx(cell(1))*dy(cell(2))*lag_params%charwidth
                if (cyl_coord) Vol = dx(cell(1))*dy(cell(2))*y_cc(cell(2))*2._wp*pi
            else
                Vol = dx(cell(1))*dy(cell(2))*dz(cell(3))
            end if

            !< Characteristic cell length
            chardist = sqrt(dx(cell(1))*dy(cell(2)))
            if (num_dims == 3) chardist = (dx(cell(1))*dy(cell(2))*dz(cell(3)))**(1._wp/3._wp)

            vol_frac = volpart/Vol

            ! !< Update void fraction field
            ! addFun1 = vol_frac
            ! $:GPU_ATOMIC(atomic='update')
            ! updatedvar(1)%sf(cell(1), cell(2), cell(3)) = &
            !     updatedvar(1)%sf(cell(1), cell(2), cell(3)) &
            !     + real(addFun1, kind=stp)

            fp_x = -lbk_f_p(l, 1)
            fp_y = -lbk_f_p(l, 2)
            fp_z = -lbk_f_p(l, 3)

            vp_x = lbk_vel(l, 1, 2)
            vp_y = lbk_vel(l, 2, 2)
            vp_z = lbk_vel(l, 3, 2)

            sigma = lag_params%epsilonb*chardist
            ! sigma = max(sigma, 0.5_wp * chardist)

            ! alpha_cut = 1.e-4_wp
            rc = sigma
            ! rc = sigma * sqrt(-2.0_wp*log(alpha_cut))
            ncx = max(1, int(rc/dx(cell(1)))) !min(3,ceiling(rc/dx(cell(1))))
            ncy = max(1, int(rc/dx(cell(1)))) !ceiling(rc/dy(cell(2))) !min(3,ceiling(rc/dy(cell(2))))
            ncz = 0
            if (num_dims == 3) ncz = max(1, int(rc/dx(cell(1)))) !ceiling(rc/dz(cell(3))) !min(3,ceiling(rc/dz(cell(3)))) !p

            i = cell(1)
            j = cell(2)
            k = cell(3)

            func_sum = 0._wp
            $:GPU_LOOP(collapse=3,private='[ix,jy,kz,r2,cellaux]')
            do ix = i - ncx, i + ncx
                do jy = j - ncy, j + ncy
                    do kz = k - ncz, k + ncz

                        cellaux(1) = ix
                        cellaux(2) = jy
                        cellaux(3) = kz
                        ! Relocate cells for particles intersecting symmetric boundaries
                        if (any((/bc_x%beg, bc_x%end, bc_y%beg, bc_y%end, bc_z%beg, bc_z%end/) == BC_REFLECTIVE)) then
                            call s_shift_cell_symmetric_bc(cellaux, cell)
                        end if

                        Vol_loc = dx(cellaux(1))*dy(cellaux(2))
                        if (num_dims == 3) Vol_loc = Vol_loc*dz(cellaux(3))
                        r2 = (x_cc(cellaux(1)) - center(1))**2 + (y_cc(cellaux(2)) - center(2))**2
                        if (num_dims == 3) r2 = r2 + (z_cc(cellaux(3)) - center(3))**2
                        g = exp(-r2/(2*sigma**2))/(sigma*sqrt(2*pi))**d
                        func_sum = func_sum + g*Vol_loc
                    end do
                end do
            end do

            $:GPU_LOOP(collapse=3,private='[ix,jy,kz,r2,weight,g]')
            do ix = i - ncx, i + ncx
                do jy = j - ncy, j + ncy
                    do kz = k - ncz, k + ncz

                        cellaux(1) = ix
                        cellaux(2) = jy
                        cellaux(3) = kz
                        ! Relocate cells for particles intersecting symmetric boundaries
                        if (any((/bc_x%beg, bc_x%end, bc_y%beg, bc_y%end, bc_z%beg, bc_z%end/) == BC_REFLECTIVE)) then
                            call s_shift_cell_symmetric_bc(cellaux, cell)
                        end if

                        r2 = (x_cc(cellaux(1)) - center(1))**2 + (y_cc(cellaux(2)) - center(2))**2
                        if (num_dims == 3) r2 = r2 + (z_cc(cellaux(3)) - center(3))**2
                        ! Vol_loc = dx(cellaux(1))*dy(cellaux(2))
                        ! if (num_dims==3) Vol_loc = Vol_loc * dz(cellaux(3))
                        g = exp(-r2/(2*sigma**2))/(sigma*sqrt(2*pi))**d
                        weight = g/func_sum

                        !Update volume fraction field
                        addFun1 = weight*volpart
                        $:GPU_ATOMIC(atomic='update')
                        updatedvar(1)%sf(cellaux(1), cellaux(2), cellaux(3)) = &
                            updatedvar(1)%sf(cellaux(1), cellaux(2), cellaux(3)) &
                            + real(addFun1, kind=stp)

                        if (lag_params%solver_approach == 2) then

                            !Update x-momentum source term
                            addFun2_x = weight*fp_x
                            $:GPU_ATOMIC(atomic='update')
                            updatedvar(2)%sf(cellaux(1), cellaux(2), cellaux(3)) = &
                                updatedvar(2)%sf(cellaux(1), cellaux(2), cellaux(3)) &
                                + real(addFun2_x, kind=stp)

                            !Update y-momentum source term
                            addFun2_y = weight*fp_y
                            $:GPU_ATOMIC(atomic='update')
                            updatedvar(3)%sf(cellaux(1), cellaux(2), cellaux(3)) = &
                                updatedvar(3)%sf(cellaux(1), cellaux(2), cellaux(3)) &
                                + real(addFun2_y, kind=stp)

                            if (num_dims == 3) then
                                !Update z-momentum source term
                                addFun2_z = weight*fp_z
                                $:GPU_ATOMIC(atomic='update')
                                updatedvar(4)%sf(cellaux(1), cellaux(2), cellaux(3)) = &
                                    updatedvar(4)%sf(cellaux(1), cellaux(2), cellaux(3)) &
                                    + real(addFun2_z, kind=stp)
                            end if
                            !Update energy source term
                            addFun_E = 0._wp
                            $:GPU_ATOMIC(atomic='update')
                            updatedvar(5)%sf(cellaux(1), cellaux(2), cellaux(3)) = &
                                updatedvar(5)%sf(cellaux(1), cellaux(2), cellaux(3)) &
                                + real(addFun_E, kind=stp)

                            !Update particle momentum field(x)
                            addFun_alphap_vp_x = weight*volpart*vp_x
                            $:GPU_ATOMIC(atomic='update')
                            updatedvar(6)%sf(cellaux(1), cellaux(2), cellaux(3)) = &
                                updatedvar(6)%sf(cellaux(1), cellaux(2), cellaux(3)) &
                                + real(addFun_alphap_vp_x, kind=stp)
                            !Update particle momentum field(y)
                            addFun_alphap_vp_y = weight*volpart*vp_y
                            $:GPU_ATOMIC(atomic='update')
                            updatedvar(7)%sf(cellaux(1), cellaux(2), cellaux(3)) = &
                                updatedvar(7)%sf(cellaux(1), cellaux(2), cellaux(3)) &
                                + real(addFun_alphap_vp_y, kind=stp)
                            if (num_dims == 3) then
                                !Update particle momentum field(z)
                                addFun_alphap_vp_z = weight*volpart*vp_z
                                $:GPU_ATOMIC(atomic='update')
                                updatedvar(8)%sf(cellaux(1), cellaux(2), cellaux(3)) = &
                                    updatedvar(8)%sf(cellaux(1), cellaux(2), cellaux(3)) &
                                    + real(addFun_alphap_vp_z, kind=stp)
                            end if
                        end if
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_gaussian

    !> This subroutine relocates the current cell, if it intersects a symmetric boundary.
            !! @param cell Cell of the current particle
            !! @param cellaux Cell to map the particle effect in.
    subroutine s_shift_cell_symmetric_bc(cellaux, cell)
        $:GPU_ROUTINE(function_name='s_shift_cell_symmetric_bc', &
            & parallelism='[seq]', cray_inline=True)

        integer, dimension(3), intent(inout) :: cellaux
        integer, dimension(3), intent(in) :: cell

        ! x-dir
        if (bc_x%beg == BC_REFLECTIVE .and. (cell(1) <= mapCells - 1)) then
            cellaux(1) = abs(cellaux(1)) - 1
        end if
        if (bc_x%end == BC_REFLECTIVE .and. (cell(1) >= m + 1 - mapCells)) then
            cellaux(1) = cellaux(1) - (2*(cellaux(1) - m) - 1)
        end if

        !y-dir
        if (bc_y%beg == BC_REFLECTIVE .and. (cell(2) <= mapCells - 1)) then
            cellaux(2) = abs(cellaux(2)) - 1
        end if
        if (bc_y%end == BC_REFLECTIVE .and. (cell(2) >= n + 1 - mapCells)) then
            cellaux(2) = cellaux(2) - (2*(cellaux(2) - n) - 1)
        end if

        if (p > 0) then
            !z-dir
            if (bc_z%beg == BC_REFLECTIVE .and. (cell(3) <= mapCells - 1)) then
                cellaux(3) = abs(cellaux(3)) - 1
            end if
            if (bc_z%end == BC_REFLECTIVE .and. (cell(3) >= p + 1 - mapCells)) then
                cellaux(3) = cellaux(3) - (2*(cellaux(3) - p) - 1)
            end if
        end if

    end subroutine s_shift_cell_symmetric_bc

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

        do dir = 1, num_dims
            if (.not. ieee_is_finite(force(dir))) then
                force(dir) = 0._wp
            end if
        end do

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

        if (abs(pos(1) - x_cc(i)) <= eps .and. &
            abs(pos(2) - y_cc(j)) <= eps .and. &
            abs(pos(3) - z_cc(k)) <= eps) then
            val = field_vf(field_index)%sf(i, j, k)
            return
        end if

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

        ! if (num_dims == 3) then
        !     local_min = minval(field_vf(field_index)%sf(i - npts:i + npts, j - npts:j + npts, k - npts:k + npts))
        !     local_max = maxval(field_vf(field_index)%sf(i - npts:i + npts, j - npts:j + npts, k - npts:k + npts))
        ! else
        !     local_min = minval(field_vf(field_index)%sf(i - npts:i + npts, j - npts:j + npts, k))
        !     local_max = maxval(field_vf(field_index)%sf(i - npts:i + npts, j - npts:j + npts, k))
        ! end if

        ! val = max(local_min, min(local_max, val))

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
            if (Knp > 0.01) then
                fKn = 1.0_wp/(1.0_wp + Knp*(2.514_wp + 0.8_wp*exp(-0.55_wp/Knp)))
            else
                fKn = 1.0_wp/(1.0_wp + Knp*(2.514_wp + 0.8_wp*exp(-0.55_wp/0.01)))

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
                CM = 1.65_wp + 0.65_wp*tanh(4_wp*mp - 3.4_wp)
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

            cd_loth = (24.0_wp/re)*(1 + 0.15_wp*(re**(0.687)))*HM + &
                      0.42_wp*CM/(1 + 42500/re**(1.16*CM) + GM/sqrt(re))

        end if

        b1 = 5.81_wp*phi/((1.0 - phi)**2) + &
             0.48_wp*(phi**(1._wp/3._wp))/((1.0 - phi)**3)

        b2 = ((1.0_wp - phi)**2)*(phi**3)* &
             re*(0.95 + 0.61*(phi**3)/((1.0 - phi)*2))

        b3 = min(sqrt(20.0_wp*mp), 1.0_wp)* &
             (5.65*phi - 22.0*(phi**2) + 23.4*(phi**3))* &
             (1 + tanh((mp - (0.65 - 0.24*phi))/0.35))

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
