!>
!! @file m_particles_EL_kernels.fpp
!! @brief Contains module m_particles_EL_kernels

#:include 'macros.fpp'

!> @brief This module contains kernel functions used to map the effect of the lagrangian particles in the Eulerian framework.
module m_particles_EL_kernels

    use m_mpi_proxy       !< Message passing interface (MPI) module proxy
    use ieee_arithmetic   !< For checking NaN
    use m_helper          !< For s_prng_splitmix32
    use m_euler_lagrange  !< s_get_char_vol

    implicit none

    ! Indices of the projected particle fields (q_particles)
    integer, parameter :: alphaf_id = 1
    integer, parameter :: alphaupx_id = 2   !< x particle momentum index
    integer, parameter :: alphaupy_id = 3   !< y particle momentum index
    integer, parameter :: alphaupz_id = 4   !< z particle momentum index
    integer, parameter :: alphaup2x_id = 5  !< x particle velocity squared index
    integer, parameter :: alphaup2y_id = 6  !< y particle velocity squared index
    integer, parameter :: alphaup2z_id = 7  !< z particle velocity squared index
    integer, parameter :: Smx_id = 8
    integer, parameter :: Smy_id = 9
    integer, parameter :: Smz_id = 10
    integer, parameter :: SE_id = 11

    ! Indices of the cell fields used for the particle forces and sources (field_vars)
    integer, parameter :: dPx_id = 1           !< Spatial pressure gradient in x, y, and z
    integer, parameter :: dPy_id = 2
    integer, parameter :: dPz_id = 3
    integer, parameter :: drhox_id = 4         !< Spatial density gradient in x, y, and z
    integer, parameter :: drhoy_id = 5
    integer, parameter :: drhoz_id = 6
    integer, parameter :: dufxdx_id = 7  ! du_x/dx
    integer, parameter :: dufxdy_id = 8  ! du_x/dy
    integer, parameter :: dufxdz_id = 9  ! du_x/dz
    integer, parameter :: dufydx_id = 10  ! du_y/dx
    integer, parameter :: dufydy_id = 11  ! du_y/dy
    integer, parameter :: dufydz_id = 12  ! du_y/dz
    integer, parameter :: dufzdx_id = 13  ! du_z/dx
    integer, parameter :: dufzdy_id = 14  ! du_z/dy
    integer, parameter :: dufzdz_id = 15  ! du_z/dz
    integer, parameter :: dalphafx_id = 16     !< Spatial fluid volume fraction gradient in x, y, and z
    integer, parameter :: dalphafy_id = 17
    integer, parameter :: dalphafz_id = 18
    integer, parameter :: dalphap_upx_id = 19  !< Spatial particle momentum gradient in x, y, and z
    integer, parameter :: dalphap_upy_id = 20
    integer, parameter :: dalphap_upz_id = 21
    integer, parameter :: src_tmp_id = 22      !< Scratch for p u_l and the derivatives in the pressure source terms
    integer, parameter :: dsrc_tmp_id = 23
    integer, parameter :: nField_vars = 23

    ! du_i/dx_j is field dufxdx_id + 3*(i - 1) + j - 1 (arithmetic, not a parameter array, which device code would need declared)

    integer, parameter  :: Ncells_proj = 3                    !< Cells per direction the Gaussian kernel projects onto
    integer, parameter  :: seed_kind = selected_int_kind(18)  !< 64-bit state of s_prng_splitmix32
    real(wp), parameter :: slip_speed_min = 1.e-8_wp          !< Slip speed below which the fluctuation direction is undefined
    !> Floor on basis-vector norms (and axis-alignment test) in the fluctuation model
    real(wp), parameter :: basis_norm_min = 1.e-8_wp
    real(wp), parameter :: tiny_positive = 1.e-30_wp         !< Keeps divisions and log() finite for vanishing arguments
    real(wp), parameter :: node_coincidence_tol = 1.e-10_wp  !< Distance, in cell widths, at which a particle sits on a node
    real(wp), parameter :: mach_min = 1.e-6_wp               !< Mach floor for Loth's O(M) rarefied terms, which are 0/0 at M = 0
    !> Volume-fraction cap in the radial distribution, below its singularity at 0.64356
    real(wp), parameter :: phi_chi_max = 0.64_wp
    !> Force terms reported when a particle force is not finite, indexed by the force_status of s_get_particle_force
    character(len=*), parameter :: force_term_names(4) = [character(len=24)::'quasi-steady drag','pressure gradient', &
              & 'added mass', 'drag fluctuation']
    integer  :: mapCells_loc
    real(wp) :: alpha
    $:GPU_DECLARE(create='[mapCells_loc, alpha]')

contains

    logical function f_is_finite_gpu(val) result(is_finite)

        $:GPU_ROUTINE(function_name='f_is_finite_gpu', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in) :: val

        is_finite = val == val .and. abs(val) <= huge(val)

    end function f_is_finite_gpu

    !> Set the Gaussian kernel half-width mapCells_loc = (Ncells_proj - 1)/2 and its decay rate, chosen so the kernel falls to 1e-4
    !! one cell past the half-width.
    subroutine s_initialize_particle_kernels()

        ! M = (N-1)/2, alpha chosen so G decays by 1e-4 at M+1 cells away
        mapCells_loc = (Ncells_proj - 1)/2
        alpha = 4._wp*log(10._wp)/real(mapCells_loc + 1, wp)**2._wp
        $:GPU_UPDATE(device='[mapCells_loc, alpha]')

    end subroutine s_initialize_particle_kernels

    !> Sum of one particle's Gaussian kernel over its stencil, func_s = sum(G V), used to normalize the projection.
    subroutine s_compute_gaussian_contribution(pos, cell, func_s)

        $:GPU_ROUTINE(function_name='s_compute_gaussian_contribution',parallelism='[seq]', cray_inline=True)

        real(wp), intent(in), dimension(3) :: pos
        integer, intent(in), dimension(3)  :: cell
        real(wp), intent(out)              :: func_s
        real(wp)                           :: Vol_loc, func
        real(wp), dimension(3)             :: nodecoord, center
        integer                            :: ip, jp, kp, di, dj, dk, di_beg, di_end, dj_beg, dj_end, dk_beg, dk_end
        integer, dimension(3)              :: cellijk

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

                    call s_get_char_vol(cellijk(1), cellijk(2), cellijk(3), particle_params%charwidth, Vol_loc)

                    call s_applygaussian_aniso(center, cellijk, nodecoord, func)

                    func_s = func_s + (func*Vol_loc)
                end do
            end do
        end do

    end subroutine s_compute_gaussian_contribution

    !> Add one particle's contribution to the projected fields ind_start:ind_end with atomic updates.
    subroutine s_gaussian_atomic(rad, vel, pos, force_p, gauSum, cell, updatedvar, ind_start, ind_end)

        $:GPU_ROUTINE(function_name='s_gaussian_atomic',parallelism='[seq]', cray_inline=True)

        real(wp), intent(in) :: rad, gauSum
        real(wp), intent(in), dimension(3) :: pos, vel, force_p
        integer, intent(in), dimension(3) :: cell
        type(scalar_field), dimension(:), intent(inout) :: updatedvar
        integer, intent(in) :: ind_start, ind_end
        real(wp) :: volpart, Vol_loc, func, weight
        real(wp) :: fp_x, fp_y, fp_z, vp_x, vp_y, vp_z
        real(wp) :: addFun
        real(wp), dimension(3) :: nodecoord, center
        integer :: ip, jp, kp, di, dj, dk, di_beg, di_end, dj_beg, dj_end, dk_beg, dk_end, field_ind
        integer, dimension(3) :: cellijk

        volpart = (4._wp/3._wp)*pi*rad**3._wp

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

                    call s_get_char_vol(cellijk(1), cellijk(2), cellijk(3), particle_params%charwidth, Vol_loc)

                    call s_applygaussian_aniso(center, cellijk, nodecoord, func)

                    weight = func/gauSum

                    do field_ind = ind_start, ind_end
                        if (field_ind == alphaf_id) then
                            addFun = weight*volpart
                        else if (field_ind == alphaupx_id) then
                            addFun = weight*volpart*vp_x
                        else if (field_ind == alphaupy_id) then
                            addFun = weight*volpart*vp_y
                        else if (field_ind == alphaupz_id) then
                            addFun = weight*volpart*vp_z
                        else if (field_ind == alphaup2x_id) then
                            addFun = weight*volpart*vp_x**2
                        else if (field_ind == alphaup2y_id) then
                            addFun = weight*volpart*vp_y**2
                        else if (field_ind == alphaup2z_id) then
                            addFun = weight*volpart*vp_z**2
                        else if (field_ind == Smx_id) then
                            addFun = weight*fp_x
                        else if (field_ind == Smy_id) then
                            addFun = weight*fp_y
                        else if (field_ind == Smz_id) then
                            addFun = weight*fp_z
                        else if (field_ind == SE_id) then
                            ! Work done on the particle leaves the gas (fp is the force on the gas): -F.u_p, so drag
                            ! dissipation stays in the gas as heat
                            addFun = weight*(fp_x*vp_x + fp_y*vp_y + fp_z*vp_z)
                        end if

                        $:GPU_ATOMIC(atomic='update')
                        updatedvar(field_ind)%sf(cellijk(1), cellijk(2), cellijk(3)) = updatedvar(field_ind)%sf(cellijk(1), &
                                   & cellijk(2), cellijk(3)) + real(addFun, kind=stp)
                    end do
                end do
            end do
        end do

    end subroutine s_gaussian_atomic

    !> Gaussian kernel at nodecoord for a particle at center, with the distance in each direction scaled by the local cell width.
    subroutine s_applygaussian_aniso(center, cellaux, nodecoord, func)

        $:GPU_ROUTINE(function_name='s_applygaussian_aniso',parallelism='[seq]', cray_inline=True)

        real(wp), dimension(3), intent(in) :: center
        integer, dimension(3), intent(in)  :: cellaux
        real(wp), dimension(3), intent(in) :: nodecoord
        real(wp), intent(out)              :: func
        real(wp)                           :: arg

        arg = alpha*(((center(1) - nodecoord(1))/dx(cellaux(1)))**2._wp + ((center(2) - nodecoord(2))/dy(cellaux(2)))**2._wp)

        if (num_dims == 3) then
            arg = arg + alpha*((center(3) - nodecoord(3))/dz(cellaux(3)))**2._wp
        end if

        func = exp(-arg)

    end subroutine s_applygaussian_aniso

    !> True when the first num_dims components of v are finite.
    function f_finite_vec(v) result(finite)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), dimension(3), intent(in) :: v
        logical                            :: finite
        integer                            :: dir

        finite = .true.
        do dir = 1, num_dims
            finite = finite .and. f_is_finite_gpu(v(dir))
        end do

    end function f_finite_vec

    !> Fluid velocity, density and pressure at the particle by barycentric interpolation.
    subroutine s_interp_fluid_properties(pos, cell, q_prim_vf, wx, wy, wz, fluid_vel, fluid_rho, fluid_pres)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), dimension(3), intent(in)                  :: pos
        integer, dimension(3), intent(in)                   :: cell
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), dimension(:), intent(in)        :: wx, wy, wz
        real(wp), dimension(3), intent(out)                 :: fluid_vel
        real(wp), intent(out)                               :: fluid_rho, fluid_pres
        integer                                             :: dir, l

        fluid_rho = 0._wp
        fluid_vel = 0._wp  ! z stays 0 in 2D; used in 3-component products

        do dir = 1, num_dims
            fluid_vel(dir) = f_interp_barycentric(pos, cell, q_prim_vf, eqn_idx%mom%beg + dir - 1, wx, wy, wz)
        end do

        do l = 1, num_fluids
            fluid_rho = fluid_rho + f_interp_barycentric(pos, cell, q_prim_vf, l, wx, wy, wz)
        end do
        fluid_pres = f_interp_barycentric(pos, cell, q_prim_vf, eqn_idx%E, wx, wy, wz)

    end subroutine s_interp_fluid_properties

    !> Force on one particle from quasi-steady drag, pressure gradient, added mass and drag fluctuations, and the added mass
    !! rmass_add that joins the particle inertia. The fluctuating force advances by a full time step only when advance_fluct (first
    !! RK stage); the other stages reuse it.
    subroutine s_get_particle_force(pos, rad, vel_p, Re, gamm, seed, fqsfluct, advance_fluct, cell, q_particles, fieldvars, &
                                    & rhs_old, wx, wy, wz, force, rmass_add, new_seed, new_fqsfluct, fluid_vel, fluid_rho, cson, &
                                    & force_status, re_p, mach_p)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in)                         :: rad, Re, gamm, fluid_rho, cson
        real(wp), dimension(3), intent(in)           :: pos
        integer, dimension(3), intent(in)            :: cell
        integer(seed_kind), intent(in)               :: seed
        real(wp), dimension(3), intent(in)           :: vel_p, fqsfluct, fluid_vel
        logical, intent(in)                          :: advance_fluct
        type(scalar_field), dimension(:), intent(in) :: q_particles
        type(scalar_field), dimension(:), intent(in) :: fieldvars
        type(scalar_field), dimension(:), intent(in) :: wx, wy, wz
        type(scalar_field), dimension(:), intent(in) :: rhs_old
        real(wp), dimension(3), intent(out)          :: force, new_fqsfluct
        real(wp), intent(out)                        :: rmass_add
        integer(seed_kind), intent(out)              :: new_seed
        integer, intent(out)                         :: force_status  !< 0, or the first non-finite term (see force_term_names)
        real(wp), intent(out)                        :: re_p, mach_p  !< Particle Reynolds and Mach numbers, for diagnostics
        integer(seed_kind)                           :: seed_loc
        real(wp)                                     :: vol, alpha_f
        real(wp), dimension(3)                       :: v_rel, dp
        real(wp)                                     :: particle_diam, gas_mu, vmag
        real(wp)                                     :: slip_velocity_x, slip_velocity_y, slip_velocity_z, beta
        real(wp)                                     :: vol_frac
        integer                                      :: dir, l

        ! Added pass params
        real(wp)               :: mach, Cam, SDrho, vrel_gradrho, drhodt
        real(wp), dimension(3) :: rhoDuDt, grad_rho, fam, udot_gradu

        ! QS Fluct
        real(wp), dimension(3) :: vel_p_mean, vel2_p_mean

        force = 0._wp
        dp = 0._wp
        grad_rho = 0._wp
        fam = 0._wp
        v_rel = 0._wp
        rhoDuDt = 0._wp
        SDrho = 0._wp
        vmag = 0._wp

        udot_gradu = 0._wp

        vel_p_mean = 0._wp
        vel2_p_mean = 0._wp

        ! Interpolate the projected particle fields and the gradients to the particle
        alpha_f = f_interp_barycentric(pos, cell, q_particles, alphaf_id, wx, wy, wz)
        vol_frac = 1._wp - alpha_f

        do dir = 1, num_dims
            vel_p_mean(dir) = f_interp_barycentric(pos, cell, q_particles, alphaupx_id + dir - 1, wx, wy, wz)/max(vol_frac, &
                       & verysmall)
            vel2_p_mean(dir) = f_interp_barycentric(pos, cell, q_particles, alphaup2x_id + dir - 1, wx, wy, wz)/max(vol_frac, &
                        & verysmall)
        end do

        if (particle_params%added_mass_force > 0) then
            drhodt = 0._wp
            do l = eqn_idx%cont%beg, eqn_idx%cont%end
                drhodt = drhodt + rhs_old(l)%sf(cell(1), cell(2), cell(3))  ! mixture density rate
            end do
        end if

        do dir = 1, num_dims
            if (particle_params%pressure_gradient_force .or. particle_params%added_mass_force > 0) then
                dp(dir) = f_interp_barycentric(pos, cell, fieldvars, dPx_id + dir - 1, wx, wy, wz)
            end if
            if (particle_params%added_mass_force > 0) then
                grad_rho(dir) = f_interp_barycentric(pos, cell, fieldvars, drhox_id + dir - 1, wx, wy, wz)
                rhoDuDt(dir) = (rhs_old(eqn_idx%mom%beg + dir - 1)%sf(cell(1), cell(2), cell(3)) - fluid_vel(dir)*drhodt)/fluid_rho
                do l = 1, num_dims
                    udot_gradu(dir) = udot_gradu(dir) + fluid_vel(l)*f_interp_barycentric(pos, cell, fieldvars, &
                               & dufxdx_id + 3*(dir - 1) + l - 1, wx, wy, wz)
                end do
            end if
        end do

        v_rel = vel_p - fluid_vel

        if (particle_params%qs_force > 0 .or. particle_params%added_mass_force > 0) then
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

            gas_mu = Re
        end if

        if (particle_params%added_mass_force > 0) then
            rhoDuDt = fluid_rho*(rhoDuDt + udot_gradu)
            vrel_gradrho = dot_product(-v_rel, grad_rho)
            SDrho = drhodt + vel_p(1)*grad_rho(1) + vel_p(2)*grad_rho(2) + vel_p(3)*grad_rho(3)
            mach = vmag/cson
        end if

        ! Step 1: Force component quasi-steady (zero at zero slip, where the correlations are singular)
        if (vmag > 0._wp) then
            if (particle_params%qs_force == 1) then
                beta = QS_Gidaspow(fluid_rho, cson, gas_mu, gamm, vmag, particle_diam, vol_frac)
                force = force - beta*v_rel
            else if (particle_params%qs_force == 2) then
                beta = QS_Parmar(fluid_rho, cson, gas_mu, gamm, vmag, particle_diam, vol_frac)
                force = force - beta*v_rel
            else if (particle_params%qs_force == 3) then
                beta = QS_Osnes(fluid_rho, cson, gas_mu, gamm, vmag, particle_diam, vol_frac)
                force = force - beta*v_rel
            else
                ! No Quasi-Steady drag
            end if
        end if
        force_status = 0
        if (.not. f_finite_vec(force)) force_status = 1

        ! Step 2: Pressure Gradient Force
        if (particle_params%pressure_gradient_force) then
            vol = (4._wp/3._wp)*pi*(rad**3._wp)
            force = force - vol*dp
            if (force_status == 0 .and. .not. f_finite_vec(force)) force_status = 2
        end if

        ! Step 4: Added Mass Force
        if (particle_params%added_mass_force == 1) then
            vol = (4._wp/3._wp)*pi*(rad**3._wp)
            if (mach > 0.6_wp) then
                Cam = 1._wp + 1.8_wp*(0.6_wp**2) + 7.6_wp*(0.6_wp**4)
            else
                Cam = 1._wp + 1.8_wp*mach**2 + 7.6_wp*mach**4
            end if

            Cam = 0.5_wp*Cam*(1._wp + 0.68_wp*vol_frac**2)
            rmass_add = fluid_rho*vol*Cam

            fam = Cam*vol*(-v_rel*SDrho + rhoDuDt + fluid_vel*(vrel_gradrho))
            force = force + fam
            if (force_status == 0 .and. (.not. f_finite_vec(force) .or. .not. f_is_finite_gpu(rmass_add))) force_status = 3
        else
            rmass_add = 0._wp
        end if

        if (particle_params%qs_fluct_force) then
            ! Step 5: quasi-steady force fluctuations, advanced once per time step
            seed_loc = seed
            if (advance_fluct) then
                call s_compute_qs_fluctuations(vel_p, fluid_vel, fluid_rho, cson, gas_mu, particle_diam, vol_frac, vmag, &
                                               & vel_p_mean, vel2_p_mean, seed_loc, fqsfluct, dt, new_fqsfluct)
            else
                new_fqsfluct = fqsfluct
            end if
            force = force + new_fqsfluct
            new_seed = seed_loc
            if (force_status == 0 .and. .not. f_finite_vec(force)) force_status = 4
        end if

        re_p = 0._wp
        mach_p = 0._wp
        if (vmag > 0._wp) then
            mach_p = vmag/cson
            if (gas_mu > 0._wp) re_p = fluid_rho*vmag*particle_diam/gas_mu
        end if

    end subroutine s_get_particle_force

    !> Stochastic quasi-steady force fluctuations: an Ornstein-Uhlenbeck process for the drag and lift fluctuations, driven by the
    !! granular temperature of the projected particle velocity moments (fluctuations assumed uncorrelated), advanced over dt_loc
    !! with its exact discretization. A. N. Osnes, M. Vartdal, M. Khalloufi, J. Capecelatro, and S. Balachandar, Comprehensive
    !! quasi-steady force correlations for compressible flow through random particle suspensions, Int. J. Multiphase Flow 165,
    !! 104485 (2023). A. M. Lattanzi, V. Tavanashad, S. Subramaniam, and J. Capecelatro, Stochastic model for the hydrodynamic force
    !! in Euler- Lagrange simulations of particle-laden flows, Phys. Rev. Fluids 7, 014301 (2022).
    subroutine s_compute_qs_fluctuations(vel_p, fluid_vel, fluid_rho, cson, gas_mu, particle_diam, vol_frac, vmag, vel_p_mean, &
                                         & vel2_p_mean, seed, fqs_fluct_old, dt_loc, fqs_fluct_new)
        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), dimension(3), intent(in)  :: vel_p, fluid_vel
        real(wp), intent(in)                :: fluid_rho, cson, gas_mu, particle_diam, vol_frac, vmag
        real(wp), dimension(3), intent(in)  :: vel_p_mean, vel2_p_mean
        integer(seed_kind), intent(inout)   :: seed
        real(wp), dimension(3), intent(in)  :: fqs_fluct_old
        real(wp), intent(in)                :: dt_loc
        real(wp), dimension(3), intent(out) :: fqs_fluct_new
        real(wp)                            :: upmean, vpmean, wpmean
        real(wp)                            :: u2pmean, v2pmean, w2pmean
        real(wp)                            :: rphip, rep, rmachp, theta, chi, tF_inv, phi_chi
        real(wp)                            :: fq, bq, Fs, sigD, sigT
        real(wp)                            :: decay, noise_amp
        real(wp)                            :: CD_prime, CD_frac, sigmoid_cf, f_CF
        real(wp)                            :: Z1, Z2, cosrand, sinrand
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
        rep = fluid_rho*vmag*particle_diam/gas_mu

        ! Particle Mach number
        rmachp = vmag/cson

        ! Granular temperature from Eulerian fields \theta = (<u_p^2> - <u_p>^2) / 3
        theta = ((u2pmean + v2pmean + w2pmean) - (upmean**2 + vpmean**2 + wpmean**2))/3._wp

        if (theta <= verysmall) theta = 0._wp

        ! Fluctuating drag magnitude (Osnes Eqs 9-12)
        fq = 6.52_wp*rphip - 22.56_wp*(rphip**2) + 49.90_wp*(rphip**3)
        Fs = 3._wp*pi*gas_mu*particle_diam*(1._wp + 0.15_wp*((rep*(1._wp - rphip))**0.687_wp))*(1._wp - rphip)*vmag
        bq = min(sqrt(20._wp*rmachp), 1._wp)*0.55_wp*(rphip**0.7_wp)*(1._wp + tanh((rmachp - 0.5_wp)/0.2_wp))
        sigD = (fq + bq)*Fs

        ! Granular kinetic theory relaxation rate, with the radial distribution function of D. Ma and G. Ahmadi, An equation of
        ! state for dense rigid sphere gases, J. Chem. Phys. 84, 3449 (1986); phi is capped below its singularity at 0.64356
        phi_chi = min(rphip, phi_chi_max)
        chi = (1._wp + 2.50_wp*phi_chi + 4.5904_wp*(phi_chi**2) + 4.515439_wp*(phi_chi**3))/((1._wp - (phi_chi/0.64356_wp)**3) &
               & **0.678021_wp)
        tF_inv = (24._wp*rphip*chi/particle_diam)*sqrt(theta/pi)

        ! Unit slip velocity direction
        if (vmag > slip_speed_min) then
            avec = slip_vel/vmag

            ! Project old fluctuating force onto slip direction
            CD_prime = fqs_fluct_old(1)*avec(1) + fqs_fluct_old(2)*avec(2) + fqs_fluct_old(3)*avec(3)
            CD_frac = CD_prime/max(sigD, tiny_positive)
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

        ! Build orthogonal basis
        eunit = [1._wp, 0._wp, 0._wp]
        if (abs(avec(2)) + abs(avec(3)) <= basis_norm_min) then
            eunit = [0._wp, 1._wp, 0._wp]
        else if (abs(avec(1)) + abs(avec(3)) <= basis_norm_min) then
            eunit = [0._wp, 0._wp, 1._wp]
        end if

        ! bvec = avec x eunit
        bvec(1) = avec(2)*eunit(3) - avec(3)*eunit(2)
        bvec(2) = avec(3)*eunit(1) - avec(1)*eunit(3)
        bvec(3) = avec(1)*eunit(2) - avec(2)*eunit(1)
        denum = max(basis_norm_min, sqrt(bvec(1)**2 + bvec(2)**2 + bvec(3)**2))
        bvec = bvec/denum

        ! cvec = avec x bvec
        cvec(1) = avec(2)*bvec(3) - avec(3)*bvec(2)
        cvec(2) = avec(3)*bvec(1) - avec(1)*bvec(3)
        cvec(3) = avec(1)*bvec(2) - avec(2)*bvec(1)
        denum = max(basis_norm_min, sqrt(cvec(1)**2 + cvec(2)**2 + cvec(3)**2))
        cvec = cvec/denum

        ! Generate random numbers
        call s_prng_splitmix32(UnifRnd(1), seed)
        call s_prng_splitmix32(UnifRnd(2), seed)
        call s_prng_splitmix32(UnifRnd(3), seed)
        call s_prng_splitmix32(UnifRnd(4), seed)
        call s_prng_splitmix32(UnifRnd(5), seed)

        UnifRnd(1) = max(UnifRnd(1), tiny_positive)
        UnifRnd(3) = max(UnifRnd(3), tiny_positive)

        ! Box-Muller transform
        Z1 = sqrt(-2._wp*log(UnifRnd(1)))*cos(TwoPi*UnifRnd(2))
        Z2 = sqrt(-2._wp*log(UnifRnd(3)))*cos(TwoPi*UnifRnd(4))

        ! Random perpendicular direction
        cosrand = cos(TwoPi*UnifRnd(5))
        sinrand = sin(TwoPi*UnifRnd(5))
        dvec = bvec*cosrand + cvec*sinrand
        denum = max(basis_norm_min, sqrt(dvec(1)**2 + dvec(2)**2 + dvec(3)**2))
        dvec = dvec/denum

        ! Exact Ornstein-Uhlenbeck update over dt_loc with relaxation rate tF_inv: the stationary standard deviations are sigD
        ! (along the slip) and sigT (perpendicular); stable for any tF_inv*dt_loc, and equal to Euler-Maruyama as it goes to 0
        decay = exp(-tF_inv*dt_loc)
        noise_amp = sqrt(max(1._wp - decay**2, 0._wp))
        fqs_fluct_new(1) = decay*fqs_fluct_old(1) + noise_amp*(sigD*Z1*avec(1) + sigT*Z2*dvec(1))
        fqs_fluct_new(2) = decay*fqs_fluct_old(2) + noise_amp*(sigD*Z1*avec(2) + sigT*Z2*dvec(2))
        if (num_dims == 3) then
            fqs_fluct_new(3) = decay*fqs_fluct_old(3) + noise_amp*(sigD*Z1*avec(3) + sigT*Z2*dvec(3))
        end if

    end subroutine s_compute_qs_fluctuations

    !> Barycentric Lagrange interpolation of field_vf(field_index) to pos with the precomputed weights, limited to the range of the
    !! stencil values so it cannot overshoot near discontinuities (e.g. a negative density behind a shock).
    function f_interp_barycentric(pos, cell, field_vf, field_index, wx, wy, wz) result(val)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), dimension(3), intent(in)           :: pos
        integer, dimension(3), intent(in)            :: cell
        type(scalar_field), dimension(:), intent(in) :: field_vf
        type(scalar_field), dimension(:), intent(in) :: wx, wy, wz
        integer, intent(in)                          :: field_index
        integer                                      :: i, j, k, ix, jy, kz, npts, npts_z, N
        integer                                      :: ix_count, jy_count, kz_count
        integer                                      :: hit_x, hit_y, hit_z  ! stencil index of a coincident node, 0 = none
        real(wp)                                     :: fx, fy, fz, weight, numerator, denominator, val, tol
        real(wp)                                     :: f_node, f_min, f_max

        i = cell(1); j = cell(2); k = cell(3)

        N = particle_params%interpolation_order
        npts = N/2
        npts_z = npts
        if (num_dims == 2) npts_z = 0

        ! Detect on-node coincidence per direction (relative to local spacing)
        hit_x = 0; hit_y = 0; hit_z = 0
        ix_count = 0
        do ix = i - npts, i + npts
            ix_count = ix_count + 1
            tol = node_coincidence_tol*dx(ix)
            if (abs(pos(1) - x_cc(ix)) <= tol) hit_x = ix_count
        end do
        jy_count = 0
        do jy = j - npts, j + npts
            jy_count = jy_count + 1
            tol = node_coincidence_tol*dy(jy)
            if (abs(pos(2) - y_cc(jy)) <= tol) hit_y = jy_count
        end do
        if (num_dims == 3) then
            kz_count = 0
            do kz = k - npts_z, k + npts_z
                kz_count = kz_count + 1
                tol = node_coincidence_tol*dz(kz)
                if (abs(pos(3) - z_cc(kz)) <= tol) hit_z = kz_count
            end do
        end if

        numerator = 0._wp
        denominator = 0._wp
        f_min = huge(1._wp)
        f_max = -huge(1._wp)

        ix_count = 0
        do ix = i - npts, i + npts
            ix_count = ix_count + 1
            if (hit_x /= 0 .and. ix_count /= hit_x) cycle
            if (hit_x /= 0) then
                fx = 1._wp
            else
                fx = wx(ix_count)%sf(i, 1, 1)/(pos(1) - x_cc(ix))
            end if

            jy_count = 0
            do jy = j - npts, j + npts
                jy_count = jy_count + 1
                if (hit_y /= 0 .and. jy_count /= hit_y) cycle
                if (hit_y /= 0) then
                    fy = 1._wp
                else
                    fy = wy(jy_count)%sf(j, 1, 1)/(pos(2) - y_cc(jy))
                end if

                kz_count = 0
                do kz = k - npts_z, k + npts_z
                    kz_count = kz_count + 1
                    if (num_dims == 3) then
                        if (hit_z /= 0 .and. kz_count /= hit_z) cycle
                        if (hit_z /= 0) then
                            fz = 1._wp
                        else
                            fz = wz(kz_count)%sf(k, 1, 1)/(pos(3) - z_cc(kz))
                        end if
                    else
                        fz = 1._wp
                    end if

                    weight = fx*fy*fz
                    f_node = field_vf(field_index)%sf(ix, jy, kz)
                    numerator = numerator + weight*f_node
                    denominator = denominator + weight
                    f_min = min(f_min, f_node)
                    f_max = max(f_max, f_node)
                end do
            end do
        end do

        ! Polynomial interpolation overshoots across a jump; keep the value within the stencil range. A NaN in the stencil is
        ! left to propagate so the particle force check reports it.
        val = numerator/denominator
        if (f_is_finite_gpu(val)) val = max(f_min, min(f_max, val))  ! min/max may drop a NaN, so clamp finite values only

        if (abs(val) <= verysmall) val = 0._wp

    end function f_interp_barycentric

    !> Quasi-steady drag coefficient beta (force = -beta*(u_p - u_f)) with Reynolds- and Mach-number corrections, and a volume-
    !! fraction correction for dilute random arrays. M. Parmar, A. Haselbacher, and S. Balachandar, Improved drag correlation for
    !! spheres and application to shock-tube experiments, AIAA J. 48(6), 1273-1276 (2010). A. S. Sangani, D. Z. Zhang, and A.
    !! Prosperetti, The added mass, Basset, and viscous drag coefficients in nondilute bubbly liquids undergoing small-amplitude
    !! oscillatory motion, Phys. Fluids A 3, 2955 (1991).
    function QS_Parmar(rho, cson, mu_fluid, gamma, vmag, dp, volume_fraction) result(beta)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: rho, cson, mu_fluid, gamma, vmag, dp, volume_fraction
        real(wp)             :: rcd1, rmacr, rcd_mcr, rcd_std, rmach_rat, rcd_M1
        real(wp)             :: rcd_M2, C1, C2, C3, f1M, f2M, f3M, lrep, factor, cd, phi_corr
        real(wp)             :: beta
        real(wp)             :: rmachp, phi, rep, re

        rmachp = vmag/cson
        phi = volume_fraction
        rep = vmag*dp*rho/mu_fluid
        re = max(rep, 0.1_wp)  ! bounds the log(Re) fits below, anchored at ln Re = 6.5 to 12.2

        rmacr = 0.6_wp  ! Critical rmachp no
        rcd_mcr = (1._wp + 0.15_wp*re**(0.684_wp)) + (re/24.0_wp)*(0.513_wp/(1._wp + 483._wp/re**(0.669_wp)))
        if (rmachp <= rmacr) then
            rcd_std = (1._wp + 0.15_wp*re**(0.687_wp)) + (re/24.0_wp)*(0.42_wp/(1._wp + 42500._wp/re**(1.16_wp)))
            rmach_rat = rmachp/rmacr
            rcd1 = rcd_std + (rcd_mcr - rcd_std)*rmach_rat
        else if (rmachp <= 1.0_wp) then
            rcd_M1 = (1.0_wp + 0.118_wp*re**0.813_wp) + (re/24.0_wp)*0.69_wp/(1.0_wp + 3550.0_wp/re**0.793_wp)
            C1 = 6.48_wp
            C2 = 9.28_wp
            C3 = 12.21_wp
            f1M = -1.884_wp + 8.422_wp*rmachp - 13.70_wp*rmachp**2 + 8.162_wp*rmachp**3
            f2M = -2.228_wp + 10.35_wp*rmachp - 16.96_wp*rmachp**2 + 9.840_wp*rmachp**3
            f3M = 4.362_wp - 16.91_wp*rmachp + 19.84_wp*rmachp**2 - 6.296_wp*rmachp**3
            lrep = log(re)
            factor = f1M*(lrep - C2)*(lrep - C3)/((C1 - C2)*(C1 - C3)) + f2M*(lrep - C1)*(lrep - C3)/((C2 - C1)*(C2 - C3)) &
                          & + f3M*(lrep - C1)*(lrep - C2)/((C3 - C1)*(C3 - C2))
            rcd1 = rcd_mcr + (rcd_M1 - rcd_mcr)*factor
        else if (rmachp < 1.75_wp) then
            rcd_M1 = (1.0_wp + 0.118_wp*re**0.813_wp) + (re/24.0_wp)*0.69_wp/(1.0_wp + 3550.0_wp/re**0.793_wp)
            rcd_M2 = (1.0_wp + 0.107_wp*re**0.867_wp) + (re/24.0_wp)*0.646_wp/(1.0_wp + 861.0_wp/re**0.634_wp)
            C1 = 6.48_wp
            C2 = 8.93_wp
            C3 = 12.21_wp
            f1M = -2.963_wp + 4.392_wp*rmachp - 1.169_wp*rmachp**2 - 0.027_wp*rmachp**3 - 0.233_wp*exp((1.0_wp - rmachp)/0.011_wp)
            f2M = -6.617_wp + 12.11_wp*rmachp - 6.501_wp*rmachp**2 + 1.182_wp*rmachp**3 - 0.174_wp*exp((1.0_wp - rmachp)/0.010_wp)
            f3M = -5.866_wp + 11.57_wp*rmachp - 6.665_wp*rmachp**2 + 1.312_wp*rmachp**3 - 0.350_wp*exp((1.0_wp - rmachp)/0.012_wp)
            lrep = log(re)
            factor = f1M*(lrep - C2)*(lrep - C3)/((C1 - C2)*(C1 - C3)) + f2M*(lrep - C1)*(lrep - C3)/((C2 - C1)*(C2 - C3)) &
                          & + f3M*(lrep - C1)*(lrep - C2)/((C3 - C1)*(C3 - C2))
            rcd1 = rcd_M1 + (rcd_M2 - rcd_M1)*factor
        else
            rcd1 = (1.0_wp + 0.107_wp*re**0.867_wp) + (re/24.0_wp)*0.646_wp/(1.0_wp + 861.0_wp/re**0.634_wp)
        end if  ! rmachp

        ! Sangani's volume fraction correction for dilute random arrays Capping volume fraction at 0.5
        phi_corr = (1.0_wp + 5.94_wp*min(phi, 0.5_wp))

        cd = (24.0_wp/re)*rcd1*phi_corr

        beta = rcd1*3.0_wp*pi*mu_fluid*dp

        beta = beta*phi_corr

    end function QS_Parmar

    !> Quasi-steady drag coefficient beta as a function of Re, Ma and volume fraction. A. N. Osnes, M. Vartdal, M. Khalloufi, J.
    !! Capecelatro, and S. Balachandar, Comprehensive quasi-steady force correlations for compressible flow through random particle
    !! suspensions, Int. J. Multiphase Flow 165, 104485 (2023). E. Loth, J. T. Daspit, M. Jeong, T. Nagata, and T. Nonomura,
    !! Supersonic and hypersonic drag coefficients for a sphere, AIAA J. 59(8), 3261-3274 (2021). The rarefied (Re < 45) branch is
    !! reformulated to avoid the singularity as Ma -> 0, and the compression branch uses coefficients that keep J_M, C_M, G_M and
    !! H_M continuous across their switch points, as given in the Appendix of T. Daoud, T. Jackson, and S. Balachandar, A careful
    !! examination of closure models in Euler-Lagrange simulations of compressible multiphase flow in a planar shock particle
    !! curtain problem, Int. J. Multiphase Flow (2026), doi:10.1016/j.ijmultiphaseflow.2026.105740.
    function QS_Osnes(rho, cson, mu_fluid, gamma, vmag, dp, volume_fraction) result(beta)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: rho, cson, mu_fluid, gamma, vmag, dp, volume_fraction
        real(wp)             :: rmachp, mp, phi, re
        real(wp)             :: Knp, fKn, CD1, s, CD2, cd_loth, CM, GM, HM, b1, b2, b3, cd, sgby2, JMt
        real(wp)             :: beta

        rmachp = vmag/cson
        mp = max(rmachp, mach_min)  ! Loth's rarefied terms CD2 and JMt are both O(M)
        phi = volume_fraction
        ! No Re floor: vmag > 0 here, and every 24/re term in cd is multiplied by re/24 in beta (Stokes limit)
        re = vmag*dp*rho/mu_fluid

        ! Loth's correlation
        if (re <= 45.0_wp) then
            ! Rarefied-dominated regime
            Knp = sqrt(0.5_wp*pi*gamma)*rmachp/re
            if (Knp > 0.01_wp) then
                fKn = 1.0_wp/(1.0_wp + Knp*(2.514_wp + 0.8_wp*exp(-0.55_wp/Knp)))
            else
                fKn = 1.0_wp/(1.0_wp + Knp*(2.514_wp + 0.8_wp*exp(-0.55_wp/0.01_wp)))
            end if
            CD1 = (24.0_wp/re)*(1.0_wp + 0.15_wp*re**(0.687_wp))*fKn
            s = mp*sqrt(0.5_wp*gamma)
            sgby2 = sqrt(0.5_wp*gamma)
            ! J_M continuous at M = 1 (Daoud et al. 2026, Appendix)
            if (mp <= 1._wp) then
                JMt = 2.26_wp*(mp**4) + 0.14_wp*mp
            else
                JMt = 1.6_wp*(mp**4) + 0.25_wp*(mp**3) + 0.11_wp*(mp**2) + 0.44_wp*mp
            end if
            ! Rarefied term reformulated to avoid the singularity at M = 0 (Daoud et al. 2026, Appendix)
            CD2 = (1.0_wp + 2.0_wp*(s**2))*exp(-s**2)*mp/((sgby2**3)*sqrt(pi)) + (4.0_wp*(s**4) + 4.0_wp*(s**2) - 1.0_wp)*erf(s) &
                   & /(2.0_wp*(sgby2**4)) + (2.0_wp*(mp**3)/(3.0_wp*sgby2))*sqrt(pi)

            CD2 = CD2/(1.0_wp + (((CD2/JMt) - 1.0_wp)*sqrt(re/45.0_wp)))
            cd_loth = CD1/(1.0_wp + (mp**4)) + CD2/(1.0_wp + (mp**4))
        else
            ! Compression-dominated regime, with coefficients that keep C_M, G_M and H_M continuous at M = 1.5, 0.8 and 1
            ! (Daoud et al. 2026, Appendix)
            if (mp < 1.5_wp) then
                CM = 1.65_wp + 0.65_wp*tanh(4._wp*mp - 3.4_wp)
            else
                CM = 2.18_wp - 0.12913149918318745_wp*tanh(0.9_wp*mp - 2.7_wp)
            end if
            if (mp < 0.8_wp) then
                GM = 166.0_wp*(mp**3) + 3.29_wp*(mp**2) - 10.9_wp*mp + 20._wp
            else
                GM = 5.0_wp + 47.809331200000017_wp*(mp**(-3))
            end if
            if (mp < 1._wp) then
                HM = 0.0239_wp*(mp**3) + 0.212_wp*(mp**2) - 0.074_wp*mp + 1._wp
            else
                HM = 0.93967777777777772_wp + 1.0_wp/(3.5_wp + (mp**5))
            end if

            cd_loth = (24.0_wp/re)*(1._wp + 0.15_wp*(re**(0.687_wp)))*HM + 0.42_wp*CM/(1._wp + 42500._wp/re**(1.16_wp*CM) &
                       & + GM/sqrt(re))
        end if

        b1 = 5.81_wp*phi/((1.0_wp - phi)**2) + 0.48_wp*(phi**(1._wp/3._wp))/((1.0_wp - phi)**3)

        b2 = ((1.0_wp - phi)**2)*(phi**3)*re*(0.95_wp + 0.61_wp*(phi**3)/((1.0_wp - phi)**2))

        b3 = min(sqrt(20.0_wp*rmachp), &
                 & 1.0_wp)*(5.65_wp*phi - 22.0_wp*(phi**2) + 23.4_wp*(phi**3))*(1._wp + tanh((rmachp - (0.65_wp - 0.24_wp*phi)) &
                 & /0.35_wp))

        cd = cd_loth/(1.0_wp - phi) + b3 + (24.0_wp/re)*(1.0_wp - phi)*(b1 + b2)

        beta = 3.0_wp*pi*mu_fluid*dp*(re/24.0_wp)*cd

    end function QS_Osnes

    !> Quasi-steady drag coefficient beta from the Gidaspow model, converted from per cell volume to per particle with the particle
    !! volume fraction and volume. D. Gidaspow, Multiphase Flow and Fluidization, Academic Press (1994).
    function QS_Gidaspow(rho, cson, mu_fluid, gamma, vmag, dp, volume_fraction) result(beta)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: rho, cson, mu_fluid, gamma, vmag, dp, volume_fraction
        real(wp)             :: cd, phifRep, phif
        real(wp)             :: phi, rep
        real(wp)             :: beta

        rep = vmag*dp*rho/mu_fluid
        phi = volume_fraction
        phif = max(1._wp - volume_fraction, 0.0001_wp)

        ! No Re floor: cd*vmag must keep its Stokes limit (vmag > 0 here, so phifRep > 0)
        phifRep = phif*rep

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

        beta = beta*(pi*dp**3)/(6.0_wp*(phi + verysmall))

    end function QS_Gidaspow

end module m_particles_EL_kernels
