!>
!! @file
!! @brief Contains module m_body_forces

#:include 'macros.fpp'

#! dV = volume of cell (j, k, l)
#:def CELL_DV()
    dV = dx(j)
    if (n > 0) dV = dV*dy(k)
    if (p > 0) dV = dV*dz(l)
#:enddef

!> @brief Computes gravitational and body force source terms for the momentum equations
module m_body_forces

    use m_derived_types
    use m_global_parameters
    use m_variables_conversion
    use m_mpi_proxy
    use m_helper, only: s_prng_splitmix32, f_unit_vector, f_cross
    use m_helper_basic, only: f_is_default
    use m_thermochem, only: num_species, gas_constant, get_mixture_molecular_weight, get_mixture_specific_heat_cp_mass
    use m_nvtx

    implicit none

    private
    public :: s_compute_body_forces_rhs, s_compute_synthetic_forces_rhs, s_initialize_body_forces_module, &
        & s_finalize_body_forces_module, s_update_mean_sources

    logical, dimension(3)                   :: mpg_dir    !< Mean pressure gradient active in this direction
    real(wp), dimension(3)                  :: mflux_old  !< Previous-step bulk mass flux, for the controller
    integer                                 :: mpg_unit   !< mean_dpdx.dat history file unit
    logical                                 :: mts_on     !< Mean gas heat source (thermostat) active
    real(wp)                                :: temp_old   !< Previous-step bulk temperature, for the controller
    integer                                 :: mts_unit   !< mean_temp.dat history file unit
    integer, parameter                      :: spbf_num_freq = 8
    real(wp)                                :: spbf_amp
    real(wp)                                :: spbf_xc
    real(wp)                                :: spbf_yc
    real(wp)                                :: spbf_conv_vel
    real(wp)                                :: spbf_sigma
    real(wp), allocatable, dimension(:)     :: freq, phase
    real(wp), allocatable, dimension(:,:,:) :: rhoM
    $:GPU_DECLARE(create='[spbf_amp, spbf_xc, spbf_yc, spbf_conv_vel, spbf_sigma, freq, phase, rhoM]')

    integer :: num_synthetic_wave_numbers
    $:GPU_DECLARE(create='[num_synthetic_wave_numbers]')

    real(wp), allocatable, dimension(:) :: synthetic_k_x, synthetic_k_y, synthetic_k_z
    real(wp), allocatable, dimension(:) :: synthetic_phase, synthetic_amp
    ! Solenoidal (divergence-free) unit vector per mode: perpendicular to k in the plane of excitation.
    ! 1-D: ex = 1, ey = 0  (streamwise only)
    ! 2-D: ex = -k_y/|k|, ey = k_x/|k|  (perpendicular to k on the circle)
    ! 3-D: stored from random frame built at init
    real(wp), allocatable, dimension(:) :: synthetic_ex, synthetic_ey, synthetic_ez
    $:GPU_DECLARE(create='[synthetic_k_x, synthetic_k_y, synthetic_k_z, synthetic_phase, synthetic_amp]')
    $:GPU_DECLARE(create='[synthetic_ex, synthetic_ey, synthetic_ez]')

contains

    !> Initialize the body forces module. When synthetic_turbulence is enabled,
    !> generates random wave vectors for each energy shell on rank 0, then
    !> broadcasts to all MPI ranks and copies to GPU.
    impure subroutine s_initialize_body_forces_module

        integer                             :: s, m_wave, m_global
        integer(kind=selected_int_kind(18)) :: seed
        real(wp)                            :: rn1, rn2, k_mag
        real(wp), dimension(3)              :: khat, xi, sig, sig_tmp

        if (n > 0) then
            if (p > 0) then
                @:ALLOCATE(rhoM(-buff_size:buff_size + m, -buff_size:buff_size + n, -buff_size:buff_size + p))
            else
                @:ALLOCATE(rhoM(-buff_size:buff_size + m, -buff_size:buff_size + n, 0:0))
            end if
        else
            @:ALLOCATE(rhoM(-buff_size:buff_size + m, 0:0, 0:0))
        end if

        if (bf_spatial_support) then
            call s_initialize_body_force_with_spatial_support
        end if

        mpg_dir = mean_dpdx /= 0._wp .or. .not. f_is_default(mean_mflux)
        mts_on = mean_qdot /= 0._wp .or. .not. f_is_default(mean_temp)
        mflux_old = dflt_real
        temp_old = dflt_real
        if (any(mpg_dir)) call s_open_history(mpg_unit, 'mean_dpdx.dat')
        if (mts_on) call s_open_history(mts_unit, 'mean_temp.dat')

        if (.not. synthetic_turbulence) return
        if (synth_n_shells <= 0) return

        num_synthetic_wave_numbers = sum(synth_n_waves_per_shell(1:synth_n_shells))

        if (num_synthetic_wave_numbers == 0) return

        @:ALLOCATE(synthetic_k_x(1:num_synthetic_wave_numbers))
        @:ALLOCATE(synthetic_phase(1:num_synthetic_wave_numbers))
        @:ALLOCATE(synthetic_amp(1:num_synthetic_wave_numbers))
        @:ALLOCATE(synthetic_ex(1:num_synthetic_wave_numbers))
        @:ALLOCATE(synthetic_ey(1:num_synthetic_wave_numbers))
        @:ALLOCATE(synthetic_ez(1:num_synthetic_wave_numbers))

        if (num_dims > 1) then
            @:ALLOCATE(synthetic_k_y(1:num_synthetic_wave_numbers))
        end if

        if (num_dims == 3) then
            @:ALLOCATE(synthetic_k_z(1:num_synthetic_wave_numbers))
        end if

        ! Generate random wave vectors and phases on rank 0, then broadcast. Uses the
        ! compiler-independent integer generator s_prng_splitmix32 so the forcing is
        ! reproducible across compilers; the 3-D polarization is built perpendicular
        ! to k via a double cross product, guaranteeing a divergence-free (solenoidal) mode.
        if (proc_rank == 0) then
            seed = int(synth_seed, kind(seed))

            m_global = 0
            do s = 1, synth_n_shells
                k_mag = synth_k_shell(s)
                do m_wave = 1, synth_n_waves_per_shell(s)
                    m_global = m_global + 1

                    if (num_dims == 1) then
                        synthetic_k_x(m_global) = k_mag
                        ! Streamwise forcing only in 1-D
                        synthetic_ex(m_global) = 1._wp
                        synthetic_ey(m_global) = 0._wp
                        synthetic_ez(m_global) = 0._wp
                    else if (num_dims == 2) then
                        ! In-plane wavevector at azimuth theta; solenoidal dir perpendicular to k
                        call s_prng_splitmix32(rn1, seed)
                        rn1 = rn1*2._wp*pi
                        synthetic_k_x(m_global) = k_mag*cos(rn1)
                        synthetic_k_y(m_global) = k_mag*sin(rn1)
                        synthetic_ex(m_global) = -sin(rn1)
                        synthetic_ey(m_global) = cos(rn1)
                        synthetic_ez(m_global) = 0._wp
                    else
                        ! Random unit wavevector uniform on the sphere
                        call s_prng_splitmix32(rn1, seed)
                        call s_prng_splitmix32(rn2, seed)
                        khat = f_unit_vector(rn1, rn2)
                        ! Random reference vector projected perpendicular to k by a double
                        ! cross product: sig = khat x (xi x khat) is a unit vector with k.sig = 0
                        call s_prng_splitmix32(rn1, seed)
                        call s_prng_splitmix32(rn2, seed)
                        xi = f_unit_vector(rn1, rn2)
                        sig_tmp = f_cross(xi, khat)
                        sig_tmp = sig_tmp/max(sqrt(sum(sig_tmp**2._wp)), 1.e-10_wp)
                        sig = f_cross(khat, sig_tmp)
                        synthetic_k_x(m_global) = k_mag*khat(1)
                        synthetic_k_y(m_global) = k_mag*khat(2)
                        synthetic_k_z(m_global) = k_mag*khat(3)
                        synthetic_ex(m_global) = sig(1)
                        synthetic_ey(m_global) = sig(2)
                        synthetic_ez(m_global) = sig(3)
                    end if

                    call s_prng_splitmix32(rn1, seed)
                    synthetic_phase(m_global) = rn1*2._wp*pi

                    synthetic_amp(m_global) = synth_amp_shell(s)
                end do
            end do
        end if

        ! Broadcast from rank 0 to all ranks
        call s_mpi_send_random_number(synthetic_k_x, num_synthetic_wave_numbers)
        call s_mpi_send_random_number(synthetic_phase, num_synthetic_wave_numbers)
        call s_mpi_send_random_number(synthetic_amp, num_synthetic_wave_numbers)
        call s_mpi_send_random_number(synthetic_ex, num_synthetic_wave_numbers)
        call s_mpi_send_random_number(synthetic_ey, num_synthetic_wave_numbers)
        call s_mpi_send_random_number(synthetic_ez, num_synthetic_wave_numbers)

        if (num_dims > 1) then
            call s_mpi_send_random_number(synthetic_k_y, num_synthetic_wave_numbers)
        end if

        if (num_dims == 3) then
            call s_mpi_send_random_number(synthetic_k_z, num_synthetic_wave_numbers)
        end if

        ! Push wave data to GPU
        $:GPU_UPDATE(device='[num_synthetic_wave_numbers, synthetic_k_x, synthetic_phase, synthetic_amp]')
        $:GPU_UPDATE(device='[synthetic_ex, synthetic_ey, synthetic_ez]')

        if (num_dims > 1) then
            $:GPU_UPDATE(device='[synthetic_k_y]')
        end if

        if (num_dims == 3) then
            $:GPU_UPDATE(device='[synthetic_k_z]')
        end if

    end subroutine s_initialize_body_forces_module

    !> Initialize a body force with spatial support presented in Wei & Freund (JFM, 2005)
    impure subroutine s_initialize_body_force_with_spatial_support

        integer :: f  !< frequency iterator

        spbf_amp = spatial_bf%amp
        spbf_xc = spatial_bf%x_centroid
        spbf_yc = spatial_bf%y_centroid
        spbf_conv_vel = spatial_bf%conv_vel
        spbf_sigma = spatial_bf%sigma

        @:ALLOCATE(freq(spbf_num_freq), phase(spbf_num_freq))
        @:PREFER_GPU(freq)
        @:PREFER_GPU(phase)
        do f = 1, spbf_num_freq
            freq(f) = spatial_bf%freq(f)
            phase(f) = spatial_bf%phase(f)
        end do
        $:GPU_UPDATE(device='[spbf_amp, spbf_xc, spbf_yc, spbf_conv_vel, spbf_sigma, freq, phase]')

    end subroutine s_initialize_body_force_with_spatial_support

    !> Compute the acceleration at time t
    subroutine s_compute_acceleration(t)

        real(wp), intent(in) :: t

        #:for DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
            if (bf_${XYZ}$) then
                accel_bf(${DIR}$) = g_${XYZ}$ + k_${XYZ}$*sin(w_${XYZ}$*t - p_${XYZ}$)
            end if
        #:endfor

        $:GPU_UPDATE(device='[accel_bf]')

    end subroutine s_compute_acceleration

    !> Apply the body force of Wei & Freund (JFM, 2005)
    subroutine s_compute_body_force_with_spatial_support(t, bounds)

        real(wp), intent(in)                              :: t
        type(int_bounds_info), dimension(1:3), intent(in) :: bounds
        real(wp)                                          :: support                    !< spatial support
        real(wp)                                          :: theta_x, theta_y, pre_fac  !< auxiliary variables
        integer                                           :: f                          !< frequency iterator
        integer                                           :: j, k, l                    !< standard iterators

        ! Safety check: if convective velocity is too small, skip computation

        if (abs(spbf_conv_vel) < 1.0e-12_wp) then
            $:GPU_PARALLEL_LOOP(private='[j, k, l]', collapse=3, copyin='[bounds]')
            do l = bounds(3)%beg, bounds(3)%end
                do k = bounds(2)%beg, bounds(2)%end
                    do j = bounds(1)%beg, bounds(1)%end
                        spbf_source_x(j, k, l) = 0._wp
                        spbf_source_y(j, k, l) = 0._wp
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            return
        end if

        $:GPU_PARALLEL_LOOP(private='[support, theta_x, theta_y, pre_fac, f, j, k, l]', collapse=3, copyin='[bounds, t]')
        do l = bounds(3)%beg, bounds(3)%end
            do k = bounds(2)%beg, bounds(2)%end
                do j = bounds(1)%beg, bounds(1)%end
                    support = exp(-spbf_sigma*((x_cc(j) - spbf_xc)**2 + (y_cc(k) - spbf_yc)**2))
                    spbf_source_x(j, k, l) = 0._wp
                    spbf_source_y(j, k, l) = 0._wp
                    do f = 1, spbf_num_freq
                        pre_fac = (freq(f)/spbf_conv_vel)
                        theta_x = pre_fac*(x_cc(j) - spbf_xc - spbf_conv_vel*t) + phase(f)
                        theta_y = pre_fac*(y_cc(k) - spbf_yc) + phase(f)
                        spbf_source_x(j, k, l) = spbf_source_x(j, k, &
                                      & l) + spbf_amp*support*(pre_fac*sin(theta_x)*cos(theta_y) - 2._wp*spbf_sigma*(y_cc(k) &
                                      & - spbf_yc)*sin(theta_x)*sin(theta_y))
                        spbf_source_y(j, k, l) = spbf_source_y(j, k, &
                                      & l) - spbf_amp*support*(pre_fac*cos(theta_x)*sin(theta_y) - 2._wp*spbf_sigma*(x_cc(j) &
                                      & - spbf_xc)*sin(theta_x)*sin(theta_y))
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_body_force_with_spatial_support

    !> Compute the mixture density at each cell center
    subroutine s_compute_mixture_density(q_cons_vf, bounds)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        type(int_bounds_info), dimension(1:3), intent(in)   :: bounds
        integer                                             :: i, j, k, l  !< standard iterators

        $:GPU_PARALLEL_LOOP(private='[j, k, l]', collapse=3, copyin='[bounds]')
        do l = bounds(3)%beg, bounds(3)%end
            do k = bounds(2)%beg, bounds(2)%end
                do j = bounds(1)%beg, bounds(1)%end
                    rhoM(j, k, l) = 0._wp
                    do i = 1, num_fluids
                        rhoM(j, k, l) = rhoM(j, k, l) + q_cons_vf(eqn_idx%cont%beg + i - 1)%sf(j, k, l)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_mixture_density

    !> Compute the body force source terms for momentum and energy equations
    subroutine s_compute_body_forces_rhs(q_prim_vf, q_cons_vf, rhs_vf, bounds)

        type(scalar_field), dimension(sys_size), intent(in)    :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(in)    :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        type(int_bounds_info), dimension(1:3), intent(in)      :: bounds
        integer                                                :: i, j, k, l  !< Loop variables
        real(wp)                                               :: f_mpg

        if (bf_x .or. bf_y .or. bf_z) then
            call s_compute_acceleration(mytime)
        end if

        if (bf_spatial_support) then
            call s_compute_body_force_with_spatial_support(mytime, bounds)
        end if

        call s_compute_mixture_density(q_cons_vf, bounds)

        $:GPU_PARALLEL_LOOP(private='[i, j, k, l]', collapse=4, copyin='[bounds]')
        do i = eqn_idx%mom%beg, eqn_idx%E
            do l = bounds(3)%beg, bounds(3)%end
                do k = bounds(2)%beg, bounds(2)%end
                    do j = bounds(1)%beg, bounds(1)%end
                        rhs_vf(i)%sf(j, k, l) = 0._wp
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        if (bf_spatial_support) then
            $:GPU_PARALLEL_LOOP(private='[j, k, l]', collapse=3, copyin='[bounds]')
            do l = bounds(3)%beg, bounds(3)%end
                do k = bounds(2)%beg, bounds(2)%end
                    do j = bounds(1)%beg, bounds(1)%end
                        rhs_vf(eqn_idx%mom%beg)%sf(j, k, l) = rhs_vf(eqn_idx%mom%beg)%sf(j, k, l) + spbf_source_x(j, k, l)
                        rhs_vf(eqn_idx%mom%beg + 1)%sf(j, k, l) = rhs_vf(eqn_idx%mom%beg + 1)%sf(j, k, l) + spbf_source_y(j, k, l)
                        ! Energy work term u*f: velocity (mom/rho) dotted with the momentum source,
                        ! matching the bf_x/y/z convention below so the forcing is energy-consistent.
                        rhs_vf(eqn_idx%E)%sf(j, k, l) = rhs_vf(eqn_idx%E)%sf(j, k, l) + (q_cons_vf(eqn_idx%mom%beg)%sf(j, k, &
                               & l)/rhoM(j, k, l))*spbf_source_x(j, k, l) + (q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k, l)/rhoM(j, &
                               & k, l))*spbf_source_y(j, k, l)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        ! bf_*: rho*g; mean pressure gradient: uniform -dP/dx in fluid cells. Energy source u*f for both.
        #:for DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
            if (bf_${XYZ}$ .or. mpg_dir(${DIR}$)) then
                $:GPU_PARALLEL_LOOP(private='[j, k, l, f_mpg]', collapse=3, copyin='[bounds]')
                do l = bounds(3)%beg, bounds(3)%end
                    do k = bounds(2)%beg, bounds(2)%end
                        do j = bounds(1)%beg, bounds(1)%end
                            f_mpg = mpg_force(${DIR}$)
                            rhs_vf(eqn_idx%mom%beg + ${DIR - 1}$)%sf(j, k, l) = rhs_vf(eqn_idx%mom%beg + ${DIR - 1}$)%sf(j, k, &
                                   & l) + rhoM(j, k, l)*accel_bf(${DIR}$) + f_mpg
                            rhs_vf(eqn_idx%E)%sf(j, k, l) = rhs_vf(eqn_idx%E)%sf(j, k, &
                                   & l) + q_cons_vf(eqn_idx%mom%beg + ${DIR - 1}$)%sf(j, k, l)*(accel_bf(${DIR}$) + f_mpg/rhoM(j, &
                                   & k, l))
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

        ! Thermostat: uniform heat source per volume
        if (mts_on) then
            $:GPU_PARALLEL_LOOP(private='[j, k, l]', collapse=3, copyin='[bounds]')
            do l = bounds(3)%beg, bounds(3)%end
                do k = bounds(2)%beg, bounds(2)%end
                    do j = bounds(1)%beg, bounds(1)%end
                        rhs_vf(eqn_idx%E)%sf(j, k, l) = rhs_vf(eqn_idx%E)%sf(j, k, l) + mts_qdot
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_compute_body_forces_rhs

    !> Rank 0 opens a per-step history file in case_dir, appending on restart
    impure subroutine s_open_history(unit, fname)

        integer, intent(out)         :: unit
        character(len=*), intent(in) :: fname

        if (proc_rank /= 0) return
        if (t_step_start == 0) then
            open (newunit=unit, file=trim(case_dir) // '/' // fname, status='replace')
        else
            open (newunit=unit, file=trim(case_dir) // '/' // fname, position='append')
        end if

    end subroutine s_open_history

    !> Integral (deadbeat for G = 1 and a constant load) update of a mean source s holding a fluid-mean quantity y at y_target: s +=
    !! (G/dt) c (y_target - 2y^n + y^(n-1)), c the source-to-y capacity. y_old carries y^(n-1) between steps.
    subroutine s_mean_controller(src, y, y_old, y_target, gain, c)

        real(wp), intent(inout) :: src, y_old
        real(wp), intent(in)    :: y, y_target, gain, c

        if (f_is_default(y_target)) return
        if (f_is_default(y_old)) y_old = y
        src = src + gain/dt*c*(y_target - 2._wp*y + y_old)
        y_old = y

    end subroutine s_mean_controller

    !> Once per step, fluid-cell averages driving the mean sources (sums over fluid cells). Mean pressure gradient: bulk mass flux
    !! Q_i = sum(rho u_i dV)/sum(dV), controlled to mean_mflux(i) by -dP/dx; logged to mean_dpdx.dat. Thermostat: mass-averaged
    !! temperature T_b = sum(rho T dV)/sum(rho dV), controlled to mean_temp by the heat source with capacity rho c_v (fluid mean);
    !! logged to mean_temp.dat with the mixing-cup T_m = sum(rho u_1 T dV)/sum(rho u_1 dV) and the fluid volume.
    impure subroutine s_update_mean_sources(q_cons_vf, q_T_sf, t_step)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        type(scalar_field), intent(in)                      :: q_T_sf
        integer, intent(in)                                 :: t_step
        real(wp), dimension(1, 8)                           :: sums_loc, sums_glb
        real(wp)                                            :: q1, q2, q3, vol, rm, rT, rcv, ruT, dV, rho, T, cp, mw
        real(wp), dimension(num_species)                    :: Ys
        real(wp), dimension(3)                              :: mflux
        real(wp)                                            :: temp_b, temp_m, rcv_b
        integer                                             :: i, j, k, l

        if (.not. (any(mpg_dir) .or. mts_on)) return

        q1 = 0._wp; q2 = 0._wp; q3 = 0._wp; vol = 0._wp
        $:GPU_PARALLEL_LOOP(private='[j, k, l, dV]', collapse=3, reduction='[[q1, q2, q3, vol]]', reductionOp='[+]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    @:CELL_DV()
                    vol = vol + dV
                    q1 = q1 + q_cons_vf(eqn_idx%mom%beg)%sf(j, k, l)*dV
                    if (n > 0) q2 = q2 + q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k, l)*dV
                    if (p > 0) q3 = q3 + q_cons_vf(eqn_idx%mom%beg + 2)%sf(j, k, l)*dV
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        rm = 0._wp; rT = 0._wp; rcv = 0._wp; ruT = 0._wp
        if (mts_on) then
            ! q_T_sf holds T at t^n here: the stage-1 RHS has just converted q_cons_vf
            $:GPU_PARALLEL_LOOP(private='[i, j, k, l, dV, rho, T, cp, mw, Ys]', collapse=3, reduction='[[rm, rT, rcv, ruT]]', &
                                & reductionOp='[+]')
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        @:CELL_DV()
                        rho = q_cons_vf(eqn_idx%cont%beg)%sf(j, k, l)
                        T = q_T_sf%sf(j, k, l)
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, num_species
                            Ys(i) = q_cons_vf(eqn_idx%species%beg + i - 1)%sf(j, k, l)/rho
                        end do
                        call get_mixture_specific_heat_cp_mass(T, Ys, cp)
                        call get_mixture_molecular_weight(Ys, mw)
                        rm = rm + rho*dV
                        rT = rT + rho*T*dV
                        rcv = rcv + rho*(cp - gas_constant/mw)*dV
                        ruT = ruT + q_cons_vf(eqn_idx%mom%beg)%sf(j, k, l)*T*dV
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        sums_loc(1,:) = [q1, q2, q3, vol, rm, rT, rcv, ruT]
        call s_mpi_allreduce_vectors_sum(sums_loc, sums_glb, 1, 8)
        vol = sums_glb(1, 4)

        if (any(mpg_dir)) then
            mflux = sums_glb(1,1:3)/vol
            do i = 1, 3
                call s_mean_controller(mpg_force(i), mflux(i), mflux_old(i), mean_mflux(i), mean_dpdx_gain, 1._wp)
            end do
            $:GPU_UPDATE(device='[mpg_force]')
            if (proc_rank == 0) write (mpg_unit, '(I0, 7(1X, ES24.16))') t_step, mytime, mflux, -mpg_force
        end if

        if (mts_on) then
            temp_b = sums_glb(1, 6)/sums_glb(1, 5)
            rcv_b = sums_glb(1, 7)/vol
            call s_mean_controller(mts_qdot, temp_b, temp_old, mean_temp, mean_temp_gain, rcv_b)
            $:GPU_UPDATE(device='[mts_qdot]')
            temp_m = temp_b
            if (sums_glb(1, 1) /= 0._wp) temp_m = sums_glb(1, 8)/sums_glb(1, 1)
            if (proc_rank == 0) write (mts_unit, '(I0, 6(1X, ES24.16))') t_step, mytime, temp_b, temp_m, mts_qdot, vol
        end if

    end subroutine s_update_mean_sources

    !> Compute the synthetic turbulence RHS contribution.
    !>
    !> Each Fourier mode m has wave vector k_m and a pre-computed solenoidal
    !> (divergence-free) direction e_m perpendicular to k_m. The force field
    !> is F = rho * (sum_m A_m * cos(k_m.(x-U*t*xhat) + phi_m) * e_m) * (U/Lx) * G.
    !> This excites all velocity components correctly and is zero outside the
    !> bounding box of the Gaussian window.
    subroutine s_compute_synthetic_forces_rhs(q_prim_vf, q_cons_vf, rhs_vf)

        type(scalar_field), dimension(sys_size), intent(in)    :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(in)    :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        integer                                                :: i, j, k, l, n_mode, turb_idx, iq
        real(wp)                                               :: pos_x, pos_y, pos_z
        real(wp)                                               :: gauss_env, G_norm, f_scale
        real(wp)                                               :: a_m, phase_arg
        real(wp)                                               :: force_x, force_y, force_z
        real(wp)                                               :: rho_local, u_local_x, u_local_y, u_local_z
        real(wp)                                               :: adv_offset
        logical                                                :: in_box

        ! G_bar for exp(-pi/2 * (2x/L)^2) over [-L/2, L/2]: erf(sqrt(pi/2))/sqrt(pi/2) ~ 0.7602
        real(wp), parameter :: G_BAR = 0.760173_wp

        ! Zero momentum and energy RHS for the synthetic turbulence pass

        $:GPU_PARALLEL_LOOP(private='[i, j, k, l]', collapse=4)
        do i = eqn_idx%mom%beg, eqn_idx%E
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        rhs_vf(i)%sf(j, k, l) = 0._wp
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! Pre-compute advection offset on CPU (mytime has no GPU declaration)
        adv_offset = synth_U_inf*mytime

        do turb_idx = 1, num_turbulent_sources
            $:GPU_PARALLEL_LOOP(private='[j, k, l, n_mode, iq, pos_x, pos_y, pos_z, gauss_env, G_norm, f_scale, a_m, phase_arg, &
                                & force_x, force_y, force_z, rho_local, u_local_x, u_local_y, u_local_z, in_box]', collapse=3)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        ! Position relative to forcing-zone center
                        pos_x = x_cc(j) - turb_pos(turb_idx, 1)
                        pos_y = 0._wp
                        pos_z = 0._wp
                        if (num_dims > 1) pos_y = y_cc(k) - turb_pos(turb_idx, 2)
                        if (num_dims == 3) pos_z = z_cc(l) - turb_pos(turb_idx, 3)

                        ! Bounding-box skip
                        in_box = abs(pos_x) <= 0.5_wp*synth_L(turb_idx, 1)
                        if (num_dims > 1) in_box = in_box .and. abs(pos_y) <= 0.5_wp*synth_L(turb_idx, 2)
                        if (num_dims == 3) in_box = in_box .and. abs(pos_z) <= 0.5_wp*synth_L(turb_idx, 3)

                        if (in_box) then
                            ! Gaussian envelope, normalised by G_BAR^num_dims
                            gauss_env = exp(-0.5_wp*pi*(2._wp*pos_x/synth_L(turb_idx, 1))**2)
                            if (num_dims > 1) gauss_env = gauss_env*exp(-0.5_wp*pi*(2._wp*pos_y/synth_L(turb_idx, 2))**2)
                            if (num_dims == 3) gauss_env = gauss_env*exp(-0.5_wp*pi*(2._wp*pos_z/synth_L(turb_idx, 3))**2)
                            G_norm = gauss_env/G_BAR**num_dims

                            ! Accumulate solenoidal force components over all modes.
                            ! Each mode contributes A_m * cos(phase) in direction e_m
                            ! (e_m is perpendicular to k_m, precomputed at init).
                            force_x = 0._wp
                            force_y = 0._wp
                            force_z = 0._wp
                            do n_mode = 1, num_synthetic_wave_numbers
                                phase_arg = synthetic_k_x(n_mode)*(x_cc(j) - adv_offset) + synthetic_phase(n_mode)
                                if (num_dims > 1) phase_arg = phase_arg + synthetic_k_y(n_mode)*y_cc(k)
                                if (num_dims == 3) phase_arg = phase_arg + synthetic_k_z(n_mode)*z_cc(l)
                                a_m = synthetic_amp(n_mode)*cos(phase_arg)
                                force_x = force_x + a_m*synthetic_ex(n_mode)
                                force_y = force_y + a_m*synthetic_ey(n_mode)
                                force_z = force_z + a_m*synthetic_ez(n_mode)
                            end do

                            ! Scale: rho * (U_inf/L_x) * G_norm
                            ! Mixture density = sum of all partial densities (continuity components)
                            rho_local = 0._wp
                            do iq = eqn_idx%cont%beg, eqn_idx%cont%end
                                rho_local = rho_local + q_prim_vf(iq)%sf(j, k, l)
                            end do
                            f_scale = rho_local*(synth_U_inf/synth_L(turb_idx, 1))*G_norm

                            force_x = force_x*f_scale
                            force_y = force_y*f_scale
                            force_z = force_z*f_scale

                            ! Local velocities for energy update F.u
                            u_local_x = q_prim_vf(eqn_idx%mom%beg)%sf(j, k, l)
                            u_local_y = 0._wp
                            u_local_z = 0._wp
                            if (num_dims > 1) u_local_y = q_prim_vf(eqn_idx%mom%beg + 1)%sf(j, k, l)
                            if (num_dims == 3) u_local_z = q_prim_vf(eqn_idx%mom%end)%sf(j, k, l)

                            rhs_vf(eqn_idx%mom%beg)%sf(j, k, l) = rhs_vf(eqn_idx%mom%beg)%sf(j, k, l) + force_x
                            if (num_dims > 1) rhs_vf(eqn_idx%mom%beg + 1)%sf(j, k, l) = rhs_vf(eqn_idx%mom%beg + 1)%sf(j, k, &
                                & l) + force_y
                            if (num_dims == 3) rhs_vf(eqn_idx%mom%end)%sf(j, k, l) = rhs_vf(eqn_idx%mom%end)%sf(j, k, l) + force_z
                            rhs_vf(eqn_idx%E)%sf(j, k, l) = rhs_vf(eqn_idx%E)%sf(j, k, &
                                   & l) + force_x*u_local_x + force_y*u_local_y + force_z*u_local_z
                        end if
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do

    end subroutine s_compute_synthetic_forces_rhs

    !> Finalize the body forces module
    impure subroutine s_finalize_body_forces_module

        @:DEALLOCATE(rhoM)

        if (proc_rank == 0) then
            if (any(mpg_dir)) close (mpg_unit)
            if (mts_on) close (mts_unit)
        end if

        if (bf_spatial_support) then
            @:DEALLOCATE(freq, phase)
        end if

        if (num_synthetic_wave_numbers > 0) then
            @:DEALLOCATE(synthetic_k_x)
            @:DEALLOCATE(synthetic_phase)
            @:DEALLOCATE(synthetic_amp)
            @:DEALLOCATE(synthetic_ex)
            @:DEALLOCATE(synthetic_ey)
            @:DEALLOCATE(synthetic_ez)
            if (num_dims > 1) then
                @:DEALLOCATE(synthetic_k_y)
            end if
            if (num_dims == 3) then
                @:DEALLOCATE(synthetic_k_z)
            end if
        end if

    end subroutine s_finalize_body_forces_module

end module m_body_forces
