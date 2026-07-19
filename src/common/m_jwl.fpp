!>
!! @file
!! @brief JWL EOS with a linear mass-fraction two-material mixture closure for ideal-gas and stiffened ambients.

#:include 'macros.fpp'
#:include 'case.fpp'

!> @brief Jones-Wilkins-Lee (JWL) EOS and its linear mass-fraction two-material closure for the five-equation model (Allaire et al.,
!! JCP 2002).
!!
!! Pure JWL products (Mie-Grueneisen referenced to an isentrope), V = rho0/rho:
!!   p = A(1 - w/(R1 V)) exp(-R1 V) + B(1 - w/(R2 V)) exp(-R2 V) + w rho e
!! A, B [Pa], R1, R2, w = omega [-] are cylinder-test fits.
!! Refs: Lee/Hornig/Kury UCRL-50422 (1968); Menikoff LA-UR-15-29536 (2015).
!!
!! Two-material closure (products + ambient of mass fraction Y = alpha_rho_prod/rho), blended
!! linearly in Y (w = Y): An = Y A, Bn = Y B, omega = air_gamma + Y(omega0 - air_gamma),
!! cv = Y cv_j + (1-Y) cv_a. Exact at Y=0 (ambient) and Y=1 (pure JWL). A stiffened ambient
!! (pi_inf > 0) adds pi_hat = (1-Y) pi_inf; pi_inf = 0 recovers the ideal-gas closure
!! bit-identically. Every coefficient depends on Y alone, never rho or e, so the (rho, p, Y) -> e
!! inverse is closed-form and the sound speed is the exact Grueneisen derivative
!! c^2 = (dp/drho)_e + (p/rho^2)(dp/de)_rho. s_jwl_verify_closure sweeps the (rho, e, Y) envelope
!! at start-up and aborts on any non-positive sound speed or failed p<->e round trip.
!!
!! Stiffened-gas ambient: Le Metayer/Massoni/Saurel, Int. J. Therm. Sci. 43, 265 (2004).
module m_jwl

    use m_global_parameters

    implicit none

    private
    public :: s_initialize_jwl_module, s_finalize_jwl_module, s_jwl_mix_state_er, s_jwl_mix_energy_pr, s_jwl_mix_sound_speed, &
        & s_jwl_mix_energy_sound_speed_pr, jwl_idx

    ! Simulation builds use m_global_parameters tables.
#ifndef MFC_SIMULATION
    real(wp), allocatable, public, dimension(:) :: jwl_As, jwl_Bs, jwl_R1s, jwl_R2s, jwl_omegas, jwl_rho0s, jwl_E0s
    real(wp), allocatable, public, dimension(:) :: jwl_air_e0s, jwl_air_rho0s, jwl_air_gammas, jwl_air_pi_infs
    $:GPU_DECLARE(create='[jwl_As, jwl_Bs, jwl_R1s, jwl_R2s, jwl_omegas, jwl_rho0s, jwl_E0s]')
    $:GPU_DECLARE(create='[jwl_air_e0s, jwl_air_rho0s, jwl_air_gammas, jwl_air_pi_infs]')
#endif

    integer  :: jwl_idx                  !< JWL fluid index.
    real(wp) :: jwl_cv_prod, jwl_cv_air  !< Products/air specific heats.
    $:GPU_DECLARE(create='[jwl_idx]')
    $:GPU_DECLARE(create='[jwl_cv_prod, jwl_cv_air]')

contains

    !> Floor x to `floor`; NaNs pass through unchanged.
    subroutine s_jwl_floor(x, floor)

        $:GPU_ROUTINE(function_name='s_jwl_floor',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(inout) :: x
        real(wp), intent(in)    :: floor

        if (x == x) then
            if (x < floor) x = floor
        end if

    end subroutine s_jwl_floor

    !> Effective mixture coefficients (linear in Y). omega relaxes air_gamma -> omega0 as products fill the cell; a stiffened
    !! ambient adds pi_hat = (1-Y)*air_pi_inf. Y-only dependence keeps the pressure law linear in e, so the inverse and c2 are
    !! closed-form.
    subroutine s_jwl_mixture_coeffs(Y, A, B, omega0, air_gamma, air_pi_inf, cv_j, cv_a, An, Bn, omega, cv, pi_c, pi_hat)

        $:GPU_ROUTINE(function_name='s_jwl_mixture_coeffs',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: Y, A, B, omega0, air_gamma, air_pi_inf, cv_j, cv_a
        real(wp), intent(out) :: An, Bn, omega, cv, pi_c, pi_hat

        An = Y*A
        Bn = Y*B
        omega = air_gamma + Y*(omega0 - air_gamma)
        pi_hat = (1._wp - Y)*air_pi_inf
        pi_c = (air_gamma + 1._wp)*pi_hat
        cv = Y*cv_j + (1._wp - Y)*cv_a  ! affects T only; p and c are cv-free

    end subroutine s_jwl_mixture_coeffs

    !> Mixture state from energy: (rho, e, Y) -> (p, T, c2, c2 floor).
    subroutine s_jwl_state_er(rho, e, Y, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, pres, T, c2, c2_floor)

        $:GPU_ROUTINE(function_name='s_jwl_state_er',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, e, Y, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a
        real(wp), intent(out) :: pres, T, c2, c2_floor
        real(wp)              :: rho_s, Y_s, An, Bn, omega, cv, pi_c, pi_hat, V, exp1, exp2, coef1, coef2

        rho_s = max(rho, sgm_eps)
        Y_s = min(max(Y, 0._wp), 1._wp)

        call s_jwl_mixture_coeffs(Y_s, A, B, omega0, air_gamma, air_pi_inf, cv_j, cv_a, An, Bn, omega, cv, pi_c, pi_hat)

        V = rho0/rho_s
        exp1 = exp(-R1*V)
        exp2 = exp(-R2*V)
        coef1 = (1._wp - omega/(R1*V))*exp1
        coef2 = (1._wp - omega/(R2*V))*exp2

        pres = An*coef1 + Bn*coef2 + omega*rho_s*e - pi_c
        T = (pres + pi_c - An*exp1 - Bn*exp2)/(omega*cv*rho_s)

        ! Frozen sound speed c2 = dp/drho|e + (p/rho^2) dp/de|rho (pi_c drops out).
        c2 = exp1*An*(R1*rho0/rho_s**2 - omega/rho_s - omega/(R1*rho0)) + exp2*Bn*(R2*rho0/rho_s**2 - omega/rho_s &
                      & - omega/(R2*rho0)) + omega*(e + pres/rho_s)

        ! Raw c2 is returned so the init scan can catch non-positive values; c2_floor is the
        ! safety bound the wrappers apply (doubles as the stiffened-ambient cavitation cutoff).
        call s_jwl_floor(pres, sgm_eps)
        call s_jwl_floor(T, sgm_eps)
        c2_floor = min(air_gamma, merge(omega0, air_gamma, air_pi_inf > 0._wp))*(max(pres, sgm_eps) + pi_hat)/rho_s

    end subroutine s_jwl_state_er

    !> Analytic inverse (rho, p, Y) -> e. Coefficients depend on Y alone, so this is a single closed form.
    subroutine s_jwl_energy_pr(rho, pres, Y, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, e)

        $:GPU_ROUTINE(function_name='s_jwl_energy_pr',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, pres, Y, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a
        real(wp), intent(out) :: e
        real(wp)              :: rho_s, Y_s, An, Bn, omega, cv, pi_c, pi_hat, V, C1, C2

        rho_s = max(rho, sgm_eps)
        Y_s = min(max(Y, 0._wp), 1._wp)

        call s_jwl_mixture_coeffs(Y_s, A, B, omega0, air_gamma, air_pi_inf, cv_j, cv_a, An, Bn, omega, cv, pi_c, pi_hat)
        V = rho0/rho_s
        C1 = (1._wp - omega/(R1*V))*exp(-R1*V)
        C2 = (1._wp - omega/(R2*V))*exp(-R2*V)

        e = (pres + pi_c - An*C1 - Bn*C2)/max(omega*rho_s, sgm_eps)
        call s_jwl_floor(e, 0._wp)

    end subroutine s_jwl_energy_pr

    !> Fused (rho, p, Y) -> (e, c): one coefficient/exponential pass yields both the energy inverse and the forward sound speed, so
    !! the Riemann path gets the reconstructed energy and the wave-speed sound speed from a single call. Expressions mirror
    !! s_jwl_energy_pr and s_jwl_state_er exactly.
    subroutine s_jwl_sound_speed_pr(rho, pres, Y, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, e, c)

        $:GPU_ROUTINE(function_name='s_jwl_sound_speed_pr',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, pres, Y, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a
        real(wp), intent(out) :: e, c
        real(wp)              :: rho_s, Y_s, An, Bn, omega, cv, pi_c, pi_hat
        real(wp)              :: V, exp1, exp2, C1, C2, p_m, cs2, cs2_floor

        rho_s = max(rho, sgm_eps)
        Y_s = min(max(Y, 0._wp), 1._wp)

        call s_jwl_mixture_coeffs(Y_s, A, B, omega0, air_gamma, air_pi_inf, cv_j, cv_a, An, Bn, omega, cv, pi_c, pi_hat)
        V = rho0/rho_s
        exp1 = exp(-R1*V)
        exp2 = exp(-R2*V)
        C1 = (1._wp - omega/(R1*V))*exp1
        C2 = (1._wp - omega/(R2*V))*exp2

        e = (pres + pi_c - An*C1 - Bn*C2)/max(omega*rho_s, sgm_eps)
        call s_jwl_floor(e, 0._wp)

        p_m = An*C1 + Bn*C2 + omega*rho_s*e - pi_c
        cs2 = exp1*An*(R1*rho0/rho_s**2 - omega/rho_s - omega/(R1*rho0)) + exp2*Bn*(R2*rho0/rho_s**2 - omega/rho_s &
                       & - omega/(R2*rho0)) + omega*(e + p_m/rho_s)
        call s_jwl_floor(p_m, sgm_eps)
        cs2_floor = min(air_gamma, merge(omega0, air_gamma, air_pi_inf > 0._wp))*(max(p_m, sgm_eps) + pi_hat)/rho_s
        c = sqrt(max(cs2, cs2_floor))

    end subroutine s_jwl_sound_speed_pr

    ! Public entry points: look up fluid jidx's parameters, then evaluate the closure.

    !> Full state from energy for fluid jidx: (rho, e, Y) -> (p, T, [c]). c is optional so pressure-only callers skip the sqrt.
    subroutine s_jwl_mix_state_er(rho, e, Y, jidx, pres, T, c)

        $:GPU_ROUTINE(function_name='s_jwl_mix_state_er',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)            :: rho, e, Y
        integer, intent(in)             :: jidx
        real(wp), intent(out)           :: pres, T
        real(wp), intent(out), optional :: c
        real(wp)                        :: c2, c2_floor

        call s_jwl_state_er(rho, e, Y, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), jwl_omegas(jidx), &
                            & jwl_rho0s(jidx), jwl_air_gammas(jidx), jwl_air_pi_infs(jidx), jwl_cv_prod, jwl_cv_air, pres, T, c2, &
                            & c2_floor)
        ! The safety floor sits below any physical mixture c2, so it only engages on
        ! unphysical (e.g. cavitated) states without overriding legitimate sound speeds.
        if (present(c)) c = sqrt(max(c2, c2_floor))

    end subroutine s_jwl_mix_state_er

    !> Energy from pressure for fluid jidx: (rho, p, Y) -> e.
    subroutine s_jwl_mix_energy_pr(rho, pres, Y, jidx, e)

        $:GPU_ROUTINE(function_name='s_jwl_mix_energy_pr',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, pres, Y
        integer, intent(in)   :: jidx
        real(wp), intent(out) :: e

        call s_jwl_energy_pr(max(rho, sgm_eps), pres, min(max(Y, 0._wp), 1._wp), jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), &
                             & jwl_R2s(jidx), jwl_omegas(jidx), jwl_rho0s(jidx), jwl_air_gammas(jidx), jwl_air_pi_infs(jidx), &
                             & jwl_cv_prod, jwl_cv_air, e)

    end subroutine s_jwl_mix_energy_pr

    !> Sound speed for fluid jidx: invert to energy, then evaluate c.
    subroutine s_jwl_mix_sound_speed(rho, pres, Y, jidx, c)

        $:GPU_ROUTINE(function_name='s_jwl_mix_sound_speed',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, pres, Y
        integer, intent(in)   :: jidx
        real(wp), intent(out) :: c
        real(wp)              :: e_unused

        call s_jwl_sound_speed_pr(rho, pres, Y, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), jwl_omegas(jidx), &
                                  & jwl_rho0s(jidx), jwl_air_gammas(jidx), jwl_air_pi_infs(jidx), jwl_cv_prod, jwl_cv_air, &
                                  & e_unused, c)

    end subroutine s_jwl_mix_sound_speed

    !> Fused energy and sound speed for fluid jidx: (rho, p, Y) -> (e, c). The Riemann faces need both, so a single call shares the
    !! one coefficient/exponential pass instead of inverting for e and again for c.
    subroutine s_jwl_mix_energy_sound_speed_pr(rho, pres, Y, jidx, e, c)

        $:GPU_ROUTINE(function_name='s_jwl_mix_energy_sound_speed_pr',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, pres, Y
        integer, intent(in)   :: jidx
        real(wp), intent(out) :: e, c

        call s_jwl_sound_speed_pr(rho, pres, Y, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), jwl_omegas(jidx), &
                                  & jwl_rho0s(jidx), jwl_air_gammas(jidx), jwl_air_pi_infs(jidx), jwl_cv_prod, jwl_cv_air, e, c)

    end subroutine s_jwl_mix_energy_sound_speed_pr

    !> Initialize JWL parameter tables.
    impure subroutine s_initialize_jwl_module

        use m_mpi_common, only: s_mpi_abort
        use m_helper_basic, only: f_approx_equal, f_is_default

        integer  :: i, air_idx, gamma_src
        real(wp) :: air_e0_from_p0

        @:ALLOCATE(jwl_As(1:num_fluids), jwl_Bs(1:num_fluids), jwl_R1s(1:num_fluids), jwl_R2s(1:num_fluids), &
                   & jwl_omegas(1:num_fluids), jwl_rho0s(1:num_fluids), jwl_E0s(1:num_fluids))
        @:ALLOCATE(jwl_air_e0s(1:num_fluids), jwl_air_rho0s(1:num_fluids), jwl_air_gammas(1:num_fluids), &
                   & jwl_air_pi_infs(1:num_fluids))

        jwl_idx = 0
        air_idx = 0
        do i = 1, num_fluids
            if (fluid_pp(i)%eos == eos_jwl .and. .not. f_is_default(fluid_pp(i)%jwl_rho0)) then
                if (f_is_default(fluid_pp(i)%jwl_E0) .and. .not. f_is_default(fluid_pp(i)%jwl_Q)) then
                    fluid_pp(i)%jwl_E0 = fluid_pp(i)%jwl_rho0*fluid_pp(i)%jwl_Q
                else if (.not. f_is_default(fluid_pp(i)%jwl_E0) .and. f_is_default(fluid_pp(i)%jwl_Q)) then
                    fluid_pp(i)%jwl_Q = fluid_pp(i)%jwl_E0/fluid_pp(i)%jwl_rho0
                end if
            end if
            jwl_As(i) = fluid_pp(i)%jwl_A
            jwl_Bs(i) = fluid_pp(i)%jwl_B
            jwl_R1s(i) = fluid_pp(i)%jwl_R1
            jwl_R2s(i) = fluid_pp(i)%jwl_R2
            jwl_omegas(i) = fluid_pp(i)%jwl_omega
            jwl_rho0s(i) = fluid_pp(i)%jwl_rho0
            jwl_E0s(i) = fluid_pp(i)%jwl_E0
            jwl_air_e0s(i) = fluid_pp(i)%jwl_air_e0
            jwl_air_rho0s(i) = fluid_pp(i)%jwl_air_rho0
            jwl_air_gammas(i) = 0._wp
            jwl_air_pi_infs(i) = 0._wp
            if (fluid_pp(i)%eos == eos_jwl) then
                jwl_idx = i
            else
                if (air_idx == 0) air_idx = i
            end if
        end do

        jwl_cv_prod = 0._wp
        if (jwl_idx > 0) jwl_cv_prod = fluid_pp(jwl_idx)%cv
        if (air_idx > 0) then
            jwl_cv_air = fluid_pp(air_idx)%cv
        else
            jwl_cv_air = jwl_cv_prod  ! No ambient fluid: mass-weighted cv degenerates to cv_prod.
        end if

        if (jwl_idx > 0) then
            ! Ambient Grueneisen coefficient Gamma = 1/fluid_pp%gamma, from the non-JWL
            ! fluid (or the JWL fluid's own gamma when it is alone).
            if (air_idx > 0) then
                gamma_src = air_idx
            else
                gamma_src = jwl_idx
            end if
            if (f_is_default(fluid_pp(gamma_src)%gamma) .or. fluid_pp(gamma_src)%gamma <= 0._wp) then
                call s_mpi_abort('The JWL closure requires positive fluid_pp%gamma for the ambient-gas Grueneisen coefficient.')
            end if
            jwl_air_gammas(jwl_idx) = 1._wp/fluid_pp(gamma_src)%gamma

            ! True ambient stiffness: pi_inf = pi_inf_mfc/(fluid_pp%gamma + 1). Unset pi_inf means ideal gas.
            if (air_idx > 0) then
                if (.not. f_is_default(fluid_pp(air_idx)%pi_inf)) then
                    jwl_air_pi_infs(jwl_idx) = fluid_pp(air_idx)%pi_inf/(fluid_pp(air_idx)%gamma + 1._wp)
                end if
            end if

            ! Ambient energy from jwl_air_p0: e = (p*gamma_mfc + pi_inf_mfc)/rho.
            if (.not. f_is_default(fluid_pp(jwl_idx)%jwl_air_p0)) then
                air_e0_from_p0 = fluid_pp(jwl_idx)%jwl_air_p0*fluid_pp(gamma_src)%gamma/fluid_pp(jwl_idx)%jwl_air_rho0
                if (jwl_air_pi_infs(jwl_idx) > 0._wp) then
                    air_e0_from_p0 = (fluid_pp(jwl_idx)%jwl_air_p0*fluid_pp(gamma_src)%gamma + fluid_pp(air_idx)%pi_inf) &
                                      & /fluid_pp(jwl_idx)%jwl_air_rho0
                end if
                if (f_is_default(fluid_pp(jwl_idx)%jwl_air_e0)) then
                    jwl_air_e0s(jwl_idx) = air_e0_from_p0
                else if (.not. f_approx_equal(fluid_pp(jwl_idx)%jwl_air_e0, air_e0_from_p0, 1.e-8_wp)) then
                    call s_mpi_abort('fluid_pp%jwl_air_e0 must equal (jwl_air_p0*fluid_pp%gamma + fluid_pp%pi_inf)' &
                                     & // '/jwl_air_rho0 when both jwl_air_e0 and jwl_air_p0 are set.')
                end if
            end if

            if (jwl_rho0s(jwl_idx) <= jwl_air_rho0s(jwl_idx) .or. jwl_E0s(jwl_idx)/jwl_rho0s(jwl_idx) <= jwl_air_e0s(jwl_idx)) then
                call s_mpi_abort('The JWL closure requires increasing air-to-products reference density and energy.')
            end if

            ! Verify the assembled closure is positive-definite and invertible over the
            ! physical envelope before it is used by the solver.
            call s_jwl_verify_closure(jwl_idx)
        end if

        $:GPU_UPDATE(device='[jwl_As, jwl_Bs, jwl_R1s, jwl_R2s, jwl_omegas, jwl_rho0s, jwl_E0s]')
        $:GPU_UPDATE(device='[jwl_air_e0s, jwl_air_rho0s, jwl_air_gammas, jwl_air_pi_infs, jwl_idx]')
        $:GPU_UPDATE(device='[jwl_cv_prod, jwl_cv_air]')

    end subroutine s_initialize_jwl_module

    !> Init-time self-check: sweep the (rho, e, Y) envelope and abort if the assembled closure yields a non-positive sound speed or
    !! fails the (rho, p, Y) -> e -> p inverse.
    impure subroutine s_jwl_verify_closure(jidx)

        use m_mpi_common, only: s_mpi_abort

        integer, intent(in) :: jidx
        integer             :: ir, ie, iy
        real(wp)            :: rho_lo, rho_hi, e_lo, e_hi, ej, rho_s, e_s, pres, T, c2, c2_floor, e_inv
        character(len=128)  :: msg
        integer, parameter  :: n_scan = 100
        real(wp), parameter :: scan_rtol = 1.e-8_wp
        real(wp), parameter :: y_scan(8) = [0._wp, 0.25_wp, 0.5_wp, 0.75_wp, 0.9_wp, 0.97_wp, 0.999_wp, 1._wp]

        ej = jwl_E0s(jidx)/jwl_rho0s(jidx)
        ! Envelope covers reflected-shock states beyond the CJ point; the density cap is
        ! tighter for a stiffened ambient (JWL is only meaningful to ~2*rho0 there).
        rho_lo = 0.1_wp*jwl_air_rho0s(jidx)
        rho_hi = 4._wp*jwl_rho0s(jidx)
        if (jwl_air_pi_infs(jidx) > 0._wp) rho_hi = 2._wp*jwl_rho0s(jidx)
        e_lo = 0.5_wp*jwl_air_e0s(jidx)
        e_hi = 5._wp*ej

        do iy = 1, size(y_scan)
            do ir = 0, n_scan - 1
                rho_s = rho_lo*(rho_hi/rho_lo)**(real(ir, wp)/real(n_scan - 1, wp))
                do ie = 0, n_scan - 1
                    e_s = e_lo + (e_hi - e_lo)*real(ie, wp)/real(n_scan - 1, wp)
                    call s_jwl_state_er(rho_s, e_s, y_scan(iy), jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), &
                                        & jwl_omegas(jidx), jwl_rho0s(jidx), jwl_air_gammas(jidx), jwl_air_pi_infs(jidx), &
                                        & jwl_cv_prod, jwl_cv_air, pres, T, c2, c2_floor)
                    ! Floored c2 must always be positive; raw c2 must be positive wherever
                    ! the pressure was not floored (floored = cavitated, handled by cutoff).
                    if (max(c2, c2_floor) <= 0._wp .or. c2_floor /= c2_floor) then
                        write (msg, &
                               & '(A,ES11.4,A,ES11.4,A,F6.3)') &
                               & 'JWL closure self-check: non-positive floored sound speed at rho=', rho_s, ', e=', e_s, ', Y=', &
                               & y_scan(iy)
                        call s_mpi_abort(trim(msg) // '. Check JWL and ambient parameters.')
                    end if
                    if (pres > sgm_eps .and. (c2 <= 0._wp .or. c2 /= c2)) then
                        write (msg, '(A,ES11.4,A,ES11.4,A,F6.3)') 'JWL closure self-check: non-positive sound speed at rho=', &
                               & rho_s, ', e=', e_s, ', Y=', y_scan(iy)
                        call s_mpi_abort(trim(msg) // '. Check JWL and ambient parameters.')
                    end if
                    if (pres > sgm_eps) then
                        call s_jwl_energy_pr(rho_s, pres, y_scan(iy), jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), &
                                             & jwl_omegas(jidx), jwl_rho0s(jidx), jwl_air_gammas(jidx), jwl_air_pi_infs(jidx), &
                                             & jwl_cv_prod, jwl_cv_air, e_inv)
                        if (e_inv /= e_inv .or. abs(e_inv - e_s) > scan_rtol*max(abs(e_s), jwl_air_e0s(jidx))) then
                            write (msg, &
                                   & '(A,ES11.4,A,ES11.4,A,F6.3)') 'JWL closure self-check: energy inversion mismatch at rho=', &
                                   & rho_s, ', e=', e_s, ', Y=', y_scan(iy)
                            call s_mpi_abort(trim(msg) // '. Check JWL and ambient-gas parameters.')
                        end if
                    end if
                end do
            end do
        end do

    end subroutine s_jwl_verify_closure

    !> Deallocate the per-fluid JWL parameter tables.
    impure subroutine s_finalize_jwl_module

        @:DEALLOCATE(jwl_As, jwl_Bs, jwl_R1s, jwl_R2s, jwl_omegas, jwl_rho0s, jwl_E0s)
        @:DEALLOCATE(jwl_air_e0s, jwl_air_rho0s, jwl_air_gammas, jwl_air_pi_infs)

    end subroutine s_finalize_jwl_module

end module m_jwl
