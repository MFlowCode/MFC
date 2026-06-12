!>
!! @file
!! @brief Contains module m_jwl

#:include 'macros.fpp'
#:include 'case.fpp'

!> @brief Jones-Wilkins-Lee (JWL) equation of state and JWL/ideal-gas mixture closures. Holds the per-fluid JWL parameter tables and
!! the pressure<->energy / sound-speed routines for every supported mixture closure (jwl_mix_type). Split out of
!! m_variables_conversion so JWL physics lives in one place.
module m_jwl

    use m_global_parameters

    implicit none

    private
    public :: s_initialize_jwl_module, s_finalize_jwl_module, s_jwl_pcold, s_jwl_sound_speed_squared, &
        & s_jwl_mixture_sound_speed_squared, s_jwl_pressure_er, s_jwl_energy_pr, s_jwl_mix_pressure_er, s_jwl_mix_energy_pr, &
        & s_jwl_kuhl_pressure_er, s_jwl_kuhl_energy_pr, s_jwl_kuhl_sound_speed_squared, jwl_idx

    ! Per-fluid JWL parameter tables. In simulation these live in m_global_parameters
    ! (shared with the solver hot path); for pre/post_process they are declared here.
#ifndef MFC_SIMULATION
    real(wp), allocatable, public, dimension(:) :: jwl_As, jwl_Bs, jwl_R1s, jwl_R2s, jwl_omegas, jwl_rho0s, jwl_E0s
    real(wp), allocatable, public, dimension(:) :: jwl_air_e0s, jwl_air_rho0s, jwl_air_gammas
    $:GPU_DECLARE(create='[jwl_As, jwl_Bs, jwl_R1s, jwl_R2s, jwl_omegas, jwl_rho0s, jwl_E0s]')
    $:GPU_DECLARE(create='[jwl_air_e0s, jwl_air_rho0s, jwl_air_gammas]')
#endif

    integer  :: jwl_idx                  !< Index of the JWL fluid (0 if none)
    real(wp) :: jwl_cv_prod, jwl_cv_air  !< Products/air specific heats for the p-T-equilibrium closure (jwl_mix_type = 2)
    $:GPU_DECLARE(create='[jwl_idx]')
    $:GPU_DECLARE(create='[jwl_cv_prod, jwl_cv_air]')

contains

    !> Compute the JWL cold pressure term with relative volume V = rho0/rho.
    subroutine s_jwl_pcold(rho, A, B, R1, R2, omega, rho0, pcold)

        $:GPU_ROUTINE(function_name='s_jwl_pcold',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, A, B, R1, R2, omega, rho0
        real(wp), intent(out) :: pcold
        real(wp)              :: V

        V = rho0/max(rho, sgm_eps)
        pcold = A*(1._wp - omega/(R1*V))*exp(-R1*V) + B*(1._wp - omega/(R2*V))*exp(-R2*V)

    end subroutine s_jwl_pcold

    !> Compute d(pcold)/d(rho) for the JWL EOS.
    subroutine s_jwl_dpcold_drho(rho, A, B, R1, R2, omega, rho0, dpcold_drho)

        $:GPU_ROUTINE(function_name='s_jwl_dpcold_drho',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, A, B, R1, R2, omega, rho0
        real(wp), intent(out) :: dpcold_drho
        real(wp)              :: V, rho_safe

        rho_safe = max(rho, sgm_eps)
        V = rho0/rho_safe
        dpcold_drho = A*exp(-R1*V)*(V/rho_safe)*(R1 - omega/V - omega/(R1*V**2)) + B*exp(-R2*V)*(V/rho_safe)*(R2 - omega/V &
                            & - omega/(R2*V**2))

    end subroutine s_jwl_dpcold_drho

    !> Compute the JWL isentropic sound-speed squared using the Rocflu/Stanley form.
    subroutine s_jwl_sound_speed_squared(rho, pres, A, B, R1, R2, omega, rho0, c2)

        $:GPU_ROUTINE(function_name='s_jwl_sound_speed_squared',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, pres, A, B, R1, R2, omega, rho0
        real(wp), intent(out) :: c2
        real(wp)              :: rho_safe, exp1, exp2, v, e

        rho_safe = max(rho, sgm_eps)
        v = rho0/rho_safe
        exp1 = exp(-R1*v)
        exp2 = exp(-R2*v)

        e = (pres - A*(1._wp - omega/(R1*v))*exp1 - B*(1._wp - omega/(R2*v))*exp2)/(omega*rho_safe)
        c2 = A*exp1*(omega/(R1*v**2._wp) - R1*(1._wp - omega/(R1*v))) + B*exp2*(omega/(R2*v**2._wp) - R2*(1._wp - omega/(R2*v)))
        c2 = c2*(-rho0/(rho_safe*rho_safe)) + omega*e + omega*pres/rho_safe
        c2 = max(c2, sgm_eps)

    end subroutine s_jwl_sound_speed_squared

    !> Sound-speed squared for mixed JWL/ideal-gas states using the frozen (mass-weighted) mixture rule, c^2 = sum_k Y_k*c_k^2, with
    !! each phase sound speed evaluated at its own density rho_k = (Y_k*rho)/alpha_k. The frozen estimate is smooth and monotone
    !! between the phase values (unlike Wood's equilibrium rule, whose sharp dip at intermediate alpha collapses the HLLC wave-speed
    !! estimates and drives spurious velocity oscillations at the contact). It recovers each single-material limit and only feeds
    !! the Riemann wave-speed/dissipation estimate, not the isobaric pressure closure.
    subroutine s_jwl_mixture_sound_speed_squared(rho, pres, Y, alpha_j, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, &
        & air_gamma, c2)

        $:GPU_ROUTINE(function_name='s_jwl_mixture_sound_speed_squared',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, pres, Y, alpha_j, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma
        real(wp), intent(out) :: c2
        real(wp)              :: rho_safe, Y_safe, a_j, a_a, rho1, rho2, c2_air, c2_jwl

        rho_safe = max(rho, sgm_eps)
        Y_safe = min(max(Y, 0._wp), 1._wp)
        a_j = min(max(alpha_j, 0._wp), 1._wp)
        a_a = 1._wp - a_j

        if (a_j <= sgm_eps) then
            c2 = max((air_gamma + 1._wp)*pres/rho_safe, sgm_eps)
            return
        end if

        if (a_a <= sgm_eps) then
            call s_jwl_sound_speed_squared(rho_safe, pres, A, B, R1, R2, omega0, rho0, c2)
            c2 = max(c2, sgm_eps)
            return
        end if

        rho1 = max(Y_safe*rho_safe/a_j, sgm_eps)
        rho2 = max((1._wp - Y_safe)*rho_safe/a_a, sgm_eps)

        call s_jwl_sound_speed_squared(rho1, pres, A, B, R1, R2, omega0, rho0, c2_jwl)
        c2_air = max((air_gamma + 1._wp)*pres/rho2, sgm_eps)

        ! Frozen (mass-weighted) mixture sound speed: smooth and monotone between the phase values, avoiding Wood's interface dip
        c2 = Y_safe*max(c2_jwl, sgm_eps) + (1._wp - Y_safe)*c2_air
        c2 = max(c2, sgm_eps)

    end subroutine s_jwl_mixture_sound_speed_squared

    !> JWL/ideal-gas mixture pressure from specific internal energy and density via the closed-form isobaric
    !! (mechanical-equilibrium) closure. Each phase is evaluated at its own density rho_k = (Y_k*rho)/alpha_k, so the JWL cold curve
    !! is taken at the dense-products density and air at its own density. Both EOS are linear in pressure, so the common pressure
    !! that partitions the cell internal energy is exact (no iteration) and recovers pure air (alpha_j->0) and pure JWL (alpha_j->1)
    !! exactly.
    subroutine s_jwl_pressure_er(rho, e, Y, alpha_j, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma, pres)

        $:GPU_ROUTINE(function_name='s_jwl_pressure_er',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, e, Y, alpha_j, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma
        real(wp), intent(out) :: pres
        real(wp)              :: rho_safe, Y_safe, a_j, a_a, rho1, rhoe, pref1, Kj

        rho_safe = max(rho, sgm_eps)
        Y_safe = min(max(Y, 0._wp), 1._wp)
        a_j = min(max(alpha_j, 0._wp), 1._wp)
        a_a = 1._wp - a_j
        rhoe = rho_safe*e

        if (a_j <= sgm_eps) then
            pres = air_gamma*rhoe
            return
        end if

        if (a_a <= sgm_eps) then
            rho1 = rho_safe
        else
            rho1 = max(Y_safe*rho_safe/a_j, sgm_eps)
        end if

        pref1 = A*(1._wp - omega0*rho1/(R1*rho0))*exp(-R1*rho0/rho1) + B*(1._wp - omega0*rho1/(R2*rho0))*exp(-R2*rho0/rho1)
        Kj = 1._wp/max(omega0, sgm_eps)

        pres = (rhoe + a_j*Kj*pref1)/max(a_j*Kj + a_a/max(air_gamma, sgm_eps), sgm_eps)

    end subroutine s_jwl_pressure_er

    !> Specific internal energy from pressure and density: exact inverse of the isobaric mixture pressure in s_jwl_pressure_er,
    !! keeping primitive<->conservative consistent.
    subroutine s_jwl_energy_pr(rho, pres, Y, alpha_j, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma, e)

        $:GPU_ROUTINE(function_name='s_jwl_energy_pr',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, pres, Y, alpha_j, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma
        real(wp), intent(out) :: e
        real(wp)              :: rho_safe, Y_work, a_j, a_a, rho1, pref1, Kj

        rho_safe = max(rho, sgm_eps)
        Y_work = min(max(Y, 0._wp), 1._wp)
        a_j = min(max(alpha_j, 0._wp), 1._wp)
        a_a = 1._wp - a_j

        if (a_j <= sgm_eps) then
            e = max(pres/max(air_gamma*rho_safe, sgm_eps), 0._wp)
            return
        end if

        if (a_a <= sgm_eps) then
            rho1 = rho_safe
        else
            rho1 = max(Y_work*rho_safe/a_j, sgm_eps)
        end if

        pref1 = A*(1._wp - omega0*rho1/(R1*rho0))*exp(-R1*rho0/rho1) + B*(1._wp - omega0*rho1/(R2*rho0))*exp(-R2*rho0/rho1)
        Kj = 1._wp/max(omega0, sgm_eps)

        e = (pres*(a_j*Kj + a_a/max(air_gamma, sgm_eps)) - a_j*Kj*pref1)/rho_safe
        e = max(e, 0._wp)

    end subroutine s_jwl_energy_pr

    !> Kuhl/Khasainov temperature-form JWL/ideal-gas mixture pressure. This is the PDF-style rho*R*T closure: p = Y*(A*exp(-R1*V) +
    !! B*exp(-R2*V)) + rho*T*(Y*omega*Cv_j + (1-Y)*R_air), with e = Y*e_cold(rho) + (Y*Cv_j + (1-Y)*Cv_air)*T. Both phases share the
    !! cell density rho and alpha_j is unused.
    subroutine s_jwl_kuhl_pressure_er(rho, e, Y, alpha_j, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma, pres)

        $:GPU_ROUTINE(function_name='s_jwl_kuhl_pressure_er',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, e, Y, alpha_j, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma
        real(wp), intent(out) :: pres
        real(wp)              :: rho_safe, Y_safe, V, exp1, exp2, pbase, ecold, cv_mix, R_mix, T

        rho_safe = max(rho, sgm_eps)
        Y_safe = min(max(Y, 0._wp), 1._wp)
        V = rho0/rho_safe
        exp1 = exp(-R1*V)
        exp2 = exp(-R2*V)
        pbase = A*exp1 + B*exp2
        ecold = A/(R1*rho0)*exp1 + B/(R2*rho0)*exp2
        cv_mix = max(Y_safe*jwl_cv_prod + (1._wp - Y_safe)*jwl_cv_air, sgm_eps)
        R_mix = Y_safe*omega0*jwl_cv_prod + (1._wp - Y_safe)*air_gamma*jwl_cv_air
        T = max((e - Y_safe*ecold)/cv_mix, sgm_eps)
        pres = max(Y_safe*pbase + rho_safe*R_mix*T, sgm_eps)

    end subroutine s_jwl_kuhl_pressure_er

    !> Specific internal energy from pressure and density: exact inverse of the Kuhl/Khasainov temperature-form closure. alpha_j is
    !! unused.
    subroutine s_jwl_kuhl_energy_pr(rho, pres, Y, alpha_j, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma, e)

        $:GPU_ROUTINE(function_name='s_jwl_kuhl_energy_pr',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, pres, Y, alpha_j, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma
        real(wp), intent(out) :: e
        real(wp)              :: rho_safe, Y_safe, p_safe, V, exp1, exp2, pbase, ecold, cv_mix, R_mix, T

        rho_safe = max(rho, sgm_eps)
        Y_safe = min(max(Y, 0._wp), 1._wp)
        p_safe = max(pres, sgm_eps)
        V = rho0/rho_safe
        exp1 = exp(-R1*V)
        exp2 = exp(-R2*V)
        pbase = A*exp1 + B*exp2
        ecold = A/(R1*rho0)*exp1 + B/(R2*rho0)*exp2
        cv_mix = max(Y_safe*jwl_cv_prod + (1._wp - Y_safe)*jwl_cv_air, sgm_eps)
        R_mix = Y_safe*omega0*jwl_cv_prod + (1._wp - Y_safe)*air_gamma*jwl_cv_air
        T = max((p_safe - Y_safe*pbase)/max(rho_safe*R_mix, sgm_eps), sgm_eps)
        e = max(Y_safe*ecold + cv_mix*T, 0._wp)

    end subroutine s_jwl_kuhl_energy_pr

    !> Sound speed for the Kuhl/Khasainov temperature-form closure, treated as a Mie-Gruneisen EOS p = P_eff(rho,Y) +
    !! Gamma_eff*rho*e with constant-Y Gamma_eff = R_mix/Cv_mix.
    subroutine s_jwl_kuhl_sound_speed_squared(rho, pres, Y, alpha_j, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma, &
        & c2)

        $:GPU_ROUTINE(function_name='s_jwl_kuhl_sound_speed_squared',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, pres, Y, alpha_j, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma
        real(wp), intent(out) :: c2
        real(wp)              :: rho_safe, Y_safe, V, exp1, exp2, pbase, ecold, dpbase, decold, cv_mix, R_mix, Gamma_eff, e, dPeff

        rho_safe = max(rho, sgm_eps)
        Y_safe = min(max(Y, 0._wp), 1._wp)
        V = rho0/rho_safe
        exp1 = exp(-R1*V)
        exp2 = exp(-R2*V)
        pbase = A*exp1 + B*exp2
        ecold = A/(R1*rho0)*exp1 + B/(R2*rho0)*exp2
        dpbase = (A*R1*rho0*exp1 + B*R2*rho0*exp2)/(rho_safe*rho_safe)
        decold = pbase/(rho_safe*rho_safe)
        cv_mix = max(Y_safe*jwl_cv_prod + (1._wp - Y_safe)*jwl_cv_air, sgm_eps)
        R_mix = Y_safe*omega0*jwl_cv_prod + (1._wp - Y_safe)*air_gamma*jwl_cv_air
        Gamma_eff = R_mix/cv_mix
        e = max((pres - Y_safe*(pbase - Gamma_eff*rho_safe*ecold))/max(Gamma_eff*rho_safe, sgm_eps), 0._wp)
        dPeff = Y_safe*dpbase - Gamma_eff*Y_safe*(ecold + rho_safe*decold)
        c2 = max(dPeff + Gamma_eff*(e + pres/rho_safe), sgm_eps)

    end subroutine s_jwl_kuhl_sound_speed_squared

    !> JWL/ideal-gas mixture pressure under the *standard* thermal-and-pressure-equilibrium closure (both phases share one T and one
    !! p; Garno 2020, ref [22]). Given the cell density, energy, and products mass fraction Y, we find the products volume fraction
    !! alpha_j that balances the two phase pressures by bisection. The shared temperature is closed-form from the energy split (both
    !! caloric EOS are linear in T), so each step is cheap. We use the temperature-explicit JWL form p_j = A e^{-R1 V} + B e^{-R2 V}
    !! + omega*rho_j*Cv_j*T (algebraically identical to Garno 7a) and floor T positive, which keeps the bracket valid for any input.
    subroutine s_jwl_ptequil_pressure_er(rho, e, Y, A, B, R1, R2, omega0, rho0, air_gamma, cv_j, cv_a, pres)

        $:GPU_ROUTINE(function_name='s_jwl_ptequil_pressure_er',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, e, Y, A, B, R1, R2, omega0, rho0, air_gamma, cv_j, cv_a
        real(wp), intent(out) :: pres
        real(wp)              :: rho_s, Y_s, cv_mix, a_lo, a_hi, a_m, rj, ra, V, ecold, T, pj, pa, f_lo, f_m, pcg
        integer               :: it

        rho_s = max(rho, sgm_eps)
        Y_s = min(max(Y, 0._wp), 1._wp)

        ! Pure cells: no equilibrium to solve, evaluate the single-phase EOS directly.
        if (Y_s <= sgm_eps) then
            pres = max(air_gamma*rho_s*e, sgm_eps); return
        else if (1._wp - Y_s <= sgm_eps) then
            call s_jwl_pcold(rho_s, A, B, R1, R2, omega0, rho0, pcg)
            pres = max(pcg + omega0*rho_s*e, sgm_eps); return
        end if

        cv_mix = max(Y_s*cv_j + (1._wp - Y_s)*cv_a, sgm_eps)

        ! Bracket alpha_j in (0,1); residual p_j - p_a is monotone-decreasing, so bisect it to zero.
        a_lo = sgm_eps; a_hi = 1._wp - sgm_eps
        rj = max(Y_s*rho_s/a_lo, sgm_eps); ra = max((1._wp - Y_s)*rho_s/(1._wp - a_lo), sgm_eps); V = rho0/rj
        ecold = A/(R1*rho0)*exp(-R1*V) + B/(R2*rho0)*exp(-R2*V)
        T = max((e - Y_s*ecold)/cv_mix, sgm_eps)
        pj = A*exp(-R1*V) + B*exp(-R2*V) + omega0*rj*cv_j*T; pa = air_gamma*ra*cv_a*T
        f_lo = pj - pa; pres = 0.5_wp*(pj + pa)

        do it = 1, 60
            a_m = 0.5_wp*(a_lo + a_hi)
            rj = max(Y_s*rho_s/a_m, sgm_eps); ra = max((1._wp - Y_s)*rho_s/(1._wp - a_m), sgm_eps); V = rho0/rj
            ecold = A/(R1*rho0)*exp(-R1*V) + B/(R2*rho0)*exp(-R2*V)
            T = max((e - Y_s*ecold)/cv_mix, sgm_eps)
            pj = A*exp(-R1*V) + B*exp(-R2*V) + omega0*rj*cv_j*T; pa = air_gamma*ra*cv_a*T
            f_m = pj - pa; pres = 0.5_wp*(pj + pa)
            if (f_lo*f_m > 0._wp) then
                a_lo = a_m; f_lo = f_m
            else
                a_hi = a_m
            end if
        end do
        pres = max(pres, sgm_eps)

    end subroutine s_jwl_ptequil_pressure_er

    !> Specific internal energy for the p-T-equilibrium mixture: inverse of s_jwl_ptequil_pressure_er. With p known, the air phase
    !! fixes T = p/(air_gamma*rho_a*Cv_a) in closed form, so we bisect alpha_j on the products-pressure residual p_j(alpha_j) - p,
    !! then assemble e = (Y*Cv_j + (1-Y)*Cv_a)*T + Y*e_cold(rho_j).
    subroutine s_jwl_ptequil_energy_pr(rho, pres, Y, A, B, R1, R2, omega0, rho0, air_gamma, cv_j, cv_a, e)

        $:GPU_ROUTINE(function_name='s_jwl_ptequil_energy_pr',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, pres, Y, A, B, R1, R2, omega0, rho0, air_gamma, cv_j, cv_a
        real(wp), intent(out) :: e
        real(wp)              :: rho_s, Y_s, p_s, a_lo, a_hi, a_m, rj, ra, V, ecold, T, pj, pcg, g_lo, g_m
        integer               :: it

        rho_s = max(rho, sgm_eps)
        Y_s = min(max(Y, 0._wp), 1._wp)
        p_s = max(pres, sgm_eps)

        if (Y_s <= sgm_eps) then
            e = max(p_s/max(air_gamma*rho_s, sgm_eps), 0._wp); return
        else if (1._wp - Y_s <= sgm_eps) then
            call s_jwl_pcold(rho_s, A, B, R1, R2, omega0, rho0, pcg)
            e = max((p_s - pcg)/max(omega0*rho_s, sgm_eps), 0._wp); return
        end if

        a_lo = sgm_eps; a_hi = 1._wp - sgm_eps
        rj = max(Y_s*rho_s/a_lo, sgm_eps); ra = max((1._wp - Y_s)*rho_s/(1._wp - a_lo), sgm_eps); V = rho0/rj
        T = p_s/max(air_gamma*ra*cv_a, sgm_eps)
        pj = A*exp(-R1*V) + B*exp(-R2*V) + omega0*rj*cv_j*T
        g_lo = pj - p_s

        do it = 1, 60
            a_m = 0.5_wp*(a_lo + a_hi)
            rj = max(Y_s*rho_s/a_m, sgm_eps); ra = max((1._wp - Y_s)*rho_s/(1._wp - a_m), sgm_eps); V = rho0/rj
            T = p_s/max(air_gamma*ra*cv_a, sgm_eps)
            pj = A*exp(-R1*V) + B*exp(-R2*V) + omega0*rj*cv_j*T
            g_m = pj - p_s
            if (g_lo*g_m > 0._wp) then
                a_lo = a_m; g_lo = g_m
            else
                a_hi = a_m
            end if
        end do

        ecold = A/(R1*rho0)*exp(-R1*V) + B/(R2*rho0)*exp(-R2*V)
        e = max((Y_s*cv_j + (1._wp - Y_s)*cv_a)*T + Y_s*ecold, 0._wp)

    end subroutine s_jwl_ptequil_energy_pr

    !> Garno (2020) "Rocflu" single-fluid blended JWL/ideal-gas EOS. In mixture cells, A and B ramp linearly with mixture specific
    !! internal energy (g_e: 0 at ambient air e_a -> 1 at explosive e0 = E0/rho0), while omega ramps linearly with mixture density
    !! (g_rho: air -> products). Garno applies the adapted coefficients only below 99% products by mass, so near-pure products use
    !! the unmodified JWL EOS; near-pure air uses the ideal-gas EOS to avoid numerical trace contamination.
    subroutine s_jwl_rocflu_pressure_er(rho, e, Y, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma, pres)

        $:GPU_ROUTINE(function_name='s_jwl_rocflu_pressure_er',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, e, Y, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma
        real(wp), intent(out) :: pres
        real(wp)              :: rho_s, Y_s, V, g_rho, g_e, om, e0s, cab

        rho_s = max(rho, sgm_eps)
        Y_s = min(max(Y, 0._wp), 1._wp)

        if (Y_s <= 1.e-4_wp) then
            pres = max(air_gamma*rho_s*e, sgm_eps)
            return
        else if (Y_s >= 0.99_wp) then
            call s_jwl_pcold(rho_s, A, B, R1, R2, omega0, rho0, cab)
            pres = max(cab + omega0*rho_s*e, sgm_eps)
            return
        end if

        V = rho0/rho_s
        g_rho = min(max((rho_s - air_rho0)/max(rho0 - air_rho0, sgm_eps), 0._wp), 1._wp)
        om = air_gamma + (omega0 - air_gamma)*g_rho
        e0s = E0/max(rho0, sgm_eps)
        g_e = min(max((e - air_e0)/max(e0s - air_e0, sgm_eps), 0._wp), 1._wp)
        cab = A*(1._wp - om/(R1*V))*exp(-R1*V) + B*(1._wp - om/(R2*V))*exp(-R2*V)
        pres = max(g_e*cab + om*rho_s*e, sgm_eps)

    end subroutine s_jwl_rocflu_pressure_er

    !> Specific internal energy for the Rocflu blend: inverse of s_jwl_rocflu_pressure_er. The pure-material guards mirror pressure
    !! recovery; in the mixed branch p is linear in e while the energy blend is unsaturated, so we invert linearly then correct if
    !! g_e clamps.
    subroutine s_jwl_rocflu_energy_pr(rho, pres, Y, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma, e)

        $:GPU_ROUTINE(function_name='s_jwl_rocflu_energy_pr',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, pres, Y, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma
        real(wp), intent(out) :: e
        real(wp)              :: rho_s, Y_s, V, g_rho, om, e0s, cab, kab, g_e

        rho_s = max(rho, sgm_eps)
        Y_s = min(max(Y, 0._wp), 1._wp)

        if (Y_s <= 1.e-4_wp) then
            e = max(pres/max(air_gamma*rho_s, sgm_eps), 0._wp)
            return
        else if (Y_s >= 0.99_wp) then
            call s_jwl_pcold(rho_s, A, B, R1, R2, omega0, rho0, cab)
            e = max((pres - cab)/max(omega0*rho_s, sgm_eps), 0._wp)
            return
        end if

        V = rho0/rho_s
        g_rho = min(max((rho_s - air_rho0)/max(rho0 - air_rho0, sgm_eps), 0._wp), 1._wp)
        om = air_gamma + (omega0 - air_gamma)*g_rho
        e0s = E0/max(rho0, sgm_eps)
        cab = A*(1._wp - om/(R1*V))*exp(-R1*V) + B*(1._wp - om/(R2*V))*exp(-R2*V)
        kab = cab/max(e0s - air_e0, sgm_eps)
        e = (pres + air_e0*kab)/max(kab + om*rho_s, sgm_eps)
        g_e = (e - air_e0)/max(e0s - air_e0, sgm_eps)
        if (g_e < 0._wp) then
            e = pres/max(om*rho_s, sgm_eps)  ! blend saturated to pure air
        else if (g_e > 1._wp) then
            e = (pres - cab)/max(om*rho_s, sgm_eps)  ! blend saturated to full products
        end if
        e = max(e, 0._wp)

    end subroutine s_jwl_rocflu_energy_pr

    !> Dispatch the JWL/ideal-gas mixture pressure-from-energy to the active closure jwl_mix_type: 0 isobaric (default), 1 Kuhl
    !! additive, 2 thermal+pressure equilibrium, 3 Rocflu blend. Keeps the per-mode argument lists in one place so call sites stay
    !! one line. jidx is the JWL fluid index into the module jwl_*s arrays.
    subroutine s_jwl_mix_pressure_er(rho, e, Y, alpha_j, jidx, pres)

        $:GPU_ROUTINE(function_name='s_jwl_mix_pressure_er',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, e, Y, alpha_j
        integer, intent(in)   :: jidx
        real(wp), intent(out) :: pres

        select case (jwl_mix_type)
        case (1)
            call s_jwl_kuhl_pressure_er(rho, e, Y, alpha_j, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), &
                                        & jwl_omegas(jidx), jwl_rho0s(jidx), jwl_E0s(jidx), jwl_air_e0s(jidx), &
                                        & jwl_air_rho0s(jidx), jwl_air_gammas(jidx), pres)
        case (2)
            call s_jwl_ptequil_pressure_er(rho, e, Y, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), jwl_omegas(jidx), &
                                           & jwl_rho0s(jidx), jwl_air_gammas(jidx), jwl_cv_prod, jwl_cv_air, pres)
        case (3)
            call s_jwl_rocflu_pressure_er(rho, e, Y, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), jwl_omegas(jidx), &
                                          & jwl_rho0s(jidx), jwl_E0s(jidx), jwl_air_e0s(jidx), jwl_air_rho0s(jidx), &
                                          & jwl_air_gammas(jidx), pres)
        case default
            call s_jwl_pressure_er(rho, e, Y, alpha_j, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), &
                                   & jwl_omegas(jidx), jwl_rho0s(jidx), jwl_E0s(jidx), jwl_air_e0s(jidx), jwl_air_rho0s(jidx), &
                                   & jwl_air_gammas(jidx), pres)
        end select

    end subroutine s_jwl_mix_pressure_er

    !> Dispatch the JWL/ideal-gas mixture energy-from-pressure to the active closure jwl_mix_type (inverse of
    !! s_jwl_mix_pressure_er).
    subroutine s_jwl_mix_energy_pr(rho, pres, Y, alpha_j, jidx, e)

        $:GPU_ROUTINE(function_name='s_jwl_mix_energy_pr',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, pres, Y, alpha_j
        integer, intent(in)   :: jidx
        real(wp), intent(out) :: e

        select case (jwl_mix_type)
        case (1)
            call s_jwl_kuhl_energy_pr(rho, pres, Y, alpha_j, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), &
                                      & jwl_omegas(jidx), jwl_rho0s(jidx), jwl_E0s(jidx), jwl_air_e0s(jidx), jwl_air_rho0s(jidx), &
                                      & jwl_air_gammas(jidx), e)
        case (2)
            call s_jwl_ptequil_energy_pr(rho, pres, Y, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), &
                                         & jwl_omegas(jidx), jwl_rho0s(jidx), jwl_air_gammas(jidx), jwl_cv_prod, jwl_cv_air, e)
        case (3)
            call s_jwl_rocflu_energy_pr(rho, pres, Y, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), jwl_omegas(jidx), &
                                        & jwl_rho0s(jidx), jwl_E0s(jidx), jwl_air_e0s(jidx), jwl_air_rho0s(jidx), &
                                        & jwl_air_gammas(jidx), e)
        case default
            call s_jwl_energy_pr(rho, pres, Y, alpha_j, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), &
                                 & jwl_omegas(jidx), jwl_rho0s(jidx), jwl_E0s(jidx), jwl_air_e0s(jidx), jwl_air_rho0s(jidx), &
                                 & jwl_air_gammas(jidx), e)
        end select

    end subroutine s_jwl_mix_energy_pr

    !> Allocate and populate the per-fluid JWL parameter tables from fluid_pp, identify the JWL fluid (jwl_idx), and cache the
    !! products/air specific heats used by the p-T-equilibrium closure. Called from the variables-conversion module initializer.
    impure subroutine s_initialize_jwl_module

        integer :: i

        @:ALLOCATE(jwl_As    (1:num_fluids))
        @:ALLOCATE(jwl_Bs    (1:num_fluids))
        @:ALLOCATE(jwl_R1s   (1:num_fluids))
        @:ALLOCATE(jwl_R2s   (1:num_fluids))
        @:ALLOCATE(jwl_omegas(1:num_fluids))
        @:ALLOCATE(jwl_rho0s (1:num_fluids))
        @:ALLOCATE(jwl_E0s   (1:num_fluids))
        @:ALLOCATE(jwl_air_e0s    (1:num_fluids))
        @:ALLOCATE(jwl_air_rho0s  (1:num_fluids))
        @:ALLOCATE(jwl_air_gammas (1:num_fluids))

        jwl_idx = 0
        do i = 1, num_fluids
            jwl_As(i) = fluid_pp(i)%jwl_A
            jwl_Bs(i) = fluid_pp(i)%jwl_B
            jwl_R1s(i) = fluid_pp(i)%jwl_R1
            jwl_R2s(i) = fluid_pp(i)%jwl_R2
            jwl_omegas(i) = fluid_pp(i)%jwl_omega
            jwl_rho0s(i) = fluid_pp(i)%jwl_rho0
            jwl_E0s(i) = fluid_pp(i)%jwl_E0
            jwl_air_e0s(i) = fluid_pp(i)%jwl_air_e0
            jwl_air_rho0s(i) = fluid_pp(i)%jwl_air_rho0
            jwl_air_gammas(i) = fluid_pp(i)%jwl_air_gamma
            if (fluid_pp(i)%eos == 2) jwl_idx = i
        end do

        ! Specific heats for the p-T-equilibrium closure: products from the JWL fluid, air from the first ideal-gas fluid.
        jwl_cv_prod = 0._wp; jwl_cv_air = 0._wp
        if (jwl_idx > 0) jwl_cv_prod = fluid_pp(jwl_idx)%cv
        do i = 1, num_fluids
            if (fluid_pp(i)%eos == 1) then
                jwl_cv_air = fluid_pp(i)%cv; exit
            end if
        end do

        $:GPU_UPDATE(device='[jwl_As, jwl_Bs, jwl_R1s, jwl_R2s, jwl_omegas, jwl_rho0s, jwl_E0s]')
        $:GPU_UPDATE(device='[jwl_air_e0s, jwl_air_rho0s, jwl_air_gammas, jwl_idx]')
        $:GPU_UPDATE(device='[jwl_cv_prod, jwl_cv_air]')

    end subroutine s_initialize_jwl_module

    !> Deallocate the per-fluid JWL parameter tables.
    impure subroutine s_finalize_jwl_module

        @:DEALLOCATE(jwl_As, jwl_Bs, jwl_R1s, jwl_R2s, jwl_omegas, jwl_rho0s, jwl_E0s)
        @:DEALLOCATE(jwl_air_e0s, jwl_air_rho0s, jwl_air_gammas)

    end subroutine s_finalize_jwl_module

end module m_jwl
