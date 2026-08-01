!>
!! @file
!! @brief JWL EOS with a composition-weighted (heat-capacity) two-material mixture closure for ideal-gas and stiffened ambients.

#:include 'macros.fpp'
#:include 'case.fpp'

!> @brief Jones-Wilkins-Lee (JWL) EOS and its composition-weighted two-material closure for the five-equation model (Allaire et al.,
!! JCP 2002).
!!
!! Pure JWL products (Mie-Grueneisen referenced to an isentrope), V = rho0/rho:
!!   p = A(1 - w/(R1 V)) exp(-R1 V) + B(1 - w/(R2 V)) exp(-R2 V) + w rho e
!! A, B [Pa], R1, R2, w = omega [-] are cylinder-test fits; the first two terms are the
!! principal isentrope, the last the thermal (Grueneisen) pressure.
!! Refs: Lee/Hornig/Kury UCRL-50422 (1968); Menikoff LA-UR-15-29536 (2015).
!!
!! Two-material closure (products + ambient of mass fraction Y = alpha_rho_prod/rho).
!! One effective EOS per cell, heat-capacity weighted with w = Y cv_j/(Y cv_j + (1-Y) cv_a):
!!   An = w A, Bn = w B, omega = air_gamma + w(omega0 - air_gamma), cv = Y cv_j + (1-Y) cv_a.
!! Exact at Y=0 (ambient) and Y=1 (pure JWL). A stiffened ambient (pi_inf > 0) adds the
!! cold-stiffness offset pi_hat = (1-w) pi_inf; pi_inf = 0 recovers the ideal-gas closure
!! bit-identically. Every coefficient depends on Y alone, never rho or e.
!!
!! Because the coefficients are e-independent, the (rho, p, Y) -> e inverse is a single
!! closed form and the sound speed is the exact Grueneisen derivative
!! c^2 = (dp/drho)_e + (p/rho^2)(dp/de)_rho (no finite differencing, no e-region branches).
!! At start-up s_jwl_verify_closure sweeps the (rho, e, Y, lambda) envelope and aborts on any
!! non-positive/non-finite sound speed or failed p<->e round trip, so a bad fit fails fast.
!!
!! Stiffened-gas ambient: Le Metayer/Massoni/Saurel, Int. J. Therm. Sci. 43, 265 (2004).
!! Reaction sources (program burn, afterburn, JWL++) live in m_jwl_sources.
module m_jwl

    use m_global_parameters

    implicit none

    private
    public :: s_initialize_jwl_module, s_finalize_jwl_module, s_jwl_mix_state_er, s_jwl_mix_energy_pr, s_jwl_mix_sound_speed, &
        & s_jwl_mix_energy_sound_speed_pr, jwl_idx

    ! Simulation builds use m_global_parameters tables.
#ifndef MFC_SIMULATION
    real(wp), allocatable, public, dimension(:) :: jwl_As, jwl_Bs, jwl_R1s, jwl_R2s, jwl_omegas, jwl_rho0s, jwl_E0s
    real(wp), allocatable, public, dimension(:) :: jwl_air_e0s, jwl_air_rho0s, jwl_air_gammas, jwl_ej_rho_refs, jwl_air_pi_infs
    real(wp), allocatable, public, dimension(:) :: jwl_delta_es
    $:GPU_DECLARE(create='[jwl_As, jwl_Bs, jwl_R1s, jwl_R2s, jwl_omegas, jwl_rho0s, jwl_E0s]')
    $:GPU_DECLARE(create='[jwl_air_e0s, jwl_air_rho0s, jwl_air_gammas, jwl_ej_rho_refs, jwl_air_pi_infs]')
    $:GPU_DECLARE(create='[jwl_delta_es]')
#endif

    integer  :: jwl_idx                  !< JWL fluid index.
    real(wp) :: jwl_cv_prod, jwl_cv_air  !< Products/air specific heats.
    $:GPU_DECLARE(create='[jwl_idx]')
    $:GPU_DECLARE(create='[jwl_cv_prod, jwl_cv_air]')

    ! Fast PT-equilibrium closure constants (Jackson, JWL EOS notes June 2026, sec.10).
    ! jwl_pt_tau_gate: relative residual below which the composition-weighted warm start
    ! is already at PT equilibrium, so a single Newton correction suffices; above it the
    ! safeguarded bracketed solve is invoked. The CB->solve handoff jump is O(tau_gate) by
    ! construction, so the gate needs no material-specific tuning. jwl_pt_p_scale is the
    ! absolute pressure floor entering every relative tolerance, so tests neither loosen at
    ! p >> p_atm nor tighten to nothing near a stiffened-ambient cavitation state.
    real(wp), parameter :: jwl_pt_tau_gate = 1.e-6_wp
    real(wp), parameter :: jwl_pt_p_scale = 1.e5_wp
    ! Within jwl_pt_y_pure of a pure cell the composition-weighted closure IS exact PT, so
    ! the wrappers hand off to it there (keeping pure/single-material cells bit-identical to
    ! the pre-PT closure and keeping the interior solve away from the degenerate Y limits).
    real(wp), parameter :: jwl_pt_y_pure = 1.e-9_wp

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

    !> Effective mixture coefficients (composition-weighted; functions of Y only). Heat-capacity weighting lets omega relax
    !! air_gamma -> omega0 as products fill the cell, exact at Y=0,1. A stiffened ambient adds pi_hat = (1-w)*air_pi_inf. A future
    !! state-dependent closure (e.g. Jackson MG) would reintroduce rho/e dependence here, and with it coefficient derivatives in c2
    !! and e-region branches in the inverse -- both absent while An/Bn/omega are e-flat.
    subroutine s_jwl_weighted_composition_coeffs(Y, A, B, omega0, air_gamma, air_pi_inf, cv_j, cv_a, An, Bn, omega, cv, pi_c, &
        & pi_hat)

        $:GPU_ROUTINE(function_name='s_jwl_weighted_composition_coeffs',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: Y, A, B, omega0, air_gamma, air_pi_inf, cv_j, cv_a
        real(wp), intent(out) :: An, Bn, omega, cv, pi_c, pi_hat
        real(wp)              :: w

        w = Y*cv_j/(Y*cv_j + (1._wp - Y)*cv_a)
        An = w*A
        Bn = w*B
        omega = air_gamma + w*(omega0 - air_gamma)
        pi_hat = (1._wp - w)*air_pi_inf
        pi_c = (air_gamma + 1._wp)*pi_hat
        cv = Y*cv_j + (1._wp - Y)*cv_a  ! affects T only; p and c are cv-free

    end subroutine s_jwl_weighted_composition_coeffs

    !> Mixture state from energy: (rho, e, Y, lambda) -> (p, T, c2, c2 floor). lambda is the jwl_reactive progress (1 = fully
    !! reacted); delta_e is the reactant/product offset (0 = off). Only the thermal term uses e_eff = e + Y*(1-lambda)*delta_e,
    !! Y-scaled so pure ambient keeps its own energy; d(e_eff)/de = 1, so the inverse and c2 are exact.
    subroutine s_jwl_weighted_composition_state_er(rho, e, Y, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, &
        & lambda, delta_e, pres, T, c2, c2_floor)

        $:GPU_ROUTINE(function_name='s_jwl_weighted_composition_state_er',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, e, Y, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, lambda, delta_e
        real(wp), intent(out) :: pres, T, c2, c2_floor
        real(wp)              :: rho_s, Y_s, e_eff, An, Bn, omega, cv, pi_c, pi_hat, V, exp1, exp2, coef1, coef2

        rho_s = max(rho, sgm_eps)
        Y_s = min(max(Y, 0._wp), 1._wp)
        ! Offset belongs to the explosive, so scale by Y: pure ambient (Y=0) is untouched.
        e_eff = e + Y_s*(1._wp - min(max(lambda, 0._wp), 1._wp))*delta_e

        call s_jwl_weighted_composition_coeffs(Y_s, A, B, omega0, air_gamma, air_pi_inf, cv_j, cv_a, An, Bn, omega, cv, pi_c, &
                                               & pi_hat)

        V = rho0/rho_s
        exp1 = exp(-R1*V)
        exp2 = exp(-R2*V)
        coef1 = (1._wp - omega/(R1*V))*exp1
        coef2 = (1._wp - omega/(R2*V))*exp2

        pres = An*coef1 + Bn*coef2 + omega*rho_s*e_eff - pi_c
        T = (pres + pi_c - An*exp1 - Bn*exp2)/(omega*cv*rho_s)

        ! Frozen sound speed c2 = dp/drho|e + (p/rho^2) dp/de|rho (pi_c drops out).
        c2 = exp1*An*(R1*rho0/rho_s**2 - omega/rho_s - omega/(R1*rho0)) + exp2*Bn*(R2*rho0/rho_s**2 - omega/rho_s &
                      & - omega/(R2*rho0)) + omega*(e_eff + pres/rho_s)

        ! Raw c2 is returned so the init scan can catch non-positive values; c2_floor is the
        ! safety bound the wrappers apply (below any physical mixture c2, doubles as the
        ! stiffened-ambient cavitation cutoff).
        call s_jwl_floor(pres, sgm_eps)
        call s_jwl_floor(T, sgm_eps)
        c2_floor = min(air_gamma, merge(omega0, air_gamma, air_pi_inf > 0._wp))*(max(pres, sgm_eps) + pi_hat)/rho_s

    end subroutine s_jwl_weighted_composition_state_er

    !> Analytic inverse (rho, p, Y) -> e. Coefficients depend on Y alone, so this is a single closed form (no e-region structure);
    !! delta_e enters the pressure target as a constant.
    subroutine s_jwl_weighted_composition_energy_pr(rho, pres, Y, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, &
        & lambda, delta_e, e)

        $:GPU_ROUTINE(function_name='s_jwl_weighted_composition_energy_pr',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, pres, Y, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, lambda, delta_e
        real(wp), intent(out) :: e
        real(wp)              :: rho_s, Y_s, An, Bn, omega, cv, pi_c, pi_hat, V, C1, C2, de_shift

        rho_s = max(rho, sgm_eps)
        Y_s = min(max(Y, 0._wp), 1._wp)

        call s_jwl_weighted_composition_coeffs(Y_s, A, B, omega0, air_gamma, air_pi_inf, cv_j, cv_a, An, Bn, omega, cv, pi_c, &
                                               & pi_hat)
        V = rho0/rho_s
        C1 = (1._wp - omega/(R1*V))*exp(-R1*V)
        C2 = (1._wp - omega/(R2*V))*exp(-R2*V)

        ! Constant pressure-target shift from the thermal-term energy offset.
        de_shift = omega*rho_s*Y_s*(1._wp - min(max(lambda, 0._wp), 1._wp))*delta_e
        e = (pres + pi_c - An*C1 - Bn*C2 - de_shift)/max(omega*rho_s, sgm_eps)
        call s_jwl_floor(e, 0._wp)

    end subroutine s_jwl_weighted_composition_energy_pr

    !> Fused (rho, p, Y) -> (e, c): one coefficient/exponential pass yields both the energy inverse and the forward sound speed, so
    !! the Riemann path gets the reconstructed energy and the wave-speed sound speed from a single call. Expressions mirror
    !! s_jwl_weighted_composition_energy_pr and s_jwl_weighted_composition_state_er exactly.
    subroutine s_jwl_weighted_composition_sound_speed_pr(rho, pres, Y, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, &
        & cv_a, lambda, delta_e, e, c)

        $:GPU_ROUTINE(function_name='s_jwl_weighted_composition_sound_speed_pr',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, pres, Y, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, lambda, delta_e
        real(wp), intent(out) :: e, c
        real(wp)              :: rho_s, Y_s, e_eff, lambda_s, de_shift, An, Bn, omega, cv, pi_c, pi_hat
        real(wp)              :: V, exp1, exp2, C1, C2, p_m, cs2, cs2_floor

        rho_s = max(rho, sgm_eps)
        Y_s = min(max(Y, 0._wp), 1._wp)
        lambda_s = min(max(lambda, 0._wp), 1._wp)

        call s_jwl_weighted_composition_coeffs(Y_s, A, B, omega0, air_gamma, air_pi_inf, cv_j, cv_a, An, Bn, omega, cv, pi_c, &
                                               & pi_hat)
        V = rho0/rho_s
        exp1 = exp(-R1*V)
        exp2 = exp(-R2*V)
        C1 = (1._wp - omega/(R1*V))*exp1
        C2 = (1._wp - omega/(R2*V))*exp2

        ! Energy inverse, then forward pressure and sound speed (e -> e_eff in the thermal term).
        de_shift = omega*rho_s*Y_s*(1._wp - lambda_s)*delta_e
        e = (pres + pi_c - An*C1 - Bn*C2 - de_shift)/max(omega*rho_s, sgm_eps)
        call s_jwl_floor(e, 0._wp)
        e_eff = e + Y_s*(1._wp - lambda_s)*delta_e

        p_m = An*C1 + Bn*C2 + omega*rho_s*e_eff - pi_c
        cs2 = exp1*An*(R1*rho0/rho_s**2 - omega/rho_s - omega/(R1*rho0)) + exp2*Bn*(R2*rho0/rho_s**2 - omega/rho_s &
                       & - omega/(R2*rho0)) + omega*(e_eff + p_m/rho_s)
        call s_jwl_floor(p_m, sgm_eps)
        cs2_floor = min(air_gamma, merge(omega0, air_gamma, air_pi_inf > 0._wp))*(max(p_m, sgm_eps) + pi_hat)/rho_s
        c = sqrt(max(cs2, cs2_floor))

    end subroutine s_jwl_weighted_composition_sound_speed_pr

    ! Fast PT-equilibrium closure (Jackson, JWL EOS notes June 2026, sec.10). Products and
    ! ambient are held at a common (p, T); both are Grueneisen-caloric so T-equality and
    ! energy conservation are linear in T and solve explicitly, leaving one nonlinear scalar
    ! equation R(rho_p) = p_products - p_ambient = 0 in the products density. The
    ! composition-weighted state above is the warm start; a residual gate fires a Newton
    ! correction (or a safeguarded bracketed solve) only where the warm start is off
    ! equilibrium. The reactant offset delta_e is folded into the effective energy by the
    ! public wrappers, so these routines never see it (their argument e is already effective).
    ! Ported from the cross-validated standalone reference (jwl_standalone/src/jwl_newton.f90).

    !> JWL cold-curve energy e_cold(rho_p) and reference pressure pref (V = rho0/rho_p).
    subroutine s_jwl_pt_ecold(rho_p, A, B, R1, R2, rho0, ecold, pref)

        $:GPU_ROUTINE(function_name='s_jwl_pt_ecold',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho_p, A, B, R1, R2, rho0
        real(wp), intent(out) :: ecold, pref
        real(wp)              :: V, e1, e2

        V = rho0/rho_p
        e1 = exp(-R1*V)
        e2 = exp(-R2*V)
        pref = A*e1 + B*e2
        ecold = (A/R1*e1 + B/R2*e2)/rho0

    end subroutine s_jwl_pt_ecold

    !> Explicit common temperature T(rho_p) and dT/drho_p from thermal equilibrium and energy conservation. The exact JWL identity
    !! d(e_cold)/drho_p = pref/rho_p^2 makes dT analytic.
    subroutine s_jwl_pt_reduced_T(rho, e, Y, rho_p, A, B, R1, R2, rho0, air_pi_inf, cv_j, cv_a, T, dT)

        $:GPU_ROUTINE(function_name='s_jwl_pt_reduced_T',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, e, Y, rho_p, A, B, R1, R2, rho0, air_pi_inf, cv_j, cv_a
        real(wp), intent(out) :: T, dT
        real(wp)              :: cvm, ecd, pref

        cvm = Y*cv_j + (1._wp - Y)*cv_a
        call s_jwl_pt_ecold(rho_p, A, B, R1, R2, rho0, ecd, pref)
        T = (e - Y*ecd - air_pi_inf*(1._wp/rho - Y/rho_p))/cvm
        dT = -(Y/cvm)*(pref/rho_p**2 + air_pi_inf/rho_p**2)

    end subroutine s_jwl_pt_reduced_T

    !> Scalar residual R(rho_p) = p_products - p_ambient and its analytic derivative dR, returning the phase states (T, p_p, p_a,
    !! rho_a) for reuse. ok = .false. on an unphysical (non-positive ambient volume) trial.
    subroutine s_jwl_pt_residual(rho, e, Y, rho_p, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, R, dR, T, p_p, &
                                 & p_a, rho_a, ok)

        $:GPU_ROUTINE(function_name='s_jwl_pt_residual',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, e, Y, rho_p, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a
        real(wp), intent(out) :: R, dR, T, p_p, p_a, rho_a
        logical, intent(out)  :: ok
        real(wp)              :: u, dT, V, e1, e2, pref, dpref, drho_a

        ok = .true.
        u = 1._wp/rho - Y/rho_p  ! (1 - Y)/rho_a
        if (u <= sgm_eps .or. rho_p <= 0._wp) then
            ok = .false.
            R = huge(1._wp); dR = 1._wp
            T = 0._wp; p_p = 0._wp; p_a = 0._wp; rho_a = 0._wp
            return
        end if
        rho_a = (1._wp - Y)/u
        drho_a = -(Y/rho_p**2)*rho_a**2/(1._wp - Y)

        call s_jwl_pt_reduced_T(rho, e, Y, rho_p, A, B, R1, R2, rho0, air_pi_inf, cv_j, cv_a, T, dT)

        V = rho0/rho_p
        e1 = exp(-R1*V)
        e2 = exp(-R2*V)
        pref = A*e1 + B*e2
        dpref = (A*R1*e1 + B*R2*e2)*rho0/rho_p**2

        p_p = pref + omega0*rho_p*cv_j*T
        p_a = air_gamma*rho_a*cv_a*T - air_pi_inf

        R = p_p - p_a
        dR = dpref + omega0*cv_j*(T + rho_p*dT) - air_gamma*cv_a*(drho_a*T + rho_a*dT)

    end subroutine s_jwl_pt_residual

    !> Composition-weighted warm start backed out into a products density via the ambient EOS and volume additivity. Falls back to
    !! rho_p0 = rho outside the physical bracket.
    subroutine s_jwl_pt_cb_base(rho, e, Y, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, rho_p0)

        $:GPU_ROUTINE(function_name='s_jwl_pt_cb_base',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, e, Y, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a
        real(wp), intent(out) :: rho_p0
        real(wp)              :: p_cb, T_cb, c2_cb, c2f_cb, rho_a0, den

        ! lambda = 1, delta_e = 0: e is already the effective energy, so the CW routine must not re-apply the reactant offset.
        call s_jwl_weighted_composition_state_er(rho, e, Y, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, 1._wp, &
            & 0._wp, p_cb, T_cb, c2_cb, c2f_cb)
        rho_p0 = rho
        if (T_cb <= 0._wp .or. p_cb + air_pi_inf <= 0._wp) return
        rho_a0 = (p_cb + air_pi_inf)/(air_gamma*cv_a*T_cb)
        den = 1._wp/rho - (1._wp - Y)/rho_a0
        if (den <= sgm_eps) return
        rho_p0 = Y/den
        if (rho_p0 <= Y*rho*(1._wp + 1.e-12_wp)) rho_p0 = rho

    end subroutine s_jwl_pt_cb_base

    !> Analytic equilibrium sound speed by the implicit function theorem on R(rho_p) = 0: F_rp = dR/drho_p (analytic), drho_p/dq =
    !! -F_q/F_rp for q in {rho, e}, then c2 = dp/drho|_e + (p/rho^2) dp/de|_rho evaluated through the mixture pressure.
    subroutine s_jwl_pt_equilibrium_c2(rho, e, Y, rho_p, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, c2, ok)

        $:GPU_ROUTINE(function_name='s_jwl_pt_equilibrium_c2',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, e, Y, rho_p, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a
        real(wp), intent(out) :: c2
        logical, intent(out)  :: ok
        real(wp)              :: R, F_rp, T, p_p, p_a, rho_a, cvm
        real(wp)              :: dT_rp, dT_rho, dT_e, dra_rho
        real(wp)              :: dpp_rp, dpp_T, F_rho, F_e, drp_drho, drp_de
        real(wp)              :: V, e1, e2, pref, dpref, dp_drho, dp_de, pmix

        call s_jwl_pt_residual(rho, e, Y, rho_p, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, R, F_rp, T, p_p, &
                               & p_a, rho_a, ok)
        if (.not. ok .or. abs(F_rp) < 1.e-30_wp) then
            ok = .false.; c2 = 0._wp
            return
        end if
        cvm = Y*cv_j + (1._wp - Y)*cv_a

        V = rho0/rho_p
        e1 = exp(-R1*V)
        e2 = exp(-R2*V)
        pref = A*e1 + B*e2
        dpref = (A*R1*e1 + B*R2*e2)*rho0/rho_p**2

        dT_rp = -(Y/cvm)*(pref/rho_p**2 + air_pi_inf/rho_p**2)
        dT_rho = air_pi_inf/(cvm*rho**2)
        dT_e = 1._wp/cvm
        dra_rho = rho_a**2/((1._wp - Y)*rho**2)

        dpp_rp = dpref + omega0*cv_j*T  ! partial in rho_p at fixed T
        dpp_T = omega0*rho_p*cv_j

        ! F_rp (from the residual) is the TOTAL d/drho_p at fixed (rho, e).
        F_rho = dpp_T*dT_rho - air_gamma*cv_a*(dra_rho*T + rho_a*dT_rho)
        F_e = dpp_T*dT_e - air_gamma*cv_a*rho_a*dT_e

        drp_drho = -F_rho/F_rp
        drp_de = -F_e/F_rp

        dp_drho = dpp_rp*drp_drho + dpp_T*(dT_rp*drp_drho + dT_rho)
        dp_de = dpp_rp*drp_de + dpp_T*(dT_rp*drp_de + dT_e)

        pmix = (Y*rho/rho_p)*p_p + (1._wp - Y*rho/rho_p)*p_a
        c2 = dp_drho + (pmix/rho**2)*dp_de
        ok = .true.

    end subroutine s_jwl_pt_equilibrium_c2

    !> Safeguarded forward solve in x = ln(rho_p - Y*rho): for a stiffened ambient the root sits in a razor-thin layer just above
    !! rho_p = Y*rho and is ill-conditioned in rho_p, smooth in x. A uniform sign-change scan brackets the root nearest the CB warm
    !! start (deep-tension stiffened states admit a spurious second root), then Numerical-Recipes rtsafe converges on the bracket
    !! even where R is non-monotone. ierr: 1 no physical state, 4 no sign change, 3 cap without tolerance.
    subroutine s_jwl_pt_bracketed_solve(rho, e, Y, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, p, T, &
                                        & rho_p_out, ierr)

        $:GPU_ROUTINE(function_name='s_jwl_pt_bracketed_solve',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, e, Y, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a
        real(wp), intent(out) :: p, T, rho_p_out
        integer, intent(out)  :: ierr
        real(wp)              :: Rlo, rho_p, R, dR, p_p, p_a, rho_a
        real(wp)              :: Tl, dTl, cand, scale, x, xlo, xhi, dRdx, rho_p_lo, dxold
        real(wp)              :: relR, relbest, xbest, xl, xh, dxnew, x_cb, xprev, Rprev, bestd
        integer               :: stall, it, k
        logical               :: found, ok
        integer, parameter    :: itmax = 60, nscan = 200
        real(wp), parameter   :: tol = 1.e-11_wp, u_floor = 1.e-12_wp

        ierr = 0
        p = 0._wp; T = 0._wp; rho_p_out = rho

        if (1._wp/rho - u_floor <= 0._wp) then
            ierr = 1
            return
        end if
        rho_p_lo = Y/(1._wp/rho - u_floor)
        xlo = log(rho_p_lo - Y*rho)
        xhi = log(10._wp*rho0 - Y*rho)

        ! Upper bound must keep T > 0: shrink onto the T = 0 boundary in x.
        call s_jwl_pt_reduced_T(rho, e, Y, Y*rho + exp(xhi), A, B, R1, R2, rho0, air_pi_inf, cv_j, cv_a, Tl, dTl)
        if (Tl <= 0._wp) then
            cand = xlo
            do k = 1, 100
                x = 0.5_wp*(cand + xhi)
                call s_jwl_pt_reduced_T(rho, e, Y, Y*rho + exp(x), A, B, R1, R2, rho0, air_pi_inf, cv_j, cv_a, Tl, dTl)
                if (Tl > 0._wp) then
                    cand = x
                else
                    xhi = x
                end if
            end do
            xhi = cand
        end if
        call s_jwl_pt_reduced_T(rho, e, Y, Y*rho + exp(xlo), A, B, R1, R2, rho0, air_pi_inf, cv_j, cv_a, Tl, dTl)
        if (Tl <= 0._wp .or. xhi <= xlo) then
            ierr = 1
            return
        end if

        ! CB warm start in x selects the thermodynamically correct root (nearest x_cb).
        call s_jwl_pt_cb_base(rho, e, Y, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, rho_p)
        x_cb = 0.5_wp*(xlo + xhi)
        if (rho_p > Y*rho) then
            cand = log(rho_p - Y*rho)
            if (cand > xlo .and. cand < xhi) x_cb = cand
        end if

        ! Uniform scan in x; keep the sign-change bracket nearest x_cb.
        found = .false.; bestd = huge(1._wp)
        xprev = xlo
        call s_jwl_pt_residual(rho, e, Y, Y*rho + exp(xlo), A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, Rprev, &
                               & dR, T, p_p, p_a, rho_a, ok)
        do k = 1, nscan
            x = xlo + (xhi - xlo)*real(k, wp)/real(nscan, wp)
            call s_jwl_pt_residual(rho, e, Y, Y*rho + exp(x), A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, R, &
                                   & dR, T, p_p, p_a, rho_a, ok)
            if (ok .and. R*Rprev < 0._wp) then
                cand = 0.5_wp*(xprev + x)
                if (abs(cand - x_cb) < bestd) then
                    bestd = abs(cand - x_cb)
                    xlo = xprev; Rlo = Rprev; xhi = x
                    found = .true.
                end if
            end if
            xprev = x; Rprev = R
        end do
        if (.not. found) then
            ierr = 4
            return
        end if

        ! Faithful Numerical-Recipes rtsafe on [xl, xh], oriented R(xl) < 0 < R(xh).
        if (Rlo < 0._wp) then
            xl = xlo; xh = xhi
        else
            xl = xhi; xh = xlo
        end if
        x = 0.5_wp*(xl + xh)
        if (x_cb > min(xl, xh) .and. x_cb < max(xl, xh)) x = x_cb
        dxold = abs(xh - xl); dxnew = dxold
        call s_jwl_pt_residual(rho, e, Y, Y*rho + exp(x), A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, R, dR, &
                               & T, p_p, p_a, rho_a, ok)
        dRdx = dR*exp(x)
        xbest = x; relbest = huge(1._wp); stall = 0
        do it = 1, itmax
            ! bisect if Newton would leave [xl, xh] or is not decreasing fast enough
            if (((x - xh)*dRdx - R)*((x - xl)*dRdx - R) > 0._wp .or. abs(2._wp*R) > abs(dxold*dRdx)) then
                dxold = dxnew
                dxnew = 0.5_wp*(xh - xl)
                x = xl + dxnew
            else
                dxold = dxnew
                dxnew = R/dRdx
                x = x - dxnew
            end if
            rho_p = Y*rho + exp(x)
            call s_jwl_pt_residual(rho, e, Y, rho_p, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, R, dR, T, &
                                   & p_p, p_a, rho_a, ok)
            if (.not. ok) then  ! stepped out of the domain
                x = 0.5_wp*(xl + xh)
                call s_jwl_pt_residual(rho, e, Y, Y*rho + exp(x), A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, &
                                       & R, dR, T, p_p, p_a, rho_a, ok)
            end if
            dRdx = dR*exp(x)
            scale = max(abs(p_p), abs(p_a), jwl_pt_p_scale)
            relR = abs(R)/scale
            if (relR < relbest) then
                relbest = relR; xbest = x; stall = 0
            else
                stall = stall + 1
            end if
            if (relR < tol) exit
            ! shrink the bracket with the fresh sign
            if (R < 0._wp) then
                xl = x
            else
                xh = x
            end if
            if (abs(xh - xl) < 4.e-15_wp*(1._wp + abs(x))) exit  ! resolved in x
            ! Residual stagnation on an already-tight bracket is the stiffened-ambient
            ! floating-point cancellation floor, not a failure: p is accurate to that floor.
            if (stall >= 8 .and. abs(xh - xl) < 1.e-6_wp) exit
        end do

        ! Return the best iterate seen.
        x = xbest
        rho_p = Y*rho + exp(x)
        call s_jwl_pt_residual(rho, e, Y, rho_p, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, R, dR, T, p_p, &
                               & p_a, rho_a, ok)
        if (.not. ok .or. T <= 0._wp) then
            ierr = 1
            return
        end if
        if (abs(R)/max(abs(p_p), abs(p_a), jwl_pt_p_scale) < 1.e-6_wp) then
            ierr = 0
        else
            ierr = max(ierr, 3)
        end if
        p = (Y*rho/rho_p)*p_p + (1._wp - Y*rho/rho_p)*p_a
        rho_p_out = rho_p

    end subroutine s_jwl_pt_bracketed_solve

    !> Fast PT forward evaluator (interior 0 < Y < 1): composition-weighted warm start, exact solve gated on the warm-start
    !! residual. Every failure path falls back to a CB estimate, so it always returns a finite (p, T, c2) with a positive (floored)
    !! c2.
    subroutine s_jwl_pt_state(rho, e, Y, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, p, T, c2)

        $:GPU_ROUTINE(function_name='s_jwl_pt_state',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, e, Y, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a
        real(wp), intent(out) :: p, T, c2
        real(wp)              :: p_cb, T_cb, c2_cb, c2f_cb, rho_p0, R0, dR, Tr, p_p, p_a, rho_a
        real(wp)              :: rho_p, c2_eq, Tr2, p_p2, p_a2, rho_a2, c2_cw
        logical               :: ok
        integer               :: ierr

        ! Warm start (also the fallback). lambda = 1, delta_e = 0: e is already effective.
        call s_jwl_weighted_composition_state_er(rho, e, Y, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, 1._wp, &
            & 0._wp, p_cb, T_cb, c2_cb, c2f_cb)
        c2_cw = max(c2_cb, c2f_cb)
        p = p_cb; T = T_cb
        ! Positivity floor (pi_hat = (1 - Y)*pi_inf) so every return path yields a real sound
        ! speed for the Riemann solver, including CB returns at tension states.
        c2 = max(c2_cw, min(air_gamma, omega0)*max(p_cb + (1._wp - Y)*air_pi_inf, jwl_pt_p_scale)/rho)

        call s_jwl_pt_cb_base(rho, e, Y, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, rho_p0)
        call s_jwl_pt_residual(rho, e, Y, rho_p0, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, R0, dR, Tr, p_p, &
                               & p_a, rho_a, ok)
        ! Gate on the residual relative to the true equilibrium pressure p_p (not p_cb, whose
        ! magnitude is set by pi_inf for a stiffened ambient). At equilibrium the two partial
        ! pressures agree and the mixture pressure is that common value.
        if (ok .and. abs(R0) <= jwl_pt_tau_gate*max(abs(p_p), jwl_pt_p_scale)) then
            ! One Newton correction makes the returned pressure O(tau_gate^2) accurate; keep
            ! the pre-step state if the step leaves the physical bracket.
            p = 0.5_wp*(p_p + p_a)
            rho_p = rho_p0 - R0/dR
            call s_jwl_pt_residual(rho, e, Y, rho_p, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, R0, dR, Tr2, &
                                   & p_p2, p_a2, rho_a2, ok)
            if (ok) then
                rho_p0 = rho_p
                p = 0.5_wp*(p_p2 + p_a2)
                T = Tr2
            else
                T = Tr
            end if
            call s_jwl_pt_equilibrium_c2(rho, e, Y, rho_p0, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, c2_eq, &
                                         & ok)
            if (ok .and. c2_eq > 0._wp) then
                c2 = c2_eq
            else
                c2 = min(air_gamma, omega0)*max(p + (1._wp - Y)*air_pi_inf, jwl_pt_p_scale)/rho
            end if
            return
        end if

        call s_jwl_pt_bracketed_solve(rho, e, Y, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, p, T, rho_p, ierr)
        ! ierr 3 (converged to ~1e-6, not 1e-11) is still a genuine equilibrium; only ierr 1/4
        ! (no physical state / no sign change) or a non-positive T mean no PT equilibrium exists.
        if ((ierr == 0 .or. ierr == 3) .and. T > 0._wp) then
            call s_jwl_pt_equilibrium_c2(rho, e, Y, rho_p, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, c2_eq, ok)
            if (ok .and. c2_eq > 0._wp) then
                c2 = c2_eq
            else
                c2 = max(c2_cw, min(air_gamma, omega0)*max(p + (1._wp - Y)*air_pi_inf, jwl_pt_p_scale)/rho)
            end if
            return
        end if
        ! No PT equilibrium (e below the products' cold-curve floor, an unphysical corner):
        ! return the composition-weighted state, which is always thermodynamically valid with a
        ! positive temperature (T feeds burn-rate laws downstream and must not go negative).
        p = p_cb; T = T_cb
        c2 = max(c2_cw, min(air_gamma, omega0)*max(p_cb + (1._wp - Y)*air_pi_inf, jwl_pt_p_scale)/rho)

    end subroutine s_jwl_pt_state

    !> Inverse scalar residual R_inv(rho_p) = p_products - p_target with T fixed by the ambient law at the target pressure, and its
    !! analytic derivative.
    subroutine s_jwl_pt_residual_inv(rho, p, Y, rho_p, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, R, dR, T, ok)

        $:GPU_ROUTINE(function_name='s_jwl_pt_residual_inv',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, p, Y, rho_p, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a
        real(wp), intent(out) :: R, dR, T
        logical, intent(out)  :: ok
        real(wp)              :: u, rho_a, drho_a, dT, V, e1, e2, pref, dpref

        ok = .true.
        u = 1._wp/rho - Y/rho_p
        if (u <= sgm_eps .or. rho_p <= 0._wp) then
            ok = .false.; R = huge(1._wp); dR = 1._wp; T = 0._wp
            return
        end if
        rho_a = (1._wp - Y)/u
        drho_a = -(Y/rho_p**2)*rho_a**2/(1._wp - Y)
        T = (p + air_pi_inf)/(air_gamma*rho_a*cv_a)
        dT = -T*drho_a/rho_a
        V = rho0/rho_p
        e1 = exp(-R1*V)
        e2 = exp(-R2*V)
        pref = A*e1 + B*e2
        dpref = (A*R1*e1 + B*R2*e2)*rho0/rho_p**2
        R = pref + omega0*rho_p*cv_j*T - p
        dR = dpref + omega0*cv_j*(T + rho_p*dT)

    end subroutine s_jwl_pt_residual_inv

    !> Robust inverse solve (rho, p, Y) -> e by the same rtsafe machinery as the forward bracketed solve, in x = ln(rho_p - Y*rho).
    !! ierr /= 0 signals the caller to fall back to the closed-form composition-weighted inverse.
    subroutine s_jwl_pt_bracketed_inverse(rho, p, Y, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, e, ierr)

        $:GPU_ROUTINE(function_name='s_jwl_pt_bracketed_inverse',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, p, Y, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a
        real(wp), intent(out) :: e
        integer, intent(out)  :: ierr
        real(wp)              :: rho_p, rho_a, u, T, R, dR, dRdx, x, xl, xh, xlo, xhi
        real(wp)              :: rho_p_lo, cand, dxold, dxnew, Rlo, Rhi, scale, relR, relbest, xbest, ecd, pref
        integer               :: it, k, stall
        logical               :: ok
        integer, parameter    :: itmax = 100, nscan = 128
        real(wp), parameter   :: tol = 1.e-12_wp, u_floor = 1.e-12_wp

        ierr = 0
        e = 0._wp
        if (1._wp/rho - u_floor <= 0._wp) then
            ierr = 1
            return
        end if

        rho_p_lo = Y/(1._wp/rho - u_floor)
        xlo = log(rho_p_lo - Y*rho)
        xhi = log(10._wp*rho0 - Y*rho)
        call s_jwl_pt_residual_inv(rho, p, Y, Y*rho + exp(xlo), A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, &
                                   & Rlo, dR, T, ok)
        if (.not. ok) then
            ierr = 1
            return
        end if
        call s_jwl_pt_residual_inv(rho, p, Y, Y*rho + exp(xhi), A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, &
                                   & Rhi, dR, T, ok)
        if (.not. ok) then
            ierr = 1
            return
        end if
        if (Rlo*Rhi > 0._wp) then
            cand = xlo
            do k = 1, nscan
                x = xlo + (xhi - xlo)*real(k, wp)/real(nscan, wp)
                call s_jwl_pt_residual_inv(rho, p, Y, Y*rho + exp(x), A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, &
                                           & cv_a, R, dR, T, ok)
                if (ok .and. R*Rlo < 0._wp) then
                    xhi = x; Rhi = R; xlo = cand
                    exit
                end if
                cand = x; Rlo = R
                if (k == nscan) then  ! no sign change: fall back
                    ierr = 1
                    return
                end if
            end do
        end if

        if (Rlo < 0._wp) then
            xl = xlo; xh = xhi
        else
            xl = xhi; xh = xlo
        end if
        x = 0.5_wp*(xl + xh)
        call s_jwl_pt_residual_inv(rho, p, Y, Y*rho + exp(x), A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, R, &
                                   & dR, T, ok)
        dRdx = dR*exp(x); dxold = abs(xh - xl); dxnew = dxold
        xbest = x; relbest = huge(1._wp); stall = 0
        do it = 1, itmax
            if (((x - xh)*dRdx - R)*((x - xl)*dRdx - R) > 0._wp .or. abs(2._wp*R) > abs(dxold*dRdx)) then
                dxold = dxnew; dxnew = 0.5_wp*(xh - xl); x = xl + dxnew
            else
                dxold = dxnew; dxnew = R/dRdx; x = x - dxnew
            end if
            rho_p = Y*rho + exp(x)
            call s_jwl_pt_residual_inv(rho, p, Y, rho_p, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, R, dR, T, &
                                       & ok)
            if (.not. ok) then
                x = 0.5_wp*(xl + xh)
                call s_jwl_pt_residual_inv(rho, p, Y, Y*rho + exp(x), A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, &
                                           & cv_a, R, dR, T, ok)
            end if
            dRdx = dR*exp(x)
            scale = max(abs(p) + air_pi_inf, jwl_pt_p_scale)
            relR = abs(R)/scale
            if (relR < relbest) then
                relbest = relR; xbest = x; stall = 0
            else
                stall = stall + 1
            end if
            if (relR < tol) exit
            if (R < 0._wp) then
                xl = x
            else
                xh = x
            end if
            if (abs(xh - xl) < 4.e-15_wp*(1._wp + abs(x))) exit
            if (stall >= 8 .and. abs(xh - xl) < 1.e-6_wp) exit
        end do

        x = xbest
        rho_p = Y*rho + exp(x)
        u = 1._wp/rho - Y/rho_p
        if (u <= sgm_eps .or. T <= 0._wp) then
            ierr = 1
            return
        end if
        rho_a = (1._wp - Y)/u
        T = (p + air_pi_inf)/(air_gamma*rho_a*cv_a)
        call s_jwl_pt_residual_inv(rho, p, Y, rho_p, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, R, dR, T, ok)
        if (abs(R)/max(abs(p) + air_pi_inf, jwl_pt_p_scale) >= 1.e-6_wp) ierr = 3
        call s_jwl_pt_ecold(rho_p, A, B, R1, R2, rho0, ecd, pref)
        e = Y*(ecd + cv_j*T) + (1._wp - Y)*(cv_a*T + air_pi_inf/rho_a)

    end subroutine s_jwl_pt_bracketed_inverse

    !> Fast PT inverse (interior 0 < Y < 1): (rho, p, Y) -> effective energy e. The exact PT inverse when the solve succeeds, else
    !! the closed-form composition-weighted inverse.
    subroutine s_jwl_pt_inverse(rho, p, Y, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, e)

        $:GPU_ROUTINE(function_name='s_jwl_pt_inverse',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, p, Y, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a
        real(wp), intent(out) :: e
        integer               :: ierr

        call s_jwl_pt_bracketed_inverse(rho, p, Y, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, e, ierr)
        ! ierr 3 (converged to ~1e-6) is accepted so the inverse stays on the PT branch the
        ! forward used; only a hard bracketing failure (ierr 1) falls back to the CW inverse.
        if (ierr == 1) then
            call s_jwl_weighted_composition_energy_pr(rho, p, Y, A, B, R1, R2, omega0, rho0, air_gamma, air_pi_inf, cv_j, cv_a, &
                & 1._wp, 0._wp, e)
        end if

    end subroutine s_jwl_pt_inverse

    ! Public entry points: look up fluid jidx's parameters, then evaluate the closure.

    !> Full state from energy for fluid jidx: (rho, e, Y, [lambda]) -> (p, T, [c]). c is optional so pressure-only callers skip the
    !! sqrt. lambda (jwl_reactive reaction progress; 1 = fully reacted) defaults to 1, recovering the closure exactly for every
    !! caller that predates the reactant/product energy offset (jwl_delta_es(jidx) = 0 by default has the same effect regardless of
    !! lambda).
    subroutine s_jwl_mix_state_er(rho, e, Y, jidx, pres, T, c, lambda)

        $:GPU_ROUTINE(function_name='s_jwl_mix_state_er',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)            :: rho, e, Y
        integer, intent(in)             :: jidx
        real(wp), intent(out)           :: pres, T
        real(wp), intent(out), optional :: c
        real(wp), intent(in), optional  :: lambda
        real(wp)                        :: c2, c2_floor, lambda_l, Yc, e_eff

        lambda_l = 1._wp; if (present(lambda)) lambda_l = lambda
        Yc = min(max(Y, 0._wp), 1._wp)

        if (Yc <= jwl_pt_y_pure .or. Yc >= 1._wp - jwl_pt_y_pure) then
            ! Single-material limit: the CW routine is exact and bit-identical here, and keeps
            ! the reactant-offset/lambda handling in one place for pure-JWL cells.
            call s_jwl_weighted_composition_state_er(rho, e, Y, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), &
                & jwl_omegas(jidx), jwl_rho0s(jidx), jwl_air_gammas(jidx), jwl_air_pi_infs(jidx), jwl_cv_prod, jwl_cv_air, &
                & lambda_l, jwl_delta_es(jidx), pres, T, c2, c2_floor)
        else
            ! Fold the Garno reactant offset into the effective energy (enters reduced_T once,
            ! d/de = 1, so all PT derivatives are unchanged) and solve the PT equilibrium.
            e_eff = e + Yc*(1._wp - lambda_l)*jwl_delta_es(jidx)
            call s_jwl_pt_state(max(rho, sgm_eps), e_eff, Yc, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), &
                                & jwl_omegas(jidx), jwl_rho0s(jidx), jwl_air_gammas(jidx), jwl_air_pi_infs(jidx), jwl_cv_prod, &
                                & jwl_cv_air, pres, T, c2)
            c2_floor = sgm_eps
        end if
        ! Sound speed is optional: cons->prim callers that only need pressure skip the sqrt.
        ! The safety floor sits below any physical mixture c2, so it only engages on
        ! unphysical (e.g. cavitated) states without overriding legitimate sound speeds.
        if (present(c)) c = sqrt(max(c2, c2_floor))

    end subroutine s_jwl_mix_state_er

    !> Energy from pressure for fluid jidx: (rho, p, Y, [lambda]) -> e. lambda defaults to 1 (fully reacted).
    subroutine s_jwl_mix_energy_pr(rho, pres, Y, jidx, e, lambda)

        $:GPU_ROUTINE(function_name='s_jwl_mix_energy_pr',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)           :: rho, pres, Y
        integer, intent(in)            :: jidx
        real(wp), intent(out)          :: e
        real(wp), intent(in), optional :: lambda
        real(wp)                       :: lambda_l, Yc, rho_c, e_eff

        lambda_l = 1._wp; if (present(lambda)) lambda_l = lambda
        Yc = min(max(Y, 0._wp), 1._wp)
        rho_c = max(rho, sgm_eps)

        if (Yc <= jwl_pt_y_pure .or. Yc >= 1._wp - jwl_pt_y_pure) then
            call s_jwl_weighted_composition_energy_pr(rho_c, pres, Yc, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), &
                & jwl_omegas(jidx), jwl_rho0s(jidx), jwl_air_gammas(jidx), jwl_air_pi_infs(jidx), jwl_cv_prod, jwl_cv_air, &
                & lambda_l, jwl_delta_es(jidx), e)
        else
            call s_jwl_pt_inverse(rho_c, pres, Yc, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), jwl_omegas(jidx), &
                                  & jwl_rho0s(jidx), jwl_air_gammas(jidx), jwl_air_pi_infs(jidx), jwl_cv_prod, jwl_cv_air, e_eff)
            ! Remove the reactant offset the forward step folded into the effective energy.
            e = e_eff - Yc*(1._wp - lambda_l)*jwl_delta_es(jidx)
        end if

    end subroutine s_jwl_mix_energy_pr

    !> Sound speed for fluid jidx: invert to energy, then evaluate c. lambda defaults to 1 (fully reacted).
    subroutine s_jwl_mix_sound_speed(rho, pres, Y, jidx, c, lambda)

        $:GPU_ROUTINE(function_name='s_jwl_mix_sound_speed',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)           :: rho, pres, Y
        integer, intent(in)            :: jidx
        real(wp), intent(out)          :: c
        real(wp), intent(in), optional :: lambda
        real(wp)                       :: lambda_l, Yc, rho_c, e_eff, p_out, T_out, c2

        lambda_l = 1._wp; if (present(lambda)) lambda_l = lambda
        Yc = min(max(Y, 0._wp), 1._wp)
        rho_c = max(rho, sgm_eps)

        if (Yc <= jwl_pt_y_pure .or. Yc >= 1._wp - jwl_pt_y_pure) then
            call s_jwl_weighted_composition_sound_speed_pr(rho, pres, Y, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), &
                & jwl_R2s(jidx), jwl_omegas(jidx), jwl_rho0s(jidx), jwl_air_gammas(jidx), jwl_air_pi_infs(jidx), jwl_cv_prod, &
                & jwl_cv_air, lambda_l, jwl_delta_es(jidx), e_eff, c)
        else
            ! Invert to the effective energy, then take the equilibrium c2 from a forward solve
            ! (the reactant offset cancels: it shifts e_eff back in, so the state is the same).
            call s_jwl_pt_inverse(rho_c, pres, Yc, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), jwl_omegas(jidx), &
                                  & jwl_rho0s(jidx), jwl_air_gammas(jidx), jwl_air_pi_infs(jidx), jwl_cv_prod, jwl_cv_air, e_eff)
            call s_jwl_pt_state(rho_c, e_eff, Yc, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), jwl_omegas(jidx), &
                                & jwl_rho0s(jidx), jwl_air_gammas(jidx), jwl_air_pi_infs(jidx), jwl_cv_prod, jwl_cv_air, p_out, &
                                & T_out, c2)
            c = sqrt(max(c2, sgm_eps))
        end if

    end subroutine s_jwl_mix_sound_speed

    !> Fused energy and sound speed for fluid jidx: (rho, p, Y, [lambda]) -> (e, c). The Riemann faces need both the reconstructed
    !! energy and a sound speed, so a single call shares the one coefficient/exponential pass instead of inverting for e and again
    !! for c. lambda defaults to 1 (fully reacted).
    subroutine s_jwl_mix_energy_sound_speed_pr(rho, pres, Y, jidx, e, c, lambda)

        $:GPU_ROUTINE(function_name='s_jwl_mix_energy_sound_speed_pr',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)           :: rho, pres, Y
        integer, intent(in)            :: jidx
        real(wp), intent(out)          :: e, c
        real(wp), intent(in), optional :: lambda
        real(wp)                       :: lambda_l, Yc, rho_c, e_eff, p_out, T_out, c2

        lambda_l = 1._wp; if (present(lambda)) lambda_l = lambda
        Yc = min(max(Y, 0._wp), 1._wp)
        rho_c = max(rho, sgm_eps)

        if (Yc <= jwl_pt_y_pure .or. Yc >= 1._wp - jwl_pt_y_pure) then
            call s_jwl_weighted_composition_sound_speed_pr(rho, pres, Y, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), &
                & jwl_R2s(jidx), jwl_omegas(jidx), jwl_rho0s(jidx), jwl_air_gammas(jidx), jwl_air_pi_infs(jidx), jwl_cv_prod, &
                & jwl_cv_air, lambda_l, jwl_delta_es(jidx), e, c)
        else
            call s_jwl_pt_inverse(rho_c, pres, Yc, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), jwl_omegas(jidx), &
                                  & jwl_rho0s(jidx), jwl_air_gammas(jidx), jwl_air_pi_infs(jidx), jwl_cv_prod, jwl_cv_air, e_eff)
            call s_jwl_pt_state(rho_c, e_eff, Yc, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), jwl_omegas(jidx), &
                                & jwl_rho0s(jidx), jwl_air_gammas(jidx), jwl_air_pi_infs(jidx), jwl_cv_prod, jwl_cv_air, p_out, &
                                & T_out, c2)
            e = e_eff - Yc*(1._wp - lambda_l)*jwl_delta_es(jidx)
            c = sqrt(max(c2, sgm_eps))
        end if

    end subroutine s_jwl_mix_energy_sound_speed_pr

    !> Initialize JWL parameter tables.
    impure subroutine s_initialize_jwl_module

        use m_mpi_common, only: s_mpi_abort
        use m_helper_basic, only: f_approx_equal, f_is_default

        integer  :: i, n_jwl, n_air, air_idx, gamma_src
        real(wp) :: jwl_E0_from_Q, air_e0_from_p0

        @:ALLOCATE(jwl_As(1:num_fluids), jwl_Bs(1:num_fluids), jwl_R1s(1:num_fluids), jwl_R2s(1:num_fluids), &
                   & jwl_omegas(1:num_fluids), jwl_rho0s(1:num_fluids), jwl_E0s(1:num_fluids))
        @:ALLOCATE(jwl_air_e0s(1:num_fluids), jwl_air_rho0s(1:num_fluids), jwl_air_gammas(1:num_fluids), &
                   & jwl_ej_rho_refs(1:num_fluids), jwl_air_pi_infs(1:num_fluids))
        @:ALLOCATE(jwl_delta_es(1:num_fluids))

        jwl_idx = 0
        n_jwl = 0
        n_air = 0
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
            jwl_ej_rho_refs(i) = fluid_pp(i)%jwl_ej_rho_ref
            jwl_delta_es(i) = 0._wp
            if (.not. f_is_default(fluid_pp(i)%jwl_delta_e)) jwl_delta_es(i) = fluid_pp(i)%jwl_delta_e
            if (fluid_pp(i)%eos == eos_jwl) then
                jwl_idx = i
                n_jwl = n_jwl + 1
                if (f_is_default(fluid_pp(i)%jwl_A) .or. f_is_default(fluid_pp(i)%jwl_B) .or. f_is_default(fluid_pp(i)%jwl_R1) &
                    & .or. f_is_default(fluid_pp(i)%jwl_R2) .or. f_is_default(fluid_pp(i)%jwl_omega) &
                    & .or. f_is_default(fluid_pp(i)%jwl_rho0) .or. f_is_default(fluid_pp(i)%jwl_E0)) then
                    call s_mpi_abort('fluid_pp%eos = eos_jwl requires jwl_A, jwl_B, jwl_R1, jwl_R2, ' &
                                     & // 'jwl_omega, jwl_rho0, and either jwl_Q or jwl_E0 to be set.')
                end if
                if (.not. f_is_default(fluid_pp(i)%jwl_Q)) then
                    jwl_E0_from_Q = fluid_pp(i)%jwl_rho0*fluid_pp(i)%jwl_Q
                    if (.not. f_approx_equal(fluid_pp(i)%jwl_E0, jwl_E0_from_Q, 1.e-8_wp)) then
                        call s_mpi_abort('fluid_pp%eos = eos_jwl requires jwl_E0 = jwl_rho0*jwl_Q when both jwl_E0 and jwl_Q are set.')
                    end if
                end if
                if (f_is_default(fluid_pp(i)%jwl_air_rho0)) then
                    call s_mpi_abort('fluid_pp%eos = eos_jwl requires jwl_air_rho0 to be set.')
                end if
                if (f_is_default(fluid_pp(i)%jwl_air_e0) .and. f_is_default(fluid_pp(i)%jwl_air_p0)) then
                    call s_mpi_abort('fluid_pp%eos = eos_jwl requires either jwl_air_e0 or jwl_air_p0 to be set.')
                end if
                if (fluid_pp(i)%jwl_R1 <= 0._wp .or. fluid_pp(i)%jwl_R2 <= 0._wp .or. fluid_pp(i)%jwl_omega <= 0._wp &
                    & .or. fluid_pp(i)%jwl_rho0 <= 0._wp .or. fluid_pp(i)%jwl_E0 <= 0._wp .or. fluid_pp(i)%jwl_air_rho0 <= 0._wp) &
                    & then
                    call s_mpi_abort('JWL parameters jwl_R1, jwl_R2, jwl_omega, jwl_rho0, jwl_Q/jwl_E0, ' &
                                     & // 'and jwl_air_rho0 must be positive.')
                end if
            else
                n_air = n_air + 1
                if (air_idx == 0) air_idx = i
            end if
        end do

        if (n_jwl > 1) then
            call s_mpi_abort('At most one fluid may use eos_jwl; found more than one.')
        end if

        if (jwl_idx > 0 .and. model_eqns /= model_eqns_5eq) then
            call s_mpi_abort('eos_jwl is only supported with model_eqns_5eq.')
        end if

        jwl_cv_prod = 0._wp
        if (jwl_idx > 0) jwl_cv_prod = fluid_pp(jwl_idx)%cv
        if (air_idx > 0) then
            jwl_cv_air = fluid_pp(air_idx)%cv
        else
            jwl_cv_air = jwl_cv_prod  ! No ambient fluid: mass-weighted cv degenerates to cv_prod.
        end if

        if (jwl_idx > 0) then
            if (f_is_default(jwl_cv_prod) .or. jwl_cv_prod <= 0._wp) then
                call s_mpi_abort('The weighted-composition closure requires positive fluid_pp%cv for the JWL fluid.')
            end if
            if (num_fluids > 1 .and. n_air /= 1) then
                call s_mpi_abort('The weighted-composition closure requires exactly one non-JWL ideal-gas fluid.')
            end if
            if (air_idx > 0) then
                if (f_is_default(fluid_pp(air_idx)%cv) .or. fluid_pp(air_idx)%cv <= 0._wp) then
                    call s_mpi_abort('The weighted-composition closure requires positive fluid_pp%cv for the non-JWL air fluid.')
                end if
            end if

            ! Ambient Grueneisen coefficient Gamma = 1/fluid_pp%gamma, from the non-JWL
            ! fluid (or the JWL fluid's own gamma when it is alone).
            if (air_idx > 0) then
                gamma_src = air_idx
            else
                gamma_src = jwl_idx
            end if
            if (f_is_default(fluid_pp(gamma_src)%gamma) .or. fluid_pp(gamma_src)%gamma <= 0._wp) then
                call s_mpi_abort('The weighted-composition closure requires positive fluid_pp%gamma for the ambient-gas Grueneisen coefficient.')
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

            ! Products-energy reference density e_j = jwl_E0/jwl_ej_rho_ref (default jwl_rho0).
            if (f_is_default(fluid_pp(jwl_idx)%jwl_ej_rho_ref)) jwl_ej_rho_refs(jwl_idx) = jwl_rho0s(jwl_idx)
            if (jwl_ej_rho_refs(jwl_idx) <= 0._wp) then
                call s_mpi_abort('fluid_pp%jwl_ej_rho_ref must be positive.')
            end if

            if (jwl_rho0s(jwl_idx) <= jwl_air_rho0s(jwl_idx) .or. jwl_E0s(jwl_idx)/jwl_ej_rho_refs(jwl_idx) &
                & <= jwl_air_e0s(jwl_idx)) then
                call s_mpi_abort('The weighted-composition closure requires increasing air-to-products reference density and energy.')
            end if

            ! Verify the assembled closure is positive-definite and invertible over the
            ! physical envelope before it is used by the solver.
            call s_jwl_verify_closure(jwl_idx)
        end if

        $:GPU_UPDATE(device='[jwl_As, jwl_Bs, jwl_R1s, jwl_R2s, jwl_omegas, jwl_rho0s, jwl_E0s]')
        $:GPU_UPDATE(device='[jwl_air_e0s, jwl_air_rho0s, jwl_air_gammas, jwl_ej_rho_refs, jwl_air_pi_infs, jwl_idx]')
        $:GPU_UPDATE(device='[jwl_delta_es, jwl_cv_prod, jwl_cv_air]')

    end subroutine s_initialize_jwl_module

    !> Init-time self-check: sweep the (rho, e, Y) envelope and abort if the assembled closure yields a non-positive sound speed or
    !! fails the (rho, p, Y) -> e -> p inverse.
    impure subroutine s_jwl_verify_closure(jidx)

        use m_mpi_common, only: s_mpi_abort

        integer, intent(in) :: jidx
        integer             :: ir, ie, iy, il
        real(wp)            :: rho_lo, rho_hi, rho_hi_l, e_lo, e_hi, ej, rho_s, e_s, pres, T, c, e_inv
        real(wp)            :: e_eff, ecd_floor, pref_floor
        character(len=128)  :: msg
        integer, parameter  :: n_scan = 100
        ! Stiffened-ambient PT states carry a floating-point cancellation floor near the
        ! rho_p = Y*rho layer; matches the reference closure's own 1e-6 acceptance.
        real(wp), parameter :: scan_rtol = 1.e-6_wp
        real(wp), parameter :: y_scan(8) = [0._wp, 0.25_wp, 0.5_wp, 0.75_wp, 0.9_wp, 0.97_wp, 0.999_wp, 1._wp]
        real(wp), parameter :: lambda_scan(3) = [0._wp, 0.5_wp, 1._wp]

        ej = jwl_E0s(jidx)/jwl_ej_rho_refs(jidx)
        ! Envelope covers reflected-shock states beyond the CJ point; the density cap is
        ! tighter for a stiffened ambient (JWL is only meaningful to ~2*rho0 there).
        rho_lo = 0.1_wp*jwl_air_rho0s(jidx)
        rho_hi = 4._wp*jwl_rho0s(jidx)
        if (jwl_air_pi_infs(jidx) > 0._wp) rho_hi = 2._wp*jwl_rho0s(jidx)
        e_lo = 0.5_wp*jwl_air_e0s(jidx)
        e_hi = 5._wp*ej

        do il = 1, size(lambda_scan)
            ! The reactant/product offset is only applied to UNREACTED explosive (lambda < 1),
            ! which physically exists only up to the von Neumann compression (~2*rho0). Scanning
            ! the offset at the full 4*rho0 reflected-shock cap tests an unreachable state whose
            ! large negative e_eff can spuriously fail the c2 check, so cap the reactant sweep.
            rho_hi_l = rho_hi
            if (jwl_delta_es(jidx) /= 0._wp .and. lambda_scan(il) < 1._wp) rho_hi_l = min(rho_hi, 2._wp*jwl_rho0s(jidx))
            do iy = 1, size(y_scan)
                do ir = 0, n_scan - 1
                    rho_s = rho_lo*(rho_hi_l/rho_lo)**(real(ir, wp)/real(n_scan - 1, wp))
                    do ie = 0, n_scan - 1
                        e_s = e_lo + (e_hi - e_lo)*real(ie, wp)/real(n_scan - 1, wp)
                        ! A PT equilibrium exists only where the effective mixture energy clears
                        ! the products' cold-curve energy at maximal admissible expansion
                        ! (rho_p -> Y*rho+); below that floor no positive-temperature equilibrium
                        ! exists (products cannot be colder than their own cold curve), so the
                        ! state is outside the physical envelope and is not required to invert.
                        ! The old CW closure only appeared to invert it by evaluating JWL at the
                        ! bulk density; the PT closure correctly rejects it.
                        if (y_scan(iy) > jwl_pt_y_pure) then
                            e_eff = e_s + y_scan(iy)*(1._wp - lambda_scan(il))*jwl_delta_es(jidx)
                            call s_jwl_pt_ecold(y_scan(iy)*rho_s, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), &
                                                & jwl_rho0s(jidx), ecd_floor, pref_floor)
                            if (e_eff <= y_scan(iy)*ecd_floor) cycle
                        end if
                        ! Exercise the closure through the public wrappers so the self-check
                        ! validates the exact PT path the solver uses (interior Y) and the CW
                        ! endpoints, including the floored sound speed the Riemann solver sees.
                        call s_jwl_mix_state_er(rho_s, e_s, y_scan(iy), jidx, pres, T, c, lambda_scan(il))
                        ! Floored sound speed must always be real and positive.
                        if (c <= 0._wp .or. c /= c) then
                            write (msg, &
                                   & '(A,ES11.4,A,ES11.4,A,F6.3,A,F6.3)') &
                                   & 'JWL closure self-check: non-positive sound speed at rho=', rho_s, ', e=', e_s, ', Y=', &
                                   & y_scan(iy), ', lambda=', lambda_scan(il)
                            call s_mpi_abort(trim(msg) // '. Check JWL and ambient parameters.')
                        end if
                        if (pres > sgm_eps) then
                            call s_jwl_mix_energy_pr(rho_s, pres, y_scan(iy), jidx, e_inv, lambda_scan(il))
                            if (e_inv /= e_inv .or. abs(e_inv - e_s) > scan_rtol*max(abs(e_s), jwl_air_e0s(jidx))) then
                                write (msg, &
                                       & '(A,ES11.4,A,ES11.4,A,F6.3,A,F6.3)') &
                                       & 'JWL closure self-check: energy inversion mismatch at rho=', rho_s, ', e=', e_s, ', Y=', &
                                       & y_scan(iy), ', lambda=', lambda_scan(il)
                                call s_mpi_abort(trim(msg) // '. Check JWL and ambient-gas parameters.')
                            end if
                        end if
                    end do
                end do
            end do
        end do

    end subroutine s_jwl_verify_closure

    !> Deallocate the per-fluid JWL parameter tables.
    impure subroutine s_finalize_jwl_module

        @:DEALLOCATE(jwl_As, jwl_Bs, jwl_R1s, jwl_R2s, jwl_omegas, jwl_rho0s, jwl_E0s)
        @:DEALLOCATE(jwl_air_e0s, jwl_air_rho0s, jwl_air_gammas, jwl_ej_rho_refs, jwl_air_pi_infs)
        @:DEALLOCATE(jwl_delta_es)

    end subroutine s_finalize_jwl_module

end module m_jwl
