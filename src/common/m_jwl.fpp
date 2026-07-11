!>
!! @file
!! @brief JWL EOS with a composition-weighted (heat-capacity) two-material mixture closure for ideal-gas and stiffened ambients.

#:include 'macros.fpp'
#:include 'case.fpp'

!> @brief Jones-Wilkins-Lee (JWL) equation of state and the composition-weighted two-material closure for the five-equation model.
!!
!! 1. Pure JWL products
!! The JWL EOS describes the detonation products of a condensed explosive as a Mie-Grueneisen
!! solid referenced to an isentrope. With relative volume V = rho0/rho (rho0 the products
!! reference density) and specific internal energy e, the pressure is
!!
!!     p(rho, e) = A (1 - w/(R1 V)) exp(-R1 V)          <- cold (isentrope) term, coeff. A
!!               + B (1 - w/(R2 V)) exp(-R2 V)          <- cold (isentrope) term, coeff. B
!!               + w rho e.                             <- thermal (Grueneisen) term
!!
!! Here A, B [Pa], R1, R2 [-], and the Grueneisen coefficient w = omega [-] are fitted to
!! cylinder-test expansion data. The first two terms are the principal isentrope p_s(V); the
!! last is the thermal pressure w rho e, with w playing the role of the Grueneisen parameter
!! Gamma = V (dp/de)_V. Reference: E. L. Lee, H. C. Hornig, J. W. Kury, "Adiabatic Expansion
!! of High Explosive Detonation Products," UCRL-50422, Lawrence Radiation Lab. (1968); see
!! also R. Menikoff, "JWL Equation of State," LA-UR-15-29536, Los Alamos (2015).
!!
!! The temperature follows from the caloric relation e = p_s(V)/(w rho) + cv T, i.e.
!!     T(rho, e) = (p - A exp(-R1 V) - B exp(-R2 V)) / (w cv rho),
!! and the frozen sound speed is the exact thermodynamic derivative
!!     c^2 = (dp/drho)_e + (p/rho^2)(dp/de)_rho,
!! evaluated in closed form in s_jwl_rocflu_state_er (no finite differencing).
!!
!! 2. Two-material closure (products + ambient)
!! In a five-equation (Allaire et al., JCP 2002) simulation a cell may hold a mixture of
!! products and an ambient fluid. Rather than solve a full pressure-temperature equilibrium,
!! MFC assembles one effective EOS per cell as a function of the state (rho, e) and the products
!! mass fraction Y = (alpha rho)_products / rho (see s_jwl_rocflu_coeffs). For an ideal-gas ambient
!! the coefficients are composition (heat-capacity) weighted: with the products' heat-capacity
!! share w = Y cv_j /(Y cv_j + (1-Y) cv_a), the amplitudes are A_n = w A, B_n = w B and the
!! Grueneisen coefficient is omega = air_gamma + w (omega0 - air_gamma). Weighting by composition
!! rather than density keeps omega relaxing toward omega0 as products fill the cell (the density
!! ramp never reaches rho0 in afterburn mixing), removing the Rocflu p/c overshoot; the closure is
!! exact at Y=0 and Y=1 and matches full p-T equilibrium to ~1e-8 in the weak-pressure regime.
!! For a stiffened ambient (pi_inf > 0) the same weight w drives a cold-stiffness offset
!! pi_hat = (1-w) pi (=> pi_c = (air_gamma+1) pi_hat); with pi = 0 this reduces bit-identically to
!! the ideal-gas closure. Every coefficient varies smoothly with Y (and is independent of rho and
!! e), so the closure is continuous in mass fraction across 0 <= Y <= 1.
!!
!! 3. Stiffened-gas ambient (underwater / condensed)
!! When the ambient fluid is a stiffened gas (e.g. water), p_ambient = Gamma rho e - (Gamma+1) pi
!! with stiffness pi > 0. The closure carries this as a cold-stiffness offset pi_c that is
!! independent of both rho and e, so the analytic pressure->energy inverse and the Grueneisen
!! sound-speed identity are preserved. Stiffened-gas EOS: O. Le Metayer, J. Massoni, R. Saurel,
!! Int. J. Thermal Sciences 43, 265-276 (2004).
!!
!! 4. Analytic inverse and self-verification
!! Because A_n(e), B_n(e) are piecewise linear in e, the map (rho, p, Y) -> e inverts in closed
!! form (s_jwl_rocflu_energy_pr): evaluate the coefficients at the mid-ramp energy to obtain the
!! exact linear slope, then correct the low- and high-energy saturated branches. At start-up the
!! module sweeps the (rho, e, Y) envelope and aborts if any state gives a non-positive or
!! non-finite sound speed, or if the p<->e round trip fails, so a bad parameter set fails fast.
!!
!! Pressure-driven and afterburn reaction sources built on this EOS live in m_jwl_sources.
!! Full design notes and validation: README-JWL-EOS.md.
module m_jwl

    use m_global_parameters

    implicit none

    private
    public :: s_initialize_jwl_module, s_finalize_jwl_module, s_jwl_mix_state_er, s_jwl_mix_energy_pr, s_jwl_mix_sound_speed, &
        & jwl_idx

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

    ! Composition-weighted closure

    !> Effective mixture coefficients. Composition (heat-capacity) weighted closure -- An/Bn = w*A, w*B and omega = air_gamma +
    !! w*(omega0 - air_gamma) with w = Y*cv_j/(Y*cv_j + (1-Y)*cv_a), all independent of rho and e (mA = mB = momega = 0), so the
    !! analytic (rho, p, Y) -> e inverse and the closed-form sound speed stay exact. A stiffened ambient (air_pi_inf > 0) adds the
    !! cold-stiffness offset pi_hat = (1-w)*air_pi_inf (=> pi_c = (air_gamma+1)*pi_hat); air_pi_inf = 0 recovers the ideal-gas
    !! closure bit-identically. cv is mass-weighted; A_sat/B_sat = An/Bn (no separate Region-III inverse now that An/Bn are e-flat).
    subroutine s_jwl_rocflu_coeffs(rho, e, Y, A, B, omega0, rho0, E0, ej_rho_ref, air_e0, air_rho0, air_gamma, air_pi_inf, cv_j, &
                                   & cv_a, An, Bn, omega, cv, mA, mB, momega, pi_c, pi_hat, A_sat, B_sat)

        $:GPU_ROUTINE(function_name='s_jwl_rocflu_coeffs',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, e, Y, A, B, omega0, rho0, E0, ej_rho_ref, air_e0, air_rho0, air_gamma, air_pi_inf, cv_j, cv_a
        real(wp), intent(out) :: An, Bn, omega, cv, mA, mB, momega, pi_c, pi_hat, A_sat, B_sat
        real(wp)              :: w

        ! Weighting by heat-capacity share rather than density lets omega relax
        ! air_gamma -> omega0 as products fill the cell (the density ramp never
        ! reached rho0 in afterburn mixing, the source of the Rocflu p/c overshoot).
        w = Y*cv_j/(Y*cv_j + (1._wp - Y)*cv_a)
        An = w*A
        Bn = w*B
        mA = 0._wp
        mB = 0._wp
        omega = air_gamma + w*(omega0 - air_gamma)
        momega = 0._wp
        pi_hat = (1._wp - w)*air_pi_inf
        pi_c = (air_gamma + 1._wp)*pi_hat
        A_sat = An
        B_sat = Bn

        ! Mass-weighted heat capacity (affects T only; p and c are cv-free).
        cv = Y*cv_j + (1._wp - Y)*cv_a

    end subroutine s_jwl_rocflu_coeffs

    !> Rocflu single-fluid state-interpolated closure: (rho, e, Y, lambda) -> (p, T, c², c² floor). lambda is the jwl_reactive
    !! reaction progress (1 = fully reacted); delta_e is the reactant/product energy offset (0 disables the effect, recovering the
    !! fully-reacted closure exactly). Only the thermal term uses the shifted e_eff = e + Y*(1-lambda)*delta_e (Y-scaled so pure
    !! ambient keeps its own energy), so d(e_eff)/de = 1 and the sound-speed algebra below needs no new derivative terms beyond
    !! replacing e with e_eff.
    subroutine s_jwl_rocflu_state_er(rho, e, Y, A, B, R1, R2, omega0, rho0, E0, ej_rho_ref, air_e0, air_rho0, air_gamma, &
                                     & air_pi_inf, cv_j, cv_a, lambda, delta_e, pres, T, c2, c2_floor)

        $:GPU_ROUTINE(function_name='s_jwl_rocflu_state_er',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in) :: rho, e, Y, A, B, R1, R2, omega0, rho0, E0, ej_rho_ref, air_e0, air_rho0, air_gamma, air_pi_inf, &
             & cv_j, cv_a, lambda, delta_e
        real(wp), intent(out) :: pres, T, c2, c2_floor
        real(wp) :: rho_s, Y_s, e_eff, An, Bn, omega, cv, mA, mB, momega, pi_c, pi_hat, A_sat, B_sat, V, exp1, exp2, coef1, coef2

        rho_s = max(rho, sgm_eps)
        Y_s = min(max(Y, 0._wp), 1._wp)
        ! The reactant/product offset belongs to the explosive only, so it is scaled by
        ! the JWL mass fraction Y: pure ambient (Y=0) keeps its own energy untouched, and
        ! the offset fades out as products mix into air. d(e_eff)/de = 1 still holds
        ! (Y, lambda, delta_e are all independent of e), so the closed-form inverse and
        ! sound speed are preserved.
        e_eff = e + Y_s*(1._wp - min(max(lambda, 0._wp), 1._wp))*delta_e

        ! Effective coefficients. The blend decays exactly to the ambient law as
        ! An, Bn -> 0, so no separate pure-ambient branch is needed. An/Bn ramp in the
        ! real e (not e_eff): the Y-vs-ambient ramp (is this JWL material) and the
        ! lambda-vs-reactant offset (has the explosive reacted) are distinct axes.
        call s_jwl_rocflu_coeffs(rho_s, e, Y_s, A, B, omega0, rho0, E0, ej_rho_ref, air_e0, air_rho0, air_gamma, air_pi_inf, &
                                 & cv_j, cv_a, An, Bn, omega, cv, mA, mB, momega, pi_c, pi_hat, A_sat, B_sat)

        ! Cold-curve exponentials at relative volume V = rho0/rho.
        V = rho0/rho_s
        exp1 = exp(-R1*V)
        exp2 = exp(-R2*V)
        coef1 = (1._wp - omega/(R1*V))*exp1
        coef2 = (1._wp - omega/(R2*V))*exp2

        ! Pressure and temperature. pi_c is independent of rho and e in both branches.
        pres = An*coef1 + Bn*coef2 + omega*rho_s*e_eff - pi_c
        T = (pres + pi_c - An*exp1 - Bn*exp2)/(omega*cv*rho_s)

        ! Sound speed c2 = dp/drho|e + (p/rho^2) dp/de|rho, with the pi_c term dropping out.
        ! d(e_eff)/de = 1 at fixed lambda, so dp/de|rho is unchanged (still omega); only the
        ! two standalone thermal-energy terms in dp/drho|e pick up e -> e_eff.
        c2 = exp1*An*(R1*rho0/rho_s**2 - omega/rho_s - omega/(R1*rho0) - rho_s*momega/(R1*rho0)) + mA*pres*coef1/rho_s**2 &
                      & + exp2*Bn*(R2*rho0/rho_s**2 - omega/rho_s - omega/(R2*rho0) - rho_s*momega/(R2*rho0)) &
                      & + mB*pres*coef2/rho_s**2 + omega*(e_eff + pres/rho_s) + momega*rho_s*e_eff

        ! Pressure floor doubles as the cavitation cutoff for stiffened ambients. c2 is
        ! returned raw so the init self-check can catch non-positive values; c2_floor is
        ! the safety bound the wrappers apply (below any physical mixture c2).
        call s_jwl_floor(pres, sgm_eps)
        call s_jwl_floor(T, sgm_eps)
        c2_floor = min(air_gamma, merge(omega0, air_gamma, air_pi_inf > 0._wp))*(max(pres, sgm_eps) + pi_hat)/rho_s

    end subroutine s_jwl_rocflu_state_er

    !> Analytic inverse (rho, p, Y) -> e. An(e)/Bn(e) are piecewise linear, so evaluating coefficients at the Region-II midpoint
    !! gives an exact inverse; Regions I/III correct the saturated offsets.

    subroutine s_jwl_rocflu_energy_pr(rho, pres, Y, A, B, R1, R2, omega0, rho0, E0, ej_rho_ref, air_e0, air_rho0, air_gamma, &
                                      & air_pi_inf, cv_j, cv_a, lambda, delta_e, e)

        $:GPU_ROUTINE(function_name='s_jwl_rocflu_energy_pr',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in) :: rho, pres, Y, A, B, R1, R2, omega0, rho0, E0, ej_rho_ref, air_e0, air_rho0, air_gamma, &
             & air_pi_inf, cv_j, cv_a, lambda, delta_e
        real(wp), intent(out) :: e
        real(wp) :: rho_s, Y_s, e_j, e_eval, An, Bn, omega, cv, mA, mB, momega, pi_c, pi_hat, A_sat, B_sat, V, C1, C2, de_shift

        rho_s = max(rho, sgm_eps)
        Y_s = min(max(Y, 0._wp), 1._wp)

        e_j = E0/ej_rho_ref
        e_eval = 0.5_wp*(e_j + air_e0)
        call s_jwl_rocflu_coeffs(rho_s, e_eval, Y_s, A, B, omega0, rho0, E0, ej_rho_ref, air_e0, air_rho0, air_gamma, air_pi_inf, &
                                 & cv_j, cv_a, An, Bn, omega, cv, mA, mB, momega, pi_c, pi_hat, A_sat, B_sat)
        V = rho0/rho_s
        C1 = (1._wp - omega/(R1*V))*exp(-R1*V)
        C2 = (1._wp - omega/(R2*V))*exp(-R2*V)

        ! Pressure target shift from the thermal-term energy offset: the forward law replaces omega*rho*e with
        ! omega*rho*(e + Y*(1-lambda)*delta_e), so every branch below solves for the real e by subtracting this constant
        ! (omega, Y, lambda do not depend on e, so the same shift applies to all three regions).
        de_shift = omega*rho_s*Y_s*(1._wp - min(max(lambda, 0._wp), 1._wp))*delta_e

        ! pi_c is independent of e, so adding it back to pres keeps every branch exact.
        ! Region II: exact linear inverse for the active coefficient branch.

        e = (pres + pi_c + (mA*C1 + mB*C2)*e_eval - An*C1 - Bn*C2 - de_shift)/max(mA*C1 + mB*C2 + omega*rho_s, sgm_eps)
        if (e < air_e0) then
            ! Region I: subtract the low-energy coefficient offsets, which are nonzero in the pure-JWL branch.
            e = (pres + pi_c - (An - mA*(e_eval - air_e0))*C1 - (Bn - mB*(e_eval - air_e0))*C2 - de_shift)/max(omega*rho_s, sgm_eps)
        else if (e > e_j) then
            ! Region III: saturated coefficients -> p = A_sat*C1 + B_sat*C2 + omega*rho*e - pi_c.
            e = (pres + pi_c - A_sat*C1 - B_sat*C2 - de_shift)/max(omega*rho_s, sgm_eps)
        end if
        call s_jwl_floor(e, 0._wp)

    end subroutine s_jwl_rocflu_energy_pr

    !> Fused (rho, p, Y) -> c: one coefficient/exponential evaluation shared between the energy inverse and the forward sound speed
    !! (the Riemann path calls this three times per face, so avoiding the second closure pass roughly halves the JWL EOS cost).
    !! Expressions mirror s_jwl_rocflu_energy_pr and s_jwl_rocflu_state_er exactly.

    subroutine s_jwl_rocflu_sound_speed_pr(rho, pres, Y, A, B, R1, R2, omega0, rho0, E0, ej_rho_ref, air_e0, air_rho0, air_gamma, &
                                           & air_pi_inf, cv_j, cv_a, lambda, delta_e, c)

        $:GPU_ROUTINE(function_name='s_jwl_rocflu_sound_speed_pr',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in) :: rho, pres, Y, A, B, R1, R2, omega0, rho0, E0, ej_rho_ref, air_e0, air_rho0, air_gamma, &
             & air_pi_inf, cv_j, cv_a, lambda, delta_e
        real(wp), intent(out) :: c
        real(wp) :: rho_s, Y_s, e_j, e_eval, e, e_eff, lambda_s, de_shift, An, Bn, omega, cv, mA, mB, momega, pi_c, pi_hat, &
             & A_sat, B_sat
        real(wp) :: V, exp1, exp2, C1, C2, p_m, cs2, cs2_floor

        rho_s = max(rho, sgm_eps)
        Y_s = min(max(Y, 0._wp), 1._wp)
        lambda_s = min(max(lambda, 0._wp), 1._wp)
        e_j = E0/ej_rho_ref

        ! Coefficients at the Region-II midpoint and the cold-curve exponentials, evaluated
        ! once and shared by the inverse and the forward sound speed below.
        e_eval = 0.5_wp*(e_j + air_e0)
        call s_jwl_rocflu_coeffs(rho_s, e_eval, Y_s, A, B, omega0, rho0, E0, ej_rho_ref, air_e0, air_rho0, air_gamma, air_pi_inf, &
                                 & cv_j, cv_a, An, Bn, omega, cv, mA, mB, momega, pi_c, pi_hat, A_sat, B_sat)
        V = rho0/rho_s
        exp1 = exp(-R1*V)
        exp2 = exp(-R2*V)
        C1 = (1._wp - omega/(R1*V))*exp1
        C2 = (1._wp - omega/(R2*V))*exp2

        ! Energy inverse (same branches and delta_e correction as s_jwl_rocflu_energy_pr).
        de_shift = omega*rho_s*Y_s*(1._wp - lambda_s)*delta_e
        e = (pres + pi_c + (mA*C1 + mB*C2)*e_eval - An*C1 - Bn*C2 - de_shift)/max(mA*C1 + mB*C2 + omega*rho_s, sgm_eps)
        if (e < air_e0) then
            e = (pres + pi_c - (An - mA*(e_eval - air_e0))*C1 - (Bn - mB*(e_eval - air_e0))*C2 - de_shift)/max(omega*rho_s, sgm_eps)
        else if (e > e_j) then
            e = (pres + pi_c - A_sat*C1 - B_sat*C2 - de_shift)/max(omega*rho_s, sgm_eps)
        end if
        call s_jwl_floor(e, 0._wp)
        e_eff = e + Y_s*(1._wp - lambda_s)*delta_e

        ! Coefficients are composition-weighted (independent of rho and e), so the values from
        ! s_jwl_rocflu_coeffs above are already final -- no re-blend at the recovered energy.

        ! Forward pressure and sound speed (same expressions as s_jwl_rocflu_state_er, e -> e_eff in the thermal term).
        p_m = An*C1 + Bn*C2 + omega*rho_s*e_eff - pi_c
        cs2 = exp1*An*(R1*rho0/rho_s**2 - omega/rho_s - omega/(R1*rho0) - rho_s*momega/(R1*rho0)) + mA*p_m*C1/rho_s**2 &
                       & + exp2*Bn*(R2*rho0/rho_s**2 - omega/rho_s - omega/(R2*rho0) - rho_s*momega/(R2*rho0)) &
                       & + mB*p_m*C2/rho_s**2 + omega*(e_eff + p_m/rho_s) + momega*rho_s*e_eff
        call s_jwl_floor(p_m, sgm_eps)
        cs2_floor = min(air_gamma, merge(omega0, air_gamma, air_pi_inf > 0._wp))*(max(p_m, sgm_eps) + pi_hat)/rho_s
        c = sqrt(max(cs2, cs2_floor))

    end subroutine s_jwl_rocflu_sound_speed_pr

    ! Public entry points: look up fluid jidx's parameters, then evaluate the closure.

    !> Full state from energy for fluid jidx: (rho, e, Y, [lambda]) -> (p, T, c). lambda (jwl_reactive reaction progress; 1 = fully
    !! reacted) defaults to 1, recovering the closure exactly for every caller that predates the reactant/product energy offset
    !! (jwl_delta_es(jidx) = 0 by default has the same effect regardless of lambda).
    subroutine s_jwl_mix_state_er(rho, e, Y, jidx, pres, T, c, lambda)

        $:GPU_ROUTINE(function_name='s_jwl_mix_state_er',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)           :: rho, e, Y
        integer, intent(in)            :: jidx
        real(wp), intent(out)          :: pres, T, c
        real(wp), intent(in), optional :: lambda
        real(wp)                       :: c2, c2_floor, lambda_l

        lambda_l = 1._wp; if (present(lambda)) lambda_l = lambda

        call s_jwl_rocflu_state_er(rho, e, Y, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), jwl_omegas(jidx), &
                                   & jwl_rho0s(jidx), jwl_E0s(jidx), jwl_ej_rho_refs(jidx), jwl_air_e0s(jidx), &
                                   & jwl_air_rho0s(jidx), jwl_air_gammas(jidx), jwl_air_pi_infs(jidx), jwl_cv_prod, jwl_cv_air, &
                                   & lambda_l, jwl_delta_es(jidx), pres, T, c2, c2_floor)
        ! The safety floor sits below any physical mixture c2, so it only engages on
        ! unphysical (e.g. cavitated) states without overriding legitimate sound speeds.
        c = sqrt(max(c2, c2_floor))

    end subroutine s_jwl_mix_state_er

    !> Energy from pressure for fluid jidx: (rho, p, Y, [lambda]) -> e. lambda defaults to 1 (fully reacted).
    subroutine s_jwl_mix_energy_pr(rho, pres, Y, jidx, e, lambda)

        $:GPU_ROUTINE(function_name='s_jwl_mix_energy_pr',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)           :: rho, pres, Y
        integer, intent(in)            :: jidx
        real(wp), intent(out)          :: e
        real(wp), intent(in), optional :: lambda
        real(wp)                       :: lambda_l

        lambda_l = 1._wp; if (present(lambda)) lambda_l = lambda

        call s_jwl_rocflu_energy_pr(max(rho, sgm_eps), pres, min(max(Y, 0._wp), 1._wp), jwl_As(jidx), jwl_Bs(jidx), &
                                    & jwl_R1s(jidx), jwl_R2s(jidx), jwl_omegas(jidx), jwl_rho0s(jidx), jwl_E0s(jidx), &
                                    & jwl_ej_rho_refs(jidx), jwl_air_e0s(jidx), jwl_air_rho0s(jidx), jwl_air_gammas(jidx), &
                                    & jwl_air_pi_infs(jidx), jwl_cv_prod, jwl_cv_air, lambda_l, jwl_delta_es(jidx), e)

    end subroutine s_jwl_mix_energy_pr

    !> Sound speed for fluid jidx: invert to energy, then evaluate c. lambda defaults to 1 (fully reacted).
    subroutine s_jwl_mix_sound_speed(rho, pres, Y, jidx, c, lambda)

        $:GPU_ROUTINE(function_name='s_jwl_mix_sound_speed',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)           :: rho, pres, Y
        integer, intent(in)            :: jidx
        real(wp), intent(out)          :: c
        real(wp), intent(in), optional :: lambda
        real(wp)                       :: lambda_l

        lambda_l = 1._wp; if (present(lambda)) lambda_l = lambda

        call s_jwl_rocflu_sound_speed_pr(rho, pres, Y, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), &
                                         & jwl_omegas(jidx), jwl_rho0s(jidx), jwl_E0s(jidx), jwl_ej_rho_refs(jidx), &
                                         & jwl_air_e0s(jidx), jwl_air_rho0s(jidx), jwl_air_gammas(jidx), jwl_air_pi_infs(jidx), &
                                         & jwl_cv_prod, jwl_cv_air, lambda_l, jwl_delta_es(jidx), c)

    end subroutine s_jwl_mix_sound_speed

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
            jwl_delta_es(i) = fluid_pp(i)%jwl_delta_e
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
                call s_mpi_abort('The Rocflu closure requires positive fluid_pp%cv for the JWL fluid.')
            end if
            if (num_fluids > 1 .and. n_air /= 1) then
                call s_mpi_abort('The Rocflu closure requires exactly one non-JWL ideal-gas fluid.')
            end if
            if (air_idx > 0) then
                if (f_is_default(fluid_pp(air_idx)%cv) .or. fluid_pp(air_idx)%cv <= 0._wp) then
                    call s_mpi_abort('The Rocflu closure requires positive fluid_pp%cv for the non-JWL air fluid.')
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
                call s_mpi_abort('The Rocflu closure requires positive fluid_pp%gamma for the ambient-gas Grueneisen coefficient.')
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
                call s_mpi_abort('The Rocflu closure requires increasing air-to-products reference density and energy.')
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
        real(wp)            :: rho_lo, rho_hi, rho_hi_l, e_lo, e_hi, ej, rho_s, e_s, pres, T, c2, c2_floor, e_inv
        character(len=128)  :: msg
        integer, parameter  :: n_scan = 100
        real(wp), parameter :: scan_rtol = 1.e-8_wp
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
                        call s_jwl_rocflu_state_er(rho_s, e_s, y_scan(iy), jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), &
                                                   & jwl_R2s(jidx), jwl_omegas(jidx), jwl_rho0s(jidx), jwl_E0s(jidx), &
                                                   & jwl_ej_rho_refs(jidx), jwl_air_e0s(jidx), jwl_air_rho0s(jidx), &
                                                   & jwl_air_gammas(jidx), jwl_air_pi_infs(jidx), jwl_cv_prod, jwl_cv_air, &
                                                   & lambda_scan(il), jwl_delta_es(jidx), pres, T, c2, c2_floor)
                        ! Floored c2 must always be positive; raw c2 must be positive wherever
                        ! the pressure was not floored (floored = cavitated, handled by cutoff).
                        if (max(c2, c2_floor) <= 0._wp .or. c2_floor /= c2_floor) then
                            write (msg, &
                                   & '(A,ES11.4,A,ES11.4,A,F6.3,A,F6.3)') &
                                   & 'JWL closure self-check: non-positive floored sound speed at rho=', rho_s, ', e=', e_s, &
                                   & ', Y=', y_scan(iy), ', lambda=', lambda_scan(il)
                            call s_mpi_abort(trim(msg) // '. Check JWL and ambient parameters.')
                        end if
                        if (pres > sgm_eps .and. (c2 <= 0._wp .or. c2 /= c2)) then
                            write (msg, &
                                   & '(A,ES11.4,A,ES11.4,A,F6.3,A,F6.3)') &
                                   & 'JWL closure self-check: non-positive sound speed at rho=', rho_s, ', e=', e_s, ', Y=', &
                                   & y_scan(iy), ', lambda=', lambda_scan(il)
                            call s_mpi_abort(trim(msg) // '. Check JWL and ambient parameters.')
                        end if
                        if (pres > sgm_eps) then
                            call s_jwl_rocflu_energy_pr(rho_s, pres, y_scan(iy), jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), &
                                                        & jwl_R2s(jidx), jwl_omegas(jidx), jwl_rho0s(jidx), jwl_E0s(jidx), &
                                                        & jwl_ej_rho_refs(jidx), jwl_air_e0s(jidx), jwl_air_rho0s(jidx), &
                                                        & jwl_air_gammas(jidx), jwl_air_pi_infs(jidx), jwl_cv_prod, jwl_cv_air, &
                                                        & lambda_scan(il), jwl_delta_es(jidx), e_inv)
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
