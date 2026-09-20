!>
!! @file
!! @brief Contains module m_eos

#:include 'macros.fpp'
#:include 'case.fpp'
#:include 'generated_eos.fpp'

!> @brief Equations of state in Gamma/Pi form, rho e = Gamma(rho) p + Pi(rho).
!!
!! Stiffened and ideal gas keep constant coefficients, resolved once at start-up. The
!! state-dependent families (Mie-Gruneisen, JWL, Vinet) evaluate theirs per cell from a
!! reference curve.
!!
!! The mixture closure rules that combine the phases -- Wood's law, the six-equation mean, the
!! bubbly branch -- are not themselves equations of state, but s_compute_mixture_coefficients,
!! s_compute_speed_of_sound and their _dt/_avg variants live here rather than in
!! m_variables_conversion, and must stay here. They are the hot-path callers of the phase
!! chain (s_phase_coefficients -> s_eos_coefficients -> s_reference_curve), and on NVHPC that
!! inlining only happens within a single file: the cross-file inliner refuses any device routine
!! with a subroutine call in its call tree. Splitting them from the chain costs ~25% of grind
!! time on NVHPC and nothing on the other backends, so it fails quietly. The solver kernels never
!! inlined these four in the first place, which is why the module boundary is drawn above them and
!! not below. See docs/documentation/gpuParallelization.md, "Module boundaries and NVHPC inlining".
!!
!! This module is a leaf: it directly uses only m_derived_types, m_constants, and
!! m_global_parameters_common. Adding an EOS family means one case in s_reference_curve.
module m_eos

    use m_derived_types
    use m_global_parameters_common
    use m_mpi_common, only: s_prohibit_abort
    use m_constants, only: eos_ideal_gas, eos_mie_gruneisen, eos_jwl, eos_vinet, eos_rk4_steps, ode_isentrope, &
        & ode_reference_temperature, sgm_eps, verysmall, dflt_real

    implicit none

    private

    public :: s_compute_mixture_coefficients, s_compute_mixture_coefficients_dt, s_compute_speed_of_sound, &
        & s_compute_speed_of_sound_avg, s_initialize_eos_module, s_finalize_eos_module, f_pressure, f_bulk_modulus, &
        & f_relativistic_enthalpy, f_isentrope_exponent, f_isentrope_pressure, f_sg_thermal, f_mixture_temperature, &
        & f_is_state_dependent, s_phase_coefficients, s_phase_pressure_on_isentrope, s_phase_temperature, &
        & s_phase_density_on_isentrope, s_phase_internal_energy, s_phase_bulk_modulus

contains

    !> Resolve every fluid's EOS coefficients once, before any conversion runs.
    impure subroutine s_initialize_eos_module()

        integer :: i
        logical :: state_dependent  !< Whether this case's fluids need a density-dependent EOS

        @:ALLOCATE(gammas (1:num_fluids))
        @:ALLOCATE(eoss (1:num_fluids))
        @:ALLOCATE(isentrope_n (1:num_fluids))
        @:ALLOCATE(pi_infs(1:num_fluids))
        @:ALLOCATE(isentrope_B(1:num_fluids))
        @:ALLOCATE(cvs    (1:num_fluids))
        @:ALLOCATE(qvs    (1:num_fluids))
        @:ALLOCATE(qvps    (1:num_fluids))

        state_dependent = .false.
        do i = 1, num_fluids
            gammas(i) = fluid_pp(i)%gamma
            isentrope_n(i) = f_isentrope_exponent(gammas(i))

            ! Each EOS supplies its own coefficients. Resolved once here, not per cell: a branch in the mixture loop costs
            ! registers in the Riemann kernels. An EOS whose coefficients depend on state must move to per-cell evaluation.
            select case (fluid_pp(i)%eos)
            case (eos_ideal_gas)
                pi_infs(i) = 0._wp
            case default
                pi_infs(i) = fluid_pp(i)%pi_inf
            end select
            isentrope_B(i) = f_isentrope_pressure(pi_infs(i), gammas(i))
            cvs(i) = fluid_pp(i)%cv
            qvs(i) = fluid_pp(i)%qv
            qvps(i) = fluid_pp(i)%qvp
            eoss(i) = fluid_pp(i)%eos
            ! Every fluid's single-source coefficients, mu_max among them: where a cubic Hugoniot fit turns over.
            ! mu(u_p) peaks where c0 = s2 u_p^2 + 2 s3 u_p^3, and past it no shock state exists, so the Newton below
            ! would wander. Solved once here, on the host.
            @:EOS_INIT_COEFFS(i)
            ! One reference state and Gruneisen closure for every family; the user-facing names keep their prefix.
            @:EOS_INIT_REFERENCE_STATE(i)
            if (f_is_state_dependent(i)) state_dependent = .true.
        end do
        ! Baked in at build time (every build, see case.py), so a case whose EOS family differs from the
        ! build's would silently run the wrong branch. The namelist still carries fluid_pp%eos, so check
        ! the two agree.
        @:PROHIBIT(state_dependent .neqv. any_state_dependent_eos, &
                   & "This case's equations of state do not match the ones this binary was built for. Rebuild with this case.")
        $:GPU_UPDATE(device='[gammas, isentrope_n, pi_infs, isentrope_B, cvs, qvs, qvps, eoss, eos_coeffs]')

    end subroutine s_initialize_eos_module

    !> Deallocate the fluid property arrays allocated in s_initialize_eos_module.
    impure subroutine s_finalize_eos_module()

        @:DEALLOCATE(gammas, isentrope_n, pi_infs, isentrope_B, cvs, qvs, qvps, eoss)

    end subroutine s_finalize_eos_module

    !> The reference curve of a state-dependent EOS at rho: p_ref, e_ref, their d/drho, and Gamma_G with its d/drho. A new family
    !! adds one case here and nothing else.
    subroutine s_reference_curve(rho, i, p_ref, e_ref, dp_drho, de_drho, G0, dG0)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)  :: rho
        integer, intent(in)   :: i
        real(wp), intent(out) :: p_ref, e_ref, dp_drho, de_drho, G0, dG0
        real(wp)              :: mu, d, V, ea, eb, up, us, dus, dup_dmu, x, ex, dp_dmu, de_dmu
        integer               :: iter

        mu = rho/eos_coeffs(i)%rho0 - 1._wp
        ! Past the fit's turnover there is no shock state to find; clamp rather than let the Newton below wander
        ! off and return a silently wrong pressure. mu_max is huge for the linear fit, so this is a no-op there.
        ! Bounded here so the step finishes and the host-side check in s_write_run_time_information can report it;
        ! past the turnover there is no shock state and the Newton below would wander.
        if (eoss(i) == eos_mie_gruneisen .and. mu > eos_coeffs(i)%mu_max) mu = eos_coeffs(i)%mu_max
        select case (eoss(i))
        case (eos_mie_gruneisen)
            ! Hugoniot reference u_s = c0 + s u_p + s2 u_p^2 + s3 u_p^3, with p_H = rho0 u_s u_p and the Hugoniot
            ! energy e_H = p_H mu/(2 rho0 (1 + mu)); linear on release. Pole at mu = 1/(s - 1) for the linear fit;
            ! the validator refuses initial states outside the EOS.
            if (mu < 0._wp) then
                p_ref = eos_coeffs(i)%rho0*eos_coeffs(i)%c0**2*mu
                dp_dmu = eos_coeffs(i)%rho0*eos_coeffs(i)%c0**2
            else if (eos_coeffs(i)%s2 == 0._wp .and. eos_coeffs(i)%s3 == 0._wp) then
                d = 1._wp - (eos_coeffs(i)%s - 1._wp)*mu
                p_ref = eos_coeffs(i)%rho0*eos_coeffs(i)%c0**2*mu*(1._wp + mu)/(d*d)
                dp_dmu = eos_coeffs(i)%rho0*eos_coeffs(i)%c0**2*((1._wp + 2._wp*mu)*d + 2._wp*(eos_coeffs(i)%s - 1._wp)*mu*(1._wp &
                                    & + mu))/(d*d*d)
            else
                ! u_p solves u_s(u_p) mu = u_p (1 + mu): Newton from the linear fit, then implicit differentiation
                up = eos_coeffs(i)%c0*mu/(1._wp - (eos_coeffs(i)%s - 1._wp)*mu)
                $:GPU_LOOP(parallelism='[seq]')
                do iter = 1, 8
                    us = eos_coeffs(i)%c0 + up*(eos_coeffs(i)%s + up*(eos_coeffs(i)%s2 + up*eos_coeffs(i)%s3))
                    dus = eos_coeffs(i)%s + up*(2._wp*eos_coeffs(i)%s2 + 3._wp*eos_coeffs(i)%s3*up)
                    up = up - (us*mu - up*(1._wp + mu))/(dus*mu - (1._wp + mu))
                end do
                us = eos_coeffs(i)%c0 + up*(eos_coeffs(i)%s + up*(eos_coeffs(i)%s2 + up*eos_coeffs(i)%s3))
                dus = eos_coeffs(i)%s + up*(2._wp*eos_coeffs(i)%s2 + 3._wp*eos_coeffs(i)%s3*up)
                dup_dmu = (up - us)/(dus*mu - (1._wp + mu))
                p_ref = eos_coeffs(i)%rho0*us*up
                dp_dmu = eos_coeffs(i)%rho0*(dus*up + us)*dup_dmu
            end if
            e_ref = p_ref*mu/(2._wp*eos_coeffs(i)%rho0*(1._wp + mu))
            de_dmu = (dp_dmu*mu*(1._wp + mu) + p_ref)/(2._wp*eos_coeffs(i)%rho0*(1._wp + mu)**2)
            dp_drho = dp_dmu/eos_coeffs(i)%rho0
            de_drho = de_dmu/eos_coeffs(i)%rho0
        case (eos_jwl)
            ! JWL: p_ref = A exp(-R1 V) + B exp(-R2 V), V = rho0/rho. The curve is itself an isentrope, so de_ref = -p_ref d(1/rho).
            V = eos_coeffs(i)%rho0/rho
            ea = eos_coeffs(i)%a*exp(-eos_coeffs(i)%r1*V)
            eb = eos_coeffs(i)%b*exp(-eos_coeffs(i)%r2*V)
            p_ref = ea + eb
            e_ref = (ea/eos_coeffs(i)%r1 + eb/eos_coeffs(i)%r2)/eos_coeffs(i)%rho0
            dp_drho = (eos_coeffs(i)%rho0/rho**2)*(eos_coeffs(i)%r1*ea + eos_coeffs(i)%r2*eb)
            de_drho = p_ref/rho**2
        case (eos_vinet)
            ! Vinet cold curve: p_c = 3 K0 (1 - x)/x^2 exp(eta (1 - x)), x = (rho0/rho)^(1/3), eta = 3 (K0' - 1)/2,
            ! an isentrope like JWL (its energy integrates in closed form).
            d = 1.5_wp*(eos_coeffs(i)%k0p - 1._wp)
            x = (eos_coeffs(i)%rho0/rho)**(1._wp/3._wp)
            ex = exp(d*(1._wp - x))
            p_ref = 3._wp*eos_coeffs(i)%k0*(1._wp - x)/x**2*ex
            e_ref = 9._wp*eos_coeffs(i)%k0/(eos_coeffs(i)%rho0*d**2)*(1._wp - (1._wp - d*(1._wp - x))*ex)
            dp_drho = 3._wp*eos_coeffs(i)%k0*ex*(-1._wp/x**2 - 2._wp*(1._wp - x)/x**3 - d*(1._wp - x)/x**2)*(-x/(3._wp*rho))
            de_drho = p_ref/rho**2
        end select
        G0 = eos_coeffs(i)%gruneisen0 + eos_coeffs(i)%gruneisen_a*mu
        dG0 = eos_coeffs(i)%gruneisen_a/eos_coeffs(i)%rho0

    end subroutine s_reference_curve

    !> Whether the EOS of fluid i is a family whose coefficients vary with density.
    function f_is_state_dependent(i) result(yes)

        $:GPU_ROUTINE(function_name='f_is_state_dependent', parallelism='[seq]', cray_inline=True)

        integer, intent(in) :: i
        logical             :: yes

        yes = ${EOS_IS_STATE_DEPENDENT('i')}$

    end function f_is_state_dependent

    !> True when fluid i's reference curve is itself an isentrope (de_ref = -p_ref d(1/rho), which holds for JWL and Vinet but not
    !! for the Mie-Gruneisen Hugoniot) and its Gruneisen coefficient is constant. Those two together make the isentrope through any
    !! state closed-form, so it never has to be integrated.
    function f_has_isentropic_reference(i) result(yes)

        $:GPU_ROUTINE(function_name='f_has_isentropic_reference', parallelism='[seq]', cray_inline=True)

        integer, intent(in) :: i
        logical             :: yes

        yes = (${EOS_HAS_ISENTROPIC_REFERENCE('i')}$) .and. eos_coeffs(i)%gruneisen_a == 0._wp

    end function f_has_isentropic_reference

    !> The largest compression a cubic Hugoniot fit can represent. mu(u_p) = u_p/(u_s - u_p) rises, peaks where c0 = s2 u_p^2 + 2 s3
    !! u_p^3, and falls after; only the rising branch is a physical shock. Returns a huge value for the linear fit, which never
    !! turns over. Host-side: called once per fluid at initialization.
    impure function f_hugoniot_compression_limit(c0, s, s2, s3) result(mu_max)

        real(wp), intent(in) :: c0, s, s2, s3
        real(wp)             :: mu_max, up, f, df, us
        integer              :: iter

        if (s2 == 0._wp .and. s3 == 0._wp) then
            mu_max = huge(1._wp)
            return
        end if

        ! Newton on c0 - s2 u^2 - 2 s3 u^3 = 0, from a guess that brackets the physical range
        up = c0
        do iter = 1, 100
            f = c0 - s2*up**2 - 2._wp*s3*up**3
            df = -2._wp*s2*up - 6._wp*s3*up**2
            if (abs(df) < verysmall) exit
            up = max(up - f/df, verysmall)
        end do
        us = c0 + up*(s + up*(s2 + up*s3))
        mu_max = up/max(us - up, verysmall)

    end function f_hugoniot_compression_limit

    !> Gamma, Pi, dPi/drho and dGamma/drho of fluid i at density rho, the coefficients of rho e = Gamma p + Pi(rho). Stiffened and
    !! ideal gas keep the constants resolved at init, bit for bit.
    subroutine s_eos_coefficients(rho, i, gamma, pi_inf, dpi, dgamma)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)  :: rho
        integer, intent(in)   :: i
        real(wp), intent(out) :: gamma, pi_inf, dpi, dgamma
        real(wp)              :: p_ref, e_ref, dp_drho, de_drho, G0, dG0

        if (.not. f_is_state_dependent(i)) then
            gamma = gammas(i)
            pi_inf = pi_infs(i)
            dpi = 0._wp
            dgamma = 0._wp
            return
        end if
        call s_reference_curve(rho, i, p_ref, e_ref, dp_drho, de_drho, G0, dG0)
        gamma = 1._wp/G0
        pi_inf = rho*e_ref - p_ref/G0
        dpi = e_ref + rho*de_drho - dp_drho/G0 + p_ref*dG0/G0**2
        dgamma = -dG0/G0**2

    end subroutine s_eos_coefficients

    !> Exponent of the stiffened-gas isentrope p + B = const rho**n. Precomputed per fluid as isentrope_n.
    function f_isentrope_exponent(gamma) result(n)

        $:GPU_ROUTINE(function_name='f_isentrope_exponent', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in) :: gamma
        real(wp)             :: n

        n = 1._wp/gamma + 1._wp

    end function f_isentrope_exponent

    !> Reference pressure of that isentrope. Precomputed per fluid as isentrope_B.
    function f_isentrope_pressure(pi_inf, gamma) result(B)

        $:GPU_ROUTINE(function_name='f_isentrope_pressure', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in) :: pi_inf, gamma
        real(wp)             :: B

        B = pi_inf/(1._wp + gamma)

    end function f_isentrope_pressure

    !> Stiffened-gas thermal law p + B = (n - 1)*cv*rho*T. Pass rho to get T, or T to get rho.
    function f_sg_thermal(pres, rho_or_T, n, B, cv) result(T_or_rho)

        $:GPU_ROUTINE(function_name='f_sg_thermal', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in) :: pres, rho_or_T, n, B, cv
        real(wp)             :: T_or_rho

        T_or_rho = (pres + B)/((n - 1._wp)*cv*rho_or_T)

    end function f_sg_thermal

    !> Thermal-equilibrium mixture temperature for stiffened gas, from primitives. Algebraically identical to the conservative form
    !! in m_phase_change's s_infinite_pt_relaxation_k, T = (rho*e + p - sum(alpha_rho_i*qv_i)) / sum(alpha_rho_i*cv_i*n_i), because
    !! rho*e = gamma_mix*p + pi_inf_mix + sum(alpha_rho_i*qv_i) in MFC's stored variables.
    function f_mixture_temperature(alpha_rho_K, pres, gamma_K, pi_inf_K) result(T)

        $:GPU_ROUTINE(function_name='f_mixture_temperature', parallelism='[seq]', cray_inline=True)

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(in) :: alpha_rho_K
        #:else
            real(wp), dimension(num_fluids), intent(in) :: alpha_rho_K
        #:endif
        real(wp), intent(in) :: pres, gamma_K, pi_inf_K
        real(wp)             :: T
        real(wp)             :: mCP  !< sum of alpha_rho_i*cp_i; cp_i = n_i*cv_i
        integer              :: i

        mCP = 0._wp
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            mCP = mCP + alpha_rho_K(i)*cvs(i)*isentrope_n(i)
        end do

        T = ((gamma_K + 1._wp)*pres + pi_inf_K)/max(mCP, sgm_eps)

    end function f_mixture_temperature

    !> Coefficients of phase i at its own density alpha_rho/alpha: the per-cell dispatch when some fluid's EOS is state dependent,
    !! the constants resolved at init otherwise (bit for bit).
    subroutine s_phase_coefficients(alpha_rho, alpha, i, rho, gamma, pi_inf, dpi, dgamma)

        $:GPU_ROUTINE(function_name='s_phase_coefficients', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in)  :: alpha_rho, alpha
        integer, intent(in)   :: i
        real(wp), intent(out) :: rho, gamma, pi_inf, dpi, dgamma

        rho = max(alpha_rho, sgm_eps)/max(alpha, sgm_eps)
        if (any_state_dependent_eos) then
            call s_eos_coefficients(rho, i, gamma, pi_inf, dpi, dgamma)
        else
            gamma = gammas(i)
            pi_inf = pi_infs(i)
            dpi = 0._wp
            dgamma = 0._wp
        end if

    end subroutine s_phase_coefficients

    !> c^2 = [((Gamma + 1) p + Pi)/rho - dPi/drho - p dGamma/drho]/Gamma, the frozen speed of one phase.
    function f_c2_from_coefficients(rho, pres, gamma, pi_inf, dpi, dgamma) result(c2)

        $:GPU_ROUTINE(function_name='f_c2_from_coefficients', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in) :: rho, pres, gamma, pi_inf, dpi, dgamma
        real(wp)             :: c2

        c2 = (((gamma + 1._wp)*pres + pi_inf)/rho - dpi - pres*dgamma)/gamma

    end function f_c2_from_coefficients

    !> Frozen sound speed squared of one phase at (rho, p) from its own coefficients. These helpers are subroutines, not functions:
    !! a device function that calls a device subroutine is a pattern no other backend-tested code in MFC uses.
    subroutine s_phase_c2(rho, pres, i, c2)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)  :: rho, pres
        integer, intent(in)   :: i
        real(wp), intent(out) :: c2
        real(wp)              :: gamma, pi_inf, dpi, dgamma

        call s_eos_coefficients(rho, i, gamma, pi_inf, dpi, dgamma)
        c2 = f_c2_from_coefficients(rho, pres, gamma, pi_inf, dpi, dgamma)

    end subroutine s_phase_c2

    !> Slope of the ODE `kind` for fluid i: dp/drho = c^2 along an isentrope (x = rho, y = p), or the reference temperature dT/dV =
    !! (de_ref/dV + p_ref)/c_v - Gamma_G T/V (x = V, y = T), the Maxwell relation applied to e = e_ref + c_v (T - T_ref).
    subroutine s_ode_slope(kind, i, x, y, dydx)

        $:GPU_ROUTINE(parallelism='[seq]')

        integer, intent(in)   :: kind, i
        real(wp), intent(in)  :: x, y
        real(wp), intent(out) :: dydx
        real(wp)              :: p_ref, e_ref, dp_drho, de_drho, G0, dG0

        if (kind == ode_isentrope) then
            call s_phase_c2(x, y, i, dydx)
        else
            call s_reference_curve(1._wp/x, i, p_ref, e_ref, dp_drho, de_drho, G0, dG0)
            dydx = (p_ref - de_drho/x**2)/cvs(i) - G0*y/x
        end if

    end subroutine s_ode_slope

    !> Fixed-step classical RK4 for the ODE `kind` from (x0, y0) to x1.
    subroutine s_rk4(kind, i, x0, y0, x1, y)

        $:GPU_ROUTINE(parallelism='[seq]')

        integer, intent(in)   :: kind, i
        real(wp), intent(in)  :: x0, y0, x1
        real(wp), intent(out) :: y
        real(wp)              :: x, h, k1, k2, k3, k4
        integer               :: step

        x = x0
        y = y0
        h = (x1 - x0)/eos_rk4_steps
        $:GPU_LOOP(parallelism='[seq]')
        do step = 1, eos_rk4_steps
            call s_ode_slope(kind, i, x, y, k1)
            call s_ode_slope(kind, i, x + 0.5_wp*h, y + 0.5_wp*h*k1, k2)
            call s_ode_slope(kind, i, x + 0.5_wp*h, y + 0.5_wp*h*k2, k3)
            call s_ode_slope(kind, i, x + h, y + h*k3, k4)
            y = y + h*(k1 + 2._wp*(k2 + k3) + k4)/6._wp
            x = x + h
        end do

    end subroutine s_rk4

    !> Pressure of phase i after the isentropic density change rho -> xi rho: closed form for the constant-coefficient families,
    !! integrated for a state-dependent EOS (the star states it serves are close to rho).
    subroutine s_phase_pressure_on_isentrope(pres, rho, xi, i, p_isen)

        $:GPU_ROUTINE(function_name='s_phase_pressure_on_isentrope', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in)  :: pres, rho, xi
        integer, intent(in)   :: i
        real(wp), intent(out) :: p_isen
        real(wp)              :: p_ref_from, p_ref_to, e_ref, dp_drho, de_drho, G0, dG0

        if (.not. f_is_state_dependent(i)) then
            p_isen = (pres + isentrope_B(i))*xi**isentrope_n(i) - isentrope_B(i)
        else if (f_has_isentropic_reference(i)) then
            ! Exact: the offset from an isentropic reference obeys dDelta/Delta = Gamma drho/rho, so
            ! p - p_ref scales as (rho'/rho)**(1 + Gamma). Integrating it instead costs a decimal per
            ! doubling of the expansion and turns the pressure negative past roughly twentyfold.
            call s_reference_curve(rho, i, p_ref_from, e_ref, dp_drho, de_drho, G0, dG0)
            call s_reference_curve(xi*rho, i, p_ref_to, e_ref, dp_drho, de_drho, G0, dG0)
            p_isen = p_ref_to + (pres - p_ref_from)*xi**(1._wp + eos_coeffs(i)%gruneisen0)
        else
            call s_rk4(ode_isentrope, i, rho, pres, xi*rho, p_isen)
        end if

    end subroutine s_phase_pressure_on_isentrope

    !> Temperature of phase i at (rho, p): the stiffened-gas relation, or T_ref(rho) + (e - e_ref)/c_v.
    subroutine s_phase_temperature(rho, pres, i, T)

        $:GPU_ROUTINE(function_name='s_phase_temperature', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in)  :: rho, pres
        integer, intent(in)   :: i
        real(wp), intent(out) :: T
        real(wp)              :: p_ref, e_ref, dp_drho, de_drho, G0, dG0, T0, T_ref

        if (f_is_state_dependent(i)) then
            call s_reference_curve(rho, i, p_ref, e_ref, dp_drho, de_drho, G0, dG0)
            T0 = eos_coeffs(i)%t0
            call s_rk4(ode_reference_temperature, i, 1._wp/eos_coeffs(i)%rho0, T0, 1._wp/rho, T_ref)
            T = T_ref + (pres - p_ref)/(rho*G0*cvs(i))
        else
            T = (pres + isentrope_B(i))/((isentrope_n(i) - 1._wp)*cvs(i)*rho)
        end if

    end subroutine s_phase_temperature

    !> Density of phase i on the isentrope through (rho_from, p_from) at p_to, and c^2 there: Newton on the pressure integrator,
    !! whose slope is c^2. The relaxation's own Newton wraps this, so a few steps suffice.
    subroutine s_phase_density_on_isentrope(i, rho_from, p_from, p_to, rho_to, c2_to)

        $:GPU_ROUTINE(function_name='s_phase_density_on_isentrope', parallelism='[seq]')

        integer, intent(in)   :: i
        real(wp), intent(in)  :: rho_from, p_from, p_to
        real(wp), intent(out) :: rho_to, c2_to
        real(wp)              :: p_at, c2_at
        integer               :: iter

        rho_to = rho_from
        $:GPU_LOOP(parallelism='[seq]')
        do iter = 1, 4
            call s_phase_pressure_on_isentrope(p_from, rho_from, rho_to/rho_from, i, p_at)
            call s_phase_c2(rho_to, p_at, i, c2_at)
            rho_to = rho_to - (p_at - p_to)/c2_at
        end do
        call s_phase_c2(rho_to, p_to, i, c2_to)

    end subroutine s_phase_density_on_isentrope

    !> Internal energy per unit volume of phase i at pressure pres: alpha (Gamma p + Pi) + alpha_rho qv, with the coefficients at
    !! the phase's own density.
    subroutine s_phase_internal_energy(pres, alpha, alpha_rho, i, e_phase)

        $:GPU_ROUTINE(function_name='s_phase_internal_energy', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in)  :: pres, alpha, alpha_rho
        integer, intent(in)   :: i
        real(wp), intent(out) :: e_phase
        real(wp)              :: rho, gamma, pi_inf, dpi, dgamma

        call s_phase_coefficients(alpha_rho, alpha, i, rho, gamma, pi_inf, dpi, dgamma)
        e_phase = alpha*(gamma*pres + pi_inf) + alpha_rho*qvs(i)

    end subroutine s_phase_internal_energy

    !> Bulk modulus rho c^2 of phase i at pressure pres: f_bulk_modulus for a constant-coefficient fluid, bit for bit, minus the
    !! reference-curve terms rho (dPi/drho + p dGamma/drho)/Gamma otherwise.
    subroutine s_phase_bulk_modulus(pres, alpha, alpha_rho, i, blkmod)

        $:GPU_ROUTINE(function_name='s_phase_bulk_modulus', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in)  :: alpha_rho, alpha, pres
        integer, intent(in)   :: i
        real(wp), intent(out) :: blkmod
        real(wp)              :: rho, gamma, pi_inf, dpi, dgamma

        call s_phase_coefficients(alpha_rho, alpha, i, rho, gamma, pi_inf, dpi, dgamma)
        blkmod = f_bulk_modulus(pres, gamma, pi_inf) - rho*(dpi + pres*dgamma)/gamma

    end subroutine s_phase_bulk_modulus

    !> Pressure of a stiffened gas from its internal energy density - the inverse of s_compute_energy. Callers subtract the kinetic,
    !! magnetic and elastic energy first; none of those are equation-of-state terms.
    function f_pressure(e_int, gamma, pi_inf, qv) result(pres)

        $:GPU_ROUTINE(function_name='f_pressure', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in) :: e_int, gamma, pi_inf, qv
        real(wp)             :: pres

        pres = (e_int - pi_inf - qv)/gamma

    end function f_pressure

    !> Isentropic bulk modulus. Takes coefficients rather than a fluid index, so a mixture - whose effective gamma and pi_inf come
    !! from s_compute_mixture_coefficients - is the same call as a single fluid. Elastic callers add their own shear term.
    function f_bulk_modulus(pres, gamma, pi_inf) result(blkmod)

        $:GPU_ROUTINE(function_name='f_bulk_modulus', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in) :: pres, gamma, pi_inf
        real(wp)             :: blkmod

        blkmod = ((gamma + 1._wp)*pres + pi_inf)/gamma

    end function f_bulk_modulus

    !> Relativistic specific enthalpy, h = 1 + (Gamma + 1)p/rho. Ideal gas only: the stiffness does not appear, so a fluid with a
    !! nonzero pi_inf is not represented here (the validator refuses that combination).
    function f_relativistic_enthalpy(pres, rho, gamma) result(H)

        $:GPU_ROUTINE(function_name='f_relativistic_enthalpy', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in) :: pres, rho, gamma
        real(wp)             :: H

        H = 1._wp + (gamma + 1._wp)*pres/rho

    end function f_relativistic_enthalpy

    !> Mixture coefficients of one state. Under bubbles_euler with num_fluids == 1 the sole advection slot aliases the void fraction
    !! (eqn_idx%alf == eqn_idx%adv%end), so alpha is not a composition there and the coefficients are the liquid's. Clipping stays
    !! with callers; it differs between solvers and cannot coincide with that case, as mpp_lim requires num_fluids > 1.
    subroutine s_compute_mixture_coefficients(alpha_rho_K, alpha_K, rho_K, gamma_K, pi_inf_K, qv_K)

        $:GPU_ROUTINE(function_name='s_compute_mixture_coefficients', parallelism='[seq]', cray_inline=True)

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(in) :: alpha_rho_K, alpha_K
        #:else
            real(wp), dimension(num_fluids), intent(in) :: alpha_rho_K, alpha_K
        #:endif
        real(wp), intent(out) :: rho_K, gamma_K, pi_inf_K, qv_K
        real(wp)              :: gamma_i, pi_inf_i, dpi_i, dgamma_i
        real(wp)              :: rho_i, alpha_i, alpha_rho_i
        integer               :: i  !< Loop iterator over fluids

        ! The bubbly closure is written for one carrier liquid, which keeps its own coefficients
        ! undiluted: Gamma_l*p_l = (E - rho|u|^2/2)/(1 - alf) - Pi_inf_l, the void entering only through
        ! the (1 - alf) that s_compute_pressure applies. There is nothing to sum - the last advection
        ! slot is the void, not a material - and the checker holds num_fluids <= 2 here.
        if (bubbles_euler) then
            rho_K = alpha_rho_K(1)
            gamma_K = gammas(1)
            pi_inf_K = pi_infs(1)
            ! Energy per unit volume, as below: alpha_rho_K(1) is the liquid partial density
            qv_K = alpha_rho_K(1)*qvs(1)
        else
            rho_K = 0._wp
            gamma_K = 0._wp
            pi_inf_K = 0._wp
            qv_K = 0._wp

            ! Stiffened-gas fast path: the phase coefficients are the constants gammas/pi_infs, and calling
            ! s_phase_coefficients per fluid per cell drags the state-dependent EOS chain (reference-curve Newton loop)
            ! into every conversion and Riemann kernel even when no fluid uses it.
            ! Same arithmetic as the general branch with gamma_i = gammas(i), pi_inf_i = pi_infs(i).
            if (.not. any_state_dependent_eos) then
                $:GPU_LOOP(parallelism='[seq]')
                do i = 1, num_fluids
                    rho_K = rho_K + alpha_rho_K(i)
                    gamma_K = gamma_K + alpha_K(i)*gammas(i)
                    pi_inf_K = pi_inf_K + alpha_K(i)*pi_infs(i)
                    qv_K = qv_K + alpha_rho_K(i)*qvs(i)
                end do
            else
                $:GPU_LOOP(parallelism='[seq]')
                do i = 1, num_fluids
                    rho_K = rho_K + alpha_rho_K(i)
                    alpha_rho_i = alpha_rho_K(i)
                    alpha_i = alpha_K(i)
                    call s_phase_coefficients(alpha_rho_i, alpha_i, i, rho_i, gamma_i, pi_inf_i, dpi_i, dgamma_i)
                    gamma_K = gamma_K + alpha_K(i)*gamma_i
                    pi_inf_K = pi_inf_K + alpha_K(i)*pi_inf_i
                    qv_K = qv_K + alpha_rho_K(i)*qvs(i)
                end do
            end if
        end if

    end subroutine s_compute_mixture_coefficients

    !> Time derivative of the mixture coefficients, mirroring s_compute_mixture_coefficients.
    subroutine s_compute_mixture_coefficients_dt(dalpha_rho_dt, dadv_dt, alpha_rho, adv, drho_dt, dgamma_dt, dpi_inf_dt, dqv_dt)

        $:GPU_ROUTINE(function_name='s_compute_mixture_coefficients_dt', parallelism='[seq]', cray_inline=True)

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(in) :: dalpha_rho_dt, dadv_dt, alpha_rho, adv
        #:else
            real(wp), dimension(num_fluids), intent(in) :: dalpha_rho_dt, dadv_dt, alpha_rho, adv
        #:endif
        real(wp), intent(out) :: drho_dt, dgamma_dt, dpi_inf_dt, dqv_dt
        real(wp)              :: rho_i, gamma_i, pi_inf_i, dpi_i, dgamma_i, alpha_i, alpha_rho_i
        integer               :: i  !< Loop iterator over fluids

        dgamma_dt = 0._wp
        dpi_inf_dt = 0._wp
        dqv_dt = 0._wp

        if (num_fluids == 1 .and. bubbles_euler) then
            ! Fluid 1's coefficients are constants here, so only rho varies.
            drho_dt = dalpha_rho_dt(1)
        else
            drho_dt = 0._wp

            $:GPU_LOOP(parallelism='[seq]')
            do i = 1, num_fluids
                drho_dt = drho_dt + dalpha_rho_dt(i)
                alpha_rho_i = alpha_rho(i)
                alpha_i = adv(i)
                call s_phase_coefficients(alpha_rho_i, alpha_i, i, rho_i, gamma_i, pi_inf_i, dpi_i, dgamma_i)
                ! d(alpha X(rho_i))/dt with rho_i = alpha_rho/alpha; the alpha in dX/dt cancels
                dgamma_dt = dgamma_dt + dadv_dt(i)*gamma_i + dgamma_i*(dalpha_rho_dt(i) - rho_i*dadv_dt(i))
                dpi_inf_dt = dpi_inf_dt + dadv_dt(i)*pi_inf_i + dpi_i*(dalpha_rho_dt(i) - rho_i*dadv_dt(i))
                dqv_dt = dqv_dt + dalpha_rho_dt(i)*qvs(i)
            end do
        end if

    end subroutine s_compute_mixture_coefficients_dt

    !> Speed of sound of a thermodynamic state. Enthalpy is not an argument: for a real state H, |u|^2 and qv all cancel out of c^2
    !! = ((Gamma + 1)p + Pi)/(Gamma rho). Averaged states, whose enthalpy is a free input, use the _avg variant.
    subroutine s_compute_speed_of_sound(pres, rho, gamma, pi_inf, adv, c, alpha_rho)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in) :: pres, rho, gamma, pi_inf
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(in) :: adv
        #:else
            real(wp), dimension(num_fluids), intent(in) :: adv
        #:endif
        real(wp), intent(out) :: c
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(in), optional :: alpha_rho
        #:else
            real(wp), dimension(num_fluids), intent(in), optional :: alpha_rho
        #:endif
        real(wp) :: alf  !< Subgrid void fraction; dilute by construction
        real(wp) :: blkmod_q, alpha_q, alpha_rho_q, gamma_q, pi_inf_q
        integer  :: q

        if (chemistry) then  ! Reacting mixture sound speed
            c = sqrt((1.0_wp + 1.0_wp/gamma)*pres/rho)
        else if (relativity) then  ! Relativistic sound speed, whose enthalpy is 1 + (Gamma + 1)p/rho
            c = sqrt((1._wp + 1._wp/gamma)*pres/rho/f_relativistic_enthalpy(pres, rho, gamma))
        else
            ! Every case below is a bulk modulus over a density. The equation of state enters
            ! only through f_bulk_modulus; the cases differ in how the phases are mixed.
            if (any_state_dependent_eos .and. present(alpha_rho)) then  ! frozen mixing: each phase's modulus at its own density
                c = 0._wp
                $:GPU_LOOP(parallelism='[seq]')
                do q = 1, num_fluids
                    alpha_q = adv(q)
                    alpha_rho_q = alpha_rho(q)
                    call s_phase_bulk_modulus(pres, alpha_q, alpha_rho_q, q, blkmod_q)
                    if (alt_soundspeed) then
                        c = c + adv(q)/blkmod_q
                    else
                        c = c + adv(q)*blkmod_q
                    end if
                end do
                if (alt_soundspeed) then
                    c = 1._wp/(rho*c)
                else
                    c = c/rho
                end if
            else if (alt_soundspeed) then  ! Wood's law: volume-weighted harmonic mean
                c = 0._wp
                $:GPU_LOOP(parallelism='[seq]')
                do q = 1, num_fluids
                    gamma_q = gammas(q)
                    pi_inf_q = pi_infs(q)
                    c = c + adv(q)/f_bulk_modulus(pres, gamma_q, pi_inf_q)
                end do
                c = 1._wp/(rho*c)
            else if (model_eqns == model_eqns_6eq) then  ! volume-weighted arithmetic mean
                c = 0._wp
                $:GPU_LOOP(parallelism='[seq]')
                do q = 1, num_fluids
                    gamma_q = gammas(q)
                    pi_inf_q = pi_infs(q)
                    c = c + adv(q)*f_bulk_modulus(pres, gamma_q, pi_inf_q)
                end do
                c = c/rho
            else  ! the mixture coefficients already carry the mixing
                c = f_bulk_modulus(pres, gamma, pi_inf)/rho

                ! Subgrid bubbles: c = c_l/(1 - alf), the carrier-liquid speed with an O(alf) void
                ! correction. alf is dilute by construction; near one means a wrong index or an
                ! out-of-regime case, which the toolchain warns about at case load (#1793).
                if (model_eqns == model_eqns_5eq .and. bubbles_euler .and. .not. (mpp_lim .and. num_fluids > 1)) then
                    alf = adv(num_fluids)
                    c = c/(1._wp - alf)
                end if
            end if

            if (mixture_err .and. c < 0._wp) then
                c = 100._wp*sgm_eps
            else
                c = sqrt(c)
            end if
        end if

    end subroutine s_compute_speed_of_sound

    !> Speed of sound of an interface-averaged state. An average of two states is not a state - its enthalpy is not the one its
    !! pressure and density imply - so the caller supplies H, |u|^2 and qv. Only the enthalpy-reading branches differ from
    !! s_compute_speed_of_sound; keep the condition below in step with the branch list there.
    subroutine s_compute_speed_of_sound_avg(pres, rho, gamma, pi_inf, qv, vel_sum, H, c_c, adv, c, alpha_rho)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in) :: pres, rho, gamma, pi_inf, qv, vel_sum, H, c_c
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(in) :: adv
        #:else
            real(wp), dimension(num_fluids), intent(in) :: adv
        #:endif
        real(wp), intent(out) :: c
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(in), optional :: alpha_rho
        #:else
            real(wp), dimension(num_fluids), intent(in), optional :: alpha_rho
        #:endif

        if (chemistry) then  ! Reacting mixture sound speed
            if (avg_state == avg_state_roe .and. abs(c_c) > verysmall) then
                c = sqrt(c_c - (gamma - 1.0_wp)*(vel_sum - H))
            else
                call s_compute_speed_of_sound(pres, rho, gamma, pi_inf, adv, c, alpha_rho)
            end if
        else if (relativity) then  ! Relativistic sound speed
            c = sqrt((1._wp + 1._wp/gamma)*pres/rho/H)
        else if (alt_soundspeed .or. model_eqns == model_eqns_6eq .or. (model_eqns == model_eqns_5eq .and. bubbles_euler) &
                 & .or. any_state_dependent_eos) then
            call s_compute_speed_of_sound(pres, rho, gamma, pi_inf, adv, c, alpha_rho)
        else  ! Stiffened-gas mixture, the one branch where the averaged enthalpy survives
            c = (H - 5.e-1*vel_sum - qv/rho)/gamma

            if (mixture_err .and. c < 0._wp) then
                c = 100._wp*sgm_eps
            else
                c = sqrt(c)
            end if
        end if

    end subroutine s_compute_speed_of_sound_avg

end module m_eos
