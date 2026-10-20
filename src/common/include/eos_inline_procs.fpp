!>
!! @file eos_inline_procs.fpp
!! @brief Hot equation-of-state routines, one source for every module that calls them from a GPU kernel.
!!
!! nvfortran 26.5 inlines a routine into an OpenACC kernel only when it is defined in the same file; a routine from
!! another module stays a real device CALL (-Minline=lib inlines host-side only). These routines are called per face
!! and per cell, so they are kept in the Fypp macro EOS_INLINE_PROCS below: m_eos expands it for the public definitions
!! and each kernel module expands it again inside its own contains section for a private copy. The kernel module must rename the public ones on use
!! (@:EOS_INLINE_USE) and mark its copies private ($:EOS_INLINE_PRIVATE()), both from eos_inline.fpp.
!! Edit the routines here only.

#:def EOS_INLINE_PROCS()
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

    end subroutine s_compute_mixture_coefficients

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

    !> Total energy per unit volume, thermodynamic terms only. Callers add magnetic and elastic energy, which are not
    !! equation-of-state terms. The chemistry and relativistic branches use a different relation and stay open-coded.
    subroutine s_compute_energy(pres, alpha_rho_K, alpha_K, vel_sum, E)

        $:GPU_ROUTINE(function_name='s_compute_energy', parallelism='[seq]', cray_inline=True)

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(in) :: alpha_rho_K, alpha_K
        #:else
            real(wp), dimension(num_fluids), intent(in) :: alpha_rho_K, alpha_K
        #:endif
        real(wp), intent(in)  :: pres, vel_sum
        real(wp), intent(out) :: E
        real(wp)              :: rho, gamma, pi_inf, qv

        call s_compute_mixture_coefficients(alpha_rho_K, alpha_K, rho, gamma, pi_inf, qv)

        ! E = dyn_p + (1 - alf)(gamma p + pi_inf) + qv. Only the liquid's internal energy is diluted; qv is already an energy
        ! density.
        E = gamma*pres + pi_inf
        if (bubbles_euler) E = E*(1._wp - alpha_K(num_fluids))
        E = E + qv + 5.e-1_wp*rho*vel_sum

    end subroutine s_compute_energy

    !> Whether the EOS of fluid i is a family whose coefficients vary with density.
    function f_is_state_dependent(i) result(yes)

        $:GPU_ROUTINE(function_name='f_is_state_dependent', parallelism='[seq]', cray_inline=True)

        integer, intent(in) :: i
        logical             :: yes

        #:if defined('EOS_INLINE_NAMES')
            #:if MFC_CASE_OPTIMIZATION and not any_state_dependent_eos
                ! Module-local copy in a case-optimized build without a state-dependent EOS: s_initialize_eos_module
                ! aborts unless every fluid is state-independent, so this is a constant and the state-dependent
                ! branches of the phase routines drop out of the kernels.
                yes = .false.
            #:else
                ! Module-local copy: the family table is generated for m_eos only, so ask the public routine.
                yes = f_is_state_dependent_ext(i)
            #:endif
        #:else
            yes = ${EOS_IS_STATE_DEPENDENT('i')}$
        #:endif

    end function f_is_state_dependent

    !> Pressure of a stiffened gas from its internal energy density - the inverse of s_compute_energy. Callers subtract the kinetic,
    !! magnetic and elastic energy first; none of those are equation-of-state terms.
    function f_pressure(e_int, gamma, pi_inf, qv) result(pres)

        $:GPU_ROUTINE(function_name='f_pressure', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in) :: e_int, gamma, pi_inf, qv
        real(wp)             :: pres

        pres = (e_int - pi_inf - qv)/gamma

    end function f_pressure

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
#:enddef
