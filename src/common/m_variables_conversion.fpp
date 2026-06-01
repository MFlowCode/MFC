!>
!! @file
!! @brief Contains module m_variables_conversion

#:include 'macros.fpp'
#:include 'case.fpp'

!> @brief Conservative-to-primitive variable conversion, mixture property evaluation, and pressure computation
module m_variables_conversion

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_mpi_common, only: s_mpi_allreduce_integer_sum
    use m_helper_basic
    use m_helper
    use m_thermochem, only: num_species, get_temperature, get_pressure, gas_constant, get_mixture_molecular_weight, &
        & get_mixture_energy_mass

    implicit none

    private
    public :: s_initialize_variables_conversion_module, &
              s_initialize_pb, &
              s_initialize_mv, &
              s_convert_to_mixture_variables, &
              s_convert_mixture_to_mixture_variables, &
              s_convert_species_to_mixture_variables, &
              s_convert_species_to_mixture_variables_acc, &
              s_convert_conservative_to_primitive_variables, &
              s_convert_primitive_to_conservative_variables, &
              s_convert_primitive_to_flux_variables, &
              s_compute_pressure, &
              s_jwl_pcold, &
              s_jwl_sound_speed_squared, &
              s_jwl_energy_pr, &
              s_jwl_reactive_pressure_er, &
              s_jwl_reactive_energy_pr, &
              f_lee_tarver_rate, &
              jwl_idx, &
              s_compute_species_fraction, &
#ifndef MFC_PRE_PROCESS
    s_compute_speed_of_sound, &
              s_compute_fast_magnetosonic_speed, &
#endif
    s_finalize_variables_conversion_module

    ! In simulation, gammas, pi_infs, and qvs are already declared in m_global_variables
#ifndef MFC_SIMULATION
    integer, allocatable, public, dimension(:)  :: eos_idxs
    real(wp), allocatable, public, dimension(:) :: gammas, gs_min, pi_infs, ps_inf, cvs, qvs, qvps
    real(wp), allocatable, public, dimension(:) :: jwl_As, jwl_Bs, jwl_R1s, jwl_R2s, jwl_omegas, jwl_rho0s, jwl_E0s
    real(wp), allocatable, public, dimension(:) :: jwl_air_e0s, jwl_air_rho0s, jwl_air_gammas
    $:GPU_DECLARE(create='[eos_idxs, gammas, gs_min, pi_infs, ps_inf, cvs, qvs, qvps]')
    $:GPU_DECLARE(create='[jwl_As, jwl_Bs, jwl_R1s, jwl_R2s, jwl_omegas, jwl_rho0s, jwl_E0s]')
    $:GPU_DECLARE(create='[jwl_air_e0s, jwl_air_rho0s, jwl_air_gammas]')
#endif

    real(wp), allocatable, dimension(:)   :: Gs_vc
    integer, allocatable, dimension(:)    :: bubrs_vc
    real(wp), allocatable, dimension(:,:) :: Res_vc
    $:GPU_DECLARE(create='[bubrs_vc, Gs_vc, Res_vc]')

    integer :: is1b, is2b, is3b, is1e, is2e, is3e
    integer :: jwl_idx
    integer :: jwl_pressure_floor_count
    $:GPU_DECLARE(create='[is1b, is2b, is3b, is1e, is2e, is3e]')
    $:GPU_DECLARE(create='[jwl_idx, jwl_pressure_floor_count]')

    real(wp), allocatable, dimension(:,:,:), public :: rho_sf     !< Scalar density function
    real(wp), allocatable, dimension(:,:,:), public :: gamma_sf   !< Scalar sp. heat ratio function
    real(wp), allocatable, dimension(:,:,:), public :: pi_inf_sf  !< Scalar liquid stiffness function
    real(wp), allocatable, dimension(:,:,:), public :: qv_sf      !< Scalar liquid energy reference function

contains

    !> Dispatch to the s_convert_mixture_to_mixture_variables and s_convert_species_to_mixture_variables subroutines. Replaces a
    !! procedure pointer.
    subroutine s_convert_to_mixture_variables(q_vf, i, j, k, rho, gamma, pi_inf, qv, Re_K, G_K, G)

        type(scalar_field), dimension(sys_size), intent(in)   :: q_vf
        integer, intent(in)                                   :: i, j, k
        real(wp), intent(out), target                         :: rho, gamma, pi_inf, qv
        real(wp), optional, dimension(2), intent(out)         :: Re_K
        real(wp), optional, intent(out)                       :: G_K
        real(wp), optional, dimension(num_fluids), intent(in) :: G

        if (model_eqns == 1) then  ! Gamma/pi_inf model
            call s_convert_mixture_to_mixture_variables(q_vf, i, j, k, rho, gamma, pi_inf, qv)
        else  ! Volume fraction model
            call s_convert_species_to_mixture_variables(q_vf, i, j, k, rho, gamma, pi_inf, qv, Re_K, G_K, G)
        end if

    end subroutine s_convert_to_mixture_variables

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

        e = (v/omega)*(pres - A*(1._wp - omega/(R1*v))*exp1 - B*(1._wp - omega/(R2*v))*exp2)
        c2 = A*exp1*(omega/(R1*v**2._wp) - R1*(1._wp - omega/(R1*v))) + B*exp2*(omega/(R2*v**2._wp) - R2*(1._wp - omega/(R2*v)))
        c2 = c2*(-rho0/(rho_safe*rho_safe)) + omega*e/rho0 + omega*pres/(rho_safe*rho0)
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
        Kj = rho0/max(omega0, sgm_eps)

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
        Kj = rho0/max(omega0, sgm_eps)

        e = (pres*(a_j*Kj + a_a/max(air_gamma, sgm_eps)) - a_j*Kj*pref1)/rho_safe
        e = max(e, 0._wp)

    end subroutine s_jwl_energy_pr

    !> Reactive (progressive-burn) JWL/ideal-gas mixture pressure from specific internal energy and density. The condensed phase is
    !! a reaction-progress (lambda) blend of the unreacted-explosive JWL EOS (lambda=0) and the products JWL EOS (lambda=1), so the
    !! products' pressure and density evolve as the burn proceeds instead of being a fixed constant state. The unreacted-explosive
    !! constants are taken from the module globals jwl_unr_*. Reduces to s_jwl_pressure_er when lambda=1.
    subroutine s_jwl_reactive_pressure_er(rho, e, Y, alpha_j, lambda, A, B, R1, R2, omega0, rho0, air_gamma, pres)

        $:GPU_ROUTINE(function_name='s_jwl_reactive_pressure_er',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, e, Y, alpha_j, lambda, A, B, R1, R2, omega0, rho0, air_gamma
        real(wp), intent(out) :: pres
        real(wp)              :: rho_safe, Y_safe, a_j, a_a, rho1, rhoe, lam, pref_p, pref_u, pref1, K_p, K_u, Kj

        rho_safe = max(rho, sgm_eps)
        Y_safe = min(max(Y, 0._wp), 1._wp)
        a_j = min(max(alpha_j, 0._wp), 1._wp)
        a_a = 1._wp - a_j
        rhoe = rho_safe*e
        lam = min(max(lambda, 0._wp), 1._wp)

        if (a_j <= sgm_eps) then
            pres = air_gamma*rhoe
            return
        end if

        if (a_a <= sgm_eps) then
            rho1 = rho_safe
        else
            rho1 = max(Y_safe*rho_safe/a_j, sgm_eps)
        end if

        pref_p = A*(1._wp - omega0*rho1/(R1*rho0))*exp(-R1*rho0/rho1) + B*(1._wp - omega0*rho1/(R2*rho0))*exp(-R2*rho0/rho1)
        K_p = rho0/max(omega0, sgm_eps)
        pref_u = jwl_unr_A*(1._wp - jwl_unr_omega*rho1/(jwl_unr_R1*jwl_unr_rho0))*exp(-jwl_unr_R1*jwl_unr_rho0/rho1) &
                            & + jwl_unr_B*(1._wp - jwl_unr_omega*rho1/(jwl_unr_R2*jwl_unr_rho0))*exp(-jwl_unr_R2*jwl_unr_rho0/rho1)
        K_u = jwl_unr_rho0/max(jwl_unr_omega, sgm_eps)

        pref1 = (1._wp - lam)*pref_u + lam*pref_p
        Kj = (1._wp - lam)*K_u + lam*K_p

        pres = (rhoe + a_j*Kj*pref1)/max(a_j*Kj + a_a/max(air_gamma, sgm_eps), sgm_eps)

    end subroutine s_jwl_reactive_pressure_er

    !> Specific internal energy from pressure and density for the reactive JWL mixture: exact inverse of s_jwl_reactive_pressure_er,
    !! keeping primitive<->conservative consistent across the burn.
    subroutine s_jwl_reactive_energy_pr(rho, pres, Y, alpha_j, lambda, A, B, R1, R2, omega0, rho0, air_gamma, e)

        $:GPU_ROUTINE(function_name='s_jwl_reactive_energy_pr',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, pres, Y, alpha_j, lambda, A, B, R1, R2, omega0, rho0, air_gamma
        real(wp), intent(out) :: e
        real(wp)              :: rho_safe, Y_safe, a_j, a_a, rho1, lam, pref_p, pref_u, pref1, K_p, K_u, Kj

        rho_safe = max(rho, sgm_eps)
        Y_safe = min(max(Y, 0._wp), 1._wp)
        a_j = min(max(alpha_j, 0._wp), 1._wp)
        a_a = 1._wp - a_j
        lam = min(max(lambda, 0._wp), 1._wp)

        if (a_j <= sgm_eps) then
            e = max(pres/max(air_gamma*rho_safe, sgm_eps), 0._wp)
            return
        end if

        if (a_a <= sgm_eps) then
            rho1 = rho_safe
        else
            rho1 = max(Y_safe*rho_safe/a_j, sgm_eps)
        end if

        pref_p = A*(1._wp - omega0*rho1/(R1*rho0))*exp(-R1*rho0/rho1) + B*(1._wp - omega0*rho1/(R2*rho0))*exp(-R2*rho0/rho1)
        K_p = rho0/max(omega0, sgm_eps)
        pref_u = jwl_unr_A*(1._wp - jwl_unr_omega*rho1/(jwl_unr_R1*jwl_unr_rho0))*exp(-jwl_unr_R1*jwl_unr_rho0/rho1) &
                            & + jwl_unr_B*(1._wp - jwl_unr_omega*rho1/(jwl_unr_R2*jwl_unr_rho0))*exp(-jwl_unr_R2*jwl_unr_rho0/rho1)
        K_u = jwl_unr_rho0/max(jwl_unr_omega, sgm_eps)

        pref1 = (1._wp - lam)*pref_u + lam*pref_p
        Kj = (1._wp - lam)*K_u + lam*K_p

        e = (pres*(a_j*Kj + a_a/max(air_gamma, sgm_eps)) - a_j*Kj*pref1)/rho_safe
        e = max(e, 0._wp)

    end subroutine s_jwl_reactive_energy_pr

    !> Lee-Tarver Ignition & Growth reaction rate dlambda/dt for progressive JWL burn. Three additive terms (ignition, growth,
    !! completion) gated by the lambda-window limiters jwl_lt_figmax/fg1max/fg2min. Coefficients (jwl_lt_*) and the unreacted
    !! reference density (jwl_unr_rho0) are read from the module globals. The caller supplies pressure and (condensed-phase)
    !! density; coefficient units must be consistent with the case's pressure units.
    pure function f_lee_tarver_rate(lambda, rho, pres) result(rate)

        $:GPU_ROUTINE(function_name='f_lee_tarver_rate',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in) :: lambda, rho, pres
        real(wp)             :: rate
        real(wp)             :: lam, one_m_lam, p_pos, comp, rho0u

        lam = min(max(lambda, 0._wp), 1._wp)
        one_m_lam = max(1._wp - lam, 0._wp)
        p_pos = max(pres, 0._wp)
        rho0u = max(jwl_unr_rho0, sgm_eps)
        rate = 0._wp

        ! Ignition term: hot-spot creation driven by compression (rho/rho0 - 1 - a)
        if (lam < jwl_lt_figmax) then
            comp = rho/rho0u - 1._wp - jwl_lt_a
            if (comp > 0._wp) then
                rate = rate + jwl_lt_I*one_m_lam**jwl_lt_b*comp**jwl_lt_x
            end if
        end if

        ! Growth term: pressure-driven burn spreading from hot spots
        if (lam < jwl_lt_fg1max .and. lam > 0._wp) then
            rate = rate + jwl_lt_G1*one_m_lam**jwl_lt_c*lam**jwl_lt_d*p_pos**jwl_lt_y
        end if

        ! Completion term: fast burnout at high reaction progress
        if (lam > jwl_lt_fg2min .and. one_m_lam > 0._wp) then
            rate = rate + jwl_lt_G2*one_m_lam**jwl_lt_e*lam**jwl_lt_g*p_pos**jwl_lt_z
        end if

        rate = max(rate, 0._wp)

    end function f_lee_tarver_rate

    !> Rocflu-style temperature from pressure and density.
    subroutine s_jwl_temperature_pr(rho, pres, e, Y, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma, cv, T)

        $:GPU_ROUTINE(function_name='s_jwl_temperature_pr',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)    :: rho, pres, e, Y, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma, cv
        real(wp), intent(inout) :: T
        real(wp)                :: rho_safe, Y_safe, cv_safe, T_air, T_jwl

        rho_safe = max(rho, sgm_eps)
        Y_safe = min(max(Y, 0._wp), 1._wp)
        cv_safe = max(cv, sgm_eps)

        ! Mass-weighted blend, consistent with the blended pressure EOS: ideal-gas air
        ! at Y->0 and the JWL thermal relation at Y->1, continuous in between.
        T_air = pres/max(air_gamma*rho_safe*cv_safe, sgm_eps)
        T_jwl = (pres - A*exp(-R1*rho0/rho_safe) - B*exp(-R2*rho0/rho_safe))/max(omega0*cv_safe*rho_safe, sgm_eps)

        T = (1._wp - Y_safe)*T_air + Y_safe*T_jwl
        if (T <= 0._wp .and. Y_safe > 1.e-2_wp) T = 270._wp

    end subroutine s_jwl_temperature_pr

    !> Compute the pressure from the appropriate equation of state
    subroutine s_compute_pressure(energy, alf, dyn_p, pi_inf, gamma, rho, qv, rhoYks, pres, T, stress, mom, G, pres_mag, jwl_Y, &
                                  & jwl_alpha)

        $:GPU_ROUTINE(function_name='s_compute_pressure',parallelism='[seq]', cray_noinline=True)

        real(stp), intent(in)           :: energy, alf
        real(wp), intent(in)            :: dyn_p
        real(wp), intent(in)            :: pi_inf, gamma, rho, qv
        real(wp), intent(out)           :: pres
        real(wp), intent(inout)         :: T
        real(stp), intent(in), optional :: stress, mom
        real(wp), intent(in), optional  :: G, pres_mag
        real(wp), intent(in), optional  :: jwl_Y
        real(wp), intent(in), optional  :: jwl_alpha

        ! Chemistry
        real(wp), dimension(1:num_species), intent(in) :: rhoYks
        real(wp), dimension(1:num_species)             :: Y_rs
        real(wp)                                       :: E_e
        real(wp)                                       :: e_Per_Kg, Pdyn_Per_Kg
        real(wp)                                       :: T_guess
        real(wp)                                       :: eint, e_sp, Y_jwl, alpha_jwl
        integer                                        :: s  !< Generic loop iterator
        #:if not chemistry
            ! Depending on model_eqns and bubbles_euler, the appropriate procedure for computing pressure is targeted by the
            ! procedure pointer

            if ((jwl_idx > 0) .and. (.not. mhd) .and. (model_eqns /= 4) .and. (bubbles_euler .neqv. .true.)) then
                Y_jwl = 1._wp
                if (present(jwl_Y)) Y_jwl = jwl_Y
                alpha_jwl = 1._wp
                if (present(jwl_alpha)) alpha_jwl = jwl_alpha
                eint = energy - dyn_p
                e_sp = eint/max(rho, sgm_eps)
                call s_jwl_pressure_er(rho, e_sp, Y_jwl, alpha_jwl, jwl_As(jwl_idx), jwl_Bs(jwl_idx), jwl_R1s(jwl_idx), &
                                       & jwl_R2s(jwl_idx), jwl_omegas(jwl_idx), jwl_rho0s(jwl_idx), jwl_E0s(jwl_idx), &
                                       & jwl_air_e0s(jwl_idx), jwl_air_rho0s(jwl_idx), jwl_air_gammas(jwl_idx), pres)
                if (pres < 1._wp) then
#ifdef MFC_SIMULATION
                    $:GPU_ATOMIC(atomic='update')
                    jwl_pressure_floor_count = jwl_pressure_floor_count + 1
#endif
                    pres = 1._wp
                end if
            else if (mhd) then
                ! MHD pressure: subtract magnetic pressure from total energy
                pres = (energy - dyn_p - pi_inf - qv - pres_mag)/gamma
            else if ((model_eqns /= 4) .and. (bubbles_euler .neqv. .true.)) then
                ! Gamma/pi_inf model or five-equation model (Allaire et al. JCP 2002): p from mixture EOS
                pres = (energy - dyn_p - pi_inf - qv)/gamma
            else if ((model_eqns /= 4) .and. bubbles_euler) then
                ! Bubble-augmented pressure with void fraction correction
                pres = ((energy - dyn_p)/(1._wp - alf) - pi_inf - qv)/gamma
            else
                ! Four-equation model (Kapila et al. PoF 2001): Tait EOS inversion
                pres = (pref + pi_inf)*(energy/(rhoref*(1 - alf)))**(1/gamma + 1) - pi_inf
            end if

            if (hypoelasticity .and. present(G)) then
                ! Subtract elastic strain energy before computing pressure (hypoelastic model)
                E_e = 0._wp
                do s = eqn_idx%stress%beg, eqn_idx%stress%end
                    if (G > 0) then
                        E_e = E_e + ((stress/rho)**2._wp)/(4._wp*G)
                        ! Double for shear stresses
                        if (any(s == shear_indices)) then
                            E_e = E_e + ((stress/rho)**2._wp)/(4._wp*G)
                        end if
                    end if
                end do

                if ((jwl_idx > 0) .and. (.not. mhd) .and. (model_eqns /= 4) .and. (bubbles_euler .neqv. .true.)) then
                    Y_jwl = 1._wp
                    if (present(jwl_Y)) Y_jwl = jwl_Y
                    alpha_jwl = 1._wp
                    if (present(jwl_alpha)) alpha_jwl = jwl_alpha
                    e_sp = (energy - 0.5_wp*(mom**2._wp)/rho - E_e)/max(rho, sgm_eps)
                    call s_jwl_pressure_er(rho, e_sp, Y_jwl, alpha_jwl, jwl_As(jwl_idx), jwl_Bs(jwl_idx), jwl_R1s(jwl_idx), &
                                           & jwl_R2s(jwl_idx), jwl_omegas(jwl_idx), jwl_rho0s(jwl_idx), jwl_E0s(jwl_idx), &
                                           & jwl_air_e0s(jwl_idx), jwl_air_rho0s(jwl_idx), jwl_air_gammas(jwl_idx), pres)
                else
                    pres = (energy - 0.5_wp*(mom**2._wp)/rho - pi_inf - qv - E_e)/gamma
                end if
            end if
        #:else
            ! Reacting mixture pressure from temperature and species
            Y_rs(:) = rhoYks(:)/rho
            e_Per_Kg = energy/rho
            Pdyn_Per_Kg = dyn_p/rho

            T_guess = T

            call get_temperature(e_Per_Kg - Pdyn_Per_Kg, T_guess, Y_rs, .true., T)
            call get_pressure(rho, T, Y_rs, pres)
        #:endif

    end subroutine s_compute_pressure

    !> Convert mixture variables to density, gamma, pi_inf, and qv for the gamma/pi_inf model. Given conservative or primitive
    !! variables, transfers the density, specific heat ratio function and the liquid stiffness function from q_vf to rho, gamma and
    !! pi_inf.
    subroutine s_convert_mixture_to_mixture_variables(q_vf, i, j, k, rho, gamma, pi_inf, qv)

        type(scalar_field), dimension(sys_size), intent(in) :: q_vf
        integer, intent(in)                                 :: i, j, k
        real(wp), intent(out), target                       :: rho
        real(wp), intent(out), target                       :: gamma
        real(wp), intent(out), target                       :: pi_inf
        real(wp), intent(out), target                       :: qv

        ! Transferring the density, the specific heat ratio function and the liquid stiffness function, respectively

        rho = q_vf(1)%sf(i, j, k)
        gamma = q_vf(eqn_idx%gamma)%sf(i, j, k)
        pi_inf = q_vf(eqn_idx%pi_inf)%sf(i, j, k)
        qv = 0._wp  ! keep this value nil for now. For future adjustment

        ! Post process requires rho_sf/gamma_sf/pi_inf_sf/qv_sf to also be updated
#ifdef MFC_POST_PROCESS
        rho_sf(i, j, k) = rho
        gamma_sf(i, j, k) = gamma
        pi_inf_sf(i, j, k) = pi_inf
        qv_sf(i, j, k) = qv
#endif

    end subroutine s_convert_mixture_to_mixture_variables

    !> Convert species volume fractions and partial densities to mixture density, gamma, pi_inf, and qv. Given conservative or
    !! primitive variables, computes the density, the specific heat ratio function and the liquid stiffness function from q_vf and
    !! stores the results into rho, gamma and pi_inf.
    subroutine s_convert_species_to_mixture_variables(q_vf, k, l, r, rho, gamma, pi_inf, qv, Re_K, G_K, G)

        type(scalar_field), dimension(sys_size), intent(in)   :: q_vf
        integer, intent(in)                                   :: k, l, r
        real(wp), intent(out), target                         :: rho
        real(wp), intent(out), target                         :: gamma
        real(wp), intent(out), target                         :: pi_inf
        real(wp), intent(out), target                         :: qv
        real(wp), optional, dimension(2), intent(out)         :: Re_K
        real(wp), optional, intent(out)                       :: G_K
        real(wp), dimension(num_fluids)                       :: alpha_rho_K, alpha_K
        real(wp), optional, dimension(num_fluids), intent(in) :: G
        integer                                               :: i, j  !< Generic loop iterator
        ! Computing the density, the specific heat ratio function and the liquid stiffness function, respectively

        call s_compute_species_fraction(q_vf, k, l, r, alpha_rho_K, alpha_K)

        ! Calculating the density, the specific heat ratio function, the liquid stiffness function, and the energy reference
        ! function, respectively, from the species analogs
        if (num_fluids == 1 .and. bubbles_euler) then
            rho = alpha_rho_K(1)
            gamma = gammas(1)
            pi_inf = pi_infs(1)
            qv = qvs(1)
        else
            rho = 0._wp; gamma = 0._wp; pi_inf = 0._wp; qv = 0._wp
            do i = 1, num_fluids
                rho = rho + alpha_rho_K(i)
                gamma = gamma + alpha_K(i)*gammas(i)
                pi_inf = pi_inf + alpha_K(i)*pi_infs(i)
                qv = qv + alpha_rho_K(i)*qvs(i)
            end do
        end if

#ifdef MFC_SIMULATION
        ! Computing the shear and bulk Reynolds numbers from species analogs
        if (viscous) then
            do i = 1, 2
                Re_K(i) = dflt_real; if (Re_size(i) > 0) Re_K(i) = 0._wp

                do j = 1, Re_size(i)
                    Re_K(i) = alpha_K(Re_idx(i, j))/fluid_pp(Re_idx(i, j))%Re(i) + Re_K(i)
                end do

                Re_K(i) = 1._wp/max(Re_K(i), sgm_eps)
            end do
        end if
#endif

        if (present(G_K)) then
            G_K = 0._wp
            do i = 1, num_fluids
                G_K = G_K + alpha_K(i)*G(i)
            end do
            G_K = max(0._wp, G_K)
        end if

        ! Post process requires rho_sf/gamma_sf/pi_inf_sf/qv_sf to also be updated
#ifdef MFC_POST_PROCESS
        rho_sf(k, l, r) = rho
        gamma_sf(k, l, r) = gamma
        pi_inf_sf(k, l, r) = pi_inf
        qv_sf(k, l, r) = qv
#endif

    end subroutine s_convert_species_to_mixture_variables

    !> GPU-accelerated conversion of species volume fractions and partial densities to mixture density, gamma, pi_inf, and qv.
    subroutine s_convert_species_to_mixture_variables_acc(rho_K, gamma_K, pi_inf_K, qv_K, alpha_K, alpha_rho_K, Re_K, G_K, G)

        $:GPU_ROUTINE(function_name='s_convert_species_to_mixture_variables_acc', parallelism='[seq]', cray_noinline=True)

        real(wp), intent(out) :: rho_K, gamma_K, pi_inf_K, qv_K
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(inout)        :: alpha_rho_K, alpha_K
            real(wp), optional, dimension(3), intent(in) :: G
        #:else
            real(wp), dimension(num_fluids), intent(inout)        :: alpha_rho_K, alpha_K
            real(wp), optional, dimension(num_fluids), intent(in) :: G
        #:endif
        real(wp), dimension(2), intent(out) :: Re_K
        real(wp), optional, intent(out)     :: G_K
        real(wp)                            :: alpha_K_sum
        integer                             :: i, j  !< Generic loop iterators

        rho_K = 0._wp
        gamma_K = 0._wp
        pi_inf_K = 0._wp
        qv_K = 0._wp
        Re_K = dflt_real
        if (present(G_K)) G_K = 0._wp

#ifdef MFC_SIMULATION
        ! Constrain partial densities and volume fractions within physical bounds
        if (num_fluids == 1 .and. bubbles_euler) then
            rho_K = alpha_rho_K(1)
            gamma_K = gammas(1)
            pi_inf_K = pi_infs(1)
            qv_K = qvs(1)
        else
            if (mpp_lim) then
                alpha_K_sum = 0._wp
                do i = 1, num_fluids
                    alpha_rho_K(i) = max(0._wp, alpha_rho_K(i))
                    alpha_K(i) = min(max(0._wp, alpha_K(i)), 1._wp)
                    alpha_K_sum = alpha_K_sum + alpha_K(i)
                end do
                alpha_K = alpha_K/max(alpha_K_sum, sgm_eps)
            end if
            rho_K = 0._wp; gamma_K = 0._wp; pi_inf_K = 0._wp; qv_K = 0._wp
            do i = 1, num_fluids
                rho_K = rho_K + alpha_rho_K(i)
                gamma_K = gamma_K + alpha_K(i)*gammas(i)
                pi_inf_K = pi_inf_K + alpha_K(i)*pi_infs(i)
                qv_K = qv_K + alpha_rho_K(i)*qvs(i)
            end do
        end if

        if (present(G_K)) then
            G_K = 0._wp
            do i = 1, num_fluids
                ! TODO: change to use Gs_vc directly here? TODO: Make this change as well for GPUs
                G_K = G_K + alpha_K(i)*G(i)
            end do
            G_K = max(0._wp, G_K)
        end if

        if (viscous) then
            do i = 1, 2
                Re_K(i) = dflt_real

                if (Re_size(i) > 0) Re_K(i) = 0._wp

                do j = 1, Re_size(i)
                    Re_K(i) = alpha_K(Re_idx(i, j))/Res_vc(i, j) + Re_K(i)
                end do

                Re_K(i) = 1._wp/max(Re_K(i), sgm_eps)
            end do
        end if
#endif

    end subroutine s_convert_species_to_mixture_variables_acc

    !> Initialize the variables conversion module.
    impure subroutine s_initialize_variables_conversion_module

        integer :: i, j, num_jwl

        $:GPU_ENTER_DATA(copyin='[is1b, is1e, is2b, is2e, is3b, is3e]')

        @:ALLOCATE(eos_idxs(1:num_fluids))
        @:ALLOCATE(gammas (1:num_fluids))
        @:ALLOCATE(gs_min (1:num_fluids))
        @:ALLOCATE(pi_infs(1:num_fluids))
        @:ALLOCATE(ps_inf(1:num_fluids))
        @:ALLOCATE(cvs    (1:num_fluids))
        @:ALLOCATE(qvs    (1:num_fluids))
        @:ALLOCATE(qvps    (1:num_fluids))
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
        @:ALLOCATE(Gs_vc     (1:num_fluids))

        jwl_idx = 0
        jwl_pressure_floor_count = 0
        num_jwl = 0
        do i = 1, num_fluids
            eos_idxs(i) = fluid_pp(i)%eos
            gammas(i) = fluid_pp(i)%gamma
            gs_min(i) = 1.0_wp/gammas(i) + 1.0_wp
            pi_infs(i) = fluid_pp(i)%pi_inf
            Gs_vc(i) = fluid_pp(i)%G
            ps_inf(i) = pi_infs(i)/(1.0_wp + gammas(i))
            cvs(i) = fluid_pp(i)%cv
            qvs(i) = fluid_pp(i)%qv
            qvps(i) = fluid_pp(i)%qvp
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
            if (eos_idxs(i) == 2) then
                jwl_idx = i
                num_jwl = num_jwl + 1
            end if
        end do

        if (num_jwl > 1) call s_mpi_abort('JWL EOS currently supports one JWL fluid per case.')

        do i = 1, num_fluids
            if (eos_idxs(i) == 2) then
                if (mhd) call s_mpi_abort('JWL EOS is not currently supported with MHD.')
                if (bubbles_euler) call s_mpi_abort('JWL EOS is not currently supported with Eulerian bubbles.')
                if (chemistry) call s_mpi_abort('JWL EOS is not currently supported with chemistry.')
                if (model_eqns == 4) call s_mpi_abort('JWL EOS is not currently supported with model_eqns = 4.')
                if (jwl_R1s(i) <= 0._wp .or. jwl_R2s(i) <= 0._wp .or. jwl_omegas(i) <= 0._wp .or. jwl_rho0s(i) <= 0._wp) then
                    call s_mpi_abort('JWL EOS requires positive R1, R2, omega, and rho0 parameters.')
                end if
                if (num_fluids > 1 .and. jwl_E0s(i) <= 0._wp) then
                    call s_mpi_abort('Mixed-material JWL EOS requires positive fluid_pp%jwl_E0.')
                end if
                if (jwl_air_e0s(i) <= 0._wp .or. jwl_air_rho0s(i) <= 0._wp .or. jwl_air_gammas(i) <= 0._wp) then
                    call s_mpi_abort('JWL EOS requires positive air reference energy, density, and gamma-minus-one.')
                end if
            else if (eos_idxs(i) /= 1) then
                call s_mpi_abort('Unsupported fluid_pp%eos selector. Use 1 for stiffened gas or 2 for JWL.')
            end if
        end do

        $:GPU_UPDATE(device='[eos_idxs, gammas, gs_min, pi_infs, ps_inf, cvs, qvs, qvps, Gs_vc]')
        $:GPU_UPDATE(device='[jwl_As, jwl_Bs, jwl_R1s, jwl_R2s, jwl_omegas, jwl_rho0s, jwl_E0s]')
        $:GPU_UPDATE(device='[jwl_air_e0s, jwl_air_rho0s, jwl_air_gammas, jwl_idx, jwl_pressure_floor_count]')

#ifdef MFC_SIMULATION
        if (viscous) then
            @:ALLOCATE(Res_vc(1:2, 1:Re_size_max))
            do i = 1, 2
                do j = 1, Re_size(i)
                    Res_vc(i, j) = fluid_pp(Re_idx(i, j))%Re(i)
                end do
            end do

            $:GPU_UPDATE(device='[Res_vc, Re_idx, Re_size]')
        end if
#endif

        if (bubbles_euler) then
            @:ALLOCATE(bubrs_vc(1:nb))
            do i = 1, nb
                bubrs_vc(i) = qbmm_idx%rs(i)
            end do
            $:GPU_UPDATE(device='[bubrs_vc]')
        end if

#ifdef MFC_POST_PROCESS
        ! Allocating the density, the specific heat ratio function and the liquid stiffness function, respectively

        ! Simulation is at least 2D
        if (n > 0) then
            ! Simulation is 3D
            if (p > 0) then
                allocate (rho_sf(-buff_size:m + buff_size,-buff_size:n + buff_size,-buff_size:p + buff_size))
                allocate (gamma_sf(-buff_size:m + buff_size,-buff_size:n + buff_size,-buff_size:p + buff_size))
                allocate (pi_inf_sf(-buff_size:m + buff_size,-buff_size:n + buff_size,-buff_size:p + buff_size))
                allocate (qv_sf(-buff_size:m + buff_size,-buff_size:n + buff_size,-buff_size:p + buff_size))

                ! Simulation is 2D
            else
                allocate (rho_sf(-buff_size:m + buff_size,-buff_size:n + buff_size,0:0))
                allocate (gamma_sf(-buff_size:m + buff_size,-buff_size:n + buff_size,0:0))
                allocate (pi_inf_sf(-buff_size:m + buff_size,-buff_size:n + buff_size,0:0))
                allocate (qv_sf(-buff_size:m + buff_size,-buff_size:n + buff_size,0:0))
            end if

            ! Simulation is 1D
        else
            allocate (rho_sf(-buff_size:m + buff_size,0:0,0:0))
            allocate (gamma_sf(-buff_size:m + buff_size,0:0,0:0))
            allocate (pi_inf_sf(-buff_size:m + buff_size,0:0,0:0))
            allocate (qv_sf(-buff_size:m + buff_size,0:0,0:0))
        end if
#endif

    end subroutine s_initialize_variables_conversion_module

    !> Initialize bubble mass-vapor values at quadrature nodes from the conserved moment statistics.
    subroutine s_initialize_mv(qK_cons_vf, mv)

        type(scalar_field), dimension(sys_size), intent(in)                                     :: qK_cons_vf
        real(stp), dimension(idwint(1)%beg:,idwint(2)%beg:,idwint(3)%beg:,1:,1:), intent(inout) :: mv
        integer                                                                                 :: i, j, k, l
        real(wp)                                                                                :: mu, sig, nbub_sc

        do l = idwint(3)%beg, idwint(3)%end
            do k = idwint(2)%beg, idwint(2)%end
                do j = idwint(1)%beg, idwint(1)%end
                    nbub_sc = qK_cons_vf(eqn_idx%bub%beg)%sf(j, k, l)

                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, nb
                        mu = qK_cons_vf(eqn_idx%bub%beg + 1 + (i - 1)*nmom)%sf(j, k, l)/nbub_sc
                        sig = (qK_cons_vf(eqn_idx%bub%beg + 3 + (i - 1)*nmom)%sf(j, k, l)/nbub_sc - mu**2)**0.5_wp

                        mv(j, k, l, 1, i) = (mass_v0(i))*(mu - sig)**(3._wp)/(R0(i)**(3._wp))
                        mv(j, k, l, 2, i) = (mass_v0(i))*(mu - sig)**(3._wp)/(R0(i)**(3._wp))
                        mv(j, k, l, 3, i) = (mass_v0(i))*(mu + sig)**(3._wp)/(R0(i)**(3._wp))
                        mv(j, k, l, 4, i) = (mass_v0(i))*(mu + sig)**(3._wp)/(R0(i)**(3._wp))
                    end do
                end do
            end do
        end do

    end subroutine s_initialize_mv

    !> Initialize bubble internal pressures at quadrature nodes using isothermal relations from the Preston model.
    subroutine s_initialize_pb(qK_cons_vf, mv, pb)

        type(scalar_field), dimension(sys_size), intent(in)                                     :: qK_cons_vf
        real(stp), dimension(idwint(1)%beg:,idwint(2)%beg:,idwint(3)%beg:,1:,1:), intent(in)    :: mv
        real(stp), dimension(idwint(1)%beg:,idwint(2)%beg:,idwint(3)%beg:,1:,1:), intent(inout) :: pb
        integer                                                                                 :: i, j, k, l
        real(wp)                                                                                :: mu, sig, nbub_sc

        do l = idwint(3)%beg, idwint(3)%end
            do k = idwint(2)%beg, idwint(2)%end
                do j = idwint(1)%beg, idwint(1)%end
                    nbub_sc = qK_cons_vf(eqn_idx%bub%beg)%sf(j, k, l)

                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, nb
                        mu = qK_cons_vf(eqn_idx%bub%beg + 1 + (i - 1)*nmom)%sf(j, k, l)/nbub_sc
                        sig = (qK_cons_vf(eqn_idx%bub%beg + 3 + (i - 1)*nmom)%sf(j, k, l)/nbub_sc - mu**2)**0.5_wp

                        ! PRESTON (ISOTHERMAL)
                        pb(j, k, l, 1, i) = (pb0(i))*(R0(i)**(3._wp))*(mass_g0(i) + mv(j, k, l, 1, &
                           & i))/(mu - sig)**(3._wp)/(mass_g0(i) + mass_v0(i))
                        pb(j, k, l, 2, i) = (pb0(i))*(R0(i)**(3._wp))*(mass_g0(i) + mv(j, k, l, 2, &
                           & i))/(mu - sig)**(3._wp)/(mass_g0(i) + mass_v0(i))
                        pb(j, k, l, 3, i) = (pb0(i))*(R0(i)**(3._wp))*(mass_g0(i) + mv(j, k, l, 3, &
                           & i))/(mu + sig)**(3._wp)/(mass_g0(i) + mass_v0(i))
                        pb(j, k, l, 4, i) = (pb0(i))*(R0(i)**(3._wp))*(mass_g0(i) + mv(j, k, l, 4, &
                           & i))/(mu + sig)**(3._wp)/(mass_g0(i) + mass_v0(i))
                    end do
                end do
            end do
        end do

    end subroutine s_initialize_pb

    !> Convert conserved variables (rho*alpha, rho*u, E, alpha) to primitives (rho, u, p, alpha). Conversion depends on model_eqns:
    !! each model has different variable sets and EOS.
    subroutine s_convert_conservative_to_primitive_variables(qK_cons_vf, q_T_sf, qK_prim_vf, ibounds)

        type(scalar_field), dimension(sys_size), intent(in)    :: qK_cons_vf
        type(scalar_field), intent(inout)                      :: q_T_sf
        type(scalar_field), dimension(sys_size), intent(inout) :: qK_prim_vf
        type(int_bounds_info), dimension(1:3), intent(in)      :: ibounds

        #:if USING_AMD and not MFC_CASE_OPTIMIZATION
            real(wp), dimension(3) :: alpha_K, alpha_rho_K
            real(wp), dimension(3) :: nRtmp
            real(wp)               :: rhoYks(1:10)
        #:else
            real(wp), dimension(num_fluids) :: alpha_K, alpha_rho_K
            real(wp), dimension(nb)         :: nRtmp
            real(wp)                        :: rhoYks(1:num_species)
        #:endif
        real(wp), dimension(2) :: Re_K
        real(wp)               :: rho_K, gamma_K, pi_inf_K, qv_K, dyn_pres_K
        real(wp)               :: Y_jwl, alpha_jwl
        real(wp)               :: vftmp, nbub_sc
        real(wp)               :: G_K
        real(wp)               :: pres
        integer                :: i, j, k, l               !< Generic loop iterators
        real(wp)               :: T
        real(wp)               :: pres_mag
        real(wp)               :: Ga                       !< Lorentz factor (gamma in relativity)
        real(wp)               :: B2                       !< Magnetic field magnitude squared
        real(wp)               :: B(3)                     !< Magnetic field components
        real(wp)               :: m2                       !< Relativistic momentum magnitude squared
        real(wp)               :: S                        !< Dot product of the magnetic field and the relativistic momentum
        real(wp)               :: W, dW                    !< W := rho*v*Ga**2; f = f(W) in Newton-Raphson
        real(wp)               :: E, D                     !< Prim/Cons variables within Newton-Raphson iteration
        real(wp)               :: f, dGa_dW, dp_dW, df_dW  !< Functions within Newton-Raphson iteration
        integer                :: iter                     !< Newton-Raphson iteration counter

        $:GPU_PARALLEL_LOOP(collapse=3, private='[alpha_K, alpha_rho_K, Re_K, nRtmp, rho_K, gamma_K, pi_inf_K, qv_K, dyn_pres_K, &
                            & rhoYks, B, pres, vftmp, nbub_sc, G_K, T, pres_mag, Ga, B2, m2, S, W, dW, E, D, f, dGa_dW, dp_dW, &
                            & df_dW, iter, Y_jwl, alpha_jwl]')
        do l = ibounds(3)%beg, ibounds(3)%end
            do k = ibounds(2)%beg, ibounds(2)%end
                do j = ibounds(1)%beg, ibounds(1)%end
                    dyn_pres_K = 0._wp

                    call s_compute_species_fraction(qK_cons_vf, j, k, l, alpha_rho_K, alpha_K)

                    if (model_eqns /= 4) then
#ifdef MFC_SIMULATION
                        ! If in simulation, use acc mixture subroutines
                        if (elasticity) then
                            call s_convert_species_to_mixture_variables_acc(rho_K, gamma_K, pi_inf_K, qv_K, alpha_K, alpha_rho_K, &
                                & Re_K, G_K, Gs_vc)
                        else
                            call s_convert_species_to_mixture_variables_acc(rho_K, gamma_K, pi_inf_K, qv_K, alpha_K, alpha_rho_K, &
                                & Re_K)
                        end if
#else
                        ! If pre-processing, use non acc mixture subroutines
                        if (elasticity) then
                            call s_convert_to_mixture_variables(qK_cons_vf, j, k, l, rho_K, gamma_K, pi_inf_K, qv_K, Re_K, G_K, &
                                                                & fluid_pp(:)%G)
                        else
                            call s_convert_to_mixture_variables(qK_cons_vf, j, k, l, rho_K, gamma_K, pi_inf_K, qv_K)
                        end if
#endif
                    end if

                    ! Relativistic MHD primitive variable recovery, Mignone & Bodo A&A (2006)
                    if (relativity) then
                        if (n == 0) then
                            B(1) = Bx0
                            B(2) = qK_cons_vf(eqn_idx%B%beg)%sf(j, k, l)
                            B(3) = qK_cons_vf(eqn_idx%B%beg + 1)%sf(j, k, l)
                        else
                            B(1) = qK_cons_vf(eqn_idx%B%beg)%sf(j, k, l)
                            B(2) = qK_cons_vf(eqn_idx%B%beg + 1)%sf(j, k, l)
                            B(3) = qK_cons_vf(eqn_idx%B%beg + 2)%sf(j, k, l)
                        end if
                        B2 = B(1)**2 + B(2)**2 + B(3)**2

                        m2 = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%mom%beg, eqn_idx%mom%end
                            m2 = m2 + qK_cons_vf(i)%sf(j, k, l)**2
                        end do

                        S = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, 3
                            S = S + qK_cons_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l)*B(i)
                        end do

                        E = qK_cons_vf(eqn_idx%E)%sf(j, k, l)

                        D = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, eqn_idx%cont%end
                            D = D + qK_cons_vf(i)%sf(j, k, l)
                        end do

                        ! Newton-Raphson
                        W = E + D
                        $:GPU_LOOP(parallelism='[seq]')
                        do iter = 1, relativity_cons_to_prim_max_iter
                            ! Lorentz factor from total enthalpy and magnetic field
                            Ga = (W + B2)*W/sqrt((W + B2)**2*W**2 - (m2*W**2 + S**2*(2*W + B2)))
                            ! Thermal pressure from EOS
                            pres = (W - D*Ga)/((gamma_K + 1)*Ga**2)
                            f = W - pres + (1 - 1/(2*Ga**2))*B2 - S**2/(2*W**2) - E - D

                            ! The first equation below corrects a typo in (Mignone & Bodo, 2006) m2*W**2 -> 2*m2*W**2, which would
                            ! cancel with the 2* in other terms This corrected version is not used as the second equation
                            ! empirically converges faster. First equation is kept for further investigation. dGa_dW = -Ga**3 * (
                            ! S**2*(3*W**2+3*W*B2+B2**2) + m2*W**2 ) / (W**3 * (W+B2)**3) ! first (corrected)
                            dGa_dW = -Ga**3*(2*S**2*(3*W**2 + 3*W*B2 + B2**2) + m2*W**2)/(2*W**3*(W + B2)**3)  ! second (in paper)

                            dp_dW = (Ga*(1 + D*dGa_dW) - 2*W*dGa_dW)/((gamma_K + 1)*Ga**3)
                            df_dW = 1 - dp_dW + (B2/Ga**3)*dGa_dW + S**2/W**3

                            dW = -f/df_dW
                            W = W + dW
                            if (abs(dW) < 1.e-12_wp*W) exit  ! Relative convergence criterion
                        end do

                        ! Recalculate pressure using converged W
                        Ga = (W + B2)*W/sqrt((W + B2)**2*W**2 - (m2*W**2 + S**2*(2*W + B2)))
                        qK_prim_vf(eqn_idx%E)%sf(j, k, l) = (W - D*Ga)/((gamma_K + 1)*Ga**2)

                        ! Recover the other primitive variables
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, 3
                            qK_prim_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l) = (qK_cons_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, &
                                       & l) + (S/W)*B(i))/(W + B2)
                        end do
                        qK_prim_vf(1)%sf(j, k, l) = D/Ga  ! Hard-coded for single-component for now

                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%B%beg, eqn_idx%B%end
                            qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)
                        end do

                        cycle  ! skip all the non-relativistic conversions below
                    end if

                    if (chemistry) then
                        ! Reacting flow: recover density from species partial densities, compute mass fractions Y_k = rhoY_k / rho
                        rho_K = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%species%beg, eqn_idx%species%end
                            rho_K = rho_K + max(0._wp, qK_cons_vf(i)%sf(j, k, l))
                        end do

                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, eqn_idx%cont%end
                            qK_prim_vf(i)%sf(j, k, l) = rho_K
                        end do

                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%species%beg, eqn_idx%species%end
                            qK_prim_vf(i)%sf(j, k, l) = max(0._wp, qK_cons_vf(i)%sf(j, k, l)/rho_K)
                        end do
                    else
                        ! Non-reacting: partial densities are directly primitive (alpha_i * rho_i)
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, eqn_idx%cont%end
                            qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)
                        end do
                    end if

#ifdef MFC_SIMULATION
                    rho_K = max(rho_K, sgm_eps)
#endif

                    ! Recover velocity from momentum: u = rho*u / rho, and accumulate dynamic pressure 0.5*rho*|u|^2
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = eqn_idx%mom%beg, eqn_idx%mom%end
                        if (model_eqns /= 4) then
                            qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)/rho_K
                            dyn_pres_K = dyn_pres_K + 5.e-1_wp*qK_cons_vf(i)%sf(j, k, l)*qK_prim_vf(i)%sf(j, k, l)
                        else
                            ! Four-equation model (Kapila et al. PoF 2001): divide by total density q_cons(1)
                            qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)/qK_cons_vf(1)%sf(j, k, l)
                        end if
                    end do

                    if (chemistry) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, num_species
                            rhoYks(i) = qK_cons_vf(eqn_idx%species%beg + i - 1)%sf(j, k, l)
                        end do

                        T = q_T_sf%sf(j, k, l)
                    end if

                    if (mhd) then
                        if (n == 0) then
                            pres_mag = 0.5_wp*(Bx0**2 + qK_cons_vf(eqn_idx%B%beg)%sf(j, k, &
                                               & l)**2 + qK_cons_vf(eqn_idx%B%beg + 1)%sf(j, k, l)**2)
                        else
                            pres_mag = 0.5_wp*(qK_cons_vf(eqn_idx%B%beg)%sf(j, k, l)**2 + qK_cons_vf(eqn_idx%B%beg + 1)%sf(j, k, &
                                               & l)**2 + qK_cons_vf(eqn_idx%B%beg + 2)%sf(j, k, l)**2)
                        end if
                    else
                        pres_mag = 0._wp
                    end if

                    Y_jwl = 1._wp
                    alpha_jwl = 1._wp
                    if (jwl_idx > 0 .and. jwl_idx <= eqn_idx%cont%end) then
                        Y_jwl = qK_cons_vf(jwl_idx)%sf(j, k, l)/max(rho_K, sgm_eps)
                        alpha_jwl = qK_cons_vf(eqn_idx%adv%beg + jwl_idx - 1)%sf(j, k, l)
                    end if

                    call s_compute_pressure(qK_cons_vf(eqn_idx%E)%sf(j, k, l), qK_cons_vf(eqn_idx%alf)%sf(j, k, l), dyn_pres_K, &
                                            & pi_inf_K, gamma_K, rho_K, qv_K, rhoYks, pres, T, pres_mag=pres_mag, jwl_Y=Y_jwl, &
                                            & jwl_alpha=alpha_jwl)

                    qK_prim_vf(eqn_idx%E)%sf(j, k, l) = pres

                    if (chemistry) then
                        q_T_sf%sf(j, k, l) = T
                    end if

                    if (bubbles_euler) then
                        ! Recover bubble primitive variables: divide conserved moments by bubble number density
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, nb
                            nRtmp(i) = qK_cons_vf(bubrs_vc(i))%sf(j, k, l)
                        end do

                        vftmp = qK_cons_vf(eqn_idx%alf)%sf(j, k, l)

                        if (qbmm) then
                            ! Get nb (constant across all R0 bins)
                            nbub_sc = qK_cons_vf(eqn_idx%bub%beg)%sf(j, k, l)

                            ! Convert cons to prim
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = eqn_idx%bub%beg, eqn_idx%bub%end
                                qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)/nbub_sc
                            end do
                            ! Need to keep track of nb in the primitive variable list (converted back to true value before output)
#ifdef MFC_SIMULATION
                            qK_prim_vf(eqn_idx%bub%beg)%sf(j, k, l) = qK_cons_vf(eqn_idx%bub%beg)%sf(j, k, l)
#endif
                        else
                            if (adv_n) then
                                qK_prim_vf(eqn_idx%n)%sf(j, k, l) = qK_cons_vf(eqn_idx%n)%sf(j, k, l)
                                nbub_sc = qK_prim_vf(eqn_idx%n)%sf(j, k, l)
                            else
                                call s_comp_n_from_cons(vftmp, nRtmp, nbub_sc, weight)
                            end if

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = eqn_idx%bub%beg, eqn_idx%bub%end
                                qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)/nbub_sc
                            end do
                        end if
                    end if

                    if (mhd) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%B%beg, eqn_idx%B%end
                            qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)
                        end do
                    end if

                    if (elasticity) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%stress%beg, eqn_idx%stress%end
                            qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)/rho_K
                        end do
                    end if

                    if (hypoelasticity) then
                        if (cont_damage) G_K = G_K*max((1._wp - qK_cons_vf(eqn_idx%damage)%sf(j, k, l)), 0._wp)
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%stress%beg, eqn_idx%stress%end
                            ! subtracting elastic contribution for pressure calculation
                            if (G_K > verysmall) then
                                qK_prim_vf(eqn_idx%E)%sf(j, k, l) = qK_prim_vf(eqn_idx%E)%sf(j, k, l) - ((qK_prim_vf(i)%sf(j, k, &
                                           & l)**2._wp)/(4._wp*G_K))/gamma_K
                                ! Double for shear stresses
                                if (any(i == shear_indices)) then
                                    qK_prim_vf(eqn_idx%E)%sf(j, k, l) = qK_prim_vf(eqn_idx%E)%sf(j, k, l) - ((qK_prim_vf(i)%sf(j, &
                                               & k, l)**2._wp)/(4._wp*G_K))/gamma_K
                                end if
                            end if
                        end do
                    end if

                    if (hyperelasticity) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%xi%beg, eqn_idx%xi%end
                            qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)/rho_K
                        end do
                    end if

                    if (.not. igr .or. num_fluids > 1) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%adv%beg, eqn_idx%adv%end
                            qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)
                        end do
                    end if

                    if (surface_tension) then
                        qK_prim_vf(eqn_idx%c)%sf(j, k, l) = qK_cons_vf(eqn_idx%c)%sf(j, k, l)
                    end if

                    if (cont_damage) qK_prim_vf(eqn_idx%damage)%sf(j, k, l) = qK_cons_vf(eqn_idx%damage)%sf(j, k, l)

                    if (hyper_cleaning) qK_prim_vf(eqn_idx%psi)%sf(j, k, l) = qK_cons_vf(eqn_idx%psi)%sf(j, k, l)
#ifdef MFC_POST_PROCESS
                    if (bubbles_lagrange) qK_prim_vf(beta_idx)%sf(j, k, l) = qK_cons_vf(beta_idx)%sf(j, k, l)
#endif
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_convert_conservative_to_primitive_variables

    !> Convert primitives (rho, u, p, alpha) to conserved variables (rho*alpha, rho*u, E, alpha).
    impure subroutine s_convert_primitive_to_conservative_variables(q_prim_vf, q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(in)    :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf

        ! Density, specific heat ratio function, liquid stiffness function and dynamic pressure, as defined in the incompressible
        ! flow sense, respectively
        real(wp)                         :: rho
        real(wp)                         :: gamma
        real(wp)                         :: pi_inf
        real(wp)                         :: qv
        real(wp)                         :: dyn_pres
        real(wp)                         :: nbub, R3tmp
        real(wp), dimension(nb)          :: Rtmp
        real(wp)                         :: G
        real(wp), dimension(2)           :: Re_K
        integer                          :: i, j, k, l  !< Generic loop iterators
        real(wp), dimension(num_species) :: Ys
        real(wp)                         :: e_mix, mix_mol_weight, T
        real(wp)                         :: pres_mag
        real(wp)                         :: pcold, e_sp, Y_jwl, alpha_jwl, phase_energy_sum, rho_phase
        real(wp)                         :: Ga          !< Lorentz factor (gamma in relativity)
        real(wp)                         :: h           !< relativistic enthalpy
        real(wp)                         :: v2          !< Square of the velocity magnitude
        real(wp)                         :: B2          !< Square of the magnetic field magnitude
        real(wp)                         :: vdotB       !< Dot product of the velocity and magnetic field vectors
        real(wp)                         :: B(3)        !< Magnetic field components

        pres_mag = 0._wp

        G = 0._wp

#ifndef MFC_SIMULATION
        ! Converting the primitive variables to the conservative variables
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    ! Obtaining the density, specific heat ratio function and the liquid stiffness function, respectively
                    call s_convert_to_mixture_variables(q_prim_vf, j, k, l, rho, gamma, pi_inf, qv, Re_K, G, fluid_pp(:)%G)

                    if (.not. igr .or. num_fluids > 1) then
                        ! Transferring the advection equation(s) variable(s)
                        do i = eqn_idx%adv%beg, eqn_idx%adv%end
                            q_cons_vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)
                        end do
                    end if

                    if (relativity) then
                        if (n == 0) then
                            B(1) = Bx0
                            B(2) = q_prim_vf(eqn_idx%B%beg)%sf(j, k, l)
                            B(3) = q_prim_vf(eqn_idx%B%beg + 1)%sf(j, k, l)
                        else
                            B(1) = q_prim_vf(eqn_idx%B%beg)%sf(j, k, l)
                            B(2) = q_prim_vf(eqn_idx%B%beg + 1)%sf(j, k, l)
                            B(3) = q_prim_vf(eqn_idx%B%beg + 2)%sf(j, k, l)
                        end if

                        v2 = 0._wp
                        do i = eqn_idx%mom%beg, eqn_idx%mom%end
                            v2 = v2 + q_prim_vf(i)%sf(j, k, l)**2
                        end do
                        if (v2 >= 1._wp) call s_mpi_abort('Error: v squared > 1 in s_convert_primitive_to_conservative_variables')

                        Ga = 1._wp/sqrt(1._wp - v2)

                        h = 1._wp + (gamma + 1)*q_prim_vf(eqn_idx%E)%sf(j, k, l)/rho  ! Assume perfect gas for now

                        B2 = 0._wp
                        do i = eqn_idx%B%beg, eqn_idx%B%end
                            B2 = B2 + q_prim_vf(i)%sf(j, k, l)**2
                        end do
                        if (n == 0) B2 = B2 + Bx0**2

                        vdotB = 0._wp
                        do i = 1, 3
                            vdotB = vdotB + q_prim_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l)*B(i)
                        end do

                        do i = 1, eqn_idx%cont%end
                            q_cons_vf(i)%sf(j, k, l) = Ga*q_prim_vf(i)%sf(j, k, l)
                        end do

                        do i = eqn_idx%mom%beg, eqn_idx%mom%end
                            q_cons_vf(i)%sf(j, k, l) = (rho*h*Ga**2 + B2)*q_prim_vf(i)%sf(j, k, &
                                      & l) - vdotB*B(i - eqn_idx%mom%beg + 1)
                        end do

                        q_cons_vf(eqn_idx%E)%sf(j, k, l) = rho*h*Ga**2 - q_prim_vf(eqn_idx%E)%sf(j, k, &
                                  & l) + 0.5_wp*(B2 + v2*B2 - vdotB**2)
                        ! Remove rest energy
                        do i = 1, eqn_idx%cont%end
                            q_cons_vf(eqn_idx%E)%sf(j, k, l) = q_cons_vf(eqn_idx%E)%sf(j, k, l) - q_cons_vf(i)%sf(j, k, l)
                        end do

                        do i = eqn_idx%B%beg, eqn_idx%B%end
                            q_cons_vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)
                        end do

                        cycle  ! skip all the non-relativistic conversions below
                    end if

                    ! Transferring the continuity equation(s) variable(s)
                    do i = 1, eqn_idx%cont%end
                        q_cons_vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)
                    end do

                    ! Zeroing out the dynamic pressure since it is computed iteratively by cycling through the velocity equations
                    dyn_pres = 0._wp

                    ! Computing momenta and dynamic pressure from velocity
                    do i = eqn_idx%mom%beg, eqn_idx%mom%end
                        q_cons_vf(i)%sf(j, k, l) = rho*q_prim_vf(i)%sf(j, k, l)
                        dyn_pres = dyn_pres + q_cons_vf(i)%sf(j, k, l)*q_prim_vf(i)%sf(j, k, l)/2._wp
                    end do

                    if (chemistry) then
                        ! Reacting mixture: compute conserved energy from species mass fractions and temperature
                        do i = eqn_idx%species%beg, eqn_idx%species%end
                            Ys(i - eqn_idx%species%beg + 1) = q_prim_vf(i)%sf(j, k, l)
                            q_cons_vf(i)%sf(j, k, l) = rho*q_prim_vf(i)%sf(j, k, l)
                        end do

                        call get_mixture_molecular_weight(Ys, mix_mol_weight)
                        T = q_prim_vf(eqn_idx%E)%sf(j, k, l)*mix_mol_weight/(gas_constant*rho)
                        call get_mixture_energy_mass(T, Ys, e_mix)

                        q_cons_vf(eqn_idx%E)%sf(j, k, l) = dyn_pres + rho*e_mix
                    else
                        ! Computing the energy from the pressure
                        if ((jwl_idx > 0) .and. (.not. mhd) .and. (model_eqns /= 4) .and. (bubbles_euler .neqv. .true.)) then
                            Y_jwl = 1._wp
                            alpha_jwl = 1._wp
                            if (jwl_idx <= eqn_idx%cont%end) then
                                Y_jwl = q_prim_vf(jwl_idx)%sf(j, k, l)/max(rho, sgm_eps)
                                alpha_jwl = q_prim_vf(eqn_idx%adv%beg + jwl_idx - 1)%sf(j, k, l)
                            end if
                            call s_jwl_energy_pr(rho, q_prim_vf(eqn_idx%E)%sf(j, k, l), Y_jwl, alpha_jwl, jwl_As(jwl_idx), &
                                                 & jwl_Bs(jwl_idx), jwl_R1s(jwl_idx), jwl_R2s(jwl_idx), jwl_omegas(jwl_idx), &
                                                 & jwl_rho0s(jwl_idx), jwl_E0s(jwl_idx), jwl_air_e0s(jwl_idx), &
                                                 & jwl_air_rho0s(jwl_idx), jwl_air_gammas(jwl_idx), e_sp)
                            q_cons_vf(eqn_idx%E)%sf(j, k, l) = dyn_pres + rho*e_sp
                        else if (mhd) then
                            if (n == 0) then
                                pres_mag = 0.5_wp*(Bx0**2 + q_prim_vf(eqn_idx%B%beg)%sf(j, k, &
                                                   & l)**2 + q_prim_vf(eqn_idx%B%beg + 1)%sf(j, k, l)**2)
                            else
                                pres_mag = 0.5_wp*(q_prim_vf(eqn_idx%B%beg)%sf(j, k, l)**2 + q_prim_vf(eqn_idx%B%beg + 1)%sf(j, &
                                                   & k, l)**2 + q_prim_vf(eqn_idx%B%beg + 2)%sf(j, k, l)**2)
                            end if
                            ! MHD energy includes magnetic pressure contribution
                            q_cons_vf(eqn_idx%E)%sf(j, k, l) = gamma*q_prim_vf(eqn_idx%E)%sf(j, k, &
                                      & l) + dyn_pres + pres_mag + pi_inf + qv
                        else if ((model_eqns /= 4) .and. (bubbles_euler .neqv. .true.)) then
                            ! Five-equation model (Allaire et al. JCP 2002): E = Gamma*p + 0.5*rho*|u|^2 + pi_inf + qv
                            q_cons_vf(eqn_idx%E)%sf(j, k, l) = gamma*q_prim_vf(eqn_idx%E)%sf(j, k, l) + dyn_pres + pi_inf + qv
                        else if ((model_eqns /= 4) .and. (bubbles_euler)) then
                            ! Bubble-augmented energy with void fraction correction
                            q_cons_vf(eqn_idx%E)%sf(j, k, l) = dyn_pres + (1._wp - q_prim_vf(eqn_idx%alf)%sf(j, k, &
                                      & l))*(gamma*q_prim_vf(eqn_idx%E)%sf(j, k, l) + pi_inf)
                        else
                            ! Four-equation model (Kapila et al. PoF 2001): Tait EOS, no conserved energy variable
                            q_cons_vf(eqn_idx%E)%sf(j, k, l) = 0._wp
                        end if
                    end if

                    ! Six-equation model (Saurel et al. JCP 2009): compute per-phase internal energies
                    if (model_eqns == 3) then
                        phase_energy_sum = 0._wp
                        do i = 1, num_fluids
                            if (eos_idxs(i) == 2) then
                                rho_phase = q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, &
                                                      & l)/max(q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l), sgm_eps)
                                call s_jwl_energy_pr(rho_phase, q_prim_vf(eqn_idx%E)%sf(j, k, l), 1._wp, 1._wp, jwl_As(i), &
                                                     & jwl_Bs(i), jwl_R1s(i), jwl_R2s(i), jwl_omegas(i), jwl_rho0s(i), &
                                                     & jwl_E0s(i), jwl_air_e0s(i), jwl_air_rho0s(i), jwl_air_gammas(i), e_sp)
                                q_cons_vf(i + eqn_idx%int_en%beg - 1)%sf(j, k, l) = q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, &
                                          & l)*e_sp
                            else
                                q_cons_vf(i + eqn_idx%int_en%beg - 1)%sf(j, k, l) = q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, &
                                          & l)*(gammas(i)*q_prim_vf(eqn_idx%E)%sf(j, k, &
                                          & l) + pi_infs(i)) + q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l)*qvs(i)
                            end if
                            phase_energy_sum = phase_energy_sum + q_cons_vf(i + eqn_idx%int_en%beg - 1)%sf(j, k, l)
                        end do
                        q_cons_vf(eqn_idx%E)%sf(j, k, l) = dyn_pres + phase_energy_sum
                    end if

                    if (bubbles_euler) then
                        ! From prim: Compute nbub = (3/4pi) * \alpha / \bar{R^3}
                        do i = 1, nb
                            Rtmp(i) = q_prim_vf(qbmm_idx%rs(i))%sf(j, k, l)
                        end do

                        if (.not. qbmm) then
                            if (adv_n) then
                                q_cons_vf(eqn_idx%n)%sf(j, k, l) = q_prim_vf(eqn_idx%n)%sf(j, k, l)
                                nbub = q_prim_vf(eqn_idx%n)%sf(j, k, l)
                            else
                                call s_comp_n_from_prim(real(q_prim_vf(eqn_idx%alf)%sf(j, k, l), kind=wp), Rtmp, nbub, weight)
                            end if
                        else
                            ! Initialize R3 averaging over R0 and R directions
                            R3tmp = 0._wp
                            do i = 1, nb
                                R3tmp = R3tmp + weight(i)*0.5_wp*(Rtmp(i) + sigR)**3._wp
                                R3tmp = R3tmp + weight(i)*0.5_wp*(Rtmp(i) - sigR)**3._wp
                            end do
                            ! Initialize nb
                            nbub = 3._wp*q_prim_vf(eqn_idx%alf)%sf(j, k, l)/(4._wp*pi*R3tmp)
                        end if

                        do i = eqn_idx%bub%beg, eqn_idx%bub%end
                            q_cons_vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)*nbub
                        end do
                    end if

                    if (mhd) then
                        do i = eqn_idx%B%beg, eqn_idx%B%end
                            q_cons_vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)
                        end do
                    end if

                    if (elasticity) then
                        ! adding the elastic contribution Multiply \tau to \rho \tau
                        do i = eqn_idx%stress%beg, eqn_idx%stress%end
                            q_cons_vf(i)%sf(j, k, l) = rho*q_prim_vf(i)%sf(j, k, l)
                        end do
                    end if

                    if (hypoelasticity) then
                        if (cont_damage) G = G*max((1._wp - q_prim_vf(eqn_idx%damage)%sf(j, k, l)), 0._wp)
                        do i = eqn_idx%stress%beg, eqn_idx%stress%end
                            ! adding elastic contribution
                            if (G > verysmall) then
                                q_cons_vf(eqn_idx%E)%sf(j, k, l) = q_cons_vf(eqn_idx%E)%sf(j, k, l) + (q_prim_vf(i)%sf(j, k, &
                                          & l)**2._wp)/(4._wp*G)
                                ! Double for shear stresses
                                if (any(i == shear_indices)) then
                                    q_cons_vf(eqn_idx%E)%sf(j, k, l) = q_cons_vf(eqn_idx%E)%sf(j, k, l) + (q_prim_vf(i)%sf(j, k, &
                                              & l)**2._wp)/(4._wp*G)
                                end if
                            end if
                        end do
                    end if

                    ! using \rho xi as the conservative formulation stated in Kamrin et al. JFM 2022
                    if (hyperelasticity) then
                        ! Multiply \xi to \rho \xi
                        do i = eqn_idx%xi%beg, eqn_idx%xi%end
                            q_cons_vf(i)%sf(j, k, l) = rho*q_prim_vf(i)%sf(j, k, l)
                        end do
                    end if

                    if (surface_tension) then
                        q_cons_vf(eqn_idx%c)%sf(j, k, l) = q_prim_vf(eqn_idx%c)%sf(j, k, l)
                    end if

                    if (cont_damage) q_cons_vf(eqn_idx%damage)%sf(j, k, l) = q_prim_vf(eqn_idx%damage)%sf(j, k, l)

                    if (hyper_cleaning) q_cons_vf(eqn_idx%psi)%sf(j, k, l) = q_prim_vf(eqn_idx%psi)%sf(j, k, l)
                end do
            end do
        end do
#else
        if (proc_rank == 0) then
            call s_mpi_abort('Conversion from primitive to ' // 'conservative variables not ' // 'implemented. Exiting.')
        end if
#endif

    end subroutine s_convert_primitive_to_conservative_variables

    !> Convert primitive variables to Eulerian flux variables.
    subroutine s_convert_primitive_to_flux_variables(qK_prim_vf, FK_vf, FK_src_vf, is1, is2, is3, s2b, s3b)

        integer, intent(in)                                                                     :: s2b, s3b
        real(wp), dimension(0:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(in)                  :: qK_prim_vf
        real(wp), dimension(0:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout)               :: FK_vf
        real(wp), dimension(0:,idwbuff(2)%beg:,idwbuff(3)%beg:,eqn_idx%adv%beg:), intent(inout) :: FK_src_vf
        type(int_bounds_info), intent(in)                                                       :: is1, is2, is3

        ! Partial densities, density, velocity, pressure, energy, advection variables, the specific heat ratio and liquid stiffness
        ! functions, the shear and volume Reynolds numbers and the Weber numbers

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3)  :: alpha_rho_K
            real(wp), dimension(3)  :: alpha_K
            real(wp), dimension(3)  :: vel_K
            real(wp), dimension(10) :: Y_K
        #:else
            real(wp), dimension(num_fluids)  :: alpha_rho_K
            real(wp), dimension(num_fluids)  :: alpha_K
            real(wp), dimension(num_vels)    :: vel_K
            real(wp), dimension(num_species) :: Y_K
        #:endif
        real(wp)               :: rho_K
        real(wp)               :: vel_K_sum
        real(wp)               :: pres_K
        real(wp)               :: E_K
        real(wp)               :: gamma_K
        real(wp)               :: pi_inf_K
        real(wp)               :: qv_K
        real(wp), dimension(2) :: Re_K
        real(wp)               :: G_K
        real(wp)               :: T_K, mix_mol_weight, R_gas
        real(wp)               :: pcold, e_sp, Y_jwl, alpha_jwl
        integer                :: i, j, k, l  !< Generic loop iterators

        is1b = is1%beg; is1e = is1%end
        is2b = is2%beg; is2e = is2%end
        is3b = is3%beg; is3e = is3%end

        $:GPU_UPDATE(device='[is1b, is2b, is3b, is1e, is2e, is3e]')

        ! Computing the flux variables from the primitive variables, without accounting for the contribution of either viscosity or
        ! capillarity
#ifdef MFC_SIMULATION
        $:GPU_PARALLEL_LOOP(collapse=3, private='[alpha_rho_K, vel_K, alpha_K, Re_K, Y_K, rho_K, vel_K_sum, pres_K, E_K, gamma_K, &
                            & pi_inf_K, qv_K, G_K, T_K, mix_mol_weight, R_gas, pcold, e_sp, Y_jwl, alpha_jwl]')
        do l = is3b, is3e
            do k = is2b, is2e
                do j = is1b, is1e
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, eqn_idx%cont%end
                        alpha_rho_K(i) = qK_prim_vf(j, k, l, i)
                    end do

                    $:GPU_LOOP(parallelism='[seq]')
                    do i = eqn_idx%adv%beg, eqn_idx%adv%end
                        alpha_K(i - eqn_idx%E) = qK_prim_vf(j, k, l, i)
                    end do

                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_vels
                        vel_K(i) = qK_prim_vf(j, k, l, eqn_idx%cont%end + i)
                    end do

                    vel_K_sum = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_vels
                        vel_K_sum = vel_K_sum + vel_K(i)**2._wp
                    end do

                    pres_K = qK_prim_vf(j, k, l, eqn_idx%E)
                    if (elasticity) then
                        call s_convert_species_to_mixture_variables_acc(rho_K, gamma_K, pi_inf_K, qv_K, alpha_K, alpha_rho_K, &
                            & Re_K, G_K, Gs_vc)
                    else
                        call s_convert_species_to_mixture_variables_acc(rho_K, gamma_K, pi_inf_K, qv_K, alpha_K, alpha_rho_K, Re_K)
                    end if

                    ! Computing the energy from the pressure

                    if (chemistry) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%species%beg, eqn_idx%species%end
                            Y_K(i - eqn_idx%species%beg + 1) = qK_prim_vf(j, k, l, i)
                        end do
                        ! Computing the energy from the internal energy of the mixture
                        call get_mixture_molecular_weight(Y_k, mix_mol_weight)
                        R_gas = gas_constant/mix_mol_weight
                        T_K = pres_K/rho_K/R_gas
                        call get_mixture_energy_mass(T_K, Y_K, E_K)
                        E_K = rho_K*E_K + 5.e-1_wp*rho_K*vel_K_sum
                    else
                        ! Computing the energy from the pressure
                        if ((jwl_idx > 0) .and. (.not. mhd) .and. (model_eqns /= 4) .and. (bubbles_euler .neqv. .true.)) then
                            Y_jwl = 1._wp
                            alpha_jwl = 1._wp
                            if (jwl_idx <= eqn_idx%cont%end) then
                                Y_jwl = alpha_rho_K(jwl_idx)/max(rho_K, sgm_eps)
                                alpha_jwl = alpha_K(jwl_idx)
                            end if
                            call s_jwl_energy_pr(rho_K, pres_K, Y_jwl, alpha_jwl, jwl_As(jwl_idx), jwl_Bs(jwl_idx), &
                                                 & jwl_R1s(jwl_idx), jwl_R2s(jwl_idx), jwl_omegas(jwl_idx), jwl_rho0s(jwl_idx), &
                                                 & jwl_E0s(jwl_idx), jwl_air_e0s(jwl_idx), jwl_air_rho0s(jwl_idx), &
                                                 & jwl_air_gammas(jwl_idx), e_sp)
                            E_K = rho_K*e_sp + 5.e-1_wp*rho_K*vel_K_sum
                        else
                            E_K = gamma_K*pres_K + pi_inf_K + 5.e-1_wp*rho_K*vel_K_sum + qv_K
                        end if
                    end if

                    ! mass flux, this should be \alpha_i \rho_i u_i
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, eqn_idx%cont%end
                        FK_vf(j, k, l, i) = alpha_rho_K(i)*vel_K(dir_idx(1))
                    end do

                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_vels
                        FK_vf(j, k, l, &
                              & eqn_idx%cont%end + dir_idx(i)) = rho_K*vel_K(dir_idx(1))*vel_K(dir_idx(i)) &
                              & + pres_K*dir_flg(dir_idx(i))
                    end do

                    ! energy flux, u(E+p)
                    FK_vf(j, k, l, eqn_idx%E) = vel_K(dir_idx(1))*(E_K + pres_K)

                    ! Species advection Flux, \rho*u*Y
                    if (chemistry) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, num_species
                            FK_vf(j, k, l, i - 1 + eqn_idx%species%beg) = vel_K(dir_idx(1))*(rho_K*Y_K(i))
                        end do
                    end if

                    if (riemann_solver == 1 .or. riemann_solver == 4) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%adv%beg, eqn_idx%adv%end
                            FK_vf(j, k, l, i) = 0._wp
                            FK_src_vf(j, k, l, i) = alpha_K(i - eqn_idx%E)
                        end do
                    else
                        ! Could be bubbles_euler!
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%adv%beg, eqn_idx%adv%end
                            FK_vf(j, k, l, i) = vel_K(dir_idx(1))*alpha_K(i - eqn_idx%E)
                        end do

                        $:GPU_LOOP(parallelism='[seq]')
                        do i = eqn_idx%adv%beg, eqn_idx%adv%end
                            FK_src_vf(j, k, l, i) = vel_K(dir_idx(1))
                        end do
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
#endif

    end subroutine s_convert_primitive_to_flux_variables

    !> Compute partial densities and volume fractions
    subroutine s_compute_species_fraction(q_vf, k, l, r, alpha_rho_K, alpha_K)

        $:GPU_ROUTINE(function_name='s_compute_species_fraction', parallelism='[seq]', cray_noinline=True)
        type(scalar_field), dimension(sys_size), intent(in) :: q_vf
        integer, intent(in)                                 :: k, l, r
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(out) :: alpha_rho_K, alpha_K
        #:else
            real(wp), dimension(num_fluids), intent(out) :: alpha_rho_K, alpha_K
        #:endif
        integer  :: i
        real(wp) :: alpha_K_sum

        if (num_fluids == 1) then
            alpha_rho_K(1) = q_vf(eqn_idx%cont%beg)%sf(k, l, r)
            if (igr .or. bubbles_euler) then
                alpha_K(1) = 1._wp
            else
                alpha_K(1) = q_vf(eqn_idx%adv%beg)%sf(k, l, r)
            end if
        else
            if (igr) then
                do i = 1, num_fluids - 1
                    alpha_rho_K(i) = q_vf(i)%sf(k, l, r)
                    alpha_K(i) = q_vf(eqn_idx%adv%beg + i - 1)%sf(k, l, r)
                end do
                alpha_rho_K(num_fluids) = q_vf(num_fluids)%sf(k, l, r)
                alpha_K(num_fluids) = 1._wp - sum(alpha_K(1:num_fluids - 1))
            else
                do i = 1, num_fluids
                    alpha_rho_K(i) = q_vf(i)%sf(k, l, r)
                    alpha_K(i) = q_vf(eqn_idx%adv%beg + i - 1)%sf(k, l, r)
                end do
            end if
        end if

        if (mpp_lim) then
            alpha_K_sum = 0._wp
            do i = 1, num_fluids
                alpha_rho_K(i) = max(0._wp, alpha_rho_K(i))
                alpha_K(i) = min(max(0._wp, alpha_K(i)), 1._wp)
                alpha_K_sum = alpha_K_sum + alpha_K(i)
            end do
            alpha_K = alpha_K/max(alpha_K_sum, 1.e-16_wp)
        end if

        if (num_fluids == 1 .and. bubbles_euler) alpha_K(1) = q_vf(eqn_idx%adv%beg)%sf(k, l, r)

    end subroutine s_compute_species_fraction

    !> Deallocate fluid property arrays and post-processing fields allocated during module initialization.
    impure subroutine s_finalize_variables_conversion_module()

#ifdef MFC_SIMULATION
        integer :: jwl_pressure_floor_count_glb
#endif

        ! Deallocating the density, the specific heat ratio function and the liquid stiffness function
#ifdef MFC_POST_PROCESS
        deallocate (rho_sf, gamma_sf, pi_inf_sf, qv_sf)
#endif

#ifdef MFC_SIMULATION
        $:GPU_UPDATE(host='[jwl_pressure_floor_count]')
        call s_mpi_allreduce_integer_sum(jwl_pressure_floor_count, jwl_pressure_floor_count_glb)
        if (proc_rank == 0 .and. jwl_pressure_floor_count_glb > 0) then
            print '(A,I0,A)', 'WARNING: JWL pressure recovery applied the 1 Pa floor ', jwl_pressure_floor_count_glb, ' times.'
        end if

        @:DEALLOCATE(eos_idxs, gammas, gs_min, pi_infs, ps_inf, cvs, qvs, qvps)
        @:DEALLOCATE(jwl_As, jwl_Bs, jwl_R1s, jwl_R2s, jwl_omegas, jwl_rho0s, jwl_E0s)
        @:DEALLOCATE(jwl_air_e0s, jwl_air_rho0s, jwl_air_gammas, Gs_vc)
        if (bubbles_euler) then
            @:DEALLOCATE(bubrs_vc)
        end if
#else
        @:DEALLOCATE(eos_idxs, gammas, gs_min, pi_infs, ps_inf, cvs, qvs, qvps)
        @:DEALLOCATE(jwl_As, jwl_Bs, jwl_R1s, jwl_R2s, jwl_omegas, jwl_rho0s, jwl_E0s)
        @:DEALLOCATE(jwl_air_e0s, jwl_air_rho0s, jwl_air_gammas, Gs_vc)
        if (bubbles_euler) then
            @:DEALLOCATE(bubrs_vc)
        end if
#endif

    end subroutine s_finalize_variables_conversion_module

#ifndef MFC_PRE_PROCESS
    !> Compute the speed of sound from thermodynamic state variables, supporting multiple equation-of-state models.
    subroutine s_compute_speed_of_sound(pres, rho, gamma, pi_inf, H, adv, vel_sum, c_c, c, qv, jwl_Y, jwl_alpha)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in) :: pres
        real(wp), intent(in) :: rho, gamma, pi_inf, qv
        real(wp), intent(in) :: H
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(in) :: adv
        #:else
            real(wp), dimension(num_fluids), intent(in) :: adv
        #:endif
        real(wp), intent(in)           :: vel_sum
        real(wp), intent(in)           :: c_c
        real(wp), intent(out)          :: c
        real(wp), intent(in), optional :: jwl_Y, jwl_alpha
        real(wp)                       :: blkmod1, blkmod2
        real(wp)                       :: Y_jwl, alpha_jwl
        integer                        :: q

        if (chemistry) then  ! Reacting mixture sound speed
            if (avg_state == 1 .and. abs(c_c) > verysmall) then
                c = sqrt(c_c - (gamma - 1.0_wp)*(vel_sum - H))
            else
                c = sqrt((1.0_wp + 1.0_wp/gamma)*pres/rho)
            end if
        else if (relativity) then  ! Relativistic sound speed
            c = sqrt((1._wp + 1._wp/gamma)*pres/rho/H)
        else
            if ((jwl_idx > 0) .and. (.not. mhd) .and. (model_eqns /= 4) .and. (bubbles_euler .neqv. .true.)) then
                Y_jwl = 1._wp
                alpha_jwl = 1._wp
                if (present(jwl_Y)) then
                    Y_jwl = jwl_Y
                else if (jwl_idx <= size(adv)) then
                    Y_jwl = adv(jwl_idx)
                end if
                if (present(jwl_alpha)) then
                    alpha_jwl = jwl_alpha
                else if (jwl_idx <= size(adv)) then
                    alpha_jwl = adv(jwl_idx)
                end if
                call s_jwl_mixture_sound_speed_squared(rho, pres, Y_jwl, alpha_jwl, jwl_As(jwl_idx), jwl_Bs(jwl_idx), &
                                                       & jwl_R1s(jwl_idx), jwl_R2s(jwl_idx), jwl_omegas(jwl_idx), &
                                                       & jwl_rho0s(jwl_idx), jwl_E0s(jwl_idx), jwl_air_e0s(jwl_idx), &
                                                       & jwl_air_rho0s(jwl_idx), jwl_air_gammas(jwl_idx), c)
                c = max(c, jwl_omegas(jwl_idx)*max(pres, 1._wp)/(max(rho, sgm_eps)*jwl_rho0s(jwl_idx)))
            else if (alt_soundspeed) then  ! Wood's mixture sound speed via bulk moduli
                blkmod1 = ((gammas(1) + 1._wp)*pres + pi_infs(1))/gammas(1)
                blkmod2 = ((gammas(2) + 1._wp)*pres + pi_infs(2))/gammas(2)
                c = (1._wp/(rho*(adv(1)/blkmod1 + adv(2)/blkmod2)))
            else if (model_eqns == 3) then  ! Six-equation model sound speed
                c = 0._wp
                $:GPU_LOOP(parallelism='[seq]')
                do q = 1, num_fluids
                    c = c + adv(q)*gs_min(q)*(pres + pi_infs(q)/(gammas(q) + 1._wp))
                end do
                c = c/rho
            else if (((model_eqns == 4) .or. (model_eqns == 2 .and. bubbles_euler))) then
                ! Sound speed for bubble mixture to order O(\alpha)

                if (mpp_lim .and. (num_fluids > 1)) then
                    c = (1._wp/gamma + 1._wp)*(pres + pi_inf/(gamma + 1._wp))/rho
                else
                    c = (1._wp/gamma + 1._wp)*(pres + pi_inf/(gamma + 1._wp))/(rho*(1._wp - adv(num_fluids)))
                end if
            else
                c = (H - 5.e-1*vel_sum - qv/rho)/gamma
            end if

            if (mixture_err .and. c < 0._wp) then
                c = 100._wp*sgm_eps
            else
                c = sqrt(c)
            end if
        end if

    end subroutine s_compute_speed_of_sound
#endif

#ifndef MFC_PRE_PROCESS
    !> Compute the fast magnetosonic wave speed from the sound speed, density, and magnetic field components.
    subroutine s_compute_fast_magnetosonic_speed(rho, c, B, norm, c_fast, h)

        $:GPU_ROUTINE(function_name='s_compute_fast_magnetosonic_speed', parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: B(3), rho, c
        real(wp), intent(in)  :: h  !< only used for relativity
        real(wp), intent(out) :: c_fast
        integer, intent(in)   :: norm
        real(wp)              :: B2, term, disc

        B2 = sum(B**2)

        if (.not. relativity) then
            term = c**2 + B2/rho
            disc = term**2 - 4*c**2*(B(norm)**2/rho)
        else
            ! Note: this is approximation for the non-relatisitic limit; accurate solution requires solving a quartic equation
            term = (c**2*(B(norm)**2 + rho*h) + B2)/(rho*h + B2)
            disc = term**2 - 4*c**2*B(norm)**2/(rho*h + B2)
        end if

#ifdef MFC_DEBUG
        if (disc < 0._wp) then
            print *, 'rho, c, Bx, By, Bz, h, term, disc:', rho, c, B(1), B(2), B(3), h, term, disc
            call s_mpi_abort('Error: negative discriminant in s_compute_fast_magnetosonic_speed')
        end if
#endif

        c_fast = sqrt(0.5_wp*(term + sqrt(disc)))

    end subroutine s_compute_fast_magnetosonic_speed
#endif
end module m_variables_conversion
