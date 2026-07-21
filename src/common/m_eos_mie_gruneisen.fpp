!>
!! @file m_eos_mie_gruneisen.fpp
!! @brief Contains module m_eos_mie_gruneisen

#:include 'macros.fpp'

!> @brief Generic Mie-Gruneisen equation of state with a constant Gruneisen coefficient, structured as a reference-curve provider
!! plus a curve-agnostic assembler. A Mie-Gruneisen material obeys p = p_ref(rho) + Gamma*rho*(e - e_ref(rho)); the reference pair
!! (p_ref, e_ref) and its density derivatives are the only thing that distinguishes one Mie-Gruneisen material from another.
!! Following the production hydrocode structure (Arienti, Morano and Shepherd, GALCIT FM99-8, 2004; LLNL singularity-eos), one
!! provider returns the reference curve for a chosen material and one assembler consumes it, so stiffened gas, ideal gas, and JWL
!! differ only in the provider. This PR ships the stiffened-gas reference curve, which recovers the existing stiffened-gas pressure
!! and sound speed exactly (Gamma = gamma_s - 1, p_ref = -pi_inf, e_ref = pi_inf/rho); the JWL principal isentrope is the next
!! provider, added with its consumer. This is a leaf module (it uses only the precision kinds) and no solver path calls it yet: it
!! is the verified thermodynamic core for the Mie-Gruneisen backend, landed ahead of the reference curve that consumes it.
module m_eos_mie_gruneisen

    use m_precision_select

    implicit none

    private
    public :: s_mg_stiffened_reference, f_mg_pressure, f_mg_internal_energy, f_mg_sound_speed_sq

contains

    !> Stiffened gas expressed as a Mie-Gruneisen reference curve. Returns the Gruneisen coefficient, the reference pressure and
    !! energy, and their density derivatives, which together make the generic assembler reproduce stiffened gas. This is the first
    !! concrete reference curve and the template the JWL principal isentrope follows. The reference pressure is density-independent
    !! for stiffened gas, so its derivative is zero; ideal gas is the pi_inf = 0 sub-case.
    pure subroutine s_mg_stiffened_reference(gamma_s, pi_inf, rho, gamma_mg, p_ref, dp_ref, e_ref, de_ref)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)  :: gamma_s, pi_inf, rho
        real(wp), intent(out) :: gamma_mg, p_ref, dp_ref, e_ref, de_ref

        gamma_mg = gamma_s - 1._wp
        p_ref = -pi_inf
        dp_ref = 0._wp
        e_ref = pi_inf/rho
        de_ref = -pi_inf/rho**2

    end subroutine s_mg_stiffened_reference

    !> Mie-Gruneisen pressure from density and specific internal energy against a reference curve: p = p_ref + Gamma*rho*(e -
    !! e_ref).
    pure function f_mg_pressure(gamma_mg, rho, e_int, p_ref, e_ref) result(pres)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in) :: gamma_mg, rho, e_int, p_ref, e_ref
        real(wp)             :: pres

        pres = p_ref + gamma_mg*rho*(e_int - e_ref)

    end function f_mg_pressure

    !> Specific internal energy from pressure, the inverse of f_mg_pressure: e = e_ref + (p - p_ref)/(Gamma*rho). The Gamma*rho
    !! divisor is positive for every admissible state (Gamma > 0, rho > 0).
    pure function f_mg_internal_energy(gamma_mg, rho, pres, p_ref, e_ref) result(e_int)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in) :: gamma_mg, rho, pres, p_ref, e_ref
        real(wp)             :: e_int

        e_int = e_ref + (pres - p_ref)/(gamma_mg*rho)

    end function f_mg_internal_energy

    !> Squared sound speed c^2 = dp/drho at constant entropy for a constant Gruneisen coefficient, from the frozen identity c^2 =
    !! dp/drho|_e + (p/rho^2) dp/de|_rho with dp/de|_rho = Gamma*rho (Arienti et al. 2004, Eq. 12; Menikoff and Plohr 1989). Written
    !! in terms of pressure so no separate internal-energy argument is needed: c^2 = dp_ref - Gamma*rho*de_ref + ((1 + Gamma)*p -
    !! p_ref)/rho, where dp_ref and de_ref are the reference-curve density derivatives dp_ref/drho and de_ref/drho. The result can
    !! be negative outside the admissible domain, so the caller tests it before taking a square root.
    pure function f_mg_sound_speed_sq(gamma_mg, rho, pres, p_ref, dp_ref, de_ref) result(c_sq)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in) :: gamma_mg, rho, pres, p_ref, dp_ref, de_ref
        real(wp)             :: c_sq

        c_sq = dp_ref - gamma_mg*rho*de_ref + ((1._wp + gamma_mg)*pres - p_ref)/rho

    end function f_mg_sound_speed_sq

end module m_eos_mie_gruneisen
