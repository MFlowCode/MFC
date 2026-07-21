!>
!! @file m_thermodynamics.fpp
!! @brief Contains module m_thermodynamics

#:include 'macros.fpp'
#:include 'case.fpp'

!> @brief Central equation-of-state interface. Every genuine thermodynamic query (pressure here; energy, temperature, sound speed,
!! and derivatives in later PRs) is answered through one module so no solver path carries its own EOS algebra. Compile-time
!! chemistry selects the adapter: the stiffened-gas analytic branch or the Pyrometheus temperature-then-query branch. This is a leaf
!! module (it uses global parameters and the thermochemistry surface, never m_variables_conversion) so no cycle forms.
module m_thermodynamics

    use m_derived_types
    use m_global_parameters
    use m_constants, only: model_eqns_4eq, model_eqns_5eq, model_eqns_6eq, avg_state_roe, verysmall, sgm_eps
    use m_thermochem, only: num_species, get_temperature, get_pressure, gas_constant, get_mixture_molecular_weight, &
        & get_mixture_energy_mass

    implicit none

    private
    public :: s_compute_pressure
    public :: s_compute_internal_energy
#ifndef MFC_PRE_PROCESS
    public :: s_compute_speed_of_sound
#endif

contains

    !> Pressure from the conserved energy and mixture properties. The #:if chemistry split selects the stiffened-gas adapter
    !! (analytic inverse of the mixture EOS) or the Pyrometheus adapter (solve temperature, then query pressure). Arithmetic is
    !! relocated unchanged from the conversion module.
    subroutine s_compute_pressure(energy, alf, dyn_p, pi_inf, gamma, rho, qv, rhoYks, pres, T, stress, mom, G, pres_mag)

        $:GPU_ROUTINE(function_name='s_compute_pressure',parallelism='[seq]', cray_noinline=True)

        real(stp), intent(in)           :: energy, alf
        real(wp), intent(in)            :: dyn_p
        real(wp), intent(in)            :: pi_inf, gamma, rho, qv
        real(wp), intent(out)           :: pres
        real(wp), intent(inout)         :: T
        real(stp), intent(in), optional :: stress, mom
        real(wp), intent(in), optional  :: G, pres_mag

        ! Chemistry
        real(wp), dimension(1:num_species), intent(in) :: rhoYks
        real(wp), dimension(1:num_species)             :: Y_rs
        real(wp)                                       :: E_e
        real(wp)                                       :: e_Per_Kg, Pdyn_Per_Kg
        real(wp)                                       :: T_guess
        integer                                        :: s  !< Generic loop iterator
        #:if not chemistry
            ! Depending on model_eqns and bubbles_euler, the appropriate procedure for computing pressure is targeted by the
            ! procedure pointer

            if (mhd) then
                ! MHD pressure: subtract magnetic pressure from total energy
                pres = (energy - dyn_p - pi_inf - qv - pres_mag)/gamma
            else if ((model_eqns /= model_eqns_4eq) .and. (bubbles_euler .neqv. .true.)) then
                ! Gamma/pi_inf model or five-equation model (Allaire et al. JCP 2002): p from mixture EOS
                pres = (energy - dyn_p - pi_inf - qv)/gamma
            else if ((model_eqns /= model_eqns_4eq) .and. bubbles_euler) then
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

                pres = (energy - 0.5_wp*(mom**2._wp)/rho - pi_inf - qv - E_e)/gamma
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

    !> Conserved energy from the primitive pressure and mixture properties: the inverse of s_compute_pressure over the same model
    !! branches. The chemistry adapter builds the mixture internal energy from temperature and species. Arithmetic is relocated
    !! unchanged from the conversion module.
    subroutine s_compute_internal_energy(pres, alf, dyn_p, pi_inf, gamma, rho, qv, Ys, energy, pres_mag)

        $:GPU_ROUTINE(function_name='s_compute_internal_energy',parallelism='[seq]', cray_noinline=True)

        real(stp), intent(in)                          :: pres, alf
        real(wp), intent(in)                           :: dyn_p, pi_inf, gamma, rho, qv
        real(stp), intent(out)                         :: energy
        real(wp), intent(in), optional                 :: pres_mag
        real(wp), dimension(1:num_species), intent(in) :: Ys
        real(wp)                                       :: e_mix, mix_mol_weight, T
        #:if not chemistry
            ! Computing the energy from the pressure
            if (mhd) then
                ! MHD energy includes magnetic pressure contribution
                energy = gamma*pres + dyn_p + pres_mag + pi_inf + qv
            else if ((model_eqns /= model_eqns_4eq) .and. (bubbles_euler .neqv. .true.)) then
                ! Five-equation model (Allaire et al. JCP 2002): E = Gamma*p + 0.5*rho*|u|^2 + pi_inf + qv
                energy = gamma*pres + dyn_p + pi_inf + qv
            else if ((model_eqns /= model_eqns_4eq) .and. bubbles_euler) then
                ! Bubble-augmented energy with void fraction correction
                energy = dyn_p + (1._wp - alf)*(gamma*pres + pi_inf)
            else
                ! Four-equation model (Kapila et al. PoF 2001): Tait EOS, no conserved energy variable
                energy = 0._wp
            end if
        #:else
            ! Reacting mixture: compute conserved energy from species mass fractions and temperature
            call get_mixture_molecular_weight(Ys, mix_mol_weight)
            T = pres*mix_mol_weight/(gas_constant*rho)
            call get_mixture_energy_mass(T, Ys, e_mix)
            energy = dyn_p + rho*e_mix
        #:endif

    end subroutine s_compute_internal_energy

#ifndef MFC_PRE_PROCESS
    !> Speed of sound from thermodynamic state, supporting the several equation-of-state models. Arithmetic is relocated unchanged
    !! from the conversion module: the non-chemistry branches hold c squared until a single final sqrt, the chemistry branch returns
    !! c already square-rooted, and the mixture_err floor is preserved.
    subroutine s_compute_speed_of_sound(pres, rho, gamma, pi_inf, H, adv, vel_sum, c_c, c, qv)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in) :: pres
        real(wp), intent(in) :: rho, gamma, pi_inf, qv
        real(wp), intent(in) :: H
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(in) :: adv
        #:else
            real(wp), dimension(num_fluids), intent(in) :: adv
        #:endif
        real(wp), intent(in)  :: vel_sum
        real(wp), intent(in)  :: c_c
        real(wp), intent(out) :: c
        real(wp)              :: blkmod1, blkmod2
        integer               :: q

        if (chemistry) then  ! Reacting mixture sound speed
            if (avg_state == avg_state_roe .and. abs(c_c) > verysmall) then
                c = sqrt(c_c - (gamma - 1.0_wp)*(vel_sum - H))
            else
                c = sqrt((1.0_wp + 1.0_wp/gamma)*pres/rho)
            end if
        else if (relativity) then  ! Relativistic sound speed
            c = sqrt((1._wp + 1._wp/gamma)*pres/rho/H)
        else
            if (alt_soundspeed) then  ! Wood's mixture sound speed via bulk moduli
                blkmod1 = ((gammas(1) + 1._wp)*pres + pi_infs(1))/gammas(1)
                blkmod2 = ((gammas(2) + 1._wp)*pres + pi_infs(2))/gammas(2)
                c = (1._wp/(rho*(adv(1)/blkmod1 + adv(2)/blkmod2)))
            else if (model_eqns == model_eqns_6eq) then  ! Six-equation model sound speed
                c = 0._wp
                $:GPU_LOOP(parallelism='[seq]')
                do q = 1, num_fluids
                    c = c + adv(q)*gs_min(q)*(pres + pi_infs(q)/(gammas(q) + 1._wp))
                end do
                c = c/rho
            else if (((model_eqns == model_eqns_4eq) .or. (model_eqns == model_eqns_5eq .and. bubbles_euler))) then
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
end module m_thermodynamics
