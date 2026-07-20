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
    use m_constants, only: model_eqns_4eq
    use m_thermochem, only: num_species, get_temperature, get_pressure

    implicit none

    private
    public :: s_compute_pressure

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

end module m_thermodynamics
