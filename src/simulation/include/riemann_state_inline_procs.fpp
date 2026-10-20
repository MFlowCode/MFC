!>
!! @file riemann_state_inline_procs.fpp
!! @brief Per-face Riemann helpers, one source for m_riemann_state and the Riemann solver modules that call them from
!! GPU kernels, kept in the Fypp macro RS_INLINE_PROCS. Same reason and usage as common/include/inline_eos_procs.fpp;
!! the rename / private helpers are RS_INLINE_USE / RS_INLINE_PRIVATE in riemann_state_inline.fpp. Edit the routines
!! here only.

#:def RS_INLINE_PROCS()
    !> Low-Mach parameter of Thornber et al. JCP (2008): the larger of the two face Mach numbers, capped at one so the correction
    !! switches itself off once the flow is no longer low speed.
    function f_low_Mach_zcoef(vel_L_rms, vel_R_rms, c_L, c_R) result(zcoef)

        $:GPU_ROUTINE(function_name='f_low_Mach_zcoef', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in) :: vel_L_rms, vel_R_rms  !< Left and right squared velocity magnitudes
        real(wp), intent(in) :: c_L, c_R              !< Left and right sound speeds
        real(wp)             :: zcoef

        zcoef = min(1._wp, max(vel_L_rms**5.e-1_wp/c_L, vel_R_rms**5.e-1_wp/c_R))

    end function f_low_Mach_zcoef

    !> The same correction for the HLLC flux, where the star state supplies the pressure jump directly and the correction scales
    !! with the mass flux through the acoustic waves instead. Zero unless low_Mach == 1.
    function f_low_Mach_pcorr_hllc(vel_L_rms, vel_R_rms, c_L, c_R, rho_L, rho_R, s_L, s_R, vel_L_norm, vel_R_norm) result(pcorr)

        $:GPU_ROUTINE(function_name='f_low_Mach_pcorr_hllc', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in) :: vel_L_rms, vel_R_rms    !< Left and right squared velocity magnitudes
        real(wp), intent(in) :: c_L, c_R                !< Left and right sound speeds
        real(wp), intent(in) :: rho_L, rho_R            !< Left and right densities
        real(wp), intent(in) :: s_L, s_R                !< Left and right wave speeds
        real(wp), intent(in) :: vel_L_norm, vel_R_norm  !< Left and right wave-normal velocities
        real(wp)             :: pcorr

        pcorr = 0._wp
        if (low_Mach == 1) then
            pcorr = rho_L*rho_R*(s_L - vel_L_norm)*(s_R - vel_R_norm)*(vel_R_norm - vel_L_norm)/(rho_R*(s_R - vel_R_norm) &
                                 & - rho_L*(s_L - vel_L_norm))*(f_low_Mach_zcoef(vel_L_rms, vel_R_rms, c_L, c_R) - 1._wp)
        end if

    end function f_low_Mach_pcorr_hllc

    !> The alternative low-Mach treatment of Thornber et al. JCP (2008) selected by low_Mach == 2: rather than correct the flux,
    !! blend the wave-normal velocities towards their mean before the wave speeds are computed, which is why this mutates its
    !! arguments and must be called ahead of s_L, s_R and s_S. The tangential velocities and vel_L/R_rms are deliberately left
    !! untouched.
    subroutine s_apply_low_Mach_velocity(vel_L_rms, vel_R_rms, c_L, c_R, vel_L_norm, vel_R_norm)

        $:GPU_ROUTINE(function_name='s_apply_low_Mach_velocity', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in)    :: vel_L_rms, vel_R_rms    !< Left and right squared velocity magnitudes
        real(wp), intent(in)    :: c_L, c_R                !< Left and right sound speeds
        real(wp), intent(inout) :: vel_L_norm, vel_R_norm  !< Left and right wave-normal velocities, blended in place
        real(wp)                :: zcoef, vel_L_tmp, vel_R_tmp

        zcoef = f_low_Mach_zcoef(vel_L_rms, vel_R_rms, c_L, c_R)

        vel_L_tmp = 5.e-1_wp*((vel_L_norm + vel_R_norm) + zcoef*(vel_L_norm - vel_R_norm))
        vel_R_tmp = 5.e-1_wp*((vel_L_norm + vel_R_norm) + zcoef*(vel_R_norm - vel_L_norm))

        vel_L_norm = vel_L_tmp
        vel_R_norm = vel_R_tmp

    end subroutine s_apply_low_Mach_velocity

    !> Interface-averaged state that the pressure-based wave-speed estimate reads. avg_state selects between the density-weighted
    !! Roe average, which costs eight square roots per face, and the plain arithmetic mean; unlike the other solver switches this
    !! one is not implied by the call site, so the dispatch stays here.
    subroutine s_compute_average_state(rho_L, rho_R, vel_L, vel_R, H_L, H_R, gamma_L, gamma_R, qv_L, qv_R, rho_avg, vel_avg_rms, &
                                       & H_avg, gamma_avg, qv_avg)

        $:GPU_ROUTINE(function_name='s_compute_average_state', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in) :: rho_L, rho_R      !< Left and right densities
        real(wp), intent(in) :: H_L, H_R          !< Left and right total enthalpies
        real(wp), intent(in) :: gamma_L, gamma_R  !< Left and right specific heat ratio functions
        real(wp), intent(in) :: qv_L, qv_R        !< Left and right reference energies
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(in) :: vel_L, vel_R
        #:else
            real(wp), dimension(num_vels), intent(in) :: vel_L, vel_R
        #:endif
        real(wp), intent(out) :: rho_avg, H_avg, gamma_avg, qv_avg
        real(wp), intent(out) :: vel_avg_rms  !< Squared magnitude of the averaged velocity, summed over all components
        integer               :: i

        vel_avg_rms = 0._wp

        if (avg_state == avg_state_roe) then
            rho_avg = sqrt(rho_L*rho_R)

            $:GPU_LOOP(parallelism='[seq]')
            do i = 1, num_vels
                vel_avg_rms = vel_avg_rms + (sqrt(rho_L)*vel_L(i) + sqrt(rho_R)*vel_R(i))**2._wp/(sqrt(rho_L) + sqrt(rho_R))**2._wp
            end do

            H_avg = (sqrt(rho_L)*H_L + sqrt(rho_R)*H_R)/(sqrt(rho_L) + sqrt(rho_R))
            gamma_avg = (sqrt(rho_L)*gamma_L + sqrt(rho_R)*gamma_R)/(sqrt(rho_L) + sqrt(rho_R))
            qv_avg = (sqrt(rho_L)*qv_L + sqrt(rho_R)*qv_R)/(sqrt(rho_L) + sqrt(rho_R))
        else
            rho_avg = 5.e-1_wp*(rho_L + rho_R)

            $:GPU_LOOP(parallelism='[seq]')
            do i = 1, num_vels
                vel_avg_rms = vel_avg_rms + (5.e-1_wp*(vel_L(i) + vel_R(i)))**2._wp
            end do

            H_avg = 5.e-1_wp*(H_L + H_R)
            gamma_avg = 5.e-1_wp*(gamma_L + gamma_R)
            qv_avg = 5.e-1_wp*(qv_L + qv_R)
        end if

    end subroutine s_compute_average_state

    !> Accumulate the hypoelastic stress contribution to the energies of the left and right Riemann states: mix the shear modulus
    !! over the fluids, add the elastic energy of each stress component (doubled for the shear components) on each side whose
    !! mixture modulus is non-negligible, then scale the returned moduli by the continuum damage state when damage is modeled
    !! (energy uses the undamaged modulus; the damaged moduli feed the callers' wave speeds). The elastic shear stresses are loaded
    !! from the state buffers by the caller, which reuses them for the stress fluxes and elastic wave speeds. The G > verysmall
    !! per-side gate is a deliberate maintainer ruling that replaces HLL's former hard-coded G > 1000 stability floor, retiring its
    !! "TODO take out if statement if stable without".
    subroutine s_compute_hypoelastic_interface_energy(nf, alpha_L, alpha_R, damage_L, damage_R, tau_e_L, tau_e_R, G_L, G_R, E_L, &
        & E_R)

        $:GPU_ROUTINE(function_name='s_compute_hypoelastic_interface_energy', parallelism='[seq]', cray_inline=True)

        integer, intent(in)                 :: nf                  !< Number of fluids to mix the shear modulus over
        real(wp), dimension(nf), intent(in) :: alpha_L, alpha_R    !< Left and right volume fractions
        real(wp), intent(in)                :: damage_L, damage_R  !< Continuum damage states (referenced only when cont_damage)
        real(wp), dimension(6), intent(in)  :: tau_e_L, tau_e_R    !< Left and right elastic shear stresses
        real(wp), intent(out)               :: G_L, G_R            !< Left and right mixture shear moduli
        real(wp), intent(inout)             :: E_L, E_R            !< Left and right state energies
        integer                             :: i                   !< Loop iterator

        G_L = 0._wp; G_R = 0._wp

        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, nf
            G_L = G_L + alpha_L(i)*Gs_rs(i)
            G_R = G_R + alpha_R(i)*Gs_rs(i)
        end do

        ! Elastic energy uses the undamaged modulus, so this loop precedes the damage scaling
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, eqn_idx%stress%end - eqn_idx%stress%beg + 1
            ! Elastic contribution to energy if G large enough
            if (G_L > verysmall) then
                E_L = E_L + (tau_e_L(i)*tau_e_L(i))/(4._wp*G_L)
                ! Double for shear stresses
                if (any(eqn_idx%stress%beg - 1 + i == shear_indices)) then
                    E_L = E_L + (tau_e_L(i)*tau_e_L(i))/(4._wp*G_L)
                end if
            end if
            if (G_R > verysmall) then
                E_R = E_R + (tau_e_R(i)*tau_e_R(i))/(4._wp*G_R)
                ! Double for shear stresses
                if (any(eqn_idx%stress%beg - 1 + i == shear_indices)) then
                    E_R = E_R + (tau_e_R(i)*tau_e_R(i))/(4._wp*G_R)
                end if
            end if
        end do

        ! Damage scaling applies only to the returned moduli (used for wave speeds)
        if (cont_damage) then
            G_L = G_L*max((1._wp - damage_L), 0._wp)
            G_R = G_R*max((1._wp - damage_R), 0._wp)
        end if

    end subroutine s_compute_hypoelastic_interface_energy

    !> Compute the advective part of the HLLC star-state momentum flux in the wave-normal direction (pressure excluded), used to
    !! assemble the geometrical source flux of the cylindrical and azimuthal sweeps.
    function f_compute_hllc_star_momentum_flux(rho_L, rho_R, vel_L_norm, vel_R_norm, s_M, s_P, s_S, xi_L, xi_R, xi_M, xi_P, &
        & dir_flg_norm) result(flux_mom)

        $:GPU_ROUTINE(function_name='f_compute_hllc_star_momentum_flux', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in) :: rho_L, rho_R            !< Left and right densities
        real(wp), intent(in) :: vel_L_norm, vel_R_norm  !< Left and right wave-normal velocities
        real(wp), intent(in) :: s_M, s_P, s_S           !< Clamped left/right and contact wave speeds
        real(wp), intent(in) :: xi_L, xi_R, xi_M, xi_P  !< Star-state compression factors and upwind selectors
        real(wp), intent(in) :: dir_flg_norm            !< Direction flag of the wave-normal direction
        real(wp)             :: flux_mom

        flux_mom = xi_M*(rho_L*(vel_L_norm*vel_L_norm + s_M*(xi_L*(dir_flg_norm*s_S + (1._wp - dir_flg_norm)*vel_L_norm) &
                         & - vel_L_norm))) + xi_P*(rho_R*(vel_R_norm*vel_R_norm + s_P*(xi_R*(dir_flg_norm*s_S + (1._wp &
                         & - dir_flg_norm)*vel_R_norm) - vel_R_norm)))

    end function f_compute_hllc_star_momentum_flux
#:enddef
