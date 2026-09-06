!>
!! @file
!! @brief Contains module m_reactive_burn

#:include 'macros.fpp'
#:include 'case.fpp'

!> @brief Condensed-phase reactive burn: a pressure-driven programmed-burn source that converts a "reactant" fluid into a "product"
!! fluid on the multi-fluid model (num_fluids=2, chemistry='F'). The two fluids share the same stiffened-gas EOS (gamma, pi_inf) and
!! differ only in their reference energy qv, so the reactant->product conversion releases (qv_reactant - qv_product) per unit mass
!! through the mixture EOS with no explicit energy source. Because the two fluids are mechanically identical the volume-fraction
!! swap is exact (the product volume fraction is the reaction progress), making this a reactive-Euler/ZND detonation model expressed
!! through the diffuse-interface framework. A shock raises the pressure above rburn%pign, the reactant burns, and the energy release
!! sustains the shock -- a self-propagating condensed-phase detonation.
module m_reactive_burn

    use m_global_parameters
    use m_variables_conversion, only: f_sg_thermal, s_compute_mixture_coefficients, f_pressure

    implicit none

    private; public :: s_compute_reactive_burn, s_reactive_burn_substep

contains

    !> Programmed-burn rate dlambda/dt for one cell state. Both the RHS source and the operator-split integrator call this, so the
    !! rate law is stated once.
    !! @param pres             Mixture pressure
    !! @param lambda           Reaction progress, i.e. the product volume fraction
    !! @param alpha_rho_react  Reactant partial density, for the optional Arrhenius factor
    !! @param alpha_react      Reactant volume fraction, for the optional Arrhenius factor
    !! @param rate             dlambda/dt; zero below the ignition pressure and once the reactant is spent
    subroutine s_burn_rate(pres, lambda, alpha_rho_react, alpha_react, rate)

        $:GPU_ROUTINE(function_name='s_burn_rate', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in)  :: pres, lambda, alpha_rho_react, alpha_react
        real(wp), intent(out) :: rate
        real(wp)              :: drive

        ! Pressure-driven programmed burn: fires only behind the shock (p > rburn%pign).
        drive = (pres - rburn%pign)/rburn%pref
        if (drive > 0._wp .and. lambda < 1._wp) then
            rate = rburn%k*(1._wp - lambda)*drive**rburn%n
            ! Optional Arrhenius dependence on the reactant phasic temperature from the stiffened-gas
            ! EOS. rburn%ta = 0, the default, leaves the pure pressure-driven rate unchanged.
            if (rburn%ta > 0._wp) then
                rate = rate*exp(-rburn%ta/f_sg_thermal(pres, alpha_rho_react/alpha_react, isentrope_n(1), isentrope_B(1), cvs(1)))
            end if
        else
            rate = 0._wp
        end if

    end subroutine s_burn_rate

    !> Add the programmed-burn reaction source to the continuity and volume-fraction RHS.
    !! @param rhs_vf     Right-hand-side accumulator (inout)
    !! @param q_cons_vf  Conserved variables (partial densities live here)
    !! @param q_prim_vf  Primitive variables (pressure and volume fractions live here)
    !! @param bounds     Interior cell bounds
    subroutine s_compute_reactive_burn(rhs_vf, q_cons_vf, q_prim_vf, bounds)

        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        type(scalar_field), dimension(sys_size), intent(in)    :: q_cons_vf, q_prim_vf
        type(int_bounds_info), dimension(1:3), intent(in)      :: bounds
        integer                                                :: x, y, z
        real(wp)                                               :: rho, pres, lambda, rate, mdot

        $:GPU_PARALLEL_LOOP(collapse=3, private='[rho, pres, lambda, rate, mdot]', copyin='[bounds]')
        do z = bounds(3)%beg, bounds(3)%end
            do y = bounds(2)%beg, bounds(2)%end
                do x = bounds(1)%beg, bounds(1)%end
                    ! reactant is fluid 1, product is fluid 2
                    rho = q_cons_vf(eqn_idx%cont%beg)%sf(x, y, z) + q_cons_vf(eqn_idx%cont%beg + 1)%sf(x, y, z)
                    pres = q_prim_vf(eqn_idx%E)%sf(x, y, z)
                    lambda = q_prim_vf(eqn_idx%adv%beg + 1)%sf(x, y, z)  ! reaction progress = product volume fraction

                    call s_burn_rate(pres, lambda, q_cons_vf(eqn_idx%cont%beg)%sf(x, y, z), q_prim_vf(eqn_idx%adv%beg)%sf(x, y, &
                                     & z), rate)
                    if (rate > 0._wp) then
                        mdot = rho*rate  ! mass reactant -> product

                        ! continuity: reactant loses mass, product gains it
                        rhs_vf(eqn_idx%cont%beg)%sf(x, y, z) = rhs_vf(eqn_idx%cont%beg)%sf(x, y, z) - mdot
                        rhs_vf(eqn_idx%cont%beg + 1)%sf(x, y, z) = rhs_vf(eqn_idx%cont%beg + 1)%sf(x, y, z) + mdot

                        ! volume fraction: exact swap (fluids share the EOS), so d(alpha)/dt = +/- rate
                        rhs_vf(eqn_idx%adv%beg)%sf(x, y, z) = rhs_vf(eqn_idx%adv%beg)%sf(x, y, z) - rate
                        rhs_vf(eqn_idx%adv%beg + 1)%sf(x, y, z) = rhs_vf(eqn_idx%adv%beg + 1)%sf(x, y, z) + rate
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_reactive_burn

    !> Operator-split alternative to s_compute_reactive_burn, used when rburn%substeps > 0. The flow is frozen and the burn ODE is
    !! integrated over one time step in equal sub-steps, so the reaction time scale is decoupled from the acoustic CFL. The mixture
    !! pressure is re-evaluated from the frozen internal energy each sub-step, which is what carries the rate's own feedback: the
    !! coefficients move as the reactant becomes product.
    !! @param q_cons_vf  Conserved variables, updated in place
    !! @param dtime      Time step to integrate across
    !! @param bounds     Interior cell bounds
    subroutine s_reactive_burn_substep(q_cons_vf, dtime, bounds)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        real(wp), intent(in)                                   :: dtime
        type(int_bounds_info), dimension(1:3), intent(in)      :: bounds
        integer                                                :: x, y, z, i, sub
        real(wp)                                               :: rho, pres, lambda, rate
        real(wp)                                               :: dt_sub, e_int, gamma_mix, pi_inf_mix, qv_mix
        real(wp)                                               :: rho_mix, dlambda, dmass

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: alpha_rho, alpha
        #:else
            real(wp), dimension(num_fluids) :: alpha_rho, alpha
        #:endif

        dt_sub = dtime/real(rburn%substeps, wp)

        $:GPU_PARALLEL_LOOP(collapse=3, private='[alpha_rho, alpha, rho, pres, lambda, rate, e_int, gamma_mix, pi_inf_mix, &
                            & qv_mix, rho_mix, dlambda, dmass, i, sub]', copyin='[bounds, dt_sub]')
        do z = bounds(3)%beg, bounds(3)%end
            do y = bounds(2)%beg, bounds(2)%end
                do x = bounds(1)%beg, bounds(1)%end
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        alpha_rho(i) = q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(x, y, z)
                        alpha(i) = q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(x, y, z)
                    end do
                    rho = alpha_rho(1) + alpha_rho(2)

                    ! Internal energy per unit volume is what the burn conserves: it moves mass and volume
                    ! between the phases, and the qv difference surfaces as pressure through the mixture EOS.
                    e_int = q_cons_vf(eqn_idx%E)%sf(x, y, z)
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = eqn_idx%mom%beg, eqn_idx%mom%end
                        e_int = e_int - 0.5_wp*q_cons_vf(i)%sf(x, y, z)**2/rho
                    end do

                    $:GPU_LOOP(parallelism='[seq]')
                    do sub = 1, rburn%substeps
                        lambda = alpha(2)
                        call s_compute_mixture_coefficients(alpha_rho, alpha, rho_mix, gamma_mix, pi_inf_mix, qv_mix)
                        pres = f_pressure(e_int, gamma_mix, pi_inf_mix, qv_mix)
                        call s_burn_rate(pres, lambda, alpha_rho(1), alpha(1), rate)
                        if (rate <= 0._wp) exit
                        ! A sub-step longer than the reaction time would carry the progress variable past one.
                        ! Stop it there and hand over the reactant's remaining mass in the same sub-step: capping
                        ! the two independently strands mass at zero volume, and the EOS divides one by the other.
                        if (rate*dt_sub >= 1._wp - lambda) then
                            dlambda = 1._wp - lambda
                            dmass = alpha_rho(1)
                        else
                            dlambda = rate*dt_sub
                            dmass = min(rho*dlambda, alpha_rho(1))
                        end if
                        alpha(1) = alpha(1) - dlambda
                        alpha(2) = alpha(2) + dlambda
                        alpha_rho(1) = alpha_rho(1) - dmass
                        alpha_rho(2) = alpha_rho(2) + dmass
                    end do

                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(x, y, z) = alpha_rho(i)
                        q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(x, y, z) = alpha(i)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_reactive_burn_substep

end module m_reactive_burn
