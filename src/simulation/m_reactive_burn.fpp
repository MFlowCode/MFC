!>
!! @file
!! @brief Contains module m_reactive_burn

#:include 'macros.fpp'
#:include 'case.fpp'

!> @brief Condensed-phase reactive burn: a pressure-driven programmed-burn source that converts a "reactant" fluid into a "product"
!! fluid on the multi-fluid model (num_fluids>=2, chemistry='F'). The two fluids share the same stiffened-gas EOS (gamma, pi_inf)
!! and differ only in their reference energy qv, so the reactant->product conversion releases (qv_reactant - qv_product) per unit
!! mass through the mixture EOS with no explicit energy source. Because the two fluids are mechanically identical the
!! volume-fraction swap is exact (the product volume fraction is the reaction progress), making this a reactive-Euler/ZND detonation
!! model expressed through the diffuse-interface framework. A shock raises the pressure above rburn_pign, the reactant burns, and
!! the energy release sustains the shock -- a self-propagating condensed-phase detonation.
module m_reactive_burn

    use m_global_parameters

    implicit none

    private; public :: s_compute_reactive_burn

contains

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
        real(wp)                                               :: rho, pres, lambda, rate, mdot, drive

        $:GPU_PARALLEL_LOOP(collapse=3, private='[rho, pres, lambda, rate, mdot, drive]', copyin='[bounds]')
        do z = bounds(3)%beg, bounds(3)%end
            do y = bounds(2)%beg, bounds(2)%end
                do x = bounds(1)%beg, bounds(1)%end
                    ! reactant is fluid 1, product is fluid 2
                    rho = q_cons_vf(eqn_idx%cont%beg)%sf(x, y, z) + q_cons_vf(eqn_idx%cont%beg + 1)%sf(x, y, z)
                    pres = q_prim_vf(eqn_idx%E)%sf(x, y, z)
                    lambda = q_prim_vf(eqn_idx%adv%beg + 1)%sf(x, y, z)  ! reaction progress = product volume fraction

                    ! pressure-driven programmed burn: fires only behind the shock (p > rburn_pign)
                    drive = (pres - rburn_pign)/rburn_pref
                    if (drive > 0._wp .and. lambda < 1._wp) then
                        rate = rburn_k*(1._wp - lambda)*drive**rburn_n  ! dlambda/dt
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

end module m_reactive_burn
