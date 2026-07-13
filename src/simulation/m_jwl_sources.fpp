!>
!! @file
!! @brief JWL reaction energy sources: kinematic program burn and afterburn.

#:include 'macros.fpp'
#:include 'case.fpp'

!> @brief Rocflu-style JWL reaction sources (after SPEC_RFLU_ModPBA and SPEC_RFLU_IntegrateChemSrcTerm, Illinois/NCSA license).
!!
!! Program burn: a front expands from the detonation point at pb_D_cj and
!! deposits jwl_Q over a band of width pb_width - exactly Q per unit explosive
!! mass, so no burn-fraction field is needed.
!!
!! Afterburn: advected progress b in [0, 1] caps the release at Y*jwl_q_ab.
!!   model 1 (mixing):    db/dt = (1-b)*(1-Y)/jwl_ab_tau
!!   model 2 (Arrhenius): db/dt = jwl_ab_A*p**jwl_ab_n*(1-b)*(1-Y)*exp(-jwl_ab_theta/T)
!!
!! Reactive burn (JWL++, Souers 2000): advected reaction progress lambda in
!! [0, 1] driven by the local pressure, dlambda/dt = jwl_G*p**jwl_b_exp*(1-lambda),
!! releasing jwl_Q as the explosive reacts. Unlike program burn the front is not
!! prescribed - it self-propagates from a hot spot in the initial condition.
!!
!! Both progress rates are linear in (1 - x), so each step applies the exact
!! frozen-coefficient solution x_new = 1 - (1-x)exp(-k dt) as an effective rate:
!! unconditionally bounded for stiff k*dt, identical to the explicit rate as
!! k*dt -> 0. This replaces the former hard (1-x)/dt clamp.
module m_jwl_sources

    use m_derived_types
    use m_global_parameters
    use m_ibm, only: ib_markers
    use m_jwl

    implicit none

    private; public :: s_compute_jwl_sources

contains

    !> Add JWL program-burn and afterburn sources to the energy (and afterburn progress) right-hand sides.
    subroutine s_compute_jwl_sources(q_cons_vf, rhs_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        real(wp) :: rho, dyn_p, e_sp, Y, b, dbdt, pres, T
        real(wp) :: q_det, r_det, r_front, lambda, dldt
        integer :: i, j, k, l
        integer :: solid_marker  !< /= 0 inside an immersed body: no reaction in solid or ghost cells

        if (prog_burn) then
            q_det = jwl_E0s(jwl_idx)/jwl_rho0s(jwl_idx)
            r_front = pb_D_cj*(mytime - pb_t_det)
            if (r_front > 0._wp) then
                $:GPU_PARALLEL_LOOP(collapse=3, private='[r_det, solid_marker]')
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            r_det = (x_cc(j) - pb_x_det)**2
                            if (n > 0) r_det = r_det + (y_cc(k) - pb_y_det)**2
                            if (p > 0) r_det = r_det + (z_cc(l) - pb_z_det)**2
                            r_det = sqrt(r_det)
                            solid_marker = 0
                            if (ib) solid_marker = ib_markers%sf(j, k, l)
                            if (solid_marker == 0 .and. r_det < r_front .and. r_det >= r_front - pb_width) then
                                rhs_vf(eqn_idx%E)%sf(j, k, l) = rhs_vf(eqn_idx%E)%sf(j, k, l) + q_cons_vf(jwl_idx)%sf(j, k, &
                                       & l)*q_det*pb_D_cj/pb_width
                            end if
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        end if

        if (jwl_afterburn) then
            $:GPU_PARALLEL_LOOP(collapse=3, private='[i, rho, dyn_p, e_sp, Y, b, dbdt, pres, T, lambda, solid_marker]')
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        rho = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, eqn_idx%cont%end
                            rho = rho + q_cons_vf(i)%sf(j, k, l)
                        end do
                        rho = max(rho, sgm_eps)
                        Y = min(max(q_cons_vf(jwl_idx)%sf(j, k, l)/rho, 0._wp), 1._wp)
                        b = min(max(q_cons_vf(eqn_idx%abn)%sf(j, k, l), 0._wp), 1._wp)

                        ! Both models vanish in the pure-products and pure-air limits
                        ! through the Y*(1 - Y) structure, and saturate as b -> 1.
                        ! No reaction inside an immersed body: solid and ghost cells hold
                        ! rebuilt (not evolved) states, so a rate there is phantom energy.
                        solid_marker = 0
                        if (ib) solid_marker = ib_markers%sf(j, k, l)
                        if (solid_marker == 0 .and. Y*(1._wp - Y)*(1._wp - b) > sgm_eps) then
                            if (jwl_ab_model == 1) then
                                dbdt = (1._wp - Y)/jwl_ab_tau
                            else
                                dyn_p = 0._wp
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = eqn_idx%mom%beg, eqn_idx%mom%end
                                    dyn_p = dyn_p + 0.5_wp*q_cons_vf(i)%sf(j, k, l)**2/rho
                                end do
                                e_sp = (q_cons_vf(eqn_idx%E)%sf(j, k, l) - dyn_p)/rho
                                ! Rate must see the same lambda-aware pressure/temperature the
                                ! solver computes (afterburn may combine with jwl_reactive).
                                lambda = 1._wp
                                if (jwl_reactive) lambda = min(max(q_cons_vf(eqn_idx%rxn)%sf(j, k, l), 0._wp), 1._wp)
                                call s_jwl_mix_state_er(rho, e_sp, Y, jwl_idx, pres, T, lambda=lambda)
                                dbdt = jwl_ab_A*pres**jwl_ab_n*(1._wp - Y)*exp(-jwl_ab_theta/T)
                            end if
                            ! Exact frozen-coefficient solution of db/dt = k(1-b) over the step
                            ! (same reasoning as the jwl_reactive rate below): b cannot overshoot
                            ! 1 within an RK substep, so the release stays capped at the
                            ! remaining Y*jwl_q_ab budget for arbitrarily stiff rates.
                            dbdt = (1._wp - b)*(1._wp - exp(-dbdt*dt))/dt
                            rhs_vf(eqn_idx%abn)%sf(j, k, l) = rhs_vf(eqn_idx%abn)%sf(j, k, l) + dbdt
                            rhs_vf(eqn_idx%E)%sf(j, k, l) = rhs_vf(eqn_idx%E)%sf(j, k, l) + q_cons_vf(jwl_idx)%sf(j, k, &
                                   & l)*jwl_q_ab*dbdt
                        end if
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        if (jwl_reactive) then
            q_det = jwl_E0s(jwl_idx)/jwl_rho0s(jwl_idx)
            $:GPU_PARALLEL_LOOP(collapse=3, private='[i, rho, dyn_p, e_sp, Y, lambda, dldt, pres, T, solid_marker]')
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        rho = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, eqn_idx%cont%end
                            rho = rho + q_cons_vf(i)%sf(j, k, l)
                        end do
                        rho = max(rho, sgm_eps)
                        Y = min(max(q_cons_vf(jwl_idx)%sf(j, k, l)/rho, 0._wp), 1._wp)
                        lambda = min(max(q_cons_vf(eqn_idx%rxn)%sf(j, k, l), 0._wp), 1._wp)

                        ! Pressure-driven reaction, active only in unreacted explosive.
                        ! No reaction inside an immersed body (see afterburn note above).
                        solid_marker = 0
                        if (ib) solid_marker = ib_markers%sf(j, k, l)
                        if (solid_marker == 0 .and. Y*(1._wp - lambda) > sgm_eps) then
                            dyn_p = 0._wp
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = eqn_idx%mom%beg, eqn_idx%mom%end
                                dyn_p = dyn_p + 0.5_wp*q_cons_vf(i)%sf(j, k, l)**2/rho
                            end do
                            e_sp = (q_cons_vf(eqn_idx%E)%sf(j, k, l) - dyn_p)/rho
                            call s_jwl_mix_state_er(rho, e_sp, Y, jwl_idx, pres, T, lambda=lambda)
                            dldt = jwl_G*pres**jwl_b_exp
                            ! Exact frozen-pressure solution of dl/dt = k(1-l) over the step,
                            ! written as an effective rate. It saturates smoothly at (1-l)/dt
                            ! for stiff k*dt and recovers the explicit rate as k*dt -> 0, so
                            ! lambda stays bounded for ANY G/p/dt: the stiff-source dt collapse
                            ! seen in multi-D reactive runs cannot occur, and G no longer needs
                            ! de-tuning below its 1D-calibrated value to survive flow focusing.
                            dldt = (1._wp - lambda)*(1._wp - exp(-dldt*dt))/dt
                            rhs_vf(eqn_idx%rxn)%sf(j, k, l) = rhs_vf(eqn_idx%rxn)%sf(j, k, l) + dldt
                            rhs_vf(eqn_idx%E)%sf(j, k, l) = rhs_vf(eqn_idx%E)%sf(j, k, l) + q_cons_vf(jwl_idx)%sf(j, k, &
                                   & l)*q_det*dldt
                        end if
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_compute_jwl_sources

end module m_jwl_sources
