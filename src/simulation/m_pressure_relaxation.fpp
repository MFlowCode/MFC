!>
!! @file
!! @brief Contains module m_pressure_relaxation

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief Pressure relaxation for the six-equation multi-component model via Newton--Raphson equilibration and volume-fraction
!! correction
module m_pressure_relaxation

    use m_derived_types
    use m_global_parameters
    use m_variables_conversion, only: s_convert_species_to_mixture_variables_kernel, f_pressure, s_phase_internal_energy, &
        & s_phase_coefficients, s_phase_density_on_isentrope, f_is_state_dependent

    implicit none

    private; public :: s_pressure_relaxation_procedure, s_report_pressure_relaxation

    !> Cell updates in which the Newton below stopped on its iteration cap rather than on the tolerance. Counted because the answer
    !! then depends on that cap, and nothing else in the code says so.
    integer  :: n_hit_cap = 0
    real(wp) :: worst_residual = 0._wp  !< Largest |f| left behind when the cap was reached

contains

    !> The main pressure relaxation procedure
    subroutine s_pressure_relaxation_procedure(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        integer                                                :: i, j, k, l

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: alpha_rho, alpha
        #:else
            real(wp), dimension(num_fluids) :: alpha_rho, alpha
        #:endif
        real(wp) :: rho, gamma, pi_inf, qv_mix
        integer  :: hit_cap, hit_cap_sum
        real(wp) :: resid, resid_max

        ! Formed here, not one call deeper: CCE OpenACC accepts a num_fluids-sized array passed to a device routine from a
        ! parallel-loop body, and rejects the same call from inside another acc routine seq.
        hit_cap_sum = 0
        resid_max = 0._wp
        $:GPU_PARALLEL_LOOP(private='[i, j, k, l, alpha_rho, alpha, rho, gamma, pi_inf, qv_mix, hit_cap, resid]', collapse=3, &
                            & reduction='[[hit_cap_sum], [resid_max]]', reductionOp='[+, max]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    if (mpp_lim) call s_correct_volume_fractions(q_cons_vf, j, k, l)

                    hit_cap = 0
                    resid = 0._wp
                    if (s_needs_pressure_relaxation(q_cons_vf, j, k, l)) then
                        call s_equilibrate_pressure(q_cons_vf, j, k, l, hit_cap, resid)
                    end if
                    hit_cap_sum = hit_cap_sum + hit_cap
                    resid_max = max(resid_max, resid)

                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        alpha_rho(i) = q_cons_vf(i)%sf(j, k, l)
                        alpha(i) = q_cons_vf(eqn_idx%E + i)%sf(j, k, l)
                    end do

                    call s_convert_species_to_mixture_variables_kernel(rho, gamma, pi_inf, qv_mix, alpha, alpha_rho)

                    call s_correct_internal_energies(q_cons_vf, j, k, l, rho, gamma, pi_inf, qv_mix)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        n_hit_cap = n_hit_cap + hit_cap_sum
        worst_residual = max(worst_residual, resid_max)

    end subroutine s_pressure_relaxation_procedure

    !> One line at the end of a run if the equilibration ever stopped on its iteration cap. Silence means every cell reached the
    !! tolerance.
    impure subroutine s_report_pressure_relaxation

        if (n_hit_cap > 0 .and. proc_rank == 0) then
            print '(A,I0,A,ES10.3,A)', ' Pressure relaxation reached its iteration cap in ', n_hit_cap, &
                & ' cell updates; worst residual left behind ', worst_residual, ' against a 1e-10 tolerance.'
        end if

    end subroutine s_report_pressure_relaxation

    !> Check if pressure relaxation is needed for this cell
    logical function s_needs_pressure_relaxation(q_cons_vf, j, k, l)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        integer, intent(in)                                 :: j, k, l
        integer                                             :: i

        s_needs_pressure_relaxation = .true.
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            if (q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l) > (1._wp - sgm_eps)) then
                s_needs_pressure_relaxation = .false.
            end if
        end do

    end function s_needs_pressure_relaxation

    !> Correct volume fractions to physical bounds
    subroutine s_correct_volume_fractions(q_cons_vf, j, k, l)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        integer, intent(in)                                    :: j, k, l
        real(wp)                                               :: sum_alpha
        integer                                                :: i

        sum_alpha = 0._wp
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            if ((q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l) < 0._wp) .or. (q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, &
                & l) < 0._wp)) then
                q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l) = 0._wp
                q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l) = 0._wp
                q_cons_vf(i + eqn_idx%int_en%beg - 1)%sf(j, k, l) = 0._wp
            end if
            if (q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l) > 1._wp) q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l) = 1._wp
            sum_alpha = sum_alpha + q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l)
        end do

        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l) = q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l)/sum_alpha
        end do

    end subroutine s_correct_volume_fractions

    !> Main pressure equilibration using Newton-Raphson
    subroutine s_equilibrate_pressure(q_cons_vf, j, k, l, hit_cap, resid)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        integer, intent(in)                                    :: j, k, l
        integer, intent(out)                                   :: hit_cap
        real(wp), intent(out)                                  :: resid
        real(wp)                                               :: pres_relax, f_pres, df_pres
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: pres_K_init, rho_K_init, rho_K_s
        #:else
            real(wp), dimension(num_fluids) :: pres_K_init, rho_K_init, rho_K_s
        #:endif
        real(wp)           :: gamma_K, pi_inf_K, dpi_K, dgamma_K, c2_K, alpha_i, alpha_rho_i, rho_i, p_i, rho_s_i
        integer, parameter :: MAX_ITER = 50
        ! Pressure relaxation convergence tolerance
        real(wp), parameter :: TOLERANCE = 1.e-10_wp
        integer             :: iter, i

        ! Initialize pressures
        pres_relax = 0._wp
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            if (q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l) > sgm_eps) then
                ! Phasic internal energy carries the formation energy: alpha_rho_k*qv_k must be
                ! removed before inverting the stiffened-gas EOS, or a nonzero qv inflates the
                ! phasic pressure by rho_k*qv_k/gamma_k (this is what breaks the reactive burn).
                alpha_rho_i = q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l)
                alpha_i = q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l)
                call s_phase_coefficients(alpha_rho_i, alpha_i, i, rho_i, gamma_K, pi_inf_K, dpi_K, dgamma_K)
                rho_K_init(i) = rho_i
                pres_K_init(i) = ((q_cons_vf(i + eqn_idx%int_en%beg - 1)%sf(j, k, l) - q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, &
                            & k, l)*qvs(i))/q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l) - pi_inf_K)/gamma_K
                if (.not. f_is_state_dependent(i)) then
                    if (pres_K_init(i) <= -(1._wp - 1.e-8_wp)*isentrope_B(i) + 1.e-8_wp) pres_K_init(i) = -(1._wp - 1.e-8_wp) &
                        & *isentrope_B(i) + 1.e-8_wp
                end if
            else
                pres_K_init(i) = 0._wp
            end if
            pres_relax = pres_relax + q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l)*pres_K_init(i)
        end do

        ! Newton-Raphson iteration
        f_pres = 1.e-9_wp
        df_pres = 1.e9_wp
        $:GPU_LOOP(parallelism='[seq]')
        do iter = 0, MAX_ITER - 1
            if (abs(f_pres) > TOLERANCE) then
                pres_relax = pres_relax - f_pres/df_pres

                ! Enforce pressure bounds
                do i = 1, num_fluids
                    if (.not. f_is_state_dependent(i)) then
                        if (pres_relax <= -(1._wp - 1.e-8_wp)*isentrope_B(i) + 1.e-8_wp) pres_relax = -(1._wp - 1.e-8_wp) &
                            & *isentrope_B(i) + 1.e-8_wp
                    end if
                end do

                ! Newton-Raphson step
                f_pres = -1._wp
                df_pres = 0._wp
                $:GPU_LOOP(parallelism='[seq]')
                do i = 1, num_fluids
                    if (q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l) > sgm_eps .and. any_state_dependent_eos) then
                        rho_i = rho_K_init(i)
                        p_i = pres_K_init(i)
                        call s_phase_density_on_isentrope(i, rho_i, p_i, pres_relax, rho_s_i, c2_K)
                        rho_K_s(i) = rho_s_i
                        f_pres = f_pres + q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l)/rho_K_s(i)
                        df_pres = df_pres - q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l)/(rho_K_s(i)**2*c2_K)
                    else if (q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l) > sgm_eps) then
                        ! Isentropic relation: rho = rho0 * (p/p0)^(1/gamma), Saurel et al. JFM (2009)
                        rho_K_s(i) = q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l)/max(q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, &
                                & k, l), &
                                & sgm_eps)*((pres_relax + isentrope_B(i))/(pres_K_init(i) + isentrope_B(i)))**(1._wp/isentrope_n(i))
                        f_pres = f_pres + q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l)/rho_K_s(i)
                        df_pres = df_pres - q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, &
                                                      & l)/(isentrope_n(i)*rho_K_s(i)*(pres_relax + isentrope_B(i)))
                    end if
                end do
            end if
        end do

        ! Written as .not. (<=) so a NaN residual is reported as a miss, like a diverged one.
        hit_cap = 0
        resid = 0._wp
        if (.not. (abs(f_pres) <= TOLERANCE)) then
            hit_cap = 1
            ! A NaN residual has to be mapped, not passed on: max() with a NaN keeps the other operand, so the
            ! reduction would report the miss as zero.
            resid = abs(f_pres)
            if (f_pres /= f_pres) resid = huge(1._wp)
        end if

        ! Update volume fractions. The Newton above often stops on the iteration cap rather than on the
        ! tolerance, and the answer then depends on that cap, so an unconverged-but-physical density is still
        ! used -- the alternative would move every six-equation result. What is refused is a density that is
        ! not a usable number: rho_K_s <= 0 or NaN, which is how one bad cell became a NaN a few steps later.
        ! The comparison is written as .not. (> 0) so a NaN takes the same path as a non-positive value.
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            if (q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l) > sgm_eps .and. .not. (rho_K_s(i) > 0._wp)) return
        end do

        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            if (q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l) > sgm_eps) q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, &
                & l) = q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l)/rho_K_s(i)
        end do

    end subroutine s_equilibrate_pressure

    !> Correct internal energies using equilibrated pressure
    subroutine s_correct_internal_energies(q_cons_vf, j, k, l, rho, gamma, pi_inf, qv_mix)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        integer, intent(in)                                    :: j, k, l
        real(wp), intent(in)                                   :: rho, gamma, pi_inf, qv_mix
        real(wp)                                               :: dyn_pres, pres_relax, alpha_i, alpha_rho_i, e_i
        integer                                                :: i

        dyn_pres = 0._wp
        $:GPU_LOOP(parallelism='[seq]')
        do i = eqn_idx%mom%beg, eqn_idx%mom%end
            dyn_pres = dyn_pres + 5.e-1_wp*q_cons_vf(i)%sf(j, k, l)*q_cons_vf(i)%sf(j, k, l)/max(rho, sgm_eps)
        end do

        pres_relax = f_pressure(q_cons_vf(eqn_idx%E)%sf(j, k, l) - dyn_pres, gamma, pi_inf, qv_mix)

        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            alpha_i = q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l)
            alpha_rho_i = q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l)
            call s_phase_internal_energy(pres_relax, alpha_i, alpha_rho_i, i, e_i)
            q_cons_vf(i + eqn_idx%int_en%beg - 1)%sf(j, k, l) = e_i
        end do

    end subroutine s_correct_internal_energies

end module m_pressure_relaxation
