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
    use m_variables_conversion, only: s_jwl_energy_pr, s_jwl_pcold, s_jwl_sound_speed_squared, jwl_idx

    implicit none

    private; public :: s_pressure_relaxation_procedure, s_initialize_pressure_relaxation_module, &
        & s_finalize_pressure_relaxation_module

    real(wp), allocatable, dimension(:,:) :: Res_pr
    $:GPU_DECLARE(create='[Res_pr]')

contains

    !> Initialize the pressure relaxation module
    impure subroutine s_initialize_pressure_relaxation_module

        integer :: i, j

        if (viscous) then
            @:ALLOCATE(Res_pr(1:2, 1:Re_size_max))
            do i = 1, 2
                do j = 1, Re_size(i)
                    Res_pr(i, j) = fluid_pp(Re_idx(i, j))%Re(i)
                end do
            end do
            $:GPU_UPDATE(device='[Res_pr, Re_idx, Re_size]')
        end if

    end subroutine s_initialize_pressure_relaxation_module

    !> Finalize the pressure relaxation module
    impure subroutine s_finalize_pressure_relaxation_module

        if (viscous) then
            @:DEALLOCATE(Res_pr)
        end if

    end subroutine s_finalize_pressure_relaxation_module

    !> The main pressure relaxation procedure
    subroutine s_pressure_relaxation_procedure(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        integer                                                :: j, k, l

        $:GPU_PARALLEL_LOOP(private='[j, k, l]', collapse=3)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    call s_relax_cell_pressure(q_cons_vf, j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_pressure_relaxation_procedure

    !> Process pressure relaxation for a single cell
    subroutine s_relax_cell_pressure(q_cons_vf, j, k, l)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        integer, intent(in)                                    :: j, k, l
        real(wp)                                               :: pres_relax

        ! Volume fraction correction
        if (mpp_lim) call s_correct_volume_fractions(q_cons_vf, j, k, l)

        ! Pressure equilibration
        if (s_needs_pressure_relaxation(q_cons_vf, j, k, l)) then
            call s_equilibrate_pressure(q_cons_vf, j, k, l, pres_relax)
        else
            call s_cell_average_pressure(q_cons_vf, j, k, l, pres_relax)
        end if

        ! Internal energy correction
        call s_correct_internal_energies(q_cons_vf, j, k, l, pres_relax)

    end subroutine s_relax_cell_pressure

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

    !> Recover a phase pressure from conserved phase volume, partial density, and internal-energy density.
    subroutine s_phase_pressure_from_energy(fluid_id, alpha, alpha_rho, alpha_energy, pres)

        $:GPU_ROUTINE(parallelism='[seq]')

        integer, intent(in)   :: fluid_id
        real(wp), intent(in)  :: alpha, alpha_rho, alpha_energy
        real(wp), intent(out) :: pres
        real(wp)              :: rho_phase, pcold, K_jwl

        if (alpha <= sgm_eps .or. alpha_rho <= sgm_eps) then
            pres = 0._wp
        else if (eos_idxs(fluid_id) == 2) then
            rho_phase = max(alpha_rho/alpha, sgm_eps)
            call s_jwl_pcold(rho_phase, jwl_As(fluid_id), jwl_Bs(fluid_id), jwl_R1s(fluid_id), jwl_R2s(fluid_id), &
                             & jwl_omegas(fluid_id), jwl_rho0s(fluid_id), pcold)
            K_jwl = jwl_rho0s(fluid_id)/max(jwl_omegas(fluid_id), sgm_eps)
            pres = pcold + alpha_energy/max(alpha*K_jwl, sgm_eps)
            pres = max(pres, 1._wp)
        else
            pres = ((alpha_energy - alpha_rho*qvs(fluid_id))/max(alpha, sgm_eps) - pi_infs(fluid_id))/gammas(fluid_id)
            if (pres <= -(1._wp - 1.e-8_wp)*ps_inf(fluid_id) + 1.e-8_wp) then
                pres = -(1._wp - 1.e-8_wp)*ps_inf(fluid_id) + 1.e-8_wp
            end if
        end if

    end subroutine s_phase_pressure_from_energy

    !> Phase internal-energy density at a common pressure.
    subroutine s_phase_energy_from_pressure(fluid_id, alpha, alpha_rho, pres, alpha_energy)

        $:GPU_ROUTINE(parallelism='[seq]')

        integer, intent(in)   :: fluid_id
        real(wp), intent(in)  :: alpha, alpha_rho, pres
        real(wp), intent(out) :: alpha_energy
        real(wp)              :: rho_phase, e_phase

        if (alpha <= sgm_eps .or. alpha_rho <= sgm_eps) then
            alpha_energy = 0._wp
        else if (eos_idxs(fluid_id) == 2) then
            rho_phase = max(alpha_rho/alpha, sgm_eps)
            call s_jwl_energy_pr(rho_phase, pres, 1._wp, 1._wp, jwl_As(fluid_id), jwl_Bs(fluid_id), jwl_R1s(fluid_id), &
                                 & jwl_R2s(fluid_id), jwl_omegas(fluid_id), jwl_rho0s(fluid_id), jwl_E0s(fluid_id), &
                                 & jwl_air_e0s(fluid_id), jwl_air_rho0s(fluid_id), jwl_air_gammas(fluid_id), e_phase)
            alpha_energy = alpha_rho*e_phase
        else
            alpha_energy = alpha*(gammas(fluid_id)*pres + pi_infs(fluid_id)) + alpha_rho*qvs(fluid_id)
        end if

    end subroutine s_phase_energy_from_pressure

    !> Density reached by isentropically moving a phase from p0 to p.
    subroutine s_phase_density_isentrope(fluid_id, rho_init, pres_init, pres, rho_s, c2_s)

        $:GPU_ROUTINE(parallelism='[seq]')

        integer, intent(in)   :: fluid_id
        real(wp), intent(in)  :: rho_init, pres_init, pres
        real(wp), intent(out) :: rho_s, c2_s
        integer, parameter    :: JWL_ISENTROPE_STEPS = 8
        integer               :: step
        real(wp)              :: dp, p0, p1, r0, k1, k2, k3, k4

        if (eos_idxs(fluid_id) == 2) then
            rho_s = max(rho_init, sgm_eps)
            p0 = max(pres_init, 1._wp)
            dp = (max(pres, 1._wp) - p0)/real(JWL_ISENTROPE_STEPS, wp)
            do step = 1, JWL_ISENTROPE_STEPS
                r0 = rho_s
                p1 = p0 + 5.e-1_wp*dp
                call s_jwl_sound_speed_squared(max(r0, sgm_eps), max(p0, 1._wp), jwl_As(fluid_id), jwl_Bs(fluid_id), &
                                               & jwl_R1s(fluid_id), jwl_R2s(fluid_id), jwl_omegas(fluid_id), jwl_rho0s(fluid_id), &
                                               & c2_s)
                k1 = 1._wp/max(c2_s, sgm_eps)
                call s_jwl_sound_speed_squared(max(r0 + 5.e-1_wp*dp*k1, sgm_eps), max(p1, 1._wp), jwl_As(fluid_id), &
                                               & jwl_Bs(fluid_id), jwl_R1s(fluid_id), jwl_R2s(fluid_id), jwl_omegas(fluid_id), &
                                               & jwl_rho0s(fluid_id), c2_s)
                k2 = 1._wp/max(c2_s, sgm_eps)
                call s_jwl_sound_speed_squared(max(r0 + 5.e-1_wp*dp*k2, sgm_eps), max(p1, 1._wp), jwl_As(fluid_id), &
                                               & jwl_Bs(fluid_id), jwl_R1s(fluid_id), jwl_R2s(fluid_id), jwl_omegas(fluid_id), &
                                               & jwl_rho0s(fluid_id), c2_s)
                k3 = 1._wp/max(c2_s, sgm_eps)
                call s_jwl_sound_speed_squared(max(r0 + dp*k3, sgm_eps), max(p0 + dp, 1._wp), jwl_As(fluid_id), jwl_Bs(fluid_id), &
                                               & jwl_R1s(fluid_id), jwl_R2s(fluid_id), jwl_omegas(fluid_id), jwl_rho0s(fluid_id), &
                                               & c2_s)
                k4 = 1._wp/max(c2_s, sgm_eps)
                rho_s = max(r0 + dp*(k1 + 2._wp*k2 + 2._wp*k3 + k4)/6._wp, sgm_eps)
                p0 = p0 + dp
            end do
            call s_jwl_sound_speed_squared(rho_s, max(pres, 1._wp), jwl_As(fluid_id), jwl_Bs(fluid_id), jwl_R1s(fluid_id), &
                                           & jwl_R2s(fluid_id), jwl_omegas(fluid_id), jwl_rho0s(fluid_id), c2_s)
        else
            rho_s = rho_init*((pres + ps_inf(fluid_id))/(pres_init + ps_inf(fluid_id)))**(1._wp/gs_min(fluid_id))
            c2_s = gs_min(fluid_id)*(pres + ps_inf(fluid_id))/max(rho_s, sgm_eps)
        end if

    end subroutine s_phase_density_isentrope

    !> Volume-fraction average of current phase pressures, used when a cell is already single-phase.
    subroutine s_cell_average_pressure(q_cons_vf, j, k, l, pres)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        integer, intent(in)                                 :: j, k, l
        real(wp), intent(out)                               :: pres
        real(wp)                                            :: alpha_i, phase_pres
        integer                                             :: i

        pres = 0._wp
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            alpha_i = q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l)
            if (alpha_i > sgm_eps) then
                call s_phase_pressure_from_energy(i, alpha_i, q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l), &
                                                  & q_cons_vf(i + eqn_idx%int_en%beg - 1)%sf(j, k, l), phase_pres)
                pres = pres + alpha_i*phase_pres
            end if
        end do

    end subroutine s_cell_average_pressure

    !> Main pressure equilibration using Newton-Raphson
    subroutine s_equilibrate_pressure(q_cons_vf, j, k, l, pres_relaxed)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        integer, intent(in)                                    :: j, k, l
        real(wp), intent(out)                                  :: pres_relaxed
        real(wp)                                               :: pres_relax, f_pres, df_pres
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: pres_K_init, rho_K_s
        #:else
            real(wp), dimension(num_fluids) :: pres_K_init, rho_K_s
        #:endif
        real(wp)           :: alpha_i, alpha_rho_i, alpha_energy_i, rho_i_init, c2_i
        integer, parameter :: MAX_ITER = 50
        ! Pressure relaxation convergence tolerance
        real(wp), parameter :: TOLERANCE = 1.e-10_wp
        integer             :: iter, i

        ! Initialize pressures
        pres_relax = 0._wp
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            alpha_i = q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l)
            alpha_rho_i = q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l)
            alpha_energy_i = q_cons_vf(i + eqn_idx%int_en%beg - 1)%sf(j, k, l)
            if (alpha_i > sgm_eps) then
                call s_phase_pressure_from_energy(i, alpha_i, alpha_rho_i, alpha_energy_i, pres_K_init(i))
            else
                pres_K_init(i) = 0._wp
            end if
            pres_relax = pres_relax + alpha_i*pres_K_init(i)
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
                    if (eos_idxs(i) == 2) then
                        pres_relax = max(pres_relax, 1._wp)
                    else if (pres_relax <= -(1._wp - 1.e-8_wp)*ps_inf(i) + 1.e-8_wp) then
                        pres_relax = -(1._wp - 1.e-8_wp)*ps_inf(i) + 1.e-8_wp
                    end if
                end do

                ! Newton-Raphson step
                f_pres = -1._wp
                df_pres = 0._wp
                $:GPU_LOOP(parallelism='[seq]')
                do i = 1, num_fluids
                    alpha_i = q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l)
                    alpha_rho_i = q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l)
                    if (alpha_i > sgm_eps) then
                        rho_i_init = alpha_rho_i/max(alpha_i, sgm_eps)
                        call s_phase_density_isentrope(i, rho_i_init, pres_K_init(i), pres_relax, rho_K_s(i), c2_i)
                        f_pres = f_pres + alpha_rho_i/rho_K_s(i)
                        df_pres = df_pres - alpha_rho_i/(rho_K_s(i)*rho_K_s(i)*max(c2_i, sgm_eps))
                    end if
                end do
            end if
        end do

        ! Update volume fractions
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            if (q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l) > sgm_eps) q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, &
                & l) = q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l)/rho_K_s(i)
        end do

        pres_relaxed = pres_relax

    end subroutine s_equilibrate_pressure

    !> Correct internal energies using equilibrated pressure
    subroutine s_correct_internal_energies(q_cons_vf, j, k, l, pres_relax)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        integer, intent(in)                                    :: j, k, l
        real(wp), intent(in)                                   :: pres_relax
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: alpha_rho, alpha
        #:else
            real(wp), dimension(num_fluids) :: alpha_rho, alpha
        #:endif
        real(wp)               :: rho, dyn_pres, gamma, pi_inf, sum_alpha, pres_eff
        real(wp), dimension(2) :: Re
        integer                :: i, q

        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            alpha_rho(i) = q_cons_vf(i)%sf(j, k, l)
            alpha(i) = q_cons_vf(eqn_idx%E + i)%sf(j, k, l)
        end do

        ! Compute mixture properties (combined bubble and standard logic)
        rho = 0._wp
        gamma = 0._wp
        pi_inf = 0._wp

        if (bubbles_euler) then
            if (mpp_lim .and. (model_eqns == 2) .and. (num_fluids > 2)) then
                $:GPU_LOOP(parallelism='[seq]')
                do i = 1, num_fluids
                    rho = rho + alpha_rho(i)
                    gamma = gamma + alpha(i)*gammas(i)
                    pi_inf = pi_inf + alpha(i)*pi_infs(i)
                end do
            else if ((model_eqns == 2) .and. (num_fluids > 2)) then
                $:GPU_LOOP(parallelism='[seq]')
                do i = 1, num_fluids - 1
                    rho = rho + alpha_rho(i)
                    gamma = gamma + alpha(i)*gammas(i)
                    pi_inf = pi_inf + alpha(i)*pi_infs(i)
                end do
            else
                rho = alpha_rho(1)
                gamma = gammas(1)
                pi_inf = pi_infs(1)
            end if
        else
            sum_alpha = 0._wp
            if (mpp_lim) then
                $:GPU_LOOP(parallelism='[seq]')
                do i = 1, num_fluids
                    alpha_rho(i) = max(0._wp, alpha_rho(i))
                    alpha(i) = min(max(0._wp, alpha(i)), 1._wp)
                    sum_alpha = sum_alpha + alpha(i)
                end do
                alpha = alpha/max(sum_alpha, sgm_eps)
            end if

            $:GPU_LOOP(parallelism='[seq]')
            do i = 1, num_fluids
                rho = rho + alpha_rho(i)
                gamma = gamma + alpha(i)*gammas(i)
                pi_inf = pi_inf + alpha(i)*pi_infs(i)
            end do

            if (viscous) then
                $:GPU_LOOP(parallelism='[seq]')
                do i = 1, 2
                    Re(i) = dflt_real
                    if (Re_size(i) > 0) Re(i) = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do q = 1, Re_size(i)
                        Re(i) = alpha(Re_idx(i, q))/Res_pr(i, q) + Re(i)
                    end do
                    Re(i) = 1._wp/max(Re(i), sgm_eps)
                end do
            end if
        end if

        ! Compute dynamic pressure and update internal energies
        dyn_pres = 0._wp
        $:GPU_LOOP(parallelism='[seq]')
        do i = eqn_idx%mom%beg, eqn_idx%mom%end
            dyn_pres = dyn_pres + 5.e-1_wp*q_cons_vf(i)%sf(j, k, l)*q_cons_vf(i)%sf(j, k, l)/max(rho, sgm_eps)
        end do

        ! For pure stiffened/ideal-gas mixtures (no JWL fluid), recover the equilibrated pressure from the conserved total
        ! energy via the exact mixture inversion. This keeps the total energy the master variable (energy-conserving) and
        ! reproduces the legacy behavior bit-for-bit. JWL mixtures cannot use the stiffened-gas mixture inversion, so they
        ! rely on the per-phase relaxed pressure computed upstream.
        if (jwl_idx == 0) then
            pres_eff = (q_cons_vf(eqn_idx%E)%sf(j, k, l) - dyn_pres - pi_inf)/gamma
        else
            pres_eff = pres_relax
        end if

        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            call s_phase_energy_from_pressure(i, q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l), &
                                              & q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l), pres_eff, &
                                              & q_cons_vf(i + eqn_idx%int_en%beg - 1)%sf(j, k, l))
        end do

    end subroutine s_correct_internal_energies

end module m_pressure_relaxation
