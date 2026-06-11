!>
!! @file
!! @brief Contains module m_riemann_solver_hllc

!> @brief HLLC Riemann solver with contact restoration, Toro et al. Shock Waves (1994)
#:include 'case.fpp'
#:include 'macros.fpp'
#:include 'inline_riemann.fpp'

module m_riemann_solver_hllc

    use m_derived_types
    use m_global_parameters
    use m_variables_conversion
    use m_bubbles
    use m_constants, only: riemann_solver_hll, riemann_solver_hllc, riemann_solver_lax_friedrichs, model_eqns_5eq, &
        & model_eqns_6eq, model_eqns_4eq, avg_state_roe, avg_state_arithmetic, wave_speeds_direct, wave_speeds_pressure
    use m_bubbles_EE
    use m_surface_tension
    use m_chemistry
    use m_thermochem, only: gas_constant, get_mixture_molecular_weight, get_mixture_specific_heat_cv_mass, &
        & get_mixture_energy_mass, get_species_specific_heats_r, get_species_enthalpies_rt, get_mixture_specific_heat_cp_mass, &
        & molecular_weights
    use m_riemann_state

    implicit none

contains

    !> HLLC Riemann solver with contact restoration, Toro et al. Shock Waves (1994)
    subroutine s_hllc_riemann_solver(qL_prim_rsx_vf, dqL_prim_dx_vf, dqL_prim_dy_vf, &

        & dqL_prim_dz_vf, qL_prim_vf, qR_prim_rsx_vf, dqR_prim_dx_vf, dqR_prim_dy_vf, dqR_prim_dz_vf, qR_prim_vf, q_prim_vf, &
            & flux_vf, flux_src_vf, flux_gsrc_vf, norm_dir, ix, iy, iz)

        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout) :: qL_prim_rsx_vf, qR_prim_rsx_vf
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), allocatable, dimension(:), intent(inout) :: qL_prim_vf, qR_prim_vf
        type(scalar_field), allocatable, dimension(:), intent(inout) :: dqL_prim_dx_vf, dqR_prim_dx_vf, dqL_prim_dy_vf, &
             & dqR_prim_dy_vf, dqL_prim_dz_vf, dqR_prim_dz_vf

        ! Intercell fluxes
        type(scalar_field), dimension(sys_size), intent(inout) :: flux_vf, flux_src_vf, flux_gsrc_vf
        integer, intent(in)                                    :: norm_dir
        type(int_bounds_info), intent(in)                      :: ix, iy, iz

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: alpha_rho_L, alpha_rho_R
            real(wp), dimension(3) :: alpha_L, alpha_R
            real(wp), dimension(3) :: alpha_lim_L, alpha_lim_R
            real(wp), dimension(3) :: vel_L, vel_R
        #:else
            real(wp), dimension(num_fluids) :: alpha_rho_L, alpha_rho_R
            real(wp), dimension(num_fluids) :: alpha_L, alpha_R
            !> Post-limiter volume fractions (alpha_L/R retain the pre-limiter loads used downstream)
            real(wp), dimension(num_fluids) :: alpha_lim_L, alpha_lim_R
            real(wp), dimension(num_dims)   :: vel_L, vel_R
        #:endif

        real(wp) :: rho_L, rho_R
        real(wp) :: pres_L, pres_R
        real(wp) :: E_L, E_R
        real(wp) :: H_L, H_R
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(10) :: Ys_L, Ys_R, Xs_L, Xs_R, Gamma_iL, Gamma_iR, Cp_iL, Cp_iR
            real(wp), dimension(10) :: Yi_avg, Phi_avg, h_iL, h_iR, h_avg_2
        #:else
            real(wp), dimension(num_species) :: Ys_L, Ys_R, Xs_L, Xs_R, Gamma_iL, Gamma_iR, Cp_iL, Cp_iR
            real(wp), dimension(num_species) :: Yi_avg, Phi_avg, h_iL, h_iR, h_avg_2
        #:endif
        real(wp)               :: Cp_avg, Cv_avg, T_avg, c_sum_Yi_Phi, eps
        real(wp)               :: T_L, T_R
        real(wp)               :: MW_L, MW_R
        real(wp)               :: R_gas_L, R_gas_R
        real(wp)               :: Cp_L, Cp_R
        real(wp)               :: Cv_L, Cv_R
        real(wp)               :: Gamm_L, Gamm_R
        real(wp)               :: Y_L, Y_R
        real(wp)               :: gamma_L, gamma_R
        real(wp)               :: pi_inf_L, pi_inf_R
        real(wp)               :: qv_L, qv_R
        real(wp)               :: c_L, c_R
        real(wp), dimension(2) :: Re_L, Re_R
        real(wp)               :: rho_avg
        real(wp)               :: H_avg
        real(wp)               :: gamma_avg
        real(wp)               :: qv_avg
        real(wp)               :: c_avg
        real(wp)               :: s_L, s_R, s_M, s_P, s_S
        real(wp)               :: xi_L, xi_R        !< Left and right wave speeds functions
        real(wp)               :: xi_L_m1, xi_R_m1  !< xi_L/R - 1, computed without cancellation
        real(wp)               :: xi_M, xi_P
        real(wp)               :: xi_MP, xi_PP
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: R0_L, R0_R
            real(wp), dimension(3) :: V0_L, V0_R
            real(wp), dimension(3) :: P0_L, P0_R
            real(wp), dimension(3) :: pbw_L, pbw_R
        #:else
            real(wp), dimension(nb) :: R0_L, R0_R
            real(wp), dimension(nb) :: V0_L, V0_R
            real(wp), dimension(nb) :: P0_L, P0_R
            real(wp), dimension(nb) :: pbw_L, pbw_R
        #:endif

        real(wp)               :: alpha_L_sum, alpha_R_sum, nbub_L, nbub_R
        real(wp)               :: ptilde_L, ptilde_R
        real(wp)               :: PbwR3Lbar, PbwR3Rbar
        real(wp)               :: R3Lbar, R3Rbar
        real(wp)               :: R3V2Lbar, R3V2Rbar
        real(wp), dimension(6) :: tau_e_L, tau_e_R
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: xi_field_L, xi_field_R
        #:else
            real(wp), dimension(num_dims) :: xi_field_L, xi_field_R
        #:endif
        real(wp) :: G_L, G_R
        real(wp) :: vel_L_rms, vel_R_rms, vel_avg_rms
        real(wp) :: vel_L_tmp, vel_R_tmp
        real(wp) :: rho_Star, E_Star, p_Star, p_K_Star, vel_K_star
        real(wp) :: pres_SL, pres_SR, Ms_L, Ms_R
        real(wp) :: flux_ene_e
        real(wp) :: zcoef, pcorr   !< low Mach number correction
        integer  :: i, j, k, l, q  !< Generic loop iterators
        ! Populating the buffers of the left and right Riemann problem states variables, based on the choice of boundary conditions

        call s_populate_riemann_states_variables_buffers(qL_prim_rsx_vf, dqL_prim_dx_vf, dqL_prim_dy_vf, dqL_prim_dz_vf, &
            & qR_prim_rsx_vf, dqR_prim_dx_vf, dqR_prim_dy_vf, dqR_prim_dz_vf, norm_dir, ix, iy, iz)

        ! Reshaping inputted data based on dimensional splitting direction

        call s_initialize_riemann_solver(flux_src_vf, norm_dir)

        #:for NORM_DIR, XYZ, STENCIL_VAR, COORDS, X_BND, Y_BND, Z_BND in &
                    [(1, 'x', 'j', '{STENCIL_IDX}, k, l', 'is1', 'is2', 'is3'), &
                     (2, 'y', 'k', 'j, {STENCIL_IDX}, l', 'is2', 'is1', 'is3'), &
                     (3, 'z', 'l', 'j, k, {STENCIL_IDX}', 'is3', 'is2', 'is1')]
            #:set SV = STENCIL_VAR
            #:set SF = lambda offs: COORDS.format(STENCIL_IDX = SV + offs)
            if (norm_dir == ${NORM_DIR}$) then
                ! 6-EQUATION MODEL WITH HLLC HLLC star-state flux with contact wave speed s_S
                if (model_eqns == model_eqns_6eq) then
                    ! 6-equation model (model_eqns=3): separate phasic internal energies
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, vel_L, vel_R, Re_L, Re_R, alpha_L, alpha_R, &
                                        & alpha_rho_L, alpha_rho_R, Ys_L, Ys_R, Xs_L, Xs_R, Gamma_iL, Gamma_iR, Cp_iL, Cp_iR, &
                                        & Yi_avg, Phi_avg, h_iL, h_iR, h_avg_2, tau_e_L, tau_e_R, flux_ene_e, xi_field_L, &
                                        & xi_field_R, pcorr, zcoef, rho_L, rho_R, pres_L, pres_R, E_L, E_R, H_L, H_R, Cp_avg, &
                                        & Cv_avg, T_avg, eps, c_sum_Yi_Phi, T_L, T_R, Y_L, Y_R, MW_L, MW_R, R_gas_L, R_gas_R, &
                                        & Cp_L, Cp_R, Cv_L, Cv_R, Gamm_L, Gamm_R, gamma_L, gamma_R, pi_inf_L, pi_inf_R, qv_L, &
                                        & qv_R, qv_avg, c_L, c_R, G_L, G_R, rho_avg, H_avg, c_avg, gamma_avg, ptilde_L, ptilde_R, &
                                        & vel_L_rms, vel_R_rms, vel_avg_rms, vel_L_tmp, vel_R_tmp, Ms_L, Ms_R, pres_SL, pres_SR, &
                                        & alpha_L_sum, alpha_R_sum, rho_Star, E_Star, p_Star, p_K_Star, vel_K_star, s_L, s_R, &
                                        & s_M, s_P, s_S, xi_M, xi_P, xi_L, xi_R, xi_L_m1, xi_R_m1, xi_MP, xi_PP]')
                    do l = ${Z_BND}$%beg, ${Z_BND}$%end
                        do k = ${Y_BND}$%beg, ${Y_BND}$%end
                            do j = ${X_BND}$%beg, ${X_BND}$%end
                                vel_L_rms = 0._wp; vel_R_rms = 0._wp
                                rho_L = 0._wp; rho_R = 0._wp
                                gamma_L = 0._wp; gamma_R = 0._wp
                                pi_inf_L = 0._wp; pi_inf_R = 0._wp
                                qv_L = 0._wp; qv_R = 0._wp
                                alpha_L_sum = 0._wp; alpha_R_sum = 0._wp

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel_L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%cont%end + i)
                                    vel_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%cont%end + i)
                                    vel_L_rms = vel_L_rms + vel_L(i)**2._wp
                                    vel_R_rms = vel_R_rms + vel_R(i)**2._wp
                                end do

                                pres_L = qL_prim_rsx_vf(${SF('')}$, eqn_idx%E)
                                pres_R = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E)

                                rho_L = 0._wp
                                gamma_L = 0._wp
                                pi_inf_L = 0._wp
                                qv_L = 0._wp

                                rho_R = 0._wp
                                gamma_R = 0._wp
                                pi_inf_R = 0._wp
                                qv_R = 0._wp

                                alpha_L_sum = 0._wp
                                alpha_R_sum = 0._wp

                                if (mpp_lim) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        qL_prim_rsx_vf(${SF('')}$, i) = max(0._wp, qL_prim_rsx_vf(${SF('')}$, i))
                                        qL_prim_rsx_vf(${SF('')}$, eqn_idx%E + i) = min(max(0._wp, qL_prim_rsx_vf(${SF('')}$, &
                                                       & eqn_idx%E + i)), 1._wp)
                                        alpha_L_sum = alpha_L_sum + qL_prim_rsx_vf(${SF('')}$, eqn_idx%E + i)
                                    end do

                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        qR_prim_rsx_vf(${SF(' + 1')}$, i) = max(0._wp, qR_prim_rsx_vf(${SF(' + 1')}$, i))
                                        qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E + i) = min(max(0._wp, &
                                                       & qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E + i)), 1._wp)
                                        alpha_R_sum = alpha_R_sum + qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E + i)
                                    end do

                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        qL_prim_rsx_vf(${SF('')}$, eqn_idx%E + i) = qL_prim_rsx_vf(${SF('')}$, &
                                                       & eqn_idx%E + i)/max(alpha_L_sum, sgm_eps)
                                        qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E + i) = qR_prim_rsx_vf(${SF(' + 1')}$, &
                                                       & eqn_idx%E + i)/max(alpha_R_sum, sgm_eps)
                                    end do
                                end if

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    alpha_rho_L(i) = qL_prim_rsx_vf(${SF('')}$, i)
                                    alpha_rho_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, i)
                                    alpha_L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%adv%beg + i - 1)
                                    alpha_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%adv%beg + i - 1)
                                end do

                                call s_accumulate_mixture_properties(num_fluids, alpha_rho_L, alpha_L, rho_L, gamma_L, pi_inf_L, &
                                                                     & qv_L)
                                call s_accumulate_mixture_properties(num_fluids, alpha_rho_R, alpha_R, rho_R, gamma_R, pi_inf_R, &
                                                                     & qv_R)

                                if (viscous) then
                                    call s_compute_interface_reynolds(alpha_L, Re_L)
                                    call s_compute_interface_reynolds(alpha_R, Re_R)
                                end if

                                E_L = gamma_L*pres_L + pi_inf_L + 5.e-1_wp*rho_L*vel_L_rms + qv_L
                                E_R = gamma_R*pres_R + pi_inf_R + 5.e-1_wp*rho_R*vel_R_rms + qv_R

                                ! ENERGY ADJUSTMENTS FOR HYPOELASTIC ENERGY
                                if (hypoelasticity) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, eqn_idx%stress%end - eqn_idx%stress%beg + 1
                                        tau_e_L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%stress%beg - 1 + i)
                                        tau_e_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%stress%beg - 1 + i)
                                    end do
                                    G_L = 0._wp; G_R = 0._wp
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        G_L = G_L + alpha_L(i)*Gs_rs(i)
                                        G_R = G_R + alpha_R(i)*Gs_rs(i)
                                    end do
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, eqn_idx%stress%end - eqn_idx%stress%beg + 1
                                        ! Elastic contribution to energy if G large enough
                                        if ((G_L > verysmall) .and. (G_R > verysmall)) then
                                            E_L = E_L + (tau_e_L(i)*tau_e_L(i))/(4._wp*G_L)
                                            E_R = E_R + (tau_e_R(i)*tau_e_R(i))/(4._wp*G_R)
                                            ! Additional terms in 2D and 3D
                                            if ((i == 2) .or. (i == 4) .or. (i == 5)) then
                                                E_L = E_L + (tau_e_L(i)*tau_e_L(i))/(4._wp*G_L)
                                                E_R = E_R + (tau_e_R(i)*tau_e_R(i))/(4._wp*G_R)
                                            end if
                                        end if
                                    end do
                                end if

                                ! Hyperelastic stress contribution: strain energy added to total energy
                                if (hyperelasticity) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        xi_field_L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%xi%beg - 1 + i)
                                        xi_field_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%xi%beg - 1 + i)
                                    end do
                                    G_L = 0._wp; G_R = 0._wp
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        ! Mixture left and right shear modulus
                                        G_L = G_L + alpha_L(i)*Gs_rs(i)
                                        G_R = G_R + alpha_R(i)*Gs_rs(i)
                                    end do
                                    ! Elastic contribution to energy if G large enough
                                    if (G_L > verysmall .and. G_R > verysmall) then
                                        E_L = E_L + G_L*qL_prim_rsx_vf(${SF('')}$, eqn_idx%xi%end + 1)
                                        E_R = E_R + G_R*qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%xi%end + 1)
                                    end if
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, b_size - 1
                                        tau_e_L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%stress%beg - 1 + i)
                                        tau_e_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%stress%beg - 1 + i)
                                    end do
                                end if

                                H_L = (E_L + pres_L)/rho_L
                                H_R = (E_R + pres_R)/rho_R

                                @:compute_average_state()

                                call s_compute_speed_of_sound(pres_L, rho_L, gamma_L, pi_inf_L, H_L, alpha_L, vel_L_rms, 0._wp, &
                                                              & c_L, qv_L)

                                call s_compute_speed_of_sound(pres_R, rho_R, gamma_R, pi_inf_R, H_R, alpha_R, vel_R_rms, 0._wp, &
                                                              & c_R, qv_R)

                                !> The computation of c_avg does not require all the variables, and therefore the non '_avg'
                                ! variables are placeholders to call the subroutine.
                                call s_compute_speed_of_sound(pres_R, rho_avg, gamma_avg, pi_inf_R, H_avg, alpha_R, vel_avg_rms, &
                                                              & 0._wp, c_avg, qv_avg)

                                if (viscous) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, 2
                                        Re_avg_rsx_vf(${SF('')}$, i) = 2._wp/(1._wp/Re_L(i) + 1._wp/Re_R(i))
                                    end do
                                end if

                                ! Low Mach correction
                                if (low_Mach == 2) then
                                    @:compute_low_Mach_correction()
                                end if

                                ! COMPUTING THE DIRECT WAVE SPEEDS
                                if (wave_speeds == wave_speeds_direct) then
                                    if (elasticity) then
                                        ! Elastic wave speed, Rodriguez et al. JCP (2019)
                                        s_L = min(vel_L(dir_idx(1)) - sqrt(c_L*c_L + (((4._wp*G_L)/3._wp) + tau_e_L(dir_idx_tau(1) &
                                                  & ))/rho_L), &
                                                  & vel_R(dir_idx(1)) - sqrt(c_R*c_R + (((4._wp*G_R)/3._wp) &
                                                  & + tau_e_R(dir_idx_tau(1)))/rho_R))
                                        s_R = max(vel_R(dir_idx(1)) + sqrt(c_R*c_R + (((4._wp*G_R)/3._wp) + tau_e_R(dir_idx_tau(1) &
                                                  & ))/rho_R), &
                                                  & vel_L(dir_idx(1)) + sqrt(c_L*c_L + (((4._wp*G_L)/3._wp) &
                                                  & + tau_e_L(dir_idx_tau(1)))/rho_L))
                                        s_S = (pres_R - tau_e_R(dir_idx_tau(1)) - pres_L + tau_e_L(dir_idx_tau(1)) &
                                               & + rho_L*vel_L(dir_idx(1))*(s_L - vel_L(dir_idx(1))) - rho_R*vel_R(dir_idx(1)) &
                                               & *(s_R - vel_R(dir_idx(1))))/(rho_L*(s_L - vel_L(dir_idx(1))) - rho_R*(s_R &
                                               & - vel_R(dir_idx(1))))
                                    else
                                        s_L = min(vel_L(dir_idx(1)) - c_L, vel_R(dir_idx(1)) - c_R)
                                        s_R = max(vel_R(dir_idx(1)) + c_R, vel_L(dir_idx(1)) + c_L)
                                        s_S = (pres_R - pres_L + rho_L*vel_L(dir_idx(1))*(s_L - vel_L(dir_idx(1))) &
                                               & - rho_R*vel_R(dir_idx(1))*(s_R - vel_R(dir_idx(1))))/(rho_L*(s_L &
                                               & - vel_L(dir_idx(1))) - rho_R*(s_R - vel_R(dir_idx(1))))
                                    end if
                                else if (wave_speeds == wave_speeds_pressure) then
                                    pres_SL = 5.e-1_wp*(pres_L + pres_R + rho_avg*c_avg*(vel_L(dir_idx(1)) - vel_R(dir_idx(1))))

                                    pres_SR = pres_SL

                                    ! Low Mach correction: Thornber et al. JCP (2008)
                                    Ms_L = max(1._wp, &
                                               & sqrt(1._wp + ((5.e-1_wp + gamma_L)/(1._wp + gamma_L))*(pres_SL/pres_L - 1._wp) &
                                               & *pres_L/((pres_L + pi_inf_L/(1._wp + gamma_L)))))
                                    Ms_R = max(1._wp, &
                                               & sqrt(1._wp + ((5.e-1_wp + gamma_R)/(1._wp + gamma_R))*(pres_SR/pres_R - 1._wp) &
                                               & *pres_R/((pres_R + pi_inf_R/(1._wp + gamma_R)))))

                                    s_L = vel_L(dir_idx(1)) - c_L*Ms_L
                                    s_R = vel_R(dir_idx(1)) + c_R*Ms_R

                                    s_S = 5.e-1_wp*((vel_L(dir_idx(1)) + vel_R(dir_idx(1))) + (pres_L - pres_R)/(rho_avg*c_avg))
                                end if

                                ! follows Einfeldt et al. s_M/P = min/max(0.,s_L/R)
                                s_M = min(0._wp, s_L); s_P = max(0._wp, s_R)

                                ! goes with q_star_L/R = xi_L/R * (variable) xi_L/R = ( ( s_L/R - u_L/R )/(s_L/R - s_star) )
                                xi_L = (s_L - vel_L(dir_idx(1)))/min(s_L - s_S, -sgm_eps)
                                xi_R = (s_R - vel_R(dir_idx(1)))/max(s_R - s_S, sgm_eps)
                                xi_L_m1 = (s_S - vel_L(dir_idx(1)))/min(s_L - s_S, -sgm_eps)
                                xi_R_m1 = (s_S - vel_R(dir_idx(1)))/max(s_R - s_S, sgm_eps)

                                ! goes with numerical star velocity in x/y/z directions xi_P/M = 0.5 +/m sgn(0.5,s_star)
                                xi_M = (5.e-1_wp + sign(0.5_wp, s_S))
                                xi_P = (5.e-1_wp - sign(0.5_wp, s_S))

                                ! goes with the numerical velocity in x/y/z directions xi_P/M (pressure) = min/max(0. sgn(1,sL/sR))
                                xi_MP = -min(0._wp, sign(1._wp, s_L))
                                xi_PP = max(0._wp, sign(1._wp, s_R))

                                E_star = xi_M*(E_L + xi_MP*(xi_L*(E_L + (s_S - vel_L(dir_idx(1)))*(rho_L*s_S + pres_L/(s_L &
                                               & - vel_L(dir_idx(1))))) - E_L)) + xi_P*(E_R + xi_PP*(xi_R*(E_R + (s_S &
                                               & - vel_R(dir_idx(1)))*(rho_R*s_S + pres_R/(s_R - vel_R(dir_idx(1))))) - E_R))
                                p_Star = xi_M*(pres_L + xi_MP*(rho_L*(s_L - vel_L(dir_idx(1)))*(s_S - vel_L(dir_idx(1))))) &
                                               & + xi_P*(pres_R + xi_PP*(rho_R*(s_R - vel_R(dir_idx(1)))*(s_S - vel_R(dir_idx(1)))))

                                rho_Star = xi_M*(rho_L*(xi_MP*xi_L + 1._wp - xi_MP)) + xi_P*(rho_R*(xi_PP*xi_R + 1._wp - xi_PP))

                                vel_K_Star = vel_L(dir_idx(1))*(1._wp - xi_MP) + xi_MP*vel_R(dir_idx(1)) + xi_MP*xi_PP*(s_S &
                                                   & - vel_R(dir_idx(1)))

                                ! Low Mach correction
                                if (low_Mach == 1) then
                                    @:compute_low_Mach_correction()
                                else
                                    pcorr = 0._wp
                                end if

                                ! COMPUTING FLUXES MASS FLUX.
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, eqn_idx%cont%end
                                    flux_rsx_vf(${SF('')}$, i) = xi_M*qL_prim_rsx_vf(${SF('')}$, &
                                                & i)*(vel_L(dir_idx(1)) + s_M*xi_L_m1) + xi_P*qR_prim_rsx_vf(${SF(' + 1')}$, &
                                                & i)*(vel_R(dir_idx(1)) + s_P*xi_R_m1)
                                end do

                                ! MOMENTUM FLUX. f = \rho u u - \sigma, q = \rho u, q_star = \xi * \rho*(s_star, v, w)
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    flux_rsx_vf(${SF('')}$, &
                                                & eqn_idx%cont%end + dir_idx(i)) = rho_Star*vel_K_Star*(dir_flg(dir_idx(i)) &
                                                & *vel_K_Star + (1._wp - dir_flg(dir_idx(i)))*(xi_M*vel_L(dir_idx(i)) &
                                                & + xi_P*vel_R(dir_idx(i)))) + dir_flg(dir_idx(i))*p_Star + (s_M/s_L)*(s_P/s_R) &
                                                & *dir_flg(dir_idx(i))*pcorr
                                end do

                                ! ENERGY FLUX. f = u*(E-\sigma), q = E, q_star = \xi*E+(s-u)(\rho s_star - \sigma/(s-u))
                                flux_rsx_vf(${SF('')}$, eqn_idx%E) = (E_star + p_Star)*vel_K_Star + (s_M/s_L)*(s_P/s_R)*pcorr*s_S

                                ! ELASTICITY. Elastic shear stress additions for the momentum and energy flux
                                if (elasticity) then
                                    flux_ene_e = 0._wp
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        ! MOMENTUM ELASTIC FLUX.
                                        flux_rsx_vf(${SF('')}$, eqn_idx%cont%end + dir_idx(i)) = flux_rsx_vf(${SF('')}$, &
                                                    & eqn_idx%cont%end + dir_idx(i)) - xi_M*tau_e_L(dir_idx_tau(i)) &
                                                    & - xi_P*tau_e_R(dir_idx_tau(i))
                                        ! ENERGY ELASTIC FLUX.
                                        flux_ene_e = flux_ene_e - xi_M*(vel_L(dir_idx(i))*tau_e_L(dir_idx_tau(i)) &
                                                                        & + s_M*(xi_L*((s_S - vel_L(i))*(tau_e_L(dir_idx_tau(i)) &
                                                                        & /(s_L - vel_L(i)))))) - xi_P*(vel_R(dir_idx(i)) &
                                                                        & *tau_e_R(dir_idx_tau(i)) + s_P*(xi_R*((s_S - vel_R(i)) &
                                                                        & *(tau_e_R(dir_idx_tau(i))/(s_R - vel_R(i))))))
                                    end do
                                    flux_rsx_vf(${SF('')}$, eqn_idx%E) = flux_rsx_vf(${SF('')}$, eqn_idx%E) + flux_ene_e
                                end if

                                ! VOLUME FRACTION FLUX.
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = eqn_idx%adv%beg, eqn_idx%adv%end
                                    flux_rsx_vf(${SF('')}$, i) = xi_M*qL_prim_rsx_vf(${SF('')}$, &
                                                & i)*s_S + xi_P*qR_prim_rsx_vf(${SF(' + 1')}$, i)*s_S
                                end do

                                ! Advection velocity source: interface velocity for volume fraction transport
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel_src_rsx_vf(${SF('')}$, &
                                                   & dir_idx(i)) = xi_M*(vel_L(dir_idx(i)) + dir_flg(dir_idx(i)) &
                                                   & *(s_S*(xi_MP*xi_L_m1 + 1) - vel_L(dir_idx(i)))) + xi_P*(vel_R(dir_idx(i)) &
                                                   & + dir_flg(dir_idx(i))*(s_S*(xi_PP*xi_R_m1 + 1) - vel_R(dir_idx(i))))
                                end do

                                ! INTERNAL ENERGIES ADVECTION FLUX. K-th pressure and velocity in preparation for the internal
                                ! energy flux
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    p_K_Star = xi_M*(xi_MP*((pres_L + pi_infs(i)/(1._wp + gammas(i)))*xi_L**(1._wp/gammas(i) &
                                                     & + 1._wp) - pi_infs(i)/(1._wp + gammas(i)) - pres_L) + pres_L) &
                                                     & + xi_P*(xi_PP*((pres_R + pi_infs(i)/(1._wp + gammas(i))) &
                                                     & *xi_R**(1._wp/gammas(i) + 1._wp) - pi_infs(i)/(1._wp + gammas(i)) - pres_R) &
                                                     & + pres_R)

                                    flux_rsx_vf(${SF('')}$, i + eqn_idx%int_en%beg - 1) = ((xi_M*qL_prim_rsx_vf(${SF('')}$, &
                                                & i + eqn_idx%adv%beg - 1) + xi_P*qR_prim_rsx_vf(${SF(' + 1')}$, &
                                                & i + eqn_idx%adv%beg - 1))*(gammas(i)*p_K_Star + pi_infs(i)) &
                                                & + (xi_M*qL_prim_rsx_vf(${SF('')}$, &
                                                & i + eqn_idx%cont%beg - 1) + xi_P*qR_prim_rsx_vf(${SF(' + 1')}$, &
                                                & i + eqn_idx%cont%beg - 1))*qvs(i))*vel_K_Star + (s_M/s_L)*(s_P/s_R) &
                                                & *pcorr*s_S*(xi_M*qL_prim_rsx_vf(${SF('')}$, &
                                                & i + eqn_idx%adv%beg - 1) + xi_P*qR_prim_rsx_vf(${SF(' + 1')}$, &
                                                & i + eqn_idx%adv%beg - 1))
                                end do

                                flux_src_rsx_vf(${SF('')}$, eqn_idx%adv%beg) = vel_src_rsx_vf(${SF('')}$, dir_idx(1))

                                ! HYPOELASTIC STRESS EVOLUTION FLUX.
                                if (hypoelasticity) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, eqn_idx%stress%end - eqn_idx%stress%beg + 1
                                        flux_rsx_vf(${SF('')}$, &
                                                    & eqn_idx%stress%beg - 1 + i) = xi_M*(s_S/(s_L - s_S))*(s_L*rho_L*tau_e_L(i) &
                                                    & - rho_L*vel_L(dir_idx(1))*tau_e_L(i)) + xi_P*(s_S/(s_R - s_S)) &
                                                    & *(s_R*rho_R*tau_e_R(i) - rho_R*vel_R(dir_idx(1))*tau_e_R(i))
                                    end do
                                end if

                                ! Hyperelastic reference map flux for material deformation tracking
                                if (hyperelasticity) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        flux_rsx_vf(${SF('')}$, &
                                                    & eqn_idx%xi%beg - 1 + i) = xi_M*(s_S/(s_L - s_S))*(s_L*rho_L*xi_field_L(i) &
                                                    & - rho_L*vel_L(dir_idx(1))*xi_field_L(i)) + xi_P*(s_S/(s_R - s_S)) &
                                                    & *(s_R*rho_R*xi_field_R(i) - rho_R*vel_R(dir_idx(1))*xi_field_R(i))
                                    end do
                                end if

                                ! COLOR FUNCTION FLUX
                                if (surface_tension) then
                                    flux_rsx_vf(${SF('')}$, eqn_idx%c) = (xi_M*qL_prim_rsx_vf(${SF('')}$, &
                                                & eqn_idx%c) + xi_P*qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%c))*s_S
                                end if

                                ! Geometrical source flux for cylindrical coordinates
                                #:if (NORM_DIR == 2)
                                    if (cyl_coord) then
                                        ! Substituting the advective flux into the inviscid geometrical source flux
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, eqn_idx%E
                                            flux_gsrc_rsx_vf(${SF('')}$, i) = flux_rsx_vf(${SF('')}$, i)
                                        end do
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = eqn_idx%int_en%beg, eqn_idx%int_en%end
                                            flux_gsrc_rsx_vf(${SF('')}$, i) = flux_rsx_vf(${SF('')}$, i)
                                        end do
                                        ! Recalculating the radial momentum geometric source flux
                                        flux_gsrc_rsx_vf(${SF('')}$, &
                                                         & eqn_idx%mom%beg - 1 + dir_idx(1)) = flux_gsrc_rsx_vf(${SF('')}$, &
                                                         & eqn_idx%mom%beg - 1 + dir_idx(1)) - p_Star
                                        ! Geometrical source of the void fraction(s) is zero
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = eqn_idx%adv%beg, eqn_idx%adv%end
                                            flux_gsrc_rsx_vf(${SF('')}$, i) = 0._wp
                                        end do
                                    end if
                                #:endif
                                #:if (NORM_DIR == 3)
                                    if (grid_geometry == 3) then
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, sys_size
                                            flux_gsrc_rsx_vf(${SF('')}$, i) = 0._wp
                                        end do
                                        flux_gsrc_rsx_vf(${SF('')}$, &
                                                         & eqn_idx%mom%beg - 1 + dir_idx(1)) = flux_gsrc_rsx_vf(${SF('')}$, &
                                                         & eqn_idx%mom%beg - 1 + dir_idx(1)) - p_Star

                                        flux_gsrc_rsx_vf(${SF('')}$, eqn_idx%mom%end) = flux_rsx_vf(${SF('')}$, eqn_idx%mom%beg + 1)
                                    end if
                                #:endif
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                else if (model_eqns == model_eqns_4eq) then
                    ! 4-equation model (model_eqns=4): single pressure, velocity equilibrium
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[i, q, alpha_rho_L, alpha_rho_R, vel_L, vel_R, alpha_L, alpha_R, &
                                        & nbub_L, nbub_R, rho_L, rho_R, pres_L, pres_R, E_L, E_R, H_L, H_R, Cp_avg, Cv_avg, &
                                        & T_avg, eps, c_sum_Yi_Phi, T_L, T_R, Y_L, Y_R, MW_L, MW_R, R_gas_L, R_gas_R, Cp_L, Cp_R, &
                                        & Gamm_L, Gamm_R, gamma_L, gamma_R, pi_inf_L, pi_inf_R, qv_L, qv_R, qv_avg, c_L, c_R, &
                                        & G_L, G_R, rho_avg, H_avg, c_avg, gamma_avg, ptilde_L, ptilde_R, vel_L_rms, vel_R_rms, &
                                        & vel_avg_rms, vel_L_tmp, vel_R_tmp, Ms_L, Ms_R, pres_SL, pres_SR, alpha_L_sum, &
                                        & alpha_R_sum, rho_Star, E_Star, p_Star, p_K_Star, vel_K_star, s_L, s_R, s_M, s_P, s_S, &
                                        & xi_M, xi_P, xi_L, xi_R, xi_L_m1, xi_R_m1, xi_MP, xi_PP, Ys_L, Ys_R, Cp_iL, Cp_iR, Xs_L, &
                                        & Xs_R, Gamma_iL, Gamma_iR, Yi_avg, Phi_avg, h_iL, h_iR, h_avg_2]')
                    do l = ${Z_BND}$%beg, ${Z_BND}$%end
                        do k = ${Y_BND}$%beg, ${Y_BND}$%end
                            do j = ${X_BND}$%beg, ${X_BND}$%end
                                vel_L_rms = 0._wp; vel_R_rms = 0._wp
                                rho_L = 0._wp; rho_R = 0._wp
                                gamma_L = 0._wp; gamma_R = 0._wp
                                pi_inf_L = 0._wp; pi_inf_R = 0._wp
                                qv_L = 0._wp; qv_R = 0._wp

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, eqn_idx%cont%end
                                    alpha_rho_L(i) = qL_prim_rsx_vf(${SF('')}$, i)
                                    alpha_rho_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, i)
                                end do

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel_L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%cont%end + i)
                                    vel_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%cont%end + i)
                                    vel_L_rms = vel_L_rms + vel_L(i)**2._wp
                                    vel_R_rms = vel_R_rms + vel_R(i)**2._wp
                                end do

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    alpha_L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%E + i)
                                    alpha_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E + i)
                                end do
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    alpha_L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%E + i)
                                    alpha_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E + i)
                                end do

                                call s_accumulate_mixture_properties(num_fluids, alpha_rho_L, alpha_L, rho_L, gamma_L, pi_inf_L, &
                                                                     & qv_L)
                                call s_accumulate_mixture_properties(num_fluids, alpha_rho_R, alpha_R, rho_R, gamma_R, pi_inf_R, &
                                                                     & qv_R)

                                pres_L = qL_prim_rsx_vf(${SF('')}$, eqn_idx%E)
                                pres_R = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E)

                                E_L = gamma_L*pres_L + pi_inf_L + 5.e-1_wp*rho_L*vel_L_rms + qv_L
                                E_R = gamma_R*pres_R + pi_inf_R + 5.e-1_wp*rho_R*vel_R_rms + qv_R

                                H_L = (E_L + pres_L)/rho_L
                                H_R = (E_R + pres_R)/rho_R

                                @:compute_average_state()

                                call s_compute_speed_of_sound(pres_L, rho_L, gamma_L, pi_inf_L, H_L, alpha_L, vel_L_rms, 0._wp, &
                                                              & c_L, qv_L)

                                call s_compute_speed_of_sound(pres_R, rho_R, gamma_R, pi_inf_R, H_R, alpha_R, vel_R_rms, 0._wp, &
                                                              & c_R, qv_R)

                                !> The computation of c_avg does not require all the variables, and therefore the non '_avg'
                                ! variables are placeholders to call the subroutine.

                                call s_compute_speed_of_sound(pres_R, rho_avg, gamma_avg, pi_inf_R, H_avg, alpha_R, vel_avg_rms, &
                                                              & 0._wp, c_avg, qv_avg)

                                if (wave_speeds == wave_speeds_direct) then
                                    s_L = min(vel_L(dir_idx(1)) - c_L, vel_R(dir_idx(1)) - c_R)
                                    s_R = max(vel_R(dir_idx(1)) + c_R, vel_L(dir_idx(1)) + c_L)

                                    s_S = (pres_R - pres_L + rho_L*vel_L(dir_idx(1))*(s_L - vel_L(dir_idx(1))) &
                                           & - rho_R*vel_R(dir_idx(1))*(s_R - vel_R(dir_idx(1))))/(rho_L*(s_L - vel_L(dir_idx(1))) &
                                           & - rho_R*(s_R - vel_R(dir_idx(1))))
                                else if (wave_speeds == wave_speeds_pressure) then
                                    pres_SL = 5.e-1_wp*(pres_L + pres_R + rho_avg*c_avg*(vel_L(dir_idx(1)) - vel_R(dir_idx(1))))

                                    pres_SR = pres_SL

                                    ! Low Mach correction: Thornber et al. JCP (2008)
                                    Ms_L = max(1._wp, &
                                               & sqrt(1._wp + ((5.e-1_wp + gamma_L)/(1._wp + gamma_L))*(pres_SL/pres_L - 1._wp) &
                                               & *pres_L/((pres_L + pi_inf_L/(1._wp + gamma_L)))))
                                    Ms_R = max(1._wp, &
                                               & sqrt(1._wp + ((5.e-1_wp + gamma_R)/(1._wp + gamma_R))*(pres_SR/pres_R - 1._wp) &
                                               & *pres_R/((pres_R + pi_inf_R/(1._wp + gamma_R)))))

                                    s_L = vel_L(dir_idx(1)) - c_L*Ms_L
                                    s_R = vel_R(dir_idx(1)) + c_R*Ms_R

                                    s_S = 5.e-1_wp*((vel_L(dir_idx(1)) + vel_R(dir_idx(1))) + (pres_L - pres_R)/(rho_avg*c_avg))
                                end if

                                ! follows Einfeldt et al. s_M/P = min/max(0.,s_L/R)
                                s_M = min(0._wp, s_L); s_P = max(0._wp, s_R)

                                ! goes with q_star_L/R = xi_L/R * (variable) xi_L/R = ( ( s_L/R - u_L/R )/(s_L/R - s_star) )
                                xi_L = (s_L - vel_L(dir_idx(1)))/min(s_L - s_S, -sgm_eps)
                                xi_R = (s_R - vel_R(dir_idx(1)))/max(s_R - s_S, sgm_eps)
                                xi_L_m1 = (s_S - vel_L(dir_idx(1)))/min(s_L - s_S, -sgm_eps)
                                xi_R_m1 = (s_S - vel_R(dir_idx(1)))/max(s_R - s_S, sgm_eps)

                                ! goes with numerical velocity in x/y/z directions xi_P/M = 0.5 +/m sgn(0.5,s_star)
                                xi_M = (5.e-1_wp + sign(5.e-1_wp, s_S))
                                xi_P = (5.e-1_wp - sign(5.e-1_wp, s_S))

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, eqn_idx%cont%end
                                    flux_rsx_vf(${SF('')}$, &
                                                & i) = xi_M*alpha_rho_L(i)*(vel_L(dir_idx(1)) + s_M*xi_L_m1) + xi_P*alpha_rho_R(i) &
                                                & *(vel_R(dir_idx(1)) + s_P*xi_R_m1)
                                end do

                                ! Momentum flux. f = \rho u u + p I, q = \rho u, q_star = \xi * \rho*(s_star, v, w)
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    flux_rsx_vf(${SF('')}$, &
                                                & eqn_idx%cont%end + dir_idx(i)) = xi_M*(rho_L*(vel_L(dir_idx(1))*vel_L(dir_idx(i) &
                                                & ) + s_M*(xi_L*(dir_flg(dir_idx(i))*s_S + (1._wp - dir_flg(dir_idx(i))) &
                                                & *vel_L(dir_idx(i))) - vel_L(dir_idx(i)))) + dir_flg(dir_idx(i))*pres_L) &
                                                & + xi_P*(rho_R*(vel_R(dir_idx(1))*vel_R(dir_idx(i)) &
                                                & + s_P*(xi_R*(dir_flg(dir_idx(i))*s_S + (1._wp - dir_flg(dir_idx(i))) &
                                                & *vel_R(dir_idx(i))) - vel_R(dir_idx(i)))) + dir_flg(dir_idx(i))*pres_R)
                                end do

                                if (bubbles_euler) then
                                    ! Put p_tilde in
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        flux_rsx_vf(${SF('')}$, eqn_idx%cont%end + dir_idx(i)) = flux_rsx_vf(${SF('')}$, &
                                                    & eqn_idx%cont%end + dir_idx(i)) + xi_M*(dir_flg(dir_idx(i))*(-1._wp*ptilde_L) &
                                                    & ) + xi_P*(dir_flg(dir_idx(i))*(-1._wp*ptilde_R))
                                    end do
                                end if

                                flux_rsx_vf(${SF('')}$, eqn_idx%E) = 0._wp

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = eqn_idx%alf, eqn_idx%alf  ! only advect the void fraction
                                    flux_rsx_vf(${SF('')}$, i) = xi_M*qL_prim_rsx_vf(${SF('')}$, &
                                                & i)*(vel_L(dir_idx(1)) + s_M*xi_L_m1) + xi_P*qR_prim_rsx_vf(${SF(' + 1')}$, &
                                                & i)*(vel_R(dir_idx(1)) + s_P*xi_R_m1)
                                end do

                                ! Advection velocity source: interface velocity for volume fraction transport
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel_src_rsx_vf(${SF('')}$, dir_idx(i)) = 0._wp
                                    ! IF ( (model_eqns == 4) .or. (num_fluids==1) ) vel_src_rs_vf(dir_idx(i))%sf(j,k,l) = 0._wp
                                end do

                                flux_src_rsx_vf(${SF('')}$, eqn_idx%adv%beg) = vel_src_rsx_vf(${SF('')}$, dir_idx(1))

                                ! Add advection flux for bubble variables
                                if (bubbles_euler) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = eqn_idx%bub%beg, eqn_idx%bub%end
                                        flux_rsx_vf(${SF('')}$, i) = xi_M*nbub_L*qL_prim_rsx_vf(${SF('')}$, &
                                                    & i)*(vel_L(dir_idx(1)) + s_M*xi_L_m1) &
                                                    & + xi_P*nbub_R*qR_prim_rsx_vf(${SF(' + 1')}$, &
                                                    & i)*(vel_R(dir_idx(1)) + s_P*xi_R_m1)
                                    end do
                                end if

                                ! Geometrical source flux for cylindrical coordinates

                                #:if (NORM_DIR == 2)
                                    if (cyl_coord) then
                                        ! Substituting the advective flux into the inviscid geometrical source flux
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, eqn_idx%E
                                            flux_gsrc_rsx_vf(${SF('')}$, i) = flux_rsx_vf(${SF('')}$, i)
                                        end do
                                        ! Recalculating the radial momentum geometric source flux
                                        flux_gsrc_rsx_vf(${SF('')}$, &
                                                         & eqn_idx%cont%end + dir_idx(1)) &
                                                         & = f_compute_hllc_star_momentum_flux(rho_L, rho_R, vel_L(dir_idx(1)), &
                                                         & vel_R(dir_idx(1)), s_M, s_P, s_S, xi_L, xi_R, xi_M, xi_P, &
                                                         & dir_flg(dir_idx(1)))
                                        ! Geometrical source of the void fraction(s) is zero
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = eqn_idx%adv%beg, eqn_idx%adv%end
                                            flux_gsrc_rsx_vf(${SF('')}$, i) = 0._wp
                                        end do
                                    end if
                                #:endif
                                #:if (NORM_DIR == 3)
                                    if (grid_geometry == 3) then
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, sys_size
                                            flux_gsrc_rsx_vf(${SF('')}$, i) = 0._wp
                                        end do
                                        flux_gsrc_rsx_vf(${SF('')}$, &
                                                         & eqn_idx%mom%beg + 1) = -f_compute_hllc_star_momentum_flux(rho_L, &
                                                         & rho_R, vel_L(dir_idx(1)), vel_R(dir_idx(1)), s_M, s_P, s_S, xi_L, &
                                                         & xi_R, xi_M, xi_P, dir_flg(dir_idx(1)))
                                        flux_gsrc_rsx_vf(${SF('')}$, eqn_idx%mom%end) = flux_rsx_vf(${SF('')}$, eqn_idx%mom%beg + 1)
                                    end if
                                #:endif
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                else if (model_eqns == model_eqns_5eq .and. bubbles_euler) then
                    ! 5-equation model with Euler-Euler bubble dynamics
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[i, q, R0_L, R0_R, V0_L, V0_R, P0_L, P0_R, pbw_L, pbw_R, vel_L, &
                                        & vel_R, rho_avg, alpha_L, alpha_R, alpha_rho_L, alpha_rho_R, h_avg, gamma_avg, Re_L, &
                                        & Re_R, pcorr, zcoef, rho_L, rho_R, pres_L, pres_R, E_L, E_R, H_L, H_R, gamma_L, gamma_R, &
                                        & pi_inf_L, pi_inf_R, qv_L, qv_R, qv_avg, c_L, c_R, c_avg, vel_L_rms, vel_R_rms, &
                                        & vel_avg_rms, vel_L_tmp, vel_R_tmp, Ms_L, Ms_R, pres_SL, pres_SR, alpha_L_sum, &
                                        & alpha_R_sum, s_L, s_R, s_M, s_P, s_S, xi_M, xi_P, xi_L, xi_R, xi_L_m1, xi_R_m1, xi_MP, &
                                        & xi_PP, nbub_L, nbub_R, PbwR3Lbar, PbwR3Rbar, R3Lbar, R3Rbar, R3V2Lbar, R3V2Rbar, Ys_L, &
                                        & Ys_R, Cp_iL, Cp_iR, Xs_L, Xs_R, Gamma_iL, Gamma_iR, Yi_avg, Phi_avg, h_iL, h_iR, h_avg_2]')
                    do l = ${Z_BND}$%beg, ${Z_BND}$%end
                        do k = ${Y_BND}$%beg, ${Y_BND}$%end
                            do j = ${X_BND}$%beg, ${X_BND}$%end
                                vel_L_rms = 0._wp; vel_R_rms = 0._wp
                                rho_L = 0._wp; rho_R = 0._wp
                                gamma_L = 0._wp; gamma_R = 0._wp
                                pi_inf_L = 0._wp; pi_inf_R = 0._wp
                                qv_L = 0._wp; qv_R = 0._wp

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    alpha_rho_L(i) = qL_prim_rsx_vf(${SF('')}$, i)
                                    alpha_rho_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, i)
                                    alpha_L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%E + i)
                                    alpha_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E + i)
                                end do

                                vel_L_rms = 0._wp; vel_R_rms = 0._wp

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel_L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%cont%end + i)
                                    vel_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%cont%end + i)
                                    vel_L_rms = vel_L_rms + vel_L(i)**2._wp
                                    vel_R_rms = vel_R_rms + vel_R(i)**2._wp
                                end do

                                ! Retain this in the refactor
                                if (mpp_lim .and. (num_fluids > 2)) then
                                    call s_accumulate_mixture_properties(num_fluids, alpha_rho_L, alpha_L, rho_L, gamma_L, &
                                                                         & pi_inf_L, qv_L)
                                    call s_accumulate_mixture_properties(num_fluids, alpha_rho_R, alpha_R, rho_R, gamma_R, &
                                                                         & pi_inf_R, qv_R)
                                else if (num_fluids > 2) then
                                    call s_accumulate_mixture_properties(num_fluids - 1, alpha_rho_L, alpha_L, rho_L, gamma_L, &
                                                                         & pi_inf_L, qv_L)
                                    call s_accumulate_mixture_properties(num_fluids - 1, alpha_rho_R, alpha_R, rho_R, gamma_R, &
                                                                         & pi_inf_R, qv_R)
                                else
                                    rho_L = qL_prim_rsx_vf(${SF('')}$, 1)
                                    gamma_L = gammas(1)
                                    pi_inf_L = pi_infs(1)
                                    qv_L = qvs(1)
                                    rho_R = qR_prim_rsx_vf(${SF(' + 1')}$, 1)
                                    gamma_R = gammas(1)
                                    pi_inf_R = pi_infs(1)
                                    qv_R = qvs(1)
                                end if

                                if (viscous) then
                                    if (num_fluids == 1) then  ! Need to consider case with num_fluids >= 2
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, 2
                                            Re_L(i) = dflt_real
                                            Re_R(i) = dflt_real

                                            if (Re_size(i) > 0) Re_L(i) = 0._wp
                                            if (Re_size(i) > 0) Re_R(i) = 0._wp

                                            $:GPU_LOOP(parallelism='[seq]')
                                            do q = 1, Re_size(i)
                                                Re_L(i) = (1._wp - qL_prim_rsx_vf(${SF('')}$, eqn_idx%E + Re_idx(i, &
                                                     & q)))/Res_gs(i, q) + Re_L(i)
                                                Re_R(i) = (1._wp - qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E + Re_idx(i, &
                                                     & q)))/Res_gs(i, q) + Re_R(i)
                                            end do

                                            Re_L(i) = 1._wp/max(Re_L(i), sgm_eps)
                                            Re_R(i) = 1._wp/max(Re_R(i), sgm_eps)
                                        end do
                                    end if
                                end if

                                pres_L = qL_prim_rsx_vf(${SF('')}$, eqn_idx%E)
                                pres_R = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E)

                                E_L = gamma_L*pres_L + pi_inf_L + 5.e-1_wp*rho_L*vel_L_rms
                                E_R = gamma_R*pres_R + pi_inf_R + 5.e-1_wp*rho_R*vel_R_rms

                                H_L = (E_L + pres_L)/rho_L
                                H_R = (E_R + pres_R)/rho_R

                                if (avg_state == avg_state_arithmetic) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, nb
                                        R0_L(i) = qL_prim_rsx_vf(${SF('')}$, rs(i))
                                        R0_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, rs(i))

                                        V0_L(i) = qL_prim_rsx_vf(${SF('')}$, vs(i))
                                        V0_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, vs(i))
                                        if (.not. polytropic .and. .not. qbmm) then
                                            P0_L(i) = qL_prim_rsx_vf(${SF('')}$, ps(i))
                                            P0_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, ps(i))
                                        end if
                                    end do

                                    if (.not. qbmm) then
                                        if (adv_n) then
                                            nbub_L = qL_prim_rsx_vf(${SF('')}$, eqn_idx%n)
                                            nbub_R = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%n)
                                        else
                                            nbub_L = 0._wp
                                            nbub_R = 0._wp
                                            $:GPU_LOOP(parallelism='[seq]')
                                            do i = 1, nb
                                                nbub_L = nbub_L + (R0_L(i)**3._wp)*weight(i)
                                                nbub_R = nbub_R + (R0_R(i)**3._wp)*weight(i)
                                            end do

                                            nbub_L = (3._wp/(4._wp*pi))*qL_prim_rsx_vf(${SF('')}$, eqn_idx%E + num_fluids)/nbub_L
                                            nbub_R = (3._wp/(4._wp*pi))*qR_prim_rsx_vf(${SF(' + 1')}$, &
                                                      & eqn_idx%E + num_fluids)/nbub_R
                                        end if
                                    else
                                        ! nb stored in 0th moment of first R0 bin in variable conversion module
                                        nbub_L = qL_prim_rsx_vf(${SF('')}$, eqn_idx%bub%beg)
                                        nbub_R = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%bub%beg)
                                    end if

                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, nb
                                        if (.not. qbmm) then
                                            pbw_L(i) = f_cpbw_KM(R0(i), R0_L(i), V0_L(i), P0_L(i))
                                            pbw_R(i) = f_cpbw_KM(R0(i), R0_R(i), V0_R(i), P0_R(i))
                                        end if
                                    end do

                                    if (qbmm) then
                                        PbwR3Lbar = mom_sp_rsx_vf(${SF('')}$, 4)
                                        PbwR3Rbar = mom_sp_rsx_vf(${SF(' + 1')}$, 4)

                                        R3Lbar = mom_sp_rsx_vf(${SF('')}$, 1)
                                        R3Rbar = mom_sp_rsx_vf(${SF(' + 1')}$, 1)

                                        R3V2Lbar = mom_sp_rsx_vf(${SF('')}$, 3)
                                        R3V2Rbar = mom_sp_rsx_vf(${SF(' + 1')}$, 3)
                                    else
                                        PbwR3Lbar = 0._wp
                                        PbwR3Rbar = 0._wp

                                        R3Lbar = 0._wp
                                        R3Rbar = 0._wp

                                        R3V2Lbar = 0._wp
                                        R3V2Rbar = 0._wp

                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, nb
                                            PbwR3Lbar = PbwR3Lbar + pbw_L(i)*(R0_L(i)**3._wp)*weight(i)
                                            PbwR3Rbar = PbwR3Rbar + pbw_R(i)*(R0_R(i)**3._wp)*weight(i)

                                            R3Lbar = R3Lbar + (R0_L(i)**3._wp)*weight(i)
                                            R3Rbar = R3Rbar + (R0_R(i)**3._wp)*weight(i)

                                            R3V2Lbar = R3V2Lbar + (R0_L(i)**3._wp)*(V0_L(i)**2._wp)*weight(i)
                                            R3V2Rbar = R3V2Rbar + (R0_R(i)**3._wp)*(V0_R(i)**2._wp)*weight(i)
                                        end do
                                    end if

                                    rho_avg = 5.e-1_wp*(rho_L + rho_R)
                                    H_avg = 5.e-1_wp*(H_L + H_R)
                                    gamma_avg = 5.e-1_wp*(gamma_L + gamma_R)
                                    qv_avg = 5.e-1_wp*(qv_L + qv_R)
                                    vel_avg_rms = 0._wp

                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        vel_avg_rms = vel_avg_rms + (5.e-1_wp*(vel_L(i) + vel_R(i)))**2._wp
                                    end do
                                end if

                                call s_compute_speed_of_sound(pres_L, rho_L, gamma_L, pi_inf_L, H_L, alpha_L, vel_L_rms, 0._wp, &
                                                              & c_L, qv_L)

                                call s_compute_speed_of_sound(pres_R, rho_R, gamma_R, pi_inf_R, H_R, alpha_R, vel_R_rms, 0._wp, &
                                                              & c_R, qv_R)

                                !> The computation of c_avg does not require all the variables, and therefore the non '_avg'
                                ! variables are placeholders to call the subroutine.
                                call s_compute_speed_of_sound(pres_R, rho_avg, gamma_avg, pi_inf_R, H_avg, alpha_R, vel_avg_rms, &
                                                              & 0._wp, c_avg, qv_avg)

                                if (viscous) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, 2
                                        Re_avg_rsx_vf(${SF('')}$, i) = 2._wp/(1._wp/Re_L(i) + 1._wp/Re_R(i))
                                    end do
                                end if

                                ! Low Mach correction
                                if (low_Mach == 2) then
                                    @:compute_low_Mach_correction()
                                end if

                                if (wave_speeds == wave_speeds_direct) then
                                    s_L = min(vel_L(dir_idx(1)) - c_L, vel_R(dir_idx(1)) - c_R)
                                    s_R = max(vel_R(dir_idx(1)) + c_R, vel_L(dir_idx(1)) + c_L)

                                    s_S = (pres_R - pres_L + rho_L*vel_L(dir_idx(1))*(s_L - vel_L(dir_idx(1))) &
                                           & - rho_R*vel_R(dir_idx(1))*(s_R - vel_R(dir_idx(1))))/(rho_L*(s_L - vel_L(dir_idx(1))) &
                                           & - rho_R*(s_R - vel_R(dir_idx(1))))
                                else if (wave_speeds == wave_speeds_pressure) then
                                    pres_SL = 5.e-1_wp*(pres_L + pres_R + rho_avg*c_avg*(vel_L(dir_idx(1)) - vel_R(dir_idx(1))))

                                    pres_SR = pres_SL

                                    ! Low Mach correction: Thornber et al. JCP (2008)
                                    Ms_L = max(1._wp, &
                                               & sqrt(1._wp + ((5.e-1_wp + gamma_L)/(1._wp + gamma_L))*(pres_SL/pres_L - 1._wp) &
                                               & *pres_L/((pres_L + pi_inf_L/(1._wp + gamma_L)))))
                                    Ms_R = max(1._wp, &
                                               & sqrt(1._wp + ((5.e-1_wp + gamma_R)/(1._wp + gamma_R))*(pres_SR/pres_R - 1._wp) &
                                               & *pres_R/((pres_R + pi_inf_R/(1._wp + gamma_R)))))

                                    s_L = vel_L(dir_idx(1)) - c_L*Ms_L
                                    s_R = vel_R(dir_idx(1)) + c_R*Ms_R

                                    s_S = 5.e-1_wp*((vel_L(dir_idx(1)) + vel_R(dir_idx(1))) + (pres_L - pres_R)/(rho_avg*c_avg))
                                end if

                                ! follows Einfeldt et al. s_M/P = min/max(0.,s_L/R)
                                s_M = min(0._wp, s_L); s_P = max(0._wp, s_R)

                                ! goes with q_star_L/R = xi_L/R * (variable) xi_L/R = ( ( s_L/R - u_L/R )/(s_L/R - s_star) )
                                xi_L = (s_L - vel_L(dir_idx(1)))/min(s_L - s_S, -sgm_eps)
                                xi_R = (s_R - vel_R(dir_idx(1)))/max(s_R - s_S, sgm_eps)
                                xi_L_m1 = (s_S - vel_L(dir_idx(1)))/min(s_L - s_S, -sgm_eps)
                                xi_R_m1 = (s_S - vel_R(dir_idx(1)))/max(s_R - s_S, sgm_eps)

                                ! goes with numerical velocity in x/y/z directions xi_P/M = 0.5 +/m sgn(0.5,s_star)
                                xi_M = (5.e-1_wp + sign(5.e-1_wp, s_S))
                                xi_P = (5.e-1_wp - sign(5.e-1_wp, s_S))

                                ! Low Mach correction
                                if (low_Mach == 1) then
                                    @:compute_low_Mach_correction()
                                else
                                    pcorr = 0._wp
                                end if

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, eqn_idx%cont%end
                                    flux_rsx_vf(${SF('')}$, i) = xi_M*qL_prim_rsx_vf(${SF('')}$, &
                                                & i)*(vel_L(dir_idx(1)) + s_M*xi_L_m1) + xi_P*qR_prim_rsx_vf(${SF(' + 1')}$, &
                                                & i)*(vel_R(dir_idx(1)) + s_P*xi_R_m1)
                                end do

                                if (bubbles_euler .and. (num_fluids > 1)) then
                                    ! Kill mass transport @ gas density
                                    flux_rsx_vf(${SF('')}$, eqn_idx%cont%end) = 0._wp
                                end if

                                ! Momentum flux. f = \rho u u + p I, q = \rho u, q_star = \xi * \rho*(s_star, v, w)

                                ! Include p_tilde

                                if (avg_state == avg_state_arithmetic) then
                                    if (alpha_L(num_fluids) < small_alf .or. R3Lbar < small_alf) then
                                        pres_L = pres_L - alpha_L(num_fluids)*pres_L
                                    else
                                        pres_L = pres_L - alpha_L(num_fluids)*(pres_L - PbwR3Lbar/R3Lbar - rho_L*R3V2Lbar/R3Lbar)
                                    end if

                                    if (alpha_R(num_fluids) < small_alf .or. R3Rbar < small_alf) then
                                        pres_R = pres_R - alpha_R(num_fluids)*pres_R
                                    else
                                        pres_R = pres_R - alpha_R(num_fluids)*(pres_R - PbwR3Rbar/R3Rbar - rho_R*R3V2Rbar/R3Rbar)
                                    end if
                                end if

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    flux_rsx_vf(${SF('')}$, &
                                                & eqn_idx%cont%end + dir_idx(i)) = xi_M*(rho_L*(vel_L(dir_idx(1))*vel_L(dir_idx(i) &
                                                & ) + s_M*(xi_L*(dir_flg(dir_idx(i))*s_S + (1._wp - dir_flg(dir_idx(i))) &
                                                & *vel_L(dir_idx(i))) - vel_L(dir_idx(i)))) + dir_flg(dir_idx(i))*(pres_L)) &
                                                & + xi_P*(rho_R*(vel_R(dir_idx(1))*vel_R(dir_idx(i)) &
                                                & + s_P*(xi_R*(dir_flg(dir_idx(i))*s_S + (1._wp - dir_flg(dir_idx(i))) &
                                                & *vel_R(dir_idx(i))) - vel_R(dir_idx(i)))) + dir_flg(dir_idx(i))*(pres_R)) &
                                                & + (s_M/s_L)*(s_P/s_R)*dir_flg(dir_idx(i))*pcorr
                                end do

                                ! Energy flux. f = u*(E+p), q = E, q_star = \xi*E+(s-u)(\rho s_star + p/(s-u))
                                flux_rsx_vf(${SF('')}$, &
                                            & eqn_idx%E) = xi_M*(vel_L(dir_idx(1))*(E_L + pres_L) + s_M*(xi_L*(E_L + (s_S &
                                            & - vel_L(dir_idx(1)))*(rho_L*s_S + (pres_L)/(s_L - vel_L(dir_idx(1))))) - E_L)) &
                                            & + xi_P*(vel_R(dir_idx(1))*(E_R + pres_R) + s_P*(xi_R*(E_R + (s_S - vel_R(dir_idx(1)) &
                                            & )*(rho_R*s_S + (pres_R)/(s_R - vel_R(dir_idx(1))))) - E_R)) + (s_M/s_L)*(s_P/s_R) &
                                            & *pcorr*s_S

                                ! Volume fraction flux
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = eqn_idx%adv%beg, eqn_idx%adv%end
                                    flux_rsx_vf(${SF('')}$, i) = xi_M*qL_prim_rsx_vf(${SF('')}$, &
                                                & i)*(vel_L(dir_idx(1)) + s_M*xi_L_m1) + xi_P*qR_prim_rsx_vf(${SF(' + 1')}$, &
                                                & i)*(vel_R(dir_idx(1)) + s_P*xi_R_m1)
                                end do

                                ! Advection velocity source: interface velocity for volume fraction transport
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel_src_rsx_vf(${SF('')}$, &
                                                   & dir_idx(i)) = xi_M*(vel_L(dir_idx(i)) + dir_flg(dir_idx(i))*s_M*xi_L_m1) &
                                                   & + xi_P*(vel_R(dir_idx(i)) + dir_flg(dir_idx(i))*s_P*xi_R_m1)

                                    ! IF ( (model_eqns == 4) .or. (num_fluids==1) ) vel_src_rs_vf(dir_idx(i))%sf(j,k,l) = 0._wp
                                end do

                                flux_src_rsx_vf(${SF('')}$, eqn_idx%adv%beg) = vel_src_rsx_vf(${SF('')}$, dir_idx(1))

                                ! Add advection flux for bubble variables
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = eqn_idx%bub%beg, eqn_idx%bub%end
                                    flux_rsx_vf(${SF('')}$, i) = xi_M*nbub_L*qL_prim_rsx_vf(${SF('')}$, &
                                                & i)*(vel_L(dir_idx(1)) + s_M*xi_L_m1) &
                                                & + xi_P*nbub_R*qR_prim_rsx_vf(${SF(' + 1')}$, i)*(vel_R(dir_idx(1)) + s_P*xi_R_m1)
                                end do

                                if (qbmm) then
                                    flux_rsx_vf(${SF('')}$, &
                                                & eqn_idx%bub%beg) = xi_M*nbub_L*(vel_L(dir_idx(1)) + s_M*xi_L_m1) &
                                                & + xi_P*nbub_R*(vel_R(dir_idx(1)) + s_P*xi_R_m1)
                                end if

                                if (adv_n) then
                                    flux_rsx_vf(${SF('')}$, &
                                                & eqn_idx%n) = xi_M*nbub_L*(vel_L(dir_idx(1)) + s_M*xi_L_m1) &
                                                & + xi_P*nbub_R*(vel_R(dir_idx(1)) + s_P*xi_R_m1)
                                end if

                                ! Geometrical source flux for cylindrical coordinates
                                #:if (NORM_DIR == 2)
                                    if (cyl_coord) then
                                        ! Substituting the advective flux into the inviscid geometrical source flux
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, eqn_idx%E
                                            flux_gsrc_rsx_vf(${SF('')}$, i) = flux_rsx_vf(${SF('')}$, i)
                                        end do
                                        ! Recalculating the radial momentum geometric source flux
                                        flux_gsrc_rsx_vf(${SF('')}$, &
                                                         & eqn_idx%cont%end + dir_idx(1)) &
                                                         & = f_compute_hllc_star_momentum_flux(rho_L, rho_R, vel_L(dir_idx(1)), &
                                                         & vel_R(dir_idx(1)), s_M, s_P, s_S, xi_L, xi_R, xi_M, xi_P, &
                                                         & dir_flg(dir_idx(1)))
                                        ! Geometrical source of the void fraction(s) is zero
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = eqn_idx%adv%beg, eqn_idx%adv%end
                                            flux_gsrc_rsx_vf(${SF('')}$, i) = 0._wp
                                        end do
                                    end if
                                #:endif
                                #:if (NORM_DIR == 3)
                                    if (grid_geometry == 3) then
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, sys_size
                                            flux_gsrc_rsx_vf(${SF('')}$, i) = 0._wp
                                        end do

                                        flux_gsrc_rsx_vf(${SF('')}$, &
                                                         & eqn_idx%mom%beg + 1) = -f_compute_hllc_star_momentum_flux(rho_L, &
                                                         & rho_R, vel_L(dir_idx(1)), vel_R(dir_idx(1)), s_M, s_P, s_S, xi_L, &
                                                         & xi_R, xi_M, xi_P, dir_flg(dir_idx(1)))
                                        flux_gsrc_rsx_vf(${SF('')}$, eqn_idx%mom%end) = flux_rsx_vf(${SF('')}$, eqn_idx%mom%beg + 1)
                                    end if
                                #:endif
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                else
                    ! 5-equation model (model_eqns=2): mixture total energy, volume fraction advection
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[i, T_L, T_R, vel_L_rms, vel_R_rms, pres_L, pres_R, rho_L, gamma_L, &
                                        & pi_inf_L, qv_L, rho_R, gamma_R, pi_inf_R, qv_R, alpha_L_sum, alpha_R_sum, E_L, E_R, &
                                        & MW_L, MW_R, R_gas_L, R_gas_R, Cp_L, Cp_R, Cv_L, Cv_R, Gamm_L, Gamm_R, Y_L, Y_R, H_L, &
                                        & H_R, qv_avg, rho_avg, gamma_avg, H_avg, c_L, c_R, c_avg, s_P, s_M, xi_P, xi_M, xi_L, &
                                        & xi_R, xi_L_m1, xi_R_m1, Ms_L, Ms_R, pres_SL, pres_SR, vel_L, vel_R, Re_L, Re_R, &
                                        & alpha_L, alpha_R, alpha_rho_L, alpha_rho_R, alpha_lim_L, alpha_lim_R, s_L, s_R, s_S, &
                                        & vel_avg_rms, pcorr, zcoef, vel_L_tmp, vel_R_tmp, Ys_L, Ys_R, Xs_L, Xs_R, Gamma_iL, &
                                        & Gamma_iR, Cp_iL, Cp_iR, tau_e_L, tau_e_R, xi_field_L, xi_field_R, Yi_avg, Phi_avg, &
                                        & h_iL, h_iR, h_avg_2, G_L, G_R]', copyin='[is1, is2, is3]')
                    do l = ${Z_BND}$%beg, ${Z_BND}$%end
                        do k = ${Y_BND}$%beg, ${Y_BND}$%end
                            do j = ${X_BND}$%beg, ${X_BND}$%end
                                vel_L_rms = 0._wp; vel_R_rms = 0._wp
                                rho_L = 0._wp; rho_R = 0._wp
                                gamma_L = 0._wp; gamma_R = 0._wp
                                pi_inf_L = 0._wp; pi_inf_R = 0._wp
                                qv_L = 0._wp; qv_R = 0._wp
                                alpha_L_sum = 0._wp; alpha_R_sum = 0._wp

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    alpha_L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%E + i)
                                    alpha_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E + i)
                                end do

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel_L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%cont%end + i)
                                    vel_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%cont%end + i)
                                    vel_L_rms = vel_L_rms + vel_L(i)**2._wp
                                    vel_R_rms = vel_R_rms + vel_R(i)**2._wp
                                end do

                                pres_L = qL_prim_rsx_vf(${SF('')}$, eqn_idx%E)
                                pres_R = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E)

                                ! Change this by splitting it into the cases present in the bubbles_euler
                                if (mpp_lim) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        qL_prim_rsx_vf(${SF('')}$, i) = max(0._wp, qL_prim_rsx_vf(${SF('')}$, i))
                                        qL_prim_rsx_vf(${SF('')}$, eqn_idx%E + i) = min(max(0._wp, qL_prim_rsx_vf(${SF('')}$, &
                                                       & eqn_idx%E + i)), 1._wp)
                                        qR_prim_rsx_vf(${SF(' + 1')}$, i) = max(0._wp, qR_prim_rsx_vf(${SF(' + 1')}$, i))
                                        qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E + i) = min(max(0._wp, &
                                                       & qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E + i)), 1._wp)
                                        alpha_L_sum = alpha_L_sum + qL_prim_rsx_vf(${SF('')}$, eqn_idx%E + i)
                                        alpha_R_sum = alpha_R_sum + qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E + i)
                                    end do

                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        qL_prim_rsx_vf(${SF('')}$, eqn_idx%E + i) = qL_prim_rsx_vf(${SF('')}$, &
                                                       & eqn_idx%E + i)/max(alpha_L_sum, sgm_eps)
                                        qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E + i) = qR_prim_rsx_vf(${SF(' + 1')}$, &
                                                       & eqn_idx%E + i)/max(alpha_R_sum, sgm_eps)
                                    end do
                                end if

                                ! Post-limiter loads for the mixture properties; alpha_L/R keep the pre-limiter loads used
                                ! downstream
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    alpha_rho_L(i) = qL_prim_rsx_vf(${SF('')}$, i)
                                    alpha_rho_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, i)
                                    alpha_lim_L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%E + i)
                                    alpha_lim_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E + i)
                                end do

                                call s_accumulate_mixture_properties(num_fluids, alpha_rho_L, alpha_lim_L, rho_L, gamma_L, &
                                                                     & pi_inf_L, qv_L)
                                call s_accumulate_mixture_properties(num_fluids, alpha_rho_R, alpha_lim_R, rho_R, gamma_R, &
                                                                     & pi_inf_R, qv_R)

                                if (viscous) then
                                    call s_compute_interface_reynolds(alpha_L, Re_L)
                                    call s_compute_interface_reynolds(alpha_R, Re_R)
                                end if

                                if (chemistry) then
                                    c_sum_Yi_Phi = 0.0_wp
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = eqn_idx%species%beg, eqn_idx%species%end
                                        Ys_L(i - eqn_idx%species%beg + 1) = qL_prim_rsx_vf(${SF('')}$, i)
                                        Ys_R(i - eqn_idx%species%beg + 1) = qR_prim_rsx_vf(${SF(' + 1')}$, i)
                                    end do

                                    call get_mixture_molecular_weight(Ys_L, MW_L)
                                    call get_mixture_molecular_weight(Ys_R, MW_R)

                                    Xs_L(:) = Ys_L(:)*MW_L/molecular_weights(:)
                                    Xs_R(:) = Ys_R(:)*MW_R/molecular_weights(:)

                                    R_gas_L = gas_constant/MW_L
                                    R_gas_R = gas_constant/MW_R

                                    T_L = pres_L/rho_L/R_gas_L
                                    T_R = pres_R/rho_R/R_gas_R

                                    call get_species_specific_heats_r(T_L, Cp_iL)
                                    call get_species_specific_heats_r(T_R, Cp_iR)

                                    if (chem_params%gamma_method == 1) then
                                        !> gamma_method = 1: Ref. Section 2.3.1 Formulation of doi:10.7907/ZKW8-ES97.
                                        Gamma_iL = Cp_iL/(Cp_iL - 1.0_wp)
                                        Gamma_iR = Cp_iR/(Cp_iR - 1.0_wp)

                                        gamma_L = sum(Xs_L(:)/(Gamma_iL(:) - 1.0_wp))
                                        gamma_R = sum(Xs_R(:)/(Gamma_iR(:) - 1.0_wp))
                                    else if (chem_params%gamma_method == 2) then
                                        !> gamma_method = 2: c_p / c_v where c_p, c_v are specific heats.
                                        call get_mixture_specific_heat_cp_mass(T_L, Ys_L, Cp_L)
                                        call get_mixture_specific_heat_cp_mass(T_R, Ys_R, Cp_R)
                                        call get_mixture_specific_heat_cv_mass(T_L, Ys_L, Cv_L)
                                        call get_mixture_specific_heat_cv_mass(T_R, Ys_R, Cv_R)

                                        Gamm_L = Cp_L/Cv_L; Gamm_R = Cp_R/Cv_R
                                        gamma_L = 1.0_wp/(Gamm_L - 1.0_wp); gamma_R = 1.0_wp/(Gamm_R - 1.0_wp)
                                    end if

                                    call get_mixture_energy_mass(T_L, Ys_L, E_L)
                                    call get_mixture_energy_mass(T_R, Ys_R, E_R)

                                    E_L = rho_L*E_L + 5.e-1*rho_L*vel_L_rms
                                    E_R = rho_R*E_R + 5.e-1*rho_R*vel_R_rms
                                    H_L = (E_L + pres_L)/rho_L
                                    H_R = (E_R + pres_R)/rho_R
                                else
                                    E_L = gamma_L*pres_L + pi_inf_L + 5.e-1*rho_L*vel_L_rms + qv_L
                                    E_R = gamma_R*pres_R + pi_inf_R + 5.e-1*rho_R*vel_R_rms + qv_R

                                    H_L = (E_L + pres_L)/rho_L
                                    H_R = (E_R + pres_R)/rho_R
                                end if

                                ! ENERGY ADJUSTMENTS FOR HYPOELASTIC ENERGY
                                if (hypoelasticity) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, eqn_idx%stress%end - eqn_idx%stress%beg + 1
                                        tau_e_L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%stress%beg - 1 + i)
                                        tau_e_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%stress%beg - 1 + i)
                                    end do
                                    G_L = 0._wp
                                    G_R = 0._wp
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        G_L = G_L + alpha_L(i)*Gs_rs(i)
                                        G_R = G_R + alpha_R(i)*Gs_rs(i)
                                    end do
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, eqn_idx%stress%end - eqn_idx%stress%beg + 1
                                        ! Elastic contribution to energy if G large enough
                                        if ((G_L > verysmall) .and. (G_R > verysmall)) then
                                            E_L = E_L + (tau_e_L(i)*tau_e_L(i))/(4._wp*G_L)
                                            E_R = E_R + (tau_e_R(i)*tau_e_R(i))/(4._wp*G_R)
                                            ! Additional terms in 2D and 3D
                                            if ((i == 2) .or. (i == 4) .or. (i == 5)) then
                                                E_L = E_L + (tau_e_L(i)*tau_e_L(i))/(4._wp*G_L)
                                                E_R = E_R + (tau_e_R(i)*tau_e_R(i))/(4._wp*G_R)
                                            end if
                                        end if
                                    end do
                                end if

                                ! Hyperelastic stress contribution: strain energy added to total energy
                                if (hyperelasticity) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        xi_field_L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%xi%beg - 1 + i)
                                        xi_field_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%xi%beg - 1 + i)
                                    end do
                                    G_L = 0._wp
                                    G_R = 0._wp
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        ! Mixture left and right shear modulus
                                        G_L = G_L + alpha_L(i)*Gs_rs(i)
                                        G_R = G_R + alpha_R(i)*Gs_rs(i)
                                    end do
                                    ! Elastic contribution to energy if G large enough
                                    if (G_L > verysmall .and. G_R > verysmall) then
                                        E_L = E_L + G_L*qL_prim_rsx_vf(${SF('')}$, eqn_idx%xi%end + 1)
                                        E_R = E_R + G_R*qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%xi%end + 1)
                                    end if
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, b_size - 1
                                        tau_e_L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%stress%beg - 1 + i)
                                        tau_e_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%stress%beg - 1 + i)
                                    end do
                                end if

                                H_L = (E_L + pres_L)/rho_L
                                H_R = (E_R + pres_R)/rho_R

                                @:compute_average_state()

                                call s_compute_speed_of_sound(pres_L, rho_L, gamma_L, pi_inf_L, H_L, alpha_L, vel_L_rms, 0._wp, &
                                                              & c_L, qv_L)

                                call s_compute_speed_of_sound(pres_R, rho_R, gamma_R, pi_inf_R, H_R, alpha_R, vel_R_rms, 0._wp, &
                                                              & c_R, qv_R)

                                !> The computation of c_avg does not require all the variables, and therefore the non '_avg'
                                !  variables are placeholders to call the subroutine.
                                call s_compute_speed_of_sound(pres_R, rho_avg, gamma_avg, pi_inf_R, H_avg, alpha_R, vel_avg_rms, &
                                                              & c_sum_Yi_Phi, c_avg, qv_avg)

                                if (viscous) then
                                    if (chemistry) then
                                        call compute_viscosity_and_inversion(T_L, Ys_L, T_R, Ys_R, Re_L(1), Re_R(1))
                                    end if
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, 2
                                        Re_avg_rsx_vf(${SF('')}$, i) = 2._wp/(1._wp/Re_L(i) + 1._wp/Re_R(i))
                                    end do
                                end if

                                ! Low Mach correction
                                if (low_Mach == 2) then
                                    @:compute_low_Mach_correction()
                                end if

                                if (wave_speeds == wave_speeds_direct) then
                                    if (elasticity) then
                                        ! Elastic wave speed, Rodriguez et al. JCP (2019)
                                        s_L = min(vel_L(dir_idx(1)) - sqrt(c_L*c_L + (((4._wp*G_L)/3._wp) + tau_e_L(dir_idx_tau(1) &
                                                  & ))/rho_L), &
                                                  & vel_R(dir_idx(1)) - sqrt(c_R*c_R + (((4._wp*G_R)/3._wp) &
                                                  & + tau_e_R(dir_idx_tau(1)))/rho_R))
                                        s_R = max(vel_R(dir_idx(1)) + sqrt(c_R*c_R + (((4._wp*G_R)/3._wp) + tau_e_R(dir_idx_tau(1) &
                                                  & ))/rho_R), &
                                                  & vel_L(dir_idx(1)) + sqrt(c_L*c_L + (((4._wp*G_L)/3._wp) &
                                                  & + tau_e_L(dir_idx_tau(1)))/rho_L))
                                        s_S = (pres_R - tau_e_R(dir_idx_tau(1)) - pres_L + tau_e_L(dir_idx_tau(1)) &
                                               & + rho_L*vel_L(dir_idx(1))*(s_L - vel_L(dir_idx(1))) - rho_R*vel_R(dir_idx(1)) &
                                               & *(s_R - vel_R(dir_idx(1))))/(rho_L*(s_L - vel_L(dir_idx(1))) - rho_R*(s_R &
                                               & - vel_R(dir_idx(1))))
                                    else
                                        s_L = min(vel_L(dir_idx(1)) - c_L, vel_R(dir_idx(1)) - c_R)
                                        s_R = max(vel_R(dir_idx(1)) + c_R, vel_L(dir_idx(1)) + c_L)
                                        s_S = (pres_R - pres_L + rho_L*vel_L(dir_idx(1))*(s_L - vel_L(dir_idx(1))) &
                                               & - rho_R*vel_R(dir_idx(1))*(s_R - vel_R(dir_idx(1))))/(rho_L*(s_L &
                                               & - vel_L(dir_idx(1))) - rho_R*(s_R - vel_R(dir_idx(1))))
                                    end if
                                else if (wave_speeds == wave_speeds_pressure) then
                                    pres_SL = 5.e-1_wp*(pres_L + pres_R + rho_avg*c_avg*(vel_L(dir_idx(1)) - vel_R(dir_idx(1))))

                                    pres_SR = pres_SL

                                    ! Low Mach correction: Thornber et al. JCP (2008)
                                    Ms_L = max(1._wp, &
                                               & sqrt(1._wp + ((5.e-1_wp + gamma_L)/(1._wp + gamma_L))*(pres_SL/pres_L - 1._wp) &
                                               & *pres_L/((pres_L + pi_inf_L/(1._wp + gamma_L)))))
                                    Ms_R = max(1._wp, &
                                               & sqrt(1._wp + ((5.e-1_wp + gamma_R)/(1._wp + gamma_R))*(pres_SR/pres_R - 1._wp) &
                                               & *pres_R/((pres_R + pi_inf_R/(1._wp + gamma_R)))))

                                    s_L = vel_L(dir_idx(1)) - c_L*Ms_L
                                    s_R = vel_R(dir_idx(1)) + c_R*Ms_R

                                    s_S = 5.e-1_wp*((vel_L(dir_idx(1)) + vel_R(dir_idx(1))) + (pres_L - pres_R)/(rho_avg*c_avg))
                                end if

                                ! follows Einfeldt et al. s_M/P = min/max(0.,s_L/R)
                                s_M = min(0._wp, s_L); s_P = max(0._wp, s_R)

                                ! goes with q_star_L/R = xi_L/R * (variable) xi_L/R = ( ( s_L/R - u_L/R )/(s_L/R - s_star) )
                                xi_L = (s_L - vel_L(dir_idx(1)))/min(s_L - s_S, -sgm_eps)
                                xi_R = (s_R - vel_R(dir_idx(1)))/max(s_R - s_S, sgm_eps)
                                ! xi_L/R - 1 = (s_S - u_L/R)/(s_L/R - s_star): avoids cancellation when xi \approx 1
                                xi_L_m1 = (s_S - vel_L(dir_idx(1)))/min(s_L - s_S, -sgm_eps)
                                xi_R_m1 = (s_S - vel_R(dir_idx(1)))/max(s_R - s_S, sgm_eps)

                                ! goes with numerical velocity in x/y/z directions xi_P/M = 0.5 +/m sgn(0.5,s_star)
                                xi_M = (5.e-1_wp + sign(5.e-1_wp, s_S))
                                xi_P = (5.e-1_wp - sign(5.e-1_wp, s_S))

                                ! Low Mach correction
                                if (low_Mach == 1) then
                                    @:compute_low_Mach_correction()
                                else
                                    pcorr = 0._wp
                                end if

                                ! COMPUTING THE HLLC FLUXES MASS FLUX.
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, eqn_idx%cont%end
                                    flux_rsx_vf(${SF('')}$, i) = xi_M*qL_prim_rsx_vf(${SF('')}$, &
                                                & i)*(vel_L(dir_idx(1)) + s_M*xi_L_m1) + xi_P*qR_prim_rsx_vf(${SF(' + 1')}$, &
                                                & i)*(vel_R(dir_idx(1)) + s_P*xi_R_m1)
                                end do

                                ! MOMENTUM FLUX. f = \rho u u - \sigma, q = \rho u, q_star = \xi * \rho*(s_star, v, w) identity:
                                ! xi*(dir_flg*s_S+(1-dir_flg)*u_i)-u_i = (dir_flg*s_L/R+(1-dir_flg)*u_i)*xi_m1
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    flux_rsx_vf(${SF('')}$, &
                                                & eqn_idx%cont%end + dir_idx(i)) = xi_M*(rho_L*(vel_L(dir_idx(1))*vel_L(dir_idx(i) &
                                                & ) + s_M*(dir_flg(dir_idx(i))*s_L + (1._wp - dir_flg(dir_idx(i))) &
                                                & *vel_L(dir_idx(i)))*xi_L_m1) + dir_flg(dir_idx(i))*(pres_L)) &
                                                & + xi_P*(rho_R*(vel_R(dir_idx(1))*vel_R(dir_idx(i)) + s_P*(dir_flg(dir_idx(i)) &
                                                & *s_R + (1._wp - dir_flg(dir_idx(i)))*vel_R(dir_idx(i)))*xi_R_m1) &
                                                & + dir_flg(dir_idx(i))*(pres_R)) + (s_M/s_L)*(s_P/s_R)*dir_flg(dir_idx(i))*pcorr
                                end do

                                ! ENERGY FLUX. f = u*(E-\sigma), q = E, q_star = \xi*E+(s-u)(\rho s_star - \sigma/(s-u))
                                ! xi*(E+expr)-E = E*xi_m1 + xi*expr avoids E*(xi-1) cancellation
                                flux_rsx_vf(${SF('')}$, &
                                            & eqn_idx%E) = xi_M*(vel_L(dir_idx(1))*(E_L + pres_L) + s_M*(E_L*xi_L_m1 + xi_L*(s_S &
                                            & - vel_L(dir_idx(1)))*(rho_L*s_S + pres_L/(s_L - vel_L(dir_idx(1)))))) &
                                            & + xi_P*(vel_R(dir_idx(1))*(E_R + pres_R) + s_P*(E_R*xi_R_m1 + xi_R*(s_S &
                                            & - vel_R(dir_idx(1)))*(rho_R*s_S + pres_R/(s_R - vel_R(dir_idx(1)))))) + (s_M/s_L) &
                                            & *(s_P/s_R)*pcorr*s_S

                                ! ELASTICITY. Elastic shear stress additions for the momentum and energy flux
                                if (elasticity) then
                                    flux_ene_e = 0._wp
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        ! MOMENTUM ELASTIC FLUX.
                                        flux_rsx_vf(${SF('')}$, eqn_idx%cont%end + dir_idx(i)) = flux_rsx_vf(${SF('')}$, &
                                                    & eqn_idx%cont%end + dir_idx(i)) - xi_M*tau_e_L(dir_idx_tau(i)) &
                                                    & - xi_P*tau_e_R(dir_idx_tau(i))
                                        ! ENERGY ELASTIC FLUX.
                                        flux_ene_e = flux_ene_e - xi_M*(vel_L(dir_idx(i))*tau_e_L(dir_idx_tau(i)) &
                                                                        & + s_M*(xi_L*((s_S - vel_L(i))*(tau_e_L(dir_idx_tau(i)) &
                                                                        & /(s_L - vel_L(i)))))) - xi_P*(vel_R(dir_idx(i)) &
                                                                        & *tau_e_R(dir_idx_tau(i)) + s_P*(xi_R*((s_S - vel_R(i)) &
                                                                        & *(tau_e_R(dir_idx_tau(i))/(s_R - vel_R(i))))))
                                    end do
                                    flux_rsx_vf(${SF('')}$, eqn_idx%E) = flux_rsx_vf(${SF('')}$, eqn_idx%E) + flux_ene_e
                                end if

                                ! HYPOELASTIC STRESS EVOLUTION FLUX.
                                if (hypoelasticity) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, eqn_idx%stress%end - eqn_idx%stress%beg + 1
                                        flux_rsx_vf(${SF('')}$, &
                                                    & eqn_idx%stress%beg - 1 + i) = xi_M*(s_S/(s_L - s_S))*(s_L*rho_L*tau_e_L(i) &
                                                    & - rho_L*vel_L(dir_idx(1))*tau_e_L(i)) + xi_P*(s_S/(s_R - s_S)) &
                                                    & *(s_R*rho_R*tau_e_R(i) - rho_R*vel_R(dir_idx(1))*tau_e_R(i))
                                    end do
                                end if

                                ! VOLUME FRACTION FLUX.
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = eqn_idx%adv%beg, eqn_idx%adv%end
                                    flux_rsx_vf(${SF('')}$, i) = xi_M*qL_prim_rsx_vf(${SF('')}$, &
                                                & i)*(vel_L(dir_idx(1)) + s_M*xi_L_m1) + xi_P*qR_prim_rsx_vf(${SF(' + 1')}$, &
                                                & i)*(vel_R(dir_idx(1)) + s_P*xi_R_m1)
                                end do

                                ! VOLUME FRACTION SOURCE FLUX.
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel_src_rsx_vf(${SF('')}$, &
                                                   & dir_idx(i)) = xi_M*(vel_L(dir_idx(i)) + dir_flg(dir_idx(i))*s_M*xi_L_m1) &
                                                   & + xi_P*(vel_R(dir_idx(i)) + dir_flg(dir_idx(i))*s_P*xi_R_m1)
                                end do

                                ! COLOR FUNCTION FLUX
                                if (surface_tension) then
                                    flux_rsx_vf(${SF('')}$, eqn_idx%c) = xi_M*qL_prim_rsx_vf(${SF('')}$, &
                                                & eqn_idx%c)*(vel_L(dir_idx(1)) + s_M*xi_L_m1) &
                                                & + xi_P*qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%c)*(vel_R(dir_idx(1)) + s_P*xi_R_m1)
                                end if

                                ! Hyperelastic reference map flux for material deformation tracking
                                if (hyperelasticity) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        flux_rsx_vf(${SF('')}$, &
                                                    & eqn_idx%xi%beg - 1 + i) = xi_M*(s_S/(s_L - s_S))*(s_L*rho_L*xi_field_L(i) &
                                                    & - rho_L*vel_L(dir_idx(1))*xi_field_L(i)) + xi_P*(s_S/(s_R - s_S)) &
                                                    & *(s_R*rho_R*xi_field_R(i) - rho_R*vel_R(dir_idx(1))*xi_field_R(i))
                                    end do
                                end if

                                flux_src_rsx_vf(${SF('')}$, eqn_idx%adv%beg) = vel_src_rsx_vf(${SF('')}$, dir_idx(1))

                                if (chemistry) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = eqn_idx%species%beg, eqn_idx%species%end
                                        Y_L = qL_prim_rsx_vf(${SF('')}$, i)
                                        Y_R = qR_prim_rsx_vf(${SF(' + 1')}$, i)

                                        flux_rsx_vf(${SF('')}$, &
                                                    & i) = xi_M*rho_L*Y_L*(vel_L(dir_idx(1)) + s_M*xi_L_m1) &
                                                    & + xi_P*rho_R*Y_R*(vel_R(dir_idx(1)) + s_P*xi_R_m1)
                                        flux_src_rsx_vf(${SF('')}$, i) = 0.0_wp
                                    end do
                                end if

                                ! Geometrical source flux for cylindrical coordinates
                                #:if (NORM_DIR == 2)
                                    if (cyl_coord) then
                                        ! Substituting the advective flux into the inviscid geometrical source flux
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, eqn_idx%E
                                            flux_gsrc_rsx_vf(${SF('')}$, i) = flux_rsx_vf(${SF('')}$, i)
                                        end do
                                        ! Recalculating the radial momentum geometric source flux
                                        flux_gsrc_rsx_vf(${SF('')}$, &
                                                         & eqn_idx%cont%end + dir_idx(1)) &
                                                         & = f_compute_hllc_star_momentum_flux(rho_L, rho_R, vel_L(dir_idx(1)), &
                                                         & vel_R(dir_idx(1)), s_M, s_P, s_S, xi_L, xi_R, xi_M, xi_P, &
                                                         & dir_flg(dir_idx(1)))
                                        ! Geometrical source of the void fraction(s) is zero
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = eqn_idx%adv%beg, eqn_idx%adv%end
                                            flux_gsrc_rsx_vf(${SF('')}$, i) = 0._wp
                                        end do
                                    end if
                                #:endif
                                #:if (NORM_DIR == 3)
                                    if (grid_geometry == 3) then
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, sys_size
                                            flux_gsrc_rsx_vf(${SF('')}$, i) = 0._wp
                                        end do

                                        flux_gsrc_rsx_vf(${SF('')}$, &
                                                         & eqn_idx%mom%beg + 1) = -f_compute_hllc_star_momentum_flux(rho_L, &
                                                         & rho_R, vel_L(dir_idx(1)), vel_R(dir_idx(1)), s_M, s_P, s_S, xi_L, &
                                                         & xi_R, xi_M, xi_P, dir_flg(dir_idx(1)))
                                        flux_gsrc_rsx_vf(${SF('')}$, eqn_idx%mom%end) = flux_rsx_vf(${SF('')}$, eqn_idx%mom%beg + 1)
                                    end if
                                #:endif
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
            end if
        #:endfor
        ! Computing HLLC flux and source flux for Euler system of equations

        if (viscous) then
            if (weno_Re_flux) then
                call s_compute_viscous_source_flux(qL_prim_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                                                   & dqL_prim_dx_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                                                   & dqL_prim_dy_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                                                   & dqL_prim_dz_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                                                   & qR_prim_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                                                   & dqR_prim_dx_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                                                   & dqR_prim_dy_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                                                   & dqR_prim_dz_vf(eqn_idx%mom%beg:eqn_idx%mom%end), flux_src_vf, q_prim_vf, &
                                                   & norm_dir, ix, iy, iz)
            else
                call s_compute_viscous_source_flux(q_prim_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                                                   & dqL_prim_dx_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                                                   & dqL_prim_dy_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                                                   & dqL_prim_dz_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                                                   & q_prim_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                                                   & dqR_prim_dx_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                                                   & dqR_prim_dy_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                                                   & dqR_prim_dz_vf(eqn_idx%mom%beg:eqn_idx%mom%end), flux_src_vf, q_prim_vf, &
                                                   & norm_dir, ix, iy, iz)
            end if
        end if

        if (surface_tension) then
            call s_compute_capillary_source_flux(vel_src_rsx_vf, flux_src_vf, norm_dir, isx, isy, isz)
        end if

        call s_finalize_riemann_solver(flux_vf, flux_src_vf, flux_gsrc_vf, norm_dir)

    end subroutine s_hllc_riemann_solver

end module m_riemann_solver_hllc
