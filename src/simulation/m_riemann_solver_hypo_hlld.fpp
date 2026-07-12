!>
!! @file
!! @brief Contains module m_riemann_solver_hypo_hlld

!> @brief Dual-pass HLLD approximate Riemann solver for hypoelastic flows, with non-conservative interface-velocity coupling
#:include 'case.fpp'
#:include 'macros.fpp'
#:include 'inline_riemann.fpp'

! Single source of truth for the per-component HLL flux on the (U_L, U_R, F_L, F_R) compact
! basis: used by the degenerate-fan fallback and the ADC blend, which must stay consistent.
! Textual inlining; codegen is identical to the materialized expression.
#:def hll_flux_component(LHS, I)
    ${LHS}$ = (S_R*F_L(${I}$) - S_L*F_R(${I}$) + S_L*S_R*(U_R(${I}$) - U_L(${I}$)))/(S_R - S_L + verysmall)
#:enddef

module m_riemann_solver_hypo_hlld

    use m_derived_types
    use m_global_parameters
    use m_variables_conversion
    use m_riemann_state

    implicit none

contains

    !> Classify the position of the interface (xi = 0) in the five-wave HLLD fan: 0 left of S_L, 1..4 for the four inner wedges, 5
    !! right of S_R. Single definition of the wave-fan geometry shared by the flux fold, the NC face-velocity export, and the
    !! axisymmetric face-state pick, which must stay consistent. In the degenerate-fan branch it is called on the fallback wave
    !! speeds and reproduces their one-sided selection.
    integer function f_hlld_wave_zone(S_L, S_Lstar, s_M, S_Rstar, S_R) result(zone)

        $:GPU_ROUTINE(function_name='f_hlld_wave_zone', parallelism='[seq]', cray_inline=True)
        real(wp), intent(in) :: S_L, S_Lstar, s_M, S_Rstar, S_R

        if (0._wp <= S_L) then
            zone = 0
        else if (0._wp <= S_Lstar) then
            zone = 1
        else if (0._wp <= s_M) then
            zone = 2
        else if (0._wp <= S_Rstar) then
            zone = 3
        else if (0._wp <= S_R) then
            zone = 4
        else
            zone = 5
        end if

    end function f_hlld_wave_zone

    !> HLLD Riemann solver resolves all 5 waves for the hypoelastic equations: 1 entropy wave, 2 shear stress waves, 2 fast waves.
    subroutine s_hypo_hlld_riemann_solver(qL_prim_rsx_vf, dqL_prim_dx_vf, dqL_prim_dy_vf, dqL_prim_dz_vf, qL_prim_vf, &
                                          & qR_prim_rsx_vf, dqR_prim_dx_vf, dqR_prim_dy_vf, dqR_prim_dz_vf, qR_prim_vf, &
                                          & q_prim_vf, flux_vf, flux_src_vf, flux_gsrc_vf, norm_dir, ix, iy, iz)

        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout) :: qL_prim_rsx_vf, qR_prim_rsx_vf
        type(scalar_field), allocatable, dimension(:), intent(inout) :: dqL_prim_dx_vf, dqR_prim_dx_vf, dqL_prim_dy_vf, &
             & dqR_prim_dy_vf, dqL_prim_dz_vf, dqR_prim_dz_vf

        type(scalar_field), allocatable, dimension(:), intent(inout) :: qL_prim_vf, qR_prim_vf
        type(scalar_field), dimension(sys_size), intent(in)          :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout)       :: flux_vf, flux_src_vf, flux_gsrc_vf
        integer, intent(in)                                          :: norm_dir
        type(int_bounds_info), intent(in)                            :: ix, iy, iz

        ! Local variables:

        ! Hypo HLLD has num_fluids==2 guard -> reads alpha_* at a hard-coded index 2.
        ! With case-opt, `num_fluids=1` makes dimension(num_fluids) length 1,
        ! so amdflang gives a compile-time out-of-bounds error. max(num_fluids, 2) fixes it

        #:if MFC_CASE_OPTIMIZATION
            real(wp), dimension(${max(num_fluids, 2)}$) :: alpha_L, alpha_R, alpha_rho_L, alpha_rho_R
        #:elif USING_AMD
            real(wp), dimension(3) :: alpha_L, alpha_R, alpha_rho_L, alpha_rho_R
        #:else
            real(wp), dimension(num_fluids) :: alpha_L, alpha_R, alpha_rho_L, alpha_rho_R
        #:endif
        type(riemann_states_vec3) :: vel
        type(riemann_states)      :: rho, pres, E, H
        type(riemann_states)      :: gamma, pi_inf, qv
        type(riemann_states)      :: vel_rms
        type(riemann_states)      :: c

        ! HLLD speeds and intermediate state variables:
        real(wp)                :: S_L, S_R, s_M, S_Lstar, S_Rstar
        real(wp)                :: pTot_L, pTot_R, rhoL_star, rhoR_star
        real(wp), dimension(14) :: U_L, U_R
        real(wp), dimension(14) :: F_L, F_R, F_hlld
        real(wp)                :: us_c, uss_c  ! selected-side U_star / U_starstar, one component at a time (register diet)
        real(wp)                :: F_HLL_c  ! per-component HLL flux for ADC blending (diet: replaces the F_HLL array)
        real(wp)                :: U_HLL_c  ! per-component HLL state for the ADC axisym trace (diet: replaces the U_HLL array)
        real(wp)                :: rho_HLL, u_n_HLL_cons, tau_nn_HLL
        real(wp)                :: u_n_HLL_trace, u_t_HLL_trace
        real(wp)                :: p_face_HLL, tau_qq_face_HLL
        integer                 :: ncomp  ! 11 for 2D/axisym, 14 for 3D Cartesian

        ! HLLD Hypo variables

        real(wp) :: C_NC, sqrtC_NC
        real(wp) :: A_L, A_R, denomA, fac_L, fac_R
        real(wp) :: u_n_L, u_t_L, u_n_R, u_t_R
        real(wp) :: u_t2_L, u_t2_R
        real(wp) :: tau_nn_L, tau_nt_L, tau_tt_L, tau_nn_R, tau_nt_R, tau_tt_R
        real(wp) :: tau_nt2_L, tau_nt2_R, tau_t2t2_L, tau_t2t2_R, tau_t1t2_L, tau_t1t2_R
        real(wp) :: tau_qq_L, tau_qq_R
        real(wp) :: G_L, G_R
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(6) :: tau_e_L, tau_e_R
        #:else
            real(wp), dimension(eqn_idx%stress%end - eqn_idx%stress%beg + 1) :: tau_e_L, tau_e_R
        #:endif

        real(wp) :: alpha1_L_star, alpha1_R_star, alpha2_L_star, alpha2_R_star
        real(wp) :: u_t_star, tau_nt_star
        real(wp) :: u_t2_star, tau_nt2_star
        real(wp) :: tau_nn_L_star, tau_nn_R_star, tau_tt_L_star, tau_tt_R_star
        real(wp) :: tau_tt_L_starstar, tau_tt_R_starstar
        real(wp) :: tau_t2t2_L_star, tau_t2t2_R_star
        real(wp) :: tau_t2t2_L_starstar, tau_t2t2_R_starstar
        real(wp) :: tau_t1t2_L_star, tau_t1t2_R_star
        real(wp) :: tau_t1t2_L_starstar, tau_t1t2_R_starstar
        real(wp) :: tau_qq_L_star, tau_qq_R_star
        real(wp) :: pTot_star
        real(wp) :: E_L_star, E_R_star
        real(wp) :: E_L_starstar, E_R_starstar
        real(wp) :: p_face, tau_qq_face
        real(wp) :: u_n_face, u_t_face
        real(wp) :: G_hat
        real(wp) :: rho_hat
        real(wp) :: tau_nn_hat, tau_nt_hat, tau_tt_hat, tau_qq_hat
        real(wp) :: tau_nt2_hat, tau_t2t2_hat, tau_t1t2_hat
        ! alpha_hat/alpha_rho_hat: same max(num_fluids, 2) reason as alpha_* above
        #:if MFC_CASE_OPTIMIZATION
            real(wp), dimension(${max(num_fluids, 2)}$)                      :: alpha_hat, alpha_rho_hat
            real(wp), dimension(eqn_idx%stress%end - eqn_idx%stress%beg + 1) :: tau_e_hat
        #:elif USING_AMD
            real(wp), dimension(3) :: alpha_hat, alpha_rho_hat
            real(wp), dimension(6) :: tau_e_hat
        #:else
            real(wp), dimension(num_fluids)                                  :: alpha_hat, alpha_rho_hat
            real(wp), dimension(eqn_idx%stress%end - eqn_idx%stress%beg + 1) :: tau_e_hat
        #:endif

        real(wp)            :: pres_hat, blkmod1_hat, blkmod2_hat, K_hat
        real(wp)            :: C_hat_1, C_hat_2
        real(wp)            :: Sigma_L, Sigma_R, dSigma, Sigma_ref
        real(wp)            :: a_L_ref, a_R_ref, a_ref
        real(wp)            :: du_t, dtau_nt, du_t2, dtau_nt2
        real(wp)            :: sensor_ptot, sensor_vt, sensor_tnt, sensor_combined
        real(wp)            :: phi
        real(wp), parameter :: ADC_power = 1.0_wp
        real(wp)            :: alpha_L_sum, alpha_R_sum
        logical             :: degenerate
        integer             :: i, j, k, l, ipass, zone

        call s_populate_riemann_states_variables_buffers(qL_prim_rsx_vf, dqL_prim_dx_vf, dqL_prim_dy_vf, dqL_prim_dz_vf, &
            & qR_prim_rsx_vf, dqR_prim_dx_vf, dqR_prim_dy_vf, dqR_prim_dz_vf, norm_dir, ix, iy, iz)

        call s_initialize_riemann_solver(flux_src_vf, norm_dir)

        #:for NORM_DIR, XYZ, STENCIL_VAR, COORDS, X_BND, Y_BND, Z_BND in &
                    [(1, 'x', 'j', '{STENCIL_IDX}, k, l', 'is1', 'is2', 'is3'), &
                     (2, 'y', 'k', 'j, {STENCIL_IDX}, l', 'is2', 'is1', 'is3'), &
                     (3, 'z', 'l', 'j, k, {STENCIL_IDX}', 'is3', 'is2', 'is1')]
            #:set SV = STENCIL_VAR
            #:set SF = lambda offs: COORDS.format(STENCIL_IDX = SV + offs)
            if (norm_dir == ${NORM_DIR}$) then
                ! Anchor-cell index pattern: the fused kernel reads both anchors (hat_L: face cell, hat_R: face cell + 1)
                ! directly from q_prim_vf; loop indices are physical, so the offset rides the stencil index.
                #:set HATIDX = SF(' + ipass - 1')
                ! The fused HLLD kernel privatizes ~150 scalars -- a private() list large enough that
                ! the single clause must be wrapped across continuation lines. The four fragments
                ! (_hlld_p1..p4) are ONLY for source readability; fypp concatenates them into one
                ! clause below. That wrapping is FOLD_DIRECTIVE's job -- its within-clause comma split
                ! exists for exactly this case -- not the fragments'.
                #:set _hlld_p1 = '[i,j,k,l,ipass,degenerate,alpha_rho_L,alpha_rho_R,vel,alpha_L,alpha_R,rho,pres,E,H,gamma,pi_inf,qv,vel_rms,c,S_L,S_R,s_M,S_Lstar,S_Rstar,pTot_L,pTot_R,rhoL_star,rhoR_star,U_L,U_R,F_L,F_R,F_hlld,us_c,uss_c,zone,F_HLL_c,U_HLL_c,rho_HLL,u_n_HLL_cons,tau_nn_HLL,u_n_HLL_trace,u_t_HLL_trace,p_face_HLL,tau_qq_face_HLL,ncomp,C_NC,sqrtC_NC,A_L,A_R,denomA,fac_L,fac_R,'
                #:set _hlld_p2 = 'u_n_L,u_t_L,u_n_R,u_t_R,u_t2_L,u_t2_R,tau_nn_L,tau_nt_L,tau_tt_L,tau_nn_R,tau_nt_R,tau_tt_R,tau_nt2_L,tau_nt2_R,tau_t2t2_L,tau_t2t2_R,tau_t1t2_L,tau_t1t2_R,tau_qq_L,tau_qq_R,G_L,G_R,tau_e_L,tau_e_R,alpha1_L_star,alpha1_R_star,alpha2_L_star,alpha2_R_star,u_t_star,tau_nt_star,u_t2_star,tau_nt2_star,tau_nn_L_star,tau_nn_R_star,tau_tt_L_star,tau_tt_R_star,tau_tt_L_starstar,tau_tt_R_starstar,'
                #:set _hlld_p3 = 'tau_t2t2_L_star,tau_t2t2_R_star,tau_t2t2_L_starstar,tau_t2t2_R_starstar,tau_t1t2_L_star,tau_t1t2_R_star,tau_t1t2_L_starstar,tau_t1t2_R_starstar,tau_qq_L_star,tau_qq_R_star,pTot_star,E_L_star,E_R_star,E_L_starstar,E_R_starstar,p_face,tau_qq_face,u_n_face,u_t_face,G_hat,rho_hat,tau_nn_hat,tau_nt_hat,tau_tt_hat,tau_qq_hat,tau_nt2_hat,tau_t2t2_hat,tau_t1t2_hat,'
                #:set _hlld_p4 = 'alpha_hat,alpha_rho_hat,tau_e_hat,pres_hat,blkmod1_hat,blkmod2_hat,K_hat,C_hat_1,C_hat_2,Sigma_L,Sigma_R,dSigma,Sigma_ref,a_L_ref,a_R_ref,a_ref,du_t,dtau_nt,du_t2,dtau_nt2,sensor_ptot,sensor_vt,sensor_tnt,sensor_combined,phi,alpha_L_sum,alpha_R_sum]'
                ! Wave-fan side table for the per-component F_hlld fold below: side name, the side's two zones,
                ! its starstar zone, and the outer/inner wave speeds. The L and R sides are mirror images.
                #:set HLLD_FAN_SIDES = [('L', 1, 2, 2, 'S_L', 'S_Lstar'), ('R', 3, 4, 3, 'S_R', 'S_Rstar')]
                ! Wave-fan component tables for the per-component fold below: (us, uss) expression pairs in
                ! F_hlld component order; uss = None means the inner (starstar) correction leaves the
                ! component unchanged (uss_c = us_c). '{S}' expands to the fan side (L/R).
                #:set HLLD_FAN_COMPS_3D = [('U_{S}(1)*fac_{S}', None), &
                    & ('U_{S}(2)*fac_{S}', None), &
                    & ('rho{S}_star*S_M', None), &
                    & ('rho{S}_star*u_t_{S}', 'rho{S}_star*u_t_star'), &
                    & ('rho{S}_star*u_t2_{S}', 'rho{S}_star*u_t2_star'), &
                    & ('E_{S}_star', 'E_{S}_starstar'), &
                    & ('alpha1_{S}_star', None), &
                    & ('alpha2_{S}_star', None), &
                    & ('rho{S}_star*tau_nn_{S}_star', None), &
                    & ('rho{S}_star*tau_nt_{S}', 'rho{S}_star*tau_nt_star'), &
                    & ('rho{S}_star*tau_nt2_{S}', 'rho{S}_star*tau_nt2_star'), &
                    & ('rho{S}_star*tau_tt_{S}_star', 'rho{S}_star*tau_tt_{S}_starstar'), &
                    & ('rho{S}_star*tau_t2t2_{S}_star', 'rho{S}_star*tau_t2t2_{S}_starstar'), &
                    & ('rho{S}_star*tau_t1t2_{S}_star', 'rho{S}_star*tau_t1t2_{S}_starstar')]
                #:set HLLD_FAN_COMPS_2D = [('U_{S}(1)*fac_{S}', None), &
                    & ('U_{S}(2)*fac_{S}', None), &
                    & ('rho{S}_star*S_M', None), &
                    & ('rho{S}_star*u_t_{S}', 'rho{S}_star*u_t_star'), &
                    & ('E_{S}_star', 'E_{S}_starstar'), &
                    & ('alpha1_{S}_star', None), &
                    & ('alpha2_{S}_star', None), &
                    & ('rho{S}_star*tau_nn_{S}_star', None), &
                    & ('rho{S}_star*tau_nt_{S}', 'rho{S}_star*tau_nt_star'), &
                    & ('rho{S}_star*tau_tt_{S}_star', 'rho{S}_star*tau_tt_{S}_starstar'), &
                    & ('rho{S}_star*tau_qq_{S}_star', None)]
                $:GPU_PARALLEL_LOOP(collapse=3, private=_hlld_p1 + _hlld_p2 + _hlld_p3 + _hlld_p4)
                do l = ${Z_BND}$%beg, ${Z_BND}$%end
                    do k = ${Y_BND}$%beg, ${Y_BND}$%end
                        do j = ${X_BND}$%beg, ${X_BND}$%end
                            ! Extract left/right primitive states

                            do i = 1, eqn_idx%cont%end
                                alpha_rho_L(i) = qL_prim_rsx_vf(${SF('')}$, i)
                                alpha_rho_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, i)
                            end do

                            ! Initialize the unused third velocity component in 2D
                            vel%L = 0._wp
                            vel%R = 0._wp

                            ! Load velocities in physical-coordinate order; map directional aliases below via dir_idx
                            do i = 1, num_vels
                                vel%L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%cont%end + i)
                                vel%R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%cont%end + i)
                            end do

                            vel_rms%L = vel%L(1)**2 + vel%L(2)**2 + vel%L(3)**2
                            vel_rms%R = vel%R(1)**2 + vel%R(2)**2 + vel%R(3)**2

                            do i = 1, num_fluids
                                alpha_L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%E + i)
                                alpha_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E + i)
                            end do

                            ! Clamp and renormalize volume fractions when mpp_lim is on (match HLL/HLLC)
                            alpha_L_sum = 0._wp
                            alpha_R_sum = 0._wp
                            if (mpp_lim) then
                                do i = 1, num_fluids
                                    alpha_rho_L(i) = max(0._wp, alpha_rho_L(i))
                                    alpha_L(i) = min(max(0._wp, alpha_L(i)), 1._wp)
                                    alpha_L_sum = alpha_L_sum + alpha_L(i)
                                end do
                                alpha_L = alpha_L/max(alpha_L_sum, sgm_eps)

                                do i = 1, num_fluids
                                    alpha_rho_R(i) = max(0._wp, alpha_rho_R(i))
                                    alpha_R(i) = min(max(0._wp, alpha_R(i)), 1._wp)
                                    alpha_R_sum = alpha_R_sum + alpha_R(i)
                                end do
                                alpha_R = alpha_R/max(alpha_R_sum, sgm_eps)
                            end if

                            pres%L = qL_prim_rsx_vf(${SF('')}$, eqn_idx%E)
                            pres%R = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E)

                            ! Hypoelasticity
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, eqn_idx%stress%end - eqn_idx%stress%beg + 1
                                tau_e_L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%stress%beg - 1 + i)
                                tau_e_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%stress%beg - 1 + i)
                            end do

                            u_t2_L = 0._wp; u_t2_R = 0._wp
                            tau_nt2_L = 0._wp; tau_nt2_R = 0._wp
                            tau_t2t2_L = 0._wp; tau_t2t2_R = 0._wp
                            tau_t1t2_L = 0._wp; tau_t1t2_R = 0._wp

                            ! Map physical-basis arrays to directional aliases via stress_perm/dir_idx
                            u_n_L = vel%L(dir_idx(1)); u_n_R = vel%R(dir_idx(1))
                            tau_nn_L = tau_e_L(stress_perm(1)); tau_nn_R = tau_e_R(stress_perm(1))
                            if (n == 0) then
                                ncomp = 11
                            else if (p == 0) then
                                ncomp = 11
                                u_t_L = vel%L(dir_idx(2)); u_t_R = vel%R(dir_idx(2))
                                tau_nt_L = tau_e_L(stress_perm(2)); tau_nt_R = tau_e_R(stress_perm(2))
                                tau_tt_L = tau_e_L(stress_perm(3)); tau_tt_R = tau_e_R(stress_perm(3))
                            else
                                ncomp = 14
                                u_t_L = vel%L(dir_idx(2)); u_t_R = vel%R(dir_idx(2))
                                u_t2_L = vel%L(dir_idx(3)); u_t2_R = vel%R(dir_idx(3))
                                tau_nt_L = tau_e_L(stress_perm(2)); tau_nt_R = tau_e_R(stress_perm(2))
                                tau_tt_L = tau_e_L(stress_perm(3)); tau_tt_R = tau_e_R(stress_perm(3))
                                tau_nt2_L = tau_e_L(stress_perm(4)); tau_nt2_R = tau_e_R(stress_perm(4))
                                tau_t1t2_L = tau_e_L(stress_perm(5)); tau_t1t2_R = tau_e_R(stress_perm(5))
                                tau_t2t2_L = tau_e_L(stress_perm(6)); tau_t2t2_R = tau_e_R(stress_perm(6))
                            end if
                            if (cyl_coord) then
                                tau_qq_L = tau_e_L(eqn_idx%stress%end - eqn_idx%stress%beg + 1)
                                tau_qq_R = tau_e_R(eqn_idx%stress%end - eqn_idx%stress%beg + 1)
                            else
                                tau_qq_L = 0._wp; tau_qq_R = 0._wp
                            end if
                            ! Total pressure (replace the usual pressure to define SM)
                            pTot_L = pres%L - tau_nn_L
                            pTot_R = pres%R - tau_nn_R

                            ! Symmetrize total pressure when it differs only by floating-point roundoff. WENO reconstruction of a
                            ! uniform field can produce slightly different L/R values at material interfaces due to different
                            ! smoothness indicators. With stiff materials (G~1e9), even 1e-12 relative pTot asymmetry creates O(1)
                            ! spurious stress through the HLLD star-state.
                            if (abs(pTot_R - pTot_L) < 1e-12_wp*max(abs(pTot_L), abs(pTot_R), 1._wp)) then
                                pTot_L = 5e-1_wp*(pTot_L + pTot_R)
                                pTot_R = pTot_L
                            end if

                            ! Sum properties of all fluid components
                            rho%L = 0._wp; gamma%L = 0._wp; pi_inf%L = 0._wp; qv%L = 0._wp
                            rho%R = 0._wp; gamma%R = 0._wp; pi_inf%R = 0._wp; qv%R = 0._wp
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                rho%L = rho%L + alpha_rho_L(i)
                                gamma%L = gamma%L + alpha_L(i)*gammas(i)
                                pi_inf%L = pi_inf%L + alpha_L(i)*pi_infs(i)
                                qv%L = qv%L + alpha_rho_L(i)*qvs(i)

                                rho%R = rho%R + alpha_rho_R(i)
                                gamma%R = gamma%R + alpha_R(i)*gammas(i)
                                pi_inf%R = pi_inf%R + alpha_R(i)*pi_infs(i)
                                qv%R = qv%R + alpha_rho_R(i)*qvs(i)
                            end do

                            G_L = 0._wp; G_R = 0._wp
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                G_L = G_L + alpha_L(i)*Gs_rs(i)
                                G_R = G_R + alpha_R(i)*Gs_rs(i)
                            end do

                            E%L = gamma%L*pres%L + pi_inf%L + 5e-1_wp*rho%L*vel_rms%L + qv%L
                            E%R = gamma%R*pres%R + pi_inf%R + 5e-1_wp*rho%R*vel_rms%R + qv%R

                            ! Freeze the thermal/kinetic enthalpy used by the EOS sound-speed call before
                            ! adding hypoelastic strain energy to the conservative total energy.
                            H%L = (E%L + pres%L)/rho%L
                            H%R = (E%R + pres%R)/rho%R

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, eqn_idx%stress%end - eqn_idx%stress%beg + 1
                                @:compute_hypo_elastic_energy(E%L, E%R, &
                                                              & (n > 0 .and. p == 0 .and. i == 2) .or. (p > 0 .and. (i == 2 &
                                                              & .or. i == 4 .or. i == 5)))
                            end do

                            ! Compute Riemann states

                            call s_compute_speed_of_sound(pres%L, rho%L, gamma%L, pi_inf%L, H%L, alpha_L, vel_rms%L, 0._wp, c%L, &
                                                          & qv%L)
                            call s_compute_speed_of_sound(pres%R, rho%R, gamma%R, pi_inf%R, H%R, alpha_R, vel_rms%R, 0._wp, c%R, &
                                                          & qv%R)

                            S_L = min(u_n_L - sqrt(max(verysmall, c%L*c%L + ((4._wp/3._wp)*G_L + tau_nn_L)/rho%L)), &
                                      & u_n_R - sqrt(max(verysmall, c%R*c%R + ((4._wp/3._wp)*G_R + tau_nn_R)/rho%R)))
                            S_R = max(u_n_R + sqrt(max(verysmall, c%R*c%R + ((4._wp/3._wp)*G_R + tau_nn_R)/rho%R)), &
                                      & u_n_L + sqrt(max(verysmall, c%L*c%L + ((4._wp/3._wp)*G_L + tau_nn_L)/rho%L)))

                            if (p > 0 .and. .not. cyl_coord) then
                                ! 3D Cartesian: 14-state compact basis
                                U_L(1) = alpha_rho_L(1); U_R(1) = alpha_rho_R(1)
                                U_L(2) = alpha_rho_L(2); U_R(2) = alpha_rho_R(2)
                                U_L(3) = rho%L*u_n_L; U_R(3) = rho%R*u_n_R
                                U_L(4) = rho%L*u_t_L; U_R(4) = rho%R*u_t_R
                                U_L(5) = rho%L*u_t2_L; U_R(5) = rho%R*u_t2_R
                                U_L(6) = E%L; U_R(6) = E%R
                                U_L(7) = alpha_L(1); U_R(7) = alpha_R(1)
                                U_L(8) = alpha_L(2); U_R(8) = alpha_R(2)
                                U_L(9) = rho%L*tau_nn_L; U_R(9) = rho%R*tau_nn_R
                                U_L(10) = rho%L*tau_nt_L; U_R(10) = rho%R*tau_nt_R
                                U_L(11) = rho%L*tau_nt2_L; U_R(11) = rho%R*tau_nt2_R
                                U_L(12) = rho%L*tau_tt_L; U_R(12) = rho%R*tau_tt_R
                                U_L(13) = rho%L*tau_t2t2_L; U_R(13) = rho%R*tau_t2t2_R
                                U_L(14) = rho%L*tau_t1t2_L; U_R(14) = rho%R*tau_t1t2_R

                                F_L(1) = U_L(1)*u_n_L; F_R(1) = U_R(1)*u_n_R
                                F_L(2) = U_L(2)*u_n_L; F_R(2) = U_R(2)*u_n_R
                                F_L(3) = rho%L*u_n_L*u_n_L + pTot_L
                                F_R(3) = rho%R*u_n_R*u_n_R + pTot_R
                                F_L(4) = rho%L*u_n_L*u_t_L - tau_nt_L
                                F_R(4) = rho%R*u_n_R*u_t_R - tau_nt_R
                                F_L(5) = rho%L*u_n_L*u_t2_L - tau_nt2_L
                                F_R(5) = rho%R*u_n_R*u_t2_R - tau_nt2_R
                                F_L(6) = (E%L + pTot_L)*u_n_L - u_t_L*tau_nt_L - u_t2_L*tau_nt2_L
                                F_R(6) = (E%R + pTot_R)*u_n_R - u_t_R*tau_nt_R - u_t2_R*tau_nt2_R
                                ! Rows 7-14 (volume fractions + stresses) carry the anchor-dependent hat coefficients and are filled
                                ! inside the per-anchor pass loop below.
                            else
                                ! 2D/axisym: 11-state compact basis (unchanged)
                                U_L(1) = alpha_rho_L(1)
                                U_L(2) = alpha_rho_L(2)
                                U_L(3) = rho%L*u_n_L
                                U_L(4) = rho%L*u_t_L
                                U_L(5) = E%L
                                U_L(6) = alpha_L(1)
                                U_L(7) = alpha_L(2)
                                U_L(8) = rho%L*tau_nn_L
                                U_L(9) = rho%L*tau_nt_L
                                U_L(10) = rho%L*tau_tt_L
                                U_L(11) = rho%L*tau_qq_L

                                U_R(1) = alpha_rho_R(1)
                                U_R(2) = alpha_rho_R(2)
                                U_R(3) = rho%R*u_n_R
                                U_R(4) = rho%R*u_t_R
                                U_R(5) = E%R
                                U_R(6) = alpha_R(1)
                                U_R(7) = alpha_R(2)
                                U_R(8) = rho%R*tau_nn_R
                                U_R(9) = rho%R*tau_nt_R
                                U_R(10) = rho%R*tau_tt_R
                                U_R(11) = rho%R*tau_qq_R

                                F_L(1) = U_L(1)*u_n_L
                                F_L(2) = U_L(2)*u_n_L
                                F_L(3) = rho%L*u_n_L*u_n_L + pTot_L
                                F_L(4) = rho%L*u_n_L*u_t_L - tau_nt_L
                                F_L(5) = (E%L + pTot_L)*u_n_L - u_t_L*tau_nt_L

                                F_R(1) = U_R(1)*u_n_R
                                F_R(2) = U_R(2)*u_n_R
                                F_R(3) = rho%R*u_n_R*u_n_R + pTot_R
                                F_R(4) = rho%R*u_n_R*u_t_R - tau_nt_R
                                F_R(5) = (E%R + pTot_R)*u_n_R - u_t_R*tau_nt_R
                                ! Rows 6-11 (volume fractions + stresses) carry the anchor-dependent hat coefficients and are filled
                                ! inside the per-anchor pass loop below.
                            end if

                            A_L = rho%L*(S_L - u_n_L)
                            A_R = rho%R*(S_R - u_n_R)
                            denomA = (A_R - A_L)

                            S_M = ((pTot_R - pTot_L) + A_L*u_n_L - A_R*u_n_R)/(A_L - A_R + verysmall)

                            ! Degenerate wave structure: denom ~ 0 or S_M not in [S_L,S_R]. The test is anchor-independent, so both
                            ! anchored solves below take the same branch.
                            degenerate = (abs(denomA) < verysmall .or. .not. (S_L - verysmall <= S_M .and. S_M <= S_R + verysmall))

                            if (.not. degenerate) then
                                ! Anchor-independent pieces of the star region (the anchor state enters only the volume-fraction and
                                ! stress star states)
                                pTot_star = pTot_L + A_L*(S_M - u_n_L)

                                ! Signed floor (HLLC-style, see m_riemann_solver_hllc.fpp) on the wave-fan denominators: guards the
                                ! degenerate S_M -> S_L / S_M -> S_R limit, sign-preserving, and a no-op away from it (|denom| >>
                                ! verysmall)
                                rhoL_star = rho%L*(S_L - u_n_L)/min(S_L - S_M, -verysmall)
                                rhoR_star = rho%R*(S_R - u_n_R)/max(S_R - S_M, verysmall)
                                fac_L = (S_L - u_n_L)/min(S_L - S_M, -verysmall)
                                fac_R = (S_R - u_n_R)/max(S_R - S_M, verysmall)

                                E_L_star = (E%L*(u_n_L - S_L) + u_n_L*pTot_L - S_M*pTot_star)/max(S_M - S_L, verysmall)
                                E_R_star = (E%R*(u_n_R - S_R) + u_n_R*pTot_R - S_M*pTot_star)/min(S_M - S_R, -verysmall)

                                if (riemann_hypo_ADC) then
                                    ! ADC sensors depend only on the L/R face states: computed once, shared by both anchored solves
                                    Sigma_L = pTot_L
                                    Sigma_R = pTot_R
                                    dSigma = Sigma_R - Sigma_L
                                    Sigma_ref = max(max(abs(Sigma_L), abs(Sigma_R)), verysmall)

                                    a_L_ref = sqrt(max(verysmall, c%L*c%L + ((4._wp/3._wp)*G_L + tau_nn_L)/rho%L))
                                    a_R_ref = sqrt(max(verysmall, c%R*c%R + ((4._wp/3._wp)*G_R + tau_nn_R)/rho%R))
                                    a_ref = max(max(a_L_ref, a_R_ref), verysmall)

                                    du_t = u_t_R - u_t_L
                                    dtau_nt = tau_nt_R - tau_nt_L
                                    du_t2 = u_t2_R - u_t2_L
                                    dtau_nt2 = tau_nt2_R - tau_nt2_L

                                    sensor_ptot = (dSigma*dSigma)/((ADC_kappa*Sigma_ref)**2 + verysmall)
                                    sensor_vt = (du_t*du_t + du_t2*du_t2)/((ADC_kappa*a_ref)**2 + verysmall)
                                    sensor_tnt = (dtau_nt*dtau_nt + dtau_nt2*dtau_nt2)/((ADC_kappa*Sigma_ref)**2 + verysmall)

                                    sensor_combined = sensor_ptot + sensor_tnt + sensor_vt

                                    phi = exp(-(sensor_combined**ADC_power))
                                end if
                            end if

                            ! Fused dual-pass: both anchored solves share everything above. The anchor (hat) state and everything it
                            ! touches is evaluated per pass: pass 1 anchors on cell j (hat_L) and writes the flux_rs* set; pass 2
                            ! anchors on cell j+1 (hat_R) and writes the flux_hatR_rs* set.
                            $:GPU_LOOP(parallelism='[seq]')
                            do ipass = 1, 2
                                do i = 1, eqn_idx%cont%end
                                    alpha_rho_hat(i) = q_prim_vf(i)%sf(${HATIDX}$)
                                end do
                                do i = 1, num_fluids
                                    alpha_hat(i) = q_prim_vf(eqn_idx%E + i)%sf(${HATIDX}$)
                                end do
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, eqn_idx%stress%end - eqn_idx%stress%beg + 1
                                    tau_e_hat(i) = q_prim_vf(eqn_idx%stress%beg - 1 + i)%sf(${HATIDX}$)
                                end do

                                ! Anchor-state directional aliases (mirrors the L/R alias block above)
                                tau_nt2_hat = 0._wp; tau_t2t2_hat = 0._wp; tau_t1t2_hat = 0._wp
                                tau_nn_hat = tau_e_hat(stress_perm(1))
                                if (n > 0) then
                                    tau_nt_hat = tau_e_hat(stress_perm(2))
                                    tau_tt_hat = tau_e_hat(stress_perm(3))
                                    if (p > 0) then
                                        tau_nt2_hat = tau_e_hat(stress_perm(4))
                                        tau_t1t2_hat = tau_e_hat(stress_perm(5))
                                        tau_t2t2_hat = tau_e_hat(stress_perm(6))
                                    end if
                                end if
                                if (cyl_coord) then
                                    tau_qq_hat = tau_e_hat(eqn_idx%stress%end - eqn_idx%stress%beg + 1)
                                else
                                    tau_qq_hat = 0._wp
                                end if

                                rho_hat = 0._wp
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    rho_hat = rho_hat + alpha_rho_hat(i)
                                end do

                                G_hat = 0._wp
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    G_hat = G_hat + alpha_hat(i)*Gs_rs(i)
                                end do

                                ! Two-component 2D only (enforced by checker restrictions)
                                K_hat = 0._wp
                                if (alt_soundspeed) then
                                    pres_hat = q_prim_vf(eqn_idx%E)%sf(${HATIDX}$)
                                    blkmod1_hat = ((gammas(1) + 1._wp)*pres_hat + pi_infs(1))/gammas(1) + (4._wp/3._wp)*Gs_rs(1)
                                    blkmod2_hat = ((gammas(2) + 1._wp)*pres_hat + pi_infs(2))/gammas(2) + (4._wp/3._wp)*Gs_rs(2)
                                    K_hat = alpha_hat(1)*alpha_hat(2)*(blkmod2_hat - blkmod1_hat)/(alpha_hat(1)*blkmod2_hat &
                                                      & + alpha_hat(2)*blkmod1_hat + verysmall)
                                end if
                                C_hat_1 = alpha_hat(1) + K_hat
                                C_hat_2 = alpha_hat(2) - K_hat

                                if (p > 0 .and. .not. cyl_coord) then
                                    ! 3D Cartesian: anchor-dependent rows (7-14) of the 14-state flux
                                    F_L(7) = U_L(7)*u_n_L - C_hat_1*u_n_L
                                    F_R(7) = U_R(7)*u_n_R - C_hat_1*u_n_R
                                    F_L(8) = U_L(8)*u_n_L - C_hat_2*u_n_L
                                    F_R(8) = U_R(8)*u_n_R - C_hat_2*u_n_R
                                    F_L(9) = U_L(9)*u_n_L - rho_hat*(4._wp/3._wp*G_hat + tau_nn_hat)*u_n_L
                                    F_R(9) = U_R(9)*u_n_R - rho_hat*(4._wp/3._wp*G_hat + tau_nn_hat)*u_n_R
                                    F_L(10) = U_L(10)*u_n_L - rho_hat*(G_hat + tau_nn_hat)*u_t_L
                                    F_R(10) = U_R(10)*u_n_R - rho_hat*(G_hat + tau_nn_hat)*u_t_R
                                    F_L(11) = U_L(11)*u_n_L - rho_hat*(G_hat + tau_nn_hat)*u_t2_L
                                    F_R(11) = U_R(11)*u_n_R - rho_hat*(G_hat + tau_nn_hat)*u_t2_R
                                    F_L(12) = U_L(12)*u_n_L + rho_hat*(2._wp/3._wp*G_hat + tau_tt_hat)*u_n_L &
                                        & - 2._wp*rho_hat*tau_nt_hat*u_t_L
                                    F_R(12) = U_R(12)*u_n_R + rho_hat*(2._wp/3._wp*G_hat + tau_tt_hat)*u_n_R &
                                        & - 2._wp*rho_hat*tau_nt_hat*u_t_R
                                    F_L(13) = U_L(13)*u_n_L + rho_hat*(2._wp/3._wp*G_hat + tau_t2t2_hat)*u_n_L &
                                        & - 2._wp*rho_hat*tau_nt2_hat*u_t2_L
                                    F_R(13) = U_R(13)*u_n_R + rho_hat*(2._wp/3._wp*G_hat + tau_t2t2_hat)*u_n_R &
                                        & - 2._wp*rho_hat*tau_nt2_hat*u_t2_R
                                    F_L(14) = U_L(14)*u_n_L + rho_hat*tau_t1t2_hat*u_n_L - rho_hat*tau_nt2_hat*u_t_L &
                                        & - rho_hat*tau_nt_hat*u_t2_L
                                    F_R(14) = U_R(14)*u_n_R + rho_hat*tau_t1t2_hat*u_n_R - rho_hat*tau_nt2_hat*u_t_R &
                                        & - rho_hat*tau_nt_hat*u_t2_R
                                else
                                    ! 2D/axisym: anchor-dependent rows (6-11) of the 11-state flux
                                    F_L(6) = U_L(6)*u_n_L - C_hat_1*u_n_L
                                    F_R(6) = U_R(6)*u_n_R - C_hat_1*u_n_R
                                    F_L(7) = U_L(7)*u_n_L - C_hat_2*u_n_L
                                    F_R(7) = U_R(7)*u_n_R - C_hat_2*u_n_R
                                    F_L(8) = U_L(8)*u_n_L - rho_hat*(4._wp/3._wp*G_hat + tau_nn_hat)*u_n_L
                                    F_R(8) = U_R(8)*u_n_R - rho_hat*(4._wp/3._wp*G_hat + tau_nn_hat)*u_n_R
                                    F_L(9) = U_L(9)*u_n_L - rho_hat*(G_hat + tau_nn_hat)*u_t_L
                                    F_R(9) = U_R(9)*u_n_R - rho_hat*(G_hat + tau_nn_hat)*u_t_R
                                    F_L(10) = U_L(10)*u_n_L + rho_hat*(2._wp/3._wp*G_hat + tau_tt_hat)*u_n_L &
                                        & - 2._wp*rho_hat*tau_nt_hat*u_t_L
                                    F_R(10) = U_R(10)*u_n_R + rho_hat*(2._wp/3._wp*G_hat + tau_tt_hat)*u_n_R &
                                        & - 2._wp*rho_hat*tau_nt_hat*u_t_R
                                    F_L(11) = U_L(11)*u_n_L + rho_hat*(2._wp/3._wp*G_hat + tau_qq_hat)*u_n_L
                                    F_R(11) = U_R(11)*u_n_R + rho_hat*(2._wp/3._wp*G_hat + tau_qq_hat)*u_n_R
                                end if

                                if (degenerate) then
                                    ! HLL (or one-sided) fallback for degenerate wave structure
                                    if (S_L < 0._wp .and. S_R > 0._wp) then
                                        do i = 1, ncomp
                                            @:hll_flux_component(F_hlld(i), i)
                                        end do
                                    else if (S_L >= 0._wp) then
                                        F_hlld(1:ncomp) = F_L(1:ncomp)
                                    else
                                        F_hlld(1:ncomp) = F_R(1:ncomp)
                                    end if
                                    ! Initialize star-state variables to safe fallback values so subsequent face-state and
                                    ! axisymmetric source logic never reads uninitialized memory regardless of wave configuration.
                                    pTot_star = 5e-1_wp*(pTot_L + pTot_R)
                                    S_Lstar = S_L
                                    S_Rstar = S_R
                                    u_t_star = 5e-1_wp*(u_t_L + u_t_R)
                                    tau_nn_L_star = tau_nn_L
                                    tau_nn_R_star = tau_nn_R
                                    tau_qq_L_star = tau_qq_L
                                    tau_qq_R_star = tau_qq_R
                                    zone = f_hlld_wave_zone(S_L, S_Lstar, s_M, S_Rstar, S_R)
                                else
                                    C_NC = rho_hat*(G_hat + tau_nn_hat)

                                    if (C_NC < verysmall) then
                                        ! Degenerate shear impedance: collapse inner waves to HLLC
                                        S_Lstar = S_M
                                        S_Rstar = S_M
                                        u_t_star = 5e-1_wp*(u_t_L + u_t_R)
                                        tau_nt_star = 5e-1_wp*(tau_nt_L + tau_nt_R)
                                        u_t2_star = 5e-1_wp*(u_t2_L + u_t2_R)
                                        tau_nt2_star = 5e-1_wp*(tau_nt2_L + tau_nt2_R)
                                    else
                                        sqrtC_NC = sqrt(C_NC)
                                        S_Lstar = S_M - sqrtC_NC/rhoL_star
                                        S_Rstar = S_M + sqrtC_NC/rhoR_star
                                        u_t_star = 5e-1_wp*((tau_nt_R - tau_nt_L)/sqrtC_NC + (u_t_R + u_t_L))
                                        tau_nt_star = 5e-1_wp*((u_t_R - u_t_L)*sqrtC_NC + (tau_nt_R + tau_nt_L))
                                        u_t2_star = 5e-1_wp*((tau_nt2_R - tau_nt2_L)/sqrtC_NC + (u_t2_R + u_t2_L))
                                        tau_nt2_star = 5e-1_wp*((u_t2_R - u_t2_L)*sqrtC_NC + (tau_nt2_R + tau_nt2_L))
                                    end if

                                    ! Outer-wave stress star states: one formula family over the components. ELC is
                                    ! the anchor-side elastic coefficient (deviatoric 4/3 vs 2/3 structure; t1t2
                                    ! carries no G term) and SGN = -1/+1 the jump sign of the component.
                                    #:for TAU, SGN, ELC in [('tau_nn', -1, '(G_hat*4._wp/3._wp + tau_nn_hat)'), &
                                        & ('tau_tt', 1, '(G_hat*2._wp/3._wp + tau_tt_hat)'), &
                                        & ('tau_t2t2', 1, '(G_hat*2._wp/3._wp + tau_t2t2_hat)'), &
                                        & ('tau_t1t2', 1, 'tau_t1t2_hat'), &
                                        & ('tau_qq', 1, '(G_hat*2._wp/3._wp + tau_qq_hat)')]
                                        #:for S, SOUT in [('L', 'S_L'), ('R', 'S_R')]
                                            ${TAU}$_${S}$_star = ${TAU}$_${S}$ ${'-' if SGN < 0 else '+'}$ (rho_hat*${ELC}$*(u_n_${S}$ - S_M))/(rho%${S}$*(u_n_${S}$ - ${SOUT}$))
                                        #:endfor
                                    #:endfor

                                    if (C_NC < verysmall) then
                                        ! Degenerate: no inner wave correction
                                        tau_tt_L_starstar = tau_tt_L_star
                                        tau_tt_R_starstar = tau_tt_R_star
                                        tau_t2t2_L_starstar = tau_t2t2_L_star
                                        tau_t2t2_R_starstar = tau_t2t2_R_star
                                        tau_t1t2_L_starstar = tau_t1t2_L_star
                                        tau_t1t2_R_starstar = tau_t1t2_R_star
                                        E_L_starstar = E_L_star
                                        E_R_starstar = E_R_star
                                    else
                                        tau_tt_L_starstar = tau_tt_L_star + 2._wp*rho_hat*tau_nt_hat/sqrtC_NC*(u_t_star - u_t_L)
                                        tau_tt_R_starstar = tau_tt_R_star - 2._wp*rho_hat*tau_nt_hat/sqrtC_NC*(u_t_star - u_t_R)
                                        tau_t2t2_L_starstar = tau_t2t2_L_star + 2._wp*rho_hat*tau_nt2_hat/sqrtC_NC*(u_t2_star &
                                            & - u_t2_L)
                                        tau_t2t2_R_starstar = tau_t2t2_R_star - 2._wp*rho_hat*tau_nt2_hat/sqrtC_NC*(u_t2_star &
                                            & - u_t2_R)
                                        tau_t1t2_L_starstar = tau_t1t2_L_star + rho_hat*(tau_nt2_hat*(u_t_star - u_t_L) &
                                            & + tau_nt_hat*(u_t2_star - u_t2_L))/sqrtC_NC
                                        tau_t1t2_R_starstar = tau_t1t2_R_star - rho_hat*(tau_nt2_hat*(u_t_star - u_t_R) &
                                            & + tau_nt_hat*(u_t2_star - u_t2_R))/sqrtC_NC
                                        E_L_starstar = E_L_star + (rhoL_star/sqrtC_NC)*((u_t_star*tau_nt_star - u_t_L*tau_nt_L) &
                                                                   & + (u_t2_star*tau_nt2_star - u_t2_L*tau_nt2_L))
                                        E_R_starstar = E_R_star - (rhoR_star/sqrtC_NC)*((u_t_star*tau_nt_star - u_t_R*tau_nt_R) &
                                                                   & + (u_t2_star*tau_nt2_star - u_t2_R*tau_nt2_R))
                                    end if

                                    alpha1_L_star = (alpha_L(1)*(S_L - u_n_L) - C_hat_1*(S_M - u_n_L))/min(S_L - S_M, -verysmall)
                                    alpha1_R_star = (alpha_R(1)*(S_R - u_n_R) - C_hat_1*(S_M - u_n_R))/max(S_R - S_M, verysmall)

                                    alpha2_L_star = (alpha_L(2)*(S_L - u_n_L) - C_hat_2*(S_M - u_n_L))/min(S_L - S_M, -verysmall)
                                    alpha2_R_star = (alpha_R(2)*(S_R - u_n_R) - C_hat_2*(S_M - u_n_R))/max(S_R - S_M, verysmall)

                                    ! HLLD flux, register-diet form: pick the wave-fan zone once (it is
                                    ! component-independent), then fold the selected side's star/starstar states into
                                    ! F_hlld one component at a time through the scalars us_c/uss_c (no fan arrays
                                    ! survive). The L and R sides are mirror images and the per-component statements
                                    ! share one shape, so the whole fold is emitted from one Fypp template driven by
                                    ! HLLD_FAN_SIDES and the HLLD_FAN_COMPS_3D/2D tables above. Per-component operation
                                    ! order matches the materialized form, so the flux is -O0 bit-identical. Do NOT
                                    ! re-expand into per-region temp arrays (or reorder the fold) without re-checking
                                    ! GPU register spill and the -O0 exactness gate.

                                    zone = f_hlld_wave_zone(S_L, S_Lstar, s_M, S_Rstar, S_R)

                                    if (zone == 0) then
                                        F_hlld(1:ncomp) = F_L(1:ncomp)
                                    else if (zone == 5) then
                                        F_hlld(1:ncomp) = F_R(1:ncomp)
                                        #:for S, ZA, ZB, ZSS, SO, SI in HLLD_FAN_SIDES
                                        else if (zone == ${ZA}$ .or. zone == ${ZB}$) then
                                            ! ${S}$ side of the fan: per component, us_c/uss_c are the selected side's
                                            ! star/starstar states (the old F_star${S}$ is folded into the first F_hlld
                                            ! statement; the starstar correction applies in zone ${ZSS}$ only, with the
                                            ! left-associative order of the materialized form preserved)
                                            if (p > 0 .and. .not. cyl_coord) then
                                                #:for CI, COMP in enumerate(HLLD_FAN_COMPS_3D)
                                                    us_c = ${COMP[0].format(S=S)}$
                                                    uss_c = ${'us_c' if COMP[1] is None else COMP[1].format(S=S)}$
                                                    F_hlld(${CI + 1}$) = F_${S}$(${CI + 1}$) + ${SO}$*(us_c - U_${S}$(${CI + 1}$))
                                                    if (zone == ${ZSS}$) F_hlld(${CI + 1}$) = F_hlld(${CI + 1}$) + ${SI}$*(uss_c &
                                                        & - us_c)
                                                #:endfor
                                            else
                                                #:for CI, COMP in enumerate(HLLD_FAN_COMPS_2D)
                                                    us_c = ${COMP[0].format(S=S)}$
                                                    uss_c = ${'us_c' if COMP[1] is None else COMP[1].format(S=S)}$
                                                    F_hlld(${CI + 1}$) = F_${S}$(${CI + 1}$) + ${SO}$*(us_c - U_${S}$(${CI + 1}$))
                                                    if (zone == ${ZSS}$) F_hlld(${CI + 1}$) = F_hlld(${CI + 1}$) + ${SI}$*(uss_c &
                                                        & - us_c)
                                                #:endfor
                                            end if
                                        #:endfor
                                    end if

                                    ! ADC blending (HLLD / HLL)

                                    if (riemann_hypo_ADC) then
                                        ! Register-diet form: the HLL flux enters per component as the scalar F_HLL_c
                                        ! instead of a materialized F_HLL array; outside the subsonic fan F_HLL equals
                                        ! F_hlld and the identity blend is kept explicitly so the arithmetic (including
                                        ! signed-zero behavior) matches the array form bit-for-bit.
                                        if (cyl_coord) then
                                            if (0._wp <= S_L) then
                                                u_n_HLL_trace = u_n_L; u_t_HLL_trace = u_t_L
                                                p_face_HLL = pres%L; tau_qq_face_HLL = tau_qq_L
                                            else if (S_R <= 0._wp) then
                                                u_n_HLL_trace = u_n_R; u_t_HLL_trace = u_t_R
                                                p_face_HLL = pres%R; tau_qq_face_HLL = tau_qq_R
                                            else
                                                u_n_HLL_trace = (S_R*u_n_L - S_L*u_n_R)/(S_R - S_L + verysmall)
                                                u_t_HLL_trace = (S_R*u_t_L - S_L*u_t_R)/(S_R - S_L + verysmall)
                                                ! Only HLL-state components 1, 2, 3, 8 and 11 feed the axisym trace
                                                U_HLL_c = (S_R*U_R(1) - S_L*U_L(1) - (F_R(1) - F_L(1)))/(S_R - S_L + verysmall)
                                                rho_HLL = U_HLL_c
                                                U_HLL_c = (S_R*U_R(2) - S_L*U_L(2) - (F_R(2) - F_L(2)))/(S_R - S_L + verysmall)
                                                rho_HLL = rho_HLL + U_HLL_c
                                                U_HLL_c = (S_R*U_R(3) - S_L*U_L(3) - (F_R(3) - F_L(3)))/(S_R - S_L + verysmall)
                                                u_n_HLL_cons = U_HLL_c/(rho_HLL + verysmall)
                                                U_HLL_c = (S_R*U_R(8) - S_L*U_L(8) - (F_R(8) - F_L(8)))/(S_R - S_L + verysmall)
                                                tau_nn_HLL = U_HLL_c/(rho_HLL + verysmall)
                                                U_HLL_c = (S_R*U_R(11) - S_L*U_L(11) - (F_R(11) - F_L(11)))/(S_R - S_L + verysmall)
                                                tau_qq_face_HLL = U_HLL_c/(rho_HLL + verysmall)
                                                ! This branch implies S_L < 0 < S_R, so component 3 of F_HLL is the
                                                ! interior HLL flux
                                                @:hll_flux_component(F_HLL_c, 3)
                                                p_face_HLL = F_HLL_c - rho_HLL*u_n_HLL_cons*u_n_HLL_cons + tau_nn_HLL
                                            end if
                                        end if

                                        ! phi is anchor-independent: computed once in the shared section above
                                        if (S_L < 0._wp .and. S_R > 0._wp) then
                                            do i = 1, ncomp
                                                @:hll_flux_component(F_HLL_c, i)
                                                F_hlld(i) = F_HLL_c + phi*(F_hlld(i) - F_HLL_c)
                                            end do
                                        else
                                            do i = 1, ncomp
                                                F_HLL_c = F_hlld(i)
                                                F_hlld(i) = F_HLL_c + phi*(F_hlld(i) - F_HLL_c)
                                            end do
                                        end if
                                    end if
                                end if

                                ! Reorder F_HLLD for output: pass 1 (hat_L) -> flux_rs*, pass 2 (hat_R) -> flux_hatR_rs*
                                #:for HATPASS, FLUX in [(1, 'flux_rsx'), (2, 'flux_hatR_rsx')]
                                    if (ipass == ${HATPASS}$) then
                                        if (p > 0 .and. .not. cyl_coord) then
                                            ! 3D Cartesian: 14-state -> physical indices
                                            ${FLUX}$_vf(${SF('')}$, 1) = F_hlld(1)
                                            ${FLUX}$_vf(${SF('')}$, 2) = F_hlld(2)
                                            ${FLUX}$_vf(${SF('')}$, eqn_idx%cont%end + dir_idx(1)) = F_hlld(3)
                                            ${FLUX}$_vf(${SF('')}$, eqn_idx%cont%end + dir_idx(2)) = F_hlld(4)
                                            ${FLUX}$_vf(${SF('')}$, eqn_idx%cont%end + dir_idx(3)) = F_hlld(5)
                                            ${FLUX}$_vf(${SF('')}$, eqn_idx%E) = F_hlld(6)
                                            ${FLUX}$_vf(${SF('')}$, eqn_idx%E + 1) = F_hlld(7)
                                            ${FLUX}$_vf(${SF('')}$, eqn_idx%E + 2) = F_hlld(8)
                                            ! Map local stress to physical stress indices
                                            if (dir_idx(1) == 1) then
                                                ${FLUX}$_vf(${SF('')}$, eqn_idx%stress%beg) = F_hlld(9)  ! tau_nn=tau_xx
                                                ${FLUX}$_vf(${SF('')}$, eqn_idx%stress%beg + 1) = F_hlld(10)  ! tau_nt1=tau_xy
                                                ${FLUX}$_vf(${SF('')}$, eqn_idx%stress%beg + 2) = F_hlld(12)  ! tau_t1t1=tau_yy
                                                ${FLUX}$_vf(${SF('')}$, eqn_idx%stress%beg + 3) = F_hlld(11)  ! tau_nt2=tau_xz
                                                ${FLUX}$_vf(${SF('')}$, eqn_idx%stress%beg + 4) = F_hlld(14)  ! tau_t1t2=tau_yz
                                                ${FLUX}$_vf(${SF('')}$, eqn_idx%stress%beg + 5) = F_hlld(13)  ! tau_t2t2=tau_zz
                                            else if (dir_idx(1) == 2) then
                                                ${FLUX}$_vf(${SF('')}$, eqn_idx%stress%beg) = F_hlld(12)  ! tau_t1t1=tau_xx
                                                ${FLUX}$_vf(${SF('')}$, eqn_idx%stress%beg + 1) = F_hlld(10)  ! tau_nt1=tau_xy
                                                ${FLUX}$_vf(${SF('')}$, eqn_idx%stress%beg + 2) = F_hlld(9)  ! tau_nn=tau_yy
                                                ${FLUX}$_vf(${SF('')}$, eqn_idx%stress%beg + 3) = F_hlld(14)  ! tau_t1t2=tau_xz
                                                ${FLUX}$_vf(${SF('')}$, eqn_idx%stress%beg + 4) = F_hlld(11)  ! tau_nt2=tau_yz
                                                ${FLUX}$_vf(${SF('')}$, eqn_idx%stress%beg + 5) = F_hlld(13)  ! tau_t2t2=tau_zz
                                            else
                                                ${FLUX}$_vf(${SF('')}$, eqn_idx%stress%beg) = F_hlld(12)  ! tau_t1t1=tau_xx
                                                ${FLUX}$_vf(${SF('')}$, eqn_idx%stress%beg + 1) = F_hlld(14)  ! tau_t1t2=tau_xy
                                                ${FLUX}$_vf(${SF('')}$, eqn_idx%stress%beg + 2) = F_hlld(13)  ! tau_t2t2=tau_yy
                                                ${FLUX}$_vf(${SF('')}$, eqn_idx%stress%beg + 3) = F_hlld(10)  ! tau_nt1=tau_xz
                                                ${FLUX}$_vf(${SF('')}$, eqn_idx%stress%beg + 4) = F_hlld(11)  ! tau_nt2=tau_yz
                                                ${FLUX}$_vf(${SF('')}$, eqn_idx%stress%beg + 5) = F_hlld(9)  ! tau_nn=tau_zz
                                            end if
                                        else
                                            ! 2D/axisym: 11-state (unchanged)
                                            ${FLUX}$_vf(${SF('')}$, 1) = F_hlld(1)
                                            ${FLUX}$_vf(${SF('')}$, 2) = F_hlld(2)
                                            ${FLUX}$_vf(${SF('')}$, eqn_idx%cont%end + dir_idx(1)) = F_hlld(3)
                                            ${FLUX}$_vf(${SF('')}$, eqn_idx%cont%end + dir_idx(2)) = F_hlld(4)
                                            ${FLUX}$_vf(${SF('')}$, eqn_idx%E) = F_hlld(5)
                                            ${FLUX}$_vf(${SF('')}$, eqn_idx%E + 1) = F_hlld(6)
                                            ${FLUX}$_vf(${SF('')}$, eqn_idx%E + 2) = F_hlld(7)
                                            ${FLUX}$_vf(${SF('')}$, eqn_idx%stress%beg + 1) = F_hlld(9)
                                            if (dir_idx(1) == 1) then
                                                ${FLUX}$_vf(${SF('')}$, eqn_idx%stress%beg) = F_hlld(8)
                                                ${FLUX}$_vf(${SF('')}$, eqn_idx%stress%beg + 2) = F_hlld(10)
                                            else
                                                ${FLUX}$_vf(${SF('')}$, eqn_idx%stress%beg) = F_hlld(10)
                                                ${FLUX}$_vf(${SF('')}$, eqn_idx%stress%beg + 2) = F_hlld(8)
                                            end if
                                            if (cyl_coord) then
                                                ${FLUX}$_vf(${SF('')}$, eqn_idx%stress%end) = F_hlld(11)
                                            end if
                                        end if
                                    end if
                                #:endfor

                                ! Export face velocities for axisym hypo source terms
                                if (grid_geometry == 2) then
                                    ! Upwind by wave-fan wedge: the inner zones ride the contact (u_n = S_M) and the
                                    ! tangential state switches L -> star -> R across the shear waves (zones 2 and 3
                                    ! share the tangential star state)
                                    if (zone == 0) then
                                        u_n_face = u_n_L; u_t_face = u_t_L
                                    else if (zone == 1) then
                                        u_n_face = S_M; u_t_face = u_t_L
                                    else if (zone <= 3) then
                                        u_n_face = S_M; u_t_face = u_t_star
                                    else if (zone == 4) then
                                        u_n_face = S_M; u_t_face = u_t_R
                                    else
                                        u_n_face = u_n_R; u_t_face = u_t_R
                                    end if
                                    ! ADC blend NC face velocities with HLL scalar traces (non-degenerate only; the degenerate
                                    ! branch does not set phi / HLL traces)
                                    if (riemann_hypo_ADC .and. .not. degenerate) then
                                        u_n_face = u_n_HLL_trace + phi*(u_n_face - u_n_HLL_trace)
                                        u_t_face = u_t_HLL_trace + phi*(u_t_face - u_t_HLL_trace)
                                    end if
                                    #:for HATPASS, NCV in [(1, 'nc_iface_vel_rsx'), (2, 'nc_iface_vel_hatR_rsx')]
                                        if (ipass == ${HATPASS}$) then
                                            if (dir_idx(1) == 1) then
                                                ${NCV}$_vf(${SF('')}$, 1) = u_n_face
                                                ${NCV}$_vf(${SF('')}$, 2) = u_t_face
                                            else
                                                ${NCV}$_vf(${SF('')}$, 1) = u_t_face
                                                ${NCV}$_vf(${SF('')}$, 2) = u_n_face
                                            end if
                                        end if
                                    #:endfor
                                end if

                                ! Radial geometric source flux for cylindrical coordinates
                                #:if (NORM_DIR == 2)
                                    if (cyl_coord) then
                                        #:for HATPASS, GSRC, FLUX in [(1, 'flux_gsrc_rsx', 'flux_rsx'), (2, &
                                                                       & 'flux_gsrc_hatR_rsx', 'flux_hatR_rsx')]
                                            if (ipass == ${HATPASS}$) then
                                                $:GPU_LOOP(parallelism='[seq]')
                                                do i = 1, sys_size
                                                    ${GSRC}$_vf(${SF('')}$, i) = ${FLUX}$_vf(${SF('')}$, i)
                                                end do
                                                ! Pure HLLD face state, upwinded by wave-fan wedge (left/right of the contact)
                                                if (zone == 0) then
                                                    p_face = pres%L; tau_qq_face = tau_qq_L
                                                else if (zone <= 2) then
                                                    p_face = pTot_star + tau_nn_L_star; tau_qq_face = tau_qq_L_star
                                                else if (zone <= 4) then
                                                    p_face = pTot_star + tau_nn_R_star; tau_qq_face = tau_qq_R_star
                                                else
                                                    p_face = pres%R; tau_qq_face = tau_qq_R
                                                end if
                                                ! ADC blend face state (non-degenerate only; the degenerate branch does not set phi
                                                ! / HLL traces)
                                                if (riemann_hypo_ADC .and. .not. degenerate) then
                                                    p_face = p_face_HLL + phi*(p_face - p_face_HLL)
                                                    tau_qq_face = tau_qq_face_HLL + phi*(tau_qq_face - tau_qq_face_HLL)
                                                end if
                                                ${GSRC}$_vf(${SF('')}$, eqn_idx%cont%end + dir_idx(1)) = ${FLUX}$_vf(${SF('')}$, &
                                                            & eqn_idx%cont%end + dir_idx(1)) - p_face + tau_qq_face
                                                $:GPU_LOOP(parallelism='[seq]')
                                                do i = eqn_idx%adv%beg, eqn_idx%adv%end
                                                    ${GSRC}$_vf(${SF('')}$, i) = 0._wp
                                                end do
                                            end if
                                        #:endfor
                                    end if
                                #:endif
                            end do

                            ! Dual-pass HLLD: all NC terms stay inside the Riemann flux (anchor-independent; written once)
                            flux_src_rsx_vf(${SF('')}$, eqn_idx%adv%beg) = 0._wp
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

        call s_finalize_riemann_solver(flux_vf, flux_src_vf, flux_gsrc_vf, norm_dir)

    end subroutine s_hypo_hlld_riemann_solver

    !> Copy the hypo interface velocities of a fused-solve pass from Riemann-space buffers to physical-space output arrays, called
    !! after the Riemann solver for each sweep direction. The plain routine handles the hat_L-anchored values (nc_iface_vel_rsx;
    !! used when hypo_nc_mode_interface is active, and as the hat_L pass of the dual-pass HLLD solve); the _hatR routine handles the
    !! hat_R-anchored values (nc_iface_vel_hatR_rsx) of the dual-pass solve. Both routines are emitted from one template so the two
    !! passes cannot drift apart.
    !! @param nc_iface_vel_vf Output: physical velocity components at interfaces
    !! @param norm_dir Sweep direction (1=x, 2=y, 3=z)
    #:for SUFFIX, NCV in [('', 'nc_iface_vel_rsx'), ('_hatR', 'nc_iface_vel_hatR_rsx')]
        subroutine s_finalize_nc_iface_vel${SUFFIX}$ (nc_iface_vel_vf, norm_dir)

            type(scalar_field), dimension(:), intent(inout) :: nc_iface_vel_vf
            integer, intent(in)                             :: norm_dir
            integer                                         :: i, j, k, l

            if (norm_dir == 2) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = 1, num_dims
                    do l = is3%beg, is3%end
                        do j = is1%beg, is1%end
                            do k = is2%beg, is2%end
                                nc_iface_vel_vf(i)%sf(k, j, l) = ${NCV}$_vf(k, j, l, i)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else if (norm_dir == 1) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = 1, num_dims
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end
                                nc_iface_vel_vf(i)%sf(j, k, l) = ${NCV}$_vf(j, k, l, i)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else if (norm_dir == 3) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = 1, num_dims
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end
                                nc_iface_vel_vf(i)%sf(l, k, j) = ${NCV}$_vf(l, k, j, i)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

        end subroutine s_finalize_nc_iface_vel${SUFFIX}$
    #:endfor
end module m_riemann_solver_hypo_hlld
