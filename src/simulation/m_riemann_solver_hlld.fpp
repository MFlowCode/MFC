!>
!! @file
!! @brief Contains module m_riemann_solver_hlld

!> @brief HLLD approximate Riemann solver for MHD, Miyoshi & Kusano JCP (2005)
#:include 'case.fpp'
#:include 'macros.fpp'

module m_riemann_solver_hlld

    use m_derived_types
    use m_global_parameters
    use m_variables_conversion
    use m_riemann_state

    implicit none

contains

    !> HLLD Riemann solver for MHD, Miyoshi & Kusano JCP (2005)
    subroutine s_hlld_riemann_solver(qL_prim_rsx_vf, dqL_prim_dx_vf, dqL_prim_dy_vf, dqL_prim_dz_vf, qL_prim_vf, qR_prim_rsx_vf, &
                                     & dqR_prim_dx_vf, dqR_prim_dy_vf, dqR_prim_dz_vf, qR_prim_vf, q_prim_vf, flux_vf, &
                                     & flux_src_vf, flux_gsrc_vf, norm_dir, ix, iy, iz)

        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout) :: qL_prim_rsx_vf, qR_prim_rsx_vf
        type(scalar_field), allocatable, dimension(:), intent(inout) :: dqL_prim_dx_vf, dqR_prim_dx_vf, dqL_prim_dy_vf, &
             & dqR_prim_dy_vf, dqL_prim_dz_vf, dqR_prim_dz_vf

        type(scalar_field), allocatable, dimension(:), intent(inout) :: qL_prim_vf, qR_prim_vf
        type(scalar_field), dimension(sys_size), intent(in)          :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout)       :: flux_vf, flux_src_vf, flux_gsrc_vf
        integer, intent(in)                                          :: norm_dir
        type(int_bounds_info), intent(in)                            :: ix, iy, iz

        ! Local variables:

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: alpha_L, alpha_R, alpha_rho_L, alpha_rho_R
        #:else
            real(wp), dimension(num_fluids) :: alpha_L, alpha_R, alpha_rho_L, alpha_rho_R
        #:endif
        type(riemann_states_vec3) :: vel
        type(riemann_states)      :: rho, pres, E, H_no_mag
        type(riemann_states)      :: gamma, pi_inf, qv
        type(riemann_states)      :: vel_rms
        type(riemann_states_vec3) :: B
        type(riemann_states)      :: c, c_fast, pres_mag

        ! HLLD speeds and intermediate state variables:
        real(wp)               :: s_L, s_R, s_M, s_starL, s_starR
        real(wp)               :: pTot_L, pTot_R, p_star, rhoL_star, rhoR_star, E_starL, E_starR
        real(wp), dimension(7) :: U_L, U_R, U_starL, U_starR, U_doubleL, U_doubleR
        real(wp), dimension(7) :: F_L, F_R, F_starL, F_starR, F_hlld

        ! Indices for U and F: (rho, rho*vel(1), rho*vel(2), rho*vel(3), By, Bz, E) Note: vel and B are permutated, so vel(1) is the
        ! normal velocity, and x is the normal direction Note: Bx is omitted as the magnetic flux is always zero in the normal
        ! direction

        real(wp) :: sqrt_rhoL_star, sqrt_rhoR_star, denom_ds, sign_Bx
        real(wp) :: vL_star, vR_star, wL_star, wR_star
        real(wp) :: v_double, w_double, By_double, Bz_double, E_doubleL, E_doubleR, E_double
        integer  :: i, j, k, l

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
                $:GPU_PARALLEL_LOOP(collapse=3, private='[alpha_rho_L, alpha_rho_R, vel, alpha_L, alpha_R, rho, pres, E, &
                                    & H_no_mag, gamma, pi_inf, qv, vel_rms, B, c, c_fast, pres_mag, U_L, U_R, U_starL, U_starR, &
                                    & U_doubleL, U_doubleR, F_L, F_R, F_starL, F_starR, F_hlld, s_L, s_R, s_M, s_starL, s_starR, &
                                    & pTot_L, pTot_R, p_star, rhoL_star, rhoR_star, E_starL, E_starR, sqrt_rhoL_star, &
                                    & sqrt_rhoR_star, denom_ds, sign_Bx, vL_star, vR_star, wL_star, wR_star, v_double, w_double, &
                                    & By_double, Bz_double, E_doubleL, E_doubleR, E_double]', copyin='[norm_dir]')
                do l = ${Z_BND}$%beg, ${Z_BND}$%end
                    do k = ${Y_BND}$%beg, ${Y_BND}$%end
                        do j = ${X_BND}$%beg, ${X_BND}$%end
                            ! (1) Extract the left/right primitive states
                            do i = 1, eqn_idx%cont%end
                                alpha_rho_L(i) = qL_prim_rsx_vf(${SF('')}$, i)
                                alpha_rho_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, i)
                            end do

                            ! NOTE: unlike HLL & HLLC, vel_L here is permutated by dir_idx for simpler logic
                            do i = 1, num_vels
                                vel%L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%cont%end + dir_idx(i))
                                vel%R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%cont%end + dir_idx(i))
                            end do

                            vel_rms%L = sum(vel%L**2._wp)
                            vel_rms%R = sum(vel%R**2._wp)

                            do i = 1, num_fluids
                                alpha_L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%E + i)
                                alpha_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E + i)
                            end do

                            pres%L = qL_prim_rsx_vf(${SF('')}$, eqn_idx%E)
                            pres%R = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E)

                            ! NOTE: unlike HLL, Bx, By, Bz are permutated by dir_idx for simpler logic
                            if (mhd) then
                                if (n == 0) then  ! 1D: constant Bx; By, Bz as variables; only in x so not permutated
                                    B%L = [Bx0, qL_prim_rsx_vf(${SF('')}$, eqn_idx%B%beg), qL_prim_rsx_vf(${SF('')}$, &
                                                               & eqn_idx%B%beg + 1)]
                                    B%R = [Bx0, qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%B%beg), qR_prim_rsx_vf(${SF(' + 1')}$, &
                                                               & eqn_idx%B%beg + 1)]
                                else  ! 2D/3D: Bx, By, Bz as variables
                                    B%L = [qL_prim_rsx_vf(${SF('')}$, eqn_idx%B%beg + dir_idx(1) - 1), qL_prim_rsx_vf(${SF('')}$, &
                                                          & eqn_idx%B%beg + dir_idx(2) - 1), qL_prim_rsx_vf(${SF('')}$, &
                                                          & eqn_idx%B%beg + dir_idx(3) - 1)]
                                    B%R = [qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%B%beg + dir_idx(1) - 1), &
                                                          & qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%B%beg + dir_idx(2) - 1), &
                                                          & qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%B%beg + dir_idx(3) - 1)]
                                end if
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

                            pres_mag%L = 0.5_wp*sum(B%L**2._wp)
                            pres_mag%R = 0.5_wp*sum(B%R**2._wp)
                            E%L = gamma%L*pres%L + pi_inf%L + 0.5_wp*rho%L*vel_rms%L + qv%L + pres_mag%L
                            E%R = gamma%R*pres%R + pi_inf%R + 0.5_wp*rho%R*vel_rms%R + qv%R + pres_mag%R  ! includes magnetic energy
                            H_no_mag%L = (E%L + pres%L - pres_mag%L)/rho%L
                            ! stagnation enthalpy here excludes magnetic energy (only used to find speed of sound)
                            H_no_mag%R = (E%R + pres%R - pres_mag%R)/rho%R

                            ! (2) Compute fast wave speeds
                            call s_compute_speed_of_sound(pres%L, rho%L, gamma%L, pi_inf%L, H_no_mag%L, alpha_L, vel_rms%L, &
                                                          & 0._wp, c%L, qv%L)
                            call s_compute_speed_of_sound(pres%R, rho%R, gamma%R, pi_inf%R, H_no_mag%R, alpha_R, vel_rms%R, &
                                                          & 0._wp, c%R, qv%R)
                            call s_compute_fast_magnetosonic_speed(rho%L, c%L, B%L, norm_dir, c_fast%L, H_no_mag%L)
                            call s_compute_fast_magnetosonic_speed(rho%R, c%R, B%R, norm_dir, c_fast%R, H_no_mag%R)

                            ! (3) Compute contact speed s_M [Miyoshi Equ. (38)]
                            s_L = min(vel%L(1) - c_fast%L, vel%R(1) - c_fast%R)
                            s_R = max(vel%R(1) + c_fast%R, vel%L(1) + c_fast%L)

                            pTot_L = pres%L + pres_mag%L
                            pTot_R = pres%R + pres_mag%R

                            s_M = (((s_R - vel%R(1))*rho%R*vel%R(1) - (s_L - vel%L(1))*rho%L*vel%L(1) - pTot_R + pTot_L)/((s_R &
                                   & - vel%R(1))*rho%R - (s_L - vel%L(1))*rho%L))

                            ! (4) Compute star state variables
                            rhoL_star = rho%L*(s_L - vel%L(1))/(s_L - s_M)
                            rhoR_star = rho%R*(s_R - vel%R(1))/(s_R - s_M)
                            p_star = pTot_L + rho%L*(s_L - vel%L(1))*(s_M - vel%L(1))/(s_L - s_M)
                            E_starL = ((s_L - vel%L(1))*E%L - pTot_L*vel%L(1) + p_star*s_M)/(s_L - s_M)
                            E_starR = ((s_R - vel%R(1))*E%R - pTot_R*vel%R(1) + p_star*s_M)/(s_R - s_M)

                            ! (5) Compute left/right state vectors and fluxes
                            U_L = [rho%L, rho%L*vel%L(1:3), B%L(2:3), E%L]
                            U_starL = [rhoL_star, rhoL_star*s_M, rhoL_star*vel%L(2:3), B%L(2:3), E_starL]
                            U_R = [rho%R, rho%R*vel%R(1:3), B%R(2:3), E%R]
                            U_starR = [rhoR_star, rhoR_star*s_M, rhoR_star*vel%R(2:3), B%R(2:3), E_starR]

                            ! Compute the left/right fluxes
                            F_L(1) = U_L(2)
                            F_L(2) = U_L(2)*vel%L(1) - B%L(1)*B%L(1) + pTot_L
                            F_L(3:4) = U_L(2)*vel%L(2:3) - B%L(1)*B%L(2:3)
                            F_L(5:6) = vel%L(1)*B%L(2:3) - vel%L(2:3)*B%L(1)
                            F_L(7) = (E%L + pTot_L)*vel%L(1) - B%L(1)*(vel%L(1)*B%L(1) + vel%L(2)*B%L(2) + vel%L(3)*B%L(3))

                            F_R(1) = U_R(2)
                            F_R(2) = U_R(2)*vel%R(1) - B%R(1)*B%R(1) + pTot_R
                            F_R(3:4) = U_R(2)*vel%R(2:3) - B%R(1)*B%R(2:3)
                            F_R(5:6) = vel%R(1)*B%R(2:3) - vel%R(2:3)*B%R(1)
                            F_R(7) = (E%R + pTot_R)*vel%R(1) - B%R(1)*(vel%R(1)*B%R(1) + vel%R(2)*B%R(2) + vel%R(3)*B%R(3))
                            ! HLLD star-state fluxes via HLL jump relation
                            F_starL = F_L + s_L*(U_starL - U_L)
                            F_starR = F_R + s_R*(U_starR - U_R)
                            ! Alfven wave speeds bounding the rotational discontinuities
                            s_starL = s_M - abs(B%L(1))/sqrt(rhoL_star)
                            s_starR = s_M + abs(B%L(1))/sqrt(rhoR_star)
                            ! HLLD double-star (intermediate) states across rotational discontinuities
                            sqrt_rhoL_star = sqrt(rhoL_star); sqrt_rhoR_star = sqrt(rhoR_star)
                            vL_star = vel%L(2); wL_star = vel%L(3)
                            vR_star = vel%R(2); wR_star = vel%R(3)

                            ! (6) Compute the double-star states [Miyoshi Eqns. (59)-(62)]
                            denom_ds = sqrt_rhoL_star + sqrt_rhoR_star
                            sign_Bx = sign(1._wp, B%L(1))
                            v_double = (sqrt_rhoL_star*vL_star + sqrt_rhoR_star*vR_star + (B%R(2) - B%L(2))*sign_Bx)/denom_ds
                            w_double = (sqrt_rhoL_star*wL_star + sqrt_rhoR_star*wR_star + (B%R(3) - B%L(3))*sign_Bx)/denom_ds
                            By_double = (sqrt_rhoL_star*B%R(2) + sqrt_rhoR_star*B%L(2) + sqrt_rhoL_star*sqrt_rhoR_star*(vR_star &
                                         & - vL_star)*sign_Bx)/denom_ds
                            Bz_double = (sqrt_rhoL_star*B%R(3) + sqrt_rhoR_star*B%L(3) + sqrt_rhoL_star*sqrt_rhoR_star*(wR_star &
                                         & - wL_star)*sign_Bx)/denom_ds

                            E_doubleL = E_starL - sqrt_rhoL_star*((vL_star*B%L(2) + wL_star*B%L(3)) - (v_double*By_double &
                                                                  & + w_double*Bz_double))*sign_Bx
                            E_doubleR = E_starR + sqrt_rhoR_star*((vR_star*B%R(2) + wR_star*B%R(3)) - (v_double*By_double &
                                                                  & + w_double*Bz_double))*sign_Bx
                            E_double = 0.5_wp*(E_doubleL + E_doubleR)

                            U_doubleL = [rhoL_star, rhoL_star*s_M, rhoL_star*v_double, rhoL_star*w_double, By_double, Bz_double, &
                                & E_double]
                            U_doubleR = [rhoR_star, rhoR_star*s_M, rhoR_star*v_double, rhoR_star*w_double, By_double, Bz_double, &
                                & E_double]

                            ! Select HLLD flux region
                            if (0.0_wp <= s_L) then
                                F_hlld = F_L
                            else if (0.0_wp <= s_starL) then
                                F_hlld = F_L + s_L*(U_starL - U_L)
                            else if (0.0_wp <= s_M) then
                                F_hlld = F_starL + s_starL*(U_doubleL - U_starL)
                            else if (0.0_wp <= s_starR) then
                                F_hlld = F_starR + s_starR*(U_doubleR - U_starR)
                            else if (0.0_wp <= s_R) then
                                F_hlld = F_R + s_R*(U_starR - U_R)
                            else
                                F_hlld = F_R
                            end if

                            ! Hybrid Riemann: overwrite HLLD with a central/Rusanov flux at a WENO-smooth face

                            ! (12) Write HLLD flux to output arrays
                            flux_rsx_vf(${SF('')}$, 1) = F_hlld(1)  ! TODO multi-component
                            ! Momentum
                            flux_rsx_vf(${SF('')}$, eqn_idx%cont%end + dir_idx(1)) = F_hlld(2)
                            flux_rsx_vf(${SF('')}$, eqn_idx%cont%end + dir_idx(2)) = F_hlld(3)
                            flux_rsx_vf(${SF('')}$, eqn_idx%cont%end + dir_idx(3)) = F_hlld(4)
                            ! Magnetic field
                            if (n == 0) then
                                flux_rsx_vf(${SF('')}$, eqn_idx%B%beg) = F_hlld(5)
                                flux_rsx_vf(${SF('')}$, eqn_idx%B%beg + 1) = F_hlld(6)
                            else
                                flux_rsx_vf(${SF('')}$, eqn_idx%B%beg + dir_idx(1) - 1) = 0._wp
                                flux_rsx_vf(${SF('')}$, eqn_idx%B%beg + dir_idx(2) - 1) = F_hlld(5)
                                flux_rsx_vf(${SF('')}$, eqn_idx%B%beg + dir_idx(3) - 1) = F_hlld(6)
                            end if
                            ! Energy
                            flux_rsx_vf(${SF('')}$, eqn_idx%E) = F_hlld(7)
                            ! Volume fractions
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = eqn_idx%adv%beg, eqn_idx%adv%end
                                flux_rsx_vf(${SF('')}$, i) = 0._wp  ! TODO multi-component (zero for now)
                            end do

                            flux_src_rsx_vf(${SF('')}$, eqn_idx%adv%beg) = 0._wp
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

        call s_finalize_riemann_solver(flux_vf, flux_src_vf, flux_gsrc_vf, norm_dir)

    end subroutine s_hlld_riemann_solver

end module m_riemann_solver_hlld
