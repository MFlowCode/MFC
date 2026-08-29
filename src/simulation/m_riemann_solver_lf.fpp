!>
!! @file
!! @brief Contains module m_riemann_solver_lf

!> @brief Lax-Friedrichs (Rusanov) approximate Riemann solver
#:include 'case.fpp'
#:include 'macros.fpp'

module m_riemann_solver_lf

    use m_derived_types
    use m_global_parameters
    use m_variables_conversion
    use m_constants, only: riemann_solver_hll, riemann_solver_hllc, riemann_solver_lax_friedrichs
    use m_thermochem, only: gas_constant, get_mixture_molecular_weight, get_mixture_specific_heat_cv_mass, &
        & get_mixture_energy_mass, get_species_specific_heats_r, get_mixture_specific_heat_cp_mass, molecular_weights
    use m_riemann_state

    implicit none

contains

    !> Lax-Friedrichs (Rusanov) approximate Riemann solver
    subroutine s_lf_riemann_solver(qL_prim_rsx_vf, dqL_prim_dx_vf, dqL_prim_dy_vf, dqL_prim_dz_vf, qL_prim_vf, qR_prim_rsx_vf, &
                                   & dqR_prim_dx_vf, dqR_prim_dy_vf, dqR_prim_dz_vf, qR_prim_vf, q_prim_vf, flux_vf, flux_src_vf, &
                                   & flux_gsrc_vf, norm_dir, ix, iy, iz)

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
            real(wp), dimension(3)    :: alpha_rho_L, alpha_rho_R
            real(wp), dimension(3)    :: vel_L, vel_R
            real(wp), dimension(3)    :: alpha_L, alpha_R
            real(wp), dimension(10)   :: Ys_L, Ys_R
            real(wp), dimension(10)   :: Cp_iL, Cp_iR, Xs_L, Xs_R, Gamma_iL, Gamma_iR
            real(wp), dimension(3, 3) :: vel_grad_L, vel_grad_R  !< Averaged velocity gradient tensor `d(vel_i)/d(coord_j)`.
        #:else
            real(wp), dimension(num_fluids)  :: alpha_rho_L, alpha_rho_R
            real(wp), dimension(num_vels)    :: vel_L, vel_R
            real(wp), dimension(num_fluids)  :: alpha_L, alpha_R
            real(wp), dimension(num_species) :: Ys_L, Ys_R
            real(wp), dimension(num_species) :: Cp_iL, Cp_iR, Xs_L, Xs_R, Gamma_iL, Gamma_iR
            !> Averaged velocity gradient tensor `d(vel_i)/d(coord_j)`.
            real(wp), dimension(num_dims, num_dims) :: vel_grad_L, vel_grad_R
        #:endif
        real(wp) :: rho_L, rho_R
        real(wp) :: pres_L, pres_R
        real(wp) :: E_L, E_R
        real(wp) :: T_L, T_R
        real(wp) :: Y_L, Y_R
        real(wp) :: MW_L, MW_R
        real(wp) :: R_gas_L, R_gas_R
        real(wp) :: Cp_L, Cp_R
        real(wp) :: Cv_L, Cv_R
        real(wp) :: Gamm_L, Gamm_R
        real(wp) :: gamma_L, gamma_R
        real(wp) :: pi_inf_L, pi_inf_R
        real(wp) :: qv_L, qv_R
        real(wp) :: c_L, c_R
        real(wp), dimension(2) :: Re_L, Re_R
        real(wp) :: s_L, s_R, s_M, s_P
        real(wp) :: ptilde_L, ptilde_R
        real(wp) :: vel_L_rms, vel_R_rms
        real(wp) :: alpha_L_sum, alpha_R_sum
        real(wp) :: pcorr  !< low Mach number correction
        integer :: i, j, k, l  !< Generic loop iterators
        integer :: Re_size_loc1, Re_size_loc2  !< host copies of Re_size; amdflang reads the declare-target original stale cross-TU
        integer, dimension(3) :: idx_right_phys  !< Physical (j,k,l) indices for right state.
        ! Populating the buffers of the left and right Riemann problem states variables, based on the choice of boundary conditions

        call s_populate_riemann_states_variables_buffers(qL_prim_rsx_vf, dqL_prim_dx_vf, dqL_prim_dy_vf, dqL_prim_dz_vf, &
            & qR_prim_rsx_vf, dqR_prim_dx_vf, dqR_prim_dy_vf, dqR_prim_dz_vf, norm_dir, ix, iy, iz)

        ! Reshaping inputted data based on dimensional splitting direction
        call s_initialize_riemann_solver(flux_src_vf, norm_dir)
        Re_size_loc1 = Re_size(1); Re_size_loc2 = Re_size(2)
        #:for NORM_DIR, XYZ, STENCIL_VAR, COORDS, X_BND, Y_BND, Z_BND in &
                    [(1, 'x', 'j', '{STENCIL_IDX}, k, l', 'is1', 'is2', 'is3'), &
                     (2, 'y', 'k', 'j, {STENCIL_IDX}, l', 'is2', 'is1', 'is3'), &
                     (3, 'z', 'l', 'j, k, {STENCIL_IDX}', 'is3', 'is2', 'is1')]
            #:set SV = STENCIL_VAR
            #:set SF = lambda offs: COORDS.format(STENCIL_IDX = SV + offs)
            if (norm_dir == ${NORM_DIR}$) then
                $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, alpha_rho_L, alpha_rho_R, vel_L, vel_R, alpha_L, alpha_R, &
                                    & Re_L, Re_R, s_L, s_R, Ys_L, Ys_R, Cp_iL, Cp_iR, Xs_L, Xs_R, Gamma_iL, Gamma_iR, pcorr, &
                                    & vel_grad_L, vel_grad_R, idx_right_phys, vel_L_rms, vel_R_rms, alpha_L_sum, alpha_R_sum, &
                                    & pres_L, pres_R, rho_L, rho_R, gamma_L, gamma_R, pi_inf_L, pi_inf_R, qv_L, qv_R, c_L, c_R, &
                                    & E_L, E_R, ptilde_L, ptilde_R, s_M, s_P, Cp_L, Cp_R, Cv_L, Cv_R, R_gas_L, R_gas_R, MW_L, &
                                    & MW_R, T_L, T_R, Y_L, Y_R]', firstprivate='[Re_size_loc1, Re_size_loc2]')
                do l = ${Z_BND}$%beg, ${Z_BND}$%end
                    do k = ${Y_BND}$%beg, ${Y_BND}$%end
                        do j = ${X_BND}$%beg, ${X_BND}$%end
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, eqn_idx%cont%end
                                alpha_rho_L(i) = qL_prim_rsx_vf(${SF('')}$, i)
                                alpha_rho_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, i)
                            end do

                            vel_L_rms = 0._wp; vel_R_rms = 0._wp

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_vels
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
                                    alpha_rho_L(i) = max(0._wp, alpha_rho_L(i))
                                    alpha_L(i) = min(max(0._wp, alpha_L(i)), 1._wp)
                                    alpha_L_sum = alpha_L_sum + alpha_L(i)
                                    alpha_rho_R(i) = max(0._wp, alpha_rho_R(i))
                                    alpha_R(i) = min(max(0._wp, alpha_R(i)), 1._wp)
                                    alpha_R_sum = alpha_R_sum + alpha_R(i)
                                end do

                                alpha_L = alpha_L/max(alpha_L_sum, sgm_eps)
                                alpha_R = alpha_R/max(alpha_R_sum, sgm_eps)
                            end if

                            call s_compute_mixture_coefficients(alpha_rho_L, alpha_L, rho_L, gamma_L, pi_inf_L, qv_L)
                            call s_compute_mixture_coefficients(alpha_rho_R, alpha_R, rho_R, gamma_R, pi_inf_R, qv_R)

                            if (viscous) then
                                call s_compute_interface_reynolds(alpha_L, Re_L, Re_size_loc1, Re_size_loc2)
                                call s_compute_interface_reynolds(alpha_R, Re_R, Re_size_loc1, Re_size_loc2)
                            end if

                            if (chemistry) then
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
                                    ! gamma_method = 1: Ref. Section 2.3.1 Formulation of doi:10.7907/ZKW8-ES97.
                                    Gamma_iL = Cp_iL/(Cp_iL - 1.0_wp)
                                    Gamma_iR = Cp_iR/(Cp_iR - 1.0_wp)

                                    gamma_L = sum(Xs_L(:)/(Gamma_iL(:) - 1.0_wp))
                                    gamma_R = sum(Xs_R(:)/(Gamma_iR(:) - 1.0_wp))
                                else if (chem_params%gamma_method == 2) then
                                    ! gamma_method = 2: c_p / c_v where c_p, c_v are specific heats.
                                    call get_mixture_specific_heat_cp_mass(T_L, Ys_L, Cp_L)
                                    call get_mixture_specific_heat_cp_mass(T_R, Ys_R, Cp_R)
                                    call get_mixture_specific_heat_cv_mass(T_L, Ys_L, Cv_L)
                                    call get_mixture_specific_heat_cv_mass(T_R, Ys_R, Cv_R)

                                    Gamm_L = Cp_L/Cv_L
                                    gamma_L = 1.0_wp/(Gamm_L - 1.0_wp)
                                    Gamm_R = Cp_R/Cv_R
                                    gamma_R = 1.0_wp/(Gamm_R - 1.0_wp)
                                end if

                                call get_mixture_energy_mass(T_L, Ys_L, E_L)
                                call get_mixture_energy_mass(T_R, Ys_R, E_R)

                                E_L = rho_L*E_L + 5.e-1*rho_L*vel_L_rms
                                E_R = rho_R*E_R + 5.e-1*rho_R*vel_R_rms
                            else
                                call s_compute_energy(pres_L, alpha_rho_L, alpha_L, vel_L_rms, E_L)
                                call s_compute_energy(pres_R, alpha_rho_R, alpha_R, vel_R_rms, E_R)
                            end if

                            call s_compute_speed_of_sound(pres_L, rho_L, gamma_L, pi_inf_L, alpha_L, c_L)

                            call s_compute_speed_of_sound(pres_R, rho_R, gamma_R, pi_inf_R, alpha_R, c_R)

                            s_L = 0._wp; s_R = 0._wp

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_dims
                                s_L = s_L + vel_L(i)**2._wp
                                s_R = s_R + vel_R(i)**2._wp
                            end do

                            s_L = sqrt(s_L)
                            s_R = sqrt(s_R)

                            s_P = max(s_L, s_R) + max(c_L, c_R)
                            s_M = -s_P

                            s_L = s_M
                            s_R = s_P

                            ! Low Mach correction
                            pcorr = f_low_Mach_pcorr_hll(vel_L_rms, vel_R_rms, c_L, c_R, rho_L, rho_R, s_M, s_P)

                            ! Mass
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, eqn_idx%cont%end
                                flux_rsx_vf(${SF('')}$, &
                                            & i) = (s_M*alpha_rho_R(i)*vel_R(norm_dir) - s_P*alpha_rho_L(i)*vel_L(norm_dir) &
                                            & + s_M*s_P*(alpha_rho_L(i) - alpha_rho_R(i)))/(s_M - s_P)
                            end do

                            ! Momentum
                            if (bubbles_euler) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_vels
                                    flux_rsx_vf(${SF('')}$, &
                                                & eqn_idx%cont%end + dir_idx(i)) = (s_M*(rho_R*vel_R(dir_idx(1))*vel_R(dir_idx(i)) &
                                                & + dir_flg(dir_idx(i))*(pres_R - ptilde_R)) - s_P*(rho_L*vel_L(dir_idx(1)) &
                                                & *vel_L(dir_idx(i)) + dir_flg(dir_idx(i))*(pres_L - ptilde_L)) &
                                                & + s_M*s_P*(rho_L*vel_L(dir_idx(i)) - rho_R*vel_R(dir_idx(i))))/(s_M - s_P) &
                                                & + (s_M/s_L)*(s_P/s_R)*pcorr*(vel_R(dir_idx(i)) - vel_L(dir_idx(i)))
                                end do
                            else
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_vels
                                    flux_rsx_vf(${SF('')}$, &
                                                & eqn_idx%cont%end + dir_idx(i)) = (s_M*(rho_R*vel_R(dir_idx(1))*vel_R(dir_idx(i)) &
                                                & + dir_flg(dir_idx(i))*pres_R) - s_P*(rho_L*vel_L(dir_idx(1))*vel_L(dir_idx(i)) &
                                                & + dir_flg(dir_idx(i))*pres_L) + s_M*s_P*(rho_L*vel_L(dir_idx(i)) &
                                                & - rho_R*vel_R(dir_idx(i))))/(s_M - s_P) + (s_M/s_L)*(s_P/s_R) &
                                                & *pcorr*(vel_R(dir_idx(i)) - vel_L(dir_idx(i)))
                                end do
                            end if

                            ! Energy
                            if (bubbles_euler) then
                                flux_rsx_vf(${SF('')}$, &
                                            & eqn_idx%E) = (s_M*vel_R(dir_idx(1))*(E_R + pres_R - ptilde_R) - s_P*vel_L(dir_idx(1) &
                                            & )*(E_L + pres_L - ptilde_L) + s_M*s_P*(E_L - E_R))/(s_M - s_P) + (s_M/s_L)*(s_P/s_R) &
                                            & *pcorr*(vel_R_rms - vel_L_rms)/2._wp
                            else
                                flux_rsx_vf(${SF('')}$, &
                                            & eqn_idx%E) = (s_M*vel_R(dir_idx(1))*(E_R + pres_R) - s_P*vel_L(dir_idx(1))*(E_L &
                                            & + pres_L) + s_M*s_P*(E_L - E_R))/(s_M - s_P) + (s_M/s_L)*(s_P/s_R)*pcorr*(vel_R_rms &
                                            & - vel_L_rms)/2._wp
                            end if

                            ! Advection flux and source: interface velocity for volume fraction transport
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = eqn_idx%adv%beg, eqn_idx%adv%end
                                flux_rsx_vf(${SF('')}$, i) = (qL_prim_rsx_vf(${SF('')}$, i) - qR_prim_rsx_vf(${SF(' + 1')}$, &
                                            & i))*s_M*s_P/(s_M - s_P)
                                flux_src_rsx_vf(${SF('')}$, i) = (s_M*qR_prim_rsx_vf(${SF(' + 1')}$, &
                                                & i) - s_P*qL_prim_rsx_vf(${SF('')}$, i))/(s_M - s_P)
                            end do

                            if (bubbles_euler) then
                                ! From HLLC: Kills mass transport @ bubble gas density
                                if (num_fluids > 1) then
                                    flux_rsx_vf(${SF('')}$, eqn_idx%cont%end) = 0._wp
                                end if
                            end if

                            if (chemistry) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = eqn_idx%species%beg, eqn_idx%species%end
                                    Y_L = qL_prim_rsx_vf(${SF('')}$, i)
                                    Y_R = qR_prim_rsx_vf(${SF(' + 1')}$, i)

                                    flux_rsx_vf(${SF('')}$, &
                                                & i) = (s_M*Y_R*rho_R*vel_R(dir_idx(1)) - s_P*Y_L*rho_L*vel_L(dir_idx(1)) &
                                                & + s_M*s_P*(Y_L*rho_L - Y_R*rho_R))/(s_M - s_P)
                                    flux_src_rsx_vf(${SF('')}$, i) = 0._wp
                                end do
                            end if

                            #:if (NORM_DIR == 2)
                                if (cyl_coord) then
                                    ! Substituting the advective flux into the inviscid geometrical source flux
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, eqn_idx%E
                                        flux_gsrc_rsx_vf(${SF('')}$, i) = flux_rsx_vf(${SF('')}$, i)
                                    end do
                                    ! Recalculating the radial momentum geometric source flux
                                    flux_gsrc_rsx_vf(${SF('')}$, eqn_idx%cont%end + 2) = flux_rsx_vf(${SF('')}$, &
                                                     & eqn_idx%cont%end + 2) - (s_M*pres_R - s_P*pres_L)/(s_M - s_P)
                                    ! Geometrical source of the void fraction(s) is zero
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = eqn_idx%adv%beg, eqn_idx%adv%end
                                        flux_gsrc_rsx_vf(${SF('')}$, i) = flux_rsx_vf(${SF('')}$, i)
                                    end do
                                end if
                            #:endif
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

        if (viscous) then
            $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, idx_right_phys, vel_grad_L, vel_grad_R, alpha_L, alpha_R, &
                                & vel_L, vel_R, Re_L, Re_R]', copyin='[norm_dir]', firstprivate='[Re_size_loc1, Re_size_loc2]')
            do l = isz%beg, isz%end
                do k = isy%beg, isy%end
                    do j = isx%beg, isx%end
                        idx_right_phys(1) = j
                        idx_right_phys(2) = k
                        idx_right_phys(3) = l
                        idx_right_phys(norm_dir) = idx_right_phys(norm_dir) + 1

                        if (norm_dir == 1) then
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                alpha_L(i) = qL_prim_rsx_vf(j, k, l, eqn_idx%E + i)
                                alpha_R(i) = qR_prim_rsx_vf(j + 1, k, l, eqn_idx%E + i)
                            end do

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_dims
                                vel_L(i) = qL_prim_rsx_vf(j, k, l, eqn_idx%mom%beg + i - 1)
                                vel_R(i) = qR_prim_rsx_vf(j + 1, k, l, eqn_idx%mom%beg + i - 1)
                            end do
                        else if (norm_dir == 2) then
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                alpha_L(i) = qL_prim_rsx_vf(j, k, l, eqn_idx%E + i)
                                alpha_R(i) = qR_prim_rsx_vf(j, k + 1, l, eqn_idx%E + i)
                            end do
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_dims
                                vel_L(i) = qL_prim_rsx_vf(j, k, l, eqn_idx%mom%beg + i - 1)
                                vel_R(i) = qR_prim_rsx_vf(j, k + 1, l, eqn_idx%mom%beg + i - 1)
                            end do
                        else
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                alpha_L(i) = qL_prim_rsx_vf(j, k, l, eqn_idx%E + i)
                                alpha_R(i) = qR_prim_rsx_vf(j, k, l + 1, eqn_idx%E + i)
                            end do

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_dims
                                vel_L(i) = qL_prim_rsx_vf(j, k, l, eqn_idx%mom%beg + i - 1)
                                vel_R(i) = qR_prim_rsx_vf(j, k, l + 1, eqn_idx%mom%beg + i - 1)
                            end do
                        end if

                        call s_compute_interface_reynolds(alpha_L, Re_L, Re_size_loc1, Re_size_loc2)
                        call s_compute_interface_reynolds(alpha_R, Re_R, Re_size_loc1, Re_size_loc2)

                        if (shear_stress) then
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_dims
                                vel_grad_L(i, 1) = (dqL_prim_dx_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l)/Re_L(1))
                                vel_grad_R(i, 1) = (dqR_prim_dx_vf(eqn_idx%mom%beg + i - 1)%sf(idx_right_phys(1), &
                                           & idx_right_phys(2), idx_right_phys(3))/Re_R(1))
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
                                    if (num_dims > 1) then
                                        vel_grad_L(i, 2) = (dqL_prim_dy_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l)/Re_L(1))
                                        vel_grad_R(i, 2) = (dqR_prim_dy_vf(eqn_idx%mom%beg + i - 1)%sf(idx_right_phys(1), &
                                                   & idx_right_phys(2), idx_right_phys(3))/Re_R(1))
                                    end if
                                    #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                                        if (num_dims > 2) then
                                            vel_grad_L(i, 3) = (dqL_prim_dz_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l)/Re_L(1))
                                            vel_grad_R(i, 3) = (dqR_prim_dz_vf(eqn_idx%mom%beg + i - 1)%sf(idx_right_phys(1), &
                                                       & idx_right_phys(2), idx_right_phys(3))/Re_R(1))
                                        end if
                                    #:endif
                                #:endif
                            end do

                            if (norm_dir == 1) then
                                flux_src_vf(eqn_idx%mom%beg)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg)%sf(j, k, &
                                            & l) - (4._wp/3._wp)*0.5_wp*(vel_grad_L(1, 1) + vel_grad_R(1, 1))
                                flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                            & l) - (4._wp/3._wp)*0.5_wp*(vel_grad_L(1, 1)*vel_L(1) + vel_grad_R(1, 1)*vel_R(1))
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
                                    if (num_dims > 1) then
                                        flux_src_vf(eqn_idx%mom%beg)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg)%sf(j, k, &
                                                    & l) - (-2._wp/3._wp)*0.5_wp*(vel_grad_L(2, 2) + vel_grad_R(2, 2))
                                        flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                    & l) - (-2._wp/3._wp)*0.5_wp*(vel_grad_L(2, 2)*vel_L(1) + vel_grad_R(2, &
                                                    & 2)*vel_R(1))

                                        flux_src_vf(eqn_idx%mom%beg + 1)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + 1)%sf(j, k, &
                                                    & l) - 0.5_wp*(vel_grad_L(1, 2) + vel_grad_R(1, 2)) - 0.5_wp*(vel_grad_L(2, &
                                                    & 1) + vel_grad_R(2, 1))
                                        flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                    & l) - 0.5_wp*(vel_grad_L(1, 2)*vel_L(2) + vel_grad_R(1, &
                                                    & 2)*vel_R(2)) - 0.5_wp*(vel_grad_L(2, 1)*vel_L(2) + vel_grad_R(2, 1)*vel_R(2))
                                        #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                                            if (num_dims > 2) then
                                                flux_src_vf(eqn_idx%mom%beg)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg)%sf(j, k, &
                                                            & l) - (-2._wp/3._wp)*0.5_wp*(vel_grad_L(3, 3) + vel_grad_R(3, 3))
                                                flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                            & l) - (-2._wp/3._wp)*0.5_wp*(vel_grad_L(3, &
                                                            & 3)*vel_L(1) + vel_grad_R(3, 3)*vel_R(1))

                                                flux_src_vf(eqn_idx%mom%beg + 2)%sf(j, k, &
                                                            & l) = flux_src_vf(eqn_idx%mom%beg + 2)%sf(j, k, &
                                                            & l) - 0.5_wp*(vel_grad_L(1, 3) + vel_grad_R(1, &
                                                            & 3)) - 0.5_wp*(vel_grad_L(3, 1) + vel_grad_R(3, 1))
                                                flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                            & l) - 0.5_wp*(vel_grad_L(1, 3)*vel_L(3) + vel_grad_R(1, &
                                                            & 3)*vel_R(3)) - 0.5_wp*(vel_grad_L(3, 1)*vel_L(3) + vel_grad_R(3, &
                                                            & 1)*vel_R(3))
                                            end if
                                        #:endif
                                    end if
                                #:endif
                            else if (norm_dir == 2) then
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
                                    flux_src_vf(eqn_idx%mom%beg + 1)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + 1)%sf(j, k, &
                                                & l) - (-2._wp/3._wp)*0.5_wp*(vel_grad_L(1, 1) + vel_grad_R(1, 1))
                                    flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                & l) - (-2._wp/3._wp)*0.5_wp*(vel_grad_L(1, 1)*vel_L(2) + vel_grad_R(1, 1)*vel_R(2))

                                    flux_src_vf(eqn_idx%mom%beg + 1)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + 1)%sf(j, k, &
                                                & l) - (4._wp/3._wp)*0.5_wp*(vel_grad_L(2, 2) + vel_grad_R(2, 2))
                                    flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                & l) - (4._wp/3._wp)*0.5_wp*(vel_grad_L(2, 2)*vel_L(2) + vel_grad_R(2, 2)*vel_R(2))

                                    flux_src_vf(eqn_idx%mom%beg)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg)%sf(j, k, &
                                                & l) - 0.5_wp*(vel_grad_L(1, 2) + vel_grad_R(1, 2)) - 0.5_wp*(vel_grad_L(2, &
                                                & 1) + vel_grad_R(2, 1))
                                    flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                & l) - 0.5_wp*(vel_grad_L(1, 2)*vel_L(1) + vel_grad_R(1, &
                                                & 2)*vel_R(1)) - 0.5_wp*(vel_grad_L(2, 1)*vel_L(1) + vel_grad_R(2, 1)*vel_R(1))
                                    #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                                        if (num_dims > 2) then
                                            flux_src_vf(eqn_idx%mom%beg + 1)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + 1)%sf(j, &
                                                        & k, l) - (-2._wp/3._wp)*0.5_wp*(vel_grad_L(3, 3) + vel_grad_R(3, 3))
                                            flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                        & l) - (-2._wp/3._wp)*0.5_wp*(vel_grad_L(3, 3)*vel_L(2) + vel_grad_R(3, &
                                                        & 3)*vel_R(2))

                                            flux_src_vf(eqn_idx%mom%beg + 2)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + 2)%sf(j, &
                                                        & k, l) - 0.5_wp*(vel_grad_L(2, 3) + vel_grad_R(2, &
                                                        & 3)) - 0.5_wp*(vel_grad_L(3, 2) + vel_grad_R(3, 2))
                                            flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                        & l) - 0.5_wp*(vel_grad_L(2, 3)*vel_L(3) + vel_grad_R(2, &
                                                        & 3)*vel_R(3)) - 0.5_wp*(vel_grad_L(3, 2)*vel_L(3) + vel_grad_R(3, &
                                                        & 2)*vel_R(3))
                                        end if
                                    #:endif
                                #:endif
                            else
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                                    flux_src_vf(eqn_idx%mom%beg + 2)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + 2)%sf(j, k, &
                                                & l) - (-2._wp/3._wp)*0.5_wp*(vel_grad_L(1, 1) + vel_grad_R(1, 1))
                                    flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                & l) - (-2._wp/3._wp)*0.5_wp*(vel_grad_L(1, 1)*vel_L(3) + vel_grad_R(1, 1)*vel_R(3))

                                    flux_src_vf(eqn_idx%mom%beg + 2)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + 2)%sf(j, k, &
                                                & l) - (-2._wp/3._wp)*0.5_wp*(vel_grad_L(2, 2) + vel_grad_R(2, 2))
                                    flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                & l) - (-2._wp/3._wp)*0.5_wp*(vel_grad_L(2, 2)*vel_L(3) + vel_grad_R(2, 2)*vel_R(3))

                                    flux_src_vf(eqn_idx%mom%beg)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg)%sf(j, k, &
                                                & l) - 0.5_wp*(vel_grad_L(1, 3) + vel_grad_R(1, 3)) - 0.5_wp*(vel_grad_L(3, &
                                                & 1) + vel_grad_R(3, 1))
                                    flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                & l) - 0.5_wp*(vel_grad_L(1, 3)*vel_L(1) + vel_grad_R(1, &
                                                & 3)*vel_R(1)) - 0.5_wp*(vel_grad_L(3, 1)*vel_L(1) + vel_grad_R(3, 1)*vel_R(1))

                                    flux_src_vf(eqn_idx%mom%beg + 2)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + 2)%sf(j, k, &
                                                & l) - (4._wp/3._wp)*0.5_wp*(vel_grad_L(3, 3) + vel_grad_R(3, 3))
                                    flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                & l) - (4._wp/3._wp)*0.5_wp*(vel_grad_L(3, 3)*vel_L(3) + vel_grad_R(3, 3)*vel_R(3))

                                    flux_src_vf(eqn_idx%mom%beg + 1)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + 1)%sf(j, k, &
                                                & l) - 0.5_wp*(vel_grad_L(2, 3) + vel_grad_R(2, 3)) - 0.5_wp*(vel_grad_L(3, &
                                                & 2) + vel_grad_R(3, 2))
                                    flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                & l) - 0.5_wp*(vel_grad_L(2, 3)*vel_L(2) + vel_grad_R(2, &
                                                & 3)*vel_R(2)) - 0.5_wp*(vel_grad_L(3, 2)*vel_L(2) + vel_grad_R(3, 2)*vel_R(2))
                                #:endif
                            end if
                        end if

                        if (bulk_stress) then
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_dims
                                vel_grad_L(i, 1) = (dqL_prim_dx_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l)/Re_L(2))
                                vel_grad_R(i, 1) = (dqR_prim_dx_vf(eqn_idx%mom%beg + i - 1)%sf(idx_right_phys(1), &
                                           & idx_right_phys(2), idx_right_phys(3))/Re_R(2))
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
                                    if (num_dims > 1) then
                                        vel_grad_L(i, 2) = (dqL_prim_dy_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l)/Re_L(2))
                                        vel_grad_R(i, 2) = (dqR_prim_dy_vf(eqn_idx%mom%beg + i - 1)%sf(idx_right_phys(1), &
                                                   & idx_right_phys(2), idx_right_phys(3))/Re_R(2))
                                    end if
                                #:endif
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                                    if (num_dims > 2) then
                                        vel_grad_L(i, 3) = (dqL_prim_dz_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l)/Re_L(2))
                                        vel_grad_R(i, 3) = (dqR_prim_dz_vf(eqn_idx%mom%beg + i - 1)%sf(idx_right_phys(1), &
                                                   & idx_right_phys(2), idx_right_phys(3))/Re_R(2))
                                    end if
                                #:endif
                            end do

                            if (norm_dir == 1) then
                                flux_src_vf(eqn_idx%mom%beg)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg)%sf(j, k, &
                                            & l) - 0.5_wp*(vel_grad_L(1, 1) + vel_grad_R(1, 1))
                                flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, l) - 0.5_wp*(vel_grad_L(1, &
                                            & 1)*vel_L(1) + vel_grad_R(1, 1)*vel_R(1))
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
                                    if (num_dims > 1) then
                                        flux_src_vf(eqn_idx%mom%beg)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg)%sf(j, k, &
                                                    & l) - 0.5_wp*(vel_grad_L(2, 2) + vel_grad_R(2, 2))
                                        flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                    & l) - 0.5_wp*(vel_grad_L(2, 2)*vel_L(1) + vel_grad_R(2, 2)*vel_R(1))

                                        #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                                            if (num_dims > 2) then
                                                flux_src_vf(eqn_idx%mom%beg)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg)%sf(j, k, &
                                                            & l) - 0.5_wp*(vel_grad_L(3, 3) + vel_grad_R(3, 3))
                                                flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                            & l) - 0.5_wp*(vel_grad_L(3, 3)*vel_L(1) + vel_grad_R(3, 3)*vel_R(1))
                                            end if
                                        #:endif
                                    end if
                                #:endif
                            else if (norm_dir == 2) then
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
                                    flux_src_vf(eqn_idx%mom%beg + 1)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + 1)%sf(j, k, &
                                                & l) - 0.5_wp*(vel_grad_L(1, 1) + vel_grad_R(1, 1))
                                    flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                & l) - 0.5_wp*(vel_grad_L(1, 1)*vel_L(2) + vel_grad_R(1, 1)*vel_R(2))

                                    flux_src_vf(eqn_idx%mom%beg + 1)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + 1)%sf(j, k, &
                                                & l) - 0.5_wp*(vel_grad_L(2, 2) + vel_grad_R(2, 2))
                                    flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                & l) - 0.5_wp*(vel_grad_L(2, 2)*vel_L(2) + vel_grad_R(2, 2)*vel_R(2))

                                    #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                                        if (num_dims > 2) then
                                            flux_src_vf(eqn_idx%mom%beg + 1)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + 1)%sf(j, &
                                                        & k, l) - 0.5_wp*(vel_grad_L(3, 3) + vel_grad_R(3, 3))
                                            flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                        & l) - 0.5_wp*(vel_grad_L(3, 3)*vel_L(2) + vel_grad_R(3, 3)*vel_R(2))
                                        end if
                                    #:endif
                                #:endif
                            else
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                                    flux_src_vf(eqn_idx%mom%beg + 2)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + 2)%sf(j, k, &
                                                & l) - 0.5_wp*(vel_grad_L(1, 1) + vel_grad_R(1, 1))
                                    flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                & l) - 0.5_wp*(vel_grad_L(1, 1)*vel_L(3) + vel_grad_R(1, 1)*vel_R(3))

                                    flux_src_vf(eqn_idx%mom%beg + 2)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + 2)%sf(j, k, &
                                                & l) - 0.5_wp*(vel_grad_L(2, 2) + vel_grad_R(2, 2))
                                    flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                & l) - 0.5_wp*(vel_grad_L(2, 2)*vel_L(3) + vel_grad_R(2, 2)*vel_R(3))

                                    flux_src_vf(eqn_idx%mom%beg + 2)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + 2)%sf(j, k, &
                                                & l) - 0.5_wp*(vel_grad_L(3, 3) + vel_grad_R(3, 3))
                                    flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                & l) - 0.5_wp*(vel_grad_L(3, 3)*vel_L(3) + vel_grad_R(3, 3)*vel_R(3))
                                #:endif
                            end if
                        end if
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        call s_finalize_riemann_solver(flux_vf, flux_src_vf, flux_gsrc_vf, norm_dir)

    end subroutine s_lf_riemann_solver

end module m_riemann_solver_lf
