!>
!! @file
!! @brief Contains module m_riemann_solvers

!> @brief Approximate and exact Riemann solvers (HLL, HLLC, HLLD, exact) for the multicomponent Navier--Stokes equations

#:include 'case.fpp'
#:include 'macros.fpp'
#:include 'inline_riemann.fpp'

module m_riemann_solvers

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_variables_conversion
    use m_bubbles
    use m_bubbles_EE
    use m_surface_tension
    use m_helper_basic
    use m_chemistry
    use m_thermochem, only: gas_constant, get_mixture_molecular_weight, get_mixture_specific_heat_cv_mass, &
        & get_mixture_energy_mass, get_species_specific_heats_r, get_species_enthalpies_rt, get_mixture_specific_heat_cp_mass

    implicit none

    private; public :: s_initialize_riemann_solvers_module, s_riemann_solver, s_hll_riemann_solver, s_hllc_riemann_solver, &
        & s_hlld_riemann_solver, s_lf_riemann_solver, s_finalize_riemann_solvers_module

    !> The cell-boundary values of the fluxes (src - source) that are computed through the chosen Riemann problem solver, and the
    !! direct evaluation of source terms, by using the left and right states given in qK_prim_rs_vf, dqK_prim_ds_vf where ds = dx,
    !! dy or dz.
    !> @{
    real(wp), allocatable, dimension(:,:,:,:) :: flux_rsx_vf, flux_src_rsx_vf
    $:GPU_DECLARE(create='[flux_rsx_vf, flux_src_rsx_vf]')
    !> @}

    !> The cell-boundary values of the geometrical source flux that are computed through the chosen Riemann problem solver by using
    !! the left and right states given in qK_prim_rs_vf. Currently 2D axisymmetric for inviscid only.
    !> @{
    real(wp), allocatable, dimension(:,:,:,:) :: flux_gsrc_rsx_vf
    $:GPU_DECLARE(create='[flux_gsrc_rsx_vf]')
    !> @}

    ! Cell-boundary velocity from Riemann solution; used for source flux

    real(wp), allocatable, dimension(:,:,:,:) :: vel_src_rsx_vf
    $:GPU_DECLARE(create='[vel_src_rsx_vf]')

    real(wp), allocatable, dimension(:,:,:,:) :: mom_sp_rsx_vf
    $:GPU_DECLARE(create='[mom_sp_rsx_vf]')

    real(wp), allocatable, dimension(:,:,:,:) :: Re_avg_rsx_vf
    $:GPU_DECLARE(create='[Re_avg_rsx_vf]')

    !> @name Indical bounds in the s1-, s2- and s3-directions
    !> @{
    type(int_bounds_info) :: is1, is2, is3
    type(int_bounds_info) :: isx, isy, isz
    !> @}

    $:GPU_DECLARE(create='[is1, is2, is3, isx, isy, isz]')

    real(wp), allocatable, dimension(:) :: Gs_rs
    $:GPU_DECLARE(create='[Gs_rs]')

    real(wp), allocatable, dimension(:,:) :: Res_gs
    $:GPU_DECLARE(create='[Res_gs]')

contains

    !> Dispatch to the subroutines that are utilized to compute the Riemann problem solution. For additional information please
    !! reference: 1) s_hll_riemann_solver 2) s_hllc_riemann_solver 3) s_lf_riemann_solver 4) s_hlld_riemann_solver
    subroutine s_riemann_solver(qL_prim_rsx_vf, dqL_prim_dx_vf, dqL_prim_dy_vf, dqL_prim_dz_vf, qL_prim_vf, qR_prim_rsx_vf, &
                                & dqR_prim_dx_vf, dqR_prim_dy_vf, dqR_prim_dz_vf, qR_prim_vf, q_prim_vf, flux_vf, flux_src_vf, &
                                & flux_gsrc_vf, norm_dir, ix, iy, iz)

        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout) :: qL_prim_rsx_vf, qR_prim_rsx_vf
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), allocatable, dimension(:), intent(inout) :: qL_prim_vf, qR_prim_vf
        type(scalar_field), allocatable, dimension(:), intent(inout) :: dqL_prim_dx_vf, dqR_prim_dx_vf, dqL_prim_dy_vf, &
             & dqR_prim_dy_vf, dqL_prim_dz_vf, dqR_prim_dz_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: flux_vf, flux_src_vf, flux_gsrc_vf
        integer, intent(in)                                    :: norm_dir
        type(int_bounds_info), intent(in)                      :: ix, iy, iz

        #:for NAME, NUM in [('hll', 1), ('hllc', 2), ('hlld', 4), ('lf', 5)]
            if (riemann_solver == ${NUM}$) then
                call s_${NAME}$_riemann_solver(qL_prim_rsx_vf, dqL_prim_dx_vf, dqL_prim_dy_vf, dqL_prim_dz_vf, qL_prim_vf, &
                                               & qR_prim_rsx_vf, dqR_prim_dx_vf, dqR_prim_dy_vf, dqR_prim_dz_vf, qR_prim_vf, &
                                               & q_prim_vf, flux_vf, flux_src_vf, flux_gsrc_vf, norm_dir, ix, iy, iz)
            end if
        #:endfor

    end subroutine s_riemann_solver

    !> Dispatch to the subroutines that are utilized to compute the viscous source fluxes for either Cartesian or cylindrical
    !! geometries. For more information please refer to: 1) s_compute_cartesian_viscous_source_flux 2)
    !! s_compute_cylindrical_viscous_source_flux
    subroutine s_compute_viscous_source_flux(velL_vf, dvelL_dx_vf, dvelL_dy_vf, dvelL_dz_vf, velR_vf, dvelR_dx_vf, dvelR_dy_vf, &

        & dvelR_dz_vf, flux_src_vf, norm_dir, ix, iy, iz)

        type(scalar_field), dimension(num_vels), intent(in) :: velL_vf, velR_vf, dvelL_dx_vf, dvelR_dx_vf, dvelL_dy_vf, &
             & dvelR_dy_vf, dvelL_dz_vf, dvelR_dz_vf

        type(scalar_field), dimension(sys_size), intent(inout) :: flux_src_vf
        integer, intent(in)                                    :: norm_dir
        type(int_bounds_info), intent(in)                      :: ix, iy, iz

        if (grid_geometry == 3) then
            call s_compute_cylindrical_viscous_source_flux(velL_vf, dvelL_dx_vf, dvelL_dy_vf, dvelL_dz_vf, velR_vf, dvelR_dx_vf, &
                & dvelR_dy_vf, dvelR_dz_vf, flux_src_vf, norm_dir, ix, iy, iz)
        else
            call s_compute_cartesian_viscous_source_flux(dvelL_dx_vf, dvelL_dy_vf, dvelL_dz_vf, dvelR_dx_vf, dvelR_dy_vf, &
                & dvelR_dz_vf, flux_src_vf, norm_dir)
        end if

    end subroutine s_compute_viscous_source_flux

    !> HLL approximate Riemann solver, Harten et al. SIAM Review (1983)
    subroutine s_hll_riemann_solver(qL_prim_rsx_vf, dqL_prim_dx_vf, dqL_prim_dy_vf, dqL_prim_dz_vf, qL_prim_vf, qR_prim_rsx_vf, &
                                    & dqR_prim_dx_vf, dqR_prim_dy_vf, dqR_prim_dz_vf, qR_prim_vf, q_prim_vf, flux_vf, &
                                    & flux_src_vf, flux_gsrc_vf, norm_dir, ix, iy, iz)

        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout) :: qL_prim_rsx_vf, qR_prim_rsx_vf
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), allocatable, dimension(:), intent(inout) :: qL_prim_vf, qR_prim_vf
        type(scalar_field), allocatable, dimension(:), intent(inout) :: dqL_prim_dx_vf, dqR_prim_dx_vf, dqL_prim_dy_vf, &
             & dqR_prim_dy_vf, dqL_prim_dz_vf, dqR_prim_dz_vf

        ! Intercell fluxes
        type(scalar_field), dimension(sys_size), intent(inout) :: flux_vf, flux_src_vf, flux_gsrc_vf
        type(riemann_states)                                   :: flux_tau
        integer, intent(in)                                    :: norm_dir
        type(int_bounds_info), intent(in)                      :: ix, iy, iz

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3)  :: alpha_rho_L, alpha_rho_R
            real(wp), dimension(3)  :: alpha_L, alpha_R
            real(wp), dimension(10) :: Ys_L, Ys_R
            real(wp), dimension(10) :: Cp_iL, Cp_iR, Xs_L, Xs_R, Gamma_iL, Gamma_iR
            real(wp), dimension(10) :: Yi_avg, Phi_avg, h_iL, h_iR, h_avg_2
        #:else
            real(wp), dimension(num_fluids)  :: alpha_rho_L, alpha_rho_R
            real(wp), dimension(num_fluids)  :: alpha_L, alpha_R
            real(wp), dimension(num_species) :: Ys_L, Ys_R
            real(wp), dimension(num_species) :: Cp_iL, Cp_iR, Xs_L, Xs_R, Gamma_iL, Gamma_iR
            real(wp), dimension(num_species) :: Yi_avg, Phi_avg, h_iL, h_iR, h_avg_2
        #:endif
        type(riemann_states_vec3) :: vel
        type(riemann_states)      :: rho, pres, E, H
        type(riemann_states)      :: T, Y, MW, R_gas, Cp, Cv, Gamm, gamma
        type(riemann_states)      :: pi_inf, qv, c, G, ptilde
        type(riemann_states)      :: vel_rms, vel_tmp
        real(wp)                  :: Cp_avg, Cv_avg, T_avg, eps, c_sum_Yi_Phi
        type(riemann_states_arr6) :: tau_e
        type(riemann_states_arr2) :: Re
        type(riemann_states_vec3) :: xi_field
        real(wp)                  :: rho_avg
        real(wp)                  :: H_avg
        real(wp)                  :: qv_avg
        real(wp)                  :: gamma_avg
        real(wp)                  :: c_avg
        real(wp)                  :: vel_avg_rms
        real(wp)                  :: s_M, s_P, s_S
        type(riemann_states)      :: s
        real(wp)                  :: xi_M, xi_P
        type(riemann_states)      :: Ms, pres_S
        real(wp)                  :: alpha_L_sum, alpha_R_sum
        real(wp)                  :: zcoef, pcorr   !< low Mach number correction
        type(riemann_states)      :: c_fast, pres_mag
        type(riemann_states_vec3) :: B
        type(riemann_states)      :: Ga             !< Gamma (Lorentz factor)
        type(riemann_states)      :: vdotB, B2
        type(riemann_states_vec3) :: b4             !< 4-magnetic field components (spatial: b4x, b4y, b4z)
        type(riemann_states_vec3) :: cm             !< Conservative momentum variables
        integer                   :: i, j, k, l, q  !< Generic loop iterators
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
                $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, q, alpha_rho_L, alpha_rho_R, vel, alpha_L, alpha_R, tau_e, &
                                    & Re, s, s_S, Ys_L, Ys_R, xi_field, Cp_iL, Cp_iR, Xs_L, Xs_R, Gamma_iL, Gamma_iR, Yi_avg, &
                                    & Phi_avg, h_iL, h_iR, h_avg_2, c_fast, pres_mag, B, Ga, vdotB, B2, b4, cm, pcorr, zcoef, &
                                    & vel_tmp, rho, pres, E, H, Cp_avg, Cv_avg, T_avg, eps, c_sum_Yi_Phi, T, Y, MW, R_gas, Cp, &
                                    & Cv, Gamm, gamma, pi_inf, qv, qv_avg, c, G, rho_avg, H_avg, c_avg, gamma_avg, ptilde, &
                                    & vel_rms, vel_avg_rms, Ms, pres_S, alpha_L_sum, alpha_R_sum, flux_tau]', copyin='[norm_dir]')
                do l = ${Z_BND}$%beg, ${Z_BND}$%end
                    do k = ${Y_BND}$%beg, ${Y_BND}$%end
                        do j = ${X_BND}$%beg, ${X_BND}$%end
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, eqn_idx%cont%end
                                alpha_rho_L(i) = qL_prim_rsx_vf(${SF('')}$, i)
                                alpha_rho_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, i)
                            end do

                            vel_rms%L = 0._wp; vel_rms%R = 0._wp

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_vels
                                vel%L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%cont%end + i)
                                vel%R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%cont%end + i)
                                vel_rms%L = vel_rms%L + vel%L(i)**2._wp
                                vel_rms%R = vel_rms%R + vel%R(i)**2._wp
                            end do

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                alpha_L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%E + i)
                                alpha_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E + i)
                            end do

                            pres%L = qL_prim_rsx_vf(${SF('')}$, eqn_idx%E)
                            pres%R = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E)

                            if (mhd) then
                                if (n == 0) then  ! 1D: constant Bx; By, Bz as variables
                                    B%L(1) = Bx0
                                    B%R(1) = Bx0
                                    B%L(2) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%B%beg)
                                    B%R(2) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%B%beg)
                                    B%L(3) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%B%beg + 1)
                                    B%R(3) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%B%beg + 1)
                                else  ! 2D/3D: Bx, By, Bz as variables
                                    B%L(1) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%B%beg)
                                    B%R(1) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%B%beg)
                                    B%L(2) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%B%beg + 1)
                                    B%R(2) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%B%beg + 1)
                                    B%L(3) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%B%beg + 2)
                                    B%R(3) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%B%beg + 2)
                                end if
                            end if

                            rho%L = 0._wp
                            gamma%L = 0._wp
                            pi_inf%L = 0._wp
                            qv%L = 0._wp

                            rho%R = 0._wp
                            gamma%R = 0._wp
                            pi_inf%R = 0._wp
                            qv%R = 0._wp

                            alpha_L_sum = 0._wp
                            alpha_R_sum = 0._wp

                            pres_mag%L = 0._wp
                            pres_mag%R = 0._wp

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

                            if (viscous) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, 2
                                    Re%L(i) = dflt_real
                                    Re%R(i) = dflt_real

                                    if (Re_size(i) > 0) Re%L(i) = 0._wp
                                    if (Re_size(i) > 0) Re%R(i) = 0._wp

                                    $:GPU_LOOP(parallelism='[seq]')
                                    do q = 1, Re_size(i)
                                        Re%L(i) = alpha_L(Re_idx(i, q))/Res_gs(i, q) + Re%L(i)
                                        Re%R(i) = alpha_R(Re_idx(i, q))/Res_gs(i, q) + Re%R(i)
                                    end do

                                    Re%L(i) = 1._wp/max(Re%L(i), sgm_eps)
                                    Re%R(i) = 1._wp/max(Re%R(i), sgm_eps)
                                end do
                            end if

                            if (chemistry) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = eqn_idx%species%beg, eqn_idx%species%end
                                    Ys_L(i - eqn_idx%species%beg + 1) = qL_prim_rsx_vf(${SF('')}$, i)
                                    Ys_R(i - eqn_idx%species%beg + 1) = qR_prim_rsx_vf(${SF(' + 1')}$, i)
                                end do

                                call get_mixture_molecular_weight(Ys_L, MW%L)
                                call get_mixture_molecular_weight(Ys_R, MW%R)
                                Xs_L(:) = Ys_L(:)*MW%L/molecular_weights(:)
                                Xs_R(:) = Ys_R(:)*MW%R/molecular_weights(:)

                                R_gas%L = gas_constant/MW%L
                                R_gas%R = gas_constant/MW%R
                                T%L = pres%L/rho%L/R_gas%L
                                T%R = pres%R/rho%R/R_gas%R

                                call get_species_specific_heats_r(T%L, Cp_iL)
                                call get_species_specific_heats_r(T%R, Cp_iR)

                                if (chem_params%gamma_method == 1) then
                                    ! gamma_method = 1: Ref. Section 2.3.1 Formulation of doi:10.7907/ZKW8-ES97.
                                    Gamma_iL = Cp_iL/(Cp_iL - 1.0_wp)
                                    Gamma_iR = Cp_iR/(Cp_iR - 1.0_wp)

                                    gamma%L = sum(Xs_L(:)/(Gamma_iL(:) - 1.0_wp))
                                    gamma%R = sum(Xs_R(:)/(Gamma_iR(:) - 1.0_wp))
                                else if (chem_params%gamma_method == 2) then
                                    ! gamma_method = 2: c_p / c_v where c_p, c_v are specific heats.
                                    call get_mixture_specific_heat_cp_mass(T%L, Ys_L, Cp%L)
                                    call get_mixture_specific_heat_cp_mass(T%R, Ys_R, Cp%R)
                                    call get_mixture_specific_heat_cv_mass(T%L, Ys_L, Cv%L)
                                    call get_mixture_specific_heat_cv_mass(T%R, Ys_R, Cv%R)

                                    Gamm%L = Cp%L/Cv%L
                                    gamma%L = 1.0_wp/(Gamm%L - 1.0_wp)
                                    Gamm%R = Cp%R/Cv%R
                                    gamma%R = 1.0_wp/(Gamm%R - 1.0_wp)
                                end if

                                call get_mixture_energy_mass(T%L, Ys_L, E%L)
                                call get_mixture_energy_mass(T%R, Ys_R, E%R)

                                E%L = rho%L*E%L + 5.e-1*rho%L*vel_rms%L
                                E%R = rho%R*E%R + 5.e-1*rho%R*vel_rms%R
                                H%L = (E%L + pres%L)/rho%L
                                H%R = (E%R + pres%R)/rho%R
                            else if (mhd .and. relativity) then
                                Ga%L = 1._wp/sqrt(1._wp - vel_rms%L)
                                Ga%R = 1._wp/sqrt(1._wp - vel_rms%R)
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                                    vdotB%L = vel%L(1)*B%L(1) + vel%L(2)*B%L(2) + vel%L(3)*B%L(3)
                                    vdotB%R = vel%R(1)*B%R(1) + vel%R(2)*B%R(2) + vel%R(3)*B%R(3)

                                    b4%L(1:3) = B%L(1:3)/Ga%L + Ga%L*vel%L(1:3)*vdotB%L
                                    b4%R(1:3) = B%R(1:3)/Ga%R + Ga%R*vel%R(1:3)*vdotB%R
                                    B2%L = B%L(1)**2._wp + B%L(2)**2._wp + B%L(3)**2._wp
                                    B2%R = B%R(1)**2._wp + B%R(2)**2._wp + B%R(3)**2._wp
                                #:endif

                                pres_mag%L = 0.5_wp*(B2%L/Ga%L**2._wp + vdotB%L**2._wp)
                                pres_mag%R = 0.5_wp*(B2%R/Ga%R**2._wp + vdotB%R**2._wp)

                                ! Hard-coded EOS
                                H%L = 1._wp + (gamma%L + 1)*pres%L/rho%L
                                H%R = 1._wp + (gamma%R + 1)*pres%R/rho%R
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                                    cm%L(1:3) = (rho%L*H%L*Ga%L**2 + B2%L)*vel%L(1:3) - vdotB%L*B%L(1:3)
                                    cm%R(1:3) = (rho%R*H%R*Ga%R**2 + B2%R)*vel%R(1:3) - vdotB%R*B%R(1:3)
                                #:endif

                                E%L = rho%L*H%L*Ga%L**2 - pres%L + 0.5_wp*(B2%L + vel_rms%L*B2%L - vdotB%L**2._wp) - rho%L*Ga%L
                                E%R = rho%R*H%R*Ga%R**2 - pres%R + 0.5_wp*(B2%R + vel_rms%R*B2%R - vdotB%R**2._wp) - rho%R*Ga%R
                            else if (mhd .and. .not. relativity) then
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                                    pres_mag%L = 0.5_wp*(B%L(1)**2._wp + B%L(2)**2._wp + B%L(3)**2._wp)
                                    pres_mag%R = 0.5_wp*(B%R(1)**2._wp + B%R(2)**2._wp + B%R(3)**2._wp)
                                #:endif
                                E%L = gamma%L*pres%L + pi_inf%L + 0.5_wp*rho%L*vel_rms%L + qv%L + pres_mag%L
                                ! includes magnetic energy
                                E%R = gamma%R*pres%R + pi_inf%R + 0.5_wp*rho%R*vel_rms%R + qv%R + pres_mag%R
                                H%L = (E%L + pres%L - pres_mag%L)/rho%L
                                ! stagnation enthalpy here excludes magnetic energy (only used to find speed of sound)
                                H%R = (E%R + pres%R - pres_mag%R)/rho%R
                            else
                                E%L = gamma%L*pres%L + pi_inf%L + 5.e-1*rho%L*vel_rms%L + qv%L
                                E%R = gamma%R*pres%R + pi_inf%R + 5.e-1*rho%R*vel_rms%R + qv%R
                                H%L = (E%L + pres%L)/rho%L
                                H%R = (E%R + pres%R)/rho%R
                            end if

                            ! elastic energy update
                            if (hypoelasticity) then
                                G%L = 0._wp; G%R = 0._wp

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    G%L = G%L + alpha_L(i)*Gs_rs(i)
                                    G%R = G%R + alpha_R(i)*Gs_rs(i)
                                end do

                                if (cont_damage) then
                                    G%L = G%L*max((1._wp - qL_prim_rsx_vf(${SF('')}$, eqn_idx%damage)), 0._wp)
                                    G%R = G%R*max((1._wp - qR_prim_rsx_vf(${SF('')}$, eqn_idx%damage)), 0._wp)
                                end if

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, eqn_idx%stress%end - eqn_idx%stress%beg + 1
                                    tau_e%L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%stress%beg - 1 + i)
                                    tau_e%R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%stress%beg - 1 + i)
                                    ! Elastic contribution to energy if G large enough TODO take out if statement if stable without
                                    if ((G%L > 1000) .and. (G%R > 1000)) then
                                        E%L = E%L + (tau_e%L(i)*tau_e%L(i))/(4._wp*G%L)
                                        E%R = E%R + (tau_e%R(i)*tau_e%R(i))/(4._wp*G%R)
                                        ! Double for shear stresses
                                        if (any(eqn_idx%stress%beg - 1 + i == shear_indices)) then
                                            E%L = E%L + (tau_e%L(i)*tau_e%L(i))/(4._wp*G%L)
                                            E%R = E%R + (tau_e%R(i)*tau_e%R(i))/(4._wp*G%R)
                                        end if
                                    end if
                                end do
                            end if

                            @:compute_average_state()

                            call s_compute_speed_of_sound(pres%L, rho%L, gamma%L, pi_inf%L, H%L, alpha_L, vel_rms%L, 0._wp, c%L, &
                                                          & qv%L)

                            call s_compute_speed_of_sound(pres%R, rho%R, gamma%R, pi_inf%R, H%R, alpha_R, vel_rms%R, 0._wp, c%R, &
                                                          & qv%R)

                            !> The computation of c_avg does not require all the variables, and therefore the non '_avg'
                            ! variables are placeholders to call the subroutine.

                            call s_compute_speed_of_sound(pres%R, rho_avg, gamma_avg, pi_inf%R, H_avg, alpha_R, vel_avg_rms, &
                                                          & c_sum_Yi_Phi, c_avg, qv_avg)

                            if (mhd) then
                                call s_compute_fast_magnetosonic_speed(rho%L, c%L, B%L, norm_dir, c_fast%L, H%L)
                                call s_compute_fast_magnetosonic_speed(rho%R, c%R, B%R, norm_dir, c_fast%R, H%R)
                            end if

                            if (viscous) then
                                if (chemistry) then
                                    call compute_viscosity_and_inversion(T%L, Ys_L, T%R, Ys_R, Re%L(1), Re%R(1))
                                end if
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, 2
                                    Re_avg_rsx_vf(${SF('')}$, i) = 2._wp/(1._wp/Re%L(i) + 1._wp/Re%R(i))
                                end do
                            end if

                            ! Wave speed estimates (wave_speeds=1: direct, wave_speeds=2: pressure-based)
                            if (wave_speeds == 1) then
                                if (mhd) then
                                    ! MHD: use fast magnetosonic speed
                                    s%L = min(vel%L(dir_idx(1)) - c_fast%L, vel%R(dir_idx(1)) - c_fast%R)
                                    s%R = max(vel%R(dir_idx(1)) + c_fast%R, vel%L(dir_idx(1)) + c_fast%L)
                                else if (hypoelasticity) then
                                    ! Elastic wave speed, Rodriguez et al. JCP (2019)
                                    s%L = min(vel%L(dir_idx(1)) - sqrt(c%L*c%L + (((4._wp*G%L)/3._wp) + tau_e%L(dir_idx_tau(1))) &
                                              & /rho%L), &
                                              & vel%R(dir_idx(1)) - sqrt(c%R*c%R + (((4._wp*G%R)/3._wp) + tau_e%R(dir_idx_tau(1))) &
                                              & /rho%R))
                                    s%R = max(vel%R(dir_idx(1)) + sqrt(c%R*c%R + (((4._wp*G%R)/3._wp) + tau_e%R(dir_idx_tau(1))) &
                                              & /rho%R), &
                                              & vel%L(dir_idx(1)) + sqrt(c%L*c%L + (((4._wp*G%L)/3._wp) + tau_e%L(dir_idx_tau(1))) &
                                              & /rho%L))
                                else if (hyperelasticity) then
                                    s%L = min(vel%L(dir_idx(1)) - sqrt(c%L*c%L + (4._wp*G%L/3._wp)/rho%L), &
                                              & vel%R(dir_idx(1)) - sqrt(c%R*c%R + (4._wp*G%R/3._wp)/rho%R))
                                    s%R = max(vel%R(dir_idx(1)) + sqrt(c%R*c%R + (4._wp*G%R/3._wp)/rho%R), &
                                              & vel%L(dir_idx(1)) + sqrt(c%L*c%L + (4._wp*G%L/3._wp)/rho%L))
                                else
                                    s%L = min(vel%L(dir_idx(1)) - c%L, vel%R(dir_idx(1)) - c%R)
                                    s%R = max(vel%R(dir_idx(1)) + c%R, vel%L(dir_idx(1)) + c%L)
                                end if

                                if (hyper_cleaning) then
                                    ! Dedner GLM divergence cleaning, Dedner et al. JCP (2002)
                                    s%L = min(s%L, -hyper_cleaning_speed)
                                    s%R = max(s%R, hyper_cleaning_speed)
                                end if

                                s_S = (pres%R - pres%L + rho%L*vel%L(dir_idx(1))*(s%L - vel%L(dir_idx(1))) &
                                       & - rho%R*vel%R(dir_idx(1))*(s%R - vel%R(dir_idx(1))))/(rho%L*(s%L - vel%L(dir_idx(1))) &
                                       & - rho%R*(s%R - vel%R(dir_idx(1))))
                            else if (wave_speeds == 2) then
                                pres_S%L = 5.e-1_wp*(pres%L + pres%R + rho_avg*c_avg*(vel%L(dir_idx(1)) - vel%R(dir_idx(1))))

                                pres_S%R = pres_S%L

                                ! Low Mach correction: Thornber et al. JCP (2008)
                                Ms%L = max(1._wp, &
                                           & sqrt(1._wp + ((5.e-1_wp + gamma%L)/(1._wp + gamma%L))*(pres_S%L/pres%L - 1._wp) &
                                           & *pres%L/((pres%L + pi_inf%L/(1._wp + gamma%L)))))
                                Ms%R = max(1._wp, &
                                           & sqrt(1._wp + ((5.e-1_wp + gamma%R)/(1._wp + gamma%R))*(pres_S%R/pres%R - 1._wp) &
                                           & *pres%R/((pres%R + pi_inf%R/(1._wp + gamma%R)))))

                                s%L = vel%L(dir_idx(1)) - c%L*Ms%L
                                s%R = vel%R(dir_idx(1)) + c%R*Ms%R

                                s_S = 5.e-1_wp*((vel%L(dir_idx(1)) + vel%R(dir_idx(1))) + (pres%L - pres%R)/(rho_avg*c_avg))
                            end if

                            s_M = min(0._wp, s%L); s_P = max(0._wp, s%R)

                            xi_M = (5.e-1_wp + sign(5.e-1_wp, s%L)) + (5.e-1_wp - sign(5.e-1_wp, s%L))*(5.e-1_wp + sign(5.e-1_wp, &
                                    & s%R))
                            xi_P = (5.e-1_wp - sign(5.e-1_wp, s%R)) + (5.e-1_wp - sign(5.e-1_wp, s%L))*(5.e-1_wp + sign(5.e-1_wp, &
                                    & s%R))

                            ! HLL intercell flux: F* = (s%R*F%L - s%L*F%R + s%L*s%R*(U%R - U%L)) / (s%R - s%L) Low Mach correction
                            if (low_Mach == 1) then
                                @:compute_low_Mach_correction()
                            else
                                pcorr = 0._wp
                            end if

                            ! Mass
                            if (.not. relativity) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, eqn_idx%cont%end
                                    flux_rsx_vf(${SF('')}$, &
                                                & i) = (s_M*alpha_rho_R(i)*vel%R(norm_dir) - s_P*alpha_rho_L(i)*vel%L(norm_dir) &
                                                & + s_M*s_P*(alpha_rho_L(i) - alpha_rho_R(i)))/(s_M - s_P)
                                end do
                            else if (relativity) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, eqn_idx%cont%end
                                    flux_rsx_vf(${SF('')}$, &
                                                & i) = (s_M*Ga%R*alpha_rho_R(i)*vel%R(norm_dir) - s_P*Ga%L*alpha_rho_L(i) &
                                                & *vel%L(norm_dir) + s_M*s_P*(Ga%L*alpha_rho_L(i) - Ga%R*alpha_rho_R(i)))/(s_M &
                                                & - s_P)
                                end do
                            end if

                            ! Momentum
                            if (mhd .and. (.not. relativity)) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, 3
                                    ! Flux of rho*v_i in the ${XYZ}$ direction = rho * v_i * v_${XYZ}$ - B_i * B_${XYZ}$ +
                                    ! delta_(${XYZ}$,i) * p_tot
                                    flux_rsx_vf(${SF('')}$, &
                                                & eqn_idx%cont%end + i) = (s_M*(rho%R*vel%R(i)*vel%R(norm_dir) - B%R(i) &
                                                & *B%R(norm_dir) + dir_flg(i)*(pres%R + pres_mag%R)) - s_P*(rho%L*vel%L(i) &
                                                & *vel%L(norm_dir) - B%L(i)*B%L(norm_dir) + dir_flg(i)*(pres%L + pres_mag%L)) &
                                                & + s_M*s_P*(rho%L*vel%L(i) - rho%R*vel%R(i)))/(s_M - s_P)
                                end do
                            else if (mhd .and. relativity) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, 3
                                    ! Flux of m_i in the ${XYZ}$ direction = m_i * v_${XYZ}$ - b_i/Gamma * B_${XYZ}$ +
                                    ! delta_(${XYZ}$,i) * p_tot
                                    flux_rsx_vf(${SF('')}$, &
                                                & eqn_idx%cont%end + i) = (s_M*(cm%R(i)*vel%R(norm_dir) - b4%R(i) &
                                                & /Ga%R*B%R(norm_dir) + dir_flg(i)*(pres%R + pres_mag%R)) - s_P*(cm%L(i) &
                                                & *vel%L(norm_dir) - b4%L(i)/Ga%L*B%L(norm_dir) + dir_flg(i)*(pres%L + pres_mag%L) &
                                                & ) + s_M*s_P*(cm%L(i) - cm%R(i)))/(s_M - s_P)
                                end do
                            else if (bubbles_euler) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_vels
                                    flux_rsx_vf(${SF('')}$, &
                                                & eqn_idx%cont%end + dir_idx(i)) = (s_M*(rho%R*vel%R(dir_idx(1))*vel%R(dir_idx(i)) &
                                                & + dir_flg(dir_idx(i))*(pres%R - ptilde%R)) - s_P*(rho%L*vel%L(dir_idx(1)) &
                                                & *vel%L(dir_idx(i)) + dir_flg(dir_idx(i))*(pres%L - ptilde%L)) &
                                                & + s_M*s_P*(rho%L*vel%L(dir_idx(i)) - rho%R*vel%R(dir_idx(i))))/(s_M - s_P) &
                                                & + (s_M/s%L)*(s_P/s%R)*pcorr*(vel%R(dir_idx(i)) - vel%L(dir_idx(i)))
                                end do
                            else if (hypoelasticity) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_vels
                                    flux_rsx_vf(${SF('')}$, &
                                                & eqn_idx%cont%end + dir_idx(i)) = (s_M*(rho%R*vel%R(dir_idx(1))*vel%R(dir_idx(i)) &
                                                & + dir_flg(dir_idx(i))*pres%R - tau_e%R(dir_idx_tau(i))) &
                                                & - s_P*(rho%L*vel%L(dir_idx(1))*vel%L(dir_idx(i)) + dir_flg(dir_idx(i))*pres%L &
                                                & - tau_e%L(dir_idx_tau(i))) + s_M*s_P*(rho%L*vel%L(dir_idx(i)) &
                                                & - rho%R*vel%R(dir_idx(i))))/(s_M - s_P)
                                end do
                            else
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_vels
                                    flux_rsx_vf(${SF('')}$, &
                                                & eqn_idx%cont%end + dir_idx(i)) = (s_M*(rho%R*vel%R(dir_idx(1))*vel%R(dir_idx(i)) &
                                                & + dir_flg(dir_idx(i))*pres%R) - s_P*(rho%L*vel%L(dir_idx(1))*vel%L(dir_idx(i)) &
                                                & + dir_flg(dir_idx(i))*pres%L) + s_M*s_P*(rho%L*vel%L(dir_idx(i)) &
                                                & - rho%R*vel%R(dir_idx(i))))/(s_M - s_P) + (s_M/s%L)*(s_P/s%R) &
                                                & *pcorr*(vel%R(dir_idx(i)) - vel%L(dir_idx(i)))
                                end do
                            end if

                            ! Energy
                            if (mhd .and. (.not. relativity)) then
                                ! energy flux = (E + p + p_mag) * v_${XYZ}$ - B_${XYZ}$ * (v_x*B_x + v_y*B_y + v_z*B_z)
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                                    flux_rsx_vf(${SF('')}$, &
                                                & eqn_idx%E) = (s_M*(vel%R(norm_dir)*(E%R + pres%R + pres_mag%R) - B%R(norm_dir) &
                                                & *(vel%R(1)*B%R(1) + vel%R(2)*B%R(2) + vel%R(3)*B%R(3))) - s_P*(vel%L(norm_dir) &
                                                & *(E%L + pres%L + pres_mag%L) - B%L(norm_dir)*(vel%L(1)*B%L(1) + vel%L(2)*B%L(2) &
                                                & + vel%L(3)*B%L(3))) + s_M*s_P*(E%L - E%R))/(s_M - s_P)
                                #:endif
                            else if (mhd .and. relativity) then
                                ! energy flux = m_${XYZ}$ - mass flux Hard-coded for single-component for now
                                flux_rsx_vf(${SF('')}$, &
                                            & eqn_idx%E) = (s_M*(cm%R(norm_dir) - Ga%R*alpha_rho_R(1)*vel%R(norm_dir)) &
                                            & - s_P*(cm%L(norm_dir) - Ga%L*alpha_rho_L(1)*vel%L(norm_dir)) + s_M*s_P*(E%L - E%R)) &
                                            & /(s_M - s_P)
                            else if (bubbles_euler) then
                                flux_rsx_vf(${SF('')}$, &
                                            & eqn_idx%E) = (s_M*vel%R(dir_idx(1))*(E%R + pres%R - ptilde%R) - s_P*vel%L(dir_idx(1) &
                                            & )*(E%L + pres%L - ptilde%L) + s_M*s_P*(E%L - E%R))/(s_M - s_P) + (s_M/s%L)*(s_P/s%R) &
                                            & *pcorr*(vel_rms%R - vel_rms%L)/2._wp
                            else if (hypoelasticity) then
                                flux_tau%L = 0._wp; flux_tau%R = 0._wp
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    flux_tau%L = flux_tau%L + tau_e%L(dir_idx_tau(i))*vel%L(dir_idx(i))
                                    flux_tau%R = flux_tau%R + tau_e%R(dir_idx_tau(i))*vel%R(dir_idx(i))
                                end do
                                flux_rsx_vf(${SF('')}$, &
                                            & eqn_idx%E) = (s_M*(vel%R(dir_idx(1))*(E%R + pres%R) - flux_tau%R) &
                                            & - s_P*(vel%L(dir_idx(1))*(E%L + pres%L) - flux_tau%L) + s_M*s_P*(E%L - E%R))/(s_M &
                                            & - s_P)
                            else
                                flux_rsx_vf(${SF('')}$, &
                                            & eqn_idx%E) = (s_M*vel%R(dir_idx(1))*(E%R + pres%R) - s_P*vel%L(dir_idx(1))*(E%L &
                                            & + pres%L) + s_M*s_P*(E%L - E%R))/(s_M - s_P) + (s_M/s%L)*(s_P/s%R)*pcorr*(vel_rms%R &
                                            & - vel_rms%L)/2._wp
                            end if

                            ! Elastic Stresses
                            if (hypoelasticity) then
                                do i = 1, eqn_idx%stress%end - eqn_idx%stress%beg + 1  ! TODO: this indexing may be slow
                                    flux_rsx_vf(${SF('')}$, &
                                                & eqn_idx%stress%beg - 1 + i) = (s_M*(rho%R*vel%R(dir_idx(1))*tau_e%R(i)) &
                                                & - s_P*(rho%L*vel%L(dir_idx(1))*tau_e%L(i)) + s_M*s_P*(rho%L*tau_e%L(i) &
                                                & - rho%R*tau_e%R(i)))/(s_M - s_P)
                                end do
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
                                    Y%L = qL_prim_rsx_vf(${SF('')}$, i)
                                    Y%R = qR_prim_rsx_vf(${SF(' + 1')}$, i)

                                    flux_rsx_vf(${SF('')}$, &
                                                & i) = (s_M*Y%R*rho%R*vel%R(dir_idx(1)) - s_P*Y%L*rho%L*vel%L(dir_idx(1)) &
                                                & + s_M*s_P*(Y%L*rho%L - Y%R*rho%R))/(s_M - s_P)
                                    flux_src_rsx_vf(${SF('')}$, i) = 0._wp
                                end do
                            end if

                            ! MHD: magnetic flux and Maxwell stress contributions
                            if (mhd) then
                                if (n == 0) then  ! 1D: d/dx flux only & Bx = Bx0 = const.
                                    ! B_y flux = v_x * B_y - v_y * Bx0 B_z flux = v_x * B_z - v_z * Bx0
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 0, 1
                                        flux_rsx_vf(j, k, l, &
                                                    & eqn_idx%B%beg + i) = (s_M*(vel%R(1)*B%R(2 + i) - vel%R(2 + i)*Bx0) &
                                                    & - s_P*(vel%L(1)*B%L(2 + i) - vel%L(2 + i)*Bx0) + s_M*s_P*(B%L(2 + i) &
                                                    & - B%R(2 + i)))/(s_M - s_P)
                                    end do
                                else  ! 2D/3D: Bx, By, Bz /= const. but zero flux component in the same direction
                                    ! B_x d/d${XYZ}$ flux = (1 - delta(x,${XYZ}$)) * (v_${XYZ}$ * B_x - v_x * B_${XYZ}$) B_y
                                    ! d/d${XYZ}$ flux = (1 - delta(y,${XYZ}$)) * (v_${XYZ}$ * B_y - v_y * B_${XYZ}$) B_z d/d${XYZ}$
                                    ! flux = (1 - delta(z,${XYZ}$)) * (v_${XYZ}$ * B_z - v_z * B_${XYZ}$)
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 0, 2
                                        flux_rsx_vf(${SF('')}$, &
                                                    & eqn_idx%B%beg + i) = (s_M*(vel%R(dir_idx(1))*B%R(i + 1) - vel%R(i + 1) &
                                                    & *B%R(norm_dir)) - s_P*(vel%L(dir_idx(1))*B%L(i + 1) - vel%L(i + 1) &
                                                    & *B%L(norm_dir)) + s_M*s_P*(B%L(i + 1) - B%R(i + 1)))/(s_M - s_P)
                                    end do

                                    if (hyper_cleaning) then
                                        ! propagate magnetic field divergence as a wave
                                        flux_rsx_vf(${SF('')}$, eqn_idx%B%beg + norm_dir - 1) = flux_rsx_vf(${SF('')}$, &
                                                    & eqn_idx%B%beg + norm_dir - 1) + (s_M*qR_prim_rsx_vf(${SF(' + 1')}$, &
                                                    & eqn_idx%psi) - s_P*qL_prim_rsx_vf(${SF('')}$, eqn_idx%psi))/(s_M - s_P)

                                        flux_rsx_vf(${SF('')}$, &
                                                    & eqn_idx%psi) = (hyper_cleaning_speed**2*(s_M*B%R(norm_dir) &
                                                    & - s_P*B%L(norm_dir)) + s_M*s_P*(qL_prim_rsx_vf(${SF('')}$, &
                                                    & eqn_idx%psi) - qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%psi)))/(s_M - s_P)
                                    else
                                        ! Without hyperbolic cleaning, make sure flux of B_normal is identically zero
                                        flux_rsx_vf(${SF('')}$, eqn_idx%B%beg + norm_dir - 1) = 0._wp
                                    end if
                                end if
                                flux_src_rsx_vf(${SF('')}$, eqn_idx%adv%beg) = 0._wp
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
                                                     & eqn_idx%cont%end + 2) - (s_M*pres%R - s_P*pres%L)/(s_M - s_P)
                                    ! Geometrical source of the void fraction(s) is zero
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = eqn_idx%adv%beg, eqn_idx%adv%end
                                        flux_gsrc_rsx_vf(${SF('')}$, i) = flux_rsx_vf(${SF('')}$, i)
                                    end do
                                end if

                                if (cyl_coord .and. hypoelasticity) then
                                    ! += tau_sigmasigma using HLL
                                    flux_gsrc_rsx_vf(${SF('')}$, eqn_idx%cont%end + 2) = flux_gsrc_rsx_vf(${SF('')}$, &
                                                     & eqn_idx%cont%end + 2) + (s_M*tau_e%R(4) - s_P*tau_e%L(4))/(s_M - s_P)

                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = eqn_idx%stress%beg, eqn_idx%stress%end
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
            if (weno_Re_flux) then
                call s_compute_viscous_source_flux(qL_prim_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                                                   & dqL_prim_dx_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                                                   & dqL_prim_dy_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                                                   & dqL_prim_dz_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                                                   & qR_prim_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                                                   & dqR_prim_dx_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                                                   & dqR_prim_dy_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                                                   & dqR_prim_dz_vf(eqn_idx%mom%beg:eqn_idx%mom%end), flux_src_vf, norm_dir, ix, &
                                                   & iy, iz)
            else
                call s_compute_viscous_source_flux(q_prim_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                                                   & dqL_prim_dx_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                                                   & dqL_prim_dy_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                                                   & dqL_prim_dz_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                                                   & q_prim_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                                                   & dqR_prim_dx_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                                                   & dqR_prim_dy_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                                                   & dqR_prim_dz_vf(eqn_idx%mom%beg:eqn_idx%mom%end), flux_src_vf, norm_dir, ix, &
                                                   & iy, iz)
            end if
        end if

        call s_finalize_riemann_solver(flux_vf, flux_src_vf, flux_gsrc_vf, norm_dir)

    end subroutine s_hll_riemann_solver

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
        type(riemann_states)                                   :: flux_tau
        integer, intent(in)                                    :: norm_dir
        type(int_bounds_info), intent(in)                      :: ix, iy, iz

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3)    :: alpha_rho_L, alpha_rho_R
            real(wp), dimension(3)    :: alpha_L, alpha_R
            real(wp), dimension(10)   :: Ys_L, Ys_R
            real(wp), dimension(10)   :: Cp_iL, Cp_iR, Xs_L, Xs_R, Gamma_iL, Gamma_iR
            real(wp), dimension(10)   :: Yi_avg, Phi_avg, h_iL, h_iR, h_avg_2
            real(wp), dimension(3, 3) :: vel_grad_L, vel_grad_R  !< Averaged velocity gradient tensor `d(vel_i)/d(coord_j)`.
        #:else
            real(wp), dimension(num_fluids)  :: alpha_rho_L, alpha_rho_R
            real(wp), dimension(num_fluids)  :: alpha_L, alpha_R
            real(wp), dimension(num_species) :: Ys_L, Ys_R
            real(wp), dimension(num_species) :: Cp_iL, Cp_iR, Xs_L, Xs_R, Gamma_iL, Gamma_iR
            real(wp), dimension(num_species) :: Yi_avg, Phi_avg, h_iL, h_iR, h_avg_2
            !> Averaged velocity gradient tensor `d(vel_i)/d(coord_j)`.
            real(wp), dimension(num_dims, num_dims) :: vel_grad_L, vel_grad_R
        #:endif
        type(riemann_states_vec3) :: vel
        type(riemann_states)      :: rho, pres, E, H
        type(riemann_states)      :: T, Y, MW, R_gas, Cp, Cv, Gamm, gamma
        type(riemann_states)      :: pi_inf, qv, c, G, ptilde
        type(riemann_states)      :: vel_rms, vel_tmp
        real(wp)                  :: Cp_avg, Cv_avg, T_avg, eps, c_sum_Yi_Phi
        type(riemann_states_arr6) :: tau_e
        type(riemann_states_arr2) :: Re
        type(riemann_states_vec3) :: xi_field
        real(wp)                  :: rho_avg
        real(wp)                  :: H_avg
        real(wp)                  :: qv_avg
        real(wp)                  :: gamma_avg
        real(wp)                  :: c_avg
        real(wp)                  :: vel_avg_rms
        real(wp)                  :: s_M, s_P, s_S
        type(riemann_states)      :: s
        real(wp)                  :: xi_M, xi_P
        type(riemann_states)      :: Ms, pres_S
        real(wp)                  :: alpha_L_sum, alpha_R_sum
        real(wp)                  :: zcoef, pcorr    !< low Mach number correction
        type(riemann_states)      :: c_fast, pres_mag
        type(riemann_states_vec3) :: B
        type(riemann_states)      :: Ga              !< Gamma (Lorentz factor)
        type(riemann_states)      :: vdotB, B2
        type(riemann_states_vec3) :: b4              !< 4-magnetic field components (spatial: b4x, b4y, b4z)
        type(riemann_states_vec3) :: cm              !< Conservative momentum variables
        integer                   :: i, j, k, l, q   !< Generic loop iterators
        integer, dimension(3)     :: idx_right_phys  !< Physical (j,k,l) indices for right state.
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
                $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, q, alpha_rho_L, alpha_rho_R, vel, alpha_L, alpha_R, tau_e, &
                                    & G, Re, rho_avg, h_avg, gamma_avg, s, s_S, Ys_L, Ys_R, xi_field, Cp_iL, Cp_iR, Xs_L, Xs_R, &
                                    & Gamma_iL, Gamma_iR, Yi_avg, Phi_avg, h_iL, h_iR, h_avg_2, c_fast, pres_mag, B, Ga, vdotB, &
                                    & B2, b4, cm, pcorr, zcoef, vel_grad_L, vel_grad_R, idx_right_phys, vel_rms, vel_avg_rms, &
                                    & vel_tmp, Ms, pres_S, alpha_L_sum, alpha_R_sum, c_avg, pres, rho, gamma, pi_inf, qv, c, E, &
                                    & H, ptilde, s_M, s_P, xi_M, xi_P, Cp_avg, Cv_avg, T_avg, eps, c_sum_Yi_Phi, Cp, Cv, Gamm, &
                                    & R_gas, MW, T, Y, flux_tau]')
                do l = ${Z_BND}$%beg, ${Z_BND}$%end
                    do k = ${Y_BND}$%beg, ${Y_BND}$%end
                        do j = ${X_BND}$%beg, ${X_BND}$%end
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, eqn_idx%cont%end
                                alpha_rho_L(i) = qL_prim_rsx_vf(${SF('')}$, i)
                                alpha_rho_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, i)
                            end do

                            vel_rms%L = 0._wp; vel_rms%R = 0._wp

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_vels
                                vel%L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%cont%end + i)
                                vel%R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%cont%end + i)
                                vel_rms%L = vel_rms%L + vel%L(i)**2._wp
                                vel_rms%R = vel_rms%R + vel%R(i)**2._wp
                            end do

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                alpha_L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%E + i)
                                alpha_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E + i)
                            end do

                            pres%L = qL_prim_rsx_vf(${SF('')}$, eqn_idx%E)
                            pres%R = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E)

                            if (mhd) then
                                if (n == 0) then  ! 1D: constant Bx; By, Bz as variables
                                    B%L(1) = Bx0
                                    B%R(1) = Bx0
                                    B%L(2) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%B%beg)
                                    B%R(2) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%B%beg)
                                    B%L(3) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%B%beg + 1)
                                    B%R(3) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%B%beg + 1)
                                else  ! 2D/3D: Bx, By, Bz as variables
                                    B%L(1) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%B%beg)
                                    B%R(1) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%B%beg)
                                    B%L(2) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%B%beg + 1)
                                    B%R(2) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%B%beg + 1)
                                    B%L(3) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%B%beg + 2)
                                    B%R(3) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%B%beg + 2)
                                end if
                            end if

                            rho%L = 0._wp
                            gamma%L = 0._wp
                            pi_inf%L = 0._wp
                            qv%L = 0._wp

                            rho%R = 0._wp
                            gamma%R = 0._wp
                            pi_inf%R = 0._wp
                            qv%R = 0._wp

                            alpha_L_sum = 0._wp
                            alpha_R_sum = 0._wp

                            pres_mag%L = 0._wp
                            pres_mag%R = 0._wp

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

                            if (viscous) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, 2
                                    Re%L(i) = dflt_real
                                    Re%R(i) = dflt_real

                                    if (Re_size(i) > 0) Re%L(i) = 0._wp
                                    if (Re_size(i) > 0) Re%R(i) = 0._wp

                                    $:GPU_LOOP(parallelism='[seq]')
                                    do q = 1, Re_size(i)
                                        Re%L(i) = alpha_L(Re_idx(i, q))/Res_gs(i, q) + Re%L(i)
                                        Re%R(i) = alpha_R(Re_idx(i, q))/Res_gs(i, q) + Re%R(i)
                                    end do

                                    Re%L(i) = 1._wp/max(Re%L(i), sgm_eps)
                                    Re%R(i) = 1._wp/max(Re%R(i), sgm_eps)
                                end do
                            end if

                            if (chemistry) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = eqn_idx%species%beg, eqn_idx%species%end
                                    Ys_L(i - eqn_idx%species%beg + 1) = qL_prim_rsx_vf(${SF('')}$, i)
                                    Ys_R(i - eqn_idx%species%beg + 1) = qR_prim_rsx_vf(${SF(' + 1')}$, i)
                                end do

                                call get_mixture_molecular_weight(Ys_L, MW%L)
                                call get_mixture_molecular_weight(Ys_R, MW%R)

                                Xs_L(:) = Ys_L(:)*MW%L/molecular_weights(:)
                                Xs_R(:) = Ys_R(:)*MW%R/molecular_weights(:)

                                R_gas%L = gas_constant/MW%L
                                R_gas%R = gas_constant/MW%R
                                T%L = pres%L/rho%L/R_gas%L
                                T%R = pres%R/rho%R/R_gas%R

                                call get_species_specific_heats_r(T%L, Cp_iL)
                                call get_species_specific_heats_r(T%R, Cp_iR)

                                if (chem_params%gamma_method == 1) then
                                    ! gamma_method = 1: Ref. Section 2.3.1 Formulation of doi:10.7907/ZKW8-ES97.
                                    Gamma_iL = Cp_iL/(Cp_iL - 1.0_wp)
                                    Gamma_iR = Cp_iR/(Cp_iR - 1.0_wp)

                                    gamma%L = sum(Xs_L(:)/(Gamma_iL(:) - 1.0_wp))
                                    gamma%R = sum(Xs_R(:)/(Gamma_iR(:) - 1.0_wp))
                                else if (chem_params%gamma_method == 2) then
                                    ! gamma_method = 2: c_p / c_v where c_p, c_v are specific heats.
                                    call get_mixture_specific_heat_cp_mass(T%L, Ys_L, Cp%L)
                                    call get_mixture_specific_heat_cp_mass(T%R, Ys_R, Cp%R)
                                    call get_mixture_specific_heat_cv_mass(T%L, Ys_L, Cv%L)
                                    call get_mixture_specific_heat_cv_mass(T%R, Ys_R, Cv%R)

                                    Gamm%L = Cp%L/Cv%L
                                    gamma%L = 1.0_wp/(Gamm%L - 1.0_wp)
                                    Gamm%R = Cp%R/Cv%R
                                    gamma%R = 1.0_wp/(Gamm%R - 1.0_wp)
                                end if

                                call get_mixture_energy_mass(T%L, Ys_L, E%L)
                                call get_mixture_energy_mass(T%R, Ys_R, E%R)

                                E%L = rho%L*E%L + 5.e-1*rho%L*vel_rms%L
                                E%R = rho%R*E%R + 5.e-1*rho%R*vel_rms%R
                                H%L = (E%L + pres%L)/rho%L
                                H%R = (E%R + pres%R)/rho%R
                            else if (mhd .and. relativity) then
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                                    Ga%L = 1._wp/sqrt(1._wp - vel_rms%L)
                                    Ga%R = 1._wp/sqrt(1._wp - vel_rms%R)
                                    vdotB%L = vel%L(1)*B%L(1) + vel%L(2)*B%L(2) + vel%L(3)*B%L(3)
                                    vdotB%R = vel%R(1)*B%R(1) + vel%R(2)*B%R(2) + vel%R(3)*B%R(3)

                                    b4%L(1:3) = B%L(1:3)/Ga%L + Ga%L*vel%L(1:3)*vdotB%L
                                    b4%R(1:3) = B%R(1:3)/Ga%R + Ga%R*vel%R(1:3)*vdotB%R
                                    B2%L = B%L(1)**2._wp + B%L(2)**2._wp + B%L(3)**2._wp
                                    B2%R = B%R(1)**2._wp + B%R(2)**2._wp + B%R(3)**2._wp

                                    pres_mag%L = 0.5_wp*(B2%L/Ga%L**2._wp + vdotB%L**2._wp)
                                    pres_mag%R = 0.5_wp*(B2%R/Ga%R**2._wp + vdotB%R**2._wp)

                                    ! Hard-coded EOS
                                    H%L = 1._wp + (gamma%L + 1)*pres%L/rho%L
                                    H%R = 1._wp + (gamma%R + 1)*pres%R/rho%R

                                    cm%L(1:3) = (rho%L*H%L*Ga%L**2 + B2%L)*vel%L(1:3) - vdotB%L*B%L(1:3)
                                    cm%R(1:3) = (rho%R*H%R*Ga%R**2 + B2%R)*vel%R(1:3) - vdotB%R*B%R(1:3)

                                    E%L = rho%L*H%L*Ga%L**2 - pres%L + 0.5_wp*(B2%L + vel_rms%L*B2%L - vdotB%L**2._wp) - rho%L*Ga%L
                                    E%R = rho%R*H%R*Ga%R**2 - pres%R + 0.5_wp*(B2%R + vel_rms%R*B2%R - vdotB%R**2._wp) - rho%R*Ga%R
                                #:endif
                            else if (mhd .and. .not. relativity) then
                                pres_mag%L = 0.5_wp*(B%L(1)**2._wp + B%L(2)**2._wp + B%L(3)**2._wp)
                                pres_mag%R = 0.5_wp*(B%R(1)**2._wp + B%R(2)**2._wp + B%R(3)**2._wp)
                                E%L = gamma%L*pres%L + pi_inf%L + 0.5_wp*rho%L*vel_rms%L + qv%L + pres_mag%L
                                ! includes magnetic energy
                                E%R = gamma%R*pres%R + pi_inf%R + 0.5_wp*rho%R*vel_rms%R + qv%R + pres_mag%R
                                H%L = (E%L + pres%L - pres_mag%L)/rho%L
                                ! stagnation enthalpy here excludes magnetic energy (only used to find speed of sound)
                                H%R = (E%R + pres%R - pres_mag%R)/rho%R
                            else
                                E%L = gamma%L*pres%L + pi_inf%L + 5.e-1*rho%L*vel_rms%L + qv%L
                                E%R = gamma%R*pres%R + pi_inf%R + 5.e-1*rho%R*vel_rms%R + qv%R
                                H%L = (E%L + pres%L)/rho%L
                                H%R = (E%R + pres%R)/rho%R
                            end if

                            ! elastic energy update
                            if (hypoelasticity) then
                                G%L = 0._wp; G%R = 0._wp

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    G%L = G%L + alpha_L(i)*Gs_rs(i)
                                    G%R = G%R + alpha_R(i)*Gs_rs(i)
                                end do

                                if (cont_damage) then
                                    G%L = G%L*max((1._wp - qL_prim_rsx_vf(${SF('')}$, eqn_idx%damage)), 0._wp)
                                    G%R = G%R*max((1._wp - qR_prim_rsx_vf(${SF('')}$, eqn_idx%damage)), 0._wp)
                                end if

                                do i = 1, eqn_idx%stress%end - eqn_idx%stress%beg + 1
                                    tau_e%L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%stress%beg - 1 + i)
                                    tau_e%R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%stress%beg - 1 + i)
                                    ! Elastic contribution to energy if G large enough TODO take out if statement if stable without
                                    if ((G%L > 1000) .and. (G%R > 1000)) then
                                        E%L = E%L + (tau_e%L(i)*tau_e%L(i))/(4._wp*G%L)
                                        E%R = E%R + (tau_e%R(i)*tau_e%R(i))/(4._wp*G%R)
                                        ! Double for shear stresses
                                        if (any(eqn_idx%stress%beg - 1 + i == shear_indices)) then
                                            E%L = E%L + (tau_e%L(i)*tau_e%L(i))/(4._wp*G%L)
                                            E%R = E%R + (tau_e%R(i)*tau_e%R(i))/(4._wp*G%R)
                                        end if
                                    end if
                                end do
                            end if

                            call s_compute_speed_of_sound(pres%L, rho%L, gamma%L, pi_inf%L, H%L, alpha_L, vel_rms%L, 0._wp, c%L, &
                                                          & qv%L)

                            call s_compute_speed_of_sound(pres%R, rho%R, gamma%R, pi_inf%R, H%R, alpha_R, vel_rms%R, 0._wp, c%R, &
                                                          & qv%R)

                            if (mhd) then
                                call s_compute_fast_magnetosonic_speed(rho%L, c%L, B%L, norm_dir, c_fast%L, H%L)
                                call s_compute_fast_magnetosonic_speed(rho%R, c%R, B%R, norm_dir, c_fast%R, H%R)
                            end if

                            s%L = 0._wp; s%R = 0._wp

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_dims
                                s%L = s%L + vel%L(i)**2._wp
                                s%R = s%R + vel%R(i)**2._wp
                            end do

                            s%L = sqrt(s%L)
                            s%R = sqrt(s%R)

                            s_P = max(s%L, s%R) + max(c%L, c%R)
                            s_M = -s_P

                            s%L = s_M
                            s%R = s_P

                            ! Low Mach correction
                            if (low_Mach == 1) then
                                @:compute_low_Mach_correction()
                            else
                                pcorr = 0._wp
                            end if

                            ! Mass
                            if (.not. relativity) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, eqn_idx%cont%end
                                    flux_rsx_vf(${SF('')}$, &
                                                & i) = (s_M*alpha_rho_R(i)*vel%R(norm_dir) - s_P*alpha_rho_L(i)*vel%L(norm_dir) &
                                                & + s_M*s_P*(alpha_rho_L(i) - alpha_rho_R(i)))/(s_M - s_P)
                                end do
                            else if (relativity) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, eqn_idx%cont%end
                                    flux_rsx_vf(${SF('')}$, &
                                                & i) = (s_M*Ga%R*alpha_rho_R(i)*vel%R(norm_dir) - s_P*Ga%L*alpha_rho_L(i) &
                                                & *vel%L(norm_dir) + s_M*s_P*(Ga%L*alpha_rho_L(i) - Ga%R*alpha_rho_R(i)))/(s_M &
                                                & - s_P)
                                end do
                            end if

                            ! Momentum
                            if (mhd .and. (.not. relativity)) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, 3
                                    ! Flux of rho*v_i in the ${XYZ}$ direction = rho * v_i * v_${XYZ}$ - B_i * B_${XYZ}$ +
                                    ! delta_(${XYZ}$,i) * p_tot
                                    flux_rsx_vf(${SF('')}$, &
                                                & eqn_idx%cont%end + i) = (s_M*(rho%R*vel%R(i)*vel%R(norm_dir) - B%R(i) &
                                                & *B%R(norm_dir) + dir_flg(i)*(pres%R + pres_mag%R)) - s_P*(rho%L*vel%L(i) &
                                                & *vel%L(norm_dir) - B%L(i)*B%L(norm_dir) + dir_flg(i)*(pres%L + pres_mag%L)) &
                                                & + s_M*s_P*(rho%L*vel%L(i) - rho%R*vel%R(i)))/(s_M - s_P)
                                end do
                            else if (mhd .and. relativity) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, 3
                                    ! Flux of m_i in the ${XYZ}$ direction = m_i * v_${XYZ}$ - b_i/Gamma * B_${XYZ}$ +
                                    ! delta_(${XYZ}$,i) * p_tot
                                    flux_rsx_vf(${SF('')}$, &
                                                & eqn_idx%cont%end + i) = (s_M*(cm%R(i)*vel%R(norm_dir) - b4%R(i) &
                                                & /Ga%R*B%R(norm_dir) + dir_flg(i)*(pres%R + pres_mag%R)) - s_P*(cm%L(i) &
                                                & *vel%L(norm_dir) - b4%L(i)/Ga%L*B%L(norm_dir) + dir_flg(i)*(pres%L + pres_mag%L) &
                                                & ) + s_M*s_P*(cm%L(i) - cm%R(i)))/(s_M - s_P)
                                end do
                            else if (bubbles_euler) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_vels
                                    flux_rsx_vf(${SF('')}$, &
                                                & eqn_idx%cont%end + dir_idx(i)) = (s_M*(rho%R*vel%R(dir_idx(1))*vel%R(dir_idx(i)) &
                                                & + dir_flg(dir_idx(i))*(pres%R - ptilde%R)) - s_P*(rho%L*vel%L(dir_idx(1)) &
                                                & *vel%L(dir_idx(i)) + dir_flg(dir_idx(i))*(pres%L - ptilde%L)) &
                                                & + s_M*s_P*(rho%L*vel%L(dir_idx(i)) - rho%R*vel%R(dir_idx(i))))/(s_M - s_P) &
                                                & + (s_M/s%L)*(s_P/s%R)*pcorr*(vel%R(dir_idx(i)) - vel%L(dir_idx(i)))
                                end do
                            else if (hypoelasticity) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_vels
                                    flux_rsx_vf(${SF('')}$, &
                                                & eqn_idx%cont%end + dir_idx(i)) = (s_M*(rho%R*vel%R(dir_idx(1))*vel%R(dir_idx(i)) &
                                                & + dir_flg(dir_idx(i))*pres%R - tau_e%R(dir_idx_tau(i))) &
                                                & - s_P*(rho%L*vel%L(dir_idx(1))*vel%L(dir_idx(i)) + dir_flg(dir_idx(i))*pres%L &
                                                & - tau_e%L(dir_idx_tau(i))) + s_M*s_P*(rho%L*vel%L(dir_idx(i)) &
                                                & - rho%R*vel%R(dir_idx(i))))/(s_M - s_P)
                                end do
                            else
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_vels
                                    flux_rsx_vf(${SF('')}$, &
                                                & eqn_idx%cont%end + dir_idx(i)) = (s_M*(rho%R*vel%R(dir_idx(1))*vel%R(dir_idx(i)) &
                                                & + dir_flg(dir_idx(i))*pres%R) - s_P*(rho%L*vel%L(dir_idx(1))*vel%L(dir_idx(i)) &
                                                & + dir_flg(dir_idx(i))*pres%L) + s_M*s_P*(rho%L*vel%L(dir_idx(i)) &
                                                & - rho%R*vel%R(dir_idx(i))))/(s_M - s_P) + (s_M/s%L)*(s_P/s%R) &
                                                & *pcorr*(vel%R(dir_idx(i)) - vel%L(dir_idx(i)))
                                end do
                            end if

                            ! Energy
                            if (mhd .and. (.not. relativity)) then
                                ! energy flux = (E + p + p_mag) * v_${XYZ}$ - B_${XYZ}$ * (v_x*B_x + v_y*B_y + v_z*B_z)
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
                                    flux_rsx_vf(${SF('')}$, &
                                                & eqn_idx%E) = (s_M*(vel%R(norm_dir)*(E%R + pres%R + pres_mag%R) - B%R(norm_dir) &
                                                & *(vel%R(1)*B%R(1) + vel%R(2)*B%R(2) + vel%R(3)*B%R(3))) - s_P*(vel%L(norm_dir) &
                                                & *(E%L + pres%L + pres_mag%L) - B%L(norm_dir)*(vel%L(1)*B%L(1) + vel%L(2)*B%L(2) &
                                                & + vel%L(3)*B%L(3))) + s_M*s_P*(E%L - E%R))/(s_M - s_P)
                                #:endif
                            else if (mhd .and. relativity) then
                                ! energy flux = m_${XYZ}$ - mass flux Hard-coded for single-component for now
                                flux_rsx_vf(${SF('')}$, &
                                            & eqn_idx%E) = (s_M*(cm%R(norm_dir) - Ga%R*alpha_rho_R(1)*vel%R(norm_dir)) &
                                            & - s_P*(cm%L(norm_dir) - Ga%L*alpha_rho_L(1)*vel%L(norm_dir)) + s_M*s_P*(E%L - E%R)) &
                                            & /(s_M - s_P)
                            else if (bubbles_euler) then
                                flux_rsx_vf(${SF('')}$, &
                                            & eqn_idx%E) = (s_M*vel%R(dir_idx(1))*(E%R + pres%R - ptilde%R) - s_P*vel%L(dir_idx(1) &
                                            & )*(E%L + pres%L - ptilde%L) + s_M*s_P*(E%L - E%R))/(s_M - s_P) + (s_M/s%L)*(s_P/s%R) &
                                            & *pcorr*(vel_rms%R - vel_rms%L)/2._wp
                            else if (hypoelasticity) then
                                flux_tau%L = 0._wp; flux_tau%R = 0._wp
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    flux_tau%L = flux_tau%L + tau_e%L(dir_idx_tau(i))*vel%L(dir_idx(i))
                                    flux_tau%R = flux_tau%R + tau_e%R(dir_idx_tau(i))*vel%R(dir_idx(i))
                                end do
                                flux_rsx_vf(${SF('')}$, &
                                            & eqn_idx%E) = (s_M*(vel%R(dir_idx(1))*(E%R + pres%R) - flux_tau%R) &
                                            & - s_P*(vel%L(dir_idx(1))*(E%L + pres%L) - flux_tau%L) + s_M*s_P*(E%L - E%R))/(s_M &
                                            & - s_P)
                            else
                                flux_rsx_vf(${SF('')}$, &
                                            & eqn_idx%E) = (s_M*vel%R(dir_idx(1))*(E%R + pres%R) - s_P*vel%L(dir_idx(1))*(E%L &
                                            & + pres%L) + s_M*s_P*(E%L - E%R))/(s_M - s_P) + (s_M/s%L)*(s_P/s%R)*pcorr*(vel_rms%R &
                                            & - vel_rms%L)/2._wp
                            end if

                            ! Elastic Stresses
                            if (hypoelasticity) then
                                do i = 1, eqn_idx%stress%end - eqn_idx%stress%beg + 1  ! TODO: this indexing may be slow
                                    flux_rsx_vf(${SF('')}$, &
                                                & eqn_idx%stress%beg - 1 + i) = (s_M*(rho%R*vel%R(dir_idx(1))*tau_e%R(i)) &
                                                & - s_P*(rho%L*vel%L(dir_idx(1))*tau_e%L(i)) + s_M*s_P*(rho%L*tau_e%L(i) &
                                                & - rho%R*tau_e%R(i)))/(s_M - s_P)
                                end do
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
                                    Y%L = qL_prim_rsx_vf(${SF('')}$, i)
                                    Y%R = qR_prim_rsx_vf(${SF(' + 1')}$, i)

                                    flux_rsx_vf(${SF('')}$, &
                                                & i) = (s_M*Y%R*rho%R*vel%R(dir_idx(1)) - s_P*Y%L*rho%L*vel%L(dir_idx(1)) &
                                                & + s_M*s_P*(Y%L*rho%L - Y%R*rho%R))/(s_M - s_P)
                                    flux_src_rsx_vf(${SF('')}$, i) = 0._wp
                                end do
                            end if

                            ! MHD: magnetic flux and Maxwell stress contributions
                            if (mhd) then
                                if (n == 0) then  ! 1D: d/dx flux only & Bx = Bx0 = const.
                                    ! B_y flux = v_x * B_y - v_y * Bx0 B_z flux = v_x * B_z - v_z * Bx0
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 0, 1
                                        flux_rsx_vf(j, k, l, &
                                                    & eqn_idx%B%beg + i) = (s_M*(vel%R(1)*B%R(2 + i) - vel%R(2 + i)*Bx0) &
                                                    & - s_P*(vel%L(1)*B%L(2 + i) - vel%L(2 + i)*Bx0) + s_M*s_P*(B%L(2 + i) &
                                                    & - B%R(2 + i)))/(s_M - s_P)
                                    end do
                                else  ! 2D/3D: Bx, By, Bz /= const. but zero flux component in the same direction
                                    ! B_x d/d${XYZ}$ flux = (1 - delta(x,${XYZ}$)) * (v_${XYZ}$ * B_x - v_x * B_${XYZ}$) B_y
                                    ! d/d${XYZ}$ flux = (1 - delta(y,${XYZ}$)) * (v_${XYZ}$ * B_y - v_y * B_${XYZ}$) B_z d/d${XYZ}$
                                    ! flux = (1 - delta(z,${XYZ}$)) * (v_${XYZ}$ * B_z - v_z * B_${XYZ}$)
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 0, 2
                                        flux_rsx_vf(${SF('')}$, &
                                                    & eqn_idx%B%beg + i) = (1 - dir_flg(i + 1))*(s_M*(vel%R(dir_idx(1))*B%R(i + 1) &
                                                    & - vel%R(i + 1)*B%R(norm_dir)) - s_P*(vel%L(dir_idx(1))*B%L(i + 1) - vel%L(i &
                                                    & + 1)*B%L(norm_dir)) + s_M*s_P*(B%L(i + 1) - B%R(i + 1)))/(s_M - s_P)
                                    end do
                                end if
                                flux_src_rsx_vf(${SF('')}$, eqn_idx%adv%beg) = 0._wp
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
                                                     & eqn_idx%cont%end + 2) - (s_M*pres%R - s_P*pres%L)/(s_M - s_P)
                                    ! Geometrical source of the void fraction(s) is zero
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = eqn_idx%adv%beg, eqn_idx%adv%end
                                        flux_gsrc_rsx_vf(${SF('')}$, i) = flux_rsx_vf(${SF('')}$, i)
                                    end do
                                end if

                                if (cyl_coord .and. hypoelasticity) then
                                    ! += tau_sigmasigma using HLL
                                    flux_gsrc_rsx_vf(${SF('')}$, eqn_idx%cont%end + 2) = flux_gsrc_rsx_vf(${SF('')}$, &
                                                     & eqn_idx%cont%end + 2) + (s_M*tau_e%R(4) - s_P*tau_e%L(4))/(s_M - s_P)

                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = eqn_idx%stress%beg, eqn_idx%stress%end
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
            $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, idx_right_phys, vel_grad_L, vel_grad_R, alpha_L, alpha_R, vel, &
                                & Re]', copyin='[norm_dir]')
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
                                vel%L(i) = qL_prim_rsx_vf(j, k, l, eqn_idx%mom%beg + i - 1)
                                vel%R(i) = qR_prim_rsx_vf(j + 1, k, l, eqn_idx%mom%beg + i - 1)
                            end do
                        else if (norm_dir == 2) then
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                alpha_L(i) = qL_prim_rsx_vf(j, k, l, eqn_idx%E + i)
                                alpha_R(i) = qR_prim_rsx_vf(j, k + 1, l, eqn_idx%E + i)
                            end do
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_dims
                                vel%L(i) = qL_prim_rsx_vf(j, k, l, eqn_idx%mom%beg + i - 1)
                                vel%R(i) = qR_prim_rsx_vf(j, k + 1, l, eqn_idx%mom%beg + i - 1)
                            end do
                        else
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                alpha_L(i) = qL_prim_rsx_vf(j, k, l, eqn_idx%E + i)
                                alpha_R(i) = qR_prim_rsx_vf(j, k, l + 1, eqn_idx%E + i)
                            end do

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_dims
                                vel%L(i) = qL_prim_rsx_vf(j, k, l, eqn_idx%mom%beg + i - 1)
                                vel%R(i) = qR_prim_rsx_vf(j, k, l + 1, eqn_idx%mom%beg + i - 1)
                            end do
                        end if

                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, 2
                            Re%L(i) = dflt_real
                            Re%R(i) = dflt_real

                            if (Re_size(i) > 0) Re%L(i) = 0._wp
                            if (Re_size(i) > 0) Re%R(i) = 0._wp

                            $:GPU_LOOP(parallelism='[seq]')
                            do q = 1, Re_size(i)
                                Re%L(i) = alpha_L(Re_idx(i, q))/Res_gs(i, q) + Re%L(i)
                                Re%R(i) = alpha_R(Re_idx(i, q))/Res_gs(i, q) + Re%R(i)
                            end do

                            Re%L(i) = 1._wp/max(Re%L(i), sgm_eps)
                            Re%R(i) = 1._wp/max(Re%R(i), sgm_eps)
                        end do

                        if (shear_stress) then
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_dims
                                vel_grad_L(i, 1) = (dqL_prim_dx_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l)/Re%L(1))
                                vel_grad_R(i, 1) = (dqR_prim_dx_vf(eqn_idx%mom%beg + i - 1)%sf(idx_right_phys(1), &
                                           & idx_right_phys(2), idx_right_phys(3))/Re%R(1))
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
                                    if (num_dims > 1) then
                                        vel_grad_L(i, 2) = (dqL_prim_dy_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l)/Re%L(1))
                                        vel_grad_R(i, 2) = (dqR_prim_dy_vf(eqn_idx%mom%beg + i - 1)%sf(idx_right_phys(1), &
                                                   & idx_right_phys(2), idx_right_phys(3))/Re%R(1))
                                    end if
                                    #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                                        if (num_dims > 2) then
                                            vel_grad_L(i, 3) = (dqL_prim_dz_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l)/Re%L(1))
                                            vel_grad_R(i, 3) = (dqR_prim_dz_vf(eqn_idx%mom%beg + i - 1)%sf(idx_right_phys(1), &
                                                       & idx_right_phys(2), idx_right_phys(3))/Re%R(1))
                                        end if
                                    #:endif
                                #:endif
                            end do

                            if (norm_dir == 1) then
                                flux_src_vf(eqn_idx%mom%beg)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg)%sf(j, k, &
                                            & l) - (4._wp/3._wp)*0.5_wp*(vel_grad_L(1, 1) + vel_grad_R(1, 1))
                                flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                            & l) - (4._wp/3._wp)*0.5_wp*(vel_grad_L(1, 1)*vel%L(1) + vel_grad_R(1, 1)*vel%R(1))
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
                                    if (num_dims > 1) then
                                        flux_src_vf(eqn_idx%mom%beg)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg)%sf(j, k, &
                                                    & l) - (-2._wp/3._wp)*0.5_wp*(vel_grad_L(2, 2) + vel_grad_R(2, 2))
                                        flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                    & l) - (-2._wp/3._wp)*0.5_wp*(vel_grad_L(2, 2)*vel%L(1) + vel_grad_R(2, &
                                                    & 2)*vel%R(1))

                                        flux_src_vf(eqn_idx%mom%beg + 1)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + 1)%sf(j, k, &
                                                    & l) - 0.5_wp*(vel_grad_L(1, 2) + vel_grad_R(1, 2)) - 0.5_wp*(vel_grad_L(2, &
                                                    & 1) + vel_grad_R(2, 1))
                                        flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                    & l) - 0.5_wp*(vel_grad_L(1, 2)*vel%L(2) + vel_grad_R(1, &
                                                    & 2)*vel%R(2)) - 0.5_wp*(vel_grad_L(2, 1)*vel%L(2) + vel_grad_R(2, 1)*vel%R(2))
                                        #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                                            if (num_dims > 2) then
                                                flux_src_vf(eqn_idx%mom%beg)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg)%sf(j, k, &
                                                            & l) - (-2._wp/3._wp)*0.5_wp*(vel_grad_L(3, 3) + vel_grad_R(3, 3))
                                                flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                            & l) - (-2._wp/3._wp)*0.5_wp*(vel_grad_L(3, &
                                                            & 3)*vel%L(1) + vel_grad_R(3, 3)*vel%R(1))

                                                flux_src_vf(eqn_idx%mom%beg + 2)%sf(j, k, &
                                                            & l) = flux_src_vf(eqn_idx%mom%beg + 2)%sf(j, k, &
                                                            & l) - 0.5_wp*(vel_grad_L(1, 3) + vel_grad_R(1, &
                                                            & 3)) - 0.5_wp*(vel_grad_L(3, 1) + vel_grad_R(3, 1))
                                                flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                            & l) - 0.5_wp*(vel_grad_L(1, 3)*vel%L(3) + vel_grad_R(1, &
                                                            & 3)*vel%R(3)) - 0.5_wp*(vel_grad_L(3, 1)*vel%L(3) + vel_grad_R(3, &
                                                            & 1)*vel%R(3))
                                            end if
                                        #:endif
                                    end if
                                #:endif
                            else if (norm_dir == 2) then
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
                                    flux_src_vf(eqn_idx%mom%beg + 1)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + 1)%sf(j, k, &
                                                & l) - (-2._wp/3._wp)*0.5_wp*(vel_grad_L(1, 1) + vel_grad_R(1, 1))
                                    flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                & l) - (-2._wp/3._wp)*0.5_wp*(vel_grad_L(1, 1)*vel%L(2) + vel_grad_R(1, 1)*vel%R(2))

                                    flux_src_vf(eqn_idx%mom%beg + 1)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + 1)%sf(j, k, &
                                                & l) - (4._wp/3._wp)*0.5_wp*(vel_grad_L(2, 2) + vel_grad_R(2, 2))
                                    flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                & l) - (4._wp/3._wp)*0.5_wp*(vel_grad_L(2, 2)*vel%L(2) + vel_grad_R(2, 2)*vel%R(2))

                                    flux_src_vf(eqn_idx%mom%beg)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg)%sf(j, k, &
                                                & l) - 0.5_wp*(vel_grad_L(1, 2) + vel_grad_R(1, 2)) - 0.5_wp*(vel_grad_L(2, &
                                                & 1) + vel_grad_R(2, 1))
                                    flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                & l) - 0.5_wp*(vel_grad_L(1, 2)*vel%L(1) + vel_grad_R(1, &
                                                & 2)*vel%R(1)) - 0.5_wp*(vel_grad_L(2, 1)*vel%L(1) + vel_grad_R(2, 1)*vel%R(1))
                                    #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                                        if (num_dims > 2) then
                                            flux_src_vf(eqn_idx%mom%beg + 1)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + 1)%sf(j, &
                                                        & k, l) - (-2._wp/3._wp)*0.5_wp*(vel_grad_L(3, 3) + vel_grad_R(3, 3))
                                            flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                        & l) - (-2._wp/3._wp)*0.5_wp*(vel_grad_L(3, 3)*vel%L(2) + vel_grad_R(3, &
                                                        & 3)*vel%R(2))

                                            flux_src_vf(eqn_idx%mom%beg + 2)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + 2)%sf(j, &
                                                        & k, l) - 0.5_wp*(vel_grad_L(2, 3) + vel_grad_R(2, &
                                                        & 3)) - 0.5_wp*(vel_grad_L(3, 2) + vel_grad_R(3, 2))
                                            flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                        & l) - 0.5_wp*(vel_grad_L(2, 3)*vel%L(3) + vel_grad_R(2, &
                                                        & 3)*vel%R(3)) - 0.5_wp*(vel_grad_L(3, 2)*vel%L(3) + vel_grad_R(3, &
                                                        & 2)*vel%R(3))
                                        end if
                                    #:endif
                                #:endif
                            else
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                                    flux_src_vf(eqn_idx%mom%beg + 2)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + 2)%sf(j, k, &
                                                & l) - (-2._wp/3._wp)*0.5_wp*(vel_grad_L(1, 1) + vel_grad_R(1, 1))
                                    flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                & l) - (-2._wp/3._wp)*0.5_wp*(vel_grad_L(1, 1)*vel%L(3) + vel_grad_R(1, 1)*vel%R(3))

                                    flux_src_vf(eqn_idx%mom%beg + 2)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + 2)%sf(j, k, &
                                                & l) - (-2._wp/3._wp)*0.5_wp*(vel_grad_L(2, 2) + vel_grad_R(2, 2))
                                    flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                & l) - (-2._wp/3._wp)*0.5_wp*(vel_grad_L(2, 2)*vel%L(3) + vel_grad_R(2, 2)*vel%R(3))

                                    flux_src_vf(eqn_idx%mom%beg)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg)%sf(j, k, &
                                                & l) - 0.5_wp*(vel_grad_L(1, 3) + vel_grad_R(1, 3)) - 0.5_wp*(vel_grad_L(3, &
                                                & 1) + vel_grad_R(3, 1))
                                    flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                & l) - 0.5_wp*(vel_grad_L(1, 3)*vel%L(1) + vel_grad_R(1, &
                                                & 3)*vel%R(1)) - 0.5_wp*(vel_grad_L(3, 1)*vel%L(1) + vel_grad_R(3, 1)*vel%R(1))

                                    flux_src_vf(eqn_idx%mom%beg + 2)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + 2)%sf(j, k, &
                                                & l) - (4._wp/3._wp)*0.5_wp*(vel_grad_L(3, 3) + vel_grad_R(3, 3))
                                    flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                & l) - (4._wp/3._wp)*0.5_wp*(vel_grad_L(3, 3)*vel%L(3) + vel_grad_R(3, 3)*vel%R(3))

                                    flux_src_vf(eqn_idx%mom%beg + 1)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + 1)%sf(j, k, &
                                                & l) - 0.5_wp*(vel_grad_L(2, 3) + vel_grad_R(2, 3)) - 0.5_wp*(vel_grad_L(3, &
                                                & 2) + vel_grad_R(3, 2))
                                    flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                & l) - 0.5_wp*(vel_grad_L(2, 3)*vel%L(2) + vel_grad_R(2, &
                                                & 3)*vel%R(2)) - 0.5_wp*(vel_grad_L(3, 2)*vel%L(2) + vel_grad_R(3, 2)*vel%R(2))
                                #:endif
                            end if
                        end if

                        if (bulk_stress) then
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_dims
                                vel_grad_L(i, 1) = (dqL_prim_dx_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l)/Re%L(2))
                                vel_grad_R(i, 1) = (dqR_prim_dx_vf(eqn_idx%mom%beg + i - 1)%sf(idx_right_phys(1), &
                                           & idx_right_phys(2), idx_right_phys(3))/Re%R(2))
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
                                    if (num_dims > 1) then
                                        vel_grad_L(i, 2) = (dqL_prim_dy_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l)/Re%L(2))
                                        vel_grad_R(i, 2) = (dqR_prim_dy_vf(eqn_idx%mom%beg + i - 1)%sf(idx_right_phys(1), &
                                                   & idx_right_phys(2), idx_right_phys(3))/Re%R(2))
                                    end if
                                #:endif
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                                    if (num_dims > 2) then
                                        vel_grad_L(i, 3) = (dqL_prim_dz_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l)/Re%L(2))
                                        vel_grad_R(i, 3) = (dqR_prim_dz_vf(eqn_idx%mom%beg + i - 1)%sf(idx_right_phys(1), &
                                                   & idx_right_phys(2), idx_right_phys(3))/Re%R(2))
                                    end if
                                #:endif
                            end do

                            if (norm_dir == 1) then
                                flux_src_vf(eqn_idx%mom%beg)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg)%sf(j, k, &
                                            & l) - 0.5_wp*(vel_grad_L(1, 1) + vel_grad_R(1, 1))
                                flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, l) - 0.5_wp*(vel_grad_L(1, &
                                            & 1)*vel%L(1) + vel_grad_R(1, 1)*vel%R(1))
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
                                    if (num_dims > 1) then
                                        flux_src_vf(eqn_idx%mom%beg)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg)%sf(j, k, &
                                                    & l) - 0.5_wp*(vel_grad_L(2, 2) + vel_grad_R(2, 2))
                                        flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                    & l) - 0.5_wp*(vel_grad_L(2, 2)*vel%L(1) + vel_grad_R(2, 2)*vel%R(1))

                                        #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                                            if (num_dims > 2) then
                                                flux_src_vf(eqn_idx%mom%beg)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg)%sf(j, k, &
                                                            & l) - 0.5_wp*(vel_grad_L(3, 3) + vel_grad_R(3, 3))
                                                flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                            & l) - 0.5_wp*(vel_grad_L(3, 3)*vel%L(1) + vel_grad_R(3, 3)*vel%R(1))
                                            end if
                                        #:endif
                                    end if
                                #:endif
                            else if (norm_dir == 2) then
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
                                    flux_src_vf(eqn_idx%mom%beg + 1)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + 1)%sf(j, k, &
                                                & l) - 0.5_wp*(vel_grad_L(1, 1) + vel_grad_R(1, 1))
                                    flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                & l) - 0.5_wp*(vel_grad_L(1, 1)*vel%L(2) + vel_grad_R(1, 1)*vel%R(2))

                                    flux_src_vf(eqn_idx%mom%beg + 1)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + 1)%sf(j, k, &
                                                & l) - 0.5_wp*(vel_grad_L(2, 2) + vel_grad_R(2, 2))
                                    flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                & l) - 0.5_wp*(vel_grad_L(2, 2)*vel%L(2) + vel_grad_R(2, 2)*vel%R(2))

                                    #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                                        if (num_dims > 2) then
                                            flux_src_vf(eqn_idx%mom%beg + 1)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + 1)%sf(j, &
                                                        & k, l) - 0.5_wp*(vel_grad_L(3, 3) + vel_grad_R(3, 3))
                                            flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                        & l) - 0.5_wp*(vel_grad_L(3, 3)*vel%L(2) + vel_grad_R(3, 3)*vel%R(2))
                                        end if
                                    #:endif
                                #:endif
                            else
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                                    flux_src_vf(eqn_idx%mom%beg + 2)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + 2)%sf(j, k, &
                                                & l) - 0.5_wp*(vel_grad_L(1, 1) + vel_grad_R(1, 1))
                                    flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                & l) - 0.5_wp*(vel_grad_L(1, 1)*vel%L(3) + vel_grad_R(1, 1)*vel%R(3))

                                    flux_src_vf(eqn_idx%mom%beg + 2)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + 2)%sf(j, k, &
                                                & l) - 0.5_wp*(vel_grad_L(2, 2) + vel_grad_R(2, 2))
                                    flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                & l) - 0.5_wp*(vel_grad_L(2, 2)*vel%L(3) + vel_grad_R(2, 2)*vel%R(3))

                                    flux_src_vf(eqn_idx%mom%beg + 2)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + 2)%sf(j, k, &
                                                & l) - 0.5_wp*(vel_grad_L(3, 3) + vel_grad_R(3, 3))
                                    flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                                & l) - 0.5_wp*(vel_grad_L(3, 3)*vel%L(3) + vel_grad_R(3, 3)*vel%R(3))
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

    !> HLLC Riemann solver with contact restoration, Toro et al. Shock Waves (1994)
    subroutine s_hllc_riemann_solver(qL_prim_rsx_vf, dqL_prim_dx_vf, dqL_prim_dy_vf, dqL_prim_dz_vf, qL_prim_vf, qR_prim_rsx_vf, &
                                     & dqR_prim_dx_vf, dqR_prim_dy_vf, dqR_prim_dz_vf, qR_prim_vf, q_prim_vf, flux_vf, &
                                     & flux_src_vf, flux_gsrc_vf, norm_dir, ix, iy, iz)

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
        #:else
            real(wp), dimension(num_fluids) :: alpha_rho_L, alpha_rho_R
            real(wp), dimension(num_fluids) :: alpha_L, alpha_R
        #:endif
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(10) :: Ys_L, Ys_R, Xs_L, Xs_R, Gamma_iL, Gamma_iR, Cp_iL, Cp_iR
            real(wp), dimension(10) :: Yi_avg, Phi_avg, h_iL, h_iR, h_avg_2
        #:else
            real(wp), dimension(num_species) :: Ys_L, Ys_R, Xs_L, Xs_R, Gamma_iL, Gamma_iR, Cp_iL, Cp_iR
            real(wp), dimension(num_species) :: Yi_avg, Phi_avg, h_iL, h_iR, h_avg_2
        #:endif
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
        type(riemann_states_vec3) :: vel
        type(riemann_states)      :: rho, pres, E, H
        type(riemann_states)      :: T, Y, MW, R_gas, Cp, Cv, Gamm, gamma
        type(riemann_states)      :: pi_inf, qv, c, G, ptilde, nbub
        type(riemann_states)      :: vel_rms, vel_tmp
        real(wp)                  :: Cp_avg, Cv_avg, T_avg, c_sum_Yi_Phi, eps
        type(riemann_states_arr2) :: Re
        type(riemann_states_arr6) :: tau_e
        type(riemann_states_vec3) :: xi_field
        real(wp)                  :: rho_avg
        real(wp)                  :: H_avg
        real(wp)                  :: gamma_avg
        real(wp)                  :: qv_avg
        real(wp)                  :: c_avg
        real(wp)                  :: vel_avg_rms
        real(wp)                  :: s_M, s_P, s_S
        type(riemann_states)      :: s
        type(riemann_states)      :: xi, xi_m1              !< Left and right wave speeds functions / xi - 1
        real(wp)                  :: xi_M, xi_P
        real(wp)                  :: xi_MP, xi_PP
        real(wp)                  :: alpha_L_sum, alpha_R_sum
        real(wp)                  :: PbwR3Lbar, PbwR3Rbar
        real(wp)                  :: R3Lbar, R3Rbar
        real(wp)                  :: R3V2Lbar, R3V2Rbar
        real(wp)                  :: rho_Star, E_Star, p_Star, p_K_Star, vel_K_star
        type(riemann_states)      :: pres_S, Ms
        real(wp)                  :: flux_ene_e
        real(wp)                  :: zcoef, pcorr           !< low Mach number correction
        integer                   :: Re_max, i, j, k, l, q  !< Generic loop iterators
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
                if (model_eqns == 3) then
                    ! 6-equation model (model_eqns=3): separate phasic internal energies
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, q, vel, Re, alpha_L, alpha_R, Ys_L, Ys_R, Xs_L, Xs_R, &
                                        & Gamma_iL, Gamma_iR, Cp_iL, Cp_iR, Yi_avg, Phi_avg, h_iL, h_iR, h_avg_2, tau_e, &
                                        & flux_ene_e, xi_field, pcorr, zcoef, rho, pres, E, H, Cp_avg, Cv_avg, T_avg, eps, &
                                        & c_sum_Yi_Phi, T, Y, MW, R_gas, Cp, Cv, Gamm, gamma, pi_inf, qv, qv_avg, c, G, rho_avg, &
                                        & H_avg, c_avg, gamma_avg, ptilde, vel_rms, vel_avg_rms, vel_tmp, Ms, pres_S, &
                                        & alpha_L_sum, alpha_R_sum, rho_Star, E_Star, p_Star, p_K_Star, vel_K_star, s, s_M, s_P, &
                                        & s_S, xi_M, xi_P, xi, xi_m1, xi_MP, xi_PP]')
                    do l = ${Z_BND}$%beg, ${Z_BND}$%end
                        do k = ${Y_BND}$%beg, ${Y_BND}$%end
                            do j = ${X_BND}$%beg, ${X_BND}$%end
                                vel_rms%L = 0._wp; vel_rms%R = 0._wp
                                rho%L = 0._wp; rho%R = 0._wp
                                gamma%L = 0._wp; gamma%R = 0._wp
                                pi_inf%L = 0._wp; pi_inf%R = 0._wp
                                qv%L = 0._wp; qv%R = 0._wp
                                alpha_L_sum = 0._wp; alpha_R_sum = 0._wp

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel%L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%cont%end + i)
                                    vel%R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%cont%end + i)
                                    vel_rms%L = vel_rms%L + vel%L(i)**2._wp
                                    vel_rms%R = vel_rms%R + vel%R(i)**2._wp
                                end do

                                pres%L = qL_prim_rsx_vf(${SF('')}$, eqn_idx%E)
                                pres%R = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E)

                                rho%L = 0._wp
                                gamma%L = 0._wp
                                pi_inf%L = 0._wp
                                qv%L = 0._wp

                                rho%R = 0._wp
                                gamma%R = 0._wp
                                pi_inf%R = 0._wp
                                qv%R = 0._wp

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
                                    rho%L = rho%L + qL_prim_rsx_vf(${SF('')}$, i)
                                    gamma%L = gamma%L + qL_prim_rsx_vf(${SF('')}$, eqn_idx%E + i)*gammas(i)
                                    pi_inf%L = pi_inf%L + qL_prim_rsx_vf(${SF('')}$, eqn_idx%E + i)*pi_infs(i)
                                    qv%L = qv%L + qL_prim_rsx_vf(${SF('')}$, i)*qvs(i)

                                    rho%R = rho%R + qR_prim_rsx_vf(${SF(' + 1')}$, i)
                                    gamma%R = gamma%R + qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E + i)*gammas(i)
                                    pi_inf%R = pi_inf%R + qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E + i)*pi_infs(i)
                                    qv%R = qv%R + qR_prim_rsx_vf(${SF(' + 1')}$, i)*qvs(i)

                                    alpha_L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%adv%beg + i - 1)
                                    alpha_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%adv%beg + i - 1)
                                end do

                                if (viscous) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, 2
                                        Re%L(i) = dflt_real
                                        Re%R(i) = dflt_real
                                        if (Re_size(i) > 0) Re%L(i) = 0._wp
                                        if (Re_size(i) > 0) Re%R(i) = 0._wp
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do q = 1, Re_size(i)
                                            Re%L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%E + Re_idx(i, q))/Res_gs(i, q) + Re%L(i)
                                            Re%R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E + Re_idx(i, q))/Res_gs(i, &
                                                 & q) + Re%R(i)
                                        end do
                                        Re%L(i) = 1._wp/max(Re%L(i), sgm_eps)
                                        Re%R(i) = 1._wp/max(Re%R(i), sgm_eps)
                                    end do
                                end if

                                E%L = gamma%L*pres%L + pi_inf%L + 5.e-1_wp*rho%L*vel_rms%L + qv%L
                                E%R = gamma%R*pres%R + pi_inf%R + 5.e-1_wp*rho%R*vel_rms%R + qv%R

                                ! ENERGY ADJUSTMENTS FOR HYPOELASTIC ENERGY
                                if (hypoelasticity) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, eqn_idx%stress%end - eqn_idx%stress%beg + 1
                                        tau_e%L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%stress%beg - 1 + i)
                                        tau_e%R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%stress%beg - 1 + i)
                                    end do
                                    G%L = 0._wp; G%R = 0._wp
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        G%L = G%L + alpha_L(i)*Gs_rs(i)
                                        G%R = G%R + alpha_R(i)*Gs_rs(i)
                                    end do
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, eqn_idx%stress%end - eqn_idx%stress%beg + 1
                                        ! Elastic contribution to energy if G large enough
                                        if ((G%L > verysmall) .and. (G%R > verysmall)) then
                                            E%L = E%L + (tau_e%L(i)*tau_e%L(i))/(4._wp*G%L)
                                            E%R = E%R + (tau_e%R(i)*tau_e%R(i))/(4._wp*G%R)
                                            ! Additional terms in 2D and 3D
                                            if ((i == 2) .or. (i == 4) .or. (i == 5)) then
                                                E%L = E%L + (tau_e%L(i)*tau_e%L(i))/(4._wp*G%L)
                                                E%R = E%R + (tau_e%R(i)*tau_e%R(i))/(4._wp*G%R)
                                            end if
                                        end if
                                    end do
                                end if

                                ! Hyperelastic stress contribution: strain energy added to total energy
                                if (hyperelasticity) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        xi_field%L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%xi%beg - 1 + i)
                                        xi_field%R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%xi%beg - 1 + i)
                                    end do
                                    G%L = 0._wp; G%R = 0._wp
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        ! Mixture left and right shear modulus
                                        G%L = G%L + alpha_L(i)*Gs_rs(i)
                                        G%R = G%R + alpha_R(i)*Gs_rs(i)
                                    end do
                                    ! Elastic contribution to energy if G large enough
                                    if (G%L > verysmall .and. G%R > verysmall) then
                                        E%L = E%L + G%L*qL_prim_rsx_vf(${SF('')}$, eqn_idx%xi%end + 1)
                                        E%R = E%R + G%R*qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%xi%end + 1)
                                    end if
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, b_size - 1
                                        tau_e%L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%stress%beg - 1 + i)
                                        tau_e%R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%stress%beg - 1 + i)
                                    end do
                                end if

                                H%L = (E%L + pres%L)/rho%L
                                H%R = (E%R + pres%R)/rho%R

                                @:compute_average_state()

                                call s_compute_speed_of_sound(pres%L, rho%L, gamma%L, pi_inf%L, H%L, alpha_L, vel_rms%L, 0._wp, &
                                                              & c%L, qv%L)

                                call s_compute_speed_of_sound(pres%R, rho%R, gamma%R, pi_inf%R, H%R, alpha_R, vel_rms%R, 0._wp, &
                                                              & c%R, qv%R)

                                !> The computation of c_avg does not require all the variables, and therefore the non '_avg'
                                ! variables are placeholders to call the subroutine.
                                call s_compute_speed_of_sound(pres%R, rho_avg, gamma_avg, pi_inf%R, H_avg, alpha_R, vel_avg_rms, &
                                                              & 0._wp, c_avg, qv_avg)

                                if (viscous) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, 2
                                        Re_avg_rsx_vf(${SF('')}$, i) = 2._wp/(1._wp/Re%L(i) + 1._wp/Re%R(i))
                                    end do
                                end if

                                ! Low Mach correction
                                if (low_Mach == 2) then
                                    @:compute_low_Mach_correction()
                                end if

                                ! COMPUTING THE DIRECT WAVE SPEEDS
                                if (wave_speeds == 1) then
                                    if (elasticity) then
                                        ! Elastic wave speed, Rodriguez et al. JCP (2019)
                                        s%L = min(vel%L(dir_idx(1)) - sqrt(c%L*c%L + (((4._wp*G%L)/3._wp) + tau_e%L(dir_idx_tau(1) &
                                                  & ))/rho%L), &
                                                  & vel%R(dir_idx(1)) - sqrt(c%R*c%R + (((4._wp*G%R)/3._wp) &
                                                  & + tau_e%R(dir_idx_tau(1)))/rho%R))
                                        s%R = max(vel%R(dir_idx(1)) + sqrt(c%R*c%R + (((4._wp*G%R)/3._wp) + tau_e%R(dir_idx_tau(1) &
                                                  & ))/rho%R), &
                                                  & vel%L(dir_idx(1)) + sqrt(c%L*c%L + (((4._wp*G%L)/3._wp) &
                                                  & + tau_e%L(dir_idx_tau(1)))/rho%L))
                                        s_S = (pres%R - tau_e%R(dir_idx_tau(1)) - pres%L + tau_e%L(dir_idx_tau(1)) &
                                               & + rho%L*vel%L(dir_idx(1))*(s%L - vel%L(dir_idx(1))) - rho%R*vel%R(dir_idx(1)) &
                                               & *(s%R - vel%R(dir_idx(1))))/(rho%L*(s%L - vel%L(dir_idx(1))) - rho%R*(s%R &
                                               & - vel%R(dir_idx(1))))
                                    else
                                        s%L = min(vel%L(dir_idx(1)) - c%L, vel%R(dir_idx(1)) - c%R)
                                        s%R = max(vel%R(dir_idx(1)) + c%R, vel%L(dir_idx(1)) + c%L)
                                        s_S = (pres%R - pres%L + rho%L*vel%L(dir_idx(1))*(s%L - vel%L(dir_idx(1))) &
                                               & - rho%R*vel%R(dir_idx(1))*(s%R - vel%R(dir_idx(1))))/(rho%L*(s%L &
                                               & - vel%L(dir_idx(1))) - rho%R*(s%R - vel%R(dir_idx(1))))
                                    end if
                                else if (wave_speeds == 2) then
                                    pres_S%L = 5.e-1_wp*(pres%L + pres%R + rho_avg*c_avg*(vel%L(dir_idx(1)) - vel%R(dir_idx(1))))

                                    pres_S%R = pres_S%L

                                    ! Low Mach correction: Thornber et al. JCP (2008)
                                    Ms%L = max(1._wp, &
                                               & sqrt(1._wp + ((5.e-1_wp + gamma%L)/(1._wp + gamma%L))*(pres_S%L/pres%L - 1._wp) &
                                               & *pres%L/((pres%L + pi_inf%L/(1._wp + gamma%L)))))
                                    Ms%R = max(1._wp, &
                                               & sqrt(1._wp + ((5.e-1_wp + gamma%R)/(1._wp + gamma%R))*(pres_S%R/pres%R - 1._wp) &
                                               & *pres%R/((pres%R + pi_inf%R/(1._wp + gamma%R)))))

                                    s%L = vel%L(dir_idx(1)) - c%L*Ms%L
                                    s%R = vel%R(dir_idx(1)) + c%R*Ms%R

                                    s_S = 5.e-1_wp*((vel%L(dir_idx(1)) + vel%R(dir_idx(1))) + (pres%L - pres%R)/(rho_avg*c_avg))
                                end if

                                ! follows Einfeldt et al. s_M/P = min/max(0.,s%L/R)
                                s_M = min(0._wp, s%L); s_P = max(0._wp, s%R)

                                ! goes with q_star_L/R = xi%L/R * (variable) xi%L/R = ( ( s%L/R - u_L/R )/(s%L/R - s_star) )
                                xi%L = (s%L - vel%L(dir_idx(1)))/min(s%L - s_S, -sgm_eps)
                                xi%R = (s%R - vel%R(dir_idx(1)))/max(s%R - s_S, sgm_eps)
                                xi_m1%L = (s_S - vel%L(dir_idx(1)))/min(s%L - s_S, -sgm_eps)
                                xi_m1%R = (s_S - vel%R(dir_idx(1)))/max(s%R - s_S, sgm_eps)

                                ! goes with numerical star velocity in x/y/z directions xi_P/M = 0.5 +/m sgn(0.5,s_star)
                                xi_M = (5.e-1_wp + sign(0.5_wp, s_S))
                                xi_P = (5.e-1_wp - sign(0.5_wp, s_S))

                                ! goes with the numerical velocity in x/y/z directions xi_P/M (pressure) = min/max(0. sgn(1,sL/sR))
                                xi_MP = -min(0._wp, sign(1._wp, s%L))
                                xi_PP = max(0._wp, sign(1._wp, s%R))

                                E_star = xi_M*(E%L + xi_MP*(xi%L*(E%L + (s_S - vel%L(dir_idx(1)))*(rho%L*s_S + pres%L/(s%L &
                                               & - vel%L(dir_idx(1))))) - E%L)) + xi_P*(E%R + xi_PP*(xi%R*(E%R + (s_S &
                                               & - vel%R(dir_idx(1)))*(rho%R*s_S + pres%R/(s%R - vel%R(dir_idx(1))))) - E%R))
                                p_Star = xi_M*(pres%L + xi_MP*(rho%L*(s%L - vel%L(dir_idx(1)))*(s_S - vel%L(dir_idx(1))))) &
                                               & + xi_P*(pres%R + xi_PP*(rho%R*(s%R - vel%R(dir_idx(1)))*(s_S - vel%R(dir_idx(1)))))

                                rho_Star = xi_M*(rho%L*(xi_MP*xi%L + 1._wp - xi_MP)) + xi_P*(rho%R*(xi_PP*xi%R + 1._wp - xi_PP))

                                vel_K_Star = vel%L(dir_idx(1))*(1._wp - xi_MP) + xi_MP*vel%R(dir_idx(1)) + xi_MP*xi_PP*(s_S &
                                                   & - vel%R(dir_idx(1)))

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
                                                & i)*(vel%L(dir_idx(1)) + s_M*xi_m1%L) + xi_P*qR_prim_rsx_vf(${SF(' + 1')}$, &
                                                & i)*(vel%R(dir_idx(1)) + s_P*xi_m1%R)
                                end do

                                ! MOMENTUM FLUX. f = \rho u u - \sigma, q = \rho u, q_star = \xi * \rho*(s_star, v, w)
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    flux_rsx_vf(${SF('')}$, &
                                                & eqn_idx%cont%end + dir_idx(i)) = rho_Star*vel_K_Star*(dir_flg(dir_idx(i)) &
                                                & *vel_K_Star + (1._wp - dir_flg(dir_idx(i)))*(xi_M*vel%L(dir_idx(i)) &
                                                & + xi_P*vel%R(dir_idx(i)))) + dir_flg(dir_idx(i))*p_Star + (s_M/s%L)*(s_P/s%R) &
                                                & *dir_flg(dir_idx(i))*pcorr
                                end do

                                ! ENERGY FLUX. f = u*(E-\sigma), q = E, q_star = \xi*E+(s-u)(\rho s_star - \sigma/(s-u))
                                flux_rsx_vf(${SF('')}$, eqn_idx%E) = (E_star + p_Star)*vel_K_Star + (s_M/s%L)*(s_P/s%R)*pcorr*s_S

                                ! ELASTICITY. Elastic shear stress additions for the momentum and energy flux
                                if (elasticity) then
                                    flux_ene_e = 0._wp
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        ! MOMENTUM ELASTIC FLUX.
                                        flux_rsx_vf(${SF('')}$, eqn_idx%cont%end + dir_idx(i)) = flux_rsx_vf(${SF('')}$, &
                                                    & eqn_idx%cont%end + dir_idx(i)) - xi_M*tau_e%L(dir_idx_tau(i)) &
                                                    & - xi_P*tau_e%R(dir_idx_tau(i))
                                        ! ENERGY ELASTIC FLUX.
                                        flux_ene_e = flux_ene_e - xi_M*(vel%L(dir_idx(i))*tau_e%L(dir_idx_tau(i)) &
                                                                        & + s_M*(xi%L*((s_S - vel%L(i))*(tau_e%L(dir_idx_tau(i)) &
                                                                        & /(s%L - vel%L(i)))))) - xi_P*(vel%R(dir_idx(i)) &
                                                                        & *tau_e%R(dir_idx_tau(i)) + s_P*(xi%R*((s_S - vel%R(i)) &
                                                                        & *(tau_e%R(dir_idx_tau(i))/(s%R - vel%R(i))))))
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
                                                   & dir_idx(i)) = xi_M*(vel%L(dir_idx(i)) + dir_flg(dir_idx(i)) &
                                                   & *(s_S*(xi_MP*xi_m1%L + 1) - vel%L(dir_idx(i)))) + xi_P*(vel%R(dir_idx(i)) &
                                                   & + dir_flg(dir_idx(i))*(s_S*(xi_PP*xi_m1%R + 1) - vel%R(dir_idx(i))))
                                end do

                                ! INTERNAL ENERGIES ADVECTION FLUX. K-th pressure and velocity in preparation for the internal
                                ! energy flux
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    p_K_Star = xi_M*(xi_MP*((pres%L + pi_infs(i)/(1._wp + gammas(i)))*xi%L**(1._wp/gammas(i) &
                                                     & + 1._wp) - pi_infs(i)/(1._wp + gammas(i)) - pres%L) + pres%L) &
                                                     & + xi_P*(xi_PP*((pres%R + pi_infs(i)/(1._wp + gammas(i))) &
                                                     & *xi%R**(1._wp/gammas(i) + 1._wp) - pi_infs(i)/(1._wp + gammas(i)) - pres%R) &
                                                     & + pres%R)

                                    flux_rsx_vf(${SF('')}$, i + eqn_idx%int_en%beg - 1) = ((xi_M*qL_prim_rsx_vf(${SF('')}$, &
                                                & i + eqn_idx%adv%beg - 1) + xi_P*qR_prim_rsx_vf(${SF(' + 1')}$, &
                                                & i + eqn_idx%adv%beg - 1))*(gammas(i)*p_K_Star + pi_infs(i)) &
                                                & + (xi_M*qL_prim_rsx_vf(${SF('')}$, &
                                                & i + eqn_idx%cont%beg - 1) + xi_P*qR_prim_rsx_vf(${SF(' + 1')}$, &
                                                & i + eqn_idx%cont%beg - 1))*qvs(i))*vel_K_Star + (s_M/s%L)*(s_P/s%R) &
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
                                                    & eqn_idx%stress%beg - 1 + i) = xi_M*(s_S/(s%L - s_S))*(s%L*rho%L*tau_e%L(i) &
                                                    & - rho%L*vel%L(dir_idx(1))*tau_e%L(i)) + xi_P*(s_S/(s%R - s_S)) &
                                                    & *(s%R*rho%R*tau_e%R(i) - rho%R*vel%R(dir_idx(1))*tau_e%R(i))
                                    end do
                                end if

                                ! Hyperelastic reference map flux for material deformation tracking
                                if (hyperelasticity) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        flux_rsx_vf(${SF('')}$, &
                                                    & eqn_idx%xi%beg - 1 + i) = xi_M*(s_S/(s%L - s_S))*(s%L*rho%L*xi_field%L(i) &
                                                    & - rho%L*vel%L(dir_idx(1))*xi_field%L(i)) + xi_P*(s_S/(s%R - s_S)) &
                                                    & *(s%R*rho%R*xi_field%R(i) - rho%R*vel%R(dir_idx(1))*xi_field%R(i))
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
                else if (model_eqns == 4) then
                    ! 4-equation model (model_eqns=4): single pressure, velocity equilibrium
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[i, q, alpha_rho_L, alpha_rho_R, vel, alpha_L, alpha_R, nbub, rho, &
                                        & pres, E, H, Cp_avg, Cv_avg, T_avg, eps, c_sum_Yi_Phi, T, Y, MW, R_gas, Cp, Gamm, gamma, &
                                        & pi_inf, qv, qv_avg, c, G, rho_avg, H_avg, c_avg, gamma_avg, ptilde, vel_rms, &
                                        & vel_avg_rms, vel_tmp, Ms, pres_S, alpha_L_sum, alpha_R_sum, rho_Star, E_Star, p_Star, &
                                        & p_K_Star, vel_K_star, s, s_M, s_P, s_S, xi_M, xi_P, xi, xi_m1, xi_MP, xi_PP, Ys_L, &
                                        & Ys_R, Cp_iL, Cp_iR, Xs_L, Xs_R, Gamma_iL, Gamma_iR, Yi_avg, Phi_avg, h_iL, h_iR, h_avg_2]')
                    do l = ${Z_BND}$%beg, ${Z_BND}$%end
                        do k = ${Y_BND}$%beg, ${Y_BND}$%end
                            do j = ${X_BND}$%beg, ${X_BND}$%end
                                vel_rms%L = 0._wp; vel_rms%R = 0._wp
                                rho%L = 0._wp; rho%R = 0._wp
                                gamma%L = 0._wp; gamma%R = 0._wp
                                pi_inf%L = 0._wp; pi_inf%R = 0._wp
                                qv%L = 0._wp; qv%R = 0._wp

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, eqn_idx%cont%end
                                    alpha_rho_L(i) = qL_prim_rsx_vf(${SF('')}$, i)
                                    alpha_rho_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, i)
                                end do

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel%L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%cont%end + i)
                                    vel%R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%cont%end + i)
                                    vel_rms%L = vel_rms%L + vel%L(i)**2._wp
                                    vel_rms%R = vel_rms%R + vel%R(i)**2._wp
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

                                pres%L = qL_prim_rsx_vf(${SF('')}$, eqn_idx%E)
                                pres%R = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E)

                                E%L = gamma%L*pres%L + pi_inf%L + 5.e-1_wp*rho%L*vel_rms%L + qv%L
                                E%R = gamma%R*pres%R + pi_inf%R + 5.e-1_wp*rho%R*vel_rms%R + qv%R

                                H%L = (E%L + pres%L)/rho%L
                                H%R = (E%R + pres%R)/rho%R

                                @:compute_average_state()

                                call s_compute_speed_of_sound(pres%L, rho%L, gamma%L, pi_inf%L, H%L, alpha_L, vel_rms%L, 0._wp, &
                                                              & c%L, qv%L)

                                call s_compute_speed_of_sound(pres%R, rho%R, gamma%R, pi_inf%R, H%R, alpha_R, vel_rms%R, 0._wp, &
                                                              & c%R, qv%R)

                                !> The computation of c_avg does not require all the variables, and therefore the non '_avg'
                                ! variables are placeholders to call the subroutine.

                                call s_compute_speed_of_sound(pres%R, rho_avg, gamma_avg, pi_inf%R, H_avg, alpha_R, vel_avg_rms, &
                                                              & 0._wp, c_avg, qv_avg)

                                if (wave_speeds == 1) then
                                    s%L = min(vel%L(dir_idx(1)) - c%L, vel%R(dir_idx(1)) - c%R)
                                    s%R = max(vel%R(dir_idx(1)) + c%R, vel%L(dir_idx(1)) + c%L)

                                    s_S = (pres%R - pres%L + rho%L*vel%L(dir_idx(1))*(s%L - vel%L(dir_idx(1))) &
                                           & - rho%R*vel%R(dir_idx(1))*(s%R - vel%R(dir_idx(1))))/(rho%L*(s%L - vel%L(dir_idx(1))) &
                                           & - rho%R*(s%R - vel%R(dir_idx(1))))
                                else if (wave_speeds == 2) then
                                    pres_S%L = 5.e-1_wp*(pres%L + pres%R + rho_avg*c_avg*(vel%L(dir_idx(1)) - vel%R(dir_idx(1))))

                                    pres_S%R = pres_S%L

                                    ! Low Mach correction: Thornber et al. JCP (2008)
                                    Ms%L = max(1._wp, &
                                               & sqrt(1._wp + ((5.e-1_wp + gamma%L)/(1._wp + gamma%L))*(pres_S%L/pres%L - 1._wp) &
                                               & *pres%L/((pres%L + pi_inf%L/(1._wp + gamma%L)))))
                                    Ms%R = max(1._wp, &
                                               & sqrt(1._wp + ((5.e-1_wp + gamma%R)/(1._wp + gamma%R))*(pres_S%R/pres%R - 1._wp) &
                                               & *pres%R/((pres%R + pi_inf%R/(1._wp + gamma%R)))))

                                    s%L = vel%L(dir_idx(1)) - c%L*Ms%L
                                    s%R = vel%R(dir_idx(1)) + c%R*Ms%R

                                    s_S = 5.e-1_wp*((vel%L(dir_idx(1)) + vel%R(dir_idx(1))) + (pres%L - pres%R)/(rho_avg*c_avg))
                                end if

                                ! follows Einfeldt et al. s_M/P = min/max(0.,s%L/R)
                                s_M = min(0._wp, s%L); s_P = max(0._wp, s%R)

                                ! goes with q_star_L/R = xi%L/R * (variable) xi%L/R = ( ( s%L/R - u_L/R )/(s%L/R - s_star) )
                                xi%L = (s%L - vel%L(dir_idx(1)))/min(s%L - s_S, -sgm_eps)
                                xi%R = (s%R - vel%R(dir_idx(1)))/max(s%R - s_S, sgm_eps)
                                xi_m1%L = (s_S - vel%L(dir_idx(1)))/min(s%L - s_S, -sgm_eps)
                                xi_m1%R = (s_S - vel%R(dir_idx(1)))/max(s%R - s_S, sgm_eps)

                                ! goes with numerical velocity in x/y/z directions xi_P/M = 0.5 +/m sgn(0.5,s_star)
                                xi_M = (5.e-1_wp + sign(5.e-1_wp, s_S))
                                xi_P = (5.e-1_wp - sign(5.e-1_wp, s_S))

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, eqn_idx%cont%end
                                    flux_rsx_vf(${SF('')}$, &
                                                & i) = xi_M*alpha_rho_L(i)*(vel%L(dir_idx(1)) + s_M*xi_m1%L) + xi_P*alpha_rho_R(i) &
                                                & *(vel%R(dir_idx(1)) + s_P*xi_m1%R)
                                end do

                                ! Momentum flux. f = \rho u u + p I, q = \rho u, q_star = \xi * \rho*(s_star, v, w)
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    flux_rsx_vf(${SF('')}$, &
                                                & eqn_idx%cont%end + dir_idx(i)) = xi_M*(rho%L*(vel%L(dir_idx(1))*vel%L(dir_idx(i) &
                                                & ) + s_M*(xi%L*(dir_flg(dir_idx(i))*s_S + (1._wp - dir_flg(dir_idx(i))) &
                                                & *vel%L(dir_idx(i))) - vel%L(dir_idx(i)))) + dir_flg(dir_idx(i))*pres%L) &
                                                & + xi_P*(rho%R*(vel%R(dir_idx(1))*vel%R(dir_idx(i)) &
                                                & + s_P*(xi%R*(dir_flg(dir_idx(i))*s_S + (1._wp - dir_flg(dir_idx(i))) &
                                                & *vel%R(dir_idx(i))) - vel%R(dir_idx(i)))) + dir_flg(dir_idx(i))*pres%R)
                                end do

                                if (bubbles_euler) then
                                    ! Put p_tilde in
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        flux_rsx_vf(${SF('')}$, eqn_idx%cont%end + dir_idx(i)) = flux_rsx_vf(${SF('')}$, &
                                                    & eqn_idx%cont%end + dir_idx(i)) + xi_M*(dir_flg(dir_idx(i))*(-1._wp*ptilde%L) &
                                                    & ) + xi_P*(dir_flg(dir_idx(i))*(-1._wp*ptilde%R))
                                    end do
                                end if

                                flux_rsx_vf(${SF('')}$, eqn_idx%E) = 0._wp

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = eqn_idx%alf, eqn_idx%alf  ! only advect the void fraction
                                    flux_rsx_vf(${SF('')}$, i) = xi_M*qL_prim_rsx_vf(${SF('')}$, &
                                                & i)*(vel%L(dir_idx(1)) + s_M*xi_m1%L) + xi_P*qR_prim_rsx_vf(${SF(' + 1')}$, &
                                                & i)*(vel%R(dir_idx(1)) + s_P*xi_m1%R)
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
                                        flux_rsx_vf(${SF('')}$, i) = xi_M*nbub%L*qL_prim_rsx_vf(${SF('')}$, &
                                                    & i)*(vel%L(dir_idx(1)) + s_M*xi_m1%L) &
                                                    & + xi_P*nbub%R*qR_prim_rsx_vf(${SF(' + 1')}$, &
                                                    & i)*(vel%R(dir_idx(1)) + s_P*xi_m1%R)
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
                                                         & eqn_idx%cont%end + dir_idx(1)) = xi_M*(rho%L*(vel%L(dir_idx(1)) &
                                                         & *vel%L(dir_idx(1)) + s_M*(xi%L*(dir_flg(dir_idx(1))*s_S + (1._wp &
                                                         & - dir_flg(dir_idx(1)))*vel%L(dir_idx(1))) - vel%L(dir_idx(1))))) &
                                                         & + xi_P*(rho%R*(vel%R(dir_idx(1))*vel%R(dir_idx(1)) &
                                                         & + s_P*(xi%R*(dir_flg(dir_idx(1))*s_S + (1._wp - dir_flg(dir_idx(1))) &
                                                         & *vel%R(dir_idx(1))) - vel%R(dir_idx(1)))))
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
                                                         & eqn_idx%mom%beg + 1) = -xi_M*(rho%L*(vel%L(dir_idx(1))*vel%L(dir_idx(1) &
                                                         & ) + s_M*(xi%L*(dir_flg(dir_idx(1))*s_S + (1._wp - dir_flg(dir_idx(1))) &
                                                         & *vel%L(dir_idx(1))) - vel%L(dir_idx(1))))) &
                                                         & - xi_P*(rho%R*(vel%R(dir_idx(1))*vel%R(dir_idx(1)) &
                                                         & + s_P*(xi%R*(dir_flg(dir_idx(1))*s_S + (1._wp - dir_flg(dir_idx(1))) &
                                                         & *vel%R(dir_idx(1))) - vel%R(dir_idx(1)))))
                                        flux_gsrc_rsx_vf(${SF('')}$, eqn_idx%mom%end) = flux_rsx_vf(${SF('')}$, eqn_idx%mom%beg + 1)
                                    end if
                                #:endif
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                else if (model_eqns == 2 .and. bubbles_euler) then
                    ! 5-equation model with Euler-Euler bubble dynamics
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[i, q, R0_L, R0_R, V0_L, V0_R, P0_L, P0_R, pbw_L, pbw_R, vel, &
                                        & rho_avg, alpha_L, alpha_R, h_avg, gamma_avg, Re, pcorr, zcoef, rho, pres, E, H, gamma, &
                                        & pi_inf, qv, qv_avg, c, c_avg, vel_rms, vel_avg_rms, vel_tmp, Ms, pres_S, alpha_L_sum, &
                                        & alpha_R_sum, s, s_M, s_P, s_S, xi_M, xi_P, xi, xi_m1, xi_MP, xi_PP, nbub, PbwR3Lbar, &
                                        & PbwR3Rbar, R3Lbar, R3Rbar, R3V2Lbar, R3V2Rbar, Ys_L, Ys_R, Cp_iL, Cp_iR, Xs_L, Xs_R, &
                                        & Gamma_iL, Gamma_iR, Yi_avg, Phi_avg, h_iL, h_iR, h_avg_2]')
                    do l = ${Z_BND}$%beg, ${Z_BND}$%end
                        do k = ${Y_BND}$%beg, ${Y_BND}$%end
                            do j = ${X_BND}$%beg, ${X_BND}$%end
                                vel_rms%L = 0._wp; vel_rms%R = 0._wp
                                rho%L = 0._wp; rho%R = 0._wp
                                gamma%L = 0._wp; gamma%R = 0._wp
                                pi_inf%L = 0._wp; pi_inf%R = 0._wp
                                qv%L = 0._wp; qv%R = 0._wp

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    alpha_L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%E + i)
                                    alpha_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E + i)
                                end do

                                vel_rms%L = 0._wp; vel_rms%R = 0._wp

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel%L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%cont%end + i)
                                    vel%R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%cont%end + i)
                                    vel_rms%L = vel_rms%L + vel%L(i)**2._wp
                                    vel_rms%R = vel_rms%R + vel%R(i)**2._wp
                                end do

                                ! Retain this in the refactor
                                if (mpp_lim .and. (num_fluids > 2)) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        rho%L = rho%L + qL_prim_rsx_vf(${SF('')}$, i)
                                        gamma%L = gamma%L + qL_prim_rsx_vf(${SF('')}$, eqn_idx%E + i)*gammas(i)
                                        pi_inf%L = pi_inf%L + qL_prim_rsx_vf(${SF('')}$, eqn_idx%E + i)*pi_infs(i)
                                        qv%L = qv%L + qL_prim_rsx_vf(${SF('')}$, i)*qvs(i)
                                        rho%R = rho%R + qR_prim_rsx_vf(${SF(' + 1')}$, i)
                                        gamma%R = gamma%R + qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E + i)*gammas(i)
                                        pi_inf%R = pi_inf%R + qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E + i)*pi_infs(i)
                                        qv%R = qv%R + qR_prim_rsx_vf(${SF(' + 1')}$, i)*qvs(i)
                                    end do
                                else if (num_fluids > 2) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids - 1
                                        rho%L = rho%L + qL_prim_rsx_vf(${SF('')}$, i)
                                        gamma%L = gamma%L + qL_prim_rsx_vf(${SF('')}$, eqn_idx%E + i)*gammas(i)
                                        pi_inf%L = pi_inf%L + qL_prim_rsx_vf(${SF('')}$, eqn_idx%E + i)*pi_infs(i)
                                        qv%L = qv%L + qL_prim_rsx_vf(${SF('')}$, i)*qvs(i)
                                        rho%R = rho%R + qR_prim_rsx_vf(${SF(' + 1')}$, i)
                                        gamma%R = gamma%R + qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E + i)*gammas(i)
                                        pi_inf%R = pi_inf%R + qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E + i)*pi_infs(i)
                                        qv%R = qv%R + qR_prim_rsx_vf(${SF(' + 1')}$, i)*qvs(i)
                                    end do
                                else
                                    rho%L = qL_prim_rsx_vf(${SF('')}$, 1)
                                    gamma%L = gammas(1)
                                    pi_inf%L = pi_infs(1)
                                    qv%L = qvs(1)
                                    rho%R = qR_prim_rsx_vf(${SF(' + 1')}$, 1)
                                    gamma%R = gammas(1)
                                    pi_inf%R = pi_infs(1)
                                    qv%R = qvs(1)
                                end if

                                if (viscous) then
                                    if (num_fluids == 1) then  ! Need to consider case with num_fluids >= 2
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, 2
                                            Re%L(i) = dflt_real
                                            Re%R(i) = dflt_real

                                            if (Re_size(i) > 0) Re%L(i) = 0._wp
                                            if (Re_size(i) > 0) Re%R(i) = 0._wp

                                            $:GPU_LOOP(parallelism='[seq]')
                                            do q = 1, Re_size(i)
                                                Re%L(i) = (1._wp - qL_prim_rsx_vf(${SF('')}$, eqn_idx%E + Re_idx(i, &
                                                     & q)))/Res_gs(i, q) + Re%L(i)
                                                Re%R(i) = (1._wp - qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E + Re_idx(i, &
                                                     & q)))/Res_gs(i, q) + Re%R(i)
                                            end do

                                            Re%L(i) = 1._wp/max(Re%L(i), sgm_eps)
                                            Re%R(i) = 1._wp/max(Re%R(i), sgm_eps)
                                        end do
                                    end if
                                end if

                                pres%L = qL_prim_rsx_vf(${SF('')}$, eqn_idx%E)
                                pres%R = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E)

                                E%L = gamma%L*pres%L + pi_inf%L + 5.e-1_wp*rho%L*vel_rms%L
                                E%R = gamma%R*pres%R + pi_inf%R + 5.e-1_wp*rho%R*vel_rms%R

                                H%L = (E%L + pres%L)/rho%L
                                H%R = (E%R + pres%R)/rho%R

                                if (avg_state == 2) then
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
                                            nbub%L = qL_prim_rsx_vf(${SF('')}$, eqn_idx%n)
                                            nbub%R = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%n)
                                        else
                                            nbub%L = 0._wp
                                            nbub%R = 0._wp
                                            $:GPU_LOOP(parallelism='[seq]')
                                            do i = 1, nb
                                                nbub%L = nbub%L + (R0_L(i)**3._wp)*weight(i)
                                                nbub%R = nbub%R + (R0_R(i)**3._wp)*weight(i)
                                            end do

                                            nbub%L = (3._wp/(4._wp*pi))*qL_prim_rsx_vf(${SF('')}$, eqn_idx%E + num_fluids)/nbub%L
                                            nbub%R = (3._wp/(4._wp*pi))*qR_prim_rsx_vf(${SF(' + 1')}$, &
                                                      & eqn_idx%E + num_fluids)/nbub%R
                                        end if
                                    else
                                        ! nb stored in 0th moment of first R0 bin in variable conversion module
                                        nbub%L = qL_prim_rsx_vf(${SF('')}$, eqn_idx%bub%beg)
                                        nbub%R = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%bub%beg)
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

                                    rho_avg = 5.e-1_wp*(rho%L + rho%R)
                                    H_avg = 5.e-1_wp*(H%L + H%R)
                                    gamma_avg = 5.e-1_wp*(gamma%L + gamma%R)
                                    qv_avg = 5.e-1_wp*(qv%L + qv%R)
                                    vel_avg_rms = 0._wp

                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        vel_avg_rms = vel_avg_rms + (5.e-1_wp*(vel%L(i) + vel%R(i)))**2._wp
                                    end do
                                end if

                                call s_compute_speed_of_sound(pres%L, rho%L, gamma%L, pi_inf%L, H%L, alpha_L, vel_rms%L, 0._wp, &
                                                              & c%L, qv%L)

                                call s_compute_speed_of_sound(pres%R, rho%R, gamma%R, pi_inf%R, H%R, alpha_R, vel_rms%R, 0._wp, &
                                                              & c%R, qv%R)

                                !> The computation of c_avg does not require all the variables, and therefore the non '_avg'
                                ! variables are placeholders to call the subroutine.
                                call s_compute_speed_of_sound(pres%R, rho_avg, gamma_avg, pi_inf%R, H_avg, alpha_R, vel_avg_rms, &
                                                              & 0._wp, c_avg, qv_avg)

                                if (viscous) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, 2
                                        Re_avg_rsx_vf(${SF('')}$, i) = 2._wp/(1._wp/Re%L(i) + 1._wp/Re%R(i))
                                    end do
                                end if

                                ! Low Mach correction
                                if (low_Mach == 2) then
                                    @:compute_low_Mach_correction()
                                end if

                                if (wave_speeds == 1) then
                                    s%L = min(vel%L(dir_idx(1)) - c%L, vel%R(dir_idx(1)) - c%R)
                                    s%R = max(vel%R(dir_idx(1)) + c%R, vel%L(dir_idx(1)) + c%L)

                                    s_S = (pres%R - pres%L + rho%L*vel%L(dir_idx(1))*(s%L - vel%L(dir_idx(1))) &
                                           & - rho%R*vel%R(dir_idx(1))*(s%R - vel%R(dir_idx(1))))/(rho%L*(s%L - vel%L(dir_idx(1))) &
                                           & - rho%R*(s%R - vel%R(dir_idx(1))))
                                else if (wave_speeds == 2) then
                                    pres_S%L = 5.e-1_wp*(pres%L + pres%R + rho_avg*c_avg*(vel%L(dir_idx(1)) - vel%R(dir_idx(1))))

                                    pres_S%R = pres_S%L

                                    ! Low Mach correction: Thornber et al. JCP (2008)
                                    Ms%L = max(1._wp, &
                                               & sqrt(1._wp + ((5.e-1_wp + gamma%L)/(1._wp + gamma%L))*(pres_S%L/pres%L - 1._wp) &
                                               & *pres%L/((pres%L + pi_inf%L/(1._wp + gamma%L)))))
                                    Ms%R = max(1._wp, &
                                               & sqrt(1._wp + ((5.e-1_wp + gamma%R)/(1._wp + gamma%R))*(pres_S%R/pres%R - 1._wp) &
                                               & *pres%R/((pres%R + pi_inf%R/(1._wp + gamma%R)))))

                                    s%L = vel%L(dir_idx(1)) - c%L*Ms%L
                                    s%R = vel%R(dir_idx(1)) + c%R*Ms%R

                                    s_S = 5.e-1_wp*((vel%L(dir_idx(1)) + vel%R(dir_idx(1))) + (pres%L - pres%R)/(rho_avg*c_avg))
                                end if

                                ! follows Einfeldt et al. s_M/P = min/max(0.,s%L/R)
                                s_M = min(0._wp, s%L); s_P = max(0._wp, s%R)

                                ! goes with q_star_L/R = xi%L/R * (variable) xi%L/R = ( ( s%L/R - u_L/R )/(s%L/R - s_star) )
                                xi%L = (s%L - vel%L(dir_idx(1)))/min(s%L - s_S, -sgm_eps)
                                xi%R = (s%R - vel%R(dir_idx(1)))/max(s%R - s_S, sgm_eps)
                                xi_m1%L = (s_S - vel%L(dir_idx(1)))/min(s%L - s_S, -sgm_eps)
                                xi_m1%R = (s_S - vel%R(dir_idx(1)))/max(s%R - s_S, sgm_eps)

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
                                                & i)*(vel%L(dir_idx(1)) + s_M*xi_m1%L) + xi_P*qR_prim_rsx_vf(${SF(' + 1')}$, &
                                                & i)*(vel%R(dir_idx(1)) + s_P*xi_m1%R)
                                end do

                                if (bubbles_euler .and. (num_fluids > 1)) then
                                    ! Kill mass transport @ gas density
                                    flux_rsx_vf(${SF('')}$, eqn_idx%cont%end) = 0._wp
                                end if

                                ! Momentum flux. f = \rho u u + p I, q = \rho u, q_star = \xi * \rho*(s_star, v, w)

                                ! Include p_tilde

                                if (avg_state == 2) then
                                    if (alpha_L(num_fluids) < small_alf .or. R3Lbar < small_alf) then
                                        pres%L = pres%L - alpha_L(num_fluids)*pres%L
                                    else
                                        pres%L = pres%L - alpha_L(num_fluids)*(pres%L - PbwR3Lbar/R3Lbar - rho%L*R3V2Lbar/R3Lbar)
                                    end if

                                    if (alpha_R(num_fluids) < small_alf .or. R3Rbar < small_alf) then
                                        pres%R = pres%R - alpha_R(num_fluids)*pres%R
                                    else
                                        pres%R = pres%R - alpha_R(num_fluids)*(pres%R - PbwR3Rbar/R3Rbar - rho%R*R3V2Rbar/R3Rbar)
                                    end if
                                end if

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    flux_rsx_vf(${SF('')}$, &
                                                & eqn_idx%cont%end + dir_idx(i)) = xi_M*(rho%L*(vel%L(dir_idx(1))*vel%L(dir_idx(i) &
                                                & ) + s_M*(xi%L*(dir_flg(dir_idx(i))*s_S + (1._wp - dir_flg(dir_idx(i))) &
                                                & *vel%L(dir_idx(i))) - vel%L(dir_idx(i)))) + dir_flg(dir_idx(i))*(pres%L)) &
                                                & + xi_P*(rho%R*(vel%R(dir_idx(1))*vel%R(dir_idx(i)) &
                                                & + s_P*(xi%R*(dir_flg(dir_idx(i))*s_S + (1._wp - dir_flg(dir_idx(i))) &
                                                & *vel%R(dir_idx(i))) - vel%R(dir_idx(i)))) + dir_flg(dir_idx(i))*(pres%R)) &
                                                & + (s_M/s%L)*(s_P/s%R)*dir_flg(dir_idx(i))*pcorr
                                end do

                                ! Energy flux. f = u*(E+p), q = E, q_star = \xi*E+(s-u)(\rho s_star + p/(s-u))
                                flux_rsx_vf(${SF('')}$, &
                                            & eqn_idx%E) = xi_M*(vel%L(dir_idx(1))*(E%L + pres%L) + s_M*(xi%L*(E%L + (s_S &
                                            & - vel%L(dir_idx(1)))*(rho%L*s_S + (pres%L)/(s%L - vel%L(dir_idx(1))))) - E%L)) &
                                            & + xi_P*(vel%R(dir_idx(1))*(E%R + pres%R) + s_P*(xi%R*(E%R + (s_S - vel%R(dir_idx(1)) &
                                            & )*(rho%R*s_S + (pres%R)/(s%R - vel%R(dir_idx(1))))) - E%R)) + (s_M/s%L)*(s_P/s%R) &
                                            & *pcorr*s_S

                                ! Volume fraction flux
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = eqn_idx%adv%beg, eqn_idx%adv%end
                                    flux_rsx_vf(${SF('')}$, i) = xi_M*qL_prim_rsx_vf(${SF('')}$, &
                                                & i)*(vel%L(dir_idx(1)) + s_M*xi_m1%L) + xi_P*qR_prim_rsx_vf(${SF(' + 1')}$, &
                                                & i)*(vel%R(dir_idx(1)) + s_P*xi_m1%R)
                                end do

                                ! Advection velocity source: interface velocity for volume fraction transport
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel_src_rsx_vf(${SF('')}$, &
                                                   & dir_idx(i)) = xi_M*(vel%L(dir_idx(i)) + dir_flg(dir_idx(i))*s_M*xi_m1%L) &
                                                   & + xi_P*(vel%R(dir_idx(i)) + dir_flg(dir_idx(i))*s_P*xi_m1%R)

                                    ! IF ( (model_eqns == 4) .or. (num_fluids==1) ) vel_src_rs_vf(dir_idx(i))%sf(j,k,l) = 0._wp
                                end do

                                flux_src_rsx_vf(${SF('')}$, eqn_idx%adv%beg) = vel_src_rsx_vf(${SF('')}$, dir_idx(1))

                                ! Add advection flux for bubble variables
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = eqn_idx%bub%beg, eqn_idx%bub%end
                                    flux_rsx_vf(${SF('')}$, i) = xi_M*nbub%L*qL_prim_rsx_vf(${SF('')}$, &
                                                & i)*(vel%L(dir_idx(1)) + s_M*xi_m1%L) &
                                                & + xi_P*nbub%R*qR_prim_rsx_vf(${SF(' + 1')}$, i)*(vel%R(dir_idx(1)) + s_P*xi_m1%R)
                                end do

                                if (qbmm) then
                                    flux_rsx_vf(${SF('')}$, &
                                                & eqn_idx%bub%beg) = xi_M*nbub%L*(vel%L(dir_idx(1)) + s_M*xi_m1%L) &
                                                & + xi_P*nbub%R*(vel%R(dir_idx(1)) + s_P*xi_m1%R)
                                end if

                                if (adv_n) then
                                    flux_rsx_vf(${SF('')}$, &
                                                & eqn_idx%n) = xi_M*nbub%L*(vel%L(dir_idx(1)) + s_M*xi_m1%L) &
                                                & + xi_P*nbub%R*(vel%R(dir_idx(1)) + s_P*xi_m1%R)
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
                                                         & eqn_idx%cont%end + dir_idx(1)) = xi_M*(rho%L*(vel%L(dir_idx(1)) &
                                                         & *vel%L(dir_idx(1)) + s_M*(xi%L*(dir_flg(dir_idx(1))*s_S + (1._wp &
                                                         & - dir_flg(dir_idx(1)))*vel%L(dir_idx(1))) - vel%L(dir_idx(1))))) &
                                                         & + xi_P*(rho%R*(vel%R(dir_idx(1))*vel%R(dir_idx(1)) &
                                                         & + s_P*(xi%R*(dir_flg(dir_idx(1))*s_S + (1._wp - dir_flg(dir_idx(1))) &
                                                         & *vel%R(dir_idx(1))) - vel%R(dir_idx(1)))))
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
                                                         & eqn_idx%mom%beg + 1) = -xi_M*(rho%L*(vel%L(dir_idx(1))*vel%L(dir_idx(1) &
                                                         & ) + s_M*(xi%L*(dir_flg(dir_idx(1))*s_S + (1._wp - dir_flg(dir_idx(1))) &
                                                         & *vel%L(dir_idx(1))) - vel%L(dir_idx(1))))) &
                                                         & - xi_P*(rho%R*(vel%R(dir_idx(1))*vel%R(dir_idx(1)) &
                                                         & + s_P*(xi%R*(dir_flg(dir_idx(1))*s_S + (1._wp - dir_flg(dir_idx(1))) &
                                                         & *vel%R(dir_idx(1))) - vel%R(dir_idx(1)))))
                                        flux_gsrc_rsx_vf(${SF('')}$, eqn_idx%mom%end) = flux_rsx_vf(${SF('')}$, eqn_idx%mom%beg + 1)
                                    end if
                                #:endif
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                else
                    ! 5-equation model (model_eqns=2): mixture total energy, volume fraction advection
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[Re_max, i, q, T, vel_rms, pres, rho, gamma, pi_inf, qv, &
                                        & alpha_L_sum, alpha_R_sum, E, MW, R_gas, Cp, Cv, Gamm, Y, H, qv_avg, rho_avg, gamma_avg, &
                                        & H_avg, c, c_avg, s_P, s_M, xi_P, xi_M, xi, xi_m1, Ms, pres_S, vel, Re, alpha_L, &
                                        & alpha_R, s, s_S, vel_avg_rms, pcorr, zcoef, vel_tmp, Ys_L, Ys_R, Xs_L, Xs_R, Gamma_iL, &
                                        & Gamma_iR, Cp_iL, Cp_iR, tau_e, xi_field, Yi_avg, Phi_avg, h_iL, h_iR, h_avg_2, G]', &
                                        & copyin='[is1, is2, is3]')
                    do l = ${Z_BND}$%beg, ${Z_BND}$%end
                        do k = ${Y_BND}$%beg, ${Y_BND}$%end
                            do j = ${X_BND}$%beg, ${X_BND}$%end
                                vel_rms%L = 0._wp; vel_rms%R = 0._wp
                                rho%L = 0._wp; rho%R = 0._wp
                                gamma%L = 0._wp; gamma%R = 0._wp
                                pi_inf%L = 0._wp; pi_inf%R = 0._wp
                                qv%L = 0._wp; qv%R = 0._wp
                                alpha_L_sum = 0._wp; alpha_R_sum = 0._wp

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    alpha_L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%E + i)
                                    alpha_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E + i)
                                end do

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel%L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%cont%end + i)
                                    vel%R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%cont%end + i)
                                    vel_rms%L = vel_rms%L + vel%L(i)**2._wp
                                    vel_rms%R = vel_rms%R + vel%R(i)**2._wp
                                end do

                                pres%L = qL_prim_rsx_vf(${SF('')}$, eqn_idx%E)
                                pres%R = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E)

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

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    rho%L = rho%L + qL_prim_rsx_vf(${SF('')}$, i)
                                    gamma%L = gamma%L + qL_prim_rsx_vf(${SF('')}$, eqn_idx%E + i)*gammas(i)
                                    pi_inf%L = pi_inf%L + qL_prim_rsx_vf(${SF('')}$, eqn_idx%E + i)*pi_infs(i)
                                    qv%L = qv%L + qL_prim_rsx_vf(${SF('')}$, i)*qvs(i)

                                    rho%R = rho%R + qR_prim_rsx_vf(${SF(' + 1')}$, i)
                                    gamma%R = gamma%R + qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E + i)*gammas(i)
                                    pi_inf%R = pi_inf%R + qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%E + i)*pi_infs(i)
                                    qv%R = qv%R + qR_prim_rsx_vf(${SF(' + 1')}$, i)*qvs(i)
                                end do

                                Re_max = 0
                                if (Re_size(1) > 0) Re_max = 1
                                if (Re_size(2) > 0) Re_max = 2

                                if (viscous) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, Re_max
                                        Re%L(i) = 0._wp
                                        Re%R(i) = 0._wp

                                        $:GPU_LOOP(parallelism='[seq]')
                                        do q = 1, Re_size(i)
                                            Re%L(i) = alpha_L(Re_idx(i, q))/Res_gs(i, q) + Re%L(i)
                                            Re%R(i) = alpha_R(Re_idx(i, q))/Res_gs(i, q) + Re%R(i)
                                        end do

                                        Re%L(i) = 1._wp/max(Re%L(i), sgm_eps)
                                        Re%R(i) = 1._wp/max(Re%R(i), sgm_eps)
                                    end do
                                end if

                                if (chemistry) then
                                    c_sum_Yi_Phi = 0.0_wp
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = eqn_idx%species%beg, eqn_idx%species%end
                                        Ys_L(i - eqn_idx%species%beg + 1) = qL_prim_rsx_vf(${SF('')}$, i)
                                        Ys_R(i - eqn_idx%species%beg + 1) = qR_prim_rsx_vf(${SF(' + 1')}$, i)
                                    end do

                                    call get_mixture_molecular_weight(Ys_L, MW%L)
                                    call get_mixture_molecular_weight(Ys_R, MW%R)

                                    Xs_L(:) = Ys_L(:)*MW%L/molecular_weights(:)
                                    Xs_R(:) = Ys_R(:)*MW%R/molecular_weights(:)

                                    R_gas%L = gas_constant/MW%L
                                    R_gas%R = gas_constant/MW%R

                                    T%L = pres%L/rho%L/R_gas%L
                                    T%R = pres%R/rho%R/R_gas%R

                                    call get_species_specific_heats_r(T%L, Cp_iL)
                                    call get_species_specific_heats_r(T%R, Cp_iR)

                                    if (chem_params%gamma_method == 1) then
                                        !> gamma_method = 1: Ref. Section 2.3.1 Formulation of doi:10.7907/ZKW8-ES97.
                                        Gamma_iL = Cp_iL/(Cp_iL - 1.0_wp)
                                        Gamma_iR = Cp_iR/(Cp_iR - 1.0_wp)

                                        gamma%L = sum(Xs_L(:)/(Gamma_iL(:) - 1.0_wp))
                                        gamma%R = sum(Xs_R(:)/(Gamma_iR(:) - 1.0_wp))
                                    else if (chem_params%gamma_method == 2) then
                                        !> gamma_method = 2: c_p / c_v where c_p, c_v are specific heats.
                                        call get_mixture_specific_heat_cp_mass(T%L, Ys_L, Cp%L)
                                        call get_mixture_specific_heat_cp_mass(T%R, Ys_R, Cp%R)
                                        call get_mixture_specific_heat_cv_mass(T%L, Ys_L, Cv%L)
                                        call get_mixture_specific_heat_cv_mass(T%R, Ys_R, Cv%R)

                                        Gamm%L = Cp%L/Cv%L; Gamm%R = Cp%R/Cv%R
                                        gamma%L = 1.0_wp/(Gamm%L - 1.0_wp); gamma%R = 1.0_wp/(Gamm%R - 1.0_wp)
                                    end if

                                    call get_mixture_energy_mass(T%L, Ys_L, E%L)
                                    call get_mixture_energy_mass(T%R, Ys_R, E%R)

                                    E%L = rho%L*E%L + 5.e-1*rho%L*vel_rms%L
                                    E%R = rho%R*E%R + 5.e-1*rho%R*vel_rms%R
                                    H%L = (E%L + pres%L)/rho%L
                                    H%R = (E%R + pres%R)/rho%R
                                else
                                    E%L = gamma%L*pres%L + pi_inf%L + 5.e-1*rho%L*vel_rms%L + qv%L
                                    E%R = gamma%R*pres%R + pi_inf%R + 5.e-1*rho%R*vel_rms%R + qv%R

                                    H%L = (E%L + pres%L)/rho%L
                                    H%R = (E%R + pres%R)/rho%R
                                end if

                                ! ENERGY ADJUSTMENTS FOR HYPOELASTIC ENERGY
                                if (hypoelasticity) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, eqn_idx%stress%end - eqn_idx%stress%beg + 1
                                        tau_e%L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%stress%beg - 1 + i)
                                        tau_e%R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%stress%beg - 1 + i)
                                    end do
                                    G%L = 0._wp
                                    G%R = 0._wp
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        G%L = G%L + alpha_L(i)*Gs_rs(i)
                                        G%R = G%R + alpha_R(i)*Gs_rs(i)
                                    end do
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, eqn_idx%stress%end - eqn_idx%stress%beg + 1
                                        ! Elastic contribution to energy if G large enough
                                        if ((G%L > verysmall) .and. (G%R > verysmall)) then
                                            E%L = E%L + (tau_e%L(i)*tau_e%L(i))/(4._wp*G%L)
                                            E%R = E%R + (tau_e%R(i)*tau_e%R(i))/(4._wp*G%R)
                                            ! Additional terms in 2D and 3D
                                            if ((i == 2) .or. (i == 4) .or. (i == 5)) then
                                                E%L = E%L + (tau_e%L(i)*tau_e%L(i))/(4._wp*G%L)
                                                E%R = E%R + (tau_e%R(i)*tau_e%R(i))/(4._wp*G%R)
                                            end if
                                        end if
                                    end do
                                end if

                                ! Hyperelastic stress contribution: strain energy added to total energy
                                if (hyperelasticity) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        xi_field%L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%xi%beg - 1 + i)
                                        xi_field%R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%xi%beg - 1 + i)
                                    end do
                                    G%L = 0._wp
                                    G%R = 0._wp
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        ! Mixture left and right shear modulus
                                        G%L = G%L + alpha_L(i)*Gs_rs(i)
                                        G%R = G%R + alpha_R(i)*Gs_rs(i)
                                    end do
                                    ! Elastic contribution to energy if G large enough
                                    if (G%L > verysmall .and. G%R > verysmall) then
                                        E%L = E%L + G%L*qL_prim_rsx_vf(${SF('')}$, eqn_idx%xi%end + 1)
                                        E%R = E%R + G%R*qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%xi%end + 1)
                                    end if
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, b_size - 1
                                        tau_e%L(i) = qL_prim_rsx_vf(${SF('')}$, eqn_idx%stress%beg - 1 + i)
                                        tau_e%R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%stress%beg - 1 + i)
                                    end do
                                end if

                                H%L = (E%L + pres%L)/rho%L
                                H%R = (E%R + pres%R)/rho%R

                                @:compute_average_state()

                                call s_compute_speed_of_sound(pres%L, rho%L, gamma%L, pi_inf%L, H%L, alpha_L, vel_rms%L, 0._wp, &
                                                              & c%L, qv%L)

                                call s_compute_speed_of_sound(pres%R, rho%R, gamma%R, pi_inf%R, H%R, alpha_R, vel_rms%R, 0._wp, &
                                                              & c%R, qv%R)

                                !> The computation of c_avg does not require all the variables, and therefore the non '_avg'
                                !  variables are placeholders to call the subroutine.
                                call s_compute_speed_of_sound(pres%R, rho_avg, gamma_avg, pi_inf%R, H_avg, alpha_R, vel_avg_rms, &
                                                              & c_sum_Yi_Phi, c_avg, qv_avg)

                                if (viscous) then
                                    if (chemistry) then
                                        call compute_viscosity_and_inversion(T%L, Ys_L, T%R, Ys_R, Re%L(1), Re%R(1))
                                    end if
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, 2
                                        Re_avg_rsx_vf(${SF('')}$, i) = 2._wp/(1._wp/Re%L(i) + 1._wp/Re%R(i))
                                    end do
                                end if

                                ! Low Mach correction
                                if (low_Mach == 2) then
                                    @:compute_low_Mach_correction()
                                end if

                                if (wave_speeds == 1) then
                                    if (elasticity) then
                                        ! Elastic wave speed, Rodriguez et al. JCP (2019)
                                        s%L = min(vel%L(dir_idx(1)) - sqrt(c%L*c%L + (((4._wp*G%L)/3._wp) + tau_e%L(dir_idx_tau(1) &
                                                  & ))/rho%L), &
                                                  & vel%R(dir_idx(1)) - sqrt(c%R*c%R + (((4._wp*G%R)/3._wp) &
                                                  & + tau_e%R(dir_idx_tau(1)))/rho%R))
                                        s%R = max(vel%R(dir_idx(1)) + sqrt(c%R*c%R + (((4._wp*G%R)/3._wp) + tau_e%R(dir_idx_tau(1) &
                                                  & ))/rho%R), &
                                                  & vel%L(dir_idx(1)) + sqrt(c%L*c%L + (((4._wp*G%L)/3._wp) &
                                                  & + tau_e%L(dir_idx_tau(1)))/rho%L))
                                        s_S = (pres%R - tau_e%R(dir_idx_tau(1)) - pres%L + tau_e%L(dir_idx_tau(1)) &
                                               & + rho%L*vel%L(dir_idx(1))*(s%L - vel%L(dir_idx(1))) - rho%R*vel%R(dir_idx(1)) &
                                               & *(s%R - vel%R(dir_idx(1))))/(rho%L*(s%L - vel%L(dir_idx(1))) - rho%R*(s%R &
                                               & - vel%R(dir_idx(1))))
                                    else
                                        s%L = min(vel%L(dir_idx(1)) - c%L, vel%R(dir_idx(1)) - c%R)
                                        s%R = max(vel%R(dir_idx(1)) + c%R, vel%L(dir_idx(1)) + c%L)
                                        s_S = (pres%R - pres%L + rho%L*vel%L(dir_idx(1))*(s%L - vel%L(dir_idx(1))) &
                                               & - rho%R*vel%R(dir_idx(1))*(s%R - vel%R(dir_idx(1))))/(rho%L*(s%L &
                                               & - vel%L(dir_idx(1))) - rho%R*(s%R - vel%R(dir_idx(1))))
                                    end if
                                else if (wave_speeds == 2) then
                                    pres_S%L = 5.e-1_wp*(pres%L + pres%R + rho_avg*c_avg*(vel%L(dir_idx(1)) - vel%R(dir_idx(1))))

                                    pres_S%R = pres_S%L

                                    ! Low Mach correction: Thornber et al. JCP (2008)
                                    Ms%L = max(1._wp, &
                                               & sqrt(1._wp + ((5.e-1_wp + gamma%L)/(1._wp + gamma%L))*(pres_S%L/pres%L - 1._wp) &
                                               & *pres%L/((pres%L + pi_inf%L/(1._wp + gamma%L)))))
                                    Ms%R = max(1._wp, &
                                               & sqrt(1._wp + ((5.e-1_wp + gamma%R)/(1._wp + gamma%R))*(pres_S%R/pres%R - 1._wp) &
                                               & *pres%R/((pres%R + pi_inf%R/(1._wp + gamma%R)))))

                                    s%L = vel%L(dir_idx(1)) - c%L*Ms%L
                                    s%R = vel%R(dir_idx(1)) + c%R*Ms%R

                                    s_S = 5.e-1_wp*((vel%L(dir_idx(1)) + vel%R(dir_idx(1))) + (pres%L - pres%R)/(rho_avg*c_avg))
                                end if

                                ! follows Einfeldt et al. s_M/P = min/max(0.,s%L/R)
                                s_M = min(0._wp, s%L); s_P = max(0._wp, s%R)

                                ! goes with q_star_L/R = xi%L/R * (variable) xi%L/R = ( ( s%L/R - u_L/R )/(s%L/R - s_star) )
                                xi%L = (s%L - vel%L(dir_idx(1)))/min(s%L - s_S, -sgm_eps)
                                xi%R = (s%R - vel%R(dir_idx(1)))/max(s%R - s_S, sgm_eps)
                                ! xi%L/R - 1 = (s_S - u_L/R)/(s%L/R - s_star): avoids cancellation when xi \approx 1
                                xi_m1%L = (s_S - vel%L(dir_idx(1)))/min(s%L - s_S, -sgm_eps)
                                xi_m1%R = (s_S - vel%R(dir_idx(1)))/max(s%R - s_S, sgm_eps)

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
                                                & i)*(vel%L(dir_idx(1)) + s_M*xi_m1%L) + xi_P*qR_prim_rsx_vf(${SF(' + 1')}$, &
                                                & i)*(vel%R(dir_idx(1)) + s_P*xi_m1%R)
                                end do

                                ! MOMENTUM FLUX. f = \rho u u - \sigma, q = \rho u, q_star = \xi * \rho*(s_star, v, w) identity:
                                ! xi*(dir_flg*s_S+(1-dir_flg)*u_i)-u_i = (dir_flg*s%L/R+(1-dir_flg)*u_i)*xi_m1
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    flux_rsx_vf(${SF('')}$, &
                                                & eqn_idx%cont%end + dir_idx(i)) = xi_M*(rho%L*(vel%L(dir_idx(1))*vel%L(dir_idx(i) &
                                                & ) + s_M*(dir_flg(dir_idx(i))*s%L + (1._wp - dir_flg(dir_idx(i))) &
                                                & *vel%L(dir_idx(i)))*xi_m1%L) + dir_flg(dir_idx(i))*(pres%L)) &
                                                & + xi_P*(rho%R*(vel%R(dir_idx(1))*vel%R(dir_idx(i)) + s_P*(dir_flg(dir_idx(i)) &
                                                & *s%R + (1._wp - dir_flg(dir_idx(i)))*vel%R(dir_idx(i)))*xi_m1%R) &
                                                & + dir_flg(dir_idx(i))*(pres%R)) + (s_M/s%L)*(s_P/s%R)*dir_flg(dir_idx(i))*pcorr
                                end do

                                ! ENERGY FLUX. f = u*(E-\sigma), q = E, q_star = \xi*E+(s-u)(\rho s_star - \sigma/(s-u))
                                ! xi*(E+expr)-E = E*xi_m1 + xi*expr avoids E*(xi-1) cancellation
                                flux_rsx_vf(${SF('')}$, &
                                            & eqn_idx%E) = xi_M*(vel%L(dir_idx(1))*(E%L + pres%L) + s_M*(E%L*xi_m1%L + xi%L*(s_S &
                                            & - vel%L(dir_idx(1)))*(rho%L*s_S + pres%L/(s%L - vel%L(dir_idx(1)))))) &
                                            & + xi_P*(vel%R(dir_idx(1))*(E%R + pres%R) + s_P*(E%R*xi_m1%R + xi%R*(s_S &
                                            & - vel%R(dir_idx(1)))*(rho%R*s_S + pres%R/(s%R - vel%R(dir_idx(1)))))) + (s_M/s%L) &
                                            & *(s_P/s%R)*pcorr*s_S

                                ! ELASTICITY. Elastic shear stress additions for the momentum and energy flux
                                if (elasticity) then
                                    flux_ene_e = 0._wp
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        ! MOMENTUM ELASTIC FLUX.
                                        flux_rsx_vf(${SF('')}$, eqn_idx%cont%end + dir_idx(i)) = flux_rsx_vf(${SF('')}$, &
                                                    & eqn_idx%cont%end + dir_idx(i)) - xi_M*tau_e%L(dir_idx_tau(i)) &
                                                    & - xi_P*tau_e%R(dir_idx_tau(i))
                                        ! ENERGY ELASTIC FLUX.
                                        flux_ene_e = flux_ene_e - xi_M*(vel%L(dir_idx(i))*tau_e%L(dir_idx_tau(i)) &
                                                                        & + s_M*(xi%L*((s_S - vel%L(i))*(tau_e%L(dir_idx_tau(i)) &
                                                                        & /(s%L - vel%L(i)))))) - xi_P*(vel%R(dir_idx(i)) &
                                                                        & *tau_e%R(dir_idx_tau(i)) + s_P*(xi%R*((s_S - vel%R(i)) &
                                                                        & *(tau_e%R(dir_idx_tau(i))/(s%R - vel%R(i))))))
                                    end do
                                    flux_rsx_vf(${SF('')}$, eqn_idx%E) = flux_rsx_vf(${SF('')}$, eqn_idx%E) + flux_ene_e
                                end if

                                ! HYPOELASTIC STRESS EVOLUTION FLUX.
                                if (hypoelasticity) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, eqn_idx%stress%end - eqn_idx%stress%beg + 1
                                        flux_rsx_vf(${SF('')}$, &
                                                    & eqn_idx%stress%beg - 1 + i) = xi_M*(s_S/(s%L - s_S))*(s%L*rho%L*tau_e%L(i) &
                                                    & - rho%L*vel%L(dir_idx(1))*tau_e%L(i)) + xi_P*(s_S/(s%R - s_S)) &
                                                    & *(s%R*rho%R*tau_e%R(i) - rho%R*vel%R(dir_idx(1))*tau_e%R(i))
                                    end do
                                end if

                                ! VOLUME FRACTION FLUX.
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = eqn_idx%adv%beg, eqn_idx%adv%end
                                    flux_rsx_vf(${SF('')}$, i) = xi_M*qL_prim_rsx_vf(${SF('')}$, &
                                                & i)*(vel%L(dir_idx(1)) + s_M*xi_m1%L) + xi_P*qR_prim_rsx_vf(${SF(' + 1')}$, &
                                                & i)*(vel%R(dir_idx(1)) + s_P*xi_m1%R)
                                end do

                                ! VOLUME FRACTION SOURCE FLUX.
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel_src_rsx_vf(${SF('')}$, &
                                                   & dir_idx(i)) = xi_M*(vel%L(dir_idx(i)) + dir_flg(dir_idx(i))*s_M*xi_m1%L) &
                                                   & + xi_P*(vel%R(dir_idx(i)) + dir_flg(dir_idx(i))*s_P*xi_m1%R)
                                end do

                                ! COLOR FUNCTION FLUX
                                if (surface_tension) then
                                    flux_rsx_vf(${SF('')}$, eqn_idx%c) = xi_M*qL_prim_rsx_vf(${SF('')}$, &
                                                & eqn_idx%c)*(vel%L(dir_idx(1)) + s_M*xi_m1%L) &
                                                & + xi_P*qR_prim_rsx_vf(${SF(' + 1')}$, eqn_idx%c)*(vel%R(dir_idx(1)) + s_P*xi_m1%R)
                                end if

                                ! Hyperelastic reference map flux for material deformation tracking
                                if (hyperelasticity) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        flux_rsx_vf(${SF('')}$, &
                                                    & eqn_idx%xi%beg - 1 + i) = xi_M*(s_S/(s%L - s_S))*(s%L*rho%L*xi_field%L(i) &
                                                    & - rho%L*vel%L(dir_idx(1))*xi_field%L(i)) + xi_P*(s_S/(s%R - s_S)) &
                                                    & *(s%R*rho%R*xi_field%R(i) - rho%R*vel%R(dir_idx(1))*xi_field%R(i))
                                    end do
                                end if

                                flux_src_rsx_vf(${SF('')}$, eqn_idx%adv%beg) = vel_src_rsx_vf(${SF('')}$, dir_idx(1))

                                if (chemistry) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = eqn_idx%species%beg, eqn_idx%species%end
                                        Y%L = qL_prim_rsx_vf(${SF('')}$, i)
                                        Y%R = qR_prim_rsx_vf(${SF(' + 1')}$, i)

                                        flux_rsx_vf(${SF('')}$, &
                                                    & i) = xi_M*rho%L*Y%L*(vel%L(dir_idx(1)) + s_M*xi_m1%L) &
                                                    & + xi_P*rho%R*Y%R*(vel%R(dir_idx(1)) + s_P*xi_m1%R)
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
                                                         & eqn_idx%cont%end + dir_idx(1)) = xi_M*(rho%L*(vel%L(dir_idx(1)) &
                                                         & *vel%L(dir_idx(1)) + s_M*(xi%L*(dir_flg(dir_idx(1))*s_S + (1._wp &
                                                         & - dir_flg(dir_idx(1)))*vel%L(dir_idx(1))) - vel%L(dir_idx(1))))) &
                                                         & + xi_P*(rho%R*(vel%R(dir_idx(1))*vel%R(dir_idx(1)) &
                                                         & + s_P*(xi%R*(dir_flg(dir_idx(1))*s_S + (1._wp - dir_flg(dir_idx(1))) &
                                                         & *vel%R(dir_idx(1))) - vel%R(dir_idx(1)))))
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
                                                         & eqn_idx%mom%beg + 1) = -xi_M*(rho%L*(vel%L(dir_idx(1))*vel%L(dir_idx(1) &
                                                         & ) + s_M*(xi%L*(dir_flg(dir_idx(1))*s_S + (1._wp - dir_flg(dir_idx(1))) &
                                                         & *vel%L(dir_idx(1))) - vel%L(dir_idx(1))))) &
                                                         & - xi_P*(rho%R*(vel%R(dir_idx(1))*vel%R(dir_idx(1)) &
                                                         & + s_P*(xi%R*(dir_flg(dir_idx(1))*s_S + (1._wp - dir_flg(dir_idx(1))) &
                                                         & *vel%R(dir_idx(1))) - vel%R(dir_idx(1)))))
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
                                                   & dqR_prim_dz_vf(eqn_idx%mom%beg:eqn_idx%mom%end), flux_src_vf, norm_dir, ix, &
                                                   & iy, iz)
            else
                call s_compute_viscous_source_flux(q_prim_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                                                   & dqL_prim_dx_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                                                   & dqL_prim_dy_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                                                   & dqL_prim_dz_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                                                   & q_prim_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                                                   & dqR_prim_dx_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                                                   & dqR_prim_dy_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                                                   & dqR_prim_dz_vf(eqn_idx%mom%beg:eqn_idx%mom%end), flux_src_vf, norm_dir, ix, &
                                                   & iy, iz)
            end if
        end if

        if (surface_tension) then
            call s_compute_capillary_source_flux(vel_src_rsx_vf, flux_src_vf, norm_dir, isx, isy, isz)
        end if

        call s_finalize_riemann_solver(flux_vf, flux_src_vf, flux_gsrc_vf, norm_dir)

    end subroutine s_hllc_riemann_solver

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
        type(riemann_states)      :: s, pTot
        real(wp)                  :: p_star, s_M, s_starL, s_starR, denom_ds, sign_Bx
        type(riemann_states)      :: rho_star, E_star, v_star, w_star, sqrt_rho_star, E_double_lr
        type(riemann_states_arr7) :: U, U_star, U_double, F, F_star
        real(wp), dimension(7)    :: F_hlld

        ! Indices for U and F: (rho, rho*vel(1), rho*vel(2), rho*vel(3), By, Bz, E) Note: vel and B are permutated, so vel(1) is the
        ! normal velocity, and x is the normal direction Note: Bx is omitted as the magnetic flux is always zero in the normal
        ! direction

        real(wp) :: v_double, w_double, By_double, Bz_double, E_double
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
                                    & H_no_mag, gamma, pi_inf, qv, vel_rms, B, c, c_fast, pres_mag, U, U_star, U_double, F, &
                                    & F_star, F_hlld, s, s_M, s_starL, s_starR, pTot, p_star, rho_star, E_star, sqrt_rho_star, &
                                    & denom_ds, sign_Bx, v_star, w_star, v_double, w_double, By_double, Bz_double, E_double_lr, &
                                    & E_double]', copyin='[norm_dir]')
                do l = ${Z_BND}$%beg, ${Z_BND}$%end
                    do k = ${Y_BND}$%beg, ${Y_BND}$%end
                        do j = ${X_BND}$%beg, ${X_BND}$%end
                            ! (1) Extract the left/right primitive states
                            do i = 1, eqn_idx%cont%end
                                alpha_rho_L(i) = qL_prim_rsx_vf(${SF('')}$, i)
                                alpha_rho_R(i) = qR_prim_rsx_vf(${SF(' + 1')}$, i)
                            end do

                            ! NOTE: unlike HLL & HLLC, vel%L here is permutated by dir_idx for simpler logic
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
                            s%L = min(vel%L(1) - c_fast%L, vel%R(1) - c_fast%R)
                            s%R = max(vel%R(1) + c_fast%R, vel%L(1) + c_fast%L)

                            pTot%L = pres%L + pres_mag%L
                            pTot%R = pres%R + pres_mag%R

                            s_M = (((s%R - vel%R(1))*rho%R*vel%R(1) - (s%L - vel%L(1))*rho%L*vel%L(1) - pTot%R + pTot%L)/((s%R &
                                   & - vel%R(1))*rho%R - (s%L - vel%L(1))*rho%L))

                            ! (4) Compute star state variables
                            rho_star%L = rho%L*(s%L - vel%L(1))/(s%L - s_M)
                            rho_star%R = rho%R*(s%R - vel%R(1))/(s%R - s_M)
                            p_star = pTot%L + rho%L*(s%L - vel%L(1))*(s_M - vel%L(1))/(s%L - s_M)
                            E_star%L = ((s%L - vel%L(1))*E%L - pTot%L*vel%L(1) + p_star*s_M)/(s%L - s_M)
                            E_star%R = ((s%R - vel%R(1))*E%R - pTot%R*vel%R(1) + p_star*s_M)/(s%R - s_M)

                            ! (5) Compute left/right state vectors and fluxes
                            U%L = [rho%L, rho%L*vel%L(1:3), B%L(2:3), E%L]
                            U_star%L = [rho_star%L, rho_star%L*s_M, rho_star%L*vel%L(2:3), B%L(2:3), E_star%L]
                            U%R = [rho%R, rho%R*vel%R(1:3), B%R(2:3), E%R]
                            U_star%R = [rho_star%R, rho_star%R*s_M, rho_star%R*vel%R(2:3), B%R(2:3), E_star%R]

                            ! Compute the left/right fluxes
                            F%L(1) = U%L(2)
                            F%L(2) = U%L(2)*vel%L(1) - B%L(1)*B%L(1) + pTot%L
                            F%L(3:4) = U%L(2)*vel%L(2:3) - B%L(1)*B%L(2:3)
                            F%L(5:6) = vel%L(1)*B%L(2:3) - vel%L(2:3)*B%L(1)
                            F%L(7) = (E%L + pTot%L)*vel%L(1) - B%L(1)*(vel%L(1)*B%L(1) + vel%L(2)*B%L(2) + vel%L(3)*B%L(3))

                            F%R(1) = U%R(2)
                            F%R(2) = U%R(2)*vel%R(1) - B%R(1)*B%R(1) + pTot%R
                            F%R(3:4) = U%R(2)*vel%R(2:3) - B%R(1)*B%R(2:3)
                            F%R(5:6) = vel%R(1)*B%R(2:3) - vel%R(2:3)*B%R(1)
                            F%R(7) = (E%R + pTot%R)*vel%R(1) - B%R(1)*(vel%R(1)*B%R(1) + vel%R(2)*B%R(2) + vel%R(3)*B%R(3))
                            ! HLLD star-state fluxes via HLL jump relation
                            F_star%L = F%L + s%L*(U_star%L - U%L)
                            F_star%R = F%R + s%R*(U_star%R - U%R)
                            ! Alfven wave speeds bounding the rotational discontinuities
                            s_starL = s_M - abs(B%L(1))/sqrt(rho_star%L)
                            s_starR = s_M + abs(B%L(1))/sqrt(rho_star%R)
                            ! HLLD double-star (intermediate) states across rotational discontinuities
                            sqrt_rho_star%L = sqrt(rho_star%L); sqrt_rho_star%R = sqrt(rho_star%R)
                            v_star%L = vel%L(2); w_star%L = vel%L(3)
                            v_star%R = vel%R(2); w_star%R = vel%R(3)

                            ! (6) Compute the double-star states [Miyoshi Eqns. (59)-(62)]
                            denom_ds = sqrt_rho_star%L + sqrt_rho_star%R
                            sign_Bx = sign(1._wp, B%L(1))
                            v_double = (sqrt_rho_star%L*v_star%L + sqrt_rho_star%R*v_star%R + (B%R(2) - B%L(2))*sign_Bx)/denom_ds
                            w_double = (sqrt_rho_star%L*w_star%L + sqrt_rho_star%R*w_star%R + (B%R(3) - B%L(3))*sign_Bx)/denom_ds
                            By_double = (sqrt_rho_star%L*B%R(2) + sqrt_rho_star%R*B%L(2) &
                                         & + sqrt_rho_star%L*sqrt_rho_star%R*(v_star%R - v_star%L)*sign_Bx)/denom_ds
                            Bz_double = (sqrt_rho_star%L*B%R(3) + sqrt_rho_star%R*B%L(3) &
                                         & + sqrt_rho_star%L*sqrt_rho_star%R*(w_star%R - w_star%L)*sign_Bx)/denom_ds

                            E_double_lr%L = E_star%L - sqrt_rho_star%L*((v_star%L*B%L(2) + w_star%L*B%L(3)) - (v_double*By_double &
                                & + w_double*Bz_double))*sign_Bx
                            E_double_lr%R = E_star%R + sqrt_rho_star%R*((v_star%R*B%R(2) + w_star%R*B%R(3)) - (v_double*By_double &
                                & + w_double*Bz_double))*sign_Bx
                            E_double = 0.5_wp*(E_double_lr%L + E_double_lr%R)

                            U_double%L = [rho_star%L, rho_star%L*s_M, rho_star%L*v_double, rho_star%L*w_double, By_double, &
                                & Bz_double, E_double]
                            U_double%R = [rho_star%R, rho_star%R*s_M, rho_star%R*v_double, rho_star%R*w_double, By_double, &
                                & Bz_double, E_double]

                            ! Select HLLD flux region
                            if (0.0_wp <= s%L) then
                                F_hlld = F%L
                            else if (0.0_wp <= s_starL) then
                                F_hlld = F%L + s%L*(U_star%L - U%L)
                            else if (0.0_wp <= s_M) then
                                F_hlld = F_star%L + s_starL*(U_double%L - U_star%L)
                            else if (0.0_wp <= s_starR) then
                                F_hlld = F_star%R + s_starR*(U_double%R - U_star%R)
                            else if (0.0_wp <= s%R) then
                                F_hlld = F%R + s%R*(U_star%R - U%R)
                            else
                                F_hlld = F%R
                            end if

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

    !> Initialize the Riemann solvers module
    impure subroutine s_initialize_riemann_solvers_module

        ! Allocating the variables that will be utilized to formulate the left, right, and average states of the Riemann problem, as
        ! well the Riemann problem solution
        integer :: i, j

        @:ALLOCATE(Gs_rs(1:num_fluids))

        do i = 1, num_fluids
            Gs_rs(i) = fluid_pp(i)%G
        end do
        $:GPU_UPDATE(device='[Gs_rs]')

        if (viscous) then
            @:ALLOCATE(Res_gs(1:2, 1:Re_size_max))
        end if

        if (viscous) then
            do i = 1, 2
                do j = 1, Re_size(i)
                    Res_gs(i, j) = fluid_pp(Re_idx(i, j))%Re(i)
                end do
            end do
            $:GPU_UPDATE(device='[Res_gs, Re_idx, Re_size]')
        end if

        $:GPU_ENTER_DATA(copyin='[is1, is2, is3, isx, isy, isz]')

        is1%beg = -1; is2%beg = 0; is3%beg = 0
        is1%end = m; is2%end = n; is3%end = p

        @:ALLOCATE(flux_rsx_vf(-1:m, -1:n, -1:p, 1:sys_size))
        @:ALLOCATE(flux_gsrc_rsx_vf(-1:m, -1:n, -1:p, 1:sys_size))
        @:ALLOCATE(flux_src_rsx_vf(-1:m, -1:n, -1:p, eqn_idx%adv%beg:sys_size))
        @:ALLOCATE(vel_src_rsx_vf(-1:m, -1:n, -1:p, 1:num_vels))
        if (qbmm) then
            @:ALLOCATE(mom_sp_rsx_vf(-1:m+1, -1:n+1, -1:p+1, 1:4))
        end if

        if (viscous) then
            @:ALLOCATE(Re_avg_rsx_vf(-1:m, -1:n, -1:p, 1:2))
        end if

    end subroutine s_initialize_riemann_solvers_module

    !> Populate the left and right Riemann state variable buffers based on boundary conditions
    subroutine s_populate_riemann_states_variables_buffers(qL_prim_rsx_vf, dqL_prim_dx_vf, dqL_prim_dy_vf, dqL_prim_dz_vf, &
        & qR_prim_rsx_vf, dqR_prim_dx_vf, dqR_prim_dy_vf, dqR_prim_dz_vf, norm_dir, ix, iy, iz)

        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout) :: qL_prim_rsx_vf, qR_prim_rsx_vf
        type(scalar_field), allocatable, dimension(:), intent(inout) :: dqL_prim_dx_vf, dqR_prim_dx_vf, dqL_prim_dy_vf, &
             & dqR_prim_dy_vf, dqL_prim_dz_vf, dqR_prim_dz_vf
        integer, intent(in)               :: norm_dir
        type(int_bounds_info), intent(in) :: ix, iy, iz
        integer                           :: i, j, k, l  !< Generic loop iterator

        if (norm_dir == 1) then
            is1 = ix; is2 = iy; is3 = iz
            dir_idx = (/1, 2, 3/); dir_flg = (/1._wp, 0._wp, 0._wp/)
        else if (norm_dir == 2) then
            is1 = iy; is2 = ix; is3 = iz
            dir_idx = (/2, 1, 3/); dir_flg = (/0._wp, 1._wp, 0._wp/)
        else
            is1 = iz; is2 = iy; is3 = ix
            dir_idx = (/3, 1, 2/); dir_flg = (/0._wp, 0._wp, 1._wp/)
        end if

        $:GPU_UPDATE(device='[is1, is2, is3]')

        if (elasticity) then
            if (norm_dir == 1) then
                dir_idx_tau = (/1, 2, 4/)
            else if (norm_dir == 2) then
                dir_idx_tau = (/3, 2, 5/)
            else
                dir_idx_tau = (/6, 4, 5/)
            end if
        end if

        isx = ix; isy = iy; isz = iz
        ! for stuff in the same module
        $:GPU_UPDATE(device='[isx, isy, isz]')
        ! for stuff in different modules
        $:GPU_UPDATE(device='[dir_idx, dir_flg, dir_idx_tau]')

        ! Population of Buffers in x-direction
        if (norm_dir == 1) then
            if (bc%x%beg == BC_RIEMANN_EXTRAP) then  ! Riemann state extrap. BC at beginning
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            qL_prim_rsx_vf(-1, k, l, i) = qR_prim_rsx_vf(0, k, l, i)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                if (viscous) then
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = eqn_idx%mom%beg, eqn_idx%mom%end
                        do l = isz%beg, isz%end
                            do k = isy%beg, isy%end
                                dqL_prim_dx_vf(i)%sf(-1, k, l) = dqR_prim_dx_vf(i)%sf(0, k, l)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()

                    if (n > 0) then
                        $:GPU_PARALLEL_LOOP(collapse=3)
                        do i = eqn_idx%mom%beg, eqn_idx%mom%end
                            do l = isz%beg, isz%end
                                do k = isy%beg, isy%end
                                    dqL_prim_dy_vf(i)%sf(-1, k, l) = dqR_prim_dy_vf(i)%sf(0, k, l)
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()

                        if (p > 0) then
                            $:GPU_PARALLEL_LOOP(collapse=3)
                            do i = eqn_idx%mom%beg, eqn_idx%mom%end
                                do l = isz%beg, isz%end
                                    do k = isy%beg, isy%end
                                        dqL_prim_dz_vf(i)%sf(-1, k, l) = dqR_prim_dz_vf(i)%sf(0, k, l)
                                    end do
                                end do
                            end do
                            $:END_GPU_PARALLEL_LOOP()
                        end if
                    end if
                end if
            end if

            if (bc%x%end == BC_RIEMANN_EXTRAP) then  ! Riemann state extrap. BC at end

                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            qR_prim_rsx_vf(m + 1, k, l, i) = qL_prim_rsx_vf(m, k, l, i)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                if (viscous) then
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = eqn_idx%mom%beg, eqn_idx%mom%end
                        do l = isz%beg, isz%end
                            do k = isy%beg, isy%end
                                dqR_prim_dx_vf(i)%sf(m + 1, k, l) = dqL_prim_dx_vf(i)%sf(m, k, l)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()

                    if (n > 0) then
                        $:GPU_PARALLEL_LOOP(collapse=3)
                        do i = eqn_idx%mom%beg, eqn_idx%mom%end
                            do l = isz%beg, isz%end
                                do k = isy%beg, isy%end
                                    dqR_prim_dy_vf(i)%sf(m + 1, k, l) = dqL_prim_dy_vf(i)%sf(m, k, l)
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()

                        if (p > 0) then
                            $:GPU_PARALLEL_LOOP(collapse=3)
                            do i = eqn_idx%mom%beg, eqn_idx%mom%end
                                do l = isz%beg, isz%end
                                    do k = isy%beg, isy%end
                                        dqR_prim_dz_vf(i)%sf(m + 1, k, l) = dqL_prim_dz_vf(i)%sf(m, k, l)
                                    end do
                                end do
                            end do
                            $:END_GPU_PARALLEL_LOOP()
                        end if
                    end if
                end if
            end if
            ! END: Population of Buffers in x-direction

            ! Population of Buffers in y-direction
        else if (norm_dir == 2) then
            if (bc%y%beg == BC_RIEMANN_EXTRAP) then  ! Riemann state extrap. BC at beginning
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            qL_prim_rsx_vf(k, -1, l, i) = qR_prim_rsx_vf(k, 0, l, i)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                if (viscous) then
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = eqn_idx%mom%beg, eqn_idx%mom%end
                        do l = isz%beg, isz%end
                            do j = isx%beg, isx%end
                                dqL_prim_dx_vf(i)%sf(j, -1, l) = dqR_prim_dx_vf(i)%sf(j, 0, l)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = eqn_idx%mom%beg, eqn_idx%mom%end
                        do l = isz%beg, isz%end
                            do j = isx%beg, isx%end
                                dqL_prim_dy_vf(i)%sf(j, -1, l) = dqR_prim_dy_vf(i)%sf(j, 0, l)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()

                    if (p > 0) then
                        $:GPU_PARALLEL_LOOP(collapse=3)
                        do i = eqn_idx%mom%beg, eqn_idx%mom%end
                            do l = isz%beg, isz%end
                                do j = isx%beg, isx%end
                                    dqL_prim_dz_vf(i)%sf(j, -1, l) = dqR_prim_dz_vf(i)%sf(j, 0, l)
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()
                    end if
                end if
            end if

            if (bc%y%end == BC_RIEMANN_EXTRAP) then  ! Riemann state extrap. BC at end

                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            qR_prim_rsx_vf(k, n + 1, l, i) = qL_prim_rsx_vf(k, n, l, i)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                if (viscous) then
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = eqn_idx%mom%beg, eqn_idx%mom%end
                        do l = isz%beg, isz%end
                            do j = isx%beg, isx%end
                                dqR_prim_dx_vf(i)%sf(j, n + 1, l) = dqL_prim_dx_vf(i)%sf(j, n, l)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = eqn_idx%mom%beg, eqn_idx%mom%end
                        do l = isz%beg, isz%end
                            do j = isx%beg, isx%end
                                dqR_prim_dy_vf(i)%sf(j, n + 1, l) = dqL_prim_dy_vf(i)%sf(j, n, l)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()

                    if (p > 0) then
                        $:GPU_PARALLEL_LOOP(collapse=3)
                        do i = eqn_idx%mom%beg, eqn_idx%mom%end
                            do l = isz%beg, isz%end
                                do j = isx%beg, isx%end
                                    dqR_prim_dz_vf(i)%sf(j, n + 1, l) = dqL_prim_dz_vf(i)%sf(j, n, l)
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()
                    end if
                end if
            end if
            ! END: Population of Buffers in y-direction

            ! Population of Buffers in z-direction
        else
            if (bc%z%beg == BC_RIEMANN_EXTRAP) then  ! Riemann state extrap. BC at beginning
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, sys_size
                    do k = is2%beg, is2%end
                        do l = is3%beg, is3%end
                            qL_prim_rsx_vf(l, k, -1, i) = qR_prim_rsx_vf(l, k, 0, i)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                if (viscous) then
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = eqn_idx%mom%beg, eqn_idx%mom%end
                        do k = isy%beg, isy%end
                            do j = isx%beg, isx%end
                                dqL_prim_dx_vf(i)%sf(j, k, -1) = dqR_prim_dx_vf(i)%sf(j, k, 0)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = eqn_idx%mom%beg, eqn_idx%mom%end
                        do k = isy%beg, isy%end
                            do j = isx%beg, isx%end
                                dqL_prim_dy_vf(i)%sf(j, k, -1) = dqR_prim_dy_vf(i)%sf(j, k, 0)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = eqn_idx%mom%beg, eqn_idx%mom%end
                        do k = isy%beg, isy%end
                            do j = isx%beg, isx%end
                                dqL_prim_dz_vf(i)%sf(j, k, -1) = dqR_prim_dz_vf(i)%sf(j, k, 0)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
            end if

            if (bc%z%end == BC_RIEMANN_EXTRAP) then  ! Riemann state extrap. BC at end

                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, sys_size
                    do k = is2%beg, is2%end
                        do l = is3%beg, is3%end
                            qR_prim_rsx_vf(l, k, p + 1, i) = qL_prim_rsx_vf(l, k, p, i)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                if (viscous) then
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = eqn_idx%mom%beg, eqn_idx%mom%end
                        do k = isy%beg, isy%end
                            do j = isx%beg, isx%end
                                dqR_prim_dx_vf(i)%sf(j, k, p + 1) = dqL_prim_dx_vf(i)%sf(j, k, p)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = eqn_idx%mom%beg, eqn_idx%mom%end
                        do k = isy%beg, isy%end
                            do j = isx%beg, isx%end
                                dqR_prim_dy_vf(i)%sf(j, k, p + 1) = dqL_prim_dy_vf(i)%sf(j, k, p)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = eqn_idx%mom%beg, eqn_idx%mom%end
                        do k = isy%beg, isy%end
                            do j = isx%beg, isx%end
                                dqR_prim_dz_vf(i)%sf(j, k, p + 1) = dqL_prim_dz_vf(i)%sf(j, k, p)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
            end if
        end if
        ! END: Population of Buffers in z-direction

    end subroutine s_populate_riemann_states_variables_buffers

    !> Set up the chosen Riemann solver algorithm for the current direction
    subroutine s_initialize_riemann_solver(flux_src_vf, norm_dir)

        type(scalar_field), dimension(sys_size), intent(inout) :: flux_src_vf
        integer, intent(in)                                    :: norm_dir
        integer                                                :: i, j, k, l  !< Generic loop iterators

        ! Reshaping Inputted Data in x-direction

        if (norm_dir == 1) then
            if (viscous .or. (surface_tension)) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = eqn_idx%mom%beg, eqn_idx%E
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end
                                flux_src_vf(i)%sf(j, k, l) = 0._wp
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            if (chem_params%diffusion) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = eqn_idx%E, eqn_idx%species%end
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end
                                if (i == eqn_idx%E .or. i >= eqn_idx%species%beg) then
                                    flux_src_vf(i)%sf(j, k, l) = 0._wp
                                end if
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            if (qbmm) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = 1, 4
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end + 1
                                mom_sp_rsx_vf(j, k, l, i) = mom_sp(i)%sf(j, k, l)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            ! Reshaping Inputted Data in y-direction
        else if (norm_dir == 2) then
            if (viscous .or. (surface_tension)) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = eqn_idx%mom%beg, eqn_idx%E
                    do l = is3%beg, is3%end
                        do j = is1%beg, is1%end
                            do k = is2%beg, is2%end
                                flux_src_vf(i)%sf(k, j, l) = 0._wp
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            if (chem_params%diffusion) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = eqn_idx%E, eqn_idx%species%end
                    do l = is3%beg, is3%end
                        do j = is1%beg, is1%end
                            do k = is2%beg, is2%end
                                if (i == eqn_idx%E .or. i >= eqn_idx%species%beg) then
                                    flux_src_vf(i)%sf(k, j, l) = 0._wp
                                end if
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            if (qbmm) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = 1, 4
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end + 1
                                mom_sp_rsx_vf(k, j, l, i) = mom_sp(i)%sf(k, j, l)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            ! Reshaping Inputted Data in z-direction
        else
            if (viscous .or. (surface_tension)) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = eqn_idx%mom%beg, eqn_idx%E
                    do j = is1%beg, is1%end
                        do k = is2%beg, is2%end
                            do l = is3%beg, is3%end
                                flux_src_vf(i)%sf(l, k, j) = 0._wp
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            if (chem_params%diffusion) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = eqn_idx%E, eqn_idx%species%end
                    do j = is1%beg, is1%end
                        do k = is2%beg, is2%end
                            do l = is3%beg, is3%end
                                if (i == eqn_idx%E .or. i >= eqn_idx%species%beg) then
                                    flux_src_vf(i)%sf(l, k, j) = 0._wp
                                end if
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            if (qbmm) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = 1, 4
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end + 1
                                mom_sp_rsx_vf(l, k, j, i) = mom_sp(i)%sf(l, k, j)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        end if

    end subroutine s_initialize_riemann_solver

    !> Compute cylindrical viscous source flux contributions for momentum and energy
    subroutine s_compute_cylindrical_viscous_source_flux(velL_vf, dvelL_dx_vf, dvelL_dy_vf, dvelL_dz_vf, velR_vf, dvelR_dx_vf, &

        & dvelR_dy_vf, dvelR_dz_vf, flux_src_vf, norm_dir, ix, iy, iz)

        type(scalar_field), dimension(num_dims), intent(in)    :: velL_vf, velR_vf
        type(scalar_field), dimension(num_dims), intent(in)    :: dvelL_dx_vf, dvelR_dx_vf
        type(scalar_field), dimension(num_dims), intent(in)    :: dvelL_dy_vf, dvelR_dy_vf
        type(scalar_field), dimension(num_dims), intent(in)    :: dvelL_dz_vf, dvelR_dz_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: flux_src_vf
        integer, intent(in)                                    :: norm_dir
        type(int_bounds_info), intent(in)                      :: ix, iy, iz

        ! Local variables

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: avg_v_int     !< Averaged interface velocity (\f$v_x, v_y, v_z\f$) (grid directions).
            real(wp), dimension(3) :: avg_dvdx_int  !< Averaged interface \f$\partial v_i/\partial x\f$ (grid dir 1).
            real(wp), dimension(3) :: avg_dvdy_int  !< Averaged interface \f$\partial v_i/\partial y\f$ (grid dir 2).
            real(wp), dimension(3) :: avg_dvdz_int  !< Averaged interface \f$\partial v_i/\partial z\f$ (grid dir 3).
            real(wp), dimension(3) :: vel_src_int   !< Interface velocity (\f$v_1,v_2,v_3\f$) (grid directions) for viscous work.

            !> Shear stress vector (\f$\sigma_{N1}, \sigma_{N2}, \sigma_{N3}\f$) on N-face (grid directions).
            real(wp), dimension(3) :: stress_vector_shear
        #:else
            real(wp), dimension(num_dims) :: avg_v_int     !< Averaged interface velocity (\f$v_x, v_y, v_z\f$) (grid directions).
            real(wp), dimension(num_dims) :: avg_dvdx_int  !< Averaged interface \f$\partial v_i/\partial x\f$ (grid dir 1).
            real(wp), dimension(num_dims) :: avg_dvdy_int  !< Averaged interface \f$\partial v_i/\partial y\f$ (grid dir 2).
            real(wp), dimension(num_dims) :: avg_dvdz_int  !< Averaged interface \f$\partial v_i/\partial z\f$ (grid dir 3).
            !> Interface velocity (\f$v_1,v_2,v_3\f$) (grid directions) for viscous work.
            real(wp), dimension(num_dims) :: vel_src_int
            !> Shear stress vector (\f$\sigma_{N1}, \sigma_{N2}, \sigma_{N3}\f$) on N-face (grid directions).
            real(wp), dimension(num_dims) :: stress_vector_shear
        #:endif
        real(wp) :: stress_normal_bulk  !< Normal bulk stress component \f$\sigma_{NN}\f$ on N-face.
        real(wp) :: Re_s, Re_b  !< Effective interface shear and bulk Reynolds numbers.
        real(wp) :: r_eff  !< Effective radius at interface for cylindrical terms.
        real(wp) :: div_v_term_const  !< Common term \f$-(2/3)(\nabla \cdot \mathbf{v}) / \text{Re}_s\f$ for shear stress diagonal.
        real(wp) :: divergence_cyl  !< Full divergence \f$\nabla \cdot \mathbf{v}\f$ in cylindrical coordinates.
        integer  :: j, k, l  !< Loop iterators for \f$x, y, z\f$ grid directions.
        integer  :: i_vel  !< Loop iterator for velocity components.
        integer  :: idx_rp(3)  !< Indices \f$(j,k,l)\f$ of 'right' point for averaging.

        $:GPU_PARALLEL_LOOP(collapse=3, private='[idx_rp, avg_v_int, avg_dvdx_int, avg_dvdy_int, avg_dvdz_int, Re_s, Re_b, &
                            & vel_src_int, r_eff, divergence_cyl, stress_vector_shear, stress_normal_bulk, div_v_term_const]')
        do l = iz%beg, iz%end
            do k = iy%beg, iy%end
                do j = ix%beg, ix%end
                    ! Determine indices for the 'right' state for averaging across the interface
                    idx_rp = [j, k, l]
                    idx_rp(norm_dir) = idx_rp(norm_dir) + 1

                    ! Average velocities and their derivatives at the interface For cylindrical: x-dir ~ axial (z_cyl), y-dir ~
                    ! radial (r_cyl), z-dir ~ azimuthal (theta_cyl)
                    $:GPU_LOOP(parallelism='[seq]')
                    do i_vel = 1, num_dims
                        avg_v_int(i_vel) = 0.5_wp*(velL_vf(i_vel)%sf(j, k, l) + velR_vf(i_vel)%sf(idx_rp(1), idx_rp(2), idx_rp(3)))

                        avg_dvdx_int(i_vel) = 0.5_wp*(dvelL_dx_vf(i_vel)%sf(j, k, l) + dvelR_dx_vf(i_vel)%sf(idx_rp(1), &
                                     & idx_rp(2), idx_rp(3)))
                        if (num_dims > 1) then
                            avg_dvdy_int(i_vel) = 0.5_wp*(dvelL_dy_vf(i_vel)%sf(j, k, l) + dvelR_dy_vf(i_vel)%sf(idx_rp(1), &
                                         & idx_rp(2), idx_rp(3)))
                        else
                            avg_dvdy_int(i_vel) = 0.0_wp
                        end if
                        if (num_dims > 2) then
                            avg_dvdz_int(i_vel) = 0.5_wp*(dvelL_dz_vf(i_vel)%sf(j, k, l) + dvelR_dz_vf(i_vel)%sf(idx_rp(1), &
                                         & idx_rp(2), idx_rp(3)))
                        else
                            avg_dvdz_int(i_vel) = 0.0_wp
                        end if
                    end do

                    ! Get Re numbers and interface velocity for viscous work
                    select case (norm_dir)
                    case (1)  ! x-face (axial face in z_cyl direction)
                        Re_s = Re_avg_rsx_vf(j, k, l, 1)
                        Re_b = Re_avg_rsx_vf(j, k, l, 2)
                        vel_src_int = vel_src_rsx_vf(j, k, l,1:num_dims)
                        r_eff = y%cc(k)
                    case (2)  ! y-face (radial face in r_cyl direction)
                        Re_s = Re_avg_rsx_vf(j, k, l, 1)
                        Re_b = Re_avg_rsx_vf(j, k, l, 2)
                        vel_src_int = vel_src_rsx_vf(j, k, l,1:num_dims)
                        r_eff = y%cb(k)
                    case (3)  ! z-face (azimuthal face in theta_cyl direction)
                        Re_s = Re_avg_rsx_vf(j, k, l, 1)
                        Re_b = Re_avg_rsx_vf(j, k, l, 2)
                        vel_src_int = vel_src_rsx_vf(j, k, l,1:num_dims)
                        r_eff = y%cc(k)
                    end select

                    ! Divergence in cylindrical coordinates (vx=vz_cyl, vy=vr_cyl, vz=vtheta_cyl)
                    #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
                        divergence_cyl = avg_dvdx_int(1) + avg_dvdy_int(2) + avg_v_int(2)/r_eff
                        if (num_dims > 2) then
                            #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                                divergence_cyl = divergence_cyl + avg_dvdz_int(3)/r_eff
                            #:endif
                        end if
                    #:endif

                    stress_vector_shear = 0.0_wp
                    stress_normal_bulk = 0.0_wp

                    if (shear_stress) then
                        div_v_term_const = -(2.0_wp/3.0_wp)*divergence_cyl/Re_s

                        select case (norm_dir)
                        case (1)  ! X-face (axial normal, z_cyl)
                            stress_vector_shear(1) = (2.0_wp*avg_dvdx_int(1))/Re_s + div_v_term_const
                            if (num_dims > 1) then
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
                                    stress_vector_shear(2) = (avg_dvdy_int(1) + avg_dvdx_int(2))/Re_s
                                #:endif
                            end if
                            if (num_dims > 2) then
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                                    stress_vector_shear(3) = (avg_dvdz_int(1)/r_eff + avg_dvdx_int(3))/Re_s
                                #:endif
                            end if
                        case (2)  ! Y-face (radial normal, r_cyl)
                            if (num_dims > 1) then
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
                                    stress_vector_shear(1) = (avg_dvdy_int(1) + avg_dvdx_int(2))/Re_s
                                    stress_vector_shear(2) = (2.0_wp*avg_dvdy_int(2))/Re_s + div_v_term_const
                                    if (num_dims > 2) then
                                        #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                                            stress_vector_shear(3) = (avg_dvdz_int(2)/r_eff - avg_v_int(3)/r_eff + avg_dvdy_int(3) &
                                                                & )/Re_s
                                        #:endif
                                    end if
                                #:endif
                            else
                                stress_vector_shear(1) = (2.0_wp*avg_dvdx_int(1))/Re_s + div_v_term_const
                            end if
                        case (3)  ! Z-face (azimuthal normal, theta_cyl)
                            if (num_dims > 2) then
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                                    stress_vector_shear(1) = (avg_dvdz_int(1)/r_eff + avg_dvdx_int(3))/Re_s
                                    stress_vector_shear(2) = (avg_dvdz_int(2)/r_eff - avg_v_int(3)/r_eff + avg_dvdy_int(3))/Re_s
                                    stress_vector_shear(3) = (2.0_wp*(avg_dvdz_int(3)/r_eff + avg_v_int(2)/r_eff))/Re_s &
                                                        & + div_v_term_const
                                #:endif
                            end if
                        end select

                        $:GPU_LOOP(parallelism='[seq]')
                        do i_vel = 1, num_dims
                            flux_src_vf(eqn_idx%mom%beg + i_vel - 1)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + i_vel - 1)%sf(j, &
                                        & k, l) - stress_vector_shear(i_vel)
                            flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                        & l) - vel_src_int(i_vel)*stress_vector_shear(i_vel)
                        end do
                    end if

                    if (bulk_stress) then
                        stress_normal_bulk = divergence_cyl/Re_b

                        flux_src_vf(eqn_idx%mom%beg + norm_dir - 1)%sf(j, k, &
                                    & l) = flux_src_vf(eqn_idx%mom%beg + norm_dir - 1)%sf(j, k, l) - stress_normal_bulk
                        flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                    & l) - vel_src_int(norm_dir)*stress_normal_bulk
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_cylindrical_viscous_source_flux

    !> Compute Cartesian viscous source flux contributions for momentum and energy
    subroutine s_compute_cartesian_viscous_source_flux(dvelL_dx_vf, dvelL_dy_vf, dvelL_dz_vf, dvelR_dx_vf, dvelR_dy_vf, &

        & dvelR_dz_vf, flux_src_vf, norm_dir)

        ! Arguments
        type(scalar_field), dimension(num_dims), intent(in)    :: dvelL_dx_vf, dvelR_dx_vf
        type(scalar_field), dimension(num_dims), intent(in)    :: dvelL_dy_vf, dvelR_dy_vf
        type(scalar_field), dimension(num_dims), intent(in)    :: dvelL_dz_vf, dvelR_dz_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: flux_src_vf
        integer, intent(in)                                    :: norm_dir

        ! Local variables

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3, 3) :: vel_grad_avg          !< Averaged velocity gradient tensor `d(vel_i)/d(coord_j)`.
            real(wp), dimension(3, 3) :: current_tau_shear     !< Current shear stress tensor.
            real(wp), dimension(3, 3) :: current_tau_bulk      !< Current bulk stress tensor.
            real(wp), dimension(3)    :: vel_src_at_interface  !< Interface velocities (u,v,w) for viscous work.
        #:else
            real(wp), dimension(num_dims, num_dims) :: vel_grad_avg  !< Averaged velocity gradient tensor `d(vel_i)/d(coord_j)`.
            real(wp), dimension(num_dims, num_dims) :: current_tau_shear  !< Current shear stress tensor.
            real(wp), dimension(num_dims, num_dims) :: current_tau_bulk  !< Current bulk stress tensor.
            real(wp), dimension(num_dims)           :: vel_src_at_interface  !< Interface velocities (u,v,w) for viscous work.
        #:endif
        integer, dimension(3) :: idx_right_phys  !< Physical (j,k,l) indices for right state.
        real(wp)              :: Re_shear        !< Interface shear Reynolds number.
        real(wp)              :: Re_bulk         !< Interface bulk Reynolds number.
        integer               :: j_loop          !< Physical x-index loop iterator.
        integer               :: k_loop          !< Physical y-index loop iterator.
        integer               :: l_loop          !< Physical z-index loop iterator.
        integer               :: i_dim           !< Generic dimension/component iterator.
        integer               :: vel_comp_idx    !< Velocity component iterator (1=u, 2=v, 3=w).
        real(wp)              :: divergence_v    !< Velocity divergence at interface.

        $:GPU_PARALLEL_LOOP(collapse=3, private='[idx_right_phys, vel_grad_avg, current_tau_shear, current_tau_bulk, &
                            & vel_src_at_interface, Re_shear, Re_bulk, divergence_v, i_dim, vel_comp_idx]')
        do l_loop = isz%beg, isz%end
            do k_loop = isy%beg, isy%end
                do j_loop = isx%beg, isx%end
                    idx_right_phys(1) = j_loop
                    idx_right_phys(2) = k_loop
                    idx_right_phys(3) = l_loop
                    idx_right_phys(norm_dir) = idx_right_phys(norm_dir) + 1

                    vel_grad_avg = 0.0_wp
                    do vel_comp_idx = 1, num_dims
                        vel_grad_avg(vel_comp_idx, 1) = 0.5_wp*(dvelL_dx_vf(vel_comp_idx)%sf(j_loop, k_loop, &
                                     & l_loop) + dvelR_dx_vf(vel_comp_idx)%sf(idx_right_phys(1), idx_right_phys(2), &
                                     & idx_right_phys(3)))
                        if (num_dims > 1) then
                            #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
                                vel_grad_avg(vel_comp_idx, 2) = 0.5_wp*(dvelL_dy_vf(vel_comp_idx)%sf(j_loop, k_loop, &
                                             & l_loop) + dvelR_dy_vf(vel_comp_idx)%sf(idx_right_phys(1), idx_right_phys(2), &
                                             & idx_right_phys(3)))
                            #:endif
                        end if
                        if (num_dims > 2) then
                            #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                                vel_grad_avg(vel_comp_idx, 3) = 0.5_wp*(dvelL_dz_vf(vel_comp_idx)%sf(j_loop, k_loop, &
                                             & l_loop) + dvelR_dz_vf(vel_comp_idx)%sf(idx_right_phys(1), idx_right_phys(2), &
                                             & idx_right_phys(3)))
                            #:endif
                        end if
                    end do

                    divergence_v = 0.0_wp
                    do i_dim = 1, num_dims
                        divergence_v = divergence_v + vel_grad_avg(i_dim, i_dim)
                    end do

                    vel_src_at_interface = 0.0_wp
                    if (norm_dir == 1) then
                        Re_shear = Re_avg_rsx_vf(j_loop, k_loop, l_loop, 1)
                        Re_bulk = Re_avg_rsx_vf(j_loop, k_loop, l_loop, 2)
                        do i_dim = 1, num_dims
                            vel_src_at_interface(i_dim) = vel_src_rsx_vf(j_loop, k_loop, l_loop, i_dim)
                        end do
                    else if (norm_dir == 2) then
                        Re_shear = Re_avg_rsx_vf(j_loop, k_loop, l_loop, 1)
                        Re_bulk = Re_avg_rsx_vf(j_loop, k_loop, l_loop, 2)
                        do i_dim = 1, num_dims
                            vel_src_at_interface(i_dim) = vel_src_rsx_vf(j_loop, k_loop, l_loop, i_dim)
                        end do
                    else
                        Re_shear = Re_avg_rsx_vf(j_loop, k_loop, l_loop, 1)
                        Re_bulk = Re_avg_rsx_vf(j_loop, k_loop, l_loop, 2)
                        do i_dim = 1, num_dims
                            vel_src_at_interface(i_dim) = vel_src_rsx_vf(j_loop, k_loop, l_loop, i_dim)
                        end do
                    end if

                    if (shear_stress) then
                        ! current_tau_shear = 0.0_wp
                        call s_calculate_shear_stress_tensor(vel_grad_avg, Re_shear, divergence_v, current_tau_shear)

                        do i_dim = 1, num_dims
                            flux_src_vf(eqn_idx%mom%beg + i_dim - 1)%sf(j_loop, k_loop, &
                                        & l_loop) = flux_src_vf(eqn_idx%mom%beg + i_dim - 1)%sf(j_loop, k_loop, &
                                        & l_loop) - current_tau_shear(norm_dir, i_dim)

                            flux_src_vf(eqn_idx%E)%sf(j_loop, k_loop, l_loop) = flux_src_vf(eqn_idx%E)%sf(j_loop, k_loop, &
                                        & l_loop) - vel_src_at_interface(i_dim)*current_tau_shear(norm_dir, i_dim)
                        end do
                    end if

                    if (bulk_stress) then
                        ! current_tau_bulk = 0.0_wp
                        call s_calculate_bulk_stress_tensor(Re_bulk, divergence_v, current_tau_bulk)

                        do i_dim = 1, num_dims
                            flux_src_vf(eqn_idx%mom%beg + i_dim - 1)%sf(j_loop, k_loop, &
                                        & l_loop) = flux_src_vf(eqn_idx%mom%beg + i_dim - 1)%sf(j_loop, k_loop, &
                                        & l_loop) - current_tau_bulk(norm_dir, i_dim)

                            flux_src_vf(eqn_idx%E)%sf(j_loop, k_loop, l_loop) = flux_src_vf(eqn_idx%E)%sf(j_loop, k_loop, &
                                        & l_loop) - vel_src_at_interface(i_dim)*current_tau_bulk(norm_dir, i_dim)
                        end do
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_cartesian_viscous_source_flux

    !> Compute shear stress tensor components
    subroutine s_calculate_shear_stress_tensor(vel_grad_avg, Re_shear, divergence_v, tau_shear_out)

        $:GPU_ROUTINE(parallelism='[seq]')

        ! Arguments
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3, 3), intent(in)  :: vel_grad_avg
            real(wp), dimension(3, 3), intent(out) :: tau_shear_out
        #:else
            real(wp), dimension(num_dims, num_dims), intent(in)  :: vel_grad_avg
            real(wp), dimension(num_dims, num_dims), intent(out) :: tau_shear_out
        #:endif
        real(wp), intent(in) :: Re_shear
        real(wp), intent(in) :: divergence_v

        ! Local variables
        integer :: i_dim  !< Loop iterator for face normal.
        integer :: j_dim  !< Loop iterator for force component direction.
        tau_shear_out = 0.0_wp

        do i_dim = 1, num_dims
            do j_dim = 1, num_dims
                tau_shear_out(i_dim, j_dim) = (vel_grad_avg(j_dim, i_dim) + vel_grad_avg(i_dim, j_dim))/Re_shear
                if (i_dim == j_dim) then
                    tau_shear_out(i_dim, j_dim) = tau_shear_out(i_dim, j_dim) - (2.0_wp/3.0_wp)*divergence_v/Re_shear
                end if
            end do
        end do

    end subroutine s_calculate_shear_stress_tensor

    !> Compute bulk stress tensor components (diagonal only)
    subroutine s_calculate_bulk_stress_tensor(Re_bulk, divergence_v, tau_bulk_out)

        $:GPU_ROUTINE(parallelism='[seq]')

        ! Arguments
        real(wp), intent(in) :: Re_bulk
        real(wp), intent(in) :: divergence_v
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3, 3), intent(out) :: tau_bulk_out
        #:else
            real(wp), dimension(num_dims, num_dims), intent(out) :: tau_bulk_out
        #:endif

        ! Local variables
        integer :: i_dim  !< Loop iterator for diagonal components.
        tau_bulk_out = 0.0_wp

        do i_dim = 1, num_dims
            tau_bulk_out(i_dim, i_dim) = divergence_v/Re_bulk
        end do

    end subroutine s_calculate_bulk_stress_tensor

    !> Deallocation and/or disassociation procedures that are needed to finalize the selected Riemann problem solver
    subroutine s_finalize_riemann_solver(flux_vf, flux_src_vf, flux_gsrc_vf, norm_dir)

        type(scalar_field), dimension(sys_size), intent(inout) :: flux_vf, flux_src_vf, flux_gsrc_vf
        integer, intent(in)                                    :: norm_dir
        integer                                                :: i, j, k, l  !< Generic loop iterators
        ! Reshaping Outputted Data in y-direction

        if (norm_dir == 2) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do l = is3%beg, is3%end
                    do j = is1%beg, is1%end
                        do k = is2%beg, is2%end
                            flux_vf(i)%sf(k, j, l) = flux_rsx_vf(k, j, l, i)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            if (cyl_coord) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do j = is1%beg, is1%end
                            do k = is2%beg, is2%end
                                flux_gsrc_vf(i)%sf(k, j, l) = flux_gsrc_rsx_vf(k, j, l, i)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            $:GPU_PARALLEL_LOOP(collapse=3)
            do l = is3%beg, is3%end
                do j = is1%beg, is1%end
                    do k = is2%beg, is2%end
                        flux_src_vf(eqn_idx%adv%beg)%sf(k, j, l) = flux_src_rsx_vf(k, j, l, eqn_idx%adv%beg)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            if (riemann_solver == 1 .or. riemann_solver == 4) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = eqn_idx%adv%beg + 1, eqn_idx%adv%end
                    do l = is3%beg, is3%end
                        do j = is1%beg, is1%end
                            do k = is2%beg, is2%end
                                flux_src_vf(i)%sf(k, j, l) = flux_src_rsx_vf(k, j, l, i)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
            ! Reshaping Outputted Data in z-direction
        else if (norm_dir == 3) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do j = is1%beg, is1%end
                    do k = is2%beg, is2%end
                        do l = is3%beg, is3%end
                            flux_vf(i)%sf(l, k, j) = flux_rsx_vf(l, k, j, i)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            if (grid_geometry == 3) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = 1, sys_size
                    do j = is1%beg, is1%end
                        do k = is2%beg, is2%end
                            do l = is3%beg, is3%end
                                flux_gsrc_vf(i)%sf(l, k, j) = flux_gsrc_rsx_vf(l, k, j, i)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            $:GPU_PARALLEL_LOOP(collapse=3)
            do j = is1%beg, is1%end
                do k = is2%beg, is2%end
                    do l = is3%beg, is3%end
                        flux_src_vf(eqn_idx%adv%beg)%sf(l, k, j) = flux_src_rsx_vf(l, k, j, eqn_idx%adv%beg)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            if (riemann_solver == 1 .or. riemann_solver == 4) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = eqn_idx%adv%beg + 1, eqn_idx%adv%end
                    do j = is1%beg, is1%end
                        do k = is2%beg, is2%end
                            do l = is3%beg, is3%end
                                flux_src_vf(i)%sf(l, k, j) = flux_src_rsx_vf(l, k, j, i)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        else if (norm_dir == 1) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do l = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = is1%beg, is1%end
                            flux_vf(i)%sf(j, k, l) = flux_rsx_vf(j, k, l, i)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            $:GPU_PARALLEL_LOOP(collapse=3)
            do l = is3%beg, is3%end
                do k = is2%beg, is2%end
                    do j = is1%beg, is1%end
                        flux_src_vf(eqn_idx%adv%beg)%sf(j, k, l) = flux_src_rsx_vf(j, k, l, eqn_idx%adv%beg)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            if (riemann_solver == 1 .or. riemann_solver == 4) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = eqn_idx%adv%beg + 1, eqn_idx%adv%end
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end
                                flux_src_vf(i)%sf(j, k, l) = flux_src_rsx_vf(j, k, l, i)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        end if

    end subroutine s_finalize_riemann_solver

    !> Module deallocation and/or disassociation procedures
    impure subroutine s_finalize_riemann_solvers_module

        if (viscous) then
            @:DEALLOCATE(Re_avg_rsx_vf)
        end if
        @:DEALLOCATE(vel_src_rsx_vf)
        @:DEALLOCATE(flux_rsx_vf)
        @:DEALLOCATE(flux_src_rsx_vf)
        @:DEALLOCATE(flux_gsrc_rsx_vf)
        if (qbmm) then
            @:DEALLOCATE(mom_sp_rsx_vf)
        end if

    end subroutine s_finalize_riemann_solvers_module

end module m_riemann_solvers
