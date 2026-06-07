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

    #:if USING_AMD
        use m_chemistry, only: molecular_weights_nonparameter
    #:endif

    implicit none

    private; public :: s_initialize_riemann_solvers_module, s_riemann_solver, s_hll_riemann_solver, s_hllc_riemann_solver, &
        & s_hlld_riemann_solver, s_hypo_hlld_riemann_solver, s_finalize_nc_iface_vel, s_lf_riemann_solver, &
        & s_finalize_riemann_solver_hatR, s_finalize_nc_iface_vel_hatR, s_finalize_riemann_solvers_module

    !> The cell-boundary values of the fluxes (src - source) that are computed through the chosen Riemann problem solver, and the
    !! direct evaluation of source terms, by using the left and right states given in qK_prim_rs_vf, dqK_prim_ds_vf where ds = dx,
    !! dy or dz.
    !> @{
    real(wp), allocatable, dimension(:,:,:,:) :: flux_rsx_vf, flux_src_rsx_vf
    real(wp), allocatable, dimension(:,:,:,:) :: flux_rsy_vf, flux_src_rsy_vf
    real(wp), allocatable, dimension(:,:,:,:) :: flux_rsz_vf, flux_src_rsz_vf
    $:GPU_DECLARE(create='[flux_rsx_vf, flux_src_rsx_vf, flux_rsy_vf, flux_src_rsy_vf, flux_rsz_vf, flux_src_rsz_vf]')
    !> @}

    !> The cell-boundary values of the geometrical source flux that are computed through the chosen Riemann problem solver by using
    !! the left and right states given in qK_prim_rs_vf. Currently 2D axisymmetric for inviscid only.
    !> @{
    real(wp), allocatable, dimension(:,:,:,:) :: flux_gsrc_rsx_vf
    real(wp), allocatable, dimension(:,:,:,:) :: flux_gsrc_rsy_vf
    real(wp), allocatable, dimension(:,:,:,:) :: flux_gsrc_rsz_vf
    $:GPU_DECLARE(create='[flux_gsrc_rsx_vf, flux_gsrc_rsy_vf, flux_gsrc_rsz_vf]')

    real(wp), allocatable, dimension(:,:,:,:) :: nc_iface_vel_rsx_vf
    real(wp), allocatable, dimension(:,:,:,:) :: nc_iface_vel_rsy_vf
    real(wp), allocatable, dimension(:,:,:,:) :: nc_iface_vel_rsz_vf
    $:GPU_DECLARE(create='[nc_iface_vel_rsx_vf, nc_iface_vel_rsy_vf, nc_iface_vel_rsz_vf]')
    !> @}

    !> Dual-pass HLLD second flux set: the hat_R-anchored fluxes (and, for axisymmetric runs, the hat_R interface velocities)
    !! written by the same fused solve that fills flux_rs* / nc_iface_vel_rs* with the hat_L-anchored values. Allocated only when
    !! hypo_nc_dual_pass.
    !> @{
    real(wp), allocatable, dimension(:,:,:,:) :: flux_hatR_rsx_vf
    real(wp), allocatable, dimension(:,:,:,:) :: flux_hatR_rsy_vf
    real(wp), allocatable, dimension(:,:,:,:) :: flux_hatR_rsz_vf
    $:GPU_DECLARE(create='[flux_hatR_rsx_vf, flux_hatR_rsy_vf, flux_hatR_rsz_vf]')

    real(wp), allocatable, dimension(:,:,:,:) :: nc_iface_vel_hatR_rsx_vf
    real(wp), allocatable, dimension(:,:,:,:) :: nc_iface_vel_hatR_rsy_vf
    real(wp), allocatable, dimension(:,:,:,:) :: nc_iface_vel_hatR_rsz_vf
    $:GPU_DECLARE(create='[nc_iface_vel_hatR_rsx_vf, nc_iface_vel_hatR_rsy_vf, nc_iface_vel_hatR_rsz_vf]')

    real(wp), allocatable, dimension(:,:,:,:) :: flux_gsrc_hatR_rsx_vf
    real(wp), allocatable, dimension(:,:,:,:) :: flux_gsrc_hatR_rsy_vf
    real(wp), allocatable, dimension(:,:,:,:) :: flux_gsrc_hatR_rsz_vf
    $:GPU_DECLARE(create='[flux_gsrc_hatR_rsx_vf, flux_gsrc_hatR_rsy_vf, flux_gsrc_hatR_rsz_vf]')
    !> @}

    ! Cell-boundary velocity from Riemann solution; used for source flux

    real(wp), allocatable, dimension(:,:,:,:) :: vel_src_rsx_vf
    real(wp), allocatable, dimension(:,:,:,:) :: vel_src_rsy_vf
    real(wp), allocatable, dimension(:,:,:,:) :: vel_src_rsz_vf
    $:GPU_DECLARE(create='[vel_src_rsx_vf, vel_src_rsy_vf, vel_src_rsz_vf]')

    real(wp), allocatable, dimension(:,:,:,:) :: mom_sp_rsx_vf
    real(wp), allocatable, dimension(:,:,:,:) :: mom_sp_rsy_vf
    real(wp), allocatable, dimension(:,:,:,:) :: mom_sp_rsz_vf
    $:GPU_DECLARE(create='[mom_sp_rsx_vf, mom_sp_rsy_vf, mom_sp_rsz_vf]')

    real(wp), allocatable, dimension(:,:,:,:) :: Re_avg_rsx_vf
    real(wp), allocatable, dimension(:,:,:,:) :: Re_avg_rsy_vf
    real(wp), allocatable, dimension(:,:,:,:) :: Re_avg_rsz_vf
    $:GPU_DECLARE(create='[Re_avg_rsx_vf, Re_avg_rsy_vf, Re_avg_rsz_vf]')

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
    subroutine s_riemann_solver(qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, dqL_prim_dx_vf, dqL_prim_dy_vf, dqL_prim_dz_vf, &

        & qL_prim_vf, qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, dqR_prim_dx_vf, dqR_prim_dy_vf, dqR_prim_dz_vf, qR_prim_vf, &
            & q_prim_vf, flux_vf, flux_src_vf, flux_gsrc_vf, norm_dir, ix, iy, iz)

        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout) :: qL_prim_rsx_vf, qL_prim_rsy_vf, &
             & qL_prim_rsz_vf, qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf
        type(scalar_field), dimension(sys_size), intent(in)          :: q_prim_vf
        type(scalar_field), allocatable, dimension(:), intent(inout) :: qL_prim_vf, qR_prim_vf
        type(scalar_field), allocatable, dimension(:), intent(inout) :: dqL_prim_dx_vf, dqR_prim_dx_vf, dqL_prim_dy_vf, &
             & dqR_prim_dy_vf, dqL_prim_dz_vf, dqR_prim_dz_vf

        type(scalar_field), dimension(sys_size), intent(inout) :: flux_vf, flux_src_vf, flux_gsrc_vf
        integer, intent(in)                                    :: norm_dir
        type(int_bounds_info), intent(in)                      :: ix, iy, iz

        if (hypo_nc_dual_pass) then
            ! Fused dual-pass: one call computes BOTH anchored flux sets (hat_L -> flux_vf via the regular finalize; hat_R into
            ! flux_hatR_rs*, finalized separately via s_finalize_riemann_solver_hatR between the two RHS assemblies).
            call s_hypo_hlld_riemann_solver(qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, dqL_prim_dx_vf, dqL_prim_dy_vf, &
                                            & dqL_prim_dz_vf, qL_prim_vf, qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, &
                                            & dqR_prim_dx_vf, dqR_prim_dy_vf, dqR_prim_dz_vf, qR_prim_vf, q_prim_vf, flux_vf, &
                                            & flux_src_vf, flux_gsrc_vf, norm_dir, ix, iy, iz)

            return
        end if

        #:for NAME, NUM in [('hll', 1), ('hllc', 2), ('hlld', 4), ('lf', 5)]
            if (riemann_solver == ${NUM}$) then
                call s_${NAME}$_riemann_solver(qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, dqL_prim_dx_vf, dqL_prim_dy_vf, &
                                               & dqL_prim_dz_vf, qL_prim_vf, qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, &
                                               & dqR_prim_dx_vf, dqR_prim_dy_vf, dqR_prim_dz_vf, qR_prim_vf, q_prim_vf, flux_vf, &
                                               & flux_src_vf, flux_gsrc_vf, norm_dir, ix, iy, iz)
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
    subroutine s_hll_riemann_solver(qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, dqL_prim_dx_vf, dqL_prim_dy_vf, &

        & dqL_prim_dz_vf, qL_prim_vf, qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, dqR_prim_dx_vf, dqR_prim_dy_vf, &
            & dqR_prim_dz_vf, qR_prim_vf, q_prim_vf, flux_vf, flux_src_vf, flux_gsrc_vf, norm_dir, ix, iy, iz)

        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout) :: qL_prim_rsx_vf, qL_prim_rsy_vf, &
             & qL_prim_rsz_vf, qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf
        type(scalar_field), dimension(sys_size), intent(in)          :: q_prim_vf
        type(scalar_field), allocatable, dimension(:), intent(inout) :: qL_prim_vf, qR_prim_vf
        type(scalar_field), allocatable, dimension(:), intent(inout) :: dqL_prim_dx_vf, dqR_prim_dx_vf, dqL_prim_dy_vf, &
             & dqR_prim_dy_vf, dqL_prim_dz_vf, dqR_prim_dz_vf

        ! Intercell fluxes
        type(scalar_field), dimension(sys_size), intent(inout) :: flux_vf, flux_src_vf, flux_gsrc_vf
        real(wp)                                               :: flux_tau_L, flux_tau_R
        integer, intent(in)                                    :: norm_dir
        type(int_bounds_info), intent(in)                      :: ix, iy, iz

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3)  :: alpha_rho_L, alpha_rho_R
            real(wp), dimension(3)  :: vel_L, vel_R
            real(wp), dimension(3)  :: alpha_L, alpha_R
            real(wp), dimension(10) :: Ys_L, Ys_R
            real(wp), dimension(10) :: Cp_iL, Cp_iR, Xs_L, Xs_R, Gamma_iL, Gamma_iR
            real(wp), dimension(10) :: Yi_avg, Phi_avg, h_iL, h_iR, h_avg_2
        #:else
            real(wp), dimension(num_fluids)  :: alpha_rho_L, alpha_rho_R
            real(wp), dimension(num_vels)    :: vel_L, vel_R
            real(wp), dimension(num_fluids)  :: alpha_L, alpha_R
            real(wp), dimension(num_species) :: Ys_L, Ys_R
            real(wp), dimension(num_species) :: Cp_iL, Cp_iR, Xs_L, Xs_R, Gamma_iL, Gamma_iR
            real(wp), dimension(num_species) :: Yi_avg, Phi_avg, h_iL, h_iR, h_avg_2
        #:endif
        real(wp)                  :: rho_L, rho_R
        real(wp)                  :: pres_L, pres_R
        real(wp)                  :: E_L, E_R
        real(wp)                  :: H_L, H_R
        real(wp)                  :: Cp_avg, Cv_avg, T_avg, eps, c_sum_Yi_Phi
        real(wp)                  :: T_L, T_R
        real(wp)                  :: Y_L, Y_R
        real(wp)                  :: MW_L, MW_R
        real(wp)                  :: R_gas_L, R_gas_R
        real(wp)                  :: Cp_L, Cp_R
        real(wp)                  :: Cv_L, Cv_R
        real(wp)                  :: Gamm_L, Gamm_R
        real(wp)                  :: gamma_L, gamma_R
        real(wp)                  :: pi_inf_L, pi_inf_R
        real(wp)                  :: qv_L, qv_R
        real(wp)                  :: c_L, c_R
        real(wp), dimension(6)    :: tau_e_L, tau_e_R
        real(wp)                  :: G_L, G_R
        real(wp), dimension(2)    :: Re_L, Re_R
        real(wp), dimension(3)    :: xi_field_L, xi_field_R
        real(wp)                  :: rho_avg
        real(wp)                  :: H_avg
        real(wp)                  :: qv_avg
        real(wp)                  :: gamma_avg
        real(wp)                  :: c_avg
        real(wp)                  :: s_L, s_R, s_M, s_P, s_S
        real(wp)                  :: xi_M, xi_P
        real(wp)                  :: ptilde_L, ptilde_R
        real(wp)                  :: vel_L_rms, vel_R_rms, vel_avg_rms
        real(wp)                  :: vel_L_tmp, vel_R_tmp
        real(wp)                  :: Ms_L, Ms_R, pres_SL, pres_SR
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

        call s_populate_riemann_states_variables_buffers(qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, dqL_prim_dx_vf, &
            & dqL_prim_dy_vf, dqL_prim_dz_vf, qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, dqR_prim_dx_vf, dqR_prim_dy_vf, &
            & dqR_prim_dz_vf, norm_dir, ix, iy, iz)

        ! Reshaping inputted data based on dimensional splitting direction
        call s_initialize_riemann_solver(flux_src_vf, norm_dir)
        #:for NORM_DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
            if (norm_dir == ${NORM_DIR}$) then
                $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, q, alpha_rho_L, alpha_rho_R, vel_L, vel_R, alpha_L, &
                                    & alpha_R, tau_e_L, tau_e_R, Re_L, Re_R, s_L, s_R, s_S, Ys_L, Ys_R, xi_field_L, xi_field_R, &
                                    & Cp_iL, Cp_iR, Xs_L, Xs_R, Gamma_iL, Gamma_iR, Yi_avg, Phi_avg, h_iL, h_iR, h_avg_2, c_fast, &
                                    & pres_mag, B, Ga, vdotB, B2, b4, cm, pcorr, zcoef, vel_L_tmp, vel_R_tmp, rho_L, rho_R, &
                                    & pres_L, pres_R, E_L, E_R, H_L, H_R, Cp_avg, Cv_avg, T_avg, eps, c_sum_Yi_Phi, T_L, T_R, &
                                    & Y_L, Y_R, MW_L, MW_R, R_gas_L, R_gas_R, Cp_L, Cp_R, Cv_L, Cv_R, Gamm_L, Gamm_R, gamma_L, &
                                    & gamma_R, pi_inf_L, pi_inf_R, qv_L, qv_R, qv_avg, c_L, c_R, G_L, G_R, rho_avg, H_avg, c_avg, &
                                    & gamma_avg, ptilde_L, ptilde_R, vel_L_rms, vel_R_rms, vel_avg_rms, Ms_L, Ms_R, pres_SL, &
                                    & pres_SR, alpha_L_sum, alpha_R_sum, flux_tau_L, flux_tau_R, s_M, s_P, xi_M, &
                                        & xi_P]', copyin='[norm_dir]')
                do l = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = is1%beg, is1%end
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, eqn_idx%cont%end
                                alpha_rho_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                alpha_rho_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)
                            end do

                            vel_L_rms = 0._wp; vel_R_rms = 0._wp

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_vels
                                vel_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%cont%end + i)
                                vel_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%cont%end + i)
                                vel_L_rms = vel_L_rms + vel_L(i)**2._wp
                                vel_R_rms = vel_R_rms + vel_R(i)**2._wp
                            end do

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                alpha_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%E + i)
                                alpha_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%E + i)
                            end do

                            pres_L = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%E)
                            pres_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%E)

                            if (mhd) then
                                if (n == 0) then  ! 1D: constant Bx; By, Bz as variables
                                    B%L(1) = Bx0
                                    B%R(1) = Bx0
                                    B%L(2) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%B%beg)
                                    B%R(2) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%B%beg)
                                    B%L(3) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%B%beg + 1)
                                    B%R(3) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%B%beg + 1)
                                else  ! 2D/3D: Bx, By, Bz as variables
                                    B%L(1) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%B%beg)
                                    B%R(1) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%B%beg)
                                    B%L(2) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%B%beg + 1)
                                    B%R(2) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%B%beg + 1)
                                    B%L(3) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%B%beg + 2)
                                    B%R(3) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%B%beg + 2)
                                end if
                            end if

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
                                rho_L = rho_L + alpha_rho_L(i)
                                gamma_L = gamma_L + alpha_L(i)*gammas(i)
                                pi_inf_L = pi_inf_L + alpha_L(i)*pi_infs(i)
                                qv_L = qv_L + alpha_rho_L(i)*qvs(i)

                                rho_R = rho_R + alpha_rho_R(i)
                                gamma_R = gamma_R + alpha_R(i)*gammas(i)
                                pi_inf_R = pi_inf_R + alpha_R(i)*pi_infs(i)
                                qv_R = qv_R + alpha_rho_R(i)*qvs(i)
                            end do

                            if (viscous) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, 2
                                    Re_L(i) = dflt_real
                                    Re_R(i) = dflt_real

                                    if (Re_size(i) > 0) Re_L(i) = 0._wp
                                    if (Re_size(i) > 0) Re_R(i) = 0._wp

                                    $:GPU_LOOP(parallelism='[seq]')
                                    do q = 1, Re_size(i)
                                        Re_L(i) = alpha_L(Re_idx(i, q))/Res_gs(i, q) + Re_L(i)
                                        Re_R(i) = alpha_R(Re_idx(i, q))/Res_gs(i, q) + Re_R(i)
                                    end do

                                    Re_L(i) = 1._wp/max(Re_L(i), sgm_eps)
                                    Re_R(i) = 1._wp/max(Re_R(i), sgm_eps)
                                end do
                            end if

                            if (chemistry) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = eqn_idx%species%beg, eqn_idx%species%end
                                    Ys_L(i - eqn_idx%species%beg + 1) = qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                    Ys_R(i - eqn_idx%species%beg + 1) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)
                                end do

                                call get_mixture_molecular_weight(Ys_L, MW_L)
                                call get_mixture_molecular_weight(Ys_R, MW_R)
                                #:if USING_AMD
                                    Xs_L(:) = Ys_L(:)*MW_L/molecular_weights_nonparameter(:)
                                    Xs_R(:) = Ys_R(:)*MW_R/molecular_weights_nonparameter(:)
                                #:else
                                    Xs_L(:) = Ys_L(:)*MW_L/molecular_weights(:)
                                    Xs_R(:) = Ys_R(:)*MW_R/molecular_weights(:)
                                #:endif

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
                                H_L = (E_L + pres_L)/rho_L
                                H_R = (E_R + pres_R)/rho_R
                            else if (mhd .and. relativity) then
                                Ga%L = 1._wp/sqrt(1._wp - vel_L_rms)
                                Ga%R = 1._wp/sqrt(1._wp - vel_R_rms)
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                                    vdotB%L = vel_L(1)*B%L(1) + vel_L(2)*B%L(2) + vel_L(3)*B%L(3)
                                    vdotB%R = vel_R(1)*B%R(1) + vel_R(2)*B%R(2) + vel_R(3)*B%R(3)

                                    b4%L(1:3) = B%L(1:3)/Ga%L + Ga%L*vel_L(1:3)*vdotB%L
                                    b4%R(1:3) = B%R(1:3)/Ga%R + Ga%R*vel_R(1:3)*vdotB%R
                                    B2%L = B%L(1)**2._wp + B%L(2)**2._wp + B%L(3)**2._wp
                                    B2%R = B%R(1)**2._wp + B%R(2)**2._wp + B%R(3)**2._wp
                                #:endif

                                pres_mag%L = 0.5_wp*(B2%L/Ga%L**2._wp + vdotB%L**2._wp)
                                pres_mag%R = 0.5_wp*(B2%R/Ga%R**2._wp + vdotB%R**2._wp)

                                ! Hard-coded EOS
                                H_L = 1._wp + (gamma_L + 1)*pres_L/rho_L
                                H_R = 1._wp + (gamma_R + 1)*pres_R/rho_R
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                                    cm%L(1:3) = (rho_L*H_L*Ga%L**2 + B2%L)*vel_L(1:3) - vdotB%L*B%L(1:3)
                                    cm%R(1:3) = (rho_R*H_R*Ga%R**2 + B2%R)*vel_R(1:3) - vdotB%R*B%R(1:3)
                                #:endif

                                E_L = rho_L*H_L*Ga%L**2 - pres_L + 0.5_wp*(B2%L + vel_L_rms*B2%L - vdotB%L**2._wp) - rho_L*Ga%L
                                E_R = rho_R*H_R*Ga%R**2 - pres_R + 0.5_wp*(B2%R + vel_R_rms*B2%R - vdotB%R**2._wp) - rho_R*Ga%R
                            else if (mhd .and. .not. relativity) then
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                                    pres_mag%L = 0.5_wp*(B%L(1)**2._wp + B%L(2)**2._wp + B%L(3)**2._wp)
                                    pres_mag%R = 0.5_wp*(B%R(1)**2._wp + B%R(2)**2._wp + B%R(3)**2._wp)
                                #:endif
                                E_L = gamma_L*pres_L + pi_inf_L + 0.5_wp*rho_L*vel_L_rms + qv_L + pres_mag%L
                                ! includes magnetic energy
                                E_R = gamma_R*pres_R + pi_inf_R + 0.5_wp*rho_R*vel_R_rms + qv_R + pres_mag%R
                                H_L = (E_L + pres_L - pres_mag%L)/rho_L
                                ! stagnation enthalpy here excludes magnetic energy (only used to find speed of sound)
                                H_R = (E_R + pres_R - pres_mag%R)/rho_R
                            else
                                E_L = gamma_L*pres_L + pi_inf_L + 5.e-1*rho_L*vel_L_rms + qv_L
                                E_R = gamma_R*pres_R + pi_inf_R + 5.e-1*rho_R*vel_R_rms + qv_R
                                H_L = (E_L + pres_L)/rho_L
                                H_R = (E_R + pres_R)/rho_R
                            end if

                            ! elastic energy update
                            if (hypoelasticity) then
                                G_L = 0._wp; G_R = 0._wp

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    G_L = G_L + alpha_L(i)*Gs_rs(i)
                                    G_R = G_R + alpha_R(i)*Gs_rs(i)
                                end do

                                if (cont_damage) then
                                    G_L = G_L*max((1._wp - qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%damage)), 0._wp)
                                    G_R = G_R*max((1._wp - qR_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%damage)), 0._wp)
                                end if

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, eqn_idx%stress%end - eqn_idx%stress%beg + 1
                                    tau_e_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%stress%beg - 1 + i)
                                    tau_e_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%stress%beg - 1 + i)
                                    ! Elastic energy (guard skips when G near zero)
                                    if (.not. hypo_energy_guard .or. ((G_L > verysmall) .and. (G_R > verysmall))) then
                                        E_L = E_L + (tau_e_L(i)*tau_e_L(i))/max(4._wp*G_L, verysmall)
                                        E_R = E_R + (tau_e_R(i)*tau_e_R(i))/max(4._wp*G_R, verysmall)
                                        ! Double for shear stresses
                                        if (any(eqn_idx%stress%beg - 1 + i == shear_indices)) then
                                            E_L = E_L + (tau_e_L(i)*tau_e_L(i))/max(4._wp*G_L, verysmall)
                                            E_R = E_R + (tau_e_R(i)*tau_e_R(i))/max(4._wp*G_R, verysmall)
                                        end if
                                    end if
                                end do
                            end if

                            @:compute_average_state()

                            call s_compute_speed_of_sound(pres_L, rho_L, gamma_L, pi_inf_L, H_L, alpha_L, vel_L_rms, 0._wp, c_L, &
                                                          & qv_L)

                            call s_compute_speed_of_sound(pres_R, rho_R, gamma_R, pi_inf_R, H_R, alpha_R, vel_R_rms, 0._wp, c_R, &
                                                          & qv_R)

                            !> The computation of c_avg does not require all the variables, and therefore the non '_avg'
                            ! variables are placeholders to call the subroutine.

                            call s_compute_speed_of_sound(pres_R, rho_avg, gamma_avg, pi_inf_R, H_avg, alpha_R, vel_avg_rms, &
                                                          & c_sum_Yi_Phi, c_avg, qv_avg)

                            if (mhd) then
                                call s_compute_fast_magnetosonic_speed(rho_L, c_L, B%L, norm_dir, c_fast%L, H_L)
                                call s_compute_fast_magnetosonic_speed(rho_R, c_R, B%R, norm_dir, c_fast%R, H_R)
                            end if

                            if (viscous) then
                                if (chemistry) then
                                    call compute_viscosity_and_inversion(T_L, Ys_L, T_R, Ys_R, Re_L(1), Re_R(1))
                                end if
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, 2
                                    Re_avg_rs${XYZ}$_vf(j, k, l, i) = 2._wp/(1._wp/Re_L(i) + 1._wp/Re_R(i))
                                end do
                            end if

                            ! Wave speed estimates (wave_speeds=1: direct, wave_speeds=2: pressure-based)
                            if (wave_speeds == 1) then
                                if (mhd) then
                                    ! MHD: use fast magnetosonic speed
                                    s_L = min(vel_L(dir_idx(1)) - c_fast%L, vel_R(dir_idx(1)) - c_fast%R)
                                    s_R = max(vel_R(dir_idx(1)) + c_fast%R, vel_L(dir_idx(1)) + c_fast%L)
                                else if (hypoelasticity) then
                                    ! Elastic wave speed, Rodriguez et al. JCP (2019)
                                    s_L = min(vel_L(dir_idx(1)) - sqrt(max(verysmall, &
                                              & c_L*c_L + (((4._wp*G_L)/3._wp) + tau_e_L(dir_idx_tau(1)))/rho_L)), &
                                              & vel_R(dir_idx(1)) - sqrt(max(verysmall, &
                                              & c_R*c_R + (((4._wp*G_R)/3._wp) + tau_e_R(dir_idx_tau(1)))/rho_R)))
                                    s_R = max(vel_R(dir_idx(1)) + sqrt(max(verysmall, &
                                              & c_R*c_R + (((4._wp*G_R)/3._wp) + tau_e_R(dir_idx_tau(1)))/rho_R)), &
                                              & vel_L(dir_idx(1)) + sqrt(max(verysmall, &
                                              & c_L*c_L + (((4._wp*G_L)/3._wp) + tau_e_L(dir_idx_tau(1)))/rho_L)))
                                else if (hyperelasticity) then
                                    s_L = min(vel_L(dir_idx(1)) - sqrt(c_L*c_L + (4._wp*G_L/3._wp)/rho_L), &
                                              & vel_R(dir_idx(1)) - sqrt(c_R*c_R + (4._wp*G_R/3._wp)/rho_R))
                                    s_R = max(vel_R(dir_idx(1)) + sqrt(c_R*c_R + (4._wp*G_R/3._wp)/rho_R), &
                                              & vel_L(dir_idx(1)) + sqrt(c_L*c_L + (4._wp*G_L/3._wp)/rho_L))
                                else
                                    s_L = min(vel_L(dir_idx(1)) - c_L, vel_R(dir_idx(1)) - c_R)
                                    s_R = max(vel_R(dir_idx(1)) + c_R, vel_L(dir_idx(1)) + c_L)
                                end if

                                if (hyper_cleaning) then
                                    ! Dedner GLM divergence cleaning, Dedner et al. JCP (2002)
                                    s_L = min(s_L, -hyper_cleaning_speed)
                                    s_R = max(s_R, hyper_cleaning_speed)
                                end if

                                s_S = (pres_R - pres_L + rho_L*vel_L(dir_idx(1))*(s_L - vel_L(dir_idx(1))) &
                                       & - rho_R*vel_R(dir_idx(1))*(s_R - vel_R(dir_idx(1))))/(rho_L*(s_L - vel_L(dir_idx(1))) &
                                       & - rho_R*(s_R - vel_R(dir_idx(1))))
                            else if (wave_speeds == 2) then
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

                            s_M = min(0._wp, s_L); s_P = max(0._wp, s_R)

                            xi_M = (5.e-1_wp + sign(5.e-1_wp, s_L)) + (5.e-1_wp - sign(5.e-1_wp, s_L))*(5.e-1_wp + sign(5.e-1_wp, &
                                    & s_R))
                            xi_P = (5.e-1_wp - sign(5.e-1_wp, s_R)) + (5.e-1_wp - sign(5.e-1_wp, s_L))*(5.e-1_wp + sign(5.e-1_wp, &
                                    & s_R))

                            ! HLL intercell flux: F* = (s_R*F_L - s_L*F_R + s_L*s_R*(U_R - U_L)) / (s_R - s_L) Low Mach correction
                            if (low_Mach == 1) then
                                @:compute_low_Mach_correction()
                            else
                                pcorr = 0._wp
                            end if

                            ! Mass
                            if (.not. relativity) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, eqn_idx%cont%end
                                    flux_rs${XYZ}$_vf(j, k, l, &
                                                      & i) = (s_M*alpha_rho_R(i)*vel_R(norm_dir) - s_P*alpha_rho_L(i) &
                                                      & *vel_L(norm_dir) + s_M*s_P*(alpha_rho_L(i) - alpha_rho_R(i)))/(s_M - s_P)
                                end do
                            else if (relativity) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, eqn_idx%cont%end
                                    flux_rs${XYZ}$_vf(j, k, l, &
                                                      & i) = (s_M*Ga%R*alpha_rho_R(i)*vel_R(norm_dir) - s_P*Ga%L*alpha_rho_L(i) &
                                                      & *vel_L(norm_dir) + s_M*s_P*(Ga%L*alpha_rho_L(i) - Ga%R*alpha_rho_R(i))) &
                                                      & /(s_M - s_P)
                                end do
                            end if

                            ! Momentum
                            if (mhd .and. (.not. relativity)) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, 3
                                    ! Flux of rho*v_i in the ${XYZ}$ direction = rho * v_i * v_${XYZ}$ - B_i * B_${XYZ}$ +
                                    ! delta_(${XYZ}$,i) * p_tot
                                    flux_rs${XYZ}$_vf(j, k, l, &
                                                      & eqn_idx%cont%end + i) = (s_M*(rho_R*vel_R(i)*vel_R(norm_dir) - B%R(i) &
                                                      & *B%R(norm_dir) + dir_flg(i)*(pres_R + pres_mag%R)) - s_P*(rho_L*vel_L(i) &
                                                      & *vel_L(norm_dir) - B%L(i)*B%L(norm_dir) + dir_flg(i)*(pres_L + pres_mag%L) &
                                                      & ) + s_M*s_P*(rho_L*vel_L(i) - rho_R*vel_R(i)))/(s_M - s_P)
                                end do
                            else if (mhd .and. relativity) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, 3
                                    ! Flux of m_i in the ${XYZ}$ direction = m_i * v_${XYZ}$ - b_i/Gamma * B_${XYZ}$ +
                                    ! delta_(${XYZ}$,i) * p_tot
                                    flux_rs${XYZ}$_vf(j, k, l, &
                                                      & eqn_idx%cont%end + i) = (s_M*(cm%R(i)*vel_R(norm_dir) - b4%R(i) &
                                                      & /Ga%R*B%R(norm_dir) + dir_flg(i)*(pres_R + pres_mag%R)) - s_P*(cm%L(i) &
                                                      & *vel_L(norm_dir) - b4%L(i)/Ga%L*B%L(norm_dir) + dir_flg(i)*(pres_L &
                                                      & + pres_mag%L)) + s_M*s_P*(cm%L(i) - cm%R(i)))/(s_M - s_P)
                                end do
                            else if (bubbles_euler) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_vels
                                    flux_rs${XYZ}$_vf(j, k, l, &
                                                      & eqn_idx%cont%end + dir_idx(i)) = (s_M*(rho_R*vel_R(dir_idx(1)) &
                                                      & *vel_R(dir_idx(i)) + dir_flg(dir_idx(i))*(pres_R - ptilde_R)) &
                                                      & - s_P*(rho_L*vel_L(dir_idx(1))*vel_L(dir_idx(i)) + dir_flg(dir_idx(i)) &
                                                      & *(pres_L - ptilde_L)) + s_M*s_P*(rho_L*vel_L(dir_idx(i)) &
                                                      & - rho_R*vel_R(dir_idx(i))))/(s_M - s_P) + (s_M/s_L)*(s_P/s_R) &
                                                      & *pcorr*(vel_R(dir_idx(i)) - vel_L(dir_idx(i)))
                                end do
                            else if (hypoelasticity) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_vels
                                    flux_rs${XYZ}$_vf(j, k, l, &
                                                      & eqn_idx%cont%end + dir_idx(i)) = (s_M*(rho_R*vel_R(dir_idx(1)) &
                                                      & *vel_R(dir_idx(i)) + dir_flg(dir_idx(i))*pres_R - tau_e_R(dir_idx_tau(i))) &
                                                      & - s_P*(rho_L*vel_L(dir_idx(1))*vel_L(dir_idx(i)) + dir_flg(dir_idx(i)) &
                                                      & *pres_L - tau_e_L(dir_idx_tau(i))) + s_M*s_P*(rho_L*vel_L(dir_idx(i)) &
                                                      & - rho_R*vel_R(dir_idx(i))))/(s_M - s_P)
                                end do
                            else
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_vels
                                    flux_rs${XYZ}$_vf(j, k, l, &
                                                      & eqn_idx%cont%end + dir_idx(i)) = (s_M*(rho_R*vel_R(dir_idx(1)) &
                                                      & *vel_R(dir_idx(i)) + dir_flg(dir_idx(i))*pres_R) &
                                                      & - s_P*(rho_L*vel_L(dir_idx(1))*vel_L(dir_idx(i)) + dir_flg(dir_idx(i)) &
                                                      & *pres_L) + s_M*s_P*(rho_L*vel_L(dir_idx(i)) - rho_R*vel_R(dir_idx(i)))) &
                                                      & /(s_M - s_P) + (s_M/s_L)*(s_P/s_R)*pcorr*(vel_R(dir_idx(i)) &
                                                      & - vel_L(dir_idx(i)))
                                end do
                            end if

                            ! Energy
                            if (mhd .and. (.not. relativity)) then
                                ! energy flux = (E + p + p_mag) * v_${XYZ}$ - B_${XYZ}$ * (v_x*B_x + v_y*B_y + v_z*B_z)
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                                    flux_rs${XYZ}$_vf(j, k, l, &
                                                      & eqn_idx%E) = (s_M*(vel_R(norm_dir)*(E_R + pres_R + pres_mag%R) &
                                                      & - B%R(norm_dir)*(vel_R(1)*B%R(1) + vel_R(2)*B%R(2) + vel_R(3)*B%R(3))) &
                                                      & - s_P*(vel_L(norm_dir)*(E_L + pres_L + pres_mag%L) - B%L(norm_dir) &
                                                      & *(vel_L(1)*B%L(1) + vel_L(2)*B%L(2) + vel_L(3)*B%L(3))) + s_M*s_P*(E_L &
                                                      & - E_R))/(s_M - s_P)
                                #:endif
                            else if (mhd .and. relativity) then
                                ! energy flux = m_${XYZ}$ - mass flux Hard-coded for single-component for now
                                flux_rs${XYZ}$_vf(j, k, l, &
                                                  & eqn_idx%E) = (s_M*(cm%R(norm_dir) - Ga%R*alpha_rho_R(1)*vel_R(norm_dir)) &
                                                  & - s_P*(cm%L(norm_dir) - Ga%L*alpha_rho_L(1)*vel_L(norm_dir)) + s_M*s_P*(E_L &
                                                  & - E_R))/(s_M - s_P)
                            else if (bubbles_euler) then
                                flux_rs${XYZ}$_vf(j, k, l, &
                                                  & eqn_idx%E) = (s_M*vel_R(dir_idx(1))*(E_R + pres_R - ptilde_R) &
                                                  & - s_P*vel_L(dir_idx(1))*(E_L + pres_L - ptilde_L) + s_M*s_P*(E_L - E_R))/(s_M &
                                                  & - s_P) + (s_M/s_L)*(s_P/s_R)*pcorr*(vel_R_rms - vel_L_rms)/2._wp
                            else if (hypoelasticity) then
                                flux_tau_L = 0._wp; flux_tau_R = 0._wp
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    flux_tau_L = flux_tau_L + tau_e_L(dir_idx_tau(i))*vel_L(dir_idx(i))
                                    flux_tau_R = flux_tau_R + tau_e_R(dir_idx_tau(i))*vel_R(dir_idx(i))
                                end do
                                flux_rs${XYZ}$_vf(j, k, l, &
                                                  & eqn_idx%E) = (s_M*(vel_R(dir_idx(1))*(E_R + pres_R) - flux_tau_R) &
                                                  & - s_P*(vel_L(dir_idx(1))*(E_L + pres_L) - flux_tau_L) + s_M*s_P*(E_L - E_R)) &
                                                  & /(s_M - s_P)
                            else
                                flux_rs${XYZ}$_vf(j, k, l, &
                                                  & eqn_idx%E) = (s_M*vel_R(dir_idx(1))*(E_R + pres_R) - s_P*vel_L(dir_idx(1)) &
                                                  & *(E_L + pres_L) + s_M*s_P*(E_L - E_R))/(s_M - s_P) + (s_M/s_L)*(s_P/s_R) &
                                                  & *pcorr*(vel_R_rms - vel_L_rms)/2._wp
                            end if

                            ! Elastic Stresses
                            if (hypoelasticity) then
                                do i = 1, eqn_idx%stress%end - eqn_idx%stress%beg + 1  ! TODO: this indexing may be slow
                                    flux_rs${XYZ}$_vf(j, k, l, &
                                                      & eqn_idx%stress%beg - 1 + i) = (s_M*(rho_R*vel_R(dir_idx(1))*tau_e_R(i)) &
                                                      & - s_P*(rho_L*vel_L(dir_idx(1))*tau_e_L(i)) + s_M*s_P*(rho_L*tau_e_L(i) &
                                                      & - rho_R*tau_e_R(i)))/(s_M - s_P)
                                end do
                            end if

                            ! Export interface velocity for NC RHS
                            if (hypo_nc_interface .or. (alt_soundspeed .and. .not. hll_u_interface)) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    if (0._wp <= s_L) then
                                        nc_iface_vel_rs${XYZ}$_vf(j, k, l, dir_idx(i)) = vel_L(dir_idx(i))
                                    else if (s_R <= 0._wp) then
                                        nc_iface_vel_rs${XYZ}$_vf(j, k, l, dir_idx(i)) = vel_R(dir_idx(i))
                                    else
                                        nc_iface_vel_rs${XYZ}$_vf(j, k, l, &
                                                                  & dir_idx(i)) = (s_R*vel_L(dir_idx(i)) - s_L*vel_R(dir_idx(i))) &
                                                                  & /(s_R - s_L)
                                    end if
                                end do
                            end if

                            if (.not. hll_u_interface) then  ! HLL Method 1: per-fluid alpha interface flux
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = eqn_idx%adv%beg, eqn_idx%adv%end
                                    if (0._wp <= s_L .or. s_R <= 0._wp) then
                                        flux_rs${XYZ}$_vf(j, k, l, i) = 0._wp
                                    else
                                        flux_rs${XYZ}$_vf(j, k, l, i) = s_M*s_P*(qR_prim_rs${XYZ}$_vf(j + 1, k, l, &
                                                          & i) - qL_prim_rs${XYZ}$_vf(j, k, l, i))/(s_R - s_L)
                                    end if
                                    if (0._wp <= s_L) then
                                        flux_src_rs${XYZ}$_vf(j, k, l, i) = qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                    else if (s_R <= 0._wp) then
                                        flux_src_rs${XYZ}$_vf(j, k, l, i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)
                                    else
                                        flux_src_rs${XYZ}$_vf(j, k, l, i) = (s_R*qL_prim_rs${XYZ}$_vf(j, k, l, &
                                                              & i) - s_L*qR_prim_rs${XYZ}$_vf(j + 1, k, l, i))/(s_R - s_L)
                                    end if
                                end do
                            else  ! HLL Method 2: shared velocity interface flux
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = eqn_idx%adv%beg, eqn_idx%adv%end
                                    if (0._wp <= s_L) then
                                        flux_rs${XYZ}$_vf(j, k, l, i) = qL_prim_rs${XYZ}$_vf(j, k, l, i)*vel_L(dir_idx(1))
                                    else if (s_R <= 0._wp) then
                                        flux_rs${XYZ}$_vf(j, k, l, i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)*vel_R(dir_idx(1))
                                    else
                                        flux_rs${XYZ}$_vf(j, k, l, i) = (s_R*qL_prim_rs${XYZ}$_vf(j, k, l, &
                                                          & i)*vel_L(dir_idx(1)) - s_L*qR_prim_rs${XYZ}$_vf(j + 1, k, l, &
                                                          & i)*vel_R(dir_idx(1)) + s_L*s_R*(qR_prim_rs${XYZ}$_vf(j + 1, k, l, &
                                                          & i) - qL_prim_rs${XYZ}$_vf(j, k, l, i)))/(s_R - s_L)
                                    end if
                                end do
                                if (0._wp <= s_L) then
                                    flux_src_rs${XYZ}$_vf(j, k, l, eqn_idx%adv%beg) = vel_L(dir_idx(1))
                                else if (s_R <= 0._wp) then
                                    flux_src_rs${XYZ}$_vf(j, k, l, eqn_idx%adv%beg) = vel_R(dir_idx(1))
                                else
                                    flux_src_rs${XYZ}$_vf(j, k, l, &
                                                          & eqn_idx%adv%beg) = (s_R*vel_L(dir_idx(1)) - s_L*vel_R(dir_idx(1))) &
                                                          & /(s_R - s_L)
                                end if
                            end if

                            if (bubbles_euler) then
                                ! From HLLC: Kills mass transport @ bubble gas density
                                if (num_fluids > 1) then
                                    flux_rs${XYZ}$_vf(j, k, l, eqn_idx%cont%end) = 0._wp
                                end if
                            end if

                            if (chemistry) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = eqn_idx%species%beg, eqn_idx%species%end
                                    Y_L = qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                    Y_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)

                                    flux_rs${XYZ}$_vf(j, k, l, &
                                                      & i) = (s_M*Y_R*rho_R*vel_R(dir_idx(1)) - s_P*Y_L*rho_L*vel_L(dir_idx(1)) &
                                                      & + s_M*s_P*(Y_L*rho_L - Y_R*rho_R))/(s_M - s_P)
                                    flux_src_rs${XYZ}$_vf(j, k, l, i) = 0._wp
                                end do
                            end if

                            ! MHD: magnetic flux and Maxwell stress contributions
                            if (mhd) then
                                if (n == 0) then  ! 1D: d/dx flux only & Bx = Bx0 = const.
                                    ! B_y flux = v_x * B_y - v_y * Bx0 B_z flux = v_x * B_z - v_z * Bx0
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 0, 1
                                        flux_rsx_vf(j, k, l, &
                                                    & eqn_idx%B%beg + i) = (s_M*(vel_R(1)*B%R(2 + i) - vel_R(2 + i)*Bx0) &
                                                    & - s_P*(vel_L(1)*B%L(2 + i) - vel_L(2 + i)*Bx0) + s_M*s_P*(B%L(2 + i) &
                                                    & - B%R(2 + i)))/(s_M - s_P)
                                    end do
                                else  ! 2D/3D: Bx, By, Bz /= const. but zero flux component in the same direction
                                    ! B_x d/d${XYZ}$ flux = (1 - delta(x,${XYZ}$)) * (v_${XYZ}$ * B_x - v_x * B_${XYZ}$) B_y
                                    ! d/d${XYZ}$ flux = (1 - delta(y,${XYZ}$)) * (v_${XYZ}$ * B_y - v_y * B_${XYZ}$) B_z d/d${XYZ}$
                                    ! flux = (1 - delta(z,${XYZ}$)) * (v_${XYZ}$ * B_z - v_z * B_${XYZ}$)
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 0, 2
                                        flux_rs${XYZ}$_vf(j, k, l, &
                                                          & eqn_idx%B%beg + i) = (s_M*(vel_R(dir_idx(1))*B%R(i + 1) - vel_R(i + 1) &
                                                          & *B%R(norm_dir)) - s_P*(vel_L(dir_idx(1))*B%L(i + 1) - vel_L(i + 1) &
                                                          & *B%L(norm_dir)) + s_M*s_P*(B%L(i + 1) - B%R(i + 1)))/(s_M - s_P)
                                    end do

                                    if (hyper_cleaning) then
                                        ! propagate magnetic field divergence as a wave
                                        flux_rs${XYZ}$_vf(j, k, l, eqn_idx%B%beg + norm_dir - 1) = flux_rs${XYZ}$_vf(j, k, l, &
                                                          & eqn_idx%B%beg + norm_dir - 1) + (s_M*qR_prim_rs${XYZ}$_vf(j + 1, k, &
                                                          & l, eqn_idx%psi) - s_P*qL_prim_rs${XYZ}$_vf(j, k, l, &
                                                          & eqn_idx%psi))/(s_M - s_P)

                                        flux_rs${XYZ}$_vf(j, k, l, &
                                                          & eqn_idx%psi) = (hyper_cleaning_speed**2*(s_M*B%R(norm_dir) &
                                                          & - s_P*B%L(norm_dir)) + s_M*s_P*(qL_prim_rs${XYZ}$_vf(j, k, l, &
                                                          & eqn_idx%psi) - qR_prim_rs${XYZ}$_vf(j + 1, k, l, &
                                                          & eqn_idx%psi)))/(s_M - s_P)
                                    else
                                        ! Without hyperbolic cleaning, make sure flux of B_normal is identically zero
                                        flux_rs${XYZ}$_vf(j, k, l, eqn_idx%B%beg + norm_dir - 1) = 0._wp
                                    end if
                                end if
                                flux_src_rs${XYZ}$_vf(j, k, l, eqn_idx%adv%beg) = 0._wp
                            end if

                            #:if (NORM_DIR == 2)
                                if (cyl_coord) then
                                    ! Substituting the advective flux into the inviscid geometrical source flux
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, eqn_idx%E
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = flux_rs${XYZ}$_vf(j, k, l, i)
                                    end do
                                    ! Recalculating the radial momentum geometric source flux
                                    flux_gsrc_rs${XYZ}$_vf(j, k, l, eqn_idx%cont%end + 2) = flux_rs${XYZ}$_vf(j, k, l, &
                                                           & eqn_idx%cont%end + 2) - (s_M*pres_R - s_P*pres_L)/(s_M - s_P)
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = eqn_idx%adv%beg, eqn_idx%adv%end
                                        if (.not. hll_u_interface) then
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = flux_rs${XYZ}$_vf(j, k, l, i)
                                        else
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = 0._wp
                                        end if
                                    end do
                                end if

                                if (cyl_coord .and. hypoelasticity) then
                                    ! += tau_sigmasigma using HLL
                                    flux_gsrc_rs${XYZ}$_vf(j, k, l, eqn_idx%cont%end + 2) = flux_gsrc_rs${XYZ}$_vf(j, k, l, &
                                                           & eqn_idx%cont%end + 2) + (s_M*tau_e_R(4) - s_P*tau_e_L(4))/(s_M - s_P)

                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = eqn_idx%stress%beg, eqn_idx%stress%end
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = flux_rs${XYZ}$_vf(j, k, l, i)
                                    end do
                                end if
                            #:endif
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

        if (viscous .or. dummy) then
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
    subroutine s_lf_riemann_solver(qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, dqL_prim_dx_vf, dqL_prim_dy_vf, &

        & dqL_prim_dz_vf, qL_prim_vf, qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, dqR_prim_dx_vf, dqR_prim_dy_vf, &
            & dqR_prim_dz_vf, qR_prim_vf, q_prim_vf, flux_vf, flux_src_vf, flux_gsrc_vf, norm_dir, ix, iy, iz)

        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout) :: qL_prim_rsx_vf, qL_prim_rsy_vf, &
             & qL_prim_rsz_vf, qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf
        type(scalar_field), dimension(sys_size), intent(in)          :: q_prim_vf
        type(scalar_field), allocatable, dimension(:), intent(inout) :: qL_prim_vf, qR_prim_vf
        type(scalar_field), allocatable, dimension(:), intent(inout) :: dqL_prim_dx_vf, dqR_prim_dx_vf, dqL_prim_dy_vf, &
             & dqR_prim_dy_vf, dqL_prim_dz_vf, dqR_prim_dz_vf

        ! Intercell fluxes
        type(scalar_field), dimension(sys_size), intent(inout) :: flux_vf, flux_src_vf, flux_gsrc_vf
        real(wp)                                               :: flux_tau_L, flux_tau_R
        integer, intent(in)                                    :: norm_dir
        type(int_bounds_info), intent(in)                      :: ix, iy, iz

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3)    :: alpha_rho_L, alpha_rho_R
            real(wp), dimension(3)    :: vel_L, vel_R
            real(wp), dimension(3)    :: alpha_L, alpha_R
            real(wp), dimension(10)   :: Ys_L, Ys_R
            real(wp), dimension(10)   :: Cp_iL, Cp_iR, Xs_L, Xs_R, Gamma_iL, Gamma_iR
            real(wp), dimension(10)   :: Yi_avg, Phi_avg, h_iL, h_iR, h_avg_2
            real(wp), dimension(3, 3) :: vel_grad_L, vel_grad_R  !< Averaged velocity gradient tensor `d(vel_i)/d(coord_j)`.
        #:else
            real(wp), dimension(num_fluids)  :: alpha_rho_L, alpha_rho_R
            real(wp), dimension(num_vels)    :: vel_L, vel_R
            real(wp), dimension(num_fluids)  :: alpha_L, alpha_R
            real(wp), dimension(num_species) :: Ys_L, Ys_R
            real(wp), dimension(num_species) :: Cp_iL, Cp_iR, Xs_L, Xs_R, Gamma_iL, Gamma_iR
            real(wp), dimension(num_species) :: Yi_avg, Phi_avg, h_iL, h_iR, h_avg_2
            !> Averaged velocity gradient tensor `d(vel_i)/d(coord_j)`.
            real(wp), dimension(num_dims, num_dims) :: vel_grad_L, vel_grad_R
        #:endif
        real(wp)                  :: rho_L, rho_R
        real(wp)                  :: pres_L, pres_R
        real(wp)                  :: E_L, E_R
        real(wp)                  :: H_L, H_R
        real(wp)                  :: Cp_avg, Cv_avg, T_avg, eps, c_sum_Yi_Phi
        real(wp)                  :: T_L, T_R
        real(wp)                  :: Y_L, Y_R
        real(wp)                  :: MW_L, MW_R
        real(wp)                  :: R_gas_L, R_gas_R
        real(wp)                  :: Cp_L, Cp_R
        real(wp)                  :: Cv_L, Cv_R
        real(wp)                  :: Gamm_L, Gamm_R
        real(wp)                  :: gamma_L, gamma_R
        real(wp)                  :: pi_inf_L, pi_inf_R
        real(wp)                  :: qv_L, qv_R
        real(wp)                  :: c_L, c_R
        real(wp), dimension(6)    :: tau_e_L, tau_e_R
        real(wp)                  :: G_L, G_R
        real(wp), dimension(2)    :: Re_L, Re_R
        real(wp), dimension(3)    :: xi_field_L, xi_field_R
        real(wp)                  :: rho_avg
        real(wp)                  :: H_avg
        real(wp)                  :: gamma_avg
        real(wp)                  :: c_avg
        real(wp)                  :: s_L, s_R, s_M, s_P, s_S
        real(wp)                  :: xi_M, xi_P
        real(wp)                  :: ptilde_L, ptilde_R
        real(wp)                  :: vel_L_rms, vel_R_rms, vel_avg_rms
        real(wp)                  :: vel_L_tmp, vel_R_tmp
        real(wp)                  :: Ms_L, Ms_R, pres_SL, pres_SR
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

        call s_populate_riemann_states_variables_buffers(qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, dqL_prim_dx_vf, &
            & dqL_prim_dy_vf, dqL_prim_dz_vf, qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, dqR_prim_dx_vf, dqR_prim_dy_vf, &
            & dqR_prim_dz_vf, norm_dir, ix, iy, iz)

        ! Reshaping inputted data based on dimensional splitting direction
        call s_initialize_riemann_solver(flux_src_vf, norm_dir)
        #:for NORM_DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
            if (norm_dir == ${NORM_DIR}$) then
                $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, q, alpha_rho_L, alpha_rho_R, vel_L, vel_R, alpha_L, &
                                    & alpha_R, tau_e_L, tau_e_R, G_L, G_R, Re_L, Re_R, rho_avg, h_avg, gamma_avg, s_L, s_R, s_S, &
                                    & Ys_L, Ys_R, xi_field_L, xi_field_R, Cp_iL, Cp_iR, Xs_L, Xs_R, Gamma_iL, Gamma_iR, Yi_avg, &
                                    & Phi_avg, h_iL, h_iR, h_avg_2, c_fast, pres_mag, B, Ga, vdotB, B2, b4, cm, pcorr, zcoef, &
                                    & vel_grad_L, vel_grad_R, idx_right_phys, vel_L_rms, vel_R_rms, vel_avg_rms, vel_L_tmp, &
                                    & vel_R_tmp, Ms_L, Ms_R, pres_SL, pres_SR, alpha_L_sum, alpha_R_sum, c_avg, pres_L, pres_R, &
                                    & rho_L, rho_R, gamma_L, gamma_R, pi_inf_L, pi_inf_R, qv_L, qv_R, c_L, c_R, E_L, E_R, H_L, &
                                    & H_R, ptilde_L, ptilde_R, s_M, s_P, xi_M, xi_P, Cp_avg, Cv_avg, T_avg, eps, c_sum_Yi_Phi, &
                                    & Cp_L, Cp_R, Cv_L, Cv_R, R_gas_L, R_gas_R, MW_L, MW_R, T_L, T_R, Y_L, Y_R]')
                do l = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = is1%beg, is1%end
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, eqn_idx%cont%end
                                alpha_rho_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                alpha_rho_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)
                            end do

                            vel_L_rms = 0._wp; vel_R_rms = 0._wp

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_vels
                                vel_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%cont%end + i)
                                vel_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%cont%end + i)
                                vel_L_rms = vel_L_rms + vel_L(i)**2._wp
                                vel_R_rms = vel_R_rms + vel_R(i)**2._wp
                            end do

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                alpha_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%E + i)
                                alpha_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%E + i)
                            end do

                            pres_L = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%E)
                            pres_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%E)

                            if (mhd) then
                                if (n == 0) then  ! 1D: constant Bx; By, Bz as variables
                                    B%L(1) = Bx0
                                    B%R(1) = Bx0
                                    B%L(2) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%B%beg)
                                    B%R(2) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%B%beg)
                                    B%L(3) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%B%beg + 1)
                                    B%R(3) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%B%beg + 1)
                                else  ! 2D/3D: Bx, By, Bz as variables
                                    B%L(1) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%B%beg)
                                    B%R(1) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%B%beg)
                                    B%L(2) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%B%beg + 1)
                                    B%R(2) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%B%beg + 1)
                                    B%L(3) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%B%beg + 2)
                                    B%R(3) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%B%beg + 2)
                                end if
                            end if

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
                                rho_L = rho_L + alpha_rho_L(i)
                                gamma_L = gamma_L + alpha_L(i)*gammas(i)
                                pi_inf_L = pi_inf_L + alpha_L(i)*pi_infs(i)
                                qv_L = qv_L + alpha_rho_L(i)*qvs(i)

                                rho_R = rho_R + alpha_rho_R(i)
                                gamma_R = gamma_R + alpha_R(i)*gammas(i)
                                pi_inf_R = pi_inf_R + alpha_R(i)*pi_infs(i)
                                qv_R = qv_R + alpha_rho_R(i)*qvs(i)
                            end do

                            if (viscous) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, 2
                                    Re_L(i) = dflt_real
                                    Re_R(i) = dflt_real

                                    if (Re_size(i) > 0) Re_L(i) = 0._wp
                                    if (Re_size(i) > 0) Re_R(i) = 0._wp

                                    $:GPU_LOOP(parallelism='[seq]')
                                    do q = 1, Re_size(i)
                                        Re_L(i) = alpha_L(Re_idx(i, q))/Res_gs(i, q) + Re_L(i)
                                        Re_R(i) = alpha_R(Re_idx(i, q))/Res_gs(i, q) + Re_R(i)
                                    end do

                                    Re_L(i) = 1._wp/max(Re_L(i), sgm_eps)
                                    Re_R(i) = 1._wp/max(Re_R(i), sgm_eps)
                                end do
                            end if

                            if (chemistry) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = eqn_idx%species%beg, eqn_idx%species%end
                                    Ys_L(i - eqn_idx%species%beg + 1) = qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                    Ys_R(i - eqn_idx%species%beg + 1) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)
                                end do

                                call get_mixture_molecular_weight(Ys_L, MW_L)
                                call get_mixture_molecular_weight(Ys_R, MW_R)

                                #:if USING_AMD
                                    Xs_L(:) = Ys_L(:)*MW_L/molecular_weights_nonparameter(:)
                                    Xs_R(:) = Ys_R(:)*MW_R/molecular_weights_nonparameter(:)
                                #:else
                                    Xs_L(:) = Ys_L(:)*MW_L/molecular_weights(:)
                                    Xs_R(:) = Ys_R(:)*MW_R/molecular_weights(:)
                                #:endif

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
                                H_L = (E_L + pres_L)/rho_L
                                H_R = (E_R + pres_R)/rho_R
                            else if (mhd .and. relativity) then
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                                    Ga%L = 1._wp/sqrt(1._wp - vel_L_rms)
                                    Ga%R = 1._wp/sqrt(1._wp - vel_R_rms)
                                    vdotB%L = vel_L(1)*B%L(1) + vel_L(2)*B%L(2) + vel_L(3)*B%L(3)
                                    vdotB%R = vel_R(1)*B%R(1) + vel_R(2)*B%R(2) + vel_R(3)*B%R(3)

                                    b4%L(1:3) = B%L(1:3)/Ga%L + Ga%L*vel_L(1:3)*vdotB%L
                                    b4%R(1:3) = B%R(1:3)/Ga%R + Ga%R*vel_R(1:3)*vdotB%R
                                    B2%L = B%L(1)**2._wp + B%L(2)**2._wp + B%L(3)**2._wp
                                    B2%R = B%R(1)**2._wp + B%R(2)**2._wp + B%R(3)**2._wp

                                    pres_mag%L = 0.5_wp*(B2%L/Ga%L**2._wp + vdotB%L**2._wp)
                                    pres_mag%R = 0.5_wp*(B2%R/Ga%R**2._wp + vdotB%R**2._wp)

                                    ! Hard-coded EOS
                                    H_L = 1._wp + (gamma_L + 1)*pres_L/rho_L
                                    H_R = 1._wp + (gamma_R + 1)*pres_R/rho_R

                                    cm%L(1:3) = (rho_L*H_L*Ga%L**2 + B2%L)*vel_L(1:3) - vdotB%L*B%L(1:3)
                                    cm%R(1:3) = (rho_R*H_R*Ga%R**2 + B2%R)*vel_R(1:3) - vdotB%R*B%R(1:3)

                                    E_L = rho_L*H_L*Ga%L**2 - pres_L + 0.5_wp*(B2%L + vel_L_rms*B2%L - vdotB%L**2._wp) - rho_L*Ga%L
                                    E_R = rho_R*H_R*Ga%R**2 - pres_R + 0.5_wp*(B2%R + vel_R_rms*B2%R - vdotB%R**2._wp) - rho_R*Ga%R
                                #:endif
                            else if (mhd .and. .not. relativity) then
                                pres_mag%L = 0.5_wp*(B%L(1)**2._wp + B%L(2)**2._wp + B%L(3)**2._wp)
                                pres_mag%R = 0.5_wp*(B%R(1)**2._wp + B%R(2)**2._wp + B%R(3)**2._wp)
                                E_L = gamma_L*pres_L + pi_inf_L + 0.5_wp*rho_L*vel_L_rms + qv_L + pres_mag%L
                                ! includes magnetic energy
                                E_R = gamma_R*pres_R + pi_inf_R + 0.5_wp*rho_R*vel_R_rms + qv_R + pres_mag%R
                                H_L = (E_L + pres_L - pres_mag%L)/rho_L
                                ! stagnation enthalpy here excludes magnetic energy (only used to find speed of sound)
                                H_R = (E_R + pres_R - pres_mag%R)/rho_R
                            else
                                E_L = gamma_L*pres_L + pi_inf_L + 5.e-1*rho_L*vel_L_rms + qv_L
                                E_R = gamma_R*pres_R + pi_inf_R + 5.e-1*rho_R*vel_R_rms + qv_R
                                H_L = (E_L + pres_L)/rho_L
                                H_R = (E_R + pres_R)/rho_R
                            end if

                            ! elastic energy update
                            if (hypoelasticity) then
                                G_L = 0._wp; G_R = 0._wp

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    G_L = G_L + alpha_L(i)*Gs_rs(i)
                                    G_R = G_R + alpha_R(i)*Gs_rs(i)
                                end do

                                if (cont_damage) then
                                    G_L = G_L*max((1._wp - qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%damage)), 0._wp)
                                    G_R = G_R*max((1._wp - qR_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%damage)), 0._wp)
                                end if

                                do i = 1, eqn_idx%stress%end - eqn_idx%stress%beg + 1
                                    tau_e_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%stress%beg - 1 + i)
                                    tau_e_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%stress%beg - 1 + i)
                                    ! Elastic energy (guard skips when G near zero)
                                    if (.not. hypo_energy_guard .or. ((G_L > verysmall) .and. (G_R > verysmall))) then
                                        E_L = E_L + (tau_e_L(i)*tau_e_L(i))/max(4._wp*G_L, verysmall)
                                        E_R = E_R + (tau_e_R(i)*tau_e_R(i))/max(4._wp*G_R, verysmall)
                                        ! Double for shear stresses
                                        if (any(eqn_idx%stress%beg - 1 + i == shear_indices)) then
                                            E_L = E_L + (tau_e_L(i)*tau_e_L(i))/max(4._wp*G_L, verysmall)
                                            E_R = E_R + (tau_e_R(i)*tau_e_R(i))/max(4._wp*G_R, verysmall)
                                        end if
                                    end if
                                end do
                            end if

                            call s_compute_speed_of_sound(pres_L, rho_L, gamma_L, pi_inf_L, H_L, alpha_L, vel_L_rms, 0._wp, c_L, &
                                                          & qv_L)

                            call s_compute_speed_of_sound(pres_R, rho_R, gamma_R, pi_inf_R, H_R, alpha_R, vel_R_rms, 0._wp, c_R, &
                                                          & qv_R)

                            if (mhd) then
                                call s_compute_fast_magnetosonic_speed(rho_L, c_L, B%L, norm_dir, c_fast%L, H_L)
                                call s_compute_fast_magnetosonic_speed(rho_R, c_R, B%R, norm_dir, c_fast%R, H_R)
                            end if

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
                            if (low_Mach == 1) then
                                @:compute_low_Mach_correction()
                            else
                                pcorr = 0._wp
                            end if

                            ! Mass
                            if (.not. relativity) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, eqn_idx%cont%end
                                    flux_rs${XYZ}$_vf(j, k, l, &
                                                      & i) = (s_M*alpha_rho_R(i)*vel_R(norm_dir) - s_P*alpha_rho_L(i) &
                                                      & *vel_L(norm_dir) + s_M*s_P*(alpha_rho_L(i) - alpha_rho_R(i)))/(s_M - s_P)
                                end do
                            else if (relativity) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, eqn_idx%cont%end
                                    flux_rs${XYZ}$_vf(j, k, l, &
                                                      & i) = (s_M*Ga%R*alpha_rho_R(i)*vel_R(norm_dir) - s_P*Ga%L*alpha_rho_L(i) &
                                                      & *vel_L(norm_dir) + s_M*s_P*(Ga%L*alpha_rho_L(i) - Ga%R*alpha_rho_R(i))) &
                                                      & /(s_M - s_P)
                                end do
                            end if

                            ! Momentum
                            if (mhd .and. (.not. relativity)) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, 3
                                    ! Flux of rho*v_i in the ${XYZ}$ direction = rho * v_i * v_${XYZ}$ - B_i * B_${XYZ}$ +
                                    ! delta_(${XYZ}$,i) * p_tot
                                    flux_rs${XYZ}$_vf(j, k, l, &
                                                      & eqn_idx%cont%end + i) = (s_M*(rho_R*vel_R(i)*vel_R(norm_dir) - B%R(i) &
                                                      & *B%R(norm_dir) + dir_flg(i)*(pres_R + pres_mag%R)) - s_P*(rho_L*vel_L(i) &
                                                      & *vel_L(norm_dir) - B%L(i)*B%L(norm_dir) + dir_flg(i)*(pres_L + pres_mag%L) &
                                                      & ) + s_M*s_P*(rho_L*vel_L(i) - rho_R*vel_R(i)))/(s_M - s_P)
                                end do
                            else if (mhd .and. relativity) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, 3
                                    ! Flux of m_i in the ${XYZ}$ direction = m_i * v_${XYZ}$ - b_i/Gamma * B_${XYZ}$ +
                                    ! delta_(${XYZ}$,i) * p_tot
                                    flux_rs${XYZ}$_vf(j, k, l, &
                                                      & eqn_idx%cont%end + i) = (s_M*(cm%R(i)*vel_R(norm_dir) - b4%R(i) &
                                                      & /Ga%R*B%R(norm_dir) + dir_flg(i)*(pres_R + pres_mag%R)) - s_P*(cm%L(i) &
                                                      & *vel_L(norm_dir) - b4%L(i)/Ga%L*B%L(norm_dir) + dir_flg(i)*(pres_L &
                                                      & + pres_mag%L)) + s_M*s_P*(cm%L(i) - cm%R(i)))/(s_M - s_P)
                                end do
                            else if (bubbles_euler) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_vels
                                    flux_rs${XYZ}$_vf(j, k, l, &
                                                      & eqn_idx%cont%end + dir_idx(i)) = (s_M*(rho_R*vel_R(dir_idx(1)) &
                                                      & *vel_R(dir_idx(i)) + dir_flg(dir_idx(i))*(pres_R - ptilde_R)) &
                                                      & - s_P*(rho_L*vel_L(dir_idx(1))*vel_L(dir_idx(i)) + dir_flg(dir_idx(i)) &
                                                      & *(pres_L - ptilde_L)) + s_M*s_P*(rho_L*vel_L(dir_idx(i)) &
                                                      & - rho_R*vel_R(dir_idx(i))))/(s_M - s_P) + (s_M/s_L)*(s_P/s_R) &
                                                      & *pcorr*(vel_R(dir_idx(i)) - vel_L(dir_idx(i)))
                                end do
                            else if (hypoelasticity) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_vels
                                    flux_rs${XYZ}$_vf(j, k, l, &
                                                      & eqn_idx%cont%end + dir_idx(i)) = (s_M*(rho_R*vel_R(dir_idx(1)) &
                                                      & *vel_R(dir_idx(i)) + dir_flg(dir_idx(i))*pres_R - tau_e_R(dir_idx_tau(i))) &
                                                      & - s_P*(rho_L*vel_L(dir_idx(1))*vel_L(dir_idx(i)) + dir_flg(dir_idx(i)) &
                                                      & *pres_L - tau_e_L(dir_idx_tau(i))) + s_M*s_P*(rho_L*vel_L(dir_idx(i)) &
                                                      & - rho_R*vel_R(dir_idx(i))))/(s_M - s_P)
                                end do
                            else
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_vels
                                    flux_rs${XYZ}$_vf(j, k, l, &
                                                      & eqn_idx%cont%end + dir_idx(i)) = (s_M*(rho_R*vel_R(dir_idx(1)) &
                                                      & *vel_R(dir_idx(i)) + dir_flg(dir_idx(i))*pres_R) &
                                                      & - s_P*(rho_L*vel_L(dir_idx(1))*vel_L(dir_idx(i)) + dir_flg(dir_idx(i)) &
                                                      & *pres_L) + s_M*s_P*(rho_L*vel_L(dir_idx(i)) - rho_R*vel_R(dir_idx(i)))) &
                                                      & /(s_M - s_P) + (s_M/s_L)*(s_P/s_R)*pcorr*(vel_R(dir_idx(i)) &
                                                      & - vel_L(dir_idx(i)))
                                end do
                            end if

                            ! Energy
                            if (mhd .and. (.not. relativity)) then
                                ! energy flux = (E + p + p_mag) * v_${XYZ}$ - B_${XYZ}$ * (v_x*B_x + v_y*B_y + v_z*B_z)
                                #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
                                    flux_rs${XYZ}$_vf(j, k, l, &
                                                      & eqn_idx%E) = (s_M*(vel_R(norm_dir)*(E_R + pres_R + pres_mag%R) &
                                                      & - B%R(norm_dir)*(vel_R(1)*B%R(1) + vel_R(2)*B%R(2) + vel_R(3)*B%R(3))) &
                                                      & - s_P*(vel_L(norm_dir)*(E_L + pres_L + pres_mag%L) - B%L(norm_dir) &
                                                      & *(vel_L(1)*B%L(1) + vel_L(2)*B%L(2) + vel_L(3)*B%L(3))) + s_M*s_P*(E_L &
                                                      & - E_R))/(s_M - s_P)
                                #:endif
                            else if (mhd .and. relativity) then
                                ! energy flux = m_${XYZ}$ - mass flux Hard-coded for single-component for now
                                flux_rs${XYZ}$_vf(j, k, l, &
                                                  & eqn_idx%E) = (s_M*(cm%R(norm_dir) - Ga%R*alpha_rho_R(1)*vel_R(norm_dir)) &
                                                  & - s_P*(cm%L(norm_dir) - Ga%L*alpha_rho_L(1)*vel_L(norm_dir)) + s_M*s_P*(E_L &
                                                  & - E_R))/(s_M - s_P)
                            else if (bubbles_euler) then
                                flux_rs${XYZ}$_vf(j, k, l, &
                                                  & eqn_idx%E) = (s_M*vel_R(dir_idx(1))*(E_R + pres_R - ptilde_R) &
                                                  & - s_P*vel_L(dir_idx(1))*(E_L + pres_L - ptilde_L) + s_M*s_P*(E_L - E_R))/(s_M &
                                                  & - s_P) + (s_M/s_L)*(s_P/s_R)*pcorr*(vel_R_rms - vel_L_rms)/2._wp
                            else if (hypoelasticity) then
                                flux_tau_L = 0._wp; flux_tau_R = 0._wp
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    flux_tau_L = flux_tau_L + tau_e_L(dir_idx_tau(i))*vel_L(dir_idx(i))
                                    flux_tau_R = flux_tau_R + tau_e_R(dir_idx_tau(i))*vel_R(dir_idx(i))
                                end do
                                flux_rs${XYZ}$_vf(j, k, l, &
                                                  & eqn_idx%E) = (s_M*(vel_R(dir_idx(1))*(E_R + pres_R) - flux_tau_R) &
                                                  & - s_P*(vel_L(dir_idx(1))*(E_L + pres_L) - flux_tau_L) + s_M*s_P*(E_L - E_R)) &
                                                  & /(s_M - s_P)
                            else
                                flux_rs${XYZ}$_vf(j, k, l, &
                                                  & eqn_idx%E) = (s_M*vel_R(dir_idx(1))*(E_R + pres_R) - s_P*vel_L(dir_idx(1)) &
                                                  & *(E_L + pres_L) + s_M*s_P*(E_L - E_R))/(s_M - s_P) + (s_M/s_L)*(s_P/s_R) &
                                                  & *pcorr*(vel_R_rms - vel_L_rms)/2._wp
                            end if

                            ! Elastic Stresses
                            if (hypoelasticity) then
                                do i = 1, eqn_idx%stress%end - eqn_idx%stress%beg + 1  ! TODO: this indexing may be slow
                                    flux_rs${XYZ}$_vf(j, k, l, &
                                                      & eqn_idx%stress%beg - 1 + i) = (s_M*(rho_R*vel_R(dir_idx(1))*tau_e_R(i)) &
                                                      & - s_P*(rho_L*vel_L(dir_idx(1))*tau_e_L(i)) + s_M*s_P*(rho_L*tau_e_L(i) &
                                                      & - rho_R*tau_e_R(i)))/(s_M - s_P)
                                end do
                            end if

                            ! Advection flux and source: interface velocity for volume fraction transport
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = eqn_idx%adv%beg, eqn_idx%adv%end
                                flux_rs${XYZ}$_vf(j, k, l, i) = (qL_prim_rs${XYZ}$_vf(j, k, l, i) - qR_prim_rs${XYZ}$_vf(j + 1, &
                                                  & k, l, i))*s_M*s_P/(s_M - s_P)
                                flux_src_rs${XYZ}$_vf(j, k, l, i) = (s_M*qR_prim_rs${XYZ}$_vf(j + 1, k, l, &
                                                      & i) - s_P*qL_prim_rs${XYZ}$_vf(j, k, l, i))/(s_M - s_P)
                            end do

                            if (bubbles_euler) then
                                ! From HLLC: Kills mass transport @ bubble gas density
                                if (num_fluids > 1) then
                                    flux_rs${XYZ}$_vf(j, k, l, eqn_idx%cont%end) = 0._wp
                                end if
                            end if

                            if (chemistry) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = eqn_idx%species%beg, eqn_idx%species%end
                                    Y_L = qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                    Y_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)

                                    flux_rs${XYZ}$_vf(j, k, l, &
                                                      & i) = (s_M*Y_R*rho_R*vel_R(dir_idx(1)) - s_P*Y_L*rho_L*vel_L(dir_idx(1)) &
                                                      & + s_M*s_P*(Y_L*rho_L - Y_R*rho_R))/(s_M - s_P)
                                    flux_src_rs${XYZ}$_vf(j, k, l, i) = 0._wp
                                end do
                            end if

                            ! MHD: magnetic flux and Maxwell stress contributions
                            if (mhd) then
                                if (n == 0) then  ! 1D: d/dx flux only & Bx = Bx0 = const.
                                    ! B_y flux = v_x * B_y - v_y * Bx0 B_z flux = v_x * B_z - v_z * Bx0
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 0, 1
                                        flux_rsx_vf(j, k, l, &
                                                    & eqn_idx%B%beg + i) = (s_M*(vel_R(1)*B%R(2 + i) - vel_R(2 + i)*Bx0) &
                                                    & - s_P*(vel_L(1)*B%L(2 + i) - vel_L(2 + i)*Bx0) + s_M*s_P*(B%L(2 + i) &
                                                    & - B%R(2 + i)))/(s_M - s_P)
                                    end do
                                else  ! 2D/3D: Bx, By, Bz /= const. but zero flux component in the same direction
                                    ! B_x d/d${XYZ}$ flux = (1 - delta(x,${XYZ}$)) * (v_${XYZ}$ * B_x - v_x * B_${XYZ}$) B_y
                                    ! d/d${XYZ}$ flux = (1 - delta(y,${XYZ}$)) * (v_${XYZ}$ * B_y - v_y * B_${XYZ}$) B_z d/d${XYZ}$
                                    ! flux = (1 - delta(z,${XYZ}$)) * (v_${XYZ}$ * B_z - v_z * B_${XYZ}$)
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 0, 2
                                        flux_rs${XYZ}$_vf(j, k, l, &
                                                          & eqn_idx%B%beg + i) = (1 - dir_flg(i + 1))*(s_M*(vel_R(dir_idx(1)) &
                                                          & *B%R(i + 1) - vel_R(i + 1)*B%R(norm_dir)) - s_P*(vel_L(dir_idx(1)) &
                                                          & *B%L(i + 1) - vel_L(i + 1)*B%L(norm_dir)) + s_M*s_P*(B%L(i + 1) &
                                                          & - B%R(i + 1)))/(s_M - s_P)
                                    end do
                                end if
                                flux_src_rs${XYZ}$_vf(j, k, l, eqn_idx%adv%beg) = 0._wp
                            end if

                            #:if (NORM_DIR == 2)
                                if (cyl_coord) then
                                    ! Substituting the advective flux into the inviscid geometrical source flux
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, eqn_idx%E
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = flux_rs${XYZ}$_vf(j, k, l, i)
                                    end do
                                    ! Recalculating the radial momentum geometric source flux
                                    flux_gsrc_rs${XYZ}$_vf(j, k, l, eqn_idx%cont%end + 2) = flux_rs${XYZ}$_vf(j, k, l, &
                                                           & eqn_idx%cont%end + 2) - (s_M*pres_R - s_P*pres_L)/(s_M - s_P)
                                    ! Geometrical source of the void fraction(s) is zero
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = eqn_idx%adv%beg, eqn_idx%adv%end
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = flux_rs${XYZ}$_vf(j, k, l, i)
                                    end do
                                end if

                                if (cyl_coord .and. hypoelasticity) then
                                    ! += tau_sigmasigma using HLL
                                    flux_gsrc_rs${XYZ}$_vf(j, k, l, eqn_idx%cont%end + 2) = flux_gsrc_rs${XYZ}$_vf(j, k, l, &
                                                           & eqn_idx%cont%end + 2) + (s_M*tau_e_R(4) - s_P*tau_e_L(4))/(s_M - s_P)

                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = eqn_idx%stress%beg, eqn_idx%stress%end
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = flux_rs${XYZ}$_vf(j, k, l, i)
                                    end do
                                end if
                            #:endif
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

        if (viscous .or. dummy) then
            $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, idx_right_phys, vel_grad_L, vel_grad_R, alpha_L, alpha_R, &
                                & vel_L, vel_R, Re_L, Re_R]', copyin='[norm_dir]')
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
                                alpha_L(i) = qL_prim_rsy_vf(k, j, l, eqn_idx%E + i)
                                alpha_R(i) = qR_prim_rsy_vf(k + 1, j, l, eqn_idx%E + i)
                            end do
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_dims
                                vel_L(i) = qL_prim_rsy_vf(k, j, l, eqn_idx%mom%beg + i - 1)
                                vel_R(i) = qR_prim_rsy_vf(k + 1, j, l, eqn_idx%mom%beg + i - 1)
                            end do
                        else
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                alpha_L(i) = qL_prim_rsz_vf(l, k, j, eqn_idx%E + i)
                                alpha_R(i) = qR_prim_rsz_vf(l + 1, k, j, eqn_idx%E + i)
                            end do

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_dims
                                vel_L(i) = qL_prim_rsz_vf(l, k, j, eqn_idx%mom%beg + i - 1)
                                vel_R(i) = qR_prim_rsz_vf(l + 1, k, j, eqn_idx%mom%beg + i - 1)
                            end do
                        end if

                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, 2
                            Re_L(i) = dflt_real
                            Re_R(i) = dflt_real

                            if (Re_size(i) > 0) Re_L(i) = 0._wp
                            if (Re_size(i) > 0) Re_R(i) = 0._wp

                            $:GPU_LOOP(parallelism='[seq]')
                            do q = 1, Re_size(i)
                                Re_L(i) = alpha_L(Re_idx(i, q))/Res_gs(i, q) + Re_L(i)
                                Re_R(i) = alpha_R(Re_idx(i, q))/Res_gs(i, q) + Re_R(i)
                            end do

                            Re_L(i) = 1._wp/max(Re_L(i), sgm_eps)
                            Re_R(i) = 1._wp/max(Re_R(i), sgm_eps)
                        end do

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

    !> HLLC Riemann solver with contact restoration, Toro et al. Shock Waves (1994)
    subroutine s_hllc_riemann_solver(qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, dqL_prim_dx_vf, dqL_prim_dy_vf, &

        & dqL_prim_dz_vf, qL_prim_vf, qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, dqR_prim_dx_vf, dqR_prim_dy_vf, &
            & dqR_prim_dz_vf, qR_prim_vf, q_prim_vf, flux_vf, flux_src_vf, flux_gsrc_vf, norm_dir, ix, iy, iz)

        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout) :: qL_prim_rsx_vf, qL_prim_rsy_vf, &
             & qL_prim_rsz_vf, qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf
        type(scalar_field), dimension(sys_size), intent(in)          :: q_prim_vf
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
            real(wp), dimension(3) :: vel_L, vel_R
        #:else
            real(wp), dimension(num_fluids) :: alpha_rho_L, alpha_rho_R
            real(wp), dimension(num_fluids) :: alpha_L, alpha_R
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
        real(wp)               :: xi_L, xi_R  !< Left and right wave speeds functions
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
        real(wp) :: flux_tau_L, flux_tau_R
        real(wp) :: zcoef, pcorr           !< low Mach number correction
        integer  :: Re_max, i, j, k, l, q  !< Generic loop iterators

        ! HLLC star-state helpers
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(20) :: U_L, U_R, U_star_L, U_star_R
            real(wp), dimension(20) :: F_L, F_R, F_star_L, F_star_R, F_HLLC
        #:else
            real(wp), dimension(sys_size) :: U_L, U_R, U_star_L, U_star_R
            real(wp), dimension(sys_size) :: F_L, F_R, F_star_L, F_star_R, F_HLLC
        #:endif
        real(wp) :: u_n_HLLC, u_t_HLLC, u_t2_HLLC
        real(wp) :: pres_tot_L, pres_tot_R
        real(wp) :: u_n_L, u_n_R, u_t_L, u_t_R
        real(wp) :: u_t2_L, u_t2_R
        real(wp) :: tau_nn_L, tau_nn_R, tau_nt_L, tau_nt_R, tau_tt_L, tau_tt_R
        real(wp) :: tau_nt2_L, tau_nt2_R, tau_t2t2_L, tau_t2t2_R, tau_t1t2_L, tau_t1t2_R
        real(wp) :: tau_qq_L, tau_qq_R
        real(wp) :: p_face, tau_qq_face
        real(wp) :: A_L, A_R, denom_A, denom_L, denom_R
        real(wp) :: fac_L, fac_R
        real(wp) :: rho_L_star, rho_R_star
        real(wp) :: u_t_star, tau_nt_star
        real(wp) :: u_t2_star, tau_nt2_star
        real(wp) :: pres_tot_star
        real(wp) :: E_L_star, E_R_star
        real(wp) :: S_Mid
        integer  :: mom_n, mom_t1, mom_t2
        integer  :: idx_phys

        ! ADC (HLL -> HLLC)
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(20) :: F_HLL
        #:else
            real(wp), dimension(sys_size) :: F_HLL
        #:endif
        real(wp)            :: u_n_HLL_trace, u_t_HLL_trace
        real(wp)            :: u_t2_HLL_trace
        real(wp)            :: p_face_HLL, tau_qq_face_HLL, tau_nn_HLL
        real(wp)            :: phi
        real(wp)            :: Sigma_L, Sigma_R, dSigma, Sigma_ref
        real(wp)            :: a_L_ref, a_R_ref, a_ref
        real(wp)            :: du_t, dtau_nt
        real(wp)            :: du_t2, dtau_nt2
        real(wp)            :: sensor_ptot, sensor_vt, sensor_tnt, sensor_combined
        real(wp), parameter :: ADC_power = 1.0_wp

        ! Populating the buffers of the left and right Riemann problem states variables, based on the choice of boundary conditions

        call s_populate_riemann_states_variables_buffers(qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, dqL_prim_dx_vf, &
            & dqL_prim_dy_vf, dqL_prim_dz_vf, qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, dqR_prim_dx_vf, dqR_prim_dy_vf, &
            & dqR_prim_dz_vf, norm_dir, ix, iy, iz)

        ! Reshaping inputted data based on dimensional splitting direction

        call s_initialize_riemann_solver(flux_src_vf, norm_dir)

        #:for NORM_DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
            if (norm_dir == ${NORM_DIR}$) then
                ! 6-EQUATION MODEL WITH HLLC HLLC star-state flux with contact wave speed s_S
                if (model_eqns == 3) then
                    ! 6-equation model (model_eqns=3): separate phasic internal energies
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, q, vel_L, vel_R, Re_L, Re_R, alpha_L, alpha_R, Ys_L, &
                                        & Ys_R, Xs_L, Xs_R, Gamma_iL, Gamma_iR, Cp_iL, Cp_iR, Yi_avg, Phi_avg, h_iL, h_iR, &
                                        & h_avg_2, tau_e_L, tau_e_R, flux_ene_e, xi_field_L, xi_field_R, pcorr, zcoef, rho_L, &
                                        & rho_R, pres_L, pres_R, E_L, E_R, H_L, H_R, Cp_avg, Cv_avg, T_avg, eps, c_sum_Yi_Phi, &
                                        & T_L, T_R, Y_L, Y_R, MW_L, MW_R, R_gas_L, R_gas_R, Cp_L, Cp_R, Cv_L, Cv_R, Gamm_L, &
                                        & Gamm_R, gamma_L, gamma_R, pi_inf_L, pi_inf_R, qv_L, qv_R, qv_avg, c_L, c_R, G_L, G_R, &
                                        & rho_avg, H_avg, c_avg, gamma_avg, ptilde_L, ptilde_R, vel_L_rms, vel_R_rms, &
                                        & vel_avg_rms, vel_L_tmp, vel_R_tmp, Ms_L, Ms_R, pres_SL, pres_SR, alpha_L_sum, &
                                        & alpha_R_sum, rho_Star, E_Star, p_Star, p_K_Star, vel_K_star, s_L, s_R, s_M, s_P, s_S, &
                                        & xi_M, xi_P, xi_L, xi_R, xi_MP, xi_PP, flux_tau_L, flux_tau_R, alpha_rho_L, alpha_rho_R]')
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end
                                vel_L_rms = 0._wp; vel_R_rms = 0._wp
                                rho_L = 0._wp; rho_R = 0._wp
                                gamma_L = 0._wp; gamma_R = 0._wp
                                pi_inf_L = 0._wp; pi_inf_R = 0._wp
                                qv_L = 0._wp; qv_R = 0._wp
                                alpha_L_sum = 0._wp; alpha_R_sum = 0._wp

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%cont%end + i)
                                    vel_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%cont%end + i)
                                    vel_L_rms = vel_L_rms + vel_L(i)**2._wp
                                    vel_R_rms = vel_R_rms + vel_R(i)**2._wp
                                end do

                                pres_L = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%E)
                                pres_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%E)

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
                                        qL_prim_rs${XYZ}$_vf(j, k, l, i) = max(0._wp, qL_prim_rs${XYZ}$_vf(j, k, l, i))
                                        qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%E + i) = min(max(0._wp, qL_prim_rs${XYZ}$_vf(j, k, &
                                                             & l, eqn_idx%E + i)), 1._wp)
                                        alpha_L_sum = alpha_L_sum + qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%E + i)
                                    end do

                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        qR_prim_rs${XYZ}$_vf(j + 1, k, l, i) = max(0._wp, qR_prim_rs${XYZ}$_vf(j + 1, k, l, i))
                                        qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%E + i) = min(max(0._wp, &
                                                             & qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%E + i)), 1._wp)
                                        alpha_R_sum = alpha_R_sum + qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%E + i)
                                    end do

                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%E + i) = qL_prim_rs${XYZ}$_vf(j, k, l, &
                                                             & eqn_idx%E + i)/max(alpha_L_sum, sgm_eps)
                                        qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%E + i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, &
                                                             & eqn_idx%E + i)/max(alpha_R_sum, sgm_eps)
                                    end do
                                end if

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    rho_L = rho_L + qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                    gamma_L = gamma_L + qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%E + i)*gammas(i)
                                    pi_inf_L = pi_inf_L + qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%E + i)*pi_infs(i)
                                    qv_L = qv_L + qL_prim_rs${XYZ}$_vf(j, k, l, i)*qvs(i)

                                    rho_R = rho_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)
                                    gamma_R = gamma_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%E + i)*gammas(i)
                                    pi_inf_R = pi_inf_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%E + i)*pi_infs(i)
                                    qv_R = qv_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)*qvs(i)

                                    alpha_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%adv%beg + i - 1)
                                    alpha_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%adv%beg + i - 1)
                                end do

                                if (viscous) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, 2
                                        Re_L(i) = dflt_real
                                        Re_R(i) = dflt_real
                                        if (Re_size(i) > 0) Re_L(i) = 0._wp
                                        if (Re_size(i) > 0) Re_R(i) = 0._wp
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do q = 1, Re_size(i)
                                            Re_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%E + Re_idx(i, q))/Res_gs(i, q) + Re_L(i)
                                            Re_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%E + Re_idx(i, q))/Res_gs(i, &
                                                 & q) + Re_R(i)
                                        end do
                                        Re_L(i) = 1._wp/max(Re_L(i), sgm_eps)
                                        Re_R(i) = 1._wp/max(Re_R(i), sgm_eps)
                                    end do
                                end if

                                E_L = gamma_L*pres_L + pi_inf_L + 5.e-1_wp*rho_L*vel_L_rms + qv_L
                                E_R = gamma_R*pres_R + pi_inf_R + 5.e-1_wp*rho_R*vel_R_rms + qv_R

                                ! ENERGY ADJUSTMENTS FOR HYPOELASTIC ENERGY
                                if (hypoelasticity) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, eqn_idx%stress%end - eqn_idx%stress%beg + 1
                                        tau_e_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%stress%beg - 1 + i)
                                        tau_e_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%stress%beg - 1 + i)
                                    end do
                                    G_L = 0._wp; G_R = 0._wp
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        G_L = G_L + alpha_L(i)*Gs_rs(i)
                                        G_R = G_R + alpha_R(i)*Gs_rs(i)
                                    end do
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, eqn_idx%stress%end - eqn_idx%stress%beg + 1
                                        ! Elastic energy (guard skips when G near zero)
                                        if (.not. hypo_energy_guard .or. ((G_L > verysmall) .and. (G_R > verysmall))) then
                                            E_L = E_L + (tau_e_L(i)*tau_e_L(i))/max(4._wp*G_L, verysmall)
                                            E_R = E_R + (tau_e_R(i)*tau_e_R(i))/max(4._wp*G_R, verysmall)
                                            ! Additional terms in 2D and 3D
                                            if ((i == 2) .or. (i == 4) .or. (i == 5)) then
                                                E_L = E_L + (tau_e_L(i)*tau_e_L(i))/max(4._wp*G_L, verysmall)
                                                E_R = E_R + (tau_e_R(i)*tau_e_R(i))/max(4._wp*G_R, verysmall)
                                            end if
                                        end if
                                    end do
                                end if

                                ! Hyperelastic stress contribution: strain energy added to total energy
                                if (hyperelasticity) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        xi_field_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%xi%beg - 1 + i)
                                        xi_field_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%xi%beg - 1 + i)
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
                                        E_L = E_L + G_L*qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%xi%end + 1)
                                        E_R = E_R + G_R*qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%xi%end + 1)
                                    end if
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, b_size - 1
                                        tau_e_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%stress%beg - 1 + i)
                                        tau_e_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%stress%beg - 1 + i)
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
                                        Re_avg_rs${XYZ}$_vf(j, k, l, i) = 2._wp/(1._wp/Re_L(i) + 1._wp/Re_R(i))
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
                                        s_L = min(vel_L(dir_idx(1)) - sqrt(max(verysmall, &
                                                  & c_L*c_L + (((4._wp*G_L)/3._wp) + tau_e_L(dir_idx_tau(1)))/rho_L)), &
                                                  & vel_R(dir_idx(1)) - sqrt(max(verysmall, &
                                                  & c_R*c_R + (((4._wp*G_R)/3._wp) + tau_e_R(dir_idx_tau(1)))/rho_R)))
                                        s_R = max(vel_R(dir_idx(1)) + sqrt(max(verysmall, &
                                                  & c_R*c_R + (((4._wp*G_R)/3._wp) + tau_e_R(dir_idx_tau(1)))/rho_R)), &
                                                  & vel_L(dir_idx(1)) + sqrt(max(verysmall, &
                                                  & c_L*c_L + (((4._wp*G_L)/3._wp) + tau_e_L(dir_idx_tau(1)))/rho_L)))
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
                                else if (wave_speeds == 2) then
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
                                xi_L = (s_L - vel_L(dir_idx(1)))/(s_L - s_S)
                                xi_R = (s_R - vel_R(dir_idx(1)))/(s_R - s_S)

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
                                    flux_rs${XYZ}$_vf(j, k, l, i) = xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, &
                                                      & i)*(vel_L(dir_idx(1)) + s_M*(xi_L - 1._wp)) + xi_P*qR_prim_rs${XYZ}$_vf(j &
                                                      & + 1, k, l, i)*(vel_R(dir_idx(1)) + s_P*(xi_R - 1._wp))
                                end do

                                ! MOMENTUM FLUX. f = \rho u u - \sigma, q = \rho u, q_star = \xi * \rho*(s_star, v, w)
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    flux_rs${XYZ}$_vf(j, k, l, &
                                                      & eqn_idx%cont%end + dir_idx(i)) = rho_Star*vel_K_Star*(dir_flg(dir_idx(i)) &
                                                      & *vel_K_Star + (1._wp - dir_flg(dir_idx(i)))*(xi_M*vel_L(dir_idx(i)) &
                                                      & + xi_P*vel_R(dir_idx(i)))) + dir_flg(dir_idx(i))*p_Star + (s_M/s_L) &
                                                      & *(s_P/s_R)*dir_flg(dir_idx(i))*pcorr
                                end do

                                ! ENERGY FLUX. f = u*(E-\sigma), q = E, q_star = \xi*E+(s-u)(\rho s_star - \sigma/(s-u))
                                flux_rs${XYZ}$_vf(j, k, l, eqn_idx%E) = (E_star + p_Star)*vel_K_Star + (s_M/s_L)*(s_P/s_R)*pcorr*s_S

                                ! ELASTICITY. Elastic shear stress additions for the momentum and energy flux
                                if (elasticity) then
                                    flux_ene_e = 0._wp
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        ! MOMENTUM ELASTIC FLUX.
                                        flux_rs${XYZ}$_vf(j, k, l, eqn_idx%cont%end + dir_idx(i)) = flux_rs${XYZ}$_vf(j, k, l, &
                                                          & eqn_idx%cont%end + dir_idx(i)) - xi_M*tau_e_L(dir_idx_tau(i)) &
                                                          & - xi_P*tau_e_R(dir_idx_tau(i))
                                        ! ENERGY ELASTIC FLUX.
                                        flux_ene_e = flux_ene_e - xi_M*(vel_L(dir_idx(i))*tau_e_L(dir_idx_tau(i)) &
                                                                        & + s_M*(xi_L*((s_S - vel_L(i))*(tau_e_L(dir_idx_tau(i)) &
                                                                        & /(s_L - vel_L(i)))))) - xi_P*(vel_R(dir_idx(i)) &
                                                                        & *tau_e_R(dir_idx_tau(i)) + s_P*(xi_R*((s_S - vel_R(i)) &
                                                                        & *(tau_e_R(dir_idx_tau(i))/(s_R - vel_R(i))))))
                                    end do
                                    flux_rs${XYZ}$_vf(j, k, l, eqn_idx%E) = flux_rs${XYZ}$_vf(j, k, l, eqn_idx%E) + flux_ene_e
                                end if

                                ! VOLUME FRACTION FLUX.
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = eqn_idx%adv%beg, eqn_idx%adv%end
                                    flux_rs${XYZ}$_vf(j, k, l, i) = xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, &
                                                      & i)*s_S + xi_P*qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)*s_S
                                end do

                                ! Advection velocity source: interface velocity for volume fraction transport
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel_src_rs${XYZ}$_vf(j, k, l, &
                                                         & dir_idx(i)) = xi_M*(vel_L(dir_idx(i)) + dir_flg(dir_idx(i)) &
                                                         & *(s_S*(xi_MP*(xi_L - 1) + 1) - vel_L(dir_idx(i)))) &
                                                         & + xi_P*(vel_R(dir_idx(i)) + dir_flg(dir_idx(i))*(s_S*(xi_PP*(xi_R - 1) &
                                                         & + 1) - vel_R(dir_idx(i))))
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

                                    flux_rs${XYZ}$_vf(j, k, l, i + eqn_idx%int_en%beg - 1) = ((xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, &
                                                      & i + eqn_idx%adv%beg - 1) + xi_P*qR_prim_rs${XYZ}$_vf(j + 1, k, l, &
                                                      & i + eqn_idx%adv%beg - 1))*(gammas(i)*p_K_Star + pi_infs(i)) &
                                                      & + (xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, &
                                                      & i + eqn_idx%cont%beg - 1) + xi_P*qR_prim_rs${XYZ}$_vf(j + 1, k, l, &
                                                      & i + eqn_idx%cont%beg - 1))*qvs(i))*vel_K_Star + (s_M/s_L)*(s_P/s_R) &
                                                      & *pcorr*s_S*(xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, &
                                                      & i + eqn_idx%adv%beg - 1) + xi_P*qR_prim_rs${XYZ}$_vf(j + 1, k, l, &
                                                      & i + eqn_idx%adv%beg - 1))
                                end do

                                flux_src_rs${XYZ}$_vf(j, k, l, eqn_idx%adv%beg) = vel_src_rs${XYZ}$_vf(j, k, l, dir_idx(1))

                                ! HYPOELASTIC STRESS EVOLUTION FLUX.
                                if (hypoelasticity) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, eqn_idx%stress%end - eqn_idx%stress%beg + 1
                                        flux_rs${XYZ}$_vf(j, k, l, &
                                                          & eqn_idx%stress%beg - 1 + i) = xi_M*(s_S/(s_L - s_S)) &
                                                          & *(s_L*rho_L*tau_e_L(i) - rho_L*vel_L(dir_idx(1))*tau_e_L(i)) &
                                                          & + xi_P*(s_S/(s_R - s_S))*(s_R*rho_R*tau_e_R(i) &
                                                          & - rho_R*vel_R(dir_idx(1))*tau_e_R(i))
                                    end do
                                end if

                                ! Hyperelastic reference map flux for material deformation tracking
                                if (hyperelasticity) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        flux_rs${XYZ}$_vf(j, k, l, &
                                                          & eqn_idx%xi%beg - 1 + i) = xi_M*(s_S/(s_L - s_S)) &
                                                          & *(s_L*rho_L*xi_field_L(i) - rho_L*vel_L(dir_idx(1))*xi_field_L(i)) &
                                                          & + xi_P*(s_S/(s_R - s_S))*(s_R*rho_R*xi_field_R(i) &
                                                          & - rho_R*vel_R(dir_idx(1))*xi_field_R(i))
                                    end do
                                end if

                                ! COLOR FUNCTION FLUX
                                if (surface_tension) then
                                    flux_rs${XYZ}$_vf(j, k, l, eqn_idx%c) = (xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, &
                                                      & eqn_idx%c) + xi_P*qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%c))*s_S
                                end if

                                ! Geometrical source flux for cylindrical coordinates
                                #:if (NORM_DIR == 2)
                                    if (cyl_coord) then
                                        ! Substituting the advective flux into the inviscid geometrical source flux
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, eqn_idx%E
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = flux_rs${XYZ}$_vf(j, k, l, i)
                                        end do
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = eqn_idx%int_en%beg, eqn_idx%int_en%end
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = flux_rs${XYZ}$_vf(j, k, l, i)
                                        end do
                                        ! Recalculating the radial momentum geometric source flux
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, &
                                                               & eqn_idx%mom%beg - 1 + dir_idx(1)) = flux_gsrc_rs${XYZ}$_vf(j, k, &
                                                               & l, eqn_idx%mom%beg - 1 + dir_idx(1)) - p_Star
                                        ! Geometrical source of the void fraction(s) is zero
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = eqn_idx%adv%beg, eqn_idx%adv%end
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = 0._wp
                                        end do
                                    end if
                                #:endif
                                #:if (NORM_DIR == 3)
                                    if (grid_geometry == 3) then
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, sys_size
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = 0._wp
                                        end do
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, &
                                                               & eqn_idx%mom%beg - 1 + dir_idx(1)) = flux_gsrc_rs${XYZ}$_vf(j, k, &
                                                               & l, eqn_idx%mom%beg - 1 + dir_idx(1)) - p_Star

                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, eqn_idx%mom%end) = flux_rs${XYZ}$_vf(j, k, l, &
                                                               & eqn_idx%mom%beg + 1)
                                    end if
                                #:endif
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                else if (model_eqns == 4) then
                    ! 4-equation model (model_eqns=4): single pressure, velocity equilibrium
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[i, q, alpha_rho_L, alpha_rho_R, vel_L, vel_R, alpha_L, alpha_R, &
                                        & nbub_L, nbub_R, rho_L, rho_R, pres_L, pres_R, E_L, E_R, H_L, H_R, Cp_avg, Cv_avg, &
                                        & T_avg, eps, c_sum_Yi_Phi, T_L, T_R, Y_L, Y_R, MW_L, MW_R, R_gas_L, R_gas_R, Cp_L, Cp_R, &
                                        & Gamm_L, Gamm_R, gamma_L, gamma_R, pi_inf_L, pi_inf_R, qv_L, qv_R, qv_avg, c_L, c_R, &
                                        & G_L, G_R, rho_avg, H_avg, c_avg, gamma_avg, ptilde_L, ptilde_R, vel_L_rms, vel_R_rms, &
                                        & vel_avg_rms, vel_L_tmp, vel_R_tmp, Ms_L, Ms_R, pres_SL, pres_SR, alpha_L_sum, &
                                        & alpha_R_sum, rho_Star, E_Star, p_Star, p_K_Star, vel_K_star, s_L, s_R, s_M, s_P, s_S, &
                                        & xi_M, xi_P, xi_L, xi_R, xi_MP, xi_PP, Ys_L, Ys_R, Cp_iL, Cp_iR, Xs_L, Xs_R, Gamma_iL, &
                                        & Gamma_iR, Yi_avg, Phi_avg, h_iL, h_iR, h_avg_2]')
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end
                                vel_L_rms = 0._wp; vel_R_rms = 0._wp
                                rho_L = 0._wp; rho_R = 0._wp
                                gamma_L = 0._wp; gamma_R = 0._wp
                                pi_inf_L = 0._wp; pi_inf_R = 0._wp
                                qv_L = 0._wp; qv_R = 0._wp

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, eqn_idx%cont%end
                                    alpha_rho_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                    alpha_rho_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)
                                end do

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%cont%end + i)
                                    vel_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%cont%end + i)
                                    vel_L_rms = vel_L_rms + vel_L(i)**2._wp
                                    vel_R_rms = vel_R_rms + vel_R(i)**2._wp
                                end do

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    alpha_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%E + i)
                                    alpha_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%E + i)
                                end do
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    alpha_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%E + i)
                                    alpha_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%E + i)
                                end do

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    rho_L = rho_L + alpha_rho_L(i)
                                    gamma_L = gamma_L + alpha_L(i)*gammas(i)
                                    pi_inf_L = pi_inf_L + alpha_L(i)*pi_infs(i)
                                    qv_L = qv_L + alpha_rho_L(i)*qvs(i)

                                    rho_R = rho_R + alpha_rho_R(i)
                                    gamma_R = gamma_R + alpha_R(i)*gammas(i)
                                    pi_inf_R = pi_inf_R + alpha_R(i)*pi_infs(i)
                                    qv_R = qv_R + alpha_rho_R(i)*qvs(i)
                                end do

                                pres_L = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%E)
                                pres_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%E)

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

                                if (wave_speeds == 1) then
                                    s_L = min(vel_L(dir_idx(1)) - c_L, vel_R(dir_idx(1)) - c_R)
                                    s_R = max(vel_R(dir_idx(1)) + c_R, vel_L(dir_idx(1)) + c_L)

                                    s_S = (pres_R - pres_L + rho_L*vel_L(dir_idx(1))*(s_L - vel_L(dir_idx(1))) &
                                           & - rho_R*vel_R(dir_idx(1))*(s_R - vel_R(dir_idx(1))))/(rho_L*(s_L - vel_L(dir_idx(1))) &
                                           & - rho_R*(s_R - vel_R(dir_idx(1))))
                                else if (wave_speeds == 2) then
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
                                xi_L = (s_L - vel_L(dir_idx(1)))/(s_L - s_S)
                                xi_R = (s_R - vel_R(dir_idx(1)))/(s_R - s_S)

                                ! goes with numerical velocity in x/y/z directions xi_P/M = 0.5 +/m sgn(0.5,s_star)
                                xi_M = (5.e-1_wp + sign(5.e-1_wp, s_S))
                                xi_P = (5.e-1_wp - sign(5.e-1_wp, s_S))

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, eqn_idx%cont%end
                                    flux_rs${XYZ}$_vf(j, k, l, &
                                                      & i) = xi_M*alpha_rho_L(i)*(vel_L(dir_idx(1)) + s_M*(xi_L - 1._wp)) &
                                                      & + xi_P*alpha_rho_R(i)*(vel_R(dir_idx(1)) + s_P*(xi_R - 1._wp))
                                end do

                                ! Momentum flux. f = \rho u u + p I, q = \rho u, q_star = \xi * \rho*(s_star, v, w)
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    flux_rs${XYZ}$_vf(j, k, l, &
                                                      & eqn_idx%cont%end + dir_idx(i)) = xi_M*(rho_L*(vel_L(dir_idx(1)) &
                                                      & *vel_L(dir_idx(i)) + s_M*(xi_L*(dir_flg(dir_idx(i))*s_S + (1._wp &
                                                      & - dir_flg(dir_idx(i)))*vel_L(dir_idx(i))) - vel_L(dir_idx(i)))) &
                                                      & + dir_flg(dir_idx(i))*pres_L) + xi_P*(rho_R*(vel_R(dir_idx(1)) &
                                                      & *vel_R(dir_idx(i)) + s_P*(xi_R*(dir_flg(dir_idx(i))*s_S + (1._wp &
                                                      & - dir_flg(dir_idx(i)))*vel_R(dir_idx(i))) - vel_R(dir_idx(i)))) &
                                                      & + dir_flg(dir_idx(i))*pres_R)
                                end do

                                if (bubbles_euler) then
                                    ! Put p_tilde in
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        flux_rs${XYZ}$_vf(j, k, l, eqn_idx%cont%end + dir_idx(i)) = flux_rs${XYZ}$_vf(j, k, l, &
                                                          & eqn_idx%cont%end + dir_idx(i)) + xi_M*(dir_flg(dir_idx(i)) &
                                                          & *(-1._wp*ptilde_L)) + xi_P*(dir_flg(dir_idx(i))*(-1._wp*ptilde_R))
                                    end do
                                end if

                                flux_rs${XYZ}$_vf(j, k, l, eqn_idx%E) = 0._wp

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = eqn_idx%alf, eqn_idx%alf  ! only advect the void fraction
                                    flux_rs${XYZ}$_vf(j, k, l, i) = xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, &
                                                      & i)*(vel_L(dir_idx(1)) + s_M*(xi_L - 1._wp)) + xi_P*qR_prim_rs${XYZ}$_vf(j &
                                                      & + 1, k, l, i)*(vel_R(dir_idx(1)) + s_P*(xi_R - 1._wp))
                                end do

                                ! Advection velocity source: interface velocity for volume fraction transport
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel_src_rs${XYZ}$_vf(j, k, l, dir_idx(i)) = 0._wp
                                    ! IF ( (model_eqns == 4) .or. (num_fluids==1) ) vel_src_rs_vf(dir_idx(i))%sf(j,k,l) = 0._wp
                                end do

                                flux_src_rs${XYZ}$_vf(j, k, l, eqn_idx%adv%beg) = vel_src_rs${XYZ}$_vf(j, k, l, dir_idx(1))

                                ! Add advection flux for bubble variables
                                if (bubbles_euler) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = eqn_idx%bub%beg, eqn_idx%bub%end
                                        flux_rs${XYZ}$_vf(j, k, l, i) = xi_M*nbub_L*qL_prim_rs${XYZ}$_vf(j, k, l, &
                                                          & i)*(vel_L(dir_idx(1)) + s_M*(xi_L - 1._wp)) &
                                                          & + xi_P*nbub_R*qR_prim_rs${XYZ}$_vf(j + 1, k, l, &
                                                          & i)*(vel_R(dir_idx(1)) + s_P*(xi_R - 1._wp))
                                    end do
                                end if

                                ! Geometrical source flux for cylindrical coordinates

                                #:if (NORM_DIR == 2)
                                    if (cyl_coord) then
                                        ! Substituting the advective flux into the inviscid geometrical source flux
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, eqn_idx%E
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = flux_rs${XYZ}$_vf(j, k, l, i)
                                        end do
                                        ! Recalculating the radial momentum geometric source flux
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, &
                                                               & eqn_idx%cont%end + dir_idx(1)) = xi_M*(rho_L*(vel_L(dir_idx(1)) &
                                                               & *vel_L(dir_idx(1)) + s_M*(xi_L*(dir_flg(dir_idx(1))*s_S + (1._wp &
                                                               & - dir_flg(dir_idx(1)))*vel_L(dir_idx(1))) - vel_L(dir_idx(1))))) &
                                                               & + xi_P*(rho_R*(vel_R(dir_idx(1))*vel_R(dir_idx(1)) &
                                                               & + s_P*(xi_R*(dir_flg(dir_idx(1))*s_S + (1._wp &
                                                               & - dir_flg(dir_idx(1)))*vel_R(dir_idx(1))) - vel_R(dir_idx(1)))))
                                        ! Geometrical source of the void fraction(s) is zero
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = eqn_idx%adv%beg, eqn_idx%adv%end
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = 0._wp
                                        end do
                                    end if
                                #:endif
                                #:if (NORM_DIR == 3)
                                    if (grid_geometry == 3) then
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, sys_size
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = 0._wp
                                        end do
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, &
                                                               & eqn_idx%mom%beg + 1) = -xi_M*(rho_L*(vel_L(dir_idx(1)) &
                                                               & *vel_L(dir_idx(1)) + s_M*(xi_L*(dir_flg(dir_idx(1))*s_S + (1._wp &
                                                               & - dir_flg(dir_idx(1)))*vel_L(dir_idx(1))) - vel_L(dir_idx(1))))) &
                                                               & - xi_P*(rho_R*(vel_R(dir_idx(1))*vel_R(dir_idx(1)) &
                                                               & + s_P*(xi_R*(dir_flg(dir_idx(1))*s_S + (1._wp &
                                                               & - dir_flg(dir_idx(1)))*vel_R(dir_idx(1))) - vel_R(dir_idx(1)))))
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, eqn_idx%mom%end) = flux_rs${XYZ}$_vf(j, k, l, &
                                                               & eqn_idx%mom%beg + 1)
                                    end if
                                #:endif
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                else if (model_eqns == 2 .and. bubbles_euler) then
                    ! 5-equation model with Euler-Euler bubble dynamics
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[i, q, R0_L, R0_R, V0_L, V0_R, P0_L, P0_R, pbw_L, pbw_R, vel_L, &
                                        & vel_R, rho_avg, alpha_L, alpha_R, h_avg, gamma_avg, Re_L, Re_R, pcorr, zcoef, rho_L, &
                                        & rho_R, pres_L, pres_R, E_L, E_R, H_L, H_R, gamma_L, gamma_R, pi_inf_L, pi_inf_R, qv_L, &
                                        & qv_R, qv_avg, c_L, c_R, c_avg, vel_L_rms, vel_R_rms, vel_avg_rms, vel_L_tmp, vel_R_tmp, &
                                        & Ms_L, Ms_R, pres_SL, pres_SR, alpha_L_sum, alpha_R_sum, s_L, s_R, s_M, s_P, s_S, xi_M, &
                                        & xi_P, xi_L, xi_R, xi_MP, xi_PP, nbub_L, nbub_R, PbwR3Lbar, PbwR3Rbar, R3Lbar, R3Rbar, &
                                        & R3V2Lbar, R3V2Rbar, Ys_L, Ys_R, Cp_iL, Cp_iR, Xs_L, Xs_R, Gamma_iL, Gamma_iR, Yi_avg, &
                                        & Phi_avg, h_iL, h_iR, h_avg_2]')
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end
                                vel_L_rms = 0._wp; vel_R_rms = 0._wp
                                rho_L = 0._wp; rho_R = 0._wp
                                gamma_L = 0._wp; gamma_R = 0._wp
                                pi_inf_L = 0._wp; pi_inf_R = 0._wp
                                qv_L = 0._wp; qv_R = 0._wp

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    alpha_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%E + i)
                                    alpha_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%E + i)
                                end do

                                vel_L_rms = 0._wp; vel_R_rms = 0._wp

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%cont%end + i)
                                    vel_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%cont%end + i)
                                    vel_L_rms = vel_L_rms + vel_L(i)**2._wp
                                    vel_R_rms = vel_R_rms + vel_R(i)**2._wp
                                end do

                                ! Retain this in the refactor
                                if (mpp_lim .and. (num_fluids > 2)) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        rho_L = rho_L + qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                        gamma_L = gamma_L + qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%E + i)*gammas(i)
                                        pi_inf_L = pi_inf_L + qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%E + i)*pi_infs(i)
                                        qv_L = qv_L + qL_prim_rs${XYZ}$_vf(j, k, l, i)*qvs(i)
                                        rho_R = rho_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)
                                        gamma_R = gamma_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%E + i)*gammas(i)
                                        pi_inf_R = pi_inf_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%E + i)*pi_infs(i)
                                        qv_R = qv_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)*qvs(i)
                                    end do
                                else if (num_fluids > 2) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids - 1
                                        rho_L = rho_L + qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                        gamma_L = gamma_L + qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%E + i)*gammas(i)
                                        pi_inf_L = pi_inf_L + qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%E + i)*pi_infs(i)
                                        qv_L = qv_L + qL_prim_rs${XYZ}$_vf(j, k, l, i)*qvs(i)
                                        rho_R = rho_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)
                                        gamma_R = gamma_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%E + i)*gammas(i)
                                        pi_inf_R = pi_inf_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%E + i)*pi_infs(i)
                                        qv_R = qv_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)*qvs(i)
                                    end do
                                else
                                    rho_L = qL_prim_rs${XYZ}$_vf(j, k, l, 1)
                                    gamma_L = gammas(1)
                                    pi_inf_L = pi_infs(1)
                                    qv_L = qvs(1)
                                    rho_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, 1)
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
                                                Re_L(i) = (1._wp - qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%E + Re_idx(i, &
                                                     & q)))/Res_gs(i, q) + Re_L(i)
                                                Re_R(i) = (1._wp - qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%E + Re_idx(i, &
                                                     & q)))/Res_gs(i, q) + Re_R(i)
                                            end do

                                            Re_L(i) = 1._wp/max(Re_L(i), sgm_eps)
                                            Re_R(i) = 1._wp/max(Re_R(i), sgm_eps)
                                        end do
                                    end if
                                end if

                                pres_L = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%E)
                                pres_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%E)

                                E_L = gamma_L*pres_L + pi_inf_L + 5.e-1_wp*rho_L*vel_L_rms
                                E_R = gamma_R*pres_R + pi_inf_R + 5.e-1_wp*rho_R*vel_R_rms

                                H_L = (E_L + pres_L)/rho_L
                                H_R = (E_R + pres_R)/rho_R

                                if (avg_state == 2) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, nb
                                        R0_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, rs(i))
                                        R0_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, rs(i))

                                        V0_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, vs(i))
                                        V0_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, vs(i))
                                        if (.not. polytropic .and. .not. qbmm) then
                                            P0_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, ps(i))
                                            P0_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, ps(i))
                                        end if
                                    end do

                                    if (.not. qbmm) then
                                        if (adv_n) then
                                            nbub_L = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%n)
                                            nbub_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%n)
                                        else
                                            nbub_L = 0._wp
                                            nbub_R = 0._wp
                                            $:GPU_LOOP(parallelism='[seq]')
                                            do i = 1, nb
                                                nbub_L = nbub_L + (R0_L(i)**3._wp)*weight(i)
                                                nbub_R = nbub_R + (R0_R(i)**3._wp)*weight(i)
                                            end do

                                            nbub_L = (3._wp/(4._wp*pi))*qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%E + num_fluids)/nbub_L
                                            nbub_R = (3._wp/(4._wp*pi))*qR_prim_rs${XYZ}$_vf(j + 1, k, l, &
                                                      & eqn_idx%E + num_fluids)/nbub_R
                                        end if
                                    else
                                        ! nb stored in 0th moment of first R0 bin in variable conversion module
                                        nbub_L = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%bub%beg)
                                        nbub_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%bub%beg)
                                    end if

                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, nb
                                        if (.not. qbmm) then
                                            pbw_L(i) = f_cpbw_KM(R0(i), R0_L(i), V0_L(i), P0_L(i))
                                            pbw_R(i) = f_cpbw_KM(R0(i), R0_R(i), V0_R(i), P0_R(i))
                                        end if
                                    end do

                                    if (qbmm) then
                                        PbwR3Lbar = mom_sp_rs${XYZ}$_vf(j, k, l, 4)
                                        PbwR3Rbar = mom_sp_rs${XYZ}$_vf(j + 1, k, l, 4)

                                        R3Lbar = mom_sp_rs${XYZ}$_vf(j, k, l, 1)
                                        R3Rbar = mom_sp_rs${XYZ}$_vf(j + 1, k, l, 1)

                                        R3V2Lbar = mom_sp_rs${XYZ}$_vf(j, k, l, 3)
                                        R3V2Rbar = mom_sp_rs${XYZ}$_vf(j + 1, k, l, 3)
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
                                        Re_avg_rs${XYZ}$_vf(j, k, l, i) = 2._wp/(1._wp/Re_L(i) + 1._wp/Re_R(i))
                                    end do
                                end if

                                ! Low Mach correction
                                if (low_Mach == 2) then
                                    @:compute_low_Mach_correction()
                                end if

                                if (wave_speeds == 1) then
                                    s_L = min(vel_L(dir_idx(1)) - c_L, vel_R(dir_idx(1)) - c_R)
                                    s_R = max(vel_R(dir_idx(1)) + c_R, vel_L(dir_idx(1)) + c_L)

                                    s_S = (pres_R - pres_L + rho_L*vel_L(dir_idx(1))*(s_L - vel_L(dir_idx(1))) &
                                           & - rho_R*vel_R(dir_idx(1))*(s_R - vel_R(dir_idx(1))))/(rho_L*(s_L - vel_L(dir_idx(1))) &
                                           & - rho_R*(s_R - vel_R(dir_idx(1))))
                                else if (wave_speeds == 2) then
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
                                xi_L = (s_L - vel_L(dir_idx(1)))/(s_L - s_S)
                                xi_R = (s_R - vel_R(dir_idx(1)))/(s_R - s_S)

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
                                    flux_rs${XYZ}$_vf(j, k, l, i) = xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, &
                                                      & i)*(vel_L(dir_idx(1)) + s_M*(xi_L - 1._wp)) + xi_P*qR_prim_rs${XYZ}$_vf(j &
                                                      & + 1, k, l, i)*(vel_R(dir_idx(1)) + s_P*(xi_R - 1._wp))
                                end do

                                if (bubbles_euler .and. (num_fluids > 1)) then
                                    ! Kill mass transport @ gas density
                                    flux_rs${XYZ}$_vf(j, k, l, eqn_idx%cont%end) = 0._wp
                                end if

                                ! Momentum flux. f = \rho u u + p I, q = \rho u, q_star = \xi * \rho*(s_star, v, w)

                                ! Include p_tilde

                                if (avg_state == 2) then
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
                                    flux_rs${XYZ}$_vf(j, k, l, &
                                                      & eqn_idx%cont%end + dir_idx(i)) = xi_M*(rho_L*(vel_L(dir_idx(1)) &
                                                      & *vel_L(dir_idx(i)) + s_M*(xi_L*(dir_flg(dir_idx(i))*s_S + (1._wp &
                                                      & - dir_flg(dir_idx(i)))*vel_L(dir_idx(i))) - vel_L(dir_idx(i)))) &
                                                      & + dir_flg(dir_idx(i))*(pres_L)) + xi_P*(rho_R*(vel_R(dir_idx(1)) &
                                                      & *vel_R(dir_idx(i)) + s_P*(xi_R*(dir_flg(dir_idx(i))*s_S + (1._wp &
                                                      & - dir_flg(dir_idx(i)))*vel_R(dir_idx(i))) - vel_R(dir_idx(i)))) &
                                                      & + dir_flg(dir_idx(i))*(pres_R)) + (s_M/s_L)*(s_P/s_R)*dir_flg(dir_idx(i)) &
                                                      & *pcorr
                                end do

                                ! Energy flux. f = u*(E+p), q = E, q_star = \xi*E+(s-u)(\rho s_star + p/(s-u))
                                flux_rs${XYZ}$_vf(j, k, l, &
                                                  & eqn_idx%E) = xi_M*(vel_L(dir_idx(1))*(E_L + pres_L) + s_M*(xi_L*(E_L + (s_S &
                                                  & - vel_L(dir_idx(1)))*(rho_L*s_S + (pres_L)/(s_L - vel_L(dir_idx(1))))) - E_L)) &
                                                  & + xi_P*(vel_R(dir_idx(1))*(E_R + pres_R) + s_P*(xi_R*(E_R + (s_S &
                                                  & - vel_R(dir_idx(1)))*(rho_R*s_S + (pres_R)/(s_R - vel_R(dir_idx(1))))) - E_R)) &
                                                  & + (s_M/s_L)*(s_P/s_R)*pcorr*s_S

                                ! Volume fraction flux
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = eqn_idx%adv%beg, eqn_idx%adv%end
                                    flux_rs${XYZ}$_vf(j, k, l, i) = xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, &
                                                      & i)*(vel_L(dir_idx(1)) + s_M*(xi_L - 1._wp)) + xi_P*qR_prim_rs${XYZ}$_vf(j &
                                                      & + 1, k, l, i)*(vel_R(dir_idx(1)) + s_P*(xi_R - 1._wp))
                                end do

                                ! Advection velocity source: interface velocity for volume fraction transport
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel_src_rs${XYZ}$_vf(j, k, l, &
                                                         & dir_idx(i)) = xi_M*(vel_L(dir_idx(i)) + dir_flg(dir_idx(i))*s_M*(xi_L &
                                                         & - 1._wp)) + xi_P*(vel_R(dir_idx(i)) + dir_flg(dir_idx(i))*s_P*(xi_R &
                                                         & - 1._wp))

                                    ! IF ( (model_eqns == 4) .or. (num_fluids==1) ) vel_src_rs_vf(dir_idx(i))%sf(j,k,l) = 0._wp
                                end do

                                flux_src_rs${XYZ}$_vf(j, k, l, eqn_idx%adv%beg) = vel_src_rs${XYZ}$_vf(j, k, l, dir_idx(1))

                                ! Add advection flux for bubble variables
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = eqn_idx%bub%beg, eqn_idx%bub%end
                                    flux_rs${XYZ}$_vf(j, k, l, i) = xi_M*nbub_L*qL_prim_rs${XYZ}$_vf(j, k, l, &
                                                      & i)*(vel_L(dir_idx(1)) + s_M*(xi_L - 1._wp)) &
                                                      & + xi_P*nbub_R*qR_prim_rs${XYZ}$_vf(j + 1, k, l, &
                                                      & i)*(vel_R(dir_idx(1)) + s_P*(xi_R - 1._wp))
                                end do

                                if (qbmm) then
                                    flux_rs${XYZ}$_vf(j, k, l, &
                                                      & eqn_idx%bub%beg) = xi_M*nbub_L*(vel_L(dir_idx(1)) + s_M*(xi_L - 1._wp)) &
                                                      & + xi_P*nbub_R*(vel_R(dir_idx(1)) + s_P*(xi_R - 1._wp))
                                end if

                                if (adv_n) then
                                    flux_rs${XYZ}$_vf(j, k, l, &
                                                      & eqn_idx%n) = xi_M*nbub_L*(vel_L(dir_idx(1)) + s_M*(xi_L - 1._wp)) &
                                                      & + xi_P*nbub_R*(vel_R(dir_idx(1)) + s_P*(xi_R - 1._wp))
                                end if

                                ! Geometrical source flux for cylindrical coordinates
                                #:if (NORM_DIR == 2)
                                    if (cyl_coord) then
                                        ! Substituting the advective flux into the inviscid geometrical source flux
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, eqn_idx%E
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = flux_rs${XYZ}$_vf(j, k, l, i)
                                        end do
                                        ! Recalculating the radial momentum geometric source flux
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, &
                                                               & eqn_idx%cont%end + dir_idx(1)) = xi_M*(rho_L*(vel_L(dir_idx(1)) &
                                                               & *vel_L(dir_idx(1)) + s_M*(xi_L*(dir_flg(dir_idx(1))*s_S + (1._wp &
                                                               & - dir_flg(dir_idx(1)))*vel_L(dir_idx(1))) - vel_L(dir_idx(1))))) &
                                                               & + xi_P*(rho_R*(vel_R(dir_idx(1))*vel_R(dir_idx(1)) &
                                                               & + s_P*(xi_R*(dir_flg(dir_idx(1))*s_S + (1._wp &
                                                               & - dir_flg(dir_idx(1)))*vel_R(dir_idx(1))) - vel_R(dir_idx(1)))))
                                        ! Geometrical source of the void fraction(s) is zero
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = eqn_idx%adv%beg, eqn_idx%adv%end
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = 0._wp
                                        end do
                                    end if
                                #:endif
                                #:if (NORM_DIR == 3)
                                    if (grid_geometry == 3) then
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, sys_size
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = 0._wp
                                        end do

                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, &
                                                               & eqn_idx%mom%beg + 1) = -xi_M*(rho_L*(vel_L(dir_idx(1)) &
                                                               & *vel_L(dir_idx(1)) + s_M*(xi_L*(dir_flg(dir_idx(1))*s_S + (1._wp &
                                                               & - dir_flg(dir_idx(1)))*vel_L(dir_idx(1))) - vel_L(dir_idx(1))))) &
                                                               & - xi_P*(rho_R*(vel_R(dir_idx(1))*vel_R(dir_idx(1)) &
                                                               & + s_P*(xi_R*(dir_flg(dir_idx(1))*s_S + (1._wp &
                                                               & - dir_flg(dir_idx(1)))*vel_R(dir_idx(1))) - vel_R(dir_idx(1)))))
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, eqn_idx%mom%end) = flux_rs${XYZ}$_vf(j, k, l, &
                                                               & eqn_idx%mom%beg + 1)
                                    end if
                                #:endif
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                else
                    ! 5-equation model (model_eqns=2): mixture total energy, volume fraction advection Private list split across
                    ! _hllc_p1/p2/p3 for Fypp line-length limits
                    #:set _hllc_p1 = '[Re_max, i, j, k, l, q, T_L, T_R, vel_L_rms, vel_R_rms, pres_L, pres_R, rho_L, gamma_L, pi_inf_L, qv_L, rho_R, gamma_R, pi_inf_R, qv_R, alpha_L_sum, alpha_R_sum, E_L, E_R, MW_L, MW_R, R_gas_L, R_gas_R, Cp_L, Cp_R, Cv_L, Cv_R, Cp_avg, Cv_avg, T_avg, eps, c_sum_Yi_Phi, Gamm_L, Gamm_R, Y_L, Y_R, H_L, H_R, qv_avg, rho_avg, gamma_avg, H_avg, c_L, c_R, c_avg, s_P, s_M, xi_P, xi_M, xi_L, xi_R, Ms_L, Ms_R, pres_SL, pres_SR, vel_L, vel_R, Re_L, Re_R, alpha_L, alpha_R, alpha_rho_L, alpha_rho_R, s_L, s_R, s_S, vel_avg_rms, pcorr, zcoef, ptilde_L, ptilde_R, vel_L_tmp, vel_R_tmp, Ys_L, Ys_R, Xs_L, Xs_R, Gamma_iL, Gamma_iR, Cp_iL, Cp_iR, tau_e_L, tau_e_R, xi_field_L, xi_field_R, Yi_avg, Phi_avg, h_iL, h_iR, h_avg_2, G_L, G_R, flux_ene_e, flux_tau_L, flux_tau_R,'
                    #:set _hllc_p2 = 'U_L, U_R, U_star_L, U_star_R, F_L, F_R, F_star_L, F_star_R, F_HLLC, u_n_HLLC, u_t_HLLC, u_t2_HLLC, pres_tot_L, pres_tot_R, u_n_L, u_n_R, u_t_L, u_t_R, u_t2_L, u_t2_R, tau_nn_L, tau_nn_R, tau_nt_L, tau_nt_R, tau_tt_L, tau_tt_R, tau_nt2_L, tau_nt2_R, tau_t2t2_L, tau_t2t2_R, tau_t1t2_L, tau_t1t2_R, tau_qq_L, tau_qq_R, p_face, tau_qq_face, A_L, A_R, denom_A, denom_L, denom_R, fac_L, fac_R, rho_L_star, rho_R_star, u_t_star, tau_nt_star, u_t2_star, tau_nt2_star, pres_tot_star, E_L_star, E_R_star, S_Mid, mom_n, mom_t1, mom_t2,'
                    #:set _hllc_p3 = 'F_HLL, u_n_HLL_trace, u_t_HLL_trace, u_t2_HLL_trace, p_face_HLL, tau_qq_face_HLL, tau_nn_HLL, phi, Sigma_L, Sigma_R, dSigma, Sigma_ref, a_L_ref, a_R_ref, a_ref, du_t, dtau_nt, du_t2, dtau_nt2, sensor_ptot, sensor_vt, sensor_tnt, sensor_combined, idx_phys]'
                    $:GPU_PARALLEL_LOOP(collapse=3, private=_hllc_p1 + _hllc_p2 + _hllc_p3, copyin='[is1, is2, is3]')
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end
                                vel_L_rms = 0._wp; vel_R_rms = 0._wp
                                rho_L = 0._wp; rho_R = 0._wp
                                gamma_L = 0._wp; gamma_R = 0._wp
                                pi_inf_L = 0._wp; pi_inf_R = 0._wp
                                qv_L = 0._wp; qv_R = 0._wp
                                alpha_L_sum = 0._wp; alpha_R_sum = 0._wp

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    alpha_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%E + i)
                                    alpha_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%E + i)
                                end do

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%cont%end + i)
                                    vel_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%cont%end + i)
                                    vel_L_rms = vel_L_rms + vel_L(i)**2._wp
                                    vel_R_rms = vel_R_rms + vel_R(i)**2._wp
                                end do

                                pres_L = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%E)
                                pres_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%E)

                                if (hypoelasticity) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, eqn_idx%stress%end - eqn_idx%stress%beg + 1
                                        tau_e_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%stress%beg - 1 + i)
                                        tau_e_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%stress%beg - 1 + i)
                                    end do

                                    ! Map physical-basis arrays to directional aliases via stress_perm/dir_idx
                                    u_n_L = vel_L(dir_idx(1)); u_n_R = vel_R(dir_idx(1))
                                    tau_nn_L = tau_e_L(stress_perm(1)); tau_nn_R = tau_e_R(stress_perm(1))
                                    if (n > 0) then
                                        u_t_L = vel_L(dir_idx(2)); u_t_R = vel_R(dir_idx(2))
                                        tau_nt_L = tau_e_L(stress_perm(2)); tau_nt_R = tau_e_R(stress_perm(2))
                                        tau_tt_L = tau_e_L(stress_perm(3)); tau_tt_R = tau_e_R(stress_perm(3))
                                    end if
                                    if (p > 0) then
                                        u_t2_L = vel_L(dir_idx(3)); u_t2_R = vel_R(dir_idx(3))
                                        tau_nt2_L = tau_e_L(stress_perm(4)); tau_nt2_R = tau_e_R(stress_perm(4))
                                        tau_t1t2_L = tau_e_L(stress_perm(5)); tau_t1t2_R = tau_e_R(stress_perm(5))
                                        tau_t2t2_L = tau_e_L(stress_perm(6)); tau_t2t2_R = tau_e_R(stress_perm(6))
                                    end if
                                    pres_tot_L = pres_L - tau_nn_L
                                    pres_tot_R = pres_R - tau_nn_R
                                    if (cyl_coord) then
                                        tau_qq_L = tau_e_L(eqn_idx%stress%end - eqn_idx%stress%beg + 1)
                                        tau_qq_R = tau_e_R(eqn_idx%stress%end - eqn_idx%stress%beg + 1)
                                    else
                                        tau_qq_L = 0._wp
                                        tau_qq_R = 0._wp
                                    end if
                                end if

                                ! Change this by splitting it into the cases present in the bubbles_euler
                                if (mpp_lim) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        qL_prim_rs${XYZ}$_vf(j, k, l, i) = max(0._wp, qL_prim_rs${XYZ}$_vf(j, k, l, i))
                                        qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%E + i) = min(max(0._wp, qL_prim_rs${XYZ}$_vf(j, k, &
                                                             & l, eqn_idx%E + i)), 1._wp)
                                        qR_prim_rs${XYZ}$_vf(j + 1, k, l, i) = max(0._wp, qR_prim_rs${XYZ}$_vf(j + 1, k, l, i))
                                        qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%E + i) = min(max(0._wp, &
                                                             & qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%E + i)), 1._wp)
                                        alpha_L_sum = alpha_L_sum + qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%E + i)
                                        alpha_R_sum = alpha_R_sum + qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%E + i)
                                    end do

                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%E + i) = qL_prim_rs${XYZ}$_vf(j, k, l, &
                                                             & eqn_idx%E + i)/max(alpha_L_sum, sgm_eps)
                                        qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%E + i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, &
                                                             & eqn_idx%E + i)/max(alpha_R_sum, sgm_eps)
                                    end do
                                end if

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    alpha_rho_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                    alpha_rho_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)

                                    rho_L = rho_L + alpha_rho_L(i)
                                    gamma_L = gamma_L + qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%E + i)*gammas(i)
                                    pi_inf_L = pi_inf_L + qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%E + i)*pi_infs(i)
                                    qv_L = qv_L + alpha_rho_L(i)*qvs(i)

                                    rho_R = rho_R + alpha_rho_R(i)
                                    gamma_R = gamma_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%E + i)*gammas(i)
                                    pi_inf_R = pi_inf_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%E + i)*pi_infs(i)
                                    qv_R = qv_R + alpha_rho_R(i)*qvs(i)
                                end do

                                Re_max = 0
                                if (Re_size(1) > 0) Re_max = 1
                                if (Re_size(2) > 0) Re_max = 2

                                if (viscous) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, Re_max
                                        Re_L(i) = 0._wp
                                        Re_R(i) = 0._wp

                                        $:GPU_LOOP(parallelism='[seq]')
                                        do q = 1, Re_size(i)
                                            Re_L(i) = alpha_L(Re_idx(i, q))/Res_gs(i, q) + Re_L(i)
                                            Re_R(i) = alpha_R(Re_idx(i, q))/Res_gs(i, q) + Re_R(i)
                                        end do

                                        Re_L(i) = 1._wp/max(Re_L(i), sgm_eps)
                                        Re_R(i) = 1._wp/max(Re_R(i), sgm_eps)
                                    end do
                                end if

                                if (chemistry) then
                                    c_sum_Yi_Phi = 0.0_wp
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = eqn_idx%species%beg, eqn_idx%species%end
                                        Ys_L(i - eqn_idx%species%beg + 1) = qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                        Ys_R(i - eqn_idx%species%beg + 1) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)
                                    end do

                                    call get_mixture_molecular_weight(Ys_L, MW_L)
                                    call get_mixture_molecular_weight(Ys_R, MW_R)

                                    #:if USING_AMD
                                        Xs_L(:) = Ys_L(:)*MW_L/molecular_weights_nonparameter(:)
                                        Xs_R(:) = Ys_R(:)*MW_R/molecular_weights_nonparameter(:)
                                    #:else
                                        Xs_L(:) = Ys_L(:)*MW_L/molecular_weights(:)
                                        Xs_R(:) = Ys_R(:)*MW_R/molecular_weights(:)
                                    #:endif

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
                                        tau_e_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%stress%beg - 1 + i)
                                        tau_e_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%stress%beg - 1 + i)
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
                                        ! Elastic energy (guard skips when G near zero)
                                        if (.not. hypo_energy_guard .or. ((G_L > verysmall) .and. (G_R > verysmall))) then
                                            E_L = E_L + (tau_e_L(i)*tau_e_L(i))/max(4._wp*G_L, verysmall)
                                            E_R = E_R + (tau_e_R(i)*tau_e_R(i))/max(4._wp*G_R, verysmall)
                                            ! Additional terms in 2D and 3D
                                            if ((i == 2) .or. (i == 4) .or. (i == 5)) then
                                                E_L = E_L + (tau_e_L(i)*tau_e_L(i))/max(4._wp*G_L, verysmall)
                                                E_R = E_R + (tau_e_R(i)*tau_e_R(i))/max(4._wp*G_R, verysmall)
                                            end if
                                        end if
                                    end do
                                end if

                                ! Hyperelastic stress contribution: strain energy added to total energy
                                if (hyperelasticity) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        xi_field_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%xi%beg - 1 + i)
                                        xi_field_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%xi%beg - 1 + i)
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
                                        E_L = E_L + G_L*qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%xi%end + 1)
                                        E_R = E_R + G_R*qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%xi%end + 1)
                                    end if
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, b_size - 1
                                        tau_e_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%stress%beg - 1 + i)
                                        tau_e_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%stress%beg - 1 + i)
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
                                        Re_avg_rs${XYZ}$_vf(j, k, l, i) = 2._wp/(1._wp/Re_L(i) + 1._wp/Re_R(i))
                                    end do
                                end if

                                ! Low Mach correction
                                if (low_Mach == 2) then
                                    @:compute_low_Mach_correction()
                                end if

                                if (wave_speeds == 1) then
                                    if (elasticity) then
                                        ! Elastic wave speed, Rodriguez et al. JCP (2019)
                                        s_L = min(vel_L(dir_idx(1)) - sqrt(max(verysmall, &
                                                  & c_L*c_L + (((4._wp*G_L)/3._wp) + tau_e_L(dir_idx_tau(1)))/rho_L)), &
                                                  & vel_R(dir_idx(1)) - sqrt(max(verysmall, &
                                                  & c_R*c_R + (((4._wp*G_R)/3._wp) + tau_e_R(dir_idx_tau(1)))/rho_R)))
                                        s_R = max(vel_R(dir_idx(1)) + sqrt(max(verysmall, &
                                                  & c_R*c_R + (((4._wp*G_R)/3._wp) + tau_e_R(dir_idx_tau(1)))/rho_R)), &
                                                  & vel_L(dir_idx(1)) + sqrt(max(verysmall, &
                                                  & c_L*c_L + (((4._wp*G_L)/3._wp) + tau_e_L(dir_idx_tau(1)))/rho_L)))
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
                                else if (wave_speeds == 2) then
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
                                xi_L = (s_L - vel_L(dir_idx(1)))/(s_L - s_S)
                                xi_R = (s_R - vel_R(dir_idx(1)))/(s_R - s_S)

                                ! goes with numerical velocity in x/y/z directions xi_P/M = 0.5 +/m sgn(0.5,s_star)
                                xi_M = (5.e-1_wp + sign(5.e-1_wp, s_S))
                                xi_P = (5.e-1_wp - sign(5.e-1_wp, s_S))

                                ! Low Mach correction
                                if (low_Mach == 1) then
                                    @:compute_low_Mach_correction()
                                else
                                    pcorr = 0._wp
                                end if

                                if (hypoelasticity) then
                                    if (n == 0) then
                                        u_t_L = 0._wp; u_t_R = 0._wp
                                        tau_nt_L = 0._wp; tau_nt_R = 0._wp
                                    end if
                                    if (p == 0) then
                                        u_t2_L = 0._wp; u_t2_R = 0._wp
                                        tau_nt2_L = 0._wp; tau_nt2_R = 0._wp
                                    end if
                                    A_L = rho_L*(s_L - vel_L(dir_idx(1)))
                                    A_R = rho_R*(s_R - vel_R(dir_idx(1)))
                                    denom_A = A_R - A_L
                                    u_t_star = (A_R*u_t_R - A_L*u_t_L + (tau_nt_R - tau_nt_L))/(denom_A + sgm_eps)
                                    tau_nt_star = (A_R*tau_nt_R - A_L*tau_nt_L)/(denom_A + sgm_eps)
                                    u_t2_star = (A_R*u_t2_R - A_L*u_t2_L + (tau_nt2_R - tau_nt2_L))/(denom_A + sgm_eps)
                                    tau_nt2_star = (A_R*tau_nt2_R - A_L*tau_nt2_L)/(denom_A + sgm_eps)
                                    pres_tot_star = pres_tot_L + A_L*(s_S - vel_L(dir_idx(1)))
                                end if

                                ! COMPUTING THE HLLC FLUXES MASS FLUX.
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, eqn_idx%cont%end
                                    flux_rs${XYZ}$_vf(j, k, l, i) = xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, &
                                                      & i)*(vel_L(dir_idx(1)) + s_M*(xi_L - 1._wp)) + xi_P*qR_prim_rs${XYZ}$_vf(j &
                                                      & + 1, k, l, i)*(vel_R(dir_idx(1)) + s_P*(xi_R - 1._wp))
                                end do

                                if (hypoelasticity) then
                                    flux_rs${XYZ}$_vf(j, k, l, &
                                                      & eqn_idx%cont%end + dir_idx(1)) = xi_M*(rho_L*(vel_L(dir_idx(1)) &
                                                      & *vel_L(dir_idx(1)) + s_M*(xi_L*s_S - vel_L(dir_idx(1)))) + pres_tot_L) &
                                                      & + xi_P*(rho_R*(vel_R(dir_idx(1))*vel_R(dir_idx(1)) + s_P*(xi_R*s_S &
                                                      & - vel_R(dir_idx(1)))) + pres_tot_R) + (s_M/s_L)*(s_P/s_R)*pcorr
                                    if (n > 0) then
                                        flux_rs${XYZ}$_vf(j, k, l, &
                                                          & eqn_idx%cont%end + dir_idx(2)) = xi_M*(rho_L*(vel_L(dir_idx(1))*u_t_L &
                                                          & + s_M*(xi_L*u_t_star - u_t_L)) - tau_nt_L) &
                                                          & + xi_P*(rho_R*(vel_R(dir_idx(1))*u_t_R + s_P*(xi_R*u_t_star - u_t_R)) &
                                                          & - tau_nt_R)
                                    end if
                                    if (p > 0) then
                                        flux_rs${XYZ}$_vf(j, k, l, &
                                                          & eqn_idx%cont%end + dir_idx(3)) = xi_M*(rho_L*(vel_L(dir_idx(1)) &
                                                          & *u_t2_L + s_M*(xi_L*u_t2_star - u_t2_L)) - tau_nt2_L) &
                                                          & + xi_P*(rho_R*(vel_R(dir_idx(1))*u_t2_R + s_P*(xi_R*u_t2_star &
                                                          & - u_t2_R)) - tau_nt2_R)
                                    end if

                                    flux_rs${XYZ}$_vf(j, k, l, &
                                                      & eqn_idx%E) = xi_M*((E_L + pres_tot_L)*vel_L(dir_idx(1)) - u_t_L*tau_nt_L &
                                                      & - u_t2_L*tau_nt2_L + s_M*(xi_L*(E_L + (s_S - vel_L(dir_idx(1))) &
                                                      & *(rho_L*s_S + pres_tot_L/(s_L - vel_L(dir_idx(1)))) + (u_t_L*tau_nt_L &
                                                      & - u_t_star*tau_nt_star)/(s_L - vel_L(dir_idx(1))) + (u_t2_L*tau_nt2_L &
                                                      & - u_t2_star*tau_nt2_star)/(s_L - vel_L(dir_idx(1)))) - E_L)) + xi_P*((E_R &
                                                      & + pres_tot_R)*vel_R(dir_idx(1)) - u_t_R*tau_nt_R - u_t2_R*tau_nt2_R &
                                                      & + s_P*(xi_R*(E_R + (s_S - vel_R(dir_idx(1)))*(rho_R*s_S + pres_tot_R/(s_R &
                                                      & - vel_R(dir_idx(1)))) + (u_t_R*tau_nt_R - u_t_star*tau_nt_star)/(s_R &
                                                      & - vel_R(dir_idx(1))) + (u_t2_R*tau_nt2_R - u_t2_star*tau_nt2_star)/(s_R &
                                                      & - vel_R(dir_idx(1)))) - E_R)) + (s_M/s_L)*(s_P/s_R)*pcorr*s_S

                                    if (n == 0) then
                                        flux_rs${XYZ}$_vf(j, k, l, &
                                                          & eqn_idx%stress%beg) = xi_M*rho_L*tau_nn_L*(vel_L(dir_idx(1)) &
                                                          & + s_M*(xi_L - 1._wp)) + xi_P*rho_R*tau_nn_R*(vel_R(dir_idx(1)) &
                                                          & + s_P*(xi_R - 1._wp))
                                    else if (p == 0) then
                                        if (dir_idx(1) == 1) then
                                            flux_rs${XYZ}$_vf(j, k, l, &
                                                              & eqn_idx%stress%beg) = xi_M*rho_L*tau_nn_L*(vel_L(dir_idx(1)) &
                                                              & + s_M*(xi_L - 1._wp)) + xi_P*rho_R*tau_nn_R*(vel_R(dir_idx(1)) &
                                                              & + s_P*(xi_R - 1._wp))
                                            flux_rs${XYZ}$_vf(j, k, l, &
                                                              & eqn_idx%stress%beg + 1) = xi_M*(rho_L*vel_L(dir_idx(1))*tau_nt_L &
                                                              & + s_M*(rho_L*xi_L*tau_nt_star - rho_L*tau_nt_L)) &
                                                              & + xi_P*(rho_R*vel_R(dir_idx(1))*tau_nt_R &
                                                              & + s_P*(rho_R*xi_R*tau_nt_star - rho_R*tau_nt_R))
                                            flux_rs${XYZ}$_vf(j, k, l, &
                                                              & eqn_idx%stress%beg + 2) = xi_M*rho_L*tau_tt_L*(vel_L(dir_idx(1)) &
                                                              & + s_M*(xi_L - 1._wp)) + xi_P*rho_R*tau_tt_R*(vel_R(dir_idx(1)) &
                                                              & + s_P*(xi_R - 1._wp))
                                        else
                                            flux_rs${XYZ}$_vf(j, k, l, &
                                                              & eqn_idx%stress%beg + 2) = xi_M*rho_L*tau_nn_L*(vel_L(dir_idx(1)) &
                                                              & + s_M*(xi_L - 1._wp)) + xi_P*rho_R*tau_nn_R*(vel_R(dir_idx(1)) &
                                                              & + s_P*(xi_R - 1._wp))
                                            flux_rs${XYZ}$_vf(j, k, l, &
                                                              & eqn_idx%stress%beg + 1) = xi_M*(rho_L*vel_L(dir_idx(1))*tau_nt_L &
                                                              & + s_M*(rho_L*xi_L*tau_nt_star - rho_L*tau_nt_L)) &
                                                              & + xi_P*(rho_R*vel_R(dir_idx(1))*tau_nt_R &
                                                              & + s_P*(rho_R*xi_R*tau_nt_star - rho_R*tau_nt_R))
                                            flux_rs${XYZ}$_vf(j, k, l, &
                                                              & eqn_idx%stress%beg) = xi_M*rho_L*tau_tt_L*(vel_L(dir_idx(1)) &
                                                              & + s_M*(xi_L - 1._wp)) + xi_P*rho_R*tau_tt_R*(vel_R(dir_idx(1)) &
                                                              & + s_P*(xi_R - 1._wp))
                                        end if
                                    else
                                        flux_rs${XYZ}$_vf(j, k, l, &
                                                          & eqn_idx%stress%beg - 1 + stress_perm(1)) &
                                                          & = xi_M*rho_L*tau_nn_L*(vel_L(dir_idx(1)) + s_M*(xi_L - 1._wp)) &
                                                          & + xi_P*rho_R*tau_nn_R*(vel_R(dir_idx(1)) + s_P*(xi_R - 1._wp))
                                        flux_rs${XYZ}$_vf(j, k, l, &
                                                          & eqn_idx%stress%beg - 1 + stress_perm(2)) &
                                                          & = xi_M*(rho_L*vel_L(dir_idx(1))*tau_nt_L &
                                                          & + s_M*(rho_L*xi_L*tau_nt_star - rho_L*tau_nt_L)) &
                                                          & + xi_P*(rho_R*vel_R(dir_idx(1))*tau_nt_R &
                                                          & + s_P*(rho_R*xi_R*tau_nt_star - rho_R*tau_nt_R))
                                        flux_rs${XYZ}$_vf(j, k, l, &
                                                          & eqn_idx%stress%beg - 1 + stress_perm(4)) &
                                                          & = xi_M*(rho_L*vel_L(dir_idx(1))*tau_nt2_L &
                                                          & + s_M*(rho_L*xi_L*tau_nt2_star - rho_L*tau_nt2_L)) &
                                                          & + xi_P*(rho_R*vel_R(dir_idx(1))*tau_nt2_R &
                                                          & + s_P*(rho_R*xi_R*tau_nt2_star - rho_R*tau_nt2_R))
                                        flux_rs${XYZ}$_vf(j, k, l, &
                                                          & eqn_idx%stress%beg - 1 + stress_perm(3)) &
                                                          & = xi_M*rho_L*tau_tt_L*(vel_L(dir_idx(1)) + s_M*(xi_L - 1._wp)) &
                                                          & + xi_P*rho_R*tau_tt_R*(vel_R(dir_idx(1)) + s_P*(xi_R - 1._wp))
                                        flux_rs${XYZ}$_vf(j, k, l, &
                                                          & eqn_idx%stress%beg - 1 + stress_perm(6)) &
                                                          & = xi_M*rho_L*tau_t2t2_L*(vel_L(dir_idx(1)) + s_M*(xi_L - 1._wp)) &
                                                          & + xi_P*rho_R*tau_t2t2_R*(vel_R(dir_idx(1)) + s_P*(xi_R - 1._wp))
                                        flux_rs${XYZ}$_vf(j, k, l, &
                                                          & eqn_idx%stress%beg - 1 + stress_perm(5)) &
                                                          & = xi_M*rho_L*tau_t1t2_L*(vel_L(dir_idx(1)) + s_M*(xi_L - 1._wp)) &
                                                          & + xi_P*rho_R*tau_t1t2_R*(vel_R(dir_idx(1)) + s_P*(xi_R - 1._wp))
                                    end if
                                    if (cyl_coord) then
                                        flux_rs${XYZ}$_vf(j, k, l, &
                                                          & eqn_idx%stress%end) = xi_M*rho_L*tau_qq_L*(vel_L(dir_idx(1)) &
                                                          & + s_M*(xi_L - 1._wp)) + xi_P*rho_R*tau_qq_R*(vel_R(dir_idx(1)) &
                                                          & + s_P*(xi_R - 1._wp))
                                    end if

                                    if (s_L >= 0._wp) then
                                        u_n_HLLC = vel_L(dir_idx(1)); u_t_HLLC = u_t_L; u_t2_HLLC = u_t2_L
                                    else if (s_R <= 0._wp) then
                                        u_n_HLLC = vel_R(dir_idx(1)); u_t_HLLC = u_t_R; u_t2_HLLC = u_t2_R
                                    else
                                        u_n_HLLC = s_S*(xi_M*xi_L + xi_P*xi_R); u_t_HLLC = u_t_star; u_t2_HLLC = u_t2_star
                                    end if
                                    nc_iface_vel_rs${XYZ}$_vf(j, k, l, dir_idx(1)) = u_n_HLLC
                                    if (n > 0) nc_iface_vel_rs${XYZ}$_vf(j, k, l, dir_idx(2)) = u_t_HLLC
                                    if (p > 0) nc_iface_vel_rs${XYZ}$_vf(j, k, l, dir_idx(3)) = u_t2_HLLC
                                else  ! not hypoelasticity
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        ! MOMENTUM FLUX.
                                        flux_rs${XYZ}$_vf(j, k, l, &
                                                          & eqn_idx%cont%end + dir_idx(i)) = xi_M*(rho_L*(vel_L(dir_idx(1)) &
                                                          & *vel_L(dir_idx(i)) + s_M*(xi_L*(dir_flg(dir_idx(i))*s_S + (1._wp &
                                                          & - dir_flg(dir_idx(i)))*vel_L(dir_idx(i))) - vel_L(dir_idx(i)))) &
                                                          & + dir_flg(dir_idx(i))*(pres_L)) + xi_P*(rho_R*(vel_R(dir_idx(1)) &
                                                          & *vel_R(dir_idx(i)) + s_P*(xi_R*(dir_flg(dir_idx(i))*s_S + (1._wp &
                                                          & - dir_flg(dir_idx(i)))*vel_R(dir_idx(i))) - vel_R(dir_idx(i)))) &
                                                          & + dir_flg(dir_idx(i))*(pres_R)) + (s_M/s_L)*(s_P/s_R) &
                                                          & *dir_flg(dir_idx(i))*pcorr
                                    end do

                                    ! ENERGY FLUX.
                                    flux_rs${XYZ}$_vf(j, k, l, &
                                                      & eqn_idx%E) = xi_M*(vel_L(dir_idx(1))*(E_L + pres_L) + s_M*(xi_L*(E_L &
                                                      & + (s_S - vel_L(dir_idx(1)))*(rho_L*s_S + pres_L/(s_L - vel_L(dir_idx(1)))) &
                                                      & ) - E_L)) + xi_P*(vel_R(dir_idx(1))*(E_R + pres_R) + s_P*(xi_R*(E_R &
                                                      & + (s_S - vel_R(dir_idx(1)))*(rho_R*s_S + pres_R/(s_R - vel_R(dir_idx(1)))) &
                                                      & ) - E_R)) + (s_M/s_L)*(s_P/s_R)*pcorr*s_S

                                    if (elasticity) then
                                        flux_ene_e = 0._wp
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, num_dims
                                            ! MOMENTUM ELASTIC FLUX.
                                            flux_rs${XYZ}$_vf(j, k, l, eqn_idx%cont%end + dir_idx(i)) = flux_rs${XYZ}$_vf(j, k, &
                                                              & l, &
                                                              & eqn_idx%cont%end + dir_idx(i)) - xi_M*tau_e_L(dir_idx_tau(i)) &
                                                              & - xi_P*tau_e_R(dir_idx_tau(i))
                                            ! ENERGY ELASTIC FLUX.
                                            flux_ene_e = flux_ene_e - xi_M*(vel_L(dir_idx(i))*tau_e_L(dir_idx_tau(i)) &
                                                                            & + s_M*(xi_L*((s_S - vel_L(i)) &
                                                                            & *(tau_e_L(dir_idx_tau(i))/(s_L - vel_L(i)))))) &
                                                                            & - xi_P*(vel_R(dir_idx(i))*tau_e_R(dir_idx_tau(i)) &
                                                                            & + s_P*(xi_R*((s_S - vel_R(i)) &
                                                                            & *(tau_e_R(dir_idx_tau(i))/(s_R - vel_R(i))))))
                                        end do
                                        flux_rs${XYZ}$_vf(j, k, l, eqn_idx%E) = flux_rs${XYZ}$_vf(j, k, l, eqn_idx%E) + flux_ene_e
                                    end if
                                end if

                                ! VOLUME FRACTION FLUX.
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = eqn_idx%adv%beg, eqn_idx%adv%end
                                    flux_rs${XYZ}$_vf(j, k, l, i) = xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, &
                                                      & i)*(vel_L(dir_idx(1)) + s_M*(xi_L - 1._wp)) + xi_P*qR_prim_rs${XYZ}$_vf(j &
                                                      & + 1, k, l, i)*(vel_R(dir_idx(1)) + s_P*(xi_R - 1._wp))
                                end do

                                ! VOLUME FRACTION SOURCE FLUX.
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel_src_rs${XYZ}$_vf(j, k, l, &
                                                         & dir_idx(i)) = xi_M*(vel_L(dir_idx(i)) + dir_flg(dir_idx(i))*s_M*(xi_L &
                                                         & - 1._wp)) + xi_P*(vel_R(dir_idx(i)) + dir_flg(dir_idx(i))*s_P*(xi_R &
                                                         & - 1._wp))
                                end do

                                ! COLOR FUNCTION FLUX
                                if (surface_tension) then
                                    flux_rs${XYZ}$_vf(j, k, l, eqn_idx%c) = xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, &
                                                      & eqn_idx%c)*(vel_L(dir_idx(1)) + s_M*(xi_L - 1._wp)) &
                                                      & + xi_P*qR_prim_rs${XYZ}$_vf(j + 1, k, l, &
                                                      & eqn_idx%c)*(vel_R(dir_idx(1)) + s_P*(xi_R - 1._wp))
                                end if

                                ! Hyperelastic reference map flux for material deformation tracking
                                if (hyperelasticity) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        flux_rs${XYZ}$_vf(j, k, l, &
                                                          & eqn_idx%xi%beg - 1 + i) = xi_M*(s_S/(s_L - s_S)) &
                                                          & *(s_L*rho_L*xi_field_L(i) - rho_L*vel_L(dir_idx(1))*xi_field_L(i)) &
                                                          & + xi_P*(s_S/(s_R - s_S))*(s_R*rho_R*xi_field_R(i) &
                                                          & - rho_R*vel_R(dir_idx(1))*xi_field_R(i))
                                    end do
                                end if

                                flux_src_rs${XYZ}$_vf(j, k, l, eqn_idx%adv%beg) = vel_src_rs${XYZ}$_vf(j, k, l, dir_idx(1))

                                if (chemistry) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = eqn_idx%species%beg, eqn_idx%species%end
                                        Y_L = qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                        Y_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)

                                        flux_rs${XYZ}$_vf(j, k, l, &
                                                          & i) = xi_M*rho_L*Y_L*(vel_L(dir_idx(1)) + s_M*(xi_L - 1._wp)) &
                                                          & + xi_P*rho_R*Y_R*(vel_R(dir_idx(1)) + s_P*(xi_R - 1._wp))
                                        flux_src_rs${XYZ}$_vf(j, k, l, i) = 0.0_wp
                                    end do
                                end if

                                ! HLLC-ADC blending for hypoelasticity
                                if (riemann_hypo_ADC .and. hypoelasticity) then
                                    ! Build U_L, U_R and F_L, F_R in local-basis layout
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        U_L(i) = alpha_rho_L(i)
                                        U_R(i) = alpha_rho_R(i)
                                        U_L(eqn_idx%adv%beg - 1 + i) = alpha_L(i)
                                        U_R(eqn_idx%adv%beg - 1 + i) = alpha_R(i)
                                        F_L(i) = alpha_rho_L(i)*u_n_L
                                        F_R(i) = alpha_rho_R(i)*u_n_R
                                        F_L(eqn_idx%adv%beg - 1 + i) = alpha_L(i)*u_n_L
                                        F_R(eqn_idx%adv%beg - 1 + i) = alpha_R(i)*u_n_R
                                    end do

                                    ! Momentum U/F in physical order via dir_idx
                                    U_L(eqn_idx%cont%end + dir_idx(1)) = rho_L*u_n_L
                                    U_R(eqn_idx%cont%end + dir_idx(1)) = rho_R*u_n_R
                                    F_L(eqn_idx%cont%end + dir_idx(1)) = rho_L*u_n_L*u_n_L + pres_tot_L
                                    F_R(eqn_idx%cont%end + dir_idx(1)) = rho_R*u_n_R*u_n_R + pres_tot_R
                                    if (n > 0) then
                                        U_L(eqn_idx%cont%end + dir_idx(2)) = rho_L*u_t_L
                                        U_R(eqn_idx%cont%end + dir_idx(2)) = rho_R*u_t_R
                                        F_L(eqn_idx%cont%end + dir_idx(2)) = rho_L*u_n_L*u_t_L - tau_nt_L
                                        F_R(eqn_idx%cont%end + dir_idx(2)) = rho_R*u_n_R*u_t_R - tau_nt_R
                                    end if
                                    if (p > 0) then
                                        U_L(eqn_idx%cont%end + dir_idx(3)) = rho_L*u_t2_L
                                        U_R(eqn_idx%cont%end + dir_idx(3)) = rho_R*u_t2_R
                                        F_L(eqn_idx%cont%end + dir_idx(3)) = rho_L*u_n_L*u_t2_L - tau_nt2_L
                                        F_R(eqn_idx%cont%end + dir_idx(3)) = rho_R*u_n_R*u_t2_R - tau_nt2_R
                                    end if

                                    U_L(eqn_idx%E) = E_L
                                    U_R(eqn_idx%E) = E_R
                                    F_L(eqn_idx%E) = (E_L + pres_tot_L)*u_n_L - u_t_L*tau_nt_L - u_t2_L*tau_nt2_L
                                    F_R(eqn_idx%E) = (E_R + pres_tot_R)*u_n_R - u_t_R*tau_nt_R - u_t2_R*tau_nt2_R

                                    ! Stress U/F in physical order via stress_perm: U = rho*tau, F = rho*u_n*tau
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, eqn_idx%stress%end - eqn_idx%stress%beg + 1 - merge(1, 0, cyl_coord)
                                        idx_phys = eqn_idx%stress%beg - 1 + stress_perm(i)
                                        U_L(idx_phys) = rho_L*tau_e_L(stress_perm(i))
                                        U_R(idx_phys) = rho_R*tau_e_R(stress_perm(i))
                                        F_L(idx_phys) = rho_L*u_n_L*tau_e_L(stress_perm(i))
                                        F_R(idx_phys) = rho_R*u_n_R*tau_e_R(stress_perm(i))
                                    end do
                                    if (cyl_coord) then
                                        U_L(eqn_idx%stress%end) = rho_L*tau_qq_L
                                        U_R(eqn_idx%stress%end) = rho_R*tau_qq_R
                                        F_L(eqn_idx%stress%end) = rho_L*u_n_L*tau_qq_L
                                        F_R(eqn_idx%stress%end) = rho_R*u_n_R*tau_qq_R
                                    end if

                                    ! Compute F_HLL (physical order) and HLL trace velocities
                                    if (s_L >= 0._wp) then
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, sys_size
                                            F_HLL(i) = F_L(i)
                                        end do
                                        u_n_HLL_trace = u_n_L; u_t_HLL_trace = u_t_L; u_t2_HLL_trace = u_t2_L
                                    else if (s_R <= 0._wp) then
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, sys_size
                                            F_HLL(i) = F_R(i)
                                        end do
                                        u_n_HLL_trace = u_n_R; u_t_HLL_trace = u_t_R; u_t2_HLL_trace = u_t2_R
                                    else
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, sys_size
                                            F_HLL(i) = (s_R*F_L(i) - s_L*F_R(i) + s_L*s_R*(U_R(i) - U_L(i)))/(s_R - s_L + verysmall)
                                        end do
                                        u_n_HLL_trace = (s_R*u_n_L - s_L*u_n_R)/(s_R - s_L + verysmall)
                                        u_t_HLL_trace = 0._wp; u_t2_HLL_trace = 0._wp
                                        if (n > 0) u_t_HLL_trace = (s_R*u_t_L - s_L*u_t_R)/(s_R - s_L + verysmall)
                                        if (p > 0) u_t2_HLL_trace = (s_R*u_t2_L - s_L*u_t2_R)/(s_R - s_L + verysmall)
                                    end if

                                    ! ADC sensor
                                    Sigma_L = pres_tot_L
                                    Sigma_R = pres_tot_R
                                    dSigma = Sigma_R - Sigma_L
                                    Sigma_ref = max(max(abs(Sigma_L), abs(Sigma_R)), verysmall)

                                    a_L_ref = sqrt(max(verysmall, c_L*c_L + ((4._wp/3._wp)*G_L + tau_nn_L)/rho_L))
                                    a_R_ref = sqrt(max(verysmall, c_R*c_R + ((4._wp/3._wp)*G_R + tau_nn_R)/rho_R))
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

                                    ! Blend all flux components: F_HLL is in physical order
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, sys_size
                                        flux_rs${XYZ}$_vf(j, k, l, i) = F_HLL(i) + phi*(flux_rs${XYZ}$_vf(j, k, l, i) - F_HLL(i))
                                    end do

                                    ! Blend interface velocities (scalar HLL traces)
                                    u_n_HLLC = u_n_HLL_trace + phi*(u_n_HLLC - u_n_HLL_trace)
                                    u_t_HLLC = u_t_HLL_trace + phi*(u_t_HLLC - u_t_HLL_trace)
                                    u_t2_HLLC = u_t2_HLL_trace + phi*(u_t2_HLLC - u_t2_HLL_trace)

                                    ! Overwrite vel_src with blended velocities
                                    vel_src_rs${XYZ}$_vf(j, k, l, dir_idx(1)) = u_n_HLLC
                                    if (n > 0) vel_src_rs${XYZ}$_vf(j, k, l, dir_idx(2)) = u_t_HLLC
                                    if (p > 0) vel_src_rs${XYZ}$_vf(j, k, l, dir_idx(3)) = u_t2_HLLC

                                    ! Update advection source flux with ADC-blended face-normal velocity
                                    flux_src_rs${XYZ}$_vf(j, k, l, eqn_idx%adv%beg) = u_n_HLLC

                                    ! Overwrite nc_iface_vel with blended velocities
                                    nc_iface_vel_rs${XYZ}$_vf(j, k, l, dir_idx(1)) = u_n_HLLC
                                    if (n > 0) nc_iface_vel_rs${XYZ}$_vf(j, k, l, dir_idx(2)) = u_t_HLLC
                                    if (p > 0) nc_iface_vel_rs${XYZ}$_vf(j, k, l, dir_idx(3)) = u_t2_HLLC
                                end if
                                ! END HLLC-ADC

                                ! Geometrical source flux for cylindrical coordinates
                                #:if (NORM_DIR == 2)
                                    if (cyl_coord .and. hypoelasticity) then
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, sys_size
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = flux_rs${XYZ}$_vf(j, k, l, i)
                                        end do
                                        if (s_L >= 0._wp) then
                                            p_face = pres_L; tau_qq_face = tau_qq_L
                                        else if (s_R <= 0._wp) then
                                            p_face = pres_R; tau_qq_face = tau_qq_R
                                        else if (s_S >= 0._wp) then
                                            p_face = pres_tot_star + tau_nn_L; tau_qq_face = tau_qq_L
                                        else
                                            p_face = pres_tot_star + tau_nn_R; tau_qq_face = tau_qq_R
                                        end if
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, eqn_idx%cont%end + dir_idx(1)) = flux_rs${XYZ}$_vf(j, k, &
                                                               & l, eqn_idx%cont%end + dir_idx(1)) - p_face + tau_qq_face
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = eqn_idx%adv%beg, eqn_idx%adv%end
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = 0._wp
                                        end do
                                    else if (cyl_coord) then
                                        ! Substituting the advective flux into the inviscid geometrical source flux
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, eqn_idx%E
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = flux_rs${XYZ}$_vf(j, k, l, i)
                                        end do
                                        ! Recalculating the radial momentum geometric source flux
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, &
                                                               & eqn_idx%cont%end + dir_idx(1)) = xi_M*(rho_L*(vel_L(dir_idx(1)) &
                                                               & *vel_L(dir_idx(1)) + s_M*(xi_L*(dir_flg(dir_idx(1))*s_S + (1._wp &
                                                               & - dir_flg(dir_idx(1)))*vel_L(dir_idx(1))) - vel_L(dir_idx(1))))) &
                                                               & + xi_P*(rho_R*(vel_R(dir_idx(1))*vel_R(dir_idx(1)) &
                                                               & + s_P*(xi_R*(dir_flg(dir_idx(1))*s_S + (1._wp &
                                                               & - dir_flg(dir_idx(1)))*vel_R(dir_idx(1))) - vel_R(dir_idx(1)))))
                                        ! Geometrical source of the void fraction(s) is zero
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = eqn_idx%adv%beg, eqn_idx%adv%end
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = 0._wp
                                        end do
                                    end if
                                #:endif
                                #:if (NORM_DIR == 3)
                                    if (grid_geometry == 3) then
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, sys_size
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = 0._wp
                                        end do

                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, &
                                                               & eqn_idx%mom%beg + 1) = -xi_M*(rho_L*(vel_L(dir_idx(1)) &
                                                               & *vel_L(dir_idx(1)) + s_M*(xi_L*(dir_flg(dir_idx(1))*s_S + (1._wp &
                                                               & - dir_flg(dir_idx(1)))*vel_L(dir_idx(1))) - vel_L(dir_idx(1))))) &
                                                               & - xi_P*(rho_R*(vel_R(dir_idx(1))*vel_R(dir_idx(1)) &
                                                               & + s_P*(xi_R*(dir_flg(dir_idx(1))*s_S + (1._wp &
                                                               & - dir_flg(dir_idx(1)))*vel_R(dir_idx(1))) - vel_R(dir_idx(1)))))
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, eqn_idx%mom%end) = flux_rs${XYZ}$_vf(j, k, l, &
                                                               & eqn_idx%mom%beg + 1)
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

        if (viscous .or. dummy) then
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
            call s_compute_capillary_source_flux(vel_src_rsx_vf, vel_src_rsy_vf, vel_src_rsz_vf, flux_src_vf, norm_dir, isx, isy, &
                                                 & isz)
        end if

        call s_finalize_riemann_solver(flux_vf, flux_src_vf, flux_gsrc_vf, norm_dir)

    end subroutine s_hllc_riemann_solver

    !> HLLD Riemann solver for MHD, Miyoshi & Kusano JCP (2005)
    subroutine s_hlld_riemann_solver(qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, dqL_prim_dx_vf, dqL_prim_dy_vf, &

        & dqL_prim_dz_vf, qL_prim_vf, qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, dqR_prim_dx_vf, dqR_prim_dy_vf, &
            & dqR_prim_dz_vf, qR_prim_vf, q_prim_vf, flux_vf, flux_src_vf, flux_gsrc_vf, norm_dir, ix, iy, iz)

        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout) :: qL_prim_rsx_vf, qL_prim_rsy_vf, &
             & qL_prim_rsz_vf, qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf

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

        call s_populate_riemann_states_variables_buffers(qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, dqL_prim_dx_vf, &
            & dqL_prim_dy_vf, dqL_prim_dz_vf, qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, dqR_prim_dx_vf, dqR_prim_dy_vf, &
            & dqR_prim_dz_vf, norm_dir, ix, iy, iz)

        call s_initialize_riemann_solver(flux_src_vf, norm_dir)

        #:for NORM_DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
            if (norm_dir == ${NORM_DIR}$) then
                $:GPU_PARALLEL_LOOP(collapse=3, private='[alpha_rho_L, alpha_rho_R, vel, alpha_L, alpha_R, rho, pres, E, &
                                    & H_no_mag, gamma, pi_inf, qv, vel_rms, B, c, c_fast, pres_mag, U_L, U_R, U_starL, U_starR, &
                                    & U_doubleL, U_doubleR, F_L, F_R, F_starL, F_starR, F_hlld, s_L, s_R, s_M, s_starL, s_starR, &
                                    & pTot_L, pTot_R, p_star, rhoL_star, rhoR_star, E_starL, E_starR, sqrt_rhoL_star, &
                                    & sqrt_rhoR_star, denom_ds, sign_Bx, vL_star, vR_star, wL_star, wR_star, v_double, w_double, &
                                    & By_double, Bz_double, E_doubleL, E_doubleR, E_double]', copyin='[norm_dir]')
                do l = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = is1%beg, is1%end
                            ! (1) Extract the left/right primitive states
                            do i = 1, eqn_idx%cont%end
                                alpha_rho_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                alpha_rho_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)
                            end do

                            ! NOTE: unlike HLL & HLLC, vel_L here is permutated by dir_idx for simpler logic
                            do i = 1, num_vels
                                vel%L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%cont%end + dir_idx(i))
                                vel%R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%cont%end + dir_idx(i))
                            end do

                            vel_rms%L = sum(vel%L**2._wp)
                            vel_rms%R = sum(vel%R**2._wp)

                            do i = 1, num_fluids
                                alpha_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%E + i)
                                alpha_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%E + i)
                            end do

                            pres%L = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%E)
                            pres%R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%E)

                            ! NOTE: unlike HLL, Bx, By, Bz are permutated by dir_idx for simpler logic
                            if (mhd) then
                                if (n == 0) then  ! 1D: constant Bx; By, Bz as variables; only in x so not permutated
                                    B%L = [Bx0, qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%B%beg), qL_prim_rs${XYZ}$_vf(j, k, l, &
                                                                     & eqn_idx%B%beg + 1)]
                                    B%R = [Bx0, qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%B%beg), qR_prim_rs${XYZ}$_vf(j + 1, k, &
                                                                     & l, eqn_idx%B%beg + 1)]
                                else  ! 2D/3D: Bx, By, Bz as variables
                                    B%L = [qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%B%beg + dir_idx(1) - 1), qL_prim_rs${XYZ}$_vf(j, &
                                                                & k, l, eqn_idx%B%beg + dir_idx(2) - 1), qL_prim_rs${XYZ}$_vf(j, &
                                                                & k, l, eqn_idx%B%beg + dir_idx(3) - 1)]
                                    B%R = [qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%B%beg + dir_idx(1) - 1), &
                                                                & qR_prim_rs${XYZ}$_vf(j + 1, k, l, &
                                                                & eqn_idx%B%beg + dir_idx(2) - 1), qR_prim_rs${XYZ}$_vf(j + 1, k, &
                                                                & l, eqn_idx%B%beg + dir_idx(3) - 1)]
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

                            ! (12) Write HLLD flux to output arrays
                            flux_rs${XYZ}$_vf(j, k, l, 1) = F_hlld(1)  ! TODO multi-component
                            ! Momentum
                            flux_rs${XYZ}$_vf(j, k, l, eqn_idx%cont%end + dir_idx(1)) = F_hlld(2)
                            flux_rs${XYZ}$_vf(j, k, l, eqn_idx%cont%end + dir_idx(2)) = F_hlld(3)
                            flux_rs${XYZ}$_vf(j, k, l, eqn_idx%cont%end + dir_idx(3)) = F_hlld(4)
                            ! Magnetic field
                            if (n == 0) then
                                flux_rs${XYZ}$_vf(j, k, l, eqn_idx%B%beg) = F_hlld(5)
                                flux_rs${XYZ}$_vf(j, k, l, eqn_idx%B%beg + 1) = F_hlld(6)
                            else
                                flux_rs${XYZ}$_vf(j, k, l, eqn_idx%B%beg + dir_idx(2) - 1) = F_hlld(5)
                                flux_rs${XYZ}$_vf(j, k, l, eqn_idx%B%beg + dir_idx(3) - 1) = F_hlld(6)
                            end if
                            ! Energy
                            flux_rs${XYZ}$_vf(j, k, l, eqn_idx%E) = F_hlld(7)
                            ! Volume fractions
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = eqn_idx%adv%beg, eqn_idx%adv%end
                                flux_rs${XYZ}$_vf(j, k, l, i) = 0._wp  ! TODO multi-component (zero for now)
                            end do

                            flux_src_rs${XYZ}$_vf(j, k, l, eqn_idx%adv%beg) = 0._wp
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

        call s_finalize_riemann_solver(flux_vf, flux_src_vf, flux_gsrc_vf, norm_dir)

    end subroutine s_hlld_riemann_solver

    !> HLLD Riemann solver resolves all 5 waves for the hypoelastic equations: 1 entropy wave, 2 shear stress waves, 2 fast waves.
    subroutine s_hypo_hlld_riemann_solver(qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, dqL_prim_dx_vf, dqL_prim_dy_vf, &
                                          & dqL_prim_dz_vf, qL_prim_vf, qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, &
                                          & dqR_prim_dx_vf, dqR_prim_dy_vf, dqR_prim_dz_vf, qR_prim_vf, q_prim_vf, flux_vf, &
                                          & flux_src_vf, flux_gsrc_vf, norm_dir, ix, iy, iz)

        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout) :: qL_prim_rsx_vf, qL_prim_rsy_vf, &
             & qL_prim_rsz_vf, qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf

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
        type(riemann_states)      :: rho, pres, E, H
        type(riemann_states)      :: gamma, pi_inf, qv
        type(riemann_states)      :: vel_rms
        type(riemann_states)      :: c

        ! HLLD speeds and intermediate state variables:
        real(wp)                :: S_L, S_R, s_M, S_Lstar, S_Rstar
        real(wp)                :: pTot_L, pTot_R, rhoL_star, rhoR_star
        real(wp), dimension(14) :: U_L, U_R, U_starL, U_starR, U_starstarL, U_starstarR
        real(wp), dimension(14) :: F_L, F_R, F_starL, F_starR, F_hlld
        real(wp), dimension(14) :: F_HLL  ! for ADC blending
        real(wp), dimension(14) :: U_HLL  ! for ADC axisym geometric flux
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
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
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
        integer             :: i, j, k, l, ipass

        call s_populate_riemann_states_variables_buffers(qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, dqL_prim_dx_vf, &
            & dqL_prim_dy_vf, dqL_prim_dz_vf, qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, dqR_prim_dx_vf, dqR_prim_dy_vf, &
            & dqR_prim_dz_vf, norm_dir, ix, iy, iz)

        call s_initialize_riemann_solver(flux_src_vf, norm_dir)

        #:for NORM_DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
            if (norm_dir == ${NORM_DIR}$) then
                ! Anchor-cell index pattern: the fused kernel reads both anchors (hat_L: cell j, hat_R: cell j+1) directly from
                ! q_prim_vf with the sweep-permuted indices the deleted q_hat fill kernel used.
                #:set HATIDX = {1: 'j + ipass - 1, k, l', 2: 'k, j + ipass - 1, l', 3: 'l, k, j + ipass - 1'}[NORM_DIR]
                #:set _hlld_p1 = '[i,j,k,l,ipass,degenerate,alpha_rho_L,alpha_rho_R,vel,alpha_L,alpha_R,rho,pres,E,H,gamma,pi_inf,qv,vel_rms,c,S_L,S_R,s_M,S_Lstar,S_Rstar,pTot_L,pTot_R,rhoL_star,rhoR_star,U_L,U_R,U_starL,U_starR,U_starstarL,U_starstarR,F_L,F_R,F_starL,F_starR,F_hlld,F_HLL,U_HLL,rho_HLL,u_n_HLL_cons,tau_nn_HLL,u_n_HLL_trace,u_t_HLL_trace,p_face_HLL,tau_qq_face_HLL,ncomp,C_NC,sqrtC_NC,A_L,A_R,denomA,fac_L,fac_R,'
                #:set _hlld_p2 = 'u_n_L,u_t_L,u_n_R,u_t_R,u_t2_L,u_t2_R,tau_nn_L,tau_nt_L,tau_tt_L,tau_nn_R,tau_nt_R,tau_tt_R,tau_nt2_L,tau_nt2_R,tau_t2t2_L,tau_t2t2_R,tau_t1t2_L,tau_t1t2_R,tau_qq_L,tau_qq_R,G_L,G_R,tau_e_L,tau_e_R,alpha1_L_star,alpha1_R_star,alpha2_L_star,alpha2_R_star,u_t_star,tau_nt_star,u_t2_star,tau_nt2_star,tau_nn_L_star,tau_nn_R_star,tau_tt_L_star,tau_tt_R_star,tau_tt_L_starstar,tau_tt_R_starstar,'
                #:set _hlld_p3 = 'tau_t2t2_L_star,tau_t2t2_R_star,tau_t2t2_L_starstar,tau_t2t2_R_starstar,tau_t1t2_L_star,tau_t1t2_R_star,tau_t1t2_L_starstar,tau_t1t2_R_starstar,tau_qq_L_star,tau_qq_R_star,pTot_star,E_L_star,E_R_star,E_L_starstar,E_R_starstar,p_face,tau_qq_face,u_n_face,u_t_face,G_hat,rho_hat,tau_nn_hat,tau_nt_hat,tau_tt_hat,tau_qq_hat,tau_nt2_hat,tau_t2t2_hat,tau_t1t2_hat,'
                #:set _hlld_p4 = 'alpha_hat,alpha_rho_hat,tau_e_hat,pres_hat,blkmod1_hat,blkmod2_hat,K_hat,C_hat_1,C_hat_2,Sigma_L,Sigma_R,dSigma,Sigma_ref,a_L_ref,a_R_ref,a_ref,du_t,dtau_nt,du_t2,dtau_nt2,sensor_ptot,sensor_vt,sensor_tnt,sensor_combined,phi,alpha_L_sum,alpha_R_sum]'
                $:GPU_PARALLEL_LOOP(collapse=3, private=_hlld_p1 + _hlld_p2 + _hlld_p3 + _hlld_p4)
                do l = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = is1%beg, is1%end
                            ! Extract left/right primitive states

                            do i = 1, eqn_idx%cont%end
                                alpha_rho_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                alpha_rho_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)
                            end do

                            ! IMP: vel%L(1:3) has (3) uninitiated for 2D
                            vel%L = 0._wp
                            vel%R = 0._wp

                            ! NOTE: unlike HLL & HLLC, vel%L here is permutated by dir_idx for simpler logic
                            do i = 1, num_vels
                                ! Don't permutate here; permutate u <-> v later at u_n_L = vel%L(1)
                                vel%L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%cont%end + i)
                                vel%R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%cont%end + i)
                            end do

                            vel_rms%L = vel%L(1)**2 + vel%L(2)**2 + vel%L(3)**2
                            vel_rms%R = vel%R(1)**2 + vel%R(2)**2 + vel%R(3)**2

                            do i = 1, num_fluids
                                alpha_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%E + i)
                                alpha_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%E + i)
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

                            pres%L = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%E)
                            pres%R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%E)

                            ! Hypoelasticity
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, eqn_idx%stress%end - eqn_idx%stress%beg + 1
                                tau_e_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, eqn_idx%stress%beg - 1 + i)
                                tau_e_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, eqn_idx%stress%beg - 1 + i)
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

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, eqn_idx%stress%end - eqn_idx%stress%beg + 1
                                ! Elastic energy (guard skips when G near zero)
                                if (.not. hypo_energy_guard .or. ((G_L > verysmall) .and. (G_R > verysmall))) then
                                    E%L = E%L + (tau_e_L(i)*tau_e_L(i))/max(4._wp*G_L, verysmall)
                                    E%R = E%R + (tau_e_R(i)*tau_e_R(i))/max(4._wp*G_R, verysmall)
                                    ! Shear terms doubled: 2D/2D-axisym i==2 only; 3D i==2,4,5
                                    if ((n > 0 .and. p == 0 .and. i == 2) .or. (p > 0 .and. (i == 2 .or. i == 4 .or. i == 5))) then
                                        E%L = E%L + (tau_e_L(i)*tau_e_L(i))/max(4._wp*G_L, verysmall)
                                        E%R = E%R + (tau_e_R(i)*tau_e_R(i))/max(4._wp*G_R, verysmall)
                                    end if
                                end if
                            end do

                            H%L = (E%L + pres%L)/rho%L
                            H%R = (E%R + pres%R)/rho%R

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

                                rhoL_star = rho%L*(S_L - u_n_L)/(S_L - S_M + verysmall)
                                rhoR_star = rho%R*(S_R - u_n_R)/(S_R - S_M + verysmall)
                                fac_L = (S_L - u_n_L)/(S_L - S_M + verysmall)
                                fac_R = (S_R - u_n_R)/(S_R - S_M + verysmall)

                                E_L_star = (E%L*(u_n_L - S_L) + u_n_L*pTot_L - S_M*pTot_star)/(S_M - S_L)
                                E_R_star = (E%R*(u_n_R - S_R) + u_n_R*pTot_R - S_M*pTot_star)/(S_M - S_R)

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
                                            F_hlld(i) = (S_R*F_L(i) - S_L*F_R(i) + S_L*S_R*(U_R(i) - U_L(i)))/(S_R - S_L &
                                                   & + verysmall)
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

                                    tau_nn_L_star = tau_nn_L - (rho_hat*(G_hat*4._wp/3._wp + tau_nn_hat)*(u_n_L - S_M)) &
                                                                & /(rho%L*(u_n_L - S_L))
                                    tau_nn_R_star = tau_nn_R - (rho_hat*(G_hat*4._wp/3._wp + tau_nn_hat)*(u_n_R - S_M)) &
                                                                & /(rho%R*(u_n_R - S_R))

                                    tau_tt_L_star = tau_tt_L + (rho_hat*(G_hat*2._wp/3._wp + tau_tt_hat)*(u_n_L - S_M)) &
                                                                & /(rho%L*(u_n_L - S_L))
                                    tau_tt_R_star = tau_tt_R + (rho_hat*(G_hat*2._wp/3._wp + tau_tt_hat)*(u_n_R - S_M)) &
                                                                & /(rho%R*(u_n_R - S_R))

                                    tau_t2t2_L_star = tau_t2t2_L + (rho_hat*(G_hat*2._wp/3._wp + tau_t2t2_hat)*(u_n_L - S_M)) &
                                                                    & /(rho%L*(u_n_L - S_L))
                                    tau_t2t2_R_star = tau_t2t2_R + (rho_hat*(G_hat*2._wp/3._wp + tau_t2t2_hat)*(u_n_R - S_M)) &
                                                                    & /(rho%R*(u_n_R - S_R))

                                    tau_t1t2_L_star = tau_t1t2_L + (rho_hat*tau_t1t2_hat*(u_n_L - S_M))/(rho%L*(u_n_L - S_L))
                                    tau_t1t2_R_star = tau_t1t2_R + (rho_hat*tau_t1t2_hat*(u_n_R - S_M))/(rho%R*(u_n_R - S_R))

                                    tau_qq_L_star = tau_qq_L + (rho_hat*(G_hat*2._wp/3._wp + tau_qq_hat)*(u_n_L - S_M)) &
                                                                & /(rho%L*(u_n_L - S_L))
                                    tau_qq_R_star = tau_qq_R + (rho_hat*(G_hat*2._wp/3._wp + tau_qq_hat)*(u_n_R - S_M)) &
                                                                & /(rho%R*(u_n_R - S_R))

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

                                    alpha1_L_star = (alpha_L(1)*(S_L - u_n_L) - C_hat_1*(S_M - u_n_L))/(S_L - S_M)
                                    alpha1_R_star = (alpha_R(1)*(S_R - u_n_R) - C_hat_1*(S_M - u_n_R))/(S_R - S_M)

                                    alpha2_L_star = (alpha_L(2)*(S_L - u_n_L) - C_hat_2*(S_M - u_n_L))/(S_L - S_M)
                                    alpha2_R_star = (alpha_R(2)*(S_R - u_n_R) - C_hat_2*(S_M - u_n_R))/(S_R - S_M)

                                    ! Compute U

                                    if (p > 0 .and. .not. cyl_coord) then
                                        ! 3D Cartesian: 14-state
                                        U_starL(1) = U_L(1)*fac_L
                                        U_starL(2) = U_L(2)*fac_L
                                        U_starL(3) = rhoL_star*S_M
                                        U_starL(4) = rhoL_star*u_t_L
                                        U_starL(5) = rhoL_star*u_t2_L
                                        U_starL(6) = E_L_star
                                        U_starL(7) = alpha1_L_star
                                        U_starL(8) = alpha2_L_star
                                        U_starL(9) = rhoL_star*tau_nn_L_star
                                        U_starL(10) = rhoL_star*tau_nt_L
                                        U_starL(11) = rhoL_star*tau_nt2_L
                                        U_starL(12) = rhoL_star*tau_tt_L_star
                                        U_starL(13) = rhoL_star*tau_t2t2_L_star
                                        U_starL(14) = rhoL_star*tau_t1t2_L_star

                                        U_starR(1) = U_R(1)*fac_R
                                        U_starR(2) = U_R(2)*fac_R
                                        U_starR(3) = rhoR_star*S_M
                                        U_starR(4) = rhoR_star*u_t_R
                                        U_starR(5) = rhoR_star*u_t2_R
                                        U_starR(6) = E_R_star
                                        U_starR(7) = alpha1_R_star
                                        U_starR(8) = alpha2_R_star
                                        U_starR(9) = rhoR_star*tau_nn_R_star
                                        U_starR(10) = rhoR_star*tau_nt_R
                                        U_starR(11) = rhoR_star*tau_nt2_R
                                        U_starR(12) = rhoR_star*tau_tt_R_star
                                        U_starR(13) = rhoR_star*tau_t2t2_R_star
                                        U_starR(14) = rhoR_star*tau_t1t2_R_star

                                        U_starstarL(1) = U_L(1)*fac_L
                                        U_starstarL(2) = U_L(2)*fac_L
                                        U_starstarL(3) = rhoL_star*S_M
                                        U_starstarL(4) = rhoL_star*u_t_star
                                        U_starstarL(5) = rhoL_star*u_t2_star
                                        U_starstarL(6) = E_L_starstar
                                        U_starstarL(7) = alpha1_L_star
                                        U_starstarL(8) = alpha2_L_star
                                        U_starstarL(9) = rhoL_star*tau_nn_L_star
                                        U_starstarL(10) = rhoL_star*tau_nt_star
                                        U_starstarL(11) = rhoL_star*tau_nt2_star
                                        U_starstarL(12) = rhoL_star*tau_tt_L_starstar
                                        U_starstarL(13) = rhoL_star*tau_t2t2_L_starstar
                                        U_starstarL(14) = rhoL_star*tau_t1t2_L_starstar

                                        U_starstarR(1) = U_R(1)*fac_R
                                        U_starstarR(2) = U_R(2)*fac_R
                                        U_starstarR(3) = rhoR_star*S_M
                                        U_starstarR(4) = rhoR_star*u_t_star
                                        U_starstarR(5) = rhoR_star*u_t2_star
                                        U_starstarR(6) = E_R_starstar
                                        U_starstarR(7) = alpha1_R_star
                                        U_starstarR(8) = alpha2_R_star
                                        U_starstarR(9) = rhoR_star*tau_nn_R_star
                                        U_starstarR(10) = rhoR_star*tau_nt_star
                                        U_starstarR(11) = rhoR_star*tau_nt2_star
                                        U_starstarR(12) = rhoR_star*tau_tt_R_starstar
                                        U_starstarR(13) = rhoR_star*tau_t2t2_R_starstar
                                        U_starstarR(14) = rhoR_star*tau_t1t2_R_starstar
                                    else
                                        ! 2D/axisym: 11-state (unchanged)
                                        U_starL(1) = U_L(1)*fac_L
                                        U_starL(2) = U_L(2)*fac_L
                                        U_starL(3) = rhoL_star*S_M
                                        U_starL(4) = rhoL_star*u_t_L
                                        U_starL(5) = E_L_star
                                        U_starL(6) = alpha1_L_star
                                        U_starL(7) = alpha2_L_star
                                        U_starL(8) = rhoL_star*tau_nn_L_star
                                        U_starL(9) = rhoL_star*tau_nt_L
                                        U_starL(10) = rhoL_star*tau_tt_L_star
                                        U_starL(11) = rhoL_star*tau_qq_L_star

                                        U_starR(1) = U_R(1)*fac_R
                                        U_starR(2) = U_R(2)*fac_R
                                        U_starR(3) = rhoR_star*S_M
                                        U_starR(4) = rhoR_star*u_t_R
                                        U_starR(5) = E_R_star
                                        U_starR(6) = alpha1_R_star
                                        U_starR(7) = alpha2_R_star
                                        U_starR(8) = rhoR_star*tau_nn_R_star
                                        U_starR(9) = rhoR_star*tau_nt_R
                                        U_starR(10) = rhoR_star*tau_tt_R_star
                                        U_starR(11) = rhoR_star*tau_qq_R_star

                                        U_starstarL(1) = U_L(1)*fac_L
                                        U_starstarL(2) = U_L(2)*fac_L
                                        U_starstarL(3) = rhoL_star*S_M
                                        U_starstarL(4) = rhoL_star*u_t_star
                                        U_starstarL(5) = E_L_starstar
                                        U_starstarL(6) = alpha1_L_star
                                        U_starstarL(7) = alpha2_L_star
                                        U_starstarL(8) = rhoL_star*tau_nn_L_star
                                        U_starstarL(9) = rhoL_star*tau_nt_star
                                        U_starstarL(10) = rhoL_star*tau_tt_L_starstar
                                        U_starstarL(11) = rhoL_star*tau_qq_L_star

                                        U_starstarR(1) = U_R(1)*fac_R
                                        U_starstarR(2) = U_R(2)*fac_R
                                        U_starstarR(3) = rhoR_star*S_M
                                        U_starstarR(4) = rhoR_star*u_t_star
                                        U_starstarR(5) = E_R_starstar
                                        U_starstarR(6) = alpha1_R_star
                                        U_starstarR(7) = alpha2_R_star
                                        U_starstarR(8) = rhoR_star*tau_nn_R_star
                                        U_starstarR(9) = rhoR_star*tau_nt_star
                                        U_starstarR(10) = rhoR_star*tau_tt_R_starstar
                                        U_starstarR(11) = rhoR_star*tau_qq_R_star
                                    end if

                                    ! Compute F and select F_HLLD

                                    F_starL(1:ncomp) = F_L(1:ncomp) + S_L*(U_starL(1:ncomp) - U_L(1:ncomp))
                                    F_starR(1:ncomp) = F_R(1:ncomp) + S_R*(U_starR(1:ncomp) - U_R(1:ncomp))

                                    if (0.0_wp <= S_L) then
                                        F_hlld(1:ncomp) = F_L(1:ncomp)
                                    else if (0.0_wp <= S_Lstar) then
                                        F_hlld(1:ncomp) = F_starL(1:ncomp)
                                    else if (0.0_wp <= s_M) then
                                        F_hlld(1:ncomp) = F_starL(1:ncomp) + S_Lstar*(U_starstarL(1:ncomp) - U_starL(1:ncomp))
                                    else if (0.0_wp <= S_Rstar) then
                                        F_hlld(1:ncomp) = F_starR(1:ncomp) + S_Rstar*(U_starstarR(1:ncomp) - U_starR(1:ncomp))
                                    else if (0.0_wp <= S_R) then
                                        F_hlld(1:ncomp) = F_starR(1:ncomp)
                                    else
                                        F_hlld(1:ncomp) = F_R(1:ncomp)
                                    end if

                                    ! ADC blending (HLLD / HLL)

                                    if (riemann_hypo_ADC) then
                                        F_HLL(1:ncomp) = F_hlld(1:ncomp)
                                        if (S_L < 0._wp .and. S_R > 0._wp) then
                                            do i = 1, ncomp
                                                F_HLL(i) = (S_R*F_L(i) - S_L*F_R(i) + S_L*S_R*(U_R(i) - U_L(i)))/(S_R - S_L &
                                                      & + verysmall)
                                            end do
                                        end if

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
                                                do i = 1, ncomp
                                                    U_HLL(i) = (S_R*U_R(i) - S_L*U_L(i) - (F_R(i) - F_L(i)))/(S_R - S_L + verysmall)
                                                end do
                                                rho_HLL = U_HLL(1) + U_HLL(2)
                                                u_n_HLL_cons = U_HLL(3)/(rho_HLL + verysmall)
                                                tau_nn_HLL = U_HLL(8)/(rho_HLL + verysmall)
                                                tau_qq_face_HLL = U_HLL(11)/(rho_HLL + verysmall)
                                                p_face_HLL = F_HLL(3) - rho_HLL*u_n_HLL_cons*u_n_HLL_cons + tau_nn_HLL
                                            end if
                                        end if

                                        ! phi is anchor-independent: computed once in the shared section above
                                        do i = 1, ncomp
                                            F_hlld(i) = F_HLL(i) + phi*(F_hlld(i) - F_HLL(i))
                                        end do
                                    end if
                                end if

                                ! Reorder F_HLLD for output: pass 1 (hat_L) -> flux_rs*, pass 2 (hat_R) -> flux_hatR_rs*
                                #:for HATPASS, FLUX in [(1, 'flux_rs' + XYZ), (2, 'flux_hatR_rs' + XYZ)]
                                    if (ipass == ${HATPASS}$) then
                                        if (p > 0 .and. .not. cyl_coord) then
                                            ! 3D Cartesian: 14-state -> physical indices
                                            ${FLUX}$_vf(j, k, l, 1) = F_hlld(1)
                                            ${FLUX}$_vf(j, k, l, 2) = F_hlld(2)
                                            ${FLUX}$_vf(j, k, l, eqn_idx%cont%end + dir_idx(1)) = F_hlld(3)
                                            ${FLUX}$_vf(j, k, l, eqn_idx%cont%end + dir_idx(2)) = F_hlld(4)
                                            ${FLUX}$_vf(j, k, l, eqn_idx%cont%end + dir_idx(3)) = F_hlld(5)
                                            ${FLUX}$_vf(j, k, l, eqn_idx%E) = F_hlld(6)
                                            ${FLUX}$_vf(j, k, l, eqn_idx%E + 1) = F_hlld(7)
                                            ${FLUX}$_vf(j, k, l, eqn_idx%E + 2) = F_hlld(8)
                                            ! Map local stress to physical stress indices
                                            if (dir_idx(1) == 1) then
                                                ${FLUX}$_vf(j, k, l, eqn_idx%stress%beg) = F_hlld(9)  ! tau_nn=tau_xx
                                                ${FLUX}$_vf(j, k, l, eqn_idx%stress%beg + 1) = F_hlld(10)  ! tau_nt1=tau_xy
                                                ${FLUX}$_vf(j, k, l, eqn_idx%stress%beg + 2) = F_hlld(12)  ! tau_t1t1=tau_yy
                                                ${FLUX}$_vf(j, k, l, eqn_idx%stress%beg + 3) = F_hlld(11)  ! tau_nt2=tau_xz
                                                ${FLUX}$_vf(j, k, l, eqn_idx%stress%beg + 4) = F_hlld(14)  ! tau_t1t2=tau_yz
                                                ${FLUX}$_vf(j, k, l, eqn_idx%stress%beg + 5) = F_hlld(13)  ! tau_t2t2=tau_zz
                                            else if (dir_idx(1) == 2) then
                                                ${FLUX}$_vf(j, k, l, eqn_idx%stress%beg) = F_hlld(12)  ! tau_t1t1=tau_xx
                                                ${FLUX}$_vf(j, k, l, eqn_idx%stress%beg + 1) = F_hlld(10)  ! tau_nt1=tau_xy
                                                ${FLUX}$_vf(j, k, l, eqn_idx%stress%beg + 2) = F_hlld(9)  ! tau_nn=tau_yy
                                                ${FLUX}$_vf(j, k, l, eqn_idx%stress%beg + 3) = F_hlld(14)  ! tau_t1t2=tau_xz
                                                ${FLUX}$_vf(j, k, l, eqn_idx%stress%beg + 4) = F_hlld(11)  ! tau_nt2=tau_yz
                                                ${FLUX}$_vf(j, k, l, eqn_idx%stress%beg + 5) = F_hlld(13)  ! tau_t2t2=tau_zz
                                            else
                                                ${FLUX}$_vf(j, k, l, eqn_idx%stress%beg) = F_hlld(12)  ! tau_t1t1=tau_xx
                                                ${FLUX}$_vf(j, k, l, eqn_idx%stress%beg + 1) = F_hlld(14)  ! tau_t1t2=tau_xy
                                                ${FLUX}$_vf(j, k, l, eqn_idx%stress%beg + 2) = F_hlld(13)  ! tau_t2t2=tau_yy
                                                ${FLUX}$_vf(j, k, l, eqn_idx%stress%beg + 3) = F_hlld(10)  ! tau_nt1=tau_xz
                                                ${FLUX}$_vf(j, k, l, eqn_idx%stress%beg + 4) = F_hlld(11)  ! tau_nt2=tau_yz
                                                ${FLUX}$_vf(j, k, l, eqn_idx%stress%beg + 5) = F_hlld(9)  ! tau_nn=tau_zz
                                            end if
                                        else
                                            ! 2D/axisym: 11-state (unchanged)
                                            ${FLUX}$_vf(j, k, l, 1) = F_hlld(1)
                                            ${FLUX}$_vf(j, k, l, 2) = F_hlld(2)
                                            ${FLUX}$_vf(j, k, l, eqn_idx%cont%end + dir_idx(1)) = F_hlld(3)
                                            ${FLUX}$_vf(j, k, l, eqn_idx%cont%end + dir_idx(2)) = F_hlld(4)
                                            ${FLUX}$_vf(j, k, l, eqn_idx%E) = F_hlld(5)
                                            ${FLUX}$_vf(j, k, l, eqn_idx%E + 1) = F_hlld(6)
                                            ${FLUX}$_vf(j, k, l, eqn_idx%E + 2) = F_hlld(7)
                                            ${FLUX}$_vf(j, k, l, eqn_idx%stress%beg + 1) = F_hlld(9)
                                            if (dir_idx(1) == 1) then
                                                ${FLUX}$_vf(j, k, l, eqn_idx%stress%beg) = F_hlld(8)
                                                ${FLUX}$_vf(j, k, l, eqn_idx%stress%beg + 2) = F_hlld(10)
                                            else
                                                ${FLUX}$_vf(j, k, l, eqn_idx%stress%beg) = F_hlld(10)
                                                ${FLUX}$_vf(j, k, l, eqn_idx%stress%beg + 2) = F_hlld(8)
                                            end if
                                            if (cyl_coord) then
                                                ${FLUX}$_vf(j, k, l, eqn_idx%stress%end) = F_hlld(11)
                                            end if
                                        end if
                                    end if
                                #:endfor

                                ! Export face velocities for axisym hypo source terms
                                if (grid_geometry == 2) then
                                    if (0._wp <= S_L) then
                                        u_n_face = u_n_L; u_t_face = u_t_L
                                    else if (0._wp <= S_Lstar) then
                                        u_n_face = S_M; u_t_face = u_t_L
                                    else if (0._wp <= S_M) then
                                        u_n_face = S_M; u_t_face = u_t_star
                                    else if (0._wp <= S_Rstar) then
                                        u_n_face = S_M; u_t_face = u_t_star
                                    else if (0._wp <= S_R) then
                                        u_n_face = S_M; u_t_face = u_t_R
                                    else
                                        u_n_face = u_n_R; u_t_face = u_t_R
                                    end if
                                    ! ADC blend NC face velocities with HLL scalar traces
                                    if (riemann_hypo_ADC) then
                                        u_n_face = u_n_HLL_trace + phi*(u_n_face - u_n_HLL_trace)
                                        u_t_face = u_t_HLL_trace + phi*(u_t_face - u_t_HLL_trace)
                                    end if
                                    #:for HATPASS, NCV in [(1, 'nc_iface_vel_rs' + XYZ), (2, 'nc_iface_vel_hatR_rs' + XYZ)]
                                        if (ipass == ${HATPASS}$) then
                                            if (dir_idx(1) == 1) then
                                                ${NCV}$_vf(j, k, l, 1) = u_n_face
                                                ${NCV}$_vf(j, k, l, 2) = u_t_face
                                            else
                                                ${NCV}$_vf(j, k, l, 1) = u_t_face
                                                ${NCV}$_vf(j, k, l, 2) = u_n_face
                                            end if
                                        end if
                                    #:endfor
                                end if

                                ! Radial geometric source flux for cylindrical coordinates
                                #:if (NORM_DIR == 2)
                                    if (cyl_coord) then
                                        #:for HATPASS, GSRC, FLUX in [(1, 'flux_gsrc_rs' + XYZ, 'flux_rs' + XYZ), (2, &
                                                                       & 'flux_gsrc_hatR_rs' + XYZ, 'flux_hatR_rs' + XYZ)]
                                            if (ipass == ${HATPASS}$) then
                                                $:GPU_LOOP(parallelism='[seq]')
                                                do i = 1, sys_size
                                                    ${GSRC}$_vf(j, k, l, i) = ${FLUX}$_vf(j, k, l, i)
                                                end do
                                                ! Pure HLLD face state
                                                if (0._wp <= S_L) then
                                                    p_face = pres%L; tau_qq_face = tau_qq_L
                                                else if (0._wp <= S_Lstar .or. 0._wp <= S_M) then
                                                    p_face = pTot_star + tau_nn_L_star; tau_qq_face = tau_qq_L_star
                                                else if (0._wp <= S_Rstar .or. 0._wp <= S_R) then
                                                    p_face = pTot_star + tau_nn_R_star; tau_qq_face = tau_qq_R_star
                                                else
                                                    p_face = pres%R; tau_qq_face = tau_qq_R
                                                end if
                                                ! ADC blend face state (HLLD ADC blends all components)
                                                if (riemann_hypo_ADC) then
                                                    p_face = p_face_HLL + phi*(p_face - p_face_HLL)
                                                    tau_qq_face = tau_qq_face_HLL + phi*(tau_qq_face - tau_qq_face_HLL)
                                                end if
                                                ${GSRC}$_vf(j, k, l, eqn_idx%cont%end + dir_idx(1)) = ${FLUX}$_vf(j, k, l, &
                                                            & eqn_idx%cont%end + dir_idx(1)) - p_face + tau_qq_face
                                                $:GPU_LOOP(parallelism='[seq]')
                                                do i = eqn_idx%adv%beg, eqn_idx%adv%end
                                                    ${GSRC}$_vf(j, k, l, i) = 0._wp
                                                end do
                                            end if
                                        #:endfor
                                    end if
                                #:endif
                            end do

                            ! Dual-pass HLLD: all NC terms stay inside the Riemann flux (anchor-independent; written once)
                            flux_src_rs${XYZ}$_vf(j, k, l, eqn_idx%adv%beg) = 0._wp
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

        call s_finalize_riemann_solver(flux_vf, flux_src_vf, flux_gsrc_vf, norm_dir)

    end subroutine s_hypo_hlld_riemann_solver

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

        @:ALLOCATE(flux_rsx_vf(is1%beg:is1%end, is2%beg:is2%end, is3%beg:is3%end, 1:sys_size))
        @:ALLOCATE(flux_gsrc_rsx_vf(is1%beg:is1%end, is2%beg:is2%end, is3%beg:is3%end, 1:sys_size))
        @:ALLOCATE(flux_src_rsx_vf(is1%beg:is1%end, is2%beg:is2%end, is3%beg:is3%end, eqn_idx%adv%beg:sys_size))
        @:ALLOCATE(vel_src_rsx_vf(is1%beg:is1%end, is2%beg:is2%end, is3%beg:is3%end, 1:num_vels))
        if (qbmm) then
            @:ALLOCATE(mom_sp_rsx_vf(is1%beg:is1%end + 1, is2%beg:is2%end, is3%beg:is3%end, 1:4))
        end if

        if (viscous) then
            @:ALLOCATE(Re_avg_rsx_vf(is1%beg:is1%end, is2%beg:is2%end, is3%beg:is3%end, 1:2))
        end if

        if (use_nc_iface_vel) then
            @:ALLOCATE(nc_iface_vel_rsx_vf(is1%beg:is1%end, is2%beg:is2%end, is3%beg:is3%end, 1:num_dims))
        end if

        if (hypo_nc_dual_pass) then
            @:ALLOCATE(flux_hatR_rsx_vf(is1%beg:is1%end, is2%beg:is2%end, is3%beg:is3%end, 1:sys_size))
            if (use_nc_iface_vel) then
                @:ALLOCATE(nc_iface_vel_hatR_rsx_vf(is1%beg:is1%end, is2%beg:is2%end, is3%beg:is3%end, 1:num_dims))
            end if
            if (cyl_coord) then
                @:ALLOCATE(flux_gsrc_hatR_rsx_vf(is1%beg:is1%end, is2%beg:is2%end, is3%beg:is3%end, 1:sys_size))
            end if
        end if

        if (n == 0) return

        is1%beg = -1; is2%beg = 0; is3%beg = 0
        is1%end = n; is2%end = m; is3%end = p

        @:ALLOCATE(flux_rsy_vf(is1%beg:is1%end, is2%beg:is2%end, is3%beg:is3%end, 1:sys_size))
        @:ALLOCATE(flux_gsrc_rsy_vf(is1%beg:is1%end, is2%beg:is2%end, is3%beg:is3%end, 1:sys_size))
        @:ALLOCATE(flux_src_rsy_vf(is1%beg:is1%end, is2%beg:is2%end, is3%beg:is3%end, eqn_idx%adv%beg:sys_size))
        @:ALLOCATE(vel_src_rsy_vf(is1%beg:is1%end, is2%beg:is2%end, is3%beg:is3%end, 1:num_vels))

        if (qbmm) then
            @:ALLOCATE(mom_sp_rsy_vf(is1%beg:is1%end + 1, is2%beg:is2%end, is3%beg:is3%end, 1:4))
        end if

        if (viscous) then
            @:ALLOCATE(Re_avg_rsy_vf(is1%beg:is1%end, is2%beg:is2%end, is3%beg:is3%end, 1:2))
        end if

        if (use_nc_iface_vel) then
            @:ALLOCATE(nc_iface_vel_rsy_vf(is1%beg:is1%end, is2%beg:is2%end, is3%beg:is3%end, 1:num_dims))
        end if

        if (hypo_nc_dual_pass) then
            @:ALLOCATE(flux_hatR_rsy_vf(is1%beg:is1%end, is2%beg:is2%end, is3%beg:is3%end, 1:sys_size))
            if (use_nc_iface_vel) then
                @:ALLOCATE(nc_iface_vel_hatR_rsy_vf(is1%beg:is1%end, is2%beg:is2%end, is3%beg:is3%end, 1:num_dims))
            end if
            if (cyl_coord) then
                @:ALLOCATE(flux_gsrc_hatR_rsy_vf(is1%beg:is1%end, is2%beg:is2%end, is3%beg:is3%end, 1:sys_size))
            end if
        end if

        if (p == 0) return

        is1%beg = -1; is2%beg = 0; is3%beg = 0
        is1%end = p; is2%end = n; is3%end = m

        @:ALLOCATE(flux_rsz_vf(is1%beg:is1%end, is2%beg:is2%end, is3%beg:is3%end, 1:sys_size))
        @:ALLOCATE(flux_gsrc_rsz_vf(is1%beg:is1%end, is2%beg:is2%end, is3%beg:is3%end, 1:sys_size))
        @:ALLOCATE(flux_src_rsz_vf(is1%beg:is1%end, is2%beg:is2%end, is3%beg:is3%end, eqn_idx%adv%beg:sys_size))
        @:ALLOCATE(vel_src_rsz_vf(is1%beg:is1%end, is2%beg:is2%end, is3%beg:is3%end, 1:num_vels))

        if (qbmm) then
            @:ALLOCATE(mom_sp_rsz_vf(is1%beg:is1%end + 1, is2%beg:is2%end, is3%beg:is3%end, 1:4))
        end if

        if (viscous) then
            @:ALLOCATE(Re_avg_rsz_vf(is1%beg:is1%end, is2%beg:is2%end, is3%beg:is3%end, 1:2))
        end if

        if (use_nc_iface_vel) then
            @:ALLOCATE(nc_iface_vel_rsz_vf(is1%beg:is1%end, is2%beg:is2%end, is3%beg:is3%end, 1:num_dims))
        end if

        if (hypo_nc_dual_pass) then
            @:ALLOCATE(flux_hatR_rsz_vf(is1%beg:is1%end, is2%beg:is2%end, is3%beg:is3%end, 1:sys_size))
            if (use_nc_iface_vel) then
                @:ALLOCATE(nc_iface_vel_hatR_rsz_vf(is1%beg:is1%end, is2%beg:is2%end, is3%beg:is3%end, 1:num_dims))
            end if
            if (cyl_coord) then
                @:ALLOCATE(flux_gsrc_hatR_rsz_vf(is1%beg:is1%end, is2%beg:is2%end, is3%beg:is3%end, 1:sys_size))
            end if
        end if

    end subroutine s_initialize_riemann_solvers_module

    !> Populate the left and right Riemann state variable buffers based on boundary conditions
    subroutine s_populate_riemann_states_variables_buffers(qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, dqL_prim_dx_vf, &

        & dqL_prim_dy_vf, dqL_prim_dz_vf, qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, dqR_prim_dx_vf, dqR_prim_dy_vf, &
            & dqR_prim_dz_vf, norm_dir, ix, iy, iz)

        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout) :: qL_prim_rsx_vf, qL_prim_rsy_vf, &
             & qL_prim_rsz_vf, qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf

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
            ! dir_idx_tau(1:3) = (nn, nt, nt2): face-normal stress row for wave speeds and momentum flux. stress_perm(1:n_stress) =
            ! full tensor permutation mapping F_HLL local basis index -> physical storage index. Local order: (nn, nt, tt, nt2,
            ! t1t2, t2t2). In 2D only entries 1-3 are used.
            if (norm_dir == 1) then
                dir_idx_tau = (/1, 2, 4/)
                stress_perm = (/1, 2, 3, 4, 5, 6/)
            else if (norm_dir == 2) then
                dir_idx_tau = (/3, 2, 5/)
                stress_perm = (/3, 2, 1, 5, 4, 6/)
            else
                dir_idx_tau = (/6, 4, 5/)
                stress_perm = (/6, 4, 1, 5, 2, 3/)
            end if
        end if

        isx = ix; isy = iy; isz = iz
        ! for stuff in the same module
        $:GPU_UPDATE(device='[isx, isy, isz]')
        ! for stuff in different modules
        $:GPU_UPDATE(device='[dir_idx, dir_flg, dir_idx_tau, stress_perm]')

        ! Population of Buffers in x-direction
        if (norm_dir == 1) then
            if (bc_x%beg == BC_RIEMANN_EXTRAP) then  ! Riemann state extrap. BC at beginning
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            qL_prim_rsx_vf(-1, k, l, i) = qR_prim_rsx_vf(0, k, l, i)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                if (viscous .or. dummy) then
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

            if (bc_x%end == BC_RIEMANN_EXTRAP) then  ! Riemann state extrap. BC at end

                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            qR_prim_rsx_vf(m + 1, k, l, i) = qL_prim_rsx_vf(m, k, l, i)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                if (viscous .or. dummy) then
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
            if (bc_y%beg == BC_RIEMANN_EXTRAP) then  ! Riemann state extrap. BC at beginning
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            qL_prim_rsy_vf(-1, k, l, i) = qR_prim_rsy_vf(0, k, l, i)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                if (viscous .or. dummy) then
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

            if (bc_y%end == BC_RIEMANN_EXTRAP) then  ! Riemann state extrap. BC at end

                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            qR_prim_rsy_vf(n + 1, k, l, i) = qL_prim_rsy_vf(n, k, l, i)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                if (viscous .or. dummy) then
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
            if (bc_z%beg == BC_RIEMANN_EXTRAP) then  ! Riemann state extrap. BC at beginning
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            qL_prim_rsz_vf(-1, k, l, i) = qR_prim_rsz_vf(0, k, l, i)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                if (viscous .or. dummy) then
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

            if (bc_z%end == BC_RIEMANN_EXTRAP) then  ! Riemann state extrap. BC at end

                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            qR_prim_rsz_vf(p + 1, k, l, i) = qL_prim_rsz_vf(p, k, l, i)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                if (viscous .or. dummy) then
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
            if (viscous .or. (surface_tension) .or. dummy) then
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
            if (viscous .or. (surface_tension) .or. dummy) then
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
                                mom_sp_rsy_vf(j, k, l, i) = mom_sp(i)%sf(k, j, l)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            ! Reshaping Inputted Data in z-direction
        else
            if (viscous .or. (surface_tension) .or. dummy) then
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
                                mom_sp_rsz_vf(j, k, l, i) = mom_sp(i)%sf(l, k, j)
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
                        r_eff = y_cc(k)
                    case (2)  ! y-face (radial face in r_cyl direction)
                        Re_s = Re_avg_rsy_vf(k, j, l, 1)
                        Re_b = Re_avg_rsy_vf(k, j, l, 2)
                        vel_src_int = vel_src_rsy_vf(k, j, l,1:num_dims)
                        r_eff = y_cb(k)
                    case (3)  ! z-face (azimuthal face in theta_cyl direction)
                        Re_s = Re_avg_rsz_vf(l, k, j, 1)
                        Re_b = Re_avg_rsz_vf(l, k, j, 2)
                        vel_src_int = vel_src_rsz_vf(l, k, j,1:num_dims)
                        r_eff = y_cc(k)
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
                        Re_shear = Re_avg_rsy_vf(k_loop, j_loop, l_loop, 1)
                        Re_bulk = Re_avg_rsy_vf(k_loop, j_loop, l_loop, 2)
                        do i_dim = 1, num_dims
                            vel_src_at_interface(i_dim) = vel_src_rsy_vf(k_loop, j_loop, l_loop, i_dim)
                        end do
                    else
                        Re_shear = Re_avg_rsz_vf(l_loop, k_loop, j_loop, 1)
                        Re_bulk = Re_avg_rsz_vf(l_loop, k_loop, j_loop, 2)
                        do i_dim = 1, num_dims
                            vel_src_at_interface(i_dim) = vel_src_rsz_vf(l_loop, k_loop, j_loop, i_dim)
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
                            flux_vf(i)%sf(k, j, l) = flux_rsy_vf(j, k, l, i)
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
                                flux_gsrc_vf(i)%sf(k, j, l) = flux_gsrc_rsy_vf(j, k, l, i)
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
                        flux_src_vf(eqn_idx%adv%beg)%sf(k, j, l) = flux_src_rsy_vf(j, k, l, eqn_idx%adv%beg)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            ! Copy the per-fluid flux_src entries when they are structurally present. HLLD writes zeros here; these entries are kept
            ! for consistency.
            if (adv_src_alpha_iface .or. adv_src_none) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = eqn_idx%adv%beg + 1, eqn_idx%adv%end
                    do l = is3%beg, is3%end
                        do j = is1%beg, is1%end
                            do k = is2%beg, is2%end
                                flux_src_vf(i)%sf(k, j, l) = flux_src_rsy_vf(j, k, l, i)
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
                            flux_vf(i)%sf(l, k, j) = flux_rsz_vf(j, k, l, i)
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
                                flux_gsrc_vf(i)%sf(l, k, j) = flux_gsrc_rsz_vf(j, k, l, i)
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
                        flux_src_vf(eqn_idx%adv%beg)%sf(l, k, j) = flux_src_rsz_vf(j, k, l, eqn_idx%adv%beg)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            ! Copy the per-fluid flux_src entries when they are structurally present. HLLD writes zeros here; these entries are kept
            ! for consistency.
            if (adv_src_alpha_iface .or. adv_src_none) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = eqn_idx%adv%beg + 1, eqn_idx%adv%end
                    do j = is1%beg, is1%end
                        do k = is2%beg, is2%end
                            do l = is3%beg, is3%end
                                flux_src_vf(i)%sf(l, k, j) = flux_src_rsz_vf(j, k, l, i)
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

            ! Copy the per-fluid flux_src entries when they are structurally present. HLLD writes zeros here; these entries are kept
            ! for consistency.
            if (adv_src_alpha_iface .or. adv_src_none) then
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

    !> Copy hypo interface velocities from Riemann-space buffers to physical-space output arrays. Called after the Riemann solver
    !! for each sweep direction when hypo_nc_interface is active.
    !! @param nc_iface_vel_vf Output: physical velocity components at interfaces
    !! @param norm_dir Sweep direction (1=x, 2=y, 3=z)
    subroutine s_finalize_nc_iface_vel(nc_iface_vel_vf, norm_dir)

        type(scalar_field), dimension(:), intent(inout) :: nc_iface_vel_vf
        integer, intent(in)                             :: norm_dir
        integer                                         :: i, j, k, l

        if (norm_dir == 2) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, num_dims
                do l = is3%beg, is3%end
                    do j = is1%beg, is1%end
                        do k = is2%beg, is2%end
                            nc_iface_vel_vf(i)%sf(k, j, l) = nc_iface_vel_rsy_vf(j, k, l, i)
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
                            nc_iface_vel_vf(i)%sf(j, k, l) = nc_iface_vel_rsx_vf(j, k, l, i)
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
                            nc_iface_vel_vf(i)%sf(l, k, j) = nc_iface_vel_rsz_vf(j, k, l, i)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_finalize_nc_iface_vel

    !> Unpermute the hat_R-anchored flux set of the fused dual-pass HLLD solve into physical-space output arrays. Mirrors
    !! s_finalize_riemann_solver for the flux_hatR_rs* buffers; flux_src is anchor-independent (already finalized with the hat_L
    !! set) and the geometric source flux only exists for the axisymmetric y-sweep.
    !! @param flux_vf Output: hat_R fluxes (overwrites the consumed hat_L values)
    !! @param flux_gsrc_vf Output: hat_R geometric source fluxes (axisymmetric y-sweep only)
    !! @param norm_dir Sweep direction (1=x, 2=y, 3=z)
    subroutine s_finalize_riemann_solver_hatR(flux_vf, flux_gsrc_vf, norm_dir)

        type(scalar_field), dimension(sys_size), intent(inout) :: flux_vf, flux_gsrc_vf
        integer, intent(in)                                    :: norm_dir
        integer                                                :: i, j, k, l  !< Generic loop iterators

        if (norm_dir == 2) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do l = is3%beg, is3%end
                    do j = is1%beg, is1%end
                        do k = is2%beg, is2%end
                            flux_vf(i)%sf(k, j, l) = flux_hatR_rsy_vf(j, k, l, i)
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
                                flux_gsrc_vf(i)%sf(k, j, l) = flux_gsrc_hatR_rsy_vf(j, k, l, i)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        else if (norm_dir == 3) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do j = is1%beg, is1%end
                    do k = is2%beg, is2%end
                        do l = is3%beg, is3%end
                            flux_vf(i)%sf(l, k, j) = flux_hatR_rsz_vf(j, k, l, i)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        else if (norm_dir == 1) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do l = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = is1%beg, is1%end
                            flux_vf(i)%sf(j, k, l) = flux_hatR_rsx_vf(j, k, l, i)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_finalize_riemann_solver_hatR

    !> Copy the hat_R-pass hypo interface velocities of the fused dual-pass HLLD solve from Riemann-space buffers to physical-space
    !! output arrays. Mirrors s_finalize_nc_iface_vel for the nc_iface_vel_hatR_rs* buffers.
    !! @param nc_iface_vel_vf Output: physical hat_R velocity components at interfaces
    !! @param norm_dir Sweep direction (1=x, 2=y, 3=z)
    subroutine s_finalize_nc_iface_vel_hatR(nc_iface_vel_vf, norm_dir)

        type(scalar_field), dimension(:), intent(inout) :: nc_iface_vel_vf
        integer, intent(in)                             :: norm_dir
        integer                                         :: i, j, k, l

        if (norm_dir == 2) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, num_dims
                do l = is3%beg, is3%end
                    do j = is1%beg, is1%end
                        do k = is2%beg, is2%end
                            nc_iface_vel_vf(i)%sf(k, j, l) = nc_iface_vel_hatR_rsy_vf(j, k, l, i)
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
                            nc_iface_vel_vf(i)%sf(j, k, l) = nc_iface_vel_hatR_rsx_vf(j, k, l, i)
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
                            nc_iface_vel_vf(i)%sf(l, k, j) = nc_iface_vel_hatR_rsz_vf(j, k, l, i)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_finalize_nc_iface_vel_hatR

    !> Module deallocation and/or disassociation procedures
    impure subroutine s_finalize_riemann_solvers_module

        if (viscous) then
            @:DEALLOCATE(Re_avg_rsx_vf)
        end if
        @:DEALLOCATE(vel_src_rsx_vf)
        @:DEALLOCATE(flux_rsx_vf)
        @:DEALLOCATE(flux_src_rsx_vf)
        @:DEALLOCATE(flux_gsrc_rsx_vf)
        if (use_nc_iface_vel) then
            @:DEALLOCATE(nc_iface_vel_rsx_vf)
        end if
        if (qbmm) then
            @:DEALLOCATE(mom_sp_rsx_vf)
        end if
        if (hypo_nc_dual_pass) then
            @:DEALLOCATE(flux_hatR_rsx_vf)
            if (use_nc_iface_vel) then
                @:DEALLOCATE(nc_iface_vel_hatR_rsx_vf)
            end if
            if (cyl_coord) then
                @:DEALLOCATE(flux_gsrc_hatR_rsx_vf)
            end if
        end if

        if (n == 0) return

        if (viscous) then
            @:DEALLOCATE(Re_avg_rsy_vf)
        end if
        @:DEALLOCATE(vel_src_rsy_vf)
        @:DEALLOCATE(flux_rsy_vf)
        @:DEALLOCATE(flux_src_rsy_vf)
        @:DEALLOCATE(flux_gsrc_rsy_vf)
        if (use_nc_iface_vel) then
            @:DEALLOCATE(nc_iface_vel_rsy_vf)
        end if
        if (qbmm) then
            @:DEALLOCATE(mom_sp_rsy_vf)
        end if
        if (hypo_nc_dual_pass) then
            @:DEALLOCATE(flux_hatR_rsy_vf)
            if (use_nc_iface_vel) then
                @:DEALLOCATE(nc_iface_vel_hatR_rsy_vf)
            end if
            if (cyl_coord) then
                @:DEALLOCATE(flux_gsrc_hatR_rsy_vf)
            end if
        end if

        if (p == 0) return

        if (viscous) then
            @:DEALLOCATE(Re_avg_rsz_vf)
        end if
        @:DEALLOCATE(vel_src_rsz_vf)
        @:DEALLOCATE(flux_rsz_vf)
        @:DEALLOCATE(flux_src_rsz_vf)
        @:DEALLOCATE(flux_gsrc_rsz_vf)
        if (use_nc_iface_vel) then
            @:DEALLOCATE(nc_iface_vel_rsz_vf)
        end if
        if (qbmm) then
            @:DEALLOCATE(mom_sp_rsz_vf)
        end if
        if (hypo_nc_dual_pass) then
            @:DEALLOCATE(flux_hatR_rsz_vf)
            if (use_nc_iface_vel) then
                @:DEALLOCATE(nc_iface_vel_hatR_rsz_vf)
            end if
            if (cyl_coord) then
                @:DEALLOCATE(flux_gsrc_hatR_rsz_vf)
            end if
        end if

    end subroutine s_finalize_riemann_solvers_module

end module m_riemann_solvers
