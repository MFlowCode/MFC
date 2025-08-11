!>
!! @file m_riemann_solvers.f90
!! @brief Contains module m_riemann_solvers

!> @brief This module features a database of approximate and exact Riemann
!!              problem solvers for the Navier-Stokes system of equations, which
!!              is supplemented by appropriate advection equations that are used
!!              to capture the material interfaces. The closure of the system is
!!              achieved by the stiffened gas equation of state and any required
!!              mixture relations. Surface tension effects are accounted for and
!!              are modeled by means of a volume force acting across the diffuse
!!              material interface region. The implementation details of viscous
!!              and capillary effects, into the Riemann solvers, may be found in
!!              Perigaud and Saurel (2005). Note that both effects are available
!!              only in the volume fraction model. At this time, the approximate
!!              and exact Riemann solvers that are listed below are available:
!!                  1) Harten-Lax-van Leer (HLL)
!!                  2) Harten-Lax-van Leer-Contact (HLLC)
!!                  3) Exact
!!                  4) Harten-Lax-van Leer Discontinuities (HLLD) - for MHD only

#:include 'case.fpp'
#:include 'macros.fpp'
#:include 'inline_riemann.fpp'

module m_riemann_solvers

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    use m_bubbles              !< To get the bubble wall pressure function

    use m_bubbles_EE

    use m_surface_tension      !< To get the capillary fluxes

    use m_helper_basic         !< Functions to compare floating point numbers

    use m_chemistry

    use m_thermochem, only: &
        gas_constant, get_mixture_molecular_weight, &
        get_mixture_specific_heat_cv_mass, get_mixture_energy_mass, &
        get_species_specific_heats_r, get_species_enthalpies_rt, &
        get_mixture_specific_heat_cp_mass

    implicit none

    private; public :: s_initialize_riemann_solvers_module, &
 s_riemann_solver, &
 s_hll_riemann_solver, &
 s_hllc_riemann_solver, &
 s_hlld_riemann_solver, &
 s_finalize_riemann_solvers_module

    !> The cell-boundary values of the fluxes (src - source) that are computed
    !! through the chosen Riemann problem solver, and the direct evaluation of
    !! source terms, by using the left and right states given in qK_prim_rs_vf,
    !! dqK_prim_ds_vf where ds = dx, dy or dz.
    !> @{

    real(wp), allocatable, dimension(:, :, :, :) :: flux_rsx_vf, flux_src_rsx_vf
    real(wp), allocatable, dimension(:, :, :, :) :: flux_rsy_vf, flux_src_rsy_vf
    real(wp), allocatable, dimension(:, :, :, :) :: flux_rsz_vf, flux_src_rsz_vf
    $:GPU_DECLARE(create='[flux_rsx_vf,flux_src_rsx_vf,flux_rsy_vf,flux_src_rsy_vf,flux_rsz_vf,flux_src_rsz_vf]')
    !> @}

    !> The cell-boundary values of the geometrical source flux that are computed
    !! through the chosen Riemann problem solver by using the left and right
    !! states given in qK_prim_rs_vf. Currently 2D axisymmetric for inviscid only.
    !> @{

    real(wp), allocatable, dimension(:, :, :, :) :: flux_gsrc_rsx_vf !<
    real(wp), allocatable, dimension(:, :, :, :) :: flux_gsrc_rsy_vf !<
    real(wp), allocatable, dimension(:, :, :, :) :: flux_gsrc_rsz_vf !<
    $:GPU_DECLARE(create='[flux_gsrc_rsx_vf,flux_gsrc_rsy_vf,flux_gsrc_rsz_vf]')
    !> @}

    ! The cell-boundary values of the velocity. vel_src_rs_vf is determined as
    ! part of Riemann problem solution and is used to evaluate the source flux.

    real(wp), allocatable, dimension(:, :, :, :) :: vel_src_rsx_vf
    real(wp), allocatable, dimension(:, :, :, :) :: vel_src_rsy_vf
    real(wp), allocatable, dimension(:, :, :, :) :: vel_src_rsz_vf
    $:GPU_DECLARE(create='[vel_src_rsx_vf,vel_src_rsy_vf,vel_src_rsz_vf]')

    real(wp), allocatable, dimension(:, :, :, :) :: mom_sp_rsx_vf
    real(wp), allocatable, dimension(:, :, :, :) :: mom_sp_rsy_vf
    real(wp), allocatable, dimension(:, :, :, :) :: mom_sp_rsz_vf
    $:GPU_DECLARE(create='[mom_sp_rsx_vf,mom_sp_rsy_vf,mom_sp_rsz_vf]')

    real(wp), allocatable, dimension(:, :, :, :) :: Re_avg_rsx_vf
    real(wp), allocatable, dimension(:, :, :, :) :: Re_avg_rsy_vf
    real(wp), allocatable, dimension(:, :, :, :) :: Re_avg_rsz_vf
    $:GPU_DECLARE(create='[Re_avg_rsx_vf,Re_avg_rsy_vf,Re_avg_rsz_vf]')

    !> @name Indical bounds in the s1-, s2- and s3-directions
    !> @{
    type(int_bounds_info) :: is1, is2, is3
    type(int_bounds_info) :: isx, isy, isz
    !> @}

    $:GPU_DECLARE(create='[is1,is2,is3,isx,isy,isz]')

    real(wp), allocatable, dimension(:) :: Gs
    $:GPU_DECLARE(create='[Gs]')

    real(wp), allocatable, dimension(:, :) :: Res
    $:GPU_DECLARE(create='[Res]')

contains

    !> Dispatch to the subroutines that are utilized to compute the
        !! Riemann problem solution. For additional information please reference:
        !!                        1) s_hll_riemann_solver
        !!                        2) s_hllc_riemann_solver
        !!                        3) s_exact_riemann_solver
        !!                        4) s_hlld_riemann_solver
        !!  @param qL_prim_vf The  left WENO-reconstructed cell-boundary values of the
        !!      cell-average primitive variables
        !!  @param qR_prim_vf The right WENO-reconstructed cell-boundary values of the
        !!      cell-average primitive variables
        !!  @param dqL_prim_dx_vf The  left WENO-reconstructed cell-boundary values of the
        !!      first-order x-dir spatial derivatives
        !!  @param dqL_prim_dy_vf The  left WENO-reconstructed cell-boundary values of the
        !!      first-order y-dir spatial derivatives
        !!  @param dqL_prim_dz_vf The  left WENO-reconstructed cell-boundary values of the
        !!      first-order z-dir spatial derivatives
        !!  @param dqR_prim_dx_vf The right WENO-reconstructed cell-boundary values of the
        !!      first-order x-dir spatial derivatives
        !!  @param dqR_prim_dy_vf The right WENO-reconstructed cell-boundary values of the
        !!      first-order y-dir spatial derivatives
        !!  @param dqR_prim_dz_vf The right WENO-reconstructed cell-boundary values of the
        !!      first-order z-dir spatial derivatives
        !!  @param gm_alphaL_vf  Left averaged gradient magnitude
        !!  @param gm_alphaR_vf Right averaged gradient magnitude
        !!  @param flux_vf Intra-cell fluxes
        !!  @param flux_src_vf Intra-cell fluxes sources
        !!  @param flux_gsrc_vf Intra-cell geometric fluxes sources
        !!  @param norm_dir Dir. splitting direction
        !!  @param ix Index bounds in the x-dir
        !!  @param iy Index bounds in the y-dir
        !!  @param iz Index bounds in the z-dir
        !!  @param q_prim_vf Cell-averaged primitive variables
    subroutine s_riemann_solver(qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, dqL_prim_dx_vf, &
                                dqL_prim_dy_vf, &
                                dqL_prim_dz_vf, &
                                qL_prim_vf, &
                                qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, dqR_prim_dx_vf, &
                                dqR_prim_dy_vf, &
                                dqR_prim_dz_vf, &
                                qR_prim_vf, &
                                q_prim_vf, &
                                flux_vf, flux_src_vf, &
                                flux_gsrc_vf, &
                                norm_dir, ix, iy, iz)

        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:), intent(INOUT) :: qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf
        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf

        type(scalar_field), allocatable, dimension(:), intent(INOUT) :: qL_prim_vf, qR_prim_vf

        type(scalar_field), &
            allocatable, dimension(:), &
            intent(INOUT) :: dqL_prim_dx_vf, dqR_prim_dx_vf, &
                             dqL_prim_dy_vf, dqR_prim_dy_vf, &
                             dqL_prim_dz_vf, dqR_prim_dz_vf

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: flux_vf, flux_src_vf, flux_gsrc_vf

        integer, intent(IN) :: norm_dir

        type(int_bounds_info), intent(IN) :: ix, iy, iz

        #:for NAME, NUM in [('hll', 1), ('hllc', 2), ('hlld', 4)]
            if (riemann_solver == ${NUM}$) then
                call s_${NAME}$_riemann_solver(qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, dqL_prim_dx_vf, &
                                               dqL_prim_dy_vf, &
                                               dqL_prim_dz_vf, &
                                               qL_prim_vf, &
                                               qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, dqR_prim_dx_vf, &
                                               dqR_prim_dy_vf, &
                                               dqR_prim_dz_vf, &
                                               qR_prim_vf, &
                                               q_prim_vf, &
                                               flux_vf, flux_src_vf, &
                                               flux_gsrc_vf, &
                                               norm_dir, ix, iy, iz)
            end if
        #:endfor

    end subroutine s_riemann_solver

    !> Dispatch to the subroutines that are utilized to compute
        !! the viscous source fluxes for either Cartesian or cylindrical geometries.
        !! For more information please refer to:
        !!      1) s_compute_cartesian_viscous_source_flux
        !!      2) s_compute_cylindrical_viscous_source_flux
    pure subroutine s_compute_viscous_source_flux(velL_vf, &
                                                  dvelL_dx_vf, &
                                                  dvelL_dy_vf, &
                                                  dvelL_dz_vf, &
                                                  velR_vf, &
                                                  dvelR_dx_vf, &
                                                  dvelR_dy_vf, &
                                                  dvelR_dz_vf, &
                                                  flux_src_vf, &
                                                  norm_dir, &
                                                  ix, iy, iz)

        type(scalar_field), &
            dimension(num_vels), &
            intent(IN) :: velL_vf, velR_vf, &
                          dvelL_dx_vf, dvelR_dx_vf, &
                          dvelL_dy_vf, dvelR_dy_vf, &
                          dvelL_dz_vf, dvelR_dz_vf

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: flux_src_vf

        integer, intent(IN) :: norm_dir

        type(int_bounds_info), intent(IN) :: ix, iy, iz

        if (grid_geometry == 3) then
            call s_compute_cylindrical_viscous_source_flux(velL_vf, &
                                                           dvelL_dx_vf, &
                                                           dvelL_dy_vf, &
                                                           dvelL_dz_vf, &
                                                           velR_vf, &
                                                           dvelR_dx_vf, &
                                                           dvelR_dy_vf, &
                                                           dvelR_dz_vf, &
                                                           flux_src_vf, &
                                                           norm_dir, &
                                                           ix, iy, iz)
        else
            call s_compute_cartesian_viscous_source_flux(dvelL_dx_vf, &
                                                         dvelL_dy_vf, &
                                                         dvelL_dz_vf, &
                                                         dvelR_dx_vf, &
                                                         dvelR_dy_vf, &
                                                         dvelR_dz_vf, &
                                                         flux_src_vf, &
                                                         norm_dir)
        end if
    end subroutine s_compute_viscous_source_flux

    subroutine s_hll_riemann_solver(qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, dqL_prim_dx_vf, &
                                    dqL_prim_dy_vf, &
                                    dqL_prim_dz_vf, &
                                    qL_prim_vf, &
                                    qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, dqR_prim_dx_vf, &
                                    dqR_prim_dy_vf, &
                                    dqR_prim_dz_vf, &
                                    qR_prim_vf, &
                                    q_prim_vf, &
                                    flux_vf, flux_src_vf, &
                                    flux_gsrc_vf, &
                                    norm_dir, ix, iy, iz)

        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:), intent(inout) :: qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf

        type(scalar_field), allocatable, dimension(:), intent(inout) :: qL_prim_vf, qR_prim_vf

        type(scalar_field), &
            allocatable, dimension(:), &
            intent(inout) :: dqL_prim_dx_vf, dqR_prim_dx_vf, &
                             dqL_prim_dy_vf, dqR_prim_dy_vf, &
                             dqL_prim_dz_vf, dqR_prim_dz_vf

        ! Intercell fluxes
        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: flux_vf, flux_src_vf, flux_gsrc_vf
        real(wp) :: flux_tau_L, flux_tau_R

        integer, intent(in) :: norm_dir
        type(int_bounds_info), intent(in) :: ix, iy, iz

        real(wp), dimension(num_fluids) :: alpha_rho_L, alpha_rho_R
        real(wp) :: rho_L, rho_R
        real(wp), dimension(num_vels) :: vel_L, vel_R
        real(wp) :: pres_L, pres_R
        real(wp) :: E_L, E_R
        real(wp) :: H_L, H_R
        real(wp), dimension(num_fluids) :: alpha_L, alpha_R
        real(wp), dimension(num_species) :: Ys_L, Ys_R
        real(wp), dimension(num_species) :: Cp_iL, Cp_iR, Xs_L, Xs_R, Gamma_iL, Gamma_iR
        real(wp), dimension(num_species) :: Yi_avg, Phi_avg, h_iL, h_iR, h_avg_2
        real(wp) :: Cp_avg, Cv_avg, T_avg, eps, c_sum_Yi_Phi
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
        real(wp), dimension(6) :: tau_e_L, tau_e_R
        real(wp) :: G_L, G_R
        real(wp), dimension(2) :: Re_L, Re_R
        real(wp), dimension(3) :: xi_field_L, xi_field_R

        real(wp) :: rho_avg
        real(wp) :: H_avg
        real(wp) :: gamma_avg
        real(wp) :: c_avg

        real(wp) :: s_L, s_R, s_M, s_P, s_S
        real(wp) :: xi_M, xi_P

        real(wp) :: ptilde_L, ptilde_R
        real(wp) :: vel_L_rms, vel_R_rms, vel_avg_rms
        real(wp) :: vel_L_tmp, vel_R_tmp
        real(wp) :: Ms_L, Ms_R, pres_SL, pres_SR
        real(wp) :: alpha_L_sum, alpha_R_sum
        real(wp) :: zcoef, pcorr !< low Mach number correction

        type(riemann_states) :: c_fast, pres_mag
        type(riemann_states_vec3) :: B

        type(riemann_states) :: Ga ! Gamma (Lorentz factor)
        type(riemann_states) :: vdotB, B2
        type(riemann_states_vec3) :: b4 ! 4-magnetic field components (spatial: b4x, b4y, b4z)
        type(riemann_states_vec3) :: cm ! Conservative momentum variables

        integer :: i, j, k, l, q !< Generic loop iterators

        ! Populating the buffers of the left and right Riemann problem
        ! states variables, based on the choice of boundary conditions
        call s_populate_riemann_states_variables_buffers( &
            qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, dqL_prim_dx_vf, &
            dqL_prim_dy_vf, &
            dqL_prim_dz_vf, &
            qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, dqR_prim_dx_vf, &
            dqR_prim_dy_vf, &
            dqR_prim_dz_vf, &
            norm_dir, ix, iy, iz)

        ! Reshaping inputted data based on dimensional splitting direction
        call s_initialize_riemann_solver( &
            flux_src_vf, &
            norm_dir)
        #:for NORM_DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]

            if (norm_dir == ${NORM_DIR}$) then
                $:GPU_PARALLEL_LOOP(collapse=3, private='[alpha_rho_L, alpha_rho_R, &
                    & vel_L, vel_R, alpha_L, alpha_R, tau_e_L, tau_e_R, &
                    & G_L, G_R, Re_L, Re_R, rho_avg, h_avg, gamma_avg, &
                    & s_L, s_R, s_S, Ys_L, Ys_R, xi_field_L, xi_field_R, &
                    & Cp_iL, Cp_iR, Xs_L, Xs_R, Gamma_iL, Gamma_iR, &
                    & Yi_avg, Phi_avg, h_iL, h_iR, h_avg_2, c_fast, &
                    & pres_mag, B, Ga, vdotB, B2, b4, cm, pcorr, &
                    & zcoef, vel_L_tmp, vel_R_tmp]')
                do l = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = is1%beg, is1%end
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, contxe
                                alpha_rho_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                alpha_rho_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)
                            end do

                            vel_L_rms = 0._wp; vel_R_rms = 0._wp

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_vels
                                vel_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, contxe + i)
                                vel_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, contxe + i)
                                vel_L_rms = vel_L_rms + vel_L(i)**2._wp
                                vel_R_rms = vel_R_rms + vel_R(i)**2._wp
                            end do

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                alpha_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)
                                alpha_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)
                            end do

                            pres_L = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx)
                            pres_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx)

                            if (mhd) then
                                if (n == 0) then ! 1D: constant Bx; By, Bz as variables
                                    B%L(1) = Bx0
                                    B%R(1) = Bx0
                                    B%L(2) = qL_prim_rs${XYZ}$_vf(j, k, l, B_idx%beg)
                                    B%R(2) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, B_idx%beg)
                                    B%L(3) = qL_prim_rs${XYZ}$_vf(j, k, l, B_idx%beg + 1)
                                    B%R(3) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, B_idx%beg + 1)
                                else ! 2D/3D: Bx, By, Bz as variables
                                    B%L(1) = qL_prim_rs${XYZ}$_vf(j, k, l, B_idx%beg)
                                    B%R(1) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, B_idx%beg)
                                    B%L(2) = qL_prim_rs${XYZ}$_vf(j, k, l, B_idx%beg + 1)
                                    B%R(2) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, B_idx%beg + 1)
                                    B%L(3) = qL_prim_rs${XYZ}$_vf(j, k, l, B_idx%beg + 2)
                                    B%R(3) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, B_idx%beg + 2)
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
                                        Re_L(i) = alpha_L(Re_idx(i, q))/Res(i, q) &
                                                  + Re_L(i)
                                        Re_R(i) = alpha_R(Re_idx(i, q))/Res(i, q) &
                                                  + Re_R(i)
                                    end do

                                    Re_L(i) = 1._wp/max(Re_L(i), sgm_eps)
                                    Re_R(i) = 1._wp/max(Re_R(i), sgm_eps)
                                end do
                            end if

                            if (chemistry) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = chemxb, chemxe
                                    Ys_L(i - chemxb + 1) = qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                    Ys_R(i - chemxb + 1) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)
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
                                H_L = (E_L + pres_L)/rho_L
                                H_R = (E_R + pres_R)/rho_R
                            elseif (mhd .and. relativity) then
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
                            elseif (mhd .and. .not. relativity) then
                                pres_mag%L = 0.5_wp*(B%L(1)**2._wp + B%L(2)**2._wp + B%L(3)**2._wp)
                                pres_mag%R = 0.5_wp*(B%R(1)**2._wp + B%R(2)**2._wp + B%R(3)**2._wp)
                                E_L = gamma_L*pres_L + pi_inf_L + 0.5_wp*rho_L*vel_L_rms + qv_L + pres_mag%L
                                E_R = gamma_R*pres_R + pi_inf_R + 0.5_wp*rho_R*vel_R_rms + qv_R + pres_mag%R ! includes magnetic energy
                                H_L = (E_L + pres_L - pres_mag%L)/rho_L
                                H_R = (E_R + pres_R - pres_mag%R)/rho_R ! stagnation enthalpy here excludes magnetic energy (only used to find speed of sound)
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
                                    G_L = G_L + alpha_L(i)*Gs(i)
                                    G_R = G_R + alpha_R(i)*Gs(i)
                                end do

                                if (cont_damage) then
                                    G_L = G_L*max((1._wp - qL_prim_rs${XYZ}$_vf(j, k, l, damage_idx)), 0._wp)
                                    G_R = G_R*max((1._wp - qR_prim_rs${XYZ}$_vf(j, k, l, damage_idx)), 0._wp)
                                end if

                                do i = 1, strxe - strxb + 1
                                    tau_e_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, strxb - 1 + i)
                                    tau_e_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, strxb - 1 + i)
                                    ! Elastic contribution to energy if G large enough
                                    !TODO take out if statement if stable without
                                    if ((G_L > 1000) .and. (G_R > 1000)) then
                                        E_L = E_L + (tau_e_L(i)*tau_e_L(i))/(4._wp*G_L)
                                        E_R = E_R + (tau_e_R(i)*tau_e_R(i))/(4._wp*G_R)
                                        ! Double for shear stresses
                                        if (any(strxb - 1 + i == shear_indices)) then
                                            E_L = E_L + (tau_e_L(i)*tau_e_L(i))/(4._wp*G_L)
                                            E_R = E_R + (tau_e_R(i)*tau_e_R(i))/(4._wp*G_R)
                                        end if
                                    end if
                                end do
                            end if

                            ! elastic energy update
                            !if ( hyperelasticity ) then
                            !    G_L = 0._wp
                            !    G_R = 0._wp
                            !
                            !    $:GPU_LOOP(parallelism='[seq]')
                            !    do i = 1, num_fluids
                            !        G_L = G_L + alpha_L(i)*Gs(i)
                            !        G_R = G_R + alpha_R(i)*Gs(i)
                            !    end do
                            !    ! Elastic contribution to energy if G large enough
                            !    if ((G_L > 1.e-3_wp) .and. (G_R > 1.e-3_wp)) then
                            !    E_L = E_L + G_L*qL_prim_rs${XYZ}$_vf(j, k, l, xiend + 1)
                            !    E_R = E_R + G_R*qR_prim_rs${XYZ}$_vf(j + 1, k, l, xiend + 1)
                            !    $:GPU_LOOP(parallelism='[seq]')
                            !    do i = 1, b_size-1
                            !        tau_e_L(i) = G_L*qL_prim_rs${XYZ}$_vf(j, k, l, strxb - 1 + i)
                            !        tau_e_R(i) = G_R*qR_prim_rs${XYZ}$_vf(j + 1, k, l, strxb - 1 + i)
                            !    end do
                            !    $:GPU_LOOP(parallelism='[seq]')
                            !    do i = 1, b_size-1
                            !        tau_e_L(i) = 0._wp
                            !        tau_e_R(i) = 0._wp
                            !    end do
                            !    $:GPU_LOOP(parallelism='[seq]')
                            !    do i = 1, num_dims
                            !        xi_field_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, xibeg - 1 + i)
                            !        xi_field_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, xibeg - 1 + i)
                            !    end do
                            !    end if
                            !end if

                            @:compute_average_state()

                            call s_compute_speed_of_sound(pres_L, rho_L, gamma_L, pi_inf_L, H_L, alpha_L, &
                                                          vel_L_rms, 0._wp, c_L)

                            call s_compute_speed_of_sound(pres_R, rho_R, gamma_R, pi_inf_R, H_R, alpha_R, &
                                                          vel_R_rms, 0._wp, c_R)

                            !> The computation of c_avg does not require all the variables, and therefore the non '_avg'
                            ! variables are placeholders to call the subroutine.

                            call s_compute_speed_of_sound(pres_R, rho_avg, gamma_avg, pi_inf_R, H_avg, alpha_R, &
                                                          vel_avg_rms, c_sum_Yi_Phi, c_avg)

                            if (mhd) then
                                call s_compute_fast_magnetosonic_speed(rho_L, c_L, B%L, norm_dir, c_fast%L, H_L)
                                call s_compute_fast_magnetosonic_speed(rho_R, c_R, B%R, norm_dir, c_fast%R, H_R)
                            end if

                            if (viscous) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, 2
                                    Re_avg_rs${XYZ}$_vf(j, k, l, i) = 2._wp/(1._wp/Re_L(i) + 1._wp/Re_R(i))
                                end do
                            end if

                            if (wave_speeds == 1) then
                                if (mhd) then
                                    s_L = min(vel_L(dir_idx(1)) - c_fast%L, vel_R(dir_idx(1)) - c_fast%R)
                                    s_R = max(vel_R(dir_idx(1)) + c_fast%R, vel_L(dir_idx(1)) + c_fast%L)
                                elseif (hypoelasticity) then
                                    s_L = min(vel_L(dir_idx(1)) - sqrt(c_L*c_L + &
                                                                       (((4._wp*G_L)/3._wp) + &
                                                                        tau_e_L(dir_idx_tau(1)))/rho_L) &
                                              , vel_R(dir_idx(1)) - sqrt(c_R*c_R + &
                                                                         (((4._wp*G_R)/3._wp) + &
                                                                          tau_e_R(dir_idx_tau(1)))/rho_R))
                                    s_R = max(vel_R(dir_idx(1)) + sqrt(c_R*c_R + &
                                                                       (((4._wp*G_R)/3._wp) + &
                                                                        tau_e_R(dir_idx_tau(1)))/rho_R) &
                                              , vel_L(dir_idx(1)) + sqrt(c_L*c_L + &
                                                                         (((4._wp*G_L)/3._wp) + &
                                                                          tau_e_L(dir_idx_tau(1)))/rho_L))
                                else if (hyperelasticity) then
                                    s_L = min(vel_L(dir_idx(1)) - sqrt(c_L*c_L + (4._wp*G_L/3._wp)/rho_L) &
                                              , vel_R(dir_idx(1)) - sqrt(c_R*c_R + (4._wp*G_R/3._wp)/rho_R))
                                    s_R = max(vel_R(dir_idx(1)) + sqrt(c_R*c_R + (4._wp*G_R/3._wp)/rho_R) &
                                              , vel_L(dir_idx(1)) + sqrt(c_L*c_L + (4._wp*G_L/3._wp)/rho_L))
                                else
                                    s_L = min(vel_L(dir_idx(1)) - c_L, vel_R(dir_idx(1)) - c_R)
                                    s_R = max(vel_R(dir_idx(1)) + c_R, vel_L(dir_idx(1)) + c_L)
                                end if

                                s_S = (pres_R - pres_L + rho_L*vel_L(dir_idx(1))* &
                                       (s_L - vel_L(dir_idx(1))) - &
                                       rho_R*vel_R(dir_idx(1))* &
                                       (s_R - vel_R(dir_idx(1)))) &
                                      /(rho_L*(s_L - vel_L(dir_idx(1))) - &
                                        rho_R*(s_R - vel_R(dir_idx(1))))
                            elseif (wave_speeds == 2) then
                                pres_SL = 5.e-1_wp*(pres_L + pres_R + rho_avg*c_avg* &
                                                    (vel_L(dir_idx(1)) - &
                                                     vel_R(dir_idx(1))))

                                pres_SR = pres_SL

                                Ms_L = max(1._wp, sqrt(1._wp + ((5.e-1_wp + gamma_L)/(1._wp + gamma_L))* &
                                                       (pres_SL/pres_L - 1._wp)*pres_L/ &
                                                       ((pres_L + pi_inf_L/(1._wp + gamma_L)))))
                                Ms_R = max(1._wp, sqrt(1._wp + ((5.e-1_wp + gamma_R)/(1._wp + gamma_R))* &
                                                       (pres_SR/pres_R - 1._wp)*pres_R/ &
                                                       ((pres_R + pi_inf_R/(1._wp + gamma_R)))))

                                s_L = vel_L(dir_idx(1)) - c_L*Ms_L
                                s_R = vel_R(dir_idx(1)) + c_R*Ms_R

                                s_S = 5.e-1_wp*((vel_L(dir_idx(1)) + vel_R(dir_idx(1))) + &
                                                (pres_L - pres_R)/ &
                                                (rho_avg*c_avg))
                            end if

                            s_M = min(0._wp, s_L); s_P = max(0._wp, s_R)

                            xi_M = (5.e-1_wp + sign(5.e-1_wp, s_L)) &
                                   + (5.e-1_wp - sign(5.e-1_wp, s_L)) &
                                   *(5.e-1_wp + sign(5.e-1_wp, s_R))
                            xi_P = (5.e-1_wp - sign(5.e-1_wp, s_R)) &
                                   + (5.e-1_wp - sign(5.e-1_wp, s_L)) &
                                   *(5.e-1_wp + sign(5.e-1_wp, s_R))

                            ! Low Mach correction
                            if (low_Mach == 1) then
                                @:compute_low_Mach_correction()
                            else
                                pcorr = 0._wp
                            end if

                            ! Mass
                            if (.not. relativity) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, contxe
                                    flux_rs${XYZ}$_vf(j, k, l, i) = &
                                        (s_M*alpha_rho_R(i)*vel_R(norm_dir) &
                                         - s_P*alpha_rho_L(i)*vel_L(norm_dir) &
                                         + s_M*s_P*(alpha_rho_L(i) &
                                                    - alpha_rho_R(i))) &
                                        /(s_M - s_P)
                                end do
                            elseif (relativity) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, contxe
                                    flux_rs${XYZ}$_vf(j, k, l, i) = &
                                        (s_M*Ga%R*alpha_rho_R(i)*vel_R(norm_dir) &
                                         - s_P*Ga%L*alpha_rho_L(i)*vel_L(norm_dir) &
                                         + s_M*s_P*(Ga%L*alpha_rho_L(i) &
                                                    - Ga%R*alpha_rho_R(i))) &
                                        /(s_M - s_P)
                                end do
                            end if

                            ! Momentum
                            if (mhd .and. (.not. relativity)) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, 3
                                    ! Flux of rho*v_i in the ${XYZ}$ direction
                                    ! = rho * v_i * v_${XYZ}$ - B_i * B_${XYZ}$ + delta_(${XYZ}$,i) * p_tot
                                    flux_rs${XYZ}$_vf(j, k, l, contxe + i) = &
                                        (s_M*(rho_R*vel_R(i)*vel_R(norm_dir) &
                                              - B%R(i)*B%R(norm_dir) &
                                              + dir_flg(i)*(pres_R + pres_mag%R)) &
                                         - s_P*(rho_L*vel_L(i)*vel_L(norm_dir) &
                                                - B%L(i)*B%L(norm_dir) &
                                                + dir_flg(i)*(pres_L + pres_mag%L)) &
                                         + s_M*s_P*(rho_L*vel_L(i) - rho_R*vel_R(i))) &
                                        /(s_M - s_P)
                                end do
                            elseif (mhd .and. relativity) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, 3
                                    ! Flux of m_i in the ${XYZ}$ direction
                                    ! = m_i * v_${XYZ}$ - b_i/Gamma * B_${XYZ}$ + delta_(${XYZ}$,i) * p_tot
                                    flux_rs${XYZ}$_vf(j, k, l, contxe + i) = &
                                        (s_M*(cm%R(i)*vel_R(norm_dir) &
                                              - b4%R(i)/Ga%R*B%R(norm_dir) &
                                              + dir_flg(i)*(pres_R + pres_mag%R)) &
                                         - s_P*(cm%L(i)*vel_L(norm_dir) &
                                                - b4%L(i)/Ga%L*B%L(norm_dir) &
                                                + dir_flg(i)*(pres_L + pres_mag%L)) &
                                         + s_M*s_P*(cm%L(i) - cm%R(i))) &
                                        /(s_M - s_P)
                                end do
                            elseif (bubbles_euler) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_vels
                                    flux_rs${XYZ}$_vf(j, k, l, contxe + dir_idx(i)) = &
                                        (s_M*(rho_R*vel_R(dir_idx(1)) &
                                              *vel_R(dir_idx(i)) &
                                              + dir_flg(dir_idx(i))*(pres_R - ptilde_R)) &
                                         - s_P*(rho_L*vel_L(dir_idx(1)) &
                                                *vel_L(dir_idx(i)) &
                                                + dir_flg(dir_idx(i))*(pres_L - ptilde_L)) &
                                         + s_M*s_P*(rho_L*vel_L(dir_idx(i)) &
                                                    - rho_R*vel_R(dir_idx(i)))) &
                                        /(s_M - s_P) &
                                        + (s_M/s_L)*(s_P/s_R)*pcorr*(vel_R(dir_idx(i)) - vel_L(dir_idx(i)))
                                end do
                            else if (hypoelasticity) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_vels
                                    flux_rs${XYZ}$_vf(j, k, l, contxe + dir_idx(i)) = &
                                        (s_M*(rho_R*vel_R(dir_idx(1)) &
                                              *vel_R(dir_idx(i)) &
                                              + dir_flg(dir_idx(i))*pres_R &
                                              - tau_e_R(dir_idx_tau(i))) &
                                         - s_P*(rho_L*vel_L(dir_idx(1)) &
                                                *vel_L(dir_idx(i)) &
                                                + dir_flg(dir_idx(i))*pres_L &
                                                - tau_e_L(dir_idx_tau(i))) &
                                         + s_M*s_P*(rho_L*vel_L(dir_idx(i)) &
                                                    - rho_R*vel_R(dir_idx(i)))) &
                                        /(s_M - s_P)
                                end do
                            else
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_vels
                                    flux_rs${XYZ}$_vf(j, k, l, contxe + dir_idx(i)) = &
                                        (s_M*(rho_R*vel_R(dir_idx(1)) &
                                              *vel_R(dir_idx(i)) &
                                              + dir_flg(dir_idx(i))*pres_R) &
                                         - s_P*(rho_L*vel_L(dir_idx(1)) &
                                                *vel_L(dir_idx(i)) &
                                                + dir_flg(dir_idx(i))*pres_L) &
                                         + s_M*s_P*(rho_L*vel_L(dir_idx(i)) &
                                                    - rho_R*vel_R(dir_idx(i)))) &
                                        /(s_M - s_P) &
                                        + (s_M/s_L)*(s_P/s_R)*pcorr*(vel_R(dir_idx(i)) - vel_L(dir_idx(i)))
                                end do
                            end if

                            ! Energy
                            if (mhd .and. (.not. relativity)) then
                                ! energy flux = (E + p + p_mag) * v_${XYZ}$ - B_${XYZ}$ * (v_x*B_x + v_y*B_y + v_z*B_z)
                                flux_rs${XYZ}$_vf(j, k, l, E_idx) = &
                                    (s_M*(vel_R(norm_dir)*(E_R + pres_R + pres_mag%R) - B%R(norm_dir)*(vel_R(1)*B%R(1) + vel_R(2)*B%R(2) + vel_R(3)*B%R(3))) &
                                     - s_P*(vel_L(norm_dir)*(E_L + pres_L + pres_mag%L) - B%L(norm_dir)*(vel_L(1)*B%L(1) + vel_L(2)*B%L(2) + vel_L(3)*B%L(3))) &
                                     + s_M*s_P*(E_L - E_R)) &
                                    /(s_M - s_P)
                            elseif (mhd .and. relativity) then
                                ! energy flux = m_${XYZ}$ - mass flux
                                ! Hard-coded for single-component for now
                                flux_rs${XYZ}$_vf(j, k, l, E_idx) = &
                                    (s_M*(cm%R(norm_dir) - Ga%R*alpha_rho_R(1)*vel_R(norm_dir)) &
                                     - s_P*(cm%L(norm_dir) - Ga%L*alpha_rho_L(1)*vel_L(norm_dir)) &
                                     + s_M*s_P*(E_L - E_R)) &
                                    /(s_M - s_P)
                            else if (bubbles_euler) then
                                flux_rs${XYZ}$_vf(j, k, l, E_idx) = &
                                    (s_M*vel_R(dir_idx(1))*(E_R + pres_R - ptilde_R) &
                                     - s_P*vel_L(dir_idx(1))*(E_L + pres_L - ptilde_L) &
                                     + s_M*s_P*(E_L - E_R)) &
                                    /(s_M - s_P) &
                                    + (s_M/s_L)*(s_P/s_R)*pcorr*(vel_R_rms - vel_L_rms)/2._wp
                            else if (hypoelasticity) then
                                flux_tau_L = 0._wp; flux_tau_R = 0._wp
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    flux_tau_L = flux_tau_L + tau_e_L(dir_idx_tau(i))*vel_L(dir_idx(i))
                                    flux_tau_R = flux_tau_R + tau_e_R(dir_idx_tau(i))*vel_R(dir_idx(i))
                                end do
                                flux_rs${XYZ}$_vf(j, k, l, E_idx) = &
                                    (s_M*(vel_R(dir_idx(1))*(E_R + pres_R) - flux_tau_R) &
                                     - s_P*(vel_L(dir_idx(1))*(E_L + pres_L) - flux_tau_L) &
                                     + s_M*s_P*(E_L - E_R))/(s_M - s_P)
                            else
                                flux_rs${XYZ}$_vf(j, k, l, E_idx) = &
                                    (s_M*vel_R(dir_idx(1))*(E_R + pres_R) &
                                     - s_P*vel_L(dir_idx(1))*(E_L + pres_L) &
                                     + s_M*s_P*(E_L - E_R)) &
                                    /(s_M - s_P) &
                                    + (s_M/s_L)*(s_P/s_R)*pcorr*(vel_R_rms - vel_L_rms)/2._wp
                            end if

                            ! Elastic Stresses
                            if (hypoelasticity) then
                                do i = 1, strxe - strxb + 1 !TODO: this indexing may be slow
                                    flux_rs${XYZ}$_vf(j, k, l, strxb - 1 + i) = &
                                        (s_M*(rho_R*vel_R(dir_idx(1)) &
                                              *tau_e_R(i)) &
                                         - s_P*(rho_L*vel_L(dir_idx(1)) &
                                                *tau_e_L(i)) &
                                         + s_M*s_P*(rho_L*tau_e_L(i) &
                                                    - rho_R*tau_e_R(i))) &
                                        /(s_M - s_P)
                                end do
                            end if

                            ! Advection
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = advxb, advxe
                                flux_rs${XYZ}$_vf(j, k, l, i) = &
                                    (qL_prim_rs${XYZ}$_vf(j, k, l, i) &
                                     - qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)) &
                                    *s_M*s_P/(s_M - s_P)
                                flux_src_rs${XYZ}$_vf(j, k, l, i) = &
                                    (s_M*qR_prim_rs${XYZ}$_vf(j + 1, k, l, i) &
                                     - s_P*qL_prim_rs${XYZ}$_vf(j, k, l, i)) &
                                    /(s_M - s_P)
                            end do

                            ! Xi field
                            !if ( hyperelasticity ) then
                            !    do i = 1, num_dims
                            !      flux_rs${XYZ}$_vf(j, k, l, xibeg - 1 + i) = &
                            !        (s_M*rho_R*vel_R(dir_idx(1))*xi_field_R(i) &
                            !         - s_P*rho_L*vel_L(dir_idx(1))*xi_field_L(i) &
                            !         + s_M*s_P*(rho_L*xi_field_L(i) &
                            !                    - rho_R*xi_field_R(i))) &
                            !        /(s_M - s_P)
                            !    end do
                            !end if

                            ! Div(U)?
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_vels
                                vel_src_rs${XYZ}$_vf(j, k, l, dir_idx(i)) = &
                                    (xi_M*(rho_L*vel_L(dir_idx(i))* &
                                           (s_L - vel_L(dir_idx(1))) - &
                                           pres_L*dir_flg(dir_idx(i))) - &
                                     xi_P*(rho_R*vel_R(dir_idx(i))* &
                                           (s_R - vel_R(dir_idx(1))) - &
                                           pres_R*dir_flg(dir_idx(i)))) &
                                    /(xi_M*rho_L*(s_L - vel_L(dir_idx(1))) - &
                                      xi_P*rho_R*(s_R - vel_R(dir_idx(1))))
                            end do

                            if (bubbles_euler) then
                                ! From HLLC: Kills mass transport @ bubble gas density
                                if (num_fluids > 1) then
                                    flux_rs${XYZ}$_vf(j, k, l, contxe) = 0._wp
                                end if
                            end if

                            if (chemistry) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = chemxb, chemxe
                                    Y_L = qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                    Y_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)

                                    flux_rs${XYZ}$_vf(j, k, l, i) = (s_M*Y_R*rho_R*vel_R(dir_idx(1)) &
                                                                     - s_P*Y_L*rho_L*vel_L(dir_idx(1)) &
                                                                     + s_M*s_P*(Y_L*rho_L - Y_R*rho_R)) &
                                                                    /(s_M - s_P)
                                    flux_src_rs${XYZ}$_vf(j, k, l, i) = 0._wp
                                end do
                            end if

                            if (mhd) then
                                if (n == 0) then ! 1D: d/dx flux only & Bx = Bx0 = const.
                                    ! B_y flux = v_x * B_y - v_y * Bx0
                                    ! B_z flux = v_x * B_z - v_z * Bx0
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 0, 1
                                        flux_rsx_vf(j, k, l, B_idx%beg + i) = (s_M*(vel_R(1)*B%R(2 + i) - vel_R(2 + i)*Bx0) &
                                                                               - s_P*(vel_L(1)*B%L(2 + i) - vel_L(2 + i)*Bx0) &
                                                                               + s_M*s_P*(B%L(2 + i) - B%R(2 + i)))/(s_M - s_P)
                                    end do
                                else ! 2D/3D: Bx, By, Bz /= const. but zero flux component in the same direction
                                    ! B_x d/d${XYZ}$ flux = (1 - delta(x,${XYZ}$)) * (v_${XYZ}$ * B_x - v_x * B_${XYZ}$)
                                    ! B_y d/d${XYZ}$ flux = (1 - delta(y,${XYZ}$)) * (v_${XYZ}$ * B_y - v_y * B_${XYZ}$)
                                    ! B_z d/d${XYZ}$ flux = (1 - delta(z,${XYZ}$)) * (v_${XYZ}$ * B_z - v_z * B_${XYZ}$)
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 0, 2
                                        flux_rs${XYZ}$_vf(j, k, l, B_idx%beg + i) = (1 - dir_flg(i + 1))*( &
                                                                                    s_M*(vel_R(dir_idx(1))*B%R(i + 1) - vel_R(i + 1)*B%R(norm_dir)) - &
                                                                                    s_P*(vel_L(dir_idx(1))*B%L(i + 1) - vel_L(i + 1)*B%L(norm_dir)) + &
                                                                                    s_M*s_P*(B%L(i + 1) - B%R(i + 1)))/(s_M - s_P)
                                    end do
                                end if
                                flux_src_rs${XYZ}$_vf(j, k, l, advxb) = 0._wp
                            end if

                            #:if (NORM_DIR == 2)
                                if (cyl_coord) then
                                    !Substituting the advective flux into the inviscid geometrical source flux
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, E_idx
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = flux_rs${XYZ}$_vf(j, k, l, i)
                                    end do
                                    ! Recalculating the radial momentum geometric source flux
                                    flux_gsrc_rs${XYZ}$_vf(j, k, l, contxe + 2) = &
                                        flux_rs${XYZ}$_vf(j, k, l, contxe + 2) &
                                        - (s_M*pres_R - s_P*pres_L)/(s_M - s_P)
                                    ! Geometrical source of the void fraction(s) is zero
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = advxb, advxe
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = flux_rs${XYZ}$_vf(j, k, l, i)
                                    end do
                                end if

                                if (cyl_coord .and. hypoelasticity) then
                                    ! += tau_sigmasigma using HLL
                                    flux_gsrc_rs${XYZ}$_vf(j, k, l, contxe + 2) = &
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, contxe + 2) + &
                                        (s_M*tau_e_R(4) - s_P*tau_e_L(4)) &
                                        /(s_M - s_P)

                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = strxb, strxe
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = flux_rs${XYZ}$_vf(j, k, l, i)
                                    end do
                                end if
                            #:endif

                        end do
                    end do
                end do
            end if

        #:endfor

        if (viscous) then
            if (weno_Re_flux) then

                call s_compute_viscous_source_flux( &
                    qL_prim_vf(momxb:momxe), &
                    dqL_prim_dx_vf(momxb:momxe), &
                    dqL_prim_dy_vf(momxb:momxe), &
                    dqL_prim_dz_vf(momxb:momxe), &
                    qR_prim_vf(momxb:momxe), &
                    dqR_prim_dx_vf(momxb:momxe), &
                    dqR_prim_dy_vf(momxb:momxe), &
                    dqR_prim_dz_vf(momxb:momxe), &
                    flux_src_vf, norm_dir, ix, iy, iz)
            else
                call s_compute_viscous_source_flux( &
                    q_prim_vf(momxb:momxe), &
                    dqL_prim_dx_vf(momxb:momxe), &
                    dqL_prim_dy_vf(momxb:momxe), &
                    dqL_prim_dz_vf(momxb:momxe), &
                    q_prim_vf(momxb:momxe), &
                    dqR_prim_dx_vf(momxb:momxe), &
                    dqR_prim_dy_vf(momxb:momxe), &
                    dqR_prim_dz_vf(momxb:momxe), &
                    flux_src_vf, norm_dir, ix, iy, iz)
            end if
        end if

        call s_finalize_riemann_solver(flux_vf, flux_src_vf, &
                                       flux_gsrc_vf, &
                                       norm_dir)

    end subroutine s_hll_riemann_solver

    !> This procedure is the implementation of the Harten, Lax,
        !!      van Leer, and contact (HLLC) approximate Riemann solver,
        !!      see Toro (1999) and Johnsen (2007). The viscous and the
        !!      surface tension effects have been included by modifying
        !!      the exact Riemann solver of Perigaud and Saurel (2005).
        !!  @param qL_prim_vf The left WENO-reconstructed cell-boundary values of the
        !!      cell-average primitive variables
        !!  @param qR_prim_vf The right WENO-reconstructed cell-boundary values of the
        !!      cell-average primitive variables
        !!  @param dqL_prim_dx_vf The left WENO-reconstructed cell-boundary values of the
        !!      first-order x-dir spatial derivatives
        !!  @param dqL_prim_dy_vf The left WENO-reconstructed cell-boundary values of the
        !!      first-order y-dir spatial derivatives
        !!  @param dqL_prim_dz_vf The left WENO-reconstructed cell-boundary values of the
        !!      first-order z-dir spatial derivatives
        !!  @param dqR_prim_dx_vf The right WENO-reconstructed cell-boundary values of the
        !!      first-order x-dir spatial derivatives
        !!  @param dqR_prim_dy_vf The right WENO-reconstructed cell-boundary values of the
        !!      first-order y-dir spatial derivatives
        !!  @param dqR_prim_dz_vf The right WENO-reconstructed cell-boundary values of the
        !!      first-order z-dir spatial derivatives
        !!  @param gm_alphaL_vf Left averaged gradient magnitude
        !!  @param gm_alphaR_vf Right averaged gradient magnitude
        !!  @param flux_vf Intra-cell fluxes
        !!  @param flux_src_vf Intra-cell fluxes sources
        !!  @param flux_gsrc_vf Intra-cell geometric fluxes sources
        !!  @param norm_dir Dir. splitting direction
        !!  @param ix Index bounds in the x-dir
        !!  @param iy Index bounds in the y-dir
        !!  @param iz Index bounds in the z-dir
        !!  @param q_prim_vf Cell-averaged primitive variables
    subroutine s_hllc_riemann_solver(qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, dqL_prim_dx_vf, &
                                     dqL_prim_dy_vf, &
                                     dqL_prim_dz_vf, &
                                     qL_prim_vf, &
                                     qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, dqR_prim_dx_vf, &
                                     dqR_prim_dy_vf, &
                                     dqR_prim_dz_vf, &
                                     qR_prim_vf, &
                                     q_prim_vf, &
                                     flux_vf, flux_src_vf, &
                                     flux_gsrc_vf, &
                                     norm_dir, ix, iy, iz)

        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:), intent(inout) :: qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), allocatable, dimension(:), intent(inout) :: qL_prim_vf, qR_prim_vf

        type(scalar_field), &
            allocatable, dimension(:), &
            intent(inout) :: dqL_prim_dx_vf, dqR_prim_dx_vf, &
                             dqL_prim_dy_vf, dqR_prim_dy_vf, &
                             dqL_prim_dz_vf, dqR_prim_dz_vf

        ! Intercell fluxes
        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: flux_vf, flux_src_vf, flux_gsrc_vf

        integer, intent(in) :: norm_dir
        type(int_bounds_info), intent(in) :: ix, iy, iz

        real(wp), dimension(num_fluids) :: alpha_rho_L, alpha_rho_R
        real(wp) :: rho_L, rho_R
        real(wp), dimension(num_dims) :: vel_L, vel_R
        real(wp) :: pres_L, pres_R
        real(wp) :: E_L, E_R
        real(wp) :: H_L, H_R
        real(wp), dimension(num_fluids) :: alpha_L, alpha_R
        real(wp), dimension(num_species) :: Ys_L, Ys_R, Xs_L, Xs_R, Gamma_iL, Gamma_iR, Cp_iL, Cp_iR
        real(wp), dimension(num_species) :: Yi_avg, Phi_avg, h_iL, h_iR, h_avg_2
        real(wp) :: Cp_avg, Cv_avg, T_avg, c_sum_Yi_Phi, eps
        real(wp) :: T_L, T_R
        real(wp) :: MW_L, MW_R
        real(wp) :: R_gas_L, R_gas_R
        real(wp) :: Cp_L, Cp_R
        real(wp) :: Cv_L, Cv_R
        real(wp) :: Gamm_L, Gamm_R
        real(wp) :: Y_L, Y_R
        real(wp) :: gamma_L, gamma_R
        real(wp) :: pi_inf_L, pi_inf_R
        real(wp) :: qv_L, qv_R
        real(wp) :: c_L, c_R
        real(wp), dimension(2) :: Re_L, Re_R

        real(wp) :: rho_avg
        real(wp) :: H_avg
        real(wp) :: gamma_avg
        real(wp) :: c_avg

        real(wp) :: s_L, s_R, s_M, s_P, s_S
        real(wp) :: xi_L, xi_R !< Left and right wave speeds functions
        real(wp) :: xi_M, xi_P
        real(wp) :: xi_MP, xi_PP

        real(wp) :: nbub_L, nbub_R
        real(wp), dimension(nb) :: R0_L, R0_R
        real(wp), dimension(nb) :: V0_L, V0_R
        real(wp), dimension(nb) :: P0_L, P0_R
        real(wp), dimension(nb) :: pbw_L, pbw_R
        real(wp) :: ptilde_L, ptilde_R

        real(wp) :: alpha_L_sum, alpha_R_sum, nbub_L_denom, nbub_R_denom

        real(wp) :: PbwR3Lbar, Pbwr3Rbar
        real(wp) :: R3Lbar, R3Rbar
        real(wp) :: R3V2Lbar, R3V2Rbar

        real(wp), dimension(6) :: tau_e_L, tau_e_R
        real(wp), dimension(num_dims) :: xi_field_L, xi_field_R
        real(wp) :: G_L, G_R

        real(wp) :: vel_L_rms, vel_R_rms, vel_avg_rms
        real(wp) :: vel_L_tmp, vel_R_tmp
        real(wp) :: rho_Star, E_Star, p_Star, p_K_Star, vel_K_star
        real(wp) :: pres_SL, pres_SR, Ms_L, Ms_R
        real(wp) :: flux_ene_e
        real(wp) :: zcoef, pcorr !< low Mach number correction

        integer :: i, j, k, l, q !< Generic loop iterators
        integer :: idx1, idxi

        ! Populating the buffers of the left and right Riemann problem
        ! states variables, based on the choice of boundary conditions

        call s_populate_riemann_states_variables_buffers( &
            qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, dqL_prim_dx_vf, &
            dqL_prim_dy_vf, &
            dqL_prim_dz_vf, &
            qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, dqR_prim_dx_vf, &
            dqR_prim_dy_vf, &
            dqR_prim_dz_vf, &
            norm_dir, ix, iy, iz)

        ! Reshaping inputted data based on dimensional splitting direction

        call s_initialize_riemann_solver( &
            flux_src_vf, &
            norm_dir)

        idx1 = 1; if (dir_idx(1) == 2) idx1 = 2; if (dir_idx(1) == 3) idx1 = 3

        #:for NORM_DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]

            if (norm_dir == ${NORM_DIR}$) then

                ! 6-EQUATION MODEL WITH HLLC
                if (model_eqns == 3) then
                    !ME3
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[vel_L, vel_R, &
                        & vel_K_Star, Re_L, Re_R, rho_avg, h_avg, &
                        & gamma_avg, s_L, s_R, s_S, vel_avg_rms, &
                        & alpha_L, alpha_R, Ys_L, Ys_R, Xs_L, Xs_R, &
                        & Gamma_iL, Gamma_iR, Cp_iL, Cp_iR, Yi_avg, &
                        & Phi_avg, h_iL, h_iR, h_avg_2, tau_e_L, &
                        & tau_e_R, G_L, G_R, flux_ene_e, xi_field_L, &
                        & xi_field_R, pcorr, zcoef, vel_L_tmp, vel_R_tmp]')
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end

                                idx1 = dir_idx(1)

                                vel_L_rms = 0._wp; vel_R_rms = 0._wp

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, contxe + i)
                                    vel_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, contxe + i)
                                    vel_L_rms = vel_L_rms + vel_L(i)**2._wp
                                    vel_R_rms = vel_R_rms + vel_R(i)**2._wp
                                end do

                                pres_L = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx)
                                pres_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx)

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
                                        qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i) = min(max(0._wp, qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)), 1._wp)
                                        alpha_L_sum = alpha_L_sum + qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)
                                    end do

                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i) = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)/max(alpha_L_sum, sgm_eps)
                                    end do

                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        qR_prim_rs${XYZ}$_vf(j + 1, k, l, i) = max(0._wp, qR_prim_rs${XYZ}$_vf(j + 1, k, l, i))
                                        qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i) = min(max(0._wp, qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)), 1._wp)
                                        alpha_R_sum = alpha_R_sum + qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)
                                    end do

                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)/max(alpha_R_sum, sgm_eps)
                                    end do
                                end if

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    rho_L = rho_L + qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                    gamma_L = gamma_L + qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)*gammas(i)
                                    pi_inf_L = pi_inf_L + qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)*pi_infs(i)
                                    qv_L = qv_L + qL_prim_rs${XYZ}$_vf(j, k, l, i)*qvs(i)

                                    rho_R = rho_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)
                                    gamma_R = gamma_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)*gammas(i)
                                    pi_inf_R = pi_inf_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)*pi_infs(i)
                                    qv_R = qv_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)*qvs(i)

                                    alpha_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, advxb + i - 1)
                                    alpha_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, advxb + i - 1)
                                end do

                                if (viscous) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, 2
                                        Re_L(i) = dflt_real

                                        if (Re_size(i) > 0) Re_L(i) = 0._wp

                                        $:GPU_LOOP(parallelism='[seq]')
                                        do q = 1, Re_size(i)
                                            Re_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + Re_idx(i, q))/Res(i, q) &
                                                      + Re_L(i)
                                        end do

                                        Re_L(i) = 1._wp/max(Re_L(i), sgm_eps)

                                    end do

                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, 2
                                        Re_R(i) = dflt_real

                                        if (Re_size(i) > 0) Re_R(i) = 0._wp

                                        $:GPU_LOOP(parallelism='[seq]')
                                        do q = 1, Re_size(i)
                                            Re_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + Re_idx(i, q))/Res(i, q) &
                                                      + Re_R(i)
                                        end do

                                        Re_R(i) = 1._wp/max(Re_R(i), sgm_eps)
                                    end do
                                end if

                                E_L = gamma_L*pres_L + pi_inf_L + 5.e-1_wp*rho_L*vel_L_rms + qv_L

                                E_R = gamma_R*pres_R + pi_inf_R + 5.e-1_wp*rho_R*vel_R_rms + qv_R

                                ! ENERGY ADJUSTMENTS FOR HYPOELASTIC ENERGY
                                if (hypoelasticity) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, strxe - strxb + 1
                                        tau_e_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, strxb - 1 + i)
                                        tau_e_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, strxb - 1 + i)
                                    end do
                                    G_L = 0._wp; G_R = 0._wp
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        G_L = G_L + alpha_L(i)*Gs(i)
                                        G_R = G_R + alpha_R(i)*Gs(i)
                                    end do
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, strxe - strxb + 1
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

                                ! ENERGY ADJUSTMENTS FOR HYPERELASTIC ENERGY
                                if (hyperelasticity) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        xi_field_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, xibeg - 1 + i)
                                        xi_field_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, xibeg - 1 + i)
                                    end do
                                    G_L = 0._wp; G_R = 0._wp; 
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        ! Mixture left and right shear modulus
                                        G_L = G_L + alpha_L(i)*Gs(i)
                                        G_R = G_R + alpha_R(i)*Gs(i)
                                    end do
                                    ! Elastic contribution to energy if G large enough
                                    if (G_L > verysmall .and. G_R > verysmall) then
                                        E_L = E_L + G_L*qL_prim_rs${XYZ}$_vf(j, k, l, xiend + 1)
                                        E_R = E_R + G_R*qR_prim_rs${XYZ}$_vf(j + 1, k, l, xiend + 1)
                                    end if
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, b_size - 1
                                        tau_e_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, strxb - 1 + i)
                                        tau_e_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, strxb - 1 + i)
                                    end do
                                end if

                                H_L = (E_L + pres_L)/rho_L
                                H_R = (E_R + pres_R)/rho_R

                                @:compute_average_state()

                                call s_compute_speed_of_sound(pres_L, rho_L, gamma_L, pi_inf_L, H_L, alpha_L, &
                                                              vel_L_rms, 0._wp, c_L)

                                call s_compute_speed_of_sound(pres_R, rho_R, gamma_R, pi_inf_R, H_R, alpha_R, &
                                                              vel_R_rms, 0._wp, c_R)

                                !> The computation of c_avg does not require all the variables, and therefore the non '_avg'
                                ! variables are placeholders to call the subroutine.
                                call s_compute_speed_of_sound(pres_R, rho_avg, gamma_avg, pi_inf_R, H_avg, alpha_R, &
                                                              vel_avg_rms, 0._wp, c_avg)

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
                                        s_L = min(vel_L(dir_idx(1)) - sqrt(c_L*c_L + &
                                                                           (((4._wp*G_L)/3._wp) + tau_e_L(dir_idx_tau(1)))/rho_L), vel_R(dir_idx(1)) - sqrt(c_R*c_R + &
                                                                                                                                                            (((4._wp*G_R)/3._wp) + tau_e_R(dir_idx_tau(1)))/rho_R))
                                        s_R = max(vel_R(dir_idx(1)) + sqrt(c_R*c_R + &
                                                                           (((4._wp*G_R)/3._wp) + tau_e_R(dir_idx_tau(1)))/rho_R), vel_L(dir_idx(1)) + sqrt(c_L*c_L + &
                                                                                                                                                            (((4._wp*G_L)/3._wp) + tau_e_L(dir_idx_tau(1)))/rho_L))
                                        s_S = (pres_R - tau_e_R(dir_idx_tau(1)) - pres_L + &
                                               tau_e_L(dir_idx_tau(1)) + rho_L*vel_L(idx1)*(s_L - vel_L(idx1)) - &
                                               rho_R*vel_R(idx1)*(s_R - vel_R(idx1)))/(rho_L*(s_L - vel_L(idx1)) - &
                                                                                       rho_R*(s_R - vel_R(idx1)))
                                    else
                                        s_L = min(vel_L(dir_idx(1)) - c_L, vel_R(dir_idx(1)) - c_R)
                                        s_R = max(vel_R(dir_idx(1)) + c_R, vel_L(dir_idx(1)) + c_L)
                                        s_S = (pres_R - pres_L + rho_L*vel_L(dir_idx(1))* &
                                               (s_L - vel_L(dir_idx(1))) - rho_R*vel_R(dir_idx(1))*(s_R - vel_R(dir_idx(1)))) &
                                              /(rho_L*(s_L - vel_L(dir_idx(1))) - rho_R*(s_R - vel_R(dir_idx(1))))

                                    end if
                                elseif (wave_speeds == 2) then
                                    pres_SL = 5.e-1_wp*(pres_L + pres_R + rho_avg*c_avg* &
                                                        (vel_L(dir_idx(1)) - &
                                                         vel_R(dir_idx(1))))

                                    pres_SR = pres_SL

                                    Ms_L = max(1._wp, sqrt(1._wp + ((5.e-1_wp + gamma_L)/(1._wp + gamma_L))* &
                                                           (pres_SL/pres_L - 1._wp)*pres_L/ &
                                                           ((pres_L + pi_inf_L/(1._wp + gamma_L)))))
                                    Ms_R = max(1._wp, sqrt(1._wp + ((5.e-1_wp + gamma_R)/(1._wp + gamma_R))* &
                                                           (pres_SR/pres_R - 1._wp)*pres_R/ &
                                                           ((pres_R + pi_inf_R/(1._wp + gamma_R)))))

                                    s_L = vel_L(dir_idx(1)) - c_L*Ms_L
                                    s_R = vel_R(dir_idx(1)) + c_R*Ms_R

                                    s_S = 5.e-1_wp*((vel_L(dir_idx(1)) + vel_R(dir_idx(1))) + &
                                                    (pres_L - pres_R)/ &
                                                    (rho_avg*c_avg))
                                end if

                                ! follows Einfeldt et al.
                                ! s_M/P = min/max(0.,s_L/R)
                                s_M = min(0._wp, s_L); s_P = max(0._wp, s_R)

                                ! goes with q_star_L/R = xi_L/R * (variable)
                                ! xi_L/R = ( ( s_L/R - u_L/R )/(s_L/R - s_star) )
                                xi_L = (s_L - vel_L(idx1))/(s_L - s_S)
                                xi_R = (s_R - vel_R(idx1))/(s_R - s_S)

                                ! goes with numerical star velocity in x/y/z directions
                                ! xi_P/M = 0.5 +/m sgn(0.5,s_star)
                                xi_M = (5.e-1_wp + sign(0.5_wp, s_S))
                                xi_P = (5.e-1_wp - sign(0.5_wp, s_S))

                                ! goes with the numerical velocity in x/y/z directions
                                ! xi_P/M (pressure) = min/max(0. sgn(1,sL/sR))
                                xi_MP = -min(0._wp, sign(1._wp, s_L))
                                xi_PP = max(0._wp, sign(1._wp, s_R))

                                E_star = xi_M*(E_L + xi_MP*(xi_L*(E_L + (s_S - vel_L(dir_idx(1)))* &
                                                                  (rho_L*s_S + pres_L/(s_L - vel_L(dir_idx(1))))) - E_L)) + &
                                         xi_P*(E_R + xi_PP*(xi_R*(E_R + (s_S - vel_R(dir_idx(1)))* &
                                                                  (rho_R*s_S + pres_R/(s_R - vel_R(dir_idx(1))))) - E_R))
                                p_Star = xi_M*(pres_L + xi_MP*(rho_L*(s_L - vel_L(dir_idx(1)))*(s_S - vel_L(dir_idx(1))))) + &
                                         xi_P*(pres_R + xi_PP*(rho_R*(s_R - vel_R(dir_idx(1)))*(s_S - vel_R(dir_idx(1)))))

                                rho_Star = xi_M*(rho_L*(xi_MP*xi_L + 1._wp - xi_MP)) + &
                                           xi_P*(rho_R*(xi_PP*xi_R + 1._wp - xi_PP))

                                vel_K_Star = vel_L(idx1)*(1._wp - xi_MP) + xi_MP*vel_R(idx1) + &
                                             xi_MP*xi_PP*(s_S - vel_R(idx1))

                                ! Low Mach correction
                                if (low_Mach == 1) then
                                    @:compute_low_Mach_correction()
                                else
                                    pcorr = 0._wp
                                end if

                                ! COMPUTING FLUXES
                                ! MASS FLUX.
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, contxe
                                    flux_rs${XYZ}$_vf(j, k, l, i) = &
                                        xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, i)*(vel_L(idx1) + s_M*(xi_L - 1._wp)) + &
                                        xi_P*qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)*(vel_R(idx1) + s_P*(xi_R - 1._wp))
                                end do

                                ! MOMENTUM FLUX.
                                ! f = \rho u u - \sigma, q = \rho u, q_star = \xi * \rho*(s_star, v, w)
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    idxi = dir_idx(i)
                                    flux_rs${XYZ}$_vf(j, k, l, contxe + idxi) = rho_Star*vel_K_Star* &
                                                                                (dir_flg(idxi)*vel_K_Star + (1._wp - dir_flg(idxi))*(xi_M*vel_L(idxi) + xi_P*vel_R(idxi))) + dir_flg(idxi)*p_Star &
                                                                                + (s_M/s_L)*(s_P/s_R)*dir_flg(idxi)*pcorr
                                end do

                                ! ENERGY FLUX.
                                ! f = u*(E-\sigma), q = E, q_star = \xi*E+(s-u)(\rho s_star - \sigma/(s-u))
                                flux_rs${XYZ}$_vf(j, k, l, E_idx) = (E_star + p_Star)*vel_K_Star &
                                                                    + (s_M/s_L)*(s_P/s_R)*pcorr*s_S

                                ! ELASTICITY. Elastic shear stress additions for the momentum and energy flux
                                if (elasticity) then
                                    flux_ene_e = 0._wp; 
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        idxi = dir_idx(i)
                                        ! MOMENTUM ELASTIC FLUX.
                                        flux_rs${XYZ}$_vf(j, k, l, contxe + idxi) = &
                                            flux_rs${XYZ}$_vf(j, k, l, contxe + idxi) &
                                            - xi_M*tau_e_L(dir_idx_tau(i)) - xi_P*tau_e_R(dir_idx_tau(i))
                                        ! ENERGY ELASTIC FLUX.
                                        flux_ene_e = flux_ene_e - &
                                                     xi_M*(vel_L(idxi)*tau_e_L(dir_idx_tau(i)) + &
                                                           s_M*(xi_L*((s_S - vel_L(i))*(tau_e_L(dir_idx_tau(i))/(s_L - vel_L(i)))))) - &
                                                     xi_P*(vel_R(idxi)*tau_e_R(dir_idx_tau(i)) + &
                                                           s_P*(xi_R*((s_S - vel_R(i))*(tau_e_R(dir_idx_tau(i))/(s_R - vel_R(i))))))
                                    end do
                                    flux_rs${XYZ}$_vf(j, k, l, E_idx) = flux_rs${XYZ}$_vf(j, k, l, E_idx) + flux_ene_e
                                end if

                                ! VOLUME FRACTION FLUX.
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = advxb, advxe
                                    flux_rs${XYZ}$_vf(j, k, l, i) = &
                                        xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, i)*s_S + &
                                        xi_P*qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)*s_S
                                end do

                                ! SOURCE TERM FOR VOLUME FRACTION ADVECTION FLUX.
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    idxi = dir_idx(i)
                                    vel_src_rs${XYZ}$_vf(j, k, l, idxi) = &
                                        xi_M*(vel_L(idxi) + dir_flg(idxi)*(s_S*(xi_MP*(xi_L - 1) + 1) - vel_L(idxi))) + &
                                        xi_P*(vel_R(idxi) + dir_flg(idxi)*(s_S*(xi_PP*(xi_R - 1) + 1) - vel_R(idxi)))
                                end do

                                ! INTERNAL ENERGIES ADVECTION FLUX.
                                ! K-th pressure and velocity in preparation for the internal energy flux
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    p_K_Star = xi_M*(xi_MP*((pres_L + pi_infs(i)/(1._wp + gammas(i)))* &
                                                            xi_L**(1._wp/gammas(i) + 1._wp) - pi_infs(i)/(1._wp + gammas(i)) - pres_L) + pres_L) + &
                                               xi_P*(xi_PP*((pres_R + pi_infs(i)/(1._wp + gammas(i)))* &
                                                            xi_R**(1._wp/gammas(i) + 1._wp) - pi_infs(i)/(1._wp + gammas(i)) - pres_R) + pres_R)

                                    flux_rs${XYZ}$_vf(j, k, l, i + intxb - 1) = &
                                        ((xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, i + advxb - 1) + xi_P*qR_prim_rs${XYZ}$_vf(j + 1, k, l, i + advxb - 1))* &
                                         (gammas(i)*p_K_Star + pi_infs(i)) + &
                                         (xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, i + contxb - 1) + xi_P*qR_prim_rs${XYZ}$_vf(j + 1, k, l, i + contxb - 1))* &
                                         qvs(i))*vel_K_Star &
                                        + (s_M/s_L)*(s_P/s_R)*pcorr*s_S*(xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, i + advxb - 1) + xi_P*qR_prim_rs${XYZ}$_vf(j + 1, k, l, i + advxb - 1))
                                end do

                                flux_src_rs${XYZ}$_vf(j, k, l, advxb) = vel_src_rs${XYZ}$_vf(j, k, l, idx1)

                                ! HYPOELASTIC STRESS EVOLUTION FLUX.
                                if (hypoelasticity) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, strxe - strxb + 1
                                        flux_rs${XYZ}$_vf(j, k, l, strxb - 1 + i) = &
                                            xi_M*(s_S/(s_L - s_S))*(s_L*rho_L*tau_e_L(i) - rho_L*vel_L(idx1)*tau_e_L(i)) + &
                                            xi_P*(s_S/(s_R - s_S))*(s_R*rho_R*tau_e_R(i) - rho_R*vel_R(idx1)*tau_e_R(i))
                                    end do
                                end if

                                ! REFERENCE MAP FLUX.
                                if (hyperelasticity) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        flux_rs${XYZ}$_vf(j, k, l, xibeg - 1 + i) = &
                                            xi_M*(s_S/(s_L - s_S))*(s_L*rho_L*xi_field_L(i) &
                                                                    - rho_L*vel_L(idx1)*xi_field_L(i)) + &
                                            xi_P*(s_S/(s_R - s_S))*(s_R*rho_R*xi_field_R(i) &
                                                                    - rho_R*vel_R(idx1)*xi_field_R(i))
                                    end do
                                end if

                                ! COLOR FUNCTION FLUX
                                if (surface_tension) then
                                    flux_rs${XYZ}$_vf(j, k, l, c_idx) = &
                                        (xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, c_idx) + &
                                         xi_P*qR_prim_rs${XYZ}$_vf(j + 1, k, l, c_idx))*s_S
                                end if

                                ! Geometrical source flux for cylindrical coordinates
                                #:if (NORM_DIR == 2)
                                    if (cyl_coord) then
                                        !Substituting the advective flux into the inviscid geometrical source flux
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, E_idx
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = flux_rs${XYZ}$_vf(j, k, l, i)
                                        end do
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = intxb, intxe
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = flux_rs${XYZ}$_vf(j, k, l, i)
                                        end do
                                        ! Recalculating the radial momentum geometric source flux
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, momxb - 1 + dir_idx(1)) = &
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, momxb - 1 + dir_idx(1)) - p_Star
                                        ! Geometrical source of the void fraction(s) is zero
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = advxb, advxe
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
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, momxb - 1 + dir_idx(1)) = &
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, momxb - 1 + dir_idx(1)) - p_Star

                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, momxe) = flux_rs${XYZ}$_vf(j, k, l, momxb + 1)
                                    end if
                                #:endif

                            end do
                        end do
                    end do

                elseif (model_eqns == 4) then
                    !ME4
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[alpha_rho_L, &
                        & alpha_rho_R, vel_L, vel_R, alpha_L, alpha_R, &
                        & rho_avg, h_avg, gamma_avg, s_L, s_R, s_S, &
                        & vel_avg_rms, nbub_L, nbub_R, ptilde_L, ptilde_R]')
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, contxe
                                    alpha_rho_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                    alpha_rho_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)
                                end do

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, contxe + i)
                                    vel_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, contxe + i)
                                end do

                                vel_L_rms = 0._wp; vel_R_rms = 0._wp
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel_L_rms = vel_L_rms + vel_L(i)**2._wp
                                    vel_R_rms = vel_R_rms + vel_R(i)**2._wp
                                end do

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    alpha_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)
                                    alpha_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)
                                end do

                                pres_L = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx)
                                pres_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx)

                                rho_L = 0._wp
                                gamma_L = 0._wp
                                pi_inf_L = 0._wp
                                qv_L = 0._wp
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    rho_L = rho_L + alpha_rho_L(i)
                                    gamma_L = gamma_L + alpha_L(i)*gammas(i)
                                    pi_inf_L = pi_inf_L + alpha_L(i)*pi_infs(i)
                                    qv_L = qv_L + alpha_rho_L(i)*qvs(i)
                                end do

                                rho_R = 0._wp
                                gamma_R = 0._wp
                                pi_inf_R = 0._wp
                                qv_R = 0._wp
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    rho_R = rho_R + alpha_rho_R(i)
                                    gamma_R = gamma_R + alpha_R(i)*gammas(i)
                                    pi_inf_R = pi_inf_R + alpha_R(i)*pi_infs(i)
                                    qv_R = qv_R + alpha_rho_R(i)*qvs(i)
                                end do

                                E_L = gamma_L*pres_L + pi_inf_L + 5.e-1_wp*rho_L*vel_L_rms + qv_L

                                E_R = gamma_R*pres_R + pi_inf_R + 5.e-1_wp*rho_R*vel_R_rms + qv_R

                                H_L = (E_L + pres_L)/rho_L
                                H_R = (E_R + pres_R)/rho_R

                                @:compute_average_state()

                                call s_compute_speed_of_sound(pres_L, rho_L, gamma_L, pi_inf_L, H_L, alpha_L, &
                                                              vel_L_rms, 0._wp, c_L)

                                call s_compute_speed_of_sound(pres_R, rho_R, gamma_R, pi_inf_R, H_R, alpha_R, &
                                                              vel_R_rms, 0._wp, c_R)

                                !> The computation of c_avg does not require all the variables, and therefore the non '_avg'
                                ! variables are placeholders to call the subroutine.

                                call s_compute_speed_of_sound(pres_R, rho_avg, gamma_avg, pi_inf_R, H_avg, alpha_R, &
                                                              vel_avg_rms, 0._wp, c_avg)

                                if (wave_speeds == 1) then
                                    s_L = min(vel_L(dir_idx(1)) - c_L, vel_R(dir_idx(1)) - c_R)
                                    s_R = max(vel_R(dir_idx(1)) + c_R, vel_L(dir_idx(1)) + c_L)

                                    s_S = (pres_R - pres_L + rho_L*vel_L(dir_idx(1))* &
                                           (s_L - vel_L(dir_idx(1))) - &
                                           rho_R*vel_R(dir_idx(1))* &
                                           (s_R - vel_R(dir_idx(1)))) &
                                          /(rho_L*(s_L - vel_L(dir_idx(1))) - &
                                            rho_R*(s_R - vel_R(dir_idx(1))))
                                elseif (wave_speeds == 2) then
                                    pres_SL = 5.e-1_wp*(pres_L + pres_R + rho_avg*c_avg* &
                                                        (vel_L(dir_idx(1)) - &
                                                         vel_R(dir_idx(1))))

                                    pres_SR = pres_SL

                                    Ms_L = max(1._wp, sqrt(1._wp + ((5.e-1_wp + gamma_L)/(1._wp + gamma_L))* &
                                                           (pres_SL/pres_L - 1._wp)*pres_L/ &
                                                           ((pres_L + pi_inf_L/(1._wp + gamma_L)))))
                                    Ms_R = max(1._wp, sqrt(1._wp + ((5.e-1_wp + gamma_R)/(1._wp + gamma_R))* &
                                                           (pres_SR/pres_R - 1._wp)*pres_R/ &
                                                           ((pres_R + pi_inf_R/(1._wp + gamma_R)))))

                                    s_L = vel_L(dir_idx(1)) - c_L*Ms_L
                                    s_R = vel_R(dir_idx(1)) + c_R*Ms_R

                                    s_S = 5.e-1_wp*((vel_L(dir_idx(1)) + vel_R(dir_idx(1))) + &
                                                    (pres_L - pres_R)/ &
                                                    (rho_avg*c_avg))
                                end if

                                ! follows Einfeldt et al.
                                ! s_M/P = min/max(0.,s_L/R)
                                s_M = min(0._wp, s_L); s_P = max(0._wp, s_R)

                                ! goes with q_star_L/R = xi_L/R * (variable)
                                ! xi_L/R = ( ( s_L/R - u_L/R )/(s_L/R - s_star) )
                                xi_L = (s_L - vel_L(dir_idx(1)))/(s_L - s_S)
                                xi_R = (s_R - vel_R(dir_idx(1)))/(s_R - s_S)

                                ! goes with numerical velocity in x/y/z directions
                                ! xi_P/M = 0.5 +/m sgn(0.5,s_star)
                                xi_M = (5.e-1_wp + sign(5.e-1_wp, s_S))
                                xi_P = (5.e-1_wp - sign(5.e-1_wp, s_S))

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, contxe
                                    flux_rs${XYZ}$_vf(j, k, l, i) = &
                                        xi_M*alpha_rho_L(i) &
                                        *(vel_L(dir_idx(1)) + s_M*(xi_L - 1._wp)) &
                                        + xi_P*alpha_rho_R(i) &
                                        *(vel_R(dir_idx(1)) + s_P*(xi_R - 1._wp))
                                end do

                                ! Momentum flux.
                                ! f = \rho u u + p I, q = \rho u, q_star = \xi * \rho*(s_star, v, w)
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    flux_rs${XYZ}$_vf(j, k, l, contxe + dir_idx(i)) = &
                                        xi_M*(rho_L*(vel_L(dir_idx(1))* &
                                                     vel_L(dir_idx(i)) + &
                                                     s_M*(xi_L*(dir_flg(dir_idx(i))*s_S + &
                                                                (1._wp - dir_flg(dir_idx(i)))* &
                                                                vel_L(dir_idx(i))) - vel_L(dir_idx(i)))) + &
                                              dir_flg(dir_idx(i))*pres_L) &
                                        + xi_P*(rho_R*(vel_R(dir_idx(1))* &
                                                       vel_R(dir_idx(i)) + &
                                                       s_P*(xi_R*(dir_flg(dir_idx(i))*s_S + &
                                                                  (1._wp - dir_flg(dir_idx(i)))* &
                                                                  vel_R(dir_idx(i))) - vel_R(dir_idx(i)))) + &
                                                dir_flg(dir_idx(i))*pres_R)
                                end do

                                if (bubbles_euler) then
                                    ! Put p_tilde in
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        flux_rs${XYZ}$_vf(j, k, l, contxe + dir_idx(i)) = &
                                            flux_rs${XYZ}$_vf(j, k, l, contxe + dir_idx(i)) + &
                                            xi_M*(dir_flg(dir_idx(i))*(-1._wp*ptilde_L)) &
                                            + xi_P*(dir_flg(dir_idx(i))*(-1._wp*ptilde_R))
                                    end do
                                end if

                                flux_rs${XYZ}$_vf(j, k, l, E_idx) = 0._wp

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = alf_idx, alf_idx !only advect the void fraction
                                    flux_rs${XYZ}$_vf(j, k, l, i) = &
                                        xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, i) &
                                        *(vel_L(dir_idx(1)) + s_M*(xi_L - 1._wp)) &
                                        + xi_P*qR_prim_rs${XYZ}$_vf(j + 1, k, l, i) &
                                        *(vel_R(dir_idx(1)) + s_P*(xi_R - 1._wp))
                                end do

                                ! Source for volume fraction advection equation
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims

                                    vel_src_rs${XYZ}$_vf(j, k, l, dir_idx(i)) = 0._wp
                                    !IF ( (model_eqns == 4) .or. (num_fluids==1) ) vel_src_rs_vf(dir_idx(i))%sf(j,k,l) = 0._wp
                                end do

                                flux_src_rs${XYZ}$_vf(j, k, l, advxb) = vel_src_rs${XYZ}$_vf(j, k, l, dir_idx(1))

                                ! Add advection flux for bubble variables
                                if (bubbles_euler) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = bubxb, bubxe
                                        flux_rs${XYZ}$_vf(j, k, l, i) = &
                                            xi_M*nbub_L*qL_prim_rs${XYZ}$_vf(j, k, l, i) &
                                            *(vel_L(dir_idx(1)) + s_M*(xi_L - 1._wp)) &
                                            + xi_P*nbub_R*qR_prim_rs${XYZ}$_vf(j + 1, k, l, i) &
                                            *(vel_R(dir_idx(1)) + s_P*(xi_R - 1._wp))
                                    end do
                                end if

                                ! Geometrical source flux for cylindrical coordinates

                                #:if (NORM_DIR == 2)
                                    if (cyl_coord) then
                                        ! Substituting the advective flux into the inviscid geometrical source flux
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, E_idx
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = flux_rs${XYZ}$_vf(j, k, l, i)
                                        end do
                                        ! Recalculating the radial momentum geometric source flux
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, contxe + dir_idx(1)) = &
                                            xi_M*(rho_L*(vel_L(dir_idx(1))* &
                                                         vel_L(dir_idx(1)) + &
                                                         s_M*(xi_L*(dir_flg(dir_idx(1))*s_S + &
                                                                    (1._wp - dir_flg(dir_idx(1)))* &
                                                                    vel_L(dir_idx(1))) - vel_L(dir_idx(1))))) &
                                            + xi_P*(rho_R*(vel_R(dir_idx(1))* &
                                                           vel_R(dir_idx(1)) + &
                                                           s_P*(xi_R*(dir_flg(dir_idx(1))*s_S + &
                                                                      (1._wp - dir_flg(dir_idx(1)))* &
                                                                      vel_R(dir_idx(1))) - vel_R(dir_idx(1)))))
                                        ! Geometrical source of the void fraction(s) is zero
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = advxb, advxe
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
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, momxb + 1) = &
                                            -xi_M*(rho_L*(vel_L(dir_idx(1))* &
                                                          vel_L(dir_idx(1)) + &
                                                          s_M*(xi_L*(dir_flg(dir_idx(1))*s_S + &
                                                                     (1._wp - dir_flg(dir_idx(1)))* &
                                                                     vel_L(dir_idx(1))) - vel_L(dir_idx(1))))) &
                                            - xi_P*(rho_R*(vel_R(dir_idx(1))* &
                                                           vel_R(dir_idx(1)) + &
                                                           s_P*(xi_R*(dir_flg(dir_idx(1))*s_S + &
                                                                      (1._wp - dir_flg(dir_idx(1)))* &
                                                                      vel_R(dir_idx(1))) - vel_R(dir_idx(1)))))
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, momxe) = flux_rs${XYZ}$_vf(j, k, l, momxb + 1)
                                    end if
                                #:endif
                            end do
                        end do
                    end do

                elseif (model_eqns == 2 .and. bubbles_euler) then
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[R0_L, R0_R, V0_L, &
                        & V0_R, P0_L, P0_R, pbw_L, pbw_R, vel_L, &
                        & vel_R, rho_avg, alpha_L, alpha_R, h_avg, &
                        & gamma_avg, s_L, s_R, s_S, nbub_L, nbub_R, &
                        & ptilde_L, ptilde_R, vel_avg_rms, Re_L, Re_R, &
                        & pcorr, zcoef, vel_L_tmp, vel_R_tmp]')
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    alpha_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)
                                    alpha_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)
                                end do

                                vel_L_rms = 0._wp; vel_R_rms = 0._wp

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, contxe + i)
                                    vel_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, contxe + i)
                                    vel_L_rms = vel_L_rms + vel_L(i)**2._wp
                                    vel_R_rms = vel_R_rms + vel_R(i)**2._wp
                                end do

                                pres_L = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx)
                                pres_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx)

                                rho_L = 0._wp
                                gamma_L = 0._wp
                                pi_inf_L = 0._wp
                                qv_L = 0._wp

                                ! Retain this in the refactor
                                if (mpp_lim .and. (num_fluids > 2)) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        rho_L = rho_L + qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                        gamma_L = gamma_L + qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)*gammas(i)
                                        pi_inf_L = pi_inf_L + qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)*pi_infs(i)
                                        qv_L = qv_L + qL_prim_rs${XYZ}$_vf(j, k, l, i)*qvs(i)
                                    end do
                                else if (num_fluids > 2) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids - 1
                                        rho_L = rho_L + qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                        gamma_L = gamma_L + qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)*gammas(i)
                                        pi_inf_L = pi_inf_L + qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)*pi_infs(i)
                                        qv_L = qv_L + qL_prim_rs${XYZ}$_vf(j, k, l, i)*qvs(i)
                                    end do
                                else
                                    rho_L = qL_prim_rs${XYZ}$_vf(j, k, l, 1)
                                    gamma_L = gammas(1)
                                    pi_inf_L = pi_infs(1)
                                    qv_L = qvs(1)
                                end if

                                rho_R = 0._wp
                                gamma_R = 0._wp
                                pi_inf_R = 0._wp
                                qv_R = 0._wp

                                if (mpp_lim .and. (num_fluids > 2)) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        rho_R = rho_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)
                                        gamma_R = gamma_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)*gammas(i)
                                        pi_inf_R = pi_inf_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)*pi_infs(i)
                                        qv_R = qv_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)*qvs(i)
                                    end do
                                else if (num_fluids > 2) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids - 1
                                        rho_R = rho_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)
                                        gamma_R = gamma_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)*gammas(i)
                                        pi_inf_R = pi_inf_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)*pi_infs(i)
                                        qv_R = qv_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)*qvs(i)
                                    end do
                                else
                                    rho_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, 1)
                                    gamma_R = gammas(1)
                                    pi_inf_R = pi_infs(1)
                                    qv_R = qvs(1)
                                end if

                                if (viscous) then
                                    if (num_fluids == 1) then ! Need to consider case with num_fluids >= 2
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, 2
                                            Re_L(i) = dflt_real

                                            if (Re_size(i) > 0) Re_L(i) = 0._wp

                                            $:GPU_LOOP(parallelism='[seq]')
                                            do q = 1, Re_size(i)
                                                Re_L(i) = (1._wp - qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + Re_idx(i, q)))/Res(i, q) &
                                                          + Re_L(i)
                                            end do

                                            Re_L(i) = 1._wp/max(Re_L(i), sgm_eps)

                                        end do

                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, 2
                                            Re_R(i) = dflt_real

                                            if (Re_size(i) > 0) Re_R(i) = 0._wp

                                            $:GPU_LOOP(parallelism='[seq]')
                                            do q = 1, Re_size(i)
                                                Re_R(i) = (1._wp - qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + Re_idx(i, q)))/Res(i, q) &
                                                          + Re_R(i)
                                            end do

                                            Re_R(i) = 1._wp/max(Re_R(i), sgm_eps)
                                        end do
                                    end if
                                end if

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
                                            nbub_L = qL_prim_rs${XYZ}$_vf(j, k, l, n_idx)
                                            nbub_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, n_idx)
                                        else
                                            nbub_L_denom = 0._wp
                                            nbub_R_denom = 0._wp
                                            $:GPU_LOOP(parallelism='[seq]')
                                            do i = 1, nb
                                                nbub_L_denom = nbub_L_denom + (R0_L(i)**3._wp)*weight(i)
                                                nbub_R_denom = nbub_R_denom + (R0_R(i)**3._wp)*weight(i)
                                            end do
                                            nbub_L = (3._wp/(4._wp*pi))*qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + num_fluids)/nbub_L_denom
                                            nbub_R = (3._wp/(4._wp*pi))*qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + num_fluids)/nbub_R_denom
                                        end if
                                    else
                                        !nb stored in 0th moment of first R0 bin in variable conversion module
                                        nbub_L = qL_prim_rs${XYZ}$_vf(j, k, l, bubxb)
                                        nbub_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, bubxb)
                                    end if

                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, nb
                                        if (.not. qbmm) then
                                            if (polytropic) then
                                                pbw_L(i) = f_cpbw_KM(R0(i), R0_L(i), V0_L(i), 0._wp)
                                                pbw_R(i) = f_cpbw_KM(R0(i), R0_R(i), V0_R(i), 0._wp)
                                            else
                                                pbw_L(i) = f_cpbw_KM(R0(i), R0_L(i), V0_L(i), P0_L(i))
                                                pbw_R(i) = f_cpbw_KM(R0(i), R0_R(i), V0_R(i), P0_R(i))
                                            end if
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

                                    if (qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + num_fluids) < small_alf .or. R3Lbar < small_alf) then
                                        ptilde_L = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + num_fluids)*pres_L
                                    else
                                        ptilde_L = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + num_fluids)*(pres_L - PbwR3Lbar/R3Lbar - &
                                                                                                      rho_L*R3V2Lbar/R3Lbar)
                                    end if

                                    if (qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + num_fluids) < small_alf .or. R3Rbar < small_alf) then
                                        ptilde_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + num_fluids)*pres_R
                                    else
                                        ptilde_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + num_fluids)*(pres_R - PbwR3Rbar/R3Rbar - &
                                                                                                          rho_R*R3V2Rbar/R3Rbar)
                                    end if

                                    if ((.not. f_approx_equal(ptilde_L, ptilde_L)) .or. (.not. f_approx_equal(ptilde_R, ptilde_R))) then
                                    end if

                                    rho_avg = 5.e-1_wp*(rho_L + rho_R)
                                    H_avg = 5.e-1_wp*(H_L + H_R)
                                    gamma_avg = 5.e-1_wp*(gamma_L + gamma_R)
                                    vel_avg_rms = 0._wp

                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        vel_avg_rms = vel_avg_rms + (5.e-1_wp*(vel_L(i) + vel_R(i)))**2._wp
                                    end do

                                end if

                                call s_compute_speed_of_sound(pres_L, rho_L, gamma_L, pi_inf_L, H_L, alpha_L, &
                                                              vel_L_rms, 0._wp, c_L)

                                call s_compute_speed_of_sound(pres_R, rho_R, gamma_R, pi_inf_R, H_R, alpha_R, &
                                                              vel_R_rms, 0._wp, c_R)

                                !> The computation of c_avg does not require all the variables, and therefore the non '_avg'
                                ! variables are placeholders to call the subroutine.
                                call s_compute_speed_of_sound(pres_R, rho_avg, gamma_avg, pi_inf_R, H_avg, alpha_R, &
                                                              vel_avg_rms, 0._wp, c_avg)

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

                                    s_S = (pres_R - pres_L + rho_L*vel_L(dir_idx(1))* &
                                           (s_L - vel_L(dir_idx(1))) - &
                                           rho_R*vel_R(dir_idx(1))* &
                                           (s_R - vel_R(dir_idx(1)))) &
                                          /(rho_L*(s_L - vel_L(dir_idx(1))) - &
                                            rho_R*(s_R - vel_R(dir_idx(1))))
                                elseif (wave_speeds == 2) then
                                    pres_SL = 5.e-1_wp*(pres_L + pres_R + rho_avg*c_avg* &
                                                        (vel_L(dir_idx(1)) - &
                                                         vel_R(dir_idx(1))))

                                    pres_SR = pres_SL

                                    Ms_L = max(1._wp, sqrt(1._wp + ((5.e-1_wp + gamma_L)/(1._wp + gamma_L))* &
                                                           (pres_SL/pres_L - 1._wp)*pres_L/ &
                                                           ((pres_L + pi_inf_L/(1._wp + gamma_L)))))
                                    Ms_R = max(1._wp, sqrt(1._wp + ((5.e-1_wp + gamma_R)/(1._wp + gamma_R))* &
                                                           (pres_SR/pres_R - 1._wp)*pres_R/ &
                                                           ((pres_R + pi_inf_R/(1._wp + gamma_R)))))

                                    s_L = vel_L(dir_idx(1)) - c_L*Ms_L
                                    s_R = vel_R(dir_idx(1)) + c_R*Ms_R

                                    s_S = 5.e-1_wp*((vel_L(dir_idx(1)) + vel_R(dir_idx(1))) + &
                                                    (pres_L - pres_R)/ &
                                                    (rho_avg*c_avg))
                                end if

                                ! follows Einfeldt et al.
                                ! s_M/P = min/max(0.,s_L/R)
                                s_M = min(0._wp, s_L); s_P = max(0._wp, s_R)

                                ! goes with q_star_L/R = xi_L/R * (variable)
                                ! xi_L/R = ( ( s_L/R - u_L/R )/(s_L/R - s_star) )
                                xi_L = (s_L - vel_L(dir_idx(1)))/(s_L - s_S)
                                xi_R = (s_R - vel_R(dir_idx(1)))/(s_R - s_S)

                                ! goes with numerical velocity in x/y/z directions
                                ! xi_P/M = 0.5 +/m sgn(0.5,s_star)
                                xi_M = (5.e-1_wp + sign(5.e-1_wp, s_S))
                                xi_P = (5.e-1_wp - sign(5.e-1_wp, s_S))

                                ! Low Mach correction
                                if (low_Mach == 1) then
                                    @:compute_low_Mach_correction()
                                else
                                    pcorr = 0._wp
                                end if

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, contxe
                                    flux_rs${XYZ}$_vf(j, k, l, i) = &
                                        xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, i) &
                                        *(vel_L(dir_idx(1)) + s_M*(xi_L - 1._wp)) &
                                        + xi_P*qR_prim_rs${XYZ}$_vf(j + 1, k, l, i) &
                                        *(vel_R(dir_idx(1)) + s_P*(xi_R - 1._wp))
                                end do

                                if (bubbles_euler .and. (num_fluids > 1)) then
                                    ! Kill mass transport @ gas density
                                    flux_rs${XYZ}$_vf(j, k, l, contxe) = 0._wp
                                end if

                                ! Momentum flux.
                                ! f = \rho u u + p I, q = \rho u, q_star = \xi * \rho*(s_star, v, w)

                                ! Include p_tilde

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    flux_rs${XYZ}$_vf(j, k, l, contxe + dir_idx(i)) = &
                                        xi_M*(rho_L*(vel_L(dir_idx(1))* &
                                                     vel_L(dir_idx(i)) + &
                                                     s_M*(xi_L*(dir_flg(dir_idx(i))*s_S + &
                                                                (1._wp - dir_flg(dir_idx(i)))* &
                                                                vel_L(dir_idx(i))) - vel_L(dir_idx(i)))) + &
                                              dir_flg(dir_idx(i))*(pres_L - ptilde_L)) &
                                        + xi_P*(rho_R*(vel_R(dir_idx(1))* &
                                                       vel_R(dir_idx(i)) + &
                                                       s_P*(xi_R*(dir_flg(dir_idx(i))*s_S + &
                                                                  (1._wp - dir_flg(dir_idx(i)))* &
                                                                  vel_R(dir_idx(i))) - vel_R(dir_idx(i)))) + &
                                                dir_flg(dir_idx(i))*(pres_R - ptilde_R)) &
                                        + (s_M/s_L)*(s_P/s_R)*dir_flg(dir_idx(i))*pcorr
                                end do

                                ! Energy flux.
                                ! f = u*(E+p), q = E, q_star = \xi*E+(s-u)(\rho s_star + p/(s-u))
                                flux_rs${XYZ}$_vf(j, k, l, E_idx) = &
                                    xi_M*(vel_L(dir_idx(1))*(E_L + pres_L - ptilde_L) + &
                                          s_M*(xi_L*(E_L + (s_S - vel_L(dir_idx(1)))* &
                                                     (rho_L*s_S + (pres_L - ptilde_L)/ &
                                                      (s_L - vel_L(dir_idx(1))))) - E_L)) &
                                    + xi_P*(vel_R(dir_idx(1))*(E_R + pres_R - ptilde_R) + &
                                            s_P*(xi_R*(E_R + (s_S - vel_R(dir_idx(1)))* &
                                                       (rho_R*s_S + (pres_R - ptilde_R)/ &
                                                        (s_R - vel_R(dir_idx(1))))) - E_R)) &
                                    + (s_M/s_L)*(s_P/s_R)*pcorr*s_S

                                ! Volume fraction flux
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = advxb, advxe
                                    flux_rs${XYZ}$_vf(j, k, l, i) = &
                                        xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, i) &
                                        *(vel_L(dir_idx(1)) + s_M*(xi_L - 1._wp)) &
                                        + xi_P*qR_prim_rs${XYZ}$_vf(j + 1, k, l, i) &
                                        *(vel_R(dir_idx(1)) + s_P*(xi_R - 1._wp))
                                end do

                                ! Source for volume fraction advection equation
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel_src_rs${XYZ}$_vf(j, k, l, dir_idx(i)) = &
                                        xi_M*(vel_L(dir_idx(i)) + &
                                              dir_flg(dir_idx(i))* &
                                              s_M*(xi_L - 1._wp)) &
                                        + xi_P*(vel_R(dir_idx(i)) + &
                                                dir_flg(dir_idx(i))* &
                                                s_P*(xi_R - 1._wp))

                                    !IF ( (model_eqns == 4) .or. (num_fluids==1) ) vel_src_rs_vf(idxi)%sf(j,k,l) = 0._wp
                                end do

                                flux_src_rs${XYZ}$_vf(j, k, l, advxb) = vel_src_rs${XYZ}$_vf(j, k, l, dir_idx(1))

                                ! Add advection flux for bubble variables
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = bubxb, bubxe
                                    flux_rs${XYZ}$_vf(j, k, l, i) = &
                                        xi_M*nbub_L*qL_prim_rs${XYZ}$_vf(j, k, l, i) &
                                        *(vel_L(dir_idx(1)) + s_M*(xi_L - 1._wp)) &
                                        + xi_P*nbub_R*qR_prim_rs${XYZ}$_vf(j + 1, k, l, i) &
                                        *(vel_R(dir_idx(1)) + s_P*(xi_R - 1._wp))
                                end do

                                if (qbmm) then
                                    flux_rs${XYZ}$_vf(j, k, l, bubxb) = &
                                        xi_M*nbub_L &
                                        *(vel_L(dir_idx(1)) + s_M*(xi_L - 1._wp)) &
                                        + xi_P*nbub_R &
                                        *(vel_R(dir_idx(1)) + s_P*(xi_R - 1._wp))
                                end if

                                if (adv_n) then
                                    flux_rs${XYZ}$_vf(j, k, l, n_idx) = &
                                        xi_M*nbub_L &
                                        *(vel_L(dir_idx(1)) + s_M*(xi_L - 1._wp)) &
                                        + xi_P*nbub_R &
                                        *(vel_R(dir_idx(1)) + s_P*(xi_R - 1._wp))
                                end if

                                ! Geometrical source flux for cylindrical coordinates
                                #:if (NORM_DIR == 2)
                                    if (cyl_coord) then
                                        ! Substituting the advective flux into the inviscid geometrical source flux
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, E_idx
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = flux_rs${XYZ}$_vf(j, k, l, i)
                                        end do
                                        ! Recalculating the radial momentum geometric source flux
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, contxe + dir_idx(1)) = &
                                            xi_M*(rho_L*(vel_L(dir_idx(1))* &
                                                         vel_L(dir_idx(1)) + &
                                                         s_M*(xi_L*(dir_flg(dir_idx(1))*s_S + &
                                                                    (1._wp - dir_flg(dir_idx(1)))* &
                                                                    vel_L(dir_idx(1))) - vel_L(dir_idx(1))))) &
                                            + xi_P*(rho_R*(vel_R(dir_idx(1))* &
                                                           vel_R(dir_idx(1)) + &
                                                           s_P*(xi_R*(dir_flg(dir_idx(1))*s_S + &
                                                                      (1._wp - dir_flg(dir_idx(1)))* &
                                                                      vel_R(dir_idx(1))) - vel_R(dir_idx(1)))))
                                        ! Geometrical source of the void fraction(s) is zero
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = advxb, advxe
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

                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, momxb + 1) = &
                                            -xi_M*(rho_L*(vel_L(dir_idx(1))* &
                                                          vel_L(dir_idx(1)) + &
                                                          s_M*(xi_L*(dir_flg(dir_idx(1))*s_S + &
                                                                     (1._wp - dir_flg(dir_idx(1)))* &
                                                                     vel_L(dir_idx(1))) - vel_L(dir_idx(1))))) &
                                            - xi_P*(rho_R*(vel_R(dir_idx(1))* &
                                                           vel_R(dir_idx(1)) + &
                                                           s_P*(xi_R*(dir_flg(dir_idx(1))*s_S + &
                                                                      (1._wp - dir_flg(dir_idx(1)))* &
                                                                      vel_R(dir_idx(1))) - vel_R(dir_idx(1)))))
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, momxe) = flux_rs${XYZ}$_vf(j, k, l, momxb + 1)

                                    end if
                                #:endif
                            end do
                        end do
                    end do
                else
                    ! 5-EQUATION MODEL WITH HLLC
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[vel_L, vel_R, &
                        & Re_L, Re_R, rho_avg, h_avg, gamma_avg, &
                        & alpha_L, alpha_R, s_L, s_R, s_S, &
                        & vel_avg_rms, pcorr, zcoef, vel_L_tmp, &
                        & vel_R_tmp, Ys_L, Ys_R, Xs_L, Xs_R, &
                        & Gamma_iL, Gamma_iR, Cp_iL, Cp_iR, tau_e_L, &
                        & tau_e_R, xi_field_L, xi_field_R, Yi_avg, &
                        & Phi_avg, h_iL, h_iR, h_avg_2]', &
                        & copyin='[is1, is2, is3]')
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end

                                !idx1 = 1; if (dir_idx(1) == 2) idx1 = 2; if (dir_idx(1) == 3) idx1 = 3

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    alpha_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)
                                    alpha_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)
                                end do

                                vel_L_rms = 0._wp; vel_R_rms = 0._wp
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    vel_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, contxe + i)
                                    vel_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, contxe + i)
                                    vel_L_rms = vel_L_rms + vel_L(i)**2._wp
                                    vel_R_rms = vel_R_rms + vel_R(i)**2._wp
                                end do

                                pres_L = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx)
                                pres_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx)

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

                                ! Change this by splitting it into the cases
                                ! present in the bubbles_euler
                                if (mpp_lim) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        qL_prim_rs${XYZ}$_vf(j, k, l, i) = max(0._wp, qL_prim_rs${XYZ}$_vf(j, k, l, i))
                                        qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i) = min(max(0._wp, qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)), 1._wp)
                                        alpha_L_sum = alpha_L_sum + qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)
                                    end do

                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i) = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)/max(alpha_L_sum, sgm_eps)
                                    end do

                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        qR_prim_rs${XYZ}$_vf(j + 1, k, l, i) = max(0._wp, qR_prim_rs${XYZ}$_vf(j + 1, k, l, i))
                                        qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i) = min(max(0._wp, qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)), 1._wp)
                                        alpha_R_sum = alpha_R_sum + qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)
                                    end do

                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)/max(alpha_R_sum, sgm_eps)
                                    end do
                                end if

                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    rho_L = rho_L + qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                    gamma_L = gamma_L + qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)*gammas(i)
                                    pi_inf_L = pi_inf_L + qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)*pi_infs(i)
                                    qv_L = qv_L + qL_prim_rs${XYZ}$_vf(j, k, l, i)*qvs(i)

                                    rho_R = rho_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)
                                    gamma_R = gamma_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)*gammas(i)
                                    pi_inf_R = pi_inf_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)*pi_infs(i)
                                    qv_R = qv_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)*qvs(i)
                                end do

                                if (viscous) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, 2
                                        Re_L(i) = dflt_real

                                        if (Re_size(i) > 0) Re_L(i) = 0._wp

                                        $:GPU_LOOP(parallelism='[seq]')
                                        do q = 1, Re_size(i)
                                            Re_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + Re_idx(i, q))/Res(i, q) &
                                                      + Re_L(i)
                                        end do

                                        Re_L(i) = 1._wp/max(Re_L(i), sgm_eps)

                                    end do

                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, 2
                                        Re_R(i) = dflt_real

                                        if (Re_size(i) > 0) Re_R(i) = 0._wp

                                        $:GPU_LOOP(parallelism='[seq]')
                                        do q = 1, Re_size(i)
                                            Re_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + Re_idx(i, q))/Res(i, q) &
                                                      + Re_R(i)
                                        end do

                                        Re_R(i) = 1._wp/max(Re_R(i), sgm_eps)
                                    end do
                                end if

                                if (chemistry) then
                                    c_sum_Yi_Phi = 0.0_wp
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = chemxb, chemxe
                                        Ys_L(i - chemxb + 1) = qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                        Ys_R(i - chemxb + 1) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)
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
                                else
                                    E_L = gamma_L*pres_L + pi_inf_L + 5.e-1*rho_L*vel_L_rms + qv_L

                                    E_R = gamma_R*pres_R + pi_inf_R + 5.e-1*rho_R*vel_R_rms + qv_R

                                    H_L = (E_L + pres_L)/rho_L
                                    H_R = (E_R + pres_R)/rho_R
                                end if

                                ! ENERGY ADJUSTMENTS FOR HYPOELASTIC ENERGY
                                if (hypoelasticity) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, strxe - strxb + 1
                                        tau_e_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, strxb - 1 + i)
                                        tau_e_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, strxb - 1 + i)
                                    end do
                                    G_L = 0._wp
                                    G_R = 0._wp
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        G_L = G_L + alpha_L(i)*Gs(i)
                                        G_R = G_R + alpha_R(i)*Gs(i)
                                    end do
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, strxe - strxb + 1
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

                                ! ENERGY ADJUSTMENTS FOR HYPERELASTIC ENERGY
                                if (hyperelasticity) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        xi_field_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, xibeg - 1 + i)
                                        xi_field_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, xibeg - 1 + i)
                                    end do
                                    G_L = 0._wp
                                    G_R = 0._wp
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_fluids
                                        ! Mixture left and right shear modulus
                                        G_L = G_L + alpha_L(i)*Gs(i)
                                        G_R = G_R + alpha_R(i)*Gs(i)
                                    end do
                                    ! Elastic contribution to energy if G large enough
                                    if (G_L > verysmall .and. G_R > verysmall) then
                                        E_L = E_L + G_L*qL_prim_rs${XYZ}$_vf(j, k, l, xiend + 1)
                                        E_R = E_R + G_R*qR_prim_rs${XYZ}$_vf(j + 1, k, l, xiend + 1)
                                    end if
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, b_size - 1
                                        tau_e_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, strxb - 1 + i)
                                        tau_e_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, strxb - 1 + i)
                                    end do
                                end if

                                H_L = (E_L + pres_L)/rho_L
                                H_R = (E_R + pres_R)/rho_R

                                @:compute_average_state()

                                call s_compute_speed_of_sound(pres_L, rho_L, gamma_L, pi_inf_L, H_L, alpha_L, &
                                                              vel_L_rms, 0._wp, c_L)

                                call s_compute_speed_of_sound(pres_R, rho_R, gamma_R, pi_inf_R, H_R, alpha_R, &
                                                              vel_R_rms, 0._wp, c_R)

                                !> The computation of c_avg does not require all the variables, and therefore the non '_avg'
                                ! variables are placeholders to call the subroutine.
                                call s_compute_speed_of_sound(pres_R, rho_avg, gamma_avg, pi_inf_R, H_avg, alpha_R, &
                                                              vel_avg_rms, c_sum_Yi_Phi, c_avg)

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
                                    if (elasticity) then
                                        s_L = min(vel_L(dir_idx(1)) - sqrt(c_L*c_L + &
                                                                           (((4._wp*G_L)/3._wp) + tau_e_L(dir_idx_tau(1)))/rho_L), vel_R(dir_idx(1)) - sqrt(c_R*c_R + &
                                                                                                                                                            (((4._wp*G_R)/3._wp) + tau_e_R(dir_idx_tau(1)))/rho_R))
                                        s_R = max(vel_R(dir_idx(1)) + sqrt(c_R*c_R + &
                                                                           (((4._wp*G_R)/3._wp) + tau_e_R(dir_idx_tau(1)))/rho_R), vel_L(dir_idx(1)) + sqrt(c_L*c_L + &
                                                                                                                                                            (((4._wp*G_L)/3._wp) + tau_e_L(dir_idx_tau(1)))/rho_L))
                                        s_S = (pres_R - tau_e_R(dir_idx_tau(1)) - pres_L + &
                                               tau_e_L(dir_idx_tau(1)) + rho_L*vel_L(idx1)*(s_L - vel_L(idx1)) - &
                                               rho_R*vel_R(idx1)*(s_R - vel_R(idx1)))/(rho_L*(s_L - vel_L(idx1)) - &
                                                                                       rho_R*(s_R - vel_R(idx1)))
                                    else
                                        s_L = min(vel_L(dir_idx(1)) - c_L, vel_R(dir_idx(1)) - c_R)
                                        s_R = max(vel_R(dir_idx(1)) + c_R, vel_L(dir_idx(1)) + c_L)
                                        s_S = (pres_R - pres_L + rho_L*vel_L(dir_idx(1))* &
                                               (s_L - vel_L(dir_idx(1))) - rho_R*vel_R(dir_idx(1))*(s_R - vel_R(dir_idx(1)))) &
                                              /(rho_L*(s_L - vel_L(dir_idx(1))) - rho_R*(s_R - vel_R(dir_idx(1))))

                                    end if
                                elseif (wave_speeds == 2) then
                                    pres_SL = 5.e-1_wp*(pres_L + pres_R + rho_avg*c_avg* &
                                                        (vel_L(idx1) - &
                                                         vel_R(idx1)))

                                    pres_SR = pres_SL

                                    Ms_L = max(1._wp, sqrt(1._wp + ((5.e-1_wp + gamma_L)/(1._wp + gamma_L))* &
                                                           (pres_SL/pres_L - 1._wp)*pres_L/ &
                                                           ((pres_L + pi_inf_L/(1._wp + gamma_L)))))
                                    Ms_R = max(1._wp, sqrt(1._wp + ((5.e-1_wp + gamma_R)/(1._wp + gamma_R))* &
                                                           (pres_SR/pres_R - 1._wp)*pres_R/ &
                                                           ((pres_R + pi_inf_R/(1._wp + gamma_R)))))

                                    s_L = vel_L(idx1) - c_L*Ms_L
                                    s_R = vel_R(idx1) + c_R*Ms_R

                                    s_S = 5.e-1_wp*((vel_L(idx1) + vel_R(idx1)) + &
                                                    (pres_L - pres_R)/ &
                                                    (rho_avg*c_avg))
                                end if

                                ! follows Einfeldt et al.
                                ! s_M/P = min/max(0.,s_L/R)
                                s_M = min(0._wp, s_L); s_P = max(0._wp, s_R)

                                ! goes with q_star_L/R = xi_L/R * (variable)
                                ! xi_L/R = ( ( s_L/R - u_L/R )/(s_L/R - s_star) )
                                xi_L = (s_L - vel_L(idx1))/(s_L - s_S)
                                xi_R = (s_R - vel_R(idx1))/(s_R - s_S)

                                ! goes with numerical velocity in x/y/z directions
                                ! xi_P/M = 0.5 +/m sgn(0.5,s_star)
                                xi_M = (5.e-1_wp + sign(5.e-1_wp, s_S))
                                xi_P = (5.e-1_wp - sign(5.e-1_wp, s_S))

                                ! Low Mach correction
                                if (low_Mach == 1) then
                                    @:compute_low_Mach_correction()
                                else
                                    pcorr = 0._wp
                                end if

                                ! COMPUTING THE HLLC FLUXES
                                ! MASS FLUX.
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, contxe
                                    flux_rs${XYZ}$_vf(j, k, l, i) = &
                                        xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, i) &
                                        *(vel_L(idx1) + s_M*(xi_L - 1._wp)) &
                                        + xi_P*qR_prim_rs${XYZ}$_vf(j + 1, k, l, i) &
                                        *(vel_R(idx1) + s_P*(xi_R - 1._wp))
                                end do

                                ! MOMENTUM FLUX.
                                ! f = \rho u u - \sigma, q = \rho u, q_star = \xi * \rho*(s_star, v, w)
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    idxi = dir_idx(i)
                                    flux_rs${XYZ}$_vf(j, k, l, contxe + idxi) = &
                                        xi_M*(rho_L*(vel_L(idx1)* &
                                                     vel_L(idxi) + &
                                                     s_M*(xi_L*(dir_flg(idxi)*s_S + &
                                                                (1._wp - dir_flg(idxi))* &
                                                                vel_L(idxi)) - vel_L(idxi))) + &
                                              dir_flg(idxi)*(pres_L)) &
                                        + xi_P*(rho_R*(vel_R(idx1)* &
                                                       vel_R(idxi) + &
                                                       s_P*(xi_R*(dir_flg(idxi)*s_S + &
                                                                  (1._wp - dir_flg(idxi))* &
                                                                  vel_R(idxi)) - vel_R(idxi))) + &
                                                dir_flg(idxi)*(pres_R)) &
                                        + (s_M/s_L)*(s_P/s_R)*dir_flg(idxi)*pcorr
                                end do

                                ! ENERGY FLUX.
                                ! f = u*(E-\sigma), q = E, q_star = \xi*E+(s-u)(\rho s_star - \sigma/(s-u))
                                flux_rs${XYZ}$_vf(j, k, l, E_idx) = &
                                    xi_M*(vel_L(idx1)*(E_L + pres_L) + &
                                          s_M*(xi_L*(E_L + (s_S - vel_L(idx1))* &
                                                     (rho_L*s_S + pres_L/ &
                                                      (s_L - vel_L(idx1)))) - E_L)) &
                                    + xi_P*(vel_R(idx1)*(E_R + pres_R) + &
                                            s_P*(xi_R*(E_R + (s_S - vel_R(idx1))* &
                                                       (rho_R*s_S + pres_R/ &
                                                        (s_R - vel_R(idx1)))) - E_R)) &
                                    + (s_M/s_L)*(s_P/s_R)*pcorr*s_S

                                ! ELASTICITY. Elastic shear stress additions for the momentum and energy flux
                                if (elasticity) then
                                    flux_ene_e = 0._wp
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        idxi = dir_idx(i)
                                        ! MOMENTUM ELASTIC FLUX.
                                        flux_rs${XYZ}$_vf(j, k, l, contxe + idxi) = &
                                            flux_rs${XYZ}$_vf(j, k, l, contxe + idxi) &
                                            - xi_M*tau_e_L(dir_idx_tau(i)) - xi_P*tau_e_R(dir_idx_tau(i))
                                        ! ENERGY ELASTIC FLUX.
                                        flux_ene_e = flux_ene_e - &
                                                     xi_M*(vel_L(idxi)*tau_e_L(dir_idx_tau(i)) + &
                                                           s_M*(xi_L*((s_S - vel_L(i))*(tau_e_L(dir_idx_tau(i))/(s_L - vel_L(i)))))) - &
                                                     xi_P*(vel_R(idxi)*tau_e_R(dir_idx_tau(i)) + &
                                                           s_P*(xi_R*((s_S - vel_R(i))*(tau_e_R(dir_idx_tau(i))/(s_R - vel_R(i))))))
                                    end do
                                    flux_rs${XYZ}$_vf(j, k, l, E_idx) = flux_rs${XYZ}$_vf(j, k, l, E_idx) + flux_ene_e
                                end if

                                ! HYPOELASTIC STRESS EVOLUTION FLUX.
                                if (hypoelasticity) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, strxe - strxb + 1
                                        flux_rs${XYZ}$_vf(j, k, l, strxb - 1 + i) = &
                                            xi_M*(s_S/(s_L - s_S))*(s_L*rho_L*tau_e_L(i) - rho_L*vel_L(idx1)*tau_e_L(i)) + &
                                            xi_P*(s_S/(s_R - s_S))*(s_R*rho_R*tau_e_R(i) - rho_R*vel_R(idx1)*tau_e_R(i))
                                    end do
                                end if

                                ! VOLUME FRACTION FLUX.
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = advxb, advxe
                                    flux_rs${XYZ}$_vf(j, k, l, i) = &
                                        xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, i) &
                                        *(vel_L(idx1) + s_M*(xi_L - 1._wp)) &
                                        + xi_P*qR_prim_rs${XYZ}$_vf(j + 1, k, l, i) &
                                        *(vel_R(idx1) + s_P*(xi_R - 1._wp))
                                end do

                                ! VOLUME FRACTION SOURCE FLUX.
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_dims
                                    idxi = dir_idx(i)
                                    vel_src_rs${XYZ}$_vf(j, k, l, idxi) = &
                                        xi_M*(vel_L(idxi) + &
                                              dir_flg(idxi)* &
                                              s_M*(xi_L - 1._wp)) &
                                        + xi_P*(vel_R(idxi) + &
                                                dir_flg(idxi)* &
                                                s_P*(xi_R - 1._wp))
                                end do

                                ! COLOR FUNCTION FLUX
                                if (surface_tension) then
                                    flux_rs${XYZ}$_vf(j, k, l, c_idx) = &
                                        xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, c_idx) &
                                        *(vel_L(idx1) + s_M*(xi_L - 1._wp)) &
                                        + xi_P*qR_prim_rs${XYZ}$_vf(j + 1, k, l, c_idx) &
                                        *(vel_R(idx1) + s_P*(xi_R - 1._wp))
                                end if

                                ! REFERENCE MAP FLUX.
                                if (hyperelasticity) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, num_dims
                                        flux_rs${XYZ}$_vf(j, k, l, xibeg - 1 + i) = &
                                            xi_M*(s_S/(s_L - s_S))*(s_L*rho_L*xi_field_L(i) &
                                                                    - rho_L*vel_L(idx1)*xi_field_L(i)) + &
                                            xi_P*(s_S/(s_R - s_S))*(s_R*rho_R*xi_field_R(i) &
                                                                    - rho_R*vel_R(idx1)*xi_field_R(i))
                                    end do
                                end if

                                flux_src_rs${XYZ}$_vf(j, k, l, advxb) = vel_src_rs${XYZ}$_vf(j, k, l, idx1)

                                if (chemistry) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = chemxb, chemxe
                                        Y_L = qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                        Y_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)

                                        flux_rs${XYZ}$_vf(j, k, l, i) = xi_M*rho_L*Y_L*(vel_L(idx1) + s_M*(xi_L - 1._wp)) &
                                                                        + xi_P*rho_R*Y_R*(vel_R(idx1) + s_P*(xi_R - 1._wp))
                                        flux_src_rs${XYZ}$_vf(j, k, l, i) = 0.0_wp
                                    end do
                                end if

                                ! Geometrical source flux for cylindrical coordinates
                                #:if (NORM_DIR == 2)
                                    if (cyl_coord) then
                                        !Substituting the advective flux into the inviscid geometrical source flux
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = 1, E_idx
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = flux_rs${XYZ}$_vf(j, k, l, i)
                                        end do
                                        ! Recalculating the radial momentum geometric source flux
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, contxe + idx1) = &
                                            xi_M*(rho_L*(vel_L(idx1)* &
                                                         vel_L(idx1) + &
                                                         s_M*(xi_L*(dir_flg(idx1)*s_S + &
                                                                    (1._wp - dir_flg(idx1))* &
                                                                    vel_L(idx1)) - vel_L(idx1)))) &
                                            + xi_P*(rho_R*(vel_R(idx1)* &
                                                           vel_R(idx1) + &
                                                           s_P*(xi_R*(dir_flg(idx1)*s_S + &
                                                                      (1._wp - dir_flg(idx1))* &
                                                                      vel_R(idx1)) - vel_R(idx1))))
                                        ! Geometrical source of the void fraction(s) is zero
                                        $:GPU_LOOP(parallelism='[seq]')
                                        do i = advxb, advxe
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

                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, momxb + 1) = &
                                            -xi_M*(rho_L*(vel_L(idx1)* &
                                                          vel_L(idx1) + &
                                                          s_M*(xi_L*(dir_flg(idx1)*s_S + &
                                                                     (1._wp - dir_flg(idx1))* &
                                                                     vel_L(idx1)) - vel_L(idx1)))) &
                                            - xi_P*(rho_R*(vel_R(idx1)* &
                                                           vel_R(idx1) + &
                                                           s_P*(xi_R*(dir_flg(idx1)*s_S + &
                                                                      (1._wp - dir_flg(idx1))* &
                                                                      vel_R(idx1)) - vel_R(idx1))))
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, momxe) = flux_rs${XYZ}$_vf(j, k, l, momxb + 1)

                                    end if
                                #:endif

                            end do
                        end do
                    end do
                end if
            end if
        #:endfor
        ! Computing HLLC flux and source flux for Euler system of equations

        if (viscous) then
            if (weno_Re_flux) then
                call s_compute_viscous_source_flux( &
                    qL_prim_vf(momxb:momxe), &
                    dqL_prim_dx_vf(momxb:momxe), &
                    dqL_prim_dy_vf(momxb:momxe), &
                    dqL_prim_dz_vf(momxb:momxe), &
                    qR_prim_vf(momxb:momxe), &
                    dqR_prim_dx_vf(momxb:momxe), &
                    dqR_prim_dy_vf(momxb:momxe), &
                    dqR_prim_dz_vf(momxb:momxe), &
                    flux_src_vf, norm_dir, ix, iy, iz)
            else
                call s_compute_viscous_source_flux( &
                    q_prim_vf(momxb:momxe), &
                    dqL_prim_dx_vf(momxb:momxe), &
                    dqL_prim_dy_vf(momxb:momxe), &
                    dqL_prim_dz_vf(momxb:momxe), &
                    q_prim_vf(momxb:momxe), &
                    dqR_prim_dx_vf(momxb:momxe), &
                    dqR_prim_dy_vf(momxb:momxe), &
                    dqR_prim_dz_vf(momxb:momxe), &
                    flux_src_vf, norm_dir, ix, iy, iz)
            end if
        end if

        if (surface_tension) then
            call s_compute_capillary_source_flux( &
                vel_src_rsx_vf, &
                vel_src_rsy_vf, &
                vel_src_rsz_vf, &
                flux_src_vf, &
                norm_dir, isx, isy, isz)
        end if

        call s_finalize_riemann_solver(flux_vf, flux_src_vf, &
                                       flux_gsrc_vf, &
                                       norm_dir)

    end subroutine s_hllc_riemann_solver

    !> HLLD Riemann solver resolves 5 of the 7 waves of MHD equations:
        !!      1 entropy wave, 2 Alfvén waves, 2 fast magnetosonic waves.
    subroutine s_hlld_riemann_solver(qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, &
                                     dqL_prim_dx_vf, dqL_prim_dy_vf, dqL_prim_dz_vf, &
                                     qL_prim_vf, &
                                     qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, &
                                     dqR_prim_dx_vf, dqR_prim_dy_vf, dqR_prim_dz_vf, &
                                     qR_prim_vf, &
                                     q_prim_vf, &
                                     flux_vf, flux_src_vf, flux_gsrc_vf, &
                                     norm_dir, ix, iy, iz)

        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:), intent(inout) :: qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, &
                                                                                                     qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf

        type(scalar_field), allocatable, dimension(:), intent(inout) :: dqL_prim_dx_vf, dqR_prim_dx_vf, &
                                                                        dqL_prim_dy_vf, dqR_prim_dy_vf, &
                                                                        dqL_prim_dz_vf, dqR_prim_dz_vf

        type(scalar_field), allocatable, dimension(:), intent(inout) :: qL_prim_vf, qR_prim_vf

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: flux_vf, flux_src_vf, flux_gsrc_vf

        integer, intent(in) :: norm_dir
        type(int_bounds_info), intent(in) :: ix, iy, iz

        ! Local variables:
        real(wp), dimension(num_fluids) :: alpha_L, alpha_R, alpha_rho_L, alpha_rho_R
        type(riemann_states_vec3) :: vel
        type(riemann_states) :: rho, pres, E, H_no_mag
        type(riemann_states) :: gamma, pi_inf, qv
        type(riemann_states) :: vel_rms

        type(riemann_states_vec3) :: B
        type(riemann_states) :: c, c_fast, pres_mag

        ! HLLD speeds and intermediate state variables:
        real(wp) :: s_L, s_R, s_M, s_starL, s_starR
        real(wp) :: pTot_L, pTot_R, p_star, rhoL_star, rhoR_star, E_starL, E_starR

        real(wp), dimension(7) :: U_L, U_R, U_starL, U_starR, U_doubleL, U_doubleR
        real(wp), dimension(7) :: F_L, F_R, F_starL, F_starR, F_hlld

        ! Indices for U and F: (rho, rho*vel(1), rho*vel(2), rho*vel(3), By, Bz, E)
        !   Note: vel and B are permutated, so vel(1) is the normal velocity, and x is the normal direction
        !   Note: Bx is omitted as the magnetic flux is always zero in the normal direction

        real(wp) :: sqrt_rhoL_star, sqrt_rhoR_star, denom_ds, sign_Bx
        real(wp) :: vL_star, vR_star, wL_star, wR_star
        real(wp) :: v_double, w_double, By_double, Bz_double, E_doubleL, E_doubleR, E_double

        integer :: i, j, k, l

        call s_populate_riemann_states_variables_buffers( &
            qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, dqL_prim_dx_vf, &
            dqL_prim_dy_vf, dqL_prim_dz_vf, &
            qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, dqR_prim_dx_vf, &
            dqR_prim_dy_vf, dqR_prim_dz_vf, &
            norm_dir, ix, iy, iz)

        call s_initialize_riemann_solver( &
            flux_src_vf, norm_dir)

        #:for NORM_DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
            if (norm_dir == ${NORM_DIR}$) then
                $:GPU_PARALLEL_LOOP(collapse=3, private='[alpha_rho_L, &
                    & alpha_rho_R, vel, alpha_L, alpha_R, rho, pres, &
                    & E, H_no_mag, gamma, pi_inf, qv, vel_rms, B, &
                    & c, c_fast, pres_mag, U_L, U_R, U_starL, &
                    & U_starR, U_doubleL, U_doubleR, F_L, F_R, &
                    & F_starL, F_starR, F_hlld]')
                do l = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = is1%beg, is1%end

                            ! (1) Extract the left/right primitive states
                            do i = 1, contxe
                                alpha_rho_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                alpha_rho_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)
                            end do

                            ! NOTE: unlike HLL & HLLC, vel_L here is permutated by dir_idx for simpler logic
                            do i = 1, num_vels
                                vel%L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, contxe + dir_idx(i))
                                vel%R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, contxe + dir_idx(i))
                            end do

                            vel_rms%L = sum(vel%L**2._wp)
                            vel_rms%R = sum(vel%R**2._wp)

                            do i = 1, num_fluids
                                alpha_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)
                                alpha_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)
                            end do

                            pres%L = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx)
                            pres%R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx)

                            ! NOTE: unlike HLL, Bx, By, Bz are permutated by dir_idx for simpler logic
                            if (mhd) then
                                if (n == 0) then ! 1D: constant Bx; By, Bz as variables; only in x so not permutated
                                    B%L = [Bx0, qL_prim_rs${XYZ}$_vf(j, k, l, B_idx%beg), qL_prim_rs${XYZ}$_vf(j, k, l, B_idx%beg + 1)]
                                    B%R = [Bx0, qR_prim_rs${XYZ}$_vf(j + 1, k, l, B_idx%beg), qR_prim_rs${XYZ}$_vf(j + 1, k, l, B_idx%beg + 1)]
                                else ! 2D/3D: Bx, By, Bz as variables
                                    B%L = [qL_prim_rs${XYZ}$_vf(j, k, l, B_idx%beg + dir_idx(1) - 1), &
                                           qL_prim_rs${XYZ}$_vf(j, k, l, B_idx%beg + dir_idx(2) - 1), &
                                           qL_prim_rs${XYZ}$_vf(j, k, l, B_idx%beg + dir_idx(3) - 1)]
                                    B%R = [qR_prim_rs${XYZ}$_vf(j + 1, k, l, B_idx%beg + dir_idx(1) - 1), &
                                           qR_prim_rs${XYZ}$_vf(j + 1, k, l, B_idx%beg + dir_idx(2) - 1), &
                                           qR_prim_rs${XYZ}$_vf(j + 1, k, l, B_idx%beg + dir_idx(3) - 1)]
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
                            E%R = gamma%R*pres%R + pi_inf%R + 0.5_wp*rho%R*vel_rms%R + qv%R + pres_mag%R ! includes magnetic energy
                            H_no_mag%L = (E%L + pres%L - pres_mag%L)/rho%L
                            H_no_mag%R = (E%R + pres%R - pres_mag%R)/rho%R ! stagnation enthalpy here excludes magnetic energy (only used to find speed of sound)

                            ! (2) Compute fast wave speeds
                            call s_compute_speed_of_sound(pres%L, rho%L, gamma%L, pi_inf%L, H_no_mag%L, alpha_L, vel_rms%L, 0._wp, c%L)
                            call s_compute_speed_of_sound(pres%R, rho%R, gamma%R, pi_inf%R, H_no_mag%R, alpha_R, vel_rms%R, 0._wp, c%R)
                            call s_compute_fast_magnetosonic_speed(rho%L, c%L, B%L, norm_dir, c_fast%L, H_no_mag%L)
                            call s_compute_fast_magnetosonic_speed(rho%R, c%R, B%R, norm_dir, c_fast%R, H_no_mag%R)

                            ! (3) Compute contact speed s_M [Miyoshi Equ. (38)]
                            s_L = min(vel%L(1) - c_fast%L, vel%R(1) - c_fast%R)
                            s_R = max(vel%R(1) + c_fast%R, vel%L(1) + c_fast%L)

                            pTot_L = pres%L + pres_mag%L
                            pTot_R = pres%R + pres_mag%R

                            s_M = (((s_R - vel%R(1))*rho%R*vel%R(1) - &
                                    (s_L - vel%L(1))*rho%L*vel%L(1) - pTot_R + pTot_L)/ &
                                   ((s_R - vel%R(1))*rho%R - (s_L - vel%L(1))*rho%L))

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
                            ! Compute the star flux using HLL relation
                            F_starL = F_L + s_L*(U_starL - U_L)
                            F_starR = F_R + s_R*(U_starR - U_R)
                            ! Compute the rotational (Alfvén) speeds
                            s_starL = s_M - abs(B%L(1))/sqrt(rhoL_star)
                            s_starR = s_M + abs(B%L(1))/sqrt(rhoR_star)
                            ! Compute the double–star states [Miyoshi Eqns. (59)-(62)]
                            sqrt_rhoL_star = sqrt(rhoL_star); sqrt_rhoR_star = sqrt(rhoR_star)
                            vL_star = vel%L(2); wL_star = vel%L(3)
                            vR_star = vel%R(2); wR_star = vel%R(3)

                            ! (6) Compute the double–star states [Miyoshi Eqns. (59)-(62)]
                            denom_ds = sqrt_rhoL_star + sqrt_rhoR_star
                            sign_Bx = sign(1._wp, B%L(1))
                            v_double = (sqrt_rhoL_star*vL_star + sqrt_rhoR_star*vR_star + (B%R(2) - B%L(2))*sign_Bx)/denom_ds
                            w_double = (sqrt_rhoL_star*wL_star + sqrt_rhoR_star*wR_star + (B%R(3) - B%L(3))*sign_Bx)/denom_ds
                            By_double = (sqrt_rhoL_star*B%R(2) + sqrt_rhoR_star*B%L(2) + sqrt_rhoL_star*sqrt_rhoR_star*(vR_star - vL_star)*sign_Bx)/denom_ds
                            Bz_double = (sqrt_rhoL_star*B%R(3) + sqrt_rhoR_star*B%L(3) + sqrt_rhoL_star*sqrt_rhoR_star*(wR_star - wL_star)*sign_Bx)/denom_ds

                            E_doubleL = E_starL - sqrt_rhoL_star*((vL_star*B%L(2) + wL_star*B%L(3)) - (v_double*By_double + w_double*Bz_double))*sign_Bx
                            E_doubleR = E_starR + sqrt_rhoR_star*((vR_star*B%R(2) + wR_star*B%R(3)) - (v_double*By_double + w_double*Bz_double))*sign_Bx
                            E_double = 0.5_wp*(E_doubleL + E_doubleR)

                            U_doubleL = [rhoL_star, rhoL_star*s_M, rhoL_star*v_double, rhoL_star*w_double, By_double, Bz_double, E_double]
                            U_doubleR = [rhoR_star, rhoR_star*s_M, rhoR_star*v_double, rhoR_star*w_double, By_double, Bz_double, E_double]

                            ! (11) Choose HLLD flux based on wave-speed regions
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

                            ! (12) Reorder and write temporary variables to the flux array
                            ! Mass
                            flux_rs${XYZ}$_vf(j, k, l, 1) = F_hlld(1) ! TODO multi-component
                            ! Momentum
                            flux_rs${XYZ}$_vf(j, k, l, [contxe + dir_idx(1), contxe + dir_idx(2), contxe + dir_idx(3)]) = F_hlld([2, 3, 4])
                            ! Magnetic field
                            if (n == 0) then
                                flux_rs${XYZ}$_vf(j, k, l, [B_idx%beg, B_idx%beg + 1]) = F_hlld([5, 6])
                            else
                                flux_rs${XYZ}$_vf(j, k, l, [B_idx%beg + dir_idx(2) - 1, B_idx%beg + dir_idx(3) - 1]) = F_hlld([5, 6])
                            end if
                            ! Energy
                            flux_rs${XYZ}$_vf(j, k, l, E_idx) = F_hlld(7)
                            ! Partial fraction
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = advxb, advxe
                                flux_rs${XYZ}$_vf(j, k, l, i) = 0._wp ! TODO multi-component (zero for now)
                            end do

                            flux_src_rs${XYZ}$_vf(j, k, l, advxb) = 0._wp
                        end do
                    end do
                end do
            end if
        #:endfor

        call s_finalize_riemann_solver(flux_vf, flux_src_vf, flux_gsrc_vf, &
                                       norm_dir)
    end subroutine s_hlld_riemann_solver

    !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    impure subroutine s_initialize_riemann_solvers_module

        ! Allocating the variables that will be utilized to formulate the
        ! left, right, and average states of the Riemann problem, as well
        ! the Riemann problem solution
        integer :: i, j

        @:ALLOCATE(Gs(1:num_fluids))

        do i = 1, num_fluids
            Gs(i) = fluid_pp(i)%G
        end do
        $:GPU_UPDATE(device='[Gs]')

        if (viscous) then
            @:ALLOCATE(Res(1:2, 1:Re_size_max))
        end if

        if (viscous) then
            do i = 1, 2
                do j = 1, Re_size(i)
                    Res(i, j) = fluid_pp(Re_idx(i, j))%Re(i)
                end do
            end do
            $:GPU_UPDATE(device='[Res,Re_idx,Re_size]')
        end if

        $:GPU_ENTER_DATA(copyin='[is1,is2,is3,isx,isy,isz]')

        is1%beg = -1; is2%beg = 0; is3%beg = 0
        is1%end = m; is2%end = n; is3%end = p

        @:ALLOCATE(flux_rsx_vf(is1%beg:is1%end, &
            is2%beg:is2%end, &
            is3%beg:is3%end, 1:sys_size))
        @:ALLOCATE(flux_gsrc_rsx_vf(is1%beg:is1%end, &
            is2%beg:is2%end, &
            is3%beg:is3%end, 1:sys_size))
        @:ALLOCATE(flux_src_rsx_vf(is1%beg:is1%end, &
            is2%beg:is2%end, &
            is3%beg:is3%end, advxb:sys_size))
        @:ALLOCATE(vel_src_rsx_vf(is1%beg:is1%end, &
            is2%beg:is2%end, &
            is3%beg:is3%end, 1:num_vels))
        if (qbmm) then
            @:ALLOCATE(mom_sp_rsx_vf(is1%beg:is1%end + 1, is2%beg:is2%end, is3%beg:is3%end, 1:4))
        end if

        if (viscous) then
            @:ALLOCATE(Re_avg_rsx_vf(is1%beg:is1%end, &
                is2%beg:is2%end, &
                is3%beg:is3%end, 1:2))
        end if

        if (n == 0) return

        is1%beg = -1; is2%beg = 0; is3%beg = 0
        is1%end = n; is2%end = m; is3%end = p

        @:ALLOCATE(flux_rsy_vf(is1%beg:is1%end, &
            is2%beg:is2%end, &
            is3%beg:is3%end, 1:sys_size))
        @:ALLOCATE(flux_gsrc_rsy_vf(is1%beg:is1%end, &
            is2%beg:is2%end, &
            is3%beg:is3%end, 1:sys_size))
        @:ALLOCATE(flux_src_rsy_vf(is1%beg:is1%end, &
            is2%beg:is2%end, &
            is3%beg:is3%end, advxb:sys_size))
        @:ALLOCATE(vel_src_rsy_vf(is1%beg:is1%end, &
            is2%beg:is2%end, &
            is3%beg:is3%end, 1:num_vels))

        if (qbmm) then
            @:ALLOCATE(mom_sp_rsy_vf(is1%beg:is1%end + 1, is2%beg:is2%end, is3%beg:is3%end, 1:4))
        end if

        if (viscous) then
            @:ALLOCATE(Re_avg_rsy_vf(is1%beg:is1%end, &
                is2%beg:is2%end, &
                is3%beg:is3%end, 1:2))
        end if

        if (p == 0) return

        is1%beg = -1; is2%beg = 0; is3%beg = 0
        is1%end = p; is2%end = n; is3%end = m

        @:ALLOCATE(flux_rsz_vf(is1%beg:is1%end, &
            is2%beg:is2%end, &
            is3%beg:is3%end, 1:sys_size))
        @:ALLOCATE(flux_gsrc_rsz_vf(is1%beg:is1%end, &
            is2%beg:is2%end, &
            is3%beg:is3%end, 1:sys_size))
        @:ALLOCATE(flux_src_rsz_vf(is1%beg:is1%end, &
            is2%beg:is2%end, &
            is3%beg:is3%end, advxb:sys_size))
        @:ALLOCATE(vel_src_rsz_vf(is1%beg:is1%end, &
            is2%beg:is2%end, &
            is3%beg:is3%end, 1:num_vels))

        if (qbmm) then
            @:ALLOCATE(mom_sp_rsz_vf(is1%beg:is1%end + 1, is2%beg:is2%end, is3%beg:is3%end, 1:4))
        end if

        if (viscous) then
            @:ALLOCATE(Re_avg_rsz_vf(is1%beg:is1%end, &
                is2%beg:is2%end, &
                is3%beg:is3%end, 1:2))
        end if

    end subroutine s_initialize_riemann_solvers_module

    !>  The purpose of this subroutine is to populate the buffers
        !!      of the left and right Riemann states variables, depending
        !!      on the boundary conditions.
        !!  @param qL_prim_vf The  left WENO-reconstructed cell-boundary values of the
        !!      cell-average primitive variables
        !!  @param qR_prim_vf The right WENO-reconstructed cell-boundary values of the
        !!      cell-average primitive variables
        !!  @param dqL_prim_dx_vf The  left WENO-reconstructed cell-boundary values of the
        !!      first-order x-dir spatial derivatives
        !!  @param dqL_prim_dy_vf The  left WENO-reconstructed cell-boundary values of the
        !!      first-order y-dir spatial derivatives
        !!  @param dqL_prim_dz_vf The  left WENO-reconstructed cell-boundary values of the
        !!      first-order z-dir spatial derivatives
        !!  @param dqR_prim_dx_vf The right WENO-reconstructed cell-boundary values of the
        !!      first-order x-dir spatial derivatives
        !!  @param dqR_prim_dy_vf The right WENO-reconstructed cell-boundary values of the
        !!      first-order y-dir spatial derivatives
        !!  @param dqR_prim_dz_vf The right WENO-reconstructed cell-boundary values of the
        !!      first-order z-dir spatial derivatives
        !!  @param gm_alphaL_vf  Left averaged gradient magnitude
        !!  @param gm_alphaR_vf Right averaged gradient magnitude
        !!  @param norm_dir Dir. splitting direction
        !!  @param ix Index bounds in the x-dir
        !!  @param iy Index bounds in the y-dir
        !!  @param iz Index bounds in the z-dir
    subroutine s_populate_riemann_states_variables_buffers( &
        qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, dqL_prim_dx_vf, &
        dqL_prim_dy_vf, &
        dqL_prim_dz_vf, &
        qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, dqR_prim_dx_vf, &
        dqR_prim_dy_vf, &
        dqR_prim_dz_vf, &
        norm_dir, ix, iy, iz)

        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:), intent(inout) :: qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf

        type(scalar_field), &
            allocatable, dimension(:), &
            intent(inout) :: dqL_prim_dx_vf, dqR_prim_dx_vf, &
                             dqL_prim_dy_vf, dqR_prim_dy_vf, &
                             dqL_prim_dz_vf, dqR_prim_dz_vf

        integer, intent(in) :: norm_dir
        type(int_bounds_info), intent(in) :: ix, iy, iz

        integer :: i, j, k, l !< Generic loop iterator

        if (norm_dir == 1) then
            is1 = ix; is2 = iy; is3 = iz
            dir_idx = (/1, 2, 3/); dir_flg = (/1._wp, 0._wp, 0._wp/)
        elseif (norm_dir == 2) then
            is1 = iy; is2 = ix; is3 = iz
            dir_idx = (/2, 1, 3/); dir_flg = (/0._wp, 1._wp, 0._wp/)
        else
            is1 = iz; is2 = iy; is3 = ix
            dir_idx = (/3, 1, 2/); dir_flg = (/0._wp, 0._wp, 1._wp/)
        end if

        $:GPU_UPDATE(device='[is1,is2,is3]')

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
        $:GPU_UPDATE(device='[isx,isy,isz]')
        ! for stuff in different modules
        $:GPU_UPDATE(device='[dir_idx,dir_flg,dir_idx_tau]')

        ! Population of Buffers in x-direction
        if (norm_dir == 1) then

            if (bc_x%beg == BC_RIEMANN_EXTRAP) then    ! Riemann state extrap. BC at beginning
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            qL_prim_rsx_vf(-1, k, l, i) = &
                                qR_prim_rsx_vf(0, k, l, i)
                        end do
                    end do
                end do

                if (viscous) then
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = momxb, momxe
                        do l = isz%beg, isz%end
                            do k = isy%beg, isy%end

                                dqL_prim_dx_vf(i)%sf(-1, k, l) = &
                                    dqR_prim_dx_vf(i)%sf(0, k, l)
                            end do
                        end do
                    end do

                    if (n > 0) then
                        $:GPU_PARALLEL_LOOP(collapse=3)
                        do i = momxb, momxe
                            do l = isz%beg, isz%end
                                do k = isy%beg, isy%end

                                    dqL_prim_dy_vf(i)%sf(-1, k, l) = &
                                        dqR_prim_dy_vf(i)%sf(0, k, l)
                                end do
                            end do
                        end do

                        if (p > 0) then
                            $:GPU_PARALLEL_LOOP(collapse=3)
                            do i = momxb, momxe
                                do l = isz%beg, isz%end
                                    do k = isy%beg, isy%end

                                        dqL_prim_dz_vf(i)%sf(-1, k, l) = &
                                            dqR_prim_dz_vf(i)%sf(0, k, l)
                                    end do
                                end do
                            end do
                        end if

                    end if

                end if

            end if

            if (bc_x%end == BC_RIEMANN_EXTRAP) then    ! Riemann state extrap. BC at end

                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            qR_prim_rsx_vf(m + 1, k, l, i) = &
                                qL_prim_rsx_vf(m, k, l, i)
                        end do
                    end do
                end do

                if (viscous) then

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = momxb, momxe
                        do l = isz%beg, isz%end
                            do k = isy%beg, isy%end

                                dqR_prim_dx_vf(i)%sf(m + 1, k, l) = &
                                    dqL_prim_dx_vf(i)%sf(m, k, l)
                            end do
                        end do
                    end do

                    if (n > 0) then
                        $:GPU_PARALLEL_LOOP(collapse=3)
                        do i = momxb, momxe
                            do l = isz%beg, isz%end
                                do k = isy%beg, isy%end

                                    dqR_prim_dy_vf(i)%sf(m + 1, k, l) = &
                                        dqL_prim_dy_vf(i)%sf(m, k, l)
                                end do
                            end do
                        end do

                        if (p > 0) then
                            $:GPU_PARALLEL_LOOP(collapse=3)
                            do i = momxb, momxe
                                do l = isz%beg, isz%end
                                    do k = isy%beg, isy%end

                                        dqR_prim_dz_vf(i)%sf(m + 1, k, l) = &
                                            dqL_prim_dz_vf(i)%sf(m, k, l)
                                    end do
                                end do
                            end do
                        end if

                    end if

                end if

            end if
            ! END: Population of Buffers in x-direction

            ! Population of Buffers in y-direction
        elseif (norm_dir == 2) then

            if (bc_y%beg == BC_RIEMANN_EXTRAP) then    ! Riemann state extrap. BC at beginning
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            qL_prim_rsy_vf(-1, k, l, i) = &
                                qR_prim_rsy_vf(0, k, l, i)
                        end do
                    end do
                end do

                if (viscous) then

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = momxb, momxe
                        do l = isz%beg, isz%end
                            do j = isx%beg, isx%end
                                dqL_prim_dx_vf(i)%sf(j, -1, l) = &
                                    dqR_prim_dx_vf(i)%sf(j, 0, l)
                            end do
                        end do
                    end do

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = momxb, momxe
                        do l = isz%beg, isz%end
                            do j = isx%beg, isx%end
                                dqL_prim_dy_vf(i)%sf(j, -1, l) = &
                                    dqR_prim_dy_vf(i)%sf(j, 0, l)
                            end do
                        end do
                    end do

                    if (p > 0) then
                        $:GPU_PARALLEL_LOOP(collapse=3)
                        do i = momxb, momxe
                            do l = isz%beg, isz%end
                                do j = isx%beg, isx%end
                                    dqL_prim_dz_vf(i)%sf(j, -1, l) = &
                                        dqR_prim_dz_vf(i)%sf(j, 0, l)
                                end do
                            end do
                        end do
                    end if

                end if

            end if

            if (bc_y%end == BC_RIEMANN_EXTRAP) then    ! Riemann state extrap. BC at end

                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            qR_prim_rsy_vf(n + 1, k, l, i) = &
                                qL_prim_rsy_vf(n, k, l, i)
                        end do
                    end do
                end do

                if (viscous) then

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = momxb, momxe
                        do l = isz%beg, isz%end
                            do j = isx%beg, isx%end
                                dqR_prim_dx_vf(i)%sf(j, n + 1, l) = &
                                    dqL_prim_dx_vf(i)%sf(j, n, l)
                            end do
                        end do
                    end do

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = momxb, momxe
                        do l = isz%beg, isz%end
                            do j = isx%beg, isx%end
                                dqR_prim_dy_vf(i)%sf(j, n + 1, l) = &
                                    dqL_prim_dy_vf(i)%sf(j, n, l)
                            end do
                        end do
                    end do

                    if (p > 0) then
                        $:GPU_PARALLEL_LOOP(collapse=3)
                        do i = momxb, momxe
                            do l = isz%beg, isz%end
                                do j = isx%beg, isx%end
                                    dqR_prim_dz_vf(i)%sf(j, n + 1, l) = &
                                        dqL_prim_dz_vf(i)%sf(j, n, l)
                                end do
                            end do
                        end do
                    end if

                end if

            end if
            ! END: Population of Buffers in y-direction

            ! Population of Buffers in z-direction
        else

            if (bc_z%beg == BC_RIEMANN_EXTRAP) then    ! Riemann state extrap. BC at beginning
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            qL_prim_rsz_vf(-1, k, l, i) = &
                                qR_prim_rsz_vf(0, k, l, i)
                        end do
                    end do
                end do

                if (viscous) then
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = momxb, momxe
                        do k = isy%beg, isy%end
                            do j = isx%beg, isx%end
                                dqL_prim_dx_vf(i)%sf(j, k, -1) = &
                                    dqR_prim_dx_vf(i)%sf(j, k, 0)
                            end do
                        end do
                    end do
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = momxb, momxe
                        do k = isy%beg, isy%end
                            do j = isx%beg, isx%end
                                dqL_prim_dy_vf(i)%sf(j, k, -1) = &
                                    dqR_prim_dy_vf(i)%sf(j, k, 0)
                            end do
                        end do
                    end do
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = momxb, momxe
                        do k = isy%beg, isy%end
                            do j = isx%beg, isx%end
                                dqL_prim_dz_vf(i)%sf(j, k, -1) = &
                                    dqR_prim_dz_vf(i)%sf(j, k, 0)
                            end do
                        end do
                    end do
                end if

            end if

            if (bc_z%end == BC_RIEMANN_EXTRAP) then    ! Riemann state extrap. BC at end

                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            qR_prim_rsz_vf(p + 1, k, l, i) = &
                                qL_prim_rsz_vf(p, k, l, i)
                        end do
                    end do
                end do

                if (viscous) then
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = momxb, momxe
                        do k = isy%beg, isy%end
                            do j = isx%beg, isx%end
                                dqR_prim_dx_vf(i)%sf(j, k, p + 1) = &
                                    dqL_prim_dx_vf(i)%sf(j, k, p)
                            end do
                        end do
                    end do

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = momxb, momxe
                        do k = isy%beg, isy%end
                            do j = isx%beg, isx%end
                                dqR_prim_dy_vf(i)%sf(j, k, p + 1) = &
                                    dqL_prim_dy_vf(i)%sf(j, k, p)
                            end do
                        end do
                    end do

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = momxb, momxe
                        do k = isy%beg, isy%end
                            do j = isx%beg, isx%end
                                dqR_prim_dz_vf(i)%sf(j, k, p + 1) = &
                                    dqL_prim_dz_vf(i)%sf(j, k, p)
                            end do
                        end do
                    end do
                end if

            end if

        end if
        ! END: Population of Buffers in z-direction

    end subroutine s_populate_riemann_states_variables_buffers

    !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures needed to configure the chosen Riemann
        !!      solver algorithm.
        !!  @param qL_prim_vf The  left WENO-reconstructed cell-boundary values of the
        !!      cell-average primitive variables
        !!  @param qR_prim_vf The right WENO-reconstructed cell-boundary values of the
        !!      cell-average primitive variables
        !!  @param flux_vf Intra-cell fluxes
        !!  @param flux_src_vf Intra-cell fluxes sources
        !!  @param flux_gsrc_vf Intra-cell geometric fluxes sources
        !!  @param norm_dir Dir. splitting direction
        !!  @param ix Index bounds in the x-dir
        !!  @param iy Index bounds in the y-dir
        !!  @param iz Index bounds in the z-dir
        !!  @param q_prim_vf Cell-averaged primitive variables
    subroutine s_initialize_riemann_solver( &
        flux_src_vf, &
        norm_dir)

        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: flux_src_vf

        integer, intent(in) :: norm_dir

        integer :: i, j, k, l ! Generic loop iterators

        ! Reshaping Inputted Data in x-direction

        if (norm_dir == 1) then

            if (viscous .or. (surface_tension)) then

                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = momxb, E_idx
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end
                                flux_src_vf(i)%sf(j, k, l) = 0._wp
                            end do
                        end do
                    end do
                end do
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
            end if

            ! Reshaping Inputted Data in y-direction
        elseif (norm_dir == 2) then

            if (viscous .or. (surface_tension)) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = momxb, E_idx
                    do l = is3%beg, is3%end
                        do j = is1%beg, is1%end
                            do k = is2%beg, is2%end
                                flux_src_vf(i)%sf(k, j, l) = 0._wp
                            end do
                        end do
                    end do
                end do
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
            end if

            ! Reshaping Inputted Data in z-direction
        else

            if (viscous .or. (surface_tension)) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = momxb, E_idx
                    do j = is1%beg, is1%end
                        do k = is2%beg, is2%end
                            do l = is3%beg, is3%end
                                flux_src_vf(i)%sf(l, k, j) = 0._wp
                            end do
                        end do
                    end do
                end do
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
            end if

        end if

    end subroutine s_initialize_riemann_solver

    !> @brief Computes cylindrical viscous source flux contributions for momentum and energy.
        !! Calculates Cartesian components of the stress tensor using averaged velocity derivatives
        !! and cylindrical geometric factors, then updates `flux_src_vf`.
        !! Assumes x-dir is axial (z_cyl), y-dir is radial (r_cyl), z-dir is azimuthal (theta_cyl for derivatives).
        !! @param[in] velL_vf Left boundary velocity ($v_x, v_y, v_z$) (num_dims scalar_field).
        !! @param[in] dvelL_dx_vf Left boundary $\partial v_i/\partial x$ (num_dims scalar_field).
        !! @param[in] dvelL_dy_vf Left boundary $\partial v_i/\partial y$ (num_dims scalar_field).
        !! @param[in] dvelL_dz_vf Left boundary $\partial v_i/\partial z$ (num_dims scalar_field).
        !! @param[in] velR_vf Right boundary velocity ($v_x, v_y, v_z$) (num_dims scalar_field).
        !! @param[in] dvelR_dx_vf Right boundary $\partial v_i/\partial x$ (num_dims scalar_field).
        !! @param[in] dvelR_dy_vf Right boundary $\partial v_i/\partial y$ (num_dims scalar_field).
        !! @param[in] dvelR_dz_vf Right boundary $\partial v_i/\partial z$ (num_dims scalar_field).
        !! @param[inout] flux_src_vf Intercell source flux array to update (sys_size scalar_field).
        !! @param[in] norm_dir Interface normal direction (1=x-face, 2=y-face, 3=z-face).
        !! @param[in] ix Global X-direction loop bounds (int_bounds_info).
        !! @param[in] iy Global Y-direction loop bounds (int_bounds_info).
        !! @param[in] iz Global Z-direction loop bounds (int_bounds_info).
    pure subroutine s_compute_cylindrical_viscous_source_flux(velL_vf, &
                                                              dvelL_dx_vf, dvelL_dy_vf, dvelL_dz_vf, &
                                                              velR_vf, &
                                                              dvelR_dx_vf, dvelR_dy_vf, dvelR_dz_vf, &
                                                              flux_src_vf, norm_dir, ix, iy, iz)

        type(scalar_field), dimension(num_dims), intent(in) :: velL_vf, velR_vf
        type(scalar_field), dimension(num_dims), intent(in) :: dvelL_dx_vf, dvelR_dx_vf
        type(scalar_field), dimension(num_dims), intent(in) :: dvelL_dy_vf, dvelR_dy_vf
        type(scalar_field), dimension(num_dims), intent(in) :: dvelL_dz_vf, dvelR_dz_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: flux_src_vf
        integer, intent(in) :: norm_dir
        type(int_bounds_info), intent(in) :: ix, iy, iz

        ! Local variables
        real(wp), dimension(num_dims) :: avg_v_int       !!< Averaged interface velocity $(v_x, v_y, v_z)$ (grid directions).
        real(wp), dimension(num_dims) :: avg_dvdx_int    !!< Averaged interface $\partial v_i/\partial x$ (grid dir 1).
        real(wp), dimension(num_dims) :: avg_dvdy_int    !!< Averaged interface $\partial v_i/\partial y$ (grid dir 2).
        real(wp), dimension(num_dims) :: avg_dvdz_int    !!< Averaged interface $\partial v_i/\partial z$ (grid dir 3).

        real(wp), dimension(num_dims) :: stress_vector_shear !!< Shear stress vector $(\sigma_{N1}, \sigma_{N2}, \sigma_{N3})$ on N-face (grid directions).
        real(wp) :: stress_normal_bulk  !!< Normal bulk stress component $\sigma_{NN}$ on N-face.

        real(wp) :: Re_s, Re_b        !!< Effective interface shear and bulk Reynolds numbers.
        real(wp), dimension(num_dims) :: vel_src_int !!< Interface velocity $(v_1,v_2,v_3)$ (grid directions) for viscous work.
        real(wp) :: r_eff             !!< Effective radius at interface for cylindrical terms.
        real(wp) :: div_v_term_const  !!< Common term $-(2/3)(\nabla \cdot \mathbf{v}) / \text{Re}_s$ for shear stress diagonal.
        real(wp) :: divergence_cyl    !!< Full divergence $\nabla \cdot \mathbf{v}$ in cylindrical coordinates.

        integer :: j, k, l           !!< Loop iterators for $x, y, z$ grid directions.
        integer :: i_vel             !!< Loop iterator for velocity components.
        integer :: idx_rp(3)         !!< Indices $(j,k,l)$ of 'right' point for averaging.

        $:GPU_PARALLEL_LOOP(collapse=3, private='[idx_rp, avg_v_int, &
            & avg_dvdx_int, avg_dvdy_int, avg_dvdz_int, Re_s, Re_b, &
            & vel_src_int, r_eff, divergence_cyl, stress_vector_shear, &
            & stress_normal_bulk, div_v_term_const]')
        do l = iz%beg, iz%end
            do k = iy%beg, iy%end
                do j = ix%beg, ix%end

                    ! Determine indices for the 'right' state for averaging across the interface
                    idx_rp = [j, k, l]
                    idx_rp(norm_dir) = idx_rp(norm_dir) + 1

                    ! Average velocities and their derivatives at the interface
                    ! For cylindrical: x-dir ~ axial (z_cyl), y-dir ~ radial (r_cyl), z-dir ~ azimuthal (theta_cyl)
                    $:GPU_LOOP(parallelism='[seq]')
                    do i_vel = 1, num_dims
                        avg_v_int(i_vel) = 0.5_wp*(velL_vf(i_vel)%sf(j, k, l) + velR_vf(i_vel)%sf(idx_rp(1), idx_rp(2), idx_rp(3)))

                        avg_dvdx_int(i_vel) = 0.5_wp*(dvelL_dx_vf(i_vel)%sf(j, k, l) + &
                                                      dvelR_dx_vf(i_vel)%sf(idx_rp(1), idx_rp(2), idx_rp(3)))
                        if (num_dims > 1) then
                            avg_dvdy_int(i_vel) = 0.5_wp*(dvelL_dy_vf(i_vel)%sf(j, k, l) + &
                                                          dvelR_dy_vf(i_vel)%sf(idx_rp(1), idx_rp(2), idx_rp(3)))
                        else
                            avg_dvdy_int(i_vel) = 0.0_wp
                        end if
                        if (num_dims > 2) then
                            avg_dvdz_int(i_vel) = 0.5_wp*(dvelL_dz_vf(i_vel)%sf(j, k, l) + &
                                                          dvelR_dz_vf(i_vel)%sf(idx_rp(1), idx_rp(2), idx_rp(3)))
                        else
                            avg_dvdz_int(i_vel) = 0.0_wp
                        end if
                    end do

                    ! Get Re numbers and interface velocity for viscous work
                    select case (norm_dir)
                    case (1) ! x-face (axial face in z_cyl direction)
                        Re_s = Re_avg_rsx_vf(j, k, l, 1)
                        Re_b = Re_avg_rsx_vf(j, k, l, 2)
                        vel_src_int = vel_src_rsx_vf(j, k, l, 1:num_dims)
                        r_eff = y_cc(k)
                    case (2) ! y-face (radial face in r_cyl direction)
                        Re_s = Re_avg_rsy_vf(k, j, l, 1)
                        Re_b = Re_avg_rsy_vf(k, j, l, 2)
                        vel_src_int = vel_src_rsy_vf(k, j, l, 1:num_dims)
                        r_eff = y_cb(k)
                    case (3) ! z-face (azimuthal face in theta_cyl direction)
                        Re_s = Re_avg_rsz_vf(l, k, j, 1)
                        Re_b = Re_avg_rsz_vf(l, k, j, 2)
                        vel_src_int = vel_src_rsz_vf(l, k, j, 1:num_dims)
                        r_eff = y_cc(k)
                    end select

                    ! Divergence in cylindrical coordinates (vx=vz_cyl, vy=vr_cyl, vz=vtheta_cyl)
                    divergence_cyl = avg_dvdx_int(1) + avg_dvdy_int(2) + avg_v_int(2)/r_eff
                    if (num_dims > 2) then
                        divergence_cyl = divergence_cyl + avg_dvdz_int(3)/r_eff
                    end if

                    stress_vector_shear = 0.0_wp
                    stress_normal_bulk = 0.0_wp

                    if (shear_stress) then
                        div_v_term_const = -(2.0_wp/3.0_wp)*divergence_cyl/Re_s

                        select case (norm_dir)
                        case (1) ! X-face (axial normal, z_cyl)
                            stress_vector_shear(1) = (2.0_wp*avg_dvdx_int(1))/Re_s + div_v_term_const
                            if (num_dims > 1) then
                                stress_vector_shear(2) = (avg_dvdy_int(1) + avg_dvdx_int(2))/Re_s
                            end if
                            if (num_dims > 2) then
                                stress_vector_shear(3) = (avg_dvdz_int(1)/r_eff + avg_dvdx_int(3))/Re_s
                            end if
                        case (2) ! Y-face (radial normal, r_cyl)
                            if (num_dims > 1) then
                                stress_vector_shear(1) = (avg_dvdy_int(1) + avg_dvdx_int(2))/Re_s
                                stress_vector_shear(2) = (2.0_wp*avg_dvdy_int(2))/Re_s + div_v_term_const
                                if (num_dims > 2) then
                                    stress_vector_shear(3) = (avg_dvdz_int(2)/r_eff - avg_v_int(3)/r_eff + avg_dvdy_int(3))/Re_s
                                end if
                            else
                                stress_vector_shear(1) = (2.0_wp*avg_dvdx_int(1))/Re_s + div_v_term_const
                            end if
                        case (3) ! Z-face (azimuthal normal, theta_cyl)
                            if (num_dims > 2) then
                                stress_vector_shear(1) = (avg_dvdz_int(1)/r_eff + avg_dvdx_int(3))/Re_s
                                stress_vector_shear(2) = (avg_dvdz_int(2)/r_eff - avg_v_int(3)/r_eff + avg_dvdy_int(3))/Re_s
                                stress_vector_shear(3) = (2.0_wp*(avg_dvdz_int(3)/r_eff + avg_v_int(2)/r_eff))/Re_s + div_v_term_const
                            end if
                        end select

                        $:GPU_LOOP(parallelism='[seq]')
                        do i_vel = 1, num_dims
                            flux_src_vf(momxb + i_vel - 1)%sf(j, k, l) = flux_src_vf(momxb + i_vel - 1)%sf(j, k, l) - stress_vector_shear(i_vel)
                            flux_src_vf(E_idx)%sf(j, k, l) = flux_src_vf(E_idx)%sf(j, k, l) - vel_src_int(i_vel)*stress_vector_shear(i_vel)
                        end do
                    end if

                    if (bulk_stress) then
                        stress_normal_bulk = divergence_cyl/Re_b

                        flux_src_vf(momxb + norm_dir - 1)%sf(j, k, l) = flux_src_vf(momxb + norm_dir - 1)%sf(j, k, l) - stress_normal_bulk
                        flux_src_vf(E_idx)%sf(j, k, l) = flux_src_vf(E_idx)%sf(j, k, l) - vel_src_int(norm_dir)*stress_normal_bulk
                    end if

                end do
            end do
        end do

    end subroutine s_compute_cylindrical_viscous_source_flux

    !> @brief Computes Cartesian viscous source flux contributions for momentum and energy.
    !! Calculates averaged velocity gradients, gets Re and interface velocities,
    !! calls helpers for shear/bulk stress, then updates `flux_src_vf`.
    !! @param[in] velL_vf Left boundary velocity (num_dims scalar_field).
    !! @param[in] dvelL_dx_vf Left boundary d(vel)/dx (num_dims scalar_field).
    !! @param[in] dvelL_dy_vf Left boundary d(vel)/dy (num_dims scalar_field).
    !! @param[in] dvelL_dz_vf Left boundary d(vel)/dz (num_dims scalar_field).
    !! @param[in] velR_vf Right boundary velocity (num_dims scalar_field).
    !! @param[in] dvelR_dx_vf Right boundary d(vel)/dx (num_dims scalar_field).
    !! @param[in] dvelR_dy_vf Right boundary d(vel)/dy (num_dims scalar_field).
    !! @param[in] dvelR_dz_vf Right boundary d(vel)/dz (num_dims scalar_field).
    !! @param[inout] flux_src_vf Intercell source flux array to update (sys_size scalar_field).
    !! @param[in] norm_dir Interface normal direction (1=x, 2=y, 3=z).
    !! @param[in] ix X-direction loop bounds (int_bounds_info).
    !! @param[in] iy Y-direction loop bounds (int_bounds_info).
    !! @param[in] iz Z-direction loop bounds (int_bounds_info).
    pure subroutine s_compute_cartesian_viscous_source_flux(dvelL_dx_vf, &
                                                            dvelL_dy_vf, &
                                                            dvelL_dz_vf, &
                                                            dvelR_dx_vf, &
                                                            dvelR_dy_vf, &
                                                            dvelR_dz_vf, &
                                                            flux_src_vf, &
                                                            norm_dir)

        ! Arguments
        type(scalar_field), dimension(num_dims), intent(in) :: dvelL_dx_vf, dvelR_dx_vf
        type(scalar_field), dimension(num_dims), intent(in) :: dvelL_dy_vf, dvelR_dy_vf
        type(scalar_field), dimension(num_dims), intent(in) :: dvelL_dz_vf, dvelR_dz_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: flux_src_vf
        integer, intent(in) :: norm_dir

        ! Local variables
        real(wp), dimension(num_dims, num_dims) :: vel_grad_avg        !< Averaged velocity gradient tensor `d(vel_i)/d(coord_j)`.
        real(wp), dimension(num_dims, num_dims) :: current_tau_shear   !< Current shear stress tensor.
        real(wp), dimension(num_dims, num_dims) :: current_tau_bulk    !< Current bulk stress tensor.
        real(wp), dimension(num_dims) :: vel_src_at_interface         !< Interface velocities (u,v,w) for viscous work.
        integer, dimension(3) :: idx_right_phys                     !< Physical (j,k,l) indices for right state.

        real(wp) :: Re_shear !< Interface shear Reynolds number.
        real(wp) :: Re_bulk  !< Interface bulk Reynolds number.

        integer :: j_loop         !< Physical x-index loop iterator.
        integer :: k_loop         !< Physical y-index loop iterator.
        integer :: l_loop         !< Physical z-index loop iterator.
        integer :: i_dim          !< Generic dimension/component iterator.
        integer :: vel_comp_idx   !< Velocity component iterator (1=u, 2=v, 3=w).

        real(wp) :: divergence_v   !< Velocity divergence at interface.

        $:GPU_PARALLEL_LOOP(collapse=3, private='[idx_right_phys, vel_grad_avg, &
            & current_tau_shear, current_tau_bulk, vel_src_at_interface, &
            & Re_shear, Re_bulk, divergence_v, i_dim, vel_comp_idx]')
        do l_loop = isz%beg, isz%end
            do k_loop = isy%beg, isy%end
                do j_loop = isx%beg, isx%end

                    idx_right_phys(1) = j_loop
                    idx_right_phys(2) = k_loop
                    idx_right_phys(3) = l_loop
                    idx_right_phys(norm_dir) = idx_right_phys(norm_dir) + 1

                    vel_grad_avg = 0.0_wp
                    do vel_comp_idx = 1, num_dims
                        vel_grad_avg(vel_comp_idx, 1) = 0.5_wp*(dvelL_dx_vf(vel_comp_idx)%sf(j_loop, k_loop, l_loop) + &
                                                                dvelR_dx_vf(vel_comp_idx)%sf(idx_right_phys(1), idx_right_phys(2), idx_right_phys(3)))
                        if (num_dims > 1) then
                            vel_grad_avg(vel_comp_idx, 2) = 0.5_wp*(dvelL_dy_vf(vel_comp_idx)%sf(j_loop, k_loop, l_loop) + &
                                                                    dvelR_dy_vf(vel_comp_idx)%sf(idx_right_phys(1), idx_right_phys(2), idx_right_phys(3)))
                        end if
                        if (num_dims > 2) then
                            vel_grad_avg(vel_comp_idx, 3) = 0.5_wp*(dvelL_dz_vf(vel_comp_idx)%sf(j_loop, k_loop, l_loop) + &
                                                                    dvelR_dz_vf(vel_comp_idx)%sf(idx_right_phys(1), idx_right_phys(2), idx_right_phys(3)))
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
                            flux_src_vf(momxb + i_dim - 1)%sf(j_loop, k_loop, l_loop) = &
                                flux_src_vf(momxb + i_dim - 1)%sf(j_loop, k_loop, l_loop) - current_tau_shear(norm_dir, i_dim)

                            flux_src_vf(E_idx)%sf(j_loop, k_loop, l_loop) = &
                                flux_src_vf(E_idx)%sf(j_loop, k_loop, l_loop) - &
                                vel_src_at_interface(i_dim)*current_tau_shear(norm_dir, i_dim)
                        end do
                    end if

                    if (bulk_stress) then
                        ! current_tau_bulk = 0.0_wp
                        call s_calculate_bulk_stress_tensor(Re_bulk, divergence_v, current_tau_bulk)

                        do i_dim = 1, num_dims
                            flux_src_vf(momxb + i_dim - 1)%sf(j_loop, k_loop, l_loop) = &
                                flux_src_vf(momxb + i_dim - 1)%sf(j_loop, k_loop, l_loop) - current_tau_bulk(norm_dir, i_dim)

                            flux_src_vf(E_idx)%sf(j_loop, k_loop, l_loop) = &
                                flux_src_vf(E_idx)%sf(j_loop, k_loop, l_loop) - &
                                vel_src_at_interface(i_dim)*current_tau_bulk(norm_dir, i_dim)
                        end do
                    end if

                end do
            end do
        end do

    end subroutine s_compute_cartesian_viscous_source_flux

    !> @brief Calculates shear stress tensor components.
    !! tau_ij_shear = ( (dui/dxj + duj/dxi) - (2/3)*(div_v)*delta_ij ) / Re_shear
    !! @param[in] vel_grad_avg Averaged velocity gradient tensor (d(vel_i)/d(coord_j)).
    !! @param[in] Re_shear Shear Reynolds number.
    !! @param[in] divergence_v Velocity divergence (du/dx + dv/dy + dw/dz).
    !! @param[out] tau_shear_out Calculated shear stress tensor (stress on i-face, j-direction).
    pure subroutine s_calculate_shear_stress_tensor(vel_grad_avg, Re_shear, divergence_v, tau_shear_out)
        $:GPU_ROUTINE(parallelism='[seq]')

        implicit none

        ! Arguments
        real(wp), dimension(num_dims, num_dims), intent(in) :: vel_grad_avg
        real(wp), intent(in) :: Re_shear
        real(wp), intent(in) :: divergence_v
        real(wp), dimension(num_dims, num_dims), intent(out) :: tau_shear_out

        ! Local variables
        integer :: i_dim !< Loop iterator for face normal.
        integer :: j_dim !< Loop iterator for force component direction.

        tau_shear_out = 0.0_wp

        do i_dim = 1, num_dims
            do j_dim = 1, num_dims
                tau_shear_out(i_dim, j_dim) = (vel_grad_avg(j_dim, i_dim) + vel_grad_avg(i_dim, j_dim))/Re_shear
                if (i_dim == j_dim) then
                    tau_shear_out(i_dim, j_dim) = tau_shear_out(i_dim, j_dim) - &
                                                  (2.0_wp/3.0_wp)*divergence_v/Re_shear
                end if
            end do
        end do

    end subroutine s_calculate_shear_stress_tensor

    !> @brief Calculates bulk stress tensor components (diagonal only).
    !! tau_ii_bulk = (div_v) / Re_bulk. Off-diagonals are zero.
    !! @param[in] Re_bulk Bulk Reynolds number.
    !! @param[in] divergence_v Velocity divergence (du/dx + dv/dy + dw/dz).
    !! @param[out] tau_bulk_out Calculated bulk stress tensor (stress on i-face, i-direction).
    pure subroutine s_calculate_bulk_stress_tensor(Re_bulk, divergence_v, tau_bulk_out)
        $:GPU_ROUTINE(parallelism='[seq]')

        implicit none

        ! Arguments
        real(wp), intent(in) :: Re_bulk
        real(wp), intent(in) :: divergence_v
        real(wp), dimension(num_dims, num_dims), intent(out) :: tau_bulk_out

        ! Local variables
        integer :: i_dim !< Loop iterator for diagonal components.

        tau_bulk_out = 0.0_wp

        do i_dim = 1, num_dims
            tau_bulk_out(i_dim, i_dim) = divergence_v/Re_bulk
        end do

    end subroutine s_calculate_bulk_stress_tensor

    !>  Deallocation and/or disassociation procedures that are
        !!      needed to finalize the selected Riemann problem solver
        !!  @param flux_vf       Intercell fluxes
        !!  @param flux_src_vf   Intercell source fluxes
        !!  @param flux_gsrc_vf  Intercell geometric source fluxes
        !!  @param norm_dir Dimensional splitting coordinate direction
    pure subroutine s_finalize_riemann_solver(flux_vf, flux_src_vf, &
                                              flux_gsrc_vf, &
                                              norm_dir)

        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: flux_vf, flux_src_vf, flux_gsrc_vf

        integer, intent(in) :: norm_dir

        integer :: i, j, k, l !< Generic loop iterators

        ! Reshaping Outputted Data in y-direction
        if (norm_dir == 2) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do l = is3%beg, is3%end
                    do j = is1%beg, is1%end
                        do k = is2%beg, is2%end
                            flux_vf(i)%sf(k, j, l) = &
                                flux_rsy_vf(j, k, l, i)
                        end do
                    end do
                end do
            end do

            if (cyl_coord) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do j = is1%beg, is1%end
                            do k = is2%beg, is2%end
                                flux_gsrc_vf(i)%sf(k, j, l) = &
                                    flux_gsrc_rsy_vf(j, k, l, i)
                            end do
                        end do
                    end do
                end do
            end if

            $:GPU_PARALLEL_LOOP(collapse=3)
            do l = is3%beg, is3%end
                do j = is1%beg, is1%end
                    do k = is2%beg, is2%end
                        flux_src_vf(advxb)%sf(k, j, l) = &
                            flux_src_rsy_vf(j, k, l, advxb)
                    end do
                end do
            end do

            if (riemann_solver == 1 .or. riemann_solver == 4) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = advxb + 1, advxe
                    do l = is3%beg, is3%end
                        do j = is1%beg, is1%end
                            do k = is2%beg, is2%end
                                flux_src_vf(i)%sf(k, j, l) = &
                                    flux_src_rsy_vf(j, k, l, i)
                            end do
                        end do
                    end do
                end do

            end if
            ! Reshaping Outputted Data in z-direction
        elseif (norm_dir == 3) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do j = is1%beg, is1%end
                    do k = is2%beg, is2%end
                        do l = is3%beg, is3%end

                            flux_vf(i)%sf(l, k, j) = &
                                flux_rsz_vf(j, k, l, i)
                        end do
                    end do
                end do
            end do
            if (grid_geometry == 3) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = 1, sys_size
                    do j = is1%beg, is1%end
                        do k = is2%beg, is2%end
                            do l = is3%beg, is3%end

                                flux_gsrc_vf(i)%sf(l, k, j) = &
                                    flux_gsrc_rsz_vf(j, k, l, i)
                            end do
                        end do
                    end do
                end do
            end if

            $:GPU_PARALLEL_LOOP(collapse=3)
            do j = is1%beg, is1%end
                do k = is2%beg, is2%end
                    do l = is3%beg, is3%end
                        flux_src_vf(advxb)%sf(l, k, j) = &
                            flux_src_rsz_vf(j, k, l, advxb)
                    end do
                end do
            end do

            if (riemann_solver == 1 .or. riemann_solver == 4) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = advxb + 1, advxe
                    do j = is1%beg, is1%end
                        do k = is2%beg, is2%end
                            do l = is3%beg, is3%end
                                flux_src_vf(i)%sf(l, k, j) = &
                                    flux_src_rsz_vf(j, k, l, i)
                            end do
                        end do
                    end do
                end do

            end if
        elseif (norm_dir == 1) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do l = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = is1%beg, is1%end
                            flux_vf(i)%sf(j, k, l) = &
                                flux_rsx_vf(j, k, l, i)
                        end do
                    end do
                end do
            end do

            $:GPU_PARALLEL_LOOP(collapse=3)
            do l = is3%beg, is3%end
                do k = is2%beg, is2%end
                    do j = is1%beg, is1%end
                        flux_src_vf(advxb)%sf(j, k, l) = &
                            flux_src_rsx_vf(j, k, l, advxb)
                    end do
                end do
            end do

            if (riemann_solver == 1 .or. riemann_solver == 4) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = advxb + 1, advxe
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end
                                flux_src_vf(i)%sf(j, k, l) = &
                                    flux_src_rsx_vf(j, k, l, i)
                            end do
                        end do
                    end do
                end do
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

        if (n == 0) return

        if (viscous) then
            @:DEALLOCATE(Re_avg_rsy_vf)
        end if
        @:DEALLOCATE(vel_src_rsy_vf)
        @:DEALLOCATE(flux_rsy_vf)
        @:DEALLOCATE(flux_src_rsy_vf)
        @:DEALLOCATE(flux_gsrc_rsy_vf)
        if (qbmm) then
            @:DEALLOCATE(mom_sp_rsy_vf)
        end if

        if (p == 0) return

        if (viscous) then
            @:DEALLOCATE(Re_avg_rsz_vf)
        end if
        @:DEALLOCATE(vel_src_rsz_vf)
        @:DEALLOCATE(flux_rsz_vf)
        @:DEALLOCATE(flux_src_rsz_vf)
        @:DEALLOCATE(flux_gsrc_rsz_vf)
        if (qbmm) then
            @:DEALLOCATE(mom_sp_rsz_vf)
        end if

    end subroutine s_finalize_riemann_solvers_module

end module m_riemann_solvers
