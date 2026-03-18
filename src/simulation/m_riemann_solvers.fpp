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

    use m_surface_tension      !< To get the capilary fluxes

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
 s_hypo_hlld_riemann_solver, &
 s_finalize_nc_iface_vel, &
 s_finalize_riemann_solvers_module

    !> The cell-boundary values of the fluxes (src - source) that are computed
    !! through the chosen Riemann problem solver, and the direct evaluation of
    !! source terms, by using the left and right states given in qK_prim_rs_vf,
    !! dqK_prim_ds_vf where ds = dx, dy or dz.
    !> @{

    real(wp), allocatable, dimension(:, :, :, :) :: flux_rsx_vf, flux_src_rsx_vf
    real(wp), allocatable, dimension(:, :, :, :) :: flux_rsy_vf, flux_src_rsy_vf
    real(wp), allocatable, dimension(:, :, :, :) :: flux_rsz_vf, flux_src_rsz_vf
    !$acc declare create( flux_rsx_vf, flux_src_rsx_vf, flux_rsy_vf,  &
    !$acc   flux_src_rsy_vf, flux_rsz_vf, flux_src_rsz_vf )
    !> @}

    !> The cell-boundary values of the geometrical source flux that are computed
    !! through the chosen Riemann problem solver by using the left and right
    !! states given in qK_prim_rs_vf. Currently 2D axisymmetric for inviscid only.
    !> @{

    real(wp), allocatable, dimension(:, :, :, :) :: flux_gsrc_rsx_vf !<
    real(wp), allocatable, dimension(:, :, :, :) :: flux_gsrc_rsy_vf !<
    real(wp), allocatable, dimension(:, :, :, :) :: flux_gsrc_rsz_vf !<
    !$acc declare create( flux_gsrc_rsx_vf, flux_gsrc_rsy_vf, flux_gsrc_rsz_vf )

    real(wp), allocatable, dimension(:, :, :, :) :: nc_iface_vel_rsx_vf
    real(wp), allocatable, dimension(:, :, :, :) :: nc_iface_vel_rsy_vf
    real(wp), allocatable, dimension(:, :, :, :) :: nc_iface_vel_rsz_vf
    !$acc declare create(nc_iface_vel_rsx_vf, nc_iface_vel_rsy_vf, nc_iface_vel_rsz_vf)
    !> @}

    ! The cell-boundary values of the velocity. vel_src_rs_vf is determined as
    ! part of Riemann problem solution and is used to evaluate the source flux.

    real(wp), allocatable, dimension(:, :, :, :) :: vel_src_rsx_vf
    real(wp), allocatable, dimension(:, :, :, :) :: vel_src_rsy_vf
    real(wp), allocatable, dimension(:, :, :, :) :: vel_src_rsz_vf
    !$acc declare create(vel_src_rsx_vf, vel_src_rsy_vf, vel_src_rsz_vf)

    real(wp), allocatable, dimension(:, :, :, :) :: mom_sp_rsx_vf
    real(wp), allocatable, dimension(:, :, :, :) :: mom_sp_rsy_vf
    real(wp), allocatable, dimension(:, :, :, :) :: mom_sp_rsz_vf
    !$acc declare create(mom_sp_rsx_vf, mom_sp_rsy_vf, mom_sp_rsz_vf)

    real(wp), allocatable, dimension(:, :, :, :) :: Re_avg_rsx_vf
    real(wp), allocatable, dimension(:, :, :, :) :: Re_avg_rsy_vf
    real(wp), allocatable, dimension(:, :, :, :) :: Re_avg_rsz_vf
    !$acc declare create(Re_avg_rsx_vf, Re_avg_rsy_vf, Re_avg_rsz_vf)

    !> @name Indical bounds in the s1-, s2- and s3-directions
    !> @{
    type(int_bounds_info) :: is1, is2, is3
    type(int_bounds_info) :: isx, isy, isz
    !> @}

    !$acc declare create(is1, is2, is3, isx, isy, isz)

    real(wp), allocatable, dimension(:) :: Gs
    !$acc declare create(Gs)

    real(wp), allocatable, dimension(:, :) :: Res
    !$acc declare create(Res)

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
                                norm_dir, ix, iy, iz, is_hat_L)

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

        logical, intent(in) :: is_hat_L

        type(int_bounds_info), intent(IN) :: ix, iy, iz

        if (hypo_nc_dual_pass) then
            call s_hypo_hlld_riemann_solver(qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, dqL_prim_dx_vf, &
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
                                            norm_dir, ix, iy, iz, is_hat_L)

            return
        end if

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
    subroutine s_compute_viscous_source_flux(velL_vf, &
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
            call s_compute_cartesian_viscous_source_flux(velL_vf, &
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

        real(wp) :: rho_HLL_ex, rhou_n_HLL_ex, rhou_t_HLL_ex, rhou_t2_HLL_ex
        real(wp) :: u_n_iface, u_t_iface, u_t2_iface
        real(wp) :: pTot_L_ex, pTot_R_ex

        integer :: i, j, k, l, q !< Generic loop iterators

        ! Populating the buffers of the left and right Riemann problem
        ! states variables, based on the choice of boundary conditions
        call s_populate_riemann_states_variables_buffers( &
            qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, dqL_prim_dx_vf, &
            dqL_prim_dy_vf, &
            dqL_prim_dz_vf, &
            qL_prim_vf, &
            qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, dqR_prim_dx_vf, &
            dqR_prim_dy_vf, &
            dqR_prim_dz_vf, &
            qR_prim_vf, &
            norm_dir, ix, iy, iz)

        ! Reshaping inputted data based on dimensional splitting direction
        call s_initialize_riemann_solver( &
            q_prim_vf, &
            flux_vf, flux_src_vf, &
            flux_gsrc_vf, &
            norm_dir, ix, iy, iz)
        #:for NORM_DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]

            if (norm_dir == ${NORM_DIR}$) then
                !$acc parallel loop collapse(3) gang vector default(present)    &
                !$acc private(alpha_rho_L, alpha_rho_R, vel_L, vel_R, alpha_L,  &
                !$acc alpha_R, tau_e_L, tau_e_R, G_L, G_R, Re_L, Re_R,          &
                !$acc rho_avg, h_avg, gamma_avg, s_L, s_R, s_S, Ys_L, Ys_R,     &
                !$acc xi_field_L, xi_field_R,                                   &
                !$acc Cp_iL, Cp_iR, Xs_L, Xs_R, Gamma_iL, Gamma_iR,             &
                !$acc Yi_avg, Phi_avg, h_iL, h_iR, h_avg_2,                     &
                !$acc c_fast, pres_mag, B, Ga, vdotB, B2, b4, cm,               &
                !$acc pcorr, zcoef, vel_L_tmp, vel_R_tmp, &
                !$acc rho_HLL_ex, rhou_n_HLL_ex, rhou_t_HLL_ex, rhou_t2_HLL_ex, u_n_iface, u_t_iface, u_t2_iface, pTot_L_ex, pTot_R_ex)
                do l = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = is1%beg, is1%end
                            !$acc loop seq
                            do i = 1, contxe
                                alpha_rho_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                alpha_rho_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)
                            end do

                            !$acc loop seq
                            do i = 1, num_vels
                                vel_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, contxe + i)
                                vel_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, contxe + i)
                            end do

                            vel_L_rms = 0._wp; vel_R_rms = 0._wp

                            !$acc loop seq
                            do i = 1, num_vels
                                vel_L_rms = vel_L_rms + vel_L(i)**2._wp
                                vel_R_rms = vel_R_rms + vel_R(i)**2._wp
                            end do

                            !$acc loop seq
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
                                !$acc loop seq
                                do i = 1, num_fluids
                                    alpha_rho_L(i) = max(0._wp, alpha_rho_L(i))
                                    alpha_L(i) = min(max(0._wp, alpha_L(i)), 1._wp)
                                    alpha_L_sum = alpha_L_sum + alpha_L(i)
                                end do

                                alpha_L = alpha_L/max(alpha_L_sum, sgm_eps)

                                !$acc loop seq
                                do i = 1, num_fluids
                                    alpha_rho_R(i) = max(0._wp, alpha_rho_R(i))
                                    alpha_R(i) = min(max(0._wp, alpha_R(i)), 1._wp)
                                    alpha_R_sum = alpha_R_sum + alpha_R(i)
                                end do

                                alpha_R = alpha_R/max(alpha_R_sum, sgm_eps)
                            end if

                            !$acc loop seq
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
                                !$acc loop seq
                                do i = 1, 2
                                    Re_L(i) = dflt_real

                                    if (Re_size(i) > 0) Re_L(i) = 0._wp

                                    !$acc loop seq
                                    do q = 1, Re_size(i)
                                        Re_L(i) = alpha_L(Re_idx(i, q))/Res(i, q) &
                                                  + Re_L(i)
                                    end do

                                    Re_L(i) = 1._wp/max(Re_L(i), sgm_eps)

                                end do

                                !$acc loop seq
                                do i = 1, 2
                                    Re_R(i) = dflt_real

                                    if (Re_size(i) > 0) Re_R(i) = 0._wp

                                    !$acc loop seq
                                    do q = 1, Re_size(i)
                                        Re_R(i) = alpha_R(Re_idx(i, q))/Res(i, q) &
                                                  + Re_R(i)
                                    end do

                                    Re_R(i) = 1._wp/max(Re_R(i), sgm_eps)
                                end do
                            end if

                            if (chemistry) then
                                !$acc loop seq
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

                                E_L = rho_L*E_L + 5e-1*rho_L*vel_L_rms
                                E_R = rho_R*E_R + 5e-1*rho_R*vel_R_rms
                                H_L = (E_L + pres_L)/rho_L
                                H_R = (E_R + pres_R)/rho_R
                            elseif (mhd .and. relativity) then
                                Ga%L = 1._wp/sqrt(1._wp - vel_L_rms)
                                Ga%R = 1._wp/sqrt(1._wp - vel_R_rms)
                                vdotB%L = vel_L(1)*B%L(1) + vel_L(2)*B%L(2) + vel_L(3)*B%L(3)
                                vdotB%R = vel_R(1)*B%R(1) + vel_R(2)*B%R(2) + vel_R(3)*B%R(3)

                                b4%L(1) = B%L(1)/Ga%L + Ga%L*vel_L(1)*vdotB%L
                                b4%L(2) = B%L(2)/Ga%L + Ga%L*vel_L(2)*vdotB%L
                                b4%L(3) = B%L(3)/Ga%L + Ga%L*vel_L(3)*vdotB%L
                                b4%R(1) = B%R(1)/Ga%R + Ga%R*vel_R(1)*vdotB%R
                                b4%R(2) = B%R(2)/Ga%R + Ga%R*vel_R(2)*vdotB%R
                                b4%R(3) = B%R(3)/Ga%R + Ga%R*vel_R(3)*vdotB%R
                                B2%L = B%L(1)**2._wp + B%L(2)**2._wp + B%L(3)**2._wp
                                B2%R = B%R(1)**2._wp + B%R(2)**2._wp + B%R(3)**2._wp

                                pres_mag%L = 0.5_wp*(B2%L/Ga%L**2._wp + vdotB%L**2._wp)
                                pres_mag%R = 0.5_wp*(B2%R/Ga%R**2._wp + vdotB%R**2._wp)

                                ! Hard-coded EOS
                                H_L = 1._wp + (gamma_L + 1)*pres_L/rho_L
                                H_R = 1._wp + (gamma_R + 1)*pres_R/rho_R

                                cm%L(1) = (rho_L*H_L*Ga%L**2 + B2%L)*vel_L(1) - vdotB%L*B%L(1)
                                cm%L(2) = (rho_L*H_L*Ga%L**2 + B2%L)*vel_L(2) - vdotB%L*B%L(2)
                                cm%L(3) = (rho_L*H_L*Ga%L**2 + B2%L)*vel_L(3) - vdotB%L*B%L(3)
                                cm%R(1) = (rho_R*H_R*Ga%R**2 + B2%R)*vel_R(1) - vdotB%R*B%R(1)
                                cm%R(2) = (rho_R*H_R*Ga%R**2 + B2%R)*vel_R(2) - vdotB%R*B%R(2)
                                cm%R(3) = (rho_R*H_R*Ga%R**2 + B2%R)*vel_R(3) - vdotB%R*B%R(3)

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
                                E_L = gamma_L*pres_L + pi_inf_L + 5e-1*rho_L*vel_L_rms + qv_L
                                E_R = gamma_R*pres_R + pi_inf_R + 5e-1*rho_R*vel_R_rms + qv_R
                                H_L = (E_L + pres_L)/rho_L
                                H_R = (E_R + pres_R)/rho_R
                            end if

                            ! elastic energy update
                            if (hypoelasticity) then
                                G_L = 0._wp; G_R = 0._wp

                                !$acc loop seq
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
                            !    !$acc loop seq
                            !    do i = 1, num_fluids
                            !        G_L = G_L + alpha_L(i)*Gs(i)
                            !        G_R = G_R + alpha_R(i)*Gs(i)
                            !    end do
                            !    ! Elastic contribution to energy if G large enough
                            !    if ((G_L > 1e-3._wp) .and. (G_R > 1e-3._wp)) then
                            !    E_L = E_L + G_L*qL_prim_rs${XYZ}$_vf(j, k, l, xiend + 1)
                            !    E_R = E_R + G_R*qR_prim_rs${XYZ}$_vf(j + 1, k, l, xiend + 1)
                            !    !$acc loop seq
                            !    do i = 1, b_size-1
                            !        tau_e_L(i) = G_L*qL_prim_rs${XYZ}$_vf(j, k, l, strxb - 1 + i)
                            !        tau_e_R(i) = G_R*qR_prim_rs${XYZ}$_vf(j + 1, k, l, strxb - 1 + i)
                            !    end do
                            !    !$acc loop seq
                            !    do i = 1, b_size-1
                            !        tau_e_L(i) = 0_wp
                            !        tau_e_R(i) = 0_wp
                            !    end do
                            !    !$acc loop seq
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
                                !$acc loop seq
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
                                pres_SL = 5e-1_wp*(pres_L + pres_R + rho_avg*c_avg* &
                                                   (vel_L(dir_idx(1)) - &
                                                    vel_R(dir_idx(1))))

                                pres_SR = pres_SL

                                Ms_L = max(1._wp, sqrt(1._wp + ((5e-1_wp + gamma_L)/(1._wp + gamma_L))* &
                                                       (pres_SL/pres_L - 1._wp)*pres_L/ &
                                                       ((pres_L + pi_inf_L/(1._wp + gamma_L)))))
                                Ms_R = max(1._wp, sqrt(1._wp + ((5e-1_wp + gamma_R)/(1._wp + gamma_R))* &
                                                       (pres_SR/pres_R - 1._wp)*pres_R/ &
                                                       ((pres_R + pi_inf_R/(1._wp + gamma_R)))))

                                s_L = vel_L(dir_idx(1)) - c_L*Ms_L
                                s_R = vel_R(dir_idx(1)) + c_R*Ms_R

                                s_S = 5e-1_wp*((vel_L(dir_idx(1)) + vel_R(dir_idx(1))) + &
                                               (pres_L - pres_R)/ &
                                               (rho_avg*c_avg))
                            end if

                            s_M = min(0._wp, s_L); s_P = max(0._wp, s_R)

                            xi_M = (5e-1_wp + sign(5e-1_wp, s_L)) &
                                   + (5e-1_wp - sign(5e-1_wp, s_L)) &
                                   *(5e-1_wp + sign(5e-1_wp, s_R))
                            xi_P = (5e-1_wp - sign(5e-1_wp, s_R)) &
                                   + (5e-1_wp - sign(5e-1_wp, s_L)) &
                                   *(5e-1_wp + sign(5e-1_wp, s_R))

                            ! Low Mach correction
                            if (low_Mach == 1) then
                                @:compute_low_Mach_correction()
                            else
                                pcorr = 0._wp
                            end if

                            ! Mass
                            if (.not. relativity) then
                                !$acc loop seq
                                do i = 1, contxe
                                    flux_rs${XYZ}$_vf(j, k, l, i) = &
                                        (s_M*alpha_rho_R(i)*vel_R(norm_dir) &
                                         - s_P*alpha_rho_L(i)*vel_L(norm_dir) &
                                         + s_M*s_P*(alpha_rho_L(i) &
                                                    - alpha_rho_R(i))) &
                                        /(s_M - s_P)
                                end do
                            elseif (relativity) then
                                !$acc loop seq
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
                                ! Flux of rho*v_x in the ${XYZ}$ direction
                                ! = rho * v_x * v_${XYZ}$ - B_x * B_${XYZ}$ + delta_(${XYZ}$,x) * p_tot
                                flux_rs${XYZ}$_vf(j, k, l, contxe + 1) = &
                                    (s_M*(rho_R*vel_R(1)*vel_R(norm_dir) &
                                          - B%R(1)*B%R(norm_dir) &
                                          + dir_flg(1)*(pres_R + pres_mag%R)) &
                                     - s_P*(rho_L*vel_L(1)*vel_L(norm_dir) &
                                            - B%L(1)*B%L(norm_dir) &
                                            + dir_flg(1)*(pres_L + pres_mag%L)) &
                                     + s_M*s_P*(rho_L*vel_L(1) - rho_R*vel_R(1))) &
                                    /(s_M - s_P)
                                ! Flux of rho*v_y in the ${XYZ}$ direction
                                ! = rho * v_y * v_${XYZ}$ - B_y * B_${XYZ}$ + delta_(${XYZ}$,y) * p_tot
                                flux_rs${XYZ}$_vf(j, k, l, contxe + 2) = &
                                    (s_M*(rho_R*vel_R(2)*vel_R(norm_dir) &
                                          - B%R(2)*B%R(norm_dir) &
                                          + dir_flg(2)*(pres_R + pres_mag%R)) &
                                     - s_P*(rho_L*vel_L(2)*vel_L(norm_dir) &
                                            - B%L(2)*B%L(norm_dir) &
                                            + dir_flg(2)*(pres_L + pres_mag%L)) &
                                     + s_M*s_P*(rho_L*vel_L(2) - rho_R*vel_R(2))) &
                                    /(s_M - s_P)
                                ! Flux of rho*v_z in the ${XYZ}$ direction
                                ! = rho * v_z * v_${XYZ}$ - B_z * B_${XYZ}$ + delta_(${XYZ}$,z) * p_tot
                                flux_rs${XYZ}$_vf(j, k, l, contxe + 3) = &
                                    (s_M*(rho_R*vel_R(3)*vel_R(norm_dir) &
                                          - B%R(3)*B%R(norm_dir) &
                                          + dir_flg(3)*(pres_R + pres_mag%R)) &
                                     - s_P*(rho_L*vel_L(3)*vel_L(norm_dir) &
                                            - B%L(3)*B%L(norm_dir) &
                                            + dir_flg(3)*(pres_L + pres_mag%L)) &
                                     + s_M*s_P*(rho_L*vel_L(3) - rho_R*vel_R(3))) &
                                    /(s_M - s_P)
                            elseif (mhd .and. relativity) then
                                ! Flux of m_x in the ${XYZ}$ direction
                                ! = m_x * v_${XYZ}$ - b_x/Gamma * B_${XYZ}$ + delta_(${XYZ}$,x) * p_tot
                                flux_rs${XYZ}$_vf(j, k, l, contxe + 1) = &
                                    (s_M*(cm%R(1)*vel_R(norm_dir) &
                                          - b4%R(1)/Ga%R*B%R(norm_dir) &
                                          + dir_flg(1)*(pres_R + pres_mag%R)) &
                                     - s_P*(cm%L(1)*vel_L(norm_dir) &
                                            - b4%L(1)/Ga%L*B%L(norm_dir) &
                                            + dir_flg(1)*(pres_L + pres_mag%L)) &
                                     + s_M*s_P*(cm%L(1) - cm%R(1))) &
                                    /(s_M - s_P)
                                ! Flux of m_y in the ${XYZ}$ direction
                                ! = rho * v_y * v_${XYZ}$ - B_y * B_${XYZ}$ + delta_(${XYZ}$,y) * p_tot
                                flux_rs${XYZ}$_vf(j, k, l, contxe + 2) = &
                                    (s_M*(cm%R(2)*vel_R(norm_dir) &
                                          - b4%R(2)/Ga%R*B%R(norm_dir) &
                                          + dir_flg(2)*(pres_R + pres_mag%R)) &
                                     - s_P*(cm%L(2)*vel_L(norm_dir) &
                                            - b4%L(2)/Ga%L*B%L(norm_dir) &
                                            + dir_flg(2)*(pres_L + pres_mag%L)) &
                                     + s_M*s_P*(cm%L(2) - cm%R(2))) &
                                    /(s_M - s_P)
                                ! Flux of m_z in the ${XYZ}$ direction
                                ! = rho * v_z * v_${XYZ}$ - B_z * B_${XYZ}$ + delta_(${XYZ}$,z) * p_tot
                                flux_rs${XYZ}$_vf(j, k, l, contxe + 3) = &
                                    (s_M*(cm%R(3)*vel_R(norm_dir) &
                                          - b4%R(3)/Ga%R*B%R(norm_dir) &
                                          + dir_flg(3)*(pres_R + pres_mag%R)) &
                                     - s_P*(cm%L(3)*vel_L(norm_dir) &
                                            - b4%L(3)/Ga%L*B%L(norm_dir) &
                                            + dir_flg(3)*(pres_L + pres_mag%L)) &
                                     + s_M*s_P*(cm%L(3) - cm%R(3))) &
                                    /(s_M - s_P)
                            elseif (bubbles_euler) then
                                !$acc loop seq
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
                                !$acc loop seq
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
                                !$acc loop seq
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
                                !TODO: simplify this so it's not split into 3
                                if (num_dims == 1) then
                                    flux_rs${XYZ}$_vf(j, k, l, E_idx) = &
                                        (s_M*(vel_R(dir_idx(1))*(E_R + pres_R) &
                                              - (tau_e_R(dir_idx_tau(1))*vel_R(dir_idx(1)))) &
                                         - s_P*(vel_L(dir_idx(1))*(E_L + pres_L) &
                                                - (tau_e_L(dir_idx_tau(1))*vel_L(dir_idx(1)))) &
                                         + s_M*s_P*(E_L - E_R)) &
                                        /(s_M - s_P)
                                else if (num_dims == 2) then
                                    flux_rs${XYZ}$_vf(j, k, l, E_idx) = &
                                        (s_M*(vel_R(dir_idx(1))*(E_R + pres_R) &
                                              - (tau_e_R(dir_idx_tau(1))*vel_R(dir_idx(1))) &
                                              - (tau_e_R(dir_idx_tau(2))*vel_R(dir_idx(2)))) &
                                         - s_P*(vel_L(dir_idx(1))*(E_L + pres_L) &
                                                - (tau_e_L(dir_idx_tau(1))*vel_L(dir_idx(1))) &
                                                - (tau_e_L(dir_idx_tau(2))*vel_L(dir_idx(2)))) &
                                         + s_M*s_P*(E_L - E_R)) &
                                        /(s_M - s_P)
                                else if (num_dims == 3) then
                                    flux_rs${XYZ}$_vf(j, k, l, E_idx) = &
                                        (s_M*(vel_R(dir_idx(1))*(E_R + pres_R) &
                                              - (tau_e_R(dir_idx_tau(1))*vel_R(dir_idx(1))) &
                                              - (tau_e_R(dir_idx_tau(2))*vel_R(dir_idx(2))) &
                                              - (tau_e_R(dir_idx_tau(3))*vel_R(dir_idx(3)))) &
                                         - s_P*(vel_L(dir_idx(1))*(E_L + pres_L) &
                                                - (tau_e_L(dir_idx_tau(1))*vel_L(dir_idx(1))) &
                                                - (tau_e_L(dir_idx_tau(2))*vel_L(dir_idx(2))) &
                                                - (tau_e_L(dir_idx_tau(3))*vel_L(dir_idx(3)))) &
                                         + s_M*s_P*(E_L - E_R)) &
                                        /(s_M - s_P)
                                end if
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

                            ! For velocity traces for NC terms, always use Method 2
                            ! (used by interface-consistent hypo RHS and HLL Method 1 K div u)
                            if (hypo_nc_interface .or. alt_soundspeed) then
                                ! Compute scalar HLL transport traces: F_HLL(U=1, F=u_i)
                                !$acc loop seq
                                do i = 1, num_dims
                                    if (0._wp <= s_L) then
                                        nc_iface_vel_rs${XYZ}$_vf(j, k, l, dir_idx(i)) = vel_L(dir_idx(i))
                                    else if (s_R <= 0._wp) then
                                        nc_iface_vel_rs${XYZ}$_vf(j, k, l, dir_idx(i)) = vel_R(dir_idx(i))
                                    else
                                        nc_iface_vel_rs${XYZ}$_vf(j, k, l, dir_idx(i)) = &
                                            (s_R*vel_L(dir_idx(i)) - s_L*vel_R(dir_idx(i)))/(s_R - s_L)
                                    end if
                                end do
                            end if

                            ! Advection: Method 1 / Method 2 split
                            if (hll_alpha_interface) then
                                ! Method 1: F(alpha)=0, source trace = F_HLL(U=1, F=alpha)
                                !$acc loop seq
                                do i = advxb, advxe
                                    if (0._wp <= s_L .or. s_R <= 0._wp) then
                                        flux_rs${XYZ}$_vf(j, k, l, i) = 0._wp
                                    else
                                        flux_rs${XYZ}$_vf(j, k, l, i) = &
                                            s_M*s_P*(qR_prim_rs${XYZ}$_vf(j + 1, k, l, i) &
                                                     - qL_prim_rs${XYZ}$_vf(j, k, l, i)) &
                                            /(s_R - s_L)
                                    end if
                                    if (0._wp <= s_L) then
                                        flux_src_rs${XYZ}$_vf(j, k, l, i) = qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                    else if (s_R <= 0._wp) then
                                        flux_src_rs${XYZ}$_vf(j, k, l, i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)
                                    else
                                        flux_src_rs${XYZ}$_vf(j, k, l, i) = &
                                            (s_R*qL_prim_rs${XYZ}$_vf(j, k, l, i) &
                                             - s_L*qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)) &
                                            /(s_R - s_L)
                                    end if
                                end do
                            else
                                ! Method 2: F(alpha)=alpha*u_n, source trace = F_HLL(U=1, F=u_n)
                                !$acc loop seq
                                do i = advxb, advxe
                                    if (0._wp <= s_L) then
                                        flux_rs${XYZ}$_vf(j, k, l, i) = &
                                            qL_prim_rs${XYZ}$_vf(j, k, l, i)*vel_L(dir_idx(1))
                                    else if (s_R <= 0._wp) then
                                        flux_rs${XYZ}$_vf(j, k, l, i) = &
                                            qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)*vel_R(dir_idx(1))
                                    else
                                        flux_rs${XYZ}$_vf(j, k, l, i) = &
                                            (s_R*qL_prim_rs${XYZ}$_vf(j, k, l, i)*vel_L(dir_idx(1)) &
                                             - s_L*qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)*vel_R(dir_idx(1)) &
                                             + s_L*s_R*(qR_prim_rs${XYZ}$_vf(j + 1, k, l, i) &
                                                        - qL_prim_rs${XYZ}$_vf(j, k, l, i))) &
                                            /(s_R - s_L)
                                    end if
                                end do
                                ! Store Psi_u (normal velocity trace) in advxb slot
                                if (0._wp <= s_L) then
                                    flux_src_rs${XYZ}$_vf(j, k, l, advxb) = vel_L(dir_idx(1))
                                else if (s_R <= 0._wp) then
                                    flux_src_rs${XYZ}$_vf(j, k, l, advxb) = vel_R(dir_idx(1))
                                else
                                    flux_src_rs${XYZ}$_vf(j, k, l, advxb) = &
                                        (s_R*vel_L(dir_idx(1)) - s_L*vel_R(dir_idx(1)))/(s_R - s_L)
                                end if
                            end if

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
                            !$acc loop seq
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
                                !$acc loop seq
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
                                    flux_rsx_vf(j, k, l, B_idx%beg) = (s_M*(vel_R(1)*B%R(2) - vel_R(2)*Bx0) &
                                                                       - s_P*(vel_L(1)*B%L(2) - vel_L(2)*Bx0) + s_M*s_P*(B%L(2) - B%R(2)))/(s_M - s_P)

                                    ! B_z flux = v_x * B_z - v_z * Bx0
                                    flux_rsx_vf(j, k, l, B_idx%beg + 1) = (s_M*(vel_R(1)*B%R(3) - vel_R(3)*Bx0) &
                                                                           - s_P*(vel_L(1)*B%L(3) - vel_L(3)*Bx0) + s_M*s_P*(B%L(3) - B%R(3)))/(s_M - s_P)

                                else ! 2D/3D: Bx, By, Bz /= const. but zero flux component in the same direction
                                    ! B_x d/d${XYZ}$ flux = (1 - delta(x,${XYZ}$)) * (v_${XYZ}$ * B_x - v_x * B_${XYZ}$)
                                    flux_rs${XYZ}$_vf(j, k, l, B_idx%beg) = (1 - dir_flg(1))*( &
                                                                            s_M*(vel_R(dir_idx(1))*B%R(1) - vel_R(1)*B%R(norm_dir)) - &
                                                                            s_P*(vel_L(dir_idx(1))*B%L(1) - vel_L(1)*B%L(norm_dir)) + &
                                                                            s_M*s_P*(B%L(1) - B%R(1)))/(s_M - s_P)

                                    ! B_y d/d${XYZ}$ flux = (1 - delta(y,${XYZ}$)) * (v_${XYZ}$ * B_y - v_y * B_${XYZ}$)
                                    flux_rs${XYZ}$_vf(j, k, l, B_idx%beg + 1) = (1 - dir_flg(2))*( &
                                                                                s_M*(vel_R(dir_idx(1))*B%R(2) - vel_R(2)*B%R(norm_dir)) - &
                                                                                s_P*(vel_L(dir_idx(1))*B%L(2) - vel_L(2)*B%L(norm_dir)) + &
                                                                                s_M*s_P*(B%L(2) - B%R(2)))/(s_M - s_P)

                                    ! B_z d/d${XYZ}$ flux = (1 - delta(z,${XYZ}$)) * (v_${XYZ}$ * B_z - v_z * B_${XYZ}$)
                                    flux_rs${XYZ}$_vf(j, k, l, B_idx%beg + 2) = (1 - dir_flg(3))*( &
                                                                                s_M*(vel_R(dir_idx(1))*B%R(3) - vel_R(3)*B%R(norm_dir)) - &
                                                                                s_P*(vel_L(dir_idx(1))*B%L(3) - vel_L(3)*B%L(norm_dir)) + &
                                                                                s_M*s_P*(B%L(3) - B%R(3)))/(s_M - s_P)

                                end if

                                flux_src_rs${XYZ}$_vf(j, k, l, advxb) = 0._wp
                            end if

                            #:if (NORM_DIR == 2)
                                if (cyl_coord) then
                                    !Substituting the advective flux into the inviscid geometrical source flux
                                    !$acc loop seq
                                    do i = 1, E_idx
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = flux_rs${XYZ}$_vf(j, k, l, i)
                                    end do
                                    ! Recalculating the radial momentum geometric source flux
                                    flux_gsrc_rs${XYZ}$_vf(j, k, l, contxe + 2) = &
                                        flux_rs${XYZ}$_vf(j, k, l, contxe + 2) &
                                        - (s_M*pres_R - s_P*pres_L)/(s_M - s_P)
                                    ! Geometrical source of the void fraction(s) is zero
                                    !$acc loop seq
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

                                    !$acc loop seq
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
                                       norm_dir, ix, iy, iz)

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

        ! --- HLLC star-state helpers / aliases (new) ---
        real(wp), dimension(sys_size) :: U_L, U_R, U_star_L, U_star_R
        real(wp), dimension(sys_size) :: F_L, F_R, F_star_L, F_star_R, F_HLLC
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
        real(wp) :: S_Mid             ! contact speed used in star formulas

        integer :: mom_n, mom_t1, mom_t2
        integer :: itnn, itnt1, itnt2, itss1, itss2, itt12

        ! --- ADC (HLL -> HLLC) ---
        real(wp), dimension(sys_size) :: F_HLL, U_HLL
        real(wp) :: u_n_HLL, u_t_HLL
        real(wp) :: u_t2_HLL
        real(wp) :: rho_HLL
        real(wp) :: p_face_HLL, tau_qq_face_HLL, tau_nn_HLL
        real(wp) :: phi

        real(wp) :: Sigma_L, Sigma_R, dSigma, Sigma_ref
        real(wp) :: a_L_ref, a_R_ref, a_ref
        real(wp) :: du_t, dtau_nt
        real(wp) :: du_t2, dtau_nt2
        real(wp) :: sensor_ptot, sensor_vt, sensor_tnt, sensor_combined

        real(wp), parameter :: ADC_power = 1.0_wp

        
        ! Populating the buffers of the left and right Riemann problem
        ! states variables, based on the choice of boundary conditions

        call s_populate_riemann_states_variables_buffers( &
            qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, dqL_prim_dx_vf, &
            dqL_prim_dy_vf, &
            dqL_prim_dz_vf, &
            qL_prim_vf, &
            qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, dqR_prim_dx_vf, &
            dqR_prim_dy_vf, &
            dqR_prim_dz_vf, &
            qR_prim_vf, &
            norm_dir, ix, iy, iz)

        ! Reshaping inputted data based on dimensional splitting direction

        call s_initialize_riemann_solver( &
            q_prim_vf, &
            flux_vf, flux_src_vf, &
            flux_gsrc_vf, &
            norm_dir, ix, iy, iz)

        idx1 = 1; if (dir_idx(1) == 2) idx1 = 2; if (dir_idx(1) == 3) idx1 = 3

        #:for NORM_DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]

            if (norm_dir == ${NORM_DIR}$) then

                ! 6-EQUATION MODEL WITH HLLC
                if (model_eqns == 3) then
                    !ME3

                    !$acc parallel loop collapse(3) gang vector default(present)                    &
                    !$acc private(vel_L, vel_R, vel_K_Star, Re_L, Re_R, rho_avg, h_avg, gamma_avg,  &
                    !$acc s_L, s_R, s_S, vel_avg_rms, alpha_L, alpha_R, Ys_L, Ys_R, Xs_L, Xs_R,     &
                    !$acc Gamma_iL, Gamma_iR, Cp_iL, Cp_iR, Yi_avg, Phi_avg, h_iL, h_iR, h_avg_2,   &
                    !$acc tau_e_L, tau_e_R, G_L, G_R, flux_ene_e, xi_field_L, xi_field_R, pcorr,    &
                    !$acc zcoef, vel_L_tmp, vel_R_tmp)
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end

                                idx1 = dir_idx(1)

                                vel_L_rms = 0._wp; vel_R_rms = 0._wp

                                !$acc loop seq
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
                                    !$acc loop seq
                                    do i = 1, num_fluids
                                        qL_prim_rs${XYZ}$_vf(j, k, l, i) = max(0._wp, qL_prim_rs${XYZ}$_vf(j, k, l, i))
                                        qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i) = min(max(0._wp, qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)), 1._wp)
                                        alpha_L_sum = alpha_L_sum + qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)
                                    end do

                                    !$acc loop seq
                                    do i = 1, num_fluids
                                        qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i) = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)/max(alpha_L_sum, sgm_eps)
                                    end do

                                    !$acc loop seq
                                    do i = 1, num_fluids
                                        qR_prim_rs${XYZ}$_vf(j + 1, k, l, i) = max(0._wp, qR_prim_rs${XYZ}$_vf(j + 1, k, l, i))
                                        qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i) = min(max(0._wp, qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)), 1._wp)
                                        alpha_R_sum = alpha_R_sum + qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)
                                    end do

                                    !$acc loop seq
                                    do i = 1, num_fluids
                                        qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)/max(alpha_R_sum, sgm_eps)
                                    end do
                                end if

                                !$acc loop seq
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
                                    !$acc loop seq
                                    do i = 1, 2
                                        Re_L(i) = dflt_real

                                        if (Re_size(i) > 0) Re_L(i) = 0._wp

                                        !$acc loop seq
                                        do q = 1, Re_size(i)
                                            Re_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + Re_idx(i, q))/Res(i, q) &
                                                      + Re_L(i)
                                        end do

                                        Re_L(i) = 1._wp/max(Re_L(i), sgm_eps)

                                    end do

                                    !$acc loop seq
                                    do i = 1, 2
                                        Re_R(i) = dflt_real

                                        if (Re_size(i) > 0) Re_R(i) = 0._wp

                                        !$acc loop seq
                                        do q = 1, Re_size(i)
                                            Re_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + Re_idx(i, q))/Res(i, q) &
                                                      + Re_R(i)
                                        end do

                                        Re_R(i) = 1._wp/max(Re_R(i), sgm_eps)
                                    end do
                                end if

                                E_L = gamma_L*pres_L + pi_inf_L + 5e-1_wp*rho_L*vel_L_rms + qv_L

                                E_R = gamma_R*pres_R + pi_inf_R + 5e-1_wp*rho_R*vel_R_rms + qv_R

                                ! ENERGY ADJUSTMENTS FOR HYPOELASTIC ENERGY
                                if (hypoelasticity) then
                                    !$acc loop seq
                                    do i = 1, strxe - strxb + 1
                                        tau_e_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, strxb - 1 + i)
                                        tau_e_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, strxb - 1 + i)
                                    end do
                                    G_L = 0_wp; G_R = 0_wp
                                    !$acc loop seq
                                    do i = 1, num_fluids
                                        G_L = G_L + alpha_L(i)*Gs(i)
                                        G_R = G_R + alpha_R(i)*Gs(i)
                                    end do
                                    !$acc loop seq
                                    do i = 1, strxe - strxb + 1
                                        ! Elastic contribution to energy if G large enough
                                        if ((G_L > verysmall) .and. (G_R > verysmall)) then
                                            E_L = E_L + (tau_e_L(i)*tau_e_L(i))/(4._wp*G_L)
                                            E_R = E_R + (tau_e_R(i)*tau_e_R(i))/(4._wp*G_R)
                                            ! Shear terms doubled: 2D/2D-axisym i==2 only; 3D i==2,4,5
                                            if ((n > 0 .and. p == 0 .and. i == 2) .or. &
                                                (p > 0 .and. (i == 2 .or. i == 4 .or. i == 5))) then
                                                E_L = E_L + (tau_e_L(i)*tau_e_L(i))/(4._wp*G_L)
                                                E_R = E_R + (tau_e_R(i)*tau_e_R(i))/(4._wp*G_R)
                                            end if
                                        end if
                                    end do
                                end if

                                ! ENERGY ADJUSTMENTS FOR HYPERELASTIC ENERGY
                                if (hyperelasticity) then
                                    !$acc loop seq
                                    do i = 1, num_dims
                                        xi_field_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, xibeg - 1 + i)
                                        xi_field_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, xibeg - 1 + i)
                                    end do
                                    G_L = 0_wp; G_R = 0_wp; 
                                    !$acc loop seq
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
                                    !$acc loop seq
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
                                    !$acc loop seq
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
                                    pres_SL = 5e-1_wp*(pres_L + pres_R + rho_avg*c_avg* &
                                                       (vel_L(dir_idx(1)) - &
                                                        vel_R(dir_idx(1))))

                                    pres_SR = pres_SL

                                    Ms_L = max(1._wp, sqrt(1._wp + ((5e-1_wp + gamma_L)/(1._wp + gamma_L))* &
                                                           (pres_SL/pres_L - 1._wp)*pres_L/ &
                                                           ((pres_L + pi_inf_L/(1._wp + gamma_L)))))
                                    Ms_R = max(1._wp, sqrt(1._wp + ((5e-1_wp + gamma_R)/(1._wp + gamma_R))* &
                                                           (pres_SR/pres_R - 1._wp)*pres_R/ &
                                                           ((pres_R + pi_inf_R/(1._wp + gamma_R)))))

                                    s_L = vel_L(dir_idx(1)) - c_L*Ms_L
                                    s_R = vel_R(dir_idx(1)) + c_R*Ms_R

                                    s_S = 5e-1_wp*((vel_L(dir_idx(1)) + vel_R(dir_idx(1))) + &
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
                                xi_M = (5e-1_wp + sign(0.5_wp, s_S))
                                xi_P = (5e-1_wp - sign(0.5_wp, s_S))

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

                                vel_K_Star = vel_L(idx1)*(1_wp - xi_MP) + xi_MP*vel_R(idx1) + &
                                             xi_MP*xi_PP*(s_S - vel_R(idx1))

                                ! Low Mach correction
                                if (low_Mach == 1) then
                                    @:compute_low_Mach_correction()
                                else
                                    pcorr = 0._wp
                                end if

                                ! COMPUTING FLUXES
                                ! MASS FLUX.
                                !$acc loop seq
                                do i = 1, contxe
                                    flux_rs${XYZ}$_vf(j, k, l, i) = &
                                        xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, i)*(vel_L(idx1) + s_M*(xi_L - 1._wp)) + &
                                        xi_P*qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)*(vel_R(idx1) + s_P*(xi_R - 1._wp))
                                end do

                                ! MOMENTUM FLUX.
                                ! f = \rho u u - \sigma, q = \rho u, q_star = \xi * \rho*(s_star, v, w)
                                !$acc loop seq
                                do i = 1, num_dims
                                    idxi = dir_idx(i)
                                    flux_rs${XYZ}$_vf(j, k, l, contxe + idxi) = rho_Star*vel_K_Star* &
                                                                                (dir_flg(idxi)*vel_K_Star + (1_wp - dir_flg(idxi))*(xi_M*vel_L(idxi) + xi_P*vel_R(idxi))) + dir_flg(idxi)*p_Star &
                                                                                + (s_M/s_L)*(s_P/s_R)*dir_flg(idxi)*pcorr
                                end do

                                ! ENERGY FLUX.
                                ! f = u*(E-\sigma), q = E, q_star = \xi*E+(s-u)(\rho s_star - \sigma/(s-u))
                                flux_rs${XYZ}$_vf(j, k, l, E_idx) = (E_star + p_Star)*vel_K_Star &
                                                                    + (s_M/s_L)*(s_P/s_R)*pcorr*s_S

                                ! ELASTICITY. Elastic shear stress additions for the momentum and energy flux
                                if (elasticity) then
                                    flux_ene_e = 0_wp; 
                                    !$acc loop seq
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
                                !$acc loop seq
                                do i = advxb, advxe
                                    flux_rs${XYZ}$_vf(j, k, l, i) = &
                                        xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, i)*s_S + &
                                        xi_P*qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)*s_S
                                end do

                                ! SOURCE TERM FOR VOLUME FRACTION ADVECTION FLUX.
                                !$acc loop seq
                                do i = 1, num_dims
                                    idxi = dir_idx(i)
                                    vel_src_rs${XYZ}$_vf(j, k, l, idxi) = &
                                        xi_M*(vel_L(idxi) + dir_flg(idxi)*(s_S*(xi_MP*(xi_L - 1) + 1) - vel_L(idxi))) + &
                                        xi_P*(vel_R(idxi) + dir_flg(idxi)*(s_S*(xi_PP*(xi_R - 1) + 1) - vel_R(idxi)))
                                end do

                                ! INTERNAL ENERGIES ADVECTION FLUX.
                                ! K-th pressure and velocity in preparation for the internal energy flux
                                !$acc loop seq
                                do i = 1, num_fluids
                                    p_K_Star = xi_M*(xi_MP*((pres_L + pi_infs(i)/(1_wp + gammas(i)))* &
                                                            xi_L**(1_wp/gammas(i) + 1_wp) - pi_infs(i)/(1_wp + gammas(i)) - pres_L) + pres_L) + &
                                               xi_P*(xi_PP*((pres_R + pi_infs(i)/(1_wp + gammas(i)))* &
                                                            xi_R**(1_wp/gammas(i) + 1_wp) - pi_infs(i)/(1_wp + gammas(i)) - pres_R) + pres_R)

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
                                    !$acc loop seq
                                    do i = 1, strxe - strxb + 1
                                        flux_rs${XYZ}$_vf(j, k, l, strxb - 1 + i) = &
                                            xi_M*(s_S/(s_L - s_S))*(s_L*rho_L*tau_e_L(i) - rho_L*vel_L(idx1)*tau_e_L(i)) + &
                                            xi_P*(s_S/(s_R - s_S))*(s_R*rho_R*tau_e_R(i) - rho_R*vel_R(idx1)*tau_e_R(i))
                                    end do
                                end if

                                ! REFERENCE MAP FLUX.
                                if (hyperelasticity) then
                                    !$acc loop seq
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
                                        !$acc loop seq
                                        do i = 1, E_idx
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = flux_rs${XYZ}$_vf(j, k, l, i)
                                        end do
                                        !$acc loop seq
                                        do i = intxb, intxe
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = flux_rs${XYZ}$_vf(j, k, l, i)
                                        end do
                                        ! Recalculating the radial momentum geometric source flux
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, momxb - 1 + dir_idx(1)) = &
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, momxb - 1 + dir_idx(1)) - p_Star
                                        ! Geometrical source of the void fraction(s) is zero
                                        !$acc loop seq
                                        do i = advxb, advxe
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = 0_wp
                                        end do
                                    end if
                                #:endif
                                #:if (NORM_DIR == 3)
                                    if (grid_geometry == 3) then
                                        !$acc loop seq
                                        do i = 1, sys_size
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = 0_wp
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
                    !$acc parallel loop collapse(3) gang vector default(present) private(alpha_rho_L, alpha_rho_R, vel_L, vel_R, alpha_L, alpha_R, &
                    !$acc rho_avg, h_avg, gamma_avg, s_L, s_R, s_S, vel_avg_rms, nbub_L, nbub_R, ptilde_L, ptilde_R)
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end

                                !$acc loop seq
                                do i = 1, contxe
                                    alpha_rho_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                    alpha_rho_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)
                                end do

                                !$acc loop seq
                                do i = 1, num_dims
                                    vel_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, contxe + i)
                                    vel_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, contxe + i)
                                end do

                                vel_L_rms = 0._wp; vel_R_rms = 0._wp
                                !$acc loop seq
                                do i = 1, num_dims
                                    vel_L_rms = vel_L_rms + vel_L(i)**2._wp
                                    vel_R_rms = vel_R_rms + vel_R(i)**2._wp
                                end do

                                !$acc loop seq
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
                                !$acc loop seq
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
                                !$acc loop seq
                                do i = 1, num_fluids
                                    rho_R = rho_R + alpha_rho_R(i)
                                    gamma_R = gamma_R + alpha_R(i)*gammas(i)
                                    pi_inf_R = pi_inf_R + alpha_R(i)*pi_infs(i)
                                    qv_R = qv_R + alpha_rho_R(i)*qvs(i)
                                end do

                                E_L = gamma_L*pres_L + pi_inf_L + 5e-1_wp*rho_L*vel_L_rms + qv_L

                                E_R = gamma_R*pres_R + pi_inf_R + 5e-1_wp*rho_R*vel_R_rms + qv_R

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
                                    pres_SL = 5e-1_wp*(pres_L + pres_R + rho_avg*c_avg* &
                                                       (vel_L(dir_idx(1)) - &
                                                        vel_R(dir_idx(1))))

                                    pres_SR = pres_SL

                                    Ms_L = max(1._wp, sqrt(1._wp + ((5e-1_wp + gamma_L)/(1._wp + gamma_L))* &
                                                           (pres_SL/pres_L - 1._wp)*pres_L/ &
                                                           ((pres_L + pi_inf_L/(1._wp + gamma_L)))))
                                    Ms_R = max(1._wp, sqrt(1._wp + ((5e-1_wp + gamma_R)/(1._wp + gamma_R))* &
                                                           (pres_SR/pres_R - 1._wp)*pres_R/ &
                                                           ((pres_R + pi_inf_R/(1._wp + gamma_R)))))

                                    s_L = vel_L(dir_idx(1)) - c_L*Ms_L
                                    s_R = vel_R(dir_idx(1)) + c_R*Ms_R

                                    s_S = 5e-1_wp*((vel_L(dir_idx(1)) + vel_R(dir_idx(1))) + &
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
                                xi_M = (5e-1_wp + sign(5e-1_wp, s_S))
                                xi_P = (5e-1_wp - sign(5e-1_wp, s_S))

                                !$acc loop seq
                                do i = 1, contxe
                                    flux_rs${XYZ}$_vf(j, k, l, i) = &
                                        xi_M*alpha_rho_L(i) &
                                        *(vel_L(dir_idx(1)) + s_M*(xi_L - 1._wp)) &
                                        + xi_P*alpha_rho_R(i) &
                                        *(vel_R(dir_idx(1)) + s_P*(xi_R - 1._wp))
                                end do

                                ! Momentum flux.
                                ! f = \rho u u + p I, q = \rho u, q_star = \xi * \rho*(s_star, v, w)
                                !$acc loop seq
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
                                    !$acc loop seq
                                    do i = 1, num_dims
                                        flux_rs${XYZ}$_vf(j, k, l, contxe + dir_idx(i)) = &
                                            flux_rs${XYZ}$_vf(j, k, l, contxe + dir_idx(i)) + &
                                            xi_M*(dir_flg(dir_idx(i))*(-1._wp*ptilde_L)) &
                                            + xi_P*(dir_flg(dir_idx(i))*(-1._wp*ptilde_R))
                                    end do
                                end if

                                flux_rs${XYZ}$_vf(j, k, l, E_idx) = 0._wp

                                !$acc loop seq
                                do i = alf_idx, alf_idx !only advect the void fraction
                                    flux_rs${XYZ}$_vf(j, k, l, i) = &
                                        xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, i) &
                                        *(vel_L(dir_idx(1)) + s_M*(xi_L - 1._wp)) &
                                        + xi_P*qR_prim_rs${XYZ}$_vf(j + 1, k, l, i) &
                                        *(vel_R(dir_idx(1)) + s_P*(xi_R - 1._wp))
                                end do

                                ! Source for volume fraction advection equation
                                !$acc loop seq
                                do i = 1, num_dims

                                    vel_src_rs${XYZ}$_vf(j, k, l, dir_idx(i)) = 0._wp
                                    !IF ( (model_eqns == 4) .or. (num_fluids==1) ) vel_src_rs_vf(dir_idx(i))%sf(j,k,l) = 0._wp
                                end do

                                flux_src_rs${XYZ}$_vf(j, k, l, advxb) = vel_src_rs${XYZ}$_vf(j, k, l, dir_idx(1))

                                ! Add advection flux for bubble variables
                                if (bubbles_euler) then
                                    !$acc loop seq
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
                                        !$acc loop seq
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
                                        !$acc loop seq
                                        do i = advxb, advxe
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = 0._wp
                                        end do
                                    end if
                                #:endif
                                #:if (NORM_DIR == 3)
                                    if (grid_geometry == 3) then
                                        !$acc loop seq
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
                    !$acc end parallel loop

                elseif (model_eqns == 2 .and. bubbles_euler) then
                    !$acc parallel loop collapse(3) gang vector default(present) private(R0_L, R0_R, V0_L, V0_R, P0_L, P0_R, pbw_L, pbw_R, vel_L, vel_R, &
                    !$acc rho_avg, alpha_L, alpha_R, h_avg, gamma_avg, s_L, s_R, s_S, nbub_L, nbub_R, ptilde_L, ptilde_R, vel_avg_rms, Re_L, Re_R, pcorr, zcoef, vel_L_tmp, vel_R_tmp)
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end

                                !$acc loop seq
                                do i = 1, num_fluids
                                    alpha_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)
                                    alpha_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)
                                end do

                                vel_L_rms = 0._wp; vel_R_rms = 0._wp

                                !$acc loop seq
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
                                    !$acc loop seq
                                    do i = 1, num_fluids
                                        rho_L = rho_L + qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                        gamma_L = gamma_L + qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)*gammas(i)
                                        pi_inf_L = pi_inf_L + qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)*pi_infs(i)
                                        qv_L = qv_L + qL_prim_rs${XYZ}$_vf(j, k, l, i)*qvs(i)
                                    end do
                                else if (num_fluids > 2) then
                                    !$acc loop seq
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
                                    !$acc loop seq
                                    do i = 1, num_fluids
                                        rho_R = rho_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)
                                        gamma_R = gamma_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)*gammas(i)
                                        pi_inf_R = pi_inf_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)*pi_infs(i)
                                        qv_R = qv_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)*qvs(i)
                                    end do
                                else if (num_fluids > 2) then
                                    !$acc loop seq
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
                                        !$acc loop seq
                                        do i = 1, 2
                                            Re_L(i) = dflt_real

                                            if (Re_size(i) > 0) Re_L(i) = 0._wp

                                            !$acc loop seq
                                            do q = 1, Re_size(i)
                                                Re_L(i) = (1._wp - qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + Re_idx(i, q)))/Res(i, q) &
                                                          + Re_L(i)
                                            end do

                                            Re_L(i) = 1._wp/max(Re_L(i), sgm_eps)

                                        end do

                                        !$acc loop seq
                                        do i = 1, 2
                                            Re_R(i) = dflt_real

                                            if (Re_size(i) > 0) Re_R(i) = 0._wp

                                            !$acc loop seq
                                            do q = 1, Re_size(i)
                                                Re_R(i) = (1._wp - qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + Re_idx(i, q)))/Res(i, q) &
                                                          + Re_R(i)
                                            end do

                                            Re_R(i) = 1._wp/max(Re_R(i), sgm_eps)
                                        end do
                                    end if
                                end if

                                E_L = gamma_L*pres_L + pi_inf_L + 5e-1_wp*rho_L*vel_L_rms

                                E_R = gamma_R*pres_R + pi_inf_R + 5e-1_wp*rho_R*vel_R_rms

                                H_L = (E_L + pres_L)/rho_L
                                H_R = (E_R + pres_R)/rho_R

                                if (avg_state == 2) then
                                    !$acc loop seq
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
                                            !$acc loop seq
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

                                    !$acc loop seq
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

                                        !$acc loop seq
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

                                    if ((ptilde_L /= ptilde_L) .or. (ptilde_R /= ptilde_R)) then
                                    end if

                                    rho_avg = 5e-1_wp*(rho_L + rho_R)
                                    H_avg = 5e-1_wp*(H_L + H_R)
                                    gamma_avg = 5e-1_wp*(gamma_L + gamma_R)
                                    vel_avg_rms = 0._wp

                                    !$acc loop seq
                                    do i = 1, num_dims
                                        vel_avg_rms = vel_avg_rms + (5e-1_wp*(vel_L(i) + vel_R(i)))**2._wp
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
                                    !$acc loop seq
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
                                    pres_SL = 5e-1_wp*(pres_L + pres_R + rho_avg*c_avg* &
                                                       (vel_L(dir_idx(1)) - &
                                                        vel_R(dir_idx(1))))

                                    pres_SR = pres_SL

                                    Ms_L = max(1._wp, sqrt(1._wp + ((5e-1_wp + gamma_L)/(1._wp + gamma_L))* &
                                                           (pres_SL/pres_L - 1._wp)*pres_L/ &
                                                           ((pres_L + pi_inf_L/(1._wp + gamma_L)))))
                                    Ms_R = max(1._wp, sqrt(1._wp + ((5e-1_wp + gamma_R)/(1._wp + gamma_R))* &
                                                           (pres_SR/pres_R - 1._wp)*pres_R/ &
                                                           ((pres_R + pi_inf_R/(1._wp + gamma_R)))))

                                    s_L = vel_L(dir_idx(1)) - c_L*Ms_L
                                    s_R = vel_R(dir_idx(1)) + c_R*Ms_R

                                    s_S = 5e-1_wp*((vel_L(dir_idx(1)) + vel_R(dir_idx(1))) + &
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
                                xi_M = (5e-1_wp + sign(5e-1_wp, s_S))
                                xi_P = (5e-1_wp - sign(5e-1_wp, s_S))

                                ! Low Mach correction
                                if (low_Mach == 1) then
                                    @:compute_low_Mach_correction()
                                else
                                    pcorr = 0._wp
                                end if

                                !$acc loop seq
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

                                !$acc loop seq
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
                                !$acc loop seq
                                do i = advxb, advxe
                                    flux_rs${XYZ}$_vf(j, k, l, i) = &
                                        xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, i) &
                                        *(vel_L(dir_idx(1)) + s_M*(xi_L - 1._wp)) &
                                        + xi_P*qR_prim_rs${XYZ}$_vf(j + 1, k, l, i) &
                                        *(vel_R(dir_idx(1)) + s_P*(xi_R - 1._wp))
                                end do

                                ! Source for volume fraction advection equation
                                !$acc loop seq
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
                                !$acc loop seq
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
                                        !$acc loop seq
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
                                        !$acc loop seq
                                        do i = advxb, advxe
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = 0._wp
                                        end do
                                    end if
                                #:endif
                                #:if (NORM_DIR == 3)
                                    if (grid_geometry == 3) then
                                        !$acc loop seq
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
                    !$acc end parallel loop
                else
                    ! 5-EQUATION MODEL WITH HLLC
                    !$acc parallel loop collapse(3) gang vector default(present) private(vel_L, vel_R, Re_L, Re_R, &
                    !$acc rho_avg, h_avg, gamma_avg, alpha_L, alpha_R, alpha_rho_L, alpha_rho_R,                    &
                    !$acc s_L, s_R, s_S, vel_avg_rms, pcorr, zcoef,                                                &
                    !$acc vel_L_tmp, vel_R_tmp, Ys_L, Ys_R, Xs_L, Xs_R, Gamma_iL, Gamma_iR, Cp_iL, Cp_iR,          &
                    !$acc tau_e_L, tau_e_R, xi_field_L, xi_field_R,                                                &
                    !$acc Yi_avg, Phi_avg, h_iL, h_iR, h_avg_2,                                                   &
                    !$acc U_L, U_R, F_L, F_R, F_HLL, U_HLL) copyin(is1,is2,is3)
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end

                                !idx1 = 1; if (dir_idx(1) == 2) idx1 = 2; if (dir_idx(1) == 3) idx1 = 3

                                !$acc loop seq
                                do i = 1, num_fluids
                                    alpha_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)
                                    alpha_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)
                                end do

                                vel_L_rms = 0._wp; vel_R_rms = 0._wp
                                !$acc loop seq
                                do i = 1, num_dims
                                    vel_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, contxe + i)
                                    vel_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, contxe + i)
                                    vel_L_rms = vel_L_rms + vel_L(i)**2._wp
                                    vel_R_rms = vel_R_rms + vel_R(i)**2._wp
                                end do

                                pres_L = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx)
                                pres_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx)

                                if (hypoelasticity) then
                                    !$acc loop seq
                                    do i = 1, strxe - strxb + 1
                                        tau_e_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, strxb - 1 + i)
                                        tau_e_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, strxb - 1 + i)
                                    end do

                                    if (n == 0) then
                                        tau_nn_L = tau_e_L(1)
                                        tau_nn_R = tau_e_R(1)
                                        u_n_L = vel_L(1)
                                        u_n_R = vel_R(1)
                                    else if (p == 0) then
                                        tau_nt_L = tau_e_L(2)
                                        tau_nt_R = tau_e_R(2)
                                        if (idx1 == 1) then
                                            tau_nn_L = tau_e_L(1)
                                            tau_nn_R = tau_e_R(1)
                                            tau_tt_L = tau_e_L(3)
                                            tau_tt_R = tau_e_R(3)
                                            u_n_L = vel_L(1)
                                            u_n_R = vel_R(1)
                                            u_t_L = vel_L(2)
                                            u_t_R = vel_R(2)
                                        else
                                            tau_nn_L = tau_e_L(3)
                                            tau_nn_R = tau_e_R(3)
                                            tau_tt_L = tau_e_L(1)
                                            tau_tt_R = tau_e_R(1)
                                            u_n_L = vel_L(2)
                                            u_n_R = vel_R(2)
                                            u_t_L = vel_L(1)
                                            u_t_R = vel_R(1)
                                        end if
                                    else
                                        ! 3D Cartesian local face basis (n, t1, t2)
                                        if (idx1 == 1) then
                                            u_n_L = vel_L(1); u_n_R = vel_R(1)
                                            u_t_L = vel_L(2); u_t_R = vel_R(2)
                                            u_t2_L = vel_L(3); u_t2_R = vel_R(3)
                                            tau_nn_L = tau_e_L(1); tau_nn_R = tau_e_R(1)
                                            tau_nt_L = tau_e_L(2); tau_nt_R = tau_e_R(2)
                                            tau_nt2_L = tau_e_L(4); tau_nt2_R = tau_e_R(4)
                                            tau_tt_L = tau_e_L(3); tau_tt_R = tau_e_R(3)
                                            tau_t2t2_L = tau_e_L(6); tau_t2t2_R = tau_e_R(6)
                                            tau_t1t2_L = tau_e_L(5); tau_t1t2_R = tau_e_R(5)
                                        else if (idx1 == 2) then
                                            u_n_L = vel_L(2); u_n_R = vel_R(2)
                                            u_t_L = vel_L(1); u_t_R = vel_R(1)
                                            u_t2_L = vel_L(3); u_t2_R = vel_R(3)
                                            tau_nn_L = tau_e_L(3); tau_nn_R = tau_e_R(3)
                                            tau_nt_L = tau_e_L(2); tau_nt_R = tau_e_R(2)
                                            tau_nt2_L = tau_e_L(5); tau_nt2_R = tau_e_R(5)
                                            tau_tt_L = tau_e_L(1); tau_tt_R = tau_e_R(1)
                                            tau_t2t2_L = tau_e_L(6); tau_t2t2_R = tau_e_R(6)
                                            tau_t1t2_L = tau_e_L(4); tau_t1t2_R = tau_e_R(4)
                                        else
                                            u_n_L = vel_L(3); u_n_R = vel_R(3)
                                            u_t_L = vel_L(1); u_t_R = vel_R(1)
                                            u_t2_L = vel_L(2); u_t2_R = vel_R(2)
                                            tau_nn_L = tau_e_L(6); tau_nn_R = tau_e_R(6)
                                            tau_nt_L = tau_e_L(4); tau_nt_R = tau_e_R(4)
                                            tau_nt2_L = tau_e_L(5); tau_nt2_R = tau_e_R(5)
                                            tau_tt_L = tau_e_L(1); tau_tt_R = tau_e_R(1)
                                            tau_t2t2_L = tau_e_L(3); tau_t2t2_R = tau_e_R(3)
                                            tau_t1t2_L = tau_e_L(2); tau_t1t2_R = tau_e_R(2)
                                        end if
                                    end if
                                    pres_tot_L = pres_L - tau_nn_L
                                    pres_tot_R = pres_R - tau_nn_R
                                    if (cyl_coord) then
                                        tau_qq_L = tau_e_L(strxe - strxb + 1)
                                        tau_qq_R = tau_e_R(strxe - strxb + 1)
                                    else
                                        tau_qq_L = 0._wp
                                        tau_qq_R = 0._wp
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

                                ! Change this by splitting it into the cases
                                ! present in the bubbles_euler
                                if (mpp_lim) then
                                    !$acc loop seq
                                    do i = 1, num_fluids
                                        qL_prim_rs${XYZ}$_vf(j, k, l, i) = max(0._wp, qL_prim_rs${XYZ}$_vf(j, k, l, i))
                                        qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i) = min(max(0._wp, qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)), 1._wp)
                                        alpha_L_sum = alpha_L_sum + qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)
                                    end do

                                    !$acc loop seq
                                    do i = 1, num_fluids
                                        qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i) = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)/max(alpha_L_sum, sgm_eps)
                                    end do

                                    !$acc loop seq
                                    do i = 1, num_fluids
                                        qR_prim_rs${XYZ}$_vf(j + 1, k, l, i) = max(0._wp, qR_prim_rs${XYZ}$_vf(j + 1, k, l, i))
                                        qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i) = min(max(0._wp, qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)), 1._wp)
                                        alpha_R_sum = alpha_R_sum + qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)
                                    end do

                                    !$acc loop seq
                                    do i = 1, num_fluids
                                        qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)/max(alpha_R_sum, sgm_eps)
                                    end do
                                end if

                                !$acc loop seq
                                do i = 1, num_fluids
                                    alpha_rho_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                    alpha_rho_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)

                                    rho_L = rho_L + alpha_rho_L(i)
                                    gamma_L = gamma_L + qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)*gammas(i)
                                    pi_inf_L = pi_inf_L + qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)*pi_infs(i)
                                    qv_L = qv_L + alpha_rho_L(i)*qvs(i)

                                    rho_R = rho_R + alpha_rho_R(i)
                                    gamma_R = gamma_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)*gammas(i)
                                    pi_inf_R = pi_inf_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)*pi_infs(i)
                                    qv_R = qv_R + alpha_rho_R(i)*qvs(i)
                                end do

                                if (viscous) then
                                    !$acc loop seq
                                    do i = 1, 2
                                        Re_L(i) = dflt_real

                                        if (Re_size(i) > 0) Re_L(i) = 0._wp

                                        !$acc loop seq
                                        do q = 1, Re_size(i)
                                            Re_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + Re_idx(i, q))/Res(i, q) &
                                                      + Re_L(i)
                                        end do

                                        Re_L(i) = 1._wp/max(Re_L(i), sgm_eps)

                                    end do

                                    !$acc loop seq
                                    do i = 1, 2
                                        Re_R(i) = dflt_real

                                        if (Re_size(i) > 0) Re_R(i) = 0._wp

                                        !$acc loop seq
                                        do q = 1, Re_size(i)
                                            Re_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + Re_idx(i, q))/Res(i, q) &
                                                      + Re_R(i)
                                        end do

                                        Re_R(i) = 1._wp/max(Re_R(i), sgm_eps)
                                    end do
                                end if

                                if (chemistry) then
                                    c_sum_Yi_Phi = 0.0_wp
                                    !$acc loop seq
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

                                    E_L = rho_L*E_L + 5e-1*rho_L*vel_L_rms
                                    E_R = rho_R*E_R + 5e-1*rho_R*vel_R_rms
                                    H_L = (E_L + pres_L)/rho_L
                                    H_R = (E_R + pres_R)/rho_R
                                else
                                    E_L = gamma_L*pres_L + pi_inf_L + 5e-1*rho_L*vel_L_rms + qv_L

                                    E_R = gamma_R*pres_R + pi_inf_R + 5e-1*rho_R*vel_R_rms + qv_R

                                    H_L = (E_L + pres_L)/rho_L
                                    H_R = (E_R + pres_R)/rho_R
                                end if

                                ! ENERGY ADJUSTMENTS FOR HYPOELASTIC ENERGY
                                if (hypoelasticity) then
                                    !$acc loop seq
                                    do i = 1, strxe - strxb + 1
                                        tau_e_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, strxb - 1 + i)
                                        tau_e_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, strxb - 1 + i)
                                    end do
                                    G_L = 0_wp
                                    G_R = 0_wp
                                    !$acc loop seq
                                    do i = 1, num_fluids
                                        G_L = G_L + alpha_L(i)*Gs(i)
                                        G_R = G_R + alpha_R(i)*Gs(i)
                                    end do
                                    !$acc loop seq
                                    do i = 1, strxe - strxb + 1
                                        ! Elastic contribution to energy if G large enough
                                        if ((G_L > verysmall) .and. (G_R > verysmall)) then
                                            E_L = E_L + (tau_e_L(i)*tau_e_L(i))/(4._wp*G_L)
                                            E_R = E_R + (tau_e_R(i)*tau_e_R(i))/(4._wp*G_R)
                                            ! Shear terms doubled: 2D/2D-axisym i==2 only; 3D i==2,4,5
                                            if ((n > 0 .and. p == 0 .and. i == 2) .or. &
                                                (p > 0 .and. (i == 2 .or. i == 4 .or. i == 5))) then
                                                E_L = E_L + (tau_e_L(i)*tau_e_L(i))/(4._wp*G_L)
                                                E_R = E_R + (tau_e_R(i)*tau_e_R(i))/(4._wp*G_R)
                                            end if
                                        end if
                                    end do
                                end if

                                ! ENERGY ADJUSTMENTS FOR HYPERELASTIC ENERGY
                                if (hyperelasticity) then
                                    !$acc loop seq
                                    do i = 1, num_dims
                                        xi_field_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, xibeg - 1 + i)
                                        xi_field_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, xibeg - 1 + i)
                                    end do
                                    G_L = 0_wp
                                    G_R = 0_wp
                                    !$acc loop seq
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
                                    !$acc loop seq
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
                                    !$acc loop seq
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
                                    pres_SL = 5e-1_wp*(pres_L + pres_R + rho_avg*c_avg* &
                                                       (vel_L(idx1) - &
                                                        vel_R(idx1)))

                                    pres_SR = pres_SL

                                    Ms_L = max(1._wp, sqrt(1._wp + ((5e-1_wp + gamma_L)/(1._wp + gamma_L))* &
                                                           (pres_SL/pres_L - 1._wp)*pres_L/ &
                                                           ((pres_L + pi_inf_L/(1._wp + gamma_L)))))
                                    Ms_R = max(1._wp, sqrt(1._wp + ((5e-1_wp + gamma_R)/(1._wp + gamma_R))* &
                                                           (pres_SR/pres_R - 1._wp)*pres_R/ &
                                                           ((pres_R + pi_inf_R/(1._wp + gamma_R)))))

                                    s_L = vel_L(idx1) - c_L*Ms_L
                                    s_R = vel_R(idx1) + c_R*Ms_R

                                    s_S = 5e-1_wp*((vel_L(idx1) + vel_R(idx1)) + &
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
                                xi_M = (5e-1_wp + sign(5e-1_wp, s_S))
                                xi_P = (5e-1_wp - sign(5e-1_wp, s_S))

                                ! Low Mach correction
                                if (low_Mach == 1) then
                                    @:compute_low_Mach_correction()
                                else
                                    pcorr = 0._wp
                                end if

                                ! Hypo HLLC precomputation
                                if (hypoelasticity) then
                                    if (n == 0) then
                                        u_t_L = 0._wp; u_t_R = 0._wp
                                        tau_nt_L = 0._wp; tau_nt_R = 0._wp
                                    end if
                                    if (p == 0) then
                                        u_t2_L = 0._wp; u_t2_R = 0._wp
                                        tau_nt2_L = 0._wp; tau_nt2_R = 0._wp
                                    end if
                                    A_L = rho_L*(s_L - vel_L(idx1))
                                    A_R = rho_R*(s_R - vel_R(idx1))
                                    denom_A = A_R - A_L
                                    u_t_star    = (A_R*u_t_R - A_L*u_t_L + (tau_nt_R - tau_nt_L))/(denom_A + sgm_eps)
                                    tau_nt_star = (A_R*tau_nt_R - A_L*tau_nt_L)/(denom_A + sgm_eps)
                                    u_t2_star    = (A_R*u_t2_R - A_L*u_t2_L + (tau_nt2_R - tau_nt2_L))/(denom_A + sgm_eps)
                                    tau_nt2_star = (A_R*tau_nt2_R - A_L*tau_nt2_L)/(denom_A + sgm_eps)
                                    pres_tot_star = pres_tot_L + A_L*(s_S - vel_L(idx1))
                                end if

                                ! MASS FLUX (unchanged for hypo)
                                !$acc loop seq
                                do i = 1, contxe
                                    flux_rs${XYZ}$_vf(j, k, l, i) = &
                                        xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, i) &
                                        *(vel_L(idx1) + s_M*(xi_L - 1._wp)) &
                                        + xi_P*qR_prim_rs${XYZ}$_vf(j + 1, k, l, i) &
                                        *(vel_R(idx1) + s_P*(xi_R - 1._wp))
                                end do

                                ! MOMENTUM FLUX
                                if (hypoelasticity) then
                                    flux_rs${XYZ}$_vf(j, k, l, contxe + idx1) = &
                                        xi_M*(rho_L*(vel_L(idx1)*vel_L(idx1) + &
                                                     s_M*(xi_L*s_S - vel_L(idx1))) + pres_tot_L) &
                                        + xi_P*(rho_R*(vel_R(idx1)*vel_R(idx1) + &
                                                       s_P*(xi_R*s_S - vel_R(idx1))) + pres_tot_R) &
                                        + (s_M/s_L)*(s_P/s_R)*pcorr
                                    if (n > 0) then
                                        idxi = dir_idx(2)
                                        flux_rs${XYZ}$_vf(j, k, l, contxe + idxi) = &
                                            xi_M*(rho_L*(vel_L(idx1)*u_t_L + &
                                                         s_M*(xi_L*u_t_star - u_t_L)) - tau_nt_L) &
                                            + xi_P*(rho_R*(vel_R(idx1)*u_t_R + &
                                                           s_P*(xi_R*u_t_star - u_t_R)) - tau_nt_R)
                                    end if
                                    if (p > 0) then
                                        idxi = dir_idx(3)
                                        flux_rs${XYZ}$_vf(j, k, l, contxe + idxi) = &
                                            xi_M*(rho_L*(vel_L(idx1)*u_t2_L + &
                                                         s_M*(xi_L*u_t2_star - u_t2_L)) - tau_nt2_L) &
                                            + xi_P*(rho_R*(vel_R(idx1)*u_t2_R + &
                                                           s_P*(xi_R*u_t2_star - u_t2_R)) - tau_nt2_R)
                                    end if
                                else
                                    !$acc loop seq
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
                                end if

                                ! ENERGY FLUX
                                if (hypoelasticity) then
                                    flux_rs${XYZ}$_vf(j, k, l, E_idx) = &
                                        xi_M*((E_L + pres_tot_L)*vel_L(idx1) &
                                              - u_t_L*tau_nt_L - u_t2_L*tau_nt2_L + &
                                              s_M*(xi_L*(E_L + (s_S - vel_L(idx1))* &
                                                         (rho_L*s_S + pres_tot_L/(s_L - vel_L(idx1))) &
                                                         + (u_t_L*tau_nt_L - u_t_star*tau_nt_star) &
                                                           /(s_L - vel_L(idx1)) &
                                                         + (u_t2_L*tau_nt2_L - u_t2_star*tau_nt2_star) &
                                                           /(s_L - vel_L(idx1))) &
                                                   - E_L)) &
                                        + xi_P*((E_R + pres_tot_R)*vel_R(idx1) &
                                                - u_t_R*tau_nt_R - u_t2_R*tau_nt2_R + &
                                                s_P*(xi_R*(E_R + (s_S - vel_R(idx1))* &
                                                           (rho_R*s_S + pres_tot_R/(s_R - vel_R(idx1))) &
                                                           + (u_t_R*tau_nt_R - u_t_star*tau_nt_star) &
                                                             /(s_R - vel_R(idx1)) &
                                                           + (u_t2_R*tau_nt2_R - u_t2_star*tau_nt2_star) &
                                                             /(s_R - vel_R(idx1))) &
                                                     - E_R)) &
                                        + (s_M/s_L)*(s_P/s_R)*pcorr*s_S
                                else
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
                                end if

                                ! ELASTICITY (non-hypo elastic stress corrections)
                                if (elasticity .and. .not. hypoelasticity) then
                                    flux_ene_e = 0_wp
                                    !$acc loop seq
                                    do i = 1, num_dims
                                        idxi = dir_idx(i)
                                        flux_rs${XYZ}$_vf(j, k, l, contxe + idxi) = &
                                            flux_rs${XYZ}$_vf(j, k, l, contxe + idxi) &
                                            - xi_M*tau_e_L(dir_idx_tau(i)) - xi_P*tau_e_R(dir_idx_tau(i))
                                        flux_ene_e = flux_ene_e - &
                                                     xi_M*(vel_L(idxi)*tau_e_L(dir_idx_tau(i)) + &
                                                           s_M*(xi_L*((s_S - vel_L(i))*(tau_e_L(dir_idx_tau(i))/(s_L - vel_L(i)))))) - &
                                                     xi_P*(vel_R(idxi)*tau_e_R(dir_idx_tau(i)) + &
                                                           s_P*(xi_R*((s_S - vel_R(i))*(tau_e_R(dir_idx_tau(i))/(s_R - vel_R(i))))))
                                    end do
                                    flux_rs${XYZ}$_vf(j, k, l, E_idx) = flux_rs${XYZ}$_vf(j, k, l, E_idx) + flux_ene_e
                                end if

                                ! HYPOELASTIC STRESS FLUX + INTERFACE VELOCITY EXPORT
                                if (hypoelasticity) then
                                    if (n == 0) then
                                        flux_rs${XYZ}$_vf(j, k, l, strxb) = &
                                            xi_M*rho_L*tau_nn_L*(vel_L(idx1) + s_M*(xi_L - 1._wp)) &
                                            + xi_P*rho_R*tau_nn_R*(vel_R(idx1) + s_P*(xi_R - 1._wp))
                                    else if (p == 0) then
                                        if (idx1 == 1) then
                                            flux_rs${XYZ}$_vf(j, k, l, strxb) = &
                                                xi_M*rho_L*tau_nn_L*(vel_L(idx1) + s_M*(xi_L - 1._wp)) &
                                                + xi_P*rho_R*tau_nn_R*(vel_R(idx1) + s_P*(xi_R - 1._wp))
                                            flux_rs${XYZ}$_vf(j, k, l, strxb + 1) = &
                                                xi_M*(rho_L*vel_L(idx1)*tau_nt_L + s_M*(rho_L*xi_L*tau_nt_star - rho_L*tau_nt_L)) &
                                                + xi_P*(rho_R*vel_R(idx1)*tau_nt_R + s_P*(rho_R*xi_R*tau_nt_star - rho_R*tau_nt_R))
                                            flux_rs${XYZ}$_vf(j, k, l, strxb + 2) = &
                                                xi_M*rho_L*tau_tt_L*(vel_L(idx1) + s_M*(xi_L - 1._wp)) &
                                                + xi_P*rho_R*tau_tt_R*(vel_R(idx1) + s_P*(xi_R - 1._wp))
                                        else
                                            flux_rs${XYZ}$_vf(j, k, l, strxb + 2) = &
                                                xi_M*rho_L*tau_nn_L*(vel_L(idx1) + s_M*(xi_L - 1._wp)) &
                                                + xi_P*rho_R*tau_nn_R*(vel_R(idx1) + s_P*(xi_R - 1._wp))
                                            flux_rs${XYZ}$_vf(j, k, l, strxb + 1) = &
                                                xi_M*(rho_L*vel_L(idx1)*tau_nt_L + s_M*(rho_L*xi_L*tau_nt_star - rho_L*tau_nt_L)) &
                                                + xi_P*(rho_R*vel_R(idx1)*tau_nt_R + s_P*(rho_R*xi_R*tau_nt_star - rho_R*tau_nt_R))
                                            flux_rs${XYZ}$_vf(j, k, l, strxb) = &
                                                xi_M*rho_L*tau_tt_L*(vel_L(idx1) + s_M*(xi_L - 1._wp)) &
                                                + xi_P*rho_R*tau_tt_R*(vel_R(idx1) + s_P*(xi_R - 1._wp))
                                        end if
                                    else
                                        ! 3D Cartesian: 6 stress components mapped via face basis
                                        if (idx1 == 1) then
                                            itnn = strxb; itnt1 = strxb + 1; itnt2 = strxb + 3
                                            itss1 = strxb + 2; itss2 = strxb + 5; itt12 = strxb + 4
                                        else if (idx1 == 2) then
                                            itnn = strxb + 2; itnt1 = strxb + 1; itnt2 = strxb + 4
                                            itss1 = strxb; itss2 = strxb + 5; itt12 = strxb + 3
                                        else
                                            itnn = strxb + 5; itnt1 = strxb + 3; itnt2 = strxb + 4
                                            itss1 = strxb; itss2 = strxb + 2; itt12 = strxb + 1
                                        end if
                                        flux_rs${XYZ}$_vf(j, k, l, itnn) = &
                                            xi_M*rho_L*tau_nn_L*(vel_L(idx1) + s_M*(xi_L - 1._wp)) &
                                            + xi_P*rho_R*tau_nn_R*(vel_R(idx1) + s_P*(xi_R - 1._wp))
                                        flux_rs${XYZ}$_vf(j, k, l, itnt1) = &
                                            xi_M*(rho_L*vel_L(idx1)*tau_nt_L + s_M*(rho_L*xi_L*tau_nt_star - rho_L*tau_nt_L)) &
                                            + xi_P*(rho_R*vel_R(idx1)*tau_nt_R + s_P*(rho_R*xi_R*tau_nt_star - rho_R*tau_nt_R))
                                        flux_rs${XYZ}$_vf(j, k, l, itnt2) = &
                                            xi_M*(rho_L*vel_L(idx1)*tau_nt2_L + s_M*(rho_L*xi_L*tau_nt2_star - rho_L*tau_nt2_L)) &
                                            + xi_P*(rho_R*vel_R(idx1)*tau_nt2_R + s_P*(rho_R*xi_R*tau_nt2_star - rho_R*tau_nt2_R))
                                        flux_rs${XYZ}$_vf(j, k, l, itss1) = &
                                            xi_M*rho_L*tau_tt_L*(vel_L(idx1) + s_M*(xi_L - 1._wp)) &
                                            + xi_P*rho_R*tau_tt_R*(vel_R(idx1) + s_P*(xi_R - 1._wp))
                                        flux_rs${XYZ}$_vf(j, k, l, itss2) = &
                                            xi_M*rho_L*tau_t2t2_L*(vel_L(idx1) + s_M*(xi_L - 1._wp)) &
                                            + xi_P*rho_R*tau_t2t2_R*(vel_R(idx1) + s_P*(xi_R - 1._wp))
                                        flux_rs${XYZ}$_vf(j, k, l, itt12) = &
                                            xi_M*rho_L*tau_t1t2_L*(vel_L(idx1) + s_M*(xi_L - 1._wp)) &
                                            + xi_P*rho_R*tau_t1t2_R*(vel_R(idx1) + s_P*(xi_R - 1._wp))
                                    end if
                                    if (cyl_coord) then
                                        flux_rs${XYZ}$_vf(j, k, l, strxe) = &
                                            xi_M*rho_L*tau_qq_L*(vel_L(idx1) + s_M*(xi_L - 1._wp)) &
                                            + xi_P*rho_R*tau_qq_R*(vel_R(idx1) + s_P*(xi_R - 1._wp))
                                    end if

                                    ! Exported interface velocities (explicit fan selection)
                                    if (s_L >= 0._wp) then
                                        u_n_HLLC = vel_L(idx1); u_t_HLLC = u_t_L; u_t2_HLLC = u_t2_L
                                    elseif (s_R <= 0._wp) then
                                        u_n_HLLC = vel_R(idx1); u_t_HLLC = u_t_R; u_t2_HLLC = u_t2_R
                                    else
                                        u_n_HLLC = s_S*(xi_M*xi_L + xi_P*xi_R); u_t_HLLC = u_t_star; u_t2_HLLC = u_t2_star
                                    end if
                                    if (p == 0) then
                                        if (idx1 == 1) then
                                            nc_iface_vel_rs${XYZ}$_vf(j, k, l, 1) = u_n_HLLC
                                            nc_iface_vel_rs${XYZ}$_vf(j, k, l, 2) = u_t_HLLC
                                        else
                                            nc_iface_vel_rs${XYZ}$_vf(j, k, l, 1) = u_t_HLLC
                                            nc_iface_vel_rs${XYZ}$_vf(j, k, l, 2) = u_n_HLLC
                                        end if
                                    else
                                        if (idx1 == 1) then
                                            nc_iface_vel_rs${XYZ}$_vf(j, k, l, 1) = u_n_HLLC
                                            nc_iface_vel_rs${XYZ}$_vf(j, k, l, 2) = u_t_HLLC
                                            nc_iface_vel_rs${XYZ}$_vf(j, k, l, 3) = u_t2_HLLC
                                        else if (idx1 == 2) then
                                            nc_iface_vel_rs${XYZ}$_vf(j, k, l, 1) = u_t_HLLC
                                            nc_iface_vel_rs${XYZ}$_vf(j, k, l, 2) = u_n_HLLC
                                            nc_iface_vel_rs${XYZ}$_vf(j, k, l, 3) = u_t2_HLLC
                                        else
                                            nc_iface_vel_rs${XYZ}$_vf(j, k, l, 1) = u_t_HLLC
                                            nc_iface_vel_rs${XYZ}$_vf(j, k, l, 2) = u_t2_HLLC
                                            nc_iface_vel_rs${XYZ}$_vf(j, k, l, 3) = u_n_HLLC
                                        end if
                                    end if
                                end if

                                ! VOLUME FRACTION FLUX (unchanged for hypo)
                                !$acc loop seq
                                do i = advxb, advxe
                                    flux_rs${XYZ}$_vf(j, k, l, i) = &
                                        xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, i) &
                                        *(vel_L(idx1) + s_M*(xi_L - 1._wp)) &
                                        + xi_P*qR_prim_rs${XYZ}$_vf(j + 1, k, l, i) &
                                        *(vel_R(idx1) + s_P*(xi_R - 1._wp))
                                end do

                                ! VOLUME FRACTION SOURCE FLUX
                                !$acc loop seq
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

                                ! REFERENCE MAP FLUX
                                if (hyperelasticity) then
                                    !$acc loop seq
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
                                    !$acc loop seq
                                    do i = chemxb, chemxe
                                        Y_L = qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                        Y_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)

                                        flux_rs${XYZ}$_vf(j, k, l, i) = xi_M*rho_L*Y_L*(vel_L(idx1) + s_M*(xi_L - 1._wp)) &
                                                                        + xi_P*rho_R*Y_R*(vel_R(idx1) + s_P*(xi_R - 1._wp))
                                        flux_src_rs${XYZ}$_vf(j, k, l, i) = 0.0_wp
                                    end do
                                end if

                                ! ===== HLLC-ADC blending for hypoelasticity =====
                                if (riemann_hypo_ADC .and. hypoelasticity) then

                                    ! --- Build U_L, U_R and F_L, F_R in local-basis layout ---
                                    !$acc loop seq
                                    do i = 1, num_fluids
                                        U_L(i) = alpha_rho_L(i)
                                        U_R(i) = alpha_rho_R(i)
                                        U_L(advxb - 1 + i) = alpha_L(i)
                                        U_R(advxb - 1 + i) = alpha_R(i)
                                        F_L(i) = alpha_rho_L(i)*u_n_L
                                        F_R(i) = alpha_rho_R(i)*u_n_R
                                        F_L(advxb - 1 + i) = alpha_L(i)*u_n_L
                                        F_R(advxb - 1 + i) = alpha_R(i)*u_n_R
                                    end do

                                    U_L(momxb) = rho_L*u_n_L
                                    U_R(momxb) = rho_R*u_n_R
                                    F_L(momxb) = rho_L*u_n_L*u_n_L + pres_tot_L
                                    F_R(momxb) = rho_R*u_n_R*u_n_R + pres_tot_R
                                    if (n > 0) then
                                        U_L(momxb + 1) = rho_L*u_t_L
                                        U_R(momxb + 1) = rho_R*u_t_R
                                        F_L(momxb + 1) = rho_L*u_n_L*u_t_L - tau_nt_L
                                        F_R(momxb + 1) = rho_R*u_n_R*u_t_R - tau_nt_R
                                    end if
                                    if (p > 0) then
                                        U_L(momxb + 2) = rho_L*u_t2_L
                                        U_R(momxb + 2) = rho_R*u_t2_R
                                        F_L(momxb + 2) = rho_L*u_n_L*u_t2_L - tau_nt2_L
                                        F_R(momxb + 2) = rho_R*u_n_R*u_t2_R - tau_nt2_R
                                    end if

                                    U_L(E_idx) = E_L
                                    U_R(E_idx) = E_R
                                    F_L(E_idx) = (E_L + pres_tot_L)*u_n_L - u_t_L*tau_nt_L - u_t2_L*tau_nt2_L
                                    F_R(E_idx) = (E_R + pres_tot_R)*u_n_R - u_t_R*tau_nt_R - u_t2_R*tau_nt2_R

                                    U_L(strxb) = rho_L*tau_nn_L
                                    U_R(strxb) = rho_R*tau_nn_R
                                    F_L(strxb) = rho_L*u_n_L*tau_nn_L
                                    F_R(strxb) = rho_R*u_n_R*tau_nn_R
                                    if (n > 0) then
                                        U_L(strxb + 1) = rho_L*tau_nt_L
                                        U_R(strxb + 1) = rho_R*tau_nt_R
                                        F_L(strxb + 1) = rho_L*u_n_L*tau_nt_L
                                        F_R(strxb + 1) = rho_R*u_n_R*tau_nt_R
                                        U_L(strxb + 2) = rho_L*tau_tt_L
                                        U_R(strxb + 2) = rho_R*tau_tt_R
                                        F_L(strxb + 2) = rho_L*u_n_L*tau_tt_L
                                        F_R(strxb + 2) = rho_R*u_n_R*tau_tt_R
                                    end if
                                    if (p > 0) then
                                        U_L(strxb + 3) = rho_L*tau_nt2_L
                                        U_R(strxb + 3) = rho_R*tau_nt2_R
                                        F_L(strxb + 3) = rho_L*u_n_L*tau_nt2_L
                                        F_R(strxb + 3) = rho_R*u_n_R*tau_nt2_R
                                        U_L(strxb + 4) = rho_L*tau_t1t2_L
                                        U_R(strxb + 4) = rho_R*tau_t1t2_R
                                        F_L(strxb + 4) = rho_L*u_n_L*tau_t1t2_L
                                        F_R(strxb + 4) = rho_R*u_n_R*tau_t1t2_R
                                        U_L(strxb + 5) = rho_L*tau_t2t2_L
                                        U_R(strxb + 5) = rho_R*tau_t2t2_R
                                        F_L(strxb + 5) = rho_L*u_n_L*tau_t2t2_L
                                        F_R(strxb + 5) = rho_R*u_n_R*tau_t2t2_R
                                    end if
                                    if (cyl_coord) then
                                        U_L(strxe) = rho_L*tau_qq_L
                                        U_R(strxe) = rho_R*tau_qq_R
                                        F_L(strxe) = rho_L*u_n_L*tau_qq_L
                                        F_R(strxe) = rho_R*u_n_R*tau_qq_R
                                    end if

                                    ! --- Compute F_HLL and U_HLL ---
                                    if (s_L >= 0._wp) then
                                        !$acc loop seq
                                        do i = 1, sys_size
                                            F_HLL(i) = F_L(i)
                                            U_HLL(i) = U_L(i)
                                        end do
                                        rho_HLL = rho_L
                                        u_n_HLL = u_n_L; u_t_HLL = u_t_L; u_t2_HLL = u_t2_L
                                    elseif (s_R <= 0._wp) then
                                        !$acc loop seq
                                        do i = 1, sys_size
                                            F_HLL(i) = F_R(i)
                                            U_HLL(i) = U_R(i)
                                        end do
                                        rho_HLL = rho_R
                                        u_n_HLL = u_n_R; u_t_HLL = u_t_R; u_t2_HLL = u_t2_R
                                    else
                                        !$acc loop seq
                                        do i = 1, sys_size
                                            F_HLL(i) = (s_R*F_L(i) - s_L*F_R(i) + s_L*s_R*(U_R(i) - U_L(i))) &
                                                       /(s_R - s_L + verysmall)
                                            U_HLL(i) = (s_R*U_R(i) - s_L*U_L(i) - (F_R(i) - F_L(i))) &
                                                       /(s_R - s_L + verysmall)
                                        end do
                                        u_n_HLL = (s_R*u_n_L - s_L*u_n_R)/(s_R - s_L + verysmall)
                                        u_t_HLL = 0._wp; u_t2_HLL = 0._wp
                                        if (n > 0) u_t_HLL = (s_R*u_t_L - s_L*u_t_R)/(s_R - s_L + verysmall)
                                        if (p > 0) u_t2_HLL = (s_R*u_t2_L - s_L*u_t2_R)/(s_R - s_L + verysmall)
                                    end if

                                    ! --- ADC sensor ---
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

                                    ! --- Blend selected flux components ---
                                    ! Partial densities: blend (same indices in local and physical basis)
                                    !$acc loop seq
                                    do i = 1, contxe
                                        flux_rs${XYZ}$_vf(j, k, l, i) = &
                                            F_HLL(i) + phi*(flux_rs${XYZ}$_vf(j, k, l, i) - F_HLL(i))
                                    end do

                                    ! Energy: blend
                                    flux_rs${XYZ}$_vf(j, k, l, E_idx) = &
                                        F_HLL(E_idx) + phi*(flux_rs${XYZ}$_vf(j, k, l, E_idx) - F_HLL(E_idx))

                                    ! Volume fractions: blend
                                    !$acc loop seq
                                    do i = advxb, advxe
                                        flux_rs${XYZ}$_vf(j, k, l, i) = &
                                            F_HLL(i) + phi*(flux_rs${XYZ}$_vf(j, k, l, i) - F_HLL(i))
                                    end do

                                    ! Stresses: blend tau_nn, tau_t1t1, tau_t2t2, tau_t1t2, tau_qq
                                    ! Keep HLLC (no blend): tau_nt1, tau_nt2
                                    if (n == 0) then
                                        ! 1D: only tau_nn (strxb in both local and physical basis)
                                        flux_rs${XYZ}$_vf(j, k, l, strxb) = &
                                            F_HLL(strxb) + phi*(flux_rs${XYZ}$_vf(j, k, l, strxb) - F_HLL(strxb))
                                    else if (p == 0) then
                                        ! 2D: tau_nn and tau_t1t1 (mapped by sweep direction)
                                        if (idx1 == 1) then
                                            flux_rs${XYZ}$_vf(j, k, l, strxb) = &
                                                F_HLL(strxb) + phi*(flux_rs${XYZ}$_vf(j, k, l, strxb) - F_HLL(strxb))
                                            flux_rs${XYZ}$_vf(j, k, l, strxb + 2) = &
                                                F_HLL(strxb + 2) + phi*(flux_rs${XYZ}$_vf(j, k, l, strxb + 2) - F_HLL(strxb + 2))
                                        else
                                            flux_rs${XYZ}$_vf(j, k, l, strxb + 2) = &
                                                F_HLL(strxb) + phi*(flux_rs${XYZ}$_vf(j, k, l, strxb + 2) - F_HLL(strxb))
                                            flux_rs${XYZ}$_vf(j, k, l, strxb) = &
                                                F_HLL(strxb + 2) + phi*(flux_rs${XYZ}$_vf(j, k, l, strxb) - F_HLL(strxb + 2))
                                        end if
                                    else
                                        ! 3D: use existing itnn/itss1/itss2/itt12 mapping
                                        flux_rs${XYZ}$_vf(j, k, l, itnn) = &
                                            F_HLL(strxb) + phi*(flux_rs${XYZ}$_vf(j, k, l, itnn) - F_HLL(strxb))
                                        flux_rs${XYZ}$_vf(j, k, l, itss1) = &
                                            F_HLL(strxb + 2) + phi*(flux_rs${XYZ}$_vf(j, k, l, itss1) - F_HLL(strxb + 2))
                                        flux_rs${XYZ}$_vf(j, k, l, itss2) = &
                                            F_HLL(strxb + 5) + phi*(flux_rs${XYZ}$_vf(j, k, l, itss2) - F_HLL(strxb + 5))
                                        flux_rs${XYZ}$_vf(j, k, l, itt12) = &
                                            F_HLL(strxb + 4) + phi*(flux_rs${XYZ}$_vf(j, k, l, itt12) - F_HLL(strxb + 4))
                                    end if
                                    if (cyl_coord) then
                                        flux_rs${XYZ}$_vf(j, k, l, strxe) = &
                                            F_HLL(strxe) + phi*(flux_rs${XYZ}$_vf(j, k, l, strxe) - F_HLL(strxe))
                                    end if

                                    ! --- Blend interface velocities ---
                                    u_n_HLLC = u_n_HLL + phi*(u_n_HLLC - u_n_HLL)
                                    u_t_HLLC = u_t_HLL + phi*(u_t_HLLC - u_t_HLL)
                                    u_t2_HLLC = u_t2_HLL + phi*(u_t2_HLLC - u_t2_HLL)

                                    ! Overwrite vel_src with blended velocities
                                    vel_src_rs${XYZ}$_vf(j, k, l, idx1) = u_n_HLLC
                                    if (n > 0) vel_src_rs${XYZ}$_vf(j, k, l, dir_idx(2)) = u_t_HLLC
                                    if (p > 0) vel_src_rs${XYZ}$_vf(j, k, l, dir_idx(3)) = u_t2_HLLC

                                    ! Update advection source flux with blended normal velocity
                                    flux_src_rs${XYZ}$_vf(j, k, l, advxb) = u_n_HLLC

                                    ! Overwrite nc_iface_vel with blended velocities
                                    if (p == 0) then
                                        if (idx1 == 1) then
                                            nc_iface_vel_rs${XYZ}$_vf(j, k, l, 1) = u_n_HLLC
                                            nc_iface_vel_rs${XYZ}$_vf(j, k, l, 2) = u_t_HLLC
                                        else
                                            nc_iface_vel_rs${XYZ}$_vf(j, k, l, 1) = u_t_HLLC
                                            nc_iface_vel_rs${XYZ}$_vf(j, k, l, 2) = u_n_HLLC
                                        end if
                                    else
                                        if (idx1 == 1) then
                                            nc_iface_vel_rs${XYZ}$_vf(j, k, l, 1) = u_n_HLLC
                                            nc_iface_vel_rs${XYZ}$_vf(j, k, l, 2) = u_t_HLLC
                                            nc_iface_vel_rs${XYZ}$_vf(j, k, l, 3) = u_t2_HLLC
                                        else if (idx1 == 2) then
                                            nc_iface_vel_rs${XYZ}$_vf(j, k, l, 1) = u_t_HLLC
                                            nc_iface_vel_rs${XYZ}$_vf(j, k, l, 2) = u_n_HLLC
                                            nc_iface_vel_rs${XYZ}$_vf(j, k, l, 3) = u_t2_HLLC
                                        else
                                            nc_iface_vel_rs${XYZ}$_vf(j, k, l, 1) = u_t_HLLC
                                            nc_iface_vel_rs${XYZ}$_vf(j, k, l, 2) = u_t2_HLLC
                                            nc_iface_vel_rs${XYZ}$_vf(j, k, l, 3) = u_n_HLLC
                                        end if
                                    end if

                                end if
                                ! ===== END HLLC-ADC =====

                                ! ! ===== Direct HLLC hypo star-state path =====
                                ! if (.false. .and. hypoelasticity) then
                                !     eps = 1e-12

                                !     ! Prim to Cons: Find U_L & U_R
                                !     ! Ordered: alpha1*rho1, (alpha2*rho2), rho, vel_n, vel_t, E, alpha1, (alpha2), rho*tau_nn, rho*tau_nt, rho*tau_tt
                                !     do i = 1, num_fluids
                                !         U_L(i) = alpha_rho_L(i)
                                !         U_R(i) = alpha_rho_R(i)

                                !         U_L(advxb - 1 + i) = alpha_L(i)
                                !         U_R(advxb - 1 + i) = alpha_R(i)
                                !     end do

                                !     U_L(momxb) = rho_L*u_n_L
                                !     U_R(momxb) = rho_R*u_n_R
                                !     if (n /= 0) then
                                !         U_L(momxb + 1) = rho_L*u_t_L
                                !         U_R(momxb + 1) = rho_R*u_t_R
                                !     end if

                                !     U_L(E_idx) = E_L
                                !     U_R(E_idx) = E_R

                                !     U_L(strxb) = rho_L*tau_nn_L
                                !     U_R(strxb) = rho_R*tau_nn_R
                                !     if (n /= 0) then
                                !         U_L(strxb + 1) = rho_L*tau_nt_L
                                !         U_L(strxb + 2) = rho_L*tau_tt_L
                                !         U_R(strxb + 1) = rho_R*tau_nt_R
                                !         U_R(strxb + 2) = rho_R*tau_tt_R
                                !     end if
                                !     if (cyl_coord) then
                                !         U_L(strxe) = rho_L*tau_qq_L
                                !         U_R(strxe) = rho_R*tau_qq_R
                                !     end if

                                !     ! Zero tangential quantities for 1D before star-state computation
                                !     if (n == 0) then
                                !         u_t_L = 0._wp; u_t_R = 0._wp
                                !         tau_nt_L = 0._wp; tau_nt_R = 0._wp
                                !     end if

                                !     ! Find prim star_L and star_R
                                !     A_L = rho_L*(S_L - u_n_L)
                                !     A_R = rho_R*(S_R - u_n_R)
                                !     denom_A = (A_R - A_L)

                                !     S_Mid = ((pres_tot_R - pres_tot_L) + A_L*u_n_L - A_R*u_n_R) / (A_L - A_R + eps) ! fallback not used

                                !     pres_tot_star = pres_tot_L + A_L*(S_Mid - u_n_L)

                                !     u_t_star    = (A_R*u_t_R - A_L*u_t_L + (tau_nt_R - tau_nt_L)) / (denom_A + eps)
                                !     tau_nt_star = (A_R*tau_nt_R - A_L*tau_nt_L) / (denom_A + eps)

                                !     denom_L = (S_L - S_Mid) + eps
                                !     denom_R = (S_R - S_Mid) + eps
                                !     rho_L_star = rho_L * (S_L - u_n_L) / denom_L
                                !     rho_R_star = rho_R * (S_R - u_n_R) / denom_R
                                !     fac_L = (S_L - u_n_L) / denom_L
                                !     fac_R = (S_R - u_n_R) / denom_R

                                !     E_L_star = ((S_L - u_n_L)*E_L - u_n_L*pres_tot_L + S_Mid*pres_tot_star - u_t_star*tau_nt_star + u_t_L*tau_nt_L) / denom_L
                                !     E_R_star = ((S_R - u_n_R)*E_R - u_n_R*pres_tot_R + S_Mid*pres_tot_star - u_t_star*tau_nt_star + u_t_R*tau_nt_R) / denom_R

                                !     ! Find U_star_L & U_star_R
                                !     do i = 1, num_fluids
                                !         U_star_L(i) = alpha_rho_L(i)*fac_L
                                !         U_star_R(i) = alpha_rho_R(i)*fac_R

                                !         U_star_L(advxb - 1 + i) = alpha_L(i)*fac_L
                                !         U_star_R(advxb - 1 + i) = alpha_R(i)*fac_R
                                !     end do

                                !     U_star_L(momxb) = rho_L_star*S_Mid
                                !     U_star_R(momxb) = rho_R_star*S_Mid
                                !     if (n /= 0) then
                                !         U_star_L(momxb + 1) = rho_L_star*u_t_star
                                !         U_star_R(momxb + 1) = rho_R_star*u_t_star
                                !     end if

                                !     U_star_L(E_idx) = E_L_star
                                !     U_star_R(E_idx) = E_R_star

                                !     U_star_L(strxb) = rho_L_star*tau_nn_L
                                !     U_star_R(strxb) = rho_R_star*tau_nn_R
                                !     if (n /= 0) then
                                !         U_star_L(strxb + 1) = rho_L_star*tau_nt_star
                                !         U_star_L(strxb + 2) = rho_L_star*tau_tt_L
                                !         U_star_R(strxb + 1) = rho_R_star*tau_nt_star
                                !         U_star_R(strxb + 2) = rho_R_star*tau_tt_R
                                !     end if
                                !     if (cyl_coord) then
                                !         U_star_L(strxe) = rho_L_star*tau_qq_L
                                !         U_star_R(strxe) = rho_R_star*tau_qq_R
                                !     end if

                                !     ! Find F_L & F_R
                                !     do i = 1, num_fluids
                                !         F_L(i) = alpha_rho_L(i)*u_n_L
                                !         F_R(i) = alpha_rho_R(i)*u_n_R

                                !         F_L(advxb - 1 + i) = alpha_L(i)*u_n_L
                                !         F_R(advxb - 1 + i) = alpha_R(i)*u_n_R
                                !     end do

                                !     F_L(momxb) = rho_L*u_n_L*u_n_L + pres_tot_L
                                !     F_R(momxb) = rho_R*u_n_R*u_n_R + pres_tot_R
                                !     if (n /= 0) then
                                !         F_L(momxb + 1) = rho_L*u_n_L*u_t_L - tau_nt_L
                                !         F_R(momxb + 1) = rho_R*u_n_R*u_t_R - tau_nt_R
                                !     end if

                                !     F_L(E_idx) = (E_L + pres_tot_L)*u_n_L - u_t_L*tau_nt_L
                                !     F_R(E_idx) = (E_R + pres_tot_R)*u_n_R - u_t_R*tau_nt_R

                                !     F_L(strxb) = rho_L*u_n_L*tau_nn_L
                                !     F_R(strxb) = rho_R*u_n_R*tau_nn_R
                                !     if (n /= 0) then
                                !         F_L(strxb + 1) = rho_L*u_n_L*tau_nt_L
                                !         F_L(strxb + 2) = rho_L*u_n_L*tau_tt_L
                                !         F_R(strxb + 1) = rho_R*u_n_R*tau_nt_R
                                !         F_R(strxb + 2) = rho_R*u_n_R*tau_tt_R
                                !     end if
                                !     if (cyl_coord) then
                                !         F_L(strxe) = rho_L*u_n_L*tau_qq_L
                                !         F_R(strxe) = rho_R*u_n_R*tau_qq_R
                                !     end if

                                !     ! Find F_star_L & F_star_R (array operation)
                                !     F_star_L = F_L + S_L * (U_star_L - U_L)
                                !     F_star_R = F_R + S_R * (U_star_R - U_R)

                                !     ! Find F_HLLC and velocity fluxes based on wave location
                                !     if (S_L >= 0d0) then
                                !         F_HLLC = F_L
                                !         u_n_HLLC = u_n_L
                                !         u_t_HLLC = u_t_L
                                !     else if (S_R <= 0d0) then
                                !         F_HLLC = F_R
                                !         u_n_HLLC = u_n_R
                                !         u_t_HLLC = u_t_R
                                !     else if (S_mid >= 0d0) then
                                !         F_HLLC = F_star_L
                                !         u_n_HLLC = S_Mid
                                !         u_t_HLLC = u_t_star
                                !     else
                                !         F_HLLC = F_star_R
                                !         u_n_HLLC = S_Mid
                                !         u_t_HLLC = u_t_star
                                !     end if

                                !     ! === ADC BLENDING ===

                                !     if (riemann_hypo_ADC) then

                                !         ! Find F_HLL and velocity fluxes for ADC
                                !         F_HLL = F_HLLC
                                !         u_n_HLL = u_n_HLLC
                                !         u_t_HLL = u_t_HLLC
                                !         if (S_L < 0d0 .and. S_R > 0d0) then
                                !             F_HLL   = (S_R*F_L - S_L*F_R + S_L*S_R*(U_R-U_L)) / (S_R-S_L + eps)
                                !             u_n_HLL = (S_R*u_n_L - S_L*u_n_R) / (S_R-S_L + eps)
                                !             u_t_HLL = (S_R*u_t_L - S_L*u_t_R) / (S_R-S_L + eps)
                                !         end if

                                !         ! Find phi
                                !         ! Total normal stress Σ = p - tau_nn on each side
                                !         Sigma_L   = pres_tot_L
                                !         Sigma_R   = pres_tot_R
                                !         dSigma    = Sigma_R - Sigma_L
                                !         Sigma_ref = max( max(abs(Sigma_L), abs(Sigma_R)), eps )

                                !         ! Directional fast speeds (normal), matching Python _a_normal:
                                !         ! a^2 = c^2 + ( (4/3)G + tau_nn ) / ρ
                                !         a_L_ref = sqrt( max( eps, c_L*c_L + (((4._wp/3._wp)*G_L + tau_nn_L)/rho_L) ) )
                                !         a_R_ref = sqrt( max( eps, c_R*c_R + (((4._wp/3._wp)*G_R + tau_nn_R)/rho_R) ) )
                                !         a_ref   = max( max(a_L_ref, a_R_ref), eps )

                                !         ! Tangential jumps
                                !         du_t    = u_t_R    - u_t_L
                                !         dtau_nt = tau_nt_R - tau_nt_L

                                !         ! Multi-sensor:
                                !         sensor_ptot = (dSigma*dSigma)   / ( (ADC_kappa*Sigma_ref)**2 + eps )
                                !         sensor_vt   = (du_t*du_t)       / ( (ADC_kappa*a_ref    )**2 + eps )
                                !         sensor_tnt  = (dtau_nt*dtau_nt) / ( (ADC_kappa*Sigma_ref)**2 + eps )

                                !         sensor_combined = sensor_ptot + sensor_tnt + sensor_vt

                                !         phi = exp( - (sensor_combined**ADC_power) )

                                !         ! Replace F_HLLC and velocity fluxes with blended fluxes
                                !         F_HLLC   = F_HLL   + phi*(F_HLLC   - F_HLL  )
                                !         u_n_HLLC = u_n_HLL + phi*(u_n_HLLC - u_n_HLL)
                                !         u_t_HLLC = u_t_HLL + phi*(u_t_HLLC - u_t_HLL)

                                !     end if

                                !     ! === END ADC BLENDING ===


                                !     ! Assign flux_rs to F_HLLC with the right indexing order
                                !     do i = 1, num_fluids
                                !         flux_rs${XYZ}$_vf(j, k, l, i) = F_HLLC(i)
                                !         flux_rs${XYZ}$_vf(j, k, l, advxb - 1 + i) = F_HLLC(advxb - 1 + i)
                                !     end do

                                !     if (n == 0) then
                                !         flux_rs${XYZ}$_vf(j, k, l, strxb) = F_HLLC(strxb)

                                !         flux_rs${XYZ}$_vf(j, k, l, momxb) = F_HLLC(momxb)
                                !     else
                                !         if (idx1 == 1) then
                                !             flux_rs${XYZ}$_vf(j, k, l, strxb)     = F_HLLC(strxb)
                                !             flux_rs${XYZ}$_vf(j, k, l, strxb + 1) = F_HLLC(strxb + 1)
                                !             flux_rs${XYZ}$_vf(j, k, l, strxb + 2) = F_HLLC(strxb + 2)

                                !             flux_rs${XYZ}$_vf(j, k, l, momxb)     = F_HLLC(momxb)
                                !             flux_rs${XYZ}$_vf(j, k, l, momxb + 1) = F_HLLC(momxb + 1)

                                !             nc_iface_vel_rs${XYZ}$_vf(j, k, l, 1) = u_n_HLLC
                                !             nc_iface_vel_rs${XYZ}$_vf(j, k, l, 2) = u_t_HLLC
                                !         else
                                !             flux_rs${XYZ}$_vf(j, k, l, strxb)     = F_HLLC(strxb + 2)
                                !             flux_rs${XYZ}$_vf(j, k, l, strxb + 1) = F_HLLC(strxb + 1)
                                !             flux_rs${XYZ}$_vf(j, k, l, strxb + 2) = F_HLLC(strxb)

                                !             flux_rs${XYZ}$_vf(j, k, l, momxb)     = F_HLLC(momxb + 1)
                                !             flux_rs${XYZ}$_vf(j, k, l, momxb + 1) = F_HLLC(momxb)

                                !             nc_iface_vel_rs${XYZ}$_vf(j, k, l, 1) = u_t_HLLC
                                !             nc_iface_vel_rs${XYZ}$_vf(j, k, l, 2) = u_n_HLLC
                                !         end if
                                !     end if

                                !     flux_rs${XYZ}$_vf(j, k, l, E_idx) = F_HLLC(E_idx)
                                !     if (cyl_coord) then
                                !         flux_rs${XYZ}$_vf(j, k, l, strxe) = F_HLLC(strxe)
                                !     end if
                                ! end if
                                ! ! ===== END HLLC hypo star-state path =====

                                ! Geometrical source flux for cylindrical coordinates
                                #:if (NORM_DIR == 2)
                                    if (cyl_coord .and. hypoelasticity) then
                                        !$acc loop seq
                                        do i = 1, sys_size
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = flux_rs${XYZ}$_vf(j, k, l, i)
                                        end do
                                        ! HLLC face state for radial momentum gsrc
                                        if (s_L >= 0._wp) then
                                            p_face = pres_L; tau_qq_face = tau_qq_L
                                        elseif (s_R <= 0._wp) then
                                            p_face = pres_R; tau_qq_face = tau_qq_R
                                        elseif (s_S >= 0._wp) then
                                            p_face = pres_tot_star + tau_nn_L; tau_qq_face = tau_qq_L
                                        else
                                            p_face = pres_tot_star + tau_nn_R; tau_qq_face = tau_qq_R
                                        end if
                                        ! ADC blending of axisym face state
                                        if (riemann_hypo_ADC) then
                                            if (s_L >= 0._wp) then
                                                p_face_HLL = pres_L; tau_qq_face_HLL = tau_qq_L
                                            elseif (s_R <= 0._wp) then
                                                p_face_HLL = pres_R; tau_qq_face_HLL = tau_qq_R
                                            else
                                                tau_nn_HLL = U_HLL(strxb)/(rho_HLL + verysmall)
                                                tau_qq_face_HLL = U_HLL(strxe)/(rho_HLL + verysmall)
                                                p_face_HLL = F_HLL(momxb) - rho_HLL*u_n_HLL*u_n_HLL + tau_nn_HLL
                                            end if
                                            p_face = p_face_HLL + phi*(p_face - p_face_HLL)
                                            tau_qq_face = tau_qq_face_HLL + phi*(tau_qq_face - tau_qq_face_HLL)
                                        end if
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, contxe + idx1) = &
                                            flux_rs${XYZ}$_vf(j, k, l, contxe + idx1) - p_face + tau_qq_face
                                        !$acc loop seq
                                        do i = advxb, advxe
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = 0._wp
                                        end do
                                    elseif (cyl_coord) then
                                        !$acc loop seq
                                        do i = 1, E_idx
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = flux_rs${XYZ}$_vf(j, k, l, i)
                                        end do
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
                                        !$acc loop seq
                                        do i = advxb, advxe
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = 0._wp
                                        end do
                                    end if
                                #:endif
                                #:if (NORM_DIR == 3)
                                    if (grid_geometry == 3) then
                                        !$acc loop seq
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
                    !$acc end parallel loop
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
            call s_compute_capilary_source_flux( &
                q_prim_vf, &
                vel_src_rsx_vf, &
                vel_src_rsy_vf, &
                vel_src_rsz_vf, &
                flux_src_vf, &
                norm_dir, isx, isy, isz)
        end if

        call s_finalize_riemann_solver(flux_vf, flux_src_vf, &
                                       flux_gsrc_vf, &
                                       norm_dir, ix, iy, iz)

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
            dqL_prim_dy_vf, dqL_prim_dz_vf, qL_prim_vf, &
            qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, dqR_prim_dx_vf, &
            dqR_prim_dy_vf, dqR_prim_dz_vf, qR_prim_vf, &
            norm_dir, ix, iy, iz)

        call s_initialize_riemann_solver( &
            q_prim_vf, flux_vf, flux_src_vf, flux_gsrc_vf, norm_dir, ix, iy, iz)

        #:for NORM_DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
            if (norm_dir == ${NORM_DIR}$) then
                !$acc parallel loop collapse(3) gang vector default(present) &
                !$acc private(alpha_rho_L, alpha_rho_R, vel, alpha_L, alpha_R, &
                !$acc rho, pres, E, H_no_mag, gamma, pi_inf, qv, vel_rms, B, c, c_fast, pres_mag, &
                !$acc U_L, U_R, U_starL, U_starR, U_doubleL, U_doubleR, F_L, F_R, F_starL, F_starR, F_hlld)
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
                                    B%L(1) = Bx0
                                    B%R(1) = Bx0
                                    B%L(2) = qL_prim_rs${XYZ}$_vf(j, k, l, B_idx%beg)
                                    B%R(2) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, B_idx%beg)
                                    B%L(3) = qL_prim_rs${XYZ}$_vf(j, k, l, B_idx%beg + 1)
                                    B%R(3) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, B_idx%beg + 1)
                                else ! 2D/3D: Bx, By, Bz as variables
                                    B%L(1) = qL_prim_rs${XYZ}$_vf(j, k, l, B_idx%beg + dir_idx(1) - 1)
                                    B%R(1) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, B_idx%beg + dir_idx(1) - 1)
                                    B%L(2) = qL_prim_rs${XYZ}$_vf(j, k, l, B_idx%beg + dir_idx(2) - 1)
                                    B%R(2) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, B_idx%beg + dir_idx(2) - 1)
                                    B%L(3) = qL_prim_rs${XYZ}$_vf(j, k, l, B_idx%beg + dir_idx(3) - 1)
                                    B%R(3) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, B_idx%beg + dir_idx(3) - 1)
                                end if
                            end if

                            ! Sum properties of all fluid components
                            rho%L = 0._wp; gamma%L = 0._wp; pi_inf%L = 0._wp; qv%L = 0._wp
                            rho%R = 0._wp; gamma%R = 0._wp; pi_inf%R = 0._wp; qv%R = 0._wp
                            !$acc loop seq
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

                            ! (5) Compute the left/right conserved state vectors
                            U_L(1) = rho%L
                            U_L(2) = rho%L*vel%L(1)
                            U_L(3) = rho%L*vel%L(2)
                            U_L(4) = rho%L*vel%L(3)
                            U_L(5) = B%L(2)
                            U_L(6) = B%L(3)
                            U_L(7) = E%L

                            U_R(1) = rho%R
                            U_R(2) = rho%R*vel%R(1)
                            U_R(3) = rho%R*vel%R(2)
                            U_R(4) = rho%R*vel%R(3)
                            U_R(5) = B%R(2)
                            U_R(6) = B%R(3)
                            U_R(7) = E%R

                            ! (6) Compute the left/right star state vectors
                            U_starL(1) = rhoL_star
                            U_starL(2) = rhoL_star*s_M
                            U_starL(3) = rhoL_star*vel%L(2)
                            U_starL(4) = rhoL_star*vel%L(3)
                            U_starL(5) = B%L(2)
                            U_starL(6) = B%L(3)
                            U_starL(7) = E_starL

                            U_starR(1) = rhoR_star
                            U_starR(2) = rhoR_star*s_M
                            U_starR(3) = rhoR_star*vel%R(2)
                            U_starR(4) = rhoR_star*vel%R(3)
                            U_starR(5) = B%R(2)
                            U_starR(6) = B%R(3)
                            U_starR(7) = E_starR

                            ! (7) Compute the left/right fluxes
                            F_L(1) = rho%L*vel%L(1)
                            F_L(2) = rho%L*vel%L(1)*vel%L(1) - B%L(1)*B%L(1) + pTot_L
                            F_L(3) = rho%L*vel%L(1)*vel%L(2) - B%L(1)*B%L(2)
                            F_L(4) = rho%L*vel%L(1)*vel%L(3) - B%L(1)*B%L(3)
                            F_L(5) = vel%L(1)*B%L(2) - vel%L(2)*B%L(1)
                            F_L(6) = vel%L(1)*B%L(3) - vel%L(3)*B%L(1)
                            F_L(7) = (E%L + pTot_L)*vel%L(1) - B%L(1)*(vel%L(1)*B%L(1) + vel%L(2)*B%L(2) + vel%L(3)*B%L(3))

                            F_R(1) = rho%R*vel%R(1)
                            F_R(2) = rho%R*vel%R(1)*vel%R(1) - B%R(1)*B%R(1) + pTot_R
                            F_R(3) = rho%R*vel%R(1)*vel%R(2) - B%R(1)*B%R(2)
                            F_R(4) = rho%R*vel%R(1)*vel%R(3) - B%R(1)*B%R(3)
                            F_R(5) = vel%R(1)*B%R(2) - vel%R(2)*B%R(1)
                            F_R(6) = vel%R(1)*B%R(3) - vel%R(3)*B%R(1)
                            F_R(7) = (E%R + pTot_R)*vel%R(1) - B%R(1)*(vel%R(1)*B%R(1) + vel%R(2)*B%R(2) + vel%R(3)*B%R(3))

                            ! (8) Compute the left/right star fluxes (note array operations)
                            F_starL = F_L + s_L*(U_starL - U_L)
                            F_starR = F_R + s_R*(U_starR - U_R)

                            ! (9) Compute the rotational (Alfvén) speeds
                            s_starL = s_M - abs(B%L(1))/sqrt(rhoL_star)
                            s_starR = s_M + abs(B%L(1))/sqrt(rhoR_star)

                            ! (10) Compute the double–star states [Miyoshi Eqns. (59)-(62)]
                            sqrt_rhoL_star = sqrt(rhoL_star)
                            sqrt_rhoR_star = sqrt(rhoR_star)
                            denom_ds = sqrt_rhoL_star + sqrt_rhoR_star
                            sign_Bx = sign(1._wp, B%L(1))
                            vL_star = vel%L(2)
                            wL_star = vel%L(3)
                            vR_star = vel%R(2)
                            wR_star = vel%R(3)
                            v_double = (sqrt_rhoL_star*vL_star + sqrt_rhoR_star*vR_star + (B%R(2) - B%L(2))*sign_Bx)/denom_ds
                            w_double = (sqrt_rhoL_star*wL_star + sqrt_rhoR_star*wR_star + (B%R(3) - B%L(3))*sign_Bx)/denom_ds
                            By_double = (sqrt_rhoL_star*B%R(2) + sqrt_rhoR_star*B%L(2) + sqrt_rhoL_star*sqrt_rhoR_star*(vR_star - vL_star)*sign_Bx)/denom_ds
                            Bz_double = (sqrt_rhoL_star*B%R(3) + sqrt_rhoR_star*B%L(3) + sqrt_rhoL_star*sqrt_rhoR_star*(wR_star - wL_star)*sign_Bx)/denom_ds

                            E_doubleL = E_starL - sqrt_rhoL_star*((vL_star*B%L(2) + wL_star*B%L(3)) - (v_double*By_double + w_double*Bz_double))*sign_Bx
                            E_doubleR = E_starR + sqrt_rhoR_star*((vR_star*B%R(2) + wR_star*B%R(3)) - (v_double*By_double + w_double*Bz_double))*sign_Bx
                            E_double = 0.5_wp*(E_doubleL + E_doubleR)

                            U_doubleL(1) = rhoL_star
                            U_doubleL(2) = rhoL_star*s_M
                            U_doubleL(3) = rhoL_star*v_double
                            U_doubleL(4) = rhoL_star*w_double
                            U_doubleL(5) = By_double
                            U_doubleL(6) = Bz_double
                            U_doubleL(7) = E_double

                            U_doubleR(1) = rhoR_star
                            U_doubleR(2) = rhoR_star*s_M
                            U_doubleR(3) = rhoR_star*v_double
                            U_doubleR(4) = rhoR_star*w_double
                            U_doubleR(5) = By_double
                            U_doubleR(6) = Bz_double
                            U_doubleR(7) = E_double

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
                            flux_rs${XYZ}$_vf(j, k, l, contxe + dir_idx(1)) = F_hlld(2)
                            flux_rs${XYZ}$_vf(j, k, l, contxe + dir_idx(2)) = F_hlld(3)
                            flux_rs${XYZ}$_vf(j, k, l, contxe + dir_idx(3)) = F_hlld(4)
                            ! Magnetic field
                            if (n == 0) then
                                flux_rs${XYZ}$_vf(j, k, l, B_idx%beg) = F_hlld(5)
                                flux_rs${XYZ}$_vf(j, k, l, B_idx%beg + 1) = F_hlld(6)
                            else
                                flux_rs${XYZ}$_vf(j, k, l, B_idx%beg + dir_idx(2) - 1) = F_hlld(5)
                                flux_rs${XYZ}$_vf(j, k, l, B_idx%beg + dir_idx(3) - 1) = F_hlld(6)
                            end if
                            ! Energy
                            flux_rs${XYZ}$_vf(j, k, l, E_idx) = F_hlld(7)
                            ! Partial fraction
                            !$acc loop seq
                            do i = advxb, advxe
                                flux_rs${XYZ}$_vf(j, k, l, i) = 0._wp ! TODO multi-component (zero for now)
                            end do

                            flux_src_rs${XYZ}$_vf(j, k, l, advxb) = 0._wp
                        end do
                    end do
                end do
                !$acc end parallel loop
            end if
        #:endfor

        call s_finalize_riemann_solver(flux_vf, flux_src_vf, flux_gsrc_vf, &
                                       norm_dir, ix, iy, iz)
    end subroutine s_hlld_riemann_solver

    !> HLLD Riemann solver resolves all 5 waves for the hypoelastic equations:
        !!      1 entropy wave, 2 shear stress waves, 2 fast waves.
    subroutine s_hypo_hlld_riemann_solver(qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, &
                                          dqL_prim_dx_vf, dqL_prim_dy_vf, dqL_prim_dz_vf, &
                                          qL_prim_vf, &
                                          qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, &
                                          dqR_prim_dx_vf, dqR_prim_dy_vf, dqR_prim_dz_vf, &
                                          qR_prim_vf, &
                                          q_prim_vf, &
                                          flux_vf, flux_src_vf, flux_gsrc_vf, &
                                          norm_dir, ix, iy, iz, is_hat_L)

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

        logical, intent(in) :: is_hat_L

        ! Local variables:

        ! q_prim_vf with indiced permutated to match qL_prim_rs${XYZ}$_vf
        ! But note that q_prim_vf and q_hat_prim_x_vf are cell-centered while qL_prim_rs${XYZ}$_vf etc. are cell-boundary values
        real(wp), dimension(idwbuff(1)%beg:idwbuff(1)%end, &
                            idwbuff(2)%beg:idwbuff(2)%end, &
                            idwbuff(3)%beg:idwbuff(3)%end, &
                            1:sys_size) :: q_hat_prim_x_vf

        real(wp), dimension(idwbuff(2)%beg:idwbuff(2)%end, &
                            idwbuff(1)%beg:idwbuff(1)%end, &
                            idwbuff(3)%beg:idwbuff(3)%end, &
                            1:sys_size) :: q_hat_prim_y_vf

        real(wp), dimension(idwbuff(3)%beg:idwbuff(3)%end, &
                            idwbuff(2)%beg:idwbuff(2)%end, &
                            idwbuff(1)%beg:idwbuff(1)%end, &
                            1:sys_size) :: q_hat_prim_z_vf

        real(wp), dimension(num_fluids) :: alpha_L, alpha_R, alpha_rho_L, alpha_rho_R
        type(riemann_states_vec3) :: vel
        type(riemann_states) :: rho, pres, E, H
        type(riemann_states) :: gamma, pi_inf, qv
        type(riemann_states) :: vel_rms

        type(riemann_states) :: c

        ! HLLD speeds and intermediate state variables:
        real(wp) :: S_L, S_R, s_M, S_Lstar, S_Rstar
        real(wp) :: pTot_L, pTot_R, rhoL_star, rhoR_star

        real(wp), dimension(14) :: U_L, U_R, U_starL, U_starR, U_starstarL, U_starstarR
        real(wp), dimension(14) :: F_L, F_R, F_starL, F_starR, F_hlld
        real(wp), dimension(14) :: F_HLL  ! for ADC blending
        integer :: ncomp  ! 11 for 2D/axisym, 14 for 3D Cartesian

        ! HLLD Hypo variables

        real(wp) :: C_NC, sqrtC_NC
        real(wp) :: A_L, A_R, denomA, fac_L, fac_R
        real(wp) :: u_n_L, u_t_L, u_n_R, u_t_R
        real(wp) :: u_t2_L, u_t2_R
        real(wp) :: tau_nn_L, tau_nt_L, tau_tt_L, tau_nn_R, tau_nt_R, tau_tt_R
        real(wp) :: tau_nt2_L, tau_nt2_R, tau_t2t2_L, tau_t2t2_R, tau_t1t2_L, tau_t1t2_R
        real(wp) :: tau_qq_L, tau_qq_R

        real(wp) :: G_L, G_R
        real(wp), dimension(strxe-strxb+1) :: tau_e_L, tau_e_R

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
        real(wp) :: u_n_hat, u_t_hat, u_t2_hat
        real(wp) :: tau_nn_hat, tau_nt_hat, tau_tt_hat, tau_qq_hat
        real(wp) :: tau_nt2_hat, tau_t2t2_hat, tau_t1t2_hat
        real(wp), dimension(num_fluids) :: alpha_hat, alpha_rho_hat
        real(wp), dimension(num_vels) :: vel_hat
        real(wp), dimension(strxe-strxb+1) :: tau_e_hat

        real(wp) :: pres_hat, blkmod1_hat, blkmod2_hat, K_hat
        real(wp) :: C_hat_1, C_hat_2

        real(wp) :: Sigma_L, Sigma_R, dSigma, Sigma_ref
        real(wp) :: a_L_ref, a_R_ref, a_ref
        real(wp) :: du_t, dtau_nt, du_t2, dtau_nt2
        real(wp) :: sensor_ptot, sensor_vt, sensor_tnt, sensor_combined
        real(wp) :: phi

        real(wp), parameter :: ADC_power = 1.0_wp

        real(wp) :: alpha_L_sum, alpha_R_sum

        integer :: i, j, k, l

        call s_populate_riemann_states_variables_buffers( &
            qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, dqL_prim_dx_vf, &
            dqL_prim_dy_vf, dqL_prim_dz_vf, qL_prim_vf, &
            qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, dqR_prim_dx_vf, &
            dqR_prim_dy_vf, dqR_prim_dz_vf, qR_prim_vf, &
            norm_dir, ix, iy, iz)

        call s_initialize_riemann_solver( &
            q_prim_vf, flux_vf, flux_src_vf, flux_gsrc_vf, norm_dir, ix, iy, iz)

        #:for NORM_DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
            if (norm_dir == ${NORM_DIR}$) then

                ! Fill q_hat_prim with cell-centered primitives at the interface location.
                !$acc data create(q_hat_prim_${XYZ}$_vf)
                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end
                                #:if NORM_DIR == 1
                                if (is_hat_L) then
                                    q_hat_prim_x_vf(j,k,l,i) = q_prim_vf(i)%sf(j  ,k,l)
                                else
                                    q_hat_prim_x_vf(j,k,l,i) = q_prim_vf(i)%sf(j+1,k,l)
                                end if
                                #:elif NORM_DIR == 2
                                if (is_hat_L) then
                                    q_hat_prim_y_vf(j,k,l,i) = q_prim_vf(i)%sf(k,j  ,l)
                                else
                                    q_hat_prim_y_vf(j,k,l,i) = q_prim_vf(i)%sf(k,j+1,l)
                                end if
                                #:else
                                if (is_hat_L) then
                                    q_hat_prim_z_vf(j,k,l,i) = q_prim_vf(i)%sf(l,k,j  )
                                else
                                    q_hat_prim_z_vf(j,k,l,i) = q_prim_vf(i)%sf(l,k,j+1)
                                end if
                                #:endif
                            end do
                        end do
                    end do
                end do

                !$acc parallel loop collapse(3) gang vector default(present) &
                !$acc private(alpha_rho_L, alpha_rho_R, vel, alpha_L, alpha_R, &
                !$acc rho, pres, E, H, gamma, pi_inf, qv, vel_rms, c, &
                !$acc U_L, U_R, U_starL, U_starR, U_starstarL, U_starstarR, &
                !$acc F_L, F_R, F_starL, F_starR, F_hlld, F_HLL, &
                !$acc tau_e_L, tau_e_R, tau_e_hat, &
                !$acc alpha_hat, alpha_rho_hat, vel_hat)
                do l = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = is1%beg, is1%end

                            ! ==================== Extract left/right primitive states ====================

                            do i = 1, contxe
                                alpha_rho_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                alpha_rho_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)

                                alpha_rho_hat(i) = q_hat_prim_${XYZ}$_vf(j, k, l, i)
                            end do

                            ! IMP: vel%L(1:3) has (3) uninitiated for 2D
                            vel%L = 0._wp
                            vel%R = 0._wp

                            ! NOTE: unlike HLL & HLLC, vel%L here is permutated by dir_idx for simpler logic
                            do i = 1, num_vels
                                vel%L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, contxe + i) ! Don't permutate here; permutate u <-> v later at u_n_L = vel%L(1)
                                vel%R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, contxe + i)

                                vel_hat(i) = q_hat_prim_${XYZ}$_vf(j, k, l, contxe + i)
                            end do

                            vel_rms%L = vel%L(1)**2._wp + vel%L(2)**2._wp + vel%L(3)**2._wp
                            vel_rms%R = vel%R(1)**2._wp + vel%R(2)**2._wp + vel%R(3)**2._wp

                            do i = 1, num_fluids
                                alpha_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)
                                alpha_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)

                                alpha_hat(i) = q_hat_prim_${XYZ}$_vf(j, k, l, E_idx + i)
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

                            pres%L = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx)
                            pres%R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx)

                            ! Hypoelasticity
                            !$acc loop seq
                            do i = 1, strxe - strxb + 1
                                tau_e_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, strxb - 1 + i)
                                tau_e_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, strxb - 1 + i)

                                tau_e_hat(i) = q_hat_prim_${XYZ}$_vf(j, k, l, strxb - 1 + i)
                            end do

                            u_t2_L = 0._wp; u_t2_R = 0._wp; u_t2_hat = 0._wp
                            tau_nt2_L = 0._wp; tau_nt2_R = 0._wp; tau_nt2_hat = 0._wp
                            tau_t2t2_L = 0._wp; tau_t2t2_R = 0._wp; tau_t2t2_hat = 0._wp
                            tau_t1t2_L = 0._wp; tau_t1t2_R = 0._wp; tau_t1t2_hat = 0._wp

                            if (n == 0) then
                                ncomp = 11
                                tau_nn_L = tau_e_L(1)
                                tau_nn_R = tau_e_R(1)
                                u_n_L = vel%L(1)
                                u_n_R = vel%R(1)
                                u_n_hat = vel_hat(1)
                                u_t_L = 0._wp; u_t_R = 0._wp; u_t_hat = 0._wp
                                tau_nt_L = 0._wp; tau_nt_R = 0._wp; tau_nt_hat = 0._wp
                                tau_tt_L = 0._wp; tau_tt_R = 0._wp; tau_tt_hat = 0._wp
                            else if (p == 0) then
                                ncomp = 11
                                tau_nt_L = tau_e_L(2)
                                tau_nt_R = tau_e_R(2)
                                tau_nt_hat = tau_e_hat(2)
                                if (dir_idx(1) == 1) then
                                    tau_nn_L = tau_e_L(1); tau_nn_R = tau_e_R(1)
                                    tau_tt_L = tau_e_L(3); tau_tt_R = tau_e_R(3)
                                    u_n_L = vel%L(1); u_n_R = vel%R(1)
                                    u_t_L = vel%L(2); u_t_R = vel%R(2)
                                    tau_nn_hat = tau_e_hat(1); tau_tt_hat = tau_e_hat(3)
                                    u_n_hat = vel_hat(1); u_t_hat = vel_hat(2)
                                else
                                    tau_nn_L = tau_e_L(3); tau_nn_R = tau_e_R(3)
                                    tau_tt_L = tau_e_L(1); tau_tt_R = tau_e_R(1)
                                    u_n_L = vel%L(2); u_n_R = vel%R(2)
                                    u_t_L = vel%L(1); u_t_R = vel%R(1)
                                    tau_nn_hat = tau_e_hat(3); tau_tt_hat = tau_e_hat(1)
                                    u_n_hat = vel_hat(2); u_t_hat = vel_hat(1)
                                end if
                            else
                                ncomp = 14
                                if (dir_idx(1) == 1) then
                                    u_n_L = vel%L(1); u_n_R = vel%R(1)
                                    u_t_L = vel%L(2); u_t_R = vel%R(2)
                                    u_t2_L = vel%L(3); u_t2_R = vel%R(3)
                                    tau_nn_L = tau_e_L(1); tau_nn_R = tau_e_R(1)
                                    tau_nt_L = tau_e_L(2); tau_nt_R = tau_e_R(2)
                                    tau_nt2_L = tau_e_L(4); tau_nt2_R = tau_e_R(4)
                                    tau_tt_L = tau_e_L(3); tau_tt_R = tau_e_R(3)
                                    tau_t2t2_L = tau_e_L(6); tau_t2t2_R = tau_e_R(6)
                                    tau_t1t2_L = tau_e_L(5); tau_t1t2_R = tau_e_R(5)
                                    u_n_hat = vel_hat(1); u_t_hat = vel_hat(2); u_t2_hat = vel_hat(3)
                                    tau_nn_hat = tau_e_hat(1); tau_nt_hat = tau_e_hat(2)
                                    tau_nt2_hat = tau_e_hat(4); tau_tt_hat = tau_e_hat(3)
                                    tau_t2t2_hat = tau_e_hat(6); tau_t1t2_hat = tau_e_hat(5)
                                else if (dir_idx(1) == 2) then
                                    u_n_L = vel%L(2); u_n_R = vel%R(2)
                                    u_t_L = vel%L(1); u_t_R = vel%R(1)
                                    u_t2_L = vel%L(3); u_t2_R = vel%R(3)
                                    tau_nn_L = tau_e_L(3); tau_nn_R = tau_e_R(3)
                                    tau_nt_L = tau_e_L(2); tau_nt_R = tau_e_R(2)
                                    tau_nt2_L = tau_e_L(5); tau_nt2_R = tau_e_R(5)
                                    tau_tt_L = tau_e_L(1); tau_tt_R = tau_e_R(1)
                                    tau_t2t2_L = tau_e_L(6); tau_t2t2_R = tau_e_R(6)
                                    tau_t1t2_L = tau_e_L(4); tau_t1t2_R = tau_e_R(4)
                                    u_n_hat = vel_hat(2); u_t_hat = vel_hat(1); u_t2_hat = vel_hat(3)
                                    tau_nn_hat = tau_e_hat(3); tau_nt_hat = tau_e_hat(2)
                                    tau_nt2_hat = tau_e_hat(5); tau_tt_hat = tau_e_hat(1)
                                    tau_t2t2_hat = tau_e_hat(6); tau_t1t2_hat = tau_e_hat(4)
                                else
                                    u_n_L = vel%L(3); u_n_R = vel%R(3)
                                    u_t_L = vel%L(1); u_t_R = vel%R(1)
                                    u_t2_L = vel%L(2); u_t2_R = vel%R(2)
                                    tau_nn_L = tau_e_L(6); tau_nn_R = tau_e_R(6)
                                    tau_nt_L = tau_e_L(4); tau_nt_R = tau_e_R(4)
                                    tau_nt2_L = tau_e_L(5); tau_nt2_R = tau_e_R(5)
                                    tau_tt_L = tau_e_L(1); tau_tt_R = tau_e_R(1)
                                    tau_t2t2_L = tau_e_L(3); tau_t2t2_R = tau_e_R(3)
                                    tau_t1t2_L = tau_e_L(2); tau_t1t2_R = tau_e_R(2)
                                    u_n_hat = vel_hat(3); u_t_hat = vel_hat(1); u_t2_hat = vel_hat(2)
                                    tau_nn_hat = tau_e_hat(6); tau_nt_hat = tau_e_hat(4)
                                    tau_nt2_hat = tau_e_hat(5); tau_tt_hat = tau_e_hat(1)
                                    tau_t2t2_hat = tau_e_hat(3); tau_t1t2_hat = tau_e_hat(2)
                                end if
                            end if
                            if (cyl_coord) then
                                tau_qq_L = tau_e_L(strxe - strxb + 1)
                                tau_qq_R = tau_e_R(strxe - strxb + 1)
                                tau_qq_hat = tau_e_hat(strxe - strxb + 1)
                            else
                                tau_qq_L = 0._wp; tau_qq_R = 0._wp; tau_qq_hat = 0._wp
                            end if
                            ! Total pressure (replace the usual pressure to define SM)
                            pTot_L = pres%L - tau_nn_L
                            pTot_R = pres%R - tau_nn_R

                            ! Sum properties of all fluid components
                            rho%L = 0._wp; gamma%L = 0._wp; pi_inf%L = 0._wp; qv%L = 0._wp
                            rho%R = 0._wp; gamma%R = 0._wp; pi_inf%R = 0._wp; qv%R = 0._wp
                            rho_hat = 0._wp
                            !$acc loop seq
                            do i = 1, num_fluids
                                rho%L = rho%L + alpha_rho_L(i)
                                gamma%L = gamma%L + alpha_L(i)*gammas(i)
                                pi_inf%L = pi_inf%L + alpha_L(i)*pi_infs(i)
                                qv%L = qv%L + alpha_rho_L(i)*qvs(i)

                                rho%R = rho%R + alpha_rho_R(i)
                                gamma%R = gamma%R + alpha_R(i)*gammas(i)
                                pi_inf%R = pi_inf%R + alpha_R(i)*pi_infs(i)
                                qv%R = qv%R + alpha_rho_R(i)*qvs(i)

                                rho_hat = rho_hat + alpha_rho_hat(i)
                            end do


                            G_L = 0._wp; G_R = 0._wp; G_hat = 0._wp
                            !$acc loop seq
                            do i = 1, num_fluids
                                G_L = G_L + alpha_L(i)*Gs(i)
                                G_R = G_R + alpha_R(i)*Gs(i)

                                G_hat = G_hat + alpha_hat(i)*Gs(i)
                            end do

                            E%L = gamma%L*pres%L + pi_inf%L + 5e-1*rho%L*vel_rms%L + qv%L
                            E%R = gamma%R*pres%R + pi_inf%R + 5e-1*rho%R*vel_rms%R + qv%R

                            !$acc loop seq
                            do i = 1, strxe - strxb + 1
                                ! Elastic contribution to energy if G large enough
                                if ((G_L > verysmall) .and. (G_R > verysmall)) then
                                    E%L = E%L + (tau_e_L(i)*tau_e_L(i))/(4._wp*G_L)
                                    E%R = E%R + (tau_e_R(i)*tau_e_R(i))/(4._wp*G_R)
                                    ! Shear terms doubled: 2D/2D-axisym i==2 only; 3D i==2,4,5
                                    if ((n > 0 .and. p == 0 .and. i == 2) .or. &
                                        (p > 0 .and. (i == 2 .or. i == 4 .or. i == 5))) then
                                        E%L = E%L + (tau_e_L(i)*tau_e_L(i))/(4._wp*G_L)
                                        E%R = E%R + (tau_e_R(i)*tau_e_R(i))/(4._wp*G_R)
                                    end if
                                end if
                            end do

                            H%L = (E%L + pres%L)/rho%L
                            H%R = (E%R + pres%R)/rho%R

                            ! ==================== Compute Riemann states ====================

                            call s_compute_speed_of_sound(pres%L, rho%L, gamma%L, pi_inf%L, H%L, alpha_L, vel_rms%L, 0._wp, c%L)
                            call s_compute_speed_of_sound(pres%R, rho%R, gamma%R, pi_inf%R, H%R, alpha_R, vel_rms%R, 0._wp, c%R)

                            S_L = min(u_n_L - sqrt(c%L*c%L + ((4._wp/3._wp)*G_L + tau_nn_L)/rho%L), u_n_R - sqrt(c%R*c%R + ((4._wp/3._wp)*G_R + tau_nn_R)/rho%R))
                            S_R = max(u_n_R + sqrt(c%R*c%R + ((4._wp/3._wp)*G_R + tau_nn_R)/rho%R), u_n_L + sqrt(c%L*c%L + ((4._wp/3._wp)*G_L + tau_nn_L)/rho%L))

                            ! Two-component 2D only (enforced by checker restrictions)

                            K_hat = 0._wp
                            if (alt_soundspeed) then
                                pres_hat = q_hat_prim_${XYZ}$_vf(j, k, l, E_idx)
                                blkmod1_hat = ((gammas(1) + 1._wp)*pres_hat + pi_infs(1))/gammas(1)
                                blkmod2_hat = ((gammas(2) + 1._wp)*pres_hat + pi_infs(2))/gammas(2)
                                K_hat = alpha_hat(1)*alpha_hat(2)*(blkmod2_hat - blkmod1_hat) &
                                        /(alpha_hat(1)*blkmod2_hat + alpha_hat(2)*blkmod1_hat + verysmall)
                            end if
                            C_hat_1 = alpha_hat(1) + K_hat
                            C_hat_2 = alpha_hat(2) - K_hat

                            if (p > 0 .and. .not. cyl_coord) then
                                ! 3D Cartesian: 14-state compact basis
                                U_L(1) = alpha_rho_L(1); U_R(1) = alpha_rho_R(1)
                                U_L(2) = alpha_rho_L(2); U_R(2) = alpha_rho_R(2)
                                U_L(3) = rho%L*u_n_L;   U_R(3) = rho%R*u_n_R
                                U_L(4) = rho%L*u_t_L;   U_R(4) = rho%R*u_t_R
                                U_L(5) = rho%L*u_t2_L;  U_R(5) = rho%R*u_t2_R
                                U_L(6) = E%L;            U_R(6) = E%R
                                U_L(7) = alpha_L(1);     U_R(7) = alpha_R(1)
                                U_L(8) = alpha_L(2);     U_R(8) = alpha_R(2)
                                U_L(9) = rho%L*tau_nn_L;     U_R(9) = rho%R*tau_nn_R
                                U_L(10) = rho%L*tau_nt_L;    U_R(10) = rho%R*tau_nt_R
                                U_L(11) = rho%L*tau_nt2_L;   U_R(11) = rho%R*tau_nt2_R
                                U_L(12) = rho%L*tau_tt_L;    U_R(12) = rho%R*tau_tt_R
                                U_L(13) = rho%L*tau_t2t2_L;  U_R(13) = rho%R*tau_t2t2_R
                                U_L(14) = rho%L*tau_t1t2_L;  U_R(14) = rho%R*tau_t1t2_R

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
                                          - 2._wp*rho_hat*tau_nt_hat*u_t_L
                                F_R(12) = U_R(12)*u_n_R + rho_hat*(2._wp/3._wp*G_hat + tau_tt_hat)*u_n_R &
                                          - 2._wp*rho_hat*tau_nt_hat*u_t_R
                                F_L(13) = U_L(13)*u_n_L + rho_hat*(2._wp/3._wp*G_hat + tau_t2t2_hat)*u_n_L &
                                          - 2._wp*rho_hat*tau_nt2_hat*u_t2_L
                                F_R(13) = U_R(13)*u_n_R + rho_hat*(2._wp/3._wp*G_hat + tau_t2t2_hat)*u_n_R &
                                          - 2._wp*rho_hat*tau_nt2_hat*u_t2_R
                                F_L(14) = U_L(14)*u_n_L + rho_hat*tau_t1t2_hat*u_n_L &
                                          - rho_hat*tau_nt2_hat*u_t_L - rho_hat*tau_nt_hat*u_t2_L
                                F_R(14) = U_R(14)*u_n_R + rho_hat*tau_t1t2_hat*u_n_R &
                                          - rho_hat*tau_nt2_hat*u_t_R - rho_hat*tau_nt_hat*u_t2_R
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
                                F_L(6) = U_L(6)*u_n_L - C_hat_1*u_n_L
                                F_L(7) = U_L(7)*u_n_L - C_hat_2*u_n_L
                                F_L(8) = U_L(8)*u_n_L - rho_hat*(4._wp/3._wp*G_hat + tau_nn_hat)*u_n_L
                                F_L(9) = U_L(9)*u_n_L - rho_hat*(G_hat + tau_nn_hat)*u_t_L
                                F_L(10) = U_L(10)*u_n_L + rho_hat*(2._wp/3._wp*G_hat + tau_tt_hat)*u_n_L - 2._wp*rho_hat*tau_nt_hat*u_t_L
                                F_L(11) = U_L(11)*u_n_L + rho_hat*(2._wp/3._wp*G_hat + tau_qq_hat)*u_n_L

                                F_R(1) = U_R(1)*u_n_R
                                F_R(2) = U_R(2)*u_n_R
                                F_R(3) = rho%R*u_n_R*u_n_R + pTot_R
                                F_R(4) = rho%R*u_n_R*u_t_R - tau_nt_R
                                F_R(5) = (E%R + pTot_R)*u_n_R - u_t_R*tau_nt_R
                                F_R(6) = U_R(6)*u_n_R - C_hat_1*u_n_R
                                F_R(7) = U_R(7)*u_n_R - C_hat_2*u_n_R
                                F_R(8) = U_R(8)*u_n_R - rho_hat*(4._wp/3._wp*G_hat + tau_nn_hat)*u_n_R
                                F_R(9) = U_R(9)*u_n_R - rho_hat*(G_hat + tau_nn_hat)*u_t_R
                                F_R(10) = U_R(10)*u_n_R + rho_hat*(2._wp/3._wp*G_hat + tau_tt_hat)*u_n_R - 2._wp*rho_hat*tau_nt_hat*u_t_R
                                F_R(11) = U_R(11)*u_n_R + rho_hat*(2._wp/3._wp*G_hat + tau_qq_hat)*u_n_R
                            end if

                            A_L = rho%L*(S_L - u_n_L)
                            A_R = rho%R*(S_R - u_n_R)
                            denomA = (A_R - A_L)

                            S_M = ((pTot_R - pTot_L) + A_L*u_n_L - A_R*u_n_R) / (A_L - A_R + verysmall)

                            ! Degenerate wave structure: denom ~ 0 or S_M not in [S_L,S_R]
                            if (abs(denomA) < verysmall .or. &
                                .not. (S_L - verysmall <= S_M .and. S_M <= S_R + verysmall)) then
                                ! HLL (or one-sided) fallback for degenerate wave structure
                                if (S_L < 0._wp .and. S_R > 0._wp) then
                                    do i = 1, ncomp
                                        F_hlld(i) = (S_R*F_L(i) - S_L*F_R(i) + S_L*S_R*(U_R(i) - U_L(i))) &
                                                    /(S_R - S_L + verysmall)
                                    end do
                                elseif (S_L >= 0._wp) then
                                    F_hlld(1:ncomp) = F_L(1:ncomp)
                                else
                                    F_hlld(1:ncomp) = F_R(1:ncomp)
                                end if
                                ! Initialize star-state variables to safe fallback values so
                                ! subsequent face-state and axisymmetric source logic never
                                ! reads uninitialized memory regardless of wave configuration.
                                pTot_star     = 5e-1_wp*(pTot_L + pTot_R)
                                S_Lstar       = S_L
                                S_Rstar       = S_R
                                u_t_star      = 5e-1_wp*(u_t_L + u_t_R)
                                tau_nn_L_star = tau_nn_L
                                tau_nn_R_star = tau_nn_R
                                tau_qq_L_star = tau_qq_L
                                tau_qq_R_star = tau_qq_R
                            else

                            pTot_star = pTot_L + A_L*(S_M - u_n_L)

                            rhoL_star = rho%L * (S_L - u_n_L) / (S_L - S_M + verysmall)
                            rhoR_star = rho%R * (S_R - u_n_R) / (S_R - S_M + verysmall)
                            fac_L = (S_L - u_n_L) / (S_L - S_M + verysmall)
                            fac_R = (S_R - u_n_R) / (S_R - S_M + verysmall)

                            C_NC = rho_hat*(G_hat + tau_nn_hat)
                            sqrtC_NC = sqrt(max(C_NC, verysmall))

                            S_Lstar = S_M - sqrtC_NC/rhoL_star
                            S_Rstar = S_M + sqrtC_NC/rhoR_star

                            u_t_star    = 0.5_wp*((tau_nt_R - tau_nt_L)/sqrtC_NC + (u_t_R + u_t_L))
                            tau_nt_star = 0.5_wp*((u_t_R - u_t_L)*sqrtC_NC + (tau_nt_R + tau_nt_L))

                            u_t2_star    = 0.5_wp*((tau_nt2_R - tau_nt2_L)/sqrtC_NC + (u_t2_R + u_t2_L))
                            tau_nt2_star = 0.5_wp*((u_t2_R - u_t2_L)*sqrtC_NC + (tau_nt2_R + tau_nt2_L))

                            tau_nn_L_star = tau_nn_L - (rho_hat*(G_hat*4._wp/3._wp + tau_nn_hat)*(u_n_L - S_M)) &
                                            /(rho%L*(u_n_L - S_L))
                            tau_nn_R_star = tau_nn_R - (rho_hat*(G_hat*4._wp/3._wp + tau_nn_hat)*(u_n_R - S_M)) &
                                            /(rho%R*(u_n_R - S_R))

                            tau_tt_L_star = tau_tt_L + (rho_hat*(G_hat*2._wp/3._wp + tau_tt_hat)*(u_n_L - S_M)) &
                                            /(rho%L*(u_n_L - S_L))
                            tau_tt_R_star = tau_tt_R + (rho_hat*(G_hat*2._wp/3._wp + tau_tt_hat)*(u_n_R - S_M)) &
                                            /(rho%R*(u_n_R - S_R))

                            tau_t2t2_L_star = tau_t2t2_L + (rho_hat*(G_hat*2._wp/3._wp + tau_t2t2_hat)*(u_n_L - S_M)) &
                                              /(rho%L*(u_n_L - S_L))
                            tau_t2t2_R_star = tau_t2t2_R + (rho_hat*(G_hat*2._wp/3._wp + tau_t2t2_hat)*(u_n_R - S_M)) &
                                              /(rho%R*(u_n_R - S_R))

                            tau_t1t2_L_star = tau_t1t2_L + (rho_hat*tau_t1t2_hat*(u_n_L - S_M)) &
                                              /(rho%L*(u_n_L - S_L))
                            tau_t1t2_R_star = tau_t1t2_R + (rho_hat*tau_t1t2_hat*(u_n_R - S_M)) &
                                              /(rho%R*(u_n_R - S_R))

                            tau_qq_L_star = tau_qq_L + (rho_hat*(G_hat*2._wp/3._wp + tau_qq_hat)*(u_n_L - S_M)) &
                                            /(rho%L*(u_n_L - S_L))
                            tau_qq_R_star = tau_qq_R + (rho_hat*(G_hat*2._wp/3._wp + tau_qq_hat)*(u_n_R - S_M)) &
                                            /(rho%R*(u_n_R - S_R))

                            tau_tt_L_starstar = tau_tt_L_star + 2._wp*rho_hat*tau_nt_hat/sqrtC_NC*(u_t_star - u_t_L)
                            tau_tt_R_starstar = tau_tt_R_star - 2._wp*rho_hat*tau_nt_hat/sqrtC_NC*(u_t_star - u_t_R)

                            tau_t2t2_L_starstar = tau_t2t2_L_star + 2._wp*rho_hat*tau_nt2_hat/sqrtC_NC*(u_t2_star - u_t2_L)
                            tau_t2t2_R_starstar = tau_t2t2_R_star - 2._wp*rho_hat*tau_nt2_hat/sqrtC_NC*(u_t2_star - u_t2_R)

                            tau_t1t2_L_starstar = tau_t1t2_L_star &
                                + rho_hat*(tau_nt2_hat*(u_t_star - u_t_L) + tau_nt_hat*(u_t2_star - u_t2_L))/sqrtC_NC
                            tau_t1t2_R_starstar = tau_t1t2_R_star &
                                - rho_hat*(tau_nt2_hat*(u_t_star - u_t_R) + tau_nt_hat*(u_t2_star - u_t2_R))/sqrtC_NC

                            E_L_star = (E%L*(u_n_L - S_L) + u_n_L*pTot_L - S_M*pTot_star)/(S_M - S_L)
                            E_R_star = (E%R*(u_n_R - S_R) + u_n_R*pTot_R - S_M*pTot_star)/(S_M - S_R)
                            E_L_starstar = E_L_star + (rhoL_star/sqrtC_NC)* &
                                ((u_t_star*tau_nt_star - u_t_L*tau_nt_L) + (u_t2_star*tau_nt2_star - u_t2_L*tau_nt2_L))
                            E_R_starstar = E_R_star - (rhoR_star/sqrtC_NC)* &
                                ((u_t_star*tau_nt_star - u_t_R*tau_nt_R) + (u_t2_star*tau_nt2_star - u_t2_R*tau_nt2_R))

                            alpha1_L_star = (alpha_L(1)*(S_L - u_n_L) - C_hat_1*(S_M - u_n_L))/(S_L - S_M)
                            alpha1_R_star = (alpha_R(1)*(S_R - u_n_R) - C_hat_1*(S_M - u_n_R))/(S_R - S_M)

                            alpha2_L_star = (alpha_L(2)*(S_L - u_n_L) - C_hat_2*(S_M - u_n_L))/(S_L - S_M)
                            alpha2_R_star = (alpha_R(2)*(S_R - u_n_R) - C_hat_2*(S_M - u_n_R))/(S_R - S_M)

                            ! ==================== Compute U ====================

                            if (p > 0 .and. .not. cyl_coord) then
                                ! 3D Cartesian: 14-state
                                U_starL(1) = U_L(1) * fac_L
                                U_starL(2) = U_L(2) * fac_L
                                U_starL(3) = rhoL_star * S_M
                                U_starL(4) = rhoL_star * u_t_L
                                U_starL(5) = rhoL_star * u_t2_L
                                U_starL(6) = E_L_star
                                U_starL(7) = alpha1_L_star
                                U_starL(8) = alpha2_L_star
                                U_starL(9) = rhoL_star * tau_nn_L_star
                                U_starL(10) = rhoL_star * tau_nt_L
                                U_starL(11) = rhoL_star * tau_nt2_L
                                U_starL(12) = rhoL_star * tau_tt_L_star
                                U_starL(13) = rhoL_star * tau_t2t2_L_star
                                U_starL(14) = rhoL_star * tau_t1t2_L_star

                                U_starR(1) = U_R(1) * fac_R
                                U_starR(2) = U_R(2) * fac_R
                                U_starR(3) = rhoR_star * S_M
                                U_starR(4) = rhoR_star * u_t_R
                                U_starR(5) = rhoR_star * u_t2_R
                                U_starR(6) = E_R_star
                                U_starR(7) = alpha1_R_star
                                U_starR(8) = alpha2_R_star
                                U_starR(9) = rhoR_star * tau_nn_R_star
                                U_starR(10) = rhoR_star * tau_nt_R
                                U_starR(11) = rhoR_star * tau_nt2_R
                                U_starR(12) = rhoR_star * tau_tt_R_star
                                U_starR(13) = rhoR_star * tau_t2t2_R_star
                                U_starR(14) = rhoR_star * tau_t1t2_R_star

                                U_starstarL(1) = U_L(1) * fac_L
                                U_starstarL(2) = U_L(2) * fac_L
                                U_starstarL(3) = rhoL_star * S_M
                                U_starstarL(4) = rhoL_star * u_t_star
                                U_starstarL(5) = rhoL_star * u_t2_star
                                U_starstarL(6) = E_L_starstar
                                U_starstarL(7) = alpha1_L_star
                                U_starstarL(8) = alpha2_L_star
                                U_starstarL(9) = rhoL_star * tau_nn_L_star
                                U_starstarL(10) = rhoL_star * tau_nt_star
                                U_starstarL(11) = rhoL_star * tau_nt2_star
                                U_starstarL(12) = rhoL_star * tau_tt_L_starstar
                                U_starstarL(13) = rhoL_star * tau_t2t2_L_starstar
                                U_starstarL(14) = rhoL_star * tau_t1t2_L_starstar

                                U_starstarR(1) = U_R(1) * fac_R
                                U_starstarR(2) = U_R(2) * fac_R
                                U_starstarR(3) = rhoR_star * S_M
                                U_starstarR(4) = rhoR_star * u_t_star
                                U_starstarR(5) = rhoR_star * u_t2_star
                                U_starstarR(6) = E_R_starstar
                                U_starstarR(7) = alpha1_R_star
                                U_starstarR(8) = alpha2_R_star
                                U_starstarR(9) = rhoR_star * tau_nn_R_star
                                U_starstarR(10) = rhoR_star * tau_nt_star
                                U_starstarR(11) = rhoR_star * tau_nt2_star
                                U_starstarR(12) = rhoR_star * tau_tt_R_starstar
                                U_starstarR(13) = rhoR_star * tau_t2t2_R_starstar
                                U_starstarR(14) = rhoR_star * tau_t1t2_R_starstar
                            else
                                ! 2D/axisym: 11-state (unchanged)
                                U_starL(1) = U_L(1) * fac_L
                                U_starL(2) = U_L(2) * fac_L
                                U_starL(3) = rhoL_star * S_M
                                U_starL(4) = rhoL_star * u_t_L
                                U_starL(5) = E_L_star
                                U_starL(6) = alpha1_L_star
                                U_starL(7) = alpha2_L_star
                                U_starL(8) = rhoL_star * tau_nn_L_star
                                U_starL(9) = rhoL_star * tau_nt_L
                                U_starL(10) = rhoL_star * tau_tt_L_star
                                U_starL(11) = rhoL_star * tau_qq_L_star

                                U_starR(1) = U_R(1) * fac_R
                                U_starR(2) = U_R(2) * fac_R
                                U_starR(3) = rhoR_star * S_M
                                U_starR(4) = rhoR_star * u_t_R
                                U_starR(5) = E_R_star
                                U_starR(6) = alpha1_R_star
                                U_starR(7) = alpha2_R_star
                                U_starR(8) = rhoR_star * tau_nn_R_star
                                U_starR(9) = rhoR_star * tau_nt_R
                                U_starR(10) = rhoR_star * tau_tt_R_star
                                U_starR(11) = rhoR_star * tau_qq_R_star

                                U_starstarL(1) = U_L(1) * fac_L
                                U_starstarL(2) = U_L(2) * fac_L
                                U_starstarL(3) = rhoL_star * S_M
                                U_starstarL(4) = rhoL_star * u_t_star
                                U_starstarL(5) = E_L_starstar
                                U_starstarL(6) = alpha1_L_star
                                U_starstarL(7) = alpha2_L_star
                                U_starstarL(8) = rhoL_star * tau_nn_L_star
                                U_starstarL(9) = rhoL_star * tau_nt_star
                                U_starstarL(10) = rhoL_star * tau_tt_L_starstar
                                U_starstarL(11) = rhoL_star * tau_qq_L_star

                                U_starstarR(1) = U_R(1) * fac_R
                                U_starstarR(2) = U_R(2) * fac_R
                                U_starstarR(3) = rhoR_star * S_M
                                U_starstarR(4) = rhoR_star * u_t_star
                                U_starstarR(5) = E_R_starstar
                                U_starstarR(6) = alpha1_R_star
                                U_starstarR(7) = alpha2_R_star
                                U_starstarR(8) = rhoR_star * tau_nn_R_star
                                U_starstarR(9) = rhoR_star * tau_nt_star
                                U_starstarR(10) = rhoR_star * tau_tt_R_starstar
                                U_starstarR(11) = rhoR_star * tau_qq_R_star
                            end if

                            ! ==================== Compute F and select F_HLLD ====================

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

                            ! ==================== ADC blending (HLLD ↔ HLL) ====================

                            if (riemann_hypo_ADC) then

                                F_HLL(1:ncomp) = F_hlld(1:ncomp)
                                if (S_L < 0._wp .and. S_R > 0._wp) then
                                    do i = 1, ncomp
                                        F_HLL(i) = ( S_R*F_L(i) - S_L*F_R(i) + S_L*S_R*(U_R(i) - U_L(i)) ) &
                                                   / (S_R - S_L + verysmall)
                                    end do
                                end if

                                Sigma_L   = pTot_L
                                Sigma_R   = pTot_R
                                dSigma    = Sigma_R - Sigma_L
                                Sigma_ref = max( max(abs(Sigma_L), abs(Sigma_R)), verysmall )

                                a_L_ref = sqrt( max( verysmall, c%L*c%L + ((4._wp/3._wp)*G_L + tau_nn_L)/rho%L ) )
                                a_R_ref = sqrt( max( verysmall, c%R*c%R + ((4._wp/3._wp)*G_R + tau_nn_R)/rho%R ) )
                                a_ref   = max( max(a_L_ref, a_R_ref), verysmall )

                                du_t    = u_t_R    - u_t_L
                                dtau_nt = tau_nt_R - tau_nt_L
                                du_t2    = u_t2_R    - u_t2_L
                                dtau_nt2 = tau_nt2_R - tau_nt2_L

                                sensor_ptot = (dSigma*dSigma)      / ( (ADC_kappa*Sigma_ref)**2 + verysmall )
                                sensor_vt   = (du_t*du_t + du_t2*du_t2) / ( (ADC_kappa*a_ref)**2 + verysmall )
                                sensor_tnt  = (dtau_nt*dtau_nt + dtau_nt2*dtau_nt2) / ( (ADC_kappa*Sigma_ref)**2 + verysmall )

                                sensor_combined = sensor_ptot + sensor_tnt + sensor_vt

                                phi = exp( - (sensor_combined**ADC_power) )

                                do i = 1, ncomp
                                    F_hlld(i) = F_HLL(i) + phi*(F_hlld(i) - F_HLL(i))
                                end do

                            end if

                            end if

                            ! ==================== Reorder F_HLLD for output ====================
                            if (p > 0 .and. .not. cyl_coord) then
                                ! 3D Cartesian: 14-state → physical indices
                                flux_rs${XYZ}$_vf(j, k, l, 1) = F_hlld(1)
                                flux_rs${XYZ}$_vf(j, k, l, 2) = F_hlld(2)
                                flux_rs${XYZ}$_vf(j, k, l, contxe + dir_idx(1)) = F_hlld(3)
                                flux_rs${XYZ}$_vf(j, k, l, contxe + dir_idx(2)) = F_hlld(4)
                                flux_rs${XYZ}$_vf(j, k, l, contxe + dir_idx(3)) = F_hlld(5)
                                flux_rs${XYZ}$_vf(j, k, l, E_idx) = F_hlld(6)
                                flux_rs${XYZ}$_vf(j, k, l, E_idx + 1) = F_hlld(7)
                                flux_rs${XYZ}$_vf(j, k, l, E_idx + 2) = F_hlld(8)
                                ! Map local stress to physical stress indices
                                if (dir_idx(1) == 1) then
                                    flux_rs${XYZ}$_vf(j, k, l, strxb    ) = F_hlld(9)  ! tau_nn=tau_xx
                                    flux_rs${XYZ}$_vf(j, k, l, strxb + 1) = F_hlld(10) ! tau_nt1=tau_xy
                                    flux_rs${XYZ}$_vf(j, k, l, strxb + 2) = F_hlld(12) ! tau_t1t1=tau_yy
                                    flux_rs${XYZ}$_vf(j, k, l, strxb + 3) = F_hlld(11) ! tau_nt2=tau_xz
                                    flux_rs${XYZ}$_vf(j, k, l, strxb + 4) = F_hlld(14) ! tau_t1t2=tau_yz
                                    flux_rs${XYZ}$_vf(j, k, l, strxb + 5) = F_hlld(13) ! tau_t2t2=tau_zz
                                else if (dir_idx(1) == 2) then
                                    flux_rs${XYZ}$_vf(j, k, l, strxb    ) = F_hlld(12) ! tau_t1t1=tau_xx
                                    flux_rs${XYZ}$_vf(j, k, l, strxb + 1) = F_hlld(10) ! tau_nt1=tau_xy
                                    flux_rs${XYZ}$_vf(j, k, l, strxb + 2) = F_hlld(9)  ! tau_nn=tau_yy
                                    flux_rs${XYZ}$_vf(j, k, l, strxb + 3) = F_hlld(14) ! tau_t1t2=tau_xz
                                    flux_rs${XYZ}$_vf(j, k, l, strxb + 4) = F_hlld(11) ! tau_nt2=tau_yz
                                    flux_rs${XYZ}$_vf(j, k, l, strxb + 5) = F_hlld(13) ! tau_t2t2=tau_zz
                                else
                                    flux_rs${XYZ}$_vf(j, k, l, strxb    ) = F_hlld(12) ! tau_t1t1=tau_xx
                                    flux_rs${XYZ}$_vf(j, k, l, strxb + 1) = F_hlld(14) ! tau_t1t2=tau_xy
                                    flux_rs${XYZ}$_vf(j, k, l, strxb + 2) = F_hlld(13) ! tau_t2t2=tau_yy
                                    flux_rs${XYZ}$_vf(j, k, l, strxb + 3) = F_hlld(10) ! tau_nt1=tau_xz
                                    flux_rs${XYZ}$_vf(j, k, l, strxb + 4) = F_hlld(11) ! tau_nt2=tau_yz
                                    flux_rs${XYZ}$_vf(j, k, l, strxb + 5) = F_hlld(9)  ! tau_nn=tau_zz
                                end if
                            else
                                ! 2D/axisym: 11-state (unchanged)
                                flux_rs${XYZ}$_vf(j, k, l, 1) = F_hlld(1)
                                flux_rs${XYZ}$_vf(j, k, l, 2) = F_hlld(2)
                                flux_rs${XYZ}$_vf(j, k, l, contxe + dir_idx(1)) = F_hlld(3)
                                flux_rs${XYZ}$_vf(j, k, l, contxe + dir_idx(2)) = F_hlld(4)
                                flux_rs${XYZ}$_vf(j, k, l, E_idx) = F_hlld(5)
                                flux_rs${XYZ}$_vf(j, k, l, E_idx + 1) = F_hlld(6)
                                flux_rs${XYZ}$_vf(j, k, l, E_idx + 2) = F_hlld(7)
                                flux_rs${XYZ}$_vf(j, k, l, strxb + 1) = F_hlld(9)
                                if (dir_idx(1) == 1) then
                                    flux_rs${XYZ}$_vf(j, k, l, strxb    ) = F_hlld(8 )
                                    flux_rs${XYZ}$_vf(j, k, l, strxb + 2) = F_hlld(10)
                                else
                                    flux_rs${XYZ}$_vf(j, k, l, strxb    ) = F_hlld(10)
                                    flux_rs${XYZ}$_vf(j, k, l, strxb + 2) = F_hlld(8 )
                                end if
                                if (cyl_coord) then
                                    flux_rs${XYZ}$_vf(j, k, l, strxe) = F_hlld(11)
                                end if
                            end if

                            ! Export face velocities for axisym hypo source terms
                            if (grid_geometry == 2) then
                                if (0._wp <= S_L) then
                                    u_n_face = u_n_L; u_t_face = u_t_L
                                elseif (0._wp <= S_Lstar) then
                                    u_n_face = S_M; u_t_face = u_t_L
                                elseif (0._wp <= S_M) then
                                    u_n_face = S_M; u_t_face = u_t_star
                                elseif (0._wp <= S_Rstar) then
                                    u_n_face = S_M; u_t_face = u_t_star
                                elseif (0._wp <= S_R) then
                                    u_n_face = S_M; u_t_face = u_t_R
                                else
                                    u_n_face = u_n_R; u_t_face = u_t_R
                                end if
                                if (dir_idx(1) == 1) then
                                    nc_iface_vel_rs${XYZ}$_vf(j, k, l, 1) = u_n_face
                                    nc_iface_vel_rs${XYZ}$_vf(j, k, l, 2) = u_t_face
                                else
                                    nc_iface_vel_rs${XYZ}$_vf(j, k, l, 1) = u_t_face
                                    nc_iface_vel_rs${XYZ}$_vf(j, k, l, 2) = u_n_face
                                end if
                            end if

                            ! Radial geometric source flux for cylindrical coordinates
                            #:if (NORM_DIR == 2)
                                if (cyl_coord) then
                                    !$acc loop seq
                                    do i = 1, sys_size
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = flux_rs${XYZ}$_vf(j, k, l, i)
                                    end do
                                    if (0._wp <= S_L) then
                                        p_face = pres%L; tau_qq_face = tau_qq_L
                                    elseif (0._wp <= S_Lstar .or. 0._wp <= S_M) then
                                        p_face = pTot_star + tau_nn_L_star; tau_qq_face = tau_qq_L_star
                                    elseif (0._wp <= S_Rstar .or. 0._wp <= S_R) then
                                        p_face = pTot_star + tau_nn_R_star; tau_qq_face = tau_qq_R_star
                                    else
                                        p_face = pres%R; tau_qq_face = tau_qq_R
                                    end if
                                    flux_gsrc_rs${XYZ}$_vf(j, k, l, contxe + dir_idx(1)) = &
                                        flux_rs${XYZ}$_vf(j, k, l, contxe + dir_idx(1)) - p_face + tau_qq_face
                                    !$acc loop seq
                                    do i = advxb, advxe
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = 0._wp
                                    end do
                                end if
                            #:endif

                            flux_src_rs${XYZ}$_vf(j, k, l, advxb) = 0._wp
                        end do
                    end do
                end do
                !$acc end parallel loop
                !$acc end data
            end if
        #:endfor

        call s_finalize_riemann_solver(flux_vf, flux_src_vf, flux_gsrc_vf, &
                                       norm_dir, ix, iy, iz)
    end subroutine s_hypo_hlld_riemann_solver

    !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_riemann_solvers_module

        ! Allocating the variables that will be utilized to formulate the
        ! left, right, and average states of the Riemann problem, as well
        ! the Riemann problem solution
        integer :: i, j

        @:ALLOCATE(Gs(1:num_fluids))

        do i = 1, num_fluids
            Gs(i) = fluid_pp(i)%G
        end do
        !$acc update device(Gs)

        if (viscous) then
            @:ALLOCATE(Res(1:2, 1:maxval(Re_size)))
        end if

        if (viscous) then
            do i = 1, 2
                do j = 1, Re_size(i)
                    Res(i, j) = fluid_pp(Re_idx(i, j))%Re(i)
                end do
            end do
            !$acc update device(Res, Re_idx, Re_size)
        end if

        !$acc enter data copyin(is1, is2, is3, isx, isy, isz)

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

        if (hypo_nc_interface .or. (hypo_nc_dual_pass .and. grid_geometry == 2) &
            .or. (riemann_solver == 1 .and. hll_alpha_interface .and. alt_soundspeed)) then
            @:ALLOCATE(nc_iface_vel_rsx_vf(is1%beg:is1%end, &
                is2%beg:is2%end, &
                is3%beg:is3%end, 1:num_dims))
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

        if (hypo_nc_interface .or. (hypo_nc_dual_pass .and. grid_geometry == 2) &
            .or. (riemann_solver == 1 .and. hll_alpha_interface .and. alt_soundspeed)) then
            @:ALLOCATE(nc_iface_vel_rsy_vf(is1%beg:is1%end, &
                is2%beg:is2%end, &
                is3%beg:is3%end, 1:num_dims))
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

        if (hypo_nc_interface .or. (hypo_nc_dual_pass .and. grid_geometry == 2) &
            .or. (riemann_solver == 1 .and. hll_alpha_interface .and. alt_soundspeed)) then
            @:ALLOCATE(nc_iface_vel_rsz_vf(is1%beg:is1%end, &
                is2%beg:is2%end, &
                is3%beg:is3%end, 1:num_dims))
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
        qL_prim_vf, &
        qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, dqR_prim_dx_vf, &
        dqR_prim_dy_vf, &
        dqR_prim_dz_vf, &
        qR_prim_vf, &
        norm_dir, ix, iy, iz)

        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:), intent(inout) :: qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf

        type(scalar_field), &
            allocatable, dimension(:), &
            intent(inout) :: dqL_prim_dx_vf, dqR_prim_dx_vf, &
                             dqL_prim_dy_vf, dqR_prim_dy_vf, &
                             dqL_prim_dz_vf, dqR_prim_dz_vf, &
                             qL_prim_vf, qR_prim_vf

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

        !$acc update device(is1, is2, is3)

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
        !$acc update device(isx, isy, isz) ! for stuff in the same module
        !$acc update device(dir_idx, dir_flg,  dir_idx_tau) ! for stuff in different modules

        ! Population of Buffers in x-direction
        if (norm_dir == 1) then

            if (bc_x%beg == BC_RIEMANN_EXTRAP) then    ! Riemann state extrap. BC at beginning
                !$acc parallel loop collapse(3) gang vector default(present)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            qL_prim_rsx_vf(-1, k, l, i) = &
                                qR_prim_rsx_vf(0, k, l, i)
                        end do
                    end do
                end do

                if (viscous) then
                    !$acc parallel loop collapse(3) gang vector default(present)
                    do i = momxb, momxe
                        do l = isz%beg, isz%end
                            do k = isy%beg, isy%end

                                dqL_prim_dx_vf(i)%sf(-1, k, l) = &
                                    dqR_prim_dx_vf(i)%sf(0, k, l)
                            end do
                        end do
                    end do

                    if (n > 0) then
                        !$acc parallel loop collapse(3) gang vector default(present)
                        do i = momxb, momxe
                            do l = isz%beg, isz%end
                                do k = isy%beg, isy%end

                                    dqL_prim_dy_vf(i)%sf(-1, k, l) = &
                                        dqR_prim_dy_vf(i)%sf(0, k, l)
                                end do
                            end do
                        end do

                        if (p > 0) then
                            !$acc parallel loop collapse(3) gang vector default(present)
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

                !$acc parallel loop collapse(3) gang vector default(present)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            qR_prim_rsx_vf(m + 1, k, l, i) = &
                                qL_prim_rsx_vf(m, k, l, i)
                        end do
                    end do
                end do

                if (viscous) then

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do i = momxb, momxe
                        do l = isz%beg, isz%end
                            do k = isy%beg, isy%end

                                dqR_prim_dx_vf(i)%sf(m + 1, k, l) = &
                                    dqL_prim_dx_vf(i)%sf(m, k, l)
                            end do
                        end do
                    end do

                    if (n > 0) then
                        !$acc parallel loop collapse(3) gang vector default(present)
                        do i = momxb, momxe
                            do l = isz%beg, isz%end
                                do k = isy%beg, isy%end

                                    dqR_prim_dy_vf(i)%sf(m + 1, k, l) = &
                                        dqL_prim_dy_vf(i)%sf(m, k, l)
                                end do
                            end do
                        end do

                        if (p > 0) then
                            !$acc parallel loop collapse(3) gang vector default(present)
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
                !$acc parallel loop collapse(3) gang vector default(present)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            qL_prim_rsy_vf(-1, k, l, i) = &
                                qR_prim_rsy_vf(0, k, l, i)
                        end do
                    end do
                end do

                if (viscous) then

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do i = momxb, momxe
                        do l = isz%beg, isz%end
                            do j = isx%beg, isx%end
                                dqL_prim_dx_vf(i)%sf(j, -1, l) = &
                                    dqR_prim_dx_vf(i)%sf(j, 0, l)
                            end do
                        end do
                    end do

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do i = momxb, momxe
                        do l = isz%beg, isz%end
                            do j = isx%beg, isx%end
                                dqL_prim_dy_vf(i)%sf(j, -1, l) = &
                                    dqR_prim_dy_vf(i)%sf(j, 0, l)
                            end do
                        end do
                    end do

                    if (p > 0) then
                        !$acc parallel loop collapse(3) gang vector default(present)
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

                !$acc parallel loop collapse(3) gang vector default(present)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            qR_prim_rsy_vf(n + 1, k, l, i) = &
                                qL_prim_rsy_vf(n, k, l, i)
                        end do
                    end do
                end do

                if (viscous) then

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do i = momxb, momxe
                        do l = isz%beg, isz%end
                            do j = isx%beg, isx%end
                                dqR_prim_dx_vf(i)%sf(j, n + 1, l) = &
                                    dqL_prim_dx_vf(i)%sf(j, n, l)
                            end do
                        end do
                    end do

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do i = momxb, momxe
                        do l = isz%beg, isz%end
                            do j = isx%beg, isx%end
                                dqR_prim_dy_vf(i)%sf(j, n + 1, l) = &
                                    dqL_prim_dy_vf(i)%sf(j, n, l)
                            end do
                        end do
                    end do

                    if (p > 0) then
                        !$acc parallel loop collapse(3) gang vector default(present)
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
                !$acc parallel loop collapse(3) gang vector default(present)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            qL_prim_rsz_vf(-1, k, l, i) = &
                                qR_prim_rsz_vf(0, k, l, i)
                        end do
                    end do
                end do

                if (viscous) then
                    !$acc parallel loop collapse(3) gang vector default(present)
                    do i = momxb, momxe
                        do k = isy%beg, isy%end
                            do j = isx%beg, isx%end
                                dqL_prim_dx_vf(i)%sf(j, k, -1) = &
                                    dqR_prim_dx_vf(i)%sf(j, k, 0)
                            end do
                        end do
                    end do
                    !$acc parallel loop collapse(3) gang vector default(present)
                    do i = momxb, momxe
                        do k = isy%beg, isy%end
                            do j = isx%beg, isx%end
                                dqL_prim_dy_vf(i)%sf(j, k, -1) = &
                                    dqR_prim_dy_vf(i)%sf(j, k, 0)
                            end do
                        end do
                    end do
                    !$acc parallel loop collapse(3) gang vector default(present)
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

                !$acc parallel loop collapse(3) gang vector default(present)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            qR_prim_rsz_vf(p + 1, k, l, i) = &
                                qL_prim_rsz_vf(p, k, l, i)
                        end do
                    end do
                end do

                if (viscous) then
                    !$acc parallel loop collapse(3) gang vector default(present)
                    do i = momxb, momxe
                        do k = isy%beg, isy%end
                            do j = isx%beg, isx%end
                                dqR_prim_dx_vf(i)%sf(j, k, p + 1) = &
                                    dqL_prim_dx_vf(i)%sf(j, k, p)
                            end do
                        end do
                    end do

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do i = momxb, momxe
                        do k = isy%beg, isy%end
                            do j = isx%beg, isx%end
                                dqR_prim_dy_vf(i)%sf(j, k, p + 1) = &
                                    dqL_prim_dy_vf(i)%sf(j, k, p)
                            end do
                        end do
                    end do

                    !$acc parallel loop collapse(3) gang vector default(present)
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
        q_prim_vf, &
        flux_vf, flux_src_vf, &
        flux_gsrc_vf, &
        norm_dir, ix, iy, iz)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: flux_vf, flux_src_vf, flux_gsrc_vf

        integer, intent(in) :: norm_dir
        type(int_bounds_info), intent(in) :: ix, iy, iz

        integer :: i, j, k, l ! Generic loop iterators

        ! Reshaping Inputted Data in x-direction

        if (norm_dir == 1) then

            if (viscous .or. (surface_tension)) then

                !$acc parallel loop collapse(4) gang vector default(present)
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

                !$acc parallel loop collapse(4) gang vector default(present)
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
                !$acc parallel loop collapse(4) gang vector default(present)
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
                !$acc parallel loop collapse(4) gang vector default(present)
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
                !$acc parallel loop collapse(4) gang vector default(present)
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
                !$acc parallel loop collapse(4) gang vector default(present)
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

    !>  The goal of this subroutine is to evaluate and account
        !!      for the contribution of viscous stresses in the source
        !!      flux for the momentum and energy.
        !!  @param velL_vf  Left, WENO reconstructed, cell-boundary values of the velocity
        !!  @param velR_vf Right, WENO reconstructed, cell-boundary values of the velocity
        !!  @param dvelL_dx_vf  Left, WENO reconstructed cell-avg. x-dir derivative of the velocity
        !!  @param dvelL_dy_vf  Left, WENO reconstructed cell-avg. y-dir derivative of the velocity
        !!  @param dvelL_dz_vf  Left, WENO reconstructed cell-avg. z-dir derivative of the velocity
        !!  @param dvelR_dx_vf Right, WENO reconstructed cell-avg. x-dir derivative of the velocity
        !!  @param dvelR_dy_vf Right, WENO reconstructed cell-avg. y-dir derivative of the velocity
        !!  @param dvelR_dz_vf Right, WENO reconstructed cell-avg. z-dir derivative of the velocity
        !!  @param flux_src_vf Intercell flux
        !!  @param norm_dir Dimensional splitting coordinate direction
        !!  @param ix Index bounds in  first coordinate direction
        !!  @param iy Index bounds in second coordinate direction
        !!  @param iz Index bounds in  third coordinate direction
    subroutine s_compute_cylindrical_viscous_source_flux(velL_vf, &
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
            dimension(num_dims), &
            intent(in) :: velL_vf, velR_vf, &
                          dvelL_dx_vf, dvelR_dx_vf, &
                          dvelL_dy_vf, dvelR_dy_vf, &
                          dvelL_dz_vf, dvelR_dz_vf

        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: flux_src_vf

        integer, intent(in) :: norm_dir

        type(int_bounds_info), intent(in) :: ix, iy, iz

        ! Arithmetic mean of the left and right, WENO-reconstructed, cell-
        ! boundary values of cell-average first-order spatial derivatives
        ! of velocity
        real(wp), dimension(num_dims) :: avg_vel
        real(wp), dimension(num_dims) :: dvel_avg_dx
        real(wp), dimension(num_dims) :: dvel_avg_dy
        real(wp), dimension(num_dims) :: dvel_avg_dz

        ! Viscous stress tensor
        real(wp), dimension(num_dims, num_dims) :: tau_Re

        ! Generic loop iterators
        integer :: i, j, k, l

        ! Viscous Stresses in z-direction
        if (norm_dir == 1) then
            if (shear_stress) then ! Shear stresses
                !$acc parallel loop collapse(3) gang vector default(present) private(avg_vel, dvel_avg_dx, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            dvel_avg_dx(1) = 5e-1_wp*(dvelL_dx_vf(1)%sf(j, k, l) &
                                                      + dvelR_dx_vf(1)%sf(j + 1, k, l))

                            tau_Re(1, 1) = (4._wp/3._wp)*dvel_avg_dx(1)/ &
                                           Re_avg_rsx_vf(j, k, l, 1)

                            flux_src_vf(momxb)%sf(j, k, l) = &
                                flux_src_vf(momxb)%sf(j, k, l) - &
                                tau_Re(1, 1)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rsx_vf(j, k, l, 1)* &
                                tau_Re(1, 1)

                        end do
                    end do
                end do
            end if

            if (bulk_stress) then ! Bulk stresses
                !$acc parallel loop collapse(3) gang vector default(present) private(avg_vel, dvel_avg_dx, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            dvel_avg_dx(1) = 5e-1_wp*(dvelL_dx_vf(1)%sf(j, k, l) &
                                                      + dvelR_dx_vf(1)%sf(j + 1, k, l))

                            tau_Re(1, 1) = dvel_avg_dx(1)/ &
                                           Re_avg_rsx_vf(j, k, l, 2)

                            flux_src_vf(momxb)%sf(j, k, l) = &
                                flux_src_vf(momxb)%sf(j, k, l) - &
                                tau_Re(1, 1)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rsx_vf(j, k, l, 1)* &
                                tau_Re(1, 1)

                        end do
                    end do
                end do
            end if

            if (n == 0) return

            if (shear_stress) then ! Shear stresses
                !$acc parallel loop collapse(3) gang vector default(present) private(avg_vel, dvel_avg_dx, dvel_avg_dy, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            avg_vel(2) = 5e-1_wp*(velL_vf(2)%sf(j, k, l) &
                                                  + velR_vf(2)%sf(j + 1, k, l))

                            !$acc loop seq
                            do i = 1, 2
                                dvel_avg_dy(i) = &
                                    5e-1_wp*(dvelL_dy_vf(i)%sf(j, k, l) &
                                             + dvelR_dy_vf(i)%sf(j + 1, k, l))
                            end do

                            dvel_avg_dx(2) = 5e-1_wp*(dvelL_dx_vf(2)%sf(j, k, l) &
                                                      + dvelR_dx_vf(2)%sf(j + 1, k, l))

                            tau_Re(1, 1) = -(2._wp/3._wp)*(dvel_avg_dy(2) + &
                                                           avg_vel(2)/y_cc(k))/ &
                                           Re_avg_rsx_vf(j, k, l, 1)

                            tau_Re(1, 2) = (dvel_avg_dy(1) + dvel_avg_dx(2))/ &
                                           Re_avg_rsx_vf(j, k, l, 1)

                            !$acc loop seq
                            do i = 1, 2
                                flux_src_vf(contxe + i)%sf(j, k, l) = &
                                    flux_src_vf(contxe + i)%sf(j, k, l) - &
                                    tau_Re(1, i)
                                flux_src_vf(E_idx)%sf(j, k, l) = &
                                    flux_src_vf(E_idx)%sf(j, k, l) - &
                                    vel_src_rsx_vf(j, k, l, i)* &
                                    tau_Re(1, i)
                            end do

                        end do
                    end do
                end do
            end if

            if (bulk_stress) then ! Bulk stresses
                !$acc parallel loop collapse(3) gang vector default(present) private(avg_vel,  dvel_avg_dy, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            avg_vel(2) = 5e-1_wp*(velL_vf(2)%sf(j, k, l) &
                                                  + velR_vf(2)%sf(j + 1, k, l))

                            dvel_avg_dy(2) = 5e-1_wp*(dvelL_dy_vf(2)%sf(j, k, l) &
                                                      + dvelR_dy_vf(2)%sf(j + 1, k, l))

                            tau_Re(1, 1) = (dvel_avg_dy(2) + &
                                            avg_vel(2)/y_cc(k))/ &
                                           Re_avg_rsx_vf(j, k, l, 2)

                            flux_src_vf(momxb)%sf(j, k, l) = &
                                flux_src_vf(momxb)%sf(j, k, l) - &
                                tau_Re(1, 1)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rsx_vf(j, k, l, 1)* &
                                tau_Re(1, 1)

                        end do
                    end do
                end do
            end if

            if (p == 0) return

            if (shear_stress) then ! Shear stresses
                !$acc parallel loop collapse(3) gang vector default(present) private(avg_vel, dvel_avg_dx, dvel_avg_dz, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            !$acc loop seq
                            do i = 1, 3, 2
                                dvel_avg_dz(i) = &
                                    5e-1_wp*(dvelL_dz_vf(i)%sf(j, k, l) &
                                             + dvelR_dz_vf(i)%sf(j + 1, k, l))
                            end do

                            dvel_avg_dx(3) = 5e-1_wp*(dvelL_dx_vf(3)%sf(j, k, l) &
                                                      + dvelR_dx_vf(3)%sf(j + 1, k, l))

                            tau_Re(1, 1) = -(2._wp/3._wp)*dvel_avg_dz(3)/y_cc(k)/ &
                                           Re_avg_rsx_vf(j, k, l, 1)

                            tau_Re(1, 3) = (dvel_avg_dz(1)/y_cc(k) + dvel_avg_dx(3))/ &
                                           Re_avg_rsx_vf(j, k, l, 1)

                            !$acc loop seq
                            do i = 1, 3, 2

                                flux_src_vf(contxe + i)%sf(j, k, l) = &
                                    flux_src_vf(contxe + i)%sf(j, k, l) - &
                                    tau_Re(1, i)

                                flux_src_vf(E_idx)%sf(j, k, l) = &
                                    flux_src_vf(E_idx)%sf(j, k, l) - &
                                    vel_src_rsx_vf(j, k, l, i)* &
                                    tau_Re(1, i)

                            end do

                        end do
                    end do
                end do
            end if

            if (bulk_stress) then ! Bulk stresses
                !$acc parallel loop collapse(3) gang vector default(present) private( avg_vel, dvel_avg_dz, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            dvel_avg_dz(3) = 5e-1_wp*(dvelL_dz_vf(3)%sf(j, k, l) &
                                                      + dvelR_dz_vf(3)%sf(j + 1, k, l))

                            tau_Re(1, 1) = dvel_avg_dz(3)/y_cc(k)/ &
                                           Re_avg_rsx_vf(j, k, l, 2)

                            flux_src_vf(momxb)%sf(j, k, l) = &
                                flux_src_vf(momxb)%sf(j, k, l) - &
                                tau_Re(1, 1)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rsx_vf(j, k, l, 1)* &
                                tau_Re(1, 1)

                        end do
                    end do
                end do
            end if
            ! END: Viscous Stresses in z-direction

            ! Viscous Stresses in r-direction
        elseif (norm_dir == 2) then

            if (shear_stress) then ! Shear stresses

                !$acc parallel loop collapse(3) gang vector default(present) private(avg_vel, dvel_avg_dx, dvel_avg_dy, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            avg_vel(2) = 5e-1_wp*(velL_vf(2)%sf(j, k, l) &
                                                  + velR_vf(2)%sf(j, k + 1, l))

                            !$acc loop seq
                            do i = 1, 2

                                dvel_avg_dx(i) = &
                                    5e-1_wp*(dvelL_dx_vf(i)%sf(j, k, l) &
                                             + dvelR_dx_vf(i)%sf(j, k + 1, l))

                                dvel_avg_dy(i) = &
                                    5e-1_wp*(dvelL_dy_vf(i)%sf(j, k, l) &
                                             + dvelR_dy_vf(i)%sf(j, k + 1, l))

                            end do

                            tau_Re(2, 1) = (dvel_avg_dy(1) + dvel_avg_dx(2))/ &
                                           Re_avg_rsy_vf(k, j, l, 1)

                            tau_Re(2, 2) = (4._wp*dvel_avg_dy(2) &
                                            - 2._wp*dvel_avg_dx(1) &
                                            - 2._wp*avg_vel(2)/y_cb(k))/ &
                                           (3._wp*Re_avg_rsy_vf(k, j, l, 1))

                            !$acc loop seq
                            do i = 1, 2

                                flux_src_vf(contxe + i)%sf(j, k, l) = &
                                    flux_src_vf(contxe + i)%sf(j, k, l) - &
                                    tau_Re(2, i)

                                flux_src_vf(E_idx)%sf(j, k, l) = &
                                    flux_src_vf(E_idx)%sf(j, k, l) - &
                                    vel_src_rsy_vf(k, j, l, i)* &
                                    tau_Re(2, i)

                            end do

                        end do
                    end do
                end do
            end if

            if (bulk_stress) then              ! Bulk stresses
                !$acc parallel loop collapse(3) gang vector default(present) private(avg_vel, dvel_avg_dx, dvel_avg_dy, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            avg_vel(2) = 5e-1_wp*(velL_vf(2)%sf(j, k, l) &
                                                  + velR_vf(2)%sf(j, k + 1, l))

                            dvel_avg_dx(1) = 5e-1_wp*(dvelL_dx_vf(1)%sf(j, k, l) &
                                                      + dvelR_dx_vf(1)%sf(j, k + 1, l))

                            dvel_avg_dy(2) = 5e-1_wp*(dvelL_dy_vf(2)%sf(j, k, l) &
                                                      + dvelR_dy_vf(2)%sf(j, k + 1, l))

                            tau_Re(2, 2) = (dvel_avg_dx(1) + dvel_avg_dy(2) + &
                                            avg_vel(2)/y_cb(k))/ &
                                           Re_avg_rsy_vf(k, j, l, 2)

                            flux_src_vf(momxb + 1)%sf(j, k, l) = &
                                flux_src_vf(momxb + 1)%sf(j, k, l) - &
                                tau_Re(2, 2)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rsy_vf(k, j, l, 2)* &
                                tau_Re(2, 2)

                        end do
                    end do
                end do
            end if

            if (p == 0) return

            if (shear_stress) then              ! Shear stresses
                !$acc parallel loop collapse(3) gang vector default(present) private(avg_vel,  dvel_avg_dy, dvel_avg_dz, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            avg_vel(3) = 5e-1_wp*(velL_vf(3)%sf(j, k, l) &
                                                  + velR_vf(3)%sf(j, k + 1, l))

                            !$acc loop seq
                            do i = 2, 3
                                dvel_avg_dz(i) = &
                                    5e-1_wp*(dvelL_dz_vf(i)%sf(j, k, l) &
                                             + dvelR_dz_vf(i)%sf(j, k + 1, l))
                            end do

                            dvel_avg_dy(3) = 5e-1_wp*(dvelL_dy_vf(3)%sf(j, k, l) &
                                                      + dvelR_dy_vf(3)%sf(j, k + 1, l))

                            tau_Re(2, 2) = -(2._wp/3._wp)*dvel_avg_dz(3)/y_cb(k)/ &
                                           Re_avg_rsy_vf(k, j, l, 1)

                            tau_Re(2, 3) = ((dvel_avg_dz(2) - avg_vel(3))/ &
                                            y_cb(k) + dvel_avg_dy(3))/ &
                                           Re_avg_rsy_vf(k, j, l, 1)

                            !$acc loop seq
                            do i = 2, 3

                                flux_src_vf(contxe + i)%sf(j, k, l) = &
                                    flux_src_vf(contxe + i)%sf(j, k, l) - &
                                    tau_Re(2, i)

                                flux_src_vf(E_idx)%sf(j, k, l) = &
                                    flux_src_vf(E_idx)%sf(j, k, l) - &
                                    vel_src_rsy_vf(k, j, l, i)* &
                                    tau_Re(2, i)

                            end do

                        end do
                    end do
                end do
            end if

            if (bulk_stress) then              ! Bulk stresses
                !$acc parallel loop collapse(3) gang vector default(present) private(avg_vel,  dvel_avg_dz, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            dvel_avg_dz(3) = 5e-1_wp*(dvelL_dz_vf(3)%sf(j, k, l) &
                                                      + dvelR_dz_vf(3)%sf(j, k + 1, l))

                            tau_Re(2, 2) = dvel_avg_dz(3)/y_cb(k)/ &
                                           Re_avg_rsy_vf(k, j, l, 2)

                            flux_src_vf(momxb + 1)%sf(j, k, l) = &
                                flux_src_vf(momxb + 1)%sf(j, k, l) - &
                                tau_Re(2, 2)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rsy_vf(k, j, l, 2)* &
                                tau_Re(2, 2)

                        end do
                    end do
                end do
            end if
            ! END: Viscous Stresses in r-direction

            ! Viscous Stresses in theta-direction
        else

            if (shear_stress) then              ! Shear stresses
                !$acc parallel loop collapse(3) gang vector default(present) private(avg_vel, dvel_avg_dx, dvel_avg_dy, dvel_avg_dz, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            !$acc loop seq
                            do i = 2, 3
                                avg_vel(i) = 5e-1_wp*(velL_vf(i)%sf(j, k, l) &
                                                      + velR_vf(i)%sf(j, k, l + 1))
                            end do

                            !$acc loop seq
                            do i = 1, 3, 2
                                dvel_avg_dx(i) = &
                                    5e-1_wp*(dvelL_dx_vf(i)%sf(j, k, l) &
                                             + dvelR_dx_vf(i)%sf(j, k, l + 1))
                            end do

                            do i = 2, 3
                                dvel_avg_dy(i) = &
                                    5e-1_wp*(dvelL_dy_vf(i)%sf(j, k, l) &
                                             + dvelR_dy_vf(i)%sf(j, k, l + 1))
                            end do

                            !$acc loop seq
                            do i = 1, 3
                                dvel_avg_dz(i) = &
                                    5e-1_wp*(dvelL_dz_vf(i)%sf(j, k, l) &
                                             + dvelR_dz_vf(i)%sf(j, k, l + 1))
                            end do

                            tau_Re(3, 1) = (dvel_avg_dz(1)/y_cc(k) + dvel_avg_dx(3))/ &
                                           Re_avg_rsz_vf(l, k, j, 1)/ &
                                           y_cc(k)

                            tau_Re(3, 2) = ((dvel_avg_dz(2) - avg_vel(3))/ &
                                            y_cc(k) + dvel_avg_dy(3))/ &
                                           Re_avg_rsz_vf(l, k, j, 1)/ &
                                           y_cc(k)

                            tau_Re(3, 3) = (4._wp*dvel_avg_dz(3)/y_cc(k) &
                                            - 2._wp*dvel_avg_dx(1) &
                                            - 2._wp*dvel_avg_dy(2) &
                                            + 4._wp*avg_vel(2)/y_cc(k))/ &
                                           (3._wp*Re_avg_rsz_vf(l, k, j, 1))/ &
                                           y_cc(k)

                            !$acc loop seq
                            do i = 1, 3
                                flux_src_vf(contxe + i)%sf(j, k, l) = &
                                    flux_src_vf(contxe + i)%sf(j, k, l) - &
                                    tau_Re(3, i)

                                flux_src_vf(E_idx)%sf(j, k, l) = &
                                    flux_src_vf(E_idx)%sf(j, k, l) - &
                                    vel_src_rsz_vf(l, k, j, i)* &
                                    tau_Re(3, i)
                            end do

                        end do
                    end do
                end do
            end if

            if (bulk_stress) then              ! Bulk stresses
                !$acc parallel loop collapse(3) gang vector default(present) private(avg_vel, dvel_avg_dx, dvel_avg_dy, dvel_avg_dz, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            avg_vel(2) = 5e-1_wp*(velL_vf(2)%sf(j, k, l) &
                                                  + velR_vf(2)%sf(j, k, l + 1))

                            dvel_avg_dx(1) = 5e-1_wp*(dvelL_dx_vf(1)%sf(j, k, l) &
                                                      + dvelR_dx_vf(1)%sf(j, k, l + 1))

                            dvel_avg_dy(2) = 5e-1_wp*(dvelL_dy_vf(2)%sf(j, k, l) &
                                                      + dvelR_dy_vf(2)%sf(j, k, l + 1))

                            dvel_avg_dz(3) = 5e-1_wp*(dvelL_dz_vf(3)%sf(j, k, l) &
                                                      + dvelR_dz_vf(3)%sf(j, k, l + 1))

                            tau_Re(3, 3) = (dvel_avg_dx(1) &
                                            + dvel_avg_dy(2) &
                                            + dvel_avg_dz(3)/y_cc(k) &
                                            + avg_vel(2)/y_cc(k))/ &
                                           Re_avg_rsz_vf(l, k, j, 2)/ &
                                           y_cc(k)

                            flux_src_vf(momxe)%sf(j, k, l) = &
                                flux_src_vf(momxe)%sf(j, k, l) - &
                                tau_Re(3, 3)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rsz_vf(l, k, j, 3)* &
                                tau_Re(3, 3)

                        end do
                    end do
                end do
            end if

        end if
        ! END: Viscous Stresses in theta-direction

    end subroutine s_compute_cylindrical_viscous_source_flux

    !>  The goal of this subroutine is to evaluate and account
        !!      for the contribution of viscous stresses in the source
        !!      flux for the momentum and energy.
        !!  @param velL_vf  Left, WENO reconstructed, cell-boundary values of the velocity
        !!  @param velR_vf Right, WENO reconstructed, cell-boundary values of the velocity
        !!  @param dvelL_dx_vf  Left, WENO reconstructed cell-avg. x-dir derivative of the velocity
        !!  @param dvelL_dy_vf  Left, WENO reconstructed cell-avg. y-dir derivative of the velocity
        !!  @param dvelL_dz_vf  Left, WENO reconstructed cell-avg. z-dir derivative of the velocity
        !!  @param dvelR_dx_vf Right, WENO reconstructed cell-avg. x-dir derivative of the velocity
        !!  @param dvelR_dy_vf Right, WENO reconstructed cell-avg. y-dir derivative of the velocity
        !!  @param dvelR_dz_vf Right, WENO reconstructed cell-avg. z-dir derivative of the velocity
        !!  @param flux_src_vf Intercell flux
        !!  @param norm_dir Dimensional splitting coordinate direction
        !!  @param ix Index bounds in  first coordinate direction
        !!  @param iy Index bounds in second coordinate direction
        !!  @param iz Index bounds in  third coordinate direction
    subroutine s_compute_cartesian_viscous_source_flux(velL_vf, &
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
            dimension(num_dims), &
            intent(in) :: velL_vf, velR_vf, &
                          dvelL_dx_vf, dvelR_dx_vf, &
                          dvelL_dy_vf, dvelR_dy_vf, &
                          dvelL_dz_vf, dvelR_dz_vf

        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: flux_src_vf

        integer, intent(in) :: norm_dir
        type(int_bounds_info), intent(in) :: ix, iy, iz

        ! Arithmetic mean of the left and right, WENO-reconstructed, cell-
        ! boundary values of cell-average first-order spatial derivatives
        ! of velocity
        real(wp), dimension(num_dims) :: dvel_avg_dx
        real(wp), dimension(num_dims) :: dvel_avg_dy
        real(wp), dimension(num_dims) :: dvel_avg_dz

        real(wp), dimension(num_dims, num_dims) :: tau_Re !< Viscous stress tensor

        integer :: i, j, k, l !< Generic loop iterators

        ! Viscous Stresses in x-direction
        if (norm_dir == 1) then

            if (shear_stress) then              ! Shear stresses
                !$acc parallel loop collapse(3) gang vector default(present) private( dvel_avg_dx, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            dvel_avg_dx(1) = 5e-1_wp*(dvelL_dx_vf(1)%sf(j, k, l) &
                                                      + dvelR_dx_vf(1)%sf(j + 1, k, l))

                            tau_Re(1, 1) = (4._wp/3._wp)*dvel_avg_dx(1)/ &
                                           Re_avg_rsx_vf(j, k, l, 1)

                            flux_src_vf(momxb)%sf(j, k, l) = &
                                flux_src_vf(momxb)%sf(j, k, l) - &
                                tau_Re(1, 1)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rsx_vf(j, k, l, 1)* &
                                tau_Re(1, 1)

                        end do
                    end do
                end do
            end if

            if (bulk_stress) then              ! Bulk stresses
                !$acc parallel loop collapse(3) gang vector default(present) private( dvel_avg_dx, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            dvel_avg_dx(1) = 5e-1_wp*(dvelL_dx_vf(1)%sf(j, k, l) &
                                                      + dvelR_dx_vf(1)%sf(j + 1, k, l))

                            tau_Re(1, 1) = dvel_avg_dx(1)/ &
                                           Re_avg_rsx_vf(j, k, l, 2)

                            flux_src_vf(momxb)%sf(j, k, l) = &
                                flux_src_vf(momxb)%sf(j, k, l) - &
                                tau_Re(1, 1)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rsx_vf(j, k, l, 1)* &
                                tau_Re(1, 1)

                        end do
                    end do
                end do
            end if

            if (n == 0) return

            if (shear_stress) then              ! Shear stresses
                !$acc parallel loop collapse(3) gang vector default(present) private(dvel_avg_dx, dvel_avg_dy, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            !$acc loop seq
                            do i = 1, 2
                                dvel_avg_dy(i) = &
                                    5e-1_wp*(dvelL_dy_vf(i)%sf(j, k, l) &
                                             + dvelR_dy_vf(i)%sf(j + 1, k, l))
                            end do

                            dvel_avg_dx(2) = 5e-1_wp*(dvelL_dx_vf(2)%sf(j, k, l) &
                                                      + dvelR_dx_vf(2)%sf(j + 1, k, l))

                            tau_Re(1, 1) = -(2._wp/3._wp)*dvel_avg_dy(2)/ &
                                           Re_avg_rsx_vf(j, k, l, 1)

                            tau_Re(1, 2) = (dvel_avg_dy(1) + dvel_avg_dx(2))/ &
                                           Re_avg_rsx_vf(j, k, l, 1)

                            !$acc loop seq
                            do i = 1, 2

                                flux_src_vf(contxe + i)%sf(j, k, l) = &
                                    flux_src_vf(contxe + i)%sf(j, k, l) - &
                                    tau_Re(1, i)

                                flux_src_vf(E_idx)%sf(j, k, l) = &
                                    flux_src_vf(E_idx)%sf(j, k, l) - &
                                    vel_src_rsx_vf(j, k, l, i)* &
                                    tau_Re(1, i)

                            end do

                        end do
                    end do
                end do
            end if

            if (bulk_stress) then              ! Bulk stresses
                !$acc parallel loop collapse(3) gang vector default(present) private( dvel_avg_dy, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            dvel_avg_dy(2) = 5e-1_wp*(dvelL_dy_vf(2)%sf(j, k, l) &
                                                      + dvelR_dy_vf(2)%sf(j + 1, k, l))

                            tau_Re(1, 1) = dvel_avg_dy(2)/ &
                                           Re_avg_rsx_vf(j, k, l, 2)

                            flux_src_vf(momxb)%sf(j, k, l) = &
                                flux_src_vf(momxb)%sf(j, k, l) - &
                                tau_Re(1, 1)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rsx_vf(j, k, l, 1)* &
                                tau_Re(1, 1)

                        end do
                    end do
                end do
            end if

            if (p == 0) return

            if (shear_stress) then              ! Shear stresses
                !$acc parallel loop collapse(3) gang vector default(present) private( dvel_avg_dx, dvel_avg_dz, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            !$acc loop seq
                            do i = 1, 3, 2
                                dvel_avg_dz(i) = &
                                    5e-1_wp*(dvelL_dz_vf(i)%sf(j, k, l) &
                                             + dvelR_dz_vf(i)%sf(j + 1, k, l))
                            end do

                            dvel_avg_dx(3) = 5e-1_wp*(dvelL_dx_vf(3)%sf(j, k, l) &
                                                      + dvelR_dx_vf(3)%sf(j + 1, k, l))

                            tau_Re(1, 1) = -(2._wp/3._wp)*dvel_avg_dz(3)/ &
                                           Re_avg_rsx_vf(j, k, l, 1)

                            tau_Re(1, 3) = (dvel_avg_dz(1) + dvel_avg_dx(3))/ &
                                           Re_avg_rsx_vf(j, k, l, 1)

                            !$acc loop seq
                            do i = 1, 3, 2
                                flux_src_vf(contxe + i)%sf(j, k, l) = &
                                    flux_src_vf(contxe + i)%sf(j, k, l) - &
                                    tau_Re(1, i)

                                flux_src_vf(E_idx)%sf(j, k, l) = &
                                    flux_src_vf(E_idx)%sf(j, k, l) - &
                                    vel_src_rsx_vf(j, k, l, i)* &
                                    tau_Re(1, i)

                            end do

                        end do
                    end do
                end do
            end if

            if (bulk_stress) then              ! Bulk stresses
                !$acc parallel loop collapse(3) gang vector default(present) private( dvel_avg_dz, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            dvel_avg_dz(3) = 5e-1_wp*(dvelL_dz_vf(3)%sf(j, k, l) &
                                                      + dvelR_dz_vf(3)%sf(j + 1, k, l))

                            tau_Re(1, 1) = dvel_avg_dz(3)/ &
                                           Re_avg_rsx_vf(j, k, l, 2)

                            flux_src_vf(momxb)%sf(j, k, l) = &
                                flux_src_vf(momxb)%sf(j, k, l) - &
                                tau_Re(1, 1)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rsx_vf(j, k, l, 1)* &
                                tau_Re(1, 1)

                        end do
                    end do
                end do
            end if
            ! END: Viscous Stresses in x-direction

            ! Viscous Stresses in y-direction
        elseif (norm_dir == 2) then

            if (shear_stress) then              ! Shear stresses
                !$acc parallel loop collapse(3) gang vector default(present) private( dvel_avg_dx, dvel_avg_dy, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            !$acc loop seq
                            do i = 1, 2

                                dvel_avg_dx(i) = &
                                    5e-1_wp*(dvelL_dx_vf(i)%sf(j, k, l) &
                                             + dvelR_dx_vf(i)%sf(j, k + 1, l))

                                dvel_avg_dy(i) = &
                                    5e-1_wp*(dvelL_dy_vf(i)%sf(j, k, l) &
                                             + dvelR_dy_vf(i)%sf(j, k + 1, l))

                            end do

                            tau_Re(2, 1) = (dvel_avg_dy(1) + dvel_avg_dx(2))/ &
                                           Re_avg_rsy_vf(k, j, l, 1)

                            tau_Re(2, 2) = (4._wp*dvel_avg_dy(2) &
                                            - 2._wp*dvel_avg_dx(1))/ &
                                           (3._wp*Re_avg_rsy_vf(k, j, l, 1))

                            !$acc loop seq
                            do i = 1, 2

                                flux_src_vf(contxe + i)%sf(j, k, l) = &
                                    flux_src_vf(contxe + i)%sf(j, k, l) - &
                                    tau_Re(2, i)

                                flux_src_vf(E_idx)%sf(j, k, l) = &
                                    flux_src_vf(E_idx)%sf(j, k, l) - &
                                    vel_src_rsy_vf(k, j, l, i)* &
                                    tau_Re(2, i)

                            end do

                        end do
                    end do
                end do
            end if

            if (bulk_stress) then              ! Bulk stresses
                !$acc parallel loop collapse(3) gang vector default(present) private( dvel_avg_dx, dvel_avg_dy, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            dvel_avg_dx(1) = 5e-1_wp*(dvelL_dx_vf(1)%sf(j, k, l) &
                                                      + dvelR_dx_vf(1)%sf(j, k + 1, l))

                            dvel_avg_dy(2) = 5e-1_wp*(dvelL_dy_vf(2)%sf(j, k, l) &
                                                      + dvelR_dy_vf(2)%sf(j, k + 1, l))

                            tau_Re(2, 2) = (dvel_avg_dx(1) + dvel_avg_dy(2))/ &
                                           Re_avg_rsy_vf(k, j, l, 2)

                            flux_src_vf(momxb + 1)%sf(j, k, l) = &
                                flux_src_vf(momxb + 1)%sf(j, k, l) - &
                                tau_Re(2, 2)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rsy_vf(k, j, l, 2)* &
                                tau_Re(2, 2)

                        end do
                    end do
                end do
            end if

            if (p == 0) return

            if (shear_stress) then              ! Shear stresses
                !$acc parallel loop collapse(3) gang vector default(present) private(  dvel_avg_dy, dvel_avg_dz, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            !$acc loop seq
                            do i = 2, 3
                                dvel_avg_dz(i) = &
                                    5e-1_wp*(dvelL_dz_vf(i)%sf(j, k, l) &
                                             + dvelR_dz_vf(i)%sf(j, k + 1, l))
                            end do

                            dvel_avg_dy(3) = 5e-1_wp*(dvelL_dy_vf(3)%sf(j, k, l) &
                                                      + dvelR_dy_vf(3)%sf(j, k + 1, l))

                            tau_Re(2, 2) = -(2._wp/3._wp)*dvel_avg_dz(3)/ &
                                           Re_avg_rsy_vf(k, j, l, 1)

                            tau_Re(2, 3) = (dvel_avg_dz(2) + dvel_avg_dy(3))/ &
                                           Re_avg_rsy_vf(k, j, l, 1)

                            !$acc loop seq
                            do i = 2, 3

                                flux_src_vf(contxe + i)%sf(j, k, l) = &
                                    flux_src_vf(contxe + i)%sf(j, k, l) - &
                                    tau_Re(2, i)

                                flux_src_vf(E_idx)%sf(j, k, l) = &
                                    flux_src_vf(E_idx)%sf(j, k, l) - &
                                    vel_src_rsy_vf(k, j, l, i)* &
                                    tau_Re(2, i)

                            end do

                        end do
                    end do
                end do
            end if

            if (bulk_stress) then              ! Bulk stresses
                !$acc parallel loop collapse(3) gang vector default(present) private( dvel_avg_dz, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            dvel_avg_dz(3) = 5e-1_wp*(dvelL_dz_vf(3)%sf(j, k, l) &
                                                      + dvelR_dz_vf(3)%sf(j, k + 1, l))

                            tau_Re(2, 2) = dvel_avg_dz(3)/ &
                                           Re_avg_rsy_vf(k, j, l, 2)

                            flux_src_vf(momxb + 1)%sf(j, k, l) = &
                                flux_src_vf(momxb + 1)%sf(j, k, l) - &
                                tau_Re(2, 2)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rsy_vf(k, j, l, 2)* &
                                tau_Re(2, 2)

                        end do
                    end do
                end do
            end if
            ! END: Viscous Stresses in y-direction

            ! Viscous Stresses in z-direction
        else

            if (shear_stress) then              ! Shear stresses
                !$acc parallel loop collapse(3) gang vector default(present) private( dvel_avg_dx, dvel_avg_dy, dvel_avg_dz, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            !$acc loop seq
                            do i = 1, 3, 2
                                dvel_avg_dx(i) = &
                                    5e-1_wp*(dvelL_dx_vf(i)%sf(j, k, l) &
                                             + dvelR_dx_vf(i)%sf(j, k, l + 1))
                            end do

                            !$acc loop seq
                            do i = 2, 3
                                dvel_avg_dy(i) = &
                                    5e-1_wp*(dvelL_dy_vf(i)%sf(j, k, l) &
                                             + dvelR_dy_vf(i)%sf(j, k, l + 1))
                            end do

                            !$acc loop seq
                            do i = 1, 3
                                dvel_avg_dz(i) = &
                                    5e-1_wp*(dvelL_dz_vf(i)%sf(j, k, l) &
                                             + dvelR_dz_vf(i)%sf(j, k, l + 1))
                            end do

                            tau_Re(3, 1) = (dvel_avg_dz(1) + dvel_avg_dx(3))/ &
                                           Re_avg_rsz_vf(l, k, j, 1)

                            tau_Re(3, 2) = (dvel_avg_dz(2) + dvel_avg_dy(3))/ &
                                           Re_avg_rsz_vf(l, k, j, 1)

                            tau_Re(3, 3) = (4._wp*dvel_avg_dz(3) &
                                            - 2._wp*dvel_avg_dx(1) &
                                            - 2._wp*dvel_avg_dy(2))/ &
                                           (3._wp*Re_avg_rsz_vf(l, k, j, 1))

                            !$acc loop seq
                            do i = 1, 3

                                flux_src_vf(contxe + i)%sf(j, k, l) = &
                                    flux_src_vf(contxe + i)%sf(j, k, l) - &
                                    tau_Re(3, i)

                                flux_src_vf(E_idx)%sf(j, k, l) = &
                                    flux_src_vf(E_idx)%sf(j, k, l) - &
                                    vel_src_rsz_vf(l, k, j, i)* &
                                    tau_Re(3, i)

                            end do

                        end do
                    end do
                end do
            end if

            if (bulk_stress) then              ! Bulk stresses
                !$acc parallel loop collapse(3) gang vector default(present) private( dvel_avg_dx, dvel_avg_dy, dvel_avg_dz, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            dvel_avg_dx(1) = 5e-1_wp*(dvelL_dx_vf(1)%sf(j, k, l) &
                                                      + dvelR_dx_vf(1)%sf(j, k, l + 1))

                            dvel_avg_dy(2) = 5e-1_wp*(dvelL_dy_vf(2)%sf(j, k, l) &
                                                      + dvelR_dy_vf(2)%sf(j, k, l + 1))

                            dvel_avg_dz(3) = 5e-1_wp*(dvelL_dz_vf(3)%sf(j, k, l) &
                                                      + dvelR_dz_vf(3)%sf(j, k, l + 1))

                            tau_Re(3, 3) = (dvel_avg_dx(1) &
                                            + dvel_avg_dy(2) &
                                            + dvel_avg_dz(3))/ &
                                           Re_avg_rsz_vf(l, k, j, 2)

                            flux_src_vf(momxe)%sf(j, k, l) = &
                                flux_src_vf(momxe)%sf(j, k, l) - &
                                tau_Re(3, 3)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rsz_vf(l, k, j, 3)* &
                                tau_Re(3, 3)

                        end do
                    end do
                end do
            end if

        end if
        ! END: Viscous Stresses in z-direction

    end subroutine s_compute_cartesian_viscous_source_flux

    !>  Deallocation and/or disassociation procedures that are
        !!      needed to finalize the selected Riemann problem solver
        !!  @param flux_vf       Intercell fluxes
        !!  @param flux_src_vf   Intercell source fluxes
        !!  @param flux_gsrc_vf  Intercell geometric source fluxes
        !!  @param norm_dir Dimensional splitting coordinate direction
        !!  @param ix   Index bounds in  first coordinate direction
        !!  @param iy   Index bounds in second coordinate direction
        !!  @param iz   Index bounds in  third coordinate direction
    subroutine s_finalize_riemann_solver(flux_vf, flux_src_vf, &
                                         flux_gsrc_vf, &
                                         norm_dir, ix, iy, iz)

        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: flux_vf, flux_src_vf, flux_gsrc_vf

        integer, intent(in) :: norm_dir
        type(int_bounds_info), intent(in) :: ix, iy, iz

        integer :: i, j, k, l !< Generic loop iterators

        ! Reshaping Outputted Data in y-direction
        if (norm_dir == 2) then
            !$acc parallel loop collapse(4) gang vector default(present)
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
                !$acc parallel loop collapse(4) gang vector default(present)
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

            !$acc parallel loop collapse(3) gang vector default(present)
            do l = is3%beg, is3%end
                do j = is1%beg, is1%end
                    do k = is2%beg, is2%end
                        flux_src_vf(advxb)%sf(k, j, l) = &
                            flux_src_rsy_vf(j, k, l, advxb)
                    end do
                end do
            end do

            if (riemann_solver == 1 .or. riemann_solver == 4) then
                !$acc parallel loop collapse(4) gang vector default(present)
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
            !$acc parallel loop collapse(4) gang vector default(present)
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
                !$acc parallel loop collapse(4) gang vector default(present)
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

            !$acc parallel loop collapse(3) gang vector default(present)
            do j = is1%beg, is1%end
                do k = is2%beg, is2%end
                    do l = is3%beg, is3%end
                        flux_src_vf(advxb)%sf(l, k, j) = &
                            flux_src_rsz_vf(j, k, l, advxb)
                    end do
                end do
            end do

            if (riemann_solver == 1 .or. riemann_solver == 4) then
                !$acc parallel loop collapse(4) gang vector default(present)
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
            !$acc parallel loop collapse(4) gang vector default(present)
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

            !$acc parallel loop collapse(3) gang vector default(present)
            do l = is3%beg, is3%end
                do k = is2%beg, is2%end
                    do j = is1%beg, is1%end
                        flux_src_vf(advxb)%sf(j, k, l) = &
                            flux_src_rsx_vf(j, k, l, advxb)
                    end do
                end do
            end do

            if ((riemann_solver == 1 .and. hll_alpha_interface) .or. riemann_solver == 4) then
                !$acc parallel loop collapse(4) gang vector default(present)
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

    !>  Copy hypo interface velocities from Riemann-space buffers to
        !!  physical-space output arrays. Called after the Riemann solver
        !!  for each sweep direction when hypo_nc_interface is active.
        !!  @param nc_iface_vel_vf Output: physical velocity components at interfaces
        !!  @param norm_dir Sweep direction (1=x, 2=y, 3=z)
    subroutine s_finalize_nc_iface_vel(nc_iface_vel_vf, norm_dir)

        type(scalar_field), dimension(:), intent(inout) :: nc_iface_vel_vf
        integer, intent(in) :: norm_dir

        integer :: i, j, k, l

        if (norm_dir == 2) then
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims
                do l = is3%beg, is3%end
                    do j = is1%beg, is1%end
                        do k = is2%beg, is2%end
                            nc_iface_vel_vf(i)%sf(k, j, l) = &
                                nc_iface_vel_rsy_vf(j, k, l, i)
                        end do
                    end do
                end do
            end do
        elseif (norm_dir == 1) then
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims
                do l = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = is1%beg, is1%end
                            nc_iface_vel_vf(i)%sf(j, k, l) = &
                                nc_iface_vel_rsx_vf(j, k, l, i)
                        end do
                    end do
                end do
            end do
        elseif (norm_dir == 3) then
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, num_dims
                do l = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = is1%beg, is1%end
                            nc_iface_vel_vf(i)%sf(l, k, j) = &
                                nc_iface_vel_rsz_vf(j, k, l, i)
                        end do
                    end do
                end do
            end do
        end if

    end subroutine s_finalize_nc_iface_vel

    !> Module deallocation and/or disassociation procedures
    subroutine s_finalize_riemann_solvers_module

        if (viscous) then
            @:DEALLOCATE(Re_avg_rsx_vf)
        end if
        @:DEALLOCATE(vel_src_rsx_vf)
        @:DEALLOCATE(flux_rsx_vf)
        @:DEALLOCATE(flux_src_rsx_vf)
        @:DEALLOCATE(flux_gsrc_rsx_vf)
        if (hypo_nc_interface .or. (hypo_nc_dual_pass .and. grid_geometry == 2) &
            .or. (riemann_solver == 1 .and. hll_alpha_interface .and. alt_soundspeed)) then
            @:DEALLOCATE(nc_iface_vel_rsx_vf)
        end if
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
        if (hypo_nc_interface .or. (hypo_nc_dual_pass .and. grid_geometry == 2) &
            .or. (riemann_solver == 1 .and. hll_alpha_interface .and. alt_soundspeed)) then
            @:DEALLOCATE(nc_iface_vel_rsy_vf)
        end if
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
        if (hypo_nc_interface .or. (hypo_nc_dual_pass .and. grid_geometry == 2) &
            .or. (riemann_solver == 1 .and. hll_alpha_interface .and. alt_soundspeed)) then
            @:DEALLOCATE(nc_iface_vel_rsz_vf)
        end if
        if (qbmm) then
            @:DEALLOCATE(mom_sp_rsz_vf)
        end if

    end subroutine s_finalize_riemann_solvers_module

end module m_riemann_solvers
