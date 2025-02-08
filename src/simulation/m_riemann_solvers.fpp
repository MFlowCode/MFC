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

        real(wp), dimension(startx:, starty:, startz:, 1:), intent(INOUT) :: qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf
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

        if (riemann_solver == 1) then
            call s_hll_riemann_solver(qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, dqL_prim_dx_vf, &
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
        elseif (riemann_solver == 2) then
            call s_hllc_riemann_solver(qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, dqL_prim_dx_vf, &
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
            dimension(num_dims), &
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

        real(wp), dimension(startx:, starty:, startz:, 1:), intent(inout) :: qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf
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
        real(wp) :: Ms_L, Ms_R, pres_SL, pres_SR
        real(wp) :: alpha_L_sum, alpha_R_sum

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
                !$acc Yi_avg, Phi_avg, h_iL, h_iR, h_avg_2)
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

                            rho_R = 0._wp
                            gamma_R = 0._wp
                            pi_inf_R = 0._wp
                            qv_R = 0._wp

                            alpha_L_sum = 0._wp
                            alpha_R_sum = 0._wp

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

                                Xs_L(:) = Ys_L(:)*MW_L/mol_weights(:)
                                Xs_R(:) = Ys_R(:)*MW_R/mol_weights(:)

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

                                do i = 1, strxe - strxb + 1
                                    tau_e_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, strxb - 1 + i)
                                    tau_e_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, strxb - 1 + i)
                                    ! Elastic contribution to energy if G large enough
                                    !TODO take out if statement if stable without
                                    if ((G_L > 1000) .and. (G_R > 1000)) then
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
                            !    if ((G_L > 1e-3_wp) .and. (G_R > 1e-3_wp)) then
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

                            ! Enthalpy with elastic energy
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

                            if (wave_speeds == 1) then
                                if (hypoelasticity) then
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
                                    s_L = min(vel_L(dir_idx(1)) - sqrt(c_L*c_L + (4_wp*G_L/3_wp)/rho_L) &
                                              , vel_R(dir_idx(1)) - sqrt(c_R*c_R + (4_wp*G_R/3_wp)/rho_R))
                                    s_R = max(vel_R(dir_idx(1)) + sqrt(c_R*c_R + (4_wp*G_R/3_wp)/rho_R) &
                                              , vel_L(dir_idx(1)) + sqrt(c_L*c_L + (4_wp*G_L/3_wp)/rho_L))
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

                            ! Mass
                            !$acc loop seq
                            do i = 1, contxe
                                flux_rs${XYZ}$_vf(j, k, l, i) = &
                                    (s_M*alpha_rho_R(i)*vel_R(dir_idx(1)) &
                                     - s_P*alpha_rho_L(i)*vel_L(dir_idx(1)) &
                                     + s_M*s_P*(alpha_rho_L(i) &
                                                - alpha_rho_R(i))) &
                                    /(s_M - s_P)
                            end do

                            ! Momentum
                            if (bubbles_euler) then
                                !$acc loop seq
                                do i = 1, num_dims
                                    flux_rs${XYZ}$_vf(j, k, l, contxe + dir_idx(i)) = &
                                        (s_M*(rho_R*vel_R(dir_idx(1)) &
                                              *vel_R(dir_idx(i)) &
                                              + dir_flg(dir_idx(i))*(pres_R - ptilde_R)) &
                                         - s_P*(rho_L*vel_L(dir_idx(1)) &
                                                *vel_L(dir_idx(i)) &
                                                + dir_flg(dir_idx(i))*(pres_L - ptilde_L)) &
                                         + s_M*s_P*(rho_L*vel_L(dir_idx(i)) &
                                                    - rho_R*vel_R(dir_idx(i)))) &
                                        /(s_M - s_P)
                                end do
                            else if (hypoelasticity) then
                                !$acc loop seq
                                do i = 1, num_dims
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
                                do i = 1, num_dims
                                    flux_rs${XYZ}$_vf(j, k, l, contxe + dir_idx(i)) = &
                                        (s_M*(rho_R*vel_R(dir_idx(1)) &
                                              *vel_R(dir_idx(i)) &
                                              + dir_flg(dir_idx(i))*pres_R) &
                                         - s_P*(rho_L*vel_L(dir_idx(1)) &
                                                *vel_L(dir_idx(i)) &
                                                + dir_flg(dir_idx(i))*pres_L) &
                                         + s_M*s_P*(rho_L*vel_L(dir_idx(i)) &
                                                    - rho_R*vel_R(dir_idx(i)))) &
                                        /(s_M - s_P)
                                end do
                            end if

                            ! Energy
                            if (bubbles_euler) then
                                flux_rs${XYZ}$_vf(j, k, l, E_idx) = &
                                    (s_M*vel_R(dir_idx(1))*(E_R + pres_R - ptilde_R) &
                                     - s_P*vel_L(dir_idx(1))*(E_L + pres_L - ptilde_L) &
                                     + s_M*s_P*(E_L - E_R)) &
                                    /(s_M - s_P)
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
                                    /(s_M - s_P)
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
                            !$acc loop seq
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
                            !$acc loop seq
                            do i = 1, num_dims
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

                                    flux_rs${XYZ}$_vf(j, k, l, i) = (s_M*Y_R*rho_R*vel_R(dir_idx(norm_dir)) &
                                                                     - s_P*Y_L*rho_L*vel_L(dir_idx(norm_dir)) &
                                                                     + s_M*s_P*(Y_L*rho_L - Y_R*rho_R)) &
                                                                    /(s_M - s_P)
                                    flux_src_rs${XYZ}$_vf(j, k, l, i) = 0._wp
                                end do
                            end if

                            #:if (NORM_DIR == 2)
                                if (cyl_coord) then
                                    !Substituting the advective flux into the inviscid geometrical source flux
                                    !$acc loop seq
                                    do i = 1, E_idx
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = flux_rs${XYZ}$_vf(j, k, l, i)
                                    end do
                                    ! Recalculating the radial momentum geometric source flux
                                    flux_gsrc_rs${XYZ}$_vf(j, k, l, contxe + dir_idx(1)) = &
                                        (s_M*(rho_R*vel_R(dir_idx(1)) &
                                              *vel_R(dir_idx(1))) &
                                         - s_P*(rho_L*vel_L(dir_idx(1)) &
                                                *vel_L(dir_idx(1))) &
                                         + s_M*s_P*(rho_L*vel_L(dir_idx(1)) &
                                                    - rho_R*vel_R(dir_idx(1)))) &
                                        /(s_M - s_P)
                                    ! Geometrical source of the void fraction(s) is zero
                                    !$acc loop seq
                                    do i = advxb, advxe
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

        real(wp), dimension(startx:, starty:, startz:, 1:), intent(inout) :: qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf
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
                    !$acc tau_e_L, tau_e_R, G_L, G_R, flux_ene_e, xi_field_L, xi_field_R)
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
                                            E_L = E_L + (tau_e_L(i)*tau_e_L(i))/(4_wp*G_L)
                                            E_R = E_R + (tau_e_R(i)*tau_e_R(i))/(4_wp*G_R)
                                            ! Additional terms in 2D and 3D
                                            if ((i == 2) .or. (i == 4) .or. (i == 5)) then
                                                E_L = E_L + (tau_e_L(i)*tau_e_L(i))/(4_wp*G_L)
                                                E_R = E_R + (tau_e_R(i)*tau_e_R(i))/(4_wp*G_R)
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

                                ! COMPUTING THE DIRECT WAVE SPEEDS
                                if (wave_speeds == 1) then
                                    if (elasticity) then
                                        s_L = min(vel_L(dir_idx(1)) - sqrt(c_L*c_L + &
                                                                           (((4_wp*G_L)/3_wp) + tau_e_L(dir_idx_tau(1)))/rho_L), vel_R(dir_idx(1)) - sqrt(c_R*c_R + &
                                                                                                                                                          (((4_wp*G_R)/3_wp) + tau_e_R(dir_idx_tau(1)))/rho_R))
                                        s_R = max(vel_R(dir_idx(1)) + sqrt(c_R*c_R + &
                                                                           (((4_wp*G_R)/3_wp) + tau_e_R(dir_idx_tau(1)))/rho_R), vel_L(dir_idx(1)) + sqrt(c_L*c_L + &
                                                                                                                                                          (((4_wp*G_L)/3_wp) + tau_e_L(dir_idx_tau(1)))/rho_L))
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
                                                                                (dir_flg(idxi)*vel_K_Star + (1_wp - dir_flg(idxi))*(xi_M*vel_L(idxi) + xi_P*vel_R(idxi))) + dir_flg(idxi)*p_Star
                                end do

                                ! ENERGY FLUX.
                                ! f = u*(E-\sigma), q = E, q_star = \xi*E+(s-u)(\rho s_star - \sigma/(s-u))
                                flux_rs${XYZ}$_vf(j, k, l, E_idx) = (E_star + p_Star)*vel_K_Star

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
                                         qvs(i))*vel_K_Star
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

                                ! SURFACE TENSION FLUX. need to check
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
                    !$acc rho_avg, h_avg, gamma_avg, alpha_L, alpha_R, s_L, s_R, s_S, vel_avg_rms, pcorr, zcoef,   &
                    !$acc vel_L_tmp, vel_R_tmp, Ys_L, Ys_R, Xs_L, Xs_R, Gamma_iL, Gamma_iR, Cp_iL, Cp_iR,          &
                    !$acc tau_e_L, tau_e_R, xi_field_L, xi_field_R,                                                &
                    !$acc Yi_avg, Phi_avg, h_iL, h_iR, h_avg_2) copyin(is1,is2,is3)
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

                                    Xs_L(:) = Ys_L(:)*MW_L/mol_weights(:)
                                    Xs_R(:) = Ys_R(:)*MW_R/mol_weights(:)

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
                                            E_L = E_L + (tau_e_L(i)*tau_e_L(i))/(4_wp*G_L)
                                            E_R = E_R + (tau_e_R(i)*tau_e_R(i))/(4_wp*G_R)
                                            ! Additional terms in 2D and 3D
                                            if ((i == 2) .or. (i == 4) .or. (i == 5)) then
                                                E_L = E_L + (tau_e_L(i)*tau_e_L(i))/(4_wp*G_L)
                                                E_R = E_R + (tau_e_R(i)*tau_e_R(i))/(4_wp*G_R)
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

                                if (low_Mach == 2) then
                                    @:compute_low_Mach_correction()
                                end if

                                if (wave_speeds == 1) then
                                    if (elasticity) then
                                        s_L = min(vel_L(dir_idx(1)) - sqrt(c_L*c_L + &
                                                                           (((4_wp*G_L)/3_wp) + tau_e_L(dir_idx_tau(1)))/rho_L), vel_R(dir_idx(1)) - sqrt(c_R*c_R + &
                                                                                                                                                          (((4_wp*G_R)/3_wp) + tau_e_R(dir_idx_tau(1)))/rho_R))
                                        s_R = max(vel_R(dir_idx(1)) + sqrt(c_R*c_R + &
                                                                           (((4_wp*G_R)/3_wp) + tau_e_R(dir_idx_tau(1)))/rho_R), vel_L(dir_idx(1)) + sqrt(c_L*c_L + &
                                                                                                                                                          (((4_wp*G_L)/3_wp) + tau_e_L(dir_idx_tau(1)))/rho_L))
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

                                ! COMPUTING THE HLLC FLUXES
                                ! MASS FLUX.
                                if (low_Mach == 1) then
                                    @:compute_low_Mach_correction()
                                else
                                    pcorr = 0._wp
                                end if

                                !$acc loop seq
                                do i = 1, contxe
                                    flux_rs${XYZ}$_vf(j, k, l, i) = &
                                        xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, i) &
                                        *(vel_L(idx1) + s_M*(xi_L - 1._wp)) &
                                        + xi_P*qR_prim_rs${XYZ}$_vf(j + 1, k, l, i) &
                                        *(vel_R(idx1) + s_P*(xi_R - 1._wp))
                                end do

                                ! MOMENTUM FLUX.
                                ! f = \rho u u - \sigma, q = \rho u, q_star = \xi * \rho*(s_star, v, w)
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
                                    flux_ene_e = 0_wp
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

                                ! HYPOELASTIC STRESS EVOLUTION FLUX.
                                if (hypoelasticity) then
                                    !$acc loop seq
                                    do i = 1, strxe - strxb + 1
                                        flux_rs${XYZ}$_vf(j, k, l, strxb - 1 + i) = &
                                            xi_M*(s_S/(s_L - s_S))*(s_L*rho_L*tau_e_L(i) - rho_L*vel_L(idx1)*tau_e_L(i)) + &
                                            xi_P*(s_S/(s_R - s_S))*(s_R*rho_R*tau_e_R(i) - rho_R*vel_R(idx1)*tau_e_R(i))
                                    end do
                                end if

                                ! VOLUME FRACTION FLUX.
                                !$acc loop seq
                                do i = advxb, advxe
                                    flux_rs${XYZ}$_vf(j, k, l, i) = &
                                        xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, i) &
                                        *(vel_L(idx1) + s_M*(xi_L - 1._wp)) &
                                        + xi_P*qR_prim_rs${XYZ}$_vf(j + 1, k, l, i) &
                                        *(vel_R(idx1) + s_P*(xi_R - 1._wp))
                                end do

                                ! VOLUME FRACTION SOURCE FLUX.
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

                                ! Geometrical source flux for cylindrical coordinates
                                #:if (NORM_DIR == 2)
                                    if (cyl_coord) then
                                        !Substituting the advective flux into the inviscid geometrical source flux
                                        !$acc loop seq
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
            is3%beg:is3%end, 1:num_dims))
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
            is3%beg:is3%end, 1:num_dims))

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
            is3%beg:is3%end, 1:num_dims))

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
        qL_prim_vf, &
        qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, dqR_prim_dx_vf, &
        dqR_prim_dy_vf, &
        dqR_prim_dz_vf, &
        qR_prim_vf, &
        norm_dir, ix, iy, iz)

        real(wp), dimension(startx:, starty:, startz:, 1:), intent(inout) :: qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf

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

            if (bc_x%beg == -4) then    ! Riemann state extrap. BC at beginning
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

            if (bc_x%end == -4) then    ! Riemann state extrap. BC at end

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

            if (bc_y%beg == -4) then    ! Riemann state extrap. BC at beginning
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

            if (bc_y%end == -4) then    ! Riemann state extrap. BC at end

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

            if (bc_z%beg == -4) then    ! Riemann state extrap. BC at beginning
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

            if (bc_z%end == -4) then    ! Riemann state extrap. BC at end

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

            if (riemann_solver == 1) then
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

            if (riemann_solver == 1) then
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

            if (riemann_solver == 1) then
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

    !> Module deallocation and/or disassociation procedures
    subroutine s_finalize_riemann_solvers_module

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
