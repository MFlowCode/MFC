!>
!! @file m_riemann_solvers.f90
!! @brief Contains module m_riemann_solvers
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

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
module m_riemann_solvers

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    use m_bubbles              !< To get the bubble wall pressure function
    ! ==========================================================================

    implicit none

    private; public :: s_initialize_riemann_solvers_module, &
 s_riemann_solver, &
 s_hll_riemann_solver, &
 s_hll_riemann_solver_acc, &
 s_hllc_riemann_solver, &
 s_hllc_riemann_solver_acc, &
 s_exact_riemann_solver, &
 s_finalize_riemann_solvers_module

    abstract interface ! =======================================================

        !> Abstract interface to the subroutines that are utilized to compute the
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
        subroutine s_abstract_riemann_solver(qL_prim_rsx_vf_flat, qL_prim_rsy_vf_flat, qL_prim_rsz_vf_flat, dqL_prim_dx_vf, &
                                             dqL_prim_dy_vf, &
                                             dqL_prim_dz_vf, &
                                             gm_alphaL_vf, &
                                             qR_prim_rsx_vf_flat, qR_prim_rsy_vf_flat, qR_prim_rsz_vf_flat, dqR_prim_dx_vf, &
                                             dqR_prim_dy_vf, &
                                             dqR_prim_dz_vf, &
                                             gm_alphaR_vf, &
                                             q_prim_vf, &
                                             flux_vf, flux_src_vf, &
                                             flux_gsrc_vf, &
                                             norm_dir, ix, iy, iz)

            import :: scalar_field, bounds_info, sys_size

            real(kind(0d0)), dimension(startx:, starty:, startz:, 1:), intent(INOUT) :: qL_prim_rsx_vf_flat, qL_prim_rsy_vf_flat, qL_prim_rsz_vf_flat, qR_prim_rsx_vf_flat, qR_prim_rsy_vf_flat, qR_prim_rsz_vf_flat
            type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf

            type(scalar_field), &
                allocatable, dimension(:), &
                intent(INOUT) :: dqL_prim_dx_vf, dqR_prim_dx_vf, &
                                 dqL_prim_dy_vf, dqR_prim_dy_vf, &
                                 dqL_prim_dz_vf, dqR_prim_dz_vf, &
                                 gm_alphaL_vf, gm_alphaR_vf

            type(scalar_field), &
                dimension(sys_size), &
                intent(INOUT) :: flux_vf, flux_src_vf, flux_gsrc_vf

            integer, intent(IN) :: norm_dir

            type(bounds_info), intent(IN) :: ix, iy, iz

        end subroutine s_abstract_riemann_solver

        !>  The abstract interface to the subroutines that are used to calculate
        !!  the Roe and arithmetic average states. For more information refer to:
        !!      1) s_compute_roe_average_state
        !!      2) s_compute_arithmetic_average_state
        !!  @param i First coordinate location index
        !!  @param j Second coordinate location index
        !!  @param k Third coordinate location index
        subroutine s_compute_abstract_average_state(qL_prim_rs_vf, qR_prim_rs_vf,i, j, k)
            import :: scalar_field, bounds_info, sys_size
            integer, intent(IN) :: i, j, k
            type(scalar_field), dimension(sys_size), intent(IN) :: qL_prim_rs_vf, qR_prim_rs_vf

        end subroutine s_compute_abstract_average_state


        !> The abstract interface to the subroutines that are utilized to compute
        !! the wave speeds of the Riemann problem either directly or by the means
        !! of pressure-velocity estimates. For more information please refer to:
        !!      1) s_compute_direct_wave_speeds
        !!      2) s_compute_pressure_velocity_wave_speeds
        !!  @param i First coordinate location index
        !!  @param j Second coordinate location index
        !!  @param k Third coordinate location index
        subroutine s_compute_abstract_wave_speeds(i, j, k)

            integer, intent(IN) :: i, j, k

        end subroutine s_compute_abstract_wave_speeds

        !> The abstract interface to the subroutines that are utilized to compute
        !! the viscous source fluxes for either Cartesian or cylindrical geometries.
        !! For more information please refer to:
        !!      1) s_compute_cartesian_viscous_source_flux
        !!      2) s_compute_cylindrical_viscous_source_flux
        subroutine s_compute_abstract_viscous_source_flux(velL_vf, & ! -------------
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

            import :: scalar_field, bounds_info, num_dims, sys_size

            type(scalar_field), &
                dimension(num_dims), &
                intent(IN) ::         velL_vf, velR_vf, &
                              dvelL_dx_vf, dvelR_dx_vf, &
                              dvelL_dy_vf, dvelR_dy_vf, &
                              dvelL_dz_vf, dvelR_dz_vf

            type(scalar_field), &
                dimension(sys_size), &
                intent(INOUT) :: flux_src_vf

            integer, intent(IN) :: norm_dir

            type(bounds_info), intent(IN) :: ix, iy, iz

        end subroutine s_compute_abstract_viscous_source_flux

    end interface ! ============================================================

    type(scalar_field), allocatable, dimension(:) :: qL_prim_rs_vf
    type(scalar_field), allocatable, dimension(:) :: qR_prim_rs_vf
    type(scalar_field), allocatable, dimension(:) :: flux_rs_vf, flux_src_rs_vf
    type(scalar_field), allocatable, dimension(:) :: flux_gsrc_rs_vf !<
    type(scalar_field), allocatable, dimension(:) :: vel_src_rs_vf


    !> The left (L) and right (R) WENO-reconstructed cell-boundary values of the
    !! cell-average primitive variables that define the left and right states of
    !! the Riemann problem. Variables qK_prim_rs_vf, K = L or R, are obtained by
    !! reshaping (RS) qK_prim_vf in a coordinate direction that is normal to the
    !! cell-boundaries along which the fluxes are to be determined.
    !> @{
    type(scalar_field), allocatable, dimension(:) :: qL_prim_rsx_vf
    type(scalar_field), allocatable, dimension(:) :: qR_prim_rsx_vf

    type(scalar_field), allocatable, dimension(:) :: qL_prim_rsy_vf
    type(scalar_field), allocatable, dimension(:) :: qR_prim_rsy_vf

    type(scalar_field), allocatable, dimension(:) :: qL_prim_rsz_vf
    type(scalar_field), allocatable, dimension(:) :: qR_prim_rsz_vf

    !> @}


    !> The cell-boundary values of the fluxes (src - source) that are computed
    !! through the chosen Riemann problem solver, and the direct evaluation of
    !! source terms, by using the left and right states given in qK_prim_rs_vf,
    !! dqK_prim_ds_vf and kappaK_rs_vf, where ds = dx, dy or dz.
    !> @{
    type(scalar_field), allocatable, dimension(:) :: flux_rsx_vf, flux_src_rsx_vf
   type(scalar_field), allocatable, dimension(:) :: flux_rsy_vf, flux_src_rsy_vf
   type(scalar_field), allocatable, dimension(:) :: flux_rsz_vf, flux_src_rsz_vf

    !> @}

    type(scalar_field), allocatable, dimension(:) :: flux_gsrc_rsx_vf !<
   type(scalar_field), allocatable, dimension(:) :: flux_gsrc_rsy_vf !<
   type(scalar_field), allocatable, dimension(:) :: flux_gsrc_rsz_vf !<

    !! The cell-boundary values of the geometrical source flux that are computed
    !! through the chosen Riemann problem solver by using the left and right
    !! states given in qK_prim_rs_vf. Currently 2D axisymmetric for inviscid only.

    ! The cell-boundary values of the velocity. vel_src_rs_vf is determined as
    ! part of Riemann problem solution and is used to evaluate the source flux.
   type(scalar_field), allocatable, dimension(:) :: vel_src_rsx_vf
   type(scalar_field), allocatable, dimension(:):: vel_src_rsy_vf
   type(scalar_field), allocatable, dimension(:) :: vel_src_rsz_vf



    !> @}


    !> The cell-boundary values of the fluxes (src - source) that are computed
    !! through the chosen Riemann problem solver, and the direct evaluation of
    !! source terms, by using the left and right states given in qK_prim_rs_vf,
    !! dqK_prim_ds_vf and kappaK_rs_vf, where ds = dx, dy or dz.
    !> @{
    real(kind(0d0)), allocatable, dimension(:,:,:,:) :: flux_rsx_vf_flat, flux_src_rsx_vf_flat
   real(kind(0d0)), allocatable, dimension(:,:,:,:) :: flux_rsy_vf_flat, flux_src_rsy_vf_flat
   real(kind(0d0)), allocatable, dimension(:,:,:,:) :: flux_rsz_vf_flat, flux_src_rsz_vf_flat

    !> @}

    real(kind(0d0)), allocatable, dimension(:,:,:,:) :: flux_gsrc_rsx_vf_flat !<
   real(kind(0d0)), allocatable, dimension(:,:,:,:) :: flux_gsrc_rsy_vf_flat !<
   real(kind(0d0)), allocatable, dimension(:,:,:,:) :: flux_gsrc_rsz_vf_flat !<

    !! The cell-boundary values of the geometrical source flux that are computed
    !! through the chosen Riemann problem solver by using the left and right
    !! states given in qK_prim_rs_vf. Currently 2D axisymmetric for inviscid only.

    ! The cell-boundary values of the velocity. vel_src_rs_vf is determined as
    ! part of Riemann problem solution and is used to evaluate the source flux.
    real(kind(0d0)), allocatable, dimension(:,:,:,:) :: vel_src_rsx_vf_flat
   real(kind(0d0)), allocatable, dimension(:,:,:,:) :: vel_src_rsy_vf_flat
   real(kind(0d0)), allocatable, dimension(:,:,:,:) :: vel_src_rsz_vf_flat


    !> @name Left and right, WENO-reconstructed, cell-boundary values of cell-average
    !! partial densities, density, velocity, pressure, internal energy, energy, enthalpy, volume
    !! fractions, mass fractions, the specific heat ratio and liquid stiffness functions, speed
    !! of sound, shear and volume Reynolds numbers and the Weber numbers. These
    !! variables are left and right states of the Riemann problem obtained from
    !! qK_prim_rs_vf and kappaK_rs_vf.
    !> @{
    real(kind(0d0)), allocatable, dimension(:)   :: alpha_rho_L, alpha_rho_R
    real(kind(0d0))                              ::       rho_L, rho_R
    real(kind(0d0)), allocatable, dimension(:)   ::       vel_L, vel_R
    real(kind(0d0))                              ::      pres_L, pres_R
    real(kind(0d0))                              ::         E_L, E_R
    real(kind(0d0))                              ::         H_L, H_R
    real(kind(0d0)), allocatable, dimension(:)   ::     alpha_L, alpha_R
    real(kind(0d0))                              ::         Y_L, Y_R
    real(kind(0d0))                              ::     gamma_L, gamma_R
    real(kind(0d0))                              ::    pi_inf_L, pi_inf_R
    real(kind(0d0))                              ::         c_L, c_R
    real(kind(0d0)), dimension(2)   ::        Re_L, Re_R
    real(kind(0d0)), allocatable, dimension(:)   ::     tau_e_L, tau_e_R

!$acc declare create(alpha_rho_L, alpha_rho_R,rho_L, rho_R,vel_L, vel_R,pres_L, pres_R, &
!$acc    E_L, E_R, H_L, H_R, alpha_L, alpha_R, Y_L, Y_R, gamma_L, gamma_R,pi_inf_L, pi_inf_R, &
!$acc    c_L, c_R,Re_L, Re_R,tau_e_L, tau_e_R)

    !> @}

    !> @name Left and right, WENO-reconstructed, cell-boundary values of cell-average
    !! bubble density, radius, radial velocity, pressure, wall pressure, and modified
    !! pressure. These variables are left and right states of the Riemann problem obtained from
    !! qK_prim_rs_vf and kappaK_rs_vf.
    !> @{
    real(kind(0d0))                              ::       nbub_L, nbub_R
    real(kind(0d0)), allocatable, dimension(:)   ::         R0_L, R0_R
    real(kind(0d0)), allocatable, dimension(:)   ::         V0_L, V0_R
    real(kind(0d0)), allocatable, dimension(:)   ::         P0_L, P0_R
    real(kind(0d0)), allocatable, dimension(:)   ::        pbw_L, pbw_R
    real(kind(0d0)), allocatable, dimension(:, :) ::       moms_L, moms_R
    real(kind(0d0))                              ::     ptilde_L, ptilde_R
    !> @}
!$acc declare create(nbub_L, nbub_R, R0_L, R0_R, V0_L, V0_R, P0_L, P0_R, pbw_L, pbw_R, moms_L, moms_R, ptilde_L, ptilde_R )

    !> @name Gamma-related constants for use in exact Riemann solver (following Toro (1999) pp.153)
    !> @{
    real(kind(0d0)) :: G1_L, G1_R
    real(kind(0d0)) :: G2_L, G2_R
    real(kind(0d0)) :: G3_L, G3_R
    real(kind(0d0)) :: G4_L, G4_R
    real(kind(0d0)) :: G5_L, G5_R
    real(kind(0d0)) :: G6_L, G6_R
    real(kind(0d0)) :: G7_L, G7_R
    real(kind(0d0)) :: G8_L, G8_R
    !> @}

    !> @name Star region pressure and velocity
    !> @{
    real(kind(0d0)) :: pres_S
    real(kind(0d0)) :: vel_S
    !> @}

    !> @name Intercell solution used to calculated intercell flux
    !> @{
    real(kind(0d0)), allocatable, dimension(:)   :: alpha_rho_IC
    real(kind(0d0))                              :: rho_IC
    real(kind(0d0)), allocatable, dimension(:)   :: vel_IC
    real(kind(0d0))                              :: pres_IC
    real(kind(0d0))                              :: E_IC
    real(kind(0d0)), allocatable, dimension(:)   :: alpha_IC
    real(kind(0d0)), allocatable, dimension(:)   :: tau_e_IC
    !> @}

    !> @name Surface tension pressure contribution
    !> @{
    real(kind(0d0)) :: dpres_L, dpres_R
    !> @}
!$acc declare create(pres_S, vel_S, alpha_IC, alpha_rho_IC, vel_IC, pres_IC, E_IC, rho_IC, tau_e_IC, dpres_L, dpres_R)

    !> @name Roe or arithmetic average density, velocity, enthalpy, volume fractions,
    !! specific heat ratio function, speed of sound, shear and volume Reynolds
    !! numbers, Weber numbers and curvatures, at the cell-boundaries, computed
    !! from the left and the right states of the Riemann problem
    !> @{
    real(kind(0d0))                                 :: rho_avg
    real(kind(0d0)), allocatable, dimension(:)   :: vel_avg
    real(kind(0d0))                                 :: H_avg
    type(scalar_field), allocatable, dimension(:)   :: alpha_avg_rs_vf
    real(kind(0d0))                                 :: gamma_avg
    real(kind(0d0))                                 :: c_avg
    type(scalar_field), allocatable, dimension(:)   :: Re_avg_rs_vf
    type(scalar_field), allocatable, dimension(:)   :: Re_avg_rsx_vf
    type(scalar_field), allocatable, dimension(:)   :: Re_avg_rsy_vf
    type(scalar_field), allocatable, dimension(:)   :: Re_avg_rsz_vf
    real(kind(0d0)), allocatable, dimension(:,:,:,:)   :: Re_avg_rsx_vf_flat
    real(kind(0d0)), allocatable, dimension(:,:,:,:)   :: Re_avg_rsy_vf_flat
    real(kind(0d0)), allocatable, dimension(:,:,:,:)   :: Re_avg_rsz_vf_flat
!$acc declare create(rho_avg, vel_avg, H_avg, alpha_avg_rs_vf, gamma_avg, c_avg, Re_avg_rs_vf, Re_avg_rsx_vf, Re_avg_rsy_vf, Re_avg_rsz_vf, Re_avg_rsx_vf_flat, Re_avg_rsy_vf_flat, Re_avg_rsz_vf_flat)
    !> @}

    !> @name Left, right and star (S) region wave speeds
    !> @{
    real(kind(0d0)) :: s_L, s_R, s_S
    !> @}

    !> @name Star region variables (HLLC)
    !> @{
    real(kind(0d0)) :: rho_Star, E_Star, p_Star, p_K_Star
    !> @}

    !> Minus (M) and plus (P) wave speeds
    !> @{
    real(kind(0d0)) :: s_M, s_P
    !> @}

    !> Minus and plus wave speeds functions
    !> @{
    real(kind(0d0)) :: xi_M, xi_P
    !> @}
    real(kind(0d0)) :: xi_L, xi_R

!$acc declare create(s_L, s_R, s_S, rho_Star, E_Star, p_Star, p_K_Star, s_M, s_P, xi_M, xi_P, xi_L, xi_R)

    procedure(s_abstract_riemann_solver), &
        pointer :: s_riemann_solver => null() !<
    !! Pointer to the procedure that is utilized to calculate either the HLL,
    !! HLLC or exact intercell fluxes, based on the choice of Riemann solver

    procedure(s_compute_abstract_average_state), &
        pointer :: s_compute_average_state => null() !<
    !! Pointer to the subroutine utilized to calculate either the Roe or the
    !! arithmetic average state variables, based on the chosen average state

    procedure(s_compute_abstract_wave_speeds), &
        pointer :: s_compute_wave_speeds => null() !<
    !! Pointer to the subroutine that is utilized to compute the wave speeds of
    !! the Riemann problem either directly or by the means of pressure-velocity
    !! estimates, based on the selected method of estimation of the wave speeds

    procedure(s_compute_abstract_viscous_source_flux), &
        pointer :: s_compute_viscous_source_flux => null() !<
    !! Pointer to the subroutine that is utilized to compute the viscous source
    !! flux for either Cartesian or cylindrical geometries.



    !> @name Indical bounds in the s1-, s2- and s3-directions
    !> @{
    type(bounds_info) :: is1, is2, is3
    !> @}
!$acc declare create(qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, &
!$acc    is1, is2, is3, flux_rsx_vf, flux_src_rsx_vf, flux_rsy_vf, flux_src_rsy_vf, flux_rsz_vf, flux_src_rsz_vf, vel_src_rsx_vf, vel_src_rsy_vf, vel_src_rsz_vf, &
!$acc    flux_gsrc_rsx_vf, flux_gsrc_rsy_vf, flux_gsrc_rsz_vf)

!$acc declare create(&
!$acc    flux_rsx_vf_flat, flux_src_rsx_vf_flat, flux_rsy_vf_flat, flux_src_rsy_vf_flat, flux_rsz_vf_flat, flux_src_rsz_vf_flat, vel_src_rsx_vf_flat, vel_src_rsy_vf_flat, vel_src_rsz_vf_flat, &
!$acc    flux_gsrc_rsx_vf_flat, flux_gsrc_rsy_vf_flat, flux_gsrc_rsz_vf_flat)

    real(kind(0d0)) :: momxb, momxe
    real(kind(0d0)) :: contxb, contxe
    real(kind(0d0)) :: advxb, advxe
    real(kind(0d0)) :: bubxb, bubxe
    real(kind(0d0)) :: intxb, intxe

!$acc declare create(momxb, momxe, contxb, contxe, advxb, advxe, bubxb, bubxe, intxb, intxe)

    real(kind(0d0)),allocatable, dimension(:) :: gammas, pi_infs
!$acc declare create(gammas, pi_infs)

    real(kind(0d0)),allocatable, dimension(:) :: rs, vs, ps, ms
!$acc declare create(rs, vs, ps, ms)
contains




    subroutine s_hll_riemann_solver_acc(qL_prim_rsx_vf_flat, qL_prim_rsy_vf_flat, qL_prim_rsz_vf_flat, dqL_prim_dx_vf, & ! -------
                                    dqL_prim_dy_vf, &
                                    dqL_prim_dz_vf, &
                                    gm_alphaL_vf, &
                                    qR_prim_rsx_vf_flat, qR_prim_rsy_vf_flat, qR_prim_rsz_vf_flat, dqR_prim_dx_vf, &
                                    dqR_prim_dy_vf, &
                                    dqR_prim_dz_vf, &
                                    gm_alphaR_vf, &
                                    q_prim_vf, &
                                    flux_vf, flux_src_vf, &
                                    flux_gsrc_vf, &
                                    norm_dir, ix, iy, iz)

        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:), intent(INOUT) :: qL_prim_rsx_vf_flat, qL_prim_rsy_vf_flat, qL_prim_rsz_vf_flat, qR_prim_rsx_vf_flat, qR_prim_rsy_vf_flat, qR_prim_rsz_vf_flat
        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf

        type(scalar_field), &
            allocatable, dimension(:), &
            intent(INOUT) :: dqL_prim_dx_vf, dqR_prim_dx_vf, &
                             dqL_prim_dy_vf, dqR_prim_dy_vf, &
                             dqL_prim_dz_vf, dqR_prim_dz_vf, &
                             gm_alphaL_vf, gm_alphaR_vf

        ! Intercell fluxes
        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: flux_vf, flux_src_vf, flux_gsrc_vf

        integer, intent(IN) :: norm_dir
        type(bounds_info), intent(IN) :: ix, iy, iz

        real(kind(0d0)),dimension(10)   :: alpha_rho_L_acc, alpha_rho_R_acc
        real(kind(0d0))                              ::       rho_L_acc, rho_R_acc
        real(kind(0d0)), dimension(3)   ::       vel_L_acc, vel_R_acc
        real(kind(0d0))                              ::      pres_L_acc, pres_R_acc
        real(kind(0d0))                              ::         E_L_acc, E_R_acc
        real(kind(0d0))                              ::         H_L_acc, H_R_acc
        real(kind(0d0)), dimension(10)   ::     alpha_L_acc, alpha_R_acc
        real(kind(0d0))                              ::         Y_L_acc, Y_R_acc
        real(kind(0d0))                              ::     gamma_L_acc, gamma_R_acc
        real(kind(0d0))                              ::    pi_inf_L_acc, pi_inf_R_acc
        real(kind(0d0))                              ::         c_L_acc, c_R_acc

        real(kind(0d0))                                 :: rho_avg_acc
        real(kind(0d0)),dimension(3)   :: vel_avg_acc
        real(kind(0d0))                                 :: H_avg_acc
        real(kind(0d0))                                 :: gamma_avg_acc
        real(kind(0d0))                                 :: c_avg_acc

        real(kind(0d0))     :: s_L_acc, s_R_acc, s_M_acc, s_P_acc, s_S_acc
        real(kind(0d0)) :: xi_L_acc, xi_R_acc !< Left and right wave speeds functions
        real(kind(0d0)) :: xi_M_acc, xi_P_acc

        real(kind(0d0))                              ::       nbub_L_acc, nbub_R_acc
        real(kind(0d0))                              ::     ptilde_L_acc, ptilde_R_acc
        real(kind(0d0))  :: vel_L_rms_acc, vel_R_rms_acc, vel_avg_rms_acc
        real(kind(0d0)) :: blkmod1, blkmod2
        real(kind(0d0)) :: rho_Star_acc, E_Star_acc, p_Star_acc, p_K_Star_acc
        real(kind(0d0)) :: Ms_L, Ms_R, pres_SL, pres_SR

        integer :: i, j, k, l !< Generic loop iterators

        ! Populating the buffers of the left and right Riemann problem
        ! states variables, based on the choice of boundary conditions
        call s_populate_riemann_states_variables_buffers( &
            qL_prim_rsx_vf_flat, qL_prim_rsy_vf_flat, qL_prim_rsz_vf_flat, dqL_prim_dx_vf, &
            dqL_prim_dy_vf, &
            dqL_prim_dz_vf, &
            gm_alphaL_vf, &
            qR_prim_rsx_vf_flat, qR_prim_rsy_vf_flat, qR_prim_rsz_vf_flat, dqR_prim_dx_vf, &
            dqR_prim_dy_vf, &
            dqR_prim_dz_vf, &
            gm_alphaR_vf, &
            norm_dir, ix, iy, iz)

        ! Reshaping inputted data based on dimensional splitting direction
        call s_initialize_riemann_solver(&
                                         q_prim_vf, &
                                         flux_vf, flux_src_vf, &
                                         flux_gsrc_vf, &
                                             norm_dir, ix, iy, iz)
#:for NORM_DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]

        if (norm_dir == ${NORM_DIR}$) then
            !$acc parallel loop collapse(3) gang vector default(present) private(alpha_rho_L_acc, alpha_rho_R_acc, vel_L_acc, vel_R_acc, alpha_L_acc, alpha_R_acc, vel_avg_acc)
            do l = is3%beg, is3%end
              do k = is2%beg, is2%end
                do j = is1%beg, is1%end
                  !$acc loop seq
                  do i = 1, contxe
                    alpha_rho_L_acc(i) = qL_prim_rs${XYZ}$_vf_flat(j,     k, l, i)
                    alpha_rho_R_acc(i) = qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, i)
                  end do

                  !$acc loop seq
                  do i = 1, num_dims
                    vel_L_acc(i) = qL_prim_rs${XYZ}$_vf_flat(j,     k, l, contxe + i)
                    vel_R_acc(i) = qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, contxe + i)
                  end do

                  vel_L_rms_acc = 0d0; vel_R_rms_acc = 0d0

                  !$acc loop seq
                  do i = 1, num_dims
                    vel_L_rms_acc = vel_L_rms_acc + vel_L_acc(i)**2d0
                    vel_R_rms_acc = vel_R_rms_acc + vel_R_acc(i)**2d0
                  end do

                  vel_L_rms_acc = sqrt(vel_L_rms_acc)
                  vel_R_rms_acc = sqrt(vel_R_rms_acc)

                  !$acc loop seq
                  do i = 1, num_fluids
                    alpha_L_acc(i) = qL_prim_rs${XYZ}$_vf_flat(j,     k, l, E_idx + i)
                    alpha_R_acc(i) = qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, E_idx + i)
                  end do

                  pres_L_acc = qL_prim_rs${XYZ}$_vf_flat(j,     k, l, E_idx)
                  pres_R_acc = qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, E_idx)

                  rho_L_acc    = 0d0
                  gamma_L_acc  = 0d0
                  pi_inf_L_acc = 0d0

                  !$acc loop seq
                  do i = 1, num_fluids
                    rho_L_acc    = rho_L_acc    + alpha_rho_L_acc(i)
                    gamma_L_acc  = gamma_L_acc  + alpha_L_acc(i)*gammas(i)
                    pi_inf_L_acc = pi_inf_L_acc + alpha_L_acc(i)*pi_infs(i)
                  end do

                  rho_R_acc    = 0d0
                  gamma_R_acc  = 0d0
                  pi_inf_R_acc = 0d0

                  !$acc loop seq
                  do i = 1, num_fluids
                    rho_R_acc    = rho_R_acc    + alpha_rho_R_acc(i)
                    gamma_R_acc  = gamma_R_acc  + alpha_R_acc(i)*gammas(i)
                    pi_inf_R_acc = pi_inf_R_acc + alpha_R_acc(i)*pi_infs(i)
                  end do

                  E_L_acc = gamma_L_acc*pres_L_acc + pi_inf_L_acc + 5d-1*rho_L_acc*vel_L_rms_acc**2d0
                  E_R_acc = gamma_R_acc*pres_R_acc + pi_inf_R_acc + 5d-1*rho_R_acc*vel_R_rms_acc**2d0

                  H_L_acc = (E_L_acc + pres_L_acc)/rho_L_acc
                  H_R_acc = (E_R_acc + pres_R_acc)/rho_R_acc

                  if(avg_state == 2) then
                    rho_avg_acc = 5d-1*(rho_L_acc + rho_R_acc)

                    !$acc loop seq
                    do i = 1, num_dims
                        vel_avg_acc(i) = 5d-1*(vel_L_acc(i) + vel_R_acc(i))
                    end do

                    H_avg_acc = 5d-1*(H_L_acc + H_R_acc)

                    gamma_avg_acc = 5d-1*(gamma_L_acc + gamma_R_acc)
                  elseif(avg_state == 1) then
                    rho_avg_acc = sqrt(rho_L_acc*rho_R_acc)

                    !$acc loop seq
                    do i = 1, num_dims
                        vel_avg_acc(i) = (sqrt(rho_L_acc)*vel_L_acc(i) + sqrt(rho_R_acc)*vel_R_acc(i))/ &
                            (sqrt(rho_L_acc) + sqrt(rho_R_acc))
                    end do

                    H_avg_acc = (sqrt(rho_L_acc)*H_L_acc + sqrt(rho_R_acc)*H_R_acc)/ &
                                (sqrt(rho_L_acc) + sqrt(rho_R_acc))

                    gamma_avg_acc = (sqrt(rho_L_acc)*gamma_L_acc + sqrt(rho_R_acc)*gamma_R_acc)/ &
                                (sqrt(rho_L_acc) + sqrt(rho_R_acc))
                  end if

                  vel_avg_rms_acc = 0d0

                  !$acc loop seq
                  do i = 1, num_dims
                    vel_avg_rms_acc = vel_avg_rms_acc + vel_avg_acc(i)**2d0
                  end do

                  vel_avg_rms_acc = sqrt(vel_avg_rms_acc)

                  if (mixture_err) then
                    if ((H_avg_acc - 5d-1*vel_avg_rms_acc**2d0) < 0d0) then
                        c_avg_acc = sgm_eps
                    else
                        c_avg_acc = sqrt((H_avg_acc - 5d-1*vel_avg_rms_acc**2d0)/gamma_avg_acc)
                    end if
                  else
                      c_avg_acc = sqrt((H_avg_acc - 5d-1*vel_avg_rms_acc**2d0)/gamma_avg_acc)
                  end if

                  if (alt_soundspeed) then
                    blkmod1 = ((gammas(1) + 1d0)*pres_L_acc + &
                               pi_infs(1))/gammas(1)

                    blkmod2 = ((gammas(2) + 1d0)*pres_L_acc + &
                               pi_infs(2))/gammas(2)

                    c_L_acc = 1d0/(rho_L_acc*(alpha_L_acc(1)/blkmod1 + alpha_L_acc(2)/blkmod2))

                    blkmod1 = ((gammas(1) + 1d0)*pres_R_acc + &
                               pi_infs(1))/gammas(1)

                    blkmod2 = ((gammas(2) + 1d0)*pres_R_acc + &
                               pi_infs(2))/gammas(2)

                    c_R_acc = 1d0/(rho_R_acc*(alpha_R_acc(1)/blkmod1 + alpha_R_acc(2)/blkmod2))
                  elseif (model_eqns == 3) then
                    c_L_acc = 0d0
                    c_R_acc = 0d0

                    !$acc loop seq
                    do i = 1, num_fluids
                        c_L_acc = c_L_acc + qL_prim_rs${XYZ}$_vf_flat(j, k, l, i + advxb - 1)*(1d0/gammas(i) + 1d0)* &
                                (qL_prim_rs${XYZ}$_vf_flat(j, k, l, E_idx) + pi_infs(i)/(gammas(i) + 1d0))

                        c_R_acc = c_R_acc + qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, i + advxb - 1)*(1d0/gammas(i) + 1d0)* &
                              (qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, E_idx) + pi_infs(i)/(gammas(i) + 1d0))
                    end do

                    c_L_acc = c_L_acc/rho_L_acc
                    c_R_acc = c_R_acc/rho_R_acc
                  elseif ((model_eqns == 4) .or. (model_eqns == 2 .and. bubbles)) then
                    ! Sound speed for bubble mmixture to order O(\alpha)

                    if (mpp_lim .and. (num_fluids > 1)) then
                      c_L_acc = (1d0/gamma_L_acc + 1d0)* &
                              (pres_L_acc + pi_inf_L_acc)/rho_L_acc
                      c_R_acc = (1d0/gamma_R_acc + 1d0)* &
                            (pres_R_acc + pi_inf_R_acc)/rho_R_acc
                    else
                      c_L_acc = &
                        (1d0/gamma_L_acc + 1d0)* &
                        (pres_L_acc + pi_inf_L_acc)/ &
                        (rho_L_acc*(1d0 - alpha_L_acc(num_fluids)))
                      c_R_acc = &
                        (1d0/gamma_R_acc + 1d0)* &
                        (pres_R_acc + pi_inf_R_acc)/ &
                        (rho_R_acc*(1d0 - alpha_R_acc(num_fluids)))
                    end if
                  else
                    c_L_acc = ((H_L_acc - 5d-1*vel_L_rms_acc**2d0)/gamma_L_acc)

                    c_R_acc = ((H_R_acc - 5d-1*vel_R_rms_acc**2d0)/gamma_R_acc)
                  end if

                  if (mixture_err .and. c_L_acc < 0d0) then
                    c_L_acc = 100.d0*sgm_eps
                  else
                    c_L_acc = sqrt(c_L_acc)
                  end if

                  if (mixture_err .and. c_R_acc < 0d0) then
                    c_R_acc = 100.d0*sgm_eps
                  else
                    c_R_acc = sqrt(c_R_acc)
                  end if

                  if(wave_speeds == 1) then
                    s_L_acc = min(vel_L_acc(dir_idx(1)) - c_L_acc, vel_R_acc(dir_idx(1)) - c_R_acc)
                    s_R_acc = max(vel_R_acc(dir_idx(1)) + c_R_acc, vel_L_acc(dir_idx(1)) + c_L_acc)

                    s_S_acc = (pres_R_acc - pres_L_acc + rho_L_acc*vel_L_acc(dir_idx(1))* &
                       (s_L_acc - vel_L_acc(dir_idx(1))) - &
                       rho_R_acc*vel_R_acc(dir_idx(1))* &
                       (s_R_acc - vel_R_acc(dir_idx(1)))) &
                      /(rho_L_acc*(s_L_acc - vel_L_acc(dir_idx(1))) - &
                        rho_R_acc*(s_R_acc - vel_R_acc(dir_idx(1))))
                  elseif(wave_speeds == 2) then
                    pres_SL = 5d-1*(pres_L_acc + pres_R_acc+ rho_avg_acc*c_avg_acc* &
                        (vel_L_acc(dir_idx(1)) - &
                            vel_R_acc(dir_idx(1))))

                    pres_SR = pres_SL

                    Ms_L = max(1d0, sqrt(1d0 + ((5d-1 + gamma_L_acc)/(1d0 + gamma_L_acc))* &
                                         (pres_SL/pres_L_acc - 1d0)*pres_L_acc/ &
                                         ((pres_L_acc + pi_inf_L_acc/(1d0 + gamma_L_acc)))))
                    Ms_R = max(1d0, sqrt(1d0 + ((5d-1 + gamma_R_acc)/(1d0 + gamma_R_acc))* &
                                         (pres_SR/pres_R_acc - 1d0)*pres_R_acc/ &
                                         ((pres_R_acc + pi_inf_R_acc/(1d0 + gamma_R_acc)))))

                    s_L_acc = vel_L_acc(dir_idx(1)) - c_L_acc*Ms_L
                    s_R_acc = vel_R_acc(dir_idx(1)) + c_R_acc*Ms_R

                    s_S_acc = 5d-1*((vel_L_acc(dir_idx(1)) + vel_R_acc(dir_idx(1))) + &
                                (pres_L_acc - pres_R_acc)/ &
                                            (rho_avg_acc*c_avg_acc))
                  end if

                  s_M_acc = min(0d0, s_L_acc); s_P_acc = max(0d0, s_R_acc)

                  xi_M_acc = (5d-1 + sign(5d-1, s_L_acc)) &
                         + (5d-1 - sign(5d-1, s_L_acc)) &
                         * (5d-1 + sign(5d-1, s_R_acc))
                  xi_P_acc = (5d-1 - sign(5d-1, s_R_acc)) &
                         + (5d-1 - sign(5d-1, s_L_acc)) &
                         * (5d-1 + sign(5d-1, s_R_acc))


                  ! Mass
                  !$acc loop seq
                  do i = 1, contxe
                    flux_rs${XYZ}$_vf_flat(j, k, l, i) = &
                      (s_M_acc*alpha_rho_R_acc(i)*vel_R_acc(dir_idx(1)) &
                       - s_P_acc*alpha_rho_L_acc(i)*vel_L_acc(dir_idx(1)) &
                       + s_M_acc*s_P_acc*(alpha_rho_L_acc(i) &
                                  - alpha_rho_R_acc(i))) &
                      /(s_M_acc - s_P_acc)
                  end do

                  ! Momentum
                  if (bubbles) then
                    !$acc loop seq
                    do i = 1, num_dims
                      flux_rs${XYZ}$_vf_flat(j, k, l, contxe + dir_idx(i)) = &
                          (s_M_acc*(rho_R_acc*vel_R_acc(dir_idx(1)) &
                                *vel_R_acc(dir_idx(i)) &
                                + dir_flg(dir_idx(i))*(pres_R_acc - ptilde_R_acc)) &
                           - s_P_acc*(rho_L_acc*vel_L_acc(dir_idx(1)) &
                                  *vel_L_acc(dir_idx(i)) &
                                  + dir_flg(dir_idx(i))*(pres_L_acc - ptilde_L_acc)) &
                           + s_M_acc*s_P_acc*(rho_L_acc*vel_L_acc(dir_idx(i)) &
                                      - rho_R_acc*vel_R_acc(dir_idx(i)))) &
                          /(s_M_acc - s_P_acc)
                    end do
                  else
                    !$acc loop seq
                    do i = 1, num_dims
                      flux_rs${XYZ}$_vf_flat(j, k, l, contxe + dir_idx(i)) = &
                        (s_M_acc*(rho_R_acc*vel_R_acc(dir_idx(1)) &
                              *vel_R_acc(dir_idx(i)) &
                              + dir_flg(dir_idx(i))*pres_R_acc) &
                         - s_P_acc*(rho_L_acc*vel_L_acc(dir_idx(1)) &
                                *vel_L_acc(dir_idx(i)) &
                                + dir_flg(dir_idx(i))*pres_L_acc) &
                         + s_M_acc*s_P_acc*(rho_L_acc*vel_L_acc(dir_idx(i)) &
                                    - rho_R_acc*vel_R_acc(dir_idx(i)))) &
                        /(s_M_acc - s_P_acc)
                    end do
                  end if

                  ! Energy
                  if (bubbles) then
                    flux_rs${XYZ}$_vf_flat(j, k, l, E_idx) = &
                      (s_M_acc*vel_R_acc(dir_idx(1))*(E_R_acc + pres_R_acc- ptilde_R_acc) &
                       - s_P_acc*vel_L_acc(dir_idx(1))*(E_L_acc + pres_L_acc - ptilde_L_acc) &
                       + s_M_acc*s_P_acc*(E_L_acc - E_R_acc)) &
                      /(s_M_acc - s_P_acc)
                  else
                    flux_rs${XYZ}$_vf_flat(j, k, l, E_idx) = &
                      (s_M_acc*vel_R_acc(dir_idx(1))*(E_R_acc + pres_R_acc) &
                       - s_P_acc*vel_L_acc(dir_idx(1))*(E_L_acc + pres_L_acc) &
                       + s_M_acc*s_P_acc*(E_L_acc - E_R_acc)) &
                      /(s_M_acc - s_P_acc)
                  end if

                  ! Advection
                  !$acc loop seq
                  do i = advxb, advxe
                    flux_rs${XYZ}$_vf_flat(j, k, l, i) = &
                      (qL_prim_rsx_vf_flat(j, k, l, i) &
                       - qR_prim_rsx_vf_flat(j + 1, k, l, i)) &
                      *s_M_acc*s_P_acc/(s_M_acc - s_P_acc)
                    flux_src_rsx_vf_flat(j, k, l, i) = &
                      (s_M_acc*qR_prim_rsx_vf_flat(j + 1, k, l, i) &
                       - s_P_acc*qL_prim_rsx_vf_flat(j, k, l, i)) &
                      /(s_M_acc - s_P_acc)
                  end do

                  ! Div(U)?
                  !$acc loop seq
                  do i = 1, num_dims
                      vel_src_rsx_vf_flat(j, k, l, dir_idx(i)) = &
                          (xi_M_acc*(rho_L_acc*vel_L_acc(dir_idx(i))* &
                                 (s_L_acc - vel_L_acc(dir_idx(1))) - &
                                 pres_L_acc*dir_flg(dir_idx(i))) - &
                           xi_P_acc*(rho_R_acc*vel_R_acc(dir_idx(i))* &
                                 (s_R_acc - vel_R_acc(dir_idx(1))) - &
                                 pres_R_acc*dir_flg(dir_idx(i)))) &
                          /(xi_M_acc*rho_L_acc*(s_L_acc - vel_L_acc(dir_idx(1))) - &
                            xi_P_acc*rho_R_acc*(s_R_acc - vel_R_acc(dir_idx(1))))
                  end do

                  if (bubbles) then
                    ! From HLLC: Kills mass transport @ bubble gas density
                    if (num_fluids > 1) then
                        flux_rs${XYZ}$_vf_flat(j, k, l, contxe) = 0d0
                    end if
                  end if
                end do
              end do
            end do
          end if

#:endfor

        call s_finalize_riemann_solver(flux_vf, flux_src_vf, &
                                       flux_gsrc_vf, &
                                       norm_dir, ix, iy, iz)

    end subroutine s_hll_riemann_solver_acc

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


        subroutine s_hllc_riemann_solver_acc(qL_prim_rsx_vf_flat, qL_prim_rsy_vf_flat, qL_prim_rsz_vf_flat, dqL_prim_dx_vf, & ! ------
                                     dqL_prim_dy_vf, &
                                     dqL_prim_dz_vf, &
                                     gm_alphaL_vf, &
                                     qR_prim_rsx_vf_flat, qR_prim_rsy_vf_flat, qR_prim_rsz_vf_flat, dqR_prim_dx_vf, &
                                     dqR_prim_dy_vf, &
                                     dqR_prim_dz_vf, &
                                     gm_alphaR_vf, &
                                     q_prim_vf, &
                                     flux_vf, flux_src_vf, &
                                     flux_gsrc_vf, &
                                     norm_dir, ix, iy, iz)

        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:), intent(INOUT) :: qL_prim_rsx_vf_flat, qL_prim_rsy_vf_flat, qL_prim_rsz_vf_flat, qR_prim_rsx_vf_flat, qR_prim_rsy_vf_flat, qR_prim_rsz_vf_flat
        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf

        type(scalar_field), &
            allocatable, dimension(:), &
            intent(INOUT) :: dqL_prim_dx_vf, dqR_prim_dx_vf, &
                             dqL_prim_dy_vf, dqR_prim_dy_vf, &
                             dqL_prim_dz_vf, dqR_prim_dz_vf, &
                             gm_alphaL_vf, gm_alphaR_vf

        ! Intercell fluxes
        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: flux_vf, flux_src_vf, flux_gsrc_vf

        integer, intent(IN) :: norm_dir
        type(bounds_info), intent(IN) :: ix, iy, iz


        real(kind(0d0)),dimension(2)   :: alpha_rho_L_acc, alpha_rho_R_acc = 0d0
        real(kind(0d0))                              ::       rho_L_acc, rho_R_acc
        real(kind(0d0)), dimension(3)   :: vel_L_acc, vel_R_acc = 0d0
        real(kind(0d0))                              ::      pres_L_acc, pres_R_acc
        real(kind(0d0))                              ::         E_L_acc, E_R_acc
        real(kind(0d0))                              ::         H_L_acc, H_R_acc
        real(kind(0d0)), dimension(2)   :: alpha_L_acc, alpha_R_acc = 0d0
        real(kind(0d0))                              ::         Y_L_acc, Y_R_acc
        real(kind(0d0))                              ::     gamma_L_acc, gamma_R_acc
        real(kind(0d0))                              ::    pi_inf_L_acc, pi_inf_R_acc
        real(kind(0d0))                              ::         c_L_acc, c_R_acc

        real(kind(0d0))                                 :: rho_avg_acc
        real(kind(0d0)),dimension(3)   :: vel_avg_acc = 0d0
        real(kind(0d0))                                 :: H_avg_acc
        real(kind(0d0))                                 :: gamma_avg_acc
        real(kind(0d0))                                 :: c_avg_acc

        real(kind(0d0))  :: s_L_acc, s_R_acc, s_M_acc, s_P_acc, s_S_acc
        real(kind(0d0)) :: xi_L_acc, xi_R_acc !< Left and right wave speeds functions
        real(kind(0d0)) :: xi_M_acc, xi_P_acc


        real(kind(0d0))                              ::       nbub_L_acc, nbub_R_acc
        real(kind(0d0)), dimension(nb)  ::         R0_L_acc, R0_R_acc
        real(kind(0d0)), dimension(nb)   ::         V0_L_acc, V0_R_acc
        real(kind(0d0)), dimension(nb)   ::         P0_L_acc, P0_R_acc
        real(kind(0d0)), dimension(nb)  ::        pbw_L_acc, pbw_R_acc
        real(kind(0d0)), dimension(3,6) ::       moms_L_acc, moms_R_acc
        real(kind(0d0))                              ::     ptilde_L_acc, ptilde_R_acc

        real(kind(0d0)) :: PbwR3Lbar_acc, Pbwr3Rbar_acc
        real(kind(0d0)) :: R3Lbar_acc, R3Rbar_acc
        real(kind(0d0)) :: R3V2Lbar_acc, R3V2Rbar_acc

        real(kind(0d0))  :: vel_L_rms_acc, vel_R_rms_acc, vel_avg_rms_acc
        real(kind(0d0)) :: blkmod1_acc, blkmod2_acc
        real(kind(0d0)) :: rho_Star_acc, E_Star_acc, p_Star_acc, p_K_Star_acc
        real(kind(0d0)) :: pres_SL, pres_SR, Ms_L, Ms_R
        integer :: i, j, k, l !< Generic loop iterators
        integer :: idx1, idxi



        ! Populating the buffers of the left and right Riemann problem
        ! states variables, based on the choice of boundary conditions
        call s_populate_riemann_states_variables_buffers( &
            qL_prim_rsx_vf_flat, qL_prim_rsy_vf_flat, qL_prim_rsz_vf_flat, dqL_prim_dx_vf, &
            dqL_prim_dy_vf, &
            dqL_prim_dz_vf, &
            gm_alphaL_vf, &
            qR_prim_rsx_vf_flat, qR_prim_rsy_vf_flat, qR_prim_rsz_vf_flat, dqR_prim_dx_vf, &
            dqR_prim_dy_vf, &
            dqR_prim_dz_vf, &
            gm_alphaR_vf, &
            norm_dir, ix, iy, iz)

        ! Reshaping inputted data based on dimensional splitting direction
        call s_initialize_riemann_solver(&
                                         q_prim_vf, &
                                         flux_vf, flux_src_vf, &
                                         flux_gsrc_vf, &
                                             norm_dir, ix, iy, iz)
#:for NORM_DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]

            if (norm_dir == ${NORM_DIR}$) then
                if(model_eqns == 3) then
                    !ME3
                
                do l = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = is1%beg, is1%end

                            !$acc loop seq
                            do i = 1, contxe
                                alpha_rho_L_acc(i) = qL_prim_rs${XYZ}$_vf_flat(j,     k, l, i)
                                alpha_rho_R_acc(i) = qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, i)
                            end do

                            !$acc loop seq
                            do i = 1, num_dims
                                vel_L_acc(i) = qL_prim_rs${XYZ}$_vf_flat(j,     k, l, contxe + i)
                                vel_R_acc(i) = qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, contxe + i)
                            end do

                            vel_L_rms_acc = 0d0; vel_R_rms_acc = 0d0
                !$acc loop seq
                            do i = 1, num_dims
                                vel_L_rms_acc = vel_L_rms_acc + vel_L_acc(i)**2d0
                                vel_R_rms_acc = vel_R_rms_acc + vel_R_acc(i)**2d0
                            end do
                            vel_L_rms_acc = sqrt(vel_L_rms_acc)
                            vel_R_rms_acc = sqrt(vel_R_rms_acc)


                !$acc loop seq
                            do i = 1, num_fluids
                                alpha_L_acc(i) = qL_prim_rs${XYZ}$_vf_flat(j, k, l, E_idx + i)
                                alpha_R_acc(i) = qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, E_idx + i)
                            end do

                            pres_L_acc = qL_prim_rs${XYZ}$_vf_flat(j, k, l, E_idx)
                            pres_R_acc = qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, E_idx)

                            rho_L_acc = 0d0
                            gamma_L_acc = 0d0
                            pi_inf_L_acc = 0d0
                !$acc loop seq
                            do i = 1, num_fluids
                                rho_L_acc = rho_L_acc + alpha_rho_L_acc(i)
                                gamma_L_acc = gamma_L_acc+ alpha_L_acc(i)*gammas(i)
                                pi_inf_L_acc = pi_inf_L_acc + alpha_L_acc(i)*pi_infs(i)
                            end do

                            rho_R_acc = 0d0
                            gamma_R_acc = 0d0
                            pi_inf_R_acc = 0d0
                !$acc loop seq
                            do i = 1, num_fluids
                                rho_R_acc = rho_R_acc + alpha_rho_R_acc(i)
                                gamma_R_acc = gamma_R_acc + alpha_R_acc(i)*gammas(i)
                                pi_inf_R_acc = pi_inf_R_acc + alpha_R_acc(i)*pi_infs(i)
                            end do


                            E_L_acc = gamma_L_acc*pres_L_acc + pi_inf_L_acc + 5d-1*rho_L_acc*vel_L_rms_acc**2d0

                            E_R_acc = gamma_R_acc*pres_R_acc + pi_inf_R_acc + 5d-1*rho_R_acc*vel_R_rms_acc**2d0

                            H_L_acc = (E_L_acc + pres_L_acc)/rho_L_acc
                            H_R_acc = (E_R_acc + pres_R_acc)/rho_R_acc
                            if(avg_state == 2) then

                                rho_avg_acc = 5d-1*(rho_L_acc + rho_R_acc)
                !$acc loop seq
                                do i = 1, num_dims
                                    vel_avg_acc(i) = 5d-1*(vel_L_acc(i) + vel_R_acc(i))
                                end do

                                H_avg_acc = 5d-1*(H_L_acc + H_R_acc)

                                gamma_avg_acc = 5d-1*(gamma_L_acc + gamma_R_acc)

                            elseif(avg_state == 1) then

                                rho_avg_acc = sqrt(rho_L_acc*rho_R_acc)
                !$acc loop seq
                                do i = 1, num_dims
                                    vel_avg_acc(i) = (sqrt(rho_L_acc)*vel_L_acc(i) + sqrt(rho_R_acc)*vel_R_acc(i))/ &
                                        (sqrt(rho_L_acc) + sqrt(rho_R_acc))
                                end do

                                H_avg_acc = (sqrt(rho_L_acc)*H_L_acc + sqrt(rho_R_acc)*H_R_acc)/ &
                                    (sqrt(rho_L_acc) + sqrt(rho_R_acc))

                                gamma_avg_acc = (sqrt(rho_L_acc)*gamma_L_acc + sqrt(rho_R_acc)*gamma_R_acc)/ &
                                    (sqrt(rho_L_acc) + sqrt(rho_R_acc))
                            end if

                            vel_avg_rms_acc = 0d0
                !$acc loop seq
                            do i = 1, num_dims
                                vel_avg_rms_acc = vel_avg_rms_acc + vel_avg_acc(i)**2d0
                            end do
                            vel_avg_rms_acc = sqrt(vel_avg_rms_acc)

                            if (mixture_err) then
                                if ((H_avg_acc - 5d-1*vel_avg_rms_acc**2d0) < 0d0) then
                                    c_avg_acc = sgm_eps
                                else

                                    c_avg_acc = sqrt((H_avg_acc - 5d-1*vel_avg_rms_acc**2d0)/gamma_avg_acc)
                                end if
                            else

                                c_avg_acc = sqrt((H_avg_acc - 5d-1*vel_avg_rms_acc**2d0)/gamma_avg_acc)
                            end if

                            if (alt_soundspeed) then


                                blkmod1_acc = ((gammas(1) + 1d0)*pres_L_acc + &
                                        pi_infs(1))/gammas(1)
                                blkmod2_acc = ((gammas(2) + 1d0)*pres_L_acc + &
                                        pi_infs(2))/gammas(2)
                                c_L_acc = 1d0/(rho_L_acc*(alpha_L_acc(1)/blkmod1_acc + alpha_L_acc(2)/blkmod2_acc))

                                blkmod1_acc = ((gammas(1) + 1d0)*pres_R_acc + &
                                        pi_infs(1))/gammas(1)
                                blkmod2_acc = ((gammas(2) + 1d0)*pres_R_acc + &
                                        pi_infs(2))/gammas(2)
                                c_R_acc = 1d0/(rho_R_acc*(alpha_R_acc(1)/blkmod1_acc + alpha_R_acc(2)/blkmod2_acc))

                            else
                                c_L_acc = 0d0
                                c_R_acc = 0d0
                !$acc loop seq
                                do i = 1, num_fluids
                                    c_L_acc = c_L_acc + qL_prim_rs${XYZ}$_vf_flat(j, k, l, i + advxb - 1)*(1d0/gammas(i) + 1d0)* &
                                        (qL_prim_rs${XYZ}$_vf_flat(j, k, l, E_idx) + pi_infs(i)/(gammas(i) + 1d0))
                                    c_R_acc = c_R_acc + qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, i + advxb - 1)*(1d0/gammas(i) + 1d0)* &
                                        (qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, E_idx) + pi_infs(i)/(gammas(i) + 1d0))
                                end do
                                c_L_acc = c_L_acc/rho_L_acc
                                c_R_acc = c_R_acc/rho_R_acc
                            end if


                            if (mixture_err .and. c_L_acc < 0d0) then
                                c_L_acc = 100.d0*sgm_eps
                            else
                                c_L_acc = sqrt(c_L_acc)
                            end if
                            if (mixture_err .and. c_R_acc < 0d0) then
                                c_R_acc = 100.d0*sgm_eps
                            else
                                c_R_acc = sqrt(c_R_acc)
                            end if

                            if(wave_speeds == 1) then
                                s_L_acc = min(vel_L_acc(dir_idx(1)) - c_L_acc, vel_R_acc(dir_idx(1)) - c_R_acc)
                                s_R_acc = max(vel_R_acc(dir_idx(1)) + c_R_acc, vel_L_acc(dir_idx(1)) + c_L_acc)

                                s_S_acc = (pres_R_acc - pres_L_acc + rho_L_acc*vel_L_acc(dir_idx(1))* &
                                (s_L_acc - vel_L_acc(dir_idx(1))) - &
                                rho_R_acc*vel_R_acc(dir_idx(1))* &
                                (s_R_acc - vel_R_acc(dir_idx(1)))) &
                                /(rho_L_acc*(s_L_acc - vel_L_acc(dir_idx(1))) - &
                                    rho_R_acc*(s_R_acc - vel_R_acc(dir_idx(1))))
                            elseif(wave_speeds == 2) then
                                pres_SL = 5d-1*(pres_L_acc + pres_R_acc+ rho_avg_acc*c_avg_acc* &
                                    (vel_L_acc(dir_idx(1)) - &
                                        vel_R_acc(dir_idx(1))))

                                pres_SR = pres_SL

                                Ms_L = max(1d0, sqrt(1d0 + ((5d-1 + gamma_L_acc)/(1d0 + gamma_L_acc))* &
                                                    (pres_SL/pres_L_acc - 1d0)*pres_L_acc/ &
                                                    ((pres_L_acc + pi_inf_L_acc/(1d0 + gamma_L_acc)))))
                                Ms_R = max(1d0, sqrt(1d0 + ((5d-1 + gamma_R_acc)/(1d0 + gamma_R_acc))* &
                                                    (pres_SR/pres_R_acc - 1d0)*pres_R_acc/ &
                                                    ((pres_R_acc + pi_inf_R_acc/(1d0 + gamma_R_acc)))))

                                s_L_acc = vel_L_acc(dir_idx(1)) - c_L_acc*Ms_L
                                s_R_acc = vel_R_acc(dir_idx(1)) + c_R_acc*Ms_R

                                s_S_acc = 5d-1*((vel_L_acc(dir_idx(1)) + vel_R_acc(dir_idx(1))) + &
                                            (pres_L_acc - pres_R_acc)/ &
                                                        (rho_avg_acc*c_avg_acc))
                            end if

                            if (s_L_acc >= 0d0) then
                                p_Star_acc = pres_L_acc ! Only usefull to recalculate the radial momentum geometric source flux
                !$acc loop seq
                                do i = 1, num_fluids
                                    flux_rs${XYZ}$_vf_flat(j, k, l, i + advxb - 1) = &
                                        qL_prim_rs${XYZ}$_vf_flat(j, k, l, i + advxb - 1)*s_S_acc

                                    flux_rs${XYZ}$_vf_flat(j, k, l, i + contxb - 1) = &
                                        qL_prim_rs${XYZ}$_vf_flat(j, k, l, i + contxb - 1)*vel_L_acc(dir_idx(1))

                                    flux_rs${XYZ}$_vf_flat(j, k, l, i + intxb - 1) = &
                                        qL_prim_rs${XYZ}$_vf_flat(j, k, l, i + advxb - 1)* &
                                        (gammas(i)*pres_L_acc + pi_infs(i))*vel_L_acc(dir_idx(1))
                                end do
                !$acc loop seq
                                do i = 1, num_dims
                                    flux_rs${XYZ}$_vf_flat(j, k, l, momxb - 1 + dir_idx(i)) = &
                                        rho_L_acc*vel_L_acc(dir_idx(1))*vel_L_acc(dir_idx(i)) + dir_flg(dir_idx(i))*pres_L_acc

                                    vel_src_rs${XYZ}$_vf_flat(j, k, l, dir_idx(i)) = vel_L_acc(dir_idx(i)) + &
                                                                            dir_flg(dir_idx(i))*(s_S_acc - vel_L_acc(dir_idx(i)))
                                    ! Compute the star velocities for the non-conservative terms
                                end do
                                flux_rs${XYZ}$_vf_flat(j, k, l, E_idx) = (E_L_acc + pres_L_acc)*vel_L_acc(dir_idx(1))

                                ! Compute right solution state
                            else if (s_R_acc <= 0d0) then
                                p_Star_acc = pres_R_acc
                                ! Only usefull to recalculate the radial momentum geometric source flux
                !$acc loop seq
                                do i = 1, num_fluids
                                    flux_rs${XYZ}$_vf_flat(j, k, l, i + advxb - 1) = &
                                        qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, i + advxb - 1)*s_S_acc

                                    flux_rs${XYZ}$_vf_flat(j, k, l, i + contxb - 1) = &
                                        qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, i + contxb - 1)*vel_R_acc(dir_idx(1))

                                    flux_rs${XYZ}$_vf_flat(j, k, l, i + intxb - 1) = &
                                        qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, i + advxb - 1)* &
                                        (gammas(i)*pres_R_acc + pi_infs(i))*vel_R_acc(dir_idx(1))
                                end do
                !$acc loop seq
                                do i = 1, num_dims
                                    flux_rs${XYZ}$_vf_flat(j, k, l, momxb - 1 + dir_idx(i)) = &
                                        rho_R_acc*vel_R_acc(dir_idx(1))*vel_R_acc(dir_idx(i)) + dir_flg(dir_idx(i))*pres_R_acc

                                    vel_src_rs${XYZ}$_vf_flat(j, k, l, dir_idx(i)) = vel_R_acc(dir_idx(i)) + &
                                                                            dir_flg(dir_idx(i))*(s_S_acc - vel_R_acc(dir_idx(i)))
                                    ! Compute the star velocities for the non-conservative terms
                                end do
                                flux_rs${XYZ}$_vf_flat(j, k, l, E_idx) = (E_R_acc + pres_R_acc)*vel_R_acc(dir_idx(1))

                                ! Compute left star solution state
                            else if (s_S_acc >= 0d0) then
                                xi_L_acc = (s_L_acc - vel_L_acc(dir_idx(1)))/(s_L_acc - s_S_acc)
                                rho_Star_acc = rho_L_acc*xi_L_acc
                                E_Star_acc = xi_L_acc*(E_L_acc + (s_S_acc - vel_L_acc(dir_idx(1)))* &
                                            (rho_L_acc*s_S_acc + pres_L_acc/(s_L_acc - vel_L_acc(dir_idx(1)))))
                                p_Star_acc = rho_L_acc*(s_L_acc - vel_L_acc(dir_idx(1)))*(s_S_acc - vel_L_acc(dir_idx(1))) + pres_L_acc
                !$acc loop seq
                                do i = 1, num_fluids
                                    p_K_Star_acc = (pres_L_acc + pi_infs(i)/(1d0 + gammas(i)))* &
                                            xi_L_acc**(1d0/gammas(i) + 1d0) - pi_infs(i)/(1d0 + gammas(i))

                                    flux_rs${XYZ}$_vf_flat(j, k, l, i + advxb - 1) = &
                                        qL_prim_rs${XYZ}$_vf_flat(j, k, l, i + advxb - 1)*s_S_acc

                                    flux_rs${XYZ}$_vf_flat(j, k, l, i + contxb - 1) = &
                                        qL_prim_rs${XYZ}$_vf_flat(j, k, l, i + contxb - 1)*xi_L_acc*s_S_acc

                                    flux_rs${XYZ}$_vf_flat(j, k, l, i + intxb - 1) = &
                                        qL_prim_rs${XYZ}$_vf_flat(j, k, l, i + advxb - 1)* &
                                        (gammas(i)*p_K_Star_acc + pi_infs(i))*s_S_acc
                                end do
                !$acc loop seq
                                do i = 1, num_dims
                                    flux_rs${XYZ}$_vf_flat(j, k, l, momxb - 1 + dir_idx(i)) = &
                                        rho_Star_acc*s_S_acc*(s_S_acc*dir_flg(dir_idx(i)) + vel_L_acc(dir_idx(i))* &
                                                    (1d0 - dir_flg(dir_idx(i)))) + dir_flg(dir_idx(i))*p_Star_acc

                                    vel_src_rs${XYZ}$_vf_flat(j, k, l, dir_idx(i)) = vel_L_acc(dir_idx(i)) + &
                                                                            dir_flg(dir_idx(i))*(s_S_acc*xi_L_acc - vel_L_acc(dir_idx(i)))
                                    ! Compute the star velocities for the non-conservative terms
                                end do
                                flux_rs${XYZ}$_vf_flat(j, k, l, E_idx) = (E_Star_acc + p_Star_acc)*s_S_acc

                                ! Compute right star solution state
                            else
                                xi_R_acc = (s_R_acc - vel_R_acc(dir_idx(1)))/(s_R_acc - s_S_acc)

                                rho_Star_acc = rho_R_acc*xi_R_acc

                                E_Star_acc = xi_R_acc*(E_R_acc + (s_S_acc - vel_R_acc(dir_idx(1)))* &
                                            (rho_R_acc*s_S_acc + pres_R_acc/(s_R_acc - vel_R_acc(dir_idx(1)))))

                                p_Star_acc = rho_R_acc*(s_R_acc - vel_R_acc(dir_idx(1)))*(s_S_acc - vel_R_acc(dir_idx(1))) + pres_R_acc
                !$acc loop seq
                                do i = 1, num_fluids
                                    p_K_Star_acc = (pres_R_acc +  pi_infs(i)/(1d0 + gammas(i)))* &
                                            xi_R_acc**(1d0/gammas(i) + 1d0) - pi_infs(i)/(1d0 + gammas(i))

                                    flux_rs${XYZ}$_vf_flat(j, k, l, i + advxb - 1) = &
                                        qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, i + advxb - 1)*s_S_acc

                                    flux_rs${XYZ}$_vf_flat(j, k, l, i + contxb - 1) = &
                                        qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, i + contxb - 1)*xi_R_acc*s_S_acc

                                    flux_rs${XYZ}$_vf_flat(j, k, l, i + intxb - 1) = &
                                        qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, i + advxb - 1)* &
                                        (gammas(i)*p_K_Star_acc + pi_infs(i))*s_S_acc
                                end do
                !$acc loop seq
                                do i = 1, num_dims
                                    flux_rs${XYZ}$_vf_flat(j, k, l, momxb - 1 + dir_idx(i)) = rho_Star_acc*s_S_acc* &
                                    (s_S_acc*dir_flg(dir_idx(i)) + vel_R_acc(dir_idx(i))*(1d0 - dir_flg(dir_idx(i)))) + &
                                    dir_flg(dir_idx(i))*p_Star_acc

                                    vel_src_rs${XYZ}$_vf_flat(j, k, l, dir_idx(i)) = vel_R_acc(dir_idx(i)) + &
                                                                            dir_flg(dir_idx(i))*(s_S_acc*xi_R_acc - vel_R_acc(dir_idx(i)))
                                    ! Compute the star velocities for the non-conservative terms
                                end do

                                flux_rs${XYZ}$_vf_flat(j, k, l, E_idx) = (E_Star_acc + p_Star_acc)*s_S_acc

                            end if

                            ! Geometrical source flux for cylindrical coordinates
                            if (cyl_coord .and. norm_dir == 2) then
                                ! Substituting the advective flux into the inviscid geometrical source flux
                !$acc loop seq
                                do i = 1, E_idx
                                    flux_gsrc_rs${XYZ}$_vf_flat(j, k, l, i) = flux_rs${XYZ}$_vf_flat(j, k, l, i)
                                end do
                !$acc loop seq
                                do i = intxb, intxe
                                    flux_gsrc_rs${XYZ}$_vf_flat(j, k, l, i) = flux_rs${XYZ}$_vf_flat(j, k, l, i)
                                end do
                                ! Recalculating the radial momentum geometric source flux (substracting the pressure part)
                                flux_gsrc_rs${XYZ}$_vf_flat(j, k, l, momxb - 1 + dir_idx(1)) = &
                                    flux_gsrc_rs${XYZ}$_vf_flat(j, k, l, momxb - 1 + dir_idx(1)) - p_Star_acc
                                ! Geometrical source of the void fraction(s) is zero
                !$acc loop seq
                                do i = advxb, advxe
                                    flux_gsrc_rs${XYZ}$_vf_flat(j, k, l, i) = 0d0
                                end do
                            end if


                        end do
                    end do
                end do
                elseif(model_eqns == 4) then
                    !ME4
                !$acc parallel loop collapse(3) gang vector default(present) private(alpha_rho_L_acc, alpha_rho_R_acc, vel_L_acc, vel_R_acc, alpha_L_acc, alpha_R_acc, vel_avg_acc)
                do l = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = is1%beg, is1%end
                !$acc loop seq
                            do i = 1, contxe
                                alpha_rho_L_acc(i) = qL_prim_rs${XYZ}$_vf_flat(j, k, l, i)
                                alpha_rho_R_acc(i) = qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, i)
                            end do

                !$acc loop seq
                            do i = 1, num_dims
                                vel_L_acc(i) = qL_prim_rs${XYZ}$_vf_flat(j, k, l, contxe + i)
                                vel_R_acc(i) = qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, contxe + i)
                            end do

                            vel_L_rms_acc = 0d0; vel_R_rms_acc = 0d0
                !$acc loop seq
                            do i = 1, num_dims
                                vel_L_rms_acc = vel_L_rms_acc + vel_L_acc(i)**2d0
                                vel_R_rms_acc = vel_R_rms_acc + vel_R_acc(i)**2d0
                            end do
                            vel_L_rms_acc = sqrt(vel_L_rms_acc)
                            vel_R_rms_acc = sqrt(vel_R_rms_acc)


                !$acc loop seq
                            do i = 1, num_fluids
                                alpha_L_acc(i) = qL_prim_rs${XYZ}$_vf_flat(j, k, l, E_idx + i)
                                alpha_R_acc(i) = qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, E_idx + i)
                            end do

                            pres_L_acc = qL_prim_rs${XYZ}$_vf_flat(j, k, l, E_idx)
                            pres_R_acc = qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, E_idx)

                            rho_L_acc = 0d0
                            gamma_L_acc = 0d0
                            pi_inf_L_acc = 0d0
                !$acc loop seq
                            do i = 1, num_fluids
                                rho_L_acc = rho_L_acc + alpha_rho_L_acc(i)
                                gamma_L_acc = gamma_L_acc+ alpha_L_acc(i)*gammas(i)
                                pi_inf_L_acc = pi_inf_L_acc + alpha_L_acc(i)*pi_infs(i)
                            end do

                            rho_R_acc = 0d0
                            gamma_R_acc = 0d0
                            pi_inf_R_acc = 0d0
                !$acc loop seq
                            do i = 1, num_fluids
                                rho_R_acc = rho_R_acc + alpha_rho_R_acc(i)
                                gamma_R_acc = gamma_R_acc + alpha_R_acc(i)*gammas(i)
                                pi_inf_R_acc = pi_inf_R_acc + alpha_R_acc(i)*pi_infs(i)
                            end do


                            E_L_acc = gamma_L_acc*pres_L_acc + pi_inf_L_acc + 5d-1*rho_L_acc*vel_L_rms_acc**2d0

                            E_R_acc = gamma_R_acc*pres_R_acc + pi_inf_R_acc + 5d-1*rho_R_acc*vel_R_rms_acc**2d0

                            H_L_acc = (E_L_acc + pres_L_acc)/rho_L_acc
                            H_R_acc = (E_R_acc + pres_R_acc)/rho_R_acc
                            if(avg_state == 2) then

                                rho_avg_acc = 5d-1*(rho_L_acc + rho_R_acc)
                !$acc loop seq
                                do i = 1, num_dims
                                    vel_avg_acc(i) = 5d-1*(vel_L_acc(i) + vel_R_acc(i))
                                end do

                                H_avg_acc = 5d-1*(H_L_acc + H_R_acc)

                                gamma_avg_acc = 5d-1*(gamma_L_acc + gamma_R_acc)

                            elseif(avg_state == 1) then

                                rho_avg_acc = sqrt(rho_L_acc*rho_R_acc)
                !$acc loop seq
                                do i = 1, num_dims
                                    vel_avg_acc(i) = (sqrt(rho_L_acc)*vel_L_acc(i) + sqrt(rho_R_acc)*vel_R_acc(i))/ &
                                        (sqrt(rho_L_acc) + sqrt(rho_R_acc))
                                end do

                                H_avg_acc = (sqrt(rho_L_acc)*H_L_acc + sqrt(rho_R_acc)*H_R_acc)/ &
                                    (sqrt(rho_L_acc) + sqrt(rho_R_acc))

                                gamma_avg_acc = (sqrt(rho_L_acc)*gamma_L_acc + sqrt(rho_R_acc)*gamma_R_acc)/ &
                                    (sqrt(rho_L_acc) + sqrt(rho_R_acc))
                            end if

                            vel_avg_rms_acc = 0d0
                !$acc loop seq
                            do i = 1, num_dims
                                vel_avg_rms_acc = vel_avg_rms_acc + vel_avg_acc(i)**2d0
                            end do
                            vel_avg_rms_acc = sqrt(vel_avg_rms_acc)

                            if (mixture_err) then
                                if ((H_avg_acc - 5d-1*vel_avg_rms_acc**2d0) < 0d0) then
                                    c_avg_acc = sgm_eps
                                else

                                    c_avg_acc = sqrt((H_avg_acc - 5d-1*vel_avg_rms_acc**2d0)/gamma_avg_acc)
                                end if
                            else

                                c_avg_acc = sqrt((H_avg_acc - 5d-1*vel_avg_rms_acc**2d0)/gamma_avg_acc)
                            end if

                            if (alt_soundspeed) then


                                blkmod1_acc = ((gammas(1) + 1d0)*pres_L_acc + &
                                        pi_infs(1))/gammas(1)
                                blkmod2_acc = ((gammas(2) + 1d0)*pres_L_acc + &
                                        pi_infs(2))/gammas(2)
                                c_L_acc = 1d0/(rho_L_acc*(alpha_L_acc(1)/blkmod1_acc + alpha_L_acc(2)/blkmod2_acc))

                                blkmod1_acc = ((gammas(1) + 1d0)*pres_R_acc + &
                                        pi_infs(1))/gammas(1)
                                blkmod2_acc = ((gammas(2) + 1d0)*pres_R_acc + &
                                        pi_infs(2))/gammas(2)
                                c_R_acc = 1d0/(rho_R_acc*(alpha_R_acc(1)/blkmod1_acc + alpha_R_acc(2)/blkmod2_acc))


                            else
                                ! Sound speed for bubble mmixture to order O(\alpha)

                                if (mpp_lim .and. (num_fluids > 1)) then
                                    c_L_acc = (1d0/gamma_L_acc + 1d0)* &
                                        (pres_L_acc + pi_inf_L_acc)/rho_L_acc
                                    c_R_acc = (1d0/gamma_R_acc + 1d0)* &
                                        (pres_R_acc + pi_inf_R_acc)/rho_R_acc
                                else
                                    c_L_acc = &
                                        (1d0/gamma_L_acc + 1d0)* &
                                        (pres_L_acc + pi_inf_L_acc)/ &
                                        (rho_L_acc*(1d0 - alpha_L_acc(num_fluids)))
                                    c_R_acc = &
                                        (1d0/gamma_R_acc + 1d0)* &
                                        (pres_R_acc + pi_inf_R_acc)/ &
                                        (rho_R_acc*(1d0 - alpha_R_acc(num_fluids)))
                                end if
                            end if

                            if (mixture_err .and. c_L_acc < 0d0) then
                                c_L_acc = 100.d0*sgm_eps
                            else
                                c_L_acc = sqrt(c_L_acc)
                            end if
                            if (mixture_err .and. c_R_acc < 0d0) then
                                c_R_acc = 100.d0*sgm_eps
                            else
                                c_R_acc = sqrt(c_R_acc)
                            end if

                            if(wave_speeds == 1) then
                                s_L_acc = min(vel_L_acc(dir_idx(1)) - c_L_acc, vel_R_acc(dir_idx(1)) - c_R_acc)
                                s_R_acc = max(vel_R_acc(dir_idx(1)) + c_R_acc, vel_L_acc(dir_idx(1)) + c_L_acc)

                                s_S_acc = (pres_R_acc - pres_L_acc + rho_L_acc*vel_L_acc(dir_idx(1))* &
                                (s_L_acc - vel_L_acc(dir_idx(1))) - &
                                rho_R_acc*vel_R_acc(dir_idx(1))* &
                                (s_R_acc - vel_R_acc(dir_idx(1)))) &
                                /(rho_L_acc*(s_L_acc - vel_L_acc(dir_idx(1))) - &
                                    rho_R_acc*(s_R_acc - vel_R_acc(dir_idx(1))))
                            elseif(wave_speeds == 2) then
                                pres_SL = 5d-1*(pres_L_acc + pres_R_acc+ rho_avg_acc*c_avg_acc* &
                                    (vel_L_acc(dir_idx(1)) - &
                                        vel_R_acc(dir_idx(1))))

                                pres_SR = pres_SL

                                Ms_L = max(1d0, sqrt(1d0 + ((5d-1 + gamma_L_acc)/(1d0 + gamma_L_acc))* &
                                                    (pres_SL/pres_L_acc - 1d0)*pres_L_acc/ &
                                                    ((pres_L_acc + pi_inf_L_acc/(1d0 + gamma_L_acc)))))
                                Ms_R = max(1d0, sqrt(1d0 + ((5d-1 + gamma_R_acc)/(1d0 + gamma_R_acc))* &
                                                    (pres_SR/pres_R_acc - 1d0)*pres_R_acc/ &
                                                    ((pres_R_acc + pi_inf_R_acc/(1d0 + gamma_R_acc)))))

                                s_L_acc = vel_L_acc(dir_idx(1)) - c_L_acc*Ms_L
                                s_R_acc = vel_R_acc(dir_idx(1)) + c_R_acc*Ms_R

                                s_S_acc = 5d-1*((vel_L_acc(dir_idx(1)) + vel_R_acc(dir_idx(1))) + &
                                            (pres_L_acc - pres_R_acc)/ &
                                                        (rho_avg_acc*c_avg_acc))
                            end if
                        ! follows Einfeldt et al.
                            ! s_M/P = min/max(0.,s_L/R)
                            s_M_acc = min(0d0, s_L_acc); s_P_acc = max(0d0, s_R_acc)

                            ! goes with q_star_L/R = xi_L/R * (variable)
                            ! xi_L/R = ( ( s_L/R - u_L/R )/(s_L/R - s_star) )
                            xi_L_acc = (s_L_acc - vel_L_acc(dir_idx(1)))/(s_L_acc - s_S_acc)
                            xi_R_acc = (s_R_acc - vel_R_acc(dir_idx(1)))/(s_R_acc - s_S_acc)

                            ! goes with numerical velocity in x/y/z directions
                            ! xi_P/M = 0.5 +/m sgn(0.5,s_star)
                            xi_M_acc = (5d-1 + sign(5d-1, s_S_acc))
                            xi_P_acc = (5d-1 - sign(5d-1, s_S_acc))

                !$acc loop seq
                            do i = 1, contxe
                                flux_rs${XYZ}$_vf_flat(j, k, l, i) = &
                                    xi_M_acc*alpha_rho_L_acc(i) &
                                    *(vel_L_acc(dir_idx(1)) + s_M_acc*(xi_L_acc - 1d0)) &
                                    + xi_P_acc*alpha_rho_R_acc(i) &
                                    *(vel_R_acc(dir_idx(1)) + s_P_acc*(xi_R_acc - 1d0))
                            end do


                            ! Momentum flux.
                            ! f = \rho u u + p I, q = \rho u, q_star = \xi * \rho*(s_star, v, w)
                            if (bubbles .neqv. .true.) then
                !$acc loop seq
                                do i = 1, num_dims
                                    flux_rs${XYZ}$_vf_flat(j, k, l, contxe + dir_idx(i)) = &
                                        xi_M_acc*(rho_L_acc*(vel_L_acc(dir_idx(1))* &
                                                    vel_L_acc(dir_idx(i)) + &
                                                    s_M_acc*(xi_L_acc*(dir_flg(dir_idx(i))*s_S_acc + &
                                                                (1d0 - dir_flg(dir_idx(i)))* &
                                                                vel_L_acc(dir_idx(i))) - vel_L_acc(dir_idx(i)))) + &
                                            dir_flg(dir_idx(i))*(pres_L_acc)) &
                                        + xi_P_acc*(rho_R_acc*(vel_R_acc(dir_idx(1))* &
                                                    vel_R_acc(dir_idx(i)) + &
                                                    s_P_acc*(xi_R_acc*(dir_flg(dir_idx(i))*s_S_acc + &
                                                                (1d0 - dir_flg(dir_idx(i)))* &
                                                                vel_R_acc(dir_idx(i))) - vel_R_acc(dir_idx(i)))) + &
                                                dir_flg(dir_idx(i))*(pres_R_acc))
                                    ! if (j==0) print*, 'flux_rs_vf', flux_rs_vf(cont_idx%end+dir_idx(i))%sf(j,k,l)
                                end do
                            else
                                ! Include p_tilde
                !$acc loop seq
                                do i = 1, num_dims
                                    flux_rs${XYZ}$_vf_flat(j, k, l, contxe + dir_idx(i)) = &
                                        xi_M_acc*(rho_L_acc*(vel_L(dir_idx(1))* &
                                                    vel_L_acc(dir_idx(i)) + &
                                                    s_M_acc*(xi_L_acc*(dir_flg(dir_idx(i))*s_S_acc + &
                                                                (1d0 - dir_flg(dir_idx(i)))* &
                                                                vel_L_acc(dir_idx(i))) - vel_L_acc(dir_idx(i)))) + &
                                            dir_flg(dir_idx(i))*(pres_L_acc - ptilde_L_acc)) &
                                        + xi_P_acc*(rho_R_acc*(vel_R_acc(dir_idx(1))* &
                                                    vel_R_acc(dir_idx(i)) + &
                                                    s_P_acc*(xi_R_acc*(dir_flg(dir_idx(i))*s_S_acc + &
                                                                (1d0 - dir_flg(dir_idx(i)))* &
                                                                vel_R_acc(dir_idx(i))) - vel_R_acc(dir_idx(i)))) + &
                                                dir_flg(dir_idx(i))*(pres_R_acc - ptilde_R_acc))
                                    ! if (j==0) print*, 'flux_rs_vf', flux_rs_vf(cont_idx%end+dir_idx(i))%sf(j,k,l)
                                end do

                            end if


                            flux_rs${XYZ}$_vf_flat(j, k, l, E_idx) = 0.d0



                !$acc loop seq
                            do i = alf_idx, alf_idx !only advect the void fraction
                                flux_rs${XYZ}$_vf_flat(j, k, l, i) = &
                                    xi_M_acc*qL_prim_rs${XYZ}$_vf_flat(j, k, l, i) &
                                    *(vel_L_acc(dir_idx(1)) + s_M_acc*(xi_L_acc - 1d0)) &
                                    + xi_P_acc*qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, i) &
                                    *(vel_R_acc(dir_idx(1)) + s_P_acc*(xi_R_acc - 1d0))
                            end do

                            ! Source for volume fraction advection equation
                !$acc loop seq
                            do i = 1, num_dims

                                vel_src_rs${XYZ}$_vf_flat(j, k, l, dir_idx(i)) = 0d0
                                !IF ( (model_eqns == 4) .or. (num_fluids==1) ) vel_src_rs_vf(dir_idx(i))%sf(j,k,l) = 0d0
                            end do

                            flux_src_rs${XYZ}$_vf_flat(j, k, l, advxb) = vel_src_rs${XYZ}$_vf_flat(j, k, l, dir_idx(1))

                            ! Add advection flux for bubble variables
                            if (bubbles) then
                !$acc loop seq
                                do i = bubxb, bubxe
                                    flux_rs${XYZ}$_vf_flat(j, k, l, i) = &
                                        xi_M_acc*nbub_L_acc*qL_prim_rs${XYZ}$_vf_flat(j, k, l, i) &
                                        *(vel_L_acc(dir_idx(1)) + s_M_acc*(xi_L_acc - 1d0)) &
                                        + xi_P_acc*nbub_R_acc*qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, i) &
                                        *(vel_R_acc(dir_idx(1)) + s_P_acc*(xi_R_acc - 1d0))
                                end do
                            end if


                            ! Geometrical source flux for cylindrical coordinates

#:if (NORM_DIR == 2)
                            if (cyl_coord) then
                                ! Substituting the advective flux into the inviscid geometrical source flux
            !$acc loop seq
                                do i = 1, E_idx
                                    flux_gsrc_rs${XYZ}$_vf_flat(j, k, l, i) = flux_rs${XYZ}$_vf_flat(j, k, l, i)
                                end do
                                ! Recalculating the radial momentum geometric source flux
                                flux_gsrc_rs${XYZ}$_vf_flat(j, k, l, contxe + dir_idx(1)) = &
                                    xi_M_acc*(rho_L_acc*(vel_L_acc(dir_idx(1))* &
                                                vel_L_acc(dir_idx(1)) + &
                                                s_M_acc*(xi_L_acc*(dir_flg(dir_idx(1))*s_S_acc + &
                                                            (1d0 - dir_flg(dir_idx(1)))* &
                                                            vel_L_acc(dir_idx(1))) - vel_L_acc(dir_idx(1))))) &
                                    + xi_P_acc*(rho_R_acc*(vel_R_acc(dir_idx(1))* &
                                                vel_R_acc(dir_idx(1)) + &
                                                s_P_acc*(xi_R_acc*(dir_flg(dir_idx(1))*s_S_acc + &
                                                            (1d0 - dir_flg(dir_idx(1)))* &
                                                            vel_R_acc(dir_idx(1))) - vel_R_acc(dir_idx(1)))))
                                ! Geometrical source of the void fraction(s) is zero
            !$acc loop seq
                                do i = advxb, advxe
                                    flux_gsrc_rs${XYZ}$_vf_flat(j, k, l, i) = 0d0
                                end do
                            end if
#:endif
#:if (NORM_DIR == 3)
                            if (grid_geometry == 3) then
                                !$acc loop seq 
                                        do i = 1, sys_size
                                            flux_gsrc_rs${XYZ}$_vf_flat(j, k, l, i) = 0d0
                                        end do
                                        flux_gsrc_rs${XYZ}$_vf_flat(j, k, l, momxb + 1) = &
                                            -xi_M_acc*(rho_L_acc*(vel_L_acc(dir_idx(1))* &
                                                          vel_L_acc(dir_idx(1)) + &
                                                          s_M_acc*(xi_L_acc*(dir_flg(dir_idx(1))*s_S_acc + &
                                                                     (1d0 - dir_flg(dir_idx(1)))* &
                                                                     vel_L_acc(dir_idx(1))) - vel_L_acc(dir_idx(1))))) &
                                            - xi_P_acc*(rho_R_acc*(vel_R_acc(dir_idx(1))* &
                                                           vel_R_acc(dir_idx(1)) + &
                                                           s_P_acc*(xi_R_acc*(dir_flg(dir_idx(1))*s_S_acc + &
                                                                      (1d0 - dir_flg(dir_idx(1)))* &
                                                                      vel_R_acc(dir_idx(1))) - vel_R_acc(dir_idx(1)))))
                                        flux_gsrc_rs${XYZ}$_vf_flat(j, k, l, momxe) = flux_rs${XYZ}$_vf_flat(j, k, l, momxb + 1)
                            end if
#:endif
                        end do
                    end do
                end do
                elseif(model_eqns == 2 .and. bubbles) then
                !$acc parallel loop collapse(3) gang vector default(present) private(R0_L_acc, R0_R_acc, V0_L_acc, V0_R_acc, P0_L_acc, P0_R_acc, pbw_L_acc, pbw_R_acc, alpha_rho_L_acc, alpha_rho_R_acc, vel_L_acc, vel_R_acc, alpha_L_acc, alpha_R_acc, vel_avg_acc)
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end


                !$acc loop seq
                                do i = 1, contxe
                                    alpha_rho_L_acc(i) = qL_prim_rs${XYZ}$_vf_flat(j, k, l, i)
                                    alpha_rho_R_acc(i) = qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, i)
                                end do

                !$acc loop seq
                                do i = 1, num_dims
                                    vel_L_acc(i) = qL_prim_rs${XYZ}$_vf_flat(j, k, l, contxe + i)
                                    vel_R_acc(i) = qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, contxe + i)
                                end do

                                vel_L_rms_acc = 0d0; vel_R_rms_acc = 0d0
                !$acc loop seq
                                do i = 1, num_dims
                                    vel_L_rms_acc = vel_L_rms_acc + vel_L_acc(i)**2d0
                                    vel_R_rms_acc = vel_R_rms_acc + vel_R_acc(i)**2d0
                                end do
                                vel_L_rms_acc = sqrt(vel_L_rms_acc)
                                vel_R_rms_acc = sqrt(vel_R_rms_acc)


                !$acc loop seq
                                do i = 1, num_fluids
                                    alpha_L_acc(i) = qL_prim_rs${XYZ}$_vf_flat(j, k, l, E_idx + i)
                                    alpha_R_acc(i) = qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, E_idx + i)
                                end do

                                pres_L_acc = qL_prim_rs${XYZ}$_vf_flat(j, k, l, E_idx)
                                pres_R_acc = qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, E_idx)

                                rho_L_acc = 0d0
                                gamma_L_acc = 0d0
                                pi_inf_L_acc = 0d0

                                if(mpp_lim .and. (num_fluids > 2)) then
                    !$acc loop seq
                                    do i = 1, num_fluids
                                        rho_L_acc = rho_L_acc + alpha_rho_L_acc(i)
                                        gamma_L_acc = gamma_L_acc+ alpha_L_acc(i)*gammas(i)
                                        pi_inf_L_acc = pi_inf_L_acc + alpha_L_acc(i)*pi_infs(i)
                                    end do
                                else if(num_fluids > 2) then
                    !$acc loop seq
                                    do i = 1, num_fluids - 1
                                        rho_L_acc = rho_L_acc + alpha_rho_L_acc(i)
                                        gamma_L_acc = gamma_L_acc+ alpha_L_acc(i)*gammas(i)
                                        pi_inf_L_acc = pi_inf_L_acc + alpha_L_acc(i)*pi_infs(i)
                                    end do
                                else
                                    rho_L_acc = alpha_rho_L_acc(1)
                                    gamma_L_acc = gammas(1)
                                    pi_inf_L_acc = pi_infs(1)
                                end if

                                rho_R_acc = 0d0
                                gamma_R_acc = 0d0
                                pi_inf_R_acc = 0d0

                                if(mpp_lim .and. (num_fluids > 2)) then
                    !$acc loop seq
                                    do i = 1, num_fluids
                                        rho_R_acc = rho_R_acc + alpha_rho_R_acc(i)
                                        gamma_R_acc = gamma_R_acc+ alpha_R_acc(i)*gammas(i)
                                        pi_inf_R_acc = pi_inf_R_acc + alpha_R_acc(i)*pi_infs(i)
                                    end do
                                else if(num_fluids > 2) then
                    !$acc loop seq
                                    do i = 1, num_fluids - 1
                                        rho_R_acc = rho_R_acc + alpha_rho_R_acc(i)
                                        gamma_R_acc = gamma_R_acc+ alpha_R_acc(i)*gammas(i)
                                        pi_inf_R_acc = pi_inf_R_acc + alpha_R_acc(i)*pi_infs(i)
                                    end do
                                else
                                    rho_R_acc = alpha_rho_R_acc(1)
                                    gamma_R_acc = gammas(1)
                                    pi_inf_R_acc = pi_infs(1)
                                end if


                                E_L_acc = gamma_L_acc*pres_L_acc + pi_inf_L_acc + 5d-1*rho_L_acc*vel_L_rms_acc**2d0

                                E_R_acc = gamma_R_acc*pres_R_acc + pi_inf_R_acc + 5d-1*rho_R_acc*vel_R_rms_acc**2d0

                                H_L_acc = (E_L_acc + pres_L_acc)/rho_L_acc
                                H_R_acc = (E_R_acc + pres_R_acc)/rho_R_acc
                                if(avg_state == 2) then



                                    if (bubbles) then
!$acc loop seq
                                        do i = 1, nb
                                            R0_L_acc(i) =  qL_prim_rs${XYZ}$_vf_flat(j, k, l, rs(i) )
                                            R0_R_acc(i) =  qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, rs(i))

                                            V0_L_acc(i) = qL_prim_rs${XYZ}$_vf_flat(j, k, l, vs(i))
                                            V0_R_acc(i) = qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, vs(i))
                                            if (.not. polytropic) then
                                                P0_L_acc(i) = qL_prim_rs${XYZ}$_vf_flat(j, k, l, ps(i))
                                                P0_R_acc(i) = qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, ps(i))
                                            end if
                                        end do

                                        call s_comp_n_from_prim(alpha_L_acc(num_fluids), R0_L_acc, nbub_L_acc)
                                        call s_comp_n_from_prim(alpha_R_acc(num_fluids), R0_R_acc, nbub_R_acc)
!$acc loop seq
                                        do i = 1, nb
                                            if (.not. qbmm) then
                                                if (polytropic) then
                                                    pbw_L_acc(i) = f_cpbw_KM(R0(i), R0_L_acc(i), V0_L_acc(i), 0d0)
                                                    pbw_R_acc(i) = f_cpbw_KM(R0(i), R0_R_acc(i), V0_R_acc(i), 0d0)
                                                else
                                                    pbw_L_acc(i) = f_cpbw_KM(R0(i), R0_L_acc(i), V0_L_acc(i), P0_L_acc(i))
                                                    pbw_R_acc(i) = f_cpbw_KM(R0(i), R0_R_acc(i), V0_R_acc(i), P0_R_acc(i))
                                                end if
                                            end if
                                        end do

                                        if (qbmm) then
                                            PbwR3Lbar_acc = mom_sp(4)%sf(j, k, l)
                                            PbwR3Rbar_acc = mom_sp(4)%sf(j + 1, k, l)

                                            R3Lbar_acc = mom_sp(1)%sf(j, k, l)
                                            R3Rbar_acc = mom_sp(1)%sf(j + 1, k, l)

                                            R3V2Lbar_acc = mom_sp(3)%sf(j, k, l)
                                            R3V2Rbar_acc = mom_sp(3)%sf(j + 1, k, l)
                                        else
                                            call s_quad(pbw_L_acc*(R0_L_acc**3.d0), PbwR3Lbar_acc)
                                            call s_quad(pbw_R_acc*(R0_R_acc**3.d0), PbwR3Rbar_acc)

                                            call s_quad(R0_L_acc**3.d0, R3Lbar_acc)
                                            call s_quad(R0_R_acc**3.d0, R3Rbar_acc)

                                            call s_quad((R0_L_acc**3.d0)*(V0_L_acc**2.d0), R3V2Lbar_acc)
                                            call s_quad((R0_R_acc**3.d0)*(V0_R_acc**2.d0), R3V2Rbar_acc)
                                        end if

                                        !ptilde = \alf( pl - \bar{ pbw R^3)/\bar{R^3} - rho \bar{R^3 \Rdot^2}/\bar{R^3} )
                                        if (alpha_L_acc(num_fluids) < small_alf .or. R3Lbar_acc < small_alf) then
                                            ptilde_L_acc = alpha_L_acc(num_fluids)*pres_L_acc
                                        else
                                            ptilde_L_acc = alpha_L_acc(num_fluids)*(pres_L_acc - PbwR3Lbar_acc/R3Lbar_acc - &
                                                                            rho_L_acc*R3V2Lbar_acc/R3Lbar_acc)
                                        end if

                                        if (alpha_R_acc(num_fluids) < small_alf .or. R3Rbar_acc < small_alf) then
                                            ptilde_R_acc = alpha_R_acc(num_fluids)*pres_R_acc
                                        else
                                            ptilde_R_acc = alpha_R_acc(num_fluids)*(pres_R_acc - PbwR3Rbar_acc/R3Rbar_acc - &
                                                                            rho_R_acc*R3V2Rbar_acc/R3Rbar_acc)
                                        end if

                                        if ((ptilde_L_acc .ne. ptilde_L_acc) .or. (ptilde_R_acc .ne. ptilde_R_acc)) then
                                        end if

                                        ptil(j, k, l) = 0.5d0*(ptilde_L_acc + ptilde_R_acc)
                                    end if

                                    rho_avg_acc = 5d-1*(rho_L_acc + rho_R_acc)
                !$acc loop seq
                                    do i = 1, num_dims
                                        vel_avg_acc(i) = 5d-1*(vel_L_acc(i) + vel_R_acc(i))
                                    end do

                                    H_avg_acc = 5d-1*(H_L_acc + H_R_acc)

                                    gamma_avg_acc = 5d-1*(gamma_L_acc + gamma_R_acc)

                                elseif(avg_state == 1) then

                                    rho_avg_acc = sqrt(rho_L_acc*rho_R_acc)
                !$acc loop seq
                                    do i = 1, num_dims
                                        vel_avg_acc(i) = (sqrt(rho_L_acc)*vel_L_acc(i) + sqrt(rho_R_acc)*vel_R_acc(i))/ &
                                            (sqrt(rho_L_acc) + sqrt(rho_R_acc))
                                    end do

                                    H_avg_acc = (sqrt(rho_L_acc)*H_L_acc + sqrt(rho_R_acc)*H_R_acc)/ &
                                        (sqrt(rho_L_acc) + sqrt(rho_R_acc))

                                    gamma_avg_acc = (sqrt(rho_L_acc)*gamma_L_acc + sqrt(rho_R_acc)*gamma_R_acc)/ &
                                        (sqrt(rho_L_acc) + sqrt(rho_R_acc))
                                end if

                                vel_avg_rms_acc = 0d0
                !$acc loop seq
                                do i = 1, num_dims
                                    vel_avg_rms_acc = vel_avg_rms_acc + vel_avg_acc(i)**2d0
                                end do
                                vel_avg_rms_acc = sqrt(vel_avg_rms_acc)

                                if (mixture_err) then
                                    if ((H_avg_acc - 5d-1*vel_avg_rms_acc**2d0) < 0d0) then
                                        c_avg_acc = sgm_eps
                                    else

                                        c_avg_acc = sqrt((H_avg_acc - 5d-1*vel_avg_rms_acc**2d0)/gamma_avg_acc)
                                    end if
                                else

                                    c_avg_acc = sqrt((H_avg_acc - 5d-1*vel_avg_rms_acc**2d0)/gamma_avg_acc)
                                end if

                                if (alt_soundspeed) then


                                    blkmod1_acc = ((gammas(1) + 1d0)*pres_L_acc + &
                                            pi_infs(1))/gammas(1)
                                    blkmod2_acc = ((gammas(2) + 1d0)*pres_L_acc + &
                                            pi_infs(2))/gammas(2)
                                    c_L_acc = 1d0/(rho_L_acc*(alpha_L_acc(1)/blkmod1_acc + alpha_L_acc(2)/blkmod2_acc))

                                    blkmod1_acc = ((gammas(1) + 1d0)*pres_R_acc + &
                                            pi_infs(1))/gammas(1)
                                    blkmod2_acc = ((gammas(2) + 1d0)*pres_R_acc + &
                                            pi_infs(2))/gammas(2)
                                    c_R_acc = 1d0/(rho_R_acc*(alpha_R_acc(1)/blkmod1_acc + alpha_R_acc(2)/blkmod2_acc))

                                else
                                    ! Sound speed for bubble mmixture to order O(\alpha)

                                    if (mpp_lim .and. (num_fluids > 1)) then
                                        c_L_acc = (1d0/gamma_L_acc + 1d0)* &
                                            (pres_L_acc + pi_inf_L_acc)/rho_L_acc
                                        c_R_acc = (1d0/gamma_R_acc + 1d0)* &
                                            (pres_R_acc + pi_inf_R_acc)/rho_R_acc
                                    else
                                        c_L_acc = &
                                            (1d0/gamma_L_acc + 1d0)* &
                                            (pres_L_acc + pi_inf_L_acc)/ &
                                            (rho_L_acc*(1d0 - alpha_L_acc(num_fluids)))
                                        c_R_acc = &
                                            (1d0/gamma_R_acc + 1d0)* &
                                            (pres_R_acc + pi_inf_R_acc)/ &
                                            (rho_R_acc*(1d0 - alpha_R_acc(num_fluids)))
                                    end if
                                end if


                                if (mixture_err .and. c_L_acc < 0d0) then
                                    c_L_acc = 100.d0*sgm_eps
                                else
                                    c_L_acc = sqrt(c_L_acc)
                                end if
                                if (mixture_err .and. c_R_acc < 0d0) then
                                    c_R_acc = 100.d0*sgm_eps
                                else
                                    c_R_acc = sqrt(c_R_acc)
                                end if

                                if(wave_speeds == 1) then
                                    s_L_acc = min(vel_L_acc(dir_idx(1)) - c_L_acc, vel_R_acc(dir_idx(1)) - c_R_acc)
                                    s_R_acc = max(vel_R_acc(dir_idx(1)) + c_R_acc, vel_L_acc(dir_idx(1)) + c_L_acc)

                                    s_S_acc = (pres_R_acc - pres_L_acc + rho_L_acc*vel_L_acc(dir_idx(1))* &
                                    (s_L_acc - vel_L_acc(dir_idx(1))) - &
                                    rho_R_acc*vel_R_acc(dir_idx(1))* &
                                    (s_R_acc - vel_R_acc(dir_idx(1)))) &
                                    /(rho_L_acc*(s_L_acc - vel_L_acc(dir_idx(1))) - &
                                        rho_R_acc*(s_R_acc - vel_R_acc(dir_idx(1))))
                                elseif(wave_speeds == 2) then
                                    pres_SL = 5d-1*(pres_L_acc + pres_R_acc+ rho_avg_acc*c_avg_acc* &
                                        (vel_L_acc(dir_idx(1)) - &
                                            vel_R_acc(dir_idx(1))))

                                    pres_SR = pres_SL

                                    Ms_L = max(1d0, sqrt(1d0 + ((5d-1 + gamma_L_acc)/(1d0 + gamma_L_acc))* &
                                                        (pres_SL/pres_L_acc - 1d0)*pres_L_acc/ &
                                                        ((pres_L_acc + pi_inf_L_acc/(1d0 + gamma_L_acc)))))
                                    Ms_R = max(1d0, sqrt(1d0 + ((5d-1 + gamma_R_acc)/(1d0 + gamma_R_acc))* &
                                                        (pres_SR/pres_R_acc - 1d0)*pres_R_acc/ &
                                                        ((pres_R_acc + pi_inf_R_acc/(1d0 + gamma_R_acc)))))

                                    s_L_acc = vel_L_acc(dir_idx(1)) - c_L_acc*Ms_L
                                    s_R_acc = vel_R_acc(dir_idx(1)) + c_R_acc*Ms_R

                                    s_S_acc = 5d-1*((vel_L_acc(dir_idx(1)) + vel_R_acc(dir_idx(1))) + &
                                                (pres_L_acc - pres_R_acc)/ &
                                                            (rho_avg_acc*c_avg_acc))
                                end if




                                ! follows Einfeldt et al.
                                ! s_M/P = min/max(0.,s_L/R)
                                s_M_acc = min(0d0, s_L_acc); s_P_acc = max(0d0, s_R_acc)

                                ! goes with q_star_L/R = xi_L/R * (variable)
                                ! xi_L/R = ( ( s_L/R - u_L/R )/(s_L/R - s_star) )
                                xi_L_acc = (s_L_acc - vel_L_acc(dir_idx(1)))/(s_L_acc - s_S_acc)
                                xi_R_acc = (s_R_acc - vel_R_acc(dir_idx(1)))/(s_R_acc - s_S_acc)

                                ! goes with numerical velocity in x/y/z directions
                                ! xi_P/M = 0.5 +/m sgn(0.5,s_star)
                                xi_M_acc = (5d-1 + sign(5d-1, s_S_acc))
                                xi_P_acc = (5d-1 - sign(5d-1, s_S_acc))

                !$acc loop seq
                                do i = 1, contxe
                                    flux_rs${XYZ}$_vf_flat(j, k, l, i) = &
                                        xi_M_acc*alpha_rho_L_acc(i) &
                                        *(vel_L_acc(dir_idx(1)) + s_M_acc*(xi_L_acc - 1d0)) &
                                        + xi_P_acc*alpha_rho_R_acc(i) &
                                        *(vel_R_acc(dir_idx(1)) + s_P_acc*(xi_R_acc - 1d0))
                                end do

                                if (bubbles  .and. (num_fluids > 1)) then
                                    ! Kill mass transport @ gas density
                                    flux_rs${XYZ}$_vf_flat(j, k, l, contxe) = 0.d0
                                end if

                                ! Momentum flux.
                                ! f = \rho u u + p I, q = \rho u, q_star = \xi * \rho*(s_star, v, w)

                                    ! Include p_tilde
                !$acc loop seq
                                do i = 1, num_dims
                                    flux_rs${XYZ}$_vf_flat(j, k, l, contxe + dir_idx(i)) = &
                                        xi_M_acc*(rho_L_acc*(vel_L_acc(dir_idx(1))* &
                                                    vel_L_acc(dir_idx(i)) + &
                                                    s_M_acc*(xi_L_acc*(dir_flg(dir_idx(i))*s_S_acc + &
                                                                (1d0 - dir_flg(dir_idx(i)))* &
                                                                vel_L_acc(dir_idx(i))) - vel_L_acc(dir_idx(i)))) + &
                                            dir_flg(dir_idx(i))*(pres_L_acc - ptilde_L_acc)) &
                                        + xi_P_acc*(rho_R_acc*(vel_R_acc(dir_idx(1))* &
                                                    vel_R_acc(dir_idx(i)) + &
                                                    s_P_acc*(xi_R_acc*(dir_flg(dir_idx(i))*s_S_acc + &
                                                                (1d0 - dir_flg(dir_idx(i)))* &
                                                                vel_R_acc(dir_idx(i))) - vel_R_acc(dir_idx(i)))) + &
                                                dir_flg(dir_idx(i))*(pres_R_acc - ptilde_R_acc))
                                    ! if (j==0) print*, 'flux_rs_vf', flux_rs_vf(cont_idx%end+dir_idx(i))%sf(j,k,l)
                                end do



                                ! Energy flux.
                                ! f = u*(E+p), q = E, q_star = \xi*E+(s-u)(\rho s_star + p/(s-u))

                                flux_rs${XYZ}$_vf_flat(j, k, l, E_idx) = &
                                    xi_M_acc*(vel_L_acc(dir_idx(1))*(E_L_acc + pres_L_acc - ptilde_L_acc) + &
                                        s_M_acc*(xi_L_acc*(E_L_acc + (s_S_acc - vel_L_acc(dir_idx(1)))* &
                                                    (rho_L_acc*s_S_acc + (pres_L_acc - ptilde_L_acc)/ &
                                                    (s_L_acc - vel_L_acc(dir_idx(1))))) - E_L_acc)) &
                                    + xi_P_acc*(vel_R_acc(dir_idx(1))*(E_R_acc + pres_R_acc - ptilde_R_acc) + &
                                            s_P_acc*(xi_R_acc*(E_R_acc + (s_S_acc - vel_R_acc(dir_idx(1)))* &
                                                    (rho_R_acc*s_S_acc + (pres_R_acc - ptilde_R_acc)/ &
                                                        (s_R_acc - vel_R_acc(dir_idx(1))))) - E_R_acc))


                                ! Volume fraction flux

                !$acc loop seq
                                do i = advxb, advxe
                                    flux_rs${XYZ}$_vf_flat(j, k, l, i) = &
                                        xi_M_acc*qL_prim_rs${XYZ}$_vf_flat(j, k, l, i) &
                                        *(vel_L_acc(dir_idx(1)) + s_M_acc*(xi_L_acc - 1d0)) &
                                        + xi_P_acc*qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, i) &
                                        *(vel_R_acc(dir_idx(1)) + s_P_acc*(xi_R_acc - 1d0))
                                end do


                                ! Source for volume fraction advection equation
                !$acc loop seq
                                do i = 1, num_dims
                                    vel_src_rs${XYZ}$_vf_flat(j, k, l, dir_idx(i)) = &
                                        xi_M_acc*(vel_L_acc(dir_idx(i)) + &
                                            dir_flg(dir_idx(i))* &
                                            s_M_acc*(xi_L_acc - 1d0)) &
                                        + xi_P_acc*(vel_R_acc(dir_idx(i)) + &
                                                dir_flg(dir_idx(i))* &
                                                s_P_acc*(xi_R_acc - 1d0))

                                    !IF ( (model_eqns == 4) .or. (num_fluids==1) ) vel_src_rs_vf(dir_idx(i))%sf(j,k,l) = 0d0
                                end do

                                flux_src_rs${XYZ}$_vf_flat(j, k, l, advxb) = vel_src_rs${XYZ}$_vf_flat(j, k, l, dir_idx(1))

                                ! Add advection flux for bubble variables

                !$acc loop seq
                                do i = bubxb, bubxe
                                    flux_rs${XYZ}$_vf_flat(j, k, l, i) = &
                                        xi_M_acc*nbub_L_acc*qL_prim_rs${XYZ}$_vf_flat(j, k, l, i) &
                                        *(vel_L_acc(dir_idx(1)) + s_M_acc*(xi_L_acc - 1d0)) &
                                        + xi_P_acc*nbub_R_acc*qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, i) &
                                        *(vel_R_acc(dir_idx(1)) + s_P_acc*(xi_R_acc - 1d0))
                                end do



                                ! Geometrical source flux for cylindrical coordinates

#:if (NORM_DIR == 2)
                                if (cyl_coord) then
                                    ! Substituting the advective flux into the inviscid geometrical source flux
                !$acc loop seq
                                    do i = 1, E_idx
                                        flux_gsrc_rs${XYZ}$_vf_flat(j, k, l, i) = flux_rs${XYZ}$_vf_flat(j, k, l, i)
                                    end do
                                    ! Recalculating the radial momentum geometric source flux
                                    flux_gsrc_rs${XYZ}$_vf_flat(j, k, l, contxe + dir_idx(1)) = &
                                        xi_M_acc*(rho_L_acc*(vel_L_acc(dir_idx(1))* &
                                                    vel_L_acc(dir_idx(1)) + &
                                                    s_M_acc*(xi_L_acc*(dir_flg(dir_idx(1))*s_S_acc + &
                                                                (1d0 - dir_flg(dir_idx(1)))* &
                                                                vel_L_acc(dir_idx(1))) - vel_L_acc(dir_idx(1))))) &
                                        + xi_P_acc*(rho_R_acc*(vel_R_acc(dir_idx(1))* &
                                                    vel_R_acc(dir_idx(1)) + &
                                                    s_P_acc*(xi_R_acc*(dir_flg(dir_idx(1))*s_S_acc + &
                                                                (1d0 - dir_flg(dir_idx(1)))* &
                                                                vel_R_acc(dir_idx(1))) - vel_R_acc(dir_idx(1)))))
                                    ! Geometrical source of the void fraction(s) is zero
                !$acc loop seq
                                    do i = advxb, advxe
                                        flux_gsrc_rs${XYZ}$_vf_flat(j, k, l, i) = 0d0
                                    end do
                                end if
#:endif
#:if (NORM_DIR == 3)
                                if (grid_geometry == 3) then
                                    !$acc loop seq 
                                    do i = 1, sys_size
                                        flux_gsrc_rs${XYZ}$_vf_flat(j, k, l, i) = 0d0
                                    end do
                                    
                                    flux_gsrc_rs${XYZ}$_vf_flat(j, k, l, momxb + 1) = &
                                        -xi_M_acc*(rho_L_acc*(vel_L_acc(dir_idx(1))* &
                                                        vel_L_acc(dir_idx(1)) + &
                                                        s_M_acc*(xi_L_acc*(dir_flg(dir_idx(1))*s_S_acc + &
                                                                    (1d0 - dir_flg(dir_idx(1)))* &
                                                                    vel_L_acc(dir_idx(1))) - vel_L_acc(dir_idx(1))))) &
                                        - xi_P_acc*(rho_R_acc*(vel_R_acc(dir_idx(1))* &
                                                        vel_R_acc(dir_idx(1)) + &
                                                        s_P_acc*(xi_R_acc*(dir_flg(dir_idx(1))*s_S_acc + &
                                                                    (1d0 - dir_flg(dir_idx(1)))* &
                                                                    vel_R_acc(dir_idx(1))) - vel_R_acc(dir_idx(1)))))
                                    flux_gsrc_rs${XYZ}$_vf_flat(j, k, l, momxe) = flux_rs${XYZ}$_vf_flat(j, k, l, momxb + 1)


                                end if
#:endif
                            end do
                        end do
                    end do
                !$acc end parallel loop
                else
        !$acc parallel loop collapse(3) gang vector default(present) private(vel_L_acc, vel_R_acc)        
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end
                                idx1 = 1; if (dir_idx(1).eq.2) idx1 = 2; if (dir_idx(1).eq.3) idx1 = 3

                                vel_L_rms_acc = 0d0; vel_R_rms_acc = 0d0
        !$acc loop seq
                                do i = 1, num_dims
                                    vel_L_acc(i) = qL_prim_rs${XYZ}$_vf_flat(j, k, l, contxe + i)
                                    vel_R_acc(i) = qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, contxe + i)
                                    vel_L_rms_acc = vel_L_rms_acc + vel_L_acc(i)**2d0
                                    vel_R_rms_acc = vel_R_rms_acc + vel_R_acc(i)**2d0
                                end do
                                vel_L_rms_acc = sqrt(vel_L_rms_acc)
                                vel_R_rms_acc = sqrt(vel_R_rms_acc)

                                pres_L_acc = qL_prim_rs${XYZ}$_vf_flat(j, k, l, E_idx)
                                pres_R_acc = qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, E_idx)

                                rho_L_acc = 0d0
                                gamma_L_acc = 0d0
                                pi_inf_L_acc = 0d0
                                rho_R_acc = 0d0
                                gamma_R_acc = 0d0
                                pi_inf_R_acc = 0d0
        !$acc loop seq 
                                do i = 1, num_fluids
                                    rho_L_acc = rho_L_acc + qL_prim_rs${XYZ}$_vf_flat(j, k, l, i)
                                    gamma_L_acc = gamma_L_acc + qL_prim_rs${XYZ}$_vf_flat(j, k, l, E_idx + i)*gammas(i)
                                    pi_inf_L_acc = pi_inf_L_acc + qL_prim_rs${XYZ}$_vf_flat(j, k, l, E_idx + i)*pi_infs(i)

                                    rho_R_acc = rho_R_acc + qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, i)
                                    gamma_R_acc = gamma_R_acc + qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, E_idx + i)*gammas(i)
                                    pi_inf_R_acc = pi_inf_R_acc + qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, E_idx + i)*pi_infs(i)
                                end do        

                                E_L_acc = gamma_L_acc*pres_L_acc + pi_inf_L_acc + 5d-1*rho_L_acc*vel_L_rms_acc**2d0

                                E_R_acc = gamma_R_acc*pres_R_acc + pi_inf_R_acc + 5d-1*rho_R_acc*vel_R_rms_acc**2d0

                                H_L_acc = (E_L_acc + pres_L_acc)/rho_L_acc
                                H_R_acc = (E_R_acc + pres_R_acc)/rho_R_acc
                                if(avg_state == 2) then

                                    rho_avg_acc = 5d-1*(rho_L_acc + rho_R_acc)
                                    vel_avg_rms_acc = (5d-1*(vel_L_acc(1) + vel_R_acc(1)))**2d0
                                    if (num_dims.ge.2) then
                                        vel_avg_rms_acc = vel_avg_rms_acc + (5d-1*(vel_L_acc(2) + vel_R_acc(2)))**2d0
                                    end if
                                    if (num_dims.eq.3)  then
                                        vel_avg_rms_acc = vel_avg_rms_acc + (5d-1*(vel_L_acc(3) + vel_R_acc(3)))**2d0
                                    end if

                                    H_avg_acc = 5d-1*(H_L_acc + H_R_acc)

                                    gamma_avg_acc = 5d-1*(gamma_L_acc + gamma_R_acc)

                                elseif(avg_state == 1) then

                                    rho_avg_acc = sqrt(rho_L_acc*rho_R_acc)
                                    vel_avg_rms_acc = (sqrt(rho_L_acc)*vel_L_acc(1) + sqrt(rho_R_acc)*vel_R_acc(1))**2d0/ &
                                            (sqrt(rho_L_acc) + sqrt(rho_R_acc))**2d0

                                    if (num_dims.ge.2) then
                                    vel_avg_rms_acc = vel_avg_rms_acc + (sqrt(rho_L_acc)*vel_L_acc(2) + sqrt(rho_R_acc)*vel_R_acc(2))**2d0/ &
                                            (sqrt(rho_L_acc) + sqrt(rho_R_acc))**2d0
                                    end if
                                    if (num_dims.eq.3) then
                                    vel_avg_rms_acc = vel_avg_rms_acc + (sqrt(rho_L_acc)*vel_L_acc(3) + sqrt(rho_R_acc)*vel_R_acc(3))**2d0/ &
                                            (sqrt(rho_L_acc) + sqrt(rho_R_acc))**2d0
                                    end if
                                        
                                    H_avg_acc = (sqrt(rho_L_acc)*H_L_acc + sqrt(rho_R_acc)*H_R_acc)/ &
                                        (sqrt(rho_L_acc) + sqrt(rho_R_acc))

                                    gamma_avg_acc = (sqrt(rho_L_acc)*gamma_L_acc + sqrt(rho_R_acc)*gamma_R_acc)/ &
                                        (sqrt(rho_L_acc) + sqrt(rho_R_acc))
                                end if
                                vel_avg_rms_acc = sqrt(vel_avg_rms_acc)

                                if (mixture_err) then
                                    if ((H_avg_acc - 5d-1*vel_avg_rms_acc**2d0) < 0d0) then
                                        c_avg_acc = sgm_eps
                                    else

                                        c_avg_acc = sqrt((H_avg_acc - 5d-1*vel_avg_rms_acc**2d0)/gamma_avg_acc)
                                    end if
                                else

                                    c_avg_acc = sqrt((H_avg_acc - 5d-1*vel_avg_rms_acc**2d0)/gamma_avg_acc)
                                end if

                                if (alt_soundspeed) then


                                    blkmod1_acc = ((gammas(1) + 1d0)*pres_L_acc + &
                                                pi_infs(1))/gammas(1)
                                    blkmod2_acc = ((gammas(2) + 1d0)*pres_L_acc + &
                                                pi_infs(2))/gammas(2)
                                    c_L_acc = 1d0/(rho_L_acc*(qL_prim_rs${XYZ}$_vf_flat(j, k, l, E_idx + 1)/blkmod1_acc &
                                                            + qL_prim_rs${XYZ}$_vf_flat(j, k, l, E_idx + 2)/blkmod2_acc))

                                    blkmod1_acc = ((gammas(1) + 1d0)*pres_R_acc + &
                                                pi_infs(1))/gammas(1)
                                    blkmod2_acc = ((gammas(2) + 1d0)*pres_R_acc + &
                                                pi_infs(2))/gammas(2)
                                    c_R_acc = 1d0/(rho_R_acc*(qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, E_idx + 1)/blkmod1_acc &
                                                            + qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, e_idx + 2)/blkmod2_acc))

                                else
                                    c_L_acc = ((H_L_acc - 5d-1*vel_L_rms_acc**2d0)/gamma_L_acc)

                                    c_R_acc = ((H_R_acc - 5d-1*vel_R_rms_acc**2d0)/gamma_R_acc)
                                end if
                                    
                                if (mixture_err .and. c_L_acc < 0d0) then
                                    c_L_acc = 100.d0*sgm_eps
                                else
                                    c_L_acc = sqrt(c_L_acc)
                                end if
                                if (mixture_err .and. c_R_acc < 0d0) then
                                    c_R_acc = 100.d0*sgm_eps
                                else
                                    c_R_acc = sqrt(c_R_acc)
                                end if

                                if(wave_speeds == 1) then
                                    s_L_acc = min(vel_L_acc(idx1) - c_L_acc, vel_R_acc(idx1) - c_R_acc)
                                    s_R_acc = max(vel_R_acc(idx1) + c_R_acc, vel_L_acc(idx1) + c_L_acc)

                                    s_S_acc = (pres_R_acc - pres_L_acc + rho_L_acc*vel_L_acc(idx1)* &
                                        (s_L_acc - vel_L_acc(idx1)) - &
                                        rho_R_acc*vel_R_acc(idx1)* &
                                        (s_R_acc - vel_R_acc(idx1))) &
                                        /(rho_L_acc*(s_L_acc - vel_L_acc(idx1)) - &
                                        rho_R_acc*(s_R_acc - vel_R_acc(idx1)))
                                elseif(wave_speeds == 2) then
                                    pres_SL = 5d-1*(pres_L_acc + pres_R_acc+ rho_avg_acc*c_avg_acc* &
                                        (vel_L_acc(idx1) - &
                                            vel_R_acc(idx1)))

                                    pres_SR = pres_SL

                                    Ms_L = max(1d0, sqrt(1d0 + ((5d-1 + gamma_L_acc)/(1d0 + gamma_L_acc))* &
                                                        (pres_SL/pres_L_acc - 1d0)*pres_L_acc/ &
                                                        ((pres_L_acc + pi_inf_L_acc/(1d0 + gamma_L_acc)))))
                                    Ms_R = max(1d0, sqrt(1d0 + ((5d-1 + gamma_R_acc)/(1d0 + gamma_R_acc))* &
                                                        (pres_SR/pres_R_acc - 1d0)*pres_R_acc/ &
                                                        ((pres_R_acc + pi_inf_R_acc/(1d0 + gamma_R_acc)))))

                                    s_L_acc = vel_L_acc(idx1) - c_L_acc*Ms_L
                                    s_R_acc = vel_R_acc(idx1) + c_R_acc*Ms_R

                                    s_S_acc = 5d-1*((vel_L_acc(idx1) + vel_R_acc(idx1)) + &
                                                (pres_L_acc - pres_R_acc)/ &
                                                            (rho_avg_acc*c_avg_acc))
                                end if




                                ! follows Einfeldt et al.
                                ! s_M/P = min/max(0.,s_L/R)
                                s_M_acc = min(0d0, s_L_acc); s_P_acc = max(0d0, s_R_acc)

                                ! goes with q_star_L/R = xi_L/R * (variable)
                                ! xi_L/R = ( ( s_L/R - u_L/R )/(s_L/R - s_star) )
                                xi_L_acc = (s_L_acc - vel_L_acc(idx1))/(s_L_acc - s_S_acc)
                                xi_R_acc = (s_R_acc - vel_R_acc(idx1))/(s_R_acc - s_S_acc)

                                ! goes with numerical velocity in x/y/z directions
                                ! xi_P/M = 0.5 +/m sgn(0.5,s_star)
                                xi_M_acc = (5d-1 + sign(5d-1, s_S_acc))
                                xi_P_acc = (5d-1 - sign(5d-1, s_S_acc))

    !$acc loop seq 
                                do i = 1, contxe
                                    flux_rs${XYZ}$_vf_flat(j, k, l, i) = &
                                        xi_M_acc*qL_prim_rs${XYZ}$_vf_flat(j, k, l, i) &
                                        *(vel_L_acc(idx1) + s_M_acc*(xi_L_acc - 1d0)) &
                                        + xi_P_acc*qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, i) &
                                        *(vel_R_acc(idx1) + s_P_acc*(xi_R_acc - 1d0))
                                end do


                                ! Momentum flux.
                                ! f = \rho u u + p I, q = \rho u, q_star = \xi * \rho*(s_star, v, w)

    !$acc loop seq 
                                do i = 1, num_dims
                                    idxi = dir_idx(i)
                                    flux_rs${XYZ}$_vf_flat(j, k, l, contxe + idxi) = &
                                        xi_M_acc*(rho_L_acc*(vel_L_acc(idx1)* &
                                                    vel_L_acc(idxi) + &
                                                    s_M_acc*(xi_L_acc*(dir_flg(idxi)*s_S_acc + &
                                                                (1d0 - dir_flg(idxi))* &
                                                                vel_L_acc(idxi)) - vel_L_acc(idxi))) + &
                                                dir_flg(idxi)*(pres_L_acc)) &
                                        + xi_P_acc*(rho_R_acc*(vel_R_acc(idx1)* &
                                                        vel_R_acc(idxi) + &
                                                        s_P_acc*(xi_R_acc*(dir_flg(idxi)*s_S_acc + &
                                                                    (1d0 - dir_flg(idxi))* &
                                                                    vel_R_acc(idxi)) - vel_R_acc(idxi))) + &
                                                dir_flg(idxi)*(pres_R_acc))
                                    ! if (j==0) print*, 'flux_rs_vf', flux_rs_vf(cont_idx%end+dir_idx(i))%sf(j,k,l)
                                end do


                                ! Energy flux.
                                ! f = u*(E+p), q = E, q_star = \xi*E+(s-u)(\rho s_star + p/(s-u))

                                flux_rs${XYZ}$_vf_flat(j, k, l, E_idx) = &
                                    xi_M_acc*(vel_L_acc(idx1)*(E_L_acc + pres_L_acc) + &
                                            s_M_acc*(xi_L_acc*(E_L_acc + (s_S_acc - vel_L_acc(idx1))* &
                                                    (rho_L_acc*s_S_acc + pres_L_acc/ &
                                                        (s_L_acc - vel_L_acc(idx1)))) - E_L_acc)) &
                                    + xi_P_acc*(vel_R_acc(idx1)*(E_R_acc + pres_R_acc) + &
                                            s_P_acc*(xi_R_acc*(E_R_acc + (s_S_acc - vel_R_acc(idx1))* &
                                                        (rho_R_acc*s_S_acc + pres_R_acc/ &
                                                        (s_R_acc - vel_R_acc(idx1)))) - E_R_acc))


                                ! Volume fraction flux

    !$acc loop seq 
                                do i = advxb, advxe
                                    flux_rs${XYZ}$_vf_flat(j, k, l, i) = &
                                        xi_M_acc*qL_prim_rs${XYZ}$_vf_flat(j, k, l, i) &
                                        *(vel_L_acc(idx1) + s_M_acc*(xi_L_acc - 1d0)) &
                                        + xi_P_acc*qR_prim_rs${XYZ}$_vf_flat(j + 1, k, l, i) &
                                        *(vel_R_acc(idx1) + s_P_acc*(xi_R_acc - 1d0))
                                end do
                                

                                ! Source for volume fraction advection equation
    !$acc loop seq 
                                do i = 1, num_dims
                                    idxi = dir_idx(i)
                                    vel_src_rs${XYZ}$_vf_flat(j, k, l, idxi) = &
                                        xi_M_acc*(vel_L_acc(idxi) + &
                                                dir_flg(idxi)* &
                                                s_M_acc*(xi_L_acc - 1d0)) &
                                        + xi_P_acc*(vel_R_acc(idxi) + &
                                                dir_flg(idxi)* &
                                                s_P_acc*(xi_R_acc - 1d0))

                                    !IF ( (model_eqns == 4) .or. (num_fluids==1) ) vel_src_rs_vf(dir_idx(i))%sf(j,k,l) = 0d0
                                end do

                                flux_src_rs${XYZ}$_vf_flat(j, k, l, advxb) = vel_src_rs${XYZ}$_vf_flat(j, k, l, idx1)


                                ! Geometrical source flux for cylindrical coordinates

#:if (NORM_DIR == 2)
                                if (cyl_coord) then
                                    !Substituting the advective flux into the inviscid geometrical source flux
    !$acc loop seq 
                                    do i = 1, E_idx
                                        flux_gsrc_rs${XYZ}$_vf_flat(j, k, l, i) = flux_rs${XYZ}$_vf_flat(j, k, l, i)
                                    end do
                                    ! Recalculating the radial momentum geometric source flux
                                    flux_gsrc_rs${XYZ}$_vf_flat(j, k, l, contxe + idx1) = &
                                        xi_M_acc*(rho_L_acc*(vel_L_acc(idx1)* &
                                                    vel_L_acc(idx1) + &
                                                    s_M_acc*(xi_L_acc*(dir_flg(idx1)*s_S_acc + &
                                                                (1d0 - dir_flg(idx1))* &
                                                                vel_L_acc(idx1)) - vel_L_acc(idx1)))) &
                                        + xi_P_acc*(rho_R_acc*(vel_R_acc(idx1)* &
                                                        vel_R_acc(idx1) + &
                                                        s_P_acc*(xi_R_acc*(dir_flg(idx1)*s_S_acc + &
                                                                    (1d0 - dir_flg(idx1))* &
                                                                    vel_R_acc(idx1)) - vel_R_acc(idx1))))
                                    ! Geometrical source of the void fraction(s) is zero
    !$acc loop seq 
                                    do i = advxb, advxe
                                        flux_gsrc_rs${XYZ}$_vf_flat(j, k, l, i) = 0d0
                                    end do
                                end if
#:endif
#:if (NORM_DIR == 3)                
                                if (grid_geometry == 3) then
                                    !$acc loop seq 
                                    do i = 1, sys_size
                                        flux_gsrc_rs${XYZ}$_vf_flat(j, k, l, i) = 0d0
                                    end do

                                    flux_gsrc_rs${XYZ}$_vf_flat(j, k, l, momxb + 1) = &
                                        -xi_M_acc*(rho_L_acc*(vel_L_acc(idx1)* &
                                                        vel_L_acc(idx1) + &
                                                        s_M_acc*(xi_L_acc*(dir_flg(idx1)*s_S_acc + &
                                                                    (1d0 - dir_flg(idx1))* &
                                                                    vel_L_acc(idx1)) - vel_L_acc(idx1)))) &
                                        - xi_P_acc*(rho_R_acc*(vel_R_acc(idx1)* &
                                                        vel_R_acc(idx1) + &
                                                        s_P_acc*(xi_R_acc*(dir_flg(idx1)*s_S_acc + &
                                                                    (1d0 - dir_flg(idx1))* &
                                                                    vel_R_acc(idx1)) - vel_R_acc(idx1))))
                                    flux_gsrc_rs${XYZ}$_vf_flat(j, k, l, momxe) = flux_rs${XYZ}$_vf_flat(j, k, l, momxb + 1)

    
                                end if
#:endif
                            end do
                        end do
                    end do 
            end if
        end if
#:endfor
        ! Computing HLLC flux and source flux for Euler system of equations



        ! print*, 'xbounds are: ', is1%beg, is1%end
        ! print*, 'ybounds are: ', is2%beg, is2%end
        ! print*, 'zbounds are: ', is3%beg, is3%end


                    ! print*, 'about to get average state'

            call s_finalize_riemann_solver(flux_vf, flux_src_vf, &
                                       flux_gsrc_vf, &
                                       norm_dir, ix, iy, iz)

        end subroutine s_hllc_riemann_solver_acc




    !>  This procedure is the implementation of the exact Riemann
        !!      solver, see Toro (1999). The effects of viscosity and the
        !!      surface tension have been incorporated following the work
        !!      of Perigaud and Saurel (2005).
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
        !!  @param q_prim_vf Cell-averaged primitive variables
        !!  @param flux_vf Intra-cell fluxes
        !!  @param flux_src_vf  Intra-cell fluxes sources
        !!  @param flux_gsrc_vf Intra-cell geometric fluxes sources
        !!  @param norm_dir Dir. splitting direction
        !!  @param ix Index bounds in the x-dir
        !!  @param iy Index bounds in the y-dir
        !!  @param iz Index bounds in the z-dir
    subroutine s_exact_riemann_solver(qL_prim_rsx_vf_flat, qL_prim_rsy_vf_flat, qL_prim_rsz_vf_flat, dqL_prim_dx_vf, & ! -----
                                      dqL_prim_dy_vf, &
                                      dqL_prim_dz_vf, &
                                      gm_alphaL_vf, &
                                      qR_prim_rsx_vf_flat, qR_prim_rsy_vf_flat, qR_prim_rsz_vf_flat, dqR_prim_dx_vf, &
                                      dqR_prim_dy_vf, &
                                      dqR_prim_dz_vf, &
                                      gm_alphaR_vf, &
                                      q_prim_vf, &
                                      flux_vf, flux_src_vf, &
                                      flux_gsrc_vf, &
                                      norm_dir, ix, iy, iz)

        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:), intent(INOUT) :: qL_prim_rsx_vf_flat, qL_prim_rsy_vf_flat, qL_prim_rsz_vf_flat, qR_prim_rsx_vf_flat, qR_prim_rsy_vf_flat, qR_prim_rsz_vf_flat
        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf

        type(scalar_field), &
            allocatable, dimension(:), &
            intent(INOUT) :: dqL_prim_dx_vf, dqR_prim_dx_vf, &
                             dqL_prim_dy_vf, dqR_prim_dy_vf, &
                             dqL_prim_dz_vf, dqR_prim_dz_vf, &
                             gm_alphaL_vf, gm_alphaR_vf

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: flux_vf, flux_src_vf, flux_gsrc_vf

        integer, intent(IN) :: norm_dir
        type(bounds_info), intent(IN) :: ix, iy, iz

        integer :: i, j, k, l !< Generic loop iterators

        ! Populating the buffers of the left and right Riemann problem
        ! states variables, based on the choice of boundary conditions
        call s_populate_riemann_states_variables_buffers( &
            qL_prim_rsx_vf_flat, qL_prim_rsy_vf_flat, qL_prim_rsz_vf_flat, dqL_prim_dx_vf, &
            dqL_prim_dy_vf, &
            dqL_prim_dz_vf, &
            gm_alphaL_vf, &
            qR_prim_rsx_vf_flat, qR_prim_rsy_vf_flat, qR_prim_rsz_vf_flat, dqR_prim_dx_vf, &
            dqR_prim_dy_vf, &
            dqR_prim_dz_vf, &
            gm_alphaR_vf, &
            norm_dir, ix, iy, iz)

        ! Reshaping inputted data based on dimensional splitting direction
        call s_initialize_riemann_solver( &
                                         q_prim_vf, &
                                         flux_vf, flux_src_vf, &
                                         flux_gsrc_vf, &
                                         norm_dir, ix, iy, iz)

        ! Computing exact flux and source flux for Euler system of equations
        do l = is3%beg, is3%end
            do k = is2%beg, is2%end
                do j = is1%beg, is1%end

                    call s_compute_constant_states(j, k, l)

                    ! Check for pressure positivity condition
                    if ((G4_L*c_L + G4_R*c_R) < (vel_R(dir_idx(1)) - vel_L(dir_idx(1)))) then
                        print '(A)', 'Vacuum is generated by Riemann data. Exiting...'
                        call s_mpi_abort()
                    end if

                    call s_compute_star_region()

                    call s_compute_intercell_solution()

                    do i = 1, cont_idx%end
                        flux_rs_vf(i)%sf(j, k, l) = alpha_rho_IC(i)*vel_IC(dir_idx(1))
                    end do

                    do i = 1, num_dims
                        flux_rs_vf(cont_idx%end + dir_idx(i))%sf(j, k, l) = &
                            rho_IC*vel_IC(dir_idx(1))*vel_IC(dir_idx(i)) + dir_flg(dir_idx(i))*pres_IC
                    end do

                    flux_rs_vf(E_idx)%sf(j, k, l) = vel_IC(dir_idx(1))*(E_IC + pres_IC)

                    do i = 1, adv_idx%end - E_idx
                        flux_rs_vf(E_idx + i)%sf(j, k, l) = alpha_IC(i)*vel_IC(dir_idx(1))
                    end do

                    do i = 1, num_dims
                        vel_src_rs_vf(dir_idx(i))%sf(j, k, l) = vel_IC(dir_idx(i))
                    end do
                end do
            end do
        end do



        ! Reshaping outputted data based on dimensional splitting direction
        call s_finalize_riemann_solver(flux_vf, flux_src_vf, &
                                       flux_gsrc_vf, &
                                       norm_dir, ix, iy, iz)

    end subroutine s_exact_riemann_solver ! --------------------------------

    !>  The procedure assigns and computes the left and right
        !!      states of the Riemann problem
        !! @param j  First coordinate index
        !! @param k Second coordinate index
        !! @param l  Third coordinate index
    subroutine s_compute_constant_states(j, k, l) ! --------------------------

        integer, intent(IN) :: j, k, l

        integer :: i !< Generic loop iterator

        ! Left and Right Riemann States
        do i = 1, cont_idx%end
            alpha_rho_L(i) = qL_prim_rs_vf(i)%sf(j, k, l)
            alpha_rho_R(i) = qR_prim_rs_vf(i)%sf(j + 1, k, l)
        end do

        do i = 1, num_dims
            vel_L(i) = qL_prim_rs_vf(cont_idx%end + i)%sf(j, k, l)
            vel_R(i) = qR_prim_rs_vf(cont_idx%end + i)%sf(j + 1, k, l)
        end do

        pres_L = qL_prim_rs_vf(E_idx)%sf(j, k, l)
        pres_R = qR_prim_rs_vf(E_idx)%sf(j + 1, k, l)

        call s_convert_to_mixture_variables(qL_prim_rs_vf, &
                                            rho_L, gamma_L, &
                                            pi_inf_L, Re_L, &
                                            j, k, l)
        call s_convert_to_mixture_variables(qR_prim_rs_vf, &
                                            rho_R, gamma_R, &
                                            pi_inf_R, Re_R, &
                                             j + 1, k, l)

        E_L = gamma_L*pres_L + pi_inf_L + 5d-1*rho_L*sum(vel_L**2d0)
        E_R = gamma_R*pres_R + pi_inf_R + 5d-1*rho_R*sum(vel_R**2d0)

        H_L = (E_L + pres_L)/rho_L
        H_R = (E_R + pres_R)/rho_R

        call s_compute_mixture_sound_speeds(qL_prim_rs_vf, qR_prim_rs_vf,j, k, l)

        do i = 1, 2
            if (Re_size(i) > 0) then
                Re_avg_rs_vf(i)%sf(j, k, l) = 2d0/(1d0/Re_L(i) + 1d0/Re_R(i))
            end if
        end do

        ! Compute gamma-related constants
        G1_L = 1d0/(2d0*(gamma_L + 1d0))
        G1_R = 1d0/(2d0*(gamma_R + 1d0))
        G2_L = (1d0 + 2d0*gamma_L)/(2d0*(gamma_L + 1d0))
        G2_R = (1d0 + 2d0*gamma_R)/(2d0*(gamma_R + 1d0))
        G3_L = 2d0*(gamma_L + 1d0)
        G3_R = 2d0*(gamma_R + 1d0)
        G4_L = 2d0*gamma_L
        G4_R = 2d0*gamma_R
        G5_L = 2d0/((1d0/gamma_L) + 2d0)
        G5_R = 2d0/((1d0/gamma_R) + 2d0)
        G6_L = 1d0/(1d0 + 2d0*gamma_L)
        G6_R = 1d0/(1d0 + 2d0*gamma_R)
        G7_L = 1d0/(2d0*gamma_L)
        G7_R = 1d0/(2d0*gamma_R)
        G8_L = 1d0/gamma_L
        G8_R = 1d0/gamma_R

        ! Surface tension pressure contribution
        dpres_L = 0d0
        dpres_R = 0d0

        dpres_L = 5d-1
        dpres_R = -5d-1

    end subroutine s_compute_constant_states ! -----------------------------

    !> Compute mixture sound speed
        !! @param j  First coordinate index
        !! @param k Second coordinate index
        !! @param l  Third coordinate index
    subroutine s_compute_mixture_sound_speeds(qL_prim_rs_vf, qR_prim_rs_vf,j, k, l) ! ---------------------
!$acc routine seq
        type(scalar_field), dimension(sys_size), intent(IN) :: qL_prim_rs_vf, qR_prim_rs_vf
        integer, intent(IN) :: j, k, l

        real(kind(0d0)) :: blkmod1, blkmod2 !< Fluid bulk modulus for alternate sound speed

        integer :: i !< Generic loop iterator

        if (alt_soundspeed) then


            blkmod1 = ((gammas(1) + 1d0)*pres_L + &
                       pi_infs(1))/gammas(1)
            blkmod2 = ((gammas(2) + 1d0)*pres_L + &
                       pi_infs(2))/gammas(2)
            c_L = 1d0/(rho_L*(alpha_L(1)/blkmod1 + alpha_L(2)/blkmod2))

            blkmod1 = ((gammas(1) + 1d0)*pres_R + &
                       pi_infs(1))/gammas(1)
            blkmod2 = ((gammas(2) + 1d0)*pres_R + &
                       pi_infs(2))/gammas(2)
            c_R = 1d0/(rho_R*(alpha_R(1)/blkmod1 + alpha_R(2)/blkmod2))

        elseif (model_eqns == 3) then
            c_L = 0d0
            c_R = 0d0
            do i = 1, num_fluids
                c_L = c_L + qL_prim_rs_vf(i + advxb - 1)%sf(j, k, l)*(1d0/gammas(i) + 1d0)* &
                      (qL_prim_rs_vf(E_idx)%sf(j, k, l) + pi_infs(i)/(gammas(i) + 1d0))
                c_R = c_R + qR_prim_rs_vf(i + advxb - 1)%sf(j + 1, k, l)*(1d0/gammas(i) + 1d0)* &
                      (qR_prim_rs_vf(E_idx)%sf(j + 1, k, l) + pi_infs(i)/(gammas(i) + 1d0))
            end do
            c_L = c_L/rho_L
            c_R = c_R/rho_R
        elseif ((model_eqns == 4) .or. (model_eqns == 2 .and. bubbles)) then
            ! Sound speed for bubble mmixture to order O(\alpha)

            if (mpp_lim .and. (num_fluids > 1)) then
                c_L = (1d0/gamma_L + 1d0)* &
                      (pres_L + pi_inf_L)/rho_L
                c_R = (1d0/gamma_R + 1d0)* &
                      (pres_R + pi_inf_R)/rho_R
            else
                c_L = &
                    (1d0/gamma_L + 1d0)* &
                    (pres_L + pi_inf_L)/ &
                    (rho_L*(1d0 - alpha_L(num_fluids)))
                c_R = &
                    (1d0/gamma_R + 1d0)* &
                    (pres_R + pi_inf_R)/ &
                    (rho_R*(1d0 - alpha_R(num_fluids)))
            end if
        else

            c_L = ((H_L - 5d-1*sum(vel_L**2d0))/gamma_L)
            c_R = ((H_R - 5d-1*sum(vel_R**2d0))/gamma_R)
        end if

        if (mixture_err .and. c_L < 0d0) then
            c_L = 100.d0*sgm_eps
        else
            c_L = sqrt(c_L)
        end if
        if (mixture_err .and. c_R < 0d0) then
            c_R = 100.d0*sgm_eps
        else
            c_R = sqrt(c_R)
        end if

    end subroutine s_compute_mixture_sound_speeds ! ------------------------

    !>  The purpose of this subroutine is to compute the solution
        !!      for pressure and velocity in the star region when using
        !!      the exact Riemann solver
    subroutine s_compute_star_region() ! -----------------------------------

        real(kind(0d0)) :: change, f_L, dfdp_L, f_R, dfdp_R, pres_old, &
                           pres_start, pres_tol, vel_diff

        integer :: iter

        pres_tol = 1d-10
        change = -1d0*dflt_real

        ! Compute starting guess value for pressure
        call s_guess_pressure(pres_start)

        pres_old = pres_start
        vel_diff = vel_R(dir_idx(1)) - vel_L(dir_idx(1))

        iter = 0

        ! Solve iteratively for pressure in the star region
        do while (change > pres_tol)
            call s_pressure_function(f_L, dfdp_L, pres_old, 1)
            call s_pressure_function(f_R, dfdp_R, pres_old, 2)
            pres_S = pres_old - (f_L + f_R + vel_diff)/(dfdp_L + dfdp_R)
            if (iter > 10000) then
                print '(A)', 'Too many iterations in pressure'
                call s_mpi_abort()
            end if
            change = 2d0*abs((pres_S - pres_old)/(pres_S + pres_old))
            pres_old = pres_S
            iter = iter + 1
        end do

        ! Compute velocity in star region
        vel_S = 5d-1*(vel_L(dir_idx(1)) + vel_R(dir_idx(1)) + f_R - f_L)

    end subroutine s_compute_star_region ! ----------------------------------

    !>  The purpose of this subroutine is to evaluate the pressure
        !!      functions f_K in the exact Riemann solver
        !!  @param f_K    Pressure function
        !!  @param dfdp_K Mixture pressure derivative
        !!  @param pres   Pressure
        !!  @param side   Wave side
    subroutine s_pressure_function(f_K, dfdp_K, pres, side) ! -

        real(kind(0d0)), intent(IN) :: pres
        integer, intent(IN) :: side
        real(kind(0d0)), intent(OUT) :: f_K, dfdp_K
        real(kind(0d0)) :: gam_L, pinf_L, gam_R, pinf_R, c_SL, c_SR

        if (side == 1) then
            gam_L = (gamma_L + 1d0)/gamma_L
            pinf_L = pi_inf_L/(gamma_L + 1d0)
            if (pres + dpres_L <= pres_L) then
                ! Rarefaction wave
                c_SL = c_L*((pres + dpres_L + pinf_L)/(pres_L + pinf_L))**G1_L
                f_K = G4_L*(c_SL/c_L - 1d0)*c_L
                dfdp_K = c_SL/(gam_L*(dpres_L + pres + pinf_L))
            elseif (pres + dpres_L > pres_L) then
                ! Shock wave
                f_K = (c_L/gam_L*((pres + dpres_L)/pres_L - 1d0)*(pres_L/(pres_L + pinf_L)))/ &
                      sqrt(G2_L*((pres + dpres_L)/pres_L - 1d0)*(pres_L/(pres_L + pinf_L)) + 1d0)
                dfdp_K = 2d0*c_L/gam_L/(pres_L + pinf_L)/ &
                         sqrt(2d0*(gam_L + 1d0)*(pres + dpres_L - pres_L)/gam_L/(pres_L + pinf_L) + 4d0) - &
                         2d0*c_L*(pres + dpres_L - pres_L)*(gam_L + 1d0)/gam_L**2d0/(pres_L + pinf_L)**2d0/ &
                         sqrt(2d0*(gam_L + 1d0)*(pres + dpres_L - pres_L)/gam_L/(pres_L + pinf_L) + 4d0)**3d0
            else
                print '(A)', 'Error in evaluating left pressure function. Exiting...'
                call s_mpi_abort()
            end if
        elseif (side == 2) then
            gam_R = (gamma_R + 1d0)/gamma_R
            pinf_R = pi_inf_R/(gamma_R + 1d0)
            if (pres + dpres_R <= pres_R) then
                ! Rarefaction wave
                c_SR = c_R*((pres + dpres_R + pinf_R)/(pres_R + pinf_R))**G1_R
                f_K = G4_R*(c_SR/c_R - 1d0)*c_R
                dfdp_K = c_SR/(gam_R*(dpres_R + pres + pinf_R))
            elseif (pres + dpres_R > pres_R) then
                ! Shock wave
                f_K = (c_R/gam_R*((pres + dpres_R)/pres_R - 1d0)*(pres_R/(pres_R + pinf_R)))/ &
                      sqrt(G2_R*((pres + dpres_R)/pres_R - 1d0)*(pres_R/(pres_R + pinf_R)) + 1d0)
                dfdp_K = 2d0*c_R/gam_R/(pres_R + pinf_R)/ &
                         sqrt(2d0*(gam_R + 1d0)*(pres + dpres_R - pres_R)/gam_R/(pres_R + pinf_R) + 4d0) - &
                         2d0*c_R*(pres + dpres_R - pres_R)*(gam_R + 1d0)/gam_R**2d0/(pres_R + pinf_R)**2d0/ &
                         sqrt(2d0*(gam_R + 1d0)*(pres + dpres_R - pres_R)/gam_R/(pres_R + pinf_R) + 4d0)**3d0
            else
                print '(A)', 'Error in evaluating right pressure function. Exiting...'
                call s_mpi_abort()
            end if
        end if

    end subroutine s_pressure_function ! ------------------------------------

    !>  The purpose of this subroutine is to provide a guess value
        !!      for pressure in the star region. The choice is made
        !!      according to adaptive Riemann solver using the PVRS, TRRS.
        !!      and TSRS approximate Riemann solvers.
        !!  @param pres_start Initial and output pressure
    subroutine s_guess_pressure(pres_start) ! -------------------------------

        real(kind(0d0)), intent(INOUT) :: pres_start

        real(kind(0d0)) :: CUP, pres_max, pres_min, pres_PV, &
                           Q_max, Q_user, pres_TS, pres_TR
        real(kind(0d0)) :: pinf_L, pinf_R, A, B

        Q_user = 2d0

        ! Compute guess pressure from PVRS Riemann solver
        CUP = 25d-2*(rho_L + rho_R)*(c_L + c_R)
        pres_PV = 5d-1*(pres_L + pres_R) + 5d-1*(vel_L(dir_idx(1)) - vel_R(dir_idx(1)))*CUP
        pres_PV = max(0d0, pres_PV)
        pres_min = min(pres_L, pres_R)
        pres_max = max(pres_L, pres_R)
        Q_max = pres_max/pres_min

        if ((Q_max <= Q_user) .and. (pres_min <= pres_PV) .and. (pres_PV <= pres_max)) then

            ! Select PVRS Riemann solver
            pres_start = pres_PV

            if (pres_start /= pres_start) then
                print '(A)', 'NaN guess for pressure using PVRS'
                call s_mpi_abort()
            end if

        elseif (pres_PV < pres_min) then

            ! Select TRRS Riemann solver
            pinf_L = pi_inf_L/(gamma_L + 1d0)
            pinf_R = pi_inf_R/(gamma_R + 1d0)
            pres_TR = ((((-G4_R*c_R*(((pres_PV + dpres_R + pinf_R)/(pres_R + pinf_R))**G1_R - 1d0) - &
                          (vel_R(dir_idx(1)) - vel_L(dir_idx(1))))/c_L/G4_L) + 1d0)**G3_L)* &
                      (pres_L + pinf_L) - (dpres_L + pinf_L)
            pres_start = max(0d0, pres_TR)

            if (pres_start /= pres_start) then
                print '(A)', 'NaN guess for pressure using TRRS'
                call s_mpi_abort()
            end if

        else
            ! Select TSRS Riemann solver with pres_PV as estimate
            A = sqrt(G2_L*((pres_PV + dpres_L - pres_L)/(pres_L + pinf_L)) + 1d0)
            B = sqrt(G2_R*((pres_PV + dpres_R - pres_R)/(pres_R + pinf_R)) + 1d0)
            pres_TS = (c_L*gamma_L/(gamma_L + 1d0)/A*((pres_L - dpres_L)/(pres_L + pinf_L)) + &
                       c_R*gamma_R/(gamma_R + 1d0)/B*((pres_R - dpres_R)/(pres_R + pinf_R)) - &
                       (vel_R(dir_idx(1)) - vel_L(dir_idx(1))))/ &
                      (c_L*gamma_L/(gamma_L + 1d0)/A/(pres_L + pinf_L) + &
                       c_R*gamma_R/(gamma_R + 1d0)/B/(pres_R + pinf_R))
            pres_start = max(0d0, pres_TS)

            if (pres_start /= pres_start) then
                print '(A)', 'NaN guess for pressure using TSRS'
                call s_mpi_abort()
            end if

        end if

    end subroutine s_guess_pressure ! --------------------------------------

    !> Computes the averaged intrercell variables for the Riemann solver
    subroutine s_compute_intercell_solution() ! --------------------

        integer :: i

        real(kind(0d0)) :: c_IC
        real(kind(0d0)) :: s_HL, S_TL, c_SL, pres_SL, s_L
        real(kind(0d0)) :: s_HR, S_TR, c_SR, pres_SR, s_R

        if (0d0 <= vel_S) then
            ! IC lies to the left of the contact discontinuity
            if (pres_S + dpres_L <= pres_L) then
                ! Left rarefaction
                s_HL = vel_L(dir_idx(1)) - c_L

                if (0d0 <= s_HL) then
                    ! IC is left data state
                    do i = 1, cont_idx%end
                        alpha_rho_IC(i) = alpha_rho_L(i)
                    end do
                    rho_IC = rho_L

                    do i = 1, num_dims
                        vel_IC(i) = vel_L(i)
                    end do

                    pres_IC = pres_L
                    E_IC = E_L

                    do i = 1, num_fluids
                        alpha_IC(i) = alpha_L(i)
                    end do

                else
                    c_SL = c_L*((pres_S + dpres_L + pi_inf_L/(gamma_L + 1d0))/(pres_L + pi_inf_L/(gamma_L + 1d0)))**G1_L
                    S_TL = vel_S - c_SL

                    if (0d0 > S_TL) then
                        ! IC is star left state
                        do i = 1, cont_idx%end
                            alpha_rho_IC(i) = alpha_rho_L(i)*((pres_S + dpres_L + pi_inf_L/(gamma_L + 1d0))/ &
                                                              (pres_L + pi_inf_L/(gamma_L + 1d0)))**(gamma_L/(gamma_L + 1d0))
                        end do
                        rho_IC = rho_L*((pres_S + dpres_L + pi_inf_L/(gamma_L + 1d0))/ &
                                        (pres_L + pi_inf_L/(gamma_L + 1d0)))**(gamma_L/(gamma_L + 1d0))

                        vel_IC(dir_idx(1)) = vel_S
                        do i = 2, num_dims
                            vel_IC(dir_idx(i)) = vel_L(dir_idx(i))
                        end do

                        pres_IC = pres_S + dpres_L
                        E_IC = gamma_L*pres_IC + pi_inf_L + 5d-1*rho_IC*sum(vel_IC**2d0)

                        do i = 1, num_fluids
                            alpha_IC(i) = alpha_L(i)
                        end do

                    else
                        ! IC is inside left rarefaction
                        vel_IC(dir_idx(1)) = G5_L*(c_L + G7_L*vel_L(dir_idx(1)) + 0d0)
                        c_IC = G5_L*(c_L + G7_L*(vel_L(dir_idx(1)) - 0d0))

                        do i = 1, cont_idx%end
                            alpha_rho_IC(i) = alpha_rho_L(i)*(c_IC/c_L)**G4_L
                        end do
                        rho_IC = rho_L*(c_IC/c_L)**G4_L

                        do i = 2, num_dims
                            vel_IC(dir_idx(i)) = vel_L(dir_idx(i))
                        end do

                        pres_IC = (pres_L + pi_inf_L/(gamma_L + 1d0))*(c_IC/c_L)**G3_L - (pi_inf_L/(gamma_L + 1d0))
                        E_IC = gamma_L*pres_IC + pi_inf_L + 5d-1*rho_IC*sum(vel_IC**2d0)

                        do i = 1, num_fluids
                            alpha_IC(i) = alpha_L(i)
                        end do

                    end if
                end if
            else
                ! Left shock
                pres_SL = (pres_S + dpres_L + pi_inf_L/(gamma_L + 1d0))/(pres_L + pi_inf_L/(gamma_L + 1d0))
                s_L = vel_L(dir_idx(1)) - c_L*sqrt(G2_L*(pres_S + dpres_L - pres_L)/(pres_L + pi_inf_L/(gamma_L + 1d0)) + 1d0)

                if (0d0 <= s_L) then
                    ! IC is left data state
                    do i = 1, cont_idx%end
                        alpha_rho_IC(i) = alpha_rho_L(i)
                    end do
                    rho_IC = rho_L

                    do i = 1, num_dims
                        vel_IC(i) = vel_L(i)
                    end do

                    pres_IC = pres_L
                    E_IC = E_L

                    do i = 1, num_fluids
                        alpha_IC(i) = alpha_L(i)
                    end do

                else
                    ! IC is star left state
                    do i = 1, cont_idx%end
                        alpha_rho_IC(i) = alpha_rho_L(i)*(pres_SL + G6_L)/(pres_SL*G6_L + 1d0)
                    end do
                    rho_IC = rho_L*(pres_SL + G6_L)/(pres_SL*G6_L + 1d0)

                    vel_IC(dir_idx(1)) = vel_S
                    do i = 2, num_dims
                        vel_IC(dir_idx(i)) = vel_L(dir_idx(i))
                    end do

                    pres_IC = pres_S + dpres_L
                    E_IC = gamma_L*pres_IC + pi_inf_L + 5d-1*rho_IC*sum(vel_IC**2d0)

                    do i = 1, num_fluids
                        alpha_IC(i) = alpha_L(i)
                    end do

                end if
            end if
        else
            ! IC is to the right of the contact discontinuity
            if (pres_S + dpres_R > pres_R) then
                ! Right shock
                pres_SR = (pres_S + dpres_R + pi_inf_R/(gamma_R + 1d0))/(pres_R + pi_inf_R/(gamma_R + 1d0))
                s_R = vel_R(dir_idx(1)) + c_R*sqrt(G2_R*(pres_S + dpres_R - pres_R)/(pres_R + pi_inf_R/(gamma_R + 1d0)) + 1d0)

                if (0d0 >= s_R) then
                    ! IC is right data state
                    do i = 1, cont_idx%end
                        alpha_rho_IC(i) = alpha_rho_R(i)
                    end do
                    rho_IC = rho_R

                    do i = 1, num_dims
                        vel_IC(i) = vel_R(i)
                    end do

                    pres_IC = pres_R
                    E_IC = E_R

                    do i = 1, num_fluids
                        alpha_IC(i) = alpha_R(i)
                    end do

                else
                    ! IC is star right state
                    do i = 1, cont_idx%end
                        alpha_rho_IC(i) = alpha_rho_R(i)*(pres_SR + G6_R)/(pres_SR*G6_R + 1d0)
                    end do
                    rho_IC = rho_R*(pres_SR + G6_R)/(pres_SR*G6_R + 1d0)

                    vel_IC(dir_idx(1)) = vel_S
                    do i = 2, num_dims
                        vel_IC(dir_idx(i)) = vel_R(dir_idx(i))
                    end do

                    pres_IC = pres_S + dpres_R
                    E_IC = gamma_R*pres_IC + pi_inf_R + 5d-1*rho_IC*sum(vel_IC**2d0)

                    do i = 1, num_fluids
                        alpha_IC(i) = alpha_R(i)
                    end do

                end if
            else
                ! Right rarefaction
                s_HR = vel_R(dir_idx(1)) + c_R

                if (0d0 >= s_HR) then
                    ! IC is right data state
                    do i = 1, cont_idx%end
                        alpha_rho_IC(i) = alpha_rho_R(i)
                    end do
                    rho_IC = rho_R

                    do i = 1, num_dims
                        vel_IC(i) = vel_R(i)
                    end do

                    pres_IC = pres_R
                    E_IC = E_R

                    do i = 1, num_fluids
                        alpha_IC(i) = alpha_R(i)
                    end do

                else
                    c_SR = c_R*((pres_S + dpres_R + pi_inf_R/(gamma_R + 1d0))/(pres_R + pi_inf_R/(gamma_R + 1d0)))**G1_R
                    S_TR = vel_S + c_SR

                    if (0d0 <= S_TR) then
                        ! IC is star right state
                        do i = 1, cont_idx%end
                            alpha_rho_IC(i) = alpha_rho_R(i)*((pres_S + dpres_R + pi_inf_R/(gamma_R + 1d0))/ &
                                                              (pres_R + pi_inf_R/(gamma_R + 1d0)))**(gamma_R/(gamma_R + 1d0))
                        end do
                        rho_IC = rho_R*((pres_S + dpres_R + pi_inf_R/(gamma_R + 1d0))/(pres_R + &
                                                                                       pi_inf_R/(gamma_R + 1d0)))**(gamma_R/(gamma_R + 1d0))

                        vel_IC(dir_idx(1)) = vel_S
                        do i = 2, num_dims
                            vel_IC(dir_idx(i)) = vel_R(dir_idx(i))
                        end do

                        pres_IC = pres_S + dpres_R
                        E_IC = gamma_R*pres_IC + pi_inf_R + 5d-1*rho_IC*sum(vel_IC**2d0)

                        do i = 1, num_fluids
                            alpha_IC(i) = alpha_R(i)
                        end do

                    else
                        ! IC is inside right rarefaction
                        vel_IC(dir_idx(1)) = G5_R*(-1d0*c_R + G7_R*vel_R(dir_idx(1)) + 0d0)
                        c_IC = G5_R*(c_R - G7_R*(vel_R(dir_idx(1)) - 0d0))

                        do i = 1, cont_idx%end
                            alpha_rho_IC(i) = alpha_rho_R(i)*(c_IC/c_R)**G4_R
                        end do
                        rho_IC = rho_R*(c_IC/c_R)**G4_R

                        do i = 2, num_dims
                            vel_IC(dir_idx(i)) = vel_R(dir_idx(i))
                        end do

                        pres_IC = (pres_R + pi_inf_R/(gamma_R + 1d0))*(c_IC/c_R)**G3_R - (pi_inf_R/(gamma_R + 1d0))
                        E_IC = gamma_R*pres_IC + pi_inf_R + 5d-1*rho_IC*sum(vel_IC**2d0)

                        do i = 1, num_fluids
                            alpha_IC(i) = alpha_R(i)
                        end do

                    end if
                end if
            end if
        end if

    end subroutine s_compute_intercell_solution ! --------------------------

    !>  The procedure computes the Roe average density, velocity,
        !!      enthalpy, volume fractions, specific heat ratio function,
        !!      speed of sound, shear and volume Reynolds numbers, Weber
        !!      numbers and curvatures, at the cell-boundaries, from the
        !!      left and right states of the Riemann problem.
        !! @param j  First coordinate index
        !! @param k Second coordinate index
        !! @param l  Third coordinate index
    subroutine s_compute_roe_average_state(qL_prim_rs_vf, qR_prim_rs_vf,j, k, l) ! ---------------
!$acc routine seq
        type(scalar_field), dimension(sys_size), intent(IN) :: qL_prim_rs_vf, qR_prim_rs_vf
        integer, intent(IN) :: j, k, l

        integer :: i

        ! Left and Right Riemann Problem States ============================


        call s_convert_species_to_mixture_variables_acc( &
                                            rho_L, gamma_L, &
                                            pi_inf_L, alpha_L, alpha_rho_L, &
                                            j, k, l)
        call s_convert_species_to_mixture_variables_acc( &
                                            rho_R, gamma_R, &
                                            pi_inf_R, alpha_R, alpha_rho_R, &
                                            j + 1, k, l)

        E_L = gamma_L*pres_L + pi_inf_L + 5d-1*rho_L*sum(vel_L**2d0)
        E_R = gamma_R*pres_R + pi_inf_R + 5d-1*rho_R*sum(vel_R**2d0)

        H_L = (E_L + pres_L)/rho_L
        H_R = (E_R + pres_R)/rho_R

        call s_compute_mixture_sound_speeds(qL_prim_rs_vf, qR_prim_rs_vf,j, k, l)

        ! ==================================================================

        ! Roe Average Riemann Problem State ================================
        rho_avg = sqrt(rho_L*rho_R)

        vel_avg = (sqrt(rho_L)*vel_L + sqrt(rho_R)*vel_R)/ &
                  (sqrt(rho_L) + sqrt(rho_R))

        H_avg = (sqrt(rho_L)*H_L + sqrt(rho_R)*H_R)/ &
                (sqrt(rho_L) + sqrt(rho_R))

        gamma_avg = (sqrt(rho_L)*gamma_L + sqrt(rho_R)*gamma_R)/ &
                    (sqrt(rho_L) + sqrt(rho_R))

        if (mixture_err) then
            if ((H_avg - 5d-1*sum(vel_avg**2d0)) < 0d0) then
                c_avg = 1d-16
            else
                c_avg = sqrt((H_avg - 5d-1*sum(vel_avg**2d0))/gamma_avg)
            end if
        else
            c_avg = sqrt((H_avg - 5d-1*sum(vel_avg**2d0))/gamma_avg)
        end if


        do i = 1, 2
            if (Re_size(i) > 0) then
                Re_avg_rs_vf(i)%sf(j, k, l) = 2d0/(1d0/Re_L(i) + 1d0/Re_R(i))
            end if
        end do

        ! ==================================================================

    end subroutine s_compute_roe_average_state ! ---------------------------

    !>  This procedure calculates the arithmetic average density,
        !!      velocity, enthalpy, volume fractions, specIFic heat ratio
        !!      function, sound speed, shear and volume Reynolds numbers,
        !!      Weber numbers and the curvatures, at the cell-boundaries,
        !!      from the left and right states of the Riemann problem.
        !!  @param j  First coordinate index
        !!  @param k Second coordinate index
        !!  @param l  Third coordinate index
    subroutine s_compute_arithmetic_average_state(qL_prim_rs_vf, qR_prim_rs_vf, j, k, l) ! --------
!$acc routine seq
        type(scalar_field), dimension(sys_size), intent(IN) :: qL_prim_rs_vf, qR_prim_rs_vf

        integer, intent(IN) :: j, k, l

        integer :: i, q !< Generic loop iterator

        !ensemble-averaged bubble variables
        real(kind(0d0)) :: PbwR3Lbar, Pbwr3Rbar
        real(kind(0d0)) :: R3Lbar, R3Rbar
        real(kind(0d0)) :: R3V2Lbar, R3V2Rbar

        ! Left and Right Riemann Problem States ============================



        !call s_convert_species_to_mixture_variables_acc( &
        !                                    rho_L, gamma_L, &
        !                                    pi_inf_L, alpha_L, alpha_rho_L, &
        !                                    j, k, l)
        !call s_convert_species_to_mixture_variables_acc( &
        !                                    rho_R, gamma_R, &
        !                                    pi_inf_R, alpha_R, alpha_rho_R, &
        !                                    j + 1, k, l)



        ! Compute left/right states for bubble number density
        if (bubbles) then
            do i = 1, num_fluids
                alpha_L(i) = qL_prim_rs_vf(E_idx + i)%sf(j, k, l)
                alpha_R(i) = qR_prim_rs_vf(E_idx + i)%sf(j + 1, k, l)
            end do!

            do i = 1, nb
                R0_L(i) = qL_prim_rs_vf(bub_idx%rs(i))%sf(j, k, l)
                R0_R(i) = qR_prim_rs_vf(bub_idx%rs(i))%sf(j + 1, k, l)

                V0_L(i) = qL_prim_rs_vf(bub_idx%vs(i))%sf(j, k, l)
                V0_R(i) = qR_prim_rs_vf(bub_idx%vs(i))%sf(j + 1, k, l)
                if (.not. polytropic) then
                    P0_L(i) = qL_prim_rs_vf(bub_idx%ps(i))%sf(j, k, l)
                    P0_R(i) = qR_prim_rs_vf(bub_idx%ps(i))%sf(j + 1, k, l)
                end if
            end do

            call s_comp_n_from_prim(alpha_L(num_fluids), R0_L, nbub_L)
            call s_comp_n_from_prim(alpha_R(num_fluids), R0_R, nbub_R)

            do i = 1, nb
                if (.not. qbmm) then
                    if (polytropic) then
                        pbw_L(i) = f_cpbw_KM(R0(i), R0_L(i), V0_L(i), 0d0)
                        pbw_R(i) = f_cpbw_KM(R0(i), R0_R(i), V0_R(i), 0d0)
                    else
                        pbw_L(i) = f_cpbw_KM(R0(i), R0_L(i), V0_L(i), P0_L(i))
                        pbw_R(i) = f_cpbw_KM(R0(i), R0_R(i), V0_R(i), P0_R(i))
                    end if
                end if
            end do

            if (qbmm) then
                PbwR3Lbar = mom_sp(4)%sf(j, k, l)
                PbwR3Rbar = mom_sp(4)%sf(j + 1, k, l)

                R3Lbar = mom_sp(1)%sf(j, k, l)
                R3Rbar = mom_sp(1)%sf(j + 1, k, l)

                R3V2Lbar = mom_sp(3)%sf(j, k, l)
                R3V2Rbar = mom_sp(3)%sf(j + 1, k, l)
            else
                call s_quad(pbw_L*(R0_L**3.d0), PbwR3Lbar)
                call s_quad(pbw_R*(R0_R**3.d0), PbwR3Rbar)

                call s_quad(R0_L**3.d0, R3Lbar)
                call s_quad(R0_R**3.d0, R3Rbar)

                call s_quad((R0_L**3.d0)*(V0_L**2.d0), R3V2Lbar)
                call s_quad((R0_R**3.d0)*(V0_R**2.d0), R3V2Rbar)
            end if

            !ptilde = \alf( pl - \bar{ pbw R^3)/\bar{R^3} - rho \bar{R^3 \Rdot^2}/\bar{R^3} )
            if (alpha_L(num_fluids) < small_alf .or. R3Lbar < small_alf) then
                ptilde_L = alpha_L(num_fluids)*pres_L
            else
                ptilde_L = alpha_L(num_fluids)*(pres_L - PbwR3Lbar/R3Lbar - &
                                                rho_L*R3V2Lbar/R3Lbar)
            end if

            if (alpha_R(num_fluids) < small_alf .or. R3Rbar < small_alf) then
                ptilde_R = alpha_R(num_fluids)*pres_R
            else
                ptilde_R = alpha_R(num_fluids)*(pres_R - PbwR3Rbar/R3Rbar - &
                                                rho_R*R3V2Rbar/R3Rbar)
            end if

            if ((ptilde_L .ne. ptilde_L) .or. (ptilde_R .ne. ptilde_R)) then
                print *, 'Ptilde NaN at ', j, k, l, x_cb(j)
                print *, nbub_L, alpha_L, pres_L, PbwR3Lbar, R3Lbar, rho_L, R3V2Lbar, R3Lbar
                print *, nbub_R, alpha_R, pres_R, PbwR3Rbar, R3Rbar, rho_R, R3V2Rbar, R3Rbar
                call s_mpi_abort()
            end if

            ptil(j, k, l) = 0.5d0*(ptilde_L + ptilde_R)
        end if

        call s_compute_mixture_sound_speeds(qL_prim_rs_vf, qR_prim_rs_vf,j, k, l)

        ! ==================================================================

        ! Arithmetic Average Riemann Problem State =========================


        do i = 1, 2
            if (Re_size(i) > 0) then
                Re_avg_rs_vf(i)%sf(j, k, l) = 2d0/(1d0/Re_L(i) + 1d0/Re_R(i))
            end if
        end do

        !
    end subroutine s_compute_arithmetic_average_state ! --------------------

    !>  The direct estimation of the left, right and middle wave
        !!      speeds, proposed by Batten et al. (1997) that results in
        !!      the exact resolution of isolated shock and contact waves.
        !!  @param j  First coordinate index
        !!  @param k Second coordinate index
        !!  @param l  Third coordinate index
    subroutine s_compute_direct_wave_speeds(j, k, l) ! -----------------------
!$acc routine seq
        integer, intent(IN) :: j, k, l

        real(kind(0d0)) :: denom

        integer :: i !< Generic loop iterator


        s_L = min(vel_L(dir_idx(1)) - c_L, vel_R(dir_idx(1)) - c_R)
        s_R = max(vel_R(dir_idx(1)) + c_R, vel_L(dir_idx(1)) + c_L)

        s_S = (pres_R - pres_L + rho_L*vel_L(dir_idx(1))* &
               (s_L - vel_L(dir_idx(1))) - &
               rho_R*vel_R(dir_idx(1))* &
               (s_R - vel_R(dir_idx(1)))) &
              /(rho_L*(s_L - vel_L(dir_idx(1))) - &
                rho_R*(s_R - vel_R(dir_idx(1))))
        denom = rho_L*(s_L - vel_L(dir_idx(1))) - rho_R*(s_R - vel_R(dir_idx(1)))


    end subroutine s_compute_direct_wave_speeds ! --------------------------

    !>  Estimation of the left, right and star region wave speeds
        !!      by the approximation of the pressures and velocity in the
        !!      star regions, see Toro (1999). The pressures and velocity
        !!      are approximated by using the primitive variables Riemann
        !!      solver (PVRS) and the wave speeds are then estimated from
        !!      those approximations using the exact wave relations.
        !!  @param j  First coordinate index
        !!  @param k Second coordinate index
        !!  @param l  Third coordinate index
    subroutine s_compute_pressure_velocity_wave_speeds(j, k, l) ! ------------
!$acc routine seq
        integer, intent(IN) :: j, k, l

        ! Left and right pressures in the star region
        real(kind(0d0)) :: pres_SL, pres_SR


        ! Left and right shock Mach numbers
        real(kind(0d0)) :: Ms_L, Ms_R

        integer :: i !< Generic loop iterator


        pres_SL = 5d-1*(pres_L + pres_R + rho_avg*c_avg* &
                        (vel_L(dir_idx(1)) - &
                         vel_R(dir_idx(1))))
        pres_SR = pres_SL

        Ms_L = max(1d0, sqrt(1d0 + ((5d-1 + gamma_L)/(1d0 + gamma_L))* &
                             (pres_SL/pres_L - 1d0)*pres_L/ &
                             ((pres_L + pi_inf_L/(1d0 + gamma_L)))))
        Ms_R = max(1d0, sqrt(1d0 + ((5d-1 + gamma_R)/(1d0 + gamma_R))* &
                             (pres_SR/pres_R - 1d0)*pres_R/ &
                             ((pres_R + pi_inf_R/(1d0 + gamma_R)))))

        s_L = vel_L(dir_idx(1)) - c_L*Ms_L
        s_R = vel_R(dir_idx(1)) + c_R*Ms_R

        s_S = 5d-1*((vel_L(dir_idx(1)) + vel_R(dir_idx(1))) + &
                    (pres_L - pres_R)/ &
                    (rho_avg*c_avg))

    end subroutine s_compute_pressure_velocity_wave_speeds ! ---------------

    !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_riemann_solvers_module() ! ---------------------

        ! Allocating the variables that will be utilized to formulate the
        ! left, right, and average states of the Riemann problem, as well
        ! the Riemann problem solution
        integer :: i


        allocate(gammas(1:num_fluids))
        allocate(pi_infs(1:num_fluids))


        do i = 1, num_fluids
            gammas(i) = fluid_pp(i)%gamma
            pi_infs(i) = fluid_pp(i)%pi_inf
        end do
!$acc update device(gammas, pi_infs)

        momxb = mom_idx%beg; momxe = mom_idx%end
        contxb = cont_idx%beg; contxe = cont_idx%end
        bubxb = bub_idx%beg; bubxe = bub_idx%end
        advxb = adv_idx%beg; advxe = adv_idx%end
        intxb = internalEnergies_idx%beg; intxe = internalEnergies_idx%end
!$acc update device(momxb, momxe, contxb, contxe, bubxb, bubxe, advxb, advxe, intxb, intxe)

        if(bubbles) then
            allocate(rs(1:nb))
            allocate(vs(1:nb))
            if(.not. polytropic) then
                allocate(ps(1:nb))
                allocate(ms(1:nb))
            end if
            
            rs = bub_idx%rs
            vs = bub_idx%vs
            if(.not. polytropic) then
                ps = bub_idx%ps
                ms = bub_idx%ms
            end if

            print *, "Rad Index", rs(1)

!$acc update device(rs, vs)
            if(.not. polytropic) then
!$acc update device(ps, ms)
            end if 
            
        end if


        allocate (qL_prim_rsx_vf(1:sys_size), qR_prim_rsx_vf(1:sys_size))
        allocate (qL_prim_rsy_vf(1:sys_size), qR_prim_rsy_vf(1:sys_size))
        allocate (qL_prim_rsz_vf(1:sys_size), qR_prim_rsz_vf(1:sys_size))

        allocate (flux_rsx_vf(1:sys_size), flux_src_rsx_vf(1:sys_size))
        allocate (flux_rsy_vf(1:sys_size), flux_src_rsy_vf(1:sys_size))
        allocate (flux_rsz_vf(1:sys_size), flux_src_rsz_vf(1:sys_size))

        allocate (flux_gsrc_rsx_vf(1:sys_size))
        allocate (flux_gsrc_rsy_vf(1:sys_size))
        allocate (flux_gsrc_rsz_vf(1:sys_size))

        allocate (vel_src_rsx_vf(1:num_dims))
        allocate (vel_src_rsy_vf(1:num_dims))
        allocate (vel_src_rsz_vf(1:num_dims))

        if (any(Re_size > 0)) allocate (Re_avg_rsx_vf(1:2))
        if (any(Re_size > 0)) allocate (Re_avg_rsy_vf(1:2))
        if (any(Re_size > 0)) allocate (Re_avg_rsz_vf(1:2))

        allocate (alpha_rho_L(1:cont_idx%end), vel_L(1:num_dims))
        allocate (alpha_rho_R(1:cont_idx%end), vel_R(1:num_dims))

        allocate (vel_avg(1:num_dims))

        allocate (alpha_L(1:num_fluids))
        allocate (alpha_R(1:num_fluids))




        if (riemann_solver == 3) then
            allocate (alpha_rho_IC(1:cont_idx%end), vel_IC(1:num_dims))
            allocate (alpha_IC(1:num_fluids))
        end if

        ! Associating procedural pointer to the subroutine that will be
        ! utilized to calculate the solution of a given Riemann problem
        if (riemann_solver == 1) then
            s_riemann_solver => s_hll_riemann_solver_acc
        elseif (riemann_solver == 2) then
            s_riemann_solver => s_hllc_riemann_solver_acc
        else
            s_riemann_solver => s_exact_riemann_solver
        end if



        if (bubbles) then
            allocate (R0_L(nb), R0_R(nb))
            allocate (V0_L(nb), V0_R(nb))
            allocate (pbw_L(nb), pbw_R(nb))
            if (qbmm) then
                allocate (moms_L(nb, nmom), moms_R(nb, nmom))
            else
                if (.not. polytropic) then
                    allocate (P0_L(nb), P0_R(nb))
                end if
            end if
        end if

        ! Associating the procedural pointers to the procedures that will be
        ! utilized to compute the average state and estimate the wave speeds
        if (riemann_solver /= 3) then

            if (avg_state == 1) then
                s_compute_average_state => s_compute_roe_average_state
            else
                s_compute_average_state => s_compute_arithmetic_average_state
            end if

            if (wave_speeds == 1) then
                s_compute_wave_speeds => s_compute_direct_wave_speeds
            else
                s_compute_wave_speeds => s_compute_pressure_velocity_wave_speeds
            end if

        end if

        ! Associating procedural pointer to the subroutine that will be
        ! utilized to compute the viscous source flux
        if (grid_geometry == 3) then
            s_compute_viscous_source_flux => s_compute_cylindrical_viscous_source_flux
        else
            s_compute_viscous_source_flux => s_compute_cartesian_viscous_source_flux
        end if

        ! Associating the procedural pointer to the appropriate subroutine
        ! that will be utilized in the conversion to the mixture variables

        if (model_eqns == 1) then        ! Gamma/pi_inf model
            s_convert_to_mixture_variables => &
                s_convert_mixture_to_mixture_variables
        elseif (bubbles) then           ! Volume fraction for bubbles
            s_convert_to_mixture_variables => &
                s_convert_species_to_mixture_variables_bubbles
        else                            ! Volume fraction model
            s_convert_to_mixture_variables => &
                s_convert_species_to_mixture_variables
        end if

        is1%beg = -1; is2%beg = 0; is3%beg = 0
        is1%end = m; is2%end = n; is3%end = p


        !allocate(qL_prim_rsx_vf_flat(is1%beg:is1%end, is2%beg:is2%end, is3%beg:is3%end, 1:sys_size))
        !allocate(qR_prim_rsx_vf_flat(is1%beg + 1:is1%end + 1, is2%beg:is2%end, is3%beg:is3%end, 1:sys_size))
        allocate (flux_rsx_vf_flat(is1%beg:is1%end, &
                                       is2%beg:is2%end, &
                                       is3%beg:is3%end, 1:sys_size))
        allocate (flux_gsrc_rsx_vf_flat(is1%beg:is1%end, &
                                            is2%beg:is2%end, &
                                            is3%beg:is3%end, 1:sys_size))
        allocate (flux_src_rsx_vf_flat(is1%beg:is1%end, &
                                           is2%beg:is2%end, &
                                           is3%beg:is3%end, advxb:sys_size))
        allocate (vel_src_rsx_vf_flat(is1%beg:is1%end, &
                                               is2%beg:is2%end, &
                                               is3%beg:is3%end, 1:num_dims))
        if(Re_size(1) > 0) then
            allocate (Re_avg_rsx_vf_flat(is1%beg:is1%end, &
                                             is2%beg:is2%end, &
                                             is3%beg:is3%end, 1:2))
        end if




        if(n == 0) return

        is1%beg = -1; is2%beg = 0; is3%beg = 0
        is1%end = n; is2%end = m; is3%end = p

        !allocate(qL_prim_rsy_vf_flat(is1%beg:is1%end, is2%beg:is2%end, is3%beg:is3%end, 1:sys_size))
        !allocate(qR_prim_rsy_vf_flat(is1%beg + 1:is1%end + 1, is2%beg:is2%end, is3%beg:is3%end, 1:sys_size))
        allocate (flux_rsy_vf_flat(is1%beg:is1%end, &
                                       is2%beg:is2%end, &
                                       is3%beg:is3%end, 1:sys_size))
        allocate (flux_gsrc_rsy_vf_flat(is1%beg:is1%end, &
                                            is2%beg:is2%end, &
                                            is3%beg:is3%end, 1:sys_size))
        allocate (flux_src_rsy_vf_flat(is1%beg:is1%end, &
                                           is2%beg:is2%end, &
                                           is3%beg:is3%end, advxb:sys_size))
        allocate (vel_src_rsy_vf_flat(is1%beg:is1%end, &
                                               is2%beg:is2%end, &
                                             is3%beg:is3%end, 1:num_dims))
        if(Re_size(1) > 0) then
            allocate (Re_avg_rsy_vf_flat(is1%beg:is1%end, &
                                             is2%beg:is2%end, &
                                             is3%beg:is3%end, 1:2))
        end if


        if(p == 0) return

        is1%beg = -1; is2%beg = 0; is3%beg = 0
        is1%end = p; is2%end = n; is3%end = m


        !allocate(qL_prim_rsz_vf_flat(is1%beg:is1%end, is2%beg:is2%end, is3%beg:is3%end, 1:sys_size))
        !allocate(qR_prim_rsz_vf_flat(is1%beg + 1:is1%end + 1, is2%beg:is2%end, is3%beg:is3%end, 1:sys_size))
        allocate (flux_rsz_vf_flat(is1%beg:is1%end, &
                                       is2%beg:is2%end, &
                                       is3%beg:is3%end, 1:sys_size))
        allocate (flux_gsrc_rsz_vf_flat(is1%beg:is1%end, &
                                            is2%beg:is2%end, &
                                            is3%beg:is3%end, 1:sys_size))
        allocate (flux_src_rsz_vf_flat(is1%beg:is1%end, &
                                           is2%beg:is2%end, &
                                           is3%beg:is3%end, advxb:sys_size))
        allocate (vel_src_rsz_vf_flat(is1%beg:is1%end, &
                                               is2%beg:is2%end, &
                                               is3%beg:is3%end, 1:num_dims))
        if(Re_size(1) > 0) then
            allocate (Re_avg_rsz_vf_flat(is1%beg:is1%end, &
                                             is2%beg:is2%end, &
                                             is3%beg:is3%end, 1:2))
        end if



    end subroutine s_initialize_riemann_solvers_module ! -------------------

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
    subroutine s_populate_riemann_states_variables_buffers( & ! ------------
        qL_prim_rsx_vf_flat, qL_prim_rsy_vf_flat, qL_prim_rsz_vf_flat,  dqL_prim_dx_vf, &
        dqL_prim_dy_vf, &
        dqL_prim_dz_vf, &
        gm_alphaL_vf, &
        qR_prim_rsx_vf_flat, qR_prim_rsy_vf_flat, qR_prim_rsz_vf_flat, dqR_prim_dx_vf, &
        dqR_prim_dy_vf, &
        dqR_prim_dz_vf, &
        gm_alphaR_vf, &
        norm_dir, ix, iy, iz)

        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:), intent(INOUT) :: qL_prim_rsx_vf_flat, qL_prim_rsy_vf_flat, qL_prim_rsz_vf_flat, qR_prim_rsx_vf_flat, qR_prim_rsy_vf_flat, qR_prim_rsz_vf_flat

        type(scalar_field), &
            allocatable, dimension(:), &
            intent(INOUT) :: dqL_prim_dx_vf, dqR_prim_dx_vf, &
                             dqL_prim_dy_vf, dqR_prim_dy_vf, &
                             dqL_prim_dz_vf, dqR_prim_dz_vf, &
                             gm_alphaL_vf, gm_alphaR_vf

        integer, intent(IN) :: norm_dir

        type(bounds_info), intent(IN) :: ix, iy, iz

        integer :: i, j, k, l !< Generic loop iterator

        if (norm_dir == 1) then
            is1 = ix; is2 = iy; is3 = iz
            dir_idx = (/1, 2, 3/); dir_flg = (/1d0, 0d0, 0d0/)
        elseif (norm_dir == 2) then
            is1 = iy; is2 = ix; is3 = iz
            dir_idx = (/2, 1, 3/); dir_flg = (/0d0, 1d0, 0d0/)
        else
            is1 = iz; is2 = iy; is3 = ix
            dir_idx = (/3, 1, 2/); dir_flg = (/0d0, 0d0, 1d0/)
        end if

!$acc update device(is1, is2, is3, dir_idx, dir_flg)


        ! Population of Buffers in x-direction =============================
        if (norm_dir == 1) then

            if (bc_x%beg == -4) then    ! Riemann state extrap. BC at beginning
!$acc parallel loop collapse(3) gang vector default(present)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                    qL_prim_rsx_vf_flat(-1, k, l, i) = &
                        qR_prim_rsx_vf_flat(0, k, l, i)
                        end do
                    end do
                end do

                if (any(Re_size > 0)) then

                    do i = mom_idx%beg, mom_idx%end
                        dqL_prim_dx_vf(i)%sf(-1, &
                                             iy%beg:iy%end, &
                                             iz%beg:iz%end) = &
                            dqR_prim_dx_vf(i)%sf(0, &
                                                 iy%beg:iy%end, &
                                                 iz%beg:iz%end)
                    end do

                    if (n > 0) then

                        do i = mom_idx%beg, mom_idx%end
                            dqL_prim_dy_vf(i)%sf(-1, &
                                                 iy%beg:iy%end, &
                                                 iz%beg:iz%end) = &
                                dqR_prim_dy_vf(i)%sf(0, &
                                                     iy%beg:iy%end, &
                                                     iz%beg:iz%end)
                        end do

                        if (p > 0) then
                            do i = mom_idx%beg, mom_idx%end
                                dqL_prim_dz_vf(i)%sf(-1, &
                                                     iy%beg:iy%end, &
                                                     iz%beg:iz%end) = &
                                    dqR_prim_dz_vf(i)%sf(0, &
                                                         iy%beg:iy%end, &
                                                         iz%beg:iz%end)
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
                    qL_prim_rsx_vf_flat(m + 1, k, l, i) = &
                        qR_prim_rsx_vf_flat(m, k, l, i)
                        end do
                    end do
                end do


                if (any(Re_size > 0)) then

                    do i = mom_idx%beg, mom_idx%end
                        dqR_prim_dx_vf(i)%sf(m + 1, &
                                             iy%beg:iy%end, &
                                             iz%beg:iz%end) = &
                            dqL_prim_dx_vf(i)%sf(m, &
                                                 iy%beg:iy%end, &
                                                 iz%beg:iz%end)
                    end do

                    if (n > 0) then

                        do i = mom_idx%beg, mom_idx%end
                            dqR_prim_dy_vf(i)%sf(m + 1, &
                                                 iy%beg:iy%end, &
                                                 iz%beg:iz%end) = &
                                dqL_prim_dy_vf(i)%sf(m, &
                                                     iy%beg:iy%end, &
                                                     iz%beg:iz%end)
                        end do

                        if (p > 0) then
                            do i = mom_idx%beg, mom_idx%end
                                dqR_prim_dz_vf(i)%sf(m + 1, &
                                                     iy%beg:iy%end, &
                                                     iz%beg:iz%end) = &
                                    dqL_prim_dz_vf(i)%sf(m, &
                                                         iy%beg:iy%end, &
                                                         iz%beg:iz%end)
                            end do
                        end if

                    end if

                end if

            end if
            ! END: Population of Buffers in x-direction ========================

            ! Population of Buffers in y-direction =============================
        elseif (norm_dir == 2) then

            if (bc_y%beg == -4) then    ! Riemann state extrap. BC at beginning
!$acc parallel loop collapse(3) gang vector default(present)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                    qL_prim_rsy_vf_flat(k, -1, l, i) = &
                        qR_prim_rsy_vf_flat(k, 0, l, i)
                        end do
                    end do
                end do

                if (any(Re_size > 0)) then

                    do i = mom_idx%beg, mom_idx%end
                        dqL_prim_dx_vf(i)%sf(ix%beg:ix%end, &
                                             -1, &
                                             iz%beg:iz%end) = &
                            dqR_prim_dx_vf(i)%sf(ix%beg:ix%end, &
                                                 0, &
                                                 iz%beg:iz%end)
                        if (n > 0) then
                            dqL_prim_dy_vf(i)%sf(ix%beg:ix%end, &
                                                 -1, &
                                                 iz%beg:iz%end) = &
                                dqR_prim_dy_vf(i)%sf(ix%beg:ix%end, &
                                                     0, &
                                                     iz%beg:iz%end)
                        end if
                    end do

                    if (p > 0) then
                        do i = mom_idx%beg, mom_idx%end
                            dqL_prim_dz_vf(i)%sf(ix%beg:ix%end, &
                                                 -1, &
                                                 iz%beg:iz%end) = &
                                dqR_prim_dz_vf(i)%sf(ix%beg:ix%end, &
                                                     0, &
                                                     iz%beg:iz%end)
                        end do
                    end if

                end if

            end if

            if (bc_y%end == -4) then    ! Riemann state extrap. BC at end

!$acc parallel loop collapse(3) gang vector default(present)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                    qL_prim_rsy_vf_flat(k, n + 1, l, i) = &
                        qR_prim_rsy_vf_flat(k, n, l, i)
                        end do
                    end do
                end do

                if (any(Re_size > 0)) then

                    do i = mom_idx%beg, mom_idx%end
                        dqR_prim_dx_vf(i)%sf(ix%beg:ix%end, &
                                             n + 1, &
                                             iz%beg:iz%end) = &
                            dqL_prim_dx_vf(i)%sf(ix%beg:ix%end, &
                                                 n, &
                                                 iz%beg:iz%end)
                        if (n > 0) then
                            dqR_prim_dy_vf(i)%sf(ix%beg:ix%end, &
                                                 n + 1, &
                                                 iz%beg:iz%end) = &
                                dqL_prim_dy_vf(i)%sf(ix%beg:ix%end, &
                                                     n, &
                                                     iz%beg:iz%end)
                        end if
                    end do

                    if (p > 0) then
                        do i = mom_idx%beg, mom_idx%end
                            dqR_prim_dz_vf(i)%sf(ix%beg:ix%end, &
                                                 n + 1, &
                                                 iz%beg:iz%end) = &
                                dqL_prim_dz_vf(i)%sf(ix%beg:ix%end, &
                                                     n, &
                                                     iz%beg:iz%end)
                        end do
                    end if

                end if

            end if
            ! END: Population of Buffers in y-direction ========================

            ! Population of Buffers in z-direction =============================
        else

            if (bc_z%beg == -4) then    ! Riemann state extrap. BC at beginning
!$acc parallel loop collapse(3) gang vector default(present)
                do i = 1, sys_size
                    do k = is2%beg, is2%end
                        do l = is3%beg, is3%end
                    qL_prim_rsz_vf_flat(l, k, -1, i) = &
                        qR_prim_rsz_vf_flat(l, k, 0, i)
                        end do
                    end do
                end do

                if (any(Re_size > 0)) then
                    do i = mom_idx%beg, mom_idx%end
                        dqL_prim_dx_vf(i)%sf(ix%beg:ix%end, &
                                             iy%beg:iy%end, &
                                             -1) = &
                            dqR_prim_dx_vf(i)%sf(ix%beg:ix%end, &
                                                 iy%beg:iy%end, &
                                                 0)
                        dqL_prim_dy_vf(i)%sf(ix%beg:ix%end, &
                                             iy%beg:iy%end, &
                                             -1) = &
                            dqR_prim_dy_vf(i)%sf(ix%beg:ix%end, &
                                                 iy%beg:iy%end, &
                                                 0)
                        dqL_prim_dz_vf(i)%sf(ix%beg:ix%end, &
                                             iy%beg:iy%end, &
                                             -1) = &
                            dqR_prim_dz_vf(i)%sf(ix%beg:ix%end, &
                                                 iy%beg:iy%end, &
                                                 0)
                    end do
                end if

            end if

            if (bc_z%end == -4) then    ! Riemann state extrap. BC at end

!$acc parallel loop collapse(3) gang vector default(present)
                do i = 1, sys_size
                    do k = is2%beg, is2%end
                        do l = is3%beg, is3%end
                    qL_prim_rsz_vf_flat(l, k, p + 1, i) = &
                        qR_prim_rsz_vf_flat(l, k, p, i)
                        end do
                    end do
                end do

                if (any(Re_size > 0)) then
                    do i = mom_idx%beg, mom_idx%end
                        dqR_prim_dx_vf(i)%sf(ix%beg:ix%end, &
                                             iy%beg:iy%end, &
                                             p + 1) = &
                            dqL_prim_dx_vf(i)%sf(ix%beg:ix%end, &
                                                 iy%beg:iy%end, &
                                                 p)
                        dqR_prim_dy_vf(i)%sf(ix%beg:ix%end, &
                                             iy%beg:iy%end, &
                                             p + 1) = &
                            dqL_prim_dy_vf(i)%sf(ix%beg:ix%end, &
                                                 iy%beg:iy%end, &
                                                 p)
                        dqR_prim_dz_vf(i)%sf(ix%beg:ix%end, &
                                             iy%beg:iy%end, &
                                             p + 1) = &
                            dqL_prim_dz_vf(i)%sf(ix%beg:ix%end, &
                                                 iy%beg:iy%end, &
                                                 p)
                    end do
                end if

            end if

        end if
        ! END: Population of Buffers in z-direction ========================

    end subroutine s_populate_riemann_states_variables_buffers ! -----------

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
    subroutine s_initialize_riemann_solver(&
                                           q_prim_vf, &
                                           flux_vf, flux_src_vf, &
                                           flux_gsrc_vf, &
                                           norm_dir, ix, iy, iz)


        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: flux_vf, flux_src_vf, flux_gsrc_vf

        integer, intent(IN) :: norm_dir

        type(bounds_info), intent(IN) :: ix, iy, iz

        integer :: i, j, k, l ! Generic loop iterators


        ! Reshaping Inputted Data in x-direction ===========================

            if (norm_dir == 1) then

                if (any(Re_size > 0)) then
    !$acc parallel loop collapse(4) gang vector default(present)
                                do i = momxb, E_idx
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end
                                    flux_src_vf(i)%sf(j, k, l) = 0d0
                                end do
                            end do
                        end do
                    end do
                end if

                ! ==================================================================

                ! Reshaping Inputted Data in y-direction ===========================
            elseif (norm_dir == 2) then

                if (any(Re_size > 0)) then
    !$acc parallel loop collapse(4) gang vector default(present)
                                do i = momxb, E_idx
                    do l = is3%beg, is3%end
                        do j = is1%beg, is1%end
                            do k = is2%beg, is2%end
                                    flux_src_vf(i)%sf(k, j, l) = 0d0
                                end do
                            end do
                        end do
                    end do
                end if
                ! ==================================================================

                ! Reshaping Inputted Data in z-direction ===========================
            else

                if (any(Re_size > 0)) then
    !$acc parallel loop collapse(4) gang vector default(present)
                                do i = momxb, E_idx
                    do j = is1%beg, is1%end
                        do k = is2%beg, is2%end
                            do l = is3%beg, is3%end
                                    flux_src_vf(i)%sf(l, k, j) = 0d0
                                end do
                            end do
                        end do
                    end do
                end if

            end if

        ! ==================================================================

    end subroutine s_initialize_riemann_solver ! ---------------------------

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
    subroutine s_compute_cylindrical_viscous_source_flux(velL_vf, & ! -------------
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
            intent(IN) ::         velL_vf, velR_vf, &
                          dvelL_dx_vf, dvelR_dx_vf, &
                          dvelL_dy_vf, dvelR_dy_vf, &
                          dvelL_dz_vf, dvelR_dz_vf

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: flux_src_vf

        integer, intent(IN) :: norm_dir

        type(bounds_info), intent(IN) :: ix, iy, iz

        ! Arithmetic mean of the left and right, WENO-reconstructed, cell-
        ! boundary values of cell-average first-order spatial derivatives
        ! of velocity
        real(kind(0d0)), dimension(num_dims) :: avg_vel
        real(kind(0d0)), dimension(num_dims) :: dvel_avg_dx
        real(kind(0d0)), dimension(num_dims) :: dvel_avg_dy
        real(kind(0d0)), dimension(num_dims) :: dvel_avg_dz

        ! Viscous stress tensor
        real(kind(0d0)), dimension(num_dims, num_dims) :: tau_Re

        ! Generic loop iterators
        integer :: i, j, k, l

        ! Viscous Stresses in z-direction ==================================
        if (norm_dir == 1) then

            if (Re_size(1) > 0) then              ! Shear stresses
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end

                            dvel_avg_dx(1) = 5d-1*(dvelL_dx_vf(1)%sf(j, k, l) &
                                                   + dvelR_dx_vf(1)%sf(j + 1, k, l))

                            tau_Re(1, 1) = (4d0/3d0)*dvel_avg_dx(1)/ &
                                           Re_avg_rs_vf(1)%sf(j, k, l)

                            flux_src_vf(mom_idx%beg)%sf(j, k, l) = &
                                flux_src_vf(mom_idx%beg)%sf(j, k, l) - &
                                tau_Re(1, 1)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rs_vf(1)%sf(j, k, l)* &
                                tau_Re(1, 1)

                        end do
                    end do
                end do
            end if

            if (Re_size(2) > 0) then              ! Bulk stresses
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end

                            dvel_avg_dx(1) = 5d-1*(dvelL_dx_vf(1)%sf(j, k, l) &
                                                   + dvelR_dx_vf(1)%sf(j + 1, k, l))

                            tau_Re(1, 1) = dvel_avg_dx(1)/ &
                                           Re_avg_rs_vf(2)%sf(j, k, l)

                            flux_src_vf(mom_idx%beg)%sf(j, k, l) = &
                                flux_src_vf(mom_idx%beg)%sf(j, k, l) - &
                                tau_Re(1, 1)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rs_vf(1)%sf(j, k, l)* &
                                tau_Re(1, 1)

                        end do
                    end do
                end do
            end if

            if (n == 0) return

            if (Re_size(1) > 0) then              ! Shear stresses
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end

                            avg_vel(2) = 5d-1*(velL_vf(2)%sf(j, k, l) &
                                               + velR_vf(2)%sf(j + 1, k, l))

                            do i = 1, 2
                                dvel_avg_dy(i) = &
                                    5d-1*(dvelL_dy_vf(i)%sf(j, k, l) &
                                          + dvelR_dy_vf(i)%sf(j + 1, k, l))
                            end do

                            dvel_avg_dx(2) = 5d-1*(dvelL_dx_vf(2)%sf(j, k, l) &
                                                   + dvelR_dx_vf(2)%sf(j + 1, k, l))

                            tau_Re(1, 1) = -(2d0/3d0)*(dvel_avg_dy(2) + &
                                                       avg_vel(2)/y_cc(k))/ &
                                           Re_avg_rs_vf(1)%sf(j, k, l)

                            tau_Re(1, 2) = (dvel_avg_dy(1) + dvel_avg_dx(2))/ &
                                           Re_avg_rs_vf(1)%sf(j, k, l)

                            do i = 1, 2

                                flux_src_vf(cont_idx%end + i)%sf(j, k, l) = &
                                    flux_src_vf(cont_idx%end + i)%sf(j, k, l) - &
                                    tau_Re(1, i)

                                flux_src_vf(E_idx)%sf(j, k, l) = &
                                    flux_src_vf(E_idx)%sf(j, k, l) - &
                                    vel_src_rs_vf(i)%sf(j, k, l)* &
                                    tau_Re(1, i)

                            end do

                        end do
                    end do
                end do
            end if

            if (Re_size(2) > 0) then              ! Bulk stresses
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end

                            avg_vel(2) = 5d-1*(velL_vf(2)%sf(j, k, l) &
                                               + velR_vf(2)%sf(j + 1, k, l))

                            dvel_avg_dy(2) = 5d-1*(dvelL_dy_vf(2)%sf(j, k, l) &
                                                   + dvelR_dy_vf(2)%sf(j + 1, k, l))

                            tau_Re(1, 1) = (dvel_avg_dy(2) + &
                                            avg_vel(2)/y_cc(k))/ &
                                           Re_avg_rs_vf(2)%sf(j, k, l)

                            flux_src_vf(mom_idx%beg)%sf(j, k, l) = &
                                flux_src_vf(mom_idx%beg)%sf(j, k, l) - &
                                tau_Re(1, 1)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rs_vf(1)%sf(j, k, l)* &
                                tau_Re(1, 1)

                        end do
                    end do
                end do
            end if

            if (p == 0) return

            if (Re_size(1) > 0) then              ! Shear stresses
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end

                            do i = 1, 3, 2
                                dvel_avg_dz(i) = &
                                    5d-1*(dvelL_dz_vf(i)%sf(j, k, l) &
                                          + dvelR_dz_vf(i)%sf(j + 1, k, l))
                            end do

                            dvel_avg_dx(3) = 5d-1*(dvelL_dx_vf(3)%sf(j, k, l) &
                                                   + dvelR_dx_vf(3)%sf(j + 1, k, l))

                            tau_Re(1, 1) = -(2d0/3d0)*dvel_avg_dz(3)/y_cc(k)/ &
                                           Re_avg_rs_vf(1)%sf(j, k, l)

                            tau_Re(1, 3) = (dvel_avg_dz(1)/y_cc(k) + dvel_avg_dx(3))/ &
                                           Re_avg_rs_vf(1)%sf(j, k, l)

                            do i = 1, 3, 2

                                flux_src_vf(cont_idx%end + i)%sf(j, k, l) = &
                                    flux_src_vf(cont_idx%end + i)%sf(j, k, l) - &
                                    tau_Re(1, i)

                                flux_src_vf(E_idx)%sf(j, k, l) = &
                                    flux_src_vf(E_idx)%sf(j, k, l) - &
                                    vel_src_rs_vf(i)%sf(j, k, l)* &
                                    tau_Re(1, i)

                            end do

                        end do
                    end do
                end do
            end if

            if (Re_size(2) > 0) then              ! Bulk stresses
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end

                            dvel_avg_dz(3) = 5d-1*(dvelL_dz_vf(3)%sf(j, k, l) &
                                                   + dvelR_dz_vf(3)%sf(j + 1, k, l))

                            tau_Re(1, 1) = dvel_avg_dz(3)/y_cc(k)/ &
                                           Re_avg_rs_vf(2)%sf(j, k, l)

                            flux_src_vf(mom_idx%beg)%sf(j, k, l) = &
                                flux_src_vf(mom_idx%beg)%sf(j, k, l) - &
                                tau_Re(1, 1)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rs_vf(1)%sf(j, k, l)* &
                                tau_Re(1, 1)

                        end do
                    end do
                end do
            end if
            ! END: Viscous Stresses in z-direction =============================

            ! Viscous Stresses in r-direction ==================================
        elseif (norm_dir == 2) then

            if (Re_size(1) > 0) then              ! Shear stresses
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end

                            avg_vel(2) = 5d-1*(velL_vf(2)%sf(j, k, l) &
                                               + velR_vf(2)%sf(j, k + 1, l))

                            do i = 1, 2

                                dvel_avg_dx(i) = &
                                    5d-1*(dvelL_dx_vf(i)%sf(j, k, l) &
                                          + dvelR_dx_vf(i)%sf(j, k + 1, l))

                                dvel_avg_dy(i) = &
                                    5d-1*(dvelL_dy_vf(i)%sf(j, k, l) &
                                          + dvelR_dy_vf(i)%sf(j, k + 1, l))

                            end do

                            tau_Re(2, 1) = (dvel_avg_dy(1) + dvel_avg_dx(2))/ &
                                           Re_avg_rs_vf(1)%sf(k, j, l)

                            tau_Re(2, 2) = (4d0*dvel_avg_dy(2) &
                                            - 2d0*dvel_avg_dx(1) &
                                            - 2d0*avg_vel(2)/y_cb(k))/ &
                                           (3d0*Re_avg_rs_vf(1)%sf(k, j, l))

                            do i = 1, 2

                                flux_src_vf(cont_idx%end + i)%sf(j, k, l) = &
                                    flux_src_vf(cont_idx%end + i)%sf(j, k, l) - &
                                    tau_Re(2, i)

                                flux_src_vf(E_idx)%sf(j, k, l) = &
                                    flux_src_vf(E_idx)%sf(j, k, l) - &
                                    vel_src_rs_vf(i)%sf(k, j, l)* &
                                    tau_Re(2, i)

                            end do

                        end do
                    end do
                end do
            end if

            if (Re_size(2) > 0) then              ! Bulk stresses
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end

                            avg_vel(2) = 5d-1*(velL_vf(2)%sf(j, k, l) &
                                               + velR_vf(2)%sf(j, k + 1, l))

                            dvel_avg_dx(1) = 5d-1*(dvelL_dx_vf(1)%sf(j, k, l) &
                                                   + dvelR_dx_vf(1)%sf(j, k + 1, l))

                            dvel_avg_dy(2) = 5d-1*(dvelL_dy_vf(2)%sf(j, k, l) &
                                                   + dvelR_dy_vf(2)%sf(j, k + 1, l))

                            tau_Re(2, 2) = (dvel_avg_dx(1) + dvel_avg_dy(2) + &
                                            avg_vel(2)/y_cb(k))/ &
                                           Re_avg_rs_vf(2)%sf(k, j, l)

                            flux_src_vf(mom_idx%beg + 1)%sf(j, k, l) = &
                                flux_src_vf(mom_idx%beg + 1)%sf(j, k, l) - &
                                tau_Re(2, 2)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rs_vf(2)%sf(k, j, l)* &
                                tau_Re(2, 2)

                        end do
                    end do
                end do
            end if

            if (p == 0) return

            if (Re_size(1) > 0) then              ! Shear stresses
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end

                            avg_vel(3) = 5d-1*(velL_vf(3)%sf(j, k, l) &
                                               + velR_vf(3)%sf(j, k + 1, l))

                            do i = 2, 3
                                dvel_avg_dz(i) = &
                                    5d-1*(dvelL_dz_vf(i)%sf(j, k, l) &
                                          + dvelR_dz_vf(i)%sf(j, k + 1, l))
                            end do

                            dvel_avg_dy(3) = 5d-1*(dvelL_dy_vf(3)%sf(j, k, l) &
                                                   + dvelR_dy_vf(3)%sf(j, k + 1, l))

                            tau_Re(2, 2) = -(2d0/3d0)*dvel_avg_dz(3)/y_cb(k)/ &
                                           Re_avg_rs_vf(1)%sf(k, j, l)

                            tau_Re(2, 3) = ((dvel_avg_dz(2) - avg_vel(3))/ &
                                            y_cb(k) + dvel_avg_dy(3))/ &
                                           Re_avg_rs_vf(1)%sf(k, j, l)

                            do i = 2, 3

                                flux_src_vf(cont_idx%end + i)%sf(j, k, l) = &
                                    flux_src_vf(cont_idx%end + i)%sf(j, k, l) - &
                                    tau_Re(2, i)

                                flux_src_vf(E_idx)%sf(j, k, l) = &
                                    flux_src_vf(E_idx)%sf(j, k, l) - &
                                    vel_src_rs_vf(i)%sf(k, j, l)* &
                                    tau_Re(2, i)

                            end do

                        end do
                    end do
                end do
            end if

            if (Re_size(2) > 0) then              ! Bulk stresses
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end

                            dvel_avg_dz(3) = 5d-1*(dvelL_dz_vf(3)%sf(j, k, l) &
                                                   + dvelR_dz_vf(3)%sf(j, k + 1, l))

                            tau_Re(2, 2) = dvel_avg_dz(3)/y_cb(k)/ &
                                           Re_avg_rs_vf(2)%sf(k, j, l)

                            flux_src_vf(mom_idx%beg + 1)%sf(j, k, l) = &
                                flux_src_vf(mom_idx%beg + 1)%sf(j, k, l) - &
                                tau_Re(2, 2)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rs_vf(2)%sf(k, j, l)* &
                                tau_Re(2, 2)

                        end do
                    end do
                end do
            end if
            ! END: Viscous Stresses in r-direction =============================

            ! Viscous Stresses in theta-direction ==================================
        else

            if (Re_size(1) > 0) then              ! Shear stresses
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end

                            do i = 2, 3
                                avg_vel(i) = 5d-1*(velL_vf(i)%sf(j, k, l) &
                                                   + velR_vf(i)%sf(j, k, l + 1))
                            end do

                            do i = 1, 3, 2
                                dvel_avg_dx(i) = &
                                    5d-1*(dvelL_dx_vf(i)%sf(j, k, l) &
                                          + dvelR_dx_vf(i)%sf(j, k, l + 1))
                            end do

                            do i = 2, 3
                                dvel_avg_dy(i) = &
                                    5d-1*(dvelL_dy_vf(i)%sf(j, k, l) &
                                          + dvelR_dy_vf(i)%sf(j, k, l + 1))
                            end do

                            do i = 1, 3
                                dvel_avg_dz(i) = &
                                    5d-1*(dvelL_dz_vf(i)%sf(j, k, l) &
                                          + dvelR_dz_vf(i)%sf(j, k, l + 1))
                            end do

                            tau_Re(3, 1) = (dvel_avg_dz(1)/y_cc(k) + dvel_avg_dx(3))/ &
                                           Re_avg_rs_vf(1)%sf(l, k, j)/ &
                                           y_cc(k)

                            tau_Re(3, 2) = ((dvel_avg_dz(2) - avg_vel(3))/ &
                                            y_cc(k) + dvel_avg_dy(3))/ &
                                           Re_avg_rs_vf(1)%sf(l, k, j)/ &
                                           y_cc(k)

                            tau_Re(3, 3) = (4d0*dvel_avg_dz(3)/y_cc(k) &
                                            - 2d0*dvel_avg_dx(1) &
                                            - 2d0*dvel_avg_dy(2) &
                                            + 4d0*avg_vel(2)/y_cc(k))/ &
                                           (3d0*Re_avg_rs_vf(1)%sf(l, k, j))/ &
                                           y_cc(k)

                            do i = 1, 3

                                flux_src_vf(cont_idx%end + i)%sf(j, k, l) = &
                                    flux_src_vf(cont_idx%end + i)%sf(j, k, l) - &
                                    tau_Re(3, i)

                                flux_src_vf(E_idx)%sf(j, k, l) = &
                                    flux_src_vf(E_idx)%sf(j, k, l) - &
                                    vel_src_rs_vf(i)%sf(l, k, j)* &
                                    tau_Re(3, i)

                            end do

                        end do
                    end do
                end do
            end if

            if (Re_size(2) > 0) then              ! Bulk stresses
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end

                            avg_vel(2) = 5d-1*(velL_vf(2)%sf(j, k, l) &
                                               + velR_vf(2)%sf(j, k, l + 1))

                            dvel_avg_dx(1) = 5d-1*(dvelL_dx_vf(1)%sf(j, k, l) &
                                                   + dvelR_dx_vf(1)%sf(j, k, l + 1))

                            dvel_avg_dy(2) = 5d-1*(dvelL_dy_vf(2)%sf(j, k, l) &
                                                   + dvelR_dy_vf(2)%sf(j, k, l + 1))

                            dvel_avg_dz(3) = 5d-1*(dvelL_dz_vf(3)%sf(j, k, l) &
                                                   + dvelR_dz_vf(3)%sf(j, k, l + 1))

                            tau_Re(3, 3) = (dvel_avg_dx(1) &
                                            + dvel_avg_dy(2) &
                                            + dvel_avg_dz(3)/y_cc(k) &
                                            + avg_vel(2)/y_cc(k))/ &
                                           Re_avg_rs_vf(2)%sf(l, k, j)/ &
                                           y_cc(k)

                            flux_src_vf(mom_idx%end)%sf(j, k, l) = &
                                flux_src_vf(mom_idx%end)%sf(j, k, l) - &
                                tau_Re(3, 3)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rs_vf(3)%sf(l, k, j)* &
                                tau_Re(3, 3)

                        end do
                    end do
                end do
            end if

        end if
        ! END: Viscous Stresses in theta-direction =============================

    end subroutine s_compute_cylindrical_viscous_source_flux ! -------------------------

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
    subroutine s_compute_cartesian_viscous_source_flux(velL_vf, & ! -------------
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
            intent(IN) ::         velL_vf, velR_vf, &
                          dvelL_dx_vf, dvelR_dx_vf, &
                          dvelL_dy_vf, dvelR_dy_vf, &
                          dvelL_dz_vf, dvelR_dz_vf

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: flux_src_vf

        integer, intent(IN) :: norm_dir

        type(bounds_info), intent(IN) :: ix, iy, iz

        ! Arithmetic mean of the left and right, WENO-reconstructed, cell-
        ! boundary values of cell-average first-order spatial derivatives
        ! of velocity
        real(kind(0d0)), dimension(num_dims) :: dvel_avg_dx
        real(kind(0d0)), dimension(num_dims) :: dvel_avg_dy
        real(kind(0d0)), dimension(num_dims) :: dvel_avg_dz

        real(kind(0d0)), dimension(num_dims, num_dims) :: tau_Re !< Viscous stress tensor

        integer :: i, j, k, l !< Generic loop iterators

        ! Viscous Stresses in x-direction ==================================
        if (norm_dir == 1) then

            if (Re_size(1) > 0) then              ! Shear stresses
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end

                            dvel_avg_dx(1) = 5d-1*(dvelL_dx_vf(1)%sf(j, k, l) &
                                                   + dvelR_dx_vf(1)%sf(j + 1, k, l))

                            tau_Re(1, 1) = (4d0/3d0)*dvel_avg_dx(1)/ &
                                           Re_avg_rs_vf(1)%sf(j, k, l)

                            flux_src_vf(mom_idx%beg)%sf(j, k, l) = &
                                flux_src_vf(mom_idx%beg)%sf(j, k, l) - &
                                tau_Re(1, 1)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rs_vf(1)%sf(j, k, l)* &
                                tau_Re(1, 1)

                        end do
                    end do
                end do
            end if

            if (Re_size(2) > 0) then              ! Bulk stresses
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end

                            dvel_avg_dx(1) = 5d-1*(dvelL_dx_vf(1)%sf(j, k, l) &
                                                   + dvelR_dx_vf(1)%sf(j + 1, k, l))

                            tau_Re(1, 1) = dvel_avg_dx(1)/ &
                                           Re_avg_rs_vf(2)%sf(j, k, l)

                            flux_src_vf(mom_idx%beg)%sf(j, k, l) = &
                                flux_src_vf(mom_idx%beg)%sf(j, k, l) - &
                                tau_Re(1, 1)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rs_vf(1)%sf(j, k, l)* &
                                tau_Re(1, 1)

                        end do
                    end do
                end do
            end if

            if (n == 0) return

            if (Re_size(1) > 0) then              ! Shear stresses
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end

                            do i = 1, 2
                                dvel_avg_dy(i) = &
                                    5d-1*(dvelL_dy_vf(i)%sf(j, k, l) &
                                          + dvelR_dy_vf(i)%sf(j + 1, k, l))
                            end do

                            dvel_avg_dx(2) = 5d-1*(dvelL_dx_vf(2)%sf(j, k, l) &
                                                   + dvelR_dx_vf(2)%sf(j + 1, k, l))

                            tau_Re(1, 1) = -(2d0/3d0)*dvel_avg_dy(2)/ &
                                           Re_avg_rs_vf(1)%sf(j, k, l)

                            tau_Re(1, 2) = (dvel_avg_dy(1) + dvel_avg_dx(2))/ &
                                           Re_avg_rs_vf(1)%sf(j, k, l)

                            do i = 1, 2

                                flux_src_vf(cont_idx%end + i)%sf(j, k, l) = &
                                    flux_src_vf(cont_idx%end + i)%sf(j, k, l) - &
                                    tau_Re(1, i)

                                flux_src_vf(E_idx)%sf(j, k, l) = &
                                    flux_src_vf(E_idx)%sf(j, k, l) - &
                                    vel_src_rs_vf(i)%sf(j, k, l)* &
                                    tau_Re(1, i)

                            end do

                        end do
                    end do
                end do
            end if

            if (Re_size(2) > 0) then              ! Bulk stresses
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end

                            dvel_avg_dy(2) = 5d-1*(dvelL_dy_vf(2)%sf(j, k, l) &
                                                   + dvelR_dy_vf(2)%sf(j + 1, k, l))

                            tau_Re(1, 1) = dvel_avg_dy(2)/ &
                                           Re_avg_rs_vf(2)%sf(j, k, l)

                            flux_src_vf(mom_idx%beg)%sf(j, k, l) = &
                                flux_src_vf(mom_idx%beg)%sf(j, k, l) - &
                                tau_Re(1, 1)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rs_vf(1)%sf(j, k, l)* &
                                tau_Re(1, 1)

                        end do
                    end do
                end do
            end if

            if (p == 0) return

            if (Re_size(1) > 0) then              ! Shear stresses
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end

                            do i = 1, 3, 2
                                dvel_avg_dz(i) = &
                                    5d-1*(dvelL_dz_vf(i)%sf(j, k, l) &
                                          + dvelR_dz_vf(i)%sf(j + 1, k, l))
                            end do

                            dvel_avg_dx(3) = 5d-1*(dvelL_dx_vf(3)%sf(j, k, l) &
                                                   + dvelR_dx_vf(3)%sf(j + 1, k, l))

                            tau_Re(1, 1) = -(2d0/3d0)*dvel_avg_dz(3)/ &
                                           Re_avg_rs_vf(1)%sf(j, k, l)

                            tau_Re(1, 3) = (dvel_avg_dz(1) + dvel_avg_dx(3))/ &
                                           Re_avg_rs_vf(1)%sf(j, k, l)

                            do i = 1, 3, 2

                                flux_src_vf(cont_idx%end + i)%sf(j, k, l) = &
                                    flux_src_vf(cont_idx%end + i)%sf(j, k, l) - &
                                    tau_Re(1, i)

                                flux_src_vf(E_idx)%sf(j, k, l) = &
                                    flux_src_vf(E_idx)%sf(j, k, l) - &
                                    vel_src_rs_vf(i)%sf(j, k, l)* &
                                    tau_Re(1, i)

                            end do

                        end do
                    end do
                end do
            end if

            if (Re_size(2) > 0) then              ! Bulk stresses
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end

                            dvel_avg_dz(3) = 5d-1*(dvelL_dz_vf(3)%sf(j, k, l) &
                                                   + dvelR_dz_vf(3)%sf(j + 1, k, l))

                            tau_Re(1, 1) = dvel_avg_dz(3)/ &
                                           Re_avg_rs_vf(2)%sf(j, k, l)

                            flux_src_vf(mom_idx%beg)%sf(j, k, l) = &
                                flux_src_vf(mom_idx%beg)%sf(j, k, l) - &
                                tau_Re(1, 1)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rs_vf(1)%sf(j, k, l)* &
                                tau_Re(1, 1)

                        end do
                    end do
                end do
            end if
            ! END: Viscous Stresses in x-direction =============================

            ! Viscous Stresses in y-direction ==================================
        elseif (norm_dir == 2) then

            if (Re_size(1) > 0) then              ! Shear stresses
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end

                            do i = 1, 2

                                dvel_avg_dx(i) = &
                                    5d-1*(dvelL_dx_vf(i)%sf(j, k, l) &
                                          + dvelR_dx_vf(i)%sf(j, k + 1, l))

                                dvel_avg_dy(i) = &
                                    5d-1*(dvelL_dy_vf(i)%sf(j, k, l) &
                                          + dvelR_dy_vf(i)%sf(j, k + 1, l))

                            end do

                            tau_Re(2, 1) = (dvel_avg_dy(1) + dvel_avg_dx(2))/ &
                                           Re_avg_rs_vf(1)%sf(k, j, l)

                            tau_Re(2, 2) = (4d0*dvel_avg_dy(2) &
                                            - 2d0*dvel_avg_dx(1))/ &
                                           (3d0*Re_avg_rs_vf(1)%sf(k, j, l))

                            do i = 1, 2

                                flux_src_vf(cont_idx%end + i)%sf(j, k, l) = &
                                    flux_src_vf(cont_idx%end + i)%sf(j, k, l) - &
                                    tau_Re(2, i)

                                flux_src_vf(E_idx)%sf(j, k, l) = &
                                    flux_src_vf(E_idx)%sf(j, k, l) - &
                                    vel_src_rs_vf(i)%sf(k, j, l)* &
                                    tau_Re(2, i)

                            end do

                        end do
                    end do
                end do
            end if

            if (Re_size(2) > 0) then              ! Bulk stresses
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end

                            dvel_avg_dx(1) = 5d-1*(dvelL_dx_vf(1)%sf(j, k, l) &
                                                   + dvelR_dx_vf(1)%sf(j, k + 1, l))

                            dvel_avg_dy(2) = 5d-1*(dvelL_dy_vf(2)%sf(j, k, l) &
                                                   + dvelR_dy_vf(2)%sf(j, k + 1, l))

                            tau_Re(2, 2) = (dvel_avg_dx(1) + dvel_avg_dy(2))/ &
                                           Re_avg_rs_vf(2)%sf(k, j, l)

                            flux_src_vf(mom_idx%beg + 1)%sf(j, k, l) = &
                                flux_src_vf(mom_idx%beg + 1)%sf(j, k, l) - &
                                tau_Re(2, 2)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rs_vf(2)%sf(k, j, l)* &
                                tau_Re(2, 2)

                        end do
                    end do
                end do
            end if

            if (p == 0) return

            if (Re_size(1) > 0) then              ! Shear stresses
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end

                            do i = 2, 3
                                dvel_avg_dz(i) = &
                                    5d-1*(dvelL_dz_vf(i)%sf(j, k, l) &
                                          + dvelR_dz_vf(i)%sf(j, k + 1, l))
                            end do

                            dvel_avg_dy(3) = 5d-1*(dvelL_dy_vf(3)%sf(j, k, l) &
                                                   + dvelR_dy_vf(3)%sf(j, k + 1, l))

                            tau_Re(2, 2) = -(2d0/3d0)*dvel_avg_dz(3)/ &
                                           Re_avg_rs_vf(1)%sf(k, j, l)

                            tau_Re(2, 3) = (dvel_avg_dz(2) + dvel_avg_dy(3))/ &
                                           Re_avg_rs_vf(1)%sf(k, j, l)

                            do i = 2, 3

                                flux_src_vf(cont_idx%end + i)%sf(j, k, l) = &
                                    flux_src_vf(cont_idx%end + i)%sf(j, k, l) - &
                                    tau_Re(2, i)

                                flux_src_vf(E_idx)%sf(j, k, l) = &
                                    flux_src_vf(E_idx)%sf(j, k, l) - &
                                    vel_src_rs_vf(i)%sf(k, j, l)* &
                                    tau_Re(2, i)

                            end do

                        end do
                    end do
                end do
            end if

            if (Re_size(2) > 0) then              ! Bulk stresses
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end

                            dvel_avg_dz(3) = 5d-1*(dvelL_dz_vf(3)%sf(j, k, l) &
                                                   + dvelR_dz_vf(3)%sf(j, k + 1, l))

                            tau_Re(2, 2) = dvel_avg_dz(3)/ &
                                           Re_avg_rs_vf(2)%sf(k, j, l)

                            flux_src_vf(mom_idx%beg + 1)%sf(j, k, l) = &
                                flux_src_vf(mom_idx%beg + 1)%sf(j, k, l) - &
                                tau_Re(2, 2)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rs_vf(2)%sf(k, j, l)* &
                                tau_Re(2, 2)

                        end do
                    end do
                end do
            end if
            ! END: Viscous Stresses in y-direction =============================

            ! Viscous Stresses in z-direction ==================================
        else

            if (Re_size(1) > 0) then              ! Shear stresses
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end

                            do i = 1, 3, 2
                                dvel_avg_dx(i) = &
                                    5d-1*(dvelL_dx_vf(i)%sf(j, k, l) &
                                          + dvelR_dx_vf(i)%sf(j, k, l + 1))
                            end do

                            do i = 2, 3
                                dvel_avg_dy(i) = &
                                    5d-1*(dvelL_dy_vf(i)%sf(j, k, l) &
                                          + dvelR_dy_vf(i)%sf(j, k, l + 1))
                            end do

                            do i = 1, 3
                                dvel_avg_dz(i) = &
                                    5d-1*(dvelL_dz_vf(i)%sf(j, k, l) &
                                          + dvelR_dz_vf(i)%sf(j, k, l + 1))
                            end do

                            tau_Re(3, 1) = (dvel_avg_dz(1) + dvel_avg_dx(3))/ &
                                           Re_avg_rs_vf(1)%sf(l, k, j)

                            tau_Re(3, 2) = (dvel_avg_dz(2) + dvel_avg_dy(3))/ &
                                           Re_avg_rs_vf(1)%sf(l, k, j)

                            tau_Re(3, 3) = (4d0*dvel_avg_dz(3) &
                                            - 2d0*dvel_avg_dx(1) &
                                            - 2d0*dvel_avg_dy(2))/ &
                                           (3d0*Re_avg_rs_vf(1)%sf(l, k, j))

                            do i = 1, 3

                                flux_src_vf(cont_idx%end + i)%sf(j, k, l) = &
                                    flux_src_vf(cont_idx%end + i)%sf(j, k, l) - &
                                    tau_Re(3, i)

                                flux_src_vf(E_idx)%sf(j, k, l) = &
                                    flux_src_vf(E_idx)%sf(j, k, l) - &
                                    vel_src_rs_vf(i)%sf(l, k, j)* &
                                    tau_Re(3, i)

                            end do

                        end do
                    end do
                end do
            end if

            if (Re_size(2) > 0) then              ! Bulk stresses
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end

                            dvel_avg_dx(1) = 5d-1*(dvelL_dx_vf(1)%sf(j, k, l) &
                                                   + dvelR_dx_vf(1)%sf(j, k, l + 1))

                            dvel_avg_dy(2) = 5d-1*(dvelL_dy_vf(2)%sf(j, k, l) &
                                                   + dvelR_dy_vf(2)%sf(j, k, l + 1))

                            dvel_avg_dz(3) = 5d-1*(dvelL_dz_vf(3)%sf(j, k, l) &
                                                   + dvelR_dz_vf(3)%sf(j, k, l + 1))

                            tau_Re(3, 3) = (dvel_avg_dx(1) &
                                            + dvel_avg_dy(2) &
                                            + dvel_avg_dz(3))/ &
                                           Re_avg_rs_vf(2)%sf(l, k, j)

                            flux_src_vf(mom_idx%end)%sf(j, k, l) = &
                                flux_src_vf(mom_idx%end)%sf(j, k, l) - &
                                tau_Re(3, 3)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rs_vf(3)%sf(l, k, j)* &
                                tau_Re(3, 3)

                        end do
                    end do
                end do
            end if

        end if
        ! END: Viscous Stresses in z-direction =============================

    end subroutine s_compute_cartesian_viscous_source_flux ! -------------------------


    !>  Deallocation and/or disassociation procedures that are
        !!      needed to finalize the selected Riemann problem solver
        !!  @param flux_vf       Intercell fluxes
        !!  @param flux_src_vf   Intercell source fluxes
        !!  @param flux_gsrc_vf  Intercell geometric source fluxes
        !!  @param norm_dir Dimensional splitting coordinate direction
        !!  @param ix   Index bounds in  first coordinate direction
        !!  @param iy   Index bounds in second coordinate direction
        !!  @param iz   Index bounds in  third coordinate direction
    subroutine s_finalize_riemann_solver(flux_vf, flux_src_vf, & ! --------
                                         flux_gsrc_vf, &
                                         norm_dir, ix, iy, iz)

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: flux_vf, flux_src_vf, flux_gsrc_vf

        integer, intent(IN) :: norm_dir

        type(bounds_info), intent(IN) :: ix, iy, iz

        integer :: i, j, k, l !< Generic loop iterators




        ! Reshaping Outputted Data in y-direction ==========================


            if (norm_dir == 2) then
!$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do j = is1%beg, is1%end
                            do k = is2%beg, is2%end
                                flux_vf(i)%sf(k, j, l) = &
                                    flux_rsy_vf_flat(j, k, l, i)
                            end do
                        end do
                    end do
                end do

                if(cyl_coord) then
!$acc parallel loop collapse(4) gang vector default(present)
                    do i = 1, sys_size
                      do l = is3%beg, is3%end
                            do j = is1%beg, is1%end
                                do k = is2%beg, is2%end
                                        flux_gsrc_vf(i)%sf(k, j, l) = &
                                            flux_gsrc_rsy_vf_flat(j, k, l, i)
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
                                flux_src_rsy_vf_flat(j, k, l, advxb)
                        end do
                    end do
                end do

                if (riemann_solver == 1) then
    !$acc parallel loop collapse(4) gang vector default(present)
                    do i = advxb + 1, sys_size
                        do l = is3%beg, is3%end
                            do j = is1%beg, is1%end
                                do k = is2%beg, is2%end
                                    flux_src_vf(i)%sf(k, j, l) = &
                                        flux_src_rsy_vf_flat(j, k, l, i)
                                end do
                            end do
                        end do
                    end do

                end if
                ! ==================================================================

                ! Reshaping Outputted Data in z-direction ==========================
            elseif (norm_dir == 3) then
    !$acc parallel loop collapse(4) gang vector default(present)
               do i = 1, sys_size
                    do j = is1%beg, is1%end
                        do k = is2%beg, is2%end
                            do l = is3%beg, is3%end

                            flux_vf(i)%sf(l, k, j) = &
                                flux_rsz_vf_flat(j, k, l, i)
                            end do
                        end do
                    end do
                end do
                if(grid_geometry == 3) then
    !$acc parallel loop collapse(4) gang vector default(present)
                   do i = 1, sys_size
                        do j = is1%beg, is1%end
                            do k = is2%beg, is2%end
                                do l = is3%beg, is3%end

                                    flux_gsrc_vf(i)%sf(l, k, j) = &
                                        flux_gsrc_rsz_vf_flat(j, k, l, i)
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
                                flux_src_rsz_vf_flat(j, k, l, advxb)
                        end do
                    end do
                end do

                if (riemann_solver == 1) then
    !$acc parallel loop collapse(4) gang vector default(present)
                    do i = advxb + 1, sys_size
                        do j = is1%beg, is1%end
                            do k = is2%beg, is2%end
                                 do l = is3%beg, is3%end
                                    flux_src_vf(i)%sf(l, k, j) = &
                                        flux_src_rsz_vf_flat(j, k, l, i)
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
                                    flux_rsx_vf_flat(j, k, l, i)
                            end do
                        end do
                    end do
                end do

    !$acc parallel loop collapse(3) gang vector default(present)
                do l = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = is1%beg, is1%end
                    flux_src_vf(advxb)%sf(j, k, l) = &
                        flux_src_rsx_vf_flat(j, k, l, advxb)
                        end do
                    end do
                end do

                if (riemann_solver == 1) then
    !$acc parallel loop collapse(4) gang vector default(present)
                do i = advxb + 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end
                                flux_src_vf(i)%sf(j, k, l) = &
                                    flux_src_rsx_vf_flat(j, k, l, i)
                                end do
                            end do
                        end do
                    end do
                end if
            end if




        ! ==================================================================


        ! ==================================================================

    end subroutine s_finalize_riemann_solver ! -----------------------------

    !> Module deallocation and/or disassociation procedures
    subroutine s_finalize_riemann_solvers_module() ! -----------------------

        ! Deallocating the variables that were utilized to formulate the
        ! left, right and average states of the Riemann problem, as well
        ! the Riemann problem solution

        integer :: i

        deallocate (alpha_rho_L, vel_L)
        deallocate (alpha_rho_R, vel_R)

        deallocate (vel_avg)

        deallocate (alpha_L, alpha_R)

        if (any(Re_size > 0)) deallocate (Re_avg_rs_vf)

        if (riemann_solver == 3) then
            deallocate (alpha_rho_IC, vel_IC)
            deallocate (alpha_IC)
        end if

        if (bubbles) then
            if (qbmm) then
                deallocate (moms_L, moms_R)
            end if
            deallocate (R0_L, R0_R, pbw_L, pbw_R)
            deallocate (V0_L, V0_R)
        end if

        deallocate(gammas, pi_infs)
        ! Disassociating procedural pointer to the subroutine which was
        ! utilized to calculate the solution of a given Riemann problem
        s_riemann_solver => null()

        ! Disassociating the procedural pointers to the procedures that were
        ! utilized to compute the average state and estimate the wave speeds
        s_compute_average_state => null(); s_compute_wave_speeds => null()

        ! Disassociating procedural pointer to the subroutine which was
        ! utilized to calculate the viscous source flux
        s_compute_viscous_source_flux => null()

        ! Disassociating the pointer to the procedure that was utilized to
        ! to convert mixture or species variables to the mixture variables
        s_convert_to_mixture_variables => null()


        if(Re_size(1) > 0) then
            deallocate(Re_avg_rsx_vf)
        end if
        deallocate(vel_src_rsx_vf_flat)
        deallocate(flux_rsx_vf_flat)
        deallocate(flux_src_rsx_vf_flat)
        deallocate(flux_gsrc_rsx_vf_flat)
        !deallocate(qL_prim_rsx_vf_flat)
        !deallocate(qR_prim_rsx_vf_flat)


        if(n == 0) return


        if(Re_size(1) > 0) then
            deallocate(Re_avg_rsy_vf_flat)
        end if
        deallocate(vel_src_rsy_vf_flat)
        deallocate(flux_rsy_vf_flat)
        deallocate(flux_src_rsy_vf_flat)
        deallocate(flux_gsrc_rsy_vf_flat)
        !deallocate(qL_prim_rsy_vf_flat)
        !deallocate(qR_prim_rsy_vf_flat)





        if(p == 0) return


        if(Re_size(1) > 0) then
            deallocate(Re_avg_rsz_vf_flat)
        end if
        deallocate(vel_src_rsz_vf_flat)
        deallocate(flux_rsz_vf_flat)
        deallocate(flux_src_rsz_vf_flat)
        deallocate(flux_gsrc_rsz_vf_flat)
        !deallocate(qL_prim_rsz_vf_flat)
        !deallocate(qR_prim_rsz_vf_flat)





    end subroutine s_finalize_riemann_solvers_module ! ---------------------

end module m_riemann_solvers
