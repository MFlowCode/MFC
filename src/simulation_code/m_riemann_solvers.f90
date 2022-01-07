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
MODULE m_riemann_solvers
    
    
    ! Dependencies =============================================================
    USE m_derived_types        !< Definitions of the derived types
    
    USE m_global_parameters    !< Definitions of the global parameters
    
    USE m_mpi_proxy            !< Message passing interface (MPI) module proxy
    
    USE m_variables_conversion !< State variables type conversion procedures

    USE m_bubbles              !< To get the bubble wall pressure function
    ! ==========================================================================
    
    IMPLICIT NONE
    
    PRIVATE; PUBLIC :: s_initialize_riemann_solvers_module, &
                       s_riemann_solver                   , &
                       s_hll_riemann_solver               , &
                       s_hllc_riemann_solver              , &
                       s_exact_riemann_solver             , &
                       s_finalize_riemann_solvers_module
    
    
    ABSTRACT INTERFACE ! =======================================================
        
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
        !!  @param kappaL_vf  Left averaged vol. frac. curvatures
        !!  @param kappaR_vf Right averaged vol. frac. curvatures
        !!  @param flux_vf Intra-cell fluxes
        !!  @param flux_src_vf Intra-cell fluxes sources
        !!  @param flux_gsrc_vf Intra-cell geometric fluxes sources
        !!  @param norm_dir Dir. splitting direction
        !!  @param ix Index bounds in the x-dir
        !!  @param iy Index bounds in the y-dir
        !!  @param iz Index bounds in the z-dir
        !!  @param q_prim_vf Cell-averaged primitive variables
        SUBROUTINE s_abstract_riemann_solver( qL_prim_vf, dqL_prim_dx_vf, &
                                                          dqL_prim_dy_vf, &
                                                          dqL_prim_dz_vf, &
                                                            gm_alphaL_vf, &
                                                               kappaL_vf, &
                                              qR_prim_vf, dqR_prim_dx_vf, &
                                                          dqR_prim_dy_vf, &
                                                          dqR_prim_dz_vf, &
                                                            gm_alphaR_vf, &
                                                               kappaR_vf, &
                                                               q_prim_vf, &
                                                    flux_vf, flux_src_vf, &
                                                            flux_gsrc_vf, &
                                                      norm_dir, ix,iy,iz  )
            
            IMPORT :: scalar_field, bounds_info, sys_size
            
            TYPE(scalar_field),  &
            DIMENSION(sys_size), &
            INTENT(INOUT) :: qL_prim_vf, qR_prim_vf
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q_prim_vf
        
            
            TYPE(scalar_field),        &
            ALLOCATABLE, DIMENSION(:), &
            INTENT(INOUT) :: dqL_prim_dx_vf, dqR_prim_dx_vf, &
                             dqL_prim_dy_vf, dqR_prim_dy_vf, &
                             dqL_prim_dz_vf, dqR_prim_dz_vf, &
                               gm_alphaL_vf,   gm_alphaR_vf, &
                                  kappaL_vf,      kappaR_vf
            
            TYPE(scalar_field),  &
            DIMENSION(sys_size), &
            INTENT(INOUT) :: flux_vf, flux_src_vf, flux_gsrc_vf
            
            INTEGER, INTENT(IN) :: norm_dir
            
            TYPE(bounds_info), INTENT(IN) :: ix,iy,iz
            
        END SUBROUTINE s_abstract_riemann_solver
        
        !>  The abstract interface to the subroutines that are used to calculate
        !!  the Roe and arithmetic average states. For more information refer to:
        !!      1) s_compute_roe_average_state
        !!      2) s_compute_arithmetic_average_state
        !!  @param i First coordinate location index
        !!  @param j Second coordinate location index
        !!  @param k Third coordinate location index
        SUBROUTINE s_compute_abstract_average_state(i,j,k)
            
            INTEGER, INTENT(IN) :: i,j,k
            
        END SUBROUTINE s_compute_abstract_average_state
        
        !> The abstract interface to the subroutines that are utilized to compute
        !! the wave speeds of the Riemann problem either directly or by the means
        !! of pressure-velocity estimates. For more information please refer to:
        !!      1) s_compute_direct_wave_speeds
        !!      2) s_compute_pressure_velocity_wave_speeds
        !!  @param i First coordinate location index
        !!  @param j Second coordinate location index
        !!  @param k Third coordinate location index
        SUBROUTINE s_compute_abstract_wave_speeds(i,j,k)
            
            INTEGER, INTENT(IN) :: i,j,k
            
        END SUBROUTINE s_compute_abstract_wave_speeds
        
        !> The abstract interface to the subroutines that are utilized to compute
        !! the viscous source fluxes for either Cartesian or cylindrical geometries.
        !! For more information please refer to:
        !!      1) s_compute_cartesian_viscous_source_flux
        !!      2) s_compute_cylindrical_viscous_source_flux
        SUBROUTINE s_compute_abstract_viscous_source_flux(     velL_vf, & ! -------------
                                                           dvelL_dx_vf, & 
                                                           dvelL_dy_vf, &
                                                           dvelL_dz_vf, &
                                                               velR_vf, &
                                                           dvelR_dx_vf, &
                                                           dvelR_dy_vf, &
                                                           dvelR_dz_vf, &
                                                           flux_src_vf, &
                                                              norm_dir, &
                                                              ix,iy,iz  )
            
            IMPORT :: scalar_field, bounds_info, num_dims, sys_size
            
            TYPE(scalar_field),  &
            DIMENSION(num_dims), &
            INTENT(IN) ::         velL_vf, velR_vf, &
                          dvelL_dx_vf, dvelR_dx_vf, &
                          dvelL_dy_vf, dvelR_dy_vf, &
                          dvelL_dz_vf, dvelR_dz_vf
            
            TYPE(scalar_field),  &
            DIMENSION(sys_size), &
            INTENT(INOUT) :: flux_src_vf
            
            INTEGER, INTENT(IN) :: norm_dir
            
            TYPE(bounds_info), INTENT(IN) :: ix,iy,iz

        END SUBROUTINE s_compute_abstract_viscous_source_flux

    END INTERFACE ! ============================================================
    
    
    !> The left (L) and right (R) WENO-reconstructed cell-boundary values of the
    !! cell-average primitive variables that define the left and right states of
    !! the Riemann problem. Variables qK_prim_rs_vf, K = L or R, are obtained by
    !! reshaping (RS) qK_prim_vf in a coordinate direction that is normal to the
    !! cell-boundaries along which the fluxes are to be determined.
    !> @{
    TYPE(scalar_field), ALLOCATABLE, DIMENSION(:) :: qL_prim_rs_vf
    TYPE(scalar_field), ALLOCATABLE, DIMENSION(:) :: qR_prim_rs_vf
    TYPE(scalar_field), ALLOCATABLE, DIMENSION(:) :: q_prim_rs_vf
    !> @}

    !> The left and right, WENO-reconstructed, cell-boundary values of the cell-
    !! average curvatures of the volume fractions that define the left and right
    !! states of the Riemann problem. The variables kappaK_rs_vf are obtained by
    !! reshaping kappaK_vf in a coordinate direction that is normal to the cell-
    !! boundaries on which the fluxes are to be determined.
    !> @{
    TYPE(scalar_field), ALLOCATABLE, DIMENSION(:) :: kappaL_rs_vf, kappaR_rs_vf
    !> @}

    !> The cell-boundary values of the fluxes (src - source) that are computed
    !! through the chosen Riemann problem solver, and the direct evaluation of
    !! source terms, by using the left and right states given in qK_prim_rs_vf,
    !! dqK_prim_ds_vf and kappaK_rs_vf, where ds = dx, dy or dz.
    !> @{
    TYPE(scalar_field), ALLOCATABLE, DIMENSION(:) :: flux_rs_vf, flux_src_rs_vf
    !> @}

   
    TYPE(scalar_field), ALLOCATABLE, DIMENSION(:) :: flux_gsrc_rs_vf !<
    !! The cell-boundary values of the geometrical source flux that are computed
    !! through the chosen Riemann problem solver by using the left and right
    !! states given in qK_prim_rs_vf. Currently 2D axisymmetric for inviscid only.
    
    ! The cell-boundary values of the velocity. vel_src_rs_vf is determined as
    ! part of Riemann problem solution and is used to evaluate the source flux.
    TYPE(scalar_field), ALLOCATABLE, DIMENSION(:) :: vel_src_rs_vf
    
    !> @name Left and right, WENO-reconstructed, cell-boundary values of cell-average
    !! partial densities, density, velocity, pressure, internal energy, energy, enthalpy, volume
    !! fractions, mass fractions, the specific heat ratio and liquid stiffness functions, speed
    !! of sound, shear and volume Reynolds numbers and the Weber numbers. These
    !! variables are left and right states of the Riemann problem obtained from
    !! qK_prim_rs_vf and kappaK_rs_vf.
    !> @{
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:)   :: alpha_rho_L, alpha_rho_R
    REAL(KIND(0d0))                              ::       rho_L,       rho_R
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:)   ::       vel_L,       vel_R
    REAL(KIND(0d0))                              ::      pres_L,      pres_R
    REAL(KIND(0d0))                              ::         E_L,         E_R
    REAL(KIND(0d0))                              ::         H_L,         H_R
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:)   ::     alpha_L,     alpha_R
    REAL(KIND(0d0))                              ::         Y_L,         Y_R
    REAL(KIND(0d0))                              ::     gamma_L,     gamma_R
    REAL(KIND(0d0))                              ::    pi_inf_L,    pi_inf_R
    REAL(KIND(0d0))                              ::         c_L,         c_R
    REAL(KIND(0d0)),              DIMENSION(2)   ::        Re_L,        Re_R
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:) ::        We_L,        We_R
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:)   ::     tau_e_L,     tau_e_R

    !> @}

    !> @name Left and right, WENO-reconstructed, cell-boundary values of cell-average
    !! bubble density, radius, radial velocity, pressure, wall pressure, and modified
    !! pressure. These variables are left and right states of the Riemann problem obtained from
    !! qK_prim_rs_vf and kappaK_rs_vf.
    !> @{
    REAL(KIND(0d0))                              ::       nbub_L,     nbub_R
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:)   ::         R0_L,       R0_R
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:)   ::         V0_L,       V0_R
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:)   ::         P0_L,       P0_R
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:)   ::        pbw_L,      pbw_R
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:) ::       moms_L,     moms_R
    REAL(KIND(0d0))                              ::     ptilde_L,   ptilde_R
    !> @}

    !> @name Left and right, WENO-reconstructed, cell-boundary, low-order values of cell-average
    !! partial densities, density, velocity, pressure, internal energy, energy, enthalpy, volume
    !! fractions, mass fractions, the specific heat ratio and liquid stiffness functions, speed
    !! of sound, shear and volume Reynolds numbers and the Weber numbers. These
    !! variables are left and right states of the Riemann problem obtained from
    !! qK_prim_rs_vf and kappaK_rs_vf.
    !> @{   
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:)   :: lo_alpha_rho_L, lo_alpha_rho_R
    REAL(KIND(0d0))                              ::       lo_rho_L,       lo_rho_R
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:)   ::       lo_vel_L,       lo_vel_R
    REAL(KIND(0d0))                              ::      lo_pres_L,      lo_pres_R
    REAL(KIND(0d0))                              ::         lo_E_L,         lo_E_R
    REAL(KIND(0d0))                              ::         lo_H_L,         lo_H_R
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:)   ::     lo_alpha_L,     lo_alpha_R
    REAL(KIND(0d0))                              ::     lo_gamma_L,     lo_gamma_R
    REAL(KIND(0d0))                              ::    lo_pi_inf_L,    lo_pi_inf_R
    REAL(KIND(0d0))                              ::         lo_c_L,         lo_c_R
    REAL(KIND(0d0)),              DIMENSION(2)   ::        lo_Re_L,        lo_Re_R
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:) ::        lo_We_L,        lo_We_R
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:)   ::     lo_tau_e_L,     lo_tau_e_R
    !> @}

    !> @name Gamma-related constants for use in exact Riemann solver (following Toro (1999) pp.153)
    !> @{
    REAL(KIND(0d0)) :: G1_L, G1_R
    REAL(KIND(0d0)) :: G2_L, G2_R
    REAL(KIND(0d0)) :: G3_L, G3_R
    REAL(KIND(0d0)) :: G4_L, G4_R
    REAL(KIND(0d0)) :: G5_L, G5_R
    REAL(KIND(0d0)) :: G6_L, G6_R
    REAL(KIND(0d0)) :: G7_L, G7_R
    REAL(KIND(0d0)) :: G8_L, G8_R
    !> @}

    !> @name Star region pressure and velocity
    !> @{
    REAL(KIND(0d0)) :: pres_S
    REAL(KIND(0d0)) :: vel_S
    !> @}

    !> @name Intercell solution used to calculated intercell flux
    !> @{
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:)   :: alpha_rho_IC
    REAL(KIND(0d0))                              :: rho_IC
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:)   :: vel_IC
    REAL(KIND(0d0))                              :: pres_IC
    REAL(KIND(0d0))                              :: E_IC
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:)   :: alpha_IC
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:)   :: tau_e_IC
    !> @}

    !> @name Surface tension pressure contribution
    !> @{
    REAL(KIND(0d0)) :: dpres_We
    REAL(KIND(0d0)) :: dpres_L, dpres_R
    !> @}

    !> @name Roe or arithmetic average density, velocity, enthalpy, volume fractions,
    !! specific heat ratio function, speed of sound, shear and volume Reynolds
    !! numbers, Weber numbers and curvatures, at the cell-boundaries, computed
    !! from the left and the right states of the Riemann problem
    !> @{
    REAL(KIND(0d0))                                 :: rho_avg
    REAL(KIND(0d0))   , ALLOCATABLE, DIMENSION(:)   :: vel_avg
    REAL(KIND(0d0))                                 :: H_avg
    TYPE(scalar_field), ALLOCATABLE, DIMENSION(:)   :: alpha_avg_rs_vf
    REAL(KIND(0d0))                                 :: gamma_avg
    REAL(KIND(0d0))                                 :: c_avg
    TYPE(scalar_field), ALLOCATABLE, DIMENSION(:)   :: Re_avg_rs_vf
    TYPE(scalar_field), ALLOCATABLE, DIMENSION(:,:) :: We_avg_rs_vf
    REAL(KIND(0d0))   , ALLOCATABLE, DIMENSION(:)   :: kappa_avg
    !> @}

    !> @name Left, right and star (S) region wave speeds
    !> @{
    REAL(KIND(0d0)) :: s_L, s_R, s_S
    !> @}

    !> @name Star region variables (HLLC)
    !> @{
    REAL(KIND(0d0)) :: rho_Star, E_Star, p_Star, p_K_Star
    !> @}

    !> Minus (M) and plus (P) wave speeds
    !> @{
    REAL(KIND(0d0)) :: s_M, s_P
    !> @}

    !> Minus and plus wave speeds functions
    !> @{
    REAL(KIND(0d0)) :: xi_M, xi_P
    !> @}

    !> @name Additional variables for applying a flux limiter
    !> @{
    REAL(KIND(0d0)) :: flux_lim_func
    REAL(KIND(0d0)) :: lo_s_L, lo_s_R, lo_s_S
    REAL(KIND(0d0)) :: lo_s_M, lo_s_P
    REAL(KIND(0d0)) :: tvd_s_M, tvd_s_P
    REAL(KIND(0d0)) :: lo_xi_M, lo_xi_P
    REAL(KIND(0d0)) :: tvd_xi_M, tvd_xi_P
    !> @}

    PROCEDURE(s_abstract_riemann_solver), &
    POINTER :: s_riemann_solver => NULL() !<
    !! Pointer to the procedure that is utilized to calculate either the HLL,
    !! HLLC or exact intercell fluxes, based on the choice of Riemann solver
    
    PROCEDURE(s_compute_abstract_average_state), &
    POINTER :: s_compute_average_state => NULL() !<
    !! Pointer to the subroutine utilized to calculate either the Roe or the
    !! arithmetic average state variables, based on the chosen average state
    
    PROCEDURE(s_compute_abstract_wave_speeds), &
    POINTER :: s_compute_wave_speeds => NULL() !<
    !! Pointer to the subroutine that is utilized to compute the wave speeds of
    !! the Riemann problem either directly or by the means of pressure-velocity
    !! estimates, based on the selected method of estimation of the wave speeds
    
    PROCEDURE(s_compute_abstract_viscous_source_flux), &
    POINTER :: s_compute_viscous_source_flux => NULL() !<
    !! Pointer to the subroutine that is utilized to compute the viscous source
    !! flux for either Cartesian or cylindrical geometries.
 
    !> @name Indical bounds in the s1-, s2- and s3-directions
    !> @{
    TYPE(bounds_info) :: is1,is2,is3
    !> @}
    
    CONTAINS
        
        
        !>  The procedure features the Harten, Lax and van Leer (HLL)
        !!      approximate Riemann solver algorithm, see Toro (1999) as
        !!      well as Saurel and Abgrall (1999). Viscous and capillary
        !!      effects have been included by adapting the exact Riemann
        !!      solver of Perigaud and Saurel (2005).
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
        !!  @param kappaL_vf  Left averaged vol. frac. curvatures
        !!  @param kappaR_vf Right averaged vol. frac. curvatures
        !!  @param flux_vf Intra-cell fluxes
        !!  @param flux_src_vf  Intra-cell fluxes sources
        !!  @param flux_gsrc_vf Intra-cell geometric fluxes sources
        !!  @param norm_dir Dir. splitting direction
        !!  @param ix Index bounds in the x-dir
        !!  @param iy Index bounds in the y-dir
        !!  @param iz Index bounds in the z-dir
        !!  @param q_prim_vf Cell-averaged primitive variables
        SUBROUTINE s_hll_riemann_solver( qL_prim_vf, dqL_prim_dx_vf, & ! -------
                                                     dqL_prim_dy_vf, &
                                                     dqL_prim_dz_vf, &
                                                       gm_alphaL_vf, &
                                                          kappaL_vf, &
                                         qR_prim_vf, dqR_prim_dx_vf, &
                                                     dqR_prim_dy_vf, &
                                                     dqR_prim_dz_vf, &
                                                       gm_alphaR_vf, &
                                                          kappaR_vf, &
                                                          q_prim_vf, &
                                               flux_vf, flux_src_vf, &
                                                       flux_gsrc_vf, &
                                                 norm_dir, ix,iy,iz  )
            
            
            TYPE(scalar_field),  &
            DIMENSION(sys_size), &
            INTENT(INOUT) :: qL_prim_vf, qR_prim_vf
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q_prim_vf
            
            TYPE(scalar_field),        &
            ALLOCATABLE, DIMENSION(:), &
            INTENT(INOUT) :: dqL_prim_dx_vf, dqR_prim_dx_vf, &
                             dqL_prim_dy_vf, dqR_prim_dy_vf, &
                             dqL_prim_dz_vf, dqR_prim_dz_vf, &
                               gm_alphaL_vf,   gm_alphaR_vf, &
                                  kappaL_vf,      kappaR_vf
            
            ! Intercell fluxes
            TYPE(scalar_field),  &
            DIMENSION(sys_size), &
            INTENT(INOUT) :: flux_vf, flux_src_vf, flux_gsrc_vf
            

            INTEGER, INTENT(IN) :: norm_dir
            TYPE(bounds_info), INTENT(IN) :: ix,iy,iz
            

            INTEGER :: i,j,k,l !< Generic loop iterators
           
            ! Placeholders for high- and low-order fluxes
            REAL(KIND(0d0)) :: hi_cons, hi_cons_L, hi_cons_R
            REAL(KIND(0d0)) :: lo_cons, lo_cons_L, lo_cons_R
            REAL(KIND(0d0)) :: tvd_cons, tvd_cons_L, tvd_cons_R
            REAL(KIND(0d0)) :: hi_flux_L, hi_flux_R, lo_flux_L, lo_flux_R
            REAL(KIND(0d0)) :: tvd_flux_L, tvd_flux_R
            
            ! Populating the buffers of the left and right Riemann problem
            ! states variables, based on the choice of boundary conditions
            CALL s_populate_riemann_states_variables_buffers( &
                                  qL_prim_vf, dqL_prim_dx_vf, &
                                              dqL_prim_dy_vf, &
                                              dqL_prim_dz_vf, &
                                                gm_alphaL_vf, &
                                                   kappaL_vf, &
                                  qR_prim_vf, dqR_prim_dx_vf, &
                                              dqR_prim_dy_vf, &
                                              dqR_prim_dz_vf, &
                                                gm_alphaR_vf, &
                                                   kappaR_vf, &
                                          norm_dir, ix,iy,iz  )
            
            
            ! Reshaping inputted data based on dimensional splitting direction
            CALL s_initialize_riemann_solver( qL_prim_vf, kappaL_vf, &
                                              qR_prim_vf, kappaR_vf, &
                                                          q_prim_vf, &
                                               flux_vf, flux_src_vf, &
                                                       flux_gsrc_vf, &
                                                 norm_dir, ix,iy,iz  )
            
            
            ! Computing HLL flux and source flux for Euler system of equations
            DO l = is3%beg, is3%end
                DO k = is2%beg, is2%end
                    DO j = is1%beg, is1%end
                        
                        CALL s_compute_average_state(j,k,l)
                        
                        CALL s_compute_wave_speeds(j,k,l)

                        IF (tvd_riemann_flux) THEN
                            CALL s_compute_flux_limiter(j,k,l,flux_lim_func,norm_dir)
                        END IF
                        
                        s_M = MIN(0d0,s_L); s_P = MAX(0d0,s_R)
                        
                        xi_M = (5d-1 + SIGN(5d-1,s_L)) &
                             + (5d-1 - SIGN(5d-1,s_L)) &
                             * (5d-1 + SIGN(5d-1,s_R))
                        xi_P = (5d-1 - SIGN(5d-1,s_R)) &
                             + (5d-1 - SIGN(5d-1,s_L)) &
                             * (5d-1 + SIGN(5d-1,s_R))
                        
                        IF (tvd_wave_speeds) THEN
                            lo_s_M = MIN(0d0,lo_s_L); lo_s_P = MAX(0d0,lo_s_R)
            
                            lo_xi_M = (5d-1 + SIGN(5d-1,lo_s_L)) & 
                                + (5d-1 - SIGN(5d-1,lo_s_L)) & 
                                * (5d-1 + SIGN(5d-1,lo_s_R))
                            lo_xi_P = (5d-1 - SIGN(5d-1,lo_s_R)) & 
                                + (5d-1 - SIGN(5d-1,lo_s_L)) & 
                                * (5d-1 + SIGN(5d-1,lo_s_R))
            
                            tvd_s_M = lo_s_M + flux_lim_func*(s_M - lo_s_M)
                            tvd_s_P = lo_s_P + flux_lim_func*(s_P - lo_s_P)
                            tvd_xi_M = lo_xi_M + flux_lim_func*(xi_M - lo_xi_M)
                            tvd_xi_P = lo_xi_P + flux_lim_func*(xi_P - lo_xi_P)
                        END IF

                        IF (tvd_riemann_flux) THEN
                            ! Mass
                            DO i = 1, cont_idx%end
                
                                 lo_cons   = lo_alpha_rho_L(i) - lo_alpha_rho_R(i)
                                 hi_cons   = alpha_rho_L(i) - alpha_rho_R(i)
                                tvd_cons   = lo_cons + flux_lim_func*(hi_cons - lo_cons)
                                 lo_flux_L = lo_alpha_rho_L(i)*lo_vel_L(dir_idx(1))
                                 lo_flux_R = lo_alpha_rho_R(i)*lo_vel_R(dir_idx(1))
                                 hi_flux_L = alpha_rho_L(i) * vel_L(dir_idx(1))
                                 hi_flux_R = alpha_rho_R(i) * vel_R(dir_idx(1))
                                tvd_flux_L = lo_flux_L + flux_lim_func*(hi_flux_L - lo_flux_L)
                                tvd_flux_R = lo_flux_R + flux_lim_func*(hi_flux_R - lo_flux_R)
                
                                IF (tvd_wave_speeds) THEN
                                    flux_rs_vf(i)%sf(j,k,l) = &
                                        ( tvd_s_M*tvd_flux_R - tvd_s_P*tvd_flux_L + tvd_s_M*tvd_s_P*tvd_cons) &
                                        / (tvd_s_M - tvd_s_P)
                                ELSE
                                    flux_rs_vf(i)%sf(j,k,l) = &
                                        ( s_M * tvd_flux_R - s_P * tvd_flux_L + s_M*s_P * tvd_cons) / (s_M - s_P)
                                END IF
                            END DO
                
                            ! Momentum
                            DO i = 1, num_dims
                
                                 lo_cons   = lo_rho_L*lo_vel_L(dir_idx(i)) - lo_rho_R*lo_vel_R(dir_idx(i))
                                 hi_cons   = rho_L * vel_L(dir_idx(i)) - rho_R * vel_R(dir_idx(i))
                                tvd_cons   = lo_cons + flux_lim_func*(hi_cons - lo_cons)
                                lo_flux_L = lo_rho_L*lo_vel_L(dir_idx(1))*lo_vel_L(dir_idx(i)) + dir_flg(dir_idx(i))*lo_pres_L
                                lo_flux_R = lo_rho_R*lo_vel_R(dir_idx(1))*lo_vel_R(dir_idx(i)) + dir_flg(dir_idx(i))*lo_pres_R
                                hi_flux_L = rho_L * vel_L(dir_idx(1)) * vel_L(dir_idx(i)) + dir_flg(dir_idx(i)) * pres_L
                                hi_flux_R = rho_R * vel_R(dir_idx(1)) * vel_R(dir_idx(i)) + dir_flg(dir_idx(i)) * pres_R

                                ! added elastic shear stress term if hypoelastic modeling is on
                                IF (hypoelasticity) THEN                                    
                                    lo_flux_L = lo_flux_L - lo_tau_e_L(dir_idx_tau(i))
                                    lo_flux_R = lo_flux_R - lo_tau_e_R(dir_idx_tau(i))
                                    hi_flux_L = hi_flux_L - tau_e_L(dir_idx_tau(i))
                                    hi_flux_R = hi_flux_R - tau_e_R(dir_idx_tau(i))
                                END IF
 
                                tvd_flux_L = lo_flux_L + flux_lim_func*(hi_flux_L - lo_flux_L)
                                tvd_flux_R = lo_flux_R + flux_lim_func*(hi_flux_R - lo_flux_R)
                
                                IF (tvd_wave_speeds) THEN
                                    flux_rs_vf(cont_idx%end+dir_idx(i))%sf(j,k,l) = &
                                        ( tvd_s_M*tvd_flux_R - tvd_s_P*tvd_flux_L + tvd_s_M*tvd_s_P*tvd_cons) &
                                        / (tvd_s_M - tvd_s_P)
                                ELSE
                                    flux_rs_vf(cont_idx%end+dir_idx(i))%sf(j,k,l) = &
                                        ( s_M* tvd_flux_R - s_P* tvd_flux_L + s_M*s_P* tvd_cons) / (s_M - s_P)
                                END IF
                            END DO
                
                            ! Energy
                
                             lo_cons   = lo_E_L - lo_E_R
                             hi_cons   = E_L - E_R
                            tvd_cons   = lo_cons + flux_lim_func*(hi_cons - lo_cons)
                             lo_flux_L = lo_vel_L(dir_idx(1))*(lo_E_L + lo_pres_L)
                             lo_flux_R = lo_vel_R(dir_idx(1))*(lo_E_R + lo_pres_R)
                             hi_flux_L = vel_L(dir_idx(1))*(E_L + pres_L)
                             hi_flux_R = vel_R(dir_idx(1))*(E_R + pres_R)

                            ! added elastic shear stress term if hypoelastic modeling is on
                            IF (hypoelasticity) THEN
                                DO i = 1, num_dims
                                        lo_flux_L = lo_flux_L - (lo_tau_e_L(dir_idx_tau(i)) * lo_vel_L(i))
                                        lo_flux_R = lo_flux_R - (lo_tau_e_R(dir_idx_tau(i)) * lo_vel_R(i))
                                        hi_flux_L = hi_flux_L - (tau_e_L(dir_idx_tau(i)) * vel_L(i))
                                        hi_flux_R = hi_flux_R - (tau_e_L(dir_idx_tau(i)) * vel_L(i))
                                END DO
                            END IF

                            tvd_flux_L = lo_flux_L + flux_lim_func*(hi_flux_L - lo_flux_L)
                            tvd_flux_R = lo_flux_R + flux_lim_func*(hi_flux_R - lo_flux_R) 
                
                            IF (tvd_wave_speeds) THEN
                                flux_rs_vf(E_idx)%sf(j,k,l) = &
                                    ( tvd_s_M*tvd_flux_R - tvd_s_P*tvd_flux_L + tvd_s_M*tvd_s_P*tvd_cons) &
                                    / (tvd_s_M - tvd_s_P)
                            ELSE
                                flux_rs_vf(E_idx)%sf(j,k,l) = &
                                    ( s_M* tvd_flux_R - s_P* tvd_flux_L + s_M*s_P* tvd_cons) / (s_M - s_P)
                            END IF
                
                            ! Volume fraction
                            DO i = adv_idx%beg, adv_idx%end
                
                                 lo_cons   =  q_prim_rs_vf(i)%sf(j,k,l) -  q_prim_rs_vf(i)%sf(j+1,k,l)
                                 hi_cons   = qL_prim_rs_vf(i)%sf(j,k,l) - qR_prim_rs_vf(i)%sf(j+1,k,l)
                                tvd_cons   =  lo_cons + flux_lim_func*(hi_cons - lo_cons)
                                 lo_flux_L =  q_prim_rs_vf(i)%sf( j ,k,l)
                                 lo_flux_R =  q_prim_rs_vf(i)%sf(j+1,k,l)
                                 hi_flux_L = qL_prim_rs_vf(i)%sf( j ,k,l)
                                 hi_flux_R = qR_prim_rs_vf(i)%sf(j+1,k,l)
                                tvd_flux_L =  lo_flux_L + flux_lim_func*(hi_flux_L - lo_flux_L)
                                tvd_flux_R =  lo_flux_R + flux_lim_func*(hi_flux_R - lo_flux_R)
                
                                IF (tvd_wave_speeds) THEN
                                        flux_rs_vf(i)%sf(j,k,l) = (tvd_cons*tvd_s_M*tvd_s_P) / (tvd_s_M - tvd_s_P)
                                    flux_src_rs_vf(i)%sf(j,k,l) = (tvd_s_M*tvd_flux_R - tvd_s_P*tvd_flux_L) &
                                                    / (tvd_s_M - tvd_s_P)
                                ELSE
                                        flux_rs_vf(i)%sf(j,k,l) = (tvd_cons*s_M*s_P) / (s_M - s_P)
                                    flux_src_rs_vf(i)%sf(j,k,l) = (s_M* tvd_flux_R - s_P* tvd_flux_L) / (s_M - s_P)
                                END IF
                            END DO

                            ! Stress eq
                            IF (hypoelasticity) THEN
                                DO i = 1, (num_dims*(num_dims+1))/2

                                    lo_cons   = lo_rho_L*lo_tau_e_L(i) - lo_rho_R*lo_tau_e_L(i)
                                    hi_cons   = rho_L * tau_e_L(i) - rho_R * tau_e_R(i)
                                   tvd_cons   = lo_cons + flux_lim_func*(hi_cons - lo_cons)
                                    lo_flux_L = lo_rho_L*lo_vel_L(dir_idx(1))*lo_tau_e_L(i)
                                    lo_flux_R = lo_rho_R*lo_vel_R(dir_idx(1))*lo_tau_e_R(i)
                                    hi_flux_L = rho_L * vel_L(dir_idx(1)) * tau_e_L(i)
                                    hi_flux_R = rho_R * vel_R(dir_idx(1)) * tau_e_R(i)
                                   tvd_flux_L = lo_flux_L + flux_lim_func*(hi_flux_L - lo_flux_L)
                                   tvd_flux_R = lo_flux_R + flux_lim_func*(hi_flux_R - lo_flux_R)

                                    IF (tvd_wave_speeds) THEN
                                        flux_rs_vf(stress_idx%beg-1+i)%sf(j,k,l) = &
                                            ( tvd_s_M*tvd_flux_R - tvd_s_P*tvd_flux_L + tvd_s_M*tvd_s_P*tvd_cons) &
                                            / (tvd_s_M - tvd_s_P)
                                    ELSE
                                        flux_rs_vf(stress_idx%beg-1+i)%sf(j,k,l) = &
                                            ( s_M * tvd_flux_R - s_P * tvd_flux_L + s_M*s_P * tvd_cons) / (s_M - s_P)
                                    END IF
                                END DO
                            END IF
                
                            ! Source terms
                            DO i = 1, num_dims
                                IF (tvd_wave_speeds) THEN
                                     lo_cons_L = lo_rho_L*(lo_s_L-lo_vel_L(dir_idx(1)))
                                     lo_cons_R = lo_rho_R*(lo_s_R-lo_vel_R(dir_idx(1)))
                                     hi_cons_L = rho_L*(s_L - vel_L(dir_idx(1)))
                                     hi_cons_R = rho_R*(s_R - vel_R(dir_idx(1)))
                                     lo_flux_L = lo_rho_L*lo_vel_L(dir_idx(i)) * &
                                         (lo_s_L-lo_vel_L(dir_idx(1)))-lo_pres_L*dir_flg(dir_idx(i))
                                     lo_flux_R = lo_rho_R*lo_vel_R(dir_idx(i)) * &
                                         (lo_s_R-lo_vel_R(dir_idx(1)))-lo_pres_R*dir_flg(dir_idx(i))
                                     hi_flux_L = rho_L * vel_L(dir_idx(i)) * &
                                         (s_L - vel_L(dir_idx(1))) - pres_L*dir_flg(dir_idx(i))
                                     hi_flux_R = rho_R * vel_R(dir_idx(i)) * &
                                         (s_R - vel_R(dir_idx(1))) - pres_R*dir_flg(dir_idx(i))
                                ELSE
                                     lo_cons_L = lo_rho_L*(s_L-lo_vel_L(dir_idx(1)))
                                     lo_cons_R = lo_rho_R*(s_R-lo_vel_R(dir_idx(1)))
                                     hi_cons_L = rho_L*(s_L - vel_L(dir_idx(1)))
                                     hi_cons_R = rho_R*(s_R - vel_R(dir_idx(1)))
                                     lo_flux_L = lo_rho_L*lo_vel_L(dir_idx(i))*(s_L-lo_vel_L(dir_idx(1)))-lo_pres_L * &
                                         dir_flg(dir_idx(i))
                                     lo_flux_R = lo_rho_R*lo_vel_R(dir_idx(i))*(s_R-lo_vel_R(dir_idx(1)))-lo_pres_R * &
                                         dir_flg(dir_idx(i))
                                     hi_flux_L = rho_L * vel_L(dir_idx(i))*(s_L - vel_L(dir_idx(1))) - &
                                         pres_L*dir_flg(dir_idx(i))
                                     hi_flux_R = rho_R * vel_R(dir_idx(i))*(s_R - vel_R(dir_idx(1))) - &
                                         pres_R*dir_flg(dir_idx(i))
                                END IF
                                tvd_cons_L = lo_cons_L + flux_lim_func*(hi_cons_L - lo_cons_L)
                                tvd_cons_R = lo_cons_R + flux_lim_func*(hi_cons_R - lo_cons_R)
                                tvd_flux_L = lo_flux_L + flux_lim_func*(hi_flux_L - lo_flux_L)
                                tvd_flux_R = lo_flux_R + flux_lim_func*(hi_flux_R - lo_flux_R)
                
                                IF (tvd_wave_speeds) THEN
                                    vel_src_rs_vf(dir_idx(i))%sf(j,k,l) = &
                                        (tvd_xi_M*tvd_flux_L - tvd_xi_P*tvd_flux_R) &
                                        / (tvd_xi_M*tvd_cons_L - tvd_xi_P*tvd_cons_R)
                                ELSE
                                    vel_src_rs_vf(dir_idx(i))%sf(j,k,l) = &
                                        (xi_M* tvd_flux_L - xi_P* tvd_flux_R) / (xi_M* tvd_cons_L - xi_P*tvd_cons_R)
                                END IF
                            END DO
                        ELSE
                            ! No TVD Riemann fluxes
                            
                            ! Mass
                            DO i = 1, cont_idx%end
                                flux_rs_vf(i)%sf(j,k,l) = &
                                      ( s_M*alpha_rho_R(i)*vel_R(dir_idx(1)) &
                                      - s_P*alpha_rho_L(i)*vel_L(dir_idx(1)) &
                                      + s_M*s_P*( alpha_rho_L(i)             &
                                                - alpha_rho_R(i) ) )         &
                                      / (s_M - s_P)
                            END DO
                            
                            ! Momentum
                            IF (bubbles) THEN
                                DO i = 1, num_dims
                                    flux_rs_vf(cont_idx%end+dir_idx(i))%sf(j,k,l) = &
                                            ( s_M*( rho_R*vel_R(dir_idx(1))         &
                                                         *vel_R(dir_idx(i))         &
                                                  + dir_flg(dir_idx(i))*(pres_R-ptilde_R) )    &
                                            - s_P*( rho_L*vel_L(dir_idx(1))         &
                                                         *vel_L(dir_idx(i))         &
                                                  + dir_flg(dir_idx(i))*(pres_L-ptilde_L) )    &
                                            + s_M*s_P*( rho_L*vel_L(dir_idx(i))     &
                                                      - rho_R*vel_R(dir_idx(i)) ) ) &
                                            / (s_M - s_P)
                                END DO
                            ELSE
                                DO i = 1, num_dims
                                    flux_rs_vf(cont_idx%end+dir_idx(i))%sf(j,k,l) = &
                                            ( s_M*( rho_R*vel_R(dir_idx(1))         &
                                                         *vel_R(dir_idx(i))         &
                                                  + dir_flg(dir_idx(i))*pres_R )    &
                                            - s_P*( rho_L*vel_L(dir_idx(1))         &
                                                         *vel_L(dir_idx(i))         &
                                                  + dir_flg(dir_idx(i))*pres_L )    &
                                            + s_M*s_P*( rho_L*vel_L(dir_idx(i))     &
                                                      - rho_R*vel_R(dir_idx(i)) ) ) &
                                            / (s_M - s_P)
                                END DO
                            END IF

                            ! Energy
                            IF (bubbles) THEN
                                flux_rs_vf(E_idx)%sf(j,k,l) = &
                                        ( s_M*vel_R(dir_idx(1))*(E_R + pres_R - ptilde_R) &
                                        - s_P*vel_L(dir_idx(1))*(E_L + pres_L - ptilde_L) &
                                        + s_M*s_P*(E_L - E_R) )                &
                                        / (s_M - s_P)
                            ELSE 
                                flux_rs_vf(E_idx)%sf(j,k,l) = &
                                        ( s_M*vel_R(dir_idx(1))*(E_R + pres_R) &
                                        - s_P*vel_L(dir_idx(1))*(E_L + pres_L) &
                                        + s_M*s_P*(E_L - E_R) )                &
                                        / (s_M - s_P)
                            END IF

                            ! Advection
                            DO i = adv_idx%beg, adv_idx%end
                                flux_rs_vf(i)%sf(j,k,l) = &
                                         ( qL_prim_rs_vf(i)%sf( j ,k,l)   &
                                         - qR_prim_rs_vf(i)%sf(j+1,k,l) ) &
                                         * s_M*s_P/(s_M - s_P)
                                flux_src_rs_vf(i)%sf(j,k,l) = &
                                       ( s_M*qR_prim_rs_vf(i)%sf(j+1,k,l)   &
                                       - s_P*qL_prim_rs_vf(i)%sf( j ,k,l) ) &
                                       / (s_M - s_P)
                            END DO
                            
                            ! Div(U)?
                            DO i = 1, num_dims
                                vel_src_rs_vf(dir_idx(i))%sf(j,k,l) = &
                                             ( xi_M*(rho_L*vel_L(dir_idx(i))  * &
                                                    (s_L - vel_L(dir_idx(1))) - &
                                                  pres_L*dir_flg(dir_idx(i))) - &
                                               xi_P*(rho_R*vel_R(dir_idx(i))  * &
                                                    (s_R - vel_R(dir_idx(1))) - &
                                                  pres_R*dir_flg(dir_idx(i))) ) &
                                     / ( xi_M*rho_L*(s_L - vel_L(dir_idx(1))) - &
                                         xi_P*rho_R*(s_R - vel_R(dir_idx(1))) )
                            END DO

                            ! SHB: Does this need to be ammended?
                            IF (hypoelasticity) THEN
                                DO i = 1, (num_dims*(num_dims+1)) / 2
                                    flux_rs_vf(stress_idx%beg-1+i)%sf(j,k,l) =  &
                                        ( s_M*( rho_R*vel_R(dir_idx(1))         &
                                                     *tau_e_R(i))               &
                                        - s_P*( rho_L*vel_L(dir_idx(1))         &
                                                     *tau_e_L(i))               &
                                        + s_M*s_P*( rho_L*tau_e_L(i)            &
                                                  - rho_R*tau_e_R(i) ) )        &
                                        / (s_M - s_P)
                                END DO 
                            END IF

                            IF (bubbles) THEN

                                ! Momentum: check signs on this
                                ! DO i = 1, num_dims
                                !     flux_rs_vf(cont_idx%end+dir_idx(i))%sf(j,k,l) = &
                                !     flux_rs_vf(cont_idx%end+dir_idx(i))%sf(j,k,l) - 
                                !             ( s_M*( rho_R*vel_R(dir_idx(1))         &
                                !                          *vel_R(dir_idx(i))         &
                                !                   + dir_flg(dir_idx(i))*ptilde_R )    &
                                !             - s_P*( rho_L*vel_L(dir_idx(1))         &
                                !                          *vel_L(dir_idx(i))         &
                                !                   + dir_flg(dir_idx(i))*ptilde_L )    &
                                !               ) &
                                !             / (s_M - s_P)
                                ! END DO
                                ! Energy: check signs on this
                                ! flux_rs_vf(E_idx)%sf(j,k,l) = &
                                ! flux_rs_vf(E_idx)%sf(j,k,l) - &
                                !         ( s_M*vel_R(dir_idx(1))*(E_R + ptilde_R) &
                                !         - s_P*vel_L(dir_idx(1))*(E_L + ptilde_L) &
                                !         )                                         &
                                !         / (s_M - s_P)


                                ! Advection: does this need to be changed? 
                                !   perhaps set bub_idx advection to zero? 
                                ! DO i = adv_idx%beg, adv_idx%end
                                !     flux_rs_vf(i)%sf(j,k,l) = &
                                !              ( qL_prim_rs_vf(i)%sf( j ,k,l)   &
                                !              - qR_prim_rs_vf(i)%sf(j+1,k,l) ) &
                                !              * s_M*s_P/(s_M - s_P)
                                !     flux_src_rs_vf(i)%sf(j,k,l) = &
                                !            ( s_M*qR_prim_rs_vf(i)%sf(j+1,k,l)   &
                                !            - s_P*qL_prim_rs_vf(i)%sf( j ,k,l) ) &
                                !            / (s_M - s_P)
                                ! END DO

                                ! From HLLC: Kills mass transport @ bubble gas density
                                IF ( num_fluids > 1 ) THEN
                                    flux_rs_vf(cont_idx%end)%sf(j,k,l) = 0d0
                                END IF

                                !Advection of bubble sources (from HLLC)
                                ! DO i = bub_idx%beg,sys_size
                                !     flux_rs_vf(i)%sf(j,k,l) =   &
                                !             xi_M*nbub_L*qL_prim_rs_vf(i)%sf(j,k,l)      &
                                !             * (vel_L(dir_idx(1)) + s_M*(xi_L - 1d0)) &
                                !             + xi_P*nbub_R*qR_prim_rs_vf(i)%sf(j+1,k,l)      &
                                !             * (vel_R(dir_idx(1)) + s_P*(xi_R - 1d0))
                                ! END DO

                            END IF

                        END IF
                    END DO
                END DO
            END DO
            
            
            ! Computing the viscous and capillary source flux
            IF(ANY(Re_size > 0) .OR. hypoelasticity) THEN
                IF (weno_Re_flux) THEN
                    CALL s_compute_viscous_source_flux( &
                                   qL_prim_vf(mom_idx%beg:mom_idx%end), &
                               dqL_prim_dx_vf(mom_idx%beg:mom_idx%end), &
                               dqL_prim_dy_vf(mom_idx%beg:mom_idx%end), &
                               dqL_prim_dz_vf(mom_idx%beg:mom_idx%end), &
                                   qR_prim_vf(mom_idx%beg:mom_idx%end), &
                               dqR_prim_dx_vf(mom_idx%beg:mom_idx%end), &
                               dqR_prim_dy_vf(mom_idx%beg:mom_idx%end), &
                               dqR_prim_dz_vf(mom_idx%beg:mom_idx%end), &
                                       flux_src_vf, norm_dir, ix,iy,iz  )
                ELSE
                    CALL s_compute_viscous_source_flux( &
                                    q_prim_vf(mom_idx%beg:mom_idx%end), &
                               dqL_prim_dx_vf(mom_idx%beg:mom_idx%end), &
                               dqL_prim_dy_vf(mom_idx%beg:mom_idx%end), &
                               dqL_prim_dz_vf(mom_idx%beg:mom_idx%end), &
                                    q_prim_vf(mom_idx%beg:mom_idx%end), &
                               dqR_prim_dx_vf(mom_idx%beg:mom_idx%end), &
                               dqR_prim_dy_vf(mom_idx%beg:mom_idx%end), &
                               dqR_prim_dz_vf(mom_idx%beg:mom_idx%end), &
                                       flux_src_vf, norm_dir, ix,iy,iz  )
                END IF
            END IF
            
            IF(We_size > 0 .AND. We_riemann_flux) THEN
                CALL s_compute_capillary_source_flux( &
                             dqL_prim_dx_vf(adv_idx%beg:adv_idx%end), &
                             dqL_prim_dy_vf(adv_idx%beg:adv_idx%end), &
                             dqL_prim_dz_vf(adv_idx%beg:adv_idx%end), &
                                                    gm_alphaL_vf , &
                             dqR_prim_dx_vf(adv_idx%beg:adv_idx%end), &
                             dqR_prim_dy_vf(adv_idx%beg:adv_idx%end), &
                             dqR_prim_dz_vf(adv_idx%beg:adv_idx%end), &
                                                    gm_alphaR_vf , &
                                 flux_src_vf, norm_dir, ix,iy,iz   )
            END IF
            
            
            ! Reshaping outputted data based on dimensional splitting direction
            CALL s_finalize_riemann_solver( flux_vf, flux_src_vf, &
                                                    flux_gsrc_vf, &
                                              norm_dir, ix,iy,iz  )
            
            
        END SUBROUTINE s_hll_riemann_solver ! ----------------------------------
        
        
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
        !!  @param kappaL_vf Left averaged vol. frac. curvatures
        !!  @param kappaR_vf Right averaged vol. rac. curvatures
        !!  @param flux_vf Intra-cell fluxes
        !!  @param flux_src_vf Intra-cell fluxes sources
        !!  @param flux_gsrc_vf Intra-cell geometric fluxes sources
        !!  @param norm_dir Dir. splitting direction
        !!  @param ix Index bounds in the x-dir
        !!  @param iy Index bounds in the y-dir
        !!  @param iz Index bounds in the z-dir
        !!  @param q_prim_vf Cell-averaged primitive variables 
        SUBROUTINE s_hllc_riemann_solver( qL_prim_vf, dqL_prim_dx_vf, & ! ------
                                                      dqL_prim_dy_vf, &
                                                      dqL_prim_dz_vf, &
                                                        gm_alphaL_vf, &
                                                           kappaL_vf, &
                                          qR_prim_vf, dqR_prim_dx_vf, &
                                                      dqR_prim_dy_vf, &
                                                      dqR_prim_dz_vf, &
                                                        gm_alphaR_vf, &
                                                           kappaR_vf, &
                                                           q_prim_vf, &
                                                flux_vf, flux_src_vf, &
                                                        flux_gsrc_vf, &
                                                  norm_dir, ix,iy,iz  )

            
            TYPE(scalar_field),  &
            DIMENSION(sys_size), &
            INTENT(INOUT) :: qL_prim_vf, qR_prim_vf
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q_prim_vf
            
            TYPE(scalar_field),        &
            ALLOCATABLE, DIMENSION(:), &
            INTENT(INOUT) :: dqL_prim_dx_vf, dqR_prim_dx_vf, &
                             dqL_prim_dy_vf, dqR_prim_dy_vf, &
                             dqL_prim_dz_vf, dqR_prim_dz_vf, &
                               gm_alphaL_vf,   gm_alphaR_vf, &
                                  kappaL_vf,      kappaR_vf
            
            ! Intercell fluxes
            TYPE(scalar_field),  &
            DIMENSION(sys_size), &
            INTENT(INOUT) :: flux_vf, flux_src_vf, flux_gsrc_vf
            
            INTEGER, INTENT(IN) :: norm_dir
            TYPE(bounds_info), INTENT(IN) :: ix,iy,iz
            

            REAL(KIND(0d0)) :: xi_L, xi_R !< Left and right wave speeds functions
            

            INTEGER :: i,j,k,l !< Generic loop iterators
            
            ! Placeholders for high- and low-order fluxes
            REAL(KIND(0d0)) :: lo_xi_L, lo_xi_R
            REAL(KIND(0d0)) :: hi_flux_L, hi_flux_R, lo_flux_L, lo_flux_R
            REAL(KIND(0d0)) :: tvd_flux_L, tvd_flux_R
            
            ! Populating the buffers of the left and right Riemann problem
            ! states variables, based on the choice of boundary conditions
            CALL s_populate_riemann_states_variables_buffers( &
                                  qL_prim_vf, dqL_prim_dx_vf, &
                                              dqL_prim_dy_vf, &
                                              dqL_prim_dz_vf, &
                                                gm_alphaL_vf, &
                                                   kappaL_vf, &
                                  qR_prim_vf, dqR_prim_dx_vf, &
                                              dqR_prim_dy_vf, &
                                              dqR_prim_dz_vf, &
                                                gm_alphaR_vf, &
                                                   kappaR_vf, &
                                          norm_dir, ix,iy,iz  )
            
            ! Reshaping inputted data based on dimensional splitting direction
            CALL s_initialize_riemann_solver( qL_prim_vf, kappaL_vf, &
                                              qR_prim_vf, kappaR_vf, &
                                                          q_prim_vf, &
                                               flux_vf, flux_src_vf, &
                                                       flux_gsrc_vf, &
                                                 norm_dir, ix,iy,iz  )

            ! Computing HLLC flux and source flux for Euler system of equations

            ! print*, 'xbounds are: ', is1%beg, is1%end
            ! print*, 'ybounds are: ', is2%beg, is2%end
            ! print*, 'zbounds are: ', is3%beg, is3%end
            DO l = is3%beg, is3%end
                DO k = is2%beg, is2%end
                    DO j = is1%beg, is1%end

                        ! print*, 'about to get average state'
                        CALL s_compute_average_state(j,k,l)
                        ! print*, 'got average state'

                        CALL s_compute_wave_speeds(j,k,l)
                        ! print*, 'got wave speeds'

                        IF(model_eqns == 3) THEN

                            ! Compute left solution state
                            IF (s_L >= 0d0) THEN
                                p_Star = pres_L ! Only usefull to recalculate the radial momentum geometric source flux
                                DO i = 1, num_fluids
                                    flux_rs_vf(i+adv_idx%beg-1)%sf(j,k,l) = &
                                        qL_prim_rs_vf(i+adv_idx%beg-1)%sf( j ,k,l)*s_S
                                    
                                    flux_rs_vf(i+cont_idx%beg-1)%sf(j,k,l) = &
                                        qL_prim_rs_vf(i+cont_idx%beg-1)%sf( j ,k,l)*vel_L(dir_idx(1))
                                    
                                    flux_rs_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l) = &
                                        qL_prim_rs_vf(i+adv_idx%beg-1)%sf( j ,k,l)*&
                                        (fluid_pp(i)%gamma*pres_L+fluid_pp(i)%pi_inf)*vel_L(dir_idx(1))
                                END DO
                                DO i = 1, num_dims
                                    flux_rs_vf(mom_idx%beg-1+dir_idx(i))%sf(j,k,l) = &
                                        rho_L*vel_L(dir_idx(1))*vel_L(dir_idx(i)) + dir_flg(dir_idx(i))*pres_L
                                    
                                    vel_src_rs_vf(dir_idx(i))%sf(j,k,l) = vel_L(dir_idx(i)) + &
                                        dir_flg(dir_idx(i))*(s_S - vel_L(dir_idx(i))) 
                                    ! Compute the star velocities for the non-conservative terms
                                END DO
                                flux_rs_vf(E_idx)%sf(j,k,l) = (E_L + pres_L)*vel_L(dir_idx(1))

                            ! Compute right solution state
                            ELSE IF (s_R <= 0d0) THEN
                                p_Star = pres_R 
                                ! Only usefull to recalculate the radial momentum geometric source flux
                                DO i = 1, num_fluids
                                    flux_rs_vf(i+adv_idx%beg-1)%sf(j,k,l) = &
                                        qR_prim_rs_vf(i+adv_idx%beg-1)%sf(j+1,k,l)*s_S
                                    
                                    flux_rs_vf(i+cont_idx%beg-1)%sf(j,k,l) = &
                                        qR_prim_rs_vf(i+cont_idx%beg-1)%sf(j+1,k,l)*vel_R(dir_idx(1))
                                    
                                    flux_rs_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l) = &
                                        qR_prim_rs_vf(i+adv_idx%beg-1)%sf(j+1,k,l) * &
                                        (fluid_pp(i)%gamma*pres_R+fluid_pp(i)%pi_inf)*vel_R(dir_idx(1))
                                END DO 
                                DO i = 1, num_dims
                                    flux_rs_vf(mom_idx%beg-1+dir_idx(i))%sf(j,k,l) = &
                                        rho_R*vel_R(dir_idx(1))*vel_R(dir_idx(i)) + dir_flg(dir_idx(i))*pres_R
                                    
                                    vel_src_rs_vf(dir_idx(i))%sf(j,k,l) = vel_R(dir_idx(i)) + &
                                        dir_flg(dir_idx(i))*(s_S - vel_R(dir_idx(i))) 
                                    ! Compute the star velocities for the non-conservative terms
                                END DO
                                flux_rs_vf(E_idx)%sf(j,k,l) = (E_R + pres_R)*vel_R(dir_idx(1))

                            ! Compute left star solution state
                            ELSE IF (s_S >= 0d0) THEN
                                xi_L = (s_L - vel_L(dir_idx(1)))/(s_L - s_S)
                                rho_Star = rho_L*xi_L
                                E_Star = xi_L*(E_L + (s_S - vel_L(dir_idx(1))) * &
                                    (rho_L*s_S + pres_L / (s_L - vel_L(dir_idx(1)))))
                                p_Star = rho_L*(s_L - vel_L(dir_idx(1)))*(s_S - vel_L(dir_idx(1))) + pres_L
                                DO i = 1, num_fluids
                                    p_K_Star = (pres_L+fluid_pp(i)%pi_inf/(1d0+fluid_pp(i)%gamma)) * &
                                        xi_L**(1d0/fluid_pp(i)%gamma+1d0) - fluid_pp(i)%pi_inf/(1d0+fluid_pp(i)%gamma)
                                    
                                    flux_rs_vf(i+adv_idx%beg-1)%sf(j,k,l) = &
                                        qL_prim_rs_vf(i+adv_idx%beg-1)%sf( j ,k,l)*s_S
                                    
                                    flux_rs_vf(i+cont_idx%beg-1)%sf(j,k,l) = &
                                        qL_prim_rs_vf(i+cont_idx%beg-1)%sf( j ,k,l)*xi_L*s_S
                                    
                                    flux_rs_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l) = &
                                        qL_prim_rs_vf(i+adv_idx%beg-1)%sf( j ,k,l) * &
                                        (fluid_pp(i)%gamma*p_K_Star+fluid_pp(i)%pi_inf)*s_S
                                END DO
                                DO i = 1, num_dims
                                    flux_rs_vf(mom_idx%beg-1+dir_idx(i))%sf(j,k,l) = &
                                        rho_Star*s_S*(s_S*dir_flg(dir_idx(i))+vel_L(dir_idx(i)) * &
                                        (1d0-dir_flg(dir_idx(i)))) + dir_flg(dir_idx(i))*p_Star

                                    vel_src_rs_vf(dir_idx(i))%sf(j,k,l) = vel_L(dir_idx(i)) + &
                                        dir_flg(dir_idx(i))*(s_S*xi_L - vel_L(dir_idx(i))) 
                                    ! Compute the star velocities for the non-conservative terms
                                END DO
                                flux_rs_vf(E_idx)%sf(j,k,l) = (E_Star + p_Star)*s_S

                            ! Compute right star solution state
                            ELSE
                                xi_R = (s_R - vel_R(dir_idx(1)))/(s_R - s_S)
                                
                                rho_Star = rho_R*xi_R
                                
                                E_Star = xi_R*(E_R + (s_S - vel_R(dir_idx(1))) * &
                                    (rho_R*s_S + pres_R / (s_R - vel_R(dir_idx(1)))))
                                
                                p_Star = rho_R*(s_R - vel_R(dir_idx(1)))*(s_S - vel_R(dir_idx(1))) + pres_R
                                
                                DO i = 1, num_fluids
                                    p_K_Star = (pres_R+fluid_pp(i)%pi_inf/(1d0+fluid_pp(i)%gamma)) * &
                                        xi_R**(1d0/fluid_pp(i)%gamma+1d0) - fluid_pp(i)%pi_inf/(1d0+fluid_pp(i)%gamma)
                                    
                                    flux_rs_vf(i+adv_idx%beg-1)%sf(j,k,l) = &
                                        qR_prim_rs_vf(i+adv_idx%beg-1)%sf(j+1,k,l)*s_S
                                    
                                    flux_rs_vf(i+cont_idx%beg-1)%sf(j,k,l) = &
                                        qR_prim_rs_vf(i+cont_idx%beg-1)%sf(j+1,k,l)*xi_R*s_S
                                    
                                    flux_rs_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l) = &
                                        qR_prim_rs_vf(i+adv_idx%beg-1)%sf(j+1,k,l) * & 
                                        (fluid_pp(i)%gamma*p_K_Star+fluid_pp(i)%pi_inf)*s_S
                                END DO
                               
                                DO i = 1, num_dims
                                    flux_rs_vf(mom_idx%beg-1+dir_idx(i))%sf(j,k,l) = rho_Star*s_S * &
                                        (s_S*dir_flg(dir_idx(i))+vel_R(dir_idx(i))*(1d0-dir_flg(dir_idx(i)))) + &
                                        dir_flg(dir_idx(i))*p_Star
                                    
                                    vel_src_rs_vf(dir_idx(i))%sf(j,k,l) = vel_R(dir_idx(i)) + &
                                        dir_flg(dir_idx(i))*(s_S*xi_R - vel_R(dir_idx(i))) 
                                    ! Compute the star velocities for the non-conservative terms
                                END DO

                                flux_rs_vf(E_idx)%sf(j,k,l) = (E_Star + p_Star)*s_S

                            END IF
                       
                            ! Geometrical source flux for cylindrical coordinates
                            IF(norm_dir == 2 .AND. cyl_coord) THEN
                                ! Substituting the advective flux into the inviscid geometrical source flux
                                DO i = 1, E_idx
                                    flux_gsrc_rs_vf(i)%sf(j,k,l) = flux_rs_vf(i)%sf(j,k,l)
                                END DO
                                DO i = internalEnergies_idx%beg, internalEnergies_idx%end
                                    flux_gsrc_rs_vf(i)%sf(j,k,l) = flux_rs_vf(i)%sf(j,k,l)
                                END DO
                                ! Recalculating the radial momentum geometric source flux (substracting the pressure part)
                                flux_gsrc_rs_vf(mom_idx%beg-1+dir_idx(1))%sf(j,k,l) = &
                                    flux_gsrc_rs_vf(mom_idx%beg-1+dir_idx(1))%sf(j,k,l) - p_Star
                                ! Geometrical source of the void fraction(s) is zero
                                DO i = adv_idx%beg, adv_idx%end
                                    flux_gsrc_rs_vf(i)%sf(j,k,l) = 0d0
                                END DO
                            ELSE
                                DO i = 1, sys_size
                                    flux_gsrc_rs_vf(i)%sf(j,k,l) = 0d0
                                END DO
                            END IF

                        ELSE

                            IF (tvd_riemann_flux) THEN
                                CALL s_compute_flux_limiter(j,k,l,flux_lim_func,norm_dir)
                            END IF

                            ! follows Einfeldt et al.
                            ! s_M/P = min/max(0.,s_L/R)
                            s_M = MIN(0d0,s_L); s_P = MAX(0d0,s_R)
                            
                            ! goes with q_star_L/R = xi_L/R * (variable)
                            ! xi_L/R = ( ( s_L/R - u_L/R )/(s_L/R - s_star) )
                            xi_L = (s_L - vel_L(dir_idx(1)))/(s_L - s_S)
                            xi_R = (s_R - vel_R(dir_idx(1)))/(s_R - s_S)
                            
                            ! goes with numerical velocity in x/y/z directions
                            ! xi_P/M = 0.5 +/m sgn(0.5,s_star) 
                            xi_M = (5d-1 + SIGN(5d-1,s_S))
                            xi_P = (5d-1 - SIGN(5d-1,s_S))


                            IF (tvd_wave_speeds) THEN
                                lo_s_M = MIN(0d0,lo_s_L); lo_s_P = MAX(0d0,lo_s_R)

                                lo_xi_L = (lo_s_L - lo_vel_L(dir_idx(1)))/(lo_s_L - lo_s_S)
                                lo_xi_R = (lo_s_R - lo_vel_R(dir_idx(1)))/(lo_s_R - lo_s_S)

                                lo_xi_M = (5d-1 + SIGN(5d-1,lo_s_S))
                                lo_xi_P = (5d-1 - SIGN(5d-1,lo_s_S))

                                tvd_xi_M = lo_xi_M + flux_lim_func*(xi_M - lo_xi_M)
                                tvd_xi_P = lo_xi_P + flux_lim_func*(xi_P - lo_xi_P)
                            END IF

                            IF (tvd_riemann_flux) THEN
                                ! Mass
                                DO i = 1, cont_idx%end
                                    IF (tvd_wave_speeds) THEN
                                         lo_flux_L = lo_alpha_rho_L(i)*(lo_vel_L(dir_idx(1)) + lo_s_M*(lo_xi_L - 1d0))
                                         lo_flux_R = lo_alpha_rho_R(i)*(lo_vel_R(dir_idx(1)) + lo_s_P*(lo_xi_R - 1d0))
                                         hi_flux_L =    alpha_rho_L(i)*(   vel_L(dir_idx(1)) +    s_M*(   xi_L - 1d0))
                                         hi_flux_R =    alpha_rho_R(i)*(   vel_R(dir_idx(1)) +    s_P*(   xi_R - 1d0))
                                    ELSE
                                         lo_flux_L = lo_alpha_rho_L(i)*(lo_vel_L(dir_idx(1)) + s_M*(xi_L - 1d0))
                                         lo_flux_R = lo_alpha_rho_R(i)*(lo_vel_R(dir_idx(1)) + s_P*(xi_R - 1d0))
                                         hi_flux_L =    alpha_rho_L(i)*(   vel_L(dir_idx(1)) + s_M*(xi_L - 1d0))
                                         hi_flux_R =    alpha_rho_R(i)*(   vel_R(dir_idx(1)) + s_P*(xi_R - 1d0))
                                    END IF
                                    tvd_flux_L = lo_flux_L + flux_lim_func*(hi_flux_L - lo_flux_L)
                                    tvd_flux_R = lo_flux_R + flux_lim_func*(hi_flux_R - lo_flux_R)
                    
                                    IF (tvd_wave_speeds) THEN
                                        flux_rs_vf(i)%sf(j,k,l) = tvd_xi_M*tvd_flux_L + tvd_xi_P*tvd_flux_R
                                    ELSE
                                        flux_rs_vf(i)%sf(j,k,l) = xi_M* tvd_flux_L + xi_P* tvd_flux_R
                                    END IF
                                END DO
                    
                                ! Momentum
                                DO i = 1, num_dims
                                    IF (tvd_wave_speeds) THEN
                                         lo_flux_L = lo_rho_L*(lo_vel_L(dir_idx(1))*lo_vel_L(dir_idx(i))+       &
                                                lo_s_M*(lo_xi_L*(dir_flg(dir_idx(i))*lo_s_S+(1d0-dir_flg(dir_idx(i)))*  &
                                                lo_vel_L(dir_idx(i)))-lo_vel_L(dir_idx(i))))+dir_flg(dir_idx(i))*lo_pres_L
                                         lo_flux_R = lo_rho_R*(lo_vel_R(dir_idx(1))*lo_vel_R(dir_idx(i))+       &
                                                lo_s_P*(lo_xi_R*(dir_flg(dir_idx(i))*lo_s_S+(1d0-dir_flg(dir_idx(i)))*  &
                                                lo_vel_R(dir_idx(i)))-lo_vel_R(dir_idx(i))))+dir_flg(dir_idx(i))*lo_pres_R
                                         hi_flux_L =    rho_L*(   vel_L(dir_idx(1))*   vel_L(dir_idx(i))+       &
                                                s_M*(   xi_L*(dir_flg(dir_idx(i))*s_S+(1d0-dir_flg(dir_idx(i)))*    &
                                                   vel_L(dir_idx(i)))-   vel_L(dir_idx(i))))+dir_flg(dir_idx(i))*   pres_L
                                         hi_flux_R =    rho_R*(   vel_R(dir_idx(1))*   vel_R(dir_idx(i))+       &
                                                s_P*(   xi_R*(dir_flg(dir_idx(i))*s_S+(1d0-dir_flg(dir_idx(i)))*    &
                                                   vel_R(dir_idx(i)))-   vel_R(dir_idx(i))))+dir_flg(dir_idx(i))*   pres_R
                                    ELSE
                                         lo_flux_L = lo_rho_L*(lo_vel_L(dir_idx(1))*lo_vel_L(dir_idx(i))+       &
                                                s_M*(xi_L*(dir_flg(dir_idx(i))*s_S+(1d0-dir_flg(dir_idx(i)))*   &
                                                lo_vel_L(dir_idx(i)))-lo_vel_L(dir_idx(i))))+dir_flg(dir_idx(i))*lo_pres_L
                                         lo_flux_R = lo_rho_R*(lo_vel_R(dir_idx(1))*lo_vel_R(dir_idx(i))+       &
                                                s_P*(xi_R*(dir_flg(dir_idx(i))*s_S+(1d0-dir_flg(dir_idx(i)))*   &
                                                lo_vel_R(dir_idx(i)))-lo_vel_R(dir_idx(i))))+dir_flg(dir_idx(i))*lo_pres_R
                                         hi_flux_L =    rho_L*(   vel_L(dir_idx(1))*   vel_L(dir_idx(i))+       &
                                                s_M*(xi_L*(dir_flg(dir_idx(i))*s_S+(1d0-dir_flg(dir_idx(i)))*   &
                                                   vel_L(dir_idx(i)))-   vel_L(dir_idx(i))))+dir_flg(dir_idx(i))*   pres_L
                                         hi_flux_R =    rho_R*(   vel_R(dir_idx(1))*   vel_R(dir_idx(i))+       &
                                                s_P*(xi_R*(dir_flg(dir_idx(i))*s_S+(1d0-dir_flg(dir_idx(i)))*   &
                                                   vel_R(dir_idx(i)))-   vel_R(dir_idx(i))))+dir_flg(dir_idx(i))*   pres_R
                                    END IF
                                    tvd_flux_L = lo_flux_L + flux_lim_func*(hi_flux_L - lo_flux_L)
                                    tvd_flux_R = lo_flux_R + flux_lim_func*(hi_flux_R - lo_flux_R)
                    
                                    IF (tvd_wave_speeds) THEN
                                        flux_rs_vf(cont_idx%end+dir_idx(i))%sf(j,k,l) = tvd_xi_M*tvd_flux_L + tvd_xi_P*tvd_flux_R
                                    ELSE
                                        flux_rs_vf(cont_idx%end+dir_idx(i))%sf(j,k,l) = xi_M* tvd_flux_L + xi_P* tvd_flux_R
                                    END IF
                                END DO
                    
                                ! Energy
                    
                                IF (tvd_wave_speeds) THEN
                                     lo_flux_L = lo_vel_L(dir_idx(1))*(lo_E_L+lo_pres_L) + &
                                            lo_s_M*(lo_xi_L*(lo_E_L+(lo_s_S-       &
                                            lo_vel_L(dir_idx(1)))*(lo_rho_L*lo_s_S + &
                                            lo_pres_L/(lo_s_L-lo_vel_L(dir_idx(1)))))-lo_E_L)
                                     lo_flux_R = lo_vel_R(dir_idx(1)) * &
                                            (lo_E_R+lo_pres_R)+lo_s_P*(lo_xi_R*(lo_E_R+(lo_s_S-  &
                                            lo_vel_R(dir_idx(1)))*(lo_rho_R*lo_s_S + &
                                            lo_pres_R/(lo_s_R-lo_vel_R(dir_idx(1)))))-lo_E_R)
                                     hi_flux_L =    vel_L(dir_idx(1))*(   E_L+   pres_L)+s_M*(   xi_L*(   E_L+(s_S-  &
                                               vel_L(dir_idx(1)))*&
                                               (   rho_L*s_S+   pres_L/(s_L-   vel_L(dir_idx(1)))))-   E_L)
                                     hi_flux_R =    vel_R(dir_idx(1))*(   E_R+   pres_R)+ &
                                         s_P*(   xi_R*(   E_R+(s_S-     &
                                               vel_R(dir_idx(1)))*(   rho_R*s_S+   &
                                               pres_R/(s_R-   vel_R(dir_idx(1)))))-   E_R)
                                ELSE
                                     lo_flux_L = lo_vel_L(dir_idx(1))*(lo_E_L+lo_pres_L)+s_M*(xi_L*(lo_E_L+(s_S-        &
                                            lo_vel_L(dir_idx(1)))*(lo_rho_L*s_S+lo_pres_L / &
                                            (s_L-lo_vel_L(dir_idx(1)))))-lo_E_L)
                                     lo_flux_R = lo_vel_R(dir_idx(1))*(lo_E_R+lo_pres_R)+s_P*(xi_R*(lo_E_R+(s_S-        &
                                            lo_vel_R(dir_idx(1)))*(lo_rho_R*s_S+lo_pres_R / &
                                            (s_R-lo_vel_R(dir_idx(1)))))-lo_E_R)
                                     hi_flux_L =    vel_L(dir_idx(1))*(   E_L+   pres_L)+s_M*(xi_L*(   E_L+(s_S-        &
                                               vel_L(dir_idx(1)))*(   rho_L*s_S+   pres_L / &
                                               (s_L-   vel_L(dir_idx(1)))))-   E_L)
                                     hi_flux_R =    vel_R(dir_idx(1))*(   E_R+   pres_R)+s_P*(xi_R*(   E_R+(s_S-        &
                                               vel_R(dir_idx(1)))*(   rho_R*s_S+   pres_R / &
                                               (s_R-   vel_R(dir_idx(1)))))-   E_R)
                                END IF
                                tvd_flux_L = lo_flux_L + flux_lim_func*(hi_flux_L - lo_flux_L)
                                tvd_flux_R = lo_flux_R + flux_lim_func*(hi_flux_R - lo_flux_R)
                    
                                IF (tvd_wave_speeds) THEN
                                    flux_rs_vf(E_idx)%sf(j,k,l) = tvd_xi_M*tvd_flux_L + tvd_xi_P*tvd_flux_R
                                ELSE
                                    flux_rs_vf(E_idx)%sf(j,k,l) = xi_M* tvd_flux_L + xi_P* tvd_flux_R
                                END IF
                    
                                ! Volume fraction
                                DO i = adv_idx%beg, adv_idx%end
                                    IF (tvd_wave_speeds) THEN
                                         lo_flux_L =  q_prim_rs_vf(i)%sf( j ,k,l)*(lo_vel_L(dir_idx(1))+lo_s_M*(lo_xi_L-1d0))
                                         lo_flux_R =  q_prim_rs_vf(i)%sf(j+1,k,l)*(lo_vel_R(dir_idx(1))+lo_s_P*(lo_xi_R-1d0))
                                         hi_flux_L = qL_prim_rs_vf(i)%sf( j ,k,l)*(   vel_L(dir_idx(1))+   s_M*(   xi_L-1d0))
                                         hi_flux_R = qR_prim_rs_vf(i)%sf(j+1,k,l)*(   vel_R(dir_idx(1))+   s_P*(   xi_R-1d0))
                                    ELSE
                                         lo_flux_L =  q_prim_rs_vf(i)%sf( j ,k,l)*(lo_vel_L(dir_idx(1))+s_M*(xi_L-1d0))
                                         lo_flux_R =  q_prim_rs_vf(i)%sf(j+1,k,l)*(lo_vel_R(dir_idx(1))+s_P*(xi_R-1d0))
                                         hi_flux_L = qL_prim_rs_vf(i)%sf( j ,k,l)*(   vel_L(dir_idx(1))+s_M*(xi_L-1d0))
                                         hi_flux_R = qR_prim_rs_vf(i)%sf(j+1,k,l)*(   vel_R(dir_idx(1))+s_P*(xi_R-1d0))
                                    END IF
                                    tvd_flux_L = lo_flux_L + flux_lim_func*(hi_flux_L - lo_flux_L)
                                    tvd_flux_R = lo_flux_R + flux_lim_func*(hi_flux_R - lo_flux_R)
                    
                                    IF (tvd_wave_speeds) THEN
                                        flux_rs_vf(i)%sf(j,k,l) = tvd_xi_M*tvd_flux_L + tvd_xi_P*tvd_flux_R
                                    ELSE
                                        flux_rs_vf(i)%sf(j,k,l) = xi_M* tvd_flux_L + xi_P* tvd_flux_R
                                    END IF
                                END DO
                    
                                ! Source terms
                                DO i = 1, num_dims
                                    IF (tvd_wave_speeds) THEN
                                         lo_flux_L = lo_vel_L(dir_idx(i)) + dir_flg(dir_idx(i))*lo_s_M*(lo_xi_L-1d0)
                                         lo_flux_R = lo_vel_R(dir_idx(i)) + dir_flg(dir_idx(i))*lo_s_P*(lo_xi_R-1d0)
                                         hi_flux_L =    vel_L(dir_idx(i)) + dir_flg(dir_idx(i))*   s_M*(   xi_L-1d0)
                                         hi_flux_R =    vel_R(dir_idx(i)) + dir_flg(dir_idx(i))*   s_P*(   xi_R-1d0)
                                    ELSE
                                         lo_flux_L = lo_vel_L(dir_idx(i)) + dir_flg(dir_idx(i))*s_M*(xi_L-1d0)
                                         lo_flux_R = lo_vel_R(dir_idx(i)) + dir_flg(dir_idx(i))*s_P*(xi_R-1d0)
                                         hi_flux_L =    vel_L(dir_idx(i)) + dir_flg(dir_idx(i))*s_M*(xi_L-1d0)
                                         hi_flux_R =    vel_R(dir_idx(i)) + dir_flg(dir_idx(i))*s_P*(xi_R-1d0)
                                    END IF
                                    tvd_flux_L = lo_flux_L + flux_lim_func*(hi_flux_L - lo_flux_L)
                                    tvd_flux_R = lo_flux_R + flux_lim_func*(hi_flux_R - lo_flux_R)
                    
                                    IF (tvd_wave_speeds) THEN
                                        vel_src_rs_vf(dir_idx(i))%sf(j,k,l) = tvd_xi_M*tvd_flux_L + tvd_xi_P*tvd_flux_R
                                    ELSE
                                        vel_src_rs_vf(dir_idx(i))%sf(j,k,l) = xi_M* tvd_flux_L + xi_P* tvd_flux_R
                                    END IF
                                END DO
                            ELSE

                                DO i = 1, cont_idx%end
                                    flux_rs_vf(i)%sf(j,k,l) = &
                                           xi_M*alpha_rho_L(i)                    &
                                         * (vel_L(dir_idx(1)) + s_M*(xi_L - 1d0)) &
                                         + xi_P*alpha_rho_R(i)                    &
                                         * (vel_R(dir_idx(1)) + s_P*(xi_R - 1d0))
                                END DO

                                IF (bubbles .AND. (model_eqns == 2) .AND. (num_fluids > 1) ) THEN
                                    ! Kill mass transport @ gas density
                                    flux_rs_vf(cont_idx%end)%sf(j,k,l) = 0.d0
                                END IF

                                ! Momentum flux. 
                                ! f = \rho u u + p I, q = \rho u, q_star = \xi * \rho*(s_star, v, w)
                                IF (bubbles .NEQV. .TRUE.) THEN
                                    DO i = 1, num_dims
                                        flux_rs_vf(cont_idx%end+dir_idx(i))%sf(j,k,l) = &
                                            xi_M*(rho_L*(vel_L(dir_idx(1)) *            &
                                            vel_L(dir_idx(i)) +                         &
                                            s_M*(xi_L*(dir_flg(dir_idx(i))*s_S +        &
                                            (1d0 - dir_flg(dir_idx(i))) *               &
                                            vel_L(dir_idx(i))) - vel_L(dir_idx(i)))) +  &
                                            dir_flg(dir_idx(i))*(pres_L))      &
                                        + xi_P*(rho_R*(vel_R(dir_idx(1)) *              &
                                            vel_R(dir_idx(i)) +                         &
                                            s_P*(xi_R*(dir_flg(dir_idx(i))*s_S +        &
                                            (1d0 - dir_flg(dir_idx(i))) *               &
                                            vel_R(dir_idx(i))) - vel_R(dir_idx(i)))) +  &
                                            dir_flg(dir_idx(i))*(pres_R))
                                        ! if (j==0) print*, 'flux_rs_vf', flux_rs_vf(cont_idx%end+dir_idx(i))%sf(j,k,l)
                                    END DO
                                ELSE
                                    ! Include p_tilde
                                    DO i = 1, num_dims
                                        flux_rs_vf(cont_idx%end+dir_idx(i))%sf(j,k,l) = &
                                            xi_M*(rho_L*(vel_L(dir_idx(1)) *            &
                                            vel_L(dir_idx(i)) +                         &
                                            s_M*(xi_L*(dir_flg(dir_idx(i))*s_S +        &
                                            (1d0 - dir_flg(dir_idx(i))) *               &
                                            vel_L(dir_idx(i))) - vel_L(dir_idx(i)))) +  &
                                            dir_flg(dir_idx(i))*(pres_L-ptilde_L))      &
                                        + xi_P*(rho_R*(vel_R(dir_idx(1)) *              &
                                            vel_R(dir_idx(i)) +                         &
                                            s_P*(xi_R*(dir_flg(dir_idx(i))*s_S +        &
                                            (1d0 - dir_flg(dir_idx(i))) *               &
                                            vel_R(dir_idx(i))) - vel_R(dir_idx(i)))) +  &
                                            dir_flg(dir_idx(i))*(pres_R-ptilde_R))
                                        ! if (j==0) print*, 'flux_rs_vf', flux_rs_vf(cont_idx%end+dir_idx(i))%sf(j,k,l)
                                    END DO

                                END IF
                                ! print*, xi_P, xi_M
                                ! print*, s_M, s_P, s_S
                                ! print*, ptilde_R, ptilde_L
                                ! print*, rho_R, rho_L
                                ! print*, vel_R(1), vel_L(1)
                                ! print*, pres_R, pres_L
                                ! print*, '              '

                                ! Energy flux.
                                ! f = u*(E+p), q = E, q_star = \xi*E+(s-u)(\rho s_star + p/(s-u))
                                IF ( (model_eqns .NE. 4) .AND. bubbles ) THEN
                                    flux_rs_vf(E_idx)%sf(j,k,l) = &
                                     xi_M*(vel_L(dir_idx(1))*(E_L + pres_L-ptilde_L) +     &
                                        s_M*(xi_L*(E_L + (s_S - vel_L(dir_idx(1))) * &
                                        (rho_L*s_S + (pres_L - ptilde_L) /                        &
                                        (s_L - vel_L(dir_idx(1))))) - E_L))          &
                                   + xi_P*(vel_R(dir_idx(1))*(E_R + pres_R-ptilde_R) +     &
                                        s_P*(xi_R*(E_R + (s_S - vel_R(dir_idx(1))) * &
                                        (rho_R*s_S + (pres_R - ptilde_R) /                        &
                                        (s_R - vel_R(dir_idx(1))))) - E_R))
                                ELSE IF ( model_eqns .NE. 4 ) THEN
                                    flux_rs_vf(E_idx)%sf(j,k,l) = &
                                     xi_M*(vel_L(dir_idx(1))*(E_L + pres_L) +     &
                                        s_M*(xi_L*(E_L + (s_S - vel_L(dir_idx(1))) * &
                                        (rho_L*s_S + pres_L /                        &
                                        (s_L - vel_L(dir_idx(1))))) - E_L))          &
                                   + xi_P*(vel_R(dir_idx(1))*(E_R + pres_R) +     &
                                        s_P*(xi_R*(E_R + (s_S - vel_R(dir_idx(1))) * &
                                        (rho_R*s_S + pres_R /                        &
                                        (s_R - vel_R(dir_idx(1))))) - E_R))
                                ELSE
                                    ! Tait EOS, no energy equation
                                    flux_rs_vf(E_idx)%sf(j,k,l) = 0.d0
                                END IF


                                ! Volume fraction flux
                                IF (model_eqns == 4) THEN
                                    DO i = alf_idx, alf_idx !only advect the void fraction
                                        flux_rs_vf(i)%sf(j,k,l) =  &
                                                  xi_M*qL_prim_rs_vf(i)%sf(j,k,l)      &
                                                * (vel_L(dir_idx(1)) + s_M*(xi_L - 1d0)) &
                                                + xi_P*qR_prim_rs_vf(i)%sf(j+1,k,l)      &
                                                * (vel_R(dir_idx(1)) + s_P*(xi_R - 1d0))
                                    END DO
                                ELSE
                                    DO i = adv_idx%beg, adv_idx%end 
                                        flux_rs_vf(i)%sf(j,k,l) =  &
                                                  xi_M*qL_prim_rs_vf(i)%sf(j,k,l)      &
                                                * (vel_L(dir_idx(1)) + s_M*(xi_L - 1d0)) &
                                                + xi_P*qR_prim_rs_vf(i)%sf(j+1,k,l)      &
                                                * (vel_R(dir_idx(1)) + s_P*(xi_R - 1d0))
                                    END DO
                                END IF

                                ! Source for volume fraction advection equation
                                DO i = 1, num_dims
                                    vel_src_rs_vf(dir_idx(i))%sf(j,k,l) = &
                                        xi_M*(   vel_L(dir_idx(i)) + &
                                            dir_flg(dir_idx(i)) * &
                                            s_M*(xi_L - 1d0) ) &
                                      + xi_P*(   vel_R(dir_idx(i)) + &
                                            dir_flg(dir_idx(i)) * &
                                            s_P*(xi_R - 1d0) ) 
                                    
                                    IF ( model_eqns == 4 ) vel_src_rs_vf(dir_idx(i))%sf(j,k,l) = 0d0
                                    !IF ( (model_eqns == 4) .or. (num_fluids==1) ) vel_src_rs_vf(dir_idx(i))%sf(j,k,l) = 0d0
                                END DO


                                ! Add advection flux for bubble variables
                                IF (bubbles) THEN
                                    DO i = bub_idx%beg,bub_idx%end
                                        flux_rs_vf(i)%sf(j,k,l) =   &
                                                xi_M*nbub_L*qL_prim_rs_vf(i)%sf(j,k,l)      &
                                                * (vel_L(dir_idx(1)) + s_M*(xi_L - 1d0)) &
                                                + xi_P*nbub_R*qR_prim_rs_vf(i)%sf(j+1,k,l)      &
                                                * (vel_R(dir_idx(1)) + s_P*(xi_R - 1d0))
                                    END DO
                                END IF

                            
                                ! print*, 'mom idx', cont_idx%end+dir_idx(i)
                                ! DO i = 1,sys_size
                                !     print*, 'flux_rs_vf', i, flux_rs_vf(i)%sf(j,k,l)
                                ! END DO

                                ! Geometrical source flux for cylindrical coordinates
                                IF(norm_dir == 2 .AND. cyl_coord) THEN
                                    ! Substituting the advective flux into the inviscid geometrical source flux
                                    DO i = 1, E_idx
                                        flux_gsrc_rs_vf(i)%sf(j,k,l) = flux_rs_vf(i)%sf(j,k,l)
                                    END DO
                                    ! Recalculating the radial momentum geometric source flux
                                    flux_gsrc_rs_vf(cont_idx%end+dir_idx(1))%sf(j,k,l) = &
                                              xi_M*(rho_L*(vel_L(dir_idx(1)) *           &
                                              vel_L(dir_idx(1)) +                        &
                                              s_M*(xi_L*(dir_flg(dir_idx(1))*s_S +       &
                                              (1d0 - dir_flg(dir_idx(1))) *              &
                                              vel_L(dir_idx(1))) - vel_L(dir_idx(1)))))  &
                                            + xi_P*(rho_R*(vel_R(dir_idx(1)) *           &
                                              vel_R(dir_idx(1)) +                        &
                                              s_P*(xi_R*(dir_flg(dir_idx(1))*s_S +       &
                                              (1d0 - dir_flg(dir_idx(1))) *              &
                                              vel_R(dir_idx(1))) - vel_R(dir_idx(1)))))
                                    ! Geometrical source of the void fraction(s) is zero
                                    DO i = adv_idx%beg, adv_idx%end
                                        flux_gsrc_rs_vf(i)%sf(j,k,l) = 0d0
                                    END DO
                                ELSEIF (norm_dir == 3 .AND. grid_geometry == 3) THEN
                                    DO i = 1, sys_size
                                        flux_gsrc_rs_vf(i)%sf(j,k,l) = 0d0
                                    END DO
                                    flux_gsrc_rs_vf(mom_idx%beg+1)%sf(j,k,l) = &
                                            - xi_M*(rho_L*(vel_L(dir_idx(1)) *           &
                                              vel_L(dir_idx(1)) +                        &
                                              s_M*(xi_L*(dir_flg(dir_idx(1))*s_S +       &
                                              (1d0 - dir_flg(dir_idx(1))) *              &
                                              vel_L(dir_idx(1))) - vel_L(dir_idx(1)))))  &
                                            - xi_P*(rho_R*(vel_R(dir_idx(1)) *           &
                                              vel_R(dir_idx(1)) +                        &
                                              s_P*(xi_R*(dir_flg(dir_idx(1))*s_S +       &
                                              (1d0 - dir_flg(dir_idx(1))) *              &
                                              vel_R(dir_idx(1))) - vel_R(dir_idx(1)))))
                                    flux_gsrc_rs_vf(mom_idx%end)%sf(j,k,l) = flux_rs_vf(mom_idx%beg+1)%sf(j,k,l)
                                ELSE
                                    DO i = 1, sys_size
                                        flux_gsrc_rs_vf(i)%sf(j,k,l) = 0d0
                                    END DO
                                END IF

                            END IF
                        END IF
                    END DO
                END DO
            END DO
            
            
            ! Computing the viscous and capillary source flux
            IF(ANY(Re_size > 0)) THEN
                IF (weno_Re_flux) THEN
                    CALL s_compute_viscous_source_flux( &
                                   qL_prim_vf(mom_idx%beg:mom_idx%end), &
                               dqL_prim_dx_vf(mom_idx%beg:mom_idx%end), &
                               dqL_prim_dy_vf(mom_idx%beg:mom_idx%end), &
                               dqL_prim_dz_vf(mom_idx%beg:mom_idx%end), &
                                   qR_prim_vf(mom_idx%beg:mom_idx%end), &
                               dqR_prim_dx_vf(mom_idx%beg:mom_idx%end), &
                               dqR_prim_dy_vf(mom_idx%beg:mom_idx%end), &
                               dqR_prim_dz_vf(mom_idx%beg:mom_idx%end), &
                                       flux_src_vf, norm_dir, ix,iy,iz  )
                ELSE
                    CALL s_compute_viscous_source_flux( &
                                    q_prim_vf(mom_idx%beg:mom_idx%end), &
                               dqL_prim_dx_vf(mom_idx%beg:mom_idx%end), &
                               dqL_prim_dy_vf(mom_idx%beg:mom_idx%end), &
                               dqL_prim_dz_vf(mom_idx%beg:mom_idx%end), &
                                    q_prim_vf(mom_idx%beg:mom_idx%end), &
                               dqR_prim_dx_vf(mom_idx%beg:mom_idx%end), &
                               dqR_prim_dy_vf(mom_idx%beg:mom_idx%end), &
                               dqR_prim_dz_vf(mom_idx%beg:mom_idx%end), &
                                       flux_src_vf, norm_dir, ix,iy,iz  )
                END IF
            END IF
            
            IF(We_size > 0 .AND. We_riemann_flux) THEN
                CALL s_compute_capillary_source_flux( &
                             dqL_prim_dx_vf(adv_idx%beg:adv_idx%end), &
                             dqL_prim_dy_vf(adv_idx%beg:adv_idx%end), &
                             dqL_prim_dz_vf(adv_idx%beg:adv_idx%end), &
                                                    gm_alphaL_vf , &
                             dqR_prim_dx_vf(adv_idx%beg:adv_idx%end), &
                             dqR_prim_dy_vf(adv_idx%beg:adv_idx%end), &
                             dqR_prim_dz_vf(adv_idx%beg:adv_idx%end), &
                                                    gm_alphaR_vf , &
                                 flux_src_vf, norm_dir, ix,iy,iz   )
            END IF
            
            
            ! Reshaping outputted data based on dimensional splitting direction
            CALL s_finalize_riemann_solver( flux_vf, flux_src_vf, &
                                                    flux_gsrc_vf, &
                                              norm_dir, ix,iy,iz  )
            
            
        END SUBROUTINE s_hllc_riemann_solver ! ---------------------------------
        
        
        
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
        !!  @param kappaL_vf  Left averaged vol. frac. curvatures
        !!  @param kappaR_vf Right averaged vol. frac. curvatures
        !!  @param q_prim_vf Cell-averaged primitive variables
        !!  @param flux_vf Intra-cell fluxes
        !!  @param flux_src_vf  Intra-cell fluxes sources
        !!  @param flux_gsrc_vf Intra-cell geometric fluxes sources
        !!  @param norm_dir Dir. splitting direction
        !!  @param ix Index bounds in the x-dir
        !!  @param iy Index bounds in the y-dir
        !!  @param iz Index bounds in the z-dir 
        SUBROUTINE s_exact_riemann_solver( qL_prim_vf, dqL_prim_dx_vf, & ! -----
                                                       dqL_prim_dy_vf, &
                                                       dqL_prim_dz_vf, &
                                                         gm_alphaL_vf, &
                                                            kappaL_vf, &
                                           qR_prim_vf, dqR_prim_dx_vf, &
                                                       dqR_prim_dy_vf, &
                                                       dqR_prim_dz_vf, &
                                                         gm_alphaR_vf, &
                                                            kappaR_vf, &
                                                            q_prim_vf, &
                                                 flux_vf, flux_src_vf, &
                                                         flux_gsrc_vf, &
                                                   norm_dir, ix,iy,iz  )
            
            TYPE(scalar_field),  &
            DIMENSION(sys_size), &
            INTENT(INOUT) :: qL_prim_vf, qR_prim_vf
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q_prim_vf
            
            TYPE(scalar_field),        &
            ALLOCATABLE, DIMENSION(:), &
            INTENT(INOUT) :: dqL_prim_dx_vf, dqR_prim_dx_vf, &
                             dqL_prim_dy_vf, dqR_prim_dy_vf, &
                             dqL_prim_dz_vf, dqR_prim_dz_vf, &
                               gm_alphaL_vf,   gm_alphaR_vf, &
                                  kappaL_vf,      kappaR_vf
            
            TYPE(scalar_field),  &
            DIMENSION(sys_size), &
            INTENT(INOUT) :: flux_vf, flux_src_vf, flux_gsrc_vf
            
            INTEGER, INTENT(IN) :: norm_dir
            TYPE(bounds_info), INTENT(IN) :: ix,iy,iz
            

            INTEGER :: i,j,k,l !< Generic loop iterators
            
            
            ! Populating the buffers of the left and right Riemann problem
            ! states variables, based on the choice of boundary conditions
            CALL s_populate_riemann_states_variables_buffers( &
                                  qL_prim_vf, dqL_prim_dx_vf, &
                                              dqL_prim_dy_vf, &
                                              dqL_prim_dz_vf, &
                                                gm_alphaL_vf, &
                                                   kappaL_vf, &
                                  qR_prim_vf, dqR_prim_dx_vf, &
                                              dqR_prim_dy_vf, &
                                              dqR_prim_dz_vf, &
                                                gm_alphaR_vf, &
                                                   kappaR_vf, &
                                          norm_dir, ix,iy,iz  )
            
            
            ! Reshaping inputted data based on dimensional splitting direction
            CALL s_initialize_riemann_solver( qL_prim_vf, kappaL_vf, &
                                              qR_prim_vf, kappaR_vf, &
                                                          q_prim_vf, &
                                               flux_vf, flux_src_vf, &
                                                       flux_gsrc_vf, &
                                                 norm_dir, ix,iy,iz  )
            
            
            ! Computing exact flux and source flux for Euler system of equations
            DO l = is3%beg, is3%end
                DO k = is2%beg, is2%end
                    DO j = is1%beg, is1%end

                        CALL s_compute_constant_states(j,k,l)
            
                        ! Check for pressure positivity condition
                        IF ((G4_L*c_L + G4_R*c_R) < (vel_R(dir_idx(1))-vel_L(dir_idx(1)))) THEN
                            PRINT '(A)','Vacuum is generated by Riemann data. Exiting...'
                            CALL s_mpi_abort()
                        END IF
            
                        CALL s_compute_star_region()
            
                        CALL s_compute_intercell_solution()
            
                        DO i = 1, cont_idx%end
                            flux_rs_vf(i)%sf(j,k,l) = alpha_rho_IC(i)*vel_IC(dir_idx(1))
                        END DO
            
                        DO i = 1, num_dims
                            flux_rs_vf(cont_idx%end+dir_idx(i))%sf(j,k,l) = &
                                rho_IC*vel_IC(dir_idx(1))*vel_IC(dir_idx(i)) + dir_flg(dir_idx(i))*pres_IC
                        END DO
            
                        flux_rs_vf(E_idx)%sf(j,k,l) = vel_IC(dir_idx(1))*(E_IC + pres_IC)
            
                        DO i = 1, adv_idx%end-E_idx
                            flux_rs_vf(E_idx+i)%sf(j,k,l) = alpha_IC(i)*vel_IC(dir_idx(1))
                        END DO
            
                        DO i = 1, num_dims
                            vel_src_rs_vf(dir_idx(i))%sf(j,k,l) = vel_IC(dir_idx(i))
                        END DO
                    END DO
                END DO
            END DO
            
            
            ! Computing the viscous and capillary source flux
            IF(ANY(Re_size > 0)) THEN
                IF (weno_Re_flux) THEN
                    CALL s_compute_viscous_source_flux( &
                                   qL_prim_vf(mom_idx%beg:mom_idx%end), &
                               dqL_prim_dx_vf(mom_idx%beg:mom_idx%end), &
                               dqL_prim_dy_vf(mom_idx%beg:mom_idx%end), &
                               dqL_prim_dz_vf(mom_idx%beg:mom_idx%end), &
                                   qR_prim_vf(mom_idx%beg:mom_idx%end), &
                               dqR_prim_dx_vf(mom_idx%beg:mom_idx%end), &
                               dqR_prim_dy_vf(mom_idx%beg:mom_idx%end), &
                               dqR_prim_dz_vf(mom_idx%beg:mom_idx%end), &
                                       flux_src_vf, norm_dir, ix,iy,iz  )
                ELSE
                    CALL s_compute_viscous_source_flux( &
                                    q_prim_vf(mom_idx%beg:mom_idx%end), &
                               dqL_prim_dx_vf(mom_idx%beg:mom_idx%end), &
                               dqL_prim_dy_vf(mom_idx%beg:mom_idx%end), &
                               dqL_prim_dz_vf(mom_idx%beg:mom_idx%end), &
                                    q_prim_vf(mom_idx%beg:mom_idx%end), &
                               dqR_prim_dx_vf(mom_idx%beg:mom_idx%end), &
                               dqR_prim_dy_vf(mom_idx%beg:mom_idx%end), &
                               dqR_prim_dz_vf(mom_idx%beg:mom_idx%end), &
                                       flux_src_vf, norm_dir, ix,iy,iz  )
                END IF
            END IF
            
            IF(We_size > 0 .AND. We_riemann_flux) THEN
                CALL s_compute_capillary_source_flux( &
                             dqL_prim_dx_vf(adv_idx%beg:adv_idx%end), &
                             dqL_prim_dy_vf(adv_idx%beg:adv_idx%end), &
                             dqL_prim_dz_vf(adv_idx%beg:adv_idx%end), &
                                                    gm_alphaL_vf , &
                             dqR_prim_dx_vf(adv_idx%beg:adv_idx%end), &
                             dqR_prim_dy_vf(adv_idx%beg:adv_idx%end), &
                             dqR_prim_dz_vf(adv_idx%beg:adv_idx%end), &
                                                    gm_alphaR_vf , &
                                 flux_src_vf, norm_dir, ix,iy,iz   )
            END IF
            
            
            ! Reshaping outputted data based on dimensional splitting direction
            CALL s_finalize_riemann_solver( flux_vf, flux_src_vf, &
                                                    flux_gsrc_vf, &
                                              norm_dir, ix,iy,iz  )
            
            
        END SUBROUTINE s_exact_riemann_solver ! --------------------------------
        
        
        
        !>  The procedure assigns and computes the left and right
        !!      states of the Riemann problem
        !! @param j  First coordinate index
        !! @param k Second coordinate index
        !! @param l  Third coordinate index
        SUBROUTINE s_compute_constant_states(j,k,l) ! --------------------------

            INTEGER, INTENT(IN) :: j,k,l

            INTEGER :: i !< Generic loop iterator

            ! Left and Right Riemann States
            DO i = 1, cont_idx%end
                alpha_rho_L(i) = qL_prim_rs_vf(i)%sf( j ,k,l)
                alpha_rho_R(i) = qR_prim_rs_vf(i)%sf(j+1,k,l)
            END DO

            DO i = 1, num_dims
                vel_L(i) = qL_prim_rs_vf(cont_idx%end+i)%sf( j ,k,l)
                vel_R(i) = qR_prim_rs_vf(cont_idx%end+i)%sf(j+1,k,l)
            END DO

            pres_L = qL_prim_rs_vf(E_idx)%sf( j ,k,l)
            pres_R = qR_prim_rs_vf(E_idx)%sf(j+1,k,l)

            CALL s_convert_to_mixture_variables( qL_prim_rs_vf, &
                                                rho_L, gamma_L, &
                                                pi_inf_L, Re_L, &
                                                 We_L,  j ,k,l  )
            CALL s_convert_to_mixture_variables( qR_prim_rs_vf, &
                                                rho_R, gamma_R, &
                                                pi_inf_R, Re_R, &
                                                 We_R, j+1,k,l  )

            E_L = gamma_L*pres_L + pi_inf_L + 5d-1*rho_L*SUM(vel_L**2d0)
            E_R = gamma_R*pres_R + pi_inf_R + 5d-1*rho_R*SUM(vel_R**2d0)

            H_L = (E_L + pres_L)/rho_L
            H_R = (E_R + pres_R)/rho_R

            CALL s_compute_mixture_sound_speeds(j,k,l)

            DO i = 1, crv_size
                alpha_avg_rs_vf(crv_idx(i))%sf(j,k,l) = &
                    5d-1*(alpha_L(crv_idx(i)) + alpha_R(crv_idx(i)))
            END DO

            DO i = 1, 2
                IF (Re_size(i) > 0) THEN
                    Re_avg_rs_vf(i)%sf(j,k,l) = 2d0/(1d0/Re_L(i) + 1d0/Re_R(i))
                END IF
            END DO

            DO i = 1, We_size
                We_avg_rs_vf(We_idx(i,1),We_idx(i,2))%sf(j,k,l) = &
                             5d-1*( We_L(We_idx(i,1),We_idx(i,2)) &
                                  + We_R(We_idx(i,1),We_idx(i,2)) )
            END DO

            DO i = 1, crv_size
                kappa_avg(crv_idx(i)) = &
                    5d-1*( kappaL_rs_vf(crv_idx(i))%sf( j ,k,l) &
                         + kappaR_rs_vf(crv_idx(i))%sf(j+1,k,l) )
            END DO

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
            dpres_We = 0d0
            dpres_L = 0d0
            dpres_R = 0d0

            IF (We_size > 0 .AND. We_riemann_flux) THEN
                DO i = 1, We_size
                    dpres_We = dpres_We &
                             - 1d0/We_avg_rs_vf(We_idx(i,1),We_idx(i,2))%sf(j,k,l) * &
                                                            kappa_avg(We_idx(i,2)) * &
                                                            ( alpha_R(We_idx(i,2))   &
                                                            - alpha_L(We_idx(i,2)) )
                END DO
            END IF

            dpres_L = 5d-1*dpres_We
            dpres_R = -5d-1*dpres_We
        
        END SUBROUTINE s_compute_constant_states ! -----------------------------

        !> Compute mixture sound speed
        !! @param j  First coordinate index
        !! @param k Second coordinate index
        !! @param l  Third coordinate index
        SUBROUTINE s_compute_mixture_sound_speeds(j,k,l) ! ---------------------

            INTEGER, INTENT(IN) :: j,k,l

            REAL(KIND(0d0)) :: blkmod1, blkmod2 !< Fluid bulk modulus for alternate sound speed

            INTEGER :: i !< Generic loop iterator

            IF ( (alt_soundspeed .or. regularization) ) THEN
                DO i = 1, num_fluids
                    alpha_L(i) = qL_prim_rs_vf(E_idx+i)%sf( j ,k,l)
                    alpha_R(i) = qR_prim_rs_vf(E_idx+i)%sf(j+1,k,l)
                END DO

                blkmod1 = ((fluid_pp(1)%gamma+1d0)*pres_L + &
                    fluid_pp(1)%pi_inf)/fluid_pp(1)%gamma
                blkmod2 = ((fluid_pp(2)%gamma+1d0)*pres_L + &
                    fluid_pp(2)%pi_inf)/fluid_pp(2)%gamma
                c_L = 1d0/(rho_L*(alpha_L(1)/blkmod1 + alpha_L(2)/blkmod2))

                blkmod1 = ((fluid_pp(1)%gamma+1d0)*pres_R + &
                    fluid_pp(1)%pi_inf)/fluid_pp(1)%gamma
                blkmod2 = ((fluid_pp(2)%gamma+1d0)*pres_R + &
                    fluid_pp(2)%pi_inf)/fluid_pp(2)%gamma
                c_R = 1d0/(rho_R*(alpha_R(1)/blkmod1 + alpha_R(2)/blkmod2))

            ELSEIF(model_eqns == 3) THEN
                c_L = 0d0
                c_R = 0d0
                DO i = 1, num_fluids
                    c_L = c_L + qL_prim_rs_vf(i+adv_idx%beg-1)%sf( j ,k,l) * (1d0/fluid_pp(i)%gamma+1d0) * &
                        (qL_prim_rs_vf(E_idx)%sf( j ,k,l) + fluid_pp(i)%pi_inf/(fluid_pp(i)%gamma+1d0))
                    c_R = c_R + qR_prim_rs_vf(i+adv_idx%beg-1)%sf(j+1,k,l) * (1d0/fluid_pp(i)%gamma+1d0) * &
                        (qR_prim_rs_vf(E_idx)%sf(j+1,k,l) + fluid_pp(i)%pi_inf/(fluid_pp(i)%gamma+1d0))
                END DO
                c_L = c_L/rho_L
                c_R = c_R/rho_R
            ELSEIF ( (model_eqns == 4) .OR. (model_eqns == 2 .AND. bubbles) ) THEN 
                ! Sound speed for bubble mmixture to order O(\alpha)
                DO i = 1,num_fluids
                    alpha_L(i) = qL_prim_rs_vf(E_idx+i)%sf( j ,k,l)
                    alpha_R(i) = qR_prim_rs_vf(E_idx+i)%sf(j+1,k,l)
                END DO
                
                IF (mpp_lim .AND. (num_fluids > 1)) THEN
                    c_L =   (1d0/gamma_L + 1d0) *   &
                            (pres_L + pi_inf_L) / rho_L  
                    c_R =   (1d0/gamma_R + 1d0) *   &
                            (pres_R + pi_inf_R) / rho_R  
                ELSE
                    c_L =   & 
                            (1d0/gamma_L + 1d0) *   &
                            (pres_L + pi_inf_L) /   &
                            (rho_L*(1d0-alpha_L(num_fluids))) 
                    c_R =   & 
                            (1d0/gamma_R + 1d0) *   &
                            (pres_R + pi_inf_R) /   &
                            (rho_R*(1d0-alpha_R(num_fluids))) 
                END IF
            ELSE
                DO i = 1, num_fluids 
                    alpha_L(i) = qL_prim_rs_vf(E_idx+i)%sf( j ,k,l)
                    alpha_R(i) = qR_prim_rs_vf(E_idx+i)%sf(j+1,k,l)
                END DO
                c_L = ((H_L - 5d-1*SUM(vel_L**2d0))/gamma_L)
                c_R = ((H_R - 5d-1*SUM(vel_R**2d0))/gamma_R)
            END IF

            IF (mixture_err .AND. c_L < 0d0) THEN
                c_L = 100.d0*sgm_eps
            ELSE
                c_L = SQRT(c_L)
            END IF
            IF (mixture_err .AND. c_R < 0d0) THEN
                c_R = 100.d0*sgm_eps
            ELSE
                c_R = SQRT(c_R)
            END IF

            IF (tvd_riemann_flux) THEN
                IF (alt_soundspeed .OR. regularization) THEN
                    DO i = 1, num_fluids
                        lo_alpha_L(i) = q_prim_rs_vf(E_idx+i)%sf( j ,k,l)
                        lo_alpha_R(i) = q_prim_rs_vf(E_idx+i)%sf(j+1,k,l)
                    END DO

                    blkmod1 = ((fluid_pp(1)%gamma+1d0)*lo_pres_L + &
                        fluid_pp(1)%pi_inf)/fluid_pp(1)%gamma
                    blkmod2 = ((fluid_pp(2)%gamma+1d0)*lo_pres_L + &
                        fluid_pp(2)%pi_inf)/fluid_pp(2)%gamma
                    lo_c_L = 1d0/(lo_rho_L*(lo_alpha_L(1)/blkmod1 + lo_alpha_L(2)/blkmod2))

                    blkmod1 = ((fluid_pp(1)%gamma+1d0)*lo_pres_R + &
                        fluid_pp(1)%pi_inf)/fluid_pp(1)%gamma
                    blkmod2 = ((fluid_pp(2)%gamma+1d0)*lo_pres_R + &
                        fluid_pp(2)%pi_inf)/fluid_pp(2)%gamma
                    lo_c_R = 1d0/(lo_rho_R*(lo_alpha_R(1)/blkmod1 + lo_alpha_R(2)/blkmod2))
                ELSE
                    DO i = 1, num_fluids 
                        lo_alpha_L(i) = q_prim_rs_vf(E_idx+i)%sf( j ,k,l)
                        lo_alpha_R(i) = q_prim_rs_vf(E_idx+i)%sf(j+1,k,l)
                    END DO
                    lo_c_L = ((lo_H_L - 5d-1*SUM(lo_vel_L**2d0))/lo_gamma_L)
                    lo_c_R = ((lo_H_R - 5d-1*SUM(lo_vel_R**2d0))/lo_gamma_R)
                END IF

                IF (mixture_err .AND. lo_c_L < 0d0) THEN
                    lo_c_L = sgm_eps
                ELSE
                    lo_c_L = SQRT(lo_c_L)
                END IF
                IF (mixture_err .AND. lo_c_R < 0d0) THEN
                    lo_c_R = sgm_eps
                ELSE
                    lo_c_R = SQRT(lo_c_R)
                END IF
            END IF

        END SUBROUTINE s_compute_mixture_sound_speeds ! ------------------------



        !>  The purpose of this subroutine is to compute the solution
        !!      for pressure and velocity in the star region when using
        !!      the exact Riemann solver
        SUBROUTINE s_compute_star_region() ! -----------------------------------

            REAL(KIND(0d0)) :: change, f_L, dfdp_L, f_R, dfdp_R, pres_old, & 
                        pres_start, pres_tol, vel_diff

            INTEGER :: iter

            pres_tol = 1d-10
            change = -1d0*dflt_real

            ! Compute starting guess value for pressure
            CALL s_guess_pressure(pres_start)

            pres_old = pres_start
            vel_diff = vel_R(dir_idx(1)) - vel_L(dir_idx(1))

            iter = 0

            ! Solve iteratively for pressure in the star region
            DO WHILE (change > pres_tol)
                CALL s_pressure_function(f_L, dfdp_L, pres_old, 1)
                CALL s_pressure_function(f_R, dfdp_R, pres_old, 2)
                pres_S = pres_old - (f_L + f_R + vel_diff)/(dfdp_L + dfdp_R)
                IF (iter > 10000) THEN
                    PRINT '(A)', 'Too many iterations in pressure'
                    CALL s_mpi_abort()
                END IF
                change = 2d0*ABS((pres_S - pres_old)/(pres_S + pres_old))
                pres_old = pres_S
                iter = iter + 1
            END DO

            ! Compute velocity in star region
            vel_S = 5d-1*(vel_L(dir_idx(1)) + vel_R(dir_idx(1)) + f_R - f_L)

        END SUBROUTINE s_compute_star_region ! ----------------------------------

        !>  The purpose of this subroutine is to evaluate the pressure
        !!      functions f_K in the exact Riemann solver
        !!  @param f_K    Pressure function
        !!  @param dfdp_K Mixture pressure derivative
        !!  @param pres   Pressure
        !!  @param side   Wave side
        SUBROUTINE s_pressure_function(f_K, dfdp_K, pres, side) ! -


            REAL(KIND(0d0)), INTENT(IN) :: pres
            INTEGER, INTENT(IN) :: side
            REAL(KIND(0d0)), INTENT(OUT) :: f_K, dfdp_K
            REAL(KIND(0d0)) :: gam_L, pinf_L, gam_R, pinf_R, c_SL, c_SR

            IF (side == 1) THEN
                gam_L = (gamma_L+1d0)/gamma_L
                pinf_L = pi_inf_L/(gamma_L+1d0)
                IF (pres + dpres_L <= pres_L) THEN
                   ! Rarefaction wave
                    c_SL = c_L*((pres+dpres_L+pinf_L)/(pres_L+pinf_L))**G1_L
                    f_K = G4_L*(c_SL/c_L - 1d0)*c_L
                    dfdp_K = c_SL/(gam_L*(dpres_L+pres+pinf_L))
                ELSEIF (pres + dpres_L > pres_L) THEN
                    ! Shock wave
                    f_K = (c_L/gam_L*((pres+dpres_L)/pres_L - 1d0)*(pres_L/(pres_L+pinf_L))) / &
                            SQRT(G2_L*((pres+dpres_L)/pres_L - 1d0)*(pres_L/(pres_L+pinf_L)) + 1d0)
                    dfdp_K = 2d0*c_L/gam_L/(pres_L+pinf_L) / &
                                SQRT(2d0*(gam_L+1d0)*(pres+dpres_L-pres_L)/gam_L/(pres_L+pinf_L)+4d0) - &
                             2d0*c_L*(pres+dpres_L-pres_L)*(gam_L+1d0)/gam_L**2d0/(pres_L+pinf_L)**2d0 / &
                                SQRT(2d0*(gam_L+1d0)*(pres+dpres_L-pres_L)/gam_L/(pres_L+pinf_L)+4d0)**3d0
                ELSE
                    PRINT '(A)','Error in evaluating left pressure function. Exiting...'
                    CALL s_mpi_abort()
                END IF
            ELSEIF (side == 2) THEN
                gam_R = (gamma_R+1d0)/gamma_R
                pinf_R = pi_inf_R/(gamma_R+1d0)
                IF (pres + dpres_R <= pres_R) THEN
                    ! Rarefaction wave
                    c_SR = c_R*((pres+dpres_R+pinf_R)/(pres_R+pinf_R))**G1_R
                    f_K = G4_R*(c_SR/c_R - 1d0)*c_R
                    dfdp_K = c_SR/(gam_R*(dpres_R+pres+pinf_R))
                ELSEIF (pres + dpres_R > pres_R) THEN
                    ! Shock wave
                    f_K = (c_R/gam_R*((pres+dpres_R)/pres_R - 1d0)*(pres_R/(pres_R+pinf_R))) / &
                            SQRT(G2_R*((pres+dpres_R)/pres_R - 1d0)*(pres_R/(pres_R+pinf_R)) + 1d0)
                    dfdp_K = 2d0*c_R/gam_R/(pres_R+pinf_R) / &
                                SQRT(2d0*(gam_R+1d0)*(pres+dpres_R-pres_R)/gam_R/(pres_R+pinf_R)+4d0) - &
                             2d0*c_R*(pres+dpres_R-pres_R)*(gam_R+1d0)/gam_R**2d0/(pres_R+pinf_R)**2d0 / &
                                SQRT(2d0*(gam_R+1d0)*(pres+dpres_R-pres_R)/gam_R/(pres_R+pinf_R)+4d0)**3d0
                ELSE
                    PRINT '(A)','Error in evaluating right pressure function. Exiting...'
                    CALL s_mpi_abort()
                END IF
            END IF

        END SUBROUTINE s_pressure_function ! ------------------------------------




        !>  The purpose of this subroutine is to provide a guess value
        !!      for pressure in the star region. The choice is made
        !!      according to adaptive Riemann solver using the PVRS, TRRS.
        !!      and TSRS approximate Riemann solvers.
        !!  @param pres_start Initial and output pressure
        SUBROUTINE s_guess_pressure(pres_start) ! -------------------------------


            REAL(KIND(0d0)), INTENT(INOUT) :: pres_start

            REAL(KIND(0d0)) :: CUP, pres_max, pres_min, pres_PV, & 
                                Q_max, Q_user, pres_TS, pres_TR
            REAL(KIND(0d0)) :: pinf_L, pinf_R,A,B

            Q_user = 2d0

            ! Compute guess pressure from PVRS Riemann solver
            CUP = 25d-2*(rho_L + rho_R)*(c_L + c_R)
            pres_PV = 5d-1*(pres_L + pres_R) + 5d-1*(vel_L(dir_idx(1)) - vel_R(dir_idx(1)))*CUP
            pres_PV = MAX(0d0, pres_PV)
            pres_min = MIN(pres_L, pres_R)
            pres_max = MAX(pres_L, pres_R)
            Q_max = pres_max / pres_min 

            IF ((Q_max <= Q_user) .AND. (pres_min <= pres_PV) .AND. (pres_PV <= pres_max)) THEN

                ! Select PVRS Riemann solver
                pres_start = pres_PV

                IF (pres_start /= pres_start) THEN
                    PRINT '(A)','NaN guess for pressure using PVRS'
                    CALL s_mpi_abort()
                END IF

            ELSEIF (pres_PV < pres_min) THEN

                ! Select TRRS Riemann solver
                pinf_L = pi_inf_L/(gamma_L + 1d0)
                pinf_R = pi_inf_R/(gamma_R + 1d0)
                pres_TR = ((((-G4_R*c_R*(((pres_PV+dpres_R+pinf_R)/(pres_R+pinf_R))**G1_R - 1d0) - &
                            (vel_R(dir_idx(1)) - vel_L(dir_idx(1)))) / c_L/G4_L )+1d0)**G3_L) * &
                            (pres_L+pinf_L) - (dpres_L+pinf_L)
                pres_start = MAX(0d0,pres_TR)

                IF (pres_start /= pres_start) THEN
                    PRINT '(A)','NaN guess for pressure using TRRS'
                    CALL s_mpi_abort()
                END IF

            ELSE
                ! Select TSRS Riemann solver with pres_PV as estimate
                A = SQRT(G2_L*((pres_PV+dpres_L-pres_L)/(pres_L+pinf_L)) + 1d0)
                B = SQRT(G2_R*((pres_PV+dpres_R-pres_R)/(pres_R+pinf_R)) + 1d0)
                pres_TS = (c_L*gamma_L/(gamma_L+1d0)/A*((pres_L-dpres_L)/(pres_L+pinf_L)) + &
                           c_R*gamma_R/(gamma_R+1d0)/B*((pres_R-dpres_R)/(pres_R+pinf_R)) - &
                          (vel_R(dir_idx(1)) - vel_L(dir_idx(1)))) / &
                          (c_L*gamma_L/(gamma_L+1d0)/A/(pres_L+pinf_L) + &
                           c_R*gamma_R/(gamma_R+1d0)/B/(pres_R+pinf_R))
                pres_start = MAX(0d0,pres_TS)

                IF (pres_start /= pres_start) THEN
                    PRINT '(A)','NaN guess for pressure using TSRS'
                    CALL s_mpi_abort()
                END IF

            END IF

            
        END SUBROUTINE s_guess_pressure ! --------------------------------------


        !> Computes the averaged intrercell variables for the Riemann solver
        SUBROUTINE s_compute_intercell_solution() ! --------------------

            INTEGER :: i

            REAL(KIND(0d0)) :: c_IC
            REAL(KIND(0d0)) :: s_HL, S_TL, c_SL, pres_SL, s_L
            REAL(KIND(0d0)) :: s_HR, S_TR, c_SR, pres_SR, s_R

            IF (0d0 <= vel_S) THEN
                ! IC lies to the left of the contact discontinuity
                IF (pres_S  + dpres_L <= pres_L) THEN
                    ! Left rarefaction
                    s_HL = vel_L(dir_idx(1)) - c_L

                    IF (0d0 <= s_HL) THEN
                        ! IC is left data state
                        DO i = 1, cont_idx%end
                            alpha_rho_IC(i) = alpha_rho_L(i)
                        END DO
                        rho_IC = rho_L

                        DO i = 1, num_dims
                            vel_IC(i) = vel_L(i)
                        END DO

                        pres_IC = pres_L
                        E_IC = E_L

                        DO i = 1, num_fluids 
                            alpha_IC(i) = alpha_L(i)
                        END DO

                    ELSE
                        c_SL = c_L*((pres_S + dpres_L + pi_inf_L/(gamma_L + 1d0))/(pres_L + pi_inf_L/(gamma_L + 1d0)))**G1_L
                        S_TL = vel_S - c_SL

                        IF (0d0 > S_TL) THEN
                            ! IC is star left state
                            DO i = 1, cont_idx%end
                                alpha_rho_IC(i) = alpha_rho_L(i)*((pres_S + dpres_L + pi_inf_L/(gamma_L + 1d0))/ &
                                    (pres_L + pi_inf_L/(gamma_L + 1d0)))**(gamma_L/(gamma_L+1d0))
                            END DO
                            rho_IC = rho_L*((pres_S + dpres_L + pi_inf_L/(gamma_L + 1d0)) / &
                                (pres_L + pi_inf_L/(gamma_L + 1d0)))**(gamma_L/(gamma_L+1d0))

                            vel_IC(dir_idx(1)) = vel_S
                            DO i = 2, num_dims
                                vel_IC(dir_idx(i)) = vel_L(dir_idx(i))
                            END DO

                            pres_IC = pres_S + dpres_L
                            E_IC = gamma_L*pres_IC + pi_inf_L + 5d-1*rho_IC*SUM(vel_IC**2d0) 

                            DO i = 1, num_fluids 
                                alpha_IC(i) = alpha_L(i)
                            END DO

                        ELSE
                            ! IC is inside left rarefaction
                            vel_IC(dir_idx(1)) = G5_L*(c_L + G7_L*vel_L(dir_idx(1)) + 0d0)
                            c_IC = G5_L*(c_L + G7_L*(vel_L(dir_idx(1)) - 0d0))

                            DO i = 1, cont_idx%end
                                alpha_rho_IC(i) = alpha_rho_L(i)*(c_IC/c_L)**G4_L
                            END DO
                            rho_IC = rho_L*(c_IC/c_L)**G4_L

                            DO i = 2, num_dims
                                vel_IC(dir_idx(i)) = vel_L(dir_idx(i))
                            END DO

                            pres_IC = (pres_L + pi_inf_L/(gamma_L + 1d0))*(c_IC/c_L)**G3_L - (pi_inf_L/(gamma_L + 1d0))
                            E_IC = gamma_L*pres_IC + pi_inf_L + 5d-1*rho_IC*SUM(vel_IC**2d0) 

                            DO i = 1, num_fluids 
                                alpha_IC(i) = alpha_L(i)
                            END DO

                        END IF
                    END IF
                ELSE
                    ! Left shock
                    pres_SL = (pres_S + dpres_L + pi_inf_L/(gamma_L + 1d0))/(pres_L + pi_inf_L/(gamma_L + 1d0))
                    s_L = vel_L(dir_idx(1)) - c_L*SQRT(G2_L*(pres_S+dpres_L-pres_L)/(pres_L+pi_inf_L/(gamma_L+1d0)) +1d0)

                    IF (0d0 <= s_L) THEN
                        ! IC is left data state
                        DO i = 1, cont_idx%end
                            alpha_rho_IC(i) = alpha_rho_L(i)
                        END DO
                        rho_IC = rho_L

                        DO i = 1, num_dims
                            vel_IC(i) = vel_L(i)
                        END DO

                        pres_IC = pres_L
                        E_IC = E_L

                        DO i = 1, num_fluids 
                            alpha_IC(i) = alpha_L(i)
                        END DO

                    ELSE
                        ! IC is star left state
                        DO i = 1, cont_idx%end
                            alpha_rho_IC(i) = alpha_rho_L(i)*(pres_SL + G6_L)/(pres_SL*G6_L + 1d0)
                        END DO
                        rho_IC = rho_L*(pres_SL + G6_L)/(pres_SL*G6_L + 1d0)

                        vel_IC(dir_idx(1)) = vel_S
                        DO i = 2, num_dims
                            vel_IC(dir_idx(i)) = vel_L(dir_idx(i))
                        END DO

                        pres_IC = pres_S + dpres_L
                        E_IC = gamma_L*pres_IC + pi_inf_L + 5d-1*rho_IC*SUM(vel_IC**2d0) 

                        DO i = 1, num_fluids 
                            alpha_IC(i) = alpha_L(i)
                        END DO

                    END IF
                END IF
            ELSE
                ! IC is to the right of the contact discontinuity
                IF (pres_S + dpres_R > pres_R) THEN
                    ! Right shock
                    pres_SR = (pres_S + dpres_R + pi_inf_R/(gamma_R + 1d0))/(pres_R + pi_inf_R/(gamma_R + 1d0))
                    s_R = vel_R(dir_idx(1)) + c_R*SQRT(G2_R*(pres_S+dpres_R-pres_R)/(pres_R+pi_inf_R/(gamma_R+1d0)) +1d0)

                    IF (0d0 >= s_R) THEN
                        ! IC is right data state
                        DO i = 1, cont_idx%end
                            alpha_rho_IC(i) = alpha_rho_R(i)
                        END DO
                        rho_IC = rho_R

                        DO i = 1, num_dims
                            vel_IC(i) = vel_R(i)
                        END DO

                        pres_IC = pres_R
                        E_IC = E_R

                        DO i = 1, num_fluids 
                            alpha_IC(i) = alpha_R(i)
                        END DO

                    ELSE
                        ! IC is star right state
                        DO i = 1, cont_idx%end
                            alpha_rho_IC(i) = alpha_rho_R(i)*(pres_SR + G6_R)/(pres_SR*G6_R + 1d0)
                        END DO
                        rho_IC = rho_R*(pres_SR + G6_R)/(pres_SR*G6_R + 1d0)

                        vel_IC(dir_idx(1)) = vel_S
                        DO i = 2, num_dims
                            vel_IC(dir_idx(i)) = vel_R(dir_idx(i))
                        END DO

                        pres_IC = pres_S + dpres_R
                        E_IC = gamma_R*pres_IC + pi_inf_R + 5d-1*rho_IC*SUM(vel_IC**2d0) 

                        DO i = 1, num_fluids 
                            alpha_IC(i) = alpha_R(i)
                        END DO

                    END IF
                ELSE
                    ! Right rarefaction
                    s_HR = vel_R(dir_idx(1)) + c_R

                    IF (0d0 >= s_HR) THEN
                        ! IC is right data state
                        DO i = 1, cont_idx%end
                            alpha_rho_IC(i) = alpha_rho_R(i)
                        END DO
                        rho_IC = rho_R

                        DO i = 1, num_dims
                            vel_IC(i) = vel_R(i)
                        END DO

                        pres_IC = pres_R
                        E_IC = E_R

                        DO i = 1, num_fluids 
                            alpha_IC(i) = alpha_R(i)
                        END DO

                    ELSE
                        c_SR = c_R*((pres_S + dpres_R + pi_inf_R/(gamma_R + 1d0))/(pres_R + pi_inf_R/(gamma_R + 1d0)))**G1_R
                        S_TR = vel_S + c_SR

                        IF (0d0 <= S_TR) THEN
                            ! IC is star right state
                            DO i = 1, cont_idx%end
                                alpha_rho_IC(i) = alpha_rho_R(i)*((pres_S + dpres_R + pi_inf_R/(gamma_R + 1d0))/ &
                                    (pres_R + pi_inf_R/(gamma_R + 1d0)))**(gamma_R/(gamma_R+1d0))
                            END DO
                            rho_IC = rho_R*((pres_S + dpres_R + pi_inf_R/(gamma_R + 1d0))/(pres_R + &
                                pi_inf_R/(gamma_R + 1d0)))**(gamma_R/(gamma_R+1d0))

                            vel_IC(dir_idx(1)) = vel_S
                            DO i = 2, num_dims
                                vel_IC(dir_idx(i)) = vel_R(dir_idx(i))
                            END DO

                            pres_IC = pres_S + dpres_R
                            E_IC = gamma_R*pres_IC + pi_inf_R + 5d-1*rho_IC*SUM(vel_IC**2d0) 

                            DO i = 1, num_fluids 
                                alpha_IC(i) = alpha_R(i)
                            END DO

                        ELSE
                            ! IC is inside right rarefaction
                            vel_IC(dir_idx(1)) = G5_R*(-1d0*c_R + G7_R*vel_R(dir_idx(1)) + 0d0)
                            c_IC = G5_R*(c_R - G7_R*(vel_R(dir_idx(1)) - 0d0))

                            DO i = 1, cont_idx%end
                                alpha_rho_IC(i) = alpha_rho_R(i)*(c_IC/c_R)**G4_R
                            END DO
                            rho_IC = rho_R*(c_IC/c_R)**G4_R

                            DO i = 2, num_dims
                                vel_IC(dir_idx(i)) = vel_R(dir_idx(i))
                            END DO

                            pres_IC = (pres_R + pi_inf_R/(gamma_R + 1d0))*(c_IC/c_R)**G3_R - (pi_inf_R/(gamma_R + 1d0))
                            E_IC = gamma_R*pres_IC + pi_inf_R + 5d-1*rho_IC*SUM(vel_IC**2d0) 

                            DO i = 1, num_fluids 
                                alpha_IC(i) = alpha_R(i)
                            END DO

                        END IF
                    END IF
                END IF
            END IF

        END SUBROUTINE s_compute_intercell_solution ! --------------------------



        !>  The procedure computes the Roe average density, velocity,
        !!      enthalpy, volume fractions, specific heat ratio function,
        !!      speed of sound, shear and volume Reynolds numbers, Weber
        !!      numbers and curvatures, at the cell-boundaries, from the
        !!      left and right states of the Riemann problem.
        !! @param j  First coordinate index
        !! @param k Second coordinate index
        !! @param l  Third coordinate index
        SUBROUTINE s_compute_roe_average_state(j,k,l) ! ---------------

            INTEGER, INTENT(IN) :: j,k,l

            INTEGER :: i
            
            ! Left and Right Riemann Problem States ============================
            DO i = 1, cont_idx%end
                alpha_rho_L(i) = qL_prim_rs_vf(i)%sf( j ,k,l)
                alpha_rho_R(i) = qR_prim_rs_vf(i)%sf(j+1,k,l)
            END DO
            
            DO i = 1, num_dims
                vel_L(i) = qL_prim_rs_vf(cont_idx%end+i)%sf( j ,k,l)
                vel_R(i) = qR_prim_rs_vf(cont_idx%end+i)%sf(j+1,k,l)
            END DO
            
            pres_L = qL_prim_rs_vf(E_idx)%sf( j ,k,l)
            pres_R = qR_prim_rs_vf(E_idx)%sf(j+1,k,l)
            
            CALL s_convert_to_mixture_variables(  qL_prim_rs_vf, &
                                                 rho_L, gamma_L, &
                                                 pi_inf_L, Re_L, &
                                                  We_L,  j ,k,l  )
            CALL s_convert_to_mixture_variables(  qR_prim_rs_vf, &
                                                 rho_R, gamma_R, &
                                                 pi_inf_R, Re_R, &
                                                  We_R, j+1,k,l  )
            
            E_L = gamma_L*pres_L + pi_inf_L + 5d-1*rho_L*SUM(vel_L**2d0)
            E_R = gamma_R*pres_R + pi_inf_R + 5d-1*rho_R*SUM(vel_R**2d0)
            
            H_L = (E_L + pres_L)/rho_L
            H_R = (E_R + pres_R)/rho_R

            IF (hypoelasticity) THEN
                DO i = 1, (num_dims*(num_dims+1)) / 2
                    tau_e_L(i) = qL_prim_rs_vf(stress_idx%beg-1+i)%sf( j ,k,l)
                    tau_e_R(i) = qR_prim_rs_vf(stress_idx%beg-1+i)%sf(j+1,k,l)
                END DO
            END IF

            IF (tvd_riemann_flux) THEN
                DO i = 1, cont_idx%end
                    lo_alpha_rho_L(i) = q_prim_rs_vf(i)%sf( j ,k,l)
                    lo_alpha_rho_R(i) = q_prim_rs_vf(i)%sf(j+1,k,l)
                END DO

                DO i = 1, num_dims
                    lo_vel_L(i) = q_prim_rs_vf(cont_idx%end+i)%sf( j ,k,l)
                    lo_vel_R(i) = q_prim_rs_vf(cont_idx%end+i)%sf(j+1,k,l)
                END DO

                lo_pres_L = q_prim_rs_vf(E_idx)%sf( j ,k,l)
                lo_pres_R = q_prim_rs_vf(E_idx)%sf(j+1,k,l)

                CALL s_convert_to_mixture_variables( q_prim_rs_vf, &
                                lo_rho_L, lo_gamma_L,&
                                lo_pi_inf_L, lo_Re_L,&
                                lo_We_L,  j , k, l)
                CALL s_convert_to_mixture_variables( q_prim_rs_vf, &
                                lo_rho_R, lo_gamma_R,&
                                lo_pi_inf_R, lo_Re_R,&
                                lo_We_R, j+1, k, l)

                lo_E_L = lo_gamma_L*lo_pres_L + lo_pi_inf_L + 5d-1*lo_rho_L*SUM(lo_vel_L**2d0)
                lo_E_R = lo_gamma_R*lo_pres_R + lo_pi_inf_R + 5d-1*lo_rho_R*SUM(lo_vel_R**2d0)
                
                lo_H_L = (lo_E_L + lo_pres_L)/lo_rho_L
                lo_H_R = (lo_E_R + lo_pres_R)/lo_rho_R

                IF (hypoelasticity) THEN
                    DO i = 1, (num_dims*(num_dims+1)) / 2
                        lo_tau_e_L(i) = q_prim_rs_vf(stress_idx%beg-1+i)%sf( j ,k,l)
                        lo_tau_e_R(i) = q_prim_rs_vf(stress_idx%beg-1+i)%sf(j+1,k,l)
                    END DO
                END IF
    
            END IF
            
            CALL s_compute_mixture_sound_speeds(j,k,l)

            ! ==================================================================
            
            
            ! Roe Average Riemann Problem State ================================
            rho_avg = SQRT(rho_L*rho_R)
            
            vel_avg = (SQRT(rho_L)*vel_L + SQRT(rho_R)*vel_R) / &
                      (SQRT(rho_L) + SQRT(rho_R))
            
            H_avg = (SQRT(rho_L)*H_L + SQRT(rho_R)*H_R) / &
                    (SQRT(rho_L) + SQRT(rho_R))
            
            DO i = 1, crv_size
               alpha_avg_rs_vf(crv_idx(i))%sf(j,k,l) = &
                        5d-1*(alpha_L(crv_idx(i)) + alpha_R(crv_idx(i)))
            END DO
            
            gamma_avg = (SQRT(rho_L)*gamma_L + SQRT(rho_R)*gamma_R) / &
                        (SQRT(rho_L) + SQRT(rho_R))
            
            IF (mixture_err) THEN
                IF ((H_avg - 5d-1*SUM(vel_avg**2d0)) < 0d0) THEN
                    c_avg = sgm_eps
                ELSE
                    c_avg = SQRT((H_avg - 5d-1*SUM(vel_avg**2d0))/gamma_avg)
                END IF
            ELSE
                c_avg = SQRT((H_avg - 5d-1*SUM(vel_avg**2d0))/gamma_avg)
            END IF
            
            DO i = 1,2
               IF(Re_size(i) > 0) THEN
                  Re_avg_rs_vf(i)%sf(j,k,l) = 2d0/(1d0/Re_L(i) + 1d0/Re_R(i))
               END IF
            END DO
            
            DO i = 1, We_size
               We_avg_rs_vf(We_idx(i,1),We_idx(i,2))%sf(j,k,l) = &
                            5d-1*( We_L(We_idx(i,1),We_idx(i,2)) &
                                 + We_R(We_idx(i,1),We_idx(i,2)) )
            END DO
            
            DO i = 1, crv_size
               kappa_avg(crv_idx(i)) = &
                         5d-1*( kappaL_rs_vf(crv_idx(i))%sf( j ,k,l) &
                              + kappaR_rs_vf(crv_idx(i))%sf(j+1,k,l) )
            END DO
            ! ==================================================================
            
            
        END SUBROUTINE s_compute_roe_average_state ! ---------------------------
        
        
        !>  This procedure calculates the arithmetic average density,
        !!      velocity, enthalpy, volume fractions, specIFic heat ratio
        !!      function, sound speed, shear and volume Reynolds numbers,
        !!      Weber numbers and the curvatures, at the cell-boundaries,
        !!      from the left and right states of the Riemann problem.
        !!  @param j  First coordinate index
        !!  @param k Second coordinate index
        !!  @param l  Third coordinate index
        SUBROUTINE s_compute_arithmetic_average_state(j,k,l) ! --------
           
            INTEGER, INTENT(IN) :: j,k,l
            
            INTEGER :: i,q !< Generic loop iterator

            !ensemble-averaged bubble variables
            REAL(KIND(0d0)) :: PbwR3Lbar, Pbwr3Rbar
            REAL(KIND(0d0)) :: R3Lbar, R3Rbar
            REAL(KIND(0d0)) :: R3V2Lbar, R3V2Rbar

            ! Left and Right Riemann Problem States ============================
            DO i = 1, cont_idx%end
               alpha_rho_L(i) = qL_prim_rs_vf(i)%sf( j ,k,l)
               alpha_rho_R(i) = qR_prim_rs_vf(i)%sf(j+1,k,l)
            END DO
            
            DO i = 1, num_dims
               vel_L(i) = qL_prim_rs_vf(cont_idx%end+i)%sf( j ,k,l)
               vel_R(i) = qR_prim_rs_vf(cont_idx%end+i)%sf(j+1,k,l)
            END DO
         
            CALL s_convert_to_mixture_variables(  qL_prim_rs_vf, &
                                                 rho_L, gamma_L, &
                                                 pi_inf_L, Re_L, &
                                                  We_L,  j ,k,l  )
            CALL s_convert_to_mixture_variables(  qR_prim_rs_vf, &
                                                 rho_R, gamma_R, &
                                                 pi_inf_R, Re_R, &
                                                  We_R, j+1,k,l  )

            pres_L = qL_prim_rs_vf(E_idx)%sf( j ,k,l)
            pres_R = qR_prim_rs_vf(E_idx)%sf(j+1,k,l)

            E_L = gamma_L*pres_L + pi_inf_L + 5d-1*rho_L*SUM(vel_L**2d0)
            E_R = gamma_R*pres_R + pi_inf_R + 5d-1*rho_R*SUM(vel_R**2d0)
            
            H_L = (E_L + pres_L)/rho_L
            H_R = (E_R + pres_R)/rho_R

            IF (hypoelasticity) THEN
                DO i = 1, (num_dims*(num_dims+1)) / 2
                    tau_e_L(i) = qL_prim_rs_vf(stress_idx%beg-1+i)%sf( j ,k,l)
                    tau_e_R(i) = qR_prim_rs_vf(stress_idx%beg-1+i)%sf(j+1,k,l)
                END DO
            END IF

            ! Compute left/right states for bubble number density
            IF (bubbles) THEN
                DO i = 1,num_fluids
                    alpha_L(i) = qL_prim_rs_vf(E_idx+i)%sf( j ,k,l)
                    alpha_R(i) = qR_prim_rs_vf(E_idx+i)%sf(j+1,k,l)
                END DO

                DO i = 1,nb
                    R0_L(i) = qL_prim_rs_vf(bub_idx%rs(i))%sf(j,k,l)
                    R0_R(i) = qR_prim_rs_vf(bub_idx%rs(i))%sf(j+1,k,l)

                    V0_L(i) = qL_prim_rs_vf(bub_idx%vs(i))%sf(j,k,l)
                    V0_R(i) = qR_prim_rs_vf(bub_idx%vs(i))%sf(j+1,k,l)
                    IF (.NOT. polytropic) THEN
                        P0_L(i) = qL_prim_rs_vf(bub_idx%ps(i))%sf(j,k,l)
                        P0_R(i) = qR_prim_rs_vf(bub_idx%ps(i))%sf(j+1,k,l)
                    END IF
                END DO
                
                CALL s_comp_n_from_prim(alpha_L(num_fluids),R0_L,nbub_L)
                CALL s_comp_n_from_prim(alpha_R(num_fluids),R0_R,nbub_R)
                
                DO i = 1,nb
                    IF ( .NOT. qbmm ) THEN
                        IF (polytropic) THEN
                            pbw_L(i) = f_cpbw_KM(R0(i),R0_L(i),V0_L(i),0d0)
                            pbw_R(i) = f_cpbw_KM(R0(i),R0_R(i),V0_R(i),0d0)
                        ELSE
                            pbw_L(i) = f_cpbw_KM(R0(i),R0_L(i),V0_L(i),P0_L(i))
                            pbw_R(i) = f_cpbw_KM(R0(i),R0_R(i),V0_R(i),P0_R(i))
                        END IF
                    END IF
                END DO

                IF (qbmm) THEN
                    PbwR3Lbar = mom_sp(4)%sf(j  ,k,l)
                    PbwR3Rbar = mom_sp(4)%sf(j+1,k,l)

                    R3Lbar = mom_sp(1)%sf(j  ,k,l)
                    R3Rbar = mom_sp(1)%sf(j+1,k,l)

                    R3V2Lbar = mom_sp(3)%sf(j  ,k,l)
                    R3V2Rbar = mom_sp(3)%sf(j+1,k,l)
                ELSE
                    CALL s_quad(pbw_L*(R0_L**3.d0), PbwR3Lbar)
                    CALL s_quad(pbw_R*(R0_R**3.d0), PbwR3Rbar)

                    CALL s_quad(R0_L**3.d0, R3Lbar)
                    CALL s_quad(R0_R**3.d0, R3Rbar)
                
                    CALL s_quad((R0_L**3.d0)*(V0_L**2.d0), R3V2Lbar)
                    CALL s_quad((R0_R**3.d0)*(V0_R**2.d0), R3V2Rbar)
                END IF 

                !ptilde = \alf( pl - \bar{ pbw R^3)/\bar{R^3} - rho \bar{R^3 \Rdot^2}/\bar{R^3} ) 
                IF (alpha_L(num_fluids) < small_alf .OR. R3Lbar < small_alf) THEN
                    ptilde_L = alpha_L(num_fluids)*pres_L
                ELSE
                    ptilde_L = alpha_L(num_fluids)*(pres_L - PbwR3Lbar/R3Lbar - & 
                        rho_L*R3V2Lbar/R3Lbar )
                END IF

                IF (alpha_R(num_fluids) < small_alf .OR. R3Rbar < small_alf) THEN
                    ptilde_R = alpha_R(num_fluids)*pres_R
                ELSE
                    ptilde_R = alpha_R(num_fluids)*(pres_R - PbwR3Rbar/R3Rbar - & 
                        rho_R*R3V2Rbar/R3Rbar )
                END IF
                
                IF ( (ptilde_L .ne. ptilde_L) .or. (ptilde_R .ne. ptilde_R) ) THEN
                    PRINT*, 'Ptilde NaN at ',j,k,l,x_cb(j)
                    PRINT*, nbub_L, alpha_L, pres_L, PbwR3Lbar, R3Lbar, rho_L, R3V2Lbar, R3Lbar
                    PRINT*, nbub_R, alpha_R, pres_R, PbwR3Rbar, R3Rbar, rho_R, R3V2Rbar, R3Rbar
                    CALL s_mpi_abort()
                END IF

                ptil(j,k,l) = 0.5d0*(ptilde_L+ptilde_R)
            END IF

            IF (tvd_riemann_flux) THEN
                DO i = 1, cont_idx%end
                    lo_alpha_rho_L(i) = q_prim_rs_vf(i)%sf( j ,k,l)
                    lo_alpha_rho_R(i) = q_prim_rs_vf(i)%sf(j+1,k,l)
                END DO
    
                DO i = 1, num_dims
                    lo_vel_L(i) = q_prim_rs_vf(cont_idx%end+i)%sf( j ,k,l)
                    lo_vel_R(i) = q_prim_rs_vf(cont_idx%end+i)%sf(j+1,k,l)
                END DO
    
                lo_pres_L = q_prim_rs_vf(E_idx)%sf( j ,k,l)
                lo_pres_R = q_prim_rs_vf(E_idx)%sf(j+1,k,l)
                
                CALL s_convert_to_mixture_variables( q_prim_rs_vf, &
                                lo_rho_L, lo_gamma_L,&
                                lo_pi_inf_L, lo_Re_L,&
                                lo_We_L,  j , k, l)
                CALL s_convert_to_mixture_variables( q_prim_rs_vf, &
                                lo_rho_R, lo_gamma_R,&
                                lo_pi_inf_R, lo_Re_R,&
                                lo_We_R, j+1, k, l)
            
                lo_E_L = lo_gamma_L*lo_pres_L + lo_pi_inf_L + 5d-1*lo_rho_L*SUM(lo_vel_L**2d0)
                lo_E_R = lo_gamma_R*lo_pres_R + lo_pi_inf_R + 5d-1*lo_rho_R*SUM(lo_vel_R**2d0)
            
                lo_H_L = (lo_E_L + lo_pres_L)/lo_rho_L
                lo_H_R = (lo_E_R + lo_pres_R)/lo_rho_R

                IF (hypoelasticity) THEN
                    DO i = 1, (num_dims*(num_dims+1)) / 2
                        lo_tau_e_L(i) = q_prim_rs_vf(stress_idx%beg-1+i)%sf( j ,k,l)
                        lo_tau_e_R(i) = q_prim_rs_vf(stress_idx%beg-1+i)%sf(j+1,k,l)
                    END DO
                END IF

            END IF


            CALL s_compute_mixture_sound_speeds(j,k,l)

            ! ==================================================================
            
            
            ! Arithmetic Average Riemann Problem State =========================
            rho_avg = 5d-1*(rho_L + rho_R)
            
            vel_avg = 5d-1*(vel_L + vel_R)
            
            H_avg = 5d-1*(H_L + H_R)
            
            DO i = 1, crv_size
               alpha_avg_rs_vf(crv_idx(i))%sf(j,k,l) = &
                        5d-1*(alpha_L(crv_idx(i)) + alpha_R(crv_idx(i)))
            END DO
            
            gamma_avg = 5d-1*(gamma_L + gamma_R)
            
            IF (mixture_err) THEN
                IF ((H_avg - 5d-1*SUM(vel_avg**2d0)) < 0d0) THEN
                    c_avg = sgm_eps
                ELSE
                    c_avg = SQRT((H_avg - 5d-1*SUM(vel_avg**2d0))/gamma_avg)
                END IF
            ELSE
                c_avg = SQRT((H_avg - 5d-1*SUM(vel_avg**2d0))/gamma_avg)
            END IF
            
            DO i = 1,2
               IF(Re_size(i) > 0) THEN
                  Re_avg_rs_vf(i)%sf(j,k,l) = 2d0/(1d0/Re_L(i) + 1d0/Re_R(i))
               END IF
            END DO
            
            DO i = 1, We_size
               We_avg_rs_vf(We_idx(i,1),We_idx(i,2))%sf(j,k,l) = &
                            5d-1*( We_L(We_idx(i,1),We_idx(i,2)) &
                                 + We_R(We_idx(i,1),We_idx(i,2)) )
            END DO
            
            DO i = 1, crv_size
               kappa_avg(crv_idx(i)) = &
                         5d-1*( kappaL_rs_vf(crv_idx(i))%sf( j ,k,l) &
                              + kappaR_rs_vf(crv_idx(i))%sf(j+1,k,l) )
            END DO
            ! ==================================================================
           

            
        END SUBROUTINE s_compute_arithmetic_average_state ! --------------------
        
        
        !>  The direct estimation of the left, right and middle wave
        !!      speeds, proposed by Batten et al. (1997) that results in
        !!      the exact resolution of isolated shock and contact waves.
        !!  @param j  First coordinate index
        !!  @param k Second coordinate index
        !!  @param l  Third coordinate index       
        SUBROUTINE s_compute_direct_wave_speeds(j,k,l) ! -----------------------
           
            INTEGER, INTENT(IN) :: j,k,l
            

            REAL(KIND(0d0)) :: denom, dpres_We !< Capillary pressure
            

            INTEGER :: i !< Generic loop iterator
            
            
            dpres_We = 0d0
            
            IF (We_size > 0 .AND. We_wave_speeds .AND. We_riemann_flux) THEN
                DO i = 1, We_size
                  dpres_We = dpres_We &
                           - 1d0/We_avg_rs_vf(We_idx(i,1),We_idx(i,2))%sf(j,k,l) * &
                                                          kappa_avg(We_idx(i,2)) * &
                                                          ( alpha_R(We_idx(i,2))   &
                                                          - alpha_L(We_idx(i,2)) )
                  
                END DO
            END IF
            
            s_L = MIN(vel_L(dir_idx(1)) - c_L, vel_R(dir_idx(1)) - c_R) 
            s_R = MAX(vel_R(dir_idx(1)) + c_R, vel_L(dir_idx(1)) + c_L) 
            
            s_S = ( pres_R - pres_L - dpres_We + rho_L*vel_L(dir_idx(1))  * &
                                                (s_L - vel_L(dir_idx(1))) - &
                                                 rho_R*vel_R(dir_idx(1))  * &
                                                (s_R - vel_R(dir_idx(1))) ) &
                / ( rho_L*(s_L - vel_L(dir_idx(1))) - &
                    rho_R*(s_R - vel_R(dir_idx(1))) )
            denom = rho_L*(s_L - vel_L(dir_idx(1))) - rho_R*(s_R - vel_R(dir_idx(1))) 

            IF (tvd_wave_speeds) THEN
                lo_s_L = MIN(lo_vel_L(dir_idx(1)) - lo_c_L, vel_avg(dir_idx(1)) - c_avg)
                lo_s_R = MAX(lo_vel_R(dir_idx(1)) + lo_c_R, vel_avg(dir_idx(1)) + c_avg)
    
                lo_s_S = (lo_pres_R - lo_pres_L - dpres_We + lo_rho_L*lo_vel_L(dir_idx(1)) * &
                                    (lo_s_L - lo_vel_L(dir_idx(1))) - &
                                    lo_rho_R*lo_vel_R(dir_idx(1)) * &
                                    (lo_s_R - lo_vel_R(dir_idx(1)))) & 
                        / (lo_rho_L*(lo_s_L - lo_vel_L(dir_idx(1))) - &
                           lo_rho_R*(lo_s_R - lo_vel_R(dir_idx(1))))
            END IF

        END SUBROUTINE s_compute_direct_wave_speeds ! --------------------------
        
        
        !>  Estimation of the left, right and star region wave speeds
        !!      by the approximation of the pressures and velocity in the
        !!      star regions, see Toro (1999). The pressures and velocity
        !!      are approximated by using the primitive variables Riemann
        !!      solver (PVRS) and the wave speeds are then estimated from
        !!      those approximations using the exact wave relations.
        !!  @param j  First coordinate index
        !!  @param k Second coordinate index
        !!  @param l  Third coordinate index          
        SUBROUTINE s_compute_pressure_velocity_wave_speeds(j,k,l) ! ------------
            
            INTEGER, INTENT(IN) :: j,k,l
            
            ! Left and right pressures in the star region
            REAL(KIND(0d0)) :: pres_SL, pres_SR
            REAL(KIND(0d0)) :: lo_pres_SL, lo_pres_SR
            

            REAL(KIND(0d0)) :: dpres_We !< Capillary pressure
            
            ! Left and right shock Mach numbers
            REAL(KIND(0d0)) :: Ms_L, Ms_R
            REAL(KIND(0d0)) :: lo_Ms_L, lo_Ms_R
            
            INTEGER :: i !< Generic loop iterator
            
            
            dpres_We = 0d0
            
            IF (We_size > 0 .AND. We_wave_speeds .AND. We_riemann_flux) THEN
                DO i = 1, We_size
                  dpres_We = dpres_We &
                           - 1d0/We_avg_rs_vf(We_idx(i,1),We_idx(i,2))%sf(j,k,l) * &
                                                          kappa_avg(We_idx(i,2)) * &
                                                          ( alpha_R(We_idx(i,2))   &
                                                          - alpha_L(We_idx(i,2)) )
                  
                END DO
            END IF
            
            pres_SL = 5d-1*( pres_L + pres_R - dpres_We + rho_avg*c_avg * &
                                                    ( vel_L(dir_idx(1)) - &
                                                      vel_R(dir_idx(1)) ) )
            pres_SR = pres_SL + dpres_We
            
            Ms_L = MAX(1d0, SQRT(1d0 + ((5d-1 + gamma_L)/(1d0 + gamma_L))    * &
                                       (pres_SL/pres_L - 1d0)*pres_L         / &
                                       ((pres_L + pi_inf_L/(1d0 + gamma_L))) ) )
            Ms_R = MAX(1d0, SQRT(1d0 + ((5d-1 + gamma_R)/(1d0 + gamma_R))    * &
                                       (pres_SR/pres_R - 1d0)*pres_R         / &
                                       ((pres_R + pi_inf_R/(1d0 + gamma_R))) ) )
            
            s_L = vel_L(dir_idx(1)) - c_L*Ms_L
            s_R = vel_R(dir_idx(1)) + c_R*Ms_R
            
            s_S = 5d-1*( (vel_L(dir_idx(1)) + vel_R(dir_idx(1))) + &
                         (pres_L - pres_R + dpres_We)            / &
                         (rho_avg*c_avg)                           )
            
            IF (tvd_wave_speeds) THEN
                lo_pres_SL = 5d-1*(lo_pres_L + lo_pres_R - dpres_We + rho_avg*c_avg * &
                                        ( lo_vel_L(dir_idx(1)) - &
                                        lo_vel_R(dir_idx(1))))
                lo_pres_SR = lo_pres_SL + dpres_We
    
                lo_Ms_L = MAX(1d0, SQRT(1d0 + ((5d-1 + lo_gamma_L)/(1d0 + lo_gamma_L)) * &
                                (lo_pres_SL/lo_pres_L - 1d0)*lo_pres_L / &
                                ((lo_pres_L + lo_pi_inf_L/(1d0 + lo_gamma_L)))))
                lo_Ms_R = MAX(1d0, SQRT(1d0 + ((5d-1 + lo_gamma_R)/(1d0 + lo_gamma_R)) * &
                                (lo_pres_SR/lo_pres_R - 1d0)*lo_pres_R / &
                                ((lo_pres_R + lo_pi_inf_R/(1d0 + lo_gamma_R)))))
    
                lo_s_L = lo_vel_L(dir_idx(1)) - lo_c_L*lo_Ms_L
                lo_s_R = lo_vel_R(dir_idx(1)) + lo_c_R*lo_Ms_R
    
                lo_s_S = 5d-1*( (lo_vel_L(dir_idx(1)) + lo_vel_R(dir_idx(1))) + &
                        (lo_pres_L - lo_pres_R + dpres_We) / &
                        (rho_avg*c_avg))
            END IF
            
        END SUBROUTINE s_compute_pressure_velocity_wave_speeds ! ---------------
        
        
        
        
        !>  This subroutine computes the flux limiter function value
        !!      at the cell boundary
        !!  @param j  First coordinate index
        !!  @param k Second coordinate index
        !!  @param l  Third coordinate index     
        !!  @param flux_lim_func The flux limiter
        !!  @param norm_dir Directional splitting index
        SUBROUTINE s_compute_flux_limiter(j,k,l,flux_lim_func,norm_dir) ! ------


            INTEGER, INTENT(IN) :: j,k,l,norm_dir
            REAL(KIND(0d0)), INTENT(OUT) :: flux_lim_func
            REAL(KIND(0d0)) :: top, bottom, slope

            IF (q_prim_rs_vf(cont_idx%end+norm_dir)%sf(j,k,l) >= 0d0) THEN
                  top = q_prim_rs_vf(adv_idx%beg)%sf( j ,k,l) - &
                    q_prim_rs_vf(adv_idx%beg)%sf(j-1,k,l)
                   bottom = q_prim_rs_vf(adv_idx%beg)%sf(j+1,k,l) - &
                    q_prim_rs_vf(adv_idx%beg)%sf( j ,k,l)
            ELSE
                  top = q_prim_rs_vf(adv_idx%beg)%sf(j+2,k,l) - &
                    q_prim_rs_vf(adv_idx%beg)%sf(j+1,k,l)
                   bottom = q_prim_rs_vf(adv_idx%beg)%sf(j+1,k,l) - &
                    q_prim_rs_vf(adv_idx%beg)%sf( j ,k,l)
            END IF
        
            ! Limit the flux limiter to only be applied where the change in
            ! volume fraction is greater than machine precision so that 
            ! insignificant fluctuations do not trip the limiter
!           IF (ABS(top) < sgm_eps) top = 0d0
!           IF (ABS(bottom) < sgm_eps) bottom = 0d0
            IF (ABS(top) < 1d-8) top = 0d0
            IF (ABS(bottom) < 1d-8) bottom = 0d0

            ! If top = bottom, then cell boundary is in a smooth region of
            ! the flow and the high order flux should be used. Also ensures
            ! that areas of no change in volume fraction (0/0) use the high
            ! order flux
            IF (top == bottom) THEN
                slope = 1d0
            ELSE
                slope = (top*bottom)/MAX(bottom**2d0,sgm_eps)
            END IF

            ! Flux limiter function
            IF (flux_lim == 1) THEN ! MINMOD (MM)
                flux_lim_func = MAX(0d0,MIN(1d0,slope))
            ELSEIF (flux_lim == 2) THEN ! MUSCL (MC)
                flux_lim_func = MAX(0d0,MIN(2d0*slope,5d-1*(1d0+slope),2d0))
            ELSEIF (flux_lim == 3) THEN ! OSPRE (OP)
                flux_lim_func = (15d-1*(slope**2d0+slope))/(slope**2d0+slope+1d0)
            ELSEIF (flux_lim == 4) THEN ! SUPERBEE (SB)
                flux_lim_func = MAX(0d0,MIN(1d0,2d0*slope),MIN(slope,2d0))
            ELSEIF (flux_lim == 5) THEN ! SWEBY (SW) (beta = 1.5)
                flux_lim_func = MAX(0d0,MIN(15d-1*slope,1d0),MIN(slope,15d-1))
            ELSEIF (flux_lim == 6) THEN ! VAN ALBADA (VA)
                flux_lim_func = (slope**2d0+slope)/(slope**2d0+1d0)
            ELSEIF (flux_lim == 7) THEN ! VAN LEER (VL)
                flux_lim_func = (ABS(slope) + slope)/(1d0 + ABS(slope))
            END IF
        END SUBROUTINE s_compute_flux_limiter ! --------------------------------




        !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
        SUBROUTINE s_initialize_riemann_solvers_module() ! ---------------------
            
            
            ! Allocating the variables that will be utilized to formulate the
            ! left, right, and average states of the Riemann problem, as well
            ! the Riemann problem solution
            ALLOCATE(qL_prim_rs_vf(1:sys_size), qR_prim_rs_vf(1:sys_size))
            
            ALLOCATE(flux_rs_vf(1:sys_size), flux_src_rs_vf(1:sys_size))
            
            ALLOCATE(flux_gsrc_rs_vf(1:sys_size))
            
            ALLOCATE(vel_src_rs_vf(1:num_dims))
            
            ALLOCATE(alpha_rho_L(1:cont_idx%end), vel_L(1:num_dims))
            ALLOCATE(alpha_rho_R(1:cont_idx%end), vel_R(1:num_dims))

            ALLOCATE(vel_avg(1:num_dims))
           

            ALLOCATE(alpha_L(1:num_fluids))
            ALLOCATE(alpha_R(1:num_fluids))
            
            
            IF(We_size > 0) THEN
              
                ALLOCATE(alpha_avg_rs_vf(1:num_fluids))
                
                ALLOCATE(kappaL_rs_vf(1:num_fluids))
                ALLOCATE(kappaR_rs_vf(1:num_fluids))
                ALLOCATE(   kappa_avg(1:num_fluids))
                
                ALLOCATE(We_avg_rs_vf(1:num_fluids,1:num_fluids))
              
            END IF
            
            ALLOCATE(We_L(1:num_fluids,1:num_fluids))
            ALLOCATE(We_R(1:num_fluids,1:num_fluids))
            
            IF(ANY(Re_size > 0)) ALLOCATE(Re_avg_rs_vf(1:2))
            
            IF (tvd_riemann_flux) THEN
                ALLOCATE(q_prim_rs_vf(1:sys_size))
                ALLOCATE(lo_alpha_rho_L(1:cont_idx%end), lo_vel_L(1:num_dims))
                ALLOCATE(lo_alpha_rho_R(1:cont_idx%end), lo_vel_R(1:num_dims))
                ALLOCATE(lo_alpha_L(1:num_fluids))
                ALLOCATE(lo_alpha_R(1:num_fluids))
                ALLOCATE(lo_We_L(1:num_fluids, 1:num_fluids))
                ALLOCATE(lo_We_R(1:num_fluids, 1:num_fluids))
            END IF

            IF (riemann_solver == 3) THEN
                ALLOCATE(alpha_rho_IC(1:cont_idx%end), vel_IC(1:num_dims))
                ALLOCATE(alpha_IC(1:num_fluids))
            END IF
            
            ! Associating procedural pointer to the subroutine that will be
            ! utilized to calculate the solution of a given Riemann problem
            IF(riemann_solver == 1) THEN
                s_riemann_solver => s_hll_riemann_solver
            ELSEIF(riemann_solver == 2) THEN
                s_riemann_solver => s_hllc_riemann_solver
            ELSE
                s_riemann_solver => s_exact_riemann_solver
            END IF
           
            IF (bubbles) THEN
                ALLOCATE( R0_L(nb),R0_R(nb) )
                ALLOCATE( V0_L(nb), V0_R(nb))
                ALLOCATE( pbw_L(nb), pbw_R(nb) )
                IF (qbmm) THEN
                    ALLOCATE( moms_L(nb,nmom), moms_R(nb,nmom) )
                ELSE
                    IF (.NOT. polytropic) THEN
                        ALLOCATE( P0_L(nb), P0_R(nb) )
                    END IF 
                END IF 
            END IF
           
            IF (hypoelasticity) THEN
                ALLOCATE(tau_e_L(1:(num_dims*(num_dims+1)/2)) &
                       , tau_e_R(1:(num_dims*(num_dims+1)/2)))

                IF (tvd_riemann_flux) THEN
                    ALLOCATE(lo_tau_e_L(1:(num_dims*(num_dims+1)/2)) &
                           , lo_tau_e_R(1:(num_dims*(num_dims+1)/2)))
                END IF

                IF (riemann_solver == 3) THEN
                    ALLOCATE(tau_e_IC(1:(num_dims*(num_dims+1)/2)))
                END IF
            END IF

            ! Associating the procedural pointers to the procedures that will be
            ! utilized to compute the average state and estimate the wave speeds
            IF(riemann_solver /= 3) THEN
              
                IF(avg_state == 1) THEN
                    s_compute_average_state => s_compute_roe_average_state
                ELSE
                    s_compute_average_state => s_compute_arithmetic_average_state
                END IF
                
                IF(wave_speeds == 1) THEN
                    s_compute_wave_speeds => s_compute_direct_wave_speeds
                ELSE
                    s_compute_wave_speeds => s_compute_pressure_velocity_wave_speeds
                END IF
              
            END IF


            ! Associating procedural pointer to the subroutine that will be
            ! utilized to compute the viscous source flux
            IF (grid_geometry == 3) THEN
                s_compute_viscous_source_flux => s_compute_cylindrical_viscous_source_flux
            ELSE
                s_compute_viscous_source_flux => s_compute_cartesian_viscous_source_flux
            END IF
            
            
            ! Associating the procedural pointer to the appropriate subroutine
            ! that will be utilized in the conversion to the mixture variables
            
            IF (model_eqns == 1) THEN        ! Gamma/pi_inf model
                s_convert_to_mixture_variables => &
                             s_convert_mixture_to_mixture_variables
            ELSEIF (bubbles) THEN           ! Volume fraction for bubbles
                s_convert_to_mixture_variables => &
                             s_convert_species_to_mixture_variables_bubbles 
            ELSE                            ! Volume fraction model
                s_convert_to_mixture_variables => &
                             s_convert_species_to_mixture_variables
            END IF
            
            
        END SUBROUTINE s_initialize_riemann_solvers_module ! -------------------
        
        
        
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
        !!  @param kappaL_vf  Left averaged vol. frac. curvatures
        !!  @param kappaR_vf Right averaged vol. frac. curvatures
        !!  @param norm_dir Dir. splitting direction
        !!  @param ix Index bounds in the x-dir
        !!  @param iy Index bounds in the y-dir
        !!  @param iz Index bounds in the z-dir
        SUBROUTINE s_populate_riemann_states_variables_buffers( & ! ------------
                                    qL_prim_vf, dqL_prim_dx_vf, &
                                                dqL_prim_dy_vf, &
                                                dqL_prim_dz_vf, &
                                                  gm_alphaL_vf, &
                                                     kappaL_vf, &
                                    qR_prim_vf, dqR_prim_dx_vf, &
                                                dqR_prim_dy_vf, &
                                                dqR_prim_dz_vf, &
                                                  gm_alphaR_vf, &
                                                     kappaR_vf, &
                                            norm_dir, ix,iy,iz  )

            TYPE(scalar_field),  &
            DIMENSION(sys_size), &
            INTENT(INOUT) :: qL_prim_vf, qR_prim_vf
            
            TYPE(scalar_field),        &
            ALLOCATABLE, DIMENSION(:), &
            INTENT(INOUT) :: dqL_prim_dx_vf, dqR_prim_dx_vf, &
                             dqL_prim_dy_vf, dqR_prim_dy_vf, &
                             dqL_prim_dz_vf, dqR_prim_dz_vf, &
                               gm_alphaL_vf,   gm_alphaR_vf, &
                                  kappaL_vf,      kappaR_vf
            
            INTEGER, INTENT(IN) :: norm_dir
            
            TYPE(bounds_info), INTENT(IN) :: ix,iy,iz
            
            INTEGER :: i !< Generic loop iterator
            
            
            ! Population of Buffers in x-direction =============================
            IF(norm_dir == 1) THEN
              
                IF(bc_x%beg == -4) THEN    ! Riemann state extrap. BC at beginning
                
                    DO i = 1, sys_size
                        qL_prim_vf(i)%sf(-1,iy%beg:iy%end,iz%beg:iz%end) = &
                        qR_prim_vf(i)%sf( 0,iy%beg:iy%end,iz%beg:iz%end)
                    END DO
                
                    IF(ANY(Re_size > 0) .OR. hypoelasticity) THEN
                  
                        DO i = mom_idx%beg, mom_idx%end
                            dqL_prim_dx_vf(i)%sf(      -1      , &
                                                  iy%beg:iy%end, &
                                                  iz%beg:iz%end  ) = &
                            dqR_prim_dx_vf(i)%sf(       0      , &
                                                  iy%beg:iy%end, &
                                                  iz%beg:iz%end  )
                        END DO
                  
                        IF(n > 0) THEN
                    
                            DO i = mom_idx%beg, mom_idx%end
                                dqL_prim_dy_vf(i)%sf(      -1      , &
                                                      iy%beg:iy%end, &
                                                      iz%beg:iz%end  ) = &
                                dqR_prim_dy_vf(i)%sf(       0      , &
                                                      iy%beg:iy%end, &
                                                      iz%beg:iz%end  )
                            END DO
                    
                            IF(p > 0) THEN
                                DO i = mom_idx%beg, mom_idx%end
                                    dqL_prim_dz_vf(i)%sf(      -1      , &
                                                          iy%beg:iy%end, &
                                                          iz%beg:iz%end  ) = &
                                    dqR_prim_dz_vf(i)%sf(       0      , &
                                                          iy%beg:iy%end, &
                                                          iz%beg:iz%end  )
                                END DO
                            END IF
                    
                        END IF
                  
                    END IF
                
                    DO i = 1, crv_size
                        dqL_prim_dx_vf(E_idx+crv_idx(i))%sf(      -1      , &
                                                             iy%beg:iy%end, &
                                                             iz%beg:iz%end  ) = &
                        dqR_prim_dx_vf(E_idx+crv_idx(i))%sf(       0      , &
                                                             iy%beg:iy%end, &
                                                             iz%beg:iz%end  )
                          gm_alphaL_vf(   crv_idx(i)   )%sf(      -1      , &
                                                             iy%beg:iy%end, &
                                                             iz%beg:iz%end  ) = &
                          gm_alphaR_vf(   crv_idx(i)   )%sf(       0      , &
                                                             iy%beg:iy%end, &
                                                             iz%beg:iz%end  )
                             kappaL_vf(   crv_idx(i)   )%sf(      -1      , &
                                                             iy%beg:iy%end, &
                                                             iz%beg:iz%end  ) = &
                             kappaR_vf(   crv_idx(i)   )%sf(       0      , &
                                                             iy%beg:iy%end, &
                                                             iz%beg:iz%end  )
                        IF (n > 0) THEN
                            dqL_prim_dy_vf(E_idx+crv_idx(i))%sf(      -1      , &
                                                                 iy%beg:iy%end, &
                                                                 iz%beg:iz%end  ) = &
                            dqR_prim_dy_vf(E_idx+crv_idx(i))%sf(       0      , &
                                                                 iy%beg:iy%end, &
                                                                 iz%beg:iz%end  )
                        END IF
                    END DO
                
                    IF(p > 0) THEN
                        DO i = 1, crv_size
                            dqL_prim_dz_vf(E_idx+crv_idx(i))%sf(      -1      , &
                                                                 iy%beg:iy%end, &
                                                                 iz%beg:iz%end  ) = &
                            dqR_prim_dz_vf(E_idx+crv_idx(i))%sf(       0      , &
                                                                 iy%beg:iy%end, &
                                                                 iz%beg:iz%end  )
                        END DO
                    END IF
                
                END IF
              
                IF(bc_x%end == -4) THEN    ! Riemann state extrap. BC at end
                
                    DO i = 1, sys_size
                        qR_prim_vf(i)%sf(m+1,iy%beg:iy%end,iz%beg:iz%end) = &
                        qL_prim_vf(i)%sf( m ,iy%beg:iy%end,iz%beg:iz%end)
                    END DO
                
                    IF(ANY(Re_size > 0) .OR. hypoelasticity) THEN
                  
                        DO i = mom_idx%beg, mom_idx%end
                            dqR_prim_dx_vf(i)%sf(      m+1     , &
                                                  iy%beg:iy%end, &
                                                  iz%beg:iz%end  ) = &
                            dqL_prim_dx_vf(i)%sf(       m      , &
                                                  iy%beg:iy%end, &
                                                  iz%beg:iz%end  )
                        END DO
                  
                        IF(n > 0) THEN
                    
                            DO i = mom_idx%beg, mom_idx%end
                                dqR_prim_dy_vf(i)%sf(      m+1     , &
                                                      iy%beg:iy%end, &
                                                      iz%beg:iz%end  ) = &
                                dqL_prim_dy_vf(i)%sf(       m      , &
                                                      iy%beg:iy%end, &
                                                      iz%beg:iz%end  )
                            END DO
                    
                            IF(p > 0) THEN
                                DO i = mom_idx%beg, mom_idx%end
                                    dqR_prim_dz_vf(i)%sf(      m+1     , &
                                                          iy%beg:iy%end, &
                                                          iz%beg:iz%end  ) = &
                                    dqL_prim_dz_vf(i)%sf(       m      , &
                                                          iy%beg:iy%end, &
                                                          iz%beg:iz%end  )
                                END DO
                            END IF
                    
                        END IF
                  
                        DO i = 1, crv_size
                            dqR_prim_dx_vf(E_idx+crv_idx(i))%sf(      m+1     , &
                                                                 iy%beg:iy%end, &
                                                                 iz%beg:iz%end  ) = &
                            dqL_prim_dx_vf(E_idx+crv_idx(i))%sf(       m      , &
                                                                 iy%beg:iy%end, &
                                                                 iz%beg:iz%end  )
                              gm_alphaR_vf(   crv_idx(i)   )%sf(      m+1     , &
                                                                 iy%beg:iy%end, &
                                                                 iz%beg:iz%end  ) = &
                              gm_alphaL_vf(   crv_idx(i)   )%sf(       m      , &
                                                                 iy%beg:iy%end, &
                                                                 iz%beg:iz%end  )
                                 kappaR_vf(   crv_idx(i)   )%sf(      m+1     , &
                                                                 iy%beg:iy%end, &
                                                                 iz%beg:iz%end  ) = &
                                 kappaL_vf(   crv_idx(i)   )%sf(       m      , &
                                                                 iy%beg:iy%end, &
                                                                 iz%beg:iz%end  )
                            IF (n > 0) THEN
                                dqR_prim_dy_vf(E_idx+crv_idx(i))%sf(      m+1     , &
                                                                     iy%beg:iy%end, &
                                                                     iz%beg:iz%end  ) = &
                                dqL_prim_dy_vf(E_idx+crv_idx(i))%sf(       m      , &
                                                                     iy%beg:iy%end, &
                                                                     iz%beg:iz%end  )
                            END IF
                        END DO
                  
                        IF(p > 0) THEN
                            DO i = 1, crv_size
                                dqR_prim_dz_vf(E_idx+crv_idx(i))%sf(      m+1     , &
                                                                     iy%beg:iy%end, &
                                                                     iz%beg:iz%end  ) = &
                                dqL_prim_dz_vf(E_idx+crv_idx(i))%sf(       m      , &
                                                                     iy%beg:iy%end, &
                                                                     iz%beg:iz%end  )
                            END DO
                        END IF
                  
                    END IF
                
                END IF
            ! END: Population of Buffers in x-direction ========================
            
            
            ! Population of Buffers in y-direction =============================
            ELSEIF(norm_dir == 2) THEN
              
                IF(bc_y%beg == -4) THEN    ! Riemann state extrap. BC at beginning
                
                    DO i = 1, sys_size
                        qL_prim_vf(i)%sf(ix%beg:ix%end,-1,iz%beg:iz%end) = &
                        qR_prim_vf(i)%sf(ix%beg:ix%end, 0,iz%beg:iz%end)
                    END DO
                
                    IF(ANY(Re_size > 0)) THEN
                  
                        DO i = mom_idx%beg, mom_idx%end
                            dqL_prim_dx_vf(i)%sf( ix%beg:ix%end, &
                                                       -1      , &
                                                  iz%beg:iz%end  ) = &
                            dqR_prim_dx_vf(i)%sf( ix%beg:ix%end, &
                                                        0      , &
                                                  iz%beg:iz%end  )
                            IF (n > 0) THEN
                                dqL_prim_dy_vf(i)%sf( ix%beg:ix%end, &
                                                           -1      , &
                                                      iz%beg:iz%end  ) = &
                                dqR_prim_dy_vf(i)%sf( ix%beg:ix%end, &
                                                            0      , &
                                                      iz%beg:iz%end  )
                            END IF
                        END DO
                  
                        IF(p > 0) THEN
                            DO i = mom_idx%beg, mom_idx%end
                                dqL_prim_dz_vf(i)%sf( ix%beg:ix%end, &
                                                           -1      , &
                                                      iz%beg:iz%end  ) = &
                                dqR_prim_dz_vf(i)%sf( ix%beg:ix%end, &
                                                            0      , &
                                                      iz%beg:iz%end  )
                            END DO
                        END IF
                  
                    END IF
                
                    DO i = 1, crv_size
                        dqL_prim_dx_vf(E_idx+crv_idx(i))%sf( ix%beg:ix%end, &
                                                                  -1      , &
                                                             iz%beg:iz%end  ) = &
                        dqR_prim_dx_vf(E_idx+crv_idx(i))%sf( ix%beg:ix%end, &
                                                                   0      , &
                                                             iz%beg:iz%end  )
                          gm_alphaL_vf(   crv_idx(i)   )%sf( ix%beg:ix%end, &
                                                                  -1      , &
                                                             iz%beg:iz%end  ) = &
                          gm_alphaR_vf(   crv_idx(i)   )%sf( ix%beg:ix%end, &
                                                                   0      , &
                                                             iz%beg:iz%end  )
                             kappaL_vf(   crv_idx(i)   )%sf( ix%beg:ix%end, &
                                                                  -1      , &
                                                             iz%beg:iz%end  ) = &
                             kappaR_vf(   crv_idx(i)   )%sf( ix%beg:ix%end, &
                                                                   0      , &
                                                             iz%beg:iz%end  )
                        IF (n > 0) THEN
                            dqL_prim_dy_vf(E_idx+crv_idx(i))%sf( ix%beg:ix%end, &
                                                                      -1      , &
                                                                 iz%beg:iz%end  ) = &
                            dqR_prim_dy_vf(E_idx+crv_idx(i))%sf( ix%beg:ix%end, &
                                                                       0      , &
                                                                 iz%beg:iz%end  )
                        END IF
                    END DO
                    
                    IF(p > 0) THEN
                        DO i = 1, crv_size
                            dqL_prim_dz_vf(E_idx+crv_idx(i))%sf( ix%beg:ix%end, &
                                                                      -1      , &
                                                                 iz%beg:iz%end  ) = &
                            dqR_prim_dz_vf(E_idx+crv_idx(i))%sf( ix%beg:ix%end, &
                                                                       0      , &
                                                                 iz%beg:iz%end  )
                        END DO
                    END IF
                
                END IF
              
                IF(bc_y%end == -4) THEN    ! Riemann state extrap. BC at end
                
                    DO i = 1, sys_size
                        qR_prim_vf(i)%sf(ix%beg:ix%end,n+1,iz%beg:iz%end) = &
                        qL_prim_vf(i)%sf(ix%beg:ix%end, n ,iz%beg:iz%end)
                    END DO
                
                    IF(ANY(Re_size > 0)) THEN
                  
                        DO i = mom_idx%beg, mom_idx%end
                            dqR_prim_dx_vf(i)%sf( ix%beg:ix%end, &
                                                       n+1     , &
                                                  iz%beg:iz%end  ) = &
                            dqL_prim_dx_vf(i)%sf( ix%beg:ix%end, &
                                                        n      , &
                                                  iz%beg:iz%end  )
                            IF (n > 0) THEN
                                dqR_prim_dy_vf(i)%sf( ix%beg:ix%end, &
                                                           n+1     , &
                                                      iz%beg:iz%end  ) = &
                                dqL_prim_dy_vf(i)%sf( ix%beg:ix%end, &
                                                            n      , &
                                                      iz%beg:iz%end  )
                            END IF
                        END DO
                  
                        IF(p > 0) THEN
                            DO i = mom_idx%beg, mom_idx%end
                                dqR_prim_dz_vf(i)%sf( ix%beg:ix%end, &
                                                           n+1     , &
                                                      iz%beg:iz%end  ) = &
                                dqL_prim_dz_vf(i)%sf( ix%beg:ix%end, &
                                                            n      , &
                                                      iz%beg:iz%end  )
                            END DO
                        END IF
                  
                        DO i = 1, crv_size
                            dqR_prim_dx_vf(E_idx+crv_idx(i))%sf( ix%beg:ix%end, &
                                                                      n+1     , &
                                                                 iz%beg:iz%end  ) = &
                            dqL_prim_dx_vf(E_idx+crv_idx(i))%sf( ix%beg:ix%end, &
                                                                       n      , &
                                                                 iz%beg:iz%end  )
                              gm_alphaR_vf(   crv_idx(i)   )%sf( ix%beg:ix%end, &
                                                                      n+1     , &
                                                                 iz%beg:iz%end  ) = &
                              gm_alphaL_vf(   crv_idx(i)   )%sf( ix%beg:ix%end, &
                                                                       n      , &
                                                                 iz%beg:iz%end  )
                                 kappaR_vf(   crv_idx(i)   )%sf( ix%beg:ix%end, &
                                                                      n+1     , &
                                                                 iz%beg:iz%end  ) = &
                                 kappaL_vf(   crv_idx(i)   )%sf( ix%beg:ix%end, &
                                                                       n      , &
                                                                 iz%beg:iz%end  )
                            IF (n > 0) THEN
                                dqR_prim_dy_vf(E_idx+crv_idx(i))%sf( ix%beg:ix%end, &
                                                                          n+1     , &
                                                                     iz%beg:iz%end  ) = &
                                dqL_prim_dy_vf(E_idx+crv_idx(i))%sf( ix%beg:ix%end, &
                                                                           n      , &
                                                                     iz%beg:iz%end  )
                            END IF
                        END DO
                  
                        IF(p > 0) THEN
                            DO i = 1, crv_size
                                dqR_prim_dz_vf(E_idx+crv_idx(i))%sf( ix%beg:ix%end, &
                                                                    n+1     , &
                                                               iz%beg:iz%end  ) = &
                                dqL_prim_dz_vf(E_idx+crv_idx(i))%sf( ix%beg:ix%end, &
                                                                     n      , &
                                                               iz%beg:iz%end  )
                            END DO
                        END IF
                  
                    END IF
                
                END IF
            ! END: Population of Buffers in y-direction ========================
            
            
            ! Population of Buffers in z-direction =============================
            ELSE
              
                IF(bc_z%beg == -4) THEN    ! Riemann state extrap. BC at beginning
                
                    DO i = 1, sys_size
                        qL_prim_vf(i)%sf(ix%beg:ix%end,iy%beg:iy%end,-1) = &
                        qR_prim_vf(i)%sf(ix%beg:ix%end,iy%beg:iy%end, 0)
                    END DO
                
                    IF(ANY(Re_size > 0)) THEN
                        DO i = mom_idx%beg, mom_idx%end
                            dqL_prim_dx_vf(i)%sf( ix%beg:ix%end, &
                                                  iy%beg:iy%end, &
                                                       -1        ) = &
                            dqR_prim_dx_vf(i)%sf( ix%beg:ix%end, &
                                                  iy%beg:iy%end, &
                                                        0        )
                            dqL_prim_dy_vf(i)%sf( ix%beg:ix%end, &
                                                  iy%beg:iy%end, &
                                                       -1        ) = &
                            dqR_prim_dy_vf(i)%sf( ix%beg:ix%end, &
                                                  iy%beg:iy%end, &
                                                        0        )
                            dqL_prim_dz_vf(i)%sf( ix%beg:ix%end, &
                                                  iy%beg:iy%end, &
                                                       -1        ) = &
                            dqR_prim_dz_vf(i)%sf( ix%beg:ix%end, &
                                                  iy%beg:iy%end, &
                                                        0        )
                        END DO
                    END IF
                
                    DO i = 1, crv_size
                        dqL_prim_dx_vf(E_idx+crv_idx(i))%sf( ix%beg:ix%end, &
                                                             iy%beg:iy%end, &
                                                                  -1        ) = &
                        dqR_prim_dx_vf(E_idx+crv_idx(i))%sf( ix%beg:ix%end, &
                                                             iy%beg:iy%end, &
                                                                   0        )
                        dqL_prim_dy_vf(E_idx+crv_idx(i))%sf( ix%beg:ix%end, &
                                                             iy%beg:iy%end, &
                                                                  -1        ) = &
                        dqR_prim_dy_vf(E_idx+crv_idx(i))%sf( ix%beg:ix%end, &
                                                             iy%beg:iy%end, &
                                                                   0        )
                        dqL_prim_dz_vf(E_idx+crv_idx(i))%sf( ix%beg:ix%end, &
                                                             iy%beg:iy%end, &
                                                                  -1        ) = &
                        dqR_prim_dz_vf(E_idx+crv_idx(i))%sf( ix%beg:ix%end, &
                                                             iy%beg:iy%end, &
                                                                   0        )
                          gm_alphaL_vf(   crv_idx(i)   )%sf( ix%beg:ix%end, &
                                                             iy%beg:iy%end, &
                                                                  -1        ) = &
                          gm_alphaR_vf(   crv_idx(i)   )%sf( ix%beg:ix%end, &
                                                             iy%beg:iy%end, &
                                                                   0        )
                             kappaL_vf(   crv_idx(i)   )%sf( ix%beg:ix%end, &
                                                             iy%beg:iy%end, &
                                                                  -1        ) = &
                             kappaR_vf(   crv_idx(i)   )%sf( ix%beg:ix%end, &
                                                             iy%beg:iy%end, &
                                                                   0        )
                    END DO
                
                END IF
              
                IF(bc_z%end == -4) THEN    ! Riemann state extrap. BC at end
                
                    DO i = 1, sys_size
                        qR_prim_vf(i)%sf(ix%beg:ix%end,iy%beg:iy%end,p+1) = &
                        qL_prim_vf(i)%sf(ix%beg:ix%end,iy%beg:iy%end, p )
                    END DO
                
                    IF(ANY(Re_size > 0)) THEN
                        DO i = mom_idx%beg, mom_idx%end
                            dqR_prim_dx_vf(i)%sf( ix%beg:ix%end, &
                                                  iy%beg:iy%end, &
                                                       p+1       ) = &
                            dqL_prim_dx_vf(i)%sf( ix%beg:ix%end, &
                                                  iy%beg:iy%end, &
                                                        p        )
                            dqR_prim_dy_vf(i)%sf( ix%beg:ix%end, &
                                                  iy%beg:iy%end, &
                                                       p+1       ) = &
                            dqL_prim_dy_vf(i)%sf( ix%beg:ix%end, &
                                                  iy%beg:iy%end, &
                                                        p        )
                            dqR_prim_dz_vf(i)%sf( ix%beg:ix%end, &
                                                  iy%beg:iy%end, &
                                                       p+1       ) = &
                            dqL_prim_dz_vf(i)%sf( ix%beg:ix%end, &
                                                  iy%beg:iy%end, &
                                                        p        )
                        END DO
                    END IF
                
                    DO i = 1, crv_size
                        dqR_prim_dx_vf(E_idx+crv_idx(i))%sf( ix%beg:ix%end, &
                                                             iy%beg:iy%end, &
                                                                  p+1       ) = &
                        dqL_prim_dx_vf(E_idx+crv_idx(i))%sf( ix%beg:ix%end, &
                                                             iy%beg:iy%end, &
                                                                   p        )
                        dqR_prim_dy_vf(E_idx+crv_idx(i))%sf( ix%beg:ix%end, &
                                                             iy%beg:iy%end, &
                                                                  p+1       ) = &
                        dqL_prim_dy_vf(E_idx+crv_idx(i))%sf( ix%beg:ix%end, &
                                                             iy%beg:iy%end, &
                                                                   p        )
                        dqR_prim_dz_vf(E_idx+crv_idx(i))%sf( ix%beg:ix%end, &
                                                             iy%beg:iy%end, &
                                                                  p+1       ) = &
                        dqL_prim_dz_vf(E_idx+crv_idx(i))%sf( ix%beg:ix%end, &
                                                             iy%beg:iy%end, &
                                                                   p        )
                          gm_alphaR_vf(   crv_idx(i)   )%sf( ix%beg:ix%end, &
                                                             iy%beg:iy%end, &
                                                                  p+1       ) = &
                          gm_alphaL_vf(   crv_idx(i)   )%sf( ix%beg:ix%end, &
                                                             iy%beg:iy%end, &
                                                                   p        )
                             kappaR_vf(   crv_idx(i)   )%sf( ix%beg:ix%end, &
                                                             iy%beg:iy%end, &
                                                                  p+1       ) = &
                             kappaL_vf(   crv_idx(i)   )%sf( ix%beg:ix%end, &
                                                             iy%beg:iy%end, &
                                                                   p        )
                    END DO
                
                END IF
              
            END IF
            ! END: Population of Buffers in z-direction ========================
            
            
        END SUBROUTINE s_populate_riemann_states_variables_buffers ! -----------
        
 

        !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures needed to configure the chosen Riemann
        !!      solver algorithm.
        !!  @param qL_prim_vf The  left WENO-reconstructed cell-boundary values of the
        !!      cell-average primitive variables 
        !!  @param qR_prim_vf The right WENO-reconstructed cell-boundary values of the
        !!      cell-average primitive variables
        !!  @param kappaL_vf  Left averaged vol. frac. curvatures
        !!  @param kappaR_vf Right averaged vol. frac. curvatures
        !!  @param flux_vf Intra-cell fluxes
        !!  @param flux_src_vf Intra-cell fluxes sources
        !!  @param flux_gsrc_vf Intra-cell geometric fluxes sources
        !!  @param norm_dir Dir. splitting direction
        !!  @param ix Index bounds in the x-dir
        !!  @param iy Index bounds in the y-dir
        !!  @param iz Index bounds in the z-dir
        !!  @param q_prim_vf Cell-averaged primitive variables
        SUBROUTINE s_initialize_riemann_solver( qL_prim_vf, kappaL_vf, & ! -----
                                                qR_prim_vf, kappaR_vf, &
                                                            q_prim_vf, &
                                                 flux_vf, flux_src_vf, &
                                                         flux_gsrc_vf, &
                                                   norm_dir, ix,iy,iz  )

            TYPE(scalar_field),  &
            DIMENSION(sys_size), &
            INTENT(IN) :: qL_prim_vf, qR_prim_vf
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q_prim_vf
            
            TYPE(scalar_field),        &
            ALLOCATABLE, DIMENSION(:), &
            INTENT(IN) :: kappaL_vf, kappaR_vf
            
            TYPE(scalar_field),  &
            DIMENSION(sys_size), &
            INTENT(INOUT) :: flux_vf, flux_src_vf, flux_gsrc_vf
            
            INTEGER, INTENT(IN) :: norm_dir
            
            TYPE(bounds_info), INTENT(IN) :: ix,iy,iz
            

            INTEGER :: i,j,k ! Generic loop iterators
           
            INTEGER :: xbeg,xend,ybeg,yend,zbeg,zend
            INTEGER :: s1beg,s1end,s2beg,s2end,s3beg,s3end
            
            ! Configuring the coordinate direction indexes and flags
            IF(norm_dir == 1) THEN
               is1 = ix; is2 = iy; is3 = iz
               dir_idx = (/1,2,3/); dir_flg = (/1d0,0d0,0d0/)
            ELSEIF(norm_dir == 2) THEN
               is1 = iy; is2 = ix; is3 = iz
               dir_idx = (/2,1,3/); dir_flg = (/0d0,1d0,0d0/)
            ELSE
               is1 = iz; is2 = iy; is3 = ix
               dir_idx = (/3,1,2/); dir_flg = (/0d0,0d0,1d0/)
            END IF
           
            ! Extra indexes needed for hypoelasticity
            IF(hypoelasticity) THEN
                IF(norm_dir == 1) THEN
                    dir_idx_tau = (/1,2,4/)
                ELSEIF(norm_dir == 2) THEN
                    dir_idx_tau = (/2,3,5/)
                ELSE
                    dir_idx_tau = (/4,5,6/)
                END IF
            END IF
 
            ! Setting up special bounds for cell-average values
            xbeg = -buff_size; ybeg = 0; zbeg = 0
            IF (n > 0) ybeg = -buff_size; IF (p > 0) zbeg = -buff_size
            xend = m - xbeg; yend = n - ybeg; zend = p - zbeg

            ! Configuring the coordinate direction indexes
            IF (norm_dir == 1) THEN
                s1beg = xbeg; s1end = xend; s2beg = ybeg; s2end = yend; s3beg = zbeg; s3end = zend
            ELSEIF (norm_dir == 2) THEN
                s1beg = ybeg; s1end = yend; s2beg = xbeg; s2end = xend; s3beg = zbeg; s3end = zend
            ELSE
                s1beg = zbeg; s1end = zend; s2beg = ybeg; s2end = yend; s3beg = xbeg; s3end = xend
            END IF
            
            ! Allocating Left, Right and Average Riemann Problem States ========
            DO i = 1, sys_size
                ALLOCATE(qL_prim_rs_vf(i)%sf( is1%beg   : is1%end  , &
                                              is2%beg   : is2%end  , &
                                              is3%beg   : is3%end   ))
                ALLOCATE(qR_prim_rs_vf(i)%sf( is1%beg+1 : is1%end+1, &
                                              is2%beg   : is2%end  , &
                                              is3%beg   : is3%end   ))
                IF (tvd_riemann_flux) THEN
                    ALLOCATE(q_prim_rs_vf(i)%sf(    s1beg : s1end, &
                                                    s2beg : s2end, &
                                                    s3beg : s3end ))
                END IF
            END DO
            
            DO i = 1, crv_size
               ALLOCATE(   kappaL_rs_vf(crv_idx(i))%sf( is1%beg   : is1%end  , &
                                                        is2%beg   : is2%end  , &
                                                        is3%beg   : is3%end   ))
               ALLOCATE(   kappaR_rs_vf(crv_idx(i))%sf( is1%beg+1 : is1%end+1, &
                                                        is2%beg   : is2%end  , &
                                                        is3%beg   : is3%end   ))
               ALLOCATE(alpha_avg_rs_vf(crv_idx(i))%sf( is1%beg   : is1%end  , &
                                                        is2%beg   : is2%end  , &
                                                        is3%beg   : is3%end   ))
            END DO
            
            DO i = 1,2
                IF (Re_size(i) > 0) THEN
                    ALLOCATE(Re_avg_rs_vf(i)%sf( is1%beg : is1%end, &
                                                 is2%beg : is2%end, &
                                                 is3%beg : is3%end ))
                END IF
            END DO
            
            DO i = 1, We_size
                ALLOCATE(We_avg_rs_vf( We_idx(i,1), &
                                       We_idx(i,2)  )%sf( is1%beg : is1%end, &
                                                          is2%beg : is2%end, &
                                                          is3%beg : is3%end ))
            END DO
            ! ==================================================================
            
            
            ! Allocating Intercell Fluxes and Velocity =========================
            IF(norm_dir /= 1) THEN
               
                DO i = 1, sys_size
                    ALLOCATE(flux_rs_vf(i)%sf( is1%beg : is1%end, &
                                               is2%beg : is2%end, &
                                               is3%beg : is3%end ))
                    ALLOCATE(flux_gsrc_rs_vf(i)%sf( is1%beg : is1%end, &
                                                    is2%beg : is2%end, &
                                                    is3%beg : is3%end ))
                END DO
              
                ALLOCATE(flux_src_rs_vf(adv_idx%beg)%sf( is1%beg : is1%end, &
                                                         is2%beg : is2%end, &
                                                         is3%beg : is3%end ))
               
                    IF(riemann_solver == 1) THEN
                        DO i = adv_idx%beg+1, sys_size
                            ALLOCATE(flux_src_rs_vf(i)%sf( is1%beg : is1%end, &
                                                       is2%beg : is2%end, &
                                                       is3%beg : is3%end ))
                        END DO
                    END IF
                !ELSE
                !    DO i = adv_idx%beg,sys_size
                !        ALLOCATE(flux_src_rs_vf(i)%sf( is1%beg : is1%end, &
                !                                       is2%beg : is2%end, &
                !                                       is3%beg : is3%end ))
                !    END DO
                !END IF
            ELSE
               
                DO i = 1, sys_size
                    flux_rs_vf(i)%sf     => flux_vf(i)%sf
                    flux_src_rs_vf(i)%sf => flux_src_vf(i)%sf
                    flux_gsrc_rs_vf(i)%sf => flux_gsrc_vf(i)%sf
                END DO
               
            END IF
            
            IF(ANY(Re_size > 0) .OR. We_size > 0 .OR. hypoelasticity) THEN
                DO i = mom_idx%beg, E_idx
                    flux_src_vf(i)%sf = 0d0
                END DO
            END IF
            
            IF(riemann_solver == 1) THEN
                ALLOCATE(vel_src_rs_vf(dir_idx(1))%sf( is1%beg : is1%end, &
                                                       is2%beg : is2%end, &
                                                       is3%beg : is3%end ))
            ELSE
                vel_src_rs_vf(dir_idx(1))%sf => flux_src_rs_vf(adv_idx%beg)%sf
            END IF
            
            DO i = 2, num_dims
                ALLOCATE(vel_src_rs_vf(dir_idx(i))%sf( is1%beg : is1%end, &
                                                       is2%beg : is2%end, &
                                                       is3%beg : is3%end ))
            END DO
            ! END: Allocating Intercell Fluxes and Velocity ====================
           
                        
            ! Reshaping Inputted Data in x-direction ===========================
            IF(norm_dir == 1) THEN
               
                DO i = 1, sys_size
                    qL_prim_rs_vf(i)%sf = qL_prim_vf(i)%sf( ix%beg   : ix%end  , &
                                                            iy%beg   : iy%end  , &
                                                            iz%beg   : iz%end    )
                    qR_prim_rs_vf(i)%sf = qR_prim_vf(i)%sf( ix%beg+1 : ix%end+1, &
                                                            iy%beg   : iy%end  , &
                                                            iz%beg   : iz%end    )
                    IF (tvd_riemann_flux) THEN
                        q_prim_rs_vf(i)%sf = q_prim_vf(i)%sf(   xbeg : xend, &
                                                                ybeg : yend, &
                                                                zbeg : zend  )
                    END IF
                END DO
                DO i = 1, crv_size
                    kappaL_rs_vf(crv_idx(i))%sf = kappaL_vf(crv_idx(i))%sf( &
                                                       ix%beg   : ix%end  , &
                                                       iy%beg   : iy%end  , &
                                                       iz%beg   : iz%end    )
                    kappaR_rs_vf(crv_idx(i))%sf = kappaR_vf(crv_idx(i))%sf( &
                                                       ix%beg+1 : ix%end+1, &
                                                       iy%beg   : iy%end  , &
                                                       iz%beg   : iz%end    )
                END DO
            ! ==================================================================

            
            ! Reshaping Inputted Data in y-direction ===========================
            ELSEIF(norm_dir == 2) THEN
               
                DO i = 1, sys_size
                    DO j = ix%beg, ix%end
                        DO k = iy%beg, iy%end
                            qL_prim_rs_vf(i)%sf( k, j,:) = &
                                  qL_prim_vf(i)%sf(j, k ,iz%beg:iz%end)
                            qR_prim_rs_vf(i)%sf(k+1,j,:) = &
                                  qR_prim_vf(i)%sf(j,k+1,iz%beg:iz%end)
                        END DO
                    END DO
                    IF (tvd_riemann_flux) THEN
                        DO j = xbeg, xend
                            DO k = ybeg, yend
                                 q_prim_rs_vf(i)%sf(k,j,:) = &
                                    q_prim_vf(i)%sf(j,k,zbeg:zend)
                            END DO
                        END DO
                    END IF
                END DO
               
                DO i = 1, crv_size
                    DO j = ix%beg, ix%end
                        DO k = iy%beg, iy%end
                            kappaL_rs_vf(crv_idx(i))%sf( k ,j,:) = &
                              kappaL_vf(crv_idx(i))%sf(j, k ,iz%beg:iz%end)
                            kappaR_rs_vf(crv_idx(i))%sf(k+1,j,:) = &
                              kappaR_vf(crv_idx(i))%sf(j,k+1,iz%beg:iz%end)
                        END DO
                    END DO
                END DO
            ! ==================================================================
            
            
            ! Reshaping Inputted Data in z-direction ===========================
            ELSE
               
                DO i = 1, sys_size
                    DO j = ix%beg, ix%end
                        DO k = iz%beg, iz%end
                            qL_prim_rs_vf(i)%sf(k,:, j ) = &
                                  qL_prim_vf(i)%sf(j,iy%beg:iy%end, k )
                            qR_prim_rs_vf(i)%sf(k+1,:,j) = &
                                  qR_prim_vf(i)%sf(j,iy%beg:iy%end,k+1)
                        END DO
                    END DO
                    IF (tvd_riemann_flux) THEN
                        DO j = xbeg, xend
                            DO k = zbeg, zend
                                 q_prim_rs_vf(i)%sf(k,:,j) = &
                                    q_prim_vf(i)%sf(j,ybeg:yend,k)
                            END DO
                        END DO
                    END IF
                END DO
               
                DO i = 1, crv_size
                    DO j = ix%beg, ix%end
                        DO k = iz%beg, iz%end
                            kappaL_rs_vf(crv_idx(i))%sf( k ,:,j) = &
                              kappaL_vf(crv_idx(i))%sf(j,iy%beg:iy%end ,k )
                            kappaR_rs_vf(crv_idx(i))%sf(k+1,:,j) = &
                              kappaR_vf(crv_idx(i))%sf(j,iy%beg:iy%end,k+1)
                        END DO
                    END DO
                END DO
               
            END IF
            ! ==================================================================
            
            
        END SUBROUTINE s_initialize_riemann_solver ! ---------------------------
        
        
        
        
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
        SUBROUTINE s_compute_cylindrical_viscous_source_flux(     velL_vf, & ! -------------
                                                              dvelL_dx_vf, &
                                                              dvelL_dy_vf, &
                                                              dvelL_dz_vf, &
                                                                  velR_vf, &
                                                              dvelR_dx_vf, &
                                                              dvelR_dy_vf, &
                                                              dvelR_dz_vf, &
                                                              flux_src_vf, &
                                                                 norm_dir, &
                                                                 ix,iy,iz  )

            TYPE(scalar_field),  &
            DIMENSION(num_dims), &
            INTENT(IN) ::         velL_vf, velR_vf, &
                          dvelL_dx_vf, dvelR_dx_vf, &
                          dvelL_dy_vf, dvelR_dy_vf, &
                          dvelL_dz_vf, dvelR_dz_vf
            
            TYPE(scalar_field),  &
            DIMENSION(sys_size), &
            INTENT(INOUT) :: flux_src_vf
            
            INTEGER, INTENT(IN) :: norm_dir
            
            TYPE(bounds_info), INTENT(IN) :: ix,iy,iz
            
            ! Arithmetic mean of the left and right, WENO-reconstructed, cell-
            ! boundary values of cell-average first-order spatial derivatives
            ! of velocity
            REAL(KIND(0d0)), DIMENSION(num_dims) :: avg_vel
            REAL(KIND(0d0)), DIMENSION(num_dims) :: dvel_avg_dx
            REAL(KIND(0d0)), DIMENSION(num_dims) :: dvel_avg_dy
            REAL(KIND(0d0)), DIMENSION(num_dims) :: dvel_avg_dz
            
            ! Viscous stress tensor
            REAL(KIND(0d0)), DIMENSION(num_dims,num_dims) :: tau_Re
            
            ! Generic loop iterators
            INTEGER :: i,j,k,l
            
            
            ! Viscous Stresses in z-direction ==================================
            IF(norm_dir == 1) THEN
               
                IF(Re_size(1) > 0) THEN              ! Shear stresses
                    DO l = iz%beg, iz%end
                        DO k = iy%beg, iy%end
                            DO j = ix%beg, ix%end
                           
                                dvel_avg_dx(1) = 5d-1*( dvelL_dx_vf(1)%sf( j ,k,l) &
                                                      + dvelR_dx_vf(1)%sf(j+1,k,l) )
                           
                                tau_Re(1,1) = (4d0/3d0)*dvel_avg_dx(1) / &
                                              Re_avg_rs_vf(1)%sf(j,k,l)
                           
                                flux_src_vf(mom_idx%beg)%sf(j,k,l) = &
                                flux_src_vf(mom_idx%beg)%sf(j,k,l) - &
                                tau_Re(1,1)
                                
                                flux_src_vf(E_idx)%sf(j,k,l) = &
                                flux_src_vf(E_idx)%sf(j,k,l) - &
                                vel_src_rs_vf(1)%sf(j,k,l)   * &
                                tau_Re(1,1)
                           
                            END DO
                        END DO
                    END DO
                END IF
               
                IF(Re_size(2) > 0) THEN              ! Bulk stresses
                    DO l = iz%beg, iz%end
                        DO k = iy%beg, iy%end
                            DO j = ix%beg, ix%end
                           
                                dvel_avg_dx(1) = 5d-1*( dvelL_dx_vf(1)%sf( j, k,l) &
                                                      + dvelR_dx_vf(1)%sf(j+1,k,l) )
                           
                                tau_Re(1,1) = dvel_avg_dx(1) / &
                                              Re_avg_rs_vf(2)%sf(j,k,l)
                           
                                flux_src_vf(mom_idx%beg)%sf(j,k,l) = &
                                flux_src_vf(mom_idx%beg)%sf(j,k,l) - &
                                tau_Re(1,1)
                                
                                flux_src_vf(E_idx)%sf(j,k,l) = &
                                flux_src_vf(E_idx)%sf(j,k,l) - &
                                vel_src_rs_vf(1)%sf(j,k,l)   * &
                                tau_Re(1,1)
                           
                            END DO
                        END DO
                    END DO
                END IF
               
                IF(n == 0) RETURN
               
                IF(Re_size(1) > 0) THEN              ! Shear stresses
                    DO l = iz%beg, iz%end
                        DO k = iy%beg, iy%end
                            DO j = ix%beg, ix%end
                           
                                avg_vel(2) = 5d-1*( velL_vf(2)%sf( j ,k,l) &
                                                  + velR_vf(2)%sf(j+1,k,l) )
                           
                                DO i = 1,2
                                    dvel_avg_dy(i) = &
                                      5d-1*( dvelL_dy_vf(i)%sf( j ,k,l) &
                                           + dvelR_dy_vf(i)%sf(j+1,k,l) )
                                END DO

                                dvel_avg_dx(2) = 5d-1*( dvelL_dx_vf(2)%sf( j ,k,l) &
                                                      + dvelR_dx_vf(2)%sf(j+1,k,l) )

                                tau_Re(1,1) = -(2d0/3d0)*(dvel_avg_dy(2) + &
                                                avg_vel(2) / y_cc(k)) / &
                                                Re_avg_rs_vf(1)%sf(j,k,l)
                                
                                tau_Re(1,2) = (dvel_avg_dy(1) + dvel_avg_dx(2)) / &
                                              Re_avg_rs_vf(1)%sf(j,k,l)
                                
                                DO i = 1,2
                              
                                    flux_src_vf(cont_idx%end+i)%sf(j,k,l) = &
                                    flux_src_vf(cont_idx%end+i)%sf(j,k,l) - &
                                    tau_Re(1,i)
                                    
                                    flux_src_vf(E_idx)%sf(j,k,l) = &
                                    flux_src_vf(E_idx)%sf(j,k,l) - &
                                    vel_src_rs_vf(i)%sf(j,k,l)   * &
                                    tau_Re(1,i)
                              
                                END DO
                           
                            END DO
                        END DO
                    END DO
                END IF
               
                IF(Re_size(2) > 0) THEN              ! Bulk stresses
                    DO l = iz%beg, iz%end
                        DO k = iy%beg, iy%end
                            DO j = ix%beg, ix%end
                           
                                avg_vel(2) = 5d-1*( velL_vf(2)%sf( j ,k,l) &
                                                  + velR_vf(2)%sf(j+1,k,l) )
                           
                                dvel_avg_dy(2) = 5d-1*( dvelL_dy_vf(2)%sf( j, k,l) &
                                                      + dvelR_dy_vf(2)%sf(j+1,k,l) )
                                
                                tau_Re(1,1) = ( dvel_avg_dy(2) + &
                                               avg_vel(2) / y_cc(k)) / &
                                               Re_avg_rs_vf(2)%sf(j,k,l)
                                
                                flux_src_vf(mom_idx%beg)%sf(j,k,l) = &
                                flux_src_vf(mom_idx%beg)%sf(j,k,l) - &
                                tau_Re(1,1)
                                
                                flux_src_vf(E_idx)%sf(j,k,l) = &
                                flux_src_vf(E_idx)%sf(j,k,l) - &
                                vel_src_rs_vf(1)%sf(j,k,l)   * &
                                tau_Re(1,1)
                           
                            END DO
                        END DO
                    END DO
                END IF
               
                IF(p == 0) RETURN
               
                IF(Re_size(1) > 0) THEN              ! Shear stresses
                    DO l = iz%beg, iz%end
                        DO k = iy%beg, iy%end
                            DO j = ix%beg, ix%end
                           
                                DO i = 1,3,2
                                    dvel_avg_dz(i) = &
                                      5d-1*( dvelL_dz_vf(i)%sf( j ,k,l) &
                                           + dvelR_dz_vf(i)%sf(j+1,k,l) )
                                END DO
                           
                                dvel_avg_dx(3) = 5d-1*( dvelL_dx_vf(3)%sf( j ,k,l) &
                                                      + dvelR_dx_vf(3)%sf(j+1,k,l) )
                                
                                tau_Re(1,1) = -(2d0/3d0)*dvel_avg_dz(3) / y_cc(k) / &
                                               Re_avg_rs_vf(1)%sf(j,k,l)
                                
                                tau_Re(1,3) = (dvel_avg_dz(1)/y_cc(k) + dvel_avg_dx(3)) / &
                                              Re_avg_rs_vf(1)%sf(j,k,l)
                                
                                DO i = 1,3,2
                                   
                                    flux_src_vf(cont_idx%end+i)%sf(j,k,l) = &
                                    flux_src_vf(cont_idx%end+i)%sf(j,k,l) - &
                                    tau_Re(1,i)
                                    
                                    flux_src_vf(E_idx)%sf(j,k,l) = &
                                    flux_src_vf(E_idx)%sf(j,k,l) - &
                                    vel_src_rs_vf(i)%sf(j,k,l)   * &
                                    tau_Re(1,i)
                                   
                                END DO
                           
                            END DO
                        END DO
                    END DO
                END IF
               
                IF(Re_size(2) > 0) THEN              ! Bulk stresses
                    DO l = iz%beg, iz%end
                        DO k = iy%beg, iy%end
                            DO j = ix%beg, ix%end
                           
                                dvel_avg_dz(3) = 5d-1*( dvelL_dz_vf(3)%sf( j, k,l) &
                                                      + dvelR_dz_vf(3)%sf(j+1,k,l) )
                                
                                tau_Re(1,1) = dvel_avg_dz(3) / y_cc(k) / &
                                              Re_avg_rs_vf(2)%sf(j,k,l)

                                flux_src_vf(mom_idx%beg)%sf(j,k,l) = &
                                flux_src_vf(mom_idx%beg)%sf(j,k,l) - &
                                tau_Re(1,1)
                                
                                flux_src_vf(E_idx)%sf(j,k,l) = &
                                flux_src_vf(E_idx)%sf(j,k,l) - &
                                vel_src_rs_vf(1)%sf(j,k,l)   * &
                                tau_Re(1,1)
                           
                            END DO
                        END DO
                    END DO
                END IF
            ! END: Viscous Stresses in z-direction =============================
            
            
            ! Viscous Stresses in r-direction ==================================
            ELSEIF(norm_dir == 2) THEN
               
                IF(Re_size(1) > 0) THEN              ! Shear stresses
                    DO l = iz%beg, iz%end
                        DO k = iy%beg, iy%end
                            DO j = ix%beg, ix%end
                           
                                avg_vel(2) = 5d-1*( velL_vf(2)%sf(j, k ,l) &
                                                  + velR_vf(2)%sf(j,k+1,l) )
                           
                                DO i = 1,2
                              
                                    dvel_avg_dx(i) = &
                                            5d-1*( dvelL_dx_vf(i)%sf(j, k ,l) &
                                                 + dvelR_dx_vf(i)%sf(j,k+1,l) )
                                    
                                    dvel_avg_dy(i) = &
                                            5d-1*( dvelL_dy_vf(i)%sf(j, k ,l) &
                                                 + dvelR_dy_vf(i)%sf(j,k+1,l) )
                              
                                END DO
                           
                                tau_Re(2,1) = (dvel_avg_dy(1) + dvel_avg_dx(2)) / &
                                              Re_avg_rs_vf(1)%sf(k,j,l)

                                tau_Re(2,2) = ( 4d0*dvel_avg_dy(2)           &
                                              - 2d0*dvel_avg_dx(1)           &
                                              - 2d0*avg_vel(2) / y_cb(k) ) / &
                                              ( 3d0*Re_avg_rs_vf(1)%sf(k,j,l))
                           
                                DO i = 1,2
                              
                                    flux_src_vf(cont_idx%end+i)%sf(j,k,l) = &
                                    flux_src_vf(cont_idx%end+i)%sf(j,k,l) - &
                                    tau_Re(2,i)
                                    
                                    flux_src_vf(E_idx)%sf(j,k,l) = &
                                    flux_src_vf(E_idx)%sf(j,k,l) - &
                                    vel_src_rs_vf(i)%sf(k,j,l)   * &
                                    tau_Re(2,i)
                              
                                END DO
                           
                            END DO
                        END DO
                    END DO
                END IF
               
                IF(Re_size(2) > 0) THEN              ! Bulk stresses
                    DO l = iz%beg, iz%end
                        DO k = iy%beg, iy%end
                            DO j = ix%beg, ix%end
                           
                                avg_vel(2) = 5d-1*( velL_vf(2)%sf(j, k ,l) &
                                                  + velR_vf(2)%sf(j,k+1,l) )
                           
                                dvel_avg_dx(1) = 5d-1*( dvelL_dx_vf(1)%sf(j, k ,l) &
                                                      + dvelR_dx_vf(1)%sf(j,k+1,l) )
                           
                                dvel_avg_dy(2) = 5d-1*( dvelL_dy_vf(2)%sf(j, k ,l) &
                                                      + dvelR_dy_vf(2)%sf(j,k+1,l) )
                                
                                tau_Re(2,2) = (dvel_avg_dx(1) + dvel_avg_dy(2) + &
                                               avg_vel(2) / y_cb(k)) / &
                                              Re_avg_rs_vf(2)%sf(k,j,l)

                                flux_src_vf(mom_idx%beg+1)%sf(j,k,l) = &
                                flux_src_vf(mom_idx%beg+1)%sf(j,k,l) - &
                                tau_Re(2,2)
                                
                                flux_src_vf(E_idx)%sf(j,k,l) = &
                                flux_src_vf(E_idx)%sf(j,k,l) - &
                                vel_src_rs_vf(2)%sf(k,j,l)   * &
                                tau_Re(2,2)
                           
                            END DO
                        END DO
                    END DO
                END IF
               
                IF(p == 0) RETURN
               
                IF(Re_size(1) > 0) THEN              ! Shear stresses
                    DO l = iz%beg, iz%end
                        DO k = iy%beg, iy%end
                            DO j = ix%beg, ix%end
                           
                                avg_vel(3) = 5d-1*( velL_vf(3)%sf(j, k ,l) &
                                                  + velR_vf(3)%sf(j,k+1,l) )
                           
                                DO i = 2,3
                                    dvel_avg_dz(i) = &
                                      5d-1*( dvelL_dz_vf(i)%sf(j, k ,l) &
                                           + dvelR_dz_vf(i)%sf(j,k+1,l) )
                                END DO
                           
                                dvel_avg_dy(3) = 5d-1*( dvelL_dy_vf(3)%sf(j, k ,l) &
                                                      + dvelR_dy_vf(3)%sf(j,k+1,l) )
                                
                                tau_Re(2,2) = -(2d0/3d0)*dvel_avg_dz(3) / y_cb(k) / &
                                               Re_avg_rs_vf(1)%sf(k,j,l)
                                
                                tau_Re(2,3) = ((dvel_avg_dz(2) - avg_vel(3)) / &
                                              y_cb(k) + dvel_avg_dy(3)) / &
                                              Re_avg_rs_vf(1)%sf(k,j,l)
                                
                                DO i = 2,3
                                   
                                    flux_src_vf(cont_idx%end+i)%sf(j,k,l) = &
                                    flux_src_vf(cont_idx%end+i)%sf(j,k,l) - &
                                    tau_Re(2,i)
                                    
                                    flux_src_vf(E_idx)%sf(j,k,l) = &
                                    flux_src_vf(E_idx)%sf(j,k,l) - &
                                    vel_src_rs_vf(i)%sf(k,j,l)   * &
                                    tau_Re(2,i)
                                   
                                END DO
                           
                            END DO
                        END DO
                    END DO
                END IF
               
                IF(Re_size(2) > 0) THEN              ! Bulk stresses
                    DO l = iz%beg, iz%end
                        DO k = iy%beg, iy%end
                            DO j = ix%beg, ix%end
                           
                                dvel_avg_dz(3) = 5d-1*( dvelL_dz_vf(3)%sf(j, k ,l) &
                                                      + dvelR_dz_vf(3)%sf(j,k+1,l) )
                                
                                tau_Re(2,2) = dvel_avg_dz(3) / y_cb(k) / &
                                              Re_avg_rs_vf(2)%sf(k,j,l)
                                
                                flux_src_vf(mom_idx%beg+1)%sf(j,k,l) = &
                                flux_src_vf(mom_idx%beg+1)%sf(j,k,l) - &
                                tau_Re(2,2)
                                
                                flux_src_vf(E_idx)%sf(j,k,l) = &
                                flux_src_vf(E_idx)%sf(j,k,l) - &
                                vel_src_rs_vf(2)%sf(k,j,l)   * &
                                tau_Re(2,2)
                           
                            END DO
                        END DO
                    END DO
                END IF
            ! END: Viscous Stresses in r-direction =============================
            
            
            ! Viscous Stresses in theta-direction ==================================
            ELSE
               
                IF(Re_size(1) > 0) THEN              ! Shear stresses
                    DO l = iz%beg, iz%end
                        DO k = iy%beg, iy%end
                            DO j = ix%beg, ix%end
                           
                                DO i = 2,3
                                    avg_vel(i) = 5d-1*( velL_vf(i)%sf(j,k, l ) &
                                                      + velR_vf(i)%sf(j,k,l+1) )
                                END DO
                           
                                DO i = 1,3,2
                                    dvel_avg_dx(i) = &
                                      5d-1*( dvelL_dx_vf(i)%sf(j,k, l ) &
                                           + dvelR_dx_vf(i)%sf(j,k,l+1) )                 
                                END DO
                           
                                DO i = 2,3
                                    dvel_avg_dy(i) = &
                                      5d-1*( dvelL_dy_vf(i)%sf(j,k, l ) &
                                           + dvelR_dy_vf(i)%sf(j,k,l+1) )
                                END DO
                           
                                DO i = 1,3
                                    dvel_avg_dz(i) = &
                                      5d-1*( dvelL_dz_vf(i)%sf(j,k, l ) &
                                           + dvelR_dz_vf(i)%sf(j,k,l+1) )
                                END DO

                                tau_Re(3,1) = (dvel_avg_dz(1)/y_cc(k) + dvel_avg_dx(3)) / &
                                              Re_avg_rs_vf(1)%sf(l,k,j) / &
                                              y_cc(k)
                                
                                tau_Re(3,2) = ((dvel_avg_dz(2) - avg_vel(3)) / &
                                              y_cc(k) + dvel_avg_dy(3)) / &
                                              Re_avg_rs_vf(1)%sf(l,k,j) / &
                                              y_cc(k)
                                
                                tau_Re(3,3) = ( 4d0*dvel_avg_dz(3)/y_cc(k) &
                                              - 2d0*dvel_avg_dx(1)         &
                                              - 2d0*dvel_avg_dy(2)         &
                                              + 4d0*avg_vel(2)/y_cc(k) ) / &
                                              ( 3d0*Re_avg_rs_vf(1)%sf(l,k,j)) / &
                                              y_cc(k)

                                DO i = 1,3
                              
                                    flux_src_vf(cont_idx%end+i)%sf(j,k,l) = &
                                    flux_src_vf(cont_idx%end+i)%sf(j,k,l) - &
                                    tau_Re(3,i)
                                    
                                    flux_src_vf(E_idx)%sf(j,k,l) = &
                                    flux_src_vf(E_idx)%sf(j,k,l) - &
                                    vel_src_rs_vf(i)%sf(l,k,j)   * &
                                    tau_Re(3,i)
                              
                                END DO
                           
                            END DO
                        END DO
                    END DO
                END IF
               
                IF(Re_size(2) > 0) THEN              ! Bulk stresses
                    DO l = iz%beg, iz%end
                        DO k = iy%beg, iy%end
                            DO j = ix%beg, ix%end
                           
                                avg_vel(2) = 5d-1*( velL_vf(2)%sf(j,k, l ) &
                                                  + velR_vf(2)%sf(j,k,l+1) )
                                
                                dvel_avg_dx(1) = 5d-1*( dvelL_dx_vf(1)%sf(j,k, l ) &
                                                      + dvelR_dx_vf(1)%sf(j,k,l+1) )
                                
                                dvel_avg_dy(2) = 5d-1*( dvelL_dy_vf(2)%sf(j,k, l ) &
                                                      + dvelR_dy_vf(2)%sf(j,k,l+1) )
                                
                                dvel_avg_dz(3) = 5d-1*( dvelL_dz_vf(3)%sf(j,k, l ) &
                                                      + dvelR_dz_vf(3)%sf(j,k,l+1) )
                                
                                tau_Re(3,3) = ( dvel_avg_dx(1)           &
                                              + dvel_avg_dy(2)           &
                                              + dvel_avg_dz(3) / y_cc(k) &
                                              + avg_vel(2) / y_cc(k) ) / &
                                              Re_avg_rs_vf(2)%sf(l,k,j) / &
                                              y_cc(k)
                                
                                flux_src_vf(mom_idx%end)%sf(j,k,l) = &
                                flux_src_vf(mom_idx%end)%sf(j,k,l) - &
                                tau_Re(3,3)
                                
                                flux_src_vf(E_idx)%sf(j,k,l) = &
                                flux_src_vf(E_idx)%sf(j,k,l) - &
                                vel_src_rs_vf(3)%sf(l,k,j)   * &
                                tau_Re(3,3)
                           
                            END DO
                        END DO
                    END DO
                END IF
               
            END IF
            ! END: Viscous Stresses in theta-direction =============================
            
            
        END SUBROUTINE s_compute_cylindrical_viscous_source_flux ! -------------------------
        
       
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
        SUBROUTINE s_compute_cartesian_viscous_source_flux(     velL_vf, & ! -------------
                                                            dvelL_dx_vf, &
                                                            dvelL_dy_vf, &
                                                            dvelL_dz_vf, &
                                                                velR_vf, &
                                                            dvelR_dx_vf, &
                                                            dvelR_dy_vf, &
                                                            dvelR_dz_vf, &
                                                            flux_src_vf, &
                                                               norm_dir, &
                                                               ix,iy,iz  )

            TYPE(scalar_field),  &
            DIMENSION(num_dims), &
            INTENT(IN) ::         velL_vf, velR_vf, &
                          dvelL_dx_vf, dvelR_dx_vf, &
                          dvelL_dy_vf, dvelR_dy_vf, &
                          dvelL_dz_vf, dvelR_dz_vf
            
            TYPE(scalar_field),  &
            DIMENSION(sys_size), &
            INTENT(INOUT) :: flux_src_vf
            
            INTEGER, INTENT(IN) :: norm_dir
            
            TYPE(bounds_info), INTENT(IN) :: ix,iy,iz
            
            ! Arithmetic mean of the left and right, WENO-reconstructed, cell-
            ! boundary values of cell-average first-order spatial derivatives
            ! of velocity
            REAL(KIND(0d0)), DIMENSION(num_dims) :: dvel_avg_dx
            REAL(KIND(0d0)), DIMENSION(num_dims) :: dvel_avg_dy
            REAL(KIND(0d0)), DIMENSION(num_dims) :: dvel_avg_dz
            

            REAL(KIND(0d0)), DIMENSION(num_dims,num_dims) :: tau_Re !< Viscous stress tensor
            

            INTEGER :: i,j,k,l !< Generic loop iterators
            
            
            ! Viscous Stresses in x-direction ==================================
            IF(norm_dir == 1) THEN
               
                IF(Re_size(1) > 0) THEN              ! Shear stresses
                    DO l = iz%beg, iz%end
                        DO k = iy%beg, iy%end
                            DO j = ix%beg, ix%end
                           
                                dvel_avg_dx(1) = 5d-1*( dvelL_dx_vf(1)%sf( j ,k,l) &
                                                      + dvelR_dx_vf(1)%sf(j+1,k,l) )
                           
                                tau_Re(1,1) = (4d0/3d0)*dvel_avg_dx(1) / &
                                              Re_avg_rs_vf(1)%sf(j,k,l)
                           
                                flux_src_vf(mom_idx%beg)%sf(j,k,l) = &
                                flux_src_vf(mom_idx%beg)%sf(j,k,l) - &
                                tau_Re(1,1)
                                
                                flux_src_vf(E_idx)%sf(j,k,l) = &
                                flux_src_vf(E_idx)%sf(j,k,l) - &
                                vel_src_rs_vf(1)%sf(j,k,l)   * &
                                tau_Re(1,1)
                           
                            END DO
                        END DO
                    END DO
                END IF
               
                IF(Re_size(2) > 0) THEN              ! Bulk stresses
                    DO l = iz%beg, iz%end
                        DO k = iy%beg, iy%end
                            DO j = ix%beg, ix%end
                           
                                dvel_avg_dx(1) = 5d-1*( dvelL_dx_vf(1)%sf( j, k,l) &
                                                      + dvelR_dx_vf(1)%sf(j+1,k,l) )
                           
                                tau_Re(1,1) = dvel_avg_dx(1) / &
                                              Re_avg_rs_vf(2)%sf(j,k,l)
                           
                                flux_src_vf(mom_idx%beg)%sf(j,k,l) = &
                                flux_src_vf(mom_idx%beg)%sf(j,k,l) - &
                                tau_Re(1,1)
                                
                                flux_src_vf(E_idx)%sf(j,k,l) = &
                                flux_src_vf(E_idx)%sf(j,k,l) - &
                                vel_src_rs_vf(1)%sf(j,k,l)   * &
                                tau_Re(1,1)
                           
                            END DO
                        END DO
                    END DO
                END IF
               
                IF(n == 0) RETURN
               
                IF(Re_size(1) > 0) THEN              ! Shear stresses
                    DO l = iz%beg, iz%end
                        DO k = iy%beg, iy%end
                            DO j = ix%beg, ix%end
                           
                                DO i = 1,2
                                    dvel_avg_dy(i) = &
                                      5d-1*( dvelL_dy_vf(i)%sf( j ,k,l) &
                                           + dvelR_dy_vf(i)%sf(j+1,k,l) )
                                END DO
                           
                                dvel_avg_dx(2) = 5d-1*( dvelL_dx_vf(2)%sf( j ,k,l) &
                                                      + dvelR_dx_vf(2)%sf(j+1,k,l) )
                           
                                tau_Re(1,1) = -(2d0/3d0)*dvel_avg_dy(2) / &
                                               Re_avg_rs_vf(1)%sf(j,k,l)
                                
                                tau_Re(1,2) = (dvel_avg_dy(1) + dvel_avg_dx(2)) / &
                                              Re_avg_rs_vf(1)%sf(j,k,l)
                                
                                DO i = 1,2
                              
                                    flux_src_vf(cont_idx%end+i)%sf(j,k,l) = &
                                    flux_src_vf(cont_idx%end+i)%sf(j,k,l) - &
                                    tau_Re(1,i)
                                    
                                    flux_src_vf(E_idx)%sf(j,k,l) = &
                                    flux_src_vf(E_idx)%sf(j,k,l) - &
                                    vel_src_rs_vf(i)%sf(j,k,l)   * &
                                    tau_Re(1,i)
                              
                                END DO
                           
                            END DO
                        END DO
                    END DO
                END IF
               
                IF(Re_size(2) > 0) THEN              ! Bulk stresses
                    DO l = iz%beg, iz%end
                        DO k = iy%beg, iy%end
                            DO j = ix%beg, ix%end
                           
                                dvel_avg_dy(2) = 5d-1*( dvelL_dy_vf(2)%sf( j, k,l) &
                                                      + dvelR_dy_vf(2)%sf(j+1,k,l) )
                                
                                tau_Re(1,1) = dvel_avg_dy(2) / &
                                              Re_avg_rs_vf(2)%sf(j,k,l)
                                
                                flux_src_vf(mom_idx%beg)%sf(j,k,l) = &
                                flux_src_vf(mom_idx%beg)%sf(j,k,l) - &
                                tau_Re(1,1)
                                
                                flux_src_vf(E_idx)%sf(j,k,l) = &
                                flux_src_vf(E_idx)%sf(j,k,l) - &
                                vel_src_rs_vf(1)%sf(j,k,l)   * &
                                tau_Re(1,1)
                           
                            END DO
                        END DO
                    END DO
                END IF
               
                IF(p == 0) RETURN
               
                IF(Re_size(1) > 0) THEN              ! Shear stresses
                    DO l = iz%beg, iz%end
                        DO k = iy%beg, iy%end
                            DO j = ix%beg, ix%end
                           
                                DO i = 1,3,2
                                    dvel_avg_dz(i) = &
                                      5d-1*( dvelL_dz_vf(i)%sf( j ,k,l) &
                                           + dvelR_dz_vf(i)%sf(j+1,k,l) )
                                END DO
                           
                                dvel_avg_dx(3) = 5d-1*( dvelL_dx_vf(3)%sf( j ,k,l) &
                                                      + dvelR_dx_vf(3)%sf(j+1,k,l) )
                                
                                tau_Re(1,1) = -(2d0/3d0)*dvel_avg_dz(3) / &
                                               Re_avg_rs_vf(1)%sf(j,k,l)
                                
                                tau_Re(1,3) = (dvel_avg_dz(1) + dvel_avg_dx(3)) / &
                                              Re_avg_rs_vf(1)%sf(j,k,l)
                                
                                DO i = 1,3,2
                                   
                                    flux_src_vf(cont_idx%end+i)%sf(j,k,l) = &
                                    flux_src_vf(cont_idx%end+i)%sf(j,k,l) - &
                                    tau_Re(1,i)
                                    
                                    flux_src_vf(E_idx)%sf(j,k,l) = &
                                    flux_src_vf(E_idx)%sf(j,k,l) - &
                                    vel_src_rs_vf(i)%sf(j,k,l)   * &
                                    tau_Re(1,i)
                                   
                                END DO
                           
                            END DO
                        END DO
                    END DO
                END IF
               
                IF(Re_size(2) > 0) THEN              ! Bulk stresses
                    DO l = iz%beg, iz%end
                        DO k = iy%beg, iy%end
                            DO j = ix%beg, ix%end
                           
                                dvel_avg_dz(3) = 5d-1*( dvelL_dz_vf(3)%sf( j, k,l) &
                                                      + dvelR_dz_vf(3)%sf(j+1,k,l) )
                                
                                tau_Re(1,1) = dvel_avg_dz(3) / &
                                              Re_avg_rs_vf(2)%sf(j,k,l)
                                
                                flux_src_vf(mom_idx%beg)%sf(j,k,l) = &
                                flux_src_vf(mom_idx%beg)%sf(j,k,l) - &
                                tau_Re(1,1)
                                
                                flux_src_vf(E_idx)%sf(j,k,l) = &
                                flux_src_vf(E_idx)%sf(j,k,l) - &
                                vel_src_rs_vf(1)%sf(j,k,l)   * &
                                tau_Re(1,1)
                           
                            END DO
                        END DO
                    END DO
                END IF
            ! END: Viscous Stresses in x-direction =============================
            
            
            ! Viscous Stresses in y-direction ==================================
            ELSEIF(norm_dir == 2) THEN
               
                IF(Re_size(1) > 0) THEN              ! Shear stresses
                    DO l = iz%beg, iz%end
                        DO k = iy%beg, iy%end
                            DO j = ix%beg, ix%end
                           
                                DO i = 1,2
                              
                                    dvel_avg_dx(i) = &
                                            5d-1*( dvelL_dx_vf(i)%sf(j, k ,l) &
                                                 + dvelR_dx_vf(i)%sf(j,k+1,l) )
                                    
                                    dvel_avg_dy(i) = &
                                            5d-1*( dvelL_dy_vf(i)%sf(j, k ,l) &
                                                 + dvelR_dy_vf(i)%sf(j,k+1,l) )
                              
                                END DO
                           
                                tau_Re(2,1) = (dvel_avg_dy(1) + dvel_avg_dx(2)) / &
                                              Re_avg_rs_vf(1)%sf(k,j,l)
                                
                                tau_Re(2,2) = ( 4d0*dvel_avg_dy(2) &
                                              - 2d0*dvel_avg_dx(1) ) / &
                                              (3d0*Re_avg_rs_vf(1)%sf(k,j,l))
                           
                                DO i = 1,2
                              
                                    flux_src_vf(cont_idx%end+i)%sf(j,k,l) = &
                                    flux_src_vf(cont_idx%end+i)%sf(j,k,l) - &
                                    tau_Re(2,i)
                                    
                                    flux_src_vf(E_idx)%sf(j,k,l) = &
                                    flux_src_vf(E_idx)%sf(j,k,l) - &
                                    vel_src_rs_vf(i)%sf(k,j,l)   * &
                                    tau_Re(2,i)
                              
                                END DO
                           
                            END DO
                        END DO
                    END DO
                END IF
               
                IF(Re_size(2) > 0) THEN              ! Bulk stresses
                    DO l = iz%beg, iz%end
                        DO k = iy%beg, iy%end
                            DO j = ix%beg, ix%end
                           
                                dvel_avg_dx(1) = 5d-1*( dvelL_dx_vf(1)%sf(j, k ,l) &
                                                      + dvelR_dx_vf(1)%sf(j,k+1,l) )
                           
                                dvel_avg_dy(2) = 5d-1*( dvelL_dy_vf(2)%sf(j, k ,l) &
                                                      + dvelR_dy_vf(2)%sf(j,k+1,l) )
                                
                                tau_Re(2,2) = (dvel_avg_dx(1) + dvel_avg_dy(2)) / &
                                              Re_avg_rs_vf(2)%sf(k,j,l)
                                
                                flux_src_vf(mom_idx%beg+1)%sf(j,k,l) = &
                                flux_src_vf(mom_idx%beg+1)%sf(j,k,l) - &
                                tau_Re(2,2)
                                
                                flux_src_vf(E_idx)%sf(j,k,l) = &
                                flux_src_vf(E_idx)%sf(j,k,l) - &
                                vel_src_rs_vf(2)%sf(k,j,l)   * &
                                tau_Re(2,2)
                           
                            END DO
                        END DO
                    END DO
                END IF
               
                IF(p == 0) RETURN
               
                IF(Re_size(1) > 0) THEN              ! Shear stresses
                    DO l = iz%beg, iz%end
                        DO k = iy%beg, iy%end
                            DO j = ix%beg, ix%end
                           
                                DO i = 2,3
                                    dvel_avg_dz(i) = &
                                      5d-1*( dvelL_dz_vf(i)%sf(j, k ,l) &
                                           + dvelR_dz_vf(i)%sf(j,k+1,l) )
                                END DO
                           
                                dvel_avg_dy(3) = 5d-1*( dvelL_dy_vf(3)%sf(j, k ,l) &
                                                      + dvelR_dy_vf(3)%sf(j,k+1,l) )
                                
                                tau_Re(2,2) = -(2d0/3d0)*dvel_avg_dz(3) / &
                                               Re_avg_rs_vf(1)%sf(k,j,l)
                                
                                tau_Re(2,3) = (dvel_avg_dz(2) + dvel_avg_dy(3)) / &
                                              Re_avg_rs_vf(1)%sf(k,j,l)
                                
                                DO i = 2,3
                                   
                                    flux_src_vf(cont_idx%end+i)%sf(j,k,l) = &
                                    flux_src_vf(cont_idx%end+i)%sf(j,k,l) - &
                                    tau_Re(2,i)
                                    
                                    flux_src_vf(E_idx)%sf(j,k,l) = &
                                    flux_src_vf(E_idx)%sf(j,k,l) - &
                                    vel_src_rs_vf(i)%sf(k,j,l)   * &
                                    tau_Re(2,i)
                                   
                                END DO
                           
                            END DO
                        END DO
                    END DO
                END IF
               
                IF(Re_size(2) > 0) THEN              ! Bulk stresses
                    DO l = iz%beg, iz%end
                        DO k = iy%beg, iy%end
                            DO j = ix%beg, ix%end
                           
                                dvel_avg_dz(3) = 5d-1*( dvelL_dz_vf(3)%sf(j, k ,l) &
                                                      + dvelR_dz_vf(3)%sf(j,k+1,l) )
                                
                                tau_Re(2,2) = dvel_avg_dz(3) / &
                                              Re_avg_rs_vf(2)%sf(k,j,l)
                                
                                flux_src_vf(mom_idx%beg+1)%sf(j,k,l) = &
                                flux_src_vf(mom_idx%beg+1)%sf(j,k,l) - &
                                tau_Re(2,2)
                                
                                flux_src_vf(E_idx)%sf(j,k,l) = &
                                flux_src_vf(E_idx)%sf(j,k,l) - &
                                vel_src_rs_vf(2)%sf(k,j,l)   * &
                                tau_Re(2,2)
                           
                            END DO
                        END DO
                    END DO
                END IF
            ! END: Viscous Stresses in y-direction =============================
            
            
            ! Viscous Stresses in z-direction ==================================
            ELSE
               
                IF(Re_size(1) > 0) THEN              ! Shear stresses
                    DO l = iz%beg, iz%end
                        DO k = iy%beg, iy%end
                            DO j = ix%beg, ix%end
                           
                                DO i = 1,3,2
                                    dvel_avg_dx(i) = &
                                      5d-1*( dvelL_dx_vf(i)%sf(j,k, l ) &
                                           + dvelR_dx_vf(i)%sf(j,k,l+1) )                 
                                END DO
                           
                                DO i = 2,3
                                    dvel_avg_dy(i) = &
                                      5d-1*( dvelL_dy_vf(i)%sf(j,k, l ) &
                                           + dvelR_dy_vf(i)%sf(j,k,l+1) )
                                END DO
                           
                                DO i = 1,3
                                    dvel_avg_dz(i) = &
                                      5d-1*( dvelL_dz_vf(i)%sf(j,k, l ) &
                                           + dvelR_dz_vf(i)%sf(j,k,l+1) )
                                END DO
                           
                                tau_Re(3,1) = (dvel_avg_dz(1) + dvel_avg_dx(3)) / &
                                              Re_avg_rs_vf(1)%sf(l,k,j)
                                
                                tau_Re(3,2) = (dvel_avg_dz(2) + dvel_avg_dy(3)) / &
                                              Re_avg_rs_vf(1)%sf(l,k,j)
                                
                                tau_Re(3,3) = ( 4d0*dvel_avg_dz(3) &
                                              - 2d0*dvel_avg_dx(1) &
                                              - 2d0*dvel_avg_dy(2) ) / &
                                              (3d0*Re_avg_rs_vf(1)%sf(l,k,j))
                           
                                DO i = 1,3
                              
                                    flux_src_vf(cont_idx%end+i)%sf(j,k,l) = &
                                    flux_src_vf(cont_idx%end+i)%sf(j,k,l) - &
                                    tau_Re(3,i)
                                    
                                    flux_src_vf(E_idx)%sf(j,k,l) = &
                                    flux_src_vf(E_idx)%sf(j,k,l) - &
                                    vel_src_rs_vf(i)%sf(l,k,j)   * &
                                    tau_Re(3,i)
                              
                                END DO
                           
                            END DO
                        END DO
                    END DO
                END IF
               
                IF(Re_size(2) > 0) THEN              ! Bulk stresses
                    DO l = iz%beg, iz%end
                        DO k = iy%beg, iy%end
                            DO j = ix%beg, ix%end
                           
                                dvel_avg_dx(1) = 5d-1*( dvelL_dx_vf(1)%sf(j,k, l ) &
                                                      + dvelR_dx_vf(1)%sf(j,k,l+1) )
                                
                                dvel_avg_dy(2) = 5d-1*( dvelL_dy_vf(2)%sf(j,k, l ) &
                                                      + dvelR_dy_vf(2)%sf(j,k,l+1) )
                                
                                dvel_avg_dz(3) = 5d-1*( dvelL_dz_vf(3)%sf(j,k, l ) &
                                                      + dvelR_dz_vf(3)%sf(j,k,l+1) )
                                
                                tau_Re(3,3) = ( dvel_avg_dx(1) &
                                              + dvel_avg_dy(2) &
                                              + dvel_avg_dz(3) ) / &
                                              Re_avg_rs_vf(2)%sf(l,k,j)
                                
                                flux_src_vf(mom_idx%end)%sf(j,k,l) = &
                                flux_src_vf(mom_idx%end)%sf(j,k,l) - &
                                tau_Re(3,3)
                                
                                flux_src_vf(E_idx)%sf(j,k,l) = &
                                flux_src_vf(E_idx)%sf(j,k,l) - &
                                vel_src_rs_vf(3)%sf(l,k,j)   * &
                                tau_Re(3,3)
                           
                            END DO
                        END DO
                    END DO
                END IF
               
            END IF
            ! END: Viscous Stresses in z-direction =============================
            
            
        END SUBROUTINE s_compute_cartesian_viscous_source_flux ! -------------------------
        
        
        !>  The purpose of the subroutine is to evaluate and account
        !!      for the contribution of capillary stresses in the source
        !!      flux for the momentum and energy.
        !! @param dalphaL_dx_vf  Left WENO-reconstructed cell-bndry value of cell-avg x-dir derivative of vol. frac.       
        !! @param dalphaL_dy_vf  Left WENO-reconstructed cell-bndry value of cell-avg y-dir derivative of vol. frac.       
        !! @param dalphaL_dz_vf  Left WENO-reconstructed cell-bndry value of cell-avg z-dir derivative of vol. frac.       
        !! @param dalphaR_dx_vf Right WENO-reconstructed cell-bndry value of cell-avg x-dir derivative of vol. frac.       
        !! @param dalphaR_dy_vf Right WENO-reconstructed cell-bndry value of cell-avg y-dir derivative of vol. frac.       
        !! @param dalphaR_dz_vf Right WENO-reconstructed cell-bndry value of cell-avg z-dir derivative of vol. frac.       
        !! @param gm_alphaL_vf   Left WENO-reconstructed cell-bndry value of the gradient magnitude
        !! @param gm_alphaR_vf  Right WENO-reconstructed cell-bndry value of the gradient magnitude
        !!  @param flux_src_vf Intercell flux
        !!  @param norm_dir Dimensional splitting coordinate direction
        !!  @param ix Index bounds in  first coordinate direction
        !!  @param iy Index bounds in second coordinate direction
        !!  @param iz Index bounds in  third coordinate direction
        SUBROUTINE s_compute_capillary_source_flux( dalphaL_dx_vf, & ! ---------
                                                    dalphaL_dy_vf, &
                                                    dalphaL_dz_vf, &
                                                     gm_alphaL_vf, &
                                                    dalphaR_dx_vf, &
                                                    dalphaR_dy_vf, &
                                                    dalphaR_dz_vf, &
                                                     gm_alphaR_vf, &
                                                      flux_src_vf, &
                                                         norm_dir, &
                                                         ix,iy,iz  )
           
            TYPE(scalar_field),    &
            DIMENSION(num_fluids), &
            INTENT(IN) :: dalphaL_dx_vf, dalphaR_dx_vf, &
                          dalphaL_dy_vf, dalphaR_dy_vf, &
                          dalphaL_dz_vf, dalphaR_dz_vf, &
                           gm_alphaL_vf,  gm_alphaR_vf
            
            TYPE(scalar_field),  &
            DIMENSION(sys_size), &
            INTENT(INOUT) :: flux_src_vf
            
            INTEGER, INTENT(IN) :: norm_dir
            
            TYPE(bounds_info), INTENT(IN) :: ix,iy,iz
            
            ! Arithmetic average of the left and right, WENO-reconstructed,
            ! cell-boundary values of the cell-average first-order spatial
            ! derivative of the volume fractions and gradient magnitude of
            ! the volume fractions
            REAL(KIND(0d0)), DIMENSION(num_fluids) :: dalpha_avg_dx
            REAL(KIND(0d0)), DIMENSION(num_fluids) :: dalpha_avg_dy
            REAL(KIND(0d0)), DIMENSION(num_fluids) :: dalpha_avg_dz
            REAL(KIND(0d0)), DIMENSION(num_fluids) ::  gm_alpha_avg
            
            ! Capillary stress tensor
            REAL(KIND(0d0)), DIMENSION(num_dims,num_dims) :: tau_We
            
            ! Generic loop iterators
            INTEGER :: i,j,k,l,r
            
            
            ! Capillary Stresses in x-direction ================================
            IF(norm_dir == 1) THEN
               
                DO i = 1, We_size
                    DO r = iz%beg, iz%end
                        DO l = iy%beg, iy%end
                            DO k = ix%beg, ix%end
                           
                                DO j = 1,2
                              
                                    dalpha_avg_dx(We_idx(i,j)) = &
                                      5d-1*( dalphaL_dx_vf(We_idx(i,j))%sf( k ,l,r) &
                                           + dalphaR_dx_vf(We_idx(i,j))%sf(k+1,l,r) )
                                    
                                    dalpha_avg_dy(We_idx(i,j)) = &
                                      5d-1*( dalphaL_dy_vf(We_idx(i,j))%sf( k ,l,r) &
                                           + dalphaR_dy_vf(We_idx(i,j))%sf(k+1,l,r) )
                                    
                                    gm_alpha_avg(We_idx(i,j)) = &
                                      5d-1*( gm_alphaL_vf(We_idx(i,j))%sf( k ,l,r) &
                                           + gm_alphaR_vf(We_idx(i,j))%sf(k+1,l,r) )
                              
                                END DO
                          
                                tau_We(1,1)= &
                                                  1d0/We_avg_rs_vf( We_idx(i,1), &
                                                      We_idx(i,2)  )%sf(k,l,r) * &
                                      ( alpha_avg_rs_vf(We_idx(i,1))%sf(k,l,r) * &
                                               (     gm_alpha_avg(We_idx(i,2)) - &
                                                    dalpha_avg_dx(We_idx(i,2)) * &
                                                    dalpha_avg_dx(We_idx(i,2)) / &
                                                 MAX(gm_alpha_avg(We_idx(i,2)) , &
                                                                      sgm_eps) ) &
                                      + alpha_avg_rs_vf(We_idx(i,2))%sf(k,l,r) * &
                                               (     gm_alpha_avg(We_idx(i,1)) - &
                                                    dalpha_avg_dx(We_idx(i,1)) * &
                                                    dalpha_avg_dx(We_idx(i,1)) / &
                                                 MAX(gm_alpha_avg(We_idx(i,1)) , &
                                                                      sgm_eps) ) )
                                
                                tau_We(1,2) = &
                                                - 1d0/We_avg_rs_vf( We_idx(i,1), &
                                                      We_idx(i,2)  )%sf(k,l,r) * &
                                      ( alpha_avg_rs_vf(We_idx(i,1))%sf(k,l,r) * &
                                               (    dalpha_avg_dx(We_idx(i,2)) * &
                                                    dalpha_avg_dy(We_idx(i,2)) / &
                                                 MAX(gm_alpha_avg(We_idx(i,2)) , &
                                                                      sgm_eps) ) &
                                      + alpha_avg_rs_vf(We_idx(i,2))%sf(k,l,r) * &
                                               (    dalpha_avg_dx(We_idx(i,1)) * &
                                                    dalpha_avg_dy(We_idx(i,1)) / &
                                                 MAX(gm_alpha_avg(We_idx(i,1)) , &
                                                                      sgm_eps) ) )

                                DO j = 1,2
                                   flux_src_vf(cont_idx%end+j)%sf(k,l,r) = &
                                   flux_src_vf(cont_idx%end+j)%sf(k,l,r) - &
                                   tau_We(1,j)
                                END DO
                           
                                tau_We(1,1)= &
                                               - 1d0/We_avg_rs_vf( We_idx(i,1), &
                                                     We_idx(i,2)  )%sf(k,l,r) * &
                                     ( alpha_avg_rs_vf(We_idx(i,1))%sf(k,l,r) * &
                                              (    dalpha_avg_dx(We_idx(i,2)) * &
                                                   dalpha_avg_dx(We_idx(i,2)) / &
                                                MAX(gm_alpha_avg(We_idx(i,2)) , &
                                                                     sgm_eps) ) &
                                     + alpha_avg_rs_vf(We_idx(i,2))%sf(k,l,r) * &
                                              (    dalpha_avg_dx(We_idx(i,1)) * &
                                                   dalpha_avg_dx(We_idx(i,1)) / &
                                                MAX(gm_alpha_avg(We_idx(i,1)) , &
                                                                     sgm_eps) ) )
                           
                                DO j = 1,2
                                    flux_src_vf(E_idx)%sf(k,l,r) = &
                                    flux_src_vf(E_idx)%sf(k,l,r) - &
                                    vel_src_rs_vf(j)%sf(k,l,r)   * &
                                    tau_We(1,j)
                                END DO
                           
                            END DO
                        END DO
                    END DO
                END DO
               
                IF(p == 0) RETURN
               
                DO i = 1, We_size
                    DO r = iz%beg, iz%end
                        DO l = iy%beg, iy%end
                            DO k = ix%beg, ix%end
                           
                                DO j = 1,2
                              
                                    dalpha_avg_dx(We_idx(i,j)) = &
                                      5d-1*( dalphaL_dx_vf(We_idx(i,j))%sf( k ,l,r) &
                                           + dalphaR_dx_vf(We_idx(i,j))%sf(k+1,l,r) )
                                    
                                    dalpha_avg_dz(We_idx(i,j)) = &
                                      5d-1*( dalphaL_dz_vf(We_idx(i,j))%sf( k ,l,r) &
                                           + dalphaR_dz_vf(We_idx(i,j))%sf(k+1,l,r) )
                                    
                                    gm_alpha_avg(We_idx(i,j)) = &
                                      5d-1*( gm_alphaL_vf(We_idx(i,j))%sf( k ,l,r) &
                                           + gm_alphaR_vf(We_idx(i,j))%sf(k+1,l,r) )
                              
                                END DO
                           
                                tau_We(1,3) = &
                                                - 1d0/We_avg_rs_vf( We_idx(i,1), &
                                                      We_idx(i,2)  )%sf(k,l,r) * &
                                      ( alpha_avg_rs_vf(We_idx(i,1))%sf(k,l,r) * &
                                               (    dalpha_avg_dx(We_idx(i,2)) * &
                                                    dalpha_avg_dz(We_idx(i,2)) / &
                                                 MAX(gm_alpha_avg(We_idx(i,2)) , &
                                                                      sgm_eps) ) &
                                      + alpha_avg_rs_vf(We_idx(i,2))%sf(k,l,r) * &
                                               (    dalpha_avg_dx(We_idx(i,1)) * &
                                                    dalpha_avg_dz(We_idx(i,1)) / &
                                                 MAX(gm_alpha_avg(We_idx(i,1)) , &
                                                                      sgm_eps) ) )
                                
                                flux_src_vf(mom_idx%end)%sf(k,l,r) = &
                                flux_src_vf(mom_idx%end)%sf(k,l,r) - &
                                tau_We(1,3)
                                
                                flux_src_vf(E_idx)%sf(k,l,r) = &
                                flux_src_vf(E_idx)%sf(k,l,r) - &
                                vel_src_rs_vf(3)%sf(k,l,r)   * &
                                tau_We(1,3)
                           
                            END DO
                        END DO
                    END DO
                END DO
            ! END: Capillary Stresses in x-direction ===========================
            
            
            ! Capillary Stresses in y-direction ================================
            ELSEIF(norm_dir == 2) THEN
               
                DO i = 1, We_size
                    DO r = iz%beg, iz%end
                        DO l = iy%beg, iy%end
                            DO k = ix%beg, ix%end
                           
                                DO j = 1,2
                              
                                    dalpha_avg_dx(We_idx(i,j)) = &
                                      5d-1*( dalphaL_dx_vf(We_idx(i,j))%sf(k, l ,r) &
                                           + dalphaR_dx_vf(We_idx(i,j))%sf(k,l+1,r) )
                                    
                                    dalpha_avg_dy(We_idx(i,j)) = &
                                      5d-1*( dalphaL_dy_vf(We_idx(i,j))%sf(k, l ,r) &
                                           + dalphaR_dy_vf(We_idx(i,j))%sf(k,l+1,r) )
                                    
                                    gm_alpha_avg(We_idx(i,j)) = &
                                      5d-1*( gm_alphaL_vf(We_idx(i,j))%sf(k, l ,r) &
                                           + gm_alphaR_vf(We_idx(i,j))%sf(k,l+1,r) )
                              
                                END DO
                           
                                tau_We(2,1) = &
                                                - 1d0/We_avg_rs_vf( We_idx(i,1), &
                                                      We_idx(i,2)  )%sf(l,k,r) * &
                                      ( alpha_avg_rs_vf(We_idx(i,1))%sf(l,k,r) * &
                                               (    dalpha_avg_dy(We_idx(i,2)) * &
                                                    dalpha_avg_dx(We_idx(i,2)) / &
                                                 MAX(gm_alpha_avg(We_idx(i,2)) , &
                                                                      sgm_eps) ) &
                                      + alpha_avg_rs_vf(We_idx(i,2))%sf(l,k,r) * &
                                               (    dalpha_avg_dy(We_idx(i,1)) * &
                                                    dalpha_avg_dx(We_idx(i,1)) / &
                                                 MAX(gm_alpha_avg(We_idx(i,1)) , &
                                                                      sgm_eps) ) )
                                
                                tau_We(2,2) = &
                                                  1d0/We_avg_rs_vf( We_idx(i,1), &
                                                      We_idx(i,2)  )%sf(l,k,r) * &
                                      ( alpha_avg_rs_vf(We_idx(i,1))%sf(l,k,r) * &
                                               (     gm_alpha_avg(We_idx(i,2)) - &
                                                    dalpha_avg_dy(We_idx(i,2)) * &
                                                    dalpha_avg_dy(We_idx(i,2)) / &
                                                 MAX(gm_alpha_avg(We_idx(i,2)) , &
                                                                      sgm_eps) ) &
                                      + alpha_avg_rs_vf(We_idx(i,2))%sf(l,k,r) * &
                                               (     gm_alpha_avg(We_idx(i,1)) - &
                                                    dalpha_avg_dy(We_idx(i,1)) * &
                                                    dalpha_avg_dy(We_idx(i,1)) / &
                                                 MAX(gm_alpha_avg(We_idx(i,1)) , &
                                                                      sgm_eps) ) )

                                DO j = 1,2
                                    flux_src_vf(cont_idx%end+j)%sf(k,l,r) = &
                                    flux_src_vf(cont_idx%end+j)%sf(k,l,r) - &
                                    tau_We(2,j)
                                END DO
                           
                                tau_We(2,2) = &
                                               - 1d0/We_avg_rs_vf( We_idx(i,1), &
                                                     We_idx(i,2)  )%sf(l,k,r) * &
                                     ( alpha_avg_rs_vf(We_idx(i,1))%sf(l,k,r) * &
                                              (    dalpha_avg_dy(We_idx(i,2)) * &
                                                   dalpha_avg_dy(We_idx(i,2)) / &
                                                MAX(gm_alpha_avg(We_idx(i,2)) , &
                                                                     sgm_eps) ) &
                                     + alpha_avg_rs_vf(We_idx(i,2))%sf(l,k,r) * &
                                              (    dalpha_avg_dy(We_idx(i,1)) * &
                                                   dalpha_avg_dy(We_idx(i,1)) / &
                                                MAX(gm_alpha_avg(We_idx(i,1)) , &
                                                                     sgm_eps) ) )
                           
                                DO j = 1,2
                                    flux_src_vf(E_idx)%sf(k,l,r) = &
                                    flux_src_vf(E_idx)%sf(k,l,r) - &
                                    vel_src_rs_vf(j)%sf(l,k,r)   * &
                                    tau_We(2,j)
                                END DO
                           
                            END DO
                        END DO
                    END DO
                END DO
               
                IF(p == 0) RETURN
               
                DO i = 1, We_size
                    DO r = iz%beg, iz%end
                        DO l = iy%beg, iy%end
                            DO k = ix%beg, ix%end
                           
                                DO j = 1,2
                              
                                    dalpha_avg_dy(We_idx(i,j)) = &
                                      5d-1*( dalphaL_dy_vf(We_idx(i,j))%sf(k, l ,r) &
                                           + dalphaR_dy_vf(We_idx(i,j))%sf(k,l+1,r) )
                                    
                                    dalpha_avg_dz(We_idx(i,j)) = &
                                      5d-1*( dalphaL_dz_vf(We_idx(i,j))%sf(k, l ,r) &
                                           + dalphaR_dz_vf(We_idx(i,j))%sf(k,l+1,r) )
                                    
                                    gm_alpha_avg(We_idx(i,j)) = &
                                      5d-1*( gm_alphaL_vf(We_idx(i,j))%sf(k, l ,r) &
                                           + gm_alphaR_vf(We_idx(i,j))%sf(k,l+1,r) )
                              
                                END DO
                           
                                tau_We(2,3) = &
                                                - 1d0/We_avg_rs_vf( We_idx(i,1), &
                                                      We_idx(i,2)  )%sf(l,k,r) * &
                                      ( alpha_avg_rs_vf(We_idx(i,1))%sf(l,k,r) * &
                                               (    dalpha_avg_dy(We_idx(i,2)) * &
                                                    dalpha_avg_dz(We_idx(i,2)) / &
                                                 MAX(gm_alpha_avg(We_idx(i,2)) , &
                                                                      sgm_eps) ) &
                                      + alpha_avg_rs_vf(We_idx(i,2))%sf(l,k,r) * &
                                               (    dalpha_avg_dy(We_idx(i,1)) * &
                                                    dalpha_avg_dz(We_idx(i,1)) / &
                                                 MAX(gm_alpha_avg(We_idx(i,1)) , &
                                                                      sgm_eps) ) )
                                
                                flux_src_vf(mom_idx%end)%sf(k,l,r) = &
                                flux_src_vf(mom_idx%end)%sf(k,l,r) - &
                                tau_We(2,3)
                                
                                flux_src_vf(E_idx)%sf(k,l,r) = &
                                flux_src_vf(E_idx)%sf(k,l,r) - &
                                vel_src_rs_vf(3)%sf(l,k,r)   * &
                                tau_We(2,3)
                           
                            END DO
                        END DO
                    END DO
                END DO
            ! END: Capillary Stresses in y-direction ===========================
            
            
            ! Capillary Stresses in z-direction ================================
            ELSE
               
                DO i = 1, We_size
                    DO r = iz%beg, iz%end
                        DO l = iy%beg, iy%end
                            DO k = ix%beg, ix%end
                           
                                DO j = 1,2
                              
                                    dalpha_avg_dx(We_idx(i,j)) = &
                                      5d-1*( dalphaL_dx_vf(We_idx(i,j))%sf(k,l, r ) &
                                           + dalphaR_dx_vf(We_idx(i,j))%sf(k,l,r+1) )
                                    
                                    dalpha_avg_dy(We_idx(i,j)) = &
                                      5d-1*( dalphaL_dy_vf(We_idx(i,j))%sf(k,l, r ) &
                                           + dalphaR_dy_vf(We_idx(i,j))%sf(k,l,r+1) )
                                    
                                    dalpha_avg_dz(We_idx(i,j)) = &
                                      5d-1*( dalphaL_dz_vf(We_idx(i,j))%sf(k,l, r ) &
                                           + dalphaR_dz_vf(We_idx(i,j))%sf(k,l,r+1) )
                                    
                                    gm_alpha_avg(We_idx(i,j)) = &
                                      5d-1*( gm_alphaL_vf(We_idx(i,j))%sf(k,l, r ) &
                                           + gm_alphaR_vf(We_idx(i,j))%sf(k,l,r+1) )
                              
                                END DO
                           
                                tau_We(3,1) = &
                                                - 1d0/We_avg_rs_vf( We_idx(i,1), &
                                                      We_idx(i,2)  )%sf(r,l,k) * &
                                      ( alpha_avg_rs_vf(We_idx(i,1))%sf(r,l,k) * &
                                               (    dalpha_avg_dz(We_idx(i,2)) * &
                                                    dalpha_avg_dx(We_idx(i,2)) / &
                                                 MAX(gm_alpha_avg(We_idx(i,2)) , &
                                                                      sgm_eps) ) &
                                      + alpha_avg_rs_vf(We_idx(i,2))%sf(r,l,k) * &
                                               (    dalpha_avg_dz(We_idx(i,1)) * &
                                                    dalpha_avg_dx(We_idx(i,1)) / &
                                                 MAX(gm_alpha_avg(We_idx(i,1)) , &
                                                                      sgm_eps) ) )
                                
                                tau_We(3,2) = &
                                                - 1d0/We_avg_rs_vf( We_idx(i,1), &
                                                      We_idx(i,2)  )%sf(r,l,k) * &
                                      ( alpha_avg_rs_vf(We_idx(i,1))%sf(r,l,k) * &
                                               (    dalpha_avg_dz(We_idx(i,2)) * &
                                                    dalpha_avg_dy(We_idx(i,2)) / &
                                                 MAX(gm_alpha_avg(We_idx(i,2)) , &
                                                                      sgm_eps) ) &
                                      + alpha_avg_rs_vf(We_idx(i,2))%sf(r,l,k) * &
                                               (    dalpha_avg_dz(We_idx(i,1)) * &
                                                    dalpha_avg_dy(We_idx(i,1)) / &
                                                 MAX(gm_alpha_avg(We_idx(i,1)) , &
                                                                      sgm_eps) ) )
                                
                                tau_We(3,3) = &
                                                  1d0/We_avg_rs_vf( We_idx(i,1), &
                                                      We_idx(i,2)  )%sf(r,l,k) * &
                                      ( alpha_avg_rs_vf(We_idx(i,1))%sf(r,l,k) * &
                                               (     gm_alpha_avg(We_idx(i,2)) - &
                                                    dalpha_avg_dz(We_idx(i,2)) * &
                                                    dalpha_avg_dz(We_idx(i,2)) / &
                                                 MAX(gm_alpha_avg(We_idx(i,2)) , &
                                                                      sgm_eps) ) &
                                      + alpha_avg_rs_vf(We_idx(i,2))%sf(r,l,k) * &
                                               (     gm_alpha_avg(We_idx(i,1)) - &
                                                    dalpha_avg_dz(We_idx(i,1)) * &
                                                    dalpha_avg_dz(We_idx(i,1)) / &
                                                 MAX(gm_alpha_avg(We_idx(i,1)) , &
                                                                      sgm_eps) ) )
                                
                                DO j = 1,3
                                    flux_src_vf(cont_idx%end+j)%sf(k,l,r) = &
                                    flux_src_vf(cont_idx%end+j)%sf(k,l,r) - &
                                    tau_We(3,j)
                                END DO
                           
                                tau_We(3,3) = &
                                               - 1d0/We_avg_rs_vf( We_idx(i,1), &
                                                     We_idx(i,2)  )%sf(r,l,k) * &
                                     ( alpha_avg_rs_vf(We_idx(i,1))%sf(r,l,k) * &
                                              (    dalpha_avg_dz(We_idx(i,2)) * &
                                                   dalpha_avg_dz(We_idx(i,2)) / &
                                                MAX(gm_alpha_avg(We_idx(i,2)) , &
                                                                     sgm_eps) ) &
                                     + alpha_avg_rs_vf(We_idx(i,2))%sf(r,l,k) * &
                                              (    dalpha_avg_dz(We_idx(i,1)) * &
                                                   dalpha_avg_dz(We_idx(i,1)) / &
                                                MAX(gm_alpha_avg(We_idx(i,1)) , &
                                                                     sgm_eps) ) )
                           
                                DO j = 1,3
                                    flux_src_vf(E_idx)%sf(k,l,r) = &
                                    flux_src_vf(E_idx)%sf(k,l,r) - &
                                    vel_src_rs_vf(j)%sf(r,l,k)   * &
                                    tau_We(3,j)
                                END DO
                           
                            END DO
                        END DO
                    END DO
                END DO
               
            END IF
            ! END: Capillary Stresses in z-direction ===========================
            
            
        END SUBROUTINE s_compute_capillary_source_flux ! -----------------------
        
        
        
        !>  Deallocation and/or disassociation procedures that are
        !!      needed to finalize the selected Riemann problem solver
        !!  @param flux_vf       Intercell fluxes 
        !!  @param flux_src_vf   Intercell source fluxes 
        !!  @param flux_gsrc_vf  Intercell geometric source fluxes 
        !!  @param norm_dir Dimensional splitting coordinate direction
        !!  @param ix   Index bounds in  first coordinate direction
        !!  @param iy   Index bounds in second coordinate direction
        !!  @param iz   Index bounds in  third coordinate direction
        SUBROUTINE s_finalize_riemann_solver( flux_vf, flux_src_vf, & ! --------
                                                      flux_gsrc_vf, &
                                                norm_dir, ix,iy,iz  )
           
            TYPE(scalar_field),  &
            DIMENSION(sys_size), &
            INTENT(INOUT) :: flux_vf, flux_src_vf, flux_gsrc_vf
            
            INTEGER, INTENT(IN) :: norm_dir
            
            TYPE(bounds_info), INTENT(IN) :: ix,iy,iz
            
            INTEGER :: i,j,k !< Generic loop iterators
            
            
            ! Reshaping Outputted Data in y-direction ==========================
            IF(norm_dir == 2) THEN
               
                DO i = 1, sys_size
                    DO k = iy%beg, iy%end
                        DO j = ix%beg, ix%end
                            flux_vf(i)%sf(j,k,iz%beg:iz%end) = &
                            flux_rs_vf(i)%sf(k,j,:)
                            flux_gsrc_vf(i)%sf(j,k,iz%beg:iz%end) = &
                            flux_gsrc_rs_vf(i)%sf(k,j,:)
                        END DO
                    END DO
                END DO
               
                DO k = iy%beg, iy%end
                    DO j = ix%beg, ix%end
                        flux_src_vf(adv_idx%beg)%sf(j,k,iz%beg:iz%end) = &
                        flux_src_rs_vf(adv_idx%beg)%sf(k,j,:)
                    END DO
                END DO
               
                IF(riemann_solver == 1) THEN
                  
                    DO i = adv_idx%beg+1, sys_size
                        DO k = iy%beg, iy%end
                            DO j = ix%beg, ix%end
                                flux_src_vf(i)%sf(j,k,iz%beg:iz%end) = &
                                flux_src_rs_vf(i)%sf(k,j,:)
                            END DO
                        END DO
                    END DO
                  
                END IF
            ! ==================================================================
            
            
            ! Reshaping Outputted Data in z-direction ==========================
            ELSEIF(norm_dir == 3) THEN
               
                DO i = 1, sys_size
                    DO k = iz%beg, iz%end
                        DO j = ix%beg, ix%end
                            flux_vf(i)%sf(j,iy%beg:iy%end,k) = &
                            flux_rs_vf(i)%sf(k,:,j)
                            flux_gsrc_vf(i)%sf(j,iy%beg:iy%end,k) = &
                            flux_gsrc_rs_vf(i)%sf(k,:,j)
                        END DO
                    END DO
                END DO
               
                DO k = iz%beg, iz%end
                    DO j = ix%beg, ix%end
                        flux_src_vf(adv_idx%beg)%sf(j,iy%beg:iy%end,k) = &
                        flux_src_rs_vf(adv_idx%beg)%sf(k,:,j)
                    END DO
                END DO
               
                IF(riemann_solver == 1) THEN
                  
                    DO i = adv_idx%beg+1, sys_size
                        DO k = iz%beg, iz%end
                            DO j = ix%beg, ix%end
                                flux_src_vf(i)%sf(j,iy%beg:iy%end,k) = &
                                flux_src_rs_vf(i)%sf(k,:,j)
                            END DO
                        END DO
                    END DO
                  
                END IF
               
            END IF
            ! ==================================================================
            
            
            ! Deallocating Left, Right and Average Riemann Problem States ======
            DO i = 1, sys_size
                DEALLOCATE(qL_prim_rs_vf(i)%sf, qR_prim_rs_vf(i)%sf)
                IF (tvd_riemann_flux) DEALLOCATE(q_prim_rs_vf(i)%sf)
            END DO
            
            DO i = 1, crv_size
                DEALLOCATE(   kappaL_rs_vf(crv_idx(i))%sf)
                DEALLOCATE(   kappaR_rs_vf(crv_idx(i))%sf)
                DEALLOCATE(alpha_avg_rs_vf(crv_idx(i))%sf)
            END DO
            
            DO i = 1,2
                IF(Re_size(i) > 0) THEN
                    DEALLOCATE(Re_avg_rs_vf(i)%sf)
                END IF
            END DO
            
            DO i = 1, We_size
                DEALLOCATE(We_avg_rs_vf(We_idx(i,1),We_idx(i,2))%sf)
            END DO
            ! ==================================================================
            
            
            ! Deallocating Intercell Fluxes and Velocity =======================
            IF(norm_dir /= 1) THEN
               
                DO i = 1, sys_size
                    DEALLOCATE(flux_rs_vf(i)%sf)
                    DEALLOCATE(flux_gsrc_rs_vf(i)%sf)
                END DO
               
                IF(riemann_solver == 1) THEN
                    DO i = adv_idx%beg+1, sys_size
                        DEALLOCATE(flux_src_rs_vf(i)%sf)
                    END DO
                END IF
               
                DEALLOCATE(flux_src_rs_vf(adv_idx%beg)%sf)
               
            ELSE
               
                DO i = 1, sys_size
                    flux_rs_vf(i)%sf     => NULL()
                    flux_src_rs_vf(i)%sf => NULL()
                    flux_gsrc_rs_vf(i)%sf => NULL()
                END DO
               
            END IF
            
            IF(riemann_solver == 1) THEN
                DEALLOCATE(vel_src_rs_vf(dir_idx(1))%sf)
            ELSE
                vel_src_rs_vf(dir_idx(1))%sf => NULL()
            END IF
            
            DO i = 2, num_dims
                DEALLOCATE(vel_src_rs_vf(dir_idx(i))%sf)
            END DO
            ! ==================================================================
            
            
        END SUBROUTINE s_finalize_riemann_solver ! -----------------------------
        
        
        
        
        
        !> Module deallocation and/or disassociation procedures
        SUBROUTINE s_finalize_riemann_solvers_module() ! -----------------------

            
            
            ! Deallocating the variables that were utilized to formulate the
            ! left, right and average states of the Riemann problem, as well
            ! the Riemann problem solution
            DEALLOCATE(qL_prim_rs_vf, qR_prim_rs_vf)
            
            DEALLOCATE(flux_rs_vf, flux_src_rs_vf, flux_gsrc_rs_vf)
            
            DEALLOCATE(vel_src_rs_vf)
            
            DEALLOCATE(alpha_rho_L, vel_L)
            DEALLOCATE(alpha_rho_R, vel_R)
            
            DEALLOCATE(vel_avg)
            
            DEALLOCATE(alpha_L, alpha_R)
            
            IF(We_size > 0) THEN
                
                DEALLOCATE(kappaL_rs_vf, kappaR_rs_vf)
                DEALLOCATE(alpha_avg_rs_vf, kappa_avg, We_avg_rs_vf)
                
            END IF
            
            DEALLOCATE(We_L, We_R)
            
            IF(ANY(Re_size > 0)) DEALLOCATE(Re_avg_rs_vf)
            
            IF (tvd_riemann_flux) THEN
                DEALLOCATE(q_prim_rs_vf)
                DEALLOCATE(lo_alpha_rho_L, lo_vel_L)
                DEALLOCATE(lo_alpha_rho_R, lo_vel_R)
                DEALLOCATE(lo_alpha_L)
                DEALLOCATE(lo_alpha_R)
                DEALLOCATE(lo_We_L, lo_We_R)
            END IF

            IF (riemann_solver == 3) THEN
                DEALLOCATE(alpha_rho_IC, vel_IC)
                DEALLOCATE(alpha_IC)
            END IF
            
            IF (bubbles) THEN
                IF (qbmm) THEN
                    DEALLOCATE(moms_L, moms_R)
                END IF
                DEALLOCATE( R0_L,R0_R,pbw_L,pbw_R)
                DEALLOCATE(V0_L,V0_R)
            END IF

            ! Disassociating procedural pointer to the subroutine which was
            ! utilized to calculate the solution of a given Riemann problem
            s_riemann_solver => NULL()
            
            
            ! Disassociating the procedural pointers to the procedures that were
            ! utilized to compute the average state and estimate the wave speeds
            s_compute_average_state => NULL(); s_compute_wave_speeds => NULL()
            
           
            ! Disassociating procedural pointer to the subroutine which was
            ! utilized to calculate the viscous source flux
            s_compute_viscous_source_flux => NULL()


            ! Disassociating the pointer to the procedure that was utilized to
            ! to convert mixture or species variables to the mixture variables
            s_convert_to_mixture_variables => NULL()
            
            
        END SUBROUTINE s_finalize_riemann_solvers_module ! ---------------------
        
        
        
        
        
END MODULE m_riemann_solvers
