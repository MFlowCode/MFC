!!       __  _______________
!!      /  |/  / ____/ ____/
!!     / /|_/ / /_  / /     
!!    / /  / / __/ / /___   
!!   /_/  /_/_/    \____/   
!!                       
!!  This file is part of MFC.
!!
!! Copyright 2021
!!
!! Permission is hereby granted, free of charge, to any person 
!! obtaining a copy of this software and associated documentation 
!! files (the "Software"), to deal in the Software without 
!! restriction, including without limitation the rights to use, 
!! copy, modify, merge, publish, distribute, sublicense, 
!! and/or sell copies of the Software, and to permit persons 
!! to whom the Software is furnished to do so, subject to the 
!! following conditions:
!!
!! The above copyright notice and this permission notice shall 
!! be included in all copies or substantial portions of the Software.
!!
!! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
!! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
!! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
!! IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
!! FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION 
!! OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
!! CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 
!! THE SOFTWARE.

!>
!! @file m_cbc.f90
!! @brief Contains module m_cbc
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief The module features a large database of characteristic boundary
!!              conditions (CBC) for the Euler system of equations. This system
!!              is augmented by the appropriate advection equations utilized to
!!              capture the material interfaces. The closure is achieved by the
!!              stiffened equation of state and mixture relations. At this time,
!!              the following CBC are available:
!!                           1) Slip Wall
!!                           2) Nonreflecting Subsonic Buffer
!!                           3) Nonreflecting Subsonic Inflow
!!                           4) Nonreflecting Subsonic Outflow
!!                           5) Force-Free Subsonic Outflow
!!                           6) Constant Pressure Subsonic Outflow
!!                           7) Supersonic Inflow
!!                           8) Supersonic Outflow
!!              Please refer to Thompson (1987, 1990) for detailed descriptions.
MODULE m_cbc
    
    
    ! Dependencies =============================================================
    USE m_derived_types        !< Definitions of the derived types
    
    USE m_global_parameters    !< Definitions of the global parameters
    
    USE m_variables_conversion !< State variables type conversion procedures
    ! ==========================================================================
    
    
    IMPLICIT NONE
    
    PRIVATE; PUBLIC :: s_initialize_cbc_module, s_cbc, s_finalize_cbc_module
    
    
    ABSTRACT INTERFACE ! =======================================================
        
        !> Abstract interface to the procedures that are utilized to calculate
        !! the L variables. For additional information refer to the following:
        !!            1) s_compute_slip_wall_L
        !!            2) s_compute_nonreflecting_subsonic_buffer_L
        !!            3) s_compute_nonreflecting_subsonic_inflow_L
        !!            4) s_compute_nonreflecting_subsonic_outflow_L
        !!            5) s_compute_force_free_subsonic_outflow_L
        !!            6) s_compute_constant_pressure_subsonic_outflow_L
        !!            7) s_compute_supersonic_inflow_L
        !!            8) s_compute_supersonic_outflow_L
        !! @param dflt_int Default null integer
        SUBROUTINE s_compute_abstract_L(dflt_int)

        INTEGER, INTENT(IN) :: dflt_int
            
        END SUBROUTINE s_compute_abstract_L
        
    END INTERFACE ! ============================================================
    
    
    TYPE(scalar_field), ALLOCATABLE, DIMENSION(:) :: q_prim_rs_vf !<
    !! The cell-average primitive variables. They are obtained by reshaping (RS)
    !! q_prim_vf in the coordinate direction normal to the domain boundary along
    !! which the CBC is applied.
 
    TYPE(scalar_field), ALLOCATABLE, DIMENSION(:) :: F_rs_vf, F_src_rs_vf !<
    !! Cell-average fluxes (src - source). These are directly determined from the
    !! cell-average primitive variables, q_prims_rs_vf, and not a Riemann solver.
    
    TYPE(scalar_field), ALLOCATABLE, DIMENSION(:) :: flux_rs_vf, flux_src_rs_vf !<
    !! The cell-boundary-average of the fluxes. They are initially determined by
    !! reshaping flux_vf and flux_src_vf in a coordinate direction normal to the
    !! domain boundary along which CBC is applied. flux_rs_vf and flux_src_rs_vf
    !! are subsequently modified based on the selected CBC.
    
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:)   :: alpha_rho   !< Cell averaged partial densiy
    REAL(KIND(0d0))                              :: rho         !< Cell averaged density
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:)   :: vel         !< Cell averaged velocity
    REAL(KIND(0d0))                              :: pres        !< Cell averaged pressure
    REAL(KIND(0d0))                              :: E           !< Cell averaged energy
    REAL(KIND(0d0))                              :: H           !< Cell averaged enthalpy
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:)   :: adv         !< Cell averaged advected variables
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:)   :: mf          !< Cell averaged mass fraction
    REAL(KIND(0d0))                              :: gamma       !< Cell averaged specific heat ratio
    REAL(KIND(0d0))                              :: pi_inf      !< Cell averaged liquid stiffness
    REAL(KIND(0d0))                              :: c           !< Cell averaged speed of sound
    REAL(KIND(0d0)),              DIMENSION(2)   :: Re          !< Cell averaged Reynolds numbers
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:) :: We          !< Cell averaged Weber numbers

    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:) :: dalpha_rho_ds !< Spatial derivatives in s-dir of partial density
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:) ::       dvel_ds !< Spatial derivatives in s-dir of velocity
    REAL(KIND(0d0))                            ::      dpres_ds !< Spatial derivatives in s-dir of pressure
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:) ::       dadv_ds !< Spatial derivatives in s-dir of advection variables
    !! Note that these are only obtained in those cells on the domain boundary along which the
    !! CBC is applied by employing finite differences (FD) on the cell-average primitive variables, q_prim_rs_vf.
 
    
    REAL(KIND(0d0)),              DIMENSION(3) :: lambda !< Eigenvalues (see Thompson 1987,1990)
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:) :: L      !< L matrix (see Thompson 1987,1990)
    

    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:) :: ds !< Cell-width distribution in the s-direction
    
    
    ! CBC Coefficients =========================================================
    REAL(KIND(0d0)), TARGET, ALLOCATABLE, DIMENSION(:,:) :: fd_coef_x !< Finite diff. coefficients x-dir
    REAL(KIND(0d0)), TARGET, ALLOCATABLE, DIMENSION(:,:) :: fd_coef_y !< Finite diff. coefficients y-dir
    REAL(KIND(0d0)), TARGET, ALLOCATABLE, DIMENSION(:,:) :: fd_coef_z !< Finite diff. coefficients z-dir
    !! The first dimension identifies the location of a coefficient in the FD
    !! formula, while the last dimension denotes the location of the CBC.


    REAL(KIND(0d0)), POINTER, DIMENSION(:,:) :: fd_coef => NULL()
    
    REAL(KIND(0d0)), TARGET, ALLOCATABLE, DIMENSION(:,:,:) :: pi_coef_x !< Polynominal interpolant coefficients in x-dir
    REAL(KIND(0d0)), TARGET, ALLOCATABLE, DIMENSION(:,:,:) :: pi_coef_y !< Polynominal interpolant coefficients in y-dir
    REAL(KIND(0d0)), TARGET, ALLOCATABLE, DIMENSION(:,:,:) :: pi_coef_z !< Polynominal interpolant coefficients in z-dir
    !! The first dimension of the array identifies the polynomial, the
    !! second dimension identifies the position of its coefficients and the last
    !! dimension denotes the location of the CBC.

    REAL(KIND(0d0)), POINTER, DIMENSION(:,:,:) :: pi_coef => NULL()
    ! ==========================================================================
    
    

    PROCEDURE(s_compute_abstract_L), POINTER :: s_compute_L => NULL() !<
    !! Pointer to procedure used to calculate L variables, based on choice of CBC
    

    TYPE(bounds_info) :: is1,is2,is3 !< Indical bounds in the s1-, s2- and s3-directions
    
    
    CONTAINS
        
        !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.        
        SUBROUTINE s_initialize_cbc_module() ! ---------------------------------

            
            
            IF(              ALL((/bc_x%beg,bc_x%end/) > -5)  &
                                    .AND.                     &
                (n > 0 .AND. ALL((/bc_y%beg,bc_y%end/) > -5)) &
                                    .AND.                     &
                (p > 0 .AND. ALL((/bc_z%beg,bc_z%end/) > -5)) ) RETURN
            
            
            ! Allocating the cell-average primitive variables
            ALLOCATE(q_prim_rs_vf(1:sys_size))

            ! Allocating the cell-average and cell-boundary-average fluxes
            ALLOCATE(   F_rs_vf(1:sys_size),    F_src_rs_vf(1:sys_size))
            ALLOCATE(flux_rs_vf(1:sys_size), flux_src_rs_vf(1:sys_size))
            
            ! Allocating the cell-average partial densities, the velocity, the
            ! advected variables, the mass fractions, as well as Weber numbers
            ALLOCATE(alpha_rho(1: cont_idx%end ))
            ALLOCATE(      vel(1:   num_dims   ))
            ALLOCATE(      adv(1:adv_idx%end-E_idx))
            ALLOCATE(       mf(1: cont_idx%end ))
            
            ALLOCATE(We(1:num_fluids,1:num_fluids))
            
            ! Allocating the first-order spatial derivatives in the s-direction
            ! of the partial densities, the velocity and the advected variables
            ALLOCATE(dalpha_rho_ds(1: cont_idx%end ))
            ALLOCATE(      dvel_ds(1:  num_dims    ))
            ALLOCATE(      dadv_ds(1:adv_idx%end-E_idx))
            
            ! Allocating L, see Thompson (1987, 1990)
            ALLOCATE(L(1:adv_idx%end))
            
            ! Allocating the cell-width distribution in the s-direction
            ALLOCATE(ds(0:buff_size))
            
            
            ! Allocating/Computing CBC Coefficients in x-direction =============
            IF(ALL((/bc_x%beg,bc_x%end/) <= -5)) THEN
                
                ALLOCATE(fd_coef_x(0:buff_size,-1: 1))
                
                IF(weno_order > 1) THEN
                    ALLOCATE(pi_coef_x(0:weno_polyn-1,0:weno_order-3,-1: 1))
                END IF
                
                CALL s_compute_cbc_coefficients(1,-1)
                CALL s_compute_cbc_coefficients(1, 1)
                
            ELSEIF(bc_x%beg <= -5) THEN
                
                ALLOCATE(fd_coef_x(0:buff_size,-1:-1))
                
                IF(weno_order > 1) THEN
                    ALLOCATE(pi_coef_x(0:weno_polyn-1,0:weno_order-3,-1:-1))
                END IF
                
                CALL s_compute_cbc_coefficients(1,-1)
                
            ELSEIF(bc_x%end <= -5) THEN
                
                ALLOCATE(fd_coef_x(0:buff_size, 1: 1))
                
                IF(weno_order > 1) THEN
                    ALLOCATE(pi_coef_x(0:weno_polyn-1,0:weno_order-3, 1: 1))
                END IF
                
                CALL s_compute_cbc_coefficients(1, 1)
                
            END IF
            ! ==================================================================
            
            
            ! Allocating/Computing CBC Coefficients in y-direction =============
            IF(n > 0) THEN
                
                IF(ALL((/bc_y%beg,bc_y%end/) <= -5)) THEN
                    
                    ALLOCATE(fd_coef_y(0:buff_size,-1: 1))
                    
                    IF(weno_order > 1) THEN
                        ALLOCATE(pi_coef_y(0:weno_polyn-1,0:weno_order-3,-1: 1))
                    END IF
                    
                    CALL s_compute_cbc_coefficients(2,-1)
                    CALL s_compute_cbc_coefficients(2, 1)
                    
                ELSEIF(bc_y%beg <= -5) THEN
                    
                    ALLOCATE(fd_coef_y(0:buff_size,-1:-1))
                    
                    IF(weno_order > 1) THEN
                        ALLOCATE(pi_coef_y(0:weno_polyn-1,0:weno_order-3,-1:-1))
                    END IF
                    
                    CALL s_compute_cbc_coefficients(2,-1)
                    
                ELSEIF(bc_y%end <= -5) THEN
                    
                    ALLOCATE(fd_coef_y(0:buff_size, 1: 1))
                    
                    IF(weno_order > 1) THEN
                        ALLOCATE(pi_coef_y(0:weno_polyn-1,0:weno_order-3, 1: 1))
                    END IF
                    
                    CALL s_compute_cbc_coefficients(2, 1)
                    
                END IF
                
            END IF
            ! ==================================================================
            
            
            ! Allocating/Computing CBC Coefficients in z-direction =============
            IF(p > 0) THEN
                
                IF(ALL((/bc_z%beg,bc_z%end/) <= -5)) THEN
                    
                    ALLOCATE(fd_coef_z(0:buff_size,-1: 1))
                    
                    IF(weno_order > 1) THEN
                        ALLOCATE(pi_coef_z(0:weno_polyn-1,0:weno_order-3,-1: 1))
                    END IF
                    
                    CALL s_compute_cbc_coefficients(3,-1)
                    CALL s_compute_cbc_coefficients(3, 1)
                    
                ELSEIF(bc_z%beg <= -5) THEN
                    
                    ALLOCATE(fd_coef_z(0:buff_size,-1:-1))
                    
                    IF(weno_order > 1) THEN
                        ALLOCATE(pi_coef_z(0:weno_polyn-1,0:weno_order-3,-1:-1))
                    END IF
                    
                    CALL s_compute_cbc_coefficients(3,-1)
                    
                ELSEIF(bc_z%end <= -5) THEN
                    
                    ALLOCATE(fd_coef_z(0:buff_size, 1: 1))
                    
                    IF(weno_order > 1) THEN
                        ALLOCATE(pi_coef_z(0:weno_polyn-1,0:weno_order-3, 1: 1))
                    END IF
                    
                    CALL s_compute_cbc_coefficients(3, 1)
                    
                END IF
                
            END IF
            ! ==================================================================
            
            
            ! Associating the procedural pointer to the appropriate subroutine
            ! that will be utilized in the conversion to the mixture variables
            
            IF (model_eqns == 1) THEN        ! Gamma/pi_inf model
                s_convert_to_mixture_variables => &
                             s_convert_mixture_to_mixture_variables
            ELSEIF (bubbles) THEN                            ! Volume fraction model
                s_convert_to_mixture_variables => &
                             s_convert_species_to_mixture_variables_bubbles    
            ELSE                            ! Volume fraction model
                s_convert_to_mixture_variables => &
                             s_convert_species_to_mixture_variables
            END IF
            
            
        END SUBROUTINE s_initialize_cbc_module ! -------------------------------
        
        
        
        
        
        SUBROUTINE s_compute_cbc_coefficients(cbc_dir, cbc_loc) ! --------------
        ! Description: The purpose of this subroutine is to compute the grid
        !              dependent FD and PI coefficients, or CBC coefficients,
        !              provided the CBC coordinate direction and location.
            
            
            ! CBC coordinate direction and location
            INTEGER, INTENT(IN) :: cbc_dir, cbc_loc
            
            ! Cell-boundary locations in the s-direction
            REAL(KIND(0d0)), DIMENSION(0:buff_size+1) :: s_cb
            
            ! Generic loop iterator
            INTEGER :: i
            
            
            ! Associating CBC coefficients pointers
            CALL s_associate_cbc_coefficients_pointers(cbc_dir, cbc_loc)
            
            
            ! Determining the cell-boundary locations in the s-direction
            s_cb(0) = 0d0
            
            DO i = 0, buff_size
                s_cb(i+1) = s_cb(i) + ds(i)
            END DO
            
            
            ! Computing CBC1 Coefficients ======================================
            IF(weno_order == 1) THEN
                
                fd_coef(:,cbc_loc) =  0d0
                fd_coef(0,cbc_loc) = -2d0/(ds(0)+ds(1))
                fd_coef(1,cbc_loc) = -fd_coef(0,cbc_loc)
                
            ! ==================================================================
            
            
            ! Computing CBC2 Coefficients ======================================
            ELSEIF(weno_order == 3) THEN
                
                fd_coef(:,cbc_loc) =  0d0
                fd_coef(0,cbc_loc) = -6d0/(3d0*ds(0)+2d0*ds(1)-ds(2))
                fd_coef(1,cbc_loc) = -4d0*fd_coef(0,cbc_loc)/3d0
                fd_coef(2,cbc_loc) =      fd_coef(0,cbc_loc)/3d0
                
                pi_coef(0,0,cbc_loc) = (s_cb(0)-s_cb(1))/(s_cb(0)-s_cb(2))
                
            ! ==================================================================
            
            
            ! Computing CBC4 Coefficients ======================================
            ELSE
                
                fd_coef(:,cbc_loc) =  0d0
                fd_coef(0,cbc_loc) = -50d0/(25d0*ds(0)+2d0*ds(1) &
                                           - 1d1*ds(2)+1d1*ds(3) &
                                           - 3d0*ds(4)           )
                fd_coef(1,cbc_loc) = -48d0*fd_coef(0,cbc_loc)/25d0
                fd_coef(2,cbc_loc) =  36d0*fd_coef(0,cbc_loc)/25d0
                fd_coef(3,cbc_loc) = -16d0*fd_coef(0,cbc_loc)/25d0
                fd_coef(4,cbc_loc) =   3d0*fd_coef(0,cbc_loc)/25d0
                
                pi_coef(0,0,cbc_loc) = &
                            ((s_cb(0)-s_cb(1))*(s_cb(1)-s_cb(2))  * &
                            (s_cb(1)-s_cb(3)))/((s_cb(1)-s_cb(4)) * &
                            (s_cb(4)-s_cb(0))*(s_cb(4)-s_cb(2)))
                pi_coef(0,1,cbc_loc) = &
                            ((s_cb(1)-s_cb(0))*(s_cb(1)-s_cb(2))  * &
                            ((s_cb(1)-s_cb(3))*(s_cb(1)-s_cb(3))  - &
                            (s_cb(0)-s_cb(4))*((s_cb(3)-s_cb(1))  + &
                            (s_cb(4)-s_cb(1)))))                  / &
                            ((s_cb(0)-s_cb(3))*(s_cb(1)-s_cb(3))  * &
                            (s_cb(0)-s_cb(4))*(s_cb(1)-s_cb(4)))
                pi_coef(0,2,cbc_loc) = &
                            (s_cb(1)-s_cb(0))*((s_cb(1)-s_cb(2))  * &
                            (s_cb(1)-s_cb(3))+((s_cb(0)-s_cb(2))  + &
                            (s_cb(1)-s_cb(3)))*(s_cb(0)-s_cb(4))) / &
                            ((s_cb(2)-s_cb(0))*(s_cb(0)-s_cb(3))  * &
                            (s_cb(0)-s_cb(4)))
                pi_coef(1,0,cbc_loc) = &
                            ((s_cb(0)-s_cb(2))*(s_cb(2)-s_cb(1))  * &
                            (s_cb(2)-s_cb(3)))/((s_cb(2)-s_cb(4)) * &
                            (s_cb(4)-s_cb(0))*(s_cb(4)-s_cb(1)))
                pi_coef(1,1,cbc_loc) = &
                            ((s_cb(0)-s_cb(2))*(s_cb(1)-s_cb(2))  * &
                            ((s_cb(1)-s_cb(3))*(s_cb(2)-s_cb(3))  + &
                            (s_cb(0)-s_cb(4))*((s_cb(1)-s_cb(3))  + &
                            (s_cb(2)-s_cb(4)))))                  / &
                            ((s_cb(0)-s_cb(3))*(s_cb(1)-s_cb(3))  * &
                            (s_cb(0)-s_cb(4))*(s_cb(1)-s_cb(4)))
                pi_coef(1,2,cbc_loc) = &
                            ((s_cb(1)-s_cb(2))*(s_cb(2)-s_cb(3))  * &
                            (s_cb(2)-s_cb(4)))/((s_cb(0)-s_cb(2)) * &
                            (s_cb(0)-s_cb(3))*(s_cb(0)-s_cb(4)))
                
            END IF
            ! END: Computing CBC4 Coefficients =================================
            
            
            ! Nullifying CBC coefficients
            NULLIFY(fd_coef, pi_coef)
            
            
        END SUBROUTINE s_compute_cbc_coefficients ! ----------------------------
        
        
        
        !!  The goal of the procedure is to associate the FD and PI
        !!      coefficients, or CBC coefficients, with the appropriate
        !!      targets, based on the coordinate direction and location
        !!      of the CBC.
        !!  @param cbc_dir CBC coordinate direction
        !!  @param cbc_loc CBC coordinate location
        SUBROUTINE s_associate_cbc_coefficients_pointers(cbc_dir, cbc_loc) ! ---
            
            
            INTEGER, INTENT(IN) :: cbc_dir, cbc_loc

            INTEGER :: i !< Generic loop iterator
            
            
            ! Associating CBC Coefficients in x-direction ======================
            IF(cbc_dir == 1) THEN
                
                fd_coef => fd_coef_x; IF(weno_order > 1) pi_coef => pi_coef_x
                
                IF(cbc_loc == -1) THEN
                    DO i = 0, buff_size
                        ds(i) = dx( i )
                    END DO
                ELSE
                    DO i = 0, buff_size
                        ds(i) = dx(m-i)
                    END DO
                END IF
            ! ==================================================================
            
            
            ! Associating CBC Coefficients in y-direction ======================
            ELSEIF(cbc_dir == 2) THEN
                
                fd_coef => fd_coef_y; IF(weno_order > 1) pi_coef => pi_coef_y
                
                IF(cbc_loc == -1) THEN
                    DO i = 0, buff_size
                        ds(i) = dy( i )
                    END DO
                ELSE
                    DO i = 0, buff_size
                        ds(i) = dy(n-i)
                    END DO
                END IF
            ! ==================================================================
            
            
            ! Associating CBC Coefficients in z-direction ======================
            ELSE
                
                fd_coef => fd_coef_z; IF(weno_order > 1) pi_coef => pi_coef_z
                
                IF(cbc_loc == -1) THEN
                    DO i = 0, buff_size
                        ds(i) = dz( i )
                    END DO
                ELSE
                    DO i = 0, buff_size
                        ds(i) = dz(p-i)
                    END DO
                END IF
                
            END IF
            ! ==================================================================
            
            
        END SUBROUTINE s_associate_cbc_coefficients_pointers ! -----------------
        
        
        
        
        !>  The following is the implementation of the CBC based on
        !!      the work of Thompson (1987, 1990) on hyperbolic systems.
        !!      The CBC is indirectly applied in the computation of the
        !!      right-hand-side (RHS) near the relevant domain boundary
        !!      through the modification of the fluxes.
        !!  @param q_prim_vf Cell-average primitive variables
        !!  @param flux_vf Cell-boundary-average fluxes
        !!  @param flux_src_vf Cell-boundary-average flux sources
        !!  @param cbc_dir CBC coordinate direction
        !!  @param cbc_loc CBC coordinate location
        !!  @param ix Index bound in the first coordinate direction 
        !!  @param iy Index bound in the second coordinate direction 
        !!  @param iz Index bound in the third coordinate direction 
        SUBROUTINE s_cbc( q_prim_vf, flux_vf, flux_src_vf, & ! -----------------
                                         cbc_dir, cbc_loc, &
                                                 ix,iy,iz  )
           
            TYPE(scalar_field),  &
            DIMENSION(sys_size), &
            INTENT(IN) :: q_prim_vf
            
            TYPE(scalar_field),  &
            DIMENSION(sys_size), &
            INTENT(INOUT) :: flux_vf, flux_src_vf
            
            INTEGER, INTENT(IN) :: cbc_dir, cbc_loc
            
            TYPE(bounds_info), INTENT(IN) :: ix,iy,iz
            
            ! First-order time derivatives of the partial densities, density,
            ! velocity, pressure, advection variables, and the specific heat
            ! ratio and liquid stiffness functions
            REAL(KIND(0d0)), DIMENSION(cont_idx%end)   :: dalpha_rho_dt
            REAL(KIND(0d0))                            ::       drho_dt
            REAL(KIND(0d0)), DIMENSION(num_dims)       ::       dvel_dt
            REAL(KIND(0d0))                            ::      dpres_dt
            REAL(KIND(0d0)), DIMENSION(adv_idx%end-E_idx) ::       dadv_dt
            REAL(KIND(0d0))                            ::     dgamma_dt
            REAL(KIND(0d0))                            ::    dpi_inf_dt
            
            INTEGER :: i,j,k,r !< Generic loop iterators


            REAL(KIND(0d0)) :: blkmod1, blkmod2 !< Fluid bulk modulus for Wood mixture sound speed
            
            
            ! Reshaping of inputted data and association of the FD and PI
            ! coefficients, or CBC coefficients, respectively, hinging on
            ! selected CBC coordinate direction
            CALL s_initialize_cbc( q_prim_vf, flux_vf, flux_src_vf, &
                                                  cbc_dir, cbc_loc, &
                                                          ix,iy,iz  )
            
            CALL s_associate_cbc_coefficients_pointers(cbc_dir, cbc_loc)
            
            
            ! PI2 of flux_rs_vf and flux_src_rs_vf at j = 1/2 ==================
            IF(weno_order == 3) THEN
               
               CALL s_convert_primitive_to_flux_variables( q_prim_rs_vf, &
                                                                F_rs_vf, &
                                                            F_src_rs_vf, &
                                                            is1,is2,is3  )
               
               DO i = 1, adv_idx%end
                  flux_rs_vf(i)%sf(0,:,:) = F_rs_vf(i)%sf(0,:,:)   &
                                          + pi_coef(0,0,cbc_loc) * &
                                          ( F_rs_vf(i)%sf(1,:,:) - &
                                            F_rs_vf(i)%sf(0,:,:) )
               END DO
               
               DO i = adv_idx%beg, adv_idx%end
                  flux_src_rs_vf(i)%sf(0,:,:) = F_src_rs_vf(i)%sf(0,:,:) + &
                                              ( F_src_rs_vf(i)%sf(1,:,:) - &
                                                F_src_rs_vf(i)%sf(0,:,:) ) &
                                              * pi_coef(0,0,cbc_loc)
               END DO
            ! ==================================================================
            
            
            ! PI4 of flux_rs_vf and flux_src_rs_vf at j = 1/2, 3/2 =============
            ELSEIF(weno_order == 5) THEN
               
               CALL s_convert_primitive_to_flux_variables( q_prim_rs_vf, &
                                                                F_rs_vf, &
                                                            F_src_rs_vf, &
                                                            is1,is2,is3  )
               
               DO i = 1, adv_idx%end
                  DO j = 0,1
                     flux_rs_vf(i)%sf(j,:,:) = F_rs_vf(i)%sf(j,:,:)   &
                                             + pi_coef(j,0,cbc_loc) * &
                                             ( F_rs_vf(i)%sf(3,:,:) - &
                                               F_rs_vf(i)%sf(2,:,:) ) &
                                             + pi_coef(j,1,cbc_loc) * &
                                             ( F_rs_vf(i)%sf(2,:,:) - &
                                               F_rs_vf(i)%sf(1,:,:) ) &
                                             + pi_coef(j,2,cbc_loc) * &
                                             ( F_rs_vf(i)%sf(1,:,:) - &
                                               F_rs_vf(i)%sf(0,:,:) )
                  END DO
               END DO
               
               DO i = adv_idx%beg, adv_idx%end
                  DO j = 0,1
                     flux_src_rs_vf(i)%sf(j,:,:) = F_src_rs_vf(i)%sf(j,:,:) + &
                                                 ( F_src_rs_vf(i)%sf(3,:,:) - &
                                                   F_src_rs_vf(i)%sf(2,:,:) ) &
                                                 * pi_coef(j,0,cbc_loc)     + &
                                                 ( F_src_rs_vf(i)%sf(2,:,:) - &
                                                   F_src_rs_vf(i)%sf(1,:,:) ) &
                                                 * pi_coef(j,1,cbc_loc)     + &
                                                 ( F_src_rs_vf(i)%sf(1,:,:) - &
                                                   F_src_rs_vf(i)%sf(0,:,:) ) &
                                                 * pi_coef(j,2,cbc_loc)
                  END DO
               END DO
               
            END IF
            ! ==================================================================
            
            
            ! FD2 or FD4 of RHS at j = 0 =======================================
            DO r = is3%beg, is3%end
               DO k = is2%beg, is2%end
                  
                  ! Transferring the Primitive Variables =======================
                  DO i = 1, cont_idx%end
                     alpha_rho(i) = q_prim_rs_vf(i)%sf(0,k,r)
                  END DO
                  
                  DO i = 1, num_dims
                     vel(i) = q_prim_rs_vf(cont_idx%end+i)%sf(0,k,r)
                  END DO
                  
                  pres = q_prim_rs_vf(E_idx)%sf(0,k,r)
                  
                  CALL s_convert_to_mixture_variables( q_prim_rs_vf, &
                                                         rho, gamma, &
                                                         pi_inf, Re, &
                                                          We, 0,k,r  )
                  
                  E = gamma*pres + pi_inf + 5d-1*rho*SUM(vel**2d0)
                  
                  H = (E + pres)/rho
                  
                  DO i = 1, adv_idx%end-E_idx
                     adv(i) = q_prim_rs_vf(E_idx+i)%sf(0,k,r)
                  END DO
                  
                  mf = alpha_rho/rho
                  
                  ! Compute mixture sound speed
                  IF (alt_soundspeed .OR. regularization) THEN
                      blkmod1 = ((fluid_pp(1)%gamma +1d0)*pres + &
                          fluid_pp(1)%pi_inf)/fluid_pp(1)%gamma
                      blkmod2 = ((fluid_pp(2)%gamma +1d0)*pres + &
                          fluid_pp(2)%pi_inf)/fluid_pp(2)%gamma
                      c = (1d0/(rho*(adv(1)/blkmod1 + adv(2)/blkmod2)))
                  ELSEIF(model_eqns == 3) THEN
                      c = 0d0
                      DO i = 1, num_fluids
                          c = c + q_prim_rs_vf(i+adv_idx%beg-1)%sf(0,k,r) * (1d0/fluid_pp(i)%gamma+1d0) * &
                              (pres + fluid_pp(i)%pi_inf/(fluid_pp(i)%gamma+1d0))
                      END DO
                      c = c/rho
                  ELSE
                      c = ((H - 5d-1*SUM(vel**2d0))/gamma)
                  END IF

                  c = SQRT(c)

!                  IF (mixture_err .AND. c < 0d0) THEN
!                    c = sgm_eps
!                  ELSE
!                    c = SQRT(c)
!                  END IF

                  ! ============================================================
                  
                  ! First-Order Spatial Derivatives of Primitive Variables =====
                  dalpha_rho_ds = 0d0
                  dvel_ds       = 0d0
                  dpres_ds      = 0d0
                  dadv_ds       = 0d0
                  
                  DO j = 0, buff_size
                     
                     DO i = 1, cont_idx%end
                        dalpha_rho_ds(i) = q_prim_rs_vf(i)%sf(j,k,r) * &
                                                  fd_coef(j,cbc_loc) + &
                                                    dalpha_rho_ds(i)
                     END DO
                     
                     DO i = 1, num_dims
                        dvel_ds(i) = q_prim_rs_vf(cont_idx%end+i)%sf(j,k,r) * &
                                                         fd_coef(j,cbc_loc) + &
                                                                 dvel_ds(i)
                     END DO
                     
                     dpres_ds = q_prim_rs_vf(E_idx)%sf(j,k,r) * &
                                           fd_coef(j,cbc_loc) + &
                                                     dpres_ds
                     
                     DO i = 1, adv_idx%end-E_idx
                        dadv_ds(i) = q_prim_rs_vf(E_idx+i)%sf(j,k,r) * &
                                                  fd_coef(j,cbc_loc) + &
                                                          dadv_ds(i)
                     END DO
                     
                  END DO
                  ! ============================================================
                  
                  ! First-Order Temporal Derivatives of Primitive Variables ====
                  lambda(1) = vel(dir_idx(1)) - c
                  lambda(2) = vel(dir_idx(1))
                  lambda(3) = vel(dir_idx(1)) + c
                  
                  CALL s_compute_L(dflt_int)
                  
                  ! Be careful about the cylindrical coordinate!
                  IF(cyl_coord .AND. cbc_dir == 2 .AND. cbc_loc == 1 ) THEN  
                      dpres_dt = -5d-1*(L(adv_idx%end) + L(1)) + rho*c*c*vel(dir_idx(1)) &
                                 / y_cc(n)
                  ELSE
                      dpres_dt = -5d-1*(L(adv_idx%end) + L(1))
                  END IF

                  DO i = 1, cont_idx%end
                      dalpha_rho_dt(i) = &
                          -(L(i+1) - mf(i)*dpres_dt)/(c*c)
                  END DO
                  
                  DO i = 1, num_dims
                     dvel_dt(dir_idx(i)) = dir_flg(dir_idx(i)) * &
                              (L(1) - L(adv_idx%end))/(2d0*rho*c) + &
                                   (dir_flg(dir_idx(i)) - 1d0) * &
                                              L(mom_idx%beg+i)
                  END DO
                  
                  ! The treatment of void fraction source is unclear
                  IF(cyl_coord .AND. cbc_dir == 2 .AND. cbc_loc == 1 ) THEN  
                      DO i = 1, adv_idx%end-E_idx
                         dadv_dt(i) = -L(mom_idx%end+i) !+ adv(i) * vel(dir_idx(1))/y_cc(n)
                      END DO
                  ELSE
                      DO i = 1, adv_idx%end-E_idx
                         dadv_dt(i) = -L(mom_idx%end+i)
                      END DO
                  END IF
                  
                  drho_dt = 0d0; dgamma_dt = 0d0; dpi_inf_dt = 0d0
                  
                  IF(model_eqns == 1) THEN
                     drho_dt    = dalpha_rho_dt(1)
                     dgamma_dt  = dadv_dt(1)
                     dpi_inf_dt = dadv_dt(2)
                  ELSE
                     DO i = 1, num_fluids
                        drho_dt    = drho_dt    + dalpha_rho_dt(i)
                        dgamma_dt  = dgamma_dt  + dadv_dt(i)*fluid_pp(i)%gamma
                        dpi_inf_dt = dpi_inf_dt + dadv_dt(i)*fluid_pp(i)%pi_inf
                     END DO
                  END IF
                  ! ============================================================
                  
                  ! flux_rs_vf and flux_src_rs_vf at j = -1/2 ==================
                  DO i = 1, cont_idx%end
                     flux_rs_vf(i)%sf(-1,k,r) = flux_rs_vf(i)%sf(0,k,r) &
                                               + ds(0)*dalpha_rho_dt(i)
                  END DO
                  
                  DO i = mom_idx%beg, mom_idx%end
                     flux_rs_vf(i)%sf(-1,k,r) = flux_rs_vf(i)%sf(0,k,r) &
                                  + ds(0)*( vel(i-cont_idx%end)*drho_dt &
                                          + rho*dvel_dt(i-cont_idx%end) )
                  END DO
                  
                  flux_rs_vf(E_idx)%sf(-1,k,r) = flux_rs_vf(E_idx)%sf(0,k,r) &
                                        + ds(0)*( pres*dgamma_dt             &
                                                + gamma*dpres_dt             &
                                                + dpi_inf_dt                 &
                                                + rho*SUM(vel*dvel_dt)       &
                                                + 5d-1*drho_dt*SUM(vel**2d0) )
                  
                  IF(riemann_solver == 1) THEN
                     
                     DO i = adv_idx%beg, adv_idx%end
                        flux_rs_vf(i)%sf(-1,k,r) = 0d0
                     END DO
                     
                     DO i = adv_idx%beg, adv_idx%end
                        flux_src_rs_vf(i)%sf(-1,k,r) = &
                                1d0/MAX(ABS(vel(dir_idx(1))),sgm_eps) &
                                      *     SIGN(1d0,vel(dir_idx(1))) &
                                      * (     flux_rs_vf(i)%sf(0,k,r) &
                                        +             vel(dir_idx(1)) &
                                        * flux_src_rs_vf(i)%sf(0,k,r) &
                                        +      ds(0)*dadv_dt(i-E_idx) )
                     END DO
                     
                  ELSE
                     
                     DO i = adv_idx%beg, adv_idx%end
                        flux_rs_vf(i)%sf(-1,k,r) = flux_rs_vf(i)%sf(0,k,r) - &
                                  adv(i-E_idx)*flux_src_rs_vf(i)%sf(0,k,r) + &
                                                    ds(0)*dadv_dt(i-E_idx)
                     END DO
                     
                     DO i = adv_idx%beg, adv_idx%end
                        flux_src_rs_vf(i)%sf(-1,k,r) = 0d0
                     END DO
                     
                  END IF
                  ! END: flux_rs_vf and flux_src_rs_vf at j = -1/2 =============
                  
               END DO
            END DO
            ! END: FD2 or FD4 of RHS at j = 0 ==================================
            
            
            ! The reshaping of outputted data and disssociation of the FD and PI
            ! coefficients, or CBC coefficients, respectively, based on selected
            ! CBC coordinate direction.
            CALL s_finalize_cbc( flux_vf, flux_src_vf, &
                                     cbc_dir, cbc_loc, &
                                             ix,iy,iz  )
            
            NULLIFY(fd_coef, pi_coef)
            
            
        END SUBROUTINE s_cbc ! -------------------------------------------------
        
        
        
        !>  The L variables for the slip wall CBC, see pg. 451 of
        !!      Thompson (1990). At the slip wall (frictionless wall),
        !!      the normal component of velocity is zero at all times,
        !!      while the transverse velocities may be nonzero.        
        !!  @param dflt_int Default null integer
        SUBROUTINE s_compute_slip_wall_L(dflt_int) ! -----------------------------------

            
        INTEGER, INTENT(IN) :: dflt_int
            
            L(1) = lambda(1)*(dpres_ds - rho*c*dvel_ds(dir_idx(1)))
            
            L(2:adv_idx%end-1) = 0d0
            
            L(adv_idx%end) = L(1)
            
            
        END SUBROUTINE s_compute_slip_wall_L ! ---------------------------------
        
        
        
        
        !>  The L variables for the nonreflecting subsonic buffer CBC
        !!      see pg. 13 of Thompson (1987). The nonreflecting subsonic
        !!      buffer reduces the amplitude of any reflections caused by
        !!      outgoing waves.
        !!  @param dflt_int Default null integer
        SUBROUTINE s_compute_nonreflecting_subsonic_buffer_L(dflt_int) ! ---------------
           
            INTEGER, INTENT(IN) :: dflt_int
            

            INTEGER :: i !< Generic loop iterator
            
            
            L(1) = (5d-1 - 5d-1*SIGN(1d0,lambda(1)))*lambda(1) &
                 * (dpres_ds - rho*c*dvel_ds(dir_idx(1)))
            
            DO i = 2, mom_idx%beg
                L(i) = (5d-1 - 5d-1*SIGN(1d0,lambda(2)))*lambda(2) &
                     * (c*c*dalpha_rho_ds(i-1) - mf(i-1)*dpres_ds)
            END DO
            
            DO i = mom_idx%beg+1, mom_idx%end
                L(i) = (5d-1 - 5d-1*SIGN(1d0,lambda(2)))*lambda(2) &
                     * (dvel_ds(dir_idx(i-cont_idx%end)))
            END DO
            
            DO i = E_idx, adv_idx%end-1
                L(i) = (5d-1 - 5d-1*SIGN(1d0,lambda(2)))*lambda(2) &
                     * (dadv_ds(i-mom_idx%end))
            END DO
            
            L(adv_idx%end) = (5d-1 - 5d-1*SIGN(1d0,lambda(3)))*lambda(3) &
                        * (dpres_ds + rho*c*dvel_ds(dir_idx(1)))
            
            
        END SUBROUTINE s_compute_nonreflecting_subsonic_buffer_L ! -------------
        
        
        
        !>  The L variables for the nonreflecting subsonic inflow CBC
        !!      see pg. 455, Thompson (1990). This nonreflecting subsonic
        !!      CBC assumes an incoming flow and reduces the amplitude of
        !!      any reflections caused by outgoing waves.
        !! @param dflt_int Default null integer
        SUBROUTINE s_compute_nonreflecting_subsonic_inflow_L(dflt_int) ! ---------------
           
            INTEGER, INTENT(IN) :: dflt_int
            
            L(1) = lambda(1)*(dpres_ds - rho*c*dvel_ds(dir_idx(1)))
            
            L(2:adv_idx%end) = 0d0
            
            
        END SUBROUTINE s_compute_nonreflecting_subsonic_inflow_L ! -------------
        
        
        
        
        !>  The L variables for the nonreflecting subsonic outflow
        !!      CBC see pg. 454 of Thompson (1990). This nonreflecting
        !!      subsonic CBC presumes an outgoing flow and reduces the
        !!      amplitude of any reflections caused by outgoing waves.
        !! @param dflt_int Default null integer
        SUBROUTINE s_compute_nonreflecting_subsonic_outflow_L(dflt_int) ! --------------
           
            INTEGER, INTENT(IN) :: dflt_int
            

            INTEGER :: i !> Generic loop iterator
            
            
            L(1) = lambda(1)*(dpres_ds - rho*c*dvel_ds(dir_idx(1)))
            
            DO i = 2, mom_idx%beg
                L(i) = lambda(2)*(c*c*dalpha_rho_ds(i-1) - mf(i-1)*dpres_ds)
            END DO
            
            DO i = mom_idx%beg+1, mom_idx%end
                L(i) = lambda(2)*(dvel_ds(dir_idx(i-cont_idx%end)))
            END DO
            
            DO i = E_idx, adv_idx%end-1
                L(i) = lambda(2)*(dadv_ds(i-mom_idx%end))
            END DO
            
            ! bubble index
            L(adv_idx%end) = 0d0
            
            
        END SUBROUTINE s_compute_nonreflecting_subsonic_outflow_L ! ------------
        
        
        !>  The L variables for the force-free subsonic outflow CBC,
        !!      see pg. 454 of Thompson (1990). The force-free subsonic
        !!      outflow sets to zero the sum of all of the forces which
        !!      are acting on a fluid element for the normal coordinate
        !!      direction to the boundary. As a result, a fluid element
        !!      at the boundary is simply advected outward at the fluid
        !!      velocity.
        !! @param dflt_int Default null integer
        SUBROUTINE s_compute_force_free_subsonic_outflow_L(dflt_int) ! -----------------
           
        INTEGER, INTENT(IN) :: dflt_int
            

            INTEGER :: i !> Generic loop iterator
            
            
            L(1) = lambda(1)*(dpres_ds - rho*c*dvel_ds(dir_idx(1)))
            
            DO i = 2, mom_idx%beg
                L(i) = lambda(2)*(c*c*dalpha_rho_ds(i-1) - mf(i-1)*dpres_ds)
            END DO
            
            DO i = mom_idx%beg+1, mom_idx%end
                L(i) = lambda(2)*(dvel_ds(dir_idx(i-cont_idx%end)))
            END DO
            
            DO i = E_idx, adv_idx%end-1
                L(i) = lambda(2)*(dadv_ds(i-mom_idx%end))
            END DO
            
            L(adv_idx%end) = L(1) + 2d0*rho*c*lambda(2)*dvel_ds(dir_idx(1))
            
            
        END SUBROUTINE s_compute_force_free_subsonic_outflow_L ! ---------------
        
        
        
        !>  L variables for the constant pressure subsonic outflow
        !!      CBC see pg. 455 Thompson (1990). The constant pressure
        !!      subsonic outflow maintains a fixed pressure at the CBC
        !!      boundary in absence of any transverse effects.
        !! @param dflt_int Default null integer
        SUBROUTINE s_compute_constant_pressure_subsonic_outflow_L(dflt_int) ! ----------
           
        INTEGER, INTENT(IN) :: dflt_int
            

            INTEGER :: i !> Generic loop iterator
            
            
            L(1) = lambda(1)*(dpres_ds - rho*c*dvel_ds(dir_idx(1)))
            
            DO i = 2, mom_idx%beg
                L(i) = lambda(2)*(c*c*dalpha_rho_ds(i-1) - mf(i-1)*dpres_ds)
            END DO
            
            DO i = mom_idx%beg+1, mom_idx%end
                L(i) = lambda(2)*(dvel_ds(dir_idx(i-cont_idx%end)))
            END DO
            
            DO i = E_idx, adv_idx%end-1
                L(i) = lambda(2)*(dadv_ds(i-mom_idx%end))
            END DO
            
            L(adv_idx%end) = -L(1)
            
            
        END SUBROUTINE s_compute_constant_pressure_subsonic_outflow_L ! --------
        
        
        
        
        !>  L variables for the supersonic inflow CBC, see pg. 453
        !!      Thompson (1990). The supersonic inflow CBC is a steady
        !!      state, or nearly a steady state, CBC in which only the
        !!      transverse terms may generate a time dependence at the
        !!      inflow boundary.
        !! @param dflt_int Default null integer
        SUBROUTINE s_compute_supersonic_inflow_L(dflt_int) ! ---------------------------
           
        INTEGER, INTENT(IN) :: dflt_int
            
            L = 0d0
            
            
        END SUBROUTINE s_compute_supersonic_inflow_L ! -------------------------
        
        
        
        !>  L variables for the supersonic outflow CBC, see pg. 453
        !!      of Thompson (1990). For the supersonic outflow CBC, the
        !!      flow evolution at the boundary is determined completely
        !!      by the interior data.
        !! @param dflt_int Default null integer
        SUBROUTINE s_compute_supersonic_outflow_L(dflt_int) ! --------------------------
           
        INTEGER, INTENT(IN) :: dflt_int
            

            INTEGER :: i !< Generic loop iterator
            
            L(1) = lambda(1)*(dpres_ds - rho*c*dvel_ds(dir_idx(1)))
            
            DO i = 2, mom_idx%beg
                L(i) = lambda(2)*(c*c*dalpha_rho_ds(i-1) - mf(i-1)*dpres_ds)
            END DO
            
            DO i = mom_idx%beg+1, mom_idx%end
                L(i) = lambda(2)*(dvel_ds(dir_idx(i-cont_idx%end)))
            END DO
            
            DO i = E_idx, adv_idx%end-1
                L(i) = lambda(2)*(dadv_ds(i-mom_idx%end))
            END DO
            
            L(adv_idx%end) = lambda(3)*(dpres_ds + rho*c*dvel_ds(dir_idx(1)))
            
            
        END SUBROUTINE s_compute_supersonic_outflow_L ! ------------------------
        
        
        
        !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are required for the setup of the
        !!      selected CBC.
        !!  @param q_prim_vf Cell-average primitive variables
        !!  @param flux_vf Cell-boundary-average fluxes
        !!  @param flux_src_vf Cell-boundary-average flux sources
        !!  @param cbc_dir CBC coordinate direction
        !!  @param cbc_loc CBC coordinate location
        !!  @param ix Index bound in the first coordinate direction 
        !!  @param iy Index bound in the second coordinate direction 
        !!  @param iz Index bound in the third coordinate direction         
        SUBROUTINE s_initialize_cbc( q_prim_vf, flux_vf, flux_src_vf, & ! ------
                                                    cbc_dir, cbc_loc, &
                                                            ix,iy,iz  )
           
            TYPE(scalar_field),  &
            DIMENSION(sys_size), &
            INTENT(IN) :: q_prim_vf
            
            TYPE(scalar_field),  &
            DIMENSION(sys_size), &
            INTENT(IN) :: flux_vf, flux_src_vf
            
            INTEGER, INTENT(IN) :: cbc_dir, cbc_loc
            TYPE(bounds_info), INTENT(IN) :: ix,iy,iz
            

            INTEGER :: dj !< Indical shift based on CBC location
            

            INTEGER :: i,j,k,r !< Generic loop iterators
            
            
            ! Configuring the coordinate direction indexes and flags
            IF(cbc_dir == 1) THEN
               is1%beg = 0; is1%end = buff_size; is2 = iy; is3 = iz
               dir_idx = (/1,2,3/); dir_flg = (/1d0,0d0,0d0/)
            ELSEIF(cbc_dir == 2) THEN
               is1%beg = 0; is1%end = buff_size; is2 = ix; is3 = iz
               dir_idx = (/2,1,3/); dir_flg = (/0d0,1d0,0d0/)
            ELSE
               is1%beg = 0; is1%end = buff_size; is2 = iy; is3 = ix
               dir_idx = (/3,1,2/); dir_flg = (/0d0,0d0,1d0/)
            END IF
            
            
            ! Determining the indicial shift based on CBC location
            dj = MAX(0,cbc_loc)
            
            
            ! Allocation/Association of Primitive and Flux Variables ===========
            DO i = 1, sys_size
               ALLOCATE(q_prim_rs_vf(i)%sf(    0    : buff_size, &
                                            is2%beg :  is2%end , &
                                            is3%beg :  is3%end  ))
            END DO
            
            IF(weno_order > 1) THEN
               
               DO i = 1, adv_idx%end
                  ALLOCATE(F_rs_vf(i)%sf(    0    : buff_size, &
                                          is2%beg :  is2%end , &
                                          is3%beg :  is3%end  ))
               END DO
               
               ALLOCATE(F_src_rs_vf(adv_idx%beg)%sf(    0    : buff_size, &
                                                     is2%beg :  is2%end , &
                                                     is3%beg :  is3%end  ))
               
               IF(riemann_solver == 1) THEN
                  DO i = adv_idx%beg+1, adv_idx%end
                     ALLOCATE(F_src_rs_vf(i)%sf(    0    : buff_size, &
                                                 is2%beg :  is2%end , &
                                                 is3%beg :  is3%end  ))
                  END DO
               ELSE
                  DO i = adv_idx%beg+1, adv_idx%end
                     F_src_rs_vf(i)%sf => F_src_rs_vf(adv_idx%beg)%sf
                  END DO
               END IF
               
            END IF
            
            DO i = 1, adv_idx%end
               ALLOCATE(flux_rs_vf(i)%sf(   -1    : buff_size, &
                                          is2%beg :  is2%end , &
                                          is3%beg :  is3%end  ))
            END DO
            
            ALLOCATE(flux_src_rs_vf(adv_idx%beg)%sf(   -1    : buff_size, &
                                                     is2%beg :  is2%end , &
                                                     is3%beg :  is3%end  ))
            
            IF(riemann_solver == 1) THEN
               DO i = adv_idx%beg+1, adv_idx%end
                  ALLOCATE(flux_src_rs_vf(i)%sf(   -1    : buff_size, &
                                                 is2%beg :  is2%end , &
                                                 is3%beg :  is3%end  ))
               END DO
            ELSE
               DO i = adv_idx%beg+1, adv_idx%end
                  flux_src_rs_vf(i)%sf => flux_src_rs_vf(adv_idx%beg)%sf
               END DO
            END IF
            ! END: Allocation/Association of Primitive and Flux Variables ======
            
            
            ! Reshaping Inputted Data in x-direction ===========================
            IF(cbc_dir == 1) THEN
               
               DO i = 1, sys_size
                  DO r = iz%beg, iz%end
                     DO k = iy%beg, iy%end
                        DO j = 0, buff_size
                           q_prim_rs_vf(i)%sf(j,k,r) = &
                           q_prim_vf(i)%sf(dj*(m-2*j)+j,k,r)
                        END DO
                     END DO
                  END DO
               END DO
               
               DO r = iz%beg, iz%end
                  DO k = iy%beg, iy%end
                     DO j = 0, buff_size
                        q_prim_rs_vf(mom_idx%beg)%sf(j,k,r) = &
                        q_prim_vf(mom_idx%beg)%sf(dj*(m-2*j)+j,k,r) * &
                        SIGN(1d0,-REAL(cbc_loc,KIND(0d0)))
                     END DO
                  END DO
               END DO
               
               DO i = 1, adv_idx%end
                  DO r = iz%beg, iz%end
                     DO k = iy%beg, iy%end
                        DO j = -1, buff_size
                           flux_rs_vf(i)%sf(j,k,r) = &
                           flux_vf(i)%sf(dj*((m-1)-2*j)+j,k,r) * &
                           SIGN(1d0,-REAL(cbc_loc,KIND(0d0)))
                        END DO
                     END DO
                  END DO
               END DO
               
               DO r = iz%beg, iz%end
                  DO k = iy%beg, iy%end
                     DO j = -1, buff_size
                        flux_rs_vf(mom_idx%beg)%sf(j,k,r) = &
                        flux_vf(mom_idx%beg)%sf(dj*((m-1)-2*j)+j,k,r)
                     END DO
                  END DO
               END DO
               
               DO r = iz%beg, iz%end
                  DO k = iy%beg, iy%end
                     DO j = -1, buff_size
                        flux_src_rs_vf(adv_idx%beg)%sf(j,k,r) = &
                        flux_src_vf(adv_idx%beg)%sf(dj*((m-1)-2*j)+j,k,r)
                     END DO
                  END DO
               END DO
               
               IF(riemann_solver == 1) THEN
                  DO i = adv_idx%beg+1, adv_idx%end
                     DO r = iz%beg, iz%end
                        DO k = iy%beg, iy%end
                           DO j = -1, buff_size
                              flux_src_rs_vf(i)%sf(j,k,r) = &
                              flux_src_vf(i)%sf(dj*((m-1)-2*j)+j,k,r)
                           END DO
                        END DO
                     END DO
                  END DO
               ELSE
                  DO r = iz%beg, iz%end
                     DO k = iy%beg, iy%end
                        DO j = -1, buff_size
                           flux_src_rs_vf(adv_idx%beg)%sf(j,k,r) = &
                           flux_src_vf(adv_idx%beg)%sf(dj*((m-1)-2*j)+j,k,r) * &
                           SIGN(1d0,-REAL(cbc_loc,KIND(0d0)))
                        END DO
                     END DO
                  END DO
               END IF
            ! END: Reshaping Inputted Data in x-direction ======================
            
            
            ! Reshaping Inputted Data in y-direction ===========================
            ELSEIF(cbc_dir == 2) THEN
               
               DO i = 1, sys_size
                  DO r = iz%beg, iz%end
                     DO k = ix%beg, ix%end
                        DO j = 0, buff_size
                           q_prim_rs_vf(i)%sf(j,k,r) = &
                           q_prim_vf(i)%sf(k,dj*(n-2*j)+j,r)
                        END DO
                     END DO
                  END DO
               END DO
               
               DO r = iz%beg, iz%end
                  DO k = ix%beg, ix%end
                     DO j = 0, buff_size
                        q_prim_rs_vf(mom_idx%beg+1)%sf(j,k,r) = &
                        q_prim_vf(mom_idx%beg+1)%sf(k,dj*(n-2*j)+j,r) * &
                        SIGN(1d0,-REAL(cbc_loc,KIND(0d0)))
                     END DO
                  END DO
               END DO
               
               DO i = 1, adv_idx%end
                  DO r = iz%beg, iz%end
                     DO k = ix%beg, ix%end
                        DO j = -1, buff_size
                           flux_rs_vf(i)%sf(j,k,r) = &
                           flux_vf(i)%sf(k,dj*((n-1)-2*j)+j,r) * &
                           SIGN(1d0,-REAL(cbc_loc,KIND(0d0)))
                        END DO
                     END DO
                  END DO
               END DO
               
               DO r = iz%beg, iz%end
                  DO k = ix%beg, ix%end
                     DO j = -1, buff_size
                        flux_rs_vf(mom_idx%beg+1)%sf(j,k,r) = &
                        flux_vf(mom_idx%beg+1)%sf(k,dj*((n-1)-2*j)+j,r)
                     END DO
                  END DO
               END DO
               
               DO r = iz%beg, iz%end
                  DO k = ix%beg, ix%end
                     DO j = -1, buff_size
                        flux_src_rs_vf(adv_idx%beg)%sf(j,k,r) = &
                        flux_src_vf(adv_idx%beg)%sf(k,dj*((n-1)-2*j)+j,r)
                     END DO
                  END DO
               END DO
               
               IF(riemann_solver == 1) THEN
                  DO i = adv_idx%beg+1, adv_idx%end
                     DO r = iz%beg, iz%end
                        DO k = ix%beg, ix%end
                           DO j = -1, buff_size
                              flux_src_rs_vf(i)%sf(j,k,r) = &
                              flux_src_vf(i)%sf(k,dj*((n-1)-2*j)+j,r)
                           END DO
                        END DO
                     END DO
                  END DO
               ELSE
                  DO r = iz%beg, iz%end
                     DO k = ix%beg, ix%end
                        DO j = -1, buff_size
                           flux_src_rs_vf(adv_idx%beg)%sf(j,k,r) = &
                           flux_src_vf(adv_idx%beg)%sf(k,dj*((n-1)-2*j)+j,r) * &
                           SIGN(1d0,-REAL(cbc_loc,KIND(0d0)))
                        END DO
                     END DO
                  END DO
               END IF
            ! END: Reshaping Inputted Data in y-direction ======================
            
            
            ! Reshaping Inputted Data in z-direction ===========================
            ELSE
               
               DO i = 1, sys_size
                  DO r = ix%beg, ix%end
                     DO k = iy%beg, iy%end
                        DO j = 0, buff_size
                           q_prim_rs_vf(i)%sf(j,k,r) = &
                           q_prim_vf(i)%sf(r,k,dj*(p-2*j)+j)
                        END DO
                     END DO
                  END DO
               END DO
               
               DO r = ix%beg, ix%end
                  DO k = iy%beg, iy%end
                     DO j = 0, buff_size
                        q_prim_rs_vf(mom_idx%end)%sf(j,k,r) = &
                        q_prim_vf(mom_idx%end)%sf(r,k,dj*(p-2*j)+j) * &
                        SIGN(1d0,-REAL(cbc_loc,KIND(0d0)))
                     END DO
                  END DO
               END DO
               
               DO i = 1, adv_idx%end
                  DO r = ix%beg, ix%end
                     DO k = iy%beg, iy%end
                        DO j = -1, buff_size
                           flux_rs_vf(i)%sf(j,k,r) = &
                           flux_vf(i)%sf(r,k,dj*((p-1)-2*j)+j) * &
                           SIGN(1d0,-REAL(cbc_loc,KIND(0d0)))
                        END DO
                     END DO
                  END DO
               END DO
               
               DO r = ix%beg, ix%end
                  DO k = iy%beg, iy%end
                     DO j = -1, buff_size
                        flux_rs_vf(mom_idx%end)%sf(j,k,r) = &
                        flux_vf(mom_idx%end)%sf(r,k,dj*((p-1)-2*j)+j)
                     END DO
                  END DO
               END DO
               
               DO r = ix%beg, ix%end
                  DO k = iy%beg, iy%end
                     DO j = -1, buff_size
                        flux_src_rs_vf(adv_idx%beg)%sf(j,k,r) = &
                        flux_src_vf(adv_idx%beg)%sf(r,k,dj*((p-1)-2*j)+j)
                     END DO
                  END DO
               END DO
               
               IF(riemann_solver == 1) THEN
                  DO i = adv_idx%beg+1, adv_idx%end
                     DO r = ix%beg, ix%end
                        DO k = iy%beg, iy%end
                           DO j = -1, buff_size
                              flux_src_rs_vf(i)%sf(j,k,r) = &
                              flux_src_vf(i)%sf(r,k,dj*((p-1)-2*j)+j)
                           END DO
                        END DO
                     END DO
                  END DO
               ELSE
                  DO r = ix%beg, ix%end
                     DO k = iy%beg, iy%end
                        DO j = -1, buff_size
                           flux_src_rs_vf(adv_idx%beg)%sf(j,k,r) = &
                           flux_src_vf(adv_idx%beg)%sf(r,k,dj*((p-1)-2*j)+j) * &
                           SIGN(1d0,-REAL(cbc_loc,KIND(0d0)))
                        END DO
                     END DO
                  END DO
               END IF
               
            END IF
            ! END: Reshaping Inputted Data in z-direction ======================
            
            
            ! Association of the procedural pointer to the appropriate procedure
            ! that will be utilized in the evaluation of L variables for the CBC
            
            ! ==================================================================
            IF(     (cbc_dir == 1 .AND. cbc_loc == -1 .AND. bc_x%beg == -5)   &
               .OR. (cbc_dir == 1 .AND. cbc_loc ==  1 .AND. bc_x%end == -5)   &
               .OR. (cbc_dir == 2 .AND. cbc_loc == -1 .AND. bc_y%beg == -5)   &
               .OR. (cbc_dir == 2 .AND. cbc_loc ==  1 .AND. bc_y%end == -5)   &
               .OR. (cbc_dir == 3 .AND. cbc_loc == -1 .AND. bc_z%beg == -5)   &
               .OR. (cbc_dir == 3 .AND. cbc_loc ==  1 .AND. bc_z%end == -5) ) &
                                             THEN
               
               s_compute_L => s_compute_slip_wall_L
               
            ELSEIF( (cbc_dir == 1 .AND. cbc_loc == -1 .AND. bc_x%beg == -6)   &
               .OR. (cbc_dir == 1 .AND. cbc_loc ==  1 .AND. bc_x%end == -6)   &
               .OR. (cbc_dir == 2 .AND. cbc_loc == -1 .AND. bc_y%beg == -6)   &
               .OR. (cbc_dir == 2 .AND. cbc_loc ==  1 .AND. bc_y%end == -6)   &
               .OR. (cbc_dir == 3 .AND. cbc_loc == -1 .AND. bc_z%beg == -6)   &
               .OR. (cbc_dir == 3 .AND. cbc_loc ==  1 .AND. bc_z%end == -6) ) &
                                             THEN
               
               s_compute_L => s_compute_nonreflecting_subsonic_buffer_L
               
            ELSEIF( (cbc_dir == 1 .AND. cbc_loc == -1 .AND. bc_x%beg == -7)   &
               .OR. (cbc_dir == 1 .AND. cbc_loc ==  1 .AND. bc_x%end == -7)   &
               .OR. (cbc_dir == 2 .AND. cbc_loc == -1 .AND. bc_y%beg == -7)   &
               .OR. (cbc_dir == 2 .AND. cbc_loc ==  1 .AND. bc_y%end == -7)   &
               .OR. (cbc_dir == 3 .AND. cbc_loc == -1 .AND. bc_z%beg == -7)   &
               .OR. (cbc_dir == 3 .AND. cbc_loc ==  1 .AND. bc_z%end == -7) ) &
                                             THEN
               
               s_compute_L => s_compute_nonreflecting_subsonic_inflow_L
               
            ELSEIF( (cbc_dir == 1 .AND. cbc_loc == -1 .AND. bc_x%beg == -8)   &
               .OR. (cbc_dir == 1 .AND. cbc_loc ==  1 .AND. bc_x%end == -8)   &
               .OR. (cbc_dir == 2 .AND. cbc_loc == -1 .AND. bc_y%beg == -8)   &
               .OR. (cbc_dir == 2 .AND. cbc_loc ==  1 .AND. bc_y%end == -8)   &
               .OR. (cbc_dir == 3 .AND. cbc_loc == -1 .AND. bc_z%beg == -8)   &
               .OR. (cbc_dir == 3 .AND. cbc_loc ==  1 .AND. bc_z%end == -8) ) &
                                             THEN
               
               s_compute_L => s_compute_nonreflecting_subsonic_outflow_L
               
            ELSEIF( (cbc_dir == 1 .AND. cbc_loc == -1 .AND. bc_x%beg == -9)   &
               .OR. (cbc_dir == 1 .AND. cbc_loc ==  1 .AND. bc_x%end == -9)   &
               .OR. (cbc_dir == 2 .AND. cbc_loc == -1 .AND. bc_y%beg == -9)   &
               .OR. (cbc_dir == 2 .AND. cbc_loc ==  1 .AND. bc_y%end == -9)   &
               .OR. (cbc_dir == 3 .AND. cbc_loc == -1 .AND. bc_z%beg == -9)   &
               .OR. (cbc_dir == 3 .AND. cbc_loc ==  1 .AND. bc_z%end == -9) ) &
                                             THEN
               
               s_compute_L => s_compute_force_free_subsonic_outflow_L
               
            ELSEIF( (cbc_dir == 1 .AND. cbc_loc == -1 .AND. bc_x%beg ==-10)   &
               .OR. (cbc_dir == 1 .AND. cbc_loc ==  1 .AND. bc_x%end ==-10)   &
               .OR. (cbc_dir == 2 .AND. cbc_loc == -1 .AND. bc_y%beg ==-10)   &
               .OR. (cbc_dir == 2 .AND. cbc_loc ==  1 .AND. bc_y%end ==-10)   &
               .OR. (cbc_dir == 3 .AND. cbc_loc == -1 .AND. bc_z%beg ==-10)   &
               .OR. (cbc_dir == 3 .AND. cbc_loc ==  1 .AND. bc_z%end ==-10) ) &
                                             THEN
               
               s_compute_L => s_compute_constant_pressure_subsonic_outflow_L
               
            ELSEIF( (cbc_dir == 1 .AND. cbc_loc == -1 .AND. bc_x%beg ==-11)   &
               .OR. (cbc_dir == 1 .AND. cbc_loc ==  1 .AND. bc_x%end ==-11)   &
               .OR. (cbc_dir == 2 .AND. cbc_loc == -1 .AND. bc_y%beg ==-11)   &
               .OR. (cbc_dir == 2 .AND. cbc_loc ==  1 .AND. bc_y%end ==-11)   &
               .OR. (cbc_dir == 3 .AND. cbc_loc == -1 .AND. bc_z%beg ==-11)   &
               .OR. (cbc_dir == 3 .AND. cbc_loc ==  1 .AND. bc_z%end ==-11) ) &
                                             THEN
               
               s_compute_L => s_compute_supersonic_inflow_L
               
            ELSE
               
               s_compute_L => s_compute_supersonic_outflow_L
               
            END IF
            ! ==================================================================
            
            
        END SUBROUTINE s_initialize_cbc ! --------------------------------------
        
        
        
        !>  Deallocation and/or the disassociation procedures that
        !!      are necessary in order to finalize the CBC application
        !!  @param flux_vf Cell-boundary-average fluxes
        !!  @param flux_src_vf Cell-boundary-average flux sources
        !!  @param cbc_dir CBC coordinate direction
        !!  @param cbc_loc CBC coordinate location
        !!  @param ix Index bound in the first coordinate direction 
        !!  @param iy Index bound in the second coordinate direction 
        !!  @param iz Index bound in the third coordinate direction         
        SUBROUTINE s_finalize_cbc( flux_vf, flux_src_vf, & ! -------------------
                                       cbc_dir, cbc_loc, &
                                               ix,iy,iz  )
           
            TYPE(scalar_field),  &
            DIMENSION(sys_size), &
            INTENT(INOUT) :: flux_vf, flux_src_vf
            
            INTEGER, INTENT(IN) :: cbc_dir, cbc_loc
            TYPE(bounds_info), INTENT(IN) :: ix,iy,iz
            

            INTEGER :: dj !< Indical shift based on CBC location
            

            INTEGER :: i,j,k,r !< Generic loop iterators
            
            
            ! Determining the indicial shift based on CBC location
            dj = MAX(0,cbc_loc)
            
            
            ! Reshaping Outputted Data in x-direction ==========================
            IF(cbc_dir == 1) THEN
               
               DO i = 1, adv_idx%end
                  DO r = iz%beg, iz%end
                     DO k = iy%beg, iy%end
                        DO j = -1, buff_size
                           flux_vf(i)%sf(dj*((m-1)-2*j)+j,k,r) = &
                           flux_rs_vf(i)%sf(j,k,r) * &
                           SIGN(1d0,-REAL(cbc_loc,KIND(0d0)))
                        END DO
                     END DO
                  END DO
               END DO
               
               DO r = iz%beg, iz%end
                  DO k = iy%beg, iy%end
                     DO j = -1, buff_size
                        flux_vf(mom_idx%beg)%sf(dj*((m-1)-2*j)+j,k,r) = &
                        flux_rs_vf(mom_idx%beg)%sf(j,k,r)
                     END DO
                  END DO
               END DO
               
               DO r = iz%beg, iz%end
                  DO k = iy%beg, iy%end
                     DO j = -1, buff_size
                        flux_src_vf(adv_idx%beg)%sf(dj*((m-1)-2*j)+j,k,r) = &
                        flux_src_rs_vf(adv_idx%beg)%sf(j,k,r)
                     END DO
                  END DO
               END DO
               
               IF(riemann_solver == 1) THEN
                  DO i = adv_idx%beg+1, adv_idx%end
                     DO r = iz%beg, iz%end
                        DO k = iy%beg, iy%end
                           DO j = -1, buff_size
                              flux_src_vf(i)%sf(dj*((m-1)-2*j)+j,k,r) = &
                              flux_src_rs_vf(i)%sf(j,k,r)
                           END DO
                        END DO
                     END DO
                  END DO
               ELSE
                  DO r = iz%beg, iz%end
                     DO k = iy%beg, iy%end
                        DO j = -1, buff_size
                           flux_src_vf(adv_idx%beg)%sf(dj*((m-1)-2*j)+j,k,r) = &
                           flux_src_rs_vf(adv_idx%beg)%sf(j,k,r) * &
                           SIGN(1d0,-REAL(cbc_loc,KIND(0d0)))
                        END DO
                     END DO
                  END DO
               END IF
            ! END: Reshaping Outputted Data in x-direction =====================
            
            
            ! Reshaping Outputted Data in y-direction ==========================
            ELSEIF(cbc_dir == 2) THEN
               
               DO i = 1, adv_idx%end
                  DO r = iz%beg, iz%end
                     DO k = ix%beg, ix%end
                        DO j = -1, buff_size
                           flux_vf(i)%sf(k,dj*((n-1)-2*j)+j,r) = &
                           flux_rs_vf(i)%sf(j,k,r) * &
                           SIGN(1d0,-REAL(cbc_loc,KIND(0d0)))
                        END DO
                     END DO
                  END DO
               END DO
               
               DO r = iz%beg, iz%end
                  DO k = ix%beg, ix%end
                     DO j = -1, buff_size
                        flux_vf(mom_idx%beg+1)%sf(k,dj*((n-1)-2*j)+j,r) = &
                        flux_rs_vf(mom_idx%beg+1)%sf(j,k,r)
                     END DO
                  END DO
               END DO
               
               DO r = iz%beg, iz%end
                  DO k = ix%beg, ix%end
                     DO j = -1, buff_size
                        flux_src_vf(adv_idx%beg)%sf(k,dj*((n-1)-2*j)+j,r) = &
                        flux_src_rs_vf(adv_idx%beg)%sf(j,k,r)
                     END DO
                  END DO
               END DO
               
               IF(riemann_solver == 1) THEN
                  DO i = adv_idx%beg+1, adv_idx%end
                     DO r = iz%beg, iz%end
                        DO k = ix%beg, ix%end
                           DO j = -1, buff_size
                              flux_src_vf(i)%sf(k,dj*((n-1)-2*j)+j,r) = &
                              flux_src_rs_vf(i)%sf(j,k,r)
                           END DO
                        END DO
                     END DO
                  END DO
               ELSE
                  DO r = iz%beg, iz%end
                     DO k = ix%beg, ix%end
                        DO j = -1, buff_size
                           flux_src_vf(adv_idx%beg)%sf(k,dj*((n-1)-2*j)+j,r) = &
                           flux_src_rs_vf(adv_idx%beg)%sf(j,k,r) * &
                           SIGN(1d0,-REAL(cbc_loc,KIND(0d0)))
                        END DO
                     END DO
                  END DO
               END IF
            ! END: Reshaping Outputted Data in y-direction =====================
            
            
            ! Reshaping Outputted Data in z-direction ==========================
            ELSE
               
               DO i = 1, adv_idx%end
                  DO r = ix%beg, ix%end
                     DO k = iy%beg, iy%end
                        DO j = -1, buff_size
                           flux_vf(i)%sf(r,k,dj*((p-1)-2*j)+j) = &
                           flux_rs_vf(i)%sf(j,k,r) * &
                           SIGN(1d0,-REAL(cbc_loc,KIND(0d0)))
                        END DO
                     END DO
                  END DO
               END DO
               
               DO r = ix%beg, ix%end
                  DO k = iy%beg, iy%end
                     DO j = -1, buff_size
                        flux_vf(mom_idx%end)%sf(r,k,dj*((p-1)-2*j)+j) = &
                        flux_rs_vf(mom_idx%end)%sf(j,k,r)
                     END DO
                  END DO
               END DO
               
               DO r = ix%beg, ix%end
                  DO k = iy%beg, iy%end
                     DO j = -1, buff_size
                        flux_src_vf(adv_idx%beg)%sf(r,k,dj*((p-1)-2*j)+j) = &
                        flux_src_rs_vf(adv_idx%beg)%sf(j,k,r)
                     END DO
                  END DO
               END DO
               
               IF(riemann_solver == 1) THEN
                  DO i = adv_idx%beg+1, adv_idx%end
                     DO r = ix%beg, ix%end
                        DO k = iy%beg, iy%end
                           DO j = -1, buff_size
                              flux_src_vf(i)%sf(r,k,dj*((p-1)-2*j)+j) = &
                              flux_src_rs_vf(i)%sf(j,k,r)
                           END DO
                        END DO
                     END DO
                  END DO
               ELSE
                  DO r = ix%beg, ix%end
                     DO k = iy%beg, iy%end
                        DO j = -1, buff_size
                           flux_src_vf(adv_idx%beg)%sf(r,k,dj*((p-1)-2*j)+j) = &
                           flux_src_rs_vf(adv_idx%beg)%sf(j,k,r) * &
                           SIGN(1d0,-REAL(cbc_loc,KIND(0d0)))
                        END DO
                     END DO
                  END DO
               END IF
               
            END IF
            ! END: Reshaping Outputted Data in z-direction =====================
            
            
            ! Deallocation/Disassociation of Primitive and Flux Variables ======
            DO i = 1, sys_size
               DEALLOCATE(q_prim_rs_vf(i)%sf)
            END DO
            
            IF(weno_order > 1) THEN
               
               DO i = 1, adv_idx%end
                  DEALLOCATE(F_rs_vf(i)%sf)
               END DO
               
               DEALLOCATE(F_src_rs_vf(adv_idx%beg)%sf)
               
               IF(riemann_solver == 1) THEN
                  DO i = adv_idx%beg+1, adv_idx%end
                     DEALLOCATE(F_src_rs_vf(i)%sf)
                  END DO
               ELSE
                  DO i = adv_idx%beg+1, adv_idx%end
                     NULLIFY(F_src_rs_vf(i)%sf)
                  END DO
               END IF
               
            END IF
            
            DO i = 1, adv_idx%end
               DEALLOCATE(flux_rs_vf(i)%sf)
            END DO
            
            DEALLOCATE(flux_src_rs_vf(adv_idx%beg)%sf)
            
            IF(riemann_solver == 1) THEN
               DO i = adv_idx%beg+1, adv_idx%end
                  DEALLOCATE(flux_src_rs_vf(i)%sf)
               END DO
            ELSE
               DO i = adv_idx%beg+1, adv_idx%end
                  NULLIFY(flux_src_rs_vf(i)%sf)
               END DO
            END IF
            ! ==================================================================
            
            
            ! Nullifying procedural pointer used in evaluation of L for the CBC
            s_compute_L => NULL()
            
            
        END SUBROUTINE s_finalize_cbc ! ----------------------------------------
        
        
        
        
        
        !> Module deallocation and/or disassociation procedures
        SUBROUTINE s_finalize_cbc_module() ! -----------------------------------
            
            
            IF(              ALL((/bc_x%beg,bc_x%end/) > -5)  &
                                    .AND.                     &
                (n > 0 .AND. ALL((/bc_y%beg,bc_y%end/) > -5)) &
                                    .AND.                     &
                (p > 0 .AND. ALL((/bc_z%beg,bc_z%end/) > -5)) ) RETURN
            
            
            ! Deallocating the cell-average primitive variables
            DEALLOCATE(q_prim_rs_vf)
            
            ! Deallocating the cell-average and cell-boundary-average fluxes
            DEALLOCATE(   F_rs_vf,    F_src_rs_vf)
            DEALLOCATE(flux_rs_vf, flux_src_rs_vf)
            
            ! Deallocating the cell-average partial densities, the velocity, the
            ! advection variables, the mass fractions and also the Weber numbers
            DEALLOCATE(alpha_rho, vel, adv, mf, We)
            
            ! Deallocating the first-order spatial derivatives, in s-direction,
            ! of the partial densities, the velocity and the advected variables
            DEALLOCATE(dalpha_rho_ds, dvel_ds, dadv_ds)
            
            ! Deallocating L, see Thompson (1987, 1990)
            DEALLOCATE(L)
            
            ! Deallocating the cell-width distribution in the s-direction
            DEALLOCATE(ds)
            
            
            ! Deallocating CBC Coefficients in x-direction =====================
            IF(ANY((/bc_x%beg,bc_x%end/) <= -5)) THEN
                DEALLOCATE(fd_coef_x); IF(weno_order > 1) DEALLOCATE(pi_coef_x)
            END IF
            ! ==================================================================
            
            
            ! Deallocating CBC Coefficients in y-direction =====================
            IF(n > 0 .AND. ANY((/bc_y%beg,bc_y%end/) <= -5)) THEN
                DEALLOCATE(fd_coef_y); IF(weno_order > 1) DEALLOCATE(pi_coef_y)
            END IF
            ! ==================================================================
            
            
            ! Deallocating CBC Coefficients in z-direction =====================
            IF(p > 0 .AND. ANY((/bc_z%beg,bc_z%end/) <= -5)) THEN
                DEALLOCATE(fd_coef_z); IF(weno_order > 1) DEALLOCATE(pi_coef_z)
            END IF
            ! ==================================================================
            
            
            ! Disassociating the pointer to the procedure that was utilized to
            ! to convert mixture or species variables to the mixture variables
            s_convert_to_mixture_variables => NULL()
            
            
        END SUBROUTINE s_finalize_cbc_module ! ---------------------------------
        
        
        
        
        
END MODULE m_cbc
