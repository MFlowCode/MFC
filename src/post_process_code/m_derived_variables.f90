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
!! @file m_derived_variables.f90
!! @brief Contains module m_derived_variables
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief This module features subroutines that allow for the derivation of
!!      numerous flow variables from the conservative and primitive ones.
!!      Currently, the available derived variables include the unadvected
!!      volume fraction, specific heat ratio, liquid stiffness, speed of
!!      sound, vorticity and the numerical Schlieren function.
MODULE m_derived_variables
    
    
    ! Dependencies =============================================================
    USE m_derived_types         !< Definitions of the derived types
    
    USE m_global_parameters     !< Global parameters for the code
    
    USE m_mpi_proxy             !< Message passing interface (MPI) module proxy
    ! ==========================================================================
    
    
    IMPLICIT NONE
    
    PRIVATE; PUBLIC :: s_initialize_derived_variables_module, &
                       s_compute_finite_difference_coefficients, &
                       s_derive_unadvected_volume_fraction, &
                       s_derive_specific_heat_ratio, &
                       s_derive_liquid_stiffness, &
                       s_derive_sound_speed, &
                       s_derive_flux_limiter, &
                       s_derive_vorticity_component, &
                       s_derive_numerical_schlieren_function, &
                       s_derive_curvature, &
                       s_finalize_derived_variables_module
    
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:,:) :: gm_rho_sf !<
    !! Gradient magnitude (gm) of the density for each cell of the computational
    !! sub-domain. This variable is employed in the calculation of the numerical
    !! Schlieren function.
    
    !> @name Finite-difference (fd) coefficients in x-, y- and z-coordinate directions.
    !! Note that because sufficient boundary information is available for all the
    !! active coordinate directions, the centered family of the finite-difference
    !! schemes is used.
    !> @{
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:), PUBLIC :: fd_coeff_x
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:), PUBLIC :: fd_coeff_y
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:), PUBLIC :: fd_coeff_z
    !> @}

    INTEGER, PRIVATE :: flg  !<
    !! Flagging (flg) variable used to annotate the dimensionality of the dataset
    !! that is undergoing the post-process. A flag value of 1 indicates that the
    !! dataset is 3D, while a flag value of 0 indicates that it is not. This flg
    !! variable is necessary to avoid cycling through the third dimension of the
    !! flow variable(s) when the simulation is not 3D and the size of the buffer
    !! is non-zero. Note that a similar procedure does not have to be applied to
    !! the second dimension since in 1D, the buffer size is always zero.
    
    
    CONTAINS
        
        
        
        
        !>  Computation of parameters, allocation procedures, and/or
        !!      any other tasks needed to properly setup the module        
        SUBROUTINE s_initialize_derived_variables_module() ! ----------------------
            
            ! Allocating the gradient magnitude of the density variable provided
            ! that numerical Schlieren function is outputted during post-process
            IF(schlieren_wrt) THEN
                ALLOCATE(gm_rho_sf( -offset_x%beg : m + offset_x%end, &
                                    -offset_y%beg : n + offset_y%end, &
                                    -offset_z%beg : p + offset_z%end ))
            END IF
            
            
            ! Allocating the variables which will store the coefficients of the
            ! centered family of finite-difference schemes. Note that sufficient
            ! space is allocated so that the coefficients up to any chosen order
            ! of accuracy may be bookkept. However, if higher than fourth-order
            ! accuracy coefficients are wanted, the formulae required to compute
            ! these coefficients will have to be implemented in the subroutine
            ! s_compute_finite_difference_coefficients.
            
            ! Allocating centered finite-difference coefficients in x-direction
            IF(omega_wrt(2) .OR. omega_wrt(3) .OR. schlieren_wrt) THEN
               ALLOCATE(fd_coeff_x(    -fd_number : fd_number      , &
                                    -offset_x%beg : m+offset_x%end ) )
            END IF
            
            ! Allocating centered finite-difference coefficients in y-direction
            IF(               omega_wrt(1) .OR. omega_wrt(3)              &
                                           .OR.                           &
                               (n > 0 .AND. schlieren_wrt)                ) THEN
               ALLOCATE(fd_coeff_y(    -fd_number : fd_number      , &
                                    -offset_y%beg : n+offset_y%end ) )
            END IF
            
            ! Allocating centered finite-difference coefficients in z-direction
            IF(               omega_wrt(1) .OR. omega_wrt(2)              &
                                           .OR.                           &
                               (p > 0 .AND. schlieren_wrt)                ) THEN
               ALLOCATE(fd_coeff_z(    -fd_number : fd_number      , &
                                    -offset_z%beg : p+offset_z%end ) )
            END IF
            
            
            ! Annotating the dimensionality of the dataset undergoing the post-
            ! process. A flag value of 1 indicates that the dataset is 3D, while
            ! a flag value of 0 indicates that it is not.
            IF(p > 0) THEN
                flg = 1
            ELSE
                flg = 0
            END IF
            
            
        END SUBROUTINE s_initialize_derived_variables_module ! --------------------
        
        
        
        !> @name The purpose of this subroutine is to compute the finite-
        !!      difference coefficients for the centered schemes utilized
        !!      in computations of first order spatial derivatives in the
        !!      s-coordinate direction. The s-coordinate direction refers
        !!      to the x-, y- or z-coordinate direction, depending on the
        !!      subroutine's inputs. Note that coefficients of up to 4th
        !!      order accuracy are available.
        !!  @param q Number of cells in the s-coordinate direction
        !!  @param offset_s  Size of the ghost zone layer in the s-coordinate direction
        !!  @param s_cc Locations of the cell-centers in the s-coordinate direction
        !!  @param fd_coeff_s Finite-diff. coefficients in the s-coordinate direction 
        SUBROUTINE s_compute_finite_difference_coefficients(    q, offset_s,   &
                                                             s_cc, fd_coeff_s  )
           
            INTEGER, INTENT(IN) :: q
            TYPE(bounds_info), INTENT(IN) :: offset_s
            
            REAL(KIND(0d0)), &
            DIMENSION(-buff_size:q+buff_size), &
            INTENT(IN) :: s_cc
            
            REAL(KIND(0d0)), &
            DIMENSION(-fd_number:fd_number, -offset_s%beg:q+offset_s%end), &
            INTENT(INOUT) :: fd_coeff_s
            

            INTEGER :: i !< Generic loop iterator
            
            ! Computing the 1st order finite-difference coefficients
            IF(fd_order == 1) THEN
                DO i = -offset_s%beg, q+offset_s%end
                    fd_coeff_s(-1,i) =  0d0
                    fd_coeff_s( 0,i) = -1d0 / (s_cc(i+1) - s_cc(i))
                    fd_coeff_s( 1,i) = -fd_coeff_s(0,i)
                END DO
                
            ! Computing the 2nd order finite-difference coefficients
            ELSEIF(fd_order == 2) THEN
                DO i = -offset_s%beg, q+offset_s%end
                    fd_coeff_s(-1,i) = -1d0 / (s_cc(i+1) - s_cc(i-1))
                    fd_coeff_s( 0,i) =  0d0
                    fd_coeff_s( 1,i) = -fd_coeff_s(-1,i)
                END DO
                
            ! Computing the 4th order finite-difference coefficients
            ELSE
                DO i = -offset_s%beg, q+offset_s%end
                    fd_coeff_s(-2,i) =  1d0 / ( s_cc(i-2) - 8d0*s_cc(i-1) &
                                              - s_cc(i+2) + 8d0*s_cc(i+1) )
                    fd_coeff_s(-1,i) = -8d0*fd_coeff_s(-2,i)
                    fd_coeff_s( 0,i) =  0d0
                    fd_coeff_s( 1,i) = -fd_coeff_s(-1,i)
                    fd_coeff_s( 2,i) = -fd_coeff_s(-2,i)
                END DO
                
            END IF
            
            
        END SUBROUTINE s_compute_finite_difference_coefficients ! --------------
        
        
        
        !>  This subroutine is used together with the volume fraction
        !!      model and when called upon, it computes the values of the
        !!      unadvected volume fraction from the inputted conservative
        !!      variables, q_cons_vf. The calculated values are stored in
        !!      the derived flow quantity storage variable, q_sf.        
        !!  @param q_cons_vf Conservative variables
        !!  @param q_sf Unadvected volume fraction
        SUBROUTINE s_derive_unadvected_volume_fraction(q_cons_vf, q_sf) ! ------

            TYPE(scalar_field), &
            DIMENSION(sys_size), &
            INTENT(IN) :: q_cons_vf
            
            REAL(KIND(0d0)), &
            DIMENSION( -offset_x%beg : m+offset_x%end  , &
                       -offset_y%beg : n+offset_y%end  , &
                       -offset_z%beg : p+offset_z%end ), &
            INTENT(INOUT) :: q_sf
            

            INTEGER :: i,j,k,l !< Generic loop iterators
            
            ! Computing unadvected volume fraction from conservative variables
            DO l = -offset_z%beg, p + offset_z%end
                DO k = -offset_y%beg, n + offset_y%end
                    DO j = -offset_x%beg, m + offset_x%end
                        
                        q_sf(j,k,l) = 1d0
                        
                        DO i = adv_idx%beg, adv_idx%end
                            q_sf(j,k,l) = q_sf(j,k,l) - q_cons_vf(i)%sf(j,k,l)
                        END DO
                        
                    END DO
                END DO
            END DO
            
            
        END SUBROUTINE s_derive_unadvected_volume_fraction ! -------------------
        
        
        
        
        !>  This subroutine receives as input the specific heat ratio
        !!      function, gamma_sf, and derives from it the specific heat
        !!      ratio. The latter is stored in the derived flow quantity
        !!      storage variable, q_sf.
        !!  @param gamma_sf Specific heat ratio function
        !!  @param q_sf Specific heat ratio
        SUBROUTINE s_derive_specific_heat_ratio(gamma_sf, q_sf) ! --------------
            
            REAL(KIND(0d0)), &
            DIMENSION( -buff_size     :  m+buff_size       , &
                       -buff_size     :  n+buff_size       , &
                       -buff_size*flg : (p+buff_size)*flg ), &
            INTENT(IN) :: gamma_sf
            
            REAL(KIND(0d0)), &
            DIMENSION( -offset_x%beg : m+offset_x%end  , &
                       -offset_y%beg : n+offset_y%end  , &
                       -offset_z%beg : p+offset_z%end ), &
            INTENT(INOUT) :: q_sf
            

            INTEGER :: i,j,k !< Generic loop iterators
            
            
            ! Computing specific heat ratio from specific heat ratio function
            DO k = -offset_z%beg, p + offset_z%end
                DO j = -offset_y%beg, n + offset_y%end
                    DO i = -offset_x%beg, m + offset_x%end
                        q_sf(i,j,k) = 1d0 + 1d0 / gamma_sf(i,j,k)
                    END DO
                END DO
            END DO
            
            
        END SUBROUTINE s_derive_specific_heat_ratio ! --------------------------
        
        
        
        !>  This subroutine admits as inputs the specific heat ratio
        !!      function and the liquid stiffness function, gamma_sf and
        !!      pi_inf_sf, respectively. These are used to calculate the
        !!      values of the liquid stiffness, which are stored in the
        !!      derived flow quantity storage variable, q_sf.        
        !!  @param gamma_sf Specific heat ratio
        !!  @param pi_inf_sf Liquid stiffness function
        !!  @param q_sf Liquid stiffness
        SUBROUTINE s_derive_liquid_stiffness(gamma_sf, pi_inf_sf, q_sf) ! ------

            REAL(KIND(0d0)), &
            DIMENSION( -buff_size     :  m+buff_size       , &
                       -buff_size     :  n+buff_size       , &
                       -buff_size*flg : (p+buff_size)*flg ), &
            INTENT(IN) :: gamma_sf, pi_inf_sf
            
            REAL(KIND(0d0)), &
            DIMENSION( -offset_x%beg : m+offset_x%end  , &
                       -offset_y%beg : n+offset_y%end  , &
                       -offset_z%beg : p+offset_z%end ), &
            INTENT(INOUT) :: q_sf
            

            INTEGER :: i,j,k !< Generic loop iterators
            
            
            ! Calculating the values of the liquid stiffness from those of the
            ! specific heat ratio function and the liquid stiffness function
            DO k = -offset_z%beg, p + offset_z%end
                DO j = -offset_y%beg, n + offset_y%end
                    DO i = -offset_x%beg, m + offset_x%end
                        q_sf(i,j,k) = pi_inf_sf(i,j,k) / (gamma_sf(i,j,k) + 1d0)
                    END DO
                END DO
            END DO
            
            
        END SUBROUTINE s_derive_liquid_stiffness ! -----------------------------
        
        
        
        !> This subroutine admits as inputs the primitive variables,
        !!      the density, the specific heat ratio function and liquid
        !!      stiffness function. It then computes from those variables
        !!      the values of the speed of sound, which are stored in the
        !!      derived flow quantity storage variable, q_sf.        
        !! @param q_prim_vf Primitive variables
        !! @param rho_sf Density
        !! @param gamma_sf Specific heat ratio function
        !! @param pi_inf_sf Liquid stiffness function
        !! @param q_sf Speed of sound
        SUBROUTINE s_derive_sound_speed( q_prim_vf, rho_sf, gamma_sf, & ! ------
                                                pi_inf_sf, q_sf       )

            TYPE(scalar_field), &
            DIMENSION(sys_size), &
            INTENT(IN) :: q_prim_vf
            
            REAL(KIND(0d0)), &
            DIMENSION( -buff_size     :  m+buff_size       , &
                       -buff_size     :  n+buff_size       , &
                       -buff_size*flg : (p+buff_size)*flg ), &
            INTENT(IN) :: rho_sf, gamma_sf, pi_inf_sf
            
            REAL(KIND(0d0)), &
            DIMENSION( -offset_x%beg : m+offset_x%end  , &
                       -offset_y%beg : n+offset_y%end  , &
                       -offset_z%beg : p+offset_z%end ), &
            INTENT(INOUT) :: q_sf
            

            INTEGER :: i,j,k !< Generic loop iterators

            ! Fluid bulk modulus for alternate sound speed
            REAL(KIND(0d0)) :: blkmod1, blkmod2
            
            ! Computing speed of sound values from those of pressure, density,
            ! specific heat ratio function and the liquid stiffness function
            DO k = -offset_z%beg, p + offset_z%end
                DO j = -offset_y%beg, n + offset_y%end
                    DO i = -offset_x%beg, m + offset_x%end

                        ! Compute mixture sound speed
                        IF (alt_soundspeed .NEQV. .TRUE.) THEN
                            q_sf(i,j,k) = (((gamma_sf(i,j,k) + 1d0) * &
                                q_prim_vf(E_idx)%sf(i,j,k) + &
                                pi_inf_sf(i,j,k)) / (gamma_sf(i,j,k) * &
                                rho_sf(i,j,k)))
                        ELSE
                            blkmod1 = ((fluid_pp(1)%gamma + 1d0)*q_prim_vf(E_idx)%sf(i,j,k) + &
                                fluid_pp(1)%pi_inf)/fluid_pp(1)%gamma
                            blkmod2 = ((fluid_pp(2)%gamma + 1d0)*q_prim_vf(E_idx)%sf(i,j,k) + &
                                fluid_pp(2)%pi_inf)/fluid_pp(2)%gamma
                            q_sf(i,j,k) = (1d0/(rho_sf(i,j,k)*(q_prim_vf(adv_idx%beg)%sf(i,j,k)/blkmod1 + &
                                    (1d0-q_prim_vf(adv_idx%beg)%sf(i,j,k))/blkmod2)))
                        END IF

                        IF (mixture_err .AND. q_sf(i,j,k) < 0d0) THEN
                            q_sf(i,j,k) = 1d-16
                        ELSE
                            q_sf(i,j,k) = SQRT(q_sf(i,j,k))
                        END IF
                    END DO
                END DO
            END DO
            
            
        END SUBROUTINE s_derive_sound_speed ! ----------------------------------
        
        
        
        
        !>  This subroutine derives the flux_limiter at cell boundary
        !!      i+1/2. This is an approximation because the velocity used
        !!      to determine the upwind direction is the velocity at the
        !!      cell center i instead of the contact velocity at the cell
        !!      boundary from the Riemann solver.        
        !!  @param i Component indicator
        !!  @param q_prim_vf Primitive variables
        !!  @param q_sf Flux limiter
        SUBROUTINE s_derive_flux_limiter(i, q_prim_vf, q_sf) ! -----------------

            INTEGER, INTENT(IN) :: i
    
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q_prim_vf
    
            REAL(KIND(0d0)), DIMENSION (    -offset_x%beg : m+offset_x%end , &
                                            -offset_y%beg : n+offset_y%end , &
                                            -offset_z%beg : p+offset_z%end), &
            INTENT(INOUT) :: q_sf

            REAL(KIND(0d0)) :: top, bottom, slope !< Flux limiter calcs
            INTEGER :: j,k,l !< Generic loop iterators
    
            DO l = -offset_z%beg, p + offset_z%end
                DO k = -offset_y%beg, n + offset_y%end
                    DO j = -offset_x%beg, m + offset_x%end
                        IF (i == 1) THEN
                            IF (q_prim_vf(cont_idx%end+i)%sf(j,k,l) >= 0d0) THEN
                                top =   q_prim_vf(adv_idx%beg)%sf( j ,k,l) - & 
                                    q_prim_vf(adv_idx%beg)%sf(j-1,k,l)
                                bottom= q_prim_vf(adv_idx%beg)%sf(j+1,k,l) - &
                                    q_prim_vf(adv_idx%beg)%sf( j ,k,l)
                            ELSE
                                top =   q_prim_vf(adv_idx%beg)%sf(j+2,k,l) - & 
                                    q_prim_vf(adv_idx%beg)%sf(j+1,k,l)
                                bottom= q_prim_vf(adv_idx%beg)%sf(j+1,k,l) - &
                                    q_prim_vf(adv_idx%beg)%sf( j ,k,l)
                            END IF
                        ELSEIF (i == 2) THEN
                            IF (q_prim_vf(cont_idx%end+i)%sf(j,k,l) >= 0d0) THEN
                                top =   q_prim_vf(adv_idx%beg)%sf(j, k ,l) - & 
                                    q_prim_vf(adv_idx%beg)%sf(j,k-1,l)
                                bottom= q_prim_vf(adv_idx%beg)%sf(j,k+1,l) - &
                                    q_prim_vf(adv_idx%beg)%sf(j, k ,l)
                            ELSE
                                top =   q_prim_vf(adv_idx%beg)%sf(j,k+2,l) - & 
                                    q_prim_vf(adv_idx%beg)%sf(j,k+1,l)
                                bottom= q_prim_vf(adv_idx%beg)%sf(j,k+1,l) - &
                                    q_prim_vf(adv_idx%beg)%sf(j, k ,l)
                            END IF
                        ELSE
                            IF (q_prim_vf(cont_idx%end+i)%sf(j,k,l) >= 0d0) THEN
                                top =   q_prim_vf(adv_idx%beg)%sf(j,k, l ) - & 
                                    q_prim_vf(adv_idx%beg)%sf(j,k,l-1)
                                bottom= q_prim_vf(adv_idx%beg)%sf(j,k,l+1) - &
                                    q_prim_vf(adv_idx%beg)%sf(j,k, l )
                            ELSE
                                top =   q_prim_vf(adv_idx%beg)%sf(j,k,l+2) - & 
                                    q_prim_vf(adv_idx%beg)%sf(j,k,l+1)
                                bottom= q_prim_vf(adv_idx%beg)%sf(j,k,l+1) - &
                                    q_prim_vf(adv_idx%beg)%sf(j,k, l )
                            END IF
                        END IF
                
                        IF (ABS(top) < 1d-8) top = 0d0
                        IF (ABS(bottom) < 1d-8) bottom = 0d0
                
                        IF (top == bottom) THEN
                            slope = 1d0
                !       ELSEIF((top == 0d0 .AND. bottom /= 0d0) &
                !               .OR.            &
                !           (bottom == 0d0 .AND. top /= 0d0)) THEN
                !           slope = 0d0
                        ELSE
                            slope = (top*bottom)/(bottom**2d0+1d-16)
                        END IF
                
                        ! Flux limiter function
                        IF (flux_lim == 1) THEN ! MINMOD (MM)
                            q_sf(j,k,l) = MAX(0d0,MIN(1d0,slope))
                        ELSEIF (flux_lim == 2) THEN ! MUSCL (MC)
                            q_sf(j,k,l) = MAX(0d0,MIN(2d0*slope,5d-1*(1d0+slope),2d0))
                        ELSEIF (flux_lim == 3) THEN ! OSPRE (OP)
                            q_sf(j,k,l) = (15d-1*(slope**2d0+slope))/(slope**2d0+slope+1d0)
                        ELSEIF (flux_lim == 4) THEN ! SUPERBEE (SB)
                            q_sf(j,k,l) = MAX(0d0,MIN(1d0,2d0*slope),MIN(slope,2d0))
                        ELSEIF (flux_lim == 5) THEN ! SWEBY (SW) (beta = 1.5)
                            q_sf(j,k,l) = MAX(0d0,MIN(15d-1*slope,1d0),MIN(slope,15d-1))
                        ELSEIF (flux_lim == 6) THEN ! VAN ALBADA (VA)
                            q_sf(j,k,l) = (slope**2d0+slope)/(slope**2d0+1d0)
                        ELSEIF (flux_lim == 7) THEN ! VAN LEER (VL)
                            q_sf(j,k,l) = (ABS(slope) + slope)/(1d0 + ABS(slope))
                        END IF
                    END DO
                END DO
            END DO
        END SUBROUTINE s_derive_flux_limiter ! ---------------------------------
    
    
    
    
        !> Computes the local curvatures
        !!  @param i Fluid indicator
        !!  @param q_prim_vf Primitive variables
        !!  @param q_sf Curvature
        SUBROUTINE s_derive_curvature(i, q_prim_vf, q_sf) ! --------------------
    
            INTEGER, INTENT(IN) :: i
    
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q_prim_vf
    
            REAL(KIND(0d0)), DIMENSION (-offset_x%beg : m+offset_x%end , &
                                        -offset_y%beg : n+offset_y%end , &
                                        -offset_z%beg : p+offset_z%end), &
            INTENT(INOUT) :: q_sf

            REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:,:) :: alpha_unadv_sf
    
            REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:) :: A
            REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:) :: sol, b, AA
            REAL(KIND(0d0)) :: norm
            REAL(KIND(0d0)) :: xloc, yloc, zloc
            INTEGER :: ndim
            INTEGER :: j,k,l, jj, kk, ll, i1, i2
            INTEGER :: stencil_j_min, stencil_j_max
            INTEGER :: stencil_k_min, stencil_k_max
            INTEGER :: stencil_l_min, stencil_l_max
            TYPE(bounds_info) :: iz1

            IF (p > 0) THEN ! 3D simulation
                ALLOCATE(A(10,10))
                ALLOCATE(sol(10)); ALLOCATE(b(10)); ALLOCATE(AA(10))
            ELSE ! 2D simulation
                ALLOCATE(A(6,6))
                ALLOCATE(sol(6)); ALLOCATE(b(6)); ALLOCATE(AA(6))
            END IF

            IF ((i == num_fluids) .AND. (adv_alphan .NEQV. .TRUE.)) THEN
                ALLOCATE(alpha_unadv_sf(-offset_x%beg : m+offset_x%end , &
                                        -offset_y%beg : n+offset_y%end , &
                                        -offset_z%beg : p+offset_z%end))
                CALL s_derive_unadvected_volume_fraction(q_prim_vf, alpha_unadv_sf)
            END IF

            IF (p > 0) THEN
                iz1%beg = -offset_z%beg; iz1%end = p + offset_z%end
            ELSE
                iz1%beg = -1; iz1%end = 1
            END IF

            ! Parabolic fitting
            ndim = SIZE(sol,1)

            DO l = iz1%beg+1, iz1%end-1
                DO k = -offset_y%beg+1, n + offset_y%end-1
                    DO j = -offset_x%beg+1, m + offset_x%end-1
                        A(:,:) = 0d0
                        b(:) = 0d0
                        sol(:) = 0d0
                        AA(:) = 0d0

                        stencil_j_min = j-1; stencil_j_max = j+1
                        stencil_k_min = k-1; stencil_k_max = k+1
                        IF (p > 0) THEN
                            stencil_l_min = l-1; stencil_l_max = l+1
                        ELSE 
                            stencil_l_min = 0; stencil_l_max = 0
                        END IF

                        DO ll = stencil_l_min, stencil_l_max
                            DO kk = stencil_k_min, stencil_k_max
                                DO jj = stencil_j_min, stencil_j_max

                                    ! Ignore corner points in 3D stencil
                                    IF ((p > 0).AND.(ABS(jj-j)==1).AND.(ABS(kk-k)==1).AND.(ABS(ll-l)==1)) CYCLE

                                    ! Find distance between cell centers
                                    xloc = x_cc(jj) - x_cc(j)
                                    yloc = y_cc(kk) - y_cc(k)
                                    IF (p > 0) THEN
                                        zloc = z_cc(ll) - z_cc(l)
                                    END IF

                                    ! Compute operator
                                    AA(1) = 1d0
                                    AA(2) = xloc
                                    AA(3) = yloc
                                    AA(4) = 5d-1*xloc**2d0
                                    AA(5) = 5d-1*yloc**2d0
                                    AA(6) = xloc*yloc
                                    IF (p > 0) THEN
                                        AA(7) = zloc
                                        AA(8) = 5d-1*zloc**2d0
                                        AA(9) = yloc*zloc
                                        AA(10) = zloc*xloc
                                    END IF

                                    DO i1 = 1, ndim
                                        DO i2 = 1, ndim
                                            A(i1,i2) = A(i1,i2) + AA(i1)*AA(i2)
                                        END DO
                                    END DO

                                    ! Form RHS vector
                                    DO i1 = 1, ndim
                                        IF ((i == num_fluids) .AND. (adv_alphan .NEQV. .TRUE.)) THEN
                                            b(i1) = b(i1) + alpha_unadv_sf(jj,kk,ll)*AA(i1)
                                        ELSE
                                            b(i1) = b(i1) + q_prim_vf(E_idx+i)%sf(jj,kk,ll)*AA(i1)
                                        END IF
                                    END DO
                                END DO
                            END DO
                        END DO

                        CALL s_solve_linear_system(A,b,sol,ndim)

                        IF (p > 0) THEN
                            norm = SQRT(sol(2)**2d0 + sol(3)**2d0 + sol(7)**2d0)
                        ELSE
                            norm = SQRT(sol(2)**2d0 + sol(3)**2d0)
                        END IF
                        
                        IF (p > 0) THEN
                            q_sf(j,k,l) = -(sol(2)**2d0*sol(5) - 2d0*sol(2)*sol(3)*sol(6) + sol(2)**2d0*sol(8) + &
                                            sol(3)**2d0*sol(8) - 2d0*sol(3)*sol(7)*sol(9) + sol(3)**2d0*sol(4) + &
                                            sol(7)**2d0*sol(4) - 2d0*sol(7)*sol(2)*sol(10)+ sol(7)**2d0*sol(5)) / &
                                            MAX(norm,sgm_eps)**3d0
                        ELSE
                            q_sf(j,k,l) = -(sol(2)**2d0*sol(5) - 2d0*sol(2)*sol(3)*sol(6) + sol(3)**2d0*sol(4)) / &
                                MAX(norm,sgm_eps)**3d0
                        END IF

                    END DO
                END DO
            END DO
            
            DEALLOCATE(A,sol,b,AA)
            IF ((i == num_fluids) .AND. (adv_alphan .NEQV. .TRUE.)) DEALLOCATE(alpha_unadv_sf)
    
        END SUBROUTINE s_derive_curvature ! ------------------------------------
    
    
    
    
        !>  Computes the solution to the linear system Ax=b w/ sol = x
        !!  @param A Input matrix
        !!  @param b right-hand-side
        !!  @param sol Solution
        !!  @param ndim Problem size
        SUBROUTINE s_solve_linear_system(A,b,sol,ndim)

            INTEGER, INTENT(IN) :: ndim
            REAL(KIND(0d0)), DIMENSION(ndim,ndim), INTENT(INOUT) :: A
            REAL(KIND(0d0)), DIMENSION(ndim), INTENT(INOUT) :: b
            REAL(KIND(0d0)), DIMENSION(ndim), INTENT(OUT) :: sol
            INTEGER, DIMENSION(ndim) :: ipiv

            INTEGER :: nrhs, lda, ldb, info
            EXTERNAL DGESV

            INTEGER :: i,j,k


            ! Solve linear system using Intel MKL (Hooke)
!               nrhs = 1
!               lda = ndim
!               ldb = ndim
!    
!               CALL DGESV(ndim, nrhs, A, lda, ipiv, b, ldb, info)
!    
!               DO i = 1, ndim
!                   sol(i) = b(i)
!               END DO
!    
!               IF (info /= 0) THEN
!                   PRINT '(A)', 'Trouble solving linear system'
!                   CALL s_mpi_abort()
!               END IF

            ! Solve linear system using own linear solver (Thomson/Darter/Comet/Stampede)
                ! Forward elimination
                DO i = 1, ndim
                    ! Pivoting
                    j = i-1+MAXLOC(ABS(A(i:ndim,i)),1)
                    sol = A(i,:)
                    A(i,:) = A(j,:)
                    A(j,:) = sol
                    sol(1) = b(i)
                    b(i) = b(j)
                    b(j) = sol(1)
                    ! Elimination
                    b(i) = b(i) / A(i,i)
                    A(i,:) = A(i,:) / A(i,i)
                    DO k = i+1, ndim
                        b(k) = b(k) - A(k,i)*b(i)
                        A(k,:) = A(k,:) - A(k,i)*A(i,:)
                    END DO
                END DO

                ! Backward substitution
                DO i = ndim, 1, -1
                    sol(i) = b(i)
                    DO k = i+1, ndim
                        sol(i) = sol(i) - A(i,k)*sol(k)
                    END DO
                END DO

        END SUBROUTINE s_solve_linear_system ! -------------------------------------


        !>  This subroutine receives as inputs the indicator of the
        !!      component of the vorticity that should be outputted and
        !!      the primitive variables. From those inputs, it proceeds
        !!      to calculate values of the desired vorticity component,
        !!      which are subsequently stored in derived flow quantity
        !!      storage variable, q_sf.
        !!  @param i Vorticity component indicator
        !!  @param q_prim_vf Primitive variables
        !!  @param q_sf Vorticity component
        SUBROUTINE s_derive_vorticity_component(i, q_prim_vf, q_sf) ! ----------

            INTEGER, INTENT(IN) :: i
            
            TYPE(scalar_field), &
            DIMENSION(sys_size), &
            INTENT(IN) :: q_prim_vf
            
            REAL(KIND(0d0)), &
            DIMENSION( -offset_x%beg : m+offset_x%end  , &
                       -offset_y%beg : n+offset_y%end  , &
                       -offset_z%beg : p+offset_z%end ), &
            INTENT(INOUT) :: q_sf
            

            INTEGER :: j,k,l,r !< Generic loop iterators
            
            
            ! Computing the vorticity component in the x-coordinate direction
            IF(i == 1) THEN
                DO l = -offset_z%beg, p+offset_z%end
                    DO k = -offset_y%beg, n+offset_y%end
                        DO j = -offset_x%beg, m+offset_x%end
                            
                            q_sf(j,k,l) = 0d0
                            
                            DO r = -fd_number, fd_number
                                IF (grid_geometry == 3) THEN
                                    q_sf(j,k,l) =                            &
                                        q_sf(j,k,l) + 1d0/y_cc(k)          * &
                                        (    fd_coeff_y(r,k) * y_cc(r+k)   * &
                                        q_prim_vf( mom_idx%end )%sf(j,r+k,l) &
                                           - fd_coeff_z(r,l) *               &
                                        q_prim_vf(mom_idx%beg+1)%sf(j,k,r+l))
                                ELSE
                                    q_sf(j,k,l) =                                  &
                                          q_sf(j,k,l) + fd_coeff_y(r,k) *          &
                                          q_prim_vf( mom_idx%end )%sf( j ,r+k, l ) &
                                                      - fd_coeff_z(r,l) *          &
                                          q_prim_vf(mom_idx%beg+1)%sf( j , k ,r+l)
                                END IF
                            END DO
                            
                        END DO
                    END DO
                END DO
                
                
            ! Computing the vorticity component in the y-coordinate direction
            ELSEIF(i == 2) THEN
                DO l = -offset_z%beg, p+offset_z%end
                    DO k = -offset_y%beg, n+offset_y%end
                        DO j = -offset_x%beg, m+offset_x%end
                            
                            q_sf(j,k,l) = 0d0
                            
                            DO r = -fd_number, fd_number
                                IF (grid_geometry == 3) THEN
                                    q_sf(j,k,l) =                                  &
                                          q_sf(j,k,l) + fd_coeff_z(r,l)/y_cc(k) *  &
                                          q_prim_vf( mom_idx%beg )%sf( j , k ,r+l) &
                                                      - fd_coeff_x(r,j) *          &
                                          q_prim_vf( mom_idx%end )%sf(r+j, k , l )
                                ELSE
                                    q_sf(j,k,l) =                                  &
                                          q_sf(j,k,l) + fd_coeff_z(r,l) *          &
                                          q_prim_vf( mom_idx%beg )%sf( j , k ,r+l) &
                                                      - fd_coeff_x(r,j) *          &
                                          q_prim_vf( mom_idx%end )%sf(r+j, k , l )
                                END IF
                            END DO
                            
                        END DO
                    END DO
                END DO
                
                
            ! Computing the vorticity component in the z-coordinate direction
            ELSE
                DO l = -offset_z%beg, p+offset_z%end
                    DO k = -offset_y%beg, n+offset_y%end
                        DO j = -offset_x%beg, m+offset_x%end
                            
                            q_sf(j,k,l) = 0d0
                            
                            DO r = -fd_number, fd_number
                                q_sf(j,k,l) = &
                                      q_sf(j,k,l) + fd_coeff_x(r,j) *          &
                                      q_prim_vf(mom_idx%beg+1)%sf(r+j, k , l ) &
                                                  - fd_coeff_y(r,k) *          &
                                      q_prim_vf( mom_idx%beg )%sf( j ,r+k, l )
                            END DO
                            
                        END DO
                    END DO
                END DO
            END IF
            
            
        END SUBROUTINE s_derive_vorticity_component ! --------------------------
        
        
        
        
        !>  This subroutine gets as inputs the conservative variables
        !!      and density. From those inputs, it proceeds to calculate
        !!      the values of the numerical Schlieren function, which are
        !!      subsequently stored in the derived flow quantity storage
        !!      variable, q_sf. 
        !!  @param q_cons_vf Conservative variables
        !!  @param rho_sf Density
        !!  @param q_sf Numerical Schlieren function
        SUBROUTINE s_derive_numerical_schlieren_function(q_cons_vf,rho_sf, q_sf)

            TYPE(scalar_field), &
            DIMENSION(sys_size), &
            INTENT(IN) :: q_cons_vf
            
            REAL(KIND(0d0)), &
            DIMENSION( -buff_size     :  m+buff_size       , &
                       -buff_size     :  n+buff_size       , &
                       -buff_size*flg : (p+buff_size)*flg ), &
            INTENT(IN) :: rho_sf
            
            REAL(KIND(0d0)), &
            DIMENSION( -offset_x%beg : m+offset_x%end  , &
                       -offset_y%beg : n+offset_y%end  , &
                       -offset_z%beg : p+offset_z%end ), &
            INTENT(INOUT) :: q_sf
            

            REAL(KIND(0d0)) :: drho_dx, drho_dy, drho_dz !<
            !! Spatial derivatives of the density in the x-, y- and z-directions
            
            REAL(KIND(0d0)), DIMENSION(2) :: gm_rho_max !<
            !! Maximum value of the gradient magnitude (gm) of the density field
            !! in entire computational domain and not just the local sub-domain.
            !! The first position in the variable contains the maximum value and
            !! the second contains the rank of the processor on which it occured.
            
            REAL(KIND(0d0)) :: alpha_unadv !< Unadvected volume fraction
            

            INTEGER :: i,j,k,l !< Generic loop iterators
            
            
            ! Computing Gradient Magnitude of Density ==========================
            
            ! Contributions from the x- and y-coordinate directions
            DO l = -offset_z%beg, p+offset_z%end
                DO k = -offset_y%beg, n+offset_y%end
                    DO j = -offset_x%beg, m+offset_x%end
                        
                        drho_dx = 0d0
                        drho_dy = 0d0
                        
                        DO i = -fd_number, fd_number
                            drho_dx = drho_dx + fd_coeff_x(i,j)*rho_sf(i+j,k,l)
                            drho_dy = drho_dy + fd_coeff_y(i,k)*rho_sf(j,i+k,l)
                        END DO
                        
                        gm_rho_sf(j,k,l) = drho_dx*drho_dx + drho_dy*drho_dy
                        
                    END DO
                END DO
            END DO
            
            ! Contribution from the z-coordinate direction
            IF(p > 0) THEN
                DO l = -offset_z%beg, p+offset_z%end
                    DO k = -offset_y%beg, n+offset_y%end
                        DO j = -offset_x%beg, m+offset_x%end
                            
                            drho_dz = 0d0
                            
                            DO i = -fd_number, fd_number
                                IF (grid_geometry == 3) THEN
                                    drho_dz = drho_dz + fd_coeff_z(i,l)/y_cc(k) * &
                                                        rho_sf(j,k,i+l)
                                ELSE
                                    drho_dz = drho_dz + fd_coeff_z(i,l) * &
                                                        rho_sf(j,k,i+l)
                                END IF
                            END DO
                            
                            gm_rho_sf(j,k,l) = gm_rho_sf(j,k,l) &
                                             + drho_dz*drho_dz
                            
                        END DO
                    END DO
                END DO
            END IF
            
            ! Up until now, only the dot product of the gradient of the density
            ! field has been calculated and stored in the gradient magnitude of
            ! density variable. So now we proceed to take the square-root as to
            ! complete the desired calculation.
            gm_rho_sf = SQRT(gm_rho_sf)
            
            ! ==================================================================
            
            
            ! Determining the local maximum of the gradient magnitude of density
            ! and bookkeeping the result, along with rank of the local processor
            gm_rho_max = (/ MAXVAL(gm_rho_sf), REAL(proc_rank, KIND(0d0)) /)
            
            
            ! Comparing the local maximum gradient magnitude of the density on
            ! this processor to the those computed on the remaining processors.
            ! This allows for the global maximum to be computed and the rank of
            ! the processor on which it has occured to be recorded.
            IF(num_procs > 1) CALL s_mpi_reduce_maxloc(gm_rho_max)
            
            
            ! Computing Numerical Schlieren Function ===========================
            
            ! The form of the numerical Schlieren function depends on the choice
            ! of the multicomponent flow model. For the gamma/pi_inf model, the
            ! exponential of the negative, normalized, gradient magnitude of the
            ! density is computed. For the volume fraction model, the amplitude
            ! of the exponential's inside is also modulated with respect to the
            ! identity of the fluid in which the function is evaluated. For more
            ! information, refer to Marquina and Mulet (2003).
            
            IF(model_eqns == 1) THEN                    ! Gamma/pi_inf model
                q_sf = -gm_rho_sf/gm_rho_max(1)
                
            ELSE                                        ! Volume fraction model
                DO l = -offset_z%beg, p+offset_z%end
                    DO k = -offset_y%beg, n+offset_y%end
                        DO j = -offset_x%beg, m+offset_x%end
                            
                            q_sf(j,k,l) = 0d0
                            
                            DO i = 1, adv_idx%end-E_idx
                                q_sf(j,k,l) =                              &
                                        q_sf(j,k,l) - schlieren_alpha(i) * &
                                        q_cons_vf(i+E_idx)%sf(j,k,l)     * &
                                        gm_rho_sf(j,k,l) / gm_rho_max(1)
                            END DO
                            
                            IF(adv_alphan .NEQV. .TRUE.) THEN
                                
                                alpha_unadv = 1d0
                                
                                DO i = 1, num_fluids-1
                                    alpha_unadv = alpha_unadv                  &
                                                - q_cons_vf(i+E_idx)%sf(j,k,l)
                                END DO
                                
                                q_sf(j,k,l) = q_sf(j,k,l)                    &
                                            - schlieren_alpha(num_fluids)  * &
                                              alpha_unadv*gm_rho_sf(j,k,l) / &
                                              gm_rho_max(1)
                                
                            END IF
                            
                        END DO
                    END DO
                END DO
            END IF
            
            ! Up until now, only the inside of the exponential of the numerical
            ! Schlieren function has been evaluated and stored. Then, to finish
            ! the computation, the exponential of the inside quantity is taken.
            q_sf = EXP(q_sf)
            
            ! ==================================================================
            
            
        END SUBROUTINE s_derive_numerical_schlieren_function ! -----------------
        
        
        
        !>  Deallocation procedures for the module        
        SUBROUTINE s_finalize_derived_variables_module() ! -------------------

            
            
            ! Deallocating the variable containing the gradient magnitude of the
            ! density field provided that the numerical Schlieren function was
            ! was outputted during the post-process
            IF(schlieren_wrt) DEALLOCATE(gm_rho_sf)
            
            
            ! Deallocating the variables that might have been used to bookkeep
            ! the finite-difference coefficients in the x-, y- and z-directions
            IF(ALLOCATED(fd_coeff_x)) DEALLOCATE(fd_coeff_x)
            IF(ALLOCATED(fd_coeff_y)) DEALLOCATE(fd_coeff_y)
            IF(ALLOCATED(fd_coeff_z)) DEALLOCATE(fd_coeff_z)
            
            
        END SUBROUTINE s_finalize_derived_variables_module ! -----------------
        
        
        
        
        
END MODULE m_derived_variables
