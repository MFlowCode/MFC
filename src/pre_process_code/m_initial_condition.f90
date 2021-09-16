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
!! @file m_initial_condition.f90
!! @brief Contains module m_initial_condition
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief This module provides a platform that is analagous to constructive
!!              solid geometry techniques and in this way allows for the creation
!!              of a wide variety of initial conditions. Several 1D, 2D and 3D
!!              fundamental geometries are included that may further be combined
!!              into more complex shapes. This is achieved by carefully setting
!!              up the order in which the patches are laid out in the domain and
!!              specifying the priority that each patch has over the preceeding
!!              ones. The resulting shapes may be identified both by the values
!!              of their primitive variables and the associated patch identities.
!!              Note that the user may choose to read in and modify a preexisting
!!              initial condition. The module m_start_up.f90 is responsible for
!!             reading in the relevant data files.
MODULE m_initial_condition
    
    
    ! Dependencies =============================================================
    USE m_derived_types         ! Definitions of the derived types
    
    USE m_global_parameters     ! Global parameters for the code
    
    USE m_variables_conversion  ! Subroutines to change the state variables from
                                ! one form to another
    ! ==========================================================================
    
    
    IMPLICIT NONE
    
    
    !> Abstract interface to the two subroutines that assign the patch primitive
    !! variables, either mixture or species, depending on the subroutine, to a
    !! particular cell in the computational domain
    ABSTRACT INTERFACE
 
        !> Skeleton of s_assign_patch_mixture_primitive_variables
        !!      and s_assign_patch_species_primitive_variables
        !! @param patch_id is the patch identifier
        !! @param j (x) cell index in which the mixture or species primitive variables from the indicated patch areassigned
        !! @param k (y,th) cell index in which the mixture or species primitive variables from the indicated patch areassigned
        !! @param l (z) cell index in which the mixture or species primitive variables from the indicated patch areassigned
        SUBROUTINE s_assign_patch_xxxxx_primitive_variables( patch_id,j,k,l)
           
            INTEGER, INTENT(IN) :: patch_id
            INTEGER, INTENT(IN) :: j,k,l
            
        END SUBROUTINE s_assign_patch_xxxxx_primitive_variables
        
    END INTERFACE
    
    ! NOTE: The abstract interface allows for the declaration of a pointer to
    ! a procedure such that the choice of the model equations does not have to
    ! be queried every time the patch primitive variables are to be assigned in
    ! a cell in the computational domain.
    TYPE(scalar_field), ALLOCATABLE, DIMENSION(:) :: q_prim_vf !< primitive variables
    
    TYPE(scalar_field), ALLOCATABLE, DIMENSION(:) :: q_cons_vf !< conservative variables
    
    TYPE(scalar_field) :: alf_sum
    
    REAL(KIND(0d0)) :: x_centroid, y_centroid, z_centroid 
    REAL(KIND(0d0)) :: length_x, length_y, length_z
    REAL(KIND(0d0)) :: radius
    REAL(KIND(0d0)) :: epsilon, beta
    INTEGER         :: smooth_patch_id
    REAL(KIND(0d0)) :: smooth_coeff !< 
    !! These variables are analogous in both meaning and use to the similarly
    !! named components in the ic_patch_parameters type (see m_derived_types.f90
    !! for additional details). They are employed as a means to more concisely
    !! perform the actions necessary to lay out a particular patch on the grid.

    TYPE(bounds_info) :: x_boundary, y_boundary, z_boundary  !<
    !! These variables combine the centroid and length parameters associated with
    !! a particular patch to yield the locations of the patch boundaries in the
    !! x-, y- and z-coordinate directions. They are used as a means to concisely
    !! perform the actions necessary to lay out a particular patch on the grid.
 
    
    REAL(KIND(0d0)) :: a, b, c, d !<
    !! When a line or a plane sweep patch geometry is employed, these variables
    !! represent the coefficients associated with the equation describing the
    !! said line or plane.
    

    REAL(KIND(0d0)) :: eta !<
    !! In the case that smoothing of patch boundaries is enabled and the boundary
    !! between two adjacent patches is to be smeared out, this variable's purpose
    !! is to act as a pseudo volume fraction to indicate the contribution of each
    !! patch toward the composition of a cell's fluid state.


    INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: patch_id_fp !<
    !! Bookkepping variable used to track the patch identities (id) associated
    !! with each of the cells in the computational domain. Note that only one
    !! patch identity may be associated with any one cell.
 
    REAL(KIND(0d0)) :: cart_y, cart_z
    REAL(KIND(0d0)) :: sph_phi !<
    !! Variables to be used to hold cell locations in Cartesian coordinates if
    !! 3D simulation is using cylindrical coordinates

    
    PROCEDURE(s_assign_patch_xxxxx_primitive_variables), &
    POINTER :: s_assign_patch_primitive_variables => NULL() !<
    !! Depending on the multicomponent flow model, this variable is a pointer to
    !! either the subroutine s_assign_patch_mixture_primitive_variables, or the
    !! subroutine s_assign_patch_species_primitive_variables
 
    
    
    CONTAINS

        !>  This subroutine assigns the mixture primitive variables
        !!              of the patch designated by the patch_id, to the cell that
        !!              is designated by the indexes (j,k,l). In addition, the
        !!              variable bookkeeping the patch identities in the entire
        !!              domain is updated with the new assignment. Note that if
        !!              the smoothing of the patch's boundaries is employed, the
        !!              ensuing primitive variables in the cell will be a type of
        !!              combination of the current patch's primitive variables
        !!              with those of the smoothing patch. The specific details
        !!              of the combination may be found in Shyue's work (1998).
        !! @param patch_id the patch identifier
        !! @param j  the x-dir node index
        !! @param k  the y-dir node index
        !! @param l  the z-dir node index
        SUBROUTINE s_assign_patch_mixture_primitive_variables( patch_id,j,k,l)
            
            
            INTEGER, INTENT(IN) :: patch_id
            INTEGER, INTENT(IN) :: j,k,l
            
            REAL(KIND(0d0))                               :: rho    !< density
            REAL(KIND(0d0)), DIMENSION(E_idx-mom_idx%beg) :: vel    !< velocity
            REAL(KIND(0d0))                               :: pres   !< pressure
            REAL(KIND(0d0))                               :: gamma  !< specific heat ratio function
            
            INTEGER :: i !< generic loop operator
            
            
            ! Assigning the mixture primitive variables of a uniform state patch
            IF(patch_icpp(patch_id)%geometry /= 6) THEN
                
                ! Transferring the identity of the smoothing patch
                smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
                
                ! Density
                q_prim_vf(1)%sf(j,k,l) = &
                              eta    *patch_icpp(    patch_id   )%rho &
                        + (1d0 - eta)*patch_icpp(smooth_patch_id)%rho
                
                ! Velocity
                DO i = 1, E_idx-mom_idx%beg
                    q_prim_vf(i+1)%sf(j,k,l) = &
                        1d0 / q_prim_vf(1)%sf(j,k,l) * &
                        (     eta    *patch_icpp(    patch_id   )%rho    &
                                     *patch_icpp(    patch_id   )%vel(i) &
                        + (1d0 - eta)*patch_icpp(smooth_patch_id)%rho    &
                                     *patch_icpp(smooth_patch_id)%vel(i) )
                END DO
                
                ! Specific heat ratio function
                q_prim_vf(gamma_idx)%sf(j,k,l) = &
                              eta    *patch_icpp(    patch_id   )%gamma &
                        + (1d0 - eta)*patch_icpp(smooth_patch_id)%gamma
                
                ! Pressure
                q_prim_vf(E_idx)%sf(j,k,l) = &
                        1d0 / q_prim_vf(gamma_idx)%sf(j,k,l) * &
                        (     eta    *patch_icpp(    patch_id   )%gamma &
                                     *patch_icpp(    patch_id   )%pres  &
                        + (1d0 - eta)*patch_icpp(smooth_patch_id)%gamma &
                                     *patch_icpp(smooth_patch_id)%pres  )
                
                ! Liquid stiffness function
                q_prim_vf(pi_inf_idx)%sf(j,k,l) = &
                              eta    *patch_icpp(    patch_id   )%pi_inf &
                        + (1d0 - eta)*patch_icpp(smooth_patch_id)%pi_inf
                
            ! Assigning mixture primitive variables of isentropic vortex patch
            ELSE
                
                ! Transferring the centroid of the isentropic vortex patch, the
                ! amplitude of its disturbance and also its domain of influence
                x_centroid = patch_icpp(patch_id)%x_centroid
                y_centroid = patch_icpp(patch_id)%y_centroid
                epsilon    = patch_icpp(patch_id)%epsilon
                beta       = patch_icpp(patch_id)%beta
                
                ! Reference density, velocity, pressure and specific heat ratio
                ! function of the isentropic vortex patch
                rho    = patch_icpp(patch_id)%rho
                vel    = patch_icpp(patch_id)%vel
                pres   = patch_icpp(patch_id)%pres
                gamma  = patch_icpp(patch_id)%gamma
                
                ! Density
                q_prim_vf(1)%sf(j,k,0) = &
                             rho*( 1d0 - (rho/pres)*(epsilon/(2d0*pi)) * &
                                   (epsilon/(8d0*beta*(gamma+1d0)*pi)) * &
                      EXP( 2d0*beta*( 1d0 - (x_cc(j)-x_centroid)**2d0    &
                                          - (y_cc(k)-y_centroid)**2d0) ) &
                                 )**gamma
                
                ! Velocity
                q_prim_vf(2)%sf(j,k,0) = &
                      vel(1) - (y_cc(k)-y_centroid)*(epsilon/(2d0*pi)) * &
                          EXP( beta*( 1d0 - (x_cc(j)-x_centroid)**2d0    &
                                          - (y_cc(k)-y_centroid)**2d0) )
                q_prim_vf(3)%sf(j,k,0) = &
                      vel(2) + (x_cc(j)-x_centroid)*(epsilon/(2d0*pi)) * &
                          EXP( beta*( 1d0 - (x_cc(j)-x_centroid)**2d0    &
                                          - (y_cc(k)-y_centroid)**2d0) )
                
                ! Pressure
                q_prim_vf(4)%sf(j,k,0) = &
                            pres*( 1d0 - (rho/pres)*(epsilon/(2d0*pi)) * &
                                   (epsilon/(8d0*beta*(gamma+1d0)*pi)) * &
                      EXP( 2d0*beta*( 1d0 - (x_cc(j)-x_centroid)**2d0    &
                                          - (y_cc(k)-y_centroid)**2d0) ) &
                                 )**(gamma+1d0)
                
                ! Specific heat ratio function
                q_prim_vf(5)%sf(j,k,0) = gamma
                
                ! Liquid stiffness function
                q_prim_vf(6)%sf(j,k,0) = 0d0
                
            END IF
            
            
            ! Updating the patch identities bookkeeping variable
            IF(1d0 - eta < 1d-16) patch_id_fp(j,k,l) = patch_id
            
            
        END SUBROUTINE s_assign_patch_mixture_primitive_variables ! ------------
        
        !>  This subroutine assigns the species primitive variables. This follows
        !!  s_assign_patch_species_primitive_variables with adaptation for
        !!  ensemble-averaged bubble modeling
        !! @param patch_id the patch identifier
        !! @param j  the x-dir node index
        !! @param k  the y-dir node index
        !! @param l  the z-dir node index
        SUBROUTINE s_assign_patch_species_primitive_variables_bubbles( patch_id,j,k,l)
            
            INTEGER, INTENT(IN) :: patch_id
            INTEGER, INTENT(IN) :: j,k,l
            
            ! Density, the specific heat ratio function and the liquid stiffness
            ! function, respectively, obtained from the combination of primitive
            ! variables of the current and smoothing patches
            REAL(KIND(0d0)) :: rho          !< density
            REAL(KIND(0d0)) :: gamma        
            REAL(KIND(0d0)) :: lit_gamma    !< specific heat ratio
            REAL(KIND(0d0)) :: pi_inf       !< stiffness from SEOS
            REAL(KIND(0d0)) :: orig_rho     
            REAL(KIND(0d0)) :: orig_gamma
            REAL(KIND(0d0)) :: orig_pi_inf 
            REAL(KIND(0d0)) :: muR, muV

            REAL(KIND(0d0)), DIMENSION(sys_size) :: orig_prim_vf !<
            !! Vector to hold original values of cell for smoothing purposes
            

            INTEGER :: i  !< Generic loop iterator
            
            ! Transferring the identity of the smoothing patch
            smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id

            ! Transferring original primitive variables
            DO i = 1, sys_size
                orig_prim_vf(i) = q_prim_vf(i)%sf(j,k,l)
            END DO

            IF (mpp_lim .AND. bubbles) THEN
                !adjust volume fractions, according to modeled gas void fraction
                alf_sum%sf = 0d0
                DO i = adv_idx%beg, adv_idx%end - 1
                    alf_sum%sf = alf_sum%sf + q_prim_vf(i)%sf
                END DO

                DO i = adv_idx%beg, adv_idx%end-1
                    q_prim_vf(i)%sf = q_prim_vf(i)%sf * (1.d0 - q_prim_vf(alf_idx)%sf) &
                            / alf_sum%sf
                END DO
            END IF
                
            ! Computing Mixture Variables from Original Primitive Variables
            CALL s_convert_species_to_mixture_variables(    &
                                        q_prim_vf, j,k,l,   &
                                        orig_rho,           &
                                        orig_gamma,         &
                                        orig_pi_inf )

            ! Computing Mixture Variables of Current Patch =====================

            ! Volume fraction(s)
            DO i = adv_idx%beg, adv_idx%end
                q_prim_vf(i)%sf(j,k,l) = patch_icpp(patch_id)%alpha(i-E_idx)
            END DO

            IF (mpp_lim .AND. bubbles) THEN
                !adjust volume fractions, according to modeled gas void fraction
                alf_sum%sf = 0d0
                DO i = adv_idx%beg, adv_idx%end - 1
                    alf_sum%sf = alf_sum%sf + q_prim_vf(i)%sf
                END DO

                DO i = adv_idx%beg, adv_idx%end-1
                    q_prim_vf(i)%sf = q_prim_vf(i)%sf * (1.d0 - q_prim_vf(alf_idx)%sf) &
                            / alf_sum%sf
                END DO
            END IF
            
            ! Partial densities
            IF (model_eqns .NE. 4) THEN
                DO i = 1, cont_idx%end
                    q_prim_vf(i)%sf(j,k,l) = patch_icpp(patch_id)%alpha_rho(i)
                END DO
            END IF

            ! Bubbles variables
            IF (bubbles) THEN
                DO i = 1,nb
                    muR = R0(i)*patch_icpp(patch_id)%r0 ! = R0(i)
                    muV = patch_icpp(patch_id)%v0 ! = 0

                    IF (qbmm) THEN
                        IF (dist_type == 1) THEN
                            q_prim_vf(bub_idx%fullmom(i,0,0))%sf(j,k,l) = 1d0
                            q_prim_vf(bub_idx%fullmom(i,1,0))%sf(j,k,l) = muR
                            q_prim_vf(bub_idx%fullmom(i,0,1))%sf(j,k,l) = muV
                            q_prim_vf(bub_idx%fullmom(i,2,0))%sf(j,k,l) = muR**2d0 + sigR**2d0
                            q_prim_vf(bub_idx%fullmom(i,1,1))%sf(j,k,l) = muR*muV + rhoRV*sigR*sigV
                            q_prim_vf(bub_idx%fullmom(i,0,2))%sf(j,k,l) = muV**2d0 + sigV**2d0
                        ELSE IF (dist_type == 2) THEN
                            q_prim_vf(bub_idx%fullmom(i,0,0))%sf(j,k,l) = 1d0
                            q_prim_vf(bub_idx%fullmom(i,1,0))%sf(j,k,l) = dexp((sigR**2d0)/2d0)*muR
                            q_prim_vf(bub_idx%fullmom(i,0,1))%sf(j,k,l) = muV
                            q_prim_vf(bub_idx%fullmom(i,2,0))%sf(j,k,l) = dexp((sigR**2d0)*2d0)*(muR**2d0)
                            q_prim_vf(bub_idx%fullmom(i,1,1))%sf(j,k,l) = dexp((sigR**2d0)/2d0)*muR*muV
                            q_prim_vf(bub_idx%fullmom(i,0,2))%sf(j,k,l) = muV**2d0 + sigV**2d0
                        END IF

                        if (j==0 .and. k==0 .and. l==0) then
                            print*, 'moments @ (0,0,0): ', &
                                q_prim_vf(bub_idx%fullmom(i,0,0))%sf(j,k,l), &
                                q_prim_vf(bub_idx%fullmom(i,1,0))%sf(j,k,l), &
                                q_prim_vf(bub_idx%fullmom(i,0,1))%sf(j,k,l), &
                                q_prim_vf(bub_idx%fullmom(i,2,0))%sf(j,k,l), &
                                q_prim_vf(bub_idx%fullmom(i,1,1))%sf(j,k,l), &
                                q_prim_vf(bub_idx%fullmom(i,0,2))%sf(j,k,l)
                        end if
                    ELSE
                        q_prim_vf(bub_idx%rs(i))%sf(j,k,l) = muR
                        q_prim_vf(bub_idx%vs(i))%sf(j,k,l) = muV
                        IF ( .NOT. polytropic ) THEN
                            q_prim_vf(bub_idx%ps(i))%sf(j,k,l) = patch_icpp(patch_id)%p0
                            q_prim_vf(bub_idx%ms(i))%sf(j,k,l) = patch_icpp(patch_id)%m0
                        END IF
                    END IF
                END DO
            END IF


            ! Density and the specific heat ratio and liquid stiffness functions
            CALL s_convert_species_to_mixture_variables(    &
                                q_prim_vf, j,k,l,           &
                                patch_icpp(patch_id)%rho,   &
                                patch_icpp(patch_id)%gamma, &
                                patch_icpp(patch_id)%pi_inf )
            
            ! ==================================================================
            
            ! Computing Mixture Variables of Smoothing Patch ===================
            
            IF (model_eqns .NE. 4) THEN
                ! Partial densities
                DO i = 1, cont_idx%end
                    q_prim_vf(i)%sf(j,k,l) = patch_icpp(smooth_patch_id)%alpha_rho(i)
                END DO
            END IF
           
            ! Volume fraction(s)
            DO i = adv_idx%beg, adv_idx%end
                q_prim_vf(i)%sf(j,k,l) = patch_icpp(smooth_patch_id)%alpha(i-E_idx)
            END DO

            IF (mpp_lim .AND. bubbles) THEN
                !adjust volume fractions, according to modeled gas void fraction
                alf_sum%sf = 0d0
                DO i = adv_idx%beg, adv_idx%end - 1
                    alf_sum%sf = alf_sum%sf + q_prim_vf(i)%sf
                END DO

                DO i = adv_idx%beg, adv_idx%end-1
                    q_prim_vf(i)%sf = q_prim_vf(i)%sf * (1.d0 - q_prim_vf(alf_idx)%sf) &
                            / alf_sum%sf
                END DO
            END IF

            ! Bubbles variables
            IF (bubbles) THEN
                DO i = 1,nb
                    muR = R0(i)*patch_icpp(smooth_patch_id)%r0 ! = R0(i)
                    muV = V0(i)*patch_icpp(smooth_patch_id)%v0 ! = 0
                    IF (qbmm) THEN
                        ! Initialize the moment set
                        IF (dist_type == 1) THEN
                            q_prim_vf(bub_idx%fullmom(i,0,0))%sf(j,k,l) = 1d0
                            q_prim_vf(bub_idx%fullmom(i,1,0))%sf(j,k,l) = muR
                            q_prim_vf(bub_idx%fullmom(i,0,1))%sf(j,k,l) = muV
                            q_prim_vf(bub_idx%fullmom(i,2,0))%sf(j,k,l) = muR**2d0 + sigR**2d0
                            q_prim_vf(bub_idx%fullmom(i,1,1))%sf(j,k,l) = muR*muV + rhoRV*sigR*sigV
                            q_prim_vf(bub_idx%fullmom(i,0,2))%sf(j,k,l) = muV**2d0 + sigV**2d0
                        ELSE IF (dist_type == 2) THEN
                            q_prim_vf(bub_idx%fullmom(i,0,0))%sf(j,k,l) = 1d0
                            q_prim_vf(bub_idx%fullmom(i,1,0))%sf(j,k,l) = dexp((sigR**2d0)/2d0)*muR
                            q_prim_vf(bub_idx%fullmom(i,0,1))%sf(j,k,l) = muV
                            q_prim_vf(bub_idx%fullmom(i,2,0))%sf(j,k,l) = dexp((sigR**2d0)*2d0)*(muR**2d0)
                            q_prim_vf(bub_idx%fullmom(i,1,1))%sf(j,k,l) = dexp((sigR**2d0)/2d0)*muR*muV
                            q_prim_vf(bub_idx%fullmom(i,0,2))%sf(j,k,l) = muV**2d0 + sigV**2d0
                        END IF
                    ELSE
                        q_prim_vf(bub_idx%rs(i))%sf(j,k,l) = muR
                        q_prim_vf(bub_idx%vs(i))%sf(j,k,l) = muV
                        IF ( .NOT. polytropic ) THEN
                            q_prim_vf(bub_idx%ps(i))%sf(j,k,l) = patch_icpp(patch_id)%p0
                            q_prim_vf(bub_idx%ms(i))%sf(j,k,l) = patch_icpp(patch_id)%m0
                        END IF
                    END IF
                END DO
            END IF

           
            ! Density and the specific heat ratio and liquid stiffness functions
            CALL s_convert_species_to_mixture_variables(       &
                            q_prim_vf, j,k,l,                  &
                            patch_icpp(smooth_patch_id)%rho,   &
                            patch_icpp(smooth_patch_id)%gamma, &
                            patch_icpp(smooth_patch_id)%pi_inf )
            
            ! ==================================================================
        
            ! Pressure
            q_prim_vf(E_idx)%sf(j,k,l) = &
                    (eta * patch_icpp(patch_id)%pres  &
                     + (1d0 - eta) * orig_prim_vf(E_idx))    

            ! Volume fractions \alpha
            DO i = adv_idx%beg, adv_idx%end
                q_prim_vf(i)%sf(j,k,l) = &
                    eta *patch_icpp(patch_id)%alpha(i-E_idx) &
                    + (1d0 - eta)   *orig_prim_vf(i)
            END DO

            IF (mpp_lim .AND. bubbles) THEN
                !adjust volume fractions, according to modeled gas void fraction
                alf_sum%sf = 0d0
                DO i = adv_idx%beg, adv_idx%end - 1
                    alf_sum%sf = alf_sum%sf + q_prim_vf(i)%sf
                END DO

                DO i = adv_idx%beg, adv_idx%end-1
                    q_prim_vf(i)%sf = q_prim_vf(i)%sf * (1.d0 - q_prim_vf(alf_idx)%sf) &
                            / alf_sum%sf
                END DO
            END IF


            ! Partial densities \alpha \rho
            IF (model_eqns .NE. 4) THEN
                !mixture density is an input
                DO i = 1, cont_idx%end
                    q_prim_vf(i)%sf(j,k,l) = &
                        eta *patch_icpp(patch_id)%alpha_rho(i) &
                        + (1d0 - eta) * orig_prim_vf(i)
                END DO
            ELSE
                !get mixture density from pressure via Tait EOS
                pi_inf = fluid_pp(1)%pi_inf
                gamma = fluid_pp(1)%gamma
                lit_gamma = (1.d0+gamma)/gamma

                ! \rho = (( p_l + pi_inf)/( p_ref + pi_inf))**(1/little_gam) * rhoref(1-alf)
                q_prim_vf(1)%sf(j,k,l) = &
                    (((q_prim_vf(E_idx)%sf(j,k,l) + pi_inf)/(pref + pi_inf))**(1/lit_gamma)) * &
                    rhoref*(1-q_prim_vf(alf_idx)%sf(j,k,l))
            END IF

            ! Density and the specific heat ratio and liquid stiffness functions
            CALL s_convert_species_to_mixture_variables( q_prim_vf , j,k,l, &
                                                          rho,gamma,pi_inf  )
            
            ! Velocity
            DO i = 1, E_idx-mom_idx%beg
                q_prim_vf(i+cont_idx%end)%sf(j,k,l) = &
                        (eta * patch_icpp(patch_id)%vel(i) & 
                         + (1d0-eta)*orig_prim_vf(i+cont_idx%end))
            END DO

            ! Smoothed bubble variables
            IF (bubbles) THEN
                DO i = 1,nb
                    muR = R0(i)*patch_icpp(patch_id)%r0 ! = 1*R0(i)
                    muV = V0(i)*patch_icpp(patch_id)%v0 ! = 1*V0(i)
                    IF (qbmm) THEN
                        ! Initialize the moment set
                        IF (dist_type == 1) THEN
                            q_prim_vf(bub_idx%fullmom(i,0,0))%sf(j,k,l) = 1d0
                            q_prim_vf(bub_idx%fullmom(i,1,0))%sf(j,k,l) = muR
                            q_prim_vf(bub_idx%fullmom(i,0,1))%sf(j,k,l) = muV
                            q_prim_vf(bub_idx%fullmom(i,2,0))%sf(j,k,l) = muR**2d0 + sigR**2d0
                            q_prim_vf(bub_idx%fullmom(i,1,1))%sf(j,k,l) = muR*muV + rhoRV*sigR*sigV
                            q_prim_vf(bub_idx%fullmom(i,0,2))%sf(j,k,l) = muV**2d0 + sigV**2d0
                        ELSE IF (dist_type == 2) THEN
                            q_prim_vf(bub_idx%fullmom(i,0,0))%sf(j,k,l) = 1d0
                            q_prim_vf(bub_idx%fullmom(i,1,0))%sf(j,k,l) = dexp((sigR**2d0)/2d0)*muR
                            q_prim_vf(bub_idx%fullmom(i,0,1))%sf(j,k,l) = muV
                            q_prim_vf(bub_idx%fullmom(i,2,0))%sf(j,k,l) = dexp((sigR**2d0)*2d0)*(muR**2d0)
                            q_prim_vf(bub_idx%fullmom(i,1,1))%sf(j,k,l) = dexp((sigR**2d0)/2d0)*muR*muV
                            q_prim_vf(bub_idx%fullmom(i,0,2))%sf(j,k,l) = muV**2d0 + sigV**2d0
                        END IF
                    ELSE
                        ! q_prim_vf(bub_idx%rs(i))%sf(j,k,l) = &
                        !     (eta * R0(i)*patch_icpp(patch_id)%r0 & 
                        !     + (1d0-eta)*orig_prim_vf(bub_idx%rs(i)))
                        ! q_prim_vf(bub_idx%vs(i))%sf(j,k,l) = &
                        !     (eta * V0(i)*patch_icpp(patch_id)%v0 & 
                        !     + (1d0-eta)*orig_prim_vf(bub_idx%vs(i)))
                        q_prim_vf(bub_idx%rs(i))%sf(j,k,l) = muR
                        q_prim_vf(bub_idx%vs(i))%sf(j,k,l) = muV

                        IF ( .NOT. polytropic ) THEN
                            q_prim_vf(bub_idx%ps(i))%sf(j,k,l) = patch_icpp(patch_id)%p0
                            q_prim_vf(bub_idx%ms(i))%sf(j,k,l) = patch_icpp(patch_id)%m0
                        END IF

                    END IF
                END DO
            END IF

            IF (mpp_lim .AND. bubbles) THEN
                !adjust volume fractions, according to modeled gas void fraction
                alf_sum%sf = 0d0
                DO i = adv_idx%beg, adv_idx%end - 1
                    alf_sum%sf = alf_sum%sf + q_prim_vf(i)%sf
                END DO

                DO i = adv_idx%beg, adv_idx%end-1
                    q_prim_vf(i)%sf = q_prim_vf(i)%sf * (1.d0 - q_prim_vf(alf_idx)%sf) &
                            / alf_sum%sf
                END DO
            END IF

            IF (bubbles .AND. (.NOT. polytropic) ) THEN
                DO i = 1,nb
                    IF( q_prim_vf(bub_idx%ps(i))%sf(j,k,l) == dflt_real ) THEN
                        q_prim_vf(bub_idx%ps(i))%sf(j,k,l) = pb0(i) 
                        print*, 'setting to pb0'
                    END IF
                    IF( q_prim_vf(bub_idx%ms(i))%sf(j,k,l) == dflt_real ) THEN
                        q_prim_vf(bub_idx%ms(i))%sf(j,k,l) = mass_v0(i)
                    END IF
                END DO
            END IF
            

            ! Updating the patch identities bookkeeping variable
            IF(1d0 - eta < 1d-16) patch_id_fp(j,k,l) = patch_id            
            
        END SUBROUTINE s_assign_patch_species_primitive_variables_bubbles ! ------------ 
        
       
        !>      This subroutine assigns the species primitive variables
        !!              of the patch designated by the patch_id, to the cell that
        !!              is designated by the indexes (j,k,l). In addition, the
        !!              variable bookkeeping the patch identities in the entire
        !!              domain is updated with the new assignment. Note that if
        !!              the smoothing of the patch's boundaries is employed, the
        !!              ensuing primitive variables in the cell will be a type of
        !!              combination of the current patch's primitive variables
        !!              with those of the smoothing patch. The specific details
        !!              of the combination may be found in Shyue's work (1998).
        !! @param patch_id the patch identifier
        !! @param j  the x-dir node index
        !! @param k  the y-dir node index
        !! @param l  the z-dir node index
        SUBROUTINE s_assign_patch_species_primitive_variables( patch_id,j,k,l)
            
            INTEGER, INTENT(IN) :: patch_id
            INTEGER, INTENT(IN) :: j,k,l
            
            REAL(KIND(0d0)) :: rho
            REAL(KIND(0d0)) :: gamma
            REAL(KIND(0d0)) :: pi_inf
            REAL(KIND(0d0)) :: orig_rho
            REAL(KIND(0d0)) :: orig_gamma
            REAL(KIND(0d0)) :: orig_pi_inf !<
            !! Density, the specific heat ratio function and the liquid stiffness
            !! function, respectively, obtained from the combination of primitive
            !! variables of the current and smoothing patches
 
            REAL(KIND(0d0)), DIMENSION(sys_size) :: orig_prim_vf !<
            ! Vector to hold original values of cell for smoothing purposes
            
            INTEGER :: i !< generic loop iterator
            
            
            ! Transferring the identity of the smoothing patch
            smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id

            ! Transferring original primitive variables
            DO i = 1, sys_size
                orig_prim_vf(i) = q_prim_vf(i)%sf(j,k,l)
            END DO
                
            ! Computing Mixture Variables from Original Primitive Variables
            CALL s_convert_species_to_mixture_variables(    &
                                        q_prim_vf, j,k,l,   &
                                        orig_rho,           &
                                        orig_gamma,         &
                                        orig_pi_inf )

            ! Computing Mixture Variables of Current Patch =====================
            
            ! Partial densities
            DO i = 1, cont_idx%end
                q_prim_vf(i)%sf(j,k,l) = patch_icpp(patch_id)%alpha_rho(i)
            END DO
            
            ! Volume fraction(s)
            DO i = adv_idx%beg, adv_idx%end
                q_prim_vf(i)%sf(j,k,l) = patch_icpp(patch_id)%alpha(i-E_idx)
            END DO
            
            ! Density and the specific heat ratio and liquid stiffness functions
            CALL s_convert_species_to_mixture_variables(    &
                                q_prim_vf, j,k,l,           &
                                patch_icpp(patch_id)%rho,   &
                                patch_icpp(patch_id)%gamma, &
                                patch_icpp(patch_id)%pi_inf )
            
            ! ==================================================================
            
            ! Computing Mixture Variables of Smoothing Patch ===================
            
            ! Partial densities
            DO i = 1, cont_idx%end
                q_prim_vf(i)%sf(j,k,l) = &
                            patch_icpp(smooth_patch_id)%alpha_rho(i)
            END DO
            
            ! Volume fraction(s)
            DO i = adv_idx%beg, adv_idx%end
                q_prim_vf(i)%sf(j,k,l) = &
                            patch_icpp(smooth_patch_id)%alpha(i-E_idx)
            END DO
            
            ! Density and the specific heat ratio and liquid stiffness functions
            CALL s_convert_species_to_mixture_variables(       &
                            q_prim_vf, j,k,l,                  &
                            patch_icpp(smooth_patch_id)%rho,   &
                            patch_icpp(smooth_patch_id)%gamma, &
                            patch_icpp(smooth_patch_id)%pi_inf )
            
            ! ==================================================================
            
            ! Partial densities
            DO i = 1, cont_idx%end
                q_prim_vf(i)%sf(j,k,l) = &
                    eta *patch_icpp(patch_id)%alpha_rho(i) &
                + (1d0 - eta)   *orig_prim_vf(i)
            END DO
            DO i = adv_idx%beg, adv_idx%end
                q_prim_vf(i)%sf(j,k,l) = &
                    eta *patch_icpp(patch_id)%alpha(i-E_idx) &
                + (1d0 - eta)   *orig_prim_vf(i)
            END DO
            
            ! Density and the specific heat ratio and liquid stiffness functions
            CALL s_convert_species_to_mixture_variables( q_prim_vf , j,k,l, &
                                                          rho,gamma,pi_inf  )
            
            ! Velocity
            DO i = 1, E_idx-mom_idx%beg
                q_prim_vf(i+cont_idx%end)%sf(j,k,l) = &
                        (eta * patch_icpp(patch_id)%vel(i) & 
                         + (1d0-eta)*orig_prim_vf(i+cont_idx%end))
            END DO

            ! Pressure
            q_prim_vf(E_idx)%sf(j,k,l) = &
                    (eta * patch_icpp(patch_id)%pres  &
                     + (1d0 - eta) * orig_prim_vf(E_idx))

            ! Elastic Shear Stress
            IF (hypoelasticity) THEN
                DO i = 1, (stress_idx%end - stress_idx%beg) + 1
                    q_prim_vf(i+stress_idx%beg - 1)%sf(j,k,l) = &
                        (eta * patch_icpp(patch_id)%tau_e(i) &
                         + (1d0-eta)*orig_prim_vf(i+stress_idx%beg -1))
                END DO            
            END IF

            ! Set partial pressures to mixture pressure
            IF(model_eqns == 3) THEN
                DO i = internalEnergies_idx%beg, internalEnergies_idx%end
                    q_prim_vf(i)%sf(j,k,l) = q_prim_vf(E_idx)%sf(j,k,l)
                END DO
            ENDIF

            ! Updating the patch identities bookkeeping variable
            IF(1d0 - eta < 1d-16) patch_id_fp(j,k,l) = patch_id
           
            
        END SUBROUTINE s_assign_patch_species_primitive_variables ! ------------  
        
        !> Computation of parameters, allocation procedures, and/or
        !!              any other tasks needed to properly setup the module
        SUBROUTINE s_initialize_initial_condition_module() ! -------------------
            
            INTEGER :: i !< generic loop iterator
            
            
            ! Allocating the primitive and conservative variables
            ALLOCATE(q_prim_vf(1:sys_size))
            ALLOCATE(q_cons_vf(1:sys_size))
            
            DO i = 1, sys_size
                ALLOCATE(q_prim_vf(i)%sf(0:m,0:n,0:p))
                ALLOCATE(q_cons_vf(i)%sf(0:m,0:n,0:p))
            END DO
            ALLOCATE( alf_sum%sf(0:m,0:n,0:p) )
            
            
            ! Allocating the patch identities bookkeeping variable
            ALLOCATE(patch_id_fp(0:m,0:n,0:p))
            
            
            ! Setting default values for conservative and primitive variables so
            ! that in the case that the initial condition is wrongly laid out on
            ! the grid the simulation component will catch the problem on start-
            ! up. The conservative variables do not need to be similarly treated
            ! since they are computed directly from the primitive variables.
            DO i = 1, sys_size
                q_cons_vf(i)%sf = dflt_real
                q_prim_vf(i)%sf = dflt_real
            END DO
            
            
            ! Setting default values for patch identities bookkeeping variable.
            ! This is necessary to avoid any confusion in the assessment of the
            ! extent of application that the overwrite permissions give a patch
            ! when it is being applied in the domain.
            patch_id_fp = 0
            
            
            ! Depending on multicomponent flow model, the appropriate procedure
            ! for assignment of the patch mixture or species primitive variables
            ! to a cell in the domain is targeted by the procedure pointer
           
            IF(model_eqns == 1) THEN        ! Gamma/pi_inf model
                s_assign_patch_primitive_variables => &
                            s_assign_patch_mixture_primitive_variables
            ELSE IF ( bubbles ) THEN
                s_assign_patch_primitive_variables => &
                            s_assign_patch_species_primitive_variables_bubbles
            ELSE ! Volume fraction model
                s_assign_patch_primitive_variables => &
                            s_assign_patch_species_primitive_variables
            END IF
            
            
        END SUBROUTINE s_initialize_initial_condition_module ! -----------------
        
        
        
        !>  This subroutine peruses the patches and depending on the
        !!              type of geometry associated with a particular patch, it
        !!              calls the related subroutine to setup the said geometry
        !!              on the grid using the primitive variables included with
        !!              the patch parameters. The subroutine is complete once the
        !!              primitive variables are converted to conservative ones.
        SUBROUTINE s_generate_initial_condition() ! ----------------------------
            
            
            INTEGER :: i  !< Generic loop operator
            
            
            ! Converting the conservative variables to the primitive ones given
            ! preexisting initial condition data files were read in on start-up
            IF(old_ic) THEN
                CALL s_convert_conservative_to_primitive_variables( q_cons_vf, &
                                                                    q_prim_vf  )
            END IF
            
            
            !  3D Patch Geometries =============================================
            IF(p > 0) THEN
                
                DO i = 1, num_patches
                    
                    ! Spherical patch
                    IF(patch_icpp(i)%geometry == 8) THEN
                        CALL s_sphere(i)
                        
                    ! Cuboidal patch
                    ELSEIF(patch_icpp(i)%geometry == 9) THEN
                        CALL s_cuboid(i)
                        
                    ! Cylindrical patch
                    ELSEIF(patch_icpp(i)%geometry == 10) THEN
                        CALL s_cylinder(i)
                        
                    ! Swept plane patch
                    ELSEIF(patch_icpp(i)%geometry == 11) THEN
                        CALL s_sweep_plane(i)

                    ! Ellipsoidal patch
                    ELSEIF(patch_icpp(i)%geometry == 12) THEN
                        CALL s_ellipsoid(i)
                        
                    ! Analytical function patch for testing purposes
                    ELSEIF(patch_icpp(i)%geometry == 13) THEN
                        CALL s_3D_analytical(i)
                        
                    ! Spherical harmonic patch
                    ELSEIF(patch_icpp(i)%geometry == 14) THEN
                        CALL s_spherical_harmonic(i)
 
                    ! 3D Modified circular patch
                    ELSEIF(patch_icpp(i)%geometry == 19) THEN
                        CALL s_3dvarcircle(i)
                    
                    ENDIF
                    
                END DO
                
            ! ==================================================================
                
                
            ! 2D Patch Geometries ==============================================
            ELSEIF(n > 0) THEN
                
                DO i = 1, num_patches
                    
                    ! Circular patch
                    IF(patch_icpp(i)%geometry == 2) THEN
                        CALL s_circle(i)
                        
                    ! Rectangular patch
                    ELSEIF(patch_icpp(i)%geometry == 3) THEN
                        CALL s_rectangle(i)
                        
                    ! Swept line patch
                    ELSEIF(patch_icpp(i)%geometry == 4) THEN
                        CALL s_sweep_line(i)

                    ! Elliptical patch
                    ELSEIF(patch_icpp(i)%geometry == 5) THEN
                        CALL s_ellipse(i)
                        
                    ! Isentropic vortex patch
                    ELSEIF(patch_icpp(i)%geometry == 6) THEN
                        CALL s_isentropic_vortex(i)

                    ! Analytical function patch for testing purposes
                    ELSEIF(patch_icpp(i)%geometry == 7) THEN
                        CALL s_2D_analytical(i)
 
                    ! Spiral patch
                    ELSEIF(patch_icpp(i)%geometry == 17) THEN
                        CALL s_spiral(i)

                    ! Modified circular patch
                    ELSEIF(patch_icpp(i)%geometry == 18) THEN
                        CALL s_varcircle(i)
 

                    END IF
                    
                END DO
                
            ! ==================================================================
                
                
            ! 1D Patch Geometries ==============================================
            ELSE
                
                DO i = 1, num_patches
               
                    ! Line segment patch
                    IF (patch_icpp(i)%geometry == 1) THEN
                        CALL s_line_segment(i)
             
                    ! 1d analytical
                    ELSEIF(patch_icpp(i)%geometry == 15) THEN
                        CALL s_1d_analytical(i)

                    ! 1d bubble screen with sinusoidal pressure pulse
                    ELSEIF(patch_icpp(i)%geometry == 16) THEN
                        CALL s_1d_bubble_pulse(i)
                    END IF
                    
                END DO
                
            END IF
            ! ==================================================================
            
            IF (perturb_flow) CALL s_perturb_surrounding_flow()
            IF (perturb_sph) CALL s_perturb_sphere()
            
            ! Converting the primitive variables to the conservative ones
            CALL s_convert_primitive_to_conservative_variables( q_prim_vf, &
                                                                q_cons_vf  )
            
            
        END SUBROUTINE s_generate_initial_condition ! --------------------------
        
        
        SUBROUTINE s_convert_cylindrical_to_cartesian_coord(cyl_y,cyl_z)

            REAL(KIND(0d0)), INTENT(IN) :: cyl_y, cyl_z

            cart_y = cyl_y*SIN(cyl_z)
            cart_z = cyl_y*COS(cyl_z)

        END SUBROUTINE s_convert_cylindrical_to_cartesian_coord ! --------------

        SUBROUTINE s_convert_cylindrical_to_spherical_coord(cyl_x,cyl_y)

            REAL(KIND(0d0)), INTENT(IN) :: cyl_x, cyl_y

            sph_phi = ATAN(cyl_y/cyl_x)

        END SUBROUTINE s_convert_cylindrical_to_spherical_coord ! --------------


        SUBROUTINE s_perturb_sphere() ! ----------------------------------------

            INTEGER :: i,j,k,l !< generic loop operators

            REAL(KIND(0d0)) :: perturb_alpha
            REAL(KIND(0d0)) :: alpha_unadv
            REAL(KIND(0d0)) :: rand_real
            CALL RANDOM_SEED()

            DO k = 0, p
                DO j = 0, n
                    DO i = 0, m
                        CALL RANDOM_NUMBER(rand_real)
                        IF ((perturb_sph_fluid == num_fluids) .AND. (adv_alphan .NEQV. .TRUE.)) THEN
                            perturb_alpha = 1d0
                            DO l = adv_idx%beg, adv_idx%end
                                perturb_alpha = perturb_alpha - q_prim_vf(l)%sf(i,j,k)
                            END DO
                        ELSE
                            perturb_alpha = q_prim_vf(E_idx+perturb_sph_fluid)%sf(i,j,k)
                        END IF

                        ! Perturb partial density fields to match perturbed volume fraction fields
!                        IF ((perturb_alpha >= 25d-2) .AND. (perturb_alpha <= 75d-2)) THEN
                        IF ((perturb_alpha /= 0d0) .AND. (perturb_alpha /= 1d0)) THEN
                            IF (adv_alphan .NEQV. .TRUE.) THEN
                                ! Find new unadvected volume fraction
                                alpha_unadv = 1d0
                                DO l = adv_idx%beg, adv_idx%end
                                    alpha_unadv = alpha_unadv - q_prim_vf(l)%sf(i,j,k)
                                END DO
                                ! Derive new partial densities
                                DO l = 1, num_fluids-1
                                    q_prim_vf(l)%sf(i,j,k) = q_prim_vf(E_idx + l)%sf(i,j,k)*fluid_rho(l)
                                END DO
                                q_prim_vf(num_fluids)%sf(i,j,k) = alpha_unadv*fluid_rho(num_fluids)
                            ELSE
                                ! Derive new partial densities
                                DO l = 1, num_fluids
                                    q_prim_vf(l)%sf(i,j,k) = q_prim_vf(E_idx + l)%sf(i,j,k)*fluid_rho(l)
                                END DO
                            END IF
                        END IF
                    END DO
                END DO
            END DO
                            
        END SUBROUTINE s_perturb_sphere ! --------------------------------------




        SUBROUTINE s_perturb_surrounding_flow() ! ------------------------------

            INTEGER :: i,j,k,l !<  generic loop iterators

            REAL(KIND(0d0)) :: perturb_alpha
            REAL(KIND(0d0)) :: rand_real
            CALL RANDOM_SEED()

            ! Perturb partial density or velocity of surrounding flow by some random small amount of noise
            DO k = 0, p
                DO j = 0, n
                    DO i = 0, m
                        IF ((perturb_flow_fluid == num_fluids) .AND. (adv_alphan .NEQV. .TRUE.)) THEN
                            perturb_alpha = 1d0
                            DO l = adv_idx%beg, adv_idx%end
                                perturb_alpha = perturb_alpha - q_prim_vf(l)%sf(i,j,k)
                            END DO
                        ELSE
                            perturb_alpha = q_prim_vf(E_idx+perturb_flow_fluid)%sf(i,j,k)
                        END IF
                        ! IF (perturb_alpha == 1d0) THEN
                            ! Perturb partial density
!                            CALL RANDOM_NUMBER(rand_real)
!                            rand_real = rand_real / 1d2 / 1d3
!                            q_prim_vf(perturb_flow_fluid)%sf(i,j,k) = q_prim_vf(perturb_flow_fluid)%sf(i,j,k) + rand_real
                            ! Perturb velocity
                            CALL RANDOM_NUMBER(rand_real)
                            rand_real = rand_real * 1.d-2
                            q_prim_vf(mom_idx%beg)%sf(i,j,k) = (1.d0+rand_real)*q_prim_vf(mom_idx%beg)%sf(i,j,k)
                            q_prim_vf(mom_idx%end)%sf(i,j,k) = rand_real*q_prim_vf(mom_idx%beg)%sf(i,j,k)
                            IF (bubbles) THEN
                                q_prim_vf(alf_idx)%sf(i,j,k) = (1.d0+rand_real)*q_prim_vf(alf_idx)%sf(i,j,k)
                            END IF
                        ! END IF
                    END DO
                END DO
            END DO

        END SUBROUTINE s_perturb_surrounding_flow ! ----------------------------


        !>          The line segment patch is a 1D geometry that may be used,
        !!              for example, in creating a Riemann problem. The geometry
        !!              of the patch is well-defined when its centroid and length
        !!              in the x-coordinate direction are provided. Note that the
        !!              line segment patch DOES NOT allow for the smearing of its
        !!              boundaries.
        !! @param patch_id patch identifier
        SUBROUTINE s_line_segment(patch_id) ! ----------------------------------

            INTEGER, INTENT(IN) :: patch_id
            
            REAL(KIND(0d0)) :: pi_inf, gamma, lit_gamma

            INTEGER :: i,j  !< Generic loop operators

            pi_inf = fluid_pp(1)%pi_inf
            gamma = fluid_pp(1)%gamma
            lit_gamma = (1d0+gamma)/gamma
            
            ! Transferring the line segment's centroid and length information
            x_centroid = patch_icpp(patch_id)%x_centroid
            length_x   = patch_icpp(patch_id)%length_x
            
            
            ! Computing the beginning and end x-coordinates of the line segment
            ! based on its centroid and length
            x_boundary%beg = x_centroid - 0.5d0*length_x
            x_boundary%end = x_centroid + 0.5d0*length_x
            
            
            ! Since the line segment patch does not allow for its boundaries to
            ! be smoothed out, the pseudo volume fraction is set to 1 to ensure
            ! that only the current patch contributes to the fluid state in the
            ! cells that this patch covers.
            eta = 1d0
            
            
            ! Checking whether the line segment covers a particular cell in the
            ! domain and verifying whether the current patch has the permission
            ! to write to that cell. If both queries check out, the primitive
            ! variables of the current patch are assigned to this cell.
            DO i = 0, m
                IF(             x_boundary%beg <= x_cc(i) .AND.              &
                                x_boundary%end >= x_cc(i) .AND.              &
                     patch_icpp(patch_id)%alter_patch(patch_id_fp(i,0,0)) ) THEN
                    
                    CALL s_assign_patch_primitive_variables(patch_id, i,0,0)
                   
                    !IF ( (q_prim_vf(1)%sf(i,0,0) < 1.e-12) .AND. (model_eqns .NE. 4)) THEN
                    !    !zero density, reassign according to Tait EOS
                    !    q_prim_vf(1)%sf(i,0,0) = &
                    !        (((q_prim_vf(E_idx)%sf(i,0,0) + pi_inf)/(pref + pi_inf))**(1d0/lit_gamma)) * &
                    !        rhoref*(1d0-q_prim_vf(alf_idx)%sf(i,0,0))
                    !END IF
                END IF
            END DO
           
        END SUBROUTINE s_line_segment ! ----------------------------------------
       
        !>  The spiral patch is a 2D geometry that may be used, The geometry
        !!              of the patch is well-defined when its centroid and radius
        !!              are provided. Note that the circular patch DOES allow for
        !!              the smoothing of its boundary.
        !!  @param patch_id patch identifier
        SUBROUTINE s_spiral(patch_id) ! ----------------------------------------
           
            INTEGER, INTENT(IN) :: patch_id
            
            INTEGER :: i,j,k !< Generic loop iterators
            REAL(KIND(0d0)) :: th, thickness, nturns, mya
            REAL(KIND(0d0)) :: spiral_x_min, spiral_x_max, spiral_y_min, spiral_y_max
            
            ! Transferring the circular patch's radius, centroid, smearing patch
            ! identity and smearing coefficient information
            x_centroid      = patch_icpp(patch_id)%x_centroid
            y_centroid      = patch_icpp(patch_id)%y_centroid
            mya             = patch_icpp(patch_id)%radius
            thickness       = patch_icpp(patch_id)%length_x
            nturns          = patch_icpp(patch_id)%length_y

!
            logic_grid = 0
            DO k = 0,int(m*91*nturns)
                th = k/REAL(int(m*91d0*nturns)) * nturns*2.d0*pi
                
                spiral_x_min = minval( (/ f_r(th,0.0d0,    mya)*cos(th), &
                                          f_r(th,thickness,mya)*cos(th) /) )
                spiral_y_min = minval( (/ f_r(th,0.0d0,    mya)*sin(th), &
                                          f_r(th,thickness,mya)*sin(th) /) )

                spiral_x_max = maxval( (/ f_r(th,0.0d0,    mya)*cos(th), &
                                          f_r(th,thickness,mya)*cos(th) /) )
                spiral_y_max = maxval( (/ f_r(th,0.0d0,    mya)*sin(th), &
                                          f_r(th,thickness,mya)*sin(th) /) )

                DO j = 0,n; DO i = 0,m;
                    if( (x_cc(i) > spiral_x_min) .AND. (x_cc(i) < spiral_x_max) .AND. &
                        (y_cc(j) > spiral_y_min) .AND. (y_cc(j) < spiral_y_max)) THEN  
                        logic_grid(i,j,0) = 1
                    END IF
                END DO; END DO
            END DO

            DO j = 0, n
                DO i = 0, m
                    IF ( (logic_grid(i,j,0) == 1) ) THEN
                        CALL s_assign_patch_primitive_variables(patch_id, i,j,0)
                    END IF
                END DO
            END DO

        END SUBROUTINE s_spiral ! ---------------------------------------------- 


        !> Archimedes spiral function
        !! @param myth Angle
        !! @param offset Thickness
        !! @param a Starting position
        FUNCTION f_r(myth,offset,a)
            REAL(KIND(0d0)), INTENT(IN) :: myth, offset,a
            REAL(KIND(0d0)) :: b
            REAL(KIND(0d0)) :: f_r

            !r(th) = a + b*th

            b = 2.d0*a/(2.d0*pi) 
            f_r = a + b*myth + offset
        END FUNCTION f_r
        
        !> The circular patch is a 2D geometry that may be used, for
        !!              example, in creating a bubble or a droplet. The geometry
        !!              of the patch is well-defined when its centroid and radius
        !!              are provided. Note that the circular patch DOES allow for
        !!              the smoothing of its boundary.   
        !! @param patch_id is the patch identifier
        SUBROUTINE s_circle(patch_id) ! ----------------------------------------
            
            INTEGER, INTENT(IN) :: patch_id
            
            INTEGER :: i,j !< Generic loop iterators
             
            ! Transferring the circular patch's radius, centroid, smearing patch
            ! identity and smearing coefficient information
            x_centroid      = patch_icpp(patch_id)%x_centroid
            y_centroid      = patch_icpp(patch_id)%y_centroid
            radius          = patch_icpp(patch_id)%radius
            smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
            smooth_coeff    = patch_icpp(patch_id)%smooth_coeff
            
            
            ! Initializing the pseudo volume fraction value to 1. The value will
            ! be modified as the patch is laid out on the grid, but only in the
            ! case that smoothing of the circular patch's boundary is enabled.
            eta = 1d0
            
            
            ! Checking whether the circle covers a particular cell in the domain
            ! and verifying whether the current patch has permission to write to
            ! that cell. If both queries check out, the primitive variables of
            ! the current patch are assigned to this cell.
            DO j = 0, n
                DO i = 0, m
                    
                    IF(patch_icpp(patch_id)%smoothen) THEN
                        
                        eta = TANH( smooth_coeff / MIN(dx,dy) *   &
                            ( SQRT( (x_cc(i) - x_centroid)**2d0   &
                                  + (y_cc(j) - y_centroid)**2d0 ) &
                                  - radius ) ) * (-0.5d0) + 0.5d0
                        
                    END IF
                    
                    IF((       (x_cc(i) - x_centroid)**2d0                     &
                             + (y_cc(j) - y_centroid)**2d0 <= radius**2d0      &
                                                  .AND.                        &
                         patch_icpp(patch_id)%alter_patch(patch_id_fp(i,j,0))) &
                                                  .OR.                         &
                                patch_id_fp(i,j,0) == smooth_patch_id        ) &
                                                  THEN
                        
                        CALL s_assign_patch_primitive_variables(patch_id, i,j,0)
                        
                    END IF
                    
                END DO
            END DO
            
            
        END SUBROUTINE s_circle ! ----------------------------------------------
        

        !>             The varcircle patch is a 2D geometry that may be used
        !!             . It  generatres an annulus
        !! @param patch_id is the patch identifier
        SUBROUTINE s_varcircle(patch_id) ! ----------------------------------------
            
            
            ! Patch identifier
            INTEGER, INTENT(IN) :: patch_id
            
            ! Generic loop iterators
            INTEGER :: i,j
           
            REAL(KIND(0d0)) :: myr, thickness

            ! Transferring the circular patch's radius, centroid, smearing patch
            ! identity and smearing coefficient information
            x_centroid      = patch_icpp(patch_id)%x_centroid
            y_centroid      = patch_icpp(patch_id)%y_centroid
            radius          = patch_icpp(patch_id)%radius
            smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
            smooth_coeff    = patch_icpp(patch_id)%smooth_coeff
            thickness       = patch_icpp(patch_id)%epsilon
            
            
            ! Initializing the pseudo volume fraction value to 1. The value will
            ! be modified as the patch is laid out on the grid, but only in the
            ! case that smoothing of the circular patch's boundary is enabled.
            eta = 1d0
            
            
            ! Checking whether the circle covers a particular cell in the domain
            ! and verifying whether the current patch has permission to write to
            ! that cell. If both queries check out, the primitive variables of
            ! the current patch are assigned to this cell.
            DO j = 0, n
                DO i = 0, m
                    myr = dsqrt(  (x_cc(i) - x_centroid)**2d0 &
                                + (y_cc(j) - y_centroid)**2d0                   )

                    IF( myr <= radius+thickness/2.d0 .AND. &       
                        myr >= radius-thickness/2.d0 .AND. &
                        patch_icpp(patch_id)%alter_patch(patch_id_fp(i,j,0)) ) THEN
                        
                        CALL s_assign_patch_primitive_variables(patch_id, i,j,0)
                        
                        q_prim_vf(alf_idx)%sf(i,j,0) = patch_icpp(patch_id)%alpha(1) * &
                            dexp( -0.5d0 * ((myr - radius)**2.d0) / (thickness/3.d0)**2.d0 )
                    END IF
                    
                END DO
            END DO
            
            
        END SUBROUTINE s_varcircle ! ----------------------------------------------
       

        SUBROUTINE s_3dvarcircle(patch_id) ! ----------------------------------------
            
            ! Patch identifier
            INTEGER, INTENT(IN) :: patch_id
            
            ! Generic loop iterators
            INTEGER :: i,j,k
           
            REAL(KIND(0d0)) :: myr, thickness

            ! Transferring the circular patch's radius, centroid, smearing patch
            ! identity and smearing coefficient information
            x_centroid      = patch_icpp(patch_id)%x_centroid
            y_centroid      = patch_icpp(patch_id)%y_centroid
            z_centroid      = patch_icpp(patch_id)%z_centroid
            length_z        = patch_icpp(patch_id)%length_z 
            radius          = patch_icpp(patch_id)%radius
            smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
            smooth_coeff    = patch_icpp(patch_id)%smooth_coeff
            thickness       = patch_icpp(patch_id)%epsilon
            
            
            ! Initializing the pseudo volume fraction value to 1. The value will
            ! be modified as the patch is laid out on the grid, but only in the
            ! case that smoothing of the circular patch's boundary is enabled.
            eta = 1d0
            
            ! write for all z
 
            ! Checking whether the circle covers a particular cell in the domain
            ! and verifying whether the current patch has permission to write to
            ! that cell. If both queries check out, the primitive variables of
            ! the current patch are assigned to this cell.
            DO k =0,p
            DO j = 0, n
                DO i = 0, m
                    myr = dsqrt(  (x_cc(i) - x_centroid)**2d0 &
                                + (y_cc(j) - y_centroid)**2d0                   )

                    IF( myr <= radius+thickness/2.d0 .AND. &       
                        myr >= radius-thickness/2.d0 .AND. &
                        patch_icpp(patch_id)%alter_patch(patch_id_fp(i,j,k)) ) THEN
                        
                        CALL s_assign_patch_primitive_variables(patch_id, i,j,k)
                        
                        q_prim_vf(alf_idx)%sf(i,j,k) = patch_icpp(patch_id)%alpha(1) * &
                            dexp( -0.5d0 * ((myr - radius)**2.d0) / (thickness/3.d0)**2.d0 )
                    END IF
                    
                END DO
            END DO
            END DO
            
            
        END SUBROUTINE s_3dvarcircle ! ----------------------------------------------
        
        !>      The elliptical patch is a 2D geometry. The geometry of 
        !!      the patch is well-defined when its centroid and radii
        !!      are provided. Note that the elliptical patch DOES allow
        !!      for the smoothing of its boundary
        !! @param patch_id is the patch identifier
        SUBROUTINE s_ellipse(patch_id) ! ---------------------------------------

            INTEGER, INTENT(IN) :: patch_id

            INTEGER :: i,j !< Generic loop operators

            ! Transferring the elliptical patch's radii, centroid, smearing
            ! patch identity, and smearing coefficient information
            x_centroid = patch_icpp(patch_id)%x_centroid
            y_centroid = patch_icpp(patch_id)%y_centroid
            a = patch_icpp(patch_id)%radii(1)
            b = patch_icpp(patch_id)%radii(2)
            smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
            smooth_coeff = patch_icpp(patch_id)%smooth_coeff

            ! Initializing the pseudo volume fraction value to 1. The value 
            ! be modified as the patch is laid out on the grid, but only in 
            ! the case that smoothing of the elliptical patch's boundary is
            ! enabled. 
            eta = 1d0

            ! Checking whether the ellipse covers a particular cell in the 
            ! domain and verifying whether the current patch has permission
            ! to write to that cell. If both queries check out, the primitive
            ! variables of the current patch are assigned to this cell.
            DO j = 0, n
                DO i = 0, m

                    IF(patch_icpp(patch_id)%smoothen) THEN
                        eta = TANH(smooth_coeff / MIN(dx,dy) * &
                            ( SQRT( ((x_cc(i) - x_centroid) / a)**2d0 + &
                                    ((y_cc(j) - y_centroid) / b)**2d0 ) &
                                    - 1d0 )) * (-0.5d0) + 0.5d0
                    END IF

                    IF ((        ((x_cc(i) - x_centroid) / a)**2d0 +          &
                                 ((y_cc(j) - y_centroid) / b)**2d0 <= 1d0     &
                                                .AND.                         &
                        patch_icpp(patch_id)%alter_patch(patch_id_fp(i,j,0))) &
                                                .OR.                          &
                              patch_id_fp(i,j,0) == smooth_patch_id )         &
                                                THEN

                        CALL s_assign_patch_primitive_variables(patch_id, i,j,0)
                    END IF
                END DO
            END DO

        END SUBROUTINE s_ellipse ! ---------------------------------------------



        !>      The ellipsoidal patch is a 3D geometry. The geometry of 
        !!       the patch is well-defined when its centroid and radii
        !!       are provided. Note that the ellipsoidal patch DOES allow
        !!       for the smoothing of its boundary
        !! @param patch_id is the patch identifier
        SUBROUTINE s_ellipsoid(patch_id) ! -------------------------------------


            ! Patch identifier
            INTEGER, INTENT(IN) :: patch_id

            ! Generic loop iterators
            INTEGER :: i,j,k

            ! Transferring the ellipsoidal patch's radii, centroid, smearing
            ! patch identity, and smearing coefficient information
            x_centroid = patch_icpp(patch_id)%x_centroid
            y_centroid = patch_icpp(patch_id)%y_centroid
            z_centroid = patch_icpp(patch_id)%z_centroid
            a = patch_icpp(patch_id)%radii(1)
            b = patch_icpp(patch_id)%radii(2)
            c = patch_icpp(patch_id)%radii(3)
            smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
            smooth_coeff = patch_icpp(patch_id)%smooth_coeff

            ! Initializing the pseudo volume fraction value to 1. The value 
            ! be modified as the patch is laid out on the grid, but only in 
            ! the case that smoothing of the ellipsoidal patch's boundary is
            ! enabled. 
            eta = 1d0

            ! Checking whether the ellipsoid covers a particular cell in the 
            ! domain and verifying whether the current patch has permission
            ! to write to that cell. If both queries check out, the primitive
            ! variables of the current patch are assigned to this cell.
            DO k = 0, p
               DO j = 0, n
                  DO i = 0, m

                     IF (grid_geometry == 3) THEN
                         CALL s_convert_cylindrical_to_cartesian_coord(y_cc(j),z_cc(k))
                     ELSE
                         cart_y = y_cc(j)
                         cart_z = z_cc(k)
                     END IF

                     IF(patch_icpp(patch_id)%smoothen) THEN
                         eta = TANH(smooth_coeff / MIN(dx,dy,dz) *       &
                             ( SQRT( ((x_cc(i) - x_centroid) / a)**2d0 + &
                                     (( cart_y - y_centroid) / b)**2d0 + &
                                     (( cart_z - z_centroid) / c)**2d0 ) &
                                 - 1d0 )) * (-0.5d0) + 0.5d0
                     END IF

                     IF ((      ((  x_cc(i) - x_centroid) / a)**2d0 +          &
                                ((   cart_y - y_centroid) / b)**2d0 +          &
                                ((   cart_z - z_centroid) / c)**2d0 <= 1d0     &
                                                 .AND.                         &
                         patch_icpp(patch_id)%alter_patch(patch_id_fp(i,j,k))) &
                                                 .OR.                          &
                         patch_id_fp(i,j,k) == smooth_patch_id )               &
                                                 THEN

                         CALL s_assign_patch_primitive_variables(patch_id, i,j,k)
                     END IF
                  END DO
               END DO
            END DO

        END SUBROUTINE s_ellipsoid ! -------------------------------------------


        !>      The rectangular patch is a 2D geometry that may be used,
        !!              for example, in creating a solid boundary, or pre-/post-
        !!              shock region, in alignment with the axes of the Cartesian
        !!              coordinate system. The geometry of such a patch is well-
        !!              defined when its centroid and lengths in the x- and y-
        !!              coordinate directions are provided. Please note that the
        !!              rectangular patch DOES NOT allow for the smoothing of its
        !!              boundaries.
        !! @param patch_id is the patch identifier
        SUBROUTINE s_rectangle(patch_id) ! -------------------------------------
           
            INTEGER, INTENT(IN) :: patch_id 

            REAL(KIND(0d0)) :: pi_inf, gamma, lit_gamma !< Equation of state parameters

            INTEGER :: i,j !< generic loop iterators

            pi_inf = fluid_pp(1)%pi_inf
            gamma = fluid_pp(1)%gamma
            lit_gamma = (1d0+gamma)/gamma

            
            
            ! Transferring the rectangle's centroid and length information
            x_centroid = patch_icpp(patch_id)%x_centroid
            y_centroid = patch_icpp(patch_id)%y_centroid
            length_x   = patch_icpp(patch_id)%length_x
            length_y   = patch_icpp(patch_id)%length_y
            
            
            ! Computing the beginning and the end x- and y-coordinates of the
            ! rectangle based on its centroid and lengths
            x_boundary%beg = x_centroid - 0.5d0*length_x
            x_boundary%end = x_centroid + 0.5d0*length_x
            y_boundary%beg = y_centroid - 0.5d0*length_y
            y_boundary%end = y_centroid + 0.5d0*length_y
            
            
            ! Since the rectangular patch does not allow for its boundaries to
            ! be smoothed out, the pseudo volume fraction is set to 1 to ensure
            ! that only the current patch contributes to the fluid state in the
            ! cells that this patch covers.
            eta = 1d0
            
            
            ! Checking whether the rectangle covers a particular cell in the
            ! domain and verifying whether the current patch has the permission
            ! to write to that cell. If both queries check out, the primitive
            ! variables of the current patch are assigned to this cell.
            DO j = 0, n
                DO i = 0, m
                    IF(             x_boundary%beg <= x_cc(i) .AND.            &
                                    x_boundary%end >= x_cc(i) .AND.            &
                                    y_boundary%beg <= y_cc(j) .AND.            &
                                    y_boundary%end >= y_cc(j)                  &
                                                .AND.                          &
                        patch_icpp(patch_id)%alter_patch(patch_id_fp(i,j,0)) ) &
                                                THEN
                        
                        CALL s_assign_patch_primitive_variables(patch_id, i,j,0)

                        IF ( (q_prim_vf(1)%sf(i,j,0) < 1.e-10) .AND. (model_eqns == 4)) THEN
                            !zero density, reassign according to Tait EOS
                            q_prim_vf(1)%sf(i,j,0) = &
                                (((q_prim_vf(E_idx)%sf(i,j,0) + pi_inf)/(pref + pi_inf))**(1d0/lit_gamma)) * &
                                rhoref*(1d0-q_prim_vf(alf_idx)%sf(i,j,0))
                        END IF                        
                    END IF
                END DO
            END DO
            
            
        END SUBROUTINE s_rectangle ! -------------------------------------------
        
        
        !>  The swept line patch is a 2D geometry that may be used,
        !!      for example, in creating a solid boundary, or pre-/post-
        !!      shock region, at an angle with respect to the axes of the
        !!      Cartesian coordinate system. The geometry of the patch is
        !!      well-defined when its centroid and normal vector, aimed
        !!      in the sweep direction, are provided. Note that the sweep
        !!      line patch DOES allow the smoothing of its boundary.        
        !! @param patch_id is the patch identifier
        SUBROUTINE s_sweep_line(patch_id) ! ------------------------------------

            INTEGER, INTENT(IN) :: patch_id
            
            INTEGER :: i,j !< Generic loop operators
            
            
            ! Transferring the centroid information of the line to be swept
            x_centroid = patch_icpp(patch_id)%x_centroid
            y_centroid = patch_icpp(patch_id)%y_centroid
            smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
            smooth_coeff    = patch_icpp(patch_id)%smooth_coeff
            
            
            ! Obtaining coefficients of the equation describing the sweep line
            a = patch_icpp(patch_id)%normal(1)
            b = patch_icpp(patch_id)%normal(2)
            c = -a*x_centroid - b*y_centroid
            
            
            ! Initializing the pseudo volume fraction value to 1. The value will
            ! be modified as the patch is laid out on the grid, but only in the
            ! case that smoothing of the sweep line patch's boundary is enabled.
            eta = 1d0
            
            ! Checking whether the region swept by the line covers a particular
            ! cell in the domain and verifying whether the current patch has the
            ! permission to write to that cell. If both queries check out, the
            ! primitive variables of the current patch are written to this cell.
            DO j = 0, n
                DO i = 0, m
                    
                    IF(patch_icpp(patch_id)%smoothen) THEN
                        eta = 5d-1 + 5d-1*TANH( smooth_coeff/MIN(dx,dy)  &
                                               * (a*x_cc(i)+b*y_cc(j)+c) &
                                               / SQRT(a**2d0+b**2d0)     )
                    END IF
                    
                    IF((           a*x_cc(i) + b*y_cc(j) + c >= 0d0            &
                                                 .AND.                         &
                         patch_icpp(patch_id)%alter_patch(patch_id_fp(i,j,0))) &
                                                 .OR.                          &
                                patch_id_fp(i,j,0) == smooth_patch_id        ) &
                                                 THEN
                        CALL s_assign_patch_primitive_variables(patch_id, i,j,0)
                        
                    END IF
                    
                END DO
            END DO
            
            
        END SUBROUTINE s_sweep_line ! ------------------------------------------
        
        
        !> The isentropic vortex is a 2D geometry that may be used,
        !!              for example, to generate an isentropic flow disturbance.
        !!              Geometry of the patch is well-defined when its centroid
        !!              and radius are provided. Notice that the patch DOES NOT
        !!              allow for the smoothing of its boundary.
        !! @param patch_id is the patch identifier
        SUBROUTINE s_isentropic_vortex(patch_id) ! ----------------------------
           
            
            ! Patch identifier
            INTEGER, INTENT(IN) :: patch_id
            
            ! Generic loop iterators
            INTEGER :: i,j
            
            
            ! Transferring isentropic vortex patch's centroid and radius info
            x_centroid      = patch_icpp(patch_id)%x_centroid
            y_centroid      = patch_icpp(patch_id)%y_centroid
            radius          = patch_icpp(patch_id)%radius
            
            
            ! Since the isentropic vortex patch does not allow for its boundary
            ! to get smoothed, the pseudo volume fraction is set to 1 to ensure
            ! that only the current patch contributes to the fluid state in the
            ! cells that this patch covers.
            eta = 1d0
            
            
            ! Verifying whether the isentropic vortex includes a particular cell
            ! and verifying whether the current patch has permission to write to
            ! that cell. If both queries work out the primitive variables of the
            ! the current patch are assigned to this cell.
            DO j = 0, n
                DO i = 0, m
                    
                    IF(        (x_cc(i) - x_centroid)**2d0                     &
                             + (y_cc(j) - y_centroid)**2d0 <= radius**2d0      &
                                                  .AND.                        &
                        patch_icpp(patch_id)%alter_patch(patch_id_fp(i,j,0)) ) &
                                                  THEN
                        
                        CALL s_assign_patch_primitive_variables(  patch_id, &
                                                                     i,j,0)
                        
                    END IF
                    
                END DO
            END DO
            
            
        END SUBROUTINE s_isentropic_vortex ! -----------------------------------
        
        !>  This patch assigns the primitive variables as analytical 
        !!  functions such that the code can be verified.         
        !!  @param patch_id is the patch identifier
        SUBROUTINE s_1D_analytical(patch_id) ! ---------------------------------

            ! Patch identifier
            INTEGER, INTENT(IN) :: patch_id

            ! Placeholders for the cell boundary values
            REAL(KIND(0d0)) :: a,b,c,d, pi_inf, gamma, lit_gamma

            ! Generic loop iterators
            INTEGER :: i,j

            pi_inf = fluid_pp(1)%pi_inf
            gamma = fluid_pp(1)%gamma
            lit_gamma = (1d0+gamma)/gamma

            ! Transferring the patch's centroid and length information
            x_centroid = patch_icpp(patch_id)%x_centroid
            length_x = patch_icpp(patch_id)%length_x

            ! Computing the beginning and the end x- and y-coordinates
            ! of the patch based on its centroid and lengths
            x_boundary%beg = x_centroid - 0.5d0*length_x
            x_boundary%end = x_centroid + 0.5d0*length_x

            ! Since the patch doesn't allow for its boundaries to be
            ! smoothed out, the pseudo volume fraction is set to 1 to 
            ! ensure that only the current patch contributes to the fluid
            ! state in the cells that this patch covers.
            eta = 1d0
            
            ! Checking whether the line segment covers a particular cell in the
            ! domain and verifying whether the current patch has the permission
            ! to write to that cell. If both queries check out, the primitive
            ! variables of the current patch are assigned to this cell.
            DO i = 0, m
                IF(             x_boundary%beg <= x_cc(i) .AND.              &
                                x_boundary%end >= x_cc(i) .AND.              &
                     patch_icpp(patch_id)%alter_patch(patch_id_fp(i,0,0)) ) THEN
                   
                    CALL s_assign_patch_primitive_variables(patch_id, i,0,0)

                    !what variables to alter
                    !bump in pressure
                    q_prim_vf(E_idx)%sf(i,0,0) = q_prim_vf(E_idx)%sf(i,0,0) * &
                        ( 1d0 + 0.2d0*dexp(-1d0*((x_cb(i)-x_centroid)**2.d0)/(2.d0*0.005d0)) )

                    !bump in void fraction
                    !q_prim_vf(adv_idx%beg)%sf(i,0,0) = q_prim_vf(adv_idx%beg)%sf(i,0,0) * &
                    !    ( 1d0 + 0.2d0*exp(-1d0*((x_cb(i)-x_centroid)**2.d0)/(2.d0*0.005d0)) )

                    !bump in R(x)
                    !q_prim_vf(adv_idx%end+1)%sf(i,0,0) = q_prim_vf(adv_idx%end+1)%sf(i,0,0) * &
                    !    ( 1d0 + 0.2d0*exp(-1d0*((x_cb(i)-x_centroid)**2.d0)/(2.d0*0.005d0)) )

                    !IF (model_eqns == 4) THEN
                        !reassign density
                    !IF (num_fluids == 1) THEN
                        q_prim_vf(1)%sf(i,0,0) = &
                            (((q_prim_vf(E_idx)%sf(i,0,0) + pi_inf)/(pref + pi_inf))**(1d0/lit_gamma)) * &
                            rhoref*(1d0-q_prim_vf(alf_idx)%sf(i,0,0))
                    !END IF
                    !ELSE IF (model_eqns == 2) THEN
                        !can manually adjust density here
                        !q_prim_vf(1)%sf(i,0,0) = q_prim_vf(1)%sf(i,0,0) * &
                        !    ( 1d0 + 0.2d0*exp(-1d0*((x_cb(i)-x_centroid)**2.d0)/(2.d0*0.005d0)) )
                    !END IF
                END IF
            END DO

            
        END SUBROUTINE s_1D_analytical ! --------------------------------------- 
        

 
        SUBROUTINE s_1d_bubble_pulse(patch_id) ! ---------------------------------
        ! Description: This patch assigns the primitive variables as analytical 
        !       functions such that the code can be verified. 

            ! Patch identifier
            INTEGER, INTENT(IN) :: patch_id

            ! Placeholders for the cell boundary values
            REAL(KIND(0d0)) :: fac, a,b,c,d, pi_inf, gamma, lit_gamma

            ! Generic loop iterators
            INTEGER :: i,j

            pi_inf = fluid_pp(1)%pi_inf
            gamma = fluid_pp(1)%gamma
            lit_gamma = (1d0+gamma)/gamma

            ! Transferring the patch's centroid and length information
            x_centroid = patch_icpp(patch_id)%x_centroid
            length_x = patch_icpp(patch_id)%length_x


            ! Computing the beginning and the end x- and y-coordinates
            ! of the patch based on its centroid and lengths
            x_boundary%beg = x_centroid - 0.5d0*length_x
            x_boundary%end = x_centroid + 0.5d0*length_x

            ! Since the patch doesn't allow for its boundaries to be
            ! smoothed out, the pseudo volume fraction is set to 1 to 
            ! ensure that only the current patch contributes to the fluid
            ! state in the cells that this patch covers.
            eta = 1d0
            
            ! Checking whether the line segment covers a particular cell in the
            ! domain and verifying whether the current patch has the permission
            ! to write to that cell. If both queries check out, the primitive
            ! variables of the current patch are assigned to this cell.
            DO i = 0, m
                IF(             x_boundary%beg <= x_cc(i) .AND.              &
                                x_boundary%end >= x_cc(i) .AND.              &
                     patch_icpp(patch_id)%alter_patch(patch_id_fp(i,0,0)) ) THEN
                   
                    CALL s_assign_patch_primitive_variables(patch_id, i,0,0)

                    !what variables to alter
                    !sinusoid in pressure
                    q_prim_vf(E_idx)%sf(i,0,0) = q_prim_vf(E_idx)%sf(i,0,0) * &
                        ( 1d0 + 0.1d0*sin(-1d0*(x_cb(i)-x_centroid)*2d0*pi/length_x ) )

                    !bump in void fraction
                    !q_prim_vf(adv_idx%beg)%sf(i,0,0) = q_prim_vf(adv_idx%beg)%sf(i,0,0) * &
                    !    ( 1d0 + 0.2d0*exp(-1d0*((x_cb(i)-x_centroid)**2.d0)/(2.d0*0.005d0)) )

                    !bump in R(x)
                    !q_prim_vf(adv_idx%end+1)%sf(i,0,0) = q_prim_vf(adv_idx%end+1)%sf(i,0,0) * &
                    !    ( 1d0 + 0.2d0*exp(-1d0*((x_cb(i)-x_centroid)**2.d0)/(2.d0*0.005d0)) )

                    !IF (model_eqns == 4) THEN
                        !reassign density
                    !IF (num_fluids == 1) THEN
                        q_prim_vf(1)%sf(i,0,0) = &
                            (((q_prim_vf(E_idx)%sf(i,0,0) + pi_inf)/(pref + pi_inf))**(1d0/lit_gamma)) * &
                            rhoref*(1d0-q_prim_vf(alf_idx)%sf(i,0,0))
                    !END IF
                    !ELSE IF (model_eqns == 2) THEN
                        !can manually adjust density here
                        !q_prim_vf(1)%sf(i,0,0) = q_prim_vf(1)%sf(i,0,0) * &
                        !    ( 1d0 + 0.2d0*exp(-1d0*((x_cb(i)-x_centroid)**2.d0)/(2.d0*0.005d0)) )
                    !END IF
                END IF
            END DO

            
        END SUBROUTINE s_1D_bubble_pulse ! --------------------------------------- 

        !>  This patch assigns the primitive variables as analytical 
        !!  functions such that the code can be verified. 
        !!  @param patch_id is the patch identifier
        SUBROUTINE s_2D_analytical(patch_id) ! ---------------------------------

            INTEGER, INTENT(IN) :: patch_id

            REAL(KIND(0d0)) :: a,b,c,d !< placeholderrs for the cell boundary values
            REAL(KIND(0d0)) :: pi_inf,gamma,lit_gamma !< equation of state parameters

            INTEGER :: i,j !< generic loop iterators

            pi_inf = fluid_pp(1)%pi_inf
            gamma = fluid_pp(1)%gamma
            lit_gamma = (1d0+gamma)/gamma
            
            ! Transferring the patch's centroid and length information
            x_centroid = patch_icpp(patch_id)%x_centroid
            y_centroid = patch_icpp(patch_id)%y_centroid
            length_x = patch_icpp(patch_id)%length_x
            length_y = patch_icpp(patch_id)%length_y

            ! Computing the beginning and the end x- and y-coordinates
            ! of the patch based on its centroid and lengths
            x_boundary%beg = x_centroid - 0.5d0*length_x
            x_boundary%end = x_centroid + 0.5d0*length_x
            y_boundary%beg = y_centroid - 0.5d0*length_y
            y_boundary%end = y_centroid + 0.5d0*length_y

            ! Since the patch doesn't allow for its boundaries to be
            ! smoothed out, the pseudo volume fraction is set to 1 to 
            ! ensure that only the current patch contributes to the fluid
            ! state in the cells that this patch covers.
            eta = 1d0

            ! Checking whether the patch covers a particular cell in the 
            ! domain and verifying whether the current patch has the 
            ! permission to write to that cell. If both queries check out, 
            ! the primitive variables of the current patch are assigned
            ! to this cell.
            DO j = 0, n
                DO i = 0, m
                    IF (    x_boundary%beg <= x_cc(i) .AND. &
                            x_boundary%end >= x_cc(i) .AND. &
                            y_boundary%beg <= y_cc(j) .AND. &
                            y_boundary%end >= y_cc(j) .AND. &
                            patch_icpp(patch_id)%alter_patch(patch_id_fp(i,j,0))) THEN
                    
                        CALL s_assign_patch_primitive_variables(patch_id, i,j,0)

                        !what variables to alter
                        !x-y bump in pressure
                        !q_prim_vf(E_idx)%sf(i,j,0) = q_prim_vf(E_idx)%sf(i,j,0) * &
                        !    ( 1d0 + 0.2d0*dexp(-1d0*((x_cb(i)-x_centroid)**2.d0 + (y_cb(j)-y_centroid)**2.d0)/(2.d0*0.005d0)) )

                        !x-bump
                        q_prim_vf(E_idx)%sf(i,j,0) = q_prim_vf(E_idx)%sf(i,j,0) * &
                            ( 1d0 + 0.2d0*dexp(-1d0*((x_cb(i)-x_centroid)**2.d0)/(2.d0*0.005d0)) )

                        !bump in void fraction
                        !q_prim_vf(adv_idx%beg)%sf(i,j,0) = q_prim_vf(adv_idx%beg)%sf(i,j,0) * &
                        !    ( 1d0 + 0.2d0*exp(-1d0*((x_cb(i)-x_centroid)**2.d0 + (y_cb(j)-y_centroid)**2.d0)/(2.d0*0.005d0)) )

                        !bump in R(x)
                        !q_prim_vf(adv_idx%end+1)%sf(i,j,0) = q_prim_vf(adv_idx%end+1)%sf(i,j,0) * &
                        !    ( 1d0 + 0.2d0*exp(-1d0*((x_cb(i)-x_centroid)**2.d0 + (y_cb(j)-y_centroid)**2.d0)/(2.d0*0.005d0)) )

                        !reassign density
                        q_prim_vf(1)%sf(i,j,0) = &
                            (((q_prim_vf(E_idx)%sf(i,j,0) + pi_inf)/(pref + pi_inf))**(1d0/lit_gamma)) * &
                            rhoref*(1d0-q_prim_vf(alf_idx)%sf(i,j,0))
                       
                       ! ================================================================================

                        ! Sinusoidal initial condition for all flow variables =============================

                        ! Cell-center values
!                        a = 0d0
!                        b = 0d0
!                        c = 0d0
!                        d = 0d0
!                        q_prim_vf(adv_idx%beg)%sf(i,j,0) = SIN(x_cc(i)) * SIN(y_cc(j))
!                        q_prim_vf(1)%sf(i,j,0) = q_prim_vf(adv_idx%beg)%sf(i,j,0) * 1d0
!                        q_prim_vf(cont_idx%end)%sf(i,j,0) = (1d0 - q_prim_vf(adv_idx%beg)%sf(i,j,0)) * 1d0
!                        q_prim_vf(mom_idx%beg)%sf(i,j,0) = SIN(x_cc(i))
!                        q_prim_vf(mom_idx%end)%sf(i,j,0) = SIN(y_cc(j))
!                        q_prim_vf(E_idx)%sf(i,j,0) = 1d0

                        ! Cell-average values
!                       a = x_cc(i) - 5d-1*dx ! x-beg
!                       b = x_cc(i) + 5d-1*dx ! x-end
!                       c = y_cc(j) - 5d-1*dy ! y-beg
!                       d = y_cc(j) + 5d-1*dy ! y-end
!                       q_prim_vf(adv_idx%beg)%sf(i,j,0) = 1d0/((b-a)*(d-c)) * &
!                               (COS(a)*COS(c) - COS(a)*COS(d) - COS(b)*COS(c) + COS(b)*COS(d))
!                       q_prim_vf(1)%sf(i,j,0) = q_prim_vf(adv_idx%beg)%sf(i,j,0) * 1d0
!                       q_prim_vf(cont_idx%end)%sf(i,j,0) = (1d0 - q_prim_vf(adv_idx%beg)%sf(i,j,0)) * 1d0
!                       q_prim_vf(mom_idx%beg)%sf(i,j,0) = (COS(a) - COS(b))/(b-a)
!                       q_prim_vf(mom_idx%end)%sf(i,j,0) = (COS(c) - COS(d))/(d-c)
!                       q_prim_vf(E_idx)%sf(i,j,0) = 1d0
                        ! ================================================================================

                        ! Initial pressure profile smearing for bubble collapse case of Tiwari (2013) ====
                        !IF((       (x_cc(i))**2d0                     &
                        !         + (y_cc(j))**2d0 <= 1d0**2d0)) THEN
                        !         q_prim_vf(E_idx)%sf(i,j,0) = 1d5 / 25d0
                        !ELSE
                        !    q_prim_vf(E_idx)%sf(i,j,0) = 1d5 + 1d0/SQRT(x_cc(i)**2d0+y_cc(j)**2d0) & 
                        !                                    * ((1d5/25d0) - 1d5)
                        !END IF
                        ! ================================================================================

                    END IF
                END DO
            END DO

        END SUBROUTINE s_2D_analytical ! ---------------------------------------



        !>      This patch assigns the primitive variables as analytical 
        !!      functions such that the code can be verified. 
        !!      @param patch_id is the patch identifier
        SUBROUTINE s_3D_analytical(patch_id) ! ---------------------------------


            INTEGER, INTENT(IN) :: patch_id
            REAL(KIND(0d0)) :: pi_inf, gamma, lit_gamma !< equation of state parameters
            
            INTEGER :: i,j,k !< generic loop iterators

            pi_inf = fluid_pp(1)%pi_inf
            gamma = fluid_pp(1)%gamma
            lit_gamma = (1d0+gamma)/gamma    

            ! Transferring the patch's centroid and length information
            x_centroid = patch_icpp(patch_id)%x_centroid
            y_centroid = patch_icpp(patch_id)%y_centroid
            z_centroid = patch_icpp(patch_id)%z_centroid
            length_x   = patch_icpp(patch_id)%length_x
            length_y   = patch_icpp(patch_id)%length_y
            length_z   = patch_icpp(patch_id)%length_z
            
            ! Computing the beginning and the end x-, y- and z-coordinates of
            ! the patch based on its centroid and lengths
            x_boundary%beg = x_centroid - 0.5d0*length_x
            x_boundary%end = x_centroid + 0.5d0*length_x
            y_boundary%beg = y_centroid - 0.5d0*length_y
            y_boundary%end = y_centroid + 0.5d0*length_y
            z_boundary%beg = z_centroid - 0.5d0*length_z
            z_boundary%end = z_centroid + 0.5d0*length_z
            
            ! Since the analytical patch does not allow for its boundaries to get
            ! smoothed out, the pseudo volume fraction is set to 1 to make sure
            ! that only the current patch contributes to the fluid state in the
            ! cells that this patch covers.
            eta = 1d0
            
            ! Checking whether the patch covers a particular cell in the domain
            ! and verifying whether the current patch has permission to write to
            ! to that cell. If both queries check out, the primitive variables
            ! of the current patch are assigned to this cell.
            DO k = 0, p
               DO j = 0, n
                  DO i = 0, m

                     IF (grid_geometry == 3) THEN
                         CALL s_convert_cylindrical_to_cartesian_coord(y_cc(j),z_cc(k))
                     ELSE
                         cart_y = y_cc(j)
                         cart_z = z_cc(k)
                     END IF

                     IF(            x_boundary%beg <= x_cc(i)   .AND.            &
                                    x_boundary%end >= x_cc(i)   .AND.            &
                                    y_boundary%beg <= cart_y    .AND.            &
                                    y_boundary%end >= cart_y    .AND.            &
                                    z_boundary%beg <= cart_z    .AND.            &
                                    z_boundary%end >= cart_z                     &
                                                .AND.                            &
                         patch_icpp(patch_id)%alter_patch(patch_id_fp(i,j,k)))   &
                                                THEN

                        CALL s_assign_patch_primitive_variables(patch_id, i,j,k)

                        !gaussian ball
                        !what variables to alter
                        !bump in pressure
                        q_prim_vf(E_idx)%sf(i,j,k) = q_prim_vf(E_idx)%sf(i,j,k) * &
                            ( 1d0 + 0.2d0*exp(-1d0 * &
                            ((x_cb(i)-x_centroid)**2.d0 + (y_cb(j)-y_centroid)**2.d0 + (z_cb(k)-z_centroid)**2.d0) &
                            /(2.d0*0.005d0)) )

                        !bump in void fraction
                        q_prim_vf(adv_idx%beg)%sf(i,j,k) = q_prim_vf(adv_idx%beg)%sf(i,j,k) * &
                            ( 1d0 + 0.2d0*exp(-1d0 * &
                            ((x_cb(i)-x_centroid)**2.d0 + (y_cb(j)-y_centroid)**2.d0 + (z_cb(k)-z_centroid)**2.d0) &
                            /(2.d0*0.005d0)) )

                        !bump in R(x)
                        q_prim_vf(adv_idx%end+1)%sf(i,j,k) = q_prim_vf(adv_idx%end+1)%sf(i,j,k) * &
                            ( 1d0 + 0.2d0*exp(-1d0 * &
                            ((x_cb(i)-x_centroid)**2.d0 + (y_cb(j)-y_centroid)**2.d0 + (z_cb(k)-z_centroid)**2.d0) &
                            /(2.d0*0.005d0)) )

                        !reassign density
                        q_prim_vf(1)%sf(i,j,k) = &
                            (((q_prim_vf(E_idx)%sf(i,j,k) + pi_inf)/(pref + pi_inf))**(1d0/lit_gamma)) * &
                            rhoref*(1d0-q_prim_vf(E_idx+1)%sf(i,j,k))

                        ! ================================================================================

                        ! Constant x-velocity in cylindrical grid ========================================
!                        q_prim_vf(cont_idx%beg )%sf(i,j,k) = 1d0 
!                        q_prim_vf(cont_idx%end )%sf(i,j,k) = 0d0
!                        q_prim_vf(mom_idx%beg  )%sf(i,j,k) = 0d0
!                        q_prim_vf(mom_idx%beg+1)%sf(i,j,k) = COS(z_cc(k))
!                        q_prim_vf(mom_idx%end  )%sf(i,j,k) = -SIN(z_cc(k))
!                        q_prim_vf(E_idx        )%sf(i,j,k) = 1d0
!                        q_prim_vf(adv_idx%beg  )%sf(i,j,k) = 1d0
                        ! ================================================================================

                        ! Couette flow in cylindrical grid ===============================================
                        !q_prim_vf(cont_idx%beg )%sf(i,j,k) = 1d0 
                        !q_prim_vf(cont_idx%end )%sf(i,j,k) = 0d0
                        !q_prim_vf(mom_idx%beg  )%sf(i,j,k) = 0d0
                        !q_prim_vf(mom_idx%beg+1)%sf(i,j,k) = y_cc(j)*COS(z_cc(k))*SIN(z_cc(k))
                        !q_prim_vf(mom_idx%end  )%sf(i,j,k) = -y_cc(j)*SIN(z_cc(k))**2d0
                        !q_prim_vf(E_idx        )%sf(i,j,k) = 1d0
                        !q_prim_vf(adv_idx%beg  )%sf(i,j,k) = 1d0
                        ! ================================================================================

                     END IF

                  END DO
               END DO
            END DO
            
        END SUBROUTINE s_3D_analytical ! ---------------------------------------



        !>      This patch generates the shape of the spherical harmonics
        !!      as a perturbation to a perfect sphere
        !!      @param patch_id is the patch identifier
        SUBROUTINE s_spherical_harmonic(patch_id) ! ----------------------------

            INTEGER, INTENT(IN) :: patch_id 
            
            INTEGER :: i,j,k !< generic loop iterators

            COMPLEX(KIND(0d0)) :: cmplx_i = (0d0,1d0)
            COMPLEX(KIND(0d0)) :: H
            
            ! Transferring the patch's centroid and radius information
            x_centroid = patch_icpp(patch_id)%x_centroid
            y_centroid = patch_icpp(patch_id)%y_centroid
            z_centroid = patch_icpp(patch_id)%z_centroid
            radius     = patch_icpp(patch_id)%radius
            epsilon    = patch_icpp(patch_id)%epsilon
            beta       = patch_icpp(patch_id)%beta
            
            ! Since the analytical patch does not allow for its boundaries to get
            ! smoothed out, the pseudo volume fraction is set to 1 to make sure
            ! that only the current patch contributes to the fluid state in the
            ! cells that this patch covers.
            eta = 1d0
            
            ! Checking whether the patch covers a particular cell in the domain
            ! and verifying whether the current patch has permission to write to
            ! to that cell. If both queries check out, the primitive variables
            ! of the current patch are assigned to this cell.
            DO k = 0, p
               DO j = 0, n
                  DO i = 0, m

                     IF (grid_geometry == 3) THEN
                         CALL s_convert_cylindrical_to_cartesian_coord(y_cc(j),z_cc(k))
                     ELSE
                         cart_y = y_cc(j)
                         cart_z = z_cc(k)
                     END IF

                     IF((       (  x_cc(i) - x_centroid)**2d0                    &
                              + (   cart_y - y_centroid)**2d0                    &
                              + (   cart_z - z_centroid)**2d0 <= radius**2d0     &
                                                   .AND.                         &
                         patch_icpp(patch_id)%alter_patch(patch_id_fp(i,j,k))))  &
                                                   THEN

                        CALL s_convert_cylindrical_to_spherical_coord(x_cc(i),y_cc(j))

                        IF (epsilon  == 1d0) THEN
                            IF (beta == 0d0) THEN
                                H = 5d-1*SQRT(3d0/pi)*COS(sph_phi)
                            ELSEIF (beta == 1d0) THEN
                                H = -5d-1*SQRT(3d0/(2d0*pi))*EXP(cmplx_i*z_cc(k))*SIN(sph_phi)
                            END IF
                        ELSEIF (epsilon == 2d0) THEN
                            IF (beta == 0d0) THEN
                                H = 25d-2*SQRT(5d0/pi)*(3d0*COS(sph_phi)**2d0 - 1d0)
                            ELSEIF (beta == 1d0) THEN
                                H = -5d-1*SQRT(15d0/(2d0*pi))*EXP(cmplx_i*z_cc(k))*SIN(sph_phi)*COS(sph_phi)
                            ELSEIF (beta == 2d0) THEN
                                H = 25d-2*SQRT(15d0/(2d0*pi))*EXP(2d0*cmplx_i*z_cc(k))*SIN(sph_phi)**2d0
                            END IF
                        ELSEIF (epsilon == 3d0) THEN
                            IF (beta == 0d0) THEN
                                H = 25d-2*SQRT(7d0/pi)*(5d0*COS(sph_phi)**3d0 - 3d0*COS(sph_phi))
                            ELSEIF (beta == 1d0) THEN
                                H = -125d-3*SQRT(21d0/pi)*EXP(cmplx_i*z_cc(k))*SIN(sph_phi) * &
                                    (5d0*COS(sph_phi)**2d0 - 1d0)
                            ELSEIF (beta == 2d0) THEN
                                H = 25d-2*SQRT(105d0/(2d0*pi))*EXP(2d0*cmplx_i*z_cc(k)) * &
                                    SIN(sph_phi)**2d0*COS(sph_phi)
                            ELSEIF (beta == 3d0) THEN
                                H = -125d-3*SQRT(35d0/pi)*EXP(3d0*cmplx_i*z_cc(k))*SIN(sph_phi)**3d0
                            END IF
                        ELSEIF (epsilon == 4d0) THEN
                            IF (beta == 0d0) THEN
                                H = 3d0/16d0*SQRT(1d0/pi)*(35d0*COS(sph_phi)**4d0 - &
                                    3d1*COS(sph_phi)**2d0 + 3d0)
                            ELSEIF (beta == 1d0) THEN
                                H = -3d0/8d0*SQRT(5d0/pi)*EXP(cmplx_i*z_cc(k))* &
                                    SIN(sph_phi)*(7d0*COS(sph_phi)**3d0 - 3d0*COS(sph_phi))
                            ELSEIF (beta == 2d0) THEN
                                H = 3d0/8d0*SQRT(5d0/(2d0*pi))*EXP(2d0*cmplx_i*z_cc(k))* &
                                    SIN(sph_phi)**2d0*(7d0*COS(sph_phi)**2d0 - 1d0)
                            ELSEIF (beta == 3d0) THEN
                                H = -3d0/8d0*SQRT(35d0/pi)*EXP(3d0*cmplx_i*z_cc(k))* &
                                    SIN(sph_phi)**3d0*COS(sph_phi)
                            ELSEIF (beta == 4d0) THEN
                                H = 3d0/16d0*SQRT(35d0/(2d0*pi))*EXP(4d0*cmplx_i*z_cc(k))* &
                                    SIN(sph_phi)**4d0
                            END IF
                        ELSEIF (epsilon == 5d0) THEN
                            IF (beta == 0d0) THEN
                                H = 1d0/16d0*SQRT(11d0/pi)*(63d0*COS(sph_phi)**5d0 - &
                                    7d1*COS(sph_phi)**3d0 + 15d0*COS(sph_phi))
                            ELSEIF (beta == 1d0) THEN
                                H = -1d0/16d0*SQRT(165d0/(2d0*pi))*EXP(cmplx_i*z_cc(k))* &
                                    SIN(sph_phi)*(21d0*COS(sph_phi)**4d0 - 14d0*COS(sph_phi)**2d0 + 1d0)
                            ELSEIF (beta == 2d0) THEN
                                H = 125d-3*SQRT(1155d0/(2d0*pi))*EXP(2d0*cmplx_i*z_cc(k)) * &
                                    SIN(sph_phi)**2d0*(3d0*COS(sph_phi)**3d0 - COS(sph_phi))
                            ELSEIF (beta == 3d0) THEN
                                H = -1d0/32d0*SQRT(385d0/pi)*EXP(3d0*cmplx_i*z_cc(k) )* &
                                    SIN(sph_phi)**3d0*(9d0*COS(sph_phi)**2d0 - 1d0)
                            ELSEIF (beta == 4d0) THEN
                                H = 3d0/16d0*SQRT(385d0/(2d0*pi))*EXP(4d0*cmplx_i*z_cc(k)) * &
                                    SIN(sph_phi)**4d0*COS(sph_phi)
                            ELSEIF (beta == 5d0) THEN
                                H = -3d0/32d0*SQRT(77d0/pi)*EXP(5d0*cmplx_i*z_cc(k)) * &
                                    SIN(sph_phi)**5d0
                            END IF
                        END IF

                        q_prim_vf(adv_idx%beg)%sf(i,j,k) = 1d0-ABS(REAL(H,KIND(0d0)))

                     END IF

                  END DO
               END DO
            END DO
            
        END SUBROUTINE s_spherical_harmonic ! ----------------------------------


        !>          The spherical patch is a 3D geometry that may be used,
        !!              for example, in creating a bubble or a droplet. The patch
        !!              geometry is well-defined when its centroid and radius are
        !!              provided. Please note that the spherical patch DOES allow
        !!              for the smoothing of its boundary.
        !!      @param patch_id is the patch identifier
        SUBROUTINE s_sphere(patch_id) ! ----------------------------------------

            INTEGER, INTENT(IN) :: patch_id
            
            ! Generic loop iterators
            INTEGER :: i,j,k !< generic loop iterators


            REAL(KIND(0d0)) :: radius_pressure, pressure_bubble, pressure_inf !< 
            !! Variables to initialize the pressure field that corresponds to the 
            !! bubble-collapse test case found in Tiwari et al. (2013)
            
            
            ! Transferring spherical patch's radius, centroid, smoothing patch
            ! identity and smoothing coefficient information
            x_centroid      = patch_icpp(patch_id)%x_centroid
            y_centroid      = patch_icpp(patch_id)%y_centroid
            z_centroid      = patch_icpp(patch_id)%z_centroid
            radius          = patch_icpp(patch_id)%radius
            smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
            smooth_coeff    = patch_icpp(patch_id)%smooth_coeff
            
            
            ! Initializing the pseudo volume fraction value to 1. The value will
            ! be modified as the patch is laid out on the grid, but only in the
            ! case that smoothing of the spherical patch's boundary is enabled.
            eta = 1d0
            
            
            ! Checking whether the sphere covers a particular cell in the domain
            ! and verifying whether the current patch has permission to write to
            ! that cell. If both queries check out, the primitive variables of
            ! the current patch are assigned to this cell.
            DO k = 0, p
               DO j = 0, n
                  DO i = 0, m
                     
                     IF (grid_geometry == 3) THEN
                         CALL s_convert_cylindrical_to_cartesian_coord(y_cc(j),z_cc(k))
                     ELSE
                         cart_y = y_cc(j)
                         cart_z = z_cc(k)
                     END IF

                     IF(patch_icpp(patch_id)%smoothen) THEN
                        
                        eta = TANH( smooth_coeff / MIN(dx,dy,dz) *  &
                            ( SQRT( (  x_cc(i) - x_centroid)**2d0   &
                                  + (   cart_y - y_centroid)**2d0   &
                                  + (   cart_z - z_centroid)**2d0 ) &
                                  - radius ) ) * (-0.5d0) + 0.5d0
                        
                     END IF
                     
                     IF((       (  x_cc(i) - x_centroid)**2d0                    &
                              + (   cart_y - y_centroid)**2d0                    &
                              + (   cart_z - z_centroid)**2d0 <= radius**2d0     &
                                                   .AND.                         &
                         patch_icpp(patch_id)%alter_patch(patch_id_fp(i,j,k)))   &
                                                   .OR.                          &
                                 patch_id_fp(i,j,k)== smooth_patch_id        )   &
                                                   THEN
                        
                        CALL s_assign_patch_primitive_variables(patch_id, i,j,k)
                        
                     END IF

                     ! Initialization of the pressure field that corresponds to the bubble-collapse 
                     !! test case found in Tiwari et al. (2013)
                     ! radius_pressure = SQRT(x_cc(i)**2d0) ! 1D
                     ! radius_pressure = SQRT(x_cc(i)**2d0 + cart_y**2d0) ! 2D
                     ! radius_pressure = SQRT(x_cc(i)**2d0 + cart_y**2d0 + cart_z**2d0) ! 3D
                     ! pressure_bubble = 1.E+04
                     ! pressure_inf    = 1.E+05
                     ! q_prim_vf(E_idx)%sf(i,j,k) = pressure_inf + radius / radius_pressure * (pressure_bubble - pressure_inf)
                     ! 
                     ! IF(       ((  x_cc(i) - x_centroid)**2d0                    &
                     !          + (   cart_y - y_centroid)**2d0                    &
                     !          + (   cart_z - z_centroid)**2d0) <= radius**2d0)   &
                     !                               THEN
                     !    
                     !    q_prim_vf(E_idx)%sf(i,j,k) = pressure_bubble
                     !    
                     ! END IF

                  END DO
               END DO
            END DO
            
            
        END SUBROUTINE s_sphere ! ----------------------------------------------
        
        
        !>      The cuboidal patch is a 3D geometry that may be used, for
        !!              example, in creating a solid boundary, or pre-/post-shock
        !!              region, which is aligned with the axes of the Cartesian
        !!              coordinate system. The geometry of such a patch is well-
        !!              defined when its centroid and lengths in the x-, y- and
        !!              z-coordinate directions are provided. Please notice that
        !!              the cuboidal patch DOES NOT allow for the smearing of its
        !!              boundaries.
        !!      @param patch_id is the patch identifier
        SUBROUTINE s_cuboid(patch_id) ! ----------------------------------------
            
            INTEGER, INTENT(IN) :: patch_id
            
            INTEGER :: i,j,k !< Generic loop iterators
            
            
            ! Transferring the cuboid's centroid and length information
            x_centroid = patch_icpp(patch_id)%x_centroid
            y_centroid = patch_icpp(patch_id)%y_centroid
            z_centroid = patch_icpp(patch_id)%z_centroid
            length_x   = patch_icpp(patch_id)%length_x
            length_y   = patch_icpp(patch_id)%length_y
            length_z   = patch_icpp(patch_id)%length_z
            
            
            ! Computing the beginning and the end x-, y- and z-coordinates of
            ! the cuboid based on its centroid and lengths
            x_boundary%beg = x_centroid - 0.5d0*length_x
            x_boundary%end = x_centroid + 0.5d0*length_x
            y_boundary%beg = y_centroid - 0.5d0*length_y
            y_boundary%end = y_centroid + 0.5d0*length_y
            z_boundary%beg = z_centroid - 0.5d0*length_z
            z_boundary%end = z_centroid + 0.5d0*length_z
            
            
            ! Since the cuboidal patch does not allow for its boundaries to get
            ! smoothed out, the pseudo volume fraction is set to 1 to make sure
            ! that only the current patch contributes to the fluid state in the
            ! cells that this patch covers.
            eta = 1d0
            
            
            ! Checking whether the cuboid covers a particular cell in the domain
            ! and verifying whether the current patch has permission to write to
            ! to that cell. If both queries check out, the primitive variables
            ! of the current patch are assigned to this cell.
            DO k = 0, p
               DO j = 0, n
                  DO i = 0, m

                     IF (grid_geometry == 3) THEN
                         CALL s_convert_cylindrical_to_cartesian_coord(y_cc(j),z_cc(k))
                     ELSE
                         cart_y = y_cc(j)
                         cart_z = z_cc(k)
                     END IF

                     IF(            x_boundary%beg <= x_cc(i)   .AND.            &
                                    x_boundary%end >= x_cc(i)   .AND.            &
                                    y_boundary%beg <= cart_y    .AND.            &
                                    y_boundary%end >= cart_y    .AND.            &
                                    z_boundary%beg <= cart_z    .AND.            &
                                    z_boundary%end >= cart_z                     &
                                                .AND.                            &
                         patch_icpp(patch_id)%alter_patch(patch_id_fp(i,j,k)))   &
                                                THEN
                        
                        CALL s_assign_patch_primitive_variables(patch_id, i,j,k)
                        
                     END IF
                  END DO
               END DO
            END DO
            
            
        END SUBROUTINE s_cuboid ! ----------------------------------------------
        
        
        !>              The cylindrical patch is a 3D geometry that may be used,
        !!              for example, in setting up a cylindrical solid boundary
        !!              confinement, like a blood vessel. The geometry of this
        !!              patch is well-defined when the centroid, the radius and
        !!              the length along the cylinder's axis, parallel to the x-,
        !!              y- or z-coordinate direction, are provided. Please note
        !!              that the cylindrical patch DOES allow for the smoothing
        !!              of its lateral boundary.       
        !!      @param patch_id is the patch identifier
        SUBROUTINE s_cylinder(patch_id) ! --------------------------------------

            INTEGER, INTENT(IN) :: patch_id
            
            INTEGER :: i,j,k !< Generic loop iterators
            
            
            ! Transferring the cylindrical patch's centroid, length, radius,
            ! smoothing patch identity and smoothing coefficient information
            x_centroid      = patch_icpp(patch_id)%x_centroid
            y_centroid      = patch_icpp(patch_id)%y_centroid
            z_centroid      = patch_icpp(patch_id)%z_centroid
            length_x        = patch_icpp(patch_id)%length_x
            length_y        = patch_icpp(patch_id)%length_y
            length_z        = patch_icpp(patch_id)%length_z
            radius          = patch_icpp(patch_id)%radius
            smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
            smooth_coeff    = patch_icpp(patch_id)%smooth_coeff
            
            
            ! Computing the beginning and the end x-, y- and z-coordinates of
            ! the cylinder based on its centroid and lengths
            x_boundary%beg = x_centroid - 0.5d0*length_x
            x_boundary%end = x_centroid + 0.5d0*length_x
            y_boundary%beg = y_centroid - 0.5d0*length_y
            y_boundary%end = y_centroid + 0.5d0*length_y
            z_boundary%beg = z_centroid - 0.5d0*length_z
            z_boundary%end = z_centroid + 0.5d0*length_z
            
            
            ! Initializing the pseudo volume fraction value to 1. The value will
            ! be modified as the patch is laid out on the grid, but only in the
            ! case that smearing of the cylindrical patch's boundary is enabled.
            eta = 1d0
            
            
            ! Checking whether the cylinder covers a particular cell in the
            ! domain and verifying whether the current patch has the permission
            ! to write to that cell. If both queries check out, the primitive
            ! variables of the current patch are assigned to this cell.
            DO k = 0, p
               DO j = 0, n
                  DO i = 0, m
                     
                     IF (grid_geometry == 3) THEN
                         CALL s_convert_cylindrical_to_cartesian_coord(y_cc(j),z_cc(k))
                     ELSE
                         cart_y = y_cc(j)
                         cart_z = z_cc(k)
                     END IF

                     IF(patch_icpp(patch_id)%smoothen) THEN
                         
                         IF(length_x /= dflt_real) THEN
                            eta = TANH( smooth_coeff / MIN(dy,dz) *     &
                                ( SQRT( (   cart_y - y_centroid)**2d0   &
                                      + (   cart_z - z_centroid)**2d0 ) &
                                      - radius ) ) * (-0.5d0) + 0.5d0
                        ELSEIF(length_y /= dflt_real) THEN
                            eta = TANH( smooth_coeff / MIN(dx,dz) *   &
                                ( SQRT( (x_cc(i) - x_centroid)**2d0   &
                                    + (   cart_z - z_centroid)**2d0 ) &
                                    - radius ) ) * (-0.5d0) + 0.5d0
                        ELSE
                            eta = TANH( smooth_coeff / MIN(dx,dy) *   &
                                ( SQRT( (x_cc(i) - x_centroid)**2d0   &
                                    + (   cart_y - y_centroid)**2d0 ) &
                                    - radius ) ) * (-0.5d0) + 0.5d0
                        END IF
                        
                     END IF
                       
                     IF((((                length_x /= dflt_real         .AND. &
                            (   cart_y - y_centroid)**2d0                      &
                          + (   cart_z - z_centroid)**2d0 <= radius**2d0 .AND. &
                                         x_boundary%beg <= x_cc(i)       .AND. &
                                         x_boundary%end >= x_cc(i)          )  &
                                                   .OR.                        &
                          (                length_y /= dflt_real         .AND. &
                              (x_cc(i) - x_centroid)**2d0                      &
                          + (cart_z - z_centroid)**2d0 <= radius**2d0    .AND. &
                                         y_boundary%beg <= cart_y        .AND. &
                                         y_boundary%end >= cart_y           )  &
                                                   .OR.                        &
                          (                length_z /= dflt_real         .AND. &
                              (x_cc(i) - x_centroid)**2d0                      &
                          + (   cart_y - y_centroid)**2d0 <= radius**2d0 .AND. &
                                         z_boundary%beg <= cart_z        .AND. &
                                         z_boundary%end >= cart_z          ))  &
                                                   .AND.                       &
                          patch_icpp(patch_id)%alter_patch(patch_id_fp(i,j,k)))&
                                                   .OR.                        &
                                 patch_id_fp(i,j,k)== smooth_patch_id         )&
                                                   THEN
                         
                        CALL s_assign_patch_primitive_variables(patch_id, i,j,k)
                        
                     END IF
                     
                  END DO
               END DO
            END DO
            
            
        END SUBROUTINE s_cylinder ! --------------------------------------------
        
        
        !>      The swept plane patch is a 3D geometry that may be used,
        !!              for example, in creating a solid boundary, or pre-/post-
        !!              shock region, at an angle with respect to the axes of the
        !!              Cartesian coordinate system. The geometry of the patch is
        !!              well-defined when its centroid and normal vector, aimed
        !!              in the sweep direction, are provided. Note that the sweep
        !!              plane patch DOES allow the smoothing of its boundary.       
        !!      @param patch_id is the patch identifier
        SUBROUTINE s_sweep_plane(patch_id) ! -----------------------------------

            INTEGER, INTENT(IN) :: patch_id
            
            INTEGER :: i,j,k !< Generic loop iterators
            
            ! Transferring the centroid information of the plane to be swept
            x_centroid = patch_icpp(patch_id)%x_centroid
            y_centroid = patch_icpp(patch_id)%y_centroid
            z_centroid = patch_icpp(patch_id)%z_centroid
            smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
            smooth_coeff = patch_icpp(patch_id)%smooth_coeff
            
            
            ! Obtaining coefficients of the equation describing the sweep plane
            a = patch_icpp(patch_id)%normal(1)
            b = patch_icpp(patch_id)%normal(2)
            c = patch_icpp(patch_id)%normal(3)
            d = -a*x_centroid - b*y_centroid - c*z_centroid
            
            
            ! Initializing the pseudo volume fraction value to 1. The value will
            ! be modified as the patch is laid out on the grid, but only in the
            ! case that smearing of the sweep plane patch's boundary is enabled.
            eta = 1d0
            
            
            ! Checking whether the region swept by the plane covers a particular
            ! cell in the domain and verifying whether the current patch has the
            ! permission to write to that cell. If both queries check out, the
            ! primitive variables of the current patch are written to this cell.
            DO k = 0, p
               DO j = 0, n
                  DO i = 0, m
                     
                     IF (grid_geometry == 3) THEN
                         CALL s_convert_cylindrical_to_cartesian_coord(y_cc(j),z_cc(k))
                     ELSE
                         cart_y = y_cc(j)
                         cart_z = z_cc(k)
                     END IF

                     IF(patch_icpp(patch_id)%smoothen) THEN
                        eta = 5d-1 + 5d-1*TANH( smooth_coeff/MIN(dx,dy,dz) &
                                              *    ( a*x_cc(i) +           &
                                                     b*cart_y +            &
                                                     c*cart_z + d )        &
                                              / SQRT(a**2d0+b**2d0+c**2d0) )
                     END IF
                     
                     IF((     a*x_cc(i) + b*cart_y + c*cart_z + d >= 0d0       &
                                                  .AND.                        &
                         patch_icpp(patch_id)%alter_patch(patch_id_fp(i,j,k))) &
                                                   .OR.                        &
                                 patch_id_fp(i,j,k)== smooth_patch_id        ) &
                                                   THEN
                        
                        CALL s_assign_patch_primitive_variables(patch_id, i,j,k)
                        
                     END IF
                     
                  END DO
               END DO
            END DO
            
            
        END SUBROUTINE s_sweep_plane ! -----------------------------------------
        
        
        
        
        !>  Deallocation procedures for the module
        SUBROUTINE s_finalize_initial_condition_module() ! ---------------------
            
            INTEGER :: i !< Generic loop iterator
            
            
            ! Dellocating the primitive and conservative variables
            DO i = 1, sys_size
                DEALLOCATE(q_prim_vf(i)%sf)
                DEALLOCATE(q_cons_vf(i)%sf)
            END DO
            
            DEALLOCATE(q_prim_vf)
            DEALLOCATE(q_cons_vf)
            
            
            ! Deallocating the patch identities bookkeeping variable
            DEALLOCATE(patch_id_fp)
            
            
            ! Nullifying procedure pointer to the subroutine assigning either
            ! the patch mixture or species primitive variables to a cell in the
            ! computational domain
            s_assign_patch_primitive_variables => NULL()
            
            
        END SUBROUTINE s_finalize_initial_condition_module ! -------------------
        
        
END MODULE m_initial_condition
