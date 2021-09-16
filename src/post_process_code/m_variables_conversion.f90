!>
!! @file m_variables_conversion.f90
!! @brief Contains module m_variables_conversion
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief This module consists of subroutines used in the conversion of the
!!              conservative variables into the primitive variables. In addition,
!!              the module also contains subroutines used to compute the mixture
!!              variables. For a specific time-step undergoing the post-process,
!!              the mixture variables are stored everywhere on the grid so that
!!              they may possibly be outputted or used to compute other derived
!!              flow quantities.
MODULE m_variables_conversion
    
    
    ! Dependencies =============================================================
    USE m_derived_types         !< Definitions of the derived types
    
    USE m_global_parameters     !< Global parameters for the code
    ! ==========================================================================
    
    
    IMPLICIT NONE
    
    PRIVATE; PUBLIC :: s_initialize_variables_conversion_module, &
                       s_convert_to_mixture_variables, &
                       s_convert_mixture_to_mixture_variables, &
                       s_convert_species_to_mixture_variables_bubbles, &
                       s_convert_species_to_mixture_variables, &
                       s_convert_conservative_to_primitive_variables, &
                       s_finalize_variables_conversion_module
    
    !> Abstract interface to two subroutines designed for the transfer/conversion
    !! of the mixture/species variables to the mixture variables
    ABSTRACT INTERFACE
        
        !>  Structure of the s_convert_mixture_to_mixture_variables
        !!      and s_convert_species_to_mixture_variables subroutines
        !! @param q_cons_vf Conservative variables
        !! @param i cell index to transfer mixture variables 
        !! @param j cell index to transfer mixture variables
        !! @param k cell index to transfer mixture variables
        SUBROUTINE s_convert_xxxxx_to_mixture_variables(q_cons_vf, i,j,k)
            
            ! Importing the derived type scalar_field from m_derived_types.f90
            ! and global variable sys_size, from m_global_variables.f90, as
            ! the abstract interface does not inherently have access to them
            IMPORT :: scalar_field, sys_size
            
            TYPE(scalar_field), &
            DIMENSION(sys_size), &
            INTENT(IN) :: q_cons_vf
            
            INTEGER, INTENT(IN) :: i,j,k
            
        END SUBROUTINE s_convert_xxxxx_to_mixture_variables
        
    END INTERFACE
    
    ! NOTE: These abstract interfaces allow for the declaration of a pointer to
    ! a procedure such that the choice of the model equations does not have to
    ! be queried every time the mixture quantities are needed.
    
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:,:), PUBLIC :: rho_sf !< Scalar density function
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:,:), PUBLIC :: gamma_sf !< Scalar sp. heat ratio function
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:,:), PUBLIC :: pi_inf_sf !< Scalar liquid stiffness function
    
    PROCEDURE(s_convert_xxxxx_to_mixture_variables), &
    POINTER :: s_convert_to_mixture_variables => NULL() !<
    !! Pointer referencing the subroutine s_convert_mixture_to_mixture_variables
    !! or s_convert_species_to_mixture_variables, based on model equations choice
    
    INTEGER, PRIVATE :: flg !<
    !! Flagging (flg) variable used to annotate the dimensionality of the dataset
    !! that is undergoing the post-process. A flag value of 1 indicates that the
    !! dataset is 3D, while a flag value of 0 indicates that it is not. This flg
    !! variable is necessary to avoid cycling through the third dimension of the
    !! flow variable(s) when the simulation is not 3D and the size of the buffer
    !! is non-zero. Note that a similar procedure does not have to be applied to
    !! the second dimension since in 1D, the buffer size is always zero.
    
    
    CONTAINS
        
        
        
        
        !>  This subroutine is constructed for the gamma/pi_inf model
        !!      and provided a set of conservative variables, transfers
        !!      the density, specific heat ratio function and the liquid
        !!      stiffness function from q_cons_vf to rho_sf, gamma_sf and
        !!      pi_inf_sf.
        !! @param q_cons_vf Conservative variables
        !! @param i cell index to transfer mixture variables 
        !! @param j cell index to transfer mixture variables
        !! @param k cell index to transfer mixture variables
        SUBROUTINE s_convert_mixture_to_mixture_variables(q_cons_vf, i,j,k) ! --
            
            TYPE(scalar_field), &
            DIMENSION(sys_size), &
            INTENT(IN) :: q_cons_vf
            
            INTEGER, INTENT(IN) :: i,j,k
            
            ! Transfering the density, the specific heat ratio function and the
            ! liquid stiffness function, respectively
            rho_sf(i,j,k)    = q_cons_vf(1)%sf(i,j,k)
            gamma_sf(i,j,k)  = q_cons_vf(gamma_idx)%sf(i,j,k)
            pi_inf_sf(i,j,k) = q_cons_vf(pi_inf_idx)%sf(i,j,k)
            
        END SUBROUTINE s_convert_mixture_to_mixture_variables ! ----------------
        
        
        !>  This subroutine is designed for the volume fraction model
        !!      with sub-grid ensemble averaged bubbles (Bryngelson 2019)
        !!      and provided a set of conservative variables, calculates
        !!      the density, the specific heat ratio function and liquid
        !!      stiffness function from q_cons_vf and stores the results
        !!      into rho_sf, gamma_sf and pi_inf_sf.
        !!  @param qK_vf Conservative variables
        !!  @param j cell index to transfer mixture variables 
        !!  @param k cell index to transfer mixture variables
        !!  @param l cell index to transfer mixture variables
        SUBROUTINE s_convert_species_to_mixture_variables_bubbles (qK_vf,j,k,l)
                                                            
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: qK_vf 
            INTEGER, INTENT(IN) :: j,k,l
            
            INTEGER :: i !< Generic loop iterator

            IF (model_eqns == 4) THEN
                rho_sf    = qK_vf(1)%sf(j,k,l)
                gamma_sf  = fluid_pp(1)%gamma   
                pi_inf_sf = fluid_pp(1)%pi_inf      
            ELSE IF (model_eqns == 2 .AND. bubbles .AND. adv_alphan) THEN
                rho_sf = 0d0; gamma_sf = 0d0; pi_inf_sf = 0d0

                IF (mpp_lim .AND. num_fluids > 1) THEN
                    DO i = 1, num_fluids
                        rho_sf    = rho_sf    + qK_vf(i+E_idx)%sf(j,k,l)*qK_vf(i)%sf(j,k,l) 
                        gamma_sf  = gamma_sf  + qK_vf(i+E_idx)%sf(j,k,l)*fluid_pp(i)%gamma
                        pi_inf_sf = pi_inf_sf + qK_vf(i+E_idx)%sf(j,k,l)*fluid_pp(i)%pi_inf
                    END DO
                ELSE IF (num_fluids > 1) THEN
                    DO i = 1, num_fluids-1 !leave out bubble part of mixture
                        rho_sf    = rho_sf    + qK_vf(i+E_idx)%sf(j,k,l)*qK_vf(i)%sf(j,k,l) 
                        gamma_sf  = gamma_sf  + qK_vf(i+E_idx)%sf(j,k,l)*fluid_pp(i)%gamma
                        pi_inf_sf = pi_inf_sf + qK_vf(i+E_idx)%sf(j,k,l)*fluid_pp(i)%pi_inf
                    END DO
                ELSE
                    rho_sf    = qK_vf(1)%sf(j,k,l)
                    gamma_sf  = fluid_pp(1)%gamma
                    pi_inf_sf = fluid_pp(1)%pi_inf   
                END IF
            END IF

        END SUBROUTINE s_convert_species_to_mixture_variables_bubbles ! ----------------        
        

        !>  This subroutine is designed for the volume fraction model
        !!      and provided a set of conservative variables, calculates
        !!      the density, the specific heat ratio function and liquid
        !!      stiffness function from q_cons_vf and stores the results
        !!      into rho_sf, gamma_sf and pi_inf_sf.
        !!  @param q_cons_vf Conservative variables
        !!  @param j cell index to transfer mixture variables 
        !!  @param k cell index to transfer mixture variables
        !!  @param l cell index to transfer mixture variables
        SUBROUTINE s_convert_species_to_mixture_variables(q_cons_vf, j,k,l) ! --

            TYPE(scalar_field), &
            DIMENSION(sys_size), &
            INTENT(IN) :: q_cons_vf
            
            INTEGER, INTENT(IN) :: j,k,l
            
            INTEGER :: i !< Generic loop iterator
            
            
            ! Computing the density, the specific heat ratio function and the
            ! liquid stiffness function, respectively
            IF(adv_alphan) THEN
                IF (bubbles .NEQV. .True.) THEN
                        rho_sf(j,k,l)    = 0d0
                        gamma_sf(j,k,l)  = 0d0
                        pi_inf_sf(j,k,l) = 0d0
                
                        DO i = 1, num_fluids
                                rho_sf(j,k,l)    = rho_sf(j,k,l)                &
                                     + q_cons_vf(   i   )%sf(j,k,l)
                                gamma_sf(j,k,l)  = gamma_sf(j,k,l)              &
                                     + q_cons_vf(i+E_idx)%sf(j,k,l) &
                                     * fluid_pp(i)%gamma
                                pi_inf_sf(j,k,l) = pi_inf_sf(j,k,l)             &
                                     + q_cons_vf(i+E_idx)%sf(j,k,l) &
                                     * fluid_pp(i)%pi_inf
                        END DO
                ELSE
                        rho_sf(j,k,l)    = q_cons_vf(1)%sf(j,k,l)
                        gamma_sf(j,k,l)  = fluid_pp(1)%gamma
                        pi_inf_sf(j,k,l) = fluid_pp(1)%pi_inf
                END IF
            ELSE
                rho_sf(j,k,l)    = q_cons_vf(num_fluids)%sf(j,k,l)
                gamma_sf(j,k,l)  = fluid_pp(num_fluids)%gamma
                pi_inf_sf(j,k,l) = fluid_pp(num_fluids)%pi_inf
                
                DO i = 1, num_fluids-1
                    rho_sf(j,k,l)    = rho_sf(j,k,l)                 &
                                     + q_cons_vf(   i   )%sf(j,k,l)
                    gamma_sf(j,k,l)  = gamma_sf(j,k,l)               &
                                     + q_cons_vf(i+E_idx)%sf(j,k,l)  &
                                     * ( fluid_pp(     i    )%gamma  &
                                       - fluid_pp(num_fluids)%gamma  )
                    pi_inf_sf(j,k,l) = pi_inf_sf(j,k,l)              &
                                     + q_cons_vf(i+E_idx)%sf(j,k,l)  &
                                     * ( fluid_pp(     i    )%pi_inf &
                                       - fluid_pp(num_fluids)%pi_inf )
                END DO
                
            END IF
            
            
        END SUBROUTINE s_convert_species_to_mixture_variables ! ----------------
        
 

        !>  Computation of parameters, allocation procedures, and/or
        !!      any other tasks needed to properly setup the module
        SUBROUTINE s_initialize_variables_conversion_module() ! -------------------

            
            ! Allocating the density, the specific heat ratio function and the
            ! liquid stiffness function, respectively
            
            ! Simulation is at least 2D
            IF(n > 0) THEN
                
                ! Simulation is 3D
                IF(p > 0) THEN
                    
                    ALLOCATE(   rho_sf( -buff_size:m+buff_size, &
                                        -buff_size:n+buff_size, &
                                        -buff_size:p+buff_size ))
                    ALLOCATE( gamma_sf( -buff_size:m+buff_size, &
                                        -buff_size:n+buff_size, &
                                        -buff_size:p+buff_size ))
                    ALLOCATE(pi_inf_sf( -buff_size:m+buff_size, &
                                        -buff_size:n+buff_size, &
                                        -buff_size:p+buff_size ))
                    
                ! Simulation is 2D
                ELSE
                    
                    ALLOCATE(   rho_sf( -buff_size:m+buff_size, &
                                        -buff_size:n+buff_size, &
                                                 0:0           ))
                    ALLOCATE( gamma_sf( -buff_size:m+buff_size, &
                                        -buff_size:n+buff_size, &
                                                 0:0           ))
                    ALLOCATE(pi_inf_sf( -buff_size:m+buff_size, &
                                        -buff_size:n+buff_size, &
                                                 0:0           ))
                    
                END IF
                
            ! Simulation is 1D
            ELSE
                
                ALLOCATE(   rho_sf( -buff_size:m+buff_size, &
                                             0:0          , &
                                             0:0           ))
                ALLOCATE( gamma_sf( -buff_size:m+buff_size, &
                                             0:0          , &
                                             0:0           ))
                ALLOCATE(pi_inf_sf( -buff_size:m+buff_size, &
                                             0:0          , &
                                             0:0           ))
                
            END IF
            
            
            ! Depending on the model selection for the equations of motion, the
            ! appropriate procedure for the conversion to the mixture variables
            ! is targeted by the procedure pointer
            
            IF(model_eqns == 1) THEN        ! Gamma/pi_inf model
                s_convert_to_mixture_variables => &
                            s_convert_mixture_to_mixture_variables
            ELSEIF (bubbles) THEN           ! Volume fraction model
                s_convert_to_mixture_variables => &
                            s_convert_species_to_mixture_variables_bubbles
            ELSE                            ! Volume fraction model
                s_convert_to_mixture_variables => &
                            s_convert_species_to_mixture_variables
            END IF
            
            
            ! Annotating the dimensionality of the dataset undergoing the post-
            ! process. A flag value of 1 indicates that the dataset is 3D, while
            ! a flag value of 0 indicates that it is not.
            IF(p > 0) THEN
                flg = 1
            ELSE
                flg = 0
            END IF
            
            
        END SUBROUTINE s_initialize_variables_conversion_module ! -----------------
        
        
        
        
        !> Converts the conservative variables to the primitive ones      
        !!  @param q_cons_vf Conservative variabels
        !!  @param q_prim_vf Primitive variables
        SUBROUTINE s_convert_conservative_to_primitive_variables( q_cons_vf, &
                                                                  q_prim_vf  )

            TYPE(scalar_field), &
            DIMENSION(sys_size), &
            INTENT(IN) :: q_cons_vf
            
            TYPE(scalar_field), &
            DIMENSION(sys_size), &
            INTENT(INOUT) :: q_prim_vf
            
            ! Dynamic pressure, as defined in the incompressible flow sense
            REAL(KIND(0d0)) :: dyn_pres

            ! Bubble parameters
            REAL(KIND(0d0)) :: nbub 
            REAL(KIND(0d0)), dimension(:), allocatable :: nRtmp 
           
            INTEGER :: i,j,k,l !< Generic loop iterators
            
            allocate( nRtmp(nb) )
            
            ! Converting the conservative variables to the primitive variables
            DO l = -buff_size*flg, (p+buff_size)*flg
                DO k = -buff_size, n+buff_size
                    DO j = -buff_size, m+buff_size
                        
                        ! Obtaining the density, specific heat ratio function
                        ! and the liquid stiffness function, respectively
                        CALL s_convert_to_mixture_variables(q_cons_vf, j,k,l)
 
                        ! Transferring the continuity equation(s) variable(s)
                        DO i = 1, cont_idx%end
                            q_prim_vf(i)%sf(j,k,l) = q_cons_vf(i)%sf(j,k,l)
                        END DO
                        
                        ! Zeroing out the dynamic pressure since it is computed
                        ! iteratively by cycling through the momentum equations
                        dyn_pres = 0d0
                        
                        ! Computing velocity and dynamic pressure from momenta
                        DO i = mom_idx%beg, mom_idx%end
                            IF (model_eqns==4) THEN
                                !u = \rho u/ \rho
                                q_prim_vf(i)%sf(j,k,l) = q_cons_vf(i)%sf(j,k,l) / &
                                                     q_cons_vf(1)%sf(j,k,l)
                            ELSE
                                q_prim_vf(i)%sf(j,k,l) = q_cons_vf(i)%sf(j,k,l) / &
                                                     rho_sf(j,k,l)
                                dyn_pres = dyn_pres + q_cons_vf(i)%sf(j,k,l) * &
                                                  q_prim_vf(i)%sf(j,k,l) / 2d0
                            END IF 
                        END DO
                        
                        IF (model_eqns == 4) THEN
                            ! Computing the pressure from the energy
                            ! Tait EOS
                            ! p_l = (pl0 + pi_infty)(rho/(rho_l0(1-alf)))^gamma
                            ! - pi_infty  
                            q_prim_vf(E_idx)%sf(j,k,l) = & 
                                   (pref+fluid_pp(1)%pi_inf) * &
                                   (( & 
                                   q_cons_vf(1)%sf(j,k,l)/ &
                                   (rhoref*(1.d0-q_cons_vf(alf_idx)%sf(j,k,l))) & 
                                   ) ** (1.d0/fluid_pp(1)%gamma + 1.d0)) - fluid_pp(1)%pi_inf
                        ELSE IF ( (model_eqns .NE. 4) .AND. bubbles .NEQV. .TRUE. ) THEN
                            ! Computing the pressure from the energy
                            q_prim_vf(E_idx)%sf(j,k,l) = &
                                ( q_cons_vf(E_idx)%sf(j,k,l) &
                                - dyn_pres -pi_inf_sf(j,k,l) ) / gamma_sf(j,k,l)
                        ELSE
                            ! p = ( E/(1-alf) - 0.5 rho u u/(1-alf) - pi_inf_k )/gamma_k
                            q_prim_vf(E_idx)%sf(j,k,l) = &
                                ( (q_cons_vf(E_idx)%sf(j,k,l) &
                                                - dyn_pres) / (1.d0 - q_cons_vf(alf_idx)%sf(j,k,l)) &
                                                - pi_inf_sf(j,k,l) &
                                                ) / gamma_sf(j,k,l)
                        END IF
                        
                        ! Set partial pressures to mixture pressure
                        IF(model_eqns == 3) THEN
                            DO i = internalEnergies_idx%beg, internalEnergies_idx%end
                                q_prim_vf(i)%sf(j,k,l) = q_prim_vf(E_idx)%sf(j,k,l)
                            END DO
                        END IF
                        
                        ! Transferring the advection equation(s) variable(s)
                        DO i = adv_idx%beg, adv_idx%end
                            q_prim_vf(i)%sf(j,k,l) = q_cons_vf(i)%sf(j,k,l)
                        END DO
                        
                        ! \phi = (n\phi)/n  (n = nbub)
                        IF (bubbles) THEN
                            ! n = sqrt( 4pi/(3 alpha) * (nR)**3 )
                            DO i = 1,nb
                                nRtmp(i) = q_cons_vf(bub_idx%rs(i))%sf(j,k,l)
                            END DO
                            CALL s_comp_n_from_cons( q_cons_vf(alf_idx)%sf(j,k,l), nRtmp, nbub)
                            DO i = bub_idx%beg, sys_size
                                q_prim_vf(i)%sf(j,k,l) = q_cons_vf(i)%sf(j,k,l)/nbub
                            END DO
                        END IF
                    END DO
                END DO
            END DO
            
           DEALLOCATE( nRtmp )
 
        END SUBROUTINE s_convert_conservative_to_primitive_variables ! ---------
        
        !> Deallocation procedures for the module
        SUBROUTINE s_finalize_variables_conversion_module() ! ----------------
            
            ! Deallocating the density, the specific heat ratio function and the
            ! liquid stiffness function
            DEALLOCATE(rho_sf)
            DEALLOCATE(gamma_sf)
            DEALLOCATE(pi_inf_sf)
            
            
            ! Nullifying the procedure pointer to the subroutine transfering/
            ! computing the mixture/species variables to the mixture variables
            s_convert_to_mixture_variables => NULL()
            
        END SUBROUTINE s_finalize_variables_conversion_module ! --------------
        
        
END MODULE m_variables_conversion
