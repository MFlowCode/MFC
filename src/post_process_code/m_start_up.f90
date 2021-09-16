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
!! @file m_start_up.f90
!! @brief  Contains module m_start_up
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief This module contains the subroutines that read in and check the
!!              consistency of the user provided inputs.
MODULE m_start_up
    
    
    ! Dependencies =============================================================
    USE m_derived_types        !< Definitions of the derived types
    
    USE m_global_parameters    !< Definitions of the global parameters
    
    USE m_mpi_proxy            !< Message passing interface (MPI) module proxy

    USE m_compile_specific
    ! ==========================================================================
    
    
    IMPLICIT NONE
    
    
    CONTAINS
        

        !>  Reads the configuration file post_process.inp, in order
        !!      to populate parameters in module m_global_parameters.f90
        !!      with the user provided inputs
        SUBROUTINE s_read_input_file() ! ---------------------------------------


            CHARACTER(LEN = name_len) :: file_loc !<
            !! Generic string used to store the address of a particular file
            
            LOGICAL :: file_check !<
            !! Generic logical used for the purpose of asserting whether a file
            !! is or is not present in the designated location
            
            ! Namelist for all of the parameters to be inputed by the user
            NAMELIST /user_inputs/ case_dir, m, n, p, t_step_start,           &
                                   t_step_stop, t_step_save, model_eqns,      &
                                   num_fluids, mpp_lim, adv_alphan,           &
                                   weno_order, bc_x,                          &
                                   bc_y, bc_z, fluid_pp, format, precision,   &
                                   alpha_rho_wrt, rho_wrt, mom_wrt, vel_wrt,  &
                                   E_wrt, pres_wrt, alpha_wrt, gamma_wrt,     &
                                   heat_ratio_wrt, pi_inf_wrt, pres_inf_wrt,  &
                                   cons_vars_wrt, prim_vars_wrt, c_wrt,       &
                                   omega_wrt, schlieren_wrt, schlieren_alpha, &
                                   fd_order, mixture_err, alt_soundspeed,     &
                                   kappa_wrt, flux_lim, flux_wrt, cyl_coord,  &
                                   parallel_io, coarsen_silo, fourier_decomp, &
                                   fourier_modes,                             &
                                   rhoref, pref, bubbles, R0ref, nb,          &
                                   polytropic, thermal, Ca, Web, Re_inv,      &
                                   polydisperse, poly_sigma

            
            ! Inquiring the status of the post_process.inp file
            file_loc = 'post_process.inp'
            INQUIRE(FILE = TRIM(file_loc), EXIST = file_check)
            
            
            ! Checking whether the input file is there. If it is, the input file
            ! is read. If not, the program is terminated.
            IF(file_check) THEN
                OPEN(1, FILE = TRIM(file_loc), FORM = 'formatted', &
                        STATUS = 'old', ACTION = 'read')
                READ(1, NML = user_inputs)
                CLOSE(1)
                ! Store m,n,p into global m,n,p
                m_glb = m
                n_glb = n
                p_glb = p
            ELSE
                PRINT '(A)', 'File post_process.inp is missing. Exiting ...'
                CALL s_mpi_abort()
            END IF
            
            
        END SUBROUTINE s_read_input_file ! -------------------------------------
        
        
        
        
        !>  Checking that the user inputs make sense, i.e. that the
        !!      individual choices are compatible with the code's options
        !!      and that the combination of these choices results into a
        !!      valid configuration for the post-process        
        SUBROUTINE s_check_input_file() ! --------------------------------------

            
            CHARACTER(LEN = LEN_TRIM(case_dir)) :: file_loc !<
            !! Generic string used to store the address of a particular file
            

            LOGICAL :: dir_check !<
            !! Logical variable used to test the existence of folders
            
            INTEGER :: i  !< Generic loop iterator
            INTEGER :: bub_fac

            bub_fac = 0; 
            IF (bubbles .AND. (num_fluids == 1)) bub_fac = 1
            
            ! Checking the existence of the case folder
            case_dir = ADJUSTL(case_dir)
            
            file_loc = TRIM(case_dir) // '/.'
            
            CALL my_inquire(file_loc,dir_check)
            
            ! Constraint on the location of the case directory
            IF(dir_check .NEQV. .TRUE.) THEN
                PRINT '(A)', 'Unsupported choice for the value of '         // &
                             'case_dir. Exiting ...'
                CALL s_mpi_abort()
             
            ! Constraints on dimensionality and the number of cells for the grid
            ELSEIF(m <= 0) THEN
                PRINT '(A)', 'Unsupported choice for the value of m. '      // &
                             'Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(n < 0) THEN
                PRINT '(A)', 'Unsupported choice for the value of n. '      // &
                             'Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(p < 0) THEN
                PRINT '(A)', 'Unsupported choice for the value of p. '      // &
                             'Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(cyl_coord .AND. p > 0 .AND. MOD(p,2) /= 1) THEN
                PRINT '(A)', 'Unsupported choice for the value of p. '      // &
                             'Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(n == 0 .AND. p > 0) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for n and p. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(                    (m+1)*(n+1)*(p+1)                  &
                                               <                          &
                       2**(MIN(1,m) + MIN(1,n) + MIN(1,p))*num_procs      ) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for num_procs, m, n and p. '           // &
                             'Exiting ...'
                CALL s_mpi_abort()
                
            ! Constraints on the time-stepping parameters
            ELSEIF(t_step_start < 0) THEN
                PRINT '(A)', 'Unsupported choice for the value of '         // &
                             't_step_start. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(t_step_stop < t_step_start) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for t_step_start and t_step_stop. '    // &
                             'Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(t_step_save > t_step_stop - t_step_start) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for t_step_start, t_step_stop and '    // &
                             't_step_save. Exiting ...'
                CALL s_mpi_abort()
                
            ! Constraints on model equations and number of fluids in the flow
            ELSEIF(ALL(model_eqns /= (/1,2,3,4/))) THEN
                PRINT '(A)', 'Unsupported value of model_eqns. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(               num_fluids /= dflt_int              &
                                          .AND.                       &
                    (num_fluids < 1 .OR. num_fluids > num_fluids_max) ) THEN
                PRINT '(A)', 'Unsupported value of num_fluids. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF( (model_eqns == 1 .AND. num_fluids /= dflt_int) &
                                         .OR.                      &
                    (model_eqns == 2 .AND. num_fluids == dflt_int) &
                                         .OR.                      &
                    (model_eqns == 3 .AND. num_fluids == dflt_int)) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'model_eqns and num_fluids. '           // &
                             'Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(model_eqns == 1 .AND. adv_alphan) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'model_eqns and adv_alphan. '           // &
                             'Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(num_fluids == 1 .AND. (adv_alphan .NEQV. .TRUE.)) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'num_fluids and adv_alphan. '           // &
                             'Exiting ...'
                CALL s_mpi_abort()
                
            ! Constraints on the order of the WENO scheme
            ELSEIF(          weno_order /= 1 .AND. weno_order /= 3        &
                                             .AND.                        &
                                        weno_order /= 5                   ) THEN
                PRINT '(A)', 'Unsupported choice for the value of '         // &
                             'weno_order. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(m+1 < weno_order) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for m and weno_order. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(n > 0 .AND. n+1 < weno_order) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for n and weno_order. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(p > 0 .AND. p+1 < weno_order) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for p and weno_order. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(                    (m+1)*(n+1)*(p+1)                       &
                                               <                               &
                     weno_order**(MIN(1,m) + MIN(1,n) + MIN(1,p))*num_procs  ) &
                                             THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for num_procs, m, n, p and '           // &
                             'weno_order. Exiting ...'
                CALL s_mpi_abort()
                
            ! Constraints on the boundary conditions in the x-direction
            ELSEIF(bc_x%beg < -12 .OR. bc_x%beg > -1) THEN
                PRINT '(A)', 'Unsupported choice for the value of '         // &
                             'bc_x%beg. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(bc_x%end < -12 .OR. bc_x%end > -1) THEN
                PRINT '(A)', 'Unsupported choice for the value of '         // &
                             'bc_x%end. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(          (bc_x%beg == -1 .AND. bc_x%end /= -1)        &
                                             .OR.                         &
                             (bc_x%end == -1 .AND. bc_x%beg /= -1)        ) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for bc_x%beg and bc_x%end. '           // &
                             'Exiting ...'
                CALL s_mpi_abort()
                
            ! Constraints on the boundary conditions in the y-direction
            ELSEIF(         bc_y%beg /= dflt_int        &
                                    .AND.               &
                   ((   ( (cyl_coord .NEQV. .TRUE.)     &
                                    .OR.                &
                          (cyl_coord .AND. p == 0))     &
                                    .AND.               &
                    (bc_y%beg < -12 .OR. bc_y%beg > -1))&
                                    .OR.                &
                    (      cyl_coord .AND. p > 0        &
                                    .AND.               &
                    (bc_y%beg < -13 .OR. bc_y%beg > -1)))) THEN
                PRINT '(A)', 'Unsupported choice for the value of '         // &
                             'bc_y%beg. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(                 bc_y%end /= dflt_int               &
                                             .AND.                        &
                             (bc_y%end < -12 .OR. bc_y%end > -1)          ) THEN
                PRINT '(A)', 'Unsupported choice for the value of '         // &
                             'bc_y%end. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(         (n == 0 .AND. bc_y%beg /= dflt_int)        &
                                             .OR.                         &
                            (n  > 0 .AND. bc_y%beg == dflt_int)        ) THEN
                PRINT '(A)', 'Unsupported choice for the value of n and '   // &
                             'bc_y%beg. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(         (n == 0 .AND. bc_y%end /= dflt_int)        &
                                             .OR.                         &
                            (n  > 0 .AND. bc_y%end == dflt_int)        ) THEN
                PRINT '(A)', 'Unsupported choice for the value of n and '   // &
                             'bc_y%end. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(                          n > 0                        &
                                             .AND.                        &
                           ( (bc_y%beg == -1 .AND. bc_y%end /= -1)        &
                                             .OR.                         &
                             (bc_y%end == -1 .AND. bc_y%beg /= -1) )      ) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for n, bc_y%beg and bc_y%end. '        // &
                             'Exiting ...'
                CALL s_mpi_abort()
                
            ! Constraints on the boundary conditions in the z-direction
            ELSEIF(                 bc_z%beg /= dflt_int               &
                                             .AND.                        &
                             (bc_z%beg < -12 .OR. bc_z%beg > -1)          ) THEN
                PRINT '(A)', 'Unsupported choice for the value of '         // &
                             'bc_z%beg. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(                 bc_z%end /= dflt_int               &
                                             .AND.                        &
                             (bc_z%end < -12 .OR. bc_z%end > -1)          ) THEN
                PRINT '(A)', 'Unsupported choice for the value of '         // &
                             'bc_z%end. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(         (p == 0 .AND. bc_z%beg /= dflt_int)        &
                                             .OR.                         &
                            (p  > 0 .AND. bc_z%beg == dflt_int)        ) THEN
                PRINT '(A)', 'Unsupported choice for the value of p and '   // &
                             'bc_z%beg. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(         (p == 0 .AND. bc_z%end /= dflt_int)        &
                                             .OR.                         &
                            (p  > 0 .AND. bc_z%end == dflt_int)        ) THEN
                PRINT '(A)', 'Unsupported choice for the value of p and '   // &
                             'bc_z%end. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(                          p > 0                        &
                                             .AND.                        &
                           ( (bc_z%beg == -1 .AND. bc_z%end /= -1)        &
                                          .OR.                            &
                             (bc_z%end == -1 .AND. bc_z%beg /= -1) )      ) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for p, bc_z%beg and bc_z%end. '        // &
                             'Exiting ...'
                CALL s_mpi_abort()

            ! Constraints on Fourier decomposition options
            ELSEIF (        fourier_decomp          &
                                .AND.               &
                       (cyl_coord .NEQV. .TRUE.     &
                                .OR.                &
                               p == 0          )) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '   // &
                             'fourier_decomp and cyl_coord or value of p. '// &
                             'Exiting ...'
                CALL s_mpi_abort()
            ELSEIF (        fourier_decomp          &
                                .AND.               &
                   ( fourier_modes%beg == dflt_int  &
                                .OR.                &
                     fourier_modes%end == dflt_int )) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '   // &
                             'fourier_decomp and fourier_modes. '// &
                             'Exiting ...'
                CALL s_mpi_abort()
            END IF
            
            
            ! Constraints on the stiffened equation of state fluids parameters
            DO i = 1, num_fluids_max
                
                IF( fluid_pp(i)%gamma /= dflt_real &
                                 .AND.             &
                    fluid_pp(i)%gamma <=    0d0    ) THEN
                    PRINT '(A,I0,A)', 'Unsupported value of ' // &
                                      'fluid_pp(',i,')%'      // &
                                      'gamma. Exiting ...'
                    CALL s_mpi_abort()
                ELSEIF(         model_eqns == 1         &
                                     .AND.              &
                        fluid_pp(i)%gamma /= dflt_real ) THEN
                    PRINT '(A,I0,A)', 'Unsupported combination ' // &
                                      'of values of model_eqns ' // &
                                      'and fluid_pp(',i,')%'     // &
                                      'gamma. Exiting ...'
                    CALL s_mpi_abort()
                ELSEIF((i <= num_fluids+bub_fac .AND. fluid_pp(i)%gamma <=    0d0   )  &
                                                .OR.                           &
                       (i >  num_fluids+bub_fac .AND. fluid_pp(i)%gamma /= dflt_real)) &
                                                THEN
                    PRINT '(A,I0,A)', 'Unsupported combination ' // &
                                      'of values of num_fluids ' // &
                                      'and fluid_pp(',i,')%'     // &
                                      'gamma. Exiting ...'
                    CALL s_mpi_abort()
                ELSEIF( fluid_pp(i)%pi_inf /= dflt_real &
                                     .AND.              &
                        fluid_pp(i)%pi_inf <     0d0    ) THEN
                    PRINT '(A,I0,A)', 'Unsupported value of ' // &
                                      'fluid_pp(',i,')%'      // &
                                      'pi_inf. Exiting ...'
                    CALL s_mpi_abort()
                ELSEIF(         model_eqns == 1         &
                                     .AND.              &
                        fluid_pp(i)%pi_inf /= dflt_real ) THEN
                    PRINT '(A,I0,A)', 'Unsupported combination ' // &
                                      'of values of model_eqns ' // &
                                      'and fluid_pp(',i,')%'     // &
                                      'pi_inf. Exiting ...'
                    CALL s_mpi_abort()
                ELSEIF((i <= num_fluids+bub_fac .AND. fluid_pp(i)%pi_inf <     0d0   ) &
                                                .OR.                           &
                       (i >  num_fluids+bub_fac .AND. fluid_pp(i)%pi_inf /= dflt_real))&
                                                THEN
                    PRINT '(A,I0,A)', 'Unsupported combination ' // &
                                      'of values of num_fluids ' // &
                                      'and fluid_pp(',i,')%'     // &
                                      'pi_inf. Exiting ...'
                    CALL s_mpi_abort()
                END IF
                
            END DO
            
            
            ! Constraints on the format of the formatted database file(s)
            IF(format /= 1 .AND. format /= 2) THEN
                PRINT '(A)', 'Unsupported choice for the value of format. ' // &
                             'Exiting ...'
                CALL s_mpi_abort()
            ELSEIF( (precision /= 2) .AND. (parallel_io .NEQV. .FALSE.) ) THEN
                PRINT '(A)', 'Unsupported combination of precision and parallel IO. '         // &
                             'Please use precision == 2 when enabling parallel_io.  Exiting ...'
                CALL s_mpi_abort()

            ! Constraints on the precision of the formatted database file(s)
            ELSEIF(precision /= 1 .AND. precision /= 2) THEN
                PRINT '(A)', 'Unsupported choice for the value of '         // &
                             'precision. Exiting ...'
                CALL s_mpi_abort()

            ! Constraints on the option to coarsen the formatted database files
            ELSEIF( coarsen_silo .AND. format /= 1) THEN
                PRINT '(A)', 'Unsupported combination of values of format ' // &
                             'and coarsen_silo. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF( coarsen_silo .AND. n == 0) THEN
                PRINT '(A)', 'Unsupported combination of values of n ' // &
                             'and coarsen_silo. Exiting ...'
                CALL s_mpi_abort()
            END IF
            
            
            ! Constraints on the post-processing of the partial densities
            DO i = 1, num_fluids_max
                IF(      (  (i >  num_fluids .OR.  model_eqns == 1)       &
                                             .AND.                        &
                                       alpha_rho_wrt(i)              )    &
                                             .OR.                         &
                         (  (i <= num_fluids .AND. model_eqns == 1)       &
                                             .AND.                        &
                                       alpha_rho_wrt(i)              )    ) THEN
                    PRINT '(A,I0,A)', 'Unsupported choice of the '          // &
                                      'combination of values for '          // &
                                      'model_eqns, num_fluids and '         // &
                                      'alpha_rho_wrt(',i,'). Exiting ...'
                    CALL s_mpi_abort()
                END IF
            END DO
            
            
            ! Constraints on the post-processing of the momentum
            IF(n == 0 .AND. mom_wrt(2)) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for n and mom_wrt(2). Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(n == 0 .AND. mom_wrt(3)) THEN
                PRINT '(A)', 'Unsupported cohice of the combination of '    // &
                             'values for n and mom_wrt(3). Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(p == 0 .AND. mom_wrt(3)) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for p and mom_wrt(3). Exiting ...'
                CALL s_mpi_abort()
                
            ! Constraints on the post-processing of the velocity
            ELSEIF(n == 0 .AND. vel_wrt(2)) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for n and vel_wrt(2). Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(n == 0 .AND. vel_wrt(3)) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for n and vel_wrt(3). Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(p == 0 .AND. vel_wrt(3)) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for p and vel_wrt(3). Exiting ...'
                CALL s_mpi_abort()
            END IF

            ! Constraints on the post-processing of the flux limiter function
            IF ( n == 0 .AND. flux_wrt(2)) THEN
                PRINT '(A)', 'Unsupported choice of the combination of ' // &
                        'values for n and flux_wrt(2). Exiting ...'
                CALL s_mpi_abort()
            ELSEIF (n == 0 .AND. flux_wrt(3)) THEN
                PRINT '(A)', 'Unsupported choice of the combination of ' // &
                        'values for n and flux_wrt(3). Exiting ...'
                CALL s_mpi_abort()
            ELSEIF (p == 0 .AND. flux_wrt(3)) THEN
                PRINT '(A)', 'Unsupported choice of the combination of ' // &
                        'values for p and flux_wrt(3). Exiting ...'
                CALL s_mpi_abort()
            ELSEIF( ALL(flux_lim /= (/dflt_int,1,2,3,4,5,6,7/))) THEN
                PRINT '(A)', 'Unsupported value of flux_lim. Exiting ...'
                CALL s_mpi_abort()
            END IF
            
           
            ! Constraints on the post-processing of the volume fractions
            DO i = 1, num_fluids_max
                IF(      (  (i >  num_fluids .OR.  model_eqns == 1)       &
                                             .AND.                        &
                                          alpha_wrt(i)               )    &
                                             .OR.                         &
                         (  (i <= num_fluids .AND. model_eqns == 1)       &
                                             .AND.                        &
                                          alpha_wrt(i)               )    ) THEN
                    PRINT '(A,I0,A)', 'Unsupported choice of the '          // &
                                      'combination of values for '          // &
                                      'model_eqns, num_fluids and '         // &
                                      'alpha_wrt(',i,'). Exiting ...'
                    CALL s_mpi_abort()
                END IF
            END DO
            
            
            ! Constraints on the post-processing of the vorticity
            IF(n == 0 .AND. omega_wrt(1)) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for n and omega_wrt(1). Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(n == 0 .AND. omega_wrt(2)) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for n and omega_wrt(2). Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(n == 0 .AND. omega_wrt(3)) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for n and omega_wrt(3). Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(p == 0 .AND. omega_wrt(1)) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for p and omega_wrt(1). Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(p == 0 .AND. omega_wrt(2)) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for p and omega_wrt(2). Exiting ...'
                CALL s_mpi_abort()
                
            ! Constraints on post-processing of numerical Schlieren function
            ELSEIF(n == 0 .AND. schlieren_wrt) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for n and schlieren_wrt. Exiting ...'
                CALL s_mpi_abort()
                
            ! Constraints on post-processing combination of flow variables
            ELSEIF(           (ANY(alpha_rho_wrt) .NEQV. .TRUE.)          &
                                             .AND.                        &
                                 (ANY(mom_wrt) .NEQV. .TRUE.)             &
                                             .AND.                        &
                                 (ANY(vel_wrt) .NEQV. .TRUE.)             &
                                             .AND.                        &
                                 (ANY(flux_wrt) .NEQV. .TRUE.)            &
                                             .AND.                        &
                       (ANY((/ rho_wrt, E_wrt, pres_wrt,                  &
                               gamma_wrt, heat_ratio_wrt,                 &
                               pi_inf_wrt, pres_inf_wrt,                  &
                               cons_vars_wrt,                             &
                               prim_vars_wrt,                             &
                               c_wrt, schlieren_wrt /)) .NEQV. .TRUE.)    &
                                             .AND.                        &
                                 (ANY(alpha_wrt) .NEQV. .TRUE.)           &
                                             .AND.                        &
                                 (ANY(omega_wrt) .NEQV. .TRUE.)           ) THEN
                PRINT '(A)', 'None of the flow variables have been '        // &
                             'selected for post-process. Exiting ...'
                CALL s_mpi_abort()
            END IF
            
            
            ! Constraints on the coefficients of numerical Schlieren function
            DO i = 1, num_fluids_max
                IF(               schlieren_alpha(i) /= dflt_real         &
                                             .AND.                        &
                                   schlieren_alpha(i) <= 0d0              ) THEN
                    PRINT '(A,I0,A)', 'Unsupported choice for the value of '// &
                                      'schlieren_alpha(',i,'). Exiting ...'
                    CALL s_mpi_abort()
                ELSEIF(((i > num_fluids .OR. (schlieren_wrt .NEQV. .TRUE.)) &
                                             .AND.                        &
                                  schlieren_alpha(i) /= dflt_real       ) &
                                             .OR.                         &
                       ((      i <= num_fluids .AND. schlieren_wrt      ) &
                                             .AND.                        &
                                   schlieren_alpha(i) <= 0d0             )) THEN
                    PRINT '(A,I0,A)', 'Unsupported choice of the '          // &
                                      'combination of values for '          // &
                                      'num_fluids, schlieren_wrt and '      // &
                                      'schlieren_alpha(',i,'). Exiting ...'
                    CALL s_mpi_abort()
                END IF
            END DO
            
            
            ! Constraints on the order of the finite difference scheme
            IF(                     fd_order /= dflt_int                  &
                                             .AND.                        &
                  fd_order /= 1 .AND. fd_order /= 2 .AND. fd_order /= 4   ) THEN
                PRINT '(A)', 'Unsupported choice for the value of '         // &
                             'fd_order. Exiting ...'
                CALL s_mpi_abort()
  !          ELSEIF(               (omega_wrt(1) .NEQV. .TRUE.)            &
  !                                           .AND.                        &
  !                                (omega_wrt(2) .NEQV. .TRUE.)            &
  !                                           .AND.                        &
  !                                (omega_wrt(3) .NEQV. .TRUE.)            &
  !                                           .AND.                        &
  !                                !(schlieren_wrt .NEQV. .TRUE.)           &
  !                                !           .AND.                        &
  !                                   fd_order /= dflt_int              ) THEN
  !              PRINT '(A)', 'AA Unsupported choice of the combination of '    // &
  !                           'values for omega_wrt, schlieren_wrt and '     // &
  !                           'fd_order. Exiting ...'
  !              CALL s_mpi_abort()
            ELSEIF(         (ANY(omega_wrt) .OR. schlieren_wrt)           &
                                            .AND.                         &
                                    fd_order == dflt_int        ) THEN
                PRINT '(A)', 'BB Unsupported choice of the combination of '    // &
                             'values for omega_wrt, schlieren_wrt and '     // &
                             'fd_order. Exiting ...'
                CALL s_mpi_abort()
            END IF
            
            
        END SUBROUTINE s_check_input_file ! ------------------------------------
        
        
        
        
        
END MODULE m_start_up
