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
!! @brief Contains module m_start_up
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief This module contains subroutines that read, and check consistency
!!              of, the user provided inputs and data.
MODULE m_start_up
    
    ! Dependencies =============================================================
    USE m_derived_types          !< Definitions of the derived types

    USE m_global_parameters      !< Global parameters for the code

    USE m_mpi_proxy              !< Message passing interface (MPI) module proxy

    USE m_data_output            !< Procedures to write the grid data and the
                                 !! conservative variables to files

    USE mpi                      !< Message passing interface (MPI) module

    USE m_compile_specific
    ! ==========================================================================
    
    
    IMPLICIT NONE
    
    PRIVATE; PUBLIC :: s_initialize_start_up_module, &
                       s_read_input_file, &
                       s_check_input_file, &
                       s_read_grid_data_files, &
                       s_read_ic_data_files, &
                       s_read_serial_grid_data_files, &
                       s_read_serial_ic_data_files, &
                       s_read_parallel_grid_data_files, &
                       s_read_parallel_ic_data_files, &
                       s_check_grid_data_files, &
                       s_finalize_start_up_module 

    ABSTRACT INTERFACE ! ===================================================

        SUBROUTINE s_read_abstract_grid_data_files(dflt_int)! ----------

            INTEGER, INTENT(IN) :: dflt_int

        END SUBROUTINE s_read_abstract_grid_data_files ! ---------------

        SUBROUTINE s_read_abstract_ic_data_files(q_cons_vf) ! -----------

            IMPORT :: scalar_field, sys_size

            ! Conservative variables
            TYPE(scalar_field), &
            DIMENSION(sys_size), &
            INTENT(INOUT) :: q_cons_vf

        END SUBROUTINE s_read_abstract_ic_data_files ! -----------------

    END INTERFACE ! ========================================================


    CHARACTER(LEN = path_len + name_len) :: proc_rank_dir !<
    !! Location of the folder associated with the rank of the local processor
    
    CHARACTER(LEN = path_len + 2*name_len), PRIVATE :: t_step_dir !<
    !! Possible location of time-step folder containing preexisting grid and/or
    !! conservative variables data to be used as starting point for pre-process    

    PROCEDURE(s_read_abstract_grid_data_files), POINTER :: s_read_grid_data_files => NULL()
    PROCEDURE(s_read_abstract_ic_data_files), POINTER :: s_read_ic_data_files => NULL()
    
    CONTAINS
        
        
        
        !>  Reads the configuration file pre_process.inp, in order to
        !!      populate the parameters in module m_global_parameters.f90
        !!      with the user provided inputs        
        SUBROUTINE s_read_input_file() ! ---------------------------------------


            CHARACTER(LEN = name_len) :: file_loc  !<
            !! Generic string used to store the address of a particular file
            
            LOGICAL :: file_check !<
            !! Generic logical used for the purpose of asserting whether a file
            !! is or is not present in the designated location
             
            ! Namelist for all of the parameters to be inputed by the user
            NAMELIST /user_inputs/ case_dir, old_grid, old_ic, t_step_old, m, &
                                   n, p, x_domain, y_domain, z_domain,        &
                                   stretch_x, stretch_y, stretch_z, a_x, a_y, &
                                   a_z, x_a, y_a, z_a, x_b, y_b, z_b,         &
                                   model_eqns, num_fluids,                    &
                                   adv_alphan, mpp_lim,                       &
                                   weno_order, bc_x, bc_y, bc_z, num_patches, &
                                   hypoelasticity, patch_icpp, fluid_pp,      &
                                   precision, parallel_io,                    &
                                   perturb_flow, perturb_flow_fluid,          &
                                   perturb_sph, perturb_sph_fluid, fluid_rho, &
                                   cyl_coord, loops_x, loops_y, loops_z,      &
                                   rhoref, pref, bubbles, R0ref, nb,          &
                                   polytropic, thermal, Ca, Web, Re_inv,      &
                                   polydisperse, poly_sigma, qbmm,      &
                                   nnode, sigR, sigV, dist_type, rhoRV, R0_type
 

            ! Inquiring the status of the pre_process.inp file
            file_loc = 'pre_process.inp'
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
                PRINT '(A)', 'File pre_process.inp is missing. Exiting ...'
                CALL s_mpi_abort()
            END IF
            
            
        END SUBROUTINE s_read_input_file ! -------------------------------------
        
        
        
        
        !>  Checking that the user inputs make sense, i.e. that the
        !!      individual choices are compatible with the code's options
        !!      and that the combination of these choices results into a
        !!      valid configuration for the pre-process       
        SUBROUTINE s_check_input_file() ! --------------------------------------
                     

            CHARACTER(LEN = LEN_TRIM(case_dir)) :: file_loc !<
            !! Generic string used to store the address of a particular file
            

            LOGICAL :: dir_check !<
            !! Logical variable used to test the existence of folders
            

            INTEGER :: i !< 
            !! Generic loop iterator

            INTEGER :: bub_fac !<
            !! For allowing an extra fluid_pp if there are subgrid bubbles

            bub_fac = 0
            IF (bubbles .AND. (num_fluids == 1)) bub_fac = 1
            
            
            ! Checking the existence of the case folder
            case_dir = ADJUSTL(case_dir)
            
            file_loc = TRIM(case_dir) // '/.'
           
            CALL my_inquire(file_loc,dir_check)
            
            ! Startup checks for bubbles and bubble variables
            IF(bubbles .AND. (model_eqns .NE. 4 .AND. model_eqns .NE. 2)) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'bubbles and model_eqns. '           // &
                             'Exiting ...'
                CALL s_mpi_abort()     
            ELSEIF(bubbles .AND. polydisperse .and. (nb==1)) THEN
                PRINT '(A)', 'Polydisperse bubble dynamics requires nb > 1 ' // &
                             'Exiting ...'
                CALL s_mpi_abort()      
            ELSEIF(bubbles .AND. polydisperse .and. (mod(nb,2)==0)) THEN
                PRINT '(A)', 'nb must be odd ' // &
                             'Exiting ...'
                CALL s_mpi_abort()      
            ELSEIF(model_eqns == 4 .AND. (rhoref == dflt_real)) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'bubbles and rhoref. '           // &
                             'Exiting ...'
                CALL s_mpi_abort()              
            ELSEIF(model_eqns == 4 .AND. (pref  == dflt_real)) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'bubbles and pref. '           // &
                             'Exiting ...'
                CALL s_mpi_abort()  
            ELSEIF(model_eqns == 4 .AND. (num_fluids > 1)) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'model_eqns and num_fluids. '           // &
                             'Exiting ...'
                CALL s_mpi_abort()   
            ELSEIF(bubbles .AND. (R0ref == dflt_real)) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'bubbles and R0ref. '           // &
                             'Exiting ...'
                CALL s_mpi_abort()   
            ELSEIF(bubbles .AND. (nb == dflt_int)) THEN
                print '(a)', 'unsupported combination of values of ' // &
                             'bubbles and nb. '           // &
                             'exiting ...'
                CALL s_mpi_abort()  
            ELSEIF(bubbles .AND. (thermal > 3)) THEN
                print '(a)', 'unsupported combination of values of ' // &
                             'bubbles and thermal. '           // &
                             'exiting ...'
                CALL s_mpi_abort() 
            END IF


            ! Constraint on the location of the case directory
            IF(dir_check .NEQV. .TRUE.) THEN
                PRINT '(A)', 'Unsupported choice for the value of '         // &
                             'case_dir.'
                PRINT '(A)', 'WARNING: Ensure that compiler flags/choices in Makefiles match your compiler! '
                PRINT '(A)', 'WARNING: Ensure that preprocessor flags are enabled! '
                CALL s_mpi_abort()
                
            ! Constraints on the use of a preexisting grid and initial condition
            ELSEIF((old_grid .NEQV. .TRUE.) .AND. old_ic) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for old_grid and old_ic. Exiting ...'
                CALL s_mpi_abort()
                
            ELSEIF((old_grid .OR. old_ic) .AND. t_step_old == dflt_int) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                    'values for old_grid and old_ic and t_step_old. Exiting ...'
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
                             'Total number of cells in azimuthal direction '// &
                             'must be an even number. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(n == 0 .AND. p > 0) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for n and p. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(                     (m+1)*(n+1)*(p+1)                 &
                                                <                         &
                        2**(MIN(1,m) + MIN(1,n) + MIN(1,p))*num_procs     ) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for num_procs, m, n and p. '           // &
                             'Exiting ...'
                CALL s_mpi_abort()
                
            ! Constraints on domain boundaries locations in the x-direction
            ELSEIF(     (old_grid .AND. x_domain%beg /= dflt_real)     &
                                             .OR.                         &
                        (       (old_grid .NEQV. .TRUE.) .AND.            &
                                 x_domain%beg == dflt_real       )     ) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for old_grid and x_domain%beg. '       // &
                             'Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(     (old_grid .AND. x_domain%end /= dflt_real)     &
                                             .OR.                         &
                        (       (old_grid .NEQV. .TRUE.) .AND.            &
                                 x_domain%end == dflt_real       )     ) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for old_grid and x_domain%end. '       // &
                             'Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(                (old_grid .NEQV. .TRUE.)               &
                                             .AND.                        &
                                 x_domain%beg >= x_domain%end             ) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for old_grid, x_domain%beg and '       // &
                             'x_domain%end. Exiting ...'
                CALL s_mpi_abort()
            ELSE IF (qbmm .and. dist_type == dflt_int) THEN
                PRINT '(A)', 'Dist type must be set if using QBMM. Exiting ...'
                CALL s_mpi_abort()
            ELSE IF (qbmm .and. (dist_type .NE. 1) .and. rhoRV > 0d0) THEN
                PRINT '(A)', 'rhoRV cannot be used with dist_type \ne 1. Exiting ...'
                CALL s_mpi_abort()
            ELSE IF (polydisperse .and. R0_type == dflt_int) THEN
                PRINT '(A)', 'R0 type must be set if using Polydisperse. Exiting ...'
                CALL s_mpi_abort()
            END IF
                
            IF (cyl_coord .NEQV. .TRUE.) THEN ! Cartesian coordinates

                ! Constraints on domain boundaries locations in the y-direction
                IF( (     n == 0 .AND. y_domain%beg /= dflt_real     )        &
                                                 .OR.                         &
                        (                        n > 0                        &
                                                 .AND.                        &
                          ( (old_grid .AND. y_domain%beg /= dflt_real)     &
                                                 .OR.                         &
                            (       (old_grid .NEQV. .TRUE.) .AND.            &
                                     y_domain%beg == dflt_real       ) ) ) ) THEN
                    PRINT '(A)', 'Unsupported choice of the combination of '    // &
                                 'values for old_grid, n and y_domain%beg. '    // &
                                 'Exiting ...'
                    CALL s_mpi_abort()
                ELSEIF( (     n == 0 .AND. y_domain%end /= dflt_real     ) &
                                                 .OR.                         &
                        (                        n > 0                        &
                                                 .AND.                        &
                          ( (old_grid .AND. y_domain%end /= dflt_real)     &
                                                 .OR.                         &
                            (       (old_grid .NEQV. .TRUE.) .AND.            &
                                     y_domain%end == dflt_real       ) ) ) ) THEN
                    PRINT '(A)', 'Unsupported choice of the combination of '    // &
                                 'values for old_grid, n and y_domain%end. '    // &
                                 'Exiting ...'
                    CALL s_mpi_abort()
                ELSEIF(                          n > 0                        &
                                                 .AND.                        &
                                        (old_grid .NEQV. .TRUE.)              &
                                                 .AND.                        &
                                      y_domain%beg >= y_domain%end            ) THEN
                    PRINT '(A)', 'Unsupported choice of the combination of '    // &
                                 'values for old_grid, n, y_domain%beg and '    // &
                                 'y_domain%end. Exiting ...'
                    CALL s_mpi_abort()
                    
                ! Constraints on domain boundaries locations in the z-direction
                ELSEIF( (     p == 0 .AND. z_domain%beg /= dflt_real     ) &
                                                 .OR.                         &
                        (                        p > 0                        &
                                                 .AND.                        &
                          ( (old_grid .AND. z_domain%beg /= dflt_real)     &
                                                 .OR.                         &
                            (       (old_grid .NEQV. .TRUE.) .AND.            &
                                     z_domain%beg == dflt_real       ) ) ) ) THEN
                    PRINT '(A)', 'Unsupported choice of the combination of '    // &
                                 'values for old_grid, p and z_domain%beg. '    // &
                                 'Exiting ...'
                    CALL s_mpi_abort()
                ELSEIF( (     p == 0 .AND. z_domain%end /= dflt_real     ) &
                                                 .OR.                         &
                        (                         p > 0                       &
                                                 .AND.                        &
                          ( (old_grid .AND. z_domain%end /= dflt_real)     &
                                                 .OR.                         &
                            (       (old_grid .NEQV. .TRUE.) .AND.            &
                                     z_domain%end == dflt_real       ) ) ) ) THEN
                    PRINT '(A)', 'Unsupported choice of the combination of '    // &
                                 'values for old_grid, p and z_domain%end. '    // &
                                 'Exiting ...'
                    CALL s_mpi_abort()
                ELSEIF(                          p > 0                        &
                                                 .AND.                        &
                                        (old_grid .NEQV. .TRUE.)              &
                                                 .AND.                        &
                                      z_domain%beg >= z_domain%end            ) THEN
                    PRINT '(A)', 'Unsupported choice of the combination of '    // &
                                 'values for old_grid, p, z_domain%beg and '    // &
                                 'z_domain%end. Exiting ...'
                    CALL s_mpi_abort()
                END IF

            ELSE ! Cylindrical coordinates

                ! Constraints on domain boundaries for cylindrical coordinates
                IF(               n == 0                  &
                                   .OR.                   &
                          y_domain%beg /= 0d0             &
                                   .OR.                   &
                          y_domain%end == dflt_real       &
                                   .OR.                   &
                          y_domain%end < 0d0              &
                                   .OR.                   &
                        y_domain%beg >= y_domain%end )  THEN
                    PRINT '(A)', 'Unsupported choice of the combination of '    // &
                                 'cyl_coord and n, y_domain%beg, or         '   // &
                                 'y_domain%end. Exiting ...'
                    CALL s_mpi_abort()
                ELSEIF ( (p == 0 .AND. z_domain%beg /= dflt_real) &
                                        .OR.                      &
                         (p == 0 .AND. z_domain%end /= dflt_real)) THEN
                    PRINT '(A)', 'Unsupported choice of the combination of '    // &
                                 'cyl_coord and p, z_domain%beg, or '           // &
                                 'z_domain%end. Exiting ...'
                    CALL s_mpi_abort()
                ELSEIF( p > 0 .AND. ( z_domain%beg /= 0d0          &
                                              .OR.                 &
                                      z_domain%end /= 2d0*pi )) THEN
                    PRINT '(A)', 'Unsupported choice of the combination of '    // &
                                 'cyl_coord and p, z_domain%beg, or '           // &
                                 'z_domain%end. Exiting ...'
                    CALL s_mpi_abort()
                END IF

            END IF

            ! Constraints on the grid stretching in the x-direction
            IF(old_grid .AND. stretch_x) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for old_grid and stretch_x. '          // &
                             'Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(stretch_x .AND. a_x == dflt_real) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for stretch_x and a_x. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(stretch_x .AND. x_a == dflt_real) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for stretch_x and x_a. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(stretch_x .AND. x_b == dflt_real) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for stretch_x and x_b. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(stretch_x .AND. x_a >= x_b) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for stretch_x, x_a and x_b. '          // &
                             'Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(                        stretch_x                        &
                                             .AND.                          &
                    (a_x + LOG(COSH(a_x*(x_domain%beg - x_a)))              &
                         + LOG(COSH(a_x*(x_domain%beg - x_b)))              &
                     - 2d0*LOG(COSH( 0.5d0*a_x*(x_b - x_a) )))/a_x <= 0d0 ) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for x_domain%beg, stretch_x, a_x, '    // &
                             'x_a, and x_b. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(                        stretch_x                        &
                                             .AND.                          &
                    (a_x + LOG(COSH(a_x*(x_domain%end - x_a)))              &
                         + LOG(COSH(a_x*(x_domain%end - x_b)))              &
                     - 2d0*LOG(COSH( 0.5d0*a_x*(x_b - x_a) )))/a_x <= 0d0 ) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for x_domain%end, stretch_x, a_x, '    // &
                             'x_a, and x_b. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(loops_x < 1) THEN
                PRINT '(A)', 'Unsupported choice for the value of loops_x. '// &
                             'Exiting ...'
                CALL s_mpi_abort()
                
            ! Constraints on the grid stretching in the y-direction
            ELSEIF(old_grid .AND. stretch_y) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for old_grid and stretch_y. '          // &
                             'Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(n == 0 .AND. stretch_y) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for n and stretch_y. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(stretch_y .AND. a_y == dflt_real) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for stretch_y and a_y. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(stretch_y .AND. y_a == dflt_real) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for stretch_y and y_a. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(stretch_y .AND. y_b == dflt_real) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for stretch_y and y_b. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(stretch_y .AND. y_a >= y_b) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for stretch_y, y_a and y_b. '          // &
                             'Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(                        stretch_y                        &
                                             .AND.                          &
                    (a_y + LOG(COSH(a_y*(y_domain%beg - y_a)))              &
                         + LOG(COSH(a_y*(y_domain%beg - y_b)))              &
                     - 2d0*LOG(COSH( 0.5d0*a_y*(y_b - y_a) )))/a_y <= 0d0 ) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for y_domain%beg, stretch_y, a_y, '    // &
                             'y_a, and y_b. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(                        stretch_y                        &
                                             .AND.                          &
                    (a_y + LOG(COSH(a_y*(y_domain%end - y_a)))              &
                         + LOG(COSH(a_y*(y_domain%end - y_b)))              &
                     - 2d0*LOG(COSH( 0.5d0*a_y*(y_b - y_a) )))/a_y <= 0d0 ) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for y_domain%end, stretch_y, a_y, '    // &
                             'y_a, and y_b. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(loops_y < 1) THEN
                PRINT '(A)', 'Unsupported choice for the value of loops_y. '// &
                             'Exiting ...'
                CALL s_mpi_abort()
                
            ! Constraints on the grid stretching in the z-direction
            ELSEIF(old_grid .AND. stretch_z) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for old_grid and stretch_z. '          // &
                             'Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(cyl_coord .AND. stretch_z) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for cyl_coord and stretch_z. '         // &
                             'Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(p == 0 .AND. stretch_z) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for p and stretch_z. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(stretch_z .AND. a_z == dflt_real) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for stretch_z and a_z. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(stretch_z .AND. z_a == dflt_real) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for stretch_z and z_a. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(stretch_z .AND. z_b == dflt_real) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for stretch_z and z_b. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(stretch_z .AND. z_a >= z_b) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for stretch_z, z_a and z_b. '          // &
                             'Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(                        stretch_z                        &
                                             .AND.                          &
                    (a_z + LOG(COSH(a_z*(z_domain%beg - z_a)))              &
                         + LOG(COSH(a_z*(z_domain%beg - z_b)))              &
                     - 2d0*LOG(COSH( 0.5d0*a_z*(z_b - z_a) )))/a_z <= 0d0 ) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for z_domain%beg, stretch_z, a_z, '    // &
                             'z_a, and z_b. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(                        stretch_z                        &
                                             .AND.                          &
                    (a_z + LOG(COSH(a_z*(z_domain%end - z_a)))              &
                         + LOG(COSH(a_z*(z_domain%end - z_b)))              &
                     - 2d0*LOG(COSH( 0.5d0*a_z*(z_b - z_a) )))/a_z <= 0d0 ) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for z_domain%end, stretch_z, a_z, '    // &
                             'z_a, and z_b. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(loops_z < 1) THEN
                PRINT '(A)', 'Unsupported choice for the value of loops_z. '// &
                             'Exiting ...'
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
                    (model_eqns == 2 .AND. num_fluids == dflt_int) ) THEN
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
            ELSEIF(         weno_order /= 1 .AND. weno_order /= 3         &
                                            .AND.                         &
                                       weno_order /= 5                    ) THEN
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
            ELSEIF(                     (m+1)*(n+1)*(p+1)                      &
                                                <                              &
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
            ELSEIF(         (bc_x%beg == -1 .AND. bc_x%end /= -1)         &
                                            .OR.                          &
                            (bc_x%end == -1 .AND. bc_x%beg /= -1)         ) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for bc_x%beg and bc_x%end. '           // &
                             'Exiting ...'
                CALL s_mpi_abort()
            END IF
                
            IF (cyl_coord .NEQV. .TRUE.) THEN ! Cartesian coordinates 

                ! Constraints on the boundary conditions in the y-direction
                IF(                       bc_y%beg /= dflt_real               &
                                                .AND.                         &
                                (bc_y%beg < -12 .OR. bc_y%beg > -1)           ) THEN
                    PRINT '(A)', 'Unsupported choice for the value of '         // &
                                 'bc_y%beg. Exiting ...'
                    CALL s_mpi_abort()
                ELSEIF(                bc_y%end /= dflt_real               &
                                                .AND.                         &
                                (bc_y%end < -12 .OR. bc_y%end > -1)           ) THEN
                    PRINT '(A)', 'Unsupported choice for the value of '         // &
                                 'bc_y%end. Exiting ...'
                    CALL s_mpi_abort()
                ELSEIF(        (n == 0 .AND. bc_y%beg /= dflt_real)        &
                                                .OR.                          &
                               (n  > 0 .AND. bc_y%beg == dflt_real)        ) THEN
                    PRINT '(A)', 'Unsupported choice for the value of n and '   // &
                                 'bc_y%beg. Exiting ...'
                    CALL s_mpi_abort()
                ELSEIF(        (n == 0 .AND. bc_y%end /= dflt_real)        &
                                                .OR.                          &
                               (n  > 0 .AND. bc_y%end == dflt_real)        ) THEN
                    PRINT '(A)', 'Unsupported choice for the value of n and '   // &
                                 'bc_y%end. Exiting ...'
                    CALL s_mpi_abort()
                ELSEIF(                         n > 0                         &
                                                .AND.                         &
                              ( (bc_y%beg == -1 .AND. bc_y%end /= -1)         &
                                                .OR.                          &
                                (bc_y%end == -1 .AND. bc_y%beg /= -1) )       ) THEN
                    PRINT '(A)', 'Unsupported choice of the combination of '    // &
                                 'values for n, bc_y%beg and bc_y%end. '        // &
                                 'Exiting ...'
                    CALL s_mpi_abort()
                    
                ! Constraints on the boundary conditions in the z-direction
                ELSEIF(                bc_z%beg /= dflt_real               &
                                                .AND.                         &
                                (bc_z%beg < -12 .OR. bc_z%beg > -1)           ) THEN
                    PRINT '(A)', 'Unsupported choice for the value of '         // &
                                 'bc_z%beg. Exiting ...'
                    CALL s_mpi_abort()
                ELSEIF(                bc_z%end /= dflt_real               &
                                                .AND.                         &
                                (bc_z%end < -12 .OR. bc_z%end > -1)           ) THEN
                    PRINT '(A)', 'Unsupported choice for the value of '         // &
                                 'bc_z%end. Exiting ...'
                    CALL s_mpi_abort()
                ELSEIF(        (p == 0 .AND. bc_z%beg /= dflt_real)        &
                                                .OR.                          &
                               (p  > 0 .AND. bc_z%beg == dflt_real)        ) THEN
                    PRINT '(A)', 'Unsupported choice for the value of p and '   // &
                                 'bc_z%beg. Exiting ...'
                    CALL s_mpi_abort()
                ELSEIF(        (p == 0 .AND. bc_z%end /= dflt_real)        &
                                                .OR.                          &
                               (p  > 0 .AND. bc_z%end == dflt_real)        ) THEN
                    PRINT '(A)', 'Unsupported choice for the value of p and '   // &
                                 'bc_z%end. Exiting ...'
                    CALL s_mpi_abort()
                ELSEIF(                        p > 0                          &
                                               .AND.                          &
                             ( (bc_z%beg == -1 .AND. bc_z%end /= -1)          &
                                               .OR.                           &
                               (bc_z%end == -1 .AND. bc_z%beg /= -1) )        ) THEN
                    PRINT '(A)', 'Unsupported choice of the combination of '    // &
                                 'values for p, bc_z%beg and bc_z%end. '        // &
                                 'Exiting ...'
                    CALL s_mpi_abort()
                END IF

            ELSE ! Cylindrical coordinates

                ! Constraints on the boundary conditions in the r-direction
                IF(                    bc_y%beg /= dflt_real               &
                                                .AND.                         &
                    (           ( p  > 0 .AND. bc_y%beg /= -13)               &
                                                .OR.                          &
                                ( p == 0 .AND. bc_y%beg /= -2 )  ))             THEN
                    PRINT '(A)', 'Unsupported choice for the value of '         // &
                                 'bc_y%beg. Exiting ...'
                    CALL s_mpi_abort()
                ELSEIF(                bc_y%end /= dflt_real               &
                                                .AND.                         &
                                (bc_y%end < -12 .OR. bc_y%end > -1)           ) THEN
                    PRINT '(A)', 'Unsupported choice for the value of '         // &
                                 'bc_y%end. Exiting ...'
                    CALL s_mpi_abort()
                ELSEIF(        (n  > 0 .AND. bc_y%beg == dflt_real)        ) THEN
                    PRINT '(A)', 'Unsupported choice for the value of n and '   // &
                                 'bc_y%beg. Exiting ...'
                    CALL s_mpi_abort()
                 ELSEIF(       (n  > 0 .AND. bc_y%end == dflt_real)        ) THEN
                    PRINT '(A)', 'Unsupported choice for the value of n and '   // &
                                 'bc_y%end. Exiting ...'
                    CALL s_mpi_abort()
                    
                ! Constraints on the boundary conditions in the theta-direction
                ELSEIF(                bc_z%beg /= dflt_real               &
                                                .AND.                         &
                                (bc_z%beg /= -1 .AND. bc_z%beg /= -2)         ) THEN
                    PRINT '(A)', 'Unsupported choice for the value of '         // &
                                 'bc_z%beg. Exiting ...'
                    CALL s_mpi_abort()
                ELSEIF(                bc_z%end /= dflt_real               &
                                                .AND.                         &
                                (bc_z%end /= -1 .AND. bc_z%end /= -2)         ) THEN
                    PRINT '(A)', 'Unsupported choice for the value of '         // &
                                 'bc_z%end. Exiting ...'
                    CALL s_mpi_abort()
                ELSEIF(        (p == 0 .AND. bc_z%beg /= dflt_real)        &
                                                .OR.                          &
                               (p  > 0 .AND. bc_z%beg == dflt_real)        ) THEN
                    PRINT '(A)', 'Unsupported choice for the value of p and '   // &
                                 'bc_z%beg. Exiting ...'
                    CALL s_mpi_abort()
                ELSEIF(        (p == 0 .AND. bc_z%end /= dflt_real)        &
                                                .OR.                          &
                               (p  > 0 .AND. bc_z%end == dflt_real)        ) THEN
                    PRINT '(A)', 'Unsupported choice for the value of p and '   // &
                                 'bc_z%end. Exiting ...'
                    CALL s_mpi_abort()
                ELSEIF(                        p > 0                          &
                                               .AND.                          &
                             ( (bc_z%beg == -1 .AND. bc_z%end /= -1)          &
                                               .OR.                           &
                               (bc_z%end == -1 .AND. bc_z%beg /= -1) )        ) THEN
                    PRINT '(A)', 'Unsupported choice of the combination of '    // &
                                 'values for p, bc_z%beg and bc_z%end. '        // &
                                 'Exiting ...'
                    CALL s_mpi_abort()
                END IF

            END IF
                
            ! Constraints on number of patches making up the initial condition
            IF(    num_patches < 0 .OR. num_patches > num_patches_max .OR. &
                  (num_patches == 0 .AND. t_step_old == dflt_int)) THEN
                PRINT '(A)', 'Unsupported choice for the value of '         // &
                             'num_patches. Exiting ...'
                CALL s_mpi_abort()
            ! Constraints on perturbing the initial condition
            ELSEIF ((perturb_flow .AND. perturb_flow_fluid == dflt_int) &
                                    .OR.                                &
                    ((perturb_flow .NEQV. .TRUE.) .AND. (perturb_flow_fluid /= dflt_int))) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '        // &
                             'values for perturb_flow and perturb_flow_fluid. ' // &
                             'Exiting ...'
                CALL s_mpi_abort()
            ELSEIF ( (perturb_flow_fluid > num_fluids)                  &
                                .OR.                                    &
                     (perturb_flow_fluid < 0 .AND. perturb_flow_fluid /= dflt_int)) THEN
                PRINT '(A)', 'Unsupported choice for the value of '        //&
                             'perturb_flow_fluid. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF ((perturb_sph .AND. perturb_sph_fluid == dflt_int) &
                                    .OR.                                &
                    ((perturb_sph .NEQV. .TRUE.) .AND. (perturb_sph_fluid /= dflt_int))) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '        // &
                             'values for perturb_sph and perturb_sph_fluid. ' // &
                             'Exiting ...'
                CALL s_mpi_abort()
            ELSEIF ( (perturb_sph_fluid > num_fluids)                  &
                                .OR.                                    &
                     (perturb_sph_fluid < 0 .AND. perturb_sph_fluid /= dflt_int)) THEN
                PRINT '(A)', 'Unsupported choice for the value of '        //&
                             'perturb_sph_fluid. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF( (ANY(fluid_rho /= dflt_real)) .AND. (perturb_sph .NEQV. .TRUE.)) THEN
                PRINT '(A)', 'Unsupported choices for values of perturb_sph '  //&
                             'and fluid_rho. Exiting ...'
                CALL s_mpi_abort()
            END IF
            
            IF (perturb_sph) THEN
                DO i = 1, num_fluids
                    IF (fluid_rho(i) == dflt_real) THEN
                        PRINT '(A,I0,A)', 'Unsupported choice for value of fluid_rho(', &
                                            i, '). Exiting ...'
                        CALL s_mpi_abort()
                    END IF
                END DO
            END IF
                                            
            
            ! Constraints on the geometric initial condition patch parameters
            DO i = 1, num_patches_max
                IF(i <= num_patches) THEN
                    IF(patch_icpp(i)%geometry == 1) THEN
                        CALL s_check_line_segment_patch_geometry(i)
                    ELSEIF(patch_icpp(i)%geometry == 2) THEN
                        CALL s_check_circle_patch_geometry(i)
                    ELSEIF(patch_icpp(i)%geometry == 3) THEN
                        CALL s_check_rectangle_patch_geometry(i)
                    ELSEIF(patch_icpp(i)%geometry == 4) THEN
                        CALL s_check_line_sweep_patch_geometry(i)
                    ELSEIF(patch_icpp(i)%geometry == 5) THEN
                        CALL s_check_ellipse_patch_geometry(i)
                    ELSEIF(patch_icpp(i)%geometry == 6) THEN
                        CALL s_check_isentropic_vortex_patch_geometry(i)
                    ELSEIF(patch_icpp(i)%geometry == 7) THEN
                        CALL s_check_2D_analytical_patch_geometry(i)
                    ELSEIF(patch_icpp(i)%geometry == 8) THEN
                        CALL s_check_sphere_patch_geometry(i)
                    ELSEIF(patch_icpp(i)%geometry == 9) THEN
                        CALL s_check_cuboid_patch_geometry(i)
                    ELSEIF(patch_icpp(i)%geometry == 10) THEN
                        CALL s_check_cylinder_patch_geometry(i)
                    ELSEIF(patch_icpp(i)%geometry == 11) THEN
                        CALL s_check_plane_sweep_patch_geometry(i)
                    ELSEIF(patch_icpp(i)%geometry == 12) THEN
                        CALL s_check_ellipsoid_patch_geometry(i)
                    ELSEIF(patch_icpp(i)%geometry == 13) THEN
                        CALL s_check_3D_analytical_patch_geometry(i)
                    ELSEIF(patch_icpp(i)%geometry == 14) THEN
                        CALL s_check_spherical_harmonic_patch_geometry(i)
                    ELSEIF(patch_icpp(i)%geometry == 15) THEN
                        CALL s_check_1d_analytical_patch_geometry(i)
                    ELSEIF(patch_icpp(i)%geometry == 16) THEN
                        print*, '1d pressure sinusoid'
                    ELSEIF(patch_icpp(i)%geometry == 17) THEN
                        print*, '2d spiral'
                    ELSEIF(patch_icpp(i)%geometry == 18) THEN
                        print*, '2d var circle'
                    ELSEIF(patch_icpp(i)%geometry == 19) THEN
                        print*, '3d var circle'
                    ELSE
                        PRINT '(A,I0,A)', 'Unsupported choice of the '      // &
                                          'geometry of active patch ',i,       &
                                          ' detected. Exiting ...'
                        CALL s_mpi_abort()
                    END IF
                ELSE
                    IF(patch_icpp(i)%geometry == dflt_int) THEN
                        CALL s_check_inactive_patch_geometry(i)
                    ELSE
                        PRINT '(A,I0,A)', 'Unsupported choice of the '      // &
                                          'geometry of inactive patch ',i,     &
                                          ' detected. Exiting ...'
                        CALL s_mpi_abort()
                    END IF
                END IF
            END DO
            
            
            ! Constraints on overwrite rights initial condition patch parameters
            DO i = 1, num_patches_max
                IF(i <= num_patches) THEN
                    CALL s_check_active_patch_alteration_rights(i)
                ELSE
                    CALL s_check_inactive_patch_alteration_rights(i)
                END IF
            END DO
            
            
            ! Constraints on smoothing initial condition patch parameters
            DO i = 1, num_patches_max
                IF(    i > 1 .AND. ( patch_icpp(i)%geometry == 2 .OR.     &
                                     patch_icpp(i)%geometry == 4 .OR.     &
                                     patch_icpp(i)%geometry == 5 .OR.     &
                                     patch_icpp(i)%geometry == 8 .OR.     &
                                     patch_icpp(i)%geometry == 10 .OR.     &
                                     patch_icpp(i)%geometry == 11 .OR.     &
                                     patch_icpp(i)%geometry == 12      )   ) THEN
                    CALL s_check_supported_patch_smoothing(i)
                ELSE
                    CALL s_check_unsupported_patch_smoothing(i)
                END IF
            END DO
            
            
            ! Constraints on flow variables initial condition patch parameters
            DO i = 1, num_patches_max
                IF(i <= num_patches) THEN
                    CALL s_check_active_patch_primitive_variables(i)
                ELSE
                    CALL s_check_inactive_patch_primitive_variables(i)
                END IF
            END DO
            
            
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
            
            
        END SUBROUTINE s_check_input_file ! ------------------------------------
        
        
        
        !> This subroutine verifies that the geometric parameters of
        !!      the line segment patch have consistently been inputted by
        !!      the user.      
        !!  @param patch_id Patch identifier
        SUBROUTINE s_check_line_segment_patch_geometry(patch_id) ! -------------

            INTEGER, INTENT(IN) :: patch_id
            
            ! Constraints on the geometric parameters of the line segment patch
            IF(n > 0 .OR. patch_icpp(patch_id)%length_x <= 0d0 &
                                    .OR.                       &
                patch_icpp(patch_id)%x_centroid == dflt_real   &
                                    .OR.                       &
                                 cyl_coord                     ) THEN
                
                PRINT '(A,I0,A)', 'Inconsistency(ies) detected in '         // &
                                  'geometric parameters of line segment '   // &
                                  'patch ',patch_id,'. Exiting ...'
                CALL s_mpi_abort()
                
            END IF
            
            
        END SUBROUTINE s_check_line_segment_patch_geometry ! -------------------
        
        
        
        !>  This subroutine verifies that the geometric parameters of
        !!      the circle patch have consistently been inputted by the
        !!      user.        
        !!  @param patch_id Patch identifier
        SUBROUTINE s_check_circle_patch_geometry(patch_id) ! -------------------

            INTEGER, INTENT(IN) :: patch_id
            
            ! Constraints on the geometric parameters of the circle patch
            IF(n == 0 .OR. p > 0 .OR. patch_icpp(patch_id)%radius <= 0d0 &
                                    .OR.                                 &
                     patch_icpp(patch_id)%x_centroid == dflt_real        &
                                    .OR.                                 &
                     patch_icpp(patch_id)%y_centroid == dflt_real) THEN
                
                PRINT '(A,I0,A)', 'Inconsistency(ies) detected in '         // &
                                  'geometric parameters of circle '         // &
                                  'patch ',patch_id,'. Exiting ...'
                CALL s_mpi_abort()
                
            END IF
            
            
        END SUBROUTINE s_check_circle_patch_geometry ! -------------------------
        
        
        
        !>  This subroutine verifies that the geometric parameters of
        !!      the rectangle patch have consistently been inputted by
        !!      the user.        
        !!  @param patch_id Patch identifier
        SUBROUTINE s_check_rectangle_patch_geometry(patch_id) ! ----------------

            INTEGER, INTENT(IN) :: patch_id
            
            
            ! Constraints on the geometric parameters of the rectangle patch
            IF(                     n == 0 .OR. p > 0                     &
                                           .OR.                           &
                      patch_icpp(patch_id)%x_centroid == dflt_real        &
                                           .OR.                           &
                      patch_icpp(patch_id)%y_centroid == dflt_real        &
                                           .OR.                           &
                           patch_icpp(patch_id)%length_x <= 0d0           &
                                           .OR.                           &
                           patch_icpp(patch_id)%length_y <= 0d0           ) THEN
                
                PRINT '(A,I0,A)', 'Inconsistency(ies) detected in '         // &
                                  'geometric parameters of rectangle '      // &
                                  'patch ',patch_id,'. Exiting ...'
                CALL s_mpi_abort()
                
            END IF
            
            
        END SUBROUTINE s_check_rectangle_patch_geometry ! ----------------------
        
        
        
        
        !> This subroutine verifies that the geometric parameters of
        !!      the line sweep patch have consistently been inputted by
        !!      the user.
        !!  @param patch_id Patch identifier
        SUBROUTINE s_check_line_sweep_patch_geometry(patch_id) ! ---------------
            
            INTEGER, INTENT(IN) :: patch_id
            
            ! Constraints on the geometric parameters of the line sweep patch
            IF(                     n == 0 .OR. p > 0                     &
                                           .OR.                           &
                      patch_icpp(patch_id)%x_centroid == dflt_real        &
                                           .OR.                           &
                      patch_icpp(patch_id)%y_centroid == dflt_real        &
                                           .OR.                           &
                     patch_icpp(patch_id)%normal(1) == dflt_real       &
                                           .OR.                           &
                     patch_icpp(patch_id)%normal(2) == dflt_real       &
                                           .OR.                           &
                     patch_icpp(patch_id)%normal(3) /= dflt_real       ) THEN
                
                PRINT '(A,I0,A)', 'Inconsistency(ies) detected in '         // &
                                  'geometric parameters of line sweep '     // &
                                  'patch ',patch_id,'. Exiting ...'
                CALL s_mpi_abort()
                
            END IF
            
            
        END SUBROUTINE s_check_line_sweep_patch_geometry ! ---------------------
        
        
        
        !>  This subroutine verifies that the geometric parameters of
        !!      the ellipse patch have consistently been inputted by
        !!      the user.
        !!  @param patch_id Patch identifier
        SUBROUTINE s_check_ellipse_patch_geometry(patch_id) ! ------------------
            
            INTEGER, INTENT(IN) :: patch_id
            
            ! Constraints on the geometric parameters of the ellipse patch
            IF(                     n == 0 .OR. p > 0                     &
                                           .OR.                           &
                      patch_icpp(patch_id)%x_centroid == dflt_real        &
                                           .OR.                           &
                      patch_icpp(patch_id)%y_centroid == dflt_real        &
                                           .OR.                           &
                     patch_icpp(patch_id)%radii(1) == dflt_real       &
                                           .OR.                           &
                     patch_icpp(patch_id)%radii(2) == dflt_real       &
                                           .OR.                           &
                     patch_icpp(patch_id)%radii(3) /= dflt_real       ) THEN
                
                PRINT '(A,I0,A)', 'Inconsistency(ies) detected in '         // &
                                  'geometric parameters of ellipse '        // &
                                  'patch ',patch_id,'. Exiting ...'
                CALL s_mpi_abort()
                
            END IF
            
            
        END SUBROUTINE s_check_ellipse_patch_geometry ! ------------------------




        !>  This subroutine verifies that the geometric parameters of
        !!      the isentropic vortex patch have been entered by the user
        !!      consistently.
        !!  @param patch_id Patch identifier
        SUBROUTINE s_check_isentropic_vortex_patch_geometry(patch_id) ! --------

            
            INTEGER, INTENT(IN) :: patch_id
            
            
            ! Constraints on the isentropic vortex patch geometric parameters
            IF(           n == 0 .OR. p > 0 .OR. model_eqns == 2          &
                                           .OR.                           &
                            patch_icpp(patch_id)%radius <= 0d0            &
                                           .OR.                           &
                           patch_icpp(patch_id)%epsilon <= 0d0            &
                                           .OR.                           &
                             patch_icpp(patch_id)%beta <= 0d0             ) THEN
                
                PRINT '(A,I0,A)', 'Inconsistency(ies) detected in '         // &
                                  'geometric parameters of isentropic '     // &
                                  'vortex patch ',patch_id,'. Exiting ...'
                CALL s_mpi_abort()
                
            END IF
            
            
        END SUBROUTINE s_check_isentropic_vortex_patch_geometry ! --------------
        
        !>  This subroutine verifies that the geometric parameters of
        !!      the analytical patch have consistently been inputted by
        !!      the user.        
        !!  @param patch_id Patch identifier
        SUBROUTINE s_check_1D_analytical_patch_geometry(patch_id) ! ---------------

            INTEGER, INTENT(IN) :: patch_id

            ! Constraints on the geometric parameters of the analytical patch
            IF(                     n > 0 .OR. p > 0                    &
                                           .OR.                         &
                        (model_eqns .ne. 4 .AND. model_eqns .ne. 2)     &
                                           .OR.                         &
                      patch_icpp(patch_id)%x_centroid == dflt_real      &
                                           .OR.                         &
                           patch_icpp(patch_id)%length_x <= 0d0           ) THEN
                
                PRINT '(A,I0,A)', 'Inconsistency(ies) detected in ' // &
                        'geometric parameters of 1D analytical ' // &
                        'patch ',patch_id,'. Exiting...'
                CALL s_mpi_abort()
            END IF
        END SUBROUTINE s_check_1D_analytical_patch_geometry ! ---------------------      
        
        !>  This subroutine verifies that the geometric parameters of
        !!      the analytical patch have consistently been inputted by
        !!      the user.        
        !!  @param patch_id Patch identifier
        SUBROUTINE s_check_2D_analytical_patch_geometry(patch_id) ! ---------------

            INTEGER, INTENT(IN) :: patch_id

            ! Constraints on the geometric parameters of the analytical patch
            IF(                     n == 0 .OR. p > 0                     &
                                           .OR.                           &
                      patch_icpp(patch_id)%x_centroid == dflt_real        &
                                           .OR.                           &
                      patch_icpp(patch_id)%y_centroid == dflt_real        &
                                           .OR.                           &
                           patch_icpp(patch_id)%length_x <= 0d0           &
                                           .OR.                           &
                           patch_icpp(patch_id)%length_y <= 0d0           ) THEN
                
                PRINT '(A,I0,A)', 'Inconsistency(ies) detected in ' // &
                        'geometric parameters of 2D analytical ' // &
                        'patch ',patch_id,'. Exiting...'
                CALL s_mpi_abort()
            END IF
        END SUBROUTINE s_check_2D_analytical_patch_geometry ! ---------------------



        !>  This subroutine verifies that the geometric parameters of
        !!      the analytical patch have consistently been inputted by
        !!      the user.
        !!  @param patch_id Patch identifier
        SUBROUTINE s_check_3D_analytical_patch_geometry(patch_id) ! ---------------

            INTEGER, INTENT(IN) :: patch_id

            ! Constraints on the geometric parameters of the analytical patch
            IF(                           p == 0                          &
                                           .OR.                           &
                      patch_icpp(patch_id)%x_centroid == dflt_real        &
                                           .OR.                           &
                      patch_icpp(patch_id)%y_centroid == dflt_real        &
                                           .OR.                           &
                      patch_icpp(patch_id)%z_centroid == dflt_real        &
                                           .OR.                           &
                           patch_icpp(patch_id)%length_x <= 0d0           &
                                           .OR.                           &
                           patch_icpp(patch_id)%length_y <= 0d0           &
                                           .OR.                           &
                           patch_icpp(patch_id)%length_z <= 0d0           ) THEN
                
                PRINT '(A,I0,A)', 'Inconsistency(ies) detected in ' // &
                        'geometric parameters of 3D analytical ' // &
                        'patch ',patch_id,'. Exiting...'
                CALL s_mpi_abort()
            END IF
        END SUBROUTINE s_check_3D_analytical_patch_geometry ! ---------------------




        !> This subroutine verifies that the geometric parameters of
        !!      the sphere patch have consistently been inputted by the
        !!      user.
        !!  @param patch_id Patch identifier
        SUBROUTINE s_check_sphere_patch_geometry(patch_id) ! -------------------
            
            INTEGER, INTENT(IN) :: patch_id
            
            ! Constraints on the geometric parameters of the sphere patch
            IF(                           p == 0                          &
                                           .OR.                           &
                           patch_icpp(patch_id)%radius <= 0d0             &
                                           .OR.                           &
                      patch_icpp(patch_id)%x_centroid == dflt_real        &
                                           .OR.                           &
                      patch_icpp(patch_id)%y_centroid == dflt_real        &
                                           .OR.                           &
                      patch_icpp(patch_id)%z_centroid == dflt_real ) THEN
                
                PRINT '(A,I0,A)', 'Inconsistency(ies) detected in '         // &
                                  'geometric parameters of sphere '         // &
                                  'patch ',patch_id,'. Exiting ...'
                CALL s_mpi_abort()
                
            END IF
            
            
        END SUBROUTINE s_check_sphere_patch_geometry ! -------------------------
        
        
        
        !>  This subroutine verifies that the geometric parameters of
        !!      the spherical harmonic  patch have consistently been
        !!      inputted by the user.        
        !!  @param patch_id Patch identifier
        SUBROUTINE s_check_spherical_harmonic_patch_geometry(patch_id) ! -------
            
            INTEGER, INTENT(IN) :: patch_id
            
            
            ! Constraints on the geometric parameters of the spherical harmonic patch
            IF(                           p == 0                               &
                                           .OR.                                &
                              patch_icpp(patch_id)%radius <= 0d0               &
                                           .OR.                                &
                      patch_icpp(patch_id)%x_centroid == dflt_real             &
                                           .OR.                                &
                      patch_icpp(patch_id)%y_centroid == dflt_real             &
                                           .OR.                                &
                      patch_icpp(patch_id)%z_centroid == dflt_real             &
                                           .OR.                                &
                ALL(patch_icpp(patch_id)%epsilon /= (/1d0,2d0,3d0,4d0,5d0/))   &
                                           .OR.                                &
                             patch_icpp(patch_id)%beta < 0d0                   &
                                           .OR.                                &
                  patch_icpp(patch_id)%beta > patch_icpp(patch_id)%epsilon) THEN
                
                PRINT '(A,I0,A)', 'Inconsistency(ies) detected in '         // &
                                  'geometric parameters of spherical '      // &
                                  'harmonic patch ',patch_id,'. Exiting ...'
                CALL s_mpi_abort()
                
            END IF
            
            
        END SUBROUTINE s_check_spherical_harmonic_patch_geometry ! -------------
        
        
        
        
        !>  This subroutine verifies that the geometric parameters of
        !!      the cuboid patch have consistently been inputted by the
        !!      user.        
        !!  @param patch_id Patch identifier
        SUBROUTINE s_check_cuboid_patch_geometry(patch_id) ! -------------------
            
            
            ! Patch identifier
            INTEGER, INTENT(IN) :: patch_id
            
            
            ! Constraints on the geometric parameters of the cuboid patch
            IF(                           p == 0                          &
                                           .OR.                           &
                      patch_icpp(patch_id)%x_centroid == dflt_real        &
                                           .OR.                           &
                      patch_icpp(patch_id)%y_centroid == dflt_real        &
                                           .OR.                           &
                      patch_icpp(patch_id)%z_centroid == dflt_real        &
                                           .OR.                           &
                           patch_icpp(patch_id)%length_x <= 0d0           &
                                           .OR.                           &
                           patch_icpp(patch_id)%length_y <= 0d0           &
                                           .OR.                           &
                           patch_icpp(patch_id)%length_z <= 0d0           ) THEN
                
                PRINT '(A,I0,A)', 'Inconsistency(ies) detected in '         // &
                                  'geometric parameters of cuboid '         // &
                                  'patch ',patch_id,'. Exiting ...'
                CALL s_mpi_abort()
                
            END IF
            
            
        END SUBROUTINE s_check_cuboid_patch_geometry ! -------------------------
        
        
        
        
        !>  This subroutine verifies that the geometric parameters of
        !!      the cylinder patch have consistently been inputted by the
        !!      user.        
        !!  @param patch_id Patch identifier
        SUBROUTINE s_check_cylinder_patch_geometry(patch_id) ! -----------------

            
            
            ! Patch identifier
            INTEGER, INTENT(IN) :: patch_id
            
            
            ! Constraints on the geometric parameters of the cylinder patch
            IF(                           p == 0                          &
                                           .OR.                           &
                      patch_icpp(patch_id)%x_centroid == dflt_real        &
                                           .OR.                           &
                      patch_icpp(patch_id)%y_centroid == dflt_real        &
                                           .OR.                           &
                      patch_icpp(patch_id)%z_centroid == dflt_real        &
                                           .OR.                           &
               (       patch_icpp(patch_id)%length_x <= 0d0 .AND.         &
                       patch_icpp(patch_id)%length_y <= 0d0 .AND.         &
                       patch_icpp(patch_id)%length_z <= 0d0             ) &
                                           .OR.                           &
               (   patch_icpp(patch_id)%length_x  >      0d0     .AND.    &
                 ( patch_icpp(patch_id)%length_y /= dflt_real .OR.     &
                   patch_icpp(patch_id)%length_z /= dflt_real      ) ) &
                                           .OR.                           &
               (   patch_icpp(patch_id)%length_y  >      0d0     .AND.    &
                 ( patch_icpp(patch_id)%length_x /= dflt_real .OR.     &
                   patch_icpp(patch_id)%length_z /= dflt_real      ) ) &
                                           .OR.                           &
               (   patch_icpp(patch_id)%length_z  >      0d0     .AND.    &
                 ( patch_icpp(patch_id)%length_x /= dflt_real .OR.     &
                   patch_icpp(patch_id)%length_y /= dflt_real      ) ) &       
                                           .OR.                           &
                            patch_icpp(patch_id)%radius <= 0d0            ) THEN
                
                PRINT '(A,I0,A)', 'Inconsistency(ies) detected in '         // &
                                  'geometric parameters of cylinder '       // &
                                  'patch ',patch_id,'. Exiting ...'
                CALL s_mpi_abort()
                
            END IF
            
            
        END SUBROUTINE s_check_cylinder_patch_geometry ! -----------------------
        
        
        
        
        !>  This subroutine verifies that the geometric parameters of
        !!      the plane sweep patch have consistently been inputted by
        !!      the user.        
        !!  @param patch_id Patch identifier
        SUBROUTINE s_check_plane_sweep_patch_geometry(patch_id) ! --------------
            
            
            ! Patch identifier
            INTEGER, INTENT(IN) :: patch_id
            
            
            ! Constraints on the geometric parameters of the plane sweep patch
            IF(                           p == 0                          &
                                           .OR.                           &
                      patch_icpp(patch_id)%x_centroid == dflt_real        &
                                           .OR.                           &
                      patch_icpp(patch_id)%y_centroid == dflt_real        &
                                           .OR.                           &
                      patch_icpp(patch_id)%z_centroid == dflt_real        &
                                           .OR.                           &
                     patch_icpp(patch_id)%normal(1) == dflt_real       &
                                           .OR.                           &
                     patch_icpp(patch_id)%normal(2) == dflt_real       &
                                           .OR.                           &
                     patch_icpp(patch_id)%normal(3) == dflt_real       ) THEN
                
                PRINT '(A,I0,A)', 'Inconsistency(ies) detected in '         // &
                                  'geometric parameters of plane sweep '    // &
                                  'patch ',patch_id,'. Exiting ...'
                CALL s_mpi_abort()
                
            END IF
            
            
        END SUBROUTINE s_check_plane_sweep_patch_geometry ! --------------------
        
        
        
        !> This subroutine verifies that the geometric parameters of
        !!      the ellipsoid patch have consistently been inputted by
        !!      the user.
        !!  @param patch_id Patch identifier
        SUBROUTINE s_check_ellipsoid_patch_geometry(patch_id) ! ----------------
           
            INTEGER, INTENT(IN) :: patch_id
            
            
            ! Constraints on the geometric parameters of the ellipsoid patch
            IF(                           p == 0                          &
                                           .OR.                           &
                      patch_icpp(patch_id)%x_centroid == dflt_real        &
                                           .OR.                           &
                      patch_icpp(patch_id)%y_centroid == dflt_real        &
                                           .OR.                           &
                      patch_icpp(patch_id)%z_centroid == dflt_real        &
                                           .OR.                           &
                     patch_icpp(patch_id)%radii(1) == dflt_real       &
                                           .OR.                           &
                     patch_icpp(patch_id)%radii(2) == dflt_real       &
                                           .OR.                           &
                     patch_icpp(patch_id)%radii(3) == dflt_real       ) THEN
                
                PRINT '(A,I0,A)', 'Inconsistency(ies) detected in '         // &
                                  'geometric parameters of ellipsoid '    // &
                                  'patch ',patch_id,'. Exiting ...'
                CALL s_mpi_abort()
                
            END IF
            
            
        END SUBROUTINE s_check_ellipsoid_patch_geometry ! ----------------------



        !>  This subroutine verifies that the geometric parameters of
        !!      the inactive patch remain unaltered by the user inputs.
        !!  @param patch_id Patch identifier
        SUBROUTINE s_check_inactive_patch_geometry(patch_id) ! -----------------
            
            INTEGER, INTENT(IN) :: patch_id
            
            
            ! Constraints on the geometric parameters of the inactive patch
            IF(      patch_icpp(patch_id)%x_centroid /= dflt_real      &
                                           .OR.                           &
                     patch_icpp(patch_id)%y_centroid /= dflt_real      &
                                           .OR.                           &
                     patch_icpp(patch_id)%z_centroid /= dflt_real      &
                                           .OR.                           &
                      patch_icpp(patch_id)%length_x /= dflt_real       &
                                           .OR.                           &
                      patch_icpp(patch_id)%length_y /= dflt_real       &
                                           .OR.                           &
                      patch_icpp(patch_id)%length_z /= dflt_real       &
                                           .OR.                           &
                       patch_icpp(patch_id)%radius /= dflt_real        &
                                           .OR.                           &
                       patch_icpp(patch_id)%epsilon /= dflt_real       &
                                           .OR.                           &
                        patch_icpp(patch_id)%beta /= dflt_real         &
                                           .OR.                           &
                     patch_icpp(patch_id)%normal(1) /= dflt_real       &
                                           .OR.                           &
                     patch_icpp(patch_id)%normal(2) /= dflt_real       &
                                           .OR.                           &
                     patch_icpp(patch_id)%normal(3) /= dflt_real       &
                                           .OR.                           &
                     patch_icpp(patch_id)%radii(1) /= dflt_real        &
                                           .OR.                           &
                     patch_icpp(patch_id)%radii(2) /= dflt_real        &
                                           .OR.                           &
                     patch_icpp(patch_id)%radii(3) /= dflt_real        ) THEN
                
                PRINT '(A,I0,A)', 'Inconsistency(ies) detected in '         // &
                                  'geometric parameters of inactive '       // &
                                  'patch ',patch_id,'. Exiting ...'
                CALL s_mpi_abort()
                
            END IF
            
            
        END SUBROUTINE s_check_inactive_patch_geometry ! -----------------------
        
        
        !>  This subroutine verifies that any rights granted to the
        !!      given active patch, to overwrite the preceding active
        !!      patches, were consistently inputted by the user.       
        !!  @param patch_id Patch identifier
        SUBROUTINE s_check_active_patch_alteration_rights(patch_id) ! ----------
            
            INTEGER, INTENT(IN) :: patch_id
            
            
            ! Constraints on the alteration rights of an active patch
            IF(      patch_icpp(patch_id)%alter_patch(0) .EQV. .FALSE.    &
                                           .OR.                           &
                     ANY(patch_icpp(patch_id)%alter_patch(patch_id:))     ) THEN
                
                PRINT '(A,I0,A)', 'Inconsistency(ies) detected in '         // &
                                  'alteration rights of active '            // &
                                  'patch ',patch_id,'. Exiting ...'
                CALL s_mpi_abort()
                
            END IF
            
            
        END SUBROUTINE s_check_active_patch_alteration_rights ! ----------------
        
        
        
        !>  This subroutine verifies that the rights of the given
        !!      inactive patch, to overwrite the preceding patches,
        !!      remain unaltered by the user inputs.        
        !!  @param patch_id Patch identifier
        SUBROUTINE s_check_inactive_patch_alteration_rights(patch_id) ! --------

            
            
            ! Patch identifier
            INTEGER, INTENT(IN) :: patch_id
            
            
            ! Constraints on the alteration rights of an inactive patch
            IF(      patch_icpp(patch_id)%alter_patch(0) .EQV. .FALSE.    &
                                           .OR.                           &
                         ANY(patch_icpp(patch_id)%alter_patch(1:))        ) THEN
                
                PRINT '(A,I0,A)', 'Inconsistency(ies) detected in '         // &
                                  'alteration rights of inactive '          // &
                                  'patch ',patch_id,'. Exiting ...'
                CALL s_mpi_abort()
                
            END IF
            
            
        END SUBROUTINE s_check_inactive_patch_alteration_rights ! --------------
        
        
        !> This subroutine verifies that the smoothing parameters of
        !!      the given patch, which supports the smoothing out of its
        !!      boundaries, have consistently been inputted by the user.
        !!  @param patch_id Patch identifier
        SUBROUTINE s_check_supported_patch_smoothing(patch_id) ! ---------------
            
            INTEGER, INTENT(IN) :: patch_id
            
            
            ! Constraints on the smoothing parameters of a supported patch
            IF( (              patch_icpp(patch_id)%smoothen              &
                                           .AND.                          &
                  (  patch_icpp(patch_id)%smooth_patch_id >= patch_id     &
                                           .OR.                           &
                        patch_icpp(patch_id)%smooth_patch_id == 0         &
                                           .OR.                           &
                         patch_icpp(patch_id)%smooth_coeff <= 0d0     ) ) &
                                           .OR.                           &
                (        (patch_icpp(patch_id)%smoothen .NEQV. .TRUE.)    &
                                           .AND.                          &
                  ( patch_icpp(patch_id)%smooth_patch_id /= patch_id      &
                                           .OR.                           &
                    patch_icpp(patch_id)%smooth_coeff /= dflt_real ) ) ) THEN
                
                PRINT '(A,I0,A)', 'Inconsistency(ies) detected in '         // &
                                  'smoothing parameters of supported '      // &
                                  'patch ',patch_id,'. Exiting ...'
                CALL s_mpi_abort()
                
            END IF
            
            
        END SUBROUTINE s_check_supported_patch_smoothing ! ---------------------
        
        
        
        !> This subroutine verifies that the smoothing parameters of
        !!      the given patch, which does not support the smoothing out
        !!          of its boundaries, remain unaltered by the user inputs.
        !!  @param patch_id Patch identifier
        SUBROUTINE s_check_unsupported_patch_smoothing(patch_id) ! -------------
            
            
            ! Patch identifier
            INTEGER, INTENT(IN) :: patch_id
            
            
            ! Constraints on the smoothing parameters of an unsupported patch
            IF(                patch_icpp(patch_id)%smoothen              &
                                           .OR.                           &
                    patch_icpp(patch_id)%smooth_patch_id /= patch_id      &
                                           .OR.                           &
                    patch_icpp(patch_id)%smooth_coeff /= dflt_real     ) THEN
                
                PRINT '(A,I0,A)', 'Inconsistency(ies) detected in '         // &
                                  'smoothing parameters of unsupported '    // &
                                  'patch ',patch_id,'. Exiting ...'
                CALL s_mpi_abort()
                
            END IF
            
            
        END SUBROUTINE s_check_unsupported_patch_smoothing ! -------------------
        
        
        
        !>  This subroutine verifies that the primitive variables
        !!      associated with the given active patch are physically
        !!      consistent.        
        !!  @param patch_id Patch identifier
        SUBROUTINE s_check_active_patch_primitive_variables(patch_id) ! --------
            
            INTEGER, INTENT(IN) :: patch_id
            
            
            ! Constraints on the primitive variables of an active patch
            IF(         patch_icpp(patch_id)%vel(1) == dflt_real       &
                                           .OR.                           &
               (n == 0 .AND. patch_icpp(patch_id)%vel(2) /= dflt_real) &
                                           .OR.                           &
               (n  > 0 .AND. patch_icpp(patch_id)%vel(2) == dflt_real) &
                                           .OR.                           &
               (p == 0 .AND. patch_icpp(patch_id)%vel(3) /= dflt_real) &
                                           .OR.                           &
               (p  > 0 .AND. patch_icpp(patch_id)%vel(3) == dflt_real) &
            !                               .OR.                           &
            !                 patch_icpp(patch_id)%pres <= 0d0             &
                                           .OR.                           &
               (                      model_eqns == 1          .AND.      &
                         (  patch_icpp(patch_id)%rho    <= 0d0 .OR.       &
                            patch_icpp(patch_id)%gamma  <= 0d0 .OR.       &
                            patch_icpp(patch_id)%pi_inf <  0d0       )  ) &
                                           .OR.                           &
               (            patch_icpp(patch_id)%geometry == 5            &
                                           .AND.                          &
                              patch_icpp(patch_id)%pi_inf > 0           ) &
                                           .OR.                           &
               (                      model_eqns == 2                     &
                                           .AND.                          &
                (ANY(patch_icpp(patch_id)%alpha_rho(1:num_fluids) < 0d0)  &
                  !                         .OR.                           &
                  !ANY(patch_icpp(patch_id)%alpha(1: adv_idx%end           &
                  !                                 -   E_idx    ) < 0d0)  &
                  !                         .OR.                           &
                  !SUM(patch_icpp(patch_id)%alpha(1: adv_idx%end           &
                  !                                 -   E_idx    ))> 1d0)  
                  ))) THEN
                
                PRINT '(A,I0,A)', 'Inconsistency(ies) detected in '         // &
                                  'primitive variables of active '          // &
                                  'patch ',patch_id,'. Exiting ...'
                CALL s_mpi_abort()
                
            END IF
            
            IF(model_eqns == 2 .AND. num_fluids < num_fluids_max) THEN
                
                IF(  ANY( patch_icpp(patch_id)%alpha_rho(num_fluids+1:)   &
                                        /= dflt_real                    ) &
                                           .OR.                           &
                     ANY(   patch_icpp(patch_id)%alpha(num_fluids+1:)     &
                                        /= dflt_real                    ) &
                                           .OR.                           &
                    (           (adv_alphan .NEQV. .TRUE.)                &
                                           .AND.                          &
                     patch_icpp(patch_id)%alpha(num_fluids) /= dflt_real) &
                                           .OR.                           &
                    (                   adv_alphan                        &
                                           .AND.                          &
                     patch_icpp(patch_id)%alpha(num_fluids) == dflt_real) ) THEN
                    
                    PRINT '(A,I0,A)', 'Inconsistency(ies) detected in '     // &
                                      'primitive variables of active '      // &
                                      'patch ',patch_id,'. Exiting ...'
                    !CALL s_mpi_abort()
                    
                END IF
                
            END IF
            
            
        END SUBROUTINE s_check_active_patch_primitive_variables ! --------------
        
        
        
        
        !>  This subroutine verifies that the primitive variables
        !!      associated with the given inactive patch remain unaltered
        !!      by the user inputs.        
        !!  @param patch_id Patch identifier
        SUBROUTINE s_check_inactive_patch_primitive_variables(patch_id) ! ------
          
            INTEGER, INTENT(IN) :: patch_id
            
            
            ! Constraints on the primitive variables of an inactive patch
            IF(    ANY(patch_icpp(patch_id)%alpha_rho /= dflt_real)    &
                                           .OR.                           &
                        patch_icpp(patch_id)%rho /= dflt_real          &
                                           .OR.                           &
                      ANY(patch_icpp(patch_id)%vel /= dflt_real)       &
                                           .OR.                           &
                        patch_icpp(patch_id)%pres /= dflt_real         &
                                           .OR.                           &
                     ANY(patch_icpp(patch_id)%alpha /= dflt_real)      &
                                           .OR.                           &
                        patch_icpp(patch_id)%gamma /= dflt_real        &
                                           .OR.                           &
                       patch_icpp(patch_id)%pi_inf /= dflt_real        ) THEN
                
                PRINT '(A,I0,A)', 'Inconsistency(ies) detected in '         // &
                                  'primitive variables of inactive '        // &
                                  'patch ',patch_id,'. Exiting ...'
                CALL s_mpi_abort()
                
            END IF
            
            
        END SUBROUTINE s_check_inactive_patch_primitive_variables ! ------------
        
        
        
        !> The goal of this subroutine is to read in any preexisting
        !!      grid data as well as based on the imported grid, complete
        !!      the necessary global computational domain parameters.        
        !! @param dflt_int Default null integer
        SUBROUTINE s_read_serial_grid_data_files(dflt_int) ! ---

            
        INTEGER, INTENT(IN) :: dflt_int
            
            ! Generic string used to store the address of a particular file
            CHARACTER(LEN = LEN_TRIM(case_dir) + 3*name_len) :: file_loc
            
            ! Logical variable used to test the existence of folders
            LOGICAL :: dir_check
            
            ! Generic logical used for the purpose of asserting whether a file
            ! is or is not present in the designated location
            LOGICAL :: file_check
            
            ! Setting address of the local processor rank and time-step directory
            WRITE(proc_rank_dir, '(A,I0)') '/p', proc_rank
            proc_rank_dir = TRIM(case_dir) // TRIM(proc_rank_dir)
            
            WRITE(t_step_dir, '(A,I0)') '/', t_step_old
            t_step_dir = TRIM(proc_rank_dir) // TRIM(t_step_dir)
            
            
            ! Inquiring as to the existence of the time-step directory
            file_loc = TRIM(t_step_dir) // '/.'
            CALL my_inquire(file_loc,dir_check)
            
            ! If the time-step directory is missing, the pre-process exits
            IF(dir_check .NEQV. .TRUE.) THEN
                PRINT '(A)', 'Time-step folder ' // TRIM(t_step_dir) // &
                             ' is missing. Exiting ...'
                CALL s_mpi_abort()
            END IF
            
            
            ! Reading the Grid Data File for the x-direction ===================
            
            ! Checking whether x_cb.dat exists
            file_loc = TRIM(t_step_dir) // '/x_cb.dat'
            INQUIRE(FILE = TRIM(file_loc), EXIST = file_check)
            
            ! If it exists, x_cb.dat is read
            IF(file_check) THEN
                OPEN(1, FILE = TRIM(file_loc), FORM = 'unformatted', &
                        STATUS = 'old', ACTION = 'read')
                READ(1) x_cb
                CLOSE(1)
            ELSE
                PRINT '(A)', 'File x_cb.dat is missing in ' // &
                             TRIM(t_step_dir) // '. Exiting ...'
                CALL s_mpi_abort()
            END IF
            
            ! Computing cell-center locations
            x_cc = (x_cb(0:m) + x_cb(-1:(m-1)))/2d0
            
            ! Computing minimum cell-width
            dx = MINVAL(x_cb(0:m) - x_cb(-1:m-1))
            IF(num_procs > 1) CALL s_mpi_reduce_min(dx)
            
            ! Setting locations of domain bounds
            x_domain%beg = x_cb(-1)
            x_domain%end = x_cb( m)
            
            ! ==================================================================
            
            
            ! Reading the Grid Data File for the y-direction ===================
            
            IF(n > 0) THEN
                
                ! Checking whether y_cb.dat exists
                file_loc = TRIM(t_step_dir) // '/y_cb.dat'
                INQUIRE(FILE = TRIM(file_loc), EXIST = file_check)
                
                ! If it exists, y_cb.dat is read
                IF(file_check) THEN
                    OPEN(1, FILE = TRIM(file_loc), FORM = 'unformatted', &
                            STATUS = 'old', ACTION = 'read')
                    READ(1) y_cb
                    CLOSE(1)
                ELSE
                    PRINT '(A)', 'File y_cb.dat is missing in ' // &
                                  TRIM(t_step_dir) // '. Exiting ...'
                    CALL s_mpi_abort()
                END IF
                
                ! Computing cell-center locations
                y_cc = (y_cb(0:n) + y_cb(-1:(n-1)))/2d0
                
                ! Computing minimum cell-width
                dy = MINVAL(y_cb(0:n) - y_cb(-1:n-1))
                IF(num_procs > 1) CALL s_mpi_reduce_min(dy)
                
                ! Setting locations of domain bounds
                y_domain%beg = y_cb(-1)
                y_domain%end = y_cb( n)
                
            ! ==================================================================
                
                
            ! Reading the Grid Data File for the z-direction ===================
                
                IF(p > 0) THEN
                    
                    ! Checking whether z_cb.dat exists
                    file_loc = TRIM(t_step_dir) // '/z_cb.dat'
                    INQUIRE(FILE = TRIM(file_loc), EXIST = file_check)
                    
                    ! If it exists, z_cb.dat is read
                    IF(file_check) THEN
                        OPEN(1, FILE = TRIM(file_loc), FORM = 'unformatted', &
                                STATUS = 'old', ACTION = 'read')
                        READ(1) z_cb
                        CLOSE(1)
                    ELSE
                        PRINT '(A)', 'File z_cb.dat is missing in ' // &
                                     TRIM(t_step_dir) // '. Exiting ...'
                        CALL s_mpi_abort()
                    END IF
                    
                    ! Computing cell-center locations
                    z_cc = (z_cb(0:p) + z_cb(-1:(p-1)))/2d0
                    
                    ! Computing minimum cell-width
                    dz = MINVAL(z_cb(0:p) - z_cb(-1:p-1))
                    IF(num_procs > 1) CALL s_mpi_reduce_min(dz)
                    
                    ! Setting locations of domain bounds
                    z_domain%beg = z_cb(-1)
                    z_domain%end = z_cb( p)
                    
                END IF
                
            END IF
            
            ! ==================================================================
            
            
            ! If only the preexisting grid data files are read in and there will
            ! not be any preexisting initial condition data files imported, then
            ! the directory associated with the rank of the local processor may
            ! be cleaned to make room for the new pre-process data. In addition,
            ! the time-step directory that will contain the new grid and initial
            ! condition data are also generated.
            IF(old_ic .NEQV. .TRUE.) THEN
                CALL SYSTEM('rm -rf ' // TRIM(proc_rank_dir) // '/*')
                CALL SYSTEM('mkdir -p ' // TRIM(proc_rank_dir) // '/0')
            END IF
            
            
        END SUBROUTINE s_read_serial_grid_data_files ! --------------------------------
        
        
        
        !> Cell-boundary data are checked for consistency by looking
        !!      at the (non-)uniform cell-width distributions for all the
        !!      active coordinate directions and making sure that all of
        !!      the cell-widths are positively valued
        SUBROUTINE s_check_grid_data_files() ! -----------------
            
            
            ! Cell-boundary Data Consistency Check in x-direction ==============
            
            IF(ANY(x_cb(0:m) - x_cb(-1:m-1) <= 0d0)) THEN
                PRINT '(A)', 'x_cb.dat in ' // TRIM(t_step_dir) // &
                             ' contains non-positive cell-spacings. Exiting ...'
                CALL s_mpi_abort()
            END IF
            
            ! ==================================================================
            
            
            ! Cell-boundary Data Consistency Check in y-direction ==============
            
            IF(n > 0) THEN
                
                IF(ANY(y_cb(0:n) - y_cb(-1:n-1) <= 0d0)) THEN
                    PRINT '(A)', 'y_cb.dat in ' // TRIM(t_step_dir) // &
                                 ' contains non-positive cell-spacings. ' // &
                                 'Exiting ...'
                    CALL s_mpi_abort()
                END IF
                
            ! ==================================================================
                
                
            ! Cell-boundary Data Consistency Check in z-direction ==============
                
                IF(p > 0) THEN
                    
                    IF(ANY(z_cb(0:p) - z_cb(-1:p-1) <= 0d0)) THEN
                        PRINT '(A)', 'z_cb.dat in ' // TRIM(t_step_dir) // &
                                     ' contains non-positive cell-spacings' // &
                                     ' .Exiting ...'
                        CALL s_mpi_abort()
                    END IF
                    
                END IF
                
            END IF
            
            ! ==================================================================
            
            
        END SUBROUTINE s_check_grid_data_files ! -------------------------------
        
        
        
        !> The goal of this subroutine is to read in any preexisting
        !!      initial condition data files so that they may be used by
        !!      the pre-process as a starting point in the creation of an
        !!      all new initial condition.      
        !! @param q_cons_vf Conservative variables
        SUBROUTINE s_read_serial_ic_data_files(q_cons_vf) ! ---------------------------

            TYPE(scalar_field), &
            DIMENSION(sys_size), &
            INTENT(INOUT) :: q_cons_vf
            

            CHARACTER(LEN = LEN_TRIM(case_dir) + 3*name_len) :: file_loc !<
            ! Generic string used to store the address of a particular file
            
            CHARACTER(LEN = &
            INT(FLOOR(LOG10(REAL(sys_size, KIND(0d0))))) + 1) :: file_num !<
            !! Used to store the variable position, in character form, of the
            !! currently manipulated conservative variable file
            
            LOGICAL :: file_check !<
            !! Generic logical used for the purpose of asserting whether a file
            !! is or is not present in the designated location
            

            INTEGER :: i !< Generic loop iterator
            
            
            ! Reading the Conservative Variables Data Files ====================
            DO i = 1, sys_size
                
                ! Checking whether data file associated with variable position
                ! of the currently manipulated conservative variable exists
                WRITE(file_num, '(I0)') i
                file_loc = TRIM(t_step_dir) // '/q_cons_vf' // &
                           TRIM(file_num) // '.dat'
                INQUIRE(FILE = TRIM(file_loc), EXIST = file_check)
                
                ! If it exists, the data file is read
                IF(file_check) THEN
                    OPEN(1, FILE = TRIM(file_loc), FORM = 'unformatted', &
                            STATUS = 'old', ACTION = 'read')
                    READ(1) q_cons_vf(i)%sf
                    CLOSE(1)
                ELSE
                    PRINT '(A)', 'File q_cons_vf' // TRIM(file_num) // &
                                 '.dat is missing in ' // TRIM(t_step_dir) // &
                                 '. Exiting ...'
                    CALL s_mpi_abort()
                END IF
                
            END DO
            
            ! ==================================================================
            
            
            ! Since the preexisting grid and initial condition data files have
            ! been read in, the directory associated with the rank of the local
            ! process may be cleaned out to make room for new pre-process data.
            ! In addition, the time-step folder that will contain the new grid
            ! and initial condition data are also generated.
            CALL SYSTEM('rm -rf ' // TRIM(proc_rank_dir) // '/*')
            CALL SYSTEM('mkdir -p ' // TRIM(proc_rank_dir) // '/0')
            
            
        END SUBROUTINE s_read_serial_ic_data_files ! ----------------------------------
        
        
        
        
        !> Cell-boundary data are checked for consistency by looking
        !!      at the (non-)uniform cell-width distributions for all the
        !!      active coordinate directions and making sure that all of
        !!      the cell-widths are positively valued
        !! @param dflt_int Default null integer
        SUBROUTINE s_read_parallel_grid_data_files(dflt_int)

            INTEGER, INTENT(IN) :: dflt_int

            REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:) :: x_cb_glb, y_cb_glb, z_cb_glb

            INTEGER :: ifile, ierr, data_size
            INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status

            CHARACTER(LEN=path_len + 2*name_len) :: file_loc
            LOGICAL :: file_exist

            ALLOCATE(x_cb_glb(-1:m_glb))
            ALLOCATE(y_cb_glb(-1:n_glb))
            ALLOCATE(z_cb_glb(-1:p_glb))

            ! Read in cell boundary locations in x-direction
            file_loc = TRIM(case_dir) // '/restart_data' // TRIM(mpiiofs) // 'x_cb.dat'
            INQUIRE(FILE = TRIM(file_loc),EXIST = file_exist)

            IF (file_exist) THEN
                data_size = m_glb+2
                CALL MPI_FILE_OPEN(MPI_COMM_WORLD,file_loc,MPI_MODE_RDONLY,mpi_info_int,ifile,ierr)
                CALL MPI_FILE_READ_ALL(ifile,x_cb_glb,data_size,MPI_DOUBLE_PRECISION,status,ierr)
                CALL MPI_FILE_CLOSE(ifile,ierr)
            ELSE
                PRINT '(A)', 'File ', TRIM(file_loc), ' is missing. Exiting... '
                CALL s_mpi_abort()
            END IF

            ! Assigning local cell boundary locations
            x_cb(-1:m) = x_cb_glb((start_idx(1)-1):(start_idx(1)+m))
            ! Computing cell center locations
            x_cc(0:m) = (x_cb(0:m) + x_cb(-1:(m-1)))/2d0
            ! Computing minimum cell width
            dx = MINVAL(x_cb(0:m) - x_cb(-1:(m-1)))
            IF (num_procs > 1) CALL s_mpi_reduce_min(dx)
            ! Setting locations of domain bounds
            x_domain%beg = x_cb(-1)
            x_domain%end = x_cb( m)

            IF (n > 0) THEN
                ! Read in cell boundary locations in y-direction
                file_loc = TRIM(case_dir) // '/restart_data' // TRIM(mpiiofs) // 'y_cb.dat'
                INQUIRE(FILE = TRIM(file_loc),EXIST = file_exist)

                IF (file_exist) THEN
                    data_size = n_glb+2
                    CALL MPI_FILE_OPEN(MPI_COMM_WORLD,file_loc,MPI_MODE_RDONLY,mpi_info_int,ifile,ierr)
                    CALL MPI_FILE_READ_ALL(ifile,y_cb_glb,data_size,MPI_DOUBLE_PRECISION,status,ierr)
                    CALL MPI_FILE_CLOSE(ifile,ierr)
                ELSE
                    PRINT '(A)', 'File ', TRIM(file_loc), ' is missing. Exiting... '
                    CALL s_mpi_abort()
                END IF

                ! Assigning local cell boundary locations
                y_cb(-1:n) = y_cb_glb((start_idx(2)-1):(start_idx(2)+n))
                ! Computing cell center locations
                y_cc = (y_cb(0:n) + y_cb(-1:(n-1)))/2d0
                ! Computing minimum cell width
                dy = MINVAL(y_cb(0:n) - y_cb(-1:(n-1)))
                IF (num_procs > 1) CALL s_mpi_reduce_min(dy)
                ! Setting locations of domain bounds
                y_domain%beg = y_cb(-1)
                y_domain%end = y_cb( n)

                IF (p > 0) THEN
                    ! Read in cell boundary locations in z-direction
                    file_loc = TRIM(case_dir) // '/restart_data' // TRIM(mpiiofs) // 'z_cb.dat'
                    INQUIRE(FILE = TRIM(file_loc),EXIST = file_exist)

                    IF (file_exist) THEN
                        data_size = p_glb+2
                        CALL MPI_FILE_OPEN(MPI_COMM_WORLD,file_loc,MPI_MODE_RDONLY,mpi_info_int,ifile,ierr)
                        CALL MPI_FILE_READ_ALL(ifile,z_cb_glb,data_size,MPI_DOUBLE_PRECISION,status,ierr)
                        CALL MPI_FILE_CLOSE(ifile,ierr)
                    ELSE
                        PRINT '(A)', 'File ', TRIM(file_loc), ' is missing. Exiting... '
                        CALL s_mpi_abort()
                    END IF

                    ! Assigning local cell boundary locations
                    z_cb(-1:p) = z_cb_glb((start_idx(3)-1):(start_idx(3)+p))
                    ! Computing cell center locations
                    z_cc = (z_cb(0:p) + z_cb(-1:(p-1)))/2d0
                    ! Computing minimum cell width
                    dz = MINVAL(z_cb(0:p) - z_cb(-1:(p-1)))
                    IF (num_procs > 1) CALL s_mpi_reduce_min(dz)
                    ! Setting locations of domain bounds
                    z_domain%beg = z_cb(-1)
                    z_domain%end = z_cb( p)

                END IF
            END IF

            DEALLOCATE(x_cb_glb, y_cb_glb, z_cb_glb)


        END SUBROUTINE s_read_parallel_grid_data_files ! -----------------------



        !> The goal of this subroutine is to read in any preexisting
        !!      initial condition data files so that they may be used by
        !!      the pre-process as a starting point in the creation of an
        !!      all new initial condition.      
        !! @param q_cons_vf Conservative variables
        SUBROUTINE s_read_parallel_ic_data_files(q_cons_vf) ! ------------------

            TYPE(scalar_field), &
            DIMENSION(sys_size), &
            INTENT(INOUT) :: q_cons_vf

            INTEGER :: ifile, ierr, data_size
            INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
            INTEGER(KIND=MPI_OFFSET_KIND) :: disp
            INTEGER(KIND=MPI_OFFSET_KIND) :: m_MOK, n_MOK, p_MOK
            INTEGER(KIND=MPI_OFFSET_KIND) :: WP_MOK, var_MOK, str_MOK
            INTEGER(KIND=MPI_OFFSET_KIND) :: NVARS_MOK
            INTEGER(KIND=MPI_OFFSET_KIND) :: MOK

            CHARACTER(LEN=path_len + 2*name_len) :: file_loc
            LOGICAL :: file_exist

            INTEGER :: i

            ! Open the file to read
            WRITE(file_loc, '(I0,A)') t_step_old, '.dat'
            file_loc = TRIM(restart_dir) // TRIM(mpiiofs) // TRIM(file_loc)
            INQUIRE(FILE = TRIM(file_loc),EXIST = file_exist)

            IF (file_exist) THEN
                CALL MPI_FILE_OPEN(MPI_COMM_WORLD,file_loc,MPI_MODE_RDONLY,mpi_info_int,ifile,ierr)

                ! Initialize MPI data I/O
                CALL s_initialize_mpi_data(q_cons_vf)

                ! Size of local arrays
                data_size = (m+1)*(n+1)*(p+1)

                ! Resize some integers so MPI can read even the biggest files
                m_MOK     = INT(m_glb+1,  MPI_OFFSET_KIND)
                n_MOK     = INT(n_glb+1,  MPI_OFFSET_KIND)
                p_MOK     = INT(p_glb+1,  MPI_OFFSET_KIND)
                WP_MOK    = INT(8d0,      MPI_OFFSET_KIND)
                MOK       = INT(1d0,      MPI_OFFSET_KIND)
                str_MOK   = INT(name_len, MPI_OFFSET_KIND)
                NVARS_MOK = INT(sys_size, MPI_OFFSET_KIND)

                ! Read the data for each variable
                DO i = 1, adv_idx%end
                    var_MOK = INT(i, MPI_OFFSET_KIND)

                    ! Initial displacement to skip at beginning of file
                    disp = m_MOK*MAX(MOK,n_MOK)*MAX(MOK,p_MOK)*WP_MOK*(var_MOK-1)

                    CALL MPI_FILE_SET_VIEW(ifile,disp,MPI_DOUBLE_PRECISION,MPI_IO_DATA%view(i), &
                                'native',mpi_info_int,ierr)
                    CALL MPI_FILE_READ(ifile,MPI_IO_DATA%var(i)%sf,data_size, &
                                MPI_DOUBLE_PRECISION,status,ierr)
                END DO

                CALL s_mpi_barrier()

                CALL MPI_FILE_CLOSE(ifile,ierr)

            ELSE
                PRINT '(A)', 'File ', TRIM(file_loc), ' is missing. Exiting... '
                CALL s_mpi_abort()
            END IF
            CALL s_mpi_barrier()
            IF (proc_rank == 0) CALL SYSTEM('rm -f ' // TRIM(file_loc))

        END SUBROUTINE s_read_parallel_ic_data_files ! -------------------------





        SUBROUTINE s_initialize_start_up_module() !-----------------------------

            IF (parallel_io .NEQV. .TRUE.) THEN
                s_read_grid_data_files => s_read_serial_grid_data_files
                s_read_ic_data_files => s_read_serial_ic_data_files
            ELSE
                s_read_grid_data_files => s_read_parallel_grid_data_files
                s_read_ic_data_files => s_read_parallel_ic_data_files
            END IF

        END SUBROUTINE s_initialize_start_up_module ! --------------------------





        SUBROUTINE s_finalize_start_up_module() ! ------------------------------

            s_read_grid_data_files => NULL()
            s_read_ic_data_files => NULL()

        END SUBROUTINE s_finalize_start_up_module ! ----------------------------

END MODULE m_start_up
