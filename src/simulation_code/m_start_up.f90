!>
!! @file m_start_up.f90
!! @brief Contains module m_start_up
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief The purpose of the module is primarily to read in the files that
!!              contain the inputs, the initial condition data and the grid data
!!              that are provided by the user. The module is additionally tasked
!!              with verifying the consistency of the user inputs and completing
!!              the grid variablesThe purpose of the module is primarily to read 
!!              in the files that
!!              contain the inputs, the initial condition data and the grid data
!!              that are provided by the user. The module is additionally tasked
!!              with verifying the consistency of the user inputs and completing
!!              the grid variables.
MODULE m_start_up
    
    
    ! Dependencies =============================================================
    USE m_derived_types        !< Definitions of the derived types
    
    USE m_global_parameters    !< Definitions of the global parameters
    
    USE m_mpi_proxy            !< Message passing interface (MPI) module proxy
    
    USE m_variables_conversion !< State variables type conversion procedures

    USE m_compile_specific
    ! ==========================================================================
    
    
    IMPLICIT NONE
    
    PRIVATE; PUBLIC :: s_initialize_start_up_module, &
                       s_read_input_file, &
                       s_check_input_file, &
                       s_read_data_files, &
                       s_read_serial_data_files, &
                       s_read_parallel_data_files, &
                       s_populate_grid_variables_buffers, &
                       s_account_for_capillary_potential_energy, &
                       s_initialize_internal_energy_equations, &
                       s_finalize_start_up_module 

    ABSTRACT INTERFACE ! ===================================================

        !! @param q_cons_vf  Conservative variables
        SUBROUTINE s_read_abstract_data_files(q_cons_vf) ! -----------

            IMPORT :: scalar_field, sys_size

            TYPE(scalar_field), &
            DIMENSION(sys_size), &
            INTENT(INOUT) :: q_cons_vf

        END SUBROUTINE s_read_abstract_data_files ! -----------------

    END INTERFACE ! ========================================================

    TYPE(scalar_field), ALLOCATABLE, DIMENSION(:)  :: grad_x_vf,grad_y_vf,grad_z_vf,norm_vf,kappa_vf

    PROCEDURE(s_read_abstract_data_files), POINTER :: s_read_data_files => NULL()
    
    CONTAINS
        
        
        
        !>  The purpose of this procedure is to first verify that an
        !!      input file has been made available by the user. Provided
        !!      that this is so, the input file is then read in.
        SUBROUTINE s_read_input_file() ! ---------------------------------------
            
            
            ! Relative path to the input file provided by the user
            CHARACTER(LEN = name_len) :: file_path = './simulation.inp'
            

            LOGICAL :: file_exist !<
            !! Logical used to check the existence of the input file
            
            ! Namelist of the global parameters which may be specified by user
            NAMELIST /user_inputs/ case_dir, run_time_info, m, n, p, dt,     &
                                   t_step_start, t_step_stop, t_step_save,   &
                                   model_eqns, num_fluids, adv_alphan,       &
                                   mpp_lim, time_stepper, weno_vars,         &
                                   weno_order, weno_eps, char_decomp,        &
                                   mapped_weno, mp_weno, weno_avg,           &
                                   riemann_solver, wave_speeds, avg_state,   &
                                   commute_err, split_err, bc_x, bc_y, bc_z, &
                                   hypoelasticity,                           &
                                   fluid_pp, com_wrt, cb_wrt, probe_wrt,     &
                                   fd_order, probe, num_probes, t_step_old,  &
                                   threshold_mf, moment_order, alt_crv,      &
                                   alt_soundspeed, mixture_err, weno_Re_flux,&
                                   tvd_riemann_flux, tvd_rhs_flux,           &
                                   tvd_wave_speeds, flux_lim, We_rhs_flux,   &
                                   We_riemann_flux, We_src, null_weights,    &
                                   We_wave_speeds, lsq_deriv, precision,     & 
                                   parallel_io,                              &
                                   regularization, reg_eps, cyl_coord,       & 
                                   rhoref, pref, bubbles, bubble_model,      &
                                   R0ref, nb, Ca, Web, Re_inv,               &
                                   monopole, mono, num_mono,                 &
                                   polytropic, thermal,                      &
                                   integral, integral_wrt, num_integrals,    &
                                   polydisperse, poly_sigma, qbmm, nnode,    &
                                   R0_type, DEBUG, t_tol
            
            
            ! Checking that an input file has been provided by the user. If it
            ! has, then the input file is read in, otherwise, simulation exits.
            INQUIRE(FILE = TRIM(file_path), EXIST = file_exist)
            
            IF(file_exist) THEN
                OPEN(1, FILE   = TRIM(file_path), &
                        FORM   = 'formatted'    , &
                        ACTION = 'read'         , &
                        STATUS = 'old'            )
                READ(1, NML = user_inputs); CLOSE(1)

                ! Store BC information into global BC
                bc_x_glb%beg = bc_x%beg
                bc_x_glb%end = bc_x%end
                bc_y_glb%beg = bc_y%beg
                bc_y_glb%end = bc_y%end
                bc_z_glb%beg = bc_z%beg
                bc_z_glb%end = bc_z%end

                ! Store m,n,p into global m,n,p
                m_glb = m
                n_glb = n
                p_glb = p

            ELSE
                PRINT '(A)', TRIM(file_path) // ' is missing. Exiting ...'
                CALL s_mpi_abort()
            END IF
            
            
        END SUBROUTINE s_read_input_file ! -------------------------------------
        
        !> The goal of this procedure is to verify that each of the
        !!      user provided inputs is valid and that their combination
        !!      consitutes a meaningful configuration for the simulation.
        SUBROUTINE s_check_input_file() ! --------------------------------------
            
            
            ! Relative path to the current directory file in the case directory
            CHARACTER(LEN = path_len) :: file_path
            
            ! Logical used to check the existence of the current directory file
            LOGICAL :: file_exist
            
            ! Generic loop iterators
            INTEGER :: i,j

            INTEGER :: bub_fac !for allowing an extra fluid_pp if there are bubbles

            bub_fac = 0
            IF ( bubbles .AND. (num_fluids == 1) ) bub_fac = 1
            
            ! Logistics ========================================================
            file_path = TRIM(case_dir) // '/.'
            
            CALL my_inquire(file_path,file_exist)

            IF(file_exist .NEQV. .TRUE.) THEN
                PRINT '(A)', TRIM(file_path) // ' is missing. Exiting ...'
                CALL s_mpi_abort()
            END IF
            ! ==================================================================
            
           
            ! Computational Domain Parameters ==================================
            IF(m <= 0) THEN
                PRINT '(A)', 'Unsupported value of m. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(n < 0) THEN
                PRINT '(A)', 'Unsupported value of n. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(p < 0) THEN
                PRINT '(A)', 'Unsupported value of p. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(cyl_coord .AND. p > 0 .AND. MOD(p,2) /= 1) THEN
                PRINT '(A)', 'Unsupported value of p. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(n == 0 .AND. p > 0) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'n and p. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(dt <= 0) THEN
                PRINT '(A)', 'Unsupported value of dt. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(t_step_start < 0) THEN
                PRINT '(A)', 'Unsupported value of t_step_start. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(t_step_stop <= t_step_start) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             't_step_start and t_step_stop. '        // &
                             'Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(t_step_save > t_step_stop - t_step_start) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             't_step_start, t_step_stop and '        // &
                             't_step_save. Exiting ...'
                CALL s_mpi_abort()
            END IF
            ! ==================================================================
            
            
            ! Simulation Algorithm Parameters ==================================
            IF(ALL(model_eqns /= (/1,2,3,4/))) THEN
                PRINT '(A)', 'Unsupported value of model_eqns. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF( model_eqns == 2 .AND. bubbles .AND. bubble_model == 1  ) THEN
                PRINT '(A)', 'The 5-equation bubbly flow model requires bubble_model = 2 (Keller--Miksis)'
                CALL s_mpi_abort()
            ELSEIF( bubbles .AND. bubble_model == 3 .AND. (polytropic .NEQV. .TRUE.)  ) THEN
                PRINT '(A)', 'RP bubbles require polytropic compression'
                CALL s_mpi_abort()
            ELSEIF( cyl_coord .AND. bubbles ) THEN
                PRINT '(A)', 'Bubble models untested in cylindrical coordinates'
                CALL s_mpi_abort()
            ELSEIF( model_eqns == 3 .AND. bubbles ) THEN
                PRINT '(A)', 'Bubble models untested with 6-equation model'
                CALL s_mpi_abort()
            ELSEIF( model_eqns == 1 .AND. bubbles ) THEN
                PRINT '(A)', 'Bubble models untested with pi-gamma model'
                CALL s_mpi_abort()
            ELSEIF( model_eqns == 4 .AND. num_fluids /= 1 ) THEN
                PRINT '(A)', 'The 4-equation model implementation is not a multi-component and requires num_fluids = 1'
                CALL s_mpi_abort()
            ELSEIF( bubbles .AND. weno_vars /= 2 ) THEN
                PRINT '(A)', 'Bubble modeling requires weno_vars = 2'
                CALL s_mpi_abort()
            ELSEIF( bubbles .AND. riemann_solver /= 2 ) THEN
                PRINT '(A)', 'Bubble modeling requires riemann_solver = 2'
                CALL s_mpi_abort()
            ELSEIF( bubbles .AND. commute_err ) THEN
                PRINT '(A)', 'Bubble modeling is not compatible with commute_err'
                CALL s_mpi_abort()
            ELSEIF( bubbles .AND. split_err ) THEN
                PRINT '(A)', 'Bubble modeling is not compatible with split_err'
                CALL s_mpi_abort()
            ELSEIF( bubbles .AND. regularization ) THEN
                PRINT '(A)', 'Bubble modeling is not compatible with regularization'
                CALL s_mpi_abort()
            ELSEIF( bubbles .AND. (adv_alphan .NEQV. .TRUE.) ) THEN
                PRINT '(A)', 'Bubble modeling requires adv_alphan'
                CALL s_mpi_abort()
            ELSEIF( (bubbles  .NEQV. .TRUE.) .AND. polydisperse ) THEN
                PRINT '(A)', 'Polydisperse bubble modeling requires the bubble switch to be activated'
                CALL s_mpi_abort()
            ELSEIF( polydisperse .and. (poly_sigma == dflt_real) ) THEN
                PRINT '(A)', 'Polydisperse bubble modeling requires poly_sigma > 0'
                CALL s_mpi_abort()
            ELSEIF( qbmm.AND. (bubbles .NEQV. .TRUE.) ) THEN
                PRINT '(A)', 'QBMM requires bubbles'
                CALL s_mpi_abort()
            ELSEIF( qbmm .AND. (nnode .NE. 4) ) THEN
                PRINT '(A)', 'nnode not supported'
                CALL s_mpi_abort()
            ELSEIF(model_eqns == 3 .AND. riemann_solver /= 2) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'model_eqns (6-eq) and riemann_solver (please use riemann_solver = 2). '  // &
                             'Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(model_eqns == 3 .AND. (alt_soundspeed .OR. regularization)) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'model_eqns (6-eq) and alt_soundspeed or regularization. '  // &
                             'Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(model_eqns == 3 .AND. (tvd_riemann_flux .OR. tvd_wave_speeds)) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'model_eqns (6-eq) and tvd_riemann_flux or tvd_wave_speeds. '  // &
                             'Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(model_eqns == 3 .AND. avg_state == 1) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'model_eqns (6-eq) and Roe average (please use avg_state = 2). '  // &
                             'Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(model_eqns == 3 .AND. wave_speeds == 2) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'model_eqns (6-eq) and wave_speeds (please use wave_speeds = 1). '  // &
                             'Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(model_eqns == 3 .AND. (cyl_coord .AND. p /= 0)) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'model_eqns (6-eq) and cylindrical coordinates. '  // &
                             'Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(model_eqns == 3 .AND. (commute_err .OR. split_err)) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'model_eqns (6-eq) and commute_err or split_err. '  // &
                             'Exiting ...'
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
            ELSEIF(model_eqns == 1 .AND. mpp_lim) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'model_eqns and mpp_lim. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(num_fluids == 1 .AND. mpp_lim) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'num_fluids and mpp_lim. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(time_stepper < 1 .OR. time_stepper > 5) THEN
                IF (time_stepper /= 23 ) THEN
                    PRINT '(A)', 'Unsupported value of time_stepper. Exiting ...'
                    CALL s_mpi_abort()
                END IF
            ELSEIF(t_tol == dflt_real .AND. time_stepper == 23 ) THEN
                PRINT '(A)', 'Adaptive timestepping requires a tolerance t_tol'
                CALL s_mpi_abort()
            ELSEIF(ALL(weno_vars /= (/1,2/))) THEN
                PRINT '(A)', 'Unsupported value of weno_vars. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(ALL(weno_order /= (/1,3,5/))) THEN
                PRINT '(A)', 'Unsupported value of weno_order. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(m+1 < num_stcls_min*weno_order) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'm and weno_order. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(n+1 < MIN(1,n)*num_stcls_min*weno_order) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'n and weno_order. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(p+1 < MIN(1,p)*num_stcls_min*weno_order) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'p and weno_order. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(weno_eps <= 0d0 .OR. weno_eps > 1d-6) THEN
                PRINT '(A)', 'Unsupported value of weno_eps. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(weno_order == 1 .AND. char_decomp) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'weno_order and char_decomp. '          // &
                             'Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(weno_order == 1 .AND. mapped_weno) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'weno_order and mapped_weno. '          // &
                             'Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(weno_order /= 5 .AND. mp_weno) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'weno_order and mp_weno. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF( (model_eqns == 1 .AND. weno_avg) .OR. (model_eqns == 4 .AND. weno_avg) ) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'model_eqns and weno_avg. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(riemann_solver < 1 .OR. riemann_solver > 3) THEN
                PRINT '(A)', 'Unsupported value of riemann_solver. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(ALL(wave_speeds /= (/dflt_int,1,2/))) THEN
                PRINT '(A)', 'Unsupported value of wave_speeds. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF( (riemann_solver /= 3 .AND. wave_speeds == dflt_int) &
                                           .OR.                         &
                    (riemann_solver == 3 .AND. wave_speeds /= dflt_int) ) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'riemann_solver and wave_speeds. '      // &
                             'Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(ALL(avg_state /= (/dflt_int,1,2/))) THEN
                PRINT '(A)', 'Unsupported value of avg_state. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(char_decomp .AND. avg_state == dflt_int) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'char_decomp and avg_state. '           // &
                             'Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(riemann_solver /= 3 .AND. avg_state == dflt_int) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'riemann_solver and avg_state. '        // &
                             'Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(weno_order == 1 .AND. commute_err) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'weno_order and commute_err. '          // &
                             'Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(n == 0 .AND. split_err) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'n and split_err. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(weno_order == 1 .AND. split_err) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'weno_order and split_err. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(bc_x%beg < -12 .OR. bc_x%beg > -1) THEN
                PRINT '(A)', 'Unsupported value of bc_x%beg. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(bc_x%end < -12 .OR. bc_x%end > -1) THEN
                PRINT '(A)', 'Unsupported value of bc_x%end. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF( (bc_x%beg == -1 .AND. bc_x%end /= -1) &
                                    .OR.                  &
                    (bc_x%end == -1 .AND. bc_x%beg /= -1) ) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'bc_x%beg and bc_x%end. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF (        bc_y%beg /= dflt_int        &
                                    .AND.               &
                   ( ((cyl_coord .NEQV. .TRUE.) .AND. (bc_y%beg < -12 .OR. bc_y%beg > -1)) & 
                                            .OR.        &
                     (cyl_coord .AND. p == 0 .AND. bc_y%beg /= -2) &
                                            .OR.        &
                     (cyl_coord .AND. p > 0 .AND. bc_y%beg /= -13) )) THEN
                PRINT '(A)', 'Unsupported value of bc_y%beg. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(         bc_y%end /= dflt_int        &
                                    .AND.               &
                    (bc_y%end < -12 .OR. bc_y%end > -1) ) THEN
                PRINT '(A)', 'Unsupported value of bc_y%end. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF( (n == 0 .AND. bc_y%beg /= dflt_int) &
                                    .OR.                &
                    (n  > 0 .AND. bc_y%beg == dflt_int) ) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'n and bc_y%beg. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF( (n == 0 .AND. bc_y%end /= dflt_int) &
                                    .OR.                &
                    (n >  0 .AND. bc_y%end == dflt_int) ) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'n and bc_y%end. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF( (bc_y%beg == -1 .AND. bc_y%end /= -1) &
                                    .OR.                  &
                    (bc_y%end == -1 .AND. bc_y%beg /= -1) ) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'bc_y%beg and bc_y%end. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(         bc_z%beg /= dflt_int        &
                                    .AND.               &
                    (bc_z%beg < -12 .OR. bc_z%beg > -1) ) THEN
                PRINT '(A)', 'Unsupported value of bc_z%beg. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(         bc_z%end /= dflt_int        &
                                    .AND.               &
                    (bc_z%end < -12 .OR. bc_z%end > -1) ) THEN
                PRINT '(A)', 'Unsupported value of bc_z%end. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF( (p == 0 .AND. bc_z%beg /= dflt_int) &
                                    .OR.                &
                    (p >  0 .AND. bc_z%beg == dflt_int) ) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'p and bc_z%beg. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF( (p == 0 .AND. bc_z%end /= dflt_int) &
                                    .OR.                &
                    (p >  0 .AND. bc_z%end == dflt_int) ) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'p and bc_z%end. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF( (bc_z%beg == -1 .AND. bc_z%end /= -1) &
                                    .OR.                  &
                    (bc_z%end == -1 .AND. bc_z%beg /= -1) ) THEN
                PRINT '(A)', 'Unsupported combination of values of ' // &
                             'bc_z%beg and bc_z%end. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF( (ANY(threshold_mf /= dflt_real)) &
                    .AND.           &
                ( ALL(cb_wrt .NEQV. .TRUE.))    ) THEN
                PRINT '(A)', 'Unsupported combination of cb_wrt ' // &
                        'and threshold_mf. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF( (ANY(moment_order /= dflt_int)) &
                    .AND.           &
                ( ALL(com_wrt .NEQV. .TRUE.))   ) THEN
                PRINT '(A)', 'Unsupported combination of com_wrt ' // &
                        'and moment_order. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF( ANY(cb_wrt) .AND. (ALL(threshold_mf == dflt_real))) THEN
                PRINT '(A)', 'Unsupported combination of cb_wrt ' // &
                        'and threshold_mf. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(model_eqns == 1 .AND. alt_soundspeed) THEN
                PRINT '(A)', 'Unsupported combination of model_eqns ' // &
                        'and alt_soundspeed. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(model_eqns == 4 .AND. alt_soundspeed) THEN
                PRINT '(A)', 'Unsupported combination of model_eqns ' // &
                        'and alt_soundspeed. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(  (num_fluids /= 2 .AND. num_fluids /= 3) .AND. alt_soundspeed) THEN
                PRINT '(A)', 'Unsupported combination of num_fluids ' // &
                        'and alt_soundspeed. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(riemann_solver /= 2 .AND. alt_soundspeed) THEN
                PRINT '(A)', 'Unsupported combination of riemann_solver ' // &
                        'and alt_soundspeed. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(model_eqns == 1 .AND. regularization) THEN
                PRINT '(A)', 'Unsupported combination of model_eqns ' // &
                        'and regularization. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(model_eqns == 4 .AND. regularization) THEN
                PRINT '(A)', 'Unsupported combination of model_eqns ' // &
                        'and regularization. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(num_fluids /= 2 .AND. regularization) THEN
                PRINT '(A)', 'Unsupported combination of num_fluids ' // &
                        'and regularization. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF(riemann_solver /= 2 .AND. regularization) THEN
                PRINT '(A)', 'Unsupported combination of riemann_solver ' // &
                        'and regularization. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF (regularization .AND. reg_eps == dflt_real) THEN
                PRINT '(A)', 'Unsupported combination of regularization ' // &
                        'and value of reg_eps. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF (regularization .AND. n == 0) THEN
                PRINT '(A)', 'Unsupported combination of regularization ' // &
                        'and value of n. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF (tvd_rhs_flux .AND. tvd_riemann_flux) THEN
                PRINT '(A)', 'Unsupported combination of tvd_rhs_flux ' // &
                        'and tvd_riemann_flux. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF (tvd_wave_speeds .AND. (tvd_riemann_flux .NEQV. .TRUE.)) THEN
                PRINT '(A)', 'Unsupported combination of tvd_riemann_flux ' // &
                        'and tvd_wave_speeds. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF ( (tvd_rhs_flux .OR. tvd_riemann_flux .OR. tvd_wave_speeds) &
                                            .AND.                           &
                                (flux_lim == dflt_int) ) THEN
                PRINT '(A)', 'Unsupported combination of tvd_rhs_flux/' // &
                        'tvd_riemann_flux/tvd_wave_speeds and flux_lim. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF ( (flux_lim /= dflt_int) .AND. (tvd_rhs_flux .NEQV. .TRUE.) .AND. &
                    (tvd_riemann_flux .NEQV. .TRUE.) .AND. (tvd_wave_speeds .NEQV. .TRUE.)) THEN
                PRINT '(A)', 'Unsupported combination of tvd_rhs_flux/' // &
                        'tvd_riemann_flux/tvd_wave_speeds and flux_lim. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF (tvd_riemann_flux .AND. (riemann_solver == 3)) THEN
                PRINT '(A)', 'Unsupported combination of tvd_riemann_flux ' // &
                        'and choice of riemann_solver. Exiting ...'
            ELSEIF (ALL(flux_lim /= (/dflt_int,1,2,3,4,5,6,7/))) THEN
                PRINT '(A)', 'Unsupported value of flux_lim. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF ( (We_riemann_flux .NEQV. .TRUE.) .AND. We_wave_speeds) THEN
                PRINT '(A)', 'Unsupported combination of We_riemann_flux and ' // &
                        'We_wave_speeds. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF (We_rhs_flux .AND. We_riemann_flux) THEN
                PRINT '(A)', 'Unsupported combination of We_rhs_flux ' // &
                        'and We_riemann_flux. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF (We_src .AND. We_riemann_flux) THEN
                PRINT '(A)', 'Unsupported combination of We_src ' // &
                        'and We_riemann_flux. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF (We_rhs_flux .AND. We_src) THEN
                PRINT '(A)', 'Unsupported combination of We_rhs_flux ' // &
                        'and We_src. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF ( (We_riemann_flux .OR. We_rhs_flux .OR. We_src) &
                        .AND. (lsq_deriv .NEQV. .TRUE.)) THEN
                PRINT '(A)', 'Unsupported combination of We_riemann_flux, ' // &
                        'We_rhs_flux, We_src, and lsq_deriv. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF (lsq_deriv .AND. n == 0) THEN
                PRINT '(A)', 'Unsupported combination of lsq_deriv ' // &
                        'and value of n. Exiting ...'
                CALL s_mpi_abort()
            END IF
            ! END: Simulation Algorithm Parameters =============================
            

            ! Finite Difference Parameters =====================================
            IF(                     fd_order /= dflt_int               &
                                             .AND.                        &
                  fd_order /= 1 .AND. fd_order /= 2 .AND. fd_order /= 4   ) THEN
                PRINT '(A)', 'Unsupported choice for the value of '         // &
                             'fd_order. Exiting ...'
                CALL s_mpi_abort()
            ELSEIF (probe_wrt .AND. fd_order == dflt_int) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for probe_wrt, and fd_order. '         // &
                             'Exiting ...'
                CALL s_mpi_abort()
            ELSEIF (integral_wrt .AND. (bubbles .NEQV. .TRUE.)) THEN
                PRINT '(A)', 'Unsupported choice of the combination of '    // &
                             'values for integral_wrt, and bubbles. '         // &
                             'Exiting ...'
                CALL s_mpi_abort()
            END IF
            ! END: Finite Difference Parameters ================================
            


            ! Fluids Physical Parameters =======================================
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
                
                DO j = 1,2
                    
                    IF( fluid_pp(i)%Re(j) /= dflt_real &
                                     .AND.             &
                        fluid_pp(i)%Re(j) <=    0d0    ) THEN
                        PRINT '(A,I0,A,I0,A)', 'Unsupported value of '  // &
                                               'fluid_pp(',i,')%'       // &
                                               'Re(',j,'). Exiting ...'
                        CALL s_mpi_abort()
                    END IF
                    
                    IF(         model_eqns == 1        &
                                     .AND.             &
                        fluid_pp(i)%Re(j) /= dflt_real ) THEN
                        PRINT '(A,I0,A,I0,A)', 'Unsupported combination ' // &
                                               'of values of model_eqns ' // &
                                               'and fluid_pp(',i,')%'     // &
                                               'Re(',j,'). Exiting ...'
                        CALL s_mpi_abort()
                    END IF
                    
                    IF(         i > num_fluids         &
                                     .AND.             &
                        fluid_pp(i)%Re(j) /= dflt_real ) THEN
                        PRINT '(A,I0,A,I0,A)', 'Unsupported combination ' // &
                                               'of values of num_fluids ' // &
                                               'and fluid_pp(',i,')%'     // &
                                               'Re(',j,'). Exiting ...'
                        CALL s_mpi_abort()
                    END IF
                    
                    IF(         weno_order == 1        &
                                     .AND.             &
                           (weno_avg .NEQV. .TRUE.)    &
                                     .AND.             &
                        fluid_pp(i)%Re(j) /= dflt_real ) THEN
                        PRINT '(A,I0,A,I0,A)', 'Unsupported combination '  // &
                                               'of values of weno_order, ' // &
                                               'weno_avg and fluid_pp'     // &
                                               '(',i,')%Re(',j,'). '       // &
                                               'Exiting ...'
                        CALL s_mpi_abort()
                    END IF
                    
                END DO
                
                DO j = 1, num_fluids_max
                    
                    IF((            i == j             &
                                     .OR.              &
                        fluid_pp(i)%We(j) <=    0d0    &
                                     .OR.              &
                        fluid_pp(j)%We(i) == dflt_real)&
                                     .AND.             &
                        fluid_pp(i)%We(j) /= dflt_real ) THEN
                        PRINT '(A,I0,A,I0,A)', 'Unsupported value of '  // &
                                               'fluid_pp(',i,')%'       // &
                                               'We(',j,'). Exiting ...'
                        CALL s_mpi_abort()
                    END IF
                    
                    IF(         model_eqns == 1        &
                                     .AND.             &
                        fluid_pp(i)%We(j) /= dflt_real ) THEN
                        PRINT '(A,I0,A,I0,A)', 'Unsupported combination ' // &
                                               'of values of model_eqns ' // &
                                               'and fluid_pp(',i,')%'     // &
                                               'We(',j,'). Exiting ...'
                        CALL s_mpi_abort()
                    END IF
                    
                    IF(      MAX(i,j) > num_fluids     &
                                     .AND.             &
                        fluid_pp(i)%We(j) /= dflt_real ) THEN
                        PRINT '(A,I0,A,I0,A)', 'Unsupported combination ' // &
                                               'of values of num_fluids ' // &
                                               'and fluid_pp(',i,')%'     // &
                                               'We(',j,'). Exiting ...'
                        CALL s_mpi_abort()
                    END IF
                    
                    IF(   ANY(num_fluids == (/i,j/))   &
                                     .AND.             &
                          (adv_alphan .NEQV. .TRUE.)  &
                                     .AND.             &
                        fluid_pp(i)%We(j) /= dflt_real ) THEN
                        PRINT '(A,I0,A,I0,A)', 'Unsupported combination ' // &
                                               'of values of num_fluids ' // &
                                               'adv_alphan and fluid_pp'  // &
                                               '(',i,')%We(',j,'). '      // &
                                               'Exiting ...'
                        CALL s_mpi_abort()
                    END IF
                    
                    IF(           char_decomp          &
                                     .AND.             &
                        fluid_pp(i)%We(j) /= dflt_real ) THEN
                        PRINT '(A,I0,A,I0,A)', 'Unsupported combination '  // &
                                               'of values of char_decomp ' // &
                                               'and fluid_pp(',i,')%'      // &
                                               'We(',j,'). Exiting ...'
                        CALL s_mpi_abort()
                    END IF
                    
                    IF(         weno_order == 1        &
                                     .AND.             &
                           (weno_avg .NEQV. .TRUE.)    &
                                     .AND.             &
                        fluid_pp(i)%We(j) /= dflt_real ) THEN
                        PRINT '(A,I0,A,I0,A)', 'Unsupported combination '  // &
                                               'of values of weno_order, ' // &
                                               'weno_avg and fluid_pp'     // &
                                               '(',i,')%We(',j,'). '       // &
                                               'Exiting ...'
                        CALL s_mpi_abort()
                    END IF
                    
                    IF((We_riemann_flux .NEQV. .TRUE.) &
                                     .AND.             &
                        (We_rhs_flux .NEQV. .TRUE.)    &
                                     .AND.             &
                           (We_src .NEQV. .TRUE.)      &
                                     .AND.             &
                        fluid_pp(i)%We(j) /= dflt_real ) THEN
                        PRINT '(A,I0,A,I0,A)', 'Unsupported combination '  // &
                                               'of values of We_rhs_flux, '// &
                                               'We_riemann_flux, We_src '  // &
                                               'and fluid_pp'     // &
                                               '(',i,')%We(',j,'). '       // &
                                               'Exiting ...'
                        CALL s_mpi_abort()
                    END IF
                    
                END DO
                
            END DO
            ! END: Fluids Physical Parameters ==================================
            
            
        END SUBROUTINE s_check_input_file ! ------------------------------------
        
        
        
        
        !> The primary purpose of this procedure is to read in the
        !!              initial condition and grid data files. The cell-average
        !!              conservative variables constitute the former, while the
        !!              cell-boundary locations in x-, y- and z-directions make
        !!              up the latter. This procedure also calculates the cell-
        !!              width distributions from the cell-boundary locations.
        !! @param q_cons_vf Cell-averaged conservative variables
        SUBROUTINE s_read_serial_data_files(q_cons_vf) ! ------------------------------

            
            
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(INOUT) :: q_cons_vf
            

            CHARACTER(LEN = path_len + 2*name_len) :: t_step_dir !<
            !! Relative path to the starting time-step directory
            
            CHARACTER(LEN = path_len + 3*name_len) :: file_path !<
            !! Relative path to the grid and conservative variables data files


            LOGICAL :: file_exist !< 
            ! Logical used to check the existence of the data files
            

            INTEGER :: i !< Generic loop iterator
            
            ! Confirming that the directory from which the initial condition and
            ! the grid data files are to be read in exists and exiting otherwise
            WRITE(t_step_dir, '(A,I0,A,I0)') &
                    TRIM(case_dir) // '/p_all/p', proc_rank, '/', t_step_start
            
            file_path = TRIM(t_step_dir) // '/.'
            CALL my_inquire(file_path,file_exist)

            IF(file_exist .NEQV. .TRUE.) THEN
                PRINT '(A)', TRIM(file_path) // ' is missing. Exiting ...'
                CALL s_mpi_abort()
            END IF
            
            
            ! Cell-boundary Locations in x-direction ===========================
            file_path = TRIM(t_step_dir) // '/x_cb.dat'
            
            INQUIRE(FILE = TRIM(file_path), EXIST = file_exist)
            
            IF(file_exist) THEN
                OPEN(2, FILE   = TRIM(file_path), &
                        FORM   = 'unformatted'  , &
                        ACTION = 'read'         , &
                        STATUS = 'old'            )
                READ(2) x_cb(-1:m); CLOSE(2)
            ELSE
                PRINT '(A)', TRIM(file_path) // ' is missing. Exiting ...'
                CALL s_mpi_abort()
            END IF
            
            dx(0:m) = x_cb(0:m) - x_cb(-1:m-1)
            x_cc(0:m) = x_cb(-1:m-1) + dx(0:m)/2d0
            ! ==================================================================
            
            
            ! Cell-boundary Locations in y-direction ===========================
            IF(n > 0) THEN
                
                file_path = TRIM(t_step_dir) // '/y_cb.dat'
                
                INQUIRE(FILE = TRIM(file_path), EXIST = file_exist)
                
                IF(file_exist) THEN
                    OPEN(2, FILE   = TRIM(file_path), &
                            FORM   = 'unformatted'  , &
                            ACTION = 'read'         , &
                            STATUS = 'old'            )
                    READ(2) y_cb(-1:n); CLOSE(2)
                ELSE
                    PRINT '(A)', TRIM(file_path) // ' is missing. Exiting ...'
                    CALL s_mpi_abort()
                END IF
                
                dy(0:n) = y_cb(0:n) - y_cb(-1:n-1)
                y_cc(0:n) = y_cb(-1:n-1) + dy(0:n)/2d0
                
            END IF
            ! ==================================================================
            
            
            ! Cell-boundary Locations in z-direction ===========================
            IF(p > 0) THEN
                
                file_path = TRIM(t_step_dir) // '/z_cb.dat'
                
                INQUIRE(FILE = TRIM(file_path), EXIST = file_exist)
                
                IF(file_exist) THEN
                    OPEN(2, FILE   = TRIM(file_path), &
                            FORM   = 'unformatted'  , &
                            ACTION = 'read'         , &
                            STATUS = 'old'            )
                    READ(2) z_cb(-1:p); CLOSE(2)
                ELSE
                    PRINT '(A)', TRIM(file_path) // ' is missing. Exiting ...'
                    CALL s_mpi_abort()
                END IF
                
                dz(0:p) = z_cb(0:p) - z_cb(-1:p-1)
                z_cc(0:p) = z_cb(-1:p-1) + dz(0:p)/2d0
                
            END IF
            ! ==================================================================
            
            
            ! Cell-average Conservative Variables ==============================
            IF ((bubbles .NEQV. .True.) .AND. (hypoelasticity .NEQV. .TRUE.)) THEN
                DO i = 1, adv_idx%end
                    WRITE(file_path, '(A,I0,A)') &
                           TRIM(t_step_dir) // '/q_cons_vf', i, '.dat'
                    INQUIRE(FILE = TRIM(file_path), EXIST = file_exist)
                    IF(file_exist) THEN
                        OPEN(2, FILE   = TRIM(file_path), &
                            FORM   = 'unformatted'  , &
                            ACTION = 'read'         , &
                            STATUS = 'old'            )
                        READ(2) q_cons_vf(i)%sf(0:m,0:n,0:p); CLOSE(2)
                    ELSE
                        PRINT '(A)', TRIM(file_path) // ' is missing. Exiting ...'
                        CALL s_mpi_abort()
                    END IF
                END DO
            ELSE 
                !make sure to read bubble variables
                DO i = 1, sys_size
                    WRITE(file_path, '(A,I0,A)') &
                           TRIM(t_step_dir) // '/q_cons_vf', i, '.dat'
                    INQUIRE(FILE = TRIM(file_path), EXIST = file_exist)
                    IF(file_exist) THEN
                        OPEN(2, FILE   = TRIM(file_path), &
                            FORM   = 'unformatted'  , &
                            ACTION = 'read'         , &
                            STATUS = 'old'            )
                        READ(2) q_cons_vf(i)%sf(0:m,0:n,0:p); CLOSE(2)
                    ELSE
                        PRINT '(A)', TRIM(file_path) // ' is missing. Exiting ...'
                        CALL s_mpi_abort()
                    END IF
                END DO
            END IF
            ! ==================================================================
            
        END SUBROUTINE s_read_serial_data_files ! -------------------------------------
        
        
        
        
        !! @param q_cons_vf Conservative variables
        SUBROUTINE s_read_parallel_data_files(q_cons_vf) ! ---------------------------

            TYPE(scalar_field), &
            DIMENSION(sys_size), &
            INTENT(INOUT) :: q_cons_vf

            REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:) :: x_cb_glb, y_cb_glb, z_cb_glb

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

            ALLOCATE(x_cb_glb(-1:m_glb))
            ALLOCATE(y_cb_glb(-1:n_glb))
            ALLOCATE(z_cb_glb(-1:p_glb))

            ! Read in cell boundary locations in x-direction
            file_loc = TRIM(case_dir) // '/restart_data' // TRIM(mpiiofs) // 'x_cb.dat'
            INQUIRE(FILE = TRIM(file_loc),EXIST = file_exist)

            IF (file_exist) THEN
                data_size = m_glb+2
                CALL MPI_FILE_OPEN(MPI_COMM_WORLD,file_loc,MPI_MODE_RDONLY,mpi_info_int,ifile,ierr)
                CALL MPI_FILE_READ(ifile,x_cb_glb,data_size,MPI_DOUBLE_PRECISION,status,ierr)
                CALL MPI_FILE_CLOSE(ifile,ierr)
            ELSE
                PRINT '(A)', 'File ', TRIM(file_loc), ' is missing. Exiting...'
                CALL s_mpi_abort()
            END IF

            ! Assigning local cell boundary locations
            x_cb(-1:m) = x_cb_glb((start_idx(1)-1):(start_idx(1)+m))
            ! Computing the cell width distribution
            dx(0:m) = x_cb(0:m) - x_cb(-1:m-1)
            ! Computing the cell center locations
            x_cc(0:m) = x_cb(-1:m-1) + dx(0:m)/2d0

            IF (n > 0) THEN
                ! Read in cell boundary locations in y-direction
                file_loc = TRIM(case_dir) // '/restart_data' // TRIM(mpiiofs) // 'y_cb.dat'
                INQUIRE(FILE = TRIM(file_loc),EXIST = file_exist)
        
                IF (file_exist) THEN
                    data_size = n_glb+2
                    CALL MPI_FILE_OPEN(MPI_COMM_WORLD,file_loc,MPI_MODE_RDONLY,mpi_info_int,ifile,ierr)
                    CALL MPI_FILE_READ(ifile,y_cb_glb,data_size,MPI_DOUBLE_PRECISION,status,ierr)
                    CALL MPI_FILE_CLOSE(ifile,ierr)
                ELSE
                    PRINT '(A)', 'File ', TRIM(file_loc), ' is missing. Exiting...'
                    CALL s_mpi_abort()
                END IF
        
                ! Assigning local cell boundary locations
                y_cb(-1:n) = y_cb_glb((start_idx(2)-1):(start_idx(2)+n))
                ! Computing the cell width distribution
                dy(0:n) = y_cb(0:n) - y_cb(-1:n-1)
                ! Computing the cell center locations
                y_cc(0:n) = y_cb(-1:n-1) + dy(0:n)/2d0

                IF (p > 0) THEN
                    ! Read in cell boundary locations in z-direction
                    file_loc = TRIM(case_dir) // '/restart_data' // TRIM(mpiiofs) // 'z_cb.dat'
                    INQUIRE(FILE = TRIM(file_loc),EXIST = file_exist)
            
                    IF (file_exist) THEN
                        data_size = p_glb+2
                        CALL MPI_FILE_OPEN(MPI_COMM_WORLD,file_loc,MPI_MODE_RDONLY,mpi_info_int,ifile,ierr)
                        CALL MPI_FILE_READ(ifile,z_cb_glb,data_size,MPI_DOUBLE_PRECISION,status,ierr)
                        CALL MPI_FILE_CLOSE(ifile,ierr)
                    ELSE
                        PRINT '(A)', 'File ', TRIM(file_loc), ' is missing. Exiting...'
                        CALL s_mpi_abort()
                    END IF
            
                    ! Assigning local cell boundary locations
                    z_cb(-1:p) = z_cb_glb((start_idx(3)-1):(start_idx(3)+p))
                    ! Computing the cell width distribution
                    dz(0:p) = z_cb(0:p) - z_cb(-1:p-1)
                    ! Computing the cell center locations
                    z_cc(0:p) = z_cb(-1:p-1) + dz(0:p)/2d0

                END IF
            END IF

            ! Open the file to read conservative variables
            WRITE(file_loc, '(I0,A)') t_step_start, '.dat'
            file_loc = TRIM(case_dir) // '/restart_data' // TRIM(mpiiofs) // TRIM(file_loc)
            INQUIRE(FILE = TRIM(file_loc),EXIST = file_exist)

            IF (file_exist) THEN
                CALL MPI_FILE_OPEN(MPI_COMM_WORLD,file_loc,MPI_MODE_RDONLY,mpi_info_int,ifile,ierr)

                ! Initialize MPI data I/O
                CALL s_initialize_mpi_data(q_cons_vf)

                ! Size of local arrays
                data_size = (m+1)*(n+1)*(p+1)

                ! Resize some integers so MPI can read even the biggest file
                m_MOK     = INT(m_glb+1,     MPI_OFFSET_KIND)
                n_MOK     = INT(n_glb+1,     MPI_OFFSET_KIND)
                p_MOK     = INT(p_glb+1,     MPI_OFFSET_KIND)
                WP_MOK    = INT(8d0,      MPI_OFFSET_KIND)
                MOK       = INT(1d0,      MPI_OFFSET_KIND)
                str_MOK   = INT(name_len, MPI_OFFSET_KIND)
                NVARS_MOK = INT(sys_size, MPI_OFFSET_KIND)

                ! Read the data for each variable
                IF ((bubbles .EQV. .TRUE.) .OR. (hypoelasticity .EQV. .TRUE.)) THEN
                    DO i = 1, sys_size!adv_idx%end
                        var_MOK = INT(i, MPI_OFFSET_KIND)

                        ! Initial displacement to skip at beginning of file
                        disp = m_MOK*MAX(MOK,n_MOK)*MAX(MOK,p_MOK)*WP_MOK*(var_MOK-1)

                        CALL MPI_FILE_SET_VIEW(ifile,disp,MPI_DOUBLE_PRECISION,MPI_IO_DATA%view(i), &
                                    'native',mpi_info_int,ierr)
                        CALL MPI_FILE_READ(ifile,MPI_IO_DATA%var(i)%sf,data_size, &
                                    MPI_DOUBLE_PRECISION,status,ierr)
                    END DO
                ELSE
                    DO i = 1, adv_idx%end
                        var_MOK = INT(i, MPI_OFFSET_KIND)

                        ! Initial displacement to skip at beginning of file
                        disp = m_MOK*MAX(MOK,n_MOK)*MAX(MOK,p_MOK)*WP_MOK*(var_MOK-1)

                        CALL MPI_FILE_SET_VIEW(ifile,disp,MPI_DOUBLE_PRECISION,MPI_IO_DATA%view(i), &
                                    'native',mpi_info_int,ierr)
                        CALL MPI_FILE_READ(ifile,MPI_IO_DATA%var(i)%sf,data_size, &
                                    MPI_DOUBLE_PRECISION,status,ierr)
                    END DO
                END IF

                CALL s_mpi_barrier()

                CALL MPI_FILE_CLOSE(ifile,ierr)
            ELSE
                PRINT '(A)', 'File ', TRIM(file_loc), ' is missing. Exiting...'
                CALL s_mpi_abort()
            END IF

            DEALLOCATE(x_cb_glb, y_cb_glb, z_cb_glb)
        END SUBROUTINE s_read_parallel_data_files ! -------------------------------



        !> The purpose of this subroutine is to populate the buffers
        !!          of the grid variables, which are constituted of the cell-
        !!          boundary locations and cell-width distributions, based on
        !!          the boundary conditions.
        SUBROUTINE s_populate_grid_variables_buffers() ! -----------------------

            
            

            INTEGER :: i !< Generic loop iterator
            
            
            ! Population of Buffers in x-direction =============================
            
            ! Populating cell-width distribution buffer, at the beginning of the
            ! coordinate direction, based on the selected boundary condition. In
            ! order, these are the ghost-cell extrapolation, symmetry, periodic,
            ! and processor boundary conditions.
            IF(bc_x%beg <= -3) THEN
                DO i = 1, buff_size
                    dx(-i) = dx(0)
                END DO
            ELSEIF(bc_x%beg == -2) THEN
                DO i = 1, buff_size
                    dx(-i) = dx(i-1)
                END DO
            ELSEIF(bc_x%beg == -1) THEN
                DO i = 1, buff_size
                    dx(-i) = dx(m-(i-1))
                END DO
            ELSE
                CALL s_mpi_sendrecv_grid_variables_buffers(1,-1)
            END IF
            
            ! Computing the cell-boundary locations buffer, at the beginning of
            ! the coordinate direction, from the cell-width distribution buffer
            DO i = 1, buff_size
                x_cb(-1-i) = x_cb(-i) - dx(-i)
            END DO
            ! Computing the cell-center locations buffer, at the beginning of
            ! the coordinate direction, from the cell-width distribution buffer
            DO i = 1, buff_size
                x_cc(-i) = x_cc(1-i) - (dx(1-i) + dx(-i))/2d0
            END DO
            
            ! Populating the cell-width distribution buffer, at the end of the
            ! coordinate direction, based on desired boundary condition. These
            ! include, in order, ghost-cell extrapolation, symmetry, periodic,
            ! and processor boundary conditions.
            IF(bc_x%end <= -3) THEN
                DO i = 1, buff_size
                    dx(m+i) = dx(m)
                END DO
            ELSEIF(bc_x%end == -2) THEN
                DO i = 1, buff_size
                    dx(m+i) = dx(m-(i-1))
                END DO
            ELSEIF(bc_x%end == -1) THEN
                DO i = 1, buff_size
                    dx(m+i) = dx(i-1)
                END DO
            ELSE
                CALL s_mpi_sendrecv_grid_variables_buffers(1,1)
            END IF
            
            ! Populating the cell-boundary locations buffer, at the end of the
            ! coordinate direction, from buffer of the cell-width distribution
            DO i = 1, buff_size
                x_cb(m+i) = x_cb(m+(i-1)) + dx(m+i)
            END DO
            ! Populating the cell-center locations buffer, at the end of the
            ! coordinate direction, from buffer of the cell-width distribution
            DO i = 1, buff_size
                x_cc(m+i) = x_cc(m+(i-1)) + (dx(m+(i-1)) + dx(m+i))/2d0
            END DO
            
            ! END: Population of Buffers in x-direction ========================
            
            
            ! Population of Buffers in y-direction =============================
            
            ! Populating cell-width distribution buffer, at the beginning of the
            ! coordinate direction, based on the selected boundary condition. In
            ! order, these are the ghost-cell extrapolation, symmetry, periodic,
            ! and processor boundary conditions.
            IF(n == 0) THEN
                RETURN
            ELSEIF(bc_y%beg <= -3 .AND. bc_y%beg /= -13) THEN
                DO i = 1, buff_size
                    dy(-i) = dy(0)
                END DO
            ELSEIF(bc_y%beg == -2 .OR. bc_y%beg == -13) THEN
                DO i = 1, buff_size
                    dy(-i) = dy(i-1)
                END DO
            ELSEIF(bc_y%beg == -1) THEN
                DO i = 1, buff_size
                    dy(-i) = dy(n-(i-1))
                END DO
            ELSE
                CALL s_mpi_sendrecv_grid_variables_buffers(2,-1)
            END IF
            
            ! Computing the cell-boundary locations buffer, at the beginning of
            ! the coordinate direction, from the cell-width distribution buffer
            DO i = 1, buff_size
                y_cb(-1-i) = y_cb(-i) - dy(-i)
            END DO
            ! Computing the cell-center locations buffer, at the beginning of
            ! the coordinate direction, from the cell-width distribution buffer
            DO i = 1, buff_size
                y_cc(-i) = y_cc(1-i) - (dy(1-i) + dy(-i))/2d0
            END DO
            
            ! Populating the cell-width distribution buffer, at the end of the
            ! coordinate direction, based on desired boundary condition. These
            ! include, in order, ghost-cell extrapolation, symmetry, periodic,
            ! and processor boundary conditions.
            IF(bc_y%end <= -3) THEN
                DO i = 1, buff_size
                    dy(n+i) = dy(n)
                END DO
            ELSEIF(bc_y%end == -2) THEN
                DO i = 1, buff_size
                    dy(n+i) = dy(n-(i-1))
                END DO
            ELSEIF(bc_y%end == -1) THEN
                DO i = 1, buff_size
                    dy(n+i) = dy(i-1)
                END DO
            ELSE
                CALL s_mpi_sendrecv_grid_variables_buffers(2,1)
            END IF
            
            ! Populating the cell-boundary locations buffer, at the end of the
            ! coordinate direction, from buffer of the cell-width distribution
            DO i = 1, buff_size
                y_cb(n+i) = y_cb(n+(i-1)) + dy(n+i)
            END DO
            ! Populating the cell-center locations buffer, at the end of the
            ! coordinate direction, from buffer of the cell-width distribution
            DO i = 1, buff_size
                y_cc(n+i) = y_cc(n+(i-1)) + (dy(n+(i-1)) + dy(n+i))/2d0
            END DO
            
            ! END: Population of Buffers in y-direction ========================
            
            
            ! Population of Buffers in z-direction =============================
            
            ! Populating cell-width distribution buffer, at the beginning of the
            ! coordinate direction, based on the selected boundary condition. In
            ! order, these are the ghost-cell extrapolation, symmetry, periodic,
            ! and processor boundary conditions.
            IF(p == 0) THEN
                RETURN
            ELSEIF(bc_z%beg <= -3) THEN
                DO i = 1, buff_size
                    dz(-i) = dz(0)
                END DO
            ELSEIF(bc_z%beg == -2) THEN
                DO i = 1, buff_size
                    dz(-i) = dz(i-1)
                END DO
            ELSEIF(bc_z%beg == -1) THEN
                DO i = 1, buff_size
                    dz(-i) = dz(p-(i-1))
                END DO
            ELSE
                CALL s_mpi_sendrecv_grid_variables_buffers(3,-1)
            END IF
            
            ! Computing the cell-boundary locations buffer, at the beginning of
            ! the coordinate direction, from the cell-width distribution buffer
            DO i = 1, buff_size
                z_cb(-1-i) = z_cb(-i) - dz(-i)
            END DO
            ! Computing the cell-center locations buffer, at the beginning of
            ! the coordinate direction, from the cell-width distribution buffer
            DO i = 1, buff_size
                z_cc(-i) = z_cc(1-i) - (dz(1-i) + dz(-i))/2d0
            END DO
            
            ! Populating the cell-width distribution buffer, at the end of the
            ! coordinate direction, based on desired boundary condition. These
            ! include, in order, ghost-cell extrapolation, symmetry, periodic,
            ! and processor boundary conditions.
            IF(bc_z%end <= -3) THEN
                DO i = 1, buff_size
                    dz(p+i) = dz(p)
                END DO
            ELSEIF(bc_z%end == -2) THEN
                DO i = 1, buff_size
                    dz(p+i) = dz(p-(i-1))
                END DO
            ELSEIF(bc_z%end == -1) THEN
                DO i = 1, buff_size
                    dz(p+i) = dz(i-1)
                END DO
            ELSE
                CALL s_mpi_sendrecv_grid_variables_buffers(3,1)
            END IF
            
            ! Populating the cell-boundary locations buffer, at the end of the
            ! coordinate direction, from buffer of the cell-width distribution
            DO i = 1, buff_size
                z_cb(p+i) = z_cb(p+(i-1)) + dz(p+i)
            END DO
            ! Populating the cell-center locations buffer, at the end of the
            ! coordinate direction, from buffer of the cell-width distribution
            DO i = 1, buff_size
                z_cc(p+i) = z_cc(p+(i-1)) + (dz(p+(i-1)) + dz(p+i))/2d0
            END DO
            
            ! END: Population of Buffers in z-direction ========================
            
        END SUBROUTINE s_populate_grid_variables_buffers ! ---------------------
        
        
        
        
        
        SUBROUTINE s_account_for_capillary_potential_energy(v_vf) !-------------

            TYPE(scalar_field), DIMENSION(sys_size), INTENT(INOUT) :: v_vf
            REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:,:) :: E_We
            TYPE(bounds_info) :: ix,iy,iz
            TYPE(bounds_info) :: iz1

            ! Placeholders (_ph) for variables 
            REAL(KIND(0d0)) :: rho_ph, gamma_ph, pi_inf_ph
            REAL(KIND(0d0)), DIMENSION(2) :: Re_ph
            REAL(KIND(0d0)), DIMENSION(num_fluids,num_fluids) :: We_ph

            INTEGER :: i,j,k,l

            ix%beg = -buff_size; iy%beg = 0; iz%beg = 0
            IF(n > 0) iy%beg = -buff_size; IF(p > 0) iz%beg = -buff_size
            ix%end = m - ix%beg; iy%end = n - iy%beg; iz%end = p - iz%beg
            ALLOCATE(E_We(ix%beg:ix%end,iy%beg:iy%end,iz%beg:iz%end))

            IF (p > 0) THEN
                iz1%beg = iz%beg; iz1%end = iz%end
            ELSE
                iz1%beg = -1; iz1%end = 1
            END IF

            ! Compute volume fraction gradient magnitude for all fluids
            CALL s_compute_lsq_gradient_curvature(v_vf,grad_x_vf,grad_y_vf,grad_z_vf,norm_vf,kappa_vf)
            
            DO k = iz1%beg+1, iz1%end-1
                DO j = iy%beg+1, iy%end-1
                    DO i = ix%beg+1, ix%end-1
                        E_We(i,j,k) = 0d0
                        
                        CALL s_convert_to_mixture_variables(v_vf, rho_ph, gamma_ph, &
                                                            pi_inf_ph, Re_ph, We_ph, i,j,k)

                        ! Compute capillary potential energy
                        DO l = 1, We_size
                            E_We(i,j,k) = E_We(i,j,k) + &
                                    v_vf(E_idx+We_idx(l,1))%sf(i,j,k) * &
                                    norm_vf(We_idx(l,2))%sf(i,j,k) / We_ph(We_idx(l,1),We_idx(l,2)) + &
                                    v_vf(E_idx+We_idx(l,2))%sf(i,j,k) * &
                                    norm_vf(We_idx(l,1))%sf(i,j,k) / We_ph(We_idx(l,1),We_idx(l,2))
                        END DO
                    END DO
                END DO
            END DO
            
            ! Add capillary potential energy to conservative variable
            v_vf(E_idx)%sf(:,:,:) = v_vf(E_idx)%sf(:,:,:) + E_We(:,:,:)
        
        END SUBROUTINE s_account_for_capillary_potential_energy !---------------

        
        
        
        !> The purpose of this procedure is to initialize the
        !!      values of the internal-energy equations of each phase
        !!      from the mass of each phase, the mixture momentum and
        !!      mixture-total-energy equations.  
        !! @param v_vf conservative variables
        SUBROUTINE s_initialize_internal_energy_equations(v_vf) !---------------


            TYPE(scalar_field), DIMENSION(sys_size), INTENT(INOUT) ::     v_vf
            REAL(KIND(0d0))                                        ::      rho
            REAL(KIND(0d0))                                        :: dyn_pres
            REAL(KIND(0d0))                                        ::     E_We
            REAL(KIND(0d0))                                        ::    gamma
            REAL(KIND(0d0))                                        ::   pi_inf
            REAL(KIND(0d0)), DIMENSION(2)                          ::       Re
            REAL(KIND(0d0)), DIMENSION(num_fluids,num_fluids)      ::       We
            REAL(KIND(0d0))                                        ::     pres

            INTEGER :: i,j,k,l

            DO j = 0, m
                DO k = 0, n
                    DO l = 0, p

                        CALL s_convert_to_mixture_variables(v_vf, rho, gamma, pi_inf, Re, We, j, k, l)

                        dyn_pres = 0d0
                        DO i = mom_idx%beg, mom_idx%end
                            dyn_pres = dyn_pres + 5d-1*v_vf(i)%sf(j,k,l)*v_vf(i)%sf(j,k,l) &
                                / MAX(rho,sgm_eps)
                        END DO

                        E_We = 0d0

                        pres = (v_vf(E_idx)%sf(j,k,l) - dyn_pres - E_We - pi_inf) / gamma

                        DO i = 1, num_fluids
                            v_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l) = v_vf(i+adv_idx%beg-1)%sf(j,k,l) * &
                                (fluid_pp(i)%gamma*pres + fluid_pp(i)%pi_inf)
                        END DO

                    END DO
                END DO
            END DO

        END SUBROUTINE s_initialize_internal_energy_equations !-----------------





        SUBROUTINE s_initialize_start_up_module() !-----------------------------

            TYPE(bounds_info) :: ix, iy, iz


            INTEGER :: i !< Generic loop iterator
            
            IF (We_size > 0 .AND. (We_riemann_flux .OR. We_rhs_flux)) THEN
                ix%beg = -buff_size; iy%beg = 0; iz%beg = 0
                IF(n > 0) iy%beg = -buff_size; IF(p > 0) iz%beg = -buff_size
                ix%end = m - ix%beg; iy%end = n - iy%beg; iz%end = p - iz%beg

                ALLOCATE(grad_x_vf(sys_size))
                ALLOCATE(grad_y_vf(sys_size))
                ALLOCATE(grad_z_vf(sys_size))
                ALLOCATE(norm_vf(1:num_fluids))
                ALLOCATE(kappa_vf(1:num_fluids))

                DO i = 1, crv_size
                    ALLOCATE(grad_x_vf(E_idx+crv_idx(i))%sf(ix%beg:ix%end, &
                                                            iy%beg:iy%end, &
                                                            iz%beg:iz%end ))
                    ALLOCATE(grad_y_vf(E_idx+crv_idx(i))%sf(ix%beg:ix%end, &
                                                            iy%beg:iy%end, &
                                                            iz%beg:iz%end ))
                    ALLOCATE(grad_z_vf(E_idx+crv_idx(i))%sf(ix%beg:ix%end, &
                                                            iy%beg:iy%end, &
                                                            iz%beg:iz%end ))
                    ALLOCATE(norm_vf(crv_idx(i))%sf(ix%beg:ix%end, &
                                                    iy%beg:iy%end, &
                                                    iz%beg:iz%end ))
                    ALLOCATE(kappa_vf(crv_idx(i))%sf(ix%beg:ix%end, &
                                                     iy%beg:iy%end, &
                                                     iz%beg:iz%end ))
                END DO

            END IF

            IF (parallel_io .NEQV. .TRUE.) THEN
                s_read_data_files => s_read_serial_data_files
            ELSE
                s_read_data_files => s_read_parallel_data_files
            END IF

        END SUBROUTINE s_initialize_start_up_module ! --------------------------





        SUBROUTINE s_finalize_start_up_module() ! ------------------------------

            INTEGER :: i !< Generic loop interator
            
            IF (We_size > 0 .AND. (We_riemann_flux .OR. We_rhs_flux)) THEN

                DO i = 1, crv_size
                    DEALLOCATE(grad_x_vf(E_idx+crv_idx(i))%sf)
                    DEALLOCATE(grad_y_vf(E_idx+crv_idx(i))%sf)
                    DEALLOCATE(grad_z_vf(E_idx+crv_idx(i))%sf)
                    DEALLOCATE(norm_vf(crv_idx(i))%sf)
                    DEALLOCATE(kappa_vf(crv_idx(i))%sf)
                END DO

                DEALLOCATE(grad_x_vf,grad_y_vf,grad_z_vf)
                DEALLOCATE(norm_vf,kappa_vf)
            END IF

            s_read_data_files => NULL()

        END SUBROUTINE s_finalize_start_up_module ! ----------------------------





END MODULE m_start_up
