!>
!! @file p_main.f90
!! @brief Contains program p_main
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief The post-process restructures raw unformatted data, outputted by
!!              the simulation, into a formatted database, Silo-HDF5 or Binary,
!!              chosen by the user. The user may also specify which variables to
!!              include in the database. The choices range from any one of the
!!              primitive and conservative variables, as well as quantities that
!!              can be derived from those such as the unadvected volume fraction,
!!              specific heat ratio, liquid stiffness, speed of sound, vorticity
!!              and the numerical Schlieren function.
PROGRAM p_main
    
    
    ! Dependencies =============================================================
    USE m_derived_types         !< Definitions of the derived types
    
    USE m_fftw                  !< Module for FFTW functions
    
    USE m_global_parameters     !< Global parameters for the code
    
    USE m_mpi_proxy             !< Message passing interface (MPI) module proxy
    
    USE m_variables_conversion  !< Subroutines to change the state variables from
                                !! one form to another
    
    USE m_start_up              !< Subroutines that read in and check consistency
                                !! of the user provided inputs and data
    
    USE m_data_input            !< Procedures reading raw simulation data to fill
                                !! the conservative, primitive and grid variables
    
    USE m_data_output           !< Procedures that write the grid and chosen flow
                                !! variable(s) to the formatted database file(s)
    
    USE m_derived_variables     !< Procedures used to compute quantites derived
                                !! from the conservative and primitive variables
    ! ==========================================================================
    
    
    IMPLICIT NONE
    
    

    INTEGER :: t_step !< Iterator for the main time-stepping loop
    
    CHARACTER(LEN = name_len) :: varname !< 
    !! Generic storage for the name(s) of the flow variable(s) that will be added
    !! to the formatted database file(s)
    
    !> @name Generic loop iterator
    !> @{
    INTEGER :: i,j,k
    !> @}

    REAL(KIND(0d0)) :: total_volume !<
    !! Variable for the total volume of the second volume fraction 
    !! to later on track the evolution of the radius of a bubble over time
 
    ! Initialization of the MPI environment
    CALL s_mpi_initialize()

    ! Processor with rank 0 assigns default user input values prior to reading
    ! those in from the input file. Next, the user inputs are read in and their
    ! consistency is checked. The detection of any inconsistencies automatically
    ! leads to the termination of the post-process.
    IF(proc_rank == 0) THEN
        CALL s_assign_default_values_to_user_inputs()
        CALL s_read_input_file()
        CALL s_check_input_file()
    END IF
    
    
    ! Broadcasting the user inputs to all of the processors and performing the
    ! parallel computational domain decomposition. Neither procedure has to be
    ! carried out if the simulation is in fact not truly executed in parallel.
    CALL s_mpi_bcast_user_inputs()
    CALL s_initialize_parallel_io()
    CALL s_mpi_decompose_computational_domain()
    
    ! Computation of parameters, allocation procedures, and/or any other tasks
    ! needed to properly setup the modules
    CALL s_initialize_global_parameters_module()
    IF(num_procs > 1) CALL s_initialize_mpi_proxy_module()
    IF (fourier_decomp) CALL s_initialize_fftw_module()
    CALL s_initialize_variables_conversion_module()
    CALL s_initialize_data_input_module()
    CALL s_initialize_derived_variables_module()
    CALL s_initialize_data_output_module()
    
        ! Associate pointers for serial or paralle I/O
    IF (parallel_io .NEQV. .TRUE.) THEN
        s_read_data_files => s_read_serial_data_files
    ELSE
        s_read_data_files => s_read_parallel_data_files
    END IF
    
    ! Setting the time-step iterator to the first time step to be post-processed
    t_step = t_step_start
   
    ! Time-Marching Loop =======================================================
    DO
        
        ! Populating the grid and conservative variables
        CALL s_read_data_files(t_step)
        ! Populating the buffer regions of the grid variables
        IF(buff_size > 0) THEN
            CALL s_populate_grid_variables_buffer_regions()
        END IF
        
        ! Populating the buffer regions of the conservative variables
        IF(buff_size > 0) THEN
            CALL s_populate_conservative_variables_buffer_regions()
        END IF
       
        ! Converting the conservative variables to the primitive ones
        CALL s_convert_conservative_to_primitive_variables(q_cons_vf, q_prim_vf)
        ! Opening a new formatted database file
        CALL s_open_formatted_database_file(t_step)
        
        
        ! Adding the grid to the formatted database file
        CALL s_write_grid_to_formatted_database_file(t_step)
        
        
        ! Computing centered finite-difference coefficients in x-direction
        IF(omega_wrt(2) .OR. omega_wrt(3) .OR. schlieren_wrt) THEN
           CALL s_compute_finite_difference_coefficients( m, offset_x, x_cc, &
                                                              fd_coeff_x     )
        END IF
        
        ! Computing centered finite-difference coefficients in y-direction
        IF(omega_wrt(1) .OR. omega_wrt(3) .OR. (n > 0 .AND. schlieren_wrt)) THEN
           CALL s_compute_finite_difference_coefficients( n, offset_y, y_cc, &
                                                              fd_coeff_y     )
        END IF
        
        ! Computing centered finite-difference coefficients in z-direction
        IF(omega_wrt(1) .OR. omega_wrt(2) .OR. (p > 0 .AND. schlieren_wrt)) THEN
           CALL s_compute_finite_difference_coefficients( p, offset_z, z_cc, &
                                                              fd_coeff_z     )
        END IF
        
        ! Adding the partial densities to the formatted database file ----------
        IF((model_eqns == 2) .OR. (model_eqns == 3) .or. (model_eqns == 4)) THEN
            DO i = 1, num_fluids
                IF (alpha_rho_wrt(i) .OR. (cons_vars_wrt .OR. prim_vars_wrt)) THEN
                    
                    q_sf = q_cons_vf(i)%sf( -offset_x%beg : m + offset_x%end, &
                                            -offset_y%beg : n + offset_y%end, &
                                            -offset_z%beg : p + offset_z%end  )

                    IF (fourier_decomp) THEN
                        dft_q_sf(:,:,:) = q_sf(:,:,:)
                        DO j = fourier_modes%beg, fourier_modes%end
                            q_sf(:,:,:) = dft_q_sf(:,:,:)
                            CALL s_apply_fourier_decomposition(q_sf,j)
                            WRITE(varname, '(A,I0,A,I0)') 'alpha_rho', i, '_', j
                            CALL s_write_variable_to_formatted_database_file(varname,t_step)
                        END DO
                    ELSE
                        IF (model_eqns .NE. 4) THEN
                            WRITE(varname, '(A,I0)') 'alpha_rho', i
                        ELSE
                            WRITE(varname, '(A,I0)') 'rho', i
                        END IF
                        CALL s_write_variable_to_formatted_database_file(varname,t_step)
                    END IF
                    
                    varname(:) = ' '
                    
                END IF
            END DO
        END IF
        ! ----------------------------------------------------------------------
        
        
        ! Adding the density to the formatted database file --------------------
        IF(                           rho_wrt                             &
                                        .OR.                              &
              (model_eqns == 1 .AND. (cons_vars_wrt .OR. prim_vars_wrt))  ) THEN
            
            q_sf = rho_sf( -offset_x%beg : m + offset_x%end, &
                           -offset_y%beg : n + offset_y%end, &
                           -offset_z%beg : p + offset_z%end  )
            
            IF (fourier_decomp) THEN
                dft_q_sf(:,:,:) = q_sf(:,:,:)
                DO j = fourier_modes%beg, fourier_modes%end
                    q_sf(:,:,:) = dft_q_sf(:,:,:)
                    CALL s_apply_fourier_decomposition(q_sf,j)
                    WRITE(varname, '(A,I0)') 'rho_', j
                    CALL s_write_variable_to_formatted_database_file(varname,t_step)
                END DO
            ELSE
                WRITE(varname, '(A)') 'rho'
                CALL s_write_variable_to_formatted_database_file(varname,t_step)
            END IF
            
            varname(:) = ' '
            
        END IF
        ! ----------------------------------------------------------------------
        
        
        ! Adding the momentum to the formatted database file -------------------
        DO i = 1, E_idx-mom_idx%beg
            IF(mom_wrt(i) .OR. cons_vars_wrt) THEN
                
                q_sf = q_cons_vf(i+cont_idx%end)%sf(            &
                              -offset_x%beg : m + offset_x%end, &
                              -offset_y%beg : n + offset_y%end, &
                              -offset_z%beg : p + offset_z%end  )
                
                IF (fourier_decomp) THEN
                    dft_q_sf(:,:,:) = q_sf(:,:,:)
                    DO j = fourier_modes%beg, fourier_modes%end
                        q_sf(:,:,:) = dft_q_sf(:,:,:)
                        CALL s_apply_fourier_decomposition(q_sf,j)
                        WRITE(varname, '(A,I0,A,I0)') 'mom', i, '_', j
                        CALL s_write_variable_to_formatted_database_file(varname,t_step)
                    END DO
                ELSE
                    WRITE(varname, '(A,I0)') 'mom', i
                    CALL s_write_variable_to_formatted_database_file(varname,t_step)
                END IF
                
                varname(:) = ' '
                
            END IF
        END DO
        ! ----------------------------------------------------------------------
        
        
        ! Adding the velocity to the formatted database file -------------------
        DO i = 1, E_idx-mom_idx%beg
            IF(vel_wrt(i) .OR. prim_vars_wrt) THEN
                
                q_sf = q_prim_vf(i+cont_idx%end)%sf(            &
                              -offset_x%beg : m + offset_x%end, &
                              -offset_y%beg : n + offset_y%end, &
                              -offset_z%beg : p + offset_z%end  )
                
                IF (fourier_decomp) THEN
                    dft_q_sf(:,:,:) = q_sf(:,:,:)
                    DO j = fourier_modes%beg, fourier_modes%end
                        IF (ANY(j == (/1,2,3,5,7,9/))) THEN
                        q_sf(:,:,:) = dft_q_sf(:,:,:)
                        CALL s_apply_fourier_decomposition(q_sf,j)
                        WRITE(varname, '(A,I0,A,I0)') 'vel', i, '_', j
                        CALL s_write_variable_to_formatted_database_file(varname,t_step)
                        END IF
                    END DO
                ELSE
                    WRITE(varname, '(A,I0)') 'vel', i
                    CALL s_write_variable_to_formatted_database_file(varname,t_step)
                END IF
                
                varname(:) = ' '
                
            END IF
        END DO
        ! ----------------------------------------------------------------------

        ! Adding the flux limiter function to the formatted database file
        DO i = 1, E_idx-mom_idx%beg
            IF (flux_wrt(i)) THEN
    
                CALL s_derive_flux_limiter(i, q_prim_vf, q_sf)
    
                IF (fourier_decomp) THEN
                    dft_q_sf(:,:,:) = q_sf(:,:,:)
                    DO j = fourier_modes%beg, fourier_modes%end
                        q_sf(:,:,:) = dft_q_sf(:,:,:)
                        CALL s_apply_fourier_decomposition(q_sf,j)
                        WRITE (varname, '(A,I0,A,I0)') 'flux', i, '_', j
                        CALL s_write_variable_to_formatted_database_file(varname,t_step)
                    END DO
                ELSE
                    WRITE (varname, '(A,I0)') 'flux', i
                    CALL s_write_variable_to_formatted_database_file(varname,t_step)
                END IF
    
                varname(:) = ' '
            END IF
        END DO
        ! ----------------------------------------------------------------------

        ! Adding the curvature to the formatted database file
        DO i = 1, num_fluids
            IF (kappa_wrt(i)) THEN

                CALL s_derive_curvature(i, q_prim_vf, q_sf)

                IF (fourier_decomp) THEN
                    dft_q_sf(:,:,:) = q_sf(:,:,:)
                    DO j = fourier_modes%beg, fourier_modes%end
                        q_sf(:,:,:) = dft_q_sf(:,:,:)
                        CALL s_apply_fourier_decomposition(q_sf,j)
                        WRITE(varname, '(A,I0,A,I0)') 'kappa', i, '_', j
                        CALL s_write_variable_to_formatted_database_file(varname,t_step)
                    END DO
                ELSE
                    WRITE(varname, '(A,I0)') 'kappa', i
                    CALL s_write_variable_to_formatted_database_file(varname,t_step)
                END IF

                varname(:) = ' '
            END IF
        END DO
        ! ----------------------------------------------------------------------
        
        ! Adding the energy to the formatted database file ---------------------
        IF(E_wrt .OR. cons_vars_wrt) THEN
            
            q_sf = q_cons_vf(E_idx)%sf( -offset_x%beg : m + offset_x%end, &
                                        -offset_y%beg : n + offset_y%end, &
                                        -offset_z%beg : p + offset_z%end  )
            
            IF (fourier_decomp) THEN
                dft_q_sf(:,:,:) = q_sf(:,:,:)
                DO j = fourier_modes%beg, fourier_modes%end
                    q_sf(:,:,:) = dft_q_sf(:,:,:)
                    CALL s_apply_fourier_decomposition(q_sf,j)
                    WRITE(varname, '(A,I0)') 'E_', j
                    CALL s_write_variable_to_formatted_database_file(varname,t_step)
                END DO
            ELSE
                WRITE(varname, '(A)') 'E'
                CALL s_write_variable_to_formatted_database_file(varname,t_step)
            END IF
            
            varname(:) = ' '
            
        END IF
        ! ----------------------------------------------------------------------
        
        
        ! Adding the pressure to the formatted database file -------------------
        IF(pres_wrt .OR. prim_vars_wrt) THEN
            
            q_sf = q_prim_vf(E_idx)%sf( -offset_x%beg : m + offset_x%end, &
                                        -offset_y%beg : n + offset_y%end, &
                                        -offset_z%beg : p + offset_z%end  )
            
            IF (fourier_decomp) THEN
                dft_q_sf(:,:,:) = q_sf(:,:,:)
                DO j = fourier_modes%beg, fourier_modes%end
                    q_sf(:,:,:) = dft_q_sf(:,:,:)
                    CALL s_apply_fourier_decomposition(q_sf,j)
                    WRITE(varname, '(A,I0)') 'pres_', j
                    CALL s_write_variable_to_formatted_database_file(varname,t_step)
                END DO
            ELSE
                WRITE(varname, '(A)') 'pres'
                CALL s_write_variable_to_formatted_database_file(varname,t_step)
            END IF
            
            varname(:) = ' '
            
        END IF
        ! ----------------------------------------------------------------------
        
        ! Adding the volume fraction(s) to the formatted database file ---------
        IF( ( (model_eqns == 2) .AND. (bubbles .NEQV. .TRUE.)) &
                .OR. (model_eqns == 3) &
                ) THEN
            
            DO i = 1, num_fluids-1
                IF(alpha_wrt(i) .OR. (cons_vars_wrt .OR. prim_vars_wrt)) THEN
                    
                    q_sf = q_cons_vf(i+E_idx)%sf(                 &
                                -offset_x%beg : m + offset_x%end, &
                                -offset_y%beg : n + offset_y%end, &
                                -offset_z%beg : p + offset_z%end  )
                    
                    IF (fourier_decomp) THEN
                        dft_q_sf(:,:,:) = q_sf(:,:,:)
                        DO j = fourier_modes%beg, fourier_modes%end
                            q_sf(:,:,:) = dft_q_sf(:,:,:)
                            CALL s_apply_fourier_decomposition(q_sf,j)
                            WRITE(varname, '(A,I0,A,I0)') 'alpha', i, '_', j
                            CALL s_write_variable_to_formatted_database_file(varname,t_step)
                        END DO
                    ELSE
                        WRITE(varname, '(A,I0)') 'alpha', i
                        CALL s_write_variable_to_formatted_database_file(varname,t_step)
                    END IF
                    
                    varname(:) = ' '
                    
                END IF
            END DO
            
            IF(                   alpha_wrt(num_fluids)                   &
                                           .OR.                           &
                  (adv_alphan .AND. (cons_vars_wrt .OR. prim_vars_wrt))   ) THEN
                
                IF(adv_alphan .NEQV. .TRUE.) THEN
                    CALL s_derive_unadvected_volume_fraction(q_cons_vf, q_sf)
                ELSE
                    q_sf = q_cons_vf(adv_idx%end)%sf(              &
                                 -offset_x%beg : m + offset_x%end, &
                                 -offset_y%beg : n + offset_y%end, &
                                 -offset_z%beg : p + offset_z%end  )
                END IF
                
                IF (fourier_decomp) THEN
                    dft_q_sf(:,:,:) = q_sf(:,:,:)
                    DO j = fourier_modes%beg, fourier_modes%end
                        q_sf(:,:,:) = dft_q_sf(:,:,:)
                        CALL s_apply_fourier_decomposition(q_sf,j)
                        WRITE(varname, '(A,I0,A,I0)') 'alpha', num_fluids, '_', j
                        CALL s_write_variable_to_formatted_database_file(varname,t_step)
                    END DO
                ELSE
                    WRITE(varname, '(A,I0)') 'alpha', num_fluids
                    CALL s_write_variable_to_formatted_database_file(varname,t_step)
                END IF
                
                varname(:) = ' '
                
            END IF
            
        END IF
        ! ----------------------------------------------------------------------
        
        ! Adding specific heat ratio function to formatted database file -------
        IF(                           gamma_wrt                           &
                                        .OR.                              &
              (model_eqns == 1 .AND. (cons_vars_wrt .OR. prim_vars_wrt))  ) THEN
            
            q_sf = gamma_sf( -offset_x%beg : m + offset_x%end, &
                             -offset_y%beg : n + offset_y%end, &
                             -offset_z%beg : p + offset_z%end  )
            
            IF (fourier_decomp) THEN
                dft_q_sf(:,:,:) = q_sf(:,:,:)
                DO j = fourier_modes%beg, fourier_modes%end
                    q_sf(:,:,:) = dft_q_sf(:,:,:)
                    CALL s_apply_fourier_decomposition(q_sf,j)
                    WRITE(varname, '(A,I0)') 'gamma_', j
                    CALL s_write_variable_to_formatted_database_file(varname,t_step)
                END DO
            ELSE
                WRITE(varname, '(A)') 'gamma'
                CALL s_write_variable_to_formatted_database_file(varname,t_step)
            END IF
            
            varname(:) = ' '
            
        END IF
        ! ----------------------------------------------------------------------
        
        
        ! Adding the specific heat ratio to the formatted database file --------
        IF(heat_ratio_wrt) THEN
            
            CALL s_derive_specific_heat_ratio(gamma_sf, q_sf)
            
            IF (fourier_decomp) THEN
                dft_q_sf(:,:,:) = q_sf(:,:,:)
                DO j = fourier_modes%beg, fourier_modes%end
                    q_sf(:,:,:) = dft_q_sf(:,:,:)
                    CALL s_apply_fourier_decomposition(q_sf,j)
                    WRITE(varname, '(A,I0)') 'heat_ratio_', j
                    CALL s_write_variable_to_formatted_database_file(varname,t_step)
                END DO
            ELSE
                WRITE(varname, '(A)') 'heat_ratio'
                CALL s_write_variable_to_formatted_database_file(varname,t_step)
            END IF
            
            varname(:) = ' '
            
        END IF
        ! ----------------------------------------------------------------------
        
        
        ! Adding liquid stiffness function to formatted database file ----------
        IF(                          pi_inf_wrt                           &
                                        .OR.                              &
              (model_eqns == 1 .AND. (cons_vars_wrt .OR. prim_vars_wrt))  ) THEN
            
            q_sf = pi_inf_sf( -offset_x%beg : m + offset_x%end, &
                              -offset_y%beg : n + offset_y%end, &
                              -offset_z%beg : p + offset_z%end  )
            
            IF (fourier_decomp) THEN
                dft_q_sf(:,:,:) = q_sf(:,:,:)
                DO j = fourier_modes%beg, fourier_modes%end
                    q_sf(:,:,:) = dft_q_sf(:,:,:)
                    CALL s_apply_fourier_decomposition(q_sf,j)
                    WRITE(varname, '(A,I0)') 'pi_inf_', j
                    CALL s_write_variable_to_formatted_database_file(varname,t_step)
                END DO
            ELSE
                WRITE(varname, '(A)') 'pi_inf'
                CALL s_write_variable_to_formatted_database_file(varname,t_step)
            END IF
            
            varname(:) = ' '
            
        END IF
        ! ----------------------------------------------------------------------
        
        
        ! Adding the liquid stiffness to the formatted database file -----------
        IF(pres_inf_wrt) THEN
            
            CALL s_derive_liquid_stiffness(gamma_sf, pi_inf_sf, q_sf)
            
            IF (fourier_decomp) THEN
                dft_q_sf(:,:,:) = q_sf(:,:,:)
                DO j = fourier_modes%beg, fourier_modes%end
                    q_sf(:,:,:) = dft_q_sf(:,:,:)
                    CALL s_apply_fourier_decomposition(q_sf,j)
                    WRITE(varname, '(A,I0)') 'pres_inf_', j
                    CALL s_write_variable_to_formatted_database_file(varname,t_step)
                END DO
            ELSE
                WRITE(varname, '(A)') 'pres_inf'
                CALL s_write_variable_to_formatted_database_file(varname,t_step)
            END IF
            
            varname(:) = ' '
            
        END IF
        ! ----------------------------------------------------------------------
        
        
        ! Adding the sound speed to the formatted database file ----------------
        IF(c_wrt) THEN
            
            CALL s_derive_sound_speed( q_prim_vf, rho_sf, gamma_sf, &
                                       pi_inf_sf, q_sf )
            
            IF (fourier_decomp) THEN
                dft_q_sf(:,:,:) = q_sf(:,:,:)
                DO j = fourier_modes%beg, fourier_modes%end
                    q_sf(:,:,:) = dft_q_sf(:,:,:)
                    CALL s_apply_fourier_decomposition(q_sf,j)
                    WRITE(varname, '(A,I0)') 'c_', j
                    CALL s_write_variable_to_formatted_database_file(varname,t_step)
                END DO
            ELSE
                WRITE(varname, '(A)') 'c'
                CALL s_write_variable_to_formatted_database_file(varname,t_step)
            END IF
            
            varname(:) = ' '
            
        END IF
        ! ----------------------------------------------------------------------
        
        
        ! Adding the vorticity to the formatted database file ------------------
        IF (p > 0) THEN
            DO i = 1, E_idx-mom_idx%beg
                IF (omega_wrt(i)) THEN

                    CALL s_derive_vorticity_component(i, q_prim_vf, q_sf)

                    IF (fourier_decomp) THEN
                        dft_q_sf(:,:,:) = q_sf(:,:,:)
                        DO j = fourier_modes%beg, fourier_modes%end
                            q_sf(:,:,:) = dft_q_sf(:,:,:)
                            CALL s_apply_fourier_decomposition(q_sf,j)
                            WRITE(varname, '(A,I0,A,I0)') 'omega', i, '_', j
                            CALL s_write_variable_to_formatted_database_file(varname,t_step)
                        END DO
                    ELSE
                        WRITE(varname, '(A,I0)') 'omega', i
                        CALL s_write_variable_to_formatted_database_file(varname,t_step)
                    END IF

                    varname(:) = ' '
                END IF
            END DO
        ELSEIF (n > 0) THEN
            DO i = 1, E_idx-cont_idx%end
                IF (omega_wrt(i)) THEN

                    CALL s_derive_vorticity_component(i, q_prim_vf, q_sf)

                    IF (fourier_decomp) THEN
                        dft_q_sf(:,:,:) = q_sf(:,:,:)
                        DO j = fourier_modes%beg, fourier_modes%end
                            q_sf(:,:,:) = dft_q_sf(:,:,:)
                            CALL s_apply_fourier_decomposition(q_sf,j)
                            WRITE(varname, '(A,I0,A,I0)') 'omega', i, '_', j
                            CALL s_write_variable_to_formatted_database_file(varname,t_step)
                        END DO
                    ELSE
                        WRITE(varname, '(A,I0)') 'omega', i
                        CALL s_write_variable_to_formatted_database_file(varname,t_step)
                    END IF

                    varname(:) = ' '
                END IF
            END DO
        END IF
        ! ----------------------------------------------------------------------
        
        
        ! Adding numerical Schlieren function to formatted database file -------
        IF(schlieren_wrt) THEN
            
            CALL s_derive_numerical_schlieren_function(q_cons_vf, rho_sf, q_sf)
            
            IF (fourier_decomp) THEN
                dft_q_sf(:,:,:) = q_sf(:,:,:)
                DO j = fourier_modes%beg, fourier_modes%end
                    q_sf(:,:,:) = dft_q_sf(:,:,:)
                    CALL s_apply_fourier_decomposition(q_sf,j)
                    WRITE(varname, '(A,I0)') 'schlieren_', j
                    CALL s_write_variable_to_formatted_database_file(varname,t_step)
                END DO
            ELSE
                WRITE(varname, '(A)') 'schlieren'
                CALL s_write_variable_to_formatted_database_file(varname,t_step)
            END IF
            
            varname(:) = ' '
            
        END IF
        ! ----------------------------------------------------------------------
        
        ! Adding the volume fraction(s) to the formatted database file ---------
        IF (bubbles) THEN
            DO i = adv_idx%beg,adv_idx%end
                q_sf = q_cons_vf(i)%sf(                 &
                                -offset_x%beg : m + offset_x%end, &
                                -offset_y%beg : n + offset_y%end, &
                                -offset_z%beg : p + offset_z%end  )
                    
                IF (fourier_decomp) THEN
                    dft_q_sf(:,:,:) = q_sf(:,:,:)
                    DO j = fourier_modes%beg, fourier_modes%end
                        q_sf(:,:,:) = dft_q_sf(:,:,:)
                        CALL s_apply_fourier_decomposition(q_sf,j)
                        WRITE(varname, '(A,I0,A,I0)') 'alpha', i-E_idx, '_', j
                        CALL s_write_variable_to_formatted_database_file(varname,t_step)
                    END DO
                ELSE
                    WRITE(varname, '(A,I0)') 'alpha', i-E_idx
                    CALL s_write_variable_to_formatted_database_file(varname,t_step)
                END IF
                varname(:) = ' '
            END DO  
        END IF
 
        ! Adding the bubble variables  to the formatted database file ---------
        IF (bubbles) THEN
                !nR
                DO i = 1,nb
                        q_sf = q_cons_vf(bub_idx%rs(i))%sf(       &
                                -offset_x%beg : m + offset_x%end, &
                                -offset_y%beg : n + offset_y%end, &
                                -offset_z%beg : p + offset_z%end  )
                        WRITE(varname, '(A,I3.3)') 'nR', i
                        CALL s_write_variable_to_formatted_database_file(varname,t_step)
                        varname(:) = ' '
                END DO

                !nRdot
                DO i = 1,nb
                        q_sf = q_cons_vf(bub_idx%vs(i))%sf(       &
                                -offset_x%beg : m + offset_x%end, &
                                -offset_y%beg : n + offset_y%end, &
                                -offset_z%beg : p + offset_z%end  )
                        WRITE(varname, '(A,I3.3)') 'nV', i
                        CALL s_write_variable_to_formatted_database_file(varname,t_step)
                        varname(:) = ' '
                END DO
                IF (polytropic .NEQV. .TRUE.) THEN
                    !nP
                    DO i = 1,nb
                            q_sf = q_cons_vf(bub_idx%ps(i))%sf(       &
                                    -offset_x%beg : m + offset_x%end, &
                                    -offset_y%beg : n + offset_y%end, &
                                    -offset_z%beg : p + offset_z%end  )
                            WRITE(varname, '(A,I3.3)') 'nP', i
                            CALL s_write_variable_to_formatted_database_file(varname,t_step)
                            varname(:) = ' '
                    END DO

                    !nM
                    DO i = 1,nb
                            q_sf = q_cons_vf(bub_idx%ms(i))%sf(       &
                                    -offset_x%beg : m + offset_x%end, &
                                    -offset_y%beg : n + offset_y%end, &
                                    -offset_z%beg : p + offset_z%end  )
                            WRITE(varname, '(A,I3.3)') 'nM', i
                            CALL s_write_variable_to_formatted_database_file(varname,t_step)
                            varname(:) = ' '
                    END DO
                END IF
        END IF
        
        ! Closing the formatted database file
        CALL s_close_formatted_database_file()
        
        ! Modifies the time-step iterator so that it may reach the final time-
        ! step to be post-processed, in the case that this one is not originally
        ! attainable through constant incrementation from the first time-step.
        ! This modification is performed upon reaching the final time-step. In
        ! case that it is not needed, the post-processor is done and may exit.
        IF((t_step_stop-t_step) < t_step_save .AND. t_step_stop /= t_step) THEN
            t_step = t_step_stop - t_step_save
        ELSEIF(t_step == t_step_stop) THEN
            EXIT
        END IF
        
        ! Incrementing time-step iterator to next time-step to be post-processed
        t_step = t_step + t_step_save
        
    END DO
    ! END: Time-Marching Loop ==================================================
    
    CLOSE(11)

    ! Disassociate pointers for serial and parallel I/O
    s_read_data_files => NULL()
    
    ! Deallocation procedures for the modules
    CALL s_finalize_data_output_module()
    CALL s_finalize_derived_variables_module()
    CALL s_finalize_data_input_module()
    CALL s_finalize_variables_conversion_module()
    IF (fourier_decomp) CALL s_finalize_fftw_module()
    IF(num_procs > 1) CALL s_finalize_mpi_proxy_module()
    CALL s_finalize_global_parameters_module()
    
    
    ! Finalizing the MPI environment
    CALL s_mpi_finalize()
    
    
END PROGRAM p_main
