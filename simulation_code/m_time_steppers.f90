! MFC v3.0 - Simulation Code: m_time_steppers.f90
! Description: The following module features a variety of time-stepping schemes.
!              Currently, it includes the following Runge-Kutta (RK) algorithms:
!                                     1) 1st Order TVD RK
!                                     2) 2nd Order TVD RK
!                                     3) 3rd Order TVD RK
!                                     4) 4th Order RK
!                                     5) 5th Order RK
!              where TVD designates a total-variation-diminishing time-stepper.
! Author: Vedran Coralic
! Date: 06/08/12


MODULE m_time_steppers
    
    
    ! Dependencies =============================================================
    USE m_derived_types        ! Definitions of the derived types
    
    USE m_global_parameters    ! Definitions of the global parameters
    
    USE m_fftw                 ! Module for FFTW functions
    
    USE m_rhs                  ! Right-hand-side (RHS) evaluation procedures
    
    USE m_data_output          ! Run-time info & solution data output procedures

    USE m_bubbles               

    USE m_mpi_proxy            ! Message passing interface (MPI) module proxy
    ! ==========================================================================
    
    
    IMPLICIT NONE
    
    
    ! Cell-average conservative variables at each time-stage (TS)
    TYPE(vector_field), ALLOCATABLE, DIMENSION(:) :: q_cons_ts
    
    ! Cell-average primitive variables at the current time-stage
    TYPE(scalar_field), PRIVATE, ALLOCATABLE, DIMENSION(:) :: q_prim_vf
    
    ! Cell-average RHS variables at the current time-stage
    TYPE(scalar_field), ALLOCATABLE, DIMENSION(:) :: rhs_vf
    
    ! Cell-average primitive variables at consecutive TIMESTEPS
    TYPE(vector_field), ALLOCATABLE, DIMENSION(:) :: q_prim_ts

    ! Number of time stages in the time-stepping scheme
    INTEGER, PRIVATE :: num_ts
    
    
    CONTAINS
        
        
        
        
        
        SUBROUTINE s_initialize_time_steppers_module() ! -----------------------
        ! Description: The computation of parameters, the allocation of memory,
        !              the association of pointers and/or the execution of any
        !              other procedures that are necessary to setup the module.
            
            
            ! Indical bounds in the x-, y- and z-directions
            TYPE(bounds_info) :: ix,iy,iz
            
            ! Generic loop iterators
            INTEGER :: i,j
            
            
            ! Setting number of time-stages for selected time-stepping scheme
            IF(time_stepper == 1) THEN
                num_ts = 1
            ELSEIF(ANY(time_stepper == (/2,3/))) THEN
                num_ts = 2
            ELSEIF(time_stepper == 4) THEN
                num_ts = 3
            ELSE
                num_ts = 6
            END IF
            
            
            ! Setting the indical bounds in the x-, y- and z-directions
            ix%beg = -buff_size; ix%end = m + buff_size
            
            IF(n > 0) THEN
                
                iy%beg = -buff_size; iy%end = n + buff_size
                
                IF(p > 0) THEN
                    iz%beg = -buff_size; iz%end = p + buff_size
                ELSE
                    iz%beg = 0; iz%end = 0
                END IF
                
            ELSE
                
                iy%beg = 0; iy%end = 0
                iz%beg = 0; iz%end = 0
                
            END IF
            
            
            ! Allocating the cell-average conservative variables
            ALLOCATE(q_cons_ts(1:num_ts))
            
            DO i = 1, num_ts
                ALLOCATE(q_cons_ts(i)%vf(1:sys_size))
            END DO
            
            DO i = 1, num_ts
                DO j = 1, sys_size
                    ALLOCATE(q_cons_ts(i)%vf(j)%sf( ix%beg:ix%end, &
                                                    iy%beg:iy%end, &
                                                    iz%beg:iz%end ))
                END DO
            END DO
            
            
            ! Allocating the cell-average primitive ts variables
            IF (ANY(com_wrt) .OR. ANY(cb_wrt) .OR. probe_wrt) THEN
                ALLOCATE(q_prim_ts(0:3))

                DO i = 0, 3
                    ALLOCATE(q_prim_ts(i)%vf(1:sys_size))
                END DO

                DO i = 0, 3
                    DO j = 1, sys_size
                        ALLOCATE(q_prim_ts(i)%vf(j)%sf( ix%beg:ix%end, &
                                        iy%beg:iy%end, &
                                        iz%beg:iz%end ))
                    END DO
                END DO
            END IF


            ! Allocating the cell-average primitive variables
            ALLOCATE(q_prim_vf(1:sys_size))
            
            DO i = mom_idx%beg, E_idx
                ALLOCATE(q_prim_vf(i)%sf( ix%beg:ix%end, &
                                          iy%beg:iy%end, &
                                          iz%beg:iz%end ))
            END DO

            ! SHB: added for bubble variables
            if (bubbles) then
                DO i = bub_idx%beg,sys_size
                    ALLOCATE(q_prim_vf(i)%sf( ix%beg:ix%end, &
                                          iy%beg:iy%end, &
                                          iz%beg:iz%end ))
                END DO
            end if

            IF (model_eqns == 3) THEN
                DO i = internalEnergies_idx%beg, internalEnergies_idx%end
                    ALLOCATE(q_prim_vf(i)%sf( ix%beg:ix%end, &
                                              iy%beg:iy%end, &
                                              iz%beg:iz%end ))
                END DO
            END IF
            
            
            ! Allocating the cell-average RHS variables
            ALLOCATE(rhs_vf(1:sys_size))
            
            DO i = 1, sys_size
                ALLOCATE(rhs_vf(i)%sf(0:m,0:n,0:p))
            END DO
            
            
            ! Opening and writing the header of the run-time information file
            IF(proc_rank == 0 .AND. run_time_info) THEN
                CALL s_open_run_time_information_file()
            END IF


        END SUBROUTINE s_initialize_time_steppers_module ! ---------------------
        
        
        
        
        
        SUBROUTINE s_1st_order_tvd_rk(t_step) ! --------------------------------
        ! Description: 1st order TVD RK time-stepping algorithm
            
            
            ! Current time-step
            INTEGER, INTENT(IN) :: t_step
            
            ! Generic loop iterator
            INTEGER :: i
!            INTEGER :: j
           
            ! Stage 1 of 1 =====================================================
            DO i = 1, cont_idx%end
                q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
            END DO

            !SHB edits
            IF (adv_alphan) THEN
                DO i = adv_idx%beg, adv_idx%end
                        q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
                END DO
            ELSE 
                DO i = adv_idx%beg, sys_size
                    q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
                END DO
            END IF

            print*, 'before rhs'
            
            CALL s_compute_rhs(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, t_step)
           

            print '(a)', 'SHB: end of rhs in time steppers'

            IF(run_time_info) THEN
                CALL s_write_run_time_information(q_prim_vf, t_step)
            END IF
            print '(a)', 'SHB: end of run_time info'

            IF (ANY(com_wrt) .OR. ANY(cb_wrt) .OR. probe_wrt) THEN
                CALL s_time_step_cycling(t_step)
            END IF
           
            
            IF(t_step == t_step_stop) RETURN
           
            DO i = 1, sys_size
                q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) = &
                               q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) &
                             + dt*rhs_vf(i)%sf
            END DO

            !print '(a)', 'SHB: end of time advancing'
            IF (grid_geometry == 3) CALL s_apply_fourier_filter(q_cons_ts(1)%vf)

            IF (model_eqns == 3) CALL s_pressure_relaxation_procedure(q_cons_ts(1)%vf)

        ! Chunk of code to output RHS of advection equation when checking changes
        ! using the analytical initial condition

        ! PRINT '(A,I0)','This is m: ',m
        ! PRINT '(A)','This is alt_soundspeed: ',alt_soundspeed
        ! i = 2d0/5d0*(m+1) + (m+1)/25d0/2d0 - 5d-1
        ! j = i
        ! PRINT '(A,E15.10,A,E15.10)', 'The location of the cell center is ',x_cc(i),' , ',y_cc(j)
        ! PRINT '(A,E20.10)','This is rhs for adv: ',rhs_vf(adv_idx%end)%sf(i,j,0)

        ! i = m/2
        ! j = n/2
        ! PRINT '(A,E15.10,A,E15.10)', 'The location of the cell center is ',x_cc(i),' , ',y_cc(j)
        ! PRINT '(A,E20.10)','This is rhs for adv: ',rhs_vf(adv_idx%end)%sf(i,j,0)

        ! i = 4d0/5d0*(m+1) + (m+1)/25d0/2d0 - 5d-1
        ! j = i
        ! PRINT '(A,E15.10,A,E15.10)', 'The location of the cell center is ',x_cc(i),' , ',y_cc(j)
        ! PRINT '(A,E20.10)','This is rhs for adv: ',rhs_vf(adv_idx%end)%sf(i,j,0)
        ! CALL s_mpi_abort()
            
            DO i = 1, cont_idx%end
                q_prim_vf(i)%sf => NULL()
            END DO
           
            !SHB edits
            if (adv_alphan) then
                DO i = adv_idx%beg, adv_idx%end
                    q_prim_vf(i)%sf => NULL()
                END DO
            else
                DO i = adv_idx%beg, sys_size ! adv_idx%end
                    q_prim_vf(i)%sf => NULL()
                END DO 
            end if
            ! ==================================================================
           
            !print '(a)', 'SHB: end of time step'

            
        END SUBROUTINE s_1st_order_tvd_rk ! ------------------------------------
        
        
        
        
        
        SUBROUTINE s_2nd_order_tvd_rk(t_step) ! --------------------------------
        ! Description: 2nd order TVD RK time-stepping algorithm
            
            
            ! Current time-step
            INTEGER, INTENT(IN) :: t_step
            
            ! Generic loop iterator
            INTEGER :: i
            
            
            ! Stage 1 of 2 =====================================================
            DO i = 1, cont_idx%end
                q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
            END DO
            
            DO i = adv_idx%beg, adv_idx%end
                q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
            END DO
            
            CALL s_compute_rhs(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, t_step)
            
            IF(run_time_info) THEN
                CALL s_write_run_time_information(q_prim_vf, t_step)
            END IF
            
            IF (ANY(com_wrt) .OR. ANY(cb_wrt) .OR. probe_wrt) THEN
                CALL s_time_step_cycling(t_step)
            END IF
            
            IF(t_step == t_step_stop) RETURN
            
            DO i = 1, sys_size
                q_cons_ts(2)%vf(i)%sf(0:m,0:n,0:p) = &
                               q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) &
                             + dt*rhs_vf(i)%sf
            END DO
            
            IF (grid_geometry == 3) CALL s_apply_fourier_filter(q_cons_ts(2)%vf)

            IF (model_eqns == 3) CALL s_pressure_relaxation_procedure(q_cons_ts(2)%vf)
            ! ==================================================================
            
            
            ! Stage 2 of 2 =====================================================
            DO i = 1, cont_idx%end
                q_prim_vf(i)%sf => q_cons_ts(2)%vf(i)%sf
            END DO
            
            DO i = adv_idx%beg, adv_idx%end
                q_prim_vf(i)%sf => q_cons_ts(2)%vf(i)%sf
            END DO
            
            CALL s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, t_step)
            
            DO i = 1, sys_size
                q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) = &
                             ( q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) &
                             + q_cons_ts(2)%vf(i)%sf(0:m,0:n,0:p) &
                             + dt*rhs_vf(i)%sf ) / 2d0
            END DO
            
            IF (grid_geometry == 3) CALL s_apply_fourier_filter(q_cons_ts(1)%vf)

            IF (model_eqns == 3) CALL s_pressure_relaxation_procedure(q_cons_ts(1)%vf)
            
            DO i = 1, cont_idx%end
                q_prim_vf(i)%sf => NULL()
            END DO
            
            DO i = adv_idx%beg, adv_idx%end
                q_prim_vf(i)%sf => NULL()
            END DO
            ! ==================================================================
            
            
        END SUBROUTINE s_2nd_order_tvd_rk ! ------------------------------------
        
        
        
        
        
        SUBROUTINE s_3rd_order_tvd_rk(t_step) ! --------------------------------
        ! Description: 3rd order TVD RK time-stepping algorithm
            
            
            ! Current time-step
            INTEGER, INTENT(IN) :: t_step
            
            ! Generic loop iterator
            INTEGER :: i,j
            
            ! Stage 1 of 3 =====================================================
            DO i = 1, cont_idx%end
                q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
            END DO
           
            DO i = adv_idx%beg, adv_idx%end
                q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
            END DO
            
            CALL s_compute_rhs(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, t_step)
            
            IF(run_time_info) THEN
                CALL s_write_run_time_information(q_prim_vf, t_step)
            END IF
            
            IF (ANY(com_wrt) .OR. ANY(cb_wrt) .OR. probe_wrt) THEN
                CALL s_time_step_cycling(t_step)
            END IF
            
            IF(t_step == t_step_stop) RETURN
            
            DO i = 1, sys_size
                q_cons_ts(2)%vf(i)%sf(0:m,0:n,0:p) = &
                               q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) &
                             + dt*rhs_vf(i)%sf
            END DO
            
            IF (grid_geometry == 3) CALL s_apply_fourier_filter(q_cons_ts(2)%vf)

            IF (model_eqns == 3) CALL s_pressure_relaxation_procedure(q_cons_ts(2)%vf)

            ! ==================================================================
            

            ! Stage 2 of 3 =====================================================
            DO i = 1, cont_idx%end
                q_prim_vf(i)%sf => q_cons_ts(2)%vf(i)%sf
            END DO
            
            DO i = adv_idx%beg, adv_idx%end
                q_prim_vf(i)%sf => q_cons_ts(2)%vf(i)%sf
            END DO
            
            CALL s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, t_step)
            
            DO i = 1, sys_size
                q_cons_ts(2)%vf(i)%sf(0:m,0:n,0:p) = &
                           ( 3d0*q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) &
                           +     q_cons_ts(2)%vf(i)%sf(0:m,0:n,0:p) &
                           +     dt*rhs_vf(i)%sf ) / 4d0
            END DO
            
            IF (grid_geometry == 3) CALL s_apply_fourier_filter(q_cons_ts(2)%vf)

            IF (model_eqns == 3) CALL s_pressure_relaxation_procedure(q_cons_ts(2)%vf) 

            ! ==================================================================
            

            ! Stage 3 of 3 =====================================================
            CALL s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, t_step)
            
            DO i = 1, sys_size
                q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) = &
                           (     q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) &
                           + 2d0*q_cons_ts(2)%vf(i)%sf(0:m,0:n,0:p) &
                           + 2d0*dt*rhs_vf(i)%sf ) / 3d0
            END DO
            
            IF (grid_geometry == 3) CALL s_apply_fourier_filter(q_cons_ts(1)%vf)

            IF (model_eqns == 3) CALL s_pressure_relaxation_procedure(q_cons_ts(1)%vf)
      
            !IF ((model_eqns == 2) .and. bubbles .and. (num_fluids > 2)) &
            !    CALL s_phase_transfer(q_cons_ts(1)%vf)
            
            DO i = 1, cont_idx%end
                q_prim_vf(i)%sf => NULL()
            END DO
            
            DO i = adv_idx%beg, adv_idx%end
                q_prim_vf(i)%sf => NULL()
            END DO
            ! ==================================================================
            
            
        END SUBROUTINE s_3rd_order_tvd_rk ! ------------------------------------
        
        
        
        
        
        SUBROUTINE s_4th_order_rk(t_step) ! ------------------------------------
        ! Description: 4th order RK time-stepping algorithm
            
            
            ! Current time-step
            INTEGER, INTENT(IN) :: t_step
            
            ! Generic loop iterator
            INTEGER :: i
            
            
            ! Stage 1 of 4 =====================================================
            DO i = 1, cont_idx%end
                q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
            END DO
            
            DO i = adv_idx%beg, adv_idx%end
                q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
            END DO
            
            CALL s_compute_rhs(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, t_step)
            
            IF(run_time_info) THEN
                CALL s_write_run_time_information(q_prim_vf, t_step)
            END IF
            
            IF (ANY(com_wrt) .OR. ANY(cb_wrt) .OR. probe_wrt) THEN
                CALL s_time_step_cycling(t_step)
            END IF
            
            IF(t_step == t_step_stop) RETURN
            
            DO i = 1, sys_size
                q_cons_ts(2)%vf(i)%sf(0:m,0:n,0:p) = &
                               q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) &
                             + dt*rhs_vf(i)%sf / 2d0
                q_cons_ts(3)%vf(i)%sf(0:m,0:n,0:p) = &
                               q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) &
                             + dt*rhs_vf(i)%sf / 6d0
            END DO
            
            IF (grid_geometry == 3) THEN
                CALL s_apply_fourier_filter(q_cons_ts(2)%vf)
                CALL s_apply_fourier_filter(q_cons_ts(3)%vf)
            END IF

            IF (model_eqns == 3) THEN
                CALL s_pressure_relaxation_procedure(q_cons_ts(2)%vf)
                CALL s_pressure_relaxation_procedure(q_cons_ts(3)%vf)
            END IF
            ! ==================================================================
            
            
            ! Stage 2 of 4 =====================================================
            DO i = 1, cont_idx%end
                q_prim_vf(i)%sf => q_cons_ts(2)%vf(i)%sf
            END DO
            
            DO i = adv_idx%beg, adv_idx%end
                q_prim_vf(i)%sf => q_cons_ts(2)%vf(i)%sf
            END DO
            
            CALL s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, t_step)
            
            DO i = 1, sys_size
                q_cons_ts(2)%vf(i)%sf(0:m,0:n,0:p) = &
                               q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) &
                             + dt*rhs_vf(i)%sf / 2d0
                q_cons_ts(3)%vf(i)%sf(0:m,0:n,0:p) = &
                               q_cons_ts(3)%vf(i)%sf(0:m,0:n,0:p) &
                             + dt*rhs_vf(i)%sf / 3d0
            END DO
            
            IF (grid_geometry == 3) THEN
                CALL s_apply_fourier_filter(q_cons_ts(2)%vf)
                CALL s_apply_fourier_filter(q_cons_ts(3)%vf)
            END IF

            IF (model_eqns == 3) THEN
                CALL s_pressure_relaxation_procedure(q_cons_ts(2)%vf)
                CALL s_pressure_relaxation_procedure(q_cons_ts(3)%vf)
            END IF
            ! ==================================================================
            
            
            ! Stage 3 of 4 =====================================================
            CALL s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, t_step)
            
            DO i = 1, sys_size
                q_cons_ts(2)%vf(i)%sf(0:m,0:n,0:p) = &
                               q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) &
                             + dt*rhs_vf(i)%sf
                q_cons_ts(3)%vf(i)%sf(0:m,0:n,0:p) = &
                               q_cons_ts(3)%vf(i)%sf(0:m,0:n,0:p) &
                             + dt*rhs_vf(i)%sf / 3d0
            END DO
            
            IF (grid_geometry == 3) THEN
                CALL s_apply_fourier_filter(q_cons_ts(2)%vf)
                CALL s_apply_fourier_filter(q_cons_ts(3)%vf)
            END IF

            IF (model_eqns == 3) THEN
                CALL s_pressure_relaxation_procedure(q_cons_ts(2)%vf)
                CALL s_pressure_relaxation_procedure(q_cons_ts(3)%vf)
            END IF
            ! ==================================================================
            
            
            ! Stage 4 of 4 =====================================================
            CALL s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, t_step)
            
            DO i = 1, sys_size
                q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) = &
                               q_cons_ts(3)%vf(i)%sf(0:m,0:n,0:p) &
                             + dt*rhs_vf(i)%sf / 6d0
            END DO
            
            IF (grid_geometry == 3) CALL s_apply_fourier_filter(q_cons_ts(1)%vf)
            
            IF (model_eqns == 3) CALL s_pressure_relaxation_procedure(q_cons_ts(1)%vf)

            DO i = 1, cont_idx%end
                q_prim_vf(i)%sf => NULL()
            END DO
            
            DO i = adv_idx%beg, adv_idx%end
                q_prim_vf(i)%sf => NULL()
            END DO
            ! ==================================================================
            
            
        END SUBROUTINE s_4th_order_rk ! ----------------------------------------
        
        
        
        
        
        SUBROUTINE s_5th_order_rk(t_step) ! ------------------------------------
        ! Description: 5th order RK time-stepping algorithm
            
            
            ! Current time-step
            INTEGER, INTENT(IN) :: t_step
            
            ! Generic loop iterator
            INTEGER :: i
            
            
            ! Stage 1 of 6 =====================================================
            DO i = 1, cont_idx%end
                q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
            END DO
            
            DO i = adv_idx%beg, adv_idx%end
                q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
            END DO
            
            CALL s_compute_rhs(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, t_step)
            
            IF(run_time_info) THEN
                CALL s_write_run_time_information(q_prim_vf, t_step)
            END IF
            
            IF (ANY(com_wrt) .OR. ANY(cb_wrt) .OR. probe_wrt) THEN
                CALL s_time_step_cycling(t_step)
            END IF
            
            IF(t_step == t_step_stop) RETURN
            
            DO i = 1, sys_size
                q_cons_ts(2)%vf(i)%sf(0:m,0:n,0:p) = &
                               q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) &
                             + 2d-1*dt*rhs_vf(i)%sf
                q_cons_ts(3)%vf(i)%sf(0:m,0:n,0:p) = &
                               q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) &
                             + 75d-3*dt*rhs_vf(i)%sf
                q_cons_ts(4)%vf(i)%sf(0:m,0:n,0:p) = &
                               q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) &
                             + 3d-1*dt*rhs_vf(i)%sf
                q_cons_ts(5)%vf(i)%sf(0:m,0:n,0:p) = &
                               q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) &
                             - (11d0/54d0)*dt*rhs_vf(i)%sf
                q_cons_ts(6)%vf(i)%sf(0:m,0:n,0:p) = &
                               q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) &
                             + (1631d0/55296d0)*dt*rhs_vf(i)%sf
                q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) = &
                               q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) &
                             + (37d0/378d0)*dt*rhs_vf(i)%sf
            END DO
            
            IF (grid_geometry == 3) THEN
                CALL s_apply_fourier_filter(q_cons_ts(2)%vf)
                CALL s_apply_fourier_filter(q_cons_ts(3)%vf)
                CALL s_apply_fourier_filter(q_cons_ts(4)%vf)
                CALL s_apply_fourier_filter(q_cons_ts(5)%vf)
                CALL s_apply_fourier_filter(q_cons_ts(6)%vf)
                CALL s_apply_fourier_filter(q_cons_ts(1)%vf)
            END IF

            IF (model_eqns == 3) THEN
                CALL s_pressure_relaxation_procedure(q_cons_ts(2)%vf)
                CALL s_pressure_relaxation_procedure(q_cons_ts(3)%vf)
                CALL s_pressure_relaxation_procedure(q_cons_ts(4)%vf)
                CALL s_pressure_relaxation_procedure(q_cons_ts(5)%vf)
                CALL s_pressure_relaxation_procedure(q_cons_ts(6)%vf)
                CALL s_pressure_relaxation_procedure(q_cons_ts(1)%vf)
            END IF
            ! ==================================================================
            
            
            ! Stage 2 of 6 =====================================================
            DO i = 1, cont_idx%end
                q_prim_vf(i)%sf => q_cons_ts(2)%vf(i)%sf
            END DO
            
            DO i = adv_idx%beg, adv_idx%end
                q_prim_vf(i)%sf => q_cons_ts(2)%vf(i)%sf
            END DO
            
            CALL s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, t_step)
            
            DO i = 1, sys_size
                q_cons_ts(3)%vf(i)%sf(0:m,0:n,0:p) = &
                               q_cons_ts(3)%vf(i)%sf(0:m,0:n,0:p) &
                             + 225d-3*dt*rhs_vf(i)%sf
                q_cons_ts(4)%vf(i)%sf(0:m,0:n,0:p) = &
                               q_cons_ts(4)%vf(i)%sf(0:m,0:n,0:p) &
                             - 9d-1*dt*rhs_vf(i)%sf
                q_cons_ts(5)%vf(i)%sf(0:m,0:n,0:p) = &
                               q_cons_ts(5)%vf(i)%sf(0:m,0:n,0:p) &
                             + 25d-1*dt*rhs_vf(i)%sf
                q_cons_ts(6)%vf(i)%sf(0:m,0:n,0:p) = &
                               q_cons_ts(6)%vf(i)%sf(0:m,0:n,0:p) &
                             + (175d0/512d0)*dt*rhs_vf(i)%sf
            END DO

            IF (grid_geometry == 3) THEN
                CALL s_apply_fourier_filter(q_cons_ts(3)%vf)
                CALL s_apply_fourier_filter(q_cons_ts(4)%vf)
                CALL s_apply_fourier_filter(q_cons_ts(5)%vf)
                CALL s_apply_fourier_filter(q_cons_ts(6)%vf)
            END IF

            IF (model_eqns == 3) THEN
                CALL s_pressure_relaxation_procedure(q_cons_ts(3)%vf)
                CALL s_pressure_relaxation_procedure(q_cons_ts(4)%vf)
                CALL s_pressure_relaxation_procedure(q_cons_ts(5)%vf)
                CALL s_pressure_relaxation_procedure(q_cons_ts(6)%vf)
            END IF
            ! ==================================================================
            
            
            ! Stage 3 of 6 =====================================================
            DO i = 1, cont_idx%end
                q_prim_vf(i)%sf => q_cons_ts(3)%vf(i)%sf
            END DO
            
            DO i = adv_idx%beg, adv_idx%end
                q_prim_vf(i)%sf => q_cons_ts(3)%vf(i)%sf
            END DO
            
            CALL s_compute_rhs(q_cons_ts(3)%vf, q_prim_vf, rhs_vf, t_step)
            
            DO i = 1, sys_size
                q_cons_ts(4)%vf(i)%sf(0:m,0:n,0:p) = &
                               q_cons_ts(4)%vf(i)%sf(0:m,0:n,0:p) &
                             + 12d-1*dt*rhs_vf(i)%sf
                q_cons_ts(5)%vf(i)%sf(0:m,0:n,0:p) = &
                               q_cons_ts(5)%vf(i)%sf(0:m,0:n,0:p) &
                             - (7d1/27d0)*dt*rhs_vf(i)%sf
                q_cons_ts(6)%vf(i)%sf(0:m,0:n,0:p) = &
                               q_cons_ts(6)%vf(i)%sf(0:m,0:n,0:p) &
                             + (575d0/13824d0)*dt*rhs_vf(i)%sf
                q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) = &
                               q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) &
                             + (25d1/621d0)*dt*rhs_vf(i)%sf
            END DO

            IF (grid_geometry == 3) THEN
                CALL s_apply_fourier_filter(q_cons_ts(4)%vf)
                CALL s_apply_fourier_filter(q_cons_ts(5)%vf)
                CALL s_apply_fourier_filter(q_cons_ts(6)%vf)
                CALL s_apply_fourier_filter(q_cons_ts(1)%vf)
            END IF

            IF (model_eqns == 3) THEN
                CALL s_pressure_relaxation_procedure(q_cons_ts(4)%vf)
                CALL s_pressure_relaxation_procedure(q_cons_ts(5)%vf)
                CALL s_pressure_relaxation_procedure(q_cons_ts(6)%vf)
                CALL s_pressure_relaxation_procedure(q_cons_ts(1)%vf)
            END IF
            ! ==================================================================
            
            
            ! Stage 4 of 6 =====================================================
            DO i = 1, cont_idx%end
                q_prim_vf(i)%sf => q_cons_ts(4)%vf(i)%sf
            END DO
            
            DO i = adv_idx%beg, adv_idx%end
                q_prim_vf(i)%sf => q_cons_ts(4)%vf(i)%sf
            END DO
            
            CALL s_compute_rhs(q_cons_ts(4)%vf, q_prim_vf, rhs_vf, t_step)
            
            DO i = 1, sys_size
                q_cons_ts(5)%vf(i)%sf(0:m,0:n,0:p) = &
                               q_cons_ts(5)%vf(i)%sf(0:m,0:n,0:p) &
                             + (35d0/27d0)*dt*rhs_vf(i)%sf
                q_cons_ts(6)%vf(i)%sf(0:m,0:n,0:p) = &
                               q_cons_ts(6)%vf(i)%sf(0:m,0:n,0:p) &
                             + (44275d0/110592d0)*dt*rhs_vf(i)%sf
                q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) = &
                               q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) &
                             + (125d0/594d0)*dt*rhs_vf(i)%sf
            END DO

            IF (grid_geometry == 3) THEN
                CALL s_apply_fourier_filter(q_cons_ts(5)%vf)
                CALL s_apply_fourier_filter(q_cons_ts(6)%vf)
                CALL s_apply_fourier_filter(q_cons_ts(1)%vf)
            END IF

            IF (model_eqns == 3) THEN
                CALL s_pressure_relaxation_procedure(q_cons_ts(5)%vf)
                CALL s_pressure_relaxation_procedure(q_cons_ts(6)%vf)
                CALL s_pressure_relaxation_procedure(q_cons_ts(1)%vf)
            END IF
            ! ==================================================================
            
            
            ! Stage 5 of 6 =====================================================
            DO i = 1, cont_idx%end
                q_prim_vf(i)%sf => q_cons_ts(5)%vf(i)%sf
            END DO
            
            DO i = adv_idx%beg, adv_idx%end
                q_prim_vf(i)%sf => q_cons_ts(5)%vf(i)%sf
            END DO
            
            CALL s_compute_rhs(q_cons_ts(5)%vf, q_prim_vf, rhs_vf, t_step)
            
            DO i = 1, sys_size
                q_cons_ts(6)%vf(i)%sf(0:m,0:n,0:p) = &
                               q_cons_ts(6)%vf(i)%sf(0:m,0:n,0:p) &
                             + (253d0/4096d0)*dt*rhs_vf(i)%sf
            END DO

            IF (grid_geometry == 3) CALL s_apply_fourier_filter(q_cons_ts(6)%vf)

            IF (model_eqns == 3) CALL s_pressure_relaxation_procedure(q_cons_ts(6)%vf)
            ! ==================================================================
            
            
            ! Stage 6 of 6 =====================================================
            DO i = 1, cont_idx%end
                q_prim_vf(i)%sf => q_cons_ts(6)%vf(i)%sf
            END DO
            
            DO i = adv_idx%beg, adv_idx%end
                q_prim_vf(i)%sf => q_cons_ts(6)%vf(i)%sf
            END DO
            
            CALL s_compute_rhs(q_cons_ts(6)%vf, q_prim_vf, rhs_vf, t_step)
            
            DO i = 1, sys_size
                q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) = &
                               q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) &
                             + (512d0/1771d0)*dt*rhs_vf(i)%sf
            END DO
            
            IF (grid_geometry == 3) CALL s_apply_fourier_filter(q_cons_ts(1)%vf)

            IF (model_eqns == 3) CALL s_pressure_relaxation_procedure(q_cons_ts(1)%vf)

            DO i = 1, cont_idx%end
                q_prim_vf(i)%sf => NULL()
            END DO
            
            DO i = adv_idx%beg, adv_idx%end
                q_prim_vf(i)%sf => NULL()
            END DO
            ! ==================================================================
            
            
        END SUBROUTINE s_5th_order_rk ! ----------------------------------------
        
        
        
        
        
        SUBROUTINE s_time_step_cycling(t_step) ! ----------------------------
        ! Description: This subroutine saves the temporary q_prim_vf vector 
        !          into the q_prim_ts vector that is then used in p_main

            ! Current time-step
            INTEGER, INTENT(IN) :: t_step

            ! Generic loop iterator
            INTEGER :: i

            IF (t_step == t_step_start) THEN
                DO i = 1, sys_size
                    q_prim_ts(3)%vf(i)%sf(:,:,:) = q_prim_vf(i)%sf(:,:,:)
                END DO
            ELSEIF (t_step == t_step_start + 1) THEN
                DO i = 1, sys_size
                    q_prim_ts(2)%vf(i)%sf(:,:,:) = q_prim_vf(i)%sf(:,:,:)
                END DO
            ELSEIF (t_step == t_step_start + 2) THEN
                DO i = 1, sys_size
                    q_prim_ts(1)%vf(i)%sf(:,:,:) = q_prim_vf(i)%sf(:,:,:)
                END DO
            ELSEIF (t_step == t_step_start + 3) THEN
                DO i = 1, sys_size
                    q_prim_ts(0)%vf(i)%sf(:,:,:) = q_prim_vf(i)%sf(:,:,:)
                END DO
            ELSE ! All other timesteps
                DO i = 1, sys_size
                    q_prim_ts(3)%vf(i)%sf(:,:,:) = q_prim_ts(2)%vf(i)%sf(:,:,:)
                    q_prim_ts(2)%vf(i)%sf(:,:,:) = q_prim_ts(1)%vf(i)%sf(:,:,:)
                    q_prim_ts(1)%vf(i)%sf(:,:,:) = q_prim_ts(0)%vf(i)%sf(:,:,:)
                    q_prim_ts(0)%vf(i)%sf(:,:,:) = q_prim_vf(i)%sf(:,:,:)
                END DO
            END IF


        END SUBROUTINE s_time_step_cycling ! -----------------------------------





        SUBROUTINE s_finalize_time_steppers_module() ! -------------------------
        ! Description: Module deallocation and/or disassociation procedures
            
            
            ! Generic loop iterators
            INTEGER :: i,j
            
            
            ! Deallocating the cell-average conservative variables
            DO i = 1, num_ts
                
                DO j = 1, sys_size
                    DEALLOCATE(q_cons_ts(i)%vf(j)%sf)
                END DO
                
                DEALLOCATE(q_cons_ts(i)%vf)
                
            END DO
            
            DEALLOCATE(q_cons_ts)


            ! Deallocating the cell-average primitive ts variables
            IF (ANY(com_wrt) .OR. ANY(cb_wrt) .OR. probe_wrt) THEN
                DO i = 0, 3
                    DO j = 1, sys_size
                        DEALLOCATE(q_prim_ts(i)%vf(j)%sf)
                    END DO
                    DEALLOCATE(q_prim_ts(i)%vf)
                END DO
                DEALLOCATE(q_prim_ts)
            END IF
            
            
            ! Deallocating the cell-average primitive variables
            DO i = mom_idx%beg, E_idx
                DEALLOCATE(q_prim_vf(i)%sf)
            END DO
            IF (model_eqns == 3) THEN
                DO i = internalEnergies_idx%beg, internalEnergies_idx%end
                    DEALLOCATE(q_prim_vf(i)%sf)
                END DO
            END IF
            
            DEALLOCATE(q_prim_vf)
            
            
            ! Deallocating the cell-average RHS variables
            DO i = 1, sys_size
                DEALLOCATE(rhs_vf(i)%sf)
            END DO
            
            DEALLOCATE(rhs_vf)
            
            
            ! Writing the footer of and closing the run-time information file
            IF(proc_rank == 0 .AND. run_time_info) THEN
                CALL s_close_run_time_information_file()
            END IF
            
            
        END SUBROUTINE s_finalize_time_steppers_module ! -----------------------
        
        
        
        
        
END MODULE m_time_steppers
