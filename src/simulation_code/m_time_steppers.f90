!>
!! @file m_time_steppers.f90
!! @brief Contains module m_time_steppers
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief The following module features a variety of time-stepping schemes.
!!              Currently, it includes the following Runge-Kutta (RK) algorithms:
!!                   1) 1st Order TVD RK
!!                   2) 2nd Order TVD RK
!!                   3) 3rd Order TVD RK
!!                   4) 4th Order RK
!!                   5) 5th Order RK
!!              where TVD designates a total-variation-diminishing time-stepper.
MODULE m_time_steppers
    
    
    ! Dependencies =============================================================
    USE m_derived_types        !< Definitions of the derived types
    
    USE m_global_parameters    !< Definitions of the global parameters
    
    USE m_fftw                 !< Module for FFTW functions
    
    USE m_rhs                  !< Right-hand-side (RHS) evaluation procedures
    
    USE m_data_output          !< Run-time info & solution data output procedures

    USE m_bubbles              !< Bubble dynamics routines

    USE m_mpi_proxy            !< Message passing interface (MPI) module proxy
    ! ==========================================================================
    
    
    IMPLICIT NONE
    
    

    TYPE(vector_field), ALLOCATABLE, DIMENSION(:) :: q_cons_ts !<
    !! Cell-average conservative variables at each time-stage (TS)
    
    TYPE(scalar_field), PRIVATE, ALLOCATABLE, DIMENSION(:) :: q_prim_vf !<
    !! Cell-average primitive variables at the current time-stage
    

    TYPE(scalar_field), ALLOCATABLE, DIMENSION(:) :: rhs_vf !<
    !! Cell-average RHS variables at the current time-stage
    

    TYPE(vector_field), ALLOCATABLE, DIMENSION(:) :: q_prim_ts !<
    !! Cell-average primitive variables at consecutive TIMESTEPS


    INTEGER, PRIVATE :: num_ts !<
    !! Number of time stages in the time-stepping scheme
    
    
    CONTAINS
        
        
        
        !> The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
        SUBROUTINE s_initialize_time_steppers_module() ! -----------------------
           
            

            TYPE(bounds_info) :: ix,iy,iz !<
            !! Indical bounds in the x-, y- and z-directions
            

            INTEGER :: i,j !< Generic loop iterators
            
            
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

            IF (bubbles) THEN
                DO i = bub_idx%beg,sys_size
                    ALLOCATE(q_prim_vf(i)%sf( ix%beg:ix%end, &
                                          iy%beg:iy%end, &
                                          iz%beg:iz%end ))
                END DO
            END IF

            IF (hypoelasticity) THEN
                DO i = stress_idx%beg, stress_idx%end
                    ALLOCATE(q_prim_vf(i)%sf( ix%beg:ix%end, &
                                              iy%beg:iy%end, &
                                              iz%beg:iz%end ))
                END DO
            END IF

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
        
        
        
        
        
        !> 1st order TVD RK time-stepping algorithm
        !! @param t_step Current time step
        SUBROUTINE s_1st_order_tvd_rk(t_step) ! --------------------------------

            INTEGER, INTENT(IN) :: t_step
            
            INTEGER :: i !< Generic loop iterator
           
            ! Stage 1 of 1 =====================================================
            DO i = 1, cont_idx%end
                q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
            END DO

            IF (adv_alphan) THEN
                DO i = adv_idx%beg, adv_idx%end
                        q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
                END DO
            ELSE 
                DO i = adv_idx%beg, sys_size
                    q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
                END DO
            END IF

            CALL s_compute_rhs(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, t_step)
            IF (DEBUG) PRINT*, 'got rhs'

            IF(run_time_info) THEN
                CALL s_write_run_time_information(q_prim_vf, t_step)
            END IF
            IF (DEBUG) print*, 'wrote runtime info'

            IF (ANY(com_wrt) .OR. ANY(cb_wrt) .OR. probe_wrt) THEN
                CALL s_time_step_cycling(t_step)
            END IF
           
            
            IF(t_step == t_step_stop) RETURN
           
            DO i = 1, sys_size
                q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) = &
                               q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) &
                             + dt*rhs_vf(i)%sf
            END DO


            IF (grid_geometry == 3) CALL s_apply_fourier_filter(q_cons_ts(1)%vf)

            IF (model_eqns == 3) CALL s_pressure_relaxation_procedure(q_cons_ts(1)%vf)
            
            DO i = 1, cont_idx%end
                q_prim_vf(i)%sf => NULL()
            END DO
           
            IF (adv_alphan) THEN
                DO i = adv_idx%beg, adv_idx%end
                    q_prim_vf(i)%sf => NULL()
                END DO
            ELSE
                DO i = adv_idx%beg, sys_size ! adv_idx%end
                    q_prim_vf(i)%sf => NULL()
                END DO 
            END IF
            ! ==================================================================
           
        END SUBROUTINE s_1st_order_tvd_rk ! ------------------------------------
        
        
        
        
        !> 2nd order TVD RK time-stepping algorithm
        !! @param t_step Current time-step
        SUBROUTINE s_2nd_order_tvd_rk(t_step) ! --------------------------------

            INTEGER, INTENT(IN) :: t_step
            
            INTEGER :: i !< Generic loop iterator
            
            
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
        
        
        
        
        
        !> 3rd order TVD RK time-stepping algorithm
        !! @param t_step Current time-step
        SUBROUTINE s_3rd_order_tvd_rk(t_step) ! --------------------------------

            INTEGER, INTENT(IN) :: t_step
            
            INTEGER :: i,j !< Generic loop iterator
            
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
      
            DO i = 1, cont_idx%end
                q_prim_vf(i)%sf => NULL()
            END DO
            
            DO i = adv_idx%beg, adv_idx%end
                q_prim_vf(i)%sf => NULL()
            END DO
            ! ==================================================================
            
            
        END SUBROUTINE s_3rd_order_tvd_rk ! ------------------------------------
        


        !> Adaptive SSP RK23 time-stepping algorithm
        !! @param t_step Current time-step
        SUBROUTINE s_23_order_tvd_rk(t_step) ! --------------------------------

            INTEGER, INTENT(IN) :: t_step
            REAL(KIND(0d0)) :: relerr, absval, tmp
            REAL(KIND(0d0)) :: dtmin,dtmax
            
            INTEGER :: i,j !< Generic loop iterator
            
            ! Stage 1 of 3 =====================================================
            DO i = 1, cont_idx%end
                q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
            END DO
            DO i = adv_idx%beg, adv_idx%end
                q_prim_vf(i)%sf => q_cons_ts(1)%vf(i)%sf
            END DO
            
            CALL s_compute_rhs(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, t_step)
            
            IF(run_time_info) CALL s_write_run_time_information(q_prim_vf, t_step)
            IF (ANY(com_wrt) .OR. ANY(cb_wrt) .OR. probe_wrt) &
                CALL s_time_step_cycling(t_step)
            
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

            ! RK2 estimate
            DO i = 1, sys_size
                q_cons_ts(3)%vf(i)%sf(0:m,0:n,0:p) = (              &
                                q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p)  &
                           +    q_cons_ts(2)%vf(i)%sf(0:m,0:n,0:p)  &
                           +    dt*rhs_vf(i)%sf ) / 2d0
            END DO

            IF (grid_geometry == 3) CALL s_apply_fourier_filter(q_cons_ts(3)%vf)
            IF (model_eqns == 3) CALL s_pressure_relaxation_procedure(q_cons_ts(3)%vf) 
            
            ! Stage 2 of RK3
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

            ! ==================================================================


            ! Approximate error =================================================
            ! err = (q_cons_ts(1)%vf(i)%sf - q_cons_ts(3)%vf(i)%sf) / &
            !     q_cons_ts(1)%vf(i)

            ! PRINT*, '          '
            ! DO i = 1,sys_size
            !     PRINT*, 'MAXVAL', i,  MAXVAL( ABS( q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p)  ), 1 )
            ! END DO

            ! DO i = 1,sys_size
            !     PRINT*, 'ABSERR', i,  MAXVAL( ABS( q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) - &
            !         q_cons_ts(3)%vf(i)%sf(0:m,0:n,0:p)  ), 1 )
            ! END DO

            relerr = 0d0
            DO i = 1,sys_size
                absval = MAXVAL( ABS( q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p)  ) )
                IF ( absval >= 1D-10 ) THEN
                    relerr =  MAX( relerr, MAXVAL( ABS(             &
                       (q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) -        &
                        q_cons_ts(3)%vf(i)%sf(0:m,0:n,0:p)) /       &
                        q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) ) )   &
                        )
                END IF
            END DO

            IF (num_procs > 1) THEN
                tmp = relerr
                CALL s_mpi_allreduce_max(tmp,relerr)
            END IF

            dtmin = 0.002d0 / 2d0
            dtmax = 0.002d0 * 2d0

            dt = dt*Min( Max( Sqrt(t_tol/(2d0*relerr)) , 0.3d0 ), 2d0 )
            dt = Max( Min( dtmax, dt ), dtmin )
            ! dt = 0.0015d0

            IF (proc_rank==0) PRINT*, 'RELERR:', relerr
            IF (proc_rank==0) PRINT*, 'dt/dt0:', dt/dt0
            ! IF (proc_rank==0) PRINT*, '---t/T:', mytime/finaltime

            ! ==================================================================
      
            DO i = 1, cont_idx%end
                q_prim_vf(i)%sf => NULL()
            END DO
            DO i = adv_idx%beg, adv_idx%end
                q_prim_vf(i)%sf => NULL()
            END DO
            ! ==================================================================
            
            
        END SUBROUTINE s_23_order_tvd_rk ! ------------------------------------
        
        
        !> 4th order RK time-stepping algorithm
        !! @param t_step Current time-step
        SUBROUTINE s_4th_order_rk(t_step) ! ------------------------------------

            INTEGER, INTENT(IN) :: t_step
            
            INTEGER :: i !< Generic loop iterator
            
            
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
        
        
        
        
        !> 5th order RK time-stepping algorithm
        !! @param t_step Current time-step
        SUBROUTINE s_5th_order_rk(t_step) ! ------------------------------------

            INTEGER, INTENT(IN) :: t_step
            
            INTEGER :: i !< Generic loop iterator
            
            
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
        
        
        
        !> This subroutine saves the temporary q_prim_vf vector 
        !!      into the q_prim_ts vector that is then used in p_main        
        !! @param t_step current time-step
        SUBROUTINE s_time_step_cycling(t_step) ! ----------------------------

            INTEGER, INTENT(IN) :: t_step

            INTEGER :: i !< Generic loop iterator

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




        !> Module deallocation and/or disassociation procedures
        SUBROUTINE s_finalize_time_steppers_module() ! -------------------------

            INTEGER :: i,j !< Generic loop iterators
            
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
