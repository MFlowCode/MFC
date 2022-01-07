!>
!! @file m_data_output.f90
!! @brief Contains module m_data_output

!> @brief The primary purpose of this module is to output the grid and the
!!              conservative variables data at the chosen time-step interval. In
!!              addition, this module is also in charge of outputting a run-time
!!              information file which summarizes the time-dependent behavior !of
!!              the stability criteria. The latter include the inviscid Courant–
!!              Friedrichs–Lewy (ICFL), viscous CFL (VCFL), capillary CFL (CCFL)
!!              and cell Reynolds (Rc) numbers.
MODULE m_data_output
    
    
    !  Dependencies ============================================================
    USE m_derived_types        !< Definitions of the derived types
    
    USE m_global_parameters    !< Definitions of the global parameters
    
    USE m_mpi_proxy            !< Message passing interface (MPI) module proxy
    
    USE m_variables_conversion !< State variables type conversion procedures

    USE m_compile_specific
    ! ==========================================================================
    
    
    IMPLICIT NONE
    
    PRIVATE; PUBLIC :: s_initialize_data_output_module  , &
                       s_open_run_time_information_file , &
                       s_open_com_files                 , &
                       s_open_cb_files                  , &
                       s_open_probe_files               , &
                       s_write_run_time_information     , &
                       s_write_data_files               , &
                       s_write_serial_data_files        , &
                       s_write_parallel_data_files      , &
                       s_write_com_files                , &
                       s_write_cb_files                 , &
                       s_write_probe_files              , &
                       s_close_run_time_information_file, &
                       s_close_com_files                , &
                       s_close_cb_files                 , &
                       s_close_probe_files              , &
                       s_finalize_data_output_module
    
    ABSTRACT INTERFACE ! ===================================================

        !> Write data files
        !! @param q_cons_vf Conservative variables
        !! @param t_step Current time step
        SUBROUTINE s_write_abstract_data_files(q_cons_vf, t_step)

            IMPORT :: scalar_field, sys_size


            TYPE(scalar_field), &
            DIMENSION(sys_size), &
            INTENT(IN) :: q_cons_vf

            INTEGER, INTENT(IN) :: t_step

        END SUBROUTINE s_write_abstract_data_files ! -------------------
    END INTERFACE ! ========================================================
    
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:,:) :: icfl_sf  !< ICFL stability criterion
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:,:) :: vcfl_sf  !< VCFL stability criterion
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:,:) :: ccfl_sf  !< CCFL stability criterion
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:,:) ::   Rc_sf  !< Rc stability criterion
    
    !> @name ICFL, VCFL, CCFL and Rc stability criteria extrema over all the time-steps
    !> @{
    REAL(KIND(0d0)) :: icfl_max !< ICFL criterion maximum
    REAL(KIND(0d0)) :: vcfl_max !< VCFL criterion maximum
    REAL(KIND(0d0)) :: ccfl_max !< CCFL criterion maximum
    REAL(KIND(0d0)) ::   Rc_min !< Rc criterion maximum
    !> @}

    ! @name Generic storage for flow variable(s) that are to be written to CoM data file
    !> @{
    REAL(KIND(0d0)), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:)  :: accel_mag 
    REAL(KIND(0d0)), PUBLIC, ALLOCATABLE, DIMENSION(:,:)    :: q_com
    REAL(KIND(0d0)), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:)  :: moments
    REAL(KIND(0d0)), PUBLIC, ALLOCATABLE, DIMENSION(:,:)    :: cb_mass
    REAL(KIND(0d0)), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:)  :: bounds
    REAL(KIND(0d0)), PUBLIC, ALLOCATABLE, DIMENSION(:,:)    :: cntrline
    REAL(KIND(0d0)), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:)  :: x_accel, y_accel, z_accel
    TYPE(scalar_field), ALLOCATABLE, DIMENSION(:)           :: grad_x_vf,grad_y_vf,grad_z_vf,norm_vf,kappa_vf
    REAL(KIND(0d0)), TARGET, ALLOCATABLE, DIMENSION(:,:,:)  :: energy !< Energy: 
    !! Used to write out correct E and p when We_size > 0
    !> @}

    PROCEDURE(s_write_abstract_data_files), POINTER :: s_write_data_files => NULL()

    
    CONTAINS
        

        !>  The purpose of this subroutine is to open a new or pre-
        !!          existing run-time information file and append to it the
        !!      basic header information relevant to current simulation.
        !!      In general, this requires generating a table header for
        !!      those stability criteria which will be written at every
        !!      time-step.        
        SUBROUTINE s_open_run_time_information_file() ! ------------------------


            CHARACTER(LEN = name_len) :: dir_name !< Name of the case directory

            CHARACTER(LEN = name_len) :: file_name = 'run_time.inf' !<
            !! Name of the run-time information file
            
            CHARACTER(LEN = path_len + name_len) :: file_path !<
            !! Relative path to a file in the case directory
            
            CHARACTER(LEN = 8) :: file_date !<
            !! Creation date of the run-time information file
            
            LOGICAL :: file_exist !<
            !! Logical used to check existence of run-time information file
            
            ! Opening the run-time information file
            file_path = TRIM(case_dir) // '/' // TRIM(file_name)
            
            INQUIRE(FILE = TRIM(file_path), EXIST = file_exist)
            
            OPEN(1, FILE     = TRIM(file_path), &
                    FORM     = 'formatted'    , &
                    POSITION = 'append'       , &
                    STATUS   = 'unknown'        )
            
            
            ! Generating file header for a new run-time information file
            IF(file_exist .NEQV. .TRUE.) THEN
                
                CALL SYSTEM('basename `pwd` > dir_name')
                
                file_path = TRIM(case_dir) // '/dir_name'
                
                OPEN(2, FILE   = TRIM(file_path), &
                        FORM   = 'formatted'    , &
                        STATUS = 'old'            )
                
                READ(2,'(A)') dir_name; CLOSE(2)
                
                CALL SYSTEM('rm -f dir_name')
                
                WRITE(1,'(A)') 'MFC - Case - ' // TRIM( dir_name) // &
                                                ': ' // TRIM(file_name)
                WRITE(1,'(A)')     'Description: Stability information at ' // &
                                   'each time-step of the simulation. This'
                WRITE(1,'(13X,A)') 'data is composed of the inviscid '      // &
                                   'Courant–Friedrichs–Lewy (ICFL)'
                WRITE(1,'(13X,A)') 'number, the viscous CFL (VCFL) number, '// &
                                   'the capillary CFL (CCFL)'
                WRITE(1,'(13X,A)') 'number and the cell Reynolds (Rc) '     // &
                                   'number. Please note that only'
                WRITE(1,'(13X,A)') 'those stability conditions pertinent '  // &
                                   'to the physics included in'
                WRITE(1,'(13X,A)') 'the current computation are displayed.'
                
                CALL DATE_AND_TIME(DATE = file_date)
                
                WRITE(1,'(A)') 'Date: ' // file_date(5:6) // '/' // &
                                           file_date(7:8) // '/' // &
                                           file_date(3:4)
                
            END IF
            
            WRITE(1,'(A)') ''; WRITE(1,'(A)') ''
            
            
            ! Generating table header for the stability criteria to be outputted
            IF(ANY(Re_size > 0) .AND. We_size > 0) THEN
                WRITE(1,'(A)') '== Time-steps ==== Time ===== ICFL Max == ' // &
                               'VCFL Max == CCFL Max ==== Rc Min ====='
            ELSEIF(ANY(Re_size > 0)) THEN
                WRITE(1,'(A)') '==== Time-steps ====== Time ======= ICFL '  // &
                               'Max ==== VCFL Max ====== Rc Min ======='
            ELSEIF(We_size > 0) THEN
                WRITE(1,'(A)') '======= Time-steps ========= Time '         // &
                               '========== ICFL Max ======= CCFL '          // &
                               'Max ========='
            ELSE
                WRITE(1,'(A)') '=========== Time-steps ============== Time '// &
                               '============== ICFL Max ============='
            END IF
            
            
        END SUBROUTINE s_open_run_time_information_file ! ----------------------
        
        
        
        !>  This opens a formatted data file where the root processor
        !!      can write out the CoM information        
        SUBROUTINE s_open_com_files() ! ----------------------------------------



            CHARACTER(LEN = path_len  + 3*name_len) :: file_path !<
            !! Relative path to the CoM file in the case directory 
        

            INTEGER :: i !< Generic loop iterator
        
            DO i = 1, num_fluids
                IF (com_wrt(i)) THEN
        
                    ! Generating the relative path to the CoM data file
                    WRITE(file_path,'(A,I0,A)') '/fluid',i,'_com.dat'
                    file_path = TRIM(case_dir) // TRIM(file_path)
        
                    ! Creating the formatted data file and setting up its
                    ! structure
                    OPEN(i+10, FILE = TRIM(file_path), &
                        FORM = 'formatted', &
                        POSITION = 'append', &
                        STATUS   = 'unknown')
        
                    IF (n == 0) THEN
                        IF (ANY(moment_order /= dflt_int)) THEN
                            WRITE(i+10,'(A)') '=== Non-Dimensional Time ' // &
                                     '=== Total Mass ' // &
                                     '=== x-loc ' // &
                                     '=== x-vel ' // &
                                     '=== x-accel ' // &
                                     '=== Higher Moments ==='
                        ELSE
                            WRITE(i+10,'(A)') '=== Non-Dimensional Time ' // &
                                     '=== Total Mass ' // &
                                     '=== x-loc ' // &
                                     '=== x-vel ' // &
                                     '=== x-accel ==='
                        END IF
                    ELSEIF (p == 0) THEN
                        IF (ANY(moment_order /= dflt_int)) THEN
                            WRITE(i+10,'(A)') '=== Non-Dimensional Time ' // &
                                     '=== Total Mass ' // &
                                     '=== x-loc ' // &
                                     '=== y-loc ' // &
                                     '=== x-vel ' // &
                                     '=== y-vel ' // &
                                     '=== x-accel ' // &
                                     '=== y-accel ' // &
                                     '=== Higher Moments ==='
                        ELSE
                            WRITE(i+10,'(A)') '=== Non-Dimensional Time ' // &
                                     '=== Total Mass ' // &
                                     '=== x-loc ' // &
                                     '=== y-loc ' // &
                                     '=== x-vel ' // &
                                     '=== y-vel ' // &
                                     '=== x-accel ' // &
                                     '=== y-accel ==='
                        END IF
                    ELSE
                        IF (ANY(moment_order /= dflt_int)) THEN
                            WRITE(i+10,'(A)') '=== Non-Dimensional Time ' // &
                                     '=== Total Mass ' // &
                                     '=== x-loc ' // &
                                     '=== y-loc ' // &
                                     '=== z-loc ' // &
                                     '=== x-vel ' // &
                                     '=== y-vel ' // &
                                     '=== z-vel ' // &
                                     '=== x-accel ' // &
                                     '=== y-accel ' // &
                                     '=== z-accel ' // &
                                     '=== Higher Moments ==='
                        ELSE
                            WRITE(i+10,'(A)') '=== Non-Dimensional Time ' // &
                                     '=== Total Mass ' // &
                                     '=== x-loc ' // &
                                     '=== y-loc ' // &
                                     '=== z-loc ' // &
                                     '=== x-vel ' // &
                                     '=== y-vel ' // &
                                     '=== z-vel ' // &
                                     '=== x-accel ' // &
                                     '=== y-accel ' // &
                                     '=== z-accel ==='
                        END IF
                    END IF
                END IF
            END DO

        END SUBROUTINE s_open_com_files ! --------------------------------------



        !>  This opens a formatted data file where the root processor
        !!      can write out the coherent body (cb) information
        SUBROUTINE s_open_cb_files() ! ----------------------------------------


            CHARACTER(LEN = path_len  + 3*name_len) :: file_path !<
            !! Relative path to the cb file in the case directory 
        

            INTEGER :: i  !< Generic loop iterator
        
            DO i = 1, num_fluids
                IF (cb_wrt(i)) THEN
        
                    ! Generating the relative path to the cb data file
                    WRITE(file_path,'(A,I0,A)') '/fluid',i,'_cb.dat'
                    file_path = TRIM(case_dir) // TRIM(file_path)
        
                    ! Creating the formatted data file and setting up its
                    ! structure
                    OPEN(i+20, FILE = TRIM(file_path), &
                        FORM = 'formatted', &
                        POSITION = 'append', &
                        STATUS   = 'unknown')

                    WRITE(i+20, '(A)') '=== Non-Dimensional Time ' // &
                                     '=== Coherent Body Mass ' // &
                                     '=== Coherent Body Area ' // &
                                     '=== Bounds ' // &
                                     '=== Centerline ==='
                END IF
            END DO

        END SUBROUTINE s_open_cb_files ! --------------------------------------




        !>  This opens a formatted data file where the root processor
        !!      can write out flow probe information
        SUBROUTINE s_open_probe_files() ! --------------------------------------

            CHARACTER(LEN = path_len + 3*name_len) :: file_path !<
            !! Relative path to the probe data file in the case directory


            INTEGER :: i !< Generic loop iterator

            DO i = 1, num_probes
                ! Generating the relative path to the data file
                WRITE(file_path,'(A,I0,A)') '/D/probe',i,'_prim.dat'
                file_path = TRIM(case_dir) // TRIM(file_path)

                ! Creating the formatted data file and setting up its
                ! structure
                OPEN(i+30, FILE = TRIM(file_path), &
                    FORM = 'formatted', &
                    STATUS   = 'unknown')
                    ! POSITION = 'append', &
                !WRITE(i+30,'(A,I0,A)') 'Probe ',i, ' located at:'
                !WRITE(i+30,'(A,F10.6)') 'x = ',probe(i)%x
                !WRITE(i+30,'(A,F10.6)') 'y = ',probe(i)%y
                !WRITE(i+30,'(A,F10.6)') 'z = ',probe(i)%z
                !WRITE(i+30, *)
                !WRITE(i+30,'(A)') '=== Non-Dimensional Time ' // &
                !                '=== Density ' // &
                !                '=== Velocity ' // &
                !                '=== Pressure ' // &
                !                '=== Gamma ' // &
                !                '=== Stiffness ' // &
                !                '=== Sound Speed ' // &
                !                '=== Acceleration ==='
            END DO


            IF (integral_wrt) THEN
            DO i = 1, num_integrals
                WRITE(file_path,'(A,I0,A)') '/D/integral',i,'_prim.dat'
                file_path = TRIM(case_dir) // TRIM(file_path)

                OPEN(i+70, FILE = TRIM(file_path), &
                    FORM = 'formatted', &
                    POSITION = 'append', &
                    STATUS   = 'unknown')
            END DO
            end IF


        END SUBROUTINE s_open_probe_files ! ------------------------------------



        !>  The goal of the procedure is to output to the run-time
        !!      information file the stability criteria extrema in the
        !!      entire computational domain and at the given time-step.
        !!      Moreover, the subroutine is also in charge of tracking
        !!      these stability criteria extrema over all time-steps.
        !!  @param q_prim_vf Cell-average primitive variables
        !!  @param t_step Current time step
        SUBROUTINE s_write_run_time_information(q_prim_vf, t_step) ! -----------
           
            
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q_prim_vf
            INTEGER, INTENT(IN) :: t_step
            
            REAL(KIND(0d0)), DIMENSION(cont_idx%end)          :: alpha_rho  !< Cell-avg. partial density
            REAL(KIND(0d0))                                   :: rho        !< Cell-avg. density
            REAL(KIND(0d0)), DIMENSION(num_dims)              :: vel        !< Cell-avg. velocity
            REAL(KIND(0d0))                                   :: pres       !< Cell-avg. pressure
            REAL(KIND(0d0)), DIMENSION(num_fluids)            :: alpha      !< Cell-avg. volume fraction
            REAL(KIND(0d0))                                   :: gamma      !< Cell-avg. sp. heat ratio
            REAL(KIND(0d0))                                   :: pi_inf     !< Cell-avg. liquid stiffness function
            REAL(KIND(0d0))                                   :: c          !< Cell-avg. sound speed
            REAL(KIND(0d0)), DIMENSION(2)                     :: Re         !< Cell-avg. Reynolds numbers
            REAL(KIND(0d0)), DIMENSION(num_fluids,num_fluids) :: We         !< Cell-avg. Weber numbers

            ! ICFL, VCFL, CCFL and Rc stability criteria extrema for the current
            ! time-step and located on both the local (loc) and the global (glb)
            ! computational domains
            REAL(KIND(0d0)) :: icfl_max_loc, icfl_max_glb !< ICFL stability extrema on local and global grids
            REAL(KIND(0d0)) :: vcfl_max_loc, vcfl_max_glb !< VCFL stability extrema on local and global grids
            REAL(KIND(0d0)) :: ccfl_max_loc, ccfl_max_glb !< CCFL stability extrema on local and global grids
            REAL(KIND(0d0)) ::   Rc_min_loc,   Rc_min_glb !< Rc   stability extrema on local and global grids


            REAL(KIND(0d0)) :: blkmod1, blkmod2 !<
            !! Fluid bulk modulus for Woods mixture sound speed
            

            INTEGER :: i,j,k,l !< Generic loop iterators

            INTEGER :: Nfq
            REAL(KIND(0d0)) :: fltr_dtheta   !< 
            !! Modified dtheta accounting for Fourier filtering in azimuthal direction. 
            

            ! Computing Stability Criteria at Current Time-step ================
            DO l = 0, p
               DO k = 0, n
                  DO j = 0, m
                     
                     DO i = 1, crv_size
                        alpha_rho(crv_idx(i)) = q_prim_vf(crv_idx(i))%sf(j,k,l)
                     END DO
                     

                     CALL s_convert_to_mixture_variables( q_prim_vf, rho, &
                                                          gamma, pi_inf,  &
                                                          Re, We, j,k,l   )
                     
                     DO i = 1, num_dims
                        vel(i) = q_prim_vf(cont_idx%end+i)%sf(j,k,l)
                     END DO
                     
                     pres = q_prim_vf(E_idx)%sf(j,k,l)
                     
                     ! Compute mixture sound speed
                     IF ( alt_soundspeed ) THEN
                         DO i = 1, num_fluids
                             alpha(i) = q_prim_vf(E_idx+i)%sf(j,k,l)
                         END DO
                         blkmod1 = ((fluid_pp(1)%gamma +1d0)*pres + &
                             fluid_pp(1)%pi_inf)/fluid_pp(1)%gamma
                         blkmod2 = ((fluid_pp(2)%gamma +1d0)*pres + &
                             fluid_pp(2)%pi_inf)/fluid_pp(2)%gamma
                         c = (1d0/(rho*(alpha(1)/blkmod1 + alpha(2)/blkmod2)))
                     ELSEIF(model_eqns == 3) THEN
                        c = 0d0
                        DO i = 1, num_fluids
                            c = c + q_prim_vf(i+adv_idx%beg-1)%sf(j,k,l) * (1d0/fluid_pp(i)%gamma+1d0) * &
                                (pres + fluid_pp(i)%pi_inf/(fluid_pp(i)%gamma+1d0))
                        END DO
                     ELSE
                         DO i = 1, crv_size
                             alpha(crv_idx(i))= q_prim_vf(E_idx+crv_idx(i))%sf(j,k,l)
                         END DO
                         c = (((gamma + 1d0)*pres + pi_inf)/(gamma*rho))
                     END IF

                     IF (mixture_err .AND. c < 0d0) THEN
                         c = sgm_eps
                     ELSE
                         c = SQRT(c)
                     END IF

                     IF (grid_geometry == 3) THEN
                        IF (k == 0) THEN
                            fltr_dtheta = 2d0*pi*y_cb(0)/3d0
                        ELSEIF (k <= fourier_rings) THEN
                            Nfq = MIN(FLOOR(2d0*REAL(k,KIND(0d0))*pi),(p+1)/2+1)
                            fltr_dtheta = 2d0*pi*y_cb(k-1)/REAL(Nfq,KIND(0d0))
                        ELSE
                            fltr_dtheta = y_cb(k-1)*dz(l)
                        END IF
                     END IF

                     IF(p > 0) THEN
                        !3D
                        IF (grid_geometry == 3) THEN
                            icfl_sf(j,k,l) = dt / MIN(       dx(j)/(ABS(vel(1))+c), &
                                                             dy(k)/(ABS(vel(2))+c), &
                                                       fltr_dtheta/(ABS(vel(3))+c)  )
                        ELSE
                            icfl_sf(j,k,l) = dt / MIN( dx(j)/(ABS(vel(1))+c), &
                                                       dy(k)/(ABS(vel(2))+c), &
                                                       dz(l)/(ABS(vel(3))+c)  )
                        END IF
                        
                        IF(ANY(Re_size > 0)) THEN
                           
                            IF (grid_geometry == 3) THEN
                               vcfl_sf(j,k,l) = MAXVAL(dt/Re) &
                                              / MIN(dx(j),dy(k),fltr_dtheta)**2d0
                               
                                 Rc_sf(j,k,l) = MIN(       dx(j)*(ABS(vel(1))+c),   &
                                                           dy(k)*(ABS(vel(2))+c),   &
                                                     fltr_dtheta*(ABS(vel(3))+c)  ) &
                                              / MAXVAL(1d0/Re)
                            ELSE
                               vcfl_sf(j,k,l) = MAXVAL(dt/Re) &
                                              / MIN(dx(j),dy(k),dz(l))**2d0
                               
                                 Rc_sf(j,k,l) = MIN( dx(j)*(ABS(vel(1))+c),   &
                                                     dy(k)*(ABS(vel(2))+c),   &
                                                     dz(l)*(ABS(vel(3))+c)  ) &
                                              / MAXVAL(1d0/Re)
                            END IF
                           
                        END IF
                        
                        IF(We_size > 0) THEN
                           
                           ccfl_sf(j,k,l) = 0d0
                           
                           DO i = 1, We_size
                              ccfl_sf(j,k,l) = MAX( ccfl_sf(j,k,l) , &
                                 1d0 / We(We_idx(i,1),We_idx(i,2)) / &
                                      ( MAX(alpha_rho(We_idx(i,1)) , &
                                                              0d0) / &
                                           MAX( alpha(We_idx(i,1)) , &
                                                          sgm_eps)   &
                                      + MAX(alpha_rho(We_idx(i,2)) , &
                                                              0d0) / &
                                           MAX( alpha(We_idx(i,2)) , &
                                                          sgm_eps) ) )
                           END DO
                           
                           IF (grid_geometry == 3) THEN
                               ccfl_sf(j,k,l) = SQRT( pi*ccfl_sf(j,k,l)*dt**2d0 &
                                                           / MIN(       dx(j),        &
                                                                        dy(k),        &
                                                                  fltr_dtheta  )**3d0 )
                           ELSE
                               ccfl_sf(j,k,l) = SQRT( pi*ccfl_sf(j,k,l)*dt**2d0 &
                                                           / MIN( dx(j),        &
                                                                  dy(k),        &
                                                                  dz(l)  )**3d0 )
                           END IF
                           
                        END IF
                        
                     ELSEIF(n > 0) THEN
                        !2D
                        icfl_sf(j,k,l) = dt / MIN( dx(j)/(ABS(vel(1))+c), &
                                                   dy(k)/(ABS(vel(2))+c)  )
                        
                        IF(ANY(Re_size > 0)) THEN
                           
                           vcfl_sf(j,k,l) = MAXVAL(dt/Re)/MIN(dx(j),dy(k))**2d0
                           
                             Rc_sf(j,k,l) = MIN( dx(j)*(ABS(vel(1))+c),   &
                                                 dy(k)*(ABS(vel(2))+c)  ) &
                                          / MAXVAL(1d0/Re)
                           
                        END IF
                        
                        IF(We_size > 0) THEN
                           
                           ccfl_sf(j,k,l) = 0d0
                           
                           DO i = 1, We_size
                              ccfl_sf(j,k,l) = MAX( ccfl_sf(j,k,l) , &
                                 1d0 / We(We_idx(i,1),We_idx(i,2)) / &
                                      ( MAX(alpha_rho(We_idx(i,1)) , &
                                                              0d0) / &
                                           MAX( alpha(We_idx(i,1)) , &
                                                          sgm_eps)   &
                                      + MAX(alpha_rho(We_idx(i,2)) , &
                                                              0d0) / &
                                           MAX( alpha(We_idx(i,2)) , &
                                                          sgm_eps) ) )
                           END DO
                           
                           ccfl_sf(j,k,l) = SQRT( pi*ccfl_sf(j,k,l)*dt**2d0 &
                                                    / MIN(dx(j),dy(k))**3d0 )
                           
                        END IF
                        
                     ELSE
                        !1D 
                        icfl_sf(j,k,l) = (dt/dx(j))*(ABS(vel(1))+c)
                        
                        IF(ANY(Re_size > 0)) THEN
                           
                           vcfl_sf(j,k,l) = MAXVAL(dt/Re)/dx(j)**2d0
                           
                             Rc_sf(j,k,l) = dx(j)*(ABS(vel(1))+c)/MAXVAL(1d0/Re)
                           
                        END IF
                        
                     END IF
                     
                  END DO
               END DO
            END DO
            ! END: Computing Stability Criteria at Current Time-step ===========


            
            ! Determining local stability criteria extrema at current time-step
                                 icfl_max_loc = MAXVAL(icfl_sf)
            IF(ANY(Re_size > 0)) vcfl_max_loc = MAXVAL(vcfl_sf)
            IF(    We_size > 0 ) ccfl_max_loc = MAXVAL(ccfl_sf)
            IF(ANY(Re_size > 0))   Rc_min_loc = MINVAL(  Rc_sf)
            
            
            ! Determining global stability criteria extrema at current time-step
            IF(num_procs > 1) THEN
                CALL s_mpi_reduce_stability_criteria_extrema( icfl_max_loc, &
                                                              vcfl_max_loc, &
                                                              ccfl_max_loc, &
                                                                Rc_min_loc, &
                                                              icfl_max_glb, &
                                                              vcfl_max_glb, &
                                                              ccfl_max_glb, &
                                                                Rc_min_glb  )
            ELSE
                                     icfl_max_glb = icfl_max_loc
                IF(ANY(Re_size > 0)) vcfl_max_glb = vcfl_max_loc
                IF(    We_size > 0 ) ccfl_max_glb = ccfl_max_loc
                IF(ANY(Re_size > 0))   Rc_min_glb =   Rc_min_loc
            END IF
            
            ! Determining the stability criteria extrema over all the time-steps
            IF(icfl_max_glb > icfl_max) icfl_max = icfl_max_glb
            
            IF(ANY(Re_size > 0)) THEN
                IF(vcfl_max_glb > vcfl_max) vcfl_max = vcfl_max_glb
                IF(  Rc_min_glb <   Rc_min)   Rc_min =   Rc_min_glb
            END IF
            
            IF(We_size > 0) THEN
                IF(ccfl_max_glb > ccfl_max) ccfl_max = ccfl_max_glb
            END IF
            
            
            ! Outputting global stability criteria extrema at current time-step
            IF(proc_rank == 0) THEN
                IF(ANY(Re_size > 0) .AND. We_size > 0) THEN
                    WRITE(1,'(4X,I8,4X,F10.6,4X,F9.6,4X,F9.6,4X,F9.6,4X,F10.6)') &
                                                t_step, t_step*dt, icfl_max_glb, &
                                                                   vcfl_max_glb, &
                                                                   ccfl_max_glb, &
                                                                     Rc_min_glb
                ELSEIF(ANY(Re_size > 0)) THEN
                    WRITE(1,'(6X,I8,6X,F10.6,6X,F9.6,6X,F9.6,6X,F10.6)') &
                                        t_step, t_step*dt, icfl_max_glb, &
                                                           vcfl_max_glb, &
                                                             Rc_min_glb
                ELSEIF(We_size > 0) THEN
                    WRITE(1,'(9X,I8,9X,F10.6,9X,F9.6,9X,F9.6)') &
                               t_step, t_step*dt, icfl_max_glb, &
                                                  ccfl_max_glb
                ELSE
                    WRITE(1,'(13X,I8,14X,F10.6,13X,F9.6)') &
                          t_step, t_step*dt, icfl_max_glb
                END IF

                IF (icfl_max_glb /= icfl_max_glb) THEN
                    PRINT '(A)', 'ICFL is NaN. Exiting ...'
                    ! print*, (dt/dx(:)),ABS(vel(1)),c

                    CALL s_mpi_abort()
                ELSEIF (icfl_max_glb > 1d0) THEN
                    PRINT '(A)', 'ICFL is greater than 1.0. Exiting ...'
                    PRINT*, 'icfl', icfl_max_glb
                    CALL s_mpi_abort()
                END IF
            END IF
            
            CALL s_mpi_barrier()
            
        END SUBROUTINE s_write_run_time_information ! --------------------------
        
        
        
        !>  The goal of this subroutine is to output the grid and
        !!      conservative variables data files for given time-step.        
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param t_step Current time-step
        SUBROUTINE s_write_serial_data_files(q_cons_vf, t_step) ! ---------------------

            
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q_cons_vf
            INTEGER, INTENT(IN) :: t_step
            
            CHARACTER(LEN = path_len + 2*name_len) :: t_step_dir !<
            !! Relative path to the current time-step directory


            CHARACTER(LEN = path_len + 3*name_len) :: file_path !<
            !! Relative path to the grid and conservative variables data files
            

            LOGICAL :: file_exist !<
            !! Logical used to check existence of current time-step directory
            
            CHARACTER(LEN=15) :: FMT

            INTEGER :: i, j, k, l, ii !< Generic loop iterators

            REAL(KIND(0d0)), DIMENSION(nb) :: nRtmp         !< Temporary bubble concentration
            REAL(KIND(0d0)) :: nbub                         !< Temporary bubble number density
            REAL(KIND(0d0)) :: gamma, lit_gamma, pi_inf     !< Temporary EOS params
            REAL(KIND(0d0)) :: rho                          !< Temporary density
            REAL(KIND(0d0)), DIMENSION(2)                   :: Re !< Temporary Reynolds number
            REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:)    :: We !< Temporary Weber number

            ! Creating or overwriting the time-step root directory
            WRITE(t_step_dir,'(A,I0,A,I0)') TRIM(case_dir) // '/p_all'
 
            ! Creating or overwriting the current time-step directory
            WRITE(t_step_dir,'(A,I0,A,I0)') TRIM(case_dir) // '/p_all/p', &
                                            proc_rank, '/', t_step
            
            file_path = TRIM(t_step_dir) // '/.'
            CALL my_inquire(file_path,file_exist)
            IF(file_exist) CALL SYSTEM('rm -rf ' // TRIM(t_step_dir))
            CALL SYSTEM('mkdir -p ' // TRIM(t_step_dir))
            
            ! Writing the grid data file in the x-direction
            file_path = TRIM(t_step_dir) // '/x_cb.dat'
            
            OPEN(2, FILE   = TRIM(file_path), &
                    FORM   = 'unformatted'  , &
                    STATUS = 'new'            )
            WRITE(2) x_cb(-1:m); CLOSE(2)
            
            ! Writing the grid data files in the y- and z-directions
            IF(n > 0) THEN
                
                file_path = TRIM(t_step_dir) // '/y_cb.dat'
                
                OPEN(2, FILE   = TRIM(file_path), &
                        FORM   = 'unformatted'  , &
                        STATUS = 'new'            )
                WRITE(2) y_cb(-1:n); CLOSE(2)
                
                IF(p > 0) THEN
                    
                    file_path = TRIM(t_step_dir) // '/z_cb.dat'
                    
                    OPEN(2, FILE   = TRIM(file_path), &
                            FORM   = 'unformatted'  , &
                            STATUS = 'new'            )
                    WRITE(2) z_cb(-1:p); CLOSE(2)
                    
                END IF
                
            END IF
            
            ! Writing the conservative variables data files
            DO i = 1, sys_size
                WRITE(file_path,'(A,I0,A)') TRIM(t_step_dir) // '/q_cons_vf', &
                                            i, '.dat'

                    OPEN(2, FILE   = TRIM(file_path), &
                        FORM   = 'unformatted'  , &
                        STATUS = 'new'            )
                
                    IF (i == E_idx .AND. We_size > 0 .AND. (We_riemann_flux .OR. We_rhs_flux)) THEN
                        CALL s_remove_capillary_potential_energy(q_cons_vf)
                        WRITE(2) energy(0:m,0:n,0:p); CLOSE(2)
                    ELSE
                        WRITE(2) q_cons_vf(i)%sf(0:m,0:n,0:p); CLOSE(2)
                    END IF  
            END DO

            gamma = fluid_pp(1)%gamma
            lit_gamma = 1d0/fluid_pp(1)%gamma + 1d0
            pi_inf = fluid_pp(1)%pi_inf

            IF (precision==1) THEN
                FMT="(2F30.7)"
            ELSE
                FMT="(2F40.14)"
            END IF

            ! writting an output directory
            WRITE(t_step_dir,'(A,I0,A,I0)') TRIM(case_dir) // '/D'
            file_path = TRIM(t_step_dir) // '/.'
        
            INQUIRE( FILE = TRIM(file_path), EXIST = file_exist       )
        
            IF(.NOT.file_exist) CALL SYSTEM('mkdir -p ' // TRIM(t_step_dir))

            !1D
            IF (n==0 .AND. p==0) THEN

                IF (model_eqns==2) THEN
                    DO i = 1, sys_size
        WRITE(file_path,'(A,I0,A,I2.2,A,I6.6,A)') TRIM(t_step_dir) // '/prim.', i, '.', proc_rank, '.', t_step,'.dat'

                        OPEN(2,FILE= TRIM(file_path) )
                            DO j=0,m
                                CALL s_convert_to_mixture_variables( q_cons_vf, rho, gamma, pi_inf, Re, We, j,0,0)
                                lit_gamma = 1d0/gamma + 1d0
                                
                                IF ( ((i.ge.cont_idx%beg) .AND. (i.le.cont_idx%end))    &
                                                          .OR.                          &
                                     ((i.ge.adv_idx%beg)  .AND. (i.le.adv_idx%end))     &
                                                         ) THEN 
                                    WRITE(2,FMT) x_cb(j),q_cons_vf(i)%sf(j,0,0)
                                ELSE IF (i.eq.mom_idx%beg) THEN !u
                                    WRITE(2,FMT) x_cb(j),q_cons_vf(mom_idx%beg)%sf(j,0,0)/rho
                                ELSE IF (i.eq.stress_idx%beg) THEN !tau_e
                                    WRITE(2,FMT) x_cb(j),q_cons_vf(stress_idx%beg)%sf(j,0,0)/rho
                                ELSE IF (i.eq.E_idx) THEN !p
                                    IF (model_eqns == 4) THEN
                                        !Tait pressure from density
                                        WRITE(2,FMT) x_cb(j), &
                                            (pref + pi_inf) * (                     &
                                            ( q_cons_vf(1)%sf(j,0,0)/               &
                                            (rhoref*(1.d0-q_cons_vf(4)%sf(j,0,0)))  & 
                                            ) ** lit_gamma )                        &
                                            - pi_inf
                                    ELSE IF (model_eqns == 2 .AND. (bubbles .NEQV. .TRUE.)) THEN
                                        !Stiffened gas pressure from energy
                                        WRITE(2,FMT) x_cb(j), &
                                            (                                       & 
                                            q_cons_vf(E_idx)%sf(j,0,0)  -            &
                                            0.5d0*(q_cons_vf(mom_idx%beg)%sf(j,0,0)**2.d0)/rho - &
                                            pi_inf &
                                            ) / gamma
                                    ELSE
                                        !Stiffened gas pressure from energy with bubbles
                                        WRITE(2,FMT) x_cb(j), &
                                            (                                       & 
                                            (q_cons_vf(E_idx)%sf(j,0,0)  -            &
                                            0.5d0*(q_cons_vf(mom_idx%beg)%sf(j,0,0)**2.d0)/rho) / &
                                            (1.d0 - q_cons_vf(alf_idx)%sf(j,0,0)) - &
                                            pi_inf &
                                            ) / gamma
                                    END IF
                                ELSE IF ((i >= bub_idx%beg) .AND. (i <= bub_idx%end) .AND. bubbles) THEN
                                    DO k = 1,nb
                                        nRtmp(k) = q_cons_vf(bub_idx%rs(k))%sf(j,0,0)
                                    END DO
                                    CALL s_comp_n_from_cons( q_cons_vf(alf_idx)%sf(j,0,0), nRtmp, nbub) 
                                    
                                    WRITE(2,FMT) x_cb(j),q_cons_vf(i)%sf(j,0,0)/nbub
                                END IF
                            END DO
                        CLOSE(2)
                    END DO
                END IF

                DO i = 1, sys_size    
        WRITE(file_path,'(A,I0,A,I2.2,A,I6.6,A)') TRIM(t_step_dir) // '/cons.', i, '.', proc_rank, '.', t_step,'.dat'

                    OPEN(2,FILE= TRIM(file_path) )
                        DO j=0,m
                            WRITE(2,FMT) x_cb(j),q_cons_vf(i)%sf(j,0,0)
                        END DO
                    CLOSE(2)
                END DO
            END IF


            IF (precision==1) THEN
                FMT="(3F30.7)"
            ELSE
                FMT="(3F40.14)"
            END IF

            ! 2D
            IF ( (n>0) .AND. (p==0) ) THEN
                DO i = 1,sys_size
        WRITE(file_path,'(A,I0,A,I2.2,A,I6.6,A)') TRIM(t_step_dir) // '/cons.', i, '.', proc_rank, '.', t_step,'.dat'
                    OPEN(2,FILE= TRIM(file_path) )
                        DO j=0,m
                        DO k=0,n
                            WRITE(2,FMT) x_cb(j),y_cb(k), q_cons_vf(i)%sf(j,k,0)
                        END DO
                        WRITE(2,*)
                        END DO
                    CLOSE(2)
                END DO
            END IF


            IF (precision==1) THEN
                FMT="(4F30.7)"
            ELSE
                FMT="(4F40.14)"
            END IF

            ! 3D
            IF ( p > 0) THEN
                DO i = 1,sys_size
        WRITE(file_path,'(A,I0,A,I2.2,A,I6.6,A)') TRIM(t_step_dir) // '/cons.', i, '.', proc_rank, '.', t_step,'.dat'
                    OPEN(2,FILE= TRIM(file_path) )
                        DO j=0,m
                        DO k=0,n
                        DO l=0,p
                            WRITE(2,FMT) x_cb(j),y_cb(k),z_cb(l), q_cons_vf(i)%sf(j,k,l)
                        END DO
                        WRITE(2,*)
                        END DO
                        WRITE(2,*)
                        END DO
                    CLOSE(2)
                END DO
            END IF

        END SUBROUTINE s_write_serial_data_files ! ------------------------------------
        
        
        SUBROUTINE s_remove_capillary_potential_energy(v_vf)

            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: v_vf
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
            
            ! Remove capillary potential energy from conservative variable
            energy(:,:,:) = v_vf(E_idx)%sf(:,:,:) - E_We(:,:,:)
            
        END SUBROUTINE s_remove_capillary_potential_energy
        
        
        
        !>  The goal of this subroutine is to output the grid and
        !!      conservative variables data files for given time-step.        
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param t_step Current time-step
        SUBROUTINE s_write_parallel_data_files(q_cons_vf, t_step) ! --

            TYPE(scalar_field), &
            DIMENSION(sys_size), &
            INTENT(IN) :: q_cons_vf

            INTEGER, INTENT(IN) :: t_step

            INTEGER :: ifile, ierr, data_size
            INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
            INTEGER(KIND=MPI_OFFSET_KIND) :: disp
            INTEGER(KIND=MPI_OFFSET_KIND) :: m_MOK, n_MOK, p_MOK
            INTEGER(KIND=MPI_OFFSET_KIND) :: WP_MOK, var_MOK, str_MOK
            INTEGER(KIND=MPI_OFFSET_KIND) :: NVARS_MOK
            INTEGER(KIND=MPI_OFFSET_KIND) :: MOK

            CHARACTER(LEN=path_len + 2*name_len) :: file_loc 
            LOGICAL :: file_exist


            INTEGER :: i !< Generic loop iterator

            ! Initialize MPI data I/O
            CALL s_initialize_mpi_data(q_cons_vf)

            IF (We_size > 0 .AND. (We_riemann_flux .OR. We_rhs_flux)) THEN
                CALL s_remove_capillary_potential_energy(q_cons_vf)
                MPI_IO_DATA%var(E_idx)%sf => energy(0:m,0:n,0:p)
            END IF

            ! Open the file to write all flow variables
            WRITE(file_loc, '(I0,A)') t_step, '.dat'
            file_loc = TRIM(case_dir) // '/restart_data' // TRIM(mpiiofs) // TRIM(file_loc)
            INQUIRE(FILE = TRIM(file_loc),EXIST = file_exist)
            IF (file_exist .AND. proc_rank == 0) THEN
                CALL MPI_FILE_DELETE(file_loc,mpi_info_int,ierr)
            END IF
            CALL MPI_FILE_OPEN(MPI_COMM_WORLD,file_loc,IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),&
                        mpi_info_int,ifile,ierr)

            ! Size of local arrays
            data_size = (m+1)*(n+1)*(p+1)

            ! Resize some integers so MPI can write even the biggest files
            m_MOK     = INT(m_glb+1,     MPI_OFFSET_KIND)
            n_MOK     = INT(n_glb+1,     MPI_OFFSET_KIND)
            p_MOK     = INT(p_glb+1,     MPI_OFFSET_KIND)
            WP_MOK    = INT(8d0,      MPI_OFFSET_KIND)
            MOK       = INT(1d0,      MPI_OFFSET_KIND)
            str_MOK   = INT(name_len, MPI_OFFSET_KIND)
            NVARS_MOK = INT(sys_size, MPI_OFFSET_KIND)
            
            IF (bubbles) THEN
                ! Write the data for each variable
                DO i = 1, sys_size
                    var_MOK = INT(i, MPI_OFFSET_KIND)

                    ! Initial displacement to skip at beginning of file
                    disp = m_MOK*MAX(MOK,n_MOK)*MAX(MOK,p_MOK)*WP_MOK*(var_MOK-1)

                    CALL MPI_FILE_SET_VIEW(ifile,disp,MPI_DOUBLE_PRECISION,MPI_IO_DATA%view(i), &
                                'native',mpi_info_int,ierr)
                    CALL MPI_FILE_WRITE_ALL(ifile,MPI_IO_DATA%var(i)%sf,data_size, &
                                MPI_DOUBLE_PRECISION,status,ierr)
                END DO
            ELSE
                DO i = 1, adv_idx%end
                    var_MOK = INT(i, MPI_OFFSET_KIND)

                    ! Initial displacement to skip at beginning of file
                    disp = m_MOK*MAX(MOK,n_MOK)*MAX(MOK,p_MOK)*WP_MOK*(var_MOK-1)

                    CALL MPI_FILE_SET_VIEW(ifile,disp,MPI_DOUBLE_PRECISION,MPI_IO_DATA%view(i), &
                                'native',mpi_info_int,ierr)
                    CALL MPI_FILE_WRITE_ALL(ifile,MPI_IO_DATA%var(i)%sf,data_size, &
                                MPI_DOUBLE_PRECISION,status,ierr)
                END DO
            END IF

            CALL MPI_FILE_CLOSE(ifile,ierr)

        END SUBROUTINE s_write_parallel_data_files ! ---------------------------

        !>  This writes a formatted data file where the root processor
        !!      can write out the CoM information    
        !!  @param t_step Current time-step
        !!  @param q_com Center of mass information
        !!  @param moments Higher moment information
        SUBROUTINE s_write_com_files(t_step,q_com,moments) ! -------------------

            INTEGER, INTENT(IN) :: t_step
            REAL(KIND(0d0)), DIMENSION(num_fluids,10), INTENT(IN) :: q_com
            REAL(KIND(0d0)), DIMENSION(num_fluids,2,5), INTENT(IN) :: moments


            INTEGER :: i !< Generic loop iterator
            REAL(KIND(0d0)) :: nondim_time !< Non-dimensional time

            ! Non-dimensional time calculation
            IF (t_step_old /= dflt_int) THEN
                nondim_time = REAL(t_step + t_step_old,KIND(0d0))*dt
            ELSE
                nondim_time = REAL(t_step,KIND(0d0))*dt
            END IF

            IF (n == 0) THEN ! 1D simulation
                DO i = 1, num_fluids ! Loop through fluids
                    IF (com_wrt(i)) THEN ! Writing out CoM data
                        IF (proc_rank == 0) THEN
                            WRITE(i+10, '(6X,F12.6,F24.8,F24.8,F24.8,F24.8)') &
                                nondim_time, &
                                q_com(i,1), &
                                q_com(i,2), &
                                q_com(i,5), &
                                q_com(i,8)
                                

                        END IF
                    END IF
                END DO
            ELSEIF (p == 0) THEN ! 2D simulation
                DO i = 1, num_fluids ! Loop through fluids
                    IF (com_wrt(i)) THEN ! Writing out CoM data
                        IF (proc_rank == 0) THEN
                            IF (moment_order(1) == dflt_int) THEN
                                WRITE(i+10, '(6X,F12.6,F24.8,F24.8,F24.8,F24.8,' // &
                                             'F24.8,F24.8,F24.8)') &
                                    nondim_time, &
                                    q_com(i,1), &
                                    q_com(i,2), &
                                    q_com(i,3), &
                                    q_com(i,5), &
                                    q_com(i,6), &
                                    q_com(i,8), &
                                    q_com(i,9)
                            ELSEIF (moment_order(2) == dflt_int) THEN
                                WRITE(i+10, '(6X,F12.6,F24.8,F24.8,F24.8,F24.8,' // &
                                             'F24.8,F24.8,F24.8,E24.8)') &
                                    nondim_time, &
                                    q_com(i,1), &
                                    q_com(i,2), &
                                    q_com(i,3), &
                                    q_com(i,5), &
                                    q_com(i,6), &
                                    q_com(i,8), &
                                    q_com(i,9), &
                                    moments(i,1,1)
                            ELSEIF (moment_order(3) == dflt_int) THEN
                                WRITE(i+10, '(6X,F12.6,F24.8,F24.8,F24.8,F24.8,' // &
                                             'F24.8,F24.8,F24.8,E24.8,' // &
                                             'E24.8)') &
                                    nondim_time, &
                                    q_com(i,1), &
                                    q_com(i,2), &
                                    q_com(i,3), &
                                    q_com(i,5), &
                                    q_com(i,6), &
                                    q_com(i,8), &
                                    q_com(i,9), &
                                    moments(i,1,1), &
                                    moments(i,1,2)
                            ELSEIF (moment_order(4) == dflt_int) THEN
                                WRITE(i+10, '(6X,F12.6,F24.8,F24.8,F24.8,F24.8,' // &
                                             'F24.8,F24.8,F24.8,E24.8,' // &
                                             'E24.8,E24.8)') &
                                    nondim_time, &
                                    q_com(i,1), &
                                    q_com(i,2), &
                                    q_com(i,3), &
                                    q_com(i,5), &
                                    q_com(i,6), &
                                    q_com(i,8), &
                                    q_com(i,9), &
                                    moments(i,1,1), &
                                    moments(i,1,2), &
                                    moments(i,1,3)
                            ELSEIF (moment_order(5) == dflt_int) THEN
                                WRITE(i+10, '(6X,F12.6,F24.8,F24.8,F24.8,F24.8,' // &
                                             'F24.8,F24.8,F24.8,E24.8,' // &
                                             'E24.8,E24.8,E24.8)') &
                                    nondim_time, &
                                    q_com(i,1), &
                                    q_com(i,2), &
                                    q_com(i,3), &
                                    q_com(i,5), &
                                    q_com(i,6), &
                                    q_com(i,8), &
                                    q_com(i,9), &
                                    moments(i,1,1), &
                                    moments(i,1,2), &
                                    moments(i,1,3), &
                                    moments(i,1,4)
                            ELSE
                                WRITE(i+10, '(6X,F12.6,F24.8,F24.8,F24.8,F24.8,' // &
                                             'F24.8,F24.8,F24.8,E24.8,' // &
                                             'E24.8,E24.8,E24.8,E24.8)') &
                                    nondim_time, &
                                    q_com(i,1), &
                                    q_com(i,2), &
                                    q_com(i,3), &
                                    q_com(i,5), &
                                    q_com(i,6), &
                                    q_com(i,8), &
                                    q_com(i,9), &
                                    moments(i,1,1), &
                                    moments(i,1,2), &
                                    moments(i,1,3), &
                                    moments(i,1,4), &
                                    moments(i,1,5)
                            END IF
                        END IF
                    END IF
                END DO
            ELSE ! 3D simulation
                DO i = 1, num_fluids ! Loop through fluids
                    IF (com_wrt(i)) THEN ! Writing out CoM data
                        IF (proc_rank == 0) THEN
                            IF (moment_order(1) == dflt_int) THEN
                                WRITE(i+10, '(6X,F12.6,F24.8,F24.8,F24.8,F24.8,' // &
                                             'F24.8,F24.8,F24.8,F24.8,' // &
                                             'F24.8,F24.8)') &
                                    nondim_time, &
                                    q_com(i,1), &
                                    q_com(i,2), &
                                    q_com(i,3), &
                                    q_com(i,4), &
                                    q_com(i,5), &
                                    q_com(i,6), &
                                    q_com(i,7), &
                                    q_com(i,8), &
                                    q_com(i,9), &
                                    q_com(i,10)
                            ELSEIF (moment_order(2) == dflt_int) THEN
                                WRITE(i+10, '(6X,F12.6,F24.8,F24.8,F24.8,F24.8,' // &
                                             'F24.8,F24.8,F24.8,F24.8,' // &
                                             'F24.8,F24.8,F24.8,F24.8)') &
                                    nondim_time, &
                                    q_com(i,1), &
                                    q_com(i,2), &
                                    q_com(i,3), &
                                    q_com(i,4), &
                                    q_com(i,5), &
                                    q_com(i,6), &
                                    q_com(i,7), &
                                    q_com(i,8), &
                                    q_com(i,9), &
                                    q_com(i,10), &
                                    moments(i,1,1), &
                                    moments(i,2,1)
                            ELSEIF (moment_order(3) == dflt_int) THEN
                                WRITE(i+10, '(6X,F12.6,F24.8,F24.8,F24.8,F24.8,' // &
                                             'F24.8,F24.8,F24.8,F24.8,' // &
                                             'F24.8,F24.8,F24.8,F24.8,' // &
                                             'F24.8,F24.8)') &
                                    nondim_time, &
                                    q_com(i,1), &
                                    q_com(i,2), &
                                    q_com(i,3), &
                                    q_com(i,4), &
                                    q_com(i,5), &
                                    q_com(i,6), &
                                    q_com(i,7), &
                                    q_com(i,8), &
                                    q_com(i,9), &
                                    q_com(i,10), &
                                    moments(i,1,1), &
                                    moments(i,1,2), &
                                    moments(i,2,1), &
                                    moments(i,2,2)
                            ELSEIF (moment_order(4) == dflt_int) THEN
                                WRITE(i+10, '(6X,F12.6,F24.8,F24.8,F24.8,F24.8,' // &
                                             'F24.8,F24.8,F24.8,F24.8,' // &
                                             'F24.8,F24.8,F24.8,F24.8,' // &
                                             'F24.8,F24.8,F24.8,F24.8)') &
                                    nondim_time, &
                                    q_com(i,1), &
                                    q_com(i,2), &
                                    q_com(i,3), &
                                    q_com(i,4), &
                                    q_com(i,5), &
                                    q_com(i,6), &
                                    q_com(i,7), &
                                    q_com(i,8), &
                                    q_com(i,9), &
                                    q_com(i,10), &
                                    moments(i,1,1), &
                                    moments(i,1,2), &
                                    moments(i,1,3), &
                                    moments(i,2,1), &
                                    moments(i,2,2), &
                                    moments(i,2,3)
                            ELSEIF (moment_order(5) == dflt_int) THEN
                                WRITE(i+10, '(6X,F12.6,F24.8,F24.8,F24.8,F24.8,' // &
                                             'F24.8,F24.8,F24.8,F24.8,' // &
                                             'F24.8,F24.8,F24.8,F24.8,' // &
                                             'F24.8,F24.8,F24.8,F24.8,' // &
                                             'F24.8,F24.8)') &
                                    nondim_time, &
                                    q_com(i,1), &
                                    q_com(i,2), &
                                    q_com(i,3), &
                                    q_com(i,4), &
                                    q_com(i,5), &
                                    q_com(i,6), &
                                    q_com(i,7), &
                                    q_com(i,8), &
                                    q_com(i,9), &
                                    q_com(i,10), &
                                    moments(i,1,1), &
                                    moments(i,1,2), &
                                    moments(i,1,3), &
                                    moments(i,1,4), &
                                    moments(i,2,1), &
                                    moments(i,2,2), &
                                    moments(i,2,3), &
                                    moments(i,2,4)
                            ELSE
                                WRITE(i+10, '(6X,F12.6,F24.8,F24.8,F24.8,F24.8,' // &
                                             'F24.8,F24.8,F24.8,F24.8,' // &
                                             'F24.8,F24.8,F24.8,F24.8,' // &
                                             'F24.8,F24.8,F24.8,F24.8,' // &
                                             'F24.8,F24.8,F24.8,F24.8)') &
                                    nondim_time, &
                                    q_com(i,1), &
                                    q_com(i,2), &
                                    q_com(i,3), &
                                    q_com(i,4), &
                                    q_com(i,5), &
                                    q_com(i,6), &
                                    q_com(i,7), &
                                    q_com(i,8), &
                                    q_com(i,9), &
                                    q_com(i,10), &
                                    moments(i,1,1), &
                                    moments(i,1,2), &
                                    moments(i,1,3), &
                                    moments(i,1,4), &
                                    moments(i,1,5), &
                                    moments(i,2,1), &
                                    moments(i,2,2), &
                                    moments(i,2,3), &
                                    moments(i,2,4), &
                                    moments(i,2,5)
                            END IF
                        END IF
                    END IF
                END DO
            END IF

        END SUBROUTINE s_write_com_files ! -------------------------------------



        !>  The goal of this subroutine is to output coherent body information.
        !!  @param t_step Current time-step
        !!  @param cb_mass Coherent body mass
        !!  @param bounds Coherent body boundary
        !!  @param cntrline Coherent body center line
        SUBROUTINE s_write_cb_files(t_step,cb_mass,bounds,cntrline) ! ----------

            INTEGER, INTENT(IN) :: t_step
            REAL(KIND(0d0)), DIMENSION(num_fluids,10), INTENT(IN) :: cb_mass
            REAL(KIND(0d0)), DIMENSION(num_fluids,5,6), INTENT(IN) :: bounds
            REAL(KIND(0d0)), DIMENSION(num_fluids,5), INTENT(IN) :: cntrline

            INTEGER :: i !< Generic loop iterator
            REAL(KIND(0d0)) :: nondim_time !< Non-dimensional time

            ! Non-dimensional time calculation
            IF (t_step_old /= dflt_int) THEN
                nondim_time = REAL(t_step + t_step_old,KIND(0d0))*dt
            ELSE
                nondim_time = REAL(t_step,KIND(0d0))*dt
            END IF

            DO i = 1, num_fluids ! Loop through fluids
                IF (cb_wrt(i)) THEN ! Writing out CoM data
                    IF (proc_rank == 0) THEN
                        IF (n == 0) THEN ! 1D simulation
                            IF (threshold_mf(2) == dflt_real) THEN
                                WRITE(i+20, '(6X,F12.6,F24.8,F24.8,F24.8,F24.8)') &
                                    nondim_time, &
                                    cb_mass(i,1), &
                                    cb_mass(i,6), &
                                    bounds(i,1,1), &
                                    bounds(i,1,2)
                            ELSEIF (threshold_mf(3) == dflt_real) THEN
                                WRITE(i+20, '(6X,F12.6,F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8)') &
                                    nondim_time, &
                                    cb_mass(i,1), &
                                    cb_mass(i,2), &
                                    cb_mass(i,6), &
                                    cb_mass(i,7), &
                                    bounds(i,1,1), &
                                    bounds(i,2,1), &
                                    bounds(i,1,2), &
                                    bounds(i,2,2)
                            ELSEIF (threshold_mf(4) == dflt_real) THEN
                                WRITE(i+20, '(6X,F12.6,F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8)') &
                                    nondim_time, &
                                    cb_mass(i,1), &
                                    cb_mass(i,2), &
                                    cb_mass(i,3), &
                                    cb_mass(i,6), &
                                    cb_mass(i,7), &
                                    cb_mass(i,8), &
                                    bounds(i,1,1), &
                                    bounds(i,2,1), &
                                    bounds(i,3,1), &
                                    bounds(i,1,2), &
                                    bounds(i,2,2), &
                                    bounds(i,3,2)
                            ELSEIF (threshold_mf(5) == dflt_real) THEN
                                WRITE(i+20, '(6X,F12.6,F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8)') &
                                    nondim_time, &
                                    cb_mass(i,1), &
                                    cb_mass(i,2), &
                                    cb_mass(i,3), &
                                    cb_mass(i,4), &
                                    cb_mass(i,6), &
                                    cb_mass(i,7), &
                                    cb_mass(i,8), &
                                    cb_mass(i,9), &
                                    bounds(i,1,1), &
                                    bounds(i,2,1), &
                                    bounds(i,3,1), &
                                    bounds(i,4,1), &
                                    bounds(i,1,2), &
                                    bounds(i,2,2), &
                                    bounds(i,3,2), &
                                    bounds(i,4,2)
                            ELSE
                                WRITE(i+20, '(6X,F12.6,F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8)') &
                                    nondim_time, &
                                    cb_mass(i,1), &
                                    cb_mass(i,2), &
                                    cb_mass(i,3), &
                                    cb_mass(i,4), &
                                    cb_mass(i,5), &
                                    cb_mass(i,6), &
                                    cb_mass(i,7), &
                                    cb_mass(i,8), &
                                    cb_mass(i,9), &
                                    cb_mass(i,10), &
                                    bounds(i,1,1), &
                                    bounds(i,2,1), &
                                    bounds(i,3,1), &
                                    bounds(i,4,1), &
                                    bounds(i,5,1), &
                                    bounds(i,1,2), &
                                    bounds(i,2,2), &
                                    bounds(i,3,2), &
                                    bounds(i,4,2), &
                                    bounds(i,5,2)
                            END IF
                        ELSEIF (p == 0) THEN ! 2D simulation
                            IF (threshold_mf(2) == dflt_real) THEN
                                WRITE(i+20, '(6X,F12.6,F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8)') &
                                    nondim_time, &
                                    cb_mass(i,1), &
                                    cb_mass(i,6), &
                                    bounds(i,1,1), &
                                    bounds(i,1,2), &
                                    bounds(i,1,3), &
                                    bounds(i,1,4), &
                                    cntrline(i,1)
                            ELSEIF (threshold_mf(3) == dflt_real) THEN
                                WRITE(i+20, '(6X,F12.6,F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8)') &
                                    nondim_time, &
                                    cb_mass(i,1), &
                                    cb_mass(i,2), &
                                    cb_mass(i,6), &
                                    cb_mass(i,7), &
                                    bounds(i,1,1), &
                                    bounds(i,2,1), &
                                    bounds(i,1,2), &
                                    bounds(i,2,2), &
                                    bounds(i,1,3), &
                                    bounds(i,2,3), &
                                    bounds(i,1,4), &
                                    bounds(i,2,4), &
                                    cntrline(i,1), &
                                    cntrline(i,2)
                            ELSEIF (threshold_mf(4) == dflt_real) THEN
                                WRITE(i+20, '(6X,F12.6,F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8)') &
                                    nondim_time, &
                                    cb_mass(i,1), &
                                    cb_mass(i,2), &
                                    cb_mass(i,3), &
                                    cb_mass(i,6), &
                                    cb_mass(i,7), &
                                    cb_mass(i,8), &
                                    bounds(i,1,1), &
                                    bounds(i,2,1), &
                                    bounds(i,3,1), &
                                    bounds(i,1,2), &
                                    bounds(i,2,2), &
                                    bounds(i,3,2), &
                                    bounds(i,1,3), &
                                    bounds(i,2,3), &
                                    bounds(i,3,3), &
                                    bounds(i,1,4), &
                                    bounds(i,2,4), &
                                    bounds(i,3,4), &
                                    cntrline(i,1), &
                                    cntrline(i,2), &
                                    cntrline(i,3)
                            ELSEIF (threshold_mf(5) == dflt_real) THEN
                                WRITE(i+20, '(6X,F12.6,F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8)') &
                                    nondim_time, &
                                    cb_mass(i,1), &
                                    cb_mass(i,2), &
                                    cb_mass(i,3), &
                                    cb_mass(i,4), &
                                    cb_mass(i,6), &
                                    cb_mass(i,7), &
                                    cb_mass(i,8), &
                                    cb_mass(i,9), &
                                    bounds(i,1,1), &
                                    bounds(i,2,1), &
                                    bounds(i,3,1), &
                                    bounds(i,4,1), &
                                    bounds(i,1,2), &
                                    bounds(i,2,2), &
                                    bounds(i,3,2), &
                                    bounds(i,4,2), &
                                    bounds(i,1,3), &
                                    bounds(i,2,3), &
                                    bounds(i,3,3), &
                                    bounds(i,4,3), &
                                    bounds(i,1,4), &
                                    bounds(i,2,4), &
                                    bounds(i,3,4), &
                                    bounds(i,4,4), &
                                    cntrline(i,1), &
                                    cntrline(i,2), &
                                    cntrline(i,3), &
                                    cntrline(i,4)
                            ELSE
                                WRITE(i+20, '(6X,F12.6,F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8)') &
                                    nondim_time, &
                                    cb_mass(i,1), &
                                    cb_mass(i,2), &
                                    cb_mass(i,3), &
                                    cb_mass(i,4), &
                                    cb_mass(i,5), &
                                    cb_mass(i,6), &
                                    cb_mass(i,7), &
                                    cb_mass(i,8), &
                                    cb_mass(i,9), &
                                    cb_mass(i,10), &
                                    bounds(i,1,1), &
                                    bounds(i,2,1), &
                                    bounds(i,3,1), &
                                    bounds(i,4,1), &
                                    bounds(i,5,1), &
                                    bounds(i,1,2), &
                                    bounds(i,2,2), &
                                    bounds(i,3,2), &
                                    bounds(i,4,2), &
                                    bounds(i,5,2), &
                                    bounds(i,1,3), &
                                    bounds(i,2,3), &
                                    bounds(i,3,3), &
                                    bounds(i,4,3), &
                                    bounds(i,5,3), &
                                    bounds(i,1,4), &
                                    bounds(i,2,4), &
                                    bounds(i,3,4), &
                                    bounds(i,4,4), &
                                    bounds(i,5,4), &
                                    cntrline(i,1), &
                                    cntrline(i,2), &
                                    cntrline(i,3), &
                                    cntrline(i,4), &
                                    cntrline(i,5)
                            END IF
                        ELSE ! 3D simulation
                            IF (threshold_mf(2) == dflt_real) THEN
                                WRITE(i+20, '(6X,F12.6,F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8)') &
                                    nondim_time, &
                                    cb_mass(i,1), &
                                    cb_mass(i,6), &
                                    bounds(i,1,1), &
                                    bounds(i,1,2), &
                                    bounds(i,1,3), &
                                    bounds(i,1,4), &
                                    bounds(i,1,5), &
                                    bounds(i,1,6), &
                                    cntrline(i,1)
                            ELSEIF (threshold_mf(3) == dflt_real) THEN
                                WRITE(i+20, '(6X,F12.6,F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8)') &
                                    nondim_time, &
                                    cb_mass(i,1), &
                                    cb_mass(i,2), &
                                    cb_mass(i,6), &
                                    cb_mass(i,7), &
                                    bounds(i,1,1), &
                                    bounds(i,2,1), &
                                    bounds(i,1,2), &
                                    bounds(i,2,2), &
                                    bounds(i,1,3), &
                                    bounds(i,2,3), &
                                    bounds(i,1,4), &
                                    bounds(i,2,4), &
                                    bounds(i,1,5), &
                                    bounds(i,2,5), &
                                    bounds(i,1,6), &
                                    bounds(i,2,6), &
                                    cntrline(i,1), &
                                    cntrline(i,2)
                            ELSEIF (threshold_mf(4) == dflt_real) THEN
                                WRITE(i+20, '(6X,F12.6,F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8)') &
                                    nondim_time, &
                                    cb_mass(i,1), &
                                    cb_mass(i,2), &
                                    cb_mass(i,3), &
                                    cb_mass(i,6), &
                                    cb_mass(i,7), &
                                    cb_mass(i,8), &
                                    bounds(i,1,1), &
                                    bounds(i,2,1), &
                                    bounds(i,3,1), &
                                    bounds(i,1,2), &
                                    bounds(i,2,2), &
                                    bounds(i,3,2), &
                                    bounds(i,1,3), &
                                    bounds(i,2,3), &
                                    bounds(i,3,3), &
                                    bounds(i,1,4), &
                                    bounds(i,2,4), &
                                    bounds(i,3,4), &
                                    bounds(i,1,5), &
                                    bounds(i,2,5), &
                                    bounds(i,3,5), &
                                    bounds(i,1,6), &
                                    bounds(i,2,6), &
                                    bounds(i,3,6), &
                                    cntrline(i,1), &
                                    cntrline(i,2), &
                                    cntrline(i,3)
                            ELSEIF (threshold_mf(5) == dflt_real) THEN
                                WRITE(i+20, '(6X,F12.6,F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8)') &
                                    nondim_time, &
                                    cb_mass(i,1), &
                                    cb_mass(i,2), &
                                    cb_mass(i,3), &
                                    cb_mass(i,4), &
                                    cb_mass(i,6), &
                                    cb_mass(i,7), &
                                    cb_mass(i,8), &
                                    cb_mass(i,9), &
                                    bounds(i,1,1), &
                                    bounds(i,2,1), &
                                    bounds(i,3,1), &
                                    bounds(i,4,1), &
                                    bounds(i,1,2), &
                                    bounds(i,2,2), &
                                    bounds(i,3,2), &
                                    bounds(i,4,2), &
                                    bounds(i,1,3), &
                                    bounds(i,2,3), &
                                    bounds(i,3,3), &
                                    bounds(i,4,3), &
                                    bounds(i,1,4), &
                                    bounds(i,2,4), &
                                    bounds(i,3,4), &
                                    bounds(i,4,4), &
                                    bounds(i,1,5), &
                                    bounds(i,2,5), &
                                    bounds(i,3,5), &
                                    bounds(i,4,5), &
                                    bounds(i,1,6), &
                                    bounds(i,2,6), &
                                    bounds(i,3,6), &
                                    bounds(i,4,6), &
                                    cntrline(i,1), &
                                    cntrline(i,2), &
                                    cntrline(i,3), &
                                    cntrline(i,4)
                            ELSE
                                WRITE(i+20, '(6X,F12.6,F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8,F24.8,F24.8,F24.8,' // &
                                                'F24.8)') &
                                    nondim_time, &
                                    cb_mass(i,1), &
                                    cb_mass(i,2), &
                                    cb_mass(i,3), &
                                    cb_mass(i,4), &
                                    cb_mass(i,5), &
                                    cb_mass(i,6), &
                                    cb_mass(i,7), &
                                    cb_mass(i,8), &
                                    cb_mass(i,9), &
                                    cb_mass(i,10), &
                                    bounds(i,1,1), &
                                    bounds(i,2,1), &
                                    bounds(i,3,1), &
                                    bounds(i,4,1), &
                                    bounds(i,5,1), &
                                    bounds(i,1,2), &
                                    bounds(i,2,2), &
                                    bounds(i,3,2), &
                                    bounds(i,4,2), &
                                    bounds(i,5,2), &
                                    bounds(i,1,3), &
                                    bounds(i,2,3), &
                                    bounds(i,3,3), &
                                    bounds(i,4,3), &
                                    bounds(i,5,3), &
                                    bounds(i,1,4), &
                                    bounds(i,2,4), &
                                    bounds(i,3,4), &
                                    bounds(i,4,4), &
                                    bounds(i,5,4), &
                                    bounds(i,1,5), &
                                    bounds(i,2,5), &
                                    bounds(i,3,5), &
                                    bounds(i,4,5), &
                                    bounds(i,5,5), &
                                    bounds(i,1,6), &
                                    bounds(i,2,6), &
                                    bounds(i,3,6), &
                                    bounds(i,4,6), &
                                    bounds(i,5,6), &
                                    cntrline(i,1), &
                                    cntrline(i,2), &
                                    cntrline(i,3), &
                                    cntrline(i,4), &
                                    cntrline(i,5)
                            END IF
                        END IF
                    END IF
                END IF
            END DO

        END SUBROUTINE s_write_cb_files ! -------------------------------------



        !>  This writes a formatted data file for the flow probe information
        !!  @param t_step Current time-step
        !!  @param q_cons_vf Conservative variables
        !!  @param accel_mag Acceleration magnitude information
        SUBROUTINE s_write_probe_files(t_step,q_cons_vf,accel_mag) ! -----------

            INTEGER, INTENT(IN) :: t_step
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q_cons_vf
            REAL(KIND(0d0)), DIMENSION(0:m,0:n,0:p), INTENT(IN) :: accel_mag

            REAL(KIND(0d0)), DIMENSION(-1:m) :: distx
            REAL(KIND(0d0)), DIMENSION(-1:n) :: disty
            REAL(KIND(0d0)), DIMENSION(-1:p) :: distz

            ! The cell-averaged partial densities, density, velocity, pressure,
            ! volume fractions, specific heat ratio function, liquid stiffness
            ! function, and sound speed.
            REAL(KIND(0d0))                                   :: lit_gamma, nbub
            REAL(KIND(0d0))                                   :: rho
            REAL(KIND(0d0)), DIMENSION(num_dims)              :: vel
            REAL(KIND(0d0))                                   :: pres
            REAL(KIND(0d0))                                   :: ptilde
            REAL(KIND(0d0))                                   :: ptot
            REAL(KIND(0d0))                                   :: alf
            REAL(KIND(0d0))                                   :: alfgr
            REAL(KIND(0d0)), DIMENSION(num_fluids)            :: alpha
            REAL(KIND(0d0))                                   :: gamma
            REAL(KIND(0d0))                                   :: pi_inf
            REAL(KIND(0d0))                                   :: c
            REAL(KIND(0d0))                                   :: M00, M10, M01, M20, M11, M02
            REAL(KIND(0d0))                                   :: varR, varV
            REAL(KIND(0d0)), DIMENSION(Nb)                    :: nR, R, nRdot, Rdot
            REAL(KIND(0d0))                                   :: accel
            REAL(KIND(0d0))                                   :: int_pres
            REAL(KIND(0d0))                                   :: max_pres
            REAL(KIND(0d0)), DIMENSION(2)             :: Re
            REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:)      :: We
            
            INTEGER :: i,j,k,l,s !< Generic loop iterator

            REAL(KIND(0d0)) :: nondim_time !< Non-dimensional time

            REAL(KIND(0d0)) :: tmp !<
            !! Temporary variable to store quantity for mpi_allreduce

            REAL(KIND(0d0)) :: blkmod1, blkmod2 !<
            !! Fluid bulk modulus for Woods mixture sound speed

            INTEGER :: npts !< Number of included integral points
            REAL(KIND(0d0)) :: rad, thickness !< For integral quantities 
            LOGICAL :: trigger !< For integral quantities

            ! Non-dimensional time calculation
            IF (time_stepper == 23) THEN
                nondim_time = mytime
            ELSE
                IF (t_step_old /= dflt_int) THEN
                    nondim_time = REAL(t_step + t_step_old,KIND(0d0))*dt
                ELSE
                    nondim_time = REAL(t_step,KIND(0d0))*dt !*1.d-5/10.0761131451d0
                END IF
            END IF

            DO i = 1, num_probes
                ! Zeroing out flow variables for all processors
                rho = 0d0
                DO s = 1, num_dims
                    vel(s) = 0d0
                END DO
                pres = 0d0
                gamma = 0d0
                pi_inf = 0d0
                c = 0d0
                accel = 0d0
                nR = 0d0; R = 0d0
                nRdot = 0d0; Rdot = 0d0
                nbub = 0d0
                M00 = 0d0
                M10 = 0d0
                M01 = 0d0
                M20 = 0d0
                M11 = 0d0
                M02 = 0d0
                varR = 0d0; varV = 0d0
                alf = 0d0

                ! Find probe location in terms of indices on a
                ! specific processor
                IF (n == 0) THEN ! 1D simulation
                    IF ( (probe(i)%x >= x_cb(-1)) .AND. (probe(i)%x <= x_cb(m)) ) THEN
                        DO s = -1, m
                            distx(s) = x_cb(s) - probe(i)%x
                            IF (distx(s) < 0d0) distx(s) = 1000d0
                        END DO
                        j = MINLOC(distx,1)
                        IF (j == 1) j = 2 ! Pick first point if probe is at edge
                        k = 0
                        l = 0

                        ! Computing/Sharing necessary state variables
                        CALL s_convert_to_mixture_variables( q_cons_vf, rho, &
                                             gamma, pi_inf, &
                                             Re, We, j-2,k,l)
                        DO s = 1, num_dims
                            vel(s) = q_cons_vf(cont_idx%end+s)%sf(j-2,k,l)/rho
                        END DO

                        IF (model_eqns == 4) THEN
                            lit_gamma = 1d0/fluid_pp(1)%gamma + 1d0
                            
                            !Tait pressure from density
                            pres = (pref + pi_inf) * (                     &
                                ( q_cons_vf(1)%sf(j-2,k,l)/               &
                                (rhoref*(1.d0-q_cons_vf(4)%sf(j-2,k,l)))  & 
                                ) ** lit_gamma )                        &
                                - pi_inf
                        ELSE IF (model_eqns == 2 .AND. (bubbles .NEQV. .TRUE.)) THEN
                            !Stiffened gas pressure from energy
                            pres = (                                       & 
                                q_cons_vf(E_idx)%sf(j-2,k,l)  -            &
                                0.5d0*(q_cons_vf(2)%sf(j-2,k,l)**2.d0)/q_cons_vf(1)%sf(j-2,k,l) - &
                                pi_inf &
                                ) / gamma
                        ELSE
                            !Stiffened gas pressure from energy with bubbles
                            pres = (                                       & 
                                (q_cons_vf(E_idx)%sf(j-2,k,l)  -            &
                                0.5d0*(q_cons_vf(mom_idx%beg)%sf(j-2,k,l)**2.d0)/rho) / &
                                (1.d0 - q_cons_vf(alf_idx)%sf(j-2,k,l)) - &
                                pi_inf &
                                ) / gamma
                        end IF

                        IF (bubbles) THEN
                            alf = q_cons_vf(alf_idx)%sf(j-2,k,l)
                            IF (num_fluids == 3) THEN
                                alfgr = q_cons_vf(alf_idx-1)%sf(j-2,k,l)
                            END IF
                            DO s = 1,nb
                                nR(s)   = q_cons_vf(bub_idx%rs(s))%sf(j-2,k,l)
                                nRdot(s)= q_cons_vf(bub_idx%vs(s))%sf(j-2,k,l)
                            END DO
                            CALL s_comp_n_from_cons( alf, nR, nbub)
                            IF (DEBUG) print*, 'In probe, nbub: ', nbub

                            IF (qbmm) THEN
                                M00 = q_cons_vf(bub_idx%moms(1,1))%sf(j-2,k,l)/nbub
                                M10 = q_cons_vf(bub_idx%moms(1,2))%sf(j-2,k,l)/nbub
                                M01 = q_cons_vf(bub_idx%moms(1,3))%sf(j-2,k,l)/nbub
                                M20 = q_cons_vf(bub_idx%moms(1,4))%sf(j-2,k,l)/nbub
                                M11 = q_cons_vf(bub_idx%moms(1,5))%sf(j-2,k,l)/nbub
                                M02 = q_cons_vf(bub_idx%moms(1,6))%sf(j-2,k,l)/nbub

                                M10 = M10/M00
                                M01 = M01/M00
                                M20 = M20/M00
                                M11 = M11/M00
                                M02 = M02/M00

                                varR = M20 - M10**2d0
                                varV = M02 - M01**2d0
                            END IF
                            R(:) = nR(:)/nbub                        
                            Rdot(:) = nRdot(:)/nbub                        
                        
                            ptilde = ptil(j-2,k,l)
                            ptot = pres - ptilde
                        END IF

                        ! Compute mixture sound speed
                        IF (alt_soundspeed .OR. regularization) THEN
                            DO s = 1, num_fluids
                                alpha(s) = q_cons_vf(E_idx+s)%sf(j-2,k,l)
                            END DO
                            blkmod1 = ((fluid_pp(1)%gamma +1d0)*pres + &
                                fluid_pp(1)%pi_inf)/fluid_pp(1)%gamma
                            blkmod2 = ((fluid_pp(2)%gamma +1d0)*pres + &
                                fluid_pp(2)%pi_inf)/fluid_pp(2)%gamma
                            c = (1d0/(rho*(alpha(1)/blkmod1 + alpha(2)/blkmod2)))
                        ELSE
                            c = (((gamma + 1d0)*pres + pi_inf)/(gamma*rho))
                        END IF

                        IF (mixture_err .AND. c < 0d0) THEN
                            c = sgm_eps
                        ELSE
                            c = SQRT(c)
                        END IF

                        accel = accel_mag(j-2,k,l)
                    END IF
                ELSEIF (p == 0) THEN ! 2D simulation
                    IF ( (probe(i)%x >= x_cb(-1)) .AND. (probe(i)%x <= x_cb(m)) ) THEN
                        IF ( (probe(i)%y >= y_cb(-1)) .AND. (probe(i)%y <= y_cb(n)) ) THEN
                            DO s = -1, m
                                distx(s) = x_cb(s) - probe(i)%x
                                IF (distx(s) < 0d0) distx(s) = 1000d0
                            END DO
                            DO s = -1, n
                                disty(s) = y_cb(s) - probe(i)%y
                                IF (disty(s) < 0d0) disty(s) = 1000d0
                            END DO
                            j = MINLOC(distx,1)
                            k = MINLOC(disty,1)
                            IF (j == 1) j = 2 ! Pick first point if probe is at edge
                            IF (k == 1) k = 2 ! Pick first point if probe is at edge
                            l = 0

                            ! Computing/Sharing necessary state variables
                            CALL s_convert_to_mixture_variables( q_cons_vf, rho, &
                                                 gamma, pi_inf, &
                                                 Re, We, j-2,k-2,l)
                            DO s = 1, num_dims
                                vel(s) = q_cons_vf(cont_idx%end+s)%sf(j-2,k-2,l)/rho
                            END DO

                            IF (model_eqns == 4) THEN
                                lit_gamma = 1d0/fluid_pp(1)%gamma + 1d0
                            
                                !Tait pressure from density
                                pres = (pref + pi_inf) * (                     &
                                    ( q_cons_vf(1)%sf(j-2,k-2,l)/               &
                                    (rhoref*(1.d0-q_cons_vf(4)%sf(j-2,k-2,l)))  & 
                                    ) ** lit_gamma )                        &
                                    - pi_inf
                            ELSE IF (model_eqns == 2 .AND. (bubbles .NEQV. .TRUE.)) THEN
                                !Stiffened gas pressure from energy
                                pres = (                                       & 
                                    q_cons_vf(E_idx)%sf(j-2,k-2,l)  -            &
                                    0.5d0*( (q_cons_vf(2)%sf(j-2,k-2,l)**2.d0 + &
                                    q_cons_vf(3)%sf(j-2,k-2,l)**2.d0)/q_cons_vf(1)%sf(j-2,k-2,l)) - &
                                    pi_inf &
                                    ) / gamma
                            ELSE
                                !Stiffened gas pressure from energy with bubbles
                                pres = (                                       & 
                                    (q_cons_vf(E_idx)%sf(j-2,k-2,l)  -            &
                                    0.5d0*(q_cons_vf(2)%sf(j-2,k-2,l)**2.d0 + q_cons_vf(3)%sf(j-2,k-2,l)**2.d0 ) &
                                    / q_cons_vf(1)%sf(j-2,k-2,l)) / &
                                    (1.d0 - q_cons_vf(alf_idx)%sf(j-2,k-2,l)) - &
                                    pi_inf &
                                    ) / gamma
                            end IF

                            IF (bubbles) THEN
                                alf = q_cons_vf(alf_idx)%sf(j-2,k-2,l)
                                DO s = 1,nb
                                    nR(s)   = q_cons_vf(bub_idx%rs(s))%sf(j-2,k-2,l)
                                    nRdot(s)= q_cons_vf(bub_idx%vs(s))%sf(j-2,k-2,l)
                                END DO
                                CALL s_comp_n_from_cons( alf, nR, nbub)
                                
                                R(:) = nR(:)/nbub                        
                                Rdot(:) = nRdot(:)/nbub                        
                            end IF

                            ! Compute mixture sound speed
                            IF (alt_soundspeed .OR. regularization) THEN
                                DO s = 1, num_fluids
                                    alpha(s) = q_cons_vf(E_idx+s)%sf(j-2,k-2,l)
                                END DO
                                blkmod1 = ((fluid_pp(1)%gamma +1d0)*pres + &
                                    fluid_pp(1)%pi_inf)/fluid_pp(1)%gamma
                                blkmod2 = ((fluid_pp(2)%gamma +1d0)*pres + &
                                    fluid_pp(2)%pi_inf)/fluid_pp(2)%gamma
                                c = (1d0/(rho*(alpha(1)/blkmod1 + alpha(2)/blkmod2)))
                            ELSE
                                c = (((gamma + 1d0)*pres + pi_inf)/(gamma*rho))
                            END IF

                            IF (mixture_err .AND. c < 0d0) THEN
                                c = sgm_eps
                            ELSE
                                c = SQRT(c)
                            END IF

                            accel = accel_mag(j-2,k-2,l)
                        END IF
                    END IF
                ELSE ! 3D simulation
                    IF ( (probe(i)%x >= x_cb(-1)) .AND. (probe(i)%x <= x_cb(m)) ) THEN
                        IF ( (probe(i)%y >= y_cb(-1)) .AND. (probe(i)%y <= y_cb(n)) ) THEN
                            IF ( (probe(i)%z >= z_cb(-1)) .AND. (probe(i)%z <= z_cb(p)) ) THEN
                                DO s = -1, m
                                    distx(s) = x_cb(s) - probe(i)%x
                                    IF (distx(s) < 0d0) distx(s) = 1000d0
                                END DO
                                DO s = -1, n
                                    disty(s) = y_cb(s) - probe(i)%y
                                    IF (disty(s) < 0d0) disty(s) = 1000d0
                                END DO
                                DO s = -1, p
                                    distz(s) = z_cb(s) - probe(i)%z
                                    IF (distz(s) < 0d0) distz(s) = 1000d0
                                END DO
                                j = MINLOC(distx,1)
                                k = MINLOC(disty,1)
                                l = MINLOC(distz,1)
                                IF (j == 1) j = 2 ! Pick first point if probe is at edge
                                IF (k == 1) k = 2 ! Pick first point if probe is at edge
                                IF (l == 1) l = 2 ! Pick first point if probe is at edge
        
                                ! Computing/Sharing necessary state variables
                                CALL s_convert_to_mixture_variables( q_cons_vf, rho, &
                                                     gamma, pi_inf, &
                                                     Re, We, j-2,k-2,l-2)
                                DO s = 1, num_dims
                                    vel(s) = q_cons_vf(cont_idx%end+s)%sf(j-2,k-2,l-2)/rho
                                END DO

                                pres = (q_cons_vf(E_idx)%sf(j-2,k-2,l-2) - 0.5d0*rho*DOT_PRODUCT(vel,vel)-pi_inf)/gamma

                                ! Compute mixture sound speed
                                IF (alt_soundspeed .OR. regularization) THEN
                                    DO s = 1, num_fluids
                                        alpha(s) = q_cons_vf(E_idx+s)%sf(j-2,k-2,l-2)
                                    END DO
                                    blkmod1 = ((fluid_pp(1)%gamma +1d0)*pres + &
                                        fluid_pp(1)%pi_inf)/fluid_pp(1)%gamma
                                    blkmod2 = ((fluid_pp(2)%gamma +1d0)*pres + &
                                        fluid_pp(2)%pi_inf)/fluid_pp(2)%gamma
                                    c = (1d0/(rho*(alpha(1)/blkmod1 + alpha(2)/blkmod2)))
                                ELSE
                                    c = (((gamma + 1d0)*pres + pi_inf)/(gamma*rho))
                                END IF

                                IF (mixture_err .AND. c < 0d0) THEN
                                    c = sgm_eps
                                ELSE
                                    c = SQRT(c)
                                END IF

                                accel = accel_mag(j-2,k-2,l-2)
                            END IF
                        END IF
                    END IF
                END IF

                IF (num_procs > 1) THEN
                    tmp = rho
                    CALL s_mpi_allreduce_sum(tmp,rho)
                    DO s = 1, num_dims
                        tmp = vel(s)
                        CALL s_mpi_allreduce_sum(tmp,vel(s))
                    END DO
                    tmp = pres
                    CALL s_mpi_allreduce_sum(tmp,pres)
                    tmp = gamma
                    CALL s_mpi_allreduce_sum(tmp,gamma)
                    tmp = pi_inf
                    CALL s_mpi_allreduce_sum(tmp,pi_inf)
                    tmp = c
                    CALL s_mpi_allreduce_sum(tmp,c)
                    tmp = accel
                    CALL s_mpi_allreduce_sum(tmp,accel)

                    IF (bubbles) THEN
                        tmp = alf
                        CALL s_mpi_allreduce_sum(tmp,alf)
                        tmp = alfgr
                        CALL s_mpi_allreduce_sum(tmp,alfgr)
                        tmp = nbub
                        CALL s_mpi_allreduce_sum(tmp,nbub)
                        tmp = nR(1)
                        CALL s_mpi_allreduce_sum(tmp,nR(1))
                        tmp = nRdot(1)
                        CALL s_mpi_allreduce_sum(tmp,nRdot(1))
                        tmp = M00
                        CALL s_mpi_allreduce_sum(tmp,M00)
                        tmp = R(1)
                        CALL s_mpi_allreduce_sum(tmp,R(1))
                        tmp = Rdot(1)
                        CALL s_mpi_allreduce_sum(tmp,Rdot(1))
                        tmp = ptilde
                        CALL s_mpi_allreduce_sum(tmp,ptilde)
                        tmp = ptot
                        CALL s_mpi_allreduce_sum(tmp,ptot)

                        IF (qbmm) THEN
                            tmp = varR
                            CALL s_mpi_allreduce_sum(tmp,varR)
                            tmp = varV
                            CALL s_mpi_allreduce_sum(tmp,varV)

                            tmp = M10
                            CALL s_mpi_allreduce_sum(tmp,M10)
                            tmp = M01
                            CALL s_mpi_allreduce_sum(tmp,M01)
                            tmp = M20
                            CALL s_mpi_allreduce_sum(tmp,M20)
                            tmp = M02
                            CALL s_mpi_allreduce_sum(tmp,M02)

                        END IF
                    END IF
                END IF

                IF (proc_rank == 0) THEN
                    IF (n == 0) THEN
                        IF (bubbles .AND. (num_fluids <= 2)) THEN
                            IF (qbmm) THEN
                                WRITE(i+30,'(6x,f12.6,14f28.16)') &
                                    nondim_time, &
                                    rho, &
                                    vel(1), &
                                    pres, &
                                    alf, &
                                    R(1), &
                                    Rdot(1), &
                                    nR(1), &
                                    nRdot(1), &
                                    varR, &
                                    varV, &
                                    M10, &
                                    M01, &
                                    M20, &
                                    M02
                            ELSE
                                WRITE(i+30,'(6x,f12.6,8f24.8)') &
                                    nondim_time, &
                                    rho, &
                                    vel(1), &
                                    pres, &
                                    alf, &
                                    R(1), &
                                    Rdot(1), &
                                    nR(1), &
                                    nRdot(1)
                                    ! ptilde, &
                                    ! ptot
                            END IF
                        ELSE IF (bubbles .AND. (num_fluids ==3)) THEN
                            WRITE(i+30,'(6x,f12.6,f24.8,f24.8,f24.8,f24.8,f24.8,' // &
                                           'f24.8,f24.8,f24.8,f24.8,f24.8, f24.8)') &
                                nondim_time, &
                                rho, &
                                vel(1), &
                                pres, &
                                alf, &
                                alfgr, &
                                nR(1), &
                                nRdot(1), &
                                R(1), &
                                Rdot(1), &
                                ptilde, &
                                ptot
                        ELSE IF (bubbles .AND. num_fluids ==4) THEN
                             WRITE(i+30,'(6x,f12.6,f24.8,f24.8,f24.8,f24.8,' // &
                                           'f24.8,f24.8,f24.8,f24.8,f24.8,f24.8,f24.8,f24.8,f24.8)') &
                                nondim_time, &
                                q_cons_vf(1)%sf(j-2,0,0), &
                                q_cons_vf(2)%sf(j-2,0,0), &
                                q_cons_vf(3)%sf(j-2,0,0), &
                                q_cons_vf(4)%sf(j-2,0,0), &
                                q_cons_vf(5)%sf(j-2,0,0), &
                                q_cons_vf(6)%sf(j-2,0,0), &
                                q_cons_vf(7)%sf(j-2,0,0), &
                                q_cons_vf(8)%sf(j-2,0,0), &
                                q_cons_vf(9)%sf(j-2,0,0), &
                                q_cons_vf(10)%sf(j-2,0,0), &
                                nbub,&
                                R(1),&
                                Rdot(1)
                        ELSE
                            WRITE(i+30,'(6X,F12.6,F24.8,F24.8,F24.8)') &
                                nondim_time, &
                                rho, &
                                vel(1), &
                                pres
                        END IF
                    ELSEIF (p == 0) THEN
                        IF (bubbles) THEN
                            WRITE(i+30,'(6X,10F24.8)') &
                                nondim_time, &
                                rho, &
                                vel(1), &
                                vel(2), &
                                pres, &
                                alf, &
                                nR(1), &
                                nRdot(1), &
                                R(1), &
                                Rdot(1)
                        ELSE
                            WRITE(i+30,'(6X,F12.6,F24.8,F24.8,F24.8)') &
                                nondim_time, &
                                rho, &
                                vel(1), &
                                pres
                        END IF
                    ELSE
                        WRITE(i+30,'(6X,F12.6,F24.8,F24.8,F24.8,F24.8,' // &
                                           'F24.8,F24.8,F24.8,F24.8,' // &
                                           'F24.8)') &
                            nondim_time, &
                            rho, &
                            vel(1), &
                            vel(2), &
                            vel(3), &
                            pres, &
                            gamma, &
                            pi_inf, &
                            c, &
                            accel
                    END IF
                END IF
            END DO

        IF (integral_wrt .AND. bubbles) THEN
            IF (n == 0) THEN ! 1D simulation
                DO i = 1, num_integrals
                    int_pres = 0d0
                    max_pres = 0d0
                    k = 0; l = 0
                    npts = 0
                    DO j = 1,m
                        pres = 0d0
                        DO s = 1, num_dims
                            vel(s) = 0d0
                        END DO
                        rho = 0d0
                        pres = 0d0
                        gamma = 0d0
                        pi_inf = 0d0

                        IF ( (integral(i)%xmin <= x_cb(j)) .AND. (integral(i)%xmax >= x_cb(j)) ) THEN
                            npts = npts + 1
                            CALL s_convert_to_mixture_variables( q_cons_vf, rho, &
                                             gamma, pi_inf, &
                                             Re, We, j,k,l)
                            DO s = 1, num_dims
                                vel(s) = q_cons_vf(cont_idx%end+s)%sf(j,k,l)/rho
                            END DO

                            pres = (                                       & 
                                (q_cons_vf(E_idx)%sf(j,k,l)  -            &
                                0.5d0*(q_cons_vf(mom_idx%beg)%sf(j,k,l)**2.d0)/rho) / &
                                (1.d0 - q_cons_vf(alf_idx)%sf(j,k,l)) - &
                                pi_inf &
                                ) / gamma
                            int_pres = int_pres + (pres-1.d0)**2.d0
                        END IF
                    END DO
                    int_pres = dsqrt(int_pres / (1.d0*npts))
    
                    IF (num_procs > 1) THEN
                        tmp = int_pres
                        CALL s_mpi_allreduce_sum(tmp,int_pres)
                    END IF

                    IF (proc_rank == 0) THEN
                        IF (bubbles .AND. (num_fluids <= 2)) THEN
                            WRITE(i+70,'(6x,f12.6,f24.8)') &
                                nondim_time, int_pres
                        end IF
                    END IF
                END DO
            ELSEIF (p == 0) THEN
                IF (num_integrals .NE. 3) THEN
                    PRINT '(A)', 'Incorrect number of integrals'
                    CALL s_mpi_abort()
                END IF

                rad = integral(1)%xmax
                thickness = integral(1)%xmin

                DO i = 1,num_integrals
                    int_pres = 0d0
                    max_pres = 0d0 
                    l = 0
                    npts = 0
                    DO j = 1,m
                        DO k = 1,n
                            trigger = .FALSE.
                            IF (i==1) THEN
                                !inner portion
                                IF ( dsqrt(x_cb(j)**2.d0 + y_cb(k)**2.d0) < (rad - 0.5d0*thickness) ) &
                                    trigger = .TRUE.
                            ELSEIF (i==2) THEN
                                !net region
                                IF ( dsqrt(x_cb(j)**2.d0 + y_cb(k)**2.d0) > (rad - 0.5d0*thickness) .AND. &
                                     dsqrt(x_cb(j)**2.d0 + y_cb(k)**2.d0) < (rad + 0.5d0*thickness) ) &
                                    trigger = .TRUE.
                            ELSEIF (i==3) THEN
                                !everything else
                                IF ( dsqrt(x_cb(j)**2.d0 + y_cb(k)**2.d0) > (rad + 0.5d0*thickness) ) &
                                    trigger = .TRUE.
                            END IF

                            pres = 0d0
                            DO s = 1, num_dims
                                vel(s) = 0d0
                            END DO
                            rho = 0d0
                            pres = 0d0
                            gamma = 0d0
                            pi_inf = 0d0

                            IF (trigger) THEN
                                npts = npts + 1
                                CALL s_convert_to_mixture_variables( q_cons_vf, rho, &
                                    gamma, pi_inf, &
                                    Re, We, j,k,l)
                                DO s = 1, num_dims
                                    vel(s) = q_cons_vf(cont_idx%end+s)%sf(j,k,l)/rho
                                END DO

                                pres = (                                       & 
                                    (q_cons_vf(E_idx)%sf(j,k,l)  -            &
                                    0.5d0*(q_cons_vf(mom_idx%beg)%sf(j,k,l)**2.d0)/rho) / &
                                    (1.d0 - q_cons_vf(alf_idx)%sf(j,k,l)) - &
                                    pi_inf &
                                    ) / gamma
                                int_pres = int_pres + abs(pres-1.d0)
                                max_pres = max(max_pres,abs(pres-1.d0))
                            END IF

                       END DO
                    END DO

                    IF (npts > 0) THEN
                        int_pres = int_pres/(1.d0*npts)
                    ELSE
                        int_pres = 0.d0
                    end IF
                    
                    IF (num_procs > 1) THEN
                        tmp = int_pres
                        CALL s_mpi_allreduce_sum(tmp,int_pres)

                        tmp = max_pres
                        CALL s_mpi_allreduce_max(tmp,max_pres)
                    END IF

                    IF (proc_rank == 0) THEN
                        IF (bubbles .AND. (num_fluids <= 2)) THEN
                            WRITE(i+70,'(6x,f12.6,f24.8,f24.8)') &
                                nondim_time, int_pres, max_pres
                        end IF
                    END IF
                END DO
            END IF
        END IF


        END SUBROUTINE s_write_probe_files ! -----------------------------------



        !>  The goal of this subroutine is to write to the run-time
        !!      information file basic footer information applicable to
        !!      the current computation and to close the file when done.
        !!      The footer contains the stability criteria extrema over
        !!      all of the time-steps and the simulation run-time.
        SUBROUTINE s_close_run_time_information_file() ! -----------------------
           
            
            REAL(KIND(0d0)) :: run_time !< Run-time of the simulation
            
            
            ! Writing the footer of and closing the run-time information file
            WRITE(1,'(A)') '----------------------------------------' // &
                           '----------------------------------------'
            WRITE(1,'(A)') ''
            
                                 WRITE(1, '(A,F9.6)') 'ICFL Max: ', icfl_max
            IF(ANY(Re_size > 0)) WRITE(1, '(A,F9.6)') 'VCFL Max: ', vcfl_max
            IF(    We_size > 0 ) WRITE(1, '(A,F9.6)') 'CCFL Max: ', ccfl_max
            IF(ANY(Re_size > 0)) WRITE(1,'(A,F10.6)')   'Rc Min: ',   Rc_min
            
            CALL CPU_TIME(run_time)
            
            WRITE(1,'(A)') ''
            WRITE(1,'(A,I0,A)') 'Run-time: ', INT(ANINT(run_time)), 's'
            WRITE(1,'(A)') '========================================' // &
                           '========================================'
            CLOSE(1)
            
            
        END SUBROUTINE s_close_run_time_information_file ! ---------------------
        
        
        
        
        !> Closes communication files 
        SUBROUTINE s_close_com_files() ! ---------------------------------------

            INTEGER :: i !< Generic loop iterator

            DO i = 1, num_fluids
                IF (com_wrt(i)) THEN
                    CLOSE(i+10)
                END IF
            END DO

        END SUBROUTINE s_close_com_files ! -------------------------------------




        !> Closes coherent body files
        SUBROUTINE s_close_cb_files() ! ---------------------------------------

            INTEGER :: i !< Generic loop iterator

            DO i = 1, num_fluids
                IF (cb_wrt(i)) THEN
                    CLOSE(i+20)
                END IF
            END DO

        END SUBROUTINE s_close_cb_files ! -------------------------------------




        !> Closes probe files
        SUBROUTINE s_close_probe_files() ! -------------------------------------


            INTEGER :: i !< Generic loop iterator

            DO i = 1, num_probes
                CLOSE(i+30)
            END DO

        END SUBROUTINE s_close_probe_files ! -----------------------------------




        !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
        SUBROUTINE s_initialize_data_output_module() ! -------------------------
           
            TYPE(bounds_info) :: ix, iy, iz

            INTEGER :: i !< Generic loop iterator
            
            ! Allocating/initializing ICFL, VCFL, CCFL and Rc stability criteria
                                 ALLOCATE(icfl_sf(0:m,0:n,0:p)); icfl_max = 0d0
            IF(ANY(Re_size > 0)) ALLOCATE(vcfl_sf(0:m,0:n,0:p)); vcfl_max = 0d0
            IF(    We_size > 0 ) ALLOCATE(ccfl_sf(0:m,0:n,0:p)); ccfl_max = 0d0
            IF(ANY(Re_size > 0)) ALLOCATE(  Rc_sf(0:m,0:n,0:p));   Rc_min = 1d3
            
            
            ! Associating the procedural pointer to the appropriate subroutine
            ! that will be utilized in the conversion to the mixture variables
           
            IF (model_eqns == 1) THEN        ! Gamma/pi_inf model
                s_convert_to_mixture_variables => &
                             s_convert_mixture_to_mixture_variables
            ELSEIF (bubbles) THEN           ! Volume fraction for bubbles
                s_convert_to_mixture_variables => &
                             s_convert_species_to_mixture_variables_bubbles
            ELSE                            ! Volume fraction model
                s_convert_to_mixture_variables => &
                             s_convert_species_to_mixture_variables
            END IF

            ! Allocating the generic storage for the flow variable(s) that are
            ! going to be written to the CoM data files
            IF (ANY(com_wrt)) THEN
                ! num_fluids, mass, x-loc, y-loc, z-loc, x-vel, y-vel, z-vel, x-acc, y-acc, z-acc
                ALLOCATE(q_com(num_fluids,10))
                ! num_fluids, 2 lateral directions, 5 higher moment orders
                ALLOCATE(moments(num_fluids,2,5))
            END IF

            IF (ANY(cb_wrt)) THEN
                ! num_fluids, mass for 5 threshold mass fractions, area for 5 threshold mass fractions
                ALLOCATE(cb_mass(num_fluids,10))
                ! num_fluids, 5 threshold mass fractions, xmin, xmax, ymin, ymax, zmin, zmax
                ALLOCATE(bounds(num_fluids,5,6))
                ! num_fluids, centerline dimension for 5 threshold mass fractions
                ALLOCATE(cntrline(num_fluids,5))
            END IF

            IF (probe_wrt) THEN
                ALLOCATE(accel_mag(0:m,0:n,0:p))
                ALLOCATE(x_accel(0:m,0:n,0:p))
                IF (n > 0) THEN
                    ALLOCATE(y_accel(0:m,0:n,0:p))
                    IF (p > 0) THEN
                        ALLOCATE(z_accel(0:m,0:n,0:p))
                    END IF
                END IF
            END IF

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
                ALLOCATE(energy(ix%beg:ix%end,iy%beg:iy%end,iz%beg:iz%end))

            END IF

            IF (parallel_io .NEQV. .TRUE.) THEN
                s_write_data_files => s_write_serial_data_files
            ELSE
                s_write_data_files => s_write_parallel_data_files
            END IF
            
        END SUBROUTINE s_initialize_data_output_module ! -----------------------
        
        
        
        
        
        !> Module deallocation and/or disassociation procedures
        SUBROUTINE s_finalize_data_output_module() ! ---------------------------

            INTEGER :: i !< Generic loop iterator
            
            ! Deallocating the ICFL, VCFL, CCFL, and Rc stability criteria
                                 DEALLOCATE(icfl_sf)
            IF(ANY(Re_size > 0)) DEALLOCATE(vcfl_sf)
            IF(    We_size > 0 ) DEALLOCATE(ccfl_sf)
            IF(ANY(Re_size > 0)) DEALLOCATE(  Rc_sf)
            
            ! Deallocating the storage employed for the flow variables that
            ! were written to CoM and probe files
            IF (ANY(com_wrt)) THEN
                DEALLOCATE(q_com)
                DEALLOCATE(moments)
            END IF

            IF (ANY(cb_wrt)) THEN
                DEALLOCATE(cb_mass)
                DEALLOCATE(bounds)
                DEALLOCATE(cntrline)
            END IF
            IF (probe_wrt) THEN
                DEALLOCATE(accel_mag)
                DEALLOCATE(x_accel)
                IF (n > 0) THEN
                    DEALLOCATE(y_accel)
                    IF (p > 0) THEN
                        DEALLOCATE(z_accel)
                    END IF
                END IF
            END IF

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
                DEALLOCATE(energy)
            END IF

            ! Disassociating the pointer to the procedure that was utilized to
            ! to convert mixture or species variables to the mixture variables
            s_convert_to_mixture_variables => NULL()
            s_write_data_files => NULL()
            
        END SUBROUTINE s_finalize_data_output_module ! -------------------------
        
        
        
        
        
END MODULE m_data_output
