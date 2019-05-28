!>
!! @file p_main.f90
!! @brief  Quasi-conservative, shock- and interface- capturing finite-volume
!!              scheme for the multicomponent Navier-Stokes equations. The system
!!              is augmented with the relevant advection equations to capture the
!!              material interfaces and closed by the stiffened equation of state
!!              as well as any required mixture relations. The effects of surface
!!              tension are included and modeled through a volume force that acts
!!              across the diffuse material interface regions. The implementation
!!              specifics of surface tension may be found in the work by Perigaud
!!              and Saurel (2005). Note that both viscous and capillarity effects
!!              are only available in the volume fraction model.
!! @author spencer
!! @version 1.1
!! @date 1/1/1
PROGRAM p_main
    
    
    ! Dependencies =============================================================
    USE m_derived_types        !< Definitions of the derived types

    USE m_fftw                 !< Module for FFTW functions
    
    USE m_global_parameters    !< Definitions of the global parameters
    
    USE m_mpi_proxy            !< Message passing interface (MPI) module proxy
    
    USE m_start_up             !< Reading and checking procedures for the input
                               !< and the initial condition and grid data files
    
    USE m_variables_conversion !< State variables type conversion procedures
    
    USE m_weno                 !< Weighted and essentially non-oscillatory (WENO)
                               !! schemes for spatial reconstruction of variables
    
    USE m_riemann_solvers      !< Exact and approximate Riemann problem solvers
    
    USE m_cbc                  !< Characteristic boundary conditions (CBC)
    
    USE m_rhs                  !< Right-hand-side (RHS) evaluation procedures
    
    USE m_data_output          !< Run-time info & solution data output procedures

    USE m_derived_variables    !< Derived variables evaluation procedures
    
    USE m_time_steppers        !< Time-stepping algorithms

    USE m_anti_diffusion       !< Anti diffusion step of Shyue
    ! ==========================================================================
    
    
    IMPLICIT NONE
    
    INTEGER :: t_step !< Iterator for the time-stepping loop

    CALL system_clock(COUNT=cpu_start, COUNT_RATE=cpu_rate)
    
    ! Initializing MPI execution environment
    CALL s_mpi_initialize()
    
    
    ! The rank 0 processor assigns default values to the user inputs prior to
    ! reading them in from the input file. Next, the user inputs are read and
    ! their consistency is checked. The identification of any inconsistencies
    ! will result in the termination of the simulation.
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
    
    ! Computation of parameters, allocation of memory, association of pointers,
    ! and/or the execution of any other tasks needed to properly configure the
    ! modules. The preparations below DO NOT DEPEND on the grid being complete.
    CALL s_initialize_global_parameters_module()
    CALL s_initialize_mpi_proxy_module()
    IF (grid_geometry == 3) CALL s_initialize_fftw_module()
    CALL s_initialize_variables_conversion_module()
    CALL s_initialize_start_up_module()
    CALL s_initialize_riemann_solvers_module()
    CALL s_initialize_rhs_module()
    CALL s_initialize_data_output_module()
    CALL s_initialize_derived_variables_module()
    CALL s_initialize_time_steppers_module()    
    IF (anti_diffusion) CALL s_initialize_anti_diffusion_module()

    ! Associate pointers for serial or parallel I/O
    IF (parallel_io .NEQV. .TRUE.) THEN
        s_read_data_files => s_read_serial_data_files
        s_write_data_files => s_write_serial_data_files
    ELSE
        s_read_data_files => s_read_parallel_data_files
        s_write_data_files => s_write_parallel_data_files
    END IF

    ! Reading in the user provided initial condition and grid data
    CALL s_read_data_files(q_cons_ts(1)%vf)
    IF (model_eqns == 3) CALL s_initialize_internal_energy_equations(q_cons_ts(1)%vf)

    ! Populating the buffers of the grid variables using the boundary conditions
    CALL s_populate_grid_variables_buffers()
    
    CALL s_populate_variables_buffers(q_cons_ts(1)%vf)
    IF (We_size > 0 .AND. (We_riemann_flux .OR. We_rhs_flux)) CALL s_account_for_capillary_potential_energy(q_cons_ts(1)%vf)
   
    ! Computation of parameters, allocation of memory, association of pointers,
    ! and/or execution of any other tasks that are needed to properly configure
    ! the modules. The preparations below DO DEPEND on the grid being complete.
    CALL s_initialize_weno_module()
    CALL s_initialize_cbc_module()

    CALL s_initialize_derived_variables()

    ! Setting the time-step iterator to the first time-step
    t_step = t_step_start


    OPEN(unit=2,file='times.dat')
    CLOSE(2)
    
    ! Time-stepping Loop =======================================================
    DO
        IF (proc_rank ==0) PRINT*, 'Time step = ', t_step, ' of ', t_step_stop
        !IF(t_step == t_step_stop) stop 
        
        ! Total-variation-diminishing (TVD) Runge-Kutta (RK) time-steppers
        IF(time_stepper == 1) THEN
            CALL s_1st_order_tvd_rk(t_step)
        ELSEIF(time_stepper == 2) THEN
            CALL s_2nd_order_tvd_rk(t_step)
        ELSEIF(time_stepper == 3) THEN
            CALL s_3rd_order_tvd_rk(t_step)
        ELSEIF(time_stepper == 4) THEN
            CALL s_4th_order_rk(t_step)
        ELSE
            CALL s_5th_order_rk(t_step)
        END IF

        CALL s_compute_derived_variables(t_step)

        IF (anti_diffusion) CALL s_compute_anti_diffusion()

        ! Time-stepping loop controls
        IF(t_step == t_step_stop) THEN
            EXIT 
        ELSE
            t_step = t_step + 1
        END IF
       
        ! Backing up the grid and conservative variables data
        IF(MOD(t_step-t_step_start, t_step_save) == 0) THEN
            CALL s_write_data_files(q_cons_ts(1)%vf, t_step)
        END IF

        CALL system_clock(cpu_end)
        IF(proc_rank == 0) THEN
            OPEN(unit=2, file='times.dat', position = "append", status = "old")
                WRITE(2,*) t_step, REAL(cpu_end-cpu_start)/cpu_rate
            CLOSE(2)
        END IF
        
        IF (proc_rank == 0) THEN
            PRINT '(a)', ' ------------------------------------------------------ '
        END IF
    END DO
    ! ==========================================================================
    
    ! Disassociate pointers for serial and parallel I/O
    s_read_data_files => NULL()
    s_write_data_files => NULL()
    
    ! Deallocation and/or disassociation procedures for the modules
    IF (anti_diffusion) CALL s_finalize_anti_diffusion_module()
    CALL s_finalize_time_steppers_module()
    CALL s_finalize_derived_variables_module()
    CALL s_finalize_data_output_module()
    CALL s_finalize_rhs_module()
    CALL s_finalize_cbc_module()
    CALL s_finalize_riemann_solvers_module()
    CALL s_finalize_weno_module()
    CALL s_finalize_start_up_module()
    CALL s_finalize_variables_conversion_module()
    IF (grid_geometry == 3) CALL s_finalize_fftw_module()
    CALL s_finalize_mpi_proxy_module()
    CALL s_finalize_global_parameters_module()
    
    
    ! Terminating MPI execution environment
    CALL s_mpi_finalize()
    
    
END PROGRAM p_main
