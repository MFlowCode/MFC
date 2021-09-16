!>
!! @file m_mpi_proxy.f90
!! @brief Contains module m_mpi_proxy
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief The module serves as a proxy to the parameters and subroutines
!!          available in the MPI implementation's MPI module. Specifically,
!!          the purpose of the proxy is to harness basic MPI commands into
!!          more complicated procedures as to accomplish the communication
!!          goals for the simulation.
MODULE m_mpi_proxy
    
    
    ! Dependencies =============================================================
    USE mpi                    !< Message passing interface (MPI) module
    
    USE m_derived_types        !< Definitions of the derived types
    
    USE m_global_parameters    !< Definitions of the global parameters
    ! ==========================================================================
    
    
    IMPLICIT NONE
    
    
    REAL(KIND(0d0)), PRIVATE, ALLOCATABLE, DIMENSION(:) :: q_cons_buff_send !<
    !! This variable is utilized to pack and send the buffer of the cell-average
    !! conservative variables, for a single computational domain boundary at the
    !! time, to the relevant neighboring processor.
    
    REAL(KIND(0d0)), PRIVATE, ALLOCATABLE, DIMENSION(:) :: q_cons_buff_recv !<
    !! q_cons_buff_recv is utilized to receive and unpack the buffer of the cell-
    !! average conservative variables, for a single computational domain boundary
    !! at the time, from the relevant neighboring processor.
    
    !> @name Generic flags used to identify and report MPI errors
    !> @{
    INTEGER, PRIVATE :: err_code, ierr
    !> @}
    

    CONTAINS
        
        
        
        !> The subroutine intializes the MPI execution environment
        !!      and queries both the number of processors which will be
        !!      available for the job and the local processor rank.        
        SUBROUTINE s_mpi_initialize() ! ----------------------------------------

            ! Initializing the MPI environment
            CALL MPI_INIT(ierr)
            
            
            ! Checking whether the MPI environment has been properly intialized
            IF(ierr /= MPI_SUCCESS) THEN
                PRINT '(A)', 'Unable to initialize MPI environment. Exiting ...'
                CALL MPI_ABORT(MPI_COMM_WORLD, err_code, ierr)
            END IF
            
            
            ! Querying the number of processors available for the job
            CALL MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
            
            
            ! Querying the rank of the local processor
            CALL MPI_COMM_RANK(MPI_COMM_WORLD, proc_rank, ierr)
            
        END SUBROUTINE s_mpi_initialize ! --------------------------------------
        
        
        
        
        
        !> The subroutine terminates the MPI execution environment.
        SUBROUTINE s_mpi_abort() ! ---------------------------------------------

            ! Terminating the MPI environment
            CALL MPI_ABORT(MPI_COMM_WORLD, err_code, ierr)
            
        END SUBROUTINE s_mpi_abort ! -------------------------------------------
        
        
        
        
        !> The subroutine that initializes MPI data structures
        !!  @param q_cons_vf Conservative variables
        SUBROUTINE s_initialize_mpi_data(q_cons_vf) ! --------------------------

            TYPE(scalar_field), &
            DIMENSION(sys_size), &
            INTENT(IN) :: q_cons_vf

            INTEGER, DIMENSION(num_dims) :: sizes_glb, sizes_loc
            INTEGER :: ierr

            INTEGER :: i !< Generic loop iterator

            DO i = 1, sys_size
                MPI_IO_DATA%var(i)%sf => q_cons_vf(i)%sf(0:m,0:n,0:p)
            END DO

            ! Define global(g) and local(l) sizes for flow variables
            sizes_glb(1) = m_glb+1; sizes_loc(1) = m+1
            IF (n > 0) THEN
                sizes_glb(2) = n_glb+1; sizes_loc(2) = n+1
                IF (p > 0) THEN
                    sizes_glb(3) = p_glb+1; sizes_loc(3) = p+1
                END IF
            END IF

            ! Define the view for each variable
            DO i = 1, sys_size
                CALL MPI_TYPE_CREATE_SUBARRAY(num_dims,sizes_glb,sizes_loc,start_idx,&
                    MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,MPI_IO_DATA%view(i),ierr)
                CALL MPI_TYPE_COMMIT(MPI_IO_DATA%view(i),ierr)
            END DO

        END SUBROUTINE s_initialize_mpi_data ! ---------------------------------
        
        
        
        
        
        !> Halts all processes until all have reached barrier.
        SUBROUTINE s_mpi_barrier() ! -------------------------------------------

            ! Calling MPI_BARRIER
            CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

        END SUBROUTINE s_mpi_barrier ! -----------------------------------------



        !> The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
        SUBROUTINE s_initialize_mpi_proxy_module() ! ---------------------------
           
            ! Allocating q_cons_buff_send and q_cons_buff_recv. Please note that
            ! for the sake of simplicity, both variables are provided sufficient
            ! storage to hold the largest buffer in the computational domain.
            
            IF(n > 0) THEN
               
               IF(p > 0) THEN
                  ALLOCATE( q_cons_buff_send( 0 : -1+buff_size*sys_size * &
                                                         (m+2*buff_size+1) * &
                                                         (n+2*buff_size+1) * &
                                                         (p+2*buff_size+1) / &
                                                (MIN(m,n,p)+2*buff_size+1) ) )
               ELSE
                  ALLOCATE( q_cons_buff_send( 0 : -1+buff_size*sys_size * &
                                                  (MAX(m,n)+2*buff_size+1) ) )
               END IF
               
            ELSE
               
               ALLOCATE(q_cons_buff_send(0:-1+buff_size*sys_size))
               
            END IF
            
            ALLOCATE(q_cons_buff_recv(0:UBOUND(q_cons_buff_send,1)))
            
            
        END SUBROUTINE s_initialize_mpi_proxy_module ! -------------------------
        
        
        !>  Since only the processor with rank 0 reads and verifies
        !!      the consistency of user inputs, these are initially not
        !!      available to the other processors. Then, the purpose of
        !!      this subroutine is to distribute the user inputs to the
        !!      remaining processors in the communicator.       
        SUBROUTINE s_mpi_bcast_user_inputs() ! ---------------------------------

            INTEGER :: i,j !< Generic loop iterator
            
            
            ! Logistics
            CALL MPI_BCAST( case_dir     , LEN(case_dir), MPI_CHARACTER ,      &
                                                 0      , MPI_COMM_WORLD, ierr )
            CALL MPI_BCAST( run_time_info,       1      , MPI_LOGICAL   ,      &
                                                 0      , MPI_COMM_WORLD, ierr )
            CALL MPI_BCAST( t_step_old, 1   , MPI_INTEGER   ,   &
                            0   , MPI_COMM_WORLD, ierr )
            
            
            ! Computational domain parameters
            CALL MPI_BCAST(m, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(p, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(m_glb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(n_glb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(p_glb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            
            CALL MPI_BCAST(cyl_coord     , 1, MPI_LOGICAL         , &
                                           0, MPI_COMM_WORLD, ierr  )

            CALL MPI_BCAST(dt, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            
            CALL MPI_BCAST(t_step_start, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST(t_step_stop , 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST(t_step_save , 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)


            CALL MPI_BCAST( t_tol  ,        1      , &
                                MPI_DOUBLE_PRECISION,        0      , &
                                MPI_COMM_WORLD, ierr                  )
            CALL MPI_BCAST(debug, 1, MPI_LOGICAL         , &
                                           0, MPI_COMM_WORLD, ierr  )
            
            
            ! Simulation algorithm parameters
            CALL MPI_BCAST(model_eqns    , 1, MPI_INTEGER         , &
                                           0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(num_fluids    , 1, MPI_INTEGER         , &
                                           0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(adv_alphan    , 1, MPI_LOGICAL         , &
                                           0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(mpp_lim       , 1, MPI_LOGICAL         , &
                                           0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(time_stepper  , 1, MPI_INTEGER         , &
                                           0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(weno_vars     , 1, MPI_INTEGER         , &
                                           0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(weno_order    , 1, MPI_INTEGER         , &
                                           0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(weno_eps      , 1, MPI_DOUBLE_PRECISION, &
                                           0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(char_decomp   , 1, MPI_LOGICAL         , &
                                           0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(mapped_weno   , 1, MPI_LOGICAL         , &
                                           0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(mp_weno       , 1, MPI_LOGICAL         , &
                                           0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(weno_avg      , 1, MPI_LOGICAL         , &
                                           0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(weno_Re_flux  , 1, MPI_LOGICAL         , &
                                           0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(riemann_solver, 1, MPI_INTEGER         , &
                                           0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(wave_speeds   , 1, MPI_INTEGER         , &
                                           0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(avg_state     , 1, MPI_INTEGER         , &
                                           0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(commute_err   , 1, MPI_LOGICAL         , &
                                           0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(split_err     , 1, MPI_LOGICAL         , &
                                           0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(alt_crv       , 1, MPI_LOGICAL         , &
                                           0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(alt_soundspeed, 1, MPI_LOGICAL         , &
                                           0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(regularization, 1, MPI_LOGICAL         , &
                                           0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(reg_eps       , 1, MPI_DOUBLE_PRECISION, &
                                           0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(null_weights  , 1, MPI_LOGICAL         , &
                                           0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(mixture_err   , 1, MPI_LOGICAL         , &
                                           0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(tvd_riemann_flux, 1, MPI_LOGICAL         , &
                                           0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(tvd_rhs_flux, 1, MPI_LOGICAL         , &
                                           0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(tvd_wave_speeds, 1, MPI_LOGICAL         , &
                                           0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(flux_lim      , 1, MPI_INTEGER         , &
                                           0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(We_riemann_flux, 1, MPI_LOGICAL         , &
                                           0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(We_rhs_flux, 1, MPI_LOGICAL         , &
                                           0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(We_src, 1, MPI_LOGICAL         , &
                                           0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(We_wave_speeds, 1, MPI_LOGICAL         , &
                                           0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(lsq_deriv, 1, MPI_LOGICAL         , &
                                           0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(parallel_io   , 1, MPI_LOGICAL         , &
                                           0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(precision      , 1, MPI_INTEGER         , &
                                           0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(hypoelasticity, 1, MPI_LOGICAL         , &
                                            0, MPI_COMM_WORLD,ierr   )
            
            CALL MPI_BCAST(bc_x%beg, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(bc_x%end, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(bc_y%beg, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(bc_y%end, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(bc_z%beg, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(bc_z%end, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            
            
            CALL MPI_BCAST(bc_x_glb%beg, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(bc_x_glb%end, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(bc_y_glb%beg, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(bc_y_glb%end, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(bc_z_glb%beg, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(bc_z_glb%end, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            
            ! Fluids physical parameters
            DO i = 1, num_fluids_max
                CALL MPI_BCAST( fluid_pp(i)%gamma   ,        1      , &
                                MPI_DOUBLE_PRECISION,        0      , &
                                MPI_COMM_WORLD, ierr                  )
                CALL MPI_BCAST( fluid_pp(i)%pi_inf  ,        1      , &
                                MPI_DOUBLE_PRECISION,        0      , &
                                MPI_COMM_WORLD, ierr                  )
                CALL MPI_BCAST( fluid_pp(i)%Re(1)   ,        2      , &
                                MPI_DOUBLE_PRECISION,        0      , &
                                MPI_COMM_WORLD, ierr                  )
                CALL MPI_BCAST( fluid_pp(i)%We(1)   , num_fluids_max, &
                                MPI_DOUBLE_PRECISION,        0      , &
                                MPI_COMM_WORLD, ierr                  )

                CALL MPI_BCAST( fluid_pp(i)%mul0  , 1, &
                                MPI_DOUBLE_PRECISION, 0, &
                                MPI_COMM_WORLD, ierr     )
                CALL MPI_BCAST( fluid_pp(i)%ss  , 1, &
                                MPI_DOUBLE_PRECISION, 0, &
                                MPI_COMM_WORLD, ierr     )
                CALL MPI_BCAST( fluid_pp(i)%pv  , 1, &
                                MPI_DOUBLE_PRECISION, 0, &
                                MPI_COMM_WORLD, ierr     )
                CALL MPI_BCAST( fluid_pp(i)%gamma_v  , 1, &
                                MPI_DOUBLE_PRECISION, 0, &
                                MPI_COMM_WORLD, ierr     )
                CALL MPI_BCAST( fluid_pp(i)%M_v  , 1, &
                                MPI_DOUBLE_PRECISION, 0, &
                                MPI_COMM_WORLD, ierr     )
                CALL MPI_BCAST( fluid_pp(i)%mu_v  , 1, &
                                MPI_DOUBLE_PRECISION, 0, &
                                MPI_COMM_WORLD, ierr     )
                CALL MPI_BCAST( fluid_pp(i)%k_v  , 1, &
                                MPI_DOUBLE_PRECISION, 0, &
                                MPI_COMM_WORLD, ierr     )                            
                Call MPI_BCAST( fluid_pp(i)%G   , 1, &
                                MPI_DOUBLE_PRECISION, 0, &
                                MPI_COMM_WORLD, ierr    )
            END DO

            !Tait EOS
            CALL MPI_BCAST( pref,1,             &
                        MPI_DOUBLE_PRECISION,0, &
                        MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST( rhoref,1,           &
                        MPI_DOUBLE_PRECISION,0, &
                        MPI_COMM_WORLD,ierr)

            !Bubble modeling
            CALL MPI_BCAST( bubbles,1,          &
                        MPI_LOGICAL,0,          &
                        MPI_COMM_WORLD,ierr  )
            CALL MPI_BCAST( bubble_model,1,            &
                        MPI_INTEGER,0, &
                        MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST( polytropic,1,          &
                        MPI_LOGICAL,0,          &
                        MPI_COMM_WORLD,ierr  )
            CALL MPI_BCAST( thermal,1,            &
                        MPI_INTEGER,0, &
                        MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST( R0ref,1,            &
                        MPI_DOUBLE_PRECISION,0, &
                        MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST( nb,1,            &
                        MPI_INTEGER,0, &
                        MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST( R0_type,1,            &
                        MPI_INTEGER,0, &
                        MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST( Web,1,            &
                        MPI_DOUBLE_PRECISION,0, &
                        MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST( Ca,1,            &
                        MPI_DOUBLE_PRECISION,0, &
                        MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST( Re_inv,1,            &
                        MPI_DOUBLE_PRECISION,0, &
                        MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST( polydisperse,1,          &
                        MPI_LOGICAL,0,          &
                        MPI_COMM_WORLD,ierr  )
            CALL MPI_BCAST( poly_sigma,1,            &
                        MPI_DOUBLE_PRECISION,0, &
                        MPI_COMM_WORLD,ierr)

            CALL MPI_BCAST( qbmm,1,          &
                        MPI_LOGICAL,0,          &
                        MPI_COMM_WORLD,ierr  )
            CALL MPI_BCAST( nnode,1,            &
                        MPI_INTEGER,0, &
                        MPI_COMM_WORLD,ierr)

            !Acoustic monopole
            CALL MPI_BCAST( monopole,1,          &
                        MPI_LOGICAL,0,          &
                        MPI_COMM_WORLD,ierr  )
            CALL MPI_BCAST( num_mono,1,          &
                        MPI_INTEGER,0,          &
                        MPI_COMM_WORLD,ierr  )
                
            DO j = 1,num_probes_max
                DO i = 1,3
                    CALL MPI_BCAST( mono(j)%loc(i)   ,              1      , &
                                MPI_DOUBLE_PRECISION,        0      , &
                                MPI_COMM_WORLD, ierr                  )
                END DO
                CALL MPI_BCAST( mono(j)%mag   ,              1      , &
                    MPI_DOUBLE_PRECISION,        0      , &
                    MPI_COMM_WORLD, ierr                  )
                CALL MPI_BCAST( mono(j)%length   ,              1      , &
                    MPI_DOUBLE_PRECISION,        0      , &
                    MPI_COMM_WORLD, ierr                  )
                CALL MPI_BCAST( mono(j)%delay,              1      , &
                    MPI_DOUBLE_PRECISION,        0      , &
                    MPI_COMM_WORLD, ierr                  )
                CALL MPI_BCAST( mono(j)%dir   ,              1      , &
                    MPI_DOUBLE_PRECISION,        0      , &
                    MPI_COMM_WORLD, ierr                  )
                CALL MPI_BCAST( mono(j)%npulse   ,              1      , &
                    MPI_DOUBLE_PRECISION,        0      , &
                    MPI_COMM_WORLD, ierr                  )
                CALL MPI_BCAST( mono(j)%pulse   ,              1      , &
                    MPI_INTEGER,        0      , &
                    MPI_COMM_WORLD, ierr                  )
                CALL MPI_BCAST( mono(j)%support   ,              1      , &
                    MPI_INTEGER,        0      , &
                    MPI_COMM_WORLD, ierr                  )
            END DO


            CALL MPI_BCAST(com_wrt(1), num_fluids_max, MPI_LOGICAL, &
                        0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(cb_wrt(1), num_fluids_max, MPI_LOGICAL, &
                        0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(fd_order, 1, MPI_INTEGER, &
                        0, MPI_COMM_WORLD, ierr)
            
            CALL MPI_BCAST(num_probes, 1, MPI_INTEGER, &
                        0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(probe_wrt, 1, MPI_LOGICAL, &
                        0, MPI_COMM_WORLD, ierr)
            DO i = 1, num_probes_max
                CALL MPI_BCAST(probe(i)%x, 1, MPI_DOUBLE_PRECISION, &
                            0, MPI_COMM_WORLD, ierr)
                CALL MPI_BCAST(probe(i)%y, 1, MPI_DOUBLE_PRECISION, &
                            0, MPI_COMM_WORLD, ierr)
                CALL MPI_BCAST(probe(i)%z, 1, MPI_DOUBLE_PRECISION, &
                            0, MPI_COMM_WORLD, ierr)
            END DO
            
            CALL MPI_BCAST(integral_wrt, 1, MPI_LOGICAL, &
                    0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(num_integrals, 1, MPI_INTEGER, &
                    0, MPI_COMM_WORLD, ierr)
            DO i = 1, num_probes_max
                CALL MPI_BCAST(integral(i)%xmin, 1, MPI_DOUBLE_PRECISION, &
                            0, MPI_COMM_WORLD, ierr)
                CALL MPI_BCAST(integral(i)%xmax, 1, MPI_DOUBLE_PRECISION, &
                            0, MPI_COMM_WORLD, ierr)
                CALL MPI_BCAST(integral(i)%ymin, 1, MPI_DOUBLE_PRECISION, &
                            0, MPI_COMM_WORLD, ierr)
                CALL MPI_BCAST(integral(i)%ymax, 1, MPI_DOUBLE_PRECISION, &
                            0, MPI_COMM_WORLD, ierr)
                CALL MPI_BCAST(integral(i)%zmin, 1, MPI_DOUBLE_PRECISION, &
                            0, MPI_COMM_WORLD, ierr)
                CALL MPI_BCAST(integral(i)%zmax, 1, MPI_DOUBLE_PRECISION, &
                            0, MPI_COMM_WORLD, ierr)
            END DO
            
            CALL MPI_BCAST(threshold_mf(1), 5, MPI_DOUBLE_PRECISION, &
                        0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(moment_order(1), 5, MPI_INTEGER, &
                        0, MPI_COMM_WORLD, ierr)
            
            
        END SUBROUTINE s_mpi_bcast_user_inputs ! -------------------------------
        
        
        
        !>  The purpose of this procedure is to optimally decompose
        !!      the computational domain among the available processors.
        !!      This is performed by attempting to award each processor,
        !!      in each of the coordinate directions, approximately the
        !!      same number of cells, and then recomputing the affected
        !!      global parameters.
        SUBROUTINE s_mpi_decompose_computational_domain() ! --------------------
            
            INTEGER :: num_procs_x, num_procs_y, num_procs_z !<
            !! Optimal number of processors in the x-, y- and z-directions
            
            REAL(KIND(0d0)) :: tmp_num_procs_x, tmp_num_procs_y, tmp_num_procs_z !<
            !! Non-optimal number of processors in the x-, y- and z-directions
            
            REAL(KIND(0d0)) :: fct_min !<
            !! Processor factorization (fct) minimization parameter
            
            INTEGER :: MPI_COMM_CART !<
            !! Cartesian processor topology communicator
            
            INTEGER :: rem_cells !<
            !! Remaining number of cells, in a particular coordinate direction,
            !! after the majority is divided up among the available processors
 

            INTEGER :: i,j !< Generic loop iterators
            
            IF (num_procs == 1 .AND. parallel_io) THEN
                DO i = 1, num_dims
                    start_idx(i) = 0
                END DO
                RETURN
            END IF
            
            ! 3D Cartesian Processor Topology ==================================
            IF(n > 0) THEN
                
                IF(p > 0) THEN
                    
                    IF (cyl_coord .AND. p > 0) THEN
                    ! Implement pencil processor blocking if using cylindrical coordinates so
                    ! that all cells in azimuthal direction are stored on a single processor.
                    ! This is necessary for efficient application of Fourier filter near axis.
                        
                        ! Initial values of the processor factorization optimization
                        num_procs_x =  1
                        num_procs_y = num_procs
                        num_procs_z =  1
                        ierr        = -1
                        
                        ! Computing minimization variable for these initial values
                        tmp_num_procs_x = num_procs_x
                        tmp_num_procs_y = num_procs_y
                        tmp_num_procs_z = num_procs_z
                        fct_min         = 10d0*ABS((m+1)/tmp_num_procs_x &
                                                  -(n+1)/tmp_num_procs_y)
                        
                        ! Searching for optimal computational domain distribution
                        DO i = 1, num_procs
                            
                            IF(        MOD(num_procs,i) == 0        &
                                               .AND.                &
                                (m+1)/i >= num_stcls_min*weno_order ) THEN
                                
                                tmp_num_procs_x = i
                                tmp_num_procs_y = num_procs / i
                                
                                IF( fct_min >= ABS((m+1)/tmp_num_procs_x  &
                                                 -(n+1)/tmp_num_procs_y) &
                                                    .AND.                &
                                            (n+1)/tmp_num_procs_y        &
                                                     >=                  &
                                           num_stcls_min*weno_order      ) THEN
                                    
                                    num_procs_x = i
                                    num_procs_y = num_procs/i
                                    fct_min     = ABS((m+1)/tmp_num_procs_x &
                                                     -(n+1)/tmp_num_procs_y)
                                    ierr        = 0
                                    
                                END IF
                                
                            END IF
                            
                        END DO

                    ELSE

                        ! Initial estimate of optimal processor topology
                        num_procs_x =  1
                        num_procs_y =  1
                        num_procs_z =  num_procs
                        ierr        = -1
                        
                        ! Benchmarking the quality of this initial guess
                        tmp_num_procs_x = num_procs_x
                        tmp_num_procs_y = num_procs_y
                        tmp_num_procs_z = num_procs_z
                        fct_min         = 10d0*ABS((m+1)/tmp_num_procs_x  &
                                                  -(n+1)/tmp_num_procs_y) &
                                        + 10d0*ABS((n+1)/tmp_num_procs_y  &
                                                  -(p+1)/tmp_num_procs_z)
                        
                        ! Optimization of the initial processor topology
                        DO i = 1, num_procs
                            
                            IF(        MOD(num_procs,i) == 0        &
                                               .AND.                &
                                (m+1)/i >= num_stcls_min*weno_order ) THEN
                                
                                DO j = 1, num_procs/i
                                    
                                    IF(       MOD(num_procs/i,j) == 0       &
                                                       .AND.                &
                                        (n+1)/j >= num_stcls_min*weno_order ) THEN
                                        
                                        tmp_num_procs_x = i
                                        tmp_num_procs_y = j
                                        tmp_num_procs_z = num_procs/(i*j)
                                        
                                        IF( fct_min >= ABS((m+1)/tmp_num_procs_x    &
                                                         -(n+1)/tmp_num_procs_y)   &
                                                    + ABS((n+1)/tmp_num_procs_y    &
                                                         -(p+1)/tmp_num_procs_z)   &
                                                            .AND.                  &
                                                    (p+1)/tmp_num_procs_z          &
                                                             >=                    &
                                                   num_stcls_min*weno_order      ) &
                                                            THEN
                                           
                                           num_procs_x = i
                                           num_procs_y = j
                                           num_procs_z = num_procs/(i*j)
                                           fct_min     = ABS((m+1)/tmp_num_procs_x &
                                                            -(n+1)/tmp_num_procs_y)&
                                                       + ABS((n+1)/tmp_num_procs_y &
                                                            -(p+1)/tmp_num_procs_z)
                                           ierr        = 0
                                           
                                        END IF
                                        
                                    END IF
                                    
                                END DO
                                
                            END IF
                            
                        END DO
                    
                    END IF
                    
                    ! Verifying that a valid decomposition of the computational
                    ! domain has been established. If not, the simulation exits.
                    IF(proc_rank == 0 .AND. ierr == -1) THEN
                        PRINT '(A)', 'Unsupported combination of values '  // &
                                     'of num_procs, m, n, p and '          // &
                                     'weno_order. Exiting ...'
                        CALL s_mpi_abort()
                    END IF
                    
                    ! Creating new communicator using the Cartesian topology
                    CALL MPI_CART_CREATE( MPI_COMM_WORLD, 3, (/ num_procs_x, &
                                          num_procs_y, num_procs_z /),       &
                                          (/ .TRUE., .TRUE., .TRUE. /),      &
                                          .FALSE., MPI_COMM_CART, ierr       )
                    
                    ! Finding the Cartesian coordinates of the local process
                    CALL MPI_CART_COORDS( MPI_COMM_CART, proc_rank, 3, &
                                          proc_coords, ierr            )
            ! END: 3D Cartesian Processor Topology =============================
                    
                    
            ! Global Parameters for z-direction ================================
                    
                    ! Number of remaining cells
                    rem_cells = MOD(p+1,num_procs_z)
                    
                    ! Optimal number of cells per processor
                    p = (p+1)/num_procs_z - 1
                    
                    ! Distributing the remaining cells
                    DO i = 1, rem_cells
                        IF(proc_coords(3) == i-1) THEN
                            p = p + 1; EXIT
                        END IF
                    END DO
                    
                    ! Boundary condition at the beginning
                    IF(proc_coords(3) > 0 .OR. bc_z%beg == -1) THEN
                        proc_coords(3) = proc_coords(3) - 1
                        CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, &
                                            bc_z%beg, ierr              )
                        proc_coords(3) = proc_coords(3) + 1
                    END IF
                    
                    ! Boundary condition at the end
                    IF(proc_coords(3) < num_procs_z-1 .OR. bc_z%end == -1) THEN
                        proc_coords(3) = proc_coords(3) + 1
                        CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, &
                                            bc_z%end, ierr              )
                        proc_coords(3) = proc_coords(3) - 1
                    END IF
                    
                    IF (parallel_io) THEN
                        IF (proc_coords(3) < rem_cells) THEN
                            start_idx(3) = (p+1) * proc_coords(3)
                        ELSE
                            start_idx(3) = (p+1) * proc_coords(3) + rem_cells
                        END IF
                    END IF
            ! ==================================================================
                    
                    
            ! 2D Cartesian Processor Topology ==================================
                ELSE
                    
                    ! Initial estimate of optimal processor topology
                    num_procs_x =  1
                    num_procs_y =  num_procs
                    ierr        = -1
                    
                    ! Benchmarking the quality of this initial guess
                    tmp_num_procs_x = num_procs_x
                    tmp_num_procs_y = num_procs_y
                    fct_min         = 10d0*ABS((m+1)/tmp_num_procs_x &
                                              -(n+1)/tmp_num_procs_y)
                    
                    ! Optimization of the initial processor topology
                    DO i = 1, num_procs
                        
                        IF(        MOD(num_procs,i) == 0        &
                                           .AND.                &
                            (m+1)/i >= num_stcls_min*weno_order ) THEN
                            
                            tmp_num_procs_x = i
                            tmp_num_procs_y = num_procs/i
                            
                            IF( fct_min >= ABS((m+1)/tmp_num_procs_x  &
                                             -(n+1)/tmp_num_procs_y) &
                                                .AND.                &
                                        (n+1)/tmp_num_procs_y        &
                                                 >=                  &
                                       num_stcls_min*weno_order      ) THEN
                                
                                num_procs_x = i
                                num_procs_y = num_procs/i
                                fct_min     = ABS((m+1)/tmp_num_procs_x &
                                                 -(n+1)/tmp_num_procs_y)
                                ierr        = 0
                                
                            END IF
                            
                        END IF
                        
                    END DO

                    ! Verifying that a valid decomposition of the computational
                    ! domain has been established. If not, the simulation exits.
                    IF(proc_rank == 0 .AND. ierr == -1) THEN
                        PRINT '(A)', 'Unsupported combination of values ' // &
                                     'of num_procs, m, n and '            // &
                                     'weno_order. Exiting ...'
                        CALL s_mpi_abort()
                    END IF

                    ! Creating new communicator using the Cartesian topology
                    CALL MPI_CART_CREATE( MPI_COMM_WORLD, 2, (/ num_procs_x, &
                                          num_procs_y /), (/ .TRUE.,         &
                                          .TRUE. /), .FALSE., MPI_COMM_CART, &
                                          ierr                               )
                    
                    ! Finding the Cartesian coordinates of the local process
                    CALL MPI_CART_COORDS( MPI_COMM_CART, proc_rank, 2, &
                                          proc_coords, ierr            )
                    
                END IF
            ! END: 2D Cartesian Processor Topology =============================
                
                
            ! Global Parameters for y-direction ================================
                
                ! Number of remaining cells
                rem_cells = MOD(n+1,num_procs_y)
                
                ! Optimal number of cells per processor
                n = (n+1)/num_procs_y - 1
                
                ! Distributing the remaining cells
                DO i = 1, rem_cells
                    IF(proc_coords(2) == i-1) THEN
                        n = n + 1; EXIT
                    END IF
                END DO
                
                ! Boundary condition at the beginning
                IF(proc_coords(2) > 0 .OR. bc_y%beg == -1) THEN
                    proc_coords(2) = proc_coords(2) - 1
                    CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, &
                                        bc_y%beg, ierr              )
                    proc_coords(2) = proc_coords(2) + 1
                END IF
                
                ! Boundary condition at the end
                IF(proc_coords(2) < num_procs_y-1 .OR. bc_y%end == -1) THEN
                    proc_coords(2) = proc_coords(2) + 1
                    CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, &
                                        bc_y%end, ierr              )
                    proc_coords(2) = proc_coords(2) - 1
                END IF
                
                IF (parallel_io) THEN
                    IF (proc_coords(2) < rem_cells) THEN
                        start_idx(2) = (n+1) * proc_coords(2)
                    ELSE
                        start_idx(2) = (n+1) * proc_coords(2) + rem_cells
                    END IF
                END IF
            ! ==================================================================
                
                
            ! 1D Cartesian Processor Topology ==================================
            ELSE
                
                ! Optimal processor topology
                num_procs_x = num_procs
                
                ! Creating new communicator using the Cartesian topology
                CALL MPI_CART_CREATE( MPI_COMM_WORLD, 1, (/ num_procs_x /), &
                                      (/ .TRUE. /), .FALSE., MPI_COMM_CART, &
                                      ierr                                  )
                
                ! Finding the Cartesian coordinates of the local process
                CALL MPI_CART_COORDS( MPI_COMM_CART, proc_rank, 1, &
                                      proc_coords, ierr            )
                
            END IF
            ! ==================================================================
            
            
            ! Global Parameters for x-direction ================================
            
            ! Number of remaining cells
            rem_cells = MOD(m+1,num_procs_x)
            
            ! Optimal number of cells per processor
            m = (m+1)/num_procs_x - 1
            
            ! Distributing the remaining cells
            DO i = 1, rem_cells
                IF(proc_coords(1) == i-1) THEN
                    m = m + 1; EXIT
                END IF
            END DO
            
            ! Boundary condition at the beginning
            IF(proc_coords(1) > 0 .OR. bc_x%beg == -1) THEN
                proc_coords(1) = proc_coords(1) - 1
                CALL MPI_CART_RANK(MPI_COMM_CART, proc_coords, bc_x%beg, ierr)
                proc_coords(1) = proc_coords(1) + 1
            END IF
            
            ! Boundary condition at the end
            IF(proc_coords(1) < num_procs_x-1 .OR. bc_x%end == -1) THEN
                proc_coords(1) = proc_coords(1) + 1
                CALL MPI_CART_RANK(MPI_COMM_CART, proc_coords, bc_x%end, ierr)
                proc_coords(1) = proc_coords(1) - 1
            END IF
            
            IF (parallel_io) THEN
                IF (proc_coords(1) < rem_cells) THEN
                    start_idx(1) = (m+1) * proc_coords(1)
                ELSE
                    start_idx(1) = (m+1) * proc_coords(1) + rem_cells
                END IF
            END IF
            ! ==================================================================
            
            
        END SUBROUTINE s_mpi_decompose_computational_domain ! ------------------
        
        
        
        !>  The goal of this procedure is to populate the buffers of
        !!      the grid variables by communicating with the neighboring
        !!      processors. Note that only the buffers of the cell-width
        !!      distributions are handled in such a way. This is because
        !!      the buffers of cell-boundary locations may be calculated
        !!      directly from those of the cell-width distributions.        
        !!  @param mpi_dir MPI communication coordinate direction
        !!  @param pbc_loc Processor boundary condition (PBC) location 
        SUBROUTINE s_mpi_sendrecv_grid_variables_buffers(mpi_dir, pbc_loc) ! ---

            INTEGER, INTENT(IN) :: mpi_dir
            INTEGER, INTENT(IN) :: pbc_loc
            
            
            ! MPI Communication in x-direction =================================
            IF(mpi_dir == 1) THEN
                
                IF(pbc_loc == -1) THEN      ! PBC at the beginning
                    
                    IF(bc_x%end >= 0) THEN      ! PBC at the beginning and end
                        
                        ! Send/receive buffer to/from bc_x%end/bc_x%beg
                        CALL MPI_SENDRECV(                              &
                                dx(m-buff_size+1), buff_size,           &
                                MPI_DOUBLE_PRECISION, bc_x%end, 0,      &
                                dx(-buff_size), buff_size,              &
                                MPI_DOUBLE_PRECISION, bc_x%beg, 0,      &
                                MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
                        
                    ELSE                        ! PBC at the beginning only
                        
                        ! Send/receive buffer to/from bc_x%beg/bc_x%beg
                        CALL MPI_SENDRECV(                              &
                                dx(0), buff_size,                       &
                                MPI_DOUBLE_PRECISION, bc_x%beg, 1,      &
                                dx(-buff_size), buff_size,              &
                                MPI_DOUBLE_PRECISION, bc_x%beg, 0,      &
                                MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
                        
                    END IF
                    
                ELSE                        ! PBC at the end
                    
                    IF(bc_x%beg >= 0) THEN      ! PBC at the end and beginning
                        
                        ! Send/receive buffer to/from bc_x%beg/bc_x%end
                        CALL MPI_SENDRECV(                              &
                                dx(0), buff_size,                       &
                                MPI_DOUBLE_PRECISION, bc_x%beg, 1,      &
                                dx(m+1), buff_size,                     &
                                MPI_DOUBLE_PRECISION, bc_x%end, 1,      &
                                MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
                        
                    ELSE                        ! PBC at the end only
                        
                        ! Send/receive buffer to/from bc_x%end/bc_x%end
                        CALL MPI_SENDRECV(                              &
                                dx(m-buff_size+1), buff_size,           &
                                MPI_DOUBLE_PRECISION, bc_x%end, 0,      &
                                dx(m+1), buff_size,                     &
                                MPI_DOUBLE_PRECISION, bc_x%end, 1,      &
                                MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
                        
                    END IF
                    
                END IF
            ! END: MPI Communication in x-direction ============================
                
                
            ! MPI Communication in y-direction =================================
            ELSEIF(mpi_dir == 2) THEN
                
                IF(pbc_loc == -1) THEN      ! PBC at the beginning
                    
                    IF(bc_y%end >= 0) THEN      ! PBC at the beginning and end
                        
                        ! Send/receive buffer to/from bc_y%end/bc_y%beg
                        CALL MPI_SENDRECV(                              &
                                dy(n-buff_size+1), buff_size,           &
                                MPI_DOUBLE_PRECISION, bc_y%end, 0,      &
                                dy(-buff_size), buff_size,              &
                                MPI_DOUBLE_PRECISION, bc_y%beg, 0,      &
                                MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
                        
                    ELSE                        ! PBC at the beginning only
                        
                        ! Send/receive buffer to/from bc_y%beg/bc_y%beg
                        CALL MPI_SENDRECV(                              &
                                dy(0), buff_size,                       &
                                MPI_DOUBLE_PRECISION, bc_y%beg, 1,      &
                                dy(-buff_size), buff_size,              &
                                MPI_DOUBLE_PRECISION, bc_y%beg, 0,      &
                                MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
                        
                    END IF
                    
                ELSE                        ! PBC at the end
                    
                    IF(bc_y%beg >= 0) THEN      ! PBC at the end and beginning
                        
                        ! Send/receive buffer to/from bc_y%beg/bc_y%end
                        CALL MPI_SENDRECV(                              &
                                dy(0), buff_size,                       &
                                MPI_DOUBLE_PRECISION, bc_y%beg, 1,      &
                                dy(n+1), buff_size,                     &
                                MPI_DOUBLE_PRECISION, bc_y%end, 1,      &
                                MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
                        
                    ELSE                        ! PBC at the end only
                        
                        ! Send/receive buffer to/from bc_y%end/bc_y%end
                        CALL MPI_SENDRECV(                              &
                                dy(n-buff_size+1), buff_size,           &
                                MPI_DOUBLE_PRECISION, bc_y%end, 0,      &
                                dy(n+1), buff_size,                     &
                                MPI_DOUBLE_PRECISION, bc_y%end, 1,      &
                                MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
                        
                    END IF
                    
                END IF
            ! END: MPI Communication in y-direction ============================
                
                
            ! MPI Communication in z-direction =================================
            ELSE
                
                IF(pbc_loc == -1) THEN      ! PBC at the beginning
                    
                    IF(bc_z%end >= 0) THEN      ! PBC at the beginning and end
                        
                        ! Send/receive buffer to/from bc_z%end/bc_z%beg
                        CALL MPI_SENDRECV(                              &
                                dz(p-buff_size+1), buff_size,           &
                                MPI_DOUBLE_PRECISION, bc_z%end, 0,      &
                                dz(-buff_size), buff_size,              &
                                MPI_DOUBLE_PRECISION, bc_z%beg, 0,      &
                                MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
                        
                    ELSE                        ! PBC at the beginning only
                        
                        ! Send/receive buffer to/from bc_z%beg/bc_z%beg
                        CALL MPI_SENDRECV(                              &
                                dz(0), buff_size,                       &
                                MPI_DOUBLE_PRECISION, bc_z%beg, 1,      &
                                dz(-buff_size), buff_size,              &
                                MPI_DOUBLE_PRECISION, bc_z%beg, 0,      &
                                MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
                        
                    END IF
                    
                ELSE                        ! PBC at the end
                    
                    IF(bc_z%beg >= 0) THEN      ! PBC at the end and beginning
                        
                        ! Send/receive buffer to/from bc_z%beg/bc_z%end
                        CALL MPI_SENDRECV(                              &
                                dz(0), buff_size,                       &
                                MPI_DOUBLE_PRECISION, bc_z%beg, 1,      &
                                dz(p+1), buff_size,                     &
                                MPI_DOUBLE_PRECISION, bc_z%end, 1,      &
                                MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
                        
                    ELSE                        ! PBC at the end only
                        
                        ! Send/receive buffer to/from bc_z%end/bc_z%end
                        CALL MPI_SENDRECV(                              &
                                dz(p-buff_size+1), buff_size,           &
                                MPI_DOUBLE_PRECISION, bc_z%end, 0,      &
                                dz(p+1), buff_size,                     &
                                MPI_DOUBLE_PRECISION, bc_z%end, 1,      &
                                MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
                        
                    END IF
                    
                END IF
                
            END IF
            ! END: MPI Communication in z-direction ============================
            
            
        END SUBROUTINE s_mpi_sendrecv_grid_variables_buffers ! -----------------
        
        
        
        !>  The goal of this subroutine is to determine the global
        !!      extrema of the stability criteria in the computational
        !!      domain. This is performed by sifting through the local
        !!      extrema of each stability criterion. Note that each of
        !!      the local extrema is from a single process, within its
        !!      assigned section of the computational domain. Finally,
        !!      note that the global extrema values are only bookkeept
        !!      on the rank 0 processor.
        !!  @param icfl_max_loc Local maximum ICFL stability criterion
        !!  @param vcfl_max_loc Local maximum VCFL stability criterion
        !!  @param ccfl_max_loc Local maximum CCFL stability criterion
        !!  @param Rc_min_loc Local minimum Rc stability criterion
        !!  @param icfl_max_glb Global maximum ICFL stability criterion
        !!  @param vcfl_max_glb Global maximum VCFL stability criterion
        !!  @param ccfl_max_glb Global maximum CCFL stability criterion
        !!  @param Rc_min_glb Global minimum Rc stability criterion
        SUBROUTINE s_mpi_reduce_stability_criteria_extrema( icfl_max_loc, & ! --
                                                            vcfl_max_loc, &
                                                            ccfl_max_loc, &
                                                              Rc_min_loc, &
                                                            icfl_max_glb, &
                                                            vcfl_max_glb, &
                                                            ccfl_max_glb, &
                                                              Rc_min_glb  )
           
            REAL(KIND(0d0)), INTENT(IN)  :: icfl_max_loc
            REAL(KIND(0d0)), INTENT(IN)  :: vcfl_max_loc
            REAL(KIND(0d0)), INTENT(IN)  :: ccfl_max_loc
            REAL(KIND(0d0)), INTENT(IN)  ::   Rc_min_loc
            
            REAL(KIND(0d0)), INTENT(OUT) :: icfl_max_glb
            REAL(KIND(0d0)), INTENT(OUT) :: vcfl_max_glb
            REAL(KIND(0d0)), INTENT(OUT) :: ccfl_max_glb
            REAL(KIND(0d0)), INTENT(OUT) ::   Rc_min_glb
            
            
            ! Reducing local extrema of ICFL, VCFL, CCFL and Rc numbers to their
            ! global extrema and bookkeeping the results on the rank 0 processor
            CALL MPI_REDUCE( icfl_max_loc, icfl_max_glb   , 1, &
                             MPI_DOUBLE_PRECISION, MPI_MAX, 0, &
                             MPI_COMM_WORLD, ierr              )
            
            IF(ANY(Re_size > 0)) THEN
                CALL MPI_REDUCE( vcfl_max_loc, vcfl_max_glb   , 1, &
                                 MPI_DOUBLE_PRECISION, MPI_MAX, 0, &
                                 MPI_COMM_WORLD, ierr              )
                CALL MPI_REDUCE( Rc_min_loc  , Rc_min_glb     , 1, &
                                 MPI_DOUBLE_PRECISION, MPI_MIN, 0, &
                                 MPI_COMM_WORLD, ierr              )
            END IF
            
            IF(We_size > 0) THEN
                CALL MPI_REDUCE( ccfl_max_loc, ccfl_max_glb   , 1, &
                                 MPI_DOUBLE_PRECISION, MPI_MAX, 0, &
                                 MPI_COMM_WORLD, ierr              )
            END IF
            
            
        END SUBROUTINE s_mpi_reduce_stability_criteria_extrema ! ---------------
        
        
        
        !>  The following subroutine takes the input local variable
        !!      from all processors and reduces to the sum of all
        !!      values. The reduced variable is recorded back onto the 
        !!      original local variable on each processor. 
        !!  @param var_loc Some variable containing the local value which should be
        !!  reduced amongst all the processors in the communicator.
        !!  @param var_glb The globally reduced value
        SUBROUTINE s_mpi_allreduce_sum(var_loc, var_glb) ! ---------------------

            REAL(KIND(0d0)), INTENT(IN) :: var_loc
            REAL(KIND(0d0)), INTENT(OUT) :: var_glb

            ! Performing the reduction procedure
            CALL MPI_ALLREDUCE(var_loc, var_glb, 1, MPI_DOUBLE_PRECISION, &
                        MPI_SUM, MPI_COMM_WORLD, ierr)

        END SUBROUTINE s_mpi_allreduce_sum ! -----------------------------------



        !>  The following subroutine takes the input local variable
        !!      from all processors and reduces to the minimum of all
        !!      values. The reduced variable is recorded back onto the 
        !!      original local variable on each processor.
        !!  @param var_loc Some variable containing the local value which should be
        !!  reduced amongst all the processors in the communicator.
        !!  @param var_glb The globally reduced value
        SUBROUTINE s_mpi_allreduce_min(var_loc, var_glb) ! ---------------------

            REAL(KIND(0d0)), INTENT(IN) :: var_loc
            REAL(KIND(0d0)), INTENT(OUT) :: var_glb

            ! Performing the reduction procedure
            CALL MPI_ALLREDUCE(var_loc, var_glb, 1, MPI_DOUBLE_PRECISION, &
                        MPI_MIN, MPI_COMM_WORLD, ierr)

        END SUBROUTINE s_mpi_allreduce_min ! -----------------------------------



        !>  The following subroutine takes the input local variable
        !!      from all processors and reduces to the maximum of all
        !!      values. The reduced variable is recorded back onto the 
        !!      original local variable on each processor.
        !!  @param var_loc Some variable containing the local value which should be
        !!  reduced amongst all the processors in the communicator.
        !!  @param var_glb The globally reduced value
        SUBROUTINE s_mpi_allreduce_max(var_loc, var_glb) ! ---------------------

            REAL(KIND(0d0)), INTENT(IN) :: var_loc
            REAL(KIND(0d0)), INTENT(OUT) :: var_glb

            ! Performing the reduction procedure
            CALL MPI_ALLREDUCE(var_loc, var_glb, 1, MPI_DOUBLE_PRECISION, &
                        MPI_MAX, MPI_COMM_WORLD, ierr)

        END SUBROUTINE s_mpi_allreduce_max ! -----------------------------------



        !>  The goal of this procedure is to populate the buffers of
        !!      the cell-average conservative variables by communicating
        !!      with the neighboring processors.
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param mpi_dir MPI communication coordinate direction
        !!  @param pbc_loc Processor boundary condition (PBC) location
        SUBROUTINE s_mpi_sendrecv_conservative_variables_buffers( q_cons_vf, &
                                                                  mpi_dir,   &
                                                                  pbc_loc    )

            
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(INOUT) :: q_cons_vf
            INTEGER, INTENT(IN) :: mpi_dir
            INTEGER, INTENT(IN) :: pbc_loc
            

            INTEGER :: i,j,k,l,r !< Generic loop iterators
            
            
            ! MPI Communication in x-direction =================================
            IF(mpi_dir == 1) THEN
                
                IF(pbc_loc == -1) THEN      ! PBC at the beginning
                    
                    IF(bc_x%end >= 0) THEN      ! PBC at the beginning and end
                        
                        ! Packing buffer to be sent to bc_x%end
                        DO l = 0, p
                           DO k = 0, n
                              DO j = m-buff_size+1, m
                                 DO i = 1, sys_size
                                    r = (i-1) + sys_size * &
                                        ((j-m-1) + buff_size*((k+1) + (n+1)*l))
                                    q_cons_buff_send(r) = q_cons_vf(i)%sf(j,k,l)
                                 END DO
                              END DO
                           END DO
                        END DO
                        
                        ! Send/receive buffer to/from bc_x%end/bc_x%beg
                        CALL MPI_SENDRECV(                              &
                                q_cons_buff_send(0),                    &
                                buff_size*sys_size*(n+1)*(p+1),      &
                                MPI_DOUBLE_PRECISION, bc_x%end, 0,      &
                                q_cons_buff_recv(0),                    &
                                buff_size*sys_size*(n+1)*(p+1),      &
                                MPI_DOUBLE_PRECISION, bc_x%beg, 0,      &
                                MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
                        
                    ELSE                        ! PBC at the beginning only
                        
                        ! Packing buffer to be sent to bc_x%beg
                        DO l = 0, p
                           DO k = 0, n
                              DO j = 0, buff_size-1
                                 DO i = 1, sys_size
                                    r = (i-1) + sys_size * &
                                        (j + buff_size*(k + (n+1)*l))
                                    q_cons_buff_send(r) = q_cons_vf(i)%sf(j,k,l)
                                 END DO
                              END DO
                           END DO
                        END DO
                        
                        ! Send/receive buffer to/from bc_x%beg/bc_x%beg
                        CALL MPI_SENDRECV(                              &
                                q_cons_buff_send(0),                    &
                                buff_size*sys_size*(n+1)*(p+1),      &
                                MPI_DOUBLE_PRECISION, bc_x%beg, 1,      &
                                q_cons_buff_recv(0),                    &
                                buff_size*sys_size*(n+1)*(p+1),      &
                                MPI_DOUBLE_PRECISION, bc_x%beg, 0,      &
                                MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
                        
                    END IF
                    
                    ! Unpacking buffer received from bc_x%beg
                    DO l = 0, p
                        DO k = 0, n
                            DO j = -buff_size, -1
                                DO i = 1, sys_size
                                    r = (i-1) + sys_size * &
                                        (j + buff_size*((k+1) + (n+1)*l))
                                    q_cons_vf(i)%sf(j,k,l) = q_cons_buff_recv(r)
                                END DO
                            END DO
                        END DO
                    END DO
                    
                ELSE                        ! PBC at the end
                    
                    IF(bc_x%beg >= 0) THEN      ! PBC at the end and beginning
                        
                        ! Packing buffer to be sent to bc_x%beg
                        DO l = 0, p
                           DO k = 0, n
                              DO j = 0, buff_size-1
                                 DO i = 1, sys_size
                                    r = (i-1) + sys_size * &
                                        (j + buff_size*(k + (n+1)*l))
                                    q_cons_buff_send(r) = q_cons_vf(i)%sf(j,k,l)
                                 END DO
                              END DO
                           END DO
                        END DO
                        
                        ! Send/receive buffer to/from bc_x%beg/bc_x%end
                        CALL MPI_SENDRECV(                              &
                                q_cons_buff_send(0),                    &
                                buff_size*sys_size*(n+1)*(p+1),      &
                                MPI_DOUBLE_PRECISION, bc_x%beg, 1,      &
                                q_cons_buff_recv(0),                    &
                                buff_size*sys_size*(n+1)*(p+1),      &
                                MPI_DOUBLE_PRECISION, bc_x%end, 1,      &
                                MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
                        
                    ELSE                        ! PBC at the end only
                        
                        ! Packing buffer to be sent to bc_x%end
                        DO l = 0, p
                           DO k = 0, n
                              DO j = m-buff_size+1, m
                                 DO i = 1, sys_size
                                    r = (i-1) + sys_size * &
                                        ((j-m-1) + buff_size*((k+1) + (n+1)*l))
                                    q_cons_buff_send(r) = q_cons_vf(i)%sf(j,k,l)
                                 END DO
                              END DO
                           END DO
                        END DO
                        
                        ! Send/receive buffer to/from bc_x%end/bc_x%end
                        CALL MPI_SENDRECV(                              &
                                q_cons_buff_send(0),                    &
                                buff_size*sys_size*(n+1)*(p+1),      &
                                MPI_DOUBLE_PRECISION, bc_x%end, 0,      &
                                q_cons_buff_recv(0),                    &
                                buff_size*sys_size*(n+1)*(p+1),      &
                                MPI_DOUBLE_PRECISION, bc_x%end, 1,      &
                                MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
                        
                    END IF
                    
                    ! Unpacking buffer received from bc_x%end
                    DO l = 0, p
                        DO k = 0, n
                            DO j = m+1, m+buff_size
                                DO i = 1, sys_size
                                    r = (i-1) + sys_size * &
                                        ((j-m-1) + buff_size*(k + (n+1)*l))
                                    q_cons_vf(i)%sf(j,k,l) = q_cons_buff_recv(r)
                                END DO
                            END DO
                        END DO
                    END DO
                    
                END IF
            ! END: MPI Communication in x-direction ============================
                
                
            ! MPI Communication in y-direction =================================
            ELSEIF(mpi_dir == 2) THEN
                
                IF(pbc_loc == -1) THEN      ! PBC at the beginning
                    
                    IF(bc_y%end >= 0) THEN      ! PBC at the beginning and end
                        
                        ! Packing buffer to be sent to bc_y%end
                        DO l = 0, p
                           DO k = n-buff_size+1, n
                              DO j = -buff_size, m+buff_size
                                 DO i = 1, sys_size
                                    r = (i-1) + sys_size * &
                                        ((j+buff_size) + (m+2*buff_size+1) * &
                                        ((k-n+buff_size-1) + buff_size*l))
                                    q_cons_buff_send(r) = q_cons_vf(i)%sf(j,k,l)
                                 END DO
                              END DO
                           END DO
                        END DO
                        
                        ! Send/receive buffer to/from bc_y%end/bc_y%beg
                        CALL MPI_SENDRECV(                                 &
                            q_cons_buff_send(0),                           &
                            buff_size*sys_size*(m+2*buff_size+1)*(p+1), &
                            MPI_DOUBLE_PRECISION, bc_y%end, 0,             &
                            q_cons_buff_recv(0),                           &
                            buff_size*sys_size*(m+2*buff_size+1)*(p+1), &
                            MPI_DOUBLE_PRECISION, bc_y%beg, 0,             &
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr        )
                        
                    ELSE                        ! PBC at the beginning only
                        
                        ! Packing buffer to be sent to bc_y%beg
                        DO l = 0, p
                           DO k = 0, buff_size-1
                              DO j = -buff_size, m+buff_size
                                 DO i = 1, sys_size
                                    r = (i-1) + sys_size * &
                                        ((j+buff_size) + (m+2*buff_size+1) * &
                                        (k + buff_size*l))
                                    q_cons_buff_send(r) = q_cons_vf(i)%sf(j,k,l)
                                 END DO
                              END DO
                           END DO
                        END DO
                        
                        ! Send/receive buffer to/from bc_y%beg/bc_y%beg
                        CALL MPI_SENDRECV(                                 &
                            q_cons_buff_send(0),                           &
                            buff_size*sys_size*(m+2*buff_size+1)*(p+1), &
                            MPI_DOUBLE_PRECISION, bc_y%beg, 1,             &
                            q_cons_buff_recv(0),                           &
                            buff_size*sys_size*(m+2*buff_size+1)*(p+1), &
                            MPI_DOUBLE_PRECISION, bc_y%beg, 0,             &
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr        )
                        
                    END IF
                    
                    ! Unpacking buffer received from bc_y%beg
                    DO l = 0, p
                        DO k = -buff_size, -1
                            DO j = -buff_size, m+buff_size
                                DO i = 1, sys_size
                                    r = (i-1) + sys_size * &
                                        ((j+buff_size) + (m+2*buff_size+1) * &
                                        ((k+buff_size) + buff_size*l))
                                    q_cons_vf(i)%sf(j,k,l) = q_cons_buff_recv(r)
                                END DO
                            END DO
                        END DO
                    END DO
                    
                ELSE                        ! PBC at the end
                    
                    IF(bc_y%beg >= 0) THEN      ! PBC at the end and beginning
                        
                        ! Packing buffer to be sent to bc_y%beg
                        DO l = 0, p
                           DO k = 0, buff_size-1
                              DO j = -buff_size, m+buff_size
                                 DO i = 1, sys_size
                                    r = (i-1) + sys_size * &
                                        ((j+buff_size) + (m+2*buff_size+1) * &
                                        (k + buff_size*l))
                                    q_cons_buff_send(r) = q_cons_vf(i)%sf(j,k,l)
                                 END DO
                              END DO
                           END DO
                        END DO
                        
                        ! Send/receive buffer to/from bc_y%beg/bc_y%end
                        CALL MPI_SENDRECV(                                 &
                            q_cons_buff_send(0),                           &
                            buff_size*sys_size*(m+2*buff_size+1)*(p+1), &
                            MPI_DOUBLE_PRECISION, bc_y%beg, 1,             &
                            q_cons_buff_recv(0),                           &
                            buff_size*sys_size*(m+2*buff_size+1)*(p+1), &
                            MPI_DOUBLE_PRECISION, bc_y%end, 1,             &
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr        )
                        
                    ELSE                        ! PBC at the end only
                        
                        ! Packing buffer to be sent to bc_y%end
                        DO l = 0, p
                           DO k = n-buff_size+1, n
                              DO j = -buff_size, m+buff_size
                                 DO i = 1, sys_size
                                    r = (i-1) + sys_size * &
                                        ((j+buff_size) + (m+2*buff_size+1) * &
                                        ((k-n+buff_size-1) + buff_size*l))
                                    q_cons_buff_send(r) = q_cons_vf(i)%sf(j,k,l)
                                 END DO
                              END DO
                           END DO
                        END DO
                        
                        ! Send/receive buffer to/from bc_y%end/bc_y%end
                        CALL MPI_SENDRECV(                                 &
                            q_cons_buff_send(0),                           &
                            buff_size*sys_size*(m+2*buff_size+1)*(p+1), &
                            MPI_DOUBLE_PRECISION, bc_y%end, 0,             &
                            q_cons_buff_recv(0),                           &
                            buff_size*sys_size*(m+2*buff_size+1)*(p+1), &
                            MPI_DOUBLE_PRECISION, bc_y%end, 1,             &
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr        )
                        
                    END IF
                    
                    ! Unpacking buffer received form bc_y%end
                    DO l = 0, p
                        DO k = n+1, n+buff_size
                            DO j = -buff_size, m+buff_size
                                DO i = 1, sys_size
                                    r = (i-1) + sys_size * &
                                        ((j+buff_size) + (m+2*buff_size+1) * &
                                        ((k-n-1) + buff_size*l))
                                    q_cons_vf(i)%sf(j,k,l) = q_cons_buff_recv(r)
                                END DO
                            END DO
                        END DO
                    END DO
                    
                END IF
            ! END: MPI Communication in y-direction ============================
                
                
            ! MPI Communication in z-direction =================================
            ELSE
                
                IF(pbc_loc == -1) THEN      ! PBC at the beginning
                    
                    IF(bc_z%end >= 0) THEN      ! PBC at the beginning and end
                        
                        ! Packing buffer to be sent to bc_z%end
                        DO l = p-buff_size+1, p
                           DO k = -buff_size, n+buff_size
                              DO j = -buff_size, m+buff_size
                                 DO i = 1, sys_size
                                    r = (i-1) + sys_size * &
                                        ((j+buff_size) + (m+2*buff_size+1) * &
                                        ((k+buff_size) + (n+2*buff_size+1) * &
                                        (l-p+buff_size-1)))
                                    q_cons_buff_send(r) = q_cons_vf(i)%sf(j,k,l)
                                 END DO
                              END DO
                           END DO
                        END DO
                        
                        ! Send/receive buffer to/from bc_z%end/bc_z%beg
                        CALL MPI_SENDRECV(                              &
                               q_cons_buff_send(0),                     &
                               buff_size*sys_size*(m+2*buff_size+1)  &
                                                    *(n+2*buff_size+1), &
                               MPI_DOUBLE_PRECISION, bc_z%end, 0,       &
                               q_cons_buff_recv(0),                     &
                               buff_size*sys_size*(m+2*buff_size+1)  &
                                                    *(n+2*buff_size+1), &
                               MPI_DOUBLE_PRECISION, bc_z%beg, 0,       &
                               MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr  )
                        
                    ELSE                        ! PBC at the beginning only
                        
                        ! Packing buffer to be sent to bc_z%beg
                        DO l = 0, buff_size-1
                           DO k = -buff_size, n+buff_size
                              DO j = -buff_size, m+buff_size
                                 DO i = 1, sys_size
                                    r = (i-1) + sys_size * &
                                        ((j+buff_size) + (m+2*buff_size+1) * &
                                        ((k+buff_size) + (n+2*buff_size+1)*l))
                                    q_cons_buff_send(r) = q_cons_vf(i)%sf(j,k,l)
                                 END DO
                              END DO
                           END DO
                        END DO
                        
                        ! Send/receive buffer to/from bc_z%beg/bc_z%beg
                        CALL MPI_SENDRECV(                              &
                               q_cons_buff_send(0),                     &
                               buff_size*sys_size*(m+2*buff_size+1)  &
                                                    *(n+2*buff_size+1), &
                               MPI_DOUBLE_PRECISION, bc_z%beg, 1,       &
                               q_cons_buff_recv(0),                     &
                               buff_size*sys_size*(m+2*buff_size+1)  &
                                                    *(n+2*buff_size+1), &
                               MPI_DOUBLE_PRECISION, bc_z%beg, 0,       &
                               MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr  )           
                        
                    END IF
                    
                    ! Unpacking buffer from bc_z%beg
                    DO l = -buff_size, -1
                        DO k = -buff_size, n+buff_size
                            DO j = -buff_size, m+buff_size
                                DO i = 1, sys_size
                                    r = (i-1) + sys_size * &
                                        ((j+buff_size) + (m+2*buff_size+1) * &
                                        ((k+buff_size) + (n+2*buff_size+1) * &
                                        (l+buff_size)))
                                    q_cons_vf(i)%sf(j,k,l) = q_cons_buff_recv(r)
                                END DO
                            END DO
                        END DO
                    END DO
                    
                ELSE                        ! PBC at the end
                    
                    IF(bc_z%beg >= 0) THEN      ! PBC at the end and beginning
                        
                        ! Packing buffer to be sent to bc_z%beg
                        DO l = 0, buff_size-1
                           DO k = -buff_size, n+buff_size
                              DO j = -buff_size, m+buff_size
                                 DO i = 1, sys_size
                                    r = (i-1) + sys_size * &
                                        ((j+buff_size) + (m+2*buff_size+1) * &
                                        ((k+buff_size) + (n+2*buff_size+1)*l))
                                    q_cons_buff_send(r) = q_cons_vf(i)%sf(j,k,l)
                                 END DO
                              END DO
                           END DO
                        END DO
                        
                        ! Send/receive buffer to/from bc_z%beg/bc_z%end
                        CALL MPI_SENDRECV(                              &
                               q_cons_buff_send(0),                     &
                               buff_size*sys_size*(m+2*buff_size+1)  &
                                                    *(n+2*buff_size+1), &
                               MPI_DOUBLE_PRECISION, bc_z%beg, 1,       &
                               q_cons_buff_recv(0),                     &
                               buff_size*sys_size*(m+2*buff_size+1)  &
                                                    *(n+2*buff_size+1), &
                               MPI_DOUBLE_PRECISION, bc_z%end, 1,       &
                               MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr  )
                        
                    ELSE                        ! PBC at the end only
                        
                        ! Packing buffer to be sent to bc_z%end
                        DO l = p-buff_size+1, p
                           DO k = -buff_size, n+buff_size
                              DO j = -buff_size, m+buff_size
                                 DO i = 1, sys_size
                                    r = (i-1) + sys_size * &
                                        ((j+buff_size) + (m+2*buff_size+1) * &
                                        ((k+buff_size) + (n+2*buff_size+1) * &
                                        (l-p+buff_size-1)))
                                    q_cons_buff_send(r) = q_cons_vf(i)%sf(j,k,l)
                                 END DO
                              END DO
                           END DO
                        END DO
                        
                        ! Send/receive buffer to/from bc_z%end/bc_z%end
                        CALL MPI_SENDRECV(                              &
                               q_cons_buff_send(0),                     &
                               buff_size*sys_size*(m+2*buff_size+1)  &
                                                    *(n+2*buff_size+1), &
                               MPI_DOUBLE_PRECISION, bc_z%end, 0,       &
                               q_cons_buff_recv(0),                     &
                               buff_size*sys_size*(m+2*buff_size+1)  &
                                                    *(n+2*buff_size+1), &
                               MPI_DOUBLE_PRECISION, bc_z%end, 1,       &
                               MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr  )
                        
                    END IF
                    
                    ! Unpacking buffer received from bc_z%end
                    DO l = p+1, p+buff_size
                        DO k = -buff_size, n+buff_size
                            DO j = -buff_size, m+buff_size
                                DO i = 1, sys_size
                                    r = (i-1) + sys_size * &
                                        ((j+buff_size) + (m+2*buff_size+1) * &
                                        ((k+buff_size) + (n+2*buff_size+1) * &
                                        (l-p-1)))
                                    q_cons_vf(i)%sf(j,k,l) = q_cons_buff_recv(r)
                                END DO
                            END DO
                        END DO
                    END DO
                    
                END IF
                
            END IF
            ! END: MPI Communication in z-direction ============================
            
            
        END SUBROUTINE s_mpi_sendrecv_conservative_variables_buffers ! ---------
        
        
        
        
        
        !> Module deallocation and/or disassociation procedures
        SUBROUTINE s_finalize_mpi_proxy_module() ! -----------------------------

            ! Deallocating q_cons_buff_send and q_cons_buff_recv
            DEALLOCATE( q_cons_buff_send, q_cons_buff_recv  )
            
        END SUBROUTINE s_finalize_mpi_proxy_module ! ---------------------------
        
        
        
        
        
        !> The subroutine finalizes the MPI execution environment.
        SUBROUTINE s_mpi_finalize() ! ------------------------------------------

            ! Finalizing the MPI environment
            CALL MPI_FINALIZE(ierr)
            
        END SUBROUTINE s_mpi_finalize ! ----------------------------------------
        
        
        
        
        
END MODULE m_mpi_proxy
