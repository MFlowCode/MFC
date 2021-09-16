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
!! @file m_mpi_proxy.f90
!! @brief Contains module m_mpi_proxy
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief  This module serves as a proxy to the parameters and subroutines
!!              available in the MPI implementation's MPI module. Specifically,
!!              the role of the proxy is to harness basic MPI commands into more
!!              complex procedures as to achieve the required communication goals
!!              for the post-process.
MODULE m_mpi_proxy
    
    
    ! Dependencies =============================================================
    USE mpi                     !< Message passing interface (MPI) module
    
    USE m_derived_types         !< Definitions of the derived types
    
    USE m_global_parameters     !< Global parameters for the code
    ! ==========================================================================
    
    
    IMPLICIT NONE
    
    
    !> @name Buffers of the conservative variables recieved/sent from/to neighbooring
    !! processors. Note that these variables are structured as vectors rather
    !! than arrays.
    !> @{
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:) :: q_cons_buffer_in
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:) :: q_cons_buffer_out
    !> @}

    !> @name Recieve counts and displacement vector variables, respectively, used in
    !! enabling MPI to gather varying amounts of data from all processes to the
    !! root process
    !> @{
    INTEGER, ALLOCATABLE, DIMENSION(:) :: recvcounts
    INTEGER, ALLOCATABLE, DIMENSION(:) :: displs
    !> @}

    !> @name Generic flags used to identify and report MPI errors
    !> @{
    INTEGER, PRIVATE :: err_code, ierr
    !> @}
    
    CONTAINS
        
        
        
        
        !>  The subroutine intializes the MPI environment and queries
        !!      both the number of processors that will be available for
        !!      the job as well as the local processor rank.        
        SUBROUTINE s_mpi_initialize() ! ----------------------------
            
            
            ! Establishing the MPI environment
            CALL MPI_INIT(ierr)
            
            
            ! Checking whether the MPI environment has been properly intialized
            IF(ierr /= MPI_SUCCESS) THEN
                PRINT '(A)', 'Unable to initialize MPI environment. Exiting ...'
                CALL MPI_ABORT(MPI_COMM_WORLD, err_code, ierr)
            END IF
            
            
            ! Querying number of processors available for the job
            CALL MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
            
            
            ! Identifying the rank of the local processor
            CALL MPI_COMM_RANK(MPI_COMM_WORLD, proc_rank ,ierr)
            
            
        END SUBROUTINE s_mpi_initialize ! --------------------------
        
        
        
        
        
        !> The subroutine terminates the MPI execution environment.
        SUBROUTINE s_mpi_abort() ! ---------------------------------------------
            
            ! Terminating the MPI environment
            CALL MPI_ABORT(MPI_COMM_WORLD, err_code, ierr)
            
        END SUBROUTINE s_mpi_abort ! -------------------------------------------
        
        
        
        
        !> This subroutine defines local and global sizes for the data
        !> @name q_cons_vf Conservative variables
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





        !>Halts all processes until all have reached barrier.
        SUBROUTINE s_mpi_barrier() ! -------------------------------------------

            ! Calling MPI_BARRIER
            CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

        END SUBROUTINE s_mpi_barrier ! -----------------------------------------



        !>  Computation of parameters, allocation procedures, and/or
        !!      any other tasks needed to properly setup the module
        SUBROUTINE s_initialize_mpi_proxy_module() ! ------------------------------
           
            INTEGER :: i !< Generic loop iterator
            
            
            ! Allocating vectorized buffer regions of conservative variables.
            ! The length of buffer vectors are set according to the size of the
            ! largest buffer region in the sub-domain.
            IF(buff_size > 0) THEN
                
                ! Simulation is at least 2D
                IF(n > 0) THEN
                    
                    ! Simulation is 3D
                    IF(p > 0) THEN
                        
                        ALLOCATE( q_cons_buffer_in( 0 : buff_size *         &
                                                        sys_size *          &
                                                        (m+2*buff_size+1) * &
                                                        (n+2*buff_size+1) * &
                                                        (p+2*buff_size+1) / &
                                                        ( MIN(m,n,p)        &
                                                        + 2*buff_size+1 )-1 ))
                        ALLOCATE(q_cons_buffer_out( 0 : buff_size *         &
                                                        sys_size *          &
                                                        (m+2*buff_size+1) * &
                                                        (n+2*buff_size+1) * &
                                                        (p+2*buff_size+1) / &
                                                        ( MIN(m,n,p)        &
                                                        + 2*buff_size+1 )-1 ))
                        
                    ! Simulation is 2D
                    ELSE
                        
                        ALLOCATE( q_cons_buffer_in( 0 : buff_size *         &
                                                        sys_size *          &
                                                        ( MAX(m,n)          &
                                                        + 2*buff_size+1 )-1 ))
                        ALLOCATE(q_cons_buffer_out( 0 : buff_size *         &
                                                        sys_size *          &
                                                        ( MAX(m,n)          &
                                                        + 2*buff_size+1 )-1 ))
                        
                    END IF
                    
                ! Simulation is 1D
                ELSE
                    
                    ALLOCATE( q_cons_buffer_in(0:buff_size*sys_size-1))
                    ALLOCATE(q_cons_buffer_out(0:buff_size*sys_size-1))
                    
                END IF
                
                ! Initially zeroing out the vectorized buffer region variables
                ! to avoid possible underflow from any unused allocated memory
                q_cons_buffer_in  = 0d0
                q_cons_buffer_out = 0d0
                
            END IF
            
            
            ! Allocating and configuring the recieve counts and the displacement
            ! vector variables used in variable-gather communication procedures.
            ! Note that these are only needed for either multidimensional runs
            ! that utilize the Silo database file format or for 1D simulations.
            IF((format == 1 .AND. n > 0) .OR. n == 0) THEN
                
                ALLOCATE(recvcounts(0:num_procs-1))
                ALLOCATE(    displs(0:num_procs-1))
                
                IF(n == 0) THEN
                    CALL MPI_GATHER( m+1, 1, MPI_INTEGER, recvcounts(0), 1, &
                                     MPI_INTEGER, 0, MPI_COMM_WORLD, ierr   )
                ELSEIF(proc_rank == 0) THEN
                    recvcounts = 1
                END IF
                
                IF(proc_rank == 0) THEN
                    displs(0) = 0
                    
                    DO i = 1, num_procs-1
                        displs(i) = displs(i-1) + recvcounts(i-1)
                    END DO
                END IF
                
            END IF
            
            
        END SUBROUTINE s_initialize_mpi_proxy_module ! ----------------------------
        
        
        
        
        !>  Since only processor with rank 0 is in charge of reading
        !!      and checking the consistency of the user provided inputs,
        !!      these are not available to the remaining processors. This
        !!      subroutine is then in charge of broadcasting the required
        !!      information.        
        SUBROUTINE s_mpi_bcast_user_inputs() ! ---------------------------------


            INTEGER :: i !< Generic loop iterator
            
            
            ! Logistics
            CALL MPI_BCAST( case_dir, LEN(case_dir), MPI_CHARACTER ,      &
                                            0      , MPI_COMM_WORLD, ierr )
            
            
            ! Computational domain parameters
            CALL MPI_BCAST(m, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(p, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(m_glb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(n_glb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(p_glb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            
            CALL MPI_BCAST(cyl_coord, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            
            CALL MPI_BCAST(t_step_start, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST(t_step_stop , 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST(t_step_save , 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
            
            
            ! Simulation algorithm parameters
            CALL MPI_BCAST(model_eqns, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(num_fluids, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(adv_alphan, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(mpp_lim, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(weno_order, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(mixture_err, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(alt_soundspeed, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            
            CALL MPI_BCAST(bc_x%beg, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(bc_x%end, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(bc_y%beg, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(bc_y%end, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(bc_z%beg, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(bc_z%end, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            
            CALL MPI_BCAST(parallel_io, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            
            ! Fluids physical parameters
            DO i = 1, num_fluids_max
                CALL MPI_BCAST( fluid_pp(i)%gamma   , 1, &
                                MPI_DOUBLE_PRECISION, 0, &
                                MPI_COMM_WORLD, ierr     )
                CALL MPI_BCAST( fluid_pp(i)%pi_inf  , 1, &
                                MPI_DOUBLE_PRECISION, 0, &
                                MPI_COMM_WORLD, ierr     )
            END DO
            
            
            ! Formatted database file(s) structure parameters
            CALL MPI_BCAST(format, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            
            CALL MPI_BCAST(precision, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
            
            CALL MPI_BCAST( coarsen_silo  ,    1    , MPI_LOGICAL,         &
                                               0    , MPI_COMM_WORLD, ierr )

            CALL MPI_BCAST( rho_wrt       ,    1    , MPI_LOGICAL,         &
                                               0    , MPI_COMM_WORLD, ierr )
            CALL MPI_BCAST( mom_wrt(1)    ,    3    , MPI_LOGICAL,         &
                                               0    , MPI_COMM_WORLD, ierr )
            CALL MPI_BCAST( vel_wrt(1)    ,    3    , MPI_LOGICAL,         &
                                               0    , MPI_COMM_WORLD, ierr )
            CALL MPI_BCAST(flux_lim      , 1, MPI_INTEGER         , &
                                           0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(flux_wrt(1)    ,    3    , MPI_LOGICAL,         &
                                               0    , MPI_COMM_WORLD, ierr )
            CALL MPI_BCAST( kappa_wrt(1)    , num_fluids_max, MPI_LOGICAL   , &
                                                    0       , MPI_COMM_WORLD, &
                                                              ierr            )
            CALL MPI_BCAST( E_wrt         ,    1    , MPI_LOGICAL,         &
                                               0    , MPI_COMM_WORLD, ierr )
            CALL MPI_BCAST( pres_wrt      ,    1    , MPI_LOGICAL,         &
                                               0    , MPI_COMM_WORLD, ierr )
            CALL MPI_BCAST( gamma_wrt     ,    1    , MPI_LOGICAL,         &
                                               0    , MPI_COMM_WORLD, ierr )
            CALL MPI_BCAST( heat_ratio_wrt,    1    , MPI_LOGICAL,         &
                                               0    , MPI_COMM_WORLD, ierr )
            CALL MPI_BCAST( pi_inf_wrt    ,    1    , MPI_LOGICAL,         &
                                               0    , MPI_COMM_WORLD, ierr )
            CALL MPI_BCAST( pres_inf_wrt  ,    1    , MPI_LOGICAL,         &
                                               0    , MPI_COMM_WORLD, ierr )
            CALL MPI_BCAST( cons_vars_wrt ,    1    , MPI_LOGICAL,         &
                                               0    , MPI_COMM_WORLD, ierr )
            CALL MPI_BCAST( prim_vars_wrt ,    1    , MPI_LOGICAL,         &
                                               0    , MPI_COMM_WORLD, ierr )
            CALL MPI_BCAST( c_wrt         ,    1    , MPI_LOGICAL,         &
                                               0    , MPI_COMM_WORLD, ierr )
            CALL MPI_BCAST( omega_wrt(1)  ,    3    , MPI_LOGICAL,         &
                                               0    , MPI_COMM_WORLD, ierr )
            CALL MPI_BCAST( schlieren_wrt ,    1    , MPI_LOGICAL,         &
                                               0    , MPI_COMM_WORLD, ierr )
            CALL MPI_BCAST( alpha_rho_wrt(1), num_fluids_max, MPI_LOGICAL   , &
                                                    0       , MPI_COMM_WORLD, &
                                                              ierr            )
            CALL MPI_BCAST( alpha_wrt(1)    , num_fluids_max, MPI_LOGICAL   , &
                                                    0       , MPI_COMM_WORLD, &
                                                              ierr            )
            
            CALL MPI_BCAST( schlieren_alpha(1)  , num_fluids_max, &
                            MPI_DOUBLE_PRECISION,        0      , &
                            MPI_COMM_WORLD, ierr                  )
            
            CALL MPI_BCAST(fd_order, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            
            CALL MPI_BCAST( fourier_decomp,    1    , MPI_LOGICAL,         &
                                               0    , MPI_COMM_WORLD, ierr )
            CALL MPI_BCAST(fourier_modes%beg, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(fourier_modes%end, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

            ! Tait EOS
            CALL MPI_BCAST( pref,1,             &
                        MPI_DOUBLE_PRECISION,0, &
                        MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST( rhoref,1,           &
                        MPI_DOUBLE_PRECISION,0, &
                        MPI_COMM_WORLD,ierr)

            ! Bubble modeling
            CALL MPI_BCAST( bubbles,1,          &
                        MPI_LOGICAL,0,          &
                        MPI_COMM_WORLD,ierr  )
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
            CALL MPI_BCAST( polydisperse,1,          &
                        MPI_LOGICAL,0,          &
                        MPI_COMM_WORLD,ierr  )
            CALL MPI_BCAST( poly_sigma,1,            &
                        MPI_DOUBLE_PRECISION,0, &
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
        END SUBROUTINE s_mpi_bcast_user_inputs ! -------------------------------
        
        
        
        
        !>  This subroutine takes care of efficiently distributing
        !!      the computational domain among the available processors
        !!      as well as recomputing some of the global parameters so
        !!      that they reflect the configuration of sub-domain that
        !!      is overseen by the local processor.        
        SUBROUTINE s_mpi_decompose_computational_domain() ! --------------------

            
            ! # of processors in the x-, y- and z-coordinate directions
            INTEGER :: num_procs_x, num_procs_y, num_procs_z
            
            ! Temporary # of processors in x-, y- and z-coordinate directions
            ! used during the processor factorization optimization procedure
            REAL(KIND(0d0)) :: tmp_num_procs_x, tmp_num_procs_y, tmp_num_procs_z
            
            ! Processor factorization (fct) minimization parameter
            REAL(KIND(0d0)) :: fct_min
            
            ! Cartesian processor topology communicator
            INTEGER :: MPI_COMM_CART
            
            ! Number of remaining cells for a particular coordinate direction
            ! after the bulk has evenly been distributed among the available
            ! processors for that coordinate direction
            INTEGER :: rem_cells
            
            ! Generic loop iterators
            INTEGER :: i,j
            
            IF (num_procs == 1 .AND. parallel_io) THEN
                DO i = 1, num_dims
                    start_idx(i) = 0
                END DO
                RETURN
            END IF
            
            ! Performing the computational domain decomposition. The procedure
            ! is optimized by ensuring that each processor contains a close to
            ! equivalent piece of the computational domain. Note that explicit
            ! type-casting is omitted here for code legibility purposes.
            
            ! Generating 3D Cartesian Processor Topology =======================
            
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

                        ! Initial values of the processor factorization optimization
                        num_procs_x =  1
                        num_procs_y =  1
                        num_procs_z = num_procs
                        ierr        = -1
                        
                        ! Computing minimization variable for these initial values
                        tmp_num_procs_x = num_procs_x
                        tmp_num_procs_y = num_procs_y
                        tmp_num_procs_z = num_procs_z
                        fct_min         = 10d0*ABS((m+1)/tmp_num_procs_x  &
                                                  -(n+1)/tmp_num_procs_y) &
                                        + 10d0*ABS((n+1)/tmp_num_procs_y  &
                                                  -(p+1)/tmp_num_procs_z)
                        
                        ! Searching for optimal computational domain distribution
                        DO i = 1, num_procs
                            
                            IF(        MOD(num_procs,i) == 0        &
                                               .AND.                &
                                (m+1)/i >= num_stcls_min*weno_order ) THEN
                                
                                DO j = 1, (num_procs/i)
                                    
                                    IF(       MOD(num_procs/i,j) == 0       &
                                                       .AND.                &
                                        (n+1)/j >= num_stcls_min*weno_order ) THEN
                                        
                                        tmp_num_procs_x = i
                                        tmp_num_procs_y = j
                                        tmp_num_procs_z = num_procs / (i*j)
                                        
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
                    
                    ! Checking whether the decomposition of the computational
                    ! domain was successful
                    IF(proc_rank == 0 .AND. ierr == -1) THEN
                        PRINT '(A)', 'Unable to decompose computational ' // &
                                     'domain for selected number of ' // &
                                     'processors. Exiting ...'
                        CALL MPI_ABORT(MPI_COMM_WORLD, err_code, ierr)
                    END IF
                    
                    ! Creating a new communicator using Cartesian topology
                    CALL MPI_CART_CREATE( MPI_COMM_WORLD, 3, (/ num_procs_x, &
                                          num_procs_y, num_procs_z /),       &
                                          (/ .TRUE., .TRUE., .TRUE. /),      &
                                          .FALSE., MPI_COMM_CART, ierr       )
                    
                    ! Finding corresponding Cartesian coordinates of the local
                    ! processor rank in newly declared cartesian communicator
                    CALL MPI_CART_COORDS( MPI_COMM_CART, proc_rank, 3, &
                                          proc_coords, ierr            )
                    
            ! END: Generating 3D Cartesian Processor Topology ==================
                    
                    
            ! Sub-domain Global Parameters in z-direction ======================
                    
                    ! Number of remaining cells after majority is distributed
                    rem_cells = MOD(p+1, num_procs_z)
                    
                    ! Optimal number of cells per processor
                    p = (p+1) / num_procs_z - 1
                    
                    ! Distributing any remaining cells
                    DO i = 1, rem_cells
                        IF(proc_coords(3) == i-1) THEN
                            p = p + 1
                            EXIT
                        END IF
                    END DO
                    
                    ! Boundary condition at the beginning
                    IF(proc_coords(3) > 0 .OR. bc_z%beg == -1) THEN
                        proc_coords(3) = proc_coords(3) - 1
                        CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, &
                                            bc_z%beg, ierr              )
                        proc_coords(3) = proc_coords(3) + 1
                    END IF
                    
                    ! Ghost zone at the beginning
                    IF(proc_coords(3) > 0 .AND. format == 1) THEN
                        offset_z%beg = 2
                    ELSE
                        offset_z%beg = 0
                    END IF
                    
                    ! Boundary condition at the end
                    IF(proc_coords(3) < num_procs_z-1 .OR. bc_z%end == -1) THEN
                        proc_coords(3) = proc_coords(3) + 1
                        CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, &
                                            bc_z%end, ierr              )
                        proc_coords(3) = proc_coords(3) - 1
                    END IF
                    
                    ! Ghost zone at the end
                    IF(proc_coords(3) < num_procs_z-1 .AND. format == 1) THEN
                        offset_z%end = 2
                    ELSE
                        offset_z%end = 0
                    END IF
                    
                    IF (parallel_io) THEN
                        IF (proc_coords(3) < rem_cells) THEN
                            start_idx(3) = (p+1) * proc_coords(3)
                        ELSE
                            start_idx(3) = (p+1) * proc_coords(3) + rem_cells
                        END IF
                    END IF
            ! ==================================================================
                    
                    
            ! Generating 2D Cartesian Processor Topology =======================
                    
                ELSE
                    
                    ! Initial values of the processor factorization optimization
                    num_procs_x =  1
                    num_procs_y = num_procs
                    ierr        = -1
                    
                    ! Computing minimization variable for these initial values
                    tmp_num_procs_x = num_procs_x
                    tmp_num_procs_y = num_procs_y
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
                    
                    ! Checking whether the decomposition of the computational
                    ! domain was successful
                    IF(proc_rank == 0 .AND. ierr == -1) THEN
                        PRINT '(A)', 'Unable to decompose computational ' // &
                                     'domain for selected number of ' // &
                                     'processors. Exiting ...'
                        CALL MPI_ABORT(MPI_COMM_WORLD, err_code, ierr)
                    END IF
                    
                    ! Creating a new communicator using Cartesian topology
                    CALL MPI_CART_CREATE( MPI_COMM_WORLD, 2, (/ num_procs_x, &
                                          num_procs_y /), (/ .TRUE.,         &
                                          .TRUE. /), .FALSE., MPI_COMM_CART, &
                                          ierr                               )
                    
                    ! Finding corresponding Cartesian coordinates of the local
                    ! processor rank in newly declared cartesian communicator
                    CALL MPI_CART_COORDS( MPI_COMM_CART, proc_rank, 2, &
                                          proc_coords, ierr            )
                    
                END IF
                
            ! END: Generating 2D Cartesian Processor Topology ==================
                
                
            ! Sub-domain Global Parameters in y-direction ======================
                
                ! Number of remaining cells after majority has been distributed
                rem_cells = MOD(n+1, num_procs_y)
                
                ! Optimal number of cells per processor
                n = (n+1) / num_procs_y - 1
                
                ! Distributing any remaining cells
                DO i = 1, rem_cells
                    IF(proc_coords(2) == i-1) THEN
                        n = n + 1
                        EXIT
                    END IF
                END DO
                
                ! Boundary condition at the beginning
                IF(proc_coords(2) > 0 .OR. bc_y%beg == -1) THEN
                    proc_coords(2) = proc_coords(2) - 1
                    CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, bc_y%beg, &
                                        ierr                                  )
                    proc_coords(2) = proc_coords(2) + 1
                END IF
                
                ! Ghost zone at the beginning
                IF(proc_coords(2) > 0 .AND. format == 1) THEN
                    offset_y%beg = 2
                ELSE
                    offset_y%beg = 0
                END IF
                
                ! Boundary condition at the end
                IF(proc_coords(2) < num_procs_y-1 .OR. bc_y%end == -1) THEN
                    proc_coords(2) = proc_coords(2) + 1
                    CALL MPI_CART_RANK( MPI_COMM_CART, proc_coords, bc_y%end, &
                                        ierr                                  )
                    proc_coords(2) = proc_coords(2) - 1
                END IF
                
                ! Ghost zone at the end
                IF(proc_coords(2) < num_procs_y-1 .AND. format == 1) THEN
                    offset_y%end = 2
                ELSE
                    offset_y%end = 0
                END IF
                
                IF (parallel_io) THEN
                    IF (proc_coords(2) < rem_cells) THEN
                        start_idx(2) = (n+1) * proc_coords(2)
                    ELSE
                        start_idx(2) = (n+1) * proc_coords(2) + rem_cells
                    END IF
                END IF
            ! ==================================================================
                
                
            ! Generating 1D Cartesian Processor Topology =======================
                
            ELSE
                
                ! Number of processors in the coordinate direction is equal to
                ! the total number of processors available
                num_procs_x = num_procs
                
                ! Number of cells in undecomposed computational domain needed
                ! for sub-domain reassembly during formatted data output
                m_root = m
                
                ! Creating a new communicator using Cartesian topology
                CALL MPI_CART_CREATE( MPI_COMM_WORLD, 1, (/ num_procs_x /), &
                                      (/ .TRUE. /), .FALSE., MPI_COMM_CART, &
                                      ierr                                  )
                
                ! Finding the corresponding Cartesian coordinates of the local
                ! processor rank in the newly declared cartesian communicator
                CALL MPI_CART_COORDS( MPI_COMM_CART, proc_rank, 1, &
                                      proc_coords, ierr            )
                
            END IF
            
            ! ==================================================================
            
            
            ! Sub-domain Global Parameters in x-direction ======================
            
            ! Number of remaining cells after majority has been distributed
            rem_cells = MOD(m+1, num_procs_x)
            
            ! Optimal number of cells per processor
            m = (m+1) / num_procs_x - 1
            
            ! Distributing any remaining cells
            DO i = 1, rem_cells
                IF(proc_coords(1) == i-1) THEN
                    m = m + 1
                    EXIT
                END IF
            END DO
            
            ! Boundary condition at the beginning
            IF(proc_coords(1) > 0 .OR. bc_x%beg == -1) THEN
                proc_coords(1) = proc_coords(1) - 1
                CALL MPI_CART_RANK(MPI_COMM_CART, proc_coords, bc_x%beg, ierr)
                proc_coords(1) = proc_coords(1) + 1
            END IF
            
            ! Ghost zone at the beginning
            IF(proc_coords(1) > 0 .AND. format == 1 .AND. n > 0) THEN
                offset_x%beg = 2
            ELSE
                offset_x%beg = 0
            END IF
            
            ! Boundary condition at the end
            IF(proc_coords(1) < num_procs_x-1 .OR. bc_x%end == -1) THEN
                proc_coords(1) = proc_coords(1) + 1
                CALL MPI_CART_RANK(MPI_COMM_CART, proc_coords, bc_x%end, ierr)
                proc_coords(1) = proc_coords(1) - 1
            END IF
            
            ! Ghost zone at the end
            IF(proc_coords(1)< num_procs_x-1 .AND. format == 1 .AND. n > 0) THEN
                offset_x%end = 2
            ELSE
                offset_x%end = 0
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
        
        
        
        
        !>  Communicates the buffer regions associated with the grid
        !!      variables with processors in charge of the neighbooring
        !!      sub-domains. Note that only cell-width spacings feature
        !!      buffer regions so that no information relating to the
        !!      cell-boundary locations is communicated.
        !!  @param pbc_loc Processor boundary condition (PBC) location
        !!  @param sweep_coord Coordinate direction normal to the processor boundary
        SUBROUTINE s_mpi_sendrecv_grid_vars_buffer_regions(pbc_loc, sweep_coord)

            CHARACTER(LEN = 3), INTENT(IN) :: pbc_loc
            CHARACTER, INTENT(IN) :: sweep_coord
            
            ! Communications in the x-direction ================================
            
            IF(sweep_coord == 'x') THEN
                
                IF(pbc_loc == 'beg') THEN    ! Buffer region at the beginning
                    
                    ! PBC at both ends of the sub-domain
                    IF(bc_x%end >= 0) THEN
                        
                        ! Sending/receiving the data to/from bc_x%end/bc_x%beg
                        CALL MPI_SENDRECV( dx(m-buff_size+1), buff_size,      &
                                           MPI_DOUBLE_PRECISION, bc_x%end, 0, &
                                           dx(-buff_size), buff_size,         &
                                           MPI_DOUBLE_PRECISION, bc_x%beg, 0, &
                                           MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                           ierr                               )
                        
                    ! PBC only at beginning of the sub-domain
                    ELSE
                        
                        ! Sending/receiving the data to/from bc_x%beg/bc_x%beg
                        CALL MPI_SENDRECV( dx(0), buff_size,                  &
                                           MPI_DOUBLE_PRECISION, bc_x%beg, 1, &
                                           dx(-buff_size), buff_size,         &
                                           MPI_DOUBLE_PRECISION, bc_x%beg, 0, &
                                           MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                           ierr                               )
                        
                    END IF
                    
                ELSE                         ! Buffer region at the end
                    
                    ! PBC at both ends of the sub-domain
                    IF(bc_x%beg >= 0) THEN
                        
                        ! Sending/receiving the data to/from bc_x%beg/bc_x%end
                        CALL MPI_SENDRECV( dx(0), buff_size,                  &
                                           MPI_DOUBLE_PRECISION, bc_x%beg, 1, &
                                           dx(m+1), buff_size,                &
                                           MPI_DOUBLE_PRECISION, bc_x%end, 1, &
                                           MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                           ierr                               )
                        
                    ! PBC only at end of the sub-domain
                    ELSE
                        
                        ! Sending/receiving the data to/from bc_x%end/bc_x%end
                        CALL MPI_SENDRECV( dx(m-buff_size+1), buff_size,      &
                                           MPI_DOUBLE_PRECISION, bc_x%end, 0, &
                                           dx(m+1), buff_size,                &
                                           MPI_DOUBLE_PRECISION, bc_x%end, 1, &
                                           MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                           ierr                               )
                        
                    END IF
                    
                END IF
                
            ! END: Communications in the x-direction ===========================
                
                
            ! Communications in the y-direction ================================
                
            ELSEIF(sweep_coord == 'y') THEN
                
                IF(pbc_loc == 'beg') THEN    ! Buffer region at the beginning
                    
                    ! PBC at both ends of the sub-domain
                    IF(bc_y%end >= 0) THEN
                        
                        ! Sending/receiving the data to/from bc_y%end/bc_y%beg
                        CALL MPI_SENDRECV( dy(n-buff_size+1), buff_size,      &
                                           MPI_DOUBLE_PRECISION, bc_y%end, 0, &
                                           dy(-buff_size), buff_size,         &
                                           MPI_DOUBLE_PRECISION, bc_y%beg, 0, &
                                           MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                           ierr                               )
                        
                    ! PBC only at beginning of the sub-domain
                    ELSE
                        
                        ! Sending/receiving the data to/from bc_y%beg/bc_y%beg
                        CALL MPI_SENDRECV( dy(0), buff_size,                  &
                                           MPI_DOUBLE_PRECISION, bc_y%beg, 1, &
                                           dy(-buff_size), buff_size,         &
                                           MPI_DOUBLE_PRECISION, bc_y%beg, 0, &
                                           MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                           ierr                               )
                        
                    END IF
                    
                ELSE                         ! Buffer region at the end
                    
                    ! PBC at both ends of the sub-domain
                    IF(bc_y%beg >= 0) THEN
                        
                        ! Sending/receiving the data to/from bc_y%beg/bc_y%end
                        CALL MPI_SENDRECV( dy(0), buff_size,                  &
                                           MPI_DOUBLE_PRECISION, bc_y%beg, 1, &
                                           dy(n+1), buff_size,                &
                                           MPI_DOUBLE_PRECISION, bc_y%end, 1, &
                                           MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                           ierr                               )
                        
                    ! PBC only at end of the sub-domain
                    ELSE
                        
                        ! Sending/receiving the data to/from bc_y%end/bc_y%end
                        CALL MPI_SENDRECV( dy(n-buff_size+1), buff_size,      &
                                           MPI_DOUBLE_PRECISION, bc_y%end, 0, &
                                           dy(n+1), buff_size,                &
                                           MPI_DOUBLE_PRECISION, bc_y%end, 1, &
                                           MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                           ierr                               )
                        
                    END IF
                    
                END IF
                
            ! END: Communications in the y-direction ===========================
                
                
            ! Communications in the z-direction ================================
                
            ELSE
                
                IF(pbc_loc == 'beg') THEN    ! Buffer region at the beginning
                    
                    ! PBC at both ends of the sub-domain
                    IF(bc_z%end >= 0) THEN
                        
                        ! Sending/receiving the data to/from bc_z%end/bc_z%beg
                        CALL MPI_SENDRECV( dz(p-buff_size+1), buff_size,      &
                                           MPI_DOUBLE_PRECISION, bc_z%end, 0, &
                                           dz(-buff_size), buff_size,         &
                                           MPI_DOUBLE_PRECISION, bc_z%beg, 0, &
                                           MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                           ierr                               )
                        
                    ! PBC only at beginning of the sub-domain
                    ELSE
                        
                        ! Sending/receiving the data to/from bc_z%beg/bc_z%beg
                        CALL MPI_SENDRECV( dz(0), buff_size,                  &
                                           MPI_DOUBLE_PRECISION, bc_z%beg, 1, &
                                           dz(-buff_size), buff_size,         &
                                           MPI_DOUBLE_PRECISION, bc_z%beg, 0, &
                                           MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                           ierr                               )
                        
                    END IF
                    
                ELSE                         ! Buffer region at the end
                    
                    ! PBC at both ends of the sub-domain
                    IF(bc_z%beg >= 0) THEN
                        
                        ! Sending/receiving the data to/from bc_z%beg/bc_z%end
                        CALL MPI_SENDRECV( dz(0), buff_size,                  &
                                           MPI_DOUBLE_PRECISION, bc_z%beg, 1, &
                                           dz(p+1), buff_size,                &
                                           MPI_DOUBLE_PRECISION, bc_z%end, 1, &
                                           MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                           ierr                               )
                        
                    ! PBC only at end of the sub-domain
                    ELSE
                        
                        ! Sending/receiving the data to/from bc_z%end/bc_z%end
                        CALL MPI_SENDRECV( dz(p-buff_size+1), buff_size,      &
                                           MPI_DOUBLE_PRECISION, bc_z%end, 0, &
                                           dz(p+1), buff_size,                &
                                           MPI_DOUBLE_PRECISION, bc_z%end, 1, &
                                           MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                           ierr                               )
                        
                    END IF
                    
                END IF
                
            END IF
            
            ! END: Communications in the z-direction ===========================
            
            
        END SUBROUTINE s_mpi_sendrecv_grid_vars_buffer_regions ! ---------------
        
        
        
        
        !>  Communicates buffer regions associated with conservative
        !!      variables with processors in charge of the neighbooring
        !!      sub-domains
        !!  @param q_cons_vf Conservative variables
        !!  @param pbc_loc Processor boundary condition (PBC) location
        !!  @param sweep_coord Coordinate direction normal to the processor boundary
        SUBROUTINE s_mpi_sendrecv_cons_vars_buffer_regions(q_cons_vf, pbc_loc, &
                                                               sweep_coord     )

            TYPE(scalar_field), &
            DIMENSION(sys_size), &
            INTENT(INOUT) :: q_cons_vf
            
            CHARACTER(LEN = 3), INTENT(IN) :: pbc_loc
            
            CHARACTER, INTENT(IN) :: sweep_coord
            
            INTEGER :: i,j,k,l,r !< Generic loop iterators
            
            
            ! Communications in the x-direction ================================
            
            IF(sweep_coord == 'x') THEN
                
                IF(pbc_loc == 'beg') THEN    ! Buffer region at the beginning
                    
                    ! PBC at both ends of the sub-domain
                    IF(bc_x%end >= 0) THEN
                        
                        ! Packing the data to be sent to bc_x%end
                        DO l = 0, p
                            DO k = 0, n
                                DO j = m-buff_size+1, m
                                    DO i = 1, sys_size
                                        r = sys_size*(j-m+buff_size-1)   &
                                          + sys_size*buff_size*k + (i-1) &
                                          + sys_size*buff_size*(n+1)*l
                                        q_cons_buffer_out(r) = &
                                                q_cons_vf(i)%sf(j,k,l)
                                    END DO
                                END DO
                            END DO
                        END DO
                        
                        ! Sending/receiving the data to/from bc_x%end/bc_x%beg
                        CALL MPI_SENDRECV(q_cons_buffer_out(0),              &
                                          buff_size*sys_size*(n+1)*(p+1),    &
                                          MPI_DOUBLE_PRECISION, bc_x%end, 0, &
                                          q_cons_buffer_in(0),               &
                                          buff_size*sys_size*(n+1)*(p+1),    &
                                          MPI_DOUBLE_PRECISION, bc_x%beg, 0, &
                                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                          ierr                               )
                        
                    ! PBC only at beginning of the sub-domain
                    ELSE
                        
                        ! Packing the data to be sent to bc_x%beg
                        DO l = 0, p
                            DO k = 0, n
                                DO j = 0, buff_size-1
                                    DO i = 1, sys_size
                                        r = (i-1) + sys_size*j         &
                                          + sys_size*buff_size*k       &
                                          + sys_size*buff_size*(n+1)*l
                                        q_cons_buffer_out(r) = &
                                                q_cons_vf(i)%sf(j,k,l)
                                    END DO
                                END DO
                            END DO
                        END DO
                        
                        ! Sending/receiving the data to/from bc_x%beg/bc_x%beg
                        CALL MPI_SENDRECV(q_cons_buffer_out(0),              &
                                          buff_size*sys_size*(n+1)*(p+1),    &
                                          MPI_DOUBLE_PRECISION, bc_x%beg, 1, &
                                          q_cons_buffer_in(0),               &
                                          buff_size*sys_size*(n+1)*(p+1),    &
                                          MPI_DOUBLE_PRECISION, bc_x%beg, 0, &
                                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                          ierr                               )
                        
                    END IF
                    
                    ! Unpacking the data received from bc_x%beg
                    DO l = 0, p
                        DO k = 0, n
                            DO j = -buff_size, -1
                                DO i = 1, sys_size
                                    r = sys_size*(j+buff_size)       &
                                      + sys_size*buff_size*k + (i-1) &
                                      + sys_size*buff_size*(n+1)*l
                                    q_cons_vf(i)%sf(j,k,l) = q_cons_buffer_in(r)
                                END DO
                            END DO
                        END DO
                    END DO
                    
                ELSE                         ! Buffer region at the end
                    
                    ! PBC at both ends of the sub-domain
                    IF(bc_x%beg >= 0) THEN
                        
                        ! Packing the data to be sent to bc_x%beg
                        DO l = 0, p
                            DO k = 0, n
                                DO j = 0, buff_size-1
                                    DO i = 1, sys_size
                                        r = (i-1) + sys_size*j         &
                                          + sys_size*buff_size*k       &
                                          + sys_size*buff_size*(n+1)*l
                                        q_cons_buffer_out(r) = &
                                                q_cons_vf(i)%sf(j,k,l)
                                    END DO
                                END DO
                            END DO
                        END DO
                        
                        ! Sending/receiving the data to/from bc_x%beg/bc_x%end
                        CALL MPI_SENDRECV(q_cons_buffer_out(0),              &
                                          buff_size*sys_size*(n+1)*(p+1),    &
                                          MPI_DOUBLE_PRECISION, bc_x%beg, 1, &
                                          q_cons_buffer_in(0),               &
                                          buff_size*sys_size*(n+1)*(p+1),    &
                                          MPI_DOUBLE_PRECISION, bc_x%end, 1, &
                                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                          ierr                               )
                        
                    ! PBC only at end of the sub-domain
                    ELSE
                        
                        ! Packing the data to be sent to bc_x%end
                        DO l = 0, p
                            DO k = 0, n
                                DO j = m-buff_size+1, m
                                    DO i = 1, sys_size
                                        r = sys_size*(j-m+buff_size-1)   &
                                          + sys_size*buff_size*k + (i-1) &
                                          + sys_size*buff_size*(n+1)*l
                                        q_cons_buffer_out(r) = &
                                                q_cons_vf(i)%sf(j,k,l)
                                    END DO
                                END DO
                            END DO
                        END DO
                        
                        ! Sending/receiving the data to/from bc_x%end/bc_x%end
                        CALL MPI_SENDRECV(q_cons_buffer_out(0),              &
                                          buff_size*sys_size*(n+1)*(p+1),    &
                                          MPI_DOUBLE_PRECISION, bc_x%end, 0, &
                                          q_cons_buffer_in(0),               &
                                          buff_size*sys_size*(n+1)*(p+1),    &
                                          MPI_DOUBLE_PRECISION, bc_x%end, 1, &
                                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                          ierr                               )
                        
                    END IF
                    
                    ! Unpacking the data received from bc_x%end
                    DO l = 0, p
                        DO k = 0, n
                            DO j = m+1, m+buff_size
                                DO i = 1, sys_size
                                    r = (i-1) + sys_size*(j-m-1)   &
                                      + sys_size*buff_size*k       &
                                      + sys_size*buff_size*(n+1)*l
                                    q_cons_vf(i)%sf(j,k,l) = q_cons_buffer_in(r)
                                END DO
                            END DO
                        END DO
                    END DO
                    
                END IF
                
            ! END: Communications in the x-direction ===========================
                
                
            ! Communications in the y-direction ================================
                
            ELSEIF(sweep_coord == 'y') THEN
                
                IF(pbc_loc == 'beg') THEN    ! Buffer region at the beginning
                    
                    ! PBC at both ends of the sub-domain
                    IF(bc_y%end >= 0) THEN
                        
                        ! Packing the data to be sent to bc_y%end
                        DO l = 0, p
                            DO k = n-buff_size+1, n
                                DO j = -buff_size, m+buff_size
                                    DO i = 1, sys_size
                                        r = sys_size*(j+buff_size)       &
                                          + sys_size*(m+2*buff_size+1) * &
                                            (k-n+buff_size-1) + (i-1)    &
                                          + sys_size*(m+2*buff_size+1) * &
                                            buff_size*l
                                        q_cons_buffer_out(r) = &
                                                   q_cons_vf(i)%sf(j,k,l)
                                    END DO
                                END DO
                            END DO
                        END DO
                        
                        ! Sending/receiving the data to/from bc_y%end/bc_y%beg
                        CALL MPI_SENDRECV( q_cons_buffer_out(0), buff_size *  &
                                           sys_size*(m+2*buff_size+1) *       &
                                           (p+1), MPI_DOUBLE_PRECISION,       &
                                           bc_y%end, 0, q_cons_buffer_in(0),  &
                                           buff_size*sys_size *               &
                                           (m+2*buff_size+1)*(p+1),           &
                                           MPI_DOUBLE_PRECISION, bc_y%beg, 0, &
                                           MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                           ierr                               )
                        
                    ! PBC only at beginning of the sub-domain
                    ELSE
                        
                        ! Packing the data to be sent to bc_y%beg
                        DO l = 0, p
                            DO k = 0, buff_size-1
                                DO j = -buff_size, m+buff_size
                                    DO i = 1, sys_size
                                        r = sys_size*(j+buff_size)       &
                                          + sys_size*(m+2*buff_size+1)*k &
                                          + sys_size*(m+2*buff_size+1) * &
                                            buff_size*l + (i-1)
                                        q_cons_buffer_out(r) = &
                                                   q_cons_vf(i)%sf(j,k,l)
                                    END DO
                                END DO
                            END DO
                        END DO
                        
                        ! Sending/receiving the data to/from bc_y%beg/bc_y%beg
                        CALL MPI_SENDRECV( q_cons_buffer_out(0), buff_size *  &
                                           sys_size*(m+2*buff_size+1) *       &
                                           (p+1), MPI_DOUBLE_PRECISION,       &
                                           bc_y%beg, 1, q_cons_buffer_in(0),  &
                                           buff_size*sys_size *               &
                                           (m+2*buff_size+1)*(p+1),           &
                                           MPI_DOUBLE_PRECISION, bc_y%beg, 0, &
                                           MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                           ierr                               )
                        
                    END IF
                    
                    ! Unpacking the data received from bc_y%beg
                    DO l = 0, p
                        DO k = -buff_size, -1
                            DO j = -buff_size, m+buff_size
                                DO i = 1, sys_size
                                    r = (i-1) + sys_size*(j+buff_size) &
                                      + sys_size*(m+2*buff_size+1) *   &
                                        (k+buff_size) + sys_size *     &
                                        (m+2*buff_size+1)*buff_size*l
                                    q_cons_vf(i)%sf(j,k,l) = q_cons_buffer_in(r)
                                END DO
                            END DO
                        END DO
                    END DO
                    
                ELSE                         ! Buffer region at the end
                    
                    ! PBC at both ends of the sub-domain
                    IF(bc_y%beg >= 0) THEN
                        
                        ! Packing the data to be sent to bc_y%beg
                        DO l = 0, p
                            DO k = 0, buff_size-1
                                DO j = -buff_size, m+buff_size
                                    DO i = 1, sys_size
                                        r = sys_size*(j+buff_size)       &
                                          + sys_size*(m+2*buff_size+1)*k &
                                          + sys_size*(m+2*buff_size+1) * &
                                            buff_size*l + (i-1)
                                        q_cons_buffer_out(r) = &
                                                   q_cons_vf(i)%sf(j,k,l)
                                    END DO
                                END DO
                            END DO
                        END DO
                        
                        ! Sending/receiving the data to/from bc_y%beg/bc_y%end
                        CALL MPI_SENDRECV( q_cons_buffer_out(0), buff_size *  &
                                           sys_size*(m+2*buff_size+1) *       &
                                           (p+1), MPI_DOUBLE_PRECISION,       &
                                           bc_y%beg, 1, q_cons_buffer_in(0),  &
                                           buff_size*sys_size *               &
                                           (m+2*buff_size+1)*(p+1),           &
                                           MPI_DOUBLE_PRECISION, bc_y%end, 1, &
                                           MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                           ierr                               )
                        
                    ! PBC only at end of the sub-domain
                    ELSE
                        
                        ! Packing the data to be sent to bc_y%end
                        DO l = 0, p
                            DO k = n-buff_size+1, n
                                DO j = -buff_size, m+buff_size
                                    DO i = 1, sys_size
                                        r = sys_size*(j+buff_size)       &
                                          + sys_size*(m+2*buff_size+1) * &
                                            (k-n+buff_size-1) + (i-1)    &
                                          + sys_size*(m+2*buff_size+1) * &
                                            buff_size*l
                                        q_cons_buffer_out(r) = &
                                                   q_cons_vf(i)%sf(j,k,l)
                                    END DO
                                END DO
                            END DO
                        END DO
                        
                        ! Sending/receiving the data to/from bc_y%end/bc_y%end
                        CALL MPI_SENDRECV( q_cons_buffer_out(0), buff_size *  &
                                           sys_size*(m+2*buff_size+1) *       &
                                           (p+1), MPI_DOUBLE_PRECISION,       &
                                           bc_y%end, 0, q_cons_buffer_in(0),  &
                                           buff_size*sys_size *               &
                                           (m+2*buff_size+1)*(p+1),           &
                                           MPI_DOUBLE_PRECISION, bc_y%end, 1, &
                                           MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                           ierr                               )
                        
                    END IF
                    
                    ! Unpacking the data received form bc_y%end
                    DO l = 0, p
                        DO k = n+1, n+buff_size
                            DO j = -buff_size, m+buff_size
                                DO i = 1, sys_size
                                    r = (i-1) + sys_size*(j+buff_size) &
                                      + sys_size*(m+2*buff_size+1) *   &
                                        (k-n-1) + sys_size *           &
                                        (m+2*buff_size+1)*buff_size*l
                                    q_cons_vf(i)%sf(j,k,l) = q_cons_buffer_in(r)
                                END DO
                            END DO
                        END DO
                    END DO
                    
                END IF
                
            ! END: Communications in the y-direction ===========================
                
                
            ! Communications in the z-direction ================================
                
            ELSE
                
                IF(pbc_loc == 'beg') THEN    ! Buffer region at the beginning
                    
                    ! PBC at both ends of the sub-domain
                    IF(bc_z%end >= 0) THEN
                        
                        ! Packing the data to be sent to bc_z%end
                        DO l = p-buff_size+1, p
                            DO k = -buff_size, n+buff_size
                                DO j = -buff_size, m+buff_size
                                    DO i = 1, sys_size
                                        r = sys_size*(j+buff_size)       &
                                          + sys_size*(m+2*buff_size+1) * &
                                            (k+buff_size) + sys_size *   &
                                            (m+2*buff_size+1) *          &
                                            (n+2*buff_size+1) *          &
                                            (l-p+buff_size-1) + (i-1)
                                        q_cons_buffer_out(r) = &
                                                q_cons_vf(i)%sf(j,k,l)
                                    END DO
                                END DO
                            END DO
                        END DO
                        
                        ! Sending/receiving the data to/from bc_z%end/bc_z%beg
                        CALL MPI_SENDRECV( q_cons_buffer_out(0), buff_size *  &
                                           sys_size*(m+2*buff_size+1) *       &
                                           (n+2*buff_size+1),                 &
                                           MPI_DOUBLE_PRECISION, bc_z%end, 0, &
                                           q_cons_buffer_in(0), buff_size *   &
                                           sys_size*(m+2*buff_size+1) *       &
                                           (n+2*buff_size+1),                 &
                                           MPI_DOUBLE_PRECISION, bc_z%beg, 0, &
                                           MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                           ierr                               )
                        
                    ! PBC only at beginning of the sub-domain
                    ELSE
                        
                        ! Packing the data to be sent to bc_z%beg
                        DO l = 0, buff_size-1
                            DO k = -buff_size, n+buff_size
                                DO j = -buff_size, m+buff_size
                                    DO i = 1, sys_size
                                        r = sys_size*(j+buff_size)       &
                                          + sys_size*(m+2*buff_size+1) * &
                                            (k+buff_size) + (i-1)        &
                                          + sys_size*(m+2*buff_size+1) * &
                                            (n+2*buff_size+1)*l
                                        q_cons_buffer_out(r) = &
                                                   q_cons_vf(i)%sf(j,k,l)
                                    END DO
                                END DO
                            END DO
                        END DO
                        
                        ! Sending/receiving the data to/from bc_z%beg/bc_z%beg
                        CALL MPI_SENDRECV( q_cons_buffer_out(0), buff_size *  &
                                           sys_size*(m+2*buff_size+1) *       &
                                           (n+2*buff_size+1),                 &
                                           MPI_DOUBLE_PRECISION, bc_z%beg, 1, &
                                           q_cons_buffer_in(0), buff_size *   &
                                           sys_size*(m+2*buff_size+1) *       &
                                           (n+2*buff_size+1),                 &
                                           MPI_DOUBLE_PRECISION, bc_z%beg, 0, &
                                           MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                           ierr                               )
                        
                    END IF
                    
                    ! Unpacking the data from bc_z%beg
                    DO l = -buff_size, -1
                        DO k = -buff_size, n+buff_size
                            DO j = -buff_size, m+buff_size
                                DO i = 1, sys_size
                                    r = sys_size*(j+buff_size)         &
                                      + sys_size*(m+2*buff_size+1) *   &
                                        (k+buff_size) + (i-1)          &
                                      + sys_size*(m+2*buff_size+1) *   &
                                        (n+2*buff_size+1)*(l+buff_size)
                                    q_cons_vf(i)%sf(j,k,l) = q_cons_buffer_in(r)
                                END DO
                            END DO
                        END DO
                    END DO
                    
                ELSE                         ! Buffer region at the end
                    
                    ! PBC at both ends of the sub-domain
                    IF(bc_z%beg >= 0) THEN
                        
                        ! Packing the data to be sent to bc_z%beg
                        DO l = 0, buff_size-1
                            DO k = -buff_size, n+buff_size
                                DO j = -buff_size, m+buff_size
                                    DO i = 1, sys_size
                                        r = sys_size*(j+buff_size)       &
                                          + sys_size*(m+2*buff_size+1) * &
                                            (k+buff_size) + (i-1)        &
                                          + sys_size*(m+2*buff_size+1) * &
                                            (n+2*buff_size+1)*l
                                        q_cons_buffer_out(r) = &
                                                   q_cons_vf(i)%sf(j,k,l)
                                    END DO
                                END DO
                            END DO
                        END DO
                        
                        ! Sending/receiving the data to/from bc_z%beg/bc_z%end
                        CALL MPI_SENDRECV( q_cons_buffer_out(0), buff_size *  &
                                           sys_size*(m+2*buff_size+1) *       &
                                           (n+2*buff_size+1),                 &
                                           MPI_DOUBLE_PRECISION, bc_z%beg, 1, &
                                           q_cons_buffer_in(0), buff_size *   &
                                           sys_size*(m+2*buff_size+1) *       &
                                           (n+2*buff_size+1),                 &
                                           MPI_DOUBLE_PRECISION, bc_z%end, 1, &
                                           MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                           ierr                               )
                        
                    ! PBC only at end of the sub-domain
                    ELSE
                        
                        ! Packing the data to be sent to bc_z%end
                        DO l = p-buff_size+1, p
                            DO k = -buff_size, n+buff_size
                                DO j = -buff_size, m+buff_size
                                    DO i = 1, sys_size
                                        r = sys_size*(j+buff_size)       &
                                          + sys_size*(m+2*buff_size+1) * &
                                            (k+buff_size) + sys_size *   &
                                            (m+2*buff_size+1) *          &
                                            (n+2*buff_size+1) *          &
                                            (l-p+buff_size-1) + (i-1)
                                        q_cons_buffer_out(r) = &
                                                q_cons_vf(i)%sf(j,k,l)
                                    END DO
                                END DO
                            END DO
                        END DO
                        
                        ! Sending/receiving the data to/from bc_z%end/bc_z%end
                        CALL MPI_SENDRECV( q_cons_buffer_out(0), buff_size *  &
                                           sys_size*(m+2*buff_size+1) *       &
                                           (n+2*buff_size+1),                 &
                                           MPI_DOUBLE_PRECISION, bc_z%end, 0, &
                                           q_cons_buffer_in(0), buff_size *   &
                                           sys_size*(m+2*buff_size+1) *       &
                                           (n+2*buff_size+1),                 &
                                           MPI_DOUBLE_PRECISION, bc_z%end, 1, &
                                           MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                           ierr                               )
                        
                    END IF
                    
                    ! Unpacking the data received from bc_z%end
                    DO l = p+1, p+buff_size
                        DO k = -buff_size, n+buff_size
                            DO j = -buff_size, m+buff_size
                                DO i = 1, sys_size
                                    r = sys_size*(j+buff_size)       &
                                      + sys_size*(m+2*buff_size+1) * &
                                        (k+buff_size) + (i-1)        &
                                      + sys_size*(m+2*buff_size+1) * &
                                        (n+2*buff_size+1)*(l-p-1)
                                    q_cons_vf(i)%sf(j,k,l) = q_cons_buffer_in(r)
                                END DO
                            END DO
                        END DO
                    END DO
                    
                END IF
                
            END IF
            
            ! END: Communications in the z-direction ===========================
            
            
        END SUBROUTINE s_mpi_sendrecv_cons_vars_buffer_regions ! ---------------
        
        
        
        
        !>  The following subroutine takes the first element of the
        !!      2-element inputted variable and determines its maximum
        !!      value on the entire computational domain. The result is
        !!      stored back into the first element of the variable while
        !!      the rank of the processor that is in charge of the sub-
        !!      domain containing the maximum is stored into the second
        !!      element of the variable.
        !!  @param var_loc On input, this variable holds the local value and processor rank,
        !!  which are to be reduced among all the processors in communicator.
        !!  On output, this variable holds the maximum value, reduced amongst
        !!  all of the local values, and the process rank to which the value
        !!  belongs.
        SUBROUTINE s_mpi_reduce_maxloc(var_loc) ! ------------------------------

            REAL(KIND(0d0)), DIMENSION(2), INTENT(INOUT) :: var_loc
            

            REAL(KIND(0d0)), DIMENSION(2) :: var_glb  !<
            !! Temporary storage variable that holds the reduced maximum value
            !! and the rank of the processor with which the value is associated            
            
            ! Performing reduction procedure and eventually storing its result
            ! into the variable that was initially inputted into the subroutine
            CALL MPI_REDUCE( var_loc, var_glb, 1, MPI_2DOUBLE_PRECISION, &
                             MPI_MAXLOC, 0, MPI_COMM_WORLD, ierr          )
            
            CALL MPI_BCAST( var_glb, 1, MPI_2DOUBLE_PRECISION, &
                                      0, MPI_COMM_WORLD, ierr   )
            
            var_loc = var_glb
            
        END SUBROUTINE s_mpi_reduce_maxloc ! -----------------------------------
        
        
        
        
        !>  This subroutine gathers the Silo database metadata for
        !!      the spatial extents in order to boost the performance of
        !!      the multidimensional visualization.  
        !!  @param spatial_extents Spatial extents for each processor's sub-domain. First dimension
        !!  corresponds to the minimum and maximum values, respectively, while
        !!  the second dimension corresponds to the processor rank.
        SUBROUTINE s_mpi_gather_spatial_extents(spatial_extents) ! -------------

            REAL(KIND(0d0)), DIMENSION(1:,0:), INTENT(INOUT) :: spatial_extents
            
            
            ! Simulation is 3D
            IF(p > 0) THEN
                IF (grid_geometry == 3) THEN
                    ! Minimum spatial extent in the r-direction
                    CALL MPI_GATHERV( MINVAL(y_cb), 1, MPI_DOUBLE_PRECISION,      &
                                      spatial_extents(1,0), recvcounts, 6*displs, &
                                      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,    &
                                      ierr                                        )
                    
                    ! Minimum spatial extent in the theta-direction
                    CALL MPI_GATHERV( MINVAL(z_cb), 1, MPI_DOUBLE_PRECISION,      &
                                      spatial_extents(2,0), recvcounts, 6*displs, &
                                      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,    &
                                      ierr                                        )
                    
                    ! Minimum spatial extent in the z-direction
                    CALL MPI_GATHERV( MINVAL(x_cb), 1, MPI_DOUBLE_PRECISION,      &
                                      spatial_extents(3,0), recvcounts, 6*displs, &
                                      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,    &
                                      ierr                                        )
                    
                    ! Maximum spatial extent in the r-direction
                    CALL MPI_GATHERV( MAXVAL(y_cb), 1, MPI_DOUBLE_PRECISION,      &
                                      spatial_extents(4,0), recvcounts, 6*displs, &
                                      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,    &
                                      ierr                                        )
                    
                    ! Maximum spatial extent in the theta-direction
                    CALL MPI_GATHERV( MAXVAL(z_cb), 1, MPI_DOUBLE_PRECISION,      &
                                      spatial_extents(5,0), recvcounts, 6*displs, &
                                      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,    &
                                      ierr                                        )
                    
                    ! Maximum spatial extent in the z-direction
                    CALL MPI_GATHERV( MAXVAL(x_cb), 1, MPI_DOUBLE_PRECISION,      &
                                      spatial_extents(6,0), recvcounts, 6*displs, &
                                      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,    &
                                      ierr                                        )
                ELSE
                    ! Minimum spatial extent in the x-direction
                    CALL MPI_GATHERV( MINVAL(x_cb), 1, MPI_DOUBLE_PRECISION,      &
                                      spatial_extents(1,0), recvcounts, 6*displs, &
                                      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,    &
                                      ierr                                        )
                    
                    ! Minimum spatial extent in the y-direction
                    CALL MPI_GATHERV( MINVAL(y_cb), 1, MPI_DOUBLE_PRECISION,      &
                                      spatial_extents(2,0), recvcounts, 6*displs, &
                                      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,    &
                                      ierr                                        )
                    
                    ! Minimum spatial extent in the z-direction
                    CALL MPI_GATHERV( MINVAL(z_cb), 1, MPI_DOUBLE_PRECISION,      &
                                      spatial_extents(3,0), recvcounts, 6*displs, &
                                      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,    &
                                      ierr                                        )
                    
                    ! Maximum spatial extent in the x-direction
                    CALL MPI_GATHERV( MAXVAL(x_cb), 1, MPI_DOUBLE_PRECISION,      &
                                      spatial_extents(4,0), recvcounts, 6*displs, &
                                      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,    &
                                      ierr                                        )
                    
                    ! Maximum spatial extent in the y-direction
                    CALL MPI_GATHERV( MAXVAL(y_cb), 1, MPI_DOUBLE_PRECISION,      &
                                      spatial_extents(5,0), recvcounts, 6*displs, &
                                      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,    &
                                      ierr                                        )
                    
                    ! Maximum spatial extent in the z-direction
                    CALL MPI_GATHERV( MAXVAL(z_cb), 1, MPI_DOUBLE_PRECISION,      &
                                      spatial_extents(6,0), recvcounts, 6*displs, &
                                      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,    &
                                      ierr                                        )
                END IF
            ! Simulation is 2D
            ELSE
                
                ! Minimum spatial extent in the x-direction
                CALL MPI_GATHERV( MINVAL(x_cb), 1, MPI_DOUBLE_PRECISION,      &
                                  spatial_extents(1,0), recvcounts, 4*displs, &
                                  MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,    &
                                  ierr                                        )
                
                ! Minimum spatial extent in the y-direction
                CALL MPI_GATHERV( MINVAL(y_cb), 1, MPI_DOUBLE_PRECISION,      &
                                  spatial_extents(2,0), recvcounts, 4*displs, &
                                  MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,    &
                                  ierr                                        )
                
                ! Maximum spatial extent in the x-direction
                CALL MPI_GATHERV( MAXVAL(x_cb), 1, MPI_DOUBLE_PRECISION,      &
                                  spatial_extents(3,0), recvcounts, 4*displs, &
                                  MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,    &
                                  ierr                                        )
                
                ! Maximum spatial extent in the y-direction
                CALL MPI_GATHERV( MAXVAL(y_cb), 1, MPI_DOUBLE_PRECISION,      &
                                  spatial_extents(4,0), recvcounts, 4*displs, &
                                  MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,    &
                                  ierr                                        )
                
            END IF
            
            
        END SUBROUTINE s_mpi_gather_spatial_extents ! --------------------------
        
        
        
        
        !>  This subroutine collects the sub-domain cell-boundary or
        !!      cell-center locations data from all of the processors and
        !!      puts back together the grid of the entire computational
        !!      domain on the rank 0 processor. This is only done for 1D
        !!      simulations.        
        SUBROUTINE s_mpi_defragment_1d_grid_variable() ! -----------------------

            ! Silo-HDF5 database format
            IF(format == 1) THEN
                
                CALL MPI_GATHERV( x_cc(0), m+1, MPI_DOUBLE_PRECISION,      &
                                  x_root_cc(0), recvcounts, displs,        &
                                  MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, &
                                  ierr                                     )
                
            ! Binary database format
            ELSE
                
                CALL MPI_GATHERV( x_cb(0), m+1, MPI_DOUBLE_PRECISION,      &
                                  x_root_cb(0), recvcounts, displs,        &
                                  MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, &
                                  ierr                                     )
                
                IF(proc_rank == 0) x_root_cb(-1) = x_cb(-1)
                
            END IF
            
        END SUBROUTINE s_mpi_defragment_1d_grid_variable ! ---------------------
        
        
        
        
        !>  This subroutine gathers the Silo database metadata for
        !!      the flow variable's extents as to boost performance of
        !!      the multidimensional visualization.    
        !!  @param q_sf Flow variable defined on a single computational sub-domain
        !!  @param data_extents The flow variable extents on each of the processor's sub-domain.
        !!   First dimension of array corresponds to the former's minimum and
        !!  maximum values, respectively, while second dimension corresponds
        !!  to each processor's rank.
        SUBROUTINE s_mpi_gather_data_extents(q_sf, data_extents) ! -------------

            
            REAL(KIND(0d0)), DIMENSION(:,:,:), INTENT(IN) :: q_sf
            
            REAL(KIND(0d0)), &
            DIMENSION(1:2,0:num_procs-1), &
            INTENT(INOUT) :: data_extents
            
            
            ! Mimimum flow variable extent
            CALL MPI_GATHERV( MINVAL(q_sf), 1, MPI_DOUBLE_PRECISION,        &
                              data_extents(1,0), recvcounts, 2*displs,      &
                              MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
            
            ! Maximum flow variable extent
            CALL MPI_GATHERV( MAXVAL(q_sf), 1, MPI_DOUBLE_PRECISION,        &
                              data_extents(2,0), recvcounts, 2*displs,      &
                              MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
            
            
        END SUBROUTINE s_mpi_gather_data_extents ! -----------------------------
        
        
        
        !>  This subroutine gathers the sub-domain flow variable data
        !!      from all of the processors and puts it back together for
        !!      the entire computational domain on the rank 0 processor.
        !!      This is only done for 1D simulations.       
        !!  @param q_sf Flow variable defined on a single computational sub-domain
        !!  @param q_root_sf Flow variable defined on the entire computational domain
        SUBROUTINE s_mpi_defragment_1d_flow_variable(q_sf, q_root_sf) ! --------

            REAL(KIND(0d0)), &
            DIMENSION(0:m,0:0,0:0), &
            INTENT(IN) :: q_sf
            
            REAL(KIND(0d0)), &
            DIMENSION(0:m_root,0:0,0:0), &
            INTENT(INOUT) :: q_root_sf
            
            ! Gathering the sub-domain flow variable data from all the processes
            ! and putting it back together for the entire computational domain
            ! on the process with rank 0
            CALL MPI_GATHERV( q_sf(0,0,0), m+1, MPI_DOUBLE_PRECISION,       &
                              q_root_sf(0,0,0), recvcounts, displs,         &
                              MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
            
            
        END SUBROUTINE s_mpi_defragment_1d_flow_variable ! ---------------------
        
        
        
        
        !> Deallocation procedures for the module 
        SUBROUTINE s_finalize_mpi_proxy_module() ! ---------------------------
            
            ! Deallocating the conservative variables buffer vectors
            IF(buff_size > 0) THEN
                DEALLOCATE( q_cons_buffer_in)
                DEALLOCATE(q_cons_buffer_out)
            END IF
            
            
            ! Deallocating the recieve counts and the displacement vector
            ! variables used in variable-gather communication procedures
            IF((format == 1 .AND. n > 0) .OR. n == 0) THEN
                DEALLOCATE(recvcounts)
                DEALLOCATE(displs)
            END IF
            
            
        END SUBROUTINE s_finalize_mpi_proxy_module ! -------------------------
        
        
        
        
       !> Finalization of all MPI related processes 
        SUBROUTINE s_mpi_finalize() ! ------------------------------
            
            ! Terminating the MPI environment
            CALL MPI_FINALIZE(ierr)
            
        END SUBROUTINE s_mpi_finalize ! ----------------------------
        
        
        
        
        
END MODULE m_mpi_proxy
