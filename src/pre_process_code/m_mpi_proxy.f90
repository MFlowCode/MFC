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

!> @brief This module serves as a proxy to the parameters and subroutines
!!              available in the MPI implementation's MPI module. Specifically,
!!              the role of the proxy is to harness basic MPI commands into more
!!              complex procedures as to achieve the required pre-processing
!!              communication goals.
MODULE m_mpi_proxy
    
    
    ! Dependencies =============================================================
    USE mpi                     !< Message passing interface (MPI) module

    USE m_derived_types         !< Definitions of the derived types
    
    USE m_global_parameters     !< Global parameters for the code
    ! ==========================================================================
    
    
    IMPLICIT NONE
    
    INTEGER, PRIVATE :: err_code, ierr !< 
    !! Generic flags used to identify and report MPI errors
    
    
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
        
        
        
        !! @param q_cons_vf Conservative variables 
        SUBROUTINE s_initialize_mpi_data(q_cons_vf) ! --------------------------

            TYPE(scalar_field), &
            DIMENSION(sys_size), &
            INTENT(IN) :: q_cons_vf

            INTEGER, DIMENSION(num_dims) :: sizes_glb, sizes_loc
            INTEGER :: ierr

            ! Generic loop iterator
            INTEGER :: i

            DO i = 1, sys_size
                MPI_IO_DATA%var(i)%sf => q_cons_vf(i)%sf
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



        !> Since only processor with rank 0 is in charge of reading
        !!       and checking the consistency of the user provided inputs,
        !!       these are not available to the remaining processors. This
        !!       subroutine is then in charge of broadcasting the required
        !!       information.
        SUBROUTINE s_mpi_bcast_user_inputs() ! ---------------------------------
           
            
            ! Generic loop iterator
            INTEGER :: i
            
            
            ! Logistics
            CALL MPI_BCAST( case_dir  , LEN(case_dir), MPI_CHARACTER ,      &
                                              0      , MPI_COMM_WORLD, ierr )
            CALL MPI_BCAST( old_grid  ,       1      , MPI_LOGICAL   ,      &
                                              0      , MPI_COMM_WORLD, ierr )
            CALL MPI_BCAST( old_ic    ,       1      , MPI_LOGICAL   ,      &
                                              0      , MPI_COMM_WORLD, ierr )
            CALL MPI_BCAST( t_step_old,       1      , MPI_INTEGER   ,      &
                                              0      , MPI_COMM_WORLD, ierr )
            
            
            ! Computational domain parameters
            CALL MPI_BCAST(m, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(p, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(m_glb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(n_glb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(p_glb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            
            CALL MPI_BCAST(cyl_coord, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            
            CALL MPI_BCAST( x_domain%beg, 1, MPI_DOUBLE_PRECISION, &
                                          0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST( x_domain%end, 1, MPI_DOUBLE_PRECISION, &
                                          0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST( y_domain%beg, 1, MPI_DOUBLE_PRECISION, &
                                          0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST( y_domain%end, 1, MPI_DOUBLE_PRECISION, &
                                          0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST( z_domain%beg, 1, MPI_DOUBLE_PRECISION, &
                                          0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST( z_domain%end, 1, MPI_DOUBLE_PRECISION, &
                                          0, MPI_COMM_WORLD, ierr  )
            
            CALL MPI_BCAST(stretch_x, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(stretch_y, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(stretch_z, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            
            CALL MPI_BCAST(a_x, 1,MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(a_y, 1,MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(a_z, 1,MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(loops_x, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(loops_y, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(loops_z, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(x_a, 1,MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(x_b, 1,MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(y_a, 1,MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(y_b, 1,MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(z_a, 1,MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(z_b, 1,MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
            
            
            ! Simulation algorithm parameters
            CALL MPI_BCAST(model_eqns, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(num_fluids, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(adv_alphan, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(mpp_lim, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(weno_order, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            
            CALL MPI_BCAST( bc_x%beg, 1, MPI_DOUBLE_PRECISION, &
                                      0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST( bc_x%end, 1, MPI_DOUBLE_PRECISION, &
                                      0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST( bc_y%beg, 1, MPI_DOUBLE_PRECISION, &
                                      0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST( bc_y%end, 1, MPI_DOUBLE_PRECISION, &
                                      0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST( bc_z%beg, 1, MPI_DOUBLE_PRECISION, &
                                      0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST( bc_z%end, 1, MPI_DOUBLE_PRECISION, &
                                      0, MPI_COMM_WORLD, ierr  )
            CALL MPI_BCAST(hypoelasticity, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

            CALL MPI_BCAST(parallel_io, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(precision, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(perturb_flow, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(perturb_flow_fluid, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(perturb_sph, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(perturb_sph_fluid, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(fluid_rho(1)     , num_fluids_max, MPI_LOGICAL   , &
                                                    0       , MPI_COMM_WORLD, &
                                                              ierr            )
           
            ! Initial condition parameters
            CALL MPI_BCAST(num_patches, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            
            DO i = 1, num_patches_max
                CALL MPI_BCAST( patch_icpp(i)%geometry       , 1, &
                                MPI_INTEGER                  , 0, &
                                MPI_COMM_WORLD, ierr              )
                CALL MPI_BCAST( patch_icpp(i)%x_centroid     , 1, &
                                MPI_DOUBLE_PRECISION         , 0, &
                                MPI_COMM_WORLD, ierr              )
                CALL MPI_BCAST( patch_icpp(i)%y_centroid     , 1, &
                                MPI_DOUBLE_PRECISION         , 0, &
                                MPI_COMM_WORLD, ierr              )
                CALL MPI_BCAST( patch_icpp(i)%z_centroid     , 1, &
                                MPI_DOUBLE_PRECISION         , 0, &
                                MPI_COMM_WORLD, ierr              )
                CALL MPI_BCAST( patch_icpp(i)%length_x       , 1, &
                                MPI_DOUBLE_PRECISION         , 0, &
                                MPI_COMM_WORLD, ierr              )
                CALL MPI_BCAST( patch_icpp(i)%length_y       , 1, &
                                MPI_DOUBLE_PRECISION         , 0, &
                                MPI_COMM_WORLD, ierr              )
                CALL MPI_BCAST( patch_icpp(i)%length_z       , 1, &
                                MPI_DOUBLE_PRECISION         , 0, &
                                MPI_COMM_WORLD, ierr              )
                CALL MPI_BCAST( patch_icpp(i)%radius         , 1, &
                                MPI_DOUBLE_PRECISION         , 0, &
                                MPI_COMM_WORLD, ierr              )
                CALL MPI_BCAST( patch_icpp(i)%epsilon        , 1, &
                                MPI_DOUBLE_PRECISION         , 0, &
                                MPI_COMM_WORLD, ierr              )
                CALL MPI_BCAST( patch_icpp(i)%beta           , 1, &
                                MPI_DOUBLE_PRECISION         , 0, &
                                MPI_COMM_WORLD, ierr              )
                CALL MPI_BCAST( patch_icpp(i)%normal(1)      , 3, &
                                MPI_DOUBLE_PRECISION         , 0, &
                                MPI_COMM_WORLD, ierr              )
                CALL MPI_BCAST( patch_icpp(i)%radii(1)       , 3, &
                                MPI_DOUBLE_PRECISION         , 0, &
                                MPI_COMM_WORLD, ierr              )
                CALL MPI_BCAST( patch_icpp(i)%smoothen       , 1, &
                                MPI_LOGICAL                  , 0, &
                                MPI_COMM_WORLD, ierr              )
                CALL MPI_BCAST( patch_icpp(i)%smooth_patch_id, 1, &
                                MPI_INTEGER                  , 0, &
                                MPI_COMM_WORLD, ierr              )
                CALL MPI_BCAST( patch_icpp(i)%smooth_coeff   , 1, &
                                MPI_DOUBLE_PRECISION         , 0, &
                                MPI_COMM_WORLD, ierr              )
                CALL MPI_BCAST( patch_icpp(i)%rho            , 1, &
                                MPI_DOUBLE_PRECISION         , 0, &
                                MPI_COMM_WORLD, ierr              )
                CALL MPI_BCAST( patch_icpp(i)%vel(1)         , 3, &
                                MPI_DOUBLE_PRECISION         , 0, &
                                MPI_COMM_WORLD, ierr              )
                CALL MPI_BCAST( patch_icpp(i)%pres           , 1, &
                                MPI_DOUBLE_PRECISION         , 0, &
                                MPI_COMM_WORLD, ierr              )
                CALL MPI_BCAST( patch_icpp(i)%gamma          , 1, &
                                MPI_DOUBLE_PRECISION         , 0, &
                                MPI_COMM_WORLD, ierr              )
                CALL MPI_BCAST( patch_icpp(i)%pi_inf         , 1, &
                                MPI_DOUBLE_PRECISION         , 0, &
                                MPI_COMM_WORLD, ierr              )
                CALL MPI_BCAST( patch_icpp(i)%tau_e(1)       , 6, &
                                MPI_DOUBLE_PRECISION         , 0, &
                                MPI_COMM_WORLD, ierr              )
                CALL MPI_BCAST( patch_icpp(i)%alter_patch(0)             , &
                                num_patches_max , MPI_LOGICAL         , 0, &
                                MPI_COMM_WORLD, ierr                       )
                CALL MPI_BCAST( patch_icpp(i)%alpha_rho(1)               , &
                                num_fluids_max  , MPI_DOUBLE_PRECISION, 0, &
                                MPI_COMM_WORLD, ierr                       )
                CALL MPI_BCAST( patch_icpp(i)%alpha(1)                   , &
                                num_fluids_max-1, MPI_DOUBLE_PRECISION, 0, &
                                MPI_COMM_WORLD, ierr                       )

                ! Bubbles
                CALL MPI_BCAST( patch_icpp(i)%r0         , 1, &
                                MPI_DOUBLE_PRECISION         , 0, &
                                MPI_COMM_WORLD, ierr              )
                CALL MPI_BCAST( patch_icpp(i)%v0         , 1, &
                                MPI_DOUBLE_PRECISION         , 0, &
                                MPI_COMM_WORLD, ierr              )
                CALL MPI_BCAST( patch_icpp(i)%p0         , 1, &
                                MPI_DOUBLE_PRECISION         , 0, &
                                MPI_COMM_WORLD, ierr              )
                CALL MPI_BCAST( patch_icpp(i)%m0         , 1, &
                                MPI_DOUBLE_PRECISION         , 0, &
                                MPI_COMM_WORLD, ierr              )



                 
            END DO
            
            
            ! Fluids physical parameters
            DO i = 1, num_fluids_max
                CALL MPI_BCAST( fluid_pp(i)%gamma   , 1, &
                                MPI_DOUBLE_PRECISION, 0, &
                                MPI_COMM_WORLD, ierr     )
                CALL MPI_BCAST( fluid_pp(i)%pi_inf  , 1, &
                                MPI_DOUBLE_PRECISION, 0, &
                                MPI_COMM_WORLD, ierr     )

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
                CALL MPI_BCAST( fluid_pp(i)%G   , 1, &
                                MPI_DOUBLE_PRECISION, 0, &
                                MPI_COMM_WORLD, ierr    )
            END DO
            
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
            CALL MPI_BCAST( polydisperse,1,          &
                        MPI_LOGICAL,0,          &
                        MPI_COMM_WORLD,ierr  )
            CALL MPI_BCAST( poly_sigma,1,            &
                        MPI_DOUBLE_PRECISION,0, &
                        MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST( thermal,1,            &
                        MPI_INTEGER,0, &
                        MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST( R0ref,1,            &
                        MPI_DOUBLE_PRECISION,0, &
                        MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST( nb,1,            &
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


            CALL MPI_BCAST( qbmm,1,          &
                        MPI_LOGICAL,0,          &
                        MPI_COMM_WORLD,ierr  )
            CALL MPI_BCAST( nnode,1,            &
                        MPI_INTEGER,0, &
                        MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST( sigR,1,            &
                        MPI_DOUBLE_PRECISION,0, &
                        MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST( sigV,1,            &
                        MPI_DOUBLE_PRECISION,0, &
                        MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST( rhoRV,1,            &
                        MPI_DOUBLE_PRECISION,0, &
                        MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST( dist_type,1,            &
                        MPI_INTEGER,0, &
                        MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST( R0_type,1,            &
                        MPI_INTEGER,0, &
                        MPI_COMM_WORLD,ierr)



        END SUBROUTINE s_mpi_bcast_user_inputs ! -------------------------------
        
        
        
        !> Description: This subroutine takes care of efficiently distributing
        !!              the computational domain among the available processors
        !!             as well as recomputing some of the global parameters so
        !!              that they reflect the configuration of sub-domain that is
        !!              overseen by the local processor.        
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
                    
                    ! Preliminary uniform cell-width spacing
                    IF(old_grid .NEQV. .TRUE.) THEN
                        dz = (z_domain%end - z_domain%beg) / REAL(p+1,KIND(0d0))
                    END IF
                    
                    ! Optimal number of cells per processor
                    p = (p+1) / num_procs_z - 1
                    
                    ! Distributing any remaining cells
                    DO i = 1, rem_cells
                        IF(proc_coords(3) == i-1) THEN
                            p = p + 1
                            EXIT
                        END IF
                    END DO
                    
                    ! Beginning and end sub-domain boundary locations
                    IF (parallel_io .NEQV. .TRUE.) THEN
                                IF(old_grid .NEQV. .TRUE.) THEN
                                    IF(proc_coords(3) < rem_cells) THEN
                                        z_domain%beg = z_domain%beg + dz*REAL( (p+1) * &
                                                       proc_coords(3) )
                                        z_domain%end = z_domain%end - dz*REAL( (p+1) * &
                                                       (num_procs_z - proc_coords(3) - 1) &
                                                     - (num_procs_z - rem_cells) )
                                    ELSE
                                        z_domain%beg = z_domain%beg + dz*REAL( (p+1) * &
                                                       proc_coords(3) + rem_cells )
                                        z_domain%end = z_domain%end - dz*REAL( (p+1) * &
                                                       (num_procs_z - proc_coords(3) - 1) )
                                    END IF
                                END IF
                    ELSE
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
                
                ! Preliminary uniform cell-width spacing
                IF(old_grid .NEQV. .TRUE.) THEN
                    dy = (y_domain%end - y_domain%beg) / REAL(n+1, KIND(0d0))
                END IF
                
                ! Optimal number of cells per processor
                n = (n+1) / num_procs_y - 1
                
                ! Distributing any remaining cells
                DO i = 1, rem_cells
                    IF(proc_coords(2) == i-1) THEN
                        n = n + 1
                        EXIT
                    END IF
                END DO
                
                ! Beginning and end sub-domain boundary locations
                IF (parallel_io .NEQV. .TRUE.) THEN
                            IF(old_grid .NEQV. .TRUE.) THEN
                                IF(proc_coords(2) < rem_cells) THEN
                                    y_domain%beg = y_domain%beg + dy*REAL( (n+1) * &
                                                   proc_coords(2) )
                                    y_domain%end = y_domain%end - dy*REAL( (n+1) * &
                                                   (num_procs_y - proc_coords(2) - 1) &
                                                 - (num_procs_y - rem_cells) )
                                ELSE
                                    y_domain%beg = y_domain%beg + dy*REAL( (n+1) * &
                                                   proc_coords(2) + rem_cells )
                                    y_domain%end = y_domain%end - dy*REAL( (n+1) * &
                                                   (num_procs_y - proc_coords(2) - 1) )
                                END IF
                            END IF
                ELSE
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
            
            ! Preliminary uniform cell-width spacing
            IF(old_grid .NEQV. .TRUE.) THEN
                dx = (x_domain%end - x_domain%beg) / REAL(m+1, KIND(0d0))
            END IF
            
            ! Optimal number of cells per processor
            m = (m+1) / num_procs_x - 1
            
            ! Distributing any remaining cells
            DO i = 1, rem_cells
                IF(proc_coords(1) == i-1) THEN
                    m = m + 1
                    EXIT
                END IF
            END DO
            
            ! Beginning and end sub-domain boundary locations
            IF (parallel_io .NEQV. .TRUE.) THEN
                    IF(old_grid .NEQV. .TRUE.) THEN
                        IF(proc_coords(1) < rem_cells) THEN
                            x_domain%beg = x_domain%beg + dx*REAL( (m+1) * &
                                           proc_coords(1) )
                            x_domain%end = x_domain%end - dx*REAL( (m+1) * &
                                           (num_procs_x - proc_coords(1) - 1) &
                                         - (num_procs_x - rem_cells) )
                        ELSE
                            x_domain%beg = x_domain%beg + dx*REAL( (m+1) * &
                                           proc_coords(1) + rem_cells )
                            x_domain%end = x_domain%end - dx*REAL( (m+1) * &
                                           (num_procs_x - proc_coords(1) - 1) )
                        END IF
                    END IF
            ELSE
                IF (proc_coords(1) < rem_cells) THEN
                    start_idx(1) = (m+1) * proc_coords(1)
                ELSE
                    start_idx(1) = (m+1) * proc_coords(1) + rem_cells
                END IF
            END IF
            
            ! ==================================================================
           
            
            
        END SUBROUTINE s_mpi_decompose_computational_domain ! ------------------
        
        
        !>  The following subroutine takes the inputted variable and
        !!      determines its minimum value on the entire computational
        !!      domain. The result is stored back into inputted variable.        
        !!  @param var_loc holds the local value to be reduced among
        !!      all the processors in communicator. On output, the variable holds
        !!      the minimum value, reduced amongst all of the local values.
        SUBROUTINE s_mpi_reduce_min(var_loc) ! ---------------------------------

            REAL(KIND(0d0)), INTENT(INOUT) :: var_loc
            
            ! Temporary storage variable that holds the reduced minimum value
            REAL(KIND(0d0)) :: var_glb
            
            
            ! Performing reduction procedure and eventually storing its result
            ! into the variable that was initially inputted into the subroutine
            CALL MPI_REDUCE( var_loc, var_glb, 1, MPI_DOUBLE_PRECISION, &
                             MPI_MIN, 0, MPI_COMM_WORLD, ierr            )
            
            CALL MPI_BCAST( var_glb, 1, MPI_DOUBLE_PRECISION, &
                                      0, MPI_COMM_WORLD, ierr  )
            
            var_loc = var_glb
            
            
        END SUBROUTINE s_mpi_reduce_min ! --------------------------------------
        
        
        
        
        !> Finalization of all MPI related processes
        SUBROUTINE s_mpi_finalize() ! ------------------------------

            ! Terminating the MPI environment
            CALL MPI_FINALIZE(ierr)
            
        END SUBROUTINE s_mpi_finalize ! ----------------------------
        
END MODULE m_mpi_proxy
