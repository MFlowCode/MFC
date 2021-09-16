!>
!! @file m_mpi_proxy.f90
!! @brief Contains module m_mpi_proxy 
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief The purpose of this module is to behave as a surrogate to the
!!	    m_mpi_proxy.f90 modules located in the pre-process, simulation,
!!	    and post-process code components, when a platform with parallel
!!	    computation capabilities is unavailable. The use of this module
!!	    instead of the ones located in the individual components' folders
!!	    should then allow the code to be compiled and executed in serial.
!!  INSTRUCTIONS: To use this module for the serial compilation and execution
!!        of one of the MFC components on a platform with no support for parallel
!!        computations, simply replace the existing m_mpi_proxy.f90 module with this
!!        one in the selected component's directory. Additionally, the makefile
!!        located in the component's folder should be modified so that the serial
!!        compiler is used to produce the executable.
MODULE m_mpi_proxy


        ! Dependencies =============================================================
        USE m_derived_types             !< Definitions of the derived types
        
        USE m_global_parameters		!< Global parameters for the code
        ! ==========================================================================
        
        
        IMPLICIT NONE
        

        CONTAINS
                !>  Refer to the m_mpi_proxy.f90 modules in the pre-process,
                !!	simulation, and/or post-process code directories for
                !!	information about this subroutine
                SUBROUTINE s_initalize_mpi_environment() ! -----------------------------

			! Serial run only has 1 processor
			num_procs = 1
			! Local processor rank is 0
			proc_rank = 0
		
                END SUBROUTINE s_initialize_mpi_environment ! --------------------------
		
		
		
		
		!> Refer to the m_mpi_proxy.f90 modules in the simulation
		!!			   and/or post-process code directories for information
		!!			   about this subroutine
		SUBROUTINE s_prepare_mpi_proxy_module() ! ------------------------------
		
			
			PRINT '(A)', 'Requesting module preparation procedures is ' // &
						 'inconsistent with use of MPI proxy module for ' // &
						 'serial execution. Exiting ...'
			STOP
			
			
		END SUBROUTINE s_prepare_mpi_proxy_module ! ----------------------------
		
		
		
		!> Refer to the m_mpi_proxy.f90 modules in the simulation
		!!			   and/or post-process code directories for information
		!!			   about this subroutine	
		SUBROUTINE s_mpi_bcast_user_inputs() ! ---------------------------------
			
			
			PRINT '(A)', 'Broadcasting user inputs is inconsistent with ' // &
						 'use of MPI proxy module for serial execution. ' // &
						 'Exiting ...'
			STOP
			
			
		END SUBROUTINE s_mpi_bcast_user_inputs ! -------------------------------
		
		
		
		
		!> Refer to the m_mpi_proxy.f90 modules in the simulation
		!!			   and/or post-process code directories for information
		!!			   about this subroutine	
		SUBROUTINE s_mpi_decompose_computational_domain() ! --------------------
			
			
			PRINT '(A)', 'Parallel computational domain decomposition is ' // &
						 'inconsistent with use of MPI proxy module for ' // &
						 'serial execution. Exiting ...'
			STOP
			
			
		END SUBROUTINE s_mpi_decompose_computational_domain ! ------------------
		
		
		
		!> Refer to the m_mpi_proxy.f90 modules in the simulation
		!!			   and/or post-process code directories for information
		!!			   about this subroutine	
		SUBROUTINE s_mpi_sendrecv_grid_vars_buffer_regions(pbc_loc, sweep_coord)
			
			CHARACTER(LEN = *), INTENT(IN) :: pbc_loc
			CHARACTER		  , INTENT(IN) :: sweep_coord
			
			
			PRINT '(A)', 'Communication of buffer regions for grid ' // &
						 'variables is inconsistent with use of MPI ' // &
						 'proxy module for serial execution. Exiting ...'
			STOP
			
			
		END SUBROUTINE s_mpi_sendrecv_grid_vars_buffer_regions ! ---------------
		
		
		
		
		!> Refer to the m_mpi_proxy.f90 modules in the simulation
		!!			   and/or post-process code directories for information
		!!			   about this subroutine	
		SUBROUTINE s_mpi_reduce_maxloc(loc_var) ! ------------------------------
			
			
			REAL(KIND(0d0)), DIMENSION(:), INTENT(IN) :: loc_var
			
			
			PRINT '(A)', 'Multiprocessor maximum location reduction ' // &
						 'procedure is inconsistent with use of MPI proxy ' // &
						 'module for serial execution. Exiting ...'
			STOP
			
			
		END SUBROUTINE s_mpi_reduce_maxloc ! -----------------------------------
		
		
		
		!> Refer to the m_mpi_proxy.f90 modules in the simulation
		!!			   and/or post-process code directories for information
		!!			   about this subroutine	
		SUBROUTINE s_mpi_reduce_min(loc_var) ! ---------------------------------
			
			REAL(KIND(0d0)), INTENT(IN) :: loc_var
			
			
			PRINT '(A)', 'Multiprocessor minimum reduction procedure is ' // &
						 'inconsistent with use of MPI proxy module for ' // &
						 'serial execution. Exiting ...'
			STOP
			
			
		END SUBROUTINE s_mpi_reduce_min ! --------------------------------------
		
		
		
		!> Refer to the m_mpi_proxy.f90 modules in the simulation
		!!			   and/or post-process code directories for information
		!!			   about this subroutine	
		SUBROUTINE s_mpi_sendrecv_cons_vars_buffer_regions(q_cons_vp, pbc_loc, &
							   	   sweep_coord	   )
			
			TYPE(field_position), DIMENSION(:), INTENT(IN) :: q_cons_vp
			CHARACTER(LEN = *)				  , INTENT(IN) :: pbc_loc
			CHARACTER						  , INTENT(IN) :: sweep_coord
			
			
			PRINT '(A)', 'Communication of buffer regions for conservative '// &
						 'variables is inconsistent with use of MPI proxy ' // &
						 'module for serial execution. Exiting ...'
			STOP
			
			
		END SUBROUTINE s_mpi_sendrecv_cons_vars_buffer_regions ! ---------------
		
		
		
		!> Refer to the m_mpi_proxy.f90 modules in the simulation
		!!			   and/or post-process code directories for information
		!!			   about this subroutine	
		SUBROUTINE s_mpi_gather_spatial_extents(spatial_extents) ! -------------
			
			
			REAL(KIND(0d0)), DIMENSION(:,:), INTENT(IN) :: spatial_extents
			
			
			PRINT '(A)', 'Gathering of the grids spatial extents is ' // &
						 'inconsistent with use of MPI proxy module for ' // &
						 'serial execution. Exiting ...'
			STOP
			
			
		END SUBROUTINE s_mpi_gather_spatial_extents ! --------------------------
		
		
		
		
		!> Refer to the m_mpi_proxy.f90 modules in the simulation
		!!			   and/or post-process code directories for information
		!!			   about this subroutine	
                SUBROUTINE s_mpi_defragment_1d_grid_variable() ! -----------------------
			
			
			PRINT '(A)', 'Defragmentation of 1D grid variable is ' // &
						 'inconsistent with use of MPI proxy module for ' // &
						 'serial execution. Exiting ...'
			STOP
			
			
		END SUBROUTINE s_mpi_defragment_1d_grid_variable ! ---------------------
		
		
		
		!> Refer to the m_mpi_proxy.f90 modules in the simulation
		!!			   and/or post-process code directories for information
		!!			   about this subroutine			
		SUBROUTINE s_mpi_gather_data_extents(q_fp, data_extents) ! -------------
			
			
			REAL(KIND(0d0)), DIMENSION(:,:,:), INTENT(IN) :: q_fp
			REAL(KIND(0d0)), DIMENSION(:,:)	 , INTENT(IN) :: data_extents
			
			
			PRINT '(A)', 'Gathering of a flow variables extents is ' // &
						 'inconsistent with use of MPI proxy module for ' // &
						 'serial execution. Exiting ...'
			STOP
			
			
		END SUBROUTINE s_mpi_gather_data_extents ! -----------------------------
		
		
		
		!> Refer to the m_mpi_proxy.f90 modules in the simulation
		!!			   and/or post-process code directories for information
		!!			   about this subroutine	
		SUBROUTINE s_mpi_defragment_1d_flow_variable(q_fp, q_root_fp) ! --------
			
			REAL(KIND(0d0)), DIMENSION(:,:,:), INTENT(IN) :: q_fp
			REAL(KIND(0d0)), DIMENSION(:,:,:), INTENT(IN) :: q_root_fp
			
			
			PRINT '(A)', 'Defragmentation of 1D flow variable is ' // &
						 'inconsistent with use of MPI proxy module for ' // &
						 'serial execution. Exiting ...'
			STOP
			
		END SUBROUTINE s_mpi_defragment_1d_flow_variable ! ---------------------
		
		
		
		
		!> Refer to the m_mpi_proxy.f90 modules in the simulation
		!!			   and/or post-process code directories for information
		!!			   about this subroutine			
		SUBROUTINE s_deallocate_mpi_proxy_module() ! ---------------------------
			
			PRINT '(A)', 'Requesting module deallocation procedures is ' // &
						 'inconsistent with use of MPI proxy module for ' // &
						 'serial execution. Exiting ...'
			STOP
			
			
		END SUBROUTINE s_deallocate_mpi_proxy_module ! -------------------------
		
		
		
		!> Refer to the m_mpi_proxy.f90 modules in the simulation
		!!			   and/or post-process code directories for information
		!!			   about this subroutine		
		SUBROUTINE s_finalize_mpi_environment() ! ------------------------------
			
			RETURN
			
		END SUBROUTINE s_finalize_mpi_environment ! ----------------------------
		
		
		
		
		
END MODULE m_mpi_proxy
