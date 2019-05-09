! MFC v3.0 - Master Scripts: m_mpi_proxy.f90
! Description: The purpose of this module is to behave as a surrogate to the
!			   m_mpi_proxy.f90 modules located in the pre-process, simulation,
!			   and post-process code components, when a platform with parallel
!			   computation capabilities is unavailable. The use of this module
!			   instead of the ones located in the individual components' folders
!			   should then allow the code to be compiled and executed in serial.
! Author: Vedran Coralic
! Date: 06/08/12


MODULE m_mpi_proxy
	
	
	! Dependencies =============================================================
	USE m_derived_types         ! Definitions of the derived types
	
	USE m_global_parameters		! Global parameters for the code
	! ==========================================================================
	
	
	IMPLICIT NONE
	
	
	! INSTRUCTIONS: To use this module for the serial compilation and execution
	! of one of the MFC components on a platform with no support for parallel
	! computations, simply replace the existing m_mpi_proxy.f90 module with this
	! one in the selected component's directory. Additionally, the makefile
	! located in the component's folder should be modified so that the serial
	! compiler is used to produce the executable.
	
	
	CONTAINS
		
		
		
		
		
		SUBROUTINE s_initalize_mpi_environment() ! -----------------------------
		! Description: Refer to the m_mpi_proxy.f90 modules in the pre-process,
		!			   simulation, and/or post-process code directories for
		!			   information about this subroutine
			
			
			! Serial run only has 1 processor
			num_procs = 1
			
			! Local processor rank is 0
			proc_rank = 0
			
			
		END SUBROUTINE s_initialize_mpi_environment ! --------------------------
		
		
		
		
		
		SUBROUTINE s_prepare_mpi_proxy_module() ! ------------------------------
		! Description: Refer to the m_mpi_proxy.f90 modules in the simulation
		!			   and/or post-process code directories for information
		!			   about this subroutine
			
			
			PRINT '(A)', 'Requesting module preparation procedures is ' // &
						 'inconsistent with use of MPI proxy module for ' // &
						 'serial execution. Exiting ...'
			STOP
			
			
		END SUBROUTINE s_prepare_mpi_proxy_module ! ----------------------------
		
		
		
		
		
		SUBROUTINE s_mpi_bcast_user_inputs() ! ---------------------------------
		! Description: Refer to the m_mpi_proxy.f90 modules in the pre-process,
		!			   simulation, and/or post-process code directories for
		!			   information about this subroutine
			
			
			PRINT '(A)', 'Broadcasting user inputs is inconsistent with ' // &
						 'use of MPI proxy module for serial execution. ' // &
						 'Exiting ...'
			STOP
			
			
		END SUBROUTINE s_mpi_bcast_user_inputs ! -------------------------------
		
		
		
		
		
		SUBROUTINE s_mpi_decompose_computational_domain() ! --------------------
		! Description: Refer to the m_mpi_proxy.f90 modules in the pre-process,
		!			   simulation, and/or post-process code directories for
		!			   information about this subroutine
			
			
			PRINT '(A)', 'Parallel computational domain decomposition is ' // &
						 'inconsistent with use of MPI proxy module for ' // &
						 'serial execution. Exiting ...'
			STOP
			
			
		END SUBROUTINE s_mpi_decompose_computational_domain ! ------------------
		
		
		
		
		
		SUBROUTINE s_mpi_sendrecv_grid_vars_buffer_regions(pbc_loc, sweep_coord)
		! Description: Refer to the m_mpi_proxy.f90 modules in the simulation
		!			   and/or post-process code directories for information
		!			   about this subroutine
			
			
			CHARACTER(LEN = *), INTENT(IN) :: pbc_loc
			CHARACTER		  , INTENT(IN) :: sweep_coord
			
			
			PRINT '(A)', 'Communication of buffer regions for grid ' // &
						 'variables is inconsistent with use of MPI ' // &
						 'proxy module for serial execution. Exiting ...'
			STOP
			
			
		END SUBROUTINE s_mpi_sendrecv_grid_vars_buffer_regions ! ---------------
		
		
		
		
		
		SUBROUTINE s_mpi_reduce_maxloc(loc_var) ! ------------------------------
		! Description: Refer to the m_mpi_proxy.f90 module in the simulation
		!			   code directory for information about this subroutine
			
			
			REAL(KIND(0d0)), DIMENSION(:), INTENT(IN) :: loc_var
			
			
			PRINT '(A)', 'Multiprocessor maximum location reduction ' // &
						 'procedure is inconsistent with use of MPI proxy ' // &
						 'module for serial execution. Exiting ...'
			STOP
			
			
		END SUBROUTINE s_mpi_reduce_maxloc ! -----------------------------------
		
		
		
		
		
		SUBROUTINE s_mpi_reduce_min(loc_var) ! ---------------------------------
		! Description: Refer to the m_mpi_proxy.f90 module in the pre-process
		!			   code directory for information about this subroutine
			
			
			REAL(KIND(0d0)), INTENT(IN) :: loc_var
			
			
			PRINT '(A)', 'Multiprocessor minimum reduction procedure is ' // &
						 'inconsistent with use of MPI proxy module for ' // &
						 'serial execution. Exiting ...'
			STOP
			
			
		END SUBROUTINE s_mpi_reduce_min ! --------------------------------------
		
		
		
		
		
		SUBROUTINE s_mpi_sendrecv_cons_vars_buffer_regions(q_cons_vp, pbc_loc, &
														   	   sweep_coord	   )
		! Description: Refer to the m_mpi_proxy.f90 module in the post-process
		!			   code directory for information about this subroutine
			
			
			TYPE(field_position), DIMENSION(:), INTENT(IN) :: q_cons_vp
			CHARACTER(LEN = *)				  , INTENT(IN) :: pbc_loc
			CHARACTER						  , INTENT(IN) :: sweep_coord
			
			
			PRINT '(A)', 'Communication of buffer regions for conservative '// &
						 'variables is inconsistent with use of MPI proxy ' // &
						 'module for serial execution. Exiting ...'
			STOP
			
			
		END SUBROUTINE s_mpi_sendrecv_cons_vars_buffer_regions ! ---------------
		
		
		
		
		
		SUBROUTINE s_mpi_gather_spatial_extents(spatial_extents) ! -------------
		! Description: Refer to the m_mpi_proxy.f90 module in the post-process
		!			   code directory for information about this subroutine
			
			
			REAL(KIND(0d0)), DIMENSION(:,:), INTENT(IN) :: spatial_extents
			
			
			PRINT '(A)', 'Gathering of the grid\'s spatial extents is ' // &
						 'inconsistent with use of MPI proxy module for ' // &
						 'serial execution. Exiting ...'
			STOP
			
			
		END SUBROUTINE s_mpi_gather_spatial_extents ! --------------------------
		
		
		
		
		
		SUBROUTINE s_mpi_defragment_1d_grid_variable() ! -----------------------
		! Description: Refer to the m_mpi_proxy.f90 module in the post-process
		!			   code directory for information about this subroutine
			
			
			PRINT '(A)', 'Defragmentation of 1D grid variable is ' // &
						 'inconsistent with use of MPI proxy module for ' // &
						 'serial execution. Exiting ...'
			STOP
			
			
		END SUBROUTINE s_mpi_defragment_1d_grid_variable ! ---------------------
		
		
		
		
		
		SUBROUTINE s_mpi_gather_data_extents(q_fp, data_extents) ! -------------
		! Description: Refer to the m_mpi_proxy.f90 module in the post-process
		!			   code directory for information about this subroutine
			
			
			REAL(KIND(0d0)), DIMENSION(:,:,:), INTENT(IN) :: q_fp
			REAL(KIND(0d0)), DIMENSION(:,:)	 , INTENT(IN) :: data_extents
			
			
			PRINT '(A)', 'Gathering of a flow variable\'s extents is ' // &
						 'inconsistent with use of MPI proxy module for ' // &
						 'serial execution. Exiting ...'
			STOP
			
			
		END SUBROUTINE s_mpi_gather_data_extents ! -----------------------------
		
		
		
		
		
		SUBROUTINE s_mpi_defragment_1d_flow_variable(q_fp, q_root_fp) ! --------
		! Description: Refer to the m_mpi_proxy.f90 module in the post-process
		!			   code directory for information about this subroutine
			
			
			REAL(KIND(0d0)), DIMENSION(:,:,:), INTENT(IN) :: q_fp
			REAL(KIND(0d0)), DIMENSION(:,:,:), INTENT(IN) :: q_root_fp
			
			
			PRINT '(A)', 'Defragmentation of 1D flow variable is ' // &
						 'inconsistent with use of MPI proxy module for ' // &
						 'serial execution. Exiting ...'
			STOP
			
			
		END SUBROUTINE s_mpi_defragment_1d_flow_variable ! ---------------------
		
		
		
		
		
		SUBROUTINE s_deallocate_mpi_proxy_module() ! ---------------------------
		! Description: Refer to the m_mpi_proxy.f90 modules in the simulation
		!			   and/or post-process code directories for information
		!			   about this subroutine
			
			
			PRINT '(A)', 'Requesting module deallocation procedures is ' // &
						 'inconsistent with use of MPI proxy module for ' // &
						 'serial execution. Exiting ...'
			STOP
			
			
		END SUBROUTINE s_deallocate_mpi_proxy_module ! -------------------------
		
		
		
		
		
		SUBROUTINE s_finalize_mpi_environment() ! ------------------------------
		! Description: Refer to the m_mpi_proxy.f90 modules in the pre-process,
		!			   simulation, and/or post-process code directories for
		!			   information about this subroutine
			
			
			RETURN
			
			
		END SUBROUTINE s_finalize_mpi_environment ! ----------------------------
		
		
		
		
		
END MODULE m_mpi_proxy
