!>
!! @file m_mpi_proxy.f90
!! @brief Contains module m_mpi_proxy
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief The purpose of this module is to behave as a surrogate to the
!!            m_mpi_proxy.f90 modules located in the pre-process, simulation,
!!            and post-process code components, when a platform with parallel
!!            computation capabilities is unavailable. The use of this module
!!            instead of the ones located in the individual components' folders
!!            should then allow the code to be compiled and executed in serial.
!!  INSTRUCTIONS: To use this module for the serial compilation and execution
!!        of one of the MFC components on a platform with no support for parallel
!!        computations, simply replace the existing m_mpi_proxy.f90 module with this
!!        one in the selected component's directory. Additionally, the makefile
!!        located in the component's folder should be modified so that the serial
!!        compiler is used to produce the executable.
module m_mpi_proxy

    ! Dependencies =============================================================
    use m_derived_types             !< Definitions of the derived types

    use m_global_parameters                !< Global parameters for the code
    ! ==========================================================================

    implicit none

contains
    !>  Refer to the m_mpi_proxy.f90 modules in the pre-process,
                !!        simulation, and/or post-process code directories for
                !!        information about this subroutine
    subroutine s_initalize_mpi_environment() ! -----------------------------

        ! Serial run only has 1 processor
        num_procs = 1
        ! Local processor rank is 0
        proc_rank = 0

    end subroutine s_initialize_mpi_environment ! --------------------------

    !> Refer to the m_mpi_proxy.f90 modules in the simulation
                !!                           and/or post-process code directories for information
                !!                           about this subroutine
    subroutine s_prepare_mpi_proxy_module() ! ------------------------------

        print '(A)', 'Requesting module preparation procedures is '// &
            'inconsistent with use of MPI proxy module for '// &
            'serial execution. Exiting ...'
        stop

    end subroutine s_prepare_mpi_proxy_module ! ----------------------------

    !> Refer to the m_mpi_proxy.f90 modules in the simulation
                !!                           and/or post-process code directories for information
                !!                           about this subroutine
    subroutine s_mpi_bcast_user_inputs() ! ---------------------------------

        print '(A)', 'Broadcasting user inputs is inconsistent with '// &
            'use of MPI proxy module for serial execution. '// &
            'Exiting ...'
        stop

    end subroutine s_mpi_bcast_user_inputs ! -------------------------------

    !> Refer to the m_mpi_proxy.f90 modules in the simulation
                !!                           and/or post-process code directories for information
                !!                           about this subroutine
    subroutine s_mpi_decompose_computational_domain() ! --------------------

        print '(A)', 'Parallel computational domain decomposition is '// &
            'inconsistent with use of MPI proxy module for '// &
            'serial execution. Exiting ...'
        stop

    end subroutine s_mpi_decompose_computational_domain ! ------------------

    !> Refer to the m_mpi_proxy.f90 modules in the simulation
                !!                           and/or post-process code directories for information
                !!                           about this subroutine
    subroutine s_mpi_sendrecv_grid_vars_buffer_regions(pbc_loc, sweep_coord)

        character(LEN=*), intent(IN) :: pbc_loc
        character, intent(IN) :: sweep_coord

        print '(A)', 'Communication of buffer regions for grid '// &
            'variables is inconsistent with use of MPI '// &
            'proxy module for serial execution. Exiting ...'
        stop

    end subroutine s_mpi_sendrecv_grid_vars_buffer_regions ! ---------------

    !> Refer to the m_mpi_proxy.f90 modules in the simulation
                !!                           and/or post-process code directories for information
                !!                           about this subroutine
    subroutine s_mpi_reduce_maxloc(loc_var) ! ------------------------------

        real(kind(0d0)), dimension(:), intent(IN) :: loc_var

        print '(A)', 'Multiprocessor maximum location reduction '// &
            'procedure is inconsistent with use of MPI proxy '// &
            'module for serial execution. Exiting ...'
        stop

    end subroutine s_mpi_reduce_maxloc ! -----------------------------------

    !> Refer to the m_mpi_proxy.f90 modules in the simulation
                !!                           and/or post-process code directories for information
                !!                           about this subroutine
    subroutine s_mpi_reduce_min(loc_var) ! ---------------------------------

        real(kind(0d0)), intent(IN) :: loc_var

        print '(A)', 'Multiprocessor minimum reduction procedure is '// &
            'inconsistent with use of MPI proxy module for '// &
            'serial execution. Exiting ...'
        stop

    end subroutine s_mpi_reduce_min ! --------------------------------------

    !> Refer to the m_mpi_proxy.f90 modules in the simulation
                !!                           and/or post-process code directories for information
                !!                           about this subroutine
    subroutine s_mpi_sendrecv_cons_vars_buffer_regions(q_cons_vp, pbc_loc, &
                                                       sweep_coord)

        type(field_position), dimension(:), intent(IN) :: q_cons_vp
        character(LEN=*), intent(IN) :: pbc_loc
        character, intent(IN) :: sweep_coord

        print '(A)', 'Communication of buffer regions for conservative '// &
            'variables is inconsistent with use of MPI proxy '// &
            'module for serial execution. Exiting ...'
        stop

    end subroutine s_mpi_sendrecv_cons_vars_buffer_regions ! ---------------

    !> Refer to the m_mpi_proxy.f90 modules in the simulation
                !!                           and/or post-process code directories for information
                !!                           about this subroutine
    subroutine s_mpi_gather_spatial_extents(spatial_extents) ! -------------

        real(kind(0d0)), dimension(:, :), intent(IN) :: spatial_extents

        print '(A)', 'Gathering of the grids spatial extents is '// &
            'inconsistent with use of MPI proxy module for '// &
            'serial execution. Exiting ...'
        stop

    end subroutine s_mpi_gather_spatial_extents ! --------------------------

    !> Refer to the m_mpi_proxy.f90 modules in the simulation
                !!                           and/or post-process code directories for information
                !!                           about this subroutine
    subroutine s_mpi_defragment_1d_grid_variable() ! -----------------------

        print '(A)', 'Defragmentation of 1D grid variable is '// &
            'inconsistent with use of MPI proxy module for '// &
            'serial execution. Exiting ...'
        stop

    end subroutine s_mpi_defragment_1d_grid_variable ! ---------------------

    !> Refer to the m_mpi_proxy.f90 modules in the simulation
                !!                           and/or post-process code directories for information
                !!                           about this subroutine
    subroutine s_mpi_gather_data_extents(q_fp, data_extents) ! -------------

        real(kind(0d0)), dimension(:, :, :), intent(IN) :: q_fp
        real(kind(0d0)), dimension(:, :), intent(IN) :: data_extents

        print '(A)', 'Gathering of a flow variables extents is '// &
            'inconsistent with use of MPI proxy module for '// &
            'serial execution. Exiting ...'
        stop

    end subroutine s_mpi_gather_data_extents ! -----------------------------

    !> Refer to the m_mpi_proxy.f90 modules in the simulation
                !!                           and/or post-process code directories for information
                !!                           about this subroutine
    subroutine s_mpi_defragment_1d_flow_variable(q_fp, q_root_fp) ! --------

        real(kind(0d0)), dimension(:, :, :), intent(IN) :: q_fp
        real(kind(0d0)), dimension(:, :, :), intent(IN) :: q_root_fp

        print '(A)', 'Defragmentation of 1D flow variable is '// &
            'inconsistent with use of MPI proxy module for '// &
            'serial execution. Exiting ...'
        stop

    end subroutine s_mpi_defragment_1d_flow_variable ! ---------------------

    !> Refer to the m_mpi_proxy.f90 modules in the simulation
                !!                           and/or post-process code directories for information
                !!                           about this subroutine
    subroutine s_deallocate_mpi_proxy_module() ! ---------------------------

        print '(A)', 'Requesting module deallocation procedures is '// &
            'inconsistent with use of MPI proxy module for '// &
            'serial execution. Exiting ...'
        stop

    end subroutine s_deallocate_mpi_proxy_module ! -------------------------

    !> Refer to the m_mpi_proxy.f90 modules in the simulation
                !!                           and/or post-process code directories for information
                !!                           about this subroutine
    subroutine s_finalize_mpi_environment() ! ------------------------------

        return

    end subroutine s_finalize_mpi_environment ! ----------------------------

end module m_mpi_proxy
