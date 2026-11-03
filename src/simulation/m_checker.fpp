!>
!!@file
!!@brief Contains module m_checker

#:include 'macros.fpp'
#:include 'case.fpp'

!> @brief Validates simulation input parameters for consistency and supported configurations
module m_checker

    use m_global_parameters
    use m_mpi_proxy
    use m_helper
    use m_constants, only: recon_type_weno, recon_type_muscl, time_stepper_rk3, BC_RIEMANN_EXTRAP, BC_CHAR_SLIP_WALL, &
        & BC_CHAR_SUP_OUTFLOW

    implicit none

    private; public :: s_check_inputs

contains

    !> Checks compatibility of parameters in the input file. Used by the simulation stage
    impure subroutine s_check_inputs

        call s_check_inputs_compilers

        if (igr) then
            call s_check_inputs_nvidia_uvm
        else
            if (recon_type == recon_type_weno) then
                call s_check_inputs_weno
            else if (recon_type == recon_type_muscl) then
                call s_check_inputs_muscl
            end if
        end if

        ! lint: runtime-check num_procs is the MPI rank count, not a case-file value

        if (amr) then
            ! Other unsupported combinations are gated in case_validator check_amr. IGR needs np = 1: the fine sigma seed is
            ! injected from the owner's local coarse jac (not P2P-gathered), so a block crossing a rank boundary would read clamped
            ! edge values.
            ! lint: runtime-check num_procs is the MPI rank count, not a case-file value
            @:PROHIBIT(igr .and. num_procs > 1, &
                       & "amr with the IGR solver is only supported at num_procs = 1: the fine block's sigma (jac) Dirichlet seed is injected from the owner's LOCAL coarse jac (not P2P-gathered like q_cons), so a block crossing a rank boundary would read clamped edge values")
            ! lint: runtime-check num_procs is the MPI rank count, not a case-file value
            @:PROHIBIT(amr_max_level > 1 .and. ib .and. num_procs > 1, &
                       & "multi-level AMR (amr_max_level > 1) with immersed boundaries is only supported at num_procs = 1 (the fine-IB image-point stencil is not decomposition-exact across a rank seam)")
        end if

        if (ib) then
            call s_check_inputs_ib_injection
        end if

    end subroutine s_check_inputs

    !> Checks constraints on compiler options
    impure subroutine s_check_inputs_compilers

#if !defined(MFC_OpenACC) && !defined(MFC_OpenMP) && !(defined(__PGI) || defined(_CRAYFTN))
        @:PROHIBIT(rdma_mpi, "Unsupported value of rdma_mpi for the current compiler")
#endif

    end subroutine s_check_inputs_compilers

    !> Checks constraints on WENO scheme parameters
    impure subroutine s_check_inputs_weno

        character(len=5) :: numStr  !< for int to string conversion

        call s_int_to_str(num_stcls_min*weno_order, numStr)
        ! lint: runtime-check m/n/p are per-rank extents after MPI decomposition, not the case-file values
        @:PROHIBIT(m + 1 < num_stcls_min*weno_order, &
                   & "m must be greater than or equal to (num_stcls_min*weno_order - 1), whose value is " // trim(numStr))
        ! lint: runtime-check per-rank n
        @:PROHIBIT(n + 1 < min(1, n)*num_stcls_min*weno_order, &
                   & "For 2D simulation, n must be greater than or equal to (num_stcls_min*weno_order - 1), whose value is " &
                   & // trim(numStr))
        ! lint: runtime-check per-rank p
        @:PROHIBIT(p + 1 < min(1, p)*num_stcls_min*weno_order, &
                   & "For 3D simulation, p must be greater than or equal to (num_stcls_min*weno_order - 1), whose value is " &
                   & // trim(numStr))

    end subroutine s_check_inputs_weno

    !> Validate that the grid resolution is sufficient for the MUSCL reconstruction order
    impure subroutine s_check_inputs_muscl

        character(len=5) :: numStr  !< for int to string conversion

        call s_int_to_str(num_stcls_min*muscl_order, numStr)
        ! lint: runtime-check m/n/p are per-rank extents after MPI decomposition, not the case-file values
        @:PROHIBIT(m + 1 < num_stcls_min*muscl_order, &
                   & "m must be greater than or equal to (num_stcls_min*muscl_order - 1), whose value is " // trim(numStr))
        ! lint: runtime-check per-rank n
        @:PROHIBIT(n + 1 < min(1, n)*num_stcls_min*muscl_order, &
                   & "For 2D simulation, n must be greater than or equal to (num_stcls_min*muscl_order - 1), whose value is " &
                   & // trim(numStr))
        ! lint: runtime-check per-rank p
        @:PROHIBIT(p + 1 < min(1, p)*num_stcls_min*muscl_order, &
                   & "For 3D simulation, p must be greater than or equal to (num_stcls_min*muscl_order - 1), whose value is " &
                   & // trim(numStr))

    end subroutine s_check_inputs_muscl

    !> Validate NVIDIA unified virtual memory configuration parameters
    impure subroutine s_check_inputs_nvidia_uvm

#ifdef __NVCOMPILER_GPU_UNIFIED_MEM
        @:PROHIBIT(nv_uvm_igr_temps_on_gpu > 3 .or. nv_uvm_igr_temps_on_gpu < 0, &
                   & "nv_uvm_igr_temps_on_gpu must be in the range [0, 3]")
        @:PROHIBIT(nv_uvm_igr_temps_on_gpu == 3 .and. igr_iter_solver == 2, &
                   & "nv_uvm_igr_temps_on_gpu must be in the range [0, 2] for igr_iter_solver == 2")
#endif

    end subroutine s_check_inputs_nvidia_uvm

    !> Validates immersed-boundary injection, thermal, and heterogeneous surface-reaction parameters.
    !> Validates the one immersed-boundary surface parameter that cannot be checked from the case file alone.
    !!
    !! Every other relation between these parameters -- the ranges, the chemistry and inj_species
    !! combinations, the Twall window -- is between case-file values and lives in
    !! case_validator.py, where it runs before any binary starts and is unit tested. num_species
    !! is populated by Cantera at run time, so this bound has to be here.
    impure subroutine s_check_inputs_ib_injection

        integer :: i

        do i = 1, num_ibs
            ! lint: runtime-check num_species is populated by Cantera, not read from the case file
            @:PROHIBIT(patch_ib(i)%inj_species > num_species, "patch_ib inj_species must be <= num_species")
        end do

    end subroutine s_check_inputs_ib_injection

end module m_checker
