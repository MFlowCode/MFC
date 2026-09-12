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
    use m_constants, only: recon_type_weno, recon_type_muscl

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

        if (ib) then
            call s_check_inputs_ib_injection
        end if

    end subroutine s_check_inputs

    !> Checks constraints on compiler options
    impure subroutine s_check_inputs_compilers

#if !defined(MFC_OpenACC) && !(defined(__PGI) || defined(_CRAYFTN))
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
    impure subroutine s_check_inputs_ib_injection

        integer :: i

        do i = 1, num_ibs
            ! Basic parameter ranges
            @:PROHIBIT(patch_ib(i)%inj_species < 0 .or. patch_ib(i)%inj_species > num_species, &
                       & "patch_ib inj_species must be in [0,num_species]")
            @:PROHIBIT(patch_ib(i)%thermal_bc < 0 .or. patch_ib(i)%thermal_bc > 2, "patch_ib thermal_bc must be 0, 1, or 2")
            @:PROHIBIT(patch_ib(i)%surface_reaction < 0 .or. patch_ib(i)%surface_reaction > 1, &
                       & "patch_ib surface_reaction must be 0 or 1")

            ! Thermal immersed-boundary condition
            !    0 = zero-normal-gradient temperature
            !    1 = prescribed wall temperature (Twall)
            !    2 = reacting surface energy balance
            if (patch_ib(i)%thermal_bc == 1) then
                @:PROHIBIT(patch_ib(i)%Twall <= 0._wp, "patch_ib Twall must be > 0 when thermal_bc = 1")
            end if

            if (patch_ib(i)%thermal_bc == 2) then
                @:PROHIBIT(patch_ib(i)%surface_reaction /= 1, "patch_ib thermal_bc = 2 requires surface_reaction = 1")
            end if

            ! Heterogeneous surface reaction    0 = none    1 = enabled
            if (patch_ib(i)%surface_reaction == 1) then
                @:PROHIBIT(.not. chemistry, "patch_ib surface_reaction = 1 requires chemistry = T")
                @:PROHIBIT(patch_ib(i)%inj_species > 0, "patch_ib surface_reaction = 1 cannot be combined with inj_species > 0")
            end if
        end do

    end subroutine s_check_inputs_ib_injection

end module m_checker
