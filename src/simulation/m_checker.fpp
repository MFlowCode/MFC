!>
!!@file
!!@brief Contains module m_checker

#:include 'macros.fpp'
#:include 'case.fpp'

!> @brief Validates simulation input parameters for consistency and supported configurations
module m_checker

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_helper

    use m_helper_basic         !< Functions to compare floating point numbers

    implicit none

    private; public :: s_check_inputs

contains

    !> Checks compatibility of parameters in the input file.
        !! Used by the simulation stage
    impure subroutine s_check_inputs

        call s_check_inputs_compilers

        if (igr) then
            call s_check_inputs_nvidia_uvm
        else
            if (recon_type == WENO_TYPE) then
                call s_check_inputs_weno
            else if (recon_type == MUSCL_TYPE) then
                call s_check_inputs_muscl
            end if
        end if

        call s_check_inputs_time_stepping

        call s_check_inputs_hypo_branch

        @:PROHIBIT(ib_state_wrt .and. .not. ib, "ib_state_wrt requires ib to be enabled")

    end subroutine s_check_inputs

    !> Checks constraints on compiler options
    impure subroutine s_check_inputs_compilers
#if !defined(MFC_OpenACC) && !(defined(__PGI) || defined(_CRAYFTN))
        @:PROHIBIT(rdma_mpi, "Unsupported value of rdma_mpi for the current compiler")
#endif
    end subroutine s_check_inputs_compilers

    !> Checks constraints on WENO scheme parameters
    impure subroutine s_check_inputs_weno
        character(len=5) :: numStr !< for int to string conversion

        call s_int_to_str(num_stcls_min*weno_order, numStr)
        @:PROHIBIT(m + 1 < num_stcls_min*weno_order, &
            "m must be greater than or equal to (num_stcls_min*weno_order - 1), whose value is "//trim(numStr))
        @:PROHIBIT(n + 1 < min(1, n)*num_stcls_min*weno_order, &
            "For 2D simulation, n must be greater than or equal to (num_stcls_min*weno_order - 1), whose value is "//trim(numStr))
        @:PROHIBIT(p + 1 < min(1, p)*num_stcls_min*weno_order, &
            "For 3D simulation, p must be greater than or equal to (num_stcls_min*weno_order - 1), whose value is "//trim(numStr))
    end subroutine s_check_inputs_weno

    !> @brief Validates that the grid resolution is sufficient for the MUSCL reconstruction order.
    impure subroutine s_check_inputs_muscl
        character(len=5) :: numStr !< for int to string conversion

        call s_int_to_str(num_stcls_min*muscl_order, numStr)
        @:PROHIBIT(m + 1 < num_stcls_min*muscl_order, &
            "m must be greater than or equal to (num_stcls_min*muscl_order - 1), whose value is "//trim(numStr))
        @:PROHIBIT(n + 1 < min(1, n)*num_stcls_min*muscl_order, &
            "For 2D simulation, n must be greater than or equal to (num_stcls_min*muscl_order - 1), whose value is "//trim(numStr))
        @:PROHIBIT(p + 1 < min(1, p)*num_stcls_min*muscl_order, &
            "For 3D simulation, p must be greater than or equal to (num_stcls_min*muscl_order - 1), whose value is "//trim(numStr))
    end subroutine s_check_inputs_muscl

    !> Checks constraints on time stepping parameters
    impure subroutine s_check_inputs_time_stepping
        if (.not. cfl_dt) then
            @:PROHIBIT(dt <= 0)
        end if
    end subroutine s_check_inputs_time_stepping

    impure subroutine s_check_inputs_nvidia_uvm
#ifdef __NVCOMPILER_GPU_UNIFIED_MEM
        @:PROHIBIT(nv_uvm_igr_temps_on_gpu > 3 .or. nv_uvm_igr_temps_on_gpu < 0, &
            "nv_uvm_igr_temps_on_gpu must be in the range [0, 3]")
        @:PROHIBIT(nv_uvm_igr_temps_on_gpu == 3 .and. igr_iter_solver == 2, &
            "nv_uvm_igr_temps_on_gpu must be in the range [0, 2] for igr_iter_solver == 2")
#endif
    end subroutine s_check_inputs_nvidia_uvm

    impure subroutine s_check_inputs_hypo_branch
        @:PROHIBIT((riemann_solver == 1) .and. hll_u_interface .and. cyl_coord .and. p > 0, &
            "HLL Method 2 is not supported for 3D cylindrical geometry")
        @:PROHIBIT(hypoelasticity .and. riemann_solver == 3, &
            "Exact Riemann (riemann_solver = 3) is not supported with hypoelasticity")
        @:PROHIBIT(hypoelasticity .and. riemann_solver == 2 .and. cyl_coord .and. p > 0, &
            "3D cylindrical hypoelastic HLLC is not supported")
        @:PROHIBIT(hypoelasticity .and. riemann_solver == 4 .and. n == 0, &
            "HLLD hypoelasticity requires at least 2D (n must be > 0)")
        @:PROHIBIT(hypoelasticity .and. riemann_solver == 4 .and. cyl_coord .and. p > 0, &
            "3D cylindrical hypoelastic HLLD is not supported")
        @:PROHIBIT(hypoelasticity .and. riemann_solver == 4 .and. num_fluids /= 2, &
            "HLLD hypoelasticity currently requires exactly 2 fluid components")
        @:PROHIBIT(riemann_solver == 4 .and. (.not. mhd) .and. (.not. hypoelasticity), &
            "HLLD is only available for MHD or hypoelasticity")
        @:PROHIBIT(riemann_hypo_ADC .and. .not. hypoelasticity, &
            "riemann_hypo_ADC requires hypoelasticity = T")
        @:PROHIBIT(riemann_hypo_ADC .and. riemann_solver /= 2 .and. riemann_solver /= 4, &
            "riemann_hypo_ADC only applies to hypo HLLC/HLLD")
        @:PROHIBIT(hypo_hll_interface_rhs .and. .not. hypoelasticity, &
            "hypo_hll_interface_rhs requires hypoelasticity = T")
        @:PROHIBIT(hypo_hll_interface_rhs .and. riemann_solver /= 1, &
            "hypo_hll_interface_rhs requires HLL Riemann solver (riemann_solver = 1)")
        @:PROHIBIT(hypo_hll_interface_rhs .and. cyl_coord .and. p > 0, &
            "3D cylindrical interface-consistent hypo RHS is not supported")
        @:PROHIBIT(alt_soundspeed .and. riemann_solver == 4 .and. .not. hypoelasticity, &
            "alt_soundspeed with HLLD requires hypoelasticity = T")
        @:PROHIBIT(alt_soundspeed .and. riemann_solver == 4 .and. num_fluids /= 2, &
            "alt_soundspeed with HLLD requires exactly 2 fluid components")
        @:PROHIBIT(alt_soundspeed .and. riemann_solver == 1 .and. (.not. hll_u_interface) .and. cyl_coord .and. p == 0, &
            "alt_soundspeed with HLL Method 1 is not supported for 2D axisymmetric geometry")
        @:PROHIBIT(alt_soundspeed .and. riemann_solver == 1 .and. cyl_coord .and. p > 0, &
            "alt_soundspeed with HLL is not currently supported for 3D cylindrical geometry")
    end subroutine s_check_inputs_hypo_branch

end module m_checker
