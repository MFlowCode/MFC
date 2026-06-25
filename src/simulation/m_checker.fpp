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
    use m_helper_basic
    use m_constants, only: recon_type_weno, recon_type_muscl, muscl_order_first_order

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

        call s_check_inputs_time_stepping

        call s_check_inputs_hll_non_conservative

        call s_check_inputs_hypo_branch

        @:PROHIBIT(ib_state_wrt .and. .not. ib, "ib_state_wrt requires ib to be enabled")
        @:PROHIBIT(many_ib_patch_parallelism .and. .not. ib, "many_ib_patch_parallelism requires ib to be enabled")

        if (num_particle_clouds > 0) then
            call s_check_inputs_particle_clouds
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
        @:PROHIBIT(m + 1 < num_stcls_min*weno_order, &
                   & "m must be greater than or equal to (num_stcls_min*weno_order - 1), whose value is " // trim(numStr))
        @:PROHIBIT(n + 1 < min(1, n)*num_stcls_min*weno_order, &
                   & "For 2D simulation, n must be greater than or equal to (num_stcls_min*weno_order - 1), whose value is " &
                   & // trim(numStr))
        @:PROHIBIT(p + 1 < min(1, p)*num_stcls_min*weno_order, &
                   & "For 3D simulation, p must be greater than or equal to (num_stcls_min*weno_order - 1), whose value is " &
                   & // trim(numStr))

    end subroutine s_check_inputs_weno

    !> Validate that the grid resolution is sufficient for the MUSCL reconstruction order
    impure subroutine s_check_inputs_muscl

        character(len=5) :: numStr  !< for int to string conversion

        call s_int_to_str(num_stcls_min*muscl_order, numStr)
        @:PROHIBIT(m + 1 < num_stcls_min*muscl_order, &
                   & "m must be greater than or equal to (num_stcls_min*muscl_order - 1), whose value is " // trim(numStr))
        @:PROHIBIT(n + 1 < min(1, n)*num_stcls_min*muscl_order, &
                   & "For 2D simulation, n must be greater than or equal to (num_stcls_min*muscl_order - 1), whose value is " &
                   & // trim(numStr))
        @:PROHIBIT(p + 1 < min(1, p)*num_stcls_min*muscl_order, &
                   & "For 3D simulation, p must be greater than or equal to (num_stcls_min*muscl_order - 1), whose value is " &
                   & // trim(numStr))
        @:PROHIBIT(muscl_order == muscl_order_first_order .and. int_comp > 0, &
                   & "int_comp requires muscl_order >= 2 (muscl_order=1 leaves the reconstruction workspace uninitialised)")

    end subroutine s_check_inputs_muscl

    !> Checks constraints on time stepping parameters
    impure subroutine s_check_inputs_time_stepping

        if (.not. cfl_dt) then
            @:PROHIBIT(dt <= 0)
        end if

    end subroutine s_check_inputs_time_stepping

    !> Validate NVIDIA unified virtual memory configuration parameters
    impure subroutine s_check_inputs_nvidia_uvm

#ifdef __NVCOMPILER_GPU_UNIFIED_MEM
        @:PROHIBIT(nv_uvm_igr_temps_on_gpu > 3 .or. nv_uvm_igr_temps_on_gpu < 0, &
                   & "nv_uvm_igr_temps_on_gpu must be in the range [0, 3]")
        @:PROHIBIT(nv_uvm_igr_temps_on_gpu == 3 .and. igr_iter_solver == 2, &
                   & "nv_uvm_igr_temps_on_gpu must be in the range [0, 2] for igr_iter_solver == 2")
#endif

    end subroutine s_check_inputs_nvidia_uvm

    impure subroutine s_check_inputs_hll_non_conservative

        @:PROHIBIT((riemann_solver == 1) .and. hll_u_interface .and. cyl_coord .and. p > 0, &
                   & "HLL Method 2 is not supported for 3D cylindrical geometry")
        @:PROHIBIT(alt_soundspeed .and. riemann_solver == 1 .and. (.not. hll_u_interface) .and. cyl_coord .and. p == 0, &
                   & "alt_soundspeed with HLL Method 1 is not supported for 2D axisymmetric geometry")
        @:PROHIBIT(alt_soundspeed .and. riemann_solver == 1 .and. cyl_coord .and. p > 0, &
                   & "alt_soundspeed with HLL is not currently supported for 3D cylindrical geometry")

    end subroutine s_check_inputs_hll_non_conservative

    impure subroutine s_check_inputs_hypo_branch

        @:PROHIBIT(hypoelasticity .and. cyl_coord .and. p > 0, "3D cylindrical hypoelasticity is not supported")

        ! Hypoelasticity solver restrictions
        @:PROHIBIT(hypoelasticity .and. riemann_solver == 3, &
                   & "Exact Riemann (riemann_solver = 3) is not supported with hypoelasticity")
        @:PROHIBIT(hypoelasticity .and. riemann_solver == 4 .and. n == 0, &
                   & "HLLD hypoelasticity requires at least 2D (n must be > 0)")
        @:PROHIBIT(hypoelasticity .and. riemann_solver == 4 .and. num_fluids /= 2, &
                   & "HLLD hypoelasticity currently requires exactly 2 fluid components")
        @:PROHIBIT(hypoelasticity .and. (riemann_solver == 1 .or. riemann_solver == 2) .and. num_fluids > 2, &
                   & "HLL/HLLC hypoelasticity supports at most 2 fluid components (blkmod uses fluid_pp(1:2) only)")
        @:PROHIBIT(hypoelasticity .and. riemann_solver == 4 .and. viscous, &
                   & "HLLD hypoelasticity does not support viscous effects (the dual-pass omits the viscous source term)")
        @:PROHIBIT(hypoelasticity .and. riemann_solver == 4 .and. surface_tension, &
                   & "HLLD hypoelasticity does not support surface tension (the dual-pass omits the surface-tension source term)")
        @:PROHIBIT(hypoelasticity .and. riemann_solver == 4 .and. cont_damage, &
                   & "HLLD hypoelasticity does not support continuum damage (the dual-pass does not damage-scale the shear modulus)")
        @:PROHIBIT(hypoelasticity .and. riemann_solver == 4 .and. bubbles_euler, &
                   & "HLLD hypoelasticity does not support Euler-Euler bubbles (the dual-pass omits the bubble source and divergence term)")
        @:PROHIBIT(hypoelasticity .and. riemann_solver == 4 .and. chemistry, "HLLD hypoelasticity does not support chemistry")
        @:PROHIBIT(riemann_solver == 4 .and. (.not. mhd) .and. (.not. hypoelasticity), &
                   & "HLLD is only available for MHD or hypoelasticity")

        ! Feature flag prerequisites
        @:PROHIBIT(riemann_hypo_ADC .and. .not. hypoelasticity, "riemann_hypo_ADC requires hypoelasticity = T")
        @:PROHIBIT(riemann_hypo_ADC .and. riemann_solver /= 2 .and. riemann_solver /= 4, &
                   & "riemann_hypo_ADC only applies to hypo HLLC/HLLD")
        @:PROHIBIT(riemann_hypo_ADC .and. (bubbles_euler .or. surface_tension .or. chemistry .or. cont_damage), &
                   & "riemann_hypo_ADC does not support bubbles, surface tension, chemistry, or continuum damage (the ADC HLL blend omits their flux components)")
        @:PROHIBIT(hypo_hll_interface_rhs .and. .not. hypoelasticity, "hypo_hll_interface_rhs requires hypoelasticity = T")
        @:PROHIBIT(hypo_hll_interface_rhs .and. riemann_solver /= 1, &
                   & "hypo_hll_interface_rhs requires HLL Riemann solver (riemann_solver = 1)")
        @:PROHIBIT(alt_soundspeed .and. riemann_solver == 4 .and. .not. hypoelasticity, &
                   & "alt_soundspeed with HLLD requires hypoelasticity = T")
        @:PROHIBIT(hypoelasticity .and. alt_soundspeed .and. num_fluids /= 2, &
                   & "hypoelastic alt_soundspeed requires exactly 2 fluid components")
        @:PROHIBIT(hypoelasticity .and. igr, "Hypoelasticity is not compatible with IGR")

    end subroutine s_check_inputs_hypo_branch

    !> Checks that each active particle cloud has a valid packing_method specified
    impure subroutine s_check_inputs_particle_clouds

        integer          :: i
        character(len=5) :: idxStr

        do i = 1, num_particle_clouds
            call s_int_to_str(i, idxStr)
            @:PROHIBIT(particle_cloud(i)%packing_method == dflt_int, &
                       & "particle_cloud("//trim(idxStr) &
                       & //")%packing_method must be specified (1 = rejection sampling, 2 = lattice)")
            @:PROHIBIT(particle_cloud(i)%packing_method /= 1 .and. particle_cloud(i)%packing_method /= 2, &
                       & "particle_cloud("//trim(idxStr) //")%packing_method must be 1 (rejection sampling) or 2 (lattice)")
        end do

    end subroutine s_check_inputs_particle_clouds

end module m_checker
