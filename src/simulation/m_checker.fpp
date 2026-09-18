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
    use m_constants, only: recon_type_weno, recon_type_muscl, eos_stiffened_gas, eos_ideal_gas, model_eqns_5eq, model_eqns_6eq

    implicit none

    private; public :: s_check_inputs

contains

    !> Checks compatibility of parameters in the input file. Used by the simulation stage
    impure subroutine s_check_inputs

        call s_check_inputs_compilers
        call s_check_inputs_conduction

        if (igr) then
            call s_check_inputs_nvidia_uvm
        else
            if (recon_type == recon_type_weno) then
                call s_check_inputs_weno
            else if (recon_type == recon_type_muscl) then
                call s_check_inputs_muscl
            end if
        end if

        if (ib .and. chemistry) then
            call s_check_inputs_ib_injection
        end if

    end subroutine s_check_inputs

    !> Checks constraints on compiler options
    impure subroutine s_check_inputs_compilers

#if !defined(MFC_OpenACC) && !(defined(__PGI) || defined(_CRAYFTN))
        @:PROHIBIT(rdma_mpi, "Unsupported value of rdma_mpi for the current compiler")
#endif

    end subroutine s_check_inputs_compilers

    !> Checks constraints on Fourier heat conduction inputs
    impure subroutine s_check_inputs_conduction

        integer :: i
        logical :: conducts

        ! Recomputed locally rather than read from the module-level heat_conduction: this runs (via
        ! s_check_input_file) before s_initialize_eqn_idx sets that variable, so it would still hold
        ! its pre-initialization value here. fluid_pp and num_fluids are already populated by the
        ! namelist read that precedes this call.

        conducts = .false.

        do i = 1, num_fluids
            ! lint: runtime-check mirrored in check_heat_conduction (case_validator.py), but duplicated here
            ! deliberately: this binary can be invoked directly against a hand-edited .inp, bypassing that
            ! gate entirely, and a negative k_therm flips conduction into anti-diffusion rather than erroring.
            @:PROHIBIT(fluid_pp(i)%k_therm < 0._wp, "fluid_pp(i)%k_therm must be non-negative")
            if (fluid_pp(i)%k_therm > 0._wp) then
                conducts = .true.
                ! lint: runtime-check same defense as above; cv <= 0 makes f_mixture_temperature divide by
                ! (a floor of) zero, producing a silently wrong temperature rather than a diagnosed error.
                @:PROHIBIT(fluid_pp(i)%cv <= 0._wp, &
                           & "fluid_pp(i)%cv must be positive when fluid_pp(i)%k_therm is set: the mixture temperature is undefined without it")
                ! lint: runtime-check same defense as above; the linear mixture-temperature closure this
                ! module uses does not hold for Mie-Gruneisen/JWL/Vinet fluids, so an unsupported EOS would
                ! silently feed a wrong temperature into the conduction flux.
                @:PROHIBIT(fluid_pp(i)%eos /= eos_stiffened_gas .and. fluid_pp(i)%eos /= eos_ideal_gas, &
                           & "heat conduction supports only the stiffened-gas and ideal-gas equations of state")
                ! lint: runtime-check same defense as above -- and here the failure mode is worse than the
                ! other three: model_eqns = 1 (gamma law) stores gamma/pi_inf, not a volume fraction, at
                ! eqn_idx%adv%beg:end (see m_global_parameters_common.fpp). m_conduction.fpp reads that range
                ! as alpha_i and clamps it to [0, 1], so under model_eqns = 1 it would silently mix gamma/pi_inf
                ! into the conductivity weighting instead of erroring -- the clamp masks the bug rather than
                ! surfacing it.
                @:PROHIBIT(model_eqns /= model_eqns_5eq .and. model_eqns /= model_eqns_6eq, &
                           & "heat conduction requires model_eqns = 2 (5-equation) or model_eqns = 3 (6-equation): fluid_pp(i)%k_therm is weighted by a volume fraction that model_eqns = 1 does not carry")
            end if
        end do

        ! lint: runtime-check load-bearing, not cosmetic: q_T_sf%sf is allocated only inside an
        ! "if (.not. igr)" block in m_time_steppers.fpp but deallocated unconditionally, so heat_conduction
        ! .and. igr would deallocate an unallocated field -- a memory-safety bug the toolchain gate alone
        ! cannot prevent if this binary is invoked directly.
        @:PROHIBIT(conducts .and. igr, "heat conduction is not supported with igr")
        ! lint: runtime-check load-bearing, not cosmetic: with chemistry, m_rhs.fpp allocates the energy
        ! flux_src slot under chemistry .and. chem_params%diffusion .and. .not. viscous, which the conduction
        ! path also allocates when heat_conduction is on, so the combination double-allocates and aborts in
        ! the allocator with no case-level diagnostic. Chemistry also carries its own mixture-averaged
        ! conduction, so the physics would double-count as well.
        @:PROHIBIT(conducts .and. chemistry, &
                   & "heat conduction is not supported with chemistry: the reacting path already carries mixture-averaged conduction through chem_params%diffusion")

    end subroutine s_check_inputs_conduction

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

    !> Validates that each burning immersed-boundary patch injects a species index within the mechanism. inj_species indexes the
    !! image-point mass-fraction array Ys_IP(1:num_species) in m_ibm; an out-of-range value is an out-of-bounds write (silent
    !! corruption). Only reachable with chemistry.
    impure subroutine s_check_inputs_ib_injection

        integer :: i

        do i = 1, num_ibs
            @:PROHIBIT(patch_ib(i)%inj_species > num_species, &
                       & "patch_ib inj_species must be <= num_species (it indexes the image-point species mass fractions; an out-of-range value writes out of bounds)")
        end do

    end subroutine s_check_inputs_ib_injection

end module m_checker
