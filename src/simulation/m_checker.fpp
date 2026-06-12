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

    implicit none

    private; public :: s_check_inputs

contains

    !> Checks compatibility of parameters in the input file. Used by the simulation stage
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

        call s_check_inputs_eos

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
        @:PROHIBIT(muscl_order == 1 .and. int_comp > 0, &
                   & "int_comp requires muscl_order >= 2 (muscl_order=1 leaves the reconstruction workspace uninitialised)")

    end subroutine s_check_inputs_muscl

    !> Checks constraints on time stepping parameters
    impure subroutine s_check_inputs_time_stepping

        if (.not. cfl_dt) then
            @:PROHIBIT(dt <= 0)
        end if

    end subroutine s_check_inputs_time_stepping

    !> Checks constraints for EOS selections and EOS/model-equation pairings.
    impure subroutine s_check_inputs_eos

        integer :: i, num_jwl

        num_jwl = 0
        do i = 1, num_fluids
            @:PROHIBIT(fluid_pp(i)%eos /= 1 .and. fluid_pp(i)%eos /= 2, &
                       & "Unsupported fluid_pp(i)%eos. Use 1 for ideal/stiffened gas or 2 for JWL.")
            if (fluid_pp(i)%eos == 1) then
                @:PROHIBIT(fluid_pp(i)%gamma <= 0._wp, "Ideal/stiffened gas EOS requires fluid_pp(i)%gamma > 0.")
                @:PROHIBIT(fluid_pp(i)%pi_inf < 0._wp, "Ideal/stiffened gas EOS requires fluid_pp(i)%pi_inf >= 0.")
            else if (fluid_pp(i)%eos == 2) then
                num_jwl = num_jwl + 1
                @:PROHIBIT(fluid_pp(i)%jwl_A <= 0._wp, "JWL EOS requires fluid_pp(i)%jwl_A > 0.")
                ! jwl_B may be negative for insensitive explosives (e.g. LX-17, PBX-9502)
                @:PROHIBIT(fluid_pp(i)%jwl_R1 <= 0._wp, "JWL EOS requires fluid_pp(i)%jwl_R1 > 0.")
                @:PROHIBIT(fluid_pp(i)%jwl_R2 <= 0._wp, "JWL EOS requires fluid_pp(i)%jwl_R2 > 0.")
                @:PROHIBIT(fluid_pp(i)%jwl_omega <= 0._wp, "JWL EOS requires fluid_pp(i)%jwl_omega > 0.")
                @:PROHIBIT(fluid_pp(i)%jwl_rho0 <= 0._wp, "JWL EOS requires fluid_pp(i)%jwl_rho0 > 0.")
                @:PROHIBIT(fluid_pp(i)%jwl_air_gamma <= 0._wp, "JWL/air closure requires fluid_pp(i)%jwl_air_gamma > 0.")
            end if
        end do

        if (num_jwl > 0) then
            @:PROHIBIT(num_jwl > 1, "JWL EOS currently supports one JWL fluid per case.")
            @:PROHIBIT(jwl_mix_type < 0 .or. jwl_mix_type > 3, &
                       & "jwl_mix_type must be 0 (isobaric), 1 (Kuhl), 2 (p-T equilibrium), or 3 (Rocflu blend).")
            @:PROHIBIT(riemann_solver == 4, "JWL EOS is not supported with the HLLD/MHD flux solver. Use HLL, HLLC, or LF.")
            @:PROHIBIT(model_eqns == 3 .and. .not. relax, &
                       & "JWL six-equation cases require relax=T so the JWL pressure-relaxation closure is applied.")
            @:PROHIBIT(num_fluids > 2, &
                       & "JWL EOS currently supports at most two fluids; the interface-compression path assumes num_fluids = 2.")
        end if

    end subroutine s_check_inputs_eos

    !> Validate NVIDIA unified virtual memory configuration parameters
    impure subroutine s_check_inputs_nvidia_uvm

#ifdef __NVCOMPILER_GPU_UNIFIED_MEM
        @:PROHIBIT(nv_uvm_igr_temps_on_gpu > 3 .or. nv_uvm_igr_temps_on_gpu < 0, &
                   & "nv_uvm_igr_temps_on_gpu must be in the range [0, 3]")
        @:PROHIBIT(nv_uvm_igr_temps_on_gpu == 3 .and. igr_iter_solver == 2, &
                   & "nv_uvm_igr_temps_on_gpu must be in the range [0, 2] for igr_iter_solver == 2")
#endif

    end subroutine s_check_inputs_nvidia_uvm

end module m_checker
