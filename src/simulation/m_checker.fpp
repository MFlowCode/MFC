!>
!!@file m_start_up.f90
!!@brief Contains module m_checker

#:include 'macros.fpp'
#:include 'case.fpp'

!> @brief The purpose of the module is to check for compatible input files
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
            call s_check_inputs_geometry_precision
        end if

        call s_check_inputs_time_stepping

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

    !> Checks constraints on geometry and precision
    impure subroutine s_check_inputs_geometry_precision
        ! Prevent spherical geometry in single precision
#ifdef MFC_SINGLE_PRECISION
        @:PROHIBIT(.not. (cyl_coord .neqv. .true. .or. (cyl_coord .and. p == 0)), "Fully 3D cylindrical grid (geometry = 3) is not supported in single precision.")
#endif
    end subroutine s_check_inputs_geometry_precision

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

end module m_checker
