!>
!!@file
!!@brief Contains module m_checker_common

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief Shared input validation checks for grid dimensions and AMD GPU compiler limits
module m_checker_common

    use m_global_parameters
    use m_mpi_proxy
    use m_helper_basic
    use m_helper

    implicit none

    private; public :: s_check_inputs_common

contains

    !> Checks compatibility of parameters in the input file. Used by all three stages
    impure subroutine s_check_inputs_common

#ifndef MFC_SIMULATION
        call s_check_total_cells
#endif
        call s_check_eos
        #:if USING_AMD
            call s_check_amd
        #:endif

    end subroutine s_check_inputs_common

    !> Restrict the per-fluid EOS selector to the currently supported adapters. Only stiffened_gas (non-chemistry) and
    !! ideal_gas_mixture (chemistry) are backed by a thermodynamics backend, and a single run uses one family for every fluid, so
    !! the unimplemented values and any intra-cell EOS mixing are both rejected here.
    impure subroutine s_check_eos

        integer :: i

        do i = 1, num_fluids
            @:PROHIBIT(chemistry .and. fluid_pp(i)%eos /= eos_ideal_gas_mixture, &
                       & "fluid_pp(:)%eos must be 'ideal_gas_mixture' for every fluid when chemistry is enabled")
            @:PROHIBIT(.not. chemistry .and. fluid_pp(i)%eos /= eos_stiffened_gas, &
                       & "fluid_pp(:)%eos selector is not supported; only 'stiffened_gas' is available " &
                       & // "(or 'ideal_gas_mixture' with a chemistry build)")
        end do

    end subroutine s_check_eos

#ifndef MFC_SIMULATION
    !> Verify that the total number of grid cells meets the minimum required by the number of dimensions and MPI ranks.
    impure subroutine s_check_total_cells

        character(len=18) :: numStr  !< for int to string conversion
        integer(kind=8)   :: min_cells

        min_cells = int(2, kind=8)**int(min(1, m) + min(1, n) + min(1, p), kind=8)*int(num_procs, kind=8)
        call s_int_to_str(2**(min(1, m) + min(1, n) + min(1, p))*num_procs, numStr)

        @:PROHIBIT(nGlobal < min_cells, &
                   & "Total number of cells must be at least (2^[number of dimensions])*num_procs, " // "which is currently " &
                   & // trim(numStr))

    end subroutine s_check_total_cells
#endif

    !> Check that simulation parameters stay within AMD GPU compiler limits when case optimization is disabled.
    impure subroutine s_check_amd

        #:if not MFC_CASE_OPTIMIZATION
            @:PROHIBIT(num_fluids > 3, "num_fluids <= 3 for AMDFLang when Case optimization is off")
            @:PROHIBIT((bubbles_euler .or. bubbles_lagrange) .and. nb > 3, "nb <= 3 for AMDFLang when Case optimization is off")
            @:PROHIBIT(chemistry .and. num_species > 10, "num_species > 10 for AMDFLang when Case optimization is off")
        #:endif

    end subroutine s_check_amd

#ifndef MFC_POST_PROCESS
#endif
end module m_checker_common
