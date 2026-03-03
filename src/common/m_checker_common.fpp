!>
!!@file
!!@brief Contains module m_checker_common

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief Shared input validation checks for grid dimensions and AMD GPU compiler limits
module m_checker_common

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_helper_basic         !< Functions to compare floating point numbers

    use m_helper

    implicit none

    private; public :: s_check_inputs_common, wp

contains

    !> Checks compatibility of parameters in the input file.
        !! Used by all three stages
    impure subroutine s_check_inputs_common

#ifndef MFC_SIMULATION
        call s_check_total_cells
#endif
        #:if USING_AMD
            call s_check_amd
        #:endif

    end subroutine s_check_inputs_common

#ifndef MFC_SIMULATION

    !> @brief Verifies that the total number of grid cells meets the minimum required by the number of dimensions and MPI ranks.
    impure subroutine s_check_total_cells
        character(len=18) :: numStr !< for int to string conversion
        integer(kind=8) :: min_cells

        min_cells = int(2, kind=8)**int(min(1, m) + min(1, n) + min(1, p), kind=8)*int(num_procs, kind=8)
        call s_int_to_str(2**(min(1, m) + min(1, n) + min(1, p))*num_procs, numStr)

        @:PROHIBIT(nGlobal < min_cells, &
            "Total number of cells must be at least (2^[number of dimensions])*num_procs, " // &
            "which is currently "//trim(numStr))
    end subroutine s_check_total_cells

#endif

    !> @brief Checks that simulation parameters stay within AMD GPU compiler limits when case optimization is disabled.
    impure subroutine s_check_amd

        #:if not MFC_CASE_OPTIMIZATION
            @:PROHIBIT(num_fluids > 3, "num_fluids <= 3 for AMDFLang when Case optimization is off")
            @:PROHIBIT((bubbles_euler .or. bubbles_lagrange) .and. nb > 3, "nb <= 3 for AMDFLang when Case optimization is off")
            @:PROHIBIT(chemistry .and. num_species /= 10, "num_species = 10 for AMDFLang when Case optimization is off")
        #:endif

    end subroutine s_check_amd

#ifndef MFC_POST_PROCESS

#endif

end module m_checker_common
