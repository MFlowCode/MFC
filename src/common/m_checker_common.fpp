!>
!!@file m_checker_common.f90
!!@brief Contains module m_checker_common

#:include 'macros.fpp'

!> @brief The purpose of the module is to check for compatible input files for.
!!              inputs common to pre-processing, post-processing and simulation
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

    end subroutine s_check_inputs_common

#ifndef MFC_SIMULATION

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

#ifndef MFC_POST_PROCESS

#endif

end module m_checker_common
