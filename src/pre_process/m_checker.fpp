!>
!!@file
!!@brief Contains module m_checker

#:include 'macros.fpp'

!> @brief Checks pre-process input file parameters for compatibility and correctness
module m_checker

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_helper_basic         !< Functions to compare floating point numbers

    use m_helper

    implicit none

    private; public :: s_check_inputs

contains

    !> Checks compatibility of parameters in the input file.
        !! Used by the pre_process stage
    impure subroutine s_check_inputs
    end subroutine s_check_inputs

end module m_checker
