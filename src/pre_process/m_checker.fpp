!>
!!@file m_checker.f90
!!@brief Contains module m_checker

#:include 'macros.fpp'

!> @brief The purpose of the module is to check for compatible input files
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

        call s_check_parallel_io

    end subroutine s_check_inputs

    !> Checks if mpi is enabled with parallel_io
    impure subroutine s_check_parallel_io
#ifndef MFC_MPI
        @:PROHIBIT(parallel_io, "MFC built with --no-mpi requires parallel_io=F")
#endif
    end subroutine s_check_parallel_io

end module m_checker
