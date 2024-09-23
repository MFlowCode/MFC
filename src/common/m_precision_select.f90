!> @file m_precision_select.f90
!> @brief Contains module m_precision_select

!> @brief This file contains the definition of floating point used in MFC
module m_precision_select
#ifdef MFC_MPI
    use mpi                    !< Message Passing Interface (MPI) module
#endif

    implicit none

    ! Define the available precision types
    integer, parameter :: single_precision = selected_real_kind(6, 37)
    integer, parameter :: double_precision = selected_real_kind(15, 307)

    ! Set the working precision (wp) to single or double precision
    integer, parameter :: wp = double_precision  ! Change this to single_precision if needed

#ifdef MFC_MPI
    ! Declare mpi_p as a module variable
    integer :: mpi_p
#else
    integer, parameter :: mpi_p = -100
#endif

contains

    ! Subroutine to initialize mpi_p based on wp
    subroutine initialize_precision()
#ifdef MFC_MPI
        if (wp == single_precision) then
            mpi_p = MPI_FLOAT
        else if (wp == double_precision) then
            mpi_p = MPI_DOUBLE_PRECISION
        else
            stop 'Unsupported precision kind.'
        end if
#endif
    end subroutine initialize_precision

end module m_precision_select