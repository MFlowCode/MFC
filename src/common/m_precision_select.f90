!> @file m_precision_select.f90
!! @brief Contains module m_precision_select

!> @brief This file contains the definition of floating point used in MFC
module m_precision_select
#ifdef MFC_MPI
    use mpi                    !< Message passing interface (MPI) module
#endif

    implicit none

    ! Define the available precision types
    integer, parameter :: single_precision = selected_real_kind(6, 37)
    integer, parameter :: double_precision = selected_real_kind(15, 307)

    ! Set the working precision (wp) to single or double precision
    integer, parameter :: wp = single_precision  ! Change to single_precision if needed

#ifdef MFC_MPI
    ! Set mpi_p based on wp using the merge intrinsic function
    integer, parameter :: mpi_p = merge(MPI_DOUBLE_PRECISION, MPI_REAL, wp == double_precision)
#else
    integer, parameter :: mpi_p = -100  ! Default value when MPI is not used
#endif

end module m_precision_select
