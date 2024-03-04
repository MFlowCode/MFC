!>
!! @file m_variables_conversion.f90
!! @brief Contains module m_variables_conversion

#:include 'macros.fpp'
#:include 'inline_conversions.fpp'
#:include '../simulation/include/case.fpp'

!> @brief This module consists of subroutines used in the conversion of the
!!              conservative variables into the primitive ones and vice versa. In
!!              addition, the module also contains the subroutines used to obtain
!!              the mixture variables and the subroutines used to compute pressure.
module m_finger_matrix_calc

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_helper
    ! ==========================================================================

    implicit none






end module m_finger_matrix_calc
