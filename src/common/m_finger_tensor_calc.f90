!>
!! @file m_variables_conversion.f90
!! @brief Contains module m_variables_conversion

#:include 'macros.fpp'
#:include 'inline_conversions.fpp'
#:include '../simulation/include/case.fpp'

!> @brief This module consists of subroutines used in the calculation of matrix 
!!              operations for the finger tensor

module m_finger_tensor_calc

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_helper
    ! ==========================================================================

    implicit none

    private; 

    public :: s_finger_tensor !variables ! name public variables for all of the subroutines

    contains 

    subroutine s_allocate_tensor(q_cons_vf,j,k,l,tensor)

        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        type(int_bounds_info), optional, intent(IN) :: j, k, l
        real(kind(0d0)), dimension(num_dims**2), intent(OUT) :: tensor

        integer :: i !< Generic loop iterators

#ifndef MFC_SIMULATION
        ! Converting the primitive variables to the conservative variables
        do i = 1, num_dims**2
           tensor(i) = q_cons_vf(stress_idx%beg+i-1)%sf(j,k,l)
        end do
#endif
    end subroutine s_allocate_tensor

    function s_calculate_determinant(tensor)
        real(kind(0d0)), dimension(num_dims**2), intent(OUT) :: tensor
        real(kind(0d0)) :: det

        det = 
        
        return det 

    end function s_calculate_determinant

    subroutine s_calculate_deviatoric(tensor,deviatoric)
        real(kind(0d0)), dimension(num_dims**2), intent(OUT) :: tensor

    end subroutine s_calculate_deviatoric

    subroutine s_calculate_atransposea(tensor,tproduct)
        real(kind(0d0)), dimension(num_dims**2), intent(OUT) :: tensor

    end subroutine s_calculate_atransposea

end module m_finger_tensor_calc
