!>
!! @file m_variables_conversion.f90
!! @brief Contains module m_variables_conversion

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

    private; public :: s_allocate_tensor, & 
              f_determinant, &
              s_calculate_deviatoric, &
              s_calculate_atransposea
    
    contains 

    subroutine s_allocate_tensor(q_cons_vf,j,k,l,tensor)

        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        integer, intent(IN) :: j, k, l
        real(kind(0d0)), dimension(num_dims**2), intent(OUT) :: tensor

        integer :: i !< Generic loop iterators

#ifndef MFC_SIMULATION
        ! Converting the primitive variables to the conservative variables
        do i = 1, num_dims**2
           tensor(i) = q_cons_vf(stress_idx%beg+i-1)%sf(j,k,l)
        end do
#endif
    end subroutine s_allocate_tensor

    function f_determinant(tensor)
        real(kind(0d0)), dimension(num_dims**2), intent(OUT) :: tensor
        real(kind(0d0)) :: f_determinant

        if (num_dims .eq. 1) then
           f_determinant = tensor(1) ! does this make sense?
        elseif (num_dims .eq. 2) then
           f_determinant = tensor(1)*tensor(4) - tensor(2)*tensor(3)
        else 
           f_determinant = tensor(1)*(tensor(2)*tensor(3))
        end if
        
    end function f_determinant

    subroutine s_calculate_deviatoric(tensor,deviatoric)
        real(kind(0d0)), dimension(num_dims**2), intent(IN) :: tensor
        real(kind(0d0)), dimension(num_dims**2), intent(OUT) :: deviatoric

            end subroutine s_calculate_deviatoric

    subroutine s_calculate_atransposea(tensor,tproduct)
        real(kind(0d0)), dimension(num_dims**2), intent(IN) :: tensor
        real(kind(0d0)), dimension(num_dims**2), intent(OUT) :: tproduct

    end subroutine s_calculate_atransposea

end module m_finger_tensor_calc
