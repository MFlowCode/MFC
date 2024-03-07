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
              s_calculate_atransposea, &
              s_calculate_elastic_energy
    
    contains 

    subroutine s_allocate_tensor(q_cons_vf,j,k,l,tensor)

        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        integer, intent(IN) :: j, k, l
        real(kind(0d0)), dimension(num_dims**2), intent(OUT) :: tensor

        integer :: i !< Generic loop iterators

        ! Converting the primitive variables to the conservative variables
        do i = 1, num_dims**2
           tensor(i) = q_cons_vf(stress_idx%beg+i-1)%sf(j,k,l)
        end do
    end subroutine s_allocate_tensor

    function f_determinant(tensor)
        real(kind(0d0)), dimension(num_dims**2), intent(IN) :: tensor
        real(kind(0d0)) :: f_determinant

        if (num_dims .eq. 1) then
           f_determinant = tensor(1) ! TODO: Mirelys: does this make sense?
        elseif (num_dims .eq. 2) then
           f_determinant = tensor(1)*tensor(4) - tensor(2)*tensor(3)
        else 
           f_determinant = tensor(1)*(tensor(5)*tensor(9) - tensor(6)*tensor(8)) & 
                           - tensor(2)*(tensor(4)*tensor(9) - tensor(6)*tensor(7)) & 
                           + tensor(3)*(tensor(4)*tensor(8) - tensor(5)*tensor(7))
        end if        
   
    end function f_determinant

    subroutine s_calculate_deviatoric(tensor,deviatoric)
        real(kind(0d0)), dimension(num_dims**2), intent(IN) :: tensor
        real(kind(0d0)), dimension(num_dims**2), intent(OUT) :: deviatoric
        real(kind(0d0)) :: trace13

        deviatoric = tensor
        trace13 = f_trace(tensor)
        trace13 = (1.0/3.0)*trace13
        deviatoric(1) = tensor(1) - trace13
        if (num_dims .eq. 2) then       
           deviatoric(4) = tensor(4) - trace13
        elseif (num_dims .eq. 3) then 
           deviatoric(5) = tensor(5) - trace13
           deviatoric(9) = tensor(9) - trace13
        end if
    end subroutine s_calculate_deviatoric

    subroutine s_calculate_atransposea(tensor,ata)
        real(kind(0d0)), dimension(num_dims**2), intent(IN) :: tensor
        real(kind(0d0)), dimension(num_dims**2), intent(OUT) :: ata

        ata(1) = tensor(1)**2 ! TODO: Mirelys: Does this make sense? 
        if (num_dims .eq. 2) then
           ata(1) = ata(1) + tensor(3)**2
           ata(2) = tensor(1)*tensor(2) + tensor(3)*tensor(4)
           ata(3) = ata(2)
           ata(4) = tensor(2)**2 + tensor(4)**2
        elseif (num_dims .eq. 3) then
           ata(1) = ata(1) + tensor(4)**2 + tensor(7)**2
           ata(5) = tensor(2) + tensor(5)**2 + tensor(8)**2
           ata(9) = tensor(3) + tensor(6)**2 + tensor(9)**2
           ata(2) = tensor(1)*tensor(2) + tensor(4)*tensor(5) + tensor(7)*tensor(8) 
           ata(3) = tensor(1)*tensor(3) + tensor(4)*tensor(6) + tensor(7)*tensor(9) 
           ata(6) = tensor(2)*tensor(3) + tensor(5)*tensor(6) + tensor(8)*tensor(9) 
           ata(4) = ata(2)
           ata(7) = ata(3)
           ata(8) = ata(4)
        end if
    end subroutine s_calculate_atransposea

    function f_elastic_energy(ghat)
        real(kind(0d0)), dimension(num_dims**2), intent(IN) :: ghat
        real(kind(0d0)), dimension(num_dims**2) :: matrix
        real(kind(0d0)) :: f_elastic_energy
        ! remove the identity and square
        matrix = ghat
        matrix(1) = matrix(1) - 1d0
        if (num_dims .eq. 2) then
            matrix(4) = matrix(4) - 1d0
        elseif (num_dims .eq. 3) then 
            matrix(5) = matrix(5) - 1d0
            matrix(9) = matrix(9) - 1d0
        end if 
        matrix(:) = matrix(:)**2
        ! compute the trace
        f_elastic_energy = f_trace(matrix)
    end subroutine s_calculate_elastic_energy

    function f_trace(tensor)
        real(kind(0d0)), dimension(num_dims**2), intent(IN) :: tensor
        real(kind(0d0)) :: f_trace
        f_trace = tensor(1)
        if (num_dims .eq. 2) then
           f_trace = f_trace + tensor(4)
        elseif (num_dims .eq. 3) then
           f_trace = f_trace + tensor(5) + tensor(9)    
        end if
    end function f_trace

end module m_finger_tensor_calc
