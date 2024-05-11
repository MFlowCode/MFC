!>
!! @file m_variables_conversion.f90
!! @brief Contains module m_variables_conversion

!> @brief This module consists of subroutines used in the calculation of matrix
!!              operations for the reference map tensor

module m_rmt_tensor_calc

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_helper
    ! ==========================================================================

    implicit none

    private; public :: s_calculate_btensor, &
 f_elastic_energy, &
 s_calculate_deviatoric

contains

    subroutine s_calculate_btensor(q_prim_vf, j, k, l, btensor)

        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf
        type(scalar_field), dimension(num_dims*(num_dims+1)/2 + 1), intent(OUT) :: btensor
        integer, intent(IN) :: j, k, l

        real(kind(0d0)), dimension(num_dims**2) :: ftensor, ftransposef, tensorb, tensor
        integer :: i !< Generic loop iterators

        ! Converting the primitive variables to the conservative variables
        do i = 1, num_dims
            tensor(i) = q_prim_vf(stress_idx%beg + i - 1)%sf(j, k, l)
        end do
        ! NOTE: btensor is symmetric, save the data space
        ! need to calculate gradxi then calculate btensor and J = det(F)
        ! store in btensor

        ! extracting the nxn tensor for the calculation
        !do i = 1, num_dims**2
        !    ftensor(i) = gradxitensor(i)%sf(j, k, l)
        !end do
        !call s_calculate_atransposea(ftensor,ftransposef)
        !call s_calculate_ainverse(ftransposef,btensor)
        !jacobian = f_determinant(ftensor)

    end subroutine s_calculate_btensor

    function f_determinant(tensor)
        real(kind(0d0)), dimension(num_dims**2), intent(IN) :: tensor
        real(kind(0d0)) :: f_determinant

        if (num_dims == 1) then
            f_determinant = tensor(1) ! TODO: Mirelys: does this make sense?
        elseif (num_dims == 2) then
            f_determinant = tensor(1)*tensor(4) - tensor(2)*tensor(3)
        else
            f_determinant = tensor(1)*(tensor(5)*tensor(9) - tensor(6)*tensor(8)) &
                            - tensor(2)*(tensor(4)*tensor(9) - tensor(6)*tensor(7)) &
                            + tensor(3)*(tensor(4)*tensor(8) - tensor(5)*tensor(7))
        end if
        ! error checking
        if (f_determinant == 0) then
            print *, 'ERROR: Determinant was zero'
            call s_mpi_abort()
        end if
    end function f_determinant

    subroutine s_calculate_deviatoric(tensor, deviatoric)
        real(kind(0d0)), dimension(num_dims**2), intent(IN) :: tensor
        real(kind(0d0)), dimension(num_dims**2), intent(OUT) :: deviatoric
        real(kind(0d0)) :: trace13
        deviatoric = tensor
        trace13 = f_trace(tensor)
        trace13 = (1.0/3.0)*trace13
        deviatoric(1) = tensor(1) - trace13
        if (num_dims == 2) then
            deviatoric(4) = tensor(4) - trace13
        elseif (num_dims == 3) then
            deviatoric(5) = tensor(5) - trace13
            deviatoric(9) = tensor(9) - trace13
        end if
    end subroutine s_calculate_deviatoric

    function f_trace(tensor)
        real(kind(0d0)), dimension(num_dims**2), intent(IN) :: tensor
        real(kind(0d0)) :: f_trace
        f_trace = tensor(1)
        if (num_dims == 2) then
            f_trace = f_trace + tensor(4)
        elseif (num_dims == 3) then
            f_trace = f_trace + tensor(5) + tensor(9)
        end if
    end function f_trace

    subroutine s_calculate_atransposea(tensor, ata)
        real(kind(0d0)), dimension(num_dims**2), intent(IN) :: tensor
        real(kind(0d0)), dimension(num_dims**2), intent(OUT) :: ata

        ata(1) = tensor(1)**2 ! TODO: Mirelys: Does this make sense?
        if (num_dims == 2) then
            ata(1) = ata(1) + tensor(3)**2
            ata(2) = tensor(1)*tensor(2) + tensor(3)*tensor(4)
            ata(3) = ata(2)
            ata(4) = tensor(2)**2 + tensor(4)**2
        elseif (num_dims == 3) then
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

    subroutine s_calculate_adjointa(tensor, dja)
        real(kind(0d0)), dimension(num_dims**2), intent(IN) :: tensor
        real(kind(0d0)), dimension(num_dims**2), intent(OUT) :: dja

        if (num_dims == 1) then
            dja(1) = 1
        elseif (num_dims == 2) then
            dja(1) = tensor(4)
            dja(2) = -tensor(3)
            dja(3) = -tensor(2)
            dja(4) = tensor(1)
        elseif (num_dims == 3) then
            dja(1) = tensor(5)*tensor(9) - tensor(6)*tensor(8)
            dja(2) = -(tensor(2)*tensor(9) - tensor(3)*tensor(8))
            dja(3) = tensor(2)*tensor(6) - tensor(3)*tensor(5)
            dja(4) = -(tensor(4)*tensor(9) - tensor(6)*tensor(7))
            dja(5) = tensor(1)*tensor(9) - tensor(3)*tensor(7)
            dja(6) = -(tensor(1)*tensor(6) - tensor(4)*tensor(3))
            dja(7) = tensor(4)*tensor(8) - tensor(5)*tensor(7)
            dja(8) = -(tensor(1)*tensor(8) - tensor(2)*tensor(7))
            dja(9) = tensor(1)*tensor(5) - tensor(2)*tensor(4)
        end if
    end subroutine s_calculate_adjointa

    subroutine s_calculate_ainverse(tensor, ainv)
        real(kind(0d0)), dimension(num_dims**2), intent(IN) :: tensor
        real(kind(0d0)), dimension(num_dims**2), intent(OUT) :: ainv
        real(kind(0d0)), dimension(num_dims**2) :: dja
        real(kind(0d0)) :: det

        call s_calculate_adjointa(tensor, dja)
        det = f_determinant(tensor)
        ainv(:) = dja(:)/det
    end subroutine s_calculate_ainverse

    ! neo-Hookean only at this time, will need to be changed later
    function f_elastic_energy(btensor, j, k, l)
        type(scalar_field), & 
            dimension(num_dims*(num_dims+1)/2 + 1), &
            intent(IN) :: btensor

        integer, intent(IN) :: j, k, l

        real(kind(0d0)), dimension(num_dims**2) :: ftransposef, tensorb
        real(kind(0d0)) :: invariant1, jacobian, f_elastic_energy
        integer :: i !< Generic loop iterators

        ! extracting the nxn tensor for the calculation
        !TODO COPY SPRATT CODE FOR SYMMETRIC TENSOR
        do i = 1, num_dims*(num_dims+1)/2
            tensorb(i) = btensor(i)%sf(j, k, l)
        end do
        tensorb(1) = btensor(1)%sf(j, k, l)
        
        jacobian = btensor(num_dims*(num_dims+1)/2 + 1)%sf(j, k, l)
        invariant1 = f_trace(tensorb)
        ! compute the invariant without the elastic modulus
        f_elastic_energy = 0.5d0*(invariant1 - 3)/jacobian
    end function f_elastic_energy


end module m_rmt_tensor_calc
