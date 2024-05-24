!>
!! @file m_hyperelastic.f90
!! @brief Contains module m_hyperelastic

!> @brief This module is used to compute source terms for hyperelastic model
module m_hyperelastic

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters
    ! ==========================================================================

    implicit none

    private; public :: s_calculate_cauchy_from_btensor

    contains

    !>  The following subroutine handles the calculation of the btensor.
        !!   The calculation of the btensor takes qprimvf.
        !! @param q_prim_vf Primitive variables
        !! @param btensor is the output
        !! calculate the grad_xi, grad_xi is a nxn tensor
        !! calculate the inverse of grad_xi to obtain F, F is a nxn tensor
        !! calculate the FFtranspose to obtain the btensor, btensor is nxn tensor
        !! btensor is symmetric, save the data space
    subroutine s_calculate_cauchy_from_btensor(btensor, q_prim_vf, ix, iy, iz)

        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_prim_vf
        type(scalar_field), dimension(b_size), intent(IN) :: btensor
        type(int_bounds_info), intent(IN) :: ix, iy, iz

        real(kind(0d0)), dimension(b_size-1) :: tensor
        real(kind(0d0)) :: trace
        integer :: i, j, k, l !< Generic loop iterators

        !$acc parallel loop collapse(3) gang vector default(present) private(trace, tensor)
        do l = iz%beg, iz%end
           do k = iy%beg, iy%end
              do j = ix%beg, ix%end
                    ! tensor is the symmetric tensor

                    !$acc loop seq
                    do i = 1, b_size - 1
                        tensor(i) = btensor(i)%sf(j, k, l) 
                    end do
                    ! calculate the trace of the tensor
                    trace = tensor(1)
                    if (num_dims == 2) then
                        trace = trace + tensor(3)
                    else
                        trace = trace + tensor(4) + tensor(6)
                    end if
                    ! calculate the deviatoric of the tensor
                    tensor(1) = tensor(1) - (1d0/3d0)*trace
                    if (num_dims == 2) then
                        tensor(3) = tensor(3) - (1d0/3d0)*trace
                    else
                        tensor(4) = tensor(4) - (1d0/3d0)*trace
                        tensor(6) = tensor(6) - (1d0/3d0)*trace
                    end if
                    ! dividing by the jacobian for neo-Hookean model
                    ! setting the tensor to the stresses for riemann solver

                    !$acc loop seq
                    do i = 1, b_size - 1
                        q_prim_vf(strxb+i)%sf(j, k, l) = tensor(i)/btensor(b_size)%sf(j, k, l)
                    end do

                end do
            end do
        end do
        !$acc end parallel loop

    end subroutine s_calculate_cauchy_from_btensor

end module m_hyperelastic
