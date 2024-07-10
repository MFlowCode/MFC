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
    subroutine s_calculate_cauchy_from_btensor(btensor, q_prim_vf, j, k, l)

        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_prim_vf
        type(scalar_field), dimension(b_size), intent(IN) :: btensor
        type(int_bounds_info), intent(IN) :: j, k, l

        real(kind(0d0)), dimension(b_size - 1) :: tensor
        real(kind(0d0)) :: trace
        integer :: i !< Generic loop iterators

        !!!$acc parallel loop collapse(3) gang vector default(present) private(trace)
        !do l = 0, p
        !    do k = 0, n
        !        do j = 0, m
                    ! tensor is the symmetric tensor & calculate the trace of the tensor
                    !trace = btensor(1)%sf(j,k,l)
                    !if (num_dims == 2) then
                    !    trace = trace + btensor(3)%sf(j,k,l)
                    !else
                    trace = btensor(1)%sf(j, k, l) + btensor(4)%sf(j, k, l) + btensor(6)%sf(j, k, l)
                    !end if
                    ! invariant calculation, saving it in the q_prim_vf field
                    !invariant1 = btensor(1)%sf(j, k, l)
                    !if (num_dims == 2) then
                    !    invariant1 = invariant1 + btensor(3)%sf(j, k, l)
                    !elseif (num_dims == 3) then
                    !    invariant1 = invariant1 + btensor(4)%sf(j, k, l) + btensor(6)%sf(j, k, l)
                    !end if

                    ! calculate the deviatoric of the tensor
                    btensor(1)%sf(j, k, l) = btensor(1)%sf(j, k, l) - (1d0/3d0)*trace
                    !if (num_dims == 2) then
                    !    btensor(3)%sf(j,k,l) = btensor(3)%sf(j,k,l) - (1d0/3d0)*trace
                    !else
                    btensor(4)%sf(j, k, l) = btensor(4)%sf(j, k, l) - (1d0/3d0)*trace
                    btensor(6)%sf(j, k, l) = btensor(6)%sf(j, k, l) - (1d0/3d0)*trace
                    !end if
                    ! dividing by the jacobian for neo-Hookean model
                    ! setting the tensor to the stresses for riemann solver

                    !$acc loop seq
                    do i = 1, b_size - 1
                        q_prim_vf(strxb + i)%sf(j, k, l) = btensor(i)%sf(j, k, l)/btensor(b_size)%sf(j, k, l)
                    end do

                    ! compute the invariant without the elastic modulus
                    ! if (btensor(b_size)%sf(j,k,l) .gt. 0d0) then
                    q_prim_vf(xiend + 1)%sf(j, k, l) = 0.5d0*(trace - 3.0d0)/btensor(b_size)%sf(j, k, l)
                    ! else
                    !     q_prim_vf(xiend+1)%sf(j,k,l) = 1d-12
                    ! end if
        !        end do
        !    end do
        !end do
        !!$acc end parallel loop

    end subroutine s_calculate_cauchy_from_btensor

end module m_hyperelastic
