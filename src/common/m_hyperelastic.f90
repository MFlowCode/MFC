!>
!! @file m_xi_tensor_calc.f90
!! @brief Contains module m_hyperelastic

!> @brief This module consists of subroutines used in the calculation 
!!              of the cauchy tensor

module m_hyperelastic

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters
    ! ==========================================================================

    implicit none

    private; public ::  s_compute_cauchy_solver, & 
 s_initialize_hyperelastic_module, &
 s_finalize_hyperelastic_module

    !> @name Abstract interface for creating function pointers
    !> @{
    abstract interface

        !> @name Abstract subroutine for the infinite relaxation solver
        !> @{
        subroutine s_abstract_hyperelastic_solver(btensor, q_prim_vf, G, j, k, l)

            import :: scalar_field, sys_size, b_size
            type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
            type(scalar_field), dimension(b_size), intent(in) :: btensor
            real(kind(0d0)), intent(in) :: G
            integer, intent(in) :: j, k, l
             
        end subroutine
        !> @}

    end interface
    !> @}

    procedure(s_abstract_hyperelastic_solver), pointer :: s_compute_cauchy_solver => null()

contains

     subroutine  s_initialize_hyperelastic_module()

        ! Associating procedural pointer to the subroutine that will be
        ! utilized to calculate the solution of a given Riemann problem
        !if (hyper_model == 1) then
            s_compute_cauchy_solver => s_neoHookean_cauchy_solver
        !elseif (riemann_solver == 2) then
        !    s_compute_cauchy_solver => s_Mooney_Rivlin_cauchy_solver
        !end if

     end subroutine

     !>  The following subroutine handles the calculation of the btensor.
        !!   The calculation of the btensor takes qprimvf.
        !! @param q_prim_vf Primitive variables
        !! @param btensor is the output
        !! calculate the grad_xi, grad_xi is a nxn tensor
        !! calculate the inverse of grad_xi to obtain F, F is a nxn tensor
        !! calculate the FFtranspose to obtain the btensor, btensor is nxn tensor
        !! btensor is symmetric, save the data space
     subroutine s_neoHookean_cauchy_solver(btensor, q_prim_vf, G, j, k, l)
#ifdef MFC_SIMULATION
        !$acc routine seq
#endif
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), dimension(b_size), intent(in) :: btensor
        integer, intent(in) :: j, k, l
        real(kind(0d0)), intent(in) :: G

        real(kind(0d0)), dimension(b_size - 1) :: tensor
        real(kind(0d0)) :: trace
        real(kind(0d0)) :: f13 = 1d0/3d0
        integer :: i !< Generic loop iterators

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
        btensor(1)%sf(j, k, l) = btensor(1)%sf(j, k, l) - f13*trace
        !if (num_dims == 2) then
        !    btensor(3)%sf(j,k,l) = btensor(3)%sf(j,k,l) - (1d0/3d0)*trace
        !else
          btensor(4)%sf(j, k, l) = btensor(4)%sf(j, k, l) - f13*trace
          btensor(6)%sf(j, k, l) = btensor(6)%sf(j, k, l) - f13*trace
        !end if
        ! dividing by the jacobian for neo-Hookean model
        ! setting the tensor to the stresses for riemann solver

        !$acc loop seq
        do i = 1, b_size - 1
           q_prim_vf(strxb + i)%sf(j, k, l) = G*btensor(i)%sf(j, k, l)/btensor(b_size)%sf(j, k, l)
        end do

        ! compute the invariant without the elastic modulus
        ! if (btensor(b_size)%sf(j,k,l) .gt. 0d0) then
        !    q_prim_vf(xiend + 1)%sf(j, k, l) = 0.5d0*(trace - 3.0d0)/btensor(b_size)%sf(j, k, l)
        ! else
        !     q_prim_vf(xiend+1)%sf(j,k,l) = 1d-12
        ! end if
     end subroutine s_neoHookean_cauchy_solver

     subroutine  s_finalize_hyperelastic_module()
        ! Disassociating procedural pointer to the subroutine which was
        ! utilized to calculate the solution of a given Riemann problem
        s_compute_cauchy_solver => null()
     end subroutine

end module m_hyperelastic
