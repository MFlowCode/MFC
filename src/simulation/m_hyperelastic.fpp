!>
!! @file m_hyperelastic.f90
!! @brief Contains module m_hyperelastic

#:include 'macros.fpp'

!> @brief This module consists of subroutines used in the calculation 
!!              of the cauchy tensor

module m_hyperelastic

    ! Dependencies =============================================================

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_variables_conversion !< State variables type conversion procedures

    use m_helper

    ! ==========================================================================

    implicit none

    private; public ::  s_hyperelastic_rmt_stress_update, &
 s_initialize_hyperelastic_module, &
 s_finalize_hyperelastic_module

    !> @name Abstract interface for creating function pointers
    !> @{
    abstract interface

        !> @name Abstract subroutine for the infinite relaxation solver
        !> @{
        subroutine s_abstract_hyperelastic_solver(btensor, q_prim_vf, G, j, k, l)

            import :: scalar_field, sys_size, b_size
            type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
            type(scalar_field), dimension(b_size), intent(inout) :: btensor
            real(kind(0d0)), intent(in) :: G
            integer, intent(in) :: j, k, l
             
        end subroutine
        !> @}

    end interface
    !> @}

    procedure(s_abstract_hyperelastic_solver), pointer :: s_compute_cauchy_solver => null()

    !! The btensor at the cell-interior Gaussian quadrature points.
    !! These tensor is needed to be calculated once and make the code DRY.
    type(vector_field) :: btensor !<
    !$acc declare create(btensor)

#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), allocatable, dimension(:, :), fd_coeff_x, fd_coeff_y, fd_coeff_z)
    !$acc declare link(fd_coeff_x,fd_coeff_y,fd_coeff_z)

#else

    real(kind(0d0)), allocatable, dimension(:, :) :: fd_coeff_x
    real(kind(0d0)), allocatable, dimension(:, :) :: fd_coeff_y
    real(kind(0d0)), allocatable, dimension(:, :) :: fd_coeff_z
    !$acc declare create(fd_coeff_x,fd_coeff_y,fd_coeff_z)
    real(kind(0d0)), allocatable, dimension(:) :: Gs
    !$acc declare create(Gs)
#endif

contains

     !>  The following subroutine handles the calculation of the btensor.
        !!   The calculation of the btensor takes qprimvf.
        !! @param q_prim_vf Primitive variables
        !! @param btensor is the output
        !! calculate the grad_xi, grad_xi is a nxn tensor
        !! calculate the inverse of grad_xi to obtain F, F is a nxn tensor
        !! calculate the FFtranspose to obtain the btensor, btensor is nxn tensor
        !! btensor is symmetric, save the data space
     subroutine s_initialize_hyperelastic_module()
        integer :: i !< generic iterator
   
        @:ALLOCATE(btensor%vf(1:b_size))
        do i = 1, b_size
          @:ALLOCATE(btensor%vf(i)%sf(0:m, 0:n, 0:p))
        end do
        @:ACC_SETUP_VFs(btensor)

        @:ALLOCATE(Gs(1:num_fluids))
        do i = 1, num_fluids
            Gs(i) = fluid_pp(i)%G
        end do
        !$acc update device(Gs)

        ! Associating procedural pointer to the subroutine that will be
        ! utilized to calculate the solution of a given Riemann problem
        !if (hyper_model == 1) then
            s_compute_cauchy_solver => s_neoHookean_cauchy_solver
        !elseif (riemann_solver == 2) then
        !    s_compute_cauchy_solver => s_Mooney_Rivlin_cauchy_solver
        !end if

        @:ALLOCATE_GLOBAL(fd_coeff_x(-fd_number:fd_number, 0:m))
        if (n > 0) then
           @:ALLOCATE_GLOBAL(fd_coeff_y(-fd_number:fd_number, 0:n))
        end if
        if (p > 0) then
           @:ALLOCATE_GLOBAL(fd_coeff_z(-fd_number:fd_number, 0:p))
        end if

        ! Computing centered finite difference coefficients
        call s_compute_finite_difference_coefficients(m, x_cc, fd_coeff_x, buff_size, &
                                                        fd_number, fd_order)
        !$acc update device(fd_coeff_x)
        if (n > 0) then
          call s_compute_finite_difference_coefficients(n, y_cc, fd_coeff_y, buff_size, &
                                                           fd_number, fd_order)
        !$acc update device(fd_coeff_y)
        end if
        if (p > 0) then
            call s_compute_finite_difference_coefficients(p, z_cc, fd_coeff_z, buff_size, &
                                                          fd_number, fd_order)
        !$acc update device(fd_coeff_z)
        end if

     end subroutine s_initialize_hyperelastic_module

     !>  The following subroutine handles the calculation of the btensor.
        !!   The calculation of the btensor takes qprimvf.
        !! @param q_prim_vf Primitive variables
        !! @param btensor is the output
        !! calculate the grad_xi, grad_xi is a nxn tensor
        !! calculate the inverse of grad_xi to obtain F, F is a nxn tensor
        !! calculate the FFtranspose to obtain the btensor, btensor is nxn tensor
        !! btensor is symmetric, save the data space
     subroutine s_hyperelastic_rmt_stress_update(q_cons_vf,q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(kind(0d0)), dimension(tensor_size) :: tensora, tensorb

        real(kind(0d0)), dimension(num_fluids) :: alpha_K, alpha_rho_K
        real(kind(0d0)), dimension(2) :: Re_K
        real(kind(0d0)) :: rho_K, gamma_K, pi_inf_K, qv_K
        real(kind(0d0)) :: G_K 
        logical :: flag
        integer :: j, k, l, i, r

        !$acc parallel loop collapse(3) gang vector default(present) private(alpha_K,alpha_rho_K,rho_K,gamma_K,pi_inf_K,qv_K,G_K,Re_K, tensora, tensorb, flag)
        do l = 0, p
           do k = 0, n
              do j = 0, m        
                flag = .true.
                !$acc loop seq
                do i = 1, num_fluids
                   alpha_rho_K(i) = q_cons_vf(i)%sf(j, k, l)
                   alpha_K(i) = q_cons_vf(advxb + i - 1)%sf(j, k, l)
                end do
                ! If in simulation, use acc mixture subroutines
                call s_convert_species_to_mixture_variables_acc(rho_K, gamma_K, pi_inf_K, qv_K, alpha_K, &
                              alpha_rho_K, Re_K, j, k, l, G_K, Gs)
                rho_K = max(rho_K, sgm_eps)
                if ( G_K .le. verysmall ) G_K = 0d0

                if ( G_K .gt. verysmall ) then
                  !$acc loop seq 
                  do i = 1, tensor_size
                    tensora(i) = 0d0
                  end do
                  ! STEP 1: computing the grad_xi tensor using finite differences
                  ! grad_xi definition / organization
                  ! number for the tensor 1-3:  dxix_dx, dxiy_dx, dxiz_dx
                  ! 4-6 :                       dxix_dy, dxiy_dy, dxiz_dy
                  ! 7-9 :                       dxix_dz, dxiy_dz, dxiz_dz
                  !$acc loop seq 
                  do r = -fd_number, fd_number        
                    ! derivatives in the x-direction
                    tensora(1) = tensora(1) + q_prim_vf(xibeg)%sf(j + r, k, l)*fd_coeff_x(r, j)
                    tensora(2) = tensora(2) + q_prim_vf(xibeg+1)%sf(j + r, k, l)*fd_coeff_x(r, j)
                    tensora(3) = tensora(3) + q_prim_vf(xiend)%sf(j + r, k, l)*fd_coeff_x(r, j)
                    ! derivatives in the y-direction
                    tensora(4) = tensora(4) + q_prim_vf(xibeg)%sf(j, k + r, l)*fd_coeff_y(r, k)
                    tensora(5) = tensora(5) + q_prim_vf(xibeg+1)%sf(j, k + r, l)*fd_coeff_y(r, k)
                    tensora(6) = tensora(6) + q_prim_vf(xiend)%sf(j, k + r, l)*fd_coeff_y(r, k)
                    ! derivatives in the z-direction
                    tensora(7) = tensora(7) + q_prim_vf(xibeg)%sf(j, k, l + r)*fd_coeff_z(r, l)
                    tensora(8) = tensora(8) + q_prim_vf(xibeg+1)%sf(j, k, l + r)*fd_coeff_z(r, l)
                    tensora(9) = tensora(9) + q_prim_vf(xiend)%sf(j, k, l + r)*fd_coeff_z(r, l)
                  end do 
                  ! STEP 2a: computing the adjoint of the grad_xi tensor for the inverse
                  tensorb(1) = tensora(5)*tensora(9) - tensora(6)*tensora(8)
                  tensorb(2) = -(tensora(2)*tensora(9) - tensora(3)*tensora(8))
                  tensorb(3) = tensora(2)*tensora(6) - tensora(3)*tensora(5)
                  tensorb(4) = -(tensora(4)*tensora(9) - tensora(6)*tensora(7))
                  tensorb(5) = tensora(1)*tensora(9) - tensora(3)*tensora(7)
                  tensorb(6) = -(tensora(1)*tensora(6) - tensora(4)*tensora(3))
                  tensorb(7) = tensora(4)*tensora(8) - tensora(5)*tensora(7)
                  tensorb(8) = -(tensora(1)*tensora(8) - tensora(2)*tensora(7))
                  tensorb(9) = tensora(1)*tensora(5) - tensora(2)*tensora(4)

                  ! STEP 2b: computing the determinant of the grad_xi tensor
                  tensorb(tensor_size) = tensora(1)*(tensora(5)*tensora(9) - tensora(6)*tensora(8)) &
                                    - tensora(2)*(tensora(4)*tensora(9) - tensora(6)*tensora(7)) &
                                    + tensora(3)*(tensora(4)*tensora(8) - tensora(5)*tensora(7))

                  ! STEP 2c: computing the inverse of grad_xi tensor = F
                  ! tensorb is the adjoint, tensora becomes F
                  !$acc loop seq
                  do i = 1, tensor_size - 1
                    tensora(i) = tensorb(i)/tensorb(tensor_size)
                  end do
 
                  ! STEP 2d: computing the J = det(F) = 1/det(\grad{\xi})
                  tensorb(tensor_size) = 1d0/tensorb(tensor_size)
 
                  ! STEP 3: computing F tranpose F
                  tensorb(1) = tensora(1)**2 + tensora(2)**2 + tensora(3)**2
                  tensorb(5) = tensora(4)**2 + tensora(5)**2 + tensora(6)**2
                  tensorb(9) = tensora(7)**2 + tensora(8)**2 + tensora(9)**2
                  tensorb(2) = tensora(1)*tensora(4) + tensora(2)*tensora(5) + tensora(3)*tensora(6)
                  tensorb(3) = tensora(1)*tensora(7) + tensora(2)*tensora(8) + tensora(3)*tensora(9)
                  tensorb(6) = tensora(4)*tensora(7) + tensora(5)*tensora(8) + tensora(6)*tensora(9)
                else
                   flag = .false.
                end if

                if (flag) then
                  ! STEP 4: update the btensor, this is consistent with Riemann solvers
                  btensor%vf(1)%sf(j, k, l) = tensorb(1)
                  btensor%vf(2)%sf(j, k, l) = tensorb(2)
                  btensor%vf(3)%sf(j, k, l) = tensorb(5)
                  btensor%vf(4)%sf(j, k, l) = tensorb(3)
                  btensor%vf(5)%sf(j, k, l) = tensorb(6)
                  btensor%vf(6)%sf(j, k, l) = tensorb(9)
                  ! store the determinant at the last entry of the btensor 
                  btensor%vf(b_size)%sf(j, k, l) = tensorb(tensor_size)
   
                  ! STEP 5a: updating the Cauchy stress primitive scalar field
                  call s_compute_cauchy_solver(btensor%vf, q_prim_vf, G_K, j, k, l)
                  ! STEP 5b: updating the pressure field
                  q_prim_vf(E_idx)%sf(j, k, l) = q_prim_vf(E_idx)%sf(j, k, l) - &
                        G_K*q_prim_vf(xiend + 1)%sf(j, k, l)/gamma_K
                  ! STEP 5c: updating the Cauchy stress conservative scalar field
                  !$acc loop seq
                  do i = 1, b_size - 1
                    q_cons_vf(strxb + i - 1)%sf(j, k, l) =  & 
                      rho_K*q_prim_vf(strxb + i - 1)%sf(j, k, l)
                  end do
                end if
            end do
          end do
        end do
        !$acc end parallel loop
     end subroutine s_hyperelastic_rmt_stress_update

     !>  The following subroutine handles the calculation of the btensor.
        !!   The calculation of the btensor takes qprimvf.
        !! @param q_prim_vf Primitive variables
        !! @param btensor is the output
        !! calculate the grad_xi, grad_xi is a nxn tensor
        !! calculate the inverse of grad_xi to obtain F, F is a nxn tensor
        !! calculate the FFtranspose to obtain the btensor, btensor is nxn tensor
        !! btensor is symmetric, save the data space
     subroutine s_neoHookean_cauchy_solver(btensor, q_prim_vf, G, j, k, l)
        !$acc routine seq
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        type(scalar_field), dimension(b_size), intent(inout) :: btensor
        real(kind(0d0)), intent(in) :: G
        integer, intent(in) :: j, k, l

        real(kind(0d0)) :: trace
        real(kind(0d0)) :: f13 = 1d0/3d0
        integer :: i !< Generic loop iterators

        !TODO Make this 1D and 2D capable
        ! tensor is the symmetric tensor & calculate the trace of the tensor
        trace = btensor(1)%sf(j, k, l) + btensor(4)%sf(j, k, l) + btensor(6)%sf(j, k, l)

        ! calculate the deviatoric of the tensor
        btensor(1)%sf(j, k, l) = btensor(1)%sf(j, k, l) - f13*trace
        btensor(4)%sf(j, k, l) = btensor(4)%sf(j, k, l) - f13*trace
        btensor(6)%sf(j, k, l) = btensor(6)%sf(j, k, l) - f13*trace

        ! dividing by the jacobian for neo-Hookean model
        ! setting the tensor to the stresses for riemann solver
        !$acc loop seq
        do i = 1, b_size - 1
          q_prim_vf(strxb + i - 1)%sf(j, k, l) =  & 
            G*btensor(i)%sf(j, k, l)/btensor(b_size)%sf(j, k, l)
        end do
        ! compute the invariant without the elastic modulus
        q_prim_vf(xiend + 1)%sf(j, k, l) =  & 
           0.5d0*(trace - 3.0d0)/btensor(b_size)%sf(j, k, l)

     end subroutine s_neoHookean_cauchy_solver

    subroutine  s_finalize_hyperelastic_module()

      integer :: i !< iterator

      ! Disassociating procedural pointer to the subroutine which was
      ! utilized to calculate the solution of a given Riemann problem
      s_compute_cauchy_solver => null()

      ! Deallocating memory
      do i = 1, b_size
           @:DEALLOCATE_GLOBAL(btensor%vf(i)%sf)
      end do
      @:DEALLOCATE_GLOBAL(fd_coeff_x)
      if (n > 0) then
         @:DEALLOCATE_GLOBAL(fd_coeff_y)
         if (p > 0) then
            @:DEALLOCATE_GLOBAL(fd_coeff_z)
         end if
      end if

    end subroutine s_finalize_hyperelastic_module

end module m_hyperelastic
