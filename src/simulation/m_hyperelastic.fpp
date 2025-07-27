!>
!! @file m_hyperelastic.f90
!! @brief Contains module m_hyperelastic

#:include 'macros.fpp'

!> @brief This module consists of subroutines used in the calculation
!!              of the cauchy tensor

module m_hyperelastic

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_variables_conversion !< State variables type conversion procedures

    use m_finite_differences

    implicit none

    private; public :: s_hyperelastic_rmt_stress_update, &
 s_initialize_hyperelastic_module, &
 s_finalize_hyperelastic_module

    !! The btensor at the cell-interior Gaussian quadrature points.
    !! These tensor is needed to be calculated once and make the code DRY.
    type(vector_field) :: btensor !<
    $:GPU_DECLARE(create='[btensor]')

    real(wp), allocatable, dimension(:, :) :: fd_coeff_x
    real(wp), allocatable, dimension(:, :) :: fd_coeff_y
    real(wp), allocatable, dimension(:, :) :: fd_coeff_z
    $:GPU_DECLARE(create='[fd_coeff_x,fd_coeff_y, fd_coeff_z]')
    real(wp), allocatable, dimension(:) :: Gs
    $:GPU_DECLARE(create='[Gs]')

contains

    !>  The following subroutine handles the calculation of the btensor.
        !!   The calculation of the btensor takes qprimvf.
        !! @param q_prim_vf Primitive variables
        !! @param btensor is the output
        !! calculate the grad_xi, grad_xi is a nxn tensor
        !! calculate the inverse of grad_xi to obtain F, F is a nxn tensor
        !! calculate the FFtranspose to obtain the btensor, btensor is nxn tensor
        !! btensor is symmetric, save the data space
    impure subroutine s_initialize_hyperelastic_module
        integer :: i !< generic iterator

        @:ALLOCATE(btensor%vf(1:b_size))
        do i = 1, b_size
            @:ALLOCATE(btensor%vf(i)%sf(0:m, 0:n, 0:p))
        end do
        @:ACC_SETUP_VFs(btensor)

        @:ALLOCATE(Gs(1:num_fluids))
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            Gs(i) = fluid_pp(i)%G
        end do
        $:GPU_UPDATE(device='[Gs]')

        @:ALLOCATE(fd_coeff_x(-fd_number:fd_number, 0:m))
        if (n > 0) then
            @:ALLOCATE(fd_coeff_y(-fd_number:fd_number, 0:n))
        end if
        if (p > 0) then
            @:ALLOCATE(fd_coeff_z(-fd_number:fd_number, 0:p))
        end if

        ! Computing centered finite difference coefficients
        call s_compute_finite_difference_coefficients(m, x_cc, fd_coeff_x, buff_size, &
                                                      fd_number, fd_order)
        $:GPU_UPDATE(device='[fd_coeff_x]')
        if (n > 0) then
            call s_compute_finite_difference_coefficients(n, y_cc, fd_coeff_y, buff_size, &
                                                          fd_number, fd_order)
            $:GPU_UPDATE(device='[fd_coeff_y]')
        end if
        if (p > 0) then
            call s_compute_finite_difference_coefficients(p, z_cc, fd_coeff_z, buff_size, &
                                                          fd_number, fd_order)
            $:GPU_UPDATE(device='[fd_coeff_z]')
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
    subroutine s_hyperelastic_rmt_stress_update(q_cons_vf, q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf

        real(wp), dimension(tensor_size) :: tensora, tensorb
        real(wp), dimension(num_fluids) :: alpha_k, alpha_rho_k
        real(wp), dimension(2) :: Re
        real(wp) :: rho, gamma, pi_inf, qv
        real(wp) :: G_local
        integer :: j, k, l, i, r

        $:GPU_PARALLEL_LOOP(collapse=3, private='[alpha_K, alpha_rho_K, rho, &
            & gamma, pi_inf, qv, G_local, Re, tensora, tensorb]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        alpha_rho_k(i) = q_cons_vf(i)%sf(j, k, l)
                        alpha_k(i) = q_cons_vf(advxb + i - 1)%sf(j, k, l)
                    end do
                    ! If in simulation, use acc mixture subroutines
                    call s_convert_species_to_mixture_variables_acc(rho, gamma, pi_inf, qv, alpha_k, &
                                                                    alpha_rho_k, Re, G_local, Gs)
                    rho = max(rho, sgm_eps)
                    G_local = max(G_local, sgm_eps)
                    !if ( G_local <= verysmall ) G_K = 0._wp

                    if (G_local > verysmall) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, tensor_size
                            tensora(i) = 0._wp
                        end do
                        ! STEP 1: computing the grad_xi tensor using finite differences
                        ! grad_xi definition / organization
                        ! number for the tensor 1-3:  dxix_dx, dxiy_dx, dxiz_dx
                        ! 4-6 :                       dxix_dy, dxiy_dy, dxiz_dy
                        ! 7-9 :                       dxix_dz, dxiy_dz, dxiz_dz
                        $:GPU_LOOP(parallelism='[seq]')
                        do r = -fd_number, fd_number
                            ! derivatives in the x-direction
                            tensora(1) = tensora(1) + q_prim_vf(xibeg)%sf(j + r, k, l)*fd_coeff_x(r, j)
                            tensora(2) = tensora(2) + q_prim_vf(xibeg + 1)%sf(j + r, k, l)*fd_coeff_x(r, j)
                            tensora(3) = tensora(3) + q_prim_vf(xiend)%sf(j + r, k, l)*fd_coeff_x(r, j)
                            ! derivatives in the y-direction
                            tensora(4) = tensora(4) + q_prim_vf(xibeg)%sf(j, k + r, l)*fd_coeff_y(r, k)
                            tensora(5) = tensora(5) + q_prim_vf(xibeg + 1)%sf(j, k + r, l)*fd_coeff_y(r, k)
                            tensora(6) = tensora(6) + q_prim_vf(xiend)%sf(j, k + r, l)*fd_coeff_y(r, k)
                            ! derivatives in the z-direction
                            tensora(7) = tensora(7) + q_prim_vf(xibeg)%sf(j, k, l + r)*fd_coeff_z(r, l)
                            tensora(8) = tensora(8) + q_prim_vf(xibeg + 1)%sf(j, k, l + r)*fd_coeff_z(r, l)
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

                        if (tensorb(tensor_size) > verysmall) then
                            ! STEP 2c: computing the inverse of grad_xi tensor = F
                            ! tensorb is the adjoint, tensora becomes F
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, tensor_size - 1
                                tensora(i) = tensorb(i)/tensorb(tensor_size)
                            end do

                            ! STEP 2d: computing the J = det(F) = 1/det(\grad{\xi})
                            tensorb(tensor_size) = 1._wp/tensorb(tensor_size)

                            ! STEP 3: computing F transpose F
                            tensorb(1) = tensora(1)**2 + tensora(2)**2 + tensora(3)**2
                            tensorb(5) = tensora(4)**2 + tensora(5)**2 + tensora(6)**2
                            tensorb(9) = tensora(7)**2 + tensora(8)**2 + tensora(9)**2
                            tensorb(2) = tensora(1)*tensora(4) + tensora(2)*tensora(5) + tensora(3)*tensora(6)
                            tensorb(3) = tensora(1)*tensora(7) + tensora(2)*tensora(8) + tensora(3)*tensora(9)
                            tensorb(6) = tensora(4)*tensora(7) + tensora(5)*tensora(8) + tensora(6)*tensora(9)
                            ! STEP 4: update the btensor, this is consistent with Riemann solvers
                            #:for BIJ, TXY in [(1,1),(2,2),(3,5),(4,3),(5,6),(6,9)]
                                btensor%vf(${BIJ}$)%sf(j, k, l) = tensorb(${TXY}$)
                            #:endfor
                            ! store the determinant at the last entry of the btensor
                            btensor%vf(b_size)%sf(j, k, l) = tensorb(tensor_size)
                            ! STEP 5a: updating the Cauchy stress primitive scalar field
                            if (hyper_model == 1) then
                                call s_neoHookean_cauchy_solver(btensor%vf, q_prim_vf, G_local, j, k, l)
                            elseif (hyper_model == 2) then
                                call s_Mooney_Rivlin_cauchy_solver(btensor%vf, q_prim_vf, G_local, j, k, l)
                            end if
                            ! STEP 5b: updating the pressure field
                            q_prim_vf(E_idx)%sf(j, k, l) = q_prim_vf(E_idx)%sf(j, k, l) - &
                                                           G_local*q_prim_vf(xiend + 1)%sf(j, k, l)/gamma
                            ! STEP 5c: updating the Cauchy stress conservative scalar field
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, b_size - 1
                                q_cons_vf(strxb + i - 1)%sf(j, k, l) = &
                                    rho*q_prim_vf(strxb + i - 1)%sf(j, k, l)
                            end do
                        end if
                    end if
                end do
            end do
        end do
    end subroutine s_hyperelastic_rmt_stress_update

    !>  The following subroutine handles the calculation of the btensor.
        !!   The calculation of the btensor takes qprimvf.
        !! @param q_prim_vf Primitive variables
        !! @param btensor is the output
        !! calculate the grad_xi, grad_xi is a nxn tensor
        !! calculate the inverse of grad_xi to obtain F, F is a nxn tensor
        !! calculate the FFtranspose to obtain the btensor, btensor is nxn tensor
        !! btensor is symmetric, save the data space
    pure subroutine s_neoHookean_cauchy_solver(btensor_in, q_prim_vf, G_param, j, k, l)
        $:GPU_ROUTINE(parallelism='[seq]')
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        type(scalar_field), dimension(b_size), intent(inout) :: btensor_in
        real(wp), intent(in) :: G_param
        integer, intent(in) :: j, k, l

        real(wp) :: trace
        real(wp), parameter :: f13 = 1._wp/3._wp
        integer :: i !< Generic loop iterators

        ! tensor is the symmetric tensor & calculate the trace of the tensor
        trace = btensor_in(1)%sf(j, k, l) + btensor_in(3)%sf(j, k, l) + btensor_in(6)%sf(j, k, l)

        ! calculate the deviatoric of the tensor
        #:for IJ in [1,3,6]
            btensor_in(${IJ}$)%sf(j, k, l) = btensor_in(${IJ}$)%sf(j, k, l) - f13*trace
        #:endfor
        ! dividing by the jacobian for neo-Hookean model
        ! setting the tensor to the stresses for riemann solver
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, b_size - 1
            q_prim_vf(strxb + i - 1)%sf(j, k, l) = &
                G_param*btensor_in(i)%sf(j, k, l)/btensor_in(b_size)%sf(j, k, l)
        end do
        ! compute the invariant without the elastic modulus
        q_prim_vf(xiend + 1)%sf(j, k, l) = &
            0.5_wp*(trace - 3.0_wp)/btensor_in(b_size)%sf(j, k, l)

    end subroutine s_neoHookean_cauchy_solver

    !>  The following subroutine handles the calculation of the btensor.
        !!   The calculation of the btensor takes qprimvf.
        !! @param q_prim_vf Primitive variables
        !! @param btensor is the output
        !! calculate the grad_xi, grad_xi is a nxn tensor
        !! calculate the inverse of grad_xi to obtain F, F is a nxn tensor
        !! calculate the FFtranspose to obtain the btensor, btensor is nxn tensor
        !! btensor is symmetric, save the data space
    pure subroutine s_Mooney_Rivlin_cauchy_solver(btensor_in, q_prim_vf, G_param, j, k, l)
        $:GPU_ROUTINE(parallelism='[seq]')
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        type(scalar_field), dimension(b_size), intent(inout) :: btensor_in
        real(wp), intent(in) :: G_param
        integer, intent(in) :: j, k, l

        real(wp) :: trace
        real(wp), parameter :: f13 = 1._wp/3._wp
        integer :: i !< Generic loop iterators

        !TODO Make this 1D and 2D capable
        ! tensor is the symmetric tensor & calculate the trace of the tensor
        trace = btensor_in(1)%sf(j, k, l) + btensor_in(3)%sf(j, k, l) + btensor_in(6)%sf(j, k, l)

        ! calculate the deviatoric of the tensor
        btensor_in(1)%sf(j, k, l) = btensor_in(1)%sf(j, k, l) - f13*trace
        btensor_in(3)%sf(j, k, l) = btensor_in(3)%sf(j, k, l) - f13*trace
        btensor_in(6)%sf(j, k, l) = btensor_in(6)%sf(j, k, l) - f13*trace

        ! dividing by the jacobian for neo-Hookean model
        ! setting the tensor to the stresses for riemann solver
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, b_size - 1
            q_prim_vf(strxb + i - 1)%sf(j, k, l) = &
                G_param*btensor_in(i)%sf(j, k, l)/btensor_in(b_size)%sf(j, k, l)
        end do
        ! compute the invariant without the elastic modulus
        q_prim_vf(xiend + 1)%sf(j, k, l) = &
            0.5_wp*(trace - 3.0_wp)/btensor_in(b_size)%sf(j, k, l)

    end subroutine s_Mooney_Rivlin_cauchy_solver

    impure subroutine s_finalize_hyperelastic_module()

        integer :: i !< iterator

        ! Deallocating memory
        do i = 1, b_size
            @:DEALLOCATE(btensor%vf(i)%sf)
        end do
        @:DEALLOCATE(fd_coeff_x)
        if (n > 0) then
            @:DEALLOCATE(fd_coeff_y)
            if (p > 0) then
                @:DEALLOCATE(fd_coeff_z)
            end if
        end if

    end subroutine s_finalize_hyperelastic_module

end module m_hyperelastic
