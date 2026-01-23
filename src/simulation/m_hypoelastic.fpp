!>
!! @file m_hypoelastic.f90
!! @brief Contains module m_hypoelastic

#:include 'macros.fpp'

!> @brief This module is used to compute source terms for hypoelastic model
module m_hypoelastic

    use m_derived_types        !< Definitions of the derived types
    use m_global_parameters    !< Definitions of the global parameters
    use m_finite_differences
    use m_helper

    implicit none

    private; public :: s_initialize_hypoelastic_module, &
 s_finalize_hypoelastic_module, &
 s_compute_hypoelastic_rhs, &
 s_compute_damage_state

    real(wp), allocatable, dimension(:) :: Gs_hypo
    $:GPU_DECLARE(create='[Gs_hypo]')

    real(wp), allocatable, dimension(:, :, :) :: du_dx_hypo, du_dy_hypo, du_dz_hypo
    real(wp), allocatable, dimension(:, :, :) :: dv_dx_hypo, dv_dy_hypo, dv_dz_hypo
    real(wp), allocatable, dimension(:, :, :) :: dw_dx_hypo, dw_dy_hypo, dw_dz_hypo
    $:GPU_DECLARE(create='[du_dx_hypo,du_dy_hypo,du_dz_hypo,dv_dx_hypo,dv_dy_hypo,dv_dz_hypo,dw_dx_hypo,dw_dy_hypo,dw_dz_hypo]')

    real(wp), allocatable, dimension(:, :, :) :: rho_K_field, G_K_field
    $:GPU_DECLARE(create='[rho_K_field,G_K_field]')

    real(wp), allocatable, dimension(:, :) :: fd_coeff_x_hypo
    real(wp), allocatable, dimension(:, :) :: fd_coeff_y_hypo
    real(wp), allocatable, dimension(:, :) :: fd_coeff_z_hypo
    $:GPU_DECLARE(create='[fd_coeff_x_hypo,fd_coeff_y_hypo,fd_coeff_z_hypo]')

contains

    impure subroutine s_initialize_hypoelastic_module

        integer :: i

        @:ALLOCATE(Gs_hypo(1:num_fluids))
        @:ALLOCATE(rho_K_field(0:m,0:n,0:p), G_K_field(0:m,0:n,0:p))
        @:ALLOCATE(du_dx_hypo(0:m,0:n,0:p))
        if (n > 0) then
            @:ALLOCATE(du_dy_hypo(0:m,0:n,0:p), dv_dx_hypo(0:m,0:n,0:p), dv_dy_hypo(0:m,0:n,0:p))
            if (p > 0) then
                @:ALLOCATE(du_dz_hypo(0:m,0:n,0:p), dv_dz_hypo(0:m,0:n,0:p))
                @:ALLOCATE(dw_dx_hypo(0:m,0:n,0:p), dw_dy_hypo(0:m,0:n,0:p), dw_dz_hypo(0:m,0:n,0:p))
            end if
        end if

        do i = 1, num_fluids
            Gs_hypo(i) = fluid_pp(i)%G
        end do
        $:GPU_UPDATE(device='[Gs_hypo]')

        @:ALLOCATE(fd_coeff_x_hypo(-fd_number:fd_number, 0:m))
        if (n > 0) then
            @:ALLOCATE(fd_coeff_y_hypo(-fd_number:fd_number, 0:n))
        end if
        if (p > 0) then
            @:ALLOCATE(fd_coeff_z_hypo(-fd_number:fd_number, 0:p))
        end if

        ! Computing centered finite difference coefficients
        call s_compute_finite_difference_coefficients(m, x_cc, fd_coeff_x_hypo, buff_size, &
                                                      fd_number, fd_order)
        $:GPU_UPDATE(device='[fd_coeff_x_hypo]')
        if (n > 0) then
            call s_compute_finite_difference_coefficients(n, y_cc, fd_coeff_y_hypo, buff_size, &
                                                          fd_number, fd_order)
            $:GPU_UPDATE(device='[fd_coeff_y_hypo]')
        end if
        if (p > 0) then
            call s_compute_finite_difference_coefficients(p, z_cc, fd_coeff_z_hypo, buff_size, &
                                                          fd_number, fd_order)
            $:GPU_UPDATE(device='[fd_coeff_z_hypo]')
        end if

    end subroutine s_initialize_hypoelastic_module

    !>  The purpose of this procedure is to compute the source terms
        !!      that are needed for the elastic stress equations
        !!  @param idir Dimension splitting index
        !!  @param q_prim_vf Primitive variables
        !!  @param rhs_vf rhs variables
    subroutine s_compute_hypoelastic_rhs(idir, q_prim_vf, rhs_vf)

        integer, intent(in) :: idir
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf

        real(wp) :: rho_K, G_K

        integer :: i, k, l, q, r !< Loop variables
        integer :: ndirs  !< Number of coordinate directions

        ndirs = 1; if (n > 0) ndirs = 2; if (p > 0) ndirs = 3

        if (idir == 1) then
            ! calculate velocity gradients + rho_K and G_K
            ! TODO: re-organize these loops one by one for GPU efficiency if possible?

            $:GPU_PARALLEL_LOOP(collapse=3)
            do q = 0, p
                do l = 0, n
                    do k = 0, m
                        du_dx_hypo(k, l, q) = 0._wp
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            $:GPU_PARALLEL_LOOP(collapse=3)
            do q = 0, p
                do l = 0, n
                    do k = 0, m
                        $:GPU_LOOP(parallelism='[seq]')
                        do r = -fd_number, fd_number
                            du_dx_hypo(k, l, q) = du_dx_hypo(k, l, q) &
                                                  + q_prim_vf(momxb)%sf(k + r, l, q)*fd_coeff_x_hypo(r, k)
                        end do

                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            if (ndirs > 1) then
                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            du_dy_hypo(k, l, q) = 0._wp; dv_dx_hypo(k, l, q) = 0._wp; dv_dy_hypo(k, l, q) = 0._wp
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                $:GPU_PARALLEL_LOOP(collapse=3)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            $:GPU_LOOP(parallelism='[seq]')
                            do r = -fd_number, fd_number
                                du_dy_hypo(k, l, q) = du_dy_hypo(k, l, q) &
                                                      + q_prim_vf(momxb)%sf(k, l + r, q)*fd_coeff_y_hypo(r, l)
                                dv_dx_hypo(k, l, q) = dv_dx_hypo(k, l, q) &
                                                      + q_prim_vf(momxb + 1)%sf(k + r, l, q)*fd_coeff_x_hypo(r, k)
                                dv_dy_hypo(k, l, q) = dv_dy_hypo(k, l, q) &
                                                      + q_prim_vf(momxb + 1)%sf(k, l + r, q)*fd_coeff_y_hypo(r, l)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                ! 3D
                if (ndirs == 3) then

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                du_dz_hypo(k, l, q) = 0._wp; dv_dz_hypo(k, l, q) = 0._wp; dw_dx_hypo(k, l, q) = 0._wp; 
                                dw_dy_hypo(k, l, q) = 0._wp; dw_dz_hypo(k, l, q) = 0._wp; 
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                $:GPU_LOOP(parallelism='[seq]')
                                do r = -fd_number, fd_number
                                    du_dz_hypo(k, l, q) = du_dz_hypo(k, l, q) &
                                                          + q_prim_vf(momxb)%sf(k, l, q + r)*fd_coeff_z_hypo(r, q)
                                    dv_dz_hypo(k, l, q) = dv_dz_hypo(k, l, q) &
                                                          + q_prim_vf(momxb + 1)%sf(k, l, q + r)*fd_coeff_z_hypo(r, q)
                                    dw_dx_hypo(k, l, q) = dw_dx_hypo(k, l, q) &
                                                          + q_prim_vf(momxe)%sf(k + r, l, q)*fd_coeff_x_hypo(r, k)
                                    dw_dy_hypo(k, l, q) = dw_dy_hypo(k, l, q) &
                                                          + q_prim_vf(momxe)%sf(k, l + r, q)*fd_coeff_y_hypo(r, l)
                                    dw_dz_hypo(k, l, q) = dw_dz_hypo(k, l, q) &
                                                          + q_prim_vf(momxe)%sf(k, l, q + r)*fd_coeff_z_hypo(r, q)
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
            end if

            $:GPU_PARALLEL_LOOP(collapse=3,private='[rho_K, G_K]')
            do q = 0, p
                do l = 0, n
                    do k = 0, m
                        rho_K = 0._wp; G_K = 0._wp
                        do i = 1, num_fluids
                            rho_K = rho_K + q_prim_vf(i)%sf(k, l, q) !alpha_rho_K(1)
                            G_K = G_K + q_prim_vf(advxb - 1 + i)%sf(k, l, q)*Gs_hypo(i)  !alpha_K(1) * Gs_hypo(1)
                        end do

                        if (cont_damage) G_K = G_K*max((1._wp - q_prim_vf(damage_idx)%sf(k, l, q)), 0._wp)

                        rho_K_field(k, l, q) = rho_K
                        G_K_field(k, l, q) = G_K

                        !TODO: take this out if not needed
                        if (G_K < verysmall) then
                            G_K_field(k, l, q) = 0
                        end if
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            ! apply rhs source term to elastic stress equation
            $:GPU_PARALLEL_LOOP(collapse=3)
            do q = 0, p
                do l = 0, n
                    do k = 0, m
                        rhs_vf(strxb)%sf(k, l, q) = &
                            rhs_vf(strxb)%sf(k, l, q) + rho_K_field(k, l, q)* &
                            ((4._wp*G_K_field(k, l, q)/3._wp) + &
                             q_prim_vf(strxb)%sf(k, l, q))* &
                            du_dx_hypo(k, l, q)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

        elseif (idir == 2) then
            $:GPU_PARALLEL_LOOP(collapse=3)
            do q = 0, p
                do l = 0, n
                    do k = 0, m
                        rhs_vf(strxb)%sf(k, l, q) = rhs_vf(strxb)%sf(k, l, q) + rho_K_field(k, l, q)* &
                                                    (q_prim_vf(strxb + 1)%sf(k, l, q)*du_dy_hypo(k, l, q) + &
                                                     q_prim_vf(strxb + 1)%sf(k, l, q)*du_dy_hypo(k, l, q) - &
                                                     q_prim_vf(strxb)%sf(k, l, q)*dv_dy_hypo(k, l, q) - &
                                                     2._wp*G_K_field(k, l, q)*(1._wp/3._wp)*dv_dy_hypo(k, l, q))

                        rhs_vf(strxb + 1)%sf(k, l, q) = rhs_vf(strxb + 1)%sf(k, l, q) + rho_K_field(k, l, q)* &
                                                        (q_prim_vf(strxb + 1)%sf(k, l, q)*du_dx_hypo(k, l, q) + &
                                                         q_prim_vf(strxb)%sf(k, l, q)*dv_dx_hypo(k, l, q) - &
                                                         q_prim_vf(strxb + 1)%sf(k, l, q)*du_dx_hypo(k, l, q) + &
                                                         q_prim_vf(strxb + 2)%sf(k, l, q)*du_dy_hypo(k, l, q) + &
                                                         q_prim_vf(strxb + 1)%sf(k, l, q)*dv_dy_hypo(k, l, q) - &
                                                         q_prim_vf(strxb + 1)%sf(k, l, q)*dv_dy_hypo(k, l, q) + &
                                                         2._wp*G_K_field(k, l, q)*(1._wp/2._wp)*(du_dy_hypo(k, l, q) + &
                                                                                                 dv_dx_hypo(k, l, q)))

                        rhs_vf(strxb + 2)%sf(k, l, q) = rhs_vf(strxb + 2)%sf(k, l, q) + rho_K_field(k, l, q)* &
                                                        (q_prim_vf(strxb + 1)%sf(k, l, q)*dv_dx_hypo(k, l, q) + &
                                                         q_prim_vf(strxb + 1)%sf(k, l, q)*dv_dx_hypo(k, l, q) - &
                                                         q_prim_vf(strxb + 2)%sf(k, l, q)*du_dx_hypo(k, l, q) + &
                                                         q_prim_vf(strxb + 2)%sf(k, l, q)*dv_dy_hypo(k, l, q) + &
                                                         q_prim_vf(strxb + 2)%sf(k, l, q)*dv_dy_hypo(k, l, q) - &
                                                         q_prim_vf(strxb + 2)%sf(k, l, q)*dv_dy_hypo(k, l, q) + &
                                                         2._wp*G_K_field(k, l, q)*(dv_dy_hypo(k, l, q) - (1._wp/3._wp)* &
                                                                                   (du_dx_hypo(k, l, q) + &
                                                                                    dv_dy_hypo(k, l, q))))
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

        elseif (idir == 3) then
            $:GPU_PARALLEL_LOOP(collapse=3)
            do q = 0, p
                do l = 0, n
                    do k = 0, m
                        rhs_vf(strxb)%sf(k, l, q) = rhs_vf(strxb)%sf(k, l, q) + rho_K_field(k, l, q)* &
                                                    (q_prim_vf(strxb + 3)%sf(k, l, q)*du_dz_hypo(k, l, q) + &
                                                     q_prim_vf(strxb + 3)%sf(k, l, q)*du_dz_hypo(k, l, q) - &
                                                     q_prim_vf(strxb)%sf(k, l, q)*dw_dz_hypo(k, l, q) - &
                                                     2._wp*G_K_field(k, l, q)*(1._wp/3._wp)*dw_dz_hypo(k, l, q))

                        rhs_vf(strxb + 1)%sf(k, l, q) = rhs_vf(strxb + 1)%sf(k, l, q) + rho_K_field(k, l, q)* &
                                                        (q_prim_vf(strxb + 4)%sf(k, l, q)*du_dz_hypo(k, l, q) + &
                                                         q_prim_vf(strxb + 3)%sf(k, l, q)*dv_dz_hypo(k, l, q) - &
                                                         q_prim_vf(strxb + 1)%sf(k, l, q)*dw_dz_hypo(k, l, q))

                        rhs_vf(strxb + 2)%sf(k, l, q) = rhs_vf(strxb + 2)%sf(k, l, q) + rho_K_field(k, l, q)* &
                                                        (q_prim_vf(strxb + 4)%sf(k, l, q)*dv_dz_hypo(k, l, q) + &
                                                         q_prim_vf(strxb + 4)%sf(k, l, q)*dv_dz_hypo(k, l, q) - &
                                                         q_prim_vf(strxb + 2)%sf(k, l, q)*dw_dz_hypo(k, l, q) - &
                                                         2._wp*G_K_field(k, l, q)*(1._wp/3._wp)*dw_dz_hypo(k, l, q))

                        rhs_vf(strxb + 3)%sf(k, l, q) = rhs_vf(strxb + 3)%sf(k, l, q) + rho_K_field(k, l, q)* &
                                                        (q_prim_vf(strxb + 3)%sf(k, l, q)*du_dx_hypo(k, l, q) + &
                                                         q_prim_vf(strxb)%sf(k, l, q)*dw_dx_hypo(k, l, q) - &
                                                         q_prim_vf(strxb + 3)%sf(k, l, q)*du_dx_hypo(k, l, q) + &
                                                         q_prim_vf(strxb + 4)%sf(k, l, q)*du_dy_hypo(k, l, q) + &
                                                         q_prim_vf(strxb + 1)%sf(k, l, q)*dw_dy_hypo(k, l, q) - &
                                                         q_prim_vf(strxb + 3)%sf(k, l, q)*dv_dy_hypo(k, l, q) + &
                                                         q_prim_vf(strxb + 5)%sf(k, l, q)*du_dz_hypo(k, l, q) + &
                                                         q_prim_vf(strxb + 3)%sf(k, l, q)*dw_dz_hypo(k, l, q) - &
                                                         q_prim_vf(strxb + 3)%sf(k, l, q)*dw_dz_hypo(k, l, q) + &
                                                         2._wp*G_K_field(k, l, q)*(1._wp/2._wp)*(du_dz_hypo(k, l, q) + &
                                                                                                 dw_dx_hypo(k, l, q)))

                        rhs_vf(strxb + 4)%sf(k, l, q) = rhs_vf(strxb + 4)%sf(k, l, q) + rho_K_field(k, l, q)* &
                                                        (q_prim_vf(strxb + 3)%sf(k, l, q)*dv_dx_hypo(k, l, q) + &
                                                         q_prim_vf(strxb + 1)%sf(k, l, q)*dw_dx_hypo(k, l, q) - &
                                                         q_prim_vf(strxb + 4)%sf(k, l, q)*du_dx_hypo(k, l, q) + &
                                                         q_prim_vf(strxb + 4)%sf(k, l, q)*dv_dy_hypo(k, l, q) + &
                                                         q_prim_vf(strxb + 2)%sf(k, l, q)*dw_dy_hypo(k, l, q) - &
                                                         q_prim_vf(strxb + 4)%sf(k, l, q)*dv_dy_hypo(k, l, q) + &
                                                         q_prim_vf(strxb + 5)%sf(k, l, q)*dv_dz_hypo(k, l, q) + &
                                                         q_prim_vf(strxb + 4)%sf(k, l, q)*dw_dz_hypo(k, l, q) - &
                                                         q_prim_vf(strxb + 4)%sf(k, l, q)*dw_dz_hypo(k, l, q) + &
                                                         2._wp*G_K_field(k, l, q)*(1._wp/2._wp)*(dv_dz_hypo(k, l, q) + &
                                                                                                 dw_dy_hypo(k, l, q)))

                        rhs_vf(strxe)%sf(k, l, q) = rhs_vf(strxe)%sf(k, l, q) + rho_K_field(k, l, q)* &
                                                    (q_prim_vf(strxe - 2)%sf(k, l, q)*dw_dx_hypo(k, l, q) + &
                                                     q_prim_vf(strxe - 2)%sf(k, l, q)*dw_dx_hypo(k, l, q) - &
                                                     q_prim_vf(strxe)%sf(k, l, q)*du_dx_hypo(k, l, q) + &
                                                     q_prim_vf(strxe - 1)%sf(k, l, q)*dw_dy_hypo(k, l, q) + &
                                                     q_prim_vf(strxe - 1)%sf(k, l, q)*dw_dy_hypo(k, l, q) - &
                                                     q_prim_vf(strxe)%sf(k, l, q)*dv_dy_hypo(k, l, q) + &
                                                     q_prim_vf(strxe)%sf(k, l, q)*dw_dz_hypo(k, l, q) + &
                                                     q_prim_vf(strxe)%sf(k, l, q)*dw_dz_hypo(k, l, q) - &
                                                     q_prim_vf(strxe)%sf(k, l, q)*dw_dz_hypo(k, l, q) + &
                                                     2._wp*G_K_field(k, l, q)*(dw_dz_hypo(k, l, q) - (1._wp/3._wp)* &
                                                                               (du_dx_hypo(k, l, q) + &
                                                                                dv_dy_hypo(k, l, q) + &
                                                                                dw_dz_hypo(k, l, q))))
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        if (cyl_coord .and. idir == 2) then

            $:GPU_PARALLEL_LOOP(collapse=3)
            do q = 0, p
                do l = 0, n
                    do k = 0, m
                        ! S_xx -= rho * v/r * (tau_xx + 2/3*G)
                        rhs_vf(strxb)%sf(k, l, q) = rhs_vf(strxb)%sf(k, l, q) - &
                                                    rho_K_field(k, l, q)*q_prim_vf(momxb + 1)%sf(k, l, q)/y_cc(l)* &
                                                    (q_prim_vf(strxb)%sf(k, l, q) + (2._wp/3._wp)*G_K_field(k, l, q)) ! tau_xx + 2/3*G

                        ! S_xr -= rho * v/r * tau_xr
                        rhs_vf(strxb + 1)%sf(k, l, q) = rhs_vf(strxb + 1)%sf(k, l, q) - &
                                                        rho_K_field(k, l, q)*q_prim_vf(momxb + 1)%sf(k, l, q)/y_cc(l)* &
                                                        q_prim_vf(strxb + 1)%sf(k, l, q) ! tau_xx

                        ! S_rr -= rho * v/r * (tau_rr + 2/3*G)
                        rhs_vf(strxb + 2)%sf(k, l, q) = rhs_vf(strxb + 2)%sf(k, l, q) - &
                                                        rho_K_field(k, l, q)*q_prim_vf(momxb + 1)%sf(k, l, q)/y_cc(l)* &
                                                        (q_prim_vf(strxb + 2)%sf(k, l, q) + (2._wp/3._wp)*G_K_field(k, l, q)) ! tau_rr + 2/3*G

                        ! S_thetatheta += rho * ( -(tau_thetatheta + 2/3*G)*(du/dx + dv/dr + v/r) + 2*(tau_thetatheta + G)*v/r )
                        rhs_vf(strxb + 3)%sf(k, l, q) = rhs_vf(strxb + 3)%sf(k, l, q) + &
                                                        rho_K_field(k, l, q)*( &
                                                        -(q_prim_vf(strxb + 3)%sf(k, l, q) + (2._wp/3._wp)*G_K_field(k, l, q))* &
                                                        (du_dx_hypo(k, l, q) + dv_dy_hypo(k, l, q) + q_prim_vf(momxb + 1)%sf(k, l, q)/y_cc(l)) &
                                                        + 2._wp*(q_prim_vf(strxb + 3)%sf(k, l, q) + G_K_field(k, l, q))*q_prim_vf(momxb + 1)%sf(k, l, q)/y_cc(l))
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

        end if

    end subroutine s_compute_hypoelastic_rhs

    impure subroutine s_finalize_hypoelastic_module()

        @:DEALLOCATE(Gs_hypo)
        @:DEALLOCATE(rho_K_field, G_K_field)
        @:DEALLOCATE(du_dx_hypo)
        @:DEALLOCATE(fd_coeff_x_hypo)
        if (n > 0) then
            @:DEALLOCATE(du_dy_hypo,dv_dx_hypo,dv_dy_hypo)
            @:DEALLOCATE(fd_coeff_y_hypo)
            if (p > 0) then
                @:DEALLOCATE(du_dz_hypo, dv_dz_hypo, dw_dx_hypo, dw_dy_hypo, dw_dz_hypo)
                @:DEALLOCATE(fd_coeff_z_hypo)
            end if
        end if

    end subroutine s_finalize_hypoelastic_module

    subroutine s_compute_damage_state(q_cons_vf, rhs_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf

        real(wp) :: tau_p ! principal stress
        real(wp) :: tau_xx, tau_xy, tau_yy, tau_zz, tau_yz, tau_xz
        real(wp) :: I1, I2, I3, argument, phi, sqrt_term_1, sqrt_term_2, temp
        integer :: q, l, k

        if (n == 0) then
            l = 0; q = 0
            $:GPU_PARALLEL_LOOP()
            do k = 0, m
                rhs_vf(damage_idx)%sf(k, l, q) = (alpha_bar*max(abs(real(q_cons_vf(stress_idx%beg)%sf(k, l, q), kind=wp)) - tau_star, 0._wp))**cont_damage_s
            end do
            $:END_GPU_PARALLEL_LOOP()
        elseif (p == 0) then
            q = 0
            $:GPU_PARALLEL_LOOP(collapse=2, private='[tau_p]')
            do l = 0, n
                do k = 0, m
                    ! Maximum principal stress
                    tau_p = 0.5_wp*(q_cons_vf(stress_idx%beg)%sf(k, l, q) + &
                                    q_cons_vf(stress_idx%beg + 2)%sf(k, l, q)) + &
                            sqrt((q_cons_vf(stress_idx%beg)%sf(k, l, q) - &
                                  q_cons_vf(stress_idx%beg + 2)%sf(k, l, q))**2.0_wp + &
                                 4._wp*q_cons_vf(stress_idx%beg + 1)%sf(k, l, q)**2.0_wp)/2._wp

                    rhs_vf(damage_idx)%sf(k, l, q) = (alpha_bar*max(tau_p - tau_star, 0._wp))**cont_damage_s
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        else
            $:GPU_PARALLEL_LOOP(collapse=3, private='[tau_xx, tau_xy, tau_yy, tau_xz, tau_yz, tau_zz, I1, I2, I3, temp, sqrt_term_1, sqrt_term_2, argument, phi, tau_p]')
            do q = 0, p
                do l = 0, n
                    do k = 0, m
                        tau_xx = q_cons_vf(stress_idx%beg)%sf(k, l, q)
                        tau_xy = q_cons_vf(stress_idx%beg + 1)%sf(k, l, q)
                        tau_yy = q_cons_vf(stress_idx%beg + 2)%sf(k, l, q)
                        tau_xz = q_cons_vf(stress_idx%beg + 3)%sf(k, l, q)
                        tau_yz = q_cons_vf(stress_idx%beg + 4)%sf(k, l, q)
                        tau_zz = q_cons_vf(stress_idx%beg + 5)%sf(k, l, q)

                        ! Invariants of the stress tensor
                        I1 = tau_xx + tau_yy + tau_zz
                        I2 = tau_xx*tau_yy + tau_xx*tau_zz + tau_yy*tau_zz - &
                             (tau_xy**2.0_wp + tau_xz**2.0_wp + tau_yz**2.0_wp)
                        I3 = tau_xx*tau_yy*tau_zz + 2.0_wp*tau_xy*tau_xz*tau_yz - &
                             tau_xx*tau_yz**2.0_wp - tau_yy*tau_xz**2.0_wp - tau_zz*tau_xy**2.0_wp

                        ! Maximum principal stress
                        temp = I1**2.0_wp - 3.0_wp*I2
                        sqrt_term_1 = sqrt(max(temp, 0.0_wp))
                        if (sqrt_term_1 > verysmall) then ! Avoid 0/0
                            argument = (2.0_wp*I1*I1*I1 - 9.0_wp*I1*I2 + 27.0_wp*I3)/ &
                                       (2.0_wp*sqrt_term_1*sqrt_term_1*sqrt_term_1)
                            if (argument > 1.0_wp) argument = 1.0_wp
                            if (argument < -1.0_wp) argument = -1.0_wp
                            phi = acos(argument)
                            sqrt_term_2 = sqrt(max(I1**2.0_wp - 3.0_wp*I2, 0.0_wp))
                            tau_p = I1/3.0_wp + 2.0_wp/sqrt(3.0_wp)*sqrt_term_2*cos(phi/3.0_wp)
                        else
                            tau_p = I1/3.0_wp
                        end if

                        rhs_vf(damage_idx)%sf(k, l, q) = (alpha_bar*max(tau_p - tau_star, 0._wp))**cont_damage_s
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_compute_damage_state

end module m_hypoelastic
