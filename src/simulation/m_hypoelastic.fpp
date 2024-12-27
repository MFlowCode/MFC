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
 s_compute_hypoelastic_rhs

    real(wp), allocatable, dimension(:) :: Gs
    !$acc declare create(Gs)

    real(wp), allocatable, dimension(:, :, :) :: du_dx, du_dy, du_dz
    real(wp), allocatable, dimension(:, :, :) :: dv_dx, dv_dy, dv_dz
    real(wp), allocatable, dimension(:, :, :) :: dw_dx, dw_dy, dw_dz
    !$acc declare create(du_dx,du_dy,du_dz,dv_dx,dv_dy,dv_dz,dw_dx,dw_dy,dw_dz)

    real(wp), allocatable, dimension(:, :, :) :: rho_K_field, G_K_field
    !$acc declare create(rho_K_field, G_K_field)

    real(wp), allocatable, dimension(:, :) :: fd_coeff_x_h
    real(wp), allocatable, dimension(:, :) :: fd_coeff_y_h
    real(wp), allocatable, dimension(:, :) :: fd_coeff_z_h
    !$acc declare create(fd_coeff_x_h,fd_coeff_y_h,fd_coeff_z_h)

contains

    subroutine s_initialize_hypoelastic_module

        integer :: i, k, r

        @:ALLOCATE(Gs(1:num_fluids))
        @:ALLOCATE(rho_K_field(0:m,0:n,0:p), G_K_field(0:m,0:n,0:p))
        @:ALLOCATE(du_dx(0:m,0:n,0:p))
        if (n > 0) then
            @:ALLOCATE(du_dy(0:m,0:n,0:p), dv_dx(0:m,0:n,0:p), dv_dy(0:m,0:n,0:p))
            if (p > 0) then
                @:ALLOCATE(du_dz(0:m,0:n,0:p), dv_dz(0:m,0:n,0:p))
                @:ALLOCATE(dw_dx(0:m,0:n,0:p), dw_dy(0:m,0:n,0:p), dw_dz(0:m,0:n,0:p))
            end if
        end if

        do i = 1, num_fluids
            Gs(i) = fluid_pp(i)%G
        end do
        !$acc update device(Gs)

        @:ALLOCATE(fd_coeff_x_h(-fd_number:fd_number, 0:m))
        if (n > 0) then
            @:ALLOCATE(fd_coeff_y_h(-fd_number:fd_number, 0:n))
        end if
        if (p > 0) then
            @:ALLOCATE(fd_coeff_z_h(-fd_number:fd_number, 0:p))
        end if

        ! Computing centered finite difference coefficients
        call s_compute_finite_difference_coefficients(m, x_cc, fd_coeff_x_h, buff_size, &
                                                      fd_number, fd_order)
        !$acc update device(fd_coeff_x_h)
        if (n > 0) then
            call s_compute_finite_difference_coefficients(n, y_cc, fd_coeff_y_h, buff_size, &
                                                          fd_number, fd_order)
            !$acc update device(fd_coeff_y_h)
        end if
        if (p > 0) then
            call s_compute_finite_difference_coefficients(p, z_cc, fd_coeff_z_h, buff_size, &
                                                          fd_number, fd_order)
            !$acc update device(fd_coeff_z_h)
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

            !$acc parallel loop collapse(3) gang vector default(present)
            do q = 0, p
                do l = 0, n
                    do k = 0, m
                        du_dx(k, l, q) = 0._wp
                    end do
                end do
            end do
            !$acc end parallel loop

            !$acc parallel loop collapse(3) gang vector default(present)
            do q = 0, p
                do l = 0, n
                    do k = 0, m
                        !$acc loop seq
                        do r = -fd_number, fd_number
                            du_dx(k, l, q) = du_dx(k, l, q) &
                                             + q_prim_vf(momxb)%sf(k + r, l, q)*fd_coeff_x_h(r, k)
                        end do

                    end do
                end do
            end do
            !$acc end parallel loop

            if (ndirs > 1) then
                !$acc parallel loop collapse(3) gang vector default(present)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            du_dy(k, l, q) = 0._wp; dv_dx(k, l, q) = 0._wp; dv_dy(k, l, q) = 0._wp
                        end do
                    end do
                end do
                !$acc end parallel loop

                !$acc parallel loop collapse(3) gang vector default(present)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            !$acc loop seq
                            do r = -fd_number, fd_number
                                du_dy(k, l, q) = du_dy(k, l, q) &
                                                 + q_prim_vf(momxb)%sf(k, l + r, q)*fd_coeff_y_h(r, l)
                                dv_dx(k, l, q) = dv_dx(k, l, q) &
                                                 + q_prim_vf(momxb + 1)%sf(k + r, l, q)*fd_coeff_x_h(r, k)
                                dv_dy(k, l, q) = dv_dy(k, l, q) &
                                                 + q_prim_vf(momxb + 1)%sf(k, l + r, q)*fd_coeff_y_h(r, l)
                            end do
                        end do
                    end do
                end do
                !$acc end parallel loop

                ! 3D
                if (ndirs == 3) then

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                du_dz(k, l, q) = 0_wp; dv_dz(k, l, q) = 0_wp; dw_dx(k, l, q) = 0_wp; 
                                dw_dy(k, l, q) = 0_wp; dw_dz(k, l, q) = 0_wp; 
                            end do
                        end do
                    end do
                    !$acc end parallel loop

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                !$acc loop seq
                                do r = -fd_number, fd_number
                                    du_dz(k, l, q) = du_dz(k, l, q) &
                                                     + q_prim_vf(momxb)%sf(k, l, q + r)*fd_coeff_z_h(r, q)
                                    dv_dz(k, l, q) = dv_dz(k, l, q) &
                                                     + q_prim_vf(momxb + 1)%sf(k, l, q + r)*fd_coeff_z_h(r, q)
                                    dw_dx(k, l, q) = dw_dx(k, l, q) &
                                                     + q_prim_vf(momxe)%sf(k + r, l, q)*fd_coeff_x_h(r, k)
                                    dw_dy(k, l, q) = dw_dy(k, l, q) &
                                                     + q_prim_vf(momxe)%sf(k, l + r, q)*fd_coeff_y_h(r, l)
                                    dw_dz(k, l, q) = dw_dz(k, l, q) &
                                                     + q_prim_vf(momxe)%sf(k, l, q + r)*fd_coeff_z_h(r, q)
                                end do
                            end do
                        end do
                    end do
                    !$acc end parallel loop
                end if
            end if

            !$acc parallel loop collapse(3) gang vector default(present)
            do q = 0, p
                do l = 0, n
                    do k = 0, m
                        rho_K = 0._wp; G_K = 0._wp
                        do i = 1, num_fluids
                            rho_K = rho_K + q_prim_vf(i)%sf(k, l, q) !alpha_rho_K(1)
                            G_K = G_K + q_prim_vf(advxb - 1 + i)%sf(k, l, q)*Gs(i)  !alpha_K(1) * Gs(1)
                        end do
                        rho_K_field(k, l, q) = rho_K
                        G_K_field(k, l, q) = G_K

                        !TODO: take this out if not needed
                        if (G_K < verysmall) then
                            G_K_field(k, l, q) = 0
                        end if
                    end do
                end do
            end do

            ! apply rhs source term to elastic stress equation
            !$acc parallel loop collapse(3) gang vector default(present)
            do q = 0, p
                do l = 0, n
                    do k = 0, m
                        rhs_vf(strxb)%sf(k, l, q) = &
                            rhs_vf(strxb)%sf(k, l, q) + rho_K_field(k, l, q)* &
                            ((4._wp*G_K_field(k, l, q)/3._wp) + &
                             q_prim_vf(strxb)%sf(k, l, q))* &
                            du_dx(k, l, q)
                    end do
                end do
            end do

        elseif (idir == 2) then
            !$acc parallel loop collapse(3) gang vector default(present)
            do q = 0, p
                do l = 0, n
                    do k = 0, m
                        rhs_vf(strxb)%sf(k, l, q) = rhs_vf(strxb)%sf(k, l, q) + rho_K_field(k, l, q)* &
                                                    (q_prim_vf(strxb + 1)%sf(k, l, q)*du_dy(k, l, q) + &
                                                     q_prim_vf(strxb + 1)%sf(k, l, q)*du_dy(k, l, q) - &
                                                     q_prim_vf(strxb)%sf(k, l, q)*dv_dy(k, l, q) - &
                                                     2._wp*G_K_field(k, l, q)*(1._wp/3._wp)*dv_dy(k, l, q))

                        rhs_vf(strxb + 1)%sf(k, l, q) = rhs_vf(strxb + 1)%sf(k, l, q) + rho_K_field(k, l, q)* &
                                                        (q_prim_vf(strxb + 1)%sf(k, l, q)*du_dx(k, l, q) + &
                                                         q_prim_vf(strxb)%sf(k, l, q)*dv_dx(k, l, q) - &
                                                         q_prim_vf(strxb + 1)%sf(k, l, q)*du_dx(k, l, q) + &
                                                         q_prim_vf(strxb + 2)%sf(k, l, q)*du_dy(k, l, q) + &
                                                         q_prim_vf(strxb + 1)%sf(k, l, q)*dv_dy(k, l, q) - &
                                                         q_prim_vf(strxb + 1)%sf(k, l, q)*dv_dy(k, l, q) + &
                                                         2._wp*G_K_field(k, l, q)*(1._wp/2._wp)*(du_dy(k, l, q) + &
                                                                                                 dv_dx(k, l, q)))

                        rhs_vf(strxb + 2)%sf(k, l, q) = rhs_vf(strxb + 2)%sf(k, l, q) + rho_K_field(k, l, q)* &
                                                        (q_prim_vf(strxb + 1)%sf(k, l, q)*dv_dx(k, l, q) + &
                                                         q_prim_vf(strxb + 1)%sf(k, l, q)*dv_dx(k, l, q) - &
                                                         q_prim_vf(strxb + 2)%sf(k, l, q)*du_dx(k, l, q) + &
                                                         q_prim_vf(strxb + 2)%sf(k, l, q)*dv_dy(k, l, q) + &
                                                         q_prim_vf(strxb + 2)%sf(k, l, q)*dv_dy(k, l, q) - &
                                                         q_prim_vf(strxb + 2)%sf(k, l, q)*dv_dy(k, l, q) + &
                                                         2._wp*G_K_field(k, l, q)*(dv_dy(k, l, q) - (1._wp/3._wp)* &
                                                                                   (du_dx(k, l, q) + &
                                                                                    dv_dy(k, l, q))))
                    end do
                end do
            end do

        elseif (idir == 3) then
            !$acc parallel loop collapse(3) gang vector default(present)
            do q = 0, p
                do l = 0, n
                    do k = 0, m
                        rhs_vf(strxb)%sf(k, l, q) = rhs_vf(strxb)%sf(k, l, q) + rho_K_field(k, l, q)* &
                                                    (q_prim_vf(strxb + 3)%sf(k, l, q)*du_dz(k, l, q) + &
                                                     q_prim_vf(strxb + 3)%sf(k, l, q)*du_dz(k, l, q) - &
                                                     q_prim_vf(strxb)%sf(k, l, q)*dw_dz(k, l, q) - &
                                                     2._wp*G_K_field(k, l, q)*(1._wp/3._wp)*dw_dz(k, l, q))

                        rhs_vf(strxb + 1)%sf(k, l, q) = rhs_vf(strxb + 1)%sf(k, l, q) + rho_K_field(k, l, q)* &
                                                        (q_prim_vf(strxb + 4)%sf(k, l, q)*du_dz(k, l, q) + &
                                                         q_prim_vf(strxb + 3)%sf(k, l, q)*dv_dz(k, l, q) - &
                                                         q_prim_vf(strxb + 1)%sf(k, l, q)*dw_dz(k, l, q))

                        rhs_vf(strxb + 2)%sf(k, l, q) = rhs_vf(strxb + 2)%sf(k, l, q) + rho_K_field(k, l, q)* &
                                                        (q_prim_vf(strxb + 4)%sf(k, l, q)*dv_dz(k, l, q) + &
                                                         q_prim_vf(strxb + 4)%sf(k, l, q)*dv_dz(k, l, q) - &
                                                         q_prim_vf(strxb + 2)%sf(k, l, q)*dw_dz(k, l, q) - &
                                                         2._wp*G_K_field(k, l, q)*(1._wp/3._wp)*dw_dz(k, l, q))

                        rhs_vf(strxb + 3)%sf(k, l, q) = rhs_vf(strxb + 3)%sf(k, l, q) + rho_K_field(k, l, q)* &
                                                        (q_prim_vf(strxb + 3)%sf(k, l, q)*du_dx(k, l, q) + &
                                                         q_prim_vf(strxb)%sf(k, l, q)*dw_dx(k, l, q) - &
                                                         q_prim_vf(strxb + 3)%sf(k, l, q)*du_dx(k, l, q) + &
                                                         q_prim_vf(strxb + 4)%sf(k, l, q)*du_dy(k, l, q) + &
                                                         q_prim_vf(strxb + 1)%sf(k, l, q)*dw_dy(k, l, q) - &
                                                         q_prim_vf(strxb + 3)%sf(k, l, q)*dv_dy(k, l, q) + &
                                                         q_prim_vf(strxb + 5)%sf(k, l, q)*du_dz(k, l, q) + &
                                                         q_prim_vf(strxb + 3)%sf(k, l, q)*dw_dz(k, l, q) - &
                                                         q_prim_vf(strxb + 3)%sf(k, l, q)*dw_dz(k, l, q) + &
                                                         2._wp*G_K_field(k, l, q)*(1._wp/2._wp)*(du_dz(k, l, q) + &
                                                                                                 dw_dx(k, l, q)))

                        rhs_vf(strxb + 4)%sf(k, l, q) = rhs_vf(strxb + 4)%sf(k, l, q) + rho_K_field(k, l, q)* &
                                                        (q_prim_vf(strxb + 3)%sf(k, l, q)*dv_dx(k, l, q) + &
                                                         q_prim_vf(strxb + 1)%sf(k, l, q)*dw_dx(k, l, q) - &
                                                         q_prim_vf(strxb + 4)%sf(k, l, q)*du_dx(k, l, q) + &
                                                         q_prim_vf(strxb + 4)%sf(k, l, q)*dv_dy(k, l, q) + &
                                                         q_prim_vf(strxb + 2)%sf(k, l, q)*dw_dy(k, l, q) - &
                                                         q_prim_vf(strxb + 4)%sf(k, l, q)*dv_dy(k, l, q) + &
                                                         q_prim_vf(strxb + 5)%sf(k, l, q)*dv_dz(k, l, q) + &
                                                         q_prim_vf(strxb + 4)%sf(k, l, q)*dw_dz(k, l, q) - &
                                                         q_prim_vf(strxb + 4)%sf(k, l, q)*dw_dz(k, l, q) + &
                                                         2._wp*G_K_field(k, l, q)*(1._wp/2._wp)*(dv_dz(k, l, q) + &
                                                                                                 dw_dy(k, l, q)))

                        rhs_vf(strxe)%sf(k, l, q) = rhs_vf(strxe)%sf(k, l, q) + rho_K_field(k, l, q)* &
                                                    (q_prim_vf(strxe - 2)%sf(k, l, q)*dw_dx(k, l, q) + &
                                                     q_prim_vf(strxe - 2)%sf(k, l, q)*dw_dx(k, l, q) - &
                                                     q_prim_vf(strxe)%sf(k, l, q)*du_dx(k, l, q) + &
                                                     q_prim_vf(strxe - 1)%sf(k, l, q)*dw_dy(k, l, q) + &
                                                     q_prim_vf(strxe - 1)%sf(k, l, q)*dw_dy(k, l, q) - &
                                                     q_prim_vf(strxe)%sf(k, l, q)*dv_dy(k, l, q) + &
                                                     q_prim_vf(strxe)%sf(k, l, q)*dw_dz(k, l, q) + &
                                                     q_prim_vf(strxe)%sf(k, l, q)*dw_dz(k, l, q) - &
                                                     q_prim_vf(strxe)%sf(k, l, q)*dw_dz(k, l, q) + &
                                                     2._wp*G_K_field(k, l, q)*(dw_dz(k, l, q) - (1._wp/3._wp)* &
                                                                               (du_dx(k, l, q) + &
                                                                                dv_dy(k, l, q) + &
                                                                                dw_dz(k, l, q))))
                    end do
                end do
            end do
        end if

    end subroutine s_compute_hypoelastic_rhs

    subroutine s_finalize_hypoelastic_module()

        @:DEALLOCATE(Gs)
        @:DEALLOCATE(rho_K_field, G_K_field)
        @:DEALLOCATE(du_dx)
        @:DEALLOCATE(fd_coeff_x_h)
        if (n > 0) then
            @:DEALLOCATE(du_dy,dv_dx,dv_dy)
            @:DEALLOCATE(fd_coeff_y_h)
            if (p > 0) then
                @:DEALLOCATE(du_dz, dv_dz, dw_dx, dw_dy, dw_dz)
                @:DEALLOCATE(fd_coeff_z_h)
            end if
        end if

    end subroutine s_finalize_hypoelastic_module

end module m_hypoelastic
