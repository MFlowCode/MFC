!>
!! @file
!! @brief Contains module m_hypoelastic

#:include 'macros.fpp'

!> @brief Computes hypoelastic stress-rate source terms and damage-state evolution
module m_hypoelastic

    use m_derived_types
    use m_global_parameters
    use m_finite_differences
    use m_helper

    implicit none

    private; public :: s_initialize_hypoelastic_module, s_finalize_hypoelastic_module, &
        & s_compute_hypoelastic_rhs_finite_diff_per_sweep, s_compute_hypoelastic_rhs_iface, &
        & s_compute_hypoelastic_rhs_axisym_geom_iface, s_compute_damage_state

    real(wp), allocatable, dimension(:) :: Gs_hypo
    $:GPU_DECLARE(create='[Gs_hypo]')

    real(wp), allocatable, dimension(:,:,:) :: du_dx_hypo, du_dy_hypo, du_dz_hypo
    real(wp), allocatable, dimension(:,:,:) :: dv_dx_hypo, dv_dy_hypo, dv_dz_hypo
    real(wp), allocatable, dimension(:,:,:) :: dw_dx_hypo, dw_dy_hypo, dw_dz_hypo
    $:GPU_DECLARE(create='[du_dx_hypo, du_dy_hypo, du_dz_hypo, dv_dx_hypo, dv_dy_hypo, dv_dz_hypo, dw_dx_hypo, dw_dy_hypo, dw_dz_hypo]')

    real(wp), allocatable, dimension(:,:,:) :: rho_K_field, G_K_field
    $:GPU_DECLARE(create='[rho_K_field, G_K_field]')

    real(wp), allocatable, dimension(:,:) :: fd_coeff_x_hypo
    real(wp), allocatable, dimension(:,:) :: fd_coeff_y_hypo
    real(wp), allocatable, dimension(:,:) :: fd_coeff_z_hypo
    $:GPU_DECLARE(create='[fd_coeff_x_hypo, fd_coeff_y_hypo, fd_coeff_z_hypo]')

contains

    !> Initialize the hypoelastic module
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
        call s_compute_finite_difference_coefficients(m, x_cc, fd_coeff_x_hypo, buff_size, fd_number, fd_order)
        $:GPU_UPDATE(device='[fd_coeff_x_hypo]')
        if (n > 0) then
            call s_compute_finite_difference_coefficients(n, y_cc, fd_coeff_y_hypo, buff_size, fd_number, fd_order)
            $:GPU_UPDATE(device='[fd_coeff_y_hypo]')
        end if
        if (p > 0) then
            call s_compute_finite_difference_coefficients(p, z_cc, fd_coeff_z_hypo, buff_size, fd_number, fd_order)
            $:GPU_UPDATE(device='[fd_coeff_z_hypo]')
        end if

    end subroutine s_initialize_hypoelastic_module

    !> Legacy FD-based hypoelastic RHS (Mode 1: HLL). Uses finite-difference velocity gradients computed from cell-centered
    !! primitive variables. Called once per direction inside the dim-split loop. Supports 1D/2D/3D Cartesian and cylindrical
    !! geometry.
    !! @param idir Dimension splitting index
    !! @param q_prim_vf Primitive variables
    !! @param rhs_vf rhs variables
    subroutine s_compute_hypoelastic_rhs_finite_diff_per_sweep(idir, q_prim_vf, rhs_vf)

        integer, intent(in)                                    :: idir
        type(scalar_field), dimension(sys_size), intent(in)    :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        real(wp)                                               :: rho_K, G_K
        integer                                                :: i, k, l, q, r  !< Loop variables
        integer                                                :: ndirs          !< Number of coordinate directions

        ndirs = 1; if (n > 0) ndirs = 2; if (p > 0) ndirs = 3

        if (idir == 1) then
            ! calculate velocity gradients + rho_K and G_K TODO: re-organize these loops one by one for GPU efficiency if possible?

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
                            du_dx_hypo(k, l, q) = du_dx_hypo(k, l, q) + q_prim_vf(eqn_idx%mom%beg)%sf(k + r, l, &
                                       & q)*fd_coeff_x_hypo(r, k)
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
                                du_dy_hypo(k, l, q) = du_dy_hypo(k, l, q) + q_prim_vf(eqn_idx%mom%beg)%sf(k, l + r, &
                                           & q)*fd_coeff_y_hypo(r, l)
                                dv_dx_hypo(k, l, q) = dv_dx_hypo(k, l, q) + q_prim_vf(eqn_idx%mom%beg + 1)%sf(k + r, l, &
                                           & q)*fd_coeff_x_hypo(r, k)
                                dv_dy_hypo(k, l, q) = dv_dy_hypo(k, l, q) + q_prim_vf(eqn_idx%mom%beg + 1)%sf(k, l + r, &
                                           & q)*fd_coeff_y_hypo(r, l)
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
                                du_dz_hypo(k, l, q) = 0._wp; dv_dz_hypo(k, l, q) = 0._wp; dw_dx_hypo(k, l, q) = 0._wp
                                dw_dy_hypo(k, l, q) = 0._wp; dw_dz_hypo(k, l, q) = 0._wp
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
                                    du_dz_hypo(k, l, q) = du_dz_hypo(k, l, q) + q_prim_vf(eqn_idx%mom%beg)%sf(k, l, &
                                               & q + r)*fd_coeff_z_hypo(r, q)
                                    dv_dz_hypo(k, l, q) = dv_dz_hypo(k, l, q) + q_prim_vf(eqn_idx%mom%beg + 1)%sf(k, l, &
                                               & q + r)*fd_coeff_z_hypo(r, q)
                                    dw_dx_hypo(k, l, q) = dw_dx_hypo(k, l, q) + q_prim_vf(eqn_idx%mom%end)%sf(k + r, l, &
                                               & q)*fd_coeff_x_hypo(r, k)
                                    dw_dy_hypo(k, l, q) = dw_dy_hypo(k, l, q) + q_prim_vf(eqn_idx%mom%end)%sf(k, l + r, &
                                               & q)*fd_coeff_y_hypo(r, l)
                                    dw_dz_hypo(k, l, q) = dw_dz_hypo(k, l, q) + q_prim_vf(eqn_idx%mom%end)%sf(k, l, &
                                               & q + r)*fd_coeff_z_hypo(r, q)
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
                            rho_K = rho_K + q_prim_vf(i)%sf(k, l, q)  ! alpha_rho_K(1)
                            G_K = G_K + q_prim_vf(eqn_idx%adv%beg - 1 + i)%sf(k, l, q)*Gs_hypo(i)  ! alpha_K(1) * Gs_hypo(1)
                        end do

                        ! Continuum damage: (1-D) scales effective stiffness, D in [0,1]
                        if (cont_damage) G_K = G_K*max((1._wp - q_prim_vf(eqn_idx%damage)%sf(k, l, q)), 0._wp)

                        rho_K_field(k, l, q) = rho_K
                        G_K_field(k, l, q) = G_K

                        if (G_K < verysmall) then
                            G_K_field(k, l, q) = 0._wp
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
                        rhs_vf(eqn_idx%stress%beg)%sf(k, l, q) = rhs_vf(eqn_idx%stress%beg)%sf(k, l, q) + rho_K_field(k, l, &
                               & q)*((4._wp*G_K_field(k, l, q)/3._wp) + q_prim_vf(eqn_idx%stress%beg)%sf(k, l, q))*du_dx_hypo(k, &
                               & l, q)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        else if (idir == 2) then
            $:GPU_PARALLEL_LOOP(collapse=3)
            do q = 0, p
                do l = 0, n
                    do k = 0, m
                        rhs_vf(eqn_idx%stress%beg)%sf(k, l, q) = rhs_vf(eqn_idx%stress%beg)%sf(k, l, q) + rho_K_field(k, l, &
                               & q)*(2._wp*q_prim_vf(eqn_idx%stress%beg + 1)%sf(k, l, q)*du_dy_hypo(k, l, &
                               & q) - q_prim_vf(eqn_idx%stress%beg)%sf(k, l, q)*dv_dy_hypo(k, l, q) - (2._wp/3._wp)*G_K_field(k, &
                               & l, q)*dv_dy_hypo(k, l, q))

                        rhs_vf(eqn_idx%stress%beg + 1)%sf(k, l, q) = rhs_vf(eqn_idx%stress%beg + 1)%sf(k, l, q) + rho_K_field(k, &
                               & l, q)*(q_prim_vf(eqn_idx%stress%beg)%sf(k, l, q)*dv_dx_hypo(k, l, &
                               & q) + q_prim_vf(eqn_idx%stress%beg + 2)%sf(k, l, q)*du_dy_hypo(k, l, q) + G_K_field(k, l, &
                               & q)*(du_dy_hypo(k, l, q) + dv_dx_hypo(k, l, q)))

                        rhs_vf(eqn_idx%stress%beg + 2)%sf(k, l, q) = rhs_vf(eqn_idx%stress%beg + 2)%sf(k, l, q) + rho_K_field(k, &
                               & l, q)*(2._wp*q_prim_vf(eqn_idx%stress%beg + 1)%sf(k, l, q)*dv_dx_hypo(k, l, &
                               & q) - q_prim_vf(eqn_idx%stress%beg + 2)%sf(k, l, q)*du_dx_hypo(k, l, &
                               & q) + q_prim_vf(eqn_idx%stress%beg + 2)%sf(k, l, q)*dv_dy_hypo(k, l, q) + 2._wp*G_K_field(k, l, &
                               & q)*(dv_dy_hypo(k, l, q) - (1._wp/3._wp)*(du_dx_hypo(k, l, q) + dv_dy_hypo(k, l, q))))
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        else if (idir == 3) then
            $:GPU_PARALLEL_LOOP(collapse=3)
            do q = 0, p
                do l = 0, n
                    do k = 0, m
                        rhs_vf(eqn_idx%stress%beg)%sf(k, l, q) = rhs_vf(eqn_idx%stress%beg)%sf(k, l, q) + rho_K_field(k, l, &
                               & q)*(2._wp*q_prim_vf(eqn_idx%stress%beg + 3)%sf(k, l, q)*du_dz_hypo(k, l, &
                               & q) - q_prim_vf(eqn_idx%stress%beg)%sf(k, l, q)*dw_dz_hypo(k, l, q) - (2._wp/3._wp)*G_K_field(k, &
                               & l, q)*dw_dz_hypo(k, l, q))

                        rhs_vf(eqn_idx%stress%beg + 1)%sf(k, l, q) = rhs_vf(eqn_idx%stress%beg + 1)%sf(k, l, q) + rho_K_field(k, &
                               & l, q)*(q_prim_vf(eqn_idx%stress%beg + 4)%sf(k, l, q)*du_dz_hypo(k, l, &
                               & q) + q_prim_vf(eqn_idx%stress%beg + 3)%sf(k, l, q)*dv_dz_hypo(k, l, &
                               & q) - q_prim_vf(eqn_idx%stress%beg + 1)%sf(k, l, q)*dw_dz_hypo(k, l, q))

                        rhs_vf(eqn_idx%stress%beg + 2)%sf(k, l, q) = rhs_vf(eqn_idx%stress%beg + 2)%sf(k, l, q) + rho_K_field(k, &
                               & l, q)*(2._wp*q_prim_vf(eqn_idx%stress%beg + 4)%sf(k, l, q)*dv_dz_hypo(k, l, &
                               & q) - q_prim_vf(eqn_idx%stress%beg + 2)%sf(k, l, q)*dw_dz_hypo(k, l, &
                               & q) - (2._wp/3._wp)*G_K_field(k, l, q)*dw_dz_hypo(k, l, q))

                        rhs_vf(eqn_idx%stress%beg + 3)%sf(k, l, q) = rhs_vf(eqn_idx%stress%beg + 3)%sf(k, l, q) + rho_K_field(k, &
                               & l, q)*(q_prim_vf(eqn_idx%stress%beg)%sf(k, l, q)*dw_dx_hypo(k, l, &
                               & q) + q_prim_vf(eqn_idx%stress%beg + 4)%sf(k, l, q)*du_dy_hypo(k, l, &
                               & q) + q_prim_vf(eqn_idx%stress%beg + 1)%sf(k, l, q)*dw_dy_hypo(k, l, &
                               & q) - q_prim_vf(eqn_idx%stress%beg + 3)%sf(k, l, q)*dv_dy_hypo(k, l, &
                               & q) + q_prim_vf(eqn_idx%stress%beg + 5)%sf(k, l, q)*du_dz_hypo(k, l, q) + G_K_field(k, l, &
                               & q)*(du_dz_hypo(k, l, q) + dw_dx_hypo(k, l, q)))

                        rhs_vf(eqn_idx%stress%beg + 4)%sf(k, l, q) = rhs_vf(eqn_idx%stress%beg + 4)%sf(k, l, q) + rho_K_field(k, &
                               & l, q)*(q_prim_vf(eqn_idx%stress%beg + 3)%sf(k, l, q)*dv_dx_hypo(k, l, &
                               & q) + q_prim_vf(eqn_idx%stress%beg + 1)%sf(k, l, q)*dw_dx_hypo(k, l, &
                               & q) - q_prim_vf(eqn_idx%stress%beg + 4)%sf(k, l, q)*du_dx_hypo(k, l, &
                               & q) + q_prim_vf(eqn_idx%stress%beg + 2)%sf(k, l, q)*dw_dy_hypo(k, l, &
                               & q) + q_prim_vf(eqn_idx%stress%beg + 5)%sf(k, l, q)*dv_dz_hypo(k, l, q) + G_K_field(k, l, &
                               & q)*(dv_dz_hypo(k, l, q) + dw_dy_hypo(k, l, q)))

                        rhs_vf(eqn_idx%stress%end)%sf(k, l, q) = rhs_vf(eqn_idx%stress%end)%sf(k, l, q) + rho_K_field(k, l, &
                               & q)*(2._wp*q_prim_vf(eqn_idx%stress%end - 2)%sf(k, l, q)*dw_dx_hypo(k, l, &
                               & q) - q_prim_vf(eqn_idx%stress%end)%sf(k, l, q)*du_dx_hypo(k, l, &
                               & q) + 2._wp*q_prim_vf(eqn_idx%stress%end - 1)%sf(k, l, q)*dw_dy_hypo(k, l, &
                               & q) - q_prim_vf(eqn_idx%stress%end)%sf(k, l, q)*dv_dy_hypo(k, l, &
                               & q) + q_prim_vf(eqn_idx%stress%end)%sf(k, l, q)*dw_dz_hypo(k, l, q) + 2._wp*G_K_field(k, l, &
                               & q)*(dw_dz_hypo(k, l, q) - (1._wp/3._wp)*(du_dx_hypo(k, l, q) + dv_dy_hypo(k, l, &
                               & q) + dw_dz_hypo(k, l, q))))
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
                        rhs_vf(eqn_idx%stress%beg)%sf(k, l, q) = rhs_vf(eqn_idx%stress%beg)%sf(k, l, q) - rho_K_field(k, l, &
                               & q)*q_prim_vf(eqn_idx%mom%beg + 1)%sf(k, l, q)/y_cc(l)*(q_prim_vf(eqn_idx%stress%beg)%sf(k, l, &
                               & q) + (2._wp/3._wp)*G_K_field(k, l, q))  ! tau_xx + 2/3*G

                        ! S_xr -= rho * v/r * tau_xr
                        rhs_vf(eqn_idx%stress%beg + 1)%sf(k, l, q) = rhs_vf(eqn_idx%stress%beg + 1)%sf(k, l, q) - rho_K_field(k, &
                               & l, q)*q_prim_vf(eqn_idx%mom%beg + 1)%sf(k, l, q)/y_cc(l)*q_prim_vf(eqn_idx%stress%beg + 1)%sf(k, &
                               & l, q)  ! tau_xx

                        ! S_rr -= rho * v/r * (tau_rr + 2/3*G)
                        rhs_vf(eqn_idx%stress%beg + 2)%sf(k, l, q) = rhs_vf(eqn_idx%stress%beg + 2)%sf(k, l, q) - rho_K_field(k, &
                               & l, q)*q_prim_vf(eqn_idx%mom%beg + 1)%sf(k, l, &
                               & q)/y_cc(l)*(q_prim_vf(eqn_idx%stress%beg + 2)%sf(k, l, q) + (2._wp/3._wp)*G_K_field(k, l, &
                               & q))  ! tau_rr + 2/3*G

                        ! S_thetatheta += rho * ( -(tau_thetatheta + 2/3*G)*(du/dx + dv/dr + v/r) + 2*(tau_thetatheta + G)*v/r )
                        rhs_vf(eqn_idx%stress%beg + 3)%sf(k, l, q) = rhs_vf(eqn_idx%stress%beg + 3)%sf(k, l, q) + rho_K_field(k, &
                               & l, q)*(-(q_prim_vf(eqn_idx%stress%beg + 3)%sf(k, l, q) + (2._wp/3._wp)*G_K_field(k, l, &
                               & q))*(du_dx_hypo(k, l, q) + dv_dy_hypo(k, l, q) + q_prim_vf(eqn_idx%mom%beg + 1)%sf(k, l, &
                               & q)/y_cc(l)) + 2._wp*(q_prim_vf(eqn_idx%stress%beg + 3)%sf(k, l, q) + G_K_field(k, l, &
                               & q))*q_prim_vf(eqn_idx%mom%beg + 1)%sf(k, l, q)/y_cc(l))
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_compute_hypoelastic_rhs_finite_diff_per_sweep

    !> Interface-consistent hypoelastic RHS (Mode 2: HLL/HLLC). Uses interface velocities from the Riemann solver to compute
    !! velocity gradients. Called once after all dimensional sweeps. Supports 1D, 2D Cartesian, 2D axisymmetric, and 3D Cartesian.
    !! @param q_prim_vf Primitive variables
    !! @param rhs_vf rhs variables
    !! @param nc_iface_vel_n Interface velocities per direction
    subroutine s_compute_hypoelastic_rhs_iface(q_prim_vf, rhs_vf, nc_iface_vel_n)

        type(scalar_field), dimension(sys_size), intent(in)    :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        type(vector_field), dimension(:), intent(in)           :: nc_iface_vel_n
        real(wp)                                               :: rho_K, G_K
        real(wp)                                               :: trace, shear, shear2, diag, diag_z, offdiag, cross1, cross2
        real(wp)                                               :: txx, txy, tyy, txz, tyz, tzz
        integer                                                :: i, k, l, q, r
        integer                                                :: ndirs

        ndirs = 1; if (n > 0) ndirs = 2; if (p > 0) ndirs = 3

        $:GPU_PARALLEL_LOOP(collapse=3)
        do q = 0, p
            do l = 0, n
                do k = 0, m
                    du_dx_hypo(k, l, q) = (nc_iface_vel_n(1)%vf(1)%sf(k, l, q) - nc_iface_vel_n(1)%vf(1)%sf(k - 1, l, q))/dx(k)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        if (ndirs > 1) then
            $:GPU_PARALLEL_LOOP(collapse=3)
            do q = 0, p
                do l = 0, n
                    do k = 0, m
                        du_dy_hypo(k, l, q) = (nc_iface_vel_n(2)%vf(1)%sf(k, l, q) - nc_iface_vel_n(2)%vf(1)%sf(k, l - 1, q))/dy(l)
                        dv_dx_hypo(k, l, q) = (nc_iface_vel_n(1)%vf(2)%sf(k, l, q) - nc_iface_vel_n(1)%vf(2)%sf(k - 1, l, q))/dx(k)
                        dv_dy_hypo(k, l, q) = (nc_iface_vel_n(2)%vf(2)%sf(k, l, q) - nc_iface_vel_n(2)%vf(2)%sf(k, l - 1, q))/dy(l)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        if (ndirs == 3) then
            $:GPU_PARALLEL_LOOP(collapse=3)
            do q = 0, p
                do l = 0, n
                    do k = 0, m
                        du_dz_hypo(k, l, q) = (nc_iface_vel_n(3)%vf(1)%sf(k, l, q) - nc_iface_vel_n(3)%vf(1)%sf(k, l, q - 1))/dz(q)
                        dv_dz_hypo(k, l, q) = (nc_iface_vel_n(3)%vf(2)%sf(k, l, q) - nc_iface_vel_n(3)%vf(2)%sf(k, l, q - 1))/dz(q)
                        dw_dx_hypo(k, l, q) = (nc_iface_vel_n(1)%vf(3)%sf(k, l, q) - nc_iface_vel_n(1)%vf(3)%sf(k - 1, l, q))/dx(k)
                        dw_dy_hypo(k, l, q) = (nc_iface_vel_n(2)%vf(3)%sf(k, l, q) - nc_iface_vel_n(2)%vf(3)%sf(k, l - 1, q))/dy(l)
                        dw_dz_hypo(k, l, q) = (nc_iface_vel_n(3)%vf(3)%sf(k, l, q) - nc_iface_vel_n(3)%vf(3)%sf(k, l, q - 1))/dz(q)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        $:GPU_PARALLEL_LOOP(collapse=3,private='[rho_K, G_K]')
        do q = 0, p
            do l = 0, n
                do k = 0, m
                    rho_K = 0._wp; G_K = 0._wp
                    do i = 1, num_fluids
                        rho_K = rho_K + q_prim_vf(i)%sf(k, l, q)
                        G_K = G_K + q_prim_vf(eqn_idx%adv%beg - 1 + i)%sf(k, l, q)*Gs_hypo(i)
                    end do

                    if (cont_damage) G_K = G_K*max((1._wp - q_prim_vf(eqn_idx%damage)%sf(k, l, q)), 0._wp)

                    rho_K_field(k, l, q) = rho_K
                    G_K_field(k, l, q) = G_K

                    if (G_K < verysmall) then
                        G_K_field(k, l, q) = 0._wp
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        $:GPU_PARALLEL_LOOP(collapse=3)
        do q = 0, p
            do l = 0, n
                do k = 0, m
                    rhs_vf(eqn_idx%stress%beg)%sf(k, l, q) = rhs_vf(eqn_idx%stress%beg)%sf(k, l, q) + rho_K_field(k, l, &
                           & q)*((4._wp*G_K_field(k, l, q)/3._wp) + q_prim_vf(eqn_idx%stress%beg)%sf(k, l, q))*du_dx_hypo(k, l, q)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        if (ndirs > 1) then
            $:GPU_PARALLEL_LOOP(collapse=3)
            do q = 0, p
                do l = 0, n
                    do k = 0, m
                        rhs_vf(eqn_idx%stress%beg)%sf(k, l, q) = rhs_vf(eqn_idx%stress%beg)%sf(k, l, q) + rho_K_field(k, l, &
                               & q)*(2._wp*q_prim_vf(eqn_idx%stress%beg + 1)%sf(k, l, q)*du_dy_hypo(k, l, &
                               & q) - q_prim_vf(eqn_idx%stress%beg)%sf(k, l, q)*dv_dy_hypo(k, l, q) - (2._wp/3._wp)*G_K_field(k, &
                               & l, q)*dv_dy_hypo(k, l, q))

                        rhs_vf(eqn_idx%stress%beg + 1)%sf(k, l, q) = rhs_vf(eqn_idx%stress%beg + 1)%sf(k, l, q) + rho_K_field(k, &
                               & l, q)*(q_prim_vf(eqn_idx%stress%beg)%sf(k, l, q)*dv_dx_hypo(k, l, &
                               & q) + q_prim_vf(eqn_idx%stress%beg + 2)%sf(k, l, q)*du_dy_hypo(k, l, q) + G_K_field(k, l, &
                               & q)*(du_dy_hypo(k, l, q) + dv_dx_hypo(k, l, q)))

                        rhs_vf(eqn_idx%stress%beg + 2)%sf(k, l, q) = rhs_vf(eqn_idx%stress%beg + 2)%sf(k, l, q) + rho_K_field(k, &
                               & l, q)*(2._wp*q_prim_vf(eqn_idx%stress%beg + 1)%sf(k, l, q)*dv_dx_hypo(k, l, &
                               & q) - q_prim_vf(eqn_idx%stress%beg + 2)%sf(k, l, q)*du_dx_hypo(k, l, &
                               & q) + q_prim_vf(eqn_idx%stress%beg + 2)%sf(k, l, q)*dv_dy_hypo(k, l, q) + 2._wp*G_K_field(k, l, &
                               & q)*(dv_dy_hypo(k, l, q) - (1._wp/3._wp)*(du_dx_hypo(k, l, q) + dv_dy_hypo(k, l, q))))
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        if (ndirs == 3 .and. .not. cyl_coord) then
            $:GPU_PARALLEL_LOOP(collapse=3,private='[txx, txy, tyy, txz, tyz, tzz, trace, shear, shear2, diag, diag_z, offdiag, &
                                & cross1, cross2]')
            do q = 0, p
                do l = 0, n
                    do k = 0, m
                        txx = q_prim_vf(eqn_idx%stress%beg)%sf(k, l, q)
                        txy = q_prim_vf(eqn_idx%stress%beg + 1)%sf(k, l, q)
                        tyy = q_prim_vf(eqn_idx%stress%beg + 2)%sf(k, l, q)
                        txz = q_prim_vf(eqn_idx%stress%beg + 3)%sf(k, l, q)
                        tyz = q_prim_vf(eqn_idx%stress%beg + 4)%sf(k, l, q)
                        tzz = q_prim_vf(eqn_idx%stress%beg + 5)%sf(k, l, q)

                        ! z-direction contributions to tau_xx
                        trace = -(2._wp/3._wp*G_K_field(k, l, q) + txx)
                        shear = 2._wp*txz
                        rhs_vf(eqn_idx%stress%beg)%sf(k, l, q) = rhs_vf(eqn_idx%stress%beg)%sf(k, l, q) + rho_K_field(k, l, &
                               & q)*(trace*dw_dz_hypo(k, l, q) + shear*du_dz_hypo(k, l, q))

                        ! z-direction contributions to tau_xy
                        offdiag = -txy
                        cross1 = tyz
                        cross2 = txz
                        rhs_vf(eqn_idx%stress%beg + 1)%sf(k, l, q) = rhs_vf(eqn_idx%stress%beg + 1)%sf(k, l, q) + rho_K_field(k, &
                               & l, q)*(offdiag*dw_dz_hypo(k, l, q) + cross1*du_dz_hypo(k, l, q) + cross2*dv_dz_hypo(k, l, q))

                        ! z-direction contributions to tau_yy
                        trace = -(2._wp/3._wp*G_K_field(k, l, q) + tyy)
                        shear = 2._wp*tyz
                        rhs_vf(eqn_idx%stress%beg + 2)%sf(k, l, q) = rhs_vf(eqn_idx%stress%beg + 2)%sf(k, l, q) + rho_K_field(k, &
                               & l, q)*(trace*dw_dz_hypo(k, l, q) + shear*dv_dz_hypo(k, l, q))

                        ! tau_xz (stress%beg+3)
                        diag = G_K_field(k, l, q) + txx
                        offdiag = -txz
                        cross1 = tyz
                        cross2 = txy
                        diag_z = G_K_field(k, l, q) + tzz
                        rhs_vf(eqn_idx%stress%beg + 3)%sf(k, l, q) = rhs_vf(eqn_idx%stress%beg + 3)%sf(k, l, q) + rho_K_field(k, &
                               & l, q)*(diag*dw_dx_hypo(k, l, q) + offdiag*dv_dy_hypo(k, l, q) + cross1*du_dy_hypo(k, l, &
                               & q) + cross2*dw_dy_hypo(k, l, q) + diag_z*du_dz_hypo(k, l, q))

                        ! tau_yz (stress%beg+4)
                        offdiag = -tyz
                        cross1 = txz
                        cross2 = txy
                        diag = G_K_field(k, l, q) + tyy
                        rhs_vf(eqn_idx%stress%beg + 4)%sf(k, l, q) = rhs_vf(eqn_idx%stress%beg + 4)%sf(k, l, q) + rho_K_field(k, &
                               & l, q)*(offdiag*du_dx_hypo(k, l, q) + cross1*dv_dx_hypo(k, l, q) + cross2*dw_dx_hypo(k, l, &
                               & q) + diag*dw_dy_hypo(k, l, q) + diag_z*dv_dz_hypo(k, l, q))

                        ! tau_zz (stress%beg+5)
                        trace = -(2._wp/3._wp*G_K_field(k, l, q) + tzz)
                        shear = 2._wp*txz
                        shear2 = 2._wp*tyz
                        diag = 4._wp/3._wp*G_K_field(k, l, q) + tzz
                        rhs_vf(eqn_idx%stress%beg + 5)%sf(k, l, q) = rhs_vf(eqn_idx%stress%beg + 5)%sf(k, l, q) + rho_K_field(k, &
                               & l, q)*(trace*du_dx_hypo(k, l, q) + shear*dw_dx_hypo(k, l, q) + trace*dv_dy_hypo(k, l, &
                               & q) + shear2*dw_dy_hypo(k, l, q) + diag*dw_dz_hypo(k, l, q))
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        if (grid_geometry == 2) then
            call s_compute_hypoelastic_rhs_axisym_geom_iface(q_prim_vf, rhs_vf, nc_iface_vel_n(1)%vf, nc_iface_vel_n(2)%vf, 1.0_wp)
        end if

    end subroutine s_compute_hypoelastic_rhs_iface

    !> Axisymmetric geometric source terms for the hypoelastic stress evolution, using interface velocities. Adds the v/r and div(u)
    !! contributions that arise in cylindrical (r-z) coordinates: tau_xx, tau_xr, tau_rr get a -rho*(v/r) source; tau_thetatheta
    !! gets a combined divergence and hoop-stress source. Called from s_compute_hypoelastic_rhs_iface when grid_geometry == 2.
    !! @param q_prim_vf Primitive variables
    !! @param rhs_vf rhs variables
    !! @param nc_iface_vel_x_vf Interface velocities in x-direction
    !! @param nc_iface_vel_y_vf Interface velocities in y-direction
    !! @param weight Sub-step weighting factor
    subroutine s_compute_hypoelastic_rhs_axisym_geom_iface(q_prim_vf, rhs_vf, nc_iface_vel_x_vf, nc_iface_vel_y_vf, weight)

        type(scalar_field), dimension(sys_size), intent(in)    :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        type(scalar_field), dimension(:), intent(in)           :: nc_iface_vel_x_vf
        type(scalar_field), dimension(:), intent(in)           :: nc_iface_vel_y_vf
        real(wp), intent(in)                                   :: weight
        integer                                                :: i, k, l, q
        real(wp)                                               :: rho_K, G_K, v_over_r, divU_axi

        $:GPU_PARALLEL_LOOP(collapse=3)
        do q = 0, p
            do l = 0, n
                do k = 0, m
                    du_dx_hypo(k, l, q) = (nc_iface_vel_x_vf(1)%sf(k, l, q) - nc_iface_vel_x_vf(1)%sf(k - 1, l, q))/dx(k)

                    dv_dy_hypo(k, l, q) = (nc_iface_vel_y_vf(2)%sf(k, l, q) - nc_iface_vel_y_vf(2)%sf(k, l - 1, q))/dy(l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        $:GPU_PARALLEL_LOOP(collapse=3,private='[rho_K, G_K, v_over_r, divU_axi]')
        do q = 0, p
            do l = 0, n
                do k = 0, m
                    rho_K = 0._wp
                    G_K = 0._wp
                    do i = 1, num_fluids
                        rho_K = rho_K + q_prim_vf(i)%sf(k, l, q)
                        G_K = G_K + q_prim_vf(eqn_idx%adv%beg - 1 + i)%sf(k, l, q)*Gs_hypo(i)
                    end do

                    if (cont_damage) G_K = G_K*max(1._wp - q_prim_vf(eqn_idx%damage)%sf(k, l, q), 0._wp)

                    v_over_r = q_prim_vf(eqn_idx%mom%beg + 1)%sf(k, l, q)/y_cc(l)
                    divU_axi = du_dx_hypo(k, l, q) + dv_dy_hypo(k, l, q) + v_over_r

                    rhs_vf(eqn_idx%stress%beg)%sf(k, l, q) = rhs_vf(eqn_idx%stress%beg)%sf(k, l, &
                           & q) - weight*rho_K*v_over_r*(q_prim_vf(eqn_idx%stress%beg)%sf(k, l, q) + 2._wp*G_K/3._wp)

                    rhs_vf(eqn_idx%stress%beg + 1)%sf(k, l, q) = rhs_vf(eqn_idx%stress%beg + 1)%sf(k, l, &
                           & q) - weight*rho_K*v_over_r*q_prim_vf(eqn_idx%stress%beg + 1)%sf(k, l, q)

                    rhs_vf(eqn_idx%stress%beg + 2)%sf(k, l, q) = rhs_vf(eqn_idx%stress%beg + 2)%sf(k, l, &
                           & q) - weight*rho_K*v_over_r*(q_prim_vf(eqn_idx%stress%beg + 2)%sf(k, l, q) + 2._wp*G_K/3._wp)

                    rhs_vf(eqn_idx%stress%beg + 3)%sf(k, l, q) = rhs_vf(eqn_idx%stress%beg + 3)%sf(k, l, &
                           & q) + weight*rho_K*(-(q_prim_vf(eqn_idx%stress%beg + 3)%sf(k, l, &
                           & q) + 2._wp*G_K/3._wp)*divU_axi + 2._wp*(q_prim_vf(eqn_idx%stress%beg + 3)%sf(k, l, q) + G_K)*v_over_r)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_hypoelastic_rhs_axisym_geom_iface

    !> Finalize the hypoelastic module
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

    !> Compute the continuum damage source term from the principal stress state
    subroutine s_compute_damage_state(q_cons_vf, rhs_vf)

        type(scalar_field), dimension(sys_size), intent(in)    :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        real(wp)                                               :: tau_p  !< principal stress
        real(wp)                                               :: tau_xx, tau_xy, tau_yy, tau_zz, tau_yz, tau_xz
        real(wp)                                               :: I1, I2, I3, argument, phi, sqrt_term_1, sqrt_term_2, temp
        integer                                                :: q, l, k

        if (n == 0) then
            l = 0; q = 0
            $:GPU_PARALLEL_LOOP()
            do k = 0, m
                rhs_vf(eqn_idx%damage)%sf(k, l, q) = (alpha_bar*max(abs(real(q_cons_vf(eqn_idx%stress%beg)%sf(k, l, q), &
                       & kind=wp)) - tau_star, 0._wp))**cont_damage_s
            end do
            $:END_GPU_PARALLEL_LOOP()
        else if (p == 0) then
            q = 0
            $:GPU_PARALLEL_LOOP(collapse=2, private='[tau_p]')
            do l = 0, n
                do k = 0, m
                    ! Maximum principal stress
                    tau_p = 0.5_wp*(q_cons_vf(eqn_idx%stress%beg)%sf(k, l, q) + q_cons_vf(eqn_idx%stress%beg + 2)%sf(k, l, &
                                    & q)) + sqrt((q_cons_vf(eqn_idx%stress%beg)%sf(k, l, &
                                    & q) - q_cons_vf(eqn_idx%stress%beg + 2)%sf(k, l, &
                                    & q))**2.0_wp + 4._wp*q_cons_vf(eqn_idx%stress%beg + 1)%sf(k, l, q)**2.0_wp)/2._wp

                    rhs_vf(eqn_idx%damage)%sf(k, l, q) = (alpha_bar*max(tau_p - tau_star, 0._wp))**cont_damage_s
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        else
            $:GPU_PARALLEL_LOOP(collapse=3, private='[tau_xx, tau_xy, tau_yy, tau_xz, tau_yz, tau_zz, I1, I2, I3, temp, &
                                & sqrt_term_1, sqrt_term_2, argument, phi, tau_p]')
            do q = 0, p
                do l = 0, n
                    do k = 0, m
                        tau_xx = q_cons_vf(eqn_idx%stress%beg)%sf(k, l, q)
                        tau_xy = q_cons_vf(eqn_idx%stress%beg + 1)%sf(k, l, q)
                        tau_yy = q_cons_vf(eqn_idx%stress%beg + 2)%sf(k, l, q)
                        tau_xz = q_cons_vf(eqn_idx%stress%beg + 3)%sf(k, l, q)
                        tau_yz = q_cons_vf(eqn_idx%stress%beg + 4)%sf(k, l, q)
                        tau_zz = q_cons_vf(eqn_idx%stress%beg + 5)%sf(k, l, q)

                        ! Invariants of the stress tensor
                        I1 = tau_xx + tau_yy + tau_zz
                        I2 = tau_xx*tau_yy + tau_xx*tau_zz + tau_yy*tau_zz - (tau_xy**2.0_wp + tau_xz**2.0_wp + tau_yz**2.0_wp)
                        I3 = tau_xx*tau_yy*tau_zz + 2.0_wp*tau_xy*tau_xz*tau_yz - tau_xx*tau_yz**2.0_wp - tau_yy*tau_xz**2.0_wp &
                            & - tau_zz*tau_xy**2.0_wp

                        ! Maximum principal stress
                        temp = I1**2.0_wp - 3.0_wp*I2
                        sqrt_term_1 = sqrt(max(temp, 0.0_wp))
                        if (sqrt_term_1 > verysmall) then  ! Avoid 0/0
                            argument = (2.0_wp*I1*I1*I1 - 9.0_wp*I1*I2 + 27.0_wp*I3)/(2.0_wp*sqrt_term_1*sqrt_term_1*sqrt_term_1)
                            if (argument > 1.0_wp) argument = 1.0_wp
                            if (argument < -1.0_wp) argument = -1.0_wp
                            phi = acos(argument)
                            sqrt_term_2 = sqrt(max(I1**2.0_wp - 3.0_wp*I2, 0.0_wp))
                            tau_p = I1/3.0_wp + 2.0_wp/sqrt(3.0_wp)*sqrt_term_2*cos(phi/3.0_wp)
                        else
                            tau_p = I1/3.0_wp
                        end if

                        rhs_vf(eqn_idx%damage)%sf(k, l, q) = (alpha_bar*max(tau_p - tau_star, 0._wp))**cont_damage_s
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_compute_damage_state

end module m_hypoelastic
