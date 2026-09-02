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
    use m_variables_conversion, only: f_bulk_modulus

    implicit none

    private; public :: s_initialize_hypoelastic_module, s_finalize_hypoelastic_module, &
        & s_compute_hypoelastic_rhs_finite_diff_per_sweep, s_compute_hypoelastic_rhs_iface, &
        & s_compute_hypoelastic_rhs_axisym_geom_iface, s_compute_hypoelastic_rhs_axisym_geom_dual_pass, s_compute_damage_state, &
        & s_enforce_cont_damage_bounds

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
                               & q)/y_cc(l)*(q_prim_vf(eqn_idx%stress%beg + 2)%sf(k, l, q) + (2._wp/3._wp)*G_K_field(k, l, q))  ! tau_rr + 2/3*G

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
        integer                                                :: i, k, l, q
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
            call s_compute_hypoelastic_rhs_axisym_geom_iface(q_prim_vf, rhs_vf, nc_iface_vel_n(1)%vf, nc_iface_vel_n(2)%vf)
        end if

    end subroutine s_compute_hypoelastic_rhs_iface

    !> Axisymmetric geometric source terms for the hypoelastic stress evolution, using interface velocities. Adds the v/r and div(u)
    !! contributions that arise in cylindrical (r-z) coordinates: tau_xx, tau_xr, tau_rr get a -rho*(v/r) source; tau_thetatheta
    !! gets a combined divergence and hoop-stress source. Called from s_compute_hypoelastic_rhs_iface when grid_geometry == 2.
    !! @param q_prim_vf Primitive variables
    !! @param rhs_vf rhs variables
    !! @param nc_iface_vel_x_vf Interface velocities in x-direction
    !! @param nc_iface_vel_y_vf Interface velocities in y-direction
    subroutine s_compute_hypoelastic_rhs_axisym_geom_iface(q_prim_vf, rhs_vf, nc_iface_vel_x_vf, nc_iface_vel_y_vf)

        type(scalar_field), dimension(sys_size), intent(in)    :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        type(scalar_field), dimension(:), intent(in)           :: nc_iface_vel_x_vf
        type(scalar_field), dimension(:), intent(in)           :: nc_iface_vel_y_vf
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
                           & q) - rho_K*v_over_r*(q_prim_vf(eqn_idx%stress%beg)%sf(k, l, q) + 2._wp*G_K/3._wp)

                    rhs_vf(eqn_idx%stress%beg + 1)%sf(k, l, q) = rhs_vf(eqn_idx%stress%beg + 1)%sf(k, l, &
                           & q) - rho_K*v_over_r*q_prim_vf(eqn_idx%stress%beg + 1)%sf(k, l, q)

                    rhs_vf(eqn_idx%stress%beg + 2)%sf(k, l, q) = rhs_vf(eqn_idx%stress%beg + 2)%sf(k, l, &
                           & q) - rho_K*v_over_r*(q_prim_vf(eqn_idx%stress%beg + 2)%sf(k, l, q) + 2._wp*G_K/3._wp)

                    rhs_vf(eqn_idx%stress%beg + 3)%sf(k, l, q) = rhs_vf(eqn_idx%stress%beg + 3)%sf(k, l, &
                           & q) + rho_K*(-(q_prim_vf(eqn_idx%stress%beg + 3)%sf(k, l, &
                           & q) + 2._wp*G_K/3._wp)*divU_axi + 2._wp*(q_prim_vf(eqn_idx%stress%beg + 3)%sf(k, l, q) + G_K)*v_over_r)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_hypoelastic_rhs_axisym_geom_iface

    !> Cylindrical completion for the dual-pass (anchored HLLD) hypoelastic path. The anchored augmented fluxes already carry every
    !! axial/radial derivative term of the stress law and the stress rows of flux_gsrc are zero, so the complete remaining
    !! cylindrical physics is the cell-local v/r family below: the advective metric -q_s*v/r plus the constitutive v/r terms (and
    !! +/-K*C for the volume fractions under alt_soundspeed). The discrete C = v/r averages the cell's own two anchored radial face
    !! traces (hat_L outer face, hat_R inner face) over y_cc, so no absolute axial velocity enters and uniform axial translation
    !! gives exactly zero. Called once after the two anchored partial RHS's are summed. Continuum damage needs no handling here:
    !! HLLD + cont_damage is prohibited (m_checker.fpp).
    !! @param q_prim_vf Primitive variables
    !! @param rhs_vf rhs variables
    !! @param nc_iface_vel_y_vf hat_L-pass radial-direction interface velocities
    !! @param nc_iface_vel_y_hatR_vf hat_R-pass radial-direction interface velocities
    subroutine s_compute_hypoelastic_rhs_axisym_geom_dual_pass(q_prim_vf, rhs_vf, nc_iface_vel_y_vf, nc_iface_vel_y_hatR_vf)

        type(scalar_field), dimension(sys_size), intent(in)    :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        type(scalar_field), dimension(:), intent(in)           :: nc_iface_vel_y_vf
        type(scalar_field), dimension(:), intent(in)           :: nc_iface_vel_y_hatR_vf
        real(wp)                                               :: rho_K, G_K, K_K, C_num, pres_K, blkmod1_K, blkmod2_K
        integer                                                :: i, k, l, q

        $:GPU_PARALLEL_LOOP(collapse=3, private='[rho_K, G_K, K_K, C_num, pres_K, blkmod1_K, blkmod2_K]')
        do q = 0, p
            do l = 0, n
                do k = 0, m
                    rho_K = 0._wp
                    G_K = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        rho_K = rho_K + q_prim_vf(i)%sf(k, l, q)
                        G_K = G_K + q_prim_vf(eqn_idx%adv%beg - 1 + i)%sf(k, l, q)*Gs_hypo(i)
                    end do

                    if (G_K < verysmall) G_K = 0._wp

                    ! Cell-owned anchored radial face traces: hat_L owns the outer face l, hat_R the inner face l - 1
                    C_num = 5e-1_wp*(nc_iface_vel_y_vf(2)%sf(k, l, q) + nc_iface_vel_y_hatR_vf(2)%sf(k, l - 1, q))/y_cc(l)

                    if (alt_soundspeed) then
                        ! Same two-component K as the HLLD anchor state (see m_riemann_solver_hypo_hlld.fpp), including the
                        ! verysmall denominator regularization
                        pres_K = q_prim_vf(eqn_idx%E)%sf(k, l, q)
                        blkmod1_K = f_bulk_modulus(pres_K, gammas(1), pi_infs(1)) + (4._wp/3._wp)*Gs_hypo(1)
                        blkmod2_K = f_bulk_modulus(pres_K, gammas(2), pi_infs(2)) + (4._wp/3._wp)*Gs_hypo(2)
                        K_K = q_prim_vf(eqn_idx%adv%beg)%sf(k, l, q)*q_prim_vf(eqn_idx%adv%end)%sf(k, l, &
                                        & q)*(blkmod2_K - blkmod1_K)/(q_prim_vf(eqn_idx%adv%beg)%sf(k, l, &
                                        & q)*blkmod2_K + q_prim_vf(eqn_idx%adv%end)%sf(k, l, q)*blkmod1_K + verysmall)
                        rhs_vf(eqn_idx%adv%beg)%sf(k, l, q) = rhs_vf(eqn_idx%adv%beg)%sf(k, l, q) + K_K*C_num
                        rhs_vf(eqn_idx%adv%end)%sf(k, l, q) = rhs_vf(eqn_idx%adv%end)%sf(k, l, q) - K_K*C_num
                    end if

                    rhs_vf(eqn_idx%stress%beg)%sf(k, l, q) = rhs_vf(eqn_idx%stress%beg)%sf(k, l, &
                           & q) - rho_K*(2._wp*q_prim_vf(eqn_idx%stress%beg)%sf(k, l, q) + 2._wp*G_K/3._wp)*C_num

                    rhs_vf(eqn_idx%stress%beg + 1)%sf(k, l, q) = rhs_vf(eqn_idx%stress%beg + 1)%sf(k, l, &
                           & q) - 2._wp*rho_K*q_prim_vf(eqn_idx%stress%beg + 1)%sf(k, l, q)*C_num

                    rhs_vf(eqn_idx%stress%beg + 2)%sf(k, l, q) = rhs_vf(eqn_idx%stress%beg + 2)%sf(k, l, &
                           & q) - rho_K*(2._wp*q_prim_vf(eqn_idx%stress%beg + 2)%sf(k, l, q) + 2._wp*G_K/3._wp)*C_num

                    rhs_vf(eqn_idx%stress%beg + 3)%sf(k, l, q) = rhs_vf(eqn_idx%stress%beg + 3)%sf(k, l, &
                           & q) + (4._wp/3._wp)*rho_K*G_K*C_num
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_hypoelastic_rhs_axisym_geom_dual_pass

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

    !> Maximum eigenvalue of the symmetric 2x2 matrix [[a, b], [b, c]]
    pure function f_max_eig_sym2x2(a, b, c) result(eig_max)

        $:GPU_ROUTINE(function_name='f_max_eig_sym2x2', parallelism='[seq]')
        real(wp), intent(in) :: a, b, c
        real(wp)             :: eig_max

        eig_max = 0.5_wp*(a + c) + sqrt((0.5_wp*(a - c))**2._wp + b*b)

    end function f_max_eig_sym2x2

    !> Maximum eigenvalue of a symmetric 3x3 matrix via the trigonometric closed form on its invariants; the acos argument is
    !! clamped and hydrostatic/repeated-eigenvalue states fall back to I1/3 to avoid 0/0
    pure function f_max_eig_sym3x3(t_xx, t_xy, t_yy, t_xz, t_yz, t_zz) result(eig_max)

        $:GPU_ROUTINE(function_name='f_max_eig_sym3x3', parallelism='[seq]')
        real(wp), intent(in) :: t_xx, t_xy, t_yy, t_xz, t_yz, t_zz
        real(wp)             :: eig_max
        real(wp)             :: I1, I2, I3, sqrt_term, argument

        I1 = t_xx + t_yy + t_zz
        I2 = t_xx*t_yy + t_xx*t_zz + t_yy*t_zz - (t_xy**2._wp + t_xz**2._wp + t_yz**2._wp)
        I3 = t_xx*t_yy*t_zz + 2._wp*t_xy*t_xz*t_yz - t_xx*t_yz**2._wp - t_yy*t_xz**2._wp - t_zz*t_xy**2._wp

        sqrt_term = sqrt(max(I1*I1 - 3._wp*I2, 0._wp))
        if (sqrt_term > verysmall) then
            argument = (2._wp*I1*I1*I1 - 9._wp*I1*I2 + 27._wp*I3)/(2._wp*sqrt_term*sqrt_term*sqrt_term)
            if (argument > 1._wp) argument = 1._wp
            if (argument < -1._wp) argument = -1._wp
            eig_max = I1/3._wp + (2._wp/3._wp)*sqrt_term*cos(acos(argument)/3._wp)
        else
            eig_max = I1/3._wp
        end if

    end function f_max_eig_sym3x3

    !> Accumulate the continuum damage source: the overstress rate on the maximum principal Cauchy stress sigma = -p I + tau_e (full
    !! 3D principal set in every dimensionality), weighted by the damageable-solid partial mass
    subroutine s_compute_damage_state(q_cons_vf, q_prim_vf, rhs_vf)

        type(scalar_field), dimension(sys_size), intent(in)    :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(in)    :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        real(wp)                                               :: sigma_p  !< maximum principal Cauchy stress
        real(wp)                                               :: pres, solid_partial_density
        real(wp)                                               :: tau_xx, tau_xy, tau_yy, tau_xz, tau_yz, tau_zz
        integer                                                :: q, l, k, i

        $:GPU_PARALLEL_LOOP(collapse=3, &
                            & private='[sigma_p, pres, solid_partial_density, tau_xx, tau_xy, tau_yy, tau_xz, tau_yz, tau_zz]')
        do q = 0, p
            do l = 0, n
                do k = 0, m
                    ! Damageable-solid partial mass
                    solid_partial_density = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        if (Gs_hypo(i) > verysmall) then
                            solid_partial_density = solid_partial_density + q_cons_vf(eqn_idx%cont%beg + i - 1)%sf(k, l, q)
                        end if
                    end do

                    if (solid_partial_density > verysmall .and. q_prim_vf(eqn_idx%damage)%sf(k, l, q) < 1._wp) then
                        pres = q_prim_vf(eqn_idx%E)%sf(k, l, q)
                        tau_xx = q_prim_vf(eqn_idx%stress%beg)%sf(k, l, q)

                        if (n == 0) then
                            ! Transverse deviatoric components are -tau_xx/2 (traceless closure)
                            sigma_p = -pres + max(tau_xx, -0.5_wp*tau_xx)
                        else if (p == 0) then
                            tau_xy = q_prim_vf(eqn_idx%stress%beg + 1)%sf(k, l, q)
                            tau_yy = q_prim_vf(eqn_idx%stress%beg + 2)%sf(k, l, q)
                            if (cyl_coord) then
                                ! Out-of-plane principal component is the stored hoop stress
                                tau_zz = q_prim_vf(eqn_idx%stress%beg + 3)%sf(k, l, q)
                            else
                                ! Out-of-plane deviatoric component from the traceless closure
                                tau_zz = -(tau_xx + tau_yy)
                            end if
                            sigma_p = -pres + max(f_max_eig_sym2x2(tau_xx, tau_xy, tau_yy), tau_zz)
                        else
                            tau_xy = q_prim_vf(eqn_idx%stress%beg + 1)%sf(k, l, q)
                            tau_yy = q_prim_vf(eqn_idx%stress%beg + 2)%sf(k, l, q)
                            tau_xz = q_prim_vf(eqn_idx%stress%beg + 3)%sf(k, l, q)
                            tau_yz = q_prim_vf(eqn_idx%stress%beg + 4)%sf(k, l, q)
                            tau_zz = q_prim_vf(eqn_idx%stress%beg + 5)%sf(k, l, q)
                            sigma_p = -pres + f_max_eig_sym3x3(tau_xx, tau_xy, tau_yy, tau_xz, tau_yz, tau_zz)
                        end if

                        if (sigma_p > tau_star) then
                            rhs_vf(eqn_idx%damage)%sf(k, l, q) = rhs_vf(eqn_idx%damage)%sf(k, l, &
                                   & q) + solid_partial_density*(alpha_bar*(sigma_p - tau_star))**cont_damage_s
                        end if
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_damage_state

    !> Project the conservative continuum-damage carrier onto 0 <= U_D <= m_s.
    subroutine s_enforce_cont_damage_bounds(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        real(stp)                                              :: solid_partial_density
        integer                                                :: q, l, k, i

        $:GPU_PARALLEL_LOOP(collapse=3, private='[solid_partial_density]')
        do q = 0, p
            do l = 0, n
                do k = 0, m
                    ! Damageable-solid partial mass
                    solid_partial_density = 0._stp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        if (Gs_hypo(i) > verysmall) then
                            solid_partial_density = solid_partial_density + q_cons_vf(eqn_idx%cont%beg + i - 1)%sf(k, l, q)
                        end if
                    end do
                    solid_partial_density = max(solid_partial_density, 0._stp)
                    q_cons_vf(eqn_idx%damage)%sf(k, l, q) = min(max(q_cons_vf(eqn_idx%damage)%sf(k, l, q), 0._stp), &
                              & solid_partial_density)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_enforce_cont_damage_bounds

end module m_hypoelastic
