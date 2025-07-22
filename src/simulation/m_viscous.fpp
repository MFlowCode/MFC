!>
!! @file m_viscous.f90
!! @brief Contains module m_viscous
#:include 'macros.fpp'

!> @brief The module contains the subroutines used to compute viscous terms.
module m_viscous

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_weno

    use m_muscl                !< Monotonic Upstream-centered (MUSCL)
                               !! schemes for conservation laws

    use m_helper

    use m_finite_differences

    private; public s_get_viscous, &
 s_compute_viscous_stress_tensor, &
 s_initialize_viscous_module, &
 s_reconstruct_cell_boundary_values_visc_deriv, &
 s_finalize_viscous_module

    type(int_bounds_info) :: iv
    type(int_bounds_info) :: is1_viscous, is2_viscous, is3_viscous
    $:GPU_DECLARE(create='[is1_viscous,is2_viscous,is3_viscous,iv]')

    real(wp), allocatable, dimension(:, :) :: Res_viscous
    $:GPU_DECLARE(create='[Res_viscous]')

contains

    impure subroutine s_initialize_viscous_module

        integer :: i, j !< generic loop iterators

        @:ALLOCATE(Res_viscous(1:2, 1:Re_size_max))

        do i = 1, 2
            do j = 1, Re_size(i)
                Res_viscous(i, j) = fluid_pp(Re_idx(i, j))%Re(i)
            end do
        end do
        $:GPU_UPDATE(device='[Res_viscous,Re_idx,Re_size]')
        $:GPU_ENTER_DATA(copyin='[is1_viscous,is2_viscous,is3_viscous,iv]')

    end subroutine s_initialize_viscous_module

    !> The purpose of this subroutine is to compute the viscous
    !      stress tensor for the cells directly next to the axis in
    !      cylindrical coordinates. This is necessary to avoid the
    !      1/r singularity that arises at the cell boundary coinciding
    !      with the axis, i.e., y_cb(-1) = 0.
    !  @param q_prim_vf Cell-average primitive variables
    !  @param grad_x_vf Cell-average primitive variable derivatives, x-dir
    !  @param grad_y_vf Cell-average primitive variable derivatives, y-dir
    !  @param grad_z_vf Cell-average primitive variable derivatives, z-dir
    subroutine s_compute_viscous_stress_tensor(q_prim_vf, grad_x_vf, grad_y_vf, grad_z_vf, &
                                               tau_Re_vf, &
                                               ix, iy, iz)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), dimension(num_dims), intent(in) :: grad_x_vf, grad_y_vf, grad_z_vf
        type(scalar_field), dimension(1:sys_size), intent(inout) :: tau_Re_vf
        type(int_bounds_info), intent(in) :: ix, iy, iz

        real(wp) :: rho_visc, gamma_visc, pi_inf_visc, alpha_visc_sum  !< Mixture variables
        real(wp), dimension(2) :: Re_visc
        real(wp), dimension(num_fluids) :: alpha_visc, alpha_rho_visc

        real(wp), dimension(num_dims, num_dims) :: tau_Re

        integer :: i, j, k, l, q !< Generic loop iterator

        is1_viscous = ix; is2_viscous = iy; is3_viscous = iz

        $:GPU_UPDATE(device='[is1_viscous,is2_viscous,is3_viscous]')

        $:GPU_PARALLEL_LOOP(collapse=3)
        do l = is3_viscous%beg, is3_viscous%end
            do k = is2_viscous%beg, is2_viscous%end
                do j = is1_viscous%beg, is1_viscous%end
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = momxb, E_idx
                        tau_Re_vf(i)%sf(j, k, l) = 0._wp
                    end do
                end do
            end do
        end do
        if (shear_stress) then    ! Shear stresses
            $:GPU_PARALLEL_LOOP(collapse=3, private='[alpha_visc, &
                & alpha_rho_visc, Re_visc, tau_Re]')
            do l = is3_viscous%beg, is3_viscous%end
                do k = -1, 1
                    do j = is1_viscous%beg, is1_viscous%end

                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, num_fluids
                            alpha_rho_visc(i) = q_prim_vf(i)%sf(j, k, l)
                            if (bubbles_euler .and. num_fluids == 1) then
                                alpha_visc(i) = 1._wp - q_prim_vf(E_idx + i)%sf(j, k, l)
                            else
                                alpha_visc(i) = q_prim_vf(E_idx + i)%sf(j, k, l)
                            end if
                        end do

                        if (bubbles_euler) then
                            rho_visc = 0._wp
                            gamma_visc = 0._wp
                            pi_inf_visc = 0._wp

                            if (mpp_lim .and. (model_eqns == 2) .and. (num_fluids > 2)) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    rho_visc = rho_visc + alpha_rho_visc(i)
                                    gamma_visc = gamma_visc + alpha_visc(i)*gammas(i)
                                    pi_inf_visc = pi_inf_visc + alpha_visc(i)*pi_infs(i)
                                end do
                            else if ((model_eqns == 2) .and. (num_fluids > 2)) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids - 1
                                    rho_visc = rho_visc + alpha_rho_visc(i)
                                    gamma_visc = gamma_visc + alpha_visc(i)*gammas(i)
                                    pi_inf_visc = pi_inf_visc + alpha_visc(i)*pi_infs(i)
                                end do
                            else
                                rho_visc = alpha_rho_visc(1)
                                gamma_visc = gammas(1)
                                pi_inf_visc = pi_infs(1)
                            end if
                        else
                            rho_visc = 0._wp
                            gamma_visc = 0._wp
                            pi_inf_visc = 0._wp

                            alpha_visc_sum = 0._wp

                            if (mpp_lim) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    alpha_rho_visc(i) = max(0._wp, alpha_rho_visc(i))
                                    alpha_visc(i) = min(max(0._wp, alpha_visc(i)), 1._wp)
                                    alpha_visc_sum = alpha_visc_sum + alpha_visc(i)
                                end do

                                alpha_visc = alpha_visc/max(alpha_visc_sum, sgm_eps)

                            end if

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                rho_visc = rho_visc + alpha_rho_visc(i)
                                gamma_visc = gamma_visc + alpha_visc(i)*gammas(i)
                                pi_inf_visc = pi_inf_visc + alpha_visc(i)*pi_infs(i)
                            end do

                            if (viscous) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, 2
                                    Re_visc(i) = dflt_real

                                    if (Re_size(i) > 0) Re_visc(i) = 0._wp
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do q = 1, Re_size(i)
                                        Re_visc(i) = alpha_visc(Re_idx(i, q))/Res_viscous(i, q) &
                                                     + Re_visc(i)
                                    end do

                                    Re_visc(i) = 1._wp/max(Re_visc(i), sgm_eps)

                                end do
                            end if
                        end if

                        tau_Re(2, 1) = (grad_y_vf(1)%sf(j, k, l) + &
                                        grad_x_vf(2)%sf(j, k, l))/ &
                                       Re_visc(1)

                        tau_Re(2, 2) = (4._wp*grad_y_vf(2)%sf(j, k, l) &
                                        - 2._wp*grad_x_vf(1)%sf(j, k, l) &
                                        - 2._wp*q_prim_vf(momxb + 1)%sf(j, k, l)/y_cc(k))/ &
                                       (3._wp*Re_visc(1))
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, 2
                            tau_Re_vf(contxe + i)%sf(j, k, l) = &
                                tau_Re_vf(contxe + i)%sf(j, k, l) - &
                                tau_Re(2, i)

                            tau_Re_vf(E_idx)%sf(j, k, l) = &
                                tau_Re_vf(E_idx)%sf(j, k, l) - &
                                q_prim_vf(contxe + i)%sf(j, k, l)*tau_Re(2, i)
                        end do
                    end do
                end do
            end do
        end if

        if (bulk_stress) then    ! Bulk stresses
            $:GPU_PARALLEL_LOOP(collapse=3, private='[alpha_visc, &
                & alpha_rho_visc, Re_visc, tau_Re]')
            do l = is3_viscous%beg, is3_viscous%end
                do k = -1, 1
                    do j = is1_viscous%beg, is1_viscous%end

                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, num_fluids
                            alpha_rho_visc(i) = q_prim_vf(i)%sf(j, k, l)
                            if (bubbles_euler .and. num_fluids == 1) then
                                alpha_visc(i) = 1._wp - q_prim_vf(E_idx + i)%sf(j, k, l)
                            else
                                alpha_visc(i) = q_prim_vf(E_idx + i)%sf(j, k, l)
                            end if
                        end do

                        if (bubbles_euler) then
                            rho_visc = 0._wp
                            gamma_visc = 0._wp
                            pi_inf_visc = 0._wp

                            if (mpp_lim .and. (model_eqns == 2) .and. (num_fluids > 2)) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    rho_visc = rho_visc + alpha_rho_visc(i)
                                    gamma_visc = gamma_visc + alpha_visc(i)*gammas(i)
                                    pi_inf_visc = pi_inf_visc + alpha_visc(i)*pi_infs(i)
                                end do
                            else if ((model_eqns == 2) .and. (num_fluids > 2)) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids - 1
                                    rho_visc = rho_visc + alpha_rho_visc(i)
                                    gamma_visc = gamma_visc + alpha_visc(i)*gammas(i)
                                    pi_inf_visc = pi_inf_visc + alpha_visc(i)*pi_infs(i)
                                end do
                            else
                                rho_visc = alpha_rho_visc(1)
                                gamma_visc = gammas(1)
                                pi_inf_visc = pi_infs(1)
                            end if
                        else
                            rho_visc = 0._wp
                            gamma_visc = 0._wp
                            pi_inf_visc = 0._wp

                            alpha_visc_sum = 0._wp

                            if (mpp_lim) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    alpha_rho_visc(i) = max(0._wp, alpha_rho_visc(i))
                                    alpha_visc(i) = min(max(0._wp, alpha_visc(i)), 1._wp)
                                    alpha_visc_sum = alpha_visc_sum + alpha_visc(i)
                                end do

                                alpha_visc = alpha_visc/max(alpha_visc_sum, sgm_eps)

                            end if

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                rho_visc = rho_visc + alpha_rho_visc(i)
                                gamma_visc = gamma_visc + alpha_visc(i)*gammas(i)
                                pi_inf_visc = pi_inf_visc + alpha_visc(i)*pi_infs(i)
                            end do

                            if (viscous) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, 2
                                    Re_visc(i) = dflt_real

                                    if (Re_size(i) > 0) Re_visc(i) = 0._wp
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do q = 1, Re_size(i)
                                        Re_visc(i) = alpha_visc(Re_idx(i, q))/Res_viscous(i, q) &
                                                     + Re_visc(i)
                                    end do

                                    Re_visc(i) = 1._wp/max(Re_visc(i), sgm_eps)

                                end do
                            end if
                        end if

                        tau_Re(2, 2) = (grad_x_vf(1)%sf(j, k, l) + &
                                        grad_y_vf(2)%sf(j, k, l) + &
                                        q_prim_vf(momxb + 1)%sf(j, k, l)/y_cc(k))/ &
                                       Re_visc(2)

                        tau_Re_vf(momxb + 1)%sf(j, k, l) = &
                            tau_Re_vf(momxb + 1)%sf(j, k, l) - &
                            tau_Re(2, 2)

                        tau_Re_vf(E_idx)%sf(j, k, l) = &
                            tau_Re_vf(E_idx)%sf(j, k, l) - &
                            q_prim_vf(momxb + 1)%sf(j, k, l)*tau_Re(2, 2)

                    end do
                end do
            end do
        end if

        if (p == 0) return

        if (shear_stress) then    ! Shear stresses
            $:GPU_PARALLEL_LOOP(collapse=3, private='[alpha_visc, &
                & alpha_rho_visc, Re_visc, tau_Re]')
            do l = is3_viscous%beg, is3_viscous%end
                do k = -1, 1
                    do j = is1_viscous%beg, is1_viscous%end

                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, num_fluids
                            alpha_rho_visc(i) = q_prim_vf(i)%sf(j, k, l)
                            if (bubbles_euler .and. num_fluids == 1) then
                                alpha_visc(i) = 1._wp - q_prim_vf(E_idx + i)%sf(j, k, l)
                            else
                                alpha_visc(i) = q_prim_vf(E_idx + i)%sf(j, k, l)
                            end if
                        end do

                        if (bubbles_euler) then
                            rho_visc = 0._wp
                            gamma_visc = 0._wp
                            pi_inf_visc = 0._wp

                            if (mpp_lim .and. (model_eqns == 2) .and. (num_fluids > 2)) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    rho_visc = rho_visc + alpha_rho_visc(i)
                                    gamma_visc = gamma_visc + alpha_visc(i)*gammas(i)
                                    pi_inf_visc = pi_inf_visc + alpha_visc(i)*pi_infs(i)
                                end do
                            else if ((model_eqns == 2) .and. (num_fluids > 2)) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids - 1
                                    rho_visc = rho_visc + alpha_rho_visc(i)
                                    gamma_visc = gamma_visc + alpha_visc(i)*gammas(i)
                                    pi_inf_visc = pi_inf_visc + alpha_visc(i)*pi_infs(i)
                                end do
                            else
                                rho_visc = alpha_rho_visc(1)
                                gamma_visc = gammas(1)
                                pi_inf_visc = pi_infs(1)
                            end if
                        else
                            rho_visc = 0._wp
                            gamma_visc = 0._wp
                            pi_inf_visc = 0._wp

                            alpha_visc_sum = 0._wp

                            if (mpp_lim) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    alpha_rho_visc(i) = max(0._wp, alpha_rho_visc(i))
                                    alpha_visc(i) = min(max(0._wp, alpha_visc(i)), 1._wp)
                                    alpha_visc_sum = alpha_visc_sum + alpha_visc(i)
                                end do

                                alpha_visc = alpha_visc/max(alpha_visc_sum, sgm_eps)

                            end if

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                rho_visc = rho_visc + alpha_rho_visc(i)
                                gamma_visc = gamma_visc + alpha_visc(i)*gammas(i)
                                pi_inf_visc = pi_inf_visc + alpha_visc(i)*pi_infs(i)
                            end do

                            if (viscous) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, 2
                                    Re_visc(i) = dflt_real

                                    if (Re_size(i) > 0) Re_visc(i) = 0._wp
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do q = 1, Re_size(i)
                                        Re_visc(i) = alpha_visc(Re_idx(i, q))/Res_viscous(i, q) &
                                                     + Re_visc(i)
                                    end do

                                    Re_visc(i) = 1._wp/max(Re_visc(i), sgm_eps)

                                end do
                            end if
                        end if

                        tau_Re(2, 2) = -(2._wp/3._wp)*grad_z_vf(3)%sf(j, k, l)/y_cc(k)/ &
                                       Re_visc(1)

                        tau_Re(2, 3) = ((grad_z_vf(2)%sf(j, k, l) - &
                                         q_prim_vf(momxe)%sf(j, k, l))/ &
                                        y_cc(k) + grad_y_vf(3)%sf(j, k, l))/ &
                                       Re_visc(1)

                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 2, 3
                            tau_Re_vf(contxe + i)%sf(j, k, l) = &
                                tau_Re_vf(contxe + i)%sf(j, k, l) - &
                                tau_Re(2, i)

                            tau_Re_vf(E_idx)%sf(j, k, l) = &
                                tau_Re_vf(E_idx)%sf(j, k, l) - &
                                q_prim_vf(contxe + i)%sf(j, k, l)*tau_Re(2, i)
                        end do

                    end do
                end do
            end do
        end if

        if (bulk_stress) then    ! Bulk stresses
            $:GPU_PARALLEL_LOOP(collapse=3, private='[alpha_visc, &
                & alpha_rho_visc, Re_visc, tau_Re]')
            do l = is3_viscous%beg, is3_viscous%end
                do k = -1, 1
                    do j = is1_viscous%beg, is1_viscous%end

                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, num_fluids
                            alpha_rho_visc(i) = q_prim_vf(i)%sf(j, k, l)
                            if (bubbles_euler .and. num_fluids == 1) then
                                alpha_visc(i) = 1._wp - q_prim_vf(E_idx + i)%sf(j, k, l)
                            else
                                alpha_visc(i) = q_prim_vf(E_idx + i)%sf(j, k, l)
                            end if
                        end do

                        if (bubbles_euler) then
                            rho_visc = 0._wp
                            gamma_visc = 0._wp
                            pi_inf_visc = 0._wp

                            if (mpp_lim .and. (model_eqns == 2) .and. (num_fluids > 2)) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    rho_visc = rho_visc + alpha_rho_visc(i)
                                    gamma_visc = gamma_visc + alpha_visc(i)*gammas(i)
                                    pi_inf_visc = pi_inf_visc + alpha_visc(i)*pi_infs(i)
                                end do
                            else if ((model_eqns == 2) .and. (num_fluids > 2)) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids - 1
                                    rho_visc = rho_visc + alpha_rho_visc(i)
                                    gamma_visc = gamma_visc + alpha_visc(i)*gammas(i)
                                    pi_inf_visc = pi_inf_visc + alpha_visc(i)*pi_infs(i)
                                end do
                            else
                                rho_visc = alpha_rho_visc(1)
                                gamma_visc = gammas(1)
                                pi_inf_visc = pi_infs(1)
                            end if
                        else
                            rho_visc = 0._wp
                            gamma_visc = 0._wp
                            pi_inf_visc = 0._wp

                            alpha_visc_sum = 0._wp

                            if (mpp_lim) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    alpha_rho_visc(i) = max(0._wp, alpha_rho_visc(i))
                                    alpha_visc(i) = min(max(0._wp, alpha_visc(i)), 1._wp)
                                    alpha_visc_sum = alpha_visc_sum + alpha_visc(i)
                                end do

                                alpha_visc = alpha_visc/max(alpha_visc_sum, sgm_eps)

                            end if

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                rho_visc = rho_visc + alpha_rho_visc(i)
                                gamma_visc = gamma_visc + alpha_visc(i)*gammas(i)
                                pi_inf_visc = pi_inf_visc + alpha_visc(i)*pi_infs(i)
                            end do

                            if (viscous) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, 2
                                    Re_visc(i) = dflt_real

                                    if (Re_size(i) > 0) Re_visc(i) = 0._wp
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do q = 1, Re_size(i)
                                        Re_visc(i) = alpha_visc(Re_idx(i, q))/Res_viscous(i, q) &
                                                     + Re_visc(i)
                                    end do

                                    Re_visc(i) = 1._wp/max(Re_visc(i), sgm_eps)

                                end do
                            end if
                        end if

                        tau_Re(2, 2) = grad_z_vf(3)%sf(j, k, l)/y_cc(k)/ &
                                       Re_visc(2)

                        tau_Re_vf(momxb + 1)%sf(j, k, l) = &
                            tau_Re_vf(momxb + 1)%sf(j, k, l) - &
                            tau_Re(2, 2)

                        tau_Re_vf(E_idx)%sf(j, k, l) = &
                            tau_Re_vf(E_idx)%sf(j, k, l) - &
                            q_prim_vf(momxb + 1)%sf(j, k, l)*tau_Re(2, 2)

                    end do
                end do
            end do
        end if
    end subroutine s_compute_viscous_stress_tensor

    !>  Computes viscous terms
    !!  @param q_cons_vf Cell-averaged conservative variables
    !!  @param q_prim_vf Cell-averaged primitive variables
    !!  @param rhs_vf Cell-averaged RHS variables
    subroutine s_get_viscous(qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, &
                             dqL_prim_dx_n, dqL_prim_dy_n, dqL_prim_dz_n, &
                             qL_prim, &
                             qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, &
                             dqR_prim_dx_n, dqR_prim_dy_n, dqR_prim_dz_n, &
                             qR_prim, &
                             q_prim_qp, &
                             dq_prim_dx_qp, dq_prim_dy_qp, dq_prim_dz_qp, &
                             ix, iy, iz)

        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:), &
            intent(inout) :: qL_prim_rsx_vf, qR_prim_rsx_vf, &
                             qL_prim_rsy_vf, qR_prim_rsy_vf, &
                             qL_prim_rsz_vf, qR_prim_rsz_vf

        type(vector_field), dimension(num_dims), intent(inout) :: qL_prim, qR_prim

        type(vector_field), intent(in) :: q_prim_qp

        type(vector_field), dimension(1:num_dims), &
            intent(inout) :: dqL_prim_dx_n, dqR_prim_dx_n, &
                             dqL_prim_dy_n, dqR_prim_dy_n, &
                             dqL_prim_dz_n, dqR_prim_dz_n

        type(vector_field), dimension(1), intent(inout) :: dq_prim_dx_qp, dq_prim_dy_qp, dq_prim_dz_qp
        type(int_bounds_info), intent(in) :: ix, iy, iz

        integer :: i, j, k, l

        do i = 1, num_dims

            iv%beg = mom_idx%beg; iv%end = mom_idx%end

            $:GPU_UPDATE(device='[iv]')

            call s_reconstruct_cell_boundary_values_visc( &
                q_prim_qp%vf(iv%beg:iv%end), &
                qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, &
                qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, &
                i, qL_prim(i)%vf(iv%beg:iv%end), qR_prim(i)%vf(iv%beg:iv%end), &
                ix, iy, iz)
        end do

        if (weno_Re_flux) then
            ! Compute velocity gradient at cell centers using scalar
            ! divergence theorem
            do i = 1, num_dims
                if (i == 1) then
                    call s_apply_scalar_divergence_theorem( &
                        qL_prim(i)%vf(iv%beg:iv%end), &
                        qR_prim(i)%vf(iv%beg:iv%end), &
                        dq_prim_dx_qp(1)%vf(iv%beg:iv%end), i, &
                        ix, iy, iz, iv, dx, m, buff_size)
                elseif (i == 2) then
                    call s_apply_scalar_divergence_theorem( &
                        qL_prim(i)%vf(iv%beg:iv%end), &
                        qR_prim(i)%vf(iv%beg:iv%end), &
                        dq_prim_dy_qp(1)%vf(iv%beg:iv%end), i, &
                        ix, iy, iz, iv, dy, n, buff_size)
                else
                    call s_apply_scalar_divergence_theorem( &
                        qL_prim(i)%vf(iv%beg:iv%end), &
                        qR_prim(i)%vf(iv%beg:iv%end), &
                        dq_prim_dz_qp(1)%vf(iv%beg:iv%end), i, &
                        ix, iy, iz, iv, dz, p, buff_size)
                end if
            end do

        else ! Compute velocity gradient at cell centers using finite differences

            iv%beg = mom_idx%beg; iv%end = mom_idx%end
            $:GPU_UPDATE(device='[iv]')

            is1_viscous = ix; is2_viscous = iy; is3_viscous = iz

            $:GPU_UPDATE(device='[is1_viscous,is2_viscous,is3_viscous]')

            $:GPU_PARALLEL_LOOP(collapse=3)
            do l = is3_viscous%beg, is3_viscous%end
                do k = iy%beg, iy%end
                    do j = is1_viscous%beg + 1, is1_viscous%end
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = iv%beg, iv%end
                            dqL_prim_dx_n(1)%vf(i)%sf(j, k, l) = &
                                (q_prim_qp%vf(i)%sf(j, k, l) - &
                                 q_prim_qp%vf(i)%sf(j - 1, k, l))/ &
                                (x_cc(j) - x_cc(j - 1))
                        end do
                    end do
                end do
            end do

            $:GPU_PARALLEL_LOOP(collapse=3)
            do l = is3_viscous%beg, is3_viscous%end
                do k = is2_viscous%beg, is2_viscous%end
                    do j = is1_viscous%beg, is1_viscous%end - 1
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = iv%beg, iv%end
                            dqR_prim_dx_n(1)%vf(i)%sf(j, k, l) = &
                                (q_prim_qp%vf(i)%sf(j + 1, k, l) - &
                                 q_prim_qp%vf(i)%sf(j, k, l))/ &
                                (x_cc(j + 1) - x_cc(j))
                        end do
                    end do
                end do
            end do

            if (n > 0) then

                $:GPU_PARALLEL_LOOP(collapse=3)
                do l = is3_viscous%beg, is3_viscous%end
                    do j = is2_viscous%beg + 1, is2_viscous%end
                        do k = is1_viscous%beg, is1_viscous%end
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = iv%beg, iv%end
                                dqL_prim_dy_n(2)%vf(i)%sf(k, j, l) = &
                                    (q_prim_qp%vf(i)%sf(k, j, l) - &
                                     q_prim_qp%vf(i)%sf(k, j - 1, l))/ &
                                    (y_cc(j) - y_cc(j - 1))
                            end do
                        end do
                    end do
                end do

                $:GPU_PARALLEL_LOOP(collapse=3)
                do l = is3_viscous%beg, is3_viscous%end
                    do j = is2_viscous%beg, is2_viscous%end - 1
                        do k = is1_viscous%beg, is1_viscous%end
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = iv%beg, iv%end
                                dqR_prim_dy_n(2)%vf(i)%sf(k, j, l) = &
                                    (q_prim_qp%vf(i)%sf(k, j + 1, l) - &
                                     q_prim_qp%vf(i)%sf(k, j, l))/ &
                                    (y_cc(j + 1) - y_cc(j))
                            end do
                        end do
                    end do
                end do

                $:GPU_PARALLEL_LOOP(collapse=3)
                do l = is3_viscous%beg, is3_viscous%end
                    do j = is2_viscous%beg + 1, is2_viscous%end
                        do k = is1_viscous%beg + 1, is1_viscous%end - 1
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = iv%beg, iv%end
                                dqL_prim_dx_n(2)%vf(i)%sf(k, j, l) = &
                                    (dqL_prim_dx_n(1)%vf(i)%sf(k, j, l) + &
                                     dqR_prim_dx_n(1)%vf(i)%sf(k, j, l) + &
                                     dqL_prim_dx_n(1)%vf(i)%sf(k, j - 1, l) + &
                                     dqR_prim_dx_n(1)%vf(i)%sf(k, j - 1, l))

                                dqL_prim_dx_n(2)%vf(i)%sf(k, j, l) = 25.e-2_wp* &
                                                                     dqL_prim_dx_n(2)%vf(i)%sf(k, j, l)
                            end do
                        end do
                    end do
                end do

                $:GPU_PARALLEL_LOOP(collapse=3)
                do l = is3_viscous%beg, is3_viscous%end
                    do j = is2_viscous%beg, is2_viscous%end - 1
                        do k = is1_viscous%beg + 1, is1_viscous%end - 1
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = iv%beg, iv%end
                                dqR_prim_dx_n(2)%vf(i)%sf(k, j, l) = &
                                    (dqL_prim_dx_n(1)%vf(i)%sf(k, j + 1, l) + &
                                     dqR_prim_dx_n(1)%vf(i)%sf(k, j + 1, l) + &
                                     dqL_prim_dx_n(1)%vf(i)%sf(k, j, l) + &
                                     dqR_prim_dx_n(1)%vf(i)%sf(k, j, l))

                                dqR_prim_dx_n(2)%vf(i)%sf(k, j, l) = 25.e-2_wp* &
                                                                     dqR_prim_dx_n(2)%vf(i)%sf(k, j, l)

                            end do
                        end do
                    end do
                end do

                $:GPU_PARALLEL_LOOP(collapse=3)
                do l = is3_viscous%beg, is3_viscous%end
                    do k = is2_viscous%beg + 1, is2_viscous%end - 1
                        do j = is1_viscous%beg + 1, is1_viscous%end
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = iv%beg, iv%end
                                dqL_prim_dy_n(1)%vf(i)%sf(j, k, l) = &
                                    (dqL_prim_dy_n(2)%vf(i)%sf(j, k, l) + &
                                     dqR_prim_dy_n(2)%vf(i)%sf(j, k, l) + &
                                     dqL_prim_dy_n(2)%vf(i)%sf(j - 1, k, l) + &
                                     dqR_prim_dy_n(2)%vf(i)%sf(j - 1, k, l))

                                dqL_prim_dy_n(1)%vf(i)%sf(j, k, l) = 25.e-2_wp* &
                                                                     dqL_prim_dy_n(1)%vf(i)%sf(j, k, l)

                            end do
                        end do
                    end do
                end do

                $:GPU_PARALLEL_LOOP(collapse=3)
                do l = is3_viscous%beg, is3_viscous%end
                    do k = is2_viscous%beg + 1, is2_viscous%end - 1
                        do j = is1_viscous%beg, is1_viscous%end - 1
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = iv%beg, iv%end
                                dqR_prim_dy_n(1)%vf(i)%sf(j, k, l) = &
                                    (dqL_prim_dy_n(2)%vf(i)%sf(j + 1, k, l) + &
                                     dqR_prim_dy_n(2)%vf(i)%sf(j + 1, k, l) + &
                                     dqL_prim_dy_n(2)%vf(i)%sf(j, k, l) + &
                                     dqR_prim_dy_n(2)%vf(i)%sf(j, k, l))

                                dqR_prim_dy_n(1)%vf(i)%sf(j, k, l) = 25.e-2_wp* &
                                                                     dqR_prim_dy_n(1)%vf(i)%sf(j, k, l)

                            end do
                        end do
                    end do
                end do

                if (p > 0) then

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do j = is3_viscous%beg + 1, is3_viscous%end
                        do l = is2_viscous%beg, is2_viscous%end
                            do k = is1_viscous%beg, is1_viscous%end
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = iv%beg, iv%end

                                    dqL_prim_dz_n(3)%vf(i)%sf(k, l, j) = &
                                        (q_prim_qp%vf(i)%sf(k, l, j) - &
                                         q_prim_qp%vf(i)%sf(k, l, j - 1))/ &
                                        (z_cc(j) - z_cc(j - 1))
                                end do
                            end do
                        end do
                    end do

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do j = is3_viscous%beg, is3_viscous%end - 1
                        do l = is2_viscous%beg, is2_viscous%end
                            do k = is1_viscous%beg, is1_viscous%end
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = iv%beg, iv%end

                                    dqR_prim_dz_n(3)%vf(i)%sf(k, l, j) = &
                                        (q_prim_qp%vf(i)%sf(k, l, j + 1) - &
                                         q_prim_qp%vf(i)%sf(k, l, j))/ &
                                        (z_cc(j + 1) - z_cc(j))
                                end do
                            end do
                        end do
                    end do

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do l = is3_viscous%beg + 1, is3_viscous%end - 1
                        do k = is2_viscous%beg, is2_viscous%end
                            do j = is1_viscous%beg + 1, is1_viscous%end
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = iv%beg, iv%end

                                    dqL_prim_dz_n(1)%vf(i)%sf(j, k, l) = &
                                        (dqL_prim_dz_n(3)%vf(i)%sf(j, k, l) + &
                                         dqR_prim_dz_n(3)%vf(i)%sf(j, k, l) + &
                                         dqL_prim_dz_n(3)%vf(i)%sf(j - 1, k, l) + &
                                         dqR_prim_dz_n(3)%vf(i)%sf(j - 1, k, l))

                                    dqL_prim_dz_n(1)%vf(i)%sf(j, k, l) = 25.e-2_wp* &
                                                                         dqL_prim_dz_n(1)%vf(i)%sf(j, k, l)

                                end do
                            end do
                        end do
                    end do

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do l = is3_viscous%beg + 1, is3_viscous%end - 1
                        do k = is2_viscous%beg, is2_viscous%end
                            do j = is1_viscous%beg, is1_viscous%end - 1
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = iv%beg, iv%end

                                    dqR_prim_dz_n(1)%vf(i)%sf(j, k, l) = &
                                        (dqL_prim_dz_n(3)%vf(i)%sf(j + 1, k, l) + &
                                         dqR_prim_dz_n(3)%vf(i)%sf(j + 1, k, l) + &
                                         dqL_prim_dz_n(3)%vf(i)%sf(j, k, l) + &
                                         dqR_prim_dz_n(3)%vf(i)%sf(j, k, l))

                                    dqR_prim_dz_n(1)%vf(i)%sf(j, k, l) = 25.e-2_wp* &
                                                                         dqR_prim_dz_n(1)%vf(i)%sf(j, k, l)

                                end do
                            end do
                        end do
                    end do

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do l = is3_viscous%beg + 1, is3_viscous%end - 1
                        do j = is2_viscous%beg + 1, is2_viscous%end
                            do k = is1_viscous%beg, is1_viscous%end
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = iv%beg, iv%end

                                    dqL_prim_dz_n(2)%vf(i)%sf(k, j, l) = &
                                        (dqL_prim_dz_n(3)%vf(i)%sf(k, j, l) + &
                                         dqR_prim_dz_n(3)%vf(i)%sf(k, j, l) + &
                                         dqL_prim_dz_n(3)%vf(i)%sf(k, j - 1, l) + &
                                         dqR_prim_dz_n(3)%vf(i)%sf(k, j - 1, l))

                                    dqL_prim_dz_n(2)%vf(i)%sf(k, j, l) = 25.e-2_wp* &
                                                                         dqL_prim_dz_n(2)%vf(i)%sf(k, j, l)

                                end do
                            end do
                        end do
                    end do

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do l = is3_viscous%beg + 1, is3_viscous%end - 1
                        do j = is2_viscous%beg, is2_viscous%end - 1
                            do k = is1_viscous%beg, is1_viscous%end
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = iv%beg, iv%end

                                    dqR_prim_dz_n(2)%vf(i)%sf(k, j, l) = &
                                        (dqL_prim_dz_n(3)%vf(i)%sf(k, j + 1, l) + &
                                         dqR_prim_dz_n(3)%vf(i)%sf(k, j + 1, l) + &
                                         dqL_prim_dz_n(3)%vf(i)%sf(k, j, l) + &
                                         dqR_prim_dz_n(3)%vf(i)%sf(k, j, l))

                                    dqR_prim_dz_n(2)%vf(i)%sf(k, j, l) = 25.e-2_wp* &
                                                                         dqR_prim_dz_n(2)%vf(i)%sf(k, j, l)

                                end do
                            end do
                        end do
                    end do

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do j = is3_viscous%beg + 1, is3_viscous%end
                        do l = is2_viscous%beg + 1, is2_viscous%end - 1
                            do k = is1_viscous%beg, is1_viscous%end
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = iv%beg, iv%end

                                    dqL_prim_dy_n(3)%vf(i)%sf(k, l, j) = &
                                        (dqL_prim_dy_n(2)%vf(i)%sf(k, l, j) + &
                                         dqR_prim_dy_n(2)%vf(i)%sf(k, l, j) + &
                                         dqL_prim_dy_n(2)%vf(i)%sf(k, l, j - 1) + &
                                         dqR_prim_dy_n(2)%vf(i)%sf(k, l, j - 1))

                                    dqL_prim_dy_n(3)%vf(i)%sf(k, l, j) = 25.e-2_wp* &
                                                                         dqL_prim_dy_n(3)%vf(i)%sf(k, l, j)

                                end do
                            end do
                        end do
                    end do

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do j = is3_viscous%beg, is3_viscous%end - 1
                        do l = is2_viscous%beg + 1, is2_viscous%end - 1
                            do k = is1_viscous%beg, is1_viscous%end
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = iv%beg, iv%end

                                    dqR_prim_dy_n(3)%vf(i)%sf(k, l, j) = &
                                        (dqL_prim_dy_n(2)%vf(i)%sf(k, l, j + 1) + &
                                         dqR_prim_dy_n(2)%vf(i)%sf(k, l, j + 1) + &
                                         dqL_prim_dy_n(2)%vf(i)%sf(k, l, j) + &
                                         dqR_prim_dy_n(2)%vf(i)%sf(k, l, j))

                                    dqR_prim_dy_n(3)%vf(i)%sf(k, l, j) = 25.e-2_wp* &
                                                                         dqR_prim_dy_n(3)%vf(i)%sf(k, l, j)

                                end do
                            end do
                        end do
                    end do
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do j = is3_viscous%beg + 1, is3_viscous%end
                        do l = is2_viscous%beg, is2_viscous%end
                            do k = is1_viscous%beg + 1, is1_viscous%end - 1
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = iv%beg, iv%end

                                    dqL_prim_dx_n(3)%vf(i)%sf(k, l, j) = &
                                        (dqL_prim_dx_n(1)%vf(i)%sf(k, l, j) + &
                                         dqR_prim_dx_n(1)%vf(i)%sf(k, l, j) + &
                                         dqL_prim_dx_n(1)%vf(i)%sf(k, l, j - 1) + &
                                         dqR_prim_dx_n(1)%vf(i)%sf(k, l, j - 1))

                                    dqL_prim_dx_n(3)%vf(i)%sf(k, l, j) = 25.e-2_wp* &
                                                                         dqL_prim_dx_n(3)%vf(i)%sf(k, l, j)

                                end do
                            end do
                        end do
                    end do
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do j = is3_viscous%beg, is3_viscous%end - 1
                        do l = is2_viscous%beg, is2_viscous%end
                            do k = is1_viscous%beg + 1, is1_viscous%end - 1
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = iv%beg, iv%end
                                    dqR_prim_dx_n(3)%vf(i)%sf(k, l, j) = &
                                        (dqL_prim_dx_n(1)%vf(i)%sf(k, l, j + 1) + &
                                         dqR_prim_dx_n(1)%vf(i)%sf(k, l, j + 1) + &
                                         dqL_prim_dx_n(1)%vf(i)%sf(k, l, j) + &
                                         dqR_prim_dx_n(1)%vf(i)%sf(k, l, j))

                                    dqR_prim_dx_n(3)%vf(i)%sf(k, l, j) = 25.e-2_wp* &
                                                                         dqR_prim_dx_n(3)%vf(i)%sf(k, l, j)

                                end do
                            end do
                        end do
                    end do

                    do i = iv%beg, iv%end
                        call s_compute_fd_gradient(q_prim_qp%vf(i), &
                                                   dq_prim_dx_qp(1)%vf(i), &
                                                   dq_prim_dy_qp(1)%vf(i), &
                                                   dq_prim_dz_qp(1)%vf(i))
                    end do

                else

                    do i = iv%beg, iv%end
                        call s_compute_fd_gradient(q_prim_qp%vf(i), &
                                                   dq_prim_dx_qp(1)%vf(i), &
                                                   dq_prim_dy_qp(1)%vf(i), &
                                                   dq_prim_dy_qp(1)%vf(i))
                    end do

                end if

            else

                do i = iv%beg, iv%end
                    call s_compute_fd_gradient(q_prim_qp%vf(i), &
                                               dq_prim_dx_qp(1)%vf(i), &
                                               dq_prim_dx_qp(1)%vf(i), &
                                               dq_prim_dx_qp(1)%vf(i))
                end do

            end if

        end if

    end subroutine s_get_viscous

    subroutine s_reconstruct_cell_boundary_values_visc(v_vf, vL_x, vL_y, vL_z, vR_x, vR_y, vR_z, &
                                                       norm_dir, vL_prim_vf, vR_prim_vf, ix, iy, iz)

        type(scalar_field), dimension(iv%beg:iv%end), intent(in) :: v_vf
        type(scalar_field), dimension(iv%beg:iv%end), intent(inout) :: vL_prim_vf, vR_prim_vf

        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:), intent(inout) :: vL_x, vL_y, vL_z, vR_x, vR_y, vR_z
        integer, intent(in) :: norm_dir
        type(int_bounds_info), intent(in) :: ix, iy, iz

        integer :: recon_dir !< Coordinate direction of the WENO reconstruction

        integer :: i, j, k, l

        #:for SCHEME, TYPE in [('weno','WENO_TYPE'), ('muscl','MUSCL_TYPE')]
            if (recon_type == ${TYPE}$) then
                ! Reconstruction in s1-direction

                if (norm_dir == 1) then
                    is1_viscous = ix; is2_viscous = iy; is3_viscous = iz
                    recon_dir = 1; is1_viscous%beg = is1_viscous%beg + ${SCHEME}$_polyn
                    is1_viscous%end = is1_viscous%end - ${SCHEME}$_polyn

                elseif (norm_dir == 2) then
                    is1_viscous = iy; is2_viscous = ix; is3_viscous = iz
                    recon_dir = 2; is1_viscous%beg = is1_viscous%beg + ${SCHEME}$_polyn
                    is1_viscous%end = is1_viscous%end - ${SCHEME}$_polyn

                else
                    is1_viscous = iz; is2_viscous = iy; is3_viscous = ix
                    recon_dir = 3; is1_viscous%beg = is1_viscous%beg + ${SCHEME}$_polyn
                    is1_viscous%end = is1_viscous%end - ${SCHEME}$_polyn

                end if

                $:GPU_UPDATE(device='[is1_viscous, is2_viscous, is3_viscous, iv]')
                if (n > 0) then
                    if (p > 0) then
                        call s_${SCHEME}$ (v_vf(iv%beg:iv%end), &
                                           vL_x(:, :, :, iv%beg:iv%end), vL_y(:, :, :, iv%beg:iv%end), vL_z(:, :, :, iv%beg:iv%end), vR_x(:, :, :, iv%beg:iv%end), vR_y(:, :, :, iv%beg:iv%end), vR_z(:, :, :, iv%beg:iv%end), &
                                           recon_dir, &
                                           is1_viscous, is2_viscous, is3_viscous)
                    else
                        call s_${SCHEME}$ (v_vf(iv%beg:iv%end), &
                                           vL_x(:, :, :, iv%beg:iv%end), vL_y(:, :, :, iv%beg:iv%end), vL_z(:, :, :, :), vR_x(:, :, :, iv%beg:iv%end), vR_y(:, :, :, iv%beg:iv%end), vR_z(:, :, :, :), &
                                           recon_dir, &
                                           is1_viscous, is2_viscous, is3_viscous)
                    end if
                else
                    call s_${SCHEME}$ (v_vf(iv%beg:iv%end), &
                                       vL_x(:, :, :, iv%beg:iv%end), vL_y(:, :, :, :), vL_z(:, :, :, :), vR_x(:, :, :, iv%beg:iv%end), vR_y(:, :, :, :), vR_z(:, :, :, :), &
                                       recon_dir, &
                                       is1_viscous, is2_viscous, is3_viscous)
                end if

                if (viscous) then
                    if (weno_Re_flux) then
                        if (norm_dir == 2) then
                            $:GPU_PARALLEL_LOOP(collapse=4)
                            do i = iv%beg, iv%end
                                do l = is3_viscous%beg, is3_viscous%end
                                    do j = is1_viscous%beg, is1_viscous%end
                                        do k = is2_viscous%beg, is2_viscous%end
                                            vL_prim_vf(i)%sf(k, j, l) = vL_y(j, k, l, i)
                                            vR_prim_vf(i)%sf(k, j, l) = vR_y(j, k, l, i)
                                        end do
                                    end do
                                end do
                            end do
                        elseif (norm_dir == 3) then
                            $:GPU_PARALLEL_LOOP(collapse=4)
                            do i = iv%beg, iv%end
                                do j = is1_viscous%beg, is1_viscous%end
                                    do k = is2_viscous%beg, is2_viscous%end
                                        do l = is3_viscous%beg, is3_viscous%end
                                            vL_prim_vf(i)%sf(l, k, j) = vL_z(j, k, l, i)
                                            vR_prim_vf(i)%sf(l, k, j) = vR_z(j, k, l, i)
                                        end do
                                    end do
                                end do
                            end do
                        elseif (norm_dir == 1) then
                            $:GPU_PARALLEL_LOOP(collapse=4)
                            do i = iv%beg, iv%end
                                do l = is3_viscous%beg, is3_viscous%end
                                    do k = is2_viscous%beg, is2_viscous%end
                                        do j = is1_viscous%beg, is1_viscous%end
                                            vL_prim_vf(i)%sf(j, k, l) = vL_x(j, k, l, i)
                                            vR_prim_vf(i)%sf(j, k, l) = vR_x(j, k, l, i)
                                        end do
                                    end do
                                end do
                            end do
                        end if
                    end if
                end if
            end if
        #:endfor
    end subroutine s_reconstruct_cell_boundary_values_visc

    subroutine s_reconstruct_cell_boundary_values_visc_deriv(v_vf, vL_x, vL_y, vL_z, vR_x, vR_y, vR_z, &
                                                             norm_dir, vL_prim_vf, vR_prim_vf, ix, iy, iz)
        type(scalar_field), dimension(iv%beg:iv%end), intent(in) :: v_vf
        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, iv%beg:), intent(inout) :: vL_x, vL_y, vL_z, vR_x, vR_y, vR_z
        type(scalar_field), dimension(iv%beg:iv%end), intent(inout) :: vL_prim_vf, vR_prim_vf
        type(int_bounds_info), intent(in) :: ix, iy, iz

        integer, intent(IN) :: norm_dir

        integer :: recon_dir !< Coordinate direction of the WENO reconstruction

        integer :: i, j, k, l

        #:for SCHEME, TYPE in [('weno','WENO_TYPE'), ('muscl','MUSCL_TYPE')]
            if (recon_type == ${TYPE}$) then
                ! Reconstruction in s1-direction

                if (norm_dir == 1) then
                    is1_viscous = ix; is2_viscous = iy; is3_viscous = iz
                    recon_dir = 1; is1_viscous%beg = is1_viscous%beg + ${SCHEME}$_polyn
                    is1_viscous%end = is1_viscous%end - ${SCHEME}$_polyn

                elseif (norm_dir == 2) then
                    is1_viscous = iy; is2_viscous = ix; is3_viscous = iz
                    recon_dir = 2; is1_viscous%beg = is1_viscous%beg + ${SCHEME}$_polyn
                    is1_viscous%end = is1_viscous%end - ${SCHEME}$_polyn

                else
                    is1_viscous = iz; is2_viscous = iy; is3_viscous = ix
                    recon_dir = 3; is1_viscous%beg = is1_viscous%beg + ${SCHEME}$_polyn
                    is1_viscous%end = is1_viscous%end - ${SCHEME}$_polyn

                end if
                $:GPU_UPDATE(device='[is1_viscous, is2_viscous, is3_viscous, iv]')
                if (n > 0) then
                    if (p > 0) then

                        call s_${SCHEME}$ (v_vf(iv%beg:iv%end), &
                                           vL_x(:, :, :, iv%beg:iv%end), vL_y(:, :, :, iv%beg:iv%end), vL_z(:, :, :, iv%beg:iv%end), vR_x(:, :, :, iv%beg:iv%end), vR_y(:, :, :, iv%beg:iv%end), vR_z(:, :, :, iv%beg:iv%end), &
                                           recon_dir, &
                                           is1_viscous, is2_viscous, is3_viscous)
                    else
                        call s_${SCHEME}$ (v_vf(iv%beg:iv%end), &
                                           vL_x(:, :, :, iv%beg:iv%end), vL_y(:, :, :, iv%beg:iv%end), vL_z(:, :, :, :), vR_x(:, :, :, iv%beg:iv%end), vR_y(:, :, :, iv%beg:iv%end), vR_z(:, :, :, :), &
                                           recon_dir, &
                                           is1_viscous, is2_viscous, is3_viscous)
                    end if
                else

                    call s_${SCHEME}$ (v_vf(iv%beg:iv%end), &
                                       vL_x(:, :, :, iv%beg:iv%end), vL_y(:, :, :, :), vL_z(:, :, :, :), vR_x(:, :, :, iv%beg:iv%end), vR_y(:, :, :, :), vR_z(:, :, :, :), &
                                       recon_dir, &
                                       is1_viscous, is2_viscous, is3_viscous)
                end if

                if (viscous) then
                    if (weno_Re_flux) then
                        if (norm_dir == 2) then
                            $:GPU_PARALLEL_LOOP(collapse=4)
                            do i = iv%beg, iv%end
                                do l = is3_viscous%beg, is3_viscous%end
                                    do j = is1_viscous%beg, is1_viscous%end
                                        do k = is2_viscous%beg, is2_viscous%end
                                            vL_prim_vf(i)%sf(k, j, l) = vL_y(j, k, l, i)
                                            vR_prim_vf(i)%sf(k, j, l) = vR_y(j, k, l, i)
                                        end do
                                    end do
                                end do
                            end do
                        elseif (norm_dir == 3) then
                            $:GPU_PARALLEL_LOOP(collapse=4)
                            do i = iv%beg, iv%end
                                do j = is1_viscous%beg, is1_viscous%end
                                    do k = is2_viscous%beg, is2_viscous%end
                                        do l = is3_viscous%beg, is3_viscous%end
                                            vL_prim_vf(i)%sf(l, k, j) = vL_z(j, k, l, i)
                                            vR_prim_vf(i)%sf(l, k, j) = vR_z(j, k, l, i)
                                        end do
                                    end do
                                end do
                            end do
                        elseif (norm_dir == 1) then
                            $:GPU_PARALLEL_LOOP(collapse=4)
                            do i = iv%beg, iv%end
                                do l = is3_viscous%beg, is3_viscous%end
                                    do k = is2_viscous%beg, is2_viscous%end
                                        do j = is1_viscous%beg, is1_viscous%end
                                            vL_prim_vf(i)%sf(j, k, l) = vL_x(j, k, l, i)
                                            vR_prim_vf(i)%sf(j, k, l) = vR_x(j, k, l, i)
                                        end do
                                    end do
                                end do
                            end do
                        end if
                    end if
                end if
            end if

        #:endfor
    end subroutine s_reconstruct_cell_boundary_values_visc_deriv

    !>  The purpose of this subroutine is to employ the inputted
        !!      left and right cell-boundary integral-averaged variables
        !!      to compute the relevant cell-average first-order spatial
        !!      derivatives in the x-, y- or z-direction by means of the
        !!      scalar divergence theorem.
        !!  @param vL_vf Left cell-boundary integral averages
        !!  @param vR_vf Right cell-boundary integral averages
        !!  @param dv_ds_vf Cell-average first-order spatial derivatives
        !!  @param norm_dir Splitting coordinate direction
    subroutine s_apply_scalar_divergence_theorem(vL_vf, vR_vf, &
                                                 dv_ds_vf, &
                                                 norm_dir, &
                                                 ix, iy, iz, iv_in, &
                                                 dL, dim, buff_size_in)

        ! arrays of cell widths
        type(scalar_field), &
            dimension(iv%beg:iv%end), &
            intent(in) :: vL_vf, vR_vf

        type(scalar_field), &
            dimension(iv%beg:iv%end), &
            intent(inout) :: dv_ds_vf

        integer, intent(in) :: norm_dir
        type(int_bounds_info), intent(in) :: ix, iy, iz, iv_in
        integer, intent(in) :: dim, buff_size_in
        real(wp), dimension(-buff_size_in:dim + buff_size_in), intent(in) :: dL

        integer :: i, j, k, l !< Generic loop iterators

        is1_viscous = ix
        is2_viscous = iy
        is3_viscous = iz
        iv = iv_in

        $:GPU_UPDATE(device='[is1_viscous, is2_viscous, is3_viscous, iv]')

        ! First-Order Spatial Derivatives in x-direction
        if (norm_dir == 1) then

            ! A general application of the scalar divergence theorem that
            ! utilizes the left and right cell-boundary integral-averages,
            ! inside each cell, or an arithmetic mean of these two at the
            ! cell-boundaries, to calculate the cell-averaged first-order
            ! spatial derivatives inside the cell.

            $:GPU_PARALLEL_LOOP(collapse=3)
            do l = is3_viscous%beg, is3_viscous%end
                do k = is2_viscous%beg, is2_viscous%end
                    do j = is1_viscous%beg + 1, is1_viscous%end - 1
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = iv%beg, iv%end
                            dv_ds_vf(i)%sf(j, k, l) = &
                                1._wp/((1._wp + wa_flg)*dL(j)) &
                                *(wa_flg*vL_vf(i)%sf(j + 1, k, l) &
                                  + vR_vf(i)%sf(j, k, l) &
                                  - vL_vf(i)%sf(j, k, l) &
                                  - wa_flg*vR_vf(i)%sf(j - 1, k, l))
                        end do
                    end do
                end do
            end do

            ! END: First-Order Spatial Derivatives in x-direction

            ! First-Order Spatial Derivatives in y-direction
        elseif (norm_dir == 2) then

            ! A general application of the scalar divergence theorem that
            ! utilizes the left and right cell-boundary integral-averages,
            ! inside each cell, or an arithmetic mean of these two at the
            ! cell-boundaries, to calculate the cell-averaged first-order
            ! spatial derivatives inside the cell.

            $:GPU_PARALLEL_LOOP(collapse=3)
            do l = is3_viscous%beg, is3_viscous%end
                do k = is2_viscous%beg + 1, is2_viscous%end - 1
                    do j = is1_viscous%beg, is1_viscous%end
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = iv%beg, iv%end
                            dv_ds_vf(i)%sf(j, k, l) = &
                                1._wp/((1._wp + wa_flg)*dL(k)) &
                                *(wa_flg*vL_vf(i)%sf(j, k + 1, l) &
                                  + vR_vf(i)%sf(j, k, l) &
                                  - vL_vf(i)%sf(j, k, l) &
                                  - wa_flg*vR_vf(i)%sf(j, k - 1, l))
                        end do
                    end do
                end do
            end do

            ! END: First-Order Spatial Derivatives in y-direction

            ! First-Order Spatial Derivatives in z-direction
        else

            ! A general application of the scalar divergence theorem that
            ! utilizes the left and right cell-boundary integral-averages,
            ! inside each cell, or an arithmetic mean of these two at the
            ! cell-boundaries, to calculate the cell-averaged first-order
            ! spatial derivatives inside the cell.

            $:GPU_PARALLEL_LOOP(collapse=3)
            do l = is3_viscous%beg + 1, is3_viscous%end - 1
                do k = is2_viscous%beg, is2_viscous%end
                    do j = is1_viscous%beg, is1_viscous%end
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = iv%beg, iv%end
                            dv_ds_vf(i)%sf(j, k, l) = &
                                1._wp/((1._wp + wa_flg)*dL(l)) &
                                *(wa_flg*vL_vf(i)%sf(j, k, l + 1) &
                                  + vR_vf(i)%sf(j, k, l) &
                                  - vL_vf(i)%sf(j, k, l) &
                                  - wa_flg*vR_vf(i)%sf(j, k, l - 1))
                        end do
                    end do
                end do
            end do

        end if
        ! END: First-Order Spatial Derivatives in z-direction

    end subroutine s_apply_scalar_divergence_theorem

    !>  Computes the scalar gradient fields via finite differences
        !!  @param var Variable to compute derivative of
        !!  @param grad_x First coordinate direction component of the derivative
        !!  @param grad_y Second coordinate direction component of the derivative
        !!  @param grad_z Third coordinate direction component of the derivative
        !!  @param norm Norm of the gradient vector
    subroutine s_compute_fd_gradient(var, grad_x, grad_y, grad_z)

        type(scalar_field), intent(in) :: var
        type(scalar_field), intent(inout) :: grad_x
        type(scalar_field), intent(inout) :: grad_y
        type(scalar_field), intent(inout) :: grad_z
        type(int_bounds_info) :: ix, iy, iz

        integer :: j, k, l !< Generic loop iterators

        ix%beg = 1 - buff_size; ix%end = m + buff_size - 1
        if (n > 0) then
            iy%beg = 1 - buff_size; iy%end = n + buff_size - 1
        else
            iy%beg = 0; iy%end = 0
        end if

        if (p > 0) then
            iz%beg = 1 - buff_size; iz%end = p + buff_size - 1
        else
            iz%beg = 0; iz%end = 0
        end if

        is1_viscous = ix; is2_viscous = iy; is3_viscous = iz

        $:GPU_UPDATE(device='[is1_viscous,is2_viscous,is3_viscous]')

        $:GPU_PARALLEL_LOOP(collapse=3)
        do l = is3_viscous%beg, is3_viscous%end
            do k = is2_viscous%beg, is2_viscous%end
                do j = is1_viscous%beg, is1_viscous%end
                    grad_x%sf(j, k, l) = &
                        (var%sf(j + 1, k, l) - var%sf(j - 1, k, l))/ &
                        (x_cc(j + 1) - x_cc(j - 1))
                end do
            end do
        end do

        if (n > 0) then
            $:GPU_PARALLEL_LOOP(collapse=3)
            do l = is3_viscous%beg, is3_viscous%end
                do k = is2_viscous%beg, is2_viscous%end
                    do j = is1_viscous%beg, is1_viscous%end
                        grad_y%sf(j, k, l) = &
                            (var%sf(j, k + 1, l) - var%sf(j, k - 1, l))/ &
                            (y_cc(k + 1) - y_cc(k - 1))
                    end do
                end do
            end do
        end if

        if (p > 0) then
            $:GPU_PARALLEL_LOOP(collapse=3)
            do l = is3_viscous%beg, is3_viscous%end
                do k = is2_viscous%beg, is2_viscous%end
                    do j = is1_viscous%beg, is1_viscous%end
                        grad_z%sf(j, k, l) = &
                            (var%sf(j, k, l + 1) - var%sf(j, k, l - 1))/ &
                            (z_cc(l + 1) - z_cc(l - 1))
                    end do
                end do
            end do
        end if

        $:GPU_PARALLEL_LOOP(collapse=2)
        do l = idwbuff(3)%beg, idwbuff(3)%end
            do k = idwbuff(2)%beg, idwbuff(2)%end
                grad_x%sf(idwbuff(1)%beg, k, l) = &
                    (-3._wp*var%sf(idwbuff(1)%beg, k, l) + 4._wp*var%sf(idwbuff(1)%beg + 1, k, l) - var%sf(idwbuff(1)%beg + 2, k, l))/ &
                    (x_cc(idwbuff(1)%beg + 2) - x_cc(idwbuff(1)%beg))
                grad_x%sf(idwbuff(1)%end, k, l) = &
                    (+3._wp*var%sf(idwbuff(1)%end, k, l) - 4._wp*var%sf(idwbuff(1)%end - 1, k, l) + var%sf(idwbuff(1)%end - 2, k, l))/ &
                    (x_cc(idwbuff(1)%end) - x_cc(idwbuff(1)%end - 2))
            end do
        end do
        if (n > 0) then
            $:GPU_PARALLEL_LOOP(collapse=2)
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do j = idwbuff(1)%beg, idwbuff(1)%end
                    grad_y%sf(j, idwbuff(2)%beg, l) = &
                        (-3._wp*var%sf(j, idwbuff(2)%beg, l) + 4._wp*var%sf(j, idwbuff(2)%beg + 1, l) - var%sf(j, idwbuff(2)%beg + 2, l))/ &
                        (y_cc(idwbuff(2)%beg + 2) - y_cc(idwbuff(2)%beg))
                    grad_y%sf(j, idwbuff(2)%end, l) = &
                        (+3._wp*var%sf(j, idwbuff(2)%end, l) - 4._wp*var%sf(j, idwbuff(2)%end - 1, l) + var%sf(j, idwbuff(2)%end - 2, l))/ &
                        (y_cc(idwbuff(2)%end) - y_cc(idwbuff(2)%end - 2))
                end do
            end do
            if (p > 0) then
                $:GPU_PARALLEL_LOOP(collapse=2)
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        grad_z%sf(j, k, idwbuff(3)%beg) = &
                            (-3._wp*var%sf(j, k, idwbuff(3)%beg) + 4._wp*var%sf(j, k, idwbuff(3)%beg + 1) - var%sf(j, k, idwbuff(3)%beg + 2))/ &
                            (z_cc(idwbuff(3)%beg + 2) - z_cc(is3_viscous%beg))
                        grad_z%sf(j, k, idwbuff(3)%end) = &
                            (+3._wp*var%sf(j, k, idwbuff(3)%end) - 4._wp*var%sf(j, k, idwbuff(3)%end - 1) + var%sf(j, k, idwbuff(3)%end - 2))/ &
                            (z_cc(idwbuff(3)%end) - z_cc(idwbuff(3)%end - 2))
                    end do
                end do
            end if
        end if

        if (bc_x%beg <= BC_GHOST_EXTRAP) then
            $:GPU_PARALLEL_LOOP(collapse=2)
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    grad_x%sf(0, k, l) = (-3._wp*var%sf(0, k, l) + 4._wp*var%sf(1, k, l) - var%sf(2, k, l))/ &
                                         (x_cc(2) - x_cc(0))
                end do
            end do
        end if
        if (bc_x%end <= BC_GHOST_EXTRAP) then
            $:GPU_PARALLEL_LOOP(collapse=2)
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    grad_x%sf(m, k, l) = (3._wp*var%sf(m, k, l) - 4._wp*var%sf(m - 1, k, l) + var%sf(m - 2, k, l))/ &
                                         (x_cc(m) - x_cc(m - 2))
                end do
            end do
        end if
        if (n > 0) then
            if (bc_y%beg <= BC_GHOST_EXTRAP .and. bc_y%beg /= BC_NULL) then
                $:GPU_PARALLEL_LOOP(collapse=2)
                do l = idwbuff(3)%beg, idwbuff(3)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        grad_y%sf(j, 0, l) = (-3._wp*var%sf(j, 0, l) + 4._wp*var%sf(j, 1, l) - var%sf(j, 2, l))/ &
                                             (y_cc(2) - y_cc(0))
                    end do
                end do
            end if
            if (bc_y%end <= BC_GHOST_EXTRAP) then
                $:GPU_PARALLEL_LOOP(collapse=2)
                do l = idwbuff(3)%beg, idwbuff(3)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        grad_y%sf(j, n, l) = (3._wp*var%sf(j, n, l) - 4._wp*var%sf(j, n - 1, l) + var%sf(j, n - 2, l))/ &
                                             (y_cc(n) - y_cc(n - 2))
                    end do
                end do
            end if
            if (p > 0) then
                if (bc_z%beg <= BC_GHOST_EXTRAP) then
                    $:GPU_PARALLEL_LOOP(collapse=2)
                    do k = idwbuff(2)%beg, idwbuff(2)%end
                        do j = idwbuff(1)%beg, idwbuff(1)%end
                            grad_z%sf(j, k, 0) = &
                                (-3._wp*var%sf(j, k, 0) + 4._wp*var%sf(j, k, 1) - var%sf(j, k, 2))/ &
                                (z_cc(2) - z_cc(0))
                        end do
                    end do
                end if
                if (bc_z%end <= BC_GHOST_EXTRAP) then
                    $:GPU_PARALLEL_LOOP(collapse=2)
                    do k = idwbuff(2)%beg, idwbuff(2)%end
                        do j = idwbuff(1)%beg, idwbuff(1)%end
                            grad_z%sf(j, k, p) = &
                                (3._wp*var%sf(j, k, p) - 4._wp*var%sf(j, k, p - 1) + var%sf(j, k, p - 2))/ &
                                (z_cc(p) - z_cc(p - 2))
                        end do
                    end do
                end if
            end if
        end if

    end subroutine s_compute_fd_gradient

    impure subroutine s_finalize_viscous_module()

        @:DEALLOCATE(Res_viscous)

    end subroutine s_finalize_viscous_module

end module m_viscous
