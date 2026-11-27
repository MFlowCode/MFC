!>
!! @file
!! @brief Contains module m_viscous
#:include 'case.fpp'
#:include 'macros.fpp'

#! Fill Re_visc(1:2) at a cylindrical-axis cell. Non-Newtonian path uses a
#! Herschel-Bulkley mixture inv-Re from the local strain rate; Newtonian path is
#! the original alpha-weighted harmonic average. Component-3 strain terms are
#! gated for compile validity (num_dims <= 2) and runtime safety (p > 0).

#! REF(D, base, off): cell index list with the direction-D entry shifted by off, for the finite differences in
#! s_compact_ib_viscous_gradient
#:set SH = lambda v, o: v + (' %+d' % o if o else '')
#:set REF = lambda D, B, o: ', '.join(SH(b, o) if i + 1 == D else b for i, b in enumerate(B))

!> @brief Computes viscous stress tensors and diffusive flux contributions for the Navier--Stokes equations
module m_viscous

    use m_derived_types
    use m_global_parameters
    use m_weno
    use m_muscl
    use m_helper
    use m_finite_differences
    use m_variables_conversion, only: s_compute_mixture_coefficients
    use m_constants, only: model_eqns_5eq, recon_type_weno, recon_type_muscl
    use m_hb_function

    private; public s_get_viscous, s_compute_viscous_stress_cylindrical_boundary, s_compute_axisym_viscous_hoop, &
        & s_initialize_viscous_module, s_reconstruct_cell_boundary_values_visc_deriv, s_compact_ib_viscous_gradients, &
        & s_compact_normal_viscous_gradients, s_finalize_viscous_module, Res_viscous

    type(int_bounds_info) :: iv
    type(int_bounds_info) :: is1_viscous, is2_viscous, is3_viscous
    $:GPU_DECLARE(create='[is1_viscous, is2_viscous, is3_viscous, iv]')

    real(wp), allocatable, dimension(:,:) :: Res_viscous
    $:GPU_DECLARE(create='[Res_viscous]')

    !> Faces within this many cells of an immersed-boundary cell take the compact velocity gradient
    integer, parameter :: ib_visc_band = 2

contains

    !> Initialize the viscous module
    impure subroutine s_initialize_viscous_module

        integer :: i, j  !< generic loop iterators

        @:ALLOCATE(Res_viscous(1:2, 1:Re_size_max))

        do i = 1, 2
            do j = 1, Re_size(i)
                Res_viscous(i, j) = fluid_pp(Re_idx(i, j))%Re(i)
            end do
        end do
        $:GPU_UPDATE(device='[Res_viscous, Re_idx, Re_size]')
        $:GPU_ENTER_DATA(copyin='[is1_viscous, is2_viscous, is3_viscous, iv]')

    end subroutine s_initialize_viscous_module

    !> Compute viscous stress tensor near cylindrical axis, avoiding 1/r singularity at y_cb(-1)=0
    !> Inverse Reynolds numbers of the mixture at one cell. Non-Newtonian fluids need the local shear rate, which comes from the
    !! velocity gradients; Newtonian ones reduce to the volume-fraction-weighted harmonic mean.
    subroutine s_compute_axis_inv_re(grad_x_vf, grad_y_vf, grad_z_vf, alpha_visc, j, k, l, Re_visc)

        $:GPU_ROUTINE(function_name='s_compute_axis_inv_re', parallelism='[seq]', cray_inline=True)

        type(scalar_field), dimension(num_dims), intent(in) :: grad_x_vf, grad_y_vf, grad_z_vf
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(in) :: alpha_visc
        #:else
            real(wp), dimension(num_fluids), intent(in) :: alpha_visc
        #:endif
        integer, intent(in)                 :: j, k, l
        real(wp), dimension(2), intent(out) :: Re_visc
        real(wp)                            :: gamma_dot_c
        integer                             :: i, q

        if (any_non_newtonian) then
            gamma_dot_c = f_compute_shear_rate_from_components(grad_x_vf(1)%sf(j, k, l), grad_y_vf(2)%sf(j, k, l), 0._wp, &
                & 0.5_wp*(grad_y_vf(1)%sf(j, k, l) + grad_x_vf(2)%sf(j, k, l)), 0._wp, 0._wp)
            #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                if (p > 0) then
                    gamma_dot_c = f_compute_shear_rate_from_components(grad_x_vf(1)%sf(j, k, l), grad_y_vf(2)%sf(j, k, l), &
                        & grad_z_vf(3)%sf(j, k, l), 0.5_wp*(grad_y_vf(1)%sf(j, k, l) + grad_x_vf(2)%sf(j, k, l)), &
                        & 0.5_wp*(grad_z_vf(1)%sf(j, k, l) + grad_x_vf(3)%sf(j, k, l)), 0.5_wp*(grad_z_vf(2)%sf(j, k, &
                        & l) + grad_y_vf(3)%sf(j, k, l)))
                end if
            #:endif
            call s_compute_mixture_inv_re(alpha_visc, gamma_dot_c, Res_viscous, Re_visc)
        else
            $:GPU_LOOP(parallelism='[seq]')
            do i = 1, 2
                Re_visc(i) = dflt_real

                if (Re_size(i) > 0) Re_visc(i) = 0._wp
                $:GPU_LOOP(parallelism='[seq]')
                do q = 1, Re_size(i)
                    Re_visc(i) = alpha_visc(Re_idx(i, q))/Res_viscous(i, q) + Re_visc(i)
                end do

                Re_visc(i) = 1._wp/max(Re_visc(i), sgm_eps)
            end do
        end if

    end subroutine s_compute_axis_inv_re

    subroutine s_compute_viscous_stress_cylindrical_boundary(q_prim_vf, grad_x_vf, grad_y_vf, grad_z_vf, tau_Re_vf, ix, iy, iz)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), dimension(num_dims), intent(in) :: grad_x_vf, grad_y_vf, grad_z_vf
        type(scalar_field), dimension(1:sys_size), intent(inout) :: tau_Re_vf
        type(int_bounds_info), intent(in) :: ix, iy, iz
        real(wp) :: rho_visc, gamma_visc, pi_inf_visc, qv_visc, alpha_visc_sum  !< Mixture variables
        real(wp), dimension(2) :: Re_visc

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3)    :: alpha_visc, alpha_rho_visc
            real(wp), dimension(3, 3) :: tau_Re
        #:else
            real(wp), dimension(num_fluids)         :: alpha_visc, alpha_rho_visc
            real(wp), dimension(num_dims, num_dims) :: tau_Re
        #:endif

        integer :: i, j, k, l, q  !< Generic loop iterator

        is1_viscous = ix; is2_viscous = iy; is3_viscous = iz

        $:GPU_UPDATE(device='[is1_viscous, is2_viscous, is3_viscous]')

        $:GPU_PARALLEL_LOOP(collapse=3)
        do l = is3_viscous%beg, is3_viscous%end
            do k = is2_viscous%beg, is2_viscous%end
                do j = is1_viscous%beg, is1_viscous%end
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = eqn_idx%mom%beg, eqn_idx%E
                        tau_Re_vf(i)%sf(j, k, l) = 0._wp
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
            if (shear_stress) then  ! Shear stresses
                $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, q, rho_visc, gamma_visc, pi_inf_visc, qv_visc, &
                                    & alpha_visc_sum, alpha_visc, alpha_rho_visc, Re_visc, tau_Re]')
                do l = is3_viscous%beg, is3_viscous%end
                    do k = -1, 1
                        do j = is1_viscous%beg, is1_viscous%end
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                alpha_rho_visc(i) = q_prim_vf(i)%sf(j, k, l)
                                if (bubbles_euler .and. num_fluids == 1) then
                                    alpha_visc(i) = 1._wp - q_prim_vf(eqn_idx%E + i)%sf(j, k, l)
                                else
                                    alpha_visc(i) = q_prim_vf(eqn_idx%E + i)%sf(j, k, l)
                                end if
                            end do

                            if (mpp_lim) then
                                alpha_visc_sum = 0._wp
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    alpha_rho_visc(i) = max(0._wp, alpha_rho_visc(i))
                                    alpha_visc(i) = min(max(0._wp, alpha_visc(i)), 1._wp)
                                    alpha_visc_sum = alpha_visc_sum + alpha_visc(i)
                                end do
                                alpha_visc = alpha_visc/max(alpha_visc_sum, sgm_eps)
                            end if
                            call s_compute_mixture_coefficients(alpha_rho_visc, alpha_visc, rho_visc, gamma_visc, pi_inf_visc, &
                                                                & qv_visc)
                            if (viscous) then
                                call s_compute_axis_inv_re(grad_x_vf, grad_y_vf, grad_z_vf, alpha_visc, j, k, l, Re_visc)
                            end if

                            ! Shear stress near cylindrical axis: includes v/r hoop term
                            tau_Re(2, 1) = (grad_y_vf(1)%sf(j, k, l) + grad_x_vf(2)%sf(j, k, l))/Re_visc(1)

                            tau_Re(2, 2) = (4._wp*grad_y_vf(2)%sf(j, k, l) - 2._wp*grad_x_vf(1)%sf(j, k, &
                                   & l) - 2._wp*q_prim_vf(eqn_idx%mom%beg + 1)%sf(j, k, l)/y_cc(k))/(3._wp*Re_visc(1))
                            ! Viscous flux contribution to momentum and energy equations
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, 2
                                tau_Re_vf(eqn_idx%cont%end + i)%sf(j, k, l) = tau_Re_vf(eqn_idx%cont%end + i)%sf(j, k, &
                                          & l) - tau_Re(2, i)

                                tau_Re_vf(eqn_idx%E)%sf(j, k, l) = tau_Re_vf(eqn_idx%E)%sf(j, k, &
                                          & l) - q_prim_vf(eqn_idx%cont%end + i)%sf(j, k, l)*tau_Re(2, i)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endif

        #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
            if (bulk_stress) then  ! Bulk stresses
                $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, q, rho_visc, gamma_visc, pi_inf_visc, qv_visc, &
                                    & alpha_visc_sum, alpha_visc, alpha_rho_visc, Re_visc, tau_Re]')
                do l = is3_viscous%beg, is3_viscous%end
                    do k = -1, 1
                        do j = is1_viscous%beg, is1_viscous%end
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                alpha_rho_visc(i) = q_prim_vf(i)%sf(j, k, l)
                                if (bubbles_euler .and. num_fluids == 1) then
                                    alpha_visc(i) = 1._wp - q_prim_vf(eqn_idx%E + i)%sf(j, k, l)
                                else
                                    alpha_visc(i) = q_prim_vf(eqn_idx%E + i)%sf(j, k, l)
                                end if
                            end do

                            if (mpp_lim) then
                                alpha_visc_sum = 0._wp
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    alpha_rho_visc(i) = max(0._wp, alpha_rho_visc(i))
                                    alpha_visc(i) = min(max(0._wp, alpha_visc(i)), 1._wp)
                                    alpha_visc_sum = alpha_visc_sum + alpha_visc(i)
                                end do
                                alpha_visc = alpha_visc/max(alpha_visc_sum, sgm_eps)
                            end if
                            call s_compute_mixture_coefficients(alpha_rho_visc, alpha_visc, rho_visc, gamma_visc, pi_inf_visc, &
                                                                & qv_visc)
                            if (viscous) then
                                call s_compute_axis_inv_re(grad_x_vf, grad_y_vf, grad_z_vf, alpha_visc, j, k, l, Re_visc)
                            end if

                            tau_Re(2, 2) = (grad_x_vf(1)%sf(j, k, l) + grad_y_vf(2)%sf(j, k, &
                                   & l) + q_prim_vf(eqn_idx%mom%beg + 1)%sf(j, k, l)/y_cc(k))/Re_visc(2)

                            tau_Re_vf(eqn_idx%mom%beg + 1)%sf(j, k, l) = tau_Re_vf(eqn_idx%mom%beg + 1)%sf(j, k, l) - tau_Re(2, 2)

                            tau_Re_vf(eqn_idx%E)%sf(j, k, l) = tau_Re_vf(eqn_idx%E)%sf(j, k, &
                                      & l) - q_prim_vf(eqn_idx%mom%beg + 1)%sf(j, k, l)*tau_Re(2, 2)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endif

        if (p == 0) return
        #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
            if (shear_stress) then  ! Shear stresses
                $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, q, rho_visc, gamma_visc, pi_inf_visc, qv_visc, &
                                    & alpha_visc_sum, alpha_visc, alpha_rho_visc, Re_visc, tau_Re]')
                do l = is3_viscous%beg, is3_viscous%end
                    do k = -1, 1
                        do j = is1_viscous%beg, is1_viscous%end
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                alpha_rho_visc(i) = q_prim_vf(i)%sf(j, k, l)
                                if (bubbles_euler .and. num_fluids == 1) then
                                    alpha_visc(i) = 1._wp - q_prim_vf(eqn_idx%E + i)%sf(j, k, l)
                                else
                                    alpha_visc(i) = q_prim_vf(eqn_idx%E + i)%sf(j, k, l)
                                end if
                            end do

                            if (mpp_lim) then
                                alpha_visc_sum = 0._wp
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    alpha_rho_visc(i) = max(0._wp, alpha_rho_visc(i))
                                    alpha_visc(i) = min(max(0._wp, alpha_visc(i)), 1._wp)
                                    alpha_visc_sum = alpha_visc_sum + alpha_visc(i)
                                end do
                                alpha_visc = alpha_visc/max(alpha_visc_sum, sgm_eps)
                            end if
                            call s_compute_mixture_coefficients(alpha_rho_visc, alpha_visc, rho_visc, gamma_visc, pi_inf_visc, &
                                                                & qv_visc)
                            if (viscous) then
                                call s_compute_axis_inv_re(grad_x_vf, grad_y_vf, grad_z_vf, alpha_visc, j, k, l, Re_visc)
                            end if

                            tau_Re(2, 2) = -(2._wp/3._wp)*grad_z_vf(3)%sf(j, k, l)/y_cc(k)/Re_visc(1)

                            tau_Re(2, 3) = ((grad_z_vf(2)%sf(j, k, l) - q_prim_vf(eqn_idx%mom%end)%sf(j, k, &
                                   & l))/y_cc(k) + grad_y_vf(3)%sf(j, k, l))/Re_visc(1)

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 2, 3
                                tau_Re_vf(eqn_idx%cont%end + i)%sf(j, k, l) = tau_Re_vf(eqn_idx%cont%end + i)%sf(j, k, &
                                          & l) - tau_Re(2, i)

                                tau_Re_vf(eqn_idx%E)%sf(j, k, l) = tau_Re_vf(eqn_idx%E)%sf(j, k, &
                                          & l) - q_prim_vf(eqn_idx%cont%end + i)%sf(j, k, l)*tau_Re(2, i)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            if (bulk_stress) then  ! Bulk stresses
                $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, q, rho_visc, gamma_visc, pi_inf_visc, qv_visc, &
                                    & alpha_visc_sum, alpha_visc, alpha_rho_visc, Re_visc, tau_Re]')
                do l = is3_viscous%beg, is3_viscous%end
                    do k = -1, 1
                        do j = is1_viscous%beg, is1_viscous%end
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                alpha_rho_visc(i) = q_prim_vf(i)%sf(j, k, l)
                                if (bubbles_euler .and. num_fluids == 1) then
                                    alpha_visc(i) = 1._wp - q_prim_vf(eqn_idx%E + i)%sf(j, k, l)
                                else
                                    alpha_visc(i) = q_prim_vf(eqn_idx%E + i)%sf(j, k, l)
                                end if
                            end do

                            if (mpp_lim) then
                                alpha_visc_sum = 0._wp
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    alpha_rho_visc(i) = max(0._wp, alpha_rho_visc(i))
                                    alpha_visc(i) = min(max(0._wp, alpha_visc(i)), 1._wp)
                                    alpha_visc_sum = alpha_visc_sum + alpha_visc(i)
                                end do
                                alpha_visc = alpha_visc/max(alpha_visc_sum, sgm_eps)
                            end if
                            call s_compute_mixture_coefficients(alpha_rho_visc, alpha_visc, rho_visc, gamma_visc, pi_inf_visc, &
                                                                & qv_visc)
                            if (viscous) then
                                call s_compute_axis_inv_re(grad_x_vf, grad_y_vf, grad_z_vf, alpha_visc, j, k, l, Re_visc)
                            end if

                            tau_Re(2, 2) = grad_z_vf(3)%sf(j, k, l)/y_cc(k)/Re_visc(2)

                            tau_Re_vf(eqn_idx%mom%beg + 1)%sf(j, k, l) = tau_Re_vf(eqn_idx%mom%beg + 1)%sf(j, k, l) - tau_Re(2, 2)

                            tau_Re_vf(eqn_idx%E)%sf(j, k, l) = tau_Re_vf(eqn_idx%E)%sf(j, k, &
                                      & l) - q_prim_vf(eqn_idx%mom%beg + 1)%sf(j, k, l)*tau_Re(2, 2)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endif

    end subroutine s_compute_viscous_stress_cylindrical_boundary

    !> 2D axisymmetric hoop stress: adds -tau_theta_theta/r to the radial momentum at every cell, with tau_theta_theta = 2 mu u_r/r
    !! + (mu_b - 2 mu/3) div(u) from cell-centred gradients (3D takes it from theta faces)
    subroutine s_compute_axisym_viscous_hoop(q_prim_vf, grad_x_vf, grad_y_vf, grad_z_vf, rhs_vf)

        type(scalar_field), dimension(sys_size), intent(in)    :: q_prim_vf
        type(scalar_field), dimension(num_dims), intent(in)    :: grad_x_vf, grad_y_vf, grad_z_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        real(wp)                                               :: alpha_visc_sum, ur_r, div_v, tau_tt
        real(wp), dimension(2)                                 :: Re_visc

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: alpha_visc
        #:else
            real(wp), dimension(num_fluids) :: alpha_visc
        #:endif
        integer :: i, j, k, l

        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, alpha_visc_sum, alpha_visc, Re_visc, ur_r, div_v, tau_tt]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        if (bubbles_euler .and. num_fluids == 1) then
                            alpha_visc(i) = 1._wp - q_prim_vf(eqn_idx%E + i)%sf(j, k, l)
                        else
                            alpha_visc(i) = q_prim_vf(eqn_idx%E + i)%sf(j, k, l)
                        end if
                    end do
                    if (mpp_lim) then
                        alpha_visc_sum = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, num_fluids
                            alpha_visc(i) = min(max(0._wp, alpha_visc(i)), 1._wp)
                            alpha_visc_sum = alpha_visc_sum + alpha_visc(i)
                        end do
                        alpha_visc = alpha_visc/max(alpha_visc_sum, sgm_eps)
                    end if
                    call s_compute_axis_inv_re(grad_x_vf, grad_y_vf, grad_z_vf, alpha_visc, j, k, l, Re_visc)

                    ur_r = q_prim_vf(eqn_idx%mom%beg + 1)%sf(j, k, l)/y_cc(k)
                    div_v = grad_x_vf(1)%sf(j, k, l) + grad_y_vf(2)%sf(j, k, l) + ur_r
                    tau_tt = 0._wp
                    if (shear_stress) tau_tt = (2._wp*ur_r - (2._wp/3._wp)*div_v)/Re_visc(1)
                    if (bulk_stress) tau_tt = tau_tt + div_v/Re_visc(2)
                    rhs_vf(eqn_idx%mom%beg + 1)%sf(j, k, l) = rhs_vf(eqn_idx%mom%beg + 1)%sf(j, k, l) - tau_tt/y_cc(k)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_axisym_viscous_hoop

    !> Computes viscous terms
    subroutine s_get_viscous(qL_prim_rsx_vf, dqL_prim_dx_n, dqL_prim_dy_n, dqL_prim_dz_n, qL_prim, qR_prim_rsx_vf, dqR_prim_dx_n, &
                             & dqR_prim_dy_n, dqR_prim_dz_n, qR_prim, q_prim_qp, dq_prim_dx_qp, dq_prim_dy_qp, dq_prim_dz_qp, ix, &
                             & iy, iz)

        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout) :: qL_prim_rsx_vf, qR_prim_rsx_vf
        type(vector_field), dimension(num_dims), intent(inout) :: qL_prim, qR_prim
        type(vector_field), intent(in) :: q_prim_qp
        type(vector_field), dimension(1:num_dims), intent(inout) :: dqL_prim_dx_n, dqR_prim_dx_n, dqL_prim_dy_n, dqR_prim_dy_n, &
             & dqL_prim_dz_n, dqR_prim_dz_n

        type(vector_field), dimension(1), intent(inout) :: dq_prim_dx_qp, dq_prim_dy_qp, dq_prim_dz_qp
        type(int_bounds_info), intent(in)               :: ix, iy, iz
        integer                                         :: i, j, k, l

        iv%beg = eqn_idx%mom%beg; iv%end = eqn_idx%mom%end
        $:GPU_UPDATE(device='[iv]')

        if (weno_Re_flux) then
            do i = 1, num_dims
                call s_reconstruct_cell_boundary_values_visc(q_prim_qp%vf(iv%beg:iv%end), qL_prim_rsx_vf, qR_prim_rsx_vf, i, &
                    & qL_prim(i)%vf(iv%beg:iv%end), qR_prim(i)%vf(iv%beg:iv%end), ix, iy, iz)
            end do

            ! Compute velocity gradients via divergence theorem on cell-boundary reconstructed values
            do i = 1, num_dims
                if (i == 1) then
                    call s_apply_scalar_divergence_theorem(qL_prim(i)%vf(iv%beg:iv%end), qR_prim(i)%vf(iv%beg:iv%end), &
                                                           & dq_prim_dx_qp(1)%vf(iv%beg:iv%end), i, ix, iy, iz, iv, dx, m, &
                                                           & buff_size)
                else if (i == 2) then
                    call s_apply_scalar_divergence_theorem(qL_prim(i)%vf(iv%beg:iv%end), qR_prim(i)%vf(iv%beg:iv%end), &
                                                           & dq_prim_dy_qp(1)%vf(iv%beg:iv%end), i, ix, iy, iz, iv, dy, n, &
                                                           & buff_size)
                else
                    call s_apply_scalar_divergence_theorem(qL_prim(i)%vf(iv%beg:iv%end), qR_prim(i)%vf(iv%beg:iv%end), &
                                                           & dq_prim_dz_qp(1)%vf(iv%beg:iv%end), i, ix, iy, iz, iv, dz, p, &
                                                           & buff_size)
                end if
            end do
        else  ! Compute velocity gradients at cell centers using central finite differences
            iv%beg = eqn_idx%mom%beg; iv%end = eqn_idx%mom%end
            $:GPU_UPDATE(device='[iv]')

            is1_viscous = ix; is2_viscous = iy; is3_viscous = iz

            $:GPU_UPDATE(device='[is1_viscous, is2_viscous, is3_viscous]')

            $:GPU_PARALLEL_LOOP(collapse=3)
            do l = is3_viscous%beg, is3_viscous%end
                do k = iy%beg, iy%end
                    do j = is1_viscous%beg + 1, is1_viscous%end
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = iv%beg, iv%end
                            dqL_prim_dx_n(1)%vf(i)%sf(j, k, l) = (q_prim_qp%vf(i)%sf(j, k, l) - q_prim_qp%vf(i)%sf(j - 1, k, &
                                          & l))/(x_cc(j) - x_cc(j - 1))
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            $:GPU_PARALLEL_LOOP(collapse=3)
            do l = is3_viscous%beg, is3_viscous%end
                do k = is2_viscous%beg, is2_viscous%end
                    do j = is1_viscous%beg, is1_viscous%end - 1
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = iv%beg, iv%end
                            dqR_prim_dx_n(1)%vf(i)%sf(j, k, l) = (q_prim_qp%vf(i)%sf(j + 1, k, l) - q_prim_qp%vf(i)%sf(j, k, &
                                          & l))/(x_cc(j + 1) - x_cc(j))
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            if (n > 0) then
                #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do l = is3_viscous%beg, is3_viscous%end
                        do j = is2_viscous%beg + 1, is2_viscous%end
                            do k = is1_viscous%beg, is1_viscous%end
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = iv%beg, iv%end
                                    dqL_prim_dy_n(2)%vf(i)%sf(k, j, l) = (q_prim_qp%vf(i)%sf(k, j, l) - q_prim_qp%vf(i)%sf(k, &
                                                  & j - 1, l))/(y_cc(j) - y_cc(j - 1))
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do l = is3_viscous%beg, is3_viscous%end
                        do j = is2_viscous%beg, is2_viscous%end - 1
                            do k = is1_viscous%beg, is1_viscous%end
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = iv%beg, iv%end
                                    dqR_prim_dy_n(2)%vf(i)%sf(k, j, l) = (q_prim_qp%vf(i)%sf(k, j + 1, l) - q_prim_qp%vf(i)%sf(k, &
                                                  & j, l))/(y_cc(j + 1) - y_cc(j))
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do l = is3_viscous%beg, is3_viscous%end
                        do j = is2_viscous%beg + 1, is2_viscous%end
                            do k = is1_viscous%beg + 1, is1_viscous%end - 1
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = iv%beg, iv%end
                                    dqL_prim_dx_n(2)%vf(i)%sf(k, j, l) = (dqL_prim_dx_n(1)%vf(i)%sf(k, j, &
                                                  & l) + dqR_prim_dx_n(1)%vf(i)%sf(k, j, l) + dqL_prim_dx_n(1)%vf(i)%sf(k, j - 1, &
                                                  & l) + dqR_prim_dx_n(1)%vf(i)%sf(k, j - 1, l))

                                    dqL_prim_dx_n(2)%vf(i)%sf(k, j, l) = 25.e-2_wp*dqL_prim_dx_n(2)%vf(i)%sf(k, j, l)
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do l = is3_viscous%beg, is3_viscous%end
                        do j = is2_viscous%beg, is2_viscous%end - 1
                            do k = is1_viscous%beg + 1, is1_viscous%end - 1
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = iv%beg, iv%end
                                    dqR_prim_dx_n(2)%vf(i)%sf(k, j, l) = (dqL_prim_dx_n(1)%vf(i)%sf(k, j + 1, &
                                                  & l) + dqR_prim_dx_n(1)%vf(i)%sf(k, j + 1, l) + dqL_prim_dx_n(1)%vf(i)%sf(k, j, &
                                                  & l) + dqR_prim_dx_n(1)%vf(i)%sf(k, j, l))

                                    dqR_prim_dx_n(2)%vf(i)%sf(k, j, l) = 25.e-2_wp*dqR_prim_dx_n(2)%vf(i)%sf(k, j, l)
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do l = is3_viscous%beg, is3_viscous%end
                        do k = is2_viscous%beg + 1, is2_viscous%end - 1
                            do j = is1_viscous%beg + 1, is1_viscous%end
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = iv%beg, iv%end
                                    dqL_prim_dy_n(1)%vf(i)%sf(j, k, l) = (dqL_prim_dy_n(2)%vf(i)%sf(j, k, &
                                                  & l) + dqR_prim_dy_n(2)%vf(i)%sf(j, k, l) + dqL_prim_dy_n(2)%vf(i)%sf(j - 1, k, &
                                                  & l) + dqR_prim_dy_n(2)%vf(i)%sf(j - 1, k, l))

                                    dqL_prim_dy_n(1)%vf(i)%sf(j, k, l) = 25.e-2_wp*dqL_prim_dy_n(1)%vf(i)%sf(j, k, l)
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do l = is3_viscous%beg, is3_viscous%end
                        do k = is2_viscous%beg + 1, is2_viscous%end - 1
                            do j = is1_viscous%beg, is1_viscous%end - 1
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = iv%beg, iv%end
                                    dqR_prim_dy_n(1)%vf(i)%sf(j, k, l) = (dqL_prim_dy_n(2)%vf(i)%sf(j + 1, k, &
                                                  & l) + dqR_prim_dy_n(2)%vf(i)%sf(j + 1, k, l) + dqL_prim_dy_n(2)%vf(i)%sf(j, k, &
                                                  & l) + dqR_prim_dy_n(2)%vf(i)%sf(j, k, l))

                                    dqR_prim_dy_n(1)%vf(i)%sf(j, k, l) = 25.e-2_wp*dqR_prim_dy_n(1)%vf(i)%sf(j, k, l)
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                #:endif

                if (p > 0) then
                    #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                        $:GPU_PARALLEL_LOOP(collapse=3)
                        do j = is3_viscous%beg + 1, is3_viscous%end
                            do l = is2_viscous%beg, is2_viscous%end
                                do k = is1_viscous%beg, is1_viscous%end
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = iv%beg, iv%end
                                        dqL_prim_dz_n(3)%vf(i)%sf(k, l, j) = (q_prim_qp%vf(i)%sf(k, l, j) - q_prim_qp%vf(i)%sf(k, &
                                                      & l, j - 1))/(z_cc(j) - z_cc(j - 1))
                                    end do
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()

                        $:GPU_PARALLEL_LOOP(collapse=3)
                        do j = is3_viscous%beg, is3_viscous%end - 1
                            do l = is2_viscous%beg, is2_viscous%end
                                do k = is1_viscous%beg, is1_viscous%end
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = iv%beg, iv%end
                                        dqR_prim_dz_n(3)%vf(i)%sf(k, l, j) = (q_prim_qp%vf(i)%sf(k, l, &
                                                      & j + 1) - q_prim_qp%vf(i)%sf(k, l, j))/(z_cc(j + 1) - z_cc(j))
                                    end do
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()

                        $:GPU_PARALLEL_LOOP(collapse=3)
                        do l = is3_viscous%beg + 1, is3_viscous%end - 1
                            do k = is2_viscous%beg, is2_viscous%end
                                do j = is1_viscous%beg + 1, is1_viscous%end
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = iv%beg, iv%end
                                        dqL_prim_dz_n(1)%vf(i)%sf(j, k, l) = (dqL_prim_dz_n(3)%vf(i)%sf(j, k, &
                                                      & l) + dqR_prim_dz_n(3)%vf(i)%sf(j, k, &
                                                      & l) + dqL_prim_dz_n(3)%vf(i)%sf(j - 1, k, &
                                                      & l) + dqR_prim_dz_n(3)%vf(i)%sf(j - 1, k, l))

                                        dqL_prim_dz_n(1)%vf(i)%sf(j, k, l) = 25.e-2_wp*dqL_prim_dz_n(1)%vf(i)%sf(j, k, l)
                                    end do
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()

                        $:GPU_PARALLEL_LOOP(collapse=3)
                        do l = is3_viscous%beg + 1, is3_viscous%end - 1
                            do k = is2_viscous%beg, is2_viscous%end
                                do j = is1_viscous%beg, is1_viscous%end - 1
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = iv%beg, iv%end
                                        dqR_prim_dz_n(1)%vf(i)%sf(j, k, l) = (dqL_prim_dz_n(3)%vf(i)%sf(j + 1, k, &
                                                      & l) + dqR_prim_dz_n(3)%vf(i)%sf(j + 1, k, &
                                                      & l) + dqL_prim_dz_n(3)%vf(i)%sf(j, k, l) + dqR_prim_dz_n(3)%vf(i)%sf(j, k, &
                                                      & l))

                                        dqR_prim_dz_n(1)%vf(i)%sf(j, k, l) = 25.e-2_wp*dqR_prim_dz_n(1)%vf(i)%sf(j, k, l)
                                    end do
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()

                        $:GPU_PARALLEL_LOOP(collapse=3)
                        do l = is3_viscous%beg + 1, is3_viscous%end - 1
                            do j = is2_viscous%beg + 1, is2_viscous%end
                                do k = is1_viscous%beg, is1_viscous%end
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = iv%beg, iv%end
                                        dqL_prim_dz_n(2)%vf(i)%sf(k, j, l) = (dqL_prim_dz_n(3)%vf(i)%sf(k, j, &
                                                      & l) + dqR_prim_dz_n(3)%vf(i)%sf(k, j, l) + dqL_prim_dz_n(3)%vf(i)%sf(k, &
                                                      & j - 1, l) + dqR_prim_dz_n(3)%vf(i)%sf(k, j - 1, l))

                                        dqL_prim_dz_n(2)%vf(i)%sf(k, j, l) = 25.e-2_wp*dqL_prim_dz_n(2)%vf(i)%sf(k, j, l)
                                    end do
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()

                        $:GPU_PARALLEL_LOOP(collapse=3)
                        do l = is3_viscous%beg + 1, is3_viscous%end - 1
                            do j = is2_viscous%beg, is2_viscous%end - 1
                                do k = is1_viscous%beg, is1_viscous%end
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = iv%beg, iv%end
                                        dqR_prim_dz_n(2)%vf(i)%sf(k, j, l) = (dqL_prim_dz_n(3)%vf(i)%sf(k, j + 1, &
                                                      & l) + dqR_prim_dz_n(3)%vf(i)%sf(k, j + 1, &
                                                      & l) + dqL_prim_dz_n(3)%vf(i)%sf(k, j, l) + dqR_prim_dz_n(3)%vf(i)%sf(k, j, &
                                                      & l))

                                        dqR_prim_dz_n(2)%vf(i)%sf(k, j, l) = 25.e-2_wp*dqR_prim_dz_n(2)%vf(i)%sf(k, j, l)
                                    end do
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()

                        $:GPU_PARALLEL_LOOP(collapse=3)
                        do j = is3_viscous%beg + 1, is3_viscous%end
                            do l = is2_viscous%beg + 1, is2_viscous%end - 1
                                do k = is1_viscous%beg, is1_viscous%end
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = iv%beg, iv%end
                                        dqL_prim_dy_n(3)%vf(i)%sf(k, l, j) = (dqL_prim_dy_n(2)%vf(i)%sf(k, l, &
                                                      & j) + dqR_prim_dy_n(2)%vf(i)%sf(k, l, j) + dqL_prim_dy_n(2)%vf(i)%sf(k, l, &
                                                      & j - 1) + dqR_prim_dy_n(2)%vf(i)%sf(k, l, j - 1))

                                        dqL_prim_dy_n(3)%vf(i)%sf(k, l, j) = 25.e-2_wp*dqL_prim_dy_n(3)%vf(i)%sf(k, l, j)
                                    end do
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()

                        $:GPU_PARALLEL_LOOP(collapse=3)
                        do j = is3_viscous%beg, is3_viscous%end - 1
                            do l = is2_viscous%beg + 1, is2_viscous%end - 1
                                do k = is1_viscous%beg, is1_viscous%end
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = iv%beg, iv%end
                                        dqR_prim_dy_n(3)%vf(i)%sf(k, l, j) = (dqL_prim_dy_n(2)%vf(i)%sf(k, l, &
                                                      & j + 1) + dqR_prim_dy_n(2)%vf(i)%sf(k, l, &
                                                      & j + 1) + dqL_prim_dy_n(2)%vf(i)%sf(k, l, &
                                                      & j) + dqR_prim_dy_n(2)%vf(i)%sf(k, l, j))

                                        dqR_prim_dy_n(3)%vf(i)%sf(k, l, j) = 25.e-2_wp*dqR_prim_dy_n(3)%vf(i)%sf(k, l, j)
                                    end do
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()
                        $:GPU_PARALLEL_LOOP(collapse=3)
                        do j = is3_viscous%beg + 1, is3_viscous%end
                            do l = is2_viscous%beg, is2_viscous%end
                                do k = is1_viscous%beg + 1, is1_viscous%end - 1
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = iv%beg, iv%end
                                        dqL_prim_dx_n(3)%vf(i)%sf(k, l, j) = (dqL_prim_dx_n(1)%vf(i)%sf(k, l, &
                                                      & j) + dqR_prim_dx_n(1)%vf(i)%sf(k, l, j) + dqL_prim_dx_n(1)%vf(i)%sf(k, l, &
                                                      & j - 1) + dqR_prim_dx_n(1)%vf(i)%sf(k, l, j - 1))

                                        dqL_prim_dx_n(3)%vf(i)%sf(k, l, j) = 25.e-2_wp*dqL_prim_dx_n(3)%vf(i)%sf(k, l, j)
                                    end do
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()
                        $:GPU_PARALLEL_LOOP(collapse=3)
                        do j = is3_viscous%beg, is3_viscous%end - 1
                            do l = is2_viscous%beg, is2_viscous%end
                                do k = is1_viscous%beg + 1, is1_viscous%end - 1
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = iv%beg, iv%end
                                        dqR_prim_dx_n(3)%vf(i)%sf(k, l, j) = (dqL_prim_dx_n(1)%vf(i)%sf(k, l, &
                                                      & j + 1) + dqR_prim_dx_n(1)%vf(i)%sf(k, l, &
                                                      & j + 1) + dqL_prim_dx_n(1)%vf(i)%sf(k, l, &
                                                      & j) + dqR_prim_dx_n(1)%vf(i)%sf(k, l, j))

                                        dqR_prim_dx_n(3)%vf(i)%sf(k, l, j) = 25.e-2_wp*dqR_prim_dx_n(3)%vf(i)%sf(k, l, j)
                                    end do
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()

                        do i = iv%beg, iv%end
                            call s_compute_fd_gradient(q_prim_qp%vf(i), dq_prim_dx_qp(1)%vf(i), dq_prim_dy_qp(1)%vf(i), &
                                                       & dq_prim_dz_qp(1)%vf(i))
                        end do
                    #:endif
                else
                    do i = iv%beg, iv%end
                        call s_compute_fd_gradient(q_prim_qp%vf(i), dq_prim_dx_qp(1)%vf(i), dq_prim_dy_qp(1)%vf(i), &
                                                   & dq_prim_dy_qp(1)%vf(i))
                    end do
                end if
            else
                do i = iv%beg, iv%end
                    call s_compute_fd_gradient(q_prim_qp%vf(i), dq_prim_dx_qp(1)%vf(i), dq_prim_dx_qp(1)%vf(i), &
                                               & dq_prim_dx_qp(1)%vf(i))
                end do
            end if
        end if

    end subroutine s_get_viscous

    !> Reconstruct left and right cell-boundary values of viscous primitive variables
    subroutine s_reconstruct_cell_boundary_values_visc(v_vf, vL_x, vR_x, norm_dir, vL_prim_vf, vR_prim_vf, ix, iy, iz)

        type(scalar_field), dimension(iv%beg:iv%end), intent(in) :: v_vf
        type(scalar_field), dimension(iv%beg:iv%end), intent(inout) :: vL_prim_vf, vR_prim_vf
        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout) :: vL_x, vR_x
        integer, intent(in) :: norm_dir
        type(int_bounds_info), intent(in) :: ix, iy, iz
        integer :: recon_dir  !< Coordinate direction of the WENO reconstruction
        integer :: i, j, k, l

        #:for SCHEME, TYPE in [('weno','recon_type_weno'), ('muscl','recon_type_muscl')]
            if (recon_type == ${TYPE}$) then
                ! Reconstruction in s1-direction

                if (norm_dir == 1) then
                    is1_viscous = ix; is2_viscous = iy; is3_viscous = iz
                    recon_dir = 1; is1_viscous%beg = is1_viscous%beg + ${SCHEME}$_polyn
                    is1_viscous%end = is1_viscous%end - ${SCHEME}$_polyn
                else if (norm_dir == 2) then
                    is1_viscous = iy; is2_viscous = ix; is3_viscous = iz
                    recon_dir = 2; is1_viscous%beg = is1_viscous%beg + ${SCHEME}$_polyn
                    is1_viscous%end = is1_viscous%end - ${SCHEME}$_polyn
                else
                    is1_viscous = iz; is2_viscous = iy; is3_viscous = ix
                    recon_dir = 3; is1_viscous%beg = is1_viscous%beg + ${SCHEME}$_polyn
                    is1_viscous%end = is1_viscous%end - ${SCHEME}$_polyn
                end if

                $:GPU_UPDATE(device='[is1_viscous, is2_viscous, is3_viscous, iv]')

                call s_${SCHEME}$ (v_vf(iv%beg:iv%end), vL_x(:,:,:,iv%beg:iv%end), vR_x(:,:,:,iv%beg:iv%end), recon_dir, &
                                   & is1_viscous, is2_viscous, is3_viscous)
            end if
        #:endfor

        if (viscous) then
            if (weno_Re_flux) then
                if (norm_dir == 2) then
                    $:GPU_PARALLEL_LOOP(collapse=4)
                    do i = iv%beg, iv%end
                        do l = is3_viscous%beg, is3_viscous%end
                            do j = is1_viscous%beg, is1_viscous%end
                                do k = is2_viscous%beg, is2_viscous%end
                                    vL_prim_vf(i)%sf(k, j, l) = vL_x(k, j, l, i)
                                    vR_prim_vf(i)%sf(k, j, l) = vR_x(k, j, l, i)
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                else if (norm_dir == 3) then
                    $:GPU_PARALLEL_LOOP(collapse=4)
                    do i = iv%beg, iv%end
                        do j = is1_viscous%beg, is1_viscous%end
                            do k = is2_viscous%beg, is2_viscous%end
                                do l = is3_viscous%beg, is3_viscous%end
                                    vL_prim_vf(i)%sf(l, k, j) = vL_x(l, k, j, i)
                                    vR_prim_vf(i)%sf(l, k, j) = vR_x(l, k, j, i)
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                else if (norm_dir == 1) then
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
                    $:END_GPU_PARALLEL_LOOP()
                end if
            end if
        end if

    end subroutine s_reconstruct_cell_boundary_values_visc

    !> Reconstruct left and right cell-boundary values of viscous primitive variable derivatives
    subroutine s_reconstruct_cell_boundary_values_visc_deriv(v_vf, vL_x, vR_x, norm_dir, vL_prim_vf, vR_prim_vf, ix, iy, iz)

        type(scalar_field), dimension(iv%beg:iv%end), intent(in) :: v_vf
        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,iv%beg:), intent(inout) :: vL_x, vR_x
        type(scalar_field), dimension(iv%beg:iv%end), intent(inout) :: vL_prim_vf, vR_prim_vf
        type(int_bounds_info), intent(in) :: ix, iy, iz
        integer, intent(in) :: norm_dir
        integer :: recon_dir  !< Coordinate direction of the WENO reconstruction
        integer :: i, j, k, l

        #:for SCHEME, TYPE in [('weno','recon_type_weno'), ('muscl','recon_type_muscl')]
            if (recon_type == ${TYPE}$) then
                ! Reconstruction in s1-direction

                if (norm_dir == 1) then
                    is1_viscous = ix; is2_viscous = iy; is3_viscous = iz
                    recon_dir = 1; is1_viscous%beg = is1_viscous%beg + ${SCHEME}$_polyn
                    is1_viscous%end = is1_viscous%end - ${SCHEME}$_polyn
                else if (norm_dir == 2) then
                    is1_viscous = iy; is2_viscous = ix; is3_viscous = iz
                    recon_dir = 2; is1_viscous%beg = is1_viscous%beg + ${SCHEME}$_polyn
                    is1_viscous%end = is1_viscous%end - ${SCHEME}$_polyn
                else
                    is1_viscous = iz; is2_viscous = iy; is3_viscous = ix
                    recon_dir = 3; is1_viscous%beg = is1_viscous%beg + ${SCHEME}$_polyn
                    is1_viscous%end = is1_viscous%end - ${SCHEME}$_polyn
                end if
                $:GPU_UPDATE(device='[is1_viscous, is2_viscous, is3_viscous, iv]')

                call s_${SCHEME}$ (v_vf(iv%beg:iv%end), vL_x(:,:,:,iv%beg:iv%end), vR_x(:,:,:,iv%beg:iv%end), recon_dir, &
                                   & is1_viscous, is2_viscous, is3_viscous)
            end if
        #:endfor

        if (viscous) then
            if (weno_Re_flux) then
                if (norm_dir == 2) then
                    $:GPU_PARALLEL_LOOP(collapse=4)
                    do i = iv%beg, iv%end
                        do l = is3_viscous%beg, is3_viscous%end
                            do j = is1_viscous%beg, is1_viscous%end
                                do k = is2_viscous%beg, is2_viscous%end
                                    vL_prim_vf(i)%sf(k, j, l) = vL_x(k, j, l, i)
                                    vR_prim_vf(i)%sf(k, j, l) = vR_x(k, j, l, i)
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                else if (norm_dir == 3) then
                    $:GPU_PARALLEL_LOOP(collapse=4)
                    do i = iv%beg, iv%end
                        do j = is1_viscous%beg, is1_viscous%end
                            do k = is2_viscous%beg, is2_viscous%end
                                do l = is3_viscous%beg, is3_viscous%end
                                    vL_prim_vf(i)%sf(l, k, j) = vL_x(l, k, j, i)
                                    vR_prim_vf(i)%sf(l, k, j) = vR_x(l, k, j, i)
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                else if (norm_dir == 1) then
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
                    $:END_GPU_PARALLEL_LOOP()
                end if
            end if
        end if

    end subroutine s_reconstruct_cell_boundary_values_visc_deriv

    !> Compute cell-average spatial derivatives via the scalar divergence theorem
    subroutine s_apply_scalar_divergence_theorem(vL_vf, vR_vf, dv_ds_vf, norm_dir, ix, iy, iz, iv_in, dL, dim, buff_size_in)

        ! arrays of cell widths
        type(scalar_field), dimension(iv%beg:iv%end), intent(in)          :: vL_vf, vR_vf
        type(scalar_field), dimension(iv%beg:iv%end), intent(inout)       :: dv_ds_vf
        integer, intent(in)                                               :: norm_dir
        type(int_bounds_info), intent(in)                                 :: ix, iy, iz, iv_in
        integer, intent(in)                                               :: dim, buff_size_in
        real(wp), dimension(-buff_size_in:dim + buff_size_in), intent(in) :: dL
        integer                                                           :: i, j, k, l  !< Generic loop iterators

        is1_viscous = ix
        is2_viscous = iy
        is3_viscous = iz
        iv = iv_in

        $:GPU_UPDATE(device='[is1_viscous, is2_viscous, is3_viscous, iv]')

        ! First-Order Spatial Derivatives in x-direction
        if (norm_dir == 1) then
            ! A general application of the scalar divergence theorem that utilizes the left and right cell-boundary
            ! integral-averages, inside each cell, or an arithmetic mean of these two at the cell-boundaries, to calculate the
            ! cell-averaged first-order spatial derivatives inside the cell.

            $:GPU_PARALLEL_LOOP(collapse=3)
            do l = is3_viscous%beg, is3_viscous%end
                do k = is2_viscous%beg, is2_viscous%end
                    do j = is1_viscous%beg + 1, is1_viscous%end - 1
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = iv%beg, iv%end
                            dv_ds_vf(i)%sf(j, k, l) = 1._wp/((1._wp + wa_flg)*dL(j))*(wa_flg*vL_vf(i)%sf(j + 1, k, &
                                     & l) + vR_vf(i)%sf(j, k, l) - vL_vf(i)%sf(j, k, l) - wa_flg*vR_vf(i)%sf(j - 1, k, l))
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            ! END: First-Order Spatial Derivatives in x-direction

            ! First-Order Spatial Derivatives in y-direction
        else if (norm_dir == 2) then
            ! A general application of the scalar divergence theorem that utilizes the left and right cell-boundary
            ! integral-averages, inside each cell, or an arithmetic mean of these two at the cell-boundaries, to calculate the
            ! cell-averaged first-order spatial derivatives inside the cell.

            $:GPU_PARALLEL_LOOP(collapse=3)
            do l = is3_viscous%beg, is3_viscous%end
                do k = is2_viscous%beg + 1, is2_viscous%end - 1
                    do j = is1_viscous%beg, is1_viscous%end
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = iv%beg, iv%end
                            dv_ds_vf(i)%sf(j, k, l) = 1._wp/((1._wp + wa_flg)*dL(k))*(wa_flg*vL_vf(i)%sf(j, k + 1, &
                                     & l) + vR_vf(i)%sf(j, k, l) - vL_vf(i)%sf(j, k, l) - wa_flg*vR_vf(i)%sf(j, k - 1, l))
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            ! END: First-Order Spatial Derivatives in y-direction

            ! First-Order Spatial Derivatives in z-direction
        else
            ! A general application of the scalar divergence theorem that utilizes the left and right cell-boundary
            ! integral-averages, inside each cell, or an arithmetic mean of these two at the cell-boundaries, to calculate the
            ! cell-averaged first-order spatial derivatives inside the cell.

            $:GPU_PARALLEL_LOOP(collapse=3)
            do l = is3_viscous%beg + 1, is3_viscous%end - 1
                do k = is2_viscous%beg, is2_viscous%end
                    do j = is1_viscous%beg, is1_viscous%end
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = iv%beg, iv%end
                            dv_ds_vf(i)%sf(j, k, l) = 1._wp/((1._wp + wa_flg)*dL(l))*(wa_flg*vL_vf(i)%sf(j, k, &
                                     & l + 1) + vR_vf(i)%sf(j, k, l) - vL_vf(i)%sf(j, k, l) - wa_flg*vR_vf(i)%sf(j, k, l - 1))
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if
        ! END: First-Order Spatial Derivatives in z-direction

    end subroutine s_apply_scalar_divergence_theorem

    !> Computes the scalar gradient fields via finite differences
    subroutine s_compute_fd_gradient(var, grad_x, grad_y, grad_z)

        type(scalar_field), intent(in)    :: var
        type(scalar_field), intent(inout) :: grad_x
        type(scalar_field), intent(inout) :: grad_y
        type(scalar_field), intent(inout) :: grad_z
        type(int_bounds_info)             :: ix, iy, iz
        integer                           :: j, k, l  !< Generic loop iterators

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

        $:GPU_UPDATE(device='[is1_viscous, is2_viscous, is3_viscous]')

        $:GPU_PARALLEL_LOOP(collapse=3)
        do l = is3_viscous%beg, is3_viscous%end
            do k = is2_viscous%beg, is2_viscous%end
                do j = is1_viscous%beg, is1_viscous%end
                    grad_x%sf(j, k, l) = (var%sf(j + 1, k, l) - var%sf(j - 1, k, l))/(x_cc(j + 1) - x_cc(j - 1))
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        if (n > 0) then
            $:GPU_PARALLEL_LOOP(collapse=3)
            do l = is3_viscous%beg, is3_viscous%end
                do k = is2_viscous%beg, is2_viscous%end
                    do j = is1_viscous%beg, is1_viscous%end
                        grad_y%sf(j, k, l) = (var%sf(j, k + 1, l) - var%sf(j, k - 1, l))/(y_cc(k + 1) - y_cc(k - 1))
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        if (p > 0) then
            $:GPU_PARALLEL_LOOP(collapse=3)
            do l = is3_viscous%beg, is3_viscous%end
                do k = is2_viscous%beg, is2_viscous%end
                    do j = is1_viscous%beg, is1_viscous%end
                        grad_z%sf(j, k, l) = (var%sf(j, k, l + 1) - var%sf(j, k, l - 1))/(z_cc(l + 1) - z_cc(l - 1))
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        $:GPU_PARALLEL_LOOP(collapse=2)
        do l = idwbuff(3)%beg, idwbuff(3)%end
            do k = idwbuff(2)%beg, idwbuff(2)%end
                grad_x%sf(idwbuff(1)%beg, k, l) = (-3._wp*var%sf(idwbuff(1)%beg, k, l) + 4._wp*var%sf(idwbuff(1)%beg + 1, k, &
                          & l) - var%sf(idwbuff(1)%beg + 2, k, l))/(x_cc(idwbuff(1)%beg + 2) - x_cc(idwbuff(1)%beg))
                grad_x%sf(idwbuff(1)%end, k, l) = (+3._wp*var%sf(idwbuff(1)%end, k, l) - 4._wp*var%sf(idwbuff(1)%end - 1, k, &
                          & l) + var%sf(idwbuff(1)%end - 2, k, l))/(x_cc(idwbuff(1)%end) - x_cc(idwbuff(1)%end - 2))
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        if (n > 0) then
            $:GPU_PARALLEL_LOOP(collapse=2)
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do j = idwbuff(1)%beg, idwbuff(1)%end
                    grad_y%sf(j, idwbuff(2)%beg, l) = (-3._wp*var%sf(j, idwbuff(2)%beg, l) + 4._wp*var%sf(j, idwbuff(2)%beg + 1, &
                              & l) - var%sf(j, idwbuff(2)%beg + 2, l))/(y_cc(idwbuff(2)%beg + 2) - y_cc(idwbuff(2)%beg))
                    grad_y%sf(j, idwbuff(2)%end, l) = (+3._wp*var%sf(j, idwbuff(2)%end, l) - 4._wp*var%sf(j, idwbuff(2)%end - 1, &
                              & l) + var%sf(j, idwbuff(2)%end - 2, l))/(y_cc(idwbuff(2)%end) - y_cc(idwbuff(2)%end - 2))
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            if (p > 0) then
                $:GPU_PARALLEL_LOOP(collapse=2)
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        grad_z%sf(j, k, idwbuff(3)%beg) = (-3._wp*var%sf(j, k, idwbuff(3)%beg) + 4._wp*var%sf(j, k, &
                                  & idwbuff(3)%beg + 1) - var%sf(j, k, &
                                  & idwbuff(3)%beg + 2))/(z_cc(idwbuff(3)%beg + 2) - z_cc(idwbuff(3)%beg))
                        grad_z%sf(j, k, idwbuff(3)%end) = (+3._wp*var%sf(j, k, idwbuff(3)%end) - 4._wp*var%sf(j, k, &
                                  & idwbuff(3)%end - 1) + var%sf(j, k, &
                                  & idwbuff(3)%end - 2))/(z_cc(idwbuff(3)%end) - z_cc(idwbuff(3)%end - 2))
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        end if

        if (bc_x%beg <= BC_GHOST_EXTRAP) then
            $:GPU_PARALLEL_LOOP(collapse=2)
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    grad_x%sf(0, k, l) = (-3._wp*var%sf(0, k, l) + 4._wp*var%sf(1, k, l) - var%sf(2, k, l))/(x_cc(2) - x_cc(0))
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if
        if (bc_x%end <= BC_GHOST_EXTRAP) then
            $:GPU_PARALLEL_LOOP(collapse=2)
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    grad_x%sf(m, k, l) = (3._wp*var%sf(m, k, l) - 4._wp*var%sf(m - 1, k, l) + var%sf(m - 2, k, &
                              & l))/(x_cc(m) - x_cc(m - 2))
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if
        if (n > 0) then
            if (bc_y%beg <= BC_GHOST_EXTRAP .and. bc_y%beg /= BC_NULL) then
                $:GPU_PARALLEL_LOOP(collapse=2)
                do l = idwbuff(3)%beg, idwbuff(3)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        grad_y%sf(j, 0, l) = (-3._wp*var%sf(j, 0, l) + 4._wp*var%sf(j, 1, l) - var%sf(j, 2, l))/(y_cc(2) - y_cc(0))
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
            if (bc_y%end <= BC_GHOST_EXTRAP) then
                $:GPU_PARALLEL_LOOP(collapse=2)
                do l = idwbuff(3)%beg, idwbuff(3)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        grad_y%sf(j, n, l) = (3._wp*var%sf(j, n, l) - 4._wp*var%sf(j, n - 1, l) + var%sf(j, n - 2, &
                                  & l))/(y_cc(n) - y_cc(n - 2))
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
            if (p > 0) then
                if (bc_z%beg <= BC_GHOST_EXTRAP) then
                    $:GPU_PARALLEL_LOOP(collapse=2)
                    do k = idwbuff(2)%beg, idwbuff(2)%end
                        do j = idwbuff(1)%beg, idwbuff(1)%end
                            grad_z%sf(j, k, 0) = (-3._wp*var%sf(j, k, 0) + 4._wp*var%sf(j, k, 1) - var%sf(j, k, &
                                      & 2))/(z_cc(2) - z_cc(0))
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
                if (bc_z%end <= BC_GHOST_EXTRAP) then
                    $:GPU_PARALLEL_LOOP(collapse=2)
                    do k = idwbuff(2)%beg, idwbuff(2)%end
                        do j = idwbuff(1)%beg, idwbuff(1)%end
                            grad_z%sf(j, k, p) = (3._wp*var%sf(j, k, p) - 4._wp*var%sf(j, k, p - 1) + var%sf(j, k, &
                                      & p - 2))/(z_cc(p) - z_cc(p - 2))
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
            end if
        end if

    end subroutine s_compute_fd_gradient

    !> With weno_Re_flux next to an immersed body, replace the WENO velocity gradient on the faces of direction id within
    !! ib_visc_band cells of a body cell by the compact one of weno_Re_flux = F: the normal difference across the face, and the mean
    !! of the two cells' central tangential differences. Ghost values extend the flow only a layer or two deep (the mirror leaves a
    !! kink at the wall), and the gradient stencil of weno_Re_flux, about five cells wide, carries that into the fluid faces next to
    !! the wall. The face j+1/2 is read by the viscous flux from dqR(j) and dqL(j+1), so both are set. On an AMR fine frame the
    !! markers (ib_markers_fine) hold the frame's interior only; cells outside it read as fluid (bodies sit inside their blocks).
    subroutine s_compact_ib_viscous_gradients(q_prim_vf, markers, id, dqL_x, dqR_x, dqL_y, dqR_y, dqL_z, dqR_z)

        type(scalar_field), dimension(sys_size), intent(in)                           :: q_prim_vf
        type(integer_field), intent(in)                                               :: markers
        integer, intent(in)                                                           :: id
        type(scalar_field), dimension(eqn_idx%mom%beg:eqn_idx%mom%end), intent(inout) :: dqL_x, dqR_x, dqL_y, dqR_y, dqL_z, dqR_z
        integer                                                                       :: i

        do i = eqn_idx%mom%beg, eqn_idx%mom%end
            call s_compact_ib_viscous_gradient(q_prim_vf(i), markers, id, dqL_x(i), dqR_x(i), dqL_y(i), dqR_y(i), dqL_z(i), &
                                               & dqR_z(i))
        end do

    end subroutine s_compact_ib_viscous_gradients

    !> One velocity component of s_compact_ib_viscous_gradients. The differences are written out per direction so that a coordinate
    !! array that does not exist (z_cc in 2D) is never touched.
    subroutine s_compact_ib_viscous_gradient(q, markers, id, dqL_x, dqR_x, dqL_y, dqR_y, dqL_z, dqR_z)

        type(scalar_field), intent(in)    :: q
        type(integer_field), intent(in)   :: markers
        integer, intent(in)               :: id
        type(scalar_field), intent(inout) :: dqL_x, dqR_x, dqL_y, dqR_y, dqL_z, dqR_z
        integer                           :: j, k, l, a, b, c, jn, kn, ln, ja, kb, lc, wz
        integer                           :: mlo1, mlo2, mlo3, mhi1, mhi2, mhi3
        logical                           :: near
        real(wp)                          :: g, f1p, f1m, f2p, f2m

        ja = merge(1, 0, id == 1); kb = merge(1, 0, id == 2); lc = merge(1, 0, id == 3)
        wz = merge(ib_visc_band, 0, p > 0)
        ! readable marker cells: the array's bounds, or the installed fine frame's interior (its markers have no ghost layers and
        ! cells past the frame hold another frame's markers)
        mlo1 = lbound(markers%sf, 1); mlo2 = lbound(markers%sf, 2); mlo3 = lbound(markers%sf, 3)
        mhi1 = ubound(markers%sf, 1); mhi2 = ubound(markers%sf, 2); mhi3 = ubound(markers%sf, 3)
        if (amr_in_fine_advance) then
            mlo1 = 0; mlo2 = 0; mlo3 = 0; mhi1 = m; mhi2 = n; mhi3 = p
        end if

        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, a, b, c, jn, kn, ln, near, g, f1p, f1m, f2p, f2m]', copyin='[ja, kb, &
                            & lc, wz, mlo1, mlo2, mlo3, mhi1, mhi2, mhi3]')
        do l = -lc, p
            do k = -kb, n
                do j = -ja, m
                    jn = j + ja; kn = k + kb; ln = l + lc
                    near = .false.
                    do c = max(-wz, mlo3 - l), min(wz + lc*min(wz, 1), mhi3 - l)
                        do b = max(-ib_visc_band*min(n, 1), mlo2 - k), min((ib_visc_band + kb)*min(n, 1), mhi2 - k)
                            do a = max(-ib_visc_band, mlo1 - j), min(ib_visc_band + ja, mhi1 - j)
                                if (markers%sf(j + a, k + b, l + c) /= 0) near = .true.
                            end do
                        end do
                    end do
                    if (near) then
                        #:for DIR, X, CC in [(1, 'x', 'x_cc'), (2, 'y', 'y_cc'), (3, 'z', 'z_cc')]
                            #:set C1 = ['j', 'k', 'l']
                            #:set C2 = ['jn', 'kn', 'ln']
                            if (num_dims >= ${DIR}$) then
                                ! forward and backward differences along ${X}$ at the cell below the face and the cell above it
                                f1p = (q%sf(${REF(DIR, C1, 1)}$) - q%sf(${REF(DIR, C1, 0)}$))/(${CC}$(${C1[DIR-1]}$ + 1) &
                                       & - ${CC}$(${C1[DIR-1]}$))
                                f1m = (q%sf(${REF(DIR, C1, 0)}$) - q%sf(${REF(DIR, C1, -1)}$))/(${CC}$(${C1[DIR-1]}$) &
                                       & - ${CC}$(${C1[DIR-1]}$ - 1))
                                f2p = (q%sf(${REF(DIR, C2, 1)}$) - q%sf(${REF(DIR, C2, 0)}$))/(${CC}$(${C2[DIR-1]}$ + 1) &
                                       & - ${CC}$(${C2[DIR-1]}$))
                                f2m = (q%sf(${REF(DIR, C2, 0)}$) - q%sf(${REF(DIR, C2, -1)}$))/(${CC}$(${C2[DIR-1]}$) &
                                       & - ${CC}$(${C2[DIR-1]}$ - 1))
                                g = merge(f1p, 0.25_wp*(f1p + f1m + f2p + f2m), id == ${DIR}$)
                                dqR_${X}$%sf(j, k, l) = g
                                dqL_${X}$%sf(jn, kn, ln) = g
                            end if
                        #:endfor
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compact_ib_viscous_gradient

    !> weno_Re_flux_compact: set the face-normal velocity gradient (direction id) on every face to the difference across the face.
    !! The WENO chain (face values, cell gradient, reconstructed gradient, face mean) has zero damping at the checkerboard
    !! wavenumber, which leaves a no-slip wall's odd-mirror forcing undamped. The face j+1/2 is read by the viscous flux from dqR(j)
    !! and dqL(j+1), so both are set. Tangential gradients are unchanged.
    subroutine s_compact_normal_viscous_gradients(q_prim_vf, id, dqL_x, dqR_x, dqL_y, dqR_y, dqL_z, dqR_z)

        type(scalar_field), dimension(sys_size), intent(in)                           :: q_prim_vf
        integer, intent(in)                                                           :: id
        type(scalar_field), dimension(eqn_idx%mom%beg:eqn_idx%mom%end), intent(inout) :: dqL_x, dqR_x, dqL_y, dqR_y, dqL_z, dqR_z
        integer                                                                       :: i

        do i = eqn_idx%mom%beg, eqn_idx%mom%end
            if (id == 1) call s_compact_normal_viscous_gradient(q_prim_vf(i), id, x_cc, dqL_x(i), dqR_x(i))
            if (id == 2 .and. n > 0) call s_compact_normal_viscous_gradient(q_prim_vf(i), id, y_cc, dqL_y(i), dqR_y(i))
            if (id == 3 .and. p > 0) call s_compact_normal_viscous_gradient(q_prim_vf(i), id, z_cc, dqL_z(i), dqR_z(i))
        end do

    end subroutine s_compact_normal_viscous_gradients

    !> One velocity component of s_compact_normal_viscous_gradients; cc holds the cell centres along id.
    subroutine s_compact_normal_viscous_gradient(q, id, cc, dqL, dqR)

        type(scalar_field), intent(in)               :: q
        integer, intent(in)                          :: id
        real(wp), dimension(-buff_size:), intent(in) :: cc
        type(scalar_field), intent(inout)            :: dqL, dqR
        integer                                      :: j, k, l, ja, kb, lc, ic
        real(wp)                                     :: g

        ja = merge(1, 0, id == 1); kb = merge(1, 0, id == 2); lc = merge(1, 0, id == 3)

        $:GPU_PARALLEL_LOOP(collapse=3, private='[ic, g]', copyin='[ja, kb, lc]')
        do l = -lc, p
            do k = -kb, n
                do j = -ja, m
                    ic = j*ja + k*kb + l*lc
                    g = (q%sf(j + ja, k + kb, l + lc) - q%sf(j, k, l))/(cc(ic + 1) - cc(ic))
                    dqR%sf(j, k, l) = g
                    dqL%sf(j + ja, k + kb, l + lc) = g
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compact_normal_viscous_gradient

    !> Finalize the viscous module
    impure subroutine s_finalize_viscous_module()

        @:DEALLOCATE(Res_viscous)

    end subroutine s_finalize_viscous_module

end module m_viscous
