!>
!! @file m_re_visc.f90
!! @brief Contains module m_re_visc

#:include 'macros.fpp'

!> @brief The module contains routines that compute viscosity-related
!!        quantities for both Newtonian and non-Newtonian fluids.
!!        s_compute_re_visc returns Re_visc = 1/mu per phase.
!!        s_compute_mixture_re computes mixture Re from per-phase values.
module m_re_visc

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_hb_function          !< Herschel-Bulkley viscosity model

    implicit none

    private; public :: s_compute_re_visc, &
 s_compute_mixture_re

contains

    !> Computes velocity gradients at a single cell using finite differences.
    !!      Uses 2nd order central differences in interior, 1st order at boundaries.
    !! @param q_prim_vf Primitive variables
    !! @param j x index
    !! @param k y index
    !! @param l z index
    !! @param D_xx Output: du/dx
    !! @param D_yy Output: dv/dy
    !! @param D_zz Output: dw/dz
    !! @param D_xy Output: 0.5*(du/dy + dv/dx)
    !! @param D_xz Output: 0.5*(du/dz + dw/dx)
    !! @param D_yz Output: 0.5*(dv/dz + dw/dy)
    pure subroutine s_compute_velocity_gradients_at_cell( &
        q_prim_vf, j, k, l, D_xx, D_yy, D_zz, D_xy, D_xz, D_yz)
        $:GPU_ROUTINE(parallelism='[seq]')

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        integer, intent(in) :: j, k, l
        real(wp), intent(out) :: D_xx, D_yy, D_zz, D_xy, D_xz, D_yz

        integer :: j_lo, j_hi, k_lo, k_hi, l_lo, l_hi

        j_lo = idwbuff(1)%beg
        j_hi = idwbuff(1)%end
        k_lo = idwbuff(2)%beg
        k_hi = idwbuff(2)%end
        l_lo = idwbuff(3)%beg
        l_hi = idwbuff(3)%end

        ! Check bounds
        if (.not. ((j >= j_lo) .and. (j <= j_hi) .and. &
                   (k >= k_lo) .and. (k <= k_hi) .and. &
                   (l >= l_lo) .and. (l <= l_hi))) then
            D_xx = 0._wp; D_yy = 0._wp; D_zz = 0._wp
            D_xy = 0._wp; D_xz = 0._wp; D_yz = 0._wp
            return
        end if

        ! D_xx = du/dx
        if (j - 1 >= j_lo .and. j + 1 <= j_hi) then
            D_xx = (q_prim_vf(momxb)%sf(j + 1, k, l) - &
                    q_prim_vf(momxb)%sf(j - 1, k, l))/ &
                   (x_cc(j + 1) - x_cc(j - 1))
        else if (j + 1 <= j_hi) then
            D_xx = (q_prim_vf(momxb)%sf(j + 1, k, l) - &
                    q_prim_vf(momxb)%sf(j, k, l))/ &
                   (x_cc(j + 1) - x_cc(j))
        else if (j - 1 >= j_lo) then
            D_xx = (q_prim_vf(momxb)%sf(j, k, l) - &
                    q_prim_vf(momxb)%sf(j - 1, k, l))/ &
                   (x_cc(j) - x_cc(j - 1))
        else
            D_xx = 0._wp
        end if

        ! D_yy = dv/dy (2D and 3D only)
        if (n > 0) then
            if (k - 1 >= k_lo .and. k + 1 <= k_hi) then
                D_yy = (q_prim_vf(momxb + 1)%sf(j, k + 1, l) - &
                        q_prim_vf(momxb + 1)%sf(j, k - 1, l))/ &
                       (y_cc(k + 1) - y_cc(k - 1))
            else if (k + 1 <= k_hi) then
                D_yy = (q_prim_vf(momxb + 1)%sf(j, k + 1, l) - &
                        q_prim_vf(momxb + 1)%sf(j, k, l))/ &
                       (y_cc(k + 1) - y_cc(k))
            else if (k - 1 >= k_lo) then
                D_yy = (q_prim_vf(momxb + 1)%sf(j, k, l) - &
                        q_prim_vf(momxb + 1)%sf(j, k - 1, l))/ &
                       (y_cc(k) - y_cc(k - 1))
            else
                D_yy = 0._wp
            end if
        else
            D_yy = 0._wp
        end if

        ! D_zz = dw/dz (3D only)
        if (p > 0) then
            if (l - 1 >= l_lo .and. l + 1 <= l_hi) then
                D_zz = (q_prim_vf(momxb + 2)%sf(j, k, l + 1) - &
                        q_prim_vf(momxb + 2)%sf(j, k, l - 1))/ &
                       (z_cc(l + 1) - z_cc(l - 1))
            else if (l + 1 <= l_hi) then
                D_zz = (q_prim_vf(momxb + 2)%sf(j, k, l + 1) - &
                        q_prim_vf(momxb + 2)%sf(j, k, l))/ &
                       (z_cc(l + 1) - z_cc(l))
            else if (l - 1 >= l_lo) then
                D_zz = (q_prim_vf(momxb + 2)%sf(j, k, l) - &
                        q_prim_vf(momxb + 2)%sf(j, k, l - 1))/ &
                       (z_cc(l) - z_cc(l - 1))
            else
                D_zz = 0._wp
            end if
        else
            D_zz = 0._wp
        end if

        ! D_xy = 0.5*(du/dy + dv/dx) (2D and 3D only)
        if (n > 0) then
            if (j - 1 >= j_lo .and. j + 1 <= j_hi .and. &
                k - 1 >= k_lo .and. k + 1 <= k_hi) then
                D_xy = 0.5_wp*( &
                       (q_prim_vf(momxb)%sf(j, k + 1, l) - &
                        q_prim_vf(momxb)%sf(j, k - 1, l))/ &
                       (y_cc(k + 1) - y_cc(k - 1)) + &
                       (q_prim_vf(momxb + 1)%sf(j + 1, k, l) - &
                        q_prim_vf(momxb + 1)%sf(j - 1, k, l))/ &
                       (x_cc(j + 1) - x_cc(j - 1)))
            else
                D_xy = 0._wp
                if (k - 1 >= k_lo .and. k + 1 <= k_hi) then
                    D_xy = 0.5_wp*(q_prim_vf(momxb)%sf(j, k + 1, l) - &
                                   q_prim_vf(momxb)%sf(j, k - 1, l))/ &
                           (y_cc(k + 1) - y_cc(k - 1))
                end if
                if (j - 1 >= j_lo .and. j + 1 <= j_hi) then
                    D_xy = D_xy + 0.5_wp* &
                           (q_prim_vf(momxb + 1)%sf(j + 1, k, l) - &
                            q_prim_vf(momxb + 1)%sf(j - 1, k, l))/ &
                           (x_cc(j + 1) - x_cc(j - 1))
                end if
            end if
        else
            D_xy = 0._wp
        end if

        ! D_xz = 0.5*(du/dz + dw/dx) (3D only)
        if (p > 0) then
            if (j - 1 >= j_lo .and. j + 1 <= j_hi .and. &
                l - 1 >= l_lo .and. l + 1 <= l_hi) then
                D_xz = 0.5_wp*( &
                       (q_prim_vf(momxb)%sf(j, k, l + 1) - &
                        q_prim_vf(momxb)%sf(j, k, l - 1))/ &
                       (z_cc(l + 1) - z_cc(l - 1)) + &
                       (q_prim_vf(momxb + 2)%sf(j + 1, k, l) - &
                        q_prim_vf(momxb + 2)%sf(j - 1, k, l))/ &
                       (x_cc(j + 1) - x_cc(j - 1)))
            else
                D_xz = 0._wp
            end if
        else
            D_xz = 0._wp
        end if

        ! D_yz = 0.5*(dv/dz + dw/dy) (3D only)
        if (p > 0 .and. n > 0) then
            if (k - 1 >= k_lo .and. k + 1 <= k_hi .and. &
                l - 1 >= l_lo .and. l + 1 <= l_hi) then
                D_yz = 0.5_wp*( &
                       (q_prim_vf(momxb + 1)%sf(j, k, l + 1) - &
                        q_prim_vf(momxb + 1)%sf(j, k, l - 1))/ &
                       (z_cc(l + 1) - z_cc(l - 1)) + &
                       (q_prim_vf(momxb + 2)%sf(j, k + 1, l) - &
                        q_prim_vf(momxb + 2)%sf(j, k - 1, l))/ &
                       (y_cc(k + 1) - y_cc(k - 1)))
            else
                D_yz = 0._wp
            end if
        else
            D_yz = 0._wp
        end if

    end subroutine s_compute_velocity_gradients_at_cell

    !> Computes Re_visc per-phase for both Newtonian and non-Newtonian fluids.
    !!      Re_visc = 1/mu for each phase in each direction (shear and bulk).
    !! @param q_prim_vf Primitive variables
    !! @param alpha_visc Volume fractions
    !! @param j x index
    !! @param k y index
    !! @param l z index
    !! @param Re_visc_per_phase Output: 1/mu per fluid per direction
    !! @param grad_x_vf Optional pre-computed x-direction gradients
    !! @param grad_y_vf Optional pre-computed y-direction gradients
    !! @param grad_z_vf Optional pre-computed z-direction gradients
    pure subroutine s_compute_re_visc(q_prim_vf, alpha_visc, j, k, l, &
                                      Re_visc_per_phase, grad_x_vf, grad_y_vf, grad_z_vf)
        $:GPU_ROUTINE(parallelism='[seq]')

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(in) :: alpha_visc
        #:else
            real(wp), dimension(num_fluids), intent(in) :: alpha_visc
        #:endif
        integer, intent(in) :: j, k, l
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3, 2), intent(out) :: Re_visc_per_phase
        #:else
            real(wp), dimension(num_fluids, 2), intent(out) :: Re_visc_per_phase
        #:endif
        type(scalar_field), dimension(:), intent(in), optional :: &
            grad_x_vf, grad_y_vf, grad_z_vf

        real(wp) :: D_xx, D_yy, D_zz, D_xy, D_xz, D_yz
        real(wp) :: shear_rate, mu_fluid
        integer :: i, q, r

        ! Initialize all to default (inviscid)
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, 2
            $:GPU_LOOP(parallelism='[seq]')
            do q = 1, num_fluids
                Re_visc_per_phase(q, i) = dflt_real
            end do
        end do

        if (any_non_newtonian) then
            ! Non-Newtonian path: compute velocity gradients and shear rate
            if (present(grad_x_vf) .and. present(grad_y_vf) .and. &
                present(grad_z_vf)) then
                ! Use provided gradients
                D_xx = grad_x_vf(1)%sf(j, k, l)
                if (n > 0) then
                    D_yy = grad_y_vf(2)%sf(j, k, l)
                    D_xy = 0.5_wp*(grad_y_vf(1)%sf(j, k, l) + &
                                   grad_x_vf(2)%sf(j, k, l))
                else
                    D_yy = 0._wp; D_xy = 0._wp
                end if
                if (p > 0) then
                    D_zz = grad_z_vf(3)%sf(j, k, l)
                    D_xz = 0.5_wp*(grad_z_vf(1)%sf(j, k, l) + &
                                   grad_x_vf(3)%sf(j, k, l))
                    D_yz = 0.5_wp*(grad_z_vf(2)%sf(j, k, l) + &
                                   grad_y_vf(3)%sf(j, k, l))
                else
                    D_zz = 0._wp; D_xz = 0._wp; D_yz = 0._wp
                end if
            else
                ! Compute gradients from primitive variables
                call s_compute_velocity_gradients_at_cell( &
                    q_prim_vf, j, k, l, D_xx, D_yy, D_zz, D_xy, D_xz, D_yz)
            end if

            ! Compute shear rate
            shear_rate = f_compute_shear_rate_from_components( &
                         D_xx, D_yy, D_zz, D_xy, D_xz, D_yz)

            ! For each phase, compute Re_visc
            $:GPU_LOOP(parallelism='[seq]')
            do q = 1, num_fluids
                if (fluid_pp(q)%non_newtonian) then
                    ! Non-Newtonian: compute shear mu from HB model
                    mu_fluid = f_compute_hb_viscosity( &
                               fluid_pp(q)%tau0, fluid_pp(q)%K, &
                               fluid_pp(q)%nn, fluid_pp(q)%mu_min, &
                               fluid_pp(q)%mu_max, shear_rate, &
                               fluid_pp(q)%hb_m)
                    Re_visc_per_phase(q, 1) = 1._wp/mu_fluid
                    ! Bulk viscosity
                    if (fluid_pp(q)%mu_bulk > 0._wp) then
                        Re_visc_per_phase(q, 2) = 1._wp/fluid_pp(q)%mu_bulk
                    else
                        Re_visc_per_phase(q, 2) = dflt_real
                    end if
                else
                    ! Newtonian: return Re input values
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, 2
                        if (Re_size(i) > 0) then
                            $:GPU_LOOP(parallelism='[seq]')
                            do r = 1, Re_size(i)
                                if (Re_idx(i, r) == q) then
                                    Re_visc_per_phase(q, i) = fluid_pp(q)%Re(i)
                                    exit
                                end if
                            end do
                        end if
                    end do
                end if
            end do
        else
            ! All Newtonian: return Re input values
            $:GPU_LOOP(parallelism='[seq]')
            do i = 1, 2
                $:GPU_LOOP(parallelism='[seq]')
                do q = 1, Re_size(i)
                    Re_visc_per_phase(Re_idx(i, q), i) = &
                        fluid_pp(Re_idx(i, q))%Re(i)
                end do
            end do
        end if

    end subroutine s_compute_re_visc

    !> Computes mixture Reynolds number from per-phase values and volume fractions.
    !!      Re_mix(i) = 1 / sum_q(alpha(q) / Re_per_phase(q, i))
    !! @param alpha Volume fractions
    !! @param Re_per_phase Per-phase Re_visc = 1/mu
    !! @param Re_mix Output: mixture Re (shear and bulk)
    pure subroutine s_compute_mixture_re(alpha, Re_per_phase, Re_mix)
        $:GPU_ROUTINE(parallelism='[seq]')

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(in) :: alpha
            real(wp), dimension(3, 2), intent(in) :: Re_per_phase
        #:else
            real(wp), dimension(num_fluids), intent(in) :: alpha
            real(wp), dimension(num_fluids, 2), intent(in) :: Re_per_phase
        #:endif
        real(wp), dimension(2), intent(out) :: Re_mix

        integer :: i, q

        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, 2
            Re_mix(i) = 0._wp
            $:GPU_LOOP(parallelism='[seq]')
            do q = 1, num_fluids
                if (Re_per_phase(q, i) /= dflt_real &
                    .and. Re_per_phase(q, i) > sgm_eps) then
                    Re_mix(i) = Re_mix(i) + alpha(q)/Re_per_phase(q, i)
                end if
            end do
            Re_mix(i) = 1._wp/max(Re_mix(i), sgm_eps)
        end do

    end subroutine s_compute_mixture_re

end module m_re_visc
