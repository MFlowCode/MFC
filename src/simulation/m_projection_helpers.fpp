!>
!! @file
!! @brief Contains module m_projection_helpers

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief Stateless helpers of the pressure projection (m_projection): loop bounds and the face formulas its kernels share
module m_projection_helpers

    use m_derived_types
    use m_global_parameters

    implicit none

    private; public :: s_set_bounds, f_wall_face, f_cond, f_face_dvel_p, f_capillary_accel, f_normal

contains

    !> Loop bounds over the cells (0 to m, n, p), but from lo to hi along direction d
    pure subroutine s_set_bounds(d, lo, hi, bounds)

        integer, intent(in)                                :: d, lo, hi
        type(int_bounds_info), dimension(3), intent(inout) :: bounds

        bounds(:)%beg = 0
        bounds(1)%end = m
        bounds(2)%end = n
        bounds(3)%end = p
        bounds(d)%beg = lo
        bounds(d)%end = hi

    end subroutine s_set_bounds

    !> A wall face: its normal velocity vanishes
    pure function f_wall_face(code) result(is_wall)

        $:GPU_ROUTINE(function_name='f_wall_face', parallelism='[seq]', cray_inline=True)

        integer, intent(in) :: code
        logical             :: is_wall

        is_wall = code == BC_REFLECTIVE .or. code == BC_SLIP_WALL .or. code == BC_NO_SLIP_WALL

    end function f_wall_face

    !> Face conductance A_f/(rho_f*d_f) with the arithmetic face density
    pure function f_cond(rho_a, rho_b, area, dist) result(cond)

        $:GPU_ROUTINE(function_name='f_cond', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in) :: rho_a, rho_b, area, dist
        real(wp)             :: cond

        cond = 2._wp*area/(max(rho_a + rho_b, sgm_eps)*dist)

    end function f_cond

    !> Pressure correction of a face velocity, tau*grad(p)/rho_f with the operator's own conductance (area 1), so div(u_face)
    !! matches the solved pressure exactly
    pure function f_face_dvel_p(rho_a, rho_b, dist, p_a, p_b, tau) result(dvel)

        $:GPU_ROUTINE(function_name='f_face_dvel_p', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in) :: rho_a, rho_b, dist, p_a, p_b, tau
        real(wp)             :: dvel

        dvel = tau*f_cond(rho_a, rho_b, 1._wp, dist)*(p_b - p_a)

    end function f_face_dvel_p

    !> Well-balanced (Brackbill CSF) capillary acceleration of a face
    pure function f_capillary_accel(kappa_a, kappa_b, w_a, w_b, c_a, c_b, rho_a, rho_b, dist) result(acc)

        $:GPU_ROUTINE(function_name='f_capillary_accel', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in) :: kappa_a, kappa_b, w_a, w_b, c_a, c_b, rho_a, rho_b, dist
        real(wp)             :: acc

        acc = 0._wp
        if (w_a + w_b > 0._wp) acc = 2._wp*sigma*(w_a*kappa_a + w_b*kappa_b)/(w_a + w_b)*(c_b - c_a)/(dist*max(rho_a + rho_b, &
            & sgm_eps))

    end function f_capillary_accel

    !> Interface normal component, grad_d(c)/|grad(c)|; zero outside the interface band
    pure function f_normal(grad_d, grad_mag) result(nd)

        $:GPU_ROUTINE(function_name='f_normal', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in) :: grad_d, grad_mag
        real(wp)             :: nd

        nd = 0._wp
        if (grad_mag > capillary_cutoff) nd = grad_d/grad_mag

    end function f_normal

end module m_projection_helpers
