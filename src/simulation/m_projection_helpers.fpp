!>
!! @file
!! @brief Contains module m_projection_helpers

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief Stateless helpers of the pressure projection (m_projection): direction-d stencil geometry and the face formulas shared by
!! its kernels
module m_projection_helpers

    use m_derived_types
    use m_global_parameters

    implicit none

    private; public :: s_unit, s_bounds, f_last, f_along, f_width, f_face_dist, f_wall_face, f_cond, f_face_gf, &
        & f_capillary_accel, f_normal

contains

    !> Unit offset (oj, ok, ol) of direction d
    pure subroutine s_unit(d, oj, ok, ol)

        $:GPU_ROUTINE(function_name='s_unit', parallelism='[seq]', cray_inline=True)

        integer, intent(in)  :: d
        integer, intent(out) :: oj, ok, ol

        oj = merge(1, 0, d == 1); ok = merge(1, 0, d == 2); ol = merge(1, 0, d == 3)

    end subroutine s_unit

    !> Loop bounds over the cells (0 to m, n, p), but from lo to hi along direction d
    pure subroutine s_bounds(d, lo, hi, bounds)

        integer, intent(in)                                :: d, lo, hi
        type(int_bounds_info), dimension(3), intent(inout) :: bounds

        bounds(:)%beg = 0; bounds(1)%end = m; bounds(2)%end = n; bounds(3)%end = p
        bounds(d)%beg = lo; bounds(d)%end = hi

    end subroutine s_bounds

    !> Last cell index along direction d
    pure integer function f_last(d)

        $:GPU_ROUTINE(function_name='f_last', parallelism='[seq]', cray_inline=True)

        integer, intent(in) :: d

        f_last = p
        if (d == 1) f_last = m
        if (d == 2) f_last = n

    end function f_last

    !> Index of (j, k, l) along direction d
    pure integer function f_along(d, j, k, l)

        $:GPU_ROUTINE(function_name='f_along', parallelism='[seq]', cray_inline=True)

        integer, intent(in) :: d, j, k, l

        f_along = l
        if (d == 1) f_along = j
        if (d == 2) f_along = k

    end function f_along

    !> Width along direction d of cell (j, k, l)
    real(wp) function f_width(d, j, k, l)
        $:GPU_ROUTINE(function_name='f_width', parallelism='[seq]', cray_inline=True)

        integer, intent(in) :: d, j, k, l

        if (d == 1) then
            f_width = dx(j)
        else if (d == 2) then
            f_width = dy(k)
        else
            f_width = dz(l)
        end if

    end function f_width

    !> Centre-to-centre distance across the direction-d face between (j, k, l) and the next cell
    real(wp) function f_face_dist(d, j, k, l)
        $:GPU_ROUTINE(function_name='f_face_dist', parallelism='[seq]', cray_inline=True)

        integer, intent(in) :: d, j, k, l

        if (d == 1) then
            f_face_dist = 0.5_wp*(dx(j) + dx(j + 1))
        else if (d == 2) then
            f_face_dist = 0.5_wp*(dy(k) + dy(k + 1))
        else
            f_face_dist = 0.5_wp*(dz(l) + dz(l + 1))
        end if

    end function f_face_dist

    !> A wall face: its normal velocity vanishes
    pure logical function f_wall_face(code)

        $:GPU_ROUTINE(function_name='f_wall_face', parallelism='[seq]', cray_inline=True)

        integer, intent(in) :: code

        f_wall_face = code == BC_REFLECTIVE .or. code == BC_SLIP_WALL .or. code == BC_NO_SLIP_WALL

    end function f_wall_face

    !> Face conductance A_f/(rho_f*d_f) with the arithmetic face density: the inertia of a face volume straddling an interface,
    !! which a heavy phase resting on a light one needs to stay at rest
    pure function f_cond(ra, rb, area, dist) result(kf)

        $:GPU_ROUTINE(function_name='f_cond', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in) :: ra, rb, area, dist
        real(wp)             :: kf

        kf = 2._wp*area/(max(ra + rb, sgm_eps)*dist)

    end function f_cond

    !> Pressure correction of a face velocity, tau*grad(p)/rho_f with the operator's own conductance (area 1), so div(uf) matches
    !! the solved pressure exactly
    pure function f_face_gf(ra, rb, dist, pa, pb, tau) result(gf)

        $:GPU_ROUTINE(function_name='f_face_gf', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in) :: ra, rb, dist, pa, pb, tau
        real(wp)             :: gf

        gf = tau*f_cond(ra, rb, 1._wp, dist)*(pb - pa)

    end function f_face_gf

    !> Well-balanced (Brackbill CSF) capillary acceleration of a face: sigma*kappa_f*(c_b - c_a)/(d_f*rho_f), with rho_f the same
    !! arithmetic face density as the pressure operator, so a constant curvature is balanced exactly by a pressure jump. kappa_f is
    !! the |grad c|-weighted mean of the adjacent cells' curvature (weights wa, wb)
    pure function f_capillary_accel(ka, kb, wa, wb, ca, cb, ra, rb, dist) result(acc)

        $:GPU_ROUTINE(function_name='f_capillary_accel', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in) :: ka, kb, wa, wb, ca, cb, ra, rb, dist
        real(wp)             :: acc

        acc = 0._wp
        if (wa + wb > 0._wp) acc = 2._wp*sigma*(wa*ka + wb*kb)/(wa + wb)*(cb - ca)/(dist*max(ra + rb, sgm_eps))

    end function f_capillary_accel

    !> Interface normal component, grad_d(c)/|grad(c)|; zero outside the interface band
    pure function f_normal(gd, g) result(nd)

        $:GPU_ROUTINE(function_name='f_normal', parallelism='[seq]', cray_inline=True)

        real(wp), intent(in) :: gd, g
        real(wp)             :: nd

        nd = 0._wp
        if (g > capillary_cutoff) nd = gd/g

    end function f_normal

end module m_projection_helpers
