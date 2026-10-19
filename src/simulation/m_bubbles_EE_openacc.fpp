!>
!! @file
!! @brief Case-specialized OpenACC Euler--Euler bubble source kernels

#:include 'macros.fpp'
#:set MFC_BUBBLE_VECTOR_LENGTH = MFC_BUBBLE_VECTOR_LENGTH if defined('MFC_BUBBLE_VECTOR_LENGTH') else 128

!> @brief Provides the optimized OpenACC Euler--Euler bubble source kernel.
module m_bubbles_EE_openacc

    use m_precision_select, only: stp, wp

    implicit none

    private

#ifdef MFC_FUSED_BUBBLE_SOURCE
    public :: s_compute_bubble_EE_source_rp_one
#endif

contains

#ifdef MFC_FUSED_BUBBLE_SOURCE
    subroutine s_compute_bubble_EE_source_rp_one(m_in, n_in, p_in, q_lo, q_shape, q_cons_rho, q_cons_v, q_prim_r, q_prim_v, &
        & q_prim_alf, q_prim_p, rhs_r, rhs_v, rhs_alf, weight_one, r0_one, ca_in, web_in, re_inv_in, gam_in, eu_in, pi_in, &
        & small_alf_in, web_active, re_active)

        integer, intent(in)   :: m_in, n_in, p_in, q_lo(3), q_shape(3)
        real(stp), intent(in) :: q_cons_rho(q_lo(1):q_lo(1) + q_shape(1) - 1,q_lo(2):q_lo(2) + q_shape(2) - 1, &
             & q_lo(3):q_lo(3) + q_shape(3) - 1)
        real(stp), intent(in) :: q_cons_v(q_lo(1):q_lo(1) + q_shape(1) - 1,q_lo(2):q_lo(2) + q_shape(2) - 1, &
             & q_lo(3):q_lo(3) + q_shape(3) - 1)
        real(stp), intent(in) :: q_prim_r(q_lo(1):q_lo(1) + q_shape(1) - 1,q_lo(2):q_lo(2) + q_shape(2) - 1, &
             & q_lo(3):q_lo(3) + q_shape(3) - 1)
        real(stp), intent(in) :: q_prim_v(q_lo(1):q_lo(1) + q_shape(1) - 1,q_lo(2):q_lo(2) + q_shape(2) - 1, &
             & q_lo(3):q_lo(3) + q_shape(3) - 1)
        real(stp), intent(in) :: q_prim_alf(q_lo(1):q_lo(1) + q_shape(1) - 1,q_lo(2):q_lo(2) + q_shape(2) - 1, &
             & q_lo(3):q_lo(3) + q_shape(3) - 1)
        real(stp), intent(in) :: q_prim_p(q_lo(1):q_lo(1) + q_shape(1) - 1,q_lo(2):q_lo(2) + q_shape(2) - 1, &
             & q_lo(3):q_lo(3) + q_shape(3) - 1)
        real(stp), intent(inout) :: rhs_r(0:m_in,0:n_in,0:p_in), rhs_v(0:m_in,0:n_in,0:p_in)
        real(stp), intent(inout) :: rhs_alf(0:m_in,0:n_in,0:p_in)
        real(wp), intent(in)     :: weight_one, r0_one, ca_in, web_in, re_inv_in, gam_in, eu_in, pi_in, small_alf_in
        logical, intent(in)      :: web_active, re_active
        real(wp)                 :: radius, velocity, alpha, pressure, density, r2v, r3, nbub, rddot
        real(wp)                 :: pr0, pr, pv, pp, prho, pca, pweb, pre, pgam, peu, ratio, cpbw
        integer                  :: cell, cells, rest, span_x, span_y, j, k, l

        pr0 = real(r0_one, wp)
        pca = real(ca_in, wp)
        pweb = real(web_in, wp)
        pre = real(re_inv_in, wp)
        pgam = real(gam_in, wp)
        peu = real(eu_in, wp)
        span_x = m_in + 1
        span_y = n_in + 1
        cells = span_x*span_y*(p_in + 1)
        $:GPU_PARALLEL_LOOP(private='[cell, rest, j, k, l, radius, velocity, alpha, pressure, density, r2v, r3, nbub, rddot, pr, &
                            & pv, pp, prho, ratio, cpbw]', firstprivate='[weight_one, pr0, pca, pweb, pre, pgam, peu, pi_in, &
                            & small_alf_in, web_active, re_active, span_x, span_y]', extraAccArgs='vector_length(%d)' % MFC_BUBBLE_VECTOR_LENGTH)
        do cell = 0, cells - 1
            rest = cell/span_x
            j = cell - rest*span_x
            k = mod(rest, span_y)
            l = rest/span_y
            radius = q_prim_r(j, k, l)
            velocity = q_prim_v(j, k, l)
            alpha = q_prim_alf(j, k, l)
            pressure = q_prim_p(j, k, l)
            density = q_cons_rho(j, k, l)
            r3 = weight_one*radius**3._wp
            r2v = weight_one*radius**2._wp*velocity
            nbub = (3._wp/(4._wp*pi_in))*alpha/r3

            if (alpha >= small_alf_in) then
                r2v = 4._wp*pi_in*nbub*r2v
                rhs_alf(j, k, l) = rhs_alf(j, k, l) + r2v
                pr = real(radius, wp)
                pv = real(velocity, wp)
                pp = real(pressure, wp)
                prho = real(density, wp)
                ratio = (pr0/pr)**(3._wp*pgam)
                cpbw = pca*ratio - pca + peu
                if (web_active) cpbw = cpbw + (2._wp/(pweb*pr0))*ratio
                if (web_active) cpbw = cpbw - 2._wp/(pr*pweb)
                if (re_active) cpbw = cpbw - 4._wp*pre*pv/pr
                rddot = real((-1.5_wp*pv**2._wp + (cpbw - pp)/prho)/pr, wp)

                rhs_v(j, k, l) = rhs_v(j, k, l) + nbub*rddot
                rhs_r(j, k, l) = rhs_r(j, k, l) + q_cons_v(j, k, l)
            end if
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_bubble_EE_source_rp_one
#endif
end module m_bubbles_EE_openacc
