!>
!! @file
!! @brief Case-specialized OpenACC QBMM moment inversion kernels

#:include 'case.fpp'
#:include 'macros.fpp'
#:include 'qbmm_quadrature.fpp'

!> @brief Provides case-specialized OpenACC QBMM moment inversion kernels.
module m_qbmm_openacc

    use m_derived_types
    use m_global_parameters, only: nb
    use m_constants, only: bubble_model_keller_miksis, bubble_model_rayleigh_plesset

    implicit none

    private

#ifdef MFC_QBMM_MOM_INV_DEVICE
    #:if MFC_CASE_OPTIMIZATION
        public :: s_initialize_qbmm_openacc_kernels, s_mom_inv_openacc_specialized

        integer, parameter :: qbmm_nout = 6
        integer, parameter :: qbmm_max_px = 10
        integer, parameter :: qbmm_max_py = 8
        integer, parameter :: qbmm_max_r0 = 4
        integer, parameter :: qbmm_not_integral = -128
        integer, parameter :: nterms = ${nterms}$
        #:if nterms == 32 and not polytropic
            #! The km_np moment-source term tables and coefficient algebra, as fypp data so the
            #! specialized kernel can unroll every (r, j) term into literal-indexed scalar code.
            #:set qbmm_km_np_px_base = [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, -2, -2, -2, -2, -2, &
                & -2, -2, -3, -3, -2, -1, -1, -2, -1, -1, -2]
            #:set qbmm_km_np_py_base = [1, 2, 3, -1, 0, 1, -1, 0, 1, 0, 1, 1, -1, 0, 1, 0, -1, 1, 2, -1, 0, -1, 0, 0, -1, 0, &
                & -1, 0, 0, -1, 0, 0]
            #! Named scalar power-table cells for the km_np kernel.
            #:set KM_PXN = lambda k: ('m' + str(-k)) if k < 0 else str(k)
            #:set KM_XPOW_NAMES = [f'xpow_{KM_PXN(k)}_n{n}' for k in range(-3, 2) for n in (1, 2)]
            #:set KM_YPOW_NAMES = [f'ypow_{KM_PXN(k)}_n{n}' for k in range(-1, 6) for n in range(1, 5)]
            #! Pre-joined for the macro-call private list, since fypp splits call arguments on
            #! commas before evaluating inline ${} expressions.
            #:set KM_POW_PRIV = ', '.join(KM_XPOW_NAMES + KM_YPOW_NAMES)
            #! Named scalars for the chyqmom/quadrature arrays (moms, msum, wght, abscX,
            #! abscY, wght_pb, wght_mv, wght_ht, ht, rhs_mv_pre) in the km_np kernel.
            #:set KM_MOM_NAMES = [f'mom{i}' for i in range(1, 7)]
            #:set KM_QUAD_NAMES = [f'{p}{n}' for p in ('wght', 'abscx', 'abscy', 'wght_pb', 'wght_mv', 'wght_ht', 'ht', &
                & 'rhs_mv_pre', 'abscx_inv') for n in range(1, 5)]
            #:set KM_MSUM_NAMES = [f'msum{i}' for i in range(1, 7)]
            #:set KM_CHY_NAMES = ['chy_bu', 'chy_bv', 'chy_d20', 'chy_d11', 'chy_d02', 'chy_c20', 'chy_c11', 'chy_c02', &
                & 'chy_ra1', 'chy_ra2', 'chy_c2a', 'chy_ua1', 'chy_ua2', 'chy_vf1', 'chy_vf2', 'chy_mu2', 'chy_rb1', &
                & 'chy_rb2', 'chy_c2b', 'chy_ub1', 'chy_ub2']
            #! Hoisted reciprocals for the thermo chains. Denominators that are invariant
            #! (pv, pe_c, pe_t, pb0, r0), shared between consumers (abscx, one reciprocal
            #! serving rhs_mv_pre, ht and three epilogue sites), or exactly hoistable
            #! (r_v/r_g, c20) become multiplies, which round differently from the divisions.
            #:set KM_HOIST_NAMES = ['inv_pv', 'inv_pe_c', 'ratio_rvg', 'inv_pb0', 'inv_r0', 'inv_pe_t', 'inv_c20']
            #:set KM_SCALAR_NAMES = KM_MOM_NAMES + KM_QUAD_NAMES + KM_MSUM_NAMES + KM_CHY_NAMES + KM_HOIST_NAMES
            #:set KM_SCALAR_PRIV = ', '.join(KM_SCALAR_NAMES)
            #! (i1, i2) per moment-output row r, and the per-term sets that drive the unrolled filters.
            #:set KM_RI = [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (0, 2)]
            #:set KM_RE_SET = {16, 18, 19, 20, 21, 22, 23, 24, 25, 26, 29, 32}
            #:set KM_WEB_SET = {17, 25}
            #:set KM_PB_QUAD = {7, 8, 9, 10, 11, 22, 23, 26}
            #:set KM_MV_QUAD = {27, 28, 29}
            #:set KM_HT_QUAD = {30, 31, 32}
            #! f_coeff_km_np's select-case as literal expressions; i1/i2 interpolate as integers.
            #:set KM_COEFF = lambda j, i1, i2: { &
                & 1: f'-3._wp*{i2}/2._wp', &
                & 2: f'-{i2}/c', &
                & 3: f'{i2}/(2._wp*c*c)', &
                & 4: f'-{i2}*pres/rho', &
                & 5: f'-2._wp*{i2}*pres/(c*rho)', &
                & 6: f'-{i2}*pres/(c*c*rho)', &
                & 7: f'{i2}/rho', &
                & 8: f'2._wp*{i2}/(c*rho)', &
                & 9: f'{i2}/(c*c*rho)', &
                & 10: f'-3._wp*{i2}*qbmm_gam/(c*rho)', &
                & 11: f'-3._wp*{i2}*qbmm_gam/(c*c*rho)', &
                & 12: f'{i1}', &
                & 16: f'-4._wp*{i2}*qbmm_re_inv/rho', &
                & 17: f'-2._wp*{i2}/qbmm_web/rho', &
                & 18: f'6._wp*{i2}*qbmm_re_inv/(rho*c)', &
                & 19: f'-2._wp*{i2}*qbmm_re_inv/(rho*c*c)', &
                & 20: f'4._wp*{i2}*pres*qbmm_re_inv/(rho*rho*c)', &
                & 21: f'4._wp*{i2}*pres*qbmm_re_inv/(rho*rho*c*c)', &
                & 22: f'-4._wp*{i2}*qbmm_re_inv/(rho*rho*c)', &
                & 23: f'-4._wp*{i2}*qbmm_re_inv/(rho*rho*c*c)', &
                & 24: f'16._wp*{i2}*qbmm_re_inv*qbmm_re_inv/(rho*rho*c)', &
                & 25: f'8._wp*{i2}*qbmm_re_inv/qbmm_web/(rho*rho*c)', &
                & 26: f'-12._wp*{i2}*qbmm_gam*qbmm_re_inv/(rho*rho*c*c)', &
                & 27: f'3._wp*{i2}*qbmm_gam*qbmm_r_v*qbmm_tw/(c*rho)', &
                & 28: f'3._wp*{i2}*qbmm_gam*qbmm_r_v*qbmm_tw/(c*c*rho)', &
                & 29: f'12._wp*{i2}*qbmm_gam*qbmm_r_v*qbmm_tw*qbmm_re_inv/(rho*rho*c*c)', &
                & 30: f'3._wp*{i2}*qbmm_gam/(c*rho)', &
                & 31: f'3._wp*{i2}*qbmm_gam/(c*c*rho)', &
                & 32: f'12._wp*{i2}*qbmm_gam*qbmm_re_inv/(rho*rho*c*c)', &
                & }[j]
        #:endif

        real(wp) :: qbmm_px_exp(qbmm_max_px, nb), qbmm_py_exp(qbmm_max_py, nb), qbmm_r0_exp(qbmm_max_r0)
        integer  :: qbmm_px_n(nb), qbmm_py_n(nb), qbmm_px_int(qbmm_max_px, nb), qbmm_py_int(qbmm_max_py, nb)
        integer  :: qbmm_r0_n, qbmm_r0_int(qbmm_max_r0)
        integer  :: qbmm_px_id(qbmm_nout, nterms, nb), qbmm_py_id(qbmm_nout, nterms, nb)
        integer  :: qbmm_r0_id(qbmm_nout, nterms, nb), qbmm_live(qbmm_nout, nterms)
        integer  :: qbmm_bubmoms(nb, 6), qbmm_re_active, qbmm_web_active
        integer  :: qbmm_alf_idx, qbmm_pres_idx, qbmm_rho_idx, qbmm_nbub_idx
        real(wp) :: qbmm_weight(nb), qbmm_r0(nb), qbmm_pb0(nb)
        real(wp) :: qbmm_k_v(nb), qbmm_k_g(nb), qbmm_pe_t(nb), qbmm_re_trans_c(nb), qbmm_re_trans_t(nb)
        real(wp) :: qbmm_mass_g0(nb), qbmm_mass_v0(nb)
        real(wp) :: qbmm_gam, qbmm_pv, qbmm_r_v, qbmm_r_g, qbmm_m_v, qbmm_m_g
        real(wp) :: qbmm_phi_vg, qbmm_phi_gv, qbmm_tw, qbmm_pe_c, qbmm_re_inv, qbmm_web
        real(wp) :: qbmm_sgm_eps, qbmm_small_alf, qbmm_pi, qbmm_gammas1, qbmm_pi_infs1

        $:GPU_DECLARE(copyin='[qbmm_px_exp, qbmm_py_exp, qbmm_r0_exp, qbmm_px_n, qbmm_py_n, qbmm_px_int, qbmm_py_int, qbmm_r0_n, &
                      & qbmm_r0_int, qbmm_px_id, qbmm_py_id, qbmm_r0_id, qbmm_live, qbmm_bubmoms, qbmm_re_active, &
                      & qbmm_web_active, qbmm_alf_idx, qbmm_pres_idx, qbmm_rho_idx, qbmm_nbub_idx, qbmm_weight, qbmm_r0, &
                      & qbmm_pb0, qbmm_k_v, qbmm_k_g, qbmm_pe_t, qbmm_re_trans_c, qbmm_re_trans_t, qbmm_mass_g0, qbmm_mass_v0, &
                      & qbmm_gam, qbmm_pv, qbmm_r_v, qbmm_r_g, qbmm_m_v, qbmm_m_g, qbmm_phi_vg, qbmm_phi_gv, qbmm_tw, qbmm_pe_c, &
                      & qbmm_re_inv, qbmm_web, qbmm_sgm_eps, qbmm_small_alf, qbmm_pi, qbmm_gammas1, qbmm_pi_infs1]')
    #:endif
#endif

contains

#ifdef MFC_QBMM_MOM_INV_DEVICE
    #:if MFC_CASE_OPTIMIZATION
        integer function f_integer_exponent(value)

            real(wp), intent(in) :: value
            integer              :: exponent

            exponent = nint(value)
            if (real(exponent, wp) == value .and. abs(exponent) <= 8) then
                f_integer_exponent = exponent
            else
                f_integer_exponent = qbmm_not_integral
            end if

        end function f_integer_exponent

        subroutine s_intern_exp(value, table, count, slot)

            real(wp), intent(in)    :: value
            real(wp), intent(inout) :: table(:)
            integer, intent(inout)  :: count
            integer, intent(out)    :: slot
            integer                 :: k

            slot = 0
            do k = 1, count
                if (table(k) == value) then
                    slot = k
                    exit
                end if
            end do
            if (slot == 0) then
                count = count + 1
                slot = count
                table(slot) = value
            end if

        end subroutine s_intern_exp

        logical function f_term_is_live(model, poly, re_active, web_active, j, i1, i2)

            integer, intent(in) :: model, re_active, web_active, j, i1, i2
            logical, intent(in) :: poly

            f_term_is_live = .false.
            if (model == bubble_model_rayleigh_plesset) then
                select case (j)
                case (1:3); f_term_is_live = i2 /= 0
                case (4); f_term_is_live = i1 /= 0
                case (5); f_term_is_live = re_active == 1 .and. i2 /= 0
                case (6); f_term_is_live = web_active == 1 .and. i2 /= 0
                case (7); f_term_is_live = poly .and. i2 /= 0
                end select
                return
            end if

            if (j == 12) then
                f_term_is_live = i1 /= 0
            else if (i2 == 0) then
                return
            else if (j <= 11) then
                f_term_is_live = .true.
            else if (j <= 15) then
                f_term_is_live = poly
            else if (j == 16) then
                f_term_is_live = re_active == 1
            else if (j == 17) then
                f_term_is_live = web_active == 1
            else if (j <= 24) then
                f_term_is_live = re_active == 1
            else if (j == 25) then
                f_term_is_live = re_active == 1 .and. web_active == 1
            else if (j == 26) then
                f_term_is_live = re_active == 1
            else if (j == 27 .or. j == 28 .or. j == 30 .or. j == 31) then
                f_term_is_live = .not. poly
            else if (j == 29 .or. j == 32) then
                f_term_is_live = .not. poly .and. re_active == 1
            end if

        end function f_term_is_live

        subroutine s_initialize_qbmm_openacc_kernels(momrhs_in, bubmoms_in, weight_in, r0_in, pb0_in, gam_in, pv_in, r_v_in, &
            & r_g_in, m_v_in, m_g_in, phi_vg_in, phi_gv_in, tw_in, pe_c_in, re_inv_in, web_in, sgm_eps_in, small_alf_in, pi_in, &
            & gammas1_in, pi_infs1_in, model_in, poly_in, re_active_in, web_active_in, alf_idx_in, pres_idx_in, rho_idx_in, &
            & nbub_idx_in, k_v_in, k_g_in, pe_t_in, re_trans_c_in, re_trans_t_in, mass_g0_in, mass_v0_in)

            real(wp), intent(in)           :: momrhs_in(3,0:2,0:2,nterms, nb)
            integer, intent(in)            :: bubmoms_in(nb, 6)
            real(wp), intent(in)           :: weight_in(nb), r0_in(nb), pb0_in(nb)
            real(wp), intent(in)           :: gam_in, pv_in, r_v_in, r_g_in, m_v_in, m_g_in, phi_vg_in, phi_gv_in
            real(wp), intent(in)           :: tw_in, pe_c_in, re_inv_in, web_in, sgm_eps_in, small_alf_in, pi_in
            real(wp), intent(in)           :: gammas1_in, pi_infs1_in
            integer, intent(in)            :: model_in, re_active_in, web_active_in
            logical, intent(in)            :: poly_in
            integer, intent(in)            :: alf_idx_in, pres_idx_in, rho_idx_in, nbub_idx_in
            real(wp), intent(in), optional :: k_v_in(nb), k_g_in(nb), pe_t_in(nb), re_trans_c_in(nb)
            real(wp), intent(in), optional :: re_trans_t_in(nb), mass_g0_in(nb), mass_v0_in(nb)
            integer                        :: i1, i2, j, q, t

            qbmm_bubmoms = bubmoms_in
            qbmm_weight = weight_in; qbmm_r0 = r0_in; qbmm_pb0 = pb0_in
            qbmm_k_v = 0._wp; qbmm_k_g = 0._wp; qbmm_pe_t = 0._wp
            qbmm_re_trans_c = 0._wp; qbmm_re_trans_t = 0._wp
            qbmm_mass_g0 = 0._wp; qbmm_mass_v0 = 0._wp
            if (present(k_v_in)) then
                qbmm_k_v = k_v_in; qbmm_k_g = k_g_in; qbmm_pe_t = pe_t_in
                qbmm_re_trans_c = re_trans_c_in; qbmm_re_trans_t = re_trans_t_in
                qbmm_mass_g0 = mass_g0_in; qbmm_mass_v0 = mass_v0_in
            end if
            qbmm_gam = gam_in; qbmm_pv = pv_in; qbmm_r_v = r_v_in; qbmm_r_g = r_g_in
            qbmm_m_v = m_v_in; qbmm_m_g = m_g_in; qbmm_phi_vg = phi_vg_in; qbmm_phi_gv = phi_gv_in
            qbmm_tw = tw_in; qbmm_pe_c = pe_c_in; qbmm_re_inv = re_inv_in; qbmm_web = web_in
            qbmm_sgm_eps = sgm_eps_in; qbmm_small_alf = small_alf_in; qbmm_pi = pi_in
            qbmm_gammas1 = gammas1_in; qbmm_pi_infs1 = pi_infs1_in
            qbmm_re_active = re_active_in; qbmm_web_active = web_active_in
            qbmm_alf_idx = alf_idx_in; qbmm_pres_idx = pres_idx_in
            qbmm_rho_idx = rho_idx_in; qbmm_nbub_idx = nbub_idx_in

            qbmm_px_exp = 0._wp; qbmm_py_exp = 0._wp; qbmm_r0_exp = 0._wp
            qbmm_px_n = 1; qbmm_py_n = 1; qbmm_r0_n = 1
            qbmm_px_int = 0; qbmm_py_int = 0; qbmm_r0_int = 0
            qbmm_px_id = 1; qbmm_py_id = 1; qbmm_r0_id = 1; qbmm_live = 0

            t = 0
            do i2 = 0, 2
                do i1 = 0, 2
                    if (i1 + i2 > 2) cycle
                    t = t + 1
                    do j = 1, nterms
                        if (.not. f_term_is_live(model_in, poly_in, re_active_in, web_active_in, j, i1, i2)) cycle
                        qbmm_live(t, j) = 1
                        do q = 1, nb
                            call s_intern_exp(momrhs_in(1, i1, i2, j, q), qbmm_px_exp(:,q), qbmm_px_n(q), qbmm_px_id(t, j, q))
                            call s_intern_exp(momrhs_in(2, i1, i2, j, q), qbmm_py_exp(:,q), qbmm_py_n(q), qbmm_py_id(t, j, q))
                            call s_intern_exp(momrhs_in(3, i1, i2, j, q), qbmm_r0_exp, qbmm_r0_n, qbmm_r0_id(t, j, q))
                        end do
                    end do
                end do
            end do

            do q = 1, nb
                do j = 1, qbmm_px_n(q)
                    qbmm_px_int(j, q) = f_integer_exponent(qbmm_px_exp(j, q))
                end do
                do j = 1, qbmm_py_n(q)
                    qbmm_py_int(j, q) = f_integer_exponent(qbmm_py_exp(j, q))
                end do
            end do
            do j = 1, qbmm_r0_n
                qbmm_r0_int(j) = f_integer_exponent(qbmm_r0_exp(j))
            end do

            $:GPU_UPDATE(device='[qbmm_px_exp, qbmm_py_exp, qbmm_r0_exp, qbmm_px_n, qbmm_py_n, qbmm_px_int, qbmm_py_int, &
                         & qbmm_r0_n, qbmm_r0_int, qbmm_px_id, qbmm_py_id, qbmm_r0_id, qbmm_live, qbmm_bubmoms, qbmm_re_active, &
                         & qbmm_web_active, qbmm_alf_idx, qbmm_pres_idx, qbmm_rho_idx, qbmm_nbub_idx, qbmm_weight, qbmm_r0, &
                         & qbmm_pb0, qbmm_k_v, qbmm_k_g, qbmm_pe_t, qbmm_re_trans_c, qbmm_re_trans_t, qbmm_mass_g0, qbmm_mass_v0, &
                         & qbmm_gam, qbmm_pv, qbmm_r_v, qbmm_r_g, qbmm_m_v, qbmm_m_g, qbmm_phi_vg, qbmm_phi_gv, qbmm_tw, &
                         & qbmm_pe_c, qbmm_re_inv, qbmm_web, qbmm_sgm_eps, qbmm_small_alf, qbmm_pi, qbmm_gammas1, qbmm_pi_infs1]')

        end subroutine s_initialize_qbmm_openacc_kernels

        real(wp) function f_integral_pow(value, exponent)
            $:GPU_ROUTINE(function_name='f_integral_pow', parallelism='[seq]')

            real(wp), intent(in) :: value
            integer, intent(in)  :: exponent
            real(wp)             :: base, result
            integer              :: power

            power = abs(exponent)
            base = value
            result = 1._wp
            do while (power > 0)
                if (iand(power, 1) /= 0) result = result*base
                power = ishft(power, -1)
                if (power > 0) base = base*base
            end do
            f_integral_pow = merge(1._wp/result, result, exponent < 0)

        end function f_integral_pow

        subroutine s_hyqmom(frho, fup, fmom, sgm_eps_in)

            $:GPU_ROUTINE(function_name='s_hyqmom', parallelism='[seq]')

            real(wp), intent(inout) :: frho(2), fup(2)
            real(wp), intent(in)    :: fmom(3)
            real(wp), intent(in)    :: sgm_eps_in

            @:QBMM_HYQMOM_BODY(sgm_eps_in)

        end subroutine s_hyqmom

        subroutine s_chyqmom(momin, wght, abscX, abscY, sgm_eps_in)

            $:GPU_ROUTINE(function_name='s_chyqmom', parallelism='[seq]')

            real(wp), intent(in)    :: momin(6)
            real(wp), intent(inout) :: wght(4), abscX(4), abscY(4)
            real(wp), intent(in)    :: sgm_eps_in

            @:QBMM_CHYQMOM_BODY(sgm_eps_in)

        end subroutine s_chyqmom

        #:def QBMM_COEFF(NAME, IS_KM, IS_POLY)
            real(wp) function f_coeff_${NAME}$(j, i1, i2, pres, rho, c, re_active, web_active, re_inv, web, pv, gam, r_v, tw)
                $:GPU_ROUTINE(parallelism='[seq]')

                integer, intent(in)  :: j, i1, i2
                real(wp), intent(in) :: pres, rho, c
                integer, intent(in)  :: re_active, web_active
                real(wp), intent(in) :: re_inv, web, pv, gam, r_v, tw

                f_coeff_${NAME}$ = 0._wp
                #:if not IS_KM
                    select case (j)
                    case (1); f_coeff_${NAME}$ = -i2*pres/rho
                    case (2); f_coeff_${NAME}$ = -3._wp*i2/2._wp
                    case (3); f_coeff_${NAME}$ = i2/rho
                    case (4); f_coeff_${NAME}$ = i1
                    case (5)
                        if (re_active == 1) f_coeff_${NAME}$ = -4._wp*i2*re_inv/rho
                    case (6)
                        if (web_active == 1) f_coeff_${NAME}$ = -2._wp*i2/web/rho
                        #:if IS_POLY
                        case (7); f_coeff_${NAME}$ = i2*pv/rho
                        #:endif
                    end select
                #:else
                    select case (j)
                    case (1); f_coeff_${NAME}$ = -3._wp*i2/2._wp
                    case (2); f_coeff_${NAME}$ = -i2/c
                    case (3); f_coeff_${NAME}$ = i2/(2._wp*c*c)
                    case (4); f_coeff_${NAME}$ = -i2*pres/rho
                    case (5); f_coeff_${NAME}$ = -2._wp*i2*pres/(c*rho)
                    case (6); f_coeff_${NAME}$ = -i2*pres/(c*c*rho)
                    case (7); f_coeff_${NAME}$ = i2/rho
                    case (8); f_coeff_${NAME}$ = 2._wp*i2/(c*rho)
                    case (9); f_coeff_${NAME}$ = i2/(c*c*rho)
                    case (10); f_coeff_${NAME}$ = -3._wp*i2*gam/(c*rho)
                    case (11); f_coeff_${NAME}$ = -3._wp*i2*gam/(c*c*rho)
                    case (12); f_coeff_${NAME}$ = i1
                        #:if IS_POLY
                        case (13); f_coeff_${NAME}$ = i2*pv/rho
                        case (14); f_coeff_${NAME}$ = 2._wp*i2*pv/(c*rho)
                        case (15); f_coeff_${NAME}$ = i2*pv/(c*c*rho)
                        #:endif
                    case (16)
                        if (re_active == 1) f_coeff_${NAME}$ = -4._wp*i2*re_inv/rho
                    case (17)
                        if (web_active == 1) f_coeff_${NAME}$ = -2._wp*i2/web/rho
                    case (18)
                        if (re_active == 1) f_coeff_${NAME}$ = 6._wp*i2*re_inv/(rho*c)
                    case (19)
                        if (re_active == 1) f_coeff_${NAME}$ = -2._wp*i2*re_inv/(rho*c*c)
                    case (20)
                        if (re_active == 1) f_coeff_${NAME}$ = 4._wp*i2*pres*re_inv/(rho*rho*c)
                    case (21)
                        if (re_active == 1) f_coeff_${NAME}$ = 4._wp*i2*pres*re_inv/(rho*rho*c*c)
                    case (22)
                        if (re_active == 1) f_coeff_${NAME}$ = -4._wp*i2*re_inv/(rho*rho*c)
                    case (23)
                        if (re_active == 1) f_coeff_${NAME}$ = -4._wp*i2*re_inv/(rho*rho*c*c)
                    case (24)
                        if (re_active == 1) f_coeff_${NAME}$ = 16._wp*i2*re_inv*re_inv/(rho*rho*c)
                    case (25)
                        if (re_active == 1 .and. web_active == 1) f_coeff_${NAME}$ = 8._wp*i2*re_inv/web/(rho*rho*c)
                    case (26)
                        if (re_active == 1) f_coeff_${NAME}$ = -12._wp*i2*gam*re_inv/(rho*rho*c*c)
                        #:if not IS_POLY
                        case (27); f_coeff_${NAME}$ = 3._wp*i2*gam*r_v*tw/(c*rho)
                        case (28); f_coeff_${NAME}$ = 3._wp*i2*gam*r_v*tw/(c*c*rho)
                        case (29)
                            if (re_active == 1) f_coeff_${NAME}$ = 12._wp*i2*gam*r_v*tw*re_inv/(rho*rho*c*c)
                        case (30); f_coeff_${NAME}$ = 3._wp*i2*gam/(c*rho)
                        case (31); f_coeff_${NAME}$ = 3._wp*i2*gam/(c*c*rho)
                        case (32)
                            if (re_active == 1) f_coeff_${NAME}$ = 12._wp*i2*gam*re_inv/(rho*rho*c*c)
                        #:endif
                    end select
                #:endif

            end function f_coeff_${NAME}$
        #:enddef

        #! The km_np coefficient algebra is inlined into its specialized kernel's unrolled
        #! term loop, so it has no runtime f_coeff function.
        #:if nterms == 32 and polytropic
            $:QBMM_COEFF('km_poly', True, True)
        #:elif nterms != 32 and polytropic
            $:QBMM_COEFF('rp_poly', False, True)
        #:elif nterms != 32
            $:QBMM_COEFF('rp_np', False, False)
        #:endif

        #:def QBMM_KERNEL(NAME, IS_KM, IS_POLY)
            subroutine s_mom_inv_${NAME}$(q_cons_vf, q_prim_vf, momsp, moms3d, pb, rhs_pb, mv, rhs_mv, ix_beg, ix_end, iy_beg, &
                                          & iy_end, iz_beg, iz_end, pb_l1, pb_l2, pb_l3, qbmm_re_active, qbmm_web_active, &
                                          & qbmm_alf_idx, qbmm_pres_idx, qbmm_rho_idx, qbmm_nbub_idx, qbmm_r0_n, qbmm_gam, &
                                          & qbmm_pv, qbmm_r_v, qbmm_r_g, qbmm_m_v, qbmm_m_g, qbmm_phi_vg, qbmm_phi_gv, qbmm_tw, &
                                          & qbmm_pe_c, qbmm_re_inv, qbmm_web, qbmm_sgm_eps, qbmm_small_alf, qbmm_pi, &
                                          & qbmm_gammas1, qbmm_pi_infs1)

                integer, intent(in)               :: ix_beg, ix_end, iy_beg, iy_end, iz_beg, iz_end, pb_l1, pb_l2, pb_l3
                type(scalar_field), intent(inout) :: q_cons_vf(:), q_prim_vf(:), momsp(:)
                type(scalar_field), intent(inout) :: moms3d(0:,0:,:)
                real(stp), intent(inout)          :: pb(pb_l1:,pb_l2:,pb_l3:,:,:), mv(pb_l1:,pb_l2:,pb_l3:,:,:)
                real(wp), intent(inout)           :: rhs_pb(pb_l1:,pb_l2:,pb_l3:,:,:), rhs_mv(pb_l1:,pb_l2:,pb_l3:,:,:)
                integer, value, intent(in)        :: qbmm_re_active, qbmm_web_active, qbmm_alf_idx, qbmm_pres_idx
                integer, value, intent(in)        :: qbmm_rho_idx, qbmm_nbub_idx, qbmm_r0_n
                real(wp), value, intent(in)       :: qbmm_gam, qbmm_pv, qbmm_r_v, qbmm_r_g, qbmm_m_v, qbmm_m_g
                real(wp), value, intent(in)       :: qbmm_phi_vg, qbmm_phi_gv, qbmm_tw, qbmm_pe_c, qbmm_re_inv, qbmm_web
                real(wp), value, intent(in)       :: qbmm_sgm_eps, qbmm_small_alf, qbmm_pi, qbmm_gammas1, qbmm_pi_infs1
                real(wp)                          :: moms(6), msum(6), wght(4), abscX(4), abscY(4), wght_pb(4)
                real(wp)                          :: wght_mv(4), wght_ht(4), ht(4), rhs_mv_pre(4)

                #:if IS_KM and not IS_POLY
                    $:'real(wp) :: ' + ' &\n&, '.join(KM_XPOW_NAMES + KM_YPOW_NAMES)
                    $:'real(wp) :: ' + ' &\n&, '.join(KM_SCALAR_NAMES)
                #:else
                    real(wp) :: r0p(qbmm_max_r0), xpe(qbmm_max_px, 2), ypowers(qbmm_max_py, 4)
                #:endif
                real(wp) :: pres, rho, nbub, c, alf, momsum, drdt, drdt2, chi_vw, x_vw, rho_mw, k_mw, grad_T
                real(wp) :: quad, coeff, xp1, xp2, yp1, yp2, yp3, yp4, n_tait, B_tait
                real(wp) :: rv1, rv2, rv3, rv4, rv5, fq1, fq2, fq3, fq4, fq5
                integer  :: cell, cells, rest, span_x, span_y
                integer  :: id1, id2, id3, i1, i2, j, q, r, k, node, kx, ky, kr

                span_x = ix_end - ix_beg + 1
                span_y = iy_end - iy_beg + 1
                cells = span_x*span_y*(iz_end - iz_beg + 1)
                #:if IS_KM and not IS_POLY
                    $:GPU_PARALLEL_LOOP(private='[cell, rest, id1, id2, id3, &
                                        & ' + KM_SCALAR_PRIV + ', ' + KM_POW_PRIV + ', r, q, k, node, n_tait, B_tait, pres, rho, nbub, c, alf, momsum, drdt, drdt2, chi_vw, x_vw, rho_mw, k_mw, grad_T, quad, coeff, xp1, xp2, yp1, yp2, yp3, yp4, i1, i2, j, rv1, rv2, rv3, rv4, rv5, fq1, fq2, fq3, fq4, fq5]', firstprivate='[span_x, span_y, ix_beg, iy_beg, iz_beg]')
                #:else
                    $:GPU_PARALLEL_LOOP(private='[cell, rest, id1, id2, id3, moms, msum, wght, abscX, abscY, wght_pb, wght_mv, &
                                        & wght_ht, ht, rhs_mv_pre, r0p, xpe, ypowers, r, q, k, node, kx, ky, kr, n_tait, B_tait, &
                                        & pres, rho, nbub, c, alf, momsum, drdt, drdt2, chi_vw, x_vw, rho_mw, k_mw, grad_T, quad, &
                                        & coeff, xp1, xp2, yp1, yp2, yp3, yp4, i1, i2, j, rv1, rv2, rv3, rv4, rv5, fq1, fq2, fq3, &
                                        & fq4, fq5]', firstprivate='[span_x, span_y, ix_beg, iy_beg, iz_beg]')
                #:endif
                do cell = 0, cells - 1
                    rest = cell/span_x
                    id1 = ix_beg + cell - rest*span_x
                    id2 = iy_beg + mod(rest, span_y)
                    id3 = iz_beg + rest/span_y
                    alf = q_prim_vf(qbmm_alf_idx)%sf(id1, id2, id3)
                    if (alf > qbmm_small_alf) then
                        pres = q_prim_vf(qbmm_pres_idx)%sf(id1, id2, id3)
                        rho = q_prim_vf(qbmm_rho_idx)%sf(id1, id2, id3)
                        #:if IS_KM
                            n_tait = 1._wp/qbmm_gammas1 + 1._wp
                            B_tait = qbmm_pi_infs1*(n_tait - 1)/n_tait
                            c = n_tait*(pres + B_tait)*(1._wp - alf)/rho
                            c = merge(sqrt(c), qbmm_sgm_eps, c > 0._wp)
                        #:endif
                        #:if IS_KM and not IS_POLY
                            ! Cell-invariant reciprocals for the thermo chains; the ratio itself is hoisted bit-exactly.
                            inv_pv = 1._wp/qbmm_pv
                            inv_pe_c = 1._wp/qbmm_pe_c
                            ratio_rvg = qbmm_r_v/qbmm_r_g
                        #:endif
                        nbub = q_cons_vf(qbmm_nbub_idx)%sf(id1, id2, id3)
                        fq1 = 0._wp; fq2 = 0._wp; fq3 = 0._wp; fq4 = 0._wp; fq5 = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do q = 1, nb
                            #:if IS_KM and not IS_POLY
                                ! s_chyqmom/s_hyqmom inlined with scalar actuals; every
                                ! statement below mirrors the routines' operation order
                                ! (hyqmom's moment vector is always (1, 0, c2), so its
                                ! divisions by 1 and subtractions of 0 fold away exactly).
                                mom2 = q_prim_vf(qbmm_bubmoms(q, 2))%sf(id1, id2, id3)
                                mom3 = q_prim_vf(qbmm_bubmoms(q, 3))%sf(id1, id2, id3)
                                mom4 = q_prim_vf(qbmm_bubmoms(q, 4))%sf(id1, id2, id3)
                                mom5 = q_prim_vf(qbmm_bubmoms(q, 5))%sf(id1, id2, id3)
                                mom6 = q_prim_vf(qbmm_bubmoms(q, 6))%sf(id1, id2, id3)
                                mom1 = 1._wp
                                chy_bu = mom2/mom1; chy_bv = mom3/mom1
                                chy_d20 = mom4/mom1; chy_d11 = mom5/mom1; chy_d02 = mom6/mom1
                                chy_c20 = chy_d20 - chy_bu**2._wp; chy_c11 = chy_d11 - chy_bu*chy_bv
                                chy_c02 = chy_d02 - chy_bv**2._wp
                                chy_ra1 = 1._wp/2._wp; chy_ra2 = 1._wp/2._wp
                                chy_c2a = max(chy_c20, qbmm_sgm_eps)
                                chy_ua1 = 0._wp - sqrt(chy_c2a); chy_ua2 = 0._wp + sqrt(chy_c2a)
                                inv_c20 = 1._wp/chy_c20
                                chy_vf1 = chy_c11*chy_ua1*inv_c20; chy_vf2 = chy_c11*chy_ua2*inv_c20
                                chy_mu2 = max(0._wp, chy_c02 - (chy_ra1*chy_vf1**2._wp + chy_ra2*chy_vf2**2._wp))
                                chy_rb1 = 1._wp/2._wp; chy_rb2 = 1._wp/2._wp
                                chy_c2b = max(chy_mu2, qbmm_sgm_eps)
                                chy_ub1 = 0._wp - sqrt(chy_c2b); chy_ub2 = 0._wp + sqrt(chy_c2b)
                                wght1 = mom1*(chy_ra1*chy_rb1); wght2 = mom1*(chy_ra1*chy_rb2)
                                wght3 = mom1*(chy_ra2*chy_rb1); wght4 = mom1*(chy_ra2*chy_rb2)
                                abscx1 = chy_bu + chy_ua1; abscx2 = chy_bu + chy_ua1
                                abscx3 = chy_bu + chy_ua2; abscx4 = chy_bu + chy_ua2
                                abscy1 = chy_bv + (chy_vf1 + chy_ub1); abscy2 = chy_bv + (chy_vf1 + chy_ub2)
                                abscy3 = chy_bv + (chy_vf2 + chy_ub1); abscy4 = chy_bv + (chy_vf2 + chy_ub2)
                            #:else
                                $:GPU_LOOP(parallelism='[seq]')
                                do r = 2, 6
                                    moms(r) = q_prim_vf(qbmm_bubmoms(q, r))%sf(id1, id2, id3)
                                end do
                                moms(1) = 1._wp
                                call s_chyqmom(moms, wght, abscX, abscY, qbmm_sgm_eps)
                            #:endif

                            #:if IS_POLY
                                $:GPU_LOOP(parallelism='[seq]')
                                do j = 1, 4
                                    wght_pb(j) = wght(j)*(qbmm_pb0(q) - qbmm_pv)
                                end do
                            #:else
                                #:if IS_KM
                                    ! Reciprocal-hoisted thermo chains: invariant and shared denominators
                                    ! become multiplies against hoisted reciprocals (see KM_HOIST_NAMES);
                                    ! per-node-varying denominators keep their divisions.
                                    inv_pb0 = 1._wp/qbmm_pb0(q)
                                    inv_r0 = 1._wp/qbmm_r0(q)
                                    inv_pe_t = 1._wp/qbmm_pe_t(q)
                                    #:for N in range(1, 5)
                                        chi_vw = 1._wp/(1._wp + ratio_rvg*(pb(id1, id2, id3, ${N}$, q)*inv_pv - 1._wp))
                                        x_vw = qbmm_m_g*chi_vw/(qbmm_m_v + (qbmm_m_g - qbmm_m_v)*chi_vw)
                                        k_mw = x_vw*qbmm_k_v(q)/(x_vw + (1._wp - x_vw)*qbmm_phi_vg) + (1._wp - x_vw)*qbmm_k_g(q) &
                                                             & /(x_vw*qbmm_phi_gv + 1._wp - x_vw)
                                        rho_mw = qbmm_pv/(chi_vw*qbmm_r_v*qbmm_tw)
                                        rhs_mv_pre${N}$ = -qbmm_re_trans_c(q)*((mv(id1, id2, id3, ${N}$, q)/(mv(id1, id2, id3, &
                                                                           & ${N}$, q) + qbmm_mass_g0(q))) - chi_vw)
                                        abscx_inv${N}$ = 1._wp/abscx${N}$
                                        rhs_mv_pre${N}$ = rho_mw*rhs_mv_pre${N}$*inv_pe_c/(1._wp - chi_vw)*abscx_inv${N}$
                                        grad_T = -qbmm_re_trans_t(q)*((pb(id1, id2, id3, ${N}$, &
                                                                  & q)*inv_pb0)*(abscx${N}$*inv_r0)*(abscx${N}$*inv_r0) &
                                                                  & *(abscx${N}$*inv_r0)*(qbmm_mass_g0(q) + qbmm_mass_v0(q)) &
                                                                  & /(qbmm_mass_g0(q) + mv(id1, id2, id3, ${N}$, q)) - 1._wp)
                                        ht${N}$ = qbmm_pb0(q)*k_mw*grad_T*inv_pe_t*abscx_inv${N}$
                                        wght_pb${N}$ = wght${N}$*pb(id1, id2, id3, ${N}$, q)
                                        wght_mv${N}$ = wght${N}$*rhs_mv_pre${N}$
                                        wght_ht${N}$ = wght${N}$*ht${N}$
                                    #:endfor
                                #:else
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do j = 1, 4
                                        chi_vw = 1._wp/(1._wp + qbmm_r_v/qbmm_r_g*(pb(id1, id2, id3, j, q)/qbmm_pv - 1._wp))
                                        x_vw = qbmm_m_g*chi_vw/(qbmm_m_v + (qbmm_m_g - qbmm_m_v)*chi_vw)
                                        k_mw = x_vw*qbmm_k_v(q)/(x_vw + (1._wp - x_vw)*qbmm_phi_vg) + (1._wp - x_vw)*qbmm_k_g(q) &
                                                             & /(x_vw*qbmm_phi_gv + 1._wp - x_vw)
                                        rho_mw = qbmm_pv/(chi_vw*qbmm_r_v*qbmm_tw)
                                        rhs_mv_pre(j) = -qbmm_re_trans_c(q)*((mv(id1, id2, id3, j, q)/(mv(id1, id2, id3, j, &
                                                   & q) + qbmm_mass_g0(q))) - chi_vw)
                                        rhs_mv_pre(j) = rho_mw*rhs_mv_pre(j)/qbmm_pe_c/(1._wp - chi_vw)/abscX(j)
                                        grad_T = -qbmm_re_trans_t(q)*((pb(id1, id2, id3, j, &
                                                                  & q)/qbmm_pb0(q))*(abscX(j)/qbmm_r0(q))**3*(qbmm_mass_g0(q) &
                                                                  & + qbmm_mass_v0(q))/(qbmm_mass_g0(q) + mv(id1, id2, id3, j, &
                                                                  & q)) - 1._wp)
                                        ht(j) = qbmm_pb0(q)*k_mw*grad_T/qbmm_pe_t(q)/abscX(j)
                                        wght_pb(j) = wght(j)*pb(id1, id2, id3, j, q)
                                        wght_mv(j) = wght(j)*rhs_mv_pre(j)
                                        wght_ht(j) = wght(j)*ht(j)
                                    end do
                                #:endif
                            #:endif

                            #:if IS_KM and not IS_POLY
                                #:for K in range(-3, 2)
                                    xpow_${KM_PXN(K)}$_n1 = f_integral_pow(abscx1, ${K}$)
                                    xpow_${KM_PXN(K)}$_n2 = f_integral_pow(abscx3, ${K}$)
                                #:endfor
                                #:for K in range(-1, 6)
                                    #:for NODE in range(1, 5)
                                        ypow_${KM_PXN(K)}$_n${NODE}$ = f_integral_pow(abscy${NODE}$, ${K}$)
                                    #:endfor
                                #:endfor
                            #:else
                                r0p(1) = 1._wp
                                $:GPU_LOOP(parallelism='[seq]')
                                do k = 2, qbmm_r0_n
                                    if (qbmm_r0_int(k) == qbmm_not_integral) then
                                        r0p(k) = qbmm_r0(q)**qbmm_r0_exp(k)
                                    else
                                        r0p(k) = f_integral_pow(qbmm_r0(q), qbmm_r0_int(k))
                                    end if
                                end do
                                xpe(1, 1) = 1._wp; xpe(1, 2) = 1._wp
                                $:GPU_LOOP(parallelism='[seq]')
                                do k = 2, qbmm_px_n(q)
                                    if (qbmm_px_int(k, q) == qbmm_not_integral) then
                                        xpe(k, 1) = abscX(1)**qbmm_px_exp(k, q)
                                        xpe(k, 2) = abscX(3)**qbmm_px_exp(k, q)
                                    else
                                        xpe(k, 1) = f_integral_pow(abscX(1), qbmm_px_int(k, q))
                                        xpe(k, 2) = f_integral_pow(abscX(3), qbmm_px_int(k, q))
                                    end if
                                end do
                                $:GPU_LOOP(parallelism='[seq]')
                                do node = 1, 4
                                    ypowers(1, node) = 1._wp
                                end do
                                $:GPU_LOOP(parallelism='[seq]')
                                do k = 2, qbmm_py_n(q)
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do node = 1, 4
                                        if (qbmm_py_int(k, q) == qbmm_not_integral) then
                                            ypowers(k, node) = abscY(node)**qbmm_py_exp(k, q)
                                        else
                                            ypowers(k, node) = f_integral_pow(abscY(node), qbmm_py_int(k, q))
                                        end if
                                    end do
                                end do
                            #:endif

                            #:if IS_KM and not IS_POLY
                                ! Moment sources with every (r, j) term a compile-time literal. The power
                                ! tables, coefficient algebra and quadrature weighting resolve at fypp time,
                                ! so the term loop indexes no local arrays and its live-term filters are
                                ! branches on constants.
                                #:for RI in KM_RI
                                    #:set I1 = RI[0]
                                    #:set I2 = RI[1]
                                    #:set R = KM_RI.index(RI) + 1
                                    momsum = 0._wp
                                    #:for J in range(1, nterms + 1)
                                        #:if (I2 == 0 and J == 12 and I1 != 0) or (I2 != 0 and not (13 <= J <= 15) &
                                            & and not (J == 12 and I1 == 0))
                                            #:set KX = qbmm_km_np_px_base[J - 1] + I1
                                            #:set KY = qbmm_km_np_py_base[J - 1] + I2
                                            #:set W = 'wght_pb' if J in KM_PB_QUAD else ('wght_mv' if J in KM_MV_QUAD &
                                                  & else ('wght_ht' if J in KM_HT_QUAD else 'wght'))
                                            #:set QP = [f'{W}{q_}*xpow_{KM_PXN(KX)}_n{1 if q_ in (1, 2) else 2}*ypow_{KM_PXN(KY)}_n{q_}' for q_ in (1, &
                                                & 2, 3, 4)]
                                            #:set GUARD = ('qbmm_re_active == 1 .and. qbmm_web_active == 1' if J in KM_RE_SET and J in KM_WEB_SET &
                                                  & else 'qbmm_re_active == 1' if J in KM_RE_SET else 'qbmm_web_active == 1' if J in KM_WEB_SET else None)
                                            #:if GUARD is not None
                                                if (${GUARD}$) then
                                                #:endif
                                                momsum = momsum + (${KM_COEFF(J, I1, I2)}$)*(${QP[0]}$ + ${QP[1]}$ + ${QP[2]}$ &
                                                                   & + ${QP[3]}$)
                                                #:if GUARD is not None
                                                end if
                                            #:endif
                                        #:endif
                                    #:endfor
                                    moms3d(${I1}$, ${I2}$, q)%sf(id1, id2, id3) = nbub*momsum
                                    msum${R}$ = momsum
                                #:endfor
                            #:else
                                $:GPU_LOOP(parallelism='[seq]')
                                do r = 1, qbmm_nout
                                    select case (r)
                                    case (1); i1 = 0; i2 = 0
                                    case (2); i1 = 1; i2 = 0
                                    case (3); i1 = 2; i2 = 0
                                    case (4); i1 = 0; i2 = 1
                                    case (5); i1 = 1; i2 = 1
                                    case default; i1 = 0; i2 = 2
                                    end select
                                    momsum = 0._wp
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do j = 1, nterms
                                        if (qbmm_live(r, j) == 0) cycle
                                        kx = qbmm_px_id(r, j, q)
                                        ky = qbmm_py_id(r, j, q)
                                        kr = qbmm_r0_id(r, j, q)
                                        xp1 = xpe(kx, 1); xp2 = xpe(kx, 2)
                                        yp1 = ypowers(ky, 1); yp2 = ypowers(ky, 2)
                                        yp3 = ypowers(ky, 3); yp4 = ypowers(ky, 4)
                                        #:if not IS_KM
                                            if (j == 3) then
                                                quad = wght_pb(1)*xp1*yp1 + wght_pb(2)*xp1*yp2 + wght_pb(3)*xp2*yp3 + wght_pb(4) &
                                                               & *xp2*yp4
                                            else
                                                quad = wght(1)*xp1*yp1 + wght(2)*xp1*yp2 + wght(3)*xp2*yp3 + wght(4)*xp2*yp4
                                            end if
                                        #:else
                                            if ((j >= 7 .and. j <= 11) .or. (j >= 22 .and. j <= 23) .or. j == 26) then
                                                quad = wght_pb(1)*xp1*yp1 + wght_pb(2)*xp1*yp2 + wght_pb(3)*xp2*yp3 + wght_pb(4) &
                                                               & *xp2*yp4
                                            else
                                                quad = wght(1)*xp1*yp1 + wght(2)*xp1*yp2 + wght(3)*xp2*yp3 + wght(4)*xp2*yp4
                                            end if
                                        #:endif
                                        coeff = f_coeff_${NAME}$(j, i1, i2, pres, rho, c, qbmm_re_active, qbmm_web_active, &
                                                                 & qbmm_re_inv, qbmm_web, qbmm_pv, qbmm_gam, qbmm_r_v, qbmm_tw)
                                        momsum = momsum + coeff*r0p(kr)*quad
                                    end do
                                    moms3d(i1, i2, q)%sf(id1, id2, id3) = nbub*momsum
                                    msum(r) = momsum
                                end do
                            #:endif

                            #:if not IS_POLY
                                #:if IS_KM
                                    #:for N in range(1, 5)
                                        drdt = msum2
                                        #:if N in (1, 2)
                                            drdt2 = -1._wp/(2._wp*sqrt(merge(mom4 - mom2**2._wp, qbmm_sgm_eps, &
                                                            & mom4 - mom2**2._wp > 0._wp)))
                                        #:else
                                            drdt2 = 1._wp/(2._wp*sqrt(merge(mom4 - mom2**2._wp, qbmm_sgm_eps, &
                                                           & mom4 - mom2**2._wp > 0._wp)))
                                        #:endif
                                        drdt2 = drdt2*(msum3 - 2._wp*mom2*msum2)
                                        drdt = drdt + drdt2
                                        rhs_pb(id1, id2, id3, ${N}$, q) = (-3._wp*qbmm_gam*drdt*abscx_inv${N}$)*pb(id1, id2, id3, &
                                               & ${N}$, q)
                                        rhs_pb(id1, id2, id3, ${N}$, q) = rhs_pb(id1, id2, id3, ${N}$, &
                                               & q) + (3._wp*qbmm_gam*abscx_inv${N}$)*rhs_mv_pre${N}$*qbmm_r_v*qbmm_tw
                                        rhs_pb(id1, id2, id3, ${N}$, q) = rhs_pb(id1, id2, id3, ${N}$, &
                                               & q) + (3._wp*qbmm_gam*abscx_inv${N}$)*ht${N}$
                                        rhs_mv(id1, id2, id3, ${N}$, q) = rhs_mv_pre${N}$*(4._wp*qbmm_pi*abscx${N}$**2._wp)
                                    #:endfor
                                #:else
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do j = 1, 4
                                        drdt = msum(2)
                                        drdt2 = merge(-1._wp, 1._wp, &
                                                      & j == 1 .or. j == 2)/(2._wp*sqrt(merge(moms(4) - moms(2)**2._wp, &
                                                      & qbmm_sgm_eps, moms(4) - moms(2)**2._wp > 0._wp)))
                                        drdt2 = drdt2*(msum(3) - 2._wp*moms(2)*msum(2))
                                        drdt = drdt + drdt2
                                        rhs_pb(id1, id2, id3, j, q) = (-3._wp*qbmm_gam*drdt/abscX(j))*pb(id1, id2, id3, j, q)
                                        rhs_pb(id1, id2, id3, j, q) = rhs_pb(id1, id2, id3, j, &
                                               & q) + (3._wp*qbmm_gam/abscX(j))*rhs_mv_pre(j)*qbmm_r_v*qbmm_tw
                                        rhs_pb(id1, id2, id3, j, q) = rhs_pb(id1, id2, id3, j, q) + (3._wp*qbmm_gam/abscX(j))*ht(j)
                                        rhs_mv(id1, id2, id3, j, q) = rhs_mv_pre(j)*(4._wp*qbmm_pi*abscX(j)**2._wp)
                                    end do
                                #:endif
                            #:endif

                            rv1 = 0._wp; rv2 = 0._wp; rv3 = 0._wp; rv4 = 0._wp
                            #:if IS_KM and not IS_POLY
                                #:for NODE in range(1, 5)
                                    rv1 = rv1 + wght${NODE}$*(abscx${NODE}$**3._wp)*(abscy${NODE}$**0._wp)
                                    rv2 = rv2 + wght${NODE}$*(abscx${NODE}$**2._wp)*(abscy${NODE}$**1._wp)
                                    rv3 = rv3 + wght${NODE}$*(abscx${NODE}$**3._wp)*(abscy${NODE}$**2._wp)
                                    rv4 = rv4 + wght${NODE}$*(abscx${NODE}$**2._wp)*(abscy${NODE}$**0._wp)
                                #:endfor
                            #:else
                                $:GPU_LOOP(parallelism='[seq]')
                                do node = 1, 4
                                    rv1 = rv1 + wght(node)*(abscX(node)**3._wp)*(abscY(node)**0._wp)
                                    rv2 = rv2 + wght(node)*(abscX(node)**2._wp)*(abscY(node)**1._wp)
                                    rv3 = rv3 + wght(node)*(abscX(node)**3._wp)*(abscY(node)**2._wp)
                                    rv4 = rv4 + wght(node)*(abscX(node)**2._wp)*(abscY(node)**0._wp)
                                end do
                            #:endif
                            fq1 = fq1 + qbmm_weight(q)*(qbmm_r0(q)**0._wp)*rv1
                            fq2 = fq2 + qbmm_weight(q)*(qbmm_r0(q)**0._wp)*rv2
                            fq3 = fq3 + qbmm_weight(q)*(qbmm_r0(q)**0._wp)*rv3
                            fq4 = fq4 + qbmm_weight(q)*(qbmm_r0(q)**0._wp)*rv4
                            if (abs(qbmm_gam - 1._wp) > 1.e-4_wp) then
                                rv5 = 0._wp
                                #:if IS_KM and not IS_POLY
                                    #:for NODE in range(1, 5)
                                        rv5 = rv5 + wght_pb${NODE}$*(abscx${NODE}$**3._wp)*(abscy${NODE}$**0._wp)
                                    #:endfor
                                #:else
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do node = 1, 4
                                        #:if IS_POLY
                                            rv5 = rv5 + wght_pb(node)*(abscX(node)**(3._wp*(1._wp - qbmm_gam)))*(abscY(node)**0._wp)
                                        #:else
                                            rv5 = rv5 + wght_pb(node)*(abscX(node)**3._wp)*(abscY(node)**0._wp)
                                        #:endif
                                    end do
                                #:endif
                                #:if IS_POLY
                                    fq5 = fq5 + qbmm_weight(q)*(qbmm_r0(q)**(3._wp*qbmm_gam))*rv5
                                #:else
                                    fq5 = fq5 + qbmm_weight(q)*(qbmm_r0(q)**0._wp)*rv5
                                #:endif
                            end if
                        end do

                        momsp(1)%sf(id1, id2, id3) = fq1
                        momsp(2)%sf(id1, id2, id3) = 4._wp*qbmm_pi*nbub*fq2
                        momsp(3)%sf(id1, id2, id3) = fq3
                        if (abs(qbmm_gam - 1._wp) <= 1.e-4_wp) then
                            momsp(4)%sf(id1, id2, id3) = 1._wp
                            #:if IS_POLY
                            else
                                momsp(4)%sf(id1, id2, id3) = fq5 + qbmm_pv*fq1 - 4._wp*qbmm_re_inv*fq2 - (2._wp/qbmm_web)*fq4
                            #:else
                            else
                                momsp(4)%sf(id1, id2, id3) = fq5 - 4._wp*qbmm_re_inv*fq2 - (2._wp/qbmm_web)*fq4
                            #:endif
                        end if
                    else
                        $:GPU_LOOP(parallelism='[seq]')
                        do q = 1, nb
                            $:GPU_LOOP(parallelism='[seq]')
                            do i1 = 0, 2
                                $:GPU_LOOP(parallelism='[seq]')
                                do i2 = 0, 2
                                    moms3d(i1, i2, q)%sf(id1, id2, id3) = 0._wp
                                end do
                            end do
                        end do
                        momsp(1)%sf(id1, id2, id3) = 0._wp
                        momsp(2)%sf(id1, id2, id3) = 0._wp
                        momsp(3)%sf(id1, id2, id3) = 0._wp
                        momsp(4)%sf(id1, id2, id3) = 0._wp
                    end if
                end do
                $:END_GPU_PARALLEL_LOOP()

            end subroutine s_mom_inv_${NAME}$
        #:enddef

        #:if nterms == 32
            #:if polytropic
                $:QBMM_KERNEL('km_poly', True, True)
            #:else
                $:QBMM_KERNEL('km_np', True, False)
            #:endif
        #:else
            #:if polytropic
                $:QBMM_KERNEL('rp_poly', False, True)
            #:else
                $:QBMM_KERNEL('rp_np', False, False)
            #:endif
        #:endif

        subroutine s_mom_inv_openacc_specialized(q_cons_vf, q_prim_vf, momsp, moms3d, pb, rhs_pb, mv, rhs_mv, ix_beg, ix_end, &
            & iy_beg, iy_end, iz_beg, iz_end, pb_l1, pb_l2, pb_l3)

            integer, intent(in)               :: ix_beg, ix_end, iy_beg, iy_end, iz_beg, iz_end, pb_l1, pb_l2, pb_l3
            type(scalar_field), intent(inout) :: q_cons_vf(:), q_prim_vf(:), momsp(:)
            type(scalar_field), intent(inout) :: moms3d(0:,0:,:)
            real(stp), intent(inout)          :: pb(pb_l1:,pb_l2:,pb_l3:,:,:), mv(pb_l1:,pb_l2:,pb_l3:,:,:)
            real(wp), intent(inout)           :: rhs_pb(pb_l1:,pb_l2:,pb_l3:,:,:), rhs_mv(pb_l1:,pb_l2:,pb_l3:,:,:)

            #:if nterms == 32
                #:if polytropic
                    call s_mom_inv_km_poly(q_cons_vf, q_prim_vf, momsp, moms3d, pb, rhs_pb, mv, rhs_mv, ix_beg, ix_end, iy_beg, &
                                           & iy_end, iz_beg, iz_end, pb_l1, pb_l2, pb_l3, qbmm_re_active, qbmm_web_active, &
                                           & qbmm_alf_idx, qbmm_pres_idx, qbmm_rho_idx, qbmm_nbub_idx, qbmm_r0_n, qbmm_gam, &
                                           & qbmm_pv, qbmm_r_v, qbmm_r_g, qbmm_m_v, qbmm_m_g, qbmm_phi_vg, qbmm_phi_gv, qbmm_tw, &
                                           & qbmm_pe_c, qbmm_re_inv, qbmm_web, qbmm_sgm_eps, qbmm_small_alf, qbmm_pi, &
                                           & qbmm_gammas1, qbmm_pi_infs1)
                #:else
                    call s_mom_inv_km_np(q_cons_vf, q_prim_vf, momsp, moms3d, pb, rhs_pb, mv, rhs_mv, ix_beg, ix_end, iy_beg, &
                                         & iy_end, iz_beg, iz_end, pb_l1, pb_l2, pb_l3, qbmm_re_active, qbmm_web_active, &
                                         & qbmm_alf_idx, qbmm_pres_idx, qbmm_rho_idx, qbmm_nbub_idx, qbmm_r0_n, qbmm_gam, &
                                         & qbmm_pv, qbmm_r_v, qbmm_r_g, qbmm_m_v, qbmm_m_g, qbmm_phi_vg, qbmm_phi_gv, qbmm_tw, &
                                         & qbmm_pe_c, qbmm_re_inv, qbmm_web, qbmm_sgm_eps, qbmm_small_alf, qbmm_pi, qbmm_gammas1, &
                                         & qbmm_pi_infs1)
                #:endif
            #:else
                #:if polytropic
                    call s_mom_inv_rp_poly(q_cons_vf, q_prim_vf, momsp, moms3d, pb, rhs_pb, mv, rhs_mv, ix_beg, ix_end, iy_beg, &
                                           & iy_end, iz_beg, iz_end, pb_l1, pb_l2, pb_l3, qbmm_re_active, qbmm_web_active, &
                                           & qbmm_alf_idx, qbmm_pres_idx, qbmm_rho_idx, qbmm_nbub_idx, qbmm_r0_n, qbmm_gam, &
                                           & qbmm_pv, qbmm_r_v, qbmm_r_g, qbmm_m_v, qbmm_m_g, qbmm_phi_vg, qbmm_phi_gv, qbmm_tw, &
                                           & qbmm_pe_c, qbmm_re_inv, qbmm_web, qbmm_sgm_eps, qbmm_small_alf, qbmm_pi, &
                                           & qbmm_gammas1, qbmm_pi_infs1)
                #:else
                    call s_mom_inv_rp_np(q_cons_vf, q_prim_vf, momsp, moms3d, pb, rhs_pb, mv, rhs_mv, ix_beg, ix_end, iy_beg, &
                                         & iy_end, iz_beg, iz_end, pb_l1, pb_l2, pb_l3, qbmm_re_active, qbmm_web_active, &
                                         & qbmm_alf_idx, qbmm_pres_idx, qbmm_rho_idx, qbmm_nbub_idx, qbmm_r0_n, qbmm_gam, &
                                         & qbmm_pv, qbmm_r_v, qbmm_r_g, qbmm_m_v, qbmm_m_g, qbmm_phi_vg, qbmm_phi_gv, qbmm_tw, &
                                         & qbmm_pe_c, qbmm_re_inv, qbmm_web, qbmm_sgm_eps, qbmm_small_alf, qbmm_pi, qbmm_gammas1, &
                                         & qbmm_pi_infs1)
                #:endif
            #:endif

        end subroutine s_mom_inv_openacc_specialized
    #:endif
#endif
end module m_qbmm_openacc
