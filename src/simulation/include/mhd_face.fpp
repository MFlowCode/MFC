! Single-fluid magnetic face closure shared by the Riemann solvers and the fused column.
#:if not defined('hllc_mhd')
    #:set hllc_mhd = False
#:endif
#:set MHD_FACE_SCALARS = ['sl', 'sr', 'sm', 'sp', 'contact', 'pstar', 'ssl', 'ssr', 'denom', 'sgn', 'vd', 'wd', 'byd', 'bzd', 'bn']
#:set MHD_FACE_SIDE_SCALARS = ['rho', 'pres', 'gamma', 'pi', 'qv', 'v2', 'b2', 'vb', 'pm', 'energy', 'h', 'ga', 'c', 'term', &
    & 'disc', 'pt', 'rstar', 'estar', 'root', 'edouble', 'dstar', 'vs', 'ws', 'bys', 'bzs']
#:set MHD_FACE_PRIVATE = [f'mf_{name}' for name in MHD_FACE_SCALARS] + [f'mf_{name}_{side}' for side in ['l', &
    & 'r'] for name in MHD_FACE_SIDE_SCALARS + ['v', 'b', 'mom', 'b4', 'u', 'f', 'us', 'fs', 'ud']]

#:def mhd_face_declarations()
    real(wp) :: ${', '.join('mf_' + name for name in MHD_FACE_SCALARS)}$
    #:for SIDE in ['l', 'r']
        real(wp) :: ${', '.join(f'mf_{name}_{SIDE}' for name in MHD_FACE_SIDE_SCALARS)}$
        real(wp) :: mf_v_${SIDE}$(3), mf_b_${SIDE}$(3), mf_mom_${SIDE}$(3), mf_b4_${SIDE}$(3)
        real(wp) :: mf_u_${SIDE}$(7), mf_f_${SIDE}$(7), mf_us_${SIDE}$(7), mf_fs_${SIDE}$(7), mf_ud_${SIDE}$(7)
    #:endfor
#:enddef

#:def mhd_face(DIM, DIR, SOLVER, Q, FLUX, GAMMA, PI, QV, BX, CLEAN, MIXTURE_ERROR)
    #! Local vectors use the common solver's normal-first permutation.
    #:set PERM = [DIR] + [v for v in range(1, 4) if v != DIR]
    #:set BROW = lambda v: 5 + v if DIM == 1 else 6 + v
    #:for SIDE in ['l', 'r']
        mf_rho_${SIDE}$ = ${Q(SIDE, 1)}$
        mf_pres_${SIDE}$ = ${Q(SIDE, 5)}$
        mf_gamma_${SIDE}$ = ${Q(SIDE, 6)}$*${GAMMA}$
        mf_pi_${SIDE}$ = ${Q(SIDE, 6)}$*${PI}$
        mf_qv_${SIDE}$ = mf_rho_${SIDE}$*${QV}$
        #:for V, ORIGINAL in enumerate(PERM, 1)
            mf_v_${SIDE}$(${V}$) = ${Q(SIDE, 1 + ORIGINAL)}$
            #:if SOLVER == 4 and DIM > 1 and V == 1
                #! The Riemann fan has a single normal field, the face average
                mf_b_${SIDE}$(1) = 0.5_wp*(${Q('l', BROW(ORIGINAL))}$ + ${Q('r', BROW(ORIGINAL))}$)
            #:else
                mf_b_${SIDE}$(${V}$) = ${BX if DIM == 1 and ORIGINAL == 1 else Q(SIDE, BROW(ORIGINAL))}$
            #:endif
        #:endfor
        mf_v2_${SIDE}$ = mf_v_${SIDE}$(1)**2 + mf_v_${SIDE}$(2)**2 + mf_v_${SIDE}$(3)**2
        mf_b2_${SIDE}$ = mf_b_${SIDE}$(1)**2 + mf_b_${SIDE}$(2)**2 + mf_b_${SIDE}$(3)**2
        mf_vb_${SIDE}$ = mf_v_${SIDE}$(1)*mf_b_${SIDE}$(1) + mf_v_${SIDE}$(2)*mf_b_${SIDE}$(2) + mf_v_${SIDE}$(3)*mf_b_${SIDE}$(3)
        #:if relativity
            mf_ga_${SIDE}$ = 1._wp/sqrt(1._wp - mf_v2_${SIDE}$)
            mf_pm_${SIDE}$ = 0.5_wp*(mf_b2_${SIDE}$/mf_ga_${SIDE}$**2 + mf_vb_${SIDE}$**2)
            mf_h_${SIDE}$ = 1._wp + (mf_gamma_${SIDE}$ + 1._wp)*mf_pres_${SIDE}$/mf_rho_${SIDE}$
            #:for V in range(1, 4)
                mf_b4_${SIDE}$(${V}$) = mf_b_${SIDE}$(${V}$)/mf_ga_${SIDE}$ + mf_ga_${SIDE}$*mf_v_${SIDE}$(${V}$)*mf_vb_${SIDE}$
                mf_mom_${SIDE}$(${V}$) = (mf_rho_${SIDE}$*mf_h_${SIDE}$*mf_ga_${SIDE}$**2 + mf_b2_${SIDE}$)*mf_v_${SIDE}$(${V}$) &
                                & - mf_vb_${SIDE}$*mf_b_${SIDE}$(${V}$)
            #:endfor
            mf_energy_${SIDE}$ = mf_rho_${SIDE}$*mf_h_${SIDE}$*mf_ga_${SIDE}$**2 - mf_pres_${SIDE}$ + 0.5_wp*(mf_b2_${SIDE}$ &
                & + mf_v2_${SIDE}$*mf_b2_${SIDE}$ - mf_vb_${SIDE}$**2) - mf_rho_${SIDE}$*mf_ga_${SIDE}$
            mf_c_${SIDE}$ = sqrt((1._wp + 1._wp/mf_gamma_${SIDE}$)*mf_pres_${SIDE}$/mf_rho_${SIDE}$/mf_h_${SIDE}$)
            mf_term_${SIDE}$ = (mf_c_${SIDE}$**2*(mf_b_${SIDE}$(1)**2 + mf_rho_${SIDE}$*mf_h_${SIDE}$) + mf_b2_${SIDE}$) &
                                & /(mf_rho_${SIDE}$*mf_h_${SIDE}$ + mf_b2_${SIDE}$)
            mf_disc_${SIDE}$ = mf_term_${SIDE}$**2 - 4._wp*mf_c_${SIDE}$**2*mf_b_${SIDE}$(1)**2/(mf_rho_${SIDE}$*mf_h_${SIDE}$ &
                & + mf_b2_${SIDE}$)
        #:else
            mf_pm_${SIDE}$ = 0.5_wp*mf_b2_${SIDE}$
            mf_energy_${SIDE}$ = mf_gamma_${SIDE}$*mf_pres_${SIDE}$ + mf_pi_${SIDE}$ + 0.5_wp*mf_rho_${SIDE}$*mf_v2_${SIDE}$ &
                & + mf_qv_${SIDE}$ + mf_pm_${SIDE}$
            mf_h_${SIDE}$ = (mf_energy_${SIDE}$ + mf_pres_${SIDE}$ - mf_pm_${SIDE}$)/mf_rho_${SIDE}$
            mf_c_${SIDE}$ = (mf_h_${SIDE}$ - 0.5_wp*mf_v2_${SIDE}$ - mf_qv_${SIDE}$/mf_rho_${SIDE}$)/mf_gamma_${SIDE}$
            if (${MIXTURE_ERROR}$ .and. mf_c_${SIDE}$ < 0._wp) then
                mf_c_${SIDE}$ = 100._wp*sgm_eps
            else
                mf_c_${SIDE}$ = sqrt(mf_c_${SIDE}$)
            end if
            mf_term_${SIDE}$ = mf_c_${SIDE}$**2 + mf_b2_${SIDE}$/mf_rho_${SIDE}$
            mf_disc_${SIDE}$ = mf_term_${SIDE}$**2 - 4._wp*mf_c_${SIDE}$**2*(mf_b_${SIDE}$(1)**2/mf_rho_${SIDE}$)
        #:endif
        mf_c_${SIDE}$ = sqrt(0.5_wp*(mf_term_${SIDE}$ + sqrt(mf_disc_${SIDE}$)))
        mf_pt_${SIDE}$ = mf_pres_${SIDE}$ + mf_pm_${SIDE}$
        mf_u_${SIDE}$(1) = mf_rho_${SIDE}$${f'*mf_ga_{SIDE}' if relativity else ''}$
        mf_f_${SIDE}$(1) = mf_u_${SIDE}$(1)*mf_v_${SIDE}$(1)
        #:for V in range(1, 4)
            #:if relativity
                mf_u_${SIDE}$(${V + 1}$) = mf_mom_${SIDE}$(${V}$)
                mf_f_${SIDE}$(${V + 1}$) = mf_mom_${SIDE}$(${V}$)*mf_v_${SIDE}$(1) - mf_b4_${SIDE}$(${V}$) &
                              & /mf_ga_${SIDE}$*mf_b_${SIDE}$(1)${f' + mf_pt_{SIDE}' if V == 1 else ''}$
            #:else
                mf_u_${SIDE}$(${V + 1}$) = mf_rho_${SIDE}$*mf_v_${SIDE}$(${V}$)
                mf_f_${SIDE}$(${V + 1}$) = mf_rho_${SIDE}$*mf_v_${SIDE}$(${V}$)*mf_v_${SIDE}$(1) - mf_b_${SIDE}$(${V}$) &
                              & *mf_b_${SIDE}$(1)${f' + mf_pt_{SIDE}' if V == 1 else ''}$
            #:endif
        #:endfor
        #:for V in [2, 3]
            mf_u_${SIDE}$(${V + 3}$) = mf_b_${SIDE}$(${V}$)
            mf_f_${SIDE}$(${V + 3}$) = mf_v_${SIDE}$(1)*mf_b_${SIDE}$(${V}$) - mf_v_${SIDE}$(${V}$)*mf_b_${SIDE}$(1)
        #:endfor
        mf_u_${SIDE}$(7) = mf_energy_${SIDE}$
        #:if relativity
            mf_f_${SIDE}$(7) = mf_mom_${SIDE}$(1) - mf_ga_${SIDE}$*mf_rho_${SIDE}$*mf_v_${SIDE}$(1)
        #:else
            mf_f_${SIDE}$(7) = (mf_energy_${SIDE}$ + mf_pt_${SIDE}$)*mf_v_${SIDE}$(1) - mf_b_${SIDE}$(1)*mf_vb_${SIDE}$
        #:endif
    #:endfor
    mf_sl = min(mf_v_l(1) - mf_c_l, mf_v_r(1) - mf_c_r)
    mf_sr = max(mf_v_r(1) + mf_c_r, mf_v_l(1) + mf_c_l)
    #:if hyper_cleaning
        mf_sl = min(mf_sl, -${CLEAN}$)
        mf_sr = max(mf_sr, ${CLEAN}$)
    #:endif
    #:set OUTPUTS = [1] + [1 + v for v in PERM] + [BROW(v) for v in PERM[1:]] + [5]
    #:if SOLVER == 1
        mf_sm = min(0._wp, mf_sl); mf_sp = max(0._wp, mf_sr)
        #:for ROW, OUT in enumerate(OUTPUTS, 1)
            ${FLUX(OUT)}$ = (mf_sm*mf_f_r(${ROW}$) - mf_sp*mf_f_l(${ROW}$) + mf_sm*mf_sp*(mf_u_l(${ROW}$) - mf_u_r(${ROW}$))) &
                   & /(mf_sm - mf_sp)
        #:endfor
        ${FLUX(6)}$ = (${Q('l', 6)}$ - ${Q('r', 6)}$)*mf_sm*mf_sp/(mf_sm - mf_sp)
        #:if DIM > 1
            #:if hyper_cleaning
                ${FLUX(BROW(DIR))}$ = (mf_sm*mf_sp*(mf_b_l(1) - mf_b_r(1)) + mf_sm*${Q('r', 10)}$ - mf_sp*${Q('l', 10)}$)/(mf_sm &
                       & - mf_sp)
                ${FLUX(10)}$ = (${CLEAN}$**2*(mf_sm*mf_b_r(1) - mf_sp*mf_b_l(1)) + mf_sm*mf_sp*(${Q('l', 10)}$ - ${Q('r', 10)}$)) &
                       & /(mf_sm - mf_sp)
            #:else
                ${FLUX(BROW(DIR))}$ = 0._wp
            #:endif
        #:endif
    #:else
        #! Miyoshi and Kusano (2005) intermediate states, as in the reference HLLD solver
        mf_contact = ((mf_sr - mf_v_r(1))*mf_rho_r*mf_v_r(1) - (mf_sl - mf_v_l(1))*mf_rho_l*mf_v_l(1) - mf_pt_r + mf_pt_l) &
                      & /((mf_sr - mf_v_r(1))*mf_rho_r - (mf_sl - mf_v_l(1))*mf_rho_l)
        mf_bn = mf_b_l(1)
        mf_pstar = mf_pt_l + mf_rho_l*(mf_sl - mf_v_l(1))*(mf_contact - mf_v_l(1))
        #:for SIDE, SPEED in [('l', 'sl'), ('r', 'sr')]
            mf_rstar_${SIDE}$ = mf_rho_${SIDE}$*(mf_${SPEED}$ - mf_v_${SIDE}$(1))/(mf_${SPEED}$ - mf_contact)
            mf_dstar_${SIDE}$ = mf_rho_${SIDE}$*(mf_${SPEED}$ - mf_v_${SIDE}$(1))*(mf_${SPEED}$ - mf_contact) - mf_bn**2
            mf_vs_${SIDE}$ = mf_v_${SIDE}$(2); mf_ws_${SIDE}$ = mf_v_${SIDE}$(3)
            mf_bys_${SIDE}$ = mf_b_${SIDE}$(2); mf_bzs_${SIDE}$ = mf_b_${SIDE}$(3)
            if (abs(mf_dstar_${SIDE}$) > sgm_eps*max(mf_rho_${SIDE}$*(mf_${SPEED}$ - mf_v_${SIDE}$(1))**2, mf_bn**2)) then
                mf_vs_${SIDE}$ = mf_v_${SIDE}$(2) - mf_bn*mf_b_${SIDE}$(2)*(mf_contact - mf_v_${SIDE}$(1))/mf_dstar_${SIDE}$
                mf_ws_${SIDE}$ = mf_v_${SIDE}$(3) - mf_bn*mf_b_${SIDE}$(3)*(mf_contact - mf_v_${SIDE}$(1))/mf_dstar_${SIDE}$
                mf_bys_${SIDE}$ = mf_b_${SIDE}$(2)*(mf_rho_${SIDE}$*(mf_${SPEED}$ - mf_v_${SIDE}$(1))**2 - mf_bn**2) &
                                                & /mf_dstar_${SIDE}$
                mf_bzs_${SIDE}$ = mf_b_${SIDE}$(3)*(mf_rho_${SIDE}$*(mf_${SPEED}$ - mf_v_${SIDE}$(1))**2 - mf_bn**2) &
                                                & /mf_dstar_${SIDE}$
            end if
            mf_estar_${SIDE}$ = ((mf_${SPEED}$ - mf_v_${SIDE}$(1))*mf_energy_${SIDE}$ - mf_pt_${SIDE}$*mf_v_${SIDE}$(1) &
                                 & + mf_pstar*mf_contact + mf_bn*(mf_v_${SIDE}$(1)*mf_bn + mf_v_${SIDE}$(2)*mf_b_${SIDE}$(2) &
                                 & + mf_v_${SIDE}$(3)*mf_b_${SIDE}$(3) - mf_contact*mf_bn - mf_vs_${SIDE}$*mf_bys_${SIDE}$ &
                                 & - mf_ws_${SIDE}$*mf_bzs_${SIDE}$))/(mf_${SPEED}$ - mf_contact)
            mf_root_${SIDE}$ = sqrt(mf_rstar_${SIDE}$)
            mf_us_${SIDE}$(1) = mf_rstar_${SIDE}$
            mf_us_${SIDE}$(2) = mf_rstar_${SIDE}$*mf_contact
            mf_us_${SIDE}$(3) = mf_rstar_${SIDE}$*mf_vs_${SIDE}$
            mf_us_${SIDE}$(4) = mf_rstar_${SIDE}$*mf_ws_${SIDE}$
            mf_us_${SIDE}$(5) = mf_bys_${SIDE}$; mf_us_${SIDE}$(6) = mf_bzs_${SIDE}$
            mf_us_${SIDE}$(7) = mf_estar_${SIDE}$
            #:for ROW in range(1, 8)
                mf_fs_${SIDE}$(${ROW}$) = mf_f_${SIDE}$(${ROW}$) + mf_${SPEED}$*(mf_us_${SIDE}$(${ROW}$) - mf_u_${SIDE}$(${ROW}$))
            #:endfor
        #:endfor
        mf_ssl = mf_contact - abs(mf_bn)/mf_root_l
        mf_ssr = mf_contact + abs(mf_bn)/mf_root_r
        mf_denom = mf_root_l + mf_root_r; mf_sgn = sign(1._wp, mf_bn)
        mf_vd = (mf_root_l*mf_vs_l + mf_root_r*mf_vs_r + (mf_bys_r - mf_bys_l)*mf_sgn)/mf_denom
        mf_wd = (mf_root_l*mf_ws_l + mf_root_r*mf_ws_r + (mf_bzs_r - mf_bzs_l)*mf_sgn)/mf_denom
        mf_byd = (mf_root_l*mf_bys_r + mf_root_r*mf_bys_l + mf_root_l*mf_root_r*(mf_vs_r - mf_vs_l)*mf_sgn)/mf_denom
        mf_bzd = (mf_root_l*mf_bzs_r + mf_root_r*mf_bzs_l + mf_root_l*mf_root_r*(mf_ws_r - mf_ws_l)*mf_sgn)/mf_denom
        #:for SIDE, SIGN in [('l', '-'), ('r', '+')]
            mf_edouble_${SIDE}$ = mf_estar_${SIDE}$ ${SIGN}$ mf_root_${SIDE}$*((mf_vs_${SIDE}$*mf_bys_${SIDE}$ &
                & + mf_ws_${SIDE}$*mf_bzs_${SIDE}$) - (mf_vd*mf_byd + mf_wd*mf_bzd))*mf_sgn
        #:endfor
        #:for SIDE in ['l', 'r']
            mf_ud_${SIDE}$(1) = mf_rstar_${SIDE}$; mf_ud_${SIDE}$(2) = mf_rstar_${SIDE}$*mf_contact
            mf_ud_${SIDE}$(3) = mf_rstar_${SIDE}$*mf_vd; mf_ud_${SIDE}$(4) = mf_rstar_${SIDE}$*mf_wd
            mf_ud_${SIDE}$(5) = mf_byd; mf_ud_${SIDE}$(6) = mf_bzd; mf_ud_${SIDE}$(7) = mf_edouble_${SIDE}$
        #:endfor
        #:for ROW, OUT in enumerate(OUTPUTS, 1)
            if (0._wp <= mf_sl) then
                ${FLUX(OUT)}$ = mf_f_l(${ROW}$)
            else if (0._wp <= mf_ssl) then
                ${FLUX(OUT)}$ = mf_fs_l(${ROW}$)
            else if (0._wp <= mf_contact) then
                ${FLUX(OUT)}$ = mf_fs_l(${ROW}$) + mf_ssl*(mf_ud_l(${ROW}$) - mf_us_l(${ROW}$))
            else if (0._wp <= mf_ssr) then
                ${FLUX(OUT)}$ = mf_fs_r(${ROW}$) + mf_ssr*(mf_ud_r(${ROW}$) - mf_us_r(${ROW}$))
            else if (0._wp <= mf_sr) then
                ${FLUX(OUT)}$ = mf_fs_r(${ROW}$)
            else
                ${FLUX(OUT)}$ = mf_f_r(${ROW}$)
            end if
        #:endfor
        ${FLUX(6)}$ = 0._wp
        #:if DIM > 1
            ${FLUX(BROW(DIR))}$ = 0._wp
        #:endif
    #:endif
#:enddef
