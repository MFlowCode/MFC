!>
!! @file
!! @brief HyQMOM/CHyQMOM moment-inversion bodies shared by m_qbmm and m_qbmm_openacc

#! Numerical bodies shared by the m_qbmm and m_qbmm_openacc wrappers, parameterized on the
#! variance-floor epsilon.

#:def QBMM_HYQMOM_BODY(EPS)
    real(wp) :: bu, d2, c2

    bu = fmom(2)/fmom(1)
    d2 = fmom(3)/fmom(1)
    c2 = d2 - bu**2._wp
    frho(1) = fmom(1)/2._wp
    frho(2) = fmom(1)/2._wp
    c2 = maxval((/c2, ${EPS}$/))
    fup(1) = bu - sqrt(c2)
    fup(2) = bu + sqrt(c2)
#:enddef

#:def QBMM_CHYQMOM_BODY(EPS)
    real(wp), dimension(0:2,0:2) :: moms
    real(wp), dimension(3)       :: M1, M3
    real(wp), dimension(2)       :: myrho, myrho3, up, up3, Vf
    real(wp)                     :: bu, bv, d20, d11, d_02, c20, c11, c02
    real(wp)                     :: mu2, vp21, vp22, rho21, rho22

    moms(0, 0) = momin(1); moms(1, 0) = momin(2); moms(0, 1) = momin(3)
    moms(2, 0) = momin(4); moms(1, 1) = momin(5); moms(0, 2) = momin(6)
    bu = moms(1, 0)/moms(0, 0); bv = moms(0, 1)/moms(0, 0)
    d20 = moms(2, 0)/moms(0, 0); d11 = moms(1, 1)/moms(0, 0); d_02 = moms(0, 2)/moms(0, 0)
    c20 = d20 - bu**2._wp; c11 = d11 - bu*bv; c02 = d_02 - bv**2._wp
    M1 = (/1._wp, 0._wp, c20/)
    call s_hyqmom(myrho, up, M1, ${EPS}$)
    Vf = c11*up/c20
    mu2 = max(0._wp, c02 - sum(myrho*(Vf**2._wp)))
    M3 = (/1._wp, 0._wp, mu2/)
    call s_hyqmom(myrho3, up3, M3, ${EPS}$)
    vp21 = up3(1); vp22 = up3(2); rho21 = myrho3(1); rho22 = myrho3(2)
    wght = moms(0, 0)*[myrho(1)*rho21, myrho(1)*rho22, myrho(2)*rho21, myrho(2)*rho22]
    abscX = bu + [up(1), up(1), up(2), up(2)]
    abscY = bv + [Vf(1) + vp21, Vf(1) + vp22, Vf(2) + vp21, Vf(2) + vp22]
#:enddef
