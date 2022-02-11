!>
!! @file m_bubbles.f90
!! @brief Contains module m_bubbles
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief This module is used to compute the ensemble-averaged bubble dynamic variables
module m_bubbles

    ! Dependencies =============================================================

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    ! ==========================================================================

    implicit none

    real(kind(0.d0)) :: chi_vw  !< Bubble wall properties (Ando 2010)
    real(kind(0.d0)) :: k_mw    !< Bubble wall properties (Ando 2010)
    real(kind(0.d0)) :: rho_mw  !< Bubble wall properties (Ando 2010)
!$acc declare create(chi_vw, k_mw, rho_mw)

contains

    !>  The purpose of this procedure is to compute the source terms
        !!      that are needed for the bubble modeling
        !!  @param idir Dimension splitting index
        !!  @param q_prim_vf Primitive variables
        !!  @param q_cons_vf Conservative variables
        !!  @param mydivu Divergence of velocity
        !!  @param bub_adv_src Advection equation source due to bubble compression/expansion
        !!  @param bub_r_src   Bubble radius equation source
        !!  @param bub_v_src   Bubble velocity equation source
        !!  @param bub_p_src   Bubble pressure equation source
        !!  @param bub_m_src   Bubble mass equation source
    subroutine s_compute_bubble_source(idir, q_prim_vf, q_cons_vf, mydivu, &
                                       bub_adv_src, bub_r_src, bub_v_src, bub_p_src, bub_m_src)

        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf, q_cons_vf
        type(scalar_field), intent(IN) :: mydivu
        integer, intent(IN) :: idir

        real(kind(0d0)), dimension(0:m, 0:n, 0:p), intent(INOUT) ::  bub_adv_src
        real(kind(0d0)), dimension(0:m, 0:n, 0:p, 1:nb), intent(INOUT) ::  bub_r_src, &
                                                                          bub_v_src, &
                                                                          bub_p_src, &
                                                                          bub_m_src

        real(kind(0d0)), dimension(0:m, 0:n, 0:p) :: nbub !< Bubble number density

        real(kind(0d0)) ::  tmp1, tmp2, tmp3, tmp4, &
                           c_gas, c_liquid, &
                           Cpbw, Cpinf, Cpinf_dot, &
                           myH, myHdot, rddot, alf_gas

        real(kind(0d0))   :: pb, mv, vflux, pldot, pbdot

        real(kind(0d0)) :: n_tait, B_tait

        real(kind(0d0)), dimension(nb)  :: Rtmp, Vtmp
        real(kind(0d0))   :: myR, myV, alf, myP, myRho, R2Vav

        real(kind(0d0)), dimension(2) :: Re !< Reynolds number

        integer :: j, k, l, q, s !< Loop variables
        integer :: ndirs  !< Number of coordinate directions

        ndirs = 1; if (n > 0) ndirs = 2; if (p > 0) ndirs = 3

        if (idir == ndirs) then
            bub_adv_src = 0.d0; bub_r_src = 0.d0; bub_v_src = 0.d0
            bub_p_src = 0.d0; bub_m_src = 0.d0

            ! advection source
            do j = 0, m; do k = 0, n; do l = 0, p
                    ! = 3 \alpha \bar{R^2 V} / \bar{R^3} = 4 pi nbub \bar{R^2 V}
                    do q = 1, nb
                        Rtmp(q) = q_prim_vf(bub_idx%rs(q))%sf(j, k, l)
                        Vtmp(q) = q_prim_vf(bub_idx%vs(q))%sf(j, k, l)
                    end do

                    ! Computes n_bub number bubble density
                    call s_comp_n_from_prim(q_prim_vf(alf_idx)%sf(j, k, l), &
                                            Rtmp, nbub(j, k, l))

                    call s_quad((Rtmp**2.d0)*Vtmp, R2Vav)
                    bub_adv_src(j, k, l) = 4.d0*pi*nbub(j, k, l)*R2Vav
                end do; end do; end do

            ! bubble radius and radial velocity source
            do q = 1, nb; do j = 0, m; do k = 0, n; do l = 0, p
                    bub_r_src( j, k, l, q) = q_cons_vf(bub_idx%vs(q))%sf(j, k, l)

                    call s_convert_to_mixture_variables(q_cons_vf, myRho, n_tait, B_tait, Re, j, k, l)

                    n_tait = 1.d0/n_tait + 1.d0 !make this the usual little 'gamma'

                    myRho = q_prim_vf(1)%sf(j, k, l)
                    myP = q_prim_vf(E_idx)%sf(j, k, l)
                    alf = q_prim_vf(alf_idx)%sf(j, k, l)
                    myR = q_prim_vf(bub_idx%rs(q))%sf(j, k, l)
                    myV = q_prim_vf(bub_idx%vs(q))%sf(j, k, l)

                    if (.not. polytropic) then
                        pb = q_prim_vf(bub_idx%ps(q))%sf(j, k, l)
                        mv = q_prim_vf(bub_idx%ms(q))%sf(j, k, l)
                        call s_bwproperty(pb, q)
                        vflux = f_vflux(myR, myV, mv, q)
                        pbdot = f_bpres_dot(vflux, myR, myV, pb, mv, q)

                        bub_p_src( j, k, l, q) = nbub(j, k, l)*pbdot
                        bub_m_src( j, k, l, q) = nbub(j, k, l)*vflux*4.d0*pi*(myR**2.d0)
                    else
                        pb = 0d0; mv = 0d0; vflux = 0d0; pbdot = 0d0
                    end if

                    if (bubble_model == 1) then
                        ! Gilmore bubbles
                        Cpinf = myP - pref
                        Cpbw = f_cpbw(R0(q), myR, myV, pb)
                        myH = f_H(Cpbw, Cpinf, n_tait, B_tait)
                        c_gas = f_cgas(Cpinf, n_tait, B_tait, myH)
                        Cpinf_dot = f_cpinfdot(myRho, myP, alf, n_tait, B_tait, bub_adv_src(j, k, l), mydivu%sf(j, k, l))
                        myHdot = f_Hdot(Cpbw, Cpinf, Cpinf_dot, n_tait, B_tait, myR, myV, R0(q), pbdot)
                        rddot = f_rddot(Cpbw, myR, myV, myH, myHdot, c_gas, n_tait, B_tait)
                    else if (bubble_model == 2) then
                        ! Keller-Miksis bubbles
                        Cpinf = myP
                        Cpbw = f_cpbw_KM(R0(q), myR, myV, pb)
                        ! c_gas = dsqrt( n_tait*(Cpbw+B_tait) / myRho)
                        c_liquid = DSQRT(n_tait*(myP + B_tait)/(myRho*(1.d0 - alf)))
                        rddot = f_rddot_KM(pbdot, Cpinf, Cpbw, myRho, myR, myV, R0(q), c_liquid)
                    else if (bubble_model == 3) then
                        ! Rayleigh-Plesset bubbles
                        Cpbw = f_cpbw_KM(R0(q), myR, myV, pb)
                        rddot = f_rddot_RP(myP, myRho, myR, myV, R0(q), Cpbw)
                    end if

                    bub_v_src( j, k, l, q) = nbub(j, k, l)*rddot

                    if (alf < 1.d-11) then
                        bub_adv_src(j, k, l) = 0d0
                        bub_r_src( j, k, l, q) = 0d0
                        bub_v_src( j, k, l, q) = 0d0
                        if (.not. polytropic) then
                            bub_p_src( j, k, l, q) = 0d0
                            bub_m_src( j, k, l, q) = 0d0
                        end if
                    end if

                end do; end do; end do; end do
        end if

    end subroutine s_compute_bubble_source

    !>  Function that computes that bubble wall pressure for Gilmore bubbles
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fpb Internal bubble pressure
    function f_cpbw(fR0, fR, fV, fpb)
!$acc routine seq
        real(kind(0d0)), intent(IN) :: fR0, fR, fV, fpb

        real(kind(0d0))             :: f_cpbw

        if (polytropic) then
            f_cpbw = (Ca + 2.d0/Web/fR0)*((fR0/fR)**(3.d0*gam)) - Ca - 4.d0*Re_inv*fV/fR - 2.d0/(fR*Web)
        else
            f_cpbw = fpb - 1.d0 - 4.d0*Re_inv*fV/fR - 2.d0/(fR*Web)
        end if

    end function f_cpbw

    !>  Function that computes the bubble enthalpy
        !!  @param fCpbw Bubble wall pressure
        !!  @param fCpinf Driving bubble pressure
        !!  @param fntait Tait EOS parameter
        !!  @param fBtait Tait EOS parameter
    function f_H(fCpbw, fCpinf, fntait, fBtait)
!$acc routine seq
        real(kind(0d0)), intent(IN) :: fCpbw, fCpinf, fntait, fBtait

        real(kind(0d0))             :: tmp1, tmp2, tmp3
        real(kind(0d0))             :: f_H

        tmp1 = (fntait - 1.d0)/fntait
        tmp2 = (fCpbw/(1.d0 + fBtait) + 1.d0)**tmp1
        tmp3 = (fCpinf/(1.d0 + fBtait) + 1.d0)**tmp1

        f_H = (tmp2 - tmp3)*fntait*(1.d0 + fBtait)/(fntait - 1.d0)

    end function f_H

    !> Function that computes the sound speed for the bubble
        !! @param fCpinf Driving bubble pressure
        !! @param fntait Tait EOS parameter
        !! @param fBtait Tait EOS parameter
        !! @param fH Bubble enthalpy
    function f_cgas(fCpinf, fntait, fBtait, fH)
!$acc routine seq
        real(kind(0d0)), intent(IN) :: fCpinf, fntait, fBtait, fH

        real(kind(0d0))             :: tmp
        real(kind(0d0))             :: f_cgas

        ! get sound speed for Gilmore equations "C" -> c_gas
        tmp = (fCpinf/(1.d0 + fBtait) + 1.d0)**((fntait - 1.d0)/fntait)
        tmp = fntait*(1.d0 + fBtait)*tmp

        f_cgas = DSQRT(tmp + (fntait - 1.d0)*fH)

    end function f_cgas

    !>  Function that computes the time derivative of the driving pressure
        !!  @param fRho Local liquid density
        !!  @param fP Local pressure
        !!  @param falf Local void fraction
        !!  @param fntait Tait EOS parameter
        !!  @param fBtait Tait EOS parameter
        !!  @param advsrc Advection equation source term
        !!  @param divu Divergence of velocity
    function f_cpinfdot(fRho, fP, falf, fntait, fBtait, advsrc, divu)
!$acc routine seq
        real(kind(0d0)), intent(IN) :: fRho, fP, falf, fntait, fBtait, advsrc, divu

        real(kind(0d0))             :: c2_liquid
        real(kind(0d0))             :: f_cpinfdot

        ! get sound speed squared for liquid (only needed for pbdot)
        ! c_l^2 = gam (p+B) / (rho*(1-alf))
        if (mpp_lim) then
            c2_liquid = fntait*(fP + fBtait)/fRho
        else
            c2_liquid = fntait*(fP + fBtait)/(fRho*(1.d0 - falf))
        end if

        ! \dot{Cp_inf} = rho sound^2 (alf_src - divu)
        f_cpinfdot = fRho*c2_liquid*(advsrc - divu)

    end function f_cpinfdot

    !>  Function that computes the time derivative of the enthalpy
        !!  @param fCpbw Bubble wall pressure
        !!  @param fCpinf Driving bubble pressure
        !!  @param fCpinf_dot Time derivative of the driving pressure
        !!  @param fntait Tait EOS parameter
        !!  @param fBtait Tait EOS parameter
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fpbdot Time derivative of the internal bubble pressure
    function f_Hdot(fCpbw, fCpinf, fCpinf_dot, fntait, fBtait, fR, fV, fR0, fpbdot)
!$acc routine seq
        real(kind(0d0)), intent(IN) :: fCpbw, fCpinf, fCpinf_dot, fntait, fBtait
        real(kind(0d0)), intent(IN) :: fR, fV, fR0, fpbdot

        real(kind(0d0))             :: tmp1, tmp2
        real(kind(0d0))             :: f_Hdot

        if (polytropic) then
            tmp1 = (fR0/fR)**(3.d0*gam)
            tmp1 = -3.d0*gam*(Ca + 2d0/Web/fR0)*tmp1*fV/fR
        else
            tmp1 = fpbdot
        end if
        tmp2 = (2.d0/Web + 4.d0*Re_inv*fV)*fV/(fR**2.d0)

        f_Hdot = &
            (fCpbw/(1.d0 + fBtait) + 1.d0)**(-1.d0/fntait)*(tmp1 + tmp2) &
            - (fCpinf/(1.d0 + fBtait) + 1.d0)**(-1.d0/fntait)*fCpinf_dot

        ! Hdot = (Cpbw/(1+B) + 1)^(-1/n_tait)*(-3 gam)*(R0/R)^(3gam) V/R
        !f_Hdot = ((fCpbw/(1d0+fBtait)+1.d0)**(-1.d0/fntait))*(-3.d0)*gam * &
        !            ( (fR0/fR)**(3.d0*gam ))*(fV/fR)

        ! Hdot = Hdot - (Cpinf/(1+B) + 1)^(-1/n_tait) Cpinfdot
        !f_Hdot = f_Hdot - ((fCpinf/(1.d0+fBtait)+1.d0)**(-1.d0/fntait))*fCpinf_dot

    end function f_Hdot

    !>  Function that computes the bubble radial acceleration for Rayleigh-Plesset bubbles
        !!  @param fCp Driving pressure
        !!  @param fRho Current density
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fCpbw Boundary wall pressure
    function f_rddot_RP(fCp, fRho, fR, fV, fR0, fCpbw)
!$acc routine seq
        real(kind(0d0)), intent(IN) :: fCp, fRho, fR, fV, fR0, fCpbw
        real(kind(0d0))             :: f_rddot_RP

            !! rddot = (1/r) (  -3/2 rdot^2 + ((r0/r)^3\gamma - Cp)/rho )
            !! rddot = (1/r) (  -3/2 rdot^2 + (tmp1 - Cp)/rho )
            !! rddot = (1/r) (  tmp2 )

        f_rddot_RP = (-1.5d0*(fV**2d0) + (fCpbw - fCp)/fRho)/fR

        if (Re_inv /= dflt_real) f_rddot_RP = f_rddot_RP - 4d0*Re_inv*fv/(fr**2d0)/fRho

    end function f_rddot_RP

    !>  Function that computes the bubble radial acceleration
        !!  @param fCpbw Bubble wall pressure
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fH Current enthalpy
        !!  @param fHdot Current time derivative of the enthalpy
        !!  @param fcgas Current gas sound speed
        !!  @param fntait Tait EOS parameter
        !!  @param fBtait Tait EOS parameter
    function f_rddot(fCpbw, fR, fV, fH, fHdot, fcgas, fntait, fBtait)
!$acc routine seq
        real(kind(0d0)), intent(IN) :: fCpbw, fR, fV, fH, fHdot
        real(kind(0d0)), intent(IN) :: fcgas, fntait, fBtait

        real(kind(0d0))             :: tmp1, tmp2, tmp3
        real(kind(0d0))             :: f_rddot

        tmp1 = fV/fcgas
        tmp2 = 1.d0 + 4.d0*Re_inv/fcgas/fR*(fCpbw/(1.d0 + fBtait) + 1.d0) &
               **(-1.d0/fntait)
        tmp3 = 1.5d0*fV**2d0*(tmp1/3.d0 - 1.d0) + fH*(1.d0 + tmp1) &
               + fR*fHdot*(1.d0 - tmp1)/fcgas

        f_rddot = tmp3/(fR*(1.d0 - tmp1)*tmp2)

    end function f_rddot

    !>  Function that computes the bubble wall pressure for Keller--Miksis bubbles
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fpb Internal bubble pressure
    function f_cpbw_KM(fR0, fR, fV, fpb)
!$acc routine seq
        real(kind(0d0)), intent(IN) :: fR0, fR, fV, fpb
        real(kind(0d0))             :: f_cpbw_KM

        if (polytropic) then
            f_cpbw_KM = Ca*((fR0/fR)**(3.d0*gam)) - Ca + 1d0
            if (Web /= dflt_real) f_cpbw_KM = f_cpbw_KM + &
                                              (2.d0/(Web*fR0))*((fR0/fR)**(3.d0*gam))
        else
            f_cpbw_KM = fpb
        end if

        if (Web /= dflt_real) f_cpbw_KM = f_cpbw_KM - 2.d0/(fR*Web)
        if (Re_inv /= dflt_real) f_cpbw_KM = f_cpbw_KM - 4.d0*Re_inv*fV/fR

    end function f_cpbw_KM

    !>  Function that computes the bubble radial acceleration for Keller--Miksis bubbles
        !!  @param fpbdot Time-derivative of internal bubble pressure
        !!  @param fCp Driving pressure
        !!  @param fCpbw Bubble wall pressure
        !!  @param fRho Current density
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fC Current sound speed
    function f_rddot_KM(fpbdot, fCp, fCpbw, fRho, fR, fV, fR0, fC)
!$acc routine seq
        real(kind(0d0)), intent(IN) :: fpbdot, fCp, fCpbw
        real(kind(0d0)), intent(IN) :: fRho, fR, fV, fR0, fC

        real(kind(0d0))             :: tmp1, tmp2, cdot_star
        real(kind(0d0))             :: f_rddot_KM

        if (polytropic) then
            cdot_star = -3d0*gam*Ca*((fR0/fR)**(3d0*gam))*fV/fR
            if (Web /= dflt_real) cdot_star = cdot_star - &
                                              3d0*gam*(2d0/(Web*fR0))*((fR0/fR)**(3d0*gam))*fV/fR
        else
            cdot_star = fpbdot
        end if

        if (Web /= dflt_real) cdot_star = cdot_star + (2d0/Web)*fV/(fR**2d0)
        if (Re_inv /= dflt_real) cdot_star = cdot_star + 4d0*Re_inv*((fV/fR)**2d0)

        tmp1 = fV/fC
        tmp2 = 1.5d0*(fV**2d0)*(tmp1/3d0 - 1d0) + &
               (1d0 + tmp1)*(fCpbw - fCp)/fRho + &
               cdot_star*fR/(fRho*fC)

        if (Re_inv == dflt_real) then
            f_rddot_KM = tmp2/(fR*(1d0 - tmp1))
        else
            f_rddot_KM = tmp2/(fR*(1d0 - tmp1) + 4d0*Re_inv/(fRho*fC))
        end if

    end function f_rddot_KM

    !>  Subroutine that computes bubble wall properties for vapor bubbles
    !>  @param pb Internal bubble pressure
    !>  @param iR0 Current bubble size index
    subroutine s_bwproperty(pb, iR0)
!$acc routine seq
        real(kind(0.d0)), intent(IN) :: pb
        integer, intent(IN) :: iR0

        real(kind(0.d0)) :: x_vw

        ! mass fraction of vapor
        chi_vw = 1.d0/(1.d0 + R_v/R_n*(pb/pv - 1.d0))
        ! mole fraction of vapor & thermal conductivity of gas mixture
        x_vw = M_n*chi_vw/(M_v + (M_n - M_v)*chi_vw)
        k_mw = x_vw*k_v(iR0)/(x_vw + (1.d0 - x_vw)*phi_vn) &
               + (1.d0 - x_vw)*k_n(iR0)/(x_vw*phi_nv + 1.d0 - x_vw)
        ! gas mixture density
        rho_mw = pv/(chi_vw*R_v*Tw)

    end subroutine s_bwproperty

    !>  Function that computes the vapour flux
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fmass_v Current mass of vapour
        !!  @param iR0 Bubble size index
    function f_vflux(fR, fV, fmass_v, iR0)
!$acc routine seq
        real(kind(0.d0)), intent(IN) :: fR
        real(kind(0.d0)), intent(IN) :: fV
        real(kind(0.d0)), intent(IN) :: fmass_v
        integer, intent(IN) :: iR0

        real(kind(0.d0)) :: chi_bar
        real(kind(0.d0)) :: grad_chi
        real(kind(0.d0)) :: f_vflux

        if (thermal == 3) then !transfer
            ! constant transfer model
            chi_bar = fmass_v/(fmass_v + mass_n0(iR0))
            grad_chi = -Re_trans_c(iR0)*(chi_bar - chi_vw)
            f_vflux = rho_mw*grad_chi/Pe_c/(1.d0 - chi_vw)/fR
        else
            ! polytropic
            f_vflux = pv*fV/(R_v*Tw)
        end if

    end function f_vflux

    !>  Function that computes the time derivative of
        !!  the internal bubble pressure
        !!  @param fvflux Vapour flux
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fpb Current internal bubble pressure
        !!  @param fmass_v Current mass of vapour
        !!  @param iR0 Bubble size index
    function f_bpres_dot(fvflux, fR, fV, fpb, fmass_v, iR0)
!$acc routine seq
        real(kind(0.d0)), intent(IN) :: fvflux
        real(kind(0.d0)), intent(IN) :: fR
        real(kind(0.d0)), intent(IN) :: fV
        real(kind(0.d0)), intent(IN) :: fpb
        real(kind(0.d0)), intent(IN) :: fmass_v
        integer, intent(IN) :: iR0

        real(kind(0.d0)) :: T_bar
        real(kind(0.d0)) :: grad_T
        real(kind(0.d0)) :: tmp1, tmp2
        real(kind(0.d0)) :: f_bpres_dot

        if (thermal == 3) then
            T_bar = Tw*(fpb/pb0(iR0))*(fR/R0(iR0))**3 &
                    *(mass_n0(iR0) + mass_v0(iR0))/(mass_n0(iR0) + fmass_v)
            grad_T = -Re_trans_T(iR0)*(T_bar - Tw)
            f_bpres_dot = 3.d0*gamma_m*(-fV*fpb + fvflux*R_v*Tw &
                                        + pb0(iR0)*k_mw*grad_T/Pe_T(iR0)/fR)/fR
        else
            f_bpres_dot = -3.d0*gamma_m*fV/fR*(fpb - pv)
        end if

    end function f_bpres_dot

end module m_bubbles
