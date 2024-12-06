!>
!! @file m_bubbles.f90
!! @brief Contains module m_bubbles

#:include 'macros.fpp'

!> @brief This module contains the procedures shared by the ensemble-averaged bubble model and the lagrange bubble model.
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

    !>  Function that computes that bubble wall pressure for Gilmore bubbles
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fpb Internal bubble pressure
    function f_cpbw(fR0, fR, fV, fpb)
        !$acc routine seq
        real(kind(0d0)), intent(in) :: fR0, fR, fV, fpb

        real(kind(0d0)) :: f_cpbw

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
        real(kind(0d0)), intent(in) :: fCpbw, fCpinf, fntait, fBtait

        real(kind(0d0)) :: tmp1, tmp2, tmp3
        real(kind(0d0)) :: f_H

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
        real(kind(0d0)), intent(in) :: fCpinf, fntait, fBtait, fH

        real(kind(0d0)) :: tmp
        real(kind(0d0)) :: f_cgas

        ! get sound speed for Gilmore equations "C" -> c_gas
        tmp = (fCpinf/(1.d0 + fBtait) + 1.d0)**((fntait - 1.d0)/fntait)
        tmp = fntait*(1.d0 + fBtait)*tmp

        f_cgas = dsqrt(tmp + (fntait - 1.d0)*fH)

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
        real(kind(0d0)), intent(in) :: fRho, fP, falf, fntait, fBtait, advsrc, divu

        real(kind(0d0)) :: c2_liquid
        real(kind(0d0)) :: f_cpinfdot

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
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fpbdot Time derivative of the internal bubble pressure
    function f_Hdot(fCpbw, fCpinf, fCpinf_dot, fntait, fBtait, fR, fV, fR0, fpbdot)
        !$acc routine seq
        real(kind(0d0)), intent(in) :: fCpbw, fCpinf, fCpinf_dot, fntait, fBtait
        real(kind(0d0)), intent(in) :: fR, fV, fR0, fpbdot

        real(kind(0d0)) :: tmp1, tmp2
        real(kind(0d0)) :: f_Hdot

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
        real(kind(0d0)), intent(in) :: fCp, fRho, fR, fV, fR0, fCpbw

        real(kind(0d0)) :: f_rddot_RP

            !! rddot = (1/r) (  -3/2 rdot^2 + ((r0/r)^3\gamma - Cp)/rho )
            !! rddot = (1/r) (  -3/2 rdot^2 + (tmp1 - Cp)/rho )
            !! rddot = (1/r) (  tmp2 )

        f_rddot_RP = (-1.5d0*(fV**2d0) + (fCpbw - fCp)/fRho)/fR

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
    function f_rddot_G(fCpbw, fR, fV, fH, fHdot, fcgas, fntait, fBtait)
        !$acc routine seq
        real(kind(0d0)), intent(in) :: fCpbw, fR, fV, fH, fHdot
        real(kind(0d0)), intent(in) :: fcgas, fntait, fBtait

        real(kind(0d0)) :: tmp1, tmp2, tmp3
        real(kind(0d0)) :: f_rddot_G

        tmp1 = fV/fcgas
        tmp2 = 1.d0 + 4.d0*Re_inv/fcgas/fR*(fCpbw/(1.d0 + fBtait) + 1.d0) &
               **(-1.d0/fntait)
        tmp3 = 1.5d0*fV**2d0*(tmp1/3.d0 - 1.d0) + fH*(1.d0 + tmp1) &
               + fR*fHdot*(1.d0 - tmp1)/fcgas

        f_rddot_G = tmp3/(fR*(1.d0 - tmp1)*tmp2)

    end function f_rddot_G

    !>  Function that computes the bubble wall pressure for Keller--Miksis bubbles
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fpb Internal bubble pressure
    function f_cpbw_KM(fR0, fR, fV, fpb)
        !$acc routine seq
        real(kind(0d0)), intent(in) :: fR0, fR, fV, fpb

        real(kind(0d0)) :: f_cpbw_KM

        if (polytropic) then
            f_cpbw_KM = Ca*((fR0/fR)**(3.d0*gam)) - Ca + 1d0
            if (.not. f_is_default(Web)) f_cpbw_KM = f_cpbw_KM + &
                                                     (2.d0/(Web*fR0))*((fR0/fR)**(3.d0*gam))
        else
            f_cpbw_KM = fpb
        end if

        if (.not. f_is_default(Web)) f_cpbw_KM = f_cpbw_KM - 2.d0/(fR*Web)
        if (.not. f_is_default(Re_inv)) f_cpbw_KM = f_cpbw_KM - 4.d0*Re_inv*fV/fR

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
        real(kind(0d0)), intent(in) :: fpbdot, fCp, fCpbw
        real(kind(0d0)), intent(in) :: fRho, fR, fV, fR0, fC

        real(kind(0d0)) :: tmp1, tmp2, cdot_star
        real(kind(0d0)) :: f_rddot_KM

        if (polytropic) then
            cdot_star = -3d0*gam*Ca*((fR0/fR)**(3d0*gam))*fV/fR
            if (.not. f_is_default(Web)) cdot_star = cdot_star - &
                                                     3d0*gam*(2d0/(Web*fR0))*((fR0/fR)**(3d0*gam))*fV/fR
        else
            cdot_star = fpbdot
        end if

        if (.not. bubbles_lagrange) then
            if (.not. f_is_default(Web)) cdot_star = cdot_star + (2d0/Web)*fV/(fR**2d0)
            if (.not. f_is_default(Re_inv)) cdot_star = cdot_star + 4d0*Re_inv*((fV/fR)**2d0)
        end if

        tmp1 = fV/fC
        tmp2 = 1.5d0*(fV**2d0)*(tmp1/3d0 - 1d0) + &
               (1d0 + tmp1)*(fCpbw - fCp)/fRho + &
               cdot_star*fR/(fRho*fC)

        if (bubbles_lagrange) then
            f_rddot_KM = tmp2/(fR*(1d0 - tmp1))
            return
        end if

        if (f_is_default(Re_inv)) then
            f_rddot_KM = tmp2/(fR*(1d0 - tmp1))
        else
            f_rddot_KM = tmp2/(fR*(1d0 - tmp1) + 4d0*Re_inv/(fRho*fC))
        end if

    end function f_rddot_KM

    !>  Subroutine that computes bubble wall properties for vapor bubbles
        !!  @param pb Internal bubble pressure
        !!  @param iR0 Current bubble size index
    subroutine s_bwproperty(pb, iR0)
        !$acc routine seq
        real(kind(0.d0)), intent(in) :: pb
        integer, intent(in) :: iR0

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
        !!  @param iR0 Bubble size index (EE) or bubble identifier (EL)
        !!  @param fmass_n Current gas mass (EL)
        !!  @param fbeta_c Mass transfer coefficient (EL)
        !!  @param fconc_v Current vapor concentration (EL)
    function f_vflux(fR, fV, fmass_v, iR0, fmass_n, fbeta_c, fconc_v)
        !$acc routine seq
        real(kind(0.d0)), intent(in) :: fR
        real(kind(0.d0)), intent(in) :: fV
        real(kind(0.d0)), intent(in) :: fmass_v
        integer, intent(in) :: iR0
        real(kind(0.d0)), intent(in), optional :: fmass_n, fbeta_c, fconc_v

        real(kind(0.d0)) :: chi_bar
        real(kind(0.d0)) :: rho_mw_lag
        real(kind(0.d0)) :: grad_chi
        real(kind(0.d0)) :: f_vflux

        if (thermal == 3) then !transfer
            ! constant transfer model
            if (bubbles_lagrange) then
                chi_bar = fmass_v/(fmass_v + fmass_n)
                grad_chi = (chi_bar - fconc_v)
                rho_mw_lag = (fmass_n + fmass_v)/(4.0d0/3.0d0*pi*fR**3d0)
                f_vflux = -fbeta_c*rho_mw_lag*grad_chi/(1.d0 - fconc_v)/fR
                return
            end if
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
        !!  @param iR0 Bubble size index (EE) or bubble identifier (EL)
        !!  @param fbeta_t Mass transfer coefficient (EL)
        !!  @param fR_m Mixture gas constant (EL)
        !!  @param fgamma_m Mixture gamma (EL)
        !!  @param fconc_v Current vapor concentration (EL)
    function f_bpres_dot(fvflux, fR, fV, fpb, fmass_v, iR0, fbeta_t, fR_m, fgamma_m, fconc_v)
        !$acc routine seq
        real(kind(0.d0)), intent(in) :: fvflux
        real(kind(0.d0)), intent(in) :: fR
        real(kind(0.d0)), intent(in) :: fV
        real(kind(0.d0)), intent(in) :: fpb
        real(kind(0.d0)), intent(in) :: fmass_v
        integer, intent(in) :: iR0
        real(kind(0.d0)), intent(in), optional :: fbeta_t, fR_m, fgamma_m, fconc_v

        real(kind(0.d0)) :: T_bar
        real(kind(0.d0)) :: grad_T
        real(kind(0.d0)) :: f_bpres_dot
        real(kind(0.d0)) :: heatflux

        if (thermal == 3) then
            if (bubbles_lagrange) then
                T_bar = fpb*(4.0d0/3.0d0*pi*fR**3d0)/fR_m
                grad_T = -fbeta_t*(T_bar - Tw)
                heatflux = (fgamma_m - 1.0d0)/fgamma_m*grad_T/fR
                f_bpres_dot = 3.d0*fgamma_m*(-fV*fpb + fvflux*R_v*Tw &
                                             + heatflux)/fR
                return
            end if
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
