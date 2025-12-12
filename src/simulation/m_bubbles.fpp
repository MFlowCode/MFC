!>
!! @file m_bubbles.f90
!! @brief Contains module m_bubbles

#:include 'macros.fpp'

!> @brief This module contains procedures shared by the ensemble-averaged and volume-averaged bubble models.
module m_bubbles

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    use m_helper_basic         !< Functions to compare floating point numbers

    implicit none

    real(wp) :: chi_vw  !< Bubble wall properties (Ando 2010)
    real(wp) :: k_mw    !< Bubble wall properties (Ando 2010)
    real(wp) :: rho_mw  !< Bubble wall properties (Ando 2010)
    $:GPU_DECLARE(create='[chi_vw,k_mw,rho_mw]')

contains

    !> Function that computes the bubble radial acceleration based on bubble models
        !!  @param fRho Current density
        !!  @param fP Current driving pressure
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fpb Internal bubble pressure
        !!  @param fpbdot Time-derivative of internal bubble pressure
        !!  @param alf bubble volume fraction
        !!  @param fntait Tait EOS parameter
        !!  @param fBtait Tait EOS parameter
        !!  @param f_bub_adv_src Source for bubble volume fraction
        !!  @param f_divu Divergence of velocity
        !!  @param fCson Speed of sound from fP (EL)
    elemental function f_rddot(fRho, fP, fR, fV, fR0, fpb, fpbdot, alf, fntait, fBtait, f_bub_adv_src, f_divu, fCson)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: fRho, fP, fR, fV, fR0, fpb, fpbdot, alf
        real(wp), intent(in) :: fntait, fBtait, f_bub_adv_src, f_divu
        real(wp), intent(in) :: fCson

        real(wp) :: fCpbw, fCpinf, fCpinf_dot, fH, fHdot, c_gas, c_liquid
        real(wp) :: f_rddot

        if (bubble_model == 1) then
            ! Gilmore bubbles
            fCpinf = fP - Eu
            fCpbw = f_cpbw(fR0, fR, fV, fpb)
            fH = f_H(fCpbw, fCpinf, fntait, fBtait)
            c_gas = f_cgas(fCpinf, fntait, fBtait, fH)
            fCpinf_dot = f_cpinfdot(fRho, fP, alf, fntait, fBtait, f_bub_adv_src, f_divu)
            fHdot = f_Hdot(fCpbw, fCpinf, fCpinf_dot, fntait, fBtait, fR, fV, fR0, fpbdot)
            f_rddot = f_rddot_G(fCpbw, fR, fV, fH, fHdot, c_gas, fntait, fBtait)
        else if (bubble_model == 2) then
            ! Keller-Miksis bubbles
            fCpinf = fP
            fCpbw = f_cpbw_KM(fR0, fR, fV, fpb)
            if (bubbles_euler) then
                c_liquid = sqrt(fntait*(fP + fBtait)/(fRho*(1._wp - alf)))
            else
                c_liquid = fCson
            end if
            f_rddot = f_rddot_KM(fpbdot, fCpinf, fCpbw, fRho, fR, fV, fR0, c_liquid)
        else if (bubble_model == 3) then
            ! Rayleigh-Plesset bubbles
            fCpbw = f_cpbw_KM(fR0, fR, fV, fpb)
            f_rddot = f_rddot_RP(fP, fRho, fR, fV, fCpbw)
        else
            ! Default: No bubble dynamics
            f_rddot = 0._wp
        end if

    end function f_rddot

    !>  Function that computes that bubble wall pressure for Gilmore bubbles
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fpb Internal bubble pressure
    elemental function f_cpbw(fR0, fR, fV, fpb)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: fR0, fR, fV, fpb

        real(wp) :: f_cpbw

        if (polytropic) then
            f_cpbw = (Ca + 2._wp/Web/fR0)*((fR0/fR)**(3._wp*gam)) - Ca - 4._wp*Re_inv*fV/fR - 2._wp/(fR*Web)
        else
            f_cpbw = fpb - 1._wp - 4._wp*Re_inv*fV/fR - 2._wp/(fR*Web)
        end if

    end function f_cpbw

    !>  Function that computes the bubble enthalpy
        !!  @param fCpbw Bubble wall pressure
        !!  @param fCpinf Driving bubble pressure
        !!  @param fntait Tait EOS parameter
        !!  @param fBtait Tait EOS parameter
    elemental function f_H(fCpbw, fCpinf, fntait, fBtait)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: fCpbw, fCpinf, fntait, fBtait

        real(wp) :: tmp1, tmp2, tmp3
        real(wp) :: f_H

        tmp1 = (fntait - 1._wp)/fntait
        tmp2 = (fCpbw/(1._wp + fBtait) + 1._wp)**tmp1
        tmp3 = (fCpinf/(1._wp + fBtait) + 1._wp)**tmp1

        f_H = (tmp2 - tmp3)*fntait*(1._wp + fBtait)/(fntait - 1._wp)

    end function f_H

    !> Function that computes the sound speed for the bubble
        !! @param fCpinf Driving bubble pressure
        !! @param fntait Tait EOS parameter
        !! @param fBtait Tait EOS parameter
        !! @param fH Bubble enthalpy
    elemental function f_cgas(fCpinf, fntait, fBtait, fH)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: fCpinf, fntait, fBtait, fH

        real(wp) :: tmp
        real(wp) :: f_cgas

        ! get sound speed for Gilmore equations "C" -> c_gas
        tmp = (fCpinf/(1._wp + fBtait) + 1._wp)**((fntait - 1._wp)/fntait)
        tmp = fntait*(1._wp + fBtait)*tmp

        f_cgas = sqrt(tmp + (fntait - 1._wp)*fH)

    end function f_cgas

    !>  Function that computes the time derivative of the driving pressure
        !!  @param fRho Local liquid density
        !!  @param fP Local pressure
        !!  @param falf Local void fraction
        !!  @param fntait Tait EOS parameter
        !!  @param fBtait Tait EOS parameter
        !!  @param advsrc Advection equation source term
        !!  @param divu Divergence of velocity
    elemental function f_cpinfdot(fRho, fP, falf, fntait, fBtait, advsrc, divu)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: fRho, fP, falf, fntait, fBtait, advsrc, divu

        real(wp) :: c2_liquid
        real(wp) :: f_cpinfdot

        ! get sound speed squared for liquid (only needed for pbdot)
        ! c_l^2 = gam (p+B) / (rho*(1-alf))
        if (mpp_lim) then
            c2_liquid = fntait*(fP + fBtait)/fRho
        else
            c2_liquid = fntait*(fP + fBtait)/(fRho*(1._wp - falf))
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
    elemental function f_Hdot(fCpbw, fCpinf, fCpinf_dot, fntait, fBtait, fR, fV, fR0, fpbdot)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: fCpbw, fCpinf, fCpinf_dot, fntait, fBtait
        real(wp), intent(in) :: fR, fV, fR0, fpbdot

        real(wp) :: tmp1, tmp2
        real(wp) :: f_Hdot

        if (polytropic) then
            tmp1 = (fR0/fR)**(3._wp*gam)
            tmp1 = -3._wp*gam*(Ca + 2._wp/Web/fR0)*tmp1*fV/fR
        else
            tmp1 = fpbdot
        end if
        tmp2 = (2._wp/Web + 4._wp*Re_inv*fV)*fV/(fR**2._wp)

        f_Hdot = &
            (fCpbw/(1._wp + fBtait) + 1._wp)**(-1._wp/fntait)*(tmp1 + tmp2) &
            - (fCpinf/(1._wp + fBtait) + 1._wp)**(-1._wp/fntait)*fCpinf_dot

        ! Hdot = (Cpbw/(1+B) + 1)^(-1/n_tait)*(-3 gam)*(R0/R)^(3gam) V/R
        !f_Hdot = ((fCpbw/(1._wp+fBtait)+1._wp)**(-1._wp/fntait))*(-3._wp)*gam * &
        !            ( (fR0/fR)**(3._wp*gam ))*(fV/fR)

        ! Hdot = Hdot - (Cpinf/(1+B) + 1)^(-1/n_tait) Cpinfdot
        !f_Hdot = f_Hdot - ((fCpinf/(1._wp+fBtait)+1._wp)**(-1._wp/fntait))*fCpinf_dot

    end function f_Hdot

    !>  Function that computes the bubble radial acceleration for Rayleigh-Plesset bubbles
        !!  @param fCp Driving pressure
        !!  @param fRho Current density
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fCpbw Boundary wall pressure
    elemental function f_rddot_RP(fCp, fRho, fR, fV, fCpbw)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: fCp, fRho, fR, fV, fCpbw

        real(wp) :: f_rddot_RP

            !! rddot = (1/r) (  -3/2 rdot^2 + ((r0/r)^3\gamma - Cp)/rho )
            !! rddot = (1/r) (  -3/2 rdot^2 + (tmp1 - Cp)/rho )
            !! rddot = (1/r) (  tmp2 )

        f_rddot_RP = (-1.5_wp*(fV**2._wp) + (fCpbw - fCp)/fRho)/fR

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
    elemental function f_rddot_G(fCpbw, fR, fV, fH, fHdot, fcgas, fntait, fBtait)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: fCpbw, fR, fV, fH, fHdot
        real(wp), intent(in) :: fcgas, fntait, fBtait

        real(wp) :: tmp1, tmp2, tmp3
        real(wp) :: f_rddot_G

        tmp1 = fV/fcgas
        tmp2 = 1._wp + 4._wp*Re_inv/fcgas/fR*(fCpbw/(1._wp + fBtait) + 1._wp) &
               **(-1._wp/fntait)
        tmp3 = 1.5_wp*fV**2._wp*(tmp1/3._wp - 1._wp) + fH*(1._wp + tmp1) &
               + fR*fHdot*(1._wp - tmp1)/fcgas

        f_rddot_G = tmp3/(fR*(1._wp - tmp1)*tmp2)

    end function f_rddot_G

    !>  Function that computes the bubble wall pressure for Keller--Miksis bubbles
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fpb Internal bubble pressure
    elemental function f_cpbw_KM(fR0, fR, fV, fpb)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: fR0, fR, fV, fpb
        real(wp) :: f_cpbw_KM

        if (polytropic) then
            f_cpbw_KM = Ca*((fR0/fR)**(3._wp*gam)) - Ca + Eu
            if (.not. f_is_default(Web)) f_cpbw_KM = f_cpbw_KM + &
                                                     (2._wp/(Web*fR0))*((fR0/fR)**(3._wp*gam))
        else
            f_cpbw_KM = fpb
        end if

        if (.not. f_is_default(Web)) f_cpbw_KM = f_cpbw_KM - 2._wp/(fR*Web)
        if (.not. f_is_default(Re_inv)) f_cpbw_KM = f_cpbw_KM - 4._wp*Re_inv*fV/fR

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
    elemental function f_rddot_KM(fpbdot, fCp, fCpbw, fRho, fR, fV, fR0, fC)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: fpbdot, fCp, fCpbw
        real(wp), intent(in) :: fRho, fR, fV, fR0, fC

        real(wp) :: tmp1, tmp2, cdot_star
        real(wp) :: f_rddot_KM
        if (polytropic) then
            cdot_star = -3._wp*gam*Ca*((fR0/fR)**(3._wp*gam))*fV/fR
            if (.not. f_is_default(Web)) cdot_star = cdot_star - &
                                                     3._wp*gam*(2._wp/(Web*fR0))*((fR0/fR)**(3._wp*gam))*fV/fR
        else
            cdot_star = fpbdot
        end if

        if (.not. f_is_default(Web)) cdot_star = cdot_star + (2._wp/Web)*fV/(fR**2._wp)
        if (.not. f_is_default(Re_inv)) cdot_star = cdot_star + 4._wp*Re_inv*((fV/fR)**2._wp)

        tmp1 = fV/fC
        tmp2 = 1.5_wp*(fV**2._wp)*(tmp1/3._wp - 1._wp) + &
               (1._wp + tmp1)*(fCpbw - fCp)/fRho + &
               cdot_star*fR/(fRho*fC)

        if (f_is_default(Re_inv)) then
            f_rddot_KM = tmp2/(fR*(1._wp - tmp1))
        else
            f_rddot_KM = tmp2/(fR*(1._wp - tmp1) + 4._wp*Re_inv/(fRho*fC))
        end if

    end function f_rddot_KM

    !>  Subroutine that computes bubble wall properties for vapor bubbles
        !!  @param pb Internal bubble pressure
        !!  @param iR0 Current bubble size index
    elemental subroutine s_bwproperty(pb_in, iR0, chi_vw_out, k_mw_out, rho_mw_out)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: pb_in
        integer, intent(in) :: iR0
        real(wp), intent(out) :: chi_vw_out  !< Bubble wall properties (Ando 2010)
        real(wp), intent(out) :: k_mw_out    !< Bubble wall properties (Ando 2010)
        real(wp), intent(out) :: rho_mw_out  !< Bubble wall properties (Ando 2010)
        real(wp) :: x_vw

        ! mass fraction of vapor
        chi_vw_out = 1._wp/(1._wp + R_v/R_g*(pb_in/pv - 1._wp))
        ! mole fraction of vapor & thermal conductivity of gas mixture
        x_vw = M_g*chi_vw_out/(M_v + (M_g - M_v)*chi_vw_out)
        k_mw_out = x_vw*k_v(iR0)/(x_vw + (1._wp - x_vw)*phi_vg) &
                   + (1._wp - x_vw)*k_g(iR0)/(x_vw*phi_gv + 1._wp - x_vw)
        ! gas mixture density
        rho_mw_out = pv/(chi_vw_out*R_v*Tw)

    end subroutine s_bwproperty

    !>  Function that computes the vapour flux
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fpb
        !!  @param fmass_v Current mass of vapour
        !!  @param iR0 Bubble size index (EE) or bubble identifier (EL)
        !!  @param fmass_g Current gas mass (EL)
        !!  @param fbeta_c Mass transfer coefficient (EL)
        !!  @param fR_m Mixture gas constant (EL)
        !!  @param fgamma_m Mixture gamma (EL)
    elemental subroutine s_vflux(fR, fV, fpb, fmass_v, iR0, vflux, fmass_g, fbeta_c, fR_m, fgamma_m)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: fR
        real(wp), intent(in) :: fV
        real(wp), intent(in) :: fpb
        real(wp), intent(in) :: fmass_v
        integer, intent(in) :: iR0
        real(wp), intent(out) :: vflux
        real(wp), intent(in), optional :: fmass_g, fbeta_c
        real(wp), intent(out), optional :: fR_m, fgamma_m

        real(wp) :: chi_bar
        real(wp) :: rho_mw_lag
        real(wp) :: grad_chi
        real(wp) :: conc_v

        if (thermal == 3) then !transfer
            ! constant transfer model
            if (bubbles_lagrange) then
                ! Mixture properties (gas+vapor) in the bubble
                conc_v = fmass_v/(fmass_v + fmass_g)
                if (lag_params%massTransfer_model) then
                    conc_v = 1._wp/(1._wp + (R_v/R_g)*(fpb/pv - 1._wp))
                end if
                fR_m = (fmass_g*R_g + fmass_v*R_v)
                fgamma_m = conc_v*gam_v + (1._wp - conc_v)*gam_g

                ! Vapor flux
                chi_bar = fmass_v/(fmass_v + fmass_g)
                grad_chi = (chi_bar - conc_v)
                rho_mw_lag = (fmass_g + fmass_v)/(4._wp/3._wp*pi*fR**3._wp)
                vflux = 0._wp
                if (lag_params%massTransfer_model) then
                    vflux = -fbeta_c*rho_mw_lag*grad_chi/(1._wp - conc_v)/fR
                end if
            else
                chi_bar = fmass_v/(fmass_v + mass_g0(iR0))
                grad_chi = -Re_trans_c(iR0)*(chi_bar - chi_vw)
                vflux = rho_mw*grad_chi/Pe_c/(1._wp - chi_vw)/fR
            end if
        else
            ! polytropic
            vflux = pv*fV/(R_v*Tw)
        end if

    end subroutine s_vflux

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
    elemental function f_bpres_dot(fvflux, fR, fV, fpb, fmass_v, iR0, fbeta_t, fR_m, fgamma_m)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: fvflux
        real(wp), intent(in) :: fR
        real(wp), intent(in) :: fV
        real(wp), intent(in) :: fpb
        real(wp), intent(in) :: fmass_v
        integer, intent(in) :: iR0
        real(wp), intent(in), optional :: fbeta_t, fR_m, fgamma_m

        real(wp) :: T_bar
        real(wp) :: grad_T
        real(wp) :: f_bpres_dot
        real(wp) :: heatflux

        if (thermal == 3) then
            if (bubbles_lagrange) then
                T_bar = fpb*(4._wp/3._wp*pi*fR**3._wp)/fR_m
                grad_T = -fbeta_t*(T_bar - Tw)
                heatflux = (fgamma_m - 1._wp)/fgamma_m*grad_T/fR
                f_bpres_dot = 3._wp*fgamma_m*(-fV*fpb + fvflux*R_v*Tw &
                                              + heatflux)/fR
                return
            end if
            grad_T = -Re_trans_T(iR0)*((fpb/pb0(iR0))*(fR/R0(iR0))**3 &
                                       *(mass_g0(iR0) + mass_v0(iR0))/(mass_g0(iR0) + fmass_v) - 1._wp)
            f_bpres_dot = 3._wp*gam_m*(-fV*fpb + fvflux*R_v*Tw &
                                       + pb0(iR0)*k_mw*grad_T/Pe_T(iR0)/fR)/fR
        else
            f_bpres_dot = -3._wp*gam_m*fV/fR*(fpb - pv)
        end if

    end function f_bpres_dot

    !> Adaptive time stepping routine for subgrid bubbles
        !!  (See Heirer, E. Hairer S.P.Nørsett G. Wanner, Solving Ordinary
        !!  Differential Equations I, Chapter II.4)
        !!  @param fRho Current density
        !!  @param fP Current driving pressure
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble radial velocity
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fpb Internal bubble pressure
        !!  @param fpbdot Time-derivative of internal bubble pressure
        !!  @param alf bubble volume fraction
        !!  @param fntait Tait EOS parameter
        !!  @param fBtait Tait EOS parameter
        !!  @param f_bub_adv_src Source for bubble volume fraction
        !!  @param f_divu Divergence of velocity
        !!  @param bub_id Bubble identifier (EL)
        !!  @param fmass_v Current mass of vapour (EL)
        !!  @param fmass_g Current mass of gas (EL)
        !!  @param fbeta_c Mass transfer coefficient (EL)
        !!  @param fbeta_t Heat transfer coefficient (EL)
        !!  @param fCson Speed of sound (EL)
        !!  @param adap_dt_stop Fail-safe exit if max iteration count reached
    subroutine s_advance_step(fRho, fP, fR, fV, fR0, fpb, fpbdot, alf, &
                              fntait, fBtait, f_bub_adv_src, f_divu, &
                              bub_id, fmass_v, fmass_g, fbeta_c, &
                              fbeta_t, fCson, adap_dt_stop)
        $:GPU_ROUTINE(function_name='s_advance_step',parallelism='[seq]', &
            & cray_inline=True)

        real(wp), intent(inout) :: fR, fV, fpb, fmass_v
        real(wp), intent(in) :: fRho, fP, fR0, fpbdot, alf
        real(wp), intent(in) :: fntait, fBtait, f_bub_adv_src, f_divu
        integer, intent(in) :: bub_id
        real(wp), intent(in) :: fmass_g, fbeta_c, fbeta_t, fCson
        integer, intent(inout) :: adap_dt_stop

        real(wp), dimension(5) :: err !< Error estimates for adaptive time stepping
        real(wp) :: t_new !< Updated time step size
        real(wp) :: h0, h !< Time step size
        real(wp), dimension(4) :: myR_tmp1, myV_tmp1, myR_tmp2, myV_tmp2 !< Bubble radius, radial velocity, and radial acceleration for the inner loop
        real(wp), dimension(4) :: myPb_tmp1, myMv_tmp1, myPb_tmp2, myMv_tmp2 !< Gas pressure and vapor mass for the inner loop (EL)
        real(wp) :: fR2, fV2, fpb2, fmass_v2
        integer :: iter_count

        call s_initial_substep_h(fRho, fP, fR, fV, fR0, fpb, fpbdot, alf, &
                                 fntait, fBtait, f_bub_adv_src, f_divu, fCson, h0)
        h = h0
        ! Advancing one step
        t_new = 0._wp
        iter_count = 0
        adap_dt_stop = 0

        do
            if (t_new + h > 0.5_wp*dt) then
                h = 0.5_wp*dt - t_new
            end if

            ! Advancing one sub-step
            do while (iter_count < adap_dt_max_iters)

                iter_count = iter_count + 1

                ! Advance one sub-step
                call s_advance_substep(err(1), &
                                       fRho, fP, fR, fV, fR0, fpb, fpbdot, alf, &
                                       fntait, fBtait, f_bub_adv_src, f_divu, &
                                       bub_id, fmass_v, fmass_g, fbeta_c, &
                                       fbeta_t, fCson, h, &
                                       myR_tmp1, myV_tmp1, myPb_tmp1, myMv_tmp1)
                if (err(1) > adap_dt_tol) then
                    h = 0.25_wp*h
                    cycle
                end if

                ! Advance one sub-step by advancing two half steps
                call s_advance_substep(err(2), &
                                       fRho, fP, fR, fV, fR0, fpb, fpbdot, alf, &
                                       fntait, fBtait, f_bub_adv_src, f_divu, &
                                       bub_id, fmass_v, fmass_g, fbeta_c, &
                                       fbeta_t, fCson, 0.5_wp*h, &
                                       myR_tmp2, myV_tmp2, myPb_tmp2, myMv_tmp2)
                if (err(2) > adap_dt_tol) then
                    h = 0.25_wp*h
                    cycle
                end if

                fR2 = myR_tmp2(4); fV2 = myV_tmp2(4)
                fpb2 = myPb_tmp2(4); fmass_v2 = myMv_tmp2(4)

                call s_advance_substep(err(3), &
                                       fRho, fP, fR2, fV2, fR0, fpb2, fpbdot, alf, &
                                       fntait, fBtait, f_bub_adv_src, f_divu, &
                                       bub_id, fmass_v2, fmass_g, fbeta_c, &
                                       fbeta_t, fCson, 0.5_wp*h, &
                                       myR_tmp2, myV_tmp2, myPb_tmp2, myMv_tmp2)
                if (err(3) > adap_dt_tol) then
                    h = 0.5_wp*h
                    cycle
                end if

                err(4) = abs((myR_tmp1(4) - myR_tmp2(4))/myR_tmp1(4))
                err(5) = abs((myV_tmp1(4) - myV_tmp2(4))/myV_tmp1(4))
                if (abs(myV_tmp1(4)) < verysmall) err(5) = 0._wp

                ! Determine acceptance/rejection and update step size
                !   Rule 1: err1, err2, err3 < tol
                !   Rule 2: myR_tmp1(4) > 0._wp
                !   Rule 3: abs((myR_tmp1(4) - myR_tmp2(4))/fR) < tol
                !   Rule 4: abs((myV_tmp1(4) - myV_tmp2(4))/fV) < tol
                if ((err(1) <= adap_dt_tol) .and. (err(2) <= adap_dt_tol) .and. &
                    (err(3) <= adap_dt_tol) .and. (err(4) <= adap_dt_tol) .and. &
                    (err(5) <= adap_dt_tol) .and. myR_tmp1(4) > 0._wp) then

                    ! Accepted. Finalize the sub-step
                    t_new = t_new + h

                    ! Update R and V
                    fR = myR_tmp1(4)
                    fV = myV_tmp1(4)

                    if (bubbles_lagrange) then
                        ! Update pb and mass_v
                        fpb = myPb_tmp1(4)
                        fmass_v = myMv_tmp1(4)
                    end if

                    ! Update step size for the next sub-step
                    h = h*min(2._wp, max(0.5_wp, (adap_dt_tol/err(1))**(1._wp/3._wp)))

                    exit
                else
                    ! Rejected. Update step size for the next try on sub-step
                    if (err(2) <= adap_dt_tol) then
                        h = 0.5_wp*h
                    else
                        h = 0.25_wp*h
                    end if
                end if
            end do

            ! Exit the loop if the final time reached dt
            if (f_approx_equal(t_new, 0.5_wp*dt) .or. iter_count >= adap_dt_max_iters) exit

        end do

        if (iter_count >= adap_dt_max_iters) adap_dt_stop = 1

    end subroutine s_advance_step

    !> Choose the initial time step size for the adaptive time stepping routine
        !!  (See Heirer, E. Hairer S.P.Nørsett G. Wanner, Solving Ordinary
        !!  Differential Equations I, Chapter II.4)
        !!  @param fRho Current density
        !!  @param fP Current driving pressure
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fpb Internal bubble pressure
        !!  @param fpbdot Time-derivative of internal bubble pressure
        !!  @param alf bubble volume fraction
        !!  @param fntait Tait EOS parameter
        !!  @param fBtait Tait EOS parameter
        !!  @param f_bub_adv_src Source for bubble volume fraction
        !!  @param f_divu Divergence of velocity
        !!  @param fCson Speed of sound (EL)
        !!  @param h Time step size
    subroutine s_initial_substep_h(fRho, fP, fR, fV, fR0, fpb, fpbdot, alf, &
                                   fntait, fBtait, f_bub_adv_src, f_divu, &
                                   fCson, h)
        $:GPU_ROUTINE(function_name='s_initial_substep_h',parallelism='[seq]', &
            & cray_inline=True)

        real(wp), intent(IN) :: fRho, fP, fR, fV, fR0, fpb, fpbdot, alf
        real(wp), intent(IN) :: fntait, fBtait, f_bub_adv_src, f_divu
        real(wp), intent(IN) :: fCson
        real(wp), intent(OUT) :: h

        real(wp), dimension(2) :: h_size !< Time step size (h0, h1)
        real(wp), dimension(3) :: d_norms !< norms (d_0, d_1, d_2)
        real(wp), dimension(2) :: myR_tmp, myV_tmp, myA_tmp !< Bubble radius, radial velocity, and radial acceleration

        ! Determine the starting time step
        ! Evaluate f(x0,y0)
        myR_tmp(1) = fR
        myV_tmp(1) = fV
        myA_tmp(1) = f_rddot(fRho, fP, myR_tmp(1), myV_tmp(1), fR0, &
                             fpb, fpbdot, alf, fntait, fBtait, &
                             f_bub_adv_src, f_divu, &
                             fCson)

        ! Compute d_0 = ||y0|| and d_1 = ||f(x0,y0)||
        d_norms(1) = sqrt((myR_tmp(1)**2._wp + myV_tmp(1)**2._wp)/2._wp)
        d_norms(2) = sqrt((myV_tmp(1)**2._wp + myA_tmp(1)**2._wp)/2._wp)
        if (d_norms(1) < threshold_first_guess .or. d_norms(2) < threshold_first_guess) then
            h_size(1) = small_guess
        else
            h_size(1) = scale_guess*(d_norms(1)/d_norms(2))
        end if

        ! Evaluate f(x0+h0,y0+h0*f(x0,y0))
        myR_tmp(2) = myR_tmp(1) + h_size(1)*myV_tmp(1)
        myV_tmp(2) = myV_tmp(1) + h_size(1)*myA_tmp(1)
        myA_tmp(2) = f_rddot(fRho, fP, myR_tmp(2), myV_tmp(2), fR0, &
                             fpb, fpbdot, alf, fntait, fBtait, &
                             f_bub_adv_src, f_divu, &
                             fCson)

        ! Compute d_2 = ||f(x0+h0,y0+h0*f(x0,y0))-f(x0,y0)||/h0
        d_norms(3) = sqrt(((myV_tmp(2) - myV_tmp(1))**2._wp + (myA_tmp(2) - myA_tmp(1))**2._wp)/2._wp)/h_size(1)

        ! Set h1 = (0.01/max(d_1,d_2))^{1/(p+1)}
        !      if max(d_1,d_2) < 1.e-15_wp, h_size(2) = max(1.e-6_wp, h0*1.e-3_wp)
        if (max(d_norms(2), d_norms(3)) < threshold_second_guess) then
            h_size(2) = max(small_guess, h_size(1)*scale_first_guess)
        else
            h_size(2) = (scale_guess/max(d_norms(2), d_norms(3)))**(1._wp/3._wp)
        end if

        h = min(h_size(1)/scale_guess, h_size(2))

    end subroutine s_initial_substep_h

    !>  Integrate bubble variables over the given time step size, h, using a
        !!      third-order accurate embedded Runge–Kutta scheme.
        !!  @param err Estimated error
        !!  @param fRho Current density
        !!  @param fP Current driving pressure
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fpb Internal bubble pressure
        !!  @param fpbdot Time-derivative of internal bubble pressure
        !!  @param alf bubble volume fraction
        !!  @param fntait Tait EOS parameter
        !!  @param fBtait Tait EOS parameter
        !!  @param f_bub_adv_src Source for bubble volume fraction
        !!  @param f_divu Divergence of velocity
        !!  @param bub_id Bubble identifier (EL)
        !!  @param fmass_v Current mass of vapour (EL)
        !!  @param fmass_g Current mass of gas (EL)
        !!  @param fbeta_c Mass transfer coefficient (EL)
        !!  @param fbeta_t Heat transfer coefficient (EL)
        !!  @param fCson Speed of sound (EL)
        !!  @param h Time step size
        !!  @param myR_tmp Bubble radius at each stage
        !!  @param myV_tmp Bubble radial velocity at each stage
        !!  @param myPb_tmp Internal bubble pressure at each stage (EL)
        !!  @param myMv_tmp Mass of vapor in the bubble at each stage (EL)
    subroutine s_advance_substep(err, fRho, fP, fR, fV, fR0, fpb, fpbdot, alf, &
                                 fntait, fBtait, f_bub_adv_src, f_divu, &
                                 bub_id, fmass_v, fmass_g, fbeta_c, &
                                 fbeta_t, fCson, h, &
                                 myR_tmp, myV_tmp, myPb_tmp, myMv_tmp)
        $:GPU_ROUTINE(function_name='s_advance_substep',parallelism='[seq]', &
            & cray_inline=True)

        real(wp), intent(OUT) :: err
        real(wp), intent(IN) :: fRho, fP, fR, fV, fR0, fpb, fpbdot, alf
        real(wp), intent(IN) :: fntait, fBtait, f_bub_adv_src, f_divu, h
        integer, intent(IN) :: bub_id
        real(wp), intent(IN) :: fmass_v, fmass_g, fbeta_c, fbeta_t, fCson
        real(wp), dimension(4), intent(OUT) :: myR_tmp, myV_tmp, myPb_tmp, myMv_tmp

        real(wp), dimension(4) :: myA_tmp, mydPbdt_tmp, mydMvdt_tmp
        real(wp) :: err_R, err_V

        myPb_tmp(1:4) = fpb
        mydPbdt_tmp(1:4) = fpbdot

        ! Stage 0
        myR_tmp(1) = fR
        myV_tmp(1) = fV
        if (bubbles_lagrange) then
            myPb_tmp(1) = fpb
            myMv_tmp(1) = fmass_v
            call s_advance_EL(myR_tmp(1), myV_tmp(1), myPb_tmp(1), myMv_tmp(1), bub_id, &
                              fmass_g, fbeta_c, fbeta_t, mydPbdt_tmp(1), mydMvdt_tmp(1))
        end if
        myA_tmp(1) = f_rddot(fRho, fP, myR_tmp(1), myV_tmp(1), fR0, &
                             myPb_tmp(1), mydPbdt_tmp(1), alf, fntait, fBtait, &
                             f_bub_adv_src, f_divu, &
                             fCson)

        ! Stage 1
        myR_tmp(2) = myR_tmp(1) + h*myV_tmp(1)
        if (myR_tmp(2) < 0._wp) then
            err = adap_dt_tol + 1._wp; return
        end if
        myV_tmp(2) = myV_tmp(1) + h*myA_tmp(1)
        if (bubbles_lagrange) then
            myPb_tmp(2) = myPb_tmp(1) + h*mydPbdt_tmp(1)
            myMv_tmp(2) = myMv_tmp(1) + h*mydMvdt_tmp(1)
            call s_advance_EL(myR_tmp(2), myV_tmp(2), myPb_tmp(2), myMv_tmp(2), &
                              bub_id, fmass_g, fbeta_c, fbeta_t, mydPbdt_tmp(2), mydMvdt_tmp(2))
        end if
        myA_tmp(2) = f_rddot(fRho, fP, myR_tmp(2), myV_tmp(2), fR0, &
                             myPb_tmp(2), mydPbdt_tmp(2), alf, fntait, fBtait, &
                             f_bub_adv_src, f_divu, &
                             fCson)

        ! Stage 2
        myR_tmp(3) = myR_tmp(1) + (h/4._wp)*(myV_tmp(1) + myV_tmp(2))
        if (myR_tmp(3) < 0._wp) then
            err = adap_dt_tol + 1._wp; return
        end if
        myV_tmp(3) = myV_tmp(1) + (h/4._wp)*(myA_tmp(1) + myA_tmp(2))
        if (bubbles_lagrange) then
            myPb_tmp(3) = myPb_tmp(1) + (h/4._wp)*(mydPbdt_tmp(1) + mydPbdt_tmp(2))
            myMv_tmp(3) = myMv_tmp(1) + (h/4._wp)*(mydMvdt_tmp(1) + mydMvdt_tmp(2))
            call s_advance_EL(myR_tmp(3), myV_tmp(3), myPb_tmp(3), myMv_tmp(3), &
                              bub_id, fmass_g, fbeta_c, fbeta_t, mydPbdt_tmp(3), mydMvdt_tmp(3))
        end if
        myA_tmp(3) = f_rddot(fRho, fP, myR_tmp(3), myV_tmp(3), fR0, &
                             myPb_tmp(3), mydPbdt_tmp(3), alf, fntait, fBtait, &
                             f_bub_adv_src, f_divu, &
                             fCson)

        ! Stage 3
        myR_tmp(4) = myR_tmp(1) + (h/6._wp)*(myV_tmp(1) + myV_tmp(2) + 4._wp*myV_tmp(3))
        if (myR_tmp(4) < 0._wp) then
            err = adap_dt_tol + 1._wp; return
        end if
        myV_tmp(4) = myV_tmp(1) + (h/6._wp)*(myA_tmp(1) + myA_tmp(2) + 4._wp*myA_tmp(3))
        if (bubbles_lagrange) then
            myPb_tmp(4) = myPb_tmp(1) + (h/6._wp)*(mydPbdt_tmp(1) + mydPbdt_tmp(2) + 4._wp*mydPbdt_tmp(3))
            myMv_tmp(4) = myMv_tmp(1) + (h/6._wp)*(mydMvdt_tmp(1) + mydMvdt_tmp(2) + 4._wp*mydMvdt_tmp(3))
            call s_advance_EL(myR_tmp(4), myV_tmp(4), myPb_tmp(4), myMv_tmp(4), &
                              bub_id, fmass_g, fbeta_c, fbeta_t, mydPbdt_tmp(4), mydMvdt_tmp(4))
        end if
        myA_tmp(4) = f_rddot(fRho, fP, myR_tmp(4), myV_tmp(4), fR0, &
                             myPb_tmp(4), mydPbdt_tmp(4), alf, fntait, fBtait, &
                             f_bub_adv_src, f_divu, &
                             fCson)

        ! Estimate error
        err_R = (-5._wp*h/24._wp)*(myV_tmp(2) + myV_tmp(3) - 2._wp*myV_tmp(4)) &
                /max(abs(myR_tmp(1)), abs(myR_tmp(4)))
        err_V = (-5._wp*h/24._wp)*(myA_tmp(2) + myA_tmp(3) - 2._wp*myA_tmp(4)) &
                /max(abs(myV_tmp(1)), abs(myV_tmp(4)))
        ! Error correction for non-oscillating bubbles
        if (max(abs(myV_tmp(1)), abs(myV_tmp(4))) < 1.e-12_wp) then
            err_V = 0._wp
        end if
        if (bubbles_lagrange .and. f_approx_equal(myA_tmp(1), 0._wp) .and. f_approx_equal(myA_tmp(2), 0._wp) .and. &
            f_approx_equal(myA_tmp(3), 0._wp) .and. f_approx_equal(myA_tmp(4), 0._wp)) then
            err_V = 0._wp
        end if
        err = sqrt((err_R**2._wp + err_V**2._wp)/2._wp)

    end subroutine s_advance_substep

    !>  Changes of pressure and vapor mass in the lagrange bubbles.
        !!  @param bub_id Bubble identifier
        !!  @param fmass_g Current mass of gas
        !!  @param fbeta_c Mass transfer coefficient
        !!  @param fbeta_t Heat transfer coefficient
        !!  @param fR_tmp Bubble radius
        !!  @param fV_tmp Bubble radial velocity
        !!  @param fPb_tmp Internal bubble pressure
        !!  @param fMv_tmp Mass of vapor in the bubble
        !!  @param fdPbdt_tmp Rate of change of the internal bubble pressure
        !!  @param fdMvdt_tmp Rate of change of the mass of vapor in the bubble
    elemental subroutine s_advance_EL(fR_tmp, fV_tmp, fPb_tmp, fMv_tmp, bub_id, &
                                      fmass_g, fbeta_c, fbeta_t, fdPbdt_tmp, advance_EL)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(IN) :: fR_tmp, fV_tmp, fPb_tmp, fMv_tmp
        real(wp), intent(IN) :: fmass_g, fbeta_c, fbeta_t
        integer, intent(IN) :: bub_id
        real(wp), intent(INOUT) :: fdPbdt_tmp
        real(wp), intent(out) :: advance_EL
        real(wp) :: fVapFlux, myR_m, mygamma_m

        call s_vflux(fR_tmp, fV_tmp, fPb_tmp, fMv_tmp, bub_id, fVapFlux, fmass_g, fbeta_c, myR_m, mygamma_m)
        fdPbdt_tmp = f_bpres_dot(fVapFlux, fR_tmp, fV_tmp, fPb_tmp, fMv_tmp, bub_id, fbeta_t, myR_m, mygamma_m)
        advance_EL = 4._wp*pi*fR_tmp**2._wp*fVapFlux

    end subroutine s_advance_EL

end module m_bubbles
