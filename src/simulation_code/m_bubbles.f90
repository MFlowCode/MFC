!!       __  _______________
!!      /  |/  / ____/ ____/
!!     / /|_/ / /_  / /     
!!    / /  / / __/ / /___   
!!   /_/  /_/_/    \____/   
!!                       
!!  This file is part of MFC.
!!
!!  MFC is the legal property of its developers, whose names 
!!  are listed in the copyright file included with this source 
!!  distribution.
!!
!!  MFC is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published 
!!  by the Free Software Foundation, either version 3 of the license 
!!  or any later version.
!!
!!  MFC is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU General Public License for more details.
!!  
!!  You should have received a copy of the GNU General Public License
!!  along with MFC (LICENSE).  
!!  If not, see <http://www.gnu.org/licenses/>.

!>
!! @file m_bubbles.f90
!! @brief Contains module m_bubbles
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief This module is used to compute the ensemble-averaged bubble dynamic variables
MODULE m_bubbles

    ! Dependencies =============================================================

    USE m_derived_types        !< Definitions of the derived types

    USE m_global_parameters    !< Definitions of the global parameters
    
    USE m_mpi_proxy            !< Message passing interface (MPI) module proxy
    
    USE m_variables_conversion !< State variables type conversion procedures

    ! ==========================================================================

    IMPLICIT NONE

    REAL(KIND(0.D0)) :: chi_vw  !< Bubble wall properties (Ando 2010)
    REAL(KIND(0.D0)) :: k_mw    !< Bubble wall properties (Ando 2010)
    REAL(KIND(0.D0)) :: rho_mw  !< Bubble wall properties (Ando 2010)

    CONTAINS

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
        SUBROUTINE s_compute_bubble_source(idir, q_prim_vf, q_cons_vf, mydivu, &
                                        bub_adv_src, bub_r_src, bub_v_src, bub_p_src, bub_m_src) 
            
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q_prim_vf, q_cons_vf
            TYPE(scalar_field), INTENT(IN) :: mydivu
            INTEGER, INTENT(IN) :: idir

            REAL(KIND(0d0)), DIMENSION(     0:m,0:n,0:p), INTENT(INOUT) ::  bub_adv_src
            REAL(KIND(0d0)), DIMENSION(1:nb,0:m,0:n,0:p), INTENT(INOUT) ::  bub_r_src, &
                                                                            bub_v_src, &
                                                                            bub_p_src, &
                                                                            bub_m_src

            REAL(KIND(0d0)), DIMENSION(0:m,0:n,0:p) :: nbub !< Bubble number density
            
            REAL(KIND(0d0)) ::  tmp1, tmp2, tmp3, tmp4, &
                                c_gas, c2_liquid, &
                                Cpbw, Cpinf,Cpinf_dot, &
                                myH, myHdot, rddot, alf_gas
            
            REAL(KIND(0d0))   :: pb, mv, vflux, pldot, pbdot
            
            REAL(KIND(0d0)) :: n_tait, B_tait, gamma_gas !< EOS parameters
            
            REAL(KIND(0d0)), DIMENSION(nb)  :: Rtmp, Vtmp 
            REAL(KIND(0d0))   :: myR, myV, alf, myP, myRho, R2Vav 
            
            
            REAL(KIND(0d0)), DIMENSION(2) :: Re !< Reynolds number 
            REAL(KIND(0d0)), DIMENSION( num_fluids, &
                                        num_fluids  ) :: We !< Weber number
           
            INTEGER :: j,k,l,q,s !< Loop variables
            INTEGER :: ndirs  !< Number of coordinate directions
            

            
            ndirs = 1; IF (n > 0) ndirs = 2; IF (p > 0) ndirs = 3

            IF (idir == ndirs) THEN
                bub_adv_src = 0.d0; bub_r_src = 0.d0; bub_v_src = 0.d0
                                    bub_p_src = 0.d0; bub_m_src = 0.d0

                ! advection source
                DO j = 0,m; DO k = 0,n; DO l = 0,p
                    ! = 3 \alpha \bar{R^2 V} / \bar{R^3} = 4 pi nbub \bar{R^2 V}
                    DO q = 1,nb
                        Rtmp(q) = q_prim_vf(bub_idx%rs(q))%sf(j,k,l)
                        Vtmp(q) = q_prim_vf(bub_idx%vs(q))%sf(j,k,l)
                    END DO

                    ! Computes n_bub number bubble density
                    CALL s_comp_n_from_prim( q_prim_vf(alf_idx)%sf(j,k,l), &
                                        Rtmp, nbub(j,k,l) )
                      
                    CALL s_quad( (Rtmp**2.d0)*Vtmp, R2Vav )
                    bub_adv_src(j,k,l)  = 4.d0*pi*nbub(j,k,l)*R2Vav
                END DO; END DO; END DO

                ! bubble radius and radial velocity source
                DO q=1,nb; DO j = 0,m; DO k = 0,n; DO l = 0,p
                    bub_r_src(q,j,k,l) = q_cons_vf(bub_idx%vs(q))%sf(j,k,l)
                
                    CALL s_convert_to_mixture_variables( q_cons_vf, myRho, n_tait, B_tait, Re, We, j, k, l )
                    
                    n_tait = 1.d0/n_tait + 1.d0 !make this the usual little 'gamma'

                    IF (num_fluids == 1) THEN
                        gamma_gas  = 1.d0/fluid_pp(num_fluids+1)%gamma + 1.d0
                    ELSE 
                        gamma_gas  = 1.d0/fluid_pp(num_fluids)%gamma + 1.d0
                    END IF

                    myRho= q_prim_vf(1)%sf(j,k,l)
                    myP  = q_prim_vf(E_idx)%sf(j,k,l)
                    alf  = q_prim_vf(alf_idx)%sf(j,k,l)
                    myR  = q_prim_vf(bub_idx%rs(q))%sf(j,k,l)
                    myV  = q_prim_vf(bub_idx%vs(q))%sf(j,k,l)

                    IF (polytropic .NEQV. .TRUE.) THEN
                        pb = q_prim_vf(bub_idx%ps(q))%sf(j,k,l)
                        mv = q_prim_vf(bub_idx%ms(q))%sf(j,k,l)
                        CALL s_bwproperty( pb, q )
                        vflux = f_vflux( myR, myV, mv, q )
                        pbdot = f_bpres_dot( vflux, myR, myV, pb, mv, q )

                        bub_p_src(q,j,k,l) = nbub(j,k,l)*pbdot
                        bub_m_src(q,j,k,l) = nbub(j,k,l)*vflux*4.d0*pi*(myR**2.d0)
                    ELSE
                        pb = 0d0; mv = 0d0; vflux = 0d0; pbdot = 0d0
                    END IF

                   IF (bubble_model == 1) THEN
                        ! Gilmore bubbles
                        Cpinf       = myP - pref
                        Cpbw        = f_cpbw( R0(q), myR, myV, gamma_gas, pb )
                        myH         = f_H       ( Cpbw, Cpinf, n_tait, B_tait )
                        c_gas       = f_cgas    ( Cpinf, n_tait, B_tait, myH )
                        Cpinf_dot   = f_cpinfdot( myRho, myP, alf, n_tait, B_tait, bub_adv_src(j,k,l), mydivu%sf(j,k,l) )
                        myHdot      = f_Hdot    ( Cpbw, Cpinf, Cpinf_dot, n_tait, B_tait, gamma_gas, myR, myV, R0(q), pbdot )
                        rddot       = f_rddot   ( Cpbw, myR, myV, myH, myHdot, c_gas, n_tait, B_tait )
                    ELSE IF (bubble_model == 2) THEN
                        ! Keller-Miksis bubbles
                        Cpinf       = myP
                        Cpbw        = f_cpbw_KM( R0(q), myR, myV, gamma_gas, pb )
                        c_gas = dsqrt( n_tait*(Cpbw+B_tait) / myRho)
                        rddot       = f_rddot_KM( pbdot, gamma_gas, Cpinf, Cpbw, myRho, myR, myV, R0(q), c_gas )
                    END IF

                    bub_v_src(q,j,k,l) = nbub(j,k,l) * rddot
                    
                    IF (isnan(bub_r_src(q,j,k,l))) THEN
                        STOP 'bub_r_src NaN'
                    ELSE IF (isnan(bub_v_src(q,j,k,l))) THEN
                        STOP 'bub_v_src NaN'
                    END IF

                    IF (alf < 1.d-10) THEN
                        bub_adv_src(j,k,l) = 0d0
                        bub_r_src(q,j,k,l) = 0d0
                        bub_v_src(q,j,k,l) = 0d0
                        IF (polytropic .NEQV. .TRUE.) THEN
                           bub_p_src(q,j,k,l) = 0d0
                           bub_m_src(q,j,k,l) = 0d0
                        END IF
                    END IF
                END DO; END DO; END DO; END DO
            END IF
            
        END SUBROUTINE s_compute_bubble_source


        !>  Function that computes that bubble wall pressure for Gilmore bubbles
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fgamma_gas Polytropic gas constant
        !!  @param fpb Internal bubble pressure
        FUNCTION f_cpbw( fR0, fR, fV, fgamma_gas, fpb )

            REAL(KIND(0d0)), INTENT(IN) :: fR0, fR, fV, fgamma_gas, fpb
            
            REAL(KIND(0d0))             :: f_cpbw

            IF (polytropic) THEN
                f_cpbw = (Ca + 2.D0/Web/fR0)*((fR0/fR)**(3.d0*fgamma_gas)) - Ca - 4.D0*Re_inv*fV/fR - 2.D0/(fR*Web)
            ELSE
                f_cpbw = fpb - 1.d0 - 4.D0*Re_inv*fV/fR - 2.D0/(fR*Web)
            END IF

        END FUNCTION f_cpbw


        !>  Function that computes the bubble enthalpy
        !!  @param fCpbw Bubble wall pressure
        !!  @param fCpinf Driving bubble pressure
        !!  @param fntait Tait EOS parameter
        !!  @param fBtait Tait EOS parameter
        FUNCTION f_H( fCpbw, fCpinf, fntait, fBtait )

            REAL(KIND(0d0)), INTENT(IN) :: fCpbw, fCpinf, fntait, fBtait
            
            REAL(KIND(0d0))             :: tmp1,tmp2,tmp3
            REAL(KIND(0d0))             :: f_H
            
            tmp1 = ( fntait-1.D0 )/fntait
            tmp2 = (  fCpbw/(1.D0+fBtait)+1.D0 )**tmp1
            tmp3 = ( fCpinf/(1.D0+fBtait)+1.D0 )**tmp1
            
            f_H = ( tmp2 - tmp3 )*fntait*( 1.D0+fBtait )/( fntait-1.D0 )

        END FUNCTION f_H


        !> Function that computes the sound speed for the bubble
        !! @param fCpinf Driving bubble pressure
        !! @param fntait Tait EOS parameter
        !! @param fBtait Tait EOS parameter
        !! @param fH Bubble enthalpy
        FUNCTION f_cgas( fCpinf, fntait, fBtait, fH )

            REAL(KIND(0d0)), INTENT(IN) :: fCpinf, fntait, fBtait, fH
            
            REAL(KIND(0d0))             :: tmp 
            REAL(KIND(0d0))             :: f_cgas 

            ! get sound speed for Gilmore equations "C" -> c_gas
            tmp = ( fCpinf/(1.D0+fBtait)+1.D0 )**( (fntait-1.D0)/fntait )
            tmp = fntait*( 1.D0+fBtait )*tmp

            f_cgas= DSQRT( tmp+(fntait-1.D0)*fH )

        END FUNCTION f_cgas


        !>  Function that computes the time derivative of the driving pressure
        !!  @param fRho Local liquid density
        !!  @param fP Local pressure
        !!  @param falf Local void fraction
        !!  @param fntait Tait EOS parameter
        !!  @param fBtait Tait EOS parameter
        !!  @param advsrc Advection equation source term
        !!  @param divu Divergence of velocity
        FUNCTION f_cpinfdot( fRho, fP, falf, fntait, fBtait, advsrc, divu )

            REAL(KIND(0d0)), INTENT(IN) :: fRho, fP, falf, fntait, fBtait, advsrc, divu
            
            REAL(KIND(0d0))             :: c2_liquid
            REAL(KIND(0d0))             :: f_cpinfdot
            
            ! get sound speed squared for liquid (only needed for pbdot)
            ! c_l^2 = gam (p+B) / (rho*(1-alf))
            IF (mpp_lim) THEN
                c2_liquid =  fntait*(fP+fBtait)/fRho
            ELSE
                c2_liquid =  fntait*(fP+fBtait)/(fRho*(1.d0-falf))
            END IF

            ! \dot{Cp_inf} = rho sound^2 (alf_src - divu) 
            f_cpinfdot = fRho*c2_liquid*(advsrc - divu)

        END FUNCTION f_cpinfdot


        !>  Function that computes the time derivative of the enthalpy
        !!  @param fCpbw Bubble wall pressure
        !!  @param fCpinf Driving bubble pressure
        !!  @param fCpinf_dot Time derivative of the driving pressure
        !!  @param fntait Tait EOS parameter
        !!  @param fBtait Tait EOS parameter
        !!  @param fgamma_gas Polytropic gas constant
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fpbdot Time derivative of the internal bubble pressure
        FUNCTION f_Hdot( fCpbw, fCpinf, fCpinf_dot, fntait, fBtait, fgamma_gas, fR, fV, fR0, fpbdot )

            REAL(KIND(0d0)), INTENT(IN) :: fCpbw, fCpinf, fCpinf_dot, fntait, fBtait
            REAL(KIND(0d0)), INTENT(IN) :: fgamma_gas, fR, fV, fR0, fpbdot
            
            REAL(KIND(0d0))             :: tmp1, tmp2
            REAL(KIND(0d0))             :: f_Hdot

            IF (polytropic) THEN
                tmp1 = ( fR0/fR )**( 3.d0*fgamma_gas )
                tmp1 = -3.D0*fgamma_gas*( Ca+2d0/Web/fR0 )*tmp1*fV/fR
            ELSE
                tmp1 = fpbdot
            END IF
            tmp2 = ( 2.d0/Web + 4.d0*Re_inv*fV )*fV/(fR**2.d0)

            f_Hdot = &
                ( fCpbw /(1.D0+fBtait)+1.D0 )**( -1.D0/fntait )*( tmp1+tmp2 ) &
              - ( fCpinf/(1.D0+fBtait)+1.D0 )**( -1.D0/fntait )*fCpinf_dot

            ! Hdot = (Cpbw/(1+B) + 1)^(-1/n_tait)*(-3 gam)*(R0/R)^(3gam) V/R
            !f_Hdot = ((fCpbw/(1d0+fBtait)+1.d0)**(-1.d0/fntait))*(-3.d0)*fgamma_gas * &
            !            ( (fR0/fR)**(3.d0*fgamma_gas ))*(fV/fR)

            ! Hdot = Hdot - (Cpinf/(1+B) + 1)^(-1/n_tait) Cpinfdot
            !f_Hdot = f_Hdot - ((fCpinf/(1.d0+fBtait)+1.d0)**(-1.d0/fntait))*fCpinf_dot

        END FUNCTION f_Hdot


        !>  Function that computes the bubble radial acceleration
        !!  @param fCpbw Bubble wall pressure
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fH Current enthalpy
        !!  @param fHdot Current time derivative of the enthalpy
        !!  @param fcgas Current gas sound speed
        !!  @param fntait Tait EOS parameter
        !!  @param fBtait Tait EOS parameter
        FUNCTION f_rddot( fCpbw, fR, fV, fH, fHdot, fcgas, fntait, fBtait )

            REAL(KIND(0d0)), INTENT(IN) :: fCpbw, fR, fV, fH, fHdot
            REAL(KIND(0d0)), INTENT(IN) :: fcgas, fntait, fBtait
            
            REAL(KIND(0d0))             :: tmp1, tmp2, tmp3
            REAL(KIND(0d0))             :: f_rddot

            tmp1 = fV/fcgas
            tmp2 = 1.D0 + 4.D0*Re_inv/fcgas/fR*( fCpbw/(1.D0+fBtait)+1.D0 ) &
                        **( -1.D0/fntait )
            tmp3 = 1.5D0*fV**2d0*( tmp1/3.D0-1.D0 ) + fH*( 1.D0+tmp1 ) &
                        + fR*fHdot*( 1.D0-tmp1 )/fcgas

            f_rddot = tmp3/( fR*(1.d0-tmp1)*tmp2 )

        END FUNCTION f_rddot


        !>  Function that computes the bubble wall pressure for Keller--Miksis bubbles
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fgamma_gas Polytropic gas constant
        !!  @param fpb Internal bubble pressure
        FUNCTION f_cpbw_KM( fR0, fR, fV, fgamma_gas, fpb )

            REAL(KIND(0d0)), INTENT(IN) :: fR0, fR, fV, fgamma_gas, fpb
            REAL(KIND(0d0))             :: f_cpbw_KM
            
            IF (polytropic) THEN
                f_cpbw_KM = (Ca + 2.D0/Web/fR0)*((fR0/fR)**(3.d0*fgamma_gas)) - &
                    Ca - 4.D0*Re_inv*fV/fR - 2.D0/(fR*Web)
            ELSE
                f_cpbw_KM = fpb - 4.D0*Re_inv*fV/fR - 2.D0/(fR*Web)
            END IF

        END FUNCTION f_cpbw_KM



        !>  Function that computes the bubble radial acceleration for Keller--Miksis bubbles
        !!  @param fpbdot Time-derivative of internal bubble pressure
        !!  @param fgamma_gas Polytropic gas constant
        !!  @param fCp Driving pressure
        !!  @param fCpbw Bubble wall pressure
        !!  @param fRho Current density
        !!  @param fR Current bubble radius
        !!  @param fV Current bubble velocity
        !!  @param fR0 Equilibrium bubble radius
        !!  @param fC Current sound speed
        FUNCTION f_rddot_KM( fpbdot, fgamma_gas, fCp, fCpbw, fRho, fR, fV, fR0, fC )

            REAL(KIND(0d0)), INTENT(IN) :: fpbdot, fgamma_gas, fCp, fCpbw
            REAL(KIND(0d0)), INTENT(IN) :: fRho, fR, fV, fR0, fC

            REAL(KIND(0d0))             :: tmp1, tmp2, cdot_star
            REAL(KIND(0d0))             :: f_rddot_KM

            IF (polytropic) THEN
                tmp1 = ( fR0/fR )**( 3.d0*fgamma_gas )
                tmp1 = -3.D0*fgamma_gas*( Ca+2d0/Web/fR0 )*tmp1*fV/fR
            ELSE
                tmp1 = fpbdot
            END IF
            cdot_star = tmp1 + ( 2.d0/Web + 4.d0*Re_inv*fV )*fV/(fR**2.d0)
            tmp1 = fV/fC
            tmp2 =  1.5D0*fV**2d0*( tmp1/3.D0-1.D0 ) +                      &
                    (fCpbw - fCp)/fRho * (1.d0 - tmp1) +      &
                    cdot_star * fR/fRho/fC
                  
            f_rddot_KM = tmp2/( fR*(1.d0-tmp1) + 4.d0*Re_inv/(fRho*fC) )

        END FUNCTION f_rddot_KM


        !>  Subroutine that computes bubble wall properties for vapor bubbles
        !>  @param pb Internal bubble pressure
        !>  @param iR0 Current bubble size index
        SUBROUTINE s_bwproperty( pb, iR0 )

            REAL(KIND(0.D0)), INTENT(IN) :: pb
            INTEGER, INTENT(IN) :: iR0 

            REAL(KIND(0.D0)) :: x_vw

            ! mass fraction of vapor
            chi_vw = 1.D0/( 1.D0+R_v/R_n*(pb/pv-1.D0) )
            ! mole fraction of vapor & thermal conductivity of gas mixture
            x_vw = M_n*chi_vw/( M_v+(M_n-M_v)*chi_vw )
            k_mw = x_vw*k_v(iR0)/( x_vw+(1.D0-x_vw)*phi_vn ) &
                + ( 1.D0-x_vw )*k_n(iR0)/( x_vw*phi_nv+1.D0-x_vw )
            ! gas mixture density
            rho_mw = pv/( chi_vw*R_v*Tw )

        END SUBROUTINE s_bwproperty



        !>  Function that computes the vapour flux 
        !!  @param fR Current bubble radius 
        !!  @param fV Current bubble velocity
        !!  @param fmass_v Current mass of vapour
        !!  @param iR0 Bubble size index
        FUNCTION f_vflux( fR,fV,fmass_v, iR0 )

            REAL(KIND(0.D0)), INTENT(IN) :: fR
            REAL(KIND(0.D0)), INTENT(IN) :: fV
            REAL(KIND(0.D0)), INTENT(IN) :: fmass_v
            INTEGER, INTENT(IN) :: iR0

            REAL(KIND(0.D0)) :: chi_bar
            REAL(KIND(0.D0)) :: grad_chi
            REAL(KIND(0.D0)) :: f_vflux

            IF ( thermal==3 ) THEN !transfer
                ! constant transfer model
                chi_bar = fmass_v/( fmass_v+mass_n0(iR0) )
                grad_chi = -Re_trans_c(iR0)*( chi_bar-chi_vw )
                f_vflux = rho_mw*grad_chi/Pe_c/( 1.D0-chi_vw )/fR
            ELSE
                ! polytropic
                f_vflux = pv*fV/( R_v*Tw )
            END IF

        END FUNCTION f_vflux


        !>  Function that computes the time derivative of 
        !!  the internal bubble pressure
        !!  @param fvflux Vapour flux
        !!  @param fR Current bubble radius 
        !!  @param fV Current bubble velocity
        !!  @param fpb Current internal bubble pressure
        !!  @param fmass_v Current mass of vapour
        !!  @param iR0 Bubble size index
        FUNCTION f_bpres_dot( fvflux,fR,fV,fpb,fmass_v,iR0 )

            REAL(KIND(0.D0)), INTENT(IN) :: fvflux
            REAL(KIND(0.D0)), INTENT(IN) :: fR
            REAL(KIND(0.D0)), INTENT(IN) :: fV
            REAL(KIND(0.D0)), INTENT(IN) :: fpb
            REAL(KIND(0.D0)), INTENT(IN) :: fmass_v
            INTEGER, INTENT(IN) :: iR0

            REAL(KIND(0.D0)) :: T_bar
            REAL(KIND(0.D0)) :: grad_T
            REAL(KIND(0.D0)) :: tmp1, tmp2
            REAL(KIND(0.D0)) :: f_bpres_dot

            IF ( thermal==3 ) THEN
                T_bar = Tw*( fpb/pb0(iR0) )*( fR/R0(iR0) )**3 &
                    * ( mass_n0(iR0)+mass_v0(iR0) )/( mass_n0(iR0)+fmass_v )
                grad_T = -Re_trans_T(iR0)*( T_bar-Tw )
                f_bpres_dot = 3.D0*gamma_m*( -fV*fpb+fvflux*R_v*Tw &
                   + pb0(iR0)*k_mw*grad_T/Pe_T(iR0)/fR )/fR
            ELSE
                f_bpres_dot = -3.D0*gamma_m*fV/fR*( fpb-pv )
            END IF
            
    END FUNCTION f_bpres_dot


END MODULE m_bubbles
