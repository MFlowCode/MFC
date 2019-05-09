MODULE m_bubbles

    USE m_derived_types        ! Definitions of the derived types
    USE m_global_parameters    ! Definitions of the global parameters
    USE m_mpi_proxy            ! Message passing interface (MPI) module proxy
    USE m_variables_conversion ! State variables type conversion procedures
    
    IMPLICIT NONE

    ! bubble wall properties
    REAL(KIND(0.D0)) :: chi_vw
    REAL(KIND(0.D0)) :: k_mw
    REAL(KIND(0.D0)) :: rho_mw

    CONTAINS

        SUBROUTINE s_compute_bubble_source(idir, q_prim_vf, q_cons_vf, mydivu, &
                                        bub_adv_src, bub_r_src, bub_v_src, bub_p_src, bub_m_src) 
        ! Description: The purpose of this procedure is to compute the source term
        ! that are needed for the bubble modeling

            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q_prim_vf, q_cons_vf
            TYPE(scalar_field), intent(in) :: mydivu
            INTEGER, INTENT(IN) :: idir

            REAL(kind(0d0)), dimension(     0:m,0:n,0:p), intent(inout) ::  bub_adv_src
            REAL(kind(0d0)), dimension(1:nb,0:m,0:n,0:p), intent(inout) ::  bub_r_src, &
                                                                            bub_v_src, &
                                                                            bub_p_src, &
                                                                            bub_m_src
            REAL(KIND(0d0)), DIMENSION(0:m,0:n,0:p) :: nbub
            
            REAL(KIND(0d0)) ::  tmp1, tmp2, tmp3, tmp4, &
                                c_gas, c2_liquid, &
                                Cpbw, Cpinf,Cpinf_dot, &
                                myH, myHdot, rddot, alf_gas

            REAL(KIND(0d0)) :: n_tait, B_tait, gamma_gas
            
            REAL(KIND(0d0)), DIMENSION(nb)  :: Rtmp, Vtmp 
            REAL(KIND(0d0))   :: myR, myV, alf, myP, myRho, R2Vav 
            
            REAL(KIND(0d0)), DIMENSION(2) :: Re
            REAL(KIND(0d0)), DIMENSION( num_fluids, &
                                        num_fluids  ) :: We
            
            INTEGER :: j,k,l,q,s,ndirs
            
            REAL(KIND(0d0))   :: pb, mv, vflux, pldot, pbdot
            
            ndirs = 1; if (n > 0) ndirs = 2; if (p > 0) ndirs = 3

            if (idir == ndirs) then
                bub_adv_src = 0.d0; bub_r_src = 0.d0; bub_v_src = 0.d0
                                    bub_p_src = 0.d0; bub_m_src = 0.d0

                ! advection source
                do j = 0,m; do k = 0,n; do l = 0,p
                    ! = 3 \alpha \bar{R^2 V} / \bar{R^3} = 4 pi nbub \bar{R^2 V}
                    do q = 1,nb
                        Rtmp(q) = q_prim_vf(bub_idx%rs(q))%sf(j,k,l)
                        Vtmp(q) = q_prim_vf(bub_idx%vs(q))%sf(j,k,l)
                    end do
                    call s_comp_n_from_prim( q_prim_vf(alf_idx)%sf(j,k,l), &
                                        Rtmp, nbub(j,k,l) ) !get n_bub
                      
                    call s_quad( (Rtmp**2.d0)*Vtmp, R2Vav )
                    bub_adv_src(j,k,l)  = b_coeff(j,k,l)*(4.d0*pi*nbub(j,k,l)*R2Vav)
                end do; end do; end do

                ! bubble radius and radial velocity source
                do q=1,nb; do j = 0,m; do k = 0,n; do l = 0,p
                    bub_r_src(q,j,k,l) = b_coeff(j,k,l)*q_cons_vf(bub_idx%vs(q))%sf(j,k,l)
                
                    CALL s_convert_to_mixture_variables( q_cons_vf, myRho, n_tait, B_tait, Re, We, j, k, l )
                    
                    n_tait = 1.d0/n_tait + 1.d0 !make this the usual little 'gamma'

                    if (num_fluids == 1) then
                        gamma_gas  = 1.d0/fluid_pp(num_fluids+1)%gamma + 1.d0
                    else 
                        gamma_gas  = 1.d0/fluid_pp(num_fluids)%gamma + 1.d0
                    end if

                    myRho= q_prim_vf(1)%sf(j,k,l)
                    myP  = q_prim_vf(E_idx)%sf(j,k,l)
                    alf  = q_prim_vf(alf_idx)%sf(j,k,l)
                    myR  = q_prim_vf(bub_idx%rs(q))%sf(j,k,l)
                    myV  = q_prim_vf(bub_idx%vs(q))%sf(j,k,l)

                    if (polytropic .neqv. .TRUE.) then
                        pb = q_prim_vf(bub_idx%ps(q))%sf(j,k,l)
                        mv = q_prim_vf(bub_idx%ms(q))%sf(j,k,l)
                        call s_bwproperty( pb, q )
                        vflux = f_vflux( myR, myV, mv, q )
                        pbdot = f_bpres_dot( vflux, myR, myV, pb, mv, q )

                        bub_p_src(q,j,k,l) = nbub(j,k,l)*pbdot
                        bub_m_src(q,j,k,l) = nbub(j,k,l)*vflux*4.d0*pi*(myR**2.d0)
                    else
                        pb = 0d0; mv = 0d0; vflux = 0d0; pbdot = 0d0
                    end if

                   if (bubble_model == 1) then
                        !Gilbert bubbles
                        Cpinf       = myP - pref
                        Cpbw        = f_cpbw( R0(q), myR, myV, gamma_gas, pb )
                        myH         = f_H       ( Cpbw, Cpinf, n_tait, B_tait )
                        c_gas       = f_cgas    ( Cpinf, n_tait, B_tait, myH )
                        Cpinf_dot   = f_cpinfdot( myRho, myP, alf, n_tait, B_tait, bub_adv_src(j,k,l), mydivu%sf(j,k,l) )
                        myHdot      = f_Hdot    ( Cpbw, Cpinf, Cpinf_dot, n_tait, B_tait, gamma_gas, myR, myV, R0(q), pbdot )
                        rddot       = f_rddot   ( Cpbw, myR, myV, myH, myHdot, c_gas, n_tait, B_tait )
                    else if (bubble_model == 2) then
                        !Keller-Miksis bubbles
                        Cpinf       = myP
                        Cpbw        = f_cpbw_KM( R0(q), myR, myV, gamma_gas, pb )

                        !actually liquid sound speed
                        !c_gas = dsqrt( n_tait*(pref*Cpbw+pref)/myRho )
                        !c_gas = dsqrt( n_tait*(Cpbw+B_tait)/myRho )
                        c_gas = dsqrt( n_tait*(Cpbw+B_tait) / myRho)
                        rddot       = f_rddot_KM( pbdot, gamma_gas, Cpinf, Cpbw, myRho, myR, myV, R0(q), c_gas, n_tait, B_tait )
                    end if
                    bub_v_src(q,j,k,l) = b_coeff(j,k,l) * nbub(j,k,l) * rddot
                    if (isnan(bub_r_src(q,j,k,l))) then
                        stop 'bub_r_src NaN'
                    else if (isnan(bub_v_src(q,j,k,l))) then
                        print*,  q_prim_vf(adv_idx%beg)%sf(j,k,l) 
                        print*,  q_prim_vf(adv_idx%beg+1)%sf(j,k,l) 
                        print*,  q_prim_vf(adv_idx%beg+2)%sf(j,k,l) 
                        print*, n_tait, B_tait, myRho
                        print*, cpinf, cpbw, myh, c_gas, cpinf_dot, myHdot, rddot, nbub(j,k,l)
                        stop 'bub_v_src NaN'
                    end if

                    if  ( (alf         < 1.d-10 ) & !.or. &
                        !.and. ( rddot < 0d0 )    & 
                        !  (nbub(j,k,l) < 1.d-10) .or. &
                        !  (myR         < 1.d-10)      &
                        ) then
                        bub_adv_src(j,k,l) = 0d0
                        bub_r_src(q,j,k,l) = 0d0
                        bub_v_src(q,j,k,l) = 0d0
                        if (polytropic .neqv. .TRUE.) then
                           bub_p_src(q,j,k,l) = 0d0
                           bub_m_src(q,j,k,l) = 0d0
                        end if
                    end if
                end do; end do; end do; end do
            end if
           
            !bub_adv_src = 0d0; bub_r_src = 0d0; bub_v_src = 0d0
            !bub_p_src = 0d0; bub_m_src = 0d0
            
            if (proc_rank == 0) then
                print*, 'Source terms:'
                print*, 'max adv src', maxval(abs(bub_adv_src(:,:,:)))
                print*, 'max r src', maxval(abs(bub_r_src(:,:,:,:)))
                print*, 'max v src', maxval(abs(bub_v_src(:,:,:,:)))
                print*, 'max p src', maxval(abs(bub_p_src(:,:,:,:)))
                print*, 'max m src', maxval(abs(bub_m_src(:,:,:,:)))
            end if
        END SUBROUTINE s_compute_bubble_source


        function f_cpbw( fR0, fR, fV, fgamma_gas, fpb )
            real(kind(0d0)), intent(in) :: fR0, fR, fV, fgamma_gas, fpb
            real(kind(0d0))             :: f_cpbw

            if (polytropic) then
                f_cpbw = (Ca + 2.D0/Web/fR0)*((fR0/fR)**(3.d0*fgamma_gas)) - Ca - 4.D0*Re_inv*fV/fR - 2.D0/(fR*Web)
            else
                f_cpbw = fpb - 1.d0 - 4.D0*Re_inv*fV/fR - 2.D0/(fR*Web)
            end if
        end function f_cpbw


        function f_H( fCpbw, fCpinf, fntait, fBtait )
            real(kind(0d0)), intent(in) :: fCpbw, fCpinf, fntait, fBtait
            real(kind(0d0))             :: tmp1,tmp2,tmp3
            real(kind(0d0))             :: f_H
            
            tmp1 = ( fntait-1.D0 )/fntait
            tmp2 = (  fCpbw/(1.D0+fBtait)+1.D0 )**tmp1
            tmp3 = ( fCpinf/(1.D0+fBtait)+1.D0 )**tmp1
            
            f_H = ( tmp2 - tmp3 )*fntait*( 1.D0+fBtait )/( fntait-1.D0 )
        end function f_H


        function f_cgas( fCpinf, fntait, fBtait, fH )
            real(kind(0d0)), intent(in) :: fCpinf, fntait, fBtait, fH
            real(kind(0d0))             :: tmp 
            real(kind(0d0))             :: f_cgas 

            ! get sound speed for Gilmore equations "C" -> c_gas
            tmp = ( fCpinf/(1.D0+fBtait)+1.D0 )**( (fntait-1.D0)/fntait )
            tmp = fntait*( 1.D0+fBtait )*tmp

            f_cgas= DSQRT( tmp+(fntait-1.D0)*fH )
        end function f_cgas


        function f_cpinfdot( fRho, fP, falf, fntait, fBtait, advsrc, divu )
            real(kind(0d0)), intent(in) :: fRho, fP, falf, fntait, fBtait, advsrc, divu
            real(kind(0d0))             :: c2_liquid
            real(kind(0d0))             :: f_cpinfdot
            
            ! get sound speed squared for liquid (only needed for pbdot)
            ! c_l^2 = gam (p+B) / (rho*(1-alf))
            if (mpp_lim) then
                c2_liquid =  fntait*(fP+fBtait)/fRho
            else
                c2_liquid =  fntait*(fP+fBtait)/(fRho*(1.d0-falf))
            end if

            ! \dot{Cp_inf} = rho sound^2 (alf_src - divu) 
            f_cpinfdot = fRho*c2_liquid*(advsrc - divu)
        end function f_cpinfdot


        function f_Hdot( fCpbw, fCpinf, fCpinf_dot, fntait, fBtait, fgamma_gas, fR, fV, fR0, fpbdot )
            real(kind(0d0)), intent(in) :: fCpbw, fCpinf, fCpinf_dot, fntait, fBtait, &
                                           fgamma_gas, fR, fV, fR0, fpbdot
            real(kind(0d0))             :: tmp1, tmp2
            real(kind(0d0))             :: f_Hdot

            if (polytropic) then
                tmp1 = ( fR0/fR )**( 3.d0*fgamma_gas )
                tmp1 = -3.D0*fgamma_gas*( Ca+2d0/Web/fR0 )*tmp1*fV/fR
            else
                tmp1 = fpbdot
            end if
            tmp2 = ( 2.d0/Web + 4.d0*Re_inv*fV )*fV/(fR**2.d0)

            f_Hdot = &
                ( fCpbw /(1.D0+fBtait)+1.D0 )**( -1.D0/fntait )*( tmp1+tmp2 ) &
              - ( fCpinf/(1.D0+fBtait)+1.D0 )**( -1.D0/fntait )*fCpinf_dot

            ! Hdot = (Cpbw/(1+B) + 1)^(-1/n_tait)*(-3 gam)*(R0/R)^(3gam) V/R
            !f_Hdot = ((fCpbw/(1d0+fBtait)+1.d0)**(-1.d0/fntait))*(-3.d0)*fgamma_gas * &
            !            ( (fR0/fR)**(3.d0*fgamma_gas ))*(fV/fR)

            ! Hdot = Hdot - (Cpinf/(1+B) + 1)^(-1/n_tait) Cpinfdot
            !f_Hdot = f_Hdot - ((fCpinf/(1.d0+fBtait)+1.d0)**(-1.d0/fntait))*fCpinf_dot
        end function f_Hdot


        function f_rddot( fCpbw, fR, fV, fH, fHdot, fcgas, fntait, fBtait )
            real(kind(0d0)), intent(in) :: fCpbw, fR, fV, fH, fHdot, fcgas, fntait, fBtait
            real(kind(0d0))             :: tmp1, tmp2, tmp3
            real(kind(0d0))             :: f_rddot

            tmp1 = fV/fcgas
            tmp2 = 1.D0 + 4.D0*Re_inv/fcgas/fR*( fCpbw/(1.D0+fBtait)+1.D0 ) &
                        **( -1.D0/fntait )
            tmp3 = 1.5D0*fV**2d0*( tmp1/3.D0-1.D0 ) + fH*( 1.D0+tmp1 ) &
                        + fR*fHdot*( 1.D0-tmp1 )/fcgas

            f_rddot = tmp3/( fR*(1.d0-tmp1)*tmp2 )
        end function f_rddot


        function f_cpbw_KM( fR0, fR, fV, fgamma_gas, fpb )
            real(kind(0d0)), intent(in) :: fR0, fR, fV, fgamma_gas, fpb
            real(kind(0d0))             :: f_cpbw_KM
            
            if (polytropic) then
                !f_cpbw = (fR0/fR)**(3.d0*fgamma_gas) - 1.d0
                f_cpbw_KM = (Ca + 2.D0/Web/fR0)*((fR0/fR)**(3.d0*fgamma_gas)) - Ca - 4.D0*Re_inv*fV/fR - 2.D0/(fR*Web)
            else
                f_cpbw_KM = fpb - 4.D0*Re_inv*fV/fR - 2.D0/(fR*Web)
            end if
        end function f_cpbw_KM


        function f_rddot_KM( fpbdot, fgamma_gas, fCp, fCpbw, fRho, fR, fV, fR0, fC, fntait, fBtait )
            real(kind(0d0)), intent(in) :: fpbdot, fgamma_gas, fCp, fCpbw, fRho, fR, fV, fR0, fC, fntait, fBtait
            real(kind(0d0))             :: tmp1, tmp2, cdot_star
            real(kind(0d0))             :: f_rddot_KM

            if (polytropic) then
                tmp1 = ( fR0/fR )**( 3.d0*fgamma_gas )
                tmp1 = -3.D0*fgamma_gas*( Ca+2d0/Web/fR0 )*tmp1*fV/fR
            else
                tmp1 = fpbdot
            end if
            cdot_star = tmp1 + ( 2.d0/Web + 4.d0*Re_inv*fV )*fV/(fR**2.d0)
            tmp1 = fV/fC
            tmp2 =  1.5D0*fV**2d0*( tmp1/3.D0-1.D0 ) +                      &
                    (fCpbw - fCp)/fRho * (1.d0 - tmp1) +      &
                    cdot_star * fR/fRho/fC
                  
            f_rddot_KM = tmp2/( fR*(1.d0-tmp1) + 4.d0*Re_inv/(fRho*fC) )
        end function f_rddot_KM


        SUBROUTINE s_bwproperty( pb, iR0 )
            REAL(KIND(0.D0)), INTENT(IN) :: pb
            INTEGER, INTENT(IN) :: iR0 
            REAL(KIND(0.D0)) :: x_vw

            !vapor='y' only

            ! mass fraction of vapor
            chi_vw = 1.D0/( 1.D0+R_v/R_n*(pb/pv-1.D0) )
            ! mole fraction of vapor & thermal conductivity of gas mixture
            x_vw = M_n*chi_vw/( M_v+(M_n-M_v)*chi_vw )
            k_mw = x_vw*k_v(iR0)/( x_vw+(1.D0-x_vw)*phi_vn ) &
                + ( 1.D0-x_vw )*k_n(iR0)/( x_vw*phi_nv+1.D0-x_vw )
            ! gas mixture density
            rho_mw = pv/( chi_vw*R_v*Tw )
        END SUBROUTINE s_bwproperty


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
                ! irrerevant procedure (polytropic)
                f_vflux = pv*fV/( R_v*Tw )
            END IF
        END FUNCTION f_vflux


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

            !IF ( thermal=='transfer'.AND.model=='Preston' ) THEN
            IF ( thermal==3 ) THEN
                ! constant transfer model
                T_bar = Tw*( fpb/pb0(iR0) )*( fR/R0(iR0) )**3 &
                    * ( mass_n0(iR0)+mass_v0(iR0) )/( mass_n0(iR0)+fmass_v )
                grad_T = -Re_trans_T(iR0)*( T_bar-Tw )
                f_bpres_dot = 3.D0*gamma_m*( -fV*fpb+fvflux*R_v*Tw &
                   + pb0(iR0)*k_mw*grad_T/Pe_T(iR0)/fR )/fR
            !ELSE IF ( thermal=='transfer'.AND.model=='Sugiyama' ) THEN
                ! only for gas bubble
            !    T_bar = Tw*( fpb/pb0(iR0) )*( fR/R0(iR0) )**3
            !    tmp1 = -fV*pb - pb0(iR0)*k_mw/Pe_T(iR0)/R &
            !        * ( Re_trans_T(iR0)*(T_bar-Tw) &
            !        + Im_trans_T(iR0)*T_bar/omegaN(iR0)*3.D0*fV/fR )
            !    tmp2 = fR/3.D0/gamma_m &
            !        + pb0(iR0)*k_mw*Im_trans_T(iR0)*T_bar &
            !        / Pe_T(iR0)/R/omegaN(iR0)/fpb
            !    f_bpres_dot = tmp1/tmp2
            ELSE
                f_bpres_dot = -3.D0*gamma_m*fV/fR*( fpb-pv )
            END IF
    END FUNCTION f_bpres_dot


    SUBROUTINE s_phase_transfer(q_cons_vf) ! ----------------
        TYPE(scalar_field), DIMENSION(sys_size), INTENT(INOUT) :: q_cons_vf
        
        REAL(KIND(0d0)), DIMENSION(Nb)  ::       nRtmp
        REAL(KIND(0d0))                 ::       nbub
        REAL(KIND(0d0))                 ::       alf
        !loop iterators
        integer :: i,j,k,l

        return 
        do j = 0,m; do k = 0,n; do l = 0,p
            if (q_cons_vf(alf_idx)%sf(j,k,l) > 0.1d0) then 

                alf = q_cons_vf(alf_idx)%sf(j,k,l)
                ! n = sqrt( 4pi/(3 alpha) * (nR)**3 )
                do i = 1,nb
                    nRtmp(i) = q_cons_vf(bub_idx%rs(i))%sf(j,k,l)
                end do
                call s_comp_n_from_cons( alf, nRtmp, nbub)   


                !transfer bubble void fraction to another component
                q_cons_vf(alf_idx-1)%sf(j,k,l) = alf


                alf  = 1d-8; 
                nbub = (3.d0/(4.d0*pi)) * alf/1.d0
                q_cons_vf(alf_idx)%sf(j,k,l) = alf   !get rid of bubble void fraction

                do i = 1,nb
                    q_cons_vf(bub_idx%rs(i))%sf(j,k,l) = nbub * 1.d0 !&
                        !* ((3.d0/(4.d0*pi) * alf/nbub)**(1.d0/3.d0))
                    q_cons_vf(bub_idx%vs(i))%sf(j,k,l) = 0d0
                end do
            end if
        end do; end do; end do

    END SUBROUTINE s_phase_transfer

END MODULE m_bubbles
