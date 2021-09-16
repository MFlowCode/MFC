!!       __  _______________
!!      /  |/  / ____/ ____/
!!     / /|_/ / /_  / /     
!!    / /  / / __/ / /___   
!!   /_/  /_/_/    \____/   
!!                       
!!  This file is part of MFC.
!!
!! Copyright 2021
!!
!! Permission is hereby granted, free of charge, to any person 
!! obtaining a copy of this software and associated documentation 
!! files (the "Software"), to deal in the Software without 
!! restriction, including without limitation the rights to use, 
!! copy, modify, merge, publish, distribute, sublicense, 
!! and/or sell copies of the Software, and to permit persons 
!! to whom the Software is furnished to do so, subject to the 
!! following conditions:
!!
!! The above copyright notice and this permission notice shall 
!! be included in all copies or substantial portions of the Software.
!!
!! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
!! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
!! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
!! IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
!! FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION 
!! OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
!! CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 
!! THE SOFTWARE.

!>
!! @file m_qbmm.f90
!! @brief Contains module m_qbmm
!! @author S. Bryngelson
!! @version 1.0
!! @date MAY 28, 2020

!> @brief This module is used to compute moment inversion via qbmm
MODULE m_qbmm

    ! Dependencies =============================================================

    USE m_derived_types        !< Definitions of the derived types

    USE m_global_parameters    !< Definitions of the global parameters
    
    USE m_mpi_proxy            !< Message passing interface (MPI) module proxy
    
    USE m_variables_conversion !< State variables type conversion procedures

    ! ==========================================================================

    IMPLICIT NONE


    REAL(KIND(0d0)), PARAMETER :: verysmall = 1.D-12
    REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:,:,:,:)  :: momrhs
    INTEGER         :: nterms !< Number of rhs terms in each moment transport equations

    CONTAINS

        SUBROUTINE s_initialize_qbmm_module()

            INTEGER :: i1,i2,q

            IF (bubble_model==2) THEN 
                ! Keller-Miksis without viscosity/surface tension
                nterms = 12
            ELSE IF (bubble_model==3) THEN
                ! Rayleigh-Plesset with viscosity/surface tension
                nterms = 6
            END IF

            ALLOCATE( momrhs(0:2,0:2,nb,nterms,3) )
            momrhs = 0d0

            ! Assigns the required RHS moments for moment transport equations
            ! The rhs%(:,3) is only to be used for R0 quadrature, not for computing X/Y indices
            DO q = 1,nb
                DO i1 = 0,2; DO i2 = 0,2
                    IF ( (i1+i2)<=2 ) THEN
                        IF (bubble_model==3) THEN
                            momrhs(i1,i2,q,1,1) = -1.d0 + i1
                            momrhs(i1,i2,q,1,2) = -1.d0 + i2
                            momrhs(i1,i2,q,1,3) = 0d0

                            momrhs(i1,i2,q,2,1) = -1.d0 + i1
                            momrhs(i1,i2,q,2,2) =  1.d0 + i2
                            momrhs(i1,i2,q,2,3) = 0d0

                            momrhs(i1,i2,q,3,1) = -1.d0 + i1 - 3.d0*gam
                            momrhs(i1,i2,q,3,2) = -1.d0 + i2
                            momrhs(i1,i2,q,3,3) = 3.d0*gam

                            momrhs(i1,i2,q,4,1) = -1.d0 + i1
                            momrhs(i1,i2,q,4,2) =  1.d0 + i2
                            momrhs(i1,i2,q,4,3) = 0d0

                            IF (Re_inv .NE. dflt_real) THEN
                                ! add viscosity
                                momrhs(i1,i2,q,5,1) = -2.d0 + i1
                                momrhs(i1,i2,q,5,2) = i2
                                momrhs(i1,i2,q,5,3) = 0d0
                            END IF

                            IF (Web .NE. dflt_real) THEN
                                ! add surface tension
                                momrhs(i1,i2,q,6,1) = -2.d0 + i1
                                momrhs(i1,i2,q,6,2) = -1.d0 + i2
                                momrhs(i1,i2,q,6,3) = 0d0
                            END IF
                        ELSE IF (bubble_model==2) THEN
                            ! KM with approximation of 1/(1-V/C) = 1+V/C
                            momrhs(i1,i2,q,1,1) = -1d0 + i1
                            momrhs(i1,i2,q,1,2) =  1d0 + i2
                            momrhs(i1,i2,q,1,3) =  0d0

                            momrhs(i1,i2,q,2,1) = -1d0 + i1
                            momrhs(i1,i2,q,2,2) =  2d0 + i2
                            momrhs(i1,i2,q,2,3) =  0d0

                            momrhs(i1,i2,q,3,1) = -1d0 + i1
                            momrhs(i1,i2,q,3,2) =  3d0 + i2
                            momrhs(i1,i2,q,3,3) =  0d0

                            momrhs(i1,i2,q,4,1) = -1d0 + i1
                            momrhs(i1,i2,q,4,2) = -1d0 + i2
                            momrhs(i1,i2,q,4,3) =  0d0

                            momrhs(i1,i2,q,5,1) = -1d0 + i1
                            momrhs(i1,i2,q,5,2) =  i2
                            momrhs(i1,i2,q,5,3) =  0d0

                            momrhs(i1,i2,q,6,1) = -1d0 + i1
                            momrhs(i1,i2,q,6,2) =  1d0 + i2
                            momrhs(i1,i2,q,6,3) =  0d0

                            momrhs(i1,i2,q,7,1) = -1d0 + i1 - 3d0*gam
                            momrhs(i1,i2,q,7,2) = -1d0 + i2
                            momrhs(i1,i2,q,7,3) =  3d0*gam

                            momrhs(i1,i2,q,8,1) = -1d0 + i1 - 3d0*gam
                            momrhs(i1,i2,q,8,2) =  i2
                            momrhs(i1,i2,q,8,3) =  3d0*gam

                            momrhs(i1,i2,q,9,1) = -1d0 + i1 - 3d0*gam
                            momrhs(i1,i2,q,9,2) =  1d0 + i2
                            momrhs(i1,i2,q,9,3) =  3d0*gam

                            momrhs(i1,i2,q,10,1) = -1d0 + i1 - 3d0*gam
                            momrhs(i1,i2,q,10,2) =  i2
                            momrhs(i1,i2,q,10,3) =  3d0*gam

                            momrhs(i1,i2,q,11,1) = -1d0 + i1 - 3d0*gam
                            momrhs(i1,i2,q,11,2) =  1d0 + i2
                            momrhs(i1,i2,q,11,3) =  3d0*gam

                            momrhs(i1,i2,q,12,1) = -1d0 + i1
                            momrhs(i1,i2,q,12,2) =  1d0 + i2
                            momrhs(i1,i2,q,12,3) =  0d0
                        END IF
                    END IF
                END DO; END DO 
            END DO

        END SUBROUTINE s_initialize_qbmm_module


        SUBROUTINE s_coeff( pres,rho,c,coeffs )

            REAL(KIND(0.D0)), INTENT(IN) :: pres, rho, c
            REAL(KIND(0.D0)), DIMENSION(nterms,0:2,0:2), INTENT(OUT) :: coeffs
            INTEGER :: i1,i2

            coeffs = 0d0

            DO i1 = 0,2; DO i2 = 0,2
                IF ((i1+i2)<=2) THEN
                    IF (bubble_model==3) THEN
                        ! RPE
                        coeffs(1,i1,i2) = -1d0*i2*pres/rho
                        coeffs(2,i1,i2) = -3d0*i2/2d0
                        coeffs(3,i1,i2) = i2/rho
                        coeffs(4,i1,i2) = i1
                        IF (Re_inv/=dflt_real) coeffs(5,i1,i2) = -4d0*i2*Re_inv/rho
                        IF (  Web /=dflt_real) coeffs(6,i1,i2) = -2d0*i2/Web/rho
                    ELSE IF (bubble_model==2) THEN
                        ! KM with approximation of 1/(1-V/C) = 1+V/C
                        coeffs(1,i1,i2)  = -3d0*i2/2d0
                        coeffs(2,i1,i2)  = -i2/c
                        coeffs(3,i1,i2)  = i2/(2d0*c*c)
                        coeffs(4,i1,i2)  = -i2*pres/rho
                        coeffs(5,i1,i2)  = -2d0*i2*pres/(c*rho)
                        coeffs(6,i1,i2)  = -i2*pres/(c*c*rho)
                        coeffs(7,i1,i2)  = i2/rho
                        coeffs(8,i1,i2)  = 2d0*i2/(c*rho)
                        coeffs(9,i1,i2)  = i2/(c*c*rho)
                        coeffs(10,i1,i2) = -3d0*i2*gam/(c*rho)
                        coeffs(11,i1,i2) = -3d0*i2*gam/(c*c*rho)
                        coeffs(12,i1,i2) = i1
                    END IF
                END IF
            END DO; END DO

        END SUBROUTINE s_coeff


        SUBROUTINE s_mom_inv( q_prim_vf, momsp, moms3d, is1, is2, is3 ) 
            
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q_prim_vf
            TYPE(scalar_field), DIMENSION(nmomsp), INTENT(INOUT) :: momsp
            TYPE(scalar_field), DIMENSION(0:2,0:2,nb), INTENT(INOUT) :: moms3d
            TYPE(bounds_info), INTENT(IN) :: is1,is2,is3

            REAL(KIND(0d0)), DIMENSION(nmom) :: moms
            REAL(KIND(0d0)), DIMENSION(nb) :: Rvec
            REAL(KIND(0d0)), DIMENSION(nb,nnode) :: wght, abscX, abscY
            REAL(KIND(0d0)), DIMENSION(nterms,0:2,0:2) :: mom3d_terms, coeff
            REAL(KIND(0d0)) :: pres, rho, nbub, c, alf
            REAL(KIND(0d0)) :: n_tait, B_tait

            INTEGER :: j,k,l,q,r,s !< Loop variables
            INTEGER :: id1,id2,id3
            INTEGER :: i1,i2

            DO id3 = is3%beg, is3%end; DO id2 = is2%beg, is2%end; DO id1 = is1%beg, is1%end

                alf = q_prim_vf(alf_idx)%sf(id1,id2,id3)
                pres = q_prim_vf(E_idx)%sf(id1,id2,id3)
                rho  = q_prim_vf(cont_idx%beg)%sf(id1,id2,id3)
                IF (bubble_model==2) THEN
                    n_tait = fluid_pp(1)%gamma
                    n_tait = 1.d0/n_tait + 1.d0 !make this the usual little 'gamma'
                    B_tait = fluid_pp(1)%pi_inf
                    c = n_tait*(pres+B_tait)/(rho*(1.d0-alf))
                    IF (c > 0.d0) THEN
                        c = DSQRT(c)
                    ELSE
                        c = sgm_eps
                    END IF                    
                END IF

                CALL s_coeff(pres,rho,c,coeff)

                ! SHB: Manually adjusted pressure here for no-coupling case
                ! pres = 1d0/0.3d0

                IF (alf > small_alf) THEN
                    DO q = 1,nb
                        Rvec(q) = q_prim_vf(bub_idx%rs(q))%sf(id1,id2,id3)
                    END DO
                    CALL s_comp_n_from_prim( alf, Rvec, nbub )

                    DO q = 1,nb
                        DO r = 1,nmom
                            moms(r) = q_prim_vf(bub_idx%moms(q,r))%sf(id1,id2,id3)
                        END DO

                        ! IF(id1==0) THEN
                        !     PRINT*, 'pres: ', pres
                        !     PRINT*, 'nb : ', nbub
                        !     PRINT*, 'alf: ', alf
                        !     DO s = 1,nmom
                        !         PRINT*, 'mom: ', moms(s)
                        !     END DO
                        ! END IF

                        CALL s_chyqmom(moms,wght(q,:),abscX(q,:),abscY(q,:))

                        DO j = 1,nterms
                            DO i1 = 0,2; DO i2 = 0,2
                                IF ( (i1+i2)<=2 ) THEN
                                    mom3d_terms(j,i1,i2) = coeff(j,i1,i2)*(R0(q)**momrhs(i1,i2,q,j,3))  &
                                       * f_quad2D(abscX(q,:),abscY(q,:),wght(q,:),momrhs(i1,i2,q,j,:))
                                END IF
                            END DO; END DO
                        END DO

                        DO i1 = 0,2; DO i2 = 0,2
                            IF ( (i1+i2)<=2 ) THEN
                                moms3d(i1,i2,q)%sf(id1,id2,id3) = nbub*SUM( mom3d_terms(:,i1,i2) )
                                ! IF (moms3d(i1,i2,q)%sf(id1,id2,id3) .NE. moms3d(i1,i2,q)%sf(id1,id2,id3)) THEN
                                !     PRINT*, 'nan in mom3d', i1,i2,id1
                                !     PRINT*, 'nbu: ', nbub
                                !     PRINT*, 'alf: ', alf
                                !     PRINT*, 'moms: ', moms(:)
                                !     CALL s_mpi_abort()
                                ! END IF
                            END IF
                        END DO; END DO
                    END DO

                    momsp(1)%sf(id1,id2,id3) = f_quad(abscX,abscY,wght,3d0,0d0,0d0)
                    momsp(2)%sf(id1,id2,id3) = 4.d0*pi*nbub*f_quad(abscX,abscY,wght,2d0,1d0,0d0)
                    momsp(3)%sf(id1,id2,id3) = f_quad(abscX,abscY,wght,3d0,2d0,0d0)
                    IF (ABS(gam-1.d0) <= 1.d-4) THEN
                        ! Gam \approx 1, don't risk imaginary quadrature
                        momsp(4)%sf(id1,id2,id3) = 1.d0
                    ELSE 
                        momsp(4)%sf(id1,id2,id3) = f_quad(abscX,abscY,wght,3d0*(1d0-gam),0d0,3d0*gam)
                    END IF

                    DO i1 = 1,4
                        IF (momsp(i1)%sf(id1,id2,id3) /= momsp(i1)%sf(id1,id2,id3)) THEN
                            PRINT*, 'NaN in sp moment', i1, 'location',id1,id2,id3
                            PRINT*, 'Rs', Rvec(:)
                            PRINT*, 'alpha', alf
                            PRINT*, 'nbub', nbub
                            PRINT*, 'abscX', abscX(:,:)
                            PRINT*, 'abscY', abscY(:,:)
                            PRINT*, 'wght' , wght(:,:)
                            CALL s_mpi_abort()
                        END IF
                    END DO
                ELSE
                    DO q = 1,nb
                        DO i1 = 0,2; DO i2 = 0,2
                            moms3d(i1,i2,q)%sf(id1,id2,id3) = 0d0
                        END DO; END DO
                    END DO
                    momsp(1)%sf(id1,id2,id3) = 0d0 
                    momsp(2)%sf(id1,id2,id3) = 0d0 
                    momsp(3)%sf(id1,id2,id3) = 0d0 
                    momsp(4)%sf(id1,id2,id3) = 0d0 
                END IF

            END DO; END DO; END DO

        END SUBROUTINE s_mom_inv


        SUBROUTINE s_chyqmom(momin,wght,abscX,abscY)
            
            REAL(KIND(0d0)), DIMENSION(nnode), INTENT(INOUT) :: wght,abscX,abscY
            REAL(KIND(0d0)), DIMENSION(nmom), INTENT(IN) :: momin
    
            REAL(KIND(0d0)), DIMENSION(0:nmom,0:nmom) :: moms
            REAL(KIND(0d0)), DIMENSION(3) :: M1, M3
            REAL(KIND(0d0)), DIMENSION(2) :: myrho, myrho3, up, up3, Vf
            REAL(KIND(0d0)) :: bu, bv, d20, d11, d02, c20, c11, c02
            REAL(KIND(0d0)) :: mu2avg, mu2, vp21, vp22, rho21, rho22

            moms(0,0) = momin(1)
            moms(1,0) = momin(2)
            moms(0,1) = momin(3)
            moms(2,0) = momin(4)
            moms(1,1) = momin(5)
            moms(0,2) = momin(6)

            bu  = moms(1,0)/moms(0,0)
            bv  = moms(0,1)/moms(0,0)
            d20 = moms(2,0)/moms(0,0)
            d11 = moms(1,1)/moms(0,0)
            d02 = moms(0,2)/moms(0,0)
            c20 = d20 - bu**2d0; 
            c11 = d11 - bu*bv;
            c02 = d02 - bv**2d0; 

            M1 = (/ 1d0, 0d0, c20 /)
            CALL s_hyqmom(myrho,up,M1)
            Vf = c11*up/c20

            mu2avg = c02 - SUM(myrho(:)*(Vf(:)**2d0))
            mu2avg = MAXVAL( (/ mu2avg, 0d0 /) )
            mu2 = mu2avg 
            M3 = (/ 1d0, 0d0, mu2 /)
            CALL s_hyqmom(myrho3,up3,M3)

            vp21 = up3(1)
            vp22 = up3(2)
            rho21 = myrho3(1)
            rho22 = myrho3(2)

            wght(1) = myrho(1)*rho21 
            wght(2) = myrho(1)*rho22 
            wght(3) = myrho(2)*rho21 
            wght(4) = myrho(2)*rho22 
            wght = moms(0,0)*wght

            abscX(1) = up(1) 
            abscX(2) = up(1) 
            abscX(3) = up(2) 
            abscX(4) = up(2) 
            abscX = bu + abscX

            abscY(1) = Vf(1) + vp21 
            abscY(2) = Vf(1) + vp22 
            abscY(3) = Vf(2) + vp21 
            abscY(4) = Vf(2) + vp22 
            abscY = bv + abscY

        END SUBROUTINE s_chyqmom


        SUBROUTINE s_hyqmom(frho,fup,fmom)
            
            REAL(KIND(0d0)), DIMENSION(2), INTENT(INOUT) :: frho, fup
            REAL(KIND(0d0)), DIMENSION(3), INTENT(IN) :: fmom
            REAL(KIND(0d0)) :: bu, d2, c2

            bu = fmom(2)/fmom(1)
            d2 = fmom(3)/fmom(1)
            c2 = d2 - bu**2d0 
            frho(1) = fmom(1)/2d0; 
            frho(2) = fmom(1)/2d0; 
            c2 = MAXVAL( (/ c2, verysmall /) )
            fup(1) = bu - DSQRT(c2)
            fup(2) = bu + DSQRT(c2) 

        END SUBROUTINE s_hyqmom

        FUNCTION f_quad( abscX,abscY,wght,q,r,s )
            REAL(KIND(0.D0)), DIMENSION(nb,nnode), INTENT(IN) :: abscX, abscY, wght
            REAL(KIND(0.D0)), INTENT(IN) :: q,r,s
            REAL(KIND(0.D0)) :: f_quad_RV, f_quad
            INTEGER :: i

            f_quad = 0d0
            DO i = 1,nb
                f_quad_RV = SUM( wght(i,:)*(abscX(i,:)**q)*(abscY(i,:)**r) )
                f_quad = f_quad + weight(i)*(R0(i)**s)*f_quad_RV
            END DO

        END FUNCTION f_quad


        FUNCTION f_quad2D( abscX,abscY,wght,pow )
            REAL(KIND(0.D0)), DIMENSION(nnode), INTENT(IN) :: abscX, abscY, wght
            REAL(KIND(0.D0)), DIMENSION(3), INTENT(IN) :: pow
            REAL(KIND(0.D0)) :: f_quad2D

            f_quad2D = SUM( wght(:)*(abscX(:)**pow(1))*(abscY(:)**pow(2)) )
        END FUNCTION f_quad2D



END MODULE m_qbmm
