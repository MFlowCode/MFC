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


    CONTAINS

        SUBROUTINE s_mom_inv( q_prim_vf, momsp, moms3d, is1, is2, is3 ) 
            
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q_prim_vf
            TYPE(scalar_field), DIMENSION(nmomsp), INTENT(INOUT) :: momsp
            TYPE(scalar_field), DIMENSION(0:nnode,0:nnode,nb), INTENT(INOUT) :: moms3d
            TYPE(bounds_info), INTENT(IN) :: is1,is2,is3

            REAL(KIND(0d0)), DIMENSION(nmom) :: moms
            REAL(KIND(0d0)), DIMENSION(nb,nnode) :: wght, abscX, abscY
            REAL(KIND(0d0)), DIMENSION(nterms,0:2,0:2,nb) :: mom3d_terms

            REAL(KIND(0d0)) :: pres, gam, moms_cond

            INTEGER :: j,k,l,q,r,s !< Loop variables
            INTEGER :: id1,id2,id3
            INTEGER :: i1,i2


            IF (num_fluids == 1) THEN
                gam  = 1.d0/fluid_pp(num_fluids+1)%gamma + 1.d0
            ELSE 
                gam  = 1.d0/fluid_pp(num_fluids)%gamma + 1.d0
            END IF

            ! TODO: should these loops be over ix%beg,ix%end or is current form sufficient?
            DO id3 = is3%beg, is3%end; DO id2 = is2%beg, is2%end; DO id1 = is1%beg, is1%end
                pres = q_prim_vf(E_idx)%sf(id1,id2,id3)
                DO q = 1,nb
                    DO r = 1,nmom
                        moms(r) = q_prim_vf(bub_idx%moms(q,r))%sf(id1,id2,id3)
                    END DO
                    CALL s_chyqmom(pres,moms,wght(q,:),abscX(q,:),abscY(q,:))

                    DO j = 1,nterms
                        i1 = 1; i2 = 0
                        moms_cond = f_quad2D(abscX(q,:),abscY(q,:),wght(q,:),momrhs(i1,i2,q,j,1),momrhs(i1,i2,q,j,2))
                        mom3d_terms(j,i1,i2,q) = f_get_coeff(j,i1,i2,q,pres) * (R0(q)**(momrhs(i1,i2,q,j,3)-q)) * moms_cond
                        i1 = 0; i2 = 1
                        moms_cond = f_quad2D(abscX(q,:),abscY(q,:),wght(q,:),momrhs(i1,i2,q,j,1),momrhs(i1,i2,q,j,2))
                        mom3d_terms(j,i1,i2,q) = f_get_coeff(j,i1,i2,q,pres) * (R0(q)**(momrhs(i1,i2,q,j,3)-q)) * moms_cond
                        i1 = 2; i2 = 0
                        moms_cond = f_quad2D(abscX(q,:),abscY(q,:),wght(q,:),momrhs(i1,i2,q,j,1),momrhs(i1,i2,q,j,2))
                        mom3d_terms(j,i1,i2,q) = f_get_coeff(j,i1,i2,q,pres) * (R0(q)**(momrhs(i1,i2,q,j,3)-q)) * moms_cond
                        i1 = 1; i2 = 1
                        moms_cond = f_quad2D(abscX(q,:),abscY(q,:),wght(q,:),momrhs(i1,i2,q,j,1),momrhs(i1,i2,q,j,2))
                        mom3d_terms(j,i1,i2,q) = f_get_coeff(j,i1,i2,q,pres) * (R0(q)**(momrhs(i1,i2,q,j,3)-q)) * moms_cond
                        i1 = 0; i2 = 2
                        moms_cond = f_quad2D(abscX(q,:),abscY(q,:),wght(q,:),momrhs(i1,i2,q,j,1),momrhs(i1,i2,q,j,2))
                        mom3d_terms(j,i1,i2,q) = f_get_coeff(j,i1,i2,q,pres) * (R0(q)**(momrhs(i1,i2,q,j,3)-q)) * moms_cond
                    END DO
                    i1=1; i2=0;
                    moms3d(i1,i2,q)%sf(id1,id2,id3) = SUM( mom3d_terms(:,i1,i2,q) )
                    i1=0; i2=1;
                    moms3d(i1,i2,q)%sf(id1,id2,id3) = SUM( mom3d_terms(:,i1,i2,q) )
                    i1=2; i2=0;
                    moms3d(i1,i2,q)%sf(id1,id2,id3) = SUM( mom3d_terms(:,i1,i2,q) )
                    i1=1; i2=1;
                    moms3d(i1,i2,q)%sf(id1,id2,id3) = SUM( mom3d_terms(:,i1,i2,q) )
                    i1=0; i2=2;
                    moms3d(i1,i2,q)%sf(id1,id2,id3) = SUM( mom3d_terms(:,i1,i2,q) )

                    ! print*, 'mom3d (1,0,', q,')', moms3d(1,0,q)%sf(id1,id2,id3)
                    ! print*, 'mom3d (0,1,', q,')', moms3d(0,1,q)%sf(id1,id2,id3)
                    ! print*, 'mom3d (2,0,', q,')', moms3d(2,0,q)%sf(id1,id2,id3)
                    ! print*, 'mom3d (1,1,', q,')', moms3d(1,1,q)%sf(id1,id2,id3)
                    ! print*, 'mom3d (0,2,', q,')', moms3d(0,2,q)%sf(id1,id2,id3)
                END DO

                momsp(1)%sf(id1,id2,id3) = f_quad(abscX,abscY,wght,3d0,0d0,0d0)
                momsp(2)%sf(id1,id2,id3) = f_quad(abscX,abscY,wght,2d0,1d0,0d0)
                momsp(3)%sf(id1,id2,id3) = f_quad(abscX,abscY,wght,3d0,2d0,0d0)
                momsp(4)%sf(id1,id2,id3) = f_quad(abscX,abscY,wght,3d0*(1d0-gam),0d0,3*gam)

                ! print*, 'momsp(1) ', momsp(1)%sf(id1,id2,id3)
                ! print*, 'momsp(2) ', momsp(2)%sf(id1,id2,id3)
                ! print*, 'momsp(3) ', momsp(3)%sf(id1,id2,id3)
                ! print*, 'momsp(4) ', momsp(4)%sf(id1,id2,id3)
                ! stop

            END DO; END DO; END DO


            print*, 'finished qbmm'


            ! CALL s_get_moms( weights, abscX, abscY, momsp, mom3d, pres, is1, is2, is3)
            ! CALL s_get_momsp ( weights(:,:,j,k,l),abscX(:,:,j,k,l),abscY(:,:,j,k,l),     momsp(1:4)%sf(j,k,l) )
            ! CALL s_get_momrhs( weights(:,:,j,k,l),abscX(:,:,j,k,l),abscY(:,:,j,k,l), mom3d(:,:,:)%sf(j,k,l), pres )

        END SUBROUTINE s_mom_inv

        ! SUBROUTINE s_get_moms( wght, abscX, abscY, moms3d, momsp, pres, is1, is2, is3) 

        !     REAL(KIND(0d0)), DIMENSION(nb,nnode), INTENT(IN) :: wght, abscX, abscY
        !     TYPE(scalar_field), DIMENSION(nmom,nmom,nb), INTENT(OUT) :: moms3d
        !     TYPE(scalar_field), DIMENSION(nmomsp), INTENT(OUT) :: momsp
        !     REAL(KIND(0d0)), INTENT(IN) :: pres
        !     TYPE(bounds_info), INTENT(IN) :: is1,is2,is3

        !     REAL(KIND(0d0)), DIMENSION(nterms,nmom,nmom,nb) :: mom3d_terms
        !     REAL(KIND(0d0)), DIMENSION(nterms,0:2,0:2,nb) :: moms_cond

        !     INTEGER :: i,j,k,i1,i2

        !     DO id3 = is3%beg, is3%end; DO id2 = is2%beg, is2%end; DO id1 = is1%beg, is1%end

        !         DO k = 1,nb
        !             DO j = 1,nterms
        !                 i1 = 1; i2 = 0
        !                 moms_cond(j,i1,i2,k) = f_quad2D(abscX(k,:),abscY(k,:),wght(k,:),momrhs(i1,i2,k,j,1),momrhs(i1,i2,k,j,2))
        !                 mom3d_terms(j,i1,i2,k) = f_get_coeff(j,i1,i2,k,pres) * (R0(k)**(momrhs(i1,i2,k,j,3)-k)) * moms_cond(j,i1,i2,k)
        !                 i1 = 0; i2 = 1
        !                 moms_cond(j,i1,i2,k) = f_quad2D(abscX(k,:),abscY(k,:),wght(k,:),momrhs(i1,i2,k,j,1),momrhs(i1,i2,k,j,2))
        !                 mom3d_terms(j,i1,i2,k) = f_get_coeff(j,i1,i2,k,pres) * (R0(k)**(momrhs(i1,i2,k,j,3)-k)) * moms_cond(j,i1,i2,k)
        !                 i1 = 2; i2 = 0
        !                 moms_cond(j,i1,i2,k) = f_quad2D(abscX(k,:),abscY(k,:),wght(k,:),momrhs(i1,i2,k,j,1),momrhs(i1,i2,k,j,2))
        !                 mom3d_terms(j,i1,i2,k) = f_get_coeff(j,i1,i2,k,pres) * (R0(k)**(momrhs(i1,i2,k,j,3)-k)) * moms_cond(j,i1,i2,k)
        !                 i1 = 1; i2 = 1
        !                 moms_cond(j,i1,i2,k) = f_quad2D(abscX(k,:),abscY(k,:),wght(k,:),momrhs(i1,i2,k,j,1),momrhs(i1,i2,k,j,2))
        !                 mom3d_terms(j,i1,i2,k) = f_get_coeff(j,i1,i2,k,pres) * (R0(k)**(momrhs(i1,i2,k,j,3)-k)) * moms_cond(j,i1,i2,k)
        !                 i1 = 0; i2 = 2
        !                 moms_cond(j,i1,i2,k) = f_quad2D(abscX(k,:),abscY(k,:),wght(k,:),momrhs(i1,i2,k,j,1),momrhs(i1,i2,k,j,2))
        !                 mom3d_terms(j,i1,i2,k) = f_get_coeff(j,i1,i2,k,pres) * (R0(k)**(momrhs(i1,i2,k,j,3)-k)) * moms_cond(j,i1,i2,k)
        !             END DO
        !             i1=1;i2=0;
        !             moms3d(i1,i2,k)%sf(id1,id3,id3) = SUM( mom3d_terms(:,i1,i2,k) )
        !             i1=0;i2=1;
        !             moms3d(i1,i2,k)%sf(id1,id3,id3) = SUM( mom3d_terms(:,i1,i2,k) )
        !             i1=2;i2=0;
        !             moms3d(i1,i2,k)%sf(id1,id3,id3) = SUM( mom3d_terms(:,i1,i2,k) )
        !             i1=1;i2=1;
        !             moms3d(i1,i2,k)%sf(id1,id3,id3) = SUM( mom3d_terms(:,i1,i2,k) )
        !             i1=0;i2=2;
        !             moms3d(i1,i2,k)%sf(id1,id3,id3) = SUM( mom3d_terms(:,i1,i2,k) )
        !         END DO

        !         momsp(1)%sf(id1,id3,id3) = f_quad(abscX,abscY,wght,3d0,0d0,0d0)
        !         momsp(2)%sf(id1,id3,id3) = f_quad(abscX,abscY,wght,2d0,1d0,0d0)
        !         momsp(3)%sf(id1,id3,id3) = f_quad(abscX,abscY,wght,3d0,2d0,0d0)
        !         momsp(4)%sf(id1,id3,id3) = f_quad(abscX,abscY,wght,3d0*(1d0-gam),0d0,3*gam)
        !     END DO; END DO; END DO

        ! END SUBROUTINE s_get_momrhs


        SUBROUTINE s_chyqmom(pres,momin,wght,abscX,abscY)
            
            REAL(KIND(0d0)), DIMENSION(nnode), INTENT(INOUT) :: wght,abscX,abscY
            REAL(KIND(0d0)), DIMENSION(nmom), INTENT(IN) :: momin
            REAL(KIND(0d0)), INTENT(IN) :: pres
    
            REAL(KIND(0d0)), DIMENSION(0:nmom,0:nmom) :: moms
            REAL(KIND(0d0)), DIMENSION(3) :: M1, M3
            REAL(KIND(0d0)), DIMENSION(2) :: myrho, myrho3, up, up3, Vf
            REAL(KIND(0d0)) :: bu, bv, d20, d11, d02, c20, c11, c02
            REAL(KIND(0d0)) :: mu2avg, mu2, vp21, vp22, rho21, rho22

            ! print*, 'moms: ', momin(:)

            moms(1,0) = momin(1)
            moms(0,1) = momin(2)
            moms(2,0) = momin(3)
            moms(1,1) = momin(4)
            moms(0,2) = momin(5)

            bu  = moms(1,0) 
            bv  = moms(0,1) 
            d20 = moms(2,0) 
            d11 = moms(1,1) 
            d02 = moms(0,2) 
            c20 = d20 - bu**2d0; 
            c11 = d11 - bu*bv;
            c02 = d02 - bv**2d0; 
            M1 = (/ 1d0, 0d0, c20 /)
            CALL s_hyqmom(myrho,up,M1)
            Vf = c11*up/c20
            mu2avg = c02 - sum(myrho(:)*(Vf(:)**2d0))
            mu2avg = maxval( (/ mu2avg, 0d0 /) )
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
            fup(1) = bu - dsqrt(c2)
            fup(2) = bu + dsqrt(c2) 

        END SUBROUTINE s_hyqmom


        FUNCTION f_quad( abscX,abscY,wght,q,r,s )
            REAL(KIND(0.D0)), DIMENSION(nb,nnode), INTENT(IN) :: abscX, abscY, wght
            REAL(KIND(0.D0)), INTENT(IN) :: q,r,s
            REAL(KIND(0.D0)) :: f_quad_RV, f_quad
            INTEGER :: i,j,k

            f_quad = 0d0
            DO i = 1,nb
                f_quad_RV = sum( wght(i,:)*(abscX(i,:)**q)*(abscY(i,:)**r) )
                f_quad = f_quad + weight(i)*(R0(i)**s)*f_quad_RV
            END DO
        END FUNCTION f_quad


        FUNCTION f_quad2D( abscX,abscY,wght,q,r )
            REAL(KIND(0.D0)), DIMENSION(nnode), INTENT(IN) :: abscX, abscY, wght
            REAL(KIND(0.D0)), INTENT(IN) :: q,r
            REAL(KIND(0.D0)) :: f_quad_RV, f_quad2D
            INTEGER :: i,j,k

            f_quad2D = sum( wght(:)*(abscX(:)**q)*(abscY(:)**r) )
        END FUNCTION f_quad2D


        FUNCTION f_get_coeff( term,i1,i2,i3,pres )
            INTEGER, INTENT(IN) :: term,i1,i2,i3
            REAL(KIND(0.D0)), INTENT(IN) :: pres
            REAL(KIND(0.D0)) :: f_get_coeff

            IF (term == 1) THEN
                f_get_coeff = -i2*pres
            ELSEIF (term == 2) THEN
                f_get_coeff = -3d0*i2/2d0
            ELSEIF (term == 3) THEN
                f_get_coeff = i2
            ELSEIF (term == 4) THEN
                f_get_coeff = i1
            END IF

        END FUNCTION f_get_coeff


END MODULE m_qbmm
