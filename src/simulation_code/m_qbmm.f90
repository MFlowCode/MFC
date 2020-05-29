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

    REAL(KIND(0.D0)) :: chi_vw  !< Bubble wall properties (Ando 2010)

    CONTAINS

        SUBROUTINE s_mom_inv( q_prim_vf, weights, abscX, abscY ) 
            
            REAL(KIND(0d0)), DIMENSION(nb,nnode,0:m,0:n,0:p), INTENT(INOUT) :: weights, abscX, abscY
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q_prim_vf

            REAL(KIND(0d0)) :: pres
            REAL(KIND(0d0)), DIMENSION(nmom) :: moms
            INTEGER :: j,k,l,q,r,s !< Loop variables

            ! TODO: should these loops be over ix%beg,ix%end or is current form sufficient?
            DO j = 0,m; DO k = 0,n; DO l = 0,p
                pres = q_prim_vf(E_idx)%sf(j,k,l)
                DO q = 1,nb
                    DO r = 1,nmom
                        moms(r) = q_prim_vf(bub_idx%moms(q,r))%sf(j,k,l)
                    END DO
                    CALL s_chyqmom(pres,moms,weights(q,:,j,k,l),abscX(q,:,j,k,l),abscY(q,:,j,k,l))

                    ! print*, weights(q,:,j,k,l)
                    ! print*, abscX(q,:,j,k,l)
                    ! print*, abscY(q,:,j,k,l)

                END DO
            END DO; END DO; END DO
            print*, 'got all weights for step 1'
            
        END SUBROUTINE s_mom_inv

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

        ! FUNCTION f_bpres_dot( fvflux,iR0 )
        !     REAL(KIND(0.D0)), INTENT(IN) :: fvflux
        !     INTEGER, INTENT(IN) :: iR0
        !     REAL(KIND(0.D0)) :: f_bpres_dot
        !     f_bpres_dot = -3.D0
        ! END FUNCTION f_bpres_dot

END MODULE m_qbmm
