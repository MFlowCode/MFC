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


MODULE m_eigen

    ! Dependencies =============================================================
    USE mpi                     ! Message passing interface (MPI) module
    ! ==========================================================================

    IMPLICIT NONE

    CONTAINS 

        SUBROUTINE s_eigs(Amat,Eigs,EVec,N)

            REAL(KIND(0d0)), TARGET, DIMENSION(:,:), INTENT(IN) :: Amat
            REAL(KIND(0d0)), TARGET, DIMENSION(:), INTENT(OUT)  :: EVec, Eigs
            INTEGER, INTENT(IN) :: N

            REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:) ::  VR, A
            REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:)   :: WR, WI, Work
            REAL(KIND(0d0)) ::  DUMMY(1,1), Eigmax

            INTEGER ::  NB2, NMAX, LDA, LDVR, LWORK, INFO,LWKOPT
            INTEGER :: i, j
            INTEGER, PARAMETER :: LWMAX = 1000

            EXTERNAL         DGEEV
            INTRINSIC        DCMPLX

            Nb2=64
            LDA=N
            LDVR=N
            LWORK=LWMAX

            ALLOCATE( VR(LDVR,N), A(N,N) )
            ALLOCATE( WORK(LWORK), WR(N), WI(N) )

            A = Amat

            PRINT*, 'N', N

            CALL DGEEV('N','V',N,A,LDA,WR,WI,DUMMY,1,VR,LDVR,WORK,LWORK,INFO)
            LWKOPT = WORK(1)

            IF (INFO/=0) THEN
                PRINT*, 'Failure in DGEEV.  INFO = ', INFO
            END IF

            EVec(1:N) = VR(1,N:1:-1)
            Eigs = WR(N:1:-1)


            ! ELSE
                ! EigMax = 0d0
                ! DO j = 1,N
                !     IF( WR(j) > EigMax ) THEN
                !         EigMax = WR(j)
                !         ! EVec(:) = VR(j,:)
                !         EVec(:) = VR(:,j)
                !     END IF
                ! END DO
            ! END IF

            ! PRINT*, 'eigenvectors'
            ! PRINT*, VR(1,:)
            ! PRINT*, VR(2,:)
            ! PRINT*, VR(3,:)


            ! print*, 'in eigs'
            ! print*, 'max eig', Eigmax
            ! print*, 'max loc', MAXLOC(WR(:),DIM=1)
            ! print*, 'real eigs: ', WR(:)
            ! print*, 'eigvec associated with largeest eigenvalue', EVec(1:N)

        END SUBROUTINE s_eigs

END MODULE m_eigen
