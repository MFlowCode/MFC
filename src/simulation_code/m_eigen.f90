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
