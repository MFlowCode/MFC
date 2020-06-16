MODULE m_eigen

    ! Dependencies =============================================================
    USE m_derived_types         ! Definitions of the derived types
    
    USE m_global_parameters     ! Global parameters for the code
    
    USE m_mpi_proxy             ! Message passing interface (MPI) module proxy
    
    USE mpi                     ! Message passing interface (MPI) module
    ! ==========================================================================

    IMPLICIT NONE

    CONTAINS 

        SUBROUTINE evsolve(Amat,Eigs,EVec,myn)

            REAL(KIND(0d0)), TARGET, DIMENSION(:,:), INTENT(IN) :: Amat
            REAL(KIND(0d0)), TARGET, DIMENSION(:), INTENT(OUT)  :: Eigs, EVec
            INTEGER, INTENT(IN) :: myn

            REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:) ::  Vr, Ain
            REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:)   :: Wi, Wr, Work
            REAL(KIND(0d0)) ::  DUMMY(1,1)

            INTEGER, PARAMETER :: NIN=5 , NOUT=6
            INTEGER ::  NB2, NMAX, LDA, LDVR, LWORK, INFO,LWKOPT
            INTEGER :: i, j

            EXTERNAL         DGEEV
            INTRINSIC        DCMPLX

            Nb2=64
            Nmax=myn*4
            Lda=Nmax
            Ldvr=Nmax
            Lwork=(2+Nb2)*Nmax

            ALLOCATE( Vr(Ldvr,Nmax), Ain(Nmax,Nmax) )
            ALLOCATE( Wi(Nmax), Work(Lwork), Wr(Nmax) )

            Ain = Amat

            CALL DGEEV('No left vectors','Vectors (right)',myn, &
                Ain,LDA,WR,WI,DUMMY,1,VR,LDVR,WORK,LWORK,INFO)
            LWKOPT = WORK(1)
            IF (INFO==0) THEN
                DO j = 1,myn
                    Eigs(j) = WR(j)
                    EVec(j) = VR(1,j)
                    ! DO i = 1,2
                        ! EVec(i,j) = VR(i,j)
                    ! END DO
                END DO
            ELSE
                PRINT*, 'Failure in DGEEV.  INFO = ', INFO
            END IF

        END SUBROUTINE evsolve

END MODULE m_eigen
