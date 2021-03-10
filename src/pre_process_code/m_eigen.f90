module m_eigen

    ! Dependencies =============================================================
    use mpi                     ! Message passing interface (MPI) module
    ! ==========================================================================

    implicit none

contains

    subroutine s_eigs(Amat, Eigs, EVec, N)

        real(kind(0d0)), target, dimension(:, :), intent(IN) :: Amat
        real(kind(0d0)), target, dimension(:), intent(OUT)  :: EVec, Eigs
        integer, intent(IN) :: N

        real(kind(0d0)), allocatable, dimension(:, :) ::  VR, A
        real(kind(0d0)), allocatable, dimension(:)   :: WR, WI, Work
        real(kind(0d0)) ::  DUMMY(1, 1), Eigmax

        integer ::  NB2, NMAX, LDA, LDVR, LWORK, INFO, LWKOPT
        integer :: i, j
        integer, parameter :: LWMAX = 1000

        external DGEEV
        intrinsic DCMPLX

        Nb2 = 64
        LDA = N
        LDVR = N
        LWORK = LWMAX

        allocate (VR(LDVR, N), A(N, N))
        allocate (WORK(LWORK), WR(N), WI(N))

        A = Amat

        print *, 'N', N

        call DGEEV('N', 'V', N, A, LDA, WR, WI, DUMMY, 1, VR, LDVR, WORK, LWORK, INFO)
        LWKOPT = WORK(1)

        if (INFO /= 0) then
            print *, 'Failure in DGEEV.  INFO = ', INFO
        end if

        EVec(1:N) = VR(1, N:1:-1)
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

    end subroutine s_eigs

end module m_eigen
