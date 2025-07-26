!>
!! @file m_lapack_example.f90
!! @brief Contains module m_lapack_example

!> @brief This module demonstrates the use of LAPACK in MFC post_process.
!!        It provides example routines that show how to use LAPACK for
!!        common linear algebra operations like solving linear systems.
module m_lapack_example

    use m_global_parameters     !< Global parameters for the code
    use m_mpi_proxy             !< Message passing interface (MPI) module proxy

    implicit none

    private; public :: s_lapack_example_solve_linear_system, &
 s_lapack_example_eigenvalues

contains

    !> @brief Example subroutine demonstrating LAPACK usage for solving
    !!        a linear system Ax = b using DGESV/SGESV
    !!        This routine shows how to use LAPACK with MFC's precision system
    impure subroutine s_lapack_example_solve_linear_system()

        ! Local variables for the linear system
        integer, parameter :: n = 3  ! Size of the system
        real(wp), dimension(n, n) :: A  ! Coefficient matrix
        real(wp), dimension(n) :: b     ! Right-hand side vector
        real(wp), dimension(n) :: x     ! Solution vector

        ! LAPACK variables
        integer, dimension(n) :: ipiv   ! Pivot indices
        integer :: info                 ! Return status
        integer, parameter :: nrhs = 1  ! Number of right-hand sides

        ! Only run on the root process to avoid duplicate output
        if (proc_rank /= 0) return

        ! Set up a simple 3x3 linear system: Ax = b
        ! Example:
        !   2x + y + z = 8
        !   x + 3y + z = 10
        !   x + y + 4z = 16
        A(1, :) = [2.0_wp, 1.0_wp, 1.0_wp]
        A(2, :) = [1.0_wp, 3.0_wp, 1.0_wp]
        A(3, :) = [1.0_wp, 1.0_wp, 4.0_wp]

        b = [8.0_wp, 10.0_wp, 16.0_wp]

        print *, "LAPACK Linear System Solver Example"
        print *, "Solving the system Ax = b where:"
        print *, "A = [2 1 1; 1 3 1; 1 1 4]"
        print *, "b = [8; 10; 16]"

        ! Copy b to x (LAPACK will overwrite the right-hand side with solution)
        x = b

        ! Call appropriate LAPACK routine based on precision
#ifdef MFC_SINGLE_PRECISION
        call sgesv(n, nrhs, A, n, ipiv, x, n, info)
        print *, "Using LAPACK (SGESV)"
#else
        call dgesv(n, nrhs, A, n, ipiv, x, n, info)
        print *, "Using LAPACK (DGESV)"
#endif

        ! Check for success
        if (info == 0) then
            print *, "Linear system solved successfully!"
            print *, "Solution: x = [", x(1), ", ", x(2), ", ", x(3), "]"
            print *, "Expected: x = [1, 2, 3]"
        else if (info < 0) then
            print *, "LAPACK error: argument ", -info, " had an illegal value"
        else
            print *, "LAPACK error: matrix is singular, solution could not be computed"
        end if

        print *, "End LAPACK Example"

    end subroutine s_lapack_example_solve_linear_system

    !> @brief Example subroutine demonstrating LAPACK usage for computing
    !!        eigenvalues of a symmetric matrix using DSYEV/SSYEV
    impure subroutine s_lapack_example_eigenvalues()

        ! Local variables for eigenvalue computation
        integer, parameter :: n = 3     ! Size of the matrix
        real(wp), dimension(n, n) :: A  ! Symmetric matrix
        real(wp), dimension(n) :: w     ! Eigenvalues
        real(wp), dimension(3*n) :: work ! Work array
        integer, parameter :: lwork = 3*n          ! Size of work array
        integer :: info                 ! Return status
        character, parameter :: jobz = 'N'         ! Compute eigenvalues only
        character, parameter :: uplo = 'U'         ! Upper triangular part of A

        ! Only run on the root process to avoid duplicate output
        if (proc_rank /= 0) return

        ! Set up a simple symmetric 3x3 matrix
        A(1, :) = [4.0_wp, 1.0_wp, 1.0_wp]
        A(2, :) = [1.0_wp, 4.0_wp, 1.0_wp]
        A(3, :) = [1.0_wp, 1.0_wp, 4.0_wp]

        print *, "LAPACK Eigenvalue Example"
        print *, "Computing eigenvalues of symmetric matrix:"
        print *, "A = [4 1 1; 1 4 1; 1 1 4]"

        ! Call appropriate LAPACK routine based on precision
#ifdef MFC_SINGLE_PRECISION
        call ssyev(jobz, uplo, n, A, n, w, work, lwork, info)
        print *, "Using LAPACK (SSYEV)"
#else
        call dsyev(jobz, uplo, n, A, n, w, work, lwork, info)
        print *, "Using LAPACK (DSYEV)"
#endif

        ! Check for success
        if (info == 0) then
            print *, "Eigenvalues computed successfully!"
            print *, "Eigenvalues: [", w(1), ", ", w(2), ", ", w(3), "]"
            print *, "Expected: [2, 5, 5] (approximately)"
        else if (info < 0) then
            print *, "LAPACK error: argument ", -info, " had an illegal value"
        else
            print *, "LAPACK error: algorithm failed to converge"
        end if

        print *, "End LAPACK Eigenvalue Example"

    end subroutine s_lapack_example_eigenvalues

end module m_lapack_example
