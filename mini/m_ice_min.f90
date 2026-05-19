! ABSOLUTE MINIMUM reproducer for ifx 2025.1.1 SPIR-V ICE #5633
!
! Trigger: matmul() inside a !$omp declare target subroutine that is
!          inlined into a !$omp target teams loop kernel.
!
! ifx version: 2025.1.1 20250418   Hardware: Intel GPU Max 1100 (Ponte Vecchio)
!
! ICE (O1/O2/O3):
!   ifx -free -fiopenmp -fopenmp-targets=spir64 -O3 m_ice_min.f90
! OK (inlining disabled):
!   ifx -free -fiopenmp -fopenmp-targets=spir64 -O3 -fno-inline m_ice_min.f90

module m_min

    implicit none
    private

    integer, parameter :: wp = kind(1.0d0)

    real(wp), dimension(3, 3) :: A  ! module-level 3x3 matrix
    real(wp), dimension(3)    :: x  ! module-level vector
    real(wp), dimension(3)    :: b  ! result vector
    !$omp declare target(A, x, b)

    public :: s_run

contains

    subroutine s_run(n)
        integer, intent(in) :: n
        integer :: i

        A(1,:) = [1._wp, 0._wp, 0._wp]
        A(2,:) = [0._wp, 1._wp, 0._wp]
        A(3,:) = [0._wp, 0._wp, 1._wp]
        x = [1._wp, 2._wp, 3._wp]
        b = 0._wp
        !$omp target enter data map(to: A, x, b)

        !$omp target teams loop private(i) map(to:A,x) map(tofrom:b)
        do i = 1, n
            call s_apply(i)
        end do

        !$omp target exit data map(from: b)
    end subroutine s_run

    subroutine s_apply(k)
        !$omp declare target
        integer, intent(in) :: k
        ! matmul inside declare-target sub -- ICE trigger when inlined
        b = real(k, wp) * matmul(A, x)
    end subroutine s_apply

end module m_min

program test_min
    use m_min
    implicit none
    call s_run(4)
end program test_min
