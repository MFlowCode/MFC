! test02_dt_array_constructor.f90
! Array constructors assigned to derived type array members inside a target region.
!   v%L = [a, b, c]
! Known AMD bug: this causes HSA_STATUS_ERROR_MEMORY_APERTURE_VIOLATION.
! Expected failure on AMD flang <= 6.x; may work on 7.x+.
program test02_dt_array_constructor
    implicit none
    integer, parameter :: wp = 8
    integer, parameter :: N = 10000

    type :: vec3
        real(wp) :: L(3), R(3)
    end type

    real(wp) :: res(N)
    type(vec3) :: v
    real(wp) :: a, b, c, expected
    integer :: i, nerr

    res = 0._wp

    !$omp target teams distribute parallel do map(from:res) private(v, a, b, c)
    do i = 1, N
        a = real(i, wp)
        b = real(i, wp) * 2._wp
        c = real(i, wp) * 3._wp
        v%L = [a, b, c]
        v%R = [-a, -b, -c]
        res(i) = v%L(1)**2._wp + v%L(2)**2._wp + v%L(3)**2._wp &
               + v%R(1)**2._wp + v%R(2)**2._wp + v%R(3)**2._wp
    end do
    !$omp end target teams distribute parallel do

    nerr = 0
    do i = 1, N
        expected = real(i, wp)**2._wp * (1._wp + 4._wp + 9._wp) * 2._wp
        if (abs(res(i) - expected) > 1.e-10_wp * expected) nerr = nerr + 1
    end do

    if (nerr == 0) then
        print *, "PASS test02: array constructor on derived type member"
    else
        print *, "FAIL test02:", nerr, "errors -- array constructor on derived type member"
    end if
end program test02_dt_array_constructor
