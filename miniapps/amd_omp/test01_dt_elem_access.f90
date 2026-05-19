! test01_dt_elem_access.f90
! Derived type with real array members as private in target parallel loop.
! Basic element-by-element access: v%L(1), v%L(2), v%L(3).
! Expected: always works (this is the baseline / reference pattern).
program test01_dt_elem_access
    implicit none
    integer, parameter :: wp = 8
    integer, parameter :: N = 10000

    type :: vec3
        real(wp) :: L(3), R(3)
    end type

    real(wp) :: res(N)
    type(vec3) :: v
    real(wp) :: expected
    integer :: i, nerr

    res = 0._wp

    !$omp target teams distribute parallel do map(from:res) private(v)
    do i = 1, N
        v%L(1) = real(i, wp)
        v%L(2) = real(i, wp) * 2._wp
        v%L(3) = real(i, wp) * 3._wp
        v%R(1) = -real(i, wp)
        v%R(2) = -real(i, wp) * 2._wp
        v%R(3) = -real(i, wp) * 3._wp
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
        print *, "PASS test01: derived type element access in private target var"
    else
        print *, "FAIL test01:", nerr, "errors -- derived type element access"
    end if
end program test01_dt_elem_access
