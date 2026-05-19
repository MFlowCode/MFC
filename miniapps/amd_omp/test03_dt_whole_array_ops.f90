! test03_dt_whole_array_ops.f90
! Whole-array intrinsics on derived type array members in target region:
!   sum(v%L**2)    -- magnitude squared
!   v%L**2         -- element-wise squaring (whole-array expression)
! Known AMD bug: produces wrong results or zero.
program test03_dt_whole_array_ops
    implicit none
    integer, parameter :: wp = 8
    integer, parameter :: N = 10000

    type :: vec3
        real(wp) :: L(3), R(3)
    end type

    real(wp) :: res_sum(N), res_max(N)
    type(vec3) :: v
    real(wp) :: expected_sum, expected_max
    integer :: i, nerr_sum, nerr_max

    res_sum = 0._wp
    res_max = 0._wp

    !$omp target teams distribute parallel do map(from:res_sum,res_max) private(v)
    do i = 1, N
        v%L(1) = real(i, wp)
        v%L(2) = real(i, wp) * 2._wp
        v%L(3) = real(i, wp) * 3._wp
        res_sum(i) = sum(v%L**2._wp)
        res_max(i) = maxval(v%L)
    end do
    !$omp end target teams distribute parallel do

    nerr_sum = 0
    nerr_max = 0
    do i = 1, N
        expected_sum = real(i, wp)**2._wp * (1._wp + 4._wp + 9._wp)
        expected_max = real(i, wp) * 3._wp
        if (abs(res_sum(i) - expected_sum) > 1.e-10_wp * expected_sum) nerr_sum = nerr_sum + 1
        if (abs(res_max(i) - expected_max) > 1.e-10_wp * expected_max) nerr_max = nerr_max + 1
    end do

    if (nerr_sum == 0) then
        print *, "PASS test03a: sum(v%L**2) on private derived type member"
    else
        print *, "FAIL test03a:", nerr_sum, "errors -- sum(v%L**2) on private derived type member"
    end if

    if (nerr_max == 0) then
        print *, "PASS test03b: maxval(v%L) on private derived type member"
    else
        print *, "FAIL test03b:", nerr_max, "errors -- maxval(v%L) on private derived type member"
    end if
end program test03_dt_whole_array_ops
