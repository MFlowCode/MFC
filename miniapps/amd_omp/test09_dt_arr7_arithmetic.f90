! test09_dt_arr7_arithmetic.f90
!
! Tests whether PRIVATE riemann_states_arr7-like DT variables (L(7), R(7))
! correctly compute F_star = F + s*(U_star - U) element-by-element inside
! an OpenMP target offload parallel loop.
!
! This is the exact HLLD F_star computation pattern.
program test09_dt_arr7_arithmetic
    implicit none
    integer, parameter :: wp = 8
    integer, parameter :: N = 10000

    type :: arr7
        real(wp) :: L(7), R(7)
    end type
    type :: scalar_lr
        real(wp) :: L, R
    end type

    real(wp) :: F_star_out(N, 7)
    real(wp) :: F_in(N, 7), U_in(N, 7), U_star_in(N, 7), s_in(N)
    type(arr7) :: F, U, U_star, F_star
    type(scalar_lr) :: s_lr
    real(wp) :: expected
    integer :: i, j, nerr

    do i = 1, N
        s_in(i) = -2.5_wp * real(i, wp)
        do j = 1, 7
            F_in(i, j) = real(i*7 + j, wp) * 0.1_wp
            U_in(i, j) = real(i*7 + j, wp) * 0.2_wp
            U_star_in(i, j) = real(i*7 + j, wp) * 0.15_wp
        end do
    end do
    F_star_out = 0._wp

    !$omp target teams distribute parallel do &
    !$omp   map(to:F_in,U_in,U_star_in,s_in) map(from:F_star_out) &
    !$omp   private(F, U, U_star, F_star, s_lr, j)
    do i = 1, N
        s_lr%L = s_in(i)
        do j = 1, 7
            F%L(j) = F_in(i, j)
            U%L(j) = U_in(i, j)
            U_star%L(j) = U_star_in(i, j)
        end do
        do j = 1, 7
            F_star%L(j) = F%L(j) + s_lr%L*(U_star%L(j) - U%L(j))
        end do
        do j = 1, 7
            F_star_out(i, j) = F_star%L(j)
        end do
    end do
    !$omp end target teams distribute parallel do

    nerr = 0
    do i = 1, N
        do j = 1, 7
            expected = F_in(i, j) + s_in(i)*(U_star_in(i, j) - U_in(i, j))
            if (abs(F_star_out(i, j) - expected) > 1.e-10_wp*max(abs(expected), 1._wp)) nerr = nerr + 1
        end do
    end do

    if (nerr == 0) then
        print *, "PASS test09: private arr7 DT F_star = F + s*(U_star - U)"
    else
        print *, "FAIL test09:", nerr, "errors -- arr7 DT element arithmetic"
    end if
end program test09_dt_arr7_arithmetic
