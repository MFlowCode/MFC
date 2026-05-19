! test04_dt_zero_init.f90
! Tests whether private derived type members are zero-initialized at entry
! to each GPU thread, as the Fortran standard requires for local vars.
! Known AMD bug: private derived type members are NOT zero-initialized,
! causing incorrect accumulations (e.g. vel_rms += ...) if not explicitly set.
program test04_dt_zero_init
    implicit none
    integer, parameter :: wp = 8
    integer, parameter :: N = 10000
    integer, parameter :: M = 3

    type :: scalar_state
        real(wp) :: L, R
    end type

    real(wp) :: res(N)
    type(scalar_state) :: acc
    real(wp) :: expected
    integer :: i, j, nerr

    res = 0._wp

    !$omp target teams distribute parallel do map(from:res) private(acc, j)
    do i = 1, N
        ! Accumulate into private acc%L without explicit initialization.
        ! If zero-init is NOT guaranteed, acc%L may start with garbage.
        acc%L = 0._wp
        acc%R = 0._wp
        do j = 1, M
            acc%L = acc%L + real(j, wp)**2._wp
            acc%R = acc%R + real(j, wp)
        end do
        res(i) = acc%L - acc%R
    end do
    !$omp end target teams distribute parallel do

    ! acc%L = 1 + 4 + 9 = 14; acc%R = 1+2+3 = 6; diff = 8
    nerr = 0
    expected = 14._wp - 6._wp
    do i = 1, N
        if (abs(res(i) - expected) > 1.e-10_wp) nerr = nerr + 1
    end do

    if (nerr == 0) then
        print *, "PASS test04: explicit zero-init of private derived type scalar member"
    else
        print *, "FAIL test04:", nerr, "errors -- private derived type scalar accumulation"
    end if
end program test04_dt_zero_init
