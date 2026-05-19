! test16_array_constructor.f90
! Tests whether array constructors [a, b, c] work inside an OpenMP target loop.
! In MFC's HLLD solver we replaced:
!   B%L = [Bx0, qL(...), qL(...)]
!   U%L = [rho, rho*vel(1:3), B(2:3), E]
! with explicit element assignments because AMD flang may fail on constructors.
program test16_array_constructor
    implicit none
    integer, parameter :: wp = 8
    integer, parameter :: N = 2000

    real(wp) :: a_in(N), b_in(N), c_in(N)
    real(wp) :: v_constructor(N, 3), v_explicit(N, 3)
    integer  :: i, nerr

    do i = 1, N
        a_in(i) = real(i, wp) * 0.001_wp
        b_in(i) = real(i, wp) * 0.002_wp
        c_in(i) = real(i, wp) * 0.003_wp
    end do
    v_constructor = 0._wp; v_explicit = 0._wp

    ! Test A: array constructor inside parallel loop
    !$omp target teams distribute parallel do &
    !$omp   map(to:a_in,b_in,c_in) map(from:v_constructor)
    do i = 1, N
        v_constructor(i, :) = [a_in(i), b_in(i), c_in(i)]
    end do
    !$omp end target teams distribute parallel do

    ! Test B: explicit element assignment (the workaround)
    !$omp target teams distribute parallel do &
    !$omp   map(to:a_in,b_in,c_in) map(from:v_explicit)
    do i = 1, N
        v_explicit(i, 1) = a_in(i)
        v_explicit(i, 2) = b_in(i)
        v_explicit(i, 3) = c_in(i)
    end do
    !$omp end target teams distribute parallel do

    nerr = 0
    do i = 1, N
        if (any(abs(v_constructor(i,:) - [a_in(i),b_in(i),c_in(i)]) > 1.e-14_wp)) nerr = nerr + 1
    end do
    if (nerr == 0) then
        print *, "PASS test16a: array constructor [a,b,c] in target loop"
    else
        print *, "FAIL test16a:", nerr, "cells wrong -- array constructor broken"
        print *, "  cell 1: got", v_constructor(1,:), "expected", a_in(1), b_in(1), c_in(1)
    end if

    nerr = 0
    do i = 1, N
        if (any(abs(v_explicit(i,:) - [a_in(i),b_in(i),c_in(i)]) > 1.e-14_wp)) nerr = nerr + 1
    end do
    if (nerr == 0) then
        print *, "PASS test16b: explicit element assignment in target loop"
    else
        print *, "FAIL test16b:", nerr, "cells wrong -- explicit assignment broken"
    end if

end program test16_array_constructor
