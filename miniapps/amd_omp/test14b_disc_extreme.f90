! test14b_disc_extreme.f90
! More extreme disc test: B nearly perfectly aligned with normal (Bn^2 ~ B^2),
! plus many sign/magnitude combinations to maximize FP cancellation.
program test14b_disc_extreme
    implicit none
    integer, parameter :: wp = 8
    integer, parameter :: N = 200000

    real(wp) :: disc_min, disc_i, B2loc, term
    real(wp) :: rho, c_sq, Bn, Bt1, Bt2
    integer  :: i, n_neg, n_tiny

    disc_min = huge(1._wp)
    n_neg = 0; n_tiny = 0

    !$omp target teams distribute parallel do &
    !$omp   reduction(min:disc_min) reduction(+:n_neg) reduction(+:n_tiny) &
    !$omp   private(rho,c_sq,Bn,Bt1,Bt2,B2loc,term,disc_i)
    do i = 1, N
        ! Span a wide range to find pathological cases
        rho  = 1.e-3_wp * (1._wp + real(mod(i, 997), wp))
        c_sq = 1.e-2_wp * (1._wp + real(mod(i, 503), wp))
        ! Bn very close to total B (nearly-parallel field)
        Bn   = 1._wp + 0.001_wp * real(mod(i, 101), wp)
        Bt1  = 1.e-4_wp * real(mod(i, 17), wp)
        Bt2  = 1.e-4_wp * real(mod(i, 13), wp)

        B2loc  = Bn**2 + Bt1**2 + Bt2**2
        term   = c_sq + B2loc / rho
        disc_i = term**2 - 4._wp * c_sq * (Bn**2 / rho)

        disc_min = min(disc_min, disc_i)
        if (disc_i < 0._wp)        n_neg  = n_neg  + 1
        if (disc_i < 1.e-14_wp)    n_tiny = n_tiny + 1
    end do
    !$omp end target teams distribute parallel do

    print *, "min(disc) =", disc_min
    print *, "n_neg (disc<0)      =", n_neg
    print *, "n_tiny (disc<1e-14) =", n_tiny
    if (n_neg == 0) then
        print *, "PASS test14b: disc >= 0 always"
    else
        print *, "FAIL test14b:", n_neg, "cells with disc < 0"
    end if
end program test14b_disc_extreme
