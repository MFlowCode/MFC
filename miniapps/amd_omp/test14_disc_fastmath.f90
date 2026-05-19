! test14_disc_fastmath.f90
! Tests whether -fopenmp-target-fast causes disc to go negative on AMD GPU.
! disc = (c^2 + B2/rho)^2 - 4*c^2*(Bn^2/rho) is guaranteed >= 0 by math.
program test14_disc_fastmath
    implicit none
    integer, parameter :: wp = 8
    integer, parameter :: N = 100000

    real(wp) :: disc_min, disc_i, B2loc, term
    real(wp) :: rho_in(N), c_in(N), B1_in(N), B2_in(N), B3_in(N)
    integer  :: i, n_neg

    ! Near-degenerate: Bn^2 ~ B^2 maximises cancellation in the discriminant
    do i = 1, N
        rho_in(i) = 0.1_wp + 0.9_wp * real(mod(i, 100), wp) / 100._wp
        c_in(i)   = 0.5_wp + 0.5_wp * real(mod(i, 37),  wp) / 37._wp
        B1_in(i)  = 1._wp  + real(mod(i, 53), wp) / 53._wp
        B2_in(i)  = 1.e-6_wp * real(mod(i, 7),  wp)
        B3_in(i)  = 1.e-6_wp * real(mod(i, 11), wp)
    end do

    disc_min = huge(1._wp)
    n_neg = 0

    !$omp target teams distribute parallel do &
    !$omp   map(to:rho_in,c_in,B1_in,B2_in,B3_in) &
    !$omp   reduction(min:disc_min) reduction(+:n_neg) &
    !$omp   private(B2loc,term,disc_i)
    do i = 1, N
        B2loc  = B1_in(i)**2 + B2_in(i)**2 + B3_in(i)**2
        term   = c_in(i)**2 + B2loc / rho_in(i)
        disc_i = term**2 - 4._wp * c_in(i)**2 * (B1_in(i)**2 / rho_in(i))
        disc_min = min(disc_min, disc_i)
        if (disc_i < 0._wp) n_neg = n_neg + 1
    end do
    !$omp end target teams distribute parallel do

    print *, "min(disc) =", disc_min
    print *, "n_neg     =", n_neg
    if (n_neg == 0) then
        print *, "PASS test14: disc >= 0 for all cells"
    else
        print *, "FAIL test14:", n_neg, "cells with disc < 0 -- fast-math is the culprit"
    end if
end program test14_disc_fastmath
