! test15_vla_private.f90
! Tests whether AMD OpenMP target correctly handles private variable-length arrays.
! In non-case-optimized MFC, num_fluids is a runtime integer, making
! dimension(num_fluids) a VLA. If AMD doesn't privatize it, threads corrupt each other.
program test15_vla_private
    implicit none
    integer, parameter :: wp = 8
    integer, parameter :: NCELLS = 2500  ! 50x50, same as 2D HLLD test

    integer :: num_fluids, i, j, nerr
    real(wp) :: alpha_rho_in(NCELLS), result_sum(NCELLS)
    real(wp) :: alpha_rho_local(3)  ! worst case: declare as 3 to simulate the VLA issue

    num_fluids = 1

    do i = 1, NCELLS
        alpha_rho_in(i) = real(i, wp) * 0.001_wp
    end do

    result_sum = 0._wp

    ! Test with runtime-variable loop bound to simulate VLA behavior
    !$omp target teams distribute parallel do &
    !$omp   map(to:alpha_rho_in,num_fluids) map(from:result_sum) &
    !$omp   private(alpha_rho_local,j)
    do i = 1, NCELLS
        do j = 1, num_fluids
            alpha_rho_local(j) = alpha_rho_in(i)
        end do
        result_sum(i) = 0._wp
        do j = 1, num_fluids
            result_sum(i) = result_sum(i) + alpha_rho_local(j)
        end do
    end do
    !$omp end target teams distribute parallel do

    nerr = 0
    do i = 1, NCELLS
        if (abs(result_sum(i) - alpha_rho_in(i)) > 1.e-12_wp) nerr = nerr + 1
    end do

    if (nerr == 0) then
        print *, "PASS test15: fixed-dim(3) private works correctly"
    else
        print *, "FAIL test15:", nerr, "cells wrong"
        print *, "  cell 1: expected", alpha_rho_in(1), "got", result_sum(1)
    end if
end program test15_vla_private
