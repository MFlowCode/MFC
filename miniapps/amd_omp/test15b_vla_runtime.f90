! test15b_vla_runtime.f90
! Tests whether OpenMP target private clause works for a VLA declared at
! subroutine scope with runtime size. In MFC non-case-optimized builds,
! alpha_rho_L(num_fluids) is such a VLA listed in the GPU_PARALLEL_LOOP private clause.
program test15b_vla_runtime
    implicit none
    integer, parameter :: wp = 8
    integer, parameter :: NCELLS = 2500

    integer :: nf, i, nerr
    real(wp) :: inp(NCELLS), out(NCELLS)

    nf = 1
    do i = 1, NCELLS
        inp(i) = real(i, wp) * 0.001_wp
    end do
    out = 0._wp

    call run_test(nf, NCELLS, inp, out)

    nerr = 0
    do i = 1, NCELLS
        if (abs(out(i) - inp(i)) > 1.e-12_wp) nerr = nerr + 1
    end do

    if (nerr == 0) then
        print *, "PASS test15b: runtime-VLA private works correctly"
    else
        print *, "FAIL test15b:", nerr, "cells wrong -- VLA private broken"
        print *, "  cell 1: expected", inp(1), "got", out(1)
        print *, "  cell 2: expected", inp(2), "got", out(2)
        print *, "  cell 3: expected", inp(3), "got", out(3)
    end if

contains

    subroutine run_test(nf, N, inp, out)
        integer, intent(in)  :: nf, N
        real(wp), intent(in) :: inp(N)
        real(wp), intent(out):: out(N)
        ! VLA at subroutine scope: size is runtime variable nf
        real(wp) :: local(nf)
        integer  :: i, j

        !$omp target teams distribute parallel do &
        !$omp   map(to:inp,nf) map(from:out) private(local,j)
        do i = 1, N
            do j = 1, nf
                local(j) = inp(i)
            end do
            out(i) = 0._wp
            do j = 1, nf
                out(i) = out(i) + local(j)
            end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine

end program test15b_vla_runtime
