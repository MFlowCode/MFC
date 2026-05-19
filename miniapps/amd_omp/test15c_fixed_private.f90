! test15c_fixed_private.f90
! Same as test15b but uses fixed-size dimension(3) instead of dimension(nf).
! This is the USING_AMD workaround in MFC.
program test15c_fixed_private
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
        print *, "PASS test15c: fixed dim(3) private works correctly"
    else
        print *, "FAIL test15c:", nerr, "cells wrong"
        print *, "  cell 1: expected", inp(1), "got", out(1)
    end if

contains

    subroutine run_test(nf, N, inp, out)
        integer, intent(in)  :: nf, N
        real(wp), intent(in) :: inp(N)
        real(wp), intent(out):: out(N)
        real(wp) :: local(3)   ! fixed size -- same as USING_AMD workaround
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

end program test15c_fixed_private
