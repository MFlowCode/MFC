! test06_dt_member_as_array_arg.f90
! Tests passing a derived type array member as an actual argument to a device
! subroutine that expects a plain array of the same size.
!   type :: vec3; real(wp) :: L(3); end type
!   call device_sub(v%L, ...)   -- passes L(3) as array argument
!
! The question: does AMD correctly pass the derived type member by reference,
! or does it corrupt the data due to misaligned/incorrect memory access?
program test06_dt_member_as_array_arg
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
        call sum_of_squares(v%L, res(i))
    end do
    !$omp end target teams distribute parallel do

    nerr = 0
    do i = 1, N
        expected = real(i, wp)**2._wp * (1._wp + 4._wp + 9._wp)
        if (abs(res(i) - expected) > 1.e-10_wp * expected) nerr = nerr + 1
    end do

    if (nerr == 0) then
        print *, "PASS test06: derived type member array passed as device routine arg"
    else
        print *, "FAIL test06:", nerr, "errors -- derived type member array as device routine arg"
    end if

contains

    subroutine sum_of_squares(B, result)
        !$omp declare target
        real(wp), intent(in)  :: B(3)
        real(wp), intent(out) :: result
        result = B(1)**2._wp + B(2)**2._wp + B(3)**2._wp
    end subroutine

end program test06_dt_member_as_array_arg
