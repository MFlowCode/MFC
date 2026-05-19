! test08_dt_private_array_elems.f90
!
! Tests whether a PRIVATE derived type with array members correctly
! maintains all element values when assigned element-by-element inside
! an OpenMP target offload parallel loop.
!
! Simulates the HLLD B%L assignment pattern:
!   B%L(1) = val1; B%L(2) = val2; B%L(3) = val3
!   pres_mag = 0.5*(B%L(1)**2 + B%L(2)**2 + B%L(3)**2)
!
! Key check: if AMD has a bug where B%L(2) or B%L(3) are stale/zero after
! element-by-element assignment in a private DT, pres_mag will be wrong.
program test08_dt_private_array_elems
    implicit none
    integer, parameter :: wp = 8
    integer, parameter :: N = 10000

    type :: vec3
        real(wp) :: L(3), R(3)
    end type

    real(wp) :: pres_mag_out(N)
    real(wp) :: Bx_in(N), By_in(N), Bz_in(N)
    type(vec3) :: B
    real(wp) :: expected
    integer :: i, nerr

    do i = 1, N
        Bx_in(i) = 0.5641895835_wp * real(i, wp)
        By_in(i) = 1.0149412604_wp * real(i, wp)
        Bz_in(i) = 0.5641895835_wp * real(i, wp) * 0.7_wp
    end do
    pres_mag_out = 0._wp

    !$omp target teams distribute parallel do &
    !$omp   map(to:Bx_in,By_in,Bz_in) map(from:pres_mag_out) &
    !$omp   private(B)
    do i = 1, N
        B%L(1) = Bx_in(i)
        B%L(2) = By_in(i)
        B%L(3) = Bz_in(i)
        pres_mag_out(i) = 0.5_wp*(B%L(1)**2 + B%L(2)**2 + B%L(3)**2)
    end do
    !$omp end target teams distribute parallel do

    nerr = 0
    do i = 1, N
        expected = 0.5_wp*(Bx_in(i)**2 + By_in(i)**2 + Bz_in(i)**2)
        if (abs(pres_mag_out(i) - expected) > 1.e-10_wp*max(abs(expected), 1._wp)) nerr = nerr + 1
    end do

    if (nerr == 0) then
        print *, "PASS test08: private DT array member elem-by-elem assign and read"
    else
        print *, "FAIL test08:", nerr, "errors -- private DT array member elem access"
    end if
end program test08_dt_private_array_elems
