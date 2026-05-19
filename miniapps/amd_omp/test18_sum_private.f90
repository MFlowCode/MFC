! test18_sum_private.f90
! Tests sum() on private fixed-size arrays inside a target loop.
! Checks:
!   (a) sum(v**2) where v is a private fixed array dimension(3)
!   (b) sum(v**2) where v is a private fixed array dimension(7)
!   (c) explicit accumulation loop (workaround)
! These mirror vel_rms = sum(vel**2) and pres_mag = sum(B**2) in HLLD.
program test18_sum_private
    implicit none
    integer, parameter :: wp = 8
    integer, parameter :: N = 2000

    real(wp) :: vx(N), vy(N), vz(N)
    real(wp) :: sum3_result(N), sum7_result(N), loop3_result(N)
    real(wp) :: ref(N)
    real(wp) :: tmp3(3), tmp7(7)
    integer  :: i, k, nerr

    do i = 1, N
        vx(i) = 0.1_wp * real(i, wp)
        vy(i) = 0.2_wp * real(i, wp)
        vz(i) = 0.3_wp * real(i, wp)
        ref(i) = vx(i)**2 + vy(i)**2 + vz(i)**2
    end do

    ! Test A: sum() on private dimension(3) array
    !$omp target teams distribute parallel do &
    !$omp   map(to:vx,vy,vz) map(from:sum3_result) private(tmp3)
    do i = 1, N
        tmp3(1) = vx(i); tmp3(2) = vy(i); tmp3(3) = vz(i)
        sum3_result(i) = sum(tmp3**2)
    end do
    !$omp end target teams distribute parallel do

    ! Test B: sum() on private dimension(7) array (like U in HLLD)
    !$omp target teams distribute parallel do &
    !$omp   map(to:vx,vy,vz) map(from:sum7_result) private(tmp7)
    do i = 1, N
        tmp7(1) = vx(i); tmp7(2) = vy(i); tmp7(3) = vz(i)
        tmp7(4) = vx(i)*vy(i); tmp7(5) = vy(i)*vz(i)
        tmp7(6) = vz(i)*vx(i); tmp7(7) = 0._wp
        sum7_result(i) = sum(tmp7(1:3)**2)
    end do
    !$omp end target teams distribute parallel do

    ! Test C: explicit accumulation (workaround)
    !$omp target teams distribute parallel do &
    !$omp   map(to:vx,vy,vz) map(from:loop3_result)
    do i = 1, N
        loop3_result(i) = vx(i)**2 + vy(i)**2 + vz(i)**2
    end do
    !$omp end target teams distribute parallel do

    nerr = 0
    do i = 1, N
        if (abs(sum3_result(i) - ref(i)) > 1.e-8_wp*ref(i)) nerr = nerr + 1
    end do
    if (nerr == 0) then
        print *, "PASS test18a: sum(v**2) on private dim(3) array"
    else
        print *, "FAIL test18a:", nerr, "wrong -- sum(v**2) on private dim(3) broken"
        print *, "  cell 1: got", sum3_result(1), "ref", ref(1)
    end if

    nerr = 0
    do i = 1, N
        if (abs(sum7_result(i) - ref(i)) > 1.e-8_wp*ref(i)) nerr = nerr + 1
    end do
    if (nerr == 0) then
        print *, "PASS test18b: sum(v**2) on private dim(7) slice"
    else
        print *, "FAIL test18b:", nerr, "wrong -- sum(v**2) on dim(7) slice broken"
        print *, "  cell 1: got", sum7_result(1), "ref", ref(1)
    end if

    nerr = 0
    do i = 1, N
        if (abs(loop3_result(i) - ref(i)) > 1.e-8_wp*ref(i)) nerr = nerr + 1
    end do
    if (nerr == 0) then
        print *, "PASS test18c: explicit accumulation workaround"
    else
        print *, "FAIL test18c:", nerr, "wrong -- explicit accumulation broken"
    end if

end program test18_sum_private
