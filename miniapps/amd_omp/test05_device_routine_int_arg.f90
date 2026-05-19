! test05_device_routine_int_arg.f90
! Tests whether passing an integer literal constant vs an integer variable
! as an argument to a device (declare target) subroutine gives the same result.
!
! The specific pattern from MFC's HLLD solver:
!   call sub(B, 1, ...)      -- literal index into B(3)
!   call sub(B, norm, ...)   -- variable, equals 1 at runtime
!   inline: B(1)**2          -- direct element access, no call
!
! Known AMD bug (flang <= 6.x): literal integer argument causes wrong GPU code.
program test05_device_routine_int_arg
    implicit none
    integer, parameter :: wp = 8
    integer, parameter :: N = 10000

    real(wp) :: res_lit(N), res_var(N), res_inline(N)
    real(wp) :: B1(N), B2(N), rho(N), c_sq(N)
    real(wp) :: B2_total, term, disc
    integer :: i, norm_idx, nerr_lit, nerr_var

    do i = 1, N
        rho(i)  = 1._wp + real(i, wp) * 0.001_wp
        c_sq(i) = 0.5_wp + real(i, wp) * 0.0001_wp
        B1(i)   = 0.1_wp * real(i, wp) * 0.001_wp
        B2(i)   = 0.05_wp + real(i, wp) * 0.00001_wp
    end do

    !$omp target enter data map(to: rho, c_sq, B1, B2)
    !$omp target enter data map(alloc: res_lit, res_var, res_inline)

    ! Pass integer literal 1
    !$omp target teams distribute parallel do
    do i = 1, N
        call magnetosonic(rho(i), c_sq(i), B1(i), B2(i), 1, res_lit(i))
    end do
    !$omp end target teams distribute parallel do

    ! Pass integer variable = 1
    !$omp target teams distribute parallel do private(norm_idx)
    do i = 1, N
        norm_idx = 1
        call magnetosonic(rho(i), c_sq(i), B1(i), B2(i), norm_idx, res_var(i))
    end do
    !$omp end target teams distribute parallel do

    ! Inline: explicit B(1) reference, no function call
    !$omp target teams distribute parallel do private(B2_total, term, disc)
    do i = 1, N
        B2_total = B1(i)**2._wp + B2(i)**2._wp
        term = c_sq(i) + B2_total / rho(i)
        disc = term**2._wp - 4._wp * c_sq(i) * (B1(i)**2._wp / rho(i))
        res_inline(i) = sqrt(max(0._wp, 0.5_wp * (term + sqrt(max(0._wp, disc)))))
    end do
    !$omp end target teams distribute parallel do

    !$omp target exit data map(from: res_lit, res_var, res_inline)
    !$omp target exit data map(delete: rho, c_sq, B1, B2)

    nerr_lit = 0
    nerr_var = 0
    do i = 1, N
        if (abs(res_lit(i) - res_inline(i)) > 1.e-12_wp * abs(res_inline(i)) + 1.e-15_wp) &
            nerr_lit = nerr_lit + 1
        if (abs(res_var(i) - res_inline(i)) > 1.e-12_wp * abs(res_inline(i)) + 1.e-15_wp) &
            nerr_var = nerr_var + 1
    end do

    if (nerr_lit == 0) then
        print *, "PASS test05a: integer literal 1 as device routine arg matches inline"
    else
        print *, "FAIL test05a:", nerr_lit, "errors -- literal int arg differs from inline"
        print *, "  sample: res_lit(1)=", res_lit(1), " vs inline=", res_inline(1)
    end if

    if (nerr_var == 0) then
        print *, "PASS test05b: integer variable (=1) as device routine arg matches inline"
    else
        print *, "FAIL test05b:", nerr_var, "errors -- variable int arg differs from inline"
    end if

contains

    subroutine magnetosonic(rho, c_sq, B_n, B_t, norm_idx, c_fast)
        !$omp declare target
        real(wp), intent(in)  :: rho, c_sq, B_n, B_t
        integer,  intent(in)  :: norm_idx
        real(wp), intent(out) :: c_fast
        real(wp) :: B_arr(2), B2_total, term, disc

        B_arr(1) = B_n
        B_arr(2) = B_t
        B2_total = B_arr(1)**2._wp + B_arr(2)**2._wp
        term = c_sq + B2_total / rho
        ! norm_idx selects the normal B component (1 = B_n, 2 = B_t)
        disc = term**2._wp - 4._wp * c_sq * (B_arr(norm_idx)**2._wp / rho)
        c_fast = sqrt(max(0._wp, 0.5_wp * (term + sqrt(max(0._wp, disc)))))
    end subroutine

end program test05_device_routine_int_arg
