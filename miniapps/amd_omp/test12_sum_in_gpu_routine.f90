! test12_sum_in_gpu_routine.f90
! Tests whether sum(B**2) on a plain real(wp) array passed as
! an argument to a GPU device subroutine works correctly.
! Also tests the explicit-loop alternative.
program test12_sum_in_gpu_routine
    implicit none
    integer, parameter :: wp = 8
    integer, parameter :: N = 1000

    ! Use separate arrays for each B component to avoid slicing issues
    real(wp) :: c_fast_sum(N), c_fast_loop(N), c_fast_ref(N)
    real(wp) :: rho_in(N), c_in(N)
    real(wp) :: B1_in(N), B2_in(N), B3_in(N)
    integer :: i, nerr_sum, nerr_loop

    do i = 1, N
        rho_in(i) = 1.0_wp + 0.001_wp*real(i,wp)
        c_in(i) = 0.9_wp
        B1_in(i) = 0.5641895835_wp
        B2_in(i) = 1.0149412604_wp * (1._wp + 0.0001_wp*real(i,wp))
        B3_in(i) = 0.5641895835_wp
    end do

    c_fast_sum = 0._wp
    c_fast_loop = 0._wp

    ! Test 1: sum() version
    !$omp target teams distribute parallel do &
    !$omp   map(to:rho_in,c_in,B1_in,B2_in,B3_in) map(from:c_fast_sum)
    do i = 1, N
        call compute_c_fast_sum(rho_in(i), c_in(i), B1_in(i), B2_in(i), B3_in(i), c_fast_sum(i))
    end do
    !$omp end target teams distribute parallel do

    ! Test 2: explicit-loop version
    !$omp target teams distribute parallel do &
    !$omp   map(to:rho_in,c_in,B1_in,B2_in,B3_in) map(from:c_fast_loop)
    do i = 1, N
        call compute_c_fast_loop(rho_in(i), c_in(i), B1_in(i), B2_in(i), B3_in(i), c_fast_loop(i))
    end do
    !$omp end target teams distribute parallel do

    ! Reference on CPU
    do i = 1, N
        call compute_c_fast_loop(rho_in(i), c_in(i), B1_in(i), B2_in(i), B3_in(i), c_fast_ref(i))
    end do

    nerr_sum = 0; nerr_loop = 0
    do i = 1, N
        if (abs(c_fast_sum(i) - c_fast_ref(i)) > 1.e-10_wp * c_fast_ref(i)) nerr_sum = nerr_sum + 1
        if (abs(c_fast_loop(i) - c_fast_ref(i)) > 1.e-10_wp * c_fast_ref(i)) nerr_loop = nerr_loop + 1
    end do

    if (nerr_sum == 0) then
        print *, "PASS test12a: sum(B**2) in GPU device routine"
    else
        print *, "FAIL test12a:", nerr_sum, "errors -- sum(B**2) in GPU device routine"
        print *, "  GPU c_fast(1)=", c_fast_sum(1), " ref=", c_fast_ref(1)
    end if

    if (nerr_loop == 0) then
        print *, "PASS test12b: explicit loop B2=B1^2+B2^2+B3^2 in GPU device routine"
    else
        print *, "FAIL test12b:", nerr_loop, "errors -- explicit loop in GPU device routine"
        print *, "  GPU c_fast(1)=", c_fast_loop(1), " ref=", c_fast_ref(1)
    end if

contains

    subroutine compute_c_fast_sum(rho, c, Bx, By, Bz, c_fast)
        !$omp declare target
        real(wp), intent(in)  :: rho, c, Bx, By, Bz
        real(wp), intent(out) :: c_fast
        real(wp)              :: B(3), B2, term, disc
        B(1) = Bx; B(2) = By; B(3) = Bz
        B2 = sum(B**2)
        term = c**2 + B2/rho
        disc = term**2 - 4._wp*c**2*(B(1)**2/rho)
        disc = max(disc, 0._wp)
        c_fast = sqrt(0.5_wp*(term + sqrt(disc)))
    end subroutine

    subroutine compute_c_fast_loop(rho, c, Bx, By, Bz, c_fast)
        !$omp declare target
        real(wp), intent(in)  :: rho, c, Bx, By, Bz
        real(wp), intent(out) :: c_fast
        real(wp)              :: B2, term, disc
        B2 = Bx**2 + By**2 + Bz**2
        term = c**2 + B2/rho
        disc = term**2 - 4._wp*c**2*(Bx**2/rho)
        disc = max(disc, 0._wp)
        c_fast = sqrt(0.5_wp*(term + sqrt(disc)))
    end subroutine

end program test12_sum_in_gpu_routine
