! test13_sum_explicit_shape.f90
! Tests sum(B**2) on an explicit-shape dummy array B(3) passed as argument
! to a GPU device routine -- exactly mimicking s_compute_fast_magnetosonic_speed.
program test13_sum_explicit_shape
    implicit none
    integer, parameter :: wp = 8
    integer, parameter :: N = 1000

    real(wp) :: c_fast_sum(N), c_fast_ref(N)
    real(wp) :: rho_in(N), c_in(N)
    ! Store B as SOA (struct of arrays) to avoid slicing issues
    real(wp) :: Bx_in(N), By_in(N), Bz_in(N)
    real(wp) :: B_vec(3)
    integer :: i, nerr

    do i = 1, N
        rho_in(i) = 1.0_wp + 0.001_wp*real(i,wp)
        c_in(i) = 0.9_wp
        Bx_in(i) = 0.5641895835_wp
        By_in(i) = 1.0149412604_wp * (1._wp + 0.0001_wp*real(i,wp))
        Bz_in(i) = 0.5641895835_wp
    end do

    c_fast_sum = 0._wp

    !$omp target teams distribute parallel do &
    !$omp   map(to:rho_in,c_in,Bx_in,By_in,Bz_in) map(from:c_fast_sum) private(B_vec)
    do i = 1, N
        B_vec(1) = Bx_in(i)
        B_vec(2) = By_in(i)
        B_vec(3) = Bz_in(i)
        call compute_c_fast(rho_in(i), c_in(i), B_vec, 1, c_fast_sum(i))
    end do
    !$omp end target teams distribute parallel do

    ! Reference on CPU (same code path)
    do i = 1, N
        B_vec(1) = Bx_in(i); B_vec(2) = By_in(i); B_vec(3) = Bz_in(i)
        call compute_c_fast(rho_in(i), c_in(i), B_vec, 1, c_fast_ref(i))
    end do

    nerr = 0
    do i = 1, N
        if (abs(c_fast_sum(i) - c_fast_ref(i)) > 1.e-10_wp * c_fast_ref(i)) nerr = nerr + 1
    end do

    if (nerr == 0) then
        print *, "PASS test13: sum(B**2) with explicit-shape B(3) arg in GPU routine"
    else
        print *, "FAIL test13:", nerr, "errors -- sum(B**2) with explicit-shape B(3) arg"
        print *, "  GPU c_fast(1)=", c_fast_sum(1), " ref=", c_fast_ref(1)
    end if

contains

    ! Mirrors s_compute_fast_magnetosonic_speed with the explicit-element fix
    subroutine compute_c_fast(rho, c, B, norm, c_fast)
        !$omp declare target
        real(wp), intent(in)  :: B(3), rho, c
        integer, intent(in)   :: norm
        real(wp), intent(out) :: c_fast
        real(wp)              :: B2, term, disc

        B2 = B(1)**2 + B(2)**2 + B(3)**2
        term = c**2 + B2/rho
        disc = term**2 - 4._wp*c**2*(B(norm)**2/rho)
        disc = max(disc, 0._wp)
        c_fast = sqrt(0.5_wp*(term + sqrt(disc)))
    end subroutine

end program test13_sum_explicit_shape
