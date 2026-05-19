! test17_whole_array_ops.f90
! Tests whether whole-array arithmetic in a target loop body works on AMD GPU.
! In MFC's HLLD solver we replaced patterns like:
!   F_star%L = F%L + s%L*(U_star%L - U%L)   (whole-array)
!   F_hlld = F%L                              (whole-array copy)
! with explicit do loops.
program test17_whole_array_ops
    implicit none
    integer, parameter :: wp = 8
    integer, parameter :: N = 2000
    integer, parameter :: M = 7

    real(wp) :: A_in(N, M), B_in(N, M), s_in(N)
    real(wp) :: out_whole(N, M), out_loop(N, M), ref(N, M)
    integer  :: i, k, nerr

    do i = 1, N
        s_in(i) = 0.5_wp + 0.001_wp*real(i, wp)
        do k = 1, M
            A_in(i, k) = real(i + k, wp) * 0.01_wp
            B_in(i, k) = real(i * k, wp) * 0.001_wp
        end do
    end do

    ! Test A: whole-array op inside private arrays in target loop
    !$omp target teams distribute parallel do &
    !$omp   map(to:A_in,B_in,s_in) map(from:out_whole) &
    !$omp   private(i,k)
    do i = 1, N
        out_whole(i, :) = A_in(i, :) + s_in(i)*(B_in(i, :) - A_in(i, :))
    end do
    !$omp end target teams distribute parallel do

    ! Test B: explicit loop (our workaround)
    !$omp target teams distribute parallel do &
    !$omp   map(to:A_in,B_in,s_in) map(from:out_loop) &
    !$omp   private(i,k)
    do i = 1, N
        do k = 1, M
            out_loop(i, k) = A_in(i, k) + s_in(i)*(B_in(i, k) - A_in(i, k))
        end do
    end do
    !$omp end target teams distribute parallel do

    ! Reference on CPU
    do i = 1, N
        do k = 1, M
            ref(i, k) = A_in(i, k) + s_in(i)*(B_in(i, k) - A_in(i, k))
        end do
    end do

    nerr = 0
    do i = 1, N
        do k = 1, M
            if (abs(out_whole(i,k) - ref(i,k)) > 1.e-12_wp) nerr = nerr + 1
        end do
    end do
    if (nerr == 0) then
        print *, "PASS test17a: whole-array op A+s*(B-A) in target loop"
    else
        print *, "FAIL test17a:", nerr, "wrong -- whole-array op broken"
        print *, "  (i=1,k=1): got", out_whole(1,1), "ref", ref(1,1)
    end if

    nerr = 0
    do i = 1, N
        do k = 1, M
            if (abs(out_loop(i,k) - ref(i,k)) > 1.e-12_wp) nerr = nerr + 1
        end do
    end do
    if (nerr == 0) then
        print *, "PASS test17b: explicit loop workaround in target loop"
    else
        print *, "FAIL test17b:", nerr, "wrong -- explicit loop also broken"
    end if

end program test17_whole_array_ops
