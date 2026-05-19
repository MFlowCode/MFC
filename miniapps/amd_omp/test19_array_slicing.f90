! test19_array_slicing.f90
! Tests array slice assignments inside a target loop body.
! Patterns like:
!   F%L(3:4) = U%L(2)*vel%L(2:3) - B%L(1)*B%L(2:3)
!   F%L(5:6) = vel%L(1)*B%L(2:3) - vel%L(2:3)*B%L(1)
! were replaced with explicit element assignments in MFC's HLLD solver.
program test19_array_slicing
    implicit none
    integer, parameter :: wp = 8
    integer, parameter :: N = 2000

    real(wp) :: u_in(N), v_in(N), w_in(N), bx_in(N), by_in(N), bz_in(N)
    real(wp) :: F_slice(N, 6), F_explicit(N, 6), ref(N, 6)
    integer  :: i, nerr

    do i = 1, N
        u_in(i) = 0.5_wp + 0.001_wp*real(i, wp)
        v_in(i) = 0.2_wp - 0.0005_wp*real(i, wp)
        w_in(i) = 0.1_wp + 0.0003_wp*real(i, wp)
        bx_in(i) = 0.4_wp
        by_in(i) = 0.3_wp + 0.0001_wp*real(i, wp)
        bz_in(i) = 0.2_wp - 0.0001_wp*real(i, wp)
    end do

    ! Test A: slice assignments in target loop
    !$omp target teams distribute parallel do &
    !$omp   map(to:u_in,v_in,w_in,bx_in,by_in,bz_in) map(from:F_slice)
    do i = 1, N
        F_slice(i, 1) = u_in(i)
        F_slice(i, 2) = u_in(i)**2 + bx_in(i)**2
        F_slice(i, 3:4) = u_in(i)*[v_in(i), w_in(i)] - bx_in(i)*[by_in(i), bz_in(i)]
        F_slice(i, 5:6) = u_in(i)*[by_in(i), bz_in(i)] - [v_in(i), w_in(i)]*bx_in(i)
    end do
    !$omp end target teams distribute parallel do

    ! Test B: explicit elements (our workaround)
    !$omp target teams distribute parallel do &
    !$omp   map(to:u_in,v_in,w_in,bx_in,by_in,bz_in) map(from:F_explicit)
    do i = 1, N
        F_explicit(i, 1) = u_in(i)
        F_explicit(i, 2) = u_in(i)**2 + bx_in(i)**2
        F_explicit(i, 3) = u_in(i)*v_in(i) - bx_in(i)*by_in(i)
        F_explicit(i, 4) = u_in(i)*w_in(i) - bx_in(i)*bz_in(i)
        F_explicit(i, 5) = u_in(i)*by_in(i) - v_in(i)*bx_in(i)
        F_explicit(i, 6) = u_in(i)*bz_in(i) - w_in(i)*bx_in(i)
    end do
    !$omp end target teams distribute parallel do

    ! Reference on CPU
    do i = 1, N
        ref(i, 1) = u_in(i)
        ref(i, 2) = u_in(i)**2 + bx_in(i)**2
        ref(i, 3) = u_in(i)*v_in(i) - bx_in(i)*by_in(i)
        ref(i, 4) = u_in(i)*w_in(i) - bx_in(i)*bz_in(i)
        ref(i, 5) = u_in(i)*by_in(i) - v_in(i)*bx_in(i)
        ref(i, 6) = u_in(i)*bz_in(i) - w_in(i)*bx_in(i)
    end do

    nerr = 0
    do i = 1, N
        if (any(abs(F_slice(i,:) - ref(i,:)) > 1.e-12_wp)) nerr = nerr + 1
    end do
    if (nerr == 0) then
        print *, "PASS test19a: slice assignment F(3:4) = ... in target loop"
    else
        print *, "FAIL test19a:", nerr, "wrong -- slice assignment broken"
        print *, "  cell 1: got", F_slice(1,:)
        print *, "  cell 1: ref", ref(1,:)
    end if

    nerr = 0
    do i = 1, N
        if (any(abs(F_explicit(i,:) - ref(i,:)) > 1.e-12_wp)) nerr = nerr + 1
    end do
    if (nerr == 0) then
        print *, "PASS test19b: explicit element assignment workaround"
    else
        print *, "FAIL test19b:", nerr, "wrong -- explicit assignment also broken"
    end if

end program test19_array_slicing
