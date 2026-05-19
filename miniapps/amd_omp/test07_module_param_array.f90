! test07_module_param_array.f90
!
! Demonstrates AMD flang 7.2.0 (amdflang) device-code failure when a Fortran
! `parameter` array declared in a `use`d module is accessed as a WHOLE ARRAY
! inside an OpenMP target offload region.
!
! Root cause in MFC: `molecular_weights(:)` from m_thermochem (a `parameter`
! array, not GPU_DECLARE'd) is used as a whole-array operand inside GPU loops
! after a merge commit.  This silently drops ALL GPU kernels during
! device LTO, producing HSA_STATUS_ERROR_INVALID_SYMBOL_NAME at runtime.
!
! This miniapp has two tests:
!   PART A - whole-array access  `param_arr(:) * A(:)` in target region  →  FAILS
!   PART B - indexed access      `param_arr(i)  * A(i)` in do-loop        →  WORKS
!
! Expected results when AMD bug is present:
!   PART A: linker error (_FortranAAssign or silent kernel drop)
!   PART B: PASS, printed result is correct

module m_constants
    implicit none
    real(8), parameter :: param_arr(5) = [2.016d0, 31.998d0, 18.015d0, 28.013d0, 44.010d0]
end module m_constants

program test07_module_param_array
    use m_constants, only: param_arr
    implicit none

    integer, parameter :: N = 5
    real(8) :: A(N), out_a(N), out_b(N)
    integer :: i

    ! initialise input array on host
    do i = 1, N
        A(i) = real(i, 8)
    end do

    out_a = 0.0d0
    out_b = 0.0d0

    ! -----------------------------------------------------------------------
    ! PART A: whole-array access to module parameter inside target region
    !         AMD flang 7.2.0 fails here (silent LTO drop or link error)
    ! -----------------------------------------------------------------------
    !$omp target map(to:A) map(from:out_a)
    out_a(:) = A(:) / param_arr(:)
    !$omp end target

    write(*,'(a)') "PART A (whole-array) result:"
    do i = 1, N
        write(*,'(2x,i2,a,f12.6)') i, ": ", out_a(i)
    end do

    ! -----------------------------------------------------------------------
    ! PART B: indexed access to module parameter inside target teams loop
    !         This is the workaround used in MFC after the fix
    ! -----------------------------------------------------------------------
    !$omp target teams distribute parallel do map(to:A) map(from:out_b)
    do i = 1, N
        out_b(i) = A(i) / param_arr(i)
    end do
    !$omp end target teams distribute parallel do

    write(*,'(a)') "PART B (indexed loop) result:"
    do i = 1, N
        write(*,'(2x,i2,a,f12.6)') i, ": ", out_b(i)
    end do

    write(*,'(a)') "PASS"
end program test07_module_param_array
