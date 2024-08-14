!>
!! @file m_helper_basic.f90
!! @brief Contains module m_helper_basic

module m_helper_basic

    ! Dependencies =============================================================

    use m_derived_types        !< Definitions of the derived types
    ! ==========================================================================

    implicit none

    private; 
    public :: f_approx_equal, &
              f_is_default, &
              f_all_default

contains

    !> This procedure checks if two floating point numbers of kind(0d0) are within tolerance.
    !! @param a First number.
    !! @param b Second number.
    !! @param tol_input Relative error (default = 1d-6).
    !! @return Result of the comparison.
    logical function f_approx_equal(a, b, tol_input) result(res)
        ! Reference: https://floating-point-gui.de/errors/comparison/

        real(kind(0d0)), intent(in) :: a, b
        real(kind(0d0)), optional, intent(in) :: tol_input
        real(kind(0d0)) :: tol

        if (present(tol_input)) then
            tol = tol_input
        else
            tol = 1d-6
        end if

        if (a == b) then
            res = .true.
        else if (a == 0d0 .or. b == 0d0 .or. (abs(a) + abs(b) < tiny(a))) then
            res = (abs(a - b) < (tol*tiny(a)))
        else
            res = (abs(a - b)/min(abs(a) + abs(b), huge(a)) < tol)
        end if
    end function f_approx_equal

    !> Checks if a real(kind(0d0)) variable is of default value.
    !! @param var Variable to check.
    logical function f_is_default(var) result(res)
        !$acc routine seq
        real(kind(0d0)), intent(in) :: var

        res = f_approx_equal(var, dflt_real)
    end function f_is_default

    !> Checks if ALL elements of a real(kind(0d0)) array are of default value.
    !! @param var_array Array to check.
    logical function f_all_default(var_array) result(res)
        real(kind(0d0)), intent(in) :: var_array(:)
        logical :: res_array(size(var_array))
        integer :: i

        do i = 1, size(var_array)
            res_array(i) = f_is_default(var_array(i))
        end do

        res = all(res_array)
    end function f_all_default

end module m_helper_basic
