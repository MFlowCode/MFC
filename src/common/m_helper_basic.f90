!>
!! @file m_helper_basic.f90
!! @brief Contains module m_helper_basic

module m_helper_basic

    use m_derived_types        !< Definitions of the derived types

    implicit none

    private; 
    public :: f_approx_equal, &
              f_is_default, &
              f_all_default, &
              f_is_integer

contains

    !> This procedure checks if two floating point numbers of wp are within tolerance.
    !! @param a First number.
    !! @param b Second number.
    !! @param tol_input Relative error (default = 1e-6_wp).
    !! @return Result of the comparison.
    logical function f_approx_equal(a, b, tol_input) result(res)
        !$acc routine seq
        real(wp), intent(in) :: a, b
        real(wp), optional, intent(in) :: tol_input
        real(wp) :: tol

        if (present(tol_input)) then
            tol = tol_input
        else
            tol = 1e-6_wp
        end if

        if (a == b) then
            res = .true.
        else if (a == 0._wp .or. b == 0._wp .or. (abs(a) + abs(b) < tiny(a))) then
            res = (abs(a - b) < (tol*tiny(a)))
        else
            res = (abs(a - b)/min(abs(a) + abs(b), huge(a)) < tol)
        end if
    end function f_approx_equal

    !> Checks if a real(wp) variable is of default value.
    !! @param var Variable to check.
    logical function f_is_default(var) result(res)
        !$acc routine seq
        real(wp), intent(in) :: var

        res = f_approx_equal(var, dflt_real)
    end function f_is_default

    !> Checks if ALL elements of a real(wp) array are of default value.
    !! @param var_array Array to check.
    logical function f_all_default(var_array) result(res)
        real(wp), intent(in) :: var_array(:)
        logical :: res_array(size(var_array))
        integer :: i

        do i = 1, size(var_array)
            res_array(i) = f_is_default(var_array(i))
        end do

        res = all(res_array)
    end function f_all_default

    !> Checks if a real(wp) variable is an integer.
    !! @param var Variable to check.
    logical function f_is_integer(var) result(res)
        !$acc routine seq
        real(wp), intent(in) :: var

        res = f_approx_equal(var, real(nint(var), wp))
    end function f_is_integer

end module m_helper_basic
