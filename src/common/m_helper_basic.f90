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
              f_is_integer, &
              s_configure_coordinate_bounds

contains

    !> This procedure checks if two floating point numbers of wp are within tolerance.
    !! @param a First number.
    !! @param b Second number.
    !! @param tol_input Relative error (default = 1e-6_wp).
    !! @return Result of the comparison.
    logical pure elemental function f_approx_equal(a, b, tol_input) result(res)
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
    logical pure elemental function f_is_default(var) result(res)
        !$acc routine seq
        real(wp), intent(in) :: var

        res = f_approx_equal(var, dflt_real)
    end function f_is_default

    !> Checks if ALL elements of a real(wp) array are of default value.
    !! @param var_array Array to check.
    logical pure function f_all_default(var_array) result(res)
        real(wp), intent(in) :: var_array(:)
        ! logical :: res_array(size(var_array))
        ! integer :: i

        res = all(f_is_default(var_array))

        ! do i = 1, size(var_array)
        !     res_array(i) = f_is_default(var_array(i))
        ! end do

        ! res = all(res_array)
    end function f_all_default

    !> Checks if a real(wp) variable is an integer.
    !! @param var Variable to check.
    logical pure elemental function f_is_integer(var) result(res)
        !$acc routine seq
        real(wp), intent(in) :: var

        res = f_approx_equal(var, real(nint(var), wp))
    end function f_is_integer

    pure subroutine s_configure_coordinate_bounds(weno_polyn, buff_size, idwint, idwbuff, &
                                                  viscous, bubbles_lagrange, m, n, p, num_dims)

        integer, intent(in) :: weno_polyn, m, n, p, num_dims
        integer, intent(inout) :: buff_size
        type(int_bounds_info), dimension(3), intent(inout) :: idwint, idwbuff
        logical, intent(in) :: viscous, bubbles_lagrange

        ! Determining the number of cells that are needed in order to store
        ! sufficient boundary conditions data as to iterate the solution in
        ! the physical computational domain from one time-step iteration to
        ! the next one
        if (viscous) then
            buff_size = 2*weno_polyn + 2
        else
            buff_size = weno_polyn + 2
        end if

        ! Correction for smearing function in the lagrangian subgrid bubble model
        if (bubbles_lagrange) then
            buff_size = max(buff_size, 6)
        end if

        ! Configuring Coordinate Direction Indexes
        idwint(1)%beg = 0; idwint(2)%beg = 0; idwint(3)%beg = 0
        idwint(1)%end = m; idwint(2)%end = n; idwint(3)%end = p

        idwbuff(1)%beg = -buff_size
        if (num_dims > 1) then; idwbuff(2)%beg = -buff_size; else; idwbuff(2)%beg = 0; end if
        if (num_dims > 2) then; idwbuff(3)%beg = -buff_size; else; idwbuff(3)%beg = 0; end if

        idwbuff(1)%end = idwint(1)%end - idwbuff(1)%beg
        idwbuff(2)%end = idwint(2)%end - idwbuff(2)%beg
        idwbuff(3)%end = idwint(3)%end - idwbuff(3)%beg

    end subroutine s_configure_coordinate_bounds

end module m_helper_basic
