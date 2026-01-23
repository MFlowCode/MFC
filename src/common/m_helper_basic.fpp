
!! @file m_helper_basic.f90
!! @brief Contains module m_helper_basic

#:include 'macros.fpp'

module m_helper_basic

    use m_derived_types        !< Definitions of the derived types

    implicit none

    private; 
    public :: f_approx_equal, &
              f_approx_in_array, &
              f_is_default, &
              f_all_default, &
              f_is_integer, &
              s_configure_coordinate_bounds, &
              s_update_cell_bounds

contains

    !> This procedure checks if two floating point numbers of wp are within tolerance.
    !! @param a First number.
    !! @param b Second number.
    !! @param tol_input Relative error (default = 1.e-10_wp).
    !! @return Result of the comparison.
    logical elemental function f_approx_equal(a, b, tol_input) result(res)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: a, b
        real(wp), optional, intent(in) :: tol_input
        real(wp) :: tol

        if (present(tol_input)) then
            tol = tol_input
        else
            tol = 1.e-10_wp
        end if

        if (a == b) then
            res = .true.
        else if (a == 0._wp .or. b == 0._wp .or. (abs(a) + abs(b) < tiny(a))) then
            res = (abs(a - b) < (tol*tiny(a)))
        else
            res = (abs(a - b)/min(abs(a) + abs(b), huge(a)) < tol)
        end if
    end function f_approx_equal

    !> This procedure checks if the point numbers of wp belongs to another array are within tolerance.
    !! @param a First number.
    !! @param b Array that contains several point numbers.
    !! @param tol_input Relative error (default = 1e-10_wp).
    !! @return Result of the comparison.
    logical function f_approx_in_array(a, b, tol_input) result(res)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: a
        real(wp), intent(in) :: b(:)
        real(wp), optional, intent(in) :: tol_input
        real(wp) :: tol
        integer :: i

        res = .false.

        if (present(tol_input)) then
            tol = tol_input
        else
            tol = 1e-10_wp
        end if

        do i = 1, size(b)
            if (f_approx_equal(a, b(i), tol)) then
                res = .true.
                exit
            end if
        end do
    end function f_approx_in_array

    !> Checks if a real(wp) variable is of default value.
    !! @param var Variable to check.
    logical elemental function f_is_default(var) result(res)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: var

        res = f_approx_equal(var, dflt_real)
    end function f_is_default

    !> Checks if ALL elements of a real(wp) array are of default value.
    !! @param var_array Array to check.
    logical function f_all_default(var_array) result(res)
        real(wp), intent(in) :: var_array(:)

        res = all(f_is_default(var_array))

        !logical :: res_array(size(var_array))
        !integer :: i

        ! do i = 1, size(var_array)
        !     res_array(i) = f_is_default(var_array(i))
        ! end do

        ! res = all(res_array)
    end function f_all_default

    !> Checks if a real(wp) variable is an integer.
    !! @param var Variable to check.
    logical elemental function f_is_integer(var) result(res)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: var

        res = f_approx_equal(var, real(nint(var), wp))
    end function f_is_integer

    subroutine s_configure_coordinate_bounds(recon_type, weno_polyn, muscl_polyn, &
                                             igr_order, buff_size, idwint, idwbuff, &
                                             viscous, bubbles_lagrange, m, n, p, num_dims, igr, ib)

        integer, intent(in) :: recon_type, weno_polyn, muscl_polyn
        integer, intent(in) :: m, n, p, num_dims, igr_order
        integer, intent(inout) :: buff_size
        type(int_bounds_info), dimension(3), intent(inout) :: idwint, idwbuff
        logical, intent(in) :: viscous, bubbles_lagrange
        logical, intent(in) :: igr
        logical, intent(in) :: ib

        ! Determining the number of cells that are needed in order to store
        ! sufficient boundary conditions data as to iterate the solution in
        ! the physical computational domain from one time-step iteration to
        ! the next one
        if (igr) then
            buff_size = (igr_order - 1)/2 + 2
        elseif (recon_type == WENO_TYPE) then
            if (viscous) then
                buff_size = 2*weno_polyn + 2
            else
                buff_size = weno_polyn + 2
            end if
        elseif (recon_type == MUSCL_TYPE) then
            buff_size = muscl_polyn + 2
        end if

        ! Correction for smearing function in the lagrangian subgrid bubble model
        if (bubbles_lagrange) then
            buff_size = max(buff_size, 6)
        end if

        if (ib) then
            buff_size = max(buff_size, 10)
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

    !> Updates the min and max number of cells in each set of axes
    !! @param bounds Min ans max values to update
    !! @param m Number of cells in x-axis
    !! @param n Number of cells in y-axis
    !! @param p Number of cells in z-axis
    elemental subroutine s_update_cell_bounds(bounds, m, n, p)
        type(cell_num_bounds), intent(out) :: bounds
        integer, intent(in) :: m, n, p

        bounds%mn_max = max(m, n)
        bounds%np_max = max(n, p)
        bounds%mp_max = max(m, p)
        bounds%mnp_max = max(m, n, p)
        bounds%mn_min = min(m, n)
        bounds%np_min = min(n, p)
        bounds%mp_min = min(m, p)
        bounds%mnp_min = min(m, n, p)

    end subroutine s_update_cell_bounds

end module m_helper_basic
