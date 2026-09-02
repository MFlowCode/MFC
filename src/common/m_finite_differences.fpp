!>
!! @file
!! @brief Contains module m_finite_differences

#:include 'macros.fpp'

!> @brief Finite difference operators for computing divergence of velocity fields
module m_finite_differences

    use m_global_parameters

    implicit none

contains

    !> Compute the centered finite-difference coefficients for first-order spatial derivatives in the s-coordinate direction (x, y,
    !! or z). Supports up to 4th order accuracy.
    !! @param fd_coeff_s Finite-diff. coefficients in the s-coordinate direction
    !! @param local_buff_size Size of the local buffer
    !! @param fd_number_in Finite-difference number
    !! @param fd_order_in Finite-difference order of accuracy
    !! @param offset_s Optional offset bounds in the s-coordinate direction
    subroutine s_compute_finite_difference_coefficients(q, s_cc, fd_coeff_s, local_buff_size, fd_number_in, fd_order_in, offset_s)

        integer                                                               :: lB, lE  !< loop bounds
        integer, intent(in)                                                   :: q
        integer, intent(in)                                                   :: local_buff_size, fd_number_in, fd_order_in
        type(int_bounds_info), optional, intent(in)                           :: offset_s
        real(wp), allocatable, dimension(:,:), intent(inout)                  :: fd_coeff_s
        real(wp), dimension(-local_buff_size:q + local_buff_size), intent(in) :: s_cc
        integer                                                               :: i       !< Generic loop iterator

        ! Coefficients always extend at least fd_number_in beyond the interior on each side, so a stencil centered on a
        ! ghost-adjacent cell (e.g. an immersed boundary near a domain boundary) has a real coefficient to read instead of
        ! reading past the caller's allocation. offset_s, when given, widens this further (never narrows it) for callers
        ! that need more than fd_number_in of margin.

        if (present(offset_s)) then
            lB = -max(fd_number_in, offset_s%beg)
            lE = q + max(fd_number_in, offset_s%end)
        else
            lB = -fd_number_in
            lE = q + fd_number_in
        end if

        ! Computing the 1st order finite-difference coefficients
        if (fd_order_in == 1) then
            do i = lB, lE
                fd_coeff_s(-1, i) = 0._wp
                fd_coeff_s(0, i) = -1._wp/(s_cc(i + 1) - s_cc(i))
                fd_coeff_s(1, i) = -fd_coeff_s(0, i)
            end do

            ! Computing the 2nd order finite-difference coefficients
        else if (fd_order_in == 2) then
            do i = lB, lE
                fd_coeff_s(-1, i) = -1._wp/(s_cc(i + 1) - s_cc(i - 1))
                fd_coeff_s(0, i) = 0._wp
                fd_coeff_s(1, i) = -fd_coeff_s(-1, i)
            end do

            ! Computing the 4th order finite-difference coefficients
        else
            do i = lB, lE
                fd_coeff_s(-2, i) = 1._wp/(s_cc(i - 2) - 8._wp*s_cc(i - 1) - s_cc(i + 2) + 8._wp*s_cc(i + 1))
                fd_coeff_s(-1, i) = -8._wp*fd_coeff_s(-2, i)
                fd_coeff_s(0, i) = 0._wp
                fd_coeff_s(1, i) = -fd_coeff_s(-1, i)
                fd_coeff_s(2, i) = -fd_coeff_s(-2, i)
            end do
        end if

    end subroutine s_compute_finite_difference_coefficients

end module m_finite_differences
