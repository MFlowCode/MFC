module m_finite_differences

    use m_global_parameters

    implicit none

contains

    !>  Computes the scalar gradient fields via finite differences
    !!  @param var Variable to compute derivative of
    !!  @param grad_x First coordinate direction component of the derivative
    !!  @param grad_y Second coordinate direction component of the derivative
    !!  @param grad_z Third coordinate direction component of the derivative
    !!  @param norm Norm of the gradient vector
    subroutine s_compute_fd_gradient(var, grad_x, grad_y, grad_z, &
                                     ix, iy, iz, buff_size_in)

        type(scalar_field), intent(IN) :: var
        type(scalar_field), intent(INOUT) :: grad_x
        type(scalar_field), intent(INOUT) :: grad_y
        type(scalar_field), intent(INOUT) :: grad_z

        integer, intent(IN) :: buff_size_in

        integer :: j, k, l !< Generic loop iterators

        type(int_bounds_info) :: ix, iy, iz

        ix%beg = -buff_size_in; ix%end = m + buff_size_in; 
        if (n > 0) then
            iy%beg = -buff_size_in; iy%end = n + buff_size_in
        else
            iy%beg = -1; iy%end = 1
        end if

        if (p > 0) then
            iz%beg = -buff_size_in; iz%end = p + buff_size_in
        else
            iz%beg = -1; iz%end = 1
        end if

        !$acc update device(ix, iy, iz)

        !$acc parallel loop collapse(3) gang vector default(present)
        do l = iz%beg + 1, iz%end - 1
            do k = iy%beg + 1, iy%end - 1
                do j = ix%beg + 1, ix%end - 1
                    grad_x%sf(j, k, l) = &
                        (var%sf(j + 1, k, l) - var%sf(j - 1, k, l))/ &
                        (x_cc(j + 1) - x_cc(j - 1))
                end do
            end do
        end do

        if (n > 0) then
            !$acc parallel loop collapse(3) gang vector
            do l = iz%beg + 1, iz%end - 1
                do k = iy%beg + 1, iy%end - 1
                    do j = ix%beg + 1, ix%end - 1
                        grad_y%sf(j, k, l) = &
                            (var%sf(j, k + 1, l) - var%sf(j, k - 1, l))/ &
                            (y_cc(k + 1) - y_cc(k - 1))
                    end do
                end do
            end do
        end if

        if (p > 0) then
            !$acc parallel loop collapse(3) gang vector
            do l = iz%beg + 1, iz%end - 1
                do k = iy%beg + 1, iy%end - 1
                    do j = ix%beg + 1, ix%end - 1
                        grad_z%sf(j, k, l) = &
                            (var%sf(j, k, l + 1) - var%sf(j, k, l - 1))/ &
                            (z_cc(l + 1) - z_cc(l - 1))
                    end do
                end do
            end do
        end if

        ix%beg = -buff_size_in; ix%end = m + buff_size_in; 
        if (n > 0) then
            iy%beg = -buff_size_in; iy%end = n + buff_size_in
        else
            iy%beg = 0; iy%end = 0
        end if

        if (p > 0) then
            iz%beg = -buff_size_in; iz%end = p + buff_size_in
        else
            iz%beg = 0; iz%end = 0
        end if

        !$acc update device(ix, iy, iz)

        !$acc parallel loop collapse(2) gang vector default(present)
        do l = iz%beg, iz%end
            do k = iy%beg, iy%end
                grad_x%sf(ix%beg, k, l) = &
                    (-3d0*var%sf(ix%beg, k, l) + 4d0*var%sf(ix%beg + 1, k, l) - var%sf(ix%beg + 2, k, l))/ &
                    (x_cc(ix%beg + 2) - x_cc(ix%beg))
                grad_x%sf(ix%end, k, l) = &
                    (3d0*var%sf(ix%end, k, l) - 4d0*var%sf(ix%end - 1, k, l) + var%sf(ix%end - 2, k, l))/ &
                    (x_cc(ix%end) - x_cc(ix%end - 2))
            end do
        end do
        if (n > 0) then
            !$acc parallel loop collapse(2) gang vector default(present)
            do l = iz%beg, iz%end
                do j = ix%beg, ix%end
                    grad_y%sf(j, iy%beg, l) = &
                        (-3d0*var%sf(j, iy%beg, l) + 4d0*var%sf(j, iy%beg + 1, l) - var%sf(j, iy%beg + 2, l))/ &
                        (y_cc(iy%beg + 2) - y_cc(iy%beg))
                    grad_y%sf(j, iy%end, l) = &
                        (3d0*var%sf(j, iy%end, l) - 4d0*var%sf(j, iy%end - 1, l) + var%sf(j, iy%end - 2, l))/ &
                        (y_cc(iy%end) - y_cc(iy%end - 2))
                end do
            end do
            if (p > 0) then
                !$acc parallel loop collapse(2) gang vector default(present)
                do k = iy%beg, iy%end
                    do j = ix%beg, ix%end
                        grad_z%sf(j, k, iz%beg) = &
                            (-3d0*var%sf(j, k, iz%beg) + 4d0*var%sf(j, k, iz%beg + 1) - var%sf(j, k, iz%beg + 2))/ &
                            (z_cc(iz%beg + 2) - z_cc(iz%beg))
                        grad_z%sf(j, k, iz%end) = &
                            (3d0*var%sf(j, k, iz%end) - 4d0*var%sf(j, k, iz%end - 1) + var%sf(j, k, iz%end - 2))/ &
                            (z_cc(iz%end) - z_cc(iz%end - 2))
                    end do
                end do
            end if
        end if

        if (bc_x%beg <= -3) then
            !$acc parallel loop collapse(2) gang vector default(present)
            do l = iz%beg, iz%end
                do k = iy%beg, iy%end
                    grad_x%sf(0, k, l) = (-3d0*var%sf(0, k, l) + 4d0*var%sf(1, k, l) - var%sf(2, k, l))/ &
                                         (x_cc(2) - x_cc(0))
                end do
            end do
        end if
        if (bc_x%end <= -3) then
            !$acc parallel loop collapse(2) gang vector default(present)
            do l = iz%beg, iz%end
                do k = iy%beg, iy%end
                    grad_x%sf(m, k, l) = (3d0*var%sf(m, k, l) - 4d0*var%sf(m - 1, k, l) + var%sf(m - 2, k, l))/ &
                                         (x_cc(m) - x_cc(m - 2))
                end do
            end do
        end if
        if (n > 0) then
            if (bc_y%beg <= -3 .and. bc_y%beg /= -14) then
                !$acc parallel loop collapse(2) gang vector default(present)
                do l = iz%beg, iz%end
                    do j = ix%beg, ix%end
                        grad_y%sf(j, 0, l) = (-3d0*var%sf(j, 0, l) + 4d0*var%sf(j, 1, l) - var%sf(j, 2, l))/ &
                                             (y_cc(2) - y_cc(0))
                    end do
                end do
            end if
            if (bc_y%end <= -3) then
                !$acc parallel loop collapse(2) gang vector default(present)
                do l = iz%beg, iz%end
                    do j = ix%beg, ix%end
                        grad_y%sf(j, n, l) = (3d0*var%sf(j, n, l) - 4d0*var%sf(j, n - 1, l) + var%sf(j, n - 2, l))/ &
                                             (y_cc(n) - y_cc(n - 2))
                    end do
                end do
            end if
            if (p > 0) then
                if (bc_z%beg <= -3) then
                    !$acc parallel loop collapse(2) gang vector default(present)
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end
                            grad_z%sf(j, k, 0) = &
                                (-3d0*var%sf(j, k, 0) + 4d0*var%sf(j, k, 1) - var%sf(j, k, 2))/ &
                                (z_cc(2) - z_cc(0))
                        end do
                    end do
                end if
                if (bc_z%end <= -3) then
                    !$acc parallel loop collapse(2) gang vector default(present)
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end
                            grad_z%sf(j, k, p) = &
                                (3d0*var%sf(j, k, p) - 4d0*var%sf(j, k, p - 1) + var%sf(j, k, p - 2))/ &
                                (z_cc(p) - z_cc(p - 2))
                        end do
                    end do
                end if
            end if
        end if

    end subroutine s_compute_fd_gradient ! --------------------------------------

    subroutine s_compute_fd_laplacian(laplacian, field, ix, iy, iz)

        type(scalar_field), intent(INOUT) :: laplacian
        type(scalar_field), intent(IN) :: field
        type(int_bounds_info), intent(IN) :: ix, iy, iz

        integer :: j, k, l !< Generic loop iterators

        type(scalar_field) :: grad_x, grad_y, grad_z

        !return

        allocate (grad_x%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
        allocate (grad_y%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
        allocate (grad_z%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))

        call s_compute_fd_gradient(field, grad_x, grad_y, grad_z, ix, iy, iz, 0)

        do l = iz%beg, iz%end
            do k = iy%beg, iy%end
                do j = ix%beg + 1, ix%end - 1
                    block
                        real(kind(0d0)) :: result

                        result = &
                            (grad_x%sf(j + 1, k, l) - grad_x%sf(j - 1, k, l))/ &
                            (x_cc(j + 1) - x_cc(j - 1))

                        if (n > 0) then
                            result = result + (grad_y%sf(j, k + 1, l) - grad_y%sf(j, k - 1, l))/ &
                                     (y_cc(k + 1) - y_cc(k - 1))
                        end if

                        if (p > 0) then
                            result = result + &
                                     (grad_z%sf(j, k, l + 1) - grad_z%sf(j, k, l - 1))/ &
                                     (z_cc(l + 1) - z_cc(l - 1))
                        end if

                        laplacian%sf(j, k, l) = result
                    end block
                end do
            end do
        end do

        deallocate (grad_x%sf, grad_y%sf, grad_z%sf)

    end subroutine s_compute_fd_laplacian ! -----------------------------------

    subroutine s_compute_fd_divergence(div, fields, ix, iy, iz)

        type(scalar_field), intent(INOUT) :: div
        type(scalar_field), intent(IN) :: fields(1:3)
        type(int_bounds_info), intent(IN) :: ix, iy, iz

        integer :: x, y, z !< Generic loop iterators

        real(kind(0d0)) :: divergence

        !$acc parallel loop collapse(3) private(divergence)
        do x = ix%beg, ix%end
            do y = iy%beg, iy%end
                do z = iz%beg, iz%end

                    if (x == ix%beg) then
                        divergence = (-3d0*fields(1)%sf(x, y, z) + 4d0*fields(1)%sf(x + 1, y, z) - fields(1)%sf(x + 2, y, z))/(x_cc(x + 2) - x_cc(x))
                    else if (x == ix%end) then
                        divergence = (+3d0*fields(1)%sf(x, y, z) - 4d0*fields(1)%sf(x - 1, y, z) + fields(1)%sf(x - 2, y, z))/(x_cc(x) - x_cc(x - 2))
                    else
                        divergence = (fields(1)%sf(x + 1, y, z) - fields(1)%sf(x - 1, y, z))/(x_cc(x + 1) - x_cc(x - 1))
                    end if

                    if (n > 0) then
                        if (y == iy%beg) then
                            divergence = divergence + (-3d0*fields(2)%sf(x, y, z) + 4d0*fields(2)%sf(x, y + 1, z) - fields(2)%sf(x, y + 2, z))/(y_cc(y + 2) - y_cc(y))
                        else if (y == iy%end) then
                            divergence = divergence + (+3d0*fields(2)%sf(x, y, z) - 4d0*fields(2)%sf(x, y - 1, z) + fields(2)%sf(x, y - 2, z))/(y_cc(y) - y_cc(y - 2))
                        else
                            divergence = divergence + (fields(2)%sf(x, y + 1, z) - fields(2)%sf(x, y - 1, z))/(y_cc(y + 1) - y_cc(y - 1))
                        end if
                    end if

                    if (p > 0) then
                        if (z == iz%beg) then
                            divergence = divergence + (-3d0*fields(3)%sf(x, y, z) + 4d0*fields(3)%sf(x, y, z + 1) - fields(3)%sf(x, y, z + 2))/(z_cc(z + 2) - z_cc(z))
                        else if (z == iz%end) then
                            divergence = divergence + (+3d0*fields(3)%sf(x, y, z) - 4d0*fields(3)%sf(x, y, z - 1) + fields(2)%sf(x, y, z - 2))/(z_cc(z) - z_cc(z - 2))
                        else
                            divergence = divergence + (fields(3)%sf(x, y, z + 1) - fields(3)%sf(x, y, z - 1))/(z_cc(z + 1) - z_cc(z - 1))
                        end if
                    end if

                    div%sf(x, y, z) = div%sf(x, y, z) + divergence

                end do
            end do
        end do

    end subroutine s_compute_fd_divergence

    !>  The purpose of this subroutine is to compute the finite-
    !!      difference coefficients for the centered schemes utilized
    !!      in computations of first order spatial derivatives in the
    !!      s-coordinate direction. The s-coordinate direction refers
    !!      to the x-, y- or z-coordinate direction, depending on the
    !!      subroutine's inputs. Note that coefficients of up to 4th
    !!      order accuracy are available.
    !!  @param q Number of cells in the s-coordinate direction
    !!  @param s_cc Locations of the cell-centers in the s-coordinate direction
    !!  @param fd_coeff_s Finite-diff. coefficients in the s-coordinate direction
    subroutine s_compute_finite_difference_coefficients(q, s_cc, fd_coeff_s, buff_size, &
                                                        fd_number_in, fd_order_in, offset_s)

        integer :: lB, lE !< loop bounds
        integer, intent(IN) :: q
        integer, intent(IN) :: buff_size, fd_number_in, fd_order_in
        type(int_bounds_info), optional, intent(IN) :: offset_s
        real(kind(0d0)), allocatable, dimension(:, :), intent(INOUT) :: fd_coeff_s

        real(kind(0d0)), &
            dimension(-buff_size:q + buff_size), &
            intent(IN) :: s_cc

        integer :: i !< Generic loop iterator

        if (present(offset_s)) then
            lB = -offset_s%beg
            lE = q + offset_s%end
        else
            lB = 0
            lE = q
        end if

        if (allocated(fd_coeff_s)) deallocate (fd_coeff_s)
        allocate (fd_coeff_s(-fd_number_in:fd_number_in, lb:lE))

        ! Computing the 1st order finite-difference coefficients
        if (fd_order_in == 1) then
            do i = lB, lE
                fd_coeff_s(-1, i) = 0d0
                fd_coeff_s(0, i) = -1d0/(s_cc(i + 1) - s_cc(i))
                fd_coeff_s(1, i) = -fd_coeff_s(0, i)
            end do

            ! Computing the 2nd order finite-difference coefficients
        elseif (fd_order_in == 2) then
            do i = lB, lE
                fd_coeff_s(-1, i) = -1d0/(s_cc(i + 1) - s_cc(i - 1))
                fd_coeff_s(0, i) = 0d0
                fd_coeff_s(1, i) = -fd_coeff_s(-1, i)
            end do

            ! Computing the 4th order finite-difference coefficients
        else
            do i = lB, lE
                fd_coeff_s(-2, i) = 1d0/(s_cc(i - 2) - 8d0*s_cc(i - 1) - s_cc(i + 2) + 8d0*s_cc(i + 1))
                fd_coeff_s(-1, i) = -8d0*fd_coeff_s(-2, i)
                fd_coeff_s(0, i) = 0d0
                fd_coeff_s(1, i) = -fd_coeff_s(-1, i)
                fd_coeff_s(2, i) = -fd_coeff_s(-2, i)
            end do

        end if

    end subroutine s_compute_finite_difference_coefficients ! --------------

end module m_finite_differences
