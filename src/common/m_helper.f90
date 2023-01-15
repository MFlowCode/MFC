!>
!! @file m_helper.f90
!! @brief Contains module m_helper
module m_helper

    ! Dependencies =============================================================
    
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    ! ==========================================================================

    implicit none

    real(kind(0d0)) :: cart_y, cart_z
    real(kind(0d0)) :: sph_phi !<
    !! Variables to be used to hold cell locations in Cartesian coordinates if
    !! 3D simulation is using cylindrical coordinates

    private; public :: s_convert_cylindrical_to_cartesian_coord, &
        s_convert_cylindrical_to_spherical_coord, &
        f_r, &
        s_compute_finite_difference_coefficients, &
        s_apply_scalar_divergence_theorem, &
        s_compute_fd_gradient, &
        s_comp_n_from_cons, &
        s_comp_n_from_prim, &
        s_quad

contains

    subroutine s_convert_cylindrical_to_cartesian_coord(cyl_y, cyl_z)

        real(kind(0d0)), intent(IN) :: cyl_y, cyl_z

        cart_y = cyl_y*sin(cyl_z)
        cart_z = cyl_y*cos(cyl_z)

    end subroutine s_convert_cylindrical_to_cartesian_coord ! --------------

    subroutine s_convert_cylindrical_to_spherical_coord(cyl_x, cyl_y)

        real(kind(0d0)), intent(IN) :: cyl_x, cyl_y

        sph_phi = atan(cyl_y/cyl_x)

    end subroutine s_convert_cylindrical_to_spherical_coord ! --------------

    !> Archimedes spiral function
        !! @param myth Angle
        !! @param offset Thickness
        !! @param a Starting position
    function f_r(myth, offset, a)
        real(kind(0d0)), intent(IN) :: myth, offset, a
        real(kind(0d0)) :: b
        real(kind(0d0)) :: f_r

        !r(th) = a + b*th

        b = 2.d0*a/(2.d0*pi)
        f_r = a + b*myth + offset
    end function f_r

        !>  The purpose of this subroutine is to employ the inputted
        !!      left and right cell-boundary integral-averaged variables
        !!      to compute the relevant cell-average first-order spatial
        !!      derivatives in the x-, y- or z-direction by means of the
        !!      scalar divergence theorem.
        !!  @param vL_vf Left cell-boundary integral averages
        !!  @param vR_vf Right cell-boundary integral averages
        !!  @param dv_ds_vf Cell-average first-order spatial derivatives
        !!  @param norm_dir Splitting coordinate direction
    subroutine s_apply_scalar_divergence_theorem(vL_vf, vR_vf, & ! --------
                                                 dv_ds_vf, &
                                                 norm_dir, &
                                                 ix, iy, iz, iv, &
                                                 dxL, dyL, dzL, buff_size)

        type(int_bounds_info) :: ix, iy, iz, iv
            
        integer :: buff_size

        real(kind(0d0)), dimension(-buff_size:m + buff_size) :: dxL
        real(kind(0d0)), dimension(-buff_size:n + buff_size) :: dyL
        real(kind(0d0)), dimension(-buff_size:p + buff_size) :: dzL
        ! arrays of cell widths

        type(scalar_field), &
            dimension(iv%beg:iv%end), &
            intent(IN) :: vL_vf, vR_vf

        type(scalar_field), &
            dimension(iv%beg:iv%end), &
            intent(INOUT) :: dv_ds_vf

        integer, intent(IN) :: norm_dir

        integer :: i, j, k, l !< Generic loop iterators

        !$acc update device(ix, iy, iz, iv)

        ! First-Order Spatial Derivatives in x-direction ===================
        if (norm_dir == 1) then

            ! A general application of the scalar divergence theorem that
            ! utilizes the left and right cell-boundary integral-averages,
            ! inside each cell, or an arithmetic mean of these two at the
            ! cell-boundaries, to calculate the cell-averaged first-order
            ! spatial derivatives inside the cell.

!$acc parallel loop collapse(3) gang vector default(present)
            do l = iz%beg, iz%end
                do k = iy%beg, iy%end
                    do j = ix%beg + 1, ix%end - 1
!$acc loop seq
                        do i = iv%beg, iv%end

                            dv_ds_vf(i)%sf(j, k, l) = &
                                1d0/dxL(j) &
                                *( &
                                vR_vf(i)%sf(j, k, l) &
                                - vL_vf(i)%sf(j, k, l) &
                                )
                        end do
                    end do
                end do
            end do

            ! END: First-Order Spatial Derivatives in x-direction ==============

            ! First-Order Spatial Derivatives in y-direction ===================
        elseif (norm_dir == 2) then

            ! A general application of the scalar divergence theorem that
            ! utilizes the left and right cell-boundary integral-averages,
            ! inside each cell, or an arithmetic mean of these two at the
            ! cell-boundaries, to calculate the cell-averaged first-order
            ! spatial derivatives inside the cell.

!$acc parallel loop collapse(3) gang vector default(present)

            do l = iz%beg, iz%end
                do k = iy%beg + 1, iy%end - 1
                    do j = ix%beg, ix%end
!$acc loop seq
                        do i = iv%beg, iv%end
                            dv_ds_vf(i)%sf(j, k, l) = &
                                1d0/dyL(k) &
                                *( &
                                vR_vf(i)%sf(j, k, l) &
                                - vL_vf(i)%sf(j, k, l) &
                                )
                        end do
                    end do
                end do
            end do

            ! END: First-Order Spatial Derivatives in y-direction ==============

            ! First-Order Spatial Derivatives in z-direction ===================
        else

            ! A general application of the scalar divergence theorem that
            ! utilizes the left and right cell-boundary integral-averages,
            ! inside each cell, or an arithmetic mean of these two at the
            ! cell-boundaries, to calculate the cell-averaged first-order
            ! spatial derivatives inside the cell.

!$acc parallel loop collapse(3) gang vector default(present)
            do l = iz%beg + 1, iz%end - 1
                do k = iy%beg, iy%end
                    do j = ix%beg, ix%end
!$acc loop seq
                        do i = iv%beg, iv%end
                            dv_ds_vf(i)%sf(j, k, l) = &
                                1d0/dzL(l) &
                                *( &
                                vR_vf(i)%sf(j, k, l) &
                                - vL_vf(i)%sf(j, k, l) &
                                )
                        end do
                    end do
                end do
            end do

        end if
        ! END: First-Order Spatial Derivatives in z-direction ==============

    end subroutine s_apply_scalar_divergence_theorem ! ---------------------

    !>  Computes the scalar gradient fields via finite differences
        !!  @param var Variable to compute derivative of
        !!  @param grad_x First coordinate direction component of the derivative
        !!  @param grad_y Second coordinate direction component of the derivative
        !!  @param grad_z Third coordinate direction component of the derivative
        !!  @param norm Norm of the gradient vector
    subroutine s_compute_fd_gradient(var, grad_x, grad_y, grad_z, norm, &
                                     ix, iy, iz, buff_size)

        type(scalar_field), intent(IN) :: var
        type(scalar_field), intent(INOUT) :: grad_x
        type(scalar_field), intent(INOUT) :: grad_y
        type(scalar_field), intent(INOUT) :: grad_z
        type(scalar_field), intent(INOUT) :: norm

        integer, intent(IN) :: buff_size

        integer :: j, k, l !< Generic loop iterators

        type(int_bounds_info) :: ix, iy, iz

        ix%beg = -buff_size; ix%end = m + buff_size; 
        if (n > 0) then
            iy%beg = -buff_size; iy%end = n + buff_size
        else
            iy%beg = -1; iy%end = 1
        end if

        if (p > 0) then
            iz%beg = -buff_size; iz%end = p + buff_size
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

        ix%beg = -buff_size; ix%end = m + buff_size; 
        if (n > 0) then
            iy%beg = -buff_size; iy%end = n + buff_size
        else
            iy%beg = 0; iy%end = 0
        end if

        if (p > 0) then
            iz%beg = -buff_size; iz%end = p + buff_size
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
            if (bc_y%beg <= -3 .and. bc_y%beg /= -13) then
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
                                                                fd_number, fd_order)

        integer, intent(IN) :: q
        integer, intent(IN) :: buff_size, fd_number, fd_order

        real(kind(0d0)), &
            dimension(-buff_size:q + buff_size), &
            intent(IN) :: s_cc

        real(kind(0d0)), &
            dimension(-fd_number:fd_number, 0:q), &
            intent(INOUT) :: fd_coeff_s

        integer :: i !< Generic loop iterator

        ! Computing the 1st order finite-difference coefficients
        if (fd_order == 1) then
            do i = 0, q
                fd_coeff_s(-1, i) = 0d0
                fd_coeff_s(0, i) = -1d0/(s_cc(i + 1) - s_cc(i))
                fd_coeff_s(1, i) = -fd_coeff_s(0, i)
            end do

            ! Computing the 2nd order finite-difference coefficients
        elseif (fd_order == 2) then
            do i = 0, q
                fd_coeff_s(-1, i) = -1d0/(s_cc(i + 1) - s_cc(i - 1))
                fd_coeff_s(0, i) = 0d0
                fd_coeff_s(1, i) = -fd_coeff_s(-1, i)
            end do

            ! Computing the 4th order finite-difference coefficients
        else
            do i = 0, q
                fd_coeff_s(-2, i) = 1d0/(s_cc(i - 2) - 8d0*s_cc(i - 1) - s_cc(i + 2) + 8d0*s_cc(i + 1))
                fd_coeff_s(-1, i) = -8d0*fd_coeff_s(-2, i)
                fd_coeff_s(0, i) = 0d0
                fd_coeff_s(1, i) = -fd_coeff_s(-1, i)
                fd_coeff_s(2, i) = -fd_coeff_s(-2, i)
            end do

        end if

    end subroutine s_compute_finite_difference_coefficients ! --------------

    !> Computes the bubble number density n from the conservative variables
        !! @param vftmp is the void fraction
        !! @param nRtmp is the bubble number  density times the bubble radii
        !! @param ntmp is the output number bubble density
    subroutine s_comp_n_from_cons(vftmp, nRtmp, ntmp, weight)
    !$acc routine seq
        real(kind(0.d0)), intent(IN) :: vftmp
        real(kind(0.d0)), dimension(nb), intent(IN) :: nRtmp
        real(kind(0.d0)), intent(OUT) :: ntmp
        real(kind(0.d0)) :: nR3
        real(kind(0.d0)), dimension(nb) :: weight

        call s_quad(nRtmp**3.d0, nR3, weight)  !returns itself if NR0 = 1
        ntmp = DSQRT((4.d0*pi/3.d0)*nR3/vftmp)
        ! ntmp = 1d0

    end subroutine s_comp_n_from_cons

    !> Computes the bubble number density n from the primitive variables
        !! @param vftmp is the void fraction
        !! @param Rtmp is the  bubble radii
        !! @param ntmp is the output number bubble density
    subroutine s_comp_n_from_prim(vftmp, Rtmp, ntmp, weight)
    !$acc routine seq
        real(kind(0.d0)), intent(IN) :: vftmp
        real(kind(0.d0)), dimension(nb), intent(IN) :: Rtmp
        real(kind(0.d0)), intent(OUT) :: ntmp
        real(kind(0.d0)) :: R3
        real(kind(0.d0)), dimension(nb) :: weight

        call s_quad(Rtmp**3.d0, R3, weight)  !returns itself if NR0 = 1
        ntmp = (3.d0/(4.d0*pi))*vftmp/R3
        ! ntmp = 1d0

    end subroutine s_comp_n_from_prim

    !> Computes the quadrature for polydisperse bubble populations
        !! @param func is the bubble dynamic variables for each bin
        !! @param mom is the computed moment
    subroutine s_quad(func, mom, weight)

        real(kind(0.d0)), dimension(nb), intent(IN) :: func
        real(kind(0.d0)), intent(OUT) :: mom
        real(kind(0.d0)), dimension(nb) :: weight

        mom = dot_product(weight, func)

    end subroutine s_quad

end module m_helper
