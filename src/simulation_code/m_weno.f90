!!       __  _______________
!!      /  |/  / ____/ ____/
!!     / /|_/ / /_  / /
!!    / /  / / __/ / /___
!!   /_/  /_/_/    \____/
!!
!!  This file is part of MFC.
!!
!!  MFC is the legal property of its developers, whose names
!!  are listed in the copyright file included with this source
!!  distribution.
!!
!!  MFC is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published
!!  by the Free Software Foundation, either version 3 of the license
!!  or any later version.
!!
!!  MFC is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with MFC (LICENSE).
!!  If not, see <http://www.gnu.org/licenses/>.

!>
!! @file m_weno.f90
!! @brief Contains module m_weno
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief  Weighted essentially non-oscillatory (WENO) reconstruction scheme
!!              that is supplemented with monotonicity preserving bounds (MPWENO)
!!              and a mapping function that boosts the accuracy of the non-linear
!!              weights (WENOM). MPWENO, see Balsara and Shu (2000), prevents the
!!              reconstructed values to lay outside the range set by the stencil,
!!              while WENOM, see Henrick et al. (2005), recovers the formal order
!!              of accuracy of the reconstruction at critical points. Please note
!!              that the basic WENO approach is implemented according to the work
!!              of Jiang and Shu (1996).
module m_weno

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_variables_conversion !< State variables type conversion procedures
    ! ==========================================================================

    implicit none

    private; public :: s_initialize_weno_module, s_weno, s_finalize_weno_module

    !> @name The cell-average variables that will be WENO-reconstructed. Formerly, they
    !! are stored in v_vf. However, they are transferred to v_rs_wsL and v_rs_wsR
    !! as to be reshaped (RS) and/or characteristically decomposed. The reshaping
    !! allows the WENO procedure to be independent of the coordinate direction of
    !! the reconstruction. Lastly, notice that the left (L) and right (R) results
    !! of the characteristic decomposition are stored in custom-constructed WENO-
    !! stencils (WS) that are annexed to each position of a given scalar field.
    !> @{
    type(vector_field), allocatable, dimension(:) :: v_rs_wsL, v_rs_wsR
    !> @}

    !> @name Left and right WENO-reconstructed values of the cell-average variables.
    !! Note that the reshaped property of the variables from which these were
    !! obtained, v_rs_wsL and v_rs_wsR, is initially kept. Once the reshaping
    !! is undone, the reconstructed values are moved into vL_vf and vR_vf.
    !> @{
    type(scalar_field), allocatable, dimension(:) :: vL_rs_vf, vR_rs_vf
    !> @}

    ! WENO Coefficients ========================================================

    !> @name Polynomial coefficients at the left and right cell-boundaries (CB) and at
    !! the left and right quadrature points (QP), in the x-, y- and z-directions.
    !! Note that the first dimension of the array identifies the polynomial, the
    !! second dimension identifies the position of its coefficients and the last
    !! dimension denotes the cell-location in the relevant coordinate direction.
    !> @{
    real(kind(0d0)), target, allocatable, dimension(:, :, :) :: poly_coef_cbL_x
    real(kind(0d0)), target, allocatable, dimension(:, :, :) :: poly_coef_cbL_y
    real(kind(0d0)), target, allocatable, dimension(:, :, :) :: poly_coef_cbL_z

    real(kind(0d0)), target, allocatable, dimension(:, :, :) :: poly_coef_cbR_x
    real(kind(0d0)), target, allocatable, dimension(:, :, :) :: poly_coef_cbR_y
    real(kind(0d0)), target, allocatable, dimension(:, :, :) :: poly_coef_cbR_z

    real(kind(0d0)), pointer, dimension(:, :, :) :: poly_coef_L => null()
    real(kind(0d0)), pointer, dimension(:, :, :) :: poly_coef_R => null()
    !> @}

    !> @name The ideal weights at the left and the right cell-boundaries and at the
    !! left and the right quadrature points, in x-, y- and z-directions. Note
    !! that the first dimension of the array identifies the weight, while the
    !! last denotes the cell-location in the relevant coordinate direction.
    !> @{
    real(kind(0d0)), target, allocatable, dimension(:, :) :: d_cbL_x
    real(kind(0d0)), target, allocatable, dimension(:, :) :: d_cbL_y
    real(kind(0d0)), target, allocatable, dimension(:, :) :: d_cbL_z

    real(kind(0d0)), target, allocatable, dimension(:, :) :: d_cbR_x
    real(kind(0d0)), target, allocatable, dimension(:, :) :: d_cbR_y
    real(kind(0d0)), target, allocatable, dimension(:, :) :: d_cbR_z

    real(kind(0d0)), pointer, dimension(:, :) :: d_L => null()
    real(kind(0d0)), pointer, dimension(:, :) :: d_R => null()
    !> @}

    !> @name Smoothness indicator coefficients in the x-, y-, and z-directions. Note
    !! that the first array dimension identifies the smoothness indicator, the
    !! second identifies the position of its coefficients and the last denotes
    !! the cell-location in the relevant coordinate direction.
    !> @{
    real(kind(0d0)), target, allocatable, dimension(:, :, :) :: beta_coef_x
    real(kind(0d0)), target, allocatable, dimension(:, :, :) :: beta_coef_y
    real(kind(0d0)), target, allocatable, dimension(:, :, :) :: beta_coef_z

    real(kind(0d0)), pointer, dimension(:, :, :) :: beta_coef => null()
    !> @}

    ! END: WENO Coefficients ===================================================

    integer :: v_size !< Number of WENO-reconstructed cell-average variables

    !> @name Indical bounds in the s1-, s2- and s3-directions
    !> @{
    type(bounds_info) :: is1, is2, is3
    !> @}

contains

    !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_weno_module() ! --------------------------------

        type(bounds_info) :: ix, iy, iz !< Indical bounds in the x-, y- and z-directions

        if (weno_order == 1) return

        ! Allocating WENO-stencil for the variables to be WENO-reconstructed
        allocate (v_rs_wsL(-weno_polyn:weno_polyn))
        allocate (v_rs_wsR(-weno_polyn:weno_polyn))

        ! Allocating/Computing WENO Coefficients in x-direction ============
        ix%beg = -buff_size + weno_polyn; ix%end = m - ix%beg

        allocate (poly_coef_cbL_x(0:weno_polyn, &
                                  0:weno_polyn - 1, &
                                  ix%beg:ix%end))
        allocate (poly_coef_cbR_x(0:weno_polyn, &
                                  0:weno_polyn - 1, &
                                  ix%beg:ix%end))

        allocate (d_cbL_x(0:weno_polyn, ix%beg:ix%end))
        allocate (d_cbR_x(0:weno_polyn, ix%beg:ix%end))

        allocate (beta_coef_x(0:weno_polyn, &
                              0:2*(weno_polyn - 1), &
                              ix%beg:ix%end))

        call s_compute_weno_coefficients(1, ix)

        ! ==================================================================

        ! Allocating/Computing WENO Coefficients in y-direction ============
        if (n == 0) return

        iy%beg = -buff_size + weno_polyn; iy%end = n - iy%beg

        allocate (poly_coef_cbL_y(0:weno_polyn, &
                                  0:weno_polyn - 1, &
                                  iy%beg:iy%end))
        allocate (poly_coef_cbR_y(0:weno_polyn, &
                                  0:weno_polyn - 1, &
                                  iy%beg:iy%end))

        allocate (d_cbL_y(0:weno_polyn, iy%beg:iy%end))
        allocate (d_cbR_y(0:weno_polyn, iy%beg:iy%end))

        allocate (beta_coef_y(0:weno_polyn, &
                              0:2*(weno_polyn - 1), &
                              iy%beg:iy%end))

        call s_compute_weno_coefficients(2, iy)

        ! ==================================================================

        ! Allocating/Computing WENO Coefficients in z-direction ============
        if (p == 0) return

        iz%beg = -buff_size + weno_polyn; iz%end = p - iz%beg

        allocate (poly_coef_cbL_z(0:weno_polyn, &
                                  0:weno_polyn - 1, &
                                  iz%beg:iz%end))
        allocate (poly_coef_cbR_z(0:weno_polyn, &
                                  0:weno_polyn - 1, &
                                  iz%beg:iz%end))

        allocate (d_cbL_z(0:weno_polyn, iz%beg:iz%end))
        allocate (d_cbR_z(0:weno_polyn, iz%beg:iz%end))

        allocate (beta_coef_z(0:weno_polyn, &
                              0:2*(weno_polyn - 1), &
                              iz%beg:iz%end))

        call s_compute_weno_coefficients(3, iz)

        ! ==================================================================

    end subroutine s_initialize_weno_module ! ------------------------------

    !>  The purpose of this subroutine is to compute the grid
        !!      dependent coefficients of the WENO polynomials, ideal
        !!      weights and smoothness indicators, provided the order,
        !!      the coordinate direction and the location of the WENO
        !!      reconstruction.
        !! @param weno_dir Coordinate direction of the WENO reconstruction
        !! @param is Index bounds in the s-direction
    subroutine s_compute_weno_coefficients(weno_dir, is) ! -------

        integer, intent(IN) :: weno_dir
        type(bounds_info), intent(IN) :: is
        integer :: s

        real(kind(0d0)), pointer, dimension(:) :: s_cb => null() !<
            !! Cell-boundary locations in the s-direction

        type(bounds_info) :: bc_s !< Boundary conditions (BC) in the s-direction

        integer :: i !< Generic loop iterator

        ! Associating WENO coefficients pointers
        call s_associate_weno_coefficients_pointers(weno_dir)

        ! Determining the number of cells, the cell-boundary locations and
        ! the boundary conditions in the coordinate direction selected for
        ! the WENO reconstruction
        if (weno_dir == 1) then
            s = m; s_cb => x_cb; bc_s = bc_x
        elseif (weno_dir == 2) then
            s = n; s_cb => y_cb; bc_s = bc_y
        else
            s = p; s_cb => z_cb; bc_s = bc_z
        end if

        ! Computing WENO3 Coefficients =====================================
        if (weno_order == 3) then

            do i = is%beg - 1, is%end - 1

                poly_coef_R(0, 0, i + 1) = (s_cb(i) - s_cb(i + 1))/ &
                                           (s_cb(i) - s_cb(i + 2))
                poly_coef_R(1, 0, i + 1) = (s_cb(i) - s_cb(i + 1))/ &
                                           (s_cb(i - 1) - s_cb(i + 1))

                poly_coef_L(0, 0, i + 1) = -poly_coef_R(0, 0, i + 1)
                poly_coef_L(1, 0, i + 1) = -poly_coef_R(1, 0, i + 1)

                d_R(0, i + 1) = (s_cb(i - 1) - s_cb(i + 1))/ &
                                (s_cb(i - 1) - s_cb(i + 2))
                d_L(0, i + 1) = (s_cb(i - 1) - s_cb(i))/ &
                                (s_cb(i - 1) - s_cb(i + 2))

                d_R(1, i + 1) = 1d0 - d_R(0, i + 1)
                d_L(1, i + 1) = 1d0 - d_L(0, i + 1)

                beta_coef(0, 0, i + 1) = 4d0*(s_cb(i) - s_cb(i + 1))**2d0/ &
                                         (s_cb(i) - s_cb(i + 2))**2d0
                beta_coef(1, 0, i + 1) = 4d0*(s_cb(i) - s_cb(i + 1))**2d0/ &
                                         (s_cb(i - 1) - s_cb(i + 1))**2d0

            end do



            ! Modifying the ideal weights coefficients in the neighborhood
            ! of beginning and end Riemann state extrapolation BC to avoid
            ! any contributions from outside of the physical domain during
            ! the WENO reconstruction
            if (null_weights) then
                if (bc_s%beg == -4) then
                    d_R(1, 0) = 0d0; d_R(0, 0) = 1d0
                    d_L(1, 0) = 0d0; d_L(0, 0) = 1d0
                end if

                if (bc_s%end == -4) then
                    d_R(0, s) = 0d0; d_R(1, s) = 1d0
                    d_L(0, s) = 0d0; d_L(1, s) = 1d0
                end if
            end if
            ! END: Computing WENO3 Coefficients ================================

            ! Computing WENO5 Coefficients =====================================
        else

            do i = is%beg - 1, is%end - 1

                poly_coef_R(0, 0, i + 1) = &
                    ((s_cb(i) - s_cb(i + 1))*(s_cb(i + 1) - s_cb(i + 2)))/ &
                    ((s_cb(i) - s_cb(i + 3))*(s_cb(i + 3) - s_cb(i + 1)))
                poly_coef_R(1, 0, i + 1) = &
                    ((s_cb(i - 1) - s_cb(i + 1))*(s_cb(i + 1) - s_cb(i)))/ &
                    ((s_cb(i - 1) - s_cb(i + 2))*(s_cb(i + 2) - s_cb(i)))
                poly_coef_R(1, 1, i + 1) = &
                    ((s_cb(i) - s_cb(i + 1))*(s_cb(i + 1) - s_cb(i + 2)))/ &
                    ((s_cb(i - 1) - s_cb(i + 1))*(s_cb(i - 1) - s_cb(i + 2)))
                poly_coef_R(2, 1, i + 1) = &
                    ((s_cb(i) - s_cb(i + 1))*(s_cb(i + 1) - s_cb(i - 1)))/ &
                    ((s_cb(i - 2) - s_cb(i))*(s_cb(i - 2) - s_cb(i + 1)))
                poly_coef_L(0, 0, i + 1) = &
                    ((s_cb(i + 1) - s_cb(i))*(s_cb(i) - s_cb(i + 2)))/ &
                    ((s_cb(i) - s_cb(i + 3))*(s_cb(i + 3) - s_cb(i + 1)))
                poly_coef_L(1, 0, i + 1) = &
                    ((s_cb(i) - s_cb(i - 1))*(s_cb(i) - s_cb(i + 1)))/ &
                    ((s_cb(i - 1) - s_cb(i + 2))*(s_cb(i) - s_cb(i + 2)))
                poly_coef_L(1, 1, i + 1) = &
                    ((s_cb(i + 1) - s_cb(i))*(s_cb(i) - s_cb(i + 2)))/ &
                    ((s_cb(i - 1) - s_cb(i + 1))*(s_cb(i - 1) - s_cb(i + 2)))
                poly_coef_L(2, 1, i + 1) = &
                    ((s_cb(i - 1) - s_cb(i))*(s_cb(i) - s_cb(i + 1)))/ &
                    ((s_cb(i - 2) - s_cb(i))*(s_cb(i - 2) - s_cb(i + 1)))

                poly_coef_R(0, 1, i + 1) = &
                    ((s_cb(i) - s_cb(i + 2)) + (s_cb(i + 1) - s_cb(i + 3)))/ &
                    ((s_cb(i) - s_cb(i + 2))*(s_cb(i) - s_cb(i + 3)))* &
                    ((s_cb(i) - s_cb(i + 1)))
                poly_coef_R(2, 0, i + 1) = &
                    ((s_cb(i - 2) - s_cb(i + 1)) + (s_cb(i - 1) - s_cb(i + 1)))/ &
                    ((s_cb(i - 1) - s_cb(i + 1))*(s_cb(i + 1) - s_cb(i - 2)))* &
                    ((s_cb(i + 1) - s_cb(i)))
                poly_coef_L(0, 1, i + 1) = &
                    ((s_cb(i) - s_cb(i + 2)) + (s_cb(i) - s_cb(i + 3)))/ &
                    ((s_cb(i) - s_cb(i + 2))*(s_cb(i) - s_cb(i + 3)))* &
                    ((s_cb(i + 1) - s_cb(i)))
                poly_coef_L(2, 0, i + 1) = &
                    ((s_cb(i - 2) - s_cb(i)) + (s_cb(i - 1) - s_cb(i + 1)))/ &
                    ((s_cb(i - 2) - s_cb(i + 1))*(s_cb(i + 1) - s_cb(i - 1)))* &
                    ((s_cb(i) - s_cb(i + 1)))

                d_R(0, i + 1) = &
                    ((s_cb(i - 2) - s_cb(i + 1))*(s_cb(i + 1) - s_cb(i - 1)))/ &
                    ((s_cb(i - 2) - s_cb(i + 3))*(s_cb(i + 3) - s_cb(i - 1)))
                d_R(2, i + 1) = &
                    ((s_cb(i + 1) - s_cb(i + 2))*(s_cb(i + 1) - s_cb(i + 3)))/ &
                    ((s_cb(i - 2) - s_cb(i + 2))*(s_cb(i - 2) - s_cb(i + 3)))
                d_L(0, i + 1) = &
                    ((s_cb(i - 2) - s_cb(i))*(s_cb(i) - s_cb(i - 1)))/ &
                    ((s_cb(i - 2) - s_cb(i + 3))*(s_cb(i + 3) - s_cb(i - 1)))
                d_L(2, i + 1) = &
                    ((s_cb(i) - s_cb(i + 2))*(s_cb(i) - s_cb(i + 3)))/ &
                    ((s_cb(i - 2) - s_cb(i + 2))*(s_cb(i - 2) - s_cb(i + 3)))

                d_R(1, i + 1) = 1d0 - d_R(0, i + 1) - d_R(2, i + 1)
                d_L(1, i + 1) = 1d0 - d_L(0, i + 1) - d_L(2, i + 1)

                beta_coef(0, 0, i + 1) = &
                4d0*(s_cb(i) - s_cb(i + 1))**2d0*(10d0*(s_cb(i + 1) - &
                s_cb(i))**2d0 + (s_cb(i + 1) - s_cb(i))*(s_cb(i + 2) - &
                s_cb(i + 1)) + (s_cb(i + 2) - s_cb(i + 1))**2d0)/((s_cb(i) - &
                s_cb(i + 3))**2d0*(s_cb(i + 1) - s_cb(i + 3))**2d0)

                beta_coef(0, 1, i + 1) = &
                4d0*(s_cb(i) - s_cb(i + 1))**2d0*(19d0*(s_cb(i + 1) - &
                s_cb(i))**2d0 - (s_cb(i + 1) - s_cb(i))*(s_cb(i + 3) - &
                s_cb(i + 1)) + 2d0*(s_cb(i + 2) - s_cb(i))*((s_cb(i + 2) - &
                s_cb(i)) + (s_cb(i + 3) - s_cb(i + 1))))/((s_cb(i) - &
                s_cb(i + 2))*(s_cb(i) - s_cb(i + 3))**2d0*(s_cb(i + 3) - &
                s_cb(i + 1)))

                beta_coef(0, 2, i + 1) = &
                4d0*(s_cb(i) - s_cb(i + 1))**2d0*(10d0*(s_cb(i + 1) - &
                s_cb(i))**2d0 + (s_cb(i + 1) - s_cb(i))*((s_cb(i + 2) - &
                s_cb(i)) + (s_cb(i + 3) - s_cb(i + 1))) + ((s_cb(i + 2) - &
                s_cb(i)) + (s_cb(i + 3) - s_cb(i + 1)))**2d0)/((s_cb(i) - &
                s_cb(i + 2))**2d0*(s_cb(i) - s_cb(i + 3))**2d0)

                beta_coef(1, 0, i + 1) = &
                4d0*(s_cb(i) - s_cb(i + 1))**2d0*(10d0*(s_cb(i + 1) - &
                s_cb(i))**2d0 + (s_cb(i) - s_cb(i - 1))**2d0 + (s_cb(i) - &
                s_cb(i - 1))*(s_cb(i + 1) - s_cb(i)))/((s_cb(i - 1) - &
                s_cb(i + 2))**2d0*(s_cb(i) - s_cb(i + 2))**2d0)

                beta_coef(1, 1, i + 1) = &
                4d0*(s_cb(i) - s_cb(i + 1))**2d0*((s_cb(i) - &
                s_cb(i + 1))*((s_cb(i) - s_cb(i - 1)) + 20d0*(s_cb(i + 1) - &
                s_cb(i))) + (2d0*(s_cb(i) - s_cb(i - 1)) + (s_cb(i + 1) - &
                s_cb(i)))*(s_cb(i + 2) - s_cb(i)))/((s_cb(i + 1) - &
                s_cb(i - 1))*(s_cb(i - 1) - s_cb(i + 2))**2d0*(s_cb(i + 2) - &
                s_cb(i)))

                beta_coef(1, 2, i + 1) = &
                4d0*(s_cb(i) - s_cb(i + 1))**2d0*(10d0*(s_cb(i + 1) - &
                s_cb(i))**2d0 + (s_cb(i + 1) - s_cb(i))*(s_cb(i + 2) - &
                s_cb(i + 1)) + (s_cb(i + 2) - s_cb(i + 1))**2d0)/ &
                ((s_cb(i - 1) - s_cb(i + 1))**2d0*(s_cb(i - 1) - &
                s_cb(i + 2))**2d0)

                beta_coef(2, 0, i + 1) = &
                4d0*(s_cb(i) - s_cb(i + 1))**2d0*(12d0*(s_cb(i + 1) - &
                s_cb(i))**2d0 + ((s_cb(i) - s_cb(i - 2)) + (s_cb(i) - &
                s_cb(i - 1)))**2d0 + 3d0*((s_cb(i) - s_cb(i - 2)) + &
                (s_cb(i) - s_cb(i - 1)))*(s_cb(i + 1) - s_cb(i)))/ &
                ((s_cb(i - 2) - s_cb(i + 1))**2d0*(s_cb(i - 1) - &
                s_cb(i + 1))**2d0)

                beta_coef(2, 1, i + 1) = &
                4d0*(s_cb(i) - s_cb(i + 1))**2d0*(19d0*(s_cb(i + 1) - &
                s_cb(i))**2d0 + ((s_cb(i) - s_cb(i - 2))*(s_cb(i) - &
                s_cb(i + 1))) + 2d0*(s_cb(i + 1) - s_cb(i - 1))*((s_cb(i) - &
                s_cb(i - 2)) + (s_cb(i + 1) - s_cb(i - 1))))/((s_cb(i - 2) - &
                s_cb(i))*(s_cb(i - 2) - s_cb(i + 1))**2d0*(s_cb(i + 1) - &
                s_cb(i - 1)))

                beta_coef(2, 2, i + 1) = &
                4d0*(s_cb(i) - s_cb(i + 1))**2d0*(10d0*(s_cb(i + 1) - &
                s_cb(i))**2d0 + (s_cb(i) - s_cb(i - 1))**2d0 + (s_cb(i) - &
                s_cb(i - 1))*(s_cb(i + 1) - s_cb(i)))/((s_cb(i - 2) - &
                s_cb(i))**2d0*(s_cb(i - 2) - s_cb(i + 1))**2d0)

            end do


            ! Modifying the ideal weights coefficients in the neighborhood
            ! of beginning and end Riemann state extrapolation BC to avoid
            ! any contributions from outside of the physical domain during
            ! the WENO reconstruction
            if (null_weights) then
                if (bc_s%beg == -4) then
                    d_R(1:2, 0) = 0d0; d_R(0, 0) = 1d0
                    d_L(1:2, 0) = 0d0; d_L(0, 0) = 1d0
                    d_R(2, 1) = 0d0; d_R(:, 1) = d_R(:, 1)/sum(d_R(:, 1))
                    d_L(2, 1) = 0d0; d_L(:, 1) = d_L(:, 1)/sum(d_L(:, 1))
                end if

                if (bc_s%end == -4) then
                    d_R(0, s - 1) = 0d0; d_R(:, s - 1) = d_R(:, s - 1)/sum(d_R(:, s - 1))
                    d_L(0, s - 1) = 0d0; d_L(:, s - 1) = d_L(:, s - 1)/sum(d_L(:, s - 1))
                    d_R(0:1, s) = 0d0; d_R(2, s) = 1d0
                    d_L(0:1, s) = 0d0; d_L(2, s) = 1d0
                end if
            end if

        end if
        ! END: Computing WENO5 Coefficients ================================

        ! Nullifying WENO coefficients and cell-boundary locations pointers
        nullify (poly_coef_L, poly_coef_R, d_L, d_R, beta_coef, s_cb)

    end subroutine s_compute_weno_coefficients ! ---------------------------

    !>  The purpose of the procedure is to associate the WENO
        !!      coefficients' pointers with their appropriate targets,
        !!      based on the coordinate direction and the location of
        !!      the WENO reconstruction.
        !! @param weno_dir Coordinate direction of the WENO reconstruction
    subroutine s_associate_weno_coefficients_pointers(weno_dir)

        integer, intent(IN) :: weno_dir

        ! Associating WENO Coefficients in x-direction =====================
        if (weno_dir == 1) then

            poly_coef_L => poly_coef_cbL_x
            poly_coef_R => poly_coef_cbR_x
            d_L => d_cbL_x
            d_R => d_cbR_x
            beta_coef => beta_coef_x

        ! Associating WENO Coefficients in y-direction =====================
        elseif (weno_dir == 2) then

            poly_coef_L => poly_coef_cbL_y
            poly_coef_R => poly_coef_cbR_y
            d_L => d_cbL_y
            d_R => d_cbR_y
            beta_coef => beta_coef_y

        ! Associating WENO Coefficients in z-direction =====================
        else

            poly_coef_L => poly_coef_cbL_z
            poly_coef_R => poly_coef_cbR_z
            d_L => d_cbL_z
            d_R => d_cbR_z
            beta_coef => beta_coef_z

        end if
        ! ==================================================================

    end subroutine s_associate_weno_coefficients_pointers ! ----------------

    !>  WENO reconstruction that is improved with monotonicity
        !!      preserving bounds (MPWENO) and a mapping function that
        !!      boosts the accuracy of the non-linear weights (WENOM).
        !!      MPWENO, Balsara and Shu (2000), prevents reconstructed
        !!      values to reside outside the range set by the stencil,
        !!      while WENOM, Henrick et al. (2005), recovers the order
        !!      of accuracy of the reconstruction for critical points.
        !!      Notice that the basic WENO scheme is implemented based
        !!      on the work of Jiang and Shu (1996).
        !! @param v_vf Cell-averaged variables
        !! @param vL_vf Left WENO reconstructed cell-boundary values
        !! @param vR_vf Right WENO reconstructed cell-boundary values
        !! @param norm_dir Characteristic decommposition coordinate direction
        !! @param weno_dir Coordinate direction of the WENO reconstruction
        !! @param ix Index bounds in first coordinate direction
        !! @param iy Index bounds in second coordinate direction
        !! @param iz Index bounds in third coordinate direction
    subroutine s_weno(v_vf, vL_vf, vR_vf, & ! -------------------
                      norm_dir, weno_dir,  &
                      ix, iy, iz)

        type(scalar_field), dimension(:), intent(IN) :: v_vf
        type(scalar_field), dimension(:), intent(INOUT) :: vL_vf, vR_vf
        integer, intent(IN) :: norm_dir
        integer, intent(IN) :: weno_dir
        type(bounds_info), intent(IN) :: ix, iy, iz

        real(kind(0d0)), dimension(-weno_polyn:weno_polyn - 1) :: dvd !<
            !! Newton divided differences

        real(kind(0d0)), dimension(0:weno_polyn) ::  poly_L, poly_R  !< Left/right polynominals
        real(kind(0d0)), dimension(0:weno_polyn) :: alpha_L, alpha_R  !< Left/right nonlinear weights
        real(kind(0d0)), dimension(0:weno_polyn) :: omega_L, omega_R  !< Left/right nonlinear weights

        real(kind(0d0)), dimension(0:weno_polyn) :: beta !< Smoothness indicators

        real(kind(0d0)), dimension(-weno_polyn:weno_polyn) :: scaling_stencil, scaled_vars
        real(kind(0d0)) :: min_u, max_u

        integer :: i, j, k, l, q !< Generic loop iterators

        ! Reshaping and/or projecting onto characteristic fields inputted
        ! data and in addition associating the WENO coefficients pointers
        if (weno_order /= 1) then
            call s_initialize_weno(v_vf, vL_vf, vR_vf, &
                                   norm_dir, weno_dir, ix, iy, iz)
            call s_associate_weno_coefficients_pointers(weno_dir)
        end if

        ! WENO1 ============================================================
        if (weno_order == 1) then

            do i = 1, ubound(v_vf, 1)
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end
                            vL_vf(i)%sf(j, k, l) = v_vf(i)%sf(j, k, l)
                            vR_vf(i)%sf(j, k, l) = v_vf(i)%sf(j, k, l)
                        end do
                    end do
                end do
            end do

            ! ==================================================================

            ! WENO3 ============================================================
        elseif (weno_order == 3) then

            do i = 1, v_size
                do l = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = is1%beg, is1%end
                            ! reconstruct from left side

                            dvd(0) = v_rs_wsL(1)%vf(i)%sf(j, k, l) &
                                     - v_rs_wsL(0)%vf(i)%sf(j, k, l)
                            dvd(-1) = v_rs_wsL(0)%vf(i)%sf(j, k, l) &
                                      - v_rs_wsL(-1)%vf(i)%sf(j, k, l)

                            poly_L(0) = v_rs_wsL(0)%vf(i)%sf(j, k, l) &
                                        + poly_coef_L(0, 0, j)*dvd(0)
                            poly_L(1) = v_rs_wsL(0)%vf(i)%sf(j, k, l) &
                                        + poly_coef_L(1, 0, j)*dvd(-1)

                            beta(0) = beta_coef(0, 0, j)*dvd(0)*dvd(0) &
                                      + weno_eps
                            beta(1) = beta_coef(1, 0, j)*dvd(-1)*dvd(-1) &
                                      + weno_eps

                            alpha_L = d_L(:, j)/(beta*beta)

                            omega_L = alpha_L/sum(alpha_L)

                            ! reconstruct from right side
                            dvd(0) = v_rs_wsR(1)%vf(i)%sf(j, k, l) &
                                     - v_rs_wsR(0)%vf(i)%sf(j, k, l)
                            dvd(-1) = v_rs_wsR(0)%vf(i)%sf(j, k, l) &
                                      - v_rs_wsR(-1)%vf(i)%sf(j, k, l)

                            poly_R(0) = v_rs_wsR(0)%vf(i)%sf(j, k, l) &
                                        + poly_coef_R(0, 0, j)*dvd(0)
                            poly_R(1) = v_rs_wsR(0)%vf(i)%sf(j, k, l) &
                                        + poly_coef_R(1, 0, j)*dvd(-1)

                            beta(0) = beta_coef(0, 0, j)*dvd(0)*dvd(0) &
                                      + weno_eps
                            beta(1) = beta_coef(1, 0, j)*dvd(-1)*dvd(-1) &
                                      + weno_eps

                            alpha_R = d_R(:, j)/(beta*beta)

                            omega_R = alpha_R/sum(alpha_R)

                            if (mapped_weno) then
                                call s_map_nonlinear_weights(d_L(:, j), &
                                                             alpha_L, &
                                                             omega_L)
                                call s_map_nonlinear_weights(d_R(:, j), &
                                                             alpha_R, &
                                                             omega_R)
                            end if

                            vL_rs_vf(i)%sf(j, k, l) = sum(omega_L*poly_L)
                            vR_rs_vf(i)%sf(j, k, l) = sum(omega_R*poly_R)

                        end do
                    end do
                end do
            end do

            ! END: WENO3 =======================================================

            ! WENO5 ============================================================
        else

            do i = 1, v_size
                do l = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = is1%beg, is1%end

                            dvd(1) = v_rs_wsL(2)%vf(i)%sf(j, k, l) &
                                     - v_rs_wsL(1)%vf(i)%sf(j, k, l)
                            dvd(0) = v_rs_wsL(1)%vf(i)%sf(j, k, l) &
                                     - v_rs_wsL(0)%vf(i)%sf(j, k, l)
                            dvd(-1) = v_rs_wsL(0)%vf(i)%sf(j, k, l) &
                                      - v_rs_wsL(-1)%vf(i)%sf(j, k, l)
                            dvd(-2) = v_rs_wsL(-1)%vf(i)%sf(j, k, l) &
                                      - v_rs_wsL(-2)%vf(i)%sf(j, k, l)

                            poly_L(0) = v_rs_wsL(0)%vf(i)%sf(j, k, l) &
                                        + poly_coef_L(0, 0, j)*dvd(1) &
                                        + poly_coef_L(0, 1, j)*dvd(0)
                            poly_L(1) = v_rs_wsL(0)%vf(i)%sf(j, k, l) &
                                        + poly_coef_L(1, 0, j)*dvd(0) &
                                        + poly_coef_L(1, 1, j)*dvd(-1)
                            poly_L(2) = v_rs_wsL(0)%vf(i)%sf(j, k, l) &
                                        + poly_coef_L(2, 0, j)*dvd(-1) &
                                        + poly_coef_L(2, 1, j)*dvd(-2)

                            beta(0) = beta_coef(0, 0, j)*dvd(1)*dvd(1) &
                                      + beta_coef(0, 1, j)*dvd(1)*dvd(0) &
                                      + beta_coef(0, 2, j)*dvd(0)*dvd(0) &
                                      + weno_eps
                            beta(1) = beta_coef(1, 0, j)*dvd(0)*dvd(0) &
                                      + beta_coef(1, 1, j)*dvd(0)*dvd(-1) &
                                      + beta_coef(1, 2, j)*dvd(-1)*dvd(-1) &
                                      + weno_eps
                            beta(2) = beta_coef(2, 0, j)*dvd(-1)*dvd(-1) &
                                      + beta_coef(2, 1, j)*dvd(-1)*dvd(-2) &
                                      + beta_coef(2, 2, j)*dvd(-2)*dvd(-2) &
                                      + weno_eps

                            alpha_L = d_L(:, j)/(beta*beta)

                            omega_L = alpha_L/sum(alpha_L)

                            dvd(1) = v_rs_wsR(2)%vf(i)%sf(j, k, l) &
                                     - v_rs_wsR(1)%vf(i)%sf(j, k, l)
                            dvd(0) = v_rs_wsR(1)%vf(i)%sf(j, k, l) &
                                     - v_rs_wsR(0)%vf(i)%sf(j, k, l)
                            dvd(-1) = v_rs_wsR(0)%vf(i)%sf(j, k, l) &
                                      - v_rs_wsR(-1)%vf(i)%sf(j, k, l)
                            dvd(-2) = v_rs_wsR(-1)%vf(i)%sf(j, k, l) &
                                      - v_rs_wsR(-2)%vf(i)%sf(j, k, l)

                            poly_R(0) = v_rs_wsR(0)%vf(i)%sf(j, k, l) &
                                        + poly_coef_R(0, 0, j)*dvd(1) &
                                        + poly_coef_R(0, 1, j)*dvd(0)
                            poly_R(1) = v_rs_wsR(0)%vf(i)%sf(j, k, l) &
                                        + poly_coef_R(1, 0, j)*dvd(0) &
                                        + poly_coef_R(1, 1, j)*dvd(-1)
                            poly_R(2) = v_rs_wsR(0)%vf(i)%sf(j, k, l) &
                                        + poly_coef_R(2, 0, j)*dvd(-1) &
                                        + poly_coef_R(2, 1, j)*dvd(-2)

                            beta(0) = beta_coef(0, 0, j)*dvd(1)*dvd(1) &
                                      + beta_coef(0, 1, j)*dvd(1)*dvd(0) &
                                      + beta_coef(0, 2, j)*dvd(0)*dvd(0) &
                                      + weno_eps
                            beta(1) = beta_coef(1, 0, j)*dvd(0)*dvd(0) &
                                      + beta_coef(1, 1, j)*dvd(0)*dvd(-1) &
                                      + beta_coef(1, 2, j)*dvd(-1)*dvd(-1) &
                                      + weno_eps
                            beta(2) = beta_coef(2, 0, j)*dvd(-1)*dvd(-1) &
                                      + beta_coef(2, 1, j)*dvd(-1)*dvd(-2) &
                                      + beta_coef(2, 2, j)*dvd(-2)*dvd(-2) &
                                      + weno_eps

                            alpha_R = d_R(:, j)/(beta*beta)

                            omega_R = alpha_R/sum(alpha_R)

                            if (mapped_weno) then
                                call s_map_nonlinear_weights(d_L(:, j), &
                                                             alpha_L, &
                                                             omega_L)
                                call s_map_nonlinear_weights(d_R(:, j), &
                                                             alpha_R, &
                                                             omega_R)
                            end if

                            vL_rs_vf(i)%sf(j, k, l) = sum(omega_L*poly_L)
                            vR_rs_vf(i)%sf(j, k, l) = sum(omega_R*poly_R)

                            if (mp_weno) then
                                call s_preserve_monotonicity(i, j, k, l)
                            end if

                        end do
                    end do
                end do
            end do

        end if
        ! END: WENO5 =======================================================

        ! Reshaping and/or projecting onto physical fields the outputted
        ! data, as well as disassociating the WENO coefficients pointers
        if (weno_order /= 1) then
            call s_finalize_weno(vL_vf, vR_vf, weno_dir, ix, iy, iz)
            nullify (poly_coef_L, poly_coef_R, d_L, d_R, beta_coef)
        end if

    end subroutine s_weno ! ------------------------------------------------

    !> The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are required for the setup of the
        !!      WENO reconstruction.
        !! @param v_vf Cell-averaged variables
        !! @param vL_vf Left WENO reconstructed cell-boundary values
        !! @param vR_vf Right WENO reconstructed cell-boundary values
        !! @param norm_dir Characteristic decommposition coordinate direction
        !! @param weno_dir Coordinate direction of the WENO reconstruction
        !! @param ix Index bounds in first coordinate direction
        !! @param iy Index bounds in second coordinate direction
        !! @param iz Index bounds in third coordinate direction
    subroutine s_initialize_weno(v_vf, vL_vf, vR_vf, & ! ---------
                                 norm_dir, weno_dir, ix, iy, iz)

        type(scalar_field), dimension(:), intent(IN) :: v_vf
        type(scalar_field), dimension(:), intent(INOUT) :: vL_vf, vR_vf
        integer, intent(IN) :: norm_dir
        integer, intent(IN) :: weno_dir
        type(bounds_info), intent(IN) :: ix, iy, iz

        integer :: i, j, k, l !< Generic loop iterators

        ! Determining the number of cell-average variables which will be
        ! WENO-reconstructed and mapping their indical bounds in the x-,
        ! y- and z-directions to those in the s1-, s2- and s3-directions
        ! as to reshape the inputted data in the coordinate direction of
        ! the WENO reconstruction
        v_size = ubound(v_vf, 1)

        if (weno_dir == 1) then
            is1 = ix; is2 = iy; is3 = iz
        elseif (weno_dir == 2) then
            is1 = iy; is2 = ix; is3 = iz
        else
            is1 = iz; is2 = iy; is3 = ix
        end if

        ! Allocating the cell-average variables, which are reshaped, and/or
        ! characteristically decomposed, in the coordinate direction of the
        ! WENO reconstruction
        do i = -weno_polyn, weno_polyn

            allocate (v_rs_wsL(i)%vf(1:v_size), v_rs_wsR(i)%vf(1:v_size))

            do j = 1, v_size

                allocate (v_rs_wsL(i)%vf(j)%sf(is1%beg:is1%end, &
                                               is2%beg:is2%end, &
                                               is3%beg:is3%end))

                v_rs_wsR(i)%vf(j)%sf => v_rs_wsL(i)%vf(j)%sf
            end do

        end do

        ! Allocating the left and right WENO reconstructions of the cell-
        ! average variables that are reshaped, and/or characteristically
        ! decomposed, in the coordinate direction of WENO reconstruction
        allocate (vL_rs_vf(1:v_size), vR_rs_vf(1:v_size))

        if (weno_dir == 1) then
            do i = 1, v_size
                vL_rs_vf(i)%sf => vL_vf(i)%sf
                vR_rs_vf(i)%sf => vR_vf(i)%sf
            end do
        else
            do i = 1, v_size
                allocate (vL_rs_vf(i)%sf(is1%beg:is1%end, &
                                         is2%beg:is2%end, &
                                         is3%beg:is3%end))
                allocate (vR_rs_vf(i)%sf(is1%beg:is1%end, &
                                         is2%beg:is2%end, &
                                         is3%beg:is3%end))
            end do
        end if

        ! Reshaping/Projecting onto Characteristic Fields in x-direction ===
        if (weno_dir == 1) then

            do i = -weno_polyn, weno_polyn
                do j = 1, v_size
                    do k = ix%beg, ix%end
                        v_rs_wsL(i)%vf(j)%sf(k, :, :) = &
                            v_vf(j)%sf(i + k, iy%beg:iy%end, iz%beg:iz%end)
                    end do
                end do
            end do

            ! ==================================================================

            ! Reshaping/Projecting onto Characteristic Fields in y-direction ===
        elseif (weno_dir == 2) then

            do i = -weno_polyn, weno_polyn
                do j = 1, v_size
                    do k = ix%beg, ix%end
                        do l = iy%beg, iy%end
                            v_rs_wsL(i)%vf(j)%sf(l, k, :) = &
                                v_vf(j)%sf(k, i + l, iz%beg:iz%end)
                        end do
                    end do
                end do
            end do


            ! ==================================================================

            ! Reshaping/Projecting onto Characteristic Fields in z-direction ===
        else

            do i = -weno_polyn, weno_polyn
                do j = 1, v_size
                    do k = ix%beg, ix%end
                        do l = iz%beg, iz%end
                            v_rs_wsL(i)%vf(j)%sf(l, :, k) = &
                                v_vf(j)%sf(k, iy%beg:iy%end, i + l)
                        end do
                    end do
                end do
            end do


        end if
        ! ==================================================================

    end subroutine s_initialize_weno ! -------------------------------------


    !>  The goal of this procedure is to map the nonlinear WENO
        !!      weights to the more accurate nonlinear WENOM weights in
        !!      order to reinstate the optimal order of accuracy of the
        !!      reconstruction in the proximity of critical points, see
        !!      Henrick et al. (2005).
        !!  @param d_K Cell boundary pointer
        !!  @param alpha_K ideal weights
        !!  @param omega_K nonlinear weights
    subroutine s_map_nonlinear_weights(d_K, alpha_K, omega_K) ! ------------

        ! Ideal and nonlinear weights
        real(kind(0d0)), dimension(0:weno_polyn), intent(IN)    ::     d_K
        real(kind(0d0)), dimension(0:weno_polyn), intent(INOUT) :: alpha_K
        real(kind(0d0)), dimension(0:weno_polyn), intent(INOUT) :: omega_K

        ! Mapping the WENO nonlinear weights to the WENOM nonlinear weights
        if (minval(d_K) == 0d0 .or. maxval(d_K) == 1d0) return

        alpha_K = (d_K*(1d0 + d_K - 3d0*omega_K) + omega_K**2d0) &
                  *(omega_K/(d_K**2d0 + omega_K*(1d0 - 2d0*d_K)))

        omega_K = alpha_K/sum(alpha_K)

    end subroutine s_map_nonlinear_weights ! -------------------------------

    !>  The goal of this subroutine is to ensure that the WENO
        !!      reconstruction is monotonic. The latter is achieved by
        !!      enforcing monotonicity preserving bounds of Suresh and
        !!      Huynh (1997). The resulting MPWENO reconstruction, see
        !!      Balsara and Shu (2000), ensures that the reconstructed
        !!      values do not reside outside the range spanned by WENO
        !!      stencil.
        !!  @param i Equation number
        !!  @param j First-coordinate cell index
        !!  @param k Second-coordinate cell index
        !!  @param l Third-coordinate cell index
    subroutine s_preserve_monotonicity(i, j, k, l) ! --------------------------

        integer, intent(IN) :: i, j, k, l

        real(kind(0d0)), dimension(-1:1) :: d !< Curvature measures at the zone centers

        real(kind(0d0)) :: d_MD, d_LC !<
            !! Median (md) curvature and large curvature (LC) measures

        ! The left and right upper bounds (UL), medians, large curvatures,
        ! minima, and maxima of the WENO-reconstructed values of the cell-
        ! average variables.
        real(kind(0d0)) :: vL_UL, vR_UL
        real(kind(0d0)) :: vL_MD, vR_MD
        real(kind(0d0)) :: vL_LC, vR_LC
        real(kind(0d0)) :: vL_min, vR_min
        real(kind(0d0)) :: vL_max, vR_max

        real(kind(0d0)), parameter :: alpha = 2d0 !>
            !! Determines the maximum Courant–Friedrichs–Lewy (CFL) number that
            !! may be utilized with the scheme. In theory, for stability, a CFL
            !! number less than 1/(1+alpha) is necessary. The default value for
            !! alpha is 2.

        real(kind(0d0)), parameter :: beta = 4d0/3d0 !<
            !! Determines the amount of freedom available from utilizing a large
            !! value for the local curvature. The default value for beta is 4/3.

        ! Left Monotonicity Preserving Bound ===============================
        d(-1) = v_rs_wsL(0)%vf(i)%sf(j, k, l) &
                + v_rs_wsL(-2)%vf(i)%sf(j, k, l) &
                - v_rs_wsL(-1)%vf(i)%sf(j, k, l) &
                *2d0
        d(0) = v_rs_wsL(1)%vf(i)%sf(j, k, l) &
               + v_rs_wsL(-1)%vf(i)%sf(j, k, l) &
               - v_rs_wsL(0)%vf(i)%sf(j, k, l) &
               *2d0
        d(1) = v_rs_wsL(2)%vf(i)%sf(j, k, l) &
               + v_rs_wsL(0)%vf(i)%sf(j, k, l) &
               - v_rs_wsL(1)%vf(i)%sf(j, k, l) &
               *2d0

        d_MD = (sign(1d0, 4d0*d(-1) - d(0)) + sign(1d0, 4d0*d(0) - d(-1))) &
               *abs((sign(1d0, 4d0*d(-1) - d(0)) + sign(1d0, d(-1))) &
                    *(sign(1d0, 4d0*d(-1) - d(0)) + sign(1d0, d(0)))) &
               *min(abs(4d0*d(-1) - d(0)), abs(d(-1)), &
                    abs(4d0*d(0) - d(-1)), abs(d(0)))/8d0

        d_LC = (sign(1d0, 4d0*d(0) - d(1)) + sign(1d0, 4d0*d(1) - d(0))) &
               *abs((sign(1d0, 4d0*d(0) - d(1)) + sign(1d0, d(0))) &
                    *(sign(1d0, 4d0*d(0) - d(1)) + sign(1d0, d(1)))) &
               *min(abs(4d0*d(0) - d(1)), abs(d(0)), &
                    abs(4d0*d(1) - d(0)), abs(d(1)))/8d0

        vL_UL = v_rs_wsL(0)%vf(i)%sf(j, k, l) &
                - (v_rs_wsL(1)%vf(i)%sf(j, k, l) &
                   - v_rs_wsL(0)%vf(i)%sf(j, k, l))*alpha

        vL_MD = (v_rs_wsL(0)%vf(i)%sf(j, k, l) &
                 + v_rs_wsL(-1)%vf(i)%sf(j, k, l) &
                 - d_MD)*5d-1

        vL_LC = v_rs_wsL(0)%vf(i)%sf(j, k, l) &
                - (v_rs_wsL(1)%vf(i)%sf(j, k, l) &
                   - v_rs_wsL(0)%vf(i)%sf(j, k, l))*5d-1 + beta*d_LC

        vL_min = max(min(v_rs_wsL(0)%vf(i)%sf(j, k, l), &
                         v_rs_wsL(-1)%vf(i)%sf(j, k, l), &
                         vL_MD), &
                     min(v_rs_wsL(0)%vf(i)%sf(j, k, l), &
                         vL_UL, &
                         vL_LC))

        vL_max = min(max(v_rs_wsL(0)%vf(i)%sf(j, k, l), &
                         v_rs_wsL(-1)%vf(i)%sf(j, k, l), &
                         vL_MD), &
                     max(v_rs_wsL(0)%vf(i)%sf(j, k, l), &
                         vL_UL, &
                         vL_LC))

        vL_rs_vf(i)%sf(j, k, l) = vL_rs_vf(i)%sf(j, k, l) &
                                  + (sign(5d-1, vL_min - vL_rs_vf(i)%sf(j, k, l)) &
                                     + sign(5d-1, vL_max - vL_rs_vf(i)%sf(j, k, l))) &
                                  *min(abs(vL_min - vL_rs_vf(i)%sf(j, k, l)), &
                                       abs(vL_max - vL_rs_vf(i)%sf(j, k, l)))
        ! END: Left Monotonicity Preserving Bound ==========================

        ! Right Monotonicity Preserving Bound ==============================
        d(-1) = v_rs_wsR(0)%vf(i)%sf(j, k, l) &
                + v_rs_wsR(-2)%vf(i)%sf(j, k, l) &
                - v_rs_wsR(-1)%vf(i)%sf(j, k, l)*2d0
        d(0) = v_rs_wsR(1)%vf(i)%sf(j, k, l) &
               + v_rs_wsR(-1)%vf(i)%sf(j, k, l) &
               - v_rs_wsR(0)%vf(i)%sf(j, k, l)*2d0
        d(1) = v_rs_wsR(2)%vf(i)%sf(j, k, l) &
               + v_rs_wsR(0)%vf(i)%sf(j, k, l) &
               - v_rs_wsR(1)%vf(i)%sf(j, k, l)*2d0

        d_MD = (sign(1d0, 4d0*d(0) - d(1)) + sign(1d0, 4d0*d(1) - d(0))) &
               *abs((sign(1d0, 4d0*d(0) - d(1)) + sign(1d0, d(0))) &
                    *(sign(1d0, 4d0*d(0) - d(1)) + sign(1d0, d(1)))) &
               *min(abs(4d0*d(0) - d(1)), abs(d(0)), &
                    abs(4d0*d(1) - d(0)), abs(d(1)))/8d0

        d_LC = (sign(1d0, 4d0*d(-1) - d(0)) + sign(1d0, 4d0*d(0) - d(-1))) &
               *abs((sign(1d0, 4d0*d(-1) - d(0)) + sign(1d0, d(-1))) &
                    *(sign(1d0, 4d0*d(-1) - d(0)) + sign(1d0, d(0)))) &
               *min(abs(4d0*d(-1) - d(0)), abs(d(-1)), &
                    abs(4d0*d(0) - d(-1)), abs(d(0)))/8d0

        vR_UL = v_rs_wsR(0)%vf(i)%sf(j, k, l) &
                + (v_rs_wsR(0)%vf(i)%sf(j, k, l) &
                   - v_rs_wsR(-1)%vf(i)%sf(j, k, l))*alpha

        vR_MD = (v_rs_wsR(0)%vf(i)%sf(j, k, l) &
                 + v_rs_wsR(1)%vf(i)%sf(j, k, l) &
                 - d_MD)*5d-1

        vR_LC = v_rs_wsR(0)%vf(i)%sf(j, k, l) &
                + (v_rs_wsR(0)%vf(i)%sf(j, k, l) &
                   - v_rs_wsR(-1)%vf(i)%sf(j, k, l))*5d-1 + beta*d_LC

        vR_min = max(min(v_rs_wsR(0)%vf(i)%sf(j, k, l), &
                         v_rs_wsR(1)%vf(i)%sf(j, k, l), &
                         vR_MD), &
                     min(v_rs_wsR(0)%vf(i)%sf(j, k, l), &
                         vR_UL, &
                         vR_LC))

        vR_max = min(max(v_rs_wsR(0)%vf(i)%sf(j, k, l), &
                         v_rs_wsR(1)%vf(i)%sf(j, k, l), &
                         vR_MD), &
                     max(v_rs_wsR(0)%vf(i)%sf(j, k, l), &
                         vR_UL, &
                         vR_LC))

        vR_rs_vf(i)%sf(j, k, l) = vR_rs_vf(i)%sf(j, k, l) &
                                  + (sign(5d-1, vR_min - vR_rs_vf(i)%sf(j, k, l)) &
                                     + sign(5d-1, vR_max - vR_rs_vf(i)%sf(j, k, l))) &
                                  *min(abs(vR_min - vR_rs_vf(i)%sf(j, k, l)), &
                                       abs(vR_max - vR_rs_vf(i)%sf(j, k, l)))
        ! END: Right Monotonicity Preserving Bound =========================

    end subroutine s_preserve_monotonicity ! -------------------------------

    !>  Deallocation and/or disassociation procedures that are
        !!      necessary in order to finalize the WENO reconstruction
        !! @param vL_vf Left WENO reconstructed cell-boundary values
        !! @param vR_vf Right WENO reconstructed cell-boundary values
        !! @param weno_dir Coordinate direction of the WENO reconstruction
        !! @param ix Index bounds in first coordinate direction
        !! @param iy Index bounds in second coordinate direction
        !! @param iz Index bounds in third coordinate direction
    subroutine s_finalize_weno(vL_vf, vR_vf, & ! -----------------
                               weno_dir, ix, iy, iz)

        type(scalar_field), dimension(:), intent(INOUT) :: vL_vf, vR_vf
        integer, intent(IN) :: weno_dir
        type(bounds_info), intent(IN) :: ix, iy, iz

        integer :: i, j, k !< Generic loop iterators

        ! Reshaping/Projecting onto Physical Fields in x-direction =========
        if (weno_dir == 1) then

            ! ==================================================================

            ! Reshaping/Projecting onto Physical Fields in y-direction =========
        elseif (weno_dir == 2) then


            do i = 1, v_size
                do j = ix%beg, ix%end
                    do k = iy%beg, iy%end
                        vL_vf(i)%sf(j, k, iz%beg:iz%end) = vL_rs_vf(i)%sf(k, j, :)
                        vR_vf(i)%sf(j, k, iz%beg:iz%end) = vR_rs_vf(i)%sf(k, j, :)
                    end do
                end do
            end do

            ! ==================================================================

            ! Reshaping/Projecting onto Physical Fields in z-direction =========
        else


            do i = 1, v_size
                do j = ix%beg, ix%end
                    do k = iz%beg, iz%end
                        vL_vf(i)%sf(j, iy%beg:iy%end, k) = vL_rs_vf(i)%sf(k, :, j)
                        vR_vf(i)%sf(j, iy%beg:iy%end, k) = vR_rs_vf(i)%sf(k, :, j)
                    end do
                end do
            end do

        end if
        ! ==================================================================

        ! Deallocating the cell-average variables that were reshaped and/or
        ! characteristically decomposed in the coordinate direction of WENO
        ! reconstruction
        do i = -weno_polyn, weno_polyn
            do j = 1, v_size
                deallocate (v_rs_wsL(i)%vf(j)%sf)
                v_rs_wsR(i)%vf(j)%sf => null()
            end do
            deallocate (v_rs_wsL(i)%vf, v_rs_wsR(i)%vf)
        end do

        ! Deallocating the left and right WENO reconstructions of the cell-
        ! average variables which were reshaped, and/or characteristically
        ! decomposed, in the coordinate direction of WENO reconstruction
        if (weno_dir == 1) then
            do i = 1, v_size
                vL_rs_vf(i)%sf => null()
                vR_rs_vf(i)%sf => null()
            end do
        else
            do i = 1, v_size
                deallocate (vL_rs_vf(i)%sf)
                deallocate (vR_rs_vf(i)%sf)
            end do
        end if

        deallocate (vL_rs_vf, vR_rs_vf)

    end subroutine s_finalize_weno ! ---------------------------------------

    !>  Module deallocation and/or disassociation procedures
    subroutine s_finalize_weno_module() ! ----------------------------------

        if (weno_order == 1) return

        ! Deallocating the WENO-stencil of the WENO-reconstructed variables
        deallocate (v_rs_wsL, v_rs_wsR)

        ! Deallocating WENO coefficients in x-direction ====================
        deallocate (poly_coef_cbL_x, poly_coef_cbR_x)
        deallocate (d_cbL_x, d_cbR_x)

        deallocate (beta_coef_x)
        ! ==================================================================

        ! Deallocating WENO coefficients in y-direction ====================
        if (n == 0) return

        deallocate (poly_coef_cbL_y, poly_coef_cbR_y)
        deallocate (d_cbL_y, d_cbR_y)

        deallocate (beta_coef_y)
        ! ==================================================================

        ! Deallocating WENO coefficients in z-direction ====================
        if (p == 0) return

        deallocate (poly_coef_cbL_z, poly_coef_cbR_z)
        deallocate (d_cbL_z, d_cbR_z)


        deallocate (beta_coef_z)
        ! ==================================================================

    end subroutine s_finalize_weno_module ! --------------------------------

end module m_weno
