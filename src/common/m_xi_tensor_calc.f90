!>
!! @file m_xi_tensor_calc.f90
!! @brief Contains module m_xi_tensor_calc

!> @brief This module consists of subroutines used in the calculation of matrix
!!              operations for the reference map tensor

module m_xi_tensor_calc

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters
    ! ==========================================================================

    implicit none

    private; public :: s_compute_gradient_xi, &
 s_compute_gradient_xi1d_acc, &
 s_compute_gradient_xi2d_acc, &
 s_compute_gradient_xi3d_acc, &
 f_elastic_energy

contains

    !>  The following subroutine handles the calculation of the btensor.
        !!   The calculation of the btensor takes qprimvf.
        !! @param q_prim_vf Primitive variables
        !! @param btensor is the output
        !! calculate the grad_xi, grad_xi is a nxn tensor
        !! calculate the inverse of grad_xi to obtain F, F is a nxn tensor
        !! calculate the FFtranspose to obtain the btensor, btensor is nxn tensor
        !! btensor is symmetric, save the data space
    subroutine s_compute_gradient_xi(q_prim_vf, xb, xe, yb, ye, & !---------
                                     zb, ze, j, k, l, tensora, tensorb)

        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf
        real(kind(0d0)), dimension(tensor_size), intent(INOUT) :: tensora, tensorb
        integer, intent(IN) :: xb, xe, yb, ye, zb, ze
        integer, intent(IN) :: j, k, l

        real(kind(0d0)) :: determinant
        integer :: i
        ! STEP 1: computing the grad_xi tensor
        ! grad_xi definition / organization
        ! number for the tensor 1-3:  dxix_dx, dxiy_dx, dxiz_dx
        ! 4-6 :                       dxix_dy, dxiy_dy, dxiz_dy
        ! 7-9 :                       dxix_dz, dxiy_dz, dxiz_dz
        if (j == xb) then
            ! dxix/dx
            tensora(1) = (-25d0*q_prim_vf(xibeg)%sf(j, k, l) &
                          + 48d0*q_prim_vf(xibeg)%sf(j + 1, k, l) &
                          - 36d0*q_prim_vf(xibeg)%sf(j + 2, k, l) &
                          + 16d0*q_prim_vf(xibeg)%sf(j + 3, k, l) &
                          - 3d0*q_prim_vf(xibeg)%sf(j + 4, k, l)) &
                         /(12d0*(x_cb(j + 1) - x_cb(j)))
        else if (j == xb + 1) then
            ! dxix/dx
            tensora(1) = (-3d0*q_prim_vf(xibeg)%sf(j - 1, k, l) &
                          - 10d0*q_prim_vf(xibeg)%sf(j, k, l) &
                          + 18d0*q_prim_vf(xibeg)%sf(j + 1, k, l) &
                          - 6d0*q_prim_vf(xibeg)%sf(j + 2, k, l) &
                          + q_prim_vf(xibeg)%sf(j + 3, k, l)) &
                         /(12d0*(x_cb(j) - x_cb(j - 1)))
        else if (j == xe - 1) then
            ! dxix/dx
            tensora(1) = (3d0*q_prim_vf(xibeg)%sf(j + 1, k, l) &
                          + 10d0*q_prim_vf(xibeg)%sf(j, k, l) &
                          - 18d0*q_prim_vf(xibeg)%sf(j - 1, k, l) &
                          + 6d0*q_prim_vf(xibeg)%sf(j - 2, k, l) &
                          - q_prim_vf(xibeg)%sf(j - 3, k, l)) &
                         /(12d0*(x_cb(j) - x_cb(j - 1)))
        else if (j == xe) then
            ! dxix/dx
            tensora(1) = (25d0*q_prim_vf(xibeg)%sf(j, k, l) &
                          - 48d0*q_prim_vf(xibeg)%sf(j - 1, k, l) &
                          + 36d0*q_prim_vf(xibeg)%sf(j - 2, k, l) &
                          - 16d0*q_prim_vf(xibeg)%sf(j - 3, k, l) &
                          + 3d0*q_prim_vf(xibeg)%sf(j - 4, k, l)) &
                         /(12d0*(x_cb(j) - x_cb(j - 1)))
        else
            ! dxix/dx
            tensora(1) = (q_prim_vf(xibeg)%sf(j - 2, k, l) &
                          - 8d0*q_prim_vf(xibeg)%sf(j - 1, k, l) &
                          + 8d0*q_prim_vf(xibeg)%sf(j + 1, k, l) &
                          - q_prim_vf(xibeg)%sf(j + 2, k, l)) &
                         /(12d0*(x_cb(j) - x_cb(j - 1)))
        end if

        if (num_dims > 1) then
            if (j == xb) then
                ! dxiy / dx
                tensora(2) = (-25d0*q_prim_vf(xibeg + 1)%sf(j, k, l) &
                              + 48d0*q_prim_vf(xibeg + 1)%sf(j + 1, k, l) &
                              - 36d0*q_prim_vf(xibeg + 1)%sf(j + 2, k, l) &
                              + 16d0*q_prim_vf(xibeg + 1)%sf(j + 3, k, l) &
                              - 3d0*q_prim_vf(xibeg + 1)%sf(j + 4, k, l)) &
                             /(12d0*(x_cb(j + 1) - x_cb(j)))
            else if (j == xb + 1) then
                ! dxiy / dx
                tensora(2) = (-3d0*q_prim_vf(xibeg + 1)%sf(j - 1, k, l) &
                              - 10d0*q_prim_vf(xibeg + 1)%sf(j, k, l) &
                              + 18d0*q_prim_vf(xibeg + 1)%sf(j + 1, k, l) &
                              - 6d0*q_prim_vf(xibeg + 1)%sf(j + 2, k, l) &
                              + q_prim_vf(xibeg + 1)%sf(j + 3, k, l)) &
                             /(12d0*(x_cb(j) - x_cb(j - 1)))
            else if (j == xe - 1) then
                ! dxiy / dx
                tensora(2) = (3d0*q_prim_vf(xibeg + 1)%sf(j + 1, k, l) &
                              + 10d0*q_prim_vf(xibeg + 1)%sf(j, k, l) &
                              - 18d0*q_prim_vf(xibeg + 1)%sf(j - 1, k, l) &
                              + 6d0*q_prim_vf(xibeg + 1)%sf(j - 2, k, l) &
                              - q_prim_vf(xibeg + 1)%sf(j - 3, k, l)) &
                             /(12d0*(x_cb(j) - x_cb(j - 1)))
            else if (j == xe) then
                ! dxiy / dx
                tensora(2) = (25d0*q_prim_vf(xibeg + 1)%sf(j, k, l) &
                              - 48d0*q_prim_vf(xibeg + 1)%sf(j - 1, k, l) &
                              + 36d0*q_prim_vf(xibeg + 1)%sf(j - 2, k, l) &
                              - 16d0*q_prim_vf(xibeg + 1)%sf(j - 3, k, l) &
                              + 3d0*q_prim_vf(xibeg + 1)%sf(j - 4, k, l)) &
                             /(12d0*(x_cb(j) - x_cb(j - 1)))
            else
                ! dxiy / dx
                tensora(2) = (q_prim_vf(xibeg + 1)%sf(j - 2, k, l) &
                              - 8d0*q_prim_vf(xibeg + 1)%sf(j - 1, k, l) &
                              + 8d0*q_prim_vf(xibeg + 1)%sf(j + 1, k, l) &
                              - q_prim_vf(xibeg + 1)%sf(j + 2, k, l)) &
                             /(12d0*(x_cb(j) - x_cb(j - 1)))
            end if

            if (k == yb) then
                ! dxix / dy
                tensora(3) = (-25d0*q_prim_vf(xibeg)%sf(j, k, l) &
                              + 48d0*q_prim_vf(xibeg)%sf(j, k + 1, l) &
                              - 36d0*q_prim_vf(xibeg)%sf(j, k + 2, l) &
                              + 16d0*q_prim_vf(xibeg)%sf(j, k + 3, l) &
                              - 3d0*q_prim_vf(xibeg)%sf(j, k + 4, l)) &
                             /(12d0*(y_cb(k + 1) - y_cb(k)))
            else if (k == yb + 1) then
                ! dxix / dy
                tensora(3) = (-3d0*q_prim_vf(xibeg)%sf(j, k - 1, l) &
                              - 10d0*q_prim_vf(xibeg)%sf(j, k, l) &
                              + 18d0*q_prim_vf(xibeg)%sf(j, k + 1, l) &
                              - 6d0*q_prim_vf(xibeg)%sf(j, k + 2, l) &
                              + q_prim_vf(xibeg)%sf(j, k + 3, l)) &
                             /(12d0*(y_cb(k) - y_cb(k - 1)))
            else if (k == ye - 1) then
                ! dxix / dy
                tensora(3) = (3d0*q_prim_vf(xibeg)%sf(j, k + 1, l) &
                              + 10d0*q_prim_vf(xibeg)%sf(j, k, l) &
                              - 18d0*q_prim_vf(xibeg)%sf(j, k - 1, l) &
                              + 6d0*q_prim_vf(xibeg)%sf(j, k - 2, l) &
                              - q_prim_vf(xibeg)%sf(j, k - 3, l)) &
                             /(12d0*(y_cb(k) - y_cb(k - 1)))
            else if (k == ye) then
                ! dxix / dy
                tensora(3) = (25d0*q_prim_vf(xibeg)%sf(j, k, l) &
                              - 48d0*q_prim_vf(xibeg)%sf(j, k - 1, l) &
                              + 36d0*q_prim_vf(xibeg)%sf(j, k - 2, l) &
                              - 16d0*q_prim_vf(xibeg)%sf(j, k - 3, l) &
                              + 3d0*q_prim_vf(xibeg)%sf(j, k - 4, l)) &
                             /(12d0*(y_cb(k) - y_cb(k - 1)))
            else
                ! dxix / dy
                tensora(3) = (q_prim_vf(xibeg)%sf(j, k - 2, l) &
                              - 8d0*q_prim_vf(xibeg)%sf(j, k - 1, l) &
                              + 8d0*q_prim_vf(xibeg)%sf(j, k + 1, l) &
                              - q_prim_vf(xibeg)%sf(j, k + 2, l)) &
                             /(12d0*(y_cb(k) - y_cb(k - 1)))
            end if

            if (k == yb) then
                ! dxiy / dy
                tensora(4) = (-25d0*q_prim_vf(xibeg + 1)%sf(j, k, l) &
                              + 48d0*q_prim_vf(xibeg + 1)%sf(j, k + 1, l) &
                              - 36d0*q_prim_vf(xibeg + 1)%sf(j, k + 2, l) &
                              + 16d0*q_prim_vf(xibeg + 1)%sf(j, k + 3, l) &
                              - 3d0*q_prim_vf(xibeg + 1)%sf(j, k + 4, l)) &
                             /(12d0*(y_cb(k + 1) - y_cb(k)))
            else if (k == yb + 1) then
                ! dxiy / dy
                tensora(4) = (-3d0*q_prim_vf(xibeg + 1)%sf(j, k - 1, l) &
                              - 10d0*q_prim_vf(xibeg + 1)%sf(j, k, l) &
                              + 18d0*q_prim_vf(xibeg + 1)%sf(j, k + 1, l) &
                              - 6d0*q_prim_vf(xibeg + 1)%sf(j, k + 2, l) &
                              + q_prim_vf(xibeg + 1)%sf(j, k + 3, l)) &
                             /(12d0*(y_cb(k) - y_cb(k - 1)))
            else if (k == ye - 1) then
                ! dxiy / dy
                tensora(4) = (3d0*q_prim_vf(xibeg + 1)%sf(j, k + 1, l) &
                              + 10d0*q_prim_vf(xibeg + 1)%sf(j, k, l) &
                              - 18d0*q_prim_vf(xibeg + 1)%sf(j, k - 1, l) &
                              + 6d0*q_prim_vf(xibeg + 1)%sf(j, k - 2, l) &
                              - q_prim_vf(xibeg + 1)%sf(j, k - 3, l)) &
                             /(12d0*(y_cb(k) - y_cb(k - 1)))
            else if (k == ye) then
                ! dxiy / dy
                tensora(4) = (25d0*q_prim_vf(xibeg + 1)%sf(j, k, l) &
                              - 48d0*q_prim_vf(xibeg + 1)%sf(j, k - 1, l) &
                              + 36d0*q_prim_vf(xibeg + 1)%sf(j, k - 2, l) &
                              - 16d0*q_prim_vf(xibeg + 1)%sf(j, k - 3, l) &
                              + 3d0*q_prim_vf(xibeg + 1)%sf(j, k - 4, l)) &
                             /(12d0*(y_cb(k) - y_cb(k - 1)))
            else
                ! dxiy / dy
                tensora(4) = (q_prim_vf(xibeg + 1)%sf(j, k - 2, l) &
                              - 8d0*q_prim_vf(xibeg + 1)%sf(j, k - 1, l) &
                              + 8d0*q_prim_vf(xibeg + 1)%sf(j, k + 1, l) &
                              - q_prim_vf(xibeg + 1)%sf(j, k + 2, l)) &
                             /(12d0*(y_cb(k) - y_cb(k - 1)))
            end if

        end if

        ! 3D
        if (num_dims > 2) then
            ! using results from upper if statement to map form 2x2 to 3x3 tensor
            tensora(5) = tensora(4)
            tensora(4) = tensora(3)

            if (l == zb) then
                ! dxix / dz
                tensora(7) = (-25d0*q_prim_vf(xibeg)%sf(j, k, l) &
                              + 48d0*q_prim_vf(xibeg)%sf(j, k, l + 1) &
                              - 36d0*q_prim_vf(xibeg)%sf(j, k, l + 2) &
                              + 16d0*q_prim_vf(xibeg)%sf(j, k, l + 3) &
                              - 3d0*q_prim_vf(xibeg)%sf(j, k, l + 4)) &
                             /(12d0*(z_cb(l + 1) - z_cb(l)))
            else if (l == zb + 1) then
                ! dxix / dz
                tensora(7) = (-3d0*q_prim_vf(xibeg)%sf(j, k, l - 1) &
                              - 10d0*q_prim_vf(xibeg)%sf(j, k, l) &
                              + 18d0*q_prim_vf(xibeg)%sf(j, k, l + 1) &
                              - 6d0*q_prim_vf(xibeg)%sf(j, k, l + 2) &
                              + q_prim_vf(xibeg)%sf(j, k, l + 3)) &
                             /(12d0*(z_cb(l) - z_cb(l - 1)))
            else if (l == ze - 1) then
                ! dxix / dz
                tensora(7) = (3d0*q_prim_vf(xibeg)%sf(j, k, l + 1) &
                              + 10d0*q_prim_vf(xibeg)%sf(j, k, l) &
                              - 18d0*q_prim_vf(xibeg)%sf(j, k, l - 1) &
                              + 6d0*q_prim_vf(xibeg)%sf(j, k, l - 2) &
                              - q_prim_vf(xibeg)%sf(j, k, l - 3)) &
                             /(12d0*(z_cb(l) - z_cb(l - 1)))
            else if (l == ze) then
                ! dxix / dz
                tensora(7) = (25d0*q_prim_vf(xibeg)%sf(j, k, l) &
                              - 48d0*q_prim_vf(xibeg)%sf(j, k, l - 1) &
                              + 36d0*q_prim_vf(xibeg)%sf(j, k, l - 2) &
                              - 16d0*q_prim_vf(xibeg)%sf(j, k, l - 3) &
                              + 3d0*q_prim_vf(xibeg)%sf(j, k, l - 4)) &
                             /(12d0*(z_cb(l) - z_cb(l - 1)))
            else
                ! dxix / dz
                tensora(7) = (q_prim_vf(xibeg)%sf(j, k, l - 2) &
                              - 8d0*q_prim_vf(xibeg)%sf(j, k, l - 1) &
                              + 8d0*q_prim_vf(xibeg)%sf(j, k, l + 1) &
                              - q_prim_vf(xibeg)%sf(j, k, l + 2)) &
                             /(12d0*(z_cb(l) - z_cb(l - 1)))
            end if

            if (l == zb) then
                ! dxiy / dz
                tensora(8) = (-25d0*q_prim_vf(xibeg + 1)%sf(j, k, l) &
                              + 48d0*q_prim_vf(xibeg + 1)%sf(j, k, l + 1) &
                              - 36d0*q_prim_vf(xibeg + 1)%sf(j, k, l + 2) &
                              + 16d0*q_prim_vf(xibeg + 1)%sf(j, k, l + 3) &
                              - 3d0*q_prim_vf(xibeg + 1)%sf(j, k, l + 4)) &
                             /(12d0*(z_cb(l + 1) - z_cb(l)))
            else if (l == zb + 1) then
                ! dxiy / dz
                tensora(8) = (-3d0*q_prim_vf(xibeg + 1)%sf(j, k, l - 1) &
                              - 10d0*q_prim_vf(xibeg + 1)%sf(j, k, l) &
                              + 18d0*q_prim_vf(xibeg + 1)%sf(j, k, l + 1) &
                              - 6d0*q_prim_vf(xibeg + 1)%sf(j, k, l + 2) &
                              + q_prim_vf(xibeg + 1)%sf(j, k, l + 3)) &
                             /(12d0*(z_cb(l) - z_cb(l - 1)))
            else if (l == ze - 1) then
                ! dxiy / dz
                tensora(8) = (3d0*q_prim_vf(xibeg + 1)%sf(j, k, l + 1) &
                              + 10d0*q_prim_vf(xibeg + 1)%sf(j, k, l) &
                              - 18d0*q_prim_vf(xibeg + 1)%sf(j, k, l - 1) &
                              + 6d0*q_prim_vf(xibeg + 1)%sf(j, k, l - 2) &
                              - q_prim_vf(xibeg + 1)%sf(j, k, l - 3)) &
                             /(12d0*(z_cb(l) - z_cb(l - 1)))
            else if (l == ze) then
                ! dxiy / dz
                tensora(8) = (25d0*q_prim_vf(xibeg + 1)%sf(j, k, l) &
                              - 48d0*q_prim_vf(xibeg + 1)%sf(j, k, l - 1) &
                              + 36d0*q_prim_vf(xibeg + 1)%sf(j, k, l - 2) &
                              - 16d0*q_prim_vf(xibeg + 1)%sf(j, k, l - 3) &
                              + 3d0*q_prim_vf(xibeg + 1)%sf(j, k, l - 4)) &
                             /(12d0*(z_cb(l) - z_cb(l - 1)))
            else
                ! dxiy / dz
                tensora(8) = (q_prim_vf(xibeg + 1)%sf(j, k, l - 2) &
                              - 8d0*q_prim_vf(xibeg + 1)%sf(j, k, l - 1) &
                              + 8d0*q_prim_vf(xibeg + 1)%sf(j, k, l + 1) &
                              - q_prim_vf(xibeg + 1)%sf(j, k, l + 2)) &
                             /(12d0*(z_cb(l) - z_cb(l - 1)))
            end if

            if (j == xb) then
                ! dxiz / dx
                tensora(3) = (-25d0*q_prim_vf(xiend)%sf(j, k, l) &
                              + 48d0*q_prim_vf(xiend)%sf(j + 1, k, l) &
                              - 36d0*q_prim_vf(xiend)%sf(j + 2, k, l) &
                              + 16d0*q_prim_vf(xiend)%sf(j + 3, k, l) &
                              - 3d0*q_prim_vf(xiend)%sf(j + 4, k, l)) &
                             /(12d0*(x_cb(j + 1) - x_cb(j)))
            else if (j == xb + 1) then
                ! dxiz / dx
                tensora(3) = (-3d0*q_prim_vf(xiend)%sf(j - 1, k, l) &
                              - 10d0*q_prim_vf(xiend)%sf(j, k, l) &
                              + 18d0*q_prim_vf(xiend)%sf(j + 1, k, l) &
                              - 6d0*q_prim_vf(xiend)%sf(j + 2, k, l) &
                              + q_prim_vf(xiend)%sf(j + 3, k, l)) &
                             /(12d0*(x_cb(j) - x_cb(j - 1)))
            else if (j == xe - 1) then
                ! dxiz / dx
                tensora(3) = (3d0*q_prim_vf(xiend)%sf(j + 1, k, l) &
                              + 10d0*q_prim_vf(xiend)%sf(j, k, l) &
                              - 18d0*q_prim_vf(xiend)%sf(j - 1, k, l) &
                              + 6d0*q_prim_vf(xiend)%sf(j - 2, k, l) &
                              - q_prim_vf(xiend)%sf(j - 3, k, l)) &
                             /(12d0*(x_cb(j) - x_cb(j - 1)))
            else if (j == xe) then
                ! dxiz / dx
                tensora(3) = (25d0*q_prim_vf(xiend)%sf(j, k, l) &
                              - 48d0*q_prim_vf(xiend)%sf(j - 1, k, l) &
                              + 36d0*q_prim_vf(xiend)%sf(j - 2, k, l) &
                              - 16d0*q_prim_vf(xiend)%sf(j - 3, k, l) &
                              + 3d0*q_prim_vf(xiend)%sf(j - 4, k, l)) &
                             /(12d0*(x_cb(j) - x_cb(j - 1)))
            else
                ! dxiz / dx
                tensora(3) = (q_prim_vf(xiend)%sf(j - 2, k, l) &
                              - 8d0*q_prim_vf(xiend)%sf(j - 1, k, l) &
                              + 8d0*q_prim_vf(xiend)%sf(j + 1, k, l) &
                              - q_prim_vf(xiend)%sf(j + 2, k, l)) &
                             /(12d0*(x_cb(j) - x_cb(j - 1)))
            end if

            if (k == yb) then
                ! dxiz / dy
                tensora(6) = (-25d0*q_prim_vf(xiend)%sf(j, k, l) &
                              + 48d0*q_prim_vf(xiend)%sf(j, k + 1, l) &
                              - 36d0*q_prim_vf(xiend)%sf(j, k + 2, l) &
                              + 16d0*q_prim_vf(xiend)%sf(j, k + 3, l) &
                              - 3d0*q_prim_vf(xiend)%sf(j, k + 4, l)) &
                             /(12d0*(y_cb(k + 1) - y_cb(k)))
            else if (k == yb + 1) then
                ! dxiz / dy
                tensora(6) = (-3d0*q_prim_vf(xiend)%sf(j, k - 1, l) &
                              - 10d0*q_prim_vf(xiend)%sf(j, k, l) &
                              + 18d0*q_prim_vf(xiend)%sf(j, k + 1, l) &
                              - 6d0*q_prim_vf(xiend)%sf(j, k + 2, l) &
                              + q_prim_vf(xiend)%sf(j, k + 3, l)) &
                             /(12d0*(y_cb(k) - y_cb(k - 1)))
            else if (k == ye - 1) then
                ! dxiz / dy
                tensora(6) = (3d0*q_prim_vf(xiend)%sf(j, k + 1, l) &
                              + 10d0*q_prim_vf(xiend)%sf(j, k, l) &
                              - 18d0*q_prim_vf(xiend)%sf(j, k - 1, l) &
                              + 6d0*q_prim_vf(xiend)%sf(j, k - 2, l) &
                              - q_prim_vf(xiend)%sf(j, k - 3, l)) &
                             /(12d0*(y_cb(k) - y_cb(k - 1)))
            else if (k == ye) then
                ! dxiz / dy
                tensora(6) = (25d0*q_prim_vf(xiend)%sf(j, k, l) &
                              - 48d0*q_prim_vf(xiend)%sf(j, k - 1, l) &
                              + 36d0*q_prim_vf(xiend)%sf(j, k - 2, l) &
                              - 16d0*q_prim_vf(xiend)%sf(j, k - 3, l) &
                              + 3d0*q_prim_vf(xiend)%sf(j, k - 4, l)) &
                             /(12d0*(y_cb(k) - y_cb(k - 1)))
            else
                ! dxiz / dy
                tensora(6) = (q_prim_vf(xiend)%sf(j, k - 2, l) &
                              - 8d0*q_prim_vf(xiend)%sf(j, k - 1, l) &
                              + 8d0*q_prim_vf(xiend)%sf(j, k + 1, l) &
                              - q_prim_vf(xiend)%sf(j, k + 2, l)) &
                             /(12d0*(y_cb(k) - y_cb(k - 1)))
            end if

            if (l == zb) then
                ! dxiz / dz
                tensora(9) = (-25d0*q_prim_vf(xiend)%sf(j, k, l) &
                              + 48d0*q_prim_vf(xiend)%sf(j, k, l + 1) &
                              - 36d0*q_prim_vf(xiend)%sf(j, k, l + 2) &
                              + 16d0*q_prim_vf(xiend)%sf(j, k, l + 3) &
                              - 3d0*q_prim_vf(xiend)%sf(j, k, l + 4)) &
                             /(12d0*(z_cb(l + 1) - z_cb(l)))
            else if (l == zb + 1) then
                ! dxiz / dz
                tensora(9) = (-3d0*q_prim_vf(xiend)%sf(j, k, l - 1) &
                              - 10d0*q_prim_vf(xiend)%sf(j, k, l) &
                              + 18d0*q_prim_vf(xiend)%sf(j, k, l + 1) &
                              - 6d0*q_prim_vf(xiend)%sf(j, k, l + 2) &
                              + q_prim_vf(xiend)%sf(j, k, l + 3)) &
                             /(12d0*(z_cb(l) - z_cb(l - 1)))
            else if (l == ze - 1) then
                ! dxiz / dz
                tensora(9) = (3d0*q_prim_vf(xiend)%sf(j, k, l + 1) &
                              + 10d0*q_prim_vf(xiend)%sf(j, k, l) &
                              - 18d0*q_prim_vf(xiend)%sf(j, k, l - 1) &
                              + 6d0*q_prim_vf(xiend)%sf(j, k, l - 2) &
                              - q_prim_vf(xiend)%sf(j, k, l - 3)) &
                             /(12d0*(z_cb(l) - z_cb(l - 1)))
            else if (l == ze) then
                ! dxiz / dz
                tensora(9) = (25d0*q_prim_vf(xiend)%sf(j, k, l) &
                              - 48d0*q_prim_vf(xiend)%sf(j, k, l - 1) &
                              + 36d0*q_prim_vf(xiend)%sf(j, k, l - 2) &
                              - 16d0*q_prim_vf(xiend)%sf(j, k, l - 3) &
                              + 3d0*q_prim_vf(xiend)%sf(j, k, l - 4)) &
                             /(12d0*(z_cb(l) - z_cb(l - 1)))
            else
                ! dxiz / dz
                tensora(9) = (q_prim_vf(xiend)%sf(j, k, l - 2) &
                              - 8d0*q_prim_vf(xiend)%sf(j, k, l - 1) &
                              + 8d0*q_prim_vf(xiend)%sf(j, k, l + 1) &
                              - q_prim_vf(xiend)%sf(j, k, l + 2)) &
                             /(12d0*(z_cb(l) - z_cb(l - 1)))
            end if
        end if

        ! STEP 2a: computing the adjoint of the grad_xi tensor for the inverse
        if (num_dims == 1) then
            tensorb(1) = 1
        elseif (num_dims == 2) then
            tensorb(1) = tensora(4)
            tensorb(2) = -tensora(3)
            tensorb(3) = -tensora(2)
            tensorb(4) = tensora(1)
        elseif (num_dims == 3) then
            tensorb(1) = tensora(5)*tensora(9) - tensora(6)*tensora(8)
            tensorb(2) = -(tensora(2)*tensora(9) - tensora(3)*tensora(8))
            tensorb(3) = tensora(2)*tensora(6) - tensora(3)*tensora(5)
            tensorb(4) = -(tensora(4)*tensora(9) - tensora(6)*tensora(7))
            tensorb(5) = tensora(1)*tensora(9) - tensora(3)*tensora(7)
            tensorb(6) = -(tensora(1)*tensora(6) - tensora(4)*tensora(3))
            tensorb(7) = tensora(4)*tensora(8) - tensora(5)*tensora(7)
            tensorb(8) = -(tensora(1)*tensora(8) - tensora(2)*tensora(7))
            tensorb(9) = tensora(1)*tensora(5) - tensora(2)*tensora(4)
        end if

        ! STEP 2b: computing the determinant of the grad_xi tensor
        if (num_dims == 1) then
            determinant = tensora(1)
        elseif (num_dims == 2) then
            determinant = tensora(1)*tensora(4) - tensora(2)*tensora(3)
        else
            determinant = tensora(1)*(tensora(5)*tensora(9) - tensora(6)*tensora(8)) &
                          - tensora(2)*(tensora(4)*tensora(9) - tensora(6)*tensora(7)) &
                          + tensora(3)*(tensora(4)*tensora(8) - tensora(5)*tensora(7))
        end if
        ! error checking
        !if (determinant == 0) then
        !    print *, 'determinant :: ', determinant
        !    print *, 'ERROR: Determinant was zero'
        !    stop
        !end if
        if (determinant < 0d0 .or. determinant > 2d0) then
            print *, 'i, j, k :: ', j, ' ', k, ' ', l, ',det ::', tensorb(tensor_size)
            !    stop
        end if

        ! STEP 2c: computing the inverse of grad_xi tensor = F
        ! tensorb is the adjoint, tensora becomes the inverse
        do i = 1, tensor_size - 1
            tensora(i) = tensorb(i)/determinant
        end do

        ! STEP 3: computing F tranpose F
        tensorb(1) = tensora(1)**2
        if (num_dims == 2) then
            tensorb(1) = tensorb(1) + tensora(3)**2
            tensorb(2) = tensora(1)*tensora(2) + tensora(3)*tensora(4)
            tensorb(3) = tensorb(2)
            tensorb(4) = tensora(2)**2 + tensora(4)**2
        elseif (num_dims == 3) then
            tensorb(1) = tensora(1)**2 + tensora(2)**2 + tensora(3)**2
            tensorb(5) = tensora(4)**2 + tensora(5)**2 + tensora(6)**2
            tensorb(9) = tensora(7)**2 + tensora(8)**2 + tensora(9)**2
            tensorb(2) = tensora(1)*tensora(4) + tensora(2)*tensora(5) + tensora(3)*tensora(6)
            tensorb(3) = tensora(1)*tensora(7) + tensora(2)*tensora(8) + tensora(3)*tensora(9)
            tensorb(6) = tensora(4)*tensora(7) + tensora(5)*tensora(8) + tensora(6)*tensora(9)
            tensorb(4) = tensorb(2)
            tensorb(7) = tensorb(3)
            tensorb(8) = tensorb(6)
        end if
        ! STEP 4: store the determinant of F in the last entry of the tensor
        tensorb(tensor_size) = determinant

    end subroutine s_compute_gradient_xi

    !>  The following subroutine handles the calculation of the btensor.
        !!   The calculation of the btensor takes qprimvf.
        !! @param q_prim_vf Primitive variables
        !! @param btensor is the output
        !! calculate the grad_xi, grad_xi is a nxn tensor
        !! calculate the inverse of grad_xi to obtain F, F is a nxn tensor
        !! calculate the FFtranspose to obtain the btensor, btensor is nxn tensor
        !! btensor is symmetric, save the data space
    subroutine s_compute_gradient_xi1d_acc(q_prim_vf, ixb, ixe, iyb, iye, & !---------
                                           izb, ize, j, k, l, tensora, tensorb)
        !$acc routine seq
        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf
        real(kind(0d0)), dimension(tensor_size), intent(INOUT) :: tensora
        real(kind(0d0)), dimension(tensor_size), intent(INOUT) :: tensorb
        integer, intent(IN) :: ixb, ixe, iyb, iye, izb, ize
        integer, intent(IN) :: j, k, l
        integer :: i

        ! STEP 1: computing the grad_xi tensor
        ! grad_xi definition / organization
        ! number for the tensor 1-3:  dxix_dx, dxiy_dx, dxiz_dx
        ! 4-6 :                       dxix_dy, dxiy_dy, dxiz_dy
        ! 7-9 :                       dxix_dz, dxiy_dz, dxiz_dz
        if (j == ixb) then
            ! dxix/dx
            tensora(1) = (-25d0*q_prim_vf(xibeg)%sf(j, k, l) &
                          + 48d0*q_prim_vf(xibeg)%sf(j + 1, k, l) &
                          - 36d0*q_prim_vf(xibeg)%sf(j + 2, k, l) &
                          + 16d0*q_prim_vf(xibeg)%sf(j + 3, k, l) &
                          - 3d0*q_prim_vf(xibeg)%sf(j + 4, k, l)) &
                         /(12d0*(x_cb(j + 1) - x_cb(j)))
        else if (j == ixb + 1) then
            ! dxix/dx
            tensora(1) = (-3d0*q_prim_vf(xibeg)%sf(j - 1, k, l) &
                          - 10d0*q_prim_vf(xibeg)%sf(j, k, l) &
                          + 18d0*q_prim_vf(xibeg)%sf(j + 1, k, l) &
                          - 6d0*q_prim_vf(xibeg)%sf(j + 2, k, l) &
                          + q_prim_vf(xibeg)%sf(j + 3, k, l)) &
                         /(12d0*(x_cb(j) - x_cb(j - 1)))
        else if (j == ixe - 1) then
            ! dxix/dx
            tensora(1) = (3d0*q_prim_vf(xibeg)%sf(j + 1, k, l) &
                          + 10d0*q_prim_vf(xibeg)%sf(j, k, l) &
                          - 18d0*q_prim_vf(xibeg)%sf(j - 1, k, l) &
                          + 6d0*q_prim_vf(xibeg)%sf(j - 2, k, l) &
                          - q_prim_vf(xibeg)%sf(j - 3, k, l)) &
                         /(12d0*(x_cb(j) - x_cb(j - 1)))
        else if (j == ixe) then
            ! dxix/dx
            tensora(1) = (25d0*q_prim_vf(xibeg)%sf(j, k, l) &
                          - 48d0*q_prim_vf(xibeg)%sf(j - 1, k, l) &
                          + 36d0*q_prim_vf(xibeg)%sf(j - 2, k, l) &
                          - 16d0*q_prim_vf(xibeg)%sf(j - 3, k, l) &
                          + 3d0*q_prim_vf(xibeg)%sf(j - 4, k, l)) &
                         /(12d0*(x_cb(j) - x_cb(j - 1)))
        else
            ! dxix/dx
            tensora(1) = (q_prim_vf(xibeg)%sf(j - 2, k, l) &
                          - 8d0*q_prim_vf(xibeg)%sf(j - 1, k, l) &
                          + 8d0*q_prim_vf(xibeg)%sf(j + 1, k, l) &
                          - q_prim_vf(xibeg)%sf(j + 2, k, l)) &
                         /(12d0*(x_cb(j) - x_cb(j - 1)))
        end if

        ! STEP 2a: computing the adjoint of the grad_xi tensor for the inverse
        tensorb(1) = 1

        ! STEP 2b: computing the determinant of the grad_xi tensor
        tensorb(tensor_size) = tensora(1)

        ! STEP 2c: computing the inverse of grad_xi tensor = F
        ! tensorb is the adjoint, tensora becomes the inverse

        !$acc loop seq
        do i = 1, tensor_size - 1
            tensora(i) = tensorb(i)/tensorb(tensor_size)
        end do

        ! STEP 3: computing F tranpose F
        tensorb(1) = tensora(1)**2
        ! STEP 4: store the determinant of F in the last entry of the tensor
        !tensorb(tensor_size) = determinant

    end subroutine s_compute_gradient_xi1d_acc

    !>  The following subroutine handles the calculation of the btensor.
        !!   The calculation of the btensor takes qprimvf.
        !! @param q_prim_vf Primitive variables
        !! @param btensor is the output
        !! calculate the grad_xi, grad_xi is a nxn tensor
        !! calculate the inverse of grad_xi to obtain F, F is a nxn tensor
        !! calculate the FFtranspose to obtain the btensor, btensor is nxn tensor
        !! btensor is symmetric, save the data space
    subroutine s_compute_gradient_xi2d_acc(q_prim_vf, ixb, ixe, iyb, iye, & !---------
                                           izb, ize, j, k, l, tensora, tensorb)
        !$acc routine seq
        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf
        real(kind(0d0)), dimension(tensor_size), intent(INOUT) :: tensora
        real(kind(0d0)), dimension(tensor_size), intent(INOUT) :: tensorb
        integer, intent(IN) :: ixb, ixe, iyb, iye, izb, ize
        integer, intent(IN) :: j, k, l
        integer :: i

        ! STEP 1: computing the grad_xi tensor
        ! grad_xi definition / organization
        ! number for the tensor 1-3:  dxix_dx, dxiy_dx, dxiz_dx
        ! 4-6 :                       dxix_dy, dxiy_dy, dxiz_dy
        ! 7-9 :                       dxix_dz, dxiy_dz, dxiz_dz
        if (j == ixb) then
            ! dxix/dx
            tensora(1) = (-25d0*q_prim_vf(xibeg)%sf(j, k, l) &
                          + 48d0*q_prim_vf(xibeg)%sf(j + 1, k, l) &
                          - 36d0*q_prim_vf(xibeg)%sf(j + 2, k, l) &
                          + 16d0*q_prim_vf(xibeg)%sf(j + 3, k, l) &
                          - 3d0*q_prim_vf(xibeg)%sf(j + 4, k, l)) &
                         /(12d0*(x_cb(j + 1) - x_cb(j)))
        else if (j == ixb + 1) then
            ! dxix/dx
            tensora(1) = (-3d0*q_prim_vf(xibeg)%sf(j - 1, k, l) &
                          - 10d0*q_prim_vf(xibeg)%sf(j, k, l) &
                          + 18d0*q_prim_vf(xibeg)%sf(j + 1, k, l) &
                          - 6d0*q_prim_vf(xibeg)%sf(j + 2, k, l) &
                          + q_prim_vf(xibeg)%sf(j + 3, k, l)) &
                         /(12d0*(x_cb(j) - x_cb(j - 1)))
        else if (j == ixe - 1) then
            ! dxix/dx
            tensora(1) = (3d0*q_prim_vf(xibeg)%sf(j + 1, k, l) &
                          + 10d0*q_prim_vf(xibeg)%sf(j, k, l) &
                          - 18d0*q_prim_vf(xibeg)%sf(j - 1, k, l) &
                          + 6d0*q_prim_vf(xibeg)%sf(j - 2, k, l) &
                          - q_prim_vf(xibeg)%sf(j - 3, k, l)) &
                         /(12d0*(x_cb(j) - x_cb(j - 1)))
        else if (j == ixe) then
            ! dxix/dx
            tensora(1) = (25d0*q_prim_vf(xibeg)%sf(j, k, l) &
                          - 48d0*q_prim_vf(xibeg)%sf(j - 1, k, l) &
                          + 36d0*q_prim_vf(xibeg)%sf(j - 2, k, l) &
                          - 16d0*q_prim_vf(xibeg)%sf(j - 3, k, l) &
                          + 3d0*q_prim_vf(xibeg)%sf(j - 4, k, l)) &
                         /(12d0*(x_cb(j) - x_cb(j - 1)))
        else
            ! dxix/dx
            tensora(1) = (q_prim_vf(xibeg)%sf(j - 2, k, l) &
                          - 8d0*q_prim_vf(xibeg)%sf(j - 1, k, l) &
                          + 8d0*q_prim_vf(xibeg)%sf(j + 1, k, l) &
                          - q_prim_vf(xibeg)%sf(j + 2, k, l)) &
                         /(12d0*(x_cb(j) - x_cb(j - 1)))
        end if

        ! 2D
        if (j == ixb) then
            ! dxiy / dx
            tensora(2) = (-25d0*q_prim_vf(xibeg + 1)%sf(j, k, l) &
                          + 48d0*q_prim_vf(xibeg + 1)%sf(j + 1, k, l) &
                          - 36d0*q_prim_vf(xibeg + 1)%sf(j + 2, k, l) &
                          + 16d0*q_prim_vf(xibeg + 1)%sf(j + 3, k, l) &
                          - 3d0*q_prim_vf(xibeg + 1)%sf(j + 4, k, l)) &
                         /(12d0*(x_cb(j + 1) - x_cb(j)))
        else if (j == ixb + 1) then
            ! dxiy / dx
            tensora(2) = (-3d0*q_prim_vf(xibeg + 1)%sf(j - 1, k, l) &
                          - 10d0*q_prim_vf(xibeg + 1)%sf(j, k, l) &
                          + 18d0*q_prim_vf(xibeg + 1)%sf(j + 1, k, l) &
                          - 6d0*q_prim_vf(xibeg + 1)%sf(j + 2, k, l) &
                          + q_prim_vf(xibeg + 1)%sf(j + 3, k, l)) &
                         /(12d0*(x_cb(j) - x_cb(j - 1)))
        else if (j == ixe - 1) then
            ! dxiy / dx
            tensora(2) = (3d0*q_prim_vf(xibeg + 1)%sf(j + 1, k, l) &
                          + 10d0*q_prim_vf(xibeg + 1)%sf(j, k, l) &
                          - 18d0*q_prim_vf(xibeg + 1)%sf(j - 1, k, l) &
                          + 6d0*q_prim_vf(xibeg + 1)%sf(j - 2, k, l) &
                          - q_prim_vf(xibeg + 1)%sf(j - 3, k, l)) &
                         /(12d0*(x_cb(j) - x_cb(j - 1)))
        else if (j == ixe) then
            ! dxiy / dx
            tensora(2) = (25d0*q_prim_vf(xibeg + 1)%sf(j, k, l) &
                          - 48d0*q_prim_vf(xibeg + 1)%sf(j - 1, k, l) &
                          + 36d0*q_prim_vf(xibeg + 1)%sf(j - 2, k, l) &
                          - 16d0*q_prim_vf(xibeg + 1)%sf(j - 3, k, l) &
                          + 3d0*q_prim_vf(xibeg + 1)%sf(j - 4, k, l)) &
                         /(12d0*(x_cb(j) - x_cb(j - 1)))
        else
            ! dxiy / dx
            tensora(2) = (q_prim_vf(xibeg + 1)%sf(j - 2, k, l) &
                          - 8d0*q_prim_vf(xibeg + 1)%sf(j - 1, k, l) &
                          + 8d0*q_prim_vf(xibeg + 1)%sf(j + 1, k, l) &
                          - q_prim_vf(xibeg + 1)%sf(j + 2, k, l)) &
                         /(12d0*(x_cb(j) - x_cb(j - 1)))
        end if

        if (k == iyb) then
            ! dxix / dy
            tensora(3) = (-25d0*q_prim_vf(xibeg)%sf(j, k, l) &
                          + 48d0*q_prim_vf(xibeg)%sf(j, k + 1, l) &
                          - 36d0*q_prim_vf(xibeg)%sf(j, k + 2, l) &
                          + 16d0*q_prim_vf(xibeg)%sf(j, k + 3, l) &
                          - 3d0*q_prim_vf(xibeg)%sf(j, k + 4, l)) &
                         /(12d0*(y_cb(k + 1) - y_cb(k)))
        else if (k == iyb + 1) then
            ! dxix / dy
            tensora(3) = (-3d0*q_prim_vf(xibeg)%sf(j, k - 1, l) &
                          - 10d0*q_prim_vf(xibeg)%sf(j, k, l) &
                          + 18d0*q_prim_vf(xibeg)%sf(j, k + 1, l) &
                          - 6d0*q_prim_vf(xibeg)%sf(j, k + 2, l) &
                          + q_prim_vf(xibeg)%sf(j, k + 3, l)) &
                         /(12d0*(y_cb(k) - y_cb(k - 1)))
        else if (k == iye - 1) then
            ! dxix / dy
            tensora(3) = (3d0*q_prim_vf(xibeg)%sf(j, k + 1, l) &
                          + 10d0*q_prim_vf(xibeg)%sf(j, k, l) &
                          - 18d0*q_prim_vf(xibeg)%sf(j, k - 1, l) &
                          + 6d0*q_prim_vf(xibeg)%sf(j, k - 2, l) &
                          - q_prim_vf(xibeg)%sf(j, k - 3, l)) &
                         /(12d0*(y_cb(k) - y_cb(k - 1)))
        else if (k == iye) then
            ! dxix / dy
            tensora(3) = (25d0*q_prim_vf(xibeg)%sf(j, k, l) &
                          - 48d0*q_prim_vf(xibeg)%sf(j, k - 1, l) &
                          + 36d0*q_prim_vf(xibeg)%sf(j, k - 2, l) &
                          - 16d0*q_prim_vf(xibeg)%sf(j, k - 3, l) &
                          + 3d0*q_prim_vf(xibeg)%sf(j, k - 4, l)) &
                         /(12d0*(y_cb(k) - y_cb(k - 1)))
        else
            ! dxix / dy
            tensora(3) = (q_prim_vf(xibeg)%sf(j, k - 2, l) &
                          - 8d0*q_prim_vf(xibeg)%sf(j, k - 1, l) &
                          + 8d0*q_prim_vf(xibeg)%sf(j, k + 1, l) &
                          - q_prim_vf(xibeg)%sf(j, k + 2, l)) &
                         /(12d0*(y_cb(k) - y_cb(k - 1)))
        end if

        if (k == iyb) then
            ! dxiy / dy
            tensora(4) = (-25d0*q_prim_vf(xibeg + 1)%sf(j, k, l) &
                          + 48d0*q_prim_vf(xibeg + 1)%sf(j, k + 1, l) &
                          - 36d0*q_prim_vf(xibeg + 1)%sf(j, k + 2, l) &
                          + 16d0*q_prim_vf(xibeg + 1)%sf(j, k + 3, l) &
                          - 3d0*q_prim_vf(xibeg + 1)%sf(j, k + 4, l)) &
                         /(12d0*(y_cb(k + 1) - y_cb(k)))
        else if (k == iyb + 1) then
            ! dxiy / dy
            tensora(4) = (-3d0*q_prim_vf(xibeg + 1)%sf(j, k - 1, l) &
                          - 10d0*q_prim_vf(xibeg + 1)%sf(j, k, l) &
                          + 18d0*q_prim_vf(xibeg + 1)%sf(j, k + 1, l) &
                          - 6d0*q_prim_vf(xibeg + 1)%sf(j, k + 2, l) &
                          + q_prim_vf(xibeg + 1)%sf(j, k + 3, l)) &
                         /(12d0*(y_cb(k) - y_cb(k - 1)))
        else if (k == iye - 1) then
            ! dxiy / dy
            tensora(4) = (3d0*q_prim_vf(xibeg + 1)%sf(j, k + 1, l) &
                          + 10d0*q_prim_vf(xibeg + 1)%sf(j, k, l) &
                          - 18d0*q_prim_vf(xibeg + 1)%sf(j, k - 1, l) &
                          + 6d0*q_prim_vf(xibeg + 1)%sf(j, k - 2, l) &
                          - q_prim_vf(xibeg + 1)%sf(j, k - 3, l)) &
                         /(12d0*(y_cb(k) - y_cb(k - 1)))
        else if (k == iye) then
            ! dxiy / dy
            tensora(4) = (25d0*q_prim_vf(xibeg + 1)%sf(j, k, l) &
                          - 48d0*q_prim_vf(xibeg + 1)%sf(j, k - 1, l) &
                          + 36d0*q_prim_vf(xibeg + 1)%sf(j, k - 2, l) &
                          - 16d0*q_prim_vf(xibeg + 1)%sf(j, k - 3, l) &
                          + 3d0*q_prim_vf(xibeg + 1)%sf(j, k - 4, l)) &
                         /(12d0*(y_cb(k) - y_cb(k - 1)))
        else
            ! dxiy / dy
            tensora(4) = (q_prim_vf(xibeg + 1)%sf(j, k - 2, l) &
                          - 8d0*q_prim_vf(xibeg + 1)%sf(j, k - 1, l) &
                          + 8d0*q_prim_vf(xibeg + 1)%sf(j, k + 1, l) &
                          - q_prim_vf(xibeg + 1)%sf(j, k + 2, l)) &
                         /(12d0*(y_cb(k) - y_cb(k - 1)))
        end if

        ! STEP 2a: computing the adjoint of the grad_xi tensor for the inverse
        tensorb(1) = tensora(4)
        tensorb(2) = -tensora(3)
        tensorb(3) = -tensora(2)
        tensorb(4) = tensora(1)

        ! STEP 2b: computing the determinant of the grad_xi tensor
        tensorb(tensor_size) = tensora(1)*tensora(4) - tensora(2)*tensora(3)

        ! STEP 2c: computing the inverse of grad_xi tensor = F
        ! tensorb is the adjoint, tensora becomes the inverse
        !$acc loop seq
        do i = 1, tensor_size - 1
            tensora(i) = tensorb(i)/tensorb(tensor_size)
        end do
        ! STEP 3: computing F tranpose F
        tensorb(1) = tensora(1)**2
        tensorb(1) = tensorb(1) + tensora(3)**2
        tensorb(2) = tensora(1)*tensora(2) + tensora(3)*tensora(4)
        tensorb(3) = tensorb(2)
        tensorb(4) = tensora(2)**2 + tensora(4)**2

        ! STEP 4: store the determinant of F in the last entry of the tensor

    end subroutine s_compute_gradient_xi2d_acc

    !>  The following subroutine handles the calculation of the btensor.
        !!   The calculation of the btensor takes qprimvf.
        !! @param q_prim_vf Primitive variables
        !! @param btensor is the output
        !! calculate the grad_xi, grad_xi is a nxn tensor
        !! calculate the inverse of grad_xi to obtain F, F is a nxn tensor
        !! calculate the FFtranspose to obtain the btensor, btensor is nxn tensor
        !! btensor is symmetric, save the data space
    subroutine s_compute_gradient_xi3d_acc(q_prim_vf, ixb, ixe, iyb, iye, & !---------
                                           izb, ize, j, k, l, tensora, tensorb)
        !$acc routine seq
        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf
        real(kind(0d0)), dimension(tensor_size), intent(INOUT) :: tensora
        real(kind(0d0)), dimension(tensor_size), intent(INOUT) :: tensorb
        integer, intent(IN) :: ixb, ixe
        integer, intent(IN) :: iyb, iye
        integer, intent(IN) :: izb, ize
        integer, intent(IN) :: j, k, l

        integer :: i

        ! STEP 1: computing the grad_xi tensor
        ! grad_xi definition / organization
        ! number for the tensor 1-3:  dxix_dx, dxiy_dx, dxiz_dx
        ! 4-6 :                       dxix_dy, dxiy_dy, dxiz_dy
        ! 7-9 :                       dxix_dz, dxiy_dz, dxiz_dz

        ! 1D
        if (j == ixb) then
            ! dxix/dx
            tensora(1) = (-25d0*q_prim_vf(xibeg)%sf(j, k, l) &
                          + 48d0*q_prim_vf(xibeg)%sf(j + 1, k, l) &
                          - 36d0*q_prim_vf(xibeg)%sf(j + 2, k, l) &
                          + 16d0*q_prim_vf(xibeg)%sf(j + 3, k, l) &
                          - 3d0*q_prim_vf(xibeg)%sf(j + 4, k, l)) &
                         /(12d0*(x_cb(j + 1) - x_cb(j)))
            ! dxiy / dx
            tensora(2) = (-25d0*q_prim_vf(xibeg + 1)%sf(j, k, l) &
                          + 48d0*q_prim_vf(xibeg + 1)%sf(j + 1, k, l) &
                          - 36d0*q_prim_vf(xibeg + 1)%sf(j + 2, k, l) &
                          + 16d0*q_prim_vf(xibeg + 1)%sf(j + 3, k, l) &
                          - 3d0*q_prim_vf(xibeg + 1)%sf(j + 4, k, l)) &
                         /(12d0*(x_cb(j + 1) - x_cb(j)))
            ! dxiz / dx
            tensora(7) = (-25d0*q_prim_vf(xiend)%sf(j, k, l) &
                          + 48d0*q_prim_vf(xiend)%sf(j + 1, k, l) &
                          - 36d0*q_prim_vf(xiend)%sf(j + 2, k, l) &
                          + 16d0*q_prim_vf(xiend)%sf(j + 3, k, l) &
                          - 3d0*q_prim_vf(xiend)%sf(j + 4, k, l)) &
                         /(12d0*(x_cb(j + 1) - x_cb(j)))

        else if (j == ixb + 1) then
            ! dxix/dx
            tensora(1) = (-3d0*q_prim_vf(xibeg)%sf(j - 1, k, l) &
                          - 10d0*q_prim_vf(xibeg)%sf(j, k, l) &
                          + 18d0*q_prim_vf(xibeg)%sf(j + 1, k, l) &
                          - 6d0*q_prim_vf(xibeg)%sf(j + 2, k, l) &
                          + q_prim_vf(xibeg)%sf(j + 3, k, l)) &
                         /(12d0*(x_cb(j) - x_cb(j - 1)))
            ! dxiy / dx
            tensora(2) = (-3d0*q_prim_vf(xibeg + 1)%sf(j - 1, k, l) &
                          - 10d0*q_prim_vf(xibeg + 1)%sf(j, k, l) &
                          + 18d0*q_prim_vf(xibeg + 1)%sf(j + 1, k, l) &
                          - 6d0*q_prim_vf(xibeg + 1)%sf(j + 2, k, l) &
                          + q_prim_vf(xibeg + 1)%sf(j + 3, k, l)) &
                         /(12d0*(x_cb(j) - x_cb(j - 1)))
            ! dxiz / dx
            tensora(7) = (-3d0*q_prim_vf(xiend)%sf(j - 1, k, l) &
                          - 10d0*q_prim_vf(xiend)%sf(j, k, l) &
                          + 18d0*q_prim_vf(xiend)%sf(j + 1, k, l) &
                          - 6d0*q_prim_vf(xiend)%sf(j + 2, k, l) &
                          + q_prim_vf(xiend)%sf(j + 3, k, l)) &
                         /(12d0*(x_cb(j) - x_cb(j - 1)))

        else if (j == ixe - 1) then
            ! dxix/dx
            tensora(1) = (3d0*q_prim_vf(xibeg)%sf(j + 1, k, l) &
                          + 10d0*q_prim_vf(xibeg)%sf(j, k, l) &
                          - 18d0*q_prim_vf(xibeg)%sf(j - 1, k, l) &
                          + 6d0*q_prim_vf(xibeg)%sf(j - 2, k, l) &
                          - q_prim_vf(xibeg)%sf(j - 3, k, l)) &
                         /(12d0*(x_cb(j) - x_cb(j - 1)))
            ! dxiy / dx
            tensora(2) = (3d0*q_prim_vf(xibeg + 1)%sf(j + 1, k, l) &
                          + 10d0*q_prim_vf(xibeg + 1)%sf(j, k, l) &
                          - 18d0*q_prim_vf(xibeg + 1)%sf(j - 1, k, l) &
                          + 6d0*q_prim_vf(xibeg + 1)%sf(j - 2, k, l) &
                          - q_prim_vf(xibeg + 1)%sf(j - 3, k, l)) &
                         /(12d0*(x_cb(j) - x_cb(j - 1)))
            ! dxiz / dx
            tensora(7) = (3d0*q_prim_vf(xiend)%sf(j + 1, k, l) &
                          + 10d0*q_prim_vf(xiend)%sf(j, k, l) &
                          - 18d0*q_prim_vf(xiend)%sf(j - 1, k, l) &
                          + 6d0*q_prim_vf(xiend)%sf(j - 2, k, l) &
                          - q_prim_vf(xiend)%sf(j - 3, k, l)) &
                         /(12d0*(x_cb(j) - x_cb(j - 1)))
        else if (j == ixe) then
            ! dxix/dx
            tensora(1) = (25d0*q_prim_vf(xibeg)%sf(j, k, l) &
                          - 48d0*q_prim_vf(xibeg)%sf(j - 1, k, l) &
                          + 36d0*q_prim_vf(xibeg)%sf(j - 2, k, l) &
                          - 16d0*q_prim_vf(xibeg)%sf(j - 3, k, l) &
                          + 3d0*q_prim_vf(xibeg)%sf(j - 4, k, l)) &
                         /(12d0*(x_cb(j) - x_cb(j - 1)))
            ! dxiy / dx
            tensora(2) = (25d0*q_prim_vf(xibeg + 1)%sf(j, k, l) &
                          - 48d0*q_prim_vf(xibeg + 1)%sf(j - 1, k, l) &
                          + 36d0*q_prim_vf(xibeg + 1)%sf(j - 2, k, l) &
                          - 16d0*q_prim_vf(xibeg + 1)%sf(j - 3, k, l) &
                          + 3d0*q_prim_vf(xibeg + 1)%sf(j - 4, k, l)) &
                         /(12d0*(x_cb(j) - x_cb(j - 1)))
            ! dxiz / dx
            tensora(7) = (25d0*q_prim_vf(xiend)%sf(j, k, l) &
                          - 48d0*q_prim_vf(xiend)%sf(j - 1, k, l) &
                          + 36d0*q_prim_vf(xiend)%sf(j - 2, k, l) &
                          - 16d0*q_prim_vf(xiend)%sf(j - 3, k, l) &
                          + 3d0*q_prim_vf(xiend)%sf(j - 4, k, l)) &
                         /(12d0*(x_cb(j) - x_cb(j - 1)))
        else
            ! dxix/dx
            tensora(1) = (q_prim_vf(xibeg)%sf(j - 2, k, l) &
                          - 8d0*q_prim_vf(xibeg)%sf(j - 1, k, l) &
                          + 8d0*q_prim_vf(xibeg)%sf(j + 1, k, l) &
                          - q_prim_vf(xibeg)%sf(j + 2, k, l)) &
                         /(12d0*(x_cb(j) - x_cb(j - 1)))
            ! dxiy / dx
            tensora(2) = (q_prim_vf(xibeg + 1)%sf(j - 2, k, l) &
                          - 8d0*q_prim_vf(xibeg + 1)%sf(j - 1, k, l) &
                          + 8d0*q_prim_vf(xibeg + 1)%sf(j + 1, k, l) &
                          - q_prim_vf(xibeg + 1)%sf(j + 2, k, l)) &
                         /(12d0*(x_cb(j) - x_cb(j - 1)))
            ! dxiz / dx
            tensora(7) = (q_prim_vf(xiend)%sf(j - 2, k, l) &
                          - 8d0*q_prim_vf(xiend)%sf(j - 1, k, l) &
                          + 8d0*q_prim_vf(xiend)%sf(j + 1, k, l) &
                          - q_prim_vf(xiend)%sf(j + 2, k, l)) &
                         /(12d0*(x_cb(j) - x_cb(j - 1)))

        end if

        ! 2D
        if (k == iyb) then
            ! dxix / dy
            tensora(4) = (-25d0*q_prim_vf(xibeg)%sf(j, k, l) &
                          + 48d0*q_prim_vf(xibeg)%sf(j, k + 1, l) &
                          - 36d0*q_prim_vf(xibeg)%sf(j, k + 2, l) &
                          + 16d0*q_prim_vf(xibeg)%sf(j, k + 3, l) &
                          - 3d0*q_prim_vf(xibeg)%sf(j, k + 4, l)) &
                         /(12d0*(y_cb(k + 1) - y_cb(k)))
            ! dxiy / dy
            tensora(5) = (-25d0*q_prim_vf(xibeg + 1)%sf(j, k, l) &
                          + 48d0*q_prim_vf(xibeg + 1)%sf(j, k + 1, l) &
                          - 36d0*q_prim_vf(xibeg + 1)%sf(j, k + 2, l) &
                          + 16d0*q_prim_vf(xibeg + 1)%sf(j, k + 3, l) &
                          - 3d0*q_prim_vf(xibeg + 1)%sf(j, k + 4, l)) &
                         /(12d0*(y_cb(k + 1) - y_cb(k)))
            ! dxiz / dy
            tensora(8) = (-25d0*q_prim_vf(xiend)%sf(j, k, l) &
                          + 48d0*q_prim_vf(xiend)%sf(j, k + 1, l) &
                          - 36d0*q_prim_vf(xiend)%sf(j, k + 2, l) &
                          + 16d0*q_prim_vf(xiend)%sf(j, k + 3, l) &
                          - 3d0*q_prim_vf(xiend)%sf(j, k + 4, l)) &
                         /(12d0*(y_cb(k + 1) - y_cb(k)))

        else if (k == iyb + 1) then
            ! dxix / dy
            tensora(4) = (-3d0*q_prim_vf(xibeg)%sf(j, k - 1, l) &
                          - 10d0*q_prim_vf(xibeg)%sf(j, k, l) &
                          + 18d0*q_prim_vf(xibeg)%sf(j, k + 1, l) &
                          - 6d0*q_prim_vf(xibeg)%sf(j, k + 2, l) &
                          + q_prim_vf(xibeg)%sf(j, k + 3, l)) &
                         /(12d0*(y_cb(k) - y_cb(k - 1)))
            ! dxiy / dy
            tensora(5) = (-3d0*q_prim_vf(xibeg + 1)%sf(j, k - 1, l) &
                          - 10d0*q_prim_vf(xibeg + 1)%sf(j, k, l) &
                          + 18d0*q_prim_vf(xibeg + 1)%sf(j, k + 1, l) &
                          - 6d0*q_prim_vf(xibeg + 1)%sf(j, k + 2, l) &
                          + q_prim_vf(xibeg + 1)%sf(j, k + 3, l)) &
                         /(12d0*(y_cb(k) - y_cb(k - 1)))
            ! dxiz / dy
            tensora(8) = (-3d0*q_prim_vf(xiend)%sf(j, k - 1, l) &
                          - 10d0*q_prim_vf(xiend)%sf(j, k, l) &
                          + 18d0*q_prim_vf(xiend)%sf(j, k + 1, l) &
                          - 6d0*q_prim_vf(xiend)%sf(j, k + 2, l) &
                          + q_prim_vf(xiend)%sf(j, k + 3, l)) &
                         /(12d0*(y_cb(k) - y_cb(k - 1)))
        else if (k == iye - 1) then
            ! dxix / dy
            tensora(4) = (3d0*q_prim_vf(xibeg)%sf(j, k + 1, l) &
                          + 10d0*q_prim_vf(xibeg)%sf(j, k, l) &
                          - 18d0*q_prim_vf(xibeg)%sf(j, k - 1, l) &
                          + 6d0*q_prim_vf(xibeg)%sf(j, k - 2, l) &
                          - q_prim_vf(xibeg)%sf(j, k - 3, l)) &
                         /(12d0*(y_cb(k) - y_cb(k - 1)))
            ! dxiy / dy
            tensora(5) = (3d0*q_prim_vf(xibeg + 1)%sf(j, k + 1, l) &
                          + 10d0*q_prim_vf(xibeg + 1)%sf(j, k, l) &
                          - 18d0*q_prim_vf(xibeg + 1)%sf(j, k - 1, l) &
                          + 6d0*q_prim_vf(xibeg + 1)%sf(j, k - 2, l) &
                          - q_prim_vf(xibeg + 1)%sf(j, k - 3, l)) &
                         /(12d0*(y_cb(k) - y_cb(k - 1)))
            ! dxiz / dy
            tensora(8) = (3d0*q_prim_vf(xiend)%sf(j, k + 1, l) &
                          + 10d0*q_prim_vf(xiend)%sf(j, k, l) &
                          - 18d0*q_prim_vf(xiend)%sf(j, k - 1, l) &
                          + 6d0*q_prim_vf(xiend)%sf(j, k - 2, l) &
                          - q_prim_vf(xiend)%sf(j, k - 3, l)) &
                         /(12d0*(y_cb(k) - y_cb(k - 1)))
        else if (k == iye) then
            ! dxix / dy
            tensora(4) = (25d0*q_prim_vf(xibeg)%sf(j, k, l) &
                          - 48d0*q_prim_vf(xibeg)%sf(j, k - 1, l) &
                          + 36d0*q_prim_vf(xibeg)%sf(j, k - 2, l) &
                          - 16d0*q_prim_vf(xibeg)%sf(j, k - 3, l) &
                          + 3d0*q_prim_vf(xibeg)%sf(j, k - 4, l)) &
                         /(12d0*(y_cb(k) - y_cb(k - 1)))
            ! dxiy / dy
            tensora(5) = (25d0*q_prim_vf(xibeg + 1)%sf(j, k, l) &
                          - 48d0*q_prim_vf(xibeg + 1)%sf(j, k - 1, l) &
                          + 36d0*q_prim_vf(xibeg + 1)%sf(j, k - 2, l) &
                          - 16d0*q_prim_vf(xibeg + 1)%sf(j, k - 3, l) &
                          + 3d0*q_prim_vf(xibeg + 1)%sf(j, k - 4, l)) &
                         /(12d0*(y_cb(k) - y_cb(k - 1)))
            ! dxiz / dy
            tensora(8) = (25d0*q_prim_vf(xiend)%sf(j, k, l) &
                          - 48d0*q_prim_vf(xiend)%sf(j, k - 1, l) &
                          + 36d0*q_prim_vf(xiend)%sf(j, k - 2, l) &
                          - 16d0*q_prim_vf(xiend)%sf(j, k - 3, l) &
                          + 3d0*q_prim_vf(xiend)%sf(j, k - 4, l)) &
                         /(12d0*(y_cb(k) - y_cb(k - 1)))
        else
            ! dxix / dy
            tensora(4) = (q_prim_vf(xibeg)%sf(j, k - 2, l) &
                          - 8d0*q_prim_vf(xibeg)%sf(j, k - 1, l) &
                          + 8d0*q_prim_vf(xibeg)%sf(j, k + 1, l) &
                          - q_prim_vf(xibeg)%sf(j, k + 2, l)) &
                         /(12d0*(y_cb(k) - y_cb(k - 1)))
            ! dxiy / dy
            tensora(5) = (q_prim_vf(xibeg + 1)%sf(j, k - 2, l) &
                          - 8d0*q_prim_vf(xibeg + 1)%sf(j, k - 1, l) &
                          + 8d0*q_prim_vf(xibeg + 1)%sf(j, k + 1, l) &
                          - q_prim_vf(xibeg + 1)%sf(j, k + 2, l)) &
                         /(12d0*(y_cb(k) - y_cb(k - 1)))
            ! dxiz / dy
            tensora(8) = (q_prim_vf(xiend)%sf(j, k - 2, l) &
                          - 8d0*q_prim_vf(xiend)%sf(j, k - 1, l) &
                          + 8d0*q_prim_vf(xiend)%sf(j, k + 1, l) &
                          - q_prim_vf(xiend)%sf(j, k + 2, l)) &
                         /(12d0*(y_cb(k) - y_cb(k - 1)))
        end if

        ! 3D
        if (l == izb) then
            ! dxix / dz
            tensora(3) = (-25d0*q_prim_vf(xibeg)%sf(j, k, l) &
                          + 48d0*q_prim_vf(xibeg)%sf(j, k, l + 1) &
                          - 36d0*q_prim_vf(xibeg)%sf(j, k, l + 2) &
                          + 16d0*q_prim_vf(xibeg)%sf(j, k, l + 3) &
                          - 3d0*q_prim_vf(xibeg)%sf(j, k, l + 4)) &
                         /(12d0*(z_cb(l + 1) - z_cb(l)))
            ! dxiy / dz
            tensora(6) = (-25d0*q_prim_vf(xibeg + 1)%sf(j, k, l) &
                          + 48d0*q_prim_vf(xibeg + 1)%sf(j, k, l + 1) &
                          - 36d0*q_prim_vf(xibeg + 1)%sf(j, k, l + 2) &
                          + 16d0*q_prim_vf(xibeg + 1)%sf(j, k, l + 3) &
                          - 3d0*q_prim_vf(xibeg + 1)%sf(j, k, l + 4)) &
                         /(12d0*(z_cb(l + 1) - z_cb(l)))
            ! dxiz / dz
            tensora(9) = (-25d0*q_prim_vf(xiend)%sf(j, k, l) &
                          + 48d0*q_prim_vf(xiend)%sf(j, k, l + 1) &
                          - 36d0*q_prim_vf(xiend)%sf(j, k, l + 2) &
                          + 16d0*q_prim_vf(xiend)%sf(j, k, l + 3) &
                          - 3d0*q_prim_vf(xiend)%sf(j, k, l + 4)) &
                         /(12d0*(z_cb(l + 1) - z_cb(l)))
        else if (l == izb + 1) then
            ! dxix / dz
            tensora(3) = (-3d0*q_prim_vf(xibeg)%sf(j, k, l - 1) &
                          - 10d0*q_prim_vf(xibeg)%sf(j, k, l) &
                          + 18d0*q_prim_vf(xibeg)%sf(j, k, l + 1) &
                          - 6d0*q_prim_vf(xibeg)%sf(j, k, l + 2) &
                          + q_prim_vf(xibeg)%sf(j, k, l + 3)) &
                         /(12d0*(z_cb(l) - z_cb(l - 1)))
            ! dxiy / dz
            tensora(6) = (-3d0*q_prim_vf(xibeg + 1)%sf(j, k, l - 1) &
                          - 10d0*q_prim_vf(xibeg + 1)%sf(j, k, l) &
                          + 18d0*q_prim_vf(xibeg + 1)%sf(j, k, l + 1) &
                          - 6d0*q_prim_vf(xibeg + 1)%sf(j, k, l + 2) &
                          + q_prim_vf(xibeg + 1)%sf(j, k, l + 3)) &
                         /(12d0*(z_cb(l) - z_cb(l - 1)))
            ! dxiz / dz
            tensora(9) = (-3d0*q_prim_vf(xiend)%sf(j, k, l - 1) &
                          - 10d0*q_prim_vf(xiend)%sf(j, k, l) &
                          + 18d0*q_prim_vf(xiend)%sf(j, k, l + 1) &
                          - 6d0*q_prim_vf(xiend)%sf(j, k, l + 2) &
                          + q_prim_vf(xiend)%sf(j, k, l + 3)) &
                         /(12d0*(z_cb(l) - z_cb(l - 1)))
        else if (l == ize - 1) then
            ! dxix / dz
            tensora(3) = (3d0*q_prim_vf(xibeg)%sf(j, k, l + 1) &
                          + 10d0*q_prim_vf(xibeg)%sf(j, k, l) &
                          - 18d0*q_prim_vf(xibeg)%sf(j, k, l - 1) &
                          + 6d0*q_prim_vf(xibeg)%sf(j, k, l - 2) &
                          - q_prim_vf(xibeg)%sf(j, k, l - 3)) &
                         /(12d0*(z_cb(l) - z_cb(l - 1)))
            ! dxiy / dz
            tensora(6) = (3d0*q_prim_vf(xibeg + 1)%sf(j, k, l + 1) &
                          + 10d0*q_prim_vf(xibeg + 1)%sf(j, k, l) &
                          - 18d0*q_prim_vf(xibeg + 1)%sf(j, k, l - 1) &
                          + 6d0*q_prim_vf(xibeg + 1)%sf(j, k, l - 2) &
                          - q_prim_vf(xibeg + 1)%sf(j, k, l - 3)) &
                         /(12d0*(z_cb(l) - z_cb(l - 1)))
            ! dxiz / dz
            tensora(9) = (3d0*q_prim_vf(xiend)%sf(j, k, l + 1) &
                          + 10d0*q_prim_vf(xiend)%sf(j, k, l) &
                          - 18d0*q_prim_vf(xiend)%sf(j, k, l - 1) &
                          + 6d0*q_prim_vf(xiend)%sf(j, k, l - 2) &
                          - q_prim_vf(xiend)%sf(j, k, l - 3)) &
                         /(12d0*(z_cb(l) - z_cb(l - 1)))
        else if (l == ize) then
            ! dxix / dz
            tensora(3) = (25d0*q_prim_vf(xibeg)%sf(j, k, l) &
                          - 48d0*q_prim_vf(xibeg)%sf(j, k, l - 1) &
                          + 36d0*q_prim_vf(xibeg)%sf(j, k, l - 2) &
                          - 16d0*q_prim_vf(xibeg)%sf(j, k, l - 3) &
                          + 3d0*q_prim_vf(xibeg)%sf(j, k, l - 4)) &
                         /(12d0*(z_cb(l) - z_cb(l - 1)))
            ! dxiy / dz
            tensora(6) = (25d0*q_prim_vf(xibeg + 1)%sf(j, k, l) &
                          - 48d0*q_prim_vf(xibeg + 1)%sf(j, k, l - 1) &
                          + 36d0*q_prim_vf(xibeg + 1)%sf(j, k, l - 2) &
                          - 16d0*q_prim_vf(xibeg + 1)%sf(j, k, l - 3) &
                          + 3d0*q_prim_vf(xibeg + 1)%sf(j, k, l - 4)) &
                         /(12d0*(z_cb(l) - z_cb(l - 1)))
            ! dxiz / dz
            tensora(9) = (25d0*q_prim_vf(xiend)%sf(j, k, l) &
                          - 48d0*q_prim_vf(xiend)%sf(j, k, l - 1) &
                          + 36d0*q_prim_vf(xiend)%sf(j, k, l - 2) &
                          - 16d0*q_prim_vf(xiend)%sf(j, k, l - 3) &
                          + 3d0*q_prim_vf(xiend)%sf(j, k, l - 4)) &
                         /(12d0*(z_cb(l) - z_cb(l - 1)))
        else
            ! dxix / dz
            tensora(3) = (q_prim_vf(xibeg)%sf(j, k, l - 2) &
                          - 8d0*q_prim_vf(xibeg)%sf(j, k, l - 1) &
                          + 8d0*q_prim_vf(xibeg)%sf(j, k, l + 1) &
                          - q_prim_vf(xibeg)%sf(j, k, l + 2)) &
                         /(12d0*(z_cb(l) - z_cb(l - 1)))
            ! dxiy / dz
            tensora(6) = (q_prim_vf(xibeg + 1)%sf(j, k, l - 2) &
                          - 8d0*q_prim_vf(xibeg + 1)%sf(j, k, l - 1) &
                          + 8d0*q_prim_vf(xibeg + 1)%sf(j, k, l + 1) &
                          - q_prim_vf(xibeg + 1)%sf(j, k, l + 2)) &
                         /(12d0*(z_cb(l) - z_cb(l - 1)))
            ! dxiz / dz
            tensora(9) = (q_prim_vf(xiend)%sf(j, k, l - 2) &
                          - 8d0*q_prim_vf(xiend)%sf(j, k, l - 1) &
                          + 8d0*q_prim_vf(xiend)%sf(j, k, l + 1) &
                          - q_prim_vf(xiend)%sf(j, k, l + 2)) &
                         /(12d0*(z_cb(l) - z_cb(l - 1)))
        end if

        ! STEP 2a: computing the adjoint of the grad_xi tensor for the inverse
        tensorb(1) = tensora(5)*tensora(9) - tensora(6)*tensora(8)
        tensorb(2) = -(tensora(2)*tensora(9) - tensora(3)*tensora(8))
        tensorb(3) = tensora(2)*tensora(6) - tensora(3)*tensora(5)
        tensorb(4) = -(tensora(4)*tensora(9) - tensora(6)*tensora(7))
        tensorb(5) = tensora(1)*tensora(9) - tensora(3)*tensora(7)
        tensorb(6) = -(tensora(1)*tensora(6) - tensora(4)*tensora(3))
        tensorb(7) = tensora(4)*tensora(8) - tensora(5)*tensora(7)
        tensorb(8) = -(tensora(1)*tensora(8) - tensora(2)*tensora(7))
        tensorb(9) = tensora(1)*tensora(5) - tensora(2)*tensora(4)

        ! STEP 2b: computing the determinant of the grad_xi tensor
        tensorb(tensor_size) = tensora(1)*(tensora(5)*tensora(9) - tensora(6)*tensora(8)) &
                               - tensora(2)*(tensora(4)*tensora(9) - tensora(6)*tensora(7)) &
                               + tensora(3)*(tensora(4)*tensora(8) - tensora(5)*tensora(7))

        ! STEP 2c: computing the inverse of grad_xi tensor = F
        ! tensorb is the adjoint, tensora becomes the inverse
        ! STEP 4: store the determinant of F in the last entry of the tensor

        !$acc loop seq
        do i = 1, tensor_size - 1
            tensora(i) = tensorb(i)/tensorb(tensor_size)
        end do

        ! STEP 3: computing F tranpose F
        tensorb(1) = tensora(1)**2
        tensorb(1) = tensorb(1) + tensora(4)**2 + tensora(7)**2
        tensorb(5) = tensora(2) + tensora(5)**2 + tensora(8)**2
        tensorb(9) = tensora(3) + tensora(6)**2 + tensora(9)**2
        tensorb(2) = tensora(1)*tensora(2) + tensora(4)*tensora(5) + tensora(7)*tensora(8)
        tensorb(3) = tensora(1)*tensora(3) + tensora(4)*tensora(6) + tensora(7)*tensora(9)
        tensorb(6) = tensora(2)*tensora(3) + tensora(5)*tensora(6) + tensora(8)*tensora(9)
        tensorb(4) = tensorb(2)
        tensorb(7) = tensorb(3)
        tensorb(8) = tensorb(4)

    end subroutine s_compute_gradient_xi3d_acc

    !>  The following subroutine handles the calculation of the btensor.
        !!   The calculation of the btensor takes qprimvf.
        !! @param q_prim_vf Primitive variables
        !! @param btensor is the output
        !! calculate the grad_xi, grad_xi is a nxn tensor
        !! calculate the inverse of grad_xi to obtain F, F is a nxn tensor
        !! calculate the FFtranspose to obtain the btensor, btensor is nxn tensor
        !! btensor is symmetric, save the data space
        !! neo-Hookean only at this time, will need to be changed later
    function f_elastic_energy(btensor, j, k, l)
#ifdef MFC_SIMULATION
        !$acc routine seq
#endif
        type(scalar_field), dimension(b_size), intent(IN) :: btensor
        integer, intent(IN) :: j, k, l
        real(kind(0d0)) :: invariant1, f_elastic_energy

        f_elastic_energy = 0d0
        invariant1 = btensor(1)%sf(j, k, l)
        !if (num_dims == 2) then
        !    invariant1 = invariant1 + btensor(3)%sf(j, k, l)
        !elseif (num_dims == 3) then
        invariant1 = invariant1 + btensor(4)%sf(j, k, l) + btensor(6)%sf(j, k, l)
        !end if

        ! compute the invariant without the elastic modulus
        f_elastic_energy = 0.5d0*(invariant1 - 3.0d0)/btensor(b_size)%sf(j, k, l)

    end function f_elastic_energy

end module m_xi_tensor_calc

