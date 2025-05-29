!>
!! @file m_mhd.f90
!! @brief Contains module m_mhd

#:include 'macros.fpp'

!> @brief This module is used to compute source terms for magnetohydrodynamics
!!              Note: not applicable for 1D

module m_mhd

    use m_derived_types        !< Definitions of the derived types
    use m_global_parameters    !< Definitions of the global parameters
    use m_finite_differences
    use m_helper
    use m_mpi_common

    implicit none

    private; public :: s_initialize_mhd_powell_module, &
 s_finalize_mhd_powell_module, &
 s_compute_mhd_powell_rhs

    real(wp), allocatable, dimension(:, :, :) :: du_dx, du_dy, du_dz
    real(wp), allocatable, dimension(:, :, :) :: dv_dx, dv_dy, dv_dz
    real(wp), allocatable, dimension(:, :, :) :: dw_dx, dw_dy, dw_dz
    !$acc declare create(du_dx,du_dy,du_dz,dv_dx,dv_dy,dv_dz,dw_dx,dw_dy,dw_dz)

    real(wp), allocatable, dimension(:, :) :: fd_coeff_x_h
    real(wp), allocatable, dimension(:, :) :: fd_coeff_y_h
    real(wp), allocatable, dimension(:, :) :: fd_coeff_z_h
    !$acc declare create(fd_coeff_x_h,fd_coeff_y_h,fd_coeff_z_h)

contains

    subroutine s_initialize_mhd_powell_module

        ! Additional safety check beyond m_checker
        if (n == 0) call s_mpi_abort('Fatal Error: Powell correction is not applicable for 1D')

        @:ALLOCATE(du_dx(0:m,0:n,0:p), dv_dx(0:m,0:n,0:p), dw_dx(0:m,0:n,0:p))
        @:ALLOCATE(du_dy(0:m,0:n,0:p), dv_dy(0:m,0:n,0:p), dw_dy(0:m,0:n,0:p))
        if (p > 0) then
            @:ALLOCATE(dw_dx(0:m,0:n,0:p), dw_dy(0:m,0:n,0:p), dw_dz(0:m,0:n,0:p))
        end if

        @:ALLOCATE(fd_coeff_x_h(-fd_number:fd_number, 0:m))
        @:ALLOCATE(fd_coeff_y_h(-fd_number:fd_number, 0:n))
        if (p > 0) then
            @:ALLOCATE(fd_coeff_z_h(-fd_number:fd_number, 0:p))
        end if

        ! Computing centered finite difference coefficients
        call s_compute_finite_difference_coefficients(m, x_cc, fd_coeff_x_h, buff_size, fd_number, fd_order)
        !$acc update device(fd_coeff_x_h)
        call s_compute_finite_difference_coefficients(n, y_cc, fd_coeff_y_h, buff_size, fd_number, fd_order)
        !$acc update device(fd_coeff_y_h)
        if (p > 0) then
            call s_compute_finite_difference_coefficients(p, z_cc, fd_coeff_z_h, buff_size, fd_number, fd_order)
            !$acc update device(fd_coeff_z_h)
        end if

    end subroutine s_initialize_mhd_powell_module

    !>  Compute the Powell source term to correct the magnetic field divergence.
        !!      The Powell source term is:
        !!      S = - (divB) [ 0, Bx, By, Bz, vdotB, vx, vy, vz ]^T
        !!  @param q_prim_vf  Primitive variables
        !!  @param rhs_vf     rhs variables
    subroutine s_compute_mhd_powell_rhs(q_prim_vf, rhs_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf

        integer :: k, l, q, r
        real(wp), dimension(3) :: v, B
        real(wp) :: divB, vdotB

        !$acc parallel loop collapse(3) gang vector default(present) &
        !$acc private(v, B)
        do q = 0, p
            do l = 0, n
                do k = 0, m

                    divB = 0._wp
                    !$acc loop seq
                    do r = -fd_number, fd_number
                        divB = divB + q_prim_vf(B_idx%beg)%sf(k + r, l, q)*fd_coeff_x_h(r, k)
                    end do
                    !$acc loop seq
                    do r = -fd_number, fd_number
                        divB = divB + q_prim_vf(B_idx%beg + 1)%sf(k, l + r, q)*fd_coeff_y_h(r, l)
                    end do
                    if (p > 0) then
                        !$acc loop seq
                        do r = -fd_number, fd_number
                            divB = divB + q_prim_vf(B_idx%beg + 2)%sf(k, l, q + r)*fd_coeff_z_h(r, q)
                        end do
                    end if

                    v(1) = q_prim_vf(momxb)%sf(k, l, q)
                    v(2) = q_prim_vf(momxb + 1)%sf(k, l, q)
                    v(3) = q_prim_vf(momxb + 2)%sf(k, l, q)

                    B(1) = q_prim_vf(B_idx%beg)%sf(k, l, q)
                    B(2) = q_prim_vf(B_idx%beg + 1)%sf(k, l, q)
                    B(3) = q_prim_vf(B_idx%beg + 2)%sf(k, l, q)

                    vdotB = sum(v*B)

                    ! 1: rho -> unchanged
                    ! 2: vx  -> - (divB) * Bx
                    ! 3: vy  -> - (divB) * By
                    ! 4: vz  -> - (divB) * Bz
                    ! 5: E   -> - (divB) * (vdotB)
                    ! 6: Bx  -> - (divB) * vx
                    ! 7: By  -> - (divB) * vy
                    ! 8: Bz  -> - (divB) * vz

                    rhs_vf(momxb)%sf(k, l, q) = rhs_vf(momxb)%sf(k, l, q) - divB*B(1)
                    rhs_vf(momxb + 1)%sf(k, l, q) = rhs_vf(momxb + 1)%sf(k, l, q) - divB*B(2)
                    rhs_vf(momxb + 2)%sf(k, l, q) = rhs_vf(momxb + 2)%sf(k, l, q) - divB*B(3)

                    rhs_vf(E_idx)%sf(k, l, q) = rhs_vf(E_idx)%sf(k, l, q) - divB*vdotB

                    rhs_vf(B_idx%beg)%sf(k, l, q) = rhs_vf(B_idx%beg)%sf(k, l, q) - divB*v(1)
                    rhs_vf(B_idx%beg + 1)%sf(k, l, q) = rhs_vf(B_idx%beg + 1)%sf(k, l, q) - divB*v(2)
                    rhs_vf(B_idx%beg + 2)%sf(k, l, q) = rhs_vf(B_idx%beg + 2)%sf(k, l, q) - divB*v(3)

                end do
            end do
        end do
        !$acc end parallel loop

    end subroutine s_compute_mhd_powell_rhs

    subroutine s_finalize_mhd_powell_module

        @:DEALLOCATE(du_dx, dv_dx, dw_dx)
        @:DEALLOCATE(fd_coeff_x_h)
        @:DEALLOCATE(du_dy, dv_dy, dw_dy)
        @:DEALLOCATE(fd_coeff_y_h)
        if (p > 0) then
            @:DEALLOCATE(dw_dx, dw_dy, dw_dz)
            @:DEALLOCATE(fd_coeff_z_h)
        end if

    end subroutine s_finalize_mhd_powell_module

end module m_mhd
