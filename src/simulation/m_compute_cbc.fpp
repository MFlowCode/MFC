!>
!! @file m_compute_cbc.f90
!! @brief CBC computation module

#:include 'macros.fpp'

module m_compute_cbc
    use m_global_parameters
    implicit none

    private; public :: s_compute_slip_wall_L, &
 s_compute_nonreflecting_subsonic_buffer_L, &
 s_compute_nonreflecting_subsonic_inflow_L, &
 s_compute_nonreflecting_subsonic_outflow_L, &
 s_compute_force_free_subsonic_outflow_L, &
 s_compute_constant_pressure_subsonic_outflow_L, &
 s_compute_supersonic_inflow_L, &
 s_compute_supersonic_outflow_L

contains
    !> Base L1 calculation
    function f_base_L1(lambda, rho, c, dpres_ds, dvel_ds) result(L1)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), dimension(3), intent(in) :: lambda
        real(wp), intent(in) :: rho, c, dpres_ds
        real(wp), dimension(num_dims), intent(in) :: dvel_ds
        real(wp) :: L1
        L1 = lambda(1)*(dpres_ds - rho*c*dvel_ds(dir_idx(1)))
    end function f_base_L1

    !> Fill density L variables
    subroutine s_fill_density_L(L, lambda_factor, lambda2, c, mf, dalpha_rho_ds, dpres_ds)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), dimension(sys_size), intent(inout) :: L
        real(wp), intent(in) :: lambda_factor, lambda2, c
        real(wp), dimension(num_fluids), intent(in) :: mf, dalpha_rho_ds
        real(wp), intent(in) :: dpres_ds
        integer :: i

        ! $:GPU_LOOP(parallelism='[seq]')
        do i = 2, momxb
            L(i) = lambda_factor*lambda2*(c*c*dalpha_rho_ds(i - 1) - mf(i - 1)*dpres_ds)
        end do
    end subroutine s_fill_density_L

    !> Fill velocity L variables
    subroutine s_fill_velocity_L(L, lambda_factor, lambda2, dvel_ds)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), dimension(sys_size), intent(inout) :: L
        real(wp), intent(in) :: lambda_factor, lambda2
        real(wp), dimension(num_dims), intent(in) :: dvel_ds
        integer :: i

        ! $:GPU_LOOP(parallelism='[seq]')
        do i = momxb + 1, momxe
            L(i) = lambda_factor*lambda2*dvel_ds(dir_idx(i - contxe))
        end do
    end subroutine s_fill_velocity_L

    !> Fill advection L variables
    subroutine s_fill_advection_L(L, lambda_factor, lambda2, dadv_ds)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), dimension(sys_size), intent(inout) :: L
        real(wp), intent(in) :: lambda_factor, lambda2
        real(wp), dimension(num_fluids), intent(in) :: dadv_ds
        integer :: i

        ! $:GPU_LOOP(parallelism='[seq]')
        do i = E_idx, advxe - 1
            L(i) = lambda_factor*lambda2*dadv_ds(i - momxe)
        end do
    end subroutine s_fill_advection_L

    !> Fill chemistry L variables
    subroutine s_fill_chemistry_L(L, lambda_factor, lambda2, dYs_ds)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), dimension(sys_size), intent(inout) :: L
        real(wp), intent(in) :: lambda_factor, lambda2
        real(wp), dimension(num_species), intent(in) :: dYs_ds
        integer :: i

        if (.not. chemistry) return

        ! $:GPU_LOOP(parallelism='[seq]')
        do i = chemxb, chemxe
            L(i) = lambda_factor*lambda2*dYs_ds(i - chemxb + 1)
        end do
    end subroutine s_fill_chemistry_L

    !> Slip wall CBC (Thompson 1990, pg. 451)
    subroutine s_compute_slip_wall_L(lambda, L, rho, c, dpres_ds, dvel_ds)
        $:GPU_ROUTINE(function_name='s_compute_slip_wall_L',parallelism='[seq]', &
            & cray_inline=True)

        real(wp), dimension(3), intent(in) :: lambda
        real(wp), dimension(sys_size), intent(inout) :: L
        real(wp), intent(in) :: rho, c, dpres_ds
        real(wp), dimension(num_dims), intent(in) :: dvel_ds
        integer :: i

        L(1) = f_base_L1(lambda, rho, c, dpres_ds, dvel_ds)
        L(2:advxe - 1) = 0._wp
        L(advxe) = L(1)
    end subroutine s_compute_slip_wall_L

    !> Nonreflecting subsonic buffer CBC (Thompson 1987, pg. 13)
    subroutine s_compute_nonreflecting_subsonic_buffer_L(lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds, dYs_ds)
        $:GPU_ROUTINE(function_name='s_compute_nonreflecting_subsonic_buffer_L', &
            & parallelism='[seq]', cray_inline=True)

        real(wp), dimension(3), intent(in) :: lambda
        real(wp), dimension(sys_size), intent(inout) :: L
        real(wp), intent(in) :: rho, c
        real(wp), dimension(num_fluids), intent(in) :: mf, dalpha_rho_ds
        real(wp), intent(in) :: dpres_ds
        real(wp), dimension(num_dims), intent(in) :: dvel_ds
        real(wp), dimension(num_fluids), intent(in) :: dadv_ds
        real(wp), dimension(num_species), intent(in) :: dYs_ds
        real(wp) :: lambda_factor

        lambda_factor = (5.e-1_wp - 5.e-1_wp*sign(1._wp, lambda(1)))
        L(1) = lambda_factor*lambda(1)*(dpres_ds - rho*c*dvel_ds(dir_idx(1)))

        lambda_factor = (5.e-1_wp - 5.e-1_wp*sign(1._wp, lambda(2)))
        call s_fill_density_L(L, lambda_factor, lambda(2), c, mf, dalpha_rho_ds, dpres_ds)
        call s_fill_velocity_L(L, lambda_factor, lambda(2), dvel_ds)
        call s_fill_advection_L(L, lambda_factor, lambda(2), dadv_ds)
        call s_fill_chemistry_L(L, lambda_factor, lambda(2), dYs_ds)

        lambda_factor = (5.e-1_wp - 5.e-1_wp*sign(1._wp, lambda(3)))
        L(advxe) = lambda_factor*lambda(3)*(dpres_ds + rho*c*dvel_ds(dir_idx(1)))
    end subroutine s_compute_nonreflecting_subsonic_buffer_L

    !> Nonreflecting subsonic inflow CBC (Thompson 1990, pg. 455)
    subroutine s_compute_nonreflecting_subsonic_inflow_L(lambda, L, rho, c, dpres_ds, dvel_ds)
        $:GPU_ROUTINE(function_name='s_compute_nonreflecting_subsonic_inflow_L', &
            & parallelism='[seq]', cray_inline=True)

        real(wp), dimension(3), intent(in) :: lambda
        real(wp), dimension(sys_size), intent(inout) :: L
        real(wp), intent(in) :: rho, c, dpres_ds
        real(wp), dimension(num_dims), intent(in) :: dvel_ds

        L(1) = f_base_L1(lambda, rho, c, dpres_ds, dvel_ds)
        L(2:advxe) = 0._wp
        if (chemistry) L(chemxb:chemxe) = 0._wp
    end subroutine s_compute_nonreflecting_subsonic_inflow_L

    !> Nonreflecting subsonic outflow CBC (Thompson 1990, pg. 454)
    subroutine s_compute_nonreflecting_subsonic_outflow_L(lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds, dYs_ds)
        $:GPU_ROUTINE(function_name='s_compute_nonreflecting_subsonic_outflow_L', &
            & parallelism='[seq]', cray_inline=True)

        real(wp), dimension(3), intent(in) :: lambda
        real(wp), dimension(sys_size), intent(inout) :: L
        real(wp), intent(in) :: rho, c
        real(wp), dimension(num_fluids), intent(in) :: mf, dalpha_rho_ds
        real(wp), intent(in) :: dpres_ds
        real(wp), dimension(num_dims), intent(in) :: dvel_ds
        real(wp), dimension(num_fluids), intent(in) :: dadv_ds
        real(wp), dimension(num_species), intent(in) :: dYs_ds

        L(1) = f_base_L1(lambda, rho, c, dpres_ds, dvel_ds)
        call s_fill_density_L(L, 1._wp, lambda(2), c, mf, dalpha_rho_ds, dpres_ds)
        call s_fill_velocity_L(L, 1._wp, lambda(2), dvel_ds)
        call s_fill_advection_L(L, 1._wp, lambda(2), dadv_ds)
        call s_fill_chemistry_L(L, 1._wp, lambda(2), dYs_ds)
        L(advxe) = 0._wp
    end subroutine s_compute_nonreflecting_subsonic_outflow_L

    !> Force-free subsonic outflow CBC (Thompson 1990, pg. 454)
    subroutine s_compute_force_free_subsonic_outflow_L(lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds)
        $:GPU_ROUTINE(function_name='s_compute_force_free_subsonic_outflow_L', &
            & parallelism='[seq]', cray_inline=True)

        real(wp), dimension(3), intent(in) :: lambda
        real(wp), dimension(sys_size), intent(inout) :: L
        real(wp), intent(in) :: rho, c
        real(wp), dimension(num_fluids), intent(in) :: mf, dalpha_rho_ds
        real(wp), intent(in) :: dpres_ds
        real(wp), dimension(num_dims), intent(in) :: dvel_ds
        real(wp), dimension(num_fluids), intent(in) :: dadv_ds

        L(1) = f_base_L1(lambda, rho, c, dpres_ds, dvel_ds)
        call s_fill_density_L(L, 1._wp, lambda(2), c, mf, dalpha_rho_ds, dpres_ds)
        call s_fill_velocity_L(L, 1._wp, lambda(2), dvel_ds)
        call s_fill_advection_L(L, 1._wp, lambda(2), dadv_ds)
        L(advxe) = L(1) + 2._wp*rho*c*lambda(2)*dvel_ds(dir_idx(1))
    end subroutine s_compute_force_free_subsonic_outflow_L

    !> Constant pressure subsonic outflow CBC (Thompson 1990, pg. 455)
    subroutine s_compute_constant_pressure_subsonic_outflow_L(lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds)
        $:GPU_ROUTINE(function_name='s_compute_constant_pressure_subsonic_outflow_L', &
            & parallelism='[seq]', cray_inline=True)

        real(wp), dimension(3), intent(in) :: lambda
        real(wp), dimension(sys_size), intent(inout) :: L
        real(wp), intent(in) :: rho, c
        real(wp), dimension(num_fluids), intent(in) :: mf, dalpha_rho_ds
        real(wp), intent(in) :: dpres_ds
        real(wp), dimension(num_dims), intent(in) :: dvel_ds
        real(wp), dimension(num_fluids), intent(in) :: dadv_ds

        L(1) = f_base_L1(lambda, rho, c, dpres_ds, dvel_ds)
        call s_fill_density_L(L, 1._wp, lambda(2), c, mf, dalpha_rho_ds, dpres_ds)
        call s_fill_velocity_L(L, 1._wp, lambda(2), dvel_ds)
        call s_fill_advection_L(L, 1._wp, lambda(2), dadv_ds)
        L(advxe) = -L(1)
    end subroutine s_compute_constant_pressure_subsonic_outflow_L

    !> Supersonic inflow CBC (Thompson 1990, pg. 453)
    subroutine s_compute_supersonic_inflow_L(L)
        $:GPU_ROUTINE(function_name='s_compute_supersonic_inflow_L', &
            & parallelism='[seq]', cray_inline=True)

        real(wp), dimension(sys_size), intent(inout) :: L
        L(1:advxe) = 0._wp
        if (chemistry) L(chemxb:chemxe) = 0._wp
    end subroutine s_compute_supersonic_inflow_L

    !> Supersonic outflow CBC (Thompson 1990, pg. 453)
    subroutine s_compute_supersonic_outflow_L(lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds, dYs_ds)
        $:GPU_ROUTINE(function_name='s_compute_supersonic_outflow_L', &
            & parallelism='[seq]', cray_inline=True)

        real(wp), dimension(3), intent(in) :: lambda
        real(wp), dimension(sys_size), intent(inout) :: L
        real(wp), intent(in) :: rho, c
        real(wp), dimension(num_fluids), intent(in) :: mf, dalpha_rho_ds
        real(wp), intent(in) :: dpres_ds
        real(wp), dimension(num_dims), intent(in) :: dvel_ds
        real(wp), dimension(num_fluids), intent(in) :: dadv_ds
        real(wp), dimension(num_species), intent(in) :: dYs_ds

        L(1) = f_base_L1(lambda, rho, c, dpres_ds, dvel_ds)
        call s_fill_density_L(L, 1._wp, lambda(2), c, mf, dalpha_rho_ds, dpres_ds)
        call s_fill_velocity_L(L, 1._wp, lambda(2), dvel_ds)
        call s_fill_advection_L(L, 1._wp, lambda(2), dadv_ds)
        call s_fill_chemistry_L(L, 1._wp, lambda(2), dYs_ds)
        L(advxe) = lambda(3)*(dpres_ds + rho*c*dvel_ds(dir_idx(1)))
    end subroutine s_compute_supersonic_outflow_L
end module m_compute_cbc
