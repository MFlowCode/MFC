!>
!! @file
!! @brief CBC computation module
#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief Characteristic boundary condition (CBC) computations for subsonic inflow, outflow, and slip walls
module m_compute_cbc

    use m_global_parameters

    implicit none

    private; public :: s_compute_slip_wall_L, s_compute_nonreflecting_subsonic_buffer_L, &
        & s_compute_nonreflecting_subsonic_inflow_L, s_compute_nonreflecting_subsonic_outflow_L, &
        & s_compute_force_free_subsonic_outflow_L, s_compute_constant_pressure_subsonic_outflow_L, s_compute_supersonic_inflow_L, &
        & s_compute_supersonic_outflow_L

contains
    !> Base L1 calculation
    function f_base_L1(lambda, rho, c, dpres_ds, dvel_ds) result(L1)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), dimension(3), intent(in) :: lambda
        real(wp), intent(in)               :: rho, c, dpres_ds
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(in) :: dvel_ds
        #:else
            real(wp), dimension(num_dims), intent(in) :: dvel_ds
        #:endif
        real(wp) :: L1
        L1 = lambda(1)*(dpres_ds - rho*c*dvel_ds(dir_idx(1)))

    end function f_base_L1

    !> Fill density L variables
    subroutine s_fill_density_L(L, lambda_factor, lambda2, c, mf, dalpha_rho_ds, dpres_ds)

        $:GPU_ROUTINE(parallelism='[seq]')
        #:if USING_AMD
            real(wp), dimension(20), intent(inout) :: L
        #:else
            real(wp), dimension(sys_size), intent(inout) :: L
        #:endif
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(in) :: mf, dalpha_rho_ds
        #:else
            real(wp), dimension(num_fluids), intent(in) :: mf, dalpha_rho_ds
        #:endif
        real(wp), intent(in) :: lambda_factor, lambda2, c
        real(wp), intent(in) :: dpres_ds
        integer              :: i

        ! $:GPU_LOOP(parallelism='[seq]')
        do i = 2, eqn_idx%mom%beg
            L(i) = lambda_factor*lambda2*(c*c*dalpha_rho_ds(i - 1) - mf(i - 1)*dpres_ds)
        end do

    end subroutine s_fill_density_L

    !> Fill velocity L variables
    subroutine s_fill_velocity_L(L, lambda_factor, lambda2, dvel_ds)

        $:GPU_ROUTINE(parallelism='[seq]')
        #:if USING_AMD
            real(wp), dimension(20), intent(inout) :: L
        #:else
            real(wp), dimension(sys_size), intent(inout) :: L
        #:endif
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(in) :: dvel_ds
        #:else
            real(wp), dimension(num_dims), intent(in) :: dvel_ds
        #:endif
        real(wp), intent(in) :: lambda_factor, lambda2
        integer              :: i

        ! $:GPU_LOOP(parallelism='[seq]')
        do i = eqn_idx%mom%beg + 1, eqn_idx%mom%end
            L(i) = lambda_factor*lambda2*dvel_ds(dir_idx(i - eqn_idx%cont%end))
        end do

    end subroutine s_fill_velocity_L

    !> Fill advection L variables
    subroutine s_fill_advection_L(L, lambda_factor, lambda2, dadv_ds)

        $:GPU_ROUTINE(parallelism='[seq]')
        #:if USING_AMD
            real(wp), dimension(20), intent(inout) :: L
        #:else
            real(wp), dimension(sys_size), intent(inout) :: L
        #:endif
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(in) :: dadv_ds
        #:else
            real(wp), dimension(num_fluids), intent(in) :: dadv_ds
        #:endif
        real(wp), intent(in) :: lambda_factor, lambda2
        integer              :: i

        ! $:GPU_LOOP(parallelism='[seq]')
        do i = eqn_idx%E, eqn_idx%adv%end - 1
            L(i) = lambda_factor*lambda2*dadv_ds(i - eqn_idx%mom%end)
        end do

    end subroutine s_fill_advection_L

    !> Fill chemistry L variables
    subroutine s_fill_chemistry_L(L, lambda_factor, lambda2, dYs_ds)

        $:GPU_ROUTINE(parallelism='[seq]')
        #:if USING_AMD
            real(wp), dimension(20), intent(inout) :: L
        #:else
            real(wp), dimension(sys_size), intent(inout) :: L
        #:endif
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(10), intent(in) :: dYs_ds
        #:else
            real(wp), dimension(num_species), intent(in) :: dYs_ds
        #:endif
        real(wp), intent(in) :: lambda_factor, lambda2
        integer              :: i

        if (.not. chemistry) return

        ! $:GPU_LOOP(parallelism='[seq]')
        do i = eqn_idx%species%beg, eqn_idx%species%end
            L(i) = lambda_factor*lambda2*dYs_ds(i - eqn_idx%species%beg + 1)
        end do

    end subroutine s_fill_chemistry_L

    !> Slip wall CBC (Thompson 1990, pg. 451)
    subroutine s_compute_slip_wall_L(lambda, L, rho, c, dpres_ds, dvel_ds)

        $:GPU_ROUTINE(function_name='s_compute_slip_wall_L',parallelism='[seq]', cray_inline=True)

        real(wp), dimension(3), intent(in) :: lambda
        #:if USING_AMD
            real(wp), dimension(20), intent(inout) :: L
        #:else
            real(wp), dimension(sys_size), intent(inout) :: L
        #:endif
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(in) :: dvel_ds
        #:else
            real(wp), dimension(num_dims), intent(in) :: dvel_ds
        #:endif
        real(wp), intent(in) :: rho, c, dpres_ds

        L(1) = f_base_L1(lambda, rho, c, dpres_ds, dvel_ds)
        L(2:eqn_idx%adv%end - 1) = 0._wp
        L(eqn_idx%adv%end) = L(1)

    end subroutine s_compute_slip_wall_L

    !> Nonreflecting subsonic buffer CBC (Thompson 1987, pg. 13)
    subroutine s_compute_nonreflecting_subsonic_buffer_L(lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds, dYs_ds)

        $:GPU_ROUTINE(function_name='s_compute_nonreflecting_subsonic_buffer_L', parallelism='[seq]', cray_inline=True)

        real(wp), dimension(3), intent(in) :: lambda
        #:if USING_AMD
            real(wp), dimension(20), intent(inout) :: L
        #:else
            real(wp), dimension(sys_size), intent(inout) :: L
        #:endif
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(in)  :: mf, dalpha_rho_ds
            real(wp), dimension(3), intent(in)  :: dvel_ds
            real(wp), dimension(3), intent(in)  :: dadv_ds
            real(wp), dimension(10), intent(in) :: dYs_ds
        #:else
            real(wp), dimension(num_fluids), intent(in)  :: mf, dalpha_rho_ds
            real(wp), dimension(num_dims), intent(in)    :: dvel_ds
            real(wp), dimension(num_fluids), intent(in)  :: dadv_ds
            real(wp), dimension(num_species), intent(in) :: dYs_ds
        #:endif
        real(wp), intent(in) :: rho, c
        real(wp), intent(in) :: dpres_ds
        real(wp)             :: lambda_factor

        lambda_factor = (5.e-1_wp - 5.e-1_wp*sign(1._wp, lambda(1)))
        L(1) = lambda_factor*lambda(1)*(dpres_ds - rho*c*dvel_ds(dir_idx(1)))

        lambda_factor = (5.e-1_wp - 5.e-1_wp*sign(1._wp, lambda(2)))
        call s_fill_density_L(L, lambda_factor, lambda(2), c, mf, dalpha_rho_ds, dpres_ds)
        call s_fill_velocity_L(L, lambda_factor, lambda(2), dvel_ds)
        call s_fill_advection_L(L, lambda_factor, lambda(2), dadv_ds)
        call s_fill_chemistry_L(L, lambda_factor, lambda(2), dYs_ds)

        lambda_factor = (5.e-1_wp - 5.e-1_wp*sign(1._wp, lambda(3)))
        L(eqn_idx%adv%end) = lambda_factor*lambda(3)*(dpres_ds + rho*c*dvel_ds(dir_idx(1)))

    end subroutine s_compute_nonreflecting_subsonic_buffer_L

    !> Nonreflecting subsonic inflow CBC (Thompson 1990, pg. 455)
    subroutine s_compute_nonreflecting_subsonic_inflow_L(lambda, L, rho, c, dpres_ds, dvel_ds)

        $:GPU_ROUTINE(function_name='s_compute_nonreflecting_subsonic_inflow_L', parallelism='[seq]', cray_inline=True)

        real(wp), dimension(3), intent(in) :: lambda
        #:if USING_AMD
            real(wp), dimension(20), intent(inout) :: L
        #:else
            real(wp), dimension(sys_size), intent(inout) :: L
        #:endif
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(in) :: dvel_ds
        #:else
            real(wp), dimension(num_dims), intent(in) :: dvel_ds
        #:endif
        real(wp), intent(in) :: rho, c, dpres_ds

        L(1) = f_base_L1(lambda, rho, c, dpres_ds, dvel_ds)
        L(2:eqn_idx%adv%end) = 0._wp
        if (chemistry) L(eqn_idx%species%beg:eqn_idx%species%end) = 0._wp

    end subroutine s_compute_nonreflecting_subsonic_inflow_L

    !> Nonreflecting subsonic outflow CBC (Thompson 1990, pg. 454)
    subroutine s_compute_nonreflecting_subsonic_outflow_L(lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds, dYs_ds)

        $:GPU_ROUTINE(function_name='s_compute_nonreflecting_subsonic_outflow_L', parallelism='[seq]', cray_inline=True)

        real(wp), dimension(3), intent(in) :: lambda
        #:if USING_AMD
            real(wp), dimension(20), intent(inout) :: L
        #:else
            real(wp), dimension(sys_size), intent(inout) :: L
        #:endif
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(in)  :: mf, dalpha_rho_ds
            real(wp), dimension(3), intent(in)  :: dvel_ds
            real(wp), dimension(3), intent(in)  :: dadv_ds
            real(wp), dimension(10), intent(in) :: dYs_ds
        #:else
            real(wp), dimension(num_fluids), intent(in)  :: mf, dalpha_rho_ds
            real(wp), dimension(num_dims), intent(in)    :: dvel_ds
            real(wp), dimension(num_fluids), intent(in)  :: dadv_ds
            real(wp), dimension(num_species), intent(in) :: dYs_ds
        #:endif
        real(wp), intent(in) :: rho, c
        real(wp), intent(in) :: dpres_ds

        L(1) = f_base_L1(lambda, rho, c, dpres_ds, dvel_ds)
        call s_fill_density_L(L, 1._wp, lambda(2), c, mf, dalpha_rho_ds, dpres_ds)
        call s_fill_velocity_L(L, 1._wp, lambda(2), dvel_ds)
        call s_fill_advection_L(L, 1._wp, lambda(2), dadv_ds)
        call s_fill_chemistry_L(L, 1._wp, lambda(2), dYs_ds)
        L(eqn_idx%adv%end) = 0._wp

    end subroutine s_compute_nonreflecting_subsonic_outflow_L

    !> Force-free subsonic outflow CBC (Thompson 1990, pg. 454)
    subroutine s_compute_force_free_subsonic_outflow_L(lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds)

        $:GPU_ROUTINE(function_name='s_compute_force_free_subsonic_outflow_L', parallelism='[seq]', cray_inline=True)

        real(wp), dimension(3), intent(in) :: lambda
        #:if USING_AMD
            real(wp), dimension(20), intent(inout) :: L
        #:else
            real(wp), dimension(sys_size), intent(inout) :: L
        #:endif
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(in) :: mf, dalpha_rho_ds
            real(wp), dimension(3), intent(in) :: dvel_ds
            real(wp), dimension(3), intent(in) :: dadv_ds
        #:else
            real(wp), dimension(num_fluids), intent(in) :: mf, dalpha_rho_ds
            real(wp), dimension(num_dims), intent(in)   :: dvel_ds
            real(wp), dimension(num_fluids), intent(in) :: dadv_ds
        #:endif
        real(wp), intent(in) :: rho, c
        real(wp), intent(in) :: dpres_ds

        L(1) = f_base_L1(lambda, rho, c, dpres_ds, dvel_ds)
        call s_fill_density_L(L, 1._wp, lambda(2), c, mf, dalpha_rho_ds, dpres_ds)
        call s_fill_velocity_L(L, 1._wp, lambda(2), dvel_ds)
        call s_fill_advection_L(L, 1._wp, lambda(2), dadv_ds)
        L(eqn_idx%adv%end) = L(1) + 2._wp*rho*c*lambda(2)*dvel_ds(dir_idx(1))

    end subroutine s_compute_force_free_subsonic_outflow_L

    !> Constant pressure subsonic outflow CBC (Thompson 1990, pg. 455)
    subroutine s_compute_constant_pressure_subsonic_outflow_L(lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds)

        $:GPU_ROUTINE(function_name='s_compute_constant_pressure_subsonic_outflow_L', parallelism='[seq]', cray_inline=True)

        real(wp), dimension(3), intent(in) :: lambda
        #:if USING_AMD
            real(wp), dimension(20), intent(inout) :: L
        #:else
            real(wp), dimension(sys_size), intent(inout) :: L
        #:endif
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(in) :: mf, dalpha_rho_ds
            real(wp), dimension(3), intent(in) :: dvel_ds
            real(wp), dimension(3), intent(in) :: dadv_ds
        #:else
            real(wp), dimension(num_fluids), intent(in) :: mf, dalpha_rho_ds
            real(wp), dimension(num_dims), intent(in)   :: dvel_ds
            real(wp), dimension(num_fluids), intent(in) :: dadv_ds
        #:endif
        real(wp), intent(in) :: rho, c
        real(wp), intent(in) :: dpres_ds

        L(1) = f_base_L1(lambda, rho, c, dpres_ds, dvel_ds)
        call s_fill_density_L(L, 1._wp, lambda(2), c, mf, dalpha_rho_ds, dpres_ds)
        call s_fill_velocity_L(L, 1._wp, lambda(2), dvel_ds)
        call s_fill_advection_L(L, 1._wp, lambda(2), dadv_ds)
        L(eqn_idx%adv%end) = -L(1)

    end subroutine s_compute_constant_pressure_subsonic_outflow_L

    !> Supersonic inflow CBC (Thompson 1990, pg. 453)
    subroutine s_compute_supersonic_inflow_L(L)

        $:GPU_ROUTINE(function_name='s_compute_supersonic_inflow_L', parallelism='[seq]', cray_inline=True)
        #:if USING_AMD
            real(wp), dimension(20), intent(inout) :: L
        #:else
            real(wp), dimension(sys_size), intent(inout) :: L
        #:endif
        L(1:eqn_idx%adv%end) = 0._wp
        if (chemistry) L(eqn_idx%species%beg:eqn_idx%species%end) = 0._wp

    end subroutine s_compute_supersonic_inflow_L

    !> Supersonic outflow CBC (Thompson 1990, pg. 453)
    subroutine s_compute_supersonic_outflow_L(lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds, dYs_ds)

        $:GPU_ROUTINE(function_name='s_compute_supersonic_outflow_L', parallelism='[seq]', cray_inline=True)

        real(wp), dimension(3), intent(in) :: lambda
        #:if USING_AMD
            real(wp), dimension(20), intent(inout) :: L
        #:else
            real(wp), dimension(sys_size), intent(inout) :: L
        #:endif
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(in)  :: mf, dalpha_rho_ds
            real(wp), dimension(3), intent(in)  :: dvel_ds
            real(wp), dimension(3), intent(in)  :: dadv_ds
            real(wp), dimension(10), intent(in) :: dYs_ds
        #:else
            real(wp), dimension(num_fluids), intent(in)  :: mf, dalpha_rho_ds
            real(wp), dimension(num_dims), intent(in)    :: dvel_ds
            real(wp), dimension(num_fluids), intent(in)  :: dadv_ds
            real(wp), dimension(num_species), intent(in) :: dYs_ds
        #:endif
        real(wp), intent(in) :: rho, c
        real(wp), intent(in) :: dpres_ds

        L(1) = f_base_L1(lambda, rho, c, dpres_ds, dvel_ds)
        call s_fill_density_L(L, 1._wp, lambda(2), c, mf, dalpha_rho_ds, dpres_ds)
        call s_fill_velocity_L(L, 1._wp, lambda(2), dvel_ds)
        call s_fill_advection_L(L, 1._wp, lambda(2), dadv_ds)
        call s_fill_chemistry_L(L, 1._wp, lambda(2), dYs_ds)
        L(eqn_idx%adv%end) = lambda(3)*(dpres_ds + rho*c*dvel_ds(dir_idx(1)))

    end subroutine s_compute_supersonic_outflow_L

end module m_compute_cbc
