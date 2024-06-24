!>
!! @file m_compute_cbc.f90
!! @brief Contains module m_compute_cbc

module m_compute_cbc

    ! Dependencies =============================================================

    use m_global_parameters    !< Definitions of the global parameters

    ! ==========================================================================

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

    !>  The L variables for the slip wall CBC, see pg. 451 of
        !!      Thompson (1990). At the slip wall (frictionless wall),
        !!      the normal component of velocity is zero at all times,
        !!      while the transverse velocities may be nonzero.
    subroutine s_compute_slip_wall_L(lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds)
#ifdef CRAY_ACC_WAR
        !DIR$ INLINEALWAYS s_compute_slip_wall_L
#else
        !$acc routine seq
#endif
        real(kind(0d0)), dimension(3), intent(in) :: lambda
        real(kind(0d0)), dimension(sys_size), intent(inout) :: L
        real(kind(0d0)), intent(in) :: rho, c
        real(kind(0d0)), dimension(num_fluids), intent(in) :: mf, dalpha_rho_ds
        real(kind(0d0)), intent(in) :: dpres_ds
        real(kind(0d0)), dimension(num_dims), intent(in) :: dvel_ds
        real(kind(0d0)), dimension(num_fluids), intent(in) :: dadv_ds

        integer :: i

        L(1) = lambda(1)*(dpres_ds - rho*c*dvel_ds(dir_idx(1)))

        do i = 2, advxe
            L(i) = 0d0
        end do

        L(advxe) = L(1)

    end subroutine s_compute_slip_wall_L

    !>  The L variables for the nonreflecting subsonic buffer CBC
        !!      see pg. 13 of Thompson (1987). The nonreflecting subsonic
        !!      buffer reduces the amplitude of any reflections caused by
        !!      outgoing waves.
    subroutine s_compute_nonreflecting_subsonic_buffer_L(lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds)
#ifdef CRAY_ACC_WAR
        !DIR$ INLINEALWAYS s_compute_nonreflecting_subsonic_buffer_L
#else
        !$acc routine seq
#endif
        real(kind(0d0)), dimension(3), intent(in) :: lambda
        real(kind(0d0)), dimension(sys_size), intent(inout) :: L
        real(kind(0d0)), intent(in) :: rho, c
        real(kind(0d0)), dimension(num_fluids), intent(in) :: mf, dalpha_rho_ds
        real(kind(0d0)), intent(in) :: dpres_ds
        real(kind(0d0)), dimension(num_dims), intent(in) :: dvel_ds
        real(kind(0d0)), dimension(num_fluids), intent(in) :: dadv_ds

        integer :: i !< Generic loop iterator

        L(1) = (5d-1 - 5d-1*sign(1d0, lambda(1)))*lambda(1) &
               *(dpres_ds - rho*c*dvel_ds(dir_idx(1)))

        do i = 2, momxb
            L(i) = (5d-1 - 5d-1*sign(1d0, lambda(2)))*lambda(2) &
                   *(c*c*dalpha_rho_ds(i - 1) - mf(i - 1)*dpres_ds)
        end do

        do i = momxb + 1, momxe
            L(i) = (5d-1 - 5d-1*sign(1d0, lambda(2)))*lambda(2) &
                   *(dvel_ds(dir_idx(i - contxe)))
        end do

        do i = E_idx, advxe - 1
            L(i) = (5d-1 - 5d-1*sign(1d0, lambda(2)))*lambda(2) &
                   *(dadv_ds(i - momxe))
        end do

        L(advxe) = (5d-1 - 5d-1*sign(1d0, lambda(3)))*lambda(3) &
                   *(dpres_ds + rho*c*dvel_ds(dir_idx(1)))

    end subroutine s_compute_nonreflecting_subsonic_buffer_L
    !>  The L variables for the nonreflecting subsonic inflow CBC
        !!      see pg. 455, Thompson (1990). This nonreflecting subsonic
        !!      CBC assumes an incoming flow and reduces the amplitude of
        !!      any reflections caused by outgoing waves.
    subroutine s_compute_nonreflecting_subsonic_inflow_L(lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds)
#ifdef CRAY_ACC_WAR
        !DIR$ INLINEALWAYS ss_compute_nonreflecting_subsonic_inflow_L
#else
        !$acc routine seq
#endif
        real(kind(0d0)), dimension(3), intent(in) :: lambda
        real(kind(0d0)), dimension(sys_size), intent(inout) :: L
        real(kind(0d0)), intent(in) :: rho, c
        real(kind(0d0)), dimension(num_fluids), intent(in) :: mf, dalpha_rho_ds
        real(kind(0d0)), intent(in) :: dpres_ds
        real(kind(0d0)), dimension(num_dims), intent(in) :: dvel_ds
        real(kind(0d0)), dimension(num_fluids), intent(in) :: dadv_ds

        integer :: i

        L(1) = lambda(1)*(dpres_ds - rho*c*dvel_ds(dir_idx(1)))

        do i = 2, advxe
            L(i) = 0d0
        end do

    end subroutine s_compute_nonreflecting_subsonic_inflow_L

    !>  The L variables for the nonreflecting subsonic outflow
        !!      CBC see pg. 454 of Thompson (1990). This nonreflecting
        !!      subsonic CBC presumes an outgoing flow and reduces the
        !!      amplitude of any reflections caused by outgoing waves.
    subroutine s_compute_nonreflecting_subsonic_outflow_L(lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds)
#ifdef CRAY_ACC_WAR
        !DIR$ INLINEALWAYS s_compute_nonreflecting_subsonic_outflow_L
#else
        !$acc routine seq
#endif
        real(kind(0d0)), dimension(3), intent(in) :: lambda
        real(kind(0d0)), dimension(sys_size), intent(inout) :: L
        real(kind(0d0)), intent(in) :: rho, c
        real(kind(0d0)), dimension(num_fluids), intent(in) :: mf, dalpha_rho_ds
        real(kind(0d0)), intent(in) :: dpres_ds
        real(kind(0d0)), dimension(num_dims), intent(in) :: dvel_ds
        real(kind(0d0)), dimension(num_fluids), intent(in) :: dadv_ds

        integer :: i !> Generic loop iterator

        L(1) = lambda(1)*(dpres_ds - rho*c*dvel_ds(dir_idx(1)))

        do i = 2, momxb
            L(i) = lambda(2)*(c*c*dalpha_rho_ds(i - 1) - mf(i - 1)*dpres_ds)
        end do

        do i = momxb + 1, momxe
            L(i) = lambda(2)*(dvel_ds(dir_idx(i - contxe)))
        end do

        do i = E_idx, advxe - 1
            L(i) = lambda(2)*(dadv_ds(i - momxe))
        end do

        ! bubble index
        L(advxe) = 0d0

    end subroutine s_compute_nonreflecting_subsonic_outflow_L

    !>  The L variables for the force-free subsonic outflow CBC,
        !!      see pg. 454 of Thompson (1990). The force-free subsonic
        !!      outflow sets to zero the sum of all of the forces which
        !!      are acting on a fluid element for the normal coordinate
        !!      direction to the boundary. As a result, a fluid element
        !!      at the boundary is simply advected outward at the fluid
        !!      velocity.
    subroutine s_compute_force_free_subsonic_outflow_L(lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds)
#ifdef CRAY_ACC_WAR
        !DIR$ INLINEALWAYS s_compute_force_free_subsonic_outflow_L
#else
        !$acc routine seq
#endif
        real(kind(0d0)), dimension(3), intent(in) :: lambda
        real(kind(0d0)), dimension(sys_size), intent(inout) :: L
        real(kind(0d0)), intent(in) :: rho, c
        real(kind(0d0)), dimension(num_fluids), intent(in) :: mf, dalpha_rho_ds
        real(kind(0d0)), intent(in) :: dpres_ds
        real(kind(0d0)), dimension(num_dims), intent(in) :: dvel_ds
        real(kind(0d0)), dimension(num_fluids), intent(in) :: dadv_ds

        integer :: i !> Generic loop iterator

        L(1) = lambda(1)*(dpres_ds - rho*c*dvel_ds(dir_idx(1)))

        do i = 2, momxb
            L(i) = lambda(2)*(c*c*dalpha_rho_ds(i - 1) - mf(i - 1)*dpres_ds)
        end do

        do i = momxb + 1, momxe
            L(i) = lambda(2)*(dvel_ds(dir_idx(i - contxe)))
        end do

        do i = E_idx, advxe - 1
            L(i) = lambda(2)*(dadv_ds(i - momxe))
        end do

        L(advxe) = L(1) + 2d0*rho*c*lambda(2)*dvel_ds(dir_idx(1))

    end subroutine s_compute_force_free_subsonic_outflow_L

    !>  L variables for the constant pressure subsonic outflow
        !!      CBC see pg. 455 Thompson (1990). The constant pressure
        !!      subsonic outflow maintains a fixed pressure at the CBC
        !!      boundary in absence of any transverse effects.
    subroutine s_compute_constant_pressure_subsonic_outflow_L(lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds)
#ifdef CRAY_ACC_WAR
        !DIR$ INLINEALWAYS s_compute_constant_pressure_subsonic_outflow_L
#else
        !$acc routine seq
#endif
        real(kind(0d0)), dimension(3), intent(in) :: lambda
        real(kind(0d0)), dimension(sys_size), intent(inout) :: L
        real(kind(0d0)), intent(in) :: rho, c
        real(kind(0d0)), dimension(num_fluids), intent(in) :: mf, dalpha_rho_ds
        real(kind(0d0)), intent(in) :: dpres_ds
        real(kind(0d0)), dimension(num_dims), intent(in) :: dvel_ds
        real(kind(0d0)), dimension(num_fluids), intent(in) :: dadv_ds

        integer :: i !> Generic loop iterator

        L(1) = lambda(1)*(dpres_ds - rho*c*dvel_ds(dir_idx(1)))

        do i = 2, momxb
            L(i) = lambda(2)*(c*c*dalpha_rho_ds(i - 1) - mf(i - 1)*dpres_ds)
        end do

        do i = momxb + 1, momxe
            L(i) = lambda(2)*(dvel_ds(dir_idx(i - contxe)))
        end do

        do i = E_idx, advxe - 1
            L(i) = lambda(2)*(dadv_ds(i - momxe))
        end do

        L(advxe) = -L(1)

    end subroutine s_compute_constant_pressure_subsonic_outflow_L

    !>  L variables for the supersonic inflow CBC, see pg. 453
        !!      Thompson (1990). The supersonic inflow CBC is a steady
        !!      state, or nearly a steady state, CBC in which only the
        !!      transverse terms may generate a time dependence at the
        !!      inflow boundary.
    subroutine s_compute_supersonic_inflow_L(lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds)
#ifdef CRAY_ACC_WAR
        !DIR$ INLINEALWAYS s_compute_supersonic_inflow_L
#else
        !$acc routine seq
#endif
        real(kind(0d0)), dimension(3), intent(in) :: lambda
        real(kind(0d0)), dimension(sys_size), intent(inout) :: L
        real(kind(0d0)), intent(in) :: rho, c
        real(kind(0d0)), dimension(num_fluids), intent(in) :: mf, dalpha_rho_ds
        real(kind(0d0)), intent(in) :: dpres_ds
        real(kind(0d0)), dimension(num_dims), intent(in) :: dvel_ds
        real(kind(0d0)), dimension(num_fluids), intent(in) :: dadv_ds
        integer :: i

        do i = 1, advxe
            L(i) = 0d0
        end do

    end subroutine s_compute_supersonic_inflow_L

    !>  L variables for the supersonic outflow CBC, see pg. 453
        !!      of Thompson (1990). For the supersonic outflow CBC, the
        !!      flow evolution at the boundary is determined completely
        !!      by the interior data.
    subroutine s_compute_supersonic_outflow_L(lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds)
#ifdef CRAY_ACC_WAR
        !DIR$ INLINEALWAYS s_compute_supersonic_outflow_L
#else
        !$acc routine seq
#endif
        real(kind(0d0)), dimension(3), intent(in) :: lambda
        real(kind(0d0)), dimension(sys_size), intent(inout) :: L
        real(kind(0d0)), intent(in) :: rho, c
        real(kind(0d0)), dimension(num_fluids), intent(in) :: mf, dalpha_rho_ds
        real(kind(0d0)), intent(in) :: dpres_ds
        real(kind(0d0)), dimension(num_dims), intent(in) :: dvel_ds
        real(kind(0d0)), dimension(num_fluids), intent(in) :: dadv_ds

        integer :: i !< Generic loop iterator

        L(1) = lambda(1)*(dpres_ds - rho*c*dvel_ds(dir_idx(1)))

        do i = 2, momxb
            L(i) = lambda(2)*(c*c*dalpha_rho_ds(i - 1) - mf(i - 1)*dpres_ds)
        end do

        do i = momxb + 1, momxe
            L(i) = lambda(2)*(dvel_ds(dir_idx(i - contxe)))
        end do

        do i = E_idx, advxe - 1
            L(i) = lambda(2)*(dadv_ds(i - momxe))
        end do

        L(advxe) = lambda(3)*(dpres_ds + rho*c*dvel_ds(dir_idx(1)))

    end subroutine s_compute_supersonic_outflow_L

end module m_compute_cbc
