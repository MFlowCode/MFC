#:include 'macros.fpp'

module m_additional_forcing
    use m_derived_types

    use m_global_parameters

    use m_ibm

    use m_mpi_proxy

    use m_volume_filtering

    implicit none

    private; public :: s_initialize_additional_forcing_module, &
 s_finalize_additional_forcing_module, s_compute_periodic_forcing

    type(scalar_field), allocatable, dimension(:) :: q_periodic_force
    real(wp) :: volfrac_phi
    integer :: N_x_total_glb
    real(wp) :: avg_coeff
    real(wp) :: spatial_rho, spatial_u, spatial_eps
    real(wp), allocatable, dimension(:) :: rho_window, u_window, eps_window
    real(wp) :: sum_rho, sum_u, sum_eps
    real(wp) :: phase_rho, phase_u, phase_eps
    real(wp) :: tau_dt
    integer :: forcing_window, window_loc, window_fill

    $:GPU_DECLARE(create='[q_periodic_force, volfrac_phi, N_x_total_glb, avg_coeff, tau_dt]')
    $:GPU_DECLARE(create='[spatial_rho, spatial_u, spatial_eps, phase_rho, phase_u, phase_eps]')

contains

    subroutine s_initialize_additional_forcing_module
        integer :: i

        @:ALLOCATE(q_periodic_force(1:3))
        do i = 1, 3
            @:ALLOCATE(q_periodic_force(i)%sf(0:m, 0:n, 0:p))
            @:ACC_SETUP_SFs(q_periodic_force(i))
        end do

        volfrac_phi = num_ibs*4._wp/3._wp*pi*patch_ib(1)%radius**3/((x_domain%end - x_domain%beg)*(y_domain%end - y_domain%beg)*(z_domain%end - z_domain%beg))
        $:GPU_UPDATE(device='[volfrac_phi]')

        N_x_total_glb = (m_glb + 1)*(n_glb + 1)*(p_glb + 1)
        $:GPU_UPDATE(device='[N_x_total_glb]')

        avg_coeff = 1._wp/(real(N_x_total_glb, wp)*(1._wp - volfrac_phi))
        $:GPU_UPDATE(device='[avg_coeff]')

        forcing_window = 1 ! forcing time window size
        window_loc = 0
        window_fill = 0

        @:ALLOCATE(rho_window(forcing_window))
        @:ALLOCATE(u_window(forcing_window))
        @:ALLOCATE(eps_window(forcing_window))

        rho_window = 0.0_wp
        u_window = 0.0_wp
        eps_window = 0.0_wp

        sum_rho = 0.0_wp
        sum_u = 0.0_wp
        sum_eps = 0.0_wp

        phase_rho = 0._wp
        phase_u = 0._wp
        phase_eps = 0._wp

        tau_dt = 1._wp/(0.5_wp*dt)
        $:GPU_UPDATE(device='[tau_dt]')

    end subroutine s_initialize_additional_forcing_module

    !< compute the space and time average of quantities, compute the periodic forcing terms described in Khalloufi and Capecelatro
    subroutine s_compute_periodic_forcing(rhs_vf, q_cons_vf, q_prim_vf, t_step)
        type(scalar_field), dimension(sys_size), intent(in) :: rhs_vf
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        integer, intent(in) :: t_step
        real(wp) :: spatial_rho_glb, spatial_u_glb, spatial_eps_glb
        integer :: i, j, k

        ! zero spatial averages
        spatial_rho = 0._wp
        spatial_u = 0._wp
        spatial_eps = 0._wp
        $:GPU_UPDATE(device='[spatial_rho, spatial_u, spatial_eps]')

        ! compute spatial averages
        $:GPU_PARALLEL_LOOP(collapse=3, reduction='[[spatial_rho, spatial_u, spatial_eps]]', reductionOp='[+]')
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    spatial_rho = spatial_rho + (q_cons_vf(1)%sf(i, j, k)*fluid_indicator_function%sf(i, j, k)) ! rho
                    spatial_u = spatial_u + (q_cons_vf(2)%sf(i, j, k)*fluid_indicator_function%sf(i, j, k)) ! rho*u
                    spatial_eps = spatial_eps + ((q_cons_vf(5)%sf(i, j, k) - 0.5_wp*( &
                                                  q_cons_vf(2)%sf(i, j, k)**2 + &
                                                  q_cons_vf(3)%sf(i, j, k)**2 + &
                                                  q_cons_vf(4)%sf(i, j, k)**2)/q_cons_vf(1)%sf(i, j, k))*fluid_indicator_function%sf(i, j, k)) ! rho*e
                end do
            end do
        end do

        $:GPU_UPDATE(host='[spatial_rho, spatial_u, spatial_eps]')

        ! reduction sum across entire domain
        call s_mpi_allreduce_sum(spatial_rho, spatial_rho_glb)
        call s_mpi_allreduce_sum(spatial_u, spatial_u_glb)
        call s_mpi_allreduce_sum(spatial_eps, spatial_eps_glb)

        spatial_rho_glb = spatial_rho_glb*avg_coeff
        spatial_u_glb = spatial_u_glb*avg_coeff
        spatial_eps_glb = spatial_eps_glb*avg_coeff

        ! update time average window location
        window_loc = 1 + mod(t_step - 1, forcing_window)

        ! update time average sum
        sum_rho = sum_rho - rho_window(window_loc) + spatial_rho_glb
        sum_u = sum_u - u_window(window_loc) + spatial_u_glb
        sum_eps = sum_eps - eps_window(window_loc) + spatial_eps_glb

        ! update window arrays
        rho_window(window_loc) = spatial_rho_glb
        u_window(window_loc) = spatial_u_glb
        eps_window(window_loc) = spatial_eps_glb

        ! update number of time samples
        if (window_fill < forcing_window) window_fill = window_fill + 1

        ! compute phase averages
        phase_rho = sum_rho/real(window_fill, wp)
        phase_u = sum_u/real(window_fill, wp)
        phase_eps = sum_eps/real(window_fill, wp)
        $:GPU_UPDATE(device='[phase_rho, phase_u, phase_eps]')

        ! compute periodic forcing terms for mass, momentum, energy
        $:GPU_PARALLEL_LOOP(collapse=3)
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    ! f_rho
                    q_periodic_force(1)%sf(i, j, k) = (rho_inf_ref - phase_rho)*tau_dt

                    ! f_u
                    q_periodic_force(2)%sf(i, j, k) = (rho_inf_ref*u_inf_ref - phase_u)*tau_dt

                    ! f_E
                    q_periodic_force(3)%sf(i, j, k) = (P_inf_ref*gammas(1) - phase_eps)*tau_dt &
                                                      + q_cons_vf(2)%sf(i, j, k)*q_periodic_force(2)%sf(i, j, k)/q_cons_vf(1)%sf(i, j, k)
                end do
            end do
        end do

        ! add the forcing terms to the RHS
        $:GPU_PARALLEL_LOOP(collapse=3)
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    rhs_vf(1)%sf(i, j, k) = rhs_vf(1)%sf(i, j, k) + q_periodic_force(1)%sf(i, j, k)*fluid_indicator_function%sf(i, j, k) ! continuity
                    rhs_vf(2)%sf(i, j, k) = rhs_vf(2)%sf(i, j, k) + q_periodic_force(2)%sf(i, j, k)*fluid_indicator_function%sf(i, j, k) ! x momentum
                    rhs_vf(5)%sf(i, j, k) = rhs_vf(5)%sf(i, j, k) + q_periodic_force(3)%sf(i, j, k)*fluid_indicator_function%sf(i, j, k) ! energy
                end do
            end do
        end do

    end subroutine s_compute_periodic_forcing

    subroutine s_finalize_additional_forcing_module
        integer :: i
        do i = 1, 3
            @:DEALLOCATE(q_periodic_force(i)%sf)
        end do
        @:DEALLOCATE(q_periodic_force)
    end subroutine s_finalize_additional_forcing_module

end module m_additional_forcing
