#:include 'macros.fpp'

module m_additional_forcing
    use m_derived_types 

    use m_global_parameters

    use m_ibm

    use m_mpi_proxy 

    use m_volume_filtering

    implicit none

    private; public :: s_initialize_additional_forcing_module, & 
 s_add_periodic_forcing, s_finalize_additional_forcing_module, & 
 s_compute_periodic_forcing

    type(scalar_field), allocatable, dimension(:) :: q_periodic_force
    real(wp) :: volfrac_phi
    integer :: N_x_total_glb
    real(wp) :: spatial_rho, spatial_u
    real(wp) :: phase_rho, phase_u

    !$acc declare create(q_periodic_force, volfrac_phi, N_x_total_glb)
    !$acc declare create(spatial_rho, spatial_u, phase_rho, phase_u)

contains

    subroutine s_initialize_additional_forcing_module
        integer :: i

        @:ALLOCATE(q_periodic_force(1:3))
        do i = 1, 3
            @:ALLOCATE(q_periodic_force(i)%sf(0:m, 0:n, 0:p))
            @:ACC_SETUP_SFs(q_periodic_force(i))
        end do

        volfrac_phi = num_ibs * 4._wp/3._wp * pi * patch_ib(1)%radius**3 / ((x_domain%end - x_domain%beg)*(y_domain%end - y_domain%beg)*(z_domain%end - z_domain%beg))
        !$acc update device(volfrac_phi)

        N_x_total_glb = (m_glb + 1) * (n_glb + 1) * (p_glb + 1)
        !$acc update device(N_x_total_glb)
    end subroutine s_initialize_additional_forcing_module

    !< adds periodic forcing terms to RHS, as detailed in Khalloufi and Capecelatro
    subroutine s_add_periodic_forcing(rhs_vf)
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        integer :: i, j, k

        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    rhs_vf(1)%sf(i, j, k) = rhs_vf(1)%sf(i, j, k) + q_periodic_force(1)%sf(i, j, k) * fluid_indicator_function%sf(i, j, k) ! continuity
                    rhs_vf(2)%sf(i, j, k) = rhs_vf(2)%sf(i, j, k) + q_periodic_force(2)%sf(i, j, k) * fluid_indicator_function%sf(i, j, k) ! x momentum
                    rhs_vf(5)%sf(i, j, k) = rhs_vf(5)%sf(i, j, k) + q_periodic_force(3)%sf(i, j, k) * fluid_indicator_function%sf(i, j, k) ! energy
                end do
            end do
        end do
    end subroutine s_add_periodic_forcing

    !< compute the space and time average of quantities, compute the periodic forcing terms described in Khalloufi and Capecelatro
    subroutine s_compute_periodic_forcing(q_cons_vf, t_step)
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        integer, intent(in) :: t_step
        real(wp) :: spatial_rho_glb, spatial_u_glb
        integer :: i, j, k

        ! zero spatial averages
        spatial_rho = 0._wp
        spatial_u = 0._wp
        !$acc update device(spatial_rho, spatial_u)

        ! compute spatial averages
        !$acc parallel loop collapse(3) gang vector default(present) reduction(+:spatial_rho, spatial_u)
        do i = 0, m 
            do j = 0, n 
                do k = 0, p 
                    spatial_rho = spatial_rho + q_cons_vf(1)%sf(i, j, k) * fluid_indicator_function%sf(i, j, k) ! rho
                    spatial_u = spatial_u + q_cons_vf(2)%sf(i, j, k) * fluid_indicator_function%sf(i, j, k) ! u
                end do
            end do
        end do

        !$acc update host(spatial_rho, spatial_u)

        ! reduction sum across entire domain
        call s_mpi_allreduce_sum(spatial_rho, spatial_rho_glb)
        call s_mpi_allreduce_sum(spatial_u, spatial_u_glb)

        ! compute phase averages
        phase_rho = phase_rho + (spatial_rho_glb / real(N_x_total_glb, wp) - phase_rho) / real(t_step, wp)
        phase_u = phase_u + (spatial_u_glb / real(N_x_total_glb, wp) - phase_u) / real(t_step, wp)
        !$acc update device(phase_rho, phase_u)

        ! compute periodic forcing terms for mass, momentum, energy
        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    ! f_rho
                    q_periodic_force(1)%sf(i, j, k) = (rho_inf_ref - phase_rho/(1._wp - volfrac_phi)) / dt

                    ! f_u
                    q_periodic_force(2)%sf(i, j, k) = (rho_inf_ref*u_inf_ref - phase_u/(1._wp - volfrac_phi)) / dt

                    ! u*f_u
                    q_periodic_force(3)%sf(i, j, k) = q_cons_vf(2)%sf(i, j, k)/q_cons_vf(1)%sf(i, j, k) * q_periodic_force(2)%sf(i, j, k)
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
