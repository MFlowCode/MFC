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
 s_compute_phase_average, s_compute_periodic_forcing;

    real(wp), allocatable, dimension(:) :: q_bar ! 1:3 rho*u, 4 rho, 5 T
    type(scalar_field), allocatable, dimension(:) :: q_periodic_force
    real(wp), allocatable, dimension(:) :: q_spatial_avg
    real(wp), allocatable, dimension(:), public :: q_spatial_avg_glb ! 1:3 rho*u, 4 rho, 5 T
    real(wp) :: volfrac_phi
    integer :: N_x_total_glb

    !$acc declare create(q_bar, q_periodic_force, q_spatial_avg, q_spatial_avg_glb, volfrac_phi, N_x_total_glb)

contains

    subroutine s_initialize_additional_forcing_module
        integer :: i
        if (periodic_forcing) then 
            @:ALLOCATE(q_bar(1:5))
            @:ALLOCATE(q_periodic_force(1:8))
            do i = 1, 8 
                @:ALLOCATE(q_periodic_force(i)%sf(0:m, 0:n, 0:p))
                @:ACC_SETUP_SFs(q_periodic_force(i))
            end do
            @:ALLOCATE(q_spatial_avg(1:5))
            @:ALLOCATE(q_spatial_avg_glb(1:5))
        end if

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
                    rhs_vf(1)%sf(i, j, k) = rhs_vf(1)%sf(i, j, k) + q_periodic_force(7)%sf(i, j, k) * fluid_indicator_function%sf(i, j, k) ! continuity
                    rhs_vf(2)%sf(i, j, k) = rhs_vf(2)%sf(i, j, k) + q_periodic_force(1)%sf(i, j, k) * fluid_indicator_function%sf(i, j, k) * fluid_indicator_function%sf(i, j, k) ! x momentum
                    rhs_vf(5)%sf(i, j, k) = rhs_vf(5)%sf(i, j, k) + (q_periodic_force(4)%sf(i, j, k) + q_periodic_force(8)%sf(i, j, k)) * fluid_indicator_function%sf(i, j, k) ! energy
                end do
            end do
        end do
    end subroutine s_add_periodic_forcing

    subroutine s_compute_phase_average(q_cons_vf, t_step)
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        integer, intent(in) :: t_step
        integer :: i, j, k

        !$acc loop seq
        do i = 1, 5
            q_spatial_avg(i) = 0._wp
        end do

        ! spatial average
        !$acc parallel loop collapse(3) gang vector default(present) reduction(+:q_spatial_avg(:))
        do i = 0, m 
            do j = 0, n 
                do k = 0, p 
                    q_spatial_avg(4) = q_spatial_avg(4) + q_cons_vf(1)%sf(i, j, k) * fluid_indicator_function%sf(i, j, k)
                    q_spatial_avg(5) = q_spatial_avg(5) + (0.4_wp/287._wp * (q_cons_vf(5)%sf(i, j, k)/q_cons_vf(1)%sf(i, j, k) & 
                                        - 0.5_wp * ((q_cons_vf(2)%sf(i, j, k)/q_cons_vf(1)%sf(i, j, k))**2 & 
                                        + (q_cons_vf(3)%sf(i, j, k)/q_cons_vf(1)%sf(i, j, k))**2 & 
                                        + (q_cons_vf(4)%sf(i, j, k)/q_cons_vf(1)%sf(i, j, k))**2))) * fluid_indicator_function%sf(i, j, k)
                                        
                    q_spatial_avg(1) = q_spatial_avg(1) + (q_cons_vf(2)%sf(i, j, k)) * fluid_indicator_function%sf(i, j, k)
                    q_spatial_avg(2) = q_spatial_avg(2) + (q_cons_vf(3)%sf(i, j, k)) * fluid_indicator_function%sf(i, j, k)
                    q_spatial_avg(3) = q_spatial_avg(3) + (q_cons_vf(4)%sf(i, j, k)) * fluid_indicator_function%sf(i, j, k)
                end do
            end do
        end do

        !$acc update host(q_spatial_avg(:))

        do i = 1, 5 
            call s_mpi_allreduce_sum(q_spatial_avg(i), q_spatial_avg_glb(i))
        end do

        !$acc update device(q_spatial_avg_glb(:))

        !$acc loop seq
        do i = 1, 5 
            q_spatial_avg_glb(i) = q_spatial_avg_glb(i) / real(N_x_total_glb, wp)
        end do

        ! time average
        !$acc loop seq
        do i = 1, 5 
            q_bar(i) = ( (q_spatial_avg_glb(i) + (t_step - 1._wp)*q_bar(i)) / t_step ) 
        end do
    end subroutine s_compute_phase_average

    !< computes the periodic forcing terms described in Khalloufi and Capecelatro
    subroutine s_compute_periodic_forcing(q_cons_vf)
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf

        integer :: i, j, k

        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    ! f_u
                    q_periodic_force(1)%sf(i, j, k) = (rho_inf_ref*u_inf_ref - q_bar(1)/(1._wp - volfrac_phi)) / dt
                    q_periodic_force(2)%sf(i, j, k) = (rho_inf_ref*u_inf_ref - q_bar(2)/(1._wp - volfrac_phi)) / dt
                    q_periodic_force(3)%sf(i, j, k) = (rho_inf_ref*u_inf_ref - q_bar(3)/(1._wp - volfrac_phi)) / dt

                    ! u*f_u
                    q_periodic_force(4)%sf(i, j, k) = q_cons_vf(2)%sf(i, j, k)/q_cons_vf(1)%sf(i, j, k) * q_periodic_force(1)%sf(i, j, k)
                    q_periodic_force(5)%sf(i, j, k) = q_cons_vf(3)%sf(i, j, k)/q_cons_vf(1)%sf(i, j, k) * q_periodic_force(2)%sf(i, j, k)
                    q_periodic_force(6)%sf(i, j, k) = q_cons_vf(4)%sf(i, j, k)/q_cons_vf(1)%sf(i, j, k) * q_periodic_force(3)%sf(i, j, k)

                    ! f_rho
                    q_periodic_force(7)%sf(i, j, k) = (rho_inf_ref - q_bar(4)/(1._wp - volfrac_phi)) / dt

                    ! f_T
                    q_periodic_force(8)%sf(i, j, k) = (q_cons_vf(1)%sf(i, j, k) / 1.4_wp) * (T_inf_ref - q_bar(5)/(1._wp - volfrac_phi)) / dt
                end do 
            end do
        end do
    end subroutine s_compute_periodic_forcing

    subroutine s_finalize_additional_forcing_module
        integer :: i
        if (periodic_forcing) then
            @:DEALLOCATE(q_bar)
            do i = 1, 8
                @:DEALLOCATE(q_periodic_force(i)%sf)
            end do
            @:DEALLOCATE(q_periodic_force)
            @:DEALLOCATE(q_spatial_avg)
            @:DEALLOCATE(q_spatial_avg_glb)
        end if
    end subroutine s_finalize_additional_forcing_module

end module m_additional_forcing