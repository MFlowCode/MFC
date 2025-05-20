#:include 'macros.fpp'

module m_compute_particle_forces
    use m_derived_types 

    use m_global_parameters

    use m_ibm

    use m_mpi_proxy 

    implicit none

    private; public :: s_initialize_particle_forces_module, & 
 s_compute_drag_coefficient, s_finalize_particle_forces_module

    real(wp), allocatable, dimension(:) :: FD_calc

    !$acc declare create(FD_calc)

contains
    
    subroutine s_initialize_particle_forces_module
        if (compute_CD) then
            @:ALLOCATE(FD_calc(0:num_ibs))
        end if

    end subroutine s_initialize_particle_forces_module

    subroutine s_compute_drag_coefficient(div_pres_visc_stress)
        type(scalar_field), dimension(momxb:momxe), intent(in) :: div_pres_visc_stress
        real(wp), dimension(0:num_ibs) :: FD_global
        real(wp) :: drag_coeff 
        integer :: i, j, k

        !$acc parallel loop gang vector default(present)
        do i = 0, num_ibs 
            FD_calc(i) = 0._wp
        end do

        !$acc parallel loop collapse(3) gang vector default(present)
        do i = 0, m 
            do j = 0, n 
                do k = 0, p  
                    !$acc atomic
                    FD_calc(ib_markers%sf(i, j, k)) = FD_calc(ib_markers%sf(i, j, k)) & 
                                                    + div_pres_visc_stress(momxb)%sf(i, j, k) * dx(i) * dy(j) * dz(k)
                end do 
            end do 
        end do

        !$acc update host(FD_calc(:))

        do i = 0, num_ibs 
            call s_mpi_allreduce_sum(FD_calc(i), FD_global(i))
        end do

        drag_coeff = FD_global(1) / (0.5_wp * rho_inf_ref * (u_inf_ref**2) * pi * (patch_ib(1)%radius**2))
        if (proc_rank == 0) then 
            print *, 'C_D: ', drag_coeff
        end if

    end subroutine s_compute_drag_coefficient

    subroutine s_finalize_particle_forces_module
        if (compute_CD) then 
            @:DEALLOCATE(FD_calc)
        end if

    end subroutine s_finalize_particle_forces_module
    
end module m_compute_particle_forces
