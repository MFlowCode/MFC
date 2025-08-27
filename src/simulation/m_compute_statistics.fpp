#:include 'macros.fpp'

module m_compute_statistics
    use m_derived_types 

    use m_global_parameters

    use m_mpi_proxy 

    use m_additional_forcing

    use m_nvtx

    implicit none

    private; public :: s_initialize_statistics_module, s_finalize_statistics_module, &
    s_compute_statistics_momentum_unclosed_terms, s_update_statistics, &
    s_compute_234_order_statistics
 
    ! terms for computing 1st, 2nd, 3rd, and 4th order statistical moments
    type(scalar_field), allocatable, dimension(:) :: Msn_reynolds_stress
    type(scalar_field), allocatable, dimension(:) :: Msn_eff_visc
    type(scalar_field), allocatable, dimension(:) :: Msn_int_mom_exch

    ! 2nd, 3rd, and 4th statistical moments for unclosed terms in volume filtered momentum equation
    type(scalar_field), allocatable, dimension(:), public :: stat_reynolds_stress
    type(scalar_field), allocatable, dimension(:), public :: stat_eff_visc
    type(scalar_field), allocatable, dimension(:), public :: stat_int_mom_exch

    !$acc declare create(Msn_reynolds_stress, Msn_eff_visc, Msn_int_mom_exch)

    !$acc declare create(stat_reynolds_stress, stat_eff_visc, stat_int_mom_exch)

contains

    subroutine s_initialize_statistics_module
        integer :: i

        @:ALLOCATE(Msn_reynolds_stress(1:4))
        do i = 1, 4
            @:ALLOCATE(Msn_reynolds_stress(i)%sf(0:m, 0:n, 0:p))
            @:ACC_SETUP_SFs(Msn_reynolds_stress(i))
        end do

        @:ALLOCATE(Msn_eff_visc(1:4))
        do i = 1, 4
            @:ALLOCATE(Msn_eff_visc(i)%sf(0:m, 0:n, 0:p))
            @:ACC_SETUP_SFs(Msn_eff_visc(i))
        end do

        @:ALLOCATE(Msn_int_mom_exch(1:4))
        do i = 1, 4
            @:ALLOCATE(Msn_int_mom_exch(i)%sf(0:m, 0:n, 0:p))
            @:ACC_SETUP_SFs(Msn_int_mom_exch(i))
        end do

        @:ALLOCATE(stat_reynolds_stress(1:4))
        do i = 1, 4
            @:ALLOCATE(stat_reynolds_stress(i)%sf(0:m, 0:n, 0:p))
            @:ACC_SETUP_SFs(stat_reynolds_stress(i))
        end do

        @:ALLOCATE(stat_eff_visc(1:4))
        do i = 1, 4
            @:ALLOCATE(stat_eff_visc(i)%sf(0:m, 0:n, 0:p))
            @:ACC_SETUP_SFs(stat_eff_visc(i))
        end do

        @:ALLOCATE(stat_int_mom_exch(1:4))
        do i = 1, 4
            @:ALLOCATE(stat_int_mom_exch(i)%sf(0:m, 0:n, 0:p))
            @:ACC_SETUP_SFs(stat_int_mom_exch(i))
        end do

    end subroutine s_initialize_statistics_module

    subroutine s_compute_statistics_momentum_unclosed_terms(n_step, reynolds_stress, eff_visc, int_mom_exch)
        type(scalar_field), intent(in) :: reynolds_stress 
        type(scalar_field), intent(in) :: eff_visc
        type(scalar_field), intent(in) :: int_mom_exch
        
        integer, intent(in) :: n_step
        real(wp) :: ns 

        ns = real(n_step, wp)

        ! update M1, M2, M3, M4
        call s_update_statistics(ns, reynolds_stress, Msn_reynolds_stress)
        call s_update_statistics(ns, eff_visc, Msn_eff_visc)
        call s_update_statistics(ns, int_mom_exch, Msn_int_mom_exch)

        ! compute 2nd, 3rd, 4th order statistical moments
        if (n_step > 3) then 
            call s_compute_234_order_statistics(ns, Msn_reynolds_stress, stat_reynolds_stress) 
            call s_compute_234_order_statistics(ns, Msn_eff_visc, stat_eff_visc) 
            call s_compute_234_order_statistics(ns, Msn_int_mom_exch, stat_int_mom_exch)  
        end if

    end subroutine s_compute_statistics_momentum_unclosed_terms

    subroutine s_update_statistics(ns, q_temp, Msn)
        type(scalar_field), intent(in) :: q_temp
        type(scalar_field), dimension(1:4), intent(inout) :: Msn

        real(wp), intent(in) :: ns
        real(wp) :: delta, delta_n, delta_n2, delta_f
        integer :: i, j, k

        !$acc parallel loop collapse(3) gang vector default(present) copyin(ns) private(delta, delta_n, delta_n2, delta_f)
        do i = 0, m 
            do j = 0, n 
                do k = 0, p
                    delta = q_temp%sf(i, j, k) - Msn(1)%sf(i, j, k)
                    delta_n = delta / ns
                    delta_n2 = delta_n**2
                    delta_f = delta * delta_n * (ns - 1._wp)

                    Msn(1)%sf(i, j, k) = Msn(1)%sf(i, j, k) + delta_n
                    Msn(4)%sf(i, j, k) = Msn(4)%sf(i, j, k) + delta_f * delta_n2 * (ns**2 - 3._wp*ns + 3._wp) + 6._wp * delta_n2 * Msn(2)%sf(i, j, k) - 4._wp * delta_n * Msn(3)%sf(i, j, k)
                    Msn(3)%sf(i, j, k) = Msn(3)%sf(i, j, k) + delta_f * delta_n * (ns - 2._wp) - 3._wp * delta_n * Msn(2)%sf(i, j, k)
                    Msn(2)%sf(i, j, k) = Msn(2)%sf(i, j, k) + delta_f
                end do 
            end do 
        end do
        
    end subroutine s_update_statistics

    subroutine s_compute_234_order_statistics(ns, Msn, q_stat)
        type(scalar_field), dimension(1:4), intent(in) :: Msn
        type(scalar_field), dimension(1:4), intent(inout) :: q_stat

        real(wp), intent(in) :: ns
        integer :: i, j, k

        !$acc parallel loop collapse(3) gang vector default(present) copyin(ns)
        do i = 0, m 
            do j = 0, n 
                do k = 0, p 
                    q_stat(1)%sf(i, j, k) = Msn(1)%sf(i, j, k)
                    q_stat(2)%sf(i, j, k) = Msn(2)%sf(i, j, k) / (ns - 1._wp)
                    q_stat(3)%sf(i, j, k) = sqrt(ns - 1._wp) / (ns - 2._wp) * ns * Msn(3)%sf(i, j, k) / (Msn(2)%sf(i, j, k)**1.5)
                    q_stat(4)%sf(i, j, k) = (ns - 1._wp) / ((ns - 2._wp) * (ns - 3._wp)) * ((ns + 1._wp) * (ns * Msn(4)%sf(i, j, k) / (Msn(2)%sf(i, j, k)**2) - 3._wp) + 6._wp)
                end do 
            end do 
        end do

    end subroutine s_compute_234_order_statistics

    subroutine s_finalize_statistics_module
        integer :: i, j
        do i = 1, 4
            @:DEALLOCATE(Msn_reynolds_stress(i)%sf)
        end do
        @:DEALLOCATE(Msn_reynolds_stress)

        do i = 1, 4
            @:DEALLOCATE(Msn_eff_visc(i)%sf)
        end do
        @:DEALLOCATE(Msn_eff_visc)

        do i = 1, 4
            @:DEALLOCATE(Msn_int_mom_exch(i)%sf)
        end do
        @:DEALLOCATE(Msn_int_mom_exch)

        do i = 1, 4
            @:DEALLOCATE(stat_reynolds_stress(i)%sf)
        end do
        @:DEALLOCATE(stat_reynolds_stress)

        do i = 1, 4
            @:DEALLOCATE(stat_eff_visc(i)%sf)
        end do
        @:DEALLOCATE(stat_eff_visc)

        do i = 1, 4
            @:DEALLOCATE(stat_int_mom_exch(i)%sf)
        end do
        @:DEALLOCATE(stat_int_mom_exch)

    end subroutine s_finalize_statistics_module

end module m_compute_statistics
