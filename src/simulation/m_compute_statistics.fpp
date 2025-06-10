#:include 'macros.fpp'

module m_compute_statistics
    use m_derived_types 

    use m_global_parameters

    use m_mpi_proxy 

    implicit none

    private; public :: s_initialize_statistics_module, s_finalize_statistics_module, s_compute_s_order_statistics

    type(scalar_field), allocatable, dimension(:) :: xnbar_stat

    type(scalar_field), allocatable, dimension(:) :: delta_stat

    type(vector_field), allocatable, dimension(:) :: Msn_stat

    !$acc declare create(xnbar_stat, delta_stat, Msn_stat)

contains

    subroutine s_initialize_statistics_module
        integer :: i, j
        @:ALLOCATE(xnbar_stat(1:3))
        do i = 1, 3
            @:ALLOCATE(xnbar_stat(i)%sf(0:m, 0:n, 0:p))
            @:ACC_SETUP_SFs(xnbar_stat(i))
        end do

        @:ALLOCATE(delta_stat(1:3))
        do i = 1, 3
            @:ALLOCATE(delta_stat(i)%sf(0:m, 0:n, 0:p))
            @:ACC_SETUP_SFs(delta_stat(i))
        end do

        @:ALLOCATE(Msn_stat(1:num_dims))
        do i = 1, 3
            @:ALLOCATE(Msn_stat(i)%vf(2:4))
        end do
        do i = 1, 3
            do j = 2, 4
                @:ALLOCATE(Msn_stat(i)%vf(j)%sf(0:m, 0:n, 0:p))
            end do
            @:ACC_SETUP_VFs(Msn_stat(i))
        end do

    end subroutine s_initialize_statistics_module

    subroutine s_compute_s_order_statistics(q_temp, n_step, s_order_stat, id)
        type(scalar_field), intent(in) :: q_temp
        integer, intent(in) :: n_step
        type(scalar_field), dimension(2:4), intent(inout) :: s_order_stat
        integer, intent(in) :: id
        real(wp) :: ns
        integer :: i, j, k, ii

        ns = real(n_step, wp)

        if (n_step == 1) then
            !$acc parallel loop collapse(3) gang vector default(present)
            do i = 0, m 
                do j = 0, n
                    do k = 0, p
                        xnbar_stat(id)%sf(i, j, k) = q_temp%sf(i, j, k)
                        Msn_stat(id)%vf(2)%sf(i, j, k) = 0.0_wp
                        Msn_stat(id)%vf(3)%sf(i, j, k) = 0.0_wp
                        Msn_stat(id)%vf(4)%sf(i, j, k) = 0.0_wp
                        s_order_stat(2)%sf(i, j, k) = 0.0_wp
                        s_order_stat(3)%sf(i, j, k) = 0.0_wp
                        s_order_stat(4)%sf(i, j, k) = 0.0_wp
                    end do 
                end do
            end do
        else 
            !$acc parallel loop collapse(3) gang vector default(present) copyin(ns)
            do i = 0, m 
                do j = 0, n
                    do k = 0, p
                        delta_stat(id)%sf(i, j, k) = q_temp%sf(i, j, k) - xnbar_stat(id)%sf(i, j, k)

                        xnbar_stat(id)%sf(i, j, k) = xnbar_stat(id)%sf(i, j, k) + delta_stat(id)%sf(i, j, k)/ns

                        Msn_stat(id)%vf(4)%sf(i, j, k) = Msn_stat(id)%vf(4)%sf(i, j, k) & 
                                                + (delta_stat(id)%sf(i, j, k)**4)*(ns - 1.0_wp)*(ns**2 - 3.0_wp*ns + 3.0_wp)/(ns**3) &
                                                + 6.0_wp*(delta_stat(id)%sf(i, j, k)**2)*Msn_stat(id)%vf(2)%sf(i, j, k)/(ns**2) &
                                                - 4.0_wp*delta_stat(id)%sf(i, j, k)*Msn_stat(id)%vf(3)%sf(i, j, k)/ns

                        Msn_stat(id)%vf(3)%sf(i, j, k) = Msn_stat(id)%vf(3)%sf(i, j, k) & 
                                                + (delta_stat(id)%sf(i, j, k)**3)*(ns - 1.0_wp)*(ns - 2.0_wp)/(ns**2) & 
                                                - 3.0_wp*delta_stat(id)%sf(i, j, k)*Msn_stat(id)%vf(2)%sf(i, j, k)/ns

                        Msn_stat(id)%vf(2)%sf(i, j, k) = Msn_stat(id)%vf(2)%sf(i, j, k) &
                                                + (delta_stat(id)%sf(i, j, k)**2)*(ns - 1.0_wp)/ns

                        s_order_stat(2)%sf(i, j, k) = Msn_stat(id)%vf(2)%sf(i, j, k)/(ns - 1.0_wp)

                        s_order_stat(3)%sf(i, j, k) = sqrt(ns)*Msn_stat(id)%vf(3)%sf(i, j, k)/(Msn_stat(id)%vf(2)%sf(i, j, k)**1.5_wp)

                        s_order_stat(4)%sf(i, j, k) = ns*Msn_stat(id)%vf(4)%sf(i, j, k)/(Msn_stat(id)%vf(2)%sf(i, j, k)**2) - 3.0_wp
                    end do 
                end do
            end do
        end if

    end subroutine s_compute_s_order_statistics

    subroutine s_finalize_statistics_module
        integer :: i, j
        do i = 1, 3
            @:DEALLOCATE(xnbar_stat(i)%sf)
        end do
        @:DEALLOCATE(xnbar_stat)

        do i = 1, 3
            @:DEALLOCATE(delta_stat(i)%sf)
        end do
        @:DEALLOCATE(delta_stat)

        do i = 1, 3
            do j = 2, 4
                @:DEALLOCATE(Msn_stat(i)%vf(j)%sf)
            end do
            @:DEALLOCATE(Msn_stat(i)%vf)
        end do
        @:DEALLOCATE(Msn_stat)
    end subroutine s_finalize_statistics_module

end module m_compute_statistics