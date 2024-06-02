#:include 'macros.fpp'

module m_body_forces

    ! Dependencies =============================================================

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_variables_conversion

    use m_nvtx

#ifdef MFC_OpenACC
    use openacc
#endif
    ! ==========================================================================

    implicit none

    private; public :: s_compute_body_forces_rhs, &
 s_initialize_body_forces_module, &
 s_finalize_body_forces_module

#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :), rhoM)
    !$acc declare link(rhoM)
#else
    real(kind(0d0)), allocatable, dimension(:, :, :) :: rhoM
    !$acc declare create(rhoM)
#endif

contains

    !> This subroutine inializes the module global array of mixture
    !! densities in each grid cell
    subroutine s_initialize_body_forces_module()

        ! Simulation is at least 2D
        if (n > 0) then
            ! Simulation is 3D
            if (p > 0) then
                @:ALLOCATE_GLOBAL (rhoM(-buff_size:buff_size + m, &
                    -buff_size:buff_size + n, &
                    -buff_size:buff_size + p))
                ! Simulation is 2D
            else
                @:ALLOCATE_GLOBAL (rhoM(-buff_size:buff_size + m, &
                    -buff_size:buff_size + n, &
                    0:0))
            end if
            ! Simulation is 1D
        else
            @:ALLOCATE_GLOBAL (rhoM(-buff_size:buff_size + m, &
                0:0, &
                0:0))
        end if

    end subroutine s_initialize_body_forces_module

    !> This subroutine computes the acceleration at time t
    subroutine s_compute_acceleration(t)

        real(kind(0d0)) :: t

        if (m > 0) then
            accel_bf(1) = g_x + k_x*sin(w_x*t - p_x)
            if (n > 0) then
                accel_bf(2) = g_y + k_y*sin(w_y*t - p_y)
                if (p > 0) then
                    accel_bf(3) = g_z + k_z*sin(w_z*t - p_z)
                end if
            end if
        end if

        !$acc update device(accel_bf)

    end subroutine s_compute_acceleration

    !> This subroutine calculates teh mixture density at each cell
    !! center
    subroutine s_compute_mixture_density(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        integer :: i, j, k, l !< standard iterators

        !$acc parallel loop collapse(3) gang vector default(present)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    rhoM(j, k, l) = 0d0
                    do i = 1, num_fluids
                        rhoM(j, k, l) = rhoM(j, k, l) + &
                                        q_cons_vf(contxb + i - 1)%sf(j, k, l)
                    end do
                end do
            end do
        end do

    end subroutine s_compute_mixture_density

    !> This subroutine calculates the source term due to body forces
    !! so the system can be advanced in time
    subroutine s_compute_body_forces_rhs(q_cons_vf, q_prim_vf, rhs_vf)

        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(INOUT) :: rhs_vf

        integer :: i, j, k, l !< Loop variables

        call s_compute_acceleration(mytime)
        call s_compute_mixture_density(q_cons_vf)

        !$acc parallel loop collapse(4) gang vector default(present)
        do i = momxb, E_idx
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        rhs_vf(i)%sf(j, k, l) = 0d0
                    end do
                end do
            end do
        end do

        if (bf_x) then ! x-direction body forces

            !$acc parallel loop collapse(3) gang vector default(present)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        rhs_vf(momxb)%sf(j, k, l) = rhs_vf(momxb)%sf(j, k, l) + &
                                                    rhoM(j, k, l)*accel_bf(1)
                        rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) + &
                                                    q_cons_vf(momxb)%sf(j, k, l)*accel_bf(1)
                    end do
                end do
            end do
        end if

        if (bf_y) then ! y-direction body forces

            !$acc parallel loop collapse(3) gang vector default(present)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        rhs_vf(momxb + 1)%sf(j, k, l) = rhs_vf(momxb + 1)%sf(j, k, l) + &
                                                        rhoM(j, k, l)*accel_bf(2)
                        rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) + &
                                                    q_cons_vf(momxb + 1)%sf(j, k, l)*accel_bf(2)
                    end do
                end do
            end do
        end if

        if (bf_z) then ! z-direction body forces

            !$acc parallel loop collapse(3) gang vector default(present)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        rhs_vf(momxe)%sf(j, k, l) = rhs_vf(momxe)%sf(j, k, l) + &
                                                    (rhoM(j, k, l))*accel_bf(3)
                        rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) + &
                                                    q_cons_vf(momxe)%sf(j, k, l)*accel_bf(3)
                    end do
                end do
            end do

        end if

    end subroutine s_compute_body_forces_rhs

    subroutine s_finalize_body_forces_module()

        @:DEALLOCATE_GLOBAL(rhoM)

    end subroutine s_finalize_body_forces_module

end module m_body_forces
