!>
!! @file
!! @brief Contains module m_body_forces

#:include 'macros.fpp'

!> @brief Computes gravitational and user-defined body force source terms for the momentum equations
module m_body_forces

    use m_derived_types
    use m_global_parameters
    use m_variables_conversion
    use m_nvtx

    ! $:USE_GPU_MODULE()

    implicit none

    private
    public :: s_compute_body_forces_rhs, s_initialize_body_forces_module, s_finalize_body_forces_module

    real(wp), allocatable, dimension(:,:,:) :: rhoM
    $:GPU_DECLARE(create='[rhoM]')

contains

    !> Initialize the body forces module
    impure subroutine s_initialize_body_forces_module

        if (n > 0) then
            if (p > 0) then
                @:ALLOCATE(rhoM(-buff_size:buff_size + m, -buff_size:buff_size + n, -buff_size:buff_size + p))
            else
                @:ALLOCATE(rhoM(-buff_size:buff_size + m, -buff_size:buff_size + n, 0:0))
            end if
        else
            @:ALLOCATE(rhoM(-buff_size:buff_size + m, 0:0, 0:0))
        end if

    end subroutine s_initialize_body_forces_module

    !> Compute the acceleration at time t
    subroutine s_compute_acceleration(t)

        real(wp), intent(in) :: t

        #:for DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
            if (bf_${XYZ}$) then
                accel_bf(${DIR}$) = g_${XYZ}$ + k_${XYZ}$*sin(w_${XYZ}$*t - p_${XYZ}$)
            end if
        #:endfor

        $:GPU_UPDATE(device='[accel_bf]')

    end subroutine s_compute_acceleration

    !> Compute the mixture density at each cell center
    subroutine s_compute_mixture_density(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        integer                                             :: i, j, k, l  !< standard iterators

        $:GPU_PARALLEL_LOOP(private='[j, k, l]', collapse=3)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    rhoM(j, k, l) = 0._wp
                    do i = 1, num_fluids
                        rhoM(j, k, l) = rhoM(j, k, l) + q_cons_vf(eqn_idx%cont%beg + i - 1)%sf(j, k, l)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_mixture_density

    !> Compute the body force source terms for momentum and energy equations
    subroutine s_compute_body_forces_rhs(q_prim_vf, q_cons_vf, rhs_vf)

        type(scalar_field), dimension(sys_size), intent(in)    :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(in)    :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        integer                                                :: i, j, k, l  !< Loop variables

        call s_compute_acceleration(mytime)
        call s_compute_mixture_density(q_cons_vf)

        $:GPU_PARALLEL_LOOP(private='[i, j, k, l]', collapse=4)
        do i = eqn_idx%mom%beg, eqn_idx%E
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        rhs_vf(i)%sf(j, k, l) = 0._wp
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        if (bf_x) then  ! x-direction body forces

            $:GPU_PARALLEL_LOOP(private='[j, k, l]', collapse=3)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        rhs_vf(eqn_idx%mom%beg)%sf(j, k, l) = rhs_vf(eqn_idx%mom%beg)%sf(j, k, l) + rhoM(j, k, l)*accel_bf(1)
                        rhs_vf(eqn_idx%E)%sf(j, k, l) = rhs_vf(eqn_idx%E)%sf(j, k, l) + q_cons_vf(eqn_idx%mom%beg)%sf(j, k, &
                               & l)*accel_bf(1)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        if (bf_y) then  ! y-direction body forces

            $:GPU_PARALLEL_LOOP(private='[j, k, l]', collapse=3)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        rhs_vf(eqn_idx%mom%beg + 1)%sf(j, k, l) = rhs_vf(eqn_idx%mom%beg + 1)%sf(j, k, l) + rhoM(j, k, &
                               & l)*accel_bf(2)
                        rhs_vf(eqn_idx%E)%sf(j, k, l) = rhs_vf(eqn_idx%E)%sf(j, k, l) + q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k, &
                               & l)*accel_bf(2)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        if (bf_z) then  ! z-direction body forces

            $:GPU_PARALLEL_LOOP(private='[j, k, l]', collapse=3)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        rhs_vf(eqn_idx%mom%end)%sf(j, k, l) = rhs_vf(eqn_idx%mom%end)%sf(j, k, l) + rhoM(j, k, l)*accel_bf(3)
                        rhs_vf(eqn_idx%E)%sf(j, k, l) = rhs_vf(eqn_idx%E)%sf(j, k, l) + q_cons_vf(eqn_idx%mom%end)%sf(j, k, &
                               & l)*accel_bf(3)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_compute_body_forces_rhs

    !> Finalize the body forces module
    impure subroutine s_finalize_body_forces_module

        @:DEALLOCATE(rhoM)

    end subroutine s_finalize_body_forces_module

end module m_body_forces
