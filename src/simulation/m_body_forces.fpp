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
    public :: s_compute_body_forces_rhs, s_compute_synthetic_forces_rhs, s_initialize_body_forces_module, &
        & s_finalize_body_forces_module

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

    subroutine s_compute_turbulent_forces_rhs(q_prim_vf, q_cons_vf, rhs_vf)

        type(scalar_field), dimension(sys_size), intent(in)    :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(in)    :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        integer                                                :: i, j, k, l, q, n_mode
        real(wp)                                               :: gauss_x, two_over_Lx
        real(wp)                                               :: u_prime, rho_local
        real(wp), dimension(num_dims)                          :: u_local, force, pos
        real(wp)                                               :: phase_arg

        ! G_bar for the 1D Gaussian exp(-pi/2 * (2x/Lx)^2) on [-Lx/2, Lx/2] Analytical: G_bar = erf(sqrt(pi/2)) / sqrt(pi/2) \approx
        ! 0.7602
        real(wp), parameter :: G_BAR = 0.760173_wp

        do turb_idx = 1, num_turbulent_sources
            $:GPU_PARALLEL_LOOP(private='[i, j, k, l, q, n_mode, x_pos, y_pos, z_pos, gauss_x, u_prime, force, rho_local, &
                                & u_local, phase_arg]', collapse=3)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        ! position relative to forcing-zone center
                        #:for 'X', 'ID', 'IND' in [('x', 1, 'j'), ('y', 2, 'k'), ('z', 3, 'l')]
                            if (${ID}$ <= num_dims) pos(${ID}$) = ${X}$_cc(${IND}$) - turb_pos(turb_idx, ${ID}$)
                        #:endfor

                        ! Skip cells outside the forcing box
                        if (abs(pos(1)) > 0.5_wp*synth_L(turb_idx, 1) .or. abs(pos(2)) > 0.5_wp*synth_L(turb_idx, &
                            & 2) .or. abs(pos(3)) > 0.5_wp*synth_L(turb_idx, 3)) then
                            do i = eqn_idx%mom%beg, eqn_idx%mom%beg + num_dims - 1
                                rhs_vf(i)%sf(j, k, l) = 0._wp
                            end do
                            cycle
                        end if

                        ! Gaussian envelope
                        gauss_scaler = 1._wp
                        do q = 1, num_dims
                            gauss_scaler = gauss_scaler*exp(-0.5_wp*pi*(2._wp*pos/(synth_L(turb_idx, q)))**2)
                        end do

                        ! Synthetic fluctuation
                        u_prime = 0._wp
                        do n_mode = 1, synth_n_modes
                            phase_arg = synth_k(n_mode)*(x_cc(j) - synth_U_inf*t_step) + synth_phase(n_mode)
                            u_prime = u_prime + synth_amp(n_mode)*cos(phase_arg)
                        end do

                        ! Local density and x-velocity from primitive vars
                        rho_local = q_prim_vf(eqn_idx%cont%beg)%sf(j, k, l)
                        do q = 1, num_dims
                            u_local(q) = q_prim_vf(eqn_idx%mom%beg + q - 1)%sf(j, k, l)
                        end do

                        ! Volume force per Tangermann & Klein Eq. (2) F_syn = rho * u' / T * G / G_bar (x-component only) T = Lx /
                        ! U_inf
                        do q = 1, num_dims
                            force(i) = rho_local*u_prime*(synth_U_inf/synth_Lx)*gauss_scaler/(G_BAR**num_dims)
                        end do

                        ! Apply to momentum (x) and energy equations
                        rhs_vf(eqn_idx%mom%beg)%sf(j, k, l) = force_x

                        ! Energy source: power injected = F . u  (only x-component nonzero)
                        rhs_vf(eqn_idx%E)%sf(j, k, l) = force_x*u_local

                        ! Continuity and y/z momentum get zero
                        rhs_vf(eqn_idx%cont%beg)%sf(j, k, l) = 0._wp
                        if (n > 0) rhs_vf(eqn_idx%mom%beg + 1)%sf(j, k, l) = 0._wp
                        if (p > 0) rhs_vf(eqn_idx%mom%beg + 2)%sf(j, k, l) = 0._wp
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do

    end subroutine s_compute_synthetic_forces_rhs

    !> Finalize the body forces module
    impure subroutine s_finalize_body_forces_module

        @:DEALLOCATE(rhoM)

    end subroutine s_finalize_body_forces_module

end module m_body_forces
