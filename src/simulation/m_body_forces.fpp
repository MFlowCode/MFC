#:include 'macros.fpp'

module m_body_forces

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_variables_conversion

    use m_nvtx

#ifdef MFC_OpenACC
    use openacc
#endif

    implicit none

    private; 
    public :: s_compute_body_forces_rhs, &
              s_initialize_body_forces_module, &
              s_finalize_body_forces_module

    real(wp), allocatable, dimension(:, :, :) :: rhoM
    $:GPU_DECLARE(create='[rhoM]')

contains

    !> This subroutine initializes the module global array of mixture
    !! densities in each grid cell
    impure subroutine s_initialize_body_forces_module

        ! Simulation is at least 2D
        if (n > 0) then
            ! Simulation is 3D
            if (p > 0) then
                @:ALLOCATE (rhoM(-buff_size:buff_size + m, &
                    -buff_size:buff_size + n, &
                    -buff_size:buff_size + p))
                ! Simulation is 2D
            else
                @:ALLOCATE (rhoM(-buff_size:buff_size + m, &
                    -buff_size:buff_size + n, &
                    0:0))
            end if
            ! Simulation is 1D
        else
            @:ALLOCATE (rhoM(-buff_size:buff_size + m, &
                0:0, &
                0:0))
        end if

    end subroutine s_initialize_body_forces_module

    !> This subroutine computes the acceleration at time t
    subroutine s_compute_acceleration(t)

        real(wp), intent(in) :: t

        #:for DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
            if (bf_${XYZ}$) then
                accel_bf(${DIR}$) = g_${XYZ}$+k_${XYZ}$*sin(w_${XYZ}$*t - p_${XYZ}$)
            end if
        #:endfor

        $:GPU_UPDATE(device='[accel_bf]')

    end subroutine s_compute_acceleration

    !> This subroutine calculates the mixture density at each cell
    !! center
    !! param q_cons_vf Conservative variable
    subroutine s_compute_mixture_density(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        integer :: i, j, k, l !< standard iterators

        $:GPU_PARALLEL_LOOP(collapse=3)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    rhoM(j, k, l) = 0._wp
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
    !! @param q_cons_vf Conservative variables
    !! @param q_prim_vf Primitive variables
    subroutine s_compute_body_forces_rhs(q_prim_vf, q_cons_vf, rhs_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf

        integer :: i, j, k, l !< Loop variables

        call s_compute_acceleration(mytime)
        call s_compute_mixture_density(q_cons_vf)

        $:GPU_PARALLEL_LOOP(collapse=4)
        do i = momxb, E_idx
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        rhs_vf(i)%sf(j, k, l) = 0._wp
                    end do
                end do
            end do
        end do

        if (bf_x) then ! x-direction body forces

            $:GPU_PARALLEL_LOOP(collapse=3)
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

            $:GPU_PARALLEL_LOOP(collapse=3)
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

            $:GPU_PARALLEL_LOOP(collapse=3)
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

    impure subroutine s_finalize_body_forces_module

        @:DEALLOCATE(rhoM)

    end subroutine s_finalize_body_forces_module

end module m_body_forces
