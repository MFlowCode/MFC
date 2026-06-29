!>
!!@file
!!@brief Contains module m_load_weight

#:include 'macros.fpp'

!> @brief Diagnostic per-cell compute-cost weight field + per-rank load-imbalance metric.
module m_load_weight

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_mpi_common

    implicit none

    private
    public :: s_initialize_load_weight_module, s_finalize_load_weight_module, s_compute_load_weight, s_report_load_imbalance, &
        & load_weight

    type(scalar_field) :: load_weight  !< per-cell modeled compute-cost weight
    $:GPU_DECLARE(create='[load_weight]')

    ! Relative cost coefficients (calibrated in Task 7; base RHS cell = 1).
    real(wp), parameter :: K_bub = 50._wp  !< per local bubble (stiff adaptive ODE); refined by validation
    real(wp), parameter :: K_ib = 2._wp    !< per IB ghost/interior cell
    real(wp), parameter :: K_pc = 3._wp    !< per phase-change Newton iteration

contains

    impure subroutine s_initialize_load_weight_module

        if (.not. load_weight_wrt) return
        @:ALLOCATE(load_weight%sf(idwint(1)%beg:idwint(1)%end, idwint(2)%beg:idwint(2)%end, idwint(3)%beg:idwint(3)%end))
        @:ACC_SETUP_SFs(load_weight)

    end subroutine s_initialize_load_weight_module

    impure subroutine s_finalize_load_weight_module

        if (.not. load_weight_wrt) return
        @:DEALLOCATE(load_weight%sf)

    end subroutine s_finalize_load_weight_module

    !> Task 1 stub: base cost 1 everywhere. Contributors added in Tasks 2,4,5,6.
    impure subroutine s_compute_load_weight(q_cons_vf, q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf, q_prim_vf
        integer                                             :: j, k, l

        if (.not. load_weight_wrt) return
        $:GPU_PARALLEL_LOOP(collapse=3)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    load_weight%sf(j, k, l) = 1._wp
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_load_weight

    !> Task 1 stub: no-op. Filled in Task 3.
    impure subroutine s_report_load_imbalance

    end subroutine s_report_load_imbalance

end module m_load_weight
