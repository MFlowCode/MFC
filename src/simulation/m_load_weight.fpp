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
    use m_active_box, only: ab_active, ab_x, ab_y, ab_z

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

    !> Base cost 1 everywhere; cells outside the active box get 0 (frozen).
    impure subroutine s_compute_load_weight(q_cons_vf, q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf, q_prim_vf
        integer                                             :: j, k, l
        integer                                             :: jlo, jhi, klo, khi, llo, lhi

        if (.not. load_weight_wrt) return

        if (ab_active) then
            jlo = ab_x%beg; jhi = ab_x%end
            klo = ab_y%beg; khi = ab_y%end
            llo = ab_z%beg; lhi = ab_z%end
        else
            jlo = 0; jhi = m; klo = 0; khi = n; llo = 0; lhi = p
        end if

        $:GPU_PARALLEL_LOOP(collapse=3)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    if (j >= jlo .and. j <= jhi .and. k >= klo .and. k <= khi .and. l >= llo .and. l <= lhi) then
                        load_weight%sf(j, k, l) = 1._wp
                    else
                        load_weight%sf(j, k, l) = 0._wp
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_load_weight

    !> Print per-rank load-imbalance metric: max/mean of per-rank weight sums.
    impure subroutine s_report_load_imbalance

        real(wp) :: w_local, w_sum, w_max, w_mean, imbalance
        integer  :: j, k, l, ierr

        if (.not. load_weight_wrt) return
        w_local = 0._wp
        $:GPU_PARALLEL_LOOP(collapse=3, reduction='[[w_local]]', reductionOp='[+]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    w_local = w_local + real(load_weight%sf(j, k, l), wp)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
#ifdef MFC_MPI
        call MPI_ALLREDUCE(w_local, w_sum, 1, mpi_p, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(w_local, w_max, 1, mpi_p, MPI_MAX, MPI_COMM_WORLD, ierr)
#else
        w_sum = w_local; w_max = w_local
#endif
        w_mean = w_sum/real(num_procs, wp)
        imbalance = w_max/max(w_mean, tiny(1._wp))
        if (proc_rank == 0) then
            print '(A,F8.3,A,ES12.5,A,ES12.5)', '[load_weight] imbalance(max/mean)=', imbalance, '  w_max=', w_max, '  w_mean=', &
                & w_mean
        end if

    end subroutine s_report_load_imbalance

end module m_load_weight
