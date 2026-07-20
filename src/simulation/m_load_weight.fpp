!>
!!@file
!!@brief Contains module m_load_weight

#:include 'macros.fpp'

!> @brief Diagnostic per-cell compute-cost weight field + per-rank load-imbalance metric.
module m_load_weight

    use m_derived_types
    use m_global_parameters
    use m_constants, only: K_bub, K_ib, K_pc
    use m_mpi_proxy
    use m_mpi_common
    use m_active_box, only: ab_active, ab_x, ab_y, ab_z
    use m_bubbles_EL, only: q_beta
    use m_ibm, only: ib_markers
    use m_phase_change, only: pc_iter_count

    implicit none

    private
    public :: s_initialize_load_weight_module, s_finalize_load_weight_module, s_compute_load_weight, s_report_load_imbalance, &
        & load_weight

    type(scalar_field) :: load_weight  !< per-cell modeled compute-cost weight
    $:GPU_DECLARE(create='[load_weight]')

    ! Relative cost coefficients K_bub/K_ib/K_pc come from m_constants (shared with the AMR block-owner cost weighting).

contains

    impure subroutine s_initialize_load_weight_module

        if (.not. load_weight_wrt .and. .not. sfc_partition_wrt) return
        @:ALLOCATE(load_weight%sf(idwint(1)%beg:idwint(1)%end, idwint(2)%beg:idwint(2)%end, idwint(3)%beg:idwint(3)%end))
        @:ACC_SETUP_SFs(load_weight)

    end subroutine s_initialize_load_weight_module

    impure subroutine s_finalize_load_weight_module

        if (.not. load_weight_wrt .and. .not. sfc_partition_wrt) return
        @:DEALLOCATE(load_weight%sf)

    end subroutine s_finalize_load_weight_module

    !> Base cost 1 everywhere; cells outside the active box get 0 (frozen).
    impure subroutine s_compute_load_weight()

        integer :: j, k, l
        integer :: jlo, jhi, klo, khi, llo, lhi

        if (.not. load_weight_wrt .and. .not. sfc_partition_wrt) return

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

        ! EE bubbles: NO weight contribution. Calibration (rank_time_wrt vs load_weight_wrt,
        ! bubblescreen -n 4) showed measured RHS-time imbalance is flat (~1.02-1.03) across
        ! vf0 = 4e-5 -> 0.1 (2500x void range), while a K_bub*void term manufactured up to
        ! 1.59x modeled imbalance. The EE bubble source (s_compute_bubble_EE_source, inside
        ! s_compute_rhs) is a roughly-uniform per-cell term whose cost does NOT scale with void
        ! magnitude -- so it is not a load-imbalance driver and gets no weight. (K_bub is an EL
        ! per-bubble-ODE constant; it does not apply to EE.) If a future EE regime -- QBMM,
        ! adv_n number-density, very high resolution -- shows real EE imbalance, re-add with an
        ! EE-specific coefficient calibrated against measured rank_time.

        ! EL bubble contributor: K_bub * per-cell bubble void fraction.
        ! q_beta(1)%sf holds the liquid volume fraction (1 - alpha_bub) after s_smear_voidfraction;
        ! 1 - q_beta(1)%sf gives the smeared bubble void fraction as a per-cell count proxy.
        if (bubbles_lagrange) then
            $:GPU_PARALLEL_LOOP(collapse=3)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        if (j >= jlo .and. j <= jhi .and. k >= klo .and. k <= khi .and. l >= llo .and. l <= lhi) then
                            load_weight%sf(j, k, l) = load_weight%sf(j, k, l) + K_bub*(1.0_wp - real(q_beta(1)%sf(j, k, l), wp))
                        end if
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        ! IB contributor: K_ib per IB-marked interior cell.
        if (ib) then
            $:GPU_PARALLEL_LOOP(collapse=3)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        if (j >= jlo .and. j <= jhi .and. k >= klo .and. k <= khi .and. l >= llo .and. l <= lhi) then
                            if (ib_markers%sf(j, k, l) /= 0) then
                                load_weight%sf(j, k, l) = load_weight%sf(j, k, l) + real(K_ib, stp)
                            end if
                        end if
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        ! Phase-change contributor: K_pc * per-cell Newton-iteration count from s_infinite_relaxation_k.
        if (relax) then
            $:GPU_PARALLEL_LOOP(collapse=3)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        if (j >= jlo .and. j <= jhi .and. k >= klo .and. k <= khi .and. l >= llo .and. l <= lhi) then
                            load_weight%sf(j, k, l) = load_weight%sf(j, k, l) + real(K_pc*real(pc_iter_count(j, k, l), wp), stp)
                        end if
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

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
