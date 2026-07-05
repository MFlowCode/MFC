!>
!!@file
!!@brief Contains module m_rank_timing

#:include 'macros.fpp'

!> @brief Per-rank compute-time (RHS + phase-change relaxation) imbalance diagnostic (Tier-2 calibration support).
module m_rank_timing

    use m_derived_types
    use m_global_parameters
    use m_mpi_common

    implicit none

    private
    public :: s_rank_time_tic, s_rank_time_toc, s_report_rank_time

    !> accumulated per-rank wall time over the timed regions this report interval (rank-local)
    real(wp) :: t_rank_compute = 0._wp
    !> system_clock tick captured at the most recent tic
    integer(8) :: tic_count
    !> tic/toc nesting depth: only the outermost pair measures, so the AMR fine advance can bracket its compute segments while the
    !! s_compute_rhs pair inside becomes a no-op (no double counting)
    integer :: tic_depth = 0

contains

    !> Start a per-rank wall-clock timing region. The device sync first ensures any prior GPU work has completed, so the interval
    !! that follows measures only the timed region.
    impure subroutine s_rank_time_tic

        if (.not. rank_time_wrt) return
        tic_depth = tic_depth + 1
        if (tic_depth > 1) return
        $:GPU_WAIT()
        call system_clock(tic_count)

    end subroutine s_rank_time_tic

    !> End a per-rank wall-clock timing region and accumulate its wall duration. The device sync first forces the timed GPU kernels
    !! to complete, so the wall interval reflects real compute time (cpu_time would capture only host-side launch overhead on the
    !! GPU backend).
    impure subroutine s_rank_time_toc

        integer(8) :: toc_count, rate

        if (.not. rank_time_wrt) return
        ! unmatched toc (caller bug): clamp instead of underflowing into an uninitialized tic_count
        if (tic_depth <= 0) then
            tic_depth = 0
            return
        end if
        tic_depth = tic_depth - 1
        if (tic_depth > 0) return
        $:GPU_WAIT()
        call system_clock(toc_count, rate)
        t_rank_compute = t_rank_compute + real(toc_count - tic_count, wp)/real(rate, wp)

    end subroutine s_rank_time_toc

    !> Reduce per-rank accumulated compute time to a max/mean imbalance and print on rank 0, then reset the accumulator for the next
    !! interval.
    impure subroutine s_report_rank_time

        real(wp) :: t_max, t_sum, t_mean, imb
        integer  :: ierr

        if (.not. rank_time_wrt) return
#ifdef MFC_MPI
        call MPI_ALLREDUCE(t_rank_compute, t_max, 1, mpi_p, MPI_MAX, MPI_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(t_rank_compute, t_sum, 1, mpi_p, MPI_SUM, MPI_COMM_WORLD, ierr)
        t_mean = t_sum/real(num_procs, wp)
#else
        t_max = t_rank_compute; t_mean = t_rank_compute
#endif
        imb = t_max/max(t_mean, tiny(1._wp))
        if (proc_rank == 0) then
            print '(A,F8.3,A,ES12.5,A,ES12.5)', '[rank_time] imbalance(max/mean)= ', imb, '  t_max= ', t_max, '  t_mean= ', t_mean
        end if
        t_rank_compute = 0._wp

    end subroutine s_report_rank_time

end module m_rank_timing
