!>
!!@file
!!@brief Contains module m_rank_timing

#:include 'macros.fpp'

!> @brief Per-rank RHS compute-time imbalance diagnostic (Tier-2 calibration support).
module m_rank_timing

    use m_derived_types
    use m_global_parameters
    use m_mpi_common

    implicit none

    private
    public :: t_rank_rhs, s_report_rank_time

    real(wp) :: t_rank_rhs = 0._wp  !< accumulated RHS cpu_time this report interval (rank-local)

contains

    !> Reduce per-rank accumulated RHS time to a max/mean imbalance and print on rank 0, then reset the accumulator for the next
    !! interval.
    impure subroutine s_report_rank_time

        real(wp) :: t_max, t_sum, t_mean, imb
        integer  :: ierr

        if (.not. rank_time_wrt) return
#ifdef MFC_MPI
        call MPI_ALLREDUCE(t_rank_rhs, t_max, 1, mpi_p, MPI_MAX, MPI_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(t_rank_rhs, t_sum, 1, mpi_p, MPI_SUM, MPI_COMM_WORLD, ierr)
        t_mean = t_sum/real(num_procs, wp)
#else
        t_max = t_rank_rhs; t_mean = t_rank_rhs
#endif
        imb = t_max/max(t_mean, tiny(1._wp))
        if (proc_rank == 0) then
            print '(A,F8.3,A,ES12.5,A,ES12.5)', '[rank_time] imbalance(max/mean)= ', imb, '  t_max= ', t_max, '  t_mean= ', t_mean
        end if
        t_rank_rhs = 0._wp

    end subroutine s_report_rank_time

end module m_rank_timing
