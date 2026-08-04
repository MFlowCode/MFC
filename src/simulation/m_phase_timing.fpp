!>
!!@file
!!@brief Per-phase wall-time budget for the AMR step, gated on rank_time_wrt.
!!
!! WHY THIS EXISTS. Profiling MFC one layer at a time - GPU kernels (rocprofv3), MPI calls (a PMPI
!! shim), inter-operation gaps - gives each layer's share of a DIFFERENT denominator, so the pieces have
!! to be reconciled by inference and the reconciliation is easy to get wrong. Two changes were built and
!! reverted on such inferences before this existed. These brackets instead sum to the measured step-loop
!! wall and print the RESIDUAL, so what is unaccounted for is visible rather than assumed.
!!
!! Measured budget at 400^3/np=8/cap 32 when this landed: rhs 45.3%, regrid 13.0%, reflux 11.8%,
!! gather 10.5%, coarse base grid 5.5%, seam 3.1%, residual 7.1%. AMR machinery is ~40% of wall against
!! 5.5% for the base-grid physics.
!!
!! Each bracket does a device sync first, so a phase's time includes the GPU work it launched (host-only
!! timing would attribute launch cost to the phase and execution cost to whoever synced next).
#:include 'macros.fpp'

!> @brief Per-phase wall-time budget for the AMR step: brackets that sum to the measured step-loop wall and report the residual, so
!! time attribution needs no cross-layer inference. Gated on rank_time_wrt.
module m_phase_timing

    use m_derived_types
    use m_global_parameters
    use m_mpi_common

    implicit none

    private
    public :: s_phase_tic, s_phase_toc, s_phase_report, PH_N, PH_HALO, PH_GATHER, PH_GFILL, PH_SEAM, PH_RHS, PH_RK, PH_REFLUX, &
        & PH_REGRID, PH_L0, PH_COARSE

    integer, parameter          :: PH_HALO = 1     !< coarse cons halo exchange (hoisted, once per stage)
    integer, parameter          :: PH_GATHER = 2   !< per-block coarse-patch gather (P2P)
    integer, parameter          :: PH_GFILL = 3    !< ghost prolongation from the gathered patch
    integer, parameter          :: PH_SEAM = 4     !< fine-fine seam halo
    integer, parameter          :: PH_RHS = 5      !< fine-block s_compute_rhs
    integer, parameter          :: PH_RK = 6       !< fine-block RK update + relax + IB
    integer, parameter          :: PH_REFLUX = 7   !< reflux (p2p faces + apply)
    integer, parameter          :: PH_REGRID = 8   !< regrid / reassignment
    integer, parameter          :: PH_L0 = 9       !< L0 tile advance
    integer, parameter          :: PH_COARSE = 10  !< coarse (non-AMR) solver work
    integer, parameter          :: PH_N = 10
    character(len=8), parameter :: PH_NAME(PH_N) = [character(len=8)::'halo','gather', 'gfill', 'seam', 'rhs', 'rk', 'reflux', &
              & 'regrid', 'L0', 'coarse']

    real(wp)   :: acc(PH_N) = 0._wp
    integer(8) :: tic_c(PH_N) = 0
    integer    :: depth(PH_N) = 0
    real(wp)   :: t_wall0 = -1._wp

contains

    impure subroutine s_phase_tic(id)

        integer, intent(in) :: id
        integer(8)          :: c, rate

        if (.not. rank_time_wrt) return
        depth(id) = depth(id) + 1
        if (depth(id) > 1) return  ! outermost bracket only, so nesting cannot double count
        $:GPU_WAIT()
        call system_clock(c, rate)
        tic_c(id) = c
        if (t_wall0 < 0._wp) t_wall0 = real(c, wp)/real(rate, wp)

    end subroutine s_phase_tic

    impure subroutine s_phase_toc(id)

        integer, intent(in) :: id
        integer(8)          :: c, rate

        if (.not. rank_time_wrt) return
        if (depth(id) <= 0) then; depth(id) = 0; return; end if
        depth(id) = depth(id) - 1
        if (depth(id) > 0) return
        $:GPU_WAIT()
        call system_clock(c, rate)
        acc(id) = acc(id) + real(c - tic_c(id), wp)/real(rate, wp)

    end subroutine s_phase_toc

    !> Print the budget on rank 0. `wall` is the caller's measured step-loop wall so the RESIDUAL - the part no bracket covers - is
    !! reported instead of being silently absorbed.
    impure subroutine s_phase_report(wall)

        real(wp), intent(in) :: wall
        real(wp)             :: tot, gmax(PH_N), gsum(PH_N)
        integer              :: i, ierr

        if (.not. rank_time_wrt) return
        tot = sum(acc)
#ifdef MFC_MPI
        call MPI_ALLREDUCE(acc, gmax, PH_N, mpi_p, MPI_MAX, MPI_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(acc, gsum, PH_N, mpi_p, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
        gmax = acc; gsum = acc*real(num_procs, wp)
#endif
        if (proc_rank /= 0) return
        print '(A)', '[phase] PHASE BUDGET'
        print '(A,F10.3,A)', '[phase] step-loop wall = ', wall, ' s'
        print '(A)', '[phase] name        mean s   max s    % wall   imbalance(max/mean)'
        do i = 1, PH_N
            if (gsum(i) <= 0._wp) cycle
            print '(A,A8,F10.3,F9.3,F9.1,A,F8.3)', '[phase] ', PH_NAME(i), gsum(i)/real(num_procs, wp), gmax(i), &
                & 100._wp*(gsum(i)/real(num_procs, wp))/wall, '%', gmax(i)/max(gsum(i)/real(num_procs, wp), tiny(1._wp))
        end do
        print '(A,F10.3,F19.1,A)', '[phase] RESIDUAL', wall - sum(gsum)/real(num_procs, wp), &
            & 100._wp*(wall - sum(gsum)/real(num_procs, wp))/wall, '%'

    end subroutine s_phase_report

end module m_phase_timing
