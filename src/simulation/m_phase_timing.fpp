!>
!!@file
!!@brief Contains module m_phase_timing

#:include 'macros.fpp'

!> @brief Per-phase wall-time budget of the AMR step, reported when rank_time_wrt is on. Brackets are disjoint (regrid sub-phases
!! nest inside regrid) and drain the device first, so a phase's time includes the GPU work it launched.
module m_phase_timing

    use m_global_parameters
    use m_mpi_common

    implicit none

    private
    public :: s_phase_tic, s_phase_toc, s_phase_step_tic, s_phase_step_toc, s_phase_report, PH_HALO, PH_GATHER, PH_GFILL, &
        & PH_SEAM, PH_RHS, PH_RK, PH_REFLUX, PH_REGRID, PH_COARSE, PH_SWAP, PH_RESTR, PH_LERP, PH_RG_HALO, PH_RG_TAG, PH_RG_CLUS, &
        & PH_RG_OWN, PH_RG_MIG, PH_RG_BUILD, PH_RG_REC

    integer, parameter          :: PH_HALO = 1       !< coarse cons halo exchange (once per stage)
    integer, parameter          :: PH_GATHER = 2     !< fill waves: coarse-patch gather into every block
    integer, parameter          :: PH_GFILL = 3      !< ghost prolongation from the gathered patch
    integer, parameter          :: PH_SEAM = 4       !< fine-fine seam halo
    integer, parameter          :: PH_RHS = 5        !< fine-block RHS
    integer, parameter          :: PH_RK = 6         !< fine-block RK update + IB
    integer, parameter          :: PH_REFLUX = 7     !< level-1 reflux (faces wave + apply)
    integer, parameter          :: PH_REGRID = 8     !< regrid
    integer, parameter          :: PH_COARSE = 9     !< coarse (base-grid) solver work
    integer, parameter          :: PH_SWAP = 10      !< per-batch grid-state swap and restore
    integer, parameter          :: PH_RESTR = 11     !< post-stage restrict + reflux-to-parent fold
    integer, parameter          :: PH_LERP = 12      !< subcycled coarse time interpolation for the level-1 fills
    integer, parameter          :: PH_NTOP = 12      !< the disjoint brackets above; the regrid sub-phases below nest in PH_REGRID
    integer, parameter          :: PH_RG_HALO = 13   !< coarse halo + host copy for tagging
    integer, parameter          :: PH_RG_TAG = 14    !< tagging
    integer, parameter          :: PH_RG_CLUS = 15   !< clustering, shaping, nesting and box checks
    integer, parameter          :: PH_RG_OWN = 16    !< stash and owner assignment
    integer, parameter          :: PH_RG_MIG = 17    !< fine-state migration
    integer, parameter          :: PH_RG_BUILD = 18  !< new-block build (fills, prolongation, overlap copy)
    integer, parameter          :: PH_RG_REC = 19    !< slot reconcile, IB setup, seam check
    integer, parameter          :: PH_N = 19
    character(len=8), parameter :: PH_NAME(PH_N) = [character(len=8)::'halo','gather', 'gfill', 'seam', 'rhs', 'rk', 'reflux', &
              & 'regrid', 'coarse', 'swap', 'restr', 'lerp', ' rg.halo', ' rg.tag', ' rg.clus', ' rg.own', ' rg.mig', &
              & ' rg.build', ' rg.rec']

    real(wp)   :: acc(PH_N) = 0._wp, wall = 0._wp  !< wall is the accumulated step-loop time the phases are budgeted against
    integer(8) :: ncall(PH_N) = 0, tic_c(PH_N) = 0, wall_c = 0
    integer    :: depth(PH_N) = 0                  !< a re-entered bracket (shared routine) measures only its outermost pair

contains

    impure subroutine s_phase_tic(id)

        integer, intent(in) :: id

        if (.not. rank_time_wrt) return
        depth(id) = depth(id) + 1
        if (depth(id) > 1) return
        $:GPU_WAIT()
        call system_clock(tic_c(id))

    end subroutine s_phase_tic

    impure subroutine s_phase_toc(id)

        integer, intent(in) :: id
        integer(8)          :: c, rate

        if (.not. rank_time_wrt) return
        depth(id) = depth(id) - 1
        if (depth(id) > 0) return
        $:GPU_WAIT()
        call system_clock(c, rate)
        acc(id) = acc(id) + real(c - tic_c(id), wp)/real(rate, wp)
        ncall(id) = ncall(id) + 1

    end subroutine s_phase_toc

    impure subroutine s_phase_step_tic()

        if (rank_time_wrt) call system_clock(wall_c)

    end subroutine s_phase_step_tic

    impure subroutine s_phase_step_toc()

        integer(8) :: c, rate

        if (.not. rank_time_wrt) return
        call system_clock(c, rate)
        wall = wall + real(c - wall_c, wp)/real(rate, wp)

    end subroutine s_phase_step_toc

    !> Mean and max over ranks per phase, share of the step-loop wall, imbalance (max/mean) and calls; the residual is the wall
    !! minus the sum of the means.
    impure subroutine s_phase_report()

        real(wp)   :: gmax(PH_N), gsum(PH_N), mean
        integer(8) :: gcall(PH_N)
        integer    :: i, ierr

        if (.not. rank_time_wrt) return
#ifdef MFC_MPI
        gmax = acc
        call s_mpi_allreduce_array_max(gmax, PH_N)
        call MPI_ALLREDUCE(acc, gsum, PH_N, mpi_p, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(ncall, gcall, PH_N, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
        gmax = acc; gsum = acc; gcall = ncall
#endif
        if (proc_rank /= 0) return
        print '(A,F10.3,A)', '[phase] step-loop wall = ', wall, ' s'
        print '(A)', '[phase] name        mean s   max s    % wall   imbalance  calls/rank    ms/call'
        do i = 1, PH_N
            if (gsum(i) <= 0._wp) cycle
            mean = gsum(i)/real(num_procs, wp)
            print '(A,A8,F10.3,F9.3,F9.1,A,F8.3,I12,F11.4)', '[phase] ', PH_NAME(i), mean, gmax(i), 100._wp*mean/wall, '%', &
                & gmax(i)/mean, gcall(i)/int(num_procs, 8), 1000._wp*mean/max(real(gcall(i)/int(num_procs, 8), wp), 1._wp)
        end do
        print '(A,F10.3,F19.1,A)', '[phase] RESIDUAL', wall - sum(gsum(1:PH_NTOP))/real(num_procs, wp), &
            & 100._wp*(wall - sum(gsum(1:PH_NTOP))/real(num_procs, wp))/wall, '%'

    end subroutine s_phase_report

end module m_phase_timing
