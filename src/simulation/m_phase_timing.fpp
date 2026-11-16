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
        & PH_RG_OWN, PH_RG_MIG, PH_RG_BUILD, PH_RG_REC, s_mem_note, s_mem_note_fields, MEM_STORE, MEM_SCRATCH, MEM_PATCH, &
        & MEM_REG, MEM_IB, MEM_IBHOST

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

    !> AMR allocations by group, noted by each owner before it frees them (s_mem_note) and printed with the phase table.
    integer, parameter          :: MEM_STORE = 1    !< flat block store (cons + RK stage copy)
    integer, parameter          :: MEM_SCRATCH = 2  !< batched-advance scratch (prim, rhs) and the block bridge
    integer, parameter          :: MEM_PATCH = 3    !< coarse patch pool and its subcycling time cache
    integer, parameter          :: MEM_REG = 4      !< flux registers
    integer, parameter          :: MEM_IB = 5       !< fine IB ghost-point park
    integer, parameter          :: MEM_IBHOST = 6   !< per-slot fine IB marker copies (host only)
    integer, parameter          :: MEM_N = 6
    character(len=8), parameter :: MEM_NAME(MEM_N + 1) = [character(len=8)::'store','scratch', 'patch', 'register', 'ib park', &
              & 'ib marks', 'total']
    real(wp) :: mem_mib(MEM_N) = 0._wp

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

        real(wp)   :: gmax(PH_N), gsum(PH_N), mean, mmax(MEM_N + 1), msum(MEM_N + 1)
        integer(8) :: gcall(PH_N)
        integer    :: i, ierr

        if (.not. rank_time_wrt) return
#ifdef MFC_MPI
        gmax = acc
        call s_mpi_allreduce_array_max(gmax, PH_N)
        call MPI_ALLREDUCE(acc, gsum, PH_N, mpi_p, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(ncall, gcall, PH_N, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
        mmax = [mem_mib, sum(mem_mib)]  ! last entry: the rank's total
        call MPI_ALLREDUCE(mmax, msum, MEM_N + 1, mpi_p, MPI_SUM, MPI_COMM_WORLD, ierr)
        call s_mpi_allreduce_array_max(mmax, MEM_N + 1)
#else
        gmax = acc; gsum = acc; gcall = ncall; mmax = [mem_mib, sum(mem_mib)]; msum = mmax
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
        if (msum(MEM_N + 1) <= 0._wp) return
        print '(A)', '[amr-mem] group       mean MiB     max MiB  (per rank, at the end of the run)'
        do i = 1, MEM_N + 1
            print '(A,A8,2F12.1)', '[amr-mem] ', MEM_NAME(i), msum(i)/real(num_procs, wp), mmax(i)
        end do
        print '(A,F8.2,A)', '[amr-mem] solver work arrays are allocated ', &
            & product(real(idwbuff_alloc%end - idwbuff_alloc%beg + 1, wp))/product(real(idwbuff%end - idwbuff%beg + 1, wp)), &
            & 'x the coarse extent (largest block or batch slab)'

    end subroutine s_phase_report

    !> Add bytes to group id of the end-of-run AMR memory table.
    impure subroutine s_mem_note(id, bytes)

        integer, intent(in)         :: id
        integer(kind=8), intent(in) :: bytes

        mem_mib(id) = mem_mib(id) + real(bytes, wp)/1024._wp**2

    end subroutine s_mem_note

    !> s_mem_note for the fields of f that are allocated.
    impure subroutine s_mem_note_fields(id, f)

        integer, intent(in)            :: id
        type(scalar_field), intent(in) :: f(:)
        integer                        :: i

        do i = 1, size(f)
            if (associated(f(i)%sf)) call s_mem_note(id, storage_size(f(i)%sf, kind=8)/8_8*size(f(i)%sf, kind=8))
        end do

    end subroutine s_mem_note_fields

end module m_phase_timing
