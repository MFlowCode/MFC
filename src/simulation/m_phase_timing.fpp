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
        & PH_RGHALO, PH_RGTAG, PH_RGCLUS, PH_RGSHAPE, PH_RGMIG, PH_RGBUILD, PH_REGRID, PH_L0, PH_COARSE
    public :: PH_RBGATH, PH_RBOVL, PH_RBPUSH, PH_RBWAIT, PH_RBALLOC, PH_RBUNPK
    public :: PH_SWAP, PH_RBOWN, PH_RBUPD, PH_RBPACK, PH_RBRSV
    public :: PH_RBSEAM, PH_RBPOST, PH_RBGEO, PH_RBSLOT, PH_RBTAIL
    public :: PH_RBSEND, PH_RBFLUSH, PH_RBXCHG, PH_RBREC, PH_RBTOPO
    public :: PH_PGALL, PH_PGSEND, PH_PGRECV
    public :: PH_RFP2P, PH_RFAPP, PH_RFRECV, PH_RFWAIT
    public :: PH_RESTR, PH_RGPART, PH_RGMOVE, PH_MGWAIT
    public :: PH_MGSLOT, PH_MGPACK, PH_MGUNPK, PH_MGPUSH
    public :: PH_GWPLAN, PH_GWPACK, PH_GWWAIT
    public :: PH_RSWAVE, PH_RSREST, PH_RSRFP
    public :: PH_CVTB, PH_BHALO
    public :: s_wait_tic, s_wait_toc, WT_GATHER, WT_PGATHER, WT_SEAM, WT_REFLUX, WT_RESTR, WT_REGRID, WT_HSLOT, WT_HSHELL, &
        & WT_HOWN, WT_HUNPK, WT_HFILL

    integer, parameter :: PH_HALO = 1     !< coarse cons halo exchange (hoisted, once per stage)
    integer, parameter :: PH_GATHER = 2   !< per-block coarse-patch gather (P2P)
    integer, parameter :: PH_GFILL = 3    !< ghost prolongation from the gathered patch
    integer, parameter :: PH_SEAM = 4     !< fine-fine seam halo
    integer, parameter :: PH_RHS = 5      !< fine-block s_compute_rhs
    integer, parameter :: PH_RK = 6       !< fine-block RK update + relax + IB
    integer, parameter :: PH_REFLUX = 7   !< reflux (p2p faces + apply)
    integer, parameter :: PH_REGRID = 8   !< regrid / reassignment
    integer, parameter :: PH_L0 = 9       !< L0 tile advance
    integer, parameter :: PH_COARSE = 10  !< coarse (non-AMR) solver work
    !> Regrid sub-phases. NESTED inside PH_REGRID, so they must NOT be summed with the top-level phases - the report prints them as
    !! a separate breakdown. Added because regrid measured 42% of wall at amr_regrid_int=2 while the optimisation effort was aimed
    !! at rhs (19%).
    integer, parameter :: PH_RGHALO = 11   !< coarse cons halo before tagging
    integer, parameter :: PH_RGTAG = 12    !< tag cells
    integer, parameter :: PH_RGCLUS = 13   !< cluster tags into boxes
    integer, parameter :: PH_RGSHAPE = 14  !< shape/nest/cap/unchanged checks
    integer, parameter :: PH_RGMIG = 15    !< stash + migrate old blocks
    integer, parameter :: PH_RGBUILD = 16  !< rebuild slots (per-box gather lives here)
    !> rg:build internals. The three candidate costs inside s_amr_regrid_rebuild_slots' per-box loop. Needed because cap32-vs-cap64
    !! scaling (cost ~ N^0.39) fits NONE of them alone: the O(N^2) old-box scan predicts 29x, the per-box rendezvous 5.4x, the
    !! volume-driven H2D copy the WRONG SIGN.
    integer, parameter :: PH_RBGATH = 17  !< (a) per-box collective gather
    integer, parameter :: PH_RBOVL = 18   !< (b) interpolate + O(old_np) overlap carry-forward
    integer, parameter :: PH_RBPUSH = 19  !< (c) per-box full-slot host->device update
    !> The MPI_WAITALL inside the REGRID-path gather only (gated by amr_rg_gather, since the same routine also serves the per-step
    !! path). rb:gath MINUS this is the gather's HOST work.
    integer, parameter :: PH_RBWAIT = 20
    !> Splitting the gather's HOST half. rb:wait was measured; the rest was attributed to the per-box allocate by CODE READING only,
    !! and the byte-proportional scaling fits the unpack equally well. These two brackets discriminate.
    integer, parameter :: PH_RBALLOC = 21  !< allocate/deallocate of rbuf,reqs,srank
    integer, parameter :: PH_RBUNPK = 22   !< post-wait unpack of rbuf into amr_cg
    !> Per-block grid-state reconfiguration: s_amr_swap_to_fine + the idwint push + s_amr_restore_coarse. TOP-LEVEL (parallel to
    !! rhs), not nested. This is what level-batching removes; it was previously unbracketed, and s_amr_restore_coarse used to sit
    !! inside PH_RHS, so the rhs bracket was charged half a swap pair.
    integer, parameter :: PH_SWAP = 23
    !> Splitting the gather's remaining HOST work (rb:gath minus wait/mem/unpk). The per-box allocate and the unpack were both
    !! REFUTED by measurement (0.002-0.010 s and 0.026-0.104 s), so these four cover what is actually left.
    integer, parameter :: PH_RBOWN = 24   !< owner's own-box local unpack (s_amr_unpack_patch)
    integer, parameter :: PH_RBUPD = 25   !< owner's per-box sys_size host->device push of amr_cg
    integer, parameter :: PH_RBPACK = 26  !< non-owner host pack loop into the send pool
    integer, parameter :: PH_RBRSV = 27   !< s_amr_gsnd_reserve - includes its force-drain MPI_WAITALL and pool resize
    !> Round 2: round 1 accounted for only 46.9%% of rb:gath and 70.0%% of rg:build, leaving 18.1%% of wall unexplained inside
    !! regrid. These five cover every remaining region on that path.
    integer, parameter :: PH_RBSEAM = 28  !< s_amr_build_seam_pairs (O(nblocks^2)) called from inside the gather
    integer, parameter :: PH_RBPOST = 29  !< the nsrc count + IRECV posting loop (per-(box,source) geometry)
    integer, parameter :: PH_RBGEO = 30   !< s_set_amr_fine_geometry - per box on EVERY rank
    integer, parameter :: PH_RBSLOT = 31  !< s_amr_alloc_slot (owner only)
    integer, parameter :: PH_RBTAIL = 32  !< post-loop tail: send flush, xchg reduce, reconcile, seam topology check
    !> Round 3. Round 2 closed rg:build to 100.0%% but left 11.9%% of wall inside rb:gath unexplained after SEVEN refuted code-read
    !! candidates; the only unbracketed code left in that routine is the non-owner ISEND and two scalar geometry calls. rb:tail
    !! (8.8%% of wall, imbalance 2.6) is split into its four collectives to separate barrier skew from work.
    integer, parameter :: PH_RBSEND = 33   !< the non-owner MPI_ISEND (rendezvous-sized: 1.5-3 MB)
    integer, parameter :: PH_RBFLUSH = 34  !< s_amr_gather_send_flush - one WAITALL over all deferred sends
    integer, parameter :: PH_RBXCHG = 35   !< s_amr_reduce_xchg_flag - MPI_ALLREDUCE, i.e. a barrier
    integer, parameter :: PH_RBREC = 36    !< s_amr_reconcile_slots
    integer, parameter :: PH_RBTOPO = 37   !< s_amr_check_seam_topology
    !> THE LEVEL>=2 PATH. `s_amr_gather_coarse_patch` returns at its FIRST branch for any block with level >= 2, into
    !! `s_amr_gather_from_parent` - so every rb:* bracket above instruments only the level-1 path, which is 64 of 224 boxes. The
    !! other 160 (71%%) were never measured. That is why EIGHT successive candidates each came back at ~0.
    integer, parameter :: PH_PGALL = 38   !< s_amr_gather_from_parent (the whole level>=2 path)
    integer, parameter :: PH_PGSEND = 39  !< parent owner: s_amr_gather_from_parent_field_cons (pack + send)
    integer, parameter :: PH_PGRECV = 40  !< block owner: s_amr_recv_parent_patch
    !> Decomposing REFLUX, whose per-call cost grows 19x between the 40-80 and 80-160 windows at a CONSTANT 64 level-1 blocks
    !! (imbalance 1.17, so it is real work, not waiting on a straggler). The owner already posts ISENDs + one WAITALL; each
    !! PARTICIPATING non-owner does 6 BLOCKING MPI_RECVs per block. rf:recv's CALL COUNT therefore measures how many blocks this
    !! rank participates in - if participation grows as the refined region spreads across ranks, that is the mechanism; if it is
    !! flat and ms/call grows instead, it is not.
    integer, parameter :: PH_RFP2P = 41   !< s_amr_p2p_reflux_faces (the whole exchange)
    integer, parameter :: PH_RFAPP = 42   !< s_amr_apply_reflux (local correction)
    integer, parameter :: PH_RFRECV = 43  !< non-owner blocking-RECV branch; CALL COUNT = participation
    integer, parameter :: PH_RFWAIT = 44  !< owner's MPI_WAITALL over its posted ISENDs
    !> The post-stage per-block restrict/reflux-to-parent chain (m_time_steppers, the reverse islot loop). Same per-box blocking P2P
    !! shape as PH_REFLUX, runs once per STEP over every block on every rank, and was entirely UNBRACKETED - it sits inside the
    !! 3.7-6.3%% residual. Its exit skew becomes the next step's entry skew, so it is the candidate for super-linear growth.
    integer, parameter :: PH_RESTR = 45
    !> The p4est 'complementarity' split of the regrid: PART decides the new partition (cluster, nest, assign owners) and moves
    !! NOTHING; MOVE is the data redistribution that follows. They are fused in s_amr_regrid_stash_migrate, so migration cost cannot
    !! be priced without this boundary - and every load-balance scheme in the literature needs it.
    integer, parameter :: PH_RGPART = 46
    integer, parameter :: PH_RGMOVE = 47
    !> The WAITALL inside the migration exchange. rg:move measures 103 MiB/s effective, far below intra-node bandwidth, so it is
    !! suspected wait-bound rather than volume-bound. If mg:wait is most of rg:move, cutting migration VOLUME (SFC hysteresis) will
    !! not convert to time.
    integer, parameter :: PH_MGWAIT = 48
    !> The rg:move work split (I4b pricing): slot = s_amr_alloc_slot_stash for received replicas (contains any store GROWTH - see
    !! s_amr_st_reserve), pack/unpk = the device pack/unpack kernels + their wire-slice transfers. mg:push is DEAD since the
    !! device-side migration (the per-received-slot full push it timed is deleted); the id stays so old budgets parse.
    integer, parameter :: PH_MGSLOT = 49
    integer, parameter :: PH_MGPACK = 50
    integer, parameter :: PH_MGUNPK = 51
    integer, parameter :: PH_MGPUSH = 52
    !> The stage-fill wave's internal split (I6 pricing): plan = the two replicated list walks, pack = the device pack kernels +
    !! their copyout transfers, wait = the single WAITALL. The residual of `gather` minus these three is recv/send posting + consume
    !! bookkeeping.
    integer, parameter :: PH_GWPLAN = 53
    integer, parameter :: PH_GWPACK = 54
    integer, parameter :: PH_GWWAIT = 55
    !> restr's internal split (the np16 rung made restr the largest inter-node growth): wave = s_amr_freg_wave (the F5b wire), rest
    !! = the restrict kernels, rfp = the level>=2 reflux-to-parent applies.
    integer, parameter :: PH_RSWAVE = 56
    integer, parameter :: PH_RSREST = 57
    integer, parameter :: PH_RSRFP = 58
    !> 2a: the batched cons->prim conversion over all owned fine blocks (s_amr_convert_prim_batch, once per stage)
    integer, parameter :: PH_CVTB = 59
    !> The BASE-GRID halo exchange (s_populate_variables_buffers) inside s_compute_rhs. It sits inside PH_COARSE, which is why the
    !! uniform (amr=F) arm reported 95%% 'coarse' and no communication at all - the one number needed to say how much of AMR's 31%%
    !! communication share is AMR's own rather than the solver's baseline.
    integer, parameter          :: PH_BHALO = 60
    integer, parameter          :: PH_N = 60
    character(len=8), parameter :: PH_NAME(PH_N) = [character(len=8)::'halo','gather', 'gfill', 'seam', 'rhs', 'rk', 'reflux', &
              & 'regrid', 'L0', 'coarse', 'rg:halo', 'rg:tag', 'rg:clus', 'rg:shape', 'rg:mig', 'rg:build', 'rb:gath', 'rb:ovl', &
              & 'rb:push', 'rb:wait', 'rb:mem', 'rb:unpk', 'swap', 'rb:own', 'rb:upd', 'rb:pack', 'rb:rsv', 'rb:seam', 'rb:post', &
              & 'rb:geo', 'rb:slot', 'rb:tail', 'rb:send', 'rb:flush', 'rb:xchg', 'rb:rec', 'rb:topo', 'pg:all', 'pg:send', &
              & 'pg:recv', 'rf:p2p', 'rf:app', 'rf:recv', 'rf:wait', 'restr', 'rg:part', 'rg:move', 'mg:wait', 'mg:slot', &
              & 'mg:pack', 'mg:unpk', 'mg:push', 'gw:plan', 'gw:pack', 'gw:wait', 'rs:wave', 'rs:rest', 'rs:rfp', 'cvt:bat', &
              & 'b:halo']

    !> The bracket-free MPI-wait table. Every s_phase_tic/toc drains the device first, so a bracket's `*:wait` row holds the GPU
    !! drain as well as the MPI wait and cannot split the excess into rank skew vs host work. These accumulate MPI_Wtime around ONLY
    !! the MPI_WAITALL / blocking MPI_RECV / MPI_SENDRECV calls, with no device sync and no MPI call anywhere on their path, keyed
    !! by the family whose [phase] bracket contains the site. The base-grid SENDRECV (m_mpi_common) serves three brackets, so its
    !! accumulator is snapshotted at their tic/toc instead; sr:other is whatever of it fell outside all three.
    integer, parameter :: WT_HALO = 1, WT_BHALO = 2, WT_RGHALO = 3, WT_SROTH = 4, WT_GATHER = 5, WT_PGATHER = 6, WT_SEAM = 7, &
        & WT_REFLUX = 8, WT_RESTR = 9, WT_REGRID = 10, WT_HSLOT = 11, WT_HSHELL = 12, WT_HOWN = 13, WT_HUNPK = 14, WT_HFILL = 15, &
        & WT_N = 15
    character(len=8), parameter :: WT_NAME(WT_N + 1) = [character(len=8)::'halo','b:halo', 'rg:halo', 'sr:other', 'gather', &
              & 'pgather', 'seam', 'reflux', 'restr', 'regrid', 'h:slot', 'h:shell', 'h:own', 'h:unpk', 'h:fill', 'TOTAL']
    integer, parameter :: SR_PH(3) = [PH_HALO, PH_BHALO, PH_RGHALO], SR_WT(3) = [WT_HALO, WT_BHALO, WT_RGHALO]
    real(dp)           :: wt(WT_N) = 0._dp, wt_t0 = 0._dp, sr_t0(3) = 0._dp  !< MPI_Wtime is double; wp may be single
    integer(8)         :: wtc(WT_N) = 0, sr_n0(3) = 0
    real(wp)           :: acc(PH_N) = 0._wp
    !> Entry count per phase. Time alone cannot distinguish "this region is slow" from "this region runs far more often than
    !! assumed"; eight code-read attributions were refuted by brackets before this column existed, and the ninth candidate had no
    !! code left to blame. ms/call is what tells the two apart.
    integer(8) :: ncall(PH_N) = 0
    integer(8) :: tic_c(PH_N) = 0
    integer    :: depth(PH_N) = 0
    real(wp)   :: t_wall0 = -1._wp

contains

    impure subroutine s_phase_tic(id)

        integer, intent(in) :: id
        integer(8)          :: c, rate
        integer             :: i

        if (.not. rank_time_wrt) return
        depth(id) = depth(id) + 1
        if (depth(id) > 1) return  ! outermost bracket only, so nesting cannot double count
        do i = 1, 3
            if (id == SR_PH(i)) then; sr_t0(i) = mpi_sr_wait; sr_n0(i) = mpi_sr_calls; end if
        end do
        $:GPU_WAIT()
        call system_clock(c, rate)
        tic_c(id) = c
        ncall(id) = ncall(id) + 1
        if (t_wall0 < 0._wp) t_wall0 = real(c, wp)/real(rate, wp)

    end subroutine s_phase_tic

    impure subroutine s_phase_toc(id)

        integer, intent(in) :: id
        integer(8)          :: c, rate
        integer             :: i

        if (.not. rank_time_wrt) return
        if (depth(id) <= 0) then; depth(id) = 0; return; end if
        depth(id) = depth(id) - 1
        if (depth(id) > 0) return
        $:GPU_WAIT()
        call system_clock(c, rate)
        acc(id) = acc(id) + real(c - tic_c(id), wp)/real(rate, wp)
        do i = 1, 3
            if (id == SR_PH(i)) then
                wt(SR_WT(i)) = wt(SR_WT(i)) + (mpi_sr_wait - sr_t0(i)); wtc(SR_WT(i)) = wtc(SR_WT(i)) + (mpi_sr_calls - sr_n0(i))
            end if
        end do

    end subroutine s_phase_toc

    !> Bracket ONE MPI wait/recv/sendrecv call: s_wait_tic() immediately before it, s_wait_toc(family) immediately after. Waits do
    !! not nest, so one timestamp suffices. No device sync, no MPI call, and nothing at all when rank_time_wrt is off.
    impure subroutine s_wait_tic()

#ifdef MFC_MPI
        if (rank_time_wrt) wt_t0 = MPI_Wtime()
#endif

    end subroutine s_wait_tic

    impure subroutine s_wait_toc(id)

        integer, intent(in) :: id

#ifdef MFC_MPI
        if (.not. rank_time_wrt) return
        wt(id) = wt(id) + (MPI_Wtime() - wt_t0); wtc(id) = wtc(id) + 1
#endif

    end subroutine s_wait_toc

    !> Print the budget on rank 0. `wall` is the caller's measured step-loop wall so the RESIDUAL - the part no bracket covers - is
    !! reported instead of being silently absorbed.
    impure subroutine s_phase_report(wall)

        real(wp), intent(in) :: wall
        real(wp)             :: tot, gmax(PH_N), gsum(PH_N)
        integer(8)           :: gcall(PH_N)
        integer              :: i, ierr, ip
        !> Per-rank times for the phases whose IMBALANCE moves with simulation time. mean/max cannot say WHICH rank is slow or
        !! whether it is the one holding more work, which is what the rhs skew (1.09 -> 2.90 between the 80- and 160-step windows)
        !! actually needs. The regrid rows split the one-rank regrid straggler the `[mpiwait] regrid` row cannot (it sums four
        !! WAITALL sites): which of migrate / rebuild, and inside them pack vs wait vs gather-wait vs the flag barrier, each rank
        !! spent its regrid seconds in.
        integer, parameter :: NPR = 12
        integer, parameter :: PR_ID(NPR) = [PH_RHS, PH_REFLUX, PH_GATHER, PH_SEAM, PH_RGMIG, PH_MGPACK, PH_MGWAIT, PH_RGBUILD, &
                                    & PH_RBGATH, PH_RBWAIT, PH_PGRECV, PH_RBXCHG]
        real(wp), allocatable :: prank(:,:)
        real(dp), allocatable :: wrank(:,:)
        integer(8)            :: wcall(WT_N + 1)

        if (.not. rank_time_wrt) return
        tot = sum(acc)
#ifdef MFC_MPI
        call MPI_ALLREDUCE(acc, gmax, PH_N, mpi_p, MPI_MAX, MPI_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(acc, gsum, PH_N, mpi_p, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(ncall, gcall, PH_N, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
        gmax = acc; gsum = acc*real(num_procs, wp); gcall = ncall*int(num_procs, 8)
#endif
        allocate (prank(0:num_procs - 1,NPR))
        do i = 1, NPR
#ifdef MFC_MPI
            call MPI_GATHER(acc(PR_ID(i)), 1, mpi_p, prank(0, i), 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
#else
            prank(0, i) = acc(PR_ID(i))
#endif
        end do
        if (proc_rank == 0) then
            do i = 1, NPR
                write (*, '(A,A8,A)', advance='no') '[phase-rank] ', PH_NAME(PR_ID(i)), ' :'
                do ip = 0, num_procs - 1
                    write (*, '(F10.2)', advance='no') prank(ip, i)
                end do
                write (*, '(A)') ''
            end do
        end if
        deallocate (prank)
        wt(WT_SROTH) = mpi_sr_wait - sum(wt(SR_WT)); wtc(WT_SROTH) = mpi_sr_calls - sum(wtc(SR_WT))
        allocate (wrank(0:num_procs - 1,WT_N + 1))
        do i = 1, WT_N
#ifdef MFC_MPI
            call MPI_GATHER(wt(i), 1, MPI_DOUBLE_PRECISION, wrank(0, i), 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#else
            wrank(0, i) = wt(i)
#endif
        end do
#ifdef MFC_MPI
        call MPI_ALLREDUCE(wtc, wcall(1:WT_N), WT_N, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
        wcall(1:WT_N) = wtc
#endif
        wcall(WT_N + 1) = sum(wcall(1:WT_N))
        if (proc_rank == 0) then
            wrank(:,WT_N + 1) = sum(wrank(:,1:WT_N), 2)
            print '(A)', '[mpiwait] MPI WAIT (inside MPI_WAITALL / MPI_RECV / MPI_SENDRECV only; no device sync on this path);'
            print '(A)', &
                & '[mpiwait] the h:* rows are HOST brackets around the per-block gather consume (ledger 81), same clock, no MPI'
            print '(A)', '[mpiwait] name       mean s    max s    min s  calls/rank    ms/call  per-rank s'
            do i = 1, WT_N + 1
                if (wcall(i) == 0) cycle
                write (*, '(A,A8,3F9.3,I12,F11.4,A)', advance='no') '[mpiwait] ', WT_NAME(i), sum(wrank(:,i))/real(num_procs, &
                       & dp), maxval(wrank(:,i)), minval(wrank(:,i)), wcall(i)/int(num_procs, 8), 1000._dp*sum(wrank(:, &
                       & i))/real(wcall(i), dp), ' :'
                do ip = 0, num_procs - 1
                    write (*, '(F9.3)', advance='no') wrank(ip, i)
                end do
                write (*, '(A)') ''
            end do
        end if
        deallocate (wrank)
        if (proc_rank /= 0) return
        print '(A)', '[phase] PHASE BUDGET'
        print '(A,F10.3,A)', '[phase] step-loop wall = ', wall, ' s'
        print '(A)', '[phase] name        mean s   max s    % wall   imbalance  calls/rank    ms/call'
        do i = 1, PH_N
            if (gsum(i) <= 0._wp) cycle
            print '(A,A8,F10.3,F9.3,F9.1,A,F8.3,I12,F11.4)', '[phase] ', PH_NAME(i), gsum(i)/real(num_procs, wp), gmax(i), &
                & 100._wp*(gsum(i)/real(num_procs, wp))/wall, '%', gmax(i)/max(gsum(i)/real(num_procs, wp), tiny(1._wp)), &
                & gcall(i)/int(num_procs, 8), 1000._wp*(gsum(i)/real(num_procs, wp))/max(real(gcall(i)/int(num_procs, 8), wp), &
                & 1._wp)
        end do
        print '(A,F10.3,F19.1,A)', '[phase] RESIDUAL', wall - sum(gsum)/real(num_procs, wp), &
            & 100._wp*(wall - sum(gsum)/real(num_procs, wp))/wall, '%'

    end subroutine s_phase_report

end module m_phase_timing
