!>
!!@file
!!@brief Per-phase wall-time budget for the AMR step, gated on rank_time_wrt.
!!
!! Profiling one layer at a time (GPU kernels, MPI calls, inter-operation gaps) gives each layer's share of a different
!! denominator, so the pieces have to be reconciled by inference. These brackets instead sum to the measured step-loop wall and
!! print the residual, so what is unaccounted for is visible rather than assumed.
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
    !> Regrid sub-phases. Nested inside PH_REGRID, so they must not be summed with the top-level phases; the report prints them as a
    !! separate breakdown.
    integer, parameter :: PH_RGHALO = 11   !< coarse cons halo before tagging
    integer, parameter :: PH_RGTAG = 12    !< tag cells
    integer, parameter :: PH_RGCLUS = 13   !< cluster tags into boxes
    integer, parameter :: PH_RGSHAPE = 14  !< shape/nest/cap/unchanged checks
    integer, parameter :: PH_RGMIG = 15    !< stash + migrate old blocks
    integer, parameter :: PH_RGBUILD = 16  !< rebuild slots (per-box gather lives here)
    !> rg:build internals: the three main costs inside the per-box loop of s_amr_regrid_rebuild_slots.
    integer, parameter :: PH_RBGATH = 17  !< (a) per-box collective gather
    integer, parameter :: PH_RBOVL = 18   !< (b) interpolate + O(old_np) overlap carry-forward
    integer, parameter :: PH_RBPUSH = 19  !< (c) per-box full-slot host->device update
    !> The MPI_WAITALL inside the regrid-path gather only (gated by amr_rg_gather, since the same routine also serves the per-step
    !! path). rb:gath minus this is the gather's host work.
    integer, parameter :: PH_RBWAIT = 20
    !> Split of the gather's host half.
    integer, parameter :: PH_RBALLOC = 21  !< allocate/deallocate of rbuf,reqs,srank
    integer, parameter :: PH_RBUNPK = 22   !< post-wait unpack of rbuf into amr_cg
    !> Per-block grid-state reconfiguration: s_amr_swap_to_fine + the idwint push + s_amr_restore_coarse. Top-level (parallel to
    !! rhs), not nested.
    integer, parameter :: PH_SWAP = 23
    !> Split of the gather's remaining host work (rb:gath minus wait/mem/unpk).
    integer, parameter :: PH_RBOWN = 24   !< owner's own-box local unpack (s_amr_unpack_patch)
    integer, parameter :: PH_RBUPD = 25   !< owner's per-box sys_size host->device push of amr_cg
    integer, parameter :: PH_RBPACK = 26  !< non-owner host pack loop into the send pool
    integer, parameter :: PH_RBRSV = 27   !< s_amr_gsnd_reserve, including its force-drain MPI_WAITALL and pool resize
    !> The remaining regions on the rg:build path.
    integer, parameter :: PH_RBSEAM = 28  !< s_amr_build_seam_pairs (O(nblocks^2)) called from inside the gather
    integer, parameter :: PH_RBPOST = 29  !< the nsrc count + IRECV posting loop (per-(box,source) geometry)
    integer, parameter :: PH_RBGEO = 30   !< s_set_amr_fine_geometry, per box on every rank
    integer, parameter :: PH_RBSLOT = 31  !< s_amr_alloc_slot (owner only)
    integer, parameter :: PH_RBTAIL = 32  !< post-loop tail: send flush, xchg reduce, reconcile, seam topology check
    !> The non-owner ISEND inside the gather, and rb:tail split into its four collectives to separate barrier skew from work.
    integer, parameter :: PH_RBSEND = 33   !< the non-owner MPI_ISEND (rendezvous-sized)
    integer, parameter :: PH_RBFLUSH = 34  !< s_amr_gather_send_flush, one WAITALL over all deferred sends
    integer, parameter :: PH_RBXCHG = 35   !< s_amr_reduce_xchg_flag, an MPI_ALLREDUCE, i.e. a barrier
    integer, parameter :: PH_RBREC = 36    !< s_amr_reconcile_slots
    integer, parameter :: PH_RBTOPO = 37   !< s_amr_check_seam_topology
    !> The level>=2 path. `s_amr_gather_coarse_patch` returns at its first branch for any block with level >= 2, into
    !! `s_amr_gather_from_parent`, so every rb:* bracket above instruments only the level-1 path; these cover the rest.
    integer, parameter :: PH_PGALL = 38   !< s_amr_gather_from_parent (the whole level>=2 path)
    integer, parameter :: PH_PGSEND = 39  !< parent owner: s_amr_gather_from_parent_field_cons (pack + send)
    integer, parameter :: PH_PGRECV = 40  !< block owner: s_amr_recv_parent_patch
    !> Reflux decomposition. The owner posts ISENDs + one WAITALL; each participating non-owner does blocking MPI_RECVs per block,
    !! so rf:recv's call count measures how many blocks this rank participates in.
    integer, parameter :: PH_RFP2P = 41   !< s_amr_p2p_reflux_faces (the whole exchange)
    integer, parameter :: PH_RFAPP = 42   !< s_amr_apply_reflux (local correction)
    integer, parameter :: PH_RFRECV = 43  !< non-owner blocking-RECV branch; call count = participation
    integer, parameter :: PH_RFWAIT = 44  !< owner's MPI_WAITALL over its posted ISENDs
    !> The post-stage per-block restrict/reflux-to-parent chain (m_time_steppers, the reverse islot loop). Same per-box blocking P2P
    !! shape as PH_REFLUX; runs once per step over every block on every rank. Its exit skew becomes the next step's entry skew.
    integer, parameter :: PH_RESTR = 45
    !> The partition/move split of the regrid: PART decides the new partition (cluster, nest, assign owners) and moves nothing; MOVE
    !! is the data redistribution that follows. Both live inside s_amr_regrid_stash_migrate, so migration cost cannot be priced
    !! without this boundary.
    integer, parameter :: PH_RGPART = 46
    integer, parameter :: PH_RGMOVE = 47
    !> The WAITALL inside the migration exchange; separates wait from volume inside rg:move.
    integer, parameter :: PH_MGWAIT = 48
    !> The rg:move work split: slot = s_amr_alloc_slot_stash for received replicas (contains any store growth, see
    !! s_amr_st_reserve), pack/unpk = the device pack/unpack kernels + their wire-slice transfers. mg:push is unused (migration is
    !! device-side, so there is no per-received-slot full push); the id stays so existing budget parsers keep working.
    integer, parameter :: PH_MGSLOT = 49
    integer, parameter :: PH_MGPACK = 50
    integer, parameter :: PH_MGUNPK = 51
    integer, parameter :: PH_MGPUSH = 52
    !> The stage-fill wave's internal split: plan = the two replicated list walks, pack = the device pack kernels + their copyout
    !! transfers, wait = the single WAITALL. The residual of `gather` minus these three is recv/send posting + consume bookkeeping.
    integer, parameter :: PH_GWPLAN = 53
    integer, parameter :: PH_GWPACK = 54
    integer, parameter :: PH_GWWAIT = 55
    !> restr's internal split. wave = the level>=2 flux-register exchange (s_amr_freg_wave, called unconditionally off the subcycle
    !! path), rest = the restrict kernels, rfp = the level>=2 reflux-to-parent applies. Naming trap: two of these three rows are
    !! reflux, not restriction (wave is the level>=2 flux-register wire and rfp is the level>=2 Berger-Colella apply); only rest is
    !! the restrict kernels. AMR's true reflux total is the `reflux` row plus these two.
    integer, parameter :: PH_RSWAVE = 56
    integer, parameter :: PH_RSREST = 57
    integer, parameter :: PH_RSRFP = 58
    !> The batched cons->prim conversion over all owned fine blocks (s_amr_convert_prim_batch, once per stage)
    integer, parameter :: PH_CVTB = 59
    !> The base-grid halo exchange (s_populate_variables_buffers) inside s_compute_rhs. It sits inside PH_COARSE; this row separates
    !! the solver's baseline communication from AMR's own.
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
    !! drain as well as the MPI wait and cannot split the excess into rank skew vs host work. These accumulate MPI_Wtime around only
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
    !! assumed"; ms/call is what tells the two apart.
    integer(8) :: ncall(PH_N) = 0
    integer(8) :: tic_c(PH_N) = 0
    integer    :: depth(PH_N) = 0
    !> Observed nesting, so the budget validates itself instead of trusting a hand-kept list of top-level rows. open_ids is the
    !! stack of brackets open on this rank (each id at most once, by the depth guard, so PH_N bounds it); tier_lo/tier_hi are the
    !! shallowest and deepest depth each phase was ever opened at. A phase with tier_lo /= tier_hi is opened both inside and outside
    !! another bracket (a shared routine) and cannot be summed at either level. n_interleave counts a toc that was not the innermost
    !! open bracket, n_orphan a toc with no tic: either one means some row double counts.
    integer  :: open_ids(PH_N) = 0, n_open = 0
    integer  :: tier_lo(PH_N) = huge(1), tier_hi(PH_N) = 0
    integer  :: n_interleave = 0, n_orphan = 0
    real(wp) :: t_wall0 = -1._wp

contains

    impure subroutine s_phase_tic(id)

        integer, intent(in) :: id
        integer(8)          :: c, rate
        integer             :: i

        if (.not. rank_time_wrt) return
        depth(id) = depth(id) + 1
        if (depth(id) > 1) return  ! outermost bracket only, so nesting cannot double count
        n_open = n_open + 1; open_ids(n_open) = id
        tier_lo(id) = min(tier_lo(id), n_open); tier_hi(id) = max(tier_hi(id), n_open)
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
        if (depth(id) <= 0) then; depth(id) = 0; n_orphan = n_orphan + 1; return; end if
        depth(id) = depth(id) - 1
        if (depth(id) > 0) return
        if (open_ids(n_open) /= id) n_interleave = n_interleave + 1
        do i = n_open, 1, -1
            if (open_ids(i) == id) then
                open_ids(i:n_open - 1) = open_ids(i + 1:n_open); n_open = n_open - 1; exit
            end if
        end do
        $:GPU_WAIT()
        call system_clock(c, rate)
        acc(id) = acc(id) + real(c - tic_c(id), wp)/real(rate, wp)
        do i = 1, 3
            if (id == SR_PH(i)) then
                wt(SR_WT(i)) = wt(SR_WT(i)) + (mpi_sr_wait - sr_t0(i)); wtc(SR_WT(i)) = wtc(SR_WT(i)) + (mpi_sr_calls - sr_n0(i))
            end if
        end do

    end subroutine s_phase_toc

    !> Bracket one MPI wait/recv/sendrecv call: s_wait_tic() immediately before it, s_wait_toc(family) immediately after. Waits do
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

    !> Print the budget on rank 0. `wall` is the caller's measured step-loop wall so the residual (the part no bracket covers) is
    !! reported instead of being silently absorbed.
    impure subroutine s_phase_report(wall)

        real(wp), intent(in) :: wall
        real(wp)             :: tot, gmax(PH_N), gsum(PH_N), t1sum
        integer              :: gtlo(PH_N), gthi(PH_N), gbad(2), nbad(2)
        character(len=8)     :: tl
        integer(8)           :: gcall(PH_N)
        integer              :: i, ierr, ip
        !> Per-rank times for every phase: the wall is one rank's serial chain, so pricing a change needs every phase per rank, not
        !! the mean/max pair. Rows whose global sum is zero are not printed. Under rank_time_wrt only, like the rest of this report.
        real(wp), allocatable :: prank(:,:)
        real(dp), allocatable :: wrank(:,:)
        integer(8)            :: wcall(WT_N + 1)

        if (.not. rank_time_wrt) return
        tot = sum(acc)
#ifdef MFC_MPI
        call MPI_ALLREDUCE(acc, gmax, PH_N, mpi_p, MPI_MAX, MPI_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(acc, gsum, PH_N, mpi_p, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(ncall, gcall, PH_N, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(tier_lo, gtlo, PH_N, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(tier_hi, gthi, PH_N, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
        nbad = [n_interleave, n_orphan]
        call MPI_ALLREDUCE(nbad, gbad, 2, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
        gmax = acc; gsum = acc*real(num_procs, wp); gcall = ncall*int(num_procs, 8)
        gtlo = tier_lo; gthi = tier_hi; gbad = [n_interleave, n_orphan]
#endif
        allocate (prank(0:num_procs - 1,PH_N))
        do i = 1, PH_N
#ifdef MFC_MPI
            call MPI_GATHER(acc(i), 1, mpi_p, prank(0, i), 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
#else
            prank(0, i) = acc(i)
#endif
        end do
        if (proc_rank == 0) then
            do i = 1, PH_N
                if (gsum(i) <= 0._wp) cycle
                write (*, '(A,A8,A)', advance='no') '[phase-rank] ', PH_NAME(i), ' :'
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
            print '(A)', '[mpiwait] the h:* rows are HOST brackets around the per-block gather consume, same clock, no MPI'
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
        print '(A)', '[phase] name        mean s   max s    % wall   imbalance  calls/rank    ms/call  tier'
        do i = 1, PH_N
            if (gsum(i) <= 0._wp) cycle
            if (gtlo(i) == gthi(i)) then
                write (tl, '(A,I0)') 'T', gtlo(i)
            else
                write (tl, '(A,I0,A,I0)') 'T', gtlo(i), '-', gthi(i)
            end if
            print '(A,A8,F10.3,F9.3,F9.1,A,F8.3,I12,F11.4,2X,A)', '[phase] ', PH_NAME(i), gsum(i)/real(num_procs, wp), gmax(i), &
                & 100._wp*(gsum(i)/real(num_procs, wp))/wall, '%', gmax(i)/max(gsum(i)/real(num_procs, wp), tiny(1._wp)), &
                & gcall(i)/int(num_procs, 8), 1000._wp*(gsum(i)/real(num_procs, wp))/max(real(gcall(i)/int(num_procs, 8), wp), &
                & 1._wp), trim(tl)
        end do
        ! Only rows that were top-level on every rank sum against wall; summing nested rows would drive the residual negative.
        t1sum = sum(gsum, mask=(gthi == 1))/real(num_procs, wp)
        print '(A,F10.3,F19.1,A)', '[phase] RESIDUAL', wall - t1sum, 100._wp*(wall - t1sum)/wall, '%'
        print '(A,F10.3,A,F6.1,A)', '[phase-tier] T1 rows sum to ', t1sum, ' s =', 100._wp*t1sum/wall, &
            & ' % of wall; RESIDUAL is wall minus T1 only'
        if (any(gsum > 0._wp .and. gtlo /= gthi)) then
            write (*, '(A)', advance='no') '[phase-tier] opened at more than one depth, excluded from T1, not summable:'
            do i = 1, PH_N
                if (gsum(i) > 0._wp .and. gtlo(i) /= gthi(i)) write (*, '(1X,A)', advance='no') trim(PH_NAME(i))
            end do
            write (*, '(A)') ''
        end if
        if (any(gbad > 0)) then
            print '(A,I0,A,I0,A)', '[phase-tier] BUDGET INVALID: ', gbad(1), ' interleaved and ', gbad(2), &
                & ' orphan brackets, so some row double counts'
        else
            print '(A)', '[phase-tier] budget valid: every bracket closed innermost-first'
        end if

    end subroutine s_phase_report

end module m_phase_timing
