!>
!!@file
!!@brief Contains module m_amr_xchg_audit

#:include 'macros.fpp'

!> @brief I1a of the plan-based exchange program (docs/documentation/amr_plan_based_exchange.md): per-call-site accounting of every
!! AMR point-to-point MPI transfer. Records what is ACTUALLY sent at the MPI call itself - never re-derived from the metadata the
!! callers read - so the I2+ plan conversions have a ground-truth baseline (message counts, words, tag ranges) and a per-family
!! conservation check (global sends == global recvs) that runs at finalize. Recording is a few integer adds per message; the report
!! prints under rank_time_wrt like the other [amr-*] instruments. Per-xfer identity headers and the destination-tiling assert are
!! I1b and layer on this registry.
module m_amr_xchg_audit

    use m_precision_select
    use m_global_parameters, only: rank_time_wrt, proc_rank, num_procs
    use m_mpi_proxy, only: s_mpi_abort

#ifdef MFC_MPI
    use mpi
#endif

    implicit none

    ! Exchange families (amr_plan_based_exchange.md, "The exchange-family inventory").
    ! XA_FL0 covers the ten s_l0_* tile-routing sites outside the seven families (instrumented
    ! read-only pending the D-l0 decision); XA_F4 is migration.
    integer, parameter :: XA_F1 = 1, XA_F2 = 2, XA_F3 = 3, XA_F4 = 4, XA_F5 = 5, XA_F6 = 6, XA_F7 = 7, XA_FL0 = 8
    integer, parameter :: XA_NFAM = 8

    ! Call-site registry. One id per PHYSICAL MPI call site (fypp twins get their own ids where
    ! they carry different payloads). Names are assigned in init; the id constants are the
    ! documentation at the call sites.
    integer, parameter :: XA_F1_SND = 1         !< s_amr_gather_coarse_patch pooled ISEND
    integer, parameter :: XA_F1_RCV = 2         !< s_amr_gather_coarse_patch IRECV
    integer, parameter :: XA_F3_SND = 3         !< s_amr_gather_coarse_patch_pbmv blocking SEND
    integer, parameter :: XA_F3_RCV = 4         !< s_amr_gather_coarse_patch_pbmv IRECV
    integer, parameter :: XA_F2_SND = 5         !< s_amr_gather_from_parent pooled ISEND (cons+stor instantiations)
    integer, parameter :: XA_F2_RCV = 6         !< s_amr_gather_from_parent blocking RECV
    integer, parameter :: XA_F4_SND = 7         !< s_amr_regrid_stash_migrate ISEND
    integer, parameter :: XA_F4_RCV = 8         !< s_amr_regrid_stash_migrate IRECV
    integer, parameter :: XA_F5_FACE_SND = 9    !< reflux face ISEND (freg lo/hi, tags 2*D/2*D+1)
    integer, parameter :: XA_F5_FACE_RCV = 10   !< reflux face IRECV
    integer, parameter :: XA_F5_FREG_SND = 11   !< freg-to-parent SEND (tags 40+)
    integer, parameter :: XA_F5_FREG_RCV = 12   !< freg-to-parent RECV
    integer, parameter :: XA_F6_XY = 13         !< seam halo SENDRECV, x-side of the pair (tag 4200 out)
    integer, parameter :: XA_F6_YX = 14         !< seam halo SENDRECV, y-side of the pair (tag 4201 out)
    integer, parameter :: XA_F7A_SND = 15       !< s_restrict_fine_to_coarse ISEND
    integer, parameter :: XA_F7A_RCV = 16       !< s_restrict_fine_to_coarse RECV
    integer, parameter :: XA_F7B_SND = 17       !< s_amr_restrict_to_parent SEND
    integer, parameter :: XA_F7B_RCV = 18       !< s_amr_restrict_to_parent RECV
    integer, parameter :: XA_F7C_SND = 19       !< s_amr_scatter_pbmv ISEND
    integer, parameter :: XA_F7C_RCV = 20       !< s_amr_scatter_pbmv RECV
    integer, parameter :: XA_L0_FILL_SND = 21   !< s_l0_fill_tiles_from_coarse SEND
    integer, parameter :: XA_L0_FILL_RCV = 22   !< s_l0_fill_tiles_from_coarse RECV
    integer, parameter :: XA_L0_SCAT_SND = 23   !< s_l0_scatter_tiles_to_coarse SEND
    integer, parameter :: XA_L0_SCAT_RCV = 24   !< s_l0_scatter_tiles_to_coarse RECV
    integer, parameter :: XA_L0_RFLX_SND = 25   !< s_l0_add_reflux_to_tiles SEND
    integer, parameter :: XA_L0_RFLX_RCV = 26   !< s_l0_add_reflux_to_tiles RECV
    integer, parameter :: XA_L0_REST_SND = 27   !< s_l0_restrict_to_tiles SEND (tag 4400+k)
    integer, parameter :: XA_L0_REST_RCV = 28   !< s_l0_restrict_to_tiles RECV
    integer, parameter :: XA_L0_MIGR_SND = 29   !< s_l0_migrate_tile SEND (tag 4300)
    integer, parameter :: XA_L0_MIGR_RCV = 30   !< s_l0_migrate_tile RECV
    integer, parameter :: XA_F1W_SND = 31       !< s_amr_stage_fill_wave per-peer aggregated q ISEND (I2a)
    integer, parameter :: XA_F1W_RCV = 32       !< s_amr_stage_fill_wave per-peer aggregated q IRECV
    integer, parameter :: XA_F3W_SND = 33       !< s_amr_stage_fill_wave per-peer aggregated pb/mv ISEND
    integer, parameter :: XA_F3W_RCV = 34       !< s_amr_stage_fill_wave per-peer aggregated pb/mv IRECV
    integer, parameter :: XA_F2W_SND = 35       !< s_amr_parent_fill_wave per-peer aggregated ISEND (I3)
    integer, parameter :: XA_F2W_RCV = 36       !< s_amr_parent_fill_wave per-peer aggregated IRECV
    integer, parameter :: XA_F6W_SND = 37       !< s_amr_fine_fine_halo per-peer aggregated ISEND (I5-F6)
    integer, parameter :: XA_F6W_RCV = 38       !< s_amr_fine_fine_halo per-peer aggregated IRECV
    integer, parameter :: XA_F5W_FACE_SND = 39  !< s_amr_reflux_faces_wave ISEND (I5-F5a, zero-copy)
    integer, parameter :: XA_F5W_FACE_RCV = 40  !< s_amr_reflux_faces_wave IRECV
    integer, parameter :: XA_F5W_FREG_SND = 41  !< s_amr_freg_wave ISEND (I5-F5b)
    integer, parameter :: XA_F5W_FREG_RCV = 42  !< s_amr_freg_wave IRECV
    integer, parameter :: XA_F7W_SND = 43       !< s_amr_restrict_l1_wave per-peer aggregated ISEND (I5b)
    integer, parameter :: XA_F7W_RCV = 44       !< s_amr_restrict_l1_wave per-peer aggregated IRECV
    integer, parameter :: XA_F7BW_SND = 45      !< s_amr_restrict_parent_wave per-peer aggregated ISEND (I5b)
    integer, parameter :: XA_F7BW_RCV = 46      !< s_amr_restrict_parent_wave per-peer aggregated IRECV
    integer, parameter :: XA_NSITE = 46
    integer, parameter :: xa_fam(XA_NSITE) = [XA_F1, XA_F1, XA_F3, XA_F3, XA_F2, XA_F2, XA_F4, XA_F4, XA_F5, XA_F5, XA_F5, XA_F5, &
                                 & XA_F6, XA_F6, XA_F7, XA_F7, XA_F7, XA_F7, XA_F7, XA_F7, XA_FL0, XA_FL0, XA_FL0, XA_FL0, &
                                 & XA_FL0, XA_FL0, XA_FL0, XA_FL0, XA_FL0, XA_FL0, XA_F1, XA_F1, XA_F3, XA_F3, XA_F2, XA_F2, &
                                 & XA_F6, XA_F6, XA_F5, XA_F5, XA_F5, XA_F5, XA_F7, XA_F7, XA_F7, XA_F7]

    ! dir 1 = send, 2 = recv; a SENDRECV site records both.
    integer(8) :: xa_msgs(XA_NSITE, 2) = 0_8
    integer(8) :: xa_words(XA_NSITE, 2) = 0_8
    !> M0 ORDER ORACLE: per-site XOR folds of mix(pair, seq, key). Sender folds (unordered pair id, its nth-send-on-that-channel,
    !! key); receiver folds the same triple for its nth-receive. Under the FIFO order contract the triples coincide
    !! message-for-message, so the global BXOR of send-folds equals the global BXOR of recv-folds; any cross-rank ordering
    !! divergence misaligns key<->seq and the finalize assert fires. The mixing hash is LOAD-BEARING: a raw packed-field XOR is
    !! provably blind to pairwise transpositions (property control: amr-bench/tools/oracle_ctl.f90, verdict in
    !! logs/oracle_ctl_final.log). M1: seq is PLAN-DERIVED at the call site (each end computes the message's position in the pair's
    !! canonically ordered transfer list from replicated metadata), so no per-peer O(P) counter state exists here anymore. Sites
    !! pass peer/key/seq opt-in; unconverted sites fold nothing.
    integer(8) :: xa_ord(XA_NSITE, 2) = 0_8
    integer    :: xa_seed = -1  !< -1 unread, 0 off, 1 = corrupt one fold (canary), 2 = shift one plan seq (order-swap gate)
    !> canary latch: fire exactly once (an XOR accumulator can return to zero, so testing it would allow self-cancelling double
    !! fires)
    logical :: xa_seeded = .false.
    integer :: xa_tag_min(XA_NSITE) = huge(0)
    integer :: xa_tag_max(XA_NSITE) = -huge(0)

    ! I1b: per-xfer identity header (amr_plan_based_exchange.md "I1b implementation binding").
    ! XA_NH real(wp) words - [site, blk, bl(3), bh(3)] as exact integer-valued reals - are
    ! PREPENDED to each converted family's wire payload under MFC_DEBUG and verified at unpack,
    ! so a plan/pack disagreement (wrong slab, wrong block, crossed families) aborts at the
    ! receiver instead of silently corrupting the patch. Zero in production, so every wire
    ! count/offset adds XA_NH unconditionally and the release arithmetic is untouched.
#ifdef MFC_DEBUG
    integer, parameter :: XA_NH = 8
#else
    integer, parameter :: XA_NH = 0
#endif

    private; public :: s_xa_rec, s_xa_report, s_xa_reset, XA_NH, s_xa_hdr_pack, s_xa_hdr_check, XA_F1_SND, XA_F1_RCV, XA_F3_SND, &
        & XA_F3_RCV, XA_F2_SND, XA_F2_RCV, XA_F4_SND, XA_F4_RCV, XA_F5_FACE_SND, XA_F5_FACE_RCV, XA_F5_FREG_SND, XA_F5_FREG_RCV, &
        & XA_F6_XY, XA_F6_YX, XA_F7A_SND, XA_F7A_RCV, XA_F7B_SND, XA_F7B_RCV, XA_F7C_SND, XA_F7C_RCV, XA_L0_FILL_SND, &
        & XA_L0_FILL_RCV, XA_L0_SCAT_SND, XA_L0_SCAT_RCV, XA_L0_RFLX_SND, XA_L0_RFLX_RCV, XA_L0_REST_SND, XA_L0_REST_RCV, &
        & XA_L0_MIGR_SND, XA_L0_MIGR_RCV, XA_F1W_SND, XA_F1W_RCV, XA_F3W_SND, XA_F3W_RCV, XA_F2W_SND, XA_F2W_RCV, XA_F6W_SND, &
        & XA_F6W_RCV, XA_F5W_FACE_SND, XA_F5W_FACE_RCV, XA_F5W_FREG_SND, XA_F5W_FREG_RCV, XA_F7W_SND, XA_F7W_RCV, XA_F7BW_SND, &
        & XA_F7BW_RCV

contains

    !> Record one transfer at the MPI call site itself. idir: 1 = send, 2 = recv. nwords is the MPI count argument verbatim; tag is
    !! the tag argument verbatim.
    impure subroutine s_xa_rec(isite, idir, nwords, tag, peer, key, seq)

        integer, intent(in)           :: isite, idir, nwords, tag
        integer, intent(in), optional :: peer, key, seq
        integer(8)                    :: pid, e
        integer                       :: sq

        xa_msgs(isite, idir) = xa_msgs(isite, idir) + 1_8
        xa_words(isite, idir) = xa_words(isite, idir) + int(nwords, 8)
        xa_tag_min(isite) = min(xa_tag_min(isite), tag)
        xa_tag_max(isite) = max(xa_tag_max(isite), tag)
        if (present(peer) .and. present(key)) then
            @:ASSERT(present(seq), "keyed xa site must pass its plan-derived seq")
            ! seeded-bug gates: MFC_XA_SEED=1 corrupts exactly one fold on rank 0 (send side, first keyed
            ! message) -- proves the end-to-end wiring can fail; MFC_XA_SEED=2 shifts one plan seq on rank 0's
            ! send side -- models a sender deriving a DIFFERENT plan order, the failure M1's keyed tags exist
            ! to catch. The finalize oracle MUST abort under either seed.
            if (xa_seed < 0) then
                block
                    character(len=8) :: ev
                    integer          :: st
                    call get_environment_variable("MFC_XA_SEED", ev, status=st)
                    xa_seed = 0
                    if (st == 0 .and. ev(1:1) == '1') xa_seed = 1
                    if (st == 0 .and. ev(1:1) == '2') xa_seed = 2
                end block
            end if
            if (xa_seed == 1 .and. proc_rank == 0 .and. idir == 1 .and. .not. xa_seeded) then
                xa_ord(isite, 1) = ieor(xa_ord(isite, 1), 12345_8)  ! the canary corruption
                xa_seeded = .true.
            end if
            sq = seq
            if (xa_seed == 2 .and. proc_rank == 0 .and. idir == 1 .and. .not. xa_seeded) then
                sq = sq + 1  ! the plan-order divergence
                xa_seeded = .true.
            end if
            pid = ior(ishft(int(min(proc_rank, peer), 8), 20), int(max(proc_rank, peer), 8))
            e = f_xa_mix(ieor(f_xa_mix(ieor(pid, ishft(int(sq, 8), 44))), int(key, 8)))
            xa_ord(isite, idir) = ieor(xa_ord(isite, idir), e)
        end if

    end subroutine s_xa_rec

    pure integer(8) function f_xa_mix(x) result(z)
        integer(8), intent(in) :: x
        z = x + int(z'9E3779B97F4A7C15', 8)
        z = ieor(z, ishft(z, -30))*int(z'BF58476D1CE4E5B9', 8)
        z = ieor(z, ishft(z, -27))*int(z'94D049BB133111EB', 8)
        z = ieor(z, ishft(z, -31))

    end function f_xa_mix

    !> Write the XA_NH-word identity header into buf(1:XA_NH): the sending site id, the block the payload is for, and the slab [bl,
    !! bh] the sender packed. Call only under `if (XA_NH > 0)`.
    impure subroutine s_xa_hdr_pack(buf, isite, blk, bl, bh)

        real(wp), intent(inout) :: buf(:)
        integer, intent(in)     :: isite, blk, bl(3), bh(3)

        buf(1) = real(isite, wp); buf(2) = real(blk, wp)
        buf(3:5) = real(bl, wp); buf(6:8) = real(bh, wp)

    end subroutine s_xa_hdr_pack

    !> Verify a received header against what THIS unpack believes it is consuming. isite is the expected SENDING site id (the
    !! matched _SND constant). Aborts with both sides on mismatch - a plan/pack disagreement caught at the wire, before it corrupts
    !! the patch.
    impure subroutine s_xa_hdr_check(buf, isite, blk, bl, bh)

        real(wp), intent(in) :: buf(:)
        integer, intent(in)  :: isite, blk, bl(3), bh(3)
        integer              :: got(8), i
        character(len=256)   :: msg

        got = nint(buf(1:8))
        if (got(1) /= isite .or. got(2) /= blk .or. any(got(3:5) /= bl) .or. any(got(6:8) /= bh)) then
            write (msg, &
                   & '(A,I0,A,I0,A,3(I0,1x),A,3(I0,1x),A,I0,A,I0,A,3(I0,1x),A,3(I0,1x))') &
                   & 'amr xchg header mismatch: expected site ', isite, ' blk ', blk, ' lo ', (bl(i), i=1, 3), 'hi ', (bh(i), &
                   & i=1, 3), '| got site ', got(1), ' blk ', got(2), ' lo ', (got(i), i=3, 5), 'hi ', (got(i), i=6, 8)
            call s_mpi_abort(trim(msg))
        end if

    end subroutine s_xa_hdr_check

    !> Zero the accumulators (a future per-window use; finalize-report runs cumulative).
    impure subroutine s_xa_reset()

        xa_msgs = 0_8; xa_words = 0_8; xa_ord = 0_8
        xa_tag_min = huge(0); xa_tag_max = -huge(0)

    end subroutine s_xa_reset

    !> Finalize-time report + the per-family conservation check: global send msgs/words must equal global recv msgs/words within
    !! each family (every family's sites are internally matched; a violation means a dropped, duplicated, or misattributed
    !! transfer). Collective over MPI_COMM_WORLD - call it from a point every rank reaches (module finalize).
    impure subroutine s_xa_report()

        integer(8)                  :: fam_m(XA_NFAM, 2), fam_w(XA_NFAM, 2)
        integer(8)                  :: gm(XA_NFAM, 2), gw(XA_NFAM, 2)
        integer                     :: i, f
        character(len=3), parameter :: fam_name(XA_NFAM) = ['F1 ', 'F2 ', 'F3 ', 'F4 ', 'F5 ', 'F6 ', 'F7 ', 'L0 ']

        fam_m = 0_8; fam_w = 0_8
        do i = 1, XA_NSITE
            f = xa_fam(i)
            fam_m(f,:) = fam_m(f,:) + xa_msgs(i,:)
            fam_w(f,:) = fam_w(f,:) + xa_words(i,:)
        end do
        gm = fam_m; gw = fam_w
#ifdef MFC_MPI
        block
            integer :: ierr
            call MPI_ALLREDUCE(fam_m, gm, XA_NFAM*2, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
            call MPI_ALLREDUCE(fam_w, gw, XA_NFAM*2, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
        end block
#endif
        if (rank_time_wrt .and. proc_rank == 0) then
            do f = 1, XA_NFAM
                if (gm(f, 1) == 0_8 .and. gm(f, 2) == 0_8) cycle
                write (0, '(A,A,A,I0,A,I0,A,I0,A,I0,A)') '[amr-xa] ', fam_name(f), ' snd ', gm(f, 1), ' msgs ', gw(f, 1), &
                       & ' words | rcv ', gm(f, 2), ' msgs ', gw(f, 2), ' words'
            end do
        end if
        ! conservation: what the senders posted is what the receivers took, family by family
        do f = 1, XA_NFAM
            @:ASSERT(gm(f, 1) == gm(f, 2), "amr xchg audit: send/recv MESSAGE count mismatch in family "//fam_name(f))
            @:ASSERT(gw(f, 1) == gw(f, 2), "amr xchg audit: send/recv WORD count mismatch in family "//fam_name(f))
        end do

        ! M0 order-oracle finalize check (see the xa_ord docs above)
        block
            integer(8) :: og(XA_NSITE, 2)
            integer    :: isite, mierr2
            og = xa_ord
#ifdef MFC_MPI
            call MPI_ALLREDUCE(MPI_IN_PLACE, og, XA_NSITE*2, MPI_INTEGER8, MPI_BXOR, MPI_COMM_WORLD, mierr2)
#endif
            ! sends and receives of the same traffic live under PAIRED sites (XA_*_SND vs XA_*_RCV):
            ! aggregate by FAMILY (xa_fam) -- a per-site snd-vs-rcv compare mismatches structurally on
            ! every healthy run (found by exactly that false positive).
            block
                integer(8) :: fs(0:63), fr(0:63)
                integer    :: fam
                fs = 0_8; fr = 0_8
                do isite = 1, XA_NSITE
                    fam = xa_fam(isite)
                    fs(fam) = ieor(fs(fam), og(isite, 1))
                    fr(fam) = ieor(fr(fam), og(isite, 2))
                end do
                do fam = 0, 63
                    if (fs(fam) /= fr(fam)) then
                        print '(A,I0,A,Z16,A,Z16)', ' [amr-xa] ORDER ORACLE MISMATCH family ', fam, ' snd ', fs(fam), ' rcv ', &
                            & fr(fam)
                        call s_mpi_abort('AMR exchange order oracle: send/recv order hashes diverge.')
                    end if
                end do
            end block
        end block

    end subroutine s_xa_report

end module m_amr_xchg_audit
