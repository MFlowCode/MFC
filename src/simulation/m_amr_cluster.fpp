!>
!!@file
!!@brief Contains module m_amr_cluster

#:include 'macros.fpp'

!> @brief Berger-Rigoutsos clustering of a sparse tag list into separated block boxes, identically on every rank: per-box tag
!! signatures, the bisection, the canonical Morton order and the min-separation merge. The regrid (m_amr_regrid) builds the tag list
!! and shapes the boxes it gets back.
module m_amr_cluster

#ifdef MFC_MPI
    use mpi  !< per-node signature reductions for the rank-invariant tree walk
#endif

    use m_derived_types  ! t_box
    use m_box, only: f_morton
    use m_global_parameters
    use m_mpi_proxy, only: s_mpi_abort
    use m_amr_state
    use m_amr_distribution

    implicit none

    private
    public :: s_amr_cluster

contains

    !> Concatenated 1D tag signatures of box [blo0:bhi0], built from the tag range [ts:te] in one pass. Axis d occupies sig(off(d) :
    !! off(d) + ext(d) - 1), and sig(off(d) + t - blo0(d)) counts the in-box tagged cells at position t along d. One signature
    !! serves the trim, the in-box count and every candidate split, so the tag list is scanned once per box. nsig returns the used
    !! length.
    impure subroutine s_amr_box_sig(tags, ts, te, blo0, bhi0, sig, off, nsig)

        integer, intent(in)  :: tags(:,:), ts, te, blo0(3), bhi0(3)
        integer, intent(out) :: sig(:), off(3), nsig
        integer              :: d, t, c(3)

        nsig = 0
        do d = 1, 3
            off(d) = nsig + 1
            if (d <= num_dims) nsig = nsig + (bhi0(d) - blo0(d) + 1)
        end do
        sig(1:nsig) = 0
        do t = ts, te
            c = tags(:,t)
            if (c(1) < blo0(1) .or. c(1) > bhi0(1)) cycle
            if (c(2) < blo0(2) .or. c(2) > bhi0(2)) cycle
            if (c(3) < blo0(3) .or. c(3) > bhi0(3)) cycle
            do d = 1, num_dims
                sig(off(d) + c(d) - blo0(d)) = sig(off(d) + c(d) - blo0(d)) + 1
            end do
        end do

    end subroutine s_amr_box_sig

    !> Shrink box [blo:bhi] to the tight bbox of its tagged cells and return their count, both read off the signature of
    !! [blo0:bhi0]. Equivalent to scanning the tag list: the per-axis min/max of the contained tags are the first and last nonzero
    !! of that axis signature, and "any tagged" is "the signature sums nonzero". ok=.false. if none tagged. Collapsed dims (lo=hi=0)
    !! survive unchanged, their signature being a single bin.
    pure subroutine s_amr_trim_from_sig(sig, off, blo0, bhi0, blo, bhi, ok, ntag)

        integer, intent(in)    :: sig(:), off(3), blo0(3), bhi0(3)
        integer, intent(inout) :: blo(3), bhi(3)
        logical, intent(out)   :: ok
        integer, intent(out)   :: ntag
        integer                :: d, t, lo, hi

        ok = .false.
        ntag = 0
        do t = blo0(1), bhi0(1)
            ntag = ntag + sig(off(1) + t - blo0(1))
        end do
        if (ntag == 0) return  ! no tags in the box; every axis signature is empty too
        do d = 1, num_dims
            lo = -1; hi = -1
            do t = blo0(d), bhi0(d)
                if (sig(off(d) + t - blo0(d)) > 0) then
                    if (lo < 0) lo = t
                    hi = t
                end if
            end do
            blo(d) = lo; bhi(d) = hi
        end do
        ok = .true.

    end subroutine s_amr_trim_from_sig

    !> Berger-Rigoutsos bisection of one (already tagged-trimmed) candidate box, read off the signature of [blo0:bhi0]: pick the
    !! longest splittable axis, prefer a zero-signature hole (widest interior run), else the strongest signature inflection
    !! (Laplacian sign change). ok=.false. if no axis admits a split leaving both children >= 2 cells. Slicing the signature to the
    !! trimmed range is exact: trim shrinks only to the tags' own bbox, so no tag leaves the box. Integer-only => identical on all
    !! ranks.
    pure subroutine s_amr_find_split_sig(sig, off, blo0, blo, bhi, sax, spos, ok)

        integer, intent(in)  :: sig(:), off(3), blo0(3)
        integer, intent(in)  :: blo(3), bhi(3)
        integer, intent(out) :: sax, spos
        logical, intent(out) :: ok
        !> Minimum child extent along the split axis, i.e. the smallest box the bisection may produce. 2 is the algorithmic floor;
        !! amr_blocking_factor raises it, which is what stops the bisection over-generating (with the floor at 2 the recursion
        !! splits until the amr_max_blocks cap stops it and the min-separation merge then collapses the result back). Note this is a
        !! minimum size, not AMReX's blocking factor: AMReX coarsens the tag lattice, but coarsening a rank-local sparse tag list
        !! cannot dedup coarse cells that straddle a rank boundary without an extra exchange, so the size floor is used instead.
        integer :: min_child
        integer :: axord(3), ext(3), d, ax, t, s, b
        integer :: run, run_start, best_run, best_start, lap, prevlap, bestmag, bestpos

        min_child = max(2, amr_blocking_factor)
        ok = .false.; sax = 0; spos = 0
        ext = bhi - blo + 1
        axord = [1, 2, 3]  ! sort axes by descending extent (deterministic bubble)
        do d = 1, 2
            do ax = 1, 3 - d
                if (ext(axord(ax)) < ext(axord(ax + 1))) then
                    s = axord(ax); axord(ax) = axord(ax + 1); axord(ax + 1) = s
                end if
            end do
        end do
        do d = 1, 3
            ax = axord(d)
            if (ax > num_dims) cycle
            if (ext(ax) < 2*min_child) cycle
            b = off(ax) - blo0(ax)  ! signature of position t on this axis is sig(b + t)
            ! (1) widest interior zero run (box is trimmed => sig(blo)>0 and sig(bhi)>0, so any run is interior)
            best_run = 0; best_start = -1; run = 0; run_start = -1
            do t = blo(ax), bhi(ax)
                if (sig(b + t) == 0) then
                    if (run == 0) run_start = t
                    run = run + 1
                else
                    if (run > best_run) then; best_run = run; best_start = run_start; end if
                    run = 0
                end if
            end do
            if (best_start > blo(ax)) then
                spos = best_start
                if (spos - blo(ax) >= min_child .and. bhi(ax) - spos + 1 >= min_child) then
                    sax = ax; ok = .true.; return
                end if
            end if
            ! (2) strongest inflection: Laplacian sign change with the largest jump
            bestmag = -1; bestpos = -1; prevlap = 0
            do t = blo(ax) + 1, bhi(ax) - 1
                lap = sig(b + t - 1) - 2*sig(b + t) + sig(b + t + 1)
                if (t > blo(ax) + 1) then
                    if (((lap < 0) .neqv. (prevlap < 0)) .and. abs(lap - prevlap) > bestmag .and. t - blo(ax) >= min_child &
                        & .and. bhi(ax) - t + 1 >= min_child) then
                        bestmag = abs(lap - prevlap); bestpos = t
                    end if
                end if
                prevlap = lap
            end do
            if (bestpos > 0) then
                sax = ax; spos = bestpos; ok = .true.; return
            end if
        end do

    end subroutine s_amr_find_split_sig

    !> Cluster a sparse tag list (level-0 cell coords, tags(1:3, 1:ntag_in)) into a list of separated block boxes, identically on
    !! every rank. Caller builds the list (s_amr_local_tags / s_amr_pack_gwin_pairs); per-rank memory is O(#tagged), not O(global
    !! grid). Berger-Rigoutsos recursive bisection until each box's tag efficiency reaches amr_cluster_eff (or it is atomic / the
    !! amr_max_blocks cap is hit), then merges any two boxes whose amr_buf-padded extents come within buff_size (so no fine-fine
    !! adjacency: separated boxes stay >= buff_size apart, nearby ones collapse to one box, their bounding box). Boxes are raw
    !! tagged extents; the caller pads, clamps, size-caps each.
    impure subroutine s_amr_cluster(tags, ntag_in, boxes, nboxes, reduce)

        integer, intent(in) :: tags(:,:), ntag_in
        !> .true.: `tags` is this rank's local list and each node's signature is reduced across ranks, so the tree is driven by
        !! global counts without any rank holding the global tag list. .false.: `tags` is already replicated on every rank.
        logical, intent(in)                   :: reduce
        type(t_box), allocatable, intent(out) :: boxes(:)
        integer, intent(out)                  :: nboxes
        integer, allocatable                  :: slo(:,:), shi(:,:), alo(:,:), ahi(:,:)
        integer, allocatable                  :: sts(:), ste(:), wt(:,:)
        integer, allocatable                  :: sdep(:)  !< recursion depth carried with each stack entry
        integer, allocatable                  :: sperm(:), t2lo(:,:), t2hi(:,:)
        integer                               :: dep
        integer, allocatable                  :: sig(:)   !< concatenated per-axis tag signature of the node's box
        integer, allocatable                  :: ovr(:)   !< scratch: ranks overlapping the node's box
        integer                               :: novr
        integer                               :: blo0(3), bhi0(3), off(3), nsig

#ifdef MFC_MPI
        integer :: ierr
#endif
        integer                 :: mg, ng, pg, t
        integer                 :: cap, nacc, i, j, k, d, sax, spos, thr, ntag
        integer(8), allocatable :: akey(:)  !< Morton key of each accepted box's lo, the canonical merge order
        integer, allocatable    :: nxt(:)  !< singly-linked survivor list: removal is O(1), so the merge is O(n) not O(n^2)
        integer                 :: head, nlive
        integer, allocatable    :: bp(:), bidx(:)  !< bin back-links + current bin of each live box (incremental refile)
        integer                 :: dirty, aa, jb2, extd
        logical                 :: need_build
        integer                 :: blo3(3), bhi3(3), nbmax, rng
        integer, allocatable    :: prv(:)  !< predecessor links: a binned hit unlinks in O(1)
        integer, allocatable    :: bh(:), bc(:)  !< bin heads + per-box chains (host scratch, rebuilt per pass)
        integer                 :: ext_max, cellw, nbx, nby, nbz, nb_tot, jbest
        integer, allocatable    :: gcnt(:), gdsp(:), sbx(:,:), gbx(:,:)  !< union of the per-rank accepted boxes
        integer                 :: ntot
        !> The level-order walk. kpos/kbat index the nodes kept at the current depth; bsig concatenates the signatures of that
        !! depth's shared nodes into the single buffer the one reduction covers, with bofs/blen/boff their slices.
        integer, allocatable :: kpos(:), kbat(:), bofs(:), blen(:), boff(:,:), bsig(:)
        integer              :: ncur, nnxt, nkeep, nbat, nbuf
        !> A node is wide when its box spans more than this many ranks. Wide nodes use the batched collective (every rank overlaps
        !! them and needs the answer); narrow ones reduce among their few overlapping ranks. The threshold only has to keep the wide
        !! count at O(log P); 8 is one 2x2x2 brick of ranks, the shape a seam node has.
        integer, parameter   :: amr_cl_wide = 8
        integer, allocatable :: bnov(:), bovr(:,:), wbuf(:)
        logical, allocatable :: bwide(:)
        integer, allocatable :: pidx(:), plist(:), scnt(:), rcnt(:), sdsp2(:), rdsp2(:), soff(:), roff(:)
        integer, allocatable :: sbuf(:), rbuf(:), creq(:)
        integer              :: np2, q, rr, nsnd, nrcv, nreq2, tagc, nwb, o1  ! t is already a loop variable above
        integer(8)           :: vol  !< box volume; a global-bbox first pass can exceed 2**31 cells
        integer              :: blo(3), bhi(3), ts, te, lo, hi, tmp(3)
        logical              :: ok, force, capped, mine
        real(wp)             :: eff

        nboxes = 0
        ! In reduce mode a rank with no local tags must still walk the tree and enter every ALLREDUCE, contributing zeros;
        ! returning early here would deadlock the ranks that do have tags. An all-empty list ends the loop via the trim.
        if (.not. reduce .and. ntag_in == 0) return
        mg = m_glb; ng = 0; pg = 0
        if (n_glb > 0) ng = n_glb
        if (p_glb > 0) pg = p_glb

        cap = amr_max_fine
        allocate (slo(3, 4*cap + 8), shi(3, 4*cap + 8), alo(3, cap), ahi(3, cap))
        allocate (sts(4*cap + 8), ste(4*cap + 8), wt(3, ntag_in), sdep(4*cap + 8))
        allocate (sig(mg + ng + pg + 3))  ! bound: the three full domain extents; reused by every node
        allocate (ovr(amr_cl_wide))  ! only narrow nodes are ever enumerated, so this does not size with P
        allocate (akey(cap))  ! merge-order scratch
        ! working copy of the tag list, partitioned in place as the tree descends so each node scans only its tags
        do t = 1, ntag_in
            wt(:,t) = tags(:,t)
        end do
        ncur = 1; slo(:,1) = [0, 0, 0]; shi(:,1) = [mg, ng, pg]  ! first node trims to the global tagged bbox
        sts(1) = 1; ste(1) = ntag_in
        sdep(1) = 0
        nacc = 0; capped = .false.
        allocate (kpos(4*cap + 8), kbat(4*cap + 8), bofs(4*cap + 8), blen(4*cap + 8), boff(3, 4*cap + 8))
        allocate (bsig(4*(mg + ng + pg + 3)))
        allocate (bnov(4*cap + 8), bwide(4*cap + 8), bovr(amr_cl_wide, 4*cap + 8))
        allocate (pidx(0:max(num_procs - 1, 0)), plist(max(num_procs, 1)))
        allocate (scnt(max(num_procs, 1)), rcnt(max(num_procs, 1)), sdsp2(max(num_procs, 1)), rdsp2(max(num_procs, 1)))
        allocate (soff(max(num_procs, 1)), roff(max(num_procs, 1)))
        allocate (wbuf(1), sbuf(1), rbuf(1), creq(1))
        pidx = 0
        ! Level-order descent: every shared node of a depth rides one reduction, so the collective count is O(tree depth). A
        ! child's box lies inside its parent's, so the ranks overlapping a child are a subset of the parent's: every ancestor of a
        ! shared node is shared, every rank walks the whole shared subtree in the same order, and the per-depth batch is identical
        ! on every rank. Rank-local nodes are excluded from the batch. Nodes 1:ncur are the current depth; children are appended
        ! past ncur and shifted down when the depth closes (peak occupancy 3*ncur <= the 4*cap + 8 the arrays carry).
        do while (ncur > 0)
            nkeep = 0; nbat = 0; nbuf = 0; nnxt = 0
            ! pass 1: classify, and stash the signatures that need reducing (rank-local nodes recompute theirs in pass 2 rather
            ! than swamp the buffer)
            do i = 1, ncur
                blo0 = slo(:,i); bhi0 = shi(:,i)
                novr = f_amr_overlap_count(blo0, bhi0)  ! rank count without enumerating the set (O(P) on a machine-wide box)
                mine = (num_procs == 1) .or. f_amr_rank_overlaps(blo0, bhi0, proc_rank)
                ! A rank holds tags only inside its subdomain, so it would contribute only zeros to a narrow node it does not
                ! reach: drop the subtree (the closing ALLGATHERV carries back anything accepted inside it). Wide nodes are never
                ! dropped: they are settled by a collective every rank must enter with the identical batch, and a wide node's
                ! ancestors are all wide, so every rank reaches every wide node. A narrow node's members all walked its parent
                ! for the same reason, so the p2p pairing below is complete.
                if (reduce .and. num_procs > 1 .and. .not. mine .and. novr <= amr_cl_wide) cycle
                ! one pass over this node's tags yields the signature; trim, count and split all read it
                call s_amr_box_sig(wt, sts(i), ste(i), blo0, bhi0, sig, off, nsig)
                nkeep = nkeep + 1; kpos(nkeep) = i; kbat(nkeep) = 0
                if (reduce .and. num_procs > 1 .and. novr > 1) then
                    call s_amr_size_int(bsig, nbuf + nsig)
                    nbat = nbat + 1; kbat(nkeep) = nbat
                    bofs(nbat) = nbuf; blen(nbat) = nsig; boff(:,nbat) = off
                    bnov(nbat) = novr; bwide(nbat) = (novr > amr_cl_wide)
                    bovr(1, nbat) = -1  ! defined for wide nodes too: Fortran does not promise .or. short-circuits
                    if (.not. bwide(nbat)) then
                        call s_amr_ranks_overlapping(blo0, bhi0, ovr, novr)  ! bounded by amr_cl_wide, so never O(P)
                        bovr(1:novr,nbat) = ovr(1:novr)
                    end if
                    bsig(nbuf + 1:nbuf + nsig) = sig(1:nsig)
                    nbuf = nbuf + nsig
                end if
            end do
#ifdef MFC_MPI
            ! The depth's reduction. Wide nodes (the O(log P) shallow ones near the root, which every rank overlaps) ride one
            ! batched collective. Narrow nodes (the deep ones straddling a rank seam, where the O(P) growth in the shared set
            ! lives) reduce point-to-point among their small rank brick: each member ships its contribution to ovr(1), which sums
            ! and ships the total back. Both ends agree on message contents without negotiation: a rank's node list at a depth is
            ! a subsequence of the one globally ordered walk, enumerated in ascending j on both sides, so one aggregated message
            ! per peer per phase matches unambiguously.
            nwb = 0
            do j = 1, nbat
                if (bwide(j)) nwb = nwb + blen(j)
            end do
            if (nwb > 0) then
                call s_amr_size_int(wbuf, nwb)
                o1 = 0
                do j = 1, nbat
                    if (.not. bwide(j)) cycle
                    wbuf(o1 + 1:o1 + blen(j)) = bsig(bofs(j) + 1:bofs(j) + blen(j)); o1 = o1 + blen(j)
                end do
                call MPI_ALLREDUCE(MPI_IN_PLACE, wbuf, nwb, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
                o1 = 0
                do j = 1, nbat
                    if (.not. bwide(j)) cycle
                    bsig(bofs(j) + 1:bofs(j) + blen(j)) = wbuf(o1 + 1:o1 + blen(j)); o1 = o1 + blen(j)
                end do
            end if
            ! peers for the narrow nodes: whoever roots a node I hold, plus whoever holds a node I root
            np2 = 0
            do j = 1, nbat
                if (bwide(j)) cycle
                if (bovr(1, j) == proc_rank) then
                    do t = 2, bnov(j)
                        rr = bovr(t, j)
                        if (pidx(rr) == 0) then; np2 = np2 + 1; plist(np2) = rr; pidx(rr) = np2; end if
                    end do
                else
                    rr = bovr(1, j)
                    if (pidx(rr) == 0) then; np2 = np2 + 1; plist(np2) = rr; pidx(rr) = np2; end if
                end if
            end do
            if (np2 > 0) then
                scnt(1:np2) = 0; rcnt(1:np2) = 0
                do j = 1, nbat
                    if (bwide(j)) cycle
                    if (bovr(1, j) == proc_rank) then
                        do t = 2, bnov(j); q = pidx(bovr(t, j)); rcnt(q) = rcnt(q) + blen(j); end do
                    else
                        q = pidx(bovr(1, j)); scnt(q) = scnt(q) + blen(j)
                    end if
                end do
                sdsp2(1) = 0; rdsp2(1) = 0
                do q = 2, np2
                    sdsp2(q) = sdsp2(q - 1) + scnt(q - 1); rdsp2(q) = rdsp2(q - 1) + rcnt(q - 1)
                end do
                nsnd = sdsp2(np2) + scnt(np2); nrcv = rdsp2(np2) + rcnt(np2)
                call s_amr_size_int(sbuf, max(nsnd, 1)); call s_amr_size_int(rbuf, max(nrcv, 1))
                call s_amr_size_int(creq, 2*np2)
                ! phase A: every member ships its own contribution up to the node's root
                soff(1:np2) = sdsp2(1:np2)
                do j = 1, nbat
                    if (bwide(j) .or. bovr(1, j) == proc_rank) cycle
                    q = pidx(bovr(1, j))
                    sbuf(soff(q) + 1:soff(q) + blen(j)) = bsig(bofs(j) + 1:bofs(j) + blen(j)); soff(q) = soff(q) + blen(j)
                end do
                tagc = amr_tag_narrow + int(mod(amr_mesh_epoch, 50_8))
                nreq2 = 0
                do q = 1, np2
                    if (rcnt(q) > 0) then
                        nreq2 = nreq2 + 1
                        call MPI_IRECV(rbuf(rdsp2(q) + 1), rcnt(q), MPI_INTEGER, plist(q), tagc, MPI_COMM_WORLD, creq(nreq2), ierr)
                    end if
                end do
                do q = 1, np2
                    if (scnt(q) > 0) then
                        nreq2 = nreq2 + 1
                        call MPI_ISEND(sbuf(sdsp2(q) + 1), scnt(q), MPI_INTEGER, plist(q), tagc, MPI_COMM_WORLD, creq(nreq2), ierr)
                    end if
                end do
                if (nreq2 > 0) call MPI_WAITALL(nreq2, creq, MPI_STATUSES_IGNORE, ierr)
                ! the root sums its members in. Integer SUM is exact and order-independent, so the total is bit-identical to what
                ! a machine-wide reduction would produce; only who is in the message differs, never the arithmetic.
                roff(1:np2) = rdsp2(1:np2)
                do j = 1, nbat
                    if (bwide(j) .or. bovr(1, j) /= proc_rank) cycle
                    do t = 2, bnov(j)
                        q = pidx(bovr(t, j))
                        bsig(bofs(j) + 1:bofs(j) + blen(j)) = bsig(bofs(j) + 1:bofs(j) + blen(j)) + rbuf(roff(q) + 1:roff(q) &
                             & + blen(j))
                        roff(q) = roff(q) + blen(j)
                    end do
                end do
                ! phase B: the total goes back down. Counts mirror phase A exactly, so the buffers swap roles.
                roff(1:np2) = rdsp2(1:np2)
                do j = 1, nbat
                    if (bwide(j) .or. bovr(1, j) /= proc_rank) cycle
                    do t = 2, bnov(j)
                        q = pidx(bovr(t, j))
                        rbuf(roff(q) + 1:roff(q) + blen(j)) = bsig(bofs(j) + 1:bofs(j) + blen(j)); roff(q) = roff(q) + blen(j)
                    end do
                end do
                nreq2 = 0
                do q = 1, np2
                    if (scnt(q) > 0) then
                        nreq2 = nreq2 + 1
                        call MPI_IRECV(sbuf(sdsp2(q) + 1), scnt(q), MPI_INTEGER, plist(q), tagc + 50, MPI_COMM_WORLD, &
                                       & creq(nreq2), ierr)
                    end if
                end do
                do q = 1, np2
                    if (rcnt(q) > 0) then
                        nreq2 = nreq2 + 1
                        call MPI_ISEND(rbuf(rdsp2(q) + 1), rcnt(q), MPI_INTEGER, plist(q), tagc + 50, MPI_COMM_WORLD, &
                                       & creq(nreq2), ierr)
                    end if
                end do
                if (nreq2 > 0) call MPI_WAITALL(nreq2, creq, MPI_STATUSES_IGNORE, ierr)
                soff(1:np2) = sdsp2(1:np2)
                do j = 1, nbat
                    if (bwide(j) .or. bovr(1, j) == proc_rank) cycle
                    q = pidx(bovr(1, j))
                    bsig(bofs(j) + 1:bofs(j) + blen(j)) = sbuf(soff(q) + 1:soff(q) + blen(j)); soff(q) = soff(q) + blen(j)
                end do
                do q = 1, np2  ! clear only what was touched: a full wipe would be O(P) per depth
                    pidx(plist(q)) = 0
                end do
            end if
#endif
            ! pass 2: trim, accept or split every node kept at this depth
            do j = 1, nkeep
                i = kpos(j); blo = slo(:,i); bhi = shi(:,i); ts = sts(i); te = ste(i); dep = sdep(i)
                blo0 = blo; bhi0 = bhi
                if (kbat(j) > 0) then
                    off = boff(:,kbat(j))
                    sig(1:blen(kbat(j))) = bsig(bofs(kbat(j)) + 1:bofs(kbat(j)) + blen(kbat(j)))
                    nsig = blen(kbat(j))
                else
                    ! rank-local: no reduction was needed, so the signature is recomputed here rather than carried. Safe because
                    ! pass 2 only ever partitions a node's own wt(:, ts:te) range, which is disjoint from every other node's.
                    call s_amr_box_sig(wt, ts, te, blo0, bhi0, sig, off, nsig)
                end if
                call s_amr_trim_from_sig(sig(1:nsig), off, blo0, bhi0, blo, bhi, ok, ntag)
                if (.not. ok) cycle
                vol = 1_8
                do d = 1, num_dims; vol = vol*int(bhi(d) - blo(d) + 1, 8); end do
                eff = real(ntag, wp)/real(max(vol, 1_8), wp)
                call s_amr_find_split_sig(sig(1:nsig), off, blo0, blo, bhi, sax, spos, ok)
                ! splitting now could overflow the amr_max_blocks cap. Under the level-order walk the pending set is the rest
                ! of this depth plus the children queued so far; with the blocking-factor floor the bisection normally stays
                ! clear of the cap, so this guard is inert.
                force = (nacc + (nkeep - j) + nnxt + 1 >= cap)
                if (eff >= amr_cluster_eff .or. .not. ok .or. force) then
                    if (nacc < cap) then; nacc = nacc + 1; alo(:,nacc) = blo; ahi(:,nacc) = bhi; end if
                    if (force .and. ok .and. eff < amr_cluster_eff) capped = .true.
                else
                    ! partition wt(:, ts:te) in place: coord(sax) < spos to the front (low child), >= spos to the back (high)
                    lo = ts; hi = te
                    do while (lo <= hi)
                        if (wt(sax, lo) < spos) then
                            lo = lo + 1
                        else
                            tmp = wt(:,lo); wt(:,lo) = wt(:,hi); wt(:,hi) = tmp
                            hi = hi - 1
                        end if
                    end do
                    ! low child = [ts:lo-1], high child = [lo:te]; every parent tag lands in exactly one (box just trimmed+split)
                    slo(:,ncur + nnxt + 1) = blo; shi(:,ncur + nnxt + 1) = bhi; shi(sax, ncur + nnxt + 1) = spos - 1
                    sts(ncur + nnxt + 1) = ts; ste(ncur + nnxt + 1) = lo - 1; sdep(ncur + nnxt + 1) = dep + 1
                    slo(:,ncur + nnxt + 2) = blo; shi(:,ncur + nnxt + 2) = bhi; slo(sax, ncur + nnxt + 2) = spos
                    sts(ncur + nnxt + 2) = lo; ste(ncur + nnxt + 2) = te; sdep(ncur + nnxt + 2) = dep + 1
                    nnxt = nnxt + 2
                end if
            end do
            ! close the depth: the children become the next current level
            do i = 1, nnxt
                slo(:,i) = slo(:,ncur + i); shi(:,i) = shi(:,ncur + i)
                sts(i) = sts(ncur + i); ste(i) = ste(ncur + i); sdep(i) = sdep(ncur + i)
            end do
            ncur = nnxt
        end do
        deallocate (kpos, kbat, bofs, blen, boff, bsig, bnov, bwide, bovr, pidx, plist)
        deallocate (scnt, rcnt, sdsp2, rdsp2, soff, roff, wbuf, sbuf, rbuf, creq)

#ifdef MFC_MPI
        ! Rank-local subtrees are walked only by their owner, so each rank holds just the boxes from the subtrees it owns. Union
        ! them once here: per-box global data, 6 ints per box. Ranks contribute in rank order, which is not the order a serial
        ! traversal would accept them in; the canonical Morton sort immediately below makes the merged result independent of that.
        if (reduce .and. num_procs > 1) then
            allocate (gcnt(num_procs), gdsp(num_procs))
            call MPI_ALLGATHER(nacc, 1, MPI_INTEGER, gcnt, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
            gdsp(1) = 0
            do i = 2, num_procs
                gdsp(i) = gdsp(i - 1) + gcnt(i - 1)
            end do
            ntot = gdsp(num_procs) + gcnt(num_procs)
            allocate (sbx(6, max(nacc, 1)), gbx(6, max(ntot, 1)))
            do i = 1, nacc
                sbx(1:3,i) = alo(:,i); sbx(4:6,i) = ahi(:,i)
            end do
            gcnt = gcnt*6; gdsp = gdsp*6
            call MPI_ALLGATHERV(sbx, nacc*6, MPI_INTEGER, gbx, gcnt, gdsp, MPI_INTEGER, MPI_COMM_WORLD, ierr)
            ! The gathered list is every rank's pre-merge leaves (the bisection splits until its per-rank guard stops it and
            ! relies on the merge below to fuse them back), so it can exceed amr_max_blocks while the merged set does not.
            ! Truncating it to the cap here would drop whole ranks' leaves (the list is in rank order) and silently leave tagged
            ! cells unrefined. The accepted arrays grow to the union instead; the cap is applied to the merged set, below.
            if (ntot > size(alo, 2)) then
                deallocate (alo, ahi, akey)
                allocate (alo(3, ntot), ahi(3, ntot), akey(ntot))
            end if
            nacc = ntot
            do i = 1, nacc
                alo(:,i) = gbx(1:3,i); ahi(:,i) = gbx(4:6,i)
            end do
            deallocate (gcnt, gdsp, sbx, gbx)
        end if
#endif

        ! Canonicalise the merge input: the merge scans in list order and fuses the first too-close pair, so sorting by Morton
        ! of lo makes its output a function of the box set alone, not of the (rank-dependent) acceptance order. Accepted boxes
        ! are disjoint, so their lo corners are distinct and the key is a total order under f_morton's 21 bits/dim; the sort is
        ! stable, so a collision above that bound falls back to acceptance order on every rank alike. Morton keeps spatial
        ! neighbours adjacent, so near pairs fuse first and the fused boxes stay compact.
        do i = 1, nacc
            akey(i) = f_morton(alo(1, i), alo(2, i), alo(3, i))
        end do
        ! stable sort on an index permutation, payload applied once
        allocate (sperm(nacc), t2lo(3, nacc), t2hi(3, nacc))
        call s_amr_sort_by_key(akey, nacc, sperm)
        do i = 1, nacc
            t2lo(:,i) = alo(:,sperm(i)); t2hi(:,i) = ahi(:,sperm(i))
        end do
        alo(:,1:nacc) = t2lo; ahi(:,1:nacc) = t2hi
        deallocate (sperm, t2lo, t2hi)

        ! min-separation merge: two boxes are separated only if some active dim's gap reaches thr; else fuse to their bounding box
        thr = buff_size + 2*amr_buf
        ! The survivors must stay in the canonical Morton order the sort established, and the fusion sequence must be
        ! reproducible, because the resulting box set is what every rank must agree on (and what the goldens hold). The
        ! walk is defined as: scan in list order, fuse the first too-close pair (minimum surviving index j for each i),
        ! restart. A next-pointer list removes an absorbed box in O(1) and visits survivors in that same order, so the
        ! same pairs are tested in the same sequence and the same fusions happen.
        allocate (nxt(max(nacc, 1)))
        do i = 1, nacc - 1
            nxt(i) = i + 1
        end do
        if (nacc >= 1) nxt(nacc) = 0
        ! head must be 0 when there is nothing to merge: nacc = 0 is reachable (a regrid where no cell is
        ! tagged globally leaves nacc at its initialization, and the reduce path has no zero-tag guard), and
        ! head = 1 there would enter the walk below and read nxt(1), which was never written.
        head = merge(1, 0, nacc >= 1); nlive = nacc
        ! Binned candidate merge. Soundness of the prune: tooclose(i,j) needs every per-dim gap < thr, which
        ! bounds |alo(d,i)-alo(d,j)| by ext_max + thr - 1, so with bin width ext_max + thr every tooclose
        ! partner of i lies within the 3^d neighbouring bins of i's lo. For each i in list order the minimum
        ! surviving index j among candidates is taken, exactly the first tooclose j a linear walk meets, so the
        ! fusion sequence is unchanged. ext_max can grow when a fusion grows a box, so bins are rebuilt when it
        ! outgrows the cell width.
        allocate (prv(max(nacc, 1)))
        do i = 1, nacc
            prv(i) = i - 1
        end do
        ! Dirty-box continuation. After fusing (i, j) only box i changed, so instead of restarting the pass:
        ! (a) re-test earlier survivors against the grown box (minimum index first, exactly what a restart
        ! would find), else (b) re-test all later survivors, else (c) the chain is exhausted and the walk
        ! resumes at the survivor's live successor; everything to its left is provably clean. Bins are built
        ! once per cellw epoch and maintained incrementally: the absorbed box is unlinked, the survivor
        ! re-filed when its lo crosses a bin (lo = min of members, so it can never drop below the epoch's
        ! blo3). Extent growth past cellw - thr doubles cellw and rebuilds (amortized log(extent range)).
        allocate (bp(max(nacc, 1)), bidx(max(nacc, 1)))
        cellw = 0
        need_build = .true.
        i = head
        outer: do while (i /= 0)
            if (need_build) then
                call s_mrg_build()
                need_build = .false.
            end if
            jbest = f_mrg_qminj(i)
            if (jbest /= 0) then
                dirty = i
                call s_mrg_fuse(dirty, jbest)
                chain: do
                    extd = 1
                    do d = 1, num_dims
                        extd = max(extd, ahi(d, dirty) - alo(d, dirty) + 1)
                    end do
                    if (extd > cellw - thr) then
                        do while (extd > cellw - thr)
                            cellw = cellw*2
                        end do
                        call s_mrg_build()
                    else if (f_mrg_binof(dirty) /= bidx(dirty)) then
                        call s_mrg_binun(dirty)
                        call s_mrg_binreg(dirty)
                    end if
                    aa = f_mrg_qmina(dirty)
                    if (aa /= 0) then
                        call s_mrg_fuse(aa, dirty)
                        dirty = aa
                        cycle chain
                    end if
                    jb2 = f_mrg_qminj(dirty)
                    if (jb2 /= 0) then
                        call s_mrg_fuse(dirty, jb2)
                        cycle chain
                    end if
                    exit chain
                end do chain
                i = nxt(dirty)
            else
                i = nxt(i)
            end if
        end do outer
        if (allocated(bh)) deallocate (bh)
        if (allocated(bc)) deallocate (bc)
        deallocate (prv, bp, bidx)
        ! compact once, in list order
        k = 0
        i = head
        do while (i /= 0)
            k = k + 1
            if (k /= i) then
                alo(:,k) = alo(:,i); ahi(:,k) = ahi(:,i)
            end if
            i = nxt(i)
        end do
        nacc = nlive
        deallocate (nxt)
        if (capped .and. proc_rank == 0) print '(A,I0)', ' [amr] WARNING: tag clustering capped at amr_max_blocks = ', cap
        ! the merged set is what the block pool must hold: past the cap, boxes simply never refine (a correctness cliff, so
        ! it is named, not silent)
        if (nacc > cap) then
            if (proc_rank == 0) print '(A,I0,A,I0)', ' [amr] WARNING: merged box set truncated: ', nacc, ' boxes, keeping ', cap
            nacc = cap
        end if

        nboxes = nacc
        allocate (boxes(nboxes))
        do i = 1, nboxes
            boxes(i)%lo = alo(:,i); boxes(i)%hi = ahi(:,i)
        end do
        deallocate (slo, shi, alo, ahi, sts, ste, wt, sdep, sig, ovr, akey)

    contains

        subroutine s_mrg_build()

            integer :: ii, d2

            ext_max = 1; blo3 = huge(0); bhi3 = -huge(0)
            ii = head
            do while (ii /= 0)
                do d2 = 1, num_dims
                    ext_max = max(ext_max, ahi(d2, ii) - alo(d2, ii) + 1)
                    blo3(d2) = min(blo3(d2), alo(d2, ii)); bhi3(d2) = max(bhi3(d2), alo(d2, ii))
                end do
                ii = nxt(ii)
            end do
            nbmax = max(2, int(real(nlive)**(1.0/3.0)) + 1)*2
            cellw = max(cellw, ext_max + thr)
            do d2 = 1, num_dims
                rng = bhi3(d2) - blo3(d2) + 1
                if (rng > cellw*nbmax) cellw = (rng + nbmax - 1)/nbmax
            end do
            nbx = (bhi3(1) - blo3(1))/cellw + 1; nby = 1; nbz = 1
            if (n_glb > 0) nby = (bhi3(2) - blo3(2))/cellw + 1
            if (p_glb > 0) nbz = (bhi3(3) - blo3(3))/cellw + 1
            nb_tot = nbx*nby*nbz
            if (allocated(bh)) then
                if (size(bh) < nb_tot) deallocate (bh)
            end if
            if (.not. allocated(bh)) allocate (bh(nb_tot))
            if (.not. allocated(bc)) allocate (bc(size(nxt)))
            bh(1:nb_tot) = 0
            ii = head
            do while (ii /= 0)
                call s_mrg_binreg(ii)
                ii = nxt(ii)
            end do

        end subroutine s_mrg_build

        integer function f_mrg_binof(ii) result(bb)

            integer, intent(in) :: ii
            integer             :: b3(3)

            b3 = merge((alo(:,ii) - blo3)/cellw, 0, amr_dim)
            bb = 1 + b3(1) + nbx*(b3(2) + nby*b3(3))

        end function f_mrg_binof

        subroutine s_mrg_binreg(ii)

            integer, intent(in) :: ii
            integer             :: bb

            bb = f_mrg_binof(ii)
            bc(ii) = bh(bb)
            if (bh(bb) /= 0) bp(bh(bb)) = ii
            bp(ii) = 0
            bh(bb) = ii
            bidx(ii) = bb

        end subroutine s_mrg_binreg

        subroutine s_mrg_binun(ii)

            integer, intent(in) :: ii

            if (bp(ii) /= 0) then
                bc(bp(ii)) = bc(ii)
            else
                bh(bidx(ii)) = bc(ii)
            end if
            if (bc(ii) /= 0) bp(bc(ii)) = bp(ii)

        end subroutine s_mrg_binun

        subroutine s_mrg_fuse(x, y)

            integer, intent(in) :: x, y

            alo(:,x) = min(alo(:,x), alo(:,y)); ahi(:,x) = max(ahi(:,x), ahi(:,y))
            nxt(prv(y)) = nxt(y)
            if (nxt(y) /= 0) prv(nxt(y)) = prv(y)
            call s_mrg_binun(y)
            nlive = nlive - 1

        end subroutine s_mrg_fuse

        integer function f_mrg_qminj(ii) result(best)

            integer, intent(in) :: ii
            integer             :: b3(3), bx, by, bz, dx1, dy1, dz1, jj, d2, zl, zh, yl, yh
            logical             :: tc

            best = 0
            b3 = merge((alo(:,ii) - blo3)/cellw, 0, amr_dim); bx = b3(1); by = b3(2); bz = b3(3)
            zl = 0; zh = 0; yl = 0; yh = 0
            if (p_glb > 0) then; zl = max(0, bz - 1); zh = min(nbz - 1, bz + 1); end if
            if (n_glb > 0) then; yl = max(0, by - 1); yh = min(nby - 1, by + 1); end if
            do dz1 = zl, zh
                do dy1 = yl, yh
                    do dx1 = max(0, bx - 1), min(nbx - 1, bx + 1)
                        jj = bh(1 + dx1 + nbx*(dy1 + nby*dz1))
                        do while (jj /= 0)
                            if (jj > ii) then
                                tc = .true.
                                do d2 = 1, num_dims
                                    if (max(alo(d2, ii), alo(d2, jj)) - min(ahi(d2, ii), ahi(d2, jj)) - 1 >= thr) tc = .false.
                                end do
                                if (tc .and. (best == 0 .or. jj < best)) best = jj
                            end if
                            jj = bc(jj)
                        end do
                    end do
                end do
            end do

        end function f_mrg_qminj

        integer function f_mrg_qmina(ii) result(best)

            integer, intent(in) :: ii
            integer             :: b3(3), bx, by, bz, dx1, dy1, dz1, jj, d2, zl, zh, yl, yh
            logical             :: tc

            best = 0
            b3 = merge((alo(:,ii) - blo3)/cellw, 0, amr_dim); bx = b3(1); by = b3(2); bz = b3(3)
            zl = 0; zh = 0; yl = 0; yh = 0
            if (p_glb > 0) then; zl = max(0, bz - 1); zh = min(nbz - 1, bz + 1); end if
            if (n_glb > 0) then; yl = max(0, by - 1); yh = min(nby - 1, by + 1); end if
            do dz1 = zl, zh
                do dy1 = yl, yh
                    do dx1 = max(0, bx - 1), min(nbx - 1, bx + 1)
                        jj = bh(1 + dx1 + nbx*(dy1 + nby*dz1))
                        do while (jj /= 0)
                            if (jj < ii) then
                                tc = .true.
                                do d2 = 1, num_dims
                                    if (max(alo(d2, ii), alo(d2, jj)) - min(ahi(d2, ii), ahi(d2, jj)) - 1 >= thr) tc = .false.
                                end do
                                if (tc .and. (best == 0 .or. jj < best)) best = jj
                            end if
                            jj = bc(jj)
                        end do
                    end do
                end do
            end do

        end function f_mrg_qmina

    end subroutine s_amr_cluster

end module m_amr_cluster
