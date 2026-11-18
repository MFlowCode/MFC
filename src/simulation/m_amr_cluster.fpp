!>
!!@file
!!@brief Contains module m_amr_cluster

#:include 'macros.fpp'

!> @brief Berger-Rigoutsos clustering of a sparse tag list into separated block boxes, identically on every rank: tag signatures,
!! bisection, canonical Morton order and min-separation merge.
module m_amr_cluster

#ifdef MFC_MPI
    use mpi
#endif

    use m_derived_types  ! t_box
    use m_box, only: f_morton
    use m_global_parameters
    use m_amr_state
    use m_amr_distribution

    implicit none

    private
    public :: s_amr_cluster

    !> A node is wide when its box spans more ranks than this; wide nodes use the batched collective, narrow ones reduce among their
    !! few overlapping ranks. 8 is one 2x2x2 brick of ranks.
    integer, parameter :: amr_cl_wide = 8
    !> Scratch of the narrow point-to-point reduction (pidx is all-zero between depths)
    integer, allocatable :: pidx(:), plist(:), scnt(:), rcnt(:), sdsp2(:), rdsp2(:), soff(:), roff(:), sbuf(:), rbuf(:), creq(:), &
                                 & wbuf(:)

contains

    !> Concatenated 1D tag signatures of box [blo0:bhi0] from tags ts:te: sig(off(d) + t - blo0(d)) counts the in-box tagged cells
    !! at position t along d. One signature serves the trim, the count and every split. nsig is the used length.
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

    !> Shrink box [blo:bhi] to the tight bbox of its tagged cells (first/last nonzero of each axis signature) and return their
    !! count; ok = .false. if none tagged.
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
        if (ntag == 0) return
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

    !> Bisect a trimmed box: on the longest splittable axis prefer the widest zero-signature hole, else the strongest signature
    !! inflection (Laplacian sign change). ok = .false. if no axis admits a split leaving both children >= min_child cells.
    pure subroutine s_amr_find_split_sig(sig, off, blo0, blo, bhi, sax, spos, ok)

        integer, intent(in)  :: sig(:), off(3), blo0(3)
        integer, intent(in)  :: blo(3), bhi(3)
        integer, intent(out) :: sax, spos
        logical, intent(out) :: ok
        !> smallest box the bisection may produce; amr_blocking_factor raises it to stop over-generation
        integer :: min_child
        integer :: axord(3), ext(3), d, ax, t, s, b
        integer :: run, run_start, best_run, best_start, lap, prevlap, bestmag, bestpos

        min_child = max(2, amr_blocking_factor)
        ok = .false.; sax = 0; spos = 0
        ext = bhi - blo + 1
        axord = [1, 2, 3]  ! axes by descending extent
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
            b = off(ax) - blo0(ax)  ! signature at position t is sig(b + t)
            ! widest zero run (interior, since the box is trimmed)
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
            ! strongest inflection: Laplacian sign change with the largest jump
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

    !> Partition wt(:, ts:te) in place by coord(sax) < spos; the low child is ts:lo-1 and the high child lo:te.
    pure subroutine s_cl_partition(wt, ts, te, sax, spos, lo)

        integer, intent(inout) :: wt(:,:)
        integer, intent(in)    :: ts, te, sax, spos
        integer, intent(out)   :: lo
        integer                :: hi, tmp(3)

        lo = ts; hi = te
        do while (lo <= hi)
            if (wt(sax, lo) < spos) then
                lo = lo + 1
            else
                tmp = wt(:,lo); wt(:,lo) = wt(:,hi); wt(:,hi) = tmp
                hi = hi - 1
            end if
        end do

    end subroutine s_cl_partition

#ifdef MFC_MPI
    !> Nonblocking exchange: receive rc(q) words into rb at rd(q) and send sc(q) words from sb at sd(q), per peer plist(q).
    impure subroutine s_cl_exchange(tagc, np2, sb, sc, sd, rb, rc, rd)

        integer, intent(in)    :: tagc, np2, sc(:), sd(:), rc(:), rd(:)
        integer, intent(inout) :: sb(:), rb(:)
        integer                :: q, nreq2, ierr

        nreq2 = 0
        do q = 1, np2
            if (rc(q) > 0) then
                nreq2 = nreq2 + 1
                call MPI_IRECV(rb(rd(q) + 1), rc(q), MPI_INTEGER, plist(q), tagc, MPI_COMM_WORLD, creq(nreq2), ierr)
            end if
        end do
        do q = 1, np2
            if (sc(q) > 0) then
                nreq2 = nreq2 + 1
                call MPI_ISEND(sb(sd(q) + 1), sc(q), MPI_INTEGER, plist(q), tagc, MPI_COMM_WORLD, creq(nreq2), ierr)
            end if
        end do
        if (nreq2 > 0) call MPI_WAITALL(nreq2, creq, MPI_STATUSES_IGNORE, ierr)

    end subroutine s_cl_exchange

    !> Register rank rr as a peer of the narrow reduction if not yet listed.
    impure subroutine s_cl_add_peer(rr, np2)

        integer, intent(in)    :: rr
        integer, intent(inout) :: np2

        if (pidx(rr) == 0) then; np2 = np2 + 1; plist(np2) = rr; pidx(rr) = np2; end if

    end subroutine s_cl_add_peer

    !> Sum the depth's batched signatures across ranks. Wide nodes (every rank overlaps) ride one ALLREDUCE. Narrow nodes reduce
    !! point-to-point within their small rank brick: members send to the root bovr(1,:), which sums and sends the total back. Both
    !! ends enumerate nodes in the same global order, so one aggregated message per peer per phase matches without negotiation.
    impure subroutine s_cl_reduce_depth(nbat, bsig, bofs, blen, bnov, bwide, bovr)

        integer, intent(in)    :: nbat, bofs(:), blen(:), bnov(:), bovr(:,:)
        integer, intent(inout) :: bsig(:)
        logical, intent(in)    :: bwide(:)
        integer                :: j, t, q, rr, o1, nwb, np2, nsnd, nrcv, tagc, ierr

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
        ! peers: whoever roots a node I hold, plus whoever holds a node I root
        np2 = 0
        do j = 1, nbat
            if (bwide(j)) cycle
            if (bovr(1, j) == proc_rank) then
                do t = 2, bnov(j)
                    call s_cl_add_peer(bovr(t, j), np2)
                end do
            else
                call s_cl_add_peer(bovr(1, j), np2)
            end if
        end do
        if (np2 == 0) return
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
        ! phase A: members send their contribution to the root
        soff(1:np2) = sdsp2(1:np2)
        do j = 1, nbat
            if (bwide(j) .or. bovr(1, j) == proc_rank) cycle
            q = pidx(bovr(1, j))
            sbuf(soff(q) + 1:soff(q) + blen(j)) = bsig(bofs(j) + 1:bofs(j) + blen(j)); soff(q) = soff(q) + blen(j)
        end do
        tagc = amr_tag_narrow + int(mod(amr_mesh_epoch, 50_8))
        call s_cl_exchange(tagc, np2, sbuf, scnt, sdsp2, rbuf, rcnt, rdsp2)
        ! integer SUM is exact and order-independent, so the total equals a machine-wide reduction
        roff(1:np2) = rdsp2(1:np2)
        do j = 1, nbat
            if (bwide(j) .or. bovr(1, j) /= proc_rank) cycle
            do t = 2, bnov(j)
                q = pidx(bovr(t, j))
                bsig(bofs(j) + 1:bofs(j) + blen(j)) = bsig(bofs(j) + 1:bofs(j) + blen(j)) + rbuf(roff(q) + 1:roff(q) + blen(j))
                roff(q) = roff(q) + blen(j)
            end do
        end do
        ! phase B: the total goes back down; counts mirror phase A, so the buffers swap roles
        roff(1:np2) = rdsp2(1:np2)
        do j = 1, nbat
            if (bwide(j) .or. bovr(1, j) /= proc_rank) cycle
            do t = 2, bnov(j)
                q = pidx(bovr(t, j))
                rbuf(roff(q) + 1:roff(q) + blen(j)) = bsig(bofs(j) + 1:bofs(j) + blen(j)); roff(q) = roff(q) + blen(j)
            end do
        end do
        call s_cl_exchange(tagc + 50, np2, rbuf, rcnt, rdsp2, sbuf, scnt, sdsp2)
        soff(1:np2) = sdsp2(1:np2)
        do j = 1, nbat
            if (bwide(j) .or. bovr(1, j) == proc_rank) cycle
            q = pidx(bovr(1, j))
            bsig(bofs(j) + 1:bofs(j) + blen(j)) = sbuf(soff(q) + 1:soff(q) + blen(j)); soff(q) = soff(q) + blen(j)
        end do
        do q = 1, np2  ! clear only what was touched: a full wipe would be O(P) per depth
            pidx(plist(q)) = 0
        end do

    end subroutine s_cl_reduce_depth

    !> Union the per-rank accepted boxes (rank-local subtrees are walked only by their owner). Rank order is not serial acceptance
    !! order; the Morton sort that follows removes the dependence. alo/ahi grow to the union (the cap applies to the merged set).
    impure subroutine s_cl_union(alo, ahi, nacc)

        integer, allocatable, intent(inout) :: alo(:,:), ahi(:,:)
        integer, intent(inout)              :: nacc
        integer, allocatable                :: gcnt(:), gdsp(:), sbx(:,:), gbx(:,:)
        integer                             :: i, ntot, ierr

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
        if (ntot > size(alo, 2)) then
            deallocate (alo, ahi)
            allocate (alo(3, ntot), ahi(3, ntot))
        end if
        nacc = ntot
        do i = 1, nacc
            alo(:,i) = gbx(1:3,i); ahi(:,i) = gbx(4:6,i)
        end do
        deallocate (gcnt, gdsp, sbx, gbx)

    end subroutine s_cl_union
#endif

    !> Level-order Berger-Rigoutsos descent: accepts boxes into alo/ahi(:, 1:nacc). Every shared node of a depth rides one
    !! reduction, so the collective count is O(tree depth). A child's overlapping ranks are a subset of its parent's, so every rank
    !! walks the whole shared subtree in the same order. Children are appended past ncur and shifted down when the depth closes
    !! (peak occupancy 3*ncur <= the 4*cap + 8 the arrays carry).
    impure subroutine s_cl_bisect(tags, ntag_in, reduce, alo, ahi, nacc, capped)

        integer, intent(in)  :: tags(:,:), ntag_in
        logical, intent(in)  :: reduce
        integer, intent(out) :: alo(:,:), ahi(:,:), nacc
        logical, intent(out) :: capped
        integer, allocatable :: slo(:,:), shi(:,:), sts(:), ste(:), wt(:,:), sdep(:), sig(:), ovr(:)
        integer, allocatable :: kpos(:), kbat(:), bofs(:), blen(:), boff(:,:), bsig(:), bnov(:), bovr(:,:)
        logical, allocatable :: bwide(:)
        integer              :: cap, mg, ng, pg, t, i, j, ncur, nnxt, nkeep, nbat, nbuf, novr, nsig, ntag, dep
        integer              :: blo0(3), bhi0(3), off(3), blo(3), bhi(3), ts, te, lo, sax, spos
        integer(8)           :: vol  !< a global-bbox first pass can exceed 2**31 cells
        logical              :: ok, force, mine
        real(wp)             :: eff

        mg = m_glb; ng = 0; pg = 0
        if (n_glb > 0) ng = n_glb
        if (p_glb > 0) pg = p_glb
        cap = amr_max_blocks
        allocate (slo(3, 4*cap + 8), shi(3, 4*cap + 8), sts(4*cap + 8), ste(4*cap + 8), wt(3, ntag_in), sdep(4*cap + 8))
        allocate (sig(mg + ng + pg + 3))  ! the three full domain extents; reused by every node
        allocate (ovr(amr_cl_wide))  ! only narrow nodes are enumerated, so this does not size with P
        allocate (kpos(4*cap + 8), kbat(4*cap + 8), bofs(4*cap + 8), blen(4*cap + 8), boff(3, 4*cap + 8))
        allocate (bsig(4*(mg + ng + pg + 3)), bnov(4*cap + 8), bwide(4*cap + 8), bovr(amr_cl_wide, 4*cap + 8))
        allocate (pidx(0:max(num_procs - 1, 0)), plist(max(num_procs, 1)))
        allocate (scnt(max(num_procs, 1)), rcnt(max(num_procs, 1)), sdsp2(max(num_procs, 1)), rdsp2(max(num_procs, 1)))
        allocate (soff(max(num_procs, 1)), roff(max(num_procs, 1)))
        allocate (wbuf(1), sbuf(1), rbuf(1), creq(1))
        pidx = 0
        ! working copy of the tags, partitioned in place as the tree descends
        do t = 1, ntag_in
            wt(:,t) = tags(:,t)
        end do
        ncur = 1; slo(:,1) = [0, 0, 0]; shi(:,1) = [mg, ng, pg]  ! first node trims to the global tagged bbox
        sts(1) = 1; ste(1) = ntag_in; sdep(1) = 0
        nacc = 0; capped = .false.
        do while (ncur > 0)
            nkeep = 0; nbat = 0; nbuf = 0; nnxt = 0
            ! pass 1: classify nodes and stash the signatures that need reducing (rank-local nodes recompute theirs in pass 2)
            do i = 1, ncur
                blo0 = slo(:,i); bhi0 = shi(:,i)
                novr = f_amr_overlap_count(blo0, bhi0)
                mine = (num_procs == 1) .or. f_amr_rank_overlaps(blo0, bhi0, proc_rank)
                ! A rank holds tags only in its subdomain, so a narrow node it does not reach would add only zeros: drop the
                ! subtree. Wide nodes are never dropped (every rank must enter their collective with the identical batch).
                if (reduce .and. num_procs > 1 .and. .not. mine .and. novr <= amr_cl_wide) cycle
                call s_amr_box_sig(wt, sts(i), ste(i), blo0, bhi0, sig, off, nsig)
                nkeep = nkeep + 1; kpos(nkeep) = i; kbat(nkeep) = 0
                if (reduce .and. num_procs > 1 .and. novr > 1) then
                    call s_amr_size_int(bsig, nbuf + nsig)
                    nbat = nbat + 1; kbat(nkeep) = nbat
                    bofs(nbat) = nbuf; blen(nbat) = nsig; boff(:,nbat) = off
                    bnov(nbat) = novr; bwide(nbat) = (novr > amr_cl_wide)
                    bovr(1, nbat) = -1  ! defined for wide nodes too: Fortran does not promise .or. short-circuits
                    if (.not. bwide(nbat)) then
                        call s_amr_ranks_overlapping(blo0, bhi0, ovr, novr)
                        bovr(1:novr,nbat) = ovr(1:novr)
                    end if
                    bsig(nbuf + 1:nbuf + nsig) = sig(1:nsig)
                    nbuf = nbuf + nsig
                end if
            end do
#ifdef MFC_MPI
            call s_cl_reduce_depth(nbat, bsig, bofs, blen, bnov, bwide, bovr)
#endif
            ! pass 2: trim, accept or split every kept node
            do j = 1, nkeep
                i = kpos(j); blo = slo(:,i); bhi = shi(:,i); ts = sts(i); te = ste(i); dep = sdep(i)
                blo0 = blo; bhi0 = bhi
                if (kbat(j) > 0) then
                    off = boff(:,kbat(j))
                    sig(1:blen(kbat(j))) = bsig(bofs(kbat(j)) + 1:bofs(kbat(j)) + blen(kbat(j)))
                    nsig = blen(kbat(j))
                else
                    ! rank-local node: recompute (safe, pass 2 only partitions this node's own wt range)
                    call s_amr_box_sig(wt, ts, te, blo0, bhi0, sig, off, nsig)
                end if
                call s_amr_trim_from_sig(sig(1:nsig), off, blo0, bhi0, blo, bhi, ok, ntag)
                if (.not. ok) cycle
                vol = 1_8
                do t = 1, num_dims; vol = vol*int(bhi(t) - blo(t) + 1, 8); end do
                eff = real(ntag, wp)/real(max(vol, 1_8), wp)
                call s_amr_find_split_sig(sig(1:nsig), off, blo0, blo, bhi, sax, spos, ok)
                ! guard against overflowing the amr_max_blocks cap: pending = rest of this depth + children queued so far
                force = (nacc + (nkeep - j) + nnxt + 1 >= cap)
                if (eff >= amr_cluster_eff .or. .not. ok .or. force) then
                    if (nacc < cap) then; nacc = nacc + 1; alo(:,nacc) = blo; ahi(:,nacc) = bhi; end if
                    if (force .and. ok .and. eff < amr_cluster_eff) capped = .true.
                else
                    call s_cl_partition(wt, ts, te, sax, spos, lo)
                    slo(:,ncur + nnxt + 1) = blo; shi(:,ncur + nnxt + 1) = bhi; shi(sax, ncur + nnxt + 1) = spos - 1
                    sts(ncur + nnxt + 1) = ts; ste(ncur + nnxt + 1) = lo - 1; sdep(ncur + nnxt + 1) = dep + 1
                    slo(:,ncur + nnxt + 2) = blo; shi(:,ncur + nnxt + 2) = bhi; slo(sax, ncur + nnxt + 2) = spos
                    sts(ncur + nnxt + 2) = lo; ste(ncur + nnxt + 2) = te; sdep(ncur + nnxt + 2) = dep + 1
                    nnxt = nnxt + 2
                end if
            end do
            do i = 1, nnxt  ! the children become the next depth
                slo(:,i) = slo(:,ncur + i); shi(:,i) = shi(:,ncur + i)
                sts(i) = sts(ncur + i); ste(i) = ste(ncur + i); sdep(i) = sdep(ncur + i)
            end do
            ncur = nnxt
        end do
        deallocate (pidx, plist, scnt, rcnt, sdsp2, rdsp2, soff, roff, wbuf, sbuf, rbuf, creq)

    end subroutine s_cl_bisect

    !> Sort boxes by the Morton key of lo so the merge input depends on the box set alone, not on the rank-dependent acceptance
    !! order. Accepted boxes are disjoint, so the keys are distinct; the sort is stable for any collision.
    impure subroutine s_cl_morton_sort(alo, ahi, nacc)

        integer, intent(inout)  :: alo(:,:), ahi(:,:)
        integer, intent(in)     :: nacc
        integer(8), allocatable :: akey(:)
        integer, allocatable    :: sperm(:), t2lo(:,:), t2hi(:,:)
        integer                 :: i

        allocate (akey(nacc), sperm(nacc), t2lo(3, nacc), t2hi(3, nacc))
        do i = 1, nacc
            akey(i) = f_morton(alo(1, i), alo(2, i), alo(3, i))
        end do
        call s_amr_sort_by_key(akey, nacc, sperm)
        do i = 1, nacc
            t2lo(:,i) = alo(:,sperm(i)); t2hi(:,i) = ahi(:,sperm(i))
        end do
        alo(:,1:nacc) = t2lo; ahi(:,1:nacc) = t2hi

    end subroutine s_cl_morton_sort

    !> Min-separation merge: two boxes are separated only if some active dim's gap reaches thr; else fuse to their bounding box. The
    !! fusion sequence is the linear walk (scan in list order, fuse the first too-close pair, restart), reproduced exactly with O(1)
    !! unlinking and a binned candidate search: a too-close partner of i lies within the 3^d bins around i's lo when the bin width
    !! is ext_max + thr, and the minimum surviving index among candidates is the first one a linear walk meets. After fusing (i, j)
    !! only box i changed, so the walk re-tests earlier survivors (minimum index first), then later ones, then resumes at the
    !! survivor's successor. Bins are rebuilt when a grown box outgrows the cell width.
    impure subroutine s_cl_merge(alo, ahi, nacc)

        integer, intent(inout) :: alo(:,:), ahi(:,:), nacc
        integer, allocatable   :: nxt(:)          !< survivor list in canonical order
        integer, allocatable   :: prv(:)          !< predecessor links
        integer, allocatable   :: bp(:), bidx(:)  !< bin back-links and current bin of each live box
        integer, allocatable   :: bh(:), bc(:)    !< bin heads and per-box chains
        integer                :: head, nlive, thr, i, k, d, dirty, aa, jb2, jbest, extd
        integer                :: blo3(3), bhi3(3), nbmax, rng, ext_max, cellw, nbx, nby, nbz, nb_tot
        logical                :: need_build

        thr = buff_size + 2*amr_buf
        allocate (nxt(max(nacc, 1)), prv(max(nacc, 1)), bp(max(nacc, 1)), bidx(max(nacc, 1)))
        do i = 1, nacc - 1
            nxt(i) = i + 1
        end do
        if (nacc >= 1) nxt(nacc) = 0
        do i = 1, nacc
            prv(i) = i - 1
        end do
        ! head must be 0 when nothing is tagged (nxt(1) is then unwritten)
        head = merge(1, 0, nacc >= 1); nlive = nacc
        cellw = 0
        need_build = .true.
        i = head
        outer: do while (i /= 0)
            if (need_build) then
                call s_mrg_build()
                need_build = .false.
            end if
            jbest = f_mrg_qmin(i, .true.)
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
                    aa = f_mrg_qmin(dirty, .false.)
                    if (aa /= 0) then
                        call s_mrg_fuse(aa, dirty)
                        dirty = aa
                        cycle chain
                    end if
                    jb2 = f_mrg_qmin(dirty, .true.)
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
        k = 0  ! compact in list order
        i = head
        do while (i /= 0)
            k = k + 1
            if (k /= i) then
                alo(:,k) = alo(:,i); ahi(:,k) = ahi(:,i)
            end if
            i = nxt(i)
        end do
        nacc = nlive

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

        !> Minimum-index too-close partner of box ii among the binned candidates above (above = .true.) or below it.
        integer function f_mrg_qmin(ii, above) result(best)

            integer, intent(in) :: ii
            logical, intent(in) :: above
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
                            if (merge(jj > ii, jj < ii, above)) then
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

        end function f_mrg_qmin

    end subroutine s_cl_merge

    !> Cluster a sparse tag list tags(1:3, 1:ntag_in) of level-0 cell coords into separated block boxes, identically on every rank.
    !! Berger-Rigoutsos bisection until a box's tag efficiency reaches amr_cluster_eff (or it is atomic or the amr_max_blocks cap is
    !! hit), then merge boxes whose amr_buf-padded extents come within buff_size. Boxes are raw tagged extents; the caller pads,
    !! clamps and size-caps each. reduce: tags is this rank's local list and node signatures are reduced across ranks; otherwise
    !! tags is already replicated.
    impure subroutine s_amr_cluster(tags, ntag_in, boxes, nboxes, reduce)

        integer, intent(in)                   :: tags(:,:), ntag_in
        logical, intent(in)                   :: reduce
        type(t_box), allocatable, intent(out) :: boxes(:)
        integer, intent(out)                  :: nboxes
        integer, allocatable                  :: alo(:,:), ahi(:,:)
        integer                               :: nacc, i
        logical                               :: capped

        nboxes = 0
        ! with reduce, a rank without tags must still walk the tree and enter every collective
        if (.not. reduce .and. ntag_in == 0) return
        allocate (alo(3, amr_max_blocks), ahi(3, amr_max_blocks))
        call s_cl_bisect(tags, ntag_in, reduce, alo, ahi, nacc, capped)
#ifdef MFC_MPI
        if (reduce .and. num_procs > 1) call s_cl_union(alo, ahi, nacc)
#endif
        call s_cl_morton_sort(alo, ahi, nacc)
        call s_cl_merge(alo, ahi, nacc)
        if (capped .and. proc_rank == 0) print '(A,I0)', ' [amr] WARNING: tag clustering capped at amr_max_blocks = ', &
            & amr_max_blocks
        if (nacc > amr_max_blocks) then  ! boxes past the cap never refine
            if (proc_rank == 0) print '(A,I0,A,I0)', ' [amr] WARNING: merged box set truncated: ', nacc, ' boxes, keeping ', &
                & amr_max_blocks
            nacc = amr_max_blocks
        end if
        nboxes = nacc
        allocate (boxes(nboxes))
        do i = 1, nboxes
            boxes(i)%lo = alo(:,i); boxes(i)%hi = ahi(:,i)
        end do

    end subroutine s_amr_cluster

end module m_amr_cluster
