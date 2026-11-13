!>
!!@file
!!@brief Contains module m_amr_nest

#:include 'macros.fpp'

!> @brief Hierarchical child nesting for AMR regrid: each level's children are tagged by the density-gradient sensor on the old
!! parent-level fine solution, clustered, and replayed in rank-invariant order. Also holds the shared tag sensor.
module m_amr_nest

#ifdef MFC_MPI
    use mpi
#endif

    use m_derived_types
    use m_global_parameters
    use m_amr_cluster, only: s_amr_cluster
    use m_amr, only: amr_slots, amr_cons_st, amr_loc_of, amr_my_blk, amr_n_my, s_amr_refresh_my_blocks, amr_maxc_fit, &
        & amr_cpat_mar, s_amr_body_bbox, s_amr_expand_box_over_bodies, s_amr_tile_box, f_amr_boxes_overlap

    implicit none

    private
    public :: s_amr_regrid_nest_children, s_amr_nest_window, f_amr_tag_grad_sf

    !> Nested children were dropped for lack of block capacity (warn once per regrid)
    logical :: amr_child_dropped = .false.

contains

    !> Grow the pack buffers sidx/skb geometrically to hold nloc+extra entries, preserving the first nloc.
    impure subroutine s_amr_grow_pack(sidx, skb, nloc, extra)

        integer(8), allocatable, intent(inout) :: sidx(:)
        integer, allocatable, intent(inout)    :: skb(:)
        integer, intent(in)                    :: nloc, extra
        integer                                :: cap, newcap
        integer(8), allocatable                :: t8(:)
        integer, allocatable                   :: ti(:)

        cap = 0
        if (allocated(sidx)) cap = size(sidx)
        if (nloc + extra <= cap) return
        newcap = max(2*cap, max(nloc + extra, 1024))
        allocate (t8(newcap), ti(newcap))
        if (nloc > 0) then
            t8(1:nloc) = sidx(1:nloc)
            ti(1:nloc) = skb(1:nloc)
        end if
        call move_alloc(t8, sidx)
        call move_alloc(ti, skb)

    end subroutine s_amr_grow_pack

    !> Append this rank's tagged cells of child window [mlo:mhi] as (int8 linear index, parent kb) pairs; the encode matches the
    !! s_amr_nest_window_tags decode. gwin is not modified.
    impure subroutine s_amr_pack_gwin_pairs(gwin, mlo, mhi, mg, ng, kb, sidx, skb, nloc)

        integer, intent(in)                    :: mlo(3), mhi(3), mg, ng, kb
        logical, intent(in)                    :: gwin(mlo(1):,mlo(2):,mlo(3):)
        integer(8), allocatable, intent(inout) :: sidx(:)
        integer, allocatable, intent(inout)    :: skb(:)
        integer, intent(inout)                 :: nloc
        integer                                :: gi, gj, gk

        do gk = mlo(3), mhi(3)
            do gj = mlo(2), mhi(2)
                do gi = mlo(1), mhi(1)
                    if (.not. gwin(gi, gj, gk)) cycle
                    call s_amr_grow_pack(sidx, skb, nloc, 1)
                    nloc = nloc + 1
                    sidx(nloc) = int(gi, 8) + int(mg + 1, 8)*(int(gj, 8) + int(ng + 1, 8)*int(gk, 8))
                    skb(nloc) = kb
                end do
            end do
        end do

    end subroutine s_amr_pack_gwin_pairs

    !> Regrid phase 3b: nest child boxes inside each level-(l-1) box for l = 2..amr_max_level (parents first) and set box_level. A
    !! child's extent is the density-gradient sensor on the old parent-level fine solution, coarsened and clustered; a new region
    !! with no fine data gets a centred inset.
    !!
    !! Per level: pass 1 tags each parent's window from owned old blocks and packs (index, kb) pairs; an alltoall routes them to the
    !! parent's round-robin owner; pass 2 rebuilds the window and clusters it into children; an allgatherv plus a stable sort by kb
    !! replays the children in serial order, so the list and its truncation are rank-invariant. IB nesting is np=1 only.
    impure subroutine s_amr_regrid_nest_children(boxes, nboxes, box_level)

        type(t_box), allocatable, intent(inout) :: boxes(:)
        integer, intent(inout)                  :: nboxes, box_level(:)
        type(t_box), allocatable                :: grown(:)
        integer                                 :: lev, plo, phi, newlo, kb, ob, jch, nloc_send, ntot_g, nmych, ntot_ch, nct
        integer                                 :: mlo(3), mhi(3)
        integer, allocatable                    :: mlo_all(:,:), mhi_all(:,:), powner(:), mych(:,:), gch(:,:), chord(:)
        integer, allocatable                    :: skb(:), gkb(:), ctags(:,:)
        integer(8), allocatable                 :: sidx(:), gidx(:)
        logical, allocatable                    :: mine(:), covered(:), gwin(:,:,:)

        box_level(1:nboxes) = 1
        if (amr_max_level < 2) return
        amr_child_dropped = .false.

        ! grow to amr_max_blocks for the appends (the IB path leaves boxes at the cluster count)
        if (size(boxes) < amr_max_blocks) then
            allocate (grown(amr_max_blocks)); grown(1:nboxes) = boxes(1:nboxes); call move_alloc(grown, boxes)
        end if
        ! host-refresh the old blocks: the sensor reads the flat store on the host, stale until the stash's update
        do ob = 1, amr_num_blocks
            if (.not. amr_owns_all(ob)) cycle
            $:GPU_UPDATE(host='[amr_cons_st(:, :, :, :, amr_loc_of(ob))]')
        end do
        call s_amr_refresh_my_blocks()

        plo = 1; phi = nboxes  ! parents
        do lev = 2, amr_max_level
            if (phi < plo) exit
            newlo = nboxes + 1
            allocate (covered(plo:phi), mine(plo:phi), mlo_all(3,plo:phi), mhi_all(3,plo:phi), powner(plo:phi), mych(7, &
                      & amr_max_blocks))
            call s_amr_nest_mark_parents(boxes, plo, phi, lev, mine, covered)
            do kb = plo, phi
                powner(kb) = mod(kb - plo, max(num_procs, 1))
                call s_amr_nest_window(boxes(kb), mlo_all(:,kb), mhi_all(:,kb))
            end do

            ! pass 1: build a window only where this rank holds an overlapping old block (IB: the owner also builds one for the
            ! body)
            nloc_send = 0
            do kb = plo, phi
                mlo = mlo_all(:,kb); mhi = mhi_all(:,kb)
                if (f_amr_nest_window_empty(mlo, mhi)) cycle
                if (.not. (mine(kb) .or. (ib .and. powner(kb) == proc_rank))) cycle
                allocate (gwin(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3)))
                call s_amr_nest_tag_parent(boxes(kb), lev, mlo, mhi, gwin, covered(kb))
                call s_amr_pack_gwin_pairs(gwin, mlo, mhi, m_glb, n_glb, kb, sidx, skb, nloc_send)
                deallocate (gwin)
            end do
            call s_amr_nest_route_pairs(powner, plo, sidx, skb, nloc_send, gidx, gkb, ntot_g)

            ! pass 2: own parents only
            nmych = 0
            do kb = plo, phi
                if (powner(kb) /= proc_rank) cycle
                if (nmych + 1 > amr_max_blocks) exit
                mlo = mlo_all(:,kb); mhi = mhi_all(:,kb)
                if (f_amr_nest_window_empty(mlo, mhi)) cycle
                call s_amr_nest_window_tags(kb, mlo, mhi, gidx, gkb, ntot_g, ctags, nct)
                if (.not. covered(kb)) then
                    call s_amr_nest_inset_child(boxes(kb), lev, kb, mych, nmych)
                else if (nct > 0) then
                    call s_amr_nest_cluster_children(ctags, nct, mlo, mhi, lev, kb, mych, nmych)
                end if
                deallocate (ctags)
            end do

            call s_amr_nest_gather_children(mych, nmych, plo, phi, gch, chord, ntot_ch)
            do jch = 1, ntot_ch
                if (nboxes + 1 > amr_max_blocks) exit  ! canonical order gives the same truncation
                nboxes = nboxes + 1
                boxes(nboxes)%lo = gch(1:3,chord(jch)); boxes(nboxes)%hi = gch(4:6,chord(jch)); box_level(nboxes) = lev
            end do
            deallocate (gch, chord, powner, mych, gidx, gkb, covered, mine, mlo_all, mhi_all)
            plo = newlo; phi = nboxes
        end do
        if (nboxes >= amr_max_blocks .and. proc_rank == 0) print '(A)', &
            & ' [amr] NOTE: block pool full during multi-level nesting; some boxes were not refined further'

    end subroutine s_amr_regrid_nest_children

    !> mine(kb): this rank holds a level-(lev-1) block overlapping parent kb; covered(kb): any rank does (one LOR; every block has
    !! one owner).
    impure subroutine s_amr_nest_mark_parents(boxes, plo, phi, lev, mine, covered)

        type(t_box), intent(in) :: boxes(:)
        integer, intent(in)     :: plo, phi, lev
        logical, intent(out)    :: mine(plo:), covered(plo:)
        integer                 :: obi, ob, kb

#ifdef MFC_MPI
        integer :: ierr
#endif

        mine = .false.
        do obi = 1, amr_n_my
            ob = amr_my_blk(obi)
            if (amr_block_level(ob) /= lev - 1) cycle
            do kb = plo, phi
                if (f_amr_boxes_overlap(boxes(kb)%lo, boxes(kb)%hi, amr_region_lo_all(:,ob), amr_region_hi_all(:, &
                    & ob))) mine(kb) = .true.
            end do
        end do
        covered = mine
#ifdef MFC_MPI
        if (num_procs > 1) call MPI_ALLREDUCE(MPI_IN_PLACE, covered, phi - plo + 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
#endif

    end subroutine s_amr_nest_mark_parents

    !> The nesting window of a parent: its box inset by amr_cpat_mar, so a child's ghost prolongation reads valid parent interior.
    pure subroutine s_amr_nest_window(box, mlo, mhi)

        type(t_box), intent(in) :: box
        integer, intent(out)    :: mlo(3), mhi(3)

        mlo = box%lo + merge(amr_cpat_mar, 0, amr_dim); mhi = box%hi - merge(amr_cpat_mar, 0, amr_dim)

    end subroutine s_amr_nest_window

    !> Too small to nest a child in some active dimension.
    pure logical function f_amr_nest_window_empty(lo, hi) result(e)

        integer, intent(in) :: lo(3), hi(3)

        e = any(amr_dim .and. hi < lo)

    end function f_amr_nest_window_empty

    !> Pass 1 for one parent: tag its window from owned level-(lev-1) blocks (amr_block_level still holds the old levels). IB: also
    !! mark the body bbox with margin max(amr_buf, 4) + amr_cpat_mar, so the clamp to the window keeps a full image-point stencil.
    impure subroutine s_amr_nest_tag_parent(box, lev, mlo, mhi, gwin, covered)

        type(t_box), intent(in) :: box
        integer, intent(in)     :: lev, mlo(3), mhi(3)
        logical, intent(out)    :: gwin(mlo(1):,mlo(2):,mlo(3):)
        logical, intent(inout)  :: covered
        integer                 :: obi, ob, ib_i, bb_lo(3), bb_hi(3)
        logical                 :: any_tag

        gwin = .false.; any_tag = .false.
        do obi = 1, amr_n_my
            ob = amr_my_blk(obi)
            if (amr_block_level(ob) /= lev - 1) cycle
            if (.not. f_amr_boxes_overlap(box%lo, box%hi, amr_region_lo_all(:,ob), amr_region_hi_all(:,ob))) cycle
            call s_amr_tag_child_from_fine(ob, mlo, mhi, gwin, any_tag)
        end do
        if (.not. ib) return
        do ib_i = 1, num_ibs
            call s_amr_body_bbox(ib_i, max(amr_buf, 4) + amr_cpat_mar, bb_lo, bb_hi)
            bb_lo = max(bb_lo, mlo); bb_hi = min(bb_hi, mhi)
            if (f_amr_nest_window_empty(bb_lo, bb_hi)) cycle
            covered = .true.
            gwin(bb_lo(1):bb_hi(1),bb_lo(2):bb_hi(2),bb_lo(3):bb_hi(3)) = .true.
        end do

    end subroutine s_amr_nest_tag_parent

    !> Route each parent's (index, kb) pairs to its clustering owner with a stable counting sort by owner. Serial: the send list is
    !! the gathered list.
    impure subroutine s_amr_nest_route_pairs(powner, plo, sidx, skb, nloc, gidx, gkb, ntot)

        integer, intent(in)                    :: plo, powner(plo:)
        integer(8), allocatable, intent(inout) :: sidx(:)
        integer, allocatable, intent(inout)    :: skb(:)
        integer, intent(in)                    :: nloc
        integer(8), allocatable, intent(out)   :: gidx(:)
        integer, allocatable, intent(out)      :: gkb(:)
        integer, intent(out)                   :: ntot

#ifdef MFC_MPI
        integer                 :: i, ip, ierr
        integer, allocatable    :: rcnt(:), rdsp(:), scnt(:), sdsp(:), phead(:), pord(:), tkb(:)
        integer(8), allocatable :: tidx(:)
#endif

        if (.not. allocated(sidx)) allocate (sidx(0), skb(0))  ! this rank owned no tags at this level
#ifdef MFC_MPI
        if (num_procs > 1) then
            allocate (rcnt(num_procs), rdsp(num_procs), scnt(num_procs), sdsp(num_procs), phead(num_procs), pord(max(nloc, 1)))
            allocate (tidx(max(nloc, 1)), tkb(max(nloc, 1)))
            phead = 0
            do i = 1, nloc
                ip = powner(skb(i)) + 1; phead(ip) = phead(ip) + 1
            end do
            scnt = phead
            sdsp(1) = 0
            do ip = 2, num_procs
                sdsp(ip) = sdsp(ip - 1) + scnt(ip - 1)
            end do
            phead = sdsp + 1
            do i = 1, nloc
                ip = powner(skb(i)) + 1; pord(phead(ip)) = i; phead(ip) = phead(ip) + 1
            end do
            do i = 1, nloc
                tidx(i) = sidx(pord(i)); tkb(i) = skb(pord(i))
            end do
            call MPI_ALLTOALL(scnt, 1, MPI_INTEGER, rcnt, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
            rdsp(1) = 0
            do ip = 2, num_procs
                rdsp(ip) = rdsp(ip - 1) + rcnt(ip - 1)
            end do
            ntot = rdsp(num_procs) + rcnt(num_procs)
            allocate (gidx(max(ntot, 1)), gkb(max(ntot, 1)))
            call MPI_ALLTOALLV(tidx, scnt, sdsp, MPI_INTEGER8, gidx, rcnt, rdsp, MPI_INTEGER8, MPI_COMM_WORLD, ierr)
            call MPI_ALLTOALLV(tkb, scnt, sdsp, MPI_INTEGER, gkb, rcnt, rdsp, MPI_INTEGER, MPI_COMM_WORLD, ierr)
            deallocate (sidx, skb)
            return
        end if
#endif
        call move_alloc(sidx, gidx); call move_alloc(skb, gkb)
        ntot = nloc

    end subroutine s_amr_nest_route_pairs

    !> Pass 2 for one parent: rebuild its dense window from the routed pairs with gkb == kb (dedups replicated tags) and extract the
    !! tags in (k,j,i) order. The int8 decode matches s_amr_pack_gwin_pairs.
    impure subroutine s_amr_nest_window_tags(kb, mlo, mhi, gidx, gkb, ntot, ctags, nct)

        integer, intent(in)               :: kb, mlo(3), mhi(3), ntot, gkb(:)
        integer(8), intent(in)            :: gidx(:)
        integer, allocatable, intent(out) :: ctags(:,:)
        integer, intent(out)              :: nct
        logical, allocatable              :: gwin(:,:,:)
        integer                           :: i, gi, gj, gk
        integer(8)                        :: jrem

        allocate (gwin(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3)))
        gwin = .false.
        do i = 1, ntot
            if (gkb(i) /= kb) cycle
            gk = int(gidx(i)/(int(m_glb + 1, 8)*int(n_glb + 1, 8)))
            jrem = gidx(i) - int(gk, 8)*int(m_glb + 1, 8)*int(n_glb + 1, 8)
            gj = int(jrem/int(m_glb + 1, 8))
            gi = int(jrem - int(gj, 8)*int(m_glb + 1, 8))
            gwin(gi, gj, gk) = .true.
        end do
        nct = count(gwin)
        allocate (ctags(3, max(nct, 1)))
        nct = 0
        do gk = mlo(3), mhi(3); do gj = mlo(2), mhi(2); do gi = mlo(1), mhi(1)
            if (gwin(gi, gj, gk)) then
                nct = nct + 1
                ctags(1, nct) = gi; ctags(2, nct) = gj; ctags(3, nct) = gk
            end if
        end do; end do; end do

    end subroutine s_amr_nest_window_tags

    !> Cluster a parent's tagged cells into children, padded by amr_buf and clamped into the nesting window. IB: expand over bodies
    !! and re-clamp (the window already holds the body plus margin).
    impure subroutine s_amr_nest_cluster_children(ctags, nct, mlo, mhi, lev, kb, mych, nmych)

        integer, intent(inout)   :: ctags(:,:), mych(:,:), nmych
        integer, intent(in)      :: nct, mlo(3), mhi(3), lev, kb
        type(t_box), allocatable :: cboxes(:)
        integer                  :: ncb, kc, clo(3), chi(3)

        call s_amr_cluster(ctags, nct, cboxes, ncb, .false.)
        do kc = 1, ncb
            clo = cboxes(kc)%lo; chi = cboxes(kc)%hi
            call s_amr_nest_clamp(clo, chi, mlo, mhi, amr_buf)
            if (ib) then
                call s_amr_expand_box_over_bodies(clo, chi)
                call s_amr_nest_clamp(clo, chi, mlo, mhi, 0)
            end if
            call s_amr_nest_emit(clo, chi, lev, kb, mych, nmych)
        end do
        if (allocated(cboxes)) deallocate (cboxes)

    end subroutine s_amr_nest_cluster_children

    !> Brand-new region (no old fine data yet): a centred inset so the child still appears this regrid.
    impure subroutine s_amr_nest_inset_child(box, lev, kb, mych, nmych)

        type(t_box), intent(in) :: box
        integer, intent(in)     :: lev, kb
        integer, intent(inout)  :: mych(:,:), nmych
        integer                 :: ins(3), clo(3), chi(3)

        ins = merge(max((box%hi - box%lo + 1)/4, amr_cpat_mar), 0, amr_dim)
        clo = box%lo + ins; chi = box%hi - ins
        if (f_amr_nest_window_empty(clo, chi)) return
        call s_amr_nest_emit(clo, chi, lev, kb, mych, nmych)

    end subroutine s_amr_nest_inset_child

    !> Pad a child by `pad` and clamp it into the nesting window; collapsed dims are pinned to 0.
    pure subroutine s_amr_nest_clamp(clo, chi, mlo, mhi, pad)

        integer, intent(inout) :: clo(3), chi(3)
        integer, intent(in)    :: mlo(3), mhi(3), pad

        clo = merge(max(clo - pad, mlo), 0, amr_dim); chi = merge(min(chi + pad, mhi), 0, amr_dim)

    end subroutine s_amr_nest_clamp

    !> Tile a child to the level's slot cap, amr_maxc_fit/amr_ref_ratio**(lev-1) coarse cells, and append the tiles to this rank's
    !! list.
    impure subroutine s_amr_nest_emit(clo, chi, lev, kb, mych, nmych)

        integer, intent(in)    :: clo(3), chi(3), lev, kb
        integer, intent(inout) :: mych(:,:), nmych
        type(t_box)            :: tiles(amr_max_blocks)
        integer                :: nt, capped, it

        nt = 0; capped = 0
        call s_amr_tile_box(clo, chi, tiles, nt, amr_max_blocks, capped, amr_maxc_fit/amr_ref_ratio**(lev - 1))
        do it = 1, nt
            if (nmych + 1 > amr_max_blocks) then
                ! refinement is dropped: warn once
                if (.not. amr_child_dropped) then
                    write (0, '(A,I0,A)') ' [amr] WARNING: nested refinement dropped, box capacity amr_max_blocks = ', &
                           & amr_max_blocks, ' reached; raise amr_max_blocks'
                    amr_child_dropped = .true.
                end if
                exit
            end if
            nmych = nmych + 1
            mych(1:3,nmych) = tiles(it)%lo; mych(4:6,nmych) = tiles(it)%hi; mych(7, nmych) = kb
        end do

    end subroutine s_amr_nest_emit

    !> Assemble every rank's children (7 ints each: lo, hi, parent kb) and order them by parent kb with a stable counting sort:
    !! chord(1:ntot) indexes gch in the canonical (kb ascending, emission) order.
    impure subroutine s_amr_nest_gather_children(mych, nmych, plo, phi, gch, chord, ntot)

        integer, intent(in)               :: mych(:,:), nmych, plo, phi
        integer, allocatable, intent(out) :: gch(:,:), chord(:)
        integer, intent(out)              :: ntot
        integer, allocatable              :: chhead(:)
        integer                           :: ich, jch, kb

#ifdef MFC_MPI
        integer              :: ip, ierr
        integer, allocatable :: rcnt(:), rdsp(:)

        if (num_procs > 1) then
            allocate (rcnt(num_procs), rdsp(num_procs))
            call MPI_ALLGATHER(nmych, 1, MPI_INTEGER, rcnt, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
            rdsp(1) = 0
            do ip = 2, num_procs
                rdsp(ip) = rdsp(ip - 1) + rcnt(ip - 1)
            end do
            ntot = rdsp(num_procs) + rcnt(num_procs)
            allocate (gch(7, max(ntot, 1)))
            rcnt = rcnt*7; rdsp = rdsp*7
            call MPI_ALLGATHERV(mych, nmych*7, MPI_INTEGER, gch, rcnt, rdsp, MPI_INTEGER, MPI_COMM_WORLD, ierr)
        else
#endif
            ntot = nmych
            allocate (gch(7, max(ntot, 1))); gch(:,1:ntot) = mych(:,1:ntot)
#ifdef MFC_MPI
        end if
#endif
        allocate (chhead(plo:phi), chord(max(ntot, 1)))
        chhead = 0
        do ich = 1, ntot
            chhead(gch(7, ich)) = chhead(gch(7, ich)) + 1
        end do
        jch = 1
        do kb = plo, phi
            ich = chhead(kb); chhead(kb) = jch; jch = jch + ich
        end do
        do ich = 1, ntot
            kb = gch(7, ich)
            chord(chhead(kb)) = ich; chhead(kb) = chhead(kb) + 1
        end do

    end subroutine s_amr_nest_gather_children

    !> Sensor-on-fine child tagging: OR the density-gradient tags of an old fine block into an L0 tag grid over a parent window.
    !! Reads the flat store on the host (caller refreshes it); one-sided differences at the fine edges avoid stale ghosts.
    impure subroutine s_amr_tag_child_from_fine(ob, win_lo, win_hi, ctag, any_tag)

        integer, intent(in)    :: ob, win_lo(3), win_hi(3)
        logical, intent(inout) :: ctag(win_lo(1):,win_lo(2):,win_lo(3):)
        logical, intent(inout) :: any_tag
        integer                :: rr, ci, cj, ck, fi, fj, fk, d1, d2, d3, fm1, fm2, fm3, olo(3), lo(3), hi(3)
        real(wp)               :: g
        logical                :: tagged

        rr = amr_slots(ob)%amr_ref_ratio
        olo = amr_region_lo_all(:,ob)
        fm1 = amr_slots(ob)%m; fm2 = amr_slots(ob)%n; fm3 = amr_slots(ob)%p
        lo = merge(max(win_lo, amr_region_lo_all(:,ob)), 0, amr_dim); hi = merge(min(win_hi, amr_region_hi_all(:,ob)), 0, amr_dim)
        do ck = lo(3), hi(3)
            do cj = lo(2), hi(2)
                do ci = lo(1), hi(1)
                    tagged = .false.
                    do d3 = 0, merge(rr - 1, 0, p_glb > 0)
                        fk = (ck - olo(3))*rr + d3
                        do d2 = 0, merge(rr - 1, 0, n_glb > 0)
                            fj = (cj - olo(2))*rr + d2
                            do d1 = 0, rr - 1
                                fi = (ci - olo(1))*rr + d1
                                g = f_amr_tag_grad_st(amr_loc_of(ob), fi, fj, fk, fm1, fm2, fm3)
                                if (g > amr_tag_eps) tagged = .true.
                            end do
                        end do
                    end do
                    if (tagged) then
                        ctag(ci, cj, ck) = .true.
                        any_tag = .true.
                    end if
                end do
            end do
        end do

    end subroutine s_amr_tag_child_from_fine

    !> Total density (sum of the continuity variables) at one cell; `_st` reads the flat store, `_sf` the level-0 field.
    #:for RSFX, RSRC in [('st', 'amr_cons_st'), ('sf', '')]
        #:set RARG = 'loc' if RSRC else 'q'
        #:set RHO = 'f_amr_rho_tot_' + RSFX
        pure function f_amr_rho_tot_${RSFX}$(${RARG}$, ci, cj, ck) result(r)

            #:if RSRC
                integer, intent(in) :: loc  !< flat-store slot
            #:else
                type(scalar_field), dimension(:), intent(in) :: q
            #:endif
            integer, intent(in) :: ci, cj, ck
            real(wp)            :: r
            integer             :: f

            r = 0._wp
            do f = eqn_idx%cont%beg, eqn_idx%cont%end
                #:if RSRC
                    r = r + real(amr_cons_st(ci, cj, ck, f, loc), wp)
                #:else
                    r = r + real(q(f)%sf(ci, cj, ck), wp)
                #:endif
            end do

        end function f_amr_rho_tot_${RSFX}$

        !> Density-gradient tag indicator at one cell, max over the active dimensions. `_st` clamps the stencil to the block
        !! (fm1:fm3 = block extents) since its ghosts are stale; `_sf` reads the level-0 halo.
        pure function f_amr_tag_grad_${RSFX}$(${RARG}$, ci, cj, ck${', fm1, fm2, fm3' if RSRC else ''}$) result(g)

            #:if RSRC
                integer, intent(in) :: loc, fm1, fm2, fm3
            #:else
                type(scalar_field), dimension(:), intent(in) :: q
            #:endif
            integer, intent(in) :: ci, cj, ck
            real(wp)            :: g, r0
            integer             :: ip(3), im(3)

            ip = [ci, cj, ck] + 1; im = [ci, cj, ck] - 1
            #:if RSRC
                ip = min(ip, [fm1, fm2, fm3]); im = max(im, 0)
            #:endif
            r0 = ${RHO}$(${RARG}$, ci, cj, ck)
            g = f_amr_tag_ind(${RHO}$(${RARG}$, ip(1), cj, ck), r0, ${RHO}$(${RARG}$, im(1), cj, ck))
            if (n_glb > 0) g = max(g, f_amr_tag_ind(${RHO}$(${RARG}$, ci, ip(2), ck), r0, ${RHO}$(${RARG}$, ci, im(2), ck)))
            if (p_glb > 0) g = max(g, f_amr_tag_ind(${RHO}$(${RARG}$, ci, cj, ip(3)), r0, ${RHO}$(${RARG}$, ci, cj, im(3))))

        end function f_amr_tag_grad_${RSFX}$
    #:endfor

    !> Tag indicator from the densities at i+1, i, i-1: central difference over 2*rho_i, or with amr_tag_abs the larger absolute
    !! face jump.
    pure real(wp) function f_amr_tag_ind(rp, r0, rm) result(v)
        real(wp), intent(in) :: rp, r0, rm

        if (amr_tag_abs) then
            v = max(abs(rp - r0), abs(r0 - rm))
        else
            v = abs(rp - rm)/(2._wp*max(abs(r0), 1.e-30_wp))
        end if

    end function f_amr_tag_ind

end module m_amr_nest
