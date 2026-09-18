!>
!!@file
!!@brief Contains module m_amr_regrid

#:include 'macros.fpp'

!> @brief Dynamic regrid for the block-structured AMR level set: density-gradient tagging, box shaping (pad/clamp/tile/IB merge)
!! around the clusterer (m_amr_cluster), hierarchical child nesting, slot rebuild with cross-rank fine-state migration. Block/slot
!! state lives in m_amr (and m_global_parameters); this module only drives it.
module m_amr_regrid

#ifdef MFC_MPI
    use mpi  !< per-node signature ALLREDUCE for the rank-invariant clustering, point-to-point for the fine-state migration
#endif

    use m_derived_types  ! scalar_field, t_box
    use m_box, only: f_morton
    use m_global_parameters
    use m_constants, only: mapCells
    use m_mpi_proxy, only: s_mpi_abort
    use m_mpi_common, only: s_mpi_allreduce_min, s_mpi_allreduce_max
    use m_amr_wave
    use m_amr_cluster, only: s_amr_cluster
    use m_amr, only: amr_slots, amr_cons_st, amr_stor_st, amr_loc_of, amr_slot_live, amr_my_blk, amr_n_my, &
        & s_amr_refresh_my_blocks, amr_maxc_fit, amr_seam_pairs_dirty, amr_mesh_epoch, amr_cpat_mar, s_amr_alloc_slot, &
        & s_amr_alloc_slot_stash, s_amr_prereserve_stash, s_amr_free_slot, s_amr_reduce_xchg_flag, s_amr_reconcile_slots, &
        & s_amr_assign_block_owners, s_amr_exchange_coarse_cons_halo, s_lag_phys_to_cells, s_amr_body_bbox, &
        & s_amr_expand_box_over_bodies, s_amr_tile_box, f_amr_seam_dim, f_amr_boxes_overlap, s_set_amr_fine_geometry, &
        & s_interpolate_coarse_to_fine, s_amr_setup_ib, f_l0_slot, amr_cad_tot, amr_cad_esc, amr_cad_armed, &
        & s_amr_ranks_overlapping, f_amr_overlap_count, f_amr_rank_overlaps, s_amr_l1_fill_exchange, s_amr_l1_fill_consume, &
        & s_amr_parent_fill_exchange, s_amr_parent_fill_consume, s_amr_fill_wave_done, amr_fw_sq, amr_fw_rq, amr_fw_dev
    use m_amr_xchg_audit, only: XA_F4_SND, XA_F4_RCV
    use m_acoustic_src, only: acoustic_supp_lo, acoustic_supp_hi
    use m_active_box, only: ab_x, ab_y, ab_z, ab_active
    use m_bubbles_EL, only: s_lag_cloud_bbox_local

    implicit none

    private
    public :: s_amr_regrid, s_amr_check_seam_topology, s_amr_check_active_box_containment

    !> Lagrangian bubble-cloud exclusion support: padded global coarse-index bbox (positions + mapCells smearing + stencil headroom
    !! [+ drift margin at regrid]). Blocks and regrid boxes stay clear: a bubble inside a block loses two-way coupling (fine advance
    !! skips the EL hooks, restriction discards the coarse result under the block). Recomputed collectively each regrid; guarded
    !! rank-locally per stage.
    integer :: lag_supp_lo(3), lag_supp_hi(3)
    logical :: lag_supp_on = .false.

contains

    !> Abort on same-level seam topologies no halo reconciles (silent conservation leaks otherwise). Run whenever the block set
    !! changes (regrid, restart) on the replicated region metadata: each rank tests its own blocks against every block, so the pairs
    !! are covered once per orientation across the machine (every block has an owner) and any hit aborts everyone; O(owned x
    !! nblocks) per rank rather than an O(nblocks^2) all-pairs scan. Two cases: (a) adjacency without the exact transverse match
    !! f_amr_seam requires, reachable only via IB body-bbox expansion (clustering merges any too-close pair; tiling emits a regular
    !! grid), which the fine-fine halo can never pair; (b) same-level box intersection, reachable only via child IB body-bbox
    !! expansion, which unlike the L1 path has no overlap-merge pass, double-restricting/refluxing the shared cells.
    impure subroutine s_amr_check_seam_topology()

        integer :: ix, xb, yb, d, t
        logical :: adj, tover

        call s_amr_refresh_my_blocks()
        do ix = 1, amr_n_my
            xb = amr_my_blk(ix)
            do yb = 1, amr_num_blocks
                if (xb == yb) cycle
                if (amr_block_level(xb) /= amr_block_level(yb)) cycle
                ! same-level intersection (different levels legitimately nest; tiling emits disjoint tiles, the L1 IB pass merges
                ! overlapping boxes, but the child IB body-bbox expansion has no overlap-merge pass)
                if (f_amr_boxes_overlap(amr_region_lo_all(:,xb), amr_region_hi_all(:,xb), amr_region_lo_all(:,yb), &
                    & amr_region_hi_all(:,yb))) then
                    call s_mpi_abort("AMR: two same-level blocks INTERSECT (the child IB body-bbox expansion route can " &
                                     & // "produce this - it has no overlap-merge pass): the overlapping cells would be " &
                                     & // "restricted and refluxed twice, silently breaking conservation. Adjust the body/" &
                                     & // "regrid inputs so body-expanded child boxes merge or separate.")
                end if
                do d = 1, num_dims
                    ! relaxed adjacency: touching faces in dim d with any transverse overlap
                    adj = amr_region_lo_all(d, yb) == amr_region_hi_all(d, xb) + 1
                    if (.not. adj) cycle
                    tover = .true.
                    do t = 1, num_dims
                        if (t /= d) tover = tover .and. amr_region_lo_all(t, xb) <= amr_region_hi_all(t, &
                            & yb) .and. amr_region_lo_all(t, yb) <= amr_region_hi_all(t, xb)
                    end do
                    if (tover .and. f_amr_seam_dim(xb, yb) == 0) then
                        call s_mpi_abort("AMR: two same-level blocks touch with PARTIAL transverse overlap (IB body-bbox " &
                                         & // "expansion can produce this): the fine-fine seam halo only reconciles exact-match " &
                                         & // "faces, so the shared-face flux would silently leak. Adjust amr_block/regrid inputs " // "so body-expanded boxes merge or separate.")
                    end if
                end do
            end do
        end do

    end subroutine s_amr_check_seam_topology

    !> Abort if any box exceeds the slot cap for its level. The slot coord/field arrays are allocated once to
    !! amr_ref_ratio*amr_maxc_fit fine cells, and a level-lev block spans amr_ref_ratio**lev fine cells per coarse cell, so its
    !! coarse extent must be <= amr_maxc_fit/amr_ref_ratio**(lev-1). Every emitter enforces that via s_amr_tile_box; this checks the
    !! invariant once, where the box set is final, rather than trusting each emitter to have done it.
    !!
    !! A violation is otherwise silent and catastrophic: s_amr_build_block_coords sizes the fine coords from the block's true
    !! extent, so an over-cap box writes past x_cb, corrupting the heap on every regrid and surfacing much later as "corrupted
    !! size vs. prev_size" inside an unrelated free(), nowhere near the bug.
    !! Matches s_amr_tile_box's own floor (tc = max(tc, 1)) so a collapsed dim, whose cap divides to 0, is not flagged.
    impure subroutine s_amr_check_box_caps(boxes, nboxes, box_level)

        type(t_box), intent(in) :: boxes(:)
        integer, intent(in)     :: nboxes, box_level(:)
        integer                 :: k, d, lev, cap, span

        do k = 1, nboxes
            lev = box_level(k)
            if (lev < 1) cycle  ! level-0 tiles are sized by the tile decomposition, not this cap
            do d = 1, num_dims
                cap = max(amr_maxc_fit(d)/amr_ref_ratio**(lev - 1), 1)
                span = boxes(k)%hi(d) - boxes(k)%lo(d) + 1
                if (span > cap) then
                    if (proc_rank == 0) print '(A,I0,A,I0,A,I0,A,I0)', ' [amr] box cap violated: level ', lev, ' dim ', d, &
                        & ' span ', span, ' > cap ', cap
                    call s_mpi_abort("AMR regrid: a fine box exceeds the slot cap for its level (span and cap printed above). " &
                                     & // "The fine coord arrays are sized to amr_ref_ratio*amr_maxc_fit, so this would write " &
                                     & // "past x_cb and corrupt the heap. A box emitter did not route through s_amr_tile_box.")
                end if
            end do
        end do

    end subroutine s_amr_check_box_caps

    !> Invariant check: same-level boxes are pairwise disjoint. Guaranteed by the cluster partition + merge threshold + IB
    !! overlap-merge; relied on by the rebuild's overlap carry-forward and by per-peer unpack reordering in the exchange plans, so
    !! it is enforced here rather than assumed. All levels share the global coarse index space (see the cap formula in
    !! s_amr_check_box_caps), so the interval test is valid across parents. O(nboxes^2) host integer compares per regrid; a sorted
    !! sweep would replace it if box counts grow.
    impure subroutine s_amr_check_box_disjoint(boxes, nboxes, box_level)

        type(t_box), intent(in) :: boxes(:)
        integer, intent(in)     :: nboxes, box_level(:)
        integer                 :: k, kk

        do k = 1, nboxes
            if (box_level(k) < 1) cycle
            do kk = k + 1, nboxes
                if (box_level(kk) /= box_level(k)) cycle
                if (all(boxes(k)%lo <= boxes(kk)%hi .and. boxes(kk)%lo <= boxes(k)%hi)) then
                    if (proc_rank == 0) print '(A,I0,A,I0,A,I0)', ' [amr] same-level box overlap: level ', box_level(k), &
                        & ' boxes ', k, ' and ', kk
                    call s_mpi_abort("AMR regrid: two same-level boxes overlap (indices printed above). The overlap " &
                                     & // "carry-forward and the exchange plans both assume same-level disjointness.")
                end if
            end do
        end do

    end subroutine s_amr_check_box_disjoint

    !> True iff global level-0 cell (gi, gj, gk) lies inside any acoustic source support bbox.
    pure logical function f_in_acoustic_support(gi, gj, gk) result(insup)

        integer, intent(in) :: gi, gj, gk
        integer             :: s

        insup = .false.
        do s = 1, num_source
            if (gi >= acoustic_supp_lo(1, s) .and. gi <= acoustic_supp_hi(1, s) .and. (n_glb == 0 .or. (gj >= acoustic_supp_lo(2, &
                & s) .and. gj <= acoustic_supp_hi(2, s))) .and. (p_glb == 0 .or. (gk >= acoustic_supp_lo(3, &
                & s) .and. gk <= acoustic_supp_hi(3, s)))) then
                insup = .true.; return
            end if
        end do

    end function f_in_acoustic_support

    !> Clip a candidate regrid box (global indices) clear of every acoustic source support bbox: per overlapping source, remove the
    !! overlap along the single axis/side keeping the largest remaining extent (deterministic: lower axis, then begin side, wins
    !! ties). Only shrinks; may empty the box (hi < lo); the caller drops empties.
    impure subroutine s_amr_clip_box_from_sources(lo, hi)

        integer, intent(inout) :: lo(3), hi(3)
        integer                :: s, d, best_d, best_side, best_ext, ext_l, ext_r

        do s = 1, num_source
            if (hi(1) < lo(1) .or. hi(2) < lo(2) .or. hi(3) < lo(3)) return  ! emptied by an earlier clip
            if (.not. f_amr_boxes_overlap(lo, hi, acoustic_supp_lo(:,s), acoustic_supp_hi(:,s))) cycle
            best_d = 1; best_side = 1; best_ext = -1
            do d = 1, num_dims
                ext_l = acoustic_supp_lo(d, s) - lo(d)  ! cells kept by [lo(d), supp_lo-1]
                ext_r = hi(d) - acoustic_supp_hi(d, s)  ! cells kept by [supp_hi+1, hi(d)]
                if (ext_l > best_ext) then; best_ext = ext_l; best_d = d; best_side = 1; end if
                if (ext_r > best_ext) then; best_ext = ext_r; best_d = d; best_side = 2; end if
            end do
            if (best_side == 1) then
                hi(best_d) = acoustic_supp_lo(best_d, s) - 1
            else
                lo(best_d) = acoustic_supp_hi(best_d, s) + 1
            end if
        end do
        ! safety net: clipping removed every overlap by construction; anything left is a bug
        do s = 1, num_source
            if (hi(1) < lo(1) .or. hi(2) < lo(2) .or. hi(3) < lo(3)) return
            if (f_amr_boxes_overlap(lo, hi, acoustic_supp_lo(:,s), acoustic_supp_hi(:, &
                & s))) call s_mpi_abort('amr regrid: acoustic source exclusion clip failed (internal error)')
        end do

    end subroutine s_amr_clip_box_from_sources

    !> Recompute the global Lagrangian-cloud exclusion bbox (collective: allreduces the rank-local position extrema). pad_cells
    !! covers smearing + stencil (+ drift until the next recompute). No-op (lag_supp_on = false) when no rank holds a bubble.
    impure subroutine s_amr_compute_lag_supp(pad_cells)

        integer, intent(in)    :: pad_cells
        real(wp), dimension(3) :: pmin_loc, pmax_loc, pmin_glb, pmax_glb
        integer                :: d

        call s_lag_cloud_bbox_local(pmin_loc, pmax_loc)
        do d = 1, 3
            call s_mpi_allreduce_min(pmin_loc(d), pmin_glb(d))
            call s_mpi_allreduce_max(pmax_loc(d), pmax_glb(d))
        end do
        lag_supp_on = pmin_glb(1) <= pmax_glb(1)
        if (.not. lag_supp_on) return
        call s_lag_phys_to_cells(pmin_glb, pmax_glb, pad_cells, lag_supp_lo, lag_supp_hi)

    end subroutine s_amr_compute_lag_supp

    !> True iff global level-0 cell (gi, gj, gk) lies inside the Lagrangian-cloud exclusion bbox.
    pure logical function f_in_lag_support(gi, gj, gk) result(insup)

        integer, intent(in) :: gi, gj, gk

        insup = .false.
        if (.not. lag_supp_on) return
        insup = gi >= lag_supp_lo(1) .and. gi <= lag_supp_hi(1) .and. (n_glb == 0 .or. (gj >= lag_supp_lo(2) &
                                  & .and. gj <= lag_supp_hi(2))) .and. (p_glb == 0 .or. (gk >= lag_supp_lo(3) &
                                  & .and. gk <= lag_supp_hi(3)))

    end function f_in_lag_support

    !> Clip a candidate regrid box (global indices) clear of one support bbox: remove the overlap along the single axis/side that
    !! keeps the largest remaining extent (deterministic: lower axis, then begin side, wins ties). Only shrinks; may empty the box
    !! (hi < lo).
    pure subroutine s_amr_clip_box_from_supp(lo, hi, slo, shi)

        integer, intent(inout) :: lo(3), hi(3)
        integer, intent(in)    :: slo(3), shi(3)
        integer                :: d, best_d, best_side, best_ext, ext_l, ext_r

        if (hi(1) < lo(1) .or. hi(2) < lo(2) .or. hi(3) < lo(3)) return
        if (.not. f_amr_boxes_overlap(lo, hi, slo, shi)) return
        best_d = 1; best_side = 1; best_ext = -1
        do d = 1, num_dims
            ext_l = slo(d) - lo(d)
            ext_r = hi(d) - shi(d)
            if (ext_l > best_ext) then; best_ext = ext_l; best_d = d; best_side = 1; end if
            if (ext_r > best_ext) then; best_ext = ext_r; best_d = d; best_side = 2; end if
        end do
        if (best_side == 1) then
            hi(best_d) = slo(best_d) - 1
        else
            lo(best_d) = shi(best_d) + 1
        end if

    end subroutine s_amr_clip_box_from_supp

    !> active_box + AMR containment: every active block must sit strictly inside the active window (one-cell margin). Two reasons:
    !! the windowed coarse RK update would silently drop a reflux correction at a face cell outside the window (conservation leak),
    !! and the coarse RHS only computes fluxes inside it. The window only grows (s_grow_active_box monotone, self-disabling at full
    !! domain), so containment set at init and re-established each regrid holds between. Collective (same window/block metadata on
    !! all ranks).
    impure subroutine s_amr_check_active_box_containment()

        integer :: k
        logical :: ok

        ! ab_active is only true at num_procs == 1 (m_active_box disables itself under MPI), so ab and block indices share
        ! the same (global == local) index space

        if ((.not. amr) .or. (.not. ab_active)) return
        do k = 1, amr_num_blocks
            ! L0 tiles span the base grid by construction; the containment rule is for fine blocks
            if (amr_block_level(k) == 0) cycle
            ok = all((amr_region_lo_all(:,k) > [ab_x%beg, ab_y%beg, ab_z%beg] .and. amr_region_hi_all(:,k) < [ab_x%end, ab_y%end, &
                     & ab_z%end]) .or. .not. amr_dim)
            if (.not. ok) then
                call s_mpi_abort('amr with active_box: an AMR block is not strictly inside the active ' &
                                 & // 'window; place the initial block (with a one-cell margin) inside the ' &
                                 & // 'initial non-ambient region plus buff_size')
            end if
        end do

    end subroutine s_amr_check_active_box_containment

    !> This rank's own tagged cells as global level-0 coordinates. Each rank scans only its interior (0:m, 0:n, 0:p), which the
    !! level-0 decomposition makes disjoint, so no global cell is emitted twice and a SUM reduction over these lists counts every
    !! tagged cell exactly once. No rank holds the global tag list: per-rank memory and wire volume scale with the rank's own tag
    !! count.
    impure subroutine s_amr_local_tags(tag_grid, sidx, tags, ntag)

        logical, intent(in)               :: tag_grid(0:,0:,0:)
        integer, intent(in)               :: sidx(3)
        integer, allocatable, intent(out) :: tags(:,:)
        integer, intent(out)              :: ntag
        integer                           :: ci, cj, ck

        ntag = 0
        do ck = 0, p; do cj = 0, n; do ci = 0, m
            if (tag_grid(ci, cj, ck)) ntag = ntag + 1
        end do; end do; end do
        allocate (tags(3, max(ntag, 1)))
        ntag = 0
        do ck = 0, p; do cj = 0, n; do ci = 0, m
            if (tag_grid(ci, cj, ck)) then
                ntag = ntag + 1
                tags(1, ntag) = ci + sidx(1)
                tags(2, ntag) = 0
                tags(3, ntag) = 0
                if (n_glb > 0) tags(2, ntag) = cj + sidx(2)
                if (p_glb > 0) tags(3, ntag) = ck + sidx(3)
            end if
        end do; end do; end do

    end subroutine s_amr_local_tags

    !> Grow the per-level pack buffers sidx(:) (int8 linear index) / skb(:) (parent box id) geometrically so at least nloc+extra
    !! slots fit; preserves the first nloc entries. Amortized O(1) append for s_amr_pack_gwin_pairs.
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

    !> Pack this rank's owned tagged cells of the child window [mlo:mhi] as (linear-index, kb) pairs, appended to the per-level send
    !! arrays sidx(:) (int8 linear index) / skb(:) (parent box id). The int8 encode matches the pass-2 decode, so gathering these
    !! pairs across ranks and setting them into a per-parent dense window dedups replicated/overlapping tags and fixes the (k,j,i)
    !! extraction order, giving identical child boxes on every rank. One allgatherv per level (caller) keeps the collective count
    !! O(#levels) rather than O(#parent-boxes). gwin is read, not modified.
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

    !> Regrid: tag by relative density gradient, cluster (Berger-Rigoutsos + min-separation merge) into separated boxes, pad/clamp/
    !! size-cap each, rebuild every active slot. Each new slot prolongs from coarse then overwrites its overlap with whichever old
    !! slot(s) covered it (rank-local; a split copies from one old slot, a merge from both). Called between steps only. No-op if
    !! nothing is tagged or the box set is unchanged.
    impure subroutine s_amr_regrid(q_cons_base)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_base
        logical, allocatable                                   :: tag_grid(:,:,:)
        type(t_box), allocatable                               :: boxes(:)
        integer                                                :: sidx(3), nboxes
        integer                                                :: old_np
        ! heap, not stack: these are O(global boxes) and overflow a default stack at large box
        ! counts. Unsaved local allocatables are deallocated automatically on every return path.
        integer, allocatable :: box_level(:)
        integer, allocatable :: old_ilo(:,:), old_ext(:,:)
        integer, allocatable :: old_level(:)
        logical, allocatable :: old_owns(:)
        logical              :: same
        integer              :: i

#ifdef MFC_MPI
#endif

        allocate (box_level(amr_max_fine), old_ilo(3, amr_max_blocks), old_ext(3, amr_max_blocks), old_level(amr_max_blocks), &
                  & old_owns(amr_max_blocks))

        ! valid coarse cons ghosts at internal rank boundaries: the tag sweep reads +/-1 across seams and the rebuild prolongation
        ! reads past the new intersection (all ranks call: pairwise per-direction exchange; complete no-op at np=1).

        call s_amr_exchange_coarse_cons_halo(q_cons_base)
        do i = 1, sys_size
            $:GPU_UPDATE(host='[q_cons_base(i)%sf]')
        end do

        ! Lagrangian-cloud exclusion bbox for this regrid (collective): smearing (mapCells) + stencil headroom (2) + drift
        ! margin until the next regrid (amr_buf)
        if (bubbles_lagrange) call s_amr_compute_lag_supp(mapCells + 2 + amr_buf)

        call s_amr_regrid_tag_cells(q_cons_base, tag_grid, sidx)
        call s_amr_cad_count(tag_grid, sidx)  ! [amr-cad] cadence containment audit (counts only; report at finalize)
        call s_amr_regrid_cluster_tags(tag_grid, sidx, boxes, nboxes)
        if (nboxes == 0) return  ! nothing tagged on any rank; keep the current blocks
        call s_amr_regrid_shape_boxes(boxes, nboxes)
        if (nboxes == 0) return  ! every box was confined to the domain margin
        call s_amr_regrid_nest_children(boxes, nboxes, box_level)
        if (amr_snap > 0) call s_amr_regrid_snap_boxes(boxes, nboxes, box_level)
        call s_amr_check_box_caps(boxes, nboxes, box_level)  ! invariant: no box may exceed its level's slot cap
        call s_amr_check_box_disjoint(boxes, nboxes, box_level)  ! invariant: same-level boxes are pairwise disjoint
        call s_amr_regrid_boxes_unchanged(boxes, nboxes, box_level, same)
        if (same) return  ! identical box set and levels: keep the live slots
        call s_amr_regrid_stash_migrate(boxes, nboxes, box_level, old_np, old_ilo, old_ext, old_level, old_owns)
        call s_amr_regrid_rebuild_slots(q_cons_base, boxes, nboxes, old_np, old_ilo, old_ext, old_level, old_owns)

    end subroutine s_amr_regrid

    !> Regrid phase 1: per-cell tag field (density-gradient criterion), skipping the two global boundary cells per active dim and
    !! suppressing tags over the acoustic source supports and the Lagrangian-cloud exclusion bbox. sidx returns the rank's global
    !! start offsets for the sparse union in phase 2.
    impure subroutine s_amr_regrid_tag_cells(q_cons_base, tag_grid, sidx)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_base
        logical, allocatable, intent(inout)                 :: tag_grid(:,:,:)
        integer, intent(out)                                :: sidx(3)
        integer                                             :: tg_lo(3), tg_hi(3), ci, cj, ck
        real(wp)                                            :: r0, g

        ! 1) per-cell tag field (density-gradient criterion), skipping the two global boundary cells per active dim

        sidx = amr_sidx
        tg_lo = merge(merge(1, 0, sidx == 0), 0, amr_dim)
        tg_hi = merge(merge(amr_ext - 1, amr_ext, sidx + amr_ext == [m_glb, n_glb, p_glb]), 0, amr_dim)
        allocate (tag_grid(0:m,0:n,0:p)); tag_grid = .false.
        do ck = tg_lo(3), tg_hi(3)
            do cj = tg_lo(2), tg_hi(2)
                do ci = tg_lo(1), tg_hi(1)
                    ! total density gradient (sum of the continuity variables): degenerates to the single-fluid tagger, immune to
                    ! trace-fluid noise. Matched-density composition-only interfaces are invisible (documented limit).
                    r0 = max(abs(f_amr_rho_tot_sf(q_cons_base, ci, cj, ck)), 1.e-30_wp)
                    g = abs(f_amr_rho_tot_sf(q_cons_base, ci + 1, cj, ck) - f_amr_rho_tot_sf(q_cons_base, ci - 1, cj, ck))
                    if (n_glb > 0) g = max(g, abs(f_amr_rho_tot_sf(q_cons_base, ci, cj + 1, ck) - f_amr_rho_tot_sf(q_cons_base, &
                        & ci, cj - 1, ck)))
                    if (p_glb > 0) g = max(g, abs(f_amr_rho_tot_sf(q_cons_base, ci, cj, ck + 1) - f_amr_rho_tot_sf(q_cons_base, &
                        & ci, cj, ck - 1)))
                    ! 2*r0 normalizes the 2-cell central difference (rho at i+1..i-1); the 2 is the stencil span, not the
                    ! refinement ratio
                    if (g/(2._wp*r0) > amr_tag_eps) tag_grid(ci, cj, ck) = .true.
                    ! the acoustic source support stays coarse (its spatials are coarse cell indices): suppress tags there so
                    ! the clusterer splits around the source
                    if (acoustic_source .and. tag_grid(ci, cj, ck)) then
                        if (f_in_acoustic_support(ci + sidx(1), cj + sidx(2), ck + sidx(3))) tag_grid(ci, cj, ck) = .false.
                    end if
                    ! the Lagrangian bubble cloud stays coarse (two-way coupling lives on the coarse grid): suppress tags over its
                    ! padded bbox
                    if (bubbles_lagrange .and. tag_grid(ci, cj, ck)) then
                        if (f_in_lag_support(ci + sidx(1), cj + sidx(2), ck + sidx(3))) tag_grid(ci, cj, ck) = .false.
                    end if
                end do
            end do
        end do

    end subroutine s_amr_regrid_tag_cells

    !> [amr-cad] cadence containment audit: count this rank's level-1 tags, and how many fall outside the pre-regrid level-1
    !! coverage, i.e. a feature that evolved unrefined since the last regrid because amr_buf did not cover its drift over
    !! amr_regrid_int steps. Counts only (reported once by s_amr_cad_report); the first refinement from scratch is skipped (every
    !! tag is new by construction). Region bounds are global coarse cells; tag_grid is rank-local, so paint each region's clip with
    !! this subdomain into a local mask first.
    impure subroutine s_amr_cad_count(tag_grid, sidx)

        logical, intent(in)  :: tag_grid(0:,0:,0:)
        integer, intent(in)  :: sidx(3)
        logical, allocatable :: cov(:,:,:)
        integer              :: k, ci, cj, ck, bl(3), bh(3)
        logical              :: any_l1

        ! skip the first regrid: it populates the hierarchy from the seed block, so nearly every tag is
        ! legitimately outside the old coverage. The instrument measures steady-state containment, regrid 2 onward.

        if (.not. amr_cad_armed) then
            amr_cad_armed = .true.
            return
        end if
        any_l1 = .false.
        do k = 1, amr_num_blocks
            if (amr_block_level(k) == 1) any_l1 = .true.
        end do
        if (.not. any_l1) return
        allocate (cov(0:m,0:n,0:p)); cov = .false.
        do k = 1, amr_num_blocks
            if (amr_block_level(k) /= 1) cycle
            bl = merge(max(amr_region_lo_all(:,k) - sidx, 0), 0, amr_dim)
            bh = merge(min(amr_region_hi_all(:,k) - sidx, amr_ext), 0, amr_dim)
            if (bl(1) > bh(1) .or. bl(2) > bh(2) .or. bl(3) > bh(3)) cycle
            cov(bl(1):bh(1),bl(2):bh(2),bl(3):bh(3)) = .true.
        end do
        do ck = 0, p
            do cj = 0, n
                do ci = 0, m
                    if (tag_grid(ci, cj, ck)) then
                        amr_cad_tot = amr_cad_tot + 1
                        if (.not. cov(ci, cj, ck)) amr_cad_esc = amr_cad_esc + 1
                    end if
                end do
            end do
        end do
        deallocate (cov)

    end subroutine s_amr_cad_count

    !> Regrid phase 2: build this rank's own sparse tag list, then cluster it with per-node signature reductions into a list of
    !! separated candidate boxes (nboxes = 0 if nothing is tagged on any rank).
    impure subroutine s_amr_regrid_cluster_tags(tag_grid, sidx, boxes, nboxes)

        logical, allocatable, intent(inout)   :: tag_grid(:,:,:)
        integer, intent(in)                   :: sidx(3)
        type(t_box), allocatable, intent(out) :: boxes(:)
        integer, intent(out)                  :: nboxes
        integer, allocatable                  :: tags(:,:)
        integer                               :: ntag

        ! 2) build this rank's local tag list, then cluster it; the tree is driven by reduced signatures, so it is rank-invariant

        call s_amr_local_tags(tag_grid, sidx, tags, ntag)
        deallocate (tag_grid)
        call s_amr_cluster(tags, ntag, boxes, nboxes, .true.)
        deallocate (tags)

    end subroutine s_amr_regrid_cluster_tags

    !> Regrid hysteresis (amr_snap > 0): every new box within amr_snap coarse cells per face of a live block of the same level takes
    !! that block's box. A feature drifting a cell between regrids otherwise shifts every tile of its envelope by that cell and
    !! re-creates every block; snapped boxes are identical, and when every box snaps s_amr_regrid_boxes_unchanged skips the rebuild
    !! outright. Coverage: a new box is the tags padded by amr_buf, so a snap of <= amr_snap <= amr_buf - 2 cells (the validator's
    !! bound) keeps >= 2 cells of padding on every face; the cadence audit ([amr-cad] escaped) is the runtime check. All-or-none:
    !! the snapped set must stay pairwise disjoint per level and every level >= 2 box must stay inside a single parent box by
    !! amr_cpat_mar (the nester's window), else the whole snap is dropped and the fresh boxes stand. Replicated inputs, so every
    !! rank decides alike.
    impure subroutine s_amr_regrid_snap_boxes(boxes, nboxes, box_level)

        type(t_box), intent(inout) :: boxes(:)
        integer, intent(in)        :: nboxes, box_level(:)
        type(t_box), allocatable   :: snapped(:)
        integer                    :: k, kk, ks, npar, nsnap, mlo(3), mhi(3)
        logical                    :: ok

        allocate (snapped(nboxes)); snapped(1:nboxes) = boxes(1:nboxes)
        nsnap = 0
        do k = 1, nboxes
            do ks = l0_slot_off + 1, amr_num_blocks
                if (amr_block_level(ks) /= box_level(k)) cycle
                if (all(abs(amr_region_lo_all(:,ks) - boxes(k)%lo) <= amr_snap) .and. all(abs(amr_region_hi_all(:, &
                    & ks) - boxes(k)%hi) <= amr_snap)) then
                    if (any(amr_region_lo_all(:,ks) /= boxes(k)%lo) .or. any(amr_region_hi_all(:, &
                        & ks) /= boxes(k)%hi)) nsnap = nsnap + 1
                    snapped(k)%lo = amr_region_lo_all(:,ks); snapped(k)%hi = amr_region_hi_all(:,ks)
                    exit
                end if
            end do
        end do
        ok = nsnap > 0
        ! same-level disjointness of the snapped set
        do k = 1, nboxes
            if (.not. ok) exit
            do kk = k + 1, nboxes
                if (box_level(kk) /= box_level(k)) cycle
                if (all(snapped(k)%lo <= snapped(kk)%hi .and. snapped(kk)%lo <= snapped(k)%hi)) then
                    ok = .false.; exit
                end if
            end do
        end do
        ! proper nesting: a level >= 2 box lies inside exactly one parent-level box, inset by the nesting margin
        do k = 1, nboxes
            if (.not. ok) exit
            if (box_level(k) < 2) cycle
            npar = 0
            do kk = 1, nboxes
                if (box_level(kk) /= box_level(k) - 1) cycle
                if (.not. all(snapped(k)%lo <= snapped(kk)%hi .and. snapped(kk)%lo <= snapped(k)%hi)) cycle
                npar = npar + 1
                call s_amr_nest_window(snapped(kk), mlo, mhi)
                if (any(snapped(k)%lo < mlo) .or. any(snapped(k)%hi > mhi)) ok = .false.
            end do
            if (npar /= 1) ok = .false.
        end do
        if (ok) boxes(1:nboxes) = snapped(1:nboxes)
        if (rank_time_wrt .and. proc_rank == 0) write (0, '(A,I0,A,I0,A,L1)') '[amr-snap] boxes ', nboxes, ' snapped ', nsnap, &
            & ' applied ', ok
        deallocate (snapped)

    end subroutine s_amr_regrid_snap_boxes

    !> Clip a box to strictly inside the active window (one-cell margin) in every active dimension.
    pure subroutine s_amr_clip_to_active_box(lo, hi)

        integer, intent(inout) :: lo(3), hi(3)

        lo = max(lo, merge([ab_x%beg, ab_y%beg, ab_z%beg] + 1, lo, amr_dim))
        hi = min(hi, merge([ab_x%end, ab_y%end, ab_z%end] - 1, hi, amr_dim))

    end subroutine s_amr_clip_to_active_box

    !> Regrid phase 3: pad + clamp + size-cap each box, clip it clear of the acoustic/Lagrangian supports and the active window,
    !! expand it over immersed bodies, then tile oversized boxes (non-IB) or merge overlapping ones (IB).
    impure subroutine s_amr_regrid_shape_boxes(boxes, nboxes)

        type(t_box), allocatable, intent(inout) :: boxes(:)
        integer, intent(inout)                  :: nboxes
        integer                                 :: lo(3), hi(3), k, kk
        logical                                 :: merged

        ! 3) pad + clamp + size-cap each box (amr_maxc_fit lets each box move freely across rank boundaries); drop margin-only boxes

        k = 0
        do kk = 1, nboxes
            lo = boxes(kk)%lo; hi = boxes(kk)%hi
            lo(1) = max(lo(1) - amr_buf, buff_size); hi(1) = min(hi(1) + amr_buf, m_glb - buff_size)
            ! IB keeps the size-cap clamp (a body needs one contiguous block; splitting a body across tiles is untested); the
            ! general path leaves boxes full-size and tiles them (below) into <= amr_maxc_fit sub-blocks with a fine-fine halo
            if (ib .and. hi(1) - lo(1) + 1 > amr_maxc_fit(1)) hi(1) = lo(1) + amr_maxc_fit(1) - 1
            if (n_glb > 0) then
                lo(2) = max(lo(2) - amr_buf, buff_size); hi(2) = min(hi(2) + amr_buf, n_glb - buff_size)
                if (ib .and. hi(2) - lo(2) + 1 > amr_maxc_fit(2)) hi(2) = lo(2) + amr_maxc_fit(2) - 1
            else
                lo(2) = 0; hi(2) = 0
            end if
            if (p_glb > 0) then
                lo(3) = max(lo(3) - amr_buf, buff_size); hi(3) = min(hi(3) + amr_buf, p_glb - buff_size)
                if (ib .and. hi(3) - lo(3) + 1 > amr_maxc_fit(3)) hi(3) = lo(3) + amr_maxc_fit(3) - 1
            else
                lo(3) = 0; hi(3) = 0
            end if
            ! keep candidate boxes clear of every acoustic source support (the source acts on the coarse grid only); clipping
            ! only shrinks, so boxes stay disjoint; empties drop below
            if (acoustic_source) call s_amr_clip_box_from_sources(lo, hi)
            if (bubbles_lagrange .and. lag_supp_on) call s_amr_clip_box_from_supp(lo, hi, lag_supp_lo, lag_supp_hi)
            ! active_box: boxes stay strictly inside the active window (the windowed coarse update would drop reflux corrections
            ! at faces outside it). Tags cannot arise outside (frozen-ambient exterior), so only the amr_buf padding is ever cut,
            ! and the cut cells are ambient. np=1 only (ab_active is false under MPI).
            if (ab_active) call s_amr_clip_to_active_box(lo, hi)
            ! a fine block that partially covers an immersed body is an untested regime (ghost prolongation through body-interior
            ! cells, refluxing across the body): any box overlapping a body's bounding box expands to contain the whole body plus
            ! margin
            if (ib) call s_amr_expand_box_over_bodies(lo, hi)
            if (hi(1) < lo(1) .or. hi(2) < lo(2) .or. hi(3) < lo(3)) cycle  ! confined to the domain margin
            k = k + 1; boxes(k)%lo = lo; boxes(k)%hi = hi
        end do
        nboxes = k
        if (nboxes == 0) return

        ! max_grid_size tiling (non-IB): split any box larger than amr_maxc_fit into contiguous <= amr_maxc_fit sub-blocks so a
        ! whole block fits a rank's local solver scratch. Tiles are adjacent; the block-to-block fine-fine halo
        ! (s_amr_fine_fine_halo) makes the seams conservative and the reflux skips fine-fine faces. (IB keeps the clamp, above.)
        if (.not. ib) then
            block
                type(t_box), allocatable :: tiled(:)
                integer                  :: kk2, ntl, capt
                allocate (tiled(amr_max_blocks))
                ntl = 0; capt = 0
                do kk2 = 1, nboxes
                    call s_amr_tile_box(boxes(kk2)%lo, boxes(kk2)%hi, tiled, ntl, amr_max_fine, capt)
                end do
                if (capt == 1 .and. proc_rank == 0) print '(A,I0)', ' [amr] WARNING: tiling capped at amr_max_blocks = ', &
                    & amr_max_blocks
                deallocate (boxes); call move_alloc(tiled, boxes)
                nboxes = ntl
            end block
        end if

        if (ib) then
            ! body-containment expansion can make boxes overlap (bisection guarantees disjoint boxes, but two boxes near one body
            ! both grow over it): merge pairs closer than a 2-cell gap to a bbox until none remain. Overlapping blocks would
            ! double-restrict/reflux; a 1-cell gap with transverse overlap gives the two blocks a coincident outside coarse
            ! cell, which the batched reflux apply writes from both blocks in one kernel, an unsynchronized read-modify-write
            ! (the clusterer's min-separation merge guarantees a >= 2 gap everywhere else; this restores it after expansion).
            merged = .true.
            do while (merged)
                merged = .false.
                outer: do k = 1, nboxes - 1
                    do kk = k + 1, nboxes
                        if (boxes(k)%lo(1) <= boxes(kk)%hi(1) + 1 .and. boxes(k)%hi(1) >= boxes(kk)%lo(1) - 1 .and. (n_glb == 0 &
                            & .or. (boxes(k)%lo(2) <= boxes(kk)%hi(2) + 1 .and. boxes(k)%hi(2) >= boxes(kk)%lo(2) - 1)) &
                            & .and. (p_glb == 0 .or. (boxes(k)%lo(3) <= boxes(kk)%hi(3) + 1 .and. boxes(k)%hi(3) &
                            & >= boxes(kk)%lo(3) - 1))) then
                            boxes(k)%lo = min(boxes(k)%lo, boxes(kk)%lo)
                            boxes(k)%hi = max(boxes(k)%hi, boxes(kk)%hi)
                            boxes(kk) = boxes(nboxes); nboxes = nboxes - 1
                            if (boxes(k)%hi(1) - boxes(k)%lo(1) + 1 > amr_maxc_fit(1) .or. (n_glb > 0 .and. boxes(k)%hi(2) &
                                & - boxes(k)%lo(2) + 1 > amr_maxc_fit(2)) .or. (p_glb > 0 .and. boxes(k)%hi(3) - boxes(k)%lo(3) &
                                & + 1 > amr_maxc_fit(3))) then
                                call s_mpi_abort('amr regrid: merging body-containing blocks exceeds ' &
                                                 & // 'the per-rank block size cap')
                            end if
                            merged = .true.
                            exit outer
                        end if
                    end do
                end do outer
            end do
            ! the expansion may also have grown a box onto an acoustic source support or the Lagrangian cloud: the constraints
            ! (contain the body, exclude the source/cloud) cannot both hold, so fail closed
            if (acoustic_source .or. (bubbles_lagrange .and. lag_supp_on)) then
                do k = 1, nboxes
                    lo = boxes(k)%lo; hi = boxes(k)%hi
                    if (acoustic_source) call s_amr_clip_box_from_sources(lo, hi)
                    if (bubbles_lagrange .and. lag_supp_on) call s_amr_clip_box_from_supp(lo, hi, lag_supp_lo, lag_supp_hi)
                    if (ab_active) call s_amr_clip_to_active_box(lo, hi)
                    if (any(lo /= boxes(k)%lo) .or. any(hi /= boxes(k)%hi)) then
                        call s_mpi_abort('amr regrid: a block must contain an immersed body AND stay ' &
                                         & // 'clear of an acoustic source support / Lagrangian bubble cloud - the ' &
                                         & // 'constraints conflict; move the body, source, or cloud apart')
                    end if
                end do
            end if
        end if

    end subroutine s_amr_regrid_shape_boxes

    !> Regrid phase 3b: multi-level nesting. Hierarchically append level-l child boxes inside each level-(l-1) box for l =
    !! 2..amr_max_level, parents first (the build loop fills a parent before its child's gather-from-parent reads it), and set
    !! box_level for every box (1 for the L0->L1 boxes). Sensor-on-fine: a child's extent is the density-gradient sensor run on the
    !! parent-level fine solution (the still-live old level-(l-1) blocks, read before the stash), coarsened to L0 cells and
    !! clustered, so children track features inside the parent. A brand-new region with no old fine data gets a centred inset (the
    !! sensor takes over next regrid); a parent with a smooth fine solution gets no child.
    !!
    !! Per level: collect -> one exchange -> process. Pass 1 tags each parent's nesting window from this rank's owned old
    !! blocks and packs the tags as (linear L0 index, parent kb) pairs; one alltoall routes every parent's pairs to its
    !! clustering owner (round-robin over kb: balanced, and a pure function of kb so every rank agrees without talking); pass 2
    !! rebuilds each owned parent's dense window from the pairs (dedup and (k,j,i) order independent of arrival) and clusters
    !! it into children, emitted without the global cap because a rank cannot see the global count mid-pass; one allgatherv
    !! plus a stable sort by kb replays the children into `boxes` in (kb, emission) order, the order a serial loop over parents
    !! would produce, so the box list, its truncation at amr_max_fine and box_level are rank-invariant. With IB, nesting is np=1
    !! only (m_checker). Regions stay in L0 cell indices.
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

        ! the nesting appends into `boxes` (up to amr_max_blocks). The non-IB path already grew it via the tiling move_alloc; the
        ! IB path (merges only) leaves it at the cluster count, so grow it here or the appends overrun the allocation
        if (size(boxes) < amr_max_blocks) then
            allocate (grown(amr_max_blocks)); grown(1:nboxes) = boxes(1:nboxes); call move_alloc(grown, boxes)
        end if
        ! host-refresh the live (old) blocks' conserved state: the fine sensor reads the flat store on the host, but the stash's
        ! GPU_UPDATE(host) runs after this nesting, so the host copy is stale here. np>1: only the owner holds the fine state.
        do ob = 1, amr_num_blocks
            if (amr_block_level(ob) == 0 .or. .not. amr_owns_all(ob)) cycle
            $:GPU_UPDATE(host='[amr_cons_st(:, :, :, :, amr_loc_of(ob))]')
        end do
        call s_amr_refresh_my_blocks()

        plo = 1; phi = nboxes  ! [plo:phi] = the boxes at level lev-1, the parents to nest inside
        do lev = 2, amr_max_level
            if (phi < plo) exit  ! nothing nested at the previous level -> no deeper levels possible
            newlo = nboxes + 1
            allocate (covered(plo:phi), mine(plo:phi), mlo_all(3,plo:phi), mhi_all(3,plo:phi), powner(plo:phi), mych(7, &
                      & amr_max_fine))
            call s_amr_nest_mark_parents(boxes, plo, phi, lev, mine, covered)
            do kb = plo, phi
                powner(kb) = mod(kb - plo, max(num_procs, 1))
                call s_amr_nest_window(boxes(kb), mlo_all(:,kb), mhi_all(:,kb))
            end do

            ! pass 1: a rank builds a parent's window only when it holds an old block overlapping it (the allocate+zero would
            ! be O(parents x window volume) on every non-contributing rank); under IB the owner always builds one for the body tags
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
                if (nmych + 1 > amr_max_fine) exit  ! local buffer full (bounded by the same global cap)
                mlo = mlo_all(:,kb); mhi = mhi_all(:,kb)
                if (f_amr_nest_window_empty(mlo, mhi)) cycle
                call s_amr_nest_window_tags(kb, mlo, mhi, gidx, gkb, ntot_g, ctags, nct)
                if (.not. covered(kb)) then
                    call s_amr_nest_inset_child(boxes(kb), lev, kb, mych, nmych)  ! brand-new region: no fine tags to cluster
                else if (nct > 0) then
                    call s_amr_nest_cluster_children(ctags, nct, mlo, mhi, lev, kb, mych, nmych)
                end if
                deallocate (ctags)
            end do

            call s_amr_nest_gather_children(mych, nmych, plo, phi, gch, chord, ntot_ch)
            do jch = 1, ntot_ch
                if (nboxes + 1 > amr_max_fine) exit  ! pool full: stop nesting (canonical order => same truncation)
                nboxes = nboxes + 1
                boxes(nboxes)%lo = gch(1:3,chord(jch)); boxes(nboxes)%hi = gch(4:6,chord(jch)); box_level(nboxes) = lev
            end do
            deallocate (gch, chord, powner, mych, gidx, gkb, covered, mine, mlo_all, mhi_all)
            plo = newlo; phi = nboxes  ! the boxes just appended are the parents for the next level
        end do
        if (nboxes >= amr_max_fine .and. proc_rank == 0) print '(A)', &
            & ' [amr] NOTE: block pool full during multi-level nesting; some boxes were not refined further'

    end subroutine s_amr_regrid_nest_children

    !> mine(kb): this rank holds a level-(lev-1) block overlapping parent kb; covered(kb): any rank does. One pass over the owned
    !! blocks (not O(parents x global blocks) per rank), and `covered`, which must stay replicated, is one LOR over the parents:
    !! every block has exactly one owner, so the union over ranks of "my blocks overlapping kb" is "all blocks overlapping kb".
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

        e = hi(1) < lo(1)
        if (n_glb > 0) e = e .or. hi(2) < lo(2)
        if (p_glb > 0) e = e .or. hi(3) < lo(3)

    end function f_amr_nest_window_empty

    !> Pass 1 for one parent: tag its window from this rank's owned level-(lev-1) blocks (amr_block_level still holds the old levels
    !! here; the rebuild resets it to box_level). IB: the body region is refined at every level even where the sensor is quiet, by
    !! marking its L0-frame bbox into the window (mirrors the L1 expand in s_amr_regrid_shape_boxes). Containment margin
    !! max(amr_buf, 4) + amr_cpat_mar: clamping the tag to the window (the parent inset by amr_cpat_mar) can eat up to amr_cpat_mar
    !! of the body's stencil margin, and the parent was widened by (amr_max_level-1)*amr_cpat_mar so the window still holds the body
    !! plus max(amr_buf, 4): the C/F boundary sits a full image-point stencil off the surface, in fluid.
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
            call s_amr_body_bbox(ib_i, max(amr_buf, 4) + amr_cpat_mar, bb_lo, bb_hi)  ! global L0 cells, same frame as mlo/mhi
            bb_lo = max(bb_lo, mlo); bb_hi = min(bb_hi, mhi)
            if (f_amr_nest_window_empty(bb_lo, bb_hi)) cycle
            covered = .true.
            gwin(bb_lo(1):bb_hi(1),bb_lo(2):bb_hi(2),bb_lo(3):bb_hi(3)) = .true.
        end do

    end subroutine s_amr_nest_tag_parent

    !> Route each parent's (index, kb) pairs to the parent's clustering owner: send volume is O(this rank's tagged cells) and
    !! receive volume O(its parents' tags). The pairs are bucketed by owner with a stable counting sort (pass 1 appends in kb order
    !! and round-robin ownership interleaves the destinations). Serial: the send list is the gathered list.
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

    !> Pass 2 for one parent: rebuild its dense window from the routed pairs with gkb == kb (setting .true. once per cell dedups
    !! replicated tags) and extract the tagged cells in (k,j,i) order, so ctags does not depend on arrival order. The int8 decode
    !! matches the s_amr_pack_gwin_pairs encode; the remainder spans an xy plane, which can exceed 2**31 cells.
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

    !> Cluster a parent's tagged L0 cells into child boxes, pad by amr_buf and clamp into the nesting window. IB: a child clustered
    !! from the (widened) body tag must fully contain every overlapping body, so expand over bodies (mirrors the L1 expand in
    !! s_amr_regrid_shape_boxes) and re-clamp; the window already holds the body plus max(amr_buf, 4) (s_amr_nest_tag_parent), so
    !! the re-clamp does not cut the body's stencil.
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
        if (f_amr_nest_window_empty(clo, chi)) return  ! the inset left no interior
        call s_amr_nest_emit(clo, chi, lev, kb, mych, nmych)

    end subroutine s_amr_nest_inset_child

    !> Pad a child by `pad` and clamp it into the nesting window; collapsed dims are pinned to 0.
    pure subroutine s_amr_nest_clamp(clo, chi, mlo, mhi, pad)

        integer, intent(inout) :: clo(3), chi(3)
        integer, intent(in)    :: mlo(3), mhi(3), pad

        clo = merge(max(clo - pad, mlo), 0, amr_dim); chi = merge(min(chi + pad, mhi), 0, amr_dim)

    end subroutine s_amr_nest_clamp

    !> Tile a child to the level's slot cap and append the tiles to this rank's emission list. A level-lev block spans
    !! amr_ref_ratio**lev*(its L0 extent) fine cells while the slot holds amr_ref_ratio*amr_maxc_fit, so a child's L0 extent must be
    !! <= amr_maxc_fit/amr_ref_ratio**(lev-1), halving per level (a fixed /2 admits, at lev = 3, a box twice what the slot holds and
    !! corrupts the heap). The centred inset needs the same cap: it bounds the child as a fraction of its parent, not in absolute
    !! cells, and a parent of span 63 gives a child of span 33 against a level-2 cap of 32. A wider feature becomes adjacent
    !! sub-blocks like the L1 tiling: the level-aware fine-fine halo matches the shared seam flux and the reflux skips those faces.
    impure subroutine s_amr_nest_emit(clo, chi, lev, kb, mych, nmych)

        integer, intent(in)    :: clo(3), chi(3), lev, kb
        integer, intent(inout) :: mych(:,:), nmych
        type(t_box)            :: tiles(amr_max_blocks)
        integer                :: nt, capped, it

        nt = 0; capped = 0
        call s_amr_tile_box(clo, chi, tiles, nt, amr_max_blocks, capped, amr_maxc_fit/amr_ref_ratio**(lev - 1))
        do it = 1, nt
            if (nmych + 1 > amr_max_fine) exit
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

    ! 4) unchanged? (same count, boxes and levels as the live slots -> keep them; a rebuild would reproduce them exactly).
    ! The level must be compared too: a box that keeps its coordinates but changes refinement level would otherwise slip
    ! through with a stale amr_block_level, corrupting the level-aware coupling.
    impure subroutine s_amr_regrid_boxes_unchanged(boxes, nboxes, box_level, same)

        type(t_box), intent(in) :: boxes(:)
        integer, intent(in)     :: nboxes, box_level(:)
        logical, intent(out)    :: same
        integer                 :: k, ks

        ! regrid manages only the fine band [l0_slot_off+1 ..] of the shared pool; the level-0 L0-tile prefix (coexist) is not
        ! part of the box set, so compare against the fine block count and index slots through f_l0_slot.

        same = .false.
        if (nboxes == amr_num_blocks - l0_slot_off) then
            same = .true.
            do k = 1, nboxes
                ks = f_l0_slot(k)
                if (any(boxes(k)%lo /= amr_slots(ks)%region%lo) .or. any(boxes(k)%hi /= amr_slots(ks)%region%hi) .or. box_level(k) &
                    & /= amr_block_level(ks)) same = .false.
            end do
        end if

    end subroutine s_amr_regrid_boxes_unchanged

    !> Device pack of an owned old block's stash into its wire-pool slice (wp wire, stp store): the store is device-authoritative
    !! during the rebuild, so pack where the data lives; the slice stays on the device when the pools are device-resident
    !! (amr_fw_dev) and is copied out otherwise. Wire layout: gi fastest, then gj, gk, ii.
    impure subroutine s_amr_mig_pack_device(loc, e1, e2, e3, buf)

        integer, intent(in)                 :: loc, e1, e2, e3
        real(wp), intent(inout), contiguous :: buf(:)
        integer                             :: ii, gk, gj, gi, n1, n2, n3

        n1 = e1 + 1; n2 = e2 + 1; n3 = e3 + 1
        $:GPU_PARALLEL_LOOP(collapse=4, copyout='[buf]')
        do ii = 1, sys_size
            do gk = 0, e3
                do gj = 0, e2
                    do gi = 0, e1
                        buf(1 + gi + n1*(gj + n2*(gk + n3*(ii - 1)))) = real(amr_stor_st(gi, gj, gk, ii, loc), wp)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_mig_pack_device

    !> Device unpack of a received old block's pool slice into its stash replica (mirror of the pack).
    impure subroutine s_amr_mig_unpack_device(loc, e1, e2, e3, buf)

        integer, intent(in)              :: loc, e1, e2, e3
        real(wp), intent(in), contiguous :: buf(:)
        integer                          :: ii, gk, gj, gi, n1, n2, n3

        n1 = e1 + 1; n2 = e2 + 1; n3 = e3 + 1
        $:GPU_PARALLEL_LOOP(collapse=4, copyin='[buf]')
        do ii = 1, sys_size
            do gk = 0, e3
                do gj = 0, e2
                    do gi = 0, e1
                        amr_stor_st(gi, gj, gk, ii, loc) = real(buf(1 + gi + n1*(gj + n2*(gk + n3*(ii - 1)))), stp)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_mig_unpack_device

    !> Device cons->stor stash copy of one owned old block's fine interior (the store is device-authoritative; no host staging).
    impure subroutine s_amr_stash_copy_device(loc, e1, e2, e3)

        integer, intent(in) :: loc, e1, e2, e3
        integer             :: ii, gk, gj, gi

        $:GPU_PARALLEL_LOOP(collapse=4)
        do ii = 1, sys_size
            do gk = 0, e3
                do gj = 0, e2
                    do gi = 0, e1
                        amr_stor_st(gi, gj, gk, ii, loc) = amr_cons_st(gi, gj, gk, ii, loc)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_stash_copy_device

    !> Device overlap carry-forward: overwrite the new block's prolonged cons (slot column vloc, extents vm/vn/vp) with the covering
    !! old block's stashed fine detail (column vold, extents ve*), shifted by sh. Cells outside the old extent are skipped.
    impure subroutine s_amr_overlap_copy_device(vloc, vold, vm, vn, vp, sh, ve1, ve2, ve3)

        integer, intent(in) :: vloc, vold, vm, vn, vp, sh(3), ve1, ve2, ve3
        integer             :: i, fi, fj, fk, ofi, ofj, ofk, sh1, sh2, sh3
        logical             :: vd2, vd3

        sh1 = sh(1); sh2 = sh(2); sh3 = sh(3)
        vd2 = n_glb > 0; vd3 = p_glb > 0
        $:GPU_PARALLEL_LOOP(collapse=4, private='[ofi, ofj, ofk]')
        do i = 1, sys_size
            do fk = 0, vp
                do fj = 0, vn
                    do fi = 0, vm
                        ofk = fk + sh3; ofj = fj + sh2; ofi = fi + sh1
                        if (vd3 .and. (ofk < 0 .or. ofk > ve3)) cycle
                        if (vd2 .and. (ofj < 0 .or. ofj > ve2)) cycle
                        if (ofi < 0 .or. ofi > ve1) cycle
                        amr_cons_st(fi, fj, fk, i, vloc) = amr_stor_st(ofi, ofj, ofk, i, vold)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_overlap_copy_device

    !> Regrid phase 5: stash every live slot's fine interior (dead-between-steps q_cons_stor bounce), record the old block set
    !! (old_*), commit the new regions/levels/owners, and migrate each stashed old block point-to-point to the ranks that now own an
    !! overlapping new block.
    impure subroutine s_amr_regrid_stash_migrate(boxes, nboxes, box_level, old_np, old_ilo, old_ext, old_level, old_owns)

        type(t_box), intent(in) :: boxes(:)
        integer, intent(in)     :: nboxes, box_level(:)
        integer, intent(out)    :: old_np, old_ilo(:,:), old_ext(:,:), old_level(:)
        logical, intent(out)    :: old_owns(:)
        integer                 :: old_chi(3, amr_max_blocks), old_owner(amr_max_blocks)
        integer                 :: k, ks, kk, k2, rr, cnt, ix, lo, hi
        logical                 :: getk(amr_max_blocks), isdest(0:num_procs - 1)

        ! 5) stash every live slot's fine interior (dead-between-steps q_cons_stor bounce), keeping its old intersection origin

        ! old_* are indexed in the regrid's own dense fine-block space [1..old_np], which maps to shared-pool slot f_l0_slot(k);
        ! under coexist the level-0 L0-tile prefix [1..l0_slot_off] is not regrid-managed and must not be stashed or migrated.

        old_np = amr_num_blocks - l0_slot_off
        do k = 1, old_np
            ks = f_l0_slot(k)
            ! global block origin + extents (replicated, valid on every rank; not the owner-only isect), so the cross-rank
            ! migration below and the overlap-copy's index shift are correct even where this rank did not own the old block
            old_ilo(:,k) = amr_region_lo_all(:,ks)
            old_chi(:,k) = amr_region_hi_all(:,ks)  ! old coarse hi (for the P2P migration overlap test below)
            ! fine extent = (amr_ref_ratio**level)*footprint - 1: a level-2 block is 4x its L0 footprint, so stashing/migrating
            ! it with the level-1 factor (2x) would truncate half its fine cells.
            old_ext(1, k) = (amr_ref_ratio**amr_block_level(ks))*(amr_region_hi_all(1, ks) - amr_region_lo_all(1, ks) + 1) - 1
            old_ext(2, k) = merge((amr_ref_ratio**amr_block_level(ks))*(amr_region_hi_all(2, ks) - amr_region_lo_all(2, &
                    & ks) + 1) - 1, 0, n_glb > 0)
            old_ext(3, k) = merge((amr_ref_ratio**amr_block_level(ks))*(amr_region_hi_all(3, ks) - amr_region_lo_all(3, &
                    & ks) + 1) - 1, 0, p_glb > 0)
            old_owner(k) = amr_block_owner(ks)
            ! overlap-copy must match levels: an old L2's stash is in the 4x parent-fine frame
            old_level(k) = amr_block_level(ks)
            old_owns(k) = amr_owns_all(ks)
            if (old_owns(k)) then
                ! Device-side stash: the store is device-authoritative, so copy cons->stor where the data lives instead of staging
                ! two full-slot transfers through the host. A mid-rebuild grow's device->host round trip preserves this stash by
                ! construction (s_amr_st_reserve's contract); the host mirror of amr_stor_st stays stale, which is fine because
                ! the migration pack and the overlap carry-forward below are device kernels and no host reader of the stash
                ! exists. (The kernel lives in its own subroutine: amdflang drops target regions nested in BLOCK constructs from
                ! the device image, and the first launch then dies on HSA_STATUS_ERROR_INVALID_SYMBOL_NAME.)
                call s_amr_stash_copy_device(amr_loc_of(ks), old_ext(1, k), old_ext(2, k), old_ext(3, k))
            end if
        end do

        ! set the regions + assign owners before the migration (P2P needs the new owners) and before the owner-dependent
        ! geometry (else s_set_amr_fine_geometry sizes the whole-block owner from a stale amr_block_owner)
        ! the fine band ends at f_l0_slot(nboxes); the level-0 tile prefix below it keeps its regions, levels and owners (a plain
        ! nboxes here would overrun the prefix under coexist)
        amr_num_blocks = f_l0_slot(nboxes)
        do k = 1, nboxes
            ks = f_l0_slot(k)
            amr_region_lo_all(:,ks) = boxes(k)%lo; amr_region_hi_all(:,ks) = boxes(k)%hi
            ! box_level(k) is the refinement level assigned during the hierarchical nesting above (1 for L0->L1 boxes, l for a
            ! box nested at level l). Setting it every regrid resets a stale level when a slot is reused across levels.
            amr_block_level(ks) = box_level(k)
        end do
        ! block set changed: dirty the cached seam-pair and overlap-rank lists now; the rebuild's per-block P2P gathers
        ! (s_amr_regrid_rebuild_slots) consume the overlap lists with the new boxes, so flagging after them would be too late
        amr_seam_pairs_dirty = .true.
        amr_mesh_epoch = amr_mesh_epoch + 1
        ! Proper-nesting guard: each level>=2 block must be covered by exactly one parent-level block. f_amr_parent_block (and
        ! the gather/reflux that key off it) take the first overlap, so a fine tile straddling two parent tiles (an internal
        ! parent-level tile seam crossed by a nested feature) would silently couple to only one parent (wrong coarse BC + a
        ! conservation leak on the other). Abort fail-closed instead. Replicated boxes -> every rank aborts together.
        block
            integer :: bk, bkk, npar
            do bk = 1, nboxes
                if (box_level(bk) < 2) cycle
                npar = 0
                do bkk = 1, nboxes
                    if (box_level(bkk) == box_level(bk) - 1 .and. f_amr_boxes_overlap(boxes(bk)%lo, boxes(bk)%hi, boxes(bkk)%lo, &
                        & boxes(bkk)%hi)) npar = npar + 1
                end do
                if (npar /= 1) call s_mpi_abort('amr multi-level: a level>=2 block overlaps more than one (or no) ' &
                    & // 'parent-level block - a fine tile straddling a parent-tile seam is unsupported (gather/reflux ' &
                    & // 'couple to a single parent); reduce max_grid_size or the refined feature extent')
            end do
        end block
        amr_num_levels = maxval(box_level(1:nboxes))
        call s_amr_assign_block_owners()
        ! The partition is decided here and nothing has moved yet; everything below redistributes data.

        ! Cross-rank fine-state migration as one wave: each old owner ships its stashed fine state to every distinct new-block
        ! owner whose region overlaps that old block, and the overlap-copy then reads every covering old block regardless of who
        ! owned it. A received old block lands in a stash-only replica slot (freed by the rebuild's early-free or the reconcile;
        ! a replica never touches q_prim/rhs, and full slots across a migration-heavy regrid's replica set would exhaust device
        ! memory). Both sides enumerate old blocks ascending, so each peer pair's transfers line up. No-op at np=1.
        if (num_procs > 1) then
            call s_amr_wave_open(amr_wave, 4)
            call s_amr_wave_reset(amr_wsend); call s_amr_wave_reset(amr_wrecv)
            do kk = 1, old_np
                cnt = sys_size*(old_ext(1, kk) + 1)*(old_ext(2, kk) + 1)*(old_ext(3, kk) + 1)
                isdest = .false.
                do k2 = 1, nboxes
                    rr = amr_block_owner(f_l0_slot(k2))
                    if (f_amr_boxes_overlap(boxes(k2)%lo, boxes(k2)%hi, old_ilo(:,kk), old_chi(:,kk))) isdest(rr) = .true.
                end do
                getk(kk) = isdest(proc_rank) .and. .not. old_owns(kk)
                if (getk(kk)) call s_amr_wave_add(amr_wrecv, old_owner(kk), kk, old_ilo(:,kk), old_chi(:,kk), cnt)
                if (.not. old_owns(kk)) cycle
                do rr = 0, num_procs - 1
                    if (isdest(rr) .and. rr /= proc_rank) call s_amr_wave_add(amr_wsend, rr, kk, old_ilo(:,kk), old_chi(:,kk), cnt)
                end do
            end do
            call s_amr_wave_close(amr_wsend, amr_fw_sq, amr_fw_dev)
            call s_amr_wave_close(amr_wrecv, amr_fw_rq, amr_fw_dev)
            call s_amr_prereserve_stash(getk, old_np)
            do kk = 1, old_np
                if (getk(kk)) call s_amr_alloc_slot_stash(f_l0_slot(kk))
            end do
            call s_amr_wave_post(amr_wave, amr_wrecv, amr_fw_rq, XA_F4_RCV, amr_fw_dev)
            do ix = 1, amr_wsend%nx
                kk = amr_wsend%blk(ix)
                call s_amr_wave_slice(amr_wsend, ix, lo, hi)
                call s_amr_mig_pack_device(amr_loc_of(f_l0_slot(kk)), old_ext(1, kk), old_ext(2, kk), old_ext(3, kk), &
                                           & amr_fw_sq(lo:hi))
                call s_amr_wave_hdr_pack(amr_wsend, amr_fw_sq, ix, XA_F4_SND)
            end do
            call s_amr_wave_send(amr_wave, amr_wsend, amr_fw_sq, XA_F4_SND, amr_fw_dev)
            call s_amr_wave_wait(amr_wave)
            do ix = 1, amr_wrecv%nx
                kk = amr_wrecv%blk(ix)
                call s_amr_wave_hdr_check(amr_wrecv, amr_fw_rq, ix, XA_F4_SND)
                call s_amr_wave_slice(amr_wrecv, ix, lo, hi)
                call s_amr_mig_unpack_device(amr_loc_of(f_l0_slot(kk)), old_ext(1, kk), old_ext(2, kk), old_ext(3, kk), &
                                             & amr_fw_rq(lo:hi))
            end do
        end if

    end subroutine s_amr_regrid_stash_migrate

    !> Regrid phase 6: build each new slot (geometry, prolong from its coarse patch, overwrite the overlap from every covering
    !! stashed old block), level by level as fill waves (the level-1 wave from the coarse grid, then one parent wave per level, so
    !! every parent is built before its children read it), then reconcile the slot pool, rebuild the fine IB state and re-validate
    !! the seam topology. Old-only slots are freed right after the last owned box whose region overlaps them is built (region
    !! overlap is a superset of every per-cell stash read), so their dense indices recycle into the next allocs instead of the whole
    !! stash/replica set peaking device memory at np >= 2.
    impure subroutine s_amr_regrid_rebuild_slots(q_cons_base, boxes, nboxes, old_np, old_ilo, old_ext, old_level, old_owns)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_base
        type(t_box), intent(in)                                :: boxes(:)
        integer, intent(in)                                    :: nboxes, old_np, old_ilo(:,:), old_ext(:,:), old_level(:)
        logical, intent(in)                                    :: old_owns(:)
        integer                                                :: sh(3), k, kk, i, hh, ks, kks, lev, nh, ohi(3), pos, nlev
        integer, allocatable                                   :: last_use(:), held(:), held_hi(:,:), vpos(:)

        ! non-owner geometry for every box: amr_owns_all = F, empty footprint, -1 fine extents, the replicated state every rank
        ! must agree on (s_amr_select_slot reads it for any block); box k lives in shared-pool slot ks = f_l0_slot(k)

        do k = 1, nboxes
            ks = f_l0_slot(k)
            if (amr_block_owner(ks) == proc_rank) cycle
            amr_cur = ks
            call s_set_amr_fine_geometry(boxes(k)%lo, boxes(k)%hi)
        end do

        ! the visit order (owned boxes, level-major, ascending) and each held old block's last reader in it; nlev is replicated,
        ! so every rank enters every level's wave (a rank can be a parent-owner sender without owning a box at that level)
        call s_amr_refresh_my_blocks()
        nlev = maxval(amr_block_level(f_l0_slot(1):amr_num_blocks))
        allocate (vpos(nboxes)); vpos = 0; pos = 0
        do lev = 1, nlev
            do i = 1, amr_n_my
                k = amr_my_blk(i) - l0_slot_off
                if (k < 1) cycle  ! L0 tile prefix
                if (amr_block_level(amr_my_blk(i)) /= lev) cycle
                pos = pos + 1; vpos(k) = pos
            end do
        end do
        allocate (last_use(old_np), held(old_np), held_hi(3, old_np)); last_use = 0
        nh = 0
        do kk = 1, old_np
            if (.not. amr_slot_live(f_l0_slot(kk))) cycle
            nh = nh + 1; held(nh) = kk
            ohi = old_ilo(:,kk) + merge((old_ext(:,kk) + 1)/amr_ref_ratio**old_level(kk) - 1, 0, amr_dim)
            held_hi(:,nh) = ohi
            do i = 1, amr_n_my
                k = amr_my_blk(i) - l0_slot_off
                if (k < 1) cycle
                if (f_amr_boxes_overlap(boxes(k)%lo, boxes(k)%hi, old_ilo(:,kk), ohi)) last_use(kk) = max(last_use(kk), vpos(k))
            end do
        end do

        do lev = 1, nlev
            if (lev == 1) then
                call s_amr_l1_fill_exchange(q_cons_base, .true.)
            else
                call s_amr_parent_fill_exchange(lev, .true.)
            end if
            do i = 1, amr_n_my
                ks = amr_my_blk(i)
                k = ks - l0_slot_off
                if (k < 1) cycle
                if (amr_block_level(ks) /= lev) cycle
                amr_cur = ks
                ! free the old-only slots no later box reads (s_amr_free_slot is a no-op once dead); a slot serving as a new
                ! owned box keeps living - the reconcile decides it
                do hh = 1, nh
                    kk = held(hh)
                    if (last_use(kk) >= vpos(k)) cycle
                    kks = f_l0_slot(kk)
                    if (kks <= amr_num_blocks) then
                        if (amr_block_owner(kks) == proc_rank) cycle
                    end if
                    call s_amr_free_slot(kks)
                end do
                call s_amr_alloc_slot(ks)
                call s_set_amr_fine_geometry(boxes(k)%lo, boxes(k)%hi)
                call s_amr_select_slot(ks)
                if (lev == 1) then
                    call s_amr_l1_fill_consume(q_cons_base, ks, .true.)
                else
                    call s_amr_parent_fill_consume(ks, .true.)
                end if
                ! prolong and overlap carry-forward are both device kernels: the slot is built in place where the store is
                ! authoritative. A level>=2 block re-prolongs from its (freshly built, parents-first) parent each regrid: its
                ! stash is in the parent-fine frame, so the L0-frame shift below does not apply; the coupling keeps conservation.
                call s_interpolate_coarse_to_fine()
                if (lev >= 2) cycle
                do hh = 1, nh
                    kk = held(hh)
                    if (old_level(kk) /= 1) cycle  ! same-level overlap only (a child's stash is 4x-framed)
                    if (.not. f_amr_boxes_overlap(boxes(k)%lo, boxes(k)%hi, old_ilo(:,kk), held_hi(:,hh))) cycle
                    kks = f_l0_slot(kk)
                    sh = amr_ref_ratio*(amr_isect_lo - old_ilo(:,kk))  ! old local fine index = new local fine index + sh
                    call s_amr_overlap_copy_device(amr_loc_of(ks), amr_loc_of(kks), amr_slots(ks)%m, amr_slots(ks)%n, &
                                                   & amr_slots(ks)%p, sh, old_ext(1, kk), old_ext(2, kk), old_ext(3, kk))
                end do
            end do
            call s_amr_fill_wave_done()
        end do
        deallocate (last_use, held, held_hi, vpos)

        ! one allreduce for the whole loop; sets amr_xchg_coarse_ghosts if any block needs it
        call s_amr_reduce_xchg_flag()
        ! lazy sizing: free the transient regrid slots (old blocks this rank stashed/received but does not now own); the
        ! new-owned slots were allocated in the build loop, so this only frees; a rank keeps just its owned blocks' fine arrays
        call s_amr_reconcile_slots()
        ! rebuild every block's fine-grid IB state for the new geometry (markers/ghost points/image points recomputed from the
        ! body definitions; no state carries across regrids)
        if (ib) call s_amr_setup_ib()
        call s_amr_select_slot(1)
        call s_amr_check_seam_topology()

    end subroutine s_amr_regrid_rebuild_slots

    !> Sensor-on-fine child tagging: OR-accumulate density-gradient tags from an old fine block's solution into an L0-cell tag grid,
    !! restricted to a parent nesting window. Reads the flat store on the host (caller host-refreshes the block first; the stash's
    !! GPU_UPDATE runs later). Fine cell (fi,fj,fk) covers L0 cell (ci,cj,ck) with fi = rr*(ci-olo(1))+d etc.; the gradient uses
    !! one-sided differences at the fine-interior edges so no stale fine ghost is read. Only decides placement; conservation is
    !! enforced downstream by restrict/reflux regardless of box extent.
    impure subroutine s_amr_tag_child_from_fine(ob, win_lo, win_hi, ctag, any_tag)

        integer, intent(in)    :: ob, win_lo(3), win_hi(3)
        logical, intent(inout) :: ctag(win_lo(1):,win_lo(2):,win_lo(3):)
        logical, intent(inout) :: any_tag
        integer                :: rr, ci, cj, ck, fi, fj, fk, d1, d2, d3, fm1, fm2, fm3, olo(3), lo(3), hi(3)
        real(wp)               :: r0, g
        logical                :: tagged

        rr = amr_slots(ob)%amr_ref_ratio
        olo = amr_region_lo_all(:,ob)
        fm1 = amr_slots(ob)%m; fm2 = amr_slots(ob)%n; fm3 = amr_slots(ob)%p
        ! overlap of this old block with the parent window, in L0 cells
        lo(1) = max(win_lo(1), amr_region_lo_all(1, ob)); hi(1) = min(win_hi(1), amr_region_hi_all(1, ob))
        lo(2) = merge(max(win_lo(2), amr_region_lo_all(2, ob)), 0, n_glb > 0)
        hi(2) = merge(min(win_hi(2), amr_region_hi_all(2, ob)), 0, n_glb > 0)
        lo(3) = merge(max(win_lo(3), amr_region_lo_all(3, ob)), 0, p_glb > 0)
        hi(3) = merge(min(win_hi(3), amr_region_hi_all(3, ob)), 0, p_glb > 0)
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
                                r0 = max(abs(f_amr_rho_tot_st(amr_loc_of(ob), fi, fj, fk)), 1.e-30_wp)
                                g = abs(f_amr_rho_tot_st(amr_loc_of(ob), min(fi + 1, fm1), fj, &
                                        & fk) - f_amr_rho_tot_st(amr_loc_of(ob), max(fi - 1, 0), fj, fk))
                                if (n_glb > 0) g = max(g, abs(f_amr_rho_tot_st(amr_loc_of(ob), fi, min(fj + 1, fm2), &
                                    & fk) - f_amr_rho_tot_st(amr_loc_of(ob), fi, max(fj - 1, 0), fk)))
                                if (p_glb > 0) g = max(g, abs(f_amr_rho_tot_st(amr_loc_of(ob), fi, fj, min(fk + 1, &
                                    & fm3)) - f_amr_rho_tot_st(amr_loc_of(ob), fi, fj, max(fk - 1, 0))))
                                ! 2*r0 normalizes the 2-cell central difference; the 2 is the stencil span, not amr_ref_ratio
                                if (g/(2._wp*r0) > amr_tag_eps) tagged = .true.
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

    !> Total density (sum of the continuity variables) at one cell: the regrid tag field. Reduces to variable 1 for one fluid. Two
    !! sources, one body: `_st` reads a refined block out of the flat store, `_sf` the level-0 monolithic field.
    #:for RSFX, RSRC in [('st', 'amr_cons_st'), ('sf', '')]
        pure function f_amr_rho_tot_${RSFX}$(${'loc' if RSRC else 'q'}$, ci, cj, ck) result(r)

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
    #:endfor
end module m_amr_regrid
