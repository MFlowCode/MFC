!>
!!@file
!!@brief Contains module m_amr_regrid

#:include 'macros.fpp'

!> @brief Dynamic regrid for the block-structured AMR level set: density-gradient tagging, Berger-Rigoutsos clustering, box shaping
!! (pad/clamp/tile/IB merge), hierarchical child nesting, slot rebuild with cross-rank fine-state migration. Block/slot state lives
!! in m_amr (and m_global_parameters); this module only drives it.
module m_amr_regrid

#ifdef MFC_MPI
    use mpi  !< per-node signature ALLREDUCE for the rank-invariant clustering, point-to-point for the fine-state migration
#endif

    use m_derived_types  ! scalar_field, t_box
    use m_box, only: f_morton  ! B1: canonical merge order (shared 3D Morton key)
    use m_global_parameters
    use m_constants, only: mapCells
    use m_mpi_proxy, only: s_mpi_abort
    use m_mpi_common, only: s_mpi_allreduce_min, s_mpi_allreduce_max
    use m_phase_timing, only: s_phase_tic, s_phase_toc, PH_RGHALO, PH_RGTAG, PH_RGCLUS, PH_RGSHAPE, PH_RGMIG, PH_RGBUILD, &
        & PH_RGPART, PH_RGMOVE, PH_MGWAIT, PH_RBGATH, PH_RBOVL, PH_RBPUSH, PH_RBSLOT, PH_RBGEO, PH_RBTAIL, PH_RBFLUSH, PH_RBXCHG, &
        & PH_RBREC, PH_RBTOPO, PH_MGSLOT, PH_MGPACK, PH_MGUNPK, PH_MGPUSH
    use m_amr, only: s_amr_build_gather_plan, amr_gpl_valid, amr_slots, amr_cons_st, amr_stor_st, amr_loc_of, &
        & s_amr_gather_chunk_post, s_amr_gather_chunk_send, s_amr_gather_consume_box, amr_gath_chunk, s_amr_cov_note, &
        & amr_maxc_fit, amr_seam_pairs_dirty, amr_mesh_epoch, amr_xchg_coarse_ghosts, amr_cpat_mar, s_amr_alloc_slot, &
        & s_amr_alloc_slot_stash, s_amr_prereserve_stash, s_amr_free_slot, s_amr_reduce_xchg_flag, s_amr_reconcile_slots, &
        & s_amr_assign_block_owners, s_amr_gather_send_flush, s_amr_gather_coarse_patch_pbmv, s_amr_prolong_pbmv, &
        & s_amr_exchange_coarse_cons_halo, s_lag_phys_to_cells, s_amr_body_bbox, s_amr_expand_box_over_bodies, s_amr_tile_box, &
        & f_amr_seam_dim, f_amr_boxes_overlap, s_set_amr_fine_geometry, s_interpolate_coarse_to_fine, s_amr_setup_ib, f_l0_slot, &
        & amr_gb_tag, amr_gb_win, amr_gb_cost, amr_gb_mig, amr_mig_snd, amr_mig_blk, amr_cad_tot, amr_cad_esc, amr_cad_armed, &
        & amr_cl_maxdep, amr_cl_maxdep_leaf, amr_cl_lmax, amr_cl_ldepth, amr_cl_nodes, amr_cl_rb, amr_cl_rb_now, &
        & amr_cl_shr_nodes, amr_cl_shr_rb, amr_cl_loc_nodes, amr_cl_loc_rb, amr_cl_shr_maxdep, s_amr_ranks_overlapping, &
        & amr_cl_shr_nodes_r, amr_cl_shr_rb_r, amr_cl_loc_nodes_r, amr_cl_loc_rb_r, amr_cl_shr_maxdep_r, amr_cl_me_nodes_r, &
        & amr_cl_me_rb_r, amr_my_blk, amr_n_my, s_amr_refresh_my_blocks, s_amr_fw_szi, f_amr_overlap_count, f_amr_rank_overlaps, &
        & amr_tag_base, amr_mesh_epoch, amr_cl_wire_r, amr_gb_box
    use m_amr_xchg_audit, only: s_xa_rec, XA_F4_SND, XA_F4_RCV  ! I1a exchange accounting (migration family)
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
    !! changes (regrid, restart); O(nblocks^2) on the replicated region metadata, identical on every rank. Two cases: (a) adjacency
    !! WITHOUT the exact transverse match f_amr_seam requires - reachable only via IB body-bbox expansion (clustering merges any
    !! too-close pair; tiling emits a regular grid) - which the fine-fine halo can never pair; (b) same-level box INTERSECTION -
    !! reachable only via CHILD IB body-bbox expansion, which unlike the L1 path has no overlap-merge pass -
    !! double-restricting/refluxing the shared cells.
    impure subroutine s_amr_check_seam_topology()

        integer :: xb, yb, d, t
        logical :: adj, tover

        do xb = 1, amr_num_blocks
            do yb = 1, amr_num_blocks
                if (xb == yb) cycle
                if (amr_block_level(xb) /= amr_block_level(yb)) cycle
                ! same-level INTERSECTION (different levels legitimately nest; tiling emits disjoint tiles, the L1 IB pass merges
                ! overlapping boxes, but the CHILD IB body-bbox expansion has no overlap-merge pass)
                if (f_amr_boxes_overlap(amr_region_lo_all(:,xb), amr_region_hi_all(:,xb), amr_region_lo_all(:,yb), &
                    & amr_region_hi_all(:,yb))) then
                    call s_mpi_abort("AMR: two same-level blocks INTERSECT (the child IB body-bbox expansion route can " &
                                     & // "produce this - it has no overlap-merge pass): the overlapping cells would be " &
                                     & // "restricted and refluxed twice, silently breaking conservation. Adjust the body/" &
                                     & // "regrid inputs so body-expanded child boxes merge or separate.")
                end if
                do d = 1, num_dims
                    ! relaxed adjacency: touching faces in dim d with ANY transverse overlap
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

    !> Abort if any box exceeds the slot cap for its level. The slot coord/field arrays are allocated ONCE to
    !! amr_ref_ratio*amr_maxc_fit fine cells, and a level-lev block spans amr_ref_ratio**lev fine cells per coarse cell, so its
    !! coarse extent must be <= amr_maxc_fit/amr_ref_ratio**(lev-1). Every emitter enforces that via s_amr_tile_box; this checks the
    !! invariant once, where the box set is final, rather than trusting each emitter to have done it.
    !!
    !! It exists because a violation is otherwise SILENT AND CATASTROPHIC: s_amr_build_block_coords sizes the fine coords from the
    !! block's true extent, so an over-cap box writes past x_cb, corrupting the heap on EVERY regrid and surfacing much later as
    !! "corrupted size vs. prev_size" inside an unrelated free(). One emitter did skip the cap (cfbacebe, the brand-new-region
    !! branch) and finding it took a bounds-checked rebuild and a multi-hour hunt, because the crash site was nowhere near the bug.
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

    !> Invariant check (plan-based exchange I0): same-level boxes are pairwise DISJOINT. Guaranteed today by the cluster partition +
    !! merge threshold + IB overlap-merge; relied on by the rebuild's overlap carry-forward and by per-peer unpack reordering in the
    !! exchange plans (amr_plan_based_exchange.md) - enforced here rather than inherited as folklore. All levels share the global
    !! coarse index space (see the cap formula in s_amr_check_box_caps), so the interval test is valid across parents. O(nboxes^2)
    !! host integer compares per regrid - negligible at current box counts; replace with a sorted sweep in increment I7 if box
    !! counts grow.
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

    !> Concatenated 1D tag signatures of box [blo0:bhi0], built from the tag range [ts:te] in ONE pass. Axis d occupies sig(off(d) :
    !! off(d) + ext(d) - 1), and sig(off(d) + t - blo0(d)) counts the in-box tagged cells at position t along d. One signature
    !! serves the trim, the in-box count AND every candidate split, replacing the up-to-num_dims+1 separate rescans of the tag list
    !! the previous form needed. nsig returns the used length.
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
    !! [blo0:bhi0]. Equivalent to scanning the tag list: the per-axis MIN/MAX of the contained tags ARE the first and last nonzero
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
    !! TRIMMED range is exact: trim shrinks only to the tags' own bbox, so no tag leaves the box. Integer-only => identical on all
    !! ranks.
    pure subroutine s_amr_find_split_sig(sig, off, blo0, blo, bhi, sax, spos, ok)

        integer, intent(in)  :: sig(:), off(3), blo0(3)
        integer, intent(in)  :: blo(3), bhi(3)
        integer, intent(out) :: sax, spos
        logical, intent(out) :: ok
        !> Minimum child extent along the split axis, i.e. the smallest box the bisection may produce. 2 is the algorithmic floor;
        !! amr_blocking_factor raises it, which is what stops the bisection over-generating. Measured 2026-08-27: with the floor at
        !! 2 and amr_cluster_eff = 0.9 the recursion never converges on its own -- it splits until the amr_max_blocks cap stops it
        !! (warning on EVERY regrid at five different caps), and the min-separation merge then collapses the result back. An 8x cap
        !! bought 61%% more tree nodes for a 0.7%% change in the final box set. NOTE this is a minimum SIZE, not AMReX's blocking
        !! factor: AMReX coarsens the TAG LATTICE, which also shrinks its global tag gather. S3.1 already deleted that gather here,
        !! and coarsening a rank-local sparse list cannot dedup coarse cells that straddle a rank boundary without an extra
        !! exchange, so the size floor is both simpler and the part that actually stops the over-generation.
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
        logical                :: ovl

        do s = 1, num_source
            if (hi(1) < lo(1) .or. hi(2) < lo(2) .or. hi(3) < lo(3)) return  ! emptied by an earlier clip
            ovl = lo(1) <= acoustic_supp_hi(1, s) .and. hi(1) >= acoustic_supp_lo(1, s)
            if (n_glb > 0) ovl = ovl .and. lo(2) <= acoustic_supp_hi(2, s) .and. hi(2) >= acoustic_supp_lo(2, s)
            if (p_glb > 0) ovl = ovl .and. lo(3) <= acoustic_supp_hi(3, s) .and. hi(3) >= acoustic_supp_lo(3, s)
            if (.not. ovl) cycle
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
        ! safety net: clipping removed every overlap by construction - anything left is a bug
        do s = 1, num_source
            if (hi(1) < lo(1) .or. hi(2) < lo(2) .or. hi(3) < lo(3)) return
            ovl = lo(1) <= acoustic_supp_hi(1, s) .and. hi(1) >= acoustic_supp_lo(1, s)
            if (n_glb > 0) ovl = ovl .and. lo(2) <= acoustic_supp_hi(2, s) .and. hi(2) >= acoustic_supp_lo(2, s)
            if (p_glb > 0) ovl = ovl .and. lo(3) <= acoustic_supp_hi(3, s) .and. hi(3) >= acoustic_supp_lo(3, s)
            if (ovl) call s_mpi_abort('amr regrid: acoustic source exclusion clip failed (internal error)')
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
        logical                :: ovl

        if (hi(1) < lo(1) .or. hi(2) < lo(2) .or. hi(3) < lo(3)) return
        ovl = lo(1) <= shi(1) .and. hi(1) >= slo(1)
        if (n_glb > 0) ovl = ovl .and. lo(2) <= shi(2) .and. hi(2) >= slo(2)
        if (p_glb > 0) ovl = ovl .and. lo(3) <= shi(3) .and. hi(3) >= slo(3)
        if (.not. ovl) return
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
    !! and the coarse RHS only computes fluxes inside it. The window only GROWS (s_grow_active_box monotone, self-disabling at full
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
            ok = amr_region_lo_all(1, k) > ab_x%beg .and. amr_region_hi_all(1, k) < ab_x%end
            if (n_glb > 0) ok = ok .and. amr_region_lo_all(2, k) > ab_y%beg .and. amr_region_hi_all(2, k) < ab_y%end
            if (p_glb > 0) ok = ok .and. amr_region_lo_all(3, k) > ab_z%beg .and. amr_region_hi_all(3, k) < ab_z%end
            if (.not. ok) then
                call s_mpi_abort('amr with active_box: an AMR block is not strictly inside the active ' &
                                 & // 'window; place the initial block (with a one-cell margin) inside the ' &
                                 & // 'initial non-ambient region plus buff_size')
            end if
        end do

    end subroutine s_amr_check_active_box_containment

    !> This rank's OWN tagged cells as GLOBAL level-0 coordinates. Each rank scans only its interior (0:m, 0:n, 0:p), which the
    !! level-0 decomposition makes disjoint, so no global cell is emitted twice and a SUM reduction over these lists counts every
    !! tagged cell exactly once. This replaces the ALLGATHERV of the global tag list (W4): per-rank memory and wire volume now scale
    !! with the rank's OWN tag count instead of the global one.
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

    !> Pack this rank's OWNED tagged cells of the child window [mlo:mhi] as (linear-index, kb) pairs, appended to the per-level send
    !! arrays sidx(:) (int8 linear index) / skb(:) (parent box id). The int8 encode matches the pass-2 decode, so gathering these
    !! pairs across ranks and setting them into a per-parent dense window reproduces the old dense-window dedup (replicated/
    !! overlapping tags collapse) and the (k,j,i) extraction order exactly -> byte-identical child boxes. One allgatherv per level
    !! (caller) drops the collective count from O(#parent-boxes) to O(#levels). gwin is read, not modified.
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

    !> Cluster a SPARSE tag list (level-0 cell coords, tags(1:3, 1:ntag_in)) into a LIST of separated block boxes, identically on
    !! every rank. Caller builds the list (s_amr_local_tags / s_amr_pack_gwin_pairs); per-rank memory is O(#tagged), not O(global
    !! grid). Berger-Rigoutsos recursive bisection until each box's tag efficiency reaches amr_cluster_eff (or it is atomic / the
    !! amr_max_blocks cap is hit), then merges any two boxes whose amr_buf-padded extents come within buff_size (so no fine-fine
    !! adjacency: separated boxes stay >= buff_size apart, nearby ones collapse to one box == the legacy bounding box). Boxes are
    !! raw tagged extents; the caller pads, clamps, size-caps each.
    impure subroutine s_amr_cluster(tags, ntag_in, boxes, nboxes, reduce)

        integer, intent(in) :: tags(:,:), ntag_in
        !> .true.: `tags` is this rank's LOCAL list and each node's signature is ALLREDUCEd, so the tree is driven by global counts
        !! without any rank holding the global tag list. .false.: `tags` is already replicated on every rank.
        logical, intent(in)                   :: reduce
        type(t_box), allocatable, intent(out) :: boxes(:)
        integer, intent(out)                  :: nboxes
        integer, allocatable                  :: slo(:,:), shi(:,:), alo(:,:), ahi(:,:)
        integer, allocatable                  :: sts(:), ste(:), wt(:,:)
        integer, allocatable                  :: sdep(:)  !< S3.0a: recursion depth carried with each stack entry
        integer                               :: dep, mxdep
        integer, allocatable                  :: sig(:)   !< concatenated per-axis tag signature of the node's box
        integer, allocatable                  :: ovr(:)   !< S3.2a scratch: ranks overlapping the node's box
        integer                               :: novr
        integer                               :: blo0(3), bhi0(3), off(3), nsig

#ifdef MFC_MPI
        integer :: ierr
#endif
        integer(8)              :: nnode
        integer                 :: mg, ng, pg, t
        integer                 :: cap, nacc, i, j, k, d, sax, spos, thr, ntag
        integer(8), allocatable :: akey(:)  !< B1: Morton key of each accepted box's lo, the canonical merge order
        integer, allocatable    :: gcnt(:), gdsp(:), sbx(:,:), gbx(:,:)  !< S3.2b: union of the per-rank accepted boxes
        integer                 :: ntot
        !> S3.2b-2: the level-order walk. kpos/kbat index the nodes kept at the current depth; bsig concatenates the signatures of
        !! that depth's SHARED nodes into the single buffer the one reduction covers, with bofs/blen/boff their slices.
        integer, allocatable :: kpos(:), kbat(:), bofs(:), blen(:), boff(:,:), bsig(:)
        integer              :: ncur, nnxt, nkeep, nbat, nbuf
        !> S3.2b-2b: a node is WIDE when its box spans more than this many ranks. Wide nodes keep the batched collective (every rank
        !! overlaps them and needs the answer); narrow ones reduce among their few overlapping ranks. The threshold only has to keep
        !! the WIDE COUNT at O(log P) -- measured, a rank participates in 6/8/10 shared nodes at np8/16/32 while the shared total is
        !! 11/23/47 -- and 8 is one 2x2x2 brick of ranks, the shape a seam node actually has.
        integer, parameter   :: amr_cl_wide = 8
        integer, allocatable :: bnov(:), bovr(:,:), wbuf(:)
        logical, allocatable :: bwide(:)
        integer, allocatable :: pidx(:), plist(:), scnt(:), rcnt(:), sdsp2(:), rdsp2(:), soff(:), roff(:)
        integer, allocatable :: sbuf(:), rbuf(:), creq(:)
        integer              :: np2, q, rr, nsnd, nrcv, nreq2, tagc, nwb, o1  ! t is already a loop variable above
        integer(8)           :: bkey
        integer(8)           :: vol  !< box volume; a global-bbox first pass can exceed 2**31 cells
        integer              :: blo(3), bhi(3), ts, te, lo, hi, tmp(3), tmp2(3)
        logical              :: ok, force, capped, changed, tooclose, mine
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
        allocate (ovr(amr_cl_wide))  ! S3.2b-2b: only NARROW nodes are ever enumerated, so this no longer sizes with P
        allocate (akey(cap))  ! B1 scratch
        ! working copy of the tag list, partitioned in place as the tree descends so each node scans only its tags
        do t = 1, ntag_in
            wt(:,t) = tags(:,t)
        end do
        ncur = 1; slo(:,1) = [0, 0, 0]; shi(:,1) = [mg, ng, pg]  ! first node trims to the global tagged bbox
        sts(1) = 1; ste(1) = ntag_in
        sdep(1) = 0; mxdep = 0; nnode = 0_8
        nacc = 0; capped = .false.
        allocate (kpos(4*cap + 8), kbat(4*cap + 8), bofs(4*cap + 8), blen(4*cap + 8), boff(3, 4*cap + 8))
        allocate (bsig(4*(mg + ng + pg + 3)))
        allocate (bnov(4*cap + 8), bwide(4*cap + 8), bovr(amr_cl_wide, 4*cap + 8))
        allocate (pidx(0:max(num_procs - 1, 0)), plist(max(num_procs, 1)))
        allocate (scnt(max(num_procs, 1)), rcnt(max(num_procs, 1)), sdsp2(max(num_procs, 1)), rdsp2(max(num_procs, 1)))
        allocate (soff(max(num_procs, 1)), roff(max(num_procs, 1)))
        allocate (wbuf(1), sbuf(1), rbuf(1), creq(1))
        pidx = 0
        ! S3.2b-2: LEVEL-ORDER descent. The old stack held mixed depths, so each SHARED node paid its own global reduction, and
        ! the shared set grows with P (measured per regrid: 11/23/47/91 at np8/16/32/64) -- O(P) full-machine syncs, ~147,000 at
        ! 1e5 ranks. Walking one whole depth at a time lets every shared node at that depth ride ONE reduction, so the count
        ! becomes O(tree depth): shr_maxdep is 2*log2(P) - 3 (3/5/7/9 over those same rungs), i.e. ~30 at 1e5 ranks.
        !
        ! WHY ONE COLLECTIVE IS EVEN LEGAL once S3.2b stopped replicating the tree: a child's box lies inside its parent's, so
        ! the ranks overlapping a child are a SUBSET of those overlapping its parent, and therefore every ancestor of a shared
        ! (novr > 1) node is itself shared. No rank ever drops a shared node's ancestor, so every rank walks the whole shared
        ! subtree, in the same deterministic order -- the per-depth batch is identical in content AND order on every rank, which
        ! is exactly what a single collective needs. Rank-local nodes differ per rank and are excluded from the batch entirely.
        !
        ! Nodes 1:ncur are the current depth; children are appended past ncur and shifted down when the depth closes. Peak
        ! occupancy is ncur + 2*ncur <= 3*cap, inside the 4*cap + 8 the arrays already carry.
        do while (ncur > 0)
            nkeep = 0; nbat = 0; nbuf = 0; nnxt = 0
            ! pass 1: classify, and stash the signatures that need reducing. Rank-local nodes are NOT stashed (their
            ! signatures would swamp the buffer); they recompute in pass 2, which is one extra tag pass over a small box.
            do i = 1, ncur
                blo0 = slo(:,i); bhi0 = shi(:,i)
                ! S3.2b: a rank-local node's tags are ALL held by its one overlapping rank -- every other rank would contribute
                ! zeros, so the reduction cannot change the answer and the subtree is that rank's alone.
                ! S3.2b-2b: how many ranks the box spans, and whether THIS rank is one of them, both without enumerating the
                ! set -- the enumeration writes one entry per overlapping rank, which is O(P) on a box spanning the machine.
                novr = f_amr_overlap_count(blo0, bhi0)
                mine = (num_procs == 1) .or. f_amr_rank_overlaps(blo0, bhi0, proc_rank)
                ! counted BEFORE the drop, as the stack walk did: amr_cl_nodes/amr_cl_maxdep describe the TREE, which is the
                ! same tree whether or not this rank descends the parts it does not overlap, and [amr-tree] compares across runs
                nnode = nnode + 1_8; mxdep = max(mxdep, sdep(i))
                ! A rank holds tags only inside its own subdomain, so a node its subdomain does not reach is one it would
                ! contribute nothing but zeros to: drop the subtree and let the closing box ALLGATHERV carry back anything
                ! accepted inside it. S3.2b did this for novr == 1; the scoped exchange below extends it to every NARROW node.
                !
                ! WIDE nodes are deliberately NOT dropped, and that is load-bearing rather than conservative. They are settled by
                ! a collective over MPI_COMM_WORLD, which every rank must enter with the identical buffer length; if a rank
                ! skipped a wide node it did not overlap, its batch would be short and the reduction would mismatch -- a hang or
                ! silent corruption that appears only once some rank stops overlapping some wide box, i.e. only at scale. A wide
                ! node's ancestors are all wide (ovr only shrinks downward), so every rank reaches every wide node and the batch
                ! stays identical. A narrow node's members all walked its parent for the same reason, so p2p pairing is complete.
                if (reduce .and. num_procs > 1 .and. .not. mine .and. novr <= amr_cl_wide) cycle
                ! ONE pass over this node's tags yields the signature; trim, count and split all read it (no rescans).
                call s_amr_box_sig(wt, sts(i), ste(i), blo0, bhi0, sig, off, nsig)
                nkeep = nkeep + 1; kpos(nkeep) = i; kbat(nkeep) = 0
                amr_cl_rb = amr_cl_rb + int(nsig, 8)*4_8
                if (reduce) amr_cl_rb_now = amr_cl_rb_now + int(nsig, 8)*4_8
                if (novr > 1) then
                    amr_cl_shr_nodes = amr_cl_shr_nodes + 1_8; amr_cl_shr_rb = amr_cl_shr_rb + int(nsig, 8)*4_8
                    amr_cl_shr_maxdep = max(amr_cl_shr_maxdep, sdep(i))
                    if (reduce) then
                        amr_cl_shr_nodes_r = amr_cl_shr_nodes_r + 1_8; amr_cl_shr_rb_r = amr_cl_shr_rb_r + int(nsig, 8)*4_8
                        amr_cl_shr_maxdep_r = max(amr_cl_shr_maxdep_r, sdep(i))
                        ! S3.2a-2: under the sparse per-depth exchange this rank pays for a shared node only if the node's box
                        ! reaches into its subdomain. Everything else is somebody else's message.
                        if (mine) then
                            amr_cl_me_nodes_r = amr_cl_me_nodes_r + 1_8; amr_cl_me_rb_r = amr_cl_me_rb_r + int(nsig, 8)*4_8
                        end if
                    end if
                else
                    amr_cl_loc_nodes = amr_cl_loc_nodes + 1_8; amr_cl_loc_rb = amr_cl_loc_rb + int(nsig, 8)*4_8
                    if (reduce) then
                        amr_cl_loc_nodes_r = amr_cl_loc_nodes_r + 1_8; amr_cl_loc_rb_r = amr_cl_loc_rb_r + int(nsig, 8)*4_8
                    end if
                end if
                if (reduce .and. num_procs > 1 .and. novr > 1) then
                    call s_amr_fw_szi(bsig, nbuf + nsig)
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
            ! S3.2b-2b: the depth's reduction, SPLIT by how many ranks a node's box actually spans.
            !
            ! WIDE nodes are the shallow ones near the root. Every rank overlaps them and genuinely needs the answer, so they
            ! ride ONE batched collective per depth (that is 2a). There are only O(log P) of them, and their volume is the
            ! domain extent, which grows as P^(1/3) under weak scaling -- not as P.
            !
            ! NARROW nodes are the deep ones straddling a rank seam, and they are where the O(P) growth in the shared set lives
            ! (measured per regrid: 11/23/47/91 shared at np8/16/32/64). Their overlap set is a small rank-coordinate brick, so
            ! they reduce POINT-TO-POINT among exactly those ranks: each member ships its contribution to ovr(1), which sums and
            ! ships the total back. Per-rank received volume then follows what a rank actually overlaps -- measured 17,280 /
            ! 27,840 / 41,600 B at np8/16/32 (1.61x, 1.49x per doubling and decelerating) against 26,920 / 59,240 / 127,080 B
            ! (2.20x, 2.15x) for what an ALLREDUCE hands every rank.
            !
            ! Both ends agree on message contents with NO negotiation: a rank's node list at a depth is a SUBSEQUENCE of the one
            ! globally-ordered tree walk (ovr_child is contained in ovr_parent, so a rank that needs a child necessarily walked
            ! its parent), and both sides enumerate nodes in ascending j -- so the nodes common to a pair appear in the same
            ! relative order on both sides, and one aggregated message per peer per phase matches unambiguously.
            nwb = 0
            do j = 1, nbat
                if (bwide(j)) nwb = nwb + blen(j)
            end do
            if (nwb > 0) then
                call s_amr_fw_szi(wbuf, nwb)
                o1 = 0
                do j = 1, nbat
                    if (.not. bwide(j)) cycle
                    wbuf(o1 + 1:o1 + blen(j)) = bsig(bofs(j) + 1:bofs(j) + blen(j)); o1 = o1 + blen(j)
                end do
                call MPI_ALLREDUCE(MPI_IN_PLACE, wbuf, nwb, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
                amr_cl_wire_r = amr_cl_wire_r + int(nwb, 8)*4_8  ! a collective hands the WHOLE buffer to every rank
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
                ! received: the members' contributions I sum as a root (phase A) + the totals sent back to me (phase B)
                amr_cl_wire_r = amr_cl_wire_r + int(nrcv, 8)*4_8 + int(nsnd, 8)*4_8
                call s_amr_fw_szi(sbuf, max(nsnd, 1)); call s_amr_fw_szi(rbuf, max(nrcv, 1))
                call s_amr_fw_szi(creq, 2*np2)
                ! phase A: every member ships its own contribution up to the node's root
                soff(1:np2) = sdsp2(1:np2)
                do j = 1, nbat
                    if (bwide(j) .or. bovr(1, j) == proc_rank) cycle
                    q = pidx(bovr(1, j))
                    sbuf(soff(q) + 1:soff(q) + blen(j)) = bsig(bofs(j) + 1:bofs(j) + blen(j)); soff(q) = soff(q) + blen(j)
                end do
                tagc = amr_tag_base(4) + int(mod(amr_mesh_epoch, 50_8))
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
                ! the machine-wide reduction produced -- the change is who is in the message, never the arithmetic.
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
                    ! pass 2 only ever partitions a node's OWN wt(:, ts:te) range, which is disjoint from every other node's.
                    call s_amr_box_sig(wt, ts, te, blo0, bhi0, sig, off, nsig)
                end if
                call s_amr_trim_from_sig(sig(1:nsig), off, blo0, bhi0, blo, bhi, ok, ntag)
                if (.not. ok) cycle
                vol = 1_8
                do d = 1, num_dims; vol = vol*int(bhi(d) - blo(d) + 1, 8); end do
                eff = real(ntag, wp)/real(max(vol, 1_8), wp)
                call s_amr_find_split_sig(sig(1:nsig), off, blo0, blo, bhi, sax, spos, ok)
                ! splitting now could overflow the amr_max_blocks cap. The level-order walk changes what is still pending when
                ! this is asked, so the term counts the rest of THIS depth plus the children queued so far; B0b keeps the
                ! bisection clear of the cap, so it stays inert (measured: the capped warning on 0 of 10 regrids).
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

        ! S3.0a: record tree shape BEFORE the merge, so nacc is still the BR leaf count (the log2 denominator). Two independent
        ! maxima -- the deepest call and the largest call -- because a single max cannot say whether a deep tree was also big.
        amr_cl_nodes = amr_cl_nodes + nnode
        if (mxdep > amr_cl_maxdep) then
            amr_cl_maxdep = mxdep; amr_cl_maxdep_leaf = nacc
        end if
        if (nacc > amr_cl_lmax) then
            amr_cl_lmax = nacc; amr_cl_ldepth = mxdep
        end if

#ifdef MFC_MPI
        ! S3.2b: with rank-local subtrees walked ONLY by their owner, each rank now holds just the boxes from the subtrees it
        ! owns. Union them once here. This is per-BOX global data, which the endstate permits, and it is tiny -- 6 ints per box
        ! against the ~13,825 per-cell signature reductions per regrid it replaces. Ranks contribute in rank order, which is not
        ! the order the serial traversal accepted them in; B1's canonical Morton sort immediately below is what makes the merged
        ! result independent of that, and is the reason B1 had to land first.
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
            amr_gb_box = amr_gb_box + int(ntot, 8)*6_8*4_8  ! every rank receives the WHOLE global box list
            ! the accepted-box array is sized to the cap; B0b keeps the bisection away from it, and a run that still reaches it
            ! truncates here exactly as the serial `if (nacc < cap)` guard did.
            nacc = min(ntot, cap)
            do i = 1, nacc
                alo(:,i) = gbx(1:3,i); ahi(:,i) = gbx(4:6,i)
            end do
            deallocate (gcnt, gdsp, sbx, gbx)
        end if
#endif

        ! B1: canonicalise the merge input. The merge below scans in list order and fuses the FIRST too-close pair, so its
        ! output is a function of the order boxes were ACCEPTED -- i.e. of the traversal. Sorting by Morton of lo makes it a
        ! function of the box SET alone, which is what lets a scoped clusterer (S3.2) complete local subtrees in parallel, in
        ! a different acceptance order, and still agree across ranks. Accepted boxes are disjoint, so their lo corners are
        ! distinct and the key is a total order under f_morton's 21 bits/dim -- the same bound the block partition assumes.
        ! The sort is stable, so even a key collision above that bound would only fall back to acceptance order, never split
        ! the ranks. Morton rather than lexicographic because it keeps spatial neighbours adjacent, so the merge fuses near
        ! pairs first and the fused bounding boxes stay compact.
        do i = 1, nacc
            akey(i) = f_morton(alo(1, i), alo(2, i), alo(3, i))
        end do
        do i = 2, nacc
            bkey = akey(i); tmp = alo(:,i); tmp2 = ahi(:,i)
            j = i - 1
            do while (j >= 1)
                if (akey(j) <= bkey) exit
                akey(j + 1) = akey(j); alo(:,j + 1) = alo(:,j); ahi(:,j + 1) = ahi(:,j)
                j = j - 1
            end do
            akey(j + 1) = bkey; alo(:,j + 1) = tmp; ahi(:,j + 1) = tmp2
        end do

        ! min-separation merge: two boxes are separated only if some active dim's gap reaches thr; else fuse to their bounding box
        thr = buff_size + 2*amr_buf
        changed = .true.
        do while (changed)
            changed = .false.
            outer: do i = 1, nacc - 1
                do j = i + 1, nacc
                    tooclose = .true.
                    do d = 1, num_dims
                        if (max(alo(d, i), alo(d, j)) - min(ahi(d, i), ahi(d, j)) - 1 >= thr) tooclose = .false.
                    end do
                    if (tooclose) then
                        alo(:,i) = min(alo(:,i), alo(:,j)); ahi(:,i) = max(ahi(:,i), ahi(:,j))
                        ! B1: remove by shifting down, not by swapping with the last entry, so the list stays in the
                        ! canonical order the sort established and later fusions keep preferring spatial neighbours.
                        do k = j, nacc - 1
                            alo(:,k) = alo(:,k + 1); ahi(:,k) = ahi(:,k + 1)
                        end do
                        nacc = nacc - 1; changed = .true.
                        exit outer
                    end if
                end do
            end do outer
        end do
        if (capped .and. proc_rank == 0) print '(A,I0)', ' [amr] WARNING: tag clustering capped at amr_max_blocks = ', cap

        nboxes = nacc
        allocate (boxes(nboxes))
        do i = 1, nboxes
            boxes(i)%lo = alo(:,i); boxes(i)%hi = ahi(:,i)
        end do
        deallocate (slo, shi, alo, ahi, sts, ste, wt, sdep, sig, ovr, akey)

    end subroutine s_amr_cluster

    !> Regrid: tag by relative density gradient, cluster (Berger-Rigoutsos + min-separation merge) into separated boxes, pad/clamp/
    !! size-cap each, rebuild every active slot. Each new slot prolongs from coarse then overwrites its overlap with whichever OLD
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
        integer(8)           :: me_l(3), me_g(3)  !< S3.2a-2: this rank's shallow-phase participation, and its max over ranks

#ifdef MFC_MPI
        integer :: mierr
#endif

        allocate (box_level(amr_max_fine), old_ilo(3, amr_max_blocks), old_ext(3, amr_max_blocks), old_level(amr_max_blocks), &
                  & old_owns(amr_max_blocks))

        ! valid coarse CONS ghosts at internal rank boundaries: the tag sweep reads +/-1 across seams and the rebuild prolongation
        ! reads past the new intersection (ALL ranks call: pairwise per-direction exchange; complete no-op at np=1).

        call s_phase_tic(PH_RGHALO); call s_amr_exchange_coarse_cons_halo(q_cons_base); call s_phase_toc(PH_RGHALO)
        do i = 1, sys_size
            $:GPU_UPDATE(host='[q_cons_base(i)%sf]')
        end do

        ! Lagrangian-cloud exclusion bbox for this regrid (collective): smearing (mapCells) + stencil headroom (2) + drift
        ! margin until the next regrid (amr_buf)
        if (bubbles_lagrange) call s_amr_compute_lag_supp(mapCells + 2 + amr_buf)

        call s_phase_tic(PH_RGTAG); call s_amr_regrid_tag_cells(q_cons_base, tag_grid, sidx); call s_phase_toc(PH_RGTAG)
        call s_amr_cad_count(tag_grid, sidx)  ! [amr-cad] cadence containment audit (counts only; report at finalize)
        call s_phase_tic(PH_RGCLUS); call s_amr_regrid_cluster_tags(tag_grid, sidx, boxes, nboxes); call s_phase_toc(PH_RGCLUS)
        if (nboxes == 0) return  ! nothing tagged on any rank; keep the current blocks
        call s_phase_tic(PH_RGSHAPE); call s_amr_regrid_shape_boxes(boxes, nboxes); call s_phase_toc(PH_RGSHAPE)
        if (nboxes == 0) return  ! every box was confined to the domain margin
        call s_amr_regrid_nest_children(boxes, nboxes, box_level)
        call s_amr_check_box_caps(boxes, nboxes, box_level)  ! invariant: no box may exceed its level's slot cap
        call s_amr_check_box_disjoint(boxes, nboxes, box_level)  ! invariant: same-level boxes are pairwise disjoint
        call s_amr_regrid_boxes_unchanged(boxes, nboxes, box_level, same)
        if (same) return  ! identical box set and levels: keep the live slots
        call s_phase_tic(PH_RGMIG); call s_amr_regrid_stash_migrate(boxes, nboxes, box_level, old_np, old_ilo, old_ext, &
                         & old_level, old_owns); call s_phase_toc(PH_RGMIG)
        call s_phase_tic(PH_RGBUILD); call s_amr_regrid_rebuild_slots(q_cons_base, boxes, nboxes, old_np, old_ilo, old_ext, &
                         & old_level, old_owns); call s_phase_toc(PH_RGBUILD)

        ! TRACK S: the quantities that must stay O(1) in problem size. Reported per regrid on rank 0
        ! because wall time at one problem size cannot see them.
        ! S3.2a-2: rank_time_wrt is a namelist flag, so this branch is entered by every rank and the reduction is safe here.
        ! MAX rather than rank 0's own value: rank 0 owns a domain CORNER and overlaps the fewest shared boxes of anyone.
        if (rank_time_wrt) then
            me_l = [amr_cl_me_nodes_r, amr_cl_me_rb_r, amr_cl_wire_r]
            me_g = me_l
#ifdef MFC_MPI
            call MPI_ALLREDUCE(me_l, me_g, 3, MPI_INTEGER8, MPI_MAX, MPI_COMM_WORLD, mierr)
#endif
        end if
        if (rank_time_wrt .and. proc_rank == 0) then
            print '(A,I0,A,I0,A,I0,A,I0,A,I0)', '[amr-scope-me] me_nodes_max ', me_g(1), ' me_rb_max ', me_g(2), ' wire_max ', &
                & me_g(3), ' shr_nodes_all ', amr_cl_shr_nodes_r, ' shr_rb_all ', amr_cl_shr_rb_r
            print '(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)', '[amr-scale] nboxes ', nboxes, ' ntag_bytes ', amr_gb_tag, ' gwin_bytes ', &
                & amr_gb_win, ' cost_bytes ', amr_gb_cost, ' box_bytes ', amr_gb_box, ' cells ', int(m_glb + 1, 8)*int(n_glb + 1, &
                & 8)*int(p_glb + 1, 8)  ! int8: int32 overflows past ~1290^3
            print '(A,I0,A,I0,A,I0)', '[amr-mig] blocks_moved ', amr_mig_blk, ' sends ', amr_mig_snd, ' bytes ', amr_gb_mig
            print '(A,I0,A,I0,A,I0,A,I0,A,I0)', '[amr-scope] shr_nodes ', amr_cl_shr_nodes, ' shr_rb ', amr_cl_shr_rb, &
                & ' loc_nodes ', amr_cl_loc_nodes, ' loc_rb ', amr_cl_loc_rb, ' shr_maxdep ', amr_cl_shr_maxdep
            print '(A,I0,A,I0,A,I0,A,I0,A,I0)', '[amr-scope-r] shr_nodes ', amr_cl_shr_nodes_r, ' shr_rb ', amr_cl_shr_rb_r, &
                & ' loc_nodes ', amr_cl_loc_nodes_r, ' loc_rb ', amr_cl_loc_rb_r, ' shr_maxdep ', amr_cl_shr_maxdep_r
            print '(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)', '[amr-tree] maxdep ', amr_cl_maxdep, ' maxdep_leaf ', amr_cl_maxdep_leaf, &
                & ' lmax ', amr_cl_lmax, ' ldepth ', amr_cl_ldepth, ' nodes ', amr_cl_nodes, ' rbytes ', amr_cl_rb
        end if

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

        sidx = 0
        sidx(1) = start_idx(1)
        if (n_glb > 0) sidx(2) = start_idx(2)
        if (p_glb > 0) sidx(3) = start_idx(3)
        tg_lo = 0; tg_hi = 0
        tg_lo(1) = merge(1, 0, sidx(1) == 0); tg_hi(1) = merge(m - 1, m, sidx(1) + m == m_glb)
        if (n_glb > 0) then; tg_lo(2) = merge(1, 0, sidx(2) == 0); tg_hi(2) = merge(n - 1, n, sidx(2) + n == n_glb); end if
        if (p_glb > 0) then; tg_lo(3) = merge(1, 0, sidx(3) == 0); tg_hi(3) = merge(p - 1, p, sidx(3) + p == p_glb); end if
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
                    ! 2*r0 normalizes the 2-cell central difference (rho at i+1..i-1); the 2 is the stencil span, NOT amr_ref_ratio
                    if (g/(2._wp*r0) > amr_tag_eps) tag_grid(ci, cj, ck) = .true.
                    ! the acoustic source support stays coarse (its spatials are coarse cell indices): suppress tags there so
                    ! the clusterer splits around the source
                    if (acoustic_source .and. tag_grid(ci, cj, ck)) then
                        if (f_in_acoustic_support(ci + sidx(1), cj + sidx(2), ck + sidx(3))) tag_grid(ci, cj, ck) = .false.
                    end if
                    ! the Lagrangian bubble cloud stays coarse (two-way coupling lives on the coarse grid): suppress tags over
                    ! its padded bbox
                    if (bubbles_lagrange .and. tag_grid(ci, cj, ck)) then
                        if (f_in_lag_support(ci + sidx(1), cj + sidx(2), ck + sidx(3))) tag_grid(ci, cj, ck) = .false.
                    end if
                end do
            end do
        end do

    end subroutine s_amr_regrid_tag_cells

    !> [amr-cad] cadence containment audit: count this rank's level-1 tags, and how many fall OUTSIDE the pre-regrid level-1
    !! coverage - a feature that evolved unrefined since the last regrid because amr_buf did not cover its drift over amr_regrid_int
    !! steps. Counts only (reported once by s_amr_cov_report); the first refinement from scratch is skipped (every tag is new by
    !! construction). Region bounds are GLOBAL coarse cells; tag_grid is rank-local, so paint each region's clip with this subdomain
    !! into a local mask first.
    impure subroutine s_amr_cad_count(tag_grid, sidx)

        logical, intent(in)  :: tag_grid(0:,0:,0:)
        integer, intent(in)  :: sidx(3)
        logical, allocatable :: cov(:,:,:)
        integer              :: k, ci, cj, ck, bl(3), bh(3)
        logical              :: any_l1

        ! skip the FIRST regrid: it populates the hierarchy from the seed block, so nearly every tag is
        ! legitimately outside the old coverage (measured 47% escaped on S0 from the t=0 transient alone).
        ! The instrument measures STEADY-STATE containment - regrid 2 onward.

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
            bl = 0; bh = 0
            bl(1) = max(amr_region_lo_all(1, k) - sidx(1), 0); bh(1) = min(amr_region_hi_all(1, k) - sidx(1), m)
            if (n_glb > 0) then
                bl(2) = max(amr_region_lo_all(2, k) - sidx(2), 0); bh(2) = min(amr_region_hi_all(2, k) - sidx(2), n)
            end if
            if (p_glb > 0) then
                bl(3) = max(amr_region_lo_all(3, k) - sidx(3), 0); bh(3) = min(amr_region_hi_all(3, k) - sidx(3), p)
            end if
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

    !> Regrid phase 2: build this rank's OWN sparse tag list, then cluster it with per-node signature reductions into a list of
    !! separated candidate boxes (nboxes = 0 if nothing is tagged on any rank).
    impure subroutine s_amr_regrid_cluster_tags(tag_grid, sidx, boxes, nboxes)

        logical, allocatable, intent(inout)   :: tag_grid(:,:,:)
        integer, intent(in)                   :: sidx(3)
        type(t_box), allocatable, intent(out) :: boxes(:)
        integer, intent(out)                  :: nboxes
        integer, allocatable                  :: tags(:,:)
        integer                               :: ntag

        ! 2) build this rank's LOCAL tag list, then cluster it; the tree is driven by reduced signatures, so it stays rank-invariant

        call s_amr_local_tags(tag_grid, sidx, tags, ntag)
        deallocate (tag_grid)
        call s_amr_cluster(tags, ntag, boxes, nboxes, .true.)
        deallocate (tags)

    end subroutine s_amr_regrid_cluster_tags

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
            ! IB keeps the size-cap CLAMP (a body needs one contiguous block; splitting a body across tiles is untested); the
            ! general path leaves boxes full-size and TILES them (below) into <= amr_maxc_fit sub-blocks with a fine-fine halo
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
            ! only shrinks, so boxes stay disjoint - empties drop below
            if (acoustic_source) call s_amr_clip_box_from_sources(lo, hi)
            if (bubbles_lagrange .and. lag_supp_on) call s_amr_clip_box_from_supp(lo, hi, lag_supp_lo, lag_supp_hi)
            ! active_box: boxes stay strictly inside the active window (the windowed coarse update would drop reflux corrections
            ! at faces outside it). Tags cannot arise outside (frozen-ambient exterior), so only the amr_buf padding is ever cut
            ! - and the cut cells are ambient. np=1 only (ab_active is false under MPI).
            if (ab_active) then
                lo(1) = max(lo(1), ab_x%beg + 1); hi(1) = min(hi(1), ab_x%end - 1)
                if (n_glb > 0) then; lo(2) = max(lo(2), ab_y%beg + 1); hi(2) = min(hi(2), ab_y%end - 1); end if
                if (p_glb > 0) then; lo(3) = max(lo(3), ab_z%beg + 1); hi(3) = min(hi(3), ab_z%end - 1); end if
            end if
            ! a fine block that PARTIALLY covers an immersed body is an untested regime (ghost prolongation through body-interior
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
        ! (s_amr_fine_fine_halo) makes the seams conservative and the reflux skips fine-fine faces. (IB keeps the clamp - above.)
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
            ! double-restrict/reflux; a 1-cell gap with transverse overlap gives the two blocks a COINCIDENT outside coarse
            ! cell, which the batched reflux apply writes from both blocks in ONE kernel - an unsynchronized read-modify-write
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
            ! (contain the body, exclude the source/cloud) cannot both hold - fail closed
            if (acoustic_source .or. (bubbles_lagrange .and. lag_supp_on)) then
                do k = 1, nboxes
                    lo = boxes(k)%lo; hi = boxes(k)%hi
                    if (acoustic_source) call s_amr_clip_box_from_sources(lo, hi)
                    if (bubbles_lagrange .and. lag_supp_on) call s_amr_clip_box_from_supp(lo, hi, lag_supp_lo, lag_supp_hi)
                    if (ab_active) then
                        lo(1) = max(lo(1), ab_x%beg + 1); hi(1) = min(hi(1), ab_x%end - 1)
                        if (n_glb > 0) then; lo(2) = max(lo(2), ab_y%beg + 1); hi(2) = min(hi(2), ab_y%end - 1); end if
                        if (p_glb > 0) then; lo(3) = max(lo(3), ab_z%beg + 1); hi(3) = min(hi(3), ab_z%end - 1); end if
                    end if
                    if (any(lo /= boxes(k)%lo) .or. any(hi /= boxes(k)%hi)) then
                        call s_mpi_abort('amr regrid: a block must contain an immersed body AND stay ' &
                                         & // 'clear of an acoustic source support / Lagrangian bubble cloud - the ' &
                                         & // 'constraints conflict; move the body, source, or cloud apart')
                    end if
                end do
            end if
        end if

    end subroutine s_amr_regrid_shape_boxes

    !> Regrid phase 3b: multi-level nesting - hierarchically append level-l child boxes (sensor-on-fine, parents-first) inside each
    !! level-(l-1) box, for l = 2..amr_max_level. Sets box_level for every box (1 for the L0->L1 boxes).
    impure subroutine s_amr_regrid_nest_children(boxes, nboxes, box_level)

        type(t_box), allocatable, intent(inout) :: boxes(:)
        integer, intent(inout)                  :: nboxes
        integer, intent(inout)                  :: box_level(:)
        integer                                 :: i

        ! 3b) multi-level nesting: hierarchically append a level-l box nested inside each level-(l-1) box, for l = 2..amr_max_level.
        ! Parents-first ordering (every level-(l-1) box precedes its level-l children) so the build loop fills a parent before its
        ! child's gather-from-parent reads it. SENSOR-ON-FINE: each child's extent is the density-gradient sensor run on the
        ! parent-level FINE solution (the still-live OLD level-(l-1) blocks, read here BEFORE the step-5 stash), coarsened to
        ! L0-cell
        ! granularity and clustered - children track features inside the parent, not a fixed centre. A brand-new region with no old
        ! fine data falls back to a centred inset (sensor takes over next regrid); a parent with a smooth fine solution gets no
        ! child.
        ! Tagging only places boxes - conservation (restrict/reflux) is independent of where they sit. np=1 + non-IB (multi-level
        ! distribution / IB nesting are future work). Regions stay in L0 cell indices.

        box_level(1:nboxes) = 1
        if (amr_max_level >= 2) then
            ! the nesting loop below APPENDS level-l child boxes into `boxes` (up to amr_max_blocks total). The non-IB path already
            ! grew `boxes` to amr_max_blocks via the tiling move_alloc; the IB path (only merges, never grows) leaves `boxes` at the
            ! cluster count, so grow it here or the child appends overrun the allocation.
            if (size(boxes) < amr_max_blocks) then
                block
                    type(t_box), allocatable :: grown(:)
                    allocate (grown(amr_max_blocks))
                    grown(1:nboxes) = boxes(1:nboxes)
                    call move_alloc(grown, boxes)
                end block
            end if
            block
                integer                 :: kb, ins(3), clo(3), chi(3), lev, plo, phi, newlo, ob, obi, ncb, kc, mlo(3), mhi(3)
                integer                 :: mg, ng, pg, nct, np_lev, nloc_send, gi, gj, gk, ntot_g
                integer(8)              :: jrem  !< decode remainder spans an xy plane, which can exceed 2**31 cells
                integer, allocatable    :: ctags(:,:), skb(:), gkb(:)
                integer(8), allocatable :: sidx(:), gidx(:)
                logical, allocatable    :: gwin(:,:,:), covered(:)
                !> S3.3c: does THIS rank hold any level-(lev-1) block overlapping parent kb? Only then does it need the parent's
                !! dense window at all. `covered` stays REPLICATED (every rank must agree) and is recovered with ONE LOR reduction
                !! over the parents: the union over ranks of "my blocks overlapping kb" is exactly "all blocks overlapping kb",
                !! since each block has exactly one owner.
                logical, allocatable     :: mine(:)
                integer, allocatable     :: mlo_all(:,:), mhi_all(:,:)
                logical                  :: any_tag
                type(t_box), allocatable :: cboxes(:)
                !> S3.3a: one rank CLUSTERS each parent's window instead of every rank clustering every parent. powner is the
                !! assignment (replicated, computed with no communication); mych_* is this rank's own children, emitted without the
                !! global slot cap because a rank cannot see the global count mid-pass; gch_* is the assembled global list. The
                !! children are replayed into `boxes` in (kb, emission) order afterwards, which is exactly the order the serial loop
                !! produced them in -- so the box list, its truncation at amr_max_fine, and box_level are unchanged.
                integer, allocatable :: powner(:)
                integer, allocatable :: mych(:,:)            !< (7, n): lo(3), hi(3), kb
                integer, allocatable :: gch(:,:)
                integer              :: nmych, ntot_ch, ich, jch
                integer, allocatable :: chhead(:), chord(:)  !< counting-sort scratch for the canonical child order
                !> S3.3b: the pair exchange is TARGETED -- a rank sends each parent's tags only to that parent's owner, so send
                !! volume is O(this rank's tagged cells) and receive volume O(its assigned parents' tags), instead of every rank
                !! receiving every rank's pairs (the O(P) `gwin_bytes` term).
                integer, allocatable    :: phead(:), pord(:)
                integer(8), allocatable :: tidx(:)
                integer, allocatable    :: tkb(:)
#ifdef MFC_MPI
                integer              :: ierr, ip
                integer, allocatable :: rcnt(:), rdsp(:), scnt(:), sdsp(:)
#endif

                ! host-refresh the live (old) blocks' conserved state: the fine sensor below reads the flat store on the host,
                ! but the step-5 stash's GPU_UPDATE(host) runs AFTER this nesting - so the host copy is stale here
                do ob = 1, amr_num_blocks
                    if (amr_block_level(ob) == 0) cycle  ! L0 tiles are not regrid-managed and carry no fine sensor
                    if (.not. amr_owns_all(ob)) cycle  ! np>1: only the owner holds this old block's fine state
                    $:GPU_UPDATE(host='[amr_cons_st(:, :, :, :, amr_loc_of(ob))]')
                end do
                ! Fine-sensor tags accumulate in a GLOBAL L0 frame: at np>1 an old block is read only by its owner, but its tag
                ! footprint can fall in ANOTHER rank's subdomain. Each parent's nesting window [mlo:mhi] is small vs the global
                ! grid,
                ! so a WINDOW-LOCAL dense field gwin (per parent, below) holds each owner's tags; s_amr_pack_gwin_pairs extracts
                ! them
                ! as (linear-index, kb) pairs, one per-level allgatherv unions all parents' pairs across ranks, and pass 2 rebuilds
                ! each parent's window from them (no O(global-grid) tag field, no local slice; the clusterer consumes the sparse
                ! per-parent list directly).
                mg = m_glb; ng = 0; pg = 0
                if (n_glb > 0) ng = n_glb
                if (p_glb > 0) pg = p_glb

                plo = 1; phi = nboxes  ! [plo:phi] = the boxes at the previous level (lev-1) to nest inside
                do lev = 2, amr_max_level
                    newlo = nboxes + 1
                    ! COLLECT -> ONE COMMUNICATE -> PROCESS, per level: the per-parent cross-rank union is batched into a SINGLE
                    ! allgatherv per level, so the collective count is O(#levels) not O(#parent-boxes). Pass 1 tags each parent's
                    ! window from OWNED obs and appends this rank's tagged cells as (linear-index, parent-kb) pairs; one allgatherv
                    ! unions them; Pass 2 rebuilds each parent's dense window from the gathered pairs whose gkb==kb, reproducing the
                    ! old dense-window dedup and the (k,j,i) extraction order exactly -> each parent's ctags set (and thus its child
                    ! boxes) is byte-identical.
                    np_lev = phi - plo + 1
                    if (np_lev < 1) exit  ! nothing nested at the previous level -> no deeper levels possible
                    allocate (covered(plo:phi), mlo_all(3,plo:phi), mhi_all(3,plo:phi), mine(plo:phi))
                    covered = .false.; mine = .false.
                    ! S3.3c: ONE pass over this rank's OWNED level-(lev-1) blocks, not `do ob = 1, amr_num_blocks` inside
                    ! `do kb = plo, phi`. That nested form was O(parents x GLOBAL blocks) on every rank -- both factors scale
                    ! with P, so it was the O(P^2) term in the regrid. Here the outer loop is O(local blocks) and `covered`,
                    ! which must stay replicated, is recovered with a single LOR over the parents.
                    call s_amr_refresh_my_blocks()
                    do obi = 1, amr_n_my
                        ob = amr_my_blk(obi)
                        if (amr_block_level(ob) /= lev - 1) cycle
                        do kb = plo, phi
                            if (boxes(kb)%lo(1) > amr_region_hi_all(1, ob) .or. boxes(kb)%hi(1) < amr_region_lo_all(1, ob)) cycle
                            if (n_glb > 0) then
                                if (boxes(kb)%lo(2) > amr_region_hi_all(2, ob) .or. boxes(kb)%hi(2) < amr_region_lo_all(2, &
                                    & ob)) cycle
                            end if
                            if (p_glb > 0) then
                                if (boxes(kb)%lo(3) > amr_region_hi_all(3, ob) .or. boxes(kb)%hi(3) < amr_region_lo_all(3, &
                                    & ob)) cycle
                            end if
                            mine(kb) = .true.; covered(kb) = .true.
                        end do
                    end do
#ifdef MFC_MPI
                    if (num_procs > 1) call MPI_ALLREDUCE(MPI_IN_PLACE, covered, np_lev, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
#endif
                    nloc_send = 0
                    ! S3.3a: assign each parent a clustering owner. Round-robin over kb both balances the parents across ranks
                    ! and is a pure function of kb and num_procs, so every rank agrees without communicating. Pass 2 below then
                    ! processes ONLY its own parents, which deletes the per-parent rescan of the whole gathered pair list
                    ! (O(parents x global tags) on EVERY rank) and the replicated clustering along with it.
                    allocate (powner(plo:phi), mych(7, amr_max_fine))
                    do kb = plo, phi
                        powner(kb) = mod(kb - plo, max(num_procs, 1))
                    end do
                    nmych = 0
                    ! Pass 1: collect (no comm)
                    do kb = plo, phi
                        ! nesting window: children keep an amr_cpat_mar margin from the parent boundary so their ghost
                        ! prolongation reads valid parent interior cells
                        mlo = boxes(kb)%lo; mhi = boxes(kb)%hi
                        mlo(1) = mlo(1) + amr_cpat_mar; mhi(1) = mhi(1) - amr_cpat_mar
                        if (n_glb > 0) then; mlo(2) = mlo(2) + amr_cpat_mar; mhi(2) = mhi(2) - amr_cpat_mar; end if
                        if (p_glb > 0) then; mlo(3) = mlo(3) + amr_cpat_mar; mhi(3) = mhi(3) - amr_cpat_mar; end if
                        mlo_all(:,kb) = mlo; mhi_all(:,kb) = mhi
                        if (mhi(1) < mlo(1)) cycle  ! too small to nest a child in x
                        if (n_glb > 0 .and. mhi(2) < mlo(2)) cycle
                        if (p_glb > 0 .and. mhi(3) < mlo(3)) cycle

                        ! sensor-on-fine: tag from this rank's OWNED level-(lev-1) blocks overlapping the parent window
                        ! (amr_block_level still holds the old levels here - it is reset to box_level at step 5b, below).
                        ! S3.3c: `covered` and `mine` are already known from the pre-pass above, so this no longer scans the
                        ! global block list, and a rank with nothing to contribute allocates no dense window at all -- that
                        ! allocate+zero was O(parents x window volume) of pure waste on every non-contributing rank.
                        ! The parent's OWNER still builds a window when IB is on, because the body tags are its to add.
                        if (.not. (mine(kb) .or. (ib .and. powner(kb) == proc_rank))) cycle
                        allocate (gwin(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3)))
                        gwin = .false.; any_tag = .false.
                        do obi = 1, amr_n_my
                            ob = amr_my_blk(obi)
                            if (amr_block_level(ob) /= lev - 1) cycle
                            if (boxes(kb)%lo(1) > amr_region_hi_all(1, ob) .or. boxes(kb)%hi(1) < amr_region_lo_all(1, ob)) cycle
                            if (n_glb > 0) then
                                if (boxes(kb)%lo(2) > amr_region_hi_all(2, ob) .or. boxes(kb)%hi(2) < amr_region_lo_all(2, &
                                    & ob)) cycle
                            end if
                            if (p_glb > 0) then
                                if (boxes(kb)%lo(3) > amr_region_hi_all(3, ob) .or. boxes(kb)%hi(3) < amr_region_lo_all(3, &
                                    & ob)) cycle
                            end if
                            call s_amr_tag_child_from_fine(ob, mlo, mhi, gwin, any_tag)
                        end do
                        ! IB: always refine the body region at this level, even where the density sensor is quiet - mark the body's
                        ! L0-frame bbox into gwin so it is clustered into a child (mirrors the L1 expand at :3836). Containment
                        ! margin
                        ! = max(amr_buf, 4) + amr_cpat_mar: the child window (mlo:mhi) is the parent inset by amr_cpat_mar, and
                        ! clamping the tag to that window can eat up to amr_cpat_mar of the body's stencil margin at the
                        ! parent-adjacent side. The parent (widened in s_amr_expand_box_over_bodies by
                        ! (amr_max_level-1)*amr_cpat_mar)
                        ! now clears the body enough that this window contains the body plus max(amr_buf, 4), so the tag survives
                        ! the
                        ! inset with a full image-point stencil of fluid on every side: the body SURFACE is refined at every level
                        ! and
                        ! the C/F boundary sits a full stencil off it, in fluid.
                        if (ib) then
                            block
                                integer :: ib_i, bb_lo(3), bb_hi(3), gii, gjj, gkk
                                do ib_i = 1, num_ibs
                                    call s_amr_body_bbox(ib_i, max(amr_buf, 4) + amr_cpat_mar, bb_lo, bb_hi)
                                    ! clamp the body bbox to this parent's nesting window (s_amr_body_bbox returns GLOBAL L0
                                    ! cell indices, same frame as mlo/mhi)
                                    bb_lo = max(bb_lo, mlo); bb_hi = min(bb_hi, mhi)
                                    if (bb_hi(1) < bb_lo(1)) cycle
                                    if (n_glb > 0 .and. bb_hi(2) < bb_lo(2)) cycle
                                    if (p_glb > 0 .and. bb_hi(3) < bb_lo(3)) cycle
                                    covered(kb) = .true.
                                    do gkk = bb_lo(3), bb_hi(3)
                                        do gjj = bb_lo(2), bb_hi(2)
                                            do gii = bb_lo(1), bb_hi(1)
                                                gwin(gii, gjj, gkk) = .true.
                                            end do
                                        end do
                                    end do
                                end do
                            end block
                        end if
                        ! extract THIS rank's OWNED tagged cells as (linear-index, kb) pairs into the per-level send arrays. The
                        ! int8 linear index matches the pass-2 decode, so the gathered pairs reproduce the same window coords. gwin
                        ! is read, then freed.
                        call s_amr_pack_gwin_pairs(gwin, mlo, mhi, mg, ng, kb, sidx, skb, nloc_send)
                        deallocate (gwin)
                    end do

                    ! COMMUNICATE: one allgatherv per level (np>1)
                    if (.not. allocated(sidx)) then
                        allocate (sidx(0), skb(0))  ! this rank owned no tags at this level
                    end if
#ifdef MFC_MPI
                    if (num_procs > 1) then
                        ! S3.3b: bucket this rank's pairs by the DESTINATION owner (stable counting sort -- pass 1 appends in
                        ! kb order, and round-robin ownership interleaves the destinations), then exchange only what each rank
                        ! actually needs. Pass 2 rebuilds a DENSE window and re-extracts in (k,j,i) order, so ctags does not
                        ! depend on the order pairs arrive in and the box set is unchanged.
                        allocate (rcnt(num_procs), rdsp(num_procs), scnt(num_procs), sdsp(num_procs))
                        allocate (phead(num_procs), pord(max(nloc_send, 1)))
                        phead = 0
                        do i = 1, nloc_send
                            ip = powner(skb(i)) + 1; phead(ip) = phead(ip) + 1
                        end do
                        scnt = phead
                        sdsp(1) = 0
                        do ip = 2, num_procs
                            sdsp(ip) = sdsp(ip - 1) + scnt(ip - 1)
                        end do
                        phead = sdsp + 1
                        do i = 1, nloc_send
                            ip = powner(skb(i)) + 1; pord(phead(ip)) = i; phead(ip) = phead(ip) + 1
                        end do
                        allocate (tidx(max(nloc_send, 1)), tkb(max(nloc_send, 1)))
                        do i = 1, nloc_send
                            tidx(i) = sidx(pord(i)); tkb(i) = skb(pord(i))
                        end do
                        call MPI_ALLTOALL(scnt, 1, MPI_INTEGER, rcnt, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
                        rdsp(1) = 0
                        do ip = 2, num_procs
                            rdsp(ip) = rdsp(ip - 1) + rcnt(ip - 1)
                        end do
                        ntot_g = rdsp(num_procs) + rcnt(num_procs)
                        amr_gb_win = amr_gb_win + int(ntot_g, 8)*(8_8 + 8_8)  ! now the RECEIVED volume, i.e. O(local)
                        allocate (gidx(max(ntot_g, 1)), gkb(max(ntot_g, 1)))
                        call MPI_ALLTOALLV(tidx, scnt, sdsp, MPI_INTEGER8, gidx, rcnt, rdsp, MPI_INTEGER8, MPI_COMM_WORLD, ierr)
                        call MPI_ALLTOALLV(tkb, scnt, sdsp, MPI_INTEGER, gkb, rcnt, rdsp, MPI_INTEGER, MPI_COMM_WORLD, ierr)
                        deallocate (rcnt, rdsp, scnt, sdsp, phead, pord, tidx, tkb)
                    else
                        call move_alloc(sidx, gidx); call move_alloc(skb, gkb)
                        ntot_g = nloc_send
                    end if
#else
                    call move_alloc(sidx, gidx); call move_alloc(skb, gkb)
                    ntot_g = nloc_send
#endif
                    if (allocated(sidx)) deallocate (sidx)
                    if (allocated(skb)) deallocate (skb)

                    ! Pass 2: process (no comm). S3.3a: OWN PARENTS ONLY. The global `nboxes + 1 > amr_max_fine` guard cannot
                    ! be evaluated here any more -- a rank does not see the other ranks' children -- so children are emitted into
                    ! mych without a cap and the replay below applies the cap in the canonical order, which is the same order and
                    ! therefore the same truncation the serial loop performed.
                    do kb = plo, phi
                        if (powner(kb) /= proc_rank) cycle
                        if (nmych + 1 > amr_max_fine) exit  ! local buffer full (bounded by the same global cap)
                        mlo = mlo_all(:,kb); mhi = mhi_all(:,kb)
                        if (mhi(1) < mlo(1)) cycle  ! too small to nest a child in x
                        if (n_glb > 0 .and. mhi(2) < mlo(2)) cycle
                        if (p_glb > 0 .and. mhi(3) < mlo(3)) cycle

                        ! rebuild this parent's dense window from the gathered pairs whose gkb==kb: setting .true. once per gathered
                        ! cell reproduces the old per-parent dedup (replicated/overlapping tags collapse), and the (k,j,i) sparse
                        ! extract below matches the old scan order -> byte-identical ctags.
                        allocate (gwin(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3)))
                        gwin = .false.
                        do i = 1, ntot_g
                            if (gkb(i) /= kb) cycle
                            gk = int(gidx(i)/(int(mg + 1, 8)*int(ng + 1, 8)))
                            jrem = gidx(i) - int(gk, 8)*int(mg + 1, 8)*int(ng + 1, 8)
                            gj = int(jrem/int(mg + 1, 8))
                            gi = int(jrem - int(gj, 8)*int(mg + 1, 8))
                            gwin(gi, gj, gk) = .true.
                        end do
                        nct = 0
                        do gk = mlo(3), mhi(3); do gj = mlo(2), mhi(2); do gi = mlo(1), mhi(1)
                            if (gwin(gi, gj, gk)) nct = nct + 1
                        end do; end do; end do
                        allocate (ctags(3, max(nct, 1)))
                        nct = 0
                        do gk = mlo(3), mhi(3); do gj = mlo(2), mhi(2); do gi = mlo(1), mhi(1)
                            if (gwin(gi, gj, gk)) then
                                nct = nct + 1
                                ctags(1, nct) = gi; ctags(2, nct) = gj; ctags(3, nct) = gk
                            end if
                        end do; end do; end do
                        deallocate (gwin)
                        any_tag = nct > 0

                        ! smooth here - no child
                        if (covered(kb) .and. .not. any_tag) then; deallocate (ctags); cycle; end if

                        if (covered(kb)) then
                            ! cluster the fine-tagged L0 cells into child boxes, pad by amr_buf, clamp into the nesting window
                            call s_amr_cluster(ctags, nct, cboxes, ncb, .false.)  ! ctags is already replicated on every rank
                            deallocate (ctags)
                            do kc = 1, ncb
                                if (nboxes + 1 > amr_max_fine) exit
                                clo = cboxes(kc)%lo; chi = cboxes(kc)%hi
                                clo(1) = max(clo(1) - amr_buf, mlo(1)); chi(1) = min(chi(1) + amr_buf, mhi(1))
                                if (n_glb > 0) then
                                    clo(2) = max(clo(2) - amr_buf, mlo(2)); chi(2) = min(chi(2) + amr_buf, mhi(2))
                                else
                                    clo(2) = 0; chi(2) = 0
                                end if
                                if (p_glb > 0) then
                                    clo(3) = max(clo(3) - amr_buf, mlo(3)); chi(3) = min(chi(3) + amr_buf, mhi(3))
                                else
                                    clo(3) = 0; chi(3) = 0
                                end if
                                ! IB: a child clustered from the (widened) body tag must fully contain every overlapping body -
                                ! expand over bodies (mirrors the L1 expand at :3836), then re-clamp to the nesting window so the
                                ! child stays nested. Because the parent was widened by (amr_max_level-1)*amr_cpat_mar, its nesting
                                ! window (mlo:mhi) already contains the body plus max(amr_buf, 4), so the re-clamp does NOT cut the
                                ! body's stencil: the child CONTAINS the body bbox and the C/F boundary lands a full image-point
                                ! stencil off the surface, in fluid (surface refined, not just the interior).
                                if (ib) then
                                    call s_amr_expand_box_over_bodies(clo, chi)
                                    clo(1) = max(clo(1), mlo(1)); chi(1) = min(chi(1), mhi(1))
                                    if (n_glb > 0) then; clo(2) = max(clo(2), mlo(2)); chi(2) = min(chi(2), mhi(2)); end if
                                    if (p_glb > 0) then; clo(3) = max(clo(3), mlo(3)); chi(3) = min(chi(3), mhi(3)); end if
                                end if
                                ! slot cap: a level-lev block's fine grid spans amr_ref_ratio**lev*(its L0 extent) cells while the
                                ! slot holds amr_ref_ratio*amr_maxc_fit (max_f* = amr_ref_ratio*amr_maxc_fit - 1), so a child's L0
                                ! extent must be <= amr_maxc_fit/amr_ref_ratio**(lev-1) - HALVING once per level, not a fixed /2.
                                ! At lev = 2 that is amr_maxc_fit/2 (unchanged); at lev = 3 a fixed /2 admits a box twice what the
                                ! slot holds. VERIFIED to fail without this: m = 255, amr_max_level = 3, amr_buf = 48, np = 1 -
                                ! the fixed /2 keeps ONE oversized level-3 box where this keeps 2, and the run dies in
                                ! s_amr_free_slot ("Invalid descriptor", core dumped) on the corrupted field descriptor. It takes
                                ! all three of a big grid, depth 3 AND a wide buffer: boxes track the FEATURE, so scaling only the
                                ! grid leaves them far below either cap and both versions agree.
                                ! TILE a wider feature into adjacent sub-blocks (like the L1 tiling): the per-stage
                                ! fine-fine halo
                                ! (s_amr_fine_fine_halo, level-aware) matches the shared seam flux and the L2->L1 reflux skips those
                                ! fine-fine faces. Subcycle used to keep ONE capped child instead - under-refining a wide feature -
                                ! because s_amr_advance_children advanced children per-block with no L2-L2 halo; it now advances
                                ! siblings transposed with the level-filtered halo interposed, so both drivers tile alike.
                                block
                                    type(t_box) :: l2t(amr_max_blocks)
                                    integer     :: nl2, cpd, it
                                    nl2 = 0; cpd = 0
                                    call s_amr_tile_box(clo, chi, l2t, nl2, amr_max_blocks, cpd, &
                                                        & amr_maxc_fit/amr_ref_ratio**(lev - 1))
                                    do it = 1, nl2
                                        if (nmych + 1 > amr_max_fine) exit
                                        nmych = nmych + 1
                                        mych(1:3,nmych) = l2t(it)%lo; mych(4:6,nmych) = l2t(it)%hi; mych(7, nmych) = kb
                                    end do
                                end block
                            end do
                            if (allocated(cboxes)) deallocate (cboxes)
                        else
                            deallocate (ctags)  ! brand-new region: no fine tags to cluster
                            ! brand-new region (no old fine data yet): centred inset so the child still appears this regrid
                            ins = 0
                            ins(1) = max((boxes(kb)%hi(1) - boxes(kb)%lo(1) + 1)/4, amr_cpat_mar)
                            if (n_glb > 0) ins(2) = max((boxes(kb)%hi(2) - boxes(kb)%lo(2) + 1)/4, amr_cpat_mar)
                            if (p_glb > 0) ins(3) = max((boxes(kb)%hi(3) - boxes(kb)%lo(3) + 1)/4, amr_cpat_mar)
                            clo = boxes(kb)%lo + ins; chi = boxes(kb)%hi - ins
                            if (chi(1) < clo(1)) cycle  ! inset left no interior in x
                            if (n_glb > 0 .and. chi(2) < clo(2)) cycle
                            if (p_glb > 0 .and. chi(3) < clo(3)) cycle
                            ! Tile to the SAME slot cap as the clustered path above. The inset bounds the child as a FRACTION of
                            ! its parent, which is not the constraint that matters: the slot coord arrays are allocated once to
                            ! amr_ref_ratio*amr_maxc_fit, so the child must be bounded in ABSOLUTE cells. A parent of span 63
                            ! (an ordinary tile - s_amr_tile_box splits a wide region into 63 and 64, not 64 and 64) gives
                            ! ins = 63/4 = 15 and a child of span 33 against a level-2 cap of 32; s_amr_build_block_coords then
                            ! sizes fcb from the TRUE extent and writes one past x_cb. Span 64 gives exactly 32 and is fine, so a
                            ! one-cell difference in the parent flipped it.
                            block
                                type(t_box) :: nrt(amr_max_blocks)
                                integer     :: nnr, nrc, it2
                                nnr = 0; nrc = 0
                                call s_amr_tile_box(clo, chi, nrt, nnr, amr_max_blocks, nrc, amr_maxc_fit/amr_ref_ratio**(lev - 1))
                                do it2 = 1, nnr
                                    if (nmych + 1 > amr_max_fine) exit
                                    nmych = nmych + 1
                                    mych(1:3,nmych) = nrt(it2)%lo; mych(4:6,nmych) = nrt(it2)%hi; mych(7, nmych) = kb
                                end do
                            end block
                        end if
                    end do

                    ! S3.3a: assemble the children. Every parent has exactly ONE owner and that owner emitted its children in
                    ! order, so ONE allgatherv of the child BOXES (7 ints each: lo, hi, parent kb -- per-box global data, which
                    ! the endstate permits, ~67 KB against the 144 MB per-cell gather above) plus a STABLE sort by kb reproduces
                    ! the (kb ascending, emission) order the serial loop appended in. Box list, truncation at amr_max_fine and
                    ! box_level are therefore unchanged -- the gate for this increment is bit-identity.
#ifdef MFC_MPI
                    if (num_procs > 1) then
                        allocate (rcnt(num_procs), rdsp(num_procs))
                        call MPI_ALLGATHER(nmych, 1, MPI_INTEGER, rcnt, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
                        rdsp(1) = 0
                        do ip = 2, num_procs
                            rdsp(ip) = rdsp(ip - 1) + rcnt(ip - 1)
                        end do
                        ntot_ch = rdsp(num_procs) + rcnt(num_procs)
                        allocate (gch(7, max(ntot_ch, 1)))
                        rcnt = rcnt*7; rdsp = rdsp*7
                        call MPI_ALLGATHERV(mych, nmych*7, MPI_INTEGER, gch, rcnt, rdsp, MPI_INTEGER, MPI_COMM_WORLD, ierr)
                        amr_gb_box = amr_gb_box + int(ntot_ch, 8)*7_8*4_8  ! likewise: the whole global child list
                        deallocate (rcnt, rdsp)
                    else
                        ntot_ch = nmych
                        allocate (gch(7, max(ntot_ch, 1))); gch(:,1:ntot_ch) = mych(:,1:ntot_ch)
                    end if
#else
                    ntot_ch = nmych
                    allocate (gch(7, max(ntot_ch, 1))); gch(:,1:ntot_ch) = mych(:,1:ntot_ch)
#endif
                    ! stable counting sort by parent kb
                    allocate (chhead(plo:phi), chord(max(ntot_ch, 1)))
                    chhead = 0
                    do ich = 1, ntot_ch
                        chhead(gch(7, ich)) = chhead(gch(7, ich)) + 1
                    end do
                    jch = 1
                    do kb = plo, phi
                        ich = chhead(kb); chhead(kb) = jch; jch = jch + ich
                    end do
                    do ich = 1, ntot_ch
                        kb = gch(7, ich)
                        chord(chhead(kb)) = ich; chhead(kb) = chhead(kb) + 1
                    end do
                    do jch = 1, ntot_ch
                        ich = chord(jch)
                        if (nboxes + 1 > amr_max_fine) exit  ! pool full - stop nesting (same order, same truncation)
                        nboxes = nboxes + 1
                        boxes(nboxes)%lo = gch(1:3,ich); boxes(nboxes)%hi = gch(4:6,ich)
                        box_level(nboxes) = lev
                    end do
                    deallocate (gch, chhead, chord, powner, mych)

                    deallocate (gidx, gkb, covered, mine, mlo_all, mhi_all)  ! per-level scratch - freed every level
                    plo = newlo; phi = nboxes  ! the boxes just appended are the parents for the next level
                    if (phi < plo) exit  ! nothing nested at this level -> no deeper levels possible
                end do
                if (nboxes >= amr_max_fine .and. proc_rank == 0) print '(A)', &
                    & ' [amr] NOTE: block pool full during multi-level nesting; some boxes were not refined further'
            end block
        end if

    end subroutine s_amr_regrid_nest_children

    ! 4) unchanged? (same count, boxes AND levels as the live slots -> keep them; a rebuild would reproduce them exactly).
    ! The level must be compared too: a box that keeps its coordinates but changes refinement level would otherwise slip
    ! through with a stale amr_block_level, corrupting the level-aware coupling.
    impure subroutine s_amr_regrid_boxes_unchanged(boxes, nboxes, box_level, same)

        type(t_box), intent(in) :: boxes(:)
        integer, intent(in)     :: nboxes, box_level(:)
        logical, intent(out)    :: same
        integer                 :: k, ks

        ! regrid manages only the FINE band [l0_slot_off+1 ..] of the shared pool; the level-0 L0-tile prefix (coexist) is not
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

    !> DEVICE pack of an owned old block's stash into the migration wire buffer (wp wire, stp store): the store is
    !! device-authoritative during the rebuild, so pack where the data lives; the copyout hands the host MPI buffer back as one
    !! contiguous transfer of exactly the interior. Wire layout (gi fastest, then gj, gk, ii) matches the old host pack
    !! byte-for-byte, so the message set and [amr-xa] F4 totals are unchanged.
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

    !> DEVICE unpack of a received old block into its stash replica (mirror of the pack; the copyin stages the wire buffer).
    !! Replaces the host cast loop AND the full-slot device push it required.
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

    !> DEVICE cons->stor stash copy of one owned old block's fine interior (the store is device-authoritative; no host staging).
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

    !> DEVICE overlap carry-forward: overwrite the new block's prolonged cons (slot column vloc, extents vm/vn/vp) with the covering
    !! old block's stashed fine detail (column vold, extents ve*), shifted by sh. Guards mirror the old host loop.
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
        integer                 :: k, i, ks
        integer                 :: np_l  !< local mirror of old_np: an INTENT(OUT) dummy is not allowed in the
        !                                   BLOCK specification expressions below (F2018 restricted expressions)

        ! 5) stash every live slot's fine interior (dead-between-steps q_cons_stor bounce), keeping its old intersection origin

        ! old_* are indexed in the regrid's own dense FINE-block space [1..old_np], which maps to shared-pool slot f_l0_slot(k);
        ! under coexist the level-0 L0-tile prefix [1..l0_slot_off] is not regrid-managed and must not be stashed or migrated.

        call s_phase_tic(PH_RGPART)
        old_np = amr_num_blocks - l0_slot_off
        np_l = old_np
        do k = 1, old_np
            ks = f_l0_slot(k)
            ! GLOBAL block origin + extents (replicated, valid on every rank - not the owner-only isect), so the cross-rank
            ! migration below and the overlap-copy's index shift are correct even where this rank did not own the old block
            old_ilo(:,k) = amr_region_lo_all(:,ks)
            old_chi(:,k) = amr_region_hi_all(:,ks)  ! old COARSE hi (for the P2P migration overlap test below)
            ! fine extent = (2**level)*footprint - 1: a level-2 block is 4x its L0 footprint, so stashing/migrating it with the
            ! level-1 factor (2x) truncates half its fine cells. Level-1 blocks (2**1 = 2) are byte-identical to before.
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
                ! DEVICE-side stash: the store is device-authoritative, so copy cons->stor where the data lives instead of
                ! staging two full-slot transfers through the host (the old pull/host-copy/push). A mid-rebuild grow's
                ! device->host round trip preserves this stash by construction (s_amr_st_reserve's contract); the host
                ! mirror of amr_stor_st stays stale, which is fine - the migration pack and the overlap carry-forward
                ! below are device kernels now and no host reader of the stash remains. (The kernel lives in its own
                ! subroutine: amdflang drops target regions nested in BLOCK constructs from the device image - the host
                ! registers them and the first launch dies on HSA_STATUS_ERROR_INVALID_SYMBOL_NAME.)
                call s_amr_stash_copy_device(amr_loc_of(ks), old_ext(1, k), old_ext(2, k), old_ext(3, k))
                ! non-polytropic QBMM: the side-state bounces through pb/mv_stor exactly like q_cons (both stors are dead between
                ! steps)
                if (qbmm .and. .not. polytropic) then
                    $:GPU_UPDATE(host='[amr_slots(ks)%pb_f%sf, amr_slots(ks)%mv_f%sf]')
                    amr_slots(ks)%pb_stor%sf(0:old_ext(1, k),0:old_ext(2, k),0:old_ext(3, k),:, &
                              & :) = amr_slots(ks)%pb_f%sf(0:old_ext(1, k),0:old_ext(2, k),0:old_ext(3, k),:,:)
                    amr_slots(ks)%mv_stor%sf(0:old_ext(1, k),0:old_ext(2, k),0:old_ext(3, k),:, &
                              & :) = amr_slots(ks)%mv_f%sf(0:old_ext(1, k),0:old_ext(2, k),0:old_ext(3, k),:,:)
                end if
            end if
        end do
        ! coarse pb/mv host-current for the per-block re-prolongation below
        if (qbmm .and. .not. polytropic) then
            $:GPU_UPDATE(host='[pb_ts(1)%sf, mv_ts(1)%sf]')
        end if

        ! set the regions + assign owners BEFORE the migration (P2P needs the new owners) and before the owner-dependent
        ! geometry (else s_set_amr_fine_geometry sizes the whole-block owner from a stale amr_block_owner)
        ! the fine band ends at f_l0_slot(nboxes); the level-0 tile prefix below it keeps its regions, levels and owners (a plain
        ! nboxes here is what overran the prefix and deadlocked the first regrid under coexist)
        amr_num_blocks = f_l0_slot(nboxes)
        do k = 1, nboxes
            ks = f_l0_slot(k)
            amr_region_lo_all(:,ks) = boxes(k)%lo; amr_region_hi_all(:,ks) = boxes(k)%hi
            ! box_level(k) is the refinement level assigned during the hierarchical nesting above (1 for L0->L1 boxes, l for a
            ! box nested at level l). Setting it every regrid resets a stale level when a slot is reused across levels.
            amr_block_level(ks) = box_level(k)
        end do
        ! block set changed: dirty the cached seam-pair AND overlap-rank lists NOW - the rebuild's per-block P2P gathers
        ! (s_amr_regrid_rebuild_slots) consume the overlap lists with the NEW boxes, so flagging after them would be too late
        amr_seam_pairs_dirty = .true.
        amr_mesh_epoch = amr_mesh_epoch + 1
        ! Proper-nesting guard: each level>=2 block must be covered by EXACTLY ONE parent-level block. f_amr_parent_block (and
        ! the gather/reflux that key off it) take the FIRST overlap, so a fine tile straddling two parent tiles - an internal
        ! parent-level tile seam crossed by a nested feature - would silently couple to only one parent (wrong coarse BC + a
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
        ! The partition is decided here and NOTHING has moved yet; everything below redistributes data.
        call s_phase_toc(PH_RGPART)
        call s_phase_tic(PH_RGMOVE)

#ifdef MFC_MPI
        ! Cross-rank fine-state migration: the overlap-copy below preserves each covering old block's fine detail by reading
        ! amr_slots(kk)%q_cons_stor, but an old block may be owned by a rank OTHER than the one now owning a covering new block.
        ! POINT-TO-POINT (mirrors s_amr_gather_coarse_patch): each old owner sends its stashed fine state ONLY to the distinct
        ! new-block owners whose region overlaps that old block. A rank that did not receive old block kk never reads it - the
        ! overlap-copy's per-(k,kk) index guard skips every cell of a non-overlapping pair. No-op at np=1 (single owner, local).
        if (num_procs > 1) then
            block
                integer               :: kk, k2, ierr2, rr, nrq
                integer               :: cnt(np_l), scol(np_l), rcol(np_l)
                integer               :: nsnd, nrcv, nsreq, maxsnd, maxrcv
                logical               :: getk(np_l), isdest(0:num_procs - 1)
                real(wp), allocatable :: spack(:,:), rpack(:,:)
                integer, allocatable  :: rq(:)
                ! I4a right-sizing: pack/request pools sized to the blocks ACTUALLY sent/received, not old_np columns of the
                ! largest block each plus an O(old_np x ranks) request array - at production counts those allocated GBs per
                ! regrid for a handful of live columns. Message set, sizes, tags, and order are UNCHANGED (byte-exact gate:
                ! identical [amr-xa] F4 totals).
                nrcv = 0; maxrcv = 0
                do kk = 1, old_np
                    cnt(kk) = sys_size*(old_ext(1, kk) + 1)*(old_ext(2, kk) + 1)*(old_ext(3, kk) + 1)
                    ! I need old block kk iff I own a NEW block overlapping it (and do not already hold kk locally)
                    getk(kk) = .false.
                    rcol(kk) = 0
                    if (.not. old_owns(kk)) then
                        do k2 = 1, nboxes
                            if (amr_block_owner(f_l0_slot(k2)) == proc_rank .and. f_amr_boxes_overlap(boxes(k2)%lo, boxes(k2)%hi, &
                                & old_ilo(:,kk), old_chi(:,kk))) then
                                getk(kk) = .true.; exit
                            end if
                        end do
                    end if
                    if (getk(kk)) then
                        nrcv = nrcv + 1; rcol(kk) = nrcv; maxrcv = max(maxrcv, cnt(kk))
                    end if
                end do
                ! pre-pass over my owned blocks: which are sent anywhere, and how many sends in total (sizes rq exactly)
                nsnd = 0; nsreq = 0; maxsnd = 0
                do kk = 1, old_np
                    scol(kk) = 0
                    if (.not. old_owns(kk)) cycle
                    isdest = .false.
                    do k2 = 1, nboxes
                        rr = amr_block_owner(f_l0_slot(k2))
                        if (rr /= proc_rank .and. f_amr_boxes_overlap(boxes(k2)%lo, boxes(k2)%hi, old_ilo(:,kk), old_chi(:, &
                            & kk))) isdest(rr) = .true.
                    end do
                    if (.not. any(isdest)) cycle
                    nsnd = nsnd + 1; scol(kk) = nsnd; maxsnd = max(maxsnd, cnt(kk))
                    nsreq = nsreq + count(isdest)
                end do
                ! a received old block needs a live slot to unpack its q_cons_stor into (freed by the rebuild's early-free or
                ! the reconcile below) - STASH-ONLY: a replica never touches q_prim/rhs, and full slots across the np-scaled
                ! replica set of a migration-heavy regrid are what OOMed the W8 gate's np=4 arm
                call s_phase_tic(PH_MGSLOT)
                call s_amr_prereserve_stash(getk, old_np)
                do kk = 1, old_np
                    if (getk(kk)) call s_amr_alloc_slot_stash(f_l0_slot(kk))
                end do
                call s_phase_toc(PH_MGSLOT)
                allocate (rq(max(nsreq + nrcv, 1)), spack(max(maxsnd, 1), max(nsnd, 1)), rpack(max(maxrcv, 1), max(nrcv, 1)))
                nrq = 0
                do kk = 1, old_np  ! post receives for the old blocks I need
                    if (.not. getk(kk)) cycle
                    nrq = nrq + 1
                    call s_xa_rec(XA_F4_RCV, 2, cnt(kk), kk)
                    call MPI_IRECV(rpack(1, rcol(kk)), cnt(kk), mpi_p, old_owner(kk), kk, MPI_COMM_WORLD, rq(nrq), ierr2)
                end do
                call s_phase_tic(PH_MGPACK)
                do kk = 1, old_np  ! pack + send each old block I own to every distinct new-owner (/= me) overlapping it
                    if (scol(kk) == 0) cycle  ! not mine, or no remote destination (pre-pass above)
                    isdest = .false.
                    do k2 = 1, nboxes
                        rr = amr_block_owner(f_l0_slot(k2))
                        if (rr /= proc_rank .and. f_amr_boxes_overlap(boxes(k2)%lo, boxes(k2)%hi, old_ilo(:,kk), old_chi(:, &
                            & kk))) isdest(rr) = .true.
                    end do
                    amr_mig_blk = amr_mig_blk + 1_8
                    call s_amr_mig_pack_device(amr_loc_of(f_l0_slot(kk)), old_ext(1, kk), old_ext(2, kk), old_ext(3, kk), &
                                               & spack(1:cnt(kk),scol(kk)))
                    do rr = 0, num_procs - 1
                        if (.not. isdest(rr)) cycle
                        nrq = nrq + 1
                        amr_mig_snd = amr_mig_snd + 1_8
                        amr_gb_mig = amr_gb_mig + int(cnt(kk), 8)*8_8
                        call s_xa_rec(XA_F4_SND, 1, cnt(kk), kk)
                        call MPI_ISEND(spack(1, scol(kk)), cnt(kk), mpi_p, rr, kk, MPI_COMM_WORLD, rq(nrq), ierr2)
                    end do
                end do
                call s_phase_toc(PH_MGPACK)
                call s_phase_tic(PH_MGWAIT)
                if (nrq > 0) call MPI_WAITALL(nrq, rq, MPI_STATUSES_IGNORE, ierr2)
                call s_phase_toc(PH_MGWAIT)
                do kk = 1, old_np  ! unpack the received old blocks into their replicated q_cons_stor slots (DEVICE-direct:
                    ! the copyin stages exactly the packed interior, and the replica lands device-side where the store is
                    ! authoritative - no host cast loop and no full-slot push, and a mid-rebuild grow preserves it)
                    if (.not. getk(kk)) cycle
                    call s_phase_tic(PH_MGUNPK)
                    call s_amr_mig_unpack_device(amr_loc_of(f_l0_slot(kk)), old_ext(1, kk), old_ext(2, kk), old_ext(3, kk), &
                                                 & rpack(1:cnt(kk),rcol(kk)))
                    call s_phase_toc(PH_MGUNPK)
                end do
                deallocate (rq, spack, rpack)
            end block
        end if
#endif
        call s_phase_toc(PH_RGMOVE)  ! outside MFC_MPI so the bracket is balanced in serial builds

    end subroutine s_amr_regrid_stash_migrate

    !> Regrid phase 6: build each new slot - geometry (collective), prolong from coarse, overwrite the overlap from every covering
    !! stashed old block - then reconcile the slot pool, rebuild the fine IB state, and re-validate the seam topology.
    impure subroutine s_amr_regrid_rebuild_slots(q_cons_base, boxes, nboxes, old_np, old_ilo, old_ext, old_level, old_owns)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_base
        type(t_box), intent(in)                                :: boxes(:)
        integer, intent(in)                                    :: nboxes, old_np, old_ilo(:,:), old_ext(:,:), old_level(:)
        logical, intent(in)                                    :: old_owns(:)
        integer                                                :: sh(3), k, kk, i, fi, fj, fk, ofi, ofj, ofk, ks, kks, ohi(3)
        integer                                                :: c_lo, c_hi
        integer, allocatable                                   :: last_use(:)

        ! 6) build each new slot: geometry (collective on all ranks), prolong, then overlap-copy from every covering old slot
        ! box k lives in shared-pool slot ks = f_l0_slot(k) (identity without L0 tiles); old block kk in slot kks

        ! W8 transient: last_use(kk) = the last new box whose region overlaps old block kk, i.e. the last iteration whose
        ! overlap-copy (or pbmv copy) can read kk's stash. Holding EVERY stashed/received old block until the reconcile is
        ! what peaks device memory at np >= 2 - the replica count grows with np - so old-only slots are freed as soon as
        ! their last covering box is built (the loop below), and the freed dense indices recycle into the very next allocs.
        ! Region overlap (all regions are L0-cell coords at every level) is a superset of every per-cell stash read.

        allocate (last_use(old_np)); last_use = 0
        do kk = 1, old_np
            ohi = old_ilo(:,kk)
            ohi(1) = ohi(1) + (old_ext(1, kk) + 1)/amr_ref_ratio**old_level(kk) - 1
            if (n_glb > 0) ohi(2) = ohi(2) + (old_ext(2, kk) + 1)/amr_ref_ratio**old_level(kk) - 1
            if (p_glb > 0) ohi(3) = ohi(3) + (old_ext(3, kk) + 1)/amr_ref_ratio**old_level(kk) - 1
            do k = 1, nboxes
                if (f_amr_boxes_overlap(boxes(k)%lo, boxes(k)%hi, old_ilo(:,kk), ohi)) last_use(kk) = k
            end do
        end do

        ! gather-batching step 1: derive the whole loop's gather message set up front; step 2 executes the exchange from it
        call s_amr_build_gather_plan(nboxes)

        c_lo = 1
        do k = 1, nboxes
            ! gather-batching step 2 (amr_regrid_gather_batching.md): at each chunk boundary, pre-post the chunk's recvs and
            ! issue its sends (level>=2 sends whose parent shares the chunk are deferred to the child's consume below), so the
            ! per-box rendezvous becomes one wait per owned box against an exchange already in flight
            if (mod(k - 1, amr_gath_chunk) == 0) then
                c_lo = k; c_hi = min(k + amr_gath_chunk - 1, nboxes)
                call s_phase_tic(PH_RBGATH)
                call s_amr_gather_chunk_post(c_lo, c_hi)
                call s_amr_gather_chunk_send(q_cons_base, c_lo, c_hi)
                call s_phase_toc(PH_RBGATH)
            end if
            ! free the old-only slots consumed by iterations before this one (last_use = k - 1 fires exactly once; a slot
            ! serving as a NEW owned box keeps living - the reconcile decides it)
            do kk = 1, old_np
                if (last_use(kk) /= k - 1) cycle
                kks = f_l0_slot(kk)
                if (kks <= amr_num_blocks) then
                    if (amr_block_owner(kks) == proc_rank) cycle
                end if
                call s_amr_free_slot(kks)
            end do
            ks = f_l0_slot(k)
            amr_cur = ks
            ! owned slot needs its arrays before geometry/prolong
            call s_phase_tic(PH_RBSLOT)
            if (amr_block_owner(ks) == proc_rank) call s_amr_alloc_slot(ks)
            call s_phase_toc(PH_RBSLOT)
            call s_phase_tic(PH_RBGEO)
            call s_set_amr_fine_geometry(boxes(k)%lo, boxes(k)%hi)
            call s_phase_toc(PH_RBGEO)
            ! fine-level distribution: consume this new block's coarse patch out of the chunk exchange (all ranks - before the
            ! owner-only cycle, which is what lets the chunk request arrays be reused; q_cons_base is host-current with valid
            ! ghosts from the exchange at the top of s_amr_regrid)
            call s_phase_tic(PH_RBGATH); call s_amr_gather_consume_box(q_cons_base, k, c_lo); call s_phase_toc(PH_RBGATH)
            ! non-polytropic QBMM: gather the coarse pb/mv patch too (ALL ranks - P2P; owners re-prolong from it below)
            if (qbmm .and. .not. polytropic) call s_amr_gather_coarse_patch_pbmv(pb_ts(1)%sf, mv_ts(1)%sf, .false.)
            if (.not. amr_rank_owns_block) cycle
            call s_amr_cov_note(old_np, old_ilo, old_ext, old_level)  ! [amr-cov] rebuild-gather coverage split
            ! prolong and overlap carry-forward are both DEVICE kernels now: the slot is built entirely in place where the
            ! store is authoritative, and the per-box full-slot push (PH_RBPUSH) is gone.
            call s_phase_tic(PH_RBOVL); call s_interpolate_coarse_to_fine()
            call s_phase_toc(PH_RBOVL)
            ! every old block's stashed fine state is now replicated in amr_slots(kk)%q_cons_stor (migration above), so copy the
            ! overlap from EVERY covering old block regardless of owner - sh is the old->new LOCAL fine index shift. A level>=2
            ! block
            ! SKIPS this: old_ilo/sh are the L0 index frame, but a child's amr_isect_lo is its PARENT-fine frame, so the shift is
            ! wrong. It re-prolongs from its (freshly-built, parents-first) parent each regrid instead; the coupling keeps
            ! conservation. Detail-preserving same-level L2 migration (parent-fine overlap) is a later increment.
            call s_phase_tic(PH_RBOVL)
            if (amr_block_level(amr_cur) < 2) then
                do kk = 1, old_np
                    ! same-level overlap only (a child's stash is 4x-framed)
                    if (old_level(kk) /= amr_block_level(amr_cur)) cycle
                    kks = f_l0_slot(kk)
                    ! old LOCAL fine index = new LOCAL fine index + sh (collapsed dims sh=0)
                    sh = amr_ref_ratio*(amr_isect_lo - old_ilo(:,kk))
                    call s_amr_overlap_copy_device(amr_loc_of(ks), amr_loc_of(kks), amr_slots(ks)%m, amr_slots(ks)%n, &
                                                   & amr_slots(ks)%p, sh, old_ext(1, kk), old_ext(2, kk), old_ext(3, kk))
                end do
            end if
            call s_phase_toc(PH_RBOVL)
            ! non-polytropic QBMM: prolong the side-state from coarse (piecewise-constant), then overwrite the overlap with the
            ! old blocks' fine data (same index shift)
            if (qbmm .and. .not. polytropic) then
                call s_amr_prolong_pbmv()
                ! level>=2 re-prolongs only (the L0-frame overlap shift is wrong for a child)
                if (amr_block_level(amr_cur) < 2) then
                    do kk = 1, old_np
                        if (old_level(kk) /= amr_block_level(amr_cur)) cycle  ! same-level overlap only
                        if (.not. old_owns(kk)) cycle
                        kks = f_l0_slot(kk)
                        sh = amr_ref_ratio*(amr_isect_lo - old_ilo(:,kk))
                        do fk = 0, amr_slots(ks)%p
                            ofk = fk + sh(3)
                            if (p_glb > 0 .and. (ofk < 0 .or. ofk > old_ext(3, kk))) cycle
                            do fj = 0, amr_slots(ks)%n
                                ofj = fj + sh(2)
                                if (n_glb > 0 .and. (ofj < 0 .or. ofj > old_ext(2, kk))) cycle
                                do fi = 0, amr_slots(ks)%m
                                    ofi = fi + sh(1)
                                    if (ofi < 0 .or. ofi > old_ext(1, kk)) cycle
                                    amr_slots(ks)%pb_f%sf(fi, fj, fk,:,:) = amr_slots(kks)%pb_stor%sf(ofi, ofj, ofk,:,:)
                                    amr_slots(ks)%mv_f%sf(fi, fj, fk,:,:) = amr_slots(kks)%mv_stor%sf(ofi, ofj, ofk,:,:)
                                end do
                            end do
                        end do
                    end do
                end if
                $:GPU_UPDATE(device='[amr_slots(ks)%pb_f%sf, amr_slots(ks)%mv_f%sf]')
            end if
            ! whole-block-per-rank: no fine-fine halo; the new block's ghost shell is (re)prolonged by the next fine advance
        end do
        deallocate (last_use)
        amr_gpl_valid = .false.  ! the plan describes THIS rebuild's box loop only; per-step gathers never consult it

        ! Drain the deferred gather sends now that every box has been posted: one WAITALL per rebuild instead of
        ! a per-box rendezvous. MUST happen before the send buffers are reused or freed.
        call s_phase_tic(PH_RBTAIL)
        call s_phase_tic(PH_RBFLUSH); call s_amr_gather_send_flush(); call s_phase_toc(PH_RBFLUSH)
        ! ONE allreduce for the whole loop; sets amr_xchg_coarse_ghosts if ANY block needs it
        call s_phase_tic(PH_RBXCHG); call s_amr_reduce_xchg_flag(); call s_phase_toc(PH_RBXCHG)
        ! lazy sizing: free the transient regrid slots (old blocks this rank stashed/received but does not now own); the
        ! new-owned slots were allocated in the build loop, so this only frees - a rank keeps just its owned blocks' fine arrays
        call s_phase_tic(PH_RBREC); call s_amr_reconcile_slots(); call s_phase_toc(PH_RBREC)
        ! rebuild every block's fine-grid IB state for the NEW geometry (markers/ghost points/image points recomputed from the
        ! body definitions; no state carries across regrids)
        if (ib) call s_amr_setup_ib()
        call s_amr_select_slot(1)
        call s_phase_tic(PH_RBTOPO); call s_amr_check_seam_topology(); call s_phase_toc(PH_RBTOPO)
        call s_phase_toc(PH_RBTAIL)  ! abort on seam topologies no halo reconciles (silent leak otherwise)

    end subroutine s_amr_regrid_rebuild_slots

    !> Sensor-on-fine child tagging: OR-accumulate density-gradient tags from an OLD fine block's solution into an L0-cell tag grid,
    !! restricted to a parent nesting window. Reads amr_slots(ob)%q_cons on the HOST (caller host-refreshes the cont range first;
    !! the step-5 stash's GPU_UPDATE runs later). Fine cell (fi,fj,fk) covers L0 cell (ci,cj,ck) with fi = rr*(ci-olo(1))+d etc.;
    !! the gradient uses one-sided differences at the fine-interior edges so no stale fine ghost is read. Only decides placement -
    !! conservation is enforced downstream by restrict/reflux regardless of box extent.
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
                                ! 2*r0 normalizes the 2-cell central difference; the 2 is the stencil span, NOT amr_ref_ratio
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
