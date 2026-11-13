!>
!!@file
!!@brief Contains module m_amr_regrid

#:include 'macros.fpp'

!> @brief Dynamic regrid for the AMR level set: density-gradient tagging, clustering (m_amr_cluster), box shaping, hierarchical
!! child nesting (m_amr_nest), and slot rebuild with cross-rank fine-state migration.
module m_amr_regrid

#ifdef MFC_MPI
    use mpi
#endif

    use m_derived_types
    use m_global_parameters
    use m_constants, only: mapCells
    use m_mpi_proxy, only: s_mpi_abort
    use m_mpi_common, only: s_mpi_allreduce_min, s_mpi_allreduce_max
    use m_amr_wave
    use m_amr_state, only: amr_band_migrate
    use m_phase_timing
    use m_amr_cluster, only: s_amr_cluster
    use m_amr_nest, only: s_amr_regrid_nest_children, s_amr_nest_window, f_amr_tag_grad_sf
    use m_amr_distribution, only: s_amr_sort_by_key
    use m_amr, only: amr_slots, amr_cons_st, amr_stor_st, amr_loc_of, amr_slot_live, amr_my_blk, amr_n_my, &
        & s_amr_refresh_my_blocks, amr_maxc_fit, amr_seam_pairs_dirty, amr_mesh_epoch, amr_cpat_mar, s_amr_alloc_slot, &
        & s_amr_free_slot, s_amr_reduce_xchg_flag, s_amr_reconcile_slots, s_amr_assign_block_owners, &
        & s_amr_exchange_coarse_cons_halo, s_amr_phys_to_cells, s_amr_body_bbox, s_amr_expand_box_over_bodies, s_amr_tile_box, &
        & f_amr_seam_dim, f_amr_boxes_overlap, s_set_amr_fine_geometry, s_interpolate_coarse_to_fine, s_amr_setup_ib, &
        & s_amr_l1_fill_exchange, s_amr_l1_fill_consume, s_amr_parent_fill_exchange, s_amr_parent_fill_consume, &
        & s_amr_fill_wave_done, amr_fw_dev
    use m_acoustic_src, only: acoustic_supp_lo, acoustic_supp_hi
    use m_bubbles_EL, only: s_lag_cloud_bbox_local

    implicit none

    private
    public :: s_amr_regrid, s_amr_check_seam_topology

    !> Padded global coarse-index bbox of the Lagrangian bubble cloud; blocks and regrid boxes stay clear of it
    integer :: lag_supp_lo(3), lag_supp_hi(3)
    logical :: lag_supp_on = .false.
    !> Migration send/receive pools (freed after the regrid); amr_mig_win(:,old block) = payload start in amr_mig_rq, fine window
    !! lo(3), hi(3) (start 0 = none)
    real(wp), allocatable :: amr_mig_sq(:), amr_mig_rq(:)
    integer, allocatable  :: amr_mig_win(:,:)

contains

    !> Abort on same-level seam topologies no halo reconciles: partial-transverse face adjacency or box intersection (both arise
    !! only from IB body-bbox expansion). Each rank tests its own blocks against all blocks.
    impure subroutine s_amr_check_seam_topology()

        integer :: ix, xb, yb, d, t
        logical :: adj, tover

        call s_amr_refresh_my_blocks()
        do ix = 1, amr_n_my
            xb = amr_my_blk(ix)
            do yb = 1, amr_num_blocks
                if (xb == yb) cycle
                if (amr_block_level(xb) /= amr_block_level(yb)) cycle
                if (f_amr_boxes_overlap(amr_region_lo_all(:,xb), amr_region_hi_all(:,xb), amr_region_lo_all(:,yb), &
                    & amr_region_hi_all(:,yb))) then
                    write (0, '(A,6I6,A,6I6)') ' [amr] intersecting blocks (lo,hi): ', amr_region_lo_all(:,xb), &
                           & amr_region_hi_all(:,xb), ' |', amr_region_lo_all(:,yb), amr_region_hi_all(:,yb)
                    call s_mpi_abort("AMR: two same-level blocks INTERSECT (the child IB body-bbox expansion route can " &
                                     & // "produce this - it has no overlap-merge pass): the overlapping cells would be " &
                                     & // "restricted and refluxed twice, silently breaking conservation. Adjust the body/" &
                                     & // "regrid inputs so body-expanded child boxes merge or separate.")
                end if
                do d = 1, num_dims
                    adj = amr_region_lo_all(d, yb) == amr_region_hi_all(d, xb) + 1
                    if (.not. adj) cycle
                    tover = .true.
                    do t = 1, num_dims
                        if (t /= d) tover = tover .and. amr_region_lo_all(t, xb) <= amr_region_hi_all(t, &
                            & yb) .and. amr_region_lo_all(t, yb) <= amr_region_hi_all(t, xb)
                    end do
                    if (tover .and. f_amr_seam_dim(xb, yb) == 0) then
                        write (0, '(A,6I6,A,6I6)') ' [amr] touching blocks (lo,hi): ', amr_region_lo_all(:,xb), &
                               & amr_region_hi_all(:,xb), ' |', amr_region_lo_all(:,yb), amr_region_hi_all(:,yb)
                        call s_mpi_abort("AMR: two same-level blocks touch with PARTIAL transverse overlap (IB body-bbox " &
                                         & // "expansion can produce this): the fine-fine seam halo only reconciles exact-match " &
                                         & // "faces, so the shared-face flux would silently leak. Adjust amr_block/regrid inputs " // "so body-expanded boxes merge or separate.")
                    end if
                end do
            end do
        end do

    end subroutine s_amr_check_seam_topology

    !> Abort if a box exceeds its level's slot cap, amr_maxc_fit/amr_ref_ratio**(lev-1) coarse cells (it would overrun the fine
    !! coord arrays); floor of 1 as in s_amr_tile_box.
    impure subroutine s_amr_check_box_caps(boxes, nboxes, box_level)

        type(t_box), intent(in) :: boxes(:)
        integer, intent(in)     :: nboxes, box_level(:)
        integer                 :: k, d, lev, cap, span

        do k = 1, nboxes
            lev = box_level(k)
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

    !> Abort if two same-level boxes overlap (the overlap carry-forward and exchange plans assume disjointness). All levels share
    !! the coarse index space; an x-lo sweep keeps it cheap.
    impure subroutine s_amr_check_box_disjoint(boxes, nboxes, box_level)

        type(t_box), intent(in) :: boxes(:)
        integer, intent(in)     :: nboxes, box_level(:)
        integer(8), allocatable :: key(:)
        integer, allocatable    :: ord(:)
        integer                 :: i, j, k, kk

        allocate (ord(nboxes)); key = [(int(boxes(k)%lo(1), 8), k=1, nboxes)]
        call s_amr_sort_by_key(key, nboxes, ord)
        do i = 1, nboxes
            k = ord(i)
            if (box_level(k) < 1) cycle
            do j = i + 1, nboxes
                kk = ord(j)
                if (boxes(kk)%lo(1) > boxes(k)%hi(1)) exit
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

    !> Clip a box (global indices) clear of every acoustic source support; may empty it (hi < lo).
    impure subroutine s_amr_clip_box_from_sources(lo, hi)

        integer, intent(inout) :: lo(3), hi(3)
        integer                :: s

        do s = 1, num_source
            call s_amr_clip_box_from_supp(lo, hi, acoustic_supp_lo(:,s), acoustic_supp_hi(:,s))
        end do
        do s = 1, num_source
            if (hi(1) < lo(1) .or. hi(2) < lo(2) .or. hi(3) < lo(3)) return
            if (f_amr_boxes_overlap(lo, hi, acoustic_supp_lo(:,s), acoustic_supp_hi(:, &
                & s))) call s_mpi_abort('amr regrid: acoustic source exclusion clip failed (internal error)')
        end do

    end subroutine s_amr_clip_box_from_sources

    !> Collectively recompute the global bubble-cloud exclusion bbox, padded by pad_cells; lag_supp_on = false when no rank holds a
    !! bubble.
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
        call s_amr_phys_to_cells(pmin_glb, pmax_glb, pad_cells, lag_supp_lo, lag_supp_hi)

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

    !> Clip a box (global indices) clear of one support bbox along the axis/side keeping the most extent (ties: lower axis, begin
    !! side); may empty it.
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

    !> Clip a box clear of the acoustic source and bubble-cloud supports.
    impure subroutine s_amr_clip_box_clear(lo, hi)

        integer, intent(inout) :: lo(3), hi(3)

        if (acoustic_source) call s_amr_clip_box_from_sources(lo, hi)
        if (bubbles_lagrange .and. lag_supp_on) call s_amr_clip_box_from_supp(lo, hi, lag_supp_lo, lag_supp_hi)

    end subroutine s_amr_clip_box_clear

    !> This rank's tagged cells as global level-0 coordinates; rank interiors are disjoint, so a SUM reduction counts each cell
    !! once.
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

    !> Regrid: tag by density gradient, cluster, shape, nest and rebuild the slots, carrying old data into the overlaps. Called
    !! between steps; no-op if nothing is tagged or the box set is unchanged.
    impure subroutine s_amr_regrid(q_cons_base)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_base
        logical, allocatable                                   :: tag_grid(:,:,:)
        type(t_box), allocatable                               :: boxes(:)
        integer                                                :: sidx(3), nboxes
        integer                                                :: old_np
        ! heap: O(global boxes) arrays overflow the default stack
        integer, allocatable :: box_level(:)
        integer, allocatable :: old_ilo(:,:), old_ext(:,:)
        integer, allocatable :: old_level(:)
        logical, allocatable :: old_owns(:)
        logical              :: same
        integer              :: i

        allocate (box_level(amr_max_blocks), old_ilo(3, amr_max_blocks), old_ext(3, amr_max_blocks), old_level(amr_max_blocks), &
                  & old_owns(amr_max_blocks))

        ! valid coarse ghosts at rank boundaries: the tag sweep and rebuild prolongation read across them

        call s_phase_tic(PH_RG_HALO)
        call s_amr_exchange_coarse_cons_halo(q_cons_base)
        do i = 1, sys_size
            $:GPU_UPDATE(host='[q_cons_base(i)%sf]')
        end do
        call s_phase_toc(PH_RG_HALO)

        ! cloud pad: smearing (mapCells) + stencil (2) + drift (amr_buf)
        if (bubbles_lagrange) call s_amr_compute_lag_supp(mapCells + 2 + amr_buf)

        call s_phase_tic(PH_RG_TAG)
        call s_amr_regrid_tag_cells(q_cons_base, tag_grid, sidx)
        call s_phase_toc(PH_RG_TAG)
        call s_phase_tic(PH_RG_CLUS)
        call s_amr_regrid_cluster_tags(tag_grid, sidx, boxes, nboxes)
        call s_phase_toc(PH_RG_CLUS)
        if (nboxes == 0) return  ! nothing tagged on any rank; keep the current blocks
        call s_amr_regrid_shape_boxes(boxes, nboxes)
        if (nboxes == 0) return  ! every box was confined to the domain margin
        call s_phase_tic(PH_RG_CLUS)
        call s_amr_regrid_nest_children(boxes, nboxes, box_level)
        if (amr_snap > 0) call s_amr_regrid_snap_boxes(boxes, nboxes, box_level)
        call s_amr_check_box_caps(boxes, nboxes, box_level)  ! invariant: no box may exceed its level's slot cap
        call s_amr_check_box_disjoint(boxes, nboxes, box_level)  ! invariant: same-level boxes are pairwise disjoint
        call s_amr_report_coverage(boxes, nboxes, box_level)
        call s_amr_regrid_boxes_unchanged(boxes, nboxes, box_level, same)
        call s_phase_toc(PH_RG_CLUS)
        if (same) return  ! identical box set and levels: keep the live slots
        call s_amr_regrid_stash_migrate(boxes, nboxes, box_level, old_np, old_ilo, old_ext, old_level, old_owns)
        call s_amr_regrid_rebuild_slots(q_cons_base, boxes, nboxes, old_np, old_ilo, old_ext, old_level)

    end subroutine s_amr_regrid

    !> Print the share of the base domain covered by each level's boxes (rank 0, rank_time_wrt).
    impure subroutine s_amr_report_coverage(boxes, nboxes, box_level)

        type(t_box), intent(in) :: boxes(:)
        integer, intent(in)     :: nboxes, box_level(:)
        real(wp)                :: cover(1:max(amr_max_level, 1)), base
        integer                 :: k

        if (.not. rank_time_wrt .or. proc_rank /= 0) return
        cover = 0._wp
        do k = 1, nboxes
            cover(box_level(k)) = cover(box_level(k)) + real(product(merge(boxes(k)%hi - boxes(k)%lo + 1, 1, amr_dim)), wp)
        end do
        base = real(m_glb + 1, wp)*real(n_glb + 1, wp)*real(max(p_glb + 1, 1), wp)
        print '(A,I0,A,*(F8.4))', ' [amr-cov] boxes ', nboxes, '  level coverage (share of domain) ', cover/base

    end subroutine s_amr_report_coverage

    !> Regrid phase 1: density-gradient tag field, skipping the global boundary cells; tags over acoustic and bubble-cloud supports
    !! are suppressed. sidx = this rank's global start offsets.
    impure subroutine s_amr_regrid_tag_cells(q_cons_base, tag_grid, sidx)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_base
        logical, allocatable, intent(inout)                 :: tag_grid(:,:,:)
        integer, intent(out)                                :: sidx(3)
        integer                                             :: tg_lo(3), tg_hi(3), ci, cj, ck
        real(wp)                                            :: g

        sidx = amr_sidx
        tg_lo = merge(merge(1, 0, sidx == 0), 0, amr_dim)
        tg_hi = merge(merge(amr_ext - 1, amr_ext, sidx + amr_ext == [m_glb, n_glb, p_glb]), 0, amr_dim)
        allocate (tag_grid(0:m,0:n,0:p)); tag_grid = .false.
        do ck = tg_lo(3), tg_hi(3)
            do cj = tg_lo(2), tg_hi(2)
                do ci = tg_lo(1), tg_hi(1)
                    ! total density (sum of continuity variables): single-fluid tagger, immune to trace-fluid noise
                    g = f_amr_tag_grad_sf(q_cons_base, ci, cj, ck)
                    if (g > amr_tag_eps) tag_grid(ci, cj, ck) = .true.
                    ! acoustic source support stays coarse
                    if (acoustic_source .and. tag_grid(ci, cj, ck)) then
                        if (f_in_acoustic_support(ci + sidx(1), cj + sidx(2), ck + sidx(3))) tag_grid(ci, cj, ck) = .false.
                    end if
                    ! bubble cloud stays coarse (two-way coupling lives there)
                    if (bubbles_lagrange .and. tag_grid(ci, cj, ck)) then
                        if (f_in_lag_support(ci + sidx(1), cj + sidx(2), ck + sidx(3))) tag_grid(ci, cj, ck) = .false.
                    end if
                end do
            end do
        end do

    end subroutine s_amr_regrid_tag_cells

    !> Regrid phase 2: cluster this rank's tags (rank-invariant, via signature reductions) into separated boxes; nboxes = 0 if
    !! nothing is tagged.
    impure subroutine s_amr_regrid_cluster_tags(tag_grid, sidx, boxes, nboxes)

        logical, allocatable, intent(inout)   :: tag_grid(:,:,:)
        integer, intent(in)                   :: sidx(3)
        type(t_box), allocatable, intent(out) :: boxes(:)
        integer, intent(out)                  :: nboxes
        integer, allocatable                  :: tags(:,:)
        integer                               :: ntag

        call s_amr_local_tags(tag_grid, sidx, tags, ntag)
        deallocate (tag_grid)
        call s_amr_cluster(tags, ntag, boxes, nboxes, .true.)
        deallocate (tags)

    end subroutine s_amr_regrid_cluster_tags

    !> True iff two same-level blocks can coexist: disjoint, and either >= 2 empty coarse cells apart or face-adjacent with an exact
    !! transverse match. Mirrors s_amr_check_seam_topology.
    pure logical function f_amr_boxes_compatible(alo, ahi, blo, bhi) result(ok)

        integer, intent(in) :: alo(3), ahi(3), blo(3), bhi(3)
        integer             :: d, t, gap
        logical             :: tmatch, tover

        ok = .true.
        if (all(alo <= bhi .and. blo <= ahi)) then; ok = .false.; return; end if
        do d = 1, 3
            if (.not. amr_dim(d)) cycle
            gap = max(blo(d) - ahi(d), alo(d) - bhi(d)) - 1  ! empty cells along d (< 0: overlap)
            if (gap >= 2 .or. gap < 0) cycle
            tover = .true.; tmatch = .true.
            do t = 1, 3
                if (t == d .or. .not. amr_dim(t)) cycle
                tover = tover .and. alo(t) <= bhi(t) .and. blo(t) <= ahi(t)
                tmatch = tmatch .and. alo(t) == blo(t) .and. ahi(t) == bhi(t)
            end do
            if (tover .and. .not. (gap == 0 .and. tmatch)) then; ok = .false.; return; end if
        end do

    end function f_amr_boxes_compatible

    !> Hysteresis (amr_snap > 0): a box within amr_snap cells per face of a live same-level block takes that block's box.
    !! All-or-none: the snapped set must stay disjoint, compatible and singly nested (within amr_cpat_mar), else the fresh boxes
    !! stand.
    impure subroutine s_amr_regrid_snap_boxes(boxes, nboxes, box_level)

        type(t_box), intent(inout) :: boxes(:)
        integer, intent(in)        :: nboxes, box_level(:)
        type(t_box), allocatable   :: snapped(:)
        integer                    :: k, kk, ks, npar, nsnap, mlo(3), mhi(3)
        logical                    :: ok

        allocate (snapped(nboxes)); snapped(1:nboxes) = boxes(1:nboxes)
        nsnap = 0
        do k = 1, nboxes
            do ks = 1, amr_num_blocks
                if (amr_block_level(ks) /= box_level(k)) cycle
                if (all(abs(amr_region_lo_all(:,ks) - boxes(k)%lo) <= amr_snap) .and. all(abs(amr_region_hi_all(:, &
                    & ks) - boxes(k)%hi) <= amr_snap)) then
                    ! with IB the old block must contain the fresh box (a lagging body margin would be lost)
                    if (ib .and. .not. (all(amr_region_lo_all(:,ks) <= boxes(k)%lo) .and. all(amr_region_hi_all(:, &
                        & ks) >= boxes(k)%hi))) cycle
                    if (any(amr_region_lo_all(:,ks) /= boxes(k)%lo) .or. any(amr_region_hi_all(:, &
                        & ks) /= boxes(k)%hi)) nsnap = nsnap + 1
                    snapped(k)%lo = amr_region_lo_all(:,ks); snapped(k)%hi = amr_region_hi_all(:,ks)
                    exit
                end if
            end do
        end do
        ok = nsnap > 0
        do k = 1, nboxes
            if (.not. ok) exit
            do kk = k + 1, nboxes
                if (box_level(kk) /= box_level(k)) cycle
                if (all(snapped(k)%lo <= snapped(kk)%hi .and. snapped(kk)%lo <= snapped(k)%hi)) then
                    ok = .false.; exit
                end if
                if (.not. f_amr_boxes_compatible(snapped(k)%lo, snapped(k)%hi, snapped(kk)%lo, snapped(kk)%hi)) then
                    ok = .false.; exit
                end if
            end do
        end do
        ! a level >= 2 box lies inside exactly one parent box, inset by the nesting margin
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
        deallocate (snapped)

    end subroutine s_amr_regrid_snap_boxes

    !> Append to out(n+1:) the pieces of `box` that stay `gap` empty cells clear of [glo:ghi] and of each other, cutting dim by dim;
    !! the part inside the block's span is dropped. The caller guarantees `box` overlaps the gap-grown block.
    pure subroutine s_amr_box_carve(box, glo, ghi, gap, out, n)

        type(t_box), intent(in)    :: box
        integer, intent(in)        :: glo(3), ghi(3), gap
        type(t_box), intent(inout) :: out(:)
        integer, intent(inout)     :: n
        type(t_box)                :: rem
        integer                    :: d

        rem = box
        do d = 1, 3
            if (.not. amr_dim(d)) cycle
            if (rem%lo(d) < glo(d) - gap) then
                n = n + 1; out(n) = rem; out(n)%hi(d) = glo(d) - gap - 1
            end if
            if (rem%hi(d) > ghi(d) + gap) then
                n = n + 1; out(n) = rem; out(n)%lo(d) = ghi(d) + gap + 1
            end if
            rem%lo(d) = max(rem%lo(d), glo(d)); rem%hi(d) = min(rem%hi(d), ghi(d))
            if (rem%lo(d) > rem%hi(d)) return
        end do

    end subroutine s_amr_box_carve

    !> Abort if ranks hold different box lists (each rank sees only nearby IB bodies; a mismatch hangs the ownership collectives).
    impure subroutine s_amr_check_boxes_agree(boxes, nboxes)

        type(t_box), intent(in) :: boxes(:)
        integer, intent(in)     :: nboxes
        real(wp)                :: chk, cmin, cmax
        integer                 :: k

        chk = real(nboxes, wp)
        do k = 1, nboxes
            chk = chk + real(modulo(k*7919, 1009), wp)*real(sum(boxes(k)%lo*(/1, 131, 17161/)) + 3*sum(boxes(k)%hi*(/1, 131, &
                             & 17161/)), wp)
        end do
        call s_mpi_allreduce_min(chk, cmin)
        call s_mpi_allreduce_max(chk, cmax)
        if (cmin /= cmax) then
            call s_mpi_abort('amr regrid with ib: the ranks built different block lists, because they hold different sets of ' &
                             & // 'immersed bodies; raise ib_neighborhood_radius until every rank sees every body')
        end if

    end subroutine s_amr_check_boxes_agree

    !> Fold body group j into group i (bbox and body-extent unions) and move the last group into slot j.
    pure subroutine s_amr_merge_body_groups(glo, ghi, clo, chi, ng, i, j)

        integer, intent(inout) :: glo(:,:), ghi(:,:), clo(:,:), chi(:,:), ng
        integer, intent(in)    :: i, j

        glo(:,i) = min(glo(:,i), glo(:,j)); ghi(:,i) = max(ghi(:,i), ghi(:,j))
        clo(:,i) = min(clo(:,i), clo(:,j)); chi(:,i) = max(chi(:,i), chi(:,j))
        glo(:,j) = glo(:,ng); ghi(:,j) = ghi(:,ng); clo(:,j) = clo(:,ng); chi(:,j) = chi(:,ng); ng = ng - 1

    end subroutine s_amr_merge_body_groups

    !> One margin-padded box (glo:ghi) and body extent (clo:chi) per immersed body, identical on every rank.
    impure subroutine s_amr_gather_body_groups(glo, ghi, clo, chi, ng, mrg, nbody)

        integer, intent(inout) :: glo(:,:), ghi(:,:), clo(:,:), chi(:,:)
        integer, intent(out)   :: ng
        integer, intent(in)    :: mrg, nbody
        integer, parameter     :: BOX_UNSET = -1000000000
        integer, allocatable   :: gbx(:,:)
        integer                :: i, gid, ierr, blo(3), bhi(3), cl(3), ch(3)

        ! Every rank must build identical boxes though it sees only nearby bodies whose positions differ in the last digits: each
        ! rank fills the boxes of its bodies by global id and one MAX reduction shares them.

        allocate (gbx(12, nbody)); gbx = BOX_UNSET
        do i = 1, num_ibs
            gid = patch_ib(i)%gbl_patch_id
            if (gid < 1 .or. gid > nbody) gid = i
            call s_amr_body_bbox(i, mrg, blo, bhi); call s_amr_body_bbox(i, 0, cl, ch)
            gbx(1:3,gid) = blo; gbx(4:6,gid) = bhi; gbx(7:9,gid) = cl; gbx(10:12,gid) = ch
        end do
#ifdef MFC_MPI
        call MPI_ALLREDUCE(MPI_IN_PLACE, gbx, size(gbx), MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
#endif
        ng = 0
        do gid = 1, nbody
            if (gbx(1, gid) <= BOX_UNSET) cycle
            blo = gbx(1:3,gid); bhi = gbx(4:6,gid)
            ! blocks stay buff_size inside the domain
            if (any(amr_dim .and. (blo < buff_size .or. bhi > [m_glb, n_glb, p_glb] - buff_size))) then
                call s_mpi_abort('amr dynamic regrid with ib: the immersed body plus its containment ' &
                                 & // 'margin does not fit inside the refinable domain interior (blocks stay buff_size off the edges)')
            end if
            ng = ng + 1; glo(:,ng) = blo; ghi(:,ng) = bhi; clo(:,ng) = gbx(7:9,gid); chi(:,ng) = gbx(10:12,gid)
        end do
        deallocate (gbx)

    end subroutine s_amr_gather_body_groups

    !> Separate body groups that face each other across < gap empty cells or overlap, reject groups over the block cap, and order
    !! them canonically.
    impure subroutine s_amr_resolve_body_groups(glo, ghi, clo, chi, ng, gap)

        integer, intent(inout) :: glo(:,:), ghi(:,:), clo(:,:), chi(:,:), ng
        integer, intent(in)    :: gap
        integer                :: e(3), d, dsep, nsep, lowb, upb, need, avail_l, avail_u, tl, tu, blo(3), bhi(3)
        integer                :: i, j, g
        logical                :: merged

        ! trim the margins apart (floored at 3 cells beyond the body: image-point stencil and moving-body guard), else merge

        merged = .true.
        do while (merged)
            merged = .false.
            outer: do i = 1, ng - 1
                do j = i + 1, ng
                    e = 0
                    do d = 1, 3
                        if (amr_dim(d)) e(d) = max(glo(d, j) - ghi(d, i), glo(d, i) - ghi(d, j)) - 1
                    end do
                    nsep = count(amr_dim .and. e >= 0)  ! dimensions in which the two boxes are apart
                    if (nsep == 0) then
                        call s_amr_merge_body_groups(glo, ghi, clo, chi, ng, i, j)
                        merged = .true.; exit outer
                    else if (nsep == 1) then
                        dsep = 0
                        do d = 1, 3
                            if (amr_dim(d) .and. e(d) >= 0) dsep = d
                        end do
                        if (e(dsep) < gap) then
                            need = gap - e(dsep)
                            lowb = merge(i, j, ghi(dsep, i) < ghi(dsep, j)); upb = merge(j, i, lowb == i)
                            avail_l = max(0, ghi(dsep, lowb) - (chi(dsep, lowb) + 3))
                            avail_u = max(0, (clo(dsep, upb) - 3) - glo(dsep, upb))
                            if (avail_l + avail_u >= need) then
                                tl = min(need, avail_l); tu = need - tl
                                ghi(dsep, lowb) = ghi(dsep, lowb) - tl; glo(dsep, upb) = glo(dsep, upb) + tu
                            else
                                call s_amr_merge_body_groups(glo, ghi, clo, chi, ng, i, j)
                                merged = .true.; exit outer
                            end if
                        end if
                    end if
                end do
            end do outer
        end do
        do g = 1, ng
            if (any(amr_dim .and. ghi(:,g) - glo(:,g) + 1 > amr_maxc_fit)) then
                write (0, '(A,I0,A,I0,A,3I6,A,3I6,A,3I6,A,I0)') ' [amr] body group ', g, ' of ', ng, ' lo ', glo(:,g), ' hi ', &
                       & ghi(:,g), ' cap ', amr_maxc_fit, ' num_ibs ', num_ibs
                call s_mpi_abort('amr dynamic regrid with ib: a body (or bodies closer than the block gap) plus ' &
                                 & // 'margin exceeds the per-rank block size cap; use fewer ranks or a larger amr_maxc_fit')
            end if
        end do

        ! canonical order by low corner: patch_ib order differs per rank and the carve depends on it
        do i = 2, ng
            blo = glo(:,i); bhi = ghi(:,i); j = i - 1
            do while (j >= 1)
                if (.not. (glo(3, j) > blo(3) .or. (glo(3, j) == blo(3) .and. (glo(2, j) > blo(2) .or. (glo(2, &
                    & j) == blo(2) .and. glo(1, j) > blo(1)))))) exit
                glo(:,j + 1) = glo(:,j); ghi(:,j + 1) = ghi(:,j); j = j - 1
            end do
            glo(:,j + 1) = blo; ghi(:,j + 1) = bhi
        end do

    end subroutine s_amr_resolve_body_groups

    !> Regrid phase 3a (IB): tile the boxes around the immersed bodies. Each body (bbox + margin; close bodies share one) owns a
    !! block; boxes give up that block plus a gap ring, the remainder tiles to <= amr_maxc_fit, and body blocks are appended last.
    !! The ring keeps body blocks off other blocks, since the seam halo and reflux need exact-match transverse faces.
    impure subroutine s_amr_carve_tile_around_bodies(boxes, nboxes, out, nt, capt)

        type(t_box), intent(in)    :: boxes(:)
        integer, intent(in)        :: nboxes
        type(t_box), intent(inout) :: out(:)
        integer, intent(inout)     :: nt, capt
        ! empty coarse cells between blocks off a shared grid (one cell gives both the same reflux outside cell)
        integer, parameter       :: gap = 2
        type(t_box), allocatable :: pieces(:), rest(:)
        integer, allocatable     :: glo(:,:), ghi(:,:), clo(:,:), chi(:,:)
        integer                  :: nbody, ng, mrg, g, k, p, npc, nn
        logical, allocatable     :: hit(:)

        ! margin as in s_amr_expand_box_over_bodies: image-point stencil plus nested-level insets

        mrg = max(amr_buf, 4) + max(0, amr_max_level - 1)*amr_cpat_mar
        nbody = max(num_gbl_ibs, num_ibs, 1)
        allocate (glo(3, nbody), ghi(3, nbody), clo(3, nbody), chi(3, nbody), hit(nbody))
        call s_amr_gather_body_groups(glo, ghi, clo, chi, ng, mrg, nbody)
        call s_amr_resolve_body_groups(glo, ghi, clo, chi, ng, gap)

        allocate (pieces(amr_max_blocks), rest(amr_max_blocks))
        hit = .false.
        do k = 1, nboxes
            npc = 1; pieces(1) = boxes(k)
            do g = 1, ng
                nn = 0
                do p = 1, npc
                    if (f_amr_boxes_overlap(pieces(p)%lo, pieces(p)%hi, glo(:,g) - gap, ghi(:,g) + gap)) then
                        if (f_amr_boxes_overlap(pieces(p)%lo, pieces(p)%hi, glo(:,g), ghi(:,g))) hit(g) = .true.
                        if (nn + 2*num_dims > amr_max_blocks) call s_mpi_abort('amr regrid: too many pieces around a body')
                        call s_amr_box_carve(pieces(p), glo(:,g), ghi(:,g), gap, rest, nn)
                    else
                        nn = nn + 1; rest(nn) = pieces(p)
                    end if
                end do
                pieces(1:nn) = rest(1:nn); npc = nn
            end do
            do p = 1, npc
                call s_amr_tile_box(pieces(p)%lo, pieces(p)%hi, out, nt, amr_max_blocks, capt)
            end do
        end do

        do g = 1, ng
            if (.not. hit(g)) cycle
            if (nt >= amr_max_blocks) then
                capt = 1
            else
                nt = nt + 1; out(nt)%lo = glo(:,g); out(nt)%hi = ghi(:,g)
            end if
        end do

    end subroutine s_amr_carve_tile_around_bodies

    !> Regrid phase 3: pad and clamp each box, clip it clear of the acoustic/bubble supports, then tile oversized boxes (IB: after
    !! carving out the body blocks).
    impure subroutine s_amr_regrid_shape_boxes(boxes, nboxes)

        type(t_box), allocatable, intent(inout) :: boxes(:)
        integer, intent(inout)                  :: nboxes
        type(t_box), allocatable                :: tiled(:)
        integer                                 :: lo(3), hi(3), k, kk, ntl, capt, pad

        pad = amr_buf
        k = 0
        do kk = 1, nboxes
            lo = merge(max(boxes(kk)%lo - pad, buff_size), 0, amr_dim)
            hi = merge(min(boxes(kk)%hi + pad, [m_glb, n_glb, p_glb] - buff_size), 0, amr_dim)
            ! the supports are coarse-grid only; clipping only shrinks, so boxes stay disjoint
            call s_amr_clip_box_clear(lo, hi)
            if (any(hi < lo)) cycle  ! confined to the domain margin
            k = k + 1; boxes(k)%lo = lo; boxes(k)%hi = hi
        end do
        nboxes = k
        if (nboxes == 0) return

        ! tile boxes larger than amr_maxc_fit (IB: after carving out body blocks, so a body is never split)
        allocate (tiled(amr_max_blocks))
        ntl = 0; capt = 0
        if (ib) then
            call s_amr_carve_tile_around_bodies(boxes, nboxes, tiled, ntl, capt)
        else
            do kk = 1, nboxes
                call s_amr_tile_box(boxes(kk)%lo, boxes(kk)%hi, tiled, ntl, amr_max_blocks, capt)
            end do
        end if
        if (capt == 1 .and. proc_rank == 0) print '(A,I0)', ' [amr] WARNING: tiling capped at amr_max_blocks = ', amr_max_blocks
        deallocate (boxes); call move_alloc(tiled, boxes)
        nboxes = ntl
        if (ib) call s_amr_check_boxes_agree(boxes, nboxes)

        ! a body block cannot both contain the body and exclude a source/cloud: fail closed
        if (ib .and. (acoustic_source .or. (bubbles_lagrange .and. lag_supp_on))) then
            do k = 1, nboxes
                lo = boxes(k)%lo; hi = boxes(k)%hi
                call s_amr_clip_box_clear(lo, hi)
                if (any(lo /= boxes(k)%lo) .or. any(hi /= boxes(k)%hi)) then
                    call s_mpi_abort('amr regrid: a block must contain an immersed body AND stay ' &
                                     & // 'clear of an acoustic source support / Lagrangian bubble cloud - the ' &
                                     & // 'constraints conflict; move the body, source, or cloud apart')
                end if
            end do
        end if

    end subroutine s_amr_regrid_shape_boxes

    !> True iff the box set and levels equal the live slots (a rebuild would reproduce them).
    impure subroutine s_amr_regrid_boxes_unchanged(boxes, nboxes, box_level, same)

        type(t_box), intent(in) :: boxes(:)
        integer, intent(in)     :: nboxes, box_level(:)
        logical, intent(out)    :: same
        integer                 :: k, ks

        same = .false.
        if (nboxes == amr_num_blocks) then
            same = .true.
            do k = 1, nboxes
                if (any(boxes(k)%lo /= amr_slots(k)%region%lo) .or. any(boxes(k)%hi /= amr_slots(k)%region%hi) .or. box_level(k) &
                    & /= amr_block_level(k)) same = .false.
            end do
        end if

    end subroutine s_amr_regrid_boxes_unchanged

    !> Pack fine window [l, h] of an owned old block's stash (device store) into its wire slice buf; gi fastest, then gj, gk, ii.
    impure subroutine s_amr_mig_pack_device(loc, l1, l2, l3, h1, h2, h3, buf)

        integer, intent(in)                 :: loc, l1, l2, l3, h1, h2, h3
        real(wp), intent(inout), contiguous :: buf(:)
        integer                             :: ii, gk, gj, gi, n1, n2, n3

        n1 = h1 - l1 + 1; n2 = h2 - l2 + 1; n3 = h3 - l3 + 1
        $:GPU_PARALLEL_LOOP(collapse=4, copyout='[buf]')
        do ii = 1, sys_size
            do gk = l3, h3
                do gj = l2, h2
                    do gi = l1, h1
                        buf(1 + gi - l1 + n1*(gj - l2 + n2*(gk - l3 + n3*(ii - 1)))) = real(amr_stor_st(gi, gj, gk, ii, loc), wp)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_mig_pack_device

    !> Fine-index window [flo, fhi] in level-1 old block ilo's stash of the coarse box [wlo, whi] (0:0 in collapsed dims).
    pure subroutine s_amr_mig_window(wlo, whi, ilo, flo, fhi)

        integer, intent(in)  :: wlo(3), whi(3), ilo(3)
        integer, intent(out) :: flo(3), fhi(3)

        flo = merge(amr_ref_ratio*(wlo - ilo), 0, amr_dim)
        fhi = merge(amr_ref_ratio*(whi - ilo + 1) - 1, 0, amr_dim)

    end subroutine s_amr_mig_window

    !> Free a migration pool (device-resident when amr_fw_dev).
    impure subroutine s_amr_mig_free(pool)

        real(wp), allocatable, intent(inout) :: pool(:)

        if (.not. allocated(pool)) return
        if (amr_fw_dev) then
            $:GPU_EXIT_DATA(delete='[pool]')
        end if
        deallocate (pool)

    end subroutine s_amr_mig_free

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

    !> Overwrite the new block's prolonged cons with the covering old block's stashed detail (old fine index = new + sh) over the
    !! old window [l, h], from store column vold or the received wire slice buf.
    #:for SRC in ['st', 'wire']
        impure subroutine s_amr_overlap_copy_${SRC}$(vloc, vm, vn, vp, sh, l1, l2, l3, h1, h2, h3, &
            & #{if SRC == 'st'}#vold#{else}#buf#{endif}#)

            integer, intent(in) :: vloc, vm, vn, vp, sh(3), l1, l2, l3, h1, h2, h3

            #:if SRC == 'st'
                integer, intent(in) :: vold
            #:else
                real(wp), intent(in), contiguous :: buf(:)
            #:endif
            integer :: i, fi, fj, fk, ofi, ofj, ofk, sh1, sh2, sh3, n1, n2, n3
            logical :: vd2, vd3

            sh1 = sh(1); sh2 = sh(2); sh3 = sh(3)
            n1 = h1 - l1 + 1; n2 = h2 - l2 + 1; n3 = h3 - l3 + 1
            vd2 = n_glb > 0; vd3 = p_glb > 0
            #:if SRC == 'st'
                $:GPU_PARALLEL_LOOP(collapse=4, private='[ofi, ofj, ofk]')
            #:else
                $:GPU_PARALLEL_LOOP(collapse=4, private='[ofi, ofj, ofk]', copyin='[buf]')
            #:endif
            do i = 1, sys_size
                do fk = 0, vp
                    do fj = 0, vn
                        do fi = 0, vm
                            ofk = fk + sh3; ofj = fj + sh2; ofi = fi + sh1
                            if (vd3 .and. (ofk < l3 .or. ofk > h3)) cycle
                            if (vd2 .and. (ofj < l2 .or. ofj > h2)) cycle
                            if (ofi < l1 .or. ofi > h1) cycle
                            #:if SRC == 'st'
                                amr_cons_st(fi, fj, fk, i, vloc) = amr_stor_st(ofi, ofj, ofk, i, vold)
                            #:else
                                amr_cons_st(fi, fj, fk, i, &
                                            & vloc) = real(buf(1 + ofi - l1 + n1*(ofj - l2 + n2*(ofk - l3 + n3*(i - 1)))), stp)
                            #:endif
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

        end subroutine s_amr_overlap_copy_${SRC}$
    #:endfor

    !> Record the live slots' fine interiors on the device and the old block set (old_*; old_chi = replicated coarse hi corner).
    impure subroutine s_amr_stash_old_blocks(old_np, old_ilo, old_chi, old_ext, old_owner, old_level, old_owns)

        integer, intent(out) :: old_np, old_ilo(:,:), old_chi(:,:), old_ext(:,:), old_owner(:), old_level(:)
        logical, intent(out) :: old_owns(:)
        integer              :: k, loc_ks, old_ext1_k, old_ext2_k, old_ext3_k  !< plain scalars for device routines (NVHPC -Minline)

        old_np = amr_num_blocks
        do k = 1, old_np
            ! replicated global origin/extent (valid on every rank, unlike the owner-only isect)
            old_ilo(:,k) = amr_region_lo_all(:,k)
            old_chi(:,k) = amr_region_hi_all(:,k)
            ! fine extent = amr_ref_ratio**level*footprint - 1 (a level-2 block is 4x its footprint)
            old_ext(:,k) = merge(amr_ref_ratio**amr_block_level(k)*(old_chi(:,k) - old_ilo(:,k) + 1) - 1, 0, amr_dim)
            old_owner(k) = amr_block_owner(k)
            ! the overlap copy must match levels
            old_level(k) = amr_block_level(k)
            old_owns(k) = amr_owns_all(k)
            if (old_owns(k)) then
                ! device stash: the store is device-authoritative; s_amr_st_reserve preserves it across a mid-rebuild grow.
                ! amdflang: the kernel is its own routine because target regions nested in BLOCK drop from the device image (HSA
                ! invalid symbol)
                loc_ks = amr_loc_of(k)
                old_ext1_k = old_ext(1, k)
                old_ext2_k = old_ext(2, k)
                old_ext3_k = old_ext(3, k)
                call s_amr_stash_copy_device(loc_ks, old_ext1_k, old_ext2_k, old_ext3_k)
            end if
        end do

    end subroutine s_amr_stash_old_blocks

    !> Commit the new regions, levels and owners, after checking proper nesting.
    impure subroutine s_amr_commit_boxes(boxes, nboxes, box_level)

        type(t_box), intent(in) :: boxes(:)
        integer, intent(in)     :: nboxes, box_level(:)
        integer                 :: k, kk, npar

        ! regions and owners first: the migration needs the new owners and s_set_amr_fine_geometry reads amr_block_owner

        amr_num_blocks = nboxes
        do k = 1, nboxes
            amr_region_lo_all(:,k) = boxes(k)%lo; amr_region_hi_all(:,k) = boxes(k)%hi
            amr_block_level(k) = box_level(k)
        end do
        ! dirty the seam-pair and overlap-rank caches before the rebuild's gathers consume them
        amr_seam_pairs_dirty = .true.
        amr_mesh_epoch = amr_mesh_epoch + 1
        ! proper nesting: a level >= 2 block must overlap exactly one parent-level block (f_amr_parent_block takes the first)
        do k = 1, nboxes
            if (box_level(k) < 2) cycle
            npar = 0
            do kk = 1, nboxes
                if (box_level(kk) == box_level(k) - 1 .and. f_amr_boxes_overlap(boxes(k)%lo, boxes(k)%hi, boxes(kk)%lo, &
                    & boxes(kk)%hi)) npar = npar + 1
            end do
            if (npar /= 1) call s_mpi_abort('amr multi-level: a level>=2 block overlaps more than one (or no) ' &
                & // 'parent-level block - a fine tile straddling a parent-tile seam is unsupported (gather/reflux ' &
                & // 'couple to a single parent); reduce max_grid_size or the refined feature extent')
        end do
        amr_num_levels = maxval(box_level(1:nboxes))
        call s_amr_assign_block_owners()

    end subroutine s_amr_commit_boxes

    !> Cross-rank migration as one wave: each old owner ships the bounding box of its overlaps to every other rank owning an
    !! overlapping new level-1 block (level >= 2 re-prolongs from its parent). Both sides enumerate old blocks ascending, so the
    !! transfers line up; received windows stay in amr_mig_rq. No-op at np=1.
    impure subroutine s_amr_migrate_old_blocks(boxes, nboxes, box_level, old_np, old_ilo, old_chi, old_owner, old_level, old_owns)

        type(t_box), intent(in) :: boxes(:)
        integer, intent(in)     :: nboxes, box_level(:), old_np, old_ilo(:,:), old_chi(:,:), old_owner(:), old_level(:)
        logical, intent(in)     :: old_owns(:)
        integer                 :: h_loc               !< plain scalar for device routines (NVHPC -Minline)
        integer                 :: k2, kk, rr, cnt, ix, lo, hi, flo(3), fhi(3)
        integer, allocatable    :: wlo(:,:), whi(:,:)  !< per-rank overlap window of the current old block

        allocate (amr_mig_win(7, old_np)); amr_mig_win = 0
        if (num_procs > 1) then
            allocate (wlo(3,0:num_procs - 1), whi(3,0:num_procs - 1))
            wlo = huge(0); whi = -huge(0)
            call s_amr_wave_open(amr_wave, amr_band_migrate)
            call s_amr_wave_reset(amr_wsend); call s_amr_wave_reset(amr_wrecv)
            do kk = 1, old_np
                if (old_level(kk) /= 1) cycle
                do k2 = 1, nboxes
                    rr = amr_block_owner(k2)
                    if (box_level(k2) /= 1 .or. rr == old_owner(kk)) cycle
                    if (.not. f_amr_boxes_overlap(boxes(k2)%lo, boxes(k2)%hi, old_ilo(:,kk), old_chi(:,kk))) cycle
                    wlo(:,rr) = min(wlo(:,rr), max(boxes(k2)%lo, old_ilo(:,kk)))
                    whi(:,rr) = max(whi(:,rr), min(boxes(k2)%hi, old_chi(:,kk)))
                end do
                do rr = 0, num_procs - 1
                    if (whi(1, rr) < wlo(1, rr)) cycle
                    call s_amr_mig_window(wlo(:,rr), whi(:,rr), old_ilo(:,kk), flo, fhi)
                    cnt = sys_size*product(fhi - flo + 1)
                    if (rr == proc_rank) then
                        call s_amr_wave_add(amr_wrecv, old_owner(kk), kk, wlo(:,rr), whi(:,rr), cnt)
                    else if (old_owns(kk)) then
                        call s_amr_wave_add(amr_wsend, rr, kk, wlo(:,rr), whi(:,rr), cnt)
                    end if
                    wlo(:,rr) = huge(0); whi(:,rr) = -huge(0)
                end do
            end do
            deallocate (wlo, whi)
            call s_amr_wave_close(amr_wsend, amr_mig_sq, amr_fw_dev)
            call s_amr_wave_close(amr_wrecv, amr_mig_rq, amr_fw_dev)
            call s_amr_wave_post(amr_wave, amr_wrecv, amr_mig_rq, amr_fw_dev)
            do ix = 1, amr_wsend%nx
                kk = amr_wsend%blk(ix)
                call s_amr_wave_slice(amr_wsend, ix, lo, hi)
                call s_amr_mig_window(amr_wsend%bl(:,ix), amr_wsend%bh(:,ix), old_ilo(:,kk), flo, fhi)
                h_loc = amr_loc_of(kk)
                call s_amr_mig_pack_device(h_loc, flo(1), flo(2), flo(3), fhi(1), fhi(2), fhi(3), amr_mig_sq(lo:hi))
            end do
            call s_amr_wave_send(amr_wave, amr_wsend, amr_mig_sq, amr_fw_dev)
            call s_amr_wave_wait(amr_wave)
            call s_amr_mig_free(amr_mig_sq)
            do ix = 1, amr_wrecv%nx
                kk = amr_wrecv%blk(ix)
                call s_amr_wave_slice(amr_wrecv, ix, lo, hi)
                call s_amr_mig_window(amr_wrecv%bl(:,ix), amr_wrecv%bh(:,ix), old_ilo(:,kk), flo, fhi)
                amr_mig_win(:,kk) = [lo, flo, fhi]
            end do
        end if

    end subroutine s_amr_migrate_old_blocks

    !> Regrid phase 5: stash the old block set, commit the new one, and migrate stashed old blocks point-to-point to the ranks
    !! owning overlapping new blocks.
    impure subroutine s_amr_regrid_stash_migrate(boxes, nboxes, box_level, old_np, old_ilo, old_ext, old_level, old_owns)

        type(t_box), intent(in) :: boxes(:)
        integer, intent(in)     :: nboxes, box_level(:)
        integer, intent(out)    :: old_np, old_ilo(:,:), old_ext(:,:), old_level(:)
        logical, intent(out)    :: old_owns(:)
        integer                 :: old_chi(3, amr_max_blocks), old_owner(amr_max_blocks)

        call s_phase_tic(PH_RG_OWN)
        call s_amr_stash_old_blocks(old_np, old_ilo, old_chi, old_ext, old_owner, old_level, old_owns)
        call s_amr_commit_boxes(boxes, nboxes, box_level)
        call s_phase_toc(PH_RG_OWN)
        call s_phase_tic(PH_RG_MIG)
        call s_amr_migrate_old_blocks(boxes, nboxes, box_level, old_np, old_ilo, old_chi, old_owner, old_level, old_owns)
        call s_phase_toc(PH_RG_MIG)

    end subroutine s_amr_regrid_stash_migrate

    !> Visit order of the owned new boxes (level-major, vpos) and, per held old slot, the last visit that reads it (last_use).
    impure subroutine s_amr_plan_rebuild_visits(boxes, nboxes, old_np, old_ilo, old_ext, old_level, vpos, last_use, held, nh)

        type(t_box), intent(in)           :: boxes(:)
        integer, intent(in)               :: nboxes, old_np, old_ilo(:,:), old_ext(:,:), old_level(:)
        integer, allocatable, intent(out) :: vpos(:), last_use(:), held(:)
        integer, intent(out)              :: nh
        integer                           :: i, k, kk, lev, pos, ohi(3)

        call s_amr_refresh_my_blocks()
        allocate (vpos(nboxes)); vpos = 0; pos = 0
        do lev = 1, maxval(amr_block_level(1:amr_num_blocks))
            do i = 1, amr_n_my
                k = amr_my_blk(i)
                if (amr_block_level(k) /= lev) cycle
                pos = pos + 1; vpos(k) = pos
            end do
        end do
        allocate (last_use(old_np), held(old_np)); last_use = 0
        nh = 0
        do kk = 1, old_np
            if (.not. amr_slot_live(kk)) cycle
            nh = nh + 1; held(nh) = kk
            ohi = old_ilo(:,kk) + merge((old_ext(:,kk) + 1)/amr_ref_ratio**old_level(kk) - 1, 0, amr_dim)
            do i = 1, amr_n_my
                k = amr_my_blk(i)
                if (f_amr_boxes_overlap(boxes(k)%lo, boxes(k)%hi, old_ilo(:,kk), ohi)) last_use(kk) = max(last_use(kk), vpos(k))
            end do
        end do

    end subroutine s_amr_plan_rebuild_visits

    !> Overwrite the freshly prolonged level-1 slot ks with the overlapping old level-1 blocks: from the migration pool when another
    !! rank owned them, else from the local stash.
    impure subroutine s_amr_carry_old_into_slot(ks, boxes, old_np, old_ilo, old_ext, old_level)

        integer, intent(in)     :: ks, old_np, old_ilo(:,:), old_ext(:,:), old_level(:)
        type(t_box), intent(in) :: boxes(:)
        integer                 :: sh(3), kk, ohi(3)
        ! NVHPC -Minline: plain scalars for device routines
        integer :: loc_kks, loc_ks, m_ks, n_ks, old_ext1_kk, old_ext2_kk, old_ext3_kk, p_ks, w1, w2, w3, w4, w5, w6, lo, hi

        loc_ks = amr_loc_of(ks)
        m_ks = amr_slots(ks)%m
        n_ks = amr_slots(ks)%n
        p_ks = amr_slots(ks)%p
        do kk = 1, old_np
            if (old_level(kk) /= 1) cycle  ! a child's stash is in the 4x frame
            ohi = old_ilo(:,kk) + merge((old_ext(:,kk) + 1)/amr_ref_ratio - 1, 0, amr_dim)
            if (.not. f_amr_boxes_overlap(boxes(ks)%lo, boxes(ks)%hi, old_ilo(:,kk), ohi)) cycle
            sh = amr_ref_ratio*(amr_isect_lo - old_ilo(:,kk))  ! old local fine index = new local fine index + sh
            if (amr_mig_win(1, kk) > 0) then
                ! owned elsewhere: window in the migration receive pool
                lo = amr_mig_win(1, kk); hi = lo + sys_size*product(amr_mig_win(5:7,kk) - amr_mig_win(2:4,kk) + 1) - 1
                w1 = amr_mig_win(2, kk); w2 = amr_mig_win(3, kk); w3 = amr_mig_win(4, kk)
                w4 = amr_mig_win(5, kk); w5 = amr_mig_win(6, kk); w6 = amr_mig_win(7, kk)
                call s_amr_overlap_copy_wire(loc_ks, m_ks, n_ks, p_ks, sh, w1, w2, w3, w4, w5, w6, amr_mig_rq(lo:hi))
            else if (amr_slot_live(kk)) then
                loc_kks = amr_loc_of(kk)
                old_ext1_kk = old_ext(1, kk)
                old_ext2_kk = old_ext(2, kk)
                old_ext3_kk = old_ext(3, kk)
                call s_amr_overlap_copy_st(loc_ks, m_ks, n_ks, p_ks, sh, 0, 0, 0, old_ext1_kk, old_ext2_kk, old_ext3_kk, loc_kks)
            end if
        end do

    end subroutine s_amr_carry_old_into_slot

    !> Regrid phase 6: build each new slot level by level (prolong from the coarse/parent patch, overwrite the overlap from covering
    !! stashed old blocks), then reconcile the slot pool, rebuild the fine IB state and re-check the seam topology. Old-only slots
    !! are freed after their last reader so device memory does not peak.
    impure subroutine s_amr_regrid_rebuild_slots(q_cons_base, boxes, nboxes, old_np, old_ilo, old_ext, old_level)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_base
        type(t_box), intent(in)                                :: boxes(:)
        integer, intent(in)                                    :: nboxes, old_np, old_ilo(:,:), old_ext(:,:), old_level(:)
        integer                                                :: k, kk, i, hh, ks, kks, lev, nh, nlev
        integer, allocatable                                   :: last_use(:), held(:), vpos(:)

        call s_phase_tic(PH_RG_BUILD)
        do k = 1, nboxes
            ks = k
            if (amr_block_owner(ks) == proc_rank) cycle
            amr_cur = ks
            call s_set_amr_fine_geometry(boxes(k)%lo, boxes(k)%hi)
        end do

        ! nlev is replicated so every rank enters every wave
        call s_amr_plan_rebuild_visits(boxes, nboxes, old_np, old_ilo, old_ext, old_level, vpos, last_use, held, nh)
        nlev = maxval(amr_block_level(1:amr_num_blocks))

        do lev = 1, nlev
            if (lev == 1) then
                call s_amr_l1_fill_exchange(q_cons_base, .true.)
            else
                call s_amr_parent_fill_exchange(lev, .true.)
            end if
            do i = 1, amr_n_my
                ks = amr_my_blk(i)
                k = ks
                if (amr_block_level(ks) /= lev) cycle
                amr_cur = ks
                ! free old-only slots no later box reads (no-op if dead; slots serving new owned boxes live on)
                do hh = 1, nh
                    kk = held(hh)
                    if (last_use(kk) >= vpos(k)) cycle
                    kks = kk
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
                ! level >= 2 re-prolongs from its freshly built parent (its stash is in the parent-fine frame): no carry-forward
                call s_interpolate_coarse_to_fine()
                if (lev >= 2) cycle
                call s_amr_carry_old_into_slot(ks, boxes, old_np, old_ilo, old_ext, old_level)
            end do
            call s_amr_fill_wave_done()
        end do
        deallocate (last_use, held, vpos, amr_mig_win)
        call s_amr_mig_free(amr_mig_rq)
        call s_phase_toc(PH_RG_BUILD)
        call s_phase_tic(PH_RG_REC)

        ! sets amr_xchg_coarse_ghosts if any block needs it
        call s_amr_reduce_xchg_flag()
        ! free the transient slots of old blocks this rank does not now own
        call s_amr_reconcile_slots()
        ! rebuild fine IB state for the new geometry
        if (ib) call s_amr_setup_ib()
        call s_amr_select_slot(1)
        call s_amr_check_seam_topology()
        call s_phase_toc(PH_RG_REC)

    end subroutine s_amr_regrid_rebuild_slots

end module m_amr_regrid
