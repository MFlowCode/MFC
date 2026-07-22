!>
!!@file
!!@brief Contains module m_amr_regrid

#:include 'macros.fpp'

!> @brief Dynamic regrid for the block-structured AMR level set: density-gradient tagging, Berger-Rigoutsos clustering, box shaping
!! (pad/clamp/tile/IB merge), hierarchical child nesting, and the slot rebuild with cross-rank fine-state migration. Split out of
!! m_amr; the block/slot state stays in m_amr (and m_global_parameters) - this module only drives it.
module m_amr_regrid

#ifdef MFC_MPI
    use mpi  !< allgather/point-to-point for the rank-invariant tag union and the fine-state migration
#endif

    use m_derived_types  ! scalar_field, t_box
    use m_global_parameters
    use m_constants, only: mapCells
    use m_mpi_proxy, only: s_mpi_abort
    use m_mpi_common, only: s_mpi_allreduce_min, s_mpi_allreduce_max
    use m_amr, only: amr_slots, amr_maxc_fit, amr_seam_pairs_dirty, amr_xchg_coarse_ghosts, amr_cpat_mar, s_amr_alloc_slot, &
        & s_amr_reconcile_slots, s_amr_assign_block_owners, s_amr_gather_coarse_patch, s_amr_gather_coarse_patch_pbmv, &
        & s_amr_prolong_pbmv, s_amr_exchange_coarse_cons_halo, s_lag_phys_to_cells, s_amr_body_bbox, &
        & s_amr_expand_box_over_bodies, s_amr_tile_box, f_amr_seam_dim, f_amr_boxes_overlap, s_set_amr_fine_geometry, &
        & s_interpolate_coarse_to_fine, s_amr_setup_ib
    use m_acoustic_src, only: acoustic_supp_lo, acoustic_supp_hi
    use m_active_box, only: ab_x, ab_y, ab_z, ab_active
    use m_bubbles_EL, only: s_lag_cloud_bbox_local

    implicit none

    private
    public :: s_amr_regrid, s_amr_check_seam_topology, s_amr_check_active_box_containment

    !> Lagrangian bubble-cloud exclusion support: padded global coarse-index bbox of the cloud (positions + mapCells smearing +
    !! stencil headroom [+ drift margin at regrid]). Blocks and regrid boxes stay clear of it: a bubble inside a block loses two-way
    !! coupling (the fine advance skips the EL hooks and the coarse result under the block is discarded by restriction). Recomputed
    !! collectively at each regrid; guarded rank-locally per stage.
    integer :: lag_supp_lo(3), lag_supp_hi(3)
    logical :: lag_supp_on = .false.

contains

    !> Abort on same-level seam topologies no halo reconciles (silent conservation leaks otherwise). Run whenever the block set
    !! changes (regrid, restart); O(nblocks^2) on the replicated region metadata, identical on every rank. Two cases: (a) adjacency
    !! WITHOUT the exact transverse match f_amr_seam requires - reachable only through the IB body-bbox box expansion (clustering
    !! merges any too-close pair; tiling emits a regular grid) - which the fine-fine halo can never pair up; (b) an exact seam pair
    !! at level >= 2 under amr_subcycle, whose per-block child advance has no L2 halo (reachable via a restart mode-switch from a
    !! lockstep-produced layout; the subcycle regrid keeps one child per box); (c) same-level box INTERSECTION - reachable only
    !! through the CHILD IB body-bbox expansion, which unlike the L1 path has no overlap-merge pass - which double-restricts and
    !! double-refluxes the shared cells.
    impure subroutine s_amr_check_seam_topology()

        integer :: xb, yb, d, t
        logical :: adj, tover

        do xb = 1, amr_num_blocks
            do yb = 1, amr_num_blocks
                if (xb == yb) cycle
                if (amr_block_level(xb) /= amr_block_level(yb)) cycle
                ! same-level INTERSECTION (different levels legitimately nest; tiling emits disjoint tiles and the L1 IB
                ! pass merges overlapping boxes, but the CHILD IB body-bbox expansion has no overlap-merge pass)
                if (f_amr_boxes_overlap(amr_region_lo_all(:,xb), amr_region_hi_all(:,xb), amr_region_lo_all(:,yb), &
                    & amr_region_hi_all(:,yb))) then
                    call s_mpi_abort("AMR: two same-level blocks INTERSECT (the child IB body-bbox expansion route can " &
                                     & // "produce this - it has no overlap-merge pass): the overlapping cells would be " &
                                     & // "restricted and refluxed twice, silently breaking conservation. Adjust the body/" &
                                     & // "regrid inputs so body-expanded child boxes merge or separate.")
                end if
                if (amr_subcycle .and. amr_block_level(xb) >= 2 .and. f_amr_seam_dim(xb, yb) > 0) then
                    call s_mpi_abort("AMR: adjacent level-2+ blocks under amr_subcycle (e.g. a lockstep-written restart " &
                                     & // "resumed with amr_subcycle = T): the per-block child advance has no L2 seam halo, so the " // "shared-face flux would silently leak. Resume with amr_subcycle = F (lockstep).")
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

    !> Shrink box [blo:bhi] to the tight bounding box of the tagged cells inside it. ok=.false. if no tagged cell. Collapsed dims
    !! (lo=hi=0) survive unchanged. Deterministic (integer scan of the identical sparse tag list).
    impure subroutine s_amr_trim_box(tags, ts, te, blo, bhi, ok)

        integer, intent(in)    :: tags(:,:), ts, te
        integer, intent(inout) :: blo(3), bhi(3)
        logical, intent(out)   :: ok
        integer                :: tlo(3), thi(3), t, i, j, k

        tlo = huge(1); thi = -huge(1)
        do t = ts, te
            i = tags(1, t); j = tags(2, t); k = tags(3, t)
            if (i < blo(1) .or. i > bhi(1)) cycle
            if (j < blo(2) .or. j > bhi(2)) cycle
            if (k < blo(3) .or. k > bhi(3)) cycle
            tlo(1) = min(tlo(1), i); thi(1) = max(thi(1), i)
            tlo(2) = min(tlo(2), j); thi(2) = max(thi(2), j)
            tlo(3) = min(tlo(3), k); thi(3) = max(thi(3), k)
        end do
        ok = thi(1) >= tlo(1)
        if (ok) then; blo = tlo; bhi = thi; end if

    end subroutine s_amr_trim_box

    !> Berger-Rigoutsos bisection of one (already tagged-trimmed) candidate box on the global tag field: pick the longest splittable
    !! axis, prefer a zero-signature hole (widest interior run), else the strongest signature inflection (Laplacian sign change).
    !! ok=.false. if no axis admits a split leaving both children >= 2 cells. Integer-only => identical on all ranks.
    impure subroutine s_amr_find_split(tags, ts, te, blo, bhi, sax, spos, ok)

        integer, intent(in)  :: tags(:,:), ts, te
        integer, intent(in)  :: blo(3), bhi(3)
        integer, intent(out) :: sax, spos
        logical, intent(out) :: ok
        integer, parameter   :: min_child = 2
        integer              :: axord(3), ext(3), d, ax, t, i, j, k, s
        integer              :: run, run_start, best_run, best_start, lap, prevlap, bestmag, bestpos
        integer, allocatable :: sig(:)

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
        allocate (sig(0:maxval(bhi)))
        do d = 1, 3
            ax = axord(d)
            if (ax > num_dims) cycle
            if (ext(ax) < 2*min_child) cycle
            ! 1D signature sig(t) = count of in-box tagged cells at axis-position t (sum over the two transverse dims)
            do t = blo(ax), bhi(ax)
                sig(t) = 0
            end do
            do t = ts, te
                i = tags(1, t); j = tags(2, t); k = tags(3, t)
                if (i < blo(1) .or. i > bhi(1)) cycle
                if (j < blo(2) .or. j > bhi(2)) cycle
                if (k < blo(3) .or. k > bhi(3)) cycle
                select case (ax)
                case (1); sig(i) = sig(i) + 1
                case (2); sig(j) = sig(j) + 1
                case (3); sig(k) = sig(k) + 1
                end select
            end do
            ! (1) widest interior zero run (box is trimmed => sig(blo)>0 and sig(bhi)>0, so any run is interior)
            best_run = 0; best_start = -1; run = 0; run_start = -1
            do t = blo(ax), bhi(ax)
                if (sig(t) == 0) then
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
                    sax = ax; ok = .true.; deallocate (sig); return
                end if
            end if
            ! (2) strongest inflection: Laplacian sign change with the largest jump
            bestmag = -1; bestpos = -1; prevlap = 0
            do t = blo(ax) + 1, bhi(ax) - 1
                lap = sig(t - 1) - 2*sig(t) + sig(t + 1)
                if (t > blo(ax) + 1) then
                    if (((lap < 0) .neqv. (prevlap < 0)) .and. abs(lap - prevlap) > bestmag .and. t - blo(ax) >= min_child &
                        & .and. bhi(ax) - t + 1 >= min_child) then
                        bestmag = abs(lap - prevlap); bestpos = t
                    end if
                end if
                prevlap = lap
            end do
            if (bestpos > 0) then
                sax = ax; spos = bestpos; ok = .true.; deallocate (sig); return
            end if
        end do
        deallocate (sig)

    end subroutine s_amr_find_split

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

    !> Clip a candidate regrid box (global indices) so it does not overlap any acoustic source support bbox: per overlapping source,
    !! remove the overlap along the single axis/side that keeps the largest remaining extent (deterministic: lower axis, then begin
    !! side, wins ties). Clipping only shrinks the box and may empty it (hi < lo); the caller drops empties.
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
    !! the coarse RK update is windowed, so a reflux correction at a face cell outside the window would be silently dropped
    !! (conservation leak), and the coarse RHS only computes fluxes inside the window. The window only ever GROWS (s_grow_active_box
    !! is monotone, self-disabling at full domain), so containment established at init and re-established at each regrid holds
    !! between them. Collective (same window and block metadata on all ranks).
    impure subroutine s_amr_check_active_box_containment()

        integer :: k
        logical :: ok

        ! ab_active is only ever true at num_procs == 1 (m_active_box disables itself under
        ! MPI), so ab and block indices share the same (global == local) index space

        if ((.not. amr) .or. (.not. ab_active)) return
        do k = 1, amr_num_blocks
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

    !> Rank-invariant SPARSE tag list for level clustering (SP7a): all-gathers each rank's tagged-cell global linear indices, then
    !! decodes them into a coordinate list tags(1:3, 1:ntag). Replaces the O(global-grid) dense tag field entirely, so per-rank
    !! memory scales with the number of tagged cells. At np=1 the list is built directly from tag_grid (no allgather).
    !! Deterministic: every rank decodes the same gathered index set into the same list, so the bisection is rank-invariant.
    impure subroutine s_amr_union_gtag(tags, ntag, tag_grid, mg, ng, pg, sidx)

        integer, allocatable, intent(out) :: tags(:,:)
        integer, intent(out)              :: ntag
        logical, intent(in)               :: tag_grid(0:,0:,0:)
        integer, intent(in)               :: mg, ng, pg, sidx(3)
        integer                           :: ci, cj, ck, gi, gj, gk

#ifdef MFC_MPI
        integer                 :: i, jrem, nloc, ierr
        integer, allocatable    :: rcnt(:), rdsp(:)
        integer(8), allocatable :: locidx(:), allidx(:)

        if (num_procs > 1) then
            allocate (locidx((m + 1)*(n + 1)*(p + 1)), rcnt(num_procs), rdsp(num_procs))
            nloc = 0
            do ck = 0, p; do cj = 0, n; do ci = 0, m
                if (tag_grid(ci, cj, ck)) then
                    gi = ci + sidx(1); gj = 0; gk = 0
                    if (n_glb > 0) gj = cj + sidx(2)
                    if (p_glb > 0) gk = ck + sidx(3)
                    nloc = nloc + 1
                    locidx(nloc) = int(gi, 8) + int(mg + 1, 8)*(int(gj, 8) + int(ng + 1, 8)*int(gk, 8))
                end if
            end do; end do; end do
            call MPI_ALLGATHER(nloc, 1, MPI_INTEGER, rcnt, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
            rdsp(1) = 0
            do i = 2, num_procs; rdsp(i) = rdsp(i - 1) + rcnt(i - 1); end do
            ntag = rdsp(num_procs) + rcnt(num_procs)
            allocate (allidx(max(ntag, 1)), tags(3, max(ntag, 1)))
            call MPI_ALLGATHERV(locidx, nloc, MPI_INTEGER8, allidx, rcnt, rdsp, MPI_INTEGER8, MPI_COMM_WORLD, ierr)
            do i = 1, ntag
                gk = int(allidx(i)/(int(mg + 1, 8)*int(ng + 1, 8)))
                jrem = int(allidx(i) - int(gk, 8)*int(mg + 1, 8)*int(ng + 1, 8))
                gj = jrem/(mg + 1)
                gi = jrem - gj*(mg + 1)
                tags(1, i) = gi; tags(2, i) = gj; tags(3, i) = gk
            end do
            deallocate (locidx, rcnt, rdsp, allidx)
            return
        end if
#endif
        ! np=1: the whole global tag field is local; decode tag_grid directly into the list
        ntag = 0
        do ck = 0, p; do cj = 0, n; do ci = 0, m
            if (tag_grid(ci, cj, ck)) ntag = ntag + 1
        end do; end do; end do
        allocate (tags(3, max(ntag, 1)))
        ntag = 0
        do ck = 0, p; do cj = 0, n; do ci = 0, m
            if (tag_grid(ci, cj, ck)) then
                gi = ci + sidx(1); gj = 0; gk = 0
                if (n_glb > 0) gj = cj + sidx(2)
                if (p_glb > 0) gk = ck + sidx(3)
                ntag = ntag + 1
                tags(1, ntag) = gi; tags(2, ntag) = gj; tags(3, ntag) = gk
            end if
        end do; end do; end do

    end subroutine s_amr_union_gtag

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
    !! arrays sidx(:) (int8 linear index) / skb(:) (parent box id). The int8 encode matches the decode in s_amr_regrid's pass 2, so
    !! gathering these pairs across ranks and setting them into a per-parent dense window reproduces the old per-parent dense-window
    !! dedup (replicated/overlapping tags collapse) and the (k,j,i) extraction order exactly -> byte-identical child boxes. Batching
    !! all parents of a level into one allgatherv (in the caller) drops the collective count from O(#parent-boxes) to O(#levels).
    !! gwin is read here (not modified).
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

    !> Cluster a rank-invariant SPARSE tag list (global level-0 cell coords, tags(1:3, 1:ntag_in)) into a LIST of separated block
    !! boxes, identically on every rank. The caller builds the list (s_amr_union_gtag / s_amr_pack_gwin_pairs). Per-rank memory is
    !! O(#tagged), not O(global grid). Runs Berger-Rigoutsos recursive bisection until each box's tag efficiency reaches
    !! amr_cluster_eff (or it is atomic / the amr_max_blocks cap is reached), then merges any two boxes whose amr_buf-padded extents
    !! come within buff_size (guaranteeing no fine-fine adjacency: separated boxes stay >= buff_size apart, nearby ones collapse to
    !! a single box == the legacy bounding box). Boxes are the raw tagged extents; the caller pads, clamps and size-caps each one.
    impure subroutine s_amr_cluster(tags, ntag_in, boxes, nboxes)

        integer, intent(in)                   :: tags(:,:), ntag_in
        type(t_box), allocatable, intent(out) :: boxes(:)
        integer, intent(out)                  :: nboxes
        integer, allocatable                  :: slo(:,:), shi(:,:), alo(:,:), ahi(:,:)
        integer, allocatable                  :: sts(:), ste(:), wt(:,:)
        integer                               :: mg, ng, pg, t
        integer                               :: cap, nwork, nacc, i, j, d, sax, spos, thr, ntag, vol
        integer                               :: blo(3), bhi(3), ts, te, lo, hi, tmp(3)
        logical                               :: ok, force, capped, changed, tooclose
        real(wp)                              :: eff

        nboxes = 0
        if (ntag_in == 0) return
        mg = m_glb; ng = 0; pg = 0
        if (n_glb > 0) ng = n_glb
        if (p_glb > 0) pg = p_glb

        cap = amr_max_blocks
        allocate (slo(3, 4*cap + 8), shi(3, 4*cap + 8), alo(3, cap), ahi(3, cap))
        allocate (sts(4*cap + 8), ste(4*cap + 8), wt(3, ntag_in))
        ! working copy of the tag list, partitioned in place as the tree descends so each node scans only its tags
        do t = 1, ntag_in
            wt(:,t) = tags(:,t)
        end do
        nwork = 1; slo(:,1) = [0, 0, 0]; shi(:,1) = [mg, ng, pg]  ! first pop trims to the global tagged bbox
        sts(1) = 1; ste(1) = ntag_in
        nacc = 0; capped = .false.
        do while (nwork > 0)
            blo = slo(:,nwork); bhi = shi(:,nwork); ts = sts(nwork); te = ste(nwork); nwork = nwork - 1
            call s_amr_trim_box(wt, ts, te, blo, bhi, ok)
            if (.not. ok) cycle
            ! invariant: [ts:te] holds exactly the tags in this box, and trim only shrinks to their bbox => count is the range size
            ntag = te - ts + 1
            vol = 1
            do d = 1, num_dims; vol = vol*(bhi(d) - blo(d) + 1); end do
            eff = real(ntag, wp)/real(max(vol, 1), wp)
            call s_amr_find_split(wt, ts, te, blo, bhi, sax, spos, ok)
            force = (nacc + nwork + 1 >= cap)  ! splitting now could overflow the amr_max_blocks cap
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
                slo(:,nwork + 1) = blo; shi(:,nwork + 1) = bhi; shi(sax, nwork + 1) = spos - 1
                sts(nwork + 1) = ts; ste(nwork + 1) = lo - 1
                slo(:,nwork + 2) = blo; shi(:,nwork + 2) = bhi; slo(sax, nwork + 2) = spos
                sts(nwork + 2) = lo; ste(nwork + 2) = te
                nwork = nwork + 2
            end if
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
                        alo(:,j) = alo(:,nacc); ahi(:,j) = ahi(:,nacc)
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
        deallocate (slo, shi, alo, ahi, sts, ste, wt)

    end subroutine s_amr_cluster

    !> Regrid: tag by relative density gradient into a per-cell field, cluster (Berger-Rigoutsos + min-separation merge) into a list
    !! of separated boxes, pad/clamp/size-cap each, and rebuild every active slot. Each new box's slot prolongs from coarse then
    !! overwrites its overlap with whichever OLD slot(s) covered it (rank-local by construction; a split copies from one old slot, a
    !! merge from both). Called between steps only. No-op if nothing is tagged or the box set is unchanged.
    impure subroutine s_amr_regrid(q_cons_base)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_base
        logical, allocatable                                   :: tag_grid(:,:,:)
        type(t_box), allocatable                               :: boxes(:)
        integer                                                :: sidx(3), nboxes, box_level(amr_max_blocks)
        integer                                                :: old_np
        integer                                                :: old_ilo(3, amr_max_blocks), old_ext(3, amr_max_blocks)
        integer                                                :: old_level(amr_max_blocks)
        logical                                                :: old_owns(amr_max_blocks), same
        integer                                                :: i

        ! valid coarse CONS ghosts at internal rank boundaries: the tag sweep reads +/-1 across seams and the rebuild
        ! prolongation
        ! reads past the new intersection (ALL ranks call: pairwise per-direction exchange; complete no-op at np=1).

        call s_amr_exchange_coarse_cons_halo(q_cons_base)
        do i = 1, sys_size
            $:GPU_UPDATE(host='[q_cons_base(i)%sf]')
        end do

        ! Lagrangian-cloud exclusion bbox for this regrid (collective): smearing (mapCells) +
        ! stencil headroom (2) + drift margin until the next regrid (amr_buf)
        if (bubbles_lagrange) call s_amr_compute_lag_supp(mapCells + 2 + amr_buf)

        call s_amr_regrid_tag_cells(q_cons_base, tag_grid, sidx)
        call s_amr_regrid_cluster_tags(tag_grid, sidx, boxes, nboxes)
        if (nboxes == 0) return  ! nothing tagged on any rank; keep the current blocks
        call s_amr_regrid_shape_boxes(boxes, nboxes)
        if (nboxes == 0) return  ! every box was confined to the domain margin
        call s_amr_regrid_nest_children(boxes, nboxes, box_level)
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

        ! 1) per-cell tag field (density-gradient criterion, unchanged), skipping the two global boundary cells per active dim

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
                    ! total density gradient (sum of the continuity variables): degenerates to the single-fluid tagger
                    ! and is
                    ! immune to trace-fluid noise. Matched-density composition-only interfaces are invisible (documented
                    ! limit).
                    r0 = max(abs(f_amr_rho_tot(q_cons_base, ci, cj, ck)), 1.e-30_wp)
                    g = abs(f_amr_rho_tot(q_cons_base, ci + 1, cj, ck) - f_amr_rho_tot(q_cons_base, ci - 1, cj, ck))
                    if (n_glb > 0) g = max(g, abs(f_amr_rho_tot(q_cons_base, ci, cj + 1, ck) - f_amr_rho_tot(q_cons_base, ci, &
                        & cj - 1, ck)))
                    if (p_glb > 0) g = max(g, abs(f_amr_rho_tot(q_cons_base, ci, cj, ck + 1) - f_amr_rho_tot(q_cons_base, ci, cj, &
                        & ck - 1)))
                    ! 2*r0 normalizes the 2-cell central difference (rho at i+1..i-1); the 2 is the stencil span, NOT amr_ref_ratio
                    if (g/(2._wp*r0) > amr_tag_eps) tag_grid(ci, cj, ck) = .true.
                    ! the acoustic source support stays coarse (its spatials are coarse cell
                    ! indices): suppress tags there so the clusterer splits around the source
                    if (acoustic_source .and. tag_grid(ci, cj, ck)) then
                        if (f_in_acoustic_support(ci + sidx(1), cj + sidx(2), ck + sidx(3))) tag_grid(ci, cj, ck) = .false.
                    end if
                    ! the Lagrangian bubble cloud stays coarse (two-way coupling lives on the
                    ! coarse grid): suppress tags over its padded bbox
                    if (bubbles_lagrange .and. tag_grid(ci, cj, ck)) then
                        if (f_in_lag_support(ci + sidx(1), cj + sidx(2), ck + sidx(3))) tag_grid(ci, cj, ck) = .false.
                    end if
                end do
            end do
        end do

    end subroutine s_amr_regrid_tag_cells

    !> Regrid phase 2: union the per-rank tag fields into the rank-invariant sparse global tag list, then cluster it into a list of
    !! separated candidate boxes (nboxes = 0 if nothing is tagged on any rank).
    impure subroutine s_amr_regrid_cluster_tags(tag_grid, sidx, boxes, nboxes)

        logical, allocatable, intent(inout)   :: tag_grid(:,:,:)
        integer, intent(in)                   :: sidx(3)
        type(t_box), allocatable, intent(out) :: boxes(:)
        integer, intent(out)                  :: nboxes
        integer, allocatable                  :: tags(:,:)
        integer                               :: mg0, ng0, pg0, ntag

        ! 2) build the rank-invariant sparse global tag list, then cluster into a list of separated boxes

        mg0 = m_glb; ng0 = 0; pg0 = 0
        if (n_glb > 0) ng0 = n_glb
        if (p_glb > 0) pg0 = p_glb
        call s_amr_union_gtag(tags, ntag, tag_grid, mg0, ng0, pg0, sidx)
        deallocate (tag_grid)
        call s_amr_cluster(tags, ntag, boxes, nboxes)
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
            ! IB keeps the size-cap CLAMP (a body needs one contiguous block; splitting a body across tiles is
            ! untested); the
            ! general path leaves boxes full-size and TILES them (below) into <= amr_maxc_fit sub-blocks with a
            ! fine-fine halo
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
            ! keep candidate boxes clear of every acoustic source support (the source acts on the
            ! coarse grid only); clipping only shrinks, so boxes stay disjoint - empties drop below
            if (acoustic_source) call s_amr_clip_box_from_sources(lo, hi)
            if (bubbles_lagrange .and. lag_supp_on) call s_amr_clip_box_from_supp(lo, hi, lag_supp_lo, lag_supp_hi)
            ! active_box: boxes stay strictly inside the active window (the windowed coarse
            ! update would drop reflux corrections at faces outside it). Tags cannot arise
            ! outside (frozen-ambient exterior), so only the amr_buf padding is ever cut -
            ! and the cut cells are ambient. np=1 only (ab_active is false under MPI).
            if (ab_active) then
                lo(1) = max(lo(1), ab_x%beg + 1); hi(1) = min(hi(1), ab_x%end - 1)
                if (n_glb > 0) then; lo(2) = max(lo(2), ab_y%beg + 1); hi(2) = min(hi(2), ab_y%end - 1); end if
                if (p_glb > 0) then; lo(3) = max(lo(3), ab_z%beg + 1); hi(3) = min(hi(3), ab_z%end - 1); end if
            end if
            ! a fine block that PARTIALLY covers an immersed body is an untested regime (ghost
            ! prolongation through body-interior cells, refluxing across the body): any box that
            ! overlaps a body's bounding box is expanded to contain the whole body plus margin
            if (ib) call s_amr_expand_box_over_bodies(lo, hi)
            if (hi(1) < lo(1) .or. hi(2) < lo(2) .or. hi(3) < lo(3)) cycle  ! confined to the domain margin
            k = k + 1; boxes(k)%lo = lo; boxes(k)%hi = hi
        end do
        nboxes = k
        if (nboxes == 0) return

        ! max_grid_size tiling (non-IB): split any box larger than amr_maxc_fit into contiguous <= amr_maxc_fit
        ! sub-blocks so a
        ! whole block fits a rank's local solver scratch. Tiles are adjacent; the block-to-block fine-fine halo
        ! (s_amr_fine_
        ! fine_halo) makes the seams conservative and the reflux skips fine-fine faces. (IB keeps the clamp - see
        ! above.)
        if (.not. ib) then
            block
                type(t_box), allocatable :: tiled(:)
                integer                  :: kk2, ntl, capt
                allocate (tiled(amr_max_blocks))
                ntl = 0; capt = 0
                do kk2 = 1, nboxes
                    call s_amr_tile_box(boxes(kk2)%lo, boxes(kk2)%hi, tiled, ntl, amr_max_blocks, capt)
                end do
                if (capt == 1 .and. proc_rank == 0) print '(A,I0)', ' [amr] WARNING: tiling capped at amr_max_blocks = ', &
                    & amr_max_blocks
                deallocate (boxes); call move_alloc(tiled, boxes)
                nboxes = ntl
            end block
        end if

        if (ib) then
            ! body-containment expansion can make boxes overlap (bisection guaranteed disjoint
            ! boxes; two boxes near one body both grow over it): merge overlapping pairs to a
            ! bounding box until none remain - overlapping blocks would double-restrict/reflux
            merged = .true.
            do while (merged)
                merged = .false.
                outer: do k = 1, nboxes - 1
                    do kk = k + 1, nboxes
                        if (boxes(k)%lo(1) <= boxes(kk)%hi(1) .and. boxes(k)%hi(1) >= boxes(kk)%lo(1) .and. (n_glb == 0 &
                            & .or. (boxes(k)%lo(2) <= boxes(kk)%hi(2) .and. boxes(k)%hi(2) >= boxes(kk)%lo(2))) .and. (p_glb == 0 &
                            & .or. (boxes(k)%lo(3) <= boxes(kk)%hi(3) .and. boxes(k)%hi(3) >= boxes(kk)%lo(3)))) then
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
            ! the expansion may also have grown a box onto an acoustic source support or the
            ! Lagrangian cloud: the constraints (contain the body, exclude the source/cloud)
            ! cannot both hold - fail closed
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

        ! 3b) multi-level nesting: hierarchically append a box at level l nested inside each level-(l-1) box, for l =
        ! 2..amr_max_
        ! level. Parents-first ordering (every level-(l-1) box precedes its level-l children) so the build loop fills a
        ! parent
        ! before its child's gather-from-parent reads it. SENSOR-ON-FINE: each child's extent is the density-gradient
        ! sensor run
        ! on the parent-level FINE solution (the still-live OLD level-(l-1) blocks, read here BEFORE the step-5 stash),
        ! coarsened
        ! to L0-cell granularity and clustered - so children track features inside the parent instead of a fixed centre.
        ! A
        ! brand-new region with no old fine data falls back to a centred inset (the sensor takes over next regrid); a
        ! parent
        ! whose
        ! fine solution is smooth gets no child. Tagging only places boxes - conservation (restrict/reflux) is
        ! independent of
        ! where they sit. np=1 + non-IB (multi-level distribution / IB nesting are future work). Regions stay in L0 cell
        ! indices.

        box_level(1:nboxes) = 1
        if (amr_max_level >= 2) then
            ! the nesting loop below APPENDS level-l child boxes into `boxes` (up to amr_max_blocks total). The non-IB
            ! path already grew `boxes` to amr_max_blocks via the tiling move_alloc; the IB path (which only merges,
            ! never
            ! grows) leaves `boxes` at the cluster count, so grow it here or the child appends overrun the allocation.
            if (size(boxes) < amr_max_blocks) then
                block
                    type(t_box), allocatable :: grown(:)
                    allocate (grown(amr_max_blocks))
                    grown(1:nboxes) = boxes(1:nboxes)
                    call move_alloc(grown, boxes)
                end block
            end if
            block
                integer                  :: kb, ins(3), clo(3), chi(3), lev, plo, phi, newlo, ob, obi, ncb, kc, mlo(3), mhi(3)
                integer                  :: mg, ng, pg, nct, np_lev, nloc_send, gi, gj, gk, jrem, ntot_g
                integer, allocatable     :: ctags(:,:), skb(:), gkb(:)
                integer(8), allocatable  :: sidx(:), gidx(:)
                logical, allocatable     :: gwin(:,:,:), covered(:)
                integer, allocatable     :: mlo_all(:,:), mhi_all(:,:)
                logical                  :: any_tag
                type(t_box), allocatable :: cboxes(:)
#ifdef MFC_MPI
                integer              :: ierr, ip
                integer, allocatable :: rcnt(:), rdsp(:)
#endif

                ! host-refresh the live (old) blocks' continuity fields: the fine sensor below reads
                ! amr_slots(ob)%q_cons on the
                ! host, but the GPU_UPDATE host that the step-5 stash does runs AFTER this nesting - so the host copy is
                ! stale
                ! here
                do ob = 1, amr_num_blocks
                    if (.not. amr_owns_all(ob)) cycle  ! np>1: only the owner holds this old block's fine q_cons
                    do obi = eqn_idx%cont%beg, eqn_idx%cont%end
                        $:GPU_UPDATE(host='[amr_slots(ob)%q_cons(obi)%sf]')
                    end do
                end do
                ! Fine-sensor tags accumulate in a GLOBAL L0 frame: at np>1 an old block is read only by its owner, but
                ! its
                ! tag footprint can fall in ANOTHER rank's subdomain. Each parent's nesting window [mlo:mhi] is small vs
                ! the global grid, so a WINDOW-LOCAL dense field gwin (allocated per parent below) holds each owner's
                ! tags; s_amr_pack_gwin_pairs extracts them as (linear-index, kb) pairs, one per-level allgatherv unions all
                ! parents' pairs across ranks, and pass 2 rebuilds each parent's window from them (no O(global-grid) tag
                ! field, no local slice; the clusterer consumes the sparse per-parent list directly).
                mg = m_glb; ng = 0; pg = 0
                if (n_glb > 0) ng = n_glb
                if (p_glb > 0) pg = p_glb

                plo = 1; phi = nboxes  ! [plo:phi] = the boxes at the previous level (lev-1) to nest inside
                do lev = 2, amr_max_level
                    newlo = nboxes + 1
                    ! COLLECT -> ONE COMMUNICATE -> PROCESS, per level: the per-parent cross-rank union (the old per-(lev,kb)
                    ! allgather) is batched into a SINGLE allgatherv per level, so the collective count is O(#levels) not
                    ! O(#parent-boxes). Pass 1 tags each parent's window from OWNED obs and appends this rank's tagged cells as
                    ! (linear-index, parent-kb) pairs; one allgatherv unions them; Pass 2 rebuilds each parent's dense window from
                    ! the gathered pairs whose gkb==kb, which reproduces the old per-parent dense-window dedup and the (k,j,i)
                    ! extraction order exactly -> each parent's ctags set (and thus its child boxes) is byte-identical.
                    np_lev = phi - plo + 1
                    if (np_lev < 1) exit  ! nothing nested at the previous level -> no deeper levels possible
                    allocate (covered(plo:phi), mlo_all(3,plo:phi), mhi_all(3,plo:phi))
                    covered = .false.
                    nloc_send = 0
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

                        ! sensor-on-fine: tag from every OLD level-(lev-1) block overlapping this parent window
                        ! (amr_block_level
                        ! still holds the old levels here - it is reset to box_level at step 5b, below)
                        allocate (gwin(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3)))
                        gwin = .false.; any_tag = .false.
                        do ob = 1, amr_num_blocks
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
                            covered(kb) = .true.  ! replicated (metadata) - identical on every rank regardless of ownership
                            if (amr_owns_all(ob)) call s_amr_tag_child_from_fine(ob, mlo, mhi, gwin, any_tag)
                        end do
                        ! IB: always refine the body region at this level, even where the density sensor is quiet - mark
                        ! the
                        ! body's L0-frame bbox into gwin so it is clustered into a child (mirrors the L1 expand at
                        ! :3836).
                        ! Containment margin = max(amr_buf, 4) + amr_cpat_mar: the child window (mlo:mhi) is the parent
                        ! inset by
                        ! amr_cpat_mar, and clamping the tag to that window can eat up to amr_cpat_mar of the body's
                        ! stencil
                        ! margin at the parent-adjacent side. The parent (widened in s_amr_expand_box_over_bodies by
                        ! (amr_max_level-1)*amr_cpat_mar) now clears the body by enough that this window contains the
                        ! body plus
                        ! max(amr_buf, 4), so the tag survives the inset with a full image-point stencil of fluid on
                        ! every side:
                        ! the body SURFACE is refined at every level and the C/F boundary sits a full stencil off it, in
                        ! fluid.
                        if (ib) then
                            block
                                integer :: ib_i, bb_lo(3), bb_hi(3), gii, gjj, gkk
                                do ib_i = 1, num_ibs
                                    call s_amr_body_bbox(ib_i, max(amr_buf, 4) + amr_cpat_mar, bb_lo, bb_hi)
                                    ! clamp the body bbox to this parent's nesting window (global L0 frame -
                                    ! s_amr_body_bbox
                                    ! returns GLOBAL L0 cell indices, same frame as mlo/mhi)
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
                        ! extract THIS rank's OWNED tagged cells as (linear-index, kb) pairs into the per-level send
                        ! arrays. The int8 linear index matches the decode in pass 2, so the gathered pairs reproduce the
                        ! same window coords. gwin is read here, then freed.
                        call s_amr_pack_gwin_pairs(gwin, mlo, mhi, mg, ng, kb, sidx, skb, nloc_send)
                        deallocate (gwin)
                    end do

                    ! COMMUNICATE: one allgatherv per level (np>1)
                    if (.not. allocated(sidx)) then
                        allocate (sidx(0), skb(0))  ! this rank owned no tags at this level
                    end if
#ifdef MFC_MPI
                    if (num_procs > 1) then
                        allocate (rcnt(num_procs), rdsp(num_procs))
                        call MPI_ALLGATHER(nloc_send, 1, MPI_INTEGER, rcnt, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
                        rdsp(1) = 0
                        do ip = 2, num_procs
                            rdsp(ip) = rdsp(ip - 1) + rcnt(ip - 1)
                        end do
                        ntot_g = rdsp(num_procs) + rcnt(num_procs)
                        allocate (gidx(max(ntot_g, 1)), gkb(max(ntot_g, 1)))
                        call MPI_ALLGATHERV(sidx, nloc_send, MPI_INTEGER8, gidx, rcnt, rdsp, MPI_INTEGER8, MPI_COMM_WORLD, ierr)
                        call MPI_ALLGATHERV(skb, nloc_send, MPI_INTEGER, gkb, rcnt, rdsp, MPI_INTEGER, MPI_COMM_WORLD, ierr)
                        deallocate (rcnt, rdsp)
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

                    ! Pass 2: process (no comm)
                    do kb = plo, phi
                        if (nboxes + 1 > amr_max_blocks) exit  ! pool full - stop nesting
                        mlo = mlo_all(:,kb); mhi = mhi_all(:,kb)
                        if (mhi(1) < mlo(1)) cycle  ! too small to nest a child in x
                        if (n_glb > 0 .and. mhi(2) < mlo(2)) cycle
                        if (p_glb > 0 .and. mhi(3) < mlo(3)) cycle

                        ! rebuild this parent's dense window from the gathered pairs whose gkb==kb: setting .true. once per
                        ! gathered cell reproduces the old per-parent dedup (replicated/overlapping tags collapse), and the
                        ! (k,j,i) sparse extract below matches the old scan order -> byte-identical ctags.
                        allocate (gwin(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3)))
                        gwin = .false.
                        do i = 1, ntot_g
                            if (gkb(i) /= kb) cycle
                            gk = int(gidx(i)/(int(mg + 1, 8)*int(ng + 1, 8)))
                            jrem = int(gidx(i) - int(gk, 8)*int(mg + 1, 8)*int(ng + 1, 8))
                            gj = jrem/(mg + 1)
                            gi = jrem - gj*(mg + 1)
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
                            call s_amr_cluster(ctags, nct, cboxes, ncb)
                            deallocate (ctags)
                            do kc = 1, ncb
                                if (nboxes + 1 > amr_max_blocks) exit
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
                                ! IB: a child clustered from the (widened) body tag must fully contain every overlapping
                                ! body -
                                ! expand over bodies (mirrors the L1 expand at :3836), then re-clamp to the nesting
                                ! window so
                                ! the
                                ! child stays nested. Because the parent was widened by (amr_max_level-1)*amr_cpat_mar,
                                ! its
                                ! nesting window (mlo:mhi) already contains the body plus max(amr_buf, 4), so the
                                ! re-clamp does
                                ! NOT cut the body's stencil: the child CONTAINS the body bbox and the C/F boundary
                                ! lands a full
                                ! image-point stencil off the surface, in fluid (surface refined, not just the
                                ! interior).
                                if (ib) then
                                    call s_amr_expand_box_over_bodies(clo, chi)
                                    clo(1) = max(clo(1), mlo(1)); chi(1) = min(chi(1), mhi(1))
                                    if (n_glb > 0) then; clo(2) = max(clo(2), mlo(2)); chi(2) = min(chi(2), mhi(2)); end if
                                    if (p_glb > 0) then; clo(3) = max(clo(3), mlo(3)); chi(3) = min(chi(3), mhi(3)); end if
                                end if
                                ! slot cap: a level->=2 block's fine grid spans 4*(its L0 extent) cells while the slot
                                ! holds
                                ! 2*amr_maxc_fit fine cells, so a child's L0 extent must be <= amr_maxc_fit/2. In
                                ! LOCK-STEP a
                                ! feature wider than that TILES into adjacent <= amr_maxc_fit/2 sub-blocks (like the L1
                                ! tiling):
                                ! the per-stage fine-fine halo (s_amr_fine_fine_halo, level-aware) matches the shared
                                ! seam flux
                                ! and the L2->L1 reflux skips those fine-fine faces. SUBCYCLE advances level-2 children
                                ! per-block
                                ! (s_amr_advance_children) with no L2-L2 halo, so it keeps ONE capped child (adjacent
                                ! tiles
                                ! would
                                ! leak at their seam there - transposing that path is future work); a wide feature is
                                ! under-
                                ! refined rather than non-conservative.
                                if (amr_subcycle) then
                                    chi(1) = min(chi(1), clo(1) + amr_maxc_fit(1)/2 - 1)
                                    if (n_glb > 0) chi(2) = min(chi(2), clo(2) + amr_maxc_fit(2)/2 - 1)
                                    if (p_glb > 0) chi(3) = min(chi(3), clo(3) + amr_maxc_fit(3)/2 - 1)
                                    nboxes = nboxes + 1
                                    boxes(nboxes)%lo = clo; boxes(nboxes)%hi = chi; box_level(nboxes) = lev
                                else
                                    block
                                        type(t_box) :: l2t(amr_max_blocks)
                                        integer     :: nl2, cpd, it
                                        nl2 = 0; cpd = 0
                                        call s_amr_tile_box(clo, chi, l2t, nl2, amr_max_blocks, cpd, amr_maxc_fit/2)
                                        do it = 1, nl2
                                            if (nboxes + 1 > amr_max_blocks) exit
                                            nboxes = nboxes + 1
                                            boxes(nboxes)%lo = l2t(it)%lo; boxes(nboxes)%hi = l2t(it)%hi
                                            box_level(nboxes) = lev
                                        end do
                                    end block
                                end if
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
                            nboxes = nboxes + 1
                            boxes(nboxes)%lo = clo; boxes(nboxes)%hi = chi; box_level(nboxes) = lev
                        end if
                    end do
                    deallocate (gidx, gkb, covered, mlo_all, mhi_all)  ! per-level scratch - freed every level (no leak)
                    plo = newlo; phi = nboxes  ! the boxes just appended are the parents for the next level
                    if (phi < plo) exit  ! nothing nested at this level -> no deeper levels possible
                end do
                if (nboxes >= amr_max_blocks .and. proc_rank == 0) print '(A)', &
                    & ' [amr] NOTE: block pool full during multi-level nesting; some boxes were not refined further'
            end block
        end if

    end subroutine s_amr_regrid_nest_children

    ! 4) unchanged? (same count, boxes AND levels as the live slots -> keep them; a rebuild would reproduce them
    ! exactly
    ! anyway). The level must be compared too: a box that keeps its coordinates but changes refinement level would
    ! otherwise slip through with a stale amr_block_level, corrupting the level-aware coupling.
    impure subroutine s_amr_regrid_boxes_unchanged(boxes, nboxes, box_level, same)

        type(t_box), intent(in) :: boxes(:)
        integer, intent(in)     :: nboxes, box_level(:)
        logical, intent(out)    :: same
        integer                 :: k

        same = .false.
        if (nboxes == amr_num_blocks) then
            same = .true.
            do k = 1, nboxes
                if (any(boxes(k)%lo /= amr_slots(k)%region%lo) .or. any(boxes(k)%hi /= amr_slots(k)%region%hi) .or. box_level(k) &
                    & /= amr_block_level(k)) same = .false.
            end do
        end if

    end subroutine s_amr_regrid_boxes_unchanged

    !> Regrid phase 5: stash every live slot's fine interior (dead-between-steps q_cons_stor bounce), record the old block set
    !! (old_*), commit the new regions/levels/owners, and migrate each stashed old block point-to-point to the ranks that now own an
    !! overlapping new block.
    impure subroutine s_amr_regrid_stash_migrate(boxes, nboxes, box_level, old_np, old_ilo, old_ext, old_level, old_owns)

        type(t_box), intent(in) :: boxes(:)
        integer, intent(in)     :: nboxes, box_level(:)
        integer, intent(out)    :: old_np, old_ilo(:,:), old_ext(:,:), old_level(:)
        logical, intent(out)    :: old_owns(:)
        integer                 :: old_chi(3, amr_max_blocks), old_owner(amr_max_blocks)
        integer                 :: k, i
        integer                 :: np_l  !< local mirror of old_np: an INTENT(OUT) dummy is not allowed in the
        !                                   BLOCK specification expressions below (F2018 restricted expressions)

        ! 5) stash every live slot's fine interior (dead-between-steps q_cons_stor bounce), keeping its old intersection origin

        old_np = amr_num_blocks
        np_l = old_np
        do k = 1, old_np
            ! GLOBAL block origin + extents (replicated, valid on every rank - not the owner-only isect), so the
            ! cross-rank
            ! migration below and the overlap-copy's index shift are correct even where this rank did not own the old
            ! block
            old_ilo(:,k) = amr_region_lo_all(:,k)
            old_chi(:,k) = amr_region_hi_all(:,k)  ! old COARSE hi (for the P2P migration overlap test below)
            ! fine extent = (2**level)*footprint - 1: a level-2 block is 4x its L0 footprint, so stashing/migrating it
            ! with the
            ! level-1 factor (2x) truncates half its fine cells. Level-1 blocks (2**1 = 2) are byte-identical to before.
            old_ext(1, k) = (amr_ref_ratio**amr_block_level(k))*(amr_region_hi_all(1, k) - amr_region_lo_all(1, k) + 1) - 1
            old_ext(2, k) = merge((amr_ref_ratio**amr_block_level(k))*(amr_region_hi_all(2, k) - amr_region_lo_all(2, &
                    & k) + 1) - 1, 0, n_glb > 0)
            old_ext(3, k) = merge((amr_ref_ratio**amr_block_level(k))*(amr_region_hi_all(3, k) - amr_region_lo_all(3, &
                    & k) + 1) - 1, 0, p_glb > 0)
            old_owner(k) = amr_block_owner(k)
            ! overlap-copy must match levels: an old L2's stash is in the 4x parent-fine frame
            old_level(k) = amr_block_level(k)
            old_owns(k) = amr_owns_all(k)
            if (old_owns(k)) then
                do i = 1, sys_size
                    $:GPU_UPDATE(host='[amr_slots(k)%q_cons(i)%sf]')
                end do
                do i = 1, sys_size
                    amr_slots(k)%q_cons_stor(i)%sf(0:old_ext(1, k),0:old_ext(2, k),0:old_ext(3, &
                              & k)) = amr_slots(k)%q_cons(i)%sf(0:old_ext(1, k),0:old_ext(2, k),0:old_ext(3, k))
                end do
                ! non-polytropic QBMM: the side-state bounces through pb/mv_stor exactly like
                ! q_cons (both stors are dead between steps)
                if (qbmm .and. .not. polytropic) then
                    $:GPU_UPDATE(host='[amr_slots(k)%pb_f%sf, amr_slots(k)%mv_f%sf]')
                    amr_slots(k)%pb_stor%sf(0:old_ext(1, k),0:old_ext(2, k),0:old_ext(3, k),:, &
                              & :) = amr_slots(k)%pb_f%sf(0:old_ext(1, k),0:old_ext(2, k),0:old_ext(3, k),:,:)
                    amr_slots(k)%mv_stor%sf(0:old_ext(1, k),0:old_ext(2, k),0:old_ext(3, k),:, &
                              & :) = amr_slots(k)%mv_f%sf(0:old_ext(1, k),0:old_ext(2, k),0:old_ext(3, k),:,:)
                end if
            end if
        end do
        ! coarse pb/mv host-current for the per-block re-prolongation below
        if (qbmm .and. .not. polytropic) then
            $:GPU_UPDATE(host='[pb_ts(1)%sf, mv_ts(1)%sf]')
        end if

        ! set the regions + assign owners BEFORE the migration (P2P needs the new owners) and before the owner-dependent
        ! geometry (else s_set_amr_fine_geometry sizes the whole-block owner from a stale amr_block_owner)
        amr_num_blocks = nboxes
        do k = 1, nboxes
            amr_region_lo_all(:,k) = boxes(k)%lo; amr_region_hi_all(:,k) = boxes(k)%hi
            ! box_level(k) is the refinement level assigned during the hierarchical nesting above (1 for the L0->L1
            ! boxes, l for
            ! a box nested at level l). Setting it every regrid resets a stale level when a slot is reused across
            ! levels.
            amr_block_level(k) = box_level(k)
        end do
        ! block set changed: dirty the cached seam-pair AND overlap-rank lists NOW - the rebuild's per-block P2P gathers
        ! (s_amr_regrid_rebuild_slots) consume the overlap lists with the NEW boxes, so flagging after them would be too late
        amr_seam_pairs_dirty = .true.
        ! Proper-nesting guard: each level>=2 block must be covered by EXACTLY ONE parent-level block.
        ! f_amr_parent_block (and
        ! the gather/reflux that key off it) take the FIRST overlap, so a fine tile straddling two parent tiles - an
        ! internal
        ! parent-level tile seam crossed by a nested feature - would silently couple to only one parent (wrong coarse BC
        ! + a
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

#ifdef MFC_MPI
        ! Cross-rank fine-state migration: the overlap-copy below preserves each covering old block's fine detail by
        ! reading
        ! amr_slots(kk)%q_cons_stor, but an old block may be owned by a rank OTHER than the one now owning a covering
        ! new block.
        ! POINT-TO-POINT (mirrors s_amr_gather_coarse_patch): each old owner sends its stashed fine state ONLY to the
        ! distinct
        ! new-block owners whose region overlaps that old block. A rank that did not receive old block kk never reads it
        ! - the
        ! overlap-copy's per-(k,kk) index guard skips every cell of a non-overlapping pair. No-op at np=1 (single owner,
        ! local).
        if (num_procs > 1) then
            block
                integer               :: kk, k2, ii, gi, gj, gk, idx2, ierr2, rr, maxcnt, nrq
                integer               :: cnt(np_l)
                logical               :: getk(np_l), isdest(0:num_procs - 1)
                real(wp), allocatable :: spack(:,:), rpack(:,:)
                integer, allocatable  :: rq(:)
                maxcnt = 0
                do kk = 1, old_np
                    cnt(kk) = sys_size*(old_ext(1, kk) + 1)*(old_ext(2, kk) + 1)*(old_ext(3, kk) + 1)
                    maxcnt = max(maxcnt, cnt(kk))
                    ! I need old block kk iff I own a NEW block overlapping it (and do not already hold kk locally)
                    getk(kk) = .false.
                    if (.not. old_owns(kk)) then
                        do k2 = 1, nboxes
                            if (amr_block_owner(k2) == proc_rank .and. f_amr_boxes_overlap(boxes(k2)%lo, boxes(k2)%hi, old_ilo(:, &
                                & kk), old_chi(:,kk))) then
                                getk(kk) = .true.; exit
                            end if
                        end do
                    end if
                end do
                ! a received old block needs a live slot to unpack its q_cons_stor into (freed by the reconcile below)
                do kk = 1, old_np
                    if (getk(kk)) call s_amr_alloc_slot(kk)
                end do
                allocate (rq(old_np*num_procs), spack(max(maxcnt, 1), old_np), rpack(max(maxcnt, 1), old_np))
                nrq = 0
                do kk = 1, old_np  ! post receives for the old blocks I need
                    if (.not. getk(kk)) cycle
                    nrq = nrq + 1
                    call MPI_IRECV(rpack(1, kk), cnt(kk), mpi_p, old_owner(kk), kk, MPI_COMM_WORLD, rq(nrq), ierr2)
                end do
                do kk = 1, old_np  ! pack + send each old block I own to every distinct new-owner (/= me) overlapping it
                    if (.not. old_owns(kk)) cycle
                    isdest = .false.
                    do k2 = 1, nboxes
                        rr = amr_block_owner(k2)
                        if (rr /= proc_rank .and. f_amr_boxes_overlap(boxes(k2)%lo, boxes(k2)%hi, old_ilo(:,kk), old_chi(:, &
                            & kk))) isdest(rr) = .true.
                    end do
                    if (.not. any(isdest)) cycle
                    idx2 = 0
                    do ii = 1, sys_size
                        do gk = 0, old_ext(3, kk)
                            do gj = 0, old_ext(2, kk)
                                do gi = 0, old_ext(1, kk)
                                    idx2 = idx2 + 1
                                    spack(idx2, kk) = real(amr_slots(kk)%q_cons_stor(ii)%sf(gi, gj, gk), wp)
                                end do
                            end do
                        end do
                    end do
                    do rr = 0, num_procs - 1
                        if (.not. isdest(rr)) cycle
                        nrq = nrq + 1
                        call MPI_ISEND(spack(1, kk), cnt(kk), mpi_p, rr, kk, MPI_COMM_WORLD, rq(nrq), ierr2)
                    end do
                end do
                if (nrq > 0) call MPI_WAITALL(nrq, rq, MPI_STATUSES_IGNORE, ierr2)
                do kk = 1, old_np  ! unpack the received old blocks into their replicated q_cons_stor slots
                    if (.not. getk(kk)) cycle
                    idx2 = 0
                    do ii = 1, sys_size
                        do gk = 0, old_ext(3, kk)
                            do gj = 0, old_ext(2, kk)
                                do gi = 0, old_ext(1, kk)
                                    idx2 = idx2 + 1
                                    amr_slots(kk)%q_cons_stor(ii)%sf(gi, gj, gk) = real(rpack(idx2, kk), stp)
                                end do
                            end do
                        end do
                    end do
                end do
                deallocate (rq, spack, rpack)
            end block
        end if
#endif

    end subroutine s_amr_regrid_stash_migrate

    !> Regrid phase 6: build each new slot - geometry (collective), prolong from coarse, overwrite the overlap from every covering
    !! stashed old block - then reconcile the slot pool, rebuild the fine IB state, and re-validate the seam topology.
    impure subroutine s_amr_regrid_rebuild_slots(q_cons_base, boxes, nboxes, old_np, old_ilo, old_ext, old_level, old_owns)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_base
        type(t_box), intent(in)                                :: boxes(:)
        integer, intent(in)                                    :: nboxes, old_np, old_ilo(:,:), old_ext(:,:), old_level(:)
        logical, intent(in)                                    :: old_owns(:)
        integer                                                :: sh(3), k, kk, i, fi, fj, fk, ofi, ofj, ofk
        logical                                                :: any_xchg

        ! 6) build each new slot: geometry (collective on all ranks), prolong, then overlap-copy from every covering old slot

        any_xchg = .false.
        do k = 1, nboxes
            amr_cur = k
            ! owned slot needs its arrays before geometry/prolong
            if (amr_block_owner(k) == proc_rank) call s_amr_alloc_slot(k)
            call s_set_amr_fine_geometry(boxes(k)%lo, boxes(k)%hi)
            any_xchg = any_xchg .or. amr_xchg_coarse_ghosts
            ! fine-level distribution: gather this new block's coarse patch (collective - before the owner-only cycle;
            ! q_cons_base is host-current with valid ghosts from the exchange at the top of s_amr_regrid)
            call s_amr_gather_coarse_patch(q_cons_base, .false.)
            ! non-polytropic QBMM: gather the coarse pb/mv patch too (ALL ranks - P2P; owners re-prolong from it below)
            if (qbmm .and. .not. polytropic) call s_amr_gather_coarse_patch_pbmv(pb_ts(1)%sf, mv_ts(1)%sf, .false.)
            if (.not. amr_rank_owns_block) cycle
            call s_interpolate_coarse_to_fine()
            ! every old block's stashed fine state is now replicated in amr_slots(kk)%q_cons_stor (migration above), so
            ! copy
            ! the overlap from EVERY covering old block regardless of who owned it - sh is the old->new LOCAL fine index
            ! shift.
            ! A level>=2 block SKIPS this: old_ilo/sh are the L0 index frame, but a child's amr_isect_lo is its
            ! PARENT-fine
            ! frame,
            ! so the shift is wrong. It re-prolongs from its (freshly-built, parents-first) parent each regrid instead;
            ! the
            ! coupling
            ! keeps conservation. Detail-preserving same-level L2 migration (parent-fine overlap) is a later increment.
            if (amr_block_level(amr_cur) < 2) then
                do kk = 1, old_np
                    ! same-level overlap only (a child's stash is 4x-framed)
                    if (old_level(kk) /= amr_block_level(amr_cur)) cycle
                    ! old LOCAL fine index = new LOCAL fine index + sh (collapsed dims sh=0)
                    sh = amr_ref_ratio*(amr_isect_lo - old_ilo(:,kk))
                    do i = 1, sys_size
                        do fk = 0, amr_slots(k)%p
                            ofk = fk + sh(3)
                            if (p_glb > 0 .and. (ofk < 0 .or. ofk > old_ext(3, kk))) cycle
                            do fj = 0, amr_slots(k)%n
                                ofj = fj + sh(2)
                                if (n_glb > 0 .and. (ofj < 0 .or. ofj > old_ext(2, kk))) cycle
                                do fi = 0, amr_slots(k)%m
                                    ofi = fi + sh(1)
                                    if (ofi < 0 .or. ofi > old_ext(1, kk)) cycle
                                    amr_slots(k)%q_cons(i)%sf(fi, fj, fk) = amr_slots(kk)%q_cons_stor(i)%sf(ofi, ofj, ofk)
                                end do
                            end do
                        end do
                    end do
                end do
            end if
            do i = 1, sys_size
                $:GPU_UPDATE(device='[amr_slots(k)%q_cons(i)%sf]')
            end do
            ! non-polytropic QBMM: prolong the side-state from coarse (piecewise-constant),
            ! then overwrite the overlap with the old blocks' fine data (same index shift)
            if (qbmm .and. .not. polytropic) then
                call s_amr_prolong_pbmv()
                ! level>=2 re-prolongs only (the L0-frame overlap shift is wrong for a child)
                if (amr_block_level(amr_cur) < 2) then
                    do kk = 1, old_np
                        if (old_level(kk) /= amr_block_level(amr_cur)) cycle  ! same-level overlap only
                        if (.not. old_owns(kk)) cycle
                        sh = amr_ref_ratio*(amr_isect_lo - old_ilo(:,kk))
                        do fk = 0, amr_slots(k)%p
                            ofk = fk + sh(3)
                            if (p_glb > 0 .and. (ofk < 0 .or. ofk > old_ext(3, kk))) cycle
                            do fj = 0, amr_slots(k)%n
                                ofj = fj + sh(2)
                                if (n_glb > 0 .and. (ofj < 0 .or. ofj > old_ext(2, kk))) cycle
                                do fi = 0, amr_slots(k)%m
                                    ofi = fi + sh(1)
                                    if (ofi < 0 .or. ofi > old_ext(1, kk)) cycle
                                    amr_slots(k)%pb_f%sf(fi, fj, fk,:,:) = amr_slots(kk)%pb_stor%sf(ofi, ofj, ofk,:,:)
                                    amr_slots(k)%mv_f%sf(fi, fj, fk,:,:) = amr_slots(kk)%mv_stor%sf(ofi, ofj, ofk,:,:)
                                end do
                            end do
                        end do
                    end do
                end if
                $:GPU_UPDATE(device='[amr_slots(k)%pb_f%sf, amr_slots(k)%mv_f%sf]')
            end if
            ! whole-block-per-rank: no fine-fine halo; the new block's ghost shell is (re)prolonged by the next fine advance
        end do
        amr_xchg_coarse_ghosts = any_xchg  ! coarse halo exchanged once per step if ANY block needs it
        ! lazy sizing: free the transient regrid slots (old blocks this rank stashed/received but does not now own); the
        ! new-owned slots were allocated in the build loop, so this only frees - a rank keeps just its owned blocks'
        ! fine arrays
        call s_amr_reconcile_slots()
        ! rebuild every block's fine-grid IB state for the NEW geometry (markers/ghost points/
        ! image points recomputed from the body definitions; no state carries across regrids)
        if (ib) call s_amr_setup_ib()
        call s_amr_select_slot(1)
        call s_amr_check_seam_topology()  ! abort on seam topologies no halo reconciles (silent leak otherwise)

    end subroutine s_amr_regrid_rebuild_slots

    !> Sensor-on-fine child tagging: OR-accumulate density-gradient tags from an OLD fine block's solution into an L0-cell tag grid,
    !! restricted to a parent nesting window. Reads amr_slots(ob)%q_cons on the HOST (the caller host-refreshes the cont range
    !! first; the step-5 stash's GPU_UPDATE runs later). Fine cell (fi,fj,fk) covers L0 cell (ci,cj,ck) with fi = rr*(ci-olo(1))+d
    !! etc.; the gradient uses one-sided differences at the fine-interior edges so no stale fine ghost is read. Only decides
    !! placement - conservation is enforced downstream by restrict/reflux regardless of the box extent.
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
                                r0 = max(abs(f_amr_rho_tot(amr_slots(ob)%q_cons, fi, fj, fk)), 1.e-30_wp)
                                g = abs(f_amr_rho_tot(amr_slots(ob)%q_cons, min(fi + 1, fm1), fj, &
                                        & fk) - f_amr_rho_tot(amr_slots(ob)%q_cons, max(fi - 1, 0), fj, fk))
                                if (n_glb > 0) g = max(g, abs(f_amr_rho_tot(amr_slots(ob)%q_cons, fi, min(fj + 1, fm2), &
                                    & fk) - f_amr_rho_tot(amr_slots(ob)%q_cons, fi, max(fj - 1, 0), fk)))
                                if (p_glb > 0) g = max(g, abs(f_amr_rho_tot(amr_slots(ob)%q_cons, fi, fj, min(fk + 1, &
                                    & fm3)) - f_amr_rho_tot(amr_slots(ob)%q_cons, fi, fj, max(fk - 1, 0))))
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

    !> Total density (sum of the continuity variables) at one cell: the regrid tag field. Reduces to variable 1 for one fluid.
    pure function f_amr_rho_tot(q, ci, cj, ck) result(r)

        type(scalar_field), dimension(:), intent(in) :: q
        integer, intent(in)                          :: ci, cj, ck
        real(wp)                                     :: r
        integer                                      :: f

        r = 0._wp
        do f = eqn_idx%cont%beg, eqn_idx%cont%end
            r = r + real(q(f)%sf(ci, cj, ck), wp)
        end do

    end function f_amr_rho_tot

end module m_amr_regrid
