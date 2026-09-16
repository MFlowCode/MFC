!>
!!@file
!!@brief Contains module m_amr_distribution

#! AMD OpenMP lane: assert allocatables present on every kernel here (see OMP_DEFAULT_STR). Every conditionally allocated
#! module array a kernel here names launches only under its allocation's own condition (amr_rvw: cyl_coord; sw_jac/jac: igr;
#! amr_cg_pb/mv: do_pbmv; amr_gst_a/b: amr_subcycle; amr_prim_st/amr_bt_*: amr_prim_batch); amr_cg and amr_cons_br/stor_st are
#! allocated before first use. A kernel naming an unallocated array aborts. Keep it so.
#:set MFC_OMP_PRESENT_ALLOCATABLE = True
#:include 'macros.fpp'

!> @brief Rank decomposition, box arithmetic and block ownership (SFC cut, owner lookups, owned-block lists).
module m_amr_distribution

#ifdef MFC_MPI
    use mpi  !< MPI-IO for the parallel_io AMR restart file
#endif

    use m_derived_types  ! scalar_field, t_box, int_bounds_info
    use m_box, only: f_morton  ! shared 3D Morton key (single-sourced with m_sfc_partition)
    use m_global_parameters
    use m_constants, only: num_fluids_max, model_eqns_6eq, mapCells, K_ib, K_pc, BC_GHOST_EXTRAP
    use m_pressure_relaxation, only: s_pressure_relaxation_procedure
    use m_mpi_proxy, only: s_mpi_abort
    use m_mpi_common, only: s_mpi_allreduce_integer_min, s_mpi_allreduce_integer_max, s_mpi_allreduce_sum, s_mpi_allreduce_min, &
        & s_mpi_allreduce_max, s_mpi_allreduce_integer_sum, s_mpi_sendrecv_variables_buffers, s_mpi_allreduce_array_max
    use m_rhs, only: s_compute_rhs, q_prim_qp
    use m_variables_conversion, only: s_convert_species_to_mixture_variables_kernel, s_compute_pressure, enforce_density_floor_vc
    use m_phase_change, only: s_infinite_relaxation_k, pc_iter_count
    use m_amr_registers, only: s_amr_zero_fine_registers, s_amr_reflux_apply_faces, s_amr_parent_foot, freg, creg, &
        & s_amr_reg_prepare, f_amr_face_is_seam
    use m_rank_timing, only: s_rank_time_tic, s_rank_time_toc
    use m_phase_timing
    use m_amr_xchg_audit  ! per-call-site accounting of every AMR p2p transfer (s_xa_rec + XA_* site ids)
    use m_ibm, only: s_ibm_alloc_fine, s_ibm_setup_fine, s_ibm_swap_to_fine, s_ibm_restore_from_fine, s_ibm_correct_state, &
        & s_ibm_load_fine_markers, s_update_mib, moving_immersed_boundary_flag, num_gps, ib_markers
    use m_hypoelastic, only: s_hypoelastic_update_fd_coeffs
    use m_weno, only: s_compute_weno_coefficients
    use m_active_box, only: ab_active
    use m_bubbles_EL, only: s_lag_cloud_bbox_local
    use m_igr, only: jac, jac_old
    use m_amr_state

    implicit none

    private
    public :: f_amr_m1_seq, f_amr_m1_tag, f_amr_overlap_count, f_amr_rank_overlaps, s_amr_assign_block_owners, s_amr_body_bbox, &
        & s_amr_box_isect, s_amr_compute_isect, s_amr_expand_box_over_bodies, s_amr_m1_wave_open, s_amr_rank_coarse_range, &
        & s_amr_rank_decomp, s_amr_rank_interior, s_amr_ranks_overlapping, s_amr_refresh_lists, s_amr_refresh_my_blocks, &
        & s_amr_sfc_cut, s_amr_tile_box, s_amr_validate_decomp, s_amr_validate_owner

contains

    !> Rank r's coarse-grid decomposition (start_idx + local extent m/n/p), computed O(1) from its cartesian coords instead of the
    !! replicated amr_decomp table. The domain cart is MPI_CART_CREATE(reorder=.false., dims=[num_procs_x,num_procs_y,num_procs_z]),
    !! so ranks keep MPI_COMM_WORLD order and r's coords are row-major in r. The per-dim split is the integer split-with-remainder
    !! on cell indices, identical to s_mpi_decompose_computational_domain (independent of grid stretching, which changes coords not
    !! indices).
    pure subroutine s_amr_rank_decomp(r, sidx, ext)

        integer, intent(in)  :: r
        integer, intent(out) :: sidx(3), ext(3)
        integer              :: coords(3), gd(3), pd(3), base, rem, d

        gd(1) = m_glb; gd(2) = n_glb; gd(3) = p_glb
        pd(1) = num_procs_x; pd(2) = num_procs_y; pd(3) = num_procs_z
        ! row-major: r = cx*(py*pz) + cy*pz + cz  (py=pz=1 when that dim is not decomposed)
        coords(1) = r/(pd(2)*pd(3))
        coords(2) = mod(r, pd(2)*pd(3))/pd(3)
        coords(3) = mod(r, pd(3))
        sidx = 0; ext = 0
        do d = 1, 3
            if (d == 2 .and. n_glb == 0) cycle
            if (d == 3 .and. p_glb == 0) cycle
            base = (gd(d) + 1)/pd(d)
            rem = mod(gd(d) + 1, pd(d))
            ext(d) = base - 1 + merge(1, 0, coords(d) < rem)
            sidx(d) = coords(d)*base + min(coords(d), rem)
        end do

    end subroutine s_amr_rank_decomp

    !> Guard: the computed accessor must reproduce this rank's actual decomposition. Every rank checks its own entry, so
    !! collectively the formula is checked for all r; fires immediately if it ever drifts from s_mpi_decompose_computational_domain.
    !! O(1), always on.
    impure subroutine s_amr_validate_decomp()

        integer :: sidx(3), ext(3)
        logical :: ok

        call s_amr_rank_decomp(proc_rank, sidx, ext)
        ! nested guards, not a single .and./.or.: Fortran does not short-circuit, so start_idx(2)/start_idx(3) (start_idx is sized
        ! num_dims) would be read out of bounds in 1D/2D even though the guard is false (a bounds-checked build aborts).
        ok = (sidx(1) == start_idx(1) .and. ext(1) == m)
        if (n_glb > 0) then
            if (sidx(2) /= start_idx(2) .or. ext(2) /= n) ok = .false.
        end if
        if (p_glb > 0) then
            if (sidx(3) /= start_idx(3) .or. ext(3) /= p) ok = .false.
        end if
        if (.not. ok) then
            call s_mpi_abort('s_amr_rank_decomp does not reproduce this rank''s decomposition - computed split disagrees with ' &
                             & // 's_mpi_decompose_computational_domain')
        end if

    end subroutine s_amr_validate_decomp

    !> Closed-form inverse of the per-dim split-with-remainder used by s_amr_rank_decomp: the cart coord owning global coarse cell
    !! g. base=(gd+1)/pd, rem=mod(gd+1,pd). Exact for g in [0,gd]; callers clamp to [0,pd-1] for out-of-range g (ghost reach).
    !! base==0 (more ranks than cells) => the g<rem*(base+1) branch always holds and c=g.
    pure integer function f_amr_cell_coord(g, base, rem) result(c)

        integer, intent(in) :: g, base, rem

        if (g < rem*(base + 1)) then
            c = g/(base + 1)
        else
            c = rem + (g - rem*(base + 1))/base
        end if

    end function f_amr_cell_coord

    !> Per-dim contiguous rank-coord range [clo:chi] whose owned coarse slab intersects box [blo:bhi], clamped to [0,pd-1]. The
    !! boundary coords' coarse_range is extended by buff_size exactly at the domain edge (s_amr_rank_coarse_range), so this clamped
    !! interior-frame range reproduces both the interior (scatter) and the coarse_range (gather) intersection sets: a box reaching
    !! the ghost zone clamps to the boundary coord whose extended slab contains it, and there is no rank beyond that coord.
    !! Collapsed dims (n_glb==0 / p_glb==0) contribute coord 0.
    pure subroutine s_amr_coord_range(blo, bhi, clo, chi)

        integer, intent(in)  :: blo(3), bhi(3)
        integer, intent(out) :: clo(3), chi(3)
        integer              :: gd(3), pd(3), base, rem, d

        gd(1) = m_glb; gd(2) = n_glb; gd(3) = p_glb
        pd(1) = num_procs_x; pd(2) = num_procs_y; pd(3) = num_procs_z
        clo = 0; chi = 0
        do d = 1, 3
            if (d == 2 .and. n_glb == 0) cycle
            if (d == 3 .and. p_glb == 0) cycle
            base = (gd(d) + 1)/pd(d)
            rem = mod(gd(d) + 1, pd(d))
            clo(d) = min(max(f_amr_cell_coord(blo(d), base, rem), 0), pd(d) - 1)
            chi(d) = min(max(f_amr_cell_coord(bhi(d), base, rem), 0), pd(d) - 1)
        end do

    end subroutine s_amr_coord_range

    !> Ascending rank list overlapping coarse box [blo:bhi]. Enumerates the coord brick cx->cy->cz so r = cx*(Py*Pz)+cy*Pz+cz is
    !! monotonic => ascending, reproducing the r=0..num_procs-1 scan order. Owner not excluded (consumers keep their own owner
    !! skip). Caller sizes ranks(:) >= f_amr_overlap_count(blo,bhi).
    pure subroutine s_amr_ranks_overlapping(blo, bhi, ranks, nr)

        integer, intent(in)  :: blo(3), bhi(3)
        integer, intent(out) :: ranks(:)
        integer, intent(out) :: nr
        integer              :: clo(3), chi(3), cx, cy, cz, pyz

        call s_amr_coord_range(blo, bhi, clo, chi)
        pyz = num_procs_y*num_procs_z
        nr = 0
        do cx = clo(1), chi(1)
            do cy = clo(2), chi(2)
                do cz = clo(3), chi(3)
                    nr = nr + 1
                    ranks(nr) = cx*pyz + cy*num_procs_z + cz
                end do
            end do
        end do

    end subroutine s_amr_ranks_overlapping

    !> Does rank r's subdomain overlap coarse box [blo:bhi]? Answers membership without enumerating the overlap set (which
    !! s_amr_ranks_overlapping must do, at O(num_procs) writes on a box spanning the machine); the clusterer asks this of every node
    !! it walks, the root included.
    pure logical function f_amr_rank_overlaps(blo, bhi, r) result(hit)

        integer, intent(in) :: blo(3), bhi(3), r
        integer             :: clo(3), chi(3), c(3)

        call s_amr_coord_range(blo, bhi, clo, chi)
        c(1) = r/(num_procs_y*num_procs_z)  ! same coord -> rank map the enumeration uses: r = cx*(Py*Pz) + cy*Pz + cz
        c(2) = mod(r/num_procs_z, num_procs_y)
        c(3) = mod(r, num_procs_z)
        hit = all(c >= clo) .and. all(c <= chi)

    end function f_amr_rank_overlaps

    !> Overlap count only (allocation sizing), = product of the per-dim coord-range widths.
    pure integer function f_amr_overlap_count(blo, bhi) result(nr)

        integer, intent(in) :: blo(3), bhi(3)
        integer             :: clo(3), chi(3)

        call s_amr_coord_range(blo, bhi, clo, chi)
        nr = (chi(1) - clo(1) + 1)*(chi(2) - clo(2) + 1)*(chi(3) - clo(3) + 1)

    end function f_amr_overlap_count

    !> Rank r's contiguous owned coarse-cell range per dim from the computed decomposition (s_amr_rank_decomp): interior
    !! [start:start+ext] plus its physical-boundary ghosts (buff_size cells only where the subdomain touches the domain edge). One
    !! contiguous span so box intersections identify contributors without a per-cell scan.
    pure subroutine s_amr_rank_coarse_range(r, crlo, crhi)

        integer, intent(in)  :: r
        integer, intent(out) :: crlo(3), crhi(3)
        integer              :: sidx(3), ext(3)

        call s_amr_rank_decomp(r, sidx, ext)
        crlo = 0; crhi = 0
        crlo(1) = sidx(1); if (sidx(1) == 0) crlo(1) = -buff_size
        crhi(1) = sidx(1) + ext(1); if (crhi(1) == m_glb) crhi(1) = crhi(1) + buff_size
        if (n_glb > 0) then
            crlo(2) = sidx(2); if (sidx(2) == 0) crlo(2) = -buff_size
            crhi(2) = sidx(2) + ext(2); if (crhi(2) == n_glb) crhi(2) = crhi(2) + buff_size
        end if
        if (p_glb > 0) then
            crlo(3) = sidx(3); if (sidx(3) == 0) crlo(3) = -buff_size
            crhi(3) = sidx(3) + ext(3); if (crhi(3) == p_glb) crhi(3) = crhi(3) + buff_size
        end if

    end subroutine s_amr_rank_coarse_range

    !> Per-dim intersection of two global boxes [alo:ahi] and [blo:bhi] -> [olo:ohi] (empty when olo > ohi in some dim).
    pure subroutine s_amr_box_isect(alo, ahi, blo, bhi, olo, ohi)

        integer, intent(in)  :: alo(3), ahi(3), blo(3), bhi(3)
        integer, intent(out) :: olo(3), ohi(3)

        olo = max(alo, blo); ohi = min(ahi, bhi)

    end subroutine s_amr_box_isect

    !> Per-block measured-cost weight over each block's level-0 footprint, replicated on every rank: each rank sums the load-weight
    !! cost model (base 1 + K_ib per IB-marked cell + K_pc per phase-change Newton iteration, when that diagnostic array is live)
    !! over its owned coarse cells inside the footprint, then one MPI_ALLREDUCE(SUM) makes the vector identical everywhere. No cost
    !! signals -> cost(k) = footprint cell count exactly (pure-geometry fallback). The Lagrangian cloud is excluded from blocks by
    !! construction, so K_bub never applies here; pc_iter_count is populated only when a load-weight diagnostic writer is on (enable
    !! load_weight_wrt to make the balance phase-change-aware), guarded by allocated().
    impure subroutine s_amr_block_cost(cost)

        real(wp), intent(out) :: cost(:)
        integer               :: k, j, kk, l, lo(3), hi(3)
        real(wp)              :: c

#ifdef MFC_MPI
        integer :: ierr
#endif

        ! one host refresh per regrid: both signal fields advance on the device between regrids
        if (ib) then
            $:GPU_UPDATE(host='[ib_markers%sf]')
        end if
        if (allocated(pc_iter_count)) then
            $:GPU_UPDATE(host='[pc_iter_count]')
        end if
        do k = 1, amr_num_blocks
            ! block footprint /\ this rank's coarse subdomain, in local interior indices (empty -> no-trip loops)
            lo = 0; hi = 0
            lo(1) = max(amr_region_lo_all(1, k) - start_idx(1), 0); hi(1) = min(amr_region_hi_all(1, k) - start_idx(1), m)
            if (n_glb > 0) then
                lo(2) = max(amr_region_lo_all(2, k) - start_idx(2), 0); hi(2) = min(amr_region_hi_all(2, k) - start_idx(2), n)
            end if
            if (p_glb > 0) then
                lo(3) = max(amr_region_lo_all(3, k) - start_idx(3), 0); hi(3) = min(amr_region_hi_all(3, k) - start_idx(3), p)
            end if
            c = 0._wp
            do l = lo(3), hi(3)
                do kk = lo(2), hi(2)
                    do j = lo(1), hi(1)
                        c = c + 1._wp
                        if (ib) then
                            if (ib_markers%sf(j, kk, l) /= 0) c = c + K_ib
                        end if
                        if (allocated(pc_iter_count)) c = c + K_pc*real(pc_iter_count(j, kk, l), wp)
                    end do
                end do
            end do
            cost(k) = c
        end do
#ifdef MFC_MPI
        call MPI_ALLREDUCE(MPI_IN_PLACE, cost, amr_num_blocks, mpi_p, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

    end subroutine s_amr_block_cost

    !> Cost-weighted SFC partition of n items into num_procs contiguous Morton-key ranges. Sorts item indices by ascending Morton
    !! key (stable merge sort), walks them accumulating wt, and advances the owner rank when the running weight crosses the next
    !! even share of the total. Writes owner(1:n) (rank per item) and cut(0:num_procs-1) (running Morton upper bound per rank; ranks
    !! that receive no item inherit the predecessor's bound so cut is non-decreasing, its top the global max key). Shared by the
    !! fine-block anchor split and the L0-tile split so the two owner maps cannot drift.
    subroutine s_amr_sfc_cut(keys, wt, n, cut, owner)

        integer, intent(in)          :: n
        integer(kind=8), intent(in)  :: keys(n)
        real(wp), intent(in)         :: wt(n)
        integer(kind=8), intent(out) :: cut(0:num_procs - 1)
        integer, intent(out)         :: owner(n)
        integer                      :: ord(n), mrg(n), k, r
        integer                      :: width, lo_m, mid_m, hi_m, i_m, j_m, t_m
        real(wp)                     :: total, cum, tgt, tol

        cut = -1_8
        if (n < 1) return

        ! Sort item indices by Morton key. Bottom-up merge sort: O(n log n), stable (ties keep their original
        ! order), iterative, and a pure function of the input; every rank must produce byte-identical
        ! order or the assignment diverges and s_amr_validate_owner aborts. O(n log n) matters because the
        ! design targets boxes_per_level >> num_procs.
        do k = 1, n
            ord(k) = k
        end do
        width = 1
        do while (width < n)
            lo_m = 1
            do while (lo_m <= n - width)
                mid_m = lo_m + width - 1
                hi_m = min(lo_m + 2*width - 1, n)
                i_m = lo_m; j_m = mid_m + 1; t_m = lo_m
                do while (i_m <= mid_m .and. j_m <= hi_m)
                    ! <= keeps the left run first on ties: stability
                    if (keys(ord(i_m)) <= keys(ord(j_m))) then
                        mrg(t_m) = ord(i_m); i_m = i_m + 1
                    else
                        mrg(t_m) = ord(j_m); j_m = j_m + 1
                    end if
                    t_m = t_m + 1
                end do
                do while (i_m <= mid_m); mrg(t_m) = ord(i_m); i_m = i_m + 1; t_m = t_m + 1; end do
                do while (j_m <= hi_m); mrg(t_m) = ord(j_m); j_m = j_m + 1; t_m = t_m + 1; end do
                ord(lo_m:hi_m) = mrg(lo_m:hi_m)
                lo_m = lo_m + 2*width
            end do
            width = 2*width
        end do

        total = 0._wp
        do k = 1, n
            total = total + wt(k)
        end do

        ! chains-on-chains over the items in SFC order; advance the owner rank when the cumulative weight crosses the next even
        ! share.
        ! All-real arithmetic on replicated weights in a fixed order, so every rank computes the identical assignment.
        r = 0; cum = 0._wp
        do k = 1, n
            tgt = real(r + 1, wp)*total/real(num_procs, wp)
            ! cum is an n-term accumulation while tgt is closed form over another n-term sum, so at an exact share boundary the
            ! two differ by rounding rather than by intent, and the comparison turns on 1 ULP. With integer-valued cost terms
            ! (footprint cells, K_ib, K_pc x integer counts) the arithmetic is exact, but a fractional cost term can split
            ! identical weights unevenly (e.g. 5/3 instead of 4/4). Tolerance is the accumulated rounding bound, O(n) ULP of
            ! the target; far from a boundary it is negligible and the greedy is unchanged.
            tol = spacing(tgt)*real(n, wp)
            if (cum >= tgt - tol .and. r < num_procs - 1) r = r + 1
            owner(ord(k)) = r
            cut(r) = keys(ord(k))  ! items visited in ascending Morton key => running upper bound for rank r
            cum = cum + wt(ord(k))
        end do
        ! ranks that received no item (or trail the last-assigned rank) inherit the predecessor's bound so the search never lands on
        ! them: cut is non-decreasing and its top equals the global max key.
        do r = 1, num_procs - 1
            if (cut(r) < cut(r - 1)) cut(r) = cut(r - 1)
        end do

    end subroutine s_amr_sfc_cut

    !> Fine-level distribution map: assigns each active block a single owner rank by chains-on-chains balancing of fine-work weight
    !! in Morton order of the block's low corner (the same SFC idea m_sfc_partition uses for the base grid, at block granularity).
    !! Fine work = measured coarse-footprint cost (s_amr_block_cost) x the level's refinement factor per active dim, so blocks
    !! concentrating IB or phase-change work weigh more than equal-size quiescent ones. The cost vector is allreduced (one
    !! collective; every rank must call this), after which the assignment is deterministic and identical on every rank.
    !! s_set_amr_fine_geometry applies it as amr_rank_owns_block = (amr_block_owner(amr_cur) == proc_rank).
    !> Refresh the owned-block list if any owner write has happened since the last refresh. Rebuilt lazily rather than in
    !! s_amr_assign_block_owners because that routine is not the only writer of amr_block_owner.
    impure subroutine s_amr_refresh_my_blocks()

        integer :: b

        if (.not. amr_myblk_dirty) return
        if (allocated(amr_my_blk)) then
            if (size(amr_my_blk) < amr_max_blocks) deallocate (amr_my_blk)
        end if
        if (.not. allocated(amr_my_blk)) allocate (amr_my_blk(amr_max_blocks))
        amr_n_my = 0
        do b = 1, amr_num_blocks
            if (amr_block_owner(b) /= proc_rank) cycle
            amr_n_my = amr_n_my + 1
            amr_my_blk(amr_n_my) = b
        end do
        amr_myblk_dirty = .false.

    end subroutine s_amr_refresh_my_blocks

    !> Rebuild the epoch-keyed block lists when the mesh epoch has moved: one O(global blocks) walk per regrid rather than one per
    !! RK stage. Also caches the parent index (f_amr_parent_block is itself an O(global blocks) scan, so calling it per block per
    !! stage would be quadratic in the global block count) and the children adjacency used by s_amr_sibling_face_weights. The cache
    !! build is that quadratic walk, paid once per regrid.
    impure subroutine s_amr_refresh_lists()

        integer :: b, rlo(3), rhi(3), milo(3), mihi(3), bl(3), bh(3)
        integer :: crlo(3), crhi(3), plo(3), phi(3), pl(3), ph(3), mar, pblk, acc, nxt

        if (amr_l1r_epoch == amr_mesh_epoch .and. amr_l1r_nblk == amr_num_blocks) return
        #:for A in ['amr_l1r_blk', 'amr_l1p_blk', 'amr_fch_blk', 'amr_own_blk', 'amr_parent_blk', 'amr_child_idx']
            if (allocated(${A}$)) then
                if (size(${A}$) < amr_max_blocks) deallocate (${A}$)
            end if
            if (.not. allocated(${A}$)) allocate (${A}$ (amr_max_blocks))
        #:endfor
        if (allocated(amr_child_ptr)) then
            if (size(amr_child_ptr) < amr_max_blocks + 1) deallocate (amr_child_ptr)
        end if
        if (.not. allocated(amr_child_ptr)) allocate (amr_child_ptr(0:amr_max_blocks))
        call s_amr_rank_interior(proc_rank, milo, mihi)
        call s_amr_rank_coarse_range(proc_rank, crlo, crhi)
        mar = amr_cpat_mar
        amr_n_l1p = 0
        amr_n_l1r = 0
        amr_n_fch = 0
        amr_n_own = 0
        amr_child_ptr(0:amr_num_blocks) = 0
        do b = 1, amr_num_blocks
            if (amr_owns_all(b)) then
                amr_n_own = amr_n_own + 1
                amr_own_blk(amr_n_own) = b
            end if
            amr_parent_blk(b) = 0
            if (amr_block_level(b) >= 2) then
                pblk = f_amr_parent_block(b)
                amr_parent_blk(b) = pblk
                ! pblk == 0 (no parent found) is structurally impossible on a nested mesh, but an unguarded count would land in
                ! amr_child_ptr(0) (the base offset parent 1's query reads), turning a broken mesh into silent corruption
                if (pblk > 0) then
                    amr_child_ptr(pblk) = amr_child_ptr(pblk) + 1
                    if (amr_block_owner(pblk) == proc_rank .and. amr_block_owner(b) /= proc_rank) then
                        amr_n_fch = amr_n_fch + 1
                        amr_fch_blk(amr_n_fch) = b
                    end if
                end if
                cycle
            end if
            if (amr_block_level(b) /= 1) cycle
            if (amr_block_owner(b) == proc_rank) cycle
            rlo = 0; rhi = 0
            rlo(1) = amr_region_lo_all(1, b); rhi(1) = amr_region_hi_all(1, b)
            if (n_glb > 0) then; rlo(2) = amr_region_lo_all(2, b); rhi(2) = amr_region_hi_all(2, b); end if
            if (p_glb > 0) then; rlo(3) = amr_region_lo_all(3, b); rhi(3) = amr_region_hi_all(3, b); end if
            call s_amr_box_isect(rlo, rhi, milo, mihi, bl, bh)
            if (.not. (bl(1) > bh(1) .or. bl(2) > bh(2) .or. bl(3) > bh(3))) then
                amr_n_l1r = amr_n_l1r + 1
                amr_l1r_blk(amr_n_l1r) = b
            end if
            plo = rlo; phi = rhi
            plo(1) = plo(1) - mar; phi(1) = phi(1) + mar
            if (n_glb > 0) then; plo(2) = plo(2) - mar; phi(2) = phi(2) + mar; end if
            if (p_glb > 0) then; plo(3) = plo(3) - mar; phi(3) = phi(3) + mar; end if
            call s_amr_box_isect(plo, phi, crlo, crhi, pl, ph)
            if (.not. (pl(1) > ph(1) .or. pl(2) > ph(2) .or. pl(3) > ph(3))) then
                amr_n_l1p = amr_n_l1p + 1
                amr_l1p_blk(amr_n_l1p) = b
            end if
        end do
        ! counts -> exclusive prefix, then the ascending fill bumps each amr_child_ptr(p) back up to p's inclusive end, restoring
        ! the query invariant: children of p = amr_child_idx(amr_child_ptr(p-1)+1 : amr_child_ptr(p))
        acc = 0
        do b = 1, amr_num_blocks
            nxt = amr_child_ptr(b); amr_child_ptr(b) = acc; acc = acc + nxt
        end do
        do b = 1, amr_num_blocks
            if (amr_block_level(b) < 2) cycle
            pblk = amr_parent_blk(b)
            if (pblk == 0) cycle  ! uncounted above; a bump here would clobber a real entry
            amr_child_ptr(pblk) = amr_child_ptr(pblk) + 1
            amr_child_idx(amr_child_ptr(pblk)) = b
        end do
        amr_l1r_epoch = amr_mesh_epoch
        amr_l1r_nblk = amr_num_blocks

    end subroutine s_amr_refresh_lists

    !> Open a keyed-tag wave: clear the in-wave per-peer seq counters (touched entries only) and bump the band's generation. Must be
    !! called at wave entry on every rank (the wave call sites are rank-unconditional).
    impure subroutine s_amr_m1_wave_open(band)

        integer, intent(in) :: band
        integer             :: i, d

        if (.not. allocated(amr_tsq)) then
            allocate (amr_tsq(0:num_procs - 1,2), amr_tsq_tch(num_procs, 2))
            amr_tsq = 0; amr_n_tsq = 0
        end if
        do d = 1, 2
            do i = 1, amr_n_tsq(d)
                amr_tsq(amr_tsq_tch(i, d), d) = 0
            end do
            amr_n_tsq(d) = 0
        end do
        amr_tag_gen(band) = mod(amr_tag_gen(band) + 1, 16)

    end subroutine s_amr_m1_wave_open

    !> Next seq on the (peer, dir) channel of the open wave. dir: 1 = send, 2 = recv; the two directions of the same peer are
    !! different channels and never share a counter.
    impure integer function f_amr_m1_seq(peer, idir) result(sq)

        integer, intent(in) :: peer, idir

        if (amr_tsq(peer, idir) == 0) then
            amr_n_tsq(idir) = amr_n_tsq(idir) + 1
            amr_tsq_tch(amr_n_tsq(idir), idir) = peer
        end if
        amr_tsq(peer, idir) = amr_tsq(peer, idir) + 1
        sq = amr_tsq(peer, idir)
        @:ASSERT(sq < 4096, "keyed-tag seq overflows its 12-bit field")

    end function f_amr_m1_seq

    pure integer function f_amr_m1_tag(band, sq) result(t)

        integer, intent(in) :: band, sq

        t = amr_m1_base + band*65536 + amr_tag_gen(band)*4096 + sq

    end function f_amr_m1_tag

    impure subroutine s_amr_assign_block_owners()

        integer :: k, a, lev, maxlev, na
        ! heap, not stack: these are O(global boxes), and at large box counts the seven together would overflow a default stack
        integer, allocatable         :: aidx(:), aown(:)
        integer(kind=8), allocatable :: key(:), akey(:)
        real(wp), allocatable        :: wt(:), cost(:), awt(:)

        if (amr_num_blocks < 1) return

        allocate (aidx(amr_num_blocks), aown(amr_num_blocks), key(amr_num_blocks), akey(amr_num_blocks), wt(amr_num_blocks), &
                  & cost(amr_num_blocks), awt(amr_num_blocks))

        call s_amr_block_cost(cost)

        ! per-block own fine-work weight = footprint cost x amr_ref_ratio**(level*active dims). A level-l block is amr_ref_ratio**l
        ! finer than L0 per dim, so its work is the footprint cost x rr**(l*d). The level factor only scales blocks within a level
        ! relative to each other (each level is cut separately), but it stays because a level's boxes can differ in footprint.
        ! With no cost signals this reduces to the fine cell count (geometry only).
        do k = 1, amr_num_blocks
            wt(k) = cost(k)*real(amr_ref_ratio, wp)**amr_block_level(k)
            if (n_glb > 0) wt(k) = wt(k)*real(amr_ref_ratio, wp)**amr_block_level(k)
            if (p_glb > 0) wt(k) = wt(k)*real(amr_ref_ratio, wp)**amr_block_level(k)
            key(k) = f_morton(amr_region_lo_all(1, k), amr_region_lo_all(2, k), amr_region_lo_all(3, k))
        end do

        ! Per-level distribution: balance every level independently, each block on its own weight. A level-1 block and its
        ! descendants are assigned separately, so a deep tower does not pin its whole subtree (weight cost*rr**(l*d)) to one
        ! rank. The parent<->child gather/restrict/reflux paths are P2P, so a split tower costs messages rather than correctness.
        !
        ! One cut per level, not one mixed cut over all fine blocks: same-level boxes are disjoint and so have distinct Morton
        ! keys, which the cut-point binary search in f_amr_owner needs. Mixed, a level-2 block sharing its parent's region_lo would
        ! collide with it and the search could not tell them apart.
        !
        ! Each level's cut goes into amr_fine_cut(:, lev): fine blocks straddle tiles, so their owner is not tile-cut-derivable and
        ! f_amr_owner reads amr_fine_cut for them. amr_owner_cut mirrors level 1 only without tiles, where the two are the same
        ! authority. Under coexist amr_owner_cut holds the tile cut that s_l0_tiles_init built; overwriting it here is harmless at
        ! init (the assigner runs first) but at regrid time would clobber the tile cut.
        ! Fine blocks occupy slots (l0_slot_off, amr_num_blocks]; slots [1, l0_slot_off] are the L0 tile prefix. At init the
        ! assigner runs before s_l0_tiles_init, so those prefix slots are still uninitialized (level reads 1 and region_lo is all
        ! zeros, i.e. Morton key 0) and must be excluded: a key-0 block can only ever resolve to rank 0 (cut is non-decreasing
        ! and the search returns the first r with key <= cut(r)), so a phantom key-0 block placed on a higher rank makes
        ! s_amr_validate_owner abort.
        maxlev = maxval(amr_block_level(l0_slot_off + 1:amr_num_blocks))
        do lev = 1, maxlev
            na = 0
            do k = l0_slot_off + 1, amr_num_blocks
                if (amr_block_level(k) /= lev) cycle
                na = na + 1
                akey(na) = key(k); awt(na) = wt(k); aidx(na) = k
            end do
            if (na < 1) cycle
            call s_amr_sfc_cut(akey, awt, na, amr_fine_cut(:,lev), aown)
            do a = 1, na
                amr_block_owner(aidx(a)) = aown(a); amr_myblk_dirty = .true.
            end do
            if (lev == 1 .and. l0_slot_off == 0) amr_owner_cut = amr_fine_cut(:,1)
        end do

        call s_amr_validate_owner()
        call s_amr_report_balance(wt, maxlev)

        deallocate (aidx, aown, key, akey, wt, cost, awt)

    end subroutine s_amr_assign_block_owners

    !> Per-level and total load-balance report: max/mean assigned block weight over ranks, the metric the balancer minimises.
    !! Without it a distribution change can only be judged by end-to-end s/step, which cannot separate "balanced" from "uniformly
    !! slow".
    !!
    !! Needs no MPI: wt, amr_block_level and amr_block_owner are replicated and identical on every rank (the cost vector is
    !! allreduced in s_amr_block_cost), so every rank computes the same numbers and rank 0 prints. ratio == 1 is perfect balance;
    !! ratio == num_procs means one rank holds everything at that level. no_blocks_ranks counts ranks holding no block at this
    !! level, which is the granularity floor showing up directly: a level with fewer boxes than ranks cannot balance, however good
    !! the cut is. It is not an idleness measure: those ranks still own level-0 work (level 0 covers every rank) and may own
    !! blocks at other levels. Only m_rank_timing measures idleness.
    impure subroutine s_amr_report_balance(wt, maxlev)

        real(wp), intent(in)  :: wt(:)
        integer, intent(in)   :: maxlev
        real(wp), allocatable :: rw(:), tw(:), rc(:)
        real(wp)              :: mx, mean, cmx, cmean
        integer               :: k, lev, nb, empty

        if (.not. load_weight_wrt) return
        if (proc_rank /= 0) return

        ! heap, not automatic: these are num_procs long, and automatic arrays would put O(num_procs) on the stack (the module
        ! puts amr_owner_cut / amr_fine_cut on the heap for the same reason)
        allocate (rw(0:num_procs - 1), tw(0:num_procs - 1), rc(0:num_procs - 1))
        tw = 0._wp
        do lev = 1, maxlev
            rw = 0._wp; rc = 0._wp
            nb = 0
            do k = 1, amr_num_blocks
                if (amr_block_level(k) /= lev) cycle
                rw(amr_block_owner(k)) = rw(amr_block_owner(k)) + wt(k)
                rc(amr_block_owner(k)) = rc(amr_block_owner(k)) + 1._wp
                nb = nb + 1
            end do
            if (nb == 0) cycle
            tw = tw + rw
            mx = maxval(rw); mean = sum(rw)/real(num_procs, wp)
            empty = count(rw <= 0._wp)
            ! Box-count imbalance beside weight imbalance. cost(k) is a footprint cell count, but per-block advance cost has a
            ! large fixed component regardless of block size, so a rank's true load also tracks how many boxes it holds. Equal
            ! cells with unequal box counts would read as perfectly balanced and run skewed; printing both shows it.
            cmx = maxval(rc); cmean = sum(rc)/real(num_procs, wp)
            ! Not merge(): merge is a function, so both arms are evaluated and the mean == 0 arm would still divide by zero.
            ! no_blocks_ranks counts ranks holding no block at this level; they are not idle (they still own level-0 work and
            ! possibly other levels), they just take no share of this level's.
            ! cmean cannot be zero here: nb >= 1, so sum(rc) = nb >= 1. No guard needed.
            if (mean > 0._wp) print '(A,I0,A,I0,A,I0,A,F8.3,A,F8.3,A,I0,A,I0)', ' [amr-balance] level ', lev, ': boxes ', nb, &
                & '/ranks ', num_procs, ' max/mean ', mx/mean, ' boxes_max/mean ', cmx/cmean, ' no_blocks_ranks ', empty, ' of ', &
                & num_procs
        end do
        ! per-rank weight (fine cells), so rhs time can be regressed against actual load
        if (proc_rank == 0) then
            write (*, '(A)', advance='no') ' [amr-balance] per-rank fine_work :'
            do k = 0, num_procs - 1; write (*, '(I12)', advance='no') nint(tw(k), kind=8); end do
            write (*, '(A)') ''
        end if
        mean = sum(tw)/real(num_procs, wp)
        ! fine_work = sum of the assigned weights = the fine cells advanced per step (with no cost signals wt is exactly that).
        ! Without it an AMR-vs-uniform wall-clock ratio cannot separate "more cells advanced" from "more overhead per cell".
        if (mean > 0._wp) print '(A,F8.3,A,I0,A,I0)', ' [amr-balance] TOTAL   : max/mean ', maxval(tw)/mean, &
            & ' ranks_with_no_fine_block ', count(tw <= 0._wp), ' fine_work ', nint(sum(tw), kind=8)
        deallocate (rw, tw, rc)

    end subroutine s_amr_report_balance

    !> Owner rank of block k from the O(num_procs) SFC cut-points: binary-search k's own Morton key in the owning authority's cut. A
    !! level-0 tile resolves against amr_owner_cut (the tile cut in tiled modes); a fine block (level>=1) resolves against its own
    !! level's cut amr_fine_cut(:, level) (level 1's == amr_owner_cut in no-tile AMR). Under per-level distribution a block's owner
    !! depends on its own key and level alone. Reproduces s_amr_assign_block_owners' / the tile split's cost-weighted SFC assignment
    !! exactly.
    pure integer function f_amr_owner(k) result(r)

        integer, intent(in) :: k
        integer(kind=8)     :: mk, cut(0:num_procs - 1)
        integer             :: a, lo, hi, mid

        a = k
        if (amr_block_level(k) == 0) then
            cut = amr_owner_cut  ! tile: own Morton key vs the tile cut
        else
            cut = amr_fine_cut(:,amr_block_level(k))
        end if
        mk = f_morton(amr_region_lo_all(1, a), amr_region_lo_all(2, a), amr_region_lo_all(3, a))
        lo = 0; hi = num_procs - 1
        do while (lo < hi)
            mid = (lo + hi)/2
            if (mk <= cut(mid)) then
                hi = mid
            else
                lo = mid + 1
            end if
        end do
        r = lo

    end function f_amr_owner

    !> Guard: the SFC cut-point accessor must reproduce the stored owner table exactly.
    impure subroutine s_amr_validate_owner()

        integer :: k

        do k = 1, amr_num_blocks
            ! Skip the L0 tile prefix when its cut has not been built yet (amr_owner_cut still -1): those
            ! slots are uninitialized at assigner time and carry a stale level with Morton key 0. The
            ! tile-init call site populates both cuts and validates them there.
            if (k <= l0_slot_off .and. amr_owner_cut(num_procs - 1) < 0_8) cycle
            ! every block resolves: tiles (level 0) via amr_owner_cut (tile cut), fine blocks (level>=1) via amr_fine_cut. The
            ! caller
            ! guarantees the relevant cut is populated for the blocks present at each call site (assigner: fine cut; tile init:
            ! both).
            if (f_amr_owner(k) /= amr_block_owner(k)) &
                & call s_mpi_abort('SFC cut-point owner disagrees with amr_block_owner - cut capture or search is wrong')
        end do

    end subroutine s_amr_validate_owner

    impure subroutine s_amr_compute_isect(lo, hi)

        integer, intent(in) :: lo(3), hi(3)
        integer             :: sidx(3), ext(3), d

        sidx = 0; ext = 0
        sidx(1) = start_idx(1); ext(1) = m
        if (n_glb > 0) then; sidx(2) = start_idx(2); ext(2) = n; end if
        if (p_glb > 0) then; sidx(3) = start_idx(3); ext(3) = p; end if
        do d = 1, 3
            amr_isect_lo(d) = max(lo(d), sidx(d))
            amr_isect_hi(d) = min(hi(d), sidx(d) + ext(d))
        end do
        amr_rank_owns_block = amr_isect_lo(1) <= amr_isect_hi(1)
        if (n_glb > 0) amr_rank_owns_block = amr_rank_owns_block .and. amr_isect_lo(2) <= amr_isect_hi(2)
        if (p_glb > 0) amr_rank_owns_block = amr_rank_owns_block .and. amr_isect_lo(3) <= amr_isect_hi(3)

    end subroutine s_amr_compute_isect

    !> Rank r's coarse interior box (global) from the computed decomposition (s_amr_rank_decomp, no ghosts). Covered coarse cells
    !! are in-domain, so restriction targets are identified by interior overlap alone.
    pure subroutine s_amr_rank_interior(r, ilo, ihi)

        integer, intent(in)  :: r
        integer, intent(out) :: ilo(3), ihi(3)
        integer              :: sidx(3), ext(3)

        call s_amr_rank_decomp(r, sidx, ext)
        ilo = 0; ihi = 0
        ilo(1) = sidx(1); ihi(1) = sidx(1) + ext(1)
        if (n_glb > 0) then; ilo(2) = sidx(2); ihi(2) = sidx(2) + ext(2); end if
        if (p_glb > 0) then; ilo(3) = sidx(3); ihi(3) = sidx(3) + ext(3); end if

    end subroutine s_amr_rank_interior

    !> Expand a candidate regrid box (global indices) to fully contain every immersed body it overlaps, with a buff_size margin (the
    !! IB image-point stencils need resolved surroundings). Expansion is re-clamped to the domain interior by the caller's own
    !! guards; a body too large for the per-rank block cap aborts with a named message. The bbox reads the live centroid, so a
    !! moving body's box tracks its current position; between regrids s_amr_update_mib_fine guards containment.

    !> Margin-padded global coarse-index bounding box of immersed body i (supported analytic geometries only; aborts on others).
    !! Reads the body's live centroid, so a moving body's box tracks its current position.
    impure subroutine s_amr_body_bbox(i, mrg, blo, bhi)

        integer, intent(in)  :: i, mrg
        integer, intent(out) :: blo(3), bhi(3)
        real(wp)             :: c(3), half(3)

        c = [patch_ib(i)%x_centroid, patch_ib(i)%y_centroid, patch_ib(i)%z_centroid]
        select case (patch_ib(i)%geometry)
        case (2, 8, 10)  ! circle, sphere, cylinder: radius-bounded (cylinder length adds below)
            half = patch_ib(i)%radius
            if (patch_ib(i)%geometry == 10) then
                half(1) = max(half(1), 0.5_wp*patch_ib(i)%length_x)
                half(2) = max(half(2), 0.5_wp*patch_ib(i)%length_y)
                half(3) = max(half(3), 0.5_wp*patch_ib(i)%length_z)
            end if
        case (3, 9)  ! rectangle, box
            half = 0.5_wp*[patch_ib(i)%length_x, patch_ib(i)%length_y, patch_ib(i)%length_z]
        case default
            call s_mpi_abort('amr dynamic regrid with ib: unsupported body geometry for the ' &
                             & // 'containment bounding box (supported: circle/rectangle/sphere/box/cylinder)')
        end select
        ! physical bbox -> global coarse indices: uniform spacing only (stretched grids with ib-dynamic-regrid/Lagrangian are
        ! aborted
        ! at init; the axisymmetric half axis cell only shrinks dy(0), so the floor is still conservative)
        blo(1) = int((c(1) - half(1) - glb_bounds(1)%beg)/dx(0)) - mrg
        bhi(1) = int((c(1) + half(1) - glb_bounds(1)%beg)/dx(0)) + mrg
        blo(2) = 0; bhi(2) = 0; blo(3) = 0; bhi(3) = 0
        if (n_glb > 0) then
            blo(2) = int((c(2) - half(2) - glb_bounds(2)%beg)/dy(min(1, n))) - mrg
            bhi(2) = int((c(2) + half(2) - glb_bounds(2)%beg)/dy(min(1, n))) + mrg
        end if
        if (p_glb > 0) then
            blo(3) = int((c(3) - half(3) - glb_bounds(3)%beg)/dz(0)) - mrg
            bhi(3) = int((c(3) + half(3) - glb_bounds(3)%beg)/dz(0)) + mrg
        end if

    end subroutine s_amr_body_bbox

    impure subroutine s_amr_expand_box_over_bodies(lo, hi)

        integer, intent(inout) :: lo(3), hi(3)
        integer                :: i, d, blo(3), bhi(3), mrg
        logical                :: ovl

        ! containment margin: the IB image-point stencil reaches a few cells beyond the surface (the static-block goldens keep
        ! ~5); buff_size (floored to 10 by ib) would exceed the per-rank block cap for ordinary bodies. For amr_max_level > 1
        ! the
        ! body must survive every child nesting inset (amr_cpat_mar per level down to amr_max_level), so the parent block clears the
        ! body by that many extra cells, keeping the finest C/F boundary a full image-point stencil off the surface (refining the
        ! surface, not the interior).

        mrg = max(amr_buf, 4) + max(0, amr_max_level - 1)*amr_cpat_mar

        do i = 1, num_ibs
            call s_amr_body_bbox(i, mrg, blo, bhi)
            ! blocks must stay buff_size inside the domain: a body whose margin-padded bbox does not fit cannot be contained; fail
            ! with a named message instead of a clipped body
            if (blo(1) < buff_size .or. bhi(1) > m_glb - buff_size .or. (n_glb > 0 .and. (blo(2) < buff_size .or. bhi(2) > n_glb &
                & - buff_size)) .or. (p_glb > 0 .and. (blo(3) < buff_size .or. bhi(3) > p_glb - buff_size))) then
                call s_mpi_abort('amr dynamic regrid with ib: the immersed body plus its containment ' &
                                 & // 'margin does not fit inside the refinable domain interior (blocks stay buff_size off the edges)')
            end if
            ovl = lo(1) <= bhi(1) .and. hi(1) >= blo(1)
            if (n_glb > 0) ovl = ovl .and. lo(2) <= bhi(2) .and. hi(2) >= blo(2)
            if (p_glb > 0) ovl = ovl .and. lo(3) <= bhi(3) .and. hi(3) >= blo(3)
            if (.not. ovl) cycle
            do d = 1, num_dims
                lo(d) = min(lo(d), blo(d))
                hi(d) = max(hi(d), bhi(d))
            end do
            if (hi(1) - lo(1) + 1 > amr_maxc_fit(1) .or. (n_glb > 0 .and. hi(2) - lo(2) + 1 > amr_maxc_fit(2)) .or. (p_glb > 0 &
                & .and. hi(3) - lo(3) + 1 > amr_maxc_fit(3))) then
                call s_mpi_abort('amr dynamic regrid with ib: containing the immersed body plus margin ' &
                                 & // 'exceeds the per-rank block size cap; use fewer ranks or a larger amr_maxc_fit')
            end if
        end do

    end subroutine s_amr_expand_box_over_bodies

    !> max_grid_size tiling: split box [lo:hi] into a grid of contiguous sub-boxes each <= amr_maxc_fit per dim (the max a rank can
    !! whole-own), appending them to out(nt+1:). Tiles are adjacent (share fine seams); the block-to-block fine-fine halo makes
    !! those seams conservative. Even split: ntl = ceil(ext/amr_maxc_fit) tiles, each of size ceil(ext/ntl) <= amr_maxc_fit. Sets
    !! capped=1 and stops if the amr_max_blocks cap is hit. Collapsed dims stay [0:0].
    pure subroutine s_amr_tile_box(lo, hi, out, nt, cap, capped, tsz)

        integer, intent(in)        :: lo(3), hi(3), cap
        type(t_box), intent(inout) :: out(:)
        integer, intent(inout)     :: nt, capped
        !> per-dim tile size (default amr_maxc_fit; a level-lev caller passes amr_maxc_fit/amr_ref_ratio**(lev-1): the slot holds
        !! amr_ref_ratio*amr_maxc_fit fine cells and a level-lev block spans amr_ref_ratio**lev per coarse cell)
        integer, intent(in), optional :: tsz(3)
        integer                       :: ntl(3), s(3), t1, t2, t3, qlo(3), qhi(3), tc(3)

        tc = amr_maxc_fit; if (present(tsz)) tc = tsz
        tc = max(tc, 1)  ! a level>=2 caller passes amr_maxc_fit/2, which is 0 when a rank's fine half-extent is 1 (small subdomain
        !                  at high np); a 0 tile size would divide-by-zero below, and a 1-cell tile is the valid floor
        ntl = 1; s = 1
        ntl(1) = (hi(1) - lo(1) + tc(1))/tc(1); s(1) = (hi(1) - lo(1) + ntl(1))/ntl(1)
        if (n_glb > 0) then
            ntl(2) = (hi(2) - lo(2) + tc(2))/tc(2); s(2) = (hi(2) - lo(2) + ntl(2))/ntl(2)
        end if
        if (p_glb > 0) then
            ntl(3) = (hi(3) - lo(3) + tc(3))/tc(3); s(3) = (hi(3) - lo(3) + ntl(3))/ntl(3)
        end if
        do t3 = 0, ntl(3) - 1
            qlo(3) = 0; qhi(3) = 0
            if (p_glb > 0) then; qlo(3) = lo(3) + t3*s(3); qhi(3) = min(lo(3) + (t3 + 1)*s(3) - 1, hi(3)); end if
            do t2 = 0, ntl(2) - 1
                qlo(2) = 0; qhi(2) = 0
                if (n_glb > 0) then; qlo(2) = lo(2) + t2*s(2); qhi(2) = min(lo(2) + (t2 + 1)*s(2) - 1, hi(2)); end if
                do t1 = 0, ntl(1) - 1
                    if (nt >= cap) then; capped = 1; return; end if
                    qlo(1) = lo(1) + t1*s(1); qhi(1) = min(lo(1) + (t1 + 1)*s(1) - 1, hi(1))
                    nt = nt + 1; out(nt)%lo = qlo; out(nt)%hi = qhi
                end do
            end do
        end do

    end subroutine s_amr_tile_box

end module m_amr_distribution
