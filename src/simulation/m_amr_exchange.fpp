!>
!!@file
!!@brief Contains module m_amr_exchange

#! AMD OpenMP lane: assert allocatables present on every kernel here (see OMP_DEFAULT_STR). Every conditionally allocated
#! module array a kernel here names launches only under its allocation's own condition (amr_rvw: cyl_coord; sw_jac/jac: igr;
#! amr_cg_pb/mv: do_pbmv; amr_gst_a/b: amr_subcycle; amr_prim_st/amr_bt_*: amr_prim_batch); amr_cg and amr_cons_br/stor_st are
#! allocated before first use. A kernel naming an unallocated array aborts. Keep it so.
#:set MFC_OMP_PRESENT_ALLOCATABLE = True
#:include 'macros.fpp'

!> @brief Parent/child and fine-fine data exchange: gather plans, pack/unpack, seams, ghost fills and the fill waves.
module m_amr_exchange

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
    use m_amr_distribution

    implicit none

    private
    public :: f_amr_seam, f_amr_seam_dim, s_amr_build_gather_plan, s_amr_build_seam_pairs, s_amr_exchange_coarse_cons_halo, &
        & s_amr_fill_fine_ghosts_pbmv, s_amr_fine_fine_drain, s_amr_fine_fine_halo, s_amr_fine_fine_post, s_amr_fw_szi, &
        & s_amr_fw_szi3, s_amr_fw_szr, s_amr_gather_chunk_post, s_amr_gather_chunk_send, s_amr_gather_coarse_patch, &
        & s_amr_gather_coarse_patch_pbmv, s_amr_gather_consume_box, s_amr_gather_from_parent_field_cons, s_amr_gather_send_flush, &
        & s_amr_parent_fill_wave, s_amr_recv_parent_patch, s_amr_stage_fill_wave, s_l0_pack_unpack_block_sf, &
        & s_l0_pack_unpack_block_st

contains

    !> Fine-level distribution: assemble the current block's coarse patch on its owner. The patch covers global coarse cells
    !! region_lo-amr_cpat_mar : region_hi+amr_cpat_mar (the full reach of every prolongation/ghost-fill stencil) for all sys_size
    !! variables, stored in amr_cg in a block-local frame (cell 0 == global amr_cpat_off). Point-to-point: the owner receives the
    !! patch cells it does not hold from exactly the coarse owners that hold them; each rank's contribution is the patch intersected
    !! with its contiguous owned coarse range (s_amr_rank_coarse_range, computed from the cartesian decomposition). Non-participants
    !! send/recv nothing (no global collective). At np=1 the owner just copies its own coarse over the patch, bit-for-bit. Runtime
    !! (pull_host) packs/unpacks the overlap boxes on the device (q_coarse device-current with valid ghosts); init/regrid fills from
    !! the host (host-current with valid ghosts). Packed data is wp, cast to stp into amr_cg (identity for stp coarse),
    !! device-current on exit. Invariant: "coarse" here means the block's parent level (level l-1), not the base grid (level 0). For
    !! a level-1 block the parent is L0, but a level>=2 block folds to/from its parent block's fine array; the C<->F
    !! prolong/restrict/gather routines all operate in the parent-fine frame, not the L0 frame. Twin s_amr_gather_coarse_patch_pbmv
    !! (q<->pb/mv): same P2P skeleton (rank-range, intersection, pack/send/recv/unpack) and patch-local frame; keep them in
    !! lockstep.
    !> Make room for one more pending gather send, draining the pool first if it is full. Draining is a WAITALL, so the pool size
    !! sets how far a contributing rank may run ahead of the owners.
    impure subroutine s_amr_gsnd_reserve(slotsz)

        integer, intent(in) :: slotsz

        if (.not. allocated(amr_gsnd_pool)) then
            allocate (amr_gsnd_pool(slotsz, amr_gsnd_max), amr_gsnd_req(amr_gsnd_max))
            amr_gsnd_n = 0
        else if (size(amr_gsnd_pool, 1) < slotsz) then
            call s_amr_gather_send_flush()  ! outstanding sends reference the old buffer - complete them before resizing
            deallocate (amr_gsnd_pool)
            allocate (amr_gsnd_pool(slotsz, amr_gsnd_max))
        end if
        if (amr_gsnd_n >= amr_gsnd_max) call s_amr_gather_send_flush()

    end subroutine s_amr_gsnd_reserve

    !> Complete every pending gather send. Must be called before the send buffers are reused or the routine returns to a caller that
    !! will free them: an ISEND whose buffer is overwritten in flight silently corrupts the receiver's patch.
    impure subroutine s_amr_gather_send_flush()

        integer :: ierr

        if (amr_gsnd_n == 0) return
#ifdef MFC_MPI
        call s_wait_tic()
        call MPI_WAITALL(amr_gsnd_n, amr_gsnd_req(1:amr_gsnd_n), MPI_STATUSES_IGNORE, ierr)
        call s_wait_toc(WT_REGRID)
#endif
        amr_gsnd_n = 0

    end subroutine s_amr_gather_send_flush

    !> Build amr_korder/amr_kpos for nboxes regrid boxes (see the declaration): per level in ascending order, per-owner FIFOs of the
    !! level's boxes (ascending box id inside each), emitted round-robin over owners.
    impure subroutine s_amr_build_korder(nboxes)

        integer, intent(in)  :: nboxes
        integer              :: k, lev, r, p, maxlev
        integer, allocatable :: cnt(:), head(:), tail(:), nxt(:)

        if (allocated(amr_korder)) then
            if (size(amr_korder) < nboxes) deallocate (amr_korder, amr_kpos)
        end if
        if (.not. allocated(amr_korder)) allocate (amr_korder(max(nboxes, 1)), amr_kpos(max(nboxes, 1)))
        if (.not. amr_korder_rot) then
            do k = 1, nboxes
                amr_korder(k) = k; amr_kpos(k) = k
            end do
            return
        end if
        allocate (cnt(0:num_procs - 1), head(0:num_procs - 1), tail(0:num_procs - 1), nxt(max(nboxes, 1)))
        maxlev = 0
        do k = 1, nboxes
            maxlev = max(maxlev, amr_block_level(f_l0_slot(k)))
        end do
        p = 0
        do lev = 1, maxlev
            cnt = 0; head = 0; tail = 0
            do k = 1, nboxes
                if (amr_block_level(f_l0_slot(k)) /= lev) cycle
                r = amr_block_owner(f_l0_slot(k))
                if (cnt(r) == 0) then
                    head(r) = k
                else
                    nxt(tail(r)) = k
                end if
                tail(r) = k; nxt(k) = 0; cnt(r) = cnt(r) + 1
            end do
            do while (any(cnt > 0))
                do r = 0, num_procs - 1
                    if (cnt(r) == 0) cycle
                    k = head(r); head(r) = nxt(k); cnt(r) = cnt(r) - 1
                    p = p + 1; amr_korder(p) = k; amr_kpos(k) = p
                end do
            end do
        end do
        @:ASSERT(p == nboxes, "rebuild walk order: box count mismatch")
        deallocate (cnt, head, tail, nxt)

    end subroutine s_amr_build_korder

    !> Derive the entire rebuild gather message set up front (per level-1 box its contributor ranks and message sizes, per level>=2
    !! box its parent source and size) from the same replicated caches the per-box path reads (amr_region_*_all, amr_ovl_gather,
    !! amr_block_owner, rank coarse ranges, s_amr_parent_foot). Caller (s_amr_regrid_rebuild_slots) clears amr_gpl_valid when its
    !! box loop ends.
    impure subroutine s_amr_build_gather_plan()

        integer :: i, ks, idx, r, nsrc, mo, pblk, im, ifc, ip, km, kf, kp
        integer :: v1hi, v2hi, v3hi, plo(3), phi(3), crlo(3), crhi(3), bl(3), bh(3), w(3)

        ! the per-box path lazily rebuilds the overlap lists inside the first gather; force the same rebuild here so the plan
        ! and the boxes read identical lists

        if (amr_seam_pairs_dirty .or. amr_seam_pairs_nblk /= amr_num_blocks) call s_amr_build_seam_pairs()
        ! the migrate step installed the new regions/levels/owners and bumped the epoch, so the lists rebuild on the new mesh here
        ! (amr_own_blk reads amr_owns_all, still the previous generation until the geometry pass; nothing in the rebuild reads
        ! it, and the reconcile's epoch bump rebuilds it before the first stage does)
        call s_amr_refresh_my_blocks()
        call s_amr_refresh_lists()
        if (allocated(amr_gpk)) then
            if (size(amr_gpk) < amr_max_blocks) deallocate (amr_gpk)
        end if
        if (.not. allocated(amr_gpk)) allocate (amr_gpk(amr_max_blocks))
        ! three-cursor ascending merge; the L0 tile prefix (slots <= l0_slot_off, owned like any block) carries no regrid box
        amr_n_gpk = 0
        im = 1; ifc = 1; ip = 1
        do
            km = huge(1); kf = huge(1); kp = huge(1)
            if (im <= amr_n_my) km = amr_my_blk(im)
            if (ifc <= amr_n_fch) kf = amr_fch_blk(ifc)
            if (ip <= amr_n_l1p) kp = amr_l1p_blk(ip)
            ks = min(km, kf, kp)
            if (ks == huge(1)) exit
            if (km == ks) im = im + 1
            if (kf == ks) ifc = ifc + 1
            if (kp == ks) ip = ip + 1
            if (ks <= l0_slot_off) cycle
            amr_n_gpk = amr_n_gpk + 1
            amr_gpk(amr_n_gpk) = ks
        end do
        ! re-emit the participants in the rebuild walk order (amr_korder): the chunk loop reads amr_gpk as one run per chunk
        call s_amr_build_korder(amr_num_blocks - l0_slot_off)
        block
            logical, allocatable :: part(:)
            integer              :: pp, np_gpk
            allocate (part(amr_num_blocks)); part = .false.
            do pp = 1, amr_n_gpk
                part(amr_gpk(pp)) = .true.
            end do
            np_gpk = 0
            do pp = 1, amr_num_blocks - l0_slot_off
                ks = f_l0_slot(amr_korder(pp))
                if (part(ks)) then
                    np_gpk = np_gpk + 1; amr_gpk(np_gpk) = ks
                end if
            end do
            @:ASSERT(np_gpk == amr_n_gpk, "gather plan: walk order lost a participant")
            deallocate (part)
        end block
        mo = size(amr_ovl_gather, 1)
        if (allocated(amr_gpl_src)) then
            if (size(amr_gpl_src, 1) < mo) deallocate (amr_gpl_src, amr_gpl_sz)
        end if
        if (.not. allocated(amr_gpl_nsrc)) allocate (amr_gpl_nsrc(amr_max_blocks), amr_gpl_psrc(amr_max_blocks), &
            & amr_gpl_psz(amr_max_blocks))
        if (.not. allocated(amr_gpl_src)) allocate (amr_gpl_src(mo, amr_max_blocks), amr_gpl_sz(mo, amr_max_blocks))
        ! plan entries for the participants only: they are the only boxes the chunk post/send/consume below ever look up
        do i = 1, amr_n_gpk
            ks = amr_gpk(i)
            amr_gpl_nsrc(ks) = 0; amr_gpl_psrc(ks) = -1; amr_gpl_psz(ks) = 0
            if (amr_block_level(ks) >= 2) then
                pblk = amr_parent_blk(ks)
                if (amr_block_owner(pblk) /= amr_block_owner(ks)) then
                    call s_amr_parent_foot(ks, pblk, plo, phi)
                    w = 0
                    w(1) = (phi(1) - plo(1)) + 2*amr_cpat_mar
                    if (n_glb > 0) w(2) = (phi(2) - plo(2)) + 2*amr_cpat_mar
                    if (p_glb > 0) w(3) = (phi(3) - plo(3)) + 2*amr_cpat_mar
                    amr_gpl_psrc(ks) = amr_block_owner(pblk)
                    amr_gpl_psz(ks) = sys_size*(w(1) + 1)*(w(2) + 1)*(w(3) + 1)
                end if
            else
                ! level-1 patch box: same arithmetic as the gather's patch-frame block (collapsed dims stay 0)
                plo = 0
                plo(1) = amr_region_lo_all(1, ks) - amr_cpat_mar
                if (n_glb > 0) plo(2) = amr_region_lo_all(2, ks) - amr_cpat_mar
                if (p_glb > 0) plo(3) = amr_region_lo_all(3, ks) - amr_cpat_mar
                v1hi = (amr_region_hi_all(1, ks) - amr_region_lo_all(1, ks)) + 2*amr_cpat_mar
                v2hi = 0; v3hi = 0
                if (n_glb > 0) v2hi = (amr_region_hi_all(2, ks) - amr_region_lo_all(2, ks)) + 2*amr_cpat_mar
                if (p_glb > 0) v3hi = (amr_region_hi_all(3, ks) - amr_region_lo_all(3, ks)) + 2*amr_cpat_mar
                phi(1) = plo(1) + v1hi; phi(2) = plo(2) + v2hi; phi(3) = plo(3) + v3hi
                nsrc = 0
                do idx = 1, amr_ovl_gather_n(ks)
                    r = amr_ovl_gather(idx, ks)
                    if (r == amr_block_owner(ks)) cycle
                    call s_amr_rank_coarse_range(r, crlo, crhi)
                    call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
                    nsrc = nsrc + 1
                    amr_gpl_src(nsrc, ks) = r
                    amr_gpl_sz(nsrc, ks) = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                end do
                amr_gpl_nsrc(ks) = nsrc
            end if
        end do
        amr_gpl_valid = .true.

    end subroutine s_amr_build_gather_plan

    !> Chunked gather, post phase: pre-post every recv this rank needs for boxes [c_lo, c_hi] (level-1 contributor slices and split
    !! level>=2 parent patches) straight from the plan into the flat chunk pool, tag = slot, appended in box order so box k's
    !! requests are one contiguous run. Ownership from amr_block_owner only (amr_owns_all / amr_rank_owns_block still mirror the
    !! previous generation until the consume phase's geometry call). Contains no MPI waits; reallocating the pool here is safe
    !! because every recv posted for the previous chunk was completed inside that chunk's consume phase. The chunk's boxes this rank
    !! has a role in are amr_gpk(i0:i1); c_lo is the chunk's first box (chunk-local indexing of the request runs).
    impure subroutine s_amr_gather_chunk_post(c_lo, i0, i1)

        integer, intent(in) :: c_lo, i0, i1
        integer             :: i, ks, cb, idx, need, nreq, off, ierr

        @:ASSERT(amr_gpl_valid, "chunk gather: no plan")
        call s_phase_tic(PH_RBPOST)
        need = 0; nreq = 0
        do i = i0, i1
            ks = amr_gpk(i)
            if (amr_block_owner(ks) /= proc_rank) cycle
            ! + XA_NH per message: the exchange-audit identity header rides ahead of each payload (zero in production)
            if (amr_block_level(ks) >= 2) then
                if (amr_gpl_psrc(ks) >= 0) then
                    need = need + amr_gpl_psz(ks) + XA_NH; nreq = nreq + 1
                end if
            else
                do idx = 1, amr_gpl_nsrc(ks)
                    need = need + amr_gpl_sz(idx, ks) + XA_NH
                end do
                nreq = nreq + amr_gpl_nsrc(ks)
            end if
        end do
        if (allocated(amr_gcr_pool)) then
            if (size(amr_gcr_pool) < need) deallocate (amr_gcr_pool)
        end if
        if (need > 0 .and. .not. allocated(amr_gcr_pool)) allocate (amr_gcr_pool(need))
        if (allocated(amr_gcr_req)) then
            if (size(amr_gcr_req) < nreq) deallocate (amr_gcr_req, amr_gcr_off)
        end if
        if (nreq > 0 .and. .not. allocated(amr_gcr_req)) allocate (amr_gcr_req(nreq), amr_gcr_off(nreq))

        amr_gcr_n = 0; off = 0
        amr_gcr_r0(:) = 1; amr_gcr_nr(:) = 0; amr_gcr_sent(:) = .false.
        do i = i0, i1
            ks = amr_gpk(i)
            if (amr_block_owner(ks) /= proc_rank) cycle
            cb = amr_kpos(ks - l0_slot_off) - c_lo + 1
            amr_gcr_r0(cb) = amr_gcr_n + 1
#ifdef MFC_MPI
            if (amr_block_level(ks) >= 2) then
                if (amr_gpl_psrc(ks) >= 0) then
                    amr_gcr_n = amr_gcr_n + 1
                    amr_gcr_off(amr_gcr_n) = off
                    call s_xa_rec(XA_F2_RCV, 2, amr_gpl_psz(ks), ks)
                    call MPI_IRECV(amr_gcr_pool(off + 1), amr_gpl_psz(ks) + XA_NH, mpi_p, amr_gpl_psrc(ks), ks, MPI_COMM_WORLD, &
                                   & amr_gcr_req(amr_gcr_n), ierr)
                    off = off + amr_gpl_psz(ks) + XA_NH
                    amr_gcr_nr(cb) = 1
                end if
            else
                do idx = 1, amr_gpl_nsrc(ks)
                    amr_gcr_n = amr_gcr_n + 1
                    amr_gcr_off(amr_gcr_n) = off
                    call s_xa_rec(XA_F1_RCV, 2, amr_gpl_sz(idx, ks), ks)
                    call MPI_IRECV(amr_gcr_pool(off + 1), amr_gpl_sz(idx, ks) + XA_NH, mpi_p, amr_gpl_src(idx, ks), ks, &
                                   & MPI_COMM_WORLD, amr_gcr_req(amr_gcr_n), ierr)
                    off = off + amr_gpl_sz(idx, ks) + XA_NH
                end do
                amr_gcr_nr(cb) = amr_gpl_nsrc(ks)
            end if
#endif
        end do
        call s_phase_toc(PH_RBPOST)

    end subroutine s_amr_gather_chunk_post

    !> Chunked gather, send phase: issue this rank's sends for boxes [c_lo, c_hi]. Level-1: pack the host slice of q_coarse (host is
    !! truth during rebuild) and ISEND through the deferred pool, driven by the plan. Level>=2 split pairs: send only when the
    !! parent was consumed in an earlier chunk (pblk < f_l0_slot(c_lo), monotone slot map); a same-chunk parent's new-generation
    !! store is not built until its own consume iteration, so that send stays at the child's consume position, where parents-first
    !! ordering guarantees the parent is complete. Geometry from the replicated caches only: no s_set_amr_fine_geometry swap, no
    !! amr_cur. Walks the chunk's participants amr_gpk(i0:i1) with the per-box predicates intact: the list only drops boxes that
    !! would have cycled (a sender is a parent-owner or a level-1 contributor).
    impure subroutine s_amr_gather_chunk_send(q_coarse, c_lo, i0, i1)

        type(scalar_field), dimension(sys_size), intent(in) :: q_coarse
        integer, intent(in)                                 :: c_lo, i0, i1
        integer                                             :: ks, cb, idx, i, ii, g1, g2, g3, o1, o2, o3, boxsz, maxsz, pblk, ierr
        integer                                             :: v1hi, v2hi, v3hi, plo(3), phi(3), crlo(3), crhi(3), bl(3), bh(3)
        logical                                             :: contrib

        o1 = start_idx(1); o2 = 0; o3 = 0
        if (n_glb > 0) o2 = start_idx(2)
        if (p_glb > 0) o3 = start_idx(3)
        do ii = i0, i1
            ks = amr_gpk(ii)
            cb = amr_kpos(ks - l0_slot_off) - c_lo + 1
            if (amr_block_level(ks) >= 2) then
                if (amr_gpl_psrc(ks) < 0) cycle  ! co-located: no message
                pblk = amr_parent_blk(ks)
                if (amr_block_owner(pblk) /= proc_rank) cycle  ! not the sender
                if (amr_kpos(pblk - l0_slot_off) >= c_lo) cycle  ! same-chunk parent: send at the child's consume position
                call s_phase_tic(PH_PGSEND)
                call s_amr_gather_from_parent_field_cons(ks, pblk, amr_loc_of(pblk), .true.)
                call s_phase_toc(PH_PGSEND)
                amr_gcr_sent(cb) = .true.
            else
                if (amr_block_owner(ks) == proc_rank) cycle  ! the owner receives
                contrib = .false.
                do idx = 1, amr_gpl_nsrc(ks)
                    if (amr_gpl_src(idx, ks) == proc_rank) contrib = .true.
                end do
                if (.not. contrib) cycle
                ! same patch-frame arithmetic as the plan builder and the per-box gather
                plo = 0
                plo(1) = amr_region_lo_all(1, ks) - amr_cpat_mar
                if (n_glb > 0) plo(2) = amr_region_lo_all(2, ks) - amr_cpat_mar
                if (p_glb > 0) plo(3) = amr_region_lo_all(3, ks) - amr_cpat_mar
                v1hi = (amr_region_hi_all(1, ks) - amr_region_lo_all(1, ks)) + 2*amr_cpat_mar
                v2hi = 0; v3hi = 0
                if (n_glb > 0) v2hi = (amr_region_hi_all(2, ks) - amr_region_lo_all(2, ks)) + 2*amr_cpat_mar
                if (p_glb > 0) v3hi = (amr_region_hi_all(3, ks) - amr_region_lo_all(3, ks)) + 2*amr_cpat_mar
                phi(1) = plo(1) + v1hi; phi(2) = plo(2) + v2hi; phi(3) = plo(3) + v3hi
                call s_amr_rank_coarse_range(proc_rank, crlo, crhi)
                call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
                boxsz = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                maxsz = sys_size*(v1hi + 1)*(v2hi + 1)*(v3hi + 1)
                call s_phase_tic(PH_RBRSV)
                call s_amr_gsnd_reserve(maxsz + XA_NH)
                call s_phase_toc(PH_RBRSV)
                amr_gsnd_n = amr_gsnd_n + 1
                call s_phase_tic(PH_RBPACK)
                if (XA_NH > 0) call s_xa_hdr_pack(amr_gsnd_pool(:,amr_gsnd_n), XA_F1_SND, ks, bl, bh)
                idx = XA_NH
                do i = 1, sys_size
                    do g3 = bl(3), bh(3)
                        do g2 = bl(2), bh(2)
                            do g1 = bl(1), bh(1)
                                idx = idx + 1
                                amr_gsnd_pool(idx, amr_gsnd_n) = real(q_coarse(i)%sf(g1 - o1, g2 - o2, g3 - o3), wp)
                            end do
                        end do
                    end do
                end do
                call s_phase_toc(PH_RBPACK)
#ifdef MFC_MPI
                call s_phase_tic(PH_RBSEND)
                call s_xa_rec(XA_F1_SND, 1, boxsz, ks)
                call MPI_ISEND(amr_gsnd_pool(1, amr_gsnd_n), boxsz + XA_NH, mpi_p, amr_block_owner(ks), ks, MPI_COMM_WORLD, &
                               & amr_gsnd_req(amr_gsnd_n), ierr)
                call s_phase_toc(PH_RBSEND)
#endif
            end if
        end do

    end subroutine s_amr_gather_chunk_send

    !> Chunked gather, consume phase: the per-box gather body with the exchange already in flight. Fills amr_cg for the current box
    !! (amr_cur, geometry already set by the caller) from the own slice plus the chunk pool's pre-posted recvs. Level-1 owner:
    !! own-box host copy, one WAITALL on this box's contiguous request run, host unpack per contributor (plan order = posting
    !! order), device push. Level>=2: co-located parent = local device copy; split parent = the parent owner packs and sends here
    !! when the parent shares this chunk (amr_gcr_sent marks the ones the send phase already covered), the child owner waits and
    !! device-unpacks its single pre-posted recv. Called for the boxes this rank owns or parents (the caller's owner-cycle comes
    !! after); every owned box's requests are waited unconditionally inside its own chunk, which is what makes the request arrays
    !! reusable next chunk.
    impure subroutine s_amr_gather_consume_box(q_coarse, k, c_lo)

        type(scalar_field), dimension(sys_size), intent(in) :: q_coarse
        integer, intent(in) :: k, c_lo
        integer :: cb, idx, i, r, g1, g2, g3, o1, o2, o3, boxsz, pblk, w1, w2, w3, ierr, r0, nr, off
        integer :: v1hi, v2hi, v3hi, plo(3), phi(3), crlo(3), crhi(3), bl(3), bh(3)

        cb = amr_kpos(k) - c_lo + 1
        r0 = amr_gcr_r0(cb); nr = amr_gcr_nr(cb)

        if (amr_block_level(amr_cur) >= 2) then
            call s_phase_tic(PH_PGALL)
            pblk = amr_parent_blk(amr_cur)
            ! the deferred same-chunk send below reads the parent's store, valid only because parents-first ordering already
            ! consumed the parent; trip immediately if the ordering is ever violated
            @:ASSERT(pblk < f_l0_slot(k), "chunk gather: parent box not before child")
            if (amr_gpl_psrc(amr_cur) < 0) then
                ! co-located: the owner's local device copy (the field routine detects co-location itself)
                if (amr_block_owner(amr_cur) == proc_rank) then
                    call s_phase_tic(PH_PGSEND)
                    call s_amr_gather_from_parent_field_cons(amr_cur, pblk, amr_loc_of(pblk), .true.)
                    call s_phase_toc(PH_PGSEND)
                end if
            else if (amr_block_owner(pblk) == proc_rank) then
                ! split, parent side: a same-chunk parent could not be packed in the send phase (its store was unbuilt);
                ! parents-first ordering means it is complete now
                if (.not. amr_gcr_sent(cb)) then
                    call s_phase_tic(PH_PGSEND)
                    call s_amr_gather_from_parent_field_cons(amr_cur, pblk, amr_loc_of(pblk), .true.)
                    call s_phase_toc(PH_PGSEND)
                end if
            else if (amr_block_owner(amr_cur) == proc_rank) then
                ! split, child side: wait on the pre-posted parent patch and unpack on the device
                call s_amr_parent_foot(amr_cur, pblk, plo, phi)
                amr_cpat_off = 0
                amr_cpat_off(1) = plo(1) - amr_cpat_mar
                if (n_glb > 0) amr_cpat_off(2) = plo(2) - amr_cpat_mar
                if (p_glb > 0) amr_cpat_off(3) = plo(3) - amr_cpat_mar
                w1 = (phi(1) - plo(1)) + 2*amr_cpat_mar
                w2 = 0; w3 = 0
                if (n_glb > 0) w2 = (phi(2) - plo(2)) + 2*amr_cpat_mar
                if (p_glb > 0) w3 = (phi(3) - plo(3)) + 2*amr_cpat_mar
#ifdef MFC_MPI
                call s_phase_tic(PH_PGRECV)
                call s_wait_tic()
                call MPI_WAITALL(nr, amr_gcr_req(r0:r0 + nr - 1), MPI_STATUSES_IGNORE, ierr)
                call s_wait_toc(WT_REGRID)
                call s_phase_toc(PH_PGRECV)
                off = amr_gcr_off(r0)
                boxsz = amr_gpl_psz(amr_cur)
                if (XA_NH > 0) call s_xa_hdr_check(amr_gcr_pool(off + 1:off + XA_NH), XA_F2_SND, amr_cur, plo, phi)
                call s_amr_unpack_parent_patch_device(w1, w2, w3, amr_gcr_pool(off + XA_NH + 1:off + XA_NH + boxsz), .true.)
#endif
            end if
            call s_phase_toc(PH_PGALL)
            return
        end if

        ! level-1: same patch frame and own fill as the per-box gather; the recvs are already posted
        amr_cpat_off = 0
        amr_cpat_off(1) = amr_region_lo_all(1, amr_cur) - amr_cpat_mar
        if (n_glb > 0) amr_cpat_off(2) = amr_region_lo_all(2, amr_cur) - amr_cpat_mar
        if (p_glb > 0) amr_cpat_off(3) = amr_region_lo_all(3, amr_cur) - amr_cpat_mar
        v1hi = (amr_region_hi_all(1, amr_cur) - amr_region_lo_all(1, amr_cur)) + 2*amr_cpat_mar
        v2hi = 0; v3hi = 0
        if (n_glb > 0) v2hi = (amr_region_hi_all(2, amr_cur) - amr_region_lo_all(2, amr_cur)) + 2*amr_cpat_mar
        if (p_glb > 0) v3hi = (amr_region_hi_all(3, amr_cur) - amr_region_lo_all(3, amr_cur)) + 2*amr_cpat_mar
        plo = amr_cpat_off
        phi(1) = amr_cpat_off(1) + v1hi; phi(2) = amr_cpat_off(2) + v2hi; phi(3) = amr_cpat_off(3) + v3hi

        if (amr_block_owner(amr_cur) /= proc_rank) return  ! contributor sends were the send phase's job

        o1 = start_idx(1); o2 = 0; o3 = 0
        if (n_glb > 0) o2 = start_idx(2)
        if (p_glb > 0) o3 = start_idx(3)
        call s_amr_rank_coarse_range(proc_rank, crlo, crhi)
        call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
        call s_phase_tic(PH_RBOWN)
        call s_amr_unpack_patch(q_coarse, bl, bh, o1, o2, o3)
        call s_phase_toc(PH_RBOWN)
#ifdef MFC_MPI
        if (nr > 0) then
            call s_phase_tic(PH_RBWAIT)
            call s_wait_tic()
            call MPI_WAITALL(nr, amr_gcr_req(r0:r0 + nr - 1), MPI_STATUSES_IGNORE, ierr)
            call s_wait_toc(WT_REGRID)
            call s_phase_toc(PH_RBWAIT)
            call s_phase_tic(PH_RBUNPK)
            do idx = 1, nr
                ! plan order = posting order; recompute each contributor's slice box exactly as the plan builder did
                call s_amr_rank_coarse_range(amr_gpl_src(idx, amr_cur), crlo, crhi)
                call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
                off = amr_gcr_off(r0 + idx - 1)
                if (XA_NH > 0) call s_xa_hdr_check(amr_gcr_pool(off + 1:off + XA_NH), XA_F1_SND, amr_cur, bl, bh)
                r = XA_NH
                do i = 1, sys_size
                    do g3 = bl(3), bh(3)
                        do g2 = bl(2), bh(2)
                            do g1 = bl(1), bh(1)
                                r = r + 1
                                amr_cg(i)%sf(g1 - amr_cpat_off(1), g2 - amr_cpat_off(2), &
                                       & g3 - amr_cpat_off(3)) = real(amr_gcr_pool(off + r), stp)
                            end do
                        end do
                    end do
                end do
            end do
            call s_phase_toc(PH_RBUNPK)
        end if
#endif
        call s_phase_tic(PH_RBUPD)
        do i = 1, sys_size
            $:GPU_UPDATE(device='[amr_cg(i)%sf]')
        end do
        call s_phase_toc(PH_RBUPD)

    end subroutine s_amr_gather_consume_box

    impure subroutine s_amr_gather_coarse_patch(q_coarse, pull_host)

        type(scalar_field), dimension(sys_size), intent(in) :: q_coarse
        !> runtime callers pass .true. (coarse device-current); init/regrid pass .false. (host is truth)
        logical, intent(in)   :: pull_host
        integer               :: i, g1, g2, g3, o1, o2, o3, owner, r, idx, boxsz, maxsz, nsrc, ierr
        integer               :: v1hi, v2hi, v3hi, plo(3), phi(3), crlo(3), crhi(3), bl(3), bh(3)
        real(wp), allocatable :: rbuf(:,:)
        integer, allocatable  :: reqs(:), srank(:)

        ! multi-level: a level>=2 block's coarse side is its parent block's fine cells, not the L0 base grid q_coarse; gather
        ! amr_cg
        ! from the parent's fine array in the parent-fine frame (isect already parent-fine from s_set_amr_fine_geometry).

        if (amr_block_level(amr_cur) >= 2) then
            call s_amr_gather_from_parent(pull_host)
            return
        end if

        ! block-local patch frame (cell 0 == global region_lo-nmar; collapsed dims -> 0) + its global cell range [plo:phi]
        amr_cpat_off = 0
        amr_cpat_off(1) = amr_region_lo_all(1, amr_cur) - amr_cpat_mar
        if (n_glb > 0) amr_cpat_off(2) = amr_region_lo_all(2, amr_cur) - amr_cpat_mar
        if (p_glb > 0) amr_cpat_off(3) = amr_region_lo_all(3, amr_cur) - amr_cpat_mar
        v1hi = (amr_region_hi_all(1, amr_cur) - amr_region_lo_all(1, amr_cur)) + 2*amr_cpat_mar
        v2hi = 0; v3hi = 0
        if (n_glb > 0) v2hi = (amr_region_hi_all(2, amr_cur) - amr_region_lo_all(2, amr_cur)) + 2*amr_cpat_mar
        if (p_glb > 0) v3hi = (amr_region_hi_all(3, amr_cur) - amr_region_lo_all(3, amr_cur)) + 2*amr_cpat_mar
        plo = amr_cpat_off
        phi(1) = amr_cpat_off(1) + v1hi; phi(2) = amr_cpat_off(2) + v2hi; phi(3) = amr_cpat_off(3) + v3hi

        owner = amr_block_owner(amr_cur)
        o1 = start_idx(1); o2 = 0; o3 = 0
        if (n_glb > 0) o2 = start_idx(2)
        if (p_glb > 0) o3 = start_idx(3)
        maxsz = sys_size*(v1hi + 1)*(v2hi + 1)*(v3hi + 1)

        ! np=1: the sole owner holds every covered coarse cell, so copy q_coarse->amr_cg on-device (same index map as
        ! s_amr_unpack_patch), skipping the device->host->device round-trip. Only for pull_host; init/regrid (.not. pull_host) falls
        ! through to the host path (device copy may be stale).
        if (num_procs == 1 .and. pull_host) then
            call s_amr_rank_coarse_range(owner, crlo, crhi)
            call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
            call s_amr_gather_own_box_device(q_coarse, bl, bh, o1, o2, o3)  ! same kernel the np>1 owner path uses
            return
        end if

        ! np>1 runtime (pull_host): no full-field host pull. The owner's own-box copy, the non-owner pack, and the received-box
        ! unpacks all run on the device over only the overlap boxes, so just the contiguous wire buffers cross PCIe (MPI stays on
        ! host buffers). Init/regrid (.not. pull_host): host is truth, so the host pack/unpack paths below read it directly.

        ! block set changed: rebuild the cached overlap-rank lists (same lazy trigger as s_amr_fine_fine_halo; local, replicated)
        if (amr_seam_pairs_dirty .or. amr_seam_pairs_nblk /= amr_num_blocks) call s_amr_build_seam_pairs()

        if (proc_rank == owner) then
            ! fill the cells this rank holds locally (own box), then receive the rest from the other coarse-owners
            call s_amr_rank_coarse_range(proc_rank, crlo, crhi)
            call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
            if (pull_host) then
                ! runtime: q_coarse is device-current - copy the own box on the device (same index map/assignment as the host path)
                call s_amr_gather_own_box_device(q_coarse, bl, bh, o1, o2, o3)
            else
                call s_amr_unpack_patch(q_coarse, bl, bh, o1, o2, o3)  ! local read: q_coarse own frame -> amr_cg patch frame
            end if
            ! count + post recvs from every other rank whose owned range overlaps the patch (cached list; every listed rank
            ! overlaps by construction)
            nsrc = 0
            do idx = 1, amr_ovl_gather_n(amr_cur)
                if (amr_ovl_gather(idx, amr_cur) /= owner) nsrc = nsrc + 1
            end do
            if (nsrc > 0) then
                allocate (rbuf(maxsz + XA_NH, nsrc), reqs(nsrc), srank(nsrc))
                nsrc = 0
                do idx = 1, amr_ovl_gather_n(amr_cur)
                    r = amr_ovl_gather(idx, amr_cur)
                    if (r == owner) cycle
                    call s_amr_rank_coarse_range(r, crlo, crhi)
                    call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
                    boxsz = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                    nsrc = nsrc + 1; srank(nsrc) = r
#ifdef MFC_MPI
                    call s_xa_rec(XA_F1_RCV, 2, boxsz, amr_cur)
                    call MPI_IRECV(rbuf(1, nsrc), boxsz + XA_NH, mpi_p, r, amr_cur, MPI_COMM_WORLD, reqs(nsrc), ierr)
#endif
                end do
#ifdef MFC_MPI
                call s_wait_tic()
                call MPI_WAITALL(nsrc, reqs, MPI_STATUSES_IGNORE, ierr)
                call s_wait_toc(WT_GATHER)
#endif
                do idx = 1, nsrc
                    call s_amr_rank_coarse_range(srank(idx), crlo, crhi)
                    call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
                    if (XA_NH > 0) call s_xa_hdr_check(rbuf(:,idx), XA_F1_SND, amr_cur, bl, bh)
                    if (pull_host) then
                        ! runtime: unpack only this box's wire buffer on the device (same order/cast as the host unpack below)
                        boxsz = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                        call s_amr_unpack_box_device(bl, bh, rbuf(XA_NH + 1:XA_NH + boxsz,idx))
                        cycle
                    end if
                    ! unpack in the same (i, g3, g2, g1) order the sender packed; place at amr_cg patch-local index
                    r = XA_NH
                    do i = 1, sys_size
                        do g3 = bl(3), bh(3)
                            do g2 = bl(2), bh(2)
                                do g1 = bl(1), bh(1)
                                    r = r + 1
                                    amr_cg(i)%sf(g1 - amr_cpat_off(1), g2 - amr_cpat_off(2), g3 - amr_cpat_off(3)) = real(rbuf(r, &
                                           & idx), stp)
                                end do
                            end do
                        end do
                    end do
                end do
                deallocate (rbuf, reqs, srank)
            end if
            ! host path only: the runtime device path wrote amr_cg on the device directly (host amr_cg stays stale, as at np=1 -
            ! runtime consumers read the device copy)
            if (.not. pull_host) then
                do i = 1, sys_size
                    $:GPU_UPDATE(device='[amr_cg(i)%sf]')
                end do
            end if
        else
            ! non-owner: if my owned coarse range overlaps the patch, pack my slice (wp) and send it to the owner
            call s_amr_rank_coarse_range(proc_rank, crlo, crhi)
            call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
            if (bl(1) <= bh(1) .and. bl(2) <= bh(2) .and. bl(3) <= bh(3)) then
                boxsz = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                call s_amr_gsnd_reserve(maxsz + XA_NH)
                amr_gsnd_n = amr_gsnd_n + 1
                if (pull_host) then
                    ! runtime: pack the overlap box on the device straight into the pool slot (only the box crosses PCIe);
                    ! the slice leaves the audit header words ahead of the data (kernel untouched)
                    call s_amr_pack_box_device(q_coarse, bl, bh, o1, o2, o3, amr_gsnd_pool(XA_NH + 1:,amr_gsnd_n))
                else
                    idx = XA_NH
                    do i = 1, sys_size
                        do g3 = bl(3), bh(3)
                            do g2 = bl(2), bh(2)
                                do g1 = bl(1), bh(1)
                                    idx = idx + 1
                                    amr_gsnd_pool(idx, amr_gsnd_n) = real(q_coarse(i)%sf(g1 - o1, g2 - o2, g3 - o3), wp)
                                end do
                            end do
                        end do
                    end do
                end if
#ifdef MFC_MPI
                ! non-blocking: the owner's per-box IRECV/WAITALL orders the data, and this rank does not rendezvous on every
                ! box. Completed by s_amr_gather_send_flush (caller) or the drain in s_amr_gsnd_reserve.
                if (XA_NH > 0) call s_xa_hdr_pack(amr_gsnd_pool(:,amr_gsnd_n), XA_F1_SND, amr_cur, bl, bh)
                call s_xa_rec(XA_F1_SND, 1, boxsz, amr_cur)
                call MPI_ISEND(amr_gsnd_pool(1, amr_gsnd_n), boxsz + XA_NH, mpi_p, owner, amr_cur, MPI_COMM_WORLD, &
                               & amr_gsnd_req(amr_gsnd_n), ierr)
#endif
            end if
        end if

    end subroutine s_amr_gather_coarse_patch

    !> Non-polytropic QBMM analogue of s_amr_gather_coarse_patch: gather the current block's coarse pb/mv patch into amr_cg_pb/mv
    !! (patch frame, cell 0 == amr_cpat_off), P2P from the coarse-cell owners into the block owner. Per-cell payload = 2*nnode*nb
    !! (pb block then mv block). Single-level only (level>=2 QBMM np>=2 is checker-gated); wire is wp, cast to stp on unpack
    !! (identity for stp coarse), so at np=1 the owner copies its own coarse over the patch bit-for-bit. Twin
    !! s_amr_gather_coarse_patch (pb/mv<->q): mirrors the q_cons gather's P2P skeleton and patch-local frame; keep them in lockstep.
    impure subroutine s_amr_gather_coarse_patch_pbmv(pb_coarse, mv_coarse, pull_host)

        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(in) :: pb_coarse, mv_coarse
        !> runtime callers pass .true. (coarse device-current); init/regrid pass .false. (host is truth)
        logical, intent(in)   :: pull_host
        integer               :: q, ib_, g1, g2, g3, o1, o2, o3, owner, r, idx, boxsz, maxsz, nsrc, ierr
        integer               :: v1hi, v2hi, v3hi, plo(3), phi(3), crlo(3), crhi(3), bl(3), bh(3), cellsz
        real(wp), allocatable :: rbuf(:,:), sbuf(:)
        integer, allocatable  :: reqs(:), srank(:)

        ! single-level only: a level>=2 block's coarse side is its parent's fine pb/mv, distributed only at np=1; the checker gate
        ! keeps multi-level QBMM np>=2 fail-closed, so this must never be reached at level>=2.

        if (amr_block_level(amr_cur) >= 2) return

        cellsz = 2*nnode*nb

        ! block-local patch frame (cell 0 == global region_lo-nmar; collapsed dims -> 0) + its global cell range [plo:phi]
        amr_cpat_off = 0
        amr_cpat_off(1) = amr_region_lo_all(1, amr_cur) - amr_cpat_mar
        if (n_glb > 0) amr_cpat_off(2) = amr_region_lo_all(2, amr_cur) - amr_cpat_mar
        if (p_glb > 0) amr_cpat_off(3) = amr_region_lo_all(3, amr_cur) - amr_cpat_mar
        v1hi = (amr_region_hi_all(1, amr_cur) - amr_region_lo_all(1, amr_cur)) + 2*amr_cpat_mar
        v2hi = 0; v3hi = 0
        if (n_glb > 0) v2hi = (amr_region_hi_all(2, amr_cur) - amr_region_lo_all(2, amr_cur)) + 2*amr_cpat_mar
        if (p_glb > 0) v3hi = (amr_region_hi_all(3, amr_cur) - amr_region_lo_all(3, amr_cur)) + 2*amr_cpat_mar
        plo = amr_cpat_off
        phi(1) = amr_cpat_off(1) + v1hi; phi(2) = amr_cpat_off(2) + v2hi; phi(3) = amr_cpat_off(3) + v3hi

        owner = amr_block_owner(amr_cur)
        o1 = start_idx(1); o2 = 0; o3 = 0
        if (n_glb > 0) o2 = start_idx(2)
        if (p_glb > 0) o3 = start_idx(3)
        maxsz = cellsz*(v1hi + 1)*(v2hi + 1)*(v3hi + 1)

        ! np=1: the sole owner holds every covered coarse cell, so copy pb_coarse/mv_coarse->amr_cg_pb/mv on-device over the
        ! in-domain patch. Only for pull_host; init/regrid (.not. pull_host) falls through to the host path (device copy may be
        ! stale).
        if (num_procs == 1 .and. pull_host) then
            call s_amr_rank_coarse_range(owner, crlo, crhi)
            call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
            call s_amr_gather_own_box_pbmv_device(pb_coarse, mv_coarse, bl, bh, o1, o2, o3)  ! same kernel the np>1 owner path uses
            return
        end if

        ! np>1 runtime (pull_host): no full-field host pull. The owner's own-box copy, the non-owner pack, and the received-box
        ! unpacks all run on the device over only the overlap boxes (mirror of s_amr_gather_coarse_patch). Init/regrid
        ! (.not. pull_host): host is truth, so the host pack/unpack paths below read it directly.

        ! block set changed: rebuild the cached overlap-rank lists (same lazy trigger as s_amr_fine_fine_halo; local, replicated)
        if (amr_seam_pairs_dirty .or. amr_seam_pairs_nblk /= amr_num_blocks) call s_amr_build_seam_pairs()

        if (proc_rank == owner) then
            ! fill the cells this rank holds locally (own box), then receive the rest from the other coarse-owners
            call s_amr_rank_coarse_range(proc_rank, crlo, crhi)
            call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
            if (pull_host) then
                ! runtime: pb/mv device-current - copy the own box on the device (same index map/assignment as the host path)
                call s_amr_gather_own_box_pbmv_device(pb_coarse, mv_coarse, bl, bh, o1, o2, o3)
            else
                do ib_ = 1, nb
                    do q = 1, nnode
                        do g3 = bl(3), bh(3)
                            do g2 = bl(2), bh(2)
                                do g1 = bl(1), bh(1)
                                    amr_cg_pb(g1 - amr_cpat_off(1), g2 - amr_cpat_off(2), g3 - amr_cpat_off(3), q, &
                                              & ib_) = pb_coarse(g1 - o1, g2 - o2, g3 - o3, q, ib_)
                                    amr_cg_mv(g1 - amr_cpat_off(1), g2 - amr_cpat_off(2), g3 - amr_cpat_off(3), q, &
                                              & ib_) = mv_coarse(g1 - o1, g2 - o2, g3 - o3, q, ib_)
                                end do
                            end do
                        end do
                    end do
                end do
            end if
            ! count + post recvs from every other rank whose owned range overlaps the patch (cached list; every listed rank
            ! overlaps by construction)
            nsrc = 0
            do idx = 1, amr_ovl_gather_n(amr_cur)
                if (amr_ovl_gather(idx, amr_cur) /= owner) nsrc = nsrc + 1
            end do
            if (nsrc > 0) then
                allocate (rbuf(maxsz + XA_NH, nsrc), reqs(nsrc), srank(nsrc))
                nsrc = 0
                do idx = 1, amr_ovl_gather_n(amr_cur)
                    r = amr_ovl_gather(idx, amr_cur)
                    if (r == owner) cycle
                    call s_amr_rank_coarse_range(r, crlo, crhi)
                    call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
                    boxsz = cellsz*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                    nsrc = nsrc + 1; srank(nsrc) = r
#ifdef MFC_MPI
                    call s_xa_rec(XA_F3_RCV, 2, boxsz, amr_cur)
                    call MPI_IRECV(rbuf(1, nsrc), boxsz + XA_NH, mpi_p, r, amr_cur, MPI_COMM_WORLD, reqs(nsrc), ierr)
#endif
                end do
#ifdef MFC_MPI
                call s_wait_tic()
                call MPI_WAITALL(nsrc, reqs, MPI_STATUSES_IGNORE, ierr)
                call s_wait_toc(WT_GATHER)
#endif
                do idx = 1, nsrc
                    call s_amr_rank_coarse_range(srank(idx), crlo, crhi)
                    call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
                    if (XA_NH > 0) call s_xa_hdr_check(rbuf(:,idx), XA_F3_SND, amr_cur, bl, bh)
                    if (pull_host) then
                        ! runtime: unpack only this box's wire buffer on the device (same order/cast as the host unpack below)
                        boxsz = cellsz*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                        call s_amr_unpack_box_pbmv_device(bl, bh, rbuf(XA_NH + 1:XA_NH + boxsz,idx))
                        cycle
                    end if
                    ! unpack in the same (ib_, q, g3, g2, g1) order the sender packed - pb block then mv block
                    r = XA_NH
                    do ib_ = 1, nb
                        do q = 1, nnode
                            do g3 = bl(3), bh(3)
                                do g2 = bl(2), bh(2)
                                    do g1 = bl(1), bh(1)
                                        r = r + 1
                                        amr_cg_pb(g1 - amr_cpat_off(1), g2 - amr_cpat_off(2), g3 - amr_cpat_off(3), q, &
                                                  & ib_) = real(rbuf(r, idx), stp)
                                    end do
                                end do
                            end do
                        end do
                    end do
                    do ib_ = 1, nb
                        do q = 1, nnode
                            do g3 = bl(3), bh(3)
                                do g2 = bl(2), bh(2)
                                    do g1 = bl(1), bh(1)
                                        r = r + 1
                                        amr_cg_mv(g1 - amr_cpat_off(1), g2 - amr_cpat_off(2), g3 - amr_cpat_off(3), q, &
                                                  & ib_) = real(rbuf(r, idx), stp)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
                deallocate (rbuf, reqs, srank)
            end if
            ! host path only: the runtime device path wrote amr_cg_pb/mv on the device directly (host copies stay stale, as at np=1
            ! -
            ! runtime consumers read the device copy)
            if (.not. pull_host) then
                $:GPU_UPDATE(device='[amr_cg_pb, amr_cg_mv]')
            end if
        else
            ! non-owner: if my owned coarse range overlaps the patch, pack my slice (wp) and send it to the owner
            call s_amr_rank_coarse_range(proc_rank, crlo, crhi)
            call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
            if (bl(1) <= bh(1) .and. bl(2) <= bh(2) .and. bl(3) <= bh(3)) then
                boxsz = cellsz*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                allocate (sbuf(boxsz + XA_NH))
                if (pull_host) then
                    ! runtime: pack the overlap box on the device straight into sbuf (only the box crosses PCIe);
                    ! the slice leaves the audit header words ahead of the data (kernel untouched)
                    call s_amr_pack_box_pbmv_device(pb_coarse, mv_coarse, bl, bh, o1, o2, o3, sbuf(XA_NH + 1:))
                else
                    idx = XA_NH
                    do ib_ = 1, nb
                        do q = 1, nnode
                            do g3 = bl(3), bh(3)
                                do g2 = bl(2), bh(2)
                                    do g1 = bl(1), bh(1)
                                        idx = idx + 1; sbuf(idx) = real(pb_coarse(g1 - o1, g2 - o2, g3 - o3, q, ib_), wp)
                                    end do
                                end do
                            end do
                        end do
                    end do
                    do ib_ = 1, nb
                        do q = 1, nnode
                            do g3 = bl(3), bh(3)
                                do g2 = bl(2), bh(2)
                                    do g1 = bl(1), bh(1)
                                        idx = idx + 1; sbuf(idx) = real(mv_coarse(g1 - o1, g2 - o2, g3 - o3, q, ib_), wp)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end if
#ifdef MFC_MPI
                if (XA_NH > 0) call s_xa_hdr_pack(sbuf, XA_F3_SND, amr_cur, bl, bh)
                call s_xa_rec(XA_F3_SND, 1, boxsz, amr_cur)
                call MPI_SEND(sbuf, boxsz + XA_NH, mpi_p, owner, amr_cur, MPI_COMM_WORLD, ierr)
#endif
                deallocate (sbuf)
            end if
        end if

        ! host-consumer callers (init/regrid prolong) need the gathered patch on the host
        if (.not. pull_host) then
            $:GPU_UPDATE(host='[amr_cg_pb, amr_cg_mv]')
        end if

    end subroutine s_amr_gather_coarse_patch_pbmv

    !> Multi-level gather: fill amr_cg (the current level>=2 block's coarse patch) from its parent block's fine array, in the
    !! parent-fine cell frame (amr_isect_lo/hi already parent-fine from s_set_amr_fine_geometry). A local copy when the block's
    !! owner also owns the parent; otherwise a point-to-point transfer from the parent owner to the block owner.
    impure subroutine s_amr_gather_from_parent(pull_host)

        logical, intent(in) :: pull_host
        integer             :: pblk

        pblk = f_amr_parent_block(amr_cur)
        ! lock-step fill: gather from the parent's current fine state. pull_host stays in the signature for the level-1 path.
        ! Owner-guard at the call site: the parent slot is allocated only on its owner, and passing its store slot on any other
        ! rank would dereference an unallocated slot. So both participants enter (the parent's owner to pack and send, the block's
        ! owner to receive) and every other rank stays out. When the two coincide (np=1, or a co-located tower) this is a local
        ! copy. to_host = .not. pull_host: init/regrid (pull_host=F) feed the host prolong/self-test; runtime (pull_host=T) reads
        ! amr_cg on the device in the C/F ghost-fill, so skip the device->host copy.
        if (amr_block_owner(pblk) == proc_rank) then
            ! parent owner: local device copy when it also owns the block, otherwise pack and send.
            call s_amr_gather_from_parent_field_cons(amr_cur, pblk, amr_loc_of(pblk), .not. pull_host)
        else if (amr_rank_owns_block) then
            ! block owner only: receive. Deliberately does not take the parent field; amr_slots(pblk) is unallocated here.
            call s_amr_recv_parent_patch(pblk, .not. pull_host)
        end if

    end subroutine s_amr_gather_from_parent

    !> Gather amr_cg (the current level>=2 block's coarse patch) from a specific parent snapshot field qp, in the parent-fine cell
    !! frame (amr_isect_lo/hi already parent-fine from s_set_amr_fine_geometry). substep (qp = the parent slot's q_cons_stor (t^n
    !! bracket) then q_cons (t^{n+1} bracket)) to build the child's two ghost-lerp sources. A local copy when the block's owner also
    !! owns the parent; otherwise point-to-point from the parent owner to the block owner. Two sources, one body: the parent's
    !! conserved state (`_cons`, amr_cons_st) and its SSP-RK stage backup (`_stor`, amr_stor_st), both in the flat store keyed by
    !! the parent's slot.
    #:for GSFX, GARR in [('cons', 'amr_cons_st')]
        impure subroutine s_amr_gather_from_parent_field_${GSFX}$(cblk, pblk, qp, to_host)

            !> the child block (explicit, not amr_cur: the chunked send phase calls this before the consume phase's geometry, when
            !! amr_cur points at another box)
            integer, intent(in) :: cblk
            integer, intent(in) :: pblk
            integer, intent(in) :: qp       !< parent's flat-store slot
            logical, intent(in) :: to_host  !< host copy of amr_cg needed (init/regrid), not runtime
            integer             :: w1, w2, w3, powner, cowner, boxsz, ierr
            integer             :: plo(3), phi(3)

            ! Patch box in the parent-fine frame. Both the child owner and the parent owner must agree on it, so derive it from
            ! replicated metadata (amr_region_*_all + the global amr_ref_ratio) rather than from amr_isect_lo/hi, which is the
            ! empty footprint on a non-owner of this block. On the child owner the two agree by construction
            ! (s_set_amr_fine_geometry).

            call s_amr_parent_foot(cblk, pblk, plo, phi)
            amr_cpat_off = 0
            amr_cpat_off(1) = plo(1) - amr_cpat_mar
            if (n_glb > 0) amr_cpat_off(2) = plo(2) - amr_cpat_mar
            if (p_glb > 0) amr_cpat_off(3) = plo(3) - amr_cpat_mar
            w1 = (phi(1) - plo(1)) + 2*amr_cpat_mar
            w2 = 0; w3 = 0
            if (n_glb > 0) w2 = (phi(2) - plo(2)) + 2*amr_cpat_mar
            if (p_glb > 0) w3 = (phi(3) - plo(3)) + 2*amr_cpat_mar

            cowner = amr_block_owner(cblk); powner = amr_block_owner(pblk)
            if (powner == cowner) then
                ! co-located (always true at np=1, and under tower co-location): straight device copy.
                call s_amr_copy_parent_patch_${GSFX}$(qp, w1, w2, w3, to_host)
                return
            end if

#ifdef MFC_MPI
            ! Split ownership, parent side: exactly one destination (the block's owner) and one box, so a single message
            ! suffices, with no overlap map and no collective (non-participants send/recv nothing, as in the L0<->L1 gather).
            ! Non-blocking, via the same deferred pool the level-1 gather uses (see s_amr_gsnd_reserve), so the parent's owner
            ! does not rendezvous with the child's owner once per box. The pool owns the buffer because an ISEND requires it to
            ! stay live until completion; the drain is s_amr_gather_send_flush after the rebuild's box loop.
            boxsz = sys_size*(w1 + 1)*(w2 + 1)*(w3 + 1)
            ! guard on the plan alone: a send packed short of the plan-sized recv completes short and the consume unpacks stale
            ! pool bytes (a silent wrong answer). amr_gpl_valid is false outside the rebuild box loop, so per-step
            ! calls never consult the plan.
            if (amr_gpl_valid) then
                @:ASSERT(amr_gpl_psz(cblk) == boxsz, "gather plan: parent send size mismatch")
            end if
            call s_amr_gsnd_reserve(boxsz + XA_NH)
            amr_gsnd_n = amr_gsnd_n + 1
            ! header written on the host after the device pack lands (copyout); data at XA_NH+1 via the slice
            call s_amr_pack_parent_patch_device_${GSFX}$(qp, w1, w2, w3, amr_gsnd_pool(XA_NH + 1:,amr_gsnd_n))
            if (XA_NH > 0) call s_xa_hdr_pack(amr_gsnd_pool(:,amr_gsnd_n), XA_F2_SND, cblk, plo, phi)
            call s_xa_rec(XA_F2_SND, 1, boxsz, cblk)
            call MPI_ISEND(amr_gsnd_pool(1, amr_gsnd_n), boxsz + XA_NH, mpi_p, cowner, cblk, MPI_COMM_WORLD, &
                           & amr_gsnd_req(amr_gsnd_n), ierr)
#endif

        end subroutine s_amr_gather_from_parent_field_${GSFX}$
    #:endfor

    !> Receive side of the split-ownership parent gather: fill amr_cg from the parent's owner. Takes only pblk: the parent slot is
    !! not allocated on this rank, so the parent field must not appear in the signature. Recomputes the patch box from the same
    !! replicated metadata the sender uses, so the two agree without a handshake.
    impure subroutine s_amr_recv_parent_patch(pblk, to_host)

        integer, intent(in)   :: pblk
        logical, intent(in)   :: to_host
        integer               :: w1, w2, w3, powner, boxsz, ierr, plo(3), phi(3)
        real(wp), allocatable :: xbuf(:)

        call s_amr_parent_foot(amr_cur, pblk, plo, phi)
        amr_cpat_off = 0
        amr_cpat_off(1) = plo(1) - amr_cpat_mar
        if (n_glb > 0) amr_cpat_off(2) = plo(2) - amr_cpat_mar
        if (p_glb > 0) amr_cpat_off(3) = plo(3) - amr_cpat_mar
        w1 = (phi(1) - plo(1)) + 2*amr_cpat_mar
        w2 = 0; w3 = 0
        if (n_glb > 0) w2 = (phi(2) - plo(2)) + 2*amr_cpat_mar
        if (p_glb > 0) w3 = (phi(3) - plo(3)) + 2*amr_cpat_mar

#ifdef MFC_MPI
        powner = amr_block_owner(pblk)
        boxsz = sys_size*(w1 + 1)*(w2 + 1)*(w3 + 1)
        allocate (xbuf(boxsz + XA_NH))
        call s_xa_rec(XA_F2_RCV, 2, boxsz, amr_cur)
        call s_wait_tic()
        call MPI_RECV(xbuf, boxsz + XA_NH, mpi_p, powner, amr_cur, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call s_wait_toc(WT_PGATHER)
        if (XA_NH > 0) call s_xa_hdr_check(xbuf, XA_F2_SND, amr_cur, plo, phi)
        call s_amr_unpack_parent_patch_device(w1, w2, w3, xbuf(XA_NH + 1:XA_NH + boxsz), to_host)
        deallocate (xbuf)
#endif

    end subroutine s_amr_recv_parent_patch

    !> Device pack of the parent's fine patch into a flat buffer. Same index map as s_amr_copy_parent_patch, writing the send buffer
    !! instead of amr_cg, so the two sides of the P2P gather cannot drift apart.
    #:for GSFX, GARR in [('cons', 'amr_cons_st')]
        #:set QP = lambda ix: GARR + '(g1 + o1, g2 + o2, g3 + o3, ' + ix + ', qp)'
        impure subroutine s_amr_pack_parent_patch_device_${GSFX}$(qp, w1, w2, w3, buf)

            integer, intent(in)                 :: qp  !< parent's flat-store slot
            integer, intent(in)                 :: w1, w2, w3
            real(wp), intent(inout), contiguous :: buf(:)
            integer                             :: i, g1, g2, g3, o1, o2, o3, n1, n2, n3

            o1 = amr_cpat_off(1); o2 = amr_cpat_off(2); o3 = amr_cpat_off(3)
            n1 = w1 + 1; n2 = w2 + 1; n3 = w3 + 1
            $:GPU_PARALLEL_LOOP(collapse=4, copyout='[buf]')
            do i = 1, sys_size
                do g3 = 0, w3
                    do g2 = 0, w2
                        do g1 = 0, w1
                            buf(1 + g1 + n1*(g2 + n2*(g3 + n3*(i - 1)))) = real(${QP('i')}$, wp)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

        end subroutine s_amr_pack_parent_patch_device_${GSFX}$
    #:endfor

    !> Device unpack of a received parent patch into amr_cg. Inverse of s_amr_pack_parent_patch_device; to_host mirrors
    !! s_amr_copy_parent_patch (init/regrid host consumers need the host copy, runtime reads amr_cg on the device).
    impure subroutine s_amr_unpack_parent_patch_device(w1, w2, w3, buf, to_host)

        integer, intent(in)              :: w1, w2, w3
        real(wp), intent(in), contiguous :: buf(:)
        logical, intent(in)              :: to_host
        integer                          :: i, g1, g2, g3, n1, n2, n3

        n1 = w1 + 1; n2 = w2 + 1; n3 = w3 + 1
        $:GPU_PARALLEL_LOOP(collapse=4, copyin='[buf]')
        do i = 1, sys_size
            do g3 = 0, w3
                do g2 = 0, w2
                    do g1 = 0, w1
                        amr_cg(i)%sf(g1, g2, g3) = buf(1 + g1 + n1*(g2 + n2*(g3 + n3*(i - 1))))
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        if (to_host) then
            do i = 1, sys_size
                $:GPU_UPDATE(host='[amr_cg(i)%sf]')
            end do
        end if

    end subroutine s_amr_unpack_parent_patch_device

    !> Device kernel for s_amr_gather_from_parent: copy the parent block's fine patch into amr_cg over [amr_cpat_off : + w]. amr_cg
    !! is then synced to host for host consumers (init self-test's restrict-prolong check). Two sources, one body; see
    !! s_amr_gather_from_parent_field_cons.
    #:for GSFX, GARR in [('cons', 'amr_cons_st')]
        #:set QP = lambda ix: GARR + '(g1 + o1, g2 + o2, g3 + o3, ' + ix + ', qp)'
        impure subroutine s_amr_copy_parent_patch_${GSFX}$(qp, w1, w2, w3, to_host)

            integer, intent(in) :: qp  !< parent's flat-store slot
            integer, intent(in) :: w1, w2, w3
            !> .true. only for the init/regrid host consumers (whole-block host prolong + restrict-prolong self-test). The runtime
            !! C/F ghost-fill reads amr_cg on the device (filled by the kernel below), so no device->host copy is needed.
            logical, intent(in) :: to_host
            integer             :: i, g1, g2, g3, o1, o2, o3

            o1 = amr_cpat_off(1); o2 = amr_cpat_off(2); o3 = amr_cpat_off(3)
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do g3 = 0, w3
                    do g2 = 0, w2
                        do g1 = 0, w1
                            amr_cg(i)%sf(g1, g2, g3) = ${QP('i')}$
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            ! amr_cg is now device-current for the runtime C/F ghost-fill. Sync to host only when a host consumer follows.
            if (to_host) then
                do i = 1, sys_size
                    $:GPU_UPDATE(host='[amr_cg(i)%sf]')
                end do
            end if

        end subroutine s_amr_copy_parent_patch_${GSFX}$
    #:endfor

    !> Sub-box variants of the parent-patch pack/unpack/copy for the ring-clipped parent-fill wave (cons only: the wave ships
    !! q_cons; pb/mv runs keep the full-patch contract). Bounds are patch-local cell ranges; the buffer holds the sub-box in the
    !! same (g1 fastest, sys_size outermost) layout as the full-patch kernels, so both wire sides agree by construction.
    impure subroutine s_amr_pack_parent_box_device_cons(qp, bl, bh, buf)

        integer, intent(in)                 :: qp  !< parent's flat-store slot
        integer, intent(in)                 :: bl(3), bh(3)
        real(wp), intent(inout), contiguous :: buf(:)
        integer                             :: i, g1, g2, g3, o1, o2, o3, n1, n2, n3, l1, l2, l3, u1, u2, u3

        o1 = amr_cpat_off(1); o2 = amr_cpat_off(2); o3 = amr_cpat_off(3)
        l1 = bl(1); l2 = bl(2); l3 = bl(3); u1 = bh(1); u2 = bh(2); u3 = bh(3)
        n1 = u1 - l1 + 1; n2 = u2 - l2 + 1; n3 = u3 - l3 + 1
        $:GPU_PARALLEL_LOOP(collapse=4, copyout='[buf]')
        do i = 1, sys_size
            do g3 = l3, u3
                do g2 = l2, u2
                    do g1 = l1, u1
                        buf(1 + (g1 - l1) + n1*((g2 - l2) + n2*((g3 - l3) + n3*(i - 1)))) = real(amr_cons_st(g1 + o1, g2 + o2, &
                            & g3 + o3, i, qp), wp)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_pack_parent_box_device_cons

    impure subroutine s_amr_unpack_parent_box_device(bl, bh, buf)

        integer, intent(in)              :: bl(3), bh(3)
        real(wp), intent(in), contiguous :: buf(:)
        integer                          :: i, g1, g2, g3, n1, n2, n3, l1, l2, l3, u1, u2, u3

        l1 = bl(1); l2 = bl(2); l3 = bl(3); u1 = bh(1); u2 = bh(2); u3 = bh(3)
        n1 = u1 - l1 + 1; n2 = u2 - l2 + 1; n3 = u3 - l3 + 1
        $:GPU_PARALLEL_LOOP(collapse=4, copyin='[buf]')
        do i = 1, sys_size
            do g3 = l3, u3
                do g2 = l2, u2
                    do g1 = l1, u1
                        amr_cg(i)%sf(g1, g2, g3) = buf(1 + (g1 - l1) + n1*((g2 - l2) + n2*((g3 - l3) + n3*(i - 1))))
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_unpack_parent_box_device

    impure subroutine s_amr_copy_parent_box_cons(qp, bl, bh)

        integer, intent(in) :: qp  !< parent's flat-store slot
        integer, intent(in) :: bl(3), bh(3)
        integer             :: i, g1, g2, g3, o1, o2, o3, l1, l2, l3, u1, u2, u3

        o1 = amr_cpat_off(1); o2 = amr_cpat_off(2); o3 = amr_cpat_off(3)
        l1 = bl(1); l2 = bl(2); l3 = bl(3); u1 = bh(1); u2 = bh(2); u3 = bh(3)
        $:GPU_PARALLEL_LOOP(collapse=4)
        do i = 1, sys_size
            do g3 = l3, u3
                do g2 = l2, u2
                    do g1 = l1, u1
                        amr_cg(i)%sf(g1, g2, g3) = amr_cons_st(g1 + o1, g2 + o2, g3 + o3, i, qp)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_copy_parent_box_cons

    !> Copy this rank's own coarse cells (box [bl:bh] global, read from q_coarse at its own start-idx frame o1/o2/o3) into amr_cg in
    !! the block-local patch frame. stp -> stp, exact.
    impure subroutine s_amr_unpack_patch(q_coarse, bl, bh, o1, o2, o3)

        type(scalar_field), dimension(sys_size), intent(in) :: q_coarse
        integer, intent(in)                                 :: bl(3), bh(3), o1, o2, o3
        integer                                             :: i, g1, g2, g3

        do i = 1, sys_size
            do g3 = bl(3), bh(3)
                do g2 = bl(2), bh(2)
                    do g1 = bl(1), bh(1)
                        amr_cg(i)%sf(g1 - amr_cpat_off(1), g2 - amr_cpat_off(2), g3 - amr_cpat_off(3)) = q_coarse(i)%sf(g1 - o1, &
                               & g2 - o2, g3 - o3)
                    end do
                end do
            end do
        end do

    end subroutine s_amr_unpack_patch

    !> Runtime device analogue of s_amr_unpack_patch: copy the owner's own coarse box [bl:bh] global from q_coarse (device) into
    !! amr_cg (device) in the patch-local frame, with no host round-trip. Same index map and direct stp assignment as the host path.
    !! Twin s_amr_gather_own_box_pbmv_device (q<->pb/mv): same own-box index map; keep them in lockstep.
    impure subroutine s_amr_gather_own_box_device(q_coarse, bl, bh, o1, o2, o3)

        type(scalar_field), dimension(sys_size), intent(in) :: q_coarse
        integer, intent(in)                                 :: bl(3), bh(3), o1, o2, o3
        integer                                             :: i, g1, g2, g3, bl1, bl2, bl3, bh1, bh2, bh3, coff1, coff2, coff3

        ! scalar copies: no host array may be referenced inside the device region (nvfortran/Cray demand it present)

        bl1 = bl(1); bh1 = bh(1); bl2 = bl(2); bh2 = bh(2); bl3 = bl(3); bh3 = bh(3)
        coff1 = amr_cpat_off(1); coff2 = amr_cpat_off(2); coff3 = amr_cpat_off(3)
        $:GPU_PARALLEL_LOOP(collapse=4)
        do i = 1, sys_size
            do g3 = bl3, bh3
                do g2 = bl2, bh2
                    do g1 = bl1, bh1
                        amr_cg(i)%sf(g1 - coff1, g2 - coff2, g3 - coff3) = q_coarse(i)%sf(g1 - o1, g2 - o2, g3 - o3)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_gather_own_box_device

    !> Runtime device pack of the overlap box [bl:bh] global from q_coarse (device) into the contiguous wire buffer buf (host, via
    !! copyout); only the box crosses PCIe, not the full field. Explicit-loop linear buf indexing (g1 fastest, then g2, g3, i) and
    !! the wp cast match the host pack in s_amr_gather_coarse_patch element-for-element, so the receiver's unpack is layout- and
    !! byte-identical (same discipline as s_amr_fine_slice: no array-section syntax near the device map). Twin
    !! s_amr_pack_box_pbmv_device (q<->pb/mv): same wire linear order + wp cast; keep them in lockstep.
    impure subroutine s_amr_pack_box_device(q_coarse, bl, bh, o1, o2, o3, buf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_coarse
        integer, intent(in)                                 :: bl(3), bh(3), o1, o2, o3
        real(wp), intent(inout), contiguous                 :: buf(:)
        integer                                             :: i, g1, g2, g3, bl1, bl2, bl3, bh1, bh2, bh3, n1, n2, n3

        bl1 = bl(1); bh1 = bh(1); bl2 = bl(2); bh2 = bh(2); bl3 = bl(3); bh3 = bh(3)
        n1 = bh1 - bl1 + 1; n2 = bh2 - bl2 + 1; n3 = bh3 - bl3 + 1
        $:GPU_PARALLEL_LOOP(collapse=4, copyout='[buf]')
        do i = 1, sys_size
            do g3 = bl3, bh3
                do g2 = bl2, bh2
                    do g1 = bl1, bh1
                        buf(1 + (g1 - bl1) + n1*((g2 - bl2) + n2*((g3 - bl3) + n3*(i - 1)))) = real(q_coarse(i)%sf(g1 - o1, &
                            & g2 - o2, g3 - o3), wp)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_pack_box_device

    !> Runtime device unpack of a received overlap box [bl:bh] global from the contiguous wire buffer buf (host, via copyin) into
    !! amr_cg (device) in the patch-local frame; only the box crosses PCIe. Same linear order and stp cast as the host unpack in
    !! s_amr_gather_coarse_patch. Twin s_amr_unpack_box_pbmv_device (q<->pb/mv): same wire linear order + stp cast; keep them in
    !! lockstep.
    impure subroutine s_amr_unpack_box_device(bl, bh, buf)

        integer, intent(in)              :: bl(3), bh(3)
        real(wp), intent(in), contiguous :: buf(:)
        integer                          :: i, g1, g2, g3, bl1, bl2, bl3, bh1, bh2, bh3, n1, n2, n3, coff1, coff2, coff3

        bl1 = bl(1); bh1 = bh(1); bl2 = bl(2); bh2 = bh(2); bl3 = bl(3); bh3 = bh(3)
        n1 = bh1 - bl1 + 1; n2 = bh2 - bl2 + 1; n3 = bh3 - bl3 + 1
        coff1 = amr_cpat_off(1); coff2 = amr_cpat_off(2); coff3 = amr_cpat_off(3)
        $:GPU_PARALLEL_LOOP(collapse=4, copyin='[buf]')
        do i = 1, sys_size
            do g3 = bl3, bh3
                do g2 = bl2, bh2
                    do g1 = bl1, bh1
                        amr_cg(i)%sf(g1 - coff1, g2 - coff2, &
                               & g3 - coff3) = real(buf(1 + (g1 - bl1) + n1*((g2 - bl2) + n2*((g3 - bl3) + n3*(i - 1)))), stp)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_unpack_box_device

    !> Runtime device own-box copy for the pbmv gather: pb/mv (device) -> amr_cg_pb/mv (device) over [bl:bh] global in the
    !! patch-local frame. Same index map and direct stp assignment as the host path in s_amr_gather_coarse_patch_pbmv. Twin
    !! s_amr_gather_own_box_device (pb/mv<->q): q_cons sibling of this own-box copy; keep the index map in lockstep.
    impure subroutine s_amr_gather_own_box_pbmv_device(pb_coarse, mv_coarse, bl, bh, o1, o2, o3)

        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(in) :: pb_coarse, mv_coarse
        integer, intent(in) :: bl(3), bh(3), o1, o2, o3
        integer :: q, ib_, g1, g2, g3, bl1, bl2, bl3, bh1, bh2, bh3, coff1, coff2, coff3

        bl1 = bl(1); bh1 = bh(1); bl2 = bl(2); bh2 = bh(2); bl3 = bl(3); bh3 = bh(3)
        coff1 = amr_cpat_off(1); coff2 = amr_cpat_off(2); coff3 = amr_cpat_off(3)
        $:GPU_PARALLEL_LOOP(collapse=5)
        do ib_ = 1, nb
            do q = 1, nnode
                do g3 = bl3, bh3
                    do g2 = bl2, bh2
                        do g1 = bl1, bh1
                            amr_cg_pb(g1 - coff1, g2 - coff2, g3 - coff3, q, ib_) = pb_coarse(g1 - o1, g2 - o2, g3 - o3, q, ib_)
                            amr_cg_mv(g1 - coff1, g2 - coff2, g3 - coff3, q, ib_) = mv_coarse(g1 - o1, g2 - o2, g3 - o3, q, ib_)
                        end do
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_gather_own_box_pbmv_device

    !> Runtime device pack for the pbmv gather: pb block then mv block of the overlap box [bl:bh] global into the contiguous wire
    !! buffer buf (host, via copyout). Linear order (g1 fastest, then g2, g3, q, ib_; mv offset by half the message) and wp cast
    !! match the host pack in s_amr_gather_coarse_patch_pbmv element-for-element. Twin s_amr_pack_box_device (pb/mv<->q): q_cons
    !! sibling; keep the wire linear order + wp cast in lockstep.
    impure subroutine s_amr_pack_box_pbmv_device(pb_coarse, mv_coarse, bl, bh, o1, o2, o3, buf)

        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(in) :: pb_coarse, mv_coarse
        integer, intent(in) :: bl(3), bh(3), o1, o2, o3
        real(wp), intent(inout), contiguous :: buf(:)
        integer :: q, ib_, g1, g2, g3, bl1, bl2, bl3, bh1, bh2, bh3, n1, n2, n3, half

        bl1 = bl(1); bh1 = bh(1); bl2 = bl(2); bh2 = bh(2); bl3 = bl(3); bh3 = bh(3)
        n1 = bh1 - bl1 + 1; n2 = bh2 - bl2 + 1; n3 = bh3 - bl3 + 1
        half = n1*n2*n3*nnode*nb
        $:GPU_PARALLEL_LOOP(collapse=5, copyout='[buf]')
        do ib_ = 1, nb
            do q = 1, nnode
                do g3 = bl3, bh3
                    do g2 = bl2, bh2
                        do g1 = bl1, bh1
                            buf(1 + (g1 - bl1) + n1*((g2 - bl2) + n2*((g3 - bl3) + n3*((q - 1) + nnode*(ib_ - 1))))) &
                                & = real(pb_coarse(g1 - o1, g2 - o2, g3 - o3, q, ib_), wp)
                            buf(half + 1 + (g1 - bl1) + n1*((g2 - bl2) + n2*((g3 - bl3) + n3*((q - 1) + nnode*(ib_ - 1))))) &
                                & = real(mv_coarse(g1 - o1, g2 - o2, g3 - o3, q, ib_), wp)
                        end do
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_pack_box_pbmv_device

    !> Runtime device unpack for the pbmv gather: received wire buffer buf (host, via copyin) -> amr_cg_pb/mv (device) over [bl:bh]
    !! global in the patch-local frame. Same linear order and stp cast as the host unpack in s_amr_gather_coarse_patch_pbmv. Twin
    !! s_amr_unpack_box_device (pb/mv<->q): q_cons sibling; keep the wire linear order + stp cast in lockstep.
    impure subroutine s_amr_unpack_box_pbmv_device(bl, bh, buf)

        integer, intent(in)              :: bl(3), bh(3)
        real(wp), intent(in), contiguous :: buf(:)
        integer                          :: q, ib_, g1, g2, g3, bl1, bl2, bl3, bh1, bh2, bh3, n1, n2, n3, half, coff1, coff2, coff3

        bl1 = bl(1); bh1 = bh(1); bl2 = bl(2); bh2 = bh(2); bl3 = bl(3); bh3 = bh(3)
        n1 = bh1 - bl1 + 1; n2 = bh2 - bl2 + 1; n3 = bh3 - bl3 + 1
        half = n1*n2*n3*nnode*nb
        coff1 = amr_cpat_off(1); coff2 = amr_cpat_off(2); coff3 = amr_cpat_off(3)
        $:GPU_PARALLEL_LOOP(collapse=5, copyin='[buf]')
        do ib_ = 1, nb
            do q = 1, nnode
                do g3 = bl3, bh3
                    do g2 = bl2, bh2
                        do g1 = bl1, bh1
                            amr_cg_pb(g1 - coff1, g2 - coff2, g3 - coff3, q, &
                                      & ib_) = real(buf(1 + (g1 - bl1) + n1*((g2 - bl2) + n2*((g3 - bl3) + n3*((q - 1) &
                                      & + nnode*(ib_ - 1))))), stp)
                            amr_cg_mv(g1 - coff1, g2 - coff2, g3 - coff3, q, &
                                      & ib_) = real(buf(half + 1 + (g1 - bl1) + n1*((g2 - bl2) + n2*((g3 - bl3) + n3*((q - 1) &
                                      & + nnode*(ib_ - 1))))), stp)
                        end do
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_unpack_box_pbmv_device

    !> Non-polytropic QBMM: piecewise-constant prolongation of the fine pb/mv ghost shell from the coarse side-state (device kernel;
    !! interior untouched). The ghosts feed the widened-idwint conversions and the qbmm rhs over the shell, mirroring the q_cons
    !! ghost fill. All four arrays are assumed-shape dummies with %sf pointer-member actuals (the pb_ts pattern); raw derived-type
    !! 5D members as actuals trip nvfortran's component-section data clauses on device. Twin s_amr_fill_fine_ghosts (pb/mv<->q):
    !! q_cons sibling; keep the ghost-fill mapping in lockstep.
    impure subroutine s_amr_fill_fine_ghosts_pbmv(pb_c, mv_c, pb_t, mv_t)

        !> coarse pb/mv read from the gathered block-local patch amr_cg_pb/mv (0-based patch frame, cell 0 == amr_cpat_off): the
        !! callers run s_amr_gather_coarse_patch_pbmv on all ranks first, so np>=2 reads the correct coarse rank's side-state
        real(stp), dimension(0:,0:,0:,1:,1:), intent(in) :: pb_c, mv_c

        real(stp), dimension(amr_slots(amr_cur)%idwbuff(1)%beg:,amr_slots(amr_cur)%idwbuff(2)%beg:, &
             & amr_slots(amr_cur)%idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_t, mv_t
        integer               :: fi, fj, fk, q, ib_, ci, cj, ck, rr, lo1, lo2, lo3, ox, oy, oz
        integer               :: s, ns, ss, g, r, n1, n2, stot
        integer, dimension(6) :: sb1, se1, sb2, se2, sb3, se3, soff, scnt
        logical               :: d2, d3

        ox = amr_cpat_off(1); oy = amr_cpat_off(2); oz = amr_cpat_off(3)
        d2 = n_glb > 0; d3 = p_glb > 0
        rr = amr_slots(amr_cur)%amr_ref_ratio
        lo1 = amr_isect_lo(1); lo2 = amr_isect_lo(2); lo3 = amr_isect_lo(3)
        call s_amr_build_ghost_slabs(ns, sb1, se1, sb2, se2, sb3, se3)
        ! flat index over the concatenated disjoint slabs, one kernel instead of ns; see s_amr_fill_fine_ghosts
        soff(1) = 0
        do s = 1, ns
            scnt(s) = (se1(s) - sb1(s) + 1)*(se2(s) - sb2(s) + 1)*(se3(s) - sb3(s) + 1)
            if (s < ns) soff(s + 1) = soff(s) + scnt(s)
        end do
        stot = soff(ns) + scnt(ns)
        $:GPU_PARALLEL_LOOP(collapse=3, copyin='[sb1, se1, sb2, se2, sb3, se3, soff, scnt]', private='[s, ss, r, n1, n2, fi, fj, &
                            & fk, ci, cj, ck]')
        do ib_ = 1, nb
            do q = 1, nnode
                do g = 0, stot - 1
                    s = 1
                    do ss = 2, ns
                        if (g >= soff(ss)) s = ss
                    end do
                    r = g - soff(s)
                    n1 = se1(s) - sb1(s) + 1; n2 = se2(s) - sb2(s) + 1
                    fi = sb1(s) + mod(r, n1)
                    fj = sb2(s) + mod(r/n1, n2)
                    fk = sb3(s) + r/(n1*n2)
                    ck = 0
                    if (d3) ck = lo3 + floor(real(fk, wp)/real(rr, wp)) - oz
                    cj = 0
                    if (d2) cj = lo2 + floor(real(fj, wp)/real(rr, wp)) - oy
                    ci = lo1 + floor(real(fi, wp)/real(rr, wp)) - ox
                    pb_t(fi, fj, fk, q, ib_) = pb_c(ci, cj, ck, q, ib_)
                    mv_t(fi, fj, fk, q, ib_) = mv_c(ci, cj, ck, q, ib_)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_fill_fine_ghosts_pbmv

    !> Decompose the current fine block's ghost shell (buffered extent minus interior) into ns disjoint face slabs whose union is
    !! exactly the non-interior cells, so the ghost-fill kernels do O(surface) work instead of masking the full buffered volume. x
    !! slabs span the full transverse extent; y slabs restrict x to the interior; z slabs restrict x and y. Collapsed dims
    !! (n_glb/p_glb == 0) contribute no slabs.
    pure subroutine s_amr_build_ghost_slabs(ns, sb1, se1, sb2, se2, sb3, se3)

        integer, intent(out)               :: ns
        integer, dimension(6), intent(out) :: sb1, se1, sb2, se2, sb3, se3
        integer                            :: fm, fn, fp, b1, e1, b2, e2, b3, e3

        fm = amr_slots(amr_cur)%m; fn = amr_slots(amr_cur)%n; fp = amr_slots(amr_cur)%p
        b1 = amr_slots(amr_cur)%idwbuff(1)%beg; e1 = amr_slots(amr_cur)%idwbuff(1)%end
        b2 = amr_slots(amr_cur)%idwbuff(2)%beg; e2 = amr_slots(amr_cur)%idwbuff(2)%end
        b3 = amr_slots(amr_cur)%idwbuff(3)%beg; e3 = amr_slots(amr_cur)%idwbuff(3)%end
        ns = 2
        sb1(1) = b1; se1(1) = -1; sb1(2) = fm + 1; se1(2) = e1
        sb2(1:2) = b2; se2(1:2) = e2; sb3(1:2) = b3; se3(1:2) = e3
        if (n_glb > 0) then
            ns = 4
            sb2(3) = b2; se2(3) = -1; sb2(4) = fn + 1; se2(4) = e2
            sb1(3:4) = 0; se1(3:4) = fm; sb3(3:4) = b3; se3(3:4) = e3
        end if
        if (p_glb > 0) then
            ns = 6
            sb3(5) = b3; se3(5) = -1; sb3(6) = fp + 1; se3(6) = e3
            sb1(5:6) = 0; se1(5:6) = fm; sb2(5:6) = 0; se2(5:6) = fn
        end if

    end subroutine s_amr_build_ghost_slabs

    !> Fill the fine ghost shell by conservative-linear prolongation from q_coarse, the gathered block-local coarse patch amr_cg
    !! (fine-level distribution; the caller gathers the source first). Device kernel: reads the patch and writes the fine target in
    !! device memory. floor/modulo mapping is valid for negative fine indices (ghosts). Interior untouched. Multi-fluid volume
    !! fractions get the same sum-preserving closure as the interior prolongation (second kernel). Twin s_amr_fill_fine_ghosts_pbmv
    !! (q<->pb/mv): pb/mv sibling; keep the mapping in lockstep.
    !!
    !! The body is generated from a Fypp accessor lambda (the idiom m_riemann_solver_hlld uses for its per-direction stencil
    !! variants) so the write target is fixed at preprocessing time: a dummy referenced in any branch of a target region is still
    !! mapped, and each mapped array costs per launch. `_cons` writes the conserved store at dense local index `loc`.
    #:for SFX, TGT in [('cons', 'amr_cons_st')]
        #:set QF = lambda ix: TGT + '(fi, fj, fk, ' + ix + ', loc)'
        impure subroutine s_amr_fill_fine_ghosts_${SFX}$(q_coarse, loc)

            type(scalar_field), dimension(sys_size), intent(in) :: q_coarse
            integer, intent(in)                                 :: loc
            integer                                             :: i, fi, fj, fk, ci, cj, ck, ox, oy, oz
            integer                                             :: rr, lo1, lo2, lo3
            integer                                             :: advb, adve, bbeg, bend, bstride
            integer                                             :: s, ns
            integer                                             :: ss, g, r, n1, n2, stot
            integer, dimension(6)                               :: sb1, se1, sb2, se2, sb3, se3, soff, scnt
            logical                                             :: d2, d3, multi, shx, shy, shz, bubEE
            real(wp)                                            :: u0, sx, sy, sz, xix, xiy, xiz, av, asum

            ! q_coarse is the gathered block-local patch amr_cg (fine-level distribution); amr_isect_lo (global, == region_lo on
            ! the owner) + f/rr - amr_cpat_off is the patch-local coarse index. Fine indices are local to this block.

            ox = amr_cpat_off(1); oy = amr_cpat_off(2); oz = amr_cpat_off(3)
            d2 = n_glb > 0; d3 = p_glb > 0
            rr = amr_slots(amr_cur)%amr_ref_ratio
            lo1 = amr_isect_lo(1); lo2 = amr_isect_lo(2); lo3 = amr_isect_lo(3)
            multi = num_fluids > 1 .and. (.not. bubbles_lagrange)  ! EL alphas sum to beta, not 1: no sum-to-one closure
            advb = eqn_idx%adv%beg; adve = eqn_idx%adv%end
            bubEE = bubbles_euler; bbeg = eqn_idx%bub%beg; bend = eqn_idx%bub%end
            bstride = 1; if (bubEE) bstride = (bend - bbeg + 1)/nb
            call s_amr_build_ghost_slabs(ns, sb1, se1, sb2, se2, sb3, se3)
            ! One kernel over the concatenation of the ns face slabs instead of one kernel each. The slabs are disjoint and their
            ! union
            ! is exactly the ghost shell (s_amr_build_ghost_slabs), so every ghost cell is written exactly once and the result is
            ! independent
            ! of how the flat index is ordered. Not the padded-hull form of s_amr_capture_creg_dense_batch: the x
            ! slabs
            ! span the full transverse extent, so a hull over all slabs is the whole buffered volume and masking it would throw away
            ! the
            ! O(surface) decomposition this routine exists to get.
            soff(1) = 0
            do s = 1, ns
                scnt(s) = (se1(s) - sb1(s) + 1)*(se2(s) - sb2(s) + 1)*(se3(s) - sb3(s) + 1)
                if (s < ns) soff(s + 1) = soff(s) + scnt(s)
            end do
            stot = soff(ns) + scnt(ns)
            amr_slab_tab(1,:) = sb1; amr_slab_tab(2,:) = se1; amr_slab_tab(3,:) = sb2; amr_slab_tab(4,:) = se2
            amr_slab_tab(5,:) = sb3; amr_slab_tab(6,:) = se3; amr_slab_tab(7,:) = soff; amr_slab_tab(8,:) = scnt
            $:GPU_UPDATE(device='[amr_slab_tab]')
            $:GPU_PARALLEL_LOOP(collapse=2, private='[s, ss, r, n1, n2, fi, fj, fk, ci, cj, ck, xix, xiy, xiz, u0, sx, sy, sz]')
            do i = 1, sys_size
                do g = 0, stot - 1
                    s = 1  ! decode the flat index: ns <= 6, so a scan beats storing a per-cell slab map
                    do ss = 2, ns
                        if (g >= amr_slab_tab(7, ss)) s = ss
                    end do
                    r = g - amr_slab_tab(7, s)
                    n1 = amr_slab_tab(2, s) - amr_slab_tab(1, s) + 1; n2 = amr_slab_tab(4, s) - amr_slab_tab(3, s) + 1
                    fi = amr_slab_tab(1, s) + mod(r, n1)
                    fj = amr_slab_tab(3, s) + mod(r/n1, n2)
                    fk = amr_slab_tab(5, s) + r/(n1*n2)
                    ! the slabs cover exactly the ghost shell; multi-fluid, skip the volume fractions (closure kernel below)
                    if (.not. (multi .and. i >= advb .and. i <= adve)) then
                        ck = 0; xiz = 0._wp
                        if (d3) then
                            ck = lo3 + floor(real(fk, wp)/real(rr, wp)) - oz
                            xiz = (real(modulo(fk, rr), wp) - real(rr - 1, wp)*0.5_wp)/real(rr, wp)
                        end if
                        cj = 0; xiy = 0._wp
                        if (d2) then
                            cj = lo2 + floor(real(fj, wp)/real(rr, wp)) - oy
                            xiy = (real(modulo(fj, rr), wp) - real(rr - 1, wp)*0.5_wp)/real(rr, wp)
                        end if
                        ci = lo1 + floor(real(fi, wp)/real(rr, wp)) - ox
                        xix = (real(modulo(fi, rr), wp) - real(rr - 1, wp)*0.5_wp)/real(rr, wp)
                        u0 = real(q_coarse(i)%sf(ci, cj, ck), wp)
                        sx = minmod(real(q_coarse(i)%sf(ci + 1, cj, ck), wp) - u0, u0 - real(q_coarse(i)%sf(ci - 1, cj, ck), wp))
                        sy = 0._wp
                        if (d2) sy = minmod(real(q_coarse(i)%sf(ci, cj + 1, ck), wp) - u0, u0 - real(q_coarse(i)%sf(ci, cj - 1, &
                            & ck), wp))
                        sz = 0._wp
                        if (d3) sz = minmod(real(q_coarse(i)%sf(ci, cj, ck + 1), wp) - u0, u0 - real(q_coarse(i)%sf(ci, cj, &
                            & ck - 1), wp))
                        ! QBMM: inject the bub block piecewise-constant (child = u0) so the ghost inherits the coarse cell's
                        ! realizable 6-moment set (CHyQMOM needs variance c20 > 0; per-component minmod slopes would break
                        ! that joint constraint). Non-QBMM Euler-Euler bubbles instead floor their positive moments (nR /
                        ! npb / nmv); the signed velocity moment nV (offset 1) is skipped.
                        if (qbmm .and. i >= bbeg .and. i <= bend) then
                            sx = 0._wp; sy = 0._wp; sz = 0._wp
                        end if
                        ${QF('i')}$ = u0 + sx*xix + sy*xiy + sz*xiz
                        if (bubEE .and. .not. qbmm .and. i >= bbeg .and. i <= bend) then
                            if (mod(i - bbeg, bstride) /= 1) ${QF('i')}$ = max(real(${QF('i')}$, wp), bub_pos_frac*u0)
                        end if
                    end if
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            ! multi-fluid volume-fraction ghosts: per-cell closure mirroring s_prolong_alphas_closure (shared limiter switch over
            ! all
            ! fluids; interpolate + clamp fluids advb..adve-1; alpha_n = 1 - sum)
            if (multi) then
                ! same flat-index fusion as the prolongation loop above, over the same disjoint slabs
                $:GPU_PARALLEL_LOOP(private='[s, ss, r, n1, n2, fi, fj, fk, i, ci, cj, ck, xix, xiy, xiz, u0, sx, sy, sz, av, &
                                    & asum, shx, shy, shz]')
                do g = 0, stot - 1
                    s = 1
                    do ss = 2, ns
                        if (g >= amr_slab_tab(7, ss)) s = ss
                    end do
                    r = g - amr_slab_tab(7, s)
                    n1 = amr_slab_tab(2, s) - amr_slab_tab(1, s) + 1; n2 = amr_slab_tab(4, s) - amr_slab_tab(3, s) + 1
                    fi = amr_slab_tab(1, s) + mod(r, n1)
                    fj = amr_slab_tab(3, s) + mod(r/n1, n2)
                    fk = amr_slab_tab(5, s) + r/(n1*n2)
                    ck = 0; xiz = 0._wp
                    if (d3) then
                        ck = lo3 + floor(real(fk, wp)/real(rr, wp)) - oz
                        xiz = (real(modulo(fk, rr), wp) - real(rr - 1, wp)*0.5_wp)/real(rr, wp)
                    end if
                    cj = 0; xiy = 0._wp
                    if (d2) then
                        cj = lo2 + floor(real(fj, wp)/real(rr, wp)) - oy
                        xiy = (real(modulo(fj, rr), wp) - real(rr - 1, wp)*0.5_wp)/real(rr, wp)
                    end if
                    ci = lo1 + floor(real(fi, wp)/real(rr, wp)) - ox
                    xix = (real(modulo(fi, rr), wp) - real(rr - 1, wp)*0.5_wp)/real(rr, wp)
                    shx = .true.; shy = d2; shz = d3
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = advb, adve
                        u0 = real(q_coarse(i)%sf(ci, cj, ck), wp)
                        if ((real(q_coarse(i)%sf(ci + 1, cj, ck), wp) - u0)*(u0 - real(q_coarse(i)%sf(ci - 1, cj, ck), &
                            & wp)) <= 0._wp) shx = .false.
                        if (d2) then
                            if ((real(q_coarse(i)%sf(ci, cj + 1, ck), wp) - u0)*(u0 - real(q_coarse(i)%sf(ci, cj - 1, ck), &
                                & wp)) <= 0._wp) shy = .false.
                        end if
                        if (d3) then
                            if ((real(q_coarse(i)%sf(ci, cj, ck + 1), wp) - u0)*(u0 - real(q_coarse(i)%sf(ci, cj, ck - 1), &
                                & wp)) <= 0._wp) shz = .false.
                        end if
                    end do
                    asum = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = advb, adve - 1
                        u0 = real(q_coarse(i)%sf(ci, cj, ck), wp)
                        sx = 0._wp
                        if (shx) sx = minmod(real(q_coarse(i)%sf(ci + 1, cj, ck), wp) - u0, u0 - real(q_coarse(i)%sf(ci - 1, cj, &
                            & ck), wp))
                        sy = 0._wp
                        if (shy) sy = minmod(real(q_coarse(i)%sf(ci, cj + 1, ck), wp) - u0, u0 - real(q_coarse(i)%sf(ci, cj - 1, &
                            & ck), wp))
                        sz = 0._wp
                        if (shz) sz = minmod(real(q_coarse(i)%sf(ci, cj, ck + 1), wp) - u0, u0 - real(q_coarse(i)%sf(ci, cj, &
                            & ck - 1), wp))
                        av = min(max(u0 + sx*xix + sy*xiy + sz*xiz, 0._wp), 1._wp)
                        ${QF('i')}$ = av
                        asum = asum + av
                    end do
                    ${QF('adve')}$ = 1._wp - asum
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

        end subroutine s_amr_fill_fine_ghosts_${SFX}$
    #:endfor

    !> Exchange the coarse conservative ghost layers at internal rank boundaries (physical-boundary ghosts untouched; per direction
    !! beg then end, mirroring s_populate_variables_buffers' disblock). The solver never fills cons ghosts (only prim), so ranks
    !! whose fine ghost-fill or prolongation stencil leaves their interior need this first. All ranks must call together (pairwise
    !! exchange per internal neighbor).
    impure subroutine s_amr_exchange_coarse_cons_halo(q_cons)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons

        if (bc_x%beg >= 0) call s_mpi_sendrecv_variables_buffers(q_cons, 1, -1, sys_size)
        if (bc_x%end >= 0) call s_mpi_sendrecv_variables_buffers(q_cons, 1, 1, sys_size)
        if (n_glb > 0) then
            if (bc_y%beg >= 0) call s_mpi_sendrecv_variables_buffers(q_cons, 2, -1, sys_size)
            if (bc_y%end >= 0) call s_mpi_sendrecv_variables_buffers(q_cons, 2, 1, sys_size)
        end if
        if (p_glb > 0) then
            if (bc_z%beg >= 0) call s_mpi_sendrecv_variables_buffers(q_cons, 3, -1, sys_size)
            if (bc_z%end >= 0) call s_mpi_sendrecv_variables_buffers(q_cons, 3, 1, sys_size)
        end if

    end subroutine s_amr_exchange_coarse_cons_halo

    !> True iff dim d is globally periodic. Uses l0_periodic (periodic_bc allreduced to all ranks in s_l0_tiles_init); periodic_bc
    !! itself is captured from the original bc before MFC folds periodicity into the MPI-cart topology, but only on rank 0, so the
    !! allreduced copy is the one that is consistent on every rank (required since f_amr_seam builds a replicated seam list).
    pure logical function f_l0_dim_periodic(d) result(per)

        integer, intent(in) :: d

        per = l0_periodic(d)

    end function f_l0_dim_periodic

    !> True iff sub-block yb sits immediately above sub-block xb on xb's high face in dim d: yb's low coarse face is xb's high face
    !! + 1, and (tiling produces a regular grid) they share the transverse coarse extents exactly. Each fine-fine seam is exactly
    !! one such ordered (xb, yb) pair (the lower block is xb). For an L0-tile periodic dim there is also a wrap-seam: xb at the
    !! domain high face (region_hi == gcell) and yb at the domain low face (region_lo == 0) with matching transverse; the fine-fine
    !! halo then fills xb's high ghost from yb's low interior and vice versa, exactly the periodic connection. Gated on l0_ntile>0
    !! so the AMR fine-block path (blocks never touch the domain edge) is unaffected.
    pure logical function f_amr_seam(xb, yb, d) result(seam)

        integer, intent(in) :: xb, yb, d
        integer             :: t, gc
        logical             :: adj

        adj = amr_region_lo_all(d, yb) == amr_region_hi_all(d, xb) + 1
        if (l0_ntile > 0 .and. f_l0_dim_periodic(d)) then
            gc = merge(m_glb, merge(n_glb, p_glb, d == 2), d == 1)
            adj = adj .or. (amr_region_hi_all(d, xb) == gc .and. amr_region_lo_all(d, yb) == 0)
        end if
        seam = adj
        do t = 1, 3
            if (t /= d) seam = seam .and. amr_region_lo_all(t, xb) == amr_region_lo_all(t, yb) .and. amr_region_hi_all(t, &
                & xb) == amr_region_hi_all(t, yb)
        end do

    end function f_amr_seam

    !> Same-rank seam exchange for every same-rank seam pair on this rank, both directions in one device kernel: xb's high near-seam
    !! interior -> yb's low seam ghost, and yb's low interior -> xb's high ghost. Fusing the two directions is safe because all four
    !! slabs are disjoint (a block's seam ghost lies outside its own interior, and the two reads are from different slots than the
    !! two writes), so neither direction can observe the other's store. Both blocks of a pair are addressed by their flat-store
    !! slot, so a runtime pair index is a plain subscript. Index order matches s_amr_fine_slice. The per-pair seam dim varies, so
    !! the index decode is a runtime select. Threads are launched over the max pair extent and masked, rather than prefix-summed:
    !! pair sizes are equal under uniform tiling, so the waste is ~0 and it avoids an O(npair) per-thread offset search.
    impure subroutine s_amr_fine_seam_exchange(npair, plx, ply, pd, pxhi, pfm, ndep)

        integer, intent(in) :: npair
        integer, intent(in) :: plx(:), ply(:), pd(:), pxhi(:), pfm(:,:)  !< per-pair: slots, seam dim, high extent, fine extents
        integer, intent(in) :: ndep
        integer             :: i, a, b, t, na, nb, pr, g, gmax, lx, ly, d, xhi, cnt
        integer             :: ix1, ix2, ix3, iy1, iy2, iy3

        ! max thread extent over the pairs on this rank (transverse product x seam depth)

        gmax = 0
        do pr = 1, npair
            select case (pd(pr))
            case (1); na = pfm(2, pr) + 1; nb = pfm(3, pr) + 1
            case (2); na = pfm(1, pr) + 1; nb = pfm(3, pr) + 1
            case default; na = pfm(1, pr) + 1; nb = pfm(2, pr) + 1
            end select
            gmax = max(gmax, na*nb*ndep)
        end do

        $:GPU_PARALLEL_LOOP(collapse=3, copyin='[plx, ply, pd, pxhi, pfm]', private='[lx, ly, d, xhi, na, nb, cnt, a, b, t, ix1, &
                            & ix2, ix3, iy1, iy2, iy3]')
        do pr = 1, npair
            do i = 1, sys_size
                do g = 0, gmax - 1
                    d = pd(pr)
                    if (d == 1) then
                        na = pfm(2, pr) + 1; nb = pfm(3, pr) + 1
                    else if (d == 2) then
                        na = pfm(1, pr) + 1; nb = pfm(3, pr) + 1
                    else
                        na = pfm(1, pr) + 1; nb = pfm(2, pr) + 1
                    end if
                    cnt = na*nb*ndep
                    if (g < cnt) then
                        lx = plx(pr); ly = ply(pr); xhi = pxhi(pr)
                        a = mod(g, na); b = mod(g/na, nb); t = g/(na*nb)
                        ! (a, b) are the transverse indices, t the depth into the seam; place them per the pair's seam dim
                        if (d == 1) then
                            ix1 = xhi - ndep + 1 + t; ix2 = a; ix3 = b
                            iy1 = -ndep + t; iy2 = a; iy3 = b
                        else if (d == 2) then
                            ix1 = a; ix2 = xhi - ndep + 1 + t; ix3 = b
                            iy1 = a; iy2 = -ndep + t; iy3 = b
                        else
                            ix1 = a; ix2 = b; ix3 = xhi - ndep + 1 + t
                            iy1 = a; iy2 = b; iy3 = -ndep + t
                        end if
                        ! xb high interior -> yb low ghost, then yb low interior -> xb high ghost. All four slabs are disjoint (a
                        ! block's seam ghost lies outside its own interior, and the reads are from different slots than the
                        ! writes), so fusing the two directions cannot let one observe the other's store.
                        amr_cons_st(iy1, iy2, iy3, i, ly) = amr_cons_st(ix1, ix2, ix3, i, lx)
                        if (d == 1) then
                            ix1 = xhi + 1 + t; iy1 = t
                        else if (d == 2) then
                            ix2 = xhi + 1 + t; iy2 = t
                        else
                            ix3 = xhi + 1 + t; iy3 = t
                        end if
                        amr_cons_st(ix1, ix2, ix3, i, lx) = amr_cons_st(iy1, iy2, iy3, i, ly)
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_fine_seam_exchange

    !> Seam dimension of the ordered pair (xb, yb): the dim d in which yb is the immediate high-face neighbour of xb at matched
    !! resolution (same level), or 0 if not a same-level fine-fine seam. Face adjacency requires transverse overlap, so a pair is
    !! adjacent in at most one dim, so the last-true assignment is unambiguous.
    pure integer function f_amr_seam_dim(xb, yb) result(d)

        integer, intent(in) :: xb, yb

        d = 0
        if (xb == yb) return
        if (amr_block_level(xb) /= amr_block_level(yb)) return
        if (f_amr_seam(xb, yb, 1)) d = 1
        if (n_glb > 0) then; if (f_amr_seam(xb, yb, 2)) d = 2; end if
        if (p_glb > 0) then; if (f_amr_seam(xb, yb, 3)) d = 3; end if

    end function f_amr_seam_dim

    !> Rebuild the cached same-level seam-pair list (amr_seam_pairs) once per regrid/restart rather than every RK stage. Same (xb,
    !! yb) order on all ranks (replicated region metadata) so the paired seam transfers stay matched. Count then fill for an
    !! exact-size list (no cap, no overflow). Also rebuilds the per-block gather/scatter overlap-rank lists (amr_ovl_gather/scatter)
    !! by O(overlap) inversion of the computed decomposition (s_amr_ranks_overlapping), sized to the max overlap, with no
    !! O(num_procs) scan or table.
    !> Binary-search the Morton-sorted block lo corners (ord/mkey) for the block whose region lo equals clo, verify the full
    !! same-level seam predicate against xb, and record it in (mb, md, nm). Blocks are disjoint, so at most one block carries a
    !! given lo corner at a level; f_amr_seam_dim takes the last true dim, so a block already recorded is raised to the higher d
    !! rather than duplicated.
    pure subroutine s_amr_seam_probe(xb, d, clo, mkey, ord, nblk, mb, md, nm)

        integer, intent(in)         :: xb, d, clo(3), nblk
        integer(kind=8), intent(in) :: mkey(:)
        integer, intent(in)         :: ord(:)
        integer, intent(inout)      :: mb(3), md(3), nm
        integer                     :: yb, t, i, lo_s, hi_s, mid_s
        integer(kind=8)             :: ck
        logical                     :: ok

        ck = f_morton(clo(1), clo(2), clo(3))
        lo_s = 1; hi_s = nblk
        do while (lo_s <= hi_s)
            mid_s = (lo_s + hi_s)/2
            if (mkey(ord(mid_s)) < ck) then
                lo_s = mid_s + 1
            else
                hi_s = mid_s - 1
            end if
        end do
        ! lo_s is the first index with key >= ck; scan the equal-key run. f_morton keeps 21 bits/dim, so
        ! above 2**21 cells/dim a run can hold distinct corners; the explicit lo check rejects those.
        do i = lo_s, nblk
            if (mkey(ord(i)) /= ck) exit
            yb = ord(i)
            if (yb == xb) cycle
            if (amr_block_level(yb) /= amr_block_level(xb)) cycle
            ok = .true.
            do t = 1, 3
                if (amr_region_lo_all(t, yb) /= clo(t)) ok = .false.
                if (t /= d) then
                    if (amr_region_hi_all(t, yb) /= amr_region_hi_all(t, xb)) ok = .false.
                end if
            end do
            if (.not. ok) cycle
            do t = 1, nm
                if (mb(t) == yb) then
                    md(t) = max(md(t), d)
                    return
                end if
            end do
            nm = nm + 1
            mb(nm) = yb; md(nm) = d
            return
        end do

    end subroutine s_amr_seam_probe

    impure subroutine s_amr_build_seam_pairs()

        integer                      :: xb, d, np, k, mx, pass, nm, im, jm, tb, td, nb
        integer                      :: plo(3), phi(3), rlo(3), rhi(3)
        integer                      :: mb(3), md(3), clo(3), gc
        integer                      :: width, lo_m, mid_m, hi_m, i_m, j_m, t_m
        integer(kind=8), allocatable :: mkey(:)
        integer, allocatable         :: ord(:), mrg(:)

        if (allocated(amr_seam_pairs)) deallocate (amr_seam_pairs)

        ! Blocks are disjoint, so (level, region lo) names one uniquely, and the seam predicate fixes the
        ! neighbour's lo corner: transverse lo equal to xb's, and lo(d) = hi(d, xb) + 1. The all-pairs
        ! O(nblocks^2) scan is therefore a lookup. Morton-sort the lo corners once, then binary-search the
        ! single candidate per (block, dim) and verify the predicate on it. Emission order is xb ascending,
        ! yb ascending within xb, which the paired seam transfers depend on; a reordered list mismatches
        ! sends to receives and deadlocks.
        nb = max(amr_num_blocks, 1)
        allocate (mkey(nb), ord(nb), mrg(nb))
        do k = 1, amr_num_blocks
            mkey(k) = f_morton(amr_region_lo_all(1, k), amr_region_lo_all(2, k), amr_region_lo_all(3, k))
            ord(k) = k
        end do

        ! Bottom-up stable merge sort by Morton key (same form as s_amr_sfc_cut): a pure function of the
        ! replicated region metadata, so every rank builds the identical order.
        width = 1
        do while (width < amr_num_blocks)
            lo_m = 1
            do while (lo_m <= amr_num_blocks - width)
                mid_m = lo_m + width - 1
                hi_m = min(lo_m + 2*width - 1, amr_num_blocks)
                i_m = lo_m; j_m = mid_m + 1; t_m = lo_m
                do while (i_m <= mid_m .and. j_m <= hi_m)
                    if (mkey(ord(i_m)) <= mkey(ord(j_m))) then
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

        ! pass 1 counts, pass 2 fills: keeps amr_seam_pairs exactly sized
        do pass = 1, 2
            np = 0
            do xb = 1, amr_num_blocks
                nm = 0
                do d = 1, 3
                    if (d == 2 .and. n_glb <= 0) cycle
                    if (d == 3 .and. p_glb <= 0) cycle
                    clo = amr_region_lo_all(:,xb)
                    clo(d) = amr_region_hi_all(d, xb) + 1
                    call s_amr_seam_probe(xb, d, clo, mkey, ord, amr_num_blocks, mb, md, nm)
                    ! periodic wrap (l0 tiling only): xb on the domain high face pairs with lo(d) = 0
                    if (l0_ntile > 0) then
                        if (f_l0_dim_periodic(d)) then
                            gc = merge(m_glb, merge(n_glb, p_glb, d == 2), d == 1)
                            if (amr_region_hi_all(d, xb) == gc) then
                                clo(d) = 0
                                call s_amr_seam_probe(xb, d, clo, mkey, ord, amr_num_blocks, mb, md, nm)
                            end if
                        end if
                    end if
                end do
                ! ascending yb within xb (nm <= 3)
                do im = 2, nm
                    tb = mb(im); td = md(im); jm = im - 1
                    do while (jm >= 1)
                        if (mb(jm) <= tb) exit
                        mb(jm + 1) = mb(jm); md(jm + 1) = md(jm); jm = jm - 1
                    end do
                    mb(jm + 1) = tb; md(jm + 1) = td
                end do
                do im = 1, nm
                    np = np + 1
                    if (pass == 2) then
                        amr_seam_pairs(1, np) = xb; amr_seam_pairs(2, np) = mb(im); amr_seam_pairs(3, np) = md(im)
                    end if
                end do
            end do
            if (pass == 1) then
                amr_num_seam_pairs = np
                allocate (amr_seam_pairs(3, max(np, 1)))
            end if
        end do
        deallocate (mkey, ord, mrg)
        ! per-block P2P overlap-rank lists by O(overlap) inversion (gather: rank coarse range vs the amr_cpat_mar-padded patch box;
        ! scatter: rank interior vs the region box), rank-ascending so iterating a list gives the same MPI send/recv order as a
        ! 0..num_procs-1 scan. The clamped interior-frame coord range reproduces both frames (see s_amr_coord_range). Bounded
        ! first dim = max overlap over all blocks (dealloc-realloc each build, like amr_seam_pairs above), not num_procs: a block
        ! spans O(1) ranks. Every consumer runs behind a build_seam_pairs guard, so the arrays are always sized before they are
        ! read.
        mx = 1
        do k = 1, amr_num_blocks
            plo = 0; phi = 0; rlo = 0; rhi = 0
            plo(1) = amr_region_lo_all(1, k) - amr_cpat_mar; phi(1) = amr_region_hi_all(1, k) + amr_cpat_mar
            rlo(1) = amr_region_lo_all(1, k); rhi(1) = amr_region_hi_all(1, k)
            if (n_glb > 0) then
                plo(2) = amr_region_lo_all(2, k) - amr_cpat_mar; phi(2) = amr_region_hi_all(2, k) + amr_cpat_mar
                rlo(2) = amr_region_lo_all(2, k); rhi(2) = amr_region_hi_all(2, k)
            end if
            if (p_glb > 0) then
                plo(3) = amr_region_lo_all(3, k) - amr_cpat_mar; phi(3) = amr_region_hi_all(3, k) + amr_cpat_mar
                rlo(3) = amr_region_lo_all(3, k); rhi(3) = amr_region_hi_all(3, k)
            end if
            mx = max(mx, f_amr_overlap_count(plo, phi), f_amr_overlap_count(rlo, rhi))
        end do
        if (allocated(amr_ovl_gather)) deallocate (amr_ovl_gather)
        if (allocated(amr_ovl_scatter)) deallocate (amr_ovl_scatter)
        allocate (amr_ovl_gather(mx, amr_max_blocks), amr_ovl_scatter(mx, amr_max_blocks))
        do k = 1, amr_num_blocks
            plo = 0; phi = 0; rlo = 0; rhi = 0
            plo(1) = amr_region_lo_all(1, k) - amr_cpat_mar; phi(1) = amr_region_hi_all(1, k) + amr_cpat_mar
            rlo(1) = amr_region_lo_all(1, k); rhi(1) = amr_region_hi_all(1, k)
            if (n_glb > 0) then
                plo(2) = amr_region_lo_all(2, k) - amr_cpat_mar; phi(2) = amr_region_hi_all(2, k) + amr_cpat_mar
                rlo(2) = amr_region_lo_all(2, k); rhi(2) = amr_region_hi_all(2, k)
            end if
            if (p_glb > 0) then
                plo(3) = amr_region_lo_all(3, k) - amr_cpat_mar; phi(3) = amr_region_hi_all(3, k) + amr_cpat_mar
                rlo(3) = amr_region_lo_all(3, k); rhi(3) = amr_region_hi_all(3, k)
            end if
            call s_amr_ranks_overlapping(plo, phi, amr_ovl_gather(:,k), amr_ovl_gather_n(k))
            call s_amr_ranks_overlapping(rlo, rhi, amr_ovl_scatter(:,k), amr_ovl_scatter_n(k))
        end do
        amr_seam_pairs_nblk = amr_num_blocks
        amr_seam_pairs_dirty = .false.

    end subroutine s_amr_build_seam_pairs

    !> Block-to-block fine-fine halo (max_grid_size tiling): overwrite each sub-block's seam ghost cells (faces shared with an
    !! adjacent sub-block) with the neighbour's stage-entry fine interior, so the shared fine flux matches on both sides
    !! (coarse-prolonged seam ghosts would be non-conservative). For each seam pair (xb below, yb above, dim d) the two owners
    !! exchange the buff_size-deep near-seam interior (a wave message per peer, or a local copy when one rank owns both). Buffer is
    !! wp, cast to stp on unpack (identity for stp fields). No-op with a single block / no adjacent pairs (any untiled case, any
    !! np).
    impure subroutine s_amr_fine_fine_post()

        integer :: xb, yb, d, rX, rY, cnt, xm(3), tsz, ierr, fmul, idx
        integer :: ip, boff, tq, sq, qbase, r, sblk, sdlo, sdhi, ublk, udlo, udhi, eblk, edlo, edhi

        amr_sw_nreq = 0; amr_sw_nsame = 0
        if (.not. amr .and. l0_ntile == 0) return
        if (amr_num_blocks < 2) return

        ! iterate the cached same-level seam list (rebuilt only when the topology changes); the list has a fixed (xb, yb) order
        ! so
        ! the paired transfers match.
        if (amr_seam_pairs_dirty .or. amr_seam_pairs_nblk /= amr_num_blocks) call s_amr_build_seam_pairs()
        ! device<->host of the fine state is done per-seam inside s_amr_fine_slice, moving only the buff_size-deep near-seam slab
        ! each
        ! pack/unpack touches (not the whole block); this runs per stage
        if (allocated(amr_sw_plx)) deallocate (amr_sw_plx, amr_sw_ply, amr_sw_pd, amr_sw_pxhi, amr_sw_pfm)
        allocate (amr_sw_plx(amr_num_seam_pairs), amr_sw_ply(amr_num_seam_pairs), amr_sw_pd(amr_num_seam_pairs), &
                  & amr_sw_pxhi(amr_num_seam_pairs), amr_sw_pfm(3, amr_num_seam_pairs))
        ! Seam wave (plan-based exchange): the cross-rank pairs are one aggregated message per (peer, direction): all recvs
        ! posted, all packs into per-transfer pool slices, all sends, one WAITALL, then all unpacks. Every pair contributes
        ! one send and one recv transfer on each of its two owners; both ranks walk the same replicated pair list ascending
        ! with per-peer running offsets, so the wire layout agrees with no metadata exchange. Wire is wp with the same stp
        ! cast on unpack; under MFC_DEBUG each slab carries the identity header [site, sending slot, (d, dlo, dhi),
        ! (cnt, 0, 0)]. Reuses the fill waves' rank-indexed scratch (the waves never overlap in time). Same-rank pairs use
        ! the batched kernel.
        if (.not. allocated(amr_fw_map)) then
            allocate (amr_fw_map(0:num_procs - 1), amr_fw_nx(0:num_procs - 1), amr_fw_pq(0:num_procs - 1), &
                      & amr_fw_pp(0:num_procs - 1))
            amr_fw_map = 0; amr_fw_nx = 0; amr_fw_pq = 0; amr_fw_pp = 0
        end if
        call s_amr_m1_wave_open(5)
        amr_sw_nsame = 0
        amr_sw_snx = 0; amr_sw_snp = 0
        do idx = 1, amr_num_seam_pairs
            xb = amr_seam_pairs(1, idx); yb = amr_seam_pairs(2, idx); d = amr_seam_pairs(3, idx)
            rX = amr_block_owner(xb); rY = amr_block_owner(yb)
            if (proc_rank /= rX .and. proc_rank /= rY) cycle
            ! fine extents from the replicated region metadata (not amr_slots%m/n/p: at np>1 this rank may own only one of the
            ! pair, and the transverse size (used for the buffer count) must be valid for both). A level-L block's region is in
            ! L0-coarse cells but its own grid is rr**L finer, so fine = rr**L*(coarse extent)-1; xb, yb share the level
            ! (same-level seam). Using the L1 factor for an L2 tile would mislocate the seam slice to half the block and fill the
            ! seam ghost from the wrong cells.
            fmul = amr_ref_ratio**amr_block_level(xb)
            xm(1) = fmul*(amr_region_hi_all(1, xb) - amr_region_lo_all(1, xb) + 1) - 1
            xm(2) = merge(fmul*(amr_region_hi_all(2, xb) - amr_region_lo_all(2, xb) + 1) - 1, 0, n_glb > 0)
            xm(3) = merge(fmul*(amr_region_hi_all(3, xb) - amr_region_lo_all(3, xb) + 1) - 1, 0, p_glb > 0)
            ! transverse fine size (dims /= d); xb and yb share it (exact-match seam)
            tsz = 1
            if (d /= 1) tsz = tsz*(xm(1) + 1)
            if (d /= 2 .and. n_glb > 0) tsz = tsz*(xm(2) + 1)
            if (d /= 3 .and. p_glb > 0) tsz = tsz*(xm(3) + 1)
            cnt = sys_size*buff_size*tsz
            if (rX == rY) then  ! same rank owns both: defer to the one batched kernel below (no host buffer, no per-pair launch)
                amr_sw_nsame = amr_sw_nsame + 1
                amr_sw_plx(amr_sw_nsame) = amr_loc_of(xb); amr_sw_ply(amr_sw_nsame) = amr_loc_of(yb)
                amr_sw_pd(amr_sw_nsame) = d; amr_sw_pxhi(amr_sw_nsame) = xm(d); amr_sw_pfm(:,amr_sw_nsame) = xm
                cycle
            end if
            ! cross-rank: append this side's send transfer (the matching recv is appended in the second pair walk below)
            if (proc_rank == rX) then
                r = rY; sblk = xb; sdlo = xm(d) - buff_size + 1; sdhi = xm(d)
            else
                r = rX; sblk = yb; sdlo = 0; sdhi = buff_size - 1
            end if
            if (amr_fw_map(r) == 0) then
                amr_sw_snp = amr_sw_snp + 1
                call s_amr_fw_szi(amr_sw_sprank, amr_sw_snp); call s_amr_fw_szi(amr_sw_sqsz, amr_sw_snp)
                call s_amr_fw_szi(amr_sw_snxp, amr_sw_snp); call s_amr_fw_szi(amr_sw_sqbase, amr_sw_snp)
                amr_fw_map(r) = amr_sw_snp
                amr_sw_sprank(amr_sw_snp) = r
            end if
            amr_sw_snx = amr_sw_snx + 1
            call s_amr_fw_szi(amr_sw_sblk, amr_sw_snx); call s_amr_fw_szi3(amr_sw_sbl, amr_sw_snx)
            call s_amr_fw_szi(amr_sw_spi, amr_sw_snx); call s_amr_fw_szi(amr_sw_sqo, amr_sw_snx)
            call s_amr_fw_szi(amr_sw_spo, amr_sw_snx)
            amr_sw_sblk(amr_sw_snx) = sblk
            amr_sw_sbl(1, amr_sw_snx) = d; amr_sw_sbl(2, amr_sw_snx) = sdlo; amr_sw_sbl(3, amr_sw_snx) = sdhi
            amr_sw_spo(amr_sw_snx) = cnt
            amr_sw_spi(amr_sw_snx) = amr_fw_map(r)
            amr_sw_sqo(amr_sw_snx) = amr_fw_pq(r) + amr_fw_nx(r)*XA_NH
            amr_fw_pq(r) = amr_fw_pq(r) + cnt
            amr_fw_nx(r) = amr_fw_nx(r) + 1
        end do
        qbase = 0
        do ip = 1, amr_sw_snp
            r = amr_sw_sprank(ip)
            amr_sw_snxp(ip) = amr_fw_nx(r)
            amr_sw_sqsz(ip) = amr_fw_pq(r) + amr_fw_nx(r)*XA_NH
            amr_sw_sqbase(ip) = qbase; qbase = qbase + amr_sw_sqsz(ip)
            amr_fw_map(r) = 0; amr_fw_nx(r) = 0; amr_fw_pq(r) = 0
        end do
        call s_amr_fw_szr(amr_sw_sq, qbase)
        ! recv transfers, second walk over the same pairs: unpack destination is my block's ghost slab; the expected header
        ! is the peer's pack (its block + its interior slab bounds), derived from the same replicated metadata
        amr_sw_rnx = 0; amr_sw_rnp = 0
        do idx = 1, amr_num_seam_pairs
            xb = amr_seam_pairs(1, idx); yb = amr_seam_pairs(2, idx); d = amr_seam_pairs(3, idx)
            rX = amr_block_owner(xb); rY = amr_block_owner(yb)
            if (rX == rY) cycle
            if (proc_rank /= rX .and. proc_rank /= rY) cycle
            fmul = amr_ref_ratio**amr_block_level(xb)
            xm(1) = fmul*(amr_region_hi_all(1, xb) - amr_region_lo_all(1, xb) + 1) - 1
            xm(2) = merge(fmul*(amr_region_hi_all(2, xb) - amr_region_lo_all(2, xb) + 1) - 1, 0, n_glb > 0)
            xm(3) = merge(fmul*(amr_region_hi_all(3, xb) - amr_region_lo_all(3, xb) + 1) - 1, 0, p_glb > 0)
            tsz = 1
            if (d /= 1) tsz = tsz*(xm(1) + 1)
            if (d /= 2 .and. n_glb > 0) tsz = tsz*(xm(2) + 1)
            if (d /= 3 .and. p_glb > 0) tsz = tsz*(xm(3) + 1)
            cnt = sys_size*buff_size*tsz
            if (proc_rank == rX) then
                ! yb's low interior arrives -> xb's high ghost
                r = rY; ublk = xb; udlo = xm(d) + 1; udhi = xm(d) + buff_size
                eblk = yb; edlo = 0; edhi = buff_size - 1
            else
                ! xb's high interior arrives -> yb's low ghost
                r = rX; ublk = yb; udlo = -buff_size; udhi = -1
                eblk = xb; edlo = xm(d) - buff_size + 1; edhi = xm(d)
            end if
            if (amr_fw_map(r) == 0) then
                amr_sw_rnp = amr_sw_rnp + 1
                call s_amr_fw_szi(amr_sw_rprank, amr_sw_rnp); call s_amr_fw_szi(amr_sw_rqsz, amr_sw_rnp)
                call s_amr_fw_szi(amr_sw_rnxp, amr_sw_rnp); call s_amr_fw_szi(amr_sw_rqbase, amr_sw_rnp)
                amr_fw_map(r) = amr_sw_rnp
                amr_sw_rprank(amr_sw_rnp) = r
            end if
            amr_sw_rnx = amr_sw_rnx + 1
            call s_amr_fw_szi(amr_sw_rblk, amr_sw_rnx); call s_amr_fw_szi3(amr_sw_rbl, amr_sw_rnx)
            call s_amr_fw_szi3(amr_sw_rbh, amr_sw_rnx); call s_amr_fw_szi(amr_sw_rpi, amr_sw_rnx)
            call s_amr_fw_szi(amr_sw_rqo, amr_sw_rnx); call s_amr_fw_szi(amr_sw_rpo, amr_sw_rnx)
            amr_sw_rblk(amr_sw_rnx) = ublk
            amr_sw_rbl(1, amr_sw_rnx) = d; amr_sw_rbl(2, amr_sw_rnx) = udlo; amr_sw_rbl(3, amr_sw_rnx) = udhi
            amr_sw_rbh(1, amr_sw_rnx) = eblk; amr_sw_rbh(2, amr_sw_rnx) = edlo; amr_sw_rbh(3, amr_sw_rnx) = edhi
            amr_sw_rpo(amr_sw_rnx) = cnt
            amr_sw_rpi(amr_sw_rnx) = amr_fw_map(r)
            amr_sw_rqo(amr_sw_rnx) = amr_fw_pq(r) + amr_fw_nx(r)*XA_NH
            amr_fw_pq(r) = amr_fw_pq(r) + cnt
            amr_fw_nx(r) = amr_fw_nx(r) + 1
        end do
        qbase = 0
        do ip = 1, amr_sw_rnp
            r = amr_sw_rprank(ip)
            amr_sw_rnxp(ip) = amr_fw_nx(r)
            amr_sw_rqsz(ip) = amr_fw_pq(r) + amr_fw_nx(r)*XA_NH
            amr_sw_rqbase(ip) = qbase; qbase = qbase + amr_sw_rqsz(ip)
            amr_fw_map(r) = 0; amr_fw_nx(r) = 0; amr_fw_pq(r) = 0
        end do
        call s_amr_fw_szr(amr_sw_rq, qbase)
        amr_sw_nreq = amr_sw_snp + amr_sw_rnp
        call s_amr_fw_szi(amr_sw_req, amr_sw_nreq); call s_amr_fw_szi(amr_sw_reqw, amr_sw_nreq)

        amr_sw_nreq = 0
#ifdef MFC_MPI
        do ip = 1, amr_sw_rnp
            sq = f_amr_m1_seq(amr_sw_rprank(ip), 2); tq = f_amr_m1_tag(5, sq)
            call s_xa_rec(XA_F6W_RCV, 2, amr_sw_rqsz(ip) - amr_sw_rnxp(ip)*XA_NH, tq, peer=amr_sw_rprank(ip), &
                          & key=amr_sw_rnxp(ip), seq=sq)
            amr_sw_nreq = amr_sw_nreq + 1; amr_sw_reqw(amr_sw_nreq) = amr_sw_rqsz(ip)
            call MPI_IRECV(amr_sw_rq(amr_sw_rqbase(ip) + 1), amr_sw_rqsz(ip), mpi_p, amr_sw_rprank(ip), tq, MPI_COMM_WORLD, &
                           & amr_sw_req(amr_sw_nreq), ierr)
        end do
#endif
        do idx = 1, amr_sw_snx
            cnt = amr_sw_spo(idx)
            boff = amr_sw_sqbase(amr_sw_spi(idx)) + amr_sw_sqo(idx)
            call s_amr_fine_slice(amr_sw_sblk(idx), amr_sw_sbl(1, idx), amr_sw_sbl(2, idx), amr_sw_sbl(3, idx), &
                                  & amr_sw_sq(boff + XA_NH + 1:boff + XA_NH + cnt), 1)
            if (XA_NH > 0) call s_xa_hdr_pack(amr_sw_sq(boff + 1:boff + XA_NH), XA_F6W_SND, amr_sw_sblk(idx), amr_sw_sbl(:,idx), &
                & [cnt, 0, 0])
        end do
#ifdef MFC_MPI
        do ip = 1, amr_sw_snp
            sq = f_amr_m1_seq(amr_sw_sprank(ip), 1); tq = f_amr_m1_tag(5, sq)
            call s_xa_rec(XA_F6W_SND, 1, amr_sw_sqsz(ip) - amr_sw_snxp(ip)*XA_NH, tq, peer=amr_sw_sprank(ip), &
                          & key=amr_sw_snxp(ip), seq=sq)
            amr_sw_nreq = amr_sw_nreq + 1; amr_sw_reqw(amr_sw_nreq) = -1
            call MPI_ISEND(amr_sw_sq(amr_sw_sqbase(ip) + 1), amr_sw_sqsz(ip), mpi_p, amr_sw_sprank(ip), tq, MPI_COMM_WORLD, &
                           & amr_sw_req(amr_sw_nreq), ierr)
        end do
#endif

    end subroutine s_amr_fine_fine_post

    !> Drain the seam wave posted by s_amr_fine_fine_post: wait, unpack the cross-rank ghosts, run the same-rank pairs. Its ghost
    !! writes stay after the coarse and parent fills (the seam wins on faces, the coarse fill on edges/corners).
    impure subroutine s_amr_fine_fine_drain()

        integer :: idx, cnt, boff, ierr

        if (amr_sw_nreq == 0 .and. amr_sw_nsame == 0) return
#ifdef MFC_MPI
        if (amr_sw_nreq > 0) then
#ifdef MFC_DEBUG
            block
                integer :: st(MPI_STATUS_SIZE, amr_sw_nreq), gotw, q
                call s_wait_tic()
                call MPI_WAITALL(amr_sw_nreq, amr_sw_req, st, ierr)
                call s_wait_toc(WT_SEAM)
                do q = 1, amr_sw_nreq
                    if (amr_sw_reqw(q) < 0) cycle
                    call MPI_GET_COUNT(st(:,q), mpi_p, gotw, ierr)
                    @:ASSERT(gotw == amr_sw_reqw(q), "seam wave: received message length differs from the plan")
                end do
            end block
#else
            call s_wait_tic()
            call MPI_WAITALL(amr_sw_nreq, amr_sw_req, MPI_STATUSES_IGNORE, ierr)
            call s_wait_toc(WT_SEAM)
#endif
        end if
#endif
        do idx = 1, amr_sw_rnx
            cnt = amr_sw_rpo(idx)
            boff = amr_sw_rqbase(amr_sw_rpi(idx)) + amr_sw_rqo(idx)
            if (XA_NH > 0) call s_xa_hdr_check(amr_sw_rq(boff + 1:boff + XA_NH), XA_F6W_SND, amr_sw_rbh(1, idx), [amr_sw_rbl(1, &
                & idx), amr_sw_rbh(2, idx), amr_sw_rbh(3, idx)], [cnt, 0, 0])
            call s_amr_fine_slice(amr_sw_rblk(idx), amr_sw_rbl(1, idx), amr_sw_rbl(2, idx), amr_sw_rbl(3, idx), &
                                  & amr_sw_rq(boff + XA_NH + 1:boff + XA_NH + cnt), -1)
        end do
        ! Every same-rank pair in one launch. Two things make this safe. Fusing across pairs: the four slabs of a pair are
        ! disjoint and no pair writes another's source. Deferring past the MPI pairs above: every seam operation reads only
        ! interior cells and writes only ghost cells (the packs read [xhi-buff+1:xhi] / [0:buff-1] and the unpacks write
        ! [xhi+1:xhi+buff] / [-buff:-1]), so no seam operation can observe another's write, in either path.
        if (amr_sw_nsame > 0) call s_amr_fine_seam_exchange(amr_sw_nsame, amr_sw_plx(1:amr_sw_nsame), amr_sw_ply(1:amr_sw_nsame), &
            & amr_sw_pd(1:amr_sw_nsame), amr_sw_pxhi(1:amr_sw_nsame), amr_sw_pfm(:,1:amr_sw_nsame), buff_size)
        call s_amr_select_slot(1)

    end subroutine s_amr_fine_fine_drain

    !> The seam exchange as one call (post + drain): the early-post-off path.
    impure subroutine s_amr_fine_fine_halo()

        call s_amr_fine_fine_post()
        call s_amr_fine_fine_drain()

    end subroutine s_amr_fine_fine_halo

    !> Ring clip: decompose the hollow shell of patch [plo:phi] minus its open core [clo:chi] (same integer frame; collapsed dims
    !! pass plo=phi=clo=chi=0 so they never cut the core) into at most 6 disjoint slabs in a fixed order: x-low/x-high spanning the
    !! full transverse extent, then y-low/high restricted to the core's x-interval, then z-low/high restricted in x and y. Width<=2
    !! cores are empty (shell = whole patch, legal); the width-1 double-cover is resolved by clamping the high slab past the low
    !! one. Both sides of every clipped exchange derive the same list from replicated metadata, so the wire layout needs no
    !! handshake (see misc/amr_ledger/stepfill_ring_clip.md).
    impure subroutine s_amr_shell_slabs(plo, phi, clo, chi, ns, sb1, se1, sb2, se2, sb3, se3, cells)

        integer, intent(in)  :: plo(3), phi(3), clo(3), chi(3)
        integer, intent(out) :: ns, sb1(6), se1(6), sb2(6), se2(6), sb3(6), se3(6), cells
        integer              :: cb(3, 6), ce(3, 6), s, ss
        integer(8)           :: words, patchw, corew

        cb(:,1) = plo; ce(:,1) = [clo(1) - 1, phi(2), phi(3)]
        cb(:,2) = [max(chi(1) + 1, clo(1)), plo(2), plo(3)]; ce(:,2) = phi
        cb(:,3) = [clo(1), plo(2), plo(3)]; ce(:,3) = [chi(1), clo(2) - 1, phi(3)]
        cb(:,4) = [clo(1), max(chi(2) + 1, clo(2)), plo(3)]; ce(:,4) = [chi(1), phi(2), phi(3)]
        cb(:,5) = [clo(1), clo(2), plo(3)]; ce(:,5) = [chi(1), chi(2), clo(3) - 1]
        cb(:,6) = [clo(1), clo(2), max(chi(3) + 1, clo(3))]; ce(:,6) = [chi(1), chi(2), phi(3)]
        ns = 0; words = 0
        do s = 1, 6
            if (cb(1, s) > ce(1, s) .or. cb(2, s) > ce(2, s) .or. cb(3, s) > ce(3, s)) cycle
            ns = ns + 1
            sb1(ns) = cb(1, s); se1(ns) = ce(1, s)
            sb2(ns) = cb(2, s); se2(ns) = ce(2, s)
            sb3(ns) = cb(3, s); se3(ns) = ce(3, s)
            words = words + int(se1(ns) - sb1(ns) + 1, 8)*int(se2(ns) - sb2(ns) + 1, 8)*int(se3(ns) - sb3(ns) + 1, 8)
        end do
        ! the slabs must tile the shell exactly: pairwise disjoint, cells summing to patch - core
        do s = 1, ns - 1
            do ss = s + 1, ns
                @:ASSERT(max(sb1(s), sb1(ss)) > min(se1(s), se1(ss)) .or. max(sb2(s), sb2(ss)) > min(se2(s), &
                         & se2(ss)) .or. max(sb3(s), sb3(ss)) > min(se3(s), se3(ss)), "shell slabs: overlap")
            end do
        end do
        patchw = int(phi(1) - plo(1) + 1, 8)*int(phi(2) - plo(2) + 1, 8)*int(phi(3) - plo(3) + 1, 8)
        corew = int(max(chi(1) - clo(1) + 1, 0), 8)*int(max(chi(2) - clo(2) + 1, 0), 8)*int(max(chi(3) - clo(3) + 1, 0), 8)
        @:ASSERT(words == patchw - corew, "shell slabs: coverage mismatch")
        cells = int(words)

    end subroutine s_amr_shell_slabs

    !> Intersect the shell-slab list with box [bl:bh]: the surviving clipped slabs in the same fixed order (each exchange side
    !! derives an identical list from replicated data, so empties drop symmetrically) plus their total cell count.
    impure subroutine s_amr_shell_clip(ns, sb1, se1, sb2, se2, sb3, se3, bl, bh, ms, tb1, te1, tb2, te2, tb3, te3, cells)

        integer, intent(in)  :: ns, sb1(6), se1(6), sb2(6), se2(6), sb3(6), se3(6), bl(3), bh(3)
        integer, intent(out) :: ms, tb1(6), te1(6), tb2(6), te2(6), tb3(6), te3(6), cells
        integer              :: s, l1, u1, l2, u2, l3, u3

        ms = 0; cells = 0
        do s = 1, ns
            l1 = max(sb1(s), bl(1)); u1 = min(se1(s), bh(1))
            l2 = max(sb2(s), bl(2)); u2 = min(se2(s), bh(2))
            l3 = max(sb3(s), bl(3)); u3 = min(se3(s), bh(3))
            if (l1 > u1 .or. l2 > u2 .or. l3 > u3) cycle
            ms = ms + 1
            tb1(ms) = l1; te1(ms) = u1; tb2(ms) = l2; te2(ms) = u2; tb3(ms) = l3; te3(ms) = u3
            cells = cells + (u1 - l1 + 1)*(u2 - l2 + 1)*(u3 - l3 + 1)
        end do

    end subroutine s_amr_shell_clip

    !> Validation arm (debug builds only): flood the current patch extent of amr_cg with quiet NaN before a clipped gather writes
    !! its shell, so a consumer read of any unshipped cell (core or a missed shell slab) NaNs the ghost fill within a step.
    impure subroutine s_amr_poison_patch_device(w1, w2, w3)

        use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        integer, intent(in) :: w1, w2, w3
        integer             :: i, g1, g2, g3
        real(stp)           :: nanv

        nanv = real(ieee_value(0._wp, ieee_quiet_nan), stp)
        $:GPU_PARALLEL_LOOP(collapse=4)
        do i = 1, sys_size
            do g3 = 0, w3
                do g2 = 0, w2
                    do g1 = 0, w1
                        amr_cg(i)%sf(g1, g2, g3) = nanv
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_poison_patch_device

    !> Ring-clipped runtime own-box copy (device): the owner's shell-slab / own-box intersections [tb:te] global from q_coarse into
    !! amr_cg in the patch-local frame, with no host round-trip and one fused kernel over the slab concatenation (the ghost-fill
    !! kernel's flat-index idiom). Same index map and direct stp assignment as the host path in s_amr_gather_coarse_patch.
    impure subroutine s_amr_gather_own_shell_device(q_coarse, ms, tb1, te1, tb2, te2, tb3, te3, o1, o2, o3)

        type(scalar_field), dimension(sys_size), intent(in) :: q_coarse
        integer, intent(in)                                 :: ms, tb1(6), te1(6), tb2(6), te2(6), tb3(6), te3(6), o1, o2, o3
        integer                                             :: lb1(6), le1(6), lb2(6), le2(6), lb3(6), le3(6), soff(6), scnt(6)
        integer                                             :: i, s, ss, g, r, n1, n2, g1, g2, g3, stot, coff1, coff2, coff3

        ! scalar/local copies: no host array may be referenced inside the device region (nvfortran/Cray demand it present)

        coff1 = amr_cpat_off(1); coff2 = amr_cpat_off(2); coff3 = amr_cpat_off(3)
        soff(1) = 0
        do s = 1, ms
            lb1(s) = tb1(s); le1(s) = te1(s); lb2(s) = tb2(s); le2(s) = te2(s); lb3(s) = tb3(s); le3(s) = te3(s)
            scnt(s) = (te1(s) - tb1(s) + 1)*(te2(s) - tb2(s) + 1)*(te3(s) - tb3(s) + 1)
            if (s < ms) soff(s + 1) = soff(s) + scnt(s)
        end do
        stot = soff(ms) + scnt(ms)
        amr_slab_tab(1,:) = lb1; amr_slab_tab(2,:) = le1; amr_slab_tab(3,:) = lb2; amr_slab_tab(4,:) = le2
        amr_slab_tab(5,:) = lb3; amr_slab_tab(6,:) = le3; amr_slab_tab(7,:) = soff; amr_slab_tab(8,:) = scnt
        $:GPU_UPDATE(device='[amr_slab_tab]')
        $:GPU_PARALLEL_LOOP(collapse=2, private='[s, ss, r, n1, n2, g1, g2, g3]')
        do i = 1, sys_size
            do g = 0, stot - 1
                s = 1  ! decode the flat index: ms <= 6, so a scan beats storing a per-cell slab map
                do ss = 2, ms
                    if (g >= amr_slab_tab(7, ss)) s = ss
                end do
                r = g - amr_slab_tab(7, s)
                n1 = amr_slab_tab(2, s) - amr_slab_tab(1, s) + 1; n2 = amr_slab_tab(4, s) - amr_slab_tab(3, s) + 1
                g1 = amr_slab_tab(1, s) + mod(r, n1)
                g2 = amr_slab_tab(3, s) + mod(r/n1, n2)
                g3 = amr_slab_tab(5, s) + r/(n1*n2)
                amr_cg(i)%sf(g1 - coff1, g2 - coff2, g3 - coff3) = q_coarse(i)%sf(g1 - o1, g2 - o2, g3 - o3)
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_gather_own_shell_device

    !> Flat-index decode shared by the three fused exchange kernels (amr_device_pack): binary-search the plan prefix for the
    !! transfer t that owns flat element e, then invert its linear buffer index into (i, g3, g2, g1). One text, so the three kernels
    !! cannot drift from each other or from the per-transfer index map.
    #:def FX_DECODE()
        lo = t0; hi = t1
        do while (lo < hi)
            md = (lo + hi)/2
            if (pre(md + 1) < e) then; lo = md + 1; else; hi = md; end if
        end do
        t = lo
        r = e - pre(t) - 1
        n1 = pl(4, t); n2 = pl(5, t); n3 = pl(6, t)
        g1 = pl(1, t) + mod(r, n1); g2 = pl(2, t) + mod(r/n1, n2); g3 = pl(3, t) + mod(r/(n1*n2), n3)
        i = r/(n1*n2*n3) + 1
    #:enddef

    !> Fused-pack plan (amr_device_pack) for transfers 1:nt of a wave's send or recv table: the slab corner, its extents, its
    !! absolute payload offset in the wire pool, and the exclusive element prefix amr_fx_pre (amr_fx_pre(t+1) - amr_fx_pre(t) is
    !! transfer t's word count). Rows 8:11 are left to the F2 pack site, the only caller whose source varies per transfer.
    impure subroutine s_amr_fx_plan(bl, bh, pbase, xoff, pi, nt)

        integer, intent(in) :: bl(:,:), bh(:,:), pbase(:), xoff(:), pi(:), nt
        integer             :: t, e

        if (allocated(amr_fx_pl)) then
            if (size(amr_fx_pl, 2) < nt) deallocate (amr_fx_pl, amr_fx_pre)
        end if
        if (.not. allocated(amr_fx_pl)) allocate (amr_fx_pl(11, max(nt, 64)), amr_fx_pre(max(nt, 64) + 1))
        e = 0
        do t = 1, nt
            amr_fx_pl(1:3,t) = bl(:,t)
            amr_fx_pl(4:6,t) = bh(:,t) - bl(:,t) + 1
            amr_fx_pl(7, t) = pbase(pi(t)) + xoff(t) + XA_NH
            amr_fx_pre(t) = e
            e = e + sys_size*amr_fx_pl(4, t)*amr_fx_pl(5, t)*amr_fx_pl(6, t)
        end do
        amr_fx_pre(nt + 1) = e

    end subroutine s_amr_fx_plan

    !> Longest run t0:t1 of one box's recv transfers that is contiguous in the wire pool, so the run unpacks in one launch from one
    !! pool slice. Both wave plans lay a box's transfers from one peer back to back, so a box whose contributors span several peers
    !! fuses per peer and a run of one is a single-transfer launch.
    pure subroutine s_amr_fx_run(k, blk, nx, t0, t1)

        integer, intent(in)  :: k, blk(:), nx, t0
        integer, intent(out) :: t1

        t1 = t0
        do while (t1 < nx)
            if (blk(t1 + 1) /= k) exit
            if (amr_fx_pl(7, t1 + 1) - XA_NH /= amr_fx_pl(7, t1) + amr_fx_pre(t1 + 1) - amr_fx_pre(t1)) exit
            t1 = t1 + 1
        end do

    end subroutine s_amr_fx_run

    !> Fused F1 pack (amr_device_pack): the wave's whole send list in one launch. The flat element index is binary-searched against
    !! the plan prefix for its transfer, then decoded to (i, g3, g2, g1), the exact inverse of the linear buffer index
    !! s_amr_pack_box_device writes. Same source expression, same wp cast, same destination word, so the wire bytes are identical to
    !! the per-transfer packs.
    impure subroutine s_amr_fx_pack_box(q_coarse, t0, t1, o1, o2, o3, pl, pre, buf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_coarse
        integer, intent(in)                                 :: t0, t1, o1, o2, o3
        integer, intent(in), contiguous                     :: pl(:,:), pre(:)
        real(wp), intent(inout), contiguous                 :: buf(:)
        integer                                             :: e, t, lo, hi, md, r, n1, n2, n3, i, g1, g2, g3

        $:GPU_PARALLEL_LOOP(copyin='[pl, pre]', copyout='[buf]', private='[t, lo, hi, md, r, n1, n2, n3, i, g1, g2, g3]')
        do e = pre(t0) + 1, pre(t1 + 1)
            @:FX_DECODE()
            buf(pl(7, t) + 1 + r) = real(q_coarse(i)%sf(g1 - o1, g2 - o2, g3 - o3), wp)
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_fx_pack_box

    !> Fused F2 pack (amr_device_pack): as s_amr_fx_pack_box, but the source is the flat store at the per-transfer parent slot
    !! pl(8,:) in the per-transfer child patch frame pl(9:11,:): the s_amr_pack_parent_box_device_cons body, one launch for the
    !! whole send list.
    impure subroutine s_amr_fx_pack_parent(t0, t1, pl, pre, buf)

        integer, intent(in)                 :: t0, t1
        integer, intent(in), contiguous     :: pl(:,:), pre(:)
        real(wp), intent(inout), contiguous :: buf(:)
        integer                             :: e, t, lo, hi, md, r, n1, n2, n3, i, g1, g2, g3

        $:GPU_PARALLEL_LOOP(copyin='[pl, pre]', copyout='[buf]', private='[t, lo, hi, md, r, n1, n2, n3, i, g1, g2, g3]')
        do e = pre(t0) + 1, pre(t1 + 1)
            @:FX_DECODE()
            buf(pl(7, t) + 1 + r) = real(amr_cons_st(g1 + pl(9, t), g2 + pl(10, t), g3 + pl(11, t), i, pl(8, t)), wp)
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_fx_pack_parent

    !> Fused unpack (amr_device_pack) of the transfer run t0:t1 into the single amr_cg patch: the run is one box's transfers that
    !! are contiguous in the wire pool (the caller checks that), so buf is the pool slice starting at word base+1 and the whole run
    !! is one launch. Serves F1 (bounds global, frame c1:c3 = amr_cpat_off) and F2 (bounds patch-local, frame 0) alike; the stp cast
    !! is the assignment both per-box unpacks perform. Cross-box fusion would need a per-box patch store.
    impure subroutine s_amr_fx_unpack(t0, t1, base, c1, c2, c3, pl, pre, buf)

        integer, intent(in)              :: t0, t1, base, c1, c2, c3
        integer, intent(in), contiguous  :: pl(:,:), pre(:)
        real(wp), intent(in), contiguous :: buf(:)
        integer                          :: e, t, lo, hi, md, r, n1, n2, n3, i, g1, g2, g3

        $:GPU_PARALLEL_LOOP(copyin='[pl, pre, buf]', private='[t, lo, hi, md, r, n1, n2, n3, i, g1, g2, g3]')
        do e = pre(t0) + 1, pre(t1 + 1)
            @:FX_DECODE()
            amr_cg(i)%sf(g1 - c1, g2 - c2, g3 - c3) = real(buf(pl(7, t) - base + 1 + r), stp)
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_fx_unpack

    !> Per-stage level-1 fill as one exchange wave: derive this stage's full (box, contributor) transfer set from the replicated
    !! caches, exchange one aggregated message per (peer, family) (F1 q_cons and, under non-polytropic QBMM, the F3 pb/mv twin) with
    !! all recvs posted first, then packs, then sends, then one waitall, and finally consume owned boxes in ascending slot order
    !! through the single amr_cg patch (own-box device copy + per-slab device unpack + ghost fill). Level>=2 blocks use the
    !! parent-fill wave. Under MFC_DEBUG every slab carries the identity header, verified at consume, and each received message
    !! length is checked against the plan.
    impure subroutine s_amr_stage_fill_wave(q_cons_coarse, pb_in, mv_in)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_coarse
        real(stp), dimension(:,:,:,:,:), intent(inout) :: pb_in, mv_in
        logical :: do_pbmv
        integer :: k, r, idx, ix, ip, owner, o1, o2, o3, qsz, psz, cellsz, tq, tp, sq, nreq, qbase, pbase, ierr, kk, kk2
        integer :: v1hi, v2hi, v3hi, plo(3), phi(3), crlo(3), crhi(3), bl(3), bh(3), boff, sqtot, ie, jx
        logical :: fuse  !< amr_device_pack: fused per-family packs (the pbmv twin keeps its per-box wire contract)
        integer :: clo(3), chi(3), nsh, msl, isl, scells
        integer :: shb1(6), she1(6), shb2(6), she2(6), shb3(6), she3(6), tb1(6), te1(6), tb2(6), te2(6), tb3(6), te3(6)

        if (amr_num_blocks <= 0) return
        @:ASSERT(amr_gsnd_n == 0, "stage-fill wave: the deferred gather-send pool must be drained")

        do_pbmv = qbmm .and. .not. polytropic
        ! the F3 pb/mv twin keeps its full-box per-transfer wire contract and a different word count per cell, so the fused
        ! plan (sys_size words per cell) covers the q_cons families only
        fuse = amr_device_pack .and. .not. do_pbmv
        cellsz = 0
        if (do_pbmv) cellsz = 2*nnode*nb
        o1 = start_idx(1); o2 = 0; o3 = 0
        if (n_glb > 0) o2 = start_idx(2)
        if (p_glb > 0) o3 = start_idx(3)
        ! two bands, one wave: the second open re-clears the shared per-peer seq counters before anything has posted, so a
        ! peer's q message takes seq 1 and its pb/mv message seq 2, the same relative order on both ends
        call s_amr_m1_wave_open(3); call s_amr_m1_wave_open(4)

        call s_phase_tic(PH_GATHER)
        call s_phase_tic(PH_GWPLAN)
        ! block set changed: rebuild the cached overlap-rank lists before reading them (same lazy trigger as the per-box path)
        if (amr_seam_pairs_dirty .or. amr_seam_pairs_nblk /= amr_num_blocks) call s_amr_build_seam_pairs()
        if (.not. allocated(amr_fw_map)) then
            allocate (amr_fw_map(0:num_procs - 1), amr_fw_nx(0:num_procs - 1), amr_fw_pq(0:num_procs - 1), &
                      & amr_fw_pp(0:num_procs - 1))
            amr_fw_map = 0; amr_fw_nx = 0; amr_fw_pq = 0; amr_fw_pp = 0
        end if

        ! send side: for every level-1 box someone else owns, my coarse-range slice of its padded patch box.
        amr_fw_snx = 0; amr_fw_snp = 0
        ! the Lagrangian-overlap safety check must cover every level-1 block (mine included), so it keeps a dedicated
        ! gated scan over all blocks rather than the owned/contributor lists
        if (bubbles_lagrange) then
            do k = 1, amr_num_blocks
                if (amr_block_level(k) /= 1) cycle
                call s_amr_select_slot(k)
                call s_amr_check_lag_clear()
            end do
        end if
        ! walk the padded cached list (region +/- amr_cpat_mar vs my coarse range, this loop's exact predicate); the body
        ! keeps its own intersection + empty cycle as belt-and-braces.
        call s_amr_refresh_lists()
        do kk2 = 1, amr_n_l1p
            k = amr_l1p_blk(kk2)
            call s_amr_select_slot(k)
            owner = amr_block_owner(k)
            plo(1) = amr_region_lo_all(1, k) - amr_cpat_mar; plo(2) = 0; plo(3) = 0
            if (n_glb > 0) plo(2) = amr_region_lo_all(2, k) - amr_cpat_mar
            if (p_glb > 0) plo(3) = amr_region_lo_all(3, k) - amr_cpat_mar
            v1hi = (amr_region_hi_all(1, k) - amr_region_lo_all(1, k)) + 2*amr_cpat_mar
            v2hi = 0; v3hi = 0
            if (n_glb > 0) v2hi = (amr_region_hi_all(2, k) - amr_region_lo_all(2, k)) + 2*amr_cpat_mar
            if (p_glb > 0) v3hi = (amr_region_hi_all(3, k) - amr_region_lo_all(3, k)) + 2*amr_cpat_mar
            phi(1) = plo(1) + v1hi; phi(2) = plo(2) + v2hi; phi(3) = plo(3) + v3hi
            call s_amr_rank_coarse_range(proc_rank, crlo, crhi)
            call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
            if (bl(1) > bh(1) .or. bl(2) > bh(2) .or. bl(3) > bh(3)) cycle
            ! ring clip (runtime q-only path): consumers of amr_cg read only the patch's hollow shell, so ship only the
            ! shell's intersection with this rank's slice, as up to 6 sub-slab transfers, derived identically on both
            ! sides from replicated metadata. The pbmv gather keeps its full-box wire contract, so qbmm+non-polytropic
            ! runs stay unclipped (full slab).
            if (do_pbmv) then
                msl = 1
                tb1(1) = bl(1); te1(1) = bh(1); tb2(1) = bl(2); te2(1) = bh(2); tb3(1) = bl(3); te3(1) = bh(3)
            else
                clo = 0; chi = 0
                clo(1) = amr_region_lo_all(1, k) + 1; chi(1) = amr_region_hi_all(1, k) - 1
                if (n_glb > 0) then; clo(2) = amr_region_lo_all(2, k) + 1; chi(2) = amr_region_hi_all(2, k) - 1; end if
                if (p_glb > 0) then; clo(3) = amr_region_lo_all(3, k) + 1; chi(3) = amr_region_hi_all(3, k) - 1; end if
                call s_amr_shell_slabs(plo, phi, clo, chi, nsh, shb1, she1, shb2, she2, shb3, she3, scells)
                call s_amr_shell_clip(nsh, shb1, she1, shb2, she2, shb3, she3, bl, bh, msl, tb1, te1, tb2, te2, tb3, te3, scells)
                if (msl == 0) cycle
            end if
            do isl = 1, msl
                bl = [tb1(isl), tb2(isl), tb3(isl)]; bh = [te1(isl), te2(isl), te3(isl)]
                qsz = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                psz = cellsz*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                if (amr_fw_map(owner) == 0) then
                    amr_fw_snp = amr_fw_snp + 1
                    call s_amr_fw_szi(amr_fw_sprank, amr_fw_snp); call s_amr_fw_szi(amr_fw_sqsz, amr_fw_snp)
                    call s_amr_fw_szi(amr_fw_spsz, amr_fw_snp); call s_amr_fw_szi(amr_fw_snxp, amr_fw_snp)
                    call s_amr_fw_szi(amr_fw_sqbase, amr_fw_snp); call s_amr_fw_szi(amr_fw_spbase, amr_fw_snp)
                    amr_fw_map(owner) = amr_fw_snp
                    amr_fw_sprank(amr_fw_snp) = owner
                end if
                amr_fw_snx = amr_fw_snx + 1
                call s_amr_fw_szi(amr_fw_sblk, amr_fw_snx); call s_amr_fw_szi3(amr_fw_sbl, amr_fw_snx)
                call s_amr_fw_szi3(amr_fw_sbh, amr_fw_snx); call s_amr_fw_szi(amr_fw_spi, amr_fw_snx)
                call s_amr_fw_szi(amr_fw_sqo, amr_fw_snx); call s_amr_fw_szi(amr_fw_spo, amr_fw_snx)
                amr_fw_sblk(amr_fw_snx) = k; amr_fw_sbl(:,amr_fw_snx) = bl; amr_fw_sbh(:,amr_fw_snx) = bh
                amr_fw_spi(amr_fw_snx) = amr_fw_map(owner)
                amr_fw_sqo(amr_fw_snx) = amr_fw_pq(owner) + amr_fw_nx(owner)*XA_NH
                amr_fw_spo(amr_fw_snx) = amr_fw_pp(owner) + amr_fw_nx(owner)*XA_NH
                amr_fw_pq(owner) = amr_fw_pq(owner) + qsz
                amr_fw_pp(owner) = amr_fw_pp(owner) + psz
                amr_fw_nx(owner) = amr_fw_nx(owner) + 1
            end do
        end do
        qbase = 0; pbase = 0
        do ip = 1, amr_fw_snp
            r = amr_fw_sprank(ip)
            amr_fw_snxp(ip) = amr_fw_nx(r)
            amr_fw_sqsz(ip) = amr_fw_pq(r) + amr_fw_nx(r)*XA_NH
            amr_fw_spsz(ip) = amr_fw_pp(r) + amr_fw_nx(r)*XA_NH
            amr_fw_sqbase(ip) = qbase; qbase = qbase + amr_fw_sqsz(ip)
            amr_fw_spbase(ip) = pbase; pbase = pbase + amr_fw_spsz(ip)
            amr_fw_map(r) = 0; amr_fw_nx(r) = 0; amr_fw_pq(r) = 0; amr_fw_pp(r) = 0
        end do
        call s_amr_fw_szr(amr_fw_sq, qbase, amr_fw_dev)
        if (do_pbmv) call s_amr_fw_szr(amr_fw_sp, pbase, amr_fw_dev)
        sqtot = qbase

        ! recv side: for every level-1 box I own, each listed contributor's slice (owner excluded; the own box is a device
        ! copy at consume). Both sides enumerate boxes ascending with per-rank running offsets, so the offsets agree.
        amr_fw_rnx = 0; amr_fw_rnp = 0
        call s_amr_refresh_my_blocks()
        do kk = 1, amr_n_my
            k = amr_my_blk(kk)
            if (amr_block_level(k) /= 1) cycle
            plo(1) = amr_region_lo_all(1, k) - amr_cpat_mar; plo(2) = 0; plo(3) = 0
            if (n_glb > 0) plo(2) = amr_region_lo_all(2, k) - amr_cpat_mar
            if (p_glb > 0) plo(3) = amr_region_lo_all(3, k) - amr_cpat_mar
            v1hi = (amr_region_hi_all(1, k) - amr_region_lo_all(1, k)) + 2*amr_cpat_mar
            v2hi = 0; v3hi = 0
            if (n_glb > 0) v2hi = (amr_region_hi_all(2, k) - amr_region_lo_all(2, k)) + 2*amr_cpat_mar
            if (p_glb > 0) v3hi = (amr_region_hi_all(3, k) - amr_region_lo_all(3, k)) + 2*amr_cpat_mar
            phi(1) = plo(1) + v1hi; phi(2) = plo(2) + v2hi; phi(3) = plo(3) + v3hi
            ! ring clip: the shell is a per-box property; clip each contributor's slice against it (mirror of the
            ! send walk, so both sides derive the identical sub-slab list)
            if (.not. do_pbmv) then
                clo = 0; chi = 0
                clo(1) = amr_region_lo_all(1, k) + 1; chi(1) = amr_region_hi_all(1, k) - 1
                if (n_glb > 0) then; clo(2) = amr_region_lo_all(2, k) + 1; chi(2) = amr_region_hi_all(2, k) - 1; end if
                if (p_glb > 0) then; clo(3) = amr_region_lo_all(3, k) + 1; chi(3) = amr_region_hi_all(3, k) - 1; end if
                call s_amr_shell_slabs(plo, phi, clo, chi, nsh, shb1, she1, shb2, she2, shb3, she3, scells)
            end if
            do idx = 1, amr_ovl_gather_n(k)
                r = amr_ovl_gather(idx, k)
                if (r == proc_rank) cycle
                call s_amr_rank_coarse_range(r, crlo, crhi)
                call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
                if (do_pbmv) then
                    msl = 1
                    tb1(1) = bl(1); te1(1) = bh(1); tb2(1) = bl(2); te2(1) = bh(2); tb3(1) = bl(3); te3(1) = bh(3)
                else
                    call s_amr_shell_clip(nsh, shb1, she1, shb2, she2, shb3, she3, bl, bh, msl, tb1, te1, tb2, te2, tb3, te3, &
                                          & scells)
                    if (msl == 0) cycle
                end if
                do isl = 1, msl
                    bl = [tb1(isl), tb2(isl), tb3(isl)]; bh = [te1(isl), te2(isl), te3(isl)]
                    qsz = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                    psz = cellsz*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                    if (amr_fw_map(r) == 0) then
                        amr_fw_rnp = amr_fw_rnp + 1
                        call s_amr_fw_szi(amr_fw_rprank, amr_fw_rnp); call s_amr_fw_szi(amr_fw_rqsz, amr_fw_rnp)
                        call s_amr_fw_szi(amr_fw_rpsz, amr_fw_rnp); call s_amr_fw_szi(amr_fw_rnxp, amr_fw_rnp)
                        call s_amr_fw_szi(amr_fw_rqbase, amr_fw_rnp); call s_amr_fw_szi(amr_fw_rpbase, amr_fw_rnp)
                        amr_fw_map(r) = amr_fw_rnp
                        amr_fw_rprank(amr_fw_rnp) = r
                    end if
                    amr_fw_rnx = amr_fw_rnx + 1
                    call s_amr_fw_szi(amr_fw_rblk, amr_fw_rnx); call s_amr_fw_szi3(amr_fw_rbl, amr_fw_rnx)
                    call s_amr_fw_szi3(amr_fw_rbh, amr_fw_rnx); call s_amr_fw_szi(amr_fw_rpi, amr_fw_rnx)
                    call s_amr_fw_szi(amr_fw_rqo, amr_fw_rnx); call s_amr_fw_szi(amr_fw_rpo, amr_fw_rnx)
                    amr_fw_rblk(amr_fw_rnx) = k; amr_fw_rbl(:,amr_fw_rnx) = bl; amr_fw_rbh(:,amr_fw_rnx) = bh
                    amr_fw_rpi(amr_fw_rnx) = amr_fw_map(r)
                    amr_fw_rqo(amr_fw_rnx) = amr_fw_pq(r) + amr_fw_nx(r)*XA_NH
                    amr_fw_rpo(amr_fw_rnx) = amr_fw_pp(r) + amr_fw_nx(r)*XA_NH
                    amr_fw_pq(r) = amr_fw_pq(r) + qsz
                    amr_fw_pp(r) = amr_fw_pp(r) + psz
                    amr_fw_nx(r) = amr_fw_nx(r) + 1
                end do
            end do
        end do
        qbase = 0; pbase = 0
        do ip = 1, amr_fw_rnp
            r = amr_fw_rprank(ip)
            amr_fw_rnxp(ip) = amr_fw_nx(r)
            amr_fw_rqsz(ip) = amr_fw_pq(r) + amr_fw_nx(r)*XA_NH
            amr_fw_rpsz(ip) = amr_fw_pp(r) + amr_fw_nx(r)*XA_NH
            amr_fw_rqbase(ip) = qbase; qbase = qbase + amr_fw_rqsz(ip)
            amr_fw_rpbase(ip) = pbase; pbase = pbase + amr_fw_rpsz(ip)
            amr_fw_map(r) = 0; amr_fw_nx(r) = 0; amr_fw_pq(r) = 0; amr_fw_pp(r) = 0
        end do
        call s_amr_fw_szr(amr_fw_rq, qbase, amr_fw_dev)
        if (do_pbmv) call s_amr_fw_szr(amr_fw_rp, pbase, amr_fw_dev)
        nreq = (amr_fw_snp + amr_fw_rnp)*(1 + merge(1, 0, do_pbmv))
        call s_amr_fw_szi(amr_fw_req, nreq); call s_amr_fw_szi(amr_fw_reqw, nreq)
        call s_phase_toc(PH_GWPLAN)

        ! post all recvs, then pack all sends (device kernels into contiguous pool slices), then post all sends, then one
        ! waitall. [amr-xa] records payload words only, so the family totals are independent of the message aggregation.
        nreq = 0
#ifdef MFC_MPI
        do ip = 1, amr_fw_rnp
            sq = f_amr_m1_seq(amr_fw_rprank(ip), 2); tq = f_amr_m1_tag(3, sq)
            call s_xa_rec(XA_F1W_RCV, 2, amr_fw_rqsz(ip) - amr_fw_rnxp(ip)*XA_NH, tq, peer=amr_fw_rprank(ip), &
                          & key=amr_fw_rnxp(ip), seq=sq)
            nreq = nreq + 1; amr_fw_reqw(nreq) = amr_fw_rqsz(ip)
            if (amr_fw_dev) then
                #:call GPU_HOST_DATA(use_device_addr='[amr_fw_rq]')
                    call MPI_IRECV(amr_fw_rq(amr_fw_rqbase(ip) + 1), amr_fw_rqsz(ip), mpi_p, amr_fw_rprank(ip), tq, &
                                   & MPI_COMM_WORLD, amr_fw_req(nreq), ierr)
                #:endcall GPU_HOST_DATA
            else
                call MPI_IRECV(amr_fw_rq(amr_fw_rqbase(ip) + 1), amr_fw_rqsz(ip), mpi_p, amr_fw_rprank(ip), tq, MPI_COMM_WORLD, &
                               & amr_fw_req(nreq), ierr)
            end if
            if (do_pbmv) then
                sq = f_amr_m1_seq(amr_fw_rprank(ip), 2); tp = f_amr_m1_tag(4, sq)
                call s_xa_rec(XA_F3W_RCV, 2, amr_fw_rpsz(ip) - amr_fw_rnxp(ip)*XA_NH, tp, peer=amr_fw_rprank(ip), &
                              & key=amr_fw_rnxp(ip), seq=sq)
                nreq = nreq + 1; amr_fw_reqw(nreq) = amr_fw_rpsz(ip)
                if (amr_fw_dev) then
                    #:call GPU_HOST_DATA(use_device_addr='[amr_fw_rp]')
                        call MPI_IRECV(amr_fw_rp(amr_fw_rpbase(ip) + 1), amr_fw_rpsz(ip), mpi_p, amr_fw_rprank(ip), tp, &
                                       & MPI_COMM_WORLD, amr_fw_req(nreq), ierr)
                    #:endcall GPU_HOST_DATA
                else
                    call MPI_IRECV(amr_fw_rp(amr_fw_rpbase(ip) + 1), amr_fw_rpsz(ip), mpi_p, amr_fw_rprank(ip), tp, &
                                   & MPI_COMM_WORLD, amr_fw_req(nreq), ierr)
                end if
            end if
        end do
#endif
        call s_phase_tic(PH_GWPACK)
        if (fuse .and. amr_fw_snx > 0) then
            ! one launch for the whole send list; the debug identity headers are written after it, because the fused copyout
            ! covers the pool prefix (payload and header words) and would otherwise clobber host-written headers.
            ! Load-bearing: this copies out the whole pool prefix, and map(from:) leaves any word the kernel did not
            ! write as uninitialised device memory on the host. It is safe only because the pool is exactly tiled
            ! (sqtot == qbase), so every word in 1:sqtot is written. Any padding or alignment in the pool would ship
            ! garbage on the wire, silently: the MFC_DEBUG NaN poison covers the patch, not the pool.
            call s_amr_fx_plan(amr_fw_sbl, amr_fw_sbh, amr_fw_sqbase, amr_fw_sqo, amr_fw_spi, amr_fw_snx)
            call s_amr_fx_pack_box(q_cons_coarse, 1, amr_fw_snx, o1, o2, o3, amr_fx_pl(:,1:amr_fw_snx), &
                                   & amr_fx_pre(1:amr_fw_snx + 1), amr_fw_sq(1:sqtot))
            if (XA_NH > 0) then
                do ix = 1, amr_fw_snx
                    boff = amr_fx_pl(7, ix) - XA_NH
                    call s_xa_hdr_pack(amr_fw_sq(boff + 1:boff + XA_NH), XA_F1W_SND, amr_fw_sblk(ix), amr_fw_sbl(:,ix), &
                                       & amr_fw_sbh(:,ix))
                end do
            end if
        else
            do ix = 1, amr_fw_snx
                bl = amr_fw_sbl(:,ix); bh = amr_fw_sbh(:,ix)
                qsz = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                boff = amr_fw_sqbase(amr_fw_spi(ix)) + amr_fw_sqo(ix)
                if (XA_NH > 0) call s_xa_hdr_pack(amr_fw_sq(boff + 1:boff + XA_NH), XA_F1W_SND, amr_fw_sblk(ix), bl, bh)
                call s_amr_pack_box_device(q_cons_coarse, bl, bh, o1, o2, o3, amr_fw_sq(boff + XA_NH + 1:boff + XA_NH + qsz))
                if (do_pbmv) then
                    psz = cellsz*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                    boff = amr_fw_spbase(amr_fw_spi(ix)) + amr_fw_spo(ix)
                    if (XA_NH > 0) call s_xa_hdr_pack(amr_fw_sp(boff + 1:boff + XA_NH), XA_F3W_SND, amr_fw_sblk(ix), bl, bh)
                    call s_amr_pack_box_pbmv_device(pb_in, mv_in, bl, bh, o1, o2, o3, &
                                                    & amr_fw_sp(boff + XA_NH + 1:boff + XA_NH + psz))
                end if
            end do
        end if
        call s_phase_toc(PH_GWPACK)
#ifdef MFC_MPI
        do ip = 1, amr_fw_snp
            sq = f_amr_m1_seq(amr_fw_sprank(ip), 1); tq = f_amr_m1_tag(3, sq)
            call s_xa_rec(XA_F1W_SND, 1, amr_fw_sqsz(ip) - amr_fw_snxp(ip)*XA_NH, tq, peer=amr_fw_sprank(ip), &
                          & key=amr_fw_snxp(ip), seq=sq)
            nreq = nreq + 1; amr_fw_reqw(nreq) = -1
            if (amr_fw_dev) then
                #:call GPU_HOST_DATA(use_device_addr='[amr_fw_sq]')
                    call MPI_ISEND(amr_fw_sq(amr_fw_sqbase(ip) + 1), amr_fw_sqsz(ip), mpi_p, amr_fw_sprank(ip), tq, &
                                   & MPI_COMM_WORLD, amr_fw_req(nreq), ierr)
                #:endcall GPU_HOST_DATA
            else
                call MPI_ISEND(amr_fw_sq(amr_fw_sqbase(ip) + 1), amr_fw_sqsz(ip), mpi_p, amr_fw_sprank(ip), tq, MPI_COMM_WORLD, &
                               & amr_fw_req(nreq), ierr)
            end if
            if (do_pbmv) then
                sq = f_amr_m1_seq(amr_fw_sprank(ip), 1); tp = f_amr_m1_tag(4, sq)
                call s_xa_rec(XA_F3W_SND, 1, amr_fw_spsz(ip) - amr_fw_snxp(ip)*XA_NH, tp, peer=amr_fw_sprank(ip), &
                              & key=amr_fw_snxp(ip), seq=sq)
                nreq = nreq + 1; amr_fw_reqw(nreq) = -1
                if (amr_fw_dev) then
                    #:call GPU_HOST_DATA(use_device_addr='[amr_fw_sp]')
                        call MPI_ISEND(amr_fw_sp(amr_fw_spbase(ip) + 1), amr_fw_spsz(ip), mpi_p, amr_fw_sprank(ip), tp, &
                                       & MPI_COMM_WORLD, amr_fw_req(nreq), ierr)
                    #:endcall GPU_HOST_DATA
                else
                    call MPI_ISEND(amr_fw_sp(amr_fw_spbase(ip) + 1), amr_fw_spsz(ip), mpi_p, amr_fw_sprank(ip), tp, &
                                   & MPI_COMM_WORLD, amr_fw_req(nreq), ierr)
                end if
            end if
        end do
        call s_phase_tic(PH_GWWAIT)
        if (nreq > 0) then
#ifdef MFC_DEBUG
            block
                integer :: st(MPI_STATUS_SIZE, nreq), gotw, q
                call s_wait_tic()
                call MPI_WAITALL(nreq, amr_fw_req, st, ierr)
                call s_wait_toc(WT_GATHER)
                do q = 1, nreq
                    if (amr_fw_reqw(q) < 0) cycle
                    call MPI_GET_COUNT(st(:,q), mpi_p, gotw, ierr)
                    @:ASSERT(gotw == amr_fw_reqw(q), "stage-fill wave: received message length differs from the plan")
                end do
            end block
#else
            call s_wait_tic()
            call MPI_WAITALL(nreq, amr_fw_req, MPI_STATUSES_IGNORE, ierr)
            call s_wait_toc(WT_GATHER)
#endif
        end if
        call s_phase_toc(PH_GWWAIT)
#endif
        call s_phase_toc(PH_GATHER)

        ! consume, ascending slot order: per owned box, patch frame + own-box device copy + per-slab device unpack (recv
        ! transfers were appended box-major, so each box's slabs are the next contiguous run), then the ghost fills.
        ix = 1
        if (fuse .and. amr_fw_rnx > 0) call s_amr_fx_plan(amr_fw_rbl, amr_fw_rbh, amr_fw_rqbase, amr_fw_rqo, amr_fw_rpi, amr_fw_rnx)
        call s_amr_refresh_my_blocks()
        do kk2 = 1, amr_n_my  ! owned list; level filter kept (list carries all owned levels)
            k = amr_my_blk(kk2)
            if (amr_block_level(k) /= 1) cycle
            call s_amr_select_slot(k)
            if (.not. amr_rank_owns_block) cycle  ! belt-and-braces; list guarantees ownership
            call s_phase_tic(PH_GATHER)
            call s_wait_tic()
            amr_cpat_off = 0
            amr_cpat_off(1) = amr_region_lo_all(1, k) - amr_cpat_mar
            if (n_glb > 0) amr_cpat_off(2) = amr_region_lo_all(2, k) - amr_cpat_mar
            if (p_glb > 0) amr_cpat_off(3) = amr_region_lo_all(3, k) - amr_cpat_mar
            v1hi = (amr_region_hi_all(1, k) - amr_region_lo_all(1, k)) + 2*amr_cpat_mar
            v2hi = 0; v3hi = 0
            if (n_glb > 0) v2hi = (amr_region_hi_all(2, k) - amr_region_lo_all(2, k)) + 2*amr_cpat_mar
            if (p_glb > 0) v3hi = (amr_region_hi_all(3, k) - amr_region_lo_all(3, k)) + 2*amr_cpat_mar
            plo = amr_cpat_off
            phi(1) = plo(1) + v1hi; phi(2) = plo(2) + v2hi; phi(3) = plo(3) + v3hi
            call s_amr_rank_coarse_range(proc_rank, crlo, crhi)
            call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
            call s_wait_toc(WT_HSLOT)
            if (do_pbmv) then
                call s_amr_gather_own_box_device(q_cons_coarse, bl, bh, o1, o2, o3)
                call s_amr_gather_own_box_pbmv_device(pb_in, mv_in, bl, bh, o1, o2, o3)
            else
#ifdef MFC_DEBUG
                ! validation arm: flood the patch with NaN before the clipped writes, so a consumer read of any
                ! unshipped cell (core or a missed shell slab) NaNs the ghost fill within a step
                call s_amr_poison_patch_device(v1hi, v2hi, v3hi)
#endif
                clo = 0; chi = 0
                clo(1) = amr_region_lo_all(1, k) + 1; chi(1) = amr_region_hi_all(1, k) - 1
                if (n_glb > 0) then; clo(2) = amr_region_lo_all(2, k) + 1; chi(2) = amr_region_hi_all(2, k) - 1; end if
                if (p_glb > 0) then; clo(3) = amr_region_lo_all(3, k) + 1; chi(3) = amr_region_hi_all(3, k) - 1; end if
                call s_wait_tic()
                call s_amr_shell_slabs(plo, phi, clo, chi, nsh, shb1, she1, shb2, she2, shb3, she3, scells)
                call s_amr_shell_clip(nsh, shb1, she1, shb2, she2, shb3, she3, bl, bh, msl, tb1, te1, tb2, te2, tb3, te3, scells)
                call s_wait_toc(WT_HSHELL)
                call s_wait_tic()
                if (msl > 0) call s_amr_gather_own_shell_device(q_cons_coarse, msl, tb1, te1, tb2, te2, tb3, te3, o1, o2, o3)
                call s_wait_toc(WT_HOWN)
            end if
            call s_wait_tic()
            do while (ix <= amr_fw_rnx)
                if (amr_fw_rblk(ix) /= k) exit
                if (fuse) then
                    call s_amr_fx_run(k, amr_fw_rblk, amr_fw_rnx, ix, ie)
                    boff = amr_fx_pl(7, ix) - XA_NH
                    if (XA_NH > 0) then
                        do jx = ix, ie
                            call s_xa_hdr_check(amr_fw_rq(amr_fx_pl(7, jx) - XA_NH + 1:amr_fx_pl(7, jx)), XA_F1W_SND, k, &
                                                & amr_fw_rbl(:,jx), amr_fw_rbh(:,jx))
                        end do
                    end if
                    call s_amr_fx_unpack(ix, ie, boff, amr_cpat_off(1), amr_cpat_off(2), amr_cpat_off(3), amr_fx_pl(:, &
                                         & 1:amr_fw_rnx), amr_fx_pre(1:amr_fw_rnx + 1), amr_fw_rq(boff + 1:amr_fx_pl(7, &
                                         & ie) + amr_fx_pre(ie + 1) - amr_fx_pre(ie)))
                    ix = ie + 1
                    cycle
                end if
                bl = amr_fw_rbl(:,ix); bh = amr_fw_rbh(:,ix)
                qsz = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                boff = amr_fw_rqbase(amr_fw_rpi(ix)) + amr_fw_rqo(ix)
                if (XA_NH > 0) call s_xa_hdr_check(amr_fw_rq(boff + 1:boff + XA_NH), XA_F1W_SND, k, bl, bh)
                call s_amr_unpack_box_device(bl, bh, amr_fw_rq(boff + XA_NH + 1:boff + XA_NH + qsz))
                if (do_pbmv) then
                    psz = cellsz*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                    boff = amr_fw_rpbase(amr_fw_rpi(ix)) + amr_fw_rpo(ix)
                    if (XA_NH > 0) call s_xa_hdr_check(amr_fw_rp(boff + 1:boff + XA_NH), XA_F3W_SND, k, bl, bh)
                    call s_amr_unpack_box_pbmv_device(bl, bh, amr_fw_rp(boff + XA_NH + 1:boff + XA_NH + psz))
                end if
                ix = ix + 1
            end do
            call s_wait_toc(WT_HUNPK)
            call s_phase_toc(PH_GATHER)
            if (rank_time_wrt) call s_rank_time_tic()
            call s_phase_tic(PH_GFILL)
            call s_wait_tic()
            call s_amr_fill_fine_ghosts_cons(amr_cg, amr_loc_of(amr_cur))
            call s_wait_toc(WT_HFILL)
            call s_phase_toc(PH_GFILL)
            if (do_pbmv) call s_amr_fill_fine_ghosts_pbmv(amr_cg_pb, amr_cg_mv, amr_slots(amr_cur)%pb_f%sf, &
                & amr_slots(amr_cur)%mv_f%sf)
            if (rank_time_wrt) call s_rank_time_toc()
        end do
        @:ASSERT(ix == amr_fw_rnx + 1, "stage-fill wave: unconsumed recv transfers")

    end subroutine s_amr_stage_fill_wave

    !> The parent-fill wave's per-box transfer list in the patch-local frame: the padded patch's hollow-shell slabs (the runtime
    !! consumer is the amr_cg ghost fill, which never reads the open interior of the parent footprint [mar+1, w-mar-1], so it never
    !! ships), or the single full patch when non-polytropic QBMM keeps the full-box contract. Send walk, recv walk, and consume all
    !! derive the list here, so the wire layout cannot drift between sides.
    impure subroutine s_amr_parent_shell(w1, w2, w3, full, msl, tb1, te1, tb2, te2, tb3, te3)

        integer, intent(in)  :: w1, w2, w3
        logical, intent(in)  :: full
        integer, intent(out) :: msl, tb1(6), te1(6), tb2(6), te2(6), tb3(6), te3(6)
        integer              :: clo(3), chi(3), scells

        if (full) then
            msl = 1
            tb1(1) = 0; te1(1) = w1; tb2(1) = 0; te2(1) = w2; tb3(1) = 0; te3(1) = w3
        else
            clo = 0; chi = 0
            clo(1) = amr_cpat_mar + 1; chi(1) = w1 - amr_cpat_mar - 1
            if (n_glb > 0) then; clo(2) = amr_cpat_mar + 1; chi(2) = w2 - amr_cpat_mar - 1; end if
            if (p_glb > 0) then; clo(3) = amr_cpat_mar + 1; chi(3) = w3 - amr_cpat_mar - 1; end if
            call s_amr_shell_slabs([0, 0, 0], [w1, w2, w3], clo, chi, msl, tb1, te1, tb2, te2, tb3, te3, scells)
        end if

    end subroutine s_amr_parent_shell

    !> Per-step level-lev fill as one exchange wave: the F2 parent gather for every level-lev block in one aggregated exchange. Each
    !! split child is its s_amr_parent_shell transfer list (ring-clipped shell slabs, or one full patch under the pbmv contract)
    !! from its parent's owner to its own owner, so the plan is a pair list, not an overlap map. Same skeleton as
    !! s_amr_stage_fill_wave (whose scratch arrays it reuses; the two never overlap in time): plans from replicated metadata
    !! (f_amr_parent_block + s_amr_parent_foot + amr_block_owner only; the per-owner mirrors lag and are empty on non-owners),
    !! recvs-packs-sends-one-WAITALL, box-major consume through the single amr_cg. Called per level ascending, so a level-(lev-1)
    !! parent's own ghost fill is complete before this wave reads its interior. Co-located parent-child is a consume-phase device
    !! copy with no wire transfer. The regrid uses the chunked F2 path; init/static use the per-box s_amr_gather_from_parent.
    impure subroutine s_amr_parent_fill_wave(lev)

        integer, intent(in) :: lev
        integer             :: k, r, ix, ip, pblk, powner, cowner, boxsz, tq, sq, nreq, qbase, ierr, kk
        integer             :: w1, w2, w3, plo(3), phi(3), boff, bl(3), bh(3), sqtot, ie, jx
        integer             :: msl, isl
        integer             :: tb1(6), te1(6), tb2(6), te2(6), tb3(6), te3(6)
        logical             :: do_pbmv

        if (amr_num_blocks <= 0) return
        @:ASSERT(amr_gsnd_n == 0, "parent-fill wave: the deferred gather-send pool must be drained")
        do_pbmv = qbmm .and. .not. polytropic

        call s_amr_m1_wave_open(2)

        call s_phase_tic(PH_GATHER)
        ! the lag guard must visit every level-lev block (it checks blocks this rank does not own), so it keeps its own global
        ! scan, gated on the one configuration that needs it
        if (bubbles_lagrange) then
            do k = 1, amr_num_blocks
                if (amr_block_level(k) /= lev) cycle
                call s_amr_select_slot(k)
                call s_amr_check_lag_clear()
            end do
        end if
        call s_amr_refresh_lists()
        ! send side: every level-lev block whose parent I own but whose child-owner is another rank (amr_fch_blk narrowed to lev)
        amr_fw_snx = 0; amr_fw_snp = 0
        if (.not. allocated(amr_fw_map)) then
            allocate (amr_fw_map(0:num_procs - 1), amr_fw_nx(0:num_procs - 1), amr_fw_pq(0:num_procs - 1), &
                      & amr_fw_pp(0:num_procs - 1))
            amr_fw_map = 0; amr_fw_nx = 0; amr_fw_pq = 0; amr_fw_pp = 0
        end if
        do kk = 1, amr_n_fch
            k = amr_fch_blk(kk)
            if (amr_block_level(k) /= lev) cycle
            call s_amr_select_slot(k)
            pblk = amr_parent_blk(k)
            powner = amr_block_owner(pblk); cowner = amr_block_owner(k)
            if (powner == cowner .or. powner /= proc_rank) cycle
            call s_amr_parent_foot(k, pblk, plo, phi)
            w1 = (phi(1) - plo(1)) + 2*amr_cpat_mar
            w2 = 0; w3 = 0
            if (n_glb > 0) w2 = (phi(2) - plo(2)) + 2*amr_cpat_mar
            if (p_glb > 0) w3 = (phi(3) - plo(3)) + 2*amr_cpat_mar
            call s_amr_parent_shell(w1, w2, w3, do_pbmv, msl, tb1, te1, tb2, te2, tb3, te3)
            do isl = 1, msl
                bl = [tb1(isl), tb2(isl), tb3(isl)]; bh = [te1(isl), te2(isl), te3(isl)]
                boxsz = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                if (amr_fw_map(cowner) == 0) then
                    amr_fw_snp = amr_fw_snp + 1
                    call s_amr_fw_szi(amr_fw_sprank, amr_fw_snp); call s_amr_fw_szi(amr_fw_sqsz, amr_fw_snp)
                    call s_amr_fw_szi(amr_fw_snxp, amr_fw_snp); call s_amr_fw_szi(amr_fw_sqbase, amr_fw_snp)
                    amr_fw_map(cowner) = amr_fw_snp
                    amr_fw_sprank(amr_fw_snp) = cowner
                end if
                amr_fw_snx = amr_fw_snx + 1
                call s_amr_fw_szi(amr_fw_sblk, amr_fw_snx); call s_amr_fw_szi3(amr_fw_sbl, amr_fw_snx)
                call s_amr_fw_szi3(amr_fw_sbh, amr_fw_snx); call s_amr_fw_szi(amr_fw_spi, amr_fw_snx)
                call s_amr_fw_szi(amr_fw_sqo, amr_fw_snx)
                amr_fw_sblk(amr_fw_snx) = k; amr_fw_sbl(:,amr_fw_snx) = bl; amr_fw_sbh(:,amr_fw_snx) = bh
                amr_fw_spi(amr_fw_snx) = amr_fw_map(cowner)
                amr_fw_sqo(amr_fw_snx) = amr_fw_pq(cowner) + amr_fw_nx(cowner)*XA_NH
                amr_fw_pq(cowner) = amr_fw_pq(cowner) + boxsz
                amr_fw_nx(cowner) = amr_fw_nx(cowner) + 1
            end do
        end do
        qbase = 0
        do ip = 1, amr_fw_snp
            r = amr_fw_sprank(ip)
            amr_fw_snxp(ip) = amr_fw_nx(r)
            amr_fw_sqsz(ip) = amr_fw_pq(r) + amr_fw_nx(r)*XA_NH
            amr_fw_sqbase(ip) = qbase; qbase = qbase + amr_fw_sqsz(ip)
            amr_fw_map(r) = 0; amr_fw_nx(r) = 0; amr_fw_pq(r) = 0
        end do
        call s_amr_fw_szr(amr_fw_sq, qbase, amr_fw_dev)
        sqtot = qbase

        ! recv side: every level-lev block I own whose parent lives on another rank; the box's shell-slab transfers (or its
        ! one full-patch transfer under the pbmv contract). Both sides enumerate boxes ascending, slabs in the fixed
        ! s_amr_parent_shell order, with per-rank running offsets, so the wire layout agrees with no metadata exchange.
        amr_fw_rnx = 0; amr_fw_rnp = 0
        call s_amr_refresh_my_blocks()
        do kk = 1, amr_n_my
            k = amr_my_blk(kk)
            if (amr_block_level(k) /= lev) cycle
            pblk = amr_parent_blk(k)
            powner = amr_block_owner(pblk)
            if (powner == proc_rank) cycle
            call s_amr_parent_foot(k, pblk, plo, phi)
            w1 = (phi(1) - plo(1)) + 2*amr_cpat_mar
            w2 = 0; w3 = 0
            if (n_glb > 0) w2 = (phi(2) - plo(2)) + 2*amr_cpat_mar
            if (p_glb > 0) w3 = (phi(3) - plo(3)) + 2*amr_cpat_mar
            call s_amr_parent_shell(w1, w2, w3, do_pbmv, msl, tb1, te1, tb2, te2, tb3, te3)
            do isl = 1, msl
                bl = [tb1(isl), tb2(isl), tb3(isl)]; bh = [te1(isl), te2(isl), te3(isl)]
                boxsz = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                if (amr_fw_map(powner) == 0) then
                    amr_fw_rnp = amr_fw_rnp + 1
                    call s_amr_fw_szi(amr_fw_rprank, amr_fw_rnp); call s_amr_fw_szi(amr_fw_rqsz, amr_fw_rnp)
                    call s_amr_fw_szi(amr_fw_rnxp, amr_fw_rnp); call s_amr_fw_szi(amr_fw_rqbase, amr_fw_rnp)
                    amr_fw_map(powner) = amr_fw_rnp
                    amr_fw_rprank(amr_fw_rnp) = powner
                end if
                amr_fw_rnx = amr_fw_rnx + 1
                call s_amr_fw_szi(amr_fw_rblk, amr_fw_rnx); call s_amr_fw_szi3(amr_fw_rbl, amr_fw_rnx)
                call s_amr_fw_szi3(amr_fw_rbh, amr_fw_rnx); call s_amr_fw_szi(amr_fw_rpi, amr_fw_rnx)
                call s_amr_fw_szi(amr_fw_rqo, amr_fw_rnx)
                amr_fw_rblk(amr_fw_rnx) = k; amr_fw_rbl(:,amr_fw_rnx) = bl; amr_fw_rbh(:,amr_fw_rnx) = bh
                amr_fw_rpi(amr_fw_rnx) = amr_fw_map(powner)
                amr_fw_rqo(amr_fw_rnx) = amr_fw_pq(powner) + amr_fw_nx(powner)*XA_NH
                amr_fw_pq(powner) = amr_fw_pq(powner) + boxsz
                amr_fw_nx(powner) = amr_fw_nx(powner) + 1
            end do
        end do
        qbase = 0
        do ip = 1, amr_fw_rnp
            r = amr_fw_rprank(ip)
            amr_fw_rnxp(ip) = amr_fw_nx(r)
            amr_fw_rqsz(ip) = amr_fw_pq(r) + amr_fw_nx(r)*XA_NH
            amr_fw_rqbase(ip) = qbase; qbase = qbase + amr_fw_rqsz(ip)
            amr_fw_map(r) = 0; amr_fw_nx(r) = 0; amr_fw_pq(r) = 0
        end do
        call s_amr_fw_szr(amr_fw_rq, qbase, amr_fw_dev)
        nreq = amr_fw_snp + amr_fw_rnp
        call s_amr_fw_szi(amr_fw_req, nreq); call s_amr_fw_szi(amr_fw_reqw, nreq)

        nreq = 0
#ifdef MFC_MPI
        do ip = 1, amr_fw_rnp
            sq = f_amr_m1_seq(amr_fw_rprank(ip), 2); tq = f_amr_m1_tag(2, sq)
            call s_xa_rec(XA_F2W_RCV, 2, amr_fw_rqsz(ip) - amr_fw_rnxp(ip)*XA_NH, tq, peer=amr_fw_rprank(ip), &
                          & key=amr_fw_rnxp(ip), seq=sq)
            nreq = nreq + 1; amr_fw_reqw(nreq) = amr_fw_rqsz(ip)
            if (amr_fw_dev) then
                #:call GPU_HOST_DATA(use_device_addr='[amr_fw_rq]')
                    call MPI_IRECV(amr_fw_rq(amr_fw_rqbase(ip) + 1), amr_fw_rqsz(ip), mpi_p, amr_fw_rprank(ip), tq, &
                                   & MPI_COMM_WORLD, amr_fw_req(nreq), ierr)
                #:endcall GPU_HOST_DATA
            else
                call MPI_IRECV(amr_fw_rq(amr_fw_rqbase(ip) + 1), amr_fw_rqsz(ip), mpi_p, amr_fw_rprank(ip), tq, MPI_COMM_WORLD, &
                               & amr_fw_req(nreq), ierr)
            end if
        end do
#endif
        ! pack: the parent-patch pack kernel reads amr_cpat_off from module scope, so set the child's frame per transfer
        ! (the consume loop recomputes it per box; phases are sequential, so the global is single-writer at any time).
        ! sbl/sbh hold the transfer's patch-local slab bounds; the frame comes from the box's parent foot.
        if (amr_device_pack .and. amr_fw_snx > 0) then
            ! one launch for the whole send list: the per-transfer parent slot and child patch frame ride the plan rows (the
            ! per-box path holds them in amr_cpat_off). Debug headers go in after the fused copyout, which covers the pool
            ! prefix.
            call s_amr_fx_plan(amr_fw_sbl, amr_fw_sbh, amr_fw_sqbase, amr_fw_sqo, amr_fw_spi, amr_fw_snx)
            do ix = 1, amr_fw_snx
                k = amr_fw_sblk(ix)
                call s_amr_parent_foot(k, amr_parent_blk(k), plo, phi)
                amr_fx_pl(8, ix) = amr_loc_of(amr_parent_blk(k))
                amr_fx_pl(9, ix) = plo(1) - amr_cpat_mar
                amr_fx_pl(10, ix) = 0; amr_fx_pl(11, ix) = 0
                if (n_glb > 0) amr_fx_pl(10, ix) = plo(2) - amr_cpat_mar
                if (p_glb > 0) amr_fx_pl(11, ix) = plo(3) - amr_cpat_mar
            end do
            call s_amr_fx_pack_parent(1, amr_fw_snx, amr_fx_pl(:,1:amr_fw_snx), amr_fx_pre(1:amr_fw_snx + 1), amr_fw_sq(1:sqtot))
            if (XA_NH > 0) then
                do ix = 1, amr_fw_snx
                    boff = amr_fx_pl(7, ix) - XA_NH
                    call s_xa_hdr_pack(amr_fw_sq(boff + 1:boff + XA_NH), XA_F2W_SND, amr_fw_sblk(ix), amr_fw_sbl(:,ix), &
                                       & amr_fw_sbh(:,ix))
                end do
            end if
        else
            do ix = 1, amr_fw_snx
                k = amr_fw_sblk(ix)
                call s_amr_parent_foot(k, amr_parent_blk(k), plo, phi)
                amr_cpat_off = 0
                amr_cpat_off(1) = plo(1) - amr_cpat_mar
                if (n_glb > 0) amr_cpat_off(2) = plo(2) - amr_cpat_mar
                if (p_glb > 0) amr_cpat_off(3) = plo(3) - amr_cpat_mar
                bl = amr_fw_sbl(:,ix); bh = amr_fw_sbh(:,ix)
                boxsz = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                boff = amr_fw_sqbase(amr_fw_spi(ix)) + amr_fw_sqo(ix)
                call s_amr_pack_parent_box_device_cons(amr_loc_of(amr_parent_blk(k)), bl, bh, &
                                                       & amr_fw_sq(boff + XA_NH + 1:boff + XA_NH + boxsz))
                if (XA_NH > 0) call s_xa_hdr_pack(amr_fw_sq(boff + 1:boff + XA_NH), XA_F2W_SND, k, bl, bh)
            end do
        end if
#ifdef MFC_MPI
        do ip = 1, amr_fw_snp
            sq = f_amr_m1_seq(amr_fw_sprank(ip), 1); tq = f_amr_m1_tag(2, sq)
            call s_xa_rec(XA_F2W_SND, 1, amr_fw_sqsz(ip) - amr_fw_snxp(ip)*XA_NH, tq, peer=amr_fw_sprank(ip), &
                          & key=amr_fw_snxp(ip), seq=sq)
            nreq = nreq + 1; amr_fw_reqw(nreq) = -1
            if (amr_fw_dev) then
                #:call GPU_HOST_DATA(use_device_addr='[amr_fw_sq]')
                    call MPI_ISEND(amr_fw_sq(amr_fw_sqbase(ip) + 1), amr_fw_sqsz(ip), mpi_p, amr_fw_sprank(ip), tq, &
                                   & MPI_COMM_WORLD, amr_fw_req(nreq), ierr)
                #:endcall GPU_HOST_DATA
            else
                call MPI_ISEND(amr_fw_sq(amr_fw_sqbase(ip) + 1), amr_fw_sqsz(ip), mpi_p, amr_fw_sprank(ip), tq, MPI_COMM_WORLD, &
                               & amr_fw_req(nreq), ierr)
            end if
        end do
        if (nreq > 0) then
#ifdef MFC_DEBUG
            block
                integer :: st(MPI_STATUS_SIZE, nreq), gotw, q
                call s_wait_tic()
                call MPI_WAITALL(nreq, amr_fw_req, st, ierr)
                call s_wait_toc(WT_PGATHER)
                do q = 1, nreq
                    if (amr_fw_reqw(q) < 0) cycle
                    call MPI_GET_COUNT(st(:,q), mpi_p, gotw, ierr)
                    @:ASSERT(gotw == amr_fw_reqw(q), "parent-fill wave: received message length differs from the plan")
                end do
            end block
#else
            call s_wait_tic()
            call MPI_WAITALL(nreq, amr_fw_req, MPI_STATUSES_IGNORE, ierr)
            call s_wait_toc(WT_PGATHER)
#endif
        end if
#endif
        call s_phase_toc(PH_GATHER)

        ! consume, ascending slot order: per owned level-lev box, patch frame + the shell-slab parent copies (co-located) or
        ! the box's received transfers, then the ghost fills.
        ix = 1
        if (amr_device_pack .and. amr_fw_rnx > 0) call s_amr_fx_plan(amr_fw_rbl, amr_fw_rbh, amr_fw_rqbase, amr_fw_rqo, &
            & amr_fw_rpi, amr_fw_rnx)
        ! amr_own_blk is the amr_owns_all (multi-owner intersection) set, not amr_my_blk, whose single-owner notion would
        ! silently narrow this loop (see the list declarations). Ascending order, so the ix cursor pairing holds.
        do kk = 1, amr_n_own
            k = amr_own_blk(kk)
            if (amr_block_level(k) /= lev) cycle
            call s_amr_select_slot(k)
            if (.not. amr_rank_owns_block) cycle
            call s_phase_tic(PH_GATHER)
            call s_wait_tic()
            pblk = amr_parent_blk(k)
            call s_amr_parent_foot(k, pblk, plo, phi)
            amr_cpat_off = 0
            amr_cpat_off(1) = plo(1) - amr_cpat_mar
            if (n_glb > 0) amr_cpat_off(2) = plo(2) - amr_cpat_mar
            if (p_glb > 0) amr_cpat_off(3) = plo(3) - amr_cpat_mar
            w1 = (phi(1) - plo(1)) + 2*amr_cpat_mar
            w2 = 0; w3 = 0
            if (n_glb > 0) w2 = (phi(2) - plo(2)) + 2*amr_cpat_mar
            if (p_glb > 0) w3 = (phi(3) - plo(3)) + 2*amr_cpat_mar
            call s_wait_toc(WT_HSLOT)
#ifdef MFC_DEBUG
            ! validation arm (mirror of the stepfill clip): NaN-flood the patch before the shell writes land, so a consumer
            ! read of any unshipped cell (the clipped core or a missed slab) NaNs the ghost fill within a step
            if (.not. do_pbmv) call s_amr_poison_patch_device(w1, w2, w3)
#endif
            if (amr_block_owner(pblk) == proc_rank) then
                call s_wait_tic()
                call s_amr_parent_shell(w1, w2, w3, do_pbmv, msl, tb1, te1, tb2, te2, tb3, te3)
                call s_wait_toc(WT_HSHELL)
                call s_wait_tic()
                do isl = 1, msl
                    call s_amr_copy_parent_box_cons(amr_loc_of(pblk), [tb1(isl), tb2(isl), tb3(isl)], [te1(isl), te2(isl), &
                                                    & te3(isl)])
                end do
                call s_wait_toc(WT_HOWN)
            else
                @:ASSERT(ix <= amr_fw_rnx .and. amr_fw_rblk(ix) == k, "parent-fill wave: missing recv transfer")
                call s_wait_tic()
                do while (ix <= amr_fw_rnx)
                    if (amr_fw_rblk(ix) /= k) exit
                    if (amr_device_pack) then
                        ! the box's transfers all come from its one parent owner, so the run is the whole box: one launch
                        call s_amr_fx_run(k, amr_fw_rblk, amr_fw_rnx, ix, ie)
                        boff = amr_fx_pl(7, ix) - XA_NH
                        if (XA_NH > 0) then
                            do jx = ix, ie
                                call s_xa_hdr_check(amr_fw_rq(amr_fx_pl(7, jx) - XA_NH + 1:amr_fx_pl(7, jx)), XA_F2W_SND, k, &
                                                    & amr_fw_rbl(:,jx), amr_fw_rbh(:,jx))
                            end do
                        end if
                        call s_amr_fx_unpack(ix, ie, boff, 0, 0, 0, amr_fx_pl(:,1:amr_fw_rnx), amr_fx_pre(1:amr_fw_rnx + 1), &
                                             & amr_fw_rq(boff + 1:amr_fx_pl(7, ie) + amr_fx_pre(ie + 1) - amr_fx_pre(ie)))
                        ix = ie + 1
                        cycle
                    end if
                    bl = amr_fw_rbl(:,ix); bh = amr_fw_rbh(:,ix)
                    boxsz = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                    boff = amr_fw_rqbase(amr_fw_rpi(ix)) + amr_fw_rqo(ix)
                    if (XA_NH > 0) call s_xa_hdr_check(amr_fw_rq(boff + 1:boff + XA_NH), XA_F2W_SND, k, bl, bh)
                    call s_amr_unpack_parent_box_device(bl, bh, amr_fw_rq(boff + XA_NH + 1:boff + XA_NH + boxsz))
                    ix = ix + 1
                end do
                call s_wait_toc(WT_HUNPK)
            end if
            call s_phase_toc(PH_GATHER)
            if (rank_time_wrt) call s_rank_time_tic()
            call s_phase_tic(PH_GFILL)
            call s_wait_tic()
            call s_amr_fill_fine_ghosts_cons(amr_cg, amr_loc_of(amr_cur))
            call s_wait_toc(WT_HFILL)
            call s_phase_toc(PH_GFILL)
            if (qbmm .and. .not. polytropic) call s_amr_fill_fine_ghosts_pbmv(amr_cg_pb, amr_cg_mv, amr_slots(amr_cur)%pb_f%sf, &
                & amr_slots(amr_cur)%mv_f%sf)
            if (rank_time_wrt) call s_rank_time_toc()
        end do
        @:ASSERT(ix == amr_fw_rnx + 1, "parent-fill wave: unconsumed recv transfers")

    end subroutine s_amr_parent_fill_wave

    !> High-water sizing for the wave's plan scratch. Callers size-then-write at append time, so a grow must preserve the entries
    !! already appended this wave.
    impure subroutine s_amr_fw_szi(a, n)

        integer, allocatable, intent(inout) :: a(:)
        integer, intent(in)                 :: n
        integer, allocatable                :: tmp(:)

        if (.not. allocated(a)) then
            allocate (a(max(n, 64)))
            return
        end if
        if (size(a) >= n) return
        call move_alloc(a, tmp)
        allocate (a(max(n, 2*size(tmp))))
        a(1:size(tmp)) = tmp

    end subroutine s_amr_fw_szi

    impure subroutine s_amr_fw_szi3(a, n)

        integer, allocatable, intent(inout) :: a(:,:)
        integer, intent(in)                 :: n
        integer, allocatable                :: tmp(:,:)

        if (.not. allocated(a)) then
            allocate (a(3, max(n, 64)))
            return
        end if
        if (size(a, 2) >= n) return
        call move_alloc(a, tmp)
        allocate (a(3, max(n, 2*size(tmp, 2))))
        a(:,1:size(tmp, 2)) = tmp

    end subroutine s_amr_fw_szi3

    !> Wire pools: preserving on grow (the F5 waves append debug header slots incrementally; the other waves size once). dev =
    !! .true. keeps the pool device-resident across (re)allocation (amr_fw_dev): the old image is deleted from the device before it
    !! is freed and the new one created after; contents never survive a wave, so nothing is copied.
    impure subroutine s_amr_fw_szr(a, n, dev)

        real(wp), allocatable, intent(inout) :: a(:)
        integer, intent(in)                  :: n
        logical, intent(in), optional        :: dev
        real(wp), allocatable                :: tmp(:)
        logical                              :: on_dev

        on_dev = .false.
        if (present(dev)) on_dev = dev
        if (.not. allocated(a)) then
            allocate (a(max(n, 64)))
            if (on_dev) then
                $:GPU_ENTER_DATA(create='[a]')
            end if
            return
        end if
        if (size(a) >= n) return
        if (on_dev) then
            $:GPU_EXIT_DATA(delete='[a]')
        end if
        call move_alloc(a, tmp)
        allocate (a(max(n, 2*size(tmp))))
        a(1:size(tmp)) = tmp
        if (on_dev) then
            $:GPU_ENTER_DATA(create='[a]')
        end if

    end subroutine s_amr_fw_szr

    !> Device pack (to_buf=T) / unpack (F) of a field's interior block [o+0:o+fm] <-> the contiguous MPI buffer buf, for the P2P
    !! migration + migrated-tile scatter. Follows s_amr_fine_slice: the pack/unpack runs on the device with copyout/copyin moving
    !! only buf host<->device (no strided %sf section in a map clause; flang miscomputes those). buf index runs j fastest then k,l,i
    !! so a matching pack/unpack aligns cell-for-cell. wp buffer, cast to/from stp (identity at double), matching the fine-fine
    !! halo. Two targets, one body: `_st` packs a block out of the flat store (the migration/scatter paths), `_sf` packs the level-0
    !! monolithic q_cons_vf, which is a real scalar_field array and not in the store.
    #:for SFX, TGT in [('st', 'amr_cons_st'), ('sf', '')]
        #:set QB = (lambda ix: TGT + '(o1 + j, o2 + k, o3 + l, ' + ix + ', loc)') if TGT else (lambda ix: 'q(' + ix &
                    & + ')%sf(o1 + j, o2 + k, o3 + l)')
        impure subroutine s_l0_pack_unpack_block_${SFX}$(${'loc' if TGT else 'q'}$, o1, o2, o3, fm1, fm2, fm3, buf, to_buf)

            #:if TGT
                integer, intent(in) :: loc
            #:else
                type(scalar_field), dimension(sys_size), intent(inout) :: q
            #:endif
            integer, intent(in)                 :: o1, o2, o3, fm1, fm2, fm3
            real(wp), intent(inout), contiguous :: buf(:)
            logical, intent(in)                 :: to_buf
            integer                             :: i, j, k, l

            if (to_buf) then
                $:GPU_PARALLEL_LOOP(collapse=4, copyout='[buf]')
                do i = 1, sys_size
                    do l = 0, fm3
                        do k = 0, fm2
                            do j = 0, fm1
                                buf(1 + j + (fm1 + 1)*(k + (fm2 + 1)*(l + (fm3 + 1)*(i - 1)))) = real(${QB('i')}$, wp)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else
                $:GPU_PARALLEL_LOOP(collapse=4, copyin='[buf]')
                do i = 1, sys_size
                    do l = 0, fm3
                        do k = 0, fm2
                            do j = 0, fm1
                                ${QB('i')}$ = real(buf(1 + j + (fm1 + 1)*(k + (fm2 + 1)*(l + (fm3 + 1)*(i - 1)))), stp)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

        end subroutine s_l0_pack_unpack_block_${SFX}$
    #:endfor
end module m_amr_exchange
