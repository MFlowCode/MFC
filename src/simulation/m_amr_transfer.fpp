!>
!!@file
!!@brief Contains module m_amr_transfer

#! AMD OpenMP lane: assert allocatables present on every kernel here (see OMP_DEFAULT_STR). A conditionally allocated module
#! array a kernel names launches only under its allocation's own condition (sw_jac/jac: igr); a kernel naming an unallocated
#! array aborts. Keep it so.
#:set MFC_OMP_PRESENT_ALLOCATABLE = True
#:include 'macros.fpp'

!> @brief Prolongation, restriction and Berger-Colella reflux, including the restrict/reflux/freg waves.
module m_amr_transfer

#ifdef MFC_MPI
    use mpi  !< MPI-IO for the parallel_io AMR restart file
#endif

    use m_derived_types  ! scalar_field, t_box, int_bounds_info
    use m_global_parameters
    use m_mpi_proxy, only: s_mpi_abort
    use m_mpi_common, only: s_mpi_allreduce_integer_max
    use m_amr_registers, only: s_amr_reflux_apply_faces, s_amr_parent_foot, freg, s_amr_reg_prepare, s_amr_reflux_faces
    use m_phase_timing
    use m_amr_xchg_audit  ! per-call-site accounting of every AMR p2p transfer (s_xa_rec + XA_* site ids)
    use m_amr_state
    use m_amr_wave
    use m_amr_distribution
    use m_amr_store
    use m_amr_exchange
    use m_amr_frame

    implicit none

    private
    public :: s_amr_freg_wave, s_amr_reduce_xchg_flag, s_amr_reflux_faces_wave, s_amr_reflux_to_parent, s_amr_restrict_to_parent, &
        & s_amr_restrict_wave, s_interpolate_coarse_to_fine, s_populate_amr_fine, s_restrict_fine_to_coarse, s_set_amr_fine_geometry

contains

    !> True iff rank r refluxes some face of the current block (its subdomain from the computed decomposition, unclipped: the
    !! participation map and both exchange paths gate their freg receives on it).
    pure logical function f_amr_reflux_participates(r) result(part)

        integer, intent(in) :: r
        integer             :: sidx(3), ext(3)
        logical             :: s_lo(3), s_hi(3)

        call s_amr_rank_decomp(r, sidx, ext)
        call s_amr_reflux_faces(sidx, ext, .false., s_lo, s_hi)
        part = any(s_lo .or. s_hi)

    end function f_amr_reflux_participates

    !> The faces of the current block that rank r applies (seam-clipped); the reflux-faces wave ships exactly these, and sender and
    !! receiver derive the identical set from replicated data.
    pure subroutine s_amr_reflux_faces_for(r, s_lo, s_hi)

        integer, intent(in)  :: r
        logical, intent(out) :: s_lo(3), s_hi(3)
        integer              :: sidx(3), ext(3)

        call s_amr_rank_decomp(r, sidx, ext)
        call s_amr_reflux_faces(sidx, ext, .true., s_lo, s_hi)

    end subroutine s_amr_reflux_faces_for

    !> The per-stage level-1 reflux-face exchange as one wave. Every rank walks the level-1 slots ascending: all receives post first
    !! (zero-copy, directly into the freg host mirrors; each box owns a register slot, so no pool is needed), then the owners stage
    !! + multicast, then one waitall, then the receivers push to device. Every message carries its own keyed tag (band 0, per-pair
    !! seq in the (ascending box, ascending dim, lo-then-hi) plan order both sides derive from the same predicates), so matching
    !! does not depend on posting order or MPI non-overtaking. Under MFC_DEBUG the identity headers travel as separate 8-word
    !! companion messages, one per (box, peer) group ahead of its payloads (a prefix cannot ride a zero-copy payload); they are
    !! never recorded in [amr-xa], so the family words stay exactly comparable. The register arrays are sized up front: the apply
    !! can reallocate them, so nothing may post into freg before s_amr_reg_prepare.
    impure subroutine s_amr_reflux_faces_wave()

#ifdef MFC_MPI
        integer :: k, r, idx, ncand, nhr, nhs, j, kk2
        integer :: cand(num_procs), glo(3), ghi(3)
        logical :: s_lo(3), s_hi(3), u_lo(3), u_hi(3)
        logical :: cl(3, num_procs), ch(3, num_procs)

        if (num_procs == 1) return
        call s_amr_reg_prepare()
        call s_amr_wave_open(amr_wave, 0)
        nhr = 0; nhs = 0
        call s_amr_refresh_lists()
        call s_amr_hdr_pools(amr_n_l1p)
        ! receive side: for every level-1 block another rank owns whose faces this rank refluxes, the owner's freg faces land
        ! straight in this rank's register slot (zero-copy); faces this rank does not reflux are poisoned so a stray read shows
        do kk2 = 1, amr_n_l1p
            k = amr_l1p_blk(kk2)
            call s_amr_select_slot(k)
            if (amr_block_owner(k) == proc_rank) cycle
            if (.not. f_amr_reflux_participates(proc_rank)) cycle
            call s_amr_reflux_faces_for(proc_rank, s_lo, s_hi)
            nhr = nhr + 1
            call s_amr_size_int(amr_fw_rblk, nhr)
            amr_fw_rblk(nhr) = k
            call s_amr_hdr_req(2, nhr, amr_block_owner(k), XA_F5W_FACE_RCV, k)
            call s_amr_freg_faces_poison(s_lo, s_hi)
            call s_amr_freg_faces_req(2, amr_block_owner(k), XA_F5W_FACE_RCV, k, s_lo, s_hi)
        end do
        ! send side: every owned level-1 block's faces to each rank that refluxes them
        call s_amr_refresh_my_blocks()
        do kk2 = 1, amr_n_my
            k = amr_my_blk(kk2)
            if (amr_block_level(k) /= 1) cycle
            call s_amr_select_slot(k)
            if (amr_block_owner(k) /= proc_rank) cycle
            glo = merge(amr_region_lo - 1, 0, amr_dim); ghi = merge(amr_region_hi + 1, 0, amr_dim)
            call s_amr_ranks_overlapping(glo, ghi, cand, ncand)
            u_lo = .false.; u_hi = .false.
            do idx = 1, ncand
                r = cand(idx)
                cl(:,idx) = .false.; ch(:,idx) = .false.
                if (r == proc_rank .or. .not. f_amr_reflux_participates(r)) cycle
                call s_amr_reflux_faces_for(r, s_lo, s_hi)
                cl(:,idx) = s_lo; ch(:,idx) = s_hi
                u_lo = u_lo .or. s_lo; u_hi = u_hi .or. s_hi
            end do
            call s_amr_freg_faces_update(.false., u_lo, u_hi)
            do idx = 1, ncand
                r = cand(idx)
                if (r == proc_rank .or. .not. f_amr_reflux_participates(r)) cycle
                nhs = nhs + 1
                call s_amr_hdr_req(1, nhs, r, XA_F5W_FACE_SND, k)
                call s_amr_freg_faces_req(1, r, XA_F5W_FACE_SND, k, cl(:,idx), ch(:,idx))
            end do
        end do
        call s_amr_wave_wait(amr_wave)
        do j = 1, nhr
            k = amr_fw_rblk(j)
            call s_amr_select_slot(k)
            call s_amr_hdr_check(j, XA_F5W_FACE_SND, k)
            call s_amr_reflux_faces_for(proc_rank, s_lo, s_hi)
            call s_amr_freg_faces_update(.true., s_lo, s_hi)
        end do
#endif

    end subroutine s_amr_reflux_faces_wave

    !> The split-ownership level>=2 freg exchange as one wave, run once before the reflux fold (the registers are final after the
    !! advance, and the applies keep their per-box reverse-order position). Same zero-copy, companion-header design as the faces
    !! wave; keyed tags on band 1 (the faces wave is band 0) keep the two disjoint.
    impure subroutine s_amr_freg_wave()

#ifdef MFC_MPI
        integer  :: k, pblk, cowner, powner, nhr, nhs, j, kk2
        real(wp) :: w_lo(3), w_hi(3)

        if (num_procs == 1) return
        call s_amr_reg_prepare()
        call s_amr_refresh_lists()
        call s_amr_wave_open(amr_wave, 1)
        nhr = 0; nhs = 0
        call s_amr_hdr_pools(amr_n_fch)
        ! receive side: the freg faces of every level>=2 child of my parents that another rank owns, straight into the
        ! child's register slot (zero-copy); faces no sibling weight selects are poisoned so a stray read shows up
        do kk2 = 1, amr_n_fch
            k = amr_fch_blk(kk2)
            call s_amr_select_slot(k)
            pblk = amr_parent_blk(k)
            cowner = amr_block_owner(k); powner = amr_block_owner(pblk)
            if (cowner == powner .or. powner /= proc_rank) cycle
            call s_amr_sibling_face_weights(k, pblk, w_lo, w_hi)
            nhr = nhr + 1
            call s_amr_size_int(amr_fw_rblk, nhr)
            amr_fw_rblk(nhr) = k
            call s_amr_hdr_req(2, nhr, cowner, XA_F5W_FREG_RCV, k)
            call s_amr_freg_faces_poison(w_lo > 0._wp, w_hi > 0._wp)
            call s_amr_freg_faces_req(2, cowner, XA_F5W_FREG_RCV, k, w_lo > 0._wp, w_hi > 0._wp)
        end do
        ! send side: my owned level>=2 blocks whose parent lives elsewhere ship the faces the sibling weights select
        call s_amr_refresh_my_blocks()
        do kk2 = 1, amr_n_my
            k = amr_my_blk(kk2)
            if (amr_block_level(k) < 2) cycle
            call s_amr_select_slot(k)
            pblk = amr_parent_blk(k)
            cowner = amr_block_owner(k); powner = amr_block_owner(pblk)
            if (cowner == powner .or. cowner /= proc_rank) cycle
            call s_amr_sibling_face_weights(k, pblk, w_lo, w_hi)
            call s_amr_freg_faces_update(.false., w_lo > 0._wp, w_hi > 0._wp)
            nhs = nhs + 1
            call s_amr_hdr_req(1, nhs, powner, XA_F5W_FREG_SND, k)
            call s_amr_freg_faces_req(1, powner, XA_F5W_FREG_SND, k, w_lo > 0._wp, w_hi > 0._wp)
        end do
        call s_amr_wave_wait(amr_wave)
        do j = 1, nhr
            k = amr_fw_rblk(j)
            call s_amr_select_slot(k)
            pblk = amr_parent_blk(k)
            call s_amr_hdr_check(j, XA_F5W_FREG_SND, k)
            call s_amr_sibling_face_weights(k, pblk, w_lo, w_hi)
            call s_amr_freg_faces_update(.true., w_lo > 0._wp, w_hi > 0._wp)
        end do
#endif

    end subroutine s_amr_freg_wave

    !> Size the header pools for a face wave under the audit: one identity header per received block (nrecv at most) and per (owned
    !! block, peer) pair.
    impure subroutine s_amr_hdr_pools(nrecv)

        integer, intent(in) :: nrecv

        if (XA_NH == 0) return
        call s_amr_refresh_my_blocks()
        call s_amr_size_real(amr_fw_rq, XA_NH*max(nrecv, 1), amr_fw_dev)
        call s_amr_size_real(amr_fw_sq, XA_NH*max(amr_n_my*num_procs, 1), amr_fw_dev)

    end subroutine s_amr_hdr_pools

    !> Under the audit, the identity-header companion message that precedes block k's faces with peer: dir 1 packs and sends header
    !! j from the send pool, 2 posts the receive of header j into the receive pool. Never recorded in [amr-xa].
    impure subroutine s_amr_hdr_req(dir, j, peer, site, k)

        integer, intent(in) :: dir, j, peer, site, k

        if (XA_NH == 0) return
        if (dir == 1) then
            @:ASSERT(size(amr_fw_sq) >= XA_NH*j, "amr_fw_sq header pool sized below the wave's send count")
            call s_xa_hdr_pack(amr_fw_sq(XA_NH*(j - 1) + 1:XA_NH*j), site, k, [0, 0, 0], [0, 0, 0])
            call s_amr_wave_req(amr_wave, 1, amr_fw_sq(XA_NH*(j - 1) + 1:XA_NH*j), XA_NH, peer, site, 0, .false., rec=.false.)
        else
            call s_amr_wave_req(amr_wave, 2, amr_fw_rq(XA_NH*(j - 1) + 1:XA_NH*j), XA_NH, peer, site, 0, .false., rec=.false.)
        end if

    end subroutine s_amr_hdr_req

    !> Under the audit, check received header j against block k and the sender's site.
    impure subroutine s_amr_hdr_check(j, site, k)

        integer, intent(in) :: j, site, k

        if (XA_NH > 0) call s_xa_hdr_check(amr_fw_rq(XA_NH*(j - 1) + 1:XA_NH*j), site, k, [0, 0, 0], [0, 0, 0])

    end subroutine s_amr_hdr_check

    !> Post one zero-copy request per flagged face register of the current slot (block k): dir 1 sends to / 2 receives from peer,
    !! keyed by (block, direction, face).
    impure subroutine s_amr_freg_faces_req(dir, peer, site, k, f_lo, f_hi)

        integer, intent(in) :: dir, peer, site, k
        logical, intent(in) :: f_lo(3), f_hi(3)
        integer             :: cnt

        #:for D in [1, 2, 3]
            if (${D}$ <= num_dims) then
                cnt = size(freg(${D}$)%lo, 1)*size(freg(${D}$)%lo, 2)*size(freg(${D}$)%lo, 3)
                if (f_lo(${D}$)) call s_amr_wave_req_raw(amr_wave, dir, freg(${D}$)%lo(:,:,:,amr_reg_cur), cnt, peer, site, &
                    & k*8 + ${D}$*2)
                if (f_hi(${D}$)) call s_amr_wave_req_raw(amr_wave, dir, freg(${D}$)%hi(:,:,:,amr_reg_cur), cnt, peer, site, &
                    & k*8 + ${D}$*2 + 1)
            end if
        #:endfor

    end subroutine s_amr_freg_faces_req

    !> Move the flagged face registers of the current slot to the device (to_device) or to the host.
    impure subroutine s_amr_freg_faces_update(to_device, f_lo, f_hi)

        logical, intent(in) :: to_device, f_lo(3), f_hi(3)

        #:for D in [1, 2, 3]
            #:for S in ['lo', 'hi']
                if (${D}$ <= num_dims .and. f_${S}$(${D}$)) then
                    if (to_device) then
                        $:GPU_UPDATE(device='[freg(' + str(D) + ')%' + S + '(:, :, :, amr_reg_cur)]')
                    else
                        $:GPU_UPDATE(host='[freg(' + str(D) + ')%' + S + '(:, :, :, amr_reg_cur)]')
                    end if
                end if
            #:endfor
        #:endfor

    end subroutine s_amr_freg_faces_update

    !> Debug builds: poison the unflagged face registers of the current slot with quiet NaN before a receive lands, so a stray read
    !! of a face nobody ships NaNs within a step.
    impure subroutine s_amr_freg_faces_poison(f_lo, f_hi)

#ifdef MFC_DEBUG
        use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
#endif
        logical, intent(in) :: f_lo(3), f_hi(3)
#ifdef MFC_DEBUG
        real(wp) :: nanv

        if (amr_reg_cur <= 0) return
        nanv = ieee_value(0._wp, ieee_quiet_nan)
        #:for D in [1, 2, 3]
            #:for S in ['lo', 'hi']
                if (${D}$ <= num_dims .and. .not. f_${S}$(${D}$)) then
                    freg(${D}$)%${S}$(:,:,:,amr_reg_cur) = nanv
                    $:GPU_UPDATE(device='[freg(' + str(D) + ')%' + S + '(:, :, :, amr_reg_cur)]')
                end if
            #:endfor
        #:endfor
#endif

    end subroutine s_amr_freg_faces_poison

    !> Set the fine level's geometry (region, intersection, extents, bounds, coordinates) for the box lo:hi. Arrays are preallocated
    !! at max size; this only updates metadata and refills coords. Collective: all ranks must call together (init and regrid do); it
    !! also refreshes the allreduced amr_xchg_coarse_ghosts flag for the new box. Invariant: a level-l block's fine extent is
    !! amr_ref_ratio**l * (coarse-region width) - 1, not amr_ref_ratio*width. (amr_ref_ratio*width holds only for the level-1
    !! initial block; nested boxes compound by amr_ref_ratio per level.) Every fine-extent computation (here, the restart-reader
    !! check, the load-weight, the fmul) uses amr_ref_ratio**level; assuming amr_ref_ratio*width rejects level>=2 blocks as corrupt.
    impure subroutine s_set_amr_fine_geometry(lo, hi)

        integer, intent(in) :: lo(3), hi(3)
        integer             :: sidx(3), ext(3), fext(3), nmar, bad_loc, pblk

        amr_slots(amr_cur)%region%lo = lo; amr_slots(amr_cur)%region%hi = hi
        amr_region_lo = lo; amr_region_hi = hi  ! global mirror for m_amr_registers (no use-cycle)
        amr_region_lo_all(:,amr_cur) = lo; amr_region_hi_all(:,amr_cur) = hi

        ! a block is owned whole by amr_block_owner(k): the owner holds every fine cell, every other rank none. amr_isect_lo/hi is
        ! the block's coarse footprint on the owner (drives the coarse<->fine gather/scatter) and empty elsewhere.
        amr_rank_owns_block = (amr_block_owner(amr_cur) == proc_rank)
        pblk = 0
        if (amr_rank_owns_block) then
            amr_isect_lo = lo; amr_isect_hi = hi
            if (amr_block_level(amr_cur) >= 2) then
                ! multi-level: the coarse footprint in the parent block's fine frame (parent-fine index of L0 cell c is
                ! rr*(c - R1.lo)), the same footprint s_amr_parent_foot derives from replicated metadata. rr is the global
                ! amr_ref_ratio, not amr_slots(pblk)%amr_ref_ratio: a rank owning this block but not its parent never allocated
                ! pblk, so that field would read undefined.
                pblk = f_amr_parent_block(amr_cur)
                call s_amr_parent_foot(amr_cur, pblk, amr_isect_lo, amr_isect_hi)
            end if
        else
            amr_isect_lo = 1; amr_isect_hi = 0  ! empty footprint
        end if
        amr_isect_lo_all(:,amr_cur) = amr_isect_lo; amr_isect_hi_all(:,amr_cur) = amr_isect_hi
        amr_owns_all(amr_cur) = amr_rank_owns_block
        ! fine extents cover the whole block on the owner; -1 (empty) on non-owners
        fext = merge(amr_ref_ratio*max(amr_isect_hi - amr_isect_lo + 1, 0) - 1, 0, amr_dim)
        amr_slots(amr_cur)%m = fext(1); amr_slots(amr_cur)%n = fext(2); amr_slots(amr_cur)%p = fext(3)
        amr_slots(amr_cur)%idwbuff%beg = merge(-buff_size, 0, amr_dim); amr_slots(amr_cur)%idwbuff%end = merge(fext + buff_size, &
                  & 0, amr_dim)
        ! coord building only on ranks with fine cells (others never read their coord arrays). Every level builds the same way:
        ! replay the ancestor chain from the global L0 boundaries. The owner may hold no part of the coarse slice it refines, and
        ! (level>=2) may not own the parent at all, so neither the local coarse coords nor the parent's slot can be read here.
        if (amr_rank_owns_block) then
            #:for D, X in [(1, 'x'), (2, 'y'), (3, 'z')]
                if (amr_dim(${D}$)) call s_amr_build_block_coords(amr_cur, amr_g${X}$cb, amr_slots(amr_cur)%${X}$_cb, &
                    & amr_slots(amr_cur)%${X}$_cc, amr_slots(amr_cur)%d${X}$, ${D}$)
            #:endfor
        end if

        ! Fine ghost prolongation reads up to nmar coarse cells past each face of the intersection; if that stencil leaves any
        ! rank's interior (block near/at/across a rank boundary), the coarse cons ghosts it reads must be halo-exchanged before
        ! every fill (the solver populates only prim ghosts). All ranks agree on the flag, so the pairwise exchanges are called
        ! consistently.
        nmar = (buff_size + amr_ref_ratio - 1)/amr_ref_ratio + 1
        bad_loc = 0
        if (amr_rank_owns_block .and. any(amr_dim .and. (amr_isect_lo - sidx < nmar .or. sidx + ext - amr_isect_hi < nmar))) &
            & bad_loc = 1
        ! Accumulate, do not reduce: the caller closes the scan with s_amr_reduce_xchg_flag. Every caller loops over blocks and
        ! wants "does any block need the exchange" (the OR over blocks, not the last block's answer), in one collective rather
        ! than one per block.
        amr_xchg_bad = max(amr_xchg_bad, bad_loc)

    end subroutine s_set_amr_fine_geometry

    !> Close a geometry scan: one allreduce of the accumulated flag, then reset so the next scan starts clean. Must be called after
    !! every s_set_amr_fine_geometry loop (or single call); the flag it sets is read by the fine advance
    !! (s_amr_exchange_coarse_cons_halo).
    impure subroutine s_amr_reduce_xchg_flag()

        integer :: bad_glb

        call s_mpi_allreduce_integer_max(amr_xchg_bad, bad_glb)
        amr_xchg_coarse_ghosts = bad_glb == 1
        amr_xchg_bad = 0

    end subroutine s_amr_reduce_xchg_flag

    !> Conservative-linear prolongation for a single variable pair. Reads coarse interior/ghost from qc; writes fine interior to qf.
    !! Minmod-limited slopes.
    impure subroutine s_prolong_one_var(qc, loc, ivar)

        type(scalar_field), intent(in) :: qc
        integer, intent(in)            :: loc, ivar  !< flat-store slot and variable of the fine target
        integer                        :: fi, fj, fk, ci, cj, ck, ox, oy, oz, rrat, mm, nn, pp, il1, il2, il3
        real(wp)                       :: u0, sx, sy, sz, xix, xiy, xiz
        logical                        :: d2, d3

        ! coarse source qc is the gathered block-local patch amr_cg (fine-level distribution): amr_isect_lo is global and equals
        ! region_lo on the owner, so amr_isect_lo + f/rr - amr_cpat_off = nmar + f/rr is the patch-local coarse index.
        ! Device kernel: reads the patch's device mirror (pushed once per prolong dispatch by s_interpolate_coarse_to_fine) and
        ! writes the fine slot in place. CPU builds compile this to the identical plain loop.

        ox = amr_cpat_off(1); oy = amr_cpat_off(2); oz = amr_cpat_off(3)
        rrat = amr_slots(amr_cur)%amr_ref_ratio
        mm = amr_slots(amr_cur)%m; nn = amr_slots(amr_cur)%n; pp = amr_slots(amr_cur)%p
        il1 = amr_isect_lo(1); il2 = amr_isect_lo(2); il3 = amr_isect_lo(3)
        d2 = n_glb > 0; d3 = p_glb > 0
        $:GPU_PARALLEL_LOOP(collapse=3, private='[ci, cj, ck, xix, xiy, xiz, u0, sx, sy, sz]')
        do fk = 0, pp
            do fj = 0, nn
                do fi = 0, mm
                    ck = il3 + fk/rrat - oz; if (.not. d3) ck = 0
                    xiz = 0._wp; if (d3) xiz = (real(mod(fk, rrat), wp) - real(rrat - 1, wp)*0.5_wp)/real(rrat, wp)
                    cj = il2 + fj/rrat - oy; if (.not. d2) cj = 0
                    xiy = 0._wp; if (d2) xiy = (real(mod(fj, rrat), wp) - real(rrat - 1, wp)*0.5_wp)/real(rrat, wp)
                    ci = il1 + fi/rrat - ox
                    xix = (real(mod(fi, rrat), wp) - real(rrat - 1, wp)*0.5_wp)/real(rrat, wp)
                    u0 = real(qc%sf(ci, cj, ck), wp)
                    sx = minmod(real(qc%sf(ci + 1, cj, ck), wp) - u0, u0 - real(qc%sf(ci - 1, cj, ck), wp))
                    sy = 0._wp
                    if (d2) sy = minmod(real(qc%sf(ci, cj + 1, ck), wp) - u0, u0 - real(qc%sf(ci, cj - 1, ck), wp))
                    sz = 0._wp
                    if (d3) sz = minmod(real(qc%sf(ci, cj, ck + 1), wp) - u0, u0 - real(qc%sf(ci, cj, ck - 1), wp))
                    amr_cons_st(fi, fj, fk, ivar, loc) = u0 + sx*xix + sy*xiy + sz*xiz
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_prolong_one_var

    !> Conservative-linear prolongation: fill amr_fine interior from coarse (level-0), minmod-limited. Symmetric child offsets
    !! (+/-1/4 of a coarse cell) => the amr_ref_ratio^d children average to the coarse value. Multi-fluid volume fractions take the
    !! sum-preserving closure path instead (single-fluid runs never branch, so their prolongation is untouched).
    impure subroutine s_interpolate_coarse_to_fine()

        integer :: i

        ! the prolong kernels read the gathered patch on the device, where every fill wave assembles it

        do i = 1, sys_size
            ! Lagrangian bubbles: alphas sum to the local liquid fraction beta (not 1), so the sum-to-one closure would corrupt
            ! the EL state; each alpha prolongs plainly instead
            if (num_fluids > 1 .and. (.not. bubbles_lagrange) .and. i >= eqn_idx%adv%beg .and. i <= eqn_idx%adv%end) cycle
            if (chemistry .and. i >= eqn_idx%species%beg .and. i <= eqn_idx%species%end) cycle  ! sum/positivity closure below
            call s_prolong_one_var(amr_cg(i), amr_loc_of(amr_cur), i)
        end do
        if (num_fluids > 1 .and. (.not. bubbles_lagrange)) call s_prolong_alphas_closure(amr_cg, amr_loc_of(amr_cur))
        if (chemistry) call s_prolong_species_closure(amr_cg, amr_loc_of(amr_cur))

    end subroutine s_interpolate_coarse_to_fine

    !> Sum-preserving volume-fraction prolongation (num_fluids > 1): fluids adv%beg..adv%end-1 are interpolated with minmod slopes
    !! under a shared per-cell limiter switch (a sign change for any fluid in a dim zeroes that dim's slope for all fluids, so the
    !! closure fluid's effective slope is limited consistently) and clamped to [0,1]; the last fluid closes alpha_n = 1 -
    !! sum(others), so sum(alpha) = 1 on the fine level by construction. For two fluids the closure is also in [0,1]; for >2 fluids
    !! any residual closure undershoot is handled by mpp_lim (required by the checker). Same fine/coarse index mapping as
    !! s_prolong_one_var.
    impure subroutine s_prolong_alphas_closure(qc, loc)

        type(scalar_field), dimension(sys_size), intent(in) :: qc
        integer, intent(in)                                 :: loc
        integer                                             :: fi, fj, fk, ci, cj, ck, ox, oy, oz, i
        integer                                             :: rrat, mm, nn, pp, il1, il2, il3, advb, adve
        real(wp)                                            :: xix, xiy, xiz, u0, sx, sy, sz, av, asum
        logical                                             :: shx, shy, shz, d2, d3

        ! coarse source qc is the gathered block-local patch amr_cg (fine-level distribution): patch-frame offset.
        ! Device kernel; the shared limiter switch (s_alpha_shared_switch) is inlined verbatim.

        ox = amr_cpat_off(1); oy = amr_cpat_off(2); oz = amr_cpat_off(3)
        rrat = amr_slots(amr_cur)%amr_ref_ratio
        mm = amr_slots(amr_cur)%m; nn = amr_slots(amr_cur)%n; pp = amr_slots(amr_cur)%p
        il1 = amr_isect_lo(1); il2 = amr_isect_lo(2); il3 = amr_isect_lo(3)
        d2 = n_glb > 0; d3 = p_glb > 0
        advb = eqn_idx%adv%beg; adve = eqn_idx%adv%end
        $:GPU_PARALLEL_LOOP(collapse=3, private='[ci, cj, ck, xix, xiy, xiz, u0, sx, sy, sz, av, asum, shx, shy, shz, i]')
        do fk = 0, pp
            do fj = 0, nn
                do fi = 0, mm
                    ck = il3 + fk/rrat - oz; if (.not. d3) ck = 0
                    xiz = 0._wp; if (d3) xiz = (real(mod(fk, rrat), wp) - real(rrat - 1, wp)*0.5_wp)/real(rrat, wp)
                    cj = il2 + fj/rrat - oy; if (.not. d2) cj = 0
                    xiy = 0._wp; if (d2) xiy = (real(mod(fj, rrat), wp) - real(rrat - 1, wp)*0.5_wp)/real(rrat, wp)
                    ci = il1 + fi/rrat - ox
                    xix = (real(mod(fi, rrat), wp) - real(rrat - 1, wp)*0.5_wp)/real(rrat, wp)
                    ! shared per-cell limiter switch (inlined s_alpha_shared_switch): per dim, slopes stay on only if no
                    ! fluid's centered differences change sign there (symmetric in the fluids, incl. the closure fluid)
                    shx = .true.; shy = d2; shz = d3
                    do i = advb, adve
                        u0 = real(qc(i)%sf(ci, cj, ck), wp)
                        if ((real(qc(i)%sf(ci + 1, cj, ck), wp) - u0)*(u0 - real(qc(i)%sf(ci - 1, cj, ck), &
                            & wp)) <= 0._wp) shx = .false.
                        if (d2) then
                            if ((real(qc(i)%sf(ci, cj + 1, ck), wp) - u0)*(u0 - real(qc(i)%sf(ci, cj - 1, ck), &
                                & wp)) <= 0._wp) shy = .false.
                        end if
                        if (d3) then
                            if ((real(qc(i)%sf(ci, cj, ck + 1), wp) - u0)*(u0 - real(qc(i)%sf(ci, cj, ck - 1), &
                                & wp)) <= 0._wp) shz = .false.
                        end if
                    end do
                    asum = 0._wp
                    do i = advb, adve - 1
                        u0 = real(qc(i)%sf(ci, cj, ck), wp)
                        sx = 0._wp
                        if (shx) sx = minmod(real(qc(i)%sf(ci + 1, cj, ck), wp) - u0, u0 - real(qc(i)%sf(ci - 1, cj, ck), wp))
                        sy = 0._wp
                        if (d2 .and. shy) sy = minmod(real(qc(i)%sf(ci, cj + 1, ck), wp) - u0, u0 - real(qc(i)%sf(ci, cj - 1, &
                            & ck), wp))
                        sz = 0._wp
                        if (d3 .and. shz) sz = minmod(real(qc(i)%sf(ci, cj, ck + 1), wp) - u0, u0 - real(qc(i)%sf(ci, cj, &
                            & ck - 1), wp))
                        av = min(max(u0 + sx*xix + sy*xiy + sz*xiz, 0._wp), 1._wp)
                        amr_cons_st(fi, fj, fk, i, loc) = av
                        asum = asum + av
                    end do
                    amr_cons_st(fi, fj, fk, adve, loc) = 1._wp - asum
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_prolong_alphas_closure

    !> Species mass-fraction prolongation closure (chemistry): each partial density rho*Y_k is minmod-prolonged and clamped
    !! non-negative, then all species are rescaled so sum_k(rho*Y_k) equals the (already prolonged) continuity density at the fine
    !! cell. This keeps the fine species realizable (Y_k >= 0, and sum(Y_k) = 1 exactly under the cons->prim recovery rho = sum
    !! rho*Y_k) and consistent with the continuity variable the reaction source reads. Same index mapping as s_prolong_one_var; cont
    !! is prolonged in the main loop before this runs.
    impure subroutine s_prolong_species_closure(qc, loc)

        type(scalar_field), dimension(sys_size), intent(in) :: qc
        integer, intent(in)                                 :: loc
        integer                                             :: fi, fj, fk, ci, cj, ck, ox, oy, oz, i
        integer                                             :: rrat, mm, nn, pp, il1, il2, il3, spb, spe, cte
        real(wp)                                            :: xix, xiy, xiz, u0, sx, sy, sz, av, rsum, rscale
        logical                                             :: d2, d3

        ! coarse source qc is the gathered block-local patch amr_cg (fine-level distribution): patch-frame offset.
        ! Device kernel; the rescale re-reads only this thread's own cell, so the loop nest is safe.

        ox = amr_cpat_off(1); oy = amr_cpat_off(2); oz = amr_cpat_off(3)
        rrat = amr_slots(amr_cur)%amr_ref_ratio
        mm = amr_slots(amr_cur)%m; nn = amr_slots(amr_cur)%n; pp = amr_slots(amr_cur)%p
        il1 = amr_isect_lo(1); il2 = amr_isect_lo(2); il3 = amr_isect_lo(3)
        d2 = n_glb > 0; d3 = p_glb > 0
        spb = eqn_idx%species%beg; spe = eqn_idx%species%end; cte = eqn_idx%cont%end
        $:GPU_PARALLEL_LOOP(collapse=3, private='[ci, cj, ck, xix, xiy, xiz, u0, sx, sy, sz, av, rsum, rscale, i]')
        do fk = 0, pp
            do fj = 0, nn
                do fi = 0, mm
                    ck = il3 + fk/rrat - oz; if (.not. d3) ck = 0
                    xiz = 0._wp; if (d3) xiz = (real(mod(fk, rrat), wp) - real(rrat - 1, wp)*0.5_wp)/real(rrat, wp)
                    cj = il2 + fj/rrat - oy; if (.not. d2) cj = 0
                    xiy = 0._wp; if (d2) xiy = (real(mod(fj, rrat), wp) - real(rrat - 1, wp)*0.5_wp)/real(rrat, wp)
                    ci = il1 + fi/rrat - ox
                    xix = (real(mod(fi, rrat), wp) - real(rrat - 1, wp)*0.5_wp)/real(rrat, wp)
                    rsum = 0._wp
                    do i = spb, spe
                        u0 = real(qc(i)%sf(ci, cj, ck), wp)
                        sx = minmod(real(qc(i)%sf(ci + 1, cj, ck), wp) - u0, u0 - real(qc(i)%sf(ci - 1, cj, ck), wp))
                        sy = 0._wp
                        if (d2) sy = minmod(real(qc(i)%sf(ci, cj + 1, ck), wp) - u0, u0 - real(qc(i)%sf(ci, cj - 1, ck), wp))
                        sz = 0._wp
                        if (d3) sz = minmod(real(qc(i)%sf(ci, cj, ck + 1), wp) - u0, u0 - real(qc(i)%sf(ci, cj, ck - 1), wp))
                        av = max(u0 + sx*xix + sy*xiy + sz*xiz, 0._wp)
                        amr_cons_st(fi, fj, fk, i, loc) = av
                        rsum = rsum + av
                    end do
                    rscale = real(amr_cons_st(fi, fj, fk, cte, loc), wp)/max(rsum, 1.e-30_wp)
                    do i = spb, spe
                        amr_cons_st(fi, fj, fk, i, loc) = real(amr_cons_st(fi, fj, fk, i, loc), wp)*rscale
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_prolong_species_closure

    !> Disblock prolongation. Guard: no-op unless amr.
    impure subroutine s_populate_amr_fine(q_cons_base)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_base
        integer                                                :: islot, kk, i

        if (.not. amr) return
        ! Prolong every block (max_grid_size tiling can make several) from its coarse patch, assembled by the level-1 fill wave
        ! from the coarse owners' interiors (so no coarse-ghost halo is needed). This runs before s_initialize_gpu_vars, and the
        ! wave packs from the device, so push the ICs first. All ranks enter the wave; only owners consume and prolong.
        do i = 1, sys_size
            $:GPU_UPDATE(device='[q_cons_base(i)%sf]')
        end do
        call s_amr_l1_fill_exchange(q_cons_base, .true.)
        call s_amr_refresh_my_blocks()
        do kk = 1, amr_n_my
            islot = amr_my_blk(kk)
            if (amr_block_level(islot) /= 1) cycle
            call s_amr_select_slot(islot)
            call s_amr_l1_fill_consume(q_cons_base, islot, .true.)
            call s_interpolate_coarse_to_fine()  ! device kernel: writes the slot in place
        end do
        call s_amr_fill_wave_done()
        if (amr_max_level >= 2) call s_amr_build_static_multilevel()
        call s_amr_select_slot(f_l0_slot(1))

    end subroutine s_populate_amr_fine

    !> Build the static multi-level hierarchy (amr_regrid_int = 0): nest exactly one level-2 block inside level-1 block 1 by a fixed
    !! geometric inset (a regrid would place it by sensor-on-fine instead), prolong the parent state into it, and keep it persistent
    !! so the advance driver steps it every timestep. The restrict/reflux identity it relies on is protected by the static
    !! multi-level goldens and the runtime conservation-defect probe.
    impure subroutine s_amr_build_static_multilevel()

        integer :: L2, n1, par, inset(3)

        if (amr_max_level < 2) return
        n1 = amr_num_blocks
        if (n1 < 1) return
        ! the static hierarchy nests exactly one level-2 block; without pool room it would silently refine only to level 1 (an
        ! under-resolved but "successful" run). n1 (the level-1 tile count) is only known here, not at checker time, so abort at the
        ! point of failure. Replicated inputs -> every rank takes the same branch (collective-safe).
        if (n1 + 1 > l0_slot_off + amr_max_fine) call s_mpi_abort('amr static multi-level (amr_max_level > 1, ' &
            & // 'amr_regrid_int = 0): amr_max_blocks is too small to nest the level-2 block (need >= level-1 block count + 1); ' &
            & // 'increase amr_max_blocks')
        L2 = n1 + 1
        ! parent is the first fine block, f_l0_slot(1), not slot 1, which under coexist is the first level-0 tile. Insetting a
        ! tile instead would put the level-2 box in the wrong place and size it off the tile (a plausible-looking box that
        ! silently corrupts the run).
        par = f_l0_slot(1)
        inset = merge(max((amr_region_hi_all(:,par) - amr_region_lo_all(:,par) + 1)/4, amr_cpat_mar), 0, amr_dim)
        amr_region_lo_all(:,L2) = amr_region_lo_all(:,par) + inset
        amr_region_hi_all(:,L2) = amr_region_hi_all(:,par) - inset
        ! Guard the fixed-inset box against configs this single-block static builder cannot represent; the dynamic regrid path has
        ! the analogous checks (proper-nesting skip + amr_maxc_fit/2 clamp), but the static path bypasses them. Replicated inputs
        ! -> every rank takes the same branch (collective-safe). (a) a level-1 block smaller than 2*inset inverts the box; (b) a
        ! level-2 L0-extent > amr_maxc_fit/2 makes its parent-fine transverse extent (2*L0) overrun the creg register (allocated
        ! 0:amr_maxc_fit-1), a silent out-of-bounds device write in the L2->L1 reflux capture.
        if (any(amr_dim .and. amr_region_lo_all(:,L2) > amr_region_hi_all(:, &
            & L2))) call s_mpi_abort('amr static multi-level: level-1 block 1 is too small to nest a level-2 block (the fixed ' &
            & // 'inset inverts the box); enlarge the base amr block or reduce amr_cpat_mar')
        if (any(amr_dim .and. amr_ref_ratio*(amr_region_hi_all(:,L2) - amr_region_lo_all(:, &
            & L2) + 1) > amr_maxc_fit)) &
            & call s_mpi_abort('amr static multi-level: the nested level-2 block exceeds the per-rank scratch cap ' &
            & // '(2*L0-extent > amr_maxc_fit); static multi-level does not tile the level-2 block - use a smaller base amr ' &
            & // 'block or the dynamic regrid path (amr_regrid_int > 0)')
        amr_block_level(L2) = 2
        amr_block_owner(L2) = amr_block_owner(par); amr_myblk_dirty = .true.
        amr_num_blocks = L2; amr_num_levels = 2
        call s_amr_reconcile_slots()
        amr_cur = L2
        call s_set_amr_fine_geometry(amr_region_lo_all(:,L2), amr_region_hi_all(:,L2))
        call s_amr_reduce_xchg_flag()
        ! the parent-fill wave assembles the L2 patch from the parent block's fine array (all ranks enter; the owner consumes and
        ! prolongs in place on the device)
        call s_amr_parent_fill_exchange(2, .true.)
        call s_amr_select_slot(L2)
        if (amr_rank_owns_block) then
            call s_amr_parent_fill_consume(L2, .true.)
            call s_interpolate_coarse_to_fine()
        end if
        call s_amr_fill_wave_done()
        ! back to the first fine block: f_l0_slot(1), not slot 1, which under coexist is an L0 tile whose geometry would size
        ! s_initialize_weno_module's device tables off the wrong bounds
        call s_amr_select_slot(f_l0_slot(1))

    end subroutine s_amr_build_static_multilevel

    !> Restriction: each covered coarse cell = the average of its amr_ref_ratio^d fine children (equal weight: the grid is Cartesian
    !! and uniform, so children share a volume). Writes the caller's coarse target: in production the level-0 state q_cons_ts(1)%vf
    !! (the deliberate fold-back of fine data each step); init-time diagnostics pass a scratch buffer instead. Device kernel.
    impure subroutine s_restrict_fine_to_coarse(coarse_tgt)

        type(scalar_field), dimension(sys_size), intent(inout) :: coarse_tgt
        integer                                                :: rr, o1, o2, o3, owner, r, idx, boxsz, maxsz, nsrc, ierr
        integer                                                :: rlo(3), rhi(3), ilo(3), ihi(3), bl(3), bh(3)
        real(wp), allocatable                                  :: sbuf(:,:), rbuf(:)
        integer, allocatable                                   :: reqs(:)

        ! multi-level: a level>=2 block folds back into its parent block's fine array (the coarse side of level l is level l-1),
        ! not the L0 coarse_tgt. Same restriction kernel, targeted at the parent in the parent-fine frame. When child and parent
        ! sit on different ranks the fold is a P2P pair, so both participants must enter or the receiver never posts.

        if (amr_block_level(amr_cur) >= 2) then
            if (amr_rank_owns_block .or. amr_block_owner(f_amr_parent_block(amr_cur)) == proc_rank) then
                call s_amr_restrict_to_parent()
            end if
            return
        end if

        ! whole-block-per-rank fold-back: the block owner restricts its fine block to coarse averages over the covered cells
        ! [region_lo:region_hi] and scatters them point-to-point to the coarse-cell owners: the owner overwrites the covered cells
        ! it holds locally and sends each other coarse-owner exactly its covered slice (all sys_size in one message). Covered
        ! cells are in-domain (no ghosts), so each is owned by exactly one interior owner. At np=1 the owner owns every covered
        ! cell, sends nothing, and overwrites locally with the same child-sum.
        rr = amr_slots(amr_cur)%amr_ref_ratio
        call s_amr_region_box(amr_cur, rlo, rhi)
        owner = amr_block_owner(amr_cur)
        o1 = amr_sidx(1); o2 = amr_sidx(2); o3 = amr_sidx(3)
        maxsz = sys_size*product(rhi - rlo + 1)

        ! block set changed: rebuild the cached overlap-rank lists (same lazy trigger as s_amr_fine_fine_halo; local, replicated)
        if (amr_seam_pairs_dirty .or. amr_seam_pairs_nblk /= amr_num_blocks) call s_amr_build_seam_pairs()

        if (proc_rank == owner) then
            ! overwrite the covered cells this rank owns, then send each other coarse-owner its covered slice
            call s_amr_rank_interior(proc_rank, ilo, ihi)
            call s_amr_box_isect(rlo, rhi, ilo, ihi, bl, bh)
            ! owner-local covered cells (every covered cell at np=1): restrict fine(device) -> coarse(device) touching only those
            ! cells. Never push the whole coarse array back to the device here: that would clobber the device-advanced non-covered
            ! coarse cells with the stale host copy, a GPU-only divergence (invisible on CPU where host==device) that
            ! IGR/MHD/acoustic amplify.
            if (all(bl <= bh)) call s_amr_restrict_device_sf(coarse_tgt, amr_loc_of(amr_cur), bl, bh, rlo, rr, o1, o2, o3)
            ! cached destination list (every listed rank's interior overlaps the region by construction)
            nsrc = 0
            do idx = 1, amr_ovl_scatter_n(amr_cur)
                if (amr_ovl_scatter(idx, amr_cur) /= owner) nsrc = nsrc + 1
            end do
            if (nsrc > 0) then
                allocate (sbuf(maxsz, nsrc), reqs(nsrc))
                nsrc = 0
                do idx = 1, amr_ovl_scatter_n(amr_cur)
                    r = amr_ovl_scatter(idx, amr_cur)
                    if (r == owner) cycle
                    call s_amr_rank_interior(r, ilo, ihi)
                    call s_amr_box_isect(rlo, rhi, ilo, ihi, bl, bh)
                    nsrc = nsrc + 1
                    boxsz = sys_size*product(bh - bl + 1)
                    ! pack this destination's covered slice on the device: restrict averages straight into the wire buffer (same
                    ! child-sum order and wp values as the device overwrite above), with no full-field host pull
                    call s_amr_restrict_device_wire(sbuf(1:boxsz,nsrc), amr_loc_of(amr_cur), bl, bh, rlo, rr)
#ifdef MFC_MPI
                    call s_xa_rec(XA_F7A_SND, 1, boxsz, amr_cur)
                    call MPI_ISEND(sbuf(1, nsrc), boxsz, mpi_p, r, amr_cur, MPI_COMM_WORLD, reqs(nsrc), ierr)
#endif
                end do
#ifdef MFC_MPI
                call MPI_WAITALL(nsrc, reqs, MPI_STATUSES_IGNORE, ierr)
#endif
                deallocate (sbuf, reqs)
            end if
        else
            ! coarse-owner: if I hold covered cells, receive my slice from the owner and overwrite my local coarse
            call s_amr_rank_interior(proc_rank, ilo, ihi)
            call s_amr_box_isect(rlo, rhi, ilo, ihi, bl, bh)
            if (all(bl <= bh)) then
                boxsz = sys_size*product(bh - bl + 1)
                allocate (rbuf(boxsz))
#ifdef MFC_MPI
                call s_xa_rec(XA_F7A_RCV, 2, boxsz, amr_cur)
                call MPI_RECV(rbuf, boxsz, mpi_p, owner, amr_cur, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
#endif
                ! Device unpack of the covered box, writing only those cells (a whole-array push would clobber the device-advanced
                ! non-covered coarse cells with this rank's stale host copy). This must not be a host unpack followed by a
                ! strided GPU_UPDATE(device=) of the box: a non-contiguous 3-D array section in an OpenMP target update is copied
                ! by AMD flang as size(section) contiguous elements, so only the first row lands on the cells it names and the
                ! remainder silently overwrites neighbouring cells with stale host data (only at np >= 2 with a block whose owner
                ! holds none of its covered cells). The wire layout (ci fastest, then cj, ck, i) is exactly
                ! s_l0_pack_unpack_block's, so it unpacks s_amr_restrict_device_wire's buffer as-is.
                call s_l0_pack_unpack_block_sf(coarse_tgt, bl(1) - o1, bl(2) - o2, bl(3) - o3, bh(1) - bl(1), bh(2) - bl(2), &
                                               & bh(3) - bl(3), rbuf, .false.)
                deallocate (rbuf)
            end if
        end if

    end subroutine s_restrict_fine_to_coarse

    !> Multi-level restriction: fold the current level>=2 block's fine averages back into its parent block's fine array over the
    !! covered cells. Same child-sum kernel as the L0 fold-back, targeted at the parent in the parent-fine frame (amr_isect already
    !! parent-fine; offset 0 = the parent's local fine indexing). Local when child and parent share an owner; otherwise a P2P pair.
    impure subroutine s_amr_restrict_to_parent()

        integer               :: pblk, rr, cowner, powner, boxsz, ierr
        integer               :: plo(3), phi(3)
        real(wp), allocatable :: xbuf(:)

        pblk = f_amr_parent_block(amr_cur)
        cowner = amr_block_owner(amr_cur); powner = amr_block_owner(pblk)
        if (proc_rank /= cowner .and. proc_rank /= powner) return  ! not a participant

        ! Same replicated-metadata box as the gather, so the folding child and the receiving parent agree without a handshake.
        call s_amr_parent_foot(amr_cur, pblk, plo, phi)
        if (any(plo > phi)) return  ! empty footprint

        rr = amr_ref_ratio

        if (powner == cowner) then
            ! co-located (np=1, or a co-located tower): fold straight into the parent.
            call s_amr_restrict_device_st(amr_loc_of(pblk), amr_loc_of(amr_cur), plo, phi, plo, rr, 0, 0, 0)
            return
        end if

#ifdef MFC_MPI
        ! Split ownership: the child restricts locally and ships coarse cells (rr**num_dims fewer values than shipping its fine
        ! block), which restriction being an overwrite (not an accumulate) makes correct. Reuses the L0<->L1 scatter's pack and
        ! unpack; their wire layout (ci fastest, then cj, ck, i) is compatible.
        boxsz = sys_size*product(phi - plo + 1)
        allocate (xbuf(boxsz))
        if (proc_rank == cowner) then
            call s_amr_restrict_device_wire(xbuf, amr_loc_of(amr_cur), plo, phi, plo, rr)
            call s_xa_rec(XA_F7B_SND, 1, boxsz, amr_cur)
            call MPI_SEND(xbuf, boxsz, mpi_p, powner, amr_cur, MPI_COMM_WORLD, ierr)
        else
            call s_xa_rec(XA_F7B_RCV, 2, boxsz, amr_cur)
            call MPI_RECV(xbuf, boxsz, mpi_p, cowner, amr_cur, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
            ! Device unpack of just the covered box; never a host unpack plus a strided GPU_UPDATE (see the L0 scatter's note: AMD
            ! flang copies a non-contiguous 3-D section as contiguous elements and silently corrupts neighbouring cells).
            call s_l0_pack_unpack_block_st(amr_loc_of(pblk), plo(1), plo(2), plo(3), phi(1) - plo(1), phi(2) - plo(2), &
                                           & phi(3) - plo(3), xbuf, .false.)
        end if
        deallocate (xbuf)
#endif

    end subroutine s_amr_restrict_to_parent

    !> The lock-step restrict fold as plan-based waves, on the lock-step np>1 path: per level (finest first) the split-ownership
    !! child->parent folds are one aggregated message per (child-owner, parent-owner) peer (s_amr_restrict_parent_wave), followed by
    !! the per-box reflux-to-parent applies (their freg was already exchanged by s_amr_freg_wave); then the level-1 -> L0
    !! covered-cell scatter is one aggregated message per (owner, coarse-owner) peer (s_amr_restrict_l1_wave). Both sides derive
    !! identical transfer lists from replicated metadata (the region-box x rank-interior intersections), so the wire layout needs no
    !! handshake and the F7 family words are exact. Level order (finest first) preserves child-before-parent folding; within a level
    !! the covered targets are disjoint, and the restrict/reflux interleave is order-free because sibling-shared faces carry weight
    !! 0 (s_amr_sibling_face_weights). np=1 uses the per-box loop.
    impure subroutine s_amr_restrict_wave(coarse_tgt, dt_reflux)

        type(scalar_field), dimension(sys_size), intent(inout) :: coarse_tgt
        real(wp), intent(in)                                   :: dt_reflux
        integer                                                :: lev, k, io, ifc, ko, kf

        call s_amr_refresh_lists()
        do lev = amr_max_level, 2, -1
            call s_phase_tic(PH_RESTR)
            call s_amr_restrict_parent_wave(lev)
            call s_phase_toc(PH_RESTR)
            ! s_amr_reflux_to_parent returns unless `own_child .or. own_parent`: own_child = amr_rank_owns_block, the
            ! multi-owner amr_owns_all notion (amr_own_blk, not amr_my_blk); own_parent = the level-lev children of my parents,
            ! i.e. amr_fch_blk plus the ones I own myself, already in amr_own_blk. The two lists overlap, and a duplicate visit
            ! would reflux twice, so walk their union: both lists are ascending, and a two-cursor merge from the tails visits each
            ! block once, in descending block order (the parent registers accumulate in that order).
            io = amr_n_own; ifc = amr_n_fch
            do
                do while (io >= 1)
                    if (amr_block_level(amr_own_blk(io)) == lev) exit
                    io = io - 1
                end do
                do while (ifc >= 1)
                    if (amr_block_level(amr_fch_blk(ifc)) == lev) exit
                    ifc = ifc - 1
                end do
                ko = 0; kf = 0
                if (io >= 1) ko = amr_own_blk(io)
                if (ifc >= 1) kf = amr_fch_blk(ifc)
                k = max(ko, kf)
                if (k == 0) exit
                if (ko == k) io = io - 1
                if (kf == k) ifc = ifc - 1
                call s_amr_select_slot(k)
                call s_phase_tic(PH_RESTR)
                call s_amr_reflux_to_parent(dt_reflux)
                call s_phase_toc(PH_RESTR)
            end do
        end do
        call s_phase_tic(PH_RESTR)
        call s_amr_restrict_l1_wave(coarse_tgt)
        call s_phase_toc(PH_RESTR)
        call s_amr_select_slot(1)

    end subroutine s_amr_restrict_wave

    !> One per-level wave of the split-ownership level>=2 child->parent restrict folds (the F7B pairs). Co-located folds run inline
    !! on the child-owner (bit-for-bit the per-box kernel); every cross-rank fold ships its parent-frame covered box in one
    !! aggregated message per (child-owner, parent-owner) peer. Both sides walk the same replicated block list ascending with
    !! per-peer running offsets, so the wire layout agrees with no metadata exchange.
    impure subroutine s_amr_restrict_parent_wave(lev)

        integer, intent(in) :: lev

#ifdef MFC_MPI
        integer :: k, pblk, cowner, powner, rr, idx, lo, hi, kk
        integer :: plo(3), phi(3), bl(3), bh(3)

        rr = amr_ref_ratio
        call s_amr_wave_open(amr_wave, 7)
        ! send side: every owned level-lev block whose parent lives elsewhere ships its whole footprint to the parent's owner;
        ! a co-located parent is folded in place
        call s_amr_wave_reset(amr_wsend)
        call s_amr_refresh_my_blocks()
        call s_amr_refresh_lists()
        do kk = 1, amr_n_my
            k = amr_my_blk(kk)
            if (amr_block_level(k) /= lev) cycle
            cowner = amr_block_owner(k)
            pblk = amr_parent_blk(k)
            powner = amr_block_owner(pblk)
            call s_amr_parent_foot(k, pblk, plo, phi)
            if (any(plo > phi)) cycle
            if (cowner == powner) then
                call s_amr_restrict_device_st(amr_loc_of(pblk), amr_loc_of(k), plo, phi, plo, rr, 0, 0, 0)
                cycle
            end if
            call s_amr_wave_add(amr_wsend, powner, k, plo, phi, sys_size*product(phi - plo + 1))
        end do
        call s_amr_wave_close(amr_wsend, amr_fw_sq, amr_fw_dev)
        ! receive side: the level-lev children of my parents that another rank owns
        call s_amr_wave_reset(amr_wrecv)
        do kk = 1, amr_n_fch
            k = amr_fch_blk(kk)
            if (amr_block_level(k) /= lev) cycle
            pblk = amr_parent_blk(k)
            cowner = amr_block_owner(k); powner = amr_block_owner(pblk)
            if (cowner == powner .or. proc_rank /= powner) cycle
            call s_amr_parent_foot(k, pblk, plo, phi)
            if (any(plo > phi)) cycle
            call s_amr_wave_add(amr_wrecv, cowner, k, plo, phi, sys_size*product(phi - plo + 1))
        end do
        call s_amr_wave_close(amr_wrecv, amr_fw_rq, amr_fw_dev)
        if (amr_wsend%np + amr_wrecv%np == 0) return
        call s_amr_wave_post(amr_wave, amr_wrecv, amr_fw_rq, XA_F7BW_RCV, amr_fw_dev)
        do idx = 1, amr_wsend%nx
            call s_amr_wave_slice(amr_wsend, idx, lo, hi)
            call s_amr_restrict_device_wire(amr_fw_sq(lo:hi), amr_loc_of(amr_wsend%blk(idx)), amr_wsend%bl(:,idx), amr_wsend%bh(:, &
                                            & idx), amr_wsend%bl(:,idx), rr)
            call s_amr_wave_hdr_pack(amr_wsend, amr_fw_sq, idx, XA_F7BW_SND)
        end do
        call s_amr_wave_send(amr_wave, amr_wsend, amr_fw_sq, XA_F7BW_SND, amr_fw_dev)
        call s_amr_wave_wait(amr_wave)
        do idx = 1, amr_wrecv%nx
            call s_amr_wave_hdr_check(amr_wrecv, amr_fw_rq, idx, XA_F7BW_SND)
            call s_amr_wave_slice(amr_wrecv, idx, lo, hi)
            bl = amr_wrecv%bl(:,idx); bh = amr_wrecv%bh(:,idx)
            call s_l0_pack_unpack_block_st(amr_loc_of(amr_parent_blk(amr_wrecv%blk(idx))), bl(1), bl(2), bl(3), bh(1) - bl(1), &
                                           & bh(2) - bl(2), bh(3) - bl(3), amr_fw_rq(lo:hi), .false.)
        end do
#endif

    end subroutine s_amr_restrict_parent_wave

    !> The level-1 -> L0 covered-cell scatter (F7A) as one wave: every owned level-1 block's covered slabs for every listed
    !! coarse-owner ship in one aggregated message per peer; the owner-local covered overwrite stays grouped per block during the
    !! pack walk. The receiver plan is my-interior x region(k) over the level-1 blocks I do not own, by construction
    !! (s_amr_ranks_overlapping) exactly the sender's list membership.
    impure subroutine s_amr_restrict_l1_wave(coarse_tgt)

        type(scalar_field), dimension(sys_size), intent(inout) :: coarse_tgt

#ifdef MFC_MPI
        integer :: k, owner, rr, idx, r, lo, hi, o1, o2, o3, cur, kk
        integer :: rlo(3), rhi(3), ilo(3), ihi(3), milo(3), mihi(3), bl(3), bh(3)

        call s_amr_wave_open(amr_wave, 6)
        o1 = amr_sidx(1); o2 = amr_sidx(2); o3 = amr_sidx(3)
        ! block set changed: rebuild the cached overlap-rank lists (same lazy trigger as s_amr_fine_fine_halo; local, replicated)
        if (amr_seam_pairs_dirty .or. amr_seam_pairs_nblk /= amr_num_blocks) call s_amr_build_seam_pairs()
        call s_amr_rank_interior(proc_rank, milo, mihi)
        ! send side: for every owned level-1 block, its covered slab inside each listed coarse-owner's interior
        call s_amr_wave_reset(amr_wsend)
        call s_amr_refresh_my_blocks()
        do kk = 1, amr_n_my
            k = amr_my_blk(kk)
            if (amr_block_level(k) /= 1) cycle
            call s_amr_region_box(k, rlo, rhi)
            do idx = 1, amr_ovl_scatter_n(k)
                r = amr_ovl_scatter(idx, k)
                if (r == proc_rank) cycle
                call s_amr_rank_interior(r, ilo, ihi)
                call s_amr_box_isect(rlo, rhi, ilo, ihi, bl, bh)
                if (any(bl > bh)) cycle
                call s_amr_wave_add(amr_wsend, r, k, bl, bh, sys_size*product(bh - bl + 1))
            end do
        end do
        call s_amr_wave_close(amr_wsend, amr_fw_sq, amr_fw_dev)
        ! receive side: every level-1 block I do not own that covers my interior, from its owner (the cached l1r list is
        ! by construction exactly the senders' membership)
        call s_amr_wave_reset(amr_wrecv)
        call s_amr_refresh_lists()
        do kk = 1, amr_n_l1r
            k = amr_l1r_blk(kk)
            owner = amr_block_owner(k)
            call s_amr_region_box(k, rlo, rhi)
            call s_amr_box_isect(rlo, rhi, milo, mihi, bl, bh)
            call s_amr_wave_add(amr_wrecv, owner, k, bl, bh, sys_size*product(bh - bl + 1))
        end do
        call s_amr_wave_close(amr_wrecv, amr_fw_rq, amr_fw_dev)
        call s_amr_wave_post(amr_wave, amr_wrecv, amr_fw_rq, XA_F7W_RCV, amr_fw_dev)
        ! owner-local covered overwrites + device packs, grouped per owned block: the transfer list is k-grouped by
        ! construction, so a monotone cursor drains each block's sends inside its group
        cur = 1
        call s_amr_refresh_my_blocks()
        do kk = 1, amr_n_my
            k = amr_my_blk(kk)
            if (amr_block_level(k) /= 1) cycle
            rr = amr_slots(k)%amr_ref_ratio
            call s_amr_region_box(k, rlo, rhi)
            call s_amr_box_isect(rlo, rhi, milo, mihi, bl, bh)
            if (all(bl <= bh)) call s_amr_restrict_device_sf(coarse_tgt, amr_loc_of(k), bl, bh, rlo, rr, o1, o2, o3)
            do while (cur <= amr_wsend%nx)
                if (amr_wsend%blk(cur) /= k) exit
                call s_amr_wave_slice(amr_wsend, cur, lo, hi)
                call s_amr_restrict_device_wire(amr_fw_sq(lo:hi), amr_loc_of(k), amr_wsend%bl(:,cur), amr_wsend%bh(:,cur), rlo, rr)
                call s_amr_wave_hdr_pack(amr_wsend, amr_fw_sq, cur, XA_F7W_SND)
                cur = cur + 1
            end do
        end do
        call s_amr_wave_send(amr_wave, amr_wsend, amr_fw_sq, XA_F7W_SND, amr_fw_dev)
        call s_amr_wave_wait(amr_wave)
        do idx = 1, amr_wrecv%nx
            call s_amr_wave_hdr_check(amr_wrecv, amr_fw_rq, idx, XA_F7W_SND)
            call s_amr_wave_slice(amr_wrecv, idx, lo, hi)
            bl = amr_wrecv%bl(:,idx); bh = amr_wrecv%bh(:,idx)
            ! Device unpack of the covered box only (the strided-update flang trap; see s_restrict_fine_to_coarse)
            call s_l0_pack_unpack_block_sf(coarse_tgt, bl(1) - o1, bl(2) - o2, bl(3) - o3, bh(1) - bl(1), bh(2) - bl(2), &
                                           & bh(3) - bl(3), amr_fw_rq(lo:hi), .false.)
        end do
#endif

    end subroutine s_amr_restrict_l1_wave

    !> Sibling-seam face weights for level>=2 block kb under parent pblk: 0 on a face shared with a same-parent sibling tile
    !! (fine-fine, not c/f: refluxing there double-writes and leaks; the outside parent cell is covered by the sibling's restrict),
    !! 1 otherwise. Replicated metadata only (f_amr_parent_block + f_amr_seam read amr_region_*_all), so every rank derives the same
    !! weights; the reflux apply and the freg wave's wire-skip both call this, so what ships and what is consumed cannot drift
    !! apart.
    impure subroutine s_amr_sibling_face_weights(kb, pblk, w_lo, w_hi)

        integer, intent(in)   :: kb, pblk
        real(wp), intent(out) :: w_lo(3), w_hi(3)
        integer               :: c, y, d

        ! iterate pblk's cached children (same-parent guarantees same level) instead of scanning every block with an
        ! O(global blocks) parent lookup per candidate. Every caller sits in a routine that already refreshed the epoch-keyed
        ! lists; refreshing here would let a loop body reallocate the very list its caller iterates.

        w_lo = 1._wp; w_hi = 1._wp
        if (pblk <= 0) return  ! orphan parent: no siblings (all weights 1); the CSR would index ptr(-1)
        do c = amr_child_ptr(pblk - 1) + 1, amr_child_ptr(pblk)
            y = amr_child_idx(c)
            if (y == kb) cycle
            do d = 1, num_dims
                if (f_amr_seam(kb, y, d)) w_hi(d) = 0._wp  ! sibling just above -> shared high face
                if (f_amr_seam(y, kb, d)) w_lo(d) = 0._wp  ! sibling just below -> shared low face
            end do
        end do

    end subroutine s_amr_sibling_face_weights

    !> Multi-level reflux: apply the Berger-Colella C/F flux correction from the current level>=2 block into its parent block's
    !! cells just outside the block footprint, in the parent-fine frame (mirror of the L0 s_amr_apply_reflux targeted at the parent;
    !! "the coarse" is level l-1). State form: q_parent(outside) += dt*(F_coarse - Fbar_fine)/dxf on the low face and +=
    !! dt*(Fbar_fine - F_coarse)/dxf on the high face, where Fbar_fine is the child-averaged fine register. creg/freg key off this
    !! block's slot. Per-face parent-fine dx (stretched-grid safe).
    !!
    !! The parent's owner applies: it holds the parent field and the parent-side creg (captured over its own advance). Only freg
    !! crosses the wire, and only when the two owners differ. Both participants must reach this routine or the P2P pair deadlocks
    !! (cf. the restrict).
    impure subroutine s_amr_reflux_to_parent(dt_reflux)

        real(wp), intent(in) :: dt_reflux
        integer              :: pblk, d, olo(3), ohi(3), glo(3), ghi(3), woff(3), plo(3), phi(3)
        real(wp)             :: w_lo(3), w_hi(3), mlo(3), mhi(3)
        logical              :: own_child, own_parent

        call s_amr_refresh_lists()  ! cached parent (f_amr_parent_block is an O(global blocks) scan; this runs per block)
        pblk = amr_parent_blk(amr_cur)
        own_child = amr_rank_owns_block
        own_parent = (amr_block_owner(pblk) == proc_rank)
        if (.not. (own_child .or. own_parent)) return
        if (.not. own_parent) return
        ! max_grid_size tiling of a level>=2 feature: a face shared with an adjacent sibling tile (same parent) is fine-fine, not a
        ! c/f boundary; its "outside" parent cell is covered by the sibling's restrict, so refluxing there double-writes and leaks.
        ! Skip those faces (weight 0); the fine-fine halo already matched the shared seam flux. No siblings -> all weights 1
        ! (no-op).
        call s_amr_sibling_face_weights(amr_cur, pblk, w_lo, w_hi)
        ! parent-fine frame for the shared reflux kernel: outside cell = isect boundary +/-1, creg-local loop range 0:extent,
        ! transverse write at the isect origin, per-face parent-fine dx. Footprint from replicated metadata (s_amr_parent_foot),
        ! not amr_isect_lo/hi, which is the empty non-owner sentinel on the parent's owner; rr likewise from the global
        ! amr_ref_ratio, since amr_slots(amr_cur) need not be allocated on this rank.
        call s_amr_parent_foot(amr_cur, pblk, plo, phi)
        olo = 0; ohi = 0; glo = 0; ghi = 0; woff = 0; mlo = 1._wp; mhi = 1._wp
        do d = 1, num_dims
            olo(d) = plo(d) - 1; ohi(d) = phi(d) + 1
            ghi(d) = phi(d) - plo(d)
            woff(d) = plo(d)
        end do
        #:for D, X in [(1, 'x'), (2, 'y'), (3, 'z')]
            if (amr_dim(${D}$)) then; mlo(${D}$) = amr_slots(pblk)%d${X}$(olo(${D}$)); mhi(${D}$) &
                & = amr_slots(pblk)%d${X}$(ohi(${D}$)); end if
        #:endfor
        call s_amr_br_load_faces(amr_loc_of(pblk), olo, ohi, glo, ghi, woff, w_lo, w_hi)
        call s_amr_reflux_apply_faces(amr_cons_br, amr_reg_cur, amr_ref_ratio, dt_reflux, olo, ohi, glo, ghi, woff, w_lo, w_hi, &
                                      & mlo, mhi)
        call s_amr_br_store_faces(amr_loc_of(pblk), olo, ohi, glo, ghi, woff, w_lo, w_hi)

    end subroutine s_amr_reflux_to_parent

    !> Device-native restriction: restrict the fine block (device, the flat store) to coarse averages over the covered coarse cells
    !! [bl:bh] global (block region origin rlo, ratio rr), touching only those cells. Child-sum order: ddk, ddj, then ddi; /nchild.
    !! Three destinations, one body, so owner-local, parent-folded and scattered coarse cells match bit-for-bit: the level-0
    !! monolithic field (`_sf`, local origin o), a parent block in the flat store (`_st`), both stp-cast, or the contiguous wire
    !! buffer buf (`_wire`, host via copyout; packed ci fastest, then cj, ck, i, s_l0_pack_unpack_block's layout, in wp since the
    !! receiver casts). No whole-coarse device push: that would clobber the device-advanced non-covered coarse cells.
    #:for SFX in ['sf', 'st', 'wire']
        #:set OARGS = '' if SFX == 'wire' else ', o1, o2, o3'
        impure subroutine s_amr_restrict_device_${SFX}$(${ {'sf': 'coarse_tgt', 'st': 'ctloc', 'wire': 'buf'}[SFX] }$, loc, bl, &
            & bh, rlo, rr${OARGS}$)

            #:if SFX == 'sf'
                type(scalar_field), dimension(sys_size), intent(inout) :: coarse_tgt
            #:elif SFX == 'st'
                integer, intent(in) :: ctloc
            #:else
                real(wp), intent(inout), contiguous :: buf(:)
            #:endif
            integer, intent(in) :: loc, bl(3), bh(3), rlo(3), rr${OARGS}$
            integer             :: i, ci, cj, ck, fi0, fj0, fk0, ddi, ddj, ddk, dj_hi, dk_hi, nchild
            integer             :: bl1, bl2, bl3, bh1, bh2, bh3, rl1, rl2, rl3, n1, n2, n3
            real(wp)            :: acc

            bl1 = bl(1); bl2 = bl(2); bl3 = bl(3); bh1 = bh(1); bh2 = bh(2); bh3 = bh(3)
            rl1 = rlo(1); rl2 = rlo(2); rl3 = rlo(3)
            n1 = bh1 - bl1 + 1; n2 = bh2 - bl2 + 1; n3 = bh3 - bl3 + 1
            dj_hi = merge(rr - 1, 0, amr_dim(2)); dk_hi = merge(rr - 1, 0, amr_dim(3)); nchild = rr**num_dims
            #:if SFX == 'wire'
                $:GPU_PARALLEL_LOOP(collapse=4, private='[fi0, fj0, fk0, ddi, ddj, ddk, acc]', copyout='[buf]')
            #:else
                $:GPU_PARALLEL_LOOP(collapse=4, private='[fi0, fj0, fk0, ddi, ddj, ddk, acc]')
            #:endif
            do i = 1, sys_size
                do ck = bl3, bh3
                    do cj = bl2, bh2
                        do ci = bl1, bh1
                            fi0 = (ci - rl1)*rr; fj0 = (cj - rl2)*rr; fk0 = (ck - rl3)*rr
                            acc = 0._wp
                            do ddk = 0, dk_hi
                                do ddj = 0, dj_hi
                                    do ddi = 0, rr - 1
                                        acc = acc + real(amr_cons_st(fi0 + ddi, fj0 + ddj, fk0 + ddk, i, loc), wp)
                                    end do
                                end do
                            end do
                            #:if SFX == 'sf'
                                coarse_tgt(i)%sf(ci - o1, cj - o2, ck - o3) = real(acc/real(nchild, wp), stp)
                            #:elif SFX == 'st'
                                amr_cons_st(ci - o1, cj - o2, ck - o3, i, ctloc) = real(acc/real(nchild, wp), stp)
                            #:else
                                buf(1 + (ci - bl1) + n1*((cj - bl2) + n2*((ck - bl3) + n3*(i - 1)))) = acc/real(nchild, wp)
                            #:endif
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

        end subroutine s_amr_restrict_device_${SFX}$
    #:endfor
end module m_amr_transfer
