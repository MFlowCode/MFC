!>
!!@file
!!@brief Contains module m_amr_transfer

#! AMD OpenMP: assert allocatables present on every kernel (see OMP_DEFAULT_STR); a kernel may name only module arrays allocated
#! under its own launch condition.
#:set MFC_OMP_PRESENT_ALLOCATABLE = True
#:include 'macros.fpp'

!> @brief Prolongation, restriction and Berger-Colella reflux, including the restrict/reflux/freg waves.
module m_amr_transfer

#ifdef MFC_MPI
    use mpi
#endif

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy, only: s_mpi_abort
    use m_mpi_common, only: s_mpi_allreduce_integer_max
    use m_amr_registers, only: s_amr_reflux_apply_faces, s_amr_parent_foot, freg, s_amr_reg_prepare, s_amr_reflux_faces, &
        & f_amr_block_near
    use m_phase_timing
    use m_amr_state
    use m_amr_wave
    use m_amr_distribution
    use m_amr_store
    use m_amr_exchange
    use m_amr_ghost, only: s_amr_slab_add
    use m_amr_seam
    use m_amr_frame

    implicit none

    private
    public :: s_amr_freg_wave, s_amr_reduce_xchg_flag, s_amr_reflux_faces_wave, s_amr_reflux_to_parent, s_amr_restrict_wave, &
        & s_amr_sync_level, s_interpolate_coarse_to_fine, s_populate_amr_fine, s_set_amr_fine_geometry

contains

    !> True iff rank r refluxes some face of the current block (its unclipped subdomain); the participation map and both exchange
    !! paths gate freg receives on it.
    pure logical function f_amr_reflux_participates(r) result(part)

        integer, intent(in) :: r
        integer             :: sidx(3), ext(3)
        logical             :: s_lo(3), s_hi(3)

        call s_amr_rank_decomp(r, sidx, ext)
        call s_amr_reflux_faces(sidx, ext, .false., s_lo, s_hi)
        part = any(s_lo .or. s_hi)

    end function f_amr_reflux_participates

    !> Faces of the current block that rank r applies (seam-clipped); sender and receiver derive the identical set from replicated
    !! data.
    pure subroutine s_amr_reflux_faces_for(r, s_lo, s_hi)

        integer, intent(in)  :: r
        logical, intent(out) :: s_lo(3), s_hi(3)
        integer              :: sidx(3), ext(3)

        call s_amr_rank_decomp(r, sidx, ext)
        call s_amr_reflux_faces(sidx, ext, .true., s_lo, s_hi)

    end subroutine s_amr_reflux_faces_for

    !> Per-stage level-1 reflux-face exchange as one wave: faces are packed on the owner's device into one pool, cross as one
    !! message per peer and unpack on the device. Both ends list a peer's faces by ascending block, dim, low before high, from the
    !! same predicates, so layouts agree with no handshake. Nothing may post into freg before s_amr_reg_prepare (the apply can
    !! reallocate the registers).
    impure subroutine s_amr_reflux_faces_wave()

#ifdef MFC_MPI
        integer :: k, r, idx, ncand, kk2, t
        integer :: cand(num_procs), glo(3), ghi(3), msidx(3), mext(3)
        logical :: s_lo(3), s_hi(3)

        if (num_procs == 1) return
        call s_amr_reg_prepare()
        call s_amr_wave_open(amr_wave, amr_band_reflux)
        call s_amr_wave_reset(amr_wsend); call s_amr_wave_reset(amr_wrecv)
        call s_amr_refresh_lists()
        call s_amr_rank_decomp(proc_rank, msidx, mext)
        do kk2 = 1, amr_n_l1p
            k = amr_l1p_blk(kk2)
            if (amr_block_owner(k) == proc_rank) cycle
            if (.not. f_amr_block_near(k, msidx, mext)) cycle
            call s_amr_select_slot(k)
            if (.not. f_amr_reflux_participates(proc_rank)) cycle
            call s_amr_reflux_faces_for(proc_rank, s_lo, s_hi)
            call s_amr_face_xfers(amr_wrecv, amr_block_owner(k), k, s_lo, s_hi)
        end do
        call s_amr_refresh_my_blocks()
        do kk2 = 1, amr_n_my
            k = amr_my_blk(kk2)
            if (amr_block_level(k) /= 1 .or. amr_block_owner(k) /= proc_rank) cycle
            call s_amr_select_slot(k)
            glo = merge(amr_region_lo - 1, 0, amr_dim); ghi = merge(amr_region_hi + 1, 0, amr_dim)
            call s_amr_ranks_overlapping(glo, ghi, cand, ncand)
            do idx = 1, ncand
                r = cand(idx)
                if (r == proc_rank .or. .not. f_amr_reflux_participates(r)) cycle
                call s_amr_reflux_faces_for(r, s_lo, s_hi)
                call s_amr_face_xfers(amr_wsend, r, k, s_lo, s_hi)
            end do
        end do
        call s_amr_wave_close(amr_wsend, amr_rf_sq, .true.)
        call s_amr_wave_close(amr_wrecv, amr_rf_rq, .true.)
        call s_amr_wave_post(amr_wave, amr_wrecv, amr_rf_rq, .false.)
        if (amr_wsend%nx > 0) then
            call s_amr_freg_pack(amr_wsend, .true.)
            $:GPU_UPDATE(host='[amr_rf_sq(1:amr_wsend%words)]')
        end if
        call s_amr_wave_send(amr_wave, amr_wsend, amr_rf_sq, .false.)
        call s_amr_wave_wait(amr_wave)
        if (amr_wrecv%nx > 0) then
            $:GPU_UPDATE(device='[amr_rf_rq(1:amr_wrecv%words)]')
            call s_amr_freg_pack(amr_wrecv, .false.)
        end if
        call s_amr_select_slot(1)
#endif

    end subroutine s_amr_reflux_faces_wave

    !> Append block k's faces f_lo/f_hi with peer as transfers (bl = dim, side 1 low / 2 high, register slot).
    impure subroutine s_amr_face_xfers(w, peer, k, f_lo, f_hi)

        type(t_amr_wave_side), intent(inout) :: w
        integer, intent(in)                  :: peer, k
        logical, intent(in)                  :: f_lo(3), f_hi(3)
        integer                              :: d, nw

        do d = 1, num_dims
            nw = size(freg(d)%lo, 1)*size(freg(d)%lo, 2)*size(freg(d)%lo, 3)
            if (f_lo(d)) call s_amr_wave_add(w, peer, k, [d, 1, amr_reg_of(k)], [0, 0, 0], nw)
            if (f_hi(d)) call s_amr_wave_add(w, peer, k, [d, 2, amr_reg_of(k)], [0, 0, 0], nw)
        end do

    end subroutine s_amr_face_xfers

    !> Copy a side's freg faces into (pack) or out of the device pool, one kernel per face direction. CCE 19 packs garbage when one
    !! procedure holds several of these loop nests (as in s_amr_apply_reflux), so each direction has its own procedure.
    impure subroutine s_amr_freg_pack(w, pack)

        type(t_amr_wave_side), intent(in) :: w
        logical, intent(in)               :: pack

        #:for D in [1, 2, 3]
            if (${D}$ <= num_dims) call s_amr_freg_pack_d${D}$(w, pack)
        #:endfor

    end subroutine s_amr_freg_pack

    #:for D in [1, 2, 3]
        !> s_amr_freg_pack for the faces normal to direction ${D}$, both sides in one kernel.
        impure subroutine s_amr_freg_pack_d${D}$(w, pack)

            type(t_amr_wave_side), intent(in) :: w
            logical, intent(in)               :: pack
            integer                           :: tsl(w%nx), toff(w%nx), tsd(w%nx)
            integer                           :: t, nt, n1, n2, eq, a, b, j

            nt = 0
            do t = 1, w%nx
                if (w%bl(1, t) /= ${D}$) cycle
                nt = nt + 1; tsl(nt) = w%bl(3, t); toff(nt) = w%off(t); tsd(nt) = w%bl(2, t)
            end do
            if (nt == 0) return
            n1 = size(freg(${D}$)%lo, 2); n2 = size(freg(${D}$)%lo, 3)
            $:GPU_PARALLEL_LOOP(collapse=4, copyin='[tsl, toff, tsd]', private='[j]')
            do t = 1, nt
                do b = 0, n2 - 1
                    do a = 0, n1 - 1
                        do eq = 1, sys_size
                            j = toff(t) + eq + sys_size*(a + n1*b)
                            if (pack) then
                                if (tsd(t) == 1) then
                                    amr_rf_sq(j) = freg(${D}$)%lo(eq, a, b, tsl(t))
                                else
                                    amr_rf_sq(j) = freg(${D}$)%hi(eq, a, b, tsl(t))
                                end if
                            else if (tsd(t) == 1) then
                                freg(${D}$)%lo(eq, a, b, tsl(t)) = amr_rf_rq(j)
                            else
                                freg(${D}$)%hi(eq, a, b, tsl(t)) = amr_rf_rq(j)
                            end if
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

        end subroutine s_amr_freg_pack_d${D}$
    #:endfor

    !> Split-ownership level>=2 freg exchange as one wave, run once before the reflux fold; zero-copy like the faces wave, on tag
    !! band 1 (faces use band 0).
    impure subroutine s_amr_freg_wave()

#ifdef MFC_MPI
        integer  :: k, pblk, cowner, powner, nhr, j, kk2
        real(wp) :: w_lo(3), w_hi(3)

        if (num_procs == 1) return
        call s_amr_reg_prepare()
        call s_amr_refresh_lists()
        call s_amr_wave_open(amr_wave, amr_band_freg)
        nhr = 0
        ! receive side: freg faces of level>=2 children of my parents that another rank owns, straight into the child's
        ! register slot; faces no sibling weight selects are poisoned
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
            call s_amr_freg_faces_poison(w_lo > 0._wp, w_hi > 0._wp)
            call s_amr_freg_faces_req(2, cowner, w_lo > 0._wp, w_hi > 0._wp)
        end do
        ! send side: my level>=2 blocks whose parent lives elsewhere ship the faces the sibling weights select
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
            call s_amr_freg_faces_req(1, powner, w_lo > 0._wp, w_hi > 0._wp)
        end do
        call s_amr_wave_wait(amr_wave)
        do j = 1, nhr
            k = amr_fw_rblk(j)
            call s_amr_select_slot(k)
            pblk = amr_parent_blk(k)
            call s_amr_sibling_face_weights(k, pblk, w_lo, w_hi)
            call s_amr_freg_faces_update(.true., w_lo > 0._wp, w_hi > 0._wp)
        end do
#endif

    end subroutine s_amr_freg_wave

    !> Post one zero-copy request per flagged face register of the current slot (block k): dir 1 sends to / 2 receives from peer, in
    !! (direction, face) order.
    impure subroutine s_amr_freg_faces_req(dir, peer, f_lo, f_hi)

        integer, intent(in) :: dir, peer
        logical, intent(in) :: f_lo(3), f_hi(3)
        integer             :: cnt

        #:for D in [1, 2, 3]
            if (${D}$ <= num_dims) then
                cnt = size(freg(${D}$)%lo, 1)*size(freg(${D}$)%lo, 2)*size(freg(${D}$)%lo, 3)
                if (f_lo(${D}$)) call s_amr_wave_req_raw(amr_wave, dir, freg(${D}$)%lo(:,:,:,amr_reg_cur), cnt, peer)
                if (f_hi(${D}$)) call s_amr_wave_req_raw(amr_wave, dir, freg(${D}$)%hi(:,:,:,amr_reg_cur), cnt, peer)
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

    !> Debug builds: NaN-poison the unflagged face registers of the current slot before a receive lands, so a stray read of an
    !! unshipped face NaNs.
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

    !> Set the fine level's geometry (region, intersection, extents, bounds, coordinates) for the box lo:hi; arrays are
    !! preallocated, this updates metadata and coords. Collective. Also refreshes the amr_xchg_coarse_ghosts flag. A level-l block's
    !! fine extent is amr_ref_ratio**l * (coarse-region width) - 1, not amr_ref_ratio*width.
    impure subroutine s_set_amr_fine_geometry(lo, hi)

        integer, intent(in) :: lo(3), hi(3)
        integer             :: fext(3), nmar, bad_loc, pblk

        amr_slots(amr_cur)%region%lo = lo; amr_slots(amr_cur)%region%hi = hi
        amr_region_lo = lo; amr_region_hi = hi  ! mirror for m_amr_registers (avoids a use-cycle)
        amr_region_lo_all(:,amr_cur) = lo; amr_region_hi_all(:,amr_cur) = hi

        ! the owner holds every fine cell of the block, other ranks none; amr_isect_lo/hi is the block's coarse footprint on
        ! the owner (empty elsewhere)
        amr_rank_owns_block = (amr_block_owner(amr_cur) == proc_rank)
        pblk = 0
        if (amr_rank_owns_block) then
            amr_isect_lo = lo; amr_isect_hi = hi
            if (amr_block_level(amr_cur) >= 2) then
                ! multi-level: footprint in the parent block's fine frame (parent-fine index of L0 cell c is rr*(c - parent lo)),
                ! as s_amr_parent_foot derives it. rr is the global amr_ref_ratio: amr_slots(pblk) may be unallocated here
                pblk = f_amr_parent_block(amr_cur)
                call s_amr_parent_foot(amr_cur, pblk, amr_isect_lo, amr_isect_hi)
            end if
        else
            amr_isect_lo = 1; amr_isect_hi = 0  ! empty
        end if
        amr_isect_lo_all(:,amr_cur) = amr_isect_lo; amr_isect_hi_all(:,amr_cur) = amr_isect_hi
        amr_owns_all(amr_cur) = amr_rank_owns_block
        ! fine extents cover the whole block on the owner, -1 (empty) elsewhere
        fext = merge(amr_ref_ratio*max(amr_isect_hi - amr_isect_lo + 1, 0) - 1, 0, amr_dim)
        amr_slots(amr_cur)%m = fext(1); amr_slots(amr_cur)%n = fext(2); amr_slots(amr_cur)%p = fext(3)
        amr_slots(amr_cur)%idwbuff%beg = merge(-buff_size, 0, amr_dim); amr_slots(amr_cur)%idwbuff%end = merge(fext + buff_size, &
                  & 0, amr_dim)
        ! coords only on ranks with fine cells: replay the ancestor chain from the global L0 boundaries (the owner may hold no
        ! coarse coords nor the parent's slot)
        if (amr_rank_owns_block) then
            #:for D, X in [(1, 'x'), (2, 'y'), (3, 'z')]
                if (amr_dim(${D}$)) call s_amr_build_block_coords(amr_cur, amr_g${X}$cb, lbound(amr_g${X}$cb, 1), &
                    & amr_slots(amr_cur)%${X}$_cb, amr_slots(amr_cur)%${X}$_cc, amr_slots(amr_cur)%d${X}$, ${D}$)
            #:endfor
        end if

        ! The ghost prolongation reads up to nmar coarse cells past each face; if that leaves any rank's interior, the coarse
        ! cons ghosts must be halo-exchanged before every fill (the solver populates only prim ghosts). All ranks agree on the
        ! flag.
        nmar = (buff_size + amr_ref_ratio - 1)/amr_ref_ratio + 1
        bad_loc = 0
        if (amr_rank_owns_block .and. any(amr_dim .and. (amr_isect_lo - amr_sidx < nmar .or. amr_sidx + amr_ext &
            & - amr_isect_hi < nmar))) bad_loc = 1
        ! Accumulate (OR over blocks), do not reduce: the caller closes the scan with s_amr_reduce_xchg_flag, one collective
        ! per scan.
        amr_xchg_bad = max(amr_xchg_bad, bad_loc)

    end subroutine s_set_amr_fine_geometry

    !> Close a geometry scan: one allreduce of the accumulated flag, then reset. Call after every s_set_amr_fine_geometry loop; the
    !! flag is read by s_amr_exchange_coarse_cons_halo.
    impure subroutine s_amr_reduce_xchg_flag()

        integer :: bad_glb

        call s_mpi_allreduce_integer_max(amr_xchg_bad, bad_glb)
        amr_xchg_coarse_ghosts = bad_glb == 1
        amr_xchg_bad = 0

    end subroutine s_amr_reduce_xchg_flag

    !> Conservative prolongation of one variable from the patch qc into the fine interior: minmod-limited linear, or the ghost
    !! fill's limited quartic with amr_prolong = 2.
    impure subroutine s_prolong_one_var(qc, loc, ivar)

        type(scalar_field), intent(in) :: qc
        integer, intent(in)            :: loc, ivar                     !< flat-store slot and variable of the fine target
        integer                        :: fi, fj, fk, ci, cj, ck, ox, oy, oz, rrat, mm, nn, pp, il1, il2, il3, ky, kz, wy, wz
        real(wp)                       :: u0, xm, xp, ym, yp, zm, zp, xix, xiy, xiz
        real(wp)                       :: q5(-2:2), pl(-2:2), ln(-2:2)  !< staged x, y, z
        logical                        :: d2, d3, quart

        ! qc is the gathered patch amr_cg: amr_isect_lo + f/rr - amr_cpat_off = nmar + f/rr is the patch-local coarse index.
        ! Device kernel reading the patch's device mirror (pushed once per dispatch by s_interpolate_coarse_to_fine).

        ox = amr_cpat_off(1); oy = amr_cpat_off(2); oz = amr_cpat_off(3)
        rrat = amr_slots(amr_cur)%amr_ref_ratio
        mm = amr_slots(amr_cur)%m; nn = amr_slots(amr_cur)%n; pp = amr_slots(amr_cur)%p
        il1 = amr_isect_lo(1); il2 = amr_isect_lo(2); il3 = amr_isect_lo(3)
        d2 = n_glb > 0; d3 = p_glb > 0
        quart = amr_prolong == 2; wy = merge(2, 0, d2); wz = merge(2, 0, d3)
        $:GPU_PARALLEL_LOOP(collapse=3, private='[ci, cj, ck, xix, xiy, xiz, u0, xm, xp, ym, yp, zm, zp, ky, kz, q5, pl, ln]')
        do fk = 0, pp
            do fj = 0, nn
                do fi = 0, mm
                    ck = il3 + fk/rrat - oz; if (.not. d3) ck = 0
                    xiz = 0._wp; if (d3) xiz = (real(mod(fk, rrat), wp) - real(rrat - 1, wp)*0.5_wp)/real(rrat, wp)
                    cj = il2 + fj/rrat - oy; if (.not. d2) cj = 0
                    xiy = 0._wp; if (d2) xiy = (real(mod(fj, rrat), wp) - real(rrat - 1, wp)*0.5_wp)/real(rrat, wp)
                    ci = il1 + fi/rrat - ox
                    xix = (real(mod(fi, rrat), wp) - real(rrat - 1, wp)*0.5_wp)/real(rrat, wp)
                    if (quart) then
                        do kz = -wz, wz
                            do ky = -wy, wy
                                q5 = real(qc%sf(ci - 2:ci + 2,cj + ky, ck + kz), wp)
                                pl(ky) = f_amr_quartic_half(q5, xix)
                            end do
                            ln(kz) = pl(0)
                            if (d2) ln(kz) = f_amr_quartic_half(pl, xiy)
                        end do
                        u0 = ln(0)
                        if (d3) u0 = f_amr_quartic_half(ln, xiz)
                        amr_cons_st(fi, fj, fk, ivar, loc) = u0
                    else
                        u0 = real(qc%sf(ci, cj, ck), wp)
                        xm = real(qc%sf(ci - 1, cj, ck), wp); xp = real(qc%sf(ci + 1, cj, ck), wp)
                        ym = u0; yp = u0
                        if (d2) then
                            ym = real(qc%sf(ci, cj - 1, ck), wp); yp = real(qc%sf(ci, cj + 1, ck), wp)
                        end if
                        zm = u0; zp = u0
                        if (d3) then
                            zm = real(qc%sf(ci, cj, ck - 1), wp); zp = real(qc%sf(ci, cj, ck + 1), wp)
                        end if
                        amr_cons_st(fi, fj, fk, ivar, loc) = f_amr_lin(u0, xm, xp, ym, yp, zm, zp, xix, xiy, xiz)
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_prolong_one_var

    !> Quartic safety net (amr_prolong = 2): a parent whose children leave the admissible state (non-positive density or internal
    !! energy) is re-prolonged minmod-limited linear in every variable.
    impure subroutine s_prolong_guard(qc, loc)

        type(scalar_field), dimension(sys_size), intent(in) :: qc
        integer, intent(in)                                 :: loc
        integer                                             :: i, pi, pj, pk, a, b, c, fi, fj, fk, ci, cj, ck, rrat, ry, rz
        integer                                             :: ox, oy, oz, il1, il2, il3, nx, ny, nz, advb, adve, spb, spe
        real(wp)                                            :: u0, xm, xp, ym, yp, zm, zp, xix, xiy, xiz
        logical                                             :: d2, d3, bad, skip_alpha

        ox = amr_cpat_off(1); oy = amr_cpat_off(2); oz = amr_cpat_off(3)
        rrat = amr_slots(amr_cur)%amr_ref_ratio
        il1 = amr_isect_lo(1); il2 = amr_isect_lo(2); il3 = amr_isect_lo(3)
        d2 = n_glb > 0; d3 = p_glb > 0
        ry = merge(rrat, 1, d2); rz = merge(rrat, 1, d3)
        nx = (amr_slots(amr_cur)%m + 1)/rrat; ny = (amr_slots(amr_cur)%n + 1)/ry; nz = (amr_slots(amr_cur)%p + 1)/rz
        skip_alpha = num_fluids > 1 .and. (.not. bubbles_lagrange)
        advb = eqn_idx%adv%beg; adve = eqn_idx%adv%end; spb = eqn_idx%species%beg; spe = eqn_idx%species%end
        $:GPU_PARALLEL_LOOP(collapse=3, &
                            & private='[i, a, b, c, fi, fj, fk, ci, cj, ck, xix, xiy, xiz, u0, xm, xp, ym, yp, zm, zp, bad]')
        do pk = 0, nz - 1
            do pj = 0, ny - 1
                do pi = 0, nx - 1
                    bad = .false.
                    do c = 0, rz - 1
                        do b = 0, ry - 1
                            do a = 0, rrat - 1
                                fi = pi*rrat + a; fj = pj*ry + b; fk = pk*rz + c
                                if (.not. f_amr_cons_ok(fi, fj, fk, loc)) bad = .true.
                            end do
                        end do
                    end do
                    if (bad) then
                        ci = il1 + pi - ox; cj = 0; ck = 0
                        if (d2) cj = il2 + pj - oy
                        if (d3) ck = il3 + pk - oz
                        do c = 0, rz - 1
                            do b = 0, ry - 1
                                do a = 0, rrat - 1
                                    fi = pi*rrat + a; fj = pj*ry + b; fk = pk*rz + c
                                    xix = (real(a, wp) - real(rrat - 1, wp)*0.5_wp)/real(rrat, wp)
                                    xiy = 0._wp; if (d2) xiy = (real(b, wp) - real(rrat - 1, wp)*0.5_wp)/real(rrat, wp)
                                    xiz = 0._wp; if (d3) xiz = (real(c, wp) - real(rrat - 1, wp)*0.5_wp)/real(rrat, wp)
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = 1, sys_size
                                        if (skip_alpha .and. i >= advb .and. i <= adve) cycle
                                        if (chemistry .and. i >= spb .and. i <= spe) cycle
                                        u0 = real(qc(i)%sf(ci, cj, ck), wp)
                                        xm = real(qc(i)%sf(ci - 1, cj, ck), wp); xp = real(qc(i)%sf(ci + 1, cj, ck), wp)
                                        ym = u0; yp = u0
                                        if (d2) then
                                            ym = real(qc(i)%sf(ci, cj - 1, ck), wp); yp = real(qc(i)%sf(ci, cj + 1, ck), wp)
                                        end if
                                        zm = u0; zp = u0
                                        if (d3) then
                                            zm = real(qc(i)%sf(ci, cj, ck - 1), wp); zp = real(qc(i)%sf(ci, cj, ck + 1), wp)
                                        end if
                                        amr_cons_st(fi, fj, fk, i, loc) = f_amr_lin(u0, xm, xp, ym, yp, zm, zp, xix, xiy, xiz)
                                    end do
                                end do
                            end do
                        end do
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_prolong_guard

    !> Conservative prolongation: fill amr_fine interior from coarse (minmod-limited linear, or quartic with amr_prolong = 2);
    !! symmetric child offsets make the children average to the coarse value. Multi-fluid volume fractions take the sum- preserving
    !! closure path.
    impure subroutine s_interpolate_coarse_to_fine()

        integer :: i

        ! the prolong kernels read the gathered patch on the device

        integer :: loc_cur  !< plain scalar for device routines (NVHPC -Minline)

        do i = 1, sys_size
            ! Lagrangian bubbles: alphas sum to beta, not 1, so each alpha prolongs plainly
            if (num_fluids > 1 .and. (.not. bubbles_lagrange) .and. i >= eqn_idx%adv%beg .and. i <= eqn_idx%adv%end) cycle
            if (chemistry .and. i >= eqn_idx%species%beg .and. i <= eqn_idx%species%end) cycle  ! closure below
            loc_cur = amr_loc_of(amr_cur)
            call s_prolong_one_var(amr_cg(i), loc_cur, i)
        end do
        loc_cur = amr_loc_of(amr_cur)
        if (num_fluids > 1 .and. (.not. bubbles_lagrange)) call s_prolong_alphas_closure(amr_cg, loc_cur)
        loc_cur = amr_loc_of(amr_cur)
        if (chemistry) call s_prolong_species_closure(amr_cg, loc_cur)
        if (amr_prolong == 2) call s_prolong_guard(amr_cg, loc_cur)

    end subroutine s_interpolate_coarse_to_fine

    !> Sum-preserving volume-fraction prolongation (num_fluids > 1): fluids adv%beg..adv%end-1 are interpolated with minmod slopes
    !! under a shared per-cell limiter switch and clamped to [0,1]; the last fluid closes alpha_n = 1 - sum(others). For >2 fluids
    !! residual closure undershoot is handled by mpp_lim. Index mapping as s_prolong_one_var.
    impure subroutine s_prolong_alphas_closure(qc, loc)

        type(scalar_field), dimension(sys_size), intent(in) :: qc
        integer, intent(in)                                 :: loc
        integer                                             :: fi, fj, fk, ci, cj, ck, ox, oy, oz, i
        integer                                             :: rrat, mm, nn, pp, il1, il2, il3, advb, adve
        real(wp)                                            :: xix, xiy, xiz, u0, xm, xp, ym, yp, zm, zp, av, asum
        logical                                             :: shx, shy, shz, d2, d3

        ! qc is the gathered patch amr_cg (patch-frame offset). Device kernel; the shared limiter switch is inlined.

        ox = amr_cpat_off(1); oy = amr_cpat_off(2); oz = amr_cpat_off(3)
        rrat = amr_slots(amr_cur)%amr_ref_ratio
        mm = amr_slots(amr_cur)%m; nn = amr_slots(amr_cur)%n; pp = amr_slots(amr_cur)%p
        il1 = amr_isect_lo(1); il2 = amr_isect_lo(2); il3 = amr_isect_lo(3)
        d2 = n_glb > 0; d3 = p_glb > 0
        advb = eqn_idx%adv%beg; adve = eqn_idx%adv%end
        $:GPU_PARALLEL_LOOP(collapse=3, &
                            & private='[ci, cj, ck, xix, xiy, xiz, u0, xm, xp, ym, yp, zm, zp, av, asum, shx, shy, shz, i]')
        do fk = 0, pp
            do fj = 0, nn
                do fi = 0, mm
                    ck = il3 + fk/rrat - oz; if (.not. d3) ck = 0
                    xiz = 0._wp; if (d3) xiz = (real(mod(fk, rrat), wp) - real(rrat - 1, wp)*0.5_wp)/real(rrat, wp)
                    cj = il2 + fj/rrat - oy; if (.not. d2) cj = 0
                    xiy = 0._wp; if (d2) xiy = (real(mod(fj, rrat), wp) - real(rrat - 1, wp)*0.5_wp)/real(rrat, wp)
                    ci = il1 + fi/rrat - ox
                    xix = (real(mod(fi, rrat), wp) - real(rrat - 1, wp)*0.5_wp)/real(rrat, wp)
                    ! shared per-cell limiter switch: per dim, slopes stay on only if no fluid's centered differences change
                    ! sign there
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
                        xm = u0; xp = u0
                        if (shx) then
                            xm = real(qc(i)%sf(ci - 1, cj, ck), wp); xp = real(qc(i)%sf(ci + 1, cj, ck), wp)
                        end if
                        ym = u0; yp = u0
                        if (shy) then
                            ym = real(qc(i)%sf(ci, cj - 1, ck), wp); yp = real(qc(i)%sf(ci, cj + 1, ck), wp)
                        end if
                        zm = u0; zp = u0
                        if (shz) then
                            zm = real(qc(i)%sf(ci, cj, ck - 1), wp); zp = real(qc(i)%sf(ci, cj, ck + 1), wp)
                        end if
                        av = min(max(f_amr_lin(u0, xm, xp, ym, yp, zm, zp, xix, xiy, xiz), 0._wp), 1._wp)
                        amr_cons_st(fi, fj, fk, i, loc) = av
                        asum = asum + av
                    end do
                    amr_cons_st(fi, fj, fk, adve, loc) = 1._wp - asum
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_prolong_alphas_closure

    !> Species prolongation closure (chemistry): each rho*Y_k is minmod-prolonged and clamped non-negative, then rescaled so
    !! sum_k(rho*Y_k) equals the prolonged continuity density. Index mapping as s_prolong_one_var; cont is prolonged before this
    !! runs.
    impure subroutine s_prolong_species_closure(qc, loc)

        type(scalar_field), dimension(sys_size), intent(in) :: qc
        integer, intent(in)                                 :: loc
        integer                                             :: fi, fj, fk, ci, cj, ck, ox, oy, oz, i
        integer                                             :: rrat, mm, nn, pp, il1, il2, il3, spb, spe, cte
        real(wp)                                            :: xix, xiy, xiz, u0, xm, xp, ym, yp, zm, zp, av, rsum, rscale
        logical                                             :: d2, d3

        ! qc is the gathered patch amr_cg (patch-frame offset). Device kernel; the rescale reads only the thread's own cell.

        ox = amr_cpat_off(1); oy = amr_cpat_off(2); oz = amr_cpat_off(3)
        rrat = amr_slots(amr_cur)%amr_ref_ratio
        mm = amr_slots(amr_cur)%m; nn = amr_slots(amr_cur)%n; pp = amr_slots(amr_cur)%p
        il1 = amr_isect_lo(1); il2 = amr_isect_lo(2); il3 = amr_isect_lo(3)
        d2 = n_glb > 0; d3 = p_glb > 0
        spb = eqn_idx%species%beg; spe = eqn_idx%species%end; cte = eqn_idx%cont%end
        $:GPU_PARALLEL_LOOP(collapse=3, private='[ci, cj, ck, xix, xiy, xiz, u0, xm, xp, ym, yp, zm, zp, av, rsum, rscale, i]')
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
                        xm = real(qc(i)%sf(ci - 1, cj, ck), wp); xp = real(qc(i)%sf(ci + 1, cj, ck), wp)
                        ym = u0; yp = u0
                        if (d2) then
                            ym = real(qc(i)%sf(ci, cj - 1, ck), wp); yp = real(qc(i)%sf(ci, cj + 1, ck), wp)
                        end if
                        zm = u0; zp = u0
                        if (d3) then
                            zm = real(qc(i)%sf(ci, cj, ck - 1), wp); zp = real(qc(i)%sf(ci, cj, ck + 1), wp)
                        end if
                        av = max(f_amr_lin(u0, xm, xp, ym, yp, zm, zp, xix, xiy, xiz), 0._wp)
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
        ! Prolong every block from its coarse patch, assembled by the level-1 fill wave (no coarse-ghost halo needed). Runs
        ! before s_initialize_gpu_vars and the wave packs from the device, so push the ICs first. All ranks enter the wave;
        ! only owners prolong.
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
            call s_interpolate_coarse_to_fine()
        end do
        call s_amr_fill_wave_done()
        if (amr_max_level >= 2) call s_amr_build_static_multilevel()
        call s_amr_select_slot(1)

    end subroutine s_populate_amr_fine

    !> Build the static multi-level hierarchy (amr_regrid_int = 0): nest one level-2 block inside level-1 block 1 by a fixed
    !! geometric inset, prolong the parent state into it and keep it persistent.
    impure subroutine s_amr_build_static_multilevel()

        integer :: L2, n1, par, inset(3)

        if (amr_max_level < 2) return
        n1 = amr_num_blocks
        if (n1 < 1) return
        ! the static hierarchy nests one level-2 block; without pool room it would silently refine only to level 1, so abort
        ! here (n1 is unknown at checker time). Replicated inputs: collective-safe.
        if (n1 + 1 > amr_max_blocks) call s_mpi_abort('amr static multi-level (amr_max_level > 1, ' &
            & // 'amr_regrid_int = 0): amr_max_blocks is too small to nest the level-2 block (need >= level-1 block count + 1); ' &
            & // 'increase amr_max_blocks')
        L2 = n1 + 1
        par = 1
        inset = merge(max((amr_region_hi_all(:,par) - amr_region_lo_all(:,par) + 1)/4, amr_cpat_mar), 0, amr_dim)
        amr_region_lo_all(:,L2) = amr_region_lo_all(:,par) + inset
        amr_region_hi_all(:,L2) = amr_region_hi_all(:,par) - inset
        ! Guard the fixed-inset box against configs this builder cannot represent: (a) a level-1 block smaller than 2*inset
        ! inverts the box; (b) a level-2 L0 extent > amr_maxc_fit/2 overruns the creg register (0:amr_maxc_fit-1) in the L2->L1
        ! reflux capture. Replicated inputs: collective-safe.
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
        ! the parent-fill wave assembles the L2 patch from the parent's fine array (all ranks enter; the owner prolongs in place)
        call s_amr_parent_fill_exchange(2, .true.)
        call s_amr_select_slot(L2)
        if (amr_rank_owns_block) then
            call s_amr_parent_fill_consume(L2, .true.)
            call s_interpolate_coarse_to_fine()
        end if
        call s_amr_fill_wave_done()
        call s_amr_select_slot(1)

    end subroutine s_amr_build_static_multilevel

    !> Sync level lev >= 2 into its parents: restrict its blocks into the covered parent cells, then state-reflux each into its
    !! parent (step-integral fluxes scaled by dt_reflux; 1 when registers carry dt).
    impure subroutine s_amr_sync_level(lev, dt_reflux)

        integer, intent(in)  :: lev
        real(wp), intent(in) :: dt_reflux
        integer              :: k, io, ifc, ko, kf

        call s_phase_tic(PH_RESTR)
        call s_amr_restrict_parent_wave(lev)
        call s_phase_toc(PH_RESTR)
        ! s_amr_reflux_to_parent acts only on the parent's owner. The two block lists overlap and a duplicate visit
        ! would reflux twice, so walk their union with a two-cursor merge from the tails (descending block order: parent
        ! registers accumulate in that order).
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

    end subroutine s_amr_sync_level

    !> Lock-step restrict fold as plan-based waves : per level, finest first, the split-ownership child->parent folds
    !! (s_amr_restrict_parent_wave) then the per-box reflux-to-parent applies (freg already exchanged by s_amr_freg_wave); then the
    !! level-1 -> L0 scatter (s_amr_restrict_l1_wave). Both sides derive identical transfer lists from replicated metadata; one rank
    !! has no peers. Within a level the covered targets are disjoint.
    impure subroutine s_amr_restrict_wave(coarse_tgt, dt_reflux)

        type(scalar_field), dimension(sys_size), intent(inout) :: coarse_tgt
        real(wp), intent(in)                                   :: dt_reflux
        integer                                                :: lev

        call s_amr_refresh_lists()
        ! under subcycling the levels >= 2 were synced after each of their parent's substeps
        if (.not. amr_subcycle) then
            do lev = amr_max_level, 2, -1
                call s_amr_sync_level(lev, dt_reflux)
            end do
        end if
        call s_phase_tic(PH_RESTR)
        call s_amr_restrict_l1_wave(coarse_tgt)
        call s_phase_toc(PH_RESTR)
        call s_amr_select_slot(1)

    end subroutine s_amr_restrict_wave

    !> One per-level wave of the split-ownership level>=2 child->parent restrict folds. Co-located folds run inline on the child
    !! owner; cross-rank folds ship the parent-frame covered box, one message per (child-owner, parent-owner) peer. Both sides walk
    !! the same block list ascending with per-peer running offsets.
    impure subroutine s_amr_restrict_parent_wave(lev)

        integer, intent(in) :: lev
        integer             :: h_bl_1, h_bl_2, h_bl_3, h_ploc, h_cloc  !< plain scalars for device routines (NVHPC -Minline)
        integer             :: k, pblk, cowner, powner, rr, idx, lo, hi, kk
        integer             :: plo(3), phi(3), bl(3), bh(3)

        rr = amr_ref_ratio
        call s_amr_wave_open(amr_wave, amr_band_restrict_parent)
        ! send side: owned level-lev blocks whose parent lives elsewhere ship their whole footprint to the parent's owner; a
        ! co-located parent is folded in place
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
                h_ploc = amr_loc_of(pblk); h_cloc = amr_loc_of(k)
                call s_amr_restrict_device_st(h_ploc, h_cloc, plo, phi, plo, rr, 0, 0, 0)
                cycle
            end if
            call s_amr_wave_add(amr_wsend, powner, k, plo, phi, sys_size*product(phi - plo + 1))
        end do
        call s_amr_wave_close(amr_wsend, amr_fw_sq, amr_fw_dev)
        ! receive side: level-lev children of my parents that another rank owns
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
        call s_amr_wave_post(amr_wave, amr_wrecv, amr_fw_rq, amr_fw_dev)
        do idx = 1, amr_wsend%nx
            call s_amr_wave_slice(amr_wsend, idx, lo, hi)
            bl = amr_wsend%bl(:,idx); bh = amr_wsend%bh(:,idx); h_cloc = amr_loc_of(amr_wsend%blk(idx))
            call s_amr_restrict_device_wire(amr_fw_sq(lo:hi), h_cloc, bl, bh, bl, rr)
        end do
        call s_amr_wave_send(amr_wave, amr_wsend, amr_fw_sq, amr_fw_dev)
        call s_amr_wave_wait(amr_wave)
        do idx = 1, amr_wrecv%nx
            call s_amr_wave_slice(amr_wrecv, idx, lo, hi)
            bl = amr_wrecv%bl(:,idx); bh = amr_wrecv%bh(:,idx)
            h_bl_1 = bl(1)
            h_bl_2 = bl(2)
            h_bl_3 = bl(3)
            h_ploc = amr_loc_of(amr_parent_blk(amr_wrecv%blk(idx)))
            call s_amr_pack_unpack_block_st(h_ploc, h_bl_1, h_bl_2, h_bl_3, bh(1) - h_bl_1, bh(2) - h_bl_2, bh(3) - h_bl_3, &
                                            & amr_fw_rq(lo:hi), .false.)
        end do

    end subroutine s_amr_restrict_parent_wave

    !> Level-1 -> L0 covered-cell scatter as one wave: every owned level-1 block's covered slabs ship to each listed coarse owner in
    !! one message per peer; the owner-local overwrite stays per block during the pack walk. The receiver plan (my interior x
    !! region(k) over level-1 blocks I do not own) equals the sender's list membership.
    impure subroutine s_amr_restrict_l1_wave(coarse_tgt)

        type(scalar_field), dimension(sys_size), intent(inout) :: coarse_tgt
        integer                                                :: k, owner, rr, idx, r, lo, hi, o1, o2, o3, cur, kk
        integer                                                :: rlo(3), rhi(3), ilo(3), ihi(3), milo(3), mihi(3), bl(3), bh(3)

        call s_amr_wave_open(amr_wave, amr_band_restrict_l1)
        o1 = amr_sidx(1); o2 = amr_sidx(2); o3 = amr_sidx(3)
        ! block set changed: rebuild the cached overlap-rank lists (local, replicated)
        if (amr_seam_pairs_dirty .or. amr_seam_pairs_nblk /= amr_num_blocks) call s_amr_build_seam_pairs()
        call s_amr_rank_interior(proc_rank, milo, mihi)
        ! send side: for every owned level-1 block, its covered slab inside each listed coarse owner's interior
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
        ! receive side: every level-1 block I do not own that covers my interior, from its owner
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
        call s_amr_wave_post(amr_wave, amr_wrecv, amr_fw_rq, amr_fw_dev)
        call s_amr_restrict_l1_batched(coarse_tgt, milo, mihi)

    end subroutine s_amr_restrict_l1_wave

    !> Rest of the level-1 restrict wave (receives posted): owned-block overwrite, send packs and received unpacks as one launch
    !! each over slab rows.
    impure subroutine s_amr_restrict_l1_batched(coarse_tgt, milo, mihi)

        type(scalar_field), dimension(sys_size), intent(inout) :: coarse_tgt
        integer, intent(in)                                    :: milo(3), mihi(3)
        integer                                                :: k, kk, t, rlo(3), rhi(3), bl(3), bh(3)

        call s_amr_refresh_my_blocks()
        do kk = 1, amr_n_my
            k = amr_my_blk(kk)
            if (amr_block_level(k) /= 1) cycle
            call s_amr_region_box(k, rlo, rhi)
            call s_amr_box_isect(rlo, rhi, milo, mihi, bl, bh)
            if (all(bl <= bh)) call s_amr_slab_add(bl, bh, amr_loc_of(k), rlo, [0, 0, 0])
        end do
        call s_amr_restrict_rows_sf(coarse_tgt, amr_sidx(1), amr_sidx(2), amr_sidx(3))
        do t = 1, amr_wsend%nx
            k = amr_wsend%blk(t)
            call s_amr_region_box(k, rlo, rhi)
            call s_amr_slab_add(amr_wsend%bl(:,t), amr_wsend%bh(:,t), amr_loc_of(k), rlo, [amr_wsend%off(t), 0, 0])
        end do
        if (amr_wsend%nx > 0) call s_amr_restrict_rows_wire(amr_fw_sq(1:amr_wsend%words))
        call s_amr_wave_send(amr_wave, amr_wsend, amr_fw_sq, amr_fw_dev)
        call s_amr_wave_wait(amr_wave)
        do t = 1, amr_wrecv%nx
            call s_amr_slab_add(amr_wrecv%bl(:,t), amr_wrecv%bh(:,t), amr_wrecv%off(t), [0, 0, 0], amr_sidx)
        end do
        if (amr_wrecv%nx > 0) call s_amr_unpack_rows(amr_fw_rq(1:amr_wrecv%words), coarse_tgt)

    end subroutine s_amr_restrict_l1_batched

    !> Sibling-seam face weights for level>=2 block kb under parent pblk: 0 on a face shared with a same-parent sibling tile
    !! (fine-fine; refluxing there double-writes), 1 otherwise. Uses replicated metadata only, and both the reflux apply and the
    !! freg wire-skip call it, so shipped and consumed faces agree.
    impure subroutine s_amr_sibling_face_weights(kb, pblk, w_lo, w_hi)

        integer, intent(in)   :: kb, pblk
        real(wp), intent(out) :: w_lo(3), w_hi(3)
        integer               :: c, y, d

        ! iterate pblk's cached children (same level by construction); do not refresh the lists here, a caller may be iterating them

        w_lo = 1._wp; w_hi = 1._wp
        if (pblk <= 0) return  ! orphan parent: no siblings; the CSR would index ptr(-1)
        do c = amr_child_ptr(pblk - 1) + 1, amr_child_ptr(pblk)
            y = amr_child_idx(c)
            if (y == kb) cycle
            do d = 1, num_dims
                if (f_amr_seam(kb, y, d)) w_hi(d) = 0._wp  ! sibling above: shared high face
                if (f_amr_seam(y, kb, d)) w_lo(d) = 0._wp  ! sibling below: shared low face
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
        integer              :: loc_pblk  !< plain scalar for device routines (NVHPC -Minline)

        call s_amr_refresh_lists()  ! caches the parent (f_amr_parent_block is a global scan)
        pblk = amr_parent_blk(amr_cur)
        if (pblk <= 0) return  ! orphan: no parent to reflux into
        if (amr_block_owner(pblk) /= proc_rank) return  ! the parent's owner applies the correction
        ! A face shared with an adjacent same-parent sibling tile is fine-fine, not c/f: its outside parent cell is covered by
        ! the sibling's restrict, so skip it (weight 0). No siblings: all weights 1.
        call s_amr_sibling_face_weights(amr_cur, pblk, w_lo, w_hi)
        ! parent-fine frame for the shared reflux kernel: outside cell = isect boundary +/-1, creg-local loop range 0:extent,
        ! transverse write at the isect origin, per-face parent-fine dx. Footprint from s_amr_parent_foot (replicated), not
        ! amr_isect_lo/hi, which is empty on the parent's owner; rr is the global amr_ref_ratio since amr_slots(amr_cur) may be
        ! unallocated.
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
        loc_pblk = amr_loc_of(pblk)
        call s_amr_br_load_faces(loc_pblk, olo, ohi, glo, ghi, woff, w_lo, w_hi)
        call s_amr_reflux_apply_faces(amr_cons_br, amr_reg_cur, amr_ref_ratio, dt_reflux, olo, ohi, glo, ghi, woff, w_lo, w_hi, &
                                      & mlo, mhi)
        loc_pblk = amr_loc_of(pblk)
        call s_amr_br_store_faces(loc_pblk, olo, ohi, glo, ghi, woff, w_lo, w_hi)

    end subroutine s_amr_reflux_to_parent

    !> Device-native restriction of the fine block (flat store) to coarse averages over the covered cells [bl:bh] (global; block
    !! origin rlo, ratio rr); child-sum order ddk, ddj, ddi, /nchild. One body for two destinations so results match bit-for-bit: a
    !! parent block in the flat store (`_st`, local origin o) or the wire buffer buf (`_wire`; ci fastest, then cj, ck, i; wp).
    !! Never push the whole coarse array to the device.
    #:for SFX in ['st', 'wire']
        #:set OARGS = '' if SFX == 'wire' else ', o1, o2, o3'
        impure subroutine s_amr_restrict_device_${SFX}$(${ {'st': 'ctloc', 'wire': 'buf'}[SFX] }$, loc, bl, bh, rlo, rr${OARGS}$)

            #:if SFX == 'st'
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
                            #:if SFX == 'st'
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
