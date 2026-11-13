!>
!!@file
!!@brief Contains module m_amr_exchange

#! AMD OpenMP: assert allocatables present on every kernel (see OMP_DEFAULT_STR); a kernel may name only module arrays allocated
#! under its own launch condition.
#:set MFC_OMP_PRESENT_ALLOCATABLE = True
#:include 'macros.fpp'

!> @brief Parent/child data exchange: gather plans, pack/unpack and the fill waves.
module m_amr_exchange

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy, only: s_mpi_abort  ! @:ASSERT expands to it
    use m_amr_registers, only: s_amr_parent_foot
    use m_phase_timing
    use m_amr_state
    use m_amr_wave
    use m_amr_distribution
    use m_amr_ghost
    use m_amr_seam, only: s_amr_build_seam_pairs

    implicit none

    private
    public :: s_amr_fill_wave_done, s_amr_l1_fill_exchange, s_amr_l1_fill_consume, s_amr_parent_fill_exchange, &
        & s_amr_parent_fill_consume, s_amr_parent_fill_wave, s_amr_stage_fill_wave, s_amr_pack_unpack_block_st, &
        & s_amr_unpack_rows, s_amr_restrict_rows_sf, s_amr_restrict_rows_wire, s_amr_l1_pool_gather, s_amr_l1_pool_fill

contains

    !> Parent value for a level >= 2 ghost fill: the current state, or its blend with the substep-start backup when th < 1.
    pure real(wp) function f_amr_par_val(g1, g2, g3, i, qp, th) result(v)
        $:GPU_ROUTINE(parallelism='[seq]')
        integer, intent(in)  :: g1, g2, g3, i, qp
        real(wp), intent(in) :: th

        if (th < 1._wp) then
            v = th*real(amr_cons_st(g1, g2, g3, i, qp), wp) + (1._wp - th)*real(amr_stor_st(g1, g2, g3, i, qp), wp)
        else
            v = real(amr_cons_st(g1, g2, g3, i, qp), wp)
        end if

    end function f_amr_par_val

    !> Sub-box variants of the parent-patch pack/copy for the ring-clipped parent-fill wave; bl:bh are patch-local, the buffer
    !! layout matches the full-patch kernels.
    impure subroutine s_amr_pack_parent_box_device(qp, bl, bh, buf)

        integer, intent(in)                 :: qp  !< parent's flat-store slot
        integer, intent(in)                 :: bl(3), bh(3)
        real(wp), intent(inout), contiguous :: buf(:)
        integer                             :: i, g1, g2, g3, o1, o2, o3, n1, n2, n3, l1, l2, l3, u1, u2, u3
        real(wp)                            :: th

        o1 = amr_cpat_off(1); o2 = amr_cpat_off(2); o3 = amr_cpat_off(3)
        l1 = bl(1); l2 = bl(2); l3 = bl(3); u1 = bh(1); u2 = bh(2); u3 = bh(3)
        n1 = u1 - l1 + 1; n2 = u2 - l2 + 1; n3 = u3 - l3 + 1
        th = amr_par_theta
        $:GPU_PARALLEL_LOOP(collapse=4, copyout='[buf]', copyin='[th]')
        do i = 1, sys_size
            do g3 = l3, u3
                do g2 = l2, u2
                    do g1 = l1, u1
                        buf(1 + (g1 - l1) + n1*((g2 - l2) + n2*((g3 - l3) + n3*(i - 1)))) = f_amr_par_val(g1 + o1, g2 + o2, &
                            & g3 + o3, i, qp, th)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_pack_parent_box_device

    impure subroutine s_amr_copy_parent_box(qp, bl, bh)

        integer, intent(in) :: qp  !< parent's flat-store slot
        integer, intent(in) :: bl(3), bh(3)
        integer             :: i, g1, g2, g3, o1, o2, o3, l1, l2, l3, u1, u2, u3
        real(wp)            :: th

        o1 = amr_cpat_off(1); o2 = amr_cpat_off(2); o3 = amr_cpat_off(3)
        l1 = bl(1); l2 = bl(2); l3 = bl(3); u1 = bh(1); u2 = bh(2); u3 = bh(3)
        th = amr_par_theta
        $:GPU_PARALLEL_LOOP(collapse=4, copyin='[th]')
        do i = 1, sys_size
            do g3 = l3, u3
                do g2 = l2, u2
                    do g1 = l1, u1
                        amr_cg(i)%sf(g1, g2, g3) = f_amr_par_val(g1 + o1, g2 + o2, g3 + o3, i, qp, th)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_copy_parent_box

    !> Pack (pack = T) the global box [bl:bh] of q_coarse into the wire buffer buf, or unpack it into amr_cg at patch-local index g
    !! - o (pack = F; stp cast). One element order (i slowest, x fastest), the same wp cast as the host pack in
    !! s_amr_gather_coarse_patch.
    #:for SFX, PK in [('pack', True), ('unpack', False)]
        impure subroutine s_amr_${SFX}$_box_device(${'q_coarse, ' if PK else ''}$bl, bh, ${'o1, o2, o3' if PK else 'off'}$, buf)

            #:if PK
                type(scalar_field), dimension(sys_size), intent(in) :: q_coarse
                integer, intent(in)                                 :: bl(3), bh(3), o1, o2, o3
                real(wp), intent(inout), contiguous                 :: buf(:)
            #:else
                integer, intent(in)              :: bl(3), bh(3), off(3)
                real(wp), intent(in), contiguous :: buf(:)
            #:endif
            integer :: i, g1, g2, g3, l1, l2, l3, u1, u2, u3, n1, n2, n3
            #:if not PK
                integer :: o1, o2, o3
            #:endif

            l1 = bl(1); l2 = bl(2); l3 = bl(3); u1 = bh(1); u2 = bh(2); u3 = bh(3)
            n1 = u1 - l1 + 1; n2 = u2 - l2 + 1; n3 = u3 - l3 + 1
            #:if not PK
                o1 = off(1); o2 = off(2); o3 = off(3)
            #:endif
            #:if PK
                $:GPU_PARALLEL_LOOP(collapse=4, copyout='[buf]')
            #:else
                $:GPU_PARALLEL_LOOP(collapse=4, copyin='[buf]')
            #:endif
            do i = 1, sys_size
                do g3 = l3, u3
                    do g2 = l2, u2
                        do g1 = l1, u1
                            #:if PK
                                buf(1 + (g1 - l1) + n1*((g2 - l2) + n2*((g3 - l3) + n3*(i - 1)))) = real(q_coarse(i)%sf(g1 - o1, &
                                    & g2 - o2, g3 - o3), wp)
                            #:else
                                amr_cg(i)%sf(g1 - o1, g2 - o2, &
                                       & g3 - o3) = real(buf(1 + (g1 - l1) + n1*((g2 - l2) + n2*((g3 - l3) + n3*(i - 1)))), stp)
                            #:endif
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

        end subroutine s_amr_${SFX}$_box_device
    #:endfor

    !> Ring clip: split the hollow shell of patch [plo:phi] minus open core [clo:chi] into at most 6 disjoint slabs in fixed order
    !! (x, y within core x, z within core x and y). Collapsed dims pass 0; the width-1 double cover is clamped. Both sides derive
    !! the list from replicated metadata.
    impure subroutine s_amr_shell_slabs(plo, phi, clo, chi, ns, sb, se, cells)

        integer, intent(in)  :: plo(3), phi(3), clo(3), chi(3)
        integer, intent(out) :: ns, sb(3, 6), se(3, 6), cells
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
            if (any(cb(:,s) > ce(:,s))) cycle
            ns = ns + 1
            sb(:,ns) = cb(:,s); se(:,ns) = ce(:,s)
            words = words + product(int(ce(:,s) - cb(:,s) + 1, 8))
        end do
        ! slabs must tile the shell exactly
        do s = 1, ns - 1
            do ss = s + 1, ns
                @:ASSERT(any(max(sb(:, s), sb(:, ss)) > min(se(:, s), se(:, ss))), "shell slabs: overlap")
            end do
        end do
        patchw = int(phi(1) - plo(1) + 1, 8)*int(phi(2) - plo(2) + 1, 8)*int(phi(3) - plo(3) + 1, 8)
        corew = int(max(chi(1) - clo(1) + 1, 0), 8)*int(max(chi(2) - clo(2) + 1, 0), 8)*int(max(chi(3) - clo(3) + 1, 0), 8)
        @:ASSERT(words == patchw - corew, "shell slabs: coverage mismatch")
        cells = int(words)

    end subroutine s_amr_shell_slabs

    !> Intersect the shell slabs with box [bl:bh], keeping order; also returns the total cell count.
    impure subroutine s_amr_shell_clip(ns, sb, se, bl, bh, ms, tb, te, cells)

        integer, intent(in)  :: ns, sb(3, 6), se(3, 6), bl(3), bh(3)
        integer, intent(out) :: ms, tb(3, 6), te(3, 6), cells
        integer              :: s, l(3), u(3)

        ms = 0; cells = 0
        do s = 1, ns
            l = max(sb(:,s), bl); u = min(se(:,s), bh)
            if (any(l > u)) cycle
            ms = ms + 1
            tb(:,ms) = l; te(:,ms) = u
            cells = cells + product(u - l + 1)
        end do

    end subroutine s_amr_shell_clip

#ifdef MFC_DEBUG
    !> Debug only: NaN-flood amr_cg before a clipped gather so any read of an unshipped cell NaNs the ghost fill.
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
#endif

    !> Ring-clipped own-box copy on device: queue shell-slab/own-box intersections as slab rows and copy q_coarse into amr_cg in one
    !! fused kernel.
    impure subroutine s_amr_gather_own_rows(q_coarse, cg, o1, o2, o3)

        type(scalar_field), dimension(sys_size), intent(in)    :: q_coarse
        type(scalar_field), dimension(sys_size), intent(inout) :: cg  !< amr_cg, passed so a regrown pool maps as an argument
        integer, intent(in)                                    :: o1, o2, o3
        integer                                                :: i, s, g, r, n1, n2, g1, g2, g3, stot, nsl, lo, hi, md

        stot = f_amr_slab_push(nsl)
        if (stot == 0) return  ! no level-1 block on this rank: the store may be unallocated
        $:GPU_PARALLEL_LOOP(collapse=2, private='[s, lo, hi, md, r, n1, n2, g1, g2, g3]')
        do i = 1, sys_size
            do g = 0, stot - 1
                lo = 1; hi = nsl
                do while (lo < hi)
                    md = (lo + hi + 1)/2
                    if (amr_slab_tab(7, md) <= g) then; lo = md; else; hi = md - 1; end if
                end do
                s = lo
                r = g - amr_slab_tab(7, s)
                n1 = amr_slab_tab(4, s) - amr_slab_tab(1, s) + 1; n2 = amr_slab_tab(5, s) - amr_slab_tab(2, s) + 1
                g1 = amr_slab_tab(1, s) + mod(r, n1)
                g2 = amr_slab_tab(2, s) + mod(r/n1, n2)
                g3 = amr_slab_tab(3, s) + r/(n1*n2)
                cg(i)%sf(g1 - amr_slab_tab(13, s), g2 - amr_slab_tab(14, s), g3 - amr_slab_tab(15, s)) = q_coarse(i)%sf(g1 - o1, &
                   & g2 - o2, g3 - o3)
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_gather_own_rows

    !> Fused-pack plan for transfers 1:nt: slab corner, extents, absolute pool payload offset and exclusive element prefix
    !! amr_fx_pre. Rows 8:11 are left to the parent-store pack site.
    impure subroutine s_amr_fx_plan(w)

        type(t_amr_wave_side), intent(in) :: w
        integer                           :: t, e, nt

        nt = w%nx
        if (allocated(amr_fx_pl)) then
            if (size(amr_fx_pl, 2) < nt) deallocate (amr_fx_pl, amr_fx_pre)
        end if
        if (.not. allocated(amr_fx_pl)) allocate (amr_fx_pl(11, max(nt, 64)), amr_fx_pre(max(nt, 64) + 1))
        e = 0
        do t = 1, nt
            amr_fx_pl(1:3,t) = w%bl(:,t)
            amr_fx_pl(4:6,t) = w%bh(:,t) - w%bl(:,t) + 1
            amr_fx_pl(7, t) = w%off(t)
            amr_fx_pre(t) = e
            e = e + sys_size*amr_fx_pl(4, t)*amr_fx_pl(5, t)*amr_fx_pl(6, t)
        end do
        amr_fx_pre(nt + 1) = e

    end subroutine s_amr_fx_plan

    !> Longest run t0:t1 of one box's recv transfers that is contiguous in the wire pool (unpacks in one launch).
    pure subroutine s_amr_fx_run(k, blk, nx, t0, t1)

        integer, intent(in)  :: k, blk(:), nx, t0
        integer, intent(out) :: t1

        t1 = t0
        do while (t1 < nx)
            if (blk(t1 + 1) /= k) exit
            if (amr_fx_pl(7, t1 + 1) /= amr_fx_pl(7, t1) + amr_fx_pre(t1 + 1) - amr_fx_pre(t1)) exit
            t1 = t1 + 1
        end do

    end subroutine s_amr_fx_run

    !> Fused coarse-field pack: the wave's whole send list in one launch; inverts the linear index s_amr_pack_box_device writes, so
    !! wire bytes match the per-transfer packs.
    impure subroutine s_amr_fx_pack_box(q_coarse, t0, t1, o1, o2, o3, pl, pre, buf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_coarse
        integer, intent(in)                                 :: t0, t1, o1, o2, o3
        integer, intent(in), contiguous                     :: pl(:,:), pre(:)
        real(wp), intent(inout), contiguous                 :: buf(:)
        integer                                             :: e, t, lo, hi, md, r, n1, n2, n3, i, g1, g2, g3

        $:GPU_PARALLEL_LOOP(copyin='[pl, pre]', copyout='[buf]', private='[t, lo, hi, md, r, n1, n2, n3, i, g1, g2, g3]')
        do e = pre(t0) + 1, pre(t1 + 1)
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
            buf(pl(7, t) + 1 + r) = real(q_coarse(i)%sf(g1 - o1, g2 - o2, g3 - o3), wp)
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_fx_pack_box

    !> Fused parent-store pack: as s_amr_fx_pack_box, with the source the flat store at parent slot pl(8,:) and child frame
    !! pl(9:11,:).
    impure subroutine s_amr_fx_pack_parent(t0, t1, pl, pre, buf)

        integer, intent(in)                 :: t0, t1
        integer, intent(in), contiguous     :: pl(:,:), pre(:)
        real(wp), intent(inout), contiguous :: buf(:)
        integer                             :: e, t, lo, hi, md, r, n1, n2, n3, i, g1, g2, g3
        real(wp)                            :: th

        th = amr_par_theta
        $:GPU_PARALLEL_LOOP(copyin='[pl, pre, th]', copyout='[buf]', private='[t, lo, hi, md, r, n1, n2, n3, i, g1, g2, g3]')
        do e = pre(t0) + 1, pre(t1 + 1)
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
            buf(pl(7, t) + 1 + r) = f_amr_par_val(g1 + pl(9, t), g2 + pl(10, t), g3 + pl(11, t), i, pl(8, t), th)
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_fx_pack_parent

    !> Fused unpack of the contiguous transfer run t0:t1 into amr_cg in one launch; serves both packs via the per-transfer frame in
    !! plan rows 9:11.
    impure subroutine s_amr_fx_unpack(t0, t1, base, pl, pre, buf, cg)

        type(scalar_field), dimension(sys_size), intent(inout) :: cg  !< amr_cg, passed so a regrown pool maps as an argument
        integer, intent(in)                                    :: t0, t1, base
        integer, intent(in), contiguous                        :: pl(:,:), pre(:)
        real(wp), intent(in), contiguous                       :: buf(:)
        integer                                                :: e, t, lo, hi, md, r, n1, n2, n3, i, g1, g2, g3

        $:GPU_PARALLEL_LOOP(copyin='[pl, pre, buf]', private='[t, lo, hi, md, r, n1, n2, n3, i, g1, g2, g3]')
        do e = pre(t0) + 1, pre(t1 + 1)
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
            cg(i)%sf(g1 - pl(9, t), g2 - pl(10, t), g3 - pl(11, t)) = real(buf(pl(7, t) - base + 1 + r), stp)
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_fx_unpack

    !> Level-1 fill wave, exchange phase: each level-1 box's padded coarse patch (region +/- amr_cpat_mar) is assembled on its owner
    !! from the coarse owners, one message per peer; both sides derive the transfers from replicated metadata. Slices are
    !! ring-clipped to the hollow shell unless full (init, regrid). Consume with s_amr_l1_fill_consume per owned box ascending.
    impure subroutine s_amr_l1_fill_exchange(q_cons_coarse, full)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_coarse
        logical, intent(in)                                    :: full
        integer                                                :: k, r, idx, ix, owner, o1, o2, o3, lo, hi, kk
        integer                                                :: plo(3), phi(3), crlo(3), crhi(3), bl(3), bh(3), msl
        integer                                                :: tb(3, 6), te(3, 6)

        amr_wcur = 1
        if (amr_num_blocks <= 0) return
        o1 = amr_sidx(1); o2 = amr_sidx(2); o3 = amr_sidx(3)
        call s_amr_wave_open(amr_wave, amr_band_stage_fill)

        call s_phase_tic(PH_GATHER)
        if (amr_seam_pairs_dirty .or. amr_seam_pairs_nblk /= amr_num_blocks) call s_amr_build_seam_pairs()
        ! the Lagrangian-overlap safety check must cover every level-1 block (mine included)
        if (bubbles_lagrange) then
            do k = 1, amr_num_blocks
                if (amr_block_level(k) /= 1) cycle
                call s_amr_select_slot(k)
                call s_amr_check_lag_clear()
            end do
        end if
        call s_amr_wave_reset(amr_wsend)
        call s_amr_refresh_lists()
        do kk = 1, amr_n_l1p
            k = amr_l1p_blk(kk)
            owner = amr_block_owner(k)
            call s_amr_patch_box(k, plo, phi)
            call s_amr_rank_coarse_range(proc_rank, crlo, crhi)
            call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
            if (any(bl > bh)) cycle
            call s_amr_l1_slice_slabs(k, plo, phi, bl, bh, full, msl, tb, te)
            call s_amr_wave_add_slabs(amr_wsend, owner, k, msl, tb, te)
        end do
        call s_amr_wave_close(amr_wsend, amr_fw_sq, amr_fw_dev)

        call s_amr_wave_reset(amr_wrecv)
        call s_amr_refresh_my_blocks()
        do kk = 1, amr_n_my
            k = amr_my_blk(kk)
            if (amr_block_level(k) /= 1) cycle
            call s_amr_patch_box(k, plo, phi)
            do idx = 1, amr_ovl_gather_n(k)
                r = amr_ovl_gather(idx, k)
                if (r == proc_rank) cycle
                call s_amr_rank_coarse_range(r, crlo, crhi)
                call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
                call s_amr_l1_slice_slabs(k, plo, phi, bl, bh, full, msl, tb, te)
                call s_amr_wave_add_slabs(amr_wrecv, r, k, msl, tb, te)
            end do
        end do
        call s_amr_wave_close(amr_wrecv, amr_fw_rq, amr_fw_dev)

        ! post recvs, pack sends, post sends, one wait
        call s_amr_wave_post(amr_wave, amr_wrecv, amr_fw_rq, amr_fw_dev)
        if (amr_device_pack .and. amr_wsend%nx > 0) then
            ! the pool is exactly tiled by the side's words; padding would ship garbage
            call s_amr_fx_plan(amr_wsend)
            call s_amr_fx_pack_box(q_cons_coarse, 1, amr_wsend%nx, o1, o2, o3, amr_fx_pl(:,1:amr_wsend%nx), &
                                   & amr_fx_pre(1:amr_wsend%nx + 1), amr_fw_sq(1:amr_wsend%words))
        else
            do ix = 1, amr_wsend%nx
                call s_amr_wave_slice(amr_wsend, ix, lo, hi)
                call s_amr_pack_box_device(q_cons_coarse, amr_wsend%bl(:,ix), amr_wsend%bh(:,ix), o1, o2, o3, amr_fw_sq(lo:hi))
            end do
        end if
        call s_amr_wave_send(amr_wave, amr_wsend, amr_fw_sq, amr_fw_dev)
        call s_amr_wave_wait(amr_wave)
        if (amr_device_pack .and. amr_wrecv%nx > 0) call s_amr_fx_plan(amr_wrecv)
        call s_phase_toc(PH_GATHER)

    end subroutine s_amr_l1_fill_exchange

    !> A contributor's slice [bl, bh] of box k's patch as transfer slabs: the whole slice (full) or its hollow-shell intersection.
    impure subroutine s_amr_l1_slice_slabs(k, plo, phi, bl, bh, full, msl, tb, te)

        integer, intent(in)  :: k, plo(3), phi(3), bl(3), bh(3)
        logical, intent(in)  :: full
        integer, intent(out) :: msl, tb(3, 6), te(3, 6)
        integer              :: clo(3), chi(3), nsh, scells, sb(3, 6), se(3, 6)

        if (full) then
            msl = 1; tb(:,1) = bl; te(:,1) = bh
            return
        end if
        call s_amr_patch_core(k, clo, chi)
        call s_amr_shell_slabs(plo, phi, clo, chi, nsh, sb, se, scells)
        call s_amr_shell_clip(nsh, sb, se, bl, bh, msl, tb, te, scells)

    end subroutine s_amr_l1_slice_slabs

    !> Drain block k's received transfers at the wave cursor into amr_cg (patch origin off): one fused unpack per contiguous run
    !! under amr_device_pack, else one unpack per transfer.
    impure subroutine s_amr_fill_drain(k, off)

        integer, intent(in) :: k, off(3)
        integer             :: lo, hi, ie, jx, boff

        do while (amr_wcur <= amr_wrecv%nx)
            if (amr_wrecv%blk(amr_wcur) /= k) exit
            if (amr_device_pack) then
                call s_amr_fx_run(k, amr_wrecv%blk, amr_wrecv%nx, amr_wcur, ie)
                boff = amr_fx_pl(7, amr_wcur)
                do jx = amr_wcur, ie
                    amr_fx_pl(9:11,jx) = off
                end do
                call s_amr_fx_unpack(amr_wcur, ie, boff, amr_fx_pl(:,1:amr_wrecv%nx), amr_fx_pre(1:amr_wrecv%nx + 1), &
                                     & amr_fw_rq(boff + 1:amr_fx_pl(7, ie) + amr_fx_pre(ie + 1) - amr_fx_pre(ie)), amr_cg)
                amr_wcur = ie + 1
                cycle
            end if
            call s_amr_wave_slice(amr_wrecv, amr_wcur, lo, hi)
            call s_amr_unpack_box_device(amr_wrecv%bl(:,amr_wcur), amr_wrecv%bh(:,amr_wcur), off, amr_fw_rq(lo:hi))
            amr_wcur = amr_wcur + 1
        end do

    end subroutine s_amr_fill_drain

    !> Level-1 fill wave, consume for owned box k: set the patch frame, copy the own slice, unpack k's transfers. Consume in
    !! ascending order: transfers are box-major, so k's are the next run at the cursor.
    impure subroutine s_amr_l1_fill_consume(q_cons_coarse, k, full)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_coarse
        integer, intent(in)                                 :: k
        logical, intent(in)                                 :: full
        integer                                             :: o1, o2, o3, msl, isl
        integer                                             :: plo(3), phi(3), crlo(3), crhi(3), bl(3), bh(3)
        integer                                             :: tb(3, 6), te(3, 6)

        call s_phase_tic(PH_GATHER)
        o1 = amr_sidx(1); o2 = amr_sidx(2); o3 = amr_sidx(3)
        call s_amr_patch_box(k, plo, phi)
        amr_cpat_off = plo
#ifdef MFC_DEBUG
        ! NaN-flood before the writes land
        call s_amr_poison_patch_device(phi(1) - plo(1), phi(2) - plo(2), phi(3) - plo(3))
#endif
        call s_amr_rank_coarse_range(proc_rank, crlo, crhi)
        call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
        call s_amr_l1_slice_slabs(k, plo, phi, bl, bh, full, msl, tb, te)
        if (msl > 0 .and. all(bl <= bh)) then
            do isl = 1, msl
                call s_amr_slab_add(tb(:,isl), te(:,isl), 0, plo, plo)
            end do
            call s_amr_gather_own_rows(q_cons_coarse, amr_cg, o1, o2, o3)
        end if
        call s_amr_fill_drain(k, amr_cpat_off)
        call s_phase_toc(PH_GATHER)

    end subroutine s_amr_l1_fill_consume

    !> Every receive transfer of the open fill wave has been consumed.
    impure subroutine s_amr_fill_wave_done()

        @:ASSERT(amr_wcur == amr_wrecv%nx + 1, "fill wave: unconsumed recv transfers")

    end subroutine s_amr_fill_wave_done

    !> Per-stage level-1 fill: the wave, then every owned level-1 box's consume and ghost fill at once.
    impure subroutine s_amr_stage_fill_wave(q_cons_coarse)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_coarse

        call s_amr_l1_fill_exchange(q_cons_coarse, .false.)
        call s_amr_l1_fill_batched(q_cons_coarse)

    end subroutine s_amr_stage_fill_wave

    !> Unpack the queued slab rows (row 9: payload offset in buf) into the pooled amr_cg patches: a whole fill wave in one launch.
    impure subroutine s_amr_unpack_rows(buf, cg)

        real(wp), intent(in), contiguous                       :: buf(:)
        type(scalar_field), dimension(sys_size), intent(inout) :: cg
        integer                                                :: i, s, g, r, n1, n2, g1, g2, g3, stot, nsl, lo, hi, md

        stot = f_amr_slab_push(nsl)
        if (stot == 0) return  ! no level-1 block on this rank: the store may be unallocated
        $:GPU_PARALLEL_LOOP(collapse=2, copyin='[buf]', private='[s, lo, hi, md, r, n1, n2, g1, g2, g3]')
        do i = 1, sys_size
            do g = 0, stot - 1
                lo = 1; hi = nsl
                do while (lo < hi)
                    md = (lo + hi + 1)/2
                    if (amr_slab_tab(7, md) <= g) then; lo = md; else; hi = md - 1; end if
                end do
                s = lo
                r = g - amr_slab_tab(7, s)
                n1 = amr_slab_tab(4, s) - amr_slab_tab(1, s) + 1; n2 = amr_slab_tab(5, s) - amr_slab_tab(2, s) + 1
                g1 = amr_slab_tab(1, s) + mod(r, n1)
                g2 = amr_slab_tab(2, s) + mod(r/n1, n2)
                g3 = amr_slab_tab(3, s) + r/(n1*n2)
                cg(i)%sf(g1 - amr_slab_tab(13, s), g2 - amr_slab_tab(14, s), g3 - amr_slab_tab(15, s)) = real(buf(amr_slab_tab(9, &
                   & s) + 1 + r + amr_slab_tab(8, s)*(i - 1)), stp)
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_unpack_rows

    !> Restriction over the queued slab rows (row 9 store column, 10:12 region lo, 13 wire payload offset): children average into
    !! the coarse field (sf, frame o) or into buf (wire layout of s_amr_unpack_rows).
    #:for SFX in ['sf', 'wire']
        impure subroutine s_amr_restrict_rows_${SFX}$(${'coarse_tgt, o1, o2, o3' if SFX == 'sf' else 'buf'}$)

            #:if SFX == 'sf'
                type(scalar_field), dimension(sys_size), intent(inout) :: coarse_tgt
                integer, intent(in)                                    :: o1, o2, o3
            #:else
                real(wp), intent(inout), contiguous :: buf(:)
            #:endif
            integer :: i, s, g, r, n1, n2, ci, cj, ck, stot, nsl, lo, hi, md, rr, fi0, fj0, fk0, ddi, ddj, ddk, dj_hi, dk_hi, nchild
            real(wp) :: acc

            stot = f_amr_slab_push(nsl)
            if (stot == 0) return
            rr = amr_ref_ratio
            dj_hi = merge(rr - 1, 0, amr_dim(2)); dk_hi = merge(rr - 1, 0, amr_dim(3)); nchild = rr**num_dims
            #:if SFX == 'wire'
                $:GPU_PARALLEL_LOOP(collapse=2, private='[s, lo, hi, md, r, n1, n2, ci, cj, ck, fi0, fj0, fk0, ddi, ddj, ddk, &
                                    & acc]', copyout='[buf]')
            #:else
                $:GPU_PARALLEL_LOOP(collapse=2, private='[s, lo, hi, md, r, n1, n2, ci, cj, ck, fi0, fj0, fk0, ddi, ddj, ddk, acc]')
            #:endif
            do i = 1, sys_size
                do g = 0, stot - 1
                    lo = 1; hi = nsl
                    do while (lo < hi)
                        md = (lo + hi + 1)/2
                        if (amr_slab_tab(7, md) <= g) then; lo = md; else; hi = md - 1; end if
                    end do
                    s = lo
                    r = g - amr_slab_tab(7, s)
                    n1 = amr_slab_tab(4, s) - amr_slab_tab(1, s) + 1; n2 = amr_slab_tab(5, s) - amr_slab_tab(2, s) + 1
                    ci = amr_slab_tab(1, s) + mod(r, n1)
                    cj = amr_slab_tab(2, s) + mod(r/n1, n2)
                    ck = amr_slab_tab(3, s) + r/(n1*n2)
                    fi0 = (ci - amr_slab_tab(10, s))*rr; fj0 = (cj - amr_slab_tab(11, s))*rr; fk0 = (ck - amr_slab_tab(12, s))*rr
                    acc = 0._wp
                    do ddk = 0, dk_hi
                        do ddj = 0, dj_hi
                            do ddi = 0, rr - 1
                                acc = acc + real(amr_cons_st(fi0 + ddi, fj0 + ddj, fk0 + ddk, i, amr_slab_tab(9, s)), wp)
                            end do
                        end do
                    end do
                    #:if SFX == 'sf'
                        coarse_tgt(i)%sf(ci - o1, cj - o2, ck - o3) = real(acc/real(nchild, wp), stp)
                    #:else
                        buf(amr_slab_tab(13, s) + 1 + r + amr_slab_tab(8, s)*(i - 1)) = acc/real(nchild, wp)
                    #:endif
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

        end subroutine s_amr_restrict_rows_${SFX}$
    #:endfor

    !> Hold np coarse patches in amr_cg stacked along dim 3 (patch p at p*(amr_cpat_hi(3) + 1)); grows only, contents do not survive
    !! a fill.
    impure subroutine s_amr_cg_reserve(np)

        integer, intent(in) :: np
        integer             :: i

        if (np <= amr_cg_np) return
        amr_cg_np = max(np, 2*amr_cg_np)
        do i = 1, sys_size
            @:ACC_TEARDOWN_SFs(amr_cg(i))
            @:DEALLOCATE(amr_cg(i)%sf)
            @:ALLOCATE(amr_cg(i)%sf(0:amr_cpat_hi(1), 0:amr_cpat_hi(2), 0:amr_cg_np*(amr_cpat_hi(3) + 1) - 1))
            amr_cg(i)%sf = 0._stp  ! finite padding for the device copy
            @:ACC_SETUP_SFs(amr_cg(i))
        end do

    end subroutine s_amr_cg_reserve

    !> Level-1 consume and ghost fill for all owned boxes at once: patches pooled in amr_cg (s_amr_cg_reserve), one kernel each.
    impure subroutine s_amr_l1_fill_batched(q_cons_coarse)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_coarse

        call s_amr_l1_pool_gather(q_cons_coarse)
        call s_amr_l1_pool_fill()

    end subroutine s_amr_l1_fill_batched

    !> Consume the open level-1 fill wave into the pooled amr_cg: box k at dim-3 offset amr_cg_pz(k), amr_cg_nuse patches.
    impure subroutine s_amr_l1_pool_gather(q_cons_coarse)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_coarse
        integer                                             :: k, kk, isl, msl, t
        integer                                             :: plo(3), phi(3), crlo(3), crhi(3), bl(3), bh(3), tb(3, 6), te(3, 6)

        call s_phase_tic(PH_GATHER)
        if (.not. allocated(amr_cg_pz)) allocate (amr_cg_pz(amr_max_blocks))
        amr_cg_nuse = 0
        do kk = 1, amr_n_my
            k = amr_my_blk(kk)
            if (amr_block_level(k) /= 1) cycle
            amr_cg_pz(k) = amr_cg_nuse*(amr_cpat_hi(3) + 1); amr_cg_nuse = amr_cg_nuse + 1
        end do
        call s_amr_cg_reserve(amr_cg_nuse)
        call s_amr_rank_coarse_range(proc_rank, crlo, crhi)
        do kk = 1, amr_n_my
            k = amr_my_blk(kk)
            if (amr_block_level(k) /= 1) cycle
            call s_amr_patch_box(k, plo, phi)
            call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
            if (any(bl > bh)) cycle
            call s_amr_l1_slice_slabs(k, plo, phi, bl, bh, .false., msl, tb, te)
            do isl = 1, msl
                call s_amr_slab_add(tb(:,isl), te(:,isl), 0, plo, plo - [0, 0, amr_cg_pz(k)])
            end do
        end do
        call s_amr_gather_own_rows(q_cons_coarse, amr_cg, amr_sidx(1), amr_sidx(2), amr_sidx(3))
        if (amr_wrecv%nx > 0) then
            do t = 1, amr_wrecv%nx
                k = amr_wrecv%blk(t)
                call s_amr_patch_box(k, plo, phi)
                call s_amr_slab_add(amr_wrecv%bl(:,t), amr_wrecv%bh(:,t), amr_wrecv%off(t), plo, plo - [0, 0, amr_cg_pz(k)])
            end do
            call s_amr_unpack_rows(amr_fw_rq(1:amr_wrecv%words), amr_cg)
        end if
        amr_wcur = amr_wrecv%nx + 1
        call s_phase_toc(PH_GATHER)

    end subroutine s_amr_l1_pool_gather

    !> Ghost-fill every owned level-1 box from its pooled patch (s_amr_l1_pool_gather's layout) in one kernel.
    impure subroutine s_amr_l1_pool_fill()

        integer :: k, kk, plo(3), phi(3)

        call s_phase_tic(PH_GFILL)
        do kk = 1, amr_n_my
            k = amr_my_blk(kk)
            if (amr_block_level(k) /= 1) cycle
            call s_amr_select_slot(k)
            call s_amr_patch_box(k, plo, phi)
            amr_cpat_off = plo
            call s_amr_queue_ghost_slabs(amr_loc_of(k), amr_cg_pz(k))
        end do
        call s_amr_fill_ghost_rows(amr_cg)
        call s_phase_toc(PH_GFILL)

    end subroutine s_amr_l1_pool_fill

    !> Padded parent-patch extents of a child with parent footprint [plo, phi] (amr_cpat_mar coarse cells each side).
    pure subroutine s_amr_patch_width(plo, phi, w1, w2, w3)

        integer, intent(in)  :: plo(3), phi(3)
        integer, intent(out) :: w1, w2, w3

        w1 = phi(1) - plo(1) + 2*amr_cpat_mar
        w2 = merge(phi(2) - plo(2) + 2*amr_cpat_mar, 0, amr_dim(2)); w3 = merge(phi(3) - plo(3) + 2*amr_cpat_mar, 0, amr_dim(3))

    end subroutine s_amr_patch_width

    !> Parent-fill transfer list of a child with padded patch (w1, w2, w3), patch-local: the whole patch (full) or its hollow shell.
    !! Send, receive and consume share it.
    impure subroutine s_amr_parent_slabs(w1, w2, w3, full, msl, tb, te)

        integer, intent(in)  :: w1, w2, w3
        logical, intent(in)  :: full
        integer, intent(out) :: msl, tb(3, 6), te(3, 6)
        integer              :: clo(3), chi(3), scells

        if (full) then
            msl = 1; tb(:,1) = 0; te(:,1) = [w1, w2, w3]
            return
        end if
        clo = merge(amr_cpat_mar + f_amr_prolong_reach(), 0, amr_dim)
        chi = merge([w1, w2, w3] - amr_cpat_mar - f_amr_prolong_reach(), 0, amr_dim)
        call s_amr_shell_slabs([0, 0, 0], [w1, w2, w3], clo, chi, msl, tb, te, scells)

    end subroutine s_amr_parent_slabs

    !> A child's patch frame: the global parent-fine origin of its padded parent patch (region footprint - amr_cpat_mar).
    impure subroutine s_amr_parent_frame(k, plo, phi, w1, w2, w3)

        integer, intent(in)  :: k
        integer, intent(out) :: plo(3), phi(3), w1, w2, w3

        call s_amr_parent_foot(k, amr_parent_blk(k), plo, phi)
        call s_amr_patch_width(plo, phi, w1, w2, w3)
        amr_cpat_off = merge(plo - amr_cpat_mar, 0, amr_dim)

    end subroutine s_amr_parent_frame

    !> Level-lev parent-fill wave, exchange phase: each block's parent patch from its parent's fine array (parent-fine frame); split
    !! children ship from the parent's owner, co-located parents copy at consume. Called per level ascending so parents are
    !! complete. Consume with s_amr_parent_fill_consume per owned box.
    impure subroutine s_amr_parent_fill_exchange(lev, full)

        integer, intent(in) :: lev
        logical, intent(in) :: full
        integer             :: k, ix, pblk, powner, cowner, lo, hi, kk
        integer             :: w1, w2, w3, plo(3), phi(3), msl
        integer             :: tb(3, 6), te(3, 6)
        integer             :: loc_parent_blk_k  !< plain scalar for device routines (NVHPC -Minline)

        amr_wcur = 1
        if (amr_num_blocks <= 0) return
        call s_amr_wave_open(amr_wave, amr_band_parent_fill)

        call s_phase_tic(PH_GATHER)
        ! the Lagrangian-overlap safety check must cover every level-lev block (mine included)
        if (bubbles_lagrange) then
            do k = 1, amr_num_blocks
                if (amr_block_level(k) /= lev) cycle
                call s_amr_select_slot(k)
                call s_amr_check_lag_clear()
            end do
        end if
        call s_amr_refresh_lists()
        ! send side: every level-lev block whose parent I own but whose owner is another rank
        call s_amr_wave_reset(amr_wsend)
        do kk = 1, amr_n_fch
            k = amr_fch_blk(kk)
            if (amr_block_level(k) /= lev) cycle
            pblk = amr_parent_blk(k)
            powner = amr_block_owner(pblk); cowner = amr_block_owner(k)
            if (powner == cowner .or. powner /= proc_rank) cycle
            call s_amr_parent_frame(k, plo, phi, w1, w2, w3)
            call s_amr_parent_slabs(w1, w2, w3, full, msl, tb, te)
            call s_amr_wave_add_slabs(amr_wsend, cowner, k, msl, tb, te)
        end do
        call s_amr_wave_close(amr_wsend, amr_fw_sq, amr_fw_dev)
        ! receive side: every level-lev block I own whose parent lives on another rank
        call s_amr_wave_reset(amr_wrecv)
        call s_amr_refresh_my_blocks()
        do kk = 1, amr_n_my
            k = amr_my_blk(kk)
            if (amr_block_level(k) /= lev) cycle
            pblk = amr_parent_blk(k)
            powner = amr_block_owner(pblk)
            if (powner == proc_rank) cycle
            call s_amr_parent_frame(k, plo, phi, w1, w2, w3)
            call s_amr_parent_slabs(w1, w2, w3, full, msl, tb, te)
            call s_amr_wave_add_slabs(amr_wrecv, powner, k, msl, tb, te)
        end do
        call s_amr_wave_close(amr_wrecv, amr_fw_rq, amr_fw_dev)
        call s_amr_wave_post(amr_wave, amr_wrecv, amr_fw_rq, amr_fw_dev)
        ! packs: one fused launch over the whole send list under amr_device_pack, else one launch per transfer
        if (amr_device_pack .and. amr_wsend%nx > 0) then
            call s_amr_fx_plan(amr_wsend)
            do ix = 1, amr_wsend%nx
                k = amr_wsend%blk(ix)
                call s_amr_parent_frame(k, plo, phi, w1, w2, w3)
                amr_fx_pl(8, ix) = amr_loc_of(amr_parent_blk(k))
                amr_fx_pl(9:11,ix) = amr_cpat_off
            end do
            call s_amr_fx_pack_parent(1, amr_wsend%nx, amr_fx_pl(:,1:amr_wsend%nx), amr_fx_pre(1:amr_wsend%nx + 1), &
                                      & amr_fw_sq(1:amr_wsend%words))
        else
            do ix = 1, amr_wsend%nx
                k = amr_wsend%blk(ix)
                call s_amr_parent_frame(k, plo, phi, w1, w2, w3)
                call s_amr_wave_slice(amr_wsend, ix, lo, hi)
                loc_parent_blk_k = amr_loc_of(amr_parent_blk(k))
                call s_amr_pack_parent_box_device(loc_parent_blk_k, amr_wsend%bl(:,ix), amr_wsend%bh(:,ix), amr_fw_sq(lo:hi))
            end do
        end if
        call s_amr_wave_send(amr_wave, amr_wsend, amr_fw_sq, amr_fw_dev)
        call s_amr_wave_wait(amr_wave)
        if (amr_device_pack .and. amr_wrecv%nx > 0) call s_amr_fx_plan(amr_wrecv)
        call s_phase_toc(PH_GATHER)

    end subroutine s_amr_parent_fill_exchange

    !> The parent-fill wave, consume phase for owned box k (the current slot): set the patch frame, then a co-located parent's
    !! own-box device copies or k's received transfers (box-major: the next contiguous run at the cursor).
    impure subroutine s_amr_parent_fill_consume(k, full)

        integer, intent(in) :: k
        logical, intent(in) :: full
        integer             :: pblk, isl, msl
        integer             :: w1, w2, w3, plo(3), phi(3)
        integer             :: tb(3, 6), te(3, 6)
        integer             :: loc_pblk  !< plain scalar for device routines (NVHPC -Minline)

        call s_phase_tic(PH_GATHER)
        pblk = amr_parent_blk(k)
        call s_amr_parent_frame(k, plo, phi, w1, w2, w3)
#ifdef MFC_DEBUG
        call s_amr_poison_patch_device(w1, w2, w3)
#endif
        if (amr_block_owner(pblk) == proc_rank) then
            call s_amr_parent_slabs(w1, w2, w3, full, msl, tb, te)
            do isl = 1, msl
                loc_pblk = amr_loc_of(pblk)
                call s_amr_copy_parent_box(loc_pblk, tb(:,isl), te(:,isl))
            end do
        else
            @:ASSERT(amr_wcur <= amr_wrecv%nx .and. amr_wrecv%blk(amr_wcur) == k, "parent-fill wave: missing recv transfer")
            call s_amr_fill_drain(k, [0, 0, 0])
        end if
        call s_phase_toc(PH_GATHER)

    end subroutine s_amr_parent_fill_consume

    !> Per-stage level-lev fill: the wave, then per owned level-lev box the consume and the ghost fill.
    impure subroutine s_amr_parent_fill_wave(lev)

        integer, intent(in) :: lev
        integer             :: k, kk
        integer             :: loc_cur  !< plain scalar for device routines (NVHPC -Minline)

        call s_amr_parent_fill_exchange(lev, .false.)
        do kk = 1, amr_n_own
            k = amr_own_blk(kk)
            if (amr_block_level(k) /= lev) cycle
            call s_amr_select_slot(k)
            if (.not. amr_rank_owns_block) cycle
            call s_amr_parent_fill_consume(k, .false.)
            call s_phase_tic(PH_GFILL)
            loc_cur = amr_loc_of(amr_cur)
            call s_amr_fill_fine_ghosts(amr_cg, loc_cur)
            call s_phase_toc(PH_GFILL)
        end do
        call s_amr_fill_wave_done()

    end subroutine s_amr_parent_fill_wave

    !> Device pack (to_buf=T) / unpack (F) of the interior block [o+0:o+fm] of the flat store at column loc <-> contiguous MPI
    !! buffer buf for P2P migration; wp buffer cast to/from stp. Only buf is mapped (flang miscomputes strided %sf sections in map
    !! clauses).
    impure subroutine s_amr_pack_unpack_block_st(loc, o1, o2, o3, fm1, fm2, fm3, buf, to_buf)

        integer, intent(in)                 :: loc, o1, o2, o3, fm1, fm2, fm3
        real(wp), intent(inout), contiguous :: buf(:)
        logical, intent(in)                 :: to_buf
        integer                             :: i, j, k, l

        #:for TB in [True, False]
            #:if TB
                if (to_buf) then
                    $:GPU_PARALLEL_LOOP(collapse=4, copyout='[buf]')
                #:else
                else
                    $:GPU_PARALLEL_LOOP(collapse=4, copyin='[buf]')
                #:endif
                do i = 1, sys_size
                    do l = 0, fm3
                        do k = 0, fm2
                            do j = 0, fm1
                                #:if TB
                                    buf(1 + j + (fm1 + 1)*(k + (fm2 + 1)*(l + (fm3 + 1)*(i - 1)))) = real(amr_cons_st(o1 + j, &
                                        & o2 + k, o3 + l, i, loc), wp)
                                #:else
                                    amr_cons_st(o1 + j, o2 + k, o3 + l, i, &
                                                & loc) = real(buf(1 + j + (fm1 + 1)*(k + (fm2 + 1)*(l + (fm3 + 1)*(i - 1)))), stp)
                                #:endif
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            #:endfor
        end if

    end subroutine s_amr_pack_unpack_block_st

end module m_amr_exchange
