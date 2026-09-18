!>
!!@file
!!@brief Contains module m_amr_exchange

#! AMD OpenMP lane: assert allocatables present on every kernel here (see OMP_DEFAULT_STR). A conditionally allocated module
#! array a kernel names launches only under its allocation's own condition (sw_jac/jac: igr); a kernel naming an unallocated
#! array aborts. Keep it so.
#:set MFC_OMP_PRESENT_ALLOCATABLE = True
#:include 'macros.fpp'

!> @brief Parent/child and fine-fine data exchange: gather plans, pack/unpack, seams, ghost fills and the fill waves.
module m_amr_exchange

#ifdef MFC_MPI
    use mpi  !< MPI-IO for the parallel_io AMR restart file
#endif

    use m_derived_types  ! scalar_field, t_box, int_bounds_info
    use m_box, only: f_morton
    use m_global_parameters
    use m_mpi_proxy, only: s_mpi_abort  ! @:ASSERT expands to it
    use m_mpi_common, only: s_mpi_sendrecv_variables_buffers
    use m_amr_registers, only: s_amr_parent_foot
    use m_phase_timing
    use m_amr_xchg_audit  ! per-call-site accounting of every AMR p2p transfer (s_xa_rec + XA_* site ids)
    use m_amr_state
    use m_amr_wave
    use m_amr_distribution

    implicit none

    private
    public :: f_amr_seam, f_amr_seam_dim, s_amr_build_seam_pairs, s_amr_exchange_coarse_cons_halo, s_amr_fine_fine_drain, &
        & s_amr_fine_fine_halo, s_amr_fine_fine_post, s_amr_fill_wave_done, s_amr_l1_fill_exchange, s_amr_l1_fill_consume, &
        & s_amr_parent_fill_exchange, s_amr_parent_fill_consume, s_amr_parent_fill_wave, s_amr_stage_fill_wave, &
        & s_l0_pack_unpack_block_sf, s_l0_pack_unpack_block_st

contains

    !> Sub-box variants of the parent-patch pack/unpack/copy for the ring-clipped parent-fill wave (the wave ships q_cons). Bounds
    !! are patch-local cell ranges; the buffer holds the sub-box in the same (g1 fastest, sys_size outermost) layout as the
    !! full-patch kernels, so both wire sides agree by construction.
    impure subroutine s_amr_pack_parent_box_device(qp, bl, bh, buf)

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

    end subroutine s_amr_pack_parent_box_device

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

    impure subroutine s_amr_copy_parent_box(qp, bl, bh)

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

    end subroutine s_amr_copy_parent_box

    !> Runtime device pack of the overlap box [bl:bh] global from q_coarse (device) into the contiguous wire buffer buf (host, via
    !! copyout); only the box crosses PCIe, not the full field. Explicit-loop linear buf indexing (g1 fastest, then g2, g3, i) and
    !! the wp cast match the host pack in s_amr_gather_coarse_patch element-for-element, so the receiver's unpack is layout- and
    !! byte-identical (same discipline as s_amr_fine_slice: no array-section syntax near the device map).
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
    !! s_amr_gather_coarse_patch.
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

    !> Decompose the current fine block's ghost shell (buffered extent minus interior) into ns disjoint face slabs whose union is
    !! exactly the non-interior cells, so the ghost-fill kernels do O(surface) work instead of masking the full buffered volume. x
    !! slabs span the full transverse extent; y slabs restrict x to the interior; z slabs restrict x and y. Collapsed dims
    !! contribute no slabs.
    pure subroutine s_amr_build_ghost_slabs(ns, sb, se)

        integer, intent(out) :: ns, sb(3, 6), se(3, 6)
        integer              :: fm(3), b(3), e(3)

        fm = [amr_slots(amr_cur)%m, amr_slots(amr_cur)%n, amr_slots(amr_cur)%p]
        b = amr_slots(amr_cur)%idwbuff%beg; e = amr_slots(amr_cur)%idwbuff%end
        ns = 2*num_dims
        sb(:,1) = b; se(:,1) = [-1, e(2), e(3)]
        sb(:,2) = [fm(1) + 1, b(2), b(3)]; se(:,2) = e
        sb(:,3) = [0, b(2), b(3)]; se(:,3) = [fm(1), -1, e(3)]
        sb(:,4) = [0, fm(2) + 1, b(3)]; se(:,4) = [fm(1), e(2), e(3)]
        sb(:,5) = [0, 0, b(3)]; se(:,5) = [fm(1), fm(2), -1]
        sb(:,6) = [0, 0, fm(3) + 1]; se(:,6) = [fm(1), fm(2), e(3)]

    end subroutine s_amr_build_ghost_slabs

    !> Load the first ns slabs [sb:se] into the device slab table (rows 1:3 lo, 4:6 hi, 7 flat offset, 8 cell count) for the fused
    !! flat-index kernels, returning the total cell count.
    impure integer function f_amr_slab_tab_load(ns, sb, se) result(stot)

        integer, intent(in) :: ns, sb(3, 6), se(3, 6)
        integer             :: s

        stot = 0
        do s = 1, ns
            amr_slab_tab(1:3,s) = sb(:,s); amr_slab_tab(4:6,s) = se(:,s)
            amr_slab_tab(7, s) = stot; amr_slab_tab(8, s) = product(se(:,s) - sb(:,s) + 1)
            stot = stot + amr_slab_tab(8, s)
        end do
        $:GPU_UPDATE(device='[amr_slab_tab]')

    end function f_amr_slab_tab_load

    !> Fill the fine ghost shell by conservative-linear prolongation from q_coarse, the gathered block-local coarse patch amr_cg
    !! (fine-level distribution; the caller gathers the source first). Device kernel: reads the patch and writes the fine target in
    !! device memory. floor/modulo mapping is valid for negative fine indices (ghosts). Interior untouched. Multi-fluid volume
    !! fractions get the same sum-preserving closure as the interior prolongation (s_amr_fill_fine_ghosts_alphas). Writes the
    !! conserved store at dense local index `loc`.
    impure subroutine s_amr_fill_fine_ghosts(q_coarse, loc)

        type(scalar_field), dimension(sys_size), intent(in) :: q_coarse
        integer, intent(in)                                 :: loc
        integer                                             :: i, fi, fj, fk, ci, cj, ck, ox, oy, oz
        integer                                             :: rr, lo1, lo2, lo3
        integer                                             :: advb, adve
        integer                                             :: s, ns
        integer                                             :: ss, g, r, n1, n2, stot, sb(3, 6), se(3, 6)
        logical                                             :: d2, d3, multi
        real(wp)                                            :: u0, sx, sy, sz, xix, xiy, xiz

        ! q_coarse is the gathered block-local patch amr_cg (fine-level distribution); amr_isect_lo (global, == region_lo on
        ! the owner) + f/rr - amr_cpat_off is the patch-local coarse index. Fine indices are local to this block.

        ox = amr_cpat_off(1); oy = amr_cpat_off(2); oz = amr_cpat_off(3)
        d2 = n_glb > 0; d3 = p_glb > 0
        rr = amr_slots(amr_cur)%amr_ref_ratio
        lo1 = amr_isect_lo(1); lo2 = amr_isect_lo(2); lo3 = amr_isect_lo(3)
        multi = num_fluids > 1 .and. (.not. bubbles_lagrange)  ! EL alphas sum to beta, not 1: no sum-to-one closure
        advb = eqn_idx%adv%beg; adve = eqn_idx%adv%end
        call s_amr_build_ghost_slabs(ns, sb, se)
        ! One kernel over the concatenation of the ns face slabs instead of one kernel each. The slabs are disjoint and their union
        ! is exactly the ghost shell (s_amr_build_ghost_slabs), so every ghost cell is written exactly once and the result is
        ! independent of how the flat index is ordered. Not the padded-hull form of s_amr_capture_batch: the x slabs
        ! span the full transverse extent, so a hull over all slabs is the whole buffered volume and masking it would throw away
        ! the O(surface) decomposition this routine exists to get.
        stot = f_amr_slab_tab_load(ns, sb, se)
        $:GPU_PARALLEL_LOOP(collapse=2, private='[s, ss, r, n1, n2, fi, fj, fk, ci, cj, ck, xix, xiy, xiz, u0, sx, sy, sz]')
        do i = 1, sys_size
            do g = 0, stot - 1
                s = 1  ! decode the flat index: ns <= 6, so a scan beats storing a per-cell slab map
                do ss = 2, ns
                    if (g >= amr_slab_tab(7, ss)) s = ss
                end do
                r = g - amr_slab_tab(7, s)
                n1 = amr_slab_tab(4, s) - amr_slab_tab(1, s) + 1; n2 = amr_slab_tab(5, s) - amr_slab_tab(2, s) + 1
                fi = amr_slab_tab(1, s) + mod(r, n1)
                fj = amr_slab_tab(2, s) + mod(r/n1, n2)
                fk = amr_slab_tab(3, s) + r/(n1*n2)
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
                    if (d2) sy = minmod(real(q_coarse(i)%sf(ci, cj + 1, ck), wp) - u0, u0 - real(q_coarse(i)%sf(ci, cj - 1, ck), &
                        & wp))
                    sz = 0._wp
                    if (d3) sz = minmod(real(q_coarse(i)%sf(ci, cj, ck + 1), wp) - u0, u0 - real(q_coarse(i)%sf(ci, cj, ck - 1), &
                        & wp))
                    amr_cons_st(fi, fj, fk, i, loc) = u0 + sx*xix + sy*xiy + sz*xiz
                end if
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! multi-fluid volume-fraction ghosts: per-cell closure mirroring s_prolong_alphas_closure (own routine: nvfortran
        ! 25.x crashes on two target regions sharing this routine's privates)
        if (multi) call s_amr_fill_fine_ghosts_alphas(q_coarse, loc, stot, advb, adve)

    end subroutine s_amr_fill_fine_ghosts

    !> Volume-fraction ghosts for the fill above: shared limiter switch over all fluids; interpolate + clamp fluids advb..adve-1;
    !! alpha_n = 1 - sum. Same flat-index fusion over the same disjoint slabs (amr_slab_tab, already on device).
    impure subroutine s_amr_fill_fine_ghosts_alphas(q_coarse, loc, stot, advb, adve)

        type(scalar_field), dimension(sys_size), intent(in) :: q_coarse
        integer, intent(in)                                 :: loc, stot, advb, adve
        integer                                             :: i, fi, fj, fk, ci, cj, ck, ox, oy, oz
        integer                                             :: rr, lo1, lo2, lo3, s, ss, g, r, n1, n2, ns
        logical                                             :: d2, d3, shx, shy, shz
        real(wp)                                            :: u0, sx, sy, sz, xix, xiy, xiz, av, asum

        ox = amr_cpat_off(1); oy = amr_cpat_off(2); oz = amr_cpat_off(3)
        d2 = n_glb > 0; d3 = p_glb > 0
        rr = amr_slots(amr_cur)%amr_ref_ratio
        lo1 = amr_isect_lo(1); lo2 = amr_isect_lo(2); lo3 = amr_isect_lo(3)
        ns = 2; if (d2) ns = 4; if (d3) ns = 6
        $:GPU_PARALLEL_LOOP(private='[s, ss, r, n1, n2, fi, fj, fk, i, ci, cj, ck, xix, xiy, xiz, u0, sx, sy, sz, av, asum, shx, &
                            & shy, shz]')
        do g = 0, stot - 1
            s = 1
            do ss = 2, ns
                if (g >= amr_slab_tab(7, ss)) s = ss
            end do
            r = g - amr_slab_tab(7, s)
            n1 = amr_slab_tab(4, s) - amr_slab_tab(1, s) + 1; n2 = amr_slab_tab(5, s) - amr_slab_tab(2, s) + 1
            fi = amr_slab_tab(1, s) + mod(r, n1)
            fj = amr_slab_tab(2, s) + mod(r/n1, n2)
            fk = amr_slab_tab(3, s) + r/(n1*n2)
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
                if (shx) sx = minmod(real(q_coarse(i)%sf(ci + 1, cj, ck), wp) - u0, u0 - real(q_coarse(i)%sf(ci - 1, cj, ck), wp))
                sy = 0._wp
                if (shy) sy = minmod(real(q_coarse(i)%sf(ci, cj + 1, ck), wp) - u0, u0 - real(q_coarse(i)%sf(ci, cj - 1, ck), wp))
                sz = 0._wp
                if (shz) sz = minmod(real(q_coarse(i)%sf(ci, cj, ck + 1), wp) - u0, u0 - real(q_coarse(i)%sf(ci, cj, ck - 1), wp))
                av = min(max(u0 + sx*xix + sy*xiy + sz*xiz, 0._wp), 1._wp)
                amr_cons_st(fi, fj, fk, i, loc) = av
                asum = asum + av
            end do
            amr_cons_st(fi, fj, fk, adve, loc) = 1._wp - asum
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_fill_fine_ghosts_alphas

    !> Exchange the coarse conservative ghost layers at internal rank boundaries (physical-boundary ghosts untouched; per direction
    !! beg then end, mirroring s_populate_variables_buffers' disblock). The solver never fills cons ghosts (only prim), so ranks
    !! whose fine ghost-fill or prolongation stencil leaves their interior need this first. All ranks must call together (pairwise
    !! exchange per internal neighbor).
    impure subroutine s_amr_exchange_coarse_cons_halo(q_cons)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons

        #:for D, X in [(1, 'x'), (2, 'y'), (3, 'z')]
            if (amr_dim(${D}$)) then
                if (bc_${X}$%beg >= 0) call s_mpi_sendrecv_variables_buffers(q_cons, ${D}$, -1, sys_size)
                if (bc_${X}$%end >= 0) call s_mpi_sendrecv_variables_buffers(q_cons, ${D}$, 1, sys_size)
            end if
        #:endfor

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

    !> Rebuild the cached same-level seam-pair list (amr_seam_pairs) once per regrid/restart, in the same (xb, yb) order on all
    !! ranks so the paired seam transfers stay matched, and the per-block overlap-rank lists (amr_ovl_gather/scatter) by O(overlap)
    !! inversion of the decomposition (s_amr_ranks_overlapping).
    impure subroutine s_amr_build_seam_pairs()

        integer                      :: xb, d, np, k, mx, pass, nm, im, jm, tb, td, nb
        integer                      :: plo(3), phi(3), rlo(3), rhi(3)
        integer                      :: mb(3), md(3), clo(3), gc
        integer(kind=8), allocatable :: mkey(:)
        integer, allocatable         :: ord(:)

        if (allocated(amr_seam_pairs)) deallocate (amr_seam_pairs)

        ! Blocks are disjoint, so (level, region lo) names one uniquely, and the seam predicate fixes the
        ! neighbour's lo corner: transverse lo equal to xb's, and lo(d) = hi(d, xb) + 1. The all-pairs
        ! O(nblocks^2) scan is therefore a lookup. Morton-sort the lo corners once, then binary-search the
        ! single candidate per (block, dim) and verify the predicate on it. Emission order is xb ascending,
        ! yb ascending within xb, which the paired seam transfers depend on; a reordered list mismatches
        ! sends to receives and deadlocks.
        nb = max(amr_num_blocks, 1)
        allocate (mkey(nb), ord(nb))
        do k = 1, amr_num_blocks
            mkey(k) = f_morton(amr_region_lo_all(1, k), amr_region_lo_all(2, k), amr_region_lo_all(3, k))
        end do

        call s_amr_sort_by_key(mkey, amr_num_blocks, ord)

        ! pass 1 counts, pass 2 fills: keeps amr_seam_pairs exactly sized
        do pass = 1, 2
            np = 0
            do xb = 1, amr_num_blocks
                nm = 0
                do d = 1, num_dims
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
        deallocate (mkey, ord)
        ! per-block P2P overlap-rank lists by O(overlap) inversion (gather: rank coarse range vs the amr_cpat_mar-padded patch box;
        ! scatter: rank interior vs the region box), rank-ascending so iterating a list gives the same MPI send/recv order as a
        ! 0..num_procs-1 scan. The clamped interior-frame coord range reproduces both frames (see s_amr_coord_range). Bounded
        ! first dim = max overlap over all blocks (dealloc-realloc each build, like amr_seam_pairs above), not num_procs: a block
        ! spans O(1) ranks. Every consumer runs behind a build_seam_pairs guard, so the arrays are always sized before they are
        ! read.
        mx = 1
        do k = 1, amr_num_blocks
            call s_amr_patch_box(k, plo, phi); call s_amr_region_box(k, rlo, rhi)
            mx = max(mx, f_amr_overlap_count(plo, phi), f_amr_overlap_count(rlo, rhi))
        end do
        if (allocated(amr_ovl_gather)) deallocate (amr_ovl_gather)
        if (allocated(amr_ovl_scatter)) deallocate (amr_ovl_scatter)
        allocate (amr_ovl_gather(mx, amr_max_blocks), amr_ovl_scatter(mx, amr_max_blocks))
        do k = 1, amr_num_blocks
            call s_amr_patch_box(k, plo, phi); call s_amr_region_box(k, rlo, rhi)
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

        integer :: xb, yb, d, rX, rY, cnt, xm(3), tsz, fmul, idx, lo, hi
        integer :: r, sblk, sdlo, sdhi, ublk, udlo, udhi, eblk, edlo, edhi

        amr_sw_nsame = 0
        if (.not. amr .and. l0_ntile == 0) return
        if (amr_num_blocks < 2) return

        if (amr_seam_pairs_dirty .or. amr_seam_pairs_nblk /= amr_num_blocks) call s_amr_build_seam_pairs()
        if (allocated(amr_sw_plx)) deallocate (amr_sw_plx, amr_sw_ply, amr_sw_pd, amr_sw_pxhi, amr_sw_pfm)
        allocate (amr_sw_plx(amr_num_seam_pairs), amr_sw_ply(amr_num_seam_pairs), amr_sw_pd(amr_num_seam_pairs), &
                  & amr_sw_pxhi(amr_num_seam_pairs), amr_sw_pfm(3, amr_num_seam_pairs))
        ! Seam wave (plan-based exchange): the cross-rank pairs are one aggregated message per (peer, direction). Every pair
        ! contributes one send and one recv transfer on each of its two owners; both ranks walk the same replicated pair list
        ! ascending, so the wire layout agrees with no metadata exchange. Wire is wp with the same stp cast on unpack; under
        ! MFC_DEBUG each slab carries the identity header [site, sending slot, (d, dlo, dhi), (cnt, 0, 0)]. Same-rank pairs
        ! use the batched kernel at drain. The wave stays open across the fill waves, so it has its own handle and sides.
        call s_amr_wave_open(amr_wave_seam, 5)
        ! transfer records: blk = the slot the data belongs to, bl = (d, dlo, dhi) its slab, bh = (peer slot, its dlo, dhi)
        call s_amr_wave_reset(amr_wseam_snd); call s_amr_wave_reset(amr_wseam_rcv)
        do idx = 1, amr_num_seam_pairs
            xb = amr_seam_pairs(1, idx); yb = amr_seam_pairs(2, idx); d = amr_seam_pairs(3, idx)
            rX = amr_block_owner(xb); rY = amr_block_owner(yb)
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
            if (rX == rY) then  ! same rank owns both: defer to the one batched kernel at drain (no host buffer, no per-pair launch)
                amr_sw_nsame = amr_sw_nsame + 1
                amr_sw_plx(amr_sw_nsame) = amr_loc_of(xb); amr_sw_ply(amr_sw_nsame) = amr_loc_of(yb)
                amr_sw_pd(amr_sw_nsame) = d; amr_sw_pxhi(amr_sw_nsame) = xm(d); amr_sw_pfm(:,amr_sw_nsame) = xm
                cycle
            end if
            if (proc_rank == rX) then
                r = rY; sblk = xb; sdlo = xm(d) - buff_size + 1; sdhi = xm(d)
                ublk = xb; udlo = xm(d) + 1; udhi = xm(d) + buff_size
                eblk = yb; edlo = 0; edhi = buff_size - 1
            else
                r = rX; sblk = yb; sdlo = 0; sdhi = buff_size - 1
                ublk = yb; udlo = -buff_size; udhi = -1
                eblk = xb; edlo = xm(d) - buff_size + 1; edhi = xm(d)
            end if
            call s_amr_wave_add(amr_wseam_snd, r, sblk, [d, sdlo, sdhi], [cnt, 0, 0], cnt)
            call s_amr_wave_add(amr_wseam_rcv, r, ublk, [d, udlo, udhi], [eblk, edlo, edhi], cnt)
        end do
        call s_amr_wave_close(amr_wseam_snd, amr_sw_sq, .false.)
        call s_amr_wave_close(amr_wseam_rcv, amr_sw_rq, .false.)
        call s_amr_wave_post(amr_wave_seam, amr_wseam_rcv, amr_sw_rq, XA_F6W_RCV, .false.)
        do idx = 1, amr_wseam_snd%nx
            call s_amr_wave_slice(amr_wseam_snd, idx, lo, hi)
            call s_amr_fine_slice(amr_wseam_snd%blk(idx), amr_wseam_snd%bl(1, idx), amr_wseam_snd%bl(2, idx), amr_wseam_snd%bl(3, &
                                  & idx), amr_sw_sq(lo:hi), 1)
            call s_amr_wave_hdr_pack(amr_wseam_snd, amr_sw_sq, idx, XA_F6W_SND)
        end do
        call s_amr_wave_send(amr_wave_seam, amr_wseam_snd, amr_sw_sq, XA_F6W_SND, .false.)

    end subroutine s_amr_fine_fine_post

    impure subroutine s_amr_fine_fine_drain()

        integer :: idx, lo, hi

        if (f_amr_wave_nreq(amr_wave_seam) == 0 .and. amr_sw_nsame == 0) return
        call s_amr_wave_wait(amr_wave_seam)
        do idx = 1, amr_wseam_rcv%nx
            call s_amr_wave_slice(amr_wseam_rcv, idx, lo, hi)
            ! the header names the SENDER's slot and slab, which this side recorded in bh
            if (XA_NH > 0) call s_xa_hdr_check(amr_sw_rq(amr_wseam_rcv%off(idx) + 1:amr_wseam_rcv%off(idx) + XA_NH), XA_F6W_SND, &
                & amr_wseam_rcv%bh(1, idx), [amr_wseam_rcv%bl(1, idx), amr_wseam_rcv%bh(2, idx), amr_wseam_rcv%bh(3, idx)], &
                & [amr_wseam_rcv%cnt(idx), 0, 0])
            call s_amr_fine_slice(amr_wseam_rcv%blk(idx), amr_wseam_rcv%bl(1, idx), amr_wseam_rcv%bl(2, idx), amr_wseam_rcv%bl(3, &
                                  & idx), amr_sw_rq(lo:hi), -1)
        end do
        if (amr_sw_nsame > 0) call s_amr_fine_seam_exchange(amr_sw_nsame, amr_sw_plx(1:amr_sw_nsame), amr_sw_ply(1:amr_sw_nsame), &
            & amr_sw_pd(1:amr_sw_nsame), amr_sw_pxhi(1:amr_sw_nsame), amr_sw_pfm(:,1:amr_sw_nsame), buff_size)
        call s_amr_select_slot(1)

    end subroutine s_amr_fine_fine_drain

    !> The seam exchange as one call (post + drain), for the L0 tile stage.
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
        ! the slabs must tile the shell exactly: pairwise disjoint, cells summing to patch - core
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

    !> Intersect the shell-slab list with box [bl:bh]: the surviving clipped slabs in the same fixed order (each exchange side
    !! derives an identical list from replicated data, so empties drop symmetrically) plus their total cell count.
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
    impure subroutine s_amr_gather_own_shell_device(q_coarse, ms, tb, te, o1, o2, o3)

        type(scalar_field), dimension(sys_size), intent(in) :: q_coarse
        integer, intent(in)                                 :: ms, tb(3, 6), te(3, 6), o1, o2, o3
        integer                                             :: i, s, ss, g, r, n1, n2, g1, g2, g3, stot, coff1, coff2, coff3

        ! scalar copies: no host array may be referenced inside the device region (nvfortran/Cray demand it present)

        coff1 = amr_cpat_off(1); coff2 = amr_cpat_off(2); coff3 = amr_cpat_off(3)
        stot = f_amr_slab_tab_load(ms, tb, te)
        $:GPU_PARALLEL_LOOP(collapse=2, private='[s, ss, r, n1, n2, g1, g2, g3]')
        do i = 1, sys_size
            do g = 0, stot - 1
                s = 1  ! decode the flat index: ms <= 6, so a scan beats storing a per-cell slab map
                do ss = 2, ms
                    if (g >= amr_slab_tab(7, ss)) s = ss
                end do
                r = g - amr_slab_tab(7, s)
                n1 = amr_slab_tab(4, s) - amr_slab_tab(1, s) + 1; n2 = amr_slab_tab(5, s) - amr_slab_tab(2, s) + 1
                g1 = amr_slab_tab(1, s) + mod(r, n1)
                g2 = amr_slab_tab(2, s) + mod(r/n1, n2)
                g3 = amr_slab_tab(3, s) + r/(n1*n2)
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
    !! absolute payload offset in the wire pool (the side's transfer offsets), and the exclusive element prefix amr_fx_pre
    !! (amr_fx_pre(t+1) - amr_fx_pre(t) is transfer t's word count). Rows 8:11 are left to the F2 pack site, the only caller whose
    !! source varies per transfer.
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
            amr_fx_pl(7, t) = w%off(t) + XA_NH
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
    !! pl(8,:) in the per-transfer child patch frame pl(9:11,:): the s_amr_pack_parent_box_device body, one launch for the whole
    !! send list.
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

    !> The level-1 fill wave, exchange phase: every level-1 box's padded coarse patch (region +/- amr_cpat_mar, the reach of the
    !! prolongation and ghost-fill stencils) is assembled on its owner from the coarse owners that hold its cells, one aggregated
    !! message per peer. Send side: for every level-1 box someone else owns, my coarse-range slice of its patch box; receive side:
    !! for every level-1 box I own, each listed contributor's slice (the own slice is a device copy at consume). Both sides derive
    !! the transfer list from replicated metadata, so the wire layout needs no handshake. The per-stage ghost fill reads only the
    !! patch's hollow shell, so by default each slice is ring-clipped to up to 6 shell sub-slabs; full = .true. (init and regrid,
    !! which prolong the whole block) ships the whole slice. Consume with s_amr_l1_fill_consume per owned level-1 box ascending.
    impure subroutine s_amr_l1_fill_exchange(q_cons_coarse, full)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_coarse
        logical, intent(in)                                    :: full
        integer                                                :: k, r, idx, ix, owner, o1, o2, o3, lo, hi, kk
        integer                                                :: plo(3), phi(3), crlo(3), crhi(3), bl(3), bh(3), msl
        integer                                                :: tb(3, 6), te(3, 6)

        amr_wcur = 1
        if (amr_num_blocks <= 0) return
        o1 = amr_sidx(1); o2 = amr_sidx(2); o3 = amr_sidx(3)
        call s_amr_wave_open(amr_wave, 3)

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

        ! post all recvs, pack all sends (device kernels into contiguous pool slices), post all sends, one wait. [amr-xa] records
        ! payload words only, so the family totals are independent of the message aggregation.
        call s_amr_wave_post(amr_wave, amr_wrecv, amr_fw_rq, XA_F1W_RCV, amr_fw_dev)
        if (amr_device_pack .and. amr_wsend%nx > 0) then
            ! one launch for the whole send list; the debug identity headers are written after it because the fused copyout
            ! covers the pool prefix (payload and header words). The pool is exactly tiled by the side's words, so every word
            ! the copyout returns was written; any padding would ship garbage on the wire silently.
            call s_amr_fx_plan(amr_wsend)
            call s_amr_fx_pack_box(q_cons_coarse, 1, amr_wsend%nx, o1, o2, o3, amr_fx_pl(:,1:amr_wsend%nx), &
                                   & amr_fx_pre(1:amr_wsend%nx + 1), amr_fw_sq(1:amr_wsend%words))
            do ix = 1, amr_wsend%nx
                call s_amr_wave_hdr_pack(amr_wsend, amr_fw_sq, ix, XA_F1W_SND)
            end do
        else
            do ix = 1, amr_wsend%nx
                call s_amr_wave_slice(amr_wsend, ix, lo, hi)
                call s_amr_wave_hdr_pack(amr_wsend, amr_fw_sq, ix, XA_F1W_SND)
                call s_amr_pack_box_device(q_cons_coarse, amr_wsend%bl(:,ix), amr_wsend%bh(:,ix), o1, o2, o3, amr_fw_sq(lo:hi))
            end do
        end if
        call s_amr_wave_send(amr_wave, amr_wsend, amr_fw_sq, XA_F1W_SND, amr_fw_dev)
        call s_amr_wave_wait(amr_wave)
        if (amr_device_pack .and. amr_wrecv%nx > 0) call s_amr_fx_plan(amr_wrecv)
        call s_phase_toc(PH_GATHER)

    end subroutine s_amr_l1_fill_exchange

    !> A contributor's slice [bl, bh] of level-1 box k's patch as the wave's transfer slabs: the whole slice (full) or its
    !! intersection with the patch's hollow shell.
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

    !> The level-1 fill wave, consume phase for owned box k (the current slot): set the patch frame, device-copy the own slice, then
    !! unpack k's received transfers. Boxes must be consumed in the order the receive side was planned (owned, ascending): the
    !! transfers were appended box-major, so k's are the next contiguous run at the cursor.
    impure subroutine s_amr_l1_fill_consume(q_cons_coarse, k, full)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_coarse
        integer, intent(in)                                 :: k
        logical, intent(in)                                 :: full
        integer                                             :: o1, o2, o3, lo, hi, ie, jx, boff, msl
        integer                                             :: plo(3), phi(3), crlo(3), crhi(3), bl(3), bh(3)
        integer                                             :: tb(3, 6), te(3, 6)

        call s_phase_tic(PH_GATHER)
        o1 = amr_sidx(1); o2 = amr_sidx(2); o3 = amr_sidx(3)
        call s_amr_patch_box(k, plo, phi)
        amr_cpat_off = plo
#ifdef MFC_DEBUG
        ! NaN-flood the patch before the writes land, so a consumer read of any unshipped cell NaNs within a step
        call s_amr_poison_patch_device(phi(1) - plo(1), phi(2) - plo(2), phi(3) - plo(3))
#endif
        call s_amr_rank_coarse_range(proc_rank, crlo, crhi)
        call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
        call s_amr_l1_slice_slabs(k, plo, phi, bl, bh, full, msl, tb, te)
        if (msl > 0 .and. all(bl <= bh)) call s_amr_gather_own_shell_device(q_cons_coarse, msl, tb, te, o1, o2, o3)
        do while (amr_wcur <= amr_wrecv%nx)
            if (amr_wrecv%blk(amr_wcur) /= k) exit
            if (amr_device_pack) then
                call s_amr_fx_run(k, amr_wrecv%blk, amr_wrecv%nx, amr_wcur, ie)
                boff = amr_fx_pl(7, amr_wcur) - XA_NH
                do jx = amr_wcur, ie
                    call s_amr_wave_hdr_check(amr_wrecv, amr_fw_rq, jx, XA_F1W_SND)
                end do
                call s_amr_fx_unpack(amr_wcur, ie, boff, amr_cpat_off(1), amr_cpat_off(2), amr_cpat_off(3), amr_fx_pl(:, &
                                     & 1:amr_wrecv%nx), amr_fx_pre(1:amr_wrecv%nx + 1), amr_fw_rq(boff + 1:amr_fx_pl(7, &
                                     & ie) + amr_fx_pre(ie + 1) - amr_fx_pre(ie)))
                amr_wcur = ie + 1
                cycle
            end if
            call s_amr_wave_hdr_check(amr_wrecv, amr_fw_rq, amr_wcur, XA_F1W_SND)
            call s_amr_wave_slice(amr_wrecv, amr_wcur, lo, hi)
            call s_amr_unpack_box_device(amr_wrecv%bl(:,amr_wcur), amr_wrecv%bh(:,amr_wcur), amr_fw_rq(lo:hi))
            amr_wcur = amr_wcur + 1
        end do
        call s_phase_toc(PH_GATHER)

    end subroutine s_amr_l1_fill_consume

    !> Every receive transfer of the open fill wave has been consumed.
    impure subroutine s_amr_fill_wave_done()

        @:ASSERT(amr_wcur == amr_wrecv%nx + 1, "fill wave: unconsumed recv transfers")

    end subroutine s_amr_fill_wave_done

    !> Per-stage level-1 fill: the wave, then per owned level-1 box the consume and the ghost fill.
    impure subroutine s_amr_stage_fill_wave(q_cons_coarse)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_coarse
        integer                                                :: k, kk

        call s_amr_l1_fill_exchange(q_cons_coarse, .false.)
        do kk = 1, amr_n_my
            k = amr_my_blk(kk)
            if (amr_block_level(k) /= 1) cycle
            call s_amr_select_slot(k)
            call s_amr_l1_fill_consume(q_cons_coarse, k, .false.)
            call s_phase_tic(PH_GFILL)
            call s_amr_fill_fine_ghosts(amr_cg, amr_loc_of(amr_cur))
            call s_phase_toc(PH_GFILL)
        end do
        call s_amr_fill_wave_done()

    end subroutine s_amr_stage_fill_wave

    !> Padded parent-patch extents of a child with parent footprint [plo, phi] (amr_cpat_mar coarse cells each side).
    pure subroutine s_amr_patch_width(plo, phi, w1, w2, w3)

        integer, intent(in)  :: plo(3), phi(3)
        integer, intent(out) :: w1, w2, w3

        w1 = phi(1) - plo(1) + 2*amr_cpat_mar
        w2 = merge(phi(2) - plo(2) + 2*amr_cpat_mar, 0, amr_dim(2)); w3 = merge(phi(3) - plo(3) + 2*amr_cpat_mar, 0, amr_dim(3))

    end subroutine s_amr_patch_width

    !> The parent-fill transfer list of a child with padded parent patch (w1, w2, w3), in the patch-local frame: the whole patch
    !! (full) or its hollow shell (the per-stage ghost fill never reads the open interior of the parent footprint [mar+1, w-mar-1]).
    !! Send, receive and consume all derive the list here, so the wire layout cannot drift between sides.
    impure subroutine s_amr_parent_slabs(w1, w2, w3, full, msl, tb, te)

        integer, intent(in)  :: w1, w2, w3
        logical, intent(in)  :: full
        integer, intent(out) :: msl, tb(3, 6), te(3, 6)
        integer              :: clo(3), chi(3), scells

        if (full) then
            msl = 1; tb(:,1) = 0; te(:,1) = [w1, w2, w3]
            return
        end if
        clo = merge(amr_cpat_mar + 1, 0, amr_dim); chi = merge([w1, w2, w3] - amr_cpat_mar - 1, 0, amr_dim)
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

    !> The level-lev parent-fill wave, exchange phase: every level-lev block's parent patch from its parent's fine array, in the
    !! parent-fine frame. A split child (parent owned elsewhere) is its s_amr_parent_slabs transfer list from the parent's owner to
    !! its own owner; a co-located parent is a device copy at consume, with no wire transfer. Called per level ascending, so every
    !! level-(lev-1) parent is complete (ghost-filled per stage; prolonged at init/regrid) before this wave reads its interior.
    !! Consume with s_amr_parent_fill_consume per owned level-lev box ascending.
    impure subroutine s_amr_parent_fill_exchange(lev, full)

        integer, intent(in) :: lev
        logical, intent(in) :: full
        integer             :: k, ix, pblk, powner, cowner, lo, hi, kk
        integer             :: w1, w2, w3, plo(3), phi(3), msl
        integer             :: tb(3, 6), te(3, 6)

        amr_wcur = 1
        if (amr_num_blocks <= 0) return
        call s_amr_wave_open(amr_wave, 2)

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
        call s_amr_wave_post(amr_wave, amr_wrecv, amr_fw_rq, XA_F2W_RCV, amr_fw_dev)
        ! packs: one fused launch over the whole send list under amr_device_pack, else one launch per transfer; the debug
        ! identity headers are written after the fused launch, whose copyout covers the pool prefix
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
            do ix = 1, amr_wsend%nx
                call s_amr_wave_hdr_pack(amr_wsend, amr_fw_sq, ix, XA_F2W_SND)
            end do
        else
            do ix = 1, amr_wsend%nx
                k = amr_wsend%blk(ix)
                call s_amr_parent_frame(k, plo, phi, w1, w2, w3)
                call s_amr_wave_slice(amr_wsend, ix, lo, hi)
                call s_amr_pack_parent_box_device(amr_loc_of(amr_parent_blk(k)), amr_wsend%bl(:,ix), amr_wsend%bh(:,ix), &
                                                  & amr_fw_sq(lo:hi))
                call s_amr_wave_hdr_pack(amr_wsend, amr_fw_sq, ix, XA_F2W_SND)
            end do
        end if
        call s_amr_wave_send(amr_wave, amr_wsend, amr_fw_sq, XA_F2W_SND, amr_fw_dev)
        call s_amr_wave_wait(amr_wave)
        if (amr_device_pack .and. amr_wrecv%nx > 0) call s_amr_fx_plan(amr_wrecv)
        call s_phase_toc(PH_GATHER)

    end subroutine s_amr_parent_fill_exchange

    !> The parent-fill wave, consume phase for owned box k (the current slot): set the patch frame, then a co-located parent's
    !! own-box device copies or k's received transfers (box-major: the next contiguous run at the cursor).
    impure subroutine s_amr_parent_fill_consume(k, full)

        integer, intent(in) :: k
        logical, intent(in) :: full
        integer             :: pblk, lo, hi, boff, ie, jx, isl, msl
        integer             :: w1, w2, w3, plo(3), phi(3)
        integer             :: tb(3, 6), te(3, 6)

        call s_phase_tic(PH_GATHER)
        pblk = amr_parent_blk(k)
        call s_amr_parent_frame(k, plo, phi, w1, w2, w3)
#ifdef MFC_DEBUG
        call s_amr_poison_patch_device(w1, w2, w3)
#endif
        if (amr_block_owner(pblk) == proc_rank) then
            call s_amr_parent_slabs(w1, w2, w3, full, msl, tb, te)
            do isl = 1, msl
                call s_amr_copy_parent_box(amr_loc_of(pblk), tb(:,isl), te(:,isl))
            end do
        else
            @:ASSERT(amr_wcur <= amr_wrecv%nx .and. amr_wrecv%blk(amr_wcur) == k, "parent-fill wave: missing recv transfer")
            do while (amr_wcur <= amr_wrecv%nx)
                if (amr_wrecv%blk(amr_wcur) /= k) exit
                if (amr_device_pack) then
                    call s_amr_fx_run(k, amr_wrecv%blk, amr_wrecv%nx, amr_wcur, ie)
                    boff = amr_fx_pl(7, amr_wcur) - XA_NH
                    do jx = amr_wcur, ie
                        call s_amr_wave_hdr_check(amr_wrecv, amr_fw_rq, jx, XA_F2W_SND)
                    end do
                    call s_amr_fx_unpack(amr_wcur, ie, boff, 0, 0, 0, amr_fx_pl(:,1:amr_wrecv%nx), &
                                         & amr_fx_pre(1:amr_wrecv%nx + 1), amr_fw_rq(boff + 1:amr_fx_pl(7, &
                                         & ie) + amr_fx_pre(ie + 1) - amr_fx_pre(ie)))
                    amr_wcur = ie + 1
                    cycle
                end if
                call s_amr_wave_hdr_check(amr_wrecv, amr_fw_rq, amr_wcur, XA_F2W_SND)
                call s_amr_wave_slice(amr_wrecv, amr_wcur, lo, hi)
                call s_amr_unpack_parent_box_device(amr_wrecv%bl(:,amr_wcur), amr_wrecv%bh(:,amr_wcur), amr_fw_rq(lo:hi))
                amr_wcur = amr_wcur + 1
            end do
        end if
        call s_phase_toc(PH_GATHER)

    end subroutine s_amr_parent_fill_consume

    !> Per-stage level-lev fill: the wave, then per owned level-lev box the consume and the ghost fill.
    impure subroutine s_amr_parent_fill_wave(lev)

        integer, intent(in) :: lev
        integer             :: k, kk

        call s_amr_parent_fill_exchange(lev, .false.)
        do kk = 1, amr_n_own
            k = amr_own_blk(kk)
            if (amr_block_level(k) /= lev) cycle
            call s_amr_select_slot(k)
            if (.not. amr_rank_owns_block) cycle
            call s_amr_parent_fill_consume(k, .false.)
            call s_phase_tic(PH_GFILL)
            call s_amr_fill_fine_ghosts(amr_cg, amr_loc_of(amr_cur))
            call s_phase_toc(PH_GFILL)
        end do
        call s_amr_fill_wave_done()

    end subroutine s_amr_parent_fill_wave

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
