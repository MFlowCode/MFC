!>
!!@file
!!@brief Contains module m_amr_exchange

#! AMD OpenMP: assert allocatables present on every kernel (see OMP_DEFAULT_STR); a kernel may name only module arrays allocated
#! under its own launch condition.
#:set MFC_OMP_PRESENT_ALLOCATABLE = True
#:include 'macros.fpp'

!> @brief Parent/child and fine-fine data exchange: gather plans, pack/unpack, seams, ghost fills and the fill waves.
module m_amr_exchange

    use m_derived_types
    use m_box, only: f_morton
    use m_global_parameters
    use m_mpi_proxy, only: s_mpi_abort  ! @:ASSERT expands to it
    use m_mpi_common, only: s_mpi_sendrecv_variables_buffers
    use m_amr_registers, only: s_amr_parent_foot
    use m_phase_timing
    use m_amr_state
    use m_amr_wave
    use m_amr_distribution

    implicit none

    private
    public :: f_amr_seam, f_amr_seam_dim, s_amr_build_seam_pairs, s_amr_exchange_coarse_cons_halo, s_amr_fine_fine_drain, &
        & s_amr_fine_fine_post, s_amr_fill_wave_done, s_amr_l1_fill_exchange, s_amr_l1_fill_consume, s_amr_parent_fill_exchange, &
        & s_amr_parent_fill_consume, s_amr_parent_fill_wave, s_amr_stage_fill_wave, s_amr_pack_unpack_block_st, s_amr_slab_add, &
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

    !> Split the fine block's ghost shell into ns disjoint face slabs: x spans the full transverse extent, y is restricted to the
    !! interior in x, z in x and y.
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

    !> Append slab [lo, hi] to amr_slab_tab (rows: 1:3 lo, 4:6 hi, 7 flat offset, 8 cells, 9 store column, 10:12 coarse origin,
    !! 13:15 patch frame).
    impure subroutine s_amr_slab_add(lo, hi, loc, ilo, coff)

        integer, intent(in)  :: lo(3), hi(3), loc, ilo(3), coff(3)
        integer, allocatable :: tmp(:,:)

        if (amr_slab_n == size(amr_slab_tab, 2)) then
            tmp = amr_slab_tab(:,1:amr_slab_n)
            @:DEALLOCATE(amr_slab_tab)
            @:ALLOCATE(amr_slab_tab(1:15, 1:2*amr_slab_n))
            amr_slab_tab(:,1:amr_slab_n) = tmp
        end if
        amr_slab_n = amr_slab_n + 1
        associate (row => amr_slab_tab(:,amr_slab_n))
            row(1:3) = lo; row(4:6) = hi; row(9) = loc; row(10:12) = ilo; row(13:15) = coff
            row(7) = 0
            if (amr_slab_n > 1) row(7) = amr_slab_tab(7, amr_slab_n - 1) + amr_slab_tab(8, amr_slab_n - 1)
            row(8) = product(hi - lo + 1)
        end associate

    end subroutine s_amr_slab_add

    !> Push the slab table to the device and empty it for the next build; returns the total cell count.
    impure integer function f_amr_slab_push(nsl) result(stot)

        integer, intent(out) :: nsl

        nsl = amr_slab_n
        stot = 0
        if (nsl > 0) stot = amr_slab_tab(7, nsl) + amr_slab_tab(8, nsl)
        if (nsl > 0) then
            $:GPU_UPDATE(device='[amr_slab_tab(:, 1:nsl)]')
        end if
        amr_slab_n = 0

    end function f_amr_slab_push

    !> Flat-index decode: binary-search row 7 of the slab table for the slab holding flat cell g, then its cell (a1, a2, a3).

    !> Queue the current block's ghost slabs (store column loc, patch at pool offset pz) for s_amr_fill_ghost_rows.
    impure subroutine s_amr_queue_ghost_slabs(loc, pz)

        integer, intent(in) :: loc, pz
        integer             :: ns, sb(3, 6), se(3, 6), j

        call s_amr_build_ghost_slabs(ns, sb, se)
        do j = 1, ns
            call s_amr_slab_add(sb(:,j), se(:,j), loc, amr_isect_lo, amr_cpat_off - [0, 0, pz])
        end do

    end subroutine s_amr_queue_ghost_slabs

    !> Fill the fine ghost shell by conservative-linear prolongation from the gathered coarse patch q_coarse (amr_cg); the interior
    !! is untouched. Writes the store at local index loc.
    impure subroutine s_amr_fill_fine_ghosts(q_coarse, loc)

        type(scalar_field), dimension(sys_size), intent(in) :: q_coarse
        integer, intent(in)                                 :: loc

        call s_amr_queue_ghost_slabs(loc, 0)
        call s_amr_fill_ghost_rows(q_coarse)

    end subroutine s_amr_fill_fine_ghosts

    !> Ghost fill over all queued slabs in one kernel; each row carries its block's store column, coarse origin and patch frame.
    impure subroutine s_amr_fill_ghost_rows(q_coarse)

        type(scalar_field), dimension(sys_size), intent(in) :: q_coarse
        integer                                             :: i, fi, fj, fk, ci, cj, ck, ox, oy, oz, loc
        integer                                             :: rr, lo1, lo2, lo3
        integer                                             :: advb, adve
        integer                                             :: s, nsl, lo, hi, md
        integer                                             :: g, r, n1, n2, stot
        logical                                             :: d2, d3, multi, quart
        real(wp)                                            :: u0, sx, sy, sz, xix, xiy, xiz
        integer                                             :: ky, kz, wy, wz
        real(wp)                                            :: q5(-2:2), pl(-2:2), ln(-2:2)  !< staged x, y, z

        ! patch-local coarse index = row coarse origin (global) + f/rr - patch frame; fine indices are block-local

        d2 = n_glb > 0; d3 = p_glb > 0
        rr = amr_ref_ratio
        quart = amr_prolong == 2; wy = merge(2, 0, d2); wz = merge(2, 0, d3)
        multi = num_fluids > 1 .and. (.not. bubbles_lagrange)  ! EL alphas sum to beta, not 1: no sum-to-one closure
        advb = eqn_idx%adv%beg; adve = eqn_idx%adv%end
        stot = f_amr_slab_push(nsl)
        if (stot == 0) return  ! no level-1 block on this rank: the store may be unallocated
        $:GPU_PARALLEL_LOOP(collapse=2, private='[s, lo, hi, md, r, n1, n2, fi, fj, fk, ci, cj, ck, xix, xiy, xiz, u0, sx, sy, &
                            & sz, ky, kz, q5, pl, ln, ox, oy, oz, lo1, lo2, lo3, loc]')
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
                fi = amr_slab_tab(1, s) + mod(r, n1)
                fj = amr_slab_tab(2, s) + mod(r/n1, n2)
                fk = amr_slab_tab(3, s) + r/(n1*n2)
                loc = amr_slab_tab(9, s)
                lo1 = amr_slab_tab(10, s); lo2 = amr_slab_tab(11, s); lo3 = amr_slab_tab(12, s)
                ox = amr_slab_tab(13, s); oy = amr_slab_tab(14, s); oz = amr_slab_tab(15, s)
                ! multi-fluid: volume fractions are handled by the closure below
                if (.not. (multi .and. i >= advb .and. i <= adve)) then
                    ck = -oz; xiz = 0._wp  ! collapsed z: oz holds minus the patch slot
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
                    if (quart) then
                        do kz = -wz, wz
                            do ky = -wy, wy
                                q5 = real(q_coarse(i)%sf(ci - 2:ci + 2,cj + ky, ck + kz), wp)
                                pl(ky) = f_amr_quartic_half(q5, xix)
                            end do
                            ln(kz) = pl(0)
                            if (d2) ln(kz) = f_amr_quartic_half(pl, xiy)
                        end do
                        u0 = ln(0)
                        if (d3) u0 = f_amr_quartic_half(ln, xiz)
                        amr_cons_st(fi, fj, fk, i, loc) = u0
                    else
                        u0 = real(q_coarse(i)%sf(ci, cj, ck), wp)
                        sx = minmod(real(q_coarse(i)%sf(ci + 1, cj, ck), wp) - u0, u0 - real(q_coarse(i)%sf(ci - 1, cj, ck), wp))
                        sy = 0._wp
                        if (d2) sy = minmod(real(q_coarse(i)%sf(ci, cj + 1, ck), wp) - u0, u0 - real(q_coarse(i)%sf(ci, cj - 1, &
                            & ck), wp))
                        sz = 0._wp
                        if (d3) sz = minmod(real(q_coarse(i)%sf(ci, cj, ck + 1), wp) - u0, u0 - real(q_coarse(i)%sf(ci, cj, &
                            & ck - 1), wp))
                        amr_cons_st(fi, fj, fk, i, loc) = u0 + sx*xix + sy*xiy + sz*xiz
                    end if
                end if
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! own routine: nvfortran 25.x crashes on two target regions sharing this routine's privates
        if (multi) call s_amr_fill_fine_ghosts_alphas(q_coarse, nsl, stot, advb, adve)
        if (quart) call s_amr_guard_fine_ghosts(q_coarse, nsl, stot, multi, advb, adve)

    end subroutine s_amr_fill_ghost_rows

    !> Quartic safety net: refill inadmissible (non-positive density or internal energy) ghost cells with minmod-limited linear;
    !! multi-fluid volume fractions keep their closure values.
    impure subroutine s_amr_guard_fine_ghosts(q_coarse, nsl, stot, multi, advb, adve)

        type(scalar_field), dimension(sys_size), intent(in) :: q_coarse
        integer, intent(in)                                 :: nsl, stot, advb, adve
        logical, intent(in)                                 :: multi
        integer                                             :: i, fi, fj, fk, ci, cj, ck, ox, oy, oz, loc
        integer                                             :: rr, lo1, lo2, lo3, s, g, r, n1, n2, lo, hi, md
        integer                                             :: amr_v
        logical                                             :: d2, d3, ok
        real(wp)                                            :: u0, sx, sy, sz, xix, xiy, xiz, amr_rho, amr_ke

        d2 = n_glb > 0; d3 = p_glb > 0
        rr = amr_ref_ratio
        $:GPU_PARALLEL_LOOP(private='[s, lo, hi, md, r, n1, n2, fi, fj, fk, i, ci, cj, ck, xix, xiy, xiz, u0, sx, sy, sz, ox, oy, &
                            & oz, lo1, lo2, lo3, loc, ok, amr_v, amr_rho, amr_ke]')
        do g = 0, stot - 1
            lo = 1; hi = nsl
            do while (lo < hi)
                md = (lo + hi + 1)/2
                if (amr_slab_tab(7, md) <= g) then; lo = md; else; hi = md - 1; end if
            end do
            s = lo
            r = g - amr_slab_tab(7, s)
            n1 = amr_slab_tab(4, s) - amr_slab_tab(1, s) + 1; n2 = amr_slab_tab(5, s) - amr_slab_tab(2, s) + 1
            fi = amr_slab_tab(1, s) + mod(r, n1)
            fj = amr_slab_tab(2, s) + mod(r/n1, n2)
            fk = amr_slab_tab(3, s) + r/(n1*n2)
            loc = amr_slab_tab(9, s)
            amr_rho = 0._wp; amr_ke = 0._wp
            $:GPU_LOOP(parallelism='[seq]')
            do amr_v = eqn_idx%cont%beg, eqn_idx%cont%end
                amr_rho = amr_rho + real(amr_cons_st(fi, fj, fk, amr_v, loc), wp)
            end do
            ok = amr_rho > 0._wp
            if (ok .and. model_eqns /= 4) then
                $:GPU_LOOP(parallelism='[seq]')
                do amr_v = eqn_idx%mom%beg, eqn_idx%mom%end
                    amr_ke = amr_ke + real(amr_cons_st(fi, fj, fk, amr_v, loc), wp)**2
                end do
                ok = real(amr_cons_st(fi, fj, fk, eqn_idx%E, loc), wp) - 0.5_wp*amr_ke/amr_rho > 0._wp
            end if
            if (ok) cycle
            lo1 = amr_slab_tab(10, s); lo2 = amr_slab_tab(11, s); lo3 = amr_slab_tab(12, s)
            ox = amr_slab_tab(13, s); oy = amr_slab_tab(14, s); oz = amr_slab_tab(15, s)
            ck = -oz; xiz = 0._wp
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
            $:GPU_LOOP(parallelism='[seq]')
            do i = 1, sys_size
                if (multi .and. i >= advb .and. i <= adve) cycle
                u0 = real(q_coarse(i)%sf(ci, cj, ck), wp)
                sx = minmod(real(q_coarse(i)%sf(ci + 1, cj, ck), wp) - u0, u0 - real(q_coarse(i)%sf(ci - 1, cj, ck), wp))
                sy = 0._wp
                if (d2) sy = minmod(real(q_coarse(i)%sf(ci, cj + 1, ck), wp) - u0, u0 - real(q_coarse(i)%sf(ci, cj - 1, ck), wp))
                sz = 0._wp
                if (d3) sz = minmod(real(q_coarse(i)%sf(ci, cj, ck + 1), wp) - u0, u0 - real(q_coarse(i)%sf(ci, cj, ck - 1), wp))
                amr_cons_st(fi, fj, fk, i, loc) = u0 + sx*xix + sy*xiy + sz*xiz
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_guard_fine_ghosts

    !> Volume-fraction ghosts: shared limiter switch over all fluids, interpolate + clamp advb..adve-1, alpha_n = 1 - sum.
    impure subroutine s_amr_fill_fine_ghosts_alphas(q_coarse, nsl, stot, advb, adve)

        type(scalar_field), dimension(sys_size), intent(in) :: q_coarse
        integer, intent(in)                                 :: nsl, stot, advb, adve
        integer                                             :: i, fi, fj, fk, ci, cj, ck, ox, oy, oz, loc
        integer                                             :: rr, lo1, lo2, lo3, s, g, r, n1, n2, lo, hi, md
        logical                                             :: d2, d3, shx, shy, shz
        real(wp)                                            :: u0, sx, sy, sz, xix, xiy, xiz, av, asum

        d2 = n_glb > 0; d3 = p_glb > 0
        rr = amr_ref_ratio
        $:GPU_PARALLEL_LOOP(private='[s, lo, hi, md, r, n1, n2, fi, fj, fk, i, ci, cj, ck, xix, xiy, xiz, u0, sx, sy, sz, av, &
                            & asum, shx, shy, shz, ox, oy, oz, lo1, lo2, lo3, loc]')
        do g = 0, stot - 1
            lo = 1; hi = nsl
            do while (lo < hi)
                md = (lo + hi + 1)/2
                if (amr_slab_tab(7, md) <= g) then; lo = md; else; hi = md - 1; end if
            end do
            s = lo
            r = g - amr_slab_tab(7, s)
            n1 = amr_slab_tab(4, s) - amr_slab_tab(1, s) + 1; n2 = amr_slab_tab(5, s) - amr_slab_tab(2, s) + 1
            fi = amr_slab_tab(1, s) + mod(r, n1)
            fj = amr_slab_tab(2, s) + mod(r/n1, n2)
            fk = amr_slab_tab(3, s) + r/(n1*n2)
            loc = amr_slab_tab(9, s)
            lo1 = amr_slab_tab(10, s); lo2 = amr_slab_tab(11, s); lo3 = amr_slab_tab(12, s)
            ox = amr_slab_tab(13, s); oy = amr_slab_tab(14, s); oz = amr_slab_tab(15, s)
            ck = -oz; xiz = 0._wp
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

    !> Exchange coarse conservative ghost layers at internal rank boundaries (beg then end per direction); the solver fills only
    !! prim ghosts. Collective.
    impure subroutine s_amr_exchange_coarse_cons_halo(q_cons)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons

        #:for D, X in [(1, 'x'), (2, 'y'), (3, 'z')]
            if (amr_dim(${D}$)) then
                if (bc_${X}$%beg >= 0) call s_mpi_sendrecv_variables_buffers(q_cons, ${D}$, -1, sys_size)
                if (bc_${X}$%end >= 0) call s_mpi_sendrecv_variables_buffers(q_cons, ${D}$, 1, sys_size)
            end if
        #:endfor

    end subroutine s_amr_exchange_coarse_cons_halo

    !> True iff yb sits directly above xb on xb's high face in dim d with identical transverse coarse extents (one ordered pair per
    !! seam).
    pure logical function f_amr_seam(xb, yb, d) result(seam)

        integer, intent(in) :: xb, yb, d
        integer             :: t

        seam = amr_region_lo_all(d, yb) == amr_region_hi_all(d, xb) + 1
        do t = 1, 3
            if (t /= d) seam = seam .and. amr_region_lo_all(t, xb) == amr_region_lo_all(t, yb) .and. amr_region_hi_all(t, &
                & xb) == amr_region_hi_all(t, yb)
        end do

    end function f_amr_seam

    !> Same-rank seam exchange, both directions in one kernel: xb high interior -> yb low ghost and yb low interior -> xb high
    !! ghost. The four slabs are disjoint, so fusing is safe. Threads span the max pair extent and are masked.
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
                        ! (a, b) transverse, t depth into the seam, placed per the pair's seam dim
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

    !> Seam dim of the ordered pair (xb, yb): the d in which yb is xb's same-level high-face neighbour, else 0 (adjacent in at most
    !! one dim).
    pure integer function f_amr_seam_dim(xb, yb) result(d)

        integer, intent(in) :: xb, yb

        d = 0
        if (xb == yb) return
        if (amr_block_level(xb) /= amr_block_level(yb)) return
        if (f_amr_seam(xb, yb, 1)) d = 1
        if (n_glb > 0) then; if (f_amr_seam(xb, yb, 2)) d = 2; end if
        if (p_glb > 0) then; if (f_amr_seam(xb, yb, 3)) d = 3; end if

    end function f_amr_seam_dim

    !> Binary-search the Morton-sorted lo corners for the block at clo, verify the seam predicate against xb and record it in (mb,
    !! md, nm); a block already recorded is raised to the higher d.
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
        ! scan the equal-key run: f_morton keeps 21 bits/dim, so the lo check rejects aliased corners
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

    !> Rebuild the same-level seam-pair list once per regrid/restart in identical (xb, yb) order on all ranks, and the per-block
    !! overlap-rank lists (amr_ovl_gather/scatter).
    impure subroutine s_amr_build_seam_pairs()

        integer                      :: xb, d, np, k, mx, pass, nm, im, jm, tb, td, nb
        integer                      :: plo(3), phi(3), rlo(3), rhi(3)
        integer                      :: mb(3), md(3), clo(3)
        integer(kind=8), allocatable :: mkey(:)
        integer, allocatable         :: ord(:)

        if (allocated(amr_seam_pairs)) deallocate (amr_seam_pairs)

        ! Blocks are disjoint and the seam fixes the neighbour's lo corner, so each (block, dim) needs one Morton binary
        ! search. Emission order (xb, then yb ascending) must match on all ranks or sends and recvs mismatch and deadlock.
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
        ! per-block overlap-rank lists, rank-ascending (same MPI order as a 0..num_procs-1 scan); the first dim is the max
        ! overlap, reallocated each build. Consumers run behind a build_seam_pairs guard.
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

    !> Fine-fine halo between tiled sub-blocks: overwrite seam ghosts with the neighbour's stage-entry interior (cross-rank via the
    !! seam wave, same-rank by local copy). Wire is wp, cast to stp on unpack.
    impure subroutine s_amr_fine_fine_post()

        integer :: xb, yb, d, rX, rY, cnt, xm(3), tsz, fmul, idx, lo, hi
        integer :: r, sblk, sdlo, sdhi, ublk, udlo, udhi, eblk, edlo, edhi

        amr_sw_nsame = 0
        if (.not. amr) return
        if (amr_num_blocks < 2) return

        if (amr_seam_pairs_dirty .or. amr_seam_pairs_nblk /= amr_num_blocks) call s_amr_build_seam_pairs()
        if (allocated(amr_sw_plx)) deallocate (amr_sw_plx, amr_sw_ply, amr_sw_pd, amr_sw_pxhi, amr_sw_pfm)
        allocate (amr_sw_plx(amr_num_seam_pairs), amr_sw_ply(amr_num_seam_pairs), amr_sw_pd(amr_num_seam_pairs), &
                  & amr_sw_pxhi(amr_num_seam_pairs), amr_sw_pfm(3, amr_num_seam_pairs))
        ! Seam wave: one aggregated message per peer and direction; both owners walk the same replicated pair list, so the
        ! layout agrees. The wave stays open across the fill waves (own handle and sides).
        call s_amr_wave_open(amr_wave_seam, amr_band_seam)
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
            if (rX == rY) then  ! same rank owns both: deferred to the batched kernel at drain
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
        call s_amr_wave_post(amr_wave_seam, amr_wseam_rcv, amr_sw_rq, .false.)
        do idx = 1, amr_wseam_snd%nx
            call s_amr_wave_slice(amr_wseam_snd, idx, lo, hi)
            call s_amr_fine_slice(amr_wseam_snd%blk(idx), amr_wseam_snd%bl(1, idx), amr_wseam_snd%bl(2, idx), amr_wseam_snd%bl(3, &
                                  & idx), amr_sw_sq(lo:hi), 1)
        end do
        call s_amr_wave_send(amr_wave_seam, amr_wseam_snd, amr_sw_sq, .false.)

    end subroutine s_amr_fine_fine_post

    impure subroutine s_amr_fine_fine_drain()

        integer :: idx, lo, hi

        if (f_amr_wave_nreq(amr_wave_seam) == 0 .and. amr_sw_nsame == 0) return
        call s_amr_wave_wait(amr_wave_seam)
        do idx = 1, amr_wseam_rcv%nx
            call s_amr_wave_slice(amr_wseam_rcv, idx, lo, hi)
            call s_amr_fine_slice(amr_wseam_rcv%blk(idx), amr_wseam_rcv%bl(1, idx), amr_wseam_rcv%bl(2, idx), amr_wseam_rcv%bl(3, &
                                  & idx), amr_sw_rq(lo:hi), -1)
        end do
        if (amr_sw_nsame > 0) call s_amr_fine_seam_exchange(amr_sw_nsame, amr_sw_plx(1:amr_sw_nsame), amr_sw_ply(1:amr_sw_nsame), &
            & amr_sw_pd(1:amr_sw_nsame), amr_sw_pxhi(1:amr_sw_nsame), amr_sw_pfm(:,1:amr_sw_nsame), buff_size)
        call s_amr_select_slot(1)

    end subroutine s_amr_fine_fine_drain

    !> Flat-index decode shared by the fused kernels: binary-search the plan prefix for the transfer owning flat element e, then
    !! invert its buffer index into (i, g3, g2, g1).

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
