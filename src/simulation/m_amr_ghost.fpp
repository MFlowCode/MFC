!>
!!@file
!!@brief Contains module m_amr_ghost

#! AMD OpenMP: assert allocatables present on every kernel (see OMP_DEFAULT_STR); a kernel may name only module arrays allocated
#! under its own launch condition.
#:set MFC_OMP_PRESENT_ALLOCATABLE = True
#:include 'macros.fpp'

!> @brief Coarse-to-fine ghost fill: slab table, limited-linear or quartic fill, quartic positivity guard, volume-fraction closure.
module m_amr_ghost

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy, only: s_mpi_abort  ! @:ASSERT expands to it
    use m_amr_state

    implicit none

    private
    public :: s_amr_slab_add, f_amr_slab_push, s_amr_queue_ghost_slabs, s_amr_fill_ghost_rows, s_amr_fill_fine_ghosts

contains

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

    !> Ghost fill over all queued slabs, one cell per thread: quartic (positivity-checked) or minmod-limited linear per variable,
    !! then the volume-fraction closure. Each row carries its block's store column, coarse origin and patch frame.
    impure subroutine s_amr_fill_ghost_rows(q_coarse)

        type(scalar_field), dimension(sys_size), intent(in) :: q_coarse
        integer                                             :: i, fi, fj, fk, ci, cj, ck, ox, oy, oz, loc
        integer                                             :: rr, lo1, lo2, lo3, advb, adve
        integer                                             :: s, nsl, lo, hi, md, g, r, n1, n2, stot
        integer                                             :: ky, kz, wy, wz
        logical                                             :: d2, d3, multi, quart, lin, shx, shy, shz
        real(wp)                                            :: u0, xm, xp, ym, yp, zm, zp, xix, xiy, xiz
        real(wp)                                            :: av, asum
        real(wp)                                            :: q5(-2:2), pl(-2:2), ln(-2:2)  !< staged x, y, z

        ! patch-local coarse index = row coarse origin (global) + f/rr - patch frame; fine indices are block-local

        d2 = n_glb > 0; d3 = p_glb > 0
        rr = amr_ref_ratio
        quart = amr_prolong == 2; wy = merge(2, 0, d2); wz = merge(2, 0, d3)
        multi = num_fluids > 1 .and. (.not. bubbles_lagrange)  ! EL alphas sum to beta, not 1: no sum-to-one closure
        advb = eqn_idx%adv%beg; adve = eqn_idx%adv%end
        stot = f_amr_slab_push(nsl)
        if (stot == 0) return  ! no level-1 block on this rank: the store may be unallocated
        $:GPU_PARALLEL_LOOP(private='[s, lo, hi, md, r, n1, n2, fi, fj, fk, i, ci, cj, ck, xix, xiy, xiz, u0, xm, xp, ym, yp, zm, &
                            & zp, ky, kz, q5, pl, ln, ox, oy, oz, lo1, lo2, lo3, loc, lin, av, asum, shx, shy, shz]')
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
            lin = .not. quart
            if (quart) then
                $:GPU_LOOP(parallelism='[seq]')
                do i = 1, sys_size
                    if (multi .and. i >= advb .and. i <= adve) cycle  ! volume fractions: closure below
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
                end do
                lin = .not. f_amr_cons_ok(fi, fj, fk, loc)  ! inadmissible quartic cell: limited linear instead
            end if
            if (lin) then
                $:GPU_LOOP(parallelism='[seq]')
                do i = 1, sys_size
                    if (multi .and. i >= advb .and. i <= adve) cycle
                    u0 = real(q_coarse(i)%sf(ci, cj, ck), wp)
                    xm = real(q_coarse(i)%sf(ci - 1, cj, ck), wp); xp = real(q_coarse(i)%sf(ci + 1, cj, ck), wp)
                    ym = u0; yp = u0
                    if (d2) then
                        ym = real(q_coarse(i)%sf(ci, cj - 1, ck), wp); yp = real(q_coarse(i)%sf(ci, cj + 1, ck), wp)
                    end if
                    zm = u0; zp = u0
                    if (d3) then
                        zm = real(q_coarse(i)%sf(ci, cj, ck - 1), wp); zp = real(q_coarse(i)%sf(ci, cj, ck + 1), wp)
                    end if
                    amr_cons_st(fi, fj, fk, i, loc) = f_amr_lin(u0, xm, xp, ym, yp, zm, zp, xix, xiy, xiz)
                end do
            end if
            if (multi) then
                ! volume fractions: shared limiter switch over all fluids, interpolate + clamp advb..adve-1, alpha_n = 1 - sum
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
                    xm = u0; xp = u0
                    if (shx) then
                        xm = real(q_coarse(i)%sf(ci - 1, cj, ck), wp); xp = real(q_coarse(i)%sf(ci + 1, cj, ck), wp)
                    end if
                    ym = u0; yp = u0
                    if (shy) then
                        ym = real(q_coarse(i)%sf(ci, cj - 1, ck), wp); yp = real(q_coarse(i)%sf(ci, cj + 1, ck), wp)
                    end if
                    zm = u0; zp = u0
                    if (shz) then
                        zm = real(q_coarse(i)%sf(ci, cj, ck - 1), wp); zp = real(q_coarse(i)%sf(ci, cj, ck + 1), wp)
                    end if
                    av = min(max(f_amr_lin(u0, xm, xp, ym, yp, zm, zp, xix, xiy, xiz), 0._wp), 1._wp)
                    amr_cons_st(fi, fj, fk, i, loc) = av
                    asum = asum + av
                end do
                amr_cons_st(fi, fj, fk, adve, loc) = 1._wp - asum
            end if
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_fill_ghost_rows

end module m_amr_ghost
