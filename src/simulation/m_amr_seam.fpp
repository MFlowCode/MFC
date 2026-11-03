!>
!!@file
!!@brief Contains module m_amr_seam

#! AMD OpenMP: assert allocatables present on every kernel (see OMP_DEFAULT_STR); a kernel may name only module arrays allocated
#! under its own launch condition.
#:set MFC_OMP_PRESENT_ALLOCATABLE = True
#:include 'macros.fpp'

!> @brief Same-level fine-fine seams and the coarse conservative halo.
module m_amr_seam

    use m_derived_types
    use m_box, only: f_morton
    use m_global_parameters
    use m_mpi_proxy, only: s_mpi_abort  ! @:ASSERT expands to it
    use m_mpi_common, only: s_mpi_sendrecv_variables_buffers
    use m_amr_state
    use m_amr_wave
    use m_amr_distribution

    implicit none

    private
    public :: f_amr_seam, f_amr_seam_dim, s_amr_build_seam_pairs, s_amr_exchange_coarse_cons_halo, s_amr_fine_fine_post, &
        & s_amr_fine_fine_drain

contains

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

end module m_amr_seam
