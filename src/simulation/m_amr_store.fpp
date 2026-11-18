!>
!!@file
!!@brief Contains module m_amr_store

#! AMD OpenMP lane: assert allocatables present on every kernel here (see OMP_DEFAULT_STR). A kernel naming an unallocated
#! module array aborts, so conditionally allocated ones (sw_jac/jac: igr) are named only under their allocation condition.
#:set MFC_OMP_PRESENT_ALLOCATABLE = True
#:include 'macros.fpp'

#:set GS_ARRAYS = [(D, A) for D, X in [(1, 'x'), (2, 'y'), (3, 'z')] for A in [X + '_cb', X + '_cc', 'd' + X]]

!> @brief Flat per-block field store: dense slot indexing, reservation, alloc/free and the prim/bridge loaders.
module m_amr_store

    use m_derived_types  ! scalar_field, t_box, int_bounds_info
    use m_global_parameters
    use m_amr_state
    use m_amr_distribution

    implicit none

    !> Staging for the coordinate push of s_amr_sync_grid_state_to_device: the active dims' cb/cc/d arrays back to back.
    real(wp), allocatable :: amr_gs_buf(:)
    $:GPU_DECLARE(create='[amr_gs_buf]')

    private
    public :: s_amr_alloc_slot, s_amr_bat_member_prim, s_amr_br_load, s_amr_br_load_batch, s_amr_br_load_faces, s_amr_br_store, &
        & s_amr_br_store_faces, s_amr_copy_fine_fields, s_amr_free_slot, s_amr_alloc_pool, s_amr_free_pool, s_amr_set_mbuf, &
        & s_amr_init_swap_buffers, s_amr_free_swap_buffers, s_amr_reconcile_slots, s_amr_sync_grid_state_to_device, s_amr_dev_mem

contains

    !> Push the host global grid state (m/n/p, idwint/idwbuff, coordinate arrays) to its GPU_DECLARE'd device copies after a
    !! swap/restore. Coordinates go through amr_gs_buf as one update since each allocatable update costs several synchronous copies.
    !! No-op on CPU.
    impure subroutine s_amr_sync_grid_state_to_device()

        integer :: o, j, lb, ub

        $:GPU_UPDATE(device='[m, n, p, idwint, idwbuff]')
        if (.not. allocated(amr_gs_buf)) then
            o = 0
            #:for D, A in GS_ARRAYS
                if (amr_dim(${D}$)) o = o + size(${A}$)
            #:endfor
            @:ALLOCATE(amr_gs_buf(1:o))
        end if
        o = 0
        #:for D, A in GS_ARRAYS
            if (amr_dim(${D}$)) then
                amr_gs_buf(o + 1:o + size(${A}$)) = ${A}$; o = o + size(${A}$)
            end if
        #:endfor
        $:GPU_UPDATE(device='[amr_gs_buf(1:o)]')
        o = 0
        #:for D, A in GS_ARRAYS
            if (amr_dim(${D}$)) then
                lb = lbound(${A}$, 1); ub = ubound(${A}$, 1)
                $:GPU_PARALLEL_LOOP()
                do j = lb, ub
                    ${A}$(j) = amr_gs_buf(o + 1 + j - lb)
                end do
                $:END_GPU_PARALLEL_LOOP()
                o = o + ub - lb + 1
            end if
        #:endfor

    end subroutine s_amr_sync_grid_state_to_device

    !> Device copy amr_cons_st -> amr_stor_st over block k's buffered box, all sys_size fields (RK step-entry backup).
    impure subroutine s_amr_copy_fine_fields(k)

        integer, intent(in) :: k  !< block slot (its flat-store column)
        integer             :: i, fi, fj, fk, loc, b1, e1, b2, e2, b3, e3

        loc = amr_loc_of(k)
        b1 = amr_slots(k)%idwbuff(1)%beg; e1 = amr_slots(k)%idwbuff(1)%end
        b2 = amr_slots(k)%idwbuff(2)%beg; e2 = amr_slots(k)%idwbuff(2)%end
        b3 = amr_slots(k)%idwbuff(3)%beg; e3 = amr_slots(k)%idwbuff(3)%end
        $:GPU_PARALLEL_LOOP(collapse=4)
        do i = 1, sys_size
            do fk = b3, e3
                do fj = b2, e2
                    do fi = b1, e1
                        amr_stor_st(fi, fj, fk, i, loc) = amr_cons_st(fi, fj, fk, i, loc)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_copy_fine_fields

    !> Block-metadata pool for amr_max_blocks slots (geometry, ownership, level, SFC cuts, overlap-list counts; the 2D rank lists
    !! are sized in s_amr_build_seam_pairs) and the dense-index maps.
    impure subroutine s_amr_alloc_pool()

        allocate (amr_slots(1:amr_max_blocks))
        call s_amr_loc_index_init()
        allocate (amr_region_lo_all(3, amr_max_blocks), amr_region_hi_all(3, amr_max_blocks))
        allocate (amr_isect_lo_all(3, amr_max_blocks), amr_isect_hi_all(3, amr_max_blocks))
        allocate (amr_owns_all(amr_max_blocks), amr_block_owner(amr_max_blocks), amr_block_level(amr_max_blocks))
        allocate (amr_fine_cut(0:num_procs - 1,1:max(amr_max_level, 1))); amr_fine_cut = -1_8
        allocate (amr_ovl_gather_n(amr_max_blocks), amr_ovl_scatter_n(amr_max_blocks))
        allocate (amr_slot_live(amr_max_blocks)); amr_slot_live = .false.
        amr_region_lo_all = 0; amr_region_hi_all = 0; amr_isect_lo_all = 0; amr_isect_hi_all = 0; amr_owns_all = .false.
        amr_block_owner = 0

    end subroutine s_amr_alloc_pool

    !> Bounce buffers for the copy-based coordinate swap; shaped like the base-level *_alloc arrays, since they are assigned whole
    !! to/from x_cb etc.
    impure subroutine s_amr_init_swap_buffers()

        #:for D, X, E in [(1, 'x', 'm'), (2, 'y', 'n'), (3, 'z', 'p')]
            if (amr_dim(${D}$)) then
                allocate (sw_${X}$_cb(-1 - buff_size:${E}$_alloc + buff_size), sw_${X}$_cc(-buff_size:${E}$_alloc + buff_size), &
                          & sw_d${X}$(-buff_size:${E}$_alloc + buff_size))
            end if
        #:endfor
        if (igr) then
            @:ALLOCATE(sw_jac(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end))
            @:ALLOCATE(sw_jac_old(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end))
        end if

    end subroutine s_amr_init_swap_buffers

    !> Free what s_amr_alloc_pool allocated (plus the per-block caches sized later), and the coordinate swap buffers.
    impure subroutine s_amr_free_pool()

        deallocate (amr_slot_live)
        call s_amr_st_finalize()
        if (allocated(amr_seam_pairs)) deallocate (amr_seam_pairs)
        if (allocated(amr_ovl_gather)) deallocate (amr_ovl_gather)
        if (allocated(amr_ovl_scatter)) deallocate (amr_ovl_scatter)
        deallocate (amr_ovl_gather_n, amr_ovl_scatter_n, amr_slots, amr_block_owner, amr_block_level, amr_fine_cut)
        deallocate (amr_region_lo_all, amr_region_hi_all, amr_isect_lo_all, amr_isect_hi_all, amr_owns_all)

    end subroutine s_amr_free_pool

    !> The buffered per-slot array bounds from the max fine extents (collapsed dims 0:0).
    impure subroutine s_amr_set_mbuf()

        mbuf_lo = merge(-buff_size, 0, amr_dim); mbuf_hi = merge(max_f + buff_size, 0, amr_dim)

    end subroutine s_amr_set_mbuf

    impure subroutine s_amr_free_swap_buffers()

        #:for D, X, E in [(1, 'x', 'm'), (2, 'y', 'n'), (3, 'z', 'p')]
            if (allocated(sw_${X}$_cb)) deallocate (sw_${X}$_cb, sw_${X}$_cc, sw_d${X}$)
            if (allocated(amr_g${X}$cb)) deallocate (amr_g${X}$cb)
        #:endfor
        if (igr) then
            @:DEALLOCATE(sw_jac)
            @:DEALLOCATE(sw_jac_old)
        end if

    end subroutine s_amr_free_swap_buffers

    !> Allocate/reset the dense local-index maps (idempotent).
    impure subroutine s_amr_loc_index_init()

        if (.not. allocated(amr_loc_of)) allocate (amr_loc_of(1:amr_max_blocks))
        if (.not. allocated(amr_loc_free)) allocate (amr_loc_free(1:amr_max_blocks))
        amr_loc_of = 0; amr_loc_free = 0; amr_loc_n = 0; amr_loc_nfree = 0

    end subroutine s_amr_loc_index_init

    !> Move one local slot's store data src -> dst on the device, in place (a staging copy would transiently double the store).
    !! Ascending-source order (s_amr_compact_store) guarantees dst's old contents are consumed. The host copy goes stale, the
    !! store's normal state between rebuilds.
    impure subroutine s_amr_st_move_slot(src, dst)

        integer, intent(in) :: src, dst
        integer             :: i, j, k, l
        integer             :: blo1, bhi1, blo2, bhi2, blo3, bhi3

        ! scalar bounds: mbuf_lo/hi are host-only, and CCE OpenACC aborts on them in a default(present) kernel

        blo1 = mbuf_lo(1); bhi1 = mbuf_hi(1)
        blo2 = mbuf_lo(2); bhi2 = mbuf_hi(2)
        blo3 = mbuf_lo(3); bhi3 = mbuf_hi(3)

        #:for ST in ['amr_cons_st', 'amr_stor_st']
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do l = blo3, bhi3
                    do k = blo2, bhi2
                        do j = blo1, bhi1
                            ${ST}$(j, k, l, i, dst) = ${ST}$(j, k, l, i, src)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        #:endfor

    end subroutine s_amr_st_move_slot

    !> Re-densify the local index space after every reconcile so amr_loc_n tracks the live count instead of ratcheting up with SFC
    !! window shifts and migration slots (capacity then plateaus at the rebuild transient). Moves are device-side and in place in
    !! ascending old-index order: each destination is <= its source and below every pending source. Called at the end of
    !! s_amr_reconcile_slots, after the rebuild's overlap carry-forward (which reads amr_stor_st at old indices) has finished.
    impure subroutine s_amr_compact_store()

        integer              :: k, v, newloc
        integer, allocatable :: inv(:)

        if (amr_st_cap <= 0) return
        ! invert the map: inv(local index) = global slot, 0 if free
        allocate (inv(amr_st_cap)); inv = 0
        do k = 1, amr_max_blocks
            if (amr_loc_of(k) > 0) inv(amr_loc_of(k)) = k
        end do
        newloc = 0
        do v = 1, amr_st_cap
            if (inv(v) == 0) cycle
            newloc = newloc + 1
            amr_loc_of(inv(v)) = newloc
            if (newloc /= v) call s_amr_st_move_slot(v, newloc)
        end do
        deallocate (inv)
        amr_loc_n = newloc
        amr_loc_nfree = 0  ! every recycled index is invalid after renumbering

    end subroutine s_amr_compact_store

    !> Free and total device memory in bytes (hipMemGetInfo on the AMD GPU build); -1 for both where there is no query.
    impure subroutine s_amr_dev_mem(free_b, total_b)

        use iso_c_binding, only: c_int, c_size_t
        integer(8), intent(out) :: free_b, total_b

        #:if MFC_COMPILER == "LLVMFlang"
#if defined(MFC_OpenMP)
            interface
                function hipMemGetInfo(f, t) bind(C, name='hipMemGetInfo') result(ierr)

                    import :: c_int, c_size_t
                    integer(c_size_t), intent(out) :: f, t
                    integer(c_int)                 :: ierr

                end function hipMemGetInfo
            end interface
            integer(c_size_t) :: f, t
#endif
        #:endif

        free_b = -1_8; total_b = -1_8
        #:if MFC_COMPILER == "LLVMFlang"
#if defined(MFC_OpenMP)
            if (hipMemGetInfo(f, t) /= 0) return
            free_b = int(f, 8); total_b = int(t, 8)
#endif
        #:endif

    end subroutine s_amr_dev_mem

    !> Size the flat store for at least nloc local slots, growing geometrically and never shrinking. Growth preserves the live
    !! columns on the device (no PCIe round trip); the host mirror comes out undefined.
    impure subroutine s_amr_st_reserve(nloc)

        integer, intent(in) :: nloc
        integer             :: oldcap, newcap, i, brlo(3), brhi(3)
        integer             :: c5, i4, k3, j2, i1
        integer             :: blo1, bhi1, blo2, bhi2, blo3, bhi3
        !> device bytes kept free through a growth; without a free-memory query, the largest old store staged on the device
        integer(8), parameter           :: amr_grow_dev_bytes = 4_8*1024_8**3
        integer(8)                      :: st_col_bytes, free_b, total_b
        logical                         :: dev_stage
        real(stp), allocatable          :: tmp(:,:,:,:,:), hstage(:,:,:,:,:)
        type(scalar_field), allocatable :: tmp_br(:)  !< CCE descriptor workaround, see below

        ! Contract: the store is device-authoritative at every call and the host mirror is undefined after a growth. A caller that
        ! just wrote the store on the host (restart read) must push its slot to the device before the next s_amr_alloc_slot.

        if (nloc <= amr_st_cap) return
        oldcap = amr_st_cap
        st_col_bytes = int(mbuf_hi(1) - mbuf_lo(1) + 1, 8)*int(mbuf_hi(2) - mbuf_lo(2) + 1, 8)*int(mbuf_hi(3) - mbuf_lo(3) + 1, &
                           & 8)*int(sys_size, 8)*int(storage_size(0._stp)/8, 8)
        ! Staging one array on the device peaks at 2*newcap - oldcap extra columns: double when that fits in free device memory,
        ! else grow 1.25x (16-slot cap, +8 floor) and stage through the host if the device copy does not fit either.
        call s_amr_dev_mem(free_b, total_b)
        newcap = max(2*oldcap, nloc, 8)
        if (.not. f_fits(newcap)) newcap = max(oldcap + max(min(oldcap/4, 16), 8), nloc)
        dev_stage = f_fits(newcap)
        if (free_b < 0) dev_stage = int(oldcap, 8)*st_col_bytes <= amr_grow_dev_bytes

        blo1 = mbuf_lo(1); bhi1 = mbuf_hi(1)
        blo2 = mbuf_lo(2); bhi2 = mbuf_hi(2)
        blo3 = mbuf_lo(3); bhi3 = mbuf_hi(3)

        #:for ST in ['amr_cons_st', 'amr_stor_st']
            if (.not. dev_stage) then
                ! near-limit fallback: slow (full PCIe both ways) but its device peak is max(old, new)
                $:GPU_UPDATE(host='[' + ST + ']')
                allocate (hstage(mbuf_lo(1):mbuf_hi(1),mbuf_lo(2):mbuf_hi(2),mbuf_lo(3):mbuf_hi(3),1:sys_size,1:oldcap))
                hstage = ${ST}$(:,:,:,:,1:oldcap)
                @:DEALLOCATE(${ST}$)
                @:ALLOCATE(${ST}$(mbuf_lo(1):mbuf_hi(1), mbuf_lo(2):mbuf_hi(2), mbuf_lo(3):mbuf_hi(3), 1:sys_size, 1:newcap))
                ${ST}$(:,:,:,:,1:oldcap) = hstage
                ${ST}$(:,:,:,:,oldcap + 1:newcap) = 0._stp
                deallocate (hstage)
                $:GPU_UPDATE(device='[' + ST + ']')
            else
                if (oldcap > 0) then
                    ! stage the live columns on the device (tmp is device-mapped by @:ALLOCATE); no PCIe traffic
                    @:ALLOCATE(tmp(mbuf_lo(1):mbuf_hi(1), mbuf_lo(2):mbuf_hi(2), mbuf_lo(3):mbuf_hi(3), 1:sys_size, 1:oldcap))
                    $:GPU_PARALLEL_LOOP(collapse=4)
                    do c5 = 1, oldcap
                        do i4 = 1, sys_size
                            do k3 = blo3, bhi3
                                do j2 = blo2, bhi2
                                    do i1 = blo1, bhi1
                                        tmp(i1, j2, k3, i4, c5) = ${ST}$(i1, j2, k3, i4, c5)
                                    end do
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                    @:DEALLOCATE(${ST}$)
                end if
                @:ALLOCATE(${ST}$(mbuf_lo(1):mbuf_hi(1), mbuf_lo(2):mbuf_hi(2), mbuf_lo(3):mbuf_hi(3), 1:sys_size, 1:newcap))
                ! restore the preserved columns and zero the rest on the device; two kernels so the oldcap == 0 path never
                ! references the unallocated tmp
                if (oldcap > 0) then
                    $:GPU_PARALLEL_LOOP(collapse=4)
                    do c5 = 1, oldcap
                        do i4 = 1, sys_size
                            do k3 = blo3, bhi3
                                do j2 = blo2, bhi2
                                    do i1 = blo1, bhi1
                                        ${ST}$(i1, j2, k3, i4, c5) = tmp(i1, j2, k3, i4, c5)
                                    end do
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                    @:DEALLOCATE(tmp)
                end if
                $:GPU_PARALLEL_LOOP(collapse=4)
                do c5 = oldcap + 1, newcap
                    do i4 = 1, sys_size
                        do k3 = blo3, bhi3
                            do j2 = blo2, bhi2
                                do i1 = blo1, bhi1
                                    ${ST}$(i1, j2, k3, i4, c5) = 0._stp
                                end do
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

        amr_st_cap = newcap

        ! the bridge spans a bounded batch of blocks along the last active dimension (amr_br_batch)
        if (.not. allocated(amr_cons_br)) then
            brlo = [mbuf_lo(1), mbuf_lo(2), mbuf_lo(3)]; brhi = [mbuf_hi(1), mbuf_hi(2), mbuf_hi(3)]
            brhi(num_dims) = brlo(num_dims) + amr_br_batch*(brhi(num_dims) - brlo(num_dims) + 1) - 1
            ! CCE descriptor defect, as amr_cg: allocate a local, move_alloc, then map
            allocate (tmp_br(1:sys_size)); call move_alloc(tmp_br, amr_cons_br)
            $:GPU_ENTER_DATA(create='[amr_cons_br]')
            do i = 1, sys_size
                @:ALLOCATE(amr_cons_br(i)%sf(brlo(1):brhi(1), brlo(2):brhi(2), brlo(3):brhi(3)))
                @:ACC_SETUP_SFs(amr_cons_br(i))
            end do
        end if

    contains

        !> Device staging of a growth to cap c leaves amr_grow_dev_bytes free (false when free memory is unknown).
        logical function f_fits(c)

            integer, intent(in) :: c

            f_fits = free_b >= 0 .and. int(2*c - oldcap, 8)*st_col_bytes + amr_grow_dev_bytes <= free_b

        end function f_fits

    end subroutine s_amr_st_reserve

    !> Move block loc's conserved state between the flat store and the bridge over the whole buffered box (all slots share the mbuf
    !! extents).
    #:set BR = 'amr_cons_br(i)%sf(j, k, l)'
    #:set ST = 'amr_cons_st(j, k, l, i, loc)'
    #:for DIR in ['load', 'store']
        #:set LHS = BR if DIR == 'load' else ST
        #:set RHS = ST if DIR == 'load' else BR
        impure subroutine s_amr_br_${DIR}$(loc)

            integer, intent(in) :: loc
            integer             :: i, j, k, l
            integer             :: blo1, bhi1, blo2, bhi2, blo3, bhi3

            blo1 = mbuf_lo(1); bhi1 = mbuf_hi(1)
            blo2 = mbuf_lo(2); bhi2 = mbuf_hi(2)
            blo3 = mbuf_lo(3); bhi3 = mbuf_hi(3)

            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do l = blo3, bhi3
                    do k = blo2, bhi2
                        do j = blo1, bhi1
                            ${LHS}$ = ${RHS}$
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

        end subroutine s_amr_br_${DIR}$
    #:endfor

    !> Face-bounded bridge move for the level>=2 reflux apply: moves only the six one-cell planes (outside cell olo(d)/ohi(d) of
    !! each active dim, over the caller's transverse window) that s_amr_reflux_apply_faces read-modify-writes, skipping zero-weight
    !! planes as the apply does. Bit-identical to a whole-box move.
    #:for DIR in ['load', 'store']
        #:set BF = 'amr_cons_br(i)%sf'
        #:set SF = 'amr_cons_st'
        impure subroutine s_amr_br_${DIR}$_faces(loc, olo, ohi, glo, ghi, woff, w_lo, w_hi)

            integer, intent(in)  :: loc, olo(3), ohi(3), glo(3), ghi(3), woff(3)
            real(wp), intent(in) :: w_lo(3), w_hi(3)
            integer              :: i, g1, g2, ol, oh, wa, wb, gla, gha, glb, ghb
            real(wp)             :: wgt_lo, wgt_hi

            ! per face direction d the transverse dims are (ta, tb); bounds and plane weights are hoisted to scalars since host
            ! dummies in a default(present) region abort CCE OpenACC

            #:for D, TA, TB, IDX in [(1, 2, 3, 'oc, wa + g1, wb + g2'), (2, 1, 3, 'wa + g1, oc, wb + g2'), (3, 1, 2, &
                                      & 'wa + g1, wb + g2, oc')]
                if (amr_dim(${D}$) .and. (w_lo(${D}$) /= 0._wp .or. w_hi(${D}$) /= 0._wp)) then
                    ol = olo(${D}$); oh = ohi(${D}$); wa = woff(${TA}$); wb = woff(${TB}$)
                    gla = glo(${TA}$); gha = ghi(${TA}$); glb = glo(${TB}$); ghb = ghi(${TB}$)
                    wgt_lo = w_lo(${D}$); wgt_hi = w_hi(${D}$)
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = 1, sys_size
                        do g2 = glb, ghb
                            do g1 = gla, gha
                                #:for OC, WT in [('ol', 'wgt_lo'), ('oh', 'wgt_hi')]
                                    if (${WT}$ /= 0._wp) then
                                        #:if DIR == 'load'
                                            ${BF}$(${IDX.replace('oc', OC)}$) = ${SF}$(${IDX.replace('oc', OC)}$, i, loc)
                                        #:else
                                            ${SF}$(${IDX.replace('oc', OC)}$, i, loc) = ${BF}$(${IDX.replace('oc', OC)}$)
                                        #:endif
                                    end if
                                #:endfor
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
            #:endfor

        end subroutine s_amr_br_${DIR}$_faces
    #:endfor

    !> Batched bridge load: every batch member (amr_bat_blk/amr_bat_loc, extents amr_bat_ext) lands at offset (ibm-1)*amr_bat_w
    !! along amr_bat_sd with its ghost shell, in one kernel. Nothing is stored back: s_compute_rhs does not write its cons dummy on
    !! the fine path.
    impure subroutine s_amr_br_load_batch(nb)

        integer, intent(in) :: nb
        integer             :: ibm, i, j, k, l, loc, b1l, b1h, b2l, b2h, b3l, b3h, jj, kk, ll, bs, mb

        bs = buff_size
        b1l = amr_slots(amr_cur)%idwbuff(1)%beg; b1h = amr_slots(amr_cur)%idwbuff(1)%end
        b2l = amr_slots(amr_cur)%idwbuff(2)%beg; b2h = amr_slots(amr_cur)%idwbuff(2)%end
        b3l = amr_slots(amr_cur)%idwbuff(3)%beg; b3h = amr_slots(amr_cur)%idwbuff(3)%end
        $:GPU_PARALLEL_LOOP(collapse=5, private='[loc, jj, kk, ll, mb]', copyin='[nb, nb, b1l, b1h, b2l, b2h, b3l, b3h, bs]')
        do ibm = 1, nb
            do i = 1, sys_size
                do l = b3l, b3h
                    do k = b2l, b2h
                        do j = b1l, b1h
                            mb = min(ibm, nb)  ! the slots past the last member repeat it (finite filler for the stencils)
                            loc = amr_bat_loc(mb)
                            ! padded member (amr_bat_pad): clamp to its own buffered region (finite filler, never read by a real
                            ! stencil)
                            jj = min(j, amr_bat_mext(1, mb) + bs); kk = min(k, amr_bat_mext(2, mb) + bs); ll = min(l, &
                                     & amr_bat_mext(3, mb) + bs)
                            amr_cons_br(i)%sf(j + amr_bat_off(1, ibm), k + amr_bat_off(2, ibm), l + amr_bat_off(3, &
                                        & ibm)) = amr_cons_st(jj, kk, ll, i, loc)
                        end do
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_br_load_batch

    !> Copy batch member ibm's primitive state from the slab scratch into the block-frame scratch over its buffered extent, for the
    !! per-member IB correction.
    impure subroutine s_amr_bat_member_prim(ibm, src, dst)

        integer, intent(in)                                    :: ibm
        type(scalar_field), dimension(sys_size), intent(in)    :: src
        type(scalar_field), dimension(sys_size), intent(inout) :: dst
        integer                                                :: i, j, k, l, h, o1, o2, o3, b1l, b1h, b2l, b2h, b3l, b3h

        h = amr_bat_blk(ibm)
        o1 = amr_bat_off(1, ibm); o2 = amr_bat_off(2, ibm); o3 = amr_bat_off(3, ibm)
        b1l = amr_slots(h)%idwbuff(1)%beg; b1h = amr_slots(h)%idwbuff(1)%end
        b2l = amr_slots(h)%idwbuff(2)%beg; b2h = amr_slots(h)%idwbuff(2)%end
        b3l = amr_slots(h)%idwbuff(3)%beg; b3h = amr_slots(h)%idwbuff(3)%end
        $:GPU_PARALLEL_LOOP(collapse=4, copyin='[o1, o2, o3, b1l, b1h, b2l, b2h, b3l, b3h]')
        do i = 1, sys_size
            do l = b3l, b3h
                do k = b2l, b2h
                    do j = b1l, b1h
                        dst(i)%sf(j, k, l) = src(i)%sf(j + o1, k + o2, l + o3)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_bat_member_prim

    !> Free the flat store and dense-index maps (idempotent; called from both finalize paths).
    impure subroutine s_amr_st_finalize()

        integer :: i

        if (allocated(amr_loc_of)) deallocate (amr_loc_of, amr_loc_free)
        if (allocated(amr_gs_buf)) then
            @:DEALLOCATE(amr_gs_buf)
        end if
        #:for ST in ['amr_cons_st', 'amr_stor_st']
            if (allocated(${ST}$)) then
                @:DEALLOCATE(${ST}$)
            end if
        #:endfor
        if (allocated(amr_cons_br)) then
            do i = 1, sys_size
                @:ACC_TEARDOWN_SFs(amr_cons_br(i))
                @:DEALLOCATE(amr_cons_br(i)%sf)
            end do
            @:DEALLOCATE(amr_cons_br)
        end if
        if (allocated(amr_scr_prim)) then
            do i = 1, sys_size
                @:ACC_TEARDOWN_SFs(amr_scr_prim(i))
                @:DEALLOCATE(amr_scr_prim(i)%sf)
                @:ACC_TEARDOWN_SFs(amr_scr_rhs(i))
                @:DEALLOCATE(amr_scr_rhs(i)%sf)
            end do
            @:DEALLOCATE(amr_scr_prim)
            @:DEALLOCATE(amr_scr_rhs)
        end if
        if (allocated(amr_scr_prim_blk)) then
            do i = 1, sys_size
                @:ACC_TEARDOWN_SFs(amr_scr_prim_blk(i))
                @:DEALLOCATE(amr_scr_prim_blk(i)%sf)
            end do
            @:DEALLOCATE(amr_scr_prim_blk)
        end if
        amr_st_cap = 0

    end subroutine s_amr_st_finalize

    !> Allocate slot islot's coords and its dense store index (idempotent).
    impure subroutine s_amr_alloc_slot(islot)

        integer, intent(in) :: islot

        if (amr_slot_live(islot)) return
        ! recycle a freed local index if one is available, else extend the dense range
        if (amr_loc_nfree > 0) then
            amr_loc_of(islot) = amr_loc_free(amr_loc_nfree)
            amr_loc_nfree = amr_loc_nfree - 1
        else
            amr_loc_n = amr_loc_n + 1
            amr_loc_of(islot) = amr_loc_n
        end if
        amr_slots(islot)%amr_ref_ratio = amr_ref_ratio
        amr_slots(islot)%buff_size = buff_size
        #:for D, X in [(1, 'x'), (2, 'y'), (3, 'z')]
            if (amr_dim(${D}$)) allocate (amr_slots(islot)%${X}$_cb(-1:max_f(${D}$)), amr_slots(islot)%${X}$_cc(0:max_f(${D}$)), &
                & amr_slots(islot)%d${X}$(0:max_f(${D}$)))
        #:endfor
        amr_slot_live(islot) = .true.
        call s_amr_st_reserve(amr_loc_n)

    end subroutine s_amr_alloc_slot

    !> Free slot islot's per-block field arrays (inverse of s_amr_alloc_slot). Idempotent (no-op if not live).
    impure subroutine s_amr_free_slot(islot)

        integer, intent(in) :: islot

        if (.not. amr_slot_live(islot)) return
        if (amr_loc_of(islot) > 0) then
            amr_loc_nfree = amr_loc_nfree + 1
            amr_loc_free(amr_loc_nfree) = amr_loc_of(islot)
            amr_loc_of(islot) = 0
        end if
        #:for X in ['x', 'y', 'z']
            if (allocated(amr_slots(islot)%${X}$_cb)) deallocate (amr_slots(islot)%${X}$_cb, amr_slots(islot)%${X}$_cc, &
                & amr_slots(islot)%d${X}$)
        #:endfor
        amr_slot_live(islot) = .false.

    end subroutine s_amr_free_slot

    !> Reconcile allocated slots to current ownership: allocate every active block this rank owns, free the rest. Call after
    !! ownership is set; regrid must allocate its transient slots before this call.
    impure subroutine s_amr_reconcile_slots()

        integer :: k
        logical :: needed

        do k = 1, amr_max_blocks
            needed = k <= amr_num_blocks
            if (needed) needed = amr_block_owner(k) == proc_rank
            if (needed) then
                call s_amr_alloc_slot(k)
            else
                call s_amr_free_slot(k)
            end if
        end do
        ! every reader of a stale slot has finished by here, which is what makes the renumbering safe
        call s_amr_compact_store()
        amr_mesh_epoch = amr_mesh_epoch + 1  ! local slot indices may have been renumbered: plans that baked them are stale

    end subroutine s_amr_reconcile_slots

end module m_amr_store
