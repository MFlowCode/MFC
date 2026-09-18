!>
!!@file
!!@brief Contains module m_amr_store

#! AMD OpenMP lane: assert allocatables present on every kernel here (see OMP_DEFAULT_STR). A conditionally allocated module
#! array a kernel names launches only under its allocation's own condition (sw_jac/jac: igr); a kernel naming an unallocated
#! array aborts. Keep it so.
#:set MFC_OMP_PRESENT_ALLOCATABLE = True
#:include 'macros.fpp'

!> @brief Flat per-block field store: dense slot indexing, reservation, alloc/free and the prim/bridge loaders.
module m_amr_store

#ifdef MFC_MPI
    use mpi  !< MPI-IO for the parallel_io AMR restart file
#endif

    use m_derived_types  ! scalar_field, t_box, int_bounds_info
    use m_global_parameters
    use m_mpi_proxy, only: s_mpi_abort  ! @:ASSERT expands to it
    use m_variables_conversion, only: s_convert_species_to_mixture_variables_kernel, s_compute_pressure, enforce_density_floor_vc
    use m_phase_timing
    use m_amr_xchg_audit  ! per-call-site accounting of every AMR p2p transfer (s_xa_rec + XA_* site ids)
    use m_amr_state
    use m_amr_distribution

    implicit none

    private
    public :: s_amr_alloc_slot, s_amr_alloc_slot_stash, s_amr_bat_member_prim, s_amr_br_load, s_amr_br_load_batch, &
        & s_amr_br_load_faces, s_amr_br_store, s_amr_br_store_faces, s_amr_copy_fine_fields, s_amr_free_slot, s_amr_alloc_pool, &
        & s_amr_free_pool, s_amr_set_mbuf, s_amr_init_swap_buffers, s_amr_free_swap_buffers, s_amr_prereserve_stash, &
        & s_amr_reconcile_slots, s_amr_sync_grid_state_to_device

contains

    !> Push the (host-side) global grid state to its device copies after a swap/restore. m/n/p, idwint/idwbuff, and the coordinate
    !! arrays are GPU_DECLARE'd; kernels read the device copies. No-op on CPU.
    impure subroutine s_amr_sync_grid_state_to_device()

        $:GPU_UPDATE(device='[m, n, p, idwint, idwbuff]')
        #:for D, X, E in [(1, 'x', 'm'), (2, 'y', 'n'), (3, 'z', 'p')]
            if (amr_dim(${D}$)) then
                $:GPU_UPDATE(device='[' + X + '_cb, ' + X + '_cc, d' + X + ']')
            end if
        #:endfor

    end subroutine s_amr_sync_grid_state_to_device

    !> Device copy amr_cons_st -> amr_stor_st over [b1:e1, b2:e2, b3:e3] for all sys_size fields (RK step-entry backup).
    impure subroutine s_amr_copy_fine_fields(loc, b1, e1, b2, e2, b3, e3)

        integer, intent(in) :: loc  !< flat-store slot: source (amr_cons_st) and destination (amr_stor_st) are the same block
        integer, intent(in) :: b1, e1, b2, e2, b3, e3
        integer             :: i, fi, fj, fk

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

    !> Allocate slot islot's per-block field arrays (coords + the device-resident field vectors), sized to the max buffered block.
    !! Idempotent (no-op if already live). The global amr_cg is not per-slot and stays in init/finalize.
    !> Allocate/reset the dense local-index maps. Called from both pool-allocation sites (s_initialize_amr_module and
    !! s_l0_tiles_init) because pure-L0 mode (amr = F) returns early from the former yet still calls s_amr_alloc_slot. Idempotent so
    !! either order is safe.
    !> The block-metadata pool for amr_max_blocks slots: per-slot geometry, ownership and level tables, the SFC cuts and the
    !! per-block overlap-list counts (the 2D rank lists are sized to the computed max overlap in s_amr_build_seam_pairs).
    impure subroutine s_amr_alloc_pool()

        allocate (amr_slots(1:amr_max_blocks))
        call s_amr_loc_index_init()
        allocate (amr_region_lo_all(3, amr_max_blocks), amr_region_hi_all(3, amr_max_blocks))
        allocate (amr_isect_lo_all(3, amr_max_blocks), amr_isect_hi_all(3, amr_max_blocks))
        allocate (amr_owns_all(amr_max_blocks), amr_block_owner(amr_max_blocks), amr_block_level(amr_max_blocks))
        allocate (amr_owner_cut(0:num_procs - 1)); amr_owner_cut = -1_8
        allocate (amr_fine_cut(0:num_procs - 1,1:max(amr_max_level, 1))); amr_fine_cut = -1_8
        allocate (amr_ovl_gather_n(amr_max_blocks), amr_ovl_scatter_n(amr_max_blocks))
        allocate (amr_slot_live(amr_max_blocks)); amr_slot_live = .false.
        amr_region_lo_all = 0; amr_region_hi_all = 0; amr_isect_lo_all = 0; amr_isect_hi_all = 0; amr_owns_all = .false.
        amr_block_owner = 0

    end subroutine s_amr_alloc_pool

    !> Bounce buffers for the copy-based coordinate swap (GPU-safe; same bounds as the base-level global arrays, which are sized on
    !! *_alloc - these are whole-array assigned to/from x_cb etc., so the shapes must agree).
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
        deallocate (amr_ovl_gather_n, amr_ovl_scatter_n, amr_slots, amr_block_owner, amr_block_level, amr_owner_cut, amr_fine_cut)
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

    impure subroutine s_amr_loc_index_init()

        if (.not. allocated(amr_loc_of)) allocate (amr_loc_of(1:amr_max_blocks))
        if (.not. allocated(amr_loc_free)) allocate (amr_loc_free(1:amr_max_blocks))
        amr_loc_of = 0; amr_loc_free = 0; amr_loc_n = 0; amr_loc_nfree = 0

    end subroutine s_amr_loc_index_init

    !> Move one local slot's store data src -> dst on the device, in place within each live store array (no staging copy: a second
    !! store-sized array on the device is a transient that can exhaust device memory). s_amr_compact_store's ascending-source
    !! ordering guarantees dst's previous contents are already consumed. The host copy goes stale, which is the store's normal state
    !! between rebuilds (device-authoritative; host readers pull per slot).
    impure subroutine s_amr_st_move_slot(src, dst)

        integer, intent(in) :: src, dst
        integer             :: i, j, k, l

        #:for ST in ['amr_cons_st', 'amr_stor_st']
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do l = mbuf_lo(3), mbuf_hi(3)
                    do k = mbuf_lo(2), mbuf_hi(2)
                        do j = mbuf_lo(1), mbuf_hi(1)
                            ${ST}$(j, k, l, i, dst) = ${ST}$(j, k, l, i, src)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        #:endfor

    end subroutine s_amr_st_move_slot

    !> Re-densify the flat store's local index space after every reconcile. Slot indices are global, so at np>=2 a rank's owned set
    !! is a shifting SFC window plus received migration slots; without re-densification `amr_loc_n` would ratchet upward over the
    !! run until the device ran out of memory at fixed per-rank work. Renumbering every reconcile pins loc_n to the live count, so
    !! the capacity high-water plateaus at the rebuild transient (old + new block generations coexist mid-rebuild, ~2x live) instead
    !! of growing without bound. The allocation itself is not shrunk. The moves are device-side and in place, processed in ascending
    !! old-index order, which makes in-place safe: each destination (the rank of its source among the live indices) is <= its
    !! source, and every pending source lies above the current destination.
    !!
    !! Called at the end of s_amr_reconcile_slots, where every reader is finished: the rebuild's overlap carry-forward
    !! (which reads amr_stor_st at the old blocks' local indices) has completed, and the next rebuild has not started.
    impure subroutine s_amr_compact_store()

        integer              :: k, v, newloc
        integer, allocatable :: inv(:)

        if (amr_st_cap <= 0) return
        ! invert the map: inv(local index) = global slot, 0 if free. The walk covers the L0 tile prefix too; those slots
        ! hold store data even though the reconcile walk skips them.
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
#ifdef MFC_DEBUG
        if (newloc < amr_loc_n) write (0, '(A,I0,A,I0,A,I0,A,I0)') '[amr-compact] rank ', proc_rank, ' loc_n ', amr_loc_n, &
            & ' -> ', newloc, '  cap ', amr_st_cap
#endif
        amr_loc_n = newloc
        amr_loc_nfree = 0  ! every recycled index is invalid after renumbering

    end subroutine s_amr_compact_store

    !> Size the flat store for at least nloc local slots, growing 1.25x and never shrinking, so a run pays few reallocations no
    !! matter how the block count churns (with the index space re-densified every reconcile, growth fires only on a new
    !! rebuild-transient high-water). Growth must preserve the live blocks' fields, and those live on the device, so each array
    !! stages through a device-side temporary (two on-device copies, no PCIe round trip). No reader depends on the host mirror
    !! through a grow.
    impure subroutine s_amr_st_reserve(nloc)

        integer, intent(in) :: nloc
        integer             :: oldcap, newcap, i, brlo(3), brhi(3)
        integer             :: c5, i4, k3, j2, i1
        !> device-native staging transiently holds old + tmp columns on the device, and growth fires at the memory high-water mark,
        !! so the transient itself is budgeted: stage on-device while the extra copy stays under amr_grow_dev_bytes, and route a
        !! larger store to the slower host path, whose device peak is max(old, new). A column count would not do: a column's bytes
        !! scale with the block cap, so the same count means wildly different bytes. Known exposure: 4 GiB is 6% of a 64 GB GCD but
        !! 25% of a 16 GB card, and no portable free-memory query exists here; revisit if a small-card production target appears.
        integer(8), parameter           :: amr_grow_dev_bytes = 4_8*1024_8**3
        integer(8)                      :: st_col_bytes
        real(stp), allocatable          :: tmp(:,:,:,:,:), hstage(:,:,:,:,:)
        type(scalar_field), allocatable :: tmp_br(:)  !< CCE descriptor workaround, see below

        ! Contract: the store is device-authoritative at every call; growth preserves the device contents only, and the host
        ! mirror comes out of a growth undefined (host readers pull per slot before reading, the store's normal state between
        ! rebuilds anyway, cf. s_amr_compact_store). A caller that has just written the store on the host (restart read) must
        ! still push its slot to the device before the next s_amr_alloc_slot, or that host data is lost.

        if (nloc <= amr_st_cap) return
        oldcap = amr_st_cap
        ! grow 1.25x with the increment capped at 16 slots: a proportional increment is itself store-scaled, and at a large cap
        ! the +25% transient is what tips a near-limit device over. The +8 floor keeps early growth cheap when oldcap is tiny.
        newcap = max(oldcap + max(min(oldcap/4, 16), 8), nloc)

        #:for ST in ['amr_cons_st', 'amr_stor_st']
            st_col_bytes = int(mbuf_hi(1) - mbuf_lo(1) + 1, 8)*int(mbuf_hi(2) - mbuf_lo(2) + 1, &
                               & 8)*int(mbuf_hi(3) - mbuf_lo(3) + 1, 8)*int(sys_size, 8)*int(storage_size(0._stp)/8, 8)
            if (int(oldcap, 8)*st_col_bytes > amr_grow_dev_bytes) then
                ! near-limit fallback: the device-native staging below transiently holds old + tmp = 2*oldcap columns
                ! on the device, and growth fires exactly at the memory high-water mark. Above the threshold, take the
                ! host round trip: slow (full PCIe both ways) but its device peak is max(old, new).
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
                            do k3 = mbuf_lo(3), mbuf_hi(3)
                                do j2 = mbuf_lo(2), mbuf_hi(2)
                                    do i1 = mbuf_lo(1), mbuf_hi(1)
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
                ! restore the preserved columns and zero the rest, both on the device; the host mirror stays undefined
                ! (see the contract above - every host reader pulls its slot first). Two kernels so the zero-only path
                ! (oldcap == 0) never references the unallocated tmp.
                if (oldcap > 0) then
                    $:GPU_PARALLEL_LOOP(collapse=4)
                    do c5 = 1, oldcap
                        do i4 = 1, sys_size
                            do k3 = mbuf_lo(3), mbuf_hi(3)
                                do j2 = mbuf_lo(2), mbuf_hi(2)
                                    do i1 = mbuf_lo(1), mbuf_hi(1)
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
                        do k3 = mbuf_lo(3), mbuf_hi(3)
                            do j2 = mbuf_lo(2), mbuf_hi(2)
                                do i1 = mbuf_lo(1), mbuf_hi(1)
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

        ! the bridge spans a bounded batch of blocks along the last active dimension (amr_br_batch), so one s_compute_rhs call
        ! can advance a whole batch instead of one block; it rides the same pool lifetime
        if (.not. allocated(amr_cons_br)) then
            brlo = [mbuf_lo(1), mbuf_lo(2), mbuf_lo(3)]; brhi = [mbuf_hi(1), mbuf_hi(2), mbuf_hi(3)]
            brhi(num_dims) = brlo(num_dims) + amr_br_batch*(brhi(num_dims) - brlo(num_dims) + 1) - 1
            ! Same CCE descriptor defect as amr_cg / amr_scr_prim: a bare module-scope derived-type allocatable must be given a
            ! valid descriptor by allocating a local and handing it over with move_alloc, then mapped.
            allocate (tmp_br(1:sys_size)); call move_alloc(tmp_br, amr_cons_br)
            $:GPU_ENTER_DATA(create='[amr_cons_br]')
            do i = 1, sys_size
                @:ALLOCATE(amr_cons_br(i)%sf(brlo(1):brhi(1), brlo(2):brhi(2), brlo(3):brhi(3)))
                @:ACC_SETUP_SFs(amr_cons_br(i))
            end do
        end if

    end subroutine s_amr_st_reserve

    !> Move block loc's conserved state between the flat store and the bridge. One kernel each way over the whole buffered box:
    !! every slot's arrays carry the same mbuf extents.
    #:set BR = 'amr_cons_br(i)%sf(j, k, l)'
    #:set ST = 'amr_cons_st(j, k, l, i, loc)'
    #:for DIR in ['load', 'store']
        #:set LHS = BR if DIR == 'load' else ST
        #:set RHS = ST if DIR == 'load' else BR
        impure subroutine s_amr_br_${DIR}$(loc)

            integer, intent(in) :: loc
            integer             :: i, j, k, l

            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do l = mbuf_lo(3), mbuf_hi(3)
                    do k = mbuf_lo(2), mbuf_hi(2)
                        do j = mbuf_lo(1), mbuf_hi(1)
                            ${LHS}$ = ${RHS}$
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

        end subroutine s_amr_br_${DIR}$
    #:endfor

    !> Face-bounded bridge move, for the level>=2 reflux apply only. s_amr_br_load/_store above move the whole buffered box, but
    !! s_amr_reflux_apply_faces only read-modify-writes six one-cell-thick planes: the outside cell olo(d)/ohi(d) of each active
    !! dim, over the transverse window the caller passes. These twins move the same cells the apply touches and nothing else, so the
    !! result is bit-identical to a whole-box move: the copies are stp->stp with no conversion, the apply's arithmetic is untouched,
    !! and a plane whose weight is zero is skipped by the apply too, so it is neither read nor written either way.
    #:for DIR in ['load', 'store']
        #:set BF = 'amr_cons_br(i)%sf'
        #:set SF = 'amr_cons_st'
        impure subroutine s_amr_br_${DIR}$_faces(loc, olo, ohi, glo, ghi, woff, w_lo, w_hi)

            integer, intent(in)  :: loc, olo(3), ohi(3), glo(3), ghi(3), woff(3)
            real(wp), intent(in) :: w_lo(3), w_hi(3)
            integer              :: i, g1, g2, ol, oh, wa, wb, gla, gha, glb, ghb

            ! per face direction d the transverse dims are (ta, tb); bounds hoisted to scalars for the device region

            #:for D, TA, TB, IDX in [(1, 2, 3, 'oc, wa + g1, wb + g2'), (2, 1, 3, 'wa + g1, oc, wb + g2'), (3, 1, 2, &
                                      & 'wa + g1, wb + g2, oc')]
                if (amr_dim(${D}$) .and. (w_lo(${D}$) /= 0._wp .or. w_hi(${D}$) /= 0._wp)) then
                    ol = olo(${D}$); oh = ohi(${D}$); wa = woff(${TA}$); wb = woff(${TB}$)
                    gla = glo(${TA}$); gha = ghi(${TA}$); glb = glo(${TB}$); ghb = ghi(${TB}$)
                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = 1, sys_size
                        do g2 = glb, ghb
                            do g1 = gla, gha
                                #:for OC, WT in [('ol', 'w_lo'), ('oh', 'w_hi')]
                                    if (${WT}$(${D}$) /= 0._wp) then
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

    !> Batched bridge load: every member of the current batch (amr_bat_blk/amr_bat_loc, extents amr_bat_ext) lands in the bridge at
    !! offset (ibm-1)*amr_bat_w along amr_bat_sd with its ghost shell, in one kernel. Only the member's own buffered box is moved
    !! (the solver never reads outside it), and nothing is stored back: s_compute_rhs does not write its cons dummy on the fine path
    !! (the buffer fill it would come through is a no-op inside the fine advance, and every feature that writes it is excluded by
    !! the validator).
    impure subroutine s_amr_br_load_batch(nb)

        integer, intent(in) :: nb
        integer             :: ibm, i, j, k, l, loc, o1, o2, o3, b1l, b1h, b2l, b2h, b3l, b3h, jj, kk, ll, bs

        bs = buff_size
        o1 = 0; o2 = 0; o3 = 0
        select case (amr_bat_sd)
        case (1); o1 = amr_bat_w
        case (2); o2 = amr_bat_w
        case (3); o3 = amr_bat_w
        end select
        b1l = amr_slots(amr_cur)%idwbuff(1)%beg; b1h = amr_slots(amr_cur)%idwbuff(1)%end
        b2l = amr_slots(amr_cur)%idwbuff(2)%beg; b2h = amr_slots(amr_cur)%idwbuff(2)%end
        b3l = amr_slots(amr_cur)%idwbuff(3)%beg; b3h = amr_slots(amr_cur)%idwbuff(3)%end
        $:GPU_PARALLEL_LOOP(collapse=5, private='[loc, jj, kk, ll]', copyin='[nb, o1, o2, o3, b1l, b1h, b2l, b2h, b3l, b3h, bs]')
        do ibm = 1, nb
            do i = 1, sys_size
                do l = b3l, b3h
                    do k = b2l, b2h
                        do j = b1l, b1h
                            loc = amr_bat_loc(ibm)
                            ! a padded member (amr_bat_pad): clamp to its own buffered region, so the padding holds its
                            ! outermost ghost values (finite and physical, never read by a real cell's stencil)
                            jj = min(j, amr_bat_mext(1, ibm) + bs); kk = min(k, amr_bat_mext(2, ibm) + bs); ll = min(l, &
                                     & amr_bat_mext(3, ibm) + bs)
                            amr_cons_br(i)%sf(j + (ibm - 1)*o1, k + (ibm - 1)*o2, l + (ibm - 1)*o3) = amr_cons_st(jj, kk, ll, i, &
                                        & loc)
                        end do
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_br_load_batch

    !> Copy batch member ibm's primitive state out of the slab scratch (members stacked amr_bat_w apart along amr_bat_sd) into the
    !! block-frame scratch over the member's own buffered extent, for the per-member IB correction after the batch RK update.
    impure subroutine s_amr_bat_member_prim(ibm, src, dst)

        integer, intent(in)                                    :: ibm
        type(scalar_field), dimension(sys_size), intent(in)    :: src
        type(scalar_field), dimension(sys_size), intent(inout) :: dst
        integer                                                :: i, j, k, l, h, o1, o2, o3, b1l, b1h, b2l, b2h, b3l, b3h

        h = amr_bat_blk(ibm)
        o1 = 0; o2 = 0; o3 = 0
        select case (amr_bat_sd)
        case (1); o1 = (ibm - 1)*amr_bat_w
        case (2); o2 = (ibm - 1)*amr_bat_w
        case (3); o3 = (ibm - 1)*amr_bat_w
        end select
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

    !> Free the flat store and the dense-index maps. Mirrors s_amr_loc_index_init: called from both finalize paths, because either
    !! pool-allocation site can have created them. Idempotent.
    impure subroutine s_amr_st_finalize()

        integer :: i

        if (allocated(amr_loc_of)) deallocate (amr_loc_of, amr_loc_free)
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

    !> One-shot pre-reserve before a batch of s_amr_alloc_slot_stash calls: grows the store at most once, to the batch's exact final
    !! size, instead of once per 8-16 incrementally allocated slots. Every growth restages the whole store (s_amr_st_reserve), so a
    !! migration wave allocating tens of replica slots would pay that restage repeatedly without this. getk(k) marks dense
    !! fine-block index k (slot f_l0_slot(k)) for a coming stash alloc.
    impure subroutine s_amr_prereserve_stash(getk, nblk)

        integer, intent(in) :: nblk
        logical, intent(in) :: getk(nblk)
        integer             :: i, nneed

        nneed = 0
        do i = 1, nblk
            if (getk(i) .and. .not. amr_slot_live(f_l0_slot(i))) nneed = nneed + 1
        end do
        ! mirrors the alloc loop exactly: the first amr_loc_nfree allocs pop the recycle stack and do
        ! not raise amr_loc_n; only the remainder grow
        call s_amr_st_reserve(amr_loc_n + max(nneed - amr_loc_nfree, 0))

    end subroutine s_amr_prereserve_stash

    !> Assign slot islot a dense store index without the per-block field arrays: a migration replica only ever has its amr_stor_st
    !! slot written (receive-unpack) and read (the rebuild's overlap carry-forward), so per-slot q_prim/rhs would roughly double its
    !! device cost across the replica set of a migration-heavy regrid. s_amr_alloc_slot upgrades a stash-only slot in place when the
    !! same global slot becomes an owned block; s_amr_free_slot handles both flavors.
    impure subroutine s_amr_alloc_slot_stash(islot)

        integer, intent(in) :: islot

        if (amr_slot_live(islot)) return
        if (amr_loc_nfree > 0) then
            amr_loc_of(islot) = amr_loc_free(amr_loc_nfree)
            amr_loc_nfree = amr_loc_nfree - 1
        else
            amr_loc_n = amr_loc_n + 1
            amr_loc_of(islot) = amr_loc_n
        end if
        amr_slot_live(islot) = .true.
        call s_amr_st_reserve(amr_loc_n)

    end subroutine s_amr_alloc_slot_stash

    impure subroutine s_amr_alloc_slot(islot)

        integer, intent(in) :: islot
        integer             :: i

        ! full-vs-stash discriminator is the grid arrays: every full slot (tile or fine) has x_cb; a stash-only slot has none
        ! (fine slots do not carry q_prim; see the pooled scratch amr_scr_prim/amr_scr_rhs)

        if (amr_slot_live(islot) .and. allocated(amr_slots(islot)%x_cb)) return
        ! recycle a freed local index if one is available, else extend the dense range; a live stash-only slot upgrading to a
        ! full one keeps the index (and so the stor data) it already holds
        if (.not. amr_slot_live(islot)) then
            if (amr_loc_nfree > 0) then
                amr_loc_of(islot) = amr_loc_free(amr_loc_nfree)
                amr_loc_nfree = amr_loc_nfree - 1
            else
                amr_loc_n = amr_loc_n + 1
                amr_loc_of(islot) = amr_loc_n
            end if
        end if
        amr_slots(islot)%amr_ref_ratio = amr_ref_ratio
        amr_slots(islot)%buff_size = buff_size
        #:for D, X in [(1, 'x'), (2, 'y'), (3, 'z')]
            if (amr_dim(${D}$)) allocate (amr_slots(islot)%${X}$_cb(-1:max_f(${D}$)), amr_slots(islot)%${X}$_cc(0:max_f(${D}$)), &
                & amr_slots(islot)%d${X}$(0:max_f(${D}$)))
        #:endfor
        ! pooled scratch: fine blocks advance through the shared scratch (amr_scr_prim/amr_scr_rhs); the fused advance leaves
        ! no cross-block q_prim/rhs lifetime. L0 tile slots are the exception: all owned tiles' rhs coexist across the
        ! MPI-synchronized reflux point (s_l0_add_reflux_to_tiles between the whole-set RHS and RK passes), and a tile's q_prim
        ! written by the RHS pass is read in the later RK pass (IB correction), so tiles keep per-slot rhs always and per-slot
        ! q_prim exactly when s_compute_rhs's copy-out gate writes it (m_rhs.fpp end-of-rhs gate).
        if (islot <= l0_slot_off) then
            @:ALLOCATE(amr_slots(islot)%rhs(1:sys_size))
            if (run_time_info .or. probe_wrt .or. ib .or. bubbles_lagrange) then
                @:ALLOCATE(amr_slots(islot)%q_prim(1:sys_size))
            end if
            do i = 1, sys_size
                ! rhs is ghost-inclusive (mbuf); igr widens to -1:+1 per dim including collapsed ones (coarse rhs_vf is -1:m+1 etc.)
                if (igr) then
                    @:ALLOCATE(amr_slots(islot)%rhs(i)%sf(mbuf_lo(1):mbuf_hi(1), min(mbuf_lo(2), -1):max(mbuf_hi(2), 1), &
                               & min(mbuf_lo(3), -1):max(mbuf_hi(3), 1)))
                else
                    @:ALLOCATE(amr_slots(islot)%rhs(i)%sf(mbuf_lo(1):mbuf_hi(1), mbuf_lo(2):mbuf_hi(2), mbuf_lo(3):mbuf_hi(3)))
                end if
                @:ACC_SETUP_SFs(amr_slots(islot)%rhs(i))
                if (allocated(amr_slots(islot)%q_prim)) then
                    @:ALLOCATE(amr_slots(islot)%q_prim(i)%sf(mbuf_lo(1):mbuf_hi(1), mbuf_lo(2):mbuf_hi(2), mbuf_lo(3):mbuf_hi(3)))
                    @:ACC_SETUP_SFs(amr_slots(islot)%q_prim(i))
                end if
            end do
        end if
        amr_slot_live(islot) = .true.
        call s_amr_st_reserve(amr_loc_n)

    end subroutine s_amr_alloc_slot

    !> Free slot islot's per-block field arrays (inverse of s_amr_alloc_slot). Idempotent (no-op if not live).
    impure subroutine s_amr_free_slot(islot)

        integer, intent(in) :: islot
        integer             :: i

        if (.not. amr_slot_live(islot)) return
        if (amr_loc_of(islot) > 0) then
            amr_loc_nfree = amr_loc_nfree + 1
            amr_loc_free(amr_loc_nfree) = amr_loc_of(islot)
            amr_loc_of(islot) = 0
        end if
        ! Undo each field's ACC_SETUP_SFs (Cray descriptor + %sf copyin) before the @:DEALLOCATE. Cray 'exit data delete'
        ! decrements the ref count, so the lone @:DEALLOCATE would leave the descriptor and the ACC_SETUP %sf ref dangling; the
        ! leaked host address is later reused (e.g. by Gs_rs at restart), tripping a Cray "Error placing / already present"
        ! present-table crash (gpu-acc). A stash-only slot (s_amr_alloc_slot_stash) has none of these arrays, only the index
        ! bookkeeping above.
        if (allocated(amr_slots(islot)%q_prim)) then
            do i = 1, sys_size
                @:ACC_TEARDOWN_SFs(amr_slots(islot)%q_prim(i))
                @:DEALLOCATE(amr_slots(islot)%q_prim(i)%sf)
            end do
            @:DEALLOCATE(amr_slots(islot)%q_prim)
        end if
        if (allocated(amr_slots(islot)%rhs)) then
            do i = 1, sys_size
                @:ACC_TEARDOWN_SFs(amr_slots(islot)%rhs(i))
                @:DEALLOCATE(amr_slots(islot)%rhs(i)%sf)
            end do
            @:DEALLOCATE(amr_slots(islot)%rhs)
        end if
        #:for X in ['x', 'y', 'z']
            if (allocated(amr_slots(islot)%${X}$_cb)) deallocate (amr_slots(islot)%${X}$_cb, amr_slots(islot)%${X}$_cc, &
                & amr_slots(islot)%d${X}$)
        #:endfor
        amr_slot_live(islot) = .false.

    end subroutine s_amr_free_slot

    !> Reconcile the allocated per-slot field arrays to the current ownership: allocate every active block this rank owns, free
    !! everything else. Call after ownership is set (init/regrid/restart). A rank ends holding only its owned blocks' fine arrays
    !! (~amr_num_blocks/num_procs of the pool), not all amr_max_blocks. Regrid must alloc its transient (received/old) slots before
    !! this call, since it frees anything not currently owned.
    impure subroutine s_amr_reconcile_slots()

        integer :: k, nliv, nfr, nal, nfree_in
        logical :: needed

        nliv = 0; nfr = 0; nal = 0; nfree_in = amr_loc_nfree
        do k = 1, amr_max_blocks
            ! Skip the L0 tile prefix [1..l0_slot_off]: those slots are owned + sized by s_l0_tiles_init (rr=1 tile geometry), not
            ! by the AMR fine-block reconcile. Without this, at coexist init the tile-prefix owner defaults to 0, so rank 0 would
            ! alloc these slots here with the fine mbuf* sizing; s_amr_alloc_slot is idempotent, so s_l0_build_tile_slot could not
            ! then resize them (a tile-undersizing hazard, benign only while fine mbuf* >= tile). l0_slot_off = 0 with no tiles,
            ! so this is a no-op for pure AMR.
            if (k <= l0_slot_off) cycle
            needed = k <= amr_num_blocks
            if (needed) needed = amr_block_owner(k) == proc_rank
            if (needed) then
                if (.not. amr_slot_live(k)) nal = nal + 1
                call s_amr_alloc_slot(k)
                nliv = nliv + 1
            else
                if (amr_slot_live(k)) nfr = nfr + 1
                call s_amr_free_slot(k)
            end if
        end do
        ! stderr (survives an abort). live = what the run actually needs; loc_n = indices ever handed out;
        ! the gap between them is the leak. stack_in/out shows whether frees accumulate for the next
        ! rebuild to recycle, and 'new' counts allocs that had to extend rather than recycle.
#ifdef MFC_DEBUG
        write (0, '(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') '[amr-recon] rank ', proc_rank, ' live ', nliv, ' loc_n ', amr_loc_n, &
               & ' freed ', nfr, ' newalloc ', nal, ' stack_in ', nfree_in, ' stack_out ', amr_loc_nfree
#endif
        ! invariant: no stash-only slot survives a reconcile; every migration replica was freed (early-free or the walk above)
        ! or upgraded to a full slot by the owned-path alloc. A survivor would reach the solver with no geometry arrays.
        do k = 1, amr_max_blocks
            if (amr_slot_live(k)) then
                @:ASSERT(allocated(amr_slots(k)%x_cb), "a stash-only replica slot survived reconcile")
            end if
        end do
        ! every reader of a stale slot has finished by here (the rebuild's overlap carry-forward is done and the next rebuild
        ! has not started), which is what makes the renumbering safe; see s_amr_compact_store.
        call s_amr_compact_store()
        ! per-rank store capacity is the weak-scaling invariant (device memory = f(live local boxes)); wall time cannot see it,
        ! so report it. live == loc_n after the compaction above; cap - live is the rebuild-transient envelope.
        if (rank_time_wrt) write (0, '(A,I0,A,I0,A,I0)') '[amr-cap] rank ', proc_rank, ' live ', amr_loc_n, ' cap ', amr_st_cap
        amr_mesh_epoch = amr_mesh_epoch + 1  ! local slot indices may have been renumbered: plans that baked them are stale

    end subroutine s_amr_reconcile_slots

end module m_amr_store
