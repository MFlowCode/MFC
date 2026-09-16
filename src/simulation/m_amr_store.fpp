!>
!!@file
!!@brief Contains module m_amr_store

#! AMD OpenMP lane: assert allocatables present on every kernel here (see OMP_DEFAULT_STR). Every conditionally allocated
#! module array a kernel here names launches only under its allocation's own condition (amr_rvw: cyl_coord; sw_jac/jac: igr;
#! amr_cg_pb/mv: do_pbmv; amr_gst_a/b: amr_subcycle; amr_prim_st/amr_bt_*: amr_prim_batch); amr_cg and amr_cons_br/stor_st are
#! allocated before first use. A kernel naming an unallocated array aborts. Keep it so.
#:set MFC_OMP_PRESENT_ALLOCATABLE = True
#:include 'macros.fpp'

!> @brief Flat per-block field store: dense slot indexing, reservation, alloc/free and the prim/bridge loaders.
module m_amr_store

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
    public :: s_amr_alloc_slot, s_amr_alloc_slot_stash, s_amr_bat_member_prim, s_amr_br_load, s_amr_br_load_batch, &
        & s_amr_br_load_faces, s_amr_br_store, s_amr_br_store_faces, s_amr_convert_prim_batch, s_amr_copy_fine_fields, &
        & s_amr_free_slot, s_amr_loc_index_init, s_amr_prereserve_stash, s_amr_prim_load, s_amr_recompute_weno_coefs, &
        & s_amr_reconcile_slots, s_amr_st_finalize, s_amr_sync_grid_state_to_device

contains

    !> Recompute the WENO reconstruction coefficient arrays from the current grid globals (the fine block's after a swap, the coarse
    !! grid's after a restore). s_compute_weno_coefficients reads the live cell-boundary arrays, refreshes uniform_grid, and pushes
    !! its own device updates; the coefficient arrays are sized to m/n/p_alloc at init, which no fine range exceeds.
    impure subroutine s_amr_recompute_weno_coefs()

        type(int_bounds_info) :: is1, is2, is3

        is1%beg = -buff_size; is1%end = m + buff_size
        call s_compute_weno_coefficients(1, is1)
        if (n_glb > 0) then
            is2%beg = -buff_size; is2%end = n + buff_size
            call s_compute_weno_coefficients(2, is2)
        end if
        if (p_glb > 0) then
            is3%beg = -buff_size; is3%end = p + buff_size
            call s_compute_weno_coefficients(3, is3)
        end if

    end subroutine s_amr_recompute_weno_coefs

    !> Push the (host-side) global grid state to its device copies after a swap/restore. m/n/p, idwint/idwbuff, and the coordinate
    !! arrays are GPU_DECLARE'd; kernels read the device copies. No-op on CPU.
    impure subroutine s_amr_sync_grid_state_to_device()

        $:GPU_UPDATE(device='[m, n, p, idwint, idwbuff]')
        $:GPU_UPDATE(device='[x_cb, x_cc, dx]')
        if (n_glb > 0) then
            $:GPU_UPDATE(device='[y_cb, y_cc, dy]')
        end if
        if (p_glb > 0) then
            $:GPU_UPDATE(device='[z_cb, z_cc, dz]')
        end if

    end subroutine s_amr_sync_grid_state_to_device

    !> Device copy amr_cons_st -> amr_stor_st over [b1:e1, b2:e2, b3:e3] for all sys_size fields (RK step-entry backup). Twin
    !! s_amr_backup_pbmv (q<->pb/mv): pb/mv sibling of this step-entry backup; keep them in lockstep.
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

    !> Allocate slot islot's per-block field arrays (coords + the 6 device-resident field vectors + non-poly QBMM side-state), sized
    !! to the max buffered block. Idempotent (no-op if already live). The single QBMM RHS scratch amr_rhs_pb_f/mv_f and the global
    !! amr_cg are not per-slot and stay in init/finalize.
    !> Allocate/reset the dense local-index maps. Called from both pool-allocation sites (s_initialize_amr_module and
    !! s_l0_tiles_init) because pure-L0 mode (amr = F) returns early from the former yet still calls s_amr_alloc_slot. Idempotent so
    !! either order is safe.
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
                do l = mbuf3_lo, mbuf3_hi
                    do k = mbuf2_lo, mbuf2_hi
                        do j = mbuf1_lo, mbuf1_hi
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
        amr_st_hw = newloc  ! let the trip-wire report the post-compaction trajectory

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

        ! Trip-wire on the store trajectory. stderr, because stdout is buffered and lost on abort.

        if (nloc > amr_st_hw) then
            amr_st_hw = nloc
#ifdef MFC_DEBUG
            write (0, '(A,I0,A,I0,A,I0,A,I0)') '[amr-store] rank ', proc_rank, ' NEW high-water nloc ', nloc, ' cap ', &
                   & amr_st_cap, ' recycle-depth ', amr_loc_nfree
#endif
        end if
        if (nloc <= amr_st_cap) return
        oldcap = amr_st_cap
        ! grow 1.25x with the increment capped at 16 slots: a proportional increment is itself store-scaled, and at a large cap
        ! the +25% transient is what tips a near-limit device over. The +8 floor keeps early growth cheap when oldcap is tiny.
        newcap = max(oldcap + max(min(oldcap/4, 16), 8), nloc)

        #:for ST in ['amr_cons_st', 'amr_stor_st']
            st_col_bytes = int(mbuf1_hi - mbuf1_lo + 1, 8)*int(mbuf2_hi - mbuf2_lo + 1, 8)*int(mbuf3_hi - mbuf3_lo + 1, &
                               & 8)*int(sys_size, 8)*int(storage_size(0._stp)/8, 8)
            if (int(oldcap, 8)*st_col_bytes > amr_grow_dev_bytes) then
                ! near-limit fallback: the device-native staging below transiently holds old + tmp = 2*oldcap columns
                ! on the device, and growth fires exactly at the memory high-water mark. Above the threshold, take the
                ! host round trip: slow (full PCIe both ways) but its device peak is max(old, new).
                $:GPU_UPDATE(host='[' + ST + ']')
                allocate (hstage(mbuf1_lo:mbuf1_hi,mbuf2_lo:mbuf2_hi,mbuf3_lo:mbuf3_hi,1:sys_size,1:oldcap))
                hstage = ${ST}$(:,:,:,:,1:oldcap)
                @:DEALLOCATE(${ST}$)
                @:ALLOCATE(${ST}$(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi, 1:sys_size, 1:newcap))
                ${ST}$(:,:,:,:,1:oldcap) = hstage
                ${ST}$(:,:,:,:,oldcap + 1:newcap) = 0._stp
                deallocate (hstage)
                $:GPU_UPDATE(device='[' + ST + ']')
            else
                if (oldcap > 0) then
                    ! stage the live columns on the device (tmp is device-mapped by @:ALLOCATE); no PCIe traffic
                    @:ALLOCATE(tmp(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi, 1:sys_size, 1:oldcap))
                    $:GPU_PARALLEL_LOOP(collapse=4)
                    do c5 = 1, oldcap
                        do i4 = 1, sys_size
                            do k3 = mbuf3_lo, mbuf3_hi
                                do j2 = mbuf2_lo, mbuf2_hi
                                    do i1 = mbuf1_lo, mbuf1_hi
                                        tmp(i1, j2, k3, i4, c5) = ${ST}$(i1, j2, k3, i4, c5)
                                    end do
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                    @:DEALLOCATE(${ST}$)
                end if
                @:ALLOCATE(${ST}$(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi, 1:sys_size, 1:newcap))
                ! restore the preserved columns and zero the rest, both on the device; the host mirror stays undefined
                ! (see the contract above - every host reader pulls its slot first). Two kernels so the zero-only path
                ! (oldcap == 0) never references the unallocated tmp.
                if (oldcap > 0) then
                    $:GPU_PARALLEL_LOOP(collapse=4)
                    do c5 = 1, oldcap
                        do i4 = 1, sys_size
                            do k3 = mbuf3_lo, mbuf3_hi
                                do j2 = mbuf2_lo, mbuf2_hi
                                    do i1 = mbuf1_lo, mbuf1_hi
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
                        do k3 = mbuf3_lo, mbuf3_hi
                            do j2 = mbuf2_lo, mbuf2_hi
                                do i1 = mbuf1_lo, mbuf1_hi
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
            brlo = [mbuf1_lo, mbuf2_lo, mbuf3_lo]; brhi = [mbuf1_hi, mbuf2_hi, mbuf3_hi]
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

        ! prim landing zone + batch metadata: per-stage scratch (rewritten by every s_amr_convert_prim_batch
        ! call), so growth discards contents; no device staging round trip, unlike the stores above.
        if (amr_prim_batch) then
            if (allocated(amr_prim_st)) then
                @:DEALLOCATE(amr_prim_st)
                @:DEALLOCATE(amr_bt_lo)
                @:DEALLOCATE(amr_bt_hi)
                @:DEALLOCATE(amr_bt_on)
            end if
            allocate (amr_prim_st(mbuf1_lo:mbuf1_hi,mbuf2_lo:mbuf2_hi,mbuf3_lo:mbuf3_hi,1:num_vels + 1,1:newcap))
            allocate (amr_bt_lo(3, newcap), amr_bt_hi(3, newcap), amr_bt_on(newcap))
        end if

    end subroutine s_amr_st_reserve

    !> One batched cons->prim conversion over every owned fine block (all levels), straight from the flat cons store into the flat
    !! prim landing zone. Runs once per RK stage after the fill + seam phases; each block's store bytes there are identical to what
    !! its per-block conversion point would read (advances write only their own slots), and the kernel is per-cell with no
    !! reductions, so the result is bit-identical to the per-block path. The cell body below is pinned to
    !! s_convert_conservative_to_primitive_variables (m_variables_conversion.fpp) restricted to the amr_prim_batch gate's configs:
    !! species fractions (s_compute_species_fraction inlined against the store; igr/bubbles_euler excluded by the gate), mixture
    !! properties, velocity + dynamic pressure, and pressure. Change the conversion and this must follow.
    impure subroutine s_amr_convert_prim_batch()

        integer :: g, loc, i, j, k, l, gg
        integer :: nl, nv, b1l, b1h, b2l, b2h, b3l, b3h

        #:if USING_AMD and not MFC_CASE_OPTIMIZATION
            real(wp), dimension(3) :: alpha_K, alpha_rho_K
            real(wp)               :: rhoYks_b(1:10)
        #:else
            real(wp), dimension(num_fluids) :: alpha_K, alpha_rho_K
            real(wp)                        :: rhoYks_b(1:num_species)
        #:endif
        real(wp) :: Re_K(2)
        real(wp) :: rho_K, gamma_K, pi_inf_K, qv_K, dyn_pres_K, alpha_K_sum, pres, T, pmag

        if (amr_loc_n == 0) return
        call s_phase_tic(PH_CVTB)
        amr_bt_on(1:amr_loc_n) = .false.
        call s_amr_refresh_my_blocks()
        do gg = 1, amr_n_my
            g = amr_my_blk(gg)
            if (amr_block_level(g) < 1) cycle
            loc = amr_loc_of(g)
            if (loc <= 0) cycle
            amr_bt_on(loc) = .true.
            do i = 1, 3
                amr_bt_lo(i, loc) = amr_slots(g)%idwbuff(i)%beg
                amr_bt_hi(i, loc) = amr_slots(g)%idwbuff(i)%end
            end do
        end do
        $:GPU_UPDATE(device='[amr_bt_on, amr_bt_lo, amr_bt_hi]')
        ! bounds through local scalars, never GPU_DECLARE'd module state (the CCE-acc stale-device-bounds class)
        nl = amr_loc_n; nv = num_vels
        b1l = mbuf1_lo; b1h = mbuf1_hi; b2l = mbuf2_lo; b2h = mbuf2_hi; b3l = mbuf3_lo; b3h = mbuf3_hi
        $:GPU_PARALLEL_LOOP(collapse=4, private='[alpha_K, alpha_rho_K, Re_K, rhoYks_b, rho_K, gamma_K, pi_inf_K, qv_K, &
                            & dyn_pres_K, alpha_K_sum, pres, T, pmag]', copyin='[nl, nv, b1l, b1h, b2l, b2h, b3l, b3h]')
        do loc = 1, nl
            do l = b3l, b3h
                do k = b2l, b2h
                    do j = b1l, b1h
                        if (.not. amr_bt_on(loc)) cycle
                        if (j < amr_bt_lo(1, loc) .or. j > amr_bt_hi(1, loc) .or. k < amr_bt_lo(2, loc) .or. k > amr_bt_hi(2, &
                            & loc) .or. l < amr_bt_lo(3, loc) .or. l > amr_bt_hi(3, loc)) cycle
                        if (num_fluids == 1) then
                            alpha_rho_K(1) = amr_cons_st(j, k, l, eqn_idx%cont%beg, loc)
                            alpha_K(1) = amr_cons_st(j, k, l, eqn_idx%adv%beg, loc)
                        else
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                alpha_rho_K(i) = amr_cons_st(j, k, l, i, loc)
                                alpha_K(i) = amr_cons_st(j, k, l, eqn_idx%adv%beg + i - 1, loc)
                            end do
                        end if
                        if (mpp_lim) then
                            alpha_K_sum = 0._wp
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                alpha_rho_K(i) = max(0._wp, alpha_rho_K(i))
                                alpha_K(i) = min(max(0._wp, alpha_K(i)), 1._wp)
                                alpha_K_sum = alpha_K_sum + alpha_K(i)
                            end do
                            ! explicit loop, not array syntax: an inline whole-array expression in a target
                            ! region is a per-thread temporary on amdflang
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                alpha_K(i) = alpha_K(i)/max(alpha_K_sum, 1.e-16_wp)
                            end do
                        end if
                        call s_convert_species_to_mixture_variables_kernel(rho_K, gamma_K, pi_inf_K, qv_K, alpha_K, alpha_rho_K, &
                            & Re_K)
                        if (enforce_density_floor_vc) rho_K = max(rho_K, sgm_eps)
                        dyn_pres_K = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, nv
                            amr_prim_st(j, k, l, i, loc) = amr_cons_st(j, k, l, eqn_idx%mom%beg + i - 1, loc)/rho_K
                            dyn_pres_K = dyn_pres_K + 5.e-1_wp*amr_cons_st(j, k, l, eqn_idx%mom%beg + i - 1, loc)*amr_prim_st(j, &
                                & k, l, i, loc)
                        end do
                        pmag = 0._wp
                        call s_compute_pressure(amr_cons_st(j, k, l, eqn_idx%E, loc), amr_cons_st(j, k, l, eqn_idx%alf, loc), &
                                                & dyn_pres_K, pi_inf_K, gamma_K, rho_K, qv_K, rhoYks_b, pres, T, pres_mag=pmag)
                        amr_prim_st(j, k, l, nv + 1, loc) = pres
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        call s_phase_toc(PH_CVTB)

    end subroutine s_amr_convert_prim_batch

    !> Land the current block's batch-computed prim vars (the contiguous mom%beg..E range) from the prim store into the m_rhs
    !! conversion scratch, replacing that block's per-block conversion.
    impure subroutine s_amr_prim_load(q_prim_b, loc)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_b
        integer, intent(in)                                    :: loc
        integer                                                :: i, lb(3)

        ! One launch per var through a plain contiguous array dummy: a scalar_field-array dummy makes the
        ! target region map the derived-type descriptors per launch. The dummy rebases to 1, so offset from
        ! the actual's own lower bounds.

        do i = eqn_idx%mom%beg, eqn_idx%E
            lb = lbound(q_prim_b(i)%sf)
            call s_amr_prim_load_one(q_prim_b(i)%sf, i - eqn_idx%mom%beg + 1, loc, amr_slots(amr_cur)%idwbuff(1)%beg, &
                                     & amr_slots(amr_cur)%idwbuff(1)%end, amr_slots(amr_cur)%idwbuff(2)%beg, &
                                     & amr_slots(amr_cur)%idwbuff(2)%end, amr_slots(amr_cur)%idwbuff(3)%beg, &
                                     & amr_slots(amr_cur)%idwbuff(3)%end, 1 - lb(1), 1 - lb(2), 1 - lb(3))
        end do

    end subroutine s_amr_prim_load

    impure subroutine s_amr_prim_load_one(dst, pv, loc, j1l, j1h, j2l, j2h, j3l, j3h, o1, o2, o3)

        real(stp), dimension(:,:,:), contiguous, intent(inout) :: dst
        integer, intent(in)                                    :: pv, loc, j1l, j1h, j2l, j2h, j3l, j3h, o1, o2, o3
        integer                                                :: j, k, l

        $:GPU_PARALLEL_LOOP(collapse=3, copyin='[pv, loc, j1l, j1h, j2l, j2h, j3l, j3h, o1, o2, o3]')
        do l = j3l, j3h
            do k = j2l, j2h
                do j = j1l, j1h
                    dst(j + o1, k + o2, l + o3) = amr_prim_st(j, k, l, pv, loc)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_prim_load_one

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
                do l = mbuf3_lo, mbuf3_hi
                    do k = mbuf2_lo, mbuf2_hi
                        do j = mbuf1_lo, mbuf1_hi
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
            integer              :: i, g1, g2, ol, oh, w1, w2, w3, gl1, gh1, gl2, gh2, gl3, gh3

            gl1 = glo(1); gh1 = ghi(1); gl2 = glo(2); gh2 = ghi(2); gl3 = glo(3); gh3 = ghi(3)
            w1 = woff(1); w2 = woff(2); w3 = woff(3)

            ! x-faces: transverse (y, z)
            if (w_lo(1) /= 0._wp .or. w_hi(1) /= 0._wp) then
                ol = olo(1); oh = ohi(1)
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, sys_size
                    do g2 = gl3, gh3
                        do g1 = gl2, gh2
                            #:for OC, WT in [('ol', 'w_lo(1)'), ('oh', 'w_hi(1)')]
                                if (${WT}$ /= 0._wp) then
                                    #:if DIR == 'load'
                                        ${BF}$(${OC}$, w2 + g1, w3 + g2) = ${SF}$(${OC}$, w2 + g1, w3 + g2, i, loc)
                                    #:else
                                        ${SF}$(${OC}$, w2 + g1, w3 + g2, i, loc) = ${BF}$(${OC}$, w2 + g1, w3 + g2)
                                    #:endif
                                end if
                            #:endfor
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
            ! y-faces: transverse (x, z)
            if (n_glb > 0 .and. (w_lo(2) /= 0._wp .or. w_hi(2) /= 0._wp)) then
                ol = olo(2); oh = ohi(2)
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, sys_size
                    do g2 = gl3, gh3
                        do g1 = gl1, gh1
                            #:for OC, WT in [('ol', 'w_lo(2)'), ('oh', 'w_hi(2)')]
                                if (${WT}$ /= 0._wp) then
                                    #:if DIR == 'load'
                                        ${BF}$(w1 + g1, ${OC}$, w3 + g2) = ${SF}$(w1 + g1, ${OC}$, w3 + g2, i, loc)
                                    #:else
                                        ${SF}$(w1 + g1, ${OC}$, w3 + g2, i, loc) = ${BF}$(w1 + g1, ${OC}$, w3 + g2)
                                    #:endif
                                end if
                            #:endfor
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
            ! z-faces: transverse (x, y)
            if (p_glb > 0 .and. (w_lo(3) /= 0._wp .or. w_hi(3) /= 0._wp)) then
                ol = olo(3); oh = ohi(3)
                $:GPU_PARALLEL_LOOP(collapse=3)
                do i = 1, sys_size
                    do g2 = gl2, gh2
                        do g1 = gl1, gh1
                            #:for OC, WT in [('ol', 'w_lo(3)'), ('oh', 'w_hi(3)')]
                                if (${WT}$ /= 0._wp) then
                                    #:if DIR == 'load'
                                        ${BF}$(w1 + g1, w2 + g2, ${OC}$) = ${SF}$(w1 + g1, w2 + g2, ${OC}$, i, loc)
                                    #:else
                                        ${SF}$(w1 + g1, w2 + g2, ${OC}$, i, loc) = ${BF}$(w1 + g1, w2 + g2, ${OC}$)
                                    #:endif
                                end if
                            #:endfor
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

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
        if (allocated(amr_prim_st)) then
            @:DEALLOCATE(amr_prim_st)
            @:DEALLOCATE(amr_bt_lo)
            @:DEALLOCATE(amr_bt_hi)
            @:DEALLOCATE(amr_bt_on)
        end if
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
        allocate (amr_slots(islot)%x_cb(-1:max_f1), amr_slots(islot)%x_cc(0:max_f1), amr_slots(islot)%dx(0:max_f1))
        if (n_glb > 0) allocate (amr_slots(islot)%y_cb(-1:max_f2), amr_slots(islot)%y_cc(0:max_f2), amr_slots(islot)%dy(0:max_f2))
        if (p_glb > 0) allocate (amr_slots(islot)%z_cb(-1:max_f3), amr_slots(islot)%z_cc(0:max_f3), amr_slots(islot)%dz(0:max_f3))
        ! pooled scratch: fine blocks advance through the shared scratch (amr_scr_prim/amr_scr_rhs); the fused per-block advance
        ! leaves no cross-block q_prim/rhs lifetime. L0 tile slots are the exception: all owned tiles' rhs coexist across the
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
                    @:ALLOCATE(amr_slots(islot)%rhs(i)%sf(mbuf1_lo:mbuf1_hi, min(mbuf2_lo, -1):max(mbuf2_hi, 1), min(mbuf3_lo, &
                               & -1):max(mbuf3_hi, 1)))
                else
                    @:ALLOCATE(amr_slots(islot)%rhs(i)%sf(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi))
                end if
                @:ACC_SETUP_SFs(amr_slots(islot)%rhs(i))
                if (allocated(amr_slots(islot)%q_prim)) then
                    @:ALLOCATE(amr_slots(islot)%q_prim(i)%sf(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi))
                    @:ACC_SETUP_SFs(amr_slots(islot)%q_prim(i))
                end if
            end do
        end if
        if (qbmm .and. .not. polytropic) then
            #:for PF in ['pb_f', 'mv_f', 'pb_stor', 'mv_stor']
                @:ALLOCATE(amr_slots(islot)%${PF}$%sf(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi, 1:nnode, 1:nb))
                @:ACC_SETUP_SFs(amr_slots(islot)%${PF}$)
            #:endfor
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
        ! decrements
        ! the ref count, so the lone @:DEALLOCATE would leave the descriptor and the ACC_SETUP %sf ref dangling; the leaked host
        ! address is later reused (e.g. by Gs_rs at restart), tripping a Cray "Error placing / already present" present-table crash
        ! (gpu-acc). A stash-only slot (s_amr_alloc_slot_stash) has none of these arrays, only the index bookkeeping above.
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
        if (qbmm .and. .not. polytropic .and. associated(amr_slots(islot)%pb_f%sf)) then
            #:for PF in ['pb_f', 'mv_f', 'pb_stor', 'mv_stor']
                @:ACC_TEARDOWN_SFs(amr_slots(islot)%${PF}$)
                @:DEALLOCATE(amr_slots(islot)%${PF}$%sf)
            #:endfor
        end if
        if (allocated(amr_slots(islot)%x_cb)) deallocate (amr_slots(islot)%x_cb, amr_slots(islot)%x_cc, amr_slots(islot)%dx)
        if (allocated(amr_slots(islot)%y_cb)) deallocate (amr_slots(islot)%y_cb, amr_slots(islot)%y_cc, amr_slots(islot)%dy)
        if (allocated(amr_slots(islot)%z_cb)) deallocate (amr_slots(islot)%z_cb, amr_slots(islot)%z_cc, amr_slots(islot)%dz)
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
