!>
!!@file
!!@brief Contains module m_amr

#! AMD OpenMP lane: assert allocatables present on every kernel here (see OMP_DEFAULT_STR). Every conditionally allocated
#! module array a kernel here names launches only under its allocation's own condition (amr_rvw: cyl_coord; sw_jac/jac: igr;
#! amr_cg_pb/mv: do_pbmv; amr_gst_a/b: amr_subcycle; amr_prim_st/amr_bt_*: amr_prim_batch); amr_cg and amr_cons_br/stor_st are
#! allocated before first use. A kernel naming an unallocated array aborts. Keep it so.
#:set MFC_OMP_PRESENT_ALLOCATABLE = True
#:include 'macros.fpp'

!> @brief Block-structured AMR: up to amr_max_blocks refined blocks (2:1 or 4:1 per amr_ref_ratio), optionally nested to
!! amr_max_level, advanced with the shared solver via grid-state swap and conservatively coupled to each block's parent level (ghost
!! prolongation, Berger-Colella flux reflux, restriction) and dynamic regrid.
module m_amr

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
    use m_amr_store
    use m_amr_exchange
    use m_amr_frame
    use m_amr_transfer
    use m_amr_advance
    use m_amr_l0

    implicit none

    !> Every AMR service module is re-exported here so that "use m_amr" keeps reaching all of them.

contains

    !> Build the static refined level-1 block. No-op unless amr. Called after the level-0 grid (x_cb/dx ready) and time-steppers
    !! (sys_size/buff_size set). Per-slot fine arrays allocated lazily (s_amr_reconcile_slots) - only the blocks a rank owns.
    impure subroutine s_initialize_amr_module()

        integer                         :: i, d
        integer                         :: sidx(3), ext(3), maxc_loc(3), bad_loc, bad_glb, fit_d
        integer                         :: blk_lo(3), blk_hi(3)
        type(scalar_field), allocatable :: tmp_cg(:)

        ! shared-pool layout: tiles are a fixed level-0 prefix; fine blocks follow. Both this init and s_l0_tiles_init read this, so
        ! it runs before the amr early-return below (this routine always executes first, per m_start_up.fpp).

        if (l0_ntile > 0) then
            l0_nt = 1; l0_nt(1) = l0_ntile
            if (n_glb > 0) l0_nt(2) = l0_ntile
            if (p_glb > 0) l0_nt(3) = l0_ntile
            l0_ntiles_tot = num_procs*l0_nt(1)*l0_nt(2)*l0_nt(3)
            l0_slot_off = l0_ntiles_tot
        end if

        if (.not. amr) return

#ifdef MFC_GPU
        amr_fw_dev = rdma_mpi .and. XA_NH == 0
#endif

        ! Batched-conversion gate (see amr_prim_batch's declaration). Off: the batched conversion kernel itself is cheap and
        ! byte-identical, but the per-block prim bridge-loads that land its output in the m_rhs scratch cost more than they save
        ! on the OpenMP-offload host path. The machinery stays for a store-native consumption path that would delete those loads.
        amr_prim_batch = .false.

        ! Fine-block cap = the case amr_max_blocks; the shared pool adds the L0 tile prefix (l0_slot_off, 0 when l0_ntile=0) ahead
        ! of
        ! it, so both AMR fine blocks and any L0 tiles draw from one amr_slots allocation.
        amr_max_fine = amr_max_blocks  ! fine/regrid cap = the case budget
        amr_max_blocks = l0_slot_off + amr_max_fine  ! total shared pool (l0_slot_off=0 when no tiles -> unchanged)

        ! fixed pool of amr_max_blocks slots; init activates exactly one (amr_cur = f_l0_slot(1), the initial fine-block slot);
        ! regrid clusters into up to amr_max_blocks
        allocate (amr_slots(1:amr_max_blocks))
        call s_amr_loc_index_init()
        allocate (amr_region_lo_all(3, amr_max_blocks), amr_region_hi_all(3, amr_max_blocks))
        allocate (amr_isect_lo_all(3, amr_max_blocks), amr_isect_hi_all(3, amr_max_blocks))
        allocate (amr_owns_all(amr_max_blocks))
        allocate (amr_block_owner(amr_max_blocks))
        allocate (amr_owner_cut(0:num_procs - 1)); amr_owner_cut = -1_8
        allocate (amr_fine_cut(0:num_procs - 1,1:max(amr_max_level, 1))); amr_fine_cut = -1_8
        allocate (amr_block_level(amr_max_blocks))
        ! amr_ovl_gather/scatter (the 2D rank lists) are allocated to the computed max overlap in s_amr_build_seam_pairs; only the
        ! per-block counts are sized here.
        allocate (amr_ovl_gather_n(amr_max_blocks), amr_ovl_scatter_n(amr_max_blocks))
        amr_region_lo_all = 0; amr_region_hi_all = 0; amr_isect_lo_all = 0; amr_isect_hi_all = 0; amr_owns_all = .false.
        amr_block_owner = 0
        amr_block_level = 1  ! init default (level-1); regrid re-tags each block's level for nesting
        amr_num_levels = 1
        amr_num_blocks = f_l0_slot(1)
        amr_cur = f_l0_slot(1)

        ! fine-level load balance is capped at min(num blocks, amr_max_blocks) ranks: the SFC map spreads whole blocks, so with
        ! fewer
        ! blocks than ranks some ranks own no fine work. Warn when the pool itself is the limit (raise amr_max_blocks).
        if (proc_rank == 0 .and. num_procs > amr_max_fine) then
            print '(A,I0,A,I0,A)', ' [amr] WARNING: amr_max_blocks (', amr_max_blocks, ') < num_procs (', num_procs, &
                & '): the fine level can occupy at most amr_max_blocks ranks - raise amr_max_blocks for better fine-level balance'
        end if

        ! Every fine block advances at the coarse dt, but a level-l cell is amr_ref_ratio**l smaller, so its CFL limit is
        ! amr_ref_ratio**amr_max_level tighter than the coarse grid's. The dt (fixed, or the coarse-only cfl_dt estimate) is NOT
        ! scaled for that, so a coarse-CFL dt silently runs the finest block unstable. The true CFL is unknown at init, so warn
        ! rather than abort - a small enough dt is valid.
        if (proc_rank == 0 .and. (amr_ref_ratio > 2 .or. amr_max_level > 1)) then
            print '(A,I0,A)', &
                & ' [amr] WARNING: fine blocks advance at the coarse dt, but the ' &
                & // 'finest cell is amr_ref_ratio**amr_max_level = ', amr_ref_ratio**amr_max_level, &
                & 'x smaller - ensure dt satisfies the FINEST cell CFL (roughly the coarse-stable dt divided by that ' &
                & // 'factor), else the fine block may go unstable'
        end if

        ! Configuration advisories. These are advice, not constraints: every setting below is legal and sometimes correct, so
        ! they warn rather than abort.
        if (proc_rank == 0) then
            ! amr_regrid_int = 0 is static AMR: the block set never changes. That is a legitimate mode
            ! (and the only one supported above amr_max_level = 2), but a user who set `amr = T` expecting
            ! adaptivity gets none, silently.
            if (amr_regrid_int == 0) then
                print '(A)', &
                    & ' [amr] NOTE: amr_regrid_int = 0 - the block set is STATIC and never adapts. ' &
                    & // 'Set amr_regrid_int > 0 (4-8 is a reasonable start) for adaptive refinement.'
            end if
            ! The derived cap is the min-over-ranks local half-extent, so it shrinks as ranks are added
            ! (the wrong direction for strong scaling), and it makes the box set (hence the answer, within
            ! tolerance) depend on the rank count.
            if (amr_max_grid_size == 0 .and. num_procs > 1) then
                print '(A)', &
                    & ' [amr] NOTE: amr_max_grid_size = 0 derives the block cap from the ' &
                    & // 'decomposition, so it SHRINKS as ranks are added and the box set depends on rank ' &
                    & // 'count. Pinning it (64 measured best in 3D on MI250X, memory-bounded) was 3.0x ' &
                    & // 'faster and makes the box set rank-invariant.'
            end if
            ! Frequent regridding is dominated by the per-cell tag sweep, which is flat in box count.
            if (amr_regrid_int > 0 .and. amr_regrid_int < 4) then
                print '(A,I0,A)', ' [amr] NOTE: amr_regrid_int = ', amr_regrid_int, &
                    & ' regrids often; the tag sweep is per-CELL and flat in box count, so interval 8 ' &
                    & // 'measured 1.39x faster. Raise it unless the refined feature moves quickly.'
            end if
        end if

        ! Mirror decomposition: each rank holds the fine cells covering block /\ its own subdomain (np=1: the intersection is the
        ! whole block). buff_size is not available at checker time, so the geometric aborts below must live here.
        sidx = 0; ext = 0
        sidx(1) = start_idx(1); ext(1) = m
        if (n_glb > 0) then; sidx(2) = start_idx(2); ext(2) = n; end if
        if (p_glb > 0) then; sidx(3) = start_idx(3); ext(3) = p; end if
        call s_amr_compute_isect(amr_block_beg, amr_block_end)

        ! the fine ghost shell and reflux outside cells must stay inside the global domain (identical inputs on all ranks; every
        ! rank takes the same branch)
        if (amr_block_beg(1) < buff_size .or. amr_block_end(1) > m_glb - buff_size .or. (n_glb > 0 .and. (amr_block_beg(2) &
            & < buff_size .or. amr_block_end(2) > n_glb - buff_size)) .or. (p_glb > 0 .and. (amr_block_beg(3) < buff_size &
            & .or. amr_block_end(3) > p_glb - buff_size))) then
            call s_mpi_abort('amr block must lie at least buff_size cells inside the domain boundaries')
        end if

        ! Scratch constraint: the fine advance reuses the solver scratch (m_rhs/WENO/Riemann work arrays) and the global coordinate
        ! arrays, all sized to this rank's local grid. Fine-level distribution gives a block whole to its owner, so the whole
        ! block's fine extent (2*block-1) must fit every rank's local extent (a big block cannot be whole-owned; it must be split
        ! into <= local-half boxes). Checked on the replicated block box so all ranks agree. (np=1: local extent = global, so
        ! 2*block-1 <= m_glb always holds.) non-IB: the block is tiled into <= amr_maxc_fit sub-blocks (each fits every rank's
        ! scratch), so no cap is needed. IB keeps a single contiguous block per body,
        ! so an IB block must itself fit a rank's local half-extent.
        bad_loc = 0
        if (ib) then
            if (amr_ref_ratio*(amr_block_end(1) - amr_block_beg(1) + 1) - 1 > m) bad_loc = 1
            if (n_glb > 0 .and. amr_ref_ratio*(amr_block_end(2) - amr_block_beg(2) + 1) - 1 > n) bad_loc = 1
            if (p_glb > 0 .and. amr_ref_ratio*(amr_block_end(3) - amr_block_beg(3) + 1) - 1 > p) bad_loc = 1
        end if
        call s_mpi_allreduce_integer_max(bad_loc, bad_glb)
        if (bad_glb == 1) then
            call s_mpi_abort('amr fine extent exceeds a rank local grid (solver scratch is local-sized): an immersed-body block ' &
                             & // 'is owned whole and un-tiled, so it may cover at most about half of any rank subdomain per ' &
                             & // 'dimension; shrink the body region or use fewer ranks')
        end if

        ! max coarse block cells per dim (upper bound for any future regrid box); 1 for collapsed dims
        amr_maxc(1) = (m_glb + 1)/amr_ref_ratio
        amr_maxc(2) = 1; amr_maxc(3) = 1
        if (n_glb > 0) amr_maxc(2) = (n_glb + 1)/amr_ref_ratio
        if (p_glb > 0) amr_maxc(3) = (p_glb + 1)/amr_ref_ratio

        ! regrid size cap. Default (amr_max_grid_size == 0): min over ranks of the local half-extent (= amr_maxc at np=1), so any
        ! clamped box satisfies every rank's scratch constraint and can move freely across ranks. That cap shrinks as ranks are
        ! added, which tiles a fixed feature into more and more blocks the further you scale (per-block cost is roughly fixed
        ! regardless of block size, so the block count is what costs), and it makes the box set (and so the answer, within
        ! tolerance) depend on the rank count. Setting amr_max_grid_size > 0 pins the cap to an absolute number of coarse cells
        ! instead, like AMReX's max_grid_size: the box set is then identical at every rank count.
        amr_maxc_fit = amr_maxc
        do d = 1, num_dims
            call s_mpi_allreduce_integer_min((ext(d) + 1)/amr_ref_ratio, fit_d)
            if (amr_max_grid_size > 0) then
                ! Rank-independent cap, independent of fit_d. The fine advance still borrows this rank's solver scratch, but that
                ! scratch is sized to the cap rather than to the subdomain (idwbuff_alloc and m/n/p_alloc in m_global_parameters
                ! give it amr_ref_ratio*amr_max_grid_size - 1 fine cells plus the ghost shell), so a block at the cap fits however
                ! small the subdomain becomes. Per-rank scratch is then O(cap**num_dims), constant in rank count.
                amr_maxc_fit(d) = min(amr_maxc(d), amr_max_grid_size)
            else
                amr_maxc_fit(d) = min(amr_maxc(d), fit_d)
            end if
        end do

        ! preallocation cap for this rank's fine arrays: a block is owned whole, so any rank must hold an entire block. regrid
        ! clamps every box to amr_maxc_fit, so amr_maxc_fit (not the global-half amr_maxc) is the true max block a rank can own;
        ! sizing to it right-sizes the fine/coord arrays. At np=1 amr_maxc_fit == amr_maxc. When amr_max_grid_size > 0,
        ! amr_maxc_fit is not bounded by the local half-extent (the solver scratch is sized to the cap instead, see above), so a
        ! rank can own a block larger than half its own subdomain.
        maxc_loc = amr_maxc_fit

        ! max fine extents and buffered bounds for preallocation
        max_f1 = amr_ref_ratio*maxc_loc(1) - 1
        max_f2 = 0; max_f3 = 0
        if (n_glb > 0) max_f2 = amr_ref_ratio*maxc_loc(2) - 1
        if (p_glb > 0) max_f3 = amr_ref_ratio*maxc_loc(3) - 1

        amr_seam_pairs_dirty = .true.; amr_seam_pairs_nblk = -1  ! force a seam-list build on the first fine-fine halo
        amr_mesh_epoch = amr_mesh_epoch + 1
        mbuf1_lo = -buff_size; mbuf1_hi = max_f1 + buff_size
        mbuf2_lo = 0; mbuf2_hi = 0; mbuf3_lo = 0; mbuf3_hi = 0
        if (n_glb > 0) then; mbuf2_lo = -buff_size; mbuf2_hi = max_f2 + buff_size; end if
        if (p_glb > 0) then; mbuf3_lo = -buff_size; mbuf3_hi = max_f3 + buff_size; end if
        if (amr_batched_advance) then
            amr_br_batch = amr_bat_max
            ! stacked members share the batch leader's coordinate arrays in the non-stacked dimensions and read the coarse WENO
            ! coefficients at their stacked index: bit-identical to the per-block advance only where the grid spacing is bitwise
            ! uniform (every cell then carries the same dx and the same coefficients). Say so once when it is not.
            block
                integer :: nonuni, nonuni_glb
                nonuni = 0
                if (any(dx(0:m) /= dx(0))) nonuni = 1
                if (n_glb > 0) then; if (any(dy(0:n) /= dy(0))) nonuni = 1; end if
                if (p_glb > 0) then; if (any(dz(0:p) /= dz(0))) nonuni = 1; end if
                call s_mpi_allreduce_integer_max(nonuni, nonuni_glb)
                if (proc_rank == 0 .and. nonuni_glb == 1) print '(A)', &
                    & ' [amr] NOTE: amr_batched_advance on a grid whose cell ' &
                    & // 'spacing is not bitwise uniform: stacked blocks reuse the batch leader''s coordinate arrays, so the ' &
                    & // 'batched advance differs from the per-block one at roundoff'
            end block
        end if
        ! with tiles, s_l0_tiles_init's mbuf union below may still enlarge these; the scratch waits for it (see s_amr_scr_init)
        if (l0_ntile == 0) call s_amr_scr_init()

        ! Memory demand, reported rather than guessed. There is no portable way to ask how much device (or host) memory is
        ! available across four compilers and three offload backends, so no cap is derived from a memory budget. What is exactly
        ! known is the demand: a block costs 2 per-slot field families (q_cons, q_cons_stor; q_prim/rhs are pooled, one shared
        ! scratch pair, not per block) x sys_size arrays on the mbuf extents. Print it and let the reader compare against
        ! their hardware. The
        ! usable cap is set by the largest slot that fits, and slot volume goes as cap**num_dims, so one cap cannot serve 2D and
        ! 3D alike. Exceeding device memory aborts inside __tgt_target_data_begin_mapper, which presents as a hang (one rank dies,
        ! the rest block in MPI).
        if (proc_rank == 0) then
            block
                real(wp) :: slot_gib, cells, nfam
                cells = real(mbuf1_hi - mbuf1_lo + 1, wp)
                if (n_glb > 0) cells = cells*real(mbuf2_hi - mbuf2_lo + 1, wp)
                if (p_glb > 0) cells = cells*real(mbuf3_hi - mbuf3_lo + 1, wp)
                nfam = 2._wp
                slot_gib = cells*real(sys_size, wp)*nfam*real(storage_size(1._wp)/8, wp)/1024._wp**3
                print '(A,I0,A,I0,A,ES10.3,A,F8.3,A)', ' [amr] per-block slot: ', nint(cells), ' cells x sys_size x ', &
                    & nint(nfam), ' fields = ', cells*real(sys_size, wp)*nfam, ' words (', slot_gib, ' GiB per owned block)'
                print '(A,F9.2,A,I0,A)', ' [amr]   worst case if one rank owned every block: ', slot_gib*real(amr_max_blocks, &
                    & wp), ' GiB (amr_max_blocks = ', amr_max_blocks, '). Typical is amr_max_blocks/num_procs blocks per rank.'
            end block
        end if

        ! bounce buffers for copy-based coord swap (GPU-safe; same bounds as the base-level global arrays, which are sized on
        ! *_alloc - these are whole-array assigned to/from x_cb etc., so the shapes must agree)
        allocate (sw_x_cb(-1 - buff_size:m_alloc + buff_size))
        allocate (sw_x_cc(-buff_size:m_alloc + buff_size))
        allocate (sw_dx(-buff_size:m_alloc + buff_size))
        if (n_glb > 0) then
            allocate (sw_y_cb(-1 - buff_size:n_alloc + buff_size))
            allocate (sw_y_cc(-buff_size:n_alloc + buff_size))
            allocate (sw_dy(-buff_size:n_alloc + buff_size))
        end if
        if (p_glb > 0) then
            allocate (sw_z_cb(-1 - buff_size:p_alloc + buff_size))
            allocate (sw_z_cc(-buff_size:p_alloc + buff_size))
            allocate (sw_dz(-buff_size:p_alloc + buff_size))
        end if
        if (igr) then
            @:ALLOCATE(sw_jac(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end))
            @:ALLOCATE(sw_jac_old(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end))
        end if
        if (cyl_coord .and. n_glb > 0) then
            @:ALLOCATE(amr_rvw(0:max_f2))
        end if

        ! Grid uniformity policy. Both spacing-uniformity consumers are handled exactly: fine-block ghost-shell coordinates extend
        ! by exact parent-cell bisection (reads sw_*_cb), and the spacing-dependent WENO reconstruction coefficients are recomputed
        ! for the active grid on every swap/restore when the grid is nonuniform anywhere (stretched grids, or 2D-axisymmetric's
        ! half-width axis cell dy(0) = dy/2). Tolerance is epsilon-scaled: an absolute 1e-12 would sit below single-precision grid
        ! roundoff and classify every grid as stretched (spuriously tripping the stretched-combo gates). On uniform grids the flag
        ! stays false and behavior is bit-identical to the reuse path. The stretch_* flags are pre_process-only, so the grid itself
        ! is checked (this also catches externally generated grids).
        if (maxval(dx(0:m)) - minval(dx(0:m)) > 1.e3_wp*epsilon(1._wp)*maxval(dx(0:m))) then
            amr_weno_coef_recompute = .true.; amr_grid_stretched = .true.
        end if
        if (n_glb > 0) then
            ! interior nonuniformity is stretching; a lone dy(0) deviation is stretching only when it is NOT the axisymmetric
            ! half-width axis cell
            if (n > 0 .and. maxval(dy(1:n)) - minval(dy(1:n)) > 1.e3_wp*epsilon(1._wp)*maxval(dy(1:n))) then
                amr_weno_coef_recompute = .true.; amr_grid_stretched = .true.
            end if
            if (abs(dy(0) - dy(min(1, n))) > 1.e3_wp*epsilon(1._wp)*dy(min(1, n))) then
                amr_weno_coef_recompute = .true.
                if (.not. cyl_coord) amr_grid_stretched = .true.
            end if
        end if
        if (p_glb > 0) then
            if (maxval(dz(0:p)) - minval(dz(0:p)) > 1.e3_wp*epsilon(1._wp)*maxval(dz(0:p))) then
                amr_weno_coef_recompute = .true.; amr_grid_stretched = .true.
            end if
        end if
        if (weno_order == 1 .or. igr) amr_weno_coef_recompute = .false.  ! order 1 / IGR: no grid-dependent WENO coefficients
        ! lint: runtime-check. The batched slab installs only the leader's dx/dy/dz (not the cell boundaries), so a per-swap
        ! coefficient
        ! recompute would give members 2..nb coefficients from stale boundaries; the grid test above is the runtime authority
        if (amr_batched_advance .and. amr_weno_coef_recompute) call s_mpi_abort('amr_batched_advance requires a uniform grid: ' &
            & // 'the per-block WENO coefficient recompute is armed on this one')

        ! persistent global coarse boundaries: the fine-distribution owner rebuilds whole-block fine coordinates from these (needed
        ! once the fine level is decoupled from the coarse decomposition; harmless otherwise)
        call s_amr_build_global_cb()
        ! Fail closed on stretched grid + Lagrangian/IB-dynamic-regrid. Two independent blockers:
        !  (1) the position->global-cell-index conversions here use int((x-beg)/dx(0)), inexact on a stretched grid and
        !      rank-inconsistent (dx(0) is rank-local). Fixable by bisection-searching the global cell-boundary arrays.
        !  (2) IB/Lagrangian floor buff_size (10/6), but s_amr_recompute_weno_coefs (armed only on nonuniform grids) indexes
        !      poly_coef_cb* over -buff_size:m+buff_size while m_weno sized those arrays with a smaller buff_size at init, an
        !      out-of-bounds write in s_compute_weno_coefficients. The WENO coefficient arrays need sizing to the final buff_size
        !      (or the recompute clamped to the module's true bounds) before this gate can lift. Fix (1) alone is insufficient.
        if (amr_grid_stretched .and. (bubbles_lagrange .or. (ib .and. amr_regrid_int > 0))) then
            call s_mpi_abort('amr on a stretched grid does not support ' &
                             & // 'Lagrangian bubbles or dynamic regrid with immersed bodies: their ' &
                             & // 'position-to-cell-index conversions assume uniform spacing')
        end if

        ! per-slot field arrays are allocated by s_amr_alloc_slot / freed by s_amr_free_slot (sized to the max buffered block). The
        ! lazy owned-only reconcile that keeps a rank's fine memory ~1/num_procs of the pool follows. The QBMM RHS scratch
        ! (amr_rhs_pb_f/mv_f) is single - allocate once.
        allocate (amr_slot_live(amr_max_blocks)); amr_slot_live = .false.
        if (qbmm .and. .not. polytropic) then
            @:ALLOCATE(amr_rhs_pb_f(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi, 1:nnode, 1:nb))
            @:ALLOCATE(amr_rhs_mv_f(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi, 1:nnode, 1:nb))
        end if
        ! per-slot field arrays are allocated lazily by s_amr_reconcile_slots once ownership is known (after the block setup +
        ! s_amr_assign_block_owners below), so a rank holds only its owned blocks' fine arrays - not all amr_max_blocks slots.

        ! fine-level distribution: coarse-patch gather buffer (see decl). Sized to the largest block's coarse footprint (block
        ! coarse
        ! cells + 2*nmar halo, block-local frame). Device-mapped so the runtime ghost-fill reads it on the owner.
        amr_cpat_mar = (buff_size + amr_ref_ratio - 1)/amr_ref_ratio + 1
        amr_cpat_hi = 0
        amr_cpat_hi(1) = maxc_loc(1) - 1 + 2*amr_cpat_mar
        if (n_glb > 0) amr_cpat_hi(2) = maxc_loc(2) - 1 + 2*amr_cpat_mar
        if (p_glb > 0) amr_cpat_hi(3) = maxc_loc(3) - 1 + 2*amr_cpat_mar
        ! CCE OpenMP-offload leaves a bare module-scope derived-type (scalar_field) allocatable's descriptor uninitialized, so a
        ! direct allocate(amr_cg(1:sys_size)) aborts with lib-4425 at program start (a local scalar_field array and a
        ! GPU_DECLARE'd module one like q_prim_vf both allocate fine; only a bare module array does not). Allocate a local, which
        ! gets a valid descriptor, and hand it to the module variable via move_alloc, then map. OpenACC is unaffected but takes
        ! the same path correctly.
        allocate (tmp_cg(1:sys_size))
        @:ALLOCATE(amr_slab_tab(1:8, 1:6))
        call move_alloc(tmp_cg, amr_cg)
        $:GPU_ENTER_DATA(create='[amr_cg]')
        do i = 1, sys_size
            @:ALLOCATE(amr_cg(i)%sf(0:amr_cpat_hi(1), 0:amr_cpat_hi(2), 0:amr_cpat_hi(3)))
            amr_cg(i)%sf = 0._stp  ! padding beyond a block's valid patch extent is never read; keep it finite for the device copy
            @:ACC_SETUP_SFs(amr_cg(i))
        end do

        ! non-polytropic QBMM: gathered coarse pb/mv patch (analogue of amr_cg, same footprint + trailing (nnode, nb) dims). Plain
        ! 5D
        ! arrays (amr_rhs_pb_f idiom): the module GPU_DECLARE + @:ALLOCATE handle device mapping - no @:ACC_SETUP_SFs.
        if (qbmm .and. .not. polytropic) then
            @:ALLOCATE(amr_cg_pb(0:amr_cpat_hi(1), 0:amr_cpat_hi(2), 0:amr_cpat_hi(3), 1:nnode, 1:nb))
            @:ALLOCATE(amr_cg_mv(0:amr_cpat_hi(1), 0:amr_cpat_hi(2), 0:amr_cpat_hi(3), 1:nnode, 1:nb))
            amr_cg_pb = 0._stp; amr_cg_mv = 0._stp
        end if

        ! the coarse decomposition (each rank's coarse start_idx + local m/n/p) is a structured cartesian split, computed O(1) per
        ! rank by s_amr_rank_decomp - no replicated table, no allgather. Validate the formula against this rank's actual values.
        call s_amr_validate_decomp()

        ! per-slot fine-grid IB marker fields (static-body AMR); sized to the same max buffered fine extents as q_cons so the fine
        ! IB pipeline can resolve the body on the block
        if (ib) call s_ibm_alloc_fine(amr_max_blocks, mbuf1_lo, mbuf1_hi, mbuf2_lo, mbuf2_hi, mbuf3_lo, mbuf3_hi)

        ! set geometry (region, m/n/p, idwbuff, coordinates) for the initial block (amr_cur = f_l0_slot(1), the initial fine-block
        ! slot). Under dynamic regrid with bodies
        ! the initial block gets the same body-containment expansion regrid boxes get (the moving-body containment guard requires it
        ! from step 1); for a static block (amr_regrid_int = 0) the user's placement is authoritative. max_grid_size tiling: the
        ! initial block splits into <= amr_maxc_fit sub-blocks (at np=1 amr_maxc_fit == amr_maxc so a normal block stays a single
        ! tile - unchanged), one per slot; IB keeps a single contiguous block.
        blk_lo = amr_block_beg; blk_hi = amr_block_end
        if (ib .and. amr_regrid_int > 0) call s_amr_expand_box_over_bodies(blk_lo, blk_hi)
        block
            type(t_box), allocatable :: tiled(:)
            integer                  :: nt, capt, kk
            allocate (tiled(amr_max_blocks)); nt = 0; capt = 0
            if (ib) then
                nt = 1; tiled(1)%lo = blk_lo; tiled(1)%hi = blk_hi
            else
                call s_amr_tile_box(blk_lo, blk_hi, tiled, nt, amr_max_fine, capt)
            end if
            amr_num_blocks = f_l0_slot(nt)  ! fine blocks occupy [l0_slot_off+1 .. l0_slot_off+nt] in the shared pool
            ! set block regions first so the owner assignment (reads amr_region_*_all) runs before the owner-dependent geometry -
            ! else s_set_amr_fine_geometry would size the whole-block owner from a stale (default) amr_block_owner
            do kk = 1, nt
                amr_region_lo_all(:,f_l0_slot(kk)) = tiled(kk)%lo; amr_region_hi_all(:,f_l0_slot(kk)) = tiled(kk)%hi
            end do
            call s_amr_assign_block_owners()  ! assign each block's single owner rank (fine-dist map)
            call s_amr_reconcile_slots()  ! allocate this rank's owned initial blocks (owner-guarded geometry writes below)
            do kk = 1, nt
                amr_cur = f_l0_slot(kk)
                call s_set_amr_fine_geometry(tiled(kk)%lo, tiled(kk)%hi)
            end do
            call s_amr_reduce_xchg_flag()
            call s_amr_select_slot(f_l0_slot(1))  ! refresh the per-block mirrors (geometry loop left them on the last tile)
            deallocate (tiled)
        end block

        ! per-family tag bases sit above the per-box tag space so the two cannot collide
        block
            integer :: f
            do f = 1, size(amr_tag_base)
                amr_tag_base(f) = amr_max_blocks + 100*f
            end do
            ! keyed band space starts at the next 65536 boundary above every per-box tag (bases + their mod-100 folds)
            amr_m1_base = ((amr_tag_base(size(amr_tag_base)) + 100)/65536 + 1)*65536
        end block
#ifdef MFC_MPI
        block
            integer(kind=MPI_ADDRESS_KIND) :: tag_ub
            logical                        :: tag_ub_set
            integer                        :: ierr
            call MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, tag_ub, tag_ub_set, ierr)
            @:ASSERT(tag_ub_set, "MPI_TAG_UB attribute unavailable")
            @:ASSERT(amr_tag_base(size(amr_tag_base)) + 100 <= tag_ub, &
                     & "AMR tag space exceeds MPI_TAG_UB: amr_max_blocks is too large for this MPI's tag range")
            @:ASSERT(amr_m1_base + 8*65536 <= tag_ub, "AMR keyed-tag band space exceeds MPI_TAG_UB")
        end block
#endif

    end subroutine s_initialize_amr_module

    !> [amr-cad] report: SUM-allreduce the cadence counters and print once on rank 0. Collective: the caller (s_finalize_amr_module)
    !! runs it before the amr early-return so every rank participates (all-zero when amr is off).
    impure subroutine s_amr_cad_report()

        integer(8) :: cad(2), cadr(2)
        integer    :: ierr

        cad(1) = amr_cad_tot; cad(2) = amr_cad_esc; cadr = cad
#ifdef MFC_MPI
        call MPI_ALLREDUCE(cad, cadr, 2, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
        if (proc_rank == 0) then
            ! cadence containment: escaped > 0 means a feature outran amr_buf between regrids (see the decl)
            if (cadr(1) > 0) write (0, '(A,I0,A,I0,A,F6.3)') ' [amr-cad] L1 tags ', cadr(1), ' escaped ', cadr(2), ' frac ', &
                & real(cadr(2))/real(cadr(1))
        end if

    end subroutine s_amr_cad_report

    impure subroutine s_finalize_amr_module()

        integer :: i, islot

        ! before the amr early-return: the report's conservation allreduce is collective, and the L0 tile
        ! families can fire with amr = F. All ranks take the same path either way.

        call s_xa_report()
        call s_amr_cad_report()
        if (rank_time_wrt) then
            write (0, '(A,I0,A,I0,A,I0,A)', advance='no') '[amr-bat] rank ', proc_rank, ' batches ', sum(amr_bat_hist), &
                   & ' single ', amr_bat_hist(1), ' sizes'
            do i = 1, amr_bat_max; write (0, '(A,I0,A,I0)', advance='no') ' ', i, 'x', amr_bat_hist(i); end do
            write (0, '(A)') ''
        end if
        if (.not. amr) return
        do islot = 1, amr_max_blocks
            call s_amr_free_slot(islot)
        end do
        if (qbmm .and. .not. polytropic) then
            @:DEALLOCATE(amr_rhs_pb_f)
            @:DEALLOCATE(amr_rhs_mv_f)
            @:DEALLOCATE(amr_cg_pb)
            @:DEALLOCATE(amr_cg_mv)
        end if
        deallocate (amr_slot_live)
        call s_amr_st_finalize()
        if (allocated(amr_seam_pairs)) deallocate (amr_seam_pairs)
        if (allocated(amr_ovl_gather)) deallocate (amr_ovl_gather)
        if (allocated(amr_ovl_scatter)) deallocate (amr_ovl_scatter)
        deallocate (amr_ovl_gather_n, amr_ovl_scatter_n)
        if (allocated(amr_gpl_nsrc)) deallocate (amr_gpl_nsrc, amr_gpl_src, amr_gpl_sz, amr_gpl_psrc, amr_gpl_psz)
        if (allocated(amr_gcr_pool)) deallocate (amr_gcr_pool)
        if (allocated(amr_gcr_req)) deallocate (amr_gcr_req, amr_gcr_off)
        ! per-array guards, not grouped on a lead member: the wave-scratch arrays of one group allocate
        ! independently (spsz/rpsz are sized only by the qbmm pb/mv wave branch), so a non-qbmm np>1 run
        ! reaches here with a group partially allocated. gfortran/ifx abort on deallocating an unallocated
        ! array (amdflang silently tolerates it).
        #:for A in ['amr_fw_sblk', 'amr_fw_sbl', 'amr_fw_sbh', 'amr_fw_spi', 'amr_fw_sqo', 'amr_fw_spo', &
            'amr_fw_rblk', 'amr_fw_rbl', 'amr_fw_rbh', 'amr_fw_rpi', 'amr_fw_rqo', 'amr_fw_rpo', &
            'amr_fw_sprank', 'amr_fw_sqsz', 'amr_fw_spsz', 'amr_fw_snxp', 'amr_fw_sqbase', 'amr_fw_spbase', &
            'amr_fw_rprank', 'amr_fw_rqsz', 'amr_fw_rpsz', 'amr_fw_rnxp', 'amr_fw_rqbase', 'amr_fw_rpbase', &
            'amr_fw_map', 'amr_fw_nx', 'amr_fw_pq', 'amr_fw_pp']
            if (allocated(${A}$)) deallocate (${A}$)
        #:endfor
        #:for A in ['amr_fw_sq', 'amr_fw_sp', 'amr_fw_rq', 'amr_fw_rp']
            if (allocated(${A}$)) then
                if (amr_fw_dev) then
                    $:GPU_EXIT_DATA(delete='[' + A + ']')
                end if
                deallocate (${A}$)
            end if
        #:endfor
        if (allocated(amr_fw_req)) deallocate (amr_fw_req, amr_fw_reqw)
        #:for A in ['amr_my_blk', 'amr_l1r_blk', 'amr_l1p_blk', 'amr_fch_blk', 'amr_own_blk', 'amr_parent_blk', &
            'amr_child_ptr', 'amr_child_idx', 'amr_gpk']
            if (allocated(${A}$)) deallocate (${A}$)
        #:endfor
        do i = 1, sys_size
            @:DEALLOCATE(amr_cg(i)%sf)
        end do
        @:DEALLOCATE(amr_cg)
        @:DEALLOCATE(amr_slab_tab)
        deallocate (amr_slots)
        deallocate (amr_region_lo_all, amr_region_hi_all, amr_isect_lo_all, amr_isect_hi_all, amr_owns_all)
        if (allocated(sw_x_cb)) deallocate (sw_x_cb, sw_x_cc, sw_dx)
        if (allocated(sw_y_cb)) deallocate (sw_y_cb, sw_y_cc, sw_dy)
        if (allocated(sw_z_cb)) deallocate (sw_z_cb, sw_z_cc, sw_dz)
        if (allocated(amr_block_owner)) deallocate (amr_block_owner)
        if (allocated(amr_owner_cut)) deallocate (amr_owner_cut)
        if (allocated(amr_fine_cut)) deallocate (amr_fine_cut)
        if (allocated(amr_block_level)) deallocate (amr_block_level)
        if (allocated(amr_gxcb)) deallocate (amr_gxcb)
        if (allocated(amr_gycb)) deallocate (amr_gycb)
        if (allocated(amr_gzcb)) deallocate (amr_gzcb)
        if (igr) then
            @:DEALLOCATE(sw_jac)
            @:DEALLOCATE(sw_jac_old)
        end if
        if (cyl_coord .and. n_glb > 0) then
            @:DEALLOCATE(amr_rvw)
        end if

    end subroutine s_finalize_amr_module

end module m_amr
