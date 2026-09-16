!>
!!@file
!!@brief Contains module m_amr_advance

#! AMD OpenMP lane: assert allocatables present on every kernel here (see OMP_DEFAULT_STR). Every conditionally allocated
#! module array a kernel here names launches only under its allocation's own condition (amr_rvw: cyl_coord; sw_jac/jac: igr;
#! amr_cg_pb/mv: do_pbmv; amr_gst_a/b: amr_subcycle; amr_prim_st/amr_bt_*: amr_prim_batch); amr_cg and amr_cons_br/stor_st are
#! allocated before first use. A kernel naming an unallocated array aborts. Keep it so.
#:set MFC_OMP_PRESENT_ALLOCATABLE = True
#:include 'macros.fpp'

!> @brief Fine-block stage advance and the immersed-boundary/Lagrange fine services.
module m_amr_advance

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

    implicit none

    private
    public :: s_amr_fine_stage_advance, s_amr_fine_stage_advance_batched, s_amr_fine_stage_rhs, s_amr_fine_stage_rk, s_amr_setup_ib

contains

    !> Compute the fine-grid IB state (markers/ghost points/levelset) for every active block from the body geometry (static-body
    !! AMR). Called once after the coarse IB setup at init (regrid+IB is gated). Per slot with fine cells: swap the grid to the fine
    !! block, swap the IB globals to the slot store, run the fine IB pipeline (writing into the slot store), restore. No-op unless
    !! amr .and. ib.
    impure subroutine s_amr_setup_ib()

        integer         :: islot, save_cur
        integer(kind=8) :: my_ib_gps, nrank_ib

        if (.not. amr .or. .not. ib) return
        save_cur = amr_cur
        my_ib_gps = 0_8
        do islot = 1, amr_num_blocks
            call s_amr_select_slot(islot)
            if (.not. amr_rank_owns_block) cycle
            call s_amr_swap_to_fine()
            call s_ibm_swap_to_fine(islot, gps_on_device=.false.)
            call s_ibm_setup_fine()
            my_ib_gps = my_ib_gps + int(num_gps, 8)
            call s_ibm_restore_from_fine(islot)
            call s_amr_restore_coarse()
        end do
        call s_amr_select_slot(save_cur)

        ! The fine-IB image-point stencil is not decomposition-exact across a rank seam. If the body's fine ghost points appear on
        ! more than one rank (the body straddles a coarse/fine rank boundary), abort rather than return a wrong body-surface state.
        ! A
        ! body wholly within one rank is decomposition-exact.
        call s_mpi_allreduce_integer_sum(merge(1_8, 0_8, my_ib_gps > 0_8), nrank_ib)
        if (nrank_ib > 1_8) then
            call s_mpi_abort('amr with ib: the immersed body straddles a rank boundary, where the ' &
                             & // 'fine-IB image-point stencil is not yet decomposition-exact; keep the ' &
                             & // 'body within a single rank subdomain (use fewer ranks or reposition it).')
        end if

    end subroutine s_amr_setup_ib

    !> Apply the IB state correction on the current fine block after its RK update (static-body AMR). Mirrors the coarse per-stage
    !! s_ibm_correct_state: swap the grid + IB globals to the fine block, correct q_cons/q_prim at the fine body/ghost cells,
    !! restore. amr_cur / amr_rank_owns_block are set by the caller (the per-block advance loop). No-op unless ib.
    impure subroutine s_amr_ib_correct_fine(q_prim_b)

        !> the q_prim the block's RHS pass filled (pooled scratch for fine blocks; per-slot for L0 tiles, where other tiles' RHS
        !! work ran in between - ib is in the copy-out gate, so a tile slot always has its own q_prim when this reads it)
        type(scalar_field), dimension(1:sys_size), intent(inout) :: q_prim_b

        if (.not. ib) return
        if (.not. amr_rank_owns_block) return
        call s_amr_swap_to_fine()
        call s_ibm_swap_to_fine(amr_cur, gps_on_device=.true.)
        call s_amr_br_load(amr_loc_of(amr_cur))
        if (qbmm .and. .not. polytropic) then
            ! mirror the coarse correct-state: non-polytropic QBMM also corrects the block's own pb/mv side-state at the body ghost
            ! points (bounds match the swapped fine idwbuff)
            call s_ibm_correct_state(amr_cons_br, q_prim_b, amr_slots(amr_cur)%pb_f%sf, amr_slots(amr_cur)%mv_f%sf)
        else
            call s_ibm_correct_state(amr_cons_br, q_prim_b)
        end if
        call s_amr_br_store(amr_loc_of(amr_cur))
        call s_ibm_restore_from_fine(amr_cur)
        call s_amr_restore_coarse()

    end subroutine s_amr_ib_correct_fine

    !> Rebuild the current fine block's IB state (markers/ghost points/image points) from the moving body's position (prescribed
    !! motion, moving_ibm==1). Reuses the coarse s_update_mib recompute on the swapped-in fine slot (grid + IB globals swapped to
    !! the fine block, recompute writes into the slot store, restore) at the body's current position. No-op unless ib. Must precede
    !! s_amr_ib_correct_fine.
    impure subroutine s_amr_update_mib_fine()

        integer :: i, blo(3), bhi(3)
        logical :: ovl, inside

        if (.not. ib) return
        if (.not. amr_rank_owns_block) return
        ! A moving body must stay inside its block (a body overlapping the block edge gets silently clipped forcing, so abort
        ! instead).
        ! Under dynamic regrid the expansion contained it with margin max(amr_buf,4) and body + image-point stencil reach (2
        ! coarse cells) must remain contained between regrids; on a static block the user's placement is authoritative (validated
        ! configs sit tighter than the regrid margin), so only the body bbox itself must stay inside. Consecutive contained
        ! positions keep every sub-time interpolate contained (axis-aligned boxes are convex in the linearly-moving corners).
        if (any(patch_ib(1:num_ibs)%moving_ibm /= 0)) then
            do i = 1, num_ibs
                if (patch_ib(i)%moving_ibm == 0) cycle
                call s_amr_body_bbox(i, merge(2, 0, amr_regrid_int > 0), blo, bhi)
                ovl = blo(1) <= amr_slots(amr_cur)%region%hi(1) .and. bhi(1) >= amr_slots(amr_cur)%region%lo(1)
                if (n_glb > 0) ovl = ovl .and. blo(2) <= amr_slots(amr_cur)%region%hi(2) .and. bhi(2) &
                    & >= amr_slots(amr_cur)%region%lo(2)
                if (p_glb > 0) ovl = ovl .and. blo(3) <= amr_slots(amr_cur)%region%hi(3) .and. bhi(3) &
                    & >= amr_slots(amr_cur)%region%lo(3)
                if (.not. ovl) cycle
                inside = blo(1) >= amr_slots(amr_cur)%region%lo(1) .and. bhi(1) <= amr_slots(amr_cur)%region%hi(1)
                if (n_glb > 0) inside = inside .and. blo(2) >= amr_slots(amr_cur)%region%lo(2) .and. bhi(2) &
                    & <= amr_slots(amr_cur)%region%hi(2)
                if (p_glb > 0) inside = inside .and. blo(3) >= amr_slots(amr_cur)%region%lo(3) .and. bhi(3) &
                    & <= amr_slots(amr_cur)%region%hi(3)
                if (.not. inside) then
                    call s_mpi_abort('amr with moving ib: the body reached the fine-block boundary; ' &
                                     & // 'under dynamic regrid reduce amr_regrid_int or increase amr_buf, ' &
                                     & // 'for a static block enlarge it to contain the trajectory')
                end if
            end do
        end if
        call s_amr_swap_to_fine()
        call s_ibm_swap_to_fine(amr_cur, gps_on_device=.true.)
        call s_update_mib(num_ibs)
        call s_ibm_restore_from_fine(amr_cur)
        call s_amr_restore_coarse()

    end subroutine s_amr_update_mib_fine

    !> Device RK stage update over the fine interior: q = (c1*q + c2*q_stor + c3*dt_in*rhs)/c4 (compute in wp, store stp). Mirrors
    !! the coarse non-IGR rk_coef form in s_tvd_rk. Twin s_amr_fine_rk_update_pbmv + s_tvd_rk (m_time_steppers): same SSP-RK stage
    !! combination; keep all three in lockstep.
    impure subroutine s_amr_fine_rk_update(loc, q_rhs, c1, c2, c3, c4, dt_in)

        integer, intent(in)                                 :: loc  !< flat-store slot of the updated block
        type(scalar_field), dimension(sys_size), intent(in) :: q_rhs
        real(wp), intent(in)                                :: c1, c2, c3, c4, dt_in
        integer                                             :: i, fi, fj, fk, fm, fn, fp

        fm = amr_slots(amr_cur)%m; fn = amr_slots(amr_cur)%n; fp = amr_slots(amr_cur)%p
        $:GPU_PARALLEL_LOOP(collapse=4)
        do i = 1, sys_size
            do fk = 0, fp
                do fj = 0, fn
                    do fi = 0, fm
                        amr_cons_st(fi, fj, fk, i, loc) = (c1*real(amr_cons_st(fi, fj, fk, i, loc), wp) + c2*real(amr_stor_st(fi, &
                                    & fj, fk, i, loc), wp) + c3*dt_in*real(q_rhs(i)%sf(fi, fj, fk), wp))/c4
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_fine_rk_update

    !> Advance phase of a fine RK stage: fine RHS + RK update (+ QBMM/6eq/IB) for the current block. Owner-only. Reads the block's
    !! ghost shell (coarse prolong + fine-fine halo already applied by the fill + halo phases).
    !> One fine-block RK stage = RHS pass then RK pass. Fused wrapper: the AMR fine blocks (m_time_steppers) call this so their
    !! rhs+rk stay back-to-back. The coexist tile path (s_l0_advance_stage) instead calls the two passes directly with the
    !! reflux-delta copy-back interposed between them, so the corrected coarse rhs reaches the tile before its RK update.
    impure subroutine s_amr_fine_stage_advance(s, coefs, bc_type, q_T_sf, pb_in, rhs_pb, mv_in, rhs_mv, t_step)

        integer, intent(in)                                        :: s, t_step
        real(wp), intent(in)                                       :: coefs(4)
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        type(scalar_field), intent(inout)                          :: q_T_sf
        real(stp), dimension(:,:,:,:,:), intent(inout)             :: pb_in, mv_in
        real(wp), dimension(:,:,:,:,:), intent(inout)              :: rhs_pb, rhs_mv

        call s_amr_fine_stage_rhs(s, bc_type, q_T_sf, amr_scr_prim, amr_scr_rhs, pb_in, rhs_pb, mv_in, rhs_mv, t_step)
        call s_amr_fine_stage_rk(s, coefs, amr_scr_prim, amr_scr_rhs)

    end subroutine s_amr_fine_stage_advance

    !> Batched fine advance (amr_batched_advance): the owned fine blocks are advanced in batches of up to amr_bat_max blocks of the
    !! same (level, extent), stacked two ghost shells apart along the last active dimension in the bridge, so one s_compute_rhs call
    !! and one RK kernel cover the batch. The per-cell arithmetic is the per-block path's; blocks are independent (each advance
    !! writes only its own store column and register slots), so the grouping order is free. The stacked members read the batch
    !! leader's coordinate arrays in the non-stacked dimensions (see the init-time note in s_initialize_amr_module).
    impure subroutine s_amr_fine_stage_advance_batched(s, coefs, bc_type, q_T_sf, pb_in, rhs_pb, mv_in, rhs_mv, t_step)

        integer, intent(in)                                        :: s, t_step
        real(wp), intent(in)                                       :: coefs(4)
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        type(scalar_field), intent(inout)                          :: q_T_sf
        real(stp), dimension(:,:,:,:,:), intent(inout)             :: pb_in, mv_in
        real(wp), dimension(:,:,:,:,:), intent(inout)              :: rhs_pb, rhs_mv
        integer                                                    :: i, j, g, h, ibm, loc, nb
        logical, allocatable                                       :: done(:)
        logical                                                    :: last_batch
        real(wp)                                                   :: tb0, tb1, tb2, tb3, tb4
        character(len=32)                                          :: bfn

        call s_amr_refresh_my_blocks()
        allocate (done(amr_n_my)); done = .false.
        do i = 1, amr_n_my
            if (done(i)) cycle
            done(i) = .true.
            g = amr_my_blk(i)
            if (amr_block_level(g) == 0) cycle  ! L0 tile slots are advanced by s_l0_advance_stage
            amr_bat_n = 1; amr_bat_blk(1) = g
            do j = i + 1, amr_n_my
                if (amr_bat_n == amr_bat_max) exit
                h = amr_my_blk(j)
                if (done(j) .or. amr_block_level(h) /= amr_block_level(g)) cycle
                ! padded membership: no larger than the leader in any dim, and the padding wastes <= amr_bat_pad of its cells
                if (amr_slots(h)%m > amr_slots(g)%m .or. amr_slots(h)%n > amr_slots(g)%n .or. amr_slots(h)%p > amr_slots(g)%p) cycle
                if (real((amr_slots(g)%m + 1)*(amr_slots(g)%n + 1)*(amr_slots(g)%p + 1) - (amr_slots(h)%m + 1)*(amr_slots(h)%n &
                    & + 1)*(amr_slots(h)%p + 1), &
                    & wp) > amr_bat_pad*real((amr_slots(h)%m + 1)*(amr_slots(h)%n + 1)*(amr_slots(h)%p + 1), wp)) cycle
                amr_bat_n = amr_bat_n + 1; amr_bat_blk(amr_bat_n) = h; done(j) = .true.
            end do
            amr_bat_hist(amr_bat_n) = amr_bat_hist(amr_bat_n) + 1
            ! no fine block left undone -> this batch's restore must push the coarse grid state (the coarse stage reads it)
            last_batch = .not. any(.not. done .and. amr_block_level(amr_my_blk(1:amr_n_my)) /= 0)
            ! the batch frame: leader selected (swap, capture and RK read amr_cur / the slot's extents), members' store columns
            call s_amr_select_slot(g)
            amr_bat_ext = [amr_slots(g)%m, amr_slots(g)%n, amr_slots(g)%p]
            amr_bat_sd = num_dims
            amr_bat_w = amr_bat_ext(amr_bat_sd) + 2*buff_size + 1
            do ibm = 1, amr_bat_n
                amr_bat_loc(ibm) = amr_loc_of(amr_bat_blk(ibm))
                amr_bat_mext(:,ibm) = [amr_slots(amr_bat_blk(ibm))%m, amr_slots(amr_bat_blk(ibm))%n, amr_slots(amr_bat_blk(ibm))%p]
            end do
            $:GPU_UPDATE(device='[amr_bat_loc, amr_bat_mext]')
            if (rank_time_wrt) call s_rank_time_tic()
            ! step-entry backup for the SSP-RK combination, per member (device copy over the member's buffered extents)
            if (s == 1) then
                do ibm = 1, amr_bat_n
                    h = amr_bat_blk(ibm); loc = amr_loc_of(h)
                    call s_amr_copy_fine_fields(loc, amr_slots(h)%idwbuff(1)%beg, amr_slots(h)%idwbuff(1)%end, &
                                                & amr_slots(h)%idwbuff(2)%beg, amr_slots(h)%idwbuff(2)%end, &
                                                & amr_slots(h)%idwbuff(3)%beg, amr_slots(h)%idwbuff(3)%end)
                end do
            end if
            amr_in_fine_advance = .true.
            tb0 = f_amr_wtime()
            call s_phase_tic(PH_SWAP)
            call s_amr_swap_to_fine()  ! the leader's grid, extended into the slab (amr_bat_n > 1)
            idwint = idwbuff  ! widen the conversion range to the ghost shells (restored by s_amr_restore_coarse)
            $:GPU_UPDATE(device='[idwint]')
            call s_phase_toc(PH_SWAP)
            tb1 = f_amr_wtime()
            call s_phase_tic(PH_RHS)
            call s_amr_br_load_batch(amr_bat_n)
            ! each member's own fine markers at its slab offset, for the RHS body-cell zeroing
            if (ib) call s_ibm_load_fine_markers(amr_bat_n, amr_bat_blk(1:amr_bat_n), amr_bat_mext(:,1:amr_bat_n), amr_bat_sd, &
                & amr_bat_w)
            call s_compute_rhs(amr_cons_br, q_T_sf, amr_scr_prim, bc_type, amr_scr_rhs, pb_in, rhs_pb, mv_in, rhs_mv, t_step, s)
            call s_phase_toc(PH_RHS)
            tb2 = f_amr_wtime()
            call s_phase_tic(PH_SWAP)
            call s_amr_restore_coarse(sync_device=last_batch)
            call s_phase_toc(PH_SWAP)
            amr_in_fine_advance = .false.
            tb3 = f_amr_wtime()
            call s_phase_tic(PH_RK)
            ! IGR folds dt into its RHS, so the update multiplies by 1 there (as the per-block advance does)
            call s_amr_fine_rk_update_batch(amr_bat_n, amr_scr_rhs, coefs(1), coefs(2), coefs(3), coefs(4), merge(1._wp, dt, igr))
            if (ib .or. (model_eqns == model_eqns_6eq .and. (.not. relax))) then
                ! the per-block path runs its post-update hooks right after each block's RK update (s_amr_fine_stage_rk):
                ! the 6-equation pressure relaxation, the moving-body rebuild, the body/ghost-cell correction. Here they run once
                ! per member after the batch's update, in the member's own frame and in the per-block order; each reads only the
                ! member's own cells, so the order across members does not matter.
                ! amr_bat_n = 1 while the members are visited: s_amr_swap_to_fine extends the installed grid into the slab
                ! whenever amr_bat_n > 1, and the hooks must see the member's extents (ib_markers is sized to a block).
                nb = amr_bat_n; amr_bat_n = 1
                do ibm = 1, nb
                    call s_amr_select_slot(amr_bat_blk(ibm))
                    if (model_eqns == model_eqns_6eq .and. (.not. relax)) call s_amr_pressure_relax_fine()
                    if (ib) then
                        if (moving_immersed_boundary_flag) call s_amr_update_mib_fine()
                        call s_amr_bat_member_prim(ibm, amr_scr_prim, amr_scr_prim_blk)
                        call s_amr_ib_correct_fine(amr_scr_prim_blk)
                    end if
                end do
                amr_bat_n = nb
            end if
            call s_phase_toc(PH_RK)
            tb4 = f_amr_wtime()
            if (rank_time_wrt) then
                call s_rank_time_toc()
                if (.not. amr_bat_open) then
                    amr_bat_open = .true.
                    write (bfn, '(A,I0,A)') 'amr_batch_r', proc_rank, '.log'
                    open (newunit=amr_bat_unit, file=trim(bfn), status='replace', action='write')
                    write (amr_bat_unit, &
                           & '(A)') &
                           & '# step stage n level m n p cells_per_member t_swap t_rhs t_restore t_rk then blk:key per member'
                end if
                write (amr_bat_unit, '(I0,1X,I0,1X,I0,1X,I0,3(1X,I0),1X,I0,4(1X,ES12.5))', advance='no') t_step, s, amr_bat_n, &
                       & amr_block_level(g), amr_bat_ext(1), amr_bat_ext(2), amr_bat_ext(3), &
                       & (amr_bat_ext(1) + 1)*(amr_bat_ext(2) + 1)*(amr_bat_ext(3) + 1), tb1 - tb0, tb2 - tb1, tb3 - tb2, tb4 - tb3
                do ibm = 1, amr_bat_n
                    h = amr_bat_blk(ibm)
                    write (amr_bat_unit, '(1X,I0,A,I0)', advance='no') h, ':', f_morton(amr_region_lo_all(1, h), &
                           & amr_region_lo_all(2, h), amr_region_lo_all(3, h))
                end do
                write (amr_bat_unit, '(A)') ''
            end if
        end do
        amr_bat_n = 0
        deallocate (done)

    end subroutine s_amr_fine_stage_advance_batched

    !> RHS pass of a fine-block RK stage: step-entry backup, swap grid globals to the block, s_compute_rhs (fills amr_slots%rhs +
    !! captures the block's freg / its children's creg), restore coarse globals. Leaves the per-slot rhs ready for the RK pass (or,
    !! under coexist, for the reflux-delta copy-back before the RK pass).
    impure subroutine s_amr_fine_stage_rhs(s, bc_type, q_T_sf, q_prim_b, rhs_b, pb_in, rhs_pb, mv_in, rhs_mv, t_step)

        integer, intent(in)                                        :: s, t_step
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        type(scalar_field), intent(inout)                          :: q_T_sf
        !> the block's q_prim/rhs target: the pooled scratch for fine blocks, the slot's own arrays for L0 tiles (whose rhs must
        !! survive the whole-set RHS pass; the caller chooses - see s_l0_advance_stage_rhs)
        type(scalar_field), dimension(1:sys_size), intent(inout) :: q_prim_b, rhs_b
        real(stp), dimension(:,:,:,:,:), intent(inout)           :: pb_in, mv_in
        real(wp), dimension(:,:,:,:,:), intent(inout)            :: rhs_pb, rhs_mv

        if (.not. amr .and. l0_ntile == 0) return
        if (.not. amr_rank_owns_block) return
        if (rank_time_wrt) call s_rank_time_tic()

        ! step-entry backup for the SSP-RK combination (device copy over the current buffered extents)
        if (s == 1) then
            call s_amr_copy_fine_fields(amr_loc_of(amr_cur), amr_slots(amr_cur)%idwbuff(1)%beg, &
                                        & amr_slots(amr_cur)%idwbuff(1)%end, amr_slots(amr_cur)%idwbuff(2)%beg, &
                                        & amr_slots(amr_cur)%idwbuff(2)%end, amr_slots(amr_cur)%idwbuff(3)%beg, &
                                        & amr_slots(amr_cur)%idwbuff(3)%end)
            if (qbmm .and. .not. polytropic) call s_amr_backup_pbmv(amr_slots(amr_cur)%pb_f%sf, amr_slots(amr_cur)%mv_f%sf, &
                & amr_slots(amr_cur)%pb_stor%sf, amr_slots(amr_cur)%mv_stor%sf)
        end if

        amr_in_fine_advance = .true.
        call s_phase_tic(PH_SWAP)
        call s_amr_swap_to_fine()
        idwint = amr_slots(amr_cur)%idwbuff  ! widen the conversion range to the ghost shell (restored by s_amr_restore_coarse)
        $:GPU_UPDATE(device='[idwint]')
        call s_phase_toc(PH_SWAP)
        call s_phase_tic(PH_RHS)
        call s_amr_br_load(amr_loc_of(amr_cur))
        ! the block's own fine markers, for the RHS body-cell zeroing (the grid globals are the block's here)
        if (ib) call s_ibm_load_fine_markers(1, [amr_cur], reshape([m, n, p], [3, 1]), 1, 0)
        ! batched conversion: this block's computed prim vars (mom, E) were already produced by the stage-top batched
        ! conversion; land them and let s_compute_rhs skip its per-block conversion. L0 tile slots (level 0) are not in
        ! the batch and keep the per-block conversion.
        if (amr_prim_batch .and. amr_block_level(amr_cur) >= 1) then
            call s_amr_prim_load(q_prim_qp%vf, amr_loc_of(amr_cur))
            amr_prim_preloaded = .true.
        end if
        if (qbmm .and. .not. polytropic) then
            ! the block's own side-state and rhs scratch: the coarse pb_in/rhs_pb must not be touched at fine indices (the coarse
            ! stage consumes them after this fine stage)
            call s_compute_rhs(amr_cons_br, q_T_sf, q_prim_b, bc_type, rhs_b, amr_slots(amr_cur)%pb_f%sf, amr_rhs_pb_f, &
                               & amr_slots(amr_cur)%mv_f%sf, amr_rhs_mv_f, t_step, s)
        else
            call s_compute_rhs(amr_cons_br, q_T_sf, q_prim_b, bc_type, rhs_b, pb_in, rhs_pb, mv_in, rhs_mv, t_step, s)
        end if
        amr_prim_preloaded = .false.
        call s_amr_br_store(amr_loc_of(amr_cur))
        call s_phase_toc(PH_RHS)
        call s_phase_tic(PH_SWAP)  ! the other half of the swap pair; keep the bracket symmetric
        call s_amr_restore_coarse()
        call s_phase_toc(PH_SWAP)
        amr_in_fine_advance = .false.

    end subroutine s_amr_fine_stage_rhs

    !> RK pass of a fine-block RK stage: SSP-RK combination consuming the per-slot rhs (already reflux-corrected under coexist),
    !! then per-stage pressure relaxation / moving-IB / IB-state correction. Uses slot bounds (no grid swap needed).
    impure subroutine s_amr_fine_stage_rk(s, coefs, q_prim_b, rhs_b)

        integer, intent(in)  :: s
        real(wp), intent(in) :: coefs(4)
        !> the same q_prim/rhs pair the RHS pass of this stage filled (pooled scratch for fine blocks, per-slot for L0 tiles)
        type(scalar_field), dimension(1:sys_size), intent(inout) :: q_prim_b, rhs_b

        if (.not. amr .and. l0_ntile == 0) return
        if (.not. amr_rank_owns_block) return

        call s_phase_tic(PH_RK)
        ! RK stage update (device kernel; mirror of the coarse form. Under IGR the rhs already embeds dt, matching the coarse igr
        ! update, so the dt factor is 1)
        call s_amr_fine_rk_update(amr_loc_of(amr_cur), rhs_b, coefs(1), coefs(2), coefs(3), coefs(4), merge(1._wp, dt, igr))
        if (qbmm .and. .not. polytropic) call s_amr_fine_rk_update_pbmv(amr_slots(amr_cur)%pb_f%sf, amr_slots(amr_cur)%mv_f%sf, &
            & amr_slots(amr_cur)%pb_stor%sf, amr_slots(amr_cur)%mv_stor%sf, amr_rhs_pb_f, amr_rhs_mv_f, coefs(1), coefs(2), &
            & coefs(3), coefs(4), dt)
        ! 6-equation model: per-stage pressure relaxation on the block (before IB correct, coarse order)
        if (model_eqns == model_eqns_6eq .and. (.not. relax)) call s_amr_pressure_relax_fine()
        ! moving body: rebuild the fine-block IB state at the current (lockstep-stage) body position before the correct-state
        if (moving_immersed_boundary_flag) call s_amr_update_mib_fine()
        ! IB state correction on the fine block (mirrors the coarse per-stage correct-state; no-op unless ib)
        call s_amr_ib_correct_fine(q_prim_b)
        call s_phase_toc(PH_RK)
        if (rank_time_wrt) call s_rank_time_toc()

    end subroutine s_amr_fine_stage_rk

    !> Batched twin of s_amr_fine_rk_update: the same per-cell combination for every member of the batch, reading its rhs at the
    !! member's slab offset. One kernel per batch instead of one per block.
    impure subroutine s_amr_fine_rk_update_batch(nb, q_rhs, c1, c2, c3, c4, dt_in)

        integer, intent(in)                                 :: nb
        type(scalar_field), dimension(sys_size), intent(in) :: q_rhs
        real(wp), intent(in)                                :: c1, c2, c3, c4, dt_in
        integer                                             :: ibm, i, fi, fj, fk, fm, fn, fp, loc, o1, o2, o3

        fm = amr_bat_ext(1); fn = amr_bat_ext(2); fp = amr_bat_ext(3)
        o1 = 0; o2 = 0; o3 = 0
        select case (amr_bat_sd)
        case (1); o1 = amr_bat_w
        case (2); o2 = amr_bat_w
        case (3); o3 = amr_bat_w
        end select
        $:GPU_PARALLEL_LOOP(collapse=5, private='[loc]', copyin='[nb, fm, fn, fp, o1, o2, o3]')
        do ibm = 1, nb
            do i = 1, sys_size
                do fk = 0, fp
                    do fj = 0, fn
                        do fi = 0, fm
                            ! padding
                            if (fi > amr_bat_mext(1, ibm) .or. fj > amr_bat_mext(2, ibm) .or. fk > amr_bat_mext(3, ibm)) cycle
                            loc = amr_bat_loc(ibm)
                            amr_cons_st(fi, fj, fk, i, loc) = (c1*real(amr_cons_st(fi, fj, fk, i, loc), &
                                        & wp) + c2*real(amr_stor_st(fi, fj, fk, i, loc), &
                                        & wp) + c3*dt_in*real(q_rhs(i)%sf(fi + (ibm - 1)*o1, fj + (ibm - 1)*o2, &
                                        & fk + (ibm - 1)*o3), wp))/c4
                        end do
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_fine_rk_update_batch

end module m_amr_advance
