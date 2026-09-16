!>
!!@file
!!@brief Contains module m_amr_advance

#! AMD OpenMP lane: assert allocatables present on every kernel here (see OMP_DEFAULT_STR). Every conditionally allocated
#! module array a kernel here names launches only under its allocation's own condition (amr_rvw: cyl_coord; sw_jac/jac: igr;
#! amr_cg_pb/mv: do_pbmv; amr_gst_a/b: amr_subcycle; amr_prim_st/amr_bt_*: amr_prim_batch); amr_cg and amr_cons_br/stor_st are
#! allocated before first use. A kernel naming an unallocated array aborts. Keep it so.
#:set MFC_OMP_PRESENT_ALLOCATABLE = True
#:include 'macros.fpp'

!> @brief Fine-block stage advance, subcycling and the immersed-boundary/Lagrange fine services.
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
    public :: s_amr_advance_fine_subcycle_all, s_amr_fine_stage_advance, s_amr_fine_stage_advance_batched, s_amr_fine_stage_rhs, &
        & s_amr_fine_stage_rk, s_amr_setup_ib

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
    !! the fine block, recompute writes into the slot store, restore). For the subcycled advance pass th in [0,1], the fine
    !! substep's fraction of the coarse step: s_update_mib snapshots the body to the linear time interpolation between the coarse
    !! t^n and t^{n+1} positions, the same interpolation the subcycle applies to the fluid ghost shell. Pass th < 0 for the
    !! non-subcycled lockstep stage (uses the body's current position). No-op unless ib. Must precede s_amr_ib_correct_fine.
    impure subroutine s_amr_update_mib_fine(th)

        real(wp), intent(in) :: th
        integer              :: i, blo(3), bhi(3)
        logical              :: ovl, inside

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
        call s_update_mib(num_ibs, th)
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
                        if (moving_immersed_boundary_flag) call s_amr_update_mib_fine(-1._wp)
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
        if (moving_immersed_boundary_flag) call s_amr_update_mib_fine(-1._wp)
        ! IB state correction on the fine block (mirrors the coarse per-stage correct-state; no-op unless ib)
        call s_amr_ib_correct_fine(q_prim_b)
        call s_phase_toc(PH_RK)
        if (rank_time_wrt) call s_rank_time_toc()

    end subroutine s_amr_fine_stage_rk

    !> Per-block setup for the transposed subcycle advance (amr_subcycle): gather+prolong the selected block's two time-lerp ghost
    !! sources (parent t^n in amr_gst_a, t^{n+1} in amr_gst_b) and zero its flux registers. The collective gathers run on all ranks;
    !! the owner-only fills and register-zero are guarded. Called once per level-1 block before the transposed stage loop (which
    !! reuses the prepared ghost sources every substep).
    impure subroutine s_amr_subcycle_setup_block(q_old, q_new, pb_old, mv_old, pb_in, mv_in)

        type(scalar_field), dimension(sys_size), intent(inout)                                  :: q_old, q_new
        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(in) :: pb_old, mv_old
        real(stp), dimension(:,:,:,:,:), intent(inout)                                          :: pb_in, mv_in

        ! the two lerp sources' coarse cons halos are loop-invariant over the setup loop (both read, neither written), so the
        ! exchanges are hoisted to s_amr_advance_fine_subcycle_all

        call s_amr_check_lag_clear()  ! every rank: non-owner bubbles can reach the block across a seam

        ! fine-level distribution: gather each lerp source's coarse patch (collective, all ranks) then prolong its ghost shell on
        ! the owner. Interleaved so the single amr_cg buffer is consumed by the fill before the next gather overwrites it.
        ! non-polytropic QBMM: the pb/mv ghost shell gets the same two-source time-lerp treatment, gathered + filled interleaved
        ! with
        ! q_cons so the single amr_cg_pb/mv buffer is consumed by each fill before the next gather overwrites it. Gathers collective
        ! (all ranks, P2P); fills owner-only.
        call s_amr_gather_coarse_patch(q_old, .true.)
        call s_amr_gather_send_flush()  ! this site has blocking semantics
        if (qbmm .and. .not. polytropic) call s_amr_gather_coarse_patch_pbmv(pb_old, mv_old, .true.)
        if (amr_rank_owns_block) call s_amr_fill_fine_ghosts_gsta(amr_cg, amr_loc_of(amr_cur))
        if (amr_rank_owns_block .and. qbmm .and. .not. polytropic) call s_amr_fill_fine_ghosts_pbmv(amr_cg_pb, amr_cg_mv, &
            & amr_slots(amr_cur)%pb_ghost_a%sf, amr_slots(amr_cur)%mv_ghost_a%sf)
        call s_amr_gather_coarse_patch(q_new, .true.)
        call s_amr_gather_send_flush()  ! this site has blocking semantics
        if (qbmm .and. .not. polytropic) call s_amr_gather_coarse_patch_pbmv(pb_in, mv_in, .true.)
        if (amr_rank_owns_block) call s_amr_fill_fine_ghosts_gstb(amr_cg, amr_loc_of(amr_cur))
        if (amr_rank_owns_block .and. qbmm .and. .not. polytropic) call s_amr_fill_fine_ghosts_pbmv(amr_cg_pb, amr_cg_mv, &
            & amr_slots(amr_cur)%pb_ghost_b%sf, amr_slots(amr_cur)%mv_ghost_b%sf)
        if (.not. amr_rank_owns_block) return

        ! registers accumulate over all six stages of the transposed loop, so zero them once at setup (the stage-1 overwrite trick
        ! cannot span two substeps)
        call s_amr_zero_fine_registers()

    end subroutine s_amr_subcycle_setup_block

    !> Subcycled fine advance (amr_subcycle) over all level-1 blocks, transposed: instead of each block running its full 2x3-stage
    !! subcycle in turn, every same-level block advances stage-by-stage in lockstep with the block-to-block fine-fine seam halo
    !! (s_amr_fine_fine_halo) interposed between the ghost lerp and the RHS at each stage. That makes max_grid_size-tiled adjacent
    !! sub-blocks (which appear at np>1 when a feature exceeds a rank's slot) compute a matching shared-face flux, so the subcycle
    !! conserves at the seam. Two dt/2 SSP-RK3 substeps after the coarse step: q_old/q_new are the coarse t^n / t^{n+1} states; each
    !! stage's ghosts are the linear time interpolation at stage time theta = (substep-1 + c_s)/2 with SSP-RK3 abscissae c = [0, 1,
    !! 1/2]. Level-1 blocks drive their level-2 children per substep (s_amr_advance_children), which applies this same transposed
    !! shape at every deeper level, so L2-L2 seams are reconciled by the level-filtered halo too. The halo is a no-op with < 2
    !! adjacent same-level blocks.
    impure subroutine s_amr_advance_fine_subcycle_all(q_old, q_new, coefs, bc_type, q_T_sf, pb_old, mv_old, pb_in, rhs_pb, mv_in, &
        & rhs_mv, t_step)

        type(scalar_field), dimension(sys_size), intent(inout)                                  :: q_old, q_new
        real(wp), dimension(:,:), intent(in)                                                    :: coefs  !< rk_coef(1:3, 1:4)
        type(integer_field), dimension(1:num_dims,1:2), intent(in)                              :: bc_type
        type(scalar_field), intent(inout)                                                       :: q_T_sf
        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(in) :: pb_old, mv_old
        real(stp), dimension(:,:,:,:,:), intent(inout)                                          :: pb_in, mv_in
        real(wp), dimension(:,:,:,:,:), intent(inout)                                           :: rhs_pb, rhs_mv
        integer, intent(in)                                                                     :: t_step
        real(wp), parameter                                                                     :: c_abs(3) = [0._wp, 1._wp, 0.5_wp]
        integer                                                                                 :: islot, sub, s
        real(wp)                                                                                :: th

        if (.not. amr) return

        ! valid coarse cons ghosts on both lerp sources, once for the whole setup loop (all ranks call: pairwise halo). Neither
        ! source is written below, so this is the loop-invariant hoist described in s_amr_subcycle_setup_block. Same phase ids
        ! as the lock-step path so the two budgets read side by side.
        call s_phase_tic(PH_HALO)
        if (amr_xchg_coarse_ghosts) then
            call s_amr_exchange_coarse_cons_halo(q_old)
            call s_amr_exchange_coarse_cons_halo(q_new)
        end if
        call s_phase_toc(PH_HALO)

        ! setup: each level-1 block prepares its two time-lerp ghost sources and zeros its registers (collective; all ranks call)
        call s_phase_tic(PH_GATHER)
        do islot = 1, amr_num_blocks
            if (amr_block_level(islot) /= 1) cycle
            call s_amr_select_slot(islot)
            call s_amr_subcycle_setup_block(q_old, q_new, pb_old, mv_old, pb_in, mv_in)
        end do
        call s_phase_toc(PH_GATHER)

        do sub = 1, 2
            do s = 1, 3
                th = (real(sub - 1, wp) + c_abs(s))*0.5_wp
                ! lerp every block's ghost shell to the stage time (+ substep-entry backup) before the seam halo reads interiors
                call s_phase_tic(PH_GATHER)
                do islot = 1, amr_num_blocks
                    if (amr_block_level(islot) /= 1) cycle
                    call s_amr_select_slot(islot)
                    if (.not. amr_rank_owns_block) cycle
                    call s_amr_subtree_stage_lerp(s, th)
                end do
                call s_phase_toc(PH_GATHER)
                ! reconcile shared seam ghosts among adjacent same-level blocks so both sides compute a matching flux. Tiling can
                ! split a wide feature into adjacent sub-blocks at any rank count (amr_maxc_fit caps a box at half the global
                ! extent even at np=1), so the halo runs unconditionally; it self-no-ops when there are no seam pairs, leaving
                ! every untiled case unaffected.
                call s_phase_tic(PH_SEAM)
                call s_amr_fine_fine_halo(0)
                call s_phase_toc(PH_SEAM)
                ! RHS + RK update every block from the reconciled ghost shell
                do islot = 1, amr_num_blocks
                    if (amr_block_level(islot) /= 1) cycle
                    call s_amr_select_slot(islot)
                    if (.not. amr_rank_owns_block) cycle
                    call s_amr_subtree_stage_advance(amr_dt_fine, coefs, bc_type, q_T_sf, pb_in, rhs_pb, mv_in, rhs_mv, t_step, &
                                                     & s, th)
                end do
            end do
            ! after this substep every level-1 block is at t_b (q_cons) with t_a in q_cons_stor: level 2 subcycles within [t_a,
            ! t_b] then folds back (restrict + Berger-Colella reflux). One level-wide call, not one per parent: the level-2 seam
            ! halo inside it spans all parents, so every owner must arrive at it together. No-op for single-level.
            if (amr_max_level >= 2) call s_amr_advance_children(1, amr_dt_fine, coefs, bc_type, q_T_sf, pb_in, rhs_pb, mv_in, &
                & rhs_mv, t_step)
        end do
        call s_amr_select_slot(1)

    end subroutine s_amr_advance_fine_subcycle_all

    !> Ghost-lerp half of one subcycled fine substage for the selected block (amr_cur): time-interpolate the ghost shell to stage
    !! time th and, on substep stage 1, back up the substep-entry state. Split from the RHS half so same-level blocks can run this
    !! together and the block-to-block fine-fine seam halo can be interposed before any block reads a neighbour's interior.
    !! Owner-only (the caller guards); no numerical coupling between blocks here.
    impure subroutine s_amr_subtree_stage_lerp(s, th)

        integer, intent(in)  :: s
        real(wp), intent(in) :: th

        if (rank_time_wrt) call s_rank_time_tic()
        ! lerp the ghost shell into q_cons at the stage time (device kernel; interior untouched)
        call s_amr_lerp_fine_ghosts(amr_loc_of(amr_cur), th)
        if (qbmm .and. .not. polytropic) call s_amr_lerp_fine_ghosts_pbmv(amr_slots(amr_cur)%pb_f%sf, amr_slots(amr_cur)%mv_f%sf, &
            & amr_slots(amr_cur)%pb_ghost_a%sf, amr_slots(amr_cur)%mv_ghost_a%sf, amr_slots(amr_cur)%pb_ghost_b%sf, &
            & amr_slots(amr_cur)%mv_ghost_b%sf, th)

        ! substep-entry backup for the SSP-RK combination (device copy, interior only)
        if (s == 1) then
            call s_amr_copy_fine_fields(amr_loc_of(amr_cur), 0, amr_slots(amr_cur)%m, 0, amr_slots(amr_cur)%n, 0, &
                                        & amr_slots(amr_cur)%p)
            if (qbmm .and. .not. polytropic) call s_amr_backup_pbmv(amr_slots(amr_cur)%pb_f%sf, amr_slots(amr_cur)%mv_f%sf, &
                & amr_slots(amr_cur)%pb_stor%sf, amr_slots(amr_cur)%mv_stor%sf)
        end if
        if (rank_time_wrt) call s_rank_time_toc()

    end subroutine s_amr_subtree_stage_lerp

    !> RHS + RK-update half of one subcycled fine substage for the selected block (amr_cur): compute the fine RHS from the (already
    !! halo-reconciled) ghost shell and apply the SSP-RK stage update at the fine substep dt_sub, plus per-stage pressure relaxation
    !! and IB correction. Split from the lerp half so the fine-fine seam halo runs between them. Owner-only (caller guards).
    impure subroutine s_amr_subtree_stage_advance(dt_sub, coefs, bc_type, q_T_sf, pb_in, rhs_pb, mv_in, rhs_mv, t_step, s, th)

        real(wp), intent(in) :: dt_sub                 !< this block's substep dt (parent step / amr_ref_ratio)
        real(wp), dimension(:,:), intent(in) :: coefs  !< rk_coef(1:3, 1:4)
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        type(scalar_field), intent(inout) :: q_T_sf
        real(stp), dimension(:,:,:,:,:), intent(inout) :: pb_in, mv_in
        real(wp), dimension(:,:,:,:,:), intent(inout) :: rhs_pb, rhs_mv
        integer, intent(in) :: t_step, s
        real(wp), intent(in) :: th

        if (rank_time_wrt) call s_rank_time_tic()
        amr_in_fine_advance = .true.
        call s_amr_swap_to_fine()
        ! widen the conversion range to the ghost shell (restored by s_amr_restore_coarse)
        idwint = amr_slots(amr_cur)%idwbuff
        $:GPU_UPDATE(device='[idwint]')
        call s_amr_br_load(amr_loc_of(amr_cur))
        if (ib) call s_ibm_load_fine_markers(1, [amr_cur], reshape([m, n, p], [3, 1]), 1, 0)
        call s_phase_tic(PH_RHS)
        if (qbmm .and. .not. polytropic) then
            ! the block's own side-state and rhs scratch (the coarse arrays stay untouched)
            call s_compute_rhs(amr_cons_br, q_T_sf, amr_scr_prim, bc_type, amr_scr_rhs, amr_slots(amr_cur)%pb_f%sf, amr_rhs_pb_f, &
                               & amr_slots(amr_cur)%mv_f%sf, amr_rhs_mv_f, t_step, s)
        else
            call s_compute_rhs(amr_cons_br, q_T_sf, amr_scr_prim, bc_type, amr_scr_rhs, pb_in, rhs_pb, mv_in, rhs_mv, t_step, s)
        end if
        call s_phase_toc(PH_RHS)
        call s_amr_br_store(amr_loc_of(amr_cur))
        call s_amr_restore_coarse()
        amr_in_fine_advance = .false.
        call s_phase_tic(PH_RK)

        ! RK stage update at the fine time step (device kernel)
        call s_amr_fine_rk_update(amr_loc_of(amr_cur), amr_scr_rhs, coefs(s, 1), coefs(s, 2), coefs(s, 3), coefs(s, 4), dt_sub)
        if (qbmm .and. .not. polytropic) then
            call s_amr_fine_rk_update_pbmv(amr_slots(amr_cur)%pb_f%sf, amr_slots(amr_cur)%mv_f%sf, amr_slots(amr_cur)%pb_stor%sf, &
                                           & amr_slots(amr_cur)%mv_stor%sf, amr_rhs_pb_f, amr_rhs_mv_f, coefs(s, 1), coefs(s, 2), &
                                           & coefs(s, 3), coefs(s, 4), dt_sub)
        end if
        ! 6-equation model: per-substage pressure relaxation (instantaneous equilibration; per stage at fine dt is the same
        ! infinite-rate limit the coarse applies per stage)
        if (model_eqns == model_eqns_6eq .and. (.not. relax)) call s_amr_pressure_relax_fine()
        ! moving body: rebuild the fine-block IB state at the body's fine sub-time position (th matches the fluid-ghost lerp)
        if (moving_immersed_boundary_flag) call s_amr_update_mib_fine(th)
        ! IB state correction on the fine block after each substep RK update (no-op unless ib)
        call s_amr_ib_correct_fine(amr_scr_prim)
        call s_phase_toc(PH_RK)
        if (rank_time_wrt) call s_rank_time_toc()

    end subroutine s_amr_subtree_stage_advance

    !> Recursively subcycle every block at level plev+1, across all parents at once, within one of the parents' substeps [t_a, t_b]
    !! (duration dt_sub). Every level-plev block has just finished that substep: q_cons = parent @ t_b, q_cons_stor = parent
    !! @ t_a. Per child: gather its two ghost-lerp sources from its own parent's two snapshots (parent-fine frame), recurse into
    !! level plev+2 at dt_sub/2 (a child takes amr_ref_ratio substeps covering [t_a, t_b]), then fold back into its parent:
    !! restrict the covered cells and apply the Berger-Colella C/F flux correction (s_amr_reflux_to_parent over dt_sub, consuming
    !! the child's freg + the parent-side creg captured during this substep). The registers carry the matching per-substep time
    !! weights (freg 1/r*rk3_w, creg rk3_w), so conservation closes with no register changes.
    !!
    !! Level-wide, not per-parent. Driving one parent's whole subtree to completion before the next parent's would put the
    !! interposed s_amr_fine_fine_halo(clev) out of lockstep: the halo exchanges every level-clev seam pair, so a pair whose two
    !! blocks sit under different parents would have only one side present (co-location hides that, since both ends of such a
    !! pair then land on one rank and the exchange is a local device copy). Walking the whole level together puts every owner
    !! at the same halo. At np=1 this only re-orders independent per-parent work (each child reads solely its own parent's
    !! finished snapshots and its same-level neighbours).
    recursive subroutine s_amr_advance_children(plev, dt_sub, coefs, bc_type, q_T_sf, pb_in, rhs_pb, mv_in, rhs_mv, t_step)

        integer, intent(in)                                        :: plev
        real(wp), intent(in)                                       :: dt_sub
        real(wp), dimension(:,:), intent(in)                       :: coefs
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        type(scalar_field), intent(inout)                          :: q_T_sf
        real(stp), dimension(:,:,:,:,:), intent(inout)             :: pb_in, mv_in
        real(wp), dimension(:,:,:,:,:), intent(inout)              :: rhs_pb, rhs_mv
        integer, intent(in)                                        :: t_step
        real(wp), parameter                                        :: c_abs(3) = [0._wp, 1._wp, 0.5_wp]
        integer                                                    :: kc, pblk, clev, sub, s
        real(wp)                                                   :: th

        clev = plev + 1
        ! setup each child: its two ghost-lerp sources from its own parent's substep endpoints (parent-fine frame) + zeroed
        ! registers
        do kc = 1, amr_num_blocks
            if (amr_block_level(kc) /= clev) cycle
            call s_amr_select_slot(kc)  ! amr_cur = kc; mirrors (isect already parent-fine)
            pblk = f_amr_parent_block(kc)
            if (.not. (amr_rank_owns_block .or. amr_block_owner(pblk) == proc_rank)) cycle
            ! Two P2P pairs (parent @ t_a, then @ t_b), so both owners must arrive or the receiver never posts. The parent owner
            ! packs and sends from its own slot; the child owner receives without naming the parent field (amr_slots(pblk) is
            ! unallocated there). Co-located (np=1, or parent and child on one rank) takes the local device-copy path.
            ! Both sends carry tag amr_cur; MPI non-overtaking on a fixed (source, tag, comm) keeps t_a ahead of t_b.
            if (amr_block_owner(pblk) == proc_rank) then
                call s_amr_gather_from_parent_field_stor(amr_cur, pblk, amr_loc_of(pblk), .false.)  ! parent @ t_a (device C/F fill)
                call s_amr_gather_send_flush()  ! this site has blocking semantics; no drain follows this loop
            else
                call s_amr_recv_parent_patch(pblk, .false.)
            end if
            if (amr_rank_owns_block) call s_amr_fill_fine_ghosts_gsta(amr_cg, amr_loc_of(kc))
            if (amr_block_owner(pblk) == proc_rank) then
                call s_amr_gather_from_parent_field_cons(amr_cur, pblk, amr_loc_of(pblk), .false.)  ! parent @ t_b (device C/F fill)
                call s_amr_gather_send_flush()  ! this site has blocking semantics; no drain follows this loop
            else
                call s_amr_recv_parent_patch(pblk, .false.)
            end if
            if (.not. amr_rank_owns_block) cycle
            call s_amr_fill_fine_ghosts_gstb(amr_cg, amr_loc_of(kc))
            call s_amr_zero_fine_registers()
        end do
        ! advance the level transposed (every level-clev block through each substep together, with the level-clev seam halo
        ! interposed), exactly as s_amr_advance_fine_subcycle_all does at level 1. Advancing each child's whole subtree in turn
        ! would leave adjacent blocks unable to see each other, so their shared face would carry mismatched fluxes. The halo is
        ! level-filtered because this runs inside one of the parents' substeps, when level plev is mid-substep and must not be
        ! touched.
        do sub = 1, 2
            do s = 1, 3
                th = (real(sub - 1, wp) + c_abs(s))*0.5_wp
                do kc = 1, amr_num_blocks
                    if (amr_block_level(kc) /= clev) cycle
                    call s_amr_select_slot(kc)
                    if (.not. amr_rank_owns_block) cycle
                    call s_amr_subtree_stage_lerp(s, th)
                end do
                call s_phase_tic(PH_SEAM)
                call s_amr_fine_fine_halo(clev)
                do kc = 1, amr_num_blocks
                    call s_phase_toc(PH_SEAM)
                    if (amr_block_level(kc) /= clev) cycle
                    call s_amr_select_slot(kc)
                    if (.not. amr_rank_owns_block) cycle
                    call s_amr_subtree_stage_advance(dt_sub*0.5_wp, coefs, bc_type, q_T_sf, pb_in, rhs_pb, mv_in, rhs_mv, t_step, &
                                                     & s, th)
                end do
            end do
            ! every level-clev block is now at its own t_b with t_a in q_cons_stor: recurse into level clev+1 within this substep
            if (amr_max_level >= clev + 1) call s_amr_advance_children(clev, dt_sub*0.5_wp, coefs, bc_type, q_T_sf, pb_in, &
                & rhs_pb, mv_in, rhs_mv, t_step)
        end do
        ! fold each child back into its parent (relax the fine phase first, matching the driver's relax -> restrict order)
        do kc = 1, amr_num_blocks
            if (amr_block_level(kc) /= clev) cycle
            call s_amr_select_slot(kc)
            ! The restrict and the reflux are each a P2P pair when child and parent are on different ranks, so both participants
            ! must reach them or the receiver never posts and the pair deadlocks. Owner-only work (relax) stays behind the guard.
            if (amr_rank_owns_block) then
                if (relax) call s_amr_relax_fine()
            end if
            if (amr_rank_owns_block .or. amr_block_owner(f_amr_parent_block(kc)) == proc_rank) then
                call s_phase_tic(PH_RSRFP)
                call s_amr_restrict_to_parent()
                call s_amr_reflux_to_parent(dt_sub, .true.)
            end if
            call s_phase_toc(PH_RSRFP)
        end do

    end subroutine s_amr_advance_children

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
