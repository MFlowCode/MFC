!>
!!@file
!!@brief Contains module m_amr_advance

#! AMD OpenMP lane: assert allocatables present on every kernel here (see OMP_DEFAULT_STR). A kernel naming an unallocated
#! module array aborts, so conditionally allocated ones (sw_jac/jac: igr) are named only under their allocation condition.
#:set MFC_OMP_PRESENT_ALLOCATABLE = True
#:include 'macros.fpp'

!> @brief Fine-block stage advance and the immersed-boundary/Lagrange fine services.
module m_amr_advance

    use m_derived_types  ! scalar_field, t_box, int_bounds_info
    use m_global_parameters
    use m_constants, only: model_eqns_6eq
    use m_mpi_proxy, only: s_mpi_abort
    use m_rhs, only: s_compute_rhs
    use m_phase_timing
    use m_ibm, only: s_ibm_setup_fine, s_ibm_swap_to_fine, s_ibm_restore_from_fine, s_ibm_correct_state, s_ibm_load_fine_markers, &
        & s_update_mib, s_ibm_reset_fine_park, moving_immersed_boundary_flag, num_gps, f_ibm_fine_gps, s_ibm_tally_gps_by_body, &
        & s_ibm_check_bodies_on_one_rank
    use m_amr_state
    use m_amr_distribution
    use m_amr_store
    use m_amr_exchange
    use m_amr_frame
    use m_amr_transfer

    implicit none

    private
    public :: s_amr_fine_stage_advance_batched, s_amr_setup_ib

contains

    !> Compute the fine-grid IB state (markers/ghost points/levelset) for every active block from the body geometry (static-body
    !! AMR), once after the coarse IB setup: swap grid and IB globals to each owned block, run the fine IB pipeline, restore.
    impure subroutine s_amr_setup_ib()

        integer                      :: islot, save_cur
        integer(kind=8), allocatable :: body_gps(:)

        if (.not. amr .or. .not. ib) return
        save_cur = amr_cur
        allocate (body_gps(max(num_gbl_ibs, 1))); body_gps = 0_8
        call s_ibm_reset_fine_park()
        do islot = 1, amr_num_blocks
            call s_amr_select_slot(islot)
            if (.not. amr_rank_owns_block) cycle
            call s_amr_swap_to_fine()
            call s_ibm_swap_to_fine(islot, gps_on_device=.false.)
            call s_ibm_setup_fine()
            call s_ibm_tally_gps_by_body(body_gps)
            call s_ibm_restore_from_fine(islot)
            call s_amr_restore_coarse()
        end do
        call s_amr_select_slot(save_cur)

        ! The fine-IB image-point stencil is not decomposition-exact across a rank seam: abort if one body's fine ghost points
        ! appear on more than one rank.
        call s_ibm_check_bodies_on_one_rank(body_gps)

    end subroutine s_amr_setup_ib

    !> Apply the IB state correction on the current fine block after its RK update, mirroring the coarse s_ibm_correct_state.
    !! amr_cur is set by the caller.
    impure subroutine s_amr_ib_correct_fine(q_prim_b, pb_in, mv_in)

        !> the q_prim the block's RHS pass filled (pooled scratch)
        type(scalar_field), dimension(1:sys_size), intent(inout) :: q_prim_b
        real(stp), dimension(:,:,:,:,:), intent(inout)           :: pb_in, mv_in  !< unread: amr prohibits qbmm
        integer                                                  :: loc_cur  !< plain scalars for device routines (NVHPC -Minline)

        if (.not. ib) return
        if (.not. amr_rank_owns_block) return
        ! the swap moves the whole marker field through the host, and a block with no ghost point needs no correction
        if (f_ibm_fine_gps(amr_cur) == 0) return
        call s_amr_swap_to_fine()
        call s_ibm_swap_to_fine(amr_cur, gps_on_device=.true.)
        loc_cur = amr_loc_of(amr_cur)
        call s_amr_br_load(loc_cur)
        call s_ibm_correct_state(amr_cons_br, q_prim_b, pb_in, mv_in)
        loc_cur = amr_loc_of(amr_cur)
        call s_amr_br_store(loc_cur)
        call s_ibm_restore_from_fine(amr_cur)
        call s_amr_restore_coarse()

    end subroutine s_amr_ib_correct_fine

    !> Rebuild the current fine block's IB state from the moving body's position (prescribed motion, moving_ibm==1) by reusing the
    !! coarse s_update_mib on the swapped-in slot. Must precede s_amr_ib_correct_fine.
    impure subroutine s_amr_update_mib_fine()

        integer :: i, blo(3), bhi(3)
        logical :: inside, moving_here

        if (.not. ib) return
        if (.not. amr_rank_owns_block) return
        moving_here = .false.
        ! A moving body must stay inside its block (an edge overlap would silently clip the forcing). Under dynamic regrid the
        ! expansion left margin max(amr_buf,4) and body + image-point stencil reach (2 coarse cells) must stay contained between
        ! regrids; on a static block only the body bbox must stay inside. Consecutive contained positions keep every sub-time
        ! interpolate contained (boxes are convex in the linearly-moving corners).
        if (any(patch_ib(1:num_ibs)%moving_ibm /= 0)) then
            do i = 1, num_ibs
                if (patch_ib(i)%moving_ibm == 0) cycle
                call s_amr_body_bbox(i, merge(2, 0, amr_regrid_int > 0), blo, bhi)
                if (.not. f_amr_boxes_overlap(blo, bhi, amr_slots(amr_cur)%region%lo, amr_slots(amr_cur)%region%hi)) cycle
                inside = all((blo >= amr_slots(amr_cur)%region%lo .and. bhi <= amr_slots(amr_cur)%region%hi) .or. .not. amr_dim)
                if (.not. inside) then
                    write (0, '(A,I0,A,3ES11.3,A,3ES11.3,A,6I6,A,6I6)') ' [amr] body ', i, ' at ', patch_ib(i)%x_centroid, &
                           & patch_ib(i)%y_centroid, patch_ib(i)%z_centroid, ' vel ', patch_ib(i)%vel(1:3), ' bbox(lo,hi) ', blo, &
                           & bhi, ' block(lo,hi) ', amr_slots(amr_cur)%region%lo, amr_slots(amr_cur)%region%hi
                    call s_mpi_abort('amr with moving ib: the body reached the fine-block boundary; ' &
                                     & // 'under dynamic regrid reduce amr_regrid_int or increase amr_buf, ' &
                                     & // 'for a static block enlarge it to contain the trajectory')
                end if
                moving_here = .true.
            end do
        end if
        ! the recompute rebuilds the whole slot for every body, so skip blocks with no moving body over them
        if (.not. moving_here) return
        call s_amr_swap_to_fine()
        call s_ibm_swap_to_fine(amr_cur, gps_on_device=.true.)
        call s_update_mib(num_ibs)
        call s_ibm_restore_from_fine(amr_cur)
        call s_amr_restore_coarse()

    end subroutine s_amr_update_mib_fine

    !> Batched fine advance: owned fine blocks are advanced in batches of up to amr_bat_lim blocks of the same level and (padded)
    !! extent, stacked two ghost shells apart along the last active dimension, so one s_compute_rhs call and one RK kernel cover a
    !! batch. Blocks are independent, so grouping order is free. Members read the leader's coordinate arrays in the non-stacked
    !! dimensions (see s_amr_init_report).
    impure subroutine s_amr_fine_stage_advance_batched(s, coefs, bc_type, q_T_sf, pb_in, rhs_pb, mv_in, rhs_mv, t_step)

        integer, intent(in) :: s, t_step
        real(wp), intent(in) :: coefs(4)
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        type(scalar_field), intent(inout) :: q_T_sf
        real(stp), dimension(:,:,:,:,:), intent(inout) :: pb_in, mv_in
        real(wp), dimension(:,:,:,:,:), intent(inout) :: rhs_pb, rhs_mv
        integer :: i, g, ibm
        integer, allocatable :: my(:)
        logical, allocatable :: done(:)
        logical :: last_batch
        integer :: h_bat_blk_ibm  !< plain scalars for device routines (NVHPC -Minline)
        real(wp) :: h_coefs_1, h_coefs_2, h_coefs_3, h_coefs_4, h_merge_1_wp_dt_igr

        call s_amr_refresh_my_blocks()
        allocate (done(amr_n_my)); done = .false.
        my = amr_my_blk(1:amr_n_my)
        do i = 1, amr_n_my
            if (done(i)) cycle
            done(i) = .true.
            g = my(i)
            if (amr_adv_level > 0 .and. amr_block_level(g) /= amr_adv_level) cycle  ! subcycling: one level at a time
            call s_amr_form_batch(i, my, done)
            ! no fine block left undone: this batch's restore must push the coarse grid state (the coarse stage reads it)
            last_batch = .not. any(.not. done .and. (amr_adv_level == 0 .or. amr_block_level(my) == amr_adv_level))
            ! leader selected (swap, capture and RK read amr_cur / its extents); members' store columns
            call s_amr_select_slot(g)
            amr_bat_ext = [amr_slots(g)%m, amr_slots(g)%n, amr_slots(g)%p]
            call s_amr_bat_layout()
            do ibm = 1, amr_bat_n
                amr_bat_loc(ibm) = amr_loc_of(amr_bat_blk(ibm))
                amr_bat_mext(:,ibm) = [amr_slots(amr_bat_blk(ibm))%m, amr_slots(amr_bat_blk(ibm))%n, amr_slots(amr_bat_blk(ibm))%p]
            end do
            $:GPU_UPDATE(device='[amr_bat_loc, amr_bat_mext, amr_bat_off]')
            ! step-entry backup for the SSP-RK combination, per member
            if (s == 1) then
                do ibm = 1, amr_bat_n
                    h_bat_blk_ibm = amr_bat_blk(ibm)
                    call s_amr_copy_fine_fields(h_bat_blk_ibm)
                end do
            end if
            amr_in_fine_advance = .true.
            call s_phase_tic(PH_SWAP)
            call s_amr_swap_to_fine(widen=.true.)  ! the leader's grid, extended into the slab (amr_bat_n > 1)
            call s_phase_toc(PH_SWAP)
            call s_phase_tic(PH_RHS)
            call s_amr_br_load_batch(amr_bat_n)
            ! each member's fine markers at its slab offset, for the RHS body-cell zeroing
            if (ib) call s_ibm_load_fine_markers(amr_bat_n, amr_bat_blk(1:amr_bat_n), amr_bat_mext(:,1:amr_bat_n), amr_bat_sd, &
                & amr_bat_w)
            call s_compute_rhs(amr_cons_br, q_T_sf, amr_scr_prim, bc_type, amr_scr_rhs, pb_in, rhs_pb, mv_in, rhs_mv, t_step, s)
            call s_phase_toc(PH_RHS)
            call s_phase_tic(PH_SWAP)
            call s_amr_restore_coarse(sync_device=last_batch)
            call s_phase_toc(PH_SWAP)
            amr_in_fine_advance = .false.
            call s_phase_tic(PH_RK)
            ! IGR folds dt into its RHS, so the update multiplies by 1 there
            h_coefs_1 = coefs(1)
            h_coefs_2 = coefs(2)
            h_coefs_3 = coefs(3)
            h_coefs_4 = coefs(4)
            h_merge_1_wp_dt_igr = merge(1._wp, dt, igr)
            call s_amr_fine_rk_update_batch(amr_bat_n, amr_scr_rhs, h_coefs_1, h_coefs_2, h_coefs_3, h_coefs_4, h_merge_1_wp_dt_igr)
            if (ib .or. bodyForces .or. cont_damage .or. (model_eqns == model_eqns_6eq .and. (.not. relax))) then
                call s_amr_batch_post_rk_hooks(coefs, pb_in, mv_in)
            end if
            call s_phase_toc(PH_RK)
        end do
        amr_bat_n = 0
        deallocate (done)

    end subroutine s_amr_fine_stage_advance_batched

    !> Start a batch at leader my(i): add later same-level blocks that fit the leader's extent, wasting <= amr_bat_pad of the
    !! member's cells; marks members done.
    impure subroutine s_amr_form_batch(i, my, done)

        integer, intent(in)    :: i, my(:)
        logical, intent(inout) :: done(:)
        integer                :: j, g, h

        g = my(i)
        amr_bat_n = 1; amr_bat_blk(1) = g
        do j = i + 1, amr_n_my
            if (amr_bat_n == amr_bat_lim) exit
            h = my(j)
            if (done(j) .or. amr_block_level(h) /= amr_block_level(g)) cycle
            ! padded membership: no larger than the leader in any dim, wasting <= amr_bat_pad of the member's cells
            if (amr_slots(h)%m > amr_slots(g)%m .or. amr_slots(h)%n > amr_slots(g)%n .or. amr_slots(h)%p > amr_slots(g)%p) cycle
            if (real((amr_slots(g)%m + 1)*(amr_slots(g)%n + 1)*(amr_slots(g)%p + 1) - (amr_slots(h)%m + 1)*(amr_slots(h)%n + 1) &
                & *(amr_slots(h)%p + 1), wp) > amr_bat_pad*real((amr_slots(h)%m + 1)*(amr_slots(h)%n + 1)*(amr_slots(h)%p + 1), &
                & wp)) cycle
            amr_bat_n = amr_bat_n + 1; amr_bat_blk(amr_bat_n) = h; done(j) = .true.
        end do

    end subroutine s_amr_form_batch

    !> Per-member post-RK hooks of the coarse stage (q_cons_ts(1), m_time_steppers).
    impure subroutine s_amr_batch_post_rk_hooks(coefs, pb_in, mv_in)

        real(wp), intent(in)                           :: coefs(4)
        real(stp), dimension(:,:,:,:,:), intent(inout) :: pb_in, mv_in
        integer                                        :: ibm, nb

        ! the coarse stage's post-RK hooks (q_cons_ts(1), m_time_steppers) run here per member in its own frame and in the
        ! coarse order: 6-equation pressure relaxation, body force, moving-body rebuild, body/ghost-cell correction,
        ! continuum-damage clamp. Each reads only the member's cells, so member order is free. Unsupported hooks are
        ! rejected with amr in case_validator. amr_bat_n = 1 while visiting members: s_amr_swap_to_fine would otherwise
        ! extend the grid into the slab, and the hooks need the member's extents.

        nb = amr_bat_n; amr_bat_n = 1
        do ibm = 1, nb
            call s_amr_select_slot(amr_bat_blk(ibm))
            if (model_eqns == model_eqns_6eq .and. (.not. relax)) call s_amr_pressure_relax_fine()
            ! same ldt as the coarse stage: rk_coef(s, 3)*dt/rk_coef(s, 4)
            if (bodyForces) call s_amr_bodyforces_fine(coefs(3)*dt/coefs(4))
            if (ib) then
                if (moving_immersed_boundary_flag) call s_amr_update_mib_fine()
                call s_amr_bat_member_prim(ibm, amr_scr_prim, amr_scr_prim_blk)
                call s_amr_ib_correct_fine(amr_scr_prim_blk, pb_in, mv_in)
            end if
            if (cont_damage) call s_amr_cont_damage_fine()
        end do
        amr_bat_n = nb

    end subroutine s_amr_batch_post_rk_hooks

    !> The batch's slab layout: members stacked along the last active dim, stride amr_bat_w.
    impure subroutine s_amr_bat_layout()

        integer :: ibm

        amr_bat_sd = num_dims
        amr_bat_w = amr_bat_ext(amr_bat_sd) + 2*buff_size + 1
        do ibm = 1, amr_bat_n
            amr_bat_off(:,ibm) = 0
            amr_bat_off(amr_bat_sd, ibm) = (ibm - 1)*amr_bat_w
        end do

    end subroutine s_amr_bat_layout

    !> Device RK stage update over the batch: q = (c1*q + c2*q_stor + c3*dt_in*rhs)/c4 (compute in wp, store stp), reading each
    !! member's rhs at its slab offset. Keep in lockstep with the non-IGR rk_coef form in s_tvd_rk (m_time_steppers).
    impure subroutine s_amr_fine_rk_update_batch(nb, q_rhs, c1, c2, c3, c4, dt_in)

        integer, intent(in)                                 :: nb
        type(scalar_field), dimension(sys_size), intent(in) :: q_rhs
        real(wp), intent(in)                                :: c1, c2, c3, c4, dt_in
        integer                                             :: ibm, i, fi, fj, fk, fm, fn, fp, loc

        fm = amr_bat_ext(1); fn = amr_bat_ext(2); fp = amr_bat_ext(3)
        $:GPU_PARALLEL_LOOP(collapse=5, private='[loc]', copyin='[nb, fm, fn, fp]')
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
                                        & wp) + c3*dt_in*real(q_rhs(i)%sf(fi + amr_bat_off(1, ibm), fj + amr_bat_off(2, ibm), &
                                        & fk + amr_bat_off(3, ibm)), wp))/c4
                        end do
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_fine_rk_update_batch

end module m_amr_advance
