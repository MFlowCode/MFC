!>
!!@file
!!@brief Contains module m_amr_stage

#:include 'macros.fpp'

!> @brief The AMR hooks of the SSP-RK stage loop (m_time_steppers): the stage-top coarse halo, the fine-level stage (fills, batched
!! advance, reflux into the coarse rhs), the L0-tile coarse update and the end-of-step fold of the fine solution into the level-0
!! covered cells.
module m_amr_stage

    use m_derived_types  ! scalar_field, integer_field
    use m_global_parameters
    use m_phase_timing
    use m_rhs, only: s_compute_rhs
    use m_active_box, only: ab_active
    use m_amr_registers, only: s_amr_apply_reflux
    use m_amr

    implicit none

    private
    public :: s_amr_stage_begin, s_amr_stage_fine, s_amr_l0_stage_update, s_amr_step_fold

contains

    !> Stage top. Coexist: tiles are the authoritative store, so refresh the L0 staging buffer from the tile interiors before the
    !! coarse RHS and the fine coarse-patch fills read it. rhs_now: the coarse RHS runs at the stage top; pure-L0 skips it (the
    !! tiles run their own), and coexist defers it past the fine advance (s_amr_stage_fine), where cross-rank stage skew is best
    !! absorbed - same inputs, so byte-identical - except with chemistry, whose coarse-vs-fine q_T_sf write order must not swap.
    !! When it runs now, the AMR cons halo and the RHS's prim halo would exchange the same stage-entry state on the same faces:
    !! hoist the cons halo here and let the RHS convert over the buffered domain (pointwise, so byte-identical), halving the
    !! base-grid SENDRECVs per step; only where the cons halo carries everything the prim halo did (no q_T_sf, no igr path).
    impure subroutine s_amr_stage_begin(q_cons, rhs_now)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons
        logical, intent(out)                                   :: rhs_now

        if (amr .and. l0_ntile > 0) call s_l0_scatter_tiles_to_coarse(q_cons)
        rhs_now = l0_ntile == 0 .or. (amr .and. chemistry)
        amr_cons_ghosts_valid = rhs_now .and. amr .and. amr_xchg_coarse_ghosts .and. (.not. bubbles_lagrange) &
            & .and. (.not. chemistry) .and. (.not. igr) .and. (.not. ab_active)
        if (.not. amr_cons_ghosts_valid) return
        call s_phase_tic(PH_HALO)
        call s_amr_exchange_coarse_cons_halo(q_cons)
        call s_phase_toc(PH_HALO)

    end subroutine s_amr_stage_begin

    !> The fine-level stage, with q_cons still the coarse stage-entry state. Phase 1 fills every block's ghost shell top-down as
    !! per-level exchange waves (the level-1 wave, then one parent-gather wave per level ascending, so each level's sources are
    !! complete before its wave), after one coarse cons halo for every block's prolongation; phase 2 overwrites the fine-fine seam
    !! ghosts with the neighbours' stage-entry interiors (posted before the fills, drained after); phase 3 advances every owned
    !! block in batches of equal shape; phase 4 refluxes the level-1 blocks into the coarse rhs as one face wave plus one batched
    !! apply (both order-free: disjoint register slots, disjoint corrections by the merge invariant). A level>=2 block's coarse side
    !! is its parent, refluxed after the stage loop (s_amr_step_fold). Coexist: the deferred coarse RHS runs before the reflux so
    !! creg(L0) exists, and its rhs is zeroed so rhs_vf becomes the pure reflux-delta accumulator that s_l0_add_reflux_to_tiles
    !! routes to each covering tile (the tiles carry their own rhs).
    impure subroutine s_amr_stage_fine(s, t_step, coefs, q_cons, q_T_sf, q_prim_vf, bc_type, rhs_vf, pb_in, rhs_pb, mv_in, rhs_mv)

        integer, intent(in)                                        :: s, t_step
        real(wp), intent(in)                                       :: coefs(4)
        type(scalar_field), dimension(sys_size), intent(inout)     :: q_cons, q_prim_vf, rhs_vf
        type(scalar_field), intent(inout)                          :: q_T_sf
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        real(stp), dimension(:,:,:,:,:), intent(inout)             :: pb_in, mv_in
        real(wp), dimension(:,:,:,:,:), intent(inout)              :: rhs_pb, rhs_mv
        integer                                                    :: ilev

        if (.not. amr) return
        if (l0_ntile > 0 .and. chemistry) call s_amr_zero_rhs(rhs_vf)  ! coexist chemistry: the stage-top RHS ran; zero it now
        call s_phase_tic(PH_HALO)
        if (amr_xchg_coarse_ghosts .and. .not. amr_cons_ghosts_valid) call s_amr_exchange_coarse_cons_halo(q_cons)
        call s_phase_toc(PH_HALO)
        amr_cons_ghosts_valid = .false.
        if (amr_early_seam_post) call s_amr_fine_fine_post()
        call s_amr_stage_fill_wave(q_cons)
        do ilev = 2, amr_num_levels
            call s_amr_parent_fill_wave(ilev)
        end do
        call s_phase_tic(PH_SEAM)
        if (amr_early_seam_post) then
            call s_amr_fine_fine_drain()
        else
            call s_amr_fine_fine_halo()
        end if
        call s_phase_toc(PH_SEAM)
        call s_amr_fine_stage_advance_batched(s, coefs, bc_type, q_T_sf, pb_in, rhs_pb, mv_in, rhs_mv, t_step)
        if (l0_ntile > 0 .and. .not. chemistry) then
            call s_phase_tic(PH_COARSE)
            call s_compute_rhs(q_cons, q_T_sf, q_prim_vf, bc_type, rhs_vf, pb_in, rhs_pb, mv_in, rhs_mv, t_step, s)
            call s_phase_toc(PH_COARSE)
            call s_amr_zero_rhs(rhs_vf)
        end if
        call s_phase_tic(PH_REFLUX)
        call s_amr_reflux_faces_wave()
        call s_amr_apply_reflux(rhs_vf)
        call s_phase_toc(PH_REFLUX)
        call s_amr_select_slot(1)  ! the next stage's coarse RHS captures creg into slot 1

    end subroutine s_amr_stage_fine

    impure subroutine s_amr_zero_rhs(rhs_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        integer                                                :: i, j, k, l

        $:GPU_PARALLEL_LOOP(collapse=4)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        rhs_vf(i)%sf(j, k, l) = 0._wp
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_zero_rhs

    !> The L0-tile coarse update of the stage (in place of the monolithic RK update; byte-identical to it). Tiles own the state
    !! across stages (copied in at stage 1) and L0 is gathered from them only at output. Coexist splits the tile advance so the
    !! fixed-L0-frame reflux delta in rhs_vf reaches each tile before its RK update.
    impure subroutine s_amr_l0_stage_update(s, t_step, coefs, q_cons, q_T_sf, bc_type, rhs_vf, pb_in, rhs_pb, mv_in, rhs_mv)

        integer, intent(in)                                        :: s, t_step
        real(wp), intent(in)                                       :: coefs(4)
        type(scalar_field), dimension(sys_size), intent(inout)     :: q_cons, rhs_vf
        type(scalar_field), intent(inout)                          :: q_T_sf
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        real(stp), dimension(:,:,:,:,:), intent(inout)             :: pb_in, mv_in
        real(wp), dimension(:,:,:,:,:), intent(inout)              :: rhs_pb, rhs_mv

        call s_phase_tic(PH_L0)
        if (s == 1) then
            call s_l0_copy_coarse_to_tiles(q_cons)
            if (l0_migrate_step > 0 .and. t_step == l0_migrate_step) call s_l0_forced_remap()
            ! nested so mod() is never reached at the default interval of 0: Fortran does not short-circuit .and., and amdflang
            ! hoists the integer divide ahead of the guard (SIGFPE)
            if (l0_rebalance_interval > 0 .and. t_step > 0) then
                if (mod(t_step, l0_rebalance_interval) == 0) call s_l0_rebalance(t_step)
            end if
        end if
        if (amr) then
            call s_l0_advance_stage_rhs(s, bc_type, q_T_sf, pb_in, rhs_pb, mv_in, rhs_mv, t_step)
            call s_l0_add_reflux_to_tiles(rhs_vf)
            call s_l0_advance_stage_rk(s, coefs)
        else
            call s_l0_advance_stage(s, coefs, bc_type, q_T_sf, pb_in, rhs_pb, mv_in, rhs_mv, t_step)
        end if
        call s_phase_toc(PH_L0)

    end subroutine s_amr_l0_stage_update

    !> End of step: fold the fine solution into the level-0 covered cells (the only deliberate level-0 write) and Berger-Colella
    !! state-reflux each level>=2 block into its parent. Bottom-up: a child folds into its parent before the parent folds into L0,
    !! so the covered cells hold the finest data; finer levels live at higher slots, so the per-box loop runs in reverse (disjoint
    !! same-level blocks make that bit-identical to forward order for one level). np>1 runs the fold as per-level waves, since a
    !! per-box loop would serialize a P2P chain that scales with the global block count. Coexist: the covered cells are then routed
    !! back to the covering tiles, the authoritative store.
    impure subroutine s_amr_step_fold(q_cons)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons
        integer                                                :: islot

        call s_phase_tic(PH_RESTR)
        call s_amr_freg_wave()  ! the split-ownership level>=2 freg exchange: the registers are final after the advance
        call s_phase_toc(PH_RESTR)
        if (num_procs > 1) then
            call s_amr_restrict_wave(q_cons, dt)
        else
            do islot = amr_num_blocks, 1, -1
                if (amr_block_level(islot) == 0) cycle  ! L0 tile slots advance separately
                call s_amr_select_slot(islot)
                call s_phase_tic(PH_RESTR)
                call s_restrict_fine_to_coarse(q_cons)
                if (amr_block_level(amr_cur) >= 2) call s_amr_reflux_to_parent(dt)
                call s_phase_toc(PH_RESTR)
            end do
        end if
        call s_amr_select_slot(1)
        if (l0_ntile > 0) call s_l0_restrict_to_tiles(q_cons)

    end subroutine s_amr_step_fold

end module m_amr_stage
