!>
!!@file
!!@brief Contains module m_amr_stage

#:include 'macros.fpp'

!> @brief The AMR hooks of the SSP-RK stage loop (m_time_steppers): the stage-top coarse halo, the fine-level stage, and the
!! end-of-step fold of the fine solution into the level-0 covered cells.
module m_amr_stage

    use m_derived_types  ! scalar_field, integer_field
    use m_global_parameters
    use m_phase_timing
    use m_amr_registers, only: s_amr_apply_reflux
    use m_amr

    implicit none

    private
    public :: s_amr_stage_begin, s_amr_stage_fine, s_amr_step_fold, s_amr_subcycle_fine

    real(wp), parameter :: c_stage(3) = [0._wp, 1._wp, 0.5_wp]  !< SSP-RK3 stage times in units of the step

contains

    !> Stage top: hoist the cons halo (the AMR fill needs it) so the RHS converts over the buffered domain instead of exchanging a
    !! prim halo on the same faces (pointwise, so byte-identical). Only where the cons halo carries everything the prim halo did (no
    !! q_T_sf, no igr).
    impure subroutine s_amr_stage_begin(q_cons)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons

        amr_cons_ghosts_valid = amr .and. amr_xchg_coarse_ghosts .and. (.not. bubbles_lagrange) .and. (.not. chemistry) &
            & .and. (.not. igr)
        if (.not. amr_cons_ghosts_valid) return
        call s_phase_tic(PH_HALO)
        call s_amr_exchange_coarse_cons_halo(q_cons)
        call s_phase_toc(PH_HALO)

    end subroutine s_amr_stage_begin

    !> The fine-level stage, with q_cons still the coarse stage-entry state: (1) fill every block's ghosts top-down as per-level
    !! waves after one coarse cons halo, (2) overwrite fine-fine seam ghosts with neighbours' stage-entry interiors (posted before
    !! the fills, drained after), (3) advance all owned blocks in batches, (4) reflux level-1 blocks into the coarse rhs (one face
    !! wave plus one batched apply, order-free). A level>=2 block's coarse side is its parent, refluxed in s_amr_step_fold.
    impure subroutine s_amr_stage_fine(s, t_step, coefs, q_cons, q_T_sf, bc_type, rhs_vf, pb_in, rhs_pb, mv_in, rhs_mv)

        integer, intent(in)                                        :: s, t_step
        real(wp), intent(in)                                       :: coefs(4)
        type(scalar_field), dimension(sys_size), intent(inout)     :: q_cons, rhs_vf
        type(scalar_field), intent(inout)                          :: q_T_sf
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        real(stp), dimension(:,:,:,:,:), intent(inout)             :: pb_in, mv_in
        real(wp), dimension(:,:,:,:,:), intent(inout)              :: rhs_pb, rhs_mv

        if (.not. amr) return
        call s_phase_tic(PH_HALO)
        if (amr_xchg_coarse_ghosts .and. .not. amr_cons_ghosts_valid) call s_amr_exchange_coarse_cons_halo(q_cons)
        call s_phase_toc(PH_HALO)
        amr_cons_ghosts_valid = .false.
        call s_amr_fill_levels(q_cons, 1, amr_num_levels)
        call s_amr_fine_stage_advance_batched(s, coefs, bc_type, q_T_sf, pb_in, rhs_pb, mv_in, rhs_mv, t_step)
        call s_phase_tic(PH_REFLUX)
        call s_amr_reflux_faces_wave()
        call s_amr_apply_reflux(rhs_vf)
        call s_phase_toc(PH_REFLUX)
        call s_amr_select_slot(1)  ! the next stage's coarse RHS captures creg into slot 1

    end subroutine s_amr_stage_fine

    !> Subcycled fine advance (amr_subcycle): after level 0's step (q_old = step-start state, q_new = end state, coarse registers
    !! holding its dt-weighted flux integral) level 1 takes two substeps of dt/2. The level-1 patch rings of q_old and q_new are
    !! gathered once (s_amr_pc_gather); each stage blends them in time, fills the fine ghosts and advances level 1 with dt/2,
    !! accumulating fine registers over both substeps; level 2 (if any) then takes its two substeps (s_amr_subcycle_l2). The
    !! registers then correct the coarse state (Berger-Colella) and the caller's fold restricts.
    impure subroutine s_amr_subcycle_fine(q_new, q_old, coefs, q_T_sf, bc_type, pb_in, rhs_pb, mv_in, rhs_mv, t_step)

        type(scalar_field), dimension(sys_size), intent(inout)     :: q_new, q_old
        real(wp), intent(in)                                       :: coefs(:,:)
        type(scalar_field), intent(inout)                          :: q_T_sf
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        real(stp), dimension(:,:,:,:,:), intent(inout)             :: pb_in, mv_in
        real(wp), dimension(:,:,:,:,:), intent(inout)              :: rhs_pb, rhs_mv
        integer, intent(in)                                        :: t_step
        real(wp)                                                   :: dt0
        integer                                                    :: k, s

        amr_cons_ghosts_valid = .false.  ! level 0's last stage may have set it; the fine RHS exchanges in lock-step
        call s_amr_pc_gather(q_old, 0)
        call s_amr_pc_gather(q_new, 1)
        dt0 = dt
        dt = 0.5_wp*dt0
        $:GPU_UPDATE(device='[dt]')
        do k = 0, 1
            amr_sub_first = k == 0
            amr_adv_level = 1
            do s = 1, 3
                call s_amr_fill_levels(q_new, 1, 1, 0.5_wp*(real(k, wp) + c_stage(s)))
                call s_amr_fine_stage_advance_batched(s, coefs(s,:), bc_type, q_T_sf, pb_in, rhs_pb, mv_in, rhs_mv, t_step)
            end do
            if (amr_max_level >= 2) call s_amr_subcycle_l2(q_new, coefs, q_T_sf, bc_type, pb_in, rhs_pb, mv_in, rhs_mv, t_step, &
                & 0.5_wp*real(k + 1, wp), dt0)
        end do
        amr_adv_level = 0
        dt = dt0
        $:GPU_UPDATE(device='[dt]')
        call s_phase_tic(PH_REFLUX)
        call s_amr_reflux_faces_wave()
        call s_amr_apply_reflux(q_new)
        call s_phase_toc(PH_REFLUX)
        call s_amr_select_slot(1)

    end subroutine s_amr_subcycle_fine

    !> Level 2's two substeps of dt0/4 after a level-1 substep ending at coarse time fraction th_end: refill level-1 ghosts at
    !! th_end (level-2 parent patches can reach them), fill each level-2 stage from parents blended between their substep-start
    !! backup and end state (amr_par_theta), then sync level 2 into level 1.
    impure subroutine s_amr_subcycle_l2(q_new, coefs, q_T_sf, bc_type, pb_in, rhs_pb, mv_in, rhs_mv, t_step, th_end, dt0)

        type(scalar_field), dimension(sys_size), intent(inout)     :: q_new
        real(wp), intent(in)                                       :: coefs(:,:)
        type(scalar_field), intent(inout)                          :: q_T_sf
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        real(stp), dimension(:,:,:,:,:), intent(inout)             :: pb_in, mv_in
        real(wp), dimension(:,:,:,:,:), intent(inout)              :: rhs_pb, rhs_mv
        integer, intent(in)                                        :: t_step
        real(wp), intent(in)                                       :: th_end, dt0
        integer                                                    :: j, s

        call s_amr_fill_levels(q_new, 1, 1, th_end)
        dt = 0.25_wp*dt0
        $:GPU_UPDATE(device='[dt]')
        amr_adv_level = 2
        do j = 0, 1
            amr_sub_first = j == 0
            do s = 1, 3
                amr_par_theta = 0.5_wp*(real(j, wp) + c_stage(s))
                call s_amr_fill_levels(q_new, 2, 2)
                call s_amr_fine_stage_advance_batched(s, coefs(s,:), bc_type, q_T_sf, pb_in, rhs_pb, mv_in, rhs_mv, t_step)
            end do
        end do
        amr_par_theta = 1._wp
        dt = 0.5_wp*dt0
        $:GPU_UPDATE(device='[dt]')
        call s_phase_tic(PH_RESTR)
        call s_amr_freg_wave()
        call s_amr_sync_level(2, 1._wp)
        call s_phase_toc(PH_RESTR)

    end subroutine s_amr_subcycle_l2

    !> Ghost fill of levels lev_lo..lev_hi: the level-1 wave from the coarse q_cons (or, with theta, the patch cache at that coarse
    !! time fraction), then the parent waves. The fine-fine seam reads stage-entry interiors only, so it is posted before the fills
    !! and drained after.
    impure subroutine s_amr_fill_levels(q_cons, lev_lo, lev_hi, theta)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons
        integer, intent(in)                                    :: lev_lo, lev_hi
        real(wp), intent(in), optional                         :: theta
        integer                                                :: ilev

        call s_amr_fine_fine_post()
        if (lev_lo == 1) then
            if (present(theta)) then
                call s_amr_pc_fill(theta)
            else
                call s_amr_stage_fill_wave(q_cons)
            end if
        end if
        do ilev = max(lev_lo, 2), lev_hi
            call s_amr_parent_fill_wave(ilev)
        end do
        call s_phase_tic(PH_SEAM)
        call s_amr_fine_fine_drain()
        call s_phase_toc(PH_SEAM)

    end subroutine s_amr_fill_levels

    !> Gather every owned level-1 block's coarse patch ring from q into the pooled amr_cg and keep it as time level t (0: step
    !! start, 1: step end), once per coarse step; stages then blend in time with no wave or halo.
    impure subroutine s_amr_pc_gather(q, t)

        type(scalar_field), dimension(sys_size), intent(inout) :: q
        integer, intent(in)                                    :: t

        call s_amr_l1_fill_exchange(q, .false.)
        call s_amr_l1_pool_gather(q)
        call s_amr_pc_copy(t, 0._wp)

    end subroutine s_amr_pc_gather

    !> Level-1 ghost fill at coarse time fraction theta: blend the cached time levels into the pool, then the batched fill.
    impure subroutine s_amr_pc_fill(theta)

        real(wp), intent(in) :: theta

        call s_phase_tic(PH_LERP)
        call s_amr_pc_copy(-1, theta)
        call s_phase_toc(PH_LERP)
        call s_amr_l1_pool_fill()

    end subroutine s_amr_pc_fill

    !> The used pool (amr_cg_nuse patches) against the cache amr_pc: t >= 0 stores it as time level t (t = 0 sizes the cache); t < 0
    !! writes the pool = (1 - theta) * level 0 + theta * level 1.
    impure subroutine s_amr_pc_copy(t, theta)

        integer, intent(in)  :: t
        real(wp), intent(in) :: theta
        integer              :: i, g1, g2, g3, n1, n2, n3, nl, o, c

        n1 = amr_cpat_hi(1) + 1; n2 = amr_cpat_hi(2) + 1; n3 = amr_cg_nuse*(amr_cpat_hi(3) + 1)
        nl = sys_size*n1*n2*n3
        if (t == 0) then
            if (allocated(amr_pc)) then
                if (size(amr_pc) < 2*nl) then
                    @:DEALLOCATE(amr_pc)
                end if
            end if
            if (.not. allocated(amr_pc)) then
                @:ALLOCATE(amr_pc(1:max(4*nl, 1)))
            end if
        end if
        if (t >= 0) then
            o = t*nl
            $:GPU_PARALLEL_LOOP(collapse=4, private='[c]')
            do i = 1, sys_size
                do g3 = 0, n3 - 1
                    do g2 = 0, n2 - 1
                        do g1 = 0, n1 - 1
                            c = 1 + g1 + n1*(g2 + n2*(g3 + n3*(i - 1)))
                            amr_pc(o + c) = amr_cg(i)%sf(g1, g2, g3)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            return
        end if
        $:GPU_PARALLEL_LOOP(collapse=4, private='[c]')
        do i = 1, sys_size
            do g3 = 0, n3 - 1
                do g2 = 0, n2 - 1
                    do g1 = 0, n1 - 1
                        c = 1 + g1 + n1*(g2 + n2*(g3 + n3*(i - 1)))
                        amr_cg(i)%sf(g1, g2, g3) = (1._wp - theta)*amr_pc(c) + theta*amr_pc(nl + c)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_pc_copy

    !> End of step: fold the fine solution into the level-0 covered cells (the only deliberate level-0 write) and Berger-Colella
    !! state-reflux each level>=2 block into its parent. Bottom-up (per-level waves), so covered cells hold the finest data.
    impure subroutine s_amr_step_fold(q_cons)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons

        call s_phase_tic(PH_RESTR)
        ! level>=2 freg exchange: the registers are final after the advance (subcycling synced level 2 already)
        if (.not. amr_subcycle) call s_amr_freg_wave()
        call s_phase_toc(PH_RESTR)
        call s_amr_restrict_wave(q_cons, dt)
        call s_amr_select_slot(1)

    end subroutine s_amr_step_fold

end module m_amr_stage
