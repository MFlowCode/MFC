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
    use m_amr_registers, only: s_amr_apply_reflux, amr_rk_c
    use m_phase_change, only: s_infinite_relaxation_k
    use m_reactive_burn, only: s_reactive_burn_substep
    use m_icpp_patches, only: s_apply_icpp_patches
    use m_assign_variables, only: s_initialize_assign_variables_module
    use m_variables_conversion, only: s_convert_primitive_to_conservative_variables
    use m_mpi_common, only: s_mpi_allreduce_sum
    use m_amr

    implicit none

    private
    public :: s_amr_stage_begin, s_amr_stage_fine, s_amr_step_fold, s_amr_subcycle_fine, s_amr_relax_fine, s_amr_burn_fine, &
        & s_amr_apply_ic_patches, s_amr_fill_levels

    integer, parameter :: fine_op_relax = 1, fine_op_burn = 2  !< cell-local end-of-step operators of s_amr_fine_cell_op

    !> Whether patch_icpp was compared with the level-0 initial data, and reproduces it
    logical :: ic_checked = .false., ic_from_patches = .false.

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
        call s_amr_apply_reflux(rhs_vf, q_cons)
        call s_phase_toc(PH_REFLUX)
        call s_amr_select_slot(1)  ! the next stage's coarse RHS captures creg into slot 1

    end subroutine s_amr_stage_fine

    !> Subcycled fine advance (amr_subcycle): after level 0's step (q_old = step-start state, q_new = end state, coarse registers
    !! holding its dt-weighted flux integral) the fine levels substep recursively (s_amr_subcycle_level), accumulating their
    !! registers. The registers then correct the coarse state (Berger-Colella) and the caller's fold restricts.
    impure subroutine s_amr_subcycle_fine(q_new, q_old, coefs, q_T_sf, bc_type, pb_in, rhs_pb, mv_in, rhs_mv, t_step)

        type(scalar_field), dimension(sys_size), intent(inout)     :: q_new, q_old
        real(wp), intent(in)                                       :: coefs(:,:)
        type(scalar_field), intent(inout)                          :: q_T_sf
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        real(stp), dimension(:,:,:,:,:), intent(inout)             :: pb_in, mv_in
        real(wp), dimension(:,:,:,:,:), intent(inout)              :: rhs_pb, rhs_mv
        integer, intent(in)                                        :: t_step
        real(wp)                                                   :: dt0

        amr_cons_ghosts_valid = .false.  ! level 0's last stage may have set it; the fine RHS exchanges in lock-step
        ! the level-1 patch rings of q_old and q_new, gathered once; each level-1 stage blends them in time
        call s_amr_pc_gather(q_old, 0)
        call s_amr_pc_gather(q_new, 1)
        dt0 = dt
        call s_amr_subcycle_level(1, 0._wp, 1._wp, 0.5_wp*dt0, q_new, coefs, q_T_sf, bc_type, pb_in, rhs_pb, mv_in, rhs_mv, t_step)
        amr_adv_level = 0
        dt = dt0
        $:GPU_UPDATE(device='[dt]')
        call s_phase_tic(PH_REFLUX)
        call s_amr_reflux_faces_wave()
        call s_amr_apply_reflux(q_new, q_old)  ! alpha for the source correction: step start (exact for uniform alpha)
        call s_phase_toc(PH_REFLUX)
        call s_amr_select_slot(1)

    end subroutine s_amr_subcycle_fine

    !> Level lev's two substeps of dtl across the fraction [th0, th0 + dth] of the coarse step. Level 1 fills its ghosts from the
    !! patch cache at that coarse time; a deeper level fills from its parents blended between their substep-start backup and end
    !! state (amr_par_theta). After each substep the next level substeps across it: this level's ghosts are first refilled at the
    !! substep end (the children's parent patches can reach them), and the children are then synced back into this level.
    recursive subroutine s_amr_subcycle_level(lev, th0, dth, dtl, q_new, coefs, q_T_sf, bc_type, pb_in, rhs_pb, mv_in, rhs_mv, &
        & t_step)

        integer, intent(in)                                        :: lev, t_step
        real(wp), intent(in)                                       :: th0, dth, dtl
        type(scalar_field), dimension(sys_size), intent(inout)     :: q_new
        real(wp), intent(in)                                       :: coefs(:,:)
        type(scalar_field), intent(inout)                          :: q_T_sf
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        real(stp), dimension(:,:,:,:,:), intent(inout)             :: pb_in, mv_in
        real(wp), dimension(:,:,:,:,:), intent(inout)              :: rhs_pb, rhs_mv
        integer                                                    :: k, s

        do k = 0, 1
            amr_sub_first = k == 0
            amr_adv_level = lev
            dt = dtl
            $:GPU_UPDATE(device='[dt]')
            do s = 1, size(amr_rk_c)
                ! this substage's time: fraction th0 + dth*(k + c_s)/2 of the coarse step dtl*2/dth
                t_rhs = f_step_time(th0 + 0.5_wp*dth*(real(k, wp) + amr_rk_c(s)), t_step, 2._wp*dtl/dth)
                $:GPU_UPDATE(device='[t_rhs]')
                call s_amr_fill_sub(lev, q_new, 0.5_wp*(real(k, wp) + amr_rk_c(s)), th0, dth)
                call s_amr_fine_stage_advance_batched(s, coefs(s,:), bc_type, q_T_sf, pb_in, rhs_pb, mv_in, rhs_mv, t_step)
            end do
            if (lev >= amr_max_level) cycle
            call s_amr_fill_sub(lev, q_new, 0.5_wp*real(k + 1, wp), th0, dth)
            call s_amr_subcycle_level(lev + 1, th0 + 0.5_wp*dth*real(k, wp), 0.5_wp*dth, 0.5_wp*dtl, q_new, coefs, q_T_sf, &
                                      & bc_type, pb_in, rhs_pb, mv_in, rhs_mv, t_step)
            amr_par_theta = 1._wp
            dt = dtl
            $:GPU_UPDATE(device='[dt]')
            call s_phase_tic(PH_RESTR)
            call s_amr_freg_wave()
            call s_amr_sync_level(lev + 1, 1._wp)
            call s_phase_toc(PH_RESTR)
        end do

    end subroutine s_amr_subcycle_level

    !> Ghost fill of level lev at the fraction th of its parent's current substep, which spans [th0, th0 + dth] of the coarse step.
    impure subroutine s_amr_fill_sub(lev, q_new, th, th0, dth)

        integer, intent(in)                                    :: lev
        type(scalar_field), dimension(sys_size), intent(inout) :: q_new
        real(wp), intent(in)                                   :: th, th0, dth

        if (lev == 1) then
            call s_amr_fill_levels(q_new, 1, 1, th0 + dth*th)
        else
            amr_par_theta = th
            call s_amr_fill_levels(q_new, lev, lev)
        end if

    end subroutine s_amr_fill_sub

    !> Ghost fill of levels lev_lo..lev_hi: the level-1 wave from the coarse q_cons (or, with theta, the patch cache at that coarse
    !! time fraction), then the parent waves. The fine-fine seam reads stage-entry interiors only, so it is posted before the fills
    !! and drained after. With a block on its parent's edge (amr_nest_seam) each level's seams are instead exchanged before the next
    !! level gathers, since that block's parent patch reads them.
    impure subroutine s_amr_fill_levels(q_cons, lev_lo, lev_hi, theta)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons
        integer, intent(in)                                    :: lev_lo, lev_hi
        real(wp), intent(in), optional                         :: theta
        integer                                                :: ilev

        call s_amr_refresh_lists()
        if (.not. amr_nest_seam) call s_amr_fine_fine_post()
        do ilev = lev_lo, max(lev_lo, lev_hi)
            if (amr_nest_seam) call s_amr_fine_fine_post(ilev)
            if (ilev > 1) then
                call s_amr_parent_fill_wave(ilev)
            else if (present(theta)) then
                call s_amr_pc_fill(theta)
            else
                call s_amr_stage_fill_wave(q_cons)
            end if
            if (amr_nest_seam .or. ilev >= lev_hi) then
                call s_phase_tic(PH_SEAM)
                call s_amr_fine_fine_drain()
                call s_phase_toc(PH_SEAM)
            end if
        end do

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

        integer, intent(in) :: t
        real(wp), intent(in) :: theta
        integer :: i, g1, g2, g3, n1, n2, n3, nl, o, c
        real(wp) :: th  !< host copy of theta for the kernel (theta can arrive as an optional dummy of an inlined caller)

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
        th = theta
        $:GPU_PARALLEL_LOOP(collapse=4, private='[c]', firstprivate='[th]')
        do i = 1, sys_size
            do g3 = 0, n3 - 1
                do g2 = 0, n2 - 1
                    do g1 = 0, n1 - 1
                        c = 1 + g1 + n1*(g2 + n2*(g3 + n3*(i - 1)))
                        amr_cg(i)%sf(g1, g2, g3) = (1._wp - th)*amr_pc(c) + th*amr_pc(nl + c)
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

    !> End-of-step phase-change relaxation on the fine blocks (the coarse grid's is in m_start_up), then a restrict so the covered
    !! cells hold the average of the relaxed fine cells. The relaxation is cell-local, so block order is free.
    impure subroutine s_amr_relax_fine(q_cons)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons

        call s_amr_fine_cell_op(fine_op_relax)
        call s_amr_restrict_fold(q_cons)

    end subroutine s_amr_relax_fine

    !> End-of-step operator-split burn (rburn%substeps > 0) on the fine blocks over the full step; the coarse grid's is in
    !! m_time_steppers. Called before s_amr_step_fold, which restricts the burnt fine cells into the covered coarse cells.
    impure subroutine s_amr_burn_fine()

        call s_amr_fine_cell_op(fine_op_burn)

    end subroutine s_amr_burn_fine

    !> Apply the cell-local end-of-step operator op to every owned fine block in its own frame (block order is free).
    impure subroutine s_amr_fine_cell_op(op)

        integer, intent(in) :: op
        integer             :: kk, nb

        call s_amr_refresh_my_blocks()
        nb = amr_bat_n; amr_bat_n = 1  ! one block's extents in the frame, not a batch slab
        do kk = 1, amr_n_my
            call s_amr_select_slot(amr_my_blk(kk))
            call s_amr_fine_enter()
            select case (op)
            case (fine_op_relax)
                call s_infinite_relaxation_k(amr_cons_br)
            case (fine_op_burn)
                call s_reactive_burn_substep(amr_cons_br, dt, idwint)
            end select
            call s_amr_fine_leave()
        end do
        amr_bat_n = nb

    end subroutine s_amr_fine_cell_op

    !> Fresh start: set every owned fine block to the initial condition patches evaluated at its own cell centers, replacing the
    !! data interpolated from level 0, then restrict. Skipped when the patches alone do not reproduce the level-0 data.
    impure subroutine s_amr_apply_ic_patches(q_cons)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons
        integer                                                :: kk, nb, i

        if (.not. ic_checked) call s_amr_check_ic_patches(q_cons)
        if (.not. ic_from_patches) return
        call s_amr_refresh_my_blocks()
        nb = amr_bat_n; amr_bat_n = 1  ! one block's extents in the frame, not a batch slab
        do kk = 1, amr_n_my
            call s_amr_select_slot(amr_my_blk(kk))
            call s_amr_fine_enter()
            ! the patches run on the host and fill the interior only
            do i = 1, sys_size
                $:GPU_UPDATE(host='[amr_cons_br(i)%sf]')
            end do
            call s_amr_ic_on_frame(amr_cons_br)
            do i = 1, sys_size
                $:GPU_UPDATE(device='[amr_cons_br(i)%sf]')
            end do
            call s_amr_fine_leave()
        end do
        amr_bat_n = nb
        call s_amr_restrict_fold(q_cons)

    end subroutine s_amr_apply_ic_patches

    !> Compare the patches evaluated on level 0 with the loaded initial data. They differ when pre_process perturbed or smoothed the
    !! patches or started from old data; hard-coded and STL patches are not evaluated in the simulation.
    impure subroutine s_amr_check_ic_patches(q_cons)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons
        type(scalar_field), dimension(sys_size)                :: q_ic
        real(wp)                                               :: bad, bad_glb, tol
        integer                                                :: i

        ic_checked = .true.
        tol = 1.e3_wp*max(epsilon(1._wp), real(epsilon(1._stp), wp))
        if (precision == 1) tol = max(tol, 1.e3_wp*real(epsilon(1._sp), wp))
        bad = 1._wp
        if (num_patches > 0) then
            if (all(patch_icpp(1:num_patches)%hcid == dflt_int) .and. all(patch_icpp(1:num_patches)%geometry /= icpp_model)) then
                call s_initialize_assign_variables_module()
                bad = 0._wp
                do i = 1, sys_size
                    allocate (q_ic(i)%sf(0:m,0:n,0:p))
                    $:GPU_UPDATE(host='[q_cons(i)%sf]')
                end do
                call s_amr_ic_on_frame(q_ic)
                do i = 1, sys_size
                    if (maxval(abs(q_ic(i)%sf - real(q_cons(i)%sf(0:m,0:n,0:p), wp))) > tol*maxval(abs(q_ic(i)%sf))) bad = 1._wp
                    deallocate (q_ic(i)%sf)
                end do
            end if
        end if
        call s_mpi_allreduce_sum(bad, bad_glb)
        ic_from_patches = bad_glb < 0.5_wp
        if (proc_rank == 0 .and. .not. ic_from_patches) print '(A)', &
            & ' [amr] NOTE: patch_icpp alone does not reproduce the ' &
            & // 'initial condition, so the fine blocks start from data interpolated from level 0.'

    end subroutine s_amr_check_ic_patches

    !> Evaluate the initial condition patches on the installed grid frame into the interior of q_cons, on the host.
    impure subroutine s_amr_ic_on_frame(q_cons)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons
        type(scalar_field), dimension(sys_size)                :: q_prim
        type(ic_grid)                                          :: g
        integer                                                :: i

#ifdef MFC_MIXED_PRECISION
        integer(kind=1), allocatable, dimension(:,:,:) :: patch_id_fp
#else
        integer, allocatable, dimension(:,:,:) :: patch_id_fp
#endif

        g%m = m; g%n = n; g%p = p
        g%x_cc = x_cc; g%dx_min = minval(dx(0:m))
        if (n > 0) then
            g%y_cc = y_cc; g%dy_min = minval(dy(0:n))
        end if
        if (p > 0) then
            g%z_cc = z_cc; g%dz_min = minval(dz(0:p))
        end if
        allocate (patch_id_fp(0:m,0:n,0:p)); patch_id_fp = 0
        do i = 1, sys_size
            allocate (q_prim(i)%sf(0:m,0:n,0:p)); q_prim(i)%sf = 0._wp
        end do
        call s_apply_icpp_patches(patch_id_fp, q_prim, g)
        call s_convert_primitive_to_conservative_variables(q_prim, q_cons)
        do i = 1, sys_size
            deallocate (q_prim(i)%sf)
        end do

    end subroutine s_amr_ic_on_frame

end module m_amr_stage
