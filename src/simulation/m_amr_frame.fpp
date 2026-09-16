!>
!!@file
!!@brief Contains module m_amr_frame

#! AMD OpenMP lane: assert allocatables present on every kernel here (see OMP_DEFAULT_STR). Every conditionally allocated
#! module array a kernel here names launches only under its allocation's own condition (sw_jac/jac: igr;
#! amr_cg_pb/mv: do_pbmv; amr_prim_st/amr_bt_*: amr_prim_batch); amr_cg and amr_cons_br/stor_st are
#! allocated before first use. A kernel naming an unallocated array aborts. Keep it so.
#:set MFC_OMP_PRESENT_ALLOCATABLE = True
#:include 'macros.fpp'

!> @brief Block-frame swap (fine grid state in/out of the shared solver) and the pb/mv side-state services.
module m_amr_frame

#ifdef MFC_MPI
    use mpi  !< MPI-IO for the parallel_io AMR restart file
#endif

    use m_derived_types  ! scalar_field, t_box, int_bounds_info
    use m_global_parameters
    use m_mpi_proxy, only: s_mpi_abort  ! @:ASSERT expands to it
    use m_pressure_relaxation, only: s_pressure_relaxation_procedure
    use m_phase_timing
    use m_amr_xchg_audit  ! per-call-site accounting of every AMR p2p transfer (s_xa_rec + XA_* site ids)
    use m_hypoelastic, only: s_hypoelastic_update_fd_coeffs
    use m_active_box, only: ab_active
    use m_igr, only: jac, jac_old
    use m_amr_state
    use m_amr_distribution
    use m_amr_store
    use m_amr_exchange

    implicit none

    private
    public :: s_amr_pressure_relax_fine, s_amr_restore_coarse, s_amr_swap_to_fine

contains

    !> 6-equation model: apply the per-stage pressure relaxation to the fine block's interior (cell-local equilibration, no
    !! stencil), mirroring the coarse per-stage call. Swaps the grid so the routine's 0:m,0:n,0:p loop covers this block.
    impure subroutine s_amr_pressure_relax_fine()

        if (.not. amr_rank_owns_block) return
        call s_amr_swap_to_fine()
        call s_amr_br_load(amr_loc_of(amr_cur))
        call s_pressure_relaxation_procedure(amr_cons_br)
        call s_amr_br_store(amr_loc_of(amr_cur))
        call s_amr_restore_coarse()

    end subroutine s_amr_pressure_relax_fine

    !> Swap the global grid state to the fine block. Must be paired with s_amr_restore_coarse.
    impure subroutine s_amr_swap_to_fine()

        ! Saving on a nested swap would overwrite the sw_* bounce buffers with fine state, and the eventual restore would install
        ! fine extents as the coarse grid (silent corruption of everything after). Hence every save below is depth-guarded; the
        ! installs are not, since re-installing the same slot is idempotent.
        amr_swap_depth = amr_swap_depth + 1
        if (amr_swap_depth == 1) then
            sw_m = m; sw_n = n; sw_p = p
            sw_idwint = idwint; sw_idwbuff = idwbuff
        end if
        ! the acoustic source's precomputed spatials are coarse-grid cell indices: applying them on the fine block would inject at
        ! wrong cells (or out of bounds). The support is guaranteed not to overlap the block (checked at startup), so the fine RHS
        ! correctly skips the source.
        if (amr_swap_depth == 1) sw_acoustic_source = acoustic_source
        acoustic_source = .false.
        ! active-box windows are coarse cell indices: applying them on the swapped fine grid would window the wrong cells. Blocks
        ! are
        ! contained in the active window (init check + regrid clamp), so the fine advance legitimately treats its whole block as
        ! active.
        if (amr_swap_depth == 1) sw_ab_active = ab_active
        ab_active = .false.
        $:GPU_UPDATE(device='[ab_active]')
        m = amr_slots(amr_cur)%m; n = amr_slots(amr_cur)%n; p = amr_slots(amr_cur)%p
        idwint(1)%beg = 0; idwint(1)%end = m
        idwint(2)%beg = 0; idwint(2)%end = n
        idwint(3)%beg = 0; idwint(3)%end = p
        idwbuff = amr_slots(amr_cur)%idwbuff
        ! save coarse coords to bounce buffers, then copy fine coords into global arrays
        if (amr_swap_depth == 1) then
            sw_x_cb = x_cb; sw_x_cc = x_cc; sw_dx = dx
            if (n_glb > 0) then; sw_y_cb = y_cb; sw_y_cc = y_cc; sw_dy = dy; end if
            if (p_glb > 0) then; sw_z_cb = z_cb; sw_z_cc = z_cc; sw_dz = dz; end if
        end if
        x_cb(-1:amr_slots(amr_cur)%m) = amr_slots(amr_cur)%x_cb(-1:amr_slots(amr_cur)%m)
        x_cc(0:amr_slots(amr_cur)%m) = amr_slots(amr_cur)%x_cc(0:amr_slots(amr_cur)%m)
        dx(0:amr_slots(amr_cur)%m) = amr_slots(amr_cur)%dx(0:amr_slots(amr_cur)%m)
        if (n_glb > 0) then
            y_cb(-1:amr_slots(amr_cur)%n) = amr_slots(amr_cur)%y_cb(-1:amr_slots(amr_cur)%n)
            y_cc(0:amr_slots(amr_cur)%n) = amr_slots(amr_cur)%y_cc(0:amr_slots(amr_cur)%n)
            dy(0:amr_slots(amr_cur)%n) = amr_slots(amr_cur)%dy(0:amr_slots(amr_cur)%n)
        end if
        if (p_glb > 0) then
            z_cb(-1:amr_slots(amr_cur)%p) = amr_slots(amr_cur)%z_cb(-1:amr_slots(amr_cur)%p)
            z_cc(0:amr_slots(amr_cur)%p) = amr_slots(amr_cur)%z_cc(0:amr_slots(amr_cur)%p)
            dz(0:amr_slots(amr_cur)%p) = amr_slots(amr_cur)%dz(0:amr_slots(amr_cur)%p)
        end if
        ! Extend the fine grid into the ghost shell (s_build_level_coords only fills the interior 0:m). Ghost cells use the exact
        ! parent-cell bisection, the same formula as the interior, with floor division for negative indices. Fine-level
        ! distribution: the owner may not hold the block's coarse coordinate slice locally, so ghost parent boundaries come from the
        ! global boundaries amr_g?cb (cl is a global coarse index, region_lo + floor(jg/rr)), matching the interior build. Blocks
        ! stay buff_size inside the domain, so every ghost parent is an in-domain coarse cell with exact coords.
        block
            integer               :: jg, cl, pblk2, k, rr, pnf
            real(wp), allocatable :: cxb(:), cyb(:), czb(:), tcc(:), tdx(:)
            rr = amr_slots(amr_cur)%amr_ref_ratio
            ! ghost parent boundaries: a level>=2 block's coarse side is its parent's fine grid (indexed in the parent-fine
            ! amr_isect frame, matching the interior s_build_level_coords), not the L0 global boundaries. amr_isect_lo is a
            ! parent-fine index, so indexing amr_g?cb (sized for L0) would read out of bounds (garbage on host, NaN on the device
            ! copy). Source the parent's fine coords for level>=2, the global L0 boundaries for level 1.
            if (amr_block_level(amr_cur) >= 2) then
                ! Rebuild the parent's fine boundaries from replicated metadata; do not read amr_slots(pblk2)%x_cb. That array is
                ! allocated only on the parent's owner, and under per-level distribution this block's owner need not be it; taking
                ! lbound/ubound of an unallocated allocatable is undefined. Same ancestor replay as the interior build, so the
                ! ghost bisection and the interior agree exactly.
                pblk2 = f_amr_parent_block(amr_cur)
                pnf = amr_ref_ratio**amr_block_level(pblk2)*(amr_region_hi_all(1, pblk2) - amr_region_lo_all(1, pblk2) + 1) - 1
                allocate (cxb(-1:pnf), tcc(0:pnf), tdx(0:pnf))
                call s_amr_build_block_coords(pblk2, amr_gxcb, cxb, tcc, tdx, 1)
                deallocate (tcc, tdx)
                if (n_glb > 0) then
                    pnf = amr_ref_ratio**amr_block_level(pblk2)*(amr_region_hi_all(2, pblk2) - amr_region_lo_all(2, pblk2) + 1) - 1
                    allocate (cyb(-1:pnf), tcc(0:pnf), tdx(0:pnf))
                    call s_amr_build_block_coords(pblk2, amr_gycb, cyb, tcc, tdx, 2)
                    deallocate (tcc, tdx)
                end if
                if (p_glb > 0) then
                    pnf = amr_ref_ratio**amr_block_level(pblk2)*(amr_region_hi_all(3, pblk2) - amr_region_lo_all(3, pblk2) + 1) - 1
                    allocate (czb(-1:pnf), tcc(0:pnf), tdx(0:pnf))
                    call s_amr_build_block_coords(pblk2, amr_gzcb, czb, tcc, tdx, 3)
                    deallocate (tcc, tdx)
                end if
            else
                allocate (cxb(lbound(amr_gxcb, 1):ubound(amr_gxcb, 1))); cxb = amr_gxcb
                if (n_glb > 0) then; allocate (cyb(lbound(amr_gycb, 1):ubound(amr_gycb, 1))); cyb = amr_gycb; end if
                if (p_glb > 0) then; allocate (czb(lbound(amr_gzcb, 1):ubound(amr_gzcb, 1))); czb = amr_gzcb; end if
            end if
            do jg = amr_slots(amr_cur)%m + 1, amr_slots(amr_cur)%m + buff_size
                cl = amr_isect_lo(1) + floor(real(jg, wp)/real(rr, wp))
                k = modulo(jg, rr)
                if (k == rr - 1) then
                    x_cb(jg) = cxb(cl)
                else
                    x_cb(jg) = (real(rr - 1 - k, wp)*cxb(cl - 1) + real(k + 1, wp)*cxb(cl))/real(rr, wp)
                end if
                dx(jg) = x_cb(jg) - x_cb(jg - 1); x_cc(jg) = 0.5_wp*(x_cb(jg - 1) + x_cb(jg))
            end do
            ! unified boundary formula (matches the interior subdivision): boundary jg belongs to
            ! parent c = isect_lo + floor(jg/rr); sub-position k=modulo(jg,rr) picks the rr-way split
            do jg = -1 - buff_size, -1
                cl = amr_isect_lo(1) + floor(real(jg, wp)/real(rr, wp))
                k = modulo(jg, rr)
                if (k == rr - 1) then
                    x_cb(jg) = cxb(cl)
                else
                    x_cb(jg) = (real(rr - 1 - k, wp)*cxb(cl - 1) + real(k + 1, wp)*cxb(cl))/real(rr, wp)
                end if
            end do
            do jg = -buff_size, -1
                dx(jg) = x_cb(jg) - x_cb(jg - 1); x_cc(jg) = 0.5_wp*(x_cb(jg - 1) + x_cb(jg))
            end do
            if (n_glb > 0) then
                do jg = amr_slots(amr_cur)%n + 1, amr_slots(amr_cur)%n + buff_size
                    cl = amr_isect_lo(2) + floor(real(jg, wp)/real(rr, wp))
                    k = modulo(jg, rr)
                    if (k == rr - 1) then
                        y_cb(jg) = cyb(cl)
                    else
                        y_cb(jg) = (real(rr - 1 - k, wp)*cyb(cl - 1) + real(k + 1, wp)*cyb(cl))/real(rr, wp)
                    end if
                    dy(jg) = y_cb(jg) - y_cb(jg - 1); y_cc(jg) = 0.5_wp*(y_cb(jg - 1) + y_cb(jg))
                end do
                ! unified boundary formula (matches the interior subdivision): boundary jg belongs to
                ! parent c = isect_lo + floor(jg/rr); sub-position k=modulo(jg,rr) picks the rr-way split
                do jg = -1 - buff_size, -1
                    cl = amr_isect_lo(2) + floor(real(jg, wp)/real(rr, wp))
                    k = modulo(jg, rr)
                    if (k == rr - 1) then
                        y_cb(jg) = cyb(cl)
                    else
                        y_cb(jg) = (real(rr - 1 - k, wp)*cyb(cl - 1) + real(k + 1, wp)*cyb(cl))/real(rr, wp)
                    end if
                end do
                do jg = -buff_size, -1
                    dy(jg) = y_cb(jg) - y_cb(jg - 1); y_cc(jg) = 0.5_wp*(y_cb(jg - 1) + y_cb(jg))
                end do
            end if
            if (p_glb > 0) then
                do jg = amr_slots(amr_cur)%p + 1, amr_slots(amr_cur)%p + buff_size
                    cl = amr_isect_lo(3) + floor(real(jg, wp)/real(rr, wp))
                    k = modulo(jg, rr)
                    if (k == rr - 1) then
                        z_cb(jg) = czb(cl)
                    else
                        z_cb(jg) = (real(rr - 1 - k, wp)*czb(cl - 1) + real(k + 1, wp)*czb(cl))/real(rr, wp)
                    end if
                    dz(jg) = z_cb(jg) - z_cb(jg - 1); z_cc(jg) = 0.5_wp*(z_cb(jg - 1) + z_cb(jg))
                end do
                ! unified boundary formula (matches the interior subdivision): boundary jg belongs to
                ! parent c = isect_lo + floor(jg/rr); sub-position k=modulo(jg,rr) picks the rr-way split
                do jg = -1 - buff_size, -1
                    cl = amr_isect_lo(3) + floor(real(jg, wp)/real(rr, wp))
                    k = modulo(jg, rr)
                    if (k == rr - 1) then
                        z_cb(jg) = czb(cl)
                    else
                        z_cb(jg) = (real(rr - 1 - k, wp)*czb(cl - 1) + real(k + 1, wp)*czb(cl))/real(rr, wp)
                    end if
                end do
                do jg = -buff_size, -1
                    dz(jg) = z_cb(jg) - z_cb(jg - 1); z_cc(jg) = 0.5_wp*(z_cb(jg - 1) + z_cb(jg))
                end do
            end if
        end block
        ! batched advance: the leader's grid is installed above; extend it into the slab of amr_bat_n stacked blocks (stride
        ! amr_bat_w along amr_bat_sd), since the flux divergence reads dx/dy/dz at every slab cell. Cell boundaries (x_cb etc.)
        ! are not replicated: nothing on the batched path reads them (WENO coefficients are not recomputed on a uniform grid).
        if (amr_bat_n > 1) then
            block
                integer :: ibm, o, e
                e = amr_bat_ext(amr_bat_sd)
                do ibm = 2, amr_bat_n
                    o = (ibm - 1)*amr_bat_w
                    select case (amr_bat_sd)
                    case (1)
                        x_cc(o - buff_size:o + e + buff_size) = x_cc(-buff_size:e + buff_size)
                        dx(o - buff_size:o + e + buff_size) = dx(-buff_size:e + buff_size)
                    case (2)
                        y_cc(o - buff_size:o + e + buff_size) = y_cc(-buff_size:e + buff_size)
                        dy(o - buff_size:o + e + buff_size) = dy(-buff_size:e + buff_size)
                    case (3)
                        z_cc(o - buff_size:o + e + buff_size) = z_cc(-buff_size:e + buff_size)
                        dz(o - buff_size:o + e + buff_size) = dz(-buff_size:e + buff_size)
                    end select
                end do
                e = (amr_bat_n - 1)*amr_bat_w + e
                select case (amr_bat_sd)
                case (1); m = e
                case (2); n = e
                case (3); p = e
                end select
                idwint(amr_bat_sd)%end = e; idwbuff(amr_bat_sd)%end = e + buff_size
            end block
        end if
        ! sync the swapped extents/bounds/coordinates to the device: RHS kernels read the device copies of these GPU_DECLARE'd
        ! globals (stale coarse bounds = OOB kernels)
        call s_amr_sync_grid_state_to_device()
        ! hypoelastic stress sources use grid-spacing-dependent FD coefficients: recompute them from the (now fine) grid, else every
        ! fine velocity gradient is halved
        if (hypoelasticity) call s_hypoelastic_update_fd_coeffs()
        ! nonuniform coarse grid (stretched, or the axisymmetric axis half-cell): the per-cell WENO coefficients must be rebuilt for
        ! the block's own grid (no-op flag on uniform grids)

        ! IGR: save the coarse sigma state and seed the fine solve. jac holds this stage's converged coarse sigma (the coarse RHS
        ! ran
        ! first), so its parent values are both the best initial guess and the frozen Dirichlet ghost data for the block-local
        ! Jacobi
        ! solve (the per-iteration BC/halo populate is skipped under amr_in_fine_advance). Piecewise-constant parent injection over
        ! the full buffered fine range.
        if (igr) call s_amr_igr_swap_sigma()

    end subroutine s_amr_swap_to_fine

    !> Restore the global grid state saved by s_amr_swap_to_fine.
    impure subroutine s_amr_restore_coarse(sync_device)

        !> .false. only from the batched stage when another batch follows at once: its swap re-pushes the whole grid state before
        !! any kernel reads it, so the restore-side push is dead work there.
        logical, intent(in), optional :: sync_device

        @:ASSERT(amr_swap_depth > 0, "s_amr_restore_coarse without a matching s_amr_swap_to_fine")
        amr_swap_depth = amr_swap_depth - 1
        ! inner restore: the enclosing swap frame still wants its slot installed, and the sw_* buffers still hold the coarse state
        if (amr_swap_depth > 0) return
        m = sw_m; n = sw_n; p = sw_p
        idwint = sw_idwint; idwbuff = sw_idwbuff
        acoustic_source = sw_acoustic_source
        ab_active = sw_ab_active
        $:GPU_UPDATE(device='[ab_active]')
        ! restore full coarse coords from bounce buffers
        x_cb = sw_x_cb; x_cc = sw_x_cc; dx = sw_dx
        if (n_glb > 0) then; y_cb = sw_y_cb; y_cc = sw_y_cc; dy = sw_dy; end if
        if (p_glb > 0) then; z_cb = sw_z_cb; z_cc = sw_z_cc; dz = sw_dz; end if
        ! sync the restored coarse extents/bounds/coordinates back to the device
        if (.not. present(sync_device)) then
            call s_amr_sync_grid_state_to_device()
        else if (sync_device) then
            call s_amr_sync_grid_state_to_device()
        end if
        if (hypoelasticity) call s_hypoelastic_update_fd_coeffs()
        if (igr) call s_amr_igr_restore_sigma()

    end subroutine s_amr_restore_coarse

    !> Save the coarse jac/jac_old and seed the (already swapped-in) fine block's sigma state by piecewise-constant parent injection
    !! from the saved coarse sigma: interior = initial guess, ghost shell = frozen Dirichlet coupling data for the block-local
    !! solve.
    impure subroutine s_amr_igr_swap_sigma()

        integer :: j, k, l, ci, cj, ck, ibm, g, o1, o2, o3, mb1, me1, mb2, me2, mb3, me3
        integer :: cb1, ce1, cb2, ce2, cb3, ce3, fb1, fe1, fb2, fe2, fb3, fe3
        integer :: lo1, lo2, lo3, ox, oy, oz

        ! bounds/offsets hoisted to scalars: sw_idwbuff (and friends) are host-only module state; referencing them inside the
        ! kernels makes OpenACC's present lookup fail (OpenMP's implicit map(to) tolerates it, so only OpenACC builds crash)

        cb1 = sw_idwbuff(1)%beg; ce1 = sw_idwbuff(1)%end
        cb2 = sw_idwbuff(2)%beg; ce2 = sw_idwbuff(2)%end
        cb3 = sw_idwbuff(3)%beg; ce3 = sw_idwbuff(3)%end
        fb1 = idwbuff(1)%beg; fe1 = idwbuff(1)%end
        fb2 = idwbuff(2)%beg; fe2 = idwbuff(2)%end
        fb3 = idwbuff(3)%beg; fe3 = idwbuff(3)%end
        lo1 = amr_isect_lo(1); lo2 = amr_isect_lo(2); lo3 = amr_isect_lo(3)
        ox = start_idx(1); oy = 0; oz = 0
        if (n_glb > 0) oy = start_idx(2)
        if (p_glb > 0) oz = start_idx(3)
        ! Save the coarse sigma, outermost swap only. A nested swap must not re-save, or sw_jac would take fine state and both the
        ! seed below and s_amr_igr_restore_sigma would work from it. The seed that follows is not guarded: it reads sw_jac, which
        ! still holds the coarse state, so every nested block seeds from the correct parent.
        if (amr_swap_depth == 1) then
            $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
            do l = cb3, ce3
                do k = cb2, ce2
                    do j = cb1, ce1
                        sw_jac(j, k, l) = jac(j, k, l)
                        sw_jac_old(j, k, l) = jac_old(j, k, l)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if
        if (amr_bat_n > 1) then
            ! batched slab: each member's buffered range is seeded from its own parent; the slab bounds above are the leader's
            do ibm = 1, amr_bat_n
                g = amr_bat_blk(ibm)
                lo1 = amr_isect_lo_all(1, g); lo2 = amr_isect_lo_all(2, g); lo3 = amr_isect_lo_all(3, g)
                o1 = 0; o2 = 0; o3 = 0
                select case (amr_bat_sd)
                case (1); o1 = (ibm - 1)*amr_bat_w
                case (2); o2 = (ibm - 1)*amr_bat_w
                case default; o3 = (ibm - 1)*amr_bat_w
                end select
                mb1 = -buff_size; me1 = amr_bat_mext(1, ibm) + buff_size
                mb2 = 0; me2 = 0; mb3 = 0; me3 = 0
                if (n_glb > 0) then; mb2 = -buff_size; me2 = amr_bat_mext(2, ibm) + buff_size; end if
                if (p_glb > 0) then; mb3 = -buff_size; me3 = amr_bat_mext(3, ibm) + buff_size; end if
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, ci, cj, ck]', copyin='[lo1, lo2, lo3, o1, o2, o3, mb1, me1, &
                                    & mb2, me2, mb3, me3]')
                do l = mb3, me3
                    do k = mb2, me2
                        do j = mb1, me1
                            ci = lo1 + floor(real(j, wp)/real(amr_ref_ratio, wp)) - ox
                            cj = 0; ck = 0
                            if (n_glb > 0) cj = lo2 + floor(real(k, wp)/real(amr_ref_ratio, wp)) - oy
                            if (p_glb > 0) ck = lo3 + floor(real(l, wp)/real(amr_ref_ratio, wp)) - oz
                            ci = min(max(ci, cb1), ce1)
                            cj = min(max(cj, cb2), ce2)
                            ck = min(max(ck, cb3), ce3)
                            jac(j + o1, k + o2, l + o3) = sw_jac(ci, cj, ck)
                            jac_old(j + o1, k + o2, l + o3) = sw_jac(ci, cj, ck)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end do
            return
        end if
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, ci, cj, ck]')
        do l = fb3, fe3
            do k = fb2, fe2
                do j = fb1, fe1
                    ci = lo1 + floor(real(j, wp)/real(amr_ref_ratio, wp)) - ox
                    cj = 0; ck = 0
                    if (n_glb > 0) cj = lo2 + floor(real(k, wp)/real(amr_ref_ratio, wp)) - oy
                    if (p_glb > 0) ck = lo3 + floor(real(l, wp)/real(amr_ref_ratio, wp)) - oz
                    ci = min(max(ci, cb1), ce1)
                    cj = min(max(cj, cb2), ce2)
                    ck = min(max(ck, cb3), ce3)
                    jac(j, k, l) = sw_jac(ci, cj, ck)
                    jac_old(j, k, l) = sw_jac(ci, cj, ck)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_igr_swap_sigma

    !> Restore the coarse jac/jac_old saved by s_amr_igr_swap_sigma (bounds already restored).
    impure subroutine s_amr_igr_restore_sigma()

        integer :: j, k, l, b1, e1, b2, e2, b3, e3

        b1 = idwbuff(1)%beg; e1 = idwbuff(1)%end
        b2 = idwbuff(2)%beg; e2 = idwbuff(2)%end
        b3 = idwbuff(3)%beg; e3 = idwbuff(3)%end
        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
        do l = b3, e3
            do k = b2, e2
                do j = b1, e1
                    jac(j, k, l) = sw_jac(j, k, l)
                    jac_old(j, k, l) = sw_jac_old(j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_igr_restore_sigma

end module m_amr_frame
