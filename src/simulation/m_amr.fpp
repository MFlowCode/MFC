!>
!!@file
!!@brief Contains module m_amr

#! AMD OpenMP lane: assert allocatables present on every kernel here (see OMP_DEFAULT_STR). A conditionally allocated module
#! array a kernel names launches only under its allocation's own condition (sw_jac/jac: igr); a kernel naming an unallocated
#! array aborts. Keep it so.
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
    use m_global_parameters
    use m_mpi_proxy, only: s_mpi_abort
    use m_mpi_common, only: s_mpi_allreduce_integer_min, s_mpi_allreduce_integer_max
    use m_phase_timing
    use m_amr_xchg_audit  ! per-call-site accounting of every AMR p2p transfer (s_xa_rec + XA_* site ids)
    use m_ibm, only: s_ibm_alloc_fine
    use m_amr_state
    use m_amr_distribution
    use m_amr_wave
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

        ! shared-pool layout: tiles are a fixed level-0 prefix; fine blocks follow. Both this init and s_l0_tiles_init read this, so
        ! it runs before the amr early-return below (this routine always executes first, per m_start_up.fpp).
        amr_dim = [.true., n_glb > 0, p_glb > 0]
        amr_sidx = 0; amr_sidx(1:num_dims) = start_idx
        amr_ext = merge([m, n, p], 0, amr_dim)
        if (l0_ntile > 0) then
            l0_nt = merge(l0_ntile, 1, amr_dim)
            l0_ntiles_tot = num_procs*l0_nt(1)*l0_nt(2)*l0_nt(3)
            l0_slot_off = l0_ntiles_tot
        end if

        if (.not. amr) return

#ifdef MFC_GPU
        amr_fw_dev = rdma_mpi .and. XA_NH == 0
#endif
        amr_br_batch = amr_bat_max

        ! Fine-block cap = the case amr_max_blocks; the shared pool adds the L0 tile prefix (l0_slot_off, 0 when l0_ntile=0)
        ! ahead of it, so both AMR fine blocks and any L0 tiles draw from one amr_slots allocation.
        amr_max_fine = amr_max_blocks
        amr_max_blocks = l0_slot_off + amr_max_fine
        call s_amr_alloc_pool()
        amr_block_level = 1  ! init default (level-1); regrid re-tags each block's level for nesting
        amr_num_levels = 1
        amr_num_blocks = f_l0_slot(1)
        amr_cur = f_l0_slot(1)
        amr_seam_pairs_dirty = .true.; amr_seam_pairs_nblk = -1  ! force a seam-list build on the first fine-fine halo
        amr_mesh_epoch = amr_mesh_epoch + 1

        call s_amr_init_extents()
        call s_amr_init_report()
        ! with tiles, s_l0_tiles_init's mbuf union may still enlarge the extents; the scratch waits for it (see s_amr_scr_init)
        if (l0_ntile == 0) call s_amr_scr_init()
        call s_amr_init_swap_buffers()
        call s_amr_build_global_cb(0)  ! the fine-distribution owner rebuilds whole-block fine coordinates from these
        call s_amr_init_coarse_patch()
        ! the coarse decomposition (each rank's coarse start_idx + local m/n/p) is a structured cartesian split, computed O(1) per
        ! rank by s_amr_rank_decomp - no replicated table, no allgather. Validate the formula against this rank's actual values.
        call s_amr_validate_decomp()
        ! per-slot fine-grid IB marker fields (static-body AMR); sized to the same max buffered fine extents as q_cons so the fine
        ! IB pipeline can resolve the body on the block
        if (ib) call s_ibm_alloc_fine(amr_max_blocks, mbuf_lo(1), mbuf_hi(1), mbuf_lo(2), mbuf_hi(2), mbuf_lo(3), mbuf_hi(3))
        call s_amr_init_first_blocks()
        call s_amr_init_tags()

    end subroutine s_initialize_amr_module

    !> The block caps and the preallocation extents. Mirror decomposition: each rank holds the fine cells covering block /\ its own
    !! subdomain (np=1: the whole block). buff_size is not available at checker time, so the geometric aborts live here.
    impure subroutine s_amr_init_extents()

        integer :: d, bad_loc, bad_glb, fit_d

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
        ! block's fine extent (2*block-1) must fit every rank's local extent. non-IB: the block is tiled into <= amr_maxc_fit
        ! sub-blocks (each fits every rank's scratch), so no cap is needed. IB keeps a single contiguous block per body, so an IB
        ! block must itself fit a rank's local half-extent. Checked on the replicated block box so all ranks agree.
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
        amr_maxc = merge(([m_glb, n_glb, p_glb] + 1)/amr_ref_ratio, 1, amr_dim)

        ! regrid size cap. Default (amr_max_grid_size == 0): min over ranks of the local half-extent (= amr_maxc at np=1), so any
        ! clamped box satisfies every rank's scratch constraint and can move freely across ranks. That cap shrinks as ranks are
        ! added, which tiles a fixed feature into more and more blocks the further you scale (per-block cost is roughly fixed
        ! regardless of block size, so the block count is what costs), and it makes the box set (and so the answer, within
        ! tolerance) depend on the rank count. amr_max_grid_size > 0 pins the cap to an absolute number of coarse cells instead,
        ! like AMReX's max_grid_size: the box set is then identical at every rank count, and the solver scratch is sized to the
        ! cap rather than to the subdomain (idwbuff_alloc and m/n/p_alloc in m_global_parameters), so a block at the cap fits
        ! however small the subdomain becomes. A block is owned whole, so amr_maxc_fit (every regrid box is clamped to it) is
        ! the largest block any rank can own and sizes the fine/coord arrays.
        amr_maxc_fit = amr_maxc
        do d = 1, num_dims
            call s_mpi_allreduce_integer_min((amr_ext(d) + 1)/amr_ref_ratio, fit_d)
            if (amr_max_grid_size > 0) then
                amr_maxc_fit(d) = min(amr_maxc(d), amr_max_grid_size)
            else
                amr_maxc_fit(d) = min(amr_maxc(d), fit_d)
            end if
        end do

        ! max fine extents and buffered bounds for preallocation
        max_f = merge(amr_ref_ratio*amr_maxc_fit - 1, 0, amr_dim)
        call s_amr_set_mbuf()

    end subroutine s_amr_init_extents

    !> Rank-0 advisories (advice, not constraints: every setting is legal and sometimes correct) and the memory demand. Collective:
    !! the uniform-spacing check reduces over ranks.
    impure subroutine s_amr_init_report()

        integer  :: nonuni, nonuni_glb
        real(wp) :: slot_gib, cells, nfam

        ! stacked members share the batch leader's coordinate arrays in the non-stacked dimensions and read the coarse WENO
        ! coefficients at their stacked index, which is exact only where the grid spacing is bitwise uniform. The validator
        ! forbids stretched grids under amr; say so once when the spacing still differs at roundoff.

        nonuni = 0
        if (any(dx(0:m) /= dx(0))) nonuni = 1
        if (n_glb > 0) then; if (any(dy(0:n) /= dy(0))) nonuni = 1; end if
        if (p_glb > 0) then; if (any(dz(0:p) /= dz(0))) nonuni = 1; end if
        call s_mpi_allreduce_integer_max(nonuni, nonuni_glb)
        if (proc_rank /= 0) return

        if (nonuni_glb == 1) print '(A)', &
            & ' [amr] NOTE: the grid''s cell spacing is not bitwise uniform: stacked blocks reuse the batch leader''s ' &
            & // 'coordinate arrays, so members after the first see the leader''s spacing'
        ! the SFC map spreads whole blocks, so with fewer blocks than ranks some ranks own no fine work
        if (num_procs > amr_max_fine) print '(A,I0,A,I0,A)', ' [amr] WARNING: amr_max_blocks (', amr_max_blocks, &
            & ') < num_procs (', num_procs, &
            & '): the fine level can occupy at most amr_max_blocks ranks - raise amr_max_blocks for better fine-level balance'
        ! Every fine block advances at the coarse dt, but a level-l cell is amr_ref_ratio**l smaller, so its CFL limit is
        ! amr_ref_ratio**amr_max_level tighter than the coarse grid's; the dt (fixed, or the coarse-only cfl_dt estimate) is
        ! not scaled for that. The true CFL is unknown at init, so warn rather than abort - a small enough dt is valid.
        if (amr_ref_ratio > 2 .or. amr_max_level > 1) print '(A,I0,A)', &
            & ' [amr] WARNING: fine blocks advance at the coarse dt, but the ' &
            & // 'finest cell is amr_ref_ratio**amr_max_level = ', amr_ref_ratio**amr_max_level, &
            & 'x smaller - ensure dt satisfies the FINEST cell CFL (roughly the coarse-stable dt divided by that ' &
            & // 'factor), else the fine block may go unstable'
        ! amr_regrid_int = 0 is static AMR (the only mode above amr_max_level = 2), but a user who set amr = T expecting
        ! adaptivity gets none, silently
        if (amr_regrid_int == 0) print '(A)', &
            & ' [amr] NOTE: amr_regrid_int = 0 - the block set is STATIC and never adapts. ' &
            & // 'Set amr_regrid_int > 0 (4-8 is a reasonable start) for adaptive refinement.'
        if (amr_max_grid_size == 0 .and. num_procs > 1) print '(A)', &
            & ' [amr] NOTE: amr_max_grid_size = 0 derives the block cap from the ' &
            & // 'decomposition, so it SHRINKS as ranks are added and the box set depends on rank ' &
            & // 'count. Pinning it (64 measured best in 3D on MI250X, memory-bounded) was 3.0x ' &
            & // 'faster and makes the box set rank-invariant.'
        if (amr_regrid_int > 0 .and. amr_regrid_int < 4) print '(A,I0,A)', ' [amr] NOTE: amr_regrid_int = ', amr_regrid_int, &
            & ' regrids often; the tag sweep is per-CELL and flat in box count, so interval 8 ' &
            & // 'measured 1.39x faster. Raise it unless the refined feature moves quickly.'

        ! Memory demand, reported rather than guessed: there is no portable way to ask how much device memory is available
        ! across four compilers and three offload backends, so no cap is derived from a budget. A block costs 2 per-slot field
        ! families (q_cons, q_cons_stor; q_prim/rhs are one pooled scratch pair) x sys_size arrays on the mbuf extents. Slot volume
        ! goes as cap**num_dims, so one cap cannot serve 2D and 3D alike. Exceeding device memory aborts inside
        ! __tgt_target_data_begin_mapper, which presents as a hang (one rank dies, the rest block in MPI).
        cells = real(product(mbuf_hi - mbuf_lo + 1), wp)
        nfam = 2._wp
        slot_gib = cells*real(sys_size, wp)*nfam*real(storage_size(1._wp)/8, wp)/1024._wp**3
        print '(A,I0,A,I0,A,ES10.3,A,F8.3,A)', ' [amr] per-block slot: ', nint(cells), ' cells x sys_size x ', nint(nfam), &
            & ' fields = ', cells*real(sys_size, wp)*nfam, ' words (', slot_gib, ' GiB per owned block)'
        print '(A,F9.2,A,I0,A)', ' [amr]   worst case if one rank owned every block: ', slot_gib*real(amr_max_blocks, wp), &
            & ' GiB (amr_max_blocks = ', amr_max_blocks, '). Typical is amr_max_blocks/num_procs blocks per rank.'

    end subroutine s_amr_init_report

    !> The coarse-patch gather buffer (see amr_cg's declaration): sized to the largest block's coarse footprint (block coarse cells
    !! + 2*amr_cpat_mar halo, block-local frame), device-mapped so the runtime ghost fill reads it on the owner.
    impure subroutine s_amr_init_coarse_patch()

        type(scalar_field), allocatable :: tmp_cg(:)
        integer                         :: i

        amr_cpat_mar = (buff_size + amr_ref_ratio - 1)/amr_ref_ratio + 1
        amr_cpat_hi = merge(amr_maxc_fit - 1 + 2*amr_cpat_mar, 0, amr_dim)
        ! CCE OpenMP-offload leaves a bare module-scope derived-type (scalar_field) allocatable's descriptor uninitialized, so a
        ! direct allocate(amr_cg(1:sys_size)) aborts with lib-4425 at program start (a local scalar_field array and a
        ! GPU_DECLARE'd module one like q_prim_vf both allocate fine; only a bare module array does not). Allocate a local, which
        ! gets a valid descriptor, and hand it to the module variable via move_alloc, then map.
        allocate (tmp_cg(1:sys_size))
        @:ALLOCATE(amr_slab_tab(1:8, 1:6))
        call move_alloc(tmp_cg, amr_cg)
        $:GPU_ENTER_DATA(create='[amr_cg]')
        do i = 1, sys_size
            @:ALLOCATE(amr_cg(i)%sf(0:amr_cpat_hi(1), 0:amr_cpat_hi(2), 0:amr_cpat_hi(3)))
            amr_cg(i)%sf = 0._stp  ! padding beyond a block's valid patch extent is never read; keep it finite for the device copy
            @:ACC_SETUP_SFs(amr_cg(i))
        end do

    end subroutine s_amr_init_coarse_patch

    !> Place the initial block(s) and set their geometry (region, m/n/p, idwbuff, coordinates). Under dynamic regrid with bodies the
    !! initial block gets the same body-containment expansion regrid boxes get (the moving-body containment guard requires it from
    !! step 1); for a static block (amr_regrid_int = 0) the user's placement is authoritative. max_grid_size tiling: the initial
    !! block splits into <= amr_maxc_fit sub-blocks (at np=1 amr_maxc_fit == amr_maxc, so a normal block stays a single tile), one
    !! per slot; IB keeps a single contiguous block. Per-slot field arrays are allocated lazily by s_amr_reconcile_slots once
    !! ownership is known, so a rank holds only its owned blocks' fine arrays.
    impure subroutine s_amr_init_first_blocks()

        type(t_box), allocatable :: tiled(:)
        integer                  :: blk_lo(3), blk_hi(3), nt, capt, kk

        blk_lo = amr_block_beg; blk_hi = amr_block_end
        if (ib .and. amr_regrid_int > 0) call s_amr_expand_box_over_bodies(blk_lo, blk_hi)
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
        call s_amr_assign_block_owners()
        call s_amr_reconcile_slots()  ! allocate this rank's owned initial blocks (owner-guarded geometry writes below)
        do kk = 1, nt
            amr_cur = f_l0_slot(kk)
            call s_set_amr_fine_geometry(tiled(kk)%lo, tiled(kk)%hi)
        end do
        call s_amr_reduce_xchg_flag()
        call s_amr_select_slot(f_l0_slot(1))  ! refresh the per-block mirrors (geometry loop left them on the last tile)
        deallocate (tiled)

    end subroutine s_amr_init_first_blocks

    !> The narrow-reduction tag base sits above the per-box tag space so the two cannot collide; the keyed band space starts at the
    !! next 65536 boundary above every per-box tag.
    impure subroutine s_amr_init_tags()

#ifdef MFC_MPI
        integer(kind=MPI_ADDRESS_KIND) :: tag_ub
        logical                        :: tag_ub_set
        integer                        :: ierr
#endif

        amr_tag_narrow = amr_max_blocks + 400
        amr_m1_base = ((amr_max_blocks + 800)/65536 + 1)*65536
#ifdef MFC_MPI
        call MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, tag_ub, tag_ub_set, ierr)
        @:ASSERT(tag_ub_set, "MPI_TAG_UB attribute unavailable")
        @:ASSERT(amr_max_blocks + 800 <= tag_ub, &
                 & "AMR tag space exceeds MPI_TAG_UB: amr_max_blocks is too large for this MPI's tag range")
        @:ASSERT(amr_m1_base + 8*65536 <= tag_ub, "AMR keyed-tag band space exceeds MPI_TAG_UB")
#endif

    end subroutine s_amr_init_tags

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
        if (.not. amr) return
        do islot = 1, amr_max_blocks
            call s_amr_free_slot(islot)
        end do
        call s_amr_free_pool()
        ! gfortran/ifx abort on deallocating an unallocated array (amdflang silently tolerates it): guard each
        if (allocated(amr_fw_rblk)) deallocate (amr_fw_rblk)
        if (allocated(amr_sw_sq)) deallocate (amr_sw_sq, amr_sw_rq)
        #:for A in ['amr_fw_sq', 'amr_fw_rq']
            if (allocated(${A}$)) then
                if (amr_fw_dev) then
                    $:GPU_EXIT_DATA(delete='[' + A + ']')
                end if
                deallocate (${A}$)
            end if
        #:endfor
        #:for A in ['amr_my_blk', 'amr_l1r_blk', 'amr_l1p_blk', 'amr_fch_blk', 'amr_own_blk', 'amr_parent_blk', &
            'amr_child_ptr', 'amr_child_idx']
            if (allocated(${A}$)) deallocate (${A}$)
        #:endfor
        do i = 1, sys_size
            @:DEALLOCATE(amr_cg(i)%sf)
        end do
        @:DEALLOCATE(amr_cg)
        @:DEALLOCATE(amr_slab_tab)
        call s_amr_free_swap_buffers()

    end subroutine s_finalize_amr_module

end module m_amr
