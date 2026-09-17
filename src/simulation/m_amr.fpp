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
        amr_br_batch = amr_bat_max

        ! Fine-block cap = the case amr_max_blocks; the shared pool adds the L0 tile prefix (l0_slot_off, 0 when l0_ntile=0)
        ! ahead of it, so both AMR fine blocks and any L0 tiles draw from one amr_slots allocation.
        amr_max_fine = amr_max_blocks
        amr_max_blocks = l0_slot_off + amr_max_fine
        allocate (amr_slots(1:amr_max_blocks))
        call s_amr_loc_index_init()
        allocate (amr_region_lo_all(3, amr_max_blocks), amr_region_hi_all(3, amr_max_blocks))
        allocate (amr_isect_lo_all(3, amr_max_blocks), amr_isect_hi_all(3, amr_max_blocks))
        allocate (amr_owns_all(amr_max_blocks), amr_block_owner(amr_max_blocks), amr_block_level(amr_max_blocks))
        allocate (amr_owner_cut(0:num_procs - 1)); amr_owner_cut = -1_8
        allocate (amr_fine_cut(0:num_procs - 1,1:max(amr_max_level, 1))); amr_fine_cut = -1_8
        ! amr_ovl_gather/scatter (the 2D rank lists) are allocated to the computed max overlap in s_amr_build_seam_pairs; only the
        ! per-block counts are sized here.
        allocate (amr_ovl_gather_n(amr_max_blocks), amr_ovl_scatter_n(amr_max_blocks))
        allocate (amr_slot_live(amr_max_blocks)); amr_slot_live = .false.
        amr_region_lo_all = 0; amr_region_hi_all = 0; amr_isect_lo_all = 0; amr_isect_hi_all = 0; amr_owns_all = .false.
        amr_block_owner = 0
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
        call s_amr_build_global_cb()  ! the fine-distribution owner rebuilds whole-block fine coordinates from these
        call s_amr_init_coarse_patch()
        ! the coarse decomposition (each rank's coarse start_idx + local m/n/p) is a structured cartesian split, computed O(1) per
        ! rank by s_amr_rank_decomp - no replicated table, no allgather. Validate the formula against this rank's actual values.
        call s_amr_validate_decomp()
        ! per-slot fine-grid IB marker fields (static-body AMR); sized to the same max buffered fine extents as q_cons so the fine
        ! IB pipeline can resolve the body on the block
        if (ib) call s_ibm_alloc_fine(amr_max_blocks, mbuf1_lo, mbuf1_hi, mbuf2_lo, mbuf2_hi, mbuf3_lo, mbuf3_hi)
        call s_amr_init_first_blocks()
        call s_amr_init_tags()

    end subroutine s_initialize_amr_module

    !> The block caps and the preallocation extents. Mirror decomposition: each rank holds the fine cells covering block /\ its own
    !! subdomain (np=1: the whole block). buff_size is not available at checker time, so the geometric aborts live here.
    impure subroutine s_amr_init_extents()

        integer :: d, ext(3), bad_loc, bad_glb, fit_d

        ext = 0
        ext(1) = m
        if (n_glb > 0) ext(2) = n
        if (p_glb > 0) ext(3) = p
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
        amr_maxc(1) = (m_glb + 1)/amr_ref_ratio
        amr_maxc(2) = 1; amr_maxc(3) = 1
        if (n_glb > 0) amr_maxc(2) = (n_glb + 1)/amr_ref_ratio
        if (p_glb > 0) amr_maxc(3) = (p_glb + 1)/amr_ref_ratio

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
            call s_mpi_allreduce_integer_min((ext(d) + 1)/amr_ref_ratio, fit_d)
            if (amr_max_grid_size > 0) then
                amr_maxc_fit(d) = min(amr_maxc(d), amr_max_grid_size)
            else
                amr_maxc_fit(d) = min(amr_maxc(d), fit_d)
            end if
        end do

        ! max fine extents and buffered bounds for preallocation
        max_f1 = amr_ref_ratio*amr_maxc_fit(1) - 1
        max_f2 = 0; max_f3 = 0
        if (n_glb > 0) max_f2 = amr_ref_ratio*amr_maxc_fit(2) - 1
        if (p_glb > 0) max_f3 = amr_ref_ratio*amr_maxc_fit(3) - 1
        mbuf1_lo = -buff_size; mbuf1_hi = max_f1 + buff_size
        mbuf2_lo = 0; mbuf2_hi = 0; mbuf3_lo = 0; mbuf3_hi = 0
        if (n_glb > 0) then; mbuf2_lo = -buff_size; mbuf2_hi = max_f2 + buff_size; end if
        if (p_glb > 0) then; mbuf3_lo = -buff_size; mbuf3_hi = max_f3 + buff_size; end if

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
        cells = real(mbuf1_hi - mbuf1_lo + 1, wp)
        if (n_glb > 0) cells = cells*real(mbuf2_hi - mbuf2_lo + 1, wp)
        if (p_glb > 0) cells = cells*real(mbuf3_hi - mbuf3_lo + 1, wp)
        nfam = 2._wp
        slot_gib = cells*real(sys_size, wp)*nfam*real(storage_size(1._wp)/8, wp)/1024._wp**3
        print '(A,I0,A,I0,A,ES10.3,A,F8.3,A)', ' [amr] per-block slot: ', nint(cells), ' cells x sys_size x ', nint(nfam), &
            & ' fields = ', cells*real(sys_size, wp)*nfam, ' words (', slot_gib, ' GiB per owned block)'
        print '(A,F9.2,A,I0,A)', ' [amr]   worst case if one rank owned every block: ', slot_gib*real(amr_max_blocks, wp), &
            & ' GiB (amr_max_blocks = ', amr_max_blocks, '). Typical is amr_max_blocks/num_procs blocks per rank.'

    end subroutine s_amr_init_report

    !> Bounce buffers for the copy-based coordinate swap (GPU-safe; same bounds as the base-level global arrays, which are sized on
    !! *_alloc - these are whole-array assigned to/from x_cb etc., so the shapes must agree).
    impure subroutine s_amr_init_swap_buffers()

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

    end subroutine s_amr_init_swap_buffers

    !> The coarse-patch gather buffer (see amr_cg's declaration): sized to the largest block's coarse footprint (block coarse cells
    !! + 2*amr_cpat_mar halo, block-local frame), device-mapped so the runtime ghost fill reads it on the owner.
    impure subroutine s_amr_init_coarse_patch()

        type(scalar_field), allocatable :: tmp_cg(:)
        integer                         :: i

        amr_cpat_mar = (buff_size + amr_ref_ratio - 1)/amr_ref_ratio + 1
        amr_cpat_hi = 0
        amr_cpat_hi(1) = amr_maxc_fit(1) - 1 + 2*amr_cpat_mar
        if (n_glb > 0) amr_cpat_hi(2) = amr_maxc_fit(2) - 1 + 2*amr_cpat_mar
        if (p_glb > 0) amr_cpat_hi(3) = amr_maxc_fit(3) - 1 + 2*amr_cpat_mar
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

    !> Per-family tag bases sit above the per-box tag space so the two cannot collide; the keyed band space starts at the next 65536
    !! boundary above every per-box tag (bases + their mod-100 folds).
    impure subroutine s_amr_init_tags()

        integer :: f

#ifdef MFC_MPI
        integer(kind=MPI_ADDRESS_KIND) :: tag_ub
        logical                        :: tag_ub_set
        integer                        :: ierr
#endif

        do f = 1, size(amr_tag_base)
            amr_tag_base(f) = amr_max_blocks + 100*f
        end do
        amr_m1_base = ((amr_tag_base(size(amr_tag_base)) + 100)/65536 + 1)*65536
#ifdef MFC_MPI
        call MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, tag_ub, tag_ub_set, ierr)
        @:ASSERT(tag_ub_set, "MPI_TAG_UB attribute unavailable")
        @:ASSERT(amr_tag_base(size(amr_tag_base)) + 100 <= tag_ub, &
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
        deallocate (amr_slot_live)
        call s_amr_st_finalize()
        if (allocated(amr_seam_pairs)) deallocate (amr_seam_pairs)
        if (allocated(amr_ovl_gather)) deallocate (amr_ovl_gather)
        if (allocated(amr_ovl_scatter)) deallocate (amr_ovl_scatter)
        deallocate (amr_ovl_gather_n, amr_ovl_scatter_n)
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

    end subroutine s_finalize_amr_module

end module m_amr
