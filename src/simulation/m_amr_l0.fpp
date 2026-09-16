!>
!!@file
!!@brief Contains module m_amr_l0

#! AMD OpenMP lane: assert allocatables present on every kernel here (see OMP_DEFAULT_STR). Every conditionally allocated
#! module array a kernel here names launches only under its allocation's own condition (amr_rvw: cyl_coord; sw_jac/jac: igr;
#! amr_cg_pb/mv: do_pbmv; amr_gst_a/b: amr_subcycle; amr_prim_st/amr_bt_*: amr_prim_batch); amr_cg and amr_cons_br/stor_st are
#! allocated before first use. A kernel naming an unallocated array aborts. Keep it so.
#:set MFC_OMP_PRESENT_ALLOCATABLE = True
#:include 'macros.fpp'

!> @brief Level-0 tiling: tile slots, coarse<->tile copies, edge BCs and the tile-stage advance.
module m_amr_l0

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
    use m_amr_advance

    implicit none

    private
    public :: s_l0_add_reflux_to_tiles, s_l0_advance_stage, s_l0_advance_stage_rhs, s_l0_advance_stage_rk, &
        & s_l0_copy_coarse_to_tiles, s_l0_fill_tiles_from_coarse, s_l0_forced_remap, s_l0_rebalance, s_l0_restrict_to_tiles, &
        & s_l0_scatter_tiles_to_coarse, s_l0_tiles_finalize, s_l0_tiles_init

contains

    ! Level-0 tiling (l0_ntile > 0):
    ! Base grid tiled into l0_ntiles_tot refinement-ratio-1 blocks advanced through the shared swap-based per-block solver; the
    ! bit-identity oracle is l0_ntile=0 (monolithic). See the l0_ntiles_tot declaration above for the design.

    !> Low global-cell index of tile `it` (0-based) when `ncell` cells split into `nt` balanced tiles: the first mod(ncell,nt) tiles
    !! get one extra cell. Tile `it` spans [f_l0_lo(it) : f_l0_lo(it+1)-1]; the widest tile is ceil(ncell/nt) cells.
    pure integer function f_l0_lo(ncell, nt, it)

        integer, intent(in) :: ncell, nt, it

        f_l0_lo = it*(ncell/nt) + min(it, mod(ncell, nt))

    end function f_l0_lo

    !> True if a domain-face BC code is a physical boundary the tiling does not support. Supported: extrapolation (BC_GHOST_EXTRAP),
    !! reflective (BC_REFLECTIVE), and periodic (BC_PERIODIC); each has a cons-space tile fill matching the monolithic prim-space
    !! BC. The characteristic / slip-wall / dirichlet family (bc < BC_GHOST_EXTRAP) has no tile fill. MPI processor boundaries (bc
    !! >= 0, incl. a periodic wrap-neighbour rank at np>1) are interior seams handled by the fine-fine halo, never a physical face
    !! here.
    pure logical function f_l0_bc_unsupported(bc)

        integer, intent(in) :: bc

        f_l0_bc_unsupported = (bc < BC_GHOST_EXTRAP)

    end function f_l0_bc_unsupported

    !> Build the base-grid tiling: allocate the slot/region/seam machinery for l0_ntiles_tot rr=1 tiles covering L0, set each tile's
    !! geometry (extent, L0-slice coords, whole-tile footprint, owner) and copy the current L0 state in. Standalone (does not call
    !! s_initialize_amr_module).
    impure subroutine s_l0_tiles_init()

        integer :: nt(3), ix, iy, iz, k, r, e
        integer :: tlo(3), thi(3)
        integer :: rsidx(3), rext(3)
        integer :: ierr

        if (l0_ntile <= 0) return

        ! periodic_bc is set on rank 0 only (s_read_input_file is rank-0-guarded), so make it globally consistent: every rank must
        ! build
        ! the same wrap-seam list in f_amr_seam / apply the same periodic edge fill, else an unmatched seam transfer deadlocks.
        ! MPI_LOR:
        ! rank 0's .true. wins on all ranks (others hold the .false. default). No-op at np=1.
        l0_periodic = periodic_bc
#ifdef MFC_MPI
        call MPI_ALLREDUCE(MPI_IN_PLACE, l0_periodic, 3, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
#endif

        ! Supported physical faces (any np): extrapolation (BC_GHOST_EXTRAP), reflective (BC_REFLECTIVE), periodic (BC_PERIODIC);
        ! each
        ! has a cons-space tile fill that commutes with the cons->prim convert so a tile matches the monolithic prim-space BC
        ! bit-for-
        ! bit. The characteristic/slip/dirichlet family (bc < BC_GHOST_EXTRAP) is not handled. Validate once here (host), not
        ! per
        ! stage.
        if (f_l0_bc_unsupported(bc_x%beg) .or. f_l0_bc_unsupported(bc_x%end) .or. (n_glb > 0 .and. (f_l0_bc_unsupported(bc_y%beg) &
            & .or. f_l0_bc_unsupported(bc_y%end))) .or. (p_glb > 0 .and. (f_l0_bc_unsupported(bc_z%beg) &
            & .or. f_l0_bc_unsupported(bc_z%end)))) then
            call s_mpi_abort('l0_ntile spike: unsupported physical BC (only extrapolation, reflective, periodic are handled)')
        end if
        ! the monolithic L0 RHS is skipped for l0_ntile>0, so the global q_prim_vf it populated is stale: run-time-info and probes
        ! (which read it at stage 1) are not supported with tiling.
        if (run_time_info .or. probe_wrt) then
            call s_mpi_abort('l0_ntile spike does not support run_time_info or probe_wrt (monolithic q_prim_vf is not maintained)')
        end if

        ! rr=1 makes the swap ghost-coord bisection and the fine-fine-halo fmul (=amr_ref_ratio**level) both identity; slots
        ! inherit it via s_amr_alloc_slot. Only clobber the global in pure-L0 (amr off): under coexist the global must stay the
        ! real 2/4 so fine blocks size correctly (s_set_amr_fine_geometry etc.), and tiles get rr=1 via the per-slot override in
        ! s_l0_build_tile_slot (level-0 tiles are rr=1 regardless of the global).
        if (.not. amr) amr_ref_ratio = 1

        ! Tiles are per-rank: each rank's local chunk is split into nt(:) pieces; the global tile table (region + owner) is the
        ! union
        ! over ranks, replicated on every rank. Total = num_procs * nt(1)*nt(2)*nt(3). Each rank allocates slot data (fields +
        ! coords) only for its own tiles; the seam-pair scan and fine-fine halo see
        ! the
        ! full table and exchange cross-rank seams over MPI.
        ! l0_nt/l0_ntiles_tot/l0_slot_off are computed once by s_initialize_amr_module (which always runs first, per
        ! m_start_up.fpp) so both inits agree on the shared-pool layout; just read them here.
        nt = l0_nt

        amr_num_levels = 1
        amr_cur = 1

        if (.not. amr) then
            ! l0-only mode: this routine owns the pool (fine budget = 0)
            amr_max_fine = 0; l0_slot_off = l0_ntiles_tot
            amr_max_blocks = l0_ntiles_tot
            amr_num_blocks = l0_ntiles_tot

            ! block-metadata pool (mirror of s_initialize_amr_module's allocation)
            allocate (amr_slots(1:amr_max_blocks))
            call s_amr_loc_index_init()
            allocate (amr_region_lo_all(3, amr_max_blocks), amr_region_hi_all(3, amr_max_blocks))
            allocate (amr_isect_lo_all(3, amr_max_blocks), amr_isect_hi_all(3, amr_max_blocks))
            allocate (amr_owns_all(amr_max_blocks))
            allocate (amr_block_owner(amr_max_blocks)); amr_block_owner = 0
            allocate (amr_owner_cut(0:num_procs - 1)); amr_owner_cut = -1_8
            allocate (amr_fine_cut(0:num_procs - 1,1:max(amr_max_level, 1))); amr_fine_cut = -1_8
            allocate (amr_tile_l0_owner(amr_max_blocks)); amr_tile_l0_owner = 0
            allocate (amr_tile_cost(amr_max_blocks)); amr_tile_cost = 0._wp
            allocate (amr_tile_cost_ema(amr_max_blocks)); amr_tile_cost_ema = 0._wp
            ! L0 tiles are the base level (fine blocks are level>=1 on a tile)
            allocate (amr_block_level(amr_max_blocks)); amr_block_level = 0
            ! 2D rank lists sized to the computed max overlap in s_amr_build_seam_pairs; only the per-block counts are sized here.
            allocate (amr_ovl_gather_n(amr_max_blocks), amr_ovl_scatter_n(amr_max_blocks))
            allocate (amr_slot_live(amr_max_blocks)); amr_slot_live = .false.
            amr_region_lo_all = 0; amr_region_hi_all = 0; amr_isect_lo_all = 0; amr_isect_hi_all = 0; amr_owns_all = .false.
        else
            ! coexist mode: s_initialize_amr_module already allocated the shared pool sized l0_slot_off+amr_max_fine.
            ! Only allocate the tile-specific side tables here, and do not touch amr_slots / amr_region_* / amr_owns_all /
            ! amr_block_owner / amr_ovl_*; those are shared with AMR and already sized/allocated.
            allocate (amr_tile_l0_owner(amr_max_blocks)); amr_tile_l0_owner = 0
            allocate (amr_tile_cost(amr_max_blocks)); amr_tile_cost = 0._wp
            allocate (amr_tile_cost_ema(amr_max_blocks)); amr_tile_cost_ema = 0._wp
            ! tiles are level 0 in slots [1..l0_ntiles_tot]; set that band without disturbing the fine slots
            amr_block_level(1:l0_ntiles_tot) = 0
        end if

        ! the per-rank coarse decomposition (global origin + local extent) that the tile geometry and max-tile-extent sizing below
        ! read for every rank is computed O(1) by s_amr_rank_decomp (no table, no allgather). In l0-only mode
        ! s_initialize_amr_module
        ! did not run, so validate the formula against this rank's actual decomposition here (coexist validates it in that routine).
        if (.not. amr) call s_amr_validate_decomp()

        ! max tile extent per dim over all ranks (= widest per-rank chunk split by nt); slots + seam buffers are sized to this
        ! global
        ! max so every rank's buffers match. Chunk r has s_amr_rank_decomp ext(d)+1 cells in dim d.
        ! Coexist: the fine sizing s_initialize_amr_module just computed must survive; these are module-level and are what
        ! s_amr_alloc_slot reads, so zeroing them here would leave the shared pool sized to the tile extent, and every slot a
        ! regrid allocates afterwards would be a fine block cut to tile size, writing past its own bounds (an out-of-range device
        ! write, not a host abort). Accumulate the max of both instead. l0-only (.not. amr): s_initialize_amr_module returned
        ! early, so no fine sizing exists; start at 0.
        if (.not. amr) then
            max_f1 = 0; max_f2 = 0; max_f3 = 0
        end if
        do r = 0, num_procs - 1
            call s_amr_rank_decomp(r, rsidx, rext)
            e = (rext(1) + 1 + nt(1) - 1)/nt(1) - 1; max_f1 = max(max_f1, e)
            if (n_glb > 0) then; e = (rext(2) + 1 + nt(2) - 1)/nt(2) - 1; max_f2 = max(max_f2, e); end if
            if (p_glb > 0) then; e = (rext(3) + 1 + nt(3) - 1)/nt(3) - 1; max_f3 = max(max_f3, e); end if
        end do
        mbuf1_lo = -buff_size; mbuf1_hi = max_f1 + buff_size
        mbuf2_lo = 0; mbuf2_hi = 0; mbuf3_lo = 0; mbuf3_hi = 0
        if (n_glb > 0) then; mbuf2_lo = -buff_size; mbuf2_hi = max_f2 + buff_size; end if
        if (p_glb > 0) then; mbuf3_lo = -buff_size; mbuf3_hi = max_f3 + buff_size; end if
        call s_amr_scr_init()  ! mbuf* now final (fine/tile union under coexist); scratch must exist on every rank

        amr_seam_pairs_dirty = .true.; amr_seam_pairs_nblk = -1
        amr_mesh_epoch = amr_mesh_epoch + 1

        ! swap bounce buffers (same bounds as the L0 global coord arrays). Shared with s_initialize_amr_module (identical m/n/p
        ! sizing), so only allocate in l0-only mode to avoid a coexist double-allocate; under coexist the AMR init's buffers already
        ! serve both the fine-block and the tile swaps.
        if (.not. amr) then
            allocate (sw_x_cb(-1 - buff_size:m_alloc + buff_size), sw_x_cc(-buff_size:m_alloc + buff_size), &
                      & sw_dx(-buff_size:m_alloc + buff_size))
            if (n_glb > 0) allocate (sw_y_cb(-1 - buff_size:n_alloc + buff_size), sw_y_cc(-buff_size:n_alloc + buff_size), &
                & sw_dy(-buff_size:n_alloc + buff_size))
            if (p_glb > 0) allocate (sw_z_cb(-1 - buff_size:p_alloc + buff_size), sw_z_cc(-buff_size:p_alloc + buff_size), &
                & sw_dz(-buff_size:p_alloc + buff_size))
        end if
        call s_l0_build_extended_global_cb()  ! global L0 boundaries extended into the domain ghost shell (edge tiles need it)
        amr_cpat_mar = (buff_size + amr_ref_ratio - 1)/amr_ref_ratio + 1
        amr_xchg_coarse_ghosts = .false.  ! tiles never prolong from a coarser level

        ! per-tile geometry: level-1 rr=1 blocks. Loop is rank-major then z,y,x (x fastest) so intra-rank neighbours are
        ! contiguous; the global region scan finds cross-rank seams. A rank writes region+owner for every tile but allocates slot
        ! data only for its own (r == proc_rank). Domain-edge detection uses the global region indices (region_lo == 0 / region_hi
        ! == m_glb).
        k = 0
        do r = 0, num_procs - 1
            call s_amr_rank_decomp(r, rsidx, rext)
            do iz = 0, nt(3) - 1
                do iy = 0, nt(2) - 1
                    do ix = 0, nt(1) - 1
                        k = k + 1
                        ! global cell range of tile (r, ix, iy, iz) = rank r's chunk split by f_l0_lo (identical split on every rank
                        ! so seam transverse extents match across ranks for an even decomposition)
                        tlo(1) = rsidx(1) + f_l0_lo(rext(1) + 1, nt(1), ix)
                        thi(1) = rsidx(1) + f_l0_lo(rext(1) + 1, nt(1), ix + 1) - 1
                        tlo(2) = 0; thi(2) = 0
                        if (n_glb > 0) then
                            tlo(2) = rsidx(2) + f_l0_lo(rext(2) + 1, nt(2), iy)
                            thi(2) = rsidx(2) + f_l0_lo(rext(2) + 1, nt(2), iy + 1) - 1
                        end if
                        tlo(3) = 0; thi(3) = 0
                        if (p_glb > 0) then
                            tlo(3) = rsidx(3) + f_l0_lo(rext(3) + 1, nt(3), iz)
                            thi(3) = rsidx(3) + f_l0_lo(rext(3) + 1, nt(3), iz + 1) - 1
                        end if
                        amr_block_owner(k) = r; amr_myblk_dirty = .true.
                        amr_tile_l0_owner(k) = r  ! L0 storage owner = init owner; stays fixed under migration
                        amr_owns_all(k) = (r == proc_rank)
                        amr_region_lo_all(:,k) = tlo; amr_region_hi_all(:,k) = thi
                        amr_isect_lo_all(:,k) = tlo; amr_isect_hi_all(:,k) = thi  ! footprint = whole tile on the owner (rr=1)
                        ! slot data is not allocated here: it must follow the SFC compute owner assigned below, which need not be
                        ! this cartesian owner r. Every rank writes region+owner metadata for every tile.
                    end do
                end do
            end do
        end do

        ! Compute owner = SFC cost-split over the tiles (fills the O(num_procs) amr_owner_cut), superseding the cartesian owner
        ! set in the loop above; amr_tile_l0_owner keeps the cartesian storage assignment (where the L0 field data physically
        ! lives). Init cost is geometric (whole-tile cell count), uniform on a uniform grid. The SFC compute owner may differ from
        ! the cartesian storage owner; whether it does depends on how the cartesian split direction lines up with Morton order, so
        ! it is a property of the grid shape (a 2:1 grid at np=2 splits in y and agrees; a square grid tie-breaks to x and
        ! diverges). s_l0_copy_coarse_to_tiles routes the initial fill storage-owner -> compute-owner when they differ, so no
        ! precondition is asserted here; the rebalancer's routed migration handles post-init divergence.
        block
            integer         :: sfco(l0_ntiles_tot), kk
            integer(kind=8) :: tkey(l0_ntiles_tot)
            real(wp)        :: twt(l0_ntiles_tot)
            do kk = 1, l0_ntiles_tot
                tkey(kk) = f_morton(amr_region_lo_all(1, kk), amr_region_lo_all(2, kk), amr_region_lo_all(3, kk))
                twt(kk) = real(amr_region_hi_all(1, kk) - amr_region_lo_all(1, kk) + 1, wp)*real(amr_region_hi_all(2, &
                    & kk) - amr_region_lo_all(2, kk) + 1, wp)*real(amr_region_hi_all(3, kk) - amr_region_lo_all(3, kk) + 1, wp)
            end do
            call s_amr_sfc_cut(tkey, twt, l0_ntiles_tot, amr_owner_cut, sfco)
            do kk = 1, l0_ntiles_tot
                amr_block_owner(kk) = sfco(kk); amr_myblk_dirty = .true.
                amr_owns_all(kk) = (sfco(kk) == proc_rank)
            end do
            ! allocate slot data for the tiles this rank computes (deferred from the cartesian loop above). s_l0_build_tile_slot
            ! reads only replicated region metadata and the global amr_g?cb, so it is valid for any tile on any rank.
            do kk = 1, l0_ntiles_tot
                if (amr_block_owner(kk) == proc_rank) call s_l0_build_tile_slot(kk)
            end do
        end block
        ! validate the full picture now that tiles exist: tiles vs amr_owner_cut (tile cut, just built), fine blocks vs amr_fine_cut
        ! (the level-1 cut the assigner saved at init). Runs in coexist too (amr_fine_cut is populated by the earlier assigner
        ! call).
        call s_amr_validate_owner()

        call s_amr_select_slot(1)
        ! tiles are persistent: L0 seeds them once at the first timestep (s_l0_copy_coarse_to_tiles self-gates on this flag), then
        ! they carry their own state across stages/timesteps. No init-time device copy (q_cons device state is not live at module
        ! init).
        l0_tiles_need_fill = .true.

    end subroutine s_l0_tiles_init

    !> Global L0 cell boundaries extended into the domain ghost shell (-1-buff_size : G+buff_size), unlike s_amr_build_global_cb
    !! (-1:G). The swap rebuilds a block's ghost-shell coordinates from these, and a tile that touches the domain boundary reaches
    !! indices beyond G (AMR fine blocks never do, being buff_size inside). Sourced from the monolithic x_cb, whose ghost cells
    !! already hold the domain's ghost coordinates, so a tile's ghost coords match the monolithic grid's bit-for-bit.
    impure subroutine s_l0_build_extended_global_cb()

        integer             :: j
        real(wp), parameter :: sentinel = -huge(1._wp)

        ! Under coexist, s_amr_build_global_cb (called by s_initialize_amr_module) already allocated amr_g?cb at the non-extended
        ! bounds (-1:G) for the fine-block geometry. The tiles need the extended bounds (-1-buff:G+buff) since an edge tile
        ! reaches into the domain ghost shell; the extended array is a value-consistent superset (same x_cb source over the
        ! overlap, and the fine geometry already copied its coords into the slots), so replace it. Without this the second
        ! allocate is a fatal error on gfortran and a silent double-allocate (leaked non-extended buffer) on flang.

        if (allocated(amr_gxcb)) deallocate (amr_gxcb)
        allocate (amr_gxcb(-1 - buff_size:m_glb + buff_size)); amr_gxcb = sentinel
        do j = -1 - buff_size, m + buff_size
            amr_gxcb(start_idx(1) + j) = x_cb(j)
        end do
        call s_mpi_allreduce_array_max(amr_gxcb, m_glb + 2 + 2*buff_size)
        if (n_glb > 0) then
            if (allocated(amr_gycb)) deallocate (amr_gycb)
            allocate (amr_gycb(-1 - buff_size:n_glb + buff_size)); amr_gycb = sentinel
            do j = -1 - buff_size, n + buff_size
                amr_gycb(start_idx(2) + j) = y_cb(j)
            end do
            call s_mpi_allreduce_array_max(amr_gycb, n_glb + 2 + 2*buff_size)
        end if
        if (p_glb > 0) then
            if (allocated(amr_gzcb)) deallocate (amr_gzcb)
            allocate (amr_gzcb(-1 - buff_size:p_glb + buff_size)); amr_gzcb = sentinel
            do j = -1 - buff_size, p + buff_size
                amr_gzcb(start_idx(3) + j) = z_cb(j)
            end do
            call s_mpi_allreduce_array_max(amr_gzcb, p_glb + 2 + 2*buff_size)
        end if

    end subroutine s_l0_build_extended_global_cb

    !> Build tile k's slot on this rank from its (already-set, replicated) region metadata: allocate the field/coord arrays and set
    !! the local extents, idwbuff, and rr=1 cell coordinates sliced from the global amr_g?cb. Shared by s_l0_tiles_init (initial
    !! owned tiles) and s_l0_migrate_tile (a tile arriving on its new owner). Requires amr_gxcb/gycb/gzcb + mbuf*/max_f* already
    !! set.
    impure subroutine s_l0_build_tile_slot(k)

        integer, intent(in) :: k
        integer             :: j, tlo(3), thi(3)

        tlo = amr_region_lo_all(:,k); thi = amr_region_hi_all(:,k)
        call s_amr_alloc_slot(k)  ! sizes to mbuf*, sets slot%amr_ref_ratio = amr_ref_ratio
        ! a base-level tile is rr=1 regardless of the global refinement ratio (the global may be 2/4 for fine blocks)
        amr_slots(k)%amr_ref_ratio = 1
        amr_slots(k)%m = thi(1) - tlo(1); amr_slots(k)%n = 0; amr_slots(k)%p = 0
        if (n_glb > 0) amr_slots(k)%n = thi(2) - tlo(2)
        if (p_glb > 0) amr_slots(k)%p = thi(3) - tlo(3)
        amr_slots(k)%idwbuff(1)%beg = -buff_size; amr_slots(k)%idwbuff(1)%end = amr_slots(k)%m + buff_size
        amr_slots(k)%idwbuff(2)%beg = 0; amr_slots(k)%idwbuff(2)%end = 0
        amr_slots(k)%idwbuff(3)%beg = 0; amr_slots(k)%idwbuff(3)%end = 0
        if (n_glb > 0) then
            amr_slots(k)%idwbuff(2)%beg = -buff_size; amr_slots(k)%idwbuff(2)%end = amr_slots(k)%n + buff_size
        end if
        if (p_glb > 0) then
            amr_slots(k)%idwbuff(3)%beg = -buff_size; amr_slots(k)%idwbuff(3)%end = amr_slots(k)%p + buff_size
        end if
        ! rr=1: tile cell j (right boundary) is the global L0 boundary amr_g?cb(tlo + j); interior coords only (the swap extends
        ! the ghost shell from amr_g?cb identically).
        do j = -1, amr_slots(k)%m
            amr_slots(k)%x_cb(j) = amr_gxcb(tlo(1) + j)
        end do
        do j = 0, amr_slots(k)%m
            amr_slots(k)%dx(j) = amr_slots(k)%x_cb(j) - amr_slots(k)%x_cb(j - 1)
            amr_slots(k)%x_cc(j) = 0.5_wp*(amr_slots(k)%x_cb(j - 1) + amr_slots(k)%x_cb(j))
        end do
        if (n_glb > 0) then
            do j = -1, amr_slots(k)%n
                amr_slots(k)%y_cb(j) = amr_gycb(tlo(2) + j)
            end do
            do j = 0, amr_slots(k)%n
                amr_slots(k)%dy(j) = amr_slots(k)%y_cb(j) - amr_slots(k)%y_cb(j - 1)
                amr_slots(k)%y_cc(j) = 0.5_wp*(amr_slots(k)%y_cb(j - 1) + amr_slots(k)%y_cb(j))
            end do
        end if
        if (p_glb > 0) then
            do j = -1, amr_slots(k)%p
                amr_slots(k)%z_cb(j) = amr_gzcb(tlo(3) + j)
            end do
            do j = 0, amr_slots(k)%p
                amr_slots(k)%dz(j) = amr_slots(k)%z_cb(j) - amr_slots(k)%z_cb(j - 1)
                amr_slots(k)%z_cc(j) = 0.5_wp*(amr_slots(k)%z_cb(j - 1) + amr_slots(k)%z_cb(j))
            end do
        end if

    end subroutine s_l0_build_tile_slot

    !> Copy the current L0 interior state into every owned tile's interior (global cell tlo+j -> tile-local cell j). A tile whose
    !! compute owner is also its L0-storage owner is seeded by a local device copy (the common case, and the entire path when the
    !! SFC cut agrees with the cartesian order). When the SFC compute owner differs from the cartesian storage owner the seed is
    !! routed: the L0-storage owner device-packs its chunk and sends it to the compute owner, which unpacks into its tile slot.
    !! Exact reverse of s_l0_scatter_tiles_to_coarse, and sound because a tile is built by subdividing one rank's cartesian chunk
    !! (s_l0_tiles_init), so it never spans two L0-storage ranks.
    impure subroutine s_l0_copy_coarse_to_tiles(q_cons_vf)

        ! inout (not in): passed as the bidirectional s_l0_copy_block q_l0 dummy (intent(inout)); read-only here (L0 -> tile)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf

        ! Persistent tiles: seed from L0 exactly once. After the first fill the tiles are authoritative; re-copying would be an
        ! identity round-trip (each stage scatters tile->L0, so L0 already mirrors the tile interior at the next timestep's stage
        ! 1).

        if (.not. l0_tiles_need_fill) return

        call s_l0_fill_tiles_from_coarse(q_cons_vf)
        l0_tiles_need_fill = .false.

    end subroutine s_l0_copy_coarse_to_tiles

    !> The fill itself, without the seed gate: overwrite every owned tile interior from the L0 field. Separate from
    !! s_l0_copy_coarse_to_tiles because the coexist subcycle path round-trips through L0 every step (tiles -> L0, fine fold writes
    !! L0, L0 -> tiles), so it needs this unconditionally, while the seed must still happen exactly once.
    impure subroutine s_l0_fill_tiles_from_coarse(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        integer                                                :: k, o1, o2, o3, fm1, fm2, fm3, bown, lown, cnt, ierr
        real(wp), allocatable                                  :: buf(:)

        do k = 1, l0_ntiles_tot
            bown = amr_block_owner(k); lown = amr_tile_l0_owner(k)
            if (bown == lown) then  ! compute owner holds the L0 cells: local device copy
                if (bown /= proc_rank) cycle
                call s_l0_tile_l0_offsets(k, o1, o2, o3)
                fm1 = amr_slots(k)%m; fm2 = amr_slots(k)%n; fm3 = amr_slots(k)%p
                call s_l0_copy_block(amr_loc_of(k), q_cons_vf, o1, o2, o3, fm1, fm2, fm3, .true.)
                cycle
            end if
            ! routed seed: extents come from the replicated region (the L0 owner has no slot for this tile)
            fm1 = amr_region_hi_all(1, k) - amr_region_lo_all(1, k)
            fm2 = 0; if (n_glb > 0) fm2 = amr_region_hi_all(2, k) - amr_region_lo_all(2, k)
            fm3 = 0; if (p_glb > 0) fm3 = amr_region_hi_all(3, k) - amr_region_lo_all(3, k)
            cnt = sys_size*(fm1 + 1)*(fm2 + 1)*(fm3 + 1)
            if (proc_rank == lown) then  ! L0-storage owner: device-pack the tile's L0 chunk, send to the compute owner
                call s_l0_tile_l0_offsets(k, o1, o2, o3)
                allocate (buf(cnt))
                call s_l0_pack_unpack_block_sf(q_cons_vf, o1, o2, o3, fm1, fm2, fm3, buf, .true.)
#ifdef MFC_MPI
                call s_xa_rec(XA_L0_FILL_SND, 1, cnt, k)
                call MPI_SEND(buf, cnt, mpi_p, bown, k, MPI_COMM_WORLD, ierr)
#endif
                deallocate (buf)
            else if (proc_rank == bown) then  ! compute owner: recv, device-unpack into the tile interior
                allocate (buf(cnt))
#ifdef MFC_MPI
                call s_xa_rec(XA_L0_FILL_RCV, 2, cnt, k)
                call MPI_RECV(buf, cnt, mpi_p, lown, k, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
#endif
                call s_l0_pack_unpack_block_st(amr_loc_of(k), 0, 0, 0, fm1, fm2, fm3, buf, .false.)
                deallocate (buf)
            end if
        end do

    end subroutine s_l0_fill_tiles_from_coarse

    !> Local-index offset of tile k's global origin in the L0 field: o(d) = region_lo(d) - start_idx(d) for active dims, 0 for a
    !! collapsed dim (start_idx is sized num_dims, so start_idx(3) must not be touched in 2D).
    subroutine s_l0_tile_l0_offsets(k, o1, o2, o3)

        integer, intent(in)  :: k
        integer, intent(out) :: o1, o2, o3

        o1 = amr_region_lo_all(1, k) - start_idx(1)
        o2 = 0; if (n_glb > 0) o2 = amr_region_lo_all(2, k) - start_idx(2)
        o3 = 0; if (p_glb > 0) o3 = amr_region_lo_all(3, k) - start_idx(3)

    end subroutine s_l0_tile_l0_offsets

    !> Scatter every tile's interior back into the L0 field (tile-local cell j -> global cell tlo+j). A tile whose compute owner is
    !! also its L0-storage owner writes locally (device kernel; the common case, and the entire no-migration path). A migrated tile
    !! (owner != l0_owner) has its interior sent by the compute owner to the L0-storage owner over MPI, which writes it into L0,
    !! keeping the fixed L0 decomposition (hence output/restart) correct after migration. Ghosts are not scattered (the tile path
    !! never reads L0 ghosts). GPU-correct: the MPI branch device-packs/unpacks via s_l0_pack_unpack_block, so the receiver writes
    !! L0 on the device; it survives the GPU_UPDATE(host) that s_save_data does before writing.
    impure subroutine s_l0_scatter_tiles_to_coarse(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        integer                                                :: k, o1, o2, o3, fm1, fm2, fm3, bown, lown, cnt, ierr
        real(wp), allocatable                                  :: buf(:)

        ! Precondition: tiles are the authoritative store. Before the first seed (s_l0_copy_coarse_to_tiles) the tile slots hold
        ! uninitialized (zero) state and L0 still holds the initial condition, so there is nothing to refresh; scattering here
        ! would overwrite the IC with zeros (zero density -> NaN once the coexist L0 coarse RHS consumes it). Skip until seeded.

        if (l0_tiles_need_fill) return

        do k = 1, l0_ntiles_tot
            bown = amr_block_owner(k); lown = amr_tile_l0_owner(k)
            fm1 = amr_region_hi_all(1, k) - amr_region_lo_all(1, k)
            fm2 = 0; if (n_glb > 0) fm2 = amr_region_hi_all(2, k) - amr_region_lo_all(2, k)
            fm3 = 0; if (p_glb > 0) fm3 = amr_region_hi_all(3, k) - amr_region_lo_all(3, k)
            if (bown == lown) then  ! not migrated: local device copy
                if (bown /= proc_rank) cycle
                call s_l0_tile_l0_offsets(k, o1, o2, o3)
                call s_l0_copy_block(amr_loc_of(k), q_cons_vf, o1, o2, o3, fm1, fm2, fm3, .false.)
                cycle
            end if
            cnt = sys_size*(fm1 + 1)*(fm2 + 1)*(fm3 + 1)
            if (proc_rank == bown) then  ! compute owner: device-pack owned tile interior, send to the L0 owner
                allocate (buf(cnt))
                call s_l0_pack_unpack_block_st(amr_loc_of(k), 0, 0, 0, fm1, fm2, fm3, buf, .true.)
#ifdef MFC_MPI
                call s_xa_rec(XA_L0_SCAT_SND, 1, cnt, k)
                call MPI_SEND(buf, cnt, mpi_p, lown, k, MPI_COMM_WORLD, ierr)
#endif
                deallocate (buf)
            else if (proc_rank == lown) then  ! L0 owner: recv, device-unpack into the local L0 chunk (device write -> survives the
                call s_l0_tile_l0_offsets(k, o1, o2, o3)  ! GPU_UPDATE(host) s_save_data does before writing)
                allocate (buf(cnt))
#ifdef MFC_MPI
                call s_xa_rec(XA_L0_SCAT_RCV, 2, cnt, k)
                call MPI_RECV(buf, cnt, mpi_p, bown, k, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
#endif
                call s_l0_pack_unpack_block_sf(q_cons_vf, o1, o2, o3, fm1, fm2, fm3, buf, .false.)
                deallocate (buf)
            end if
        end do

    end subroutine s_l0_scatter_tiles_to_coarse

    !> Device add of an L0 block [o+0:o+fm] into a tile rhs interior [0:fm] (q_rhs += q_l0). Additive twin of s_l0_copy_block's
    !! to_tile branch, in wp (the rhs is computed in wp, stored stp). Slot rhs is a dummy so the kernel reads a valid mapped
    !! descriptor (indexing module amr_slots%rhs in a kernel is a null deref; see s_l0_copy_block / s_amr_fine_slice).
    impure subroutine s_l0_add_block(q_rhs, q_l0, o1, o2, o3, fm1, fm2, fm3)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_rhs
        type(scalar_field), dimension(sys_size), intent(in)    :: q_l0
        integer, intent(in)                                    :: o1, o2, o3, fm1, fm2, fm3
        integer                                                :: i, j, k, l

        $:GPU_PARALLEL_LOOP(collapse=4)
        do i = 1, sys_size
            do l = 0, fm3
                do k = 0, fm2
                    do j = 0, fm1
                        q_rhs(i)%sf(j, k, l) = real(real(q_rhs(i)%sf(j, k, l), wp) + real(q_l0(i)%sf(o1 + j, o2 + k, o3 + l), &
                              & wp), stp)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_l0_add_block

    !> Device add of the contiguous MPI buffer buf into a tile rhs interior [0:fm] (q_rhs += buf). Additive twin of
    !! s_l0_pack_unpack_block's unpack branch (same j-fastest buf layout), in wp. Slot rhs passed as a dummy (GPU-safe).
    impure subroutine s_l0_unpack_add_block(q_rhs, fm1, fm2, fm3, buf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_rhs
        integer, intent(in)                                    :: fm1, fm2, fm3
        real(wp), intent(inout), contiguous                    :: buf(:)
        integer                                                :: i, j, k, l

        $:GPU_PARALLEL_LOOP(collapse=4, copyin='[buf]')
        do i = 1, sys_size
            do l = 0, fm3
                do k = 0, fm2
                    do j = 0, fm1
                        q_rhs(i)%sf(j, k, l) = real(real(q_rhs(i)%sf(j, k, l), &
                              & wp) + buf(1 + j + (fm1 + 1)*(k + (fm2 + 1)*(l + (fm3 + 1)*(i - 1)))), stp)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_l0_unpack_add_block

    !> Coexist reflux copy-back: add the fixed-L0-frame Berger-Colella reflux delta into each tile's per-slot rhs on its (possibly
    !! migrated) compute owner, so the tile RK update sees the same c/f-face correction the monolithic coarse update would. The
    !! delta lives in rhs_delta (the L0 rhs, which s_tvd_rk zeroed before the fine loop so s_amr_apply_reflux filled it with the
    !! pure delta). Reverse of s_l0_scatter_tiles_to_coarse: source is the fixed L0-storage owner (amr_tile_l0_owner), dest is the
    !! compute owner (amr_block_owner); local when they coincide, else P2P (L0-owner packs the tile's L0 region, compute-owner adds
    !! it). Whole tile interior is routed; the delta is zero outside the c/f reflux shell, so the add is identity elsewhere.
    impure subroutine s_l0_add_reflux_to_tiles(rhs_delta)

        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_delta
        integer                                                :: k, o1, o2, o3, fm1, fm2, fm3, bown, lown, cnt, ierr
        real(wp), allocatable                                  :: buf(:)

        do k = 1, l0_ntiles_tot
            bown = amr_block_owner(k); lown = amr_tile_l0_owner(k)
            fm1 = amr_region_hi_all(1, k) - amr_region_lo_all(1, k)
            fm2 = 0; if (n_glb > 0) fm2 = amr_region_hi_all(2, k) - amr_region_lo_all(2, k)
            fm3 = 0; if (p_glb > 0) fm3 = amr_region_hi_all(3, k) - amr_region_lo_all(3, k)
            if (bown == lown) then  ! not migrated: local device add
                if (bown /= proc_rank) cycle
                call s_l0_tile_l0_offsets(k, o1, o2, o3)
                call s_l0_add_block(amr_slots(k)%rhs, rhs_delta, o1, o2, o3, fm1, fm2, fm3)
                cycle
            end if
            cnt = sys_size*(fm1 + 1)*(fm2 + 1)*(fm3 + 1)
            if (proc_rank == lown) then  ! L0 owner: device-pack the delta over this tile's L0 region, send to the compute owner
                call s_l0_tile_l0_offsets(k, o1, o2, o3)
                allocate (buf(cnt))
                call s_l0_pack_unpack_block_sf(rhs_delta, o1, o2, o3, fm1, fm2, fm3, buf, .true.)
#ifdef MFC_MPI
                call s_xa_rec(XA_L0_RFLX_SND, 1, cnt, k)
                call MPI_SEND(buf, cnt, mpi_p, bown, k, MPI_COMM_WORLD, ierr)
#endif
                deallocate (buf)
            else if (proc_rank == bown) then  ! compute owner: recv, device-add the delta into the tile rhs
                allocate (buf(cnt))
#ifdef MFC_MPI
                call s_xa_rec(XA_L0_RFLX_RCV, 2, cnt, k)
                call MPI_RECV(buf, cnt, mpi_p, lown, k, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
#endif
                call s_l0_unpack_add_block(amr_slots(k)%rhs, fm1, fm2, fm3, buf)
                deallocate (buf)
            end if
        end do

    end subroutine s_l0_add_reflux_to_tiles

    !> Coexist restrict copy-back: after the fine blocks restrict their solution into the L0 covered cells (fixed-L0-frame q_cons),
    !! overwrite each covering tile's matching cells with those restricted values on the tile's (possibly migrated) compute owner
    !! (the coexist twin of the monolithic level-0 covered-cell overwrite). Only the covered footprint moves (non-covered tile cells
    !! keep their advanced state; disjoint from the reflux shell). Per (tile, level-1 block) footprint intersection in the L0 frame:
    !! local when L0-owner == compute-owner (buffer roundtrip), else P2P (L0-owner packs the intersection, compute-owner unpacks).
    !! Reuses s_l0_pack_unpack_block with per-side offsets (its offset arg is per-call, so src L0 and dst tile offsets differ).
    impure subroutine s_l0_restrict_to_tiles(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        integer                                                :: k, b, d, bown, lown, cnt, ierr
        integer                                                :: ilo(3), ihi(3), e1, e2, e3, lo1, lo2, lo3, to1, to2, to3
        real(wp), allocatable                                  :: buf(:)
        logical                                                :: nonempty

        do k = 1, l0_ntiles_tot  ! tiles are the level-0 prefix
            bown = amr_block_owner(k); lown = amr_tile_l0_owner(k)
            do b = 1, amr_num_blocks
                ! only level-1 fine blocks restrict into L0 covered cells (level>=2 fold to parent)
                if (amr_block_level(b) /= 1) cycle
                nonempty = .true.  ! L0-frame intersection of tile k's region and fine block b's footprint
                do d = 1, 3
                    ilo(d) = max(amr_region_lo_all(d, k), amr_region_lo_all(d, b))
                    ihi(d) = min(amr_region_hi_all(d, k), amr_region_hi_all(d, b))
                    if (ilo(d) > ihi(d)) nonempty = .false.
                end do
                if (.not. nonempty) cycle
                e1 = ihi(1) - ilo(1)
                e2 = 0; if (n_glb > 0) e2 = ihi(2) - ilo(2)
                e3 = 0; if (p_glb > 0) e3 = ihi(3) - ilo(3)
                lo1 = ilo(1) - start_idx(1); to1 = ilo(1) - amr_region_lo_all(1, k)  ! L0-local (src) vs tile-local (dst) offsets
                lo2 = 0; to2 = 0
                if (n_glb > 0) then; lo2 = ilo(2) - start_idx(2); to2 = ilo(2) - amr_region_lo_all(2, k); end if
                lo3 = 0; to3 = 0
                if (p_glb > 0) then; lo3 = ilo(3) - start_idx(3); to3 = ilo(3) - amr_region_lo_all(3, k); end if
                cnt = sys_size*(e1 + 1)*(e2 + 1)*(e3 + 1)
                if (bown == lown) then  ! not migrated: local device pack (L0 region) -> unpack (tile region), same rank
                    if (bown /= proc_rank) cycle
                    allocate (buf(cnt))
                    call s_l0_pack_unpack_block_sf(q_cons_vf, lo1, lo2, lo3, e1, e2, e3, buf, .true.)
                    call s_l0_pack_unpack_block_st(amr_loc_of(k), to1, to2, to3, e1, e2, e3, buf, .false.)
                    deallocate (buf)
                else if (proc_rank == lown) then  ! L0 owner: device-pack the intersection, send to the compute owner
                    allocate (buf(cnt))
                    call s_l0_pack_unpack_block_sf(q_cons_vf, lo1, lo2, lo3, e1, e2, e3, buf, .true.)
#ifdef MFC_MPI
                    call s_xa_rec(XA_L0_REST_SND, 1, cnt, 4400 + k)
                    call MPI_SEND(buf, cnt, mpi_p, bown, 4400 + k, MPI_COMM_WORLD, ierr)
#endif
                    deallocate (buf)
                else if (proc_rank == bown) then  ! compute owner: recv, device-unpack (overwrite) into the tile covered cells
                    allocate (buf(cnt))
#ifdef MFC_MPI
                    call s_xa_rec(XA_L0_REST_RCV, 2, cnt, 4400 + k)
                    call MPI_RECV(buf, cnt, mpi_p, lown, 4400 + k, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
#endif
                    call s_l0_pack_unpack_block_st(amr_loc_of(k), to1, to2, to3, e1, e2, e3, buf, .false.)
                    deallocate (buf)
                end if
            end do
        end do

    end subroutine s_l0_restrict_to_tiles

    !> Migrate tile k from its current compute owner to new_owner: P2P-move the persistent interior state, (re)build the slot on the
    !! receiver, free it on the sender, and update the replicated owner map + seam topology. All ranks call with the same (k,
    !! new_owner). This is the load-balance migration primitive; the decision of which tile moves where is made by the caller
    !! (s_l0_forced_remap or s_l0_rebalance). Ghosts are not moved (refilled by edge-BC + fine-fine halo before the next stage).
    !! GPU-correct: interior is device-packed/unpacked via s_l0_pack_unpack_block (wp buffer, cast to/from stp).
    impure subroutine s_l0_migrate_tile(k, new_owner)

        integer, intent(in)   :: k, new_owner
        integer               :: old_owner, ni, nj, nl, cnt, ierr
        real(wp), allocatable :: buf(:)

        old_owner = amr_block_owner(k)
        if (old_owner == new_owner) return

        ni = amr_region_hi_all(1, k) - amr_region_lo_all(1, k)
        nj = 0; if (n_glb > 0) nj = amr_region_hi_all(2, k) - amr_region_lo_all(2, k)
        nl = 0; if (p_glb > 0) nl = amr_region_hi_all(3, k) - amr_region_lo_all(3, k)
        cnt = sys_size*(ni + 1)*(nj + 1)*(nl + 1)

        if (proc_rank == old_owner) then  ! device-pack + send the interior, then release the slot
            allocate (buf(cnt))
            call s_l0_pack_unpack_block_st(amr_loc_of(k), 0, 0, 0, ni, nj, nl, buf, .true.)
#ifdef MFC_MPI
            call s_xa_rec(XA_L0_MIGR_SND, 1, cnt, 4300)
            call MPI_SEND(buf, cnt, mpi_p, new_owner, 4300, MPI_COMM_WORLD, ierr)
#endif
            deallocate (buf)
            call s_amr_free_slot(k)
        else if (proc_rank == new_owner) then  ! build the slot, recv + device-unpack the interior into it
            call s_l0_build_tile_slot(k)
            allocate (buf(cnt))
#ifdef MFC_MPI
            call s_xa_rec(XA_L0_MIGR_RCV, 2, cnt, 4300)
            call MPI_RECV(buf, cnt, mpi_p, old_owner, 4300, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
#endif
            call s_l0_pack_unpack_block_st(amr_loc_of(k), 0, 0, 0, ni, nj, nl, buf, .false.)
            deallocate (buf)
        end if

        ! replicated ownership update on every rank; mark the seam topology dirty so the next halo rebuilds pair/overlap lists.
        ! The epoch bump matters most here: ownership changed with no regrid, which the consumed boolean cannot express to a
        ! cached exchange plan.
        amr_block_owner(k) = new_owner; amr_myblk_dirty = .true.
        amr_owns_all(k) = (new_owner == proc_rank)
        amr_seam_pairs_dirty = .true.
        amr_mesh_epoch = amr_mesh_epoch + 1

    end subroutine s_l0_migrate_tile

    !> Test hook: at t_step == l0_migrate_step, force-migrate the last tile (initially owned by rank num_procs-1) to rank 0,
    !! exercising the migration primitive + seam-topology rebuild. Output must stay byte-identical to the no-migration run. No-op at
    !! np=1 (the last tile already lives on rank 0). All ranks call with identical arguments.
    impure subroutine s_l0_forced_remap()

        integer :: k

        k = l0_ntiles_tot
        if (amr_block_owner(k) /= 0) call s_l0_migrate_tile(k, 0)

    end subroutine s_l0_forced_remap

    !> Closed-loop rebalancer driven by measured per-tile compute time. Each rank accumulated amr_tile_cost for its own tiles since
    !! the last rebalance; an allreduce(SUM) makes the full cost vector replicated and bit-identical on every rank (each tile has
    !! exactly one nonzero contributor, so the sum is exact), so every rank runs the identical re-cut and issues matching P2P
    !! migrations. The re-cut is an SFC cost-split of the smoothed cost, with a small relative deadband so timing noise does not
    !! cause churn. Because migration is bit-preserving and the decision touches no field data, output is byte-identical regardless
    !! of the (run-to-run nondeterministic) measured schedule; the tile path's decomposition invariance is what keeps a
    !! nondeterministic cost signal golden-safe. Costs reset after each rebalance. No-op at np=1.
    impure subroutine s_l0_rebalance(t_step)

        integer, intent(in) :: t_step
        integer             :: k, nmig, ierr
        integer             :: newo(l0_ntiles_tot)
        real(wp)            :: cost(l0_ntiles_tot), load(0:num_procs - 1), gap0, gap1, mean, tol
        real(wp), parameter :: ema_hist = 0.5_wp  ! weight on the running estimate vs this window's measurement

        if (num_procs < 2) then
            amr_tile_cost = 0._wp  ! nothing to balance; still clear the window
            return
        end if

        cost = amr_tile_cost  ! local: nonzero only for this rank's owned tiles
#ifdef MFC_MPI
        call MPI_ALLREDUCE(MPI_IN_PLACE, cost, l0_ntiles_tot, mpi_p, MPI_SUM, MPI_COMM_WORLD, ierr)  ! -> replicated, bit-identical
#endif
        ! smooth the (replicated) window cost with a per-tile EMA so GPU per-tile launch-latency noise does not drive spurious
        ! migrations; seed on the first window (ema still all-zero) with the raw measurement to avoid a cold-start bias toward 0.
        ! amr_tile_cost_ema is derived only from the replicated cost, so it stays bit-identical on every rank -> consistent
        ! decision.
        if (all(amr_tile_cost_ema(1:l0_ntiles_tot) == 0._wp)) then
            amr_tile_cost_ema(1:l0_ntiles_tot) = cost
        else
            amr_tile_cost_ema(1:l0_ntiles_tot) = ema_hist*amr_tile_cost_ema(1:l0_ntiles_tot) + (1._wp - ema_hist)*cost
        end if
        cost = amr_tile_cost_ema(1:l0_ntiles_tot)  ! decide on the smoothed cost
        newo = amr_block_owner
        load = 0._wp
        do k = 1, l0_ntiles_tot
            load(newo(k)) = load(newo(k)) + cost(k)
        end do
        mean = sum(load)/real(num_procs, wp)
        tol = 0.05_wp*mean  ! deadband: ignore imbalance below 5% of the mean load so measurement noise does not churn migrations
        gap0 = maxval(load) - minval(load)

        ! SFC weighted re-cut: Morton-sort the tiles, then cumulative-split the smoothed cost into num_procs contiguous SFC
        ! ranges, the same partition logic as s_amr_assign_block_owners' cut. Ownership stays SFC-contiguous and O(num_procs)
        ! cut-derivable, and locality is preserved. Deadband: skip the re-cut while the load gap is already within tol (no churn
        ! on sub-5% imbalance).
        if (gap0 > tol) then
            block
                integer(kind=8) :: tkey(l0_ntiles_tot), cut_try(0:num_procs - 1)
                integer         :: newo_try(l0_ntiles_tot)
                real(wp)        :: load_try(0:num_procs - 1)
                do k = 1, l0_ntiles_tot
                    tkey(k) = f_morton(amr_region_lo_all(1, k), amr_region_lo_all(2, k), amr_region_lo_all(3, k))
                end do
                ! shared SFC cut: cumulative-split the smoothed cost into num_procs contiguous Morton ranges. Because the cut is
                ! restricted to contiguous ranges, the finest correction it can make is one whole tile; when tiles-per-rank is
                ! small that quantum exceeds the gap the deadband admits, and the re-cut can return a partition worse than the
                ! current one. Evaluate into a temporary and reject only a strict worsening. Not "accept only a strict
                ! improvement": with near-uniform tile costs the Morton partition and the cartesian one have equal gaps, and
                ! rejecting those would silently skip the migration the SFC re-cut rebalance golden exists to cover (it would
                ! still pass, because migration is bit-neutral; the coverage would just be gone). amr_owner_cut must move with
                ! newo: f_amr_owner resolves tile ownership against it, so refreshing it without migrating would leave it
                ! disagreeing with amr_block_owner.
                call s_amr_sfc_cut(tkey, cost, l0_ntiles_tot, cut_try, newo_try)
                load_try = 0._wp
                do k = 1, l0_ntiles_tot
                    load_try(newo_try(k)) = load_try(newo_try(k)) + cost(k)
                end do
                if (maxval(load_try) - minval(load_try) <= gap0) then
                    amr_owner_cut = cut_try
                    newo = newo_try
                end if
            end block
        end if
        load = 0._wp
        do k = 1, l0_ntiles_tot
            load(newo(k)) = load(newo(k)) + cost(k)
        end do
        gap1 = maxval(load) - minval(load)

        nmig = 0
        do k = 1, l0_ntiles_tot
            if (newo(k) /= amr_block_owner(k)) then
                call s_l0_migrate_tile(k, newo(k))
                nmig = nmig + 1
            end if
        end do
        ! tol is printed because without it "deadband skipped the re-cut" (gap0 <= tol) and "the re-cut ran and was rejected for
        ! not improving the gap" (gap0 > tol, guard above) emit identical output (gap0 -> gap0, 0 migrations). Those are different
        ! behaviours and one of them is the guard doing its job; a verification run cannot tell them apart otherwise.
        if (proc_rank == 0) print '(A,I0,A,ES10.3,A,ES10.3,A,ES10.3,A,I0,A)', ' [l0 rebalance] t_step=', t_step, ' load-gap ', &
            & gap0, ' -> ', gap1, ' (deadband ', tol, ', ', nmig, ' migrations)'

        amr_tile_cost = 0._wp  ! reset the measurement window

    end subroutine s_l0_rebalance

    !> Device copy between a tile interior [0:fm] and the L0 field [o+0:o+fm]. to_tile=T copies L0->tile, F copies tile->L0. The
    !! slot is addressed through the flat store, which is a plain GPU_DECLARE'd module array (indexing a per-slot scalar_field array
    !! inside a kernel would be a null deref; see s_amr_fine_slice).
    impure subroutine s_l0_copy_block(loc, q_l0, o1, o2, o3, fm1, fm2, fm3, to_tile)

        integer, intent(in)                                    :: loc
        type(scalar_field), dimension(sys_size), intent(inout) :: q_l0
        integer, intent(in)                                    :: o1, o2, o3, fm1, fm2, fm3
        logical, intent(in)                                    :: to_tile
        integer                                                :: i, j, k, l

        if (to_tile) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do l = 0, fm3
                    do k = 0, fm2
                        do j = 0, fm1
                            amr_cons_st(j, k, l, i, loc) = q_l0(i)%sf(o1 + j, o2 + k, o3 + l)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        else
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do l = 0, fm3
                    do k = 0, fm2
                        do j = 0, fm1
                            q_l0(i)%sf(o1 + j, o2 + k, o3 + l) = amr_cons_st(j, k, l, i, loc)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_l0_copy_block

    !> Fill each tile's domain-edge face ghosts with the physical BC (interior-seam faces are overwritten by s_amr_fine_fine_halo
    !! afterward). Supported: extrapolation, reflective, periodic (see f_l0_bc_unsupported); other codes abort. A face (d,side) is a
    !! domain edge iff the tile touches the global boundary there.
    impure subroutine s_l0_fill_edge_bc()

        integer :: k, fm(3), gcell(3)

        gcell(1) = m_glb; gcell(2) = n_glb; gcell(3) = p_glb
        do k = 1, l0_ntiles_tot
            if (amr_block_owner(k) /= proc_rank) cycle
            fm(1) = amr_slots(k)%m; fm(2) = amr_slots(k)%n; fm(3) = amr_slots(k)%p
            call s_l0_edge_bc_tile(amr_loc_of(k), amr_region_lo_all(1, k), amr_region_hi_all(1, k), gcell(1), fm, 1, bc_x%beg, &
                                   & bc_x%end)
            if (n_glb > 0) call s_l0_edge_bc_tile(amr_loc_of(k), amr_region_lo_all(2, k), amr_region_hi_all(2, k), gcell(2), fm, &
                & 2, bc_y%beg, bc_y%end)
            if (p_glb > 0) call s_l0_edge_bc_tile(amr_loc_of(k), amr_region_lo_all(3, k), amr_region_hi_all(3, k), gcell(3), fm, &
                & 3, bc_z%beg, bc_z%end)
        end do

    end subroutine s_l0_fill_edge_bc

    !> One tile, one dimension d: extrapolate the low/high face ghost from the edge interior cell, but only where the tile touches
    !! the domain boundary (rlo==0 low / rhi==gcell high). Applied to q_cons; convert (identity-commuting for extrapolation) makes
    !! the prim ghost the monolithic path produces. The transverse loop spans the face interior only (dimension-split scheme reads
    !! no corner ghost).
    impure subroutine s_l0_edge_bc_tile(loc, rlo, rhi, gcell, fm, d, bcbeg, bcend)

        integer, intent(in) :: loc
        integer, intent(in) :: rlo, rhi, gcell, fm(3), d, bcbeg, bcend

        ! BC support is validated once at init (s_l0_tiles_init); here we only apply it at domain-edge faces. Periodicity is read
        ! from
        ! the global periodic_bc(d) (not bcbeg, which becomes a wrap-neighbour rank at a decomposed periodic boundary): a periodic
        ! dim
        ! wraps; a tile that spans it (rlo==0 .and. rhi==gcell) self-wraps here, a partial tile's periodic faces are cross-tile
        ! wrap-seams filled by s_amr_fine_fine_halo (skipped here). A non-periodic domain-edge face gets reflective (mirror + normal
        ! momentum flip) or 0th-order extrapolation per its physical bc code.

        if (l0_periodic(d)) then
            if (rlo == 0 .and. rhi == gcell) call s_l0_wrap_one(loc, d, fm)
            return
        end if
        if (rlo == 0) then
            if (bcbeg == BC_REFLECTIVE) then; call s_l0_reflect_one(loc, d, -1, fm); else; call s_l0_extrap_one(loc, d, -1, &
                & fm); end if
        end if
        if (rhi == gcell) then
            if (bcend == BC_REFLECTIVE) then; call s_l0_reflect_one(loc, d, 1, fm); else; call s_l0_extrap_one(loc, d, 1, &
                & fm); end if
        end if

    end subroutine s_l0_edge_bc_tile

    !> Extrapolate tile face ghosts in dim d, side (-1 low / +1 high): ghost cells 1..buff_size = the edge interior cell (0 or md).
    !! Transverse extents (na, nb) and md are read into scalars before the device region (no host array element in the kernel).
    impure subroutine s_l0_extrap_one(loc, d, side, fm)

        integer, intent(in) :: loc
        integer, intent(in) :: d, side, fm(3)
        integer             :: i, jg, a, b, e, gc, na, nb, md

        #:for D, TA, TB in [(1, 2, 3), (2, 1, 3), (3, 1, 2)]
            #:set SIDX = {1: 'e, a, b', 2: 'a, e, b', 3: 'a, b, e'}[D]
            #:set GIDX = {1: 'gc, a, b', 2: 'a, gc, b', 3: 'a, b, gc'}[D]
            if (d == ${D}$) then
                na = fm(${TA}$); nb = fm(${TB}$); md = fm(${D}$)
                e = merge(0, md, side == -1)
                $:GPU_PARALLEL_LOOP(collapse=3, private='[gc]')
                do i = 1, sys_size
                    do b = 0, nb
                        do a = 0, na
                            do jg = 1, buff_size
                                gc = merge(-jg, md + jg, side == -1)
                                amr_cons_st(${GIDX}$, i, loc) = amr_cons_st(${SIDX}$, i, loc)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

    end subroutine s_l0_extrap_one

    !> Reflective (symmetry) tile face ghosts in dim d, side (-1 low / +1 high): ghost cell 1..buff_size mirrors the near-edge
    !! interior (ghost -jg <- interior jg-1 low; md+jg <- md-(jg-1) high) with the normal-direction momentum (eqn_idx%mom%beg + d -
    !! 1) negated, all other conserved variables copied. Done on q_cons; negating conserved normal momentum commutes with the
    !! cons->prim convert (velocity flips, rho and mom**2, hence pressure, are unchanged), so this reproduces the monolithic
    !! prim-space s_symmetry bit-for-bit. Transverse extent is the face interior only (dimension-split reads no corner ghost),
    !! matching s_l0_extrap_one.
    impure subroutine s_l0_reflect_one(loc, d, side, fm)

        integer, intent(in) :: loc
        integer, intent(in) :: d, side, fm(3)
        integer             :: i, jg, a, b, gc, sc, na, nb, md, nrm

        #:for D, TA, TB in [(1, 2, 3), (2, 1, 3), (3, 1, 2)]
            #:set SIDX = {1: 'sc, a, b', 2: 'a, sc, b', 3: 'a, b, sc'}[D]
            #:set GIDX = {1: 'gc, a, b', 2: 'a, gc, b', 3: 'a, b, gc'}[D]
            if (d == ${D}$) then
                na = fm(${TA}$); nb = fm(${TB}$); md = fm(${D}$)
                nrm = eqn_idx%mom%beg + ${D}$ - 1
                $:GPU_PARALLEL_LOOP(collapse=3, private='[gc, sc]')
                do i = 1, sys_size
                    do b = 0, nb
                        do a = 0, na
                            do jg = 1, buff_size
                                gc = merge(-jg, md + jg, side == -1)
                                sc = merge(jg - 1, md - (jg - 1), side == -1)
                                amr_cons_st(${GIDX}$, i, loc) = merge(-amr_cons_st(${SIDX}$, i, loc), amr_cons_st(${SIDX}$, i, &
                                            & loc), i == nrm)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

    end subroutine s_l0_reflect_one

    !> Periodic self-wrap for a tile that spans dim d (its low and high faces are both the domain boundary, i.e. l0_ntile==1 in d):
    !! fill both ghost shells from the opposite-end interior of the same tile: low ghost -jg <- interior md-(jg-1), high ghost md+jg
    !! <- interior jg-1 (a pure copy, matching the monolithic prim-space s_periodic; copy commutes with cons->prim convert). Partial
    !! tiles never reach here (their periodic faces are cross-tile wrap-seams handled by s_amr_fine_fine_halo). Face interior only.
    impure subroutine s_l0_wrap_one(loc, d, fm)

        integer, intent(in) :: loc
        integer, intent(in) :: d, fm(3)
        integer             :: i, jg, a, b, glo, shi, ghi, slo, na, nb, md

        #:for D, TA, TB in [(1, 2, 3), (2, 1, 3), (3, 1, 2)]
            #:set GLO = {1: 'glo, a, b', 2: 'a, glo, b', 3: 'a, b, glo'}[D]
            #:set SHI = {1: 'shi, a, b', 2: 'a, shi, b', 3: 'a, b, shi'}[D]
            #:set GHI = {1: 'ghi, a, b', 2: 'a, ghi, b', 3: 'a, b, ghi'}[D]
            #:set SLO = {1: 'slo, a, b', 2: 'a, slo, b', 3: 'a, b, slo'}[D]
            if (d == ${D}$) then
                na = fm(${TA}$); nb = fm(${TB}$); md = fm(${D}$)
                $:GPU_PARALLEL_LOOP(collapse=3, private='[glo, shi, ghi, slo]')
                do i = 1, sys_size
                    do b = 0, nb
                        do a = 0, na
                            do jg = 1, buff_size
                                glo = -jg; shi = md - (jg - 1)
                                ghi = md + jg; slo = jg - 1
                                amr_cons_st(${GLO}$, i, loc) = amr_cons_st(${SHI}$, i, loc)
                                amr_cons_st(${GHI}$, i, loc) = amr_cons_st(${SLO}$, i, loc)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

    end subroutine s_l0_wrap_one

    !> Fill a tile's multi-dim ghost cells (>= 2 dims in the ghost region: 2D diagonal corners, 3D ghost edges + corners) from the
    !! nearest interior cell (each index clamped to [0:m]/[0:n]/[0:p]). The dimension-split face fills only set single-ghost face
    !! slabs (they are all the RHS stencil reads), leaving these unset; the cons->prim convert still visits them and an unset ghost
    !! (void fractions 0 -> gamma 0) is a 0/0 (traps under -ffpe-trap, a stray NaN otherwise). The clamp source is always an
    !! interior cell (valid, real data) and never a face ghost, so the RHS-relevant ghosts are untouched and output is
    !! bit-unchanged.
    impure subroutine s_l0_fill_ghost_corners(loc, mx, ny, pz)

        integer, intent(in) :: loc
        integer, intent(in) :: mx, ny, pz
        integer             :: i, jb, kb, lb, jc, kc, lc, ng, lo2, hi2, lo3, hi3

        lo2 = 0; hi2 = 0; if (n_glb > 0) then; lo2 = -buff_size; hi2 = ny + buff_size; end if
        lo3 = 0; hi3 = 0; if (p_glb > 0) then; lo3 = -buff_size; hi3 = pz + buff_size; end if
        $:GPU_PARALLEL_LOOP(collapse=4, private='[jc, kc, lc, ng]')
        do i = 1, sys_size
            do lb = lo3, hi3
                do kb = lo2, hi2
                    do jb = -buff_size, mx + buff_size
                        ng = 0
                        if (jb < 0 .or. jb > mx) ng = ng + 1
                        if (n_glb > 0) then; if (kb < 0 .or. kb > ny) ng = ng + 1; end if
                        if (p_glb > 0) then; if (lb < 0 .or. lb > pz) ng = ng + 1; end if
                        if (ng >= 2) then
                            jc = min(max(jb, 0), mx); kc = min(max(kb, 0), ny); lc = min(max(lb, 0), pz)
                            amr_cons_st(jb, kb, lb, i, loc) = amr_cons_st(jc, kc, lc, i, loc)
                        end if
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_l0_fill_ghost_corners

    !> Advance every tile one RK stage: fill domain-edge BC ghosts (phase 1), overwrite interior-seam ghosts with neighbour interior
    !! (phase 2, fmul=1), then advance + RK-update each tile through the shared swap-based solver (phase 3). Mirrors the AMR
    !! fine-block phase structure (m_time_steppers) with prolong/reflux dropped (a base-res tile has no coarser level).
    !> Advance every owned tile one RK stage. Fused wrapper (pure-L0): RHS pass then RK pass back-to-back. Under coexist (amr)
    !! s_tvd_rk instead calls the two passes directly with s_l0_add_reflux_to_tiles between them, so each tile's coarse rhs is
    !! Berger-Colella corrected (from the fixed-L0-frame reflux) before its RK update.
    impure subroutine s_l0_advance_stage(s, coefs, bc_type, q_T_sf, pb_in, rhs_pb, mv_in, rhs_mv, t_step)

        integer, intent(in)                                        :: s, t_step
        real(wp), intent(in)                                       :: coefs(4)
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        type(scalar_field), intent(inout)                          :: q_T_sf
        real(stp), dimension(:,:,:,:,:), intent(inout)             :: pb_in, mv_in
        real(wp), dimension(:,:,:,:,:), intent(inout)              :: rhs_pb, rhs_mv

        call s_l0_advance_stage_rhs(s, bc_type, q_T_sf, pb_in, rhs_pb, mv_in, rhs_mv, t_step)
        call s_l0_advance_stage_rk(s, coefs)

    end subroutine s_l0_advance_stage

    !> RHS pass for all owned tiles: fill domain-edge BC + interior-seam (fine-fine) halo, then s_compute_rhs each owned tile into
    !! its per-slot rhs. Leaves amr_slots(k)%rhs ready for the RK pass (or, under coexist, for the reflux-delta copy-back first).
    impure subroutine s_l0_advance_stage_rhs(s, bc_type, q_T_sf, pb_in, rhs_pb, mv_in, rhs_mv, t_step)

        integer, intent(in)                                        :: s, t_step
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        type(scalar_field), intent(inout)                          :: q_T_sf
        real(stp), dimension(:,:,:,:,:), intent(inout)             :: pb_in, mv_in
        real(wp), dimension(:,:,:,:,:), intent(inout)              :: rhs_pb, rhs_mv
        integer                                                    :: islot
        integer(8)                                                 :: tc0, tc1, crate
        logical                                                    :: measure

        ! measure per-tile compute time only when rebalancing is active (the GPU_WAIT bracketing serialises the GPU, so it is off by
        ! default). GPU-synced wall time (cpu_time would capture only host launch overhead under offload); accumulated across
        ! stages, reset at each rebalance. Timing is a pure side-channel; it never touches field data, so output stays
        ! bit-identical.

        measure = (l0_rebalance_interval > 0)

        call s_l0_fill_edge_bc()
        call s_amr_fine_fine_halo(0)
        ! Fill the multi-dim ghost cells (2D diagonal corners; 3D also the ghost edges) that the dimension-split face fills
        ! (s_l0_fill_edge_bc + s_amr_fine_fine_halo) deliberately leave unset. The RHS never reads them, but the cons->prim
        ! convert processes the whole buffered range, and an unset ghost (all void fractions 0 -> gamma 0) is a 0/0 that traps
        ! under -ffpe-trap and is a stray NaN otherwise. Each is copied from the nearest interior cell (valid state; the face
        ! ghosts the RHS reads are single-ghost and untouched, so output is unaffected).
        do islot = 1, l0_ntiles_tot
            if (amr_block_owner(islot) /= proc_rank) cycle
            call s_l0_fill_ghost_corners(amr_loc_of(islot), amr_slots(islot)%m, amr_slots(islot)%n, amr_slots(islot)%p)
        end do
        do islot = 1, l0_ntiles_tot
            if (amr_block_owner(islot) /= proc_rank) cycle  ! advance only owned tiles; remote tiles live on their owner rank
            call s_amr_select_slot(islot)
            if (measure) then
                $:GPU_WAIT()
                call system_clock(tc0)
            end if
            ! tiles fill their own rhs (it must survive the whole-set RHS pass through the reflux point to the RK pass); the
            ! per-slot q_prim exists exactly when the copy-out gate writes it; otherwise the pooled scratch takes the (unread,
            ! unwritten) dummy
            if (allocated(amr_slots(islot)%q_prim)) then
                call s_amr_fine_stage_rhs(s, bc_type, q_T_sf, amr_slots(islot)%q_prim, amr_slots(islot)%rhs, pb_in, rhs_pb, &
                                          & mv_in, rhs_mv, t_step)
            else
                call s_amr_fine_stage_rhs(s, bc_type, q_T_sf, amr_scr_prim, amr_slots(islot)%rhs, pb_in, rhs_pb, mv_in, rhs_mv, &
                                          & t_step)
            end if
            if (measure) then
                $:GPU_WAIT()
                call system_clock(tc1, crate)
                amr_tile_cost(islot) = amr_tile_cost(islot) + real(tc1 - tc0, wp)/real(crate, wp)
            end if
        end do
        call s_amr_select_slot(1)

    end subroutine s_l0_advance_stage_rhs

    !> RK pass for all owned tiles: SSP-RK update consuming each tile's per-slot rhs (already reflux-corrected under coexist).
    impure subroutine s_l0_advance_stage_rk(s, coefs)

        integer, intent(in)  :: s
        real(wp), intent(in) :: coefs(4)
        integer              :: islot

        do islot = 1, l0_ntiles_tot
            if (amr_block_owner(islot) /= proc_rank) cycle
            call s_amr_select_slot(islot)
            if (allocated(amr_slots(islot)%q_prim)) then
                call s_amr_fine_stage_rk(s, coefs, amr_slots(islot)%q_prim, amr_slots(islot)%rhs)
            else
                call s_amr_fine_stage_rk(s, coefs, amr_scr_prim, amr_slots(islot)%rhs)
            end if
        end do
        call s_amr_select_slot(1)

    end subroutine s_l0_advance_stage_rk

    !> Free the base-grid tiling allocations (mirror of s_l0_tiles_init).
    impure subroutine s_l0_tiles_finalize()

        integer :: islot

        if (l0_ntile <= 0) return
        ! Only free the shared slot pool in pure-L0 mode: under coexist (amr) s_finalize_amr_module runs first (m_start_up call
        ! order) and already freed every slot 1..amr_max_blocks and deallocated amr_slot_live, so re-running s_amr_free_slot here
        ! would read the deallocated amr_slot_live (use-after-free).
        if (.not. amr) then
            do islot = 1, amr_max_blocks
                call s_amr_free_slot(islot)
            end do
        end if
        if (allocated(amr_seam_pairs)) deallocate (amr_seam_pairs)
        ! amr_slots, amr_region_*, amr_isect_*, amr_owns_all, amr_block_owner, amr_block_level, amr_ovl_*, and
        ! amr_slot_live are shared with s_initialize_amr_module/s_finalize_amr_module: when amr, that pair owns them, so only
        ! free them here in l0-only mode to avoid a coexist double-free. amr_tile_l0_owner/amr_tile_cost/amr_tile_cost_ema are
        ! tile-only and always freed here.
        if (.not. amr) then
            deallocate (amr_slot_live)
            call s_amr_st_finalize()
            if (allocated(amr_ovl_gather)) deallocate (amr_ovl_gather)
            if (allocated(amr_ovl_scatter)) deallocate (amr_ovl_scatter)
            deallocate (amr_ovl_gather_n, amr_ovl_scatter_n)
            if (allocated(amr_gpl_nsrc)) deallocate (amr_gpl_nsrc, amr_gpl_src, amr_gpl_sz, amr_gpl_psrc, amr_gpl_psz)
            if (allocated(amr_gcr_pool)) deallocate (amr_gcr_pool)
            if (allocated(amr_gcr_req)) deallocate (amr_gcr_req, amr_gcr_off)
            deallocate (amr_slots)
            deallocate (amr_region_lo_all, amr_region_hi_all, amr_isect_lo_all, amr_isect_hi_all, amr_owns_all)
            deallocate (amr_block_owner, amr_block_level)
            if (allocated(amr_owner_cut)) deallocate (amr_owner_cut)
            if (allocated(amr_fine_cut)) deallocate (amr_fine_cut)
        end if
        deallocate (amr_tile_l0_owner, amr_tile_cost, amr_tile_cost_ema)
        if (allocated(sw_x_cb)) deallocate (sw_x_cb, sw_x_cc, sw_dx)
        if (allocated(sw_y_cb)) deallocate (sw_y_cb, sw_y_cc, sw_dy)
        if (allocated(sw_z_cb)) deallocate (sw_z_cb, sw_z_cc, sw_dz)
        if (allocated(amr_gxcb)) deallocate (amr_gxcb)
        if (allocated(amr_gycb)) deallocate (amr_gycb)
        if (allocated(amr_gzcb)) deallocate (amr_gzcb)

    end subroutine s_l0_tiles_finalize

end module m_amr_l0
