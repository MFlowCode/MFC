!>
!!@file
!!@brief Contains module m_amr_state

#! AMD OpenMP lane: assert allocatables present on every kernel here (see OMP_DEFAULT_STR). Every conditionally allocated
#! module array a kernel here names launches only under its allocation's own condition (amr_rvw: cyl_coord; sw_jac/jac: igr;
#! amr_cg_pb/mv: do_pbmv; amr_gst_a/b: amr_subcycle; amr_prim_st/amr_bt_*: amr_prim_batch); amr_cg and amr_cons_br/stor_st are
#! allocated before first use. A kernel naming an unallocated array aborts. Keep it so.
#:set MFC_OMP_PRESENT_ALLOCATABLE = True
#:include 'macros.fpp'

!> @brief Module-scope AMR state (block pool, flat field store, wave buffers) and the leaf helpers every AMR module needs.
module m_amr_state

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

    implicit none

    !> Block/slot state and fine-distribution services consumed by the regrid and restart drivers (m_amr_regrid, m_amr_restart). The
    !! state lives in this module; those modules only drive it.

    !> Realizability floor for prolonged Euler-Euler bubble positive moments (radius nR, non-polytropic partial pressure npb / vapor
    !! mass nmv): a positive fraction of the coarse parent so derived R = nR/n, pb, mv stay >= 0. Minmod keeps a positive field
    !! positive, so this fires only under floating-point edge cases (conservation defect ~0 otherwise).
    real(wp), parameter :: bub_pos_frac = 1.0e-10_wp

    !> One refined level: its own grid + conservative fields. Field arrays are device-resident (@:ALLOCATE); coords/metadata
    !! host-only.
    type t_level
        integer               :: amr_ref_ratio
        type(t_box)           :: region   !< block extent in parent (level-0) cell indices
        integer               :: m, n, p  !< this level's interior extents
        integer               :: buff_size
        type(int_bounds_info) :: idwbuff(3)
        real(wp), allocatable :: x_cb(:), x_cc(:), dx(:)
        real(wp), allocatable :: y_cb(:), y_cc(:), dy(:)
        real(wp), allocatable :: z_cb(:), z_cc(:), dz(:)
        !> conserved state lives in the flat store amr_cons_st, indexed by this slot's dense local index amr_loc_of - not here.
        !> SSP-RK stage storage lives in the flat store amr_stor_st, indexed by this slot's amr_loc_of.
        type(scalar_field), allocatable :: q_prim(:)  !< primitive stage (fine advance, same bounds)
        type(scalar_field), allocatable :: rhs(:)     !< RHS (fine interior only: 0:m, 0:n, 0:p)
        !> non-polytropic QBMM quadrature side-state on the block (nnode x nb per cell). pb/mv evolve cell-locally (their rhs reads
        !! only the local cell + the block's own moment fluxes), so the fine treatment is prolong -> advance -> restrict with no
        !! reflux; ghosts feed the widened-idwint conversions and are prolonged piecewise-constant (CHyQMOM realizability, like the
        !! moments).
        type(pres_field) :: pb_f, mv_f        !< fine pb/mv (ghost-inclusive)
        type(pres_field) :: pb_stor, mv_stor  !< SSP-RK step-entry backup (also the regrid bounce)
    end type t_level

    !> Fixed pool of refined-block slots (at init one slot is active; dynamic regrid activates up to amr_max_blocks). The working
    !! slot amr_cur (m_global_parameters) selects which slot every per-block routine operates on.
    type(t_level), allocatable :: amr_slots(:)

    !> Dense local index for live slots. `amr_slots` is indexed by global block index and allocated lazily, so live slots are sparse
    !! across 1:amr_max_blocks. The contiguous per-block field store is indexed densely by a local index instead, so it need not be
    !! sized for the whole global pool. `amr_loc_of(g)` is the local index of global slot g (0 if not live), `amr_loc_n` is the
    !! high-water mark, and freed indices are recycled through `amr_loc_free` so the dense range stays tight under regrid churn.
    integer, allocatable :: amr_loc_of(:)    !< global slot -> dense local index, 0 if not live
    integer, allocatable :: amr_loc_free(:)  !< stack of recycled local indices
    integer              :: amr_loc_n = 0    !< high-water mark of local indices handed out
    !> Last high-water mark reported by the store trip-wire. Module scope (not `save`) so the wire prints one line per new
    !! high-water instead of one per slot allocation.
    integer :: amr_st_hw = 0
    integer :: amr_loc_nfree = 0  !< depth of the recycle stack

    !> Flat per-block field store, indexed (x, y, z, var, local slot) by the dense index above. One contiguous module array rather
    !! than a per-slot vector of independently allocated scalar_fields, so one kernel can run over every live block. Every slot's
    !! arrays carry the same mbuf extents, so a single array serves them all. Sized by s_amr_st_reserve.
    real(stp), allocatable, dimension(:,:,:,:,:) :: amr_cons_st, amr_stor_st
    $:GPU_DECLARE(create='[amr_cons_st, amr_stor_st]')
    !> local slots the store is sized for; grows and never shrinks, but plateaus at the rebuild-transient high-water because
    !! s_amr_compact_store re-densifies the index space every reconcile
    integer :: amr_st_cap = 0

    !> Prim landing zone for the batched cons->prim conversion (s_amr_convert_prim_batch): the computed prim vars only, i.e. the
    !! contiguous eqn_idx%mom%beg..eqn_idx%E range (velocities + pressure), var dim 1..num_vels+1. The aliased prim vars (cont, adv,
    !! c, psi) ride the cons copy-in inside s_compute_rhs as always. Per-stage scratch: rewritten by every batch call, so store
    !! growth discards it (no staging round trip). Sized with the store in s_amr_st_reserve, only under the amr_prim_batch gate.
    real(stp), allocatable, dimension(:,:,:,:,:) :: amr_prim_st
    $:GPU_DECLARE(create='[amr_prim_st]')
    !> per-dense-slot batch metadata: participating fine slot + its idwbuff window (host-filled each batch call)
    integer, allocatable :: amr_bt_lo(:,:), amr_bt_hi(:,:)  !< (3, loc)
    logical, allocatable :: amr_bt_on(:)
    $:GPU_DECLARE(create='[amr_bt_lo, amr_bt_hi, amr_bt_on]')
    !> Batched-conversion gate, derived once at init: the batched conversion covers the plain multi-fluid configs (5/6-eq, WENO,
    !! with/without viscous); every feature that adds conversion write-set members or changes its inputs (igr, chemistry,
    !! relativity, hypoelasticity, mhd, cont_damage, ib, Lagrangian bubbles) uses the per-block conversion path instead.
    logical :: amr_prim_batch = .false.

    !> Copy bridge to the shared solver. s_compute_rhs, s_ibm_correct_state, s_pressure_relaxation_procedure and
    !! s_infinite_relaxation_k all take type(scalar_field), dimension(sys_size) and serve the monolithic path too, so the flat store
    !! cannot be handed to them, and a pointer view into the store is not attachable on the OpenMP-offload backend. One block-shaped
    !! scalar_field array bridges instead: load it from the store, call, store it back. All four dummies are intent(inout)
    !! (s_compute_rhs writes the buffer region through s_populate_variables_buffers), so both directions are required at every
    !! crossing.
    !> Device-declared because `@:ALLOCATE` expands to `allocate` followed by `GPU_ENTER_DATA(create=)`, i.e. an `omp target enter
    !! data map` on the variable itself. Mapping a module allocatable that was never `declare target` has no device descriptor to
    !! attach to: Cray CCE aborts at runtime with `lib-4425 UNRECOVERABLE library error: Unitialized descriptor for ALLOCATE
    !! statement argument`, while amdflang's runtime tolerates it and creates the mapping implicitly. The fix is the `move_alloc` +
    !! GPU_ENTER_DATA pair at the allocation site; do not add a GPU_DECLARE on top of it (see amr_scr_prim below).
    type(scalar_field), allocatable :: amr_cons_br(:)
    !> Batched advance (amr_batched_advance): the bridge spans amr_br_batch blocks along the last active dimension; member i of a
    !! batch sits at offset (i-1)*amr_bat_w there (amr_bat_w = the batch's block width + two ghost shells, m_global_parameters),
    !! carrying its own ghost shell, so consecutive blocks are separated by two ghost shells and no block's stencil can reach
    !! another's interior (the property the batched advance's correctness rests on). amr_bat_loc = the members' flat-store columns,
    !! device-resident so the batch kernels index the store without a per-launch map.
    integer :: amr_bat_loc(amr_bat_max) = 0
    $:GPU_DECLARE(create='[amr_bat_loc]')
    !> Batched-advance population histogram: run-lifetime count of the batches formed, indexed by member count. Reported at finalize
    !! as [amr-bat] under rank_time_wrt. hist(1) is the single-member count; a deck whose batches are all single-member exercises
    !! the batching frame but not the stacking.
    integer :: amr_bat_hist(amr_bat_max) = 0
    !> Per-batch timing log (rank_time_wrt): one line per batch per stage (step, stage, members, level, extents, cells per member,
    !! seconds in swap / rhs / rk, then the members' block ids and Morton keys) to amr_batch_r<rank>.log.
    integer :: amr_bat_unit = -1
    logical :: amr_bat_open = .false.  ! newunit= hands back a negative unit, so the unit's sign cannot serve as the sentinel
    !> Pooled advance scratch: the fused per-block fine advance (rhs then rk on one block, s_amr_fine_stage_advance) leaves no
    !! cross-block q_prim/rhs lifetime, so every fine block shares this one slot-shaped pair instead of carrying per-slot arrays
    !! (which would multiply the live footprint by the slot count and churn the device allocator). Same shared-scratch pattern as
    !! amr_rhs_pb_f/amr_cg. L0 tile slots are the exception and keep per-slot arrays (see s_amr_alloc_slot).
    type(scalar_field), allocatable :: amr_scr_prim(:), amr_scr_rhs(:)
    !> block-frame primitive scratch for the batched advance's per-member IB correction (allocated only with ib): the slab's prim
    !! holds the members stacked along amr_bat_sd, and s_ibm_correct_state reads a block in its own frame
    type(scalar_field), allocatable :: amr_scr_prim_blk(:)
    !> Not device-declared. A GPU_DECLARE(create=) on a module allocatable binds a present-table entry to the descriptor at program
    !! init; the `move_alloc` at the allocation site then swaps that descriptor out and every later kernel lookup misses (Cray CCE
    !! gpu-acc fails with `find_in_present_table failed`). The `move_alloc` + GPU_ENTER_DATA pair at the allocation site is what
    !! fixes the lib-4425 descriptor abort, and it is sufficient on its own; `amr_cg` uses exactly that shape with no declare.
    !> True only while the regrid path is inside s_amr_gather_coarse_patch, so the WAITALL bracket attributes to rb:wait rather than
    !! mixing in the per-step gather that shares this routine.
    !> Blocks per batched s_compute_rhs call: amr_bat_max under amr_batched_advance, else 1 (the bridge holds one block). Bounded on
    !! purpose: sizing the bridge per live block can exhaust device memory on large cases.
    integer :: amr_br_batch = 1
    integer :: amr_maxc(3)  !< max coarse block cells per dim: (m_glb+1)/2 etc.; 1 for collapsed dims

    !> Per-slot field-array sizing (module-scope, used by s_amr_alloc_slot/s_amr_free_slot): max fine cells per dim (2*maxc_loc-1)
    !! and the buffered array bounds. amr_slot_live(k) tracks whether slot k's field arrays are allocated - lazy owned-only sizing
    !! keeps a rank's fine memory ~1/num_procs of the pool.
    integer              :: max_f1, max_f2, max_f3
    integer              :: mbuf1_lo, mbuf1_hi, mbuf2_lo, mbuf2_hi, mbuf3_lo, mbuf3_hi
    logical, allocatable :: amr_slot_live(:)
    !! cached same-level adjacent-seam list (3, npairs) = (xb, yb, seam-dim), so s_amr_fine_fine_halo iterates O(#seams) instead of
    !! rescanning all O(nblocks^2) pairs every RK stage. Block topology changes only at regrid/restart, so the list is rebuilt only
    !! when amr_seam_pairs_dirty is set (or the block count changes - a tripwire).
    integer, allocatable :: amr_seam_pairs(:,:)
    integer              :: amr_num_seam_pairs, amr_seam_pairs_nblk
    logical              :: amr_seam_pairs_dirty
    !! amr_mesh_epoch lives in m_global_parameters (m_amr_registers keys its participation-map rebuild on it and cannot use
    !! m_amr); it is use-associated here and re-exported.
    !! Per-family message tag bases: amr_max_blocks + 100*f keeps this tag space disjoint from the per-box space (tags in
    !! [1..amr_max_blocks]). The wave families use the keyed tags (amr_m1_base bands below); the only user of this array is the
    !! regrid's per-box migration, amr_tag_base(4) + mod(amr_mesh_epoch, 50) in m_amr_regrid.fpp (entries 1..3, 5..7 are unused).
    !! The init MPI_TAG_UB assert is the scale tripwire; Open MPI reports 2**31 - 1 and Cray MPICH 2**29 - 1, so the tag space is
    !! not a scaling limit on either. The amr_max_blocks term can only go once no site tags per box (the migration still does).
    integer :: amr_tag_base(7) = 0
    !> Keyed wave tags: tag = amr_m1_base + band*65536 + gen*4096 + seq, checked against MPI_TAG_UB at init. band 0 = reflux-faces
    !! wave, 1 = freg wave, 2 = parent-fill wave (F2W), 3/4 = stage-fill q / pb-mv waves (F1W/F3W), 5 = fine-fine halo wave (F6W), 6
    !! = level-1 restrict wave (F7W), 7 = parent restrict wave (F7BW). gen (mod 16) bumps at wave entry on every rank (every wave
    !! call site is rank-unconditional), separating successive waves that share a band; seq is the message's position in the pair's
    !! canonically ordered transfer list (ascending block id, then dim, lo before hi), derived independently by each end from
    !! replicated metadata, so message matching does not depend on posting order.
    integer              :: amr_m1_base = 0
    integer              :: amr_tag_gen(0:7) = 0
    integer, allocatable :: amr_tsq(:,:)      !< (0:np-1, dir) in-wave per-peer seq counters; touched-reset
    integer, allocatable :: amr_tsq_tch(:,:)  !< touched peers per dir, so the reset is O(active peers), not O(P)
    integer              :: amr_n_tsq(2) = 0
    !! cached per-block P2P overlap-rank lists (rebuilt with the seam list - same dirty flag): amr_ovl_gather(:,k) = ranks whose
    !! owned coarse range (s_amr_rank_coarse_range) intersects block k's amr_cpat_mar-padded patch box (gather contributors);
    !! amr_ovl_scatter(:,k) = ranks whose coarse interior (s_amr_rank_interior) intersects block k's region box (restrict-scatter
    !! destinations). Built by O(overlap) inversion (s_amr_ranks_overlapping), rank-ascending and not owner-excluded (consumers
    !! keep their owner skip), so iterating a list gives the same MPI send/recv order as a 0..num_procs-1 scan.
    !> (max-overlap, amr_max_blocks); sized in s_amr_build_seam_pairs
    integer, allocatable :: amr_ovl_gather(:,:), amr_ovl_scatter(:,:)
    !> The block indices this rank owns, ascending. Iterating this list is O(local) where a scan of the global block list filtered
    !! by ownership would be O(global blocks) per stage. Built ascending in s_amr_assign_block_owners, the single authority for
    !! ownership (regrid, init and both restart paths all route through it), so iterating it visits blocks in the same order as an
    !! ascending global scan, which the paired MPI_SENDRECVs depend on.
    integer, allocatable :: amr_my_blk(:)
    integer              :: amr_n_my = 0
    !> The level-1 blocks this rank receives restriction for: not its own, and overlapping its interior. The predicate is level +
    !! single-owner + geometry, all fixed between regrids, so the list is rebuilt once per mesh epoch and the per-stage loop walks
    !! it instead of every block. Keyed on `amr_mesh_epoch`, not `amr_myblk_dirty`: that flag tracks owner writes only, and this
    !! predicate also depends on level and region, which m_amr_regrid.fpp/m_amr_restart.fpp change without touching it. Note the
    !! ownership notion: `amr_block_owner(k) == proc_rank` (single owner), which is not the same as `amr_rank_owns_block` (the
    !! multi-owner intersection). Loops testing the latter need their own list.
    integer, allocatable :: amr_l1r_blk(:)
    integer, allocatable :: amr_l1p_blk(:)  !< padded variant: region +/- amr_cpat_mar vs my COARSE range
    integer              :: amr_n_l1p = 0   !< (the stage-fill gather predicate; superset differs from amr_l1r)
    integer              :: amr_n_l1r = 0
    integer, allocatable :: amr_fch_blk(:)  !< foreign children of my owned parents (level >= 2, powner == me, cowner /= me):
    integer              :: amr_n_fch = 0   !! the freg-recv / parent-fill-send / restrict-parent-recv survivor superset
    integer, allocatable :: amr_own_blk(:)  !< blocks this rank intersects (amr_owns_all, the multi-owner notion, all levels)
    integer              :: amr_n_own = 0
    !> cached f_amr_parent_block (0 for level <= 1); the function is itself an O(global blocks) scan, so per-stage wave bodies
    !! calling it per block would be quadratic in the global block count
    integer, allocatable :: amr_parent_blk(:)
    integer, allocatable :: amr_child_ptr(:)   !< children of p = amr_child_idx(amr_child_ptr(p-1)+1 : amr_child_ptr(p)),
    integer, allocatable :: amr_child_idx(:)   !! ascending; amr_child_ptr(0) = 0
    integer(8)           :: amr_l1r_epoch = -1_8
    integer              :: amr_l1r_nblk = -1  !< second key half, matching s_amr_reg_prepare's (epoch, num_blocks) pair
    !> Set wherever amr_block_owner is written. s_amr_assign_block_owners is not the only writer: a tiled level-2 block inherits its
    !! parent's owner (s_amr_add_l2_tile), and the restart/migration paths assign directly. A list built only in the assigner would
    !! go stale and a loop over it would visit the wrong blocks.
    logical              :: amr_myblk_dirty = .true.
    integer, allocatable :: amr_ovl_gather_n(:), amr_ovl_scatter_n(:)  !< per-block list lengths
    !> Rebuild gather plan: the whole rebuild's gather message set, derived up front by s_amr_build_gather_plan from the replicated
    !! caches. Per level-1 slot: contributor count/ranks/message sizes (owner excluded, list order = amr_ovl_gather order). Per
    !! level>=2 slot: the parent-owner source rank (-1 when co-located, no message) and its message size. The chunked path below
    !! checks, guarded on amr_gpl_valid, that the plan reproduces the message set.
    integer, allocatable :: amr_gpl_nsrc(:), amr_gpl_src(:,:), amr_gpl_sz(:,:), amr_gpl_psrc(:), amr_gpl_psz(:)
    logical              :: amr_gpl_valid = .false.  !< true only between plan build and the end of the rebuild box loop
    !> The rebuild's participant list: the ascending union of amr_my_blk (owner: posts, consumes), amr_fch_blk (owner of a foreign
    !! child's parent: the level>=2 send) and amr_l1p_blk (level-1 contributor: the send phase and the pb/mv gather), fine band
    !! only. The consumers keep their per-box predicates; the list only drops boxes they would have cycled. A box this rank has no
    !! role in touches nothing of its own but the replicated non-owner geometry, which the rebuild fills in one plain pass. Built by
    !! s_amr_build_gather_plan from the epoch-keyed lists; valid exactly as long as amr_gpl_valid.
    integer, allocatable :: amr_gpk(:)
    integer              :: amr_n_gpk = 0
    !> Rebuild walk order: amr_korder(p) is the box visited at position p, amr_kpos(k) its inverse. Level-major (so parents-first
    !! holds), and inside a level round-robin over owners: the Morton cut makes the box id monotone in owner, so an ascending walk
    !! would give every chunk to one rank and the rebuild would run rank after rank. A pure function of replicated metadata, so
    !! every rank derives the same order. amr_korder_rot = .false. is the identity walk.
    integer, allocatable :: amr_korder(:), amr_kpos(:)
    logical, parameter   :: amr_korder_rot = .true.
    !> Chunked rebuild gather: the rebuild box loop runs in chunks of amr_gath_chunk boxes. Every owned box's recvs (level-1
    !! contributor slices and split level>=2 parent patches) are pre-posted from the plan into one flat pool, this rank's sends are
    !! issued (level>=2 only when the parent was consumed in an earlier chunk; a same-chunk parent's store is unbuilt until its own
    !! consume, so that send stays at the child's consume position), then boxes are consumed in order with a per-box wait. Requests
    !! are appended in box order, so each box's recvs are the contiguous run amr_gcr_r0 : +amr_gcr_nr-1. The pool/request arrays
    !! grow monotonically and are reusable across chunks only because the consume phase waits every owned box's requests
    !! unconditionally inside its own chunk.
    integer, parameter    :: amr_gath_chunk = 32  !< boxes per chunk: staging memory vs message batching
    real(wp), allocatable :: amr_gcr_pool(:)  !< flat recv staging for one chunk
    integer, allocatable  :: amr_gcr_req(:), amr_gcr_off(:)  !< request handle + pool offset per posted recv
    integer               :: amr_gcr_r0(amr_gath_chunk), amr_gcr_nr(amr_gath_chunk)  !< per chunk-local box: first recv, count
    logical               :: amr_gcr_sent(amr_gath_chunk)  !< chunk-local: level>=2 send already issued in the send phase
    integer               :: amr_gcr_n = 0  !< posted recvs in the current chunk
    !> Stage-fill wave (plan-based exchange): the level-1 per-stage fill's F1 q_cons + F3 pb/mv gathers as one per-(peer, family)
    !! aggregated exchange per RK stage. Transfer records are SoA flat arrays (no derived types); the wire layout of each peer
    !! message is the ascending-box concatenation of [XA_NH header | slab], which sender and receiver derive independently from the
    !! replicated caches (rank coarse ranges x patch boxes), so no metadata is exchanged. Plans are rebuilt every wave. All scratch
    !! is high-water and its contents never survive a wave; the rank-indexed build counters (amr_fw_map/nx/pq/pp) are re-zeroed for
    !! touched ranks after each build so they stay all-zero between builds.
    integer               :: amr_fw_snx = 0, amr_fw_rnx = 0, amr_fw_snp = 0, amr_fw_rnp = 0
    integer, allocatable  :: amr_fw_sblk(:), amr_fw_sbl(:,:), amr_fw_sbh(:,:), amr_fw_spi(:), amr_fw_sqo(:), amr_fw_spo(:)
    integer, allocatable  :: amr_fw_rblk(:), amr_fw_rbl(:,:), amr_fw_rbh(:,:), amr_fw_rpi(:), amr_fw_rqo(:), amr_fw_rpo(:)
    integer, allocatable  :: amr_fw_sprank(:), amr_fw_sqsz(:), amr_fw_spsz(:), amr_fw_snxp(:), amr_fw_sqbase(:), amr_fw_spbase(:)
    integer, allocatable  :: amr_fw_rprank(:), amr_fw_rqsz(:), amr_fw_rpsz(:), amr_fw_rnxp(:), amr_fw_rqbase(:), amr_fw_rpbase(:)
    integer, allocatable  :: amr_fw_map(:), amr_fw_nx(:), amr_fw_pq(:), amr_fw_pp(:)  !< rank-indexed build scratch (0:num_procs-1)
    real(wp), allocatable :: amr_fw_sq(:), amr_fw_sp(:), amr_fw_rq(:), amr_fw_rp(:)   !< wire pools (live across the ISENDs)
    !> Device-resident wire pools: with rdma_mpi the four pools live on the device and MPI sends and receives them by device
    !! address, as the base halo does. Every step-path writer and reader of a pool is a device kernel whose copyin/copyout of the
    !! slice then finds the pool present and copies nothing; otherwise each box's slice would cross PCIe twice per wave. Off when
    !! the exchange audit writes host headers.
    logical              :: amr_fw_dev = .false.
    integer, allocatable :: amr_fw_req(:), amr_fw_reqw(:)  !< requests + expected recv word counts (-1 for sends; debug check)
    !> Seam wave's private pools: the seam is posted at the top of the stage and drained after the parent fills, so it must not
    !! share the wave scratch the gather/parent waves rebuild in between. Same layout as amr_fw_*; the rank-indexed build scratch
    !! (amr_fw_map/nx/pq/pp) stays shared because plan builds never overlap.
    integer               :: amr_sw_snx = 0, amr_sw_rnx = 0, amr_sw_snp = 0, amr_sw_rnp = 0, amr_sw_nreq = 0, amr_sw_nsame = 0
    integer, allocatable  :: amr_sw_sblk(:), amr_sw_sbl(:,:), amr_sw_spi(:), amr_sw_sqo(:), amr_sw_spo(:)
    integer, allocatable  :: amr_sw_rblk(:), amr_sw_rbl(:,:), amr_sw_rbh(:,:), amr_sw_rpi(:), amr_sw_rqo(:), amr_sw_rpo(:)
    integer, allocatable  :: amr_sw_sprank(:), amr_sw_sqsz(:), amr_sw_snxp(:), amr_sw_sqbase(:)
    integer, allocatable  :: amr_sw_rprank(:), amr_sw_rqsz(:), amr_sw_rnxp(:), amr_sw_rqbase(:)
    real(wp), allocatable :: amr_sw_sq(:), amr_sw_rq(:)
    integer, allocatable  :: amr_sw_req(:), amr_sw_reqw(:)
    integer, allocatable  :: amr_sw_plx(:), amr_sw_ply(:), amr_sw_pd(:), amr_sw_pxhi(:), amr_sw_pfm(:,:)  !< same-rank pairs
    logical, parameter    :: amr_early_seam_post = .true.
    !> Fused exchange packs (amr_device_pack): one row per wave transfer (slab corner (1:3), slab extents (4:6), the transfer's
    !! absolute payload offset in the wire pool (7), and for the F2 pack the source store slot (8) and the child's patch frame
    !! (9:11)) plus the exclusive element prefix, so one kernel walks a whole family's transfer list by flat index instead of one
    !! launch per transfer. Rebuilt per wave from the amr_fw_* tables; contents never survive a wave.
    integer, allocatable :: amr_fx_pl(:,:), amr_fx_pre(:)
    !> [amr-cad] regrid-cadence containment audit: level-1 tags counted at each regrid, and how many fell outside the pre-regrid
    !! level-1 coverage (a feature that evolved unrefined since the last regrid, i.e. the tag buffer amr_buf did not cover its
    !! drift). Zero escaped validates the (amr_regrid_int, amr_buf) pair for that run; the case validator only warns (the CFL <= 1
    !! worst case is too strict for low-CFL cases). Incremented by m_amr_regrid, reported by s_amr_cad_report.
    integer(8) :: amr_cad_tot = 0, amr_cad_esc = 0
    logical    :: amr_cad_armed = .false.  !< first regrid (hierarchy population) is skipped, see s_amr_cad_count
    !> (0:num_procs-1) SFC Morton-key upper bound per rank from the cost-weighted split; owner = cut-search (f_amr_owner). O(P)
    !! storage, where the amr_block_owner table is O(global blocks).
    integer(kind=8), allocatable :: amr_owner_cut(:)
    !> (0:num_procs-1, 1:amr_max_level) companion cuts for fine-block (level>=1) owners: one independent cut per level. Each level's
    !! boxes are balanced across all ranks on their own weight, so a deep refinement tower does not pin its whole subtree to one
    !! rank. Per level rather than one mixed cut because same-level boxes are disjoint and so have distinct Morton keys, which the
    !! cut-point binary search requires; a mixed cut would let a child share its parent's region_lo and make the search ambiguous.
    !! In no-tile AMR level 1's cut also mirrors amr_owner_cut, but in coexist amr_owner_cut is overwritten by the tile cut, so the
    !! fine cuts are kept here for f_amr_owner (fine blocks straddle tiles and cannot be derived from the tile cut).
    integer(kind=8), allocatable :: amr_fine_cut(:,:)

    !> Regrid box size cap per dim (fixed for the run, identical on all ranks; 1 in collapsed dims): a box of at most min-over-ranks
    !! of (local extent + 1)/2 cells intersects every rank in at most (its extent + 1)/2 cells, so the per-rank scratch constraint
    !! 2*(isect cells) - 1 <= local extent holds by construction. Equals amr_maxc at np=1.
    integer :: amr_maxc_fit(3) = 1

    !> Swap contract (s_amr_swap_to_fine / s_amr_restore_coarse). The fine advance runs the shared solver on a fine block by
    !! swapping these coarse-grid globals to the block's values and restoring after: m/n/p, idwint/idwbuff, the nine coordinate
    !! arrays (sw_x_cb..sw_dz below), acoustic_source, ab_active; WENO/hypoelastic/IGR spacing coefficients are recomputed for the
    !! active grid rather than saved. Rule for anyone adding grid-dependent state: any module-level variable derived from
    !! m/n/p/idwint/idwbuff/coords that a kernel reads on the fine grid must be swapped here or refreshed on every fine call at its
    !! use site, and if GPU_DECLARE'd, its device copy too. A stale device copy of coarse bounds reads out of range on the fine grid
    !! under CCE OpenACC (ab_int is refreshed by a per-call GPU_UPDATE in s_compute_rhs for this reason; see m_rhs.fpp and
    !! .claude/rules/common-pitfalls.md). amr_swap_depth makes the swap re-entrant and guards against an unpaired restore.
    !> Saved coarse-level global state for swap/restore
    integer               :: sw_m, sw_n, sw_p
    type(int_bounds_info) :: sw_idwint(3), sw_idwbuff(3)
    logical               :: sw_acoustic_source
    logical               :: sw_ab_active

    !> IGR sigma-state bounce (igr only): the fine solve reuses the module jac/jac_old arrays at fine indices (the extent guard
    !! keeps fine bounds inside), so the coarse contents - jac_old is the Jacobi warm start persisting across steps - are saved here
    !! across the fine advance.
    real(wp), allocatable :: sw_jac(:,:,:), sw_jac_old(:,:,:)
    $:GPU_DECLARE(create='[sw_jac, sw_jac_old]')

    !> Per-fine-cell radial volume weight for cyl_coord restriction (axisymmetric): the fold-back must be volume-weighted and cell
    !! volume ~ radius, so a fine child is weighted by its own cell-center radius y_cc. Filled from the active block's fine y_cc
    !! each restriction, read identically by the device kernel and host scatter path so np=1 == np>=2 stays element-exact. Allocated
    !! only for cyl_coord. Its device copy is refreshed (GPU_UPDATE) only in s_restrict_fine_to_coarse; s_amr_restrict_to_parent
    !! reads it without refreshing, which is safe only because m_checker.fpp forbids cyl_coord with amr_max_level > 1. Lifting that
    !! gate makes this a stale-device bug (the parent fold needs a fresh per-block radius table). See the swap contract note above.
    real(wp), allocatable :: amr_rvw(:)
    $:GPU_DECLARE(create='[amr_rvw]')

    !> Non-polytropic QBMM fine rhs scratch, shared across slots (slots advance sequentially). Module-level raw arrays mirror the
    !! coarse rhs_pb/rhs_mv pattern: derived-type component actuals here tripped nvfortran's component-section data clauses on
    !! device.
    real(wp), allocatable :: amr_rhs_pb_f(:,:,:,:,:), amr_rhs_mv_f(:,:,:,:,:)
    $:GPU_DECLARE(create='[amr_rhs_pb_f, amr_rhs_mv_f]')

    !> Swap nesting depth. Only the outermost swap saves the coarse state into the sw_* bounce buffers, and only its matching
    !! restore puts it back; an inner swap re-installs and an inner restore is a no-op. Every nested swap site swaps to the same
    !! slot (amr_cur), so the enclosing frame's view is unchanged either way.
    integer :: amr_swap_depth = 0
    !> True when the coarse grid is nonuniform (stretched grids, or 2D-axisymmetric's half-width axis cell): the spacing-dependent
    !! WENO coefficients are then recomputed for the active grid on every block swap (the fine block's grid is itself nonuniform
    !! under stretching) and restored after. False on fully uniform grids, where the recompute is skipped.
    logical :: amr_weno_coef_recompute = .false.
    logical :: amr_grid_stretched = .false.  !< stretched coarse spacing (beyond the axisym axis half-cell; set at init)
    !> Persistent global coarse cell-boundary arrays (indices -1:X_glb), assembled once at init. The fine-distribution owner
    !! reconstructs whole-block fine coordinates from these (its fine cells cover coarse cells it does not own the coordinate slice
    !! for). Exact on any grid.
    real(wp), allocatable :: amr_gxcb(:), amr_gycb(:), amr_gzcb(:)
    real(wp), allocatable :: sw_x_cb(:), sw_x_cc(:), sw_dx(:)
    real(wp), allocatable :: sw_y_cb(:), sw_y_cc(:), sw_dy(:)
    real(wp), allocatable :: sw_z_cb(:), sw_z_cc(:), sw_dz(:)

    !> Conservation-defect baselines (level-0 interior integrals at init; per-fluid masses + energy)

    !> True (identically on all ranks) iff some rank's fine ghost-fill stencil reads its coarse ghost cells: the solver populates
    !! only prim ghosts, so the cons ghosts the fill prolongs from must be halo-exchanged first. Never true at np=1 (block faces sit
    !! >= buff_size inside the domain).
    logical :: amr_xchg_coarse_ghosts = .false.
    !> local (un-reduced) accumulator behind amr_xchg_coarse_ghosts. s_set_amr_fine_geometry ORs each block's answer in here and
    !! s_amr_reduce_xchg_flag performs one allreduce for the whole scan, rather than one collective per block.
    integer :: amr_xchg_bad = 0

    !> Per-block gathered coarse patch (fine-level distribution). The block owner may not hold the coarse cells its block refines,
    !! so before each prolongation/ghost-fill the coarse patch spanning region_lo-amr_cpat_mar : region_hi+amr_cpat_mar (the full
    !! coarse-cell reach of every prolongation stencil) is gathered here point-to-point from the coarse owners
    !! (s_amr_gather_coarse_patch). Stored in amr_cg as stp scalar_fields (a drop-in for the coarse q_cons in the prolong/ghost-fill
    !! kernels) in a block-local frame: amr_cg cell 0 is global coarse cell amr_cpat_off(d). Messages carry wp, cast to stp
    !! (identity for stp coarse), so at np=1 (owner copies its own coarse) the patch equals the local coarse read bit-for-bit. Sized
    !! to the largest block. amr_slab_tab is one device-resident slab table for the shell/ghost kernels: rows
    !! sb1,se1,sb2,se2,sb3,se3,soff,scnt over <= 6 slabs, refreshed by one GPU_UPDATE per launch.
    integer, allocatable            :: amr_slab_tab(:,:)
    type(scalar_field), allocatable :: amr_cg(:)
    integer                         :: amr_cpat_mar = 0    !< coarse-cell stencil reach = (buff_size+1)/2 + 1 (matches nmar)
    integer                         :: amr_cpat_hi(3) = 0  !< amr_cg upper local bounds per dim (0 in collapsed dims)

    !> Deferred-send pool for the per-box coarse-patch gather. The gather is called once per box; the non-owner side's sends are
    !! non-blocking and completed in batches, so a rank that only contributes data can run ahead instead of rendezvousing with the
    !! owner on each box. The pool owns the buffers because MPI_ISEND requires them to stay live until completion.
    integer, parameter    :: amr_gsnd_max = 64    !< pending sends before a forced drain (bounds pool memory)
    real(wp), allocatable :: amr_gsnd_pool(:,:)
    integer, allocatable  :: amr_gsnd_req(:)
    integer               :: amr_gsnd_n = 0
    integer               :: amr_cpat_off(3) = 0  !< global coarse index of amr_cg local cell 0 (region_lo - amr_cpat_mar)
    !> Gathered coarse pb/mv patch for non-polytropic QBMM (analogue of amr_cg): the block's coarse-side pb/mv side-state,
    !! P2P-gathered from the coarse-cell owners into the block owner in the amr_cg patch-local frame (cell 0 == amr_cpat_off). Read
    !! by the pb/mv prolong + ghost-fill so np>=2 couples to the correct coarse rank. Allocated only for non-polytropic QBMM.
    real(stp), allocatable, dimension(:,:,:,:,:) :: amr_cg_pb, amr_cg_mv
    $:GPU_DECLARE(create='[amr_cg_pb, amr_cg_mv]')

    !> Level-0 tiling (l0_ntile > 0, amr off): tile the base grid into l0_ntile**num_dims base-resolution (refinement-ratio-1)
    !! blocks and advance each through the same swap-based per-block solver the AMR fine overlay uses, with tile-tile same-level
    !! seam halos (s_amr_fine_fine_halo, fmul=1) at interior faces and the physical BC at domain-edge faces. l0_ntile>0 must be
    !! byte-identical to l0_ntile=0 (monolithic). Reuses the full amr_slots/region/owner/seam machinery; the tiles are level-1
    !! blocks with amr_ref_ratio 1. Off (l0_ntile=0) => no effect.
    integer :: l0_ntiles_tot = 0  !< total tiles = l0_ntile**num_dims (0 when off)
    integer :: l0_nt(3) = 1       !< tiles per dim (1 in collapsed dims)
    ! Unification pool layout (amr_max_fine, l0_slot_off) lives in m_global_parameters beside amr_num_blocks/amr_max_blocks, so
    ! m_amr_regrid (a separate module) sees it too.
    !> per-dim global periodicity, allreduced from periodic_bc in s_l0_tiles_init (periodic_bc is set on rank 0 only, since
    !! s_read_input_file is rank-0-guarded, so it must be made consistent for the wrap-seam decision in f_amr_seam /
    !! s_l0_edge_bc_tile, which every rank must agree on)
    logical :: l0_periodic(3) = .false.
    !> tiles are persistent: L0 seeds them once (first timestep); set at init, cleared after fill
    logical :: l0_tiles_need_fill = .false.
    !> rank owning each tile's L0 storage cells (fixed = init owner); scatter routes the compute-owner's interior back to this rank
    !! when a tile has migrated (owner != l0_owner)
    integer, allocatable :: amr_tile_l0_owner(:)
    !> accumulated measured compute time per owned tile since the last rebalance (GPU-synced wall time); allreduced to a replicated
    !! cost vector that drives s_l0_rebalance, then reset
    real(wp), allocatable :: amr_tile_cost(:)
    !> per-tile exponential moving average of the (replicated) measured cost, smoothed across rebalance windows so GPU
    !! launch-latency timing noise does not drive churn
    real(wp), allocatable :: amr_tile_cost_ema(:)

contains

    !> Wall clock for the AMR instruments (batch and regrid timings); 0 without MPI, where the instruments are not used. Keeps the
    !! serial (no-MPI) build compiling.
    impure function f_amr_wtime() result(t)

        real(wp) :: t

#ifdef MFC_MPI
        t = MPI_Wtime()
#else
        t = 0._wp
#endif

    end function f_amr_wtime

    !> Fill level-1 fcb/fcc/fdx by bisecting parent cells; pcb_lb is lbound(parent_cb, 1). Passing pcb as assumed-shape resets
    !! lbound to 1; pcb_lb + idx_offset recovers original indexing. Arrays preallocated at max size; only 0..nfine filled.
    subroutine s_build_level_coords(pcb, pcb_lb, lo, nfine, fcb, fcc, fdx)

        real(wp), intent(in)                 :: pcb(:)
        integer, intent(in)                  :: pcb_lb, lo, nfine
        real(wp), allocatable, intent(inout) :: fcb(:), fcc(:), fdx(:)
        integer                              :: fi, c, idx_offset, k, rr
        real(wp)                             :: xl, xr
        ! pcb(k) = parent_cb(k + pcb_lb - 1); to access parent_cb(j): k = j - pcb_lb + 1

        rr = amr_slots(amr_cur)%amr_ref_ratio
        idx_offset = 1 - pcb_lb
        ! fine cell fi (0..nfine) subdivides coarse cell c = lo + fi/rr into rr equal parts
        fcb(-1) = pcb(lo - 1 + idx_offset)  ! left boundary of the fine region
        do fi = 0, nfine
            c = lo + fi/rr
            xl = pcb(c - 1 + idx_offset)  ! left boundary of coarse cell c
            xr = pcb(c + idx_offset)  ! right boundary of coarse cell c
            k = mod(fi, rr)  ! fi >= 0, so mod gives the sub-position in [0, rr-1]
            if (k == rr - 1) then
                fcb(fi) = xr  ! right edge of parent cell c
            else
                fcb(fi) = (real(rr - 1 - k, wp)*xl + real(k + 1, wp)*xr)/real(rr, wp)
            end if
        end do
        do fi = 0, nfine
            fdx(fi) = fcb(fi) - fcb(fi - 1)
            fcc(fi) = 0.5_wp*(fcb(fi - 1) + fcb(fi))
        end do

    end subroutine s_build_level_coords

    !> Fine cell coordinates of block k in dimension d, rebuilt from the global coarse boundaries gcb by replaying k's ancestor
    !! chain, touching no other block's slot arrays. A level-l block's grid is l nested midpoint subdivisions of the L0 boundaries,
    !! and every box in the chain is known from replicated metadata (amr_region_*_all + the global amr_ref_ratio), so any rank can
    !! reproduce it. Bit-identical to bisecting the parent's stored coords, because that array is itself the same subdivision.
    !!
    !! Do not bisect amr_slots(parent)%x_cb instead: the parent's slot is allocated only on the parent's owner, and under
    !! per-level distribution a level>=2 block can be owned by a rank holding no part of its parent. Bisecting the unallocated
    !! array there gives garbage cell widths and NaNs a few steps later.
    impure subroutine s_amr_build_block_coords(k, gcb, fcb, fcc, fdx, d)

        integer, intent(in)                  :: k, d
        real(wp), intent(in)                 :: gcb(:)  !< global L0 cell boundaries, lbound -1
        real(wp), allocatable, intent(inout) :: fcb(:), fcc(:), fdx(:)
        integer                              :: chain(0:amr_max_level), lev, j, a, lo, nf, span, rr, plo(3), phi(3)
        real(wp), allocatable                :: cur(:), scb(:), scc(:), sdx(:)

        rr = amr_ref_ratio
        lev = amr_block_level(k)
        a = k
        do j = lev, 1, -1  ! chain(j) = k's ancestor at level j; chain(lev) = k
            chain(j) = a
            if (j > 1) a = f_amr_parent_block(a)
        end do

        cur = gcb  ! level 0: the global coarse boundaries (allocatable assignment carries lbound -1)
        do j = 1, lev
            a = chain(j)
            span = amr_region_hi_all(d, a) - amr_region_lo_all(d, a) + 1  ! L0 cells the box covers
            nf = rr**j*span - 1  ! its fine extent at level j
            if (j == 1) then
                lo = amr_region_lo_all(d, a)  ! global L0 index of the box's low corner
            else
                call s_amr_parent_foot(a, chain(j - 1), plo, phi)  ! low corner in the parent's fine frame
                lo = plo(d)
            end if
            if (j == lev) then
                call s_build_level_coords(cur, -1, lo, nf, fcb, fcc, fdx)
            else
                if (allocated(scb)) deallocate (scb, scc, sdx)
                allocate (scb(-1:nf), scc(0:nf), sdx(0:nf))
                call s_build_level_coords(cur, -1, lo, nf, scb, scc, sdx)
                cur = scb
            end if
        end do
        if (allocated(scb)) deallocate (scb, scc, sdx)

    end subroutine s_amr_build_block_coords

    !> Compute this rank's per-dim intersection of the box lo:hi with its subdomain (global indices, mirrored to amr_isect_lo/hi)
    !! and whether it holds fine cells (amr_rank_owns_block: nonempty in all active dims). Must be called with the coarse grid state
    !! in m/n/p (never from inside the fine advance).
    !> Assemble the persistent global coarse cell-boundary arrays. Each rank writes the boundaries of the cells it owns (shared
    !! inter-rank faces written identically by both neighbours) into a sentinel-filled global array; an elementwise MAX allreduce
    !! recovers the exact global array on every rank. Grid fixed for the run, so this runs once.
    impure subroutine s_amr_build_global_cb()

        integer             :: j
        real(wp), parameter :: sentinel = -huge(1._wp)

        allocate (amr_gxcb(-1:m_glb)); amr_gxcb = sentinel
        do j = -1, m
            amr_gxcb(start_idx(1) + j) = x_cb(j)
        end do
        call s_mpi_allreduce_array_max(amr_gxcb, m_glb + 2)
        if (n_glb > 0) then
            allocate (amr_gycb(-1:n_glb)); amr_gycb = sentinel
            do j = -1, n
                amr_gycb(start_idx(2) + j) = y_cb(j)
            end do
            call s_mpi_allreduce_array_max(amr_gycb, n_glb + 2)
        end if
        if (p_glb > 0) then
            allocate (amr_gzcb(-1:p_glb)); amr_gzcb = sentinel
            do j = -1, p
                amr_gzcb(start_idx(3) + j) = z_cb(j)
            end do
            call s_mpi_allreduce_array_max(amr_gzcb, p_glb + 2)
        end if

    end subroutine s_amr_build_global_cb

    !> Do two coarse-index boxes [alo:ahi] and [blo:bhi] overlap? Collapsed dims (n_glb/p_glb == 0) never disqualify.
    pure logical function f_amr_boxes_overlap(alo, ahi, blo, bhi) result(ov)

        integer, intent(in) :: alo(3), ahi(3), blo(3), bhi(3)

        ov = alo(1) <= bhi(1) .and. ahi(1) >= blo(1)
        if (n_glb > 0) ov = ov .and. alo(2) <= bhi(2) .and. ahi(2) >= blo(2)
        if (p_glb > 0) ov = ov .and. alo(3) <= bhi(3) .and. ahi(3) >= blo(3)

    end function f_amr_boxes_overlap

    !> Multi-level nesting: index of the covering level-(level(k)-1) block that block k refines (its coarse parent), or 0 when block
    !! k is level 1 (parent is the L0 base grid). Regions are in L0 cell indices at every level, so the parent is the level-below
    !! block whose box contains k's; proper nesting guarantees exactly one, and the first overlap is returned.
    pure integer function f_amr_parent_block(k) result(p)

        integer, intent(in) :: k
        integer             :: j

        p = 0
        if (amr_block_level(k) <= 1) return
        do j = 1, amr_num_blocks
            if (amr_block_level(j) == amr_block_level(k) - 1 .and. f_amr_boxes_overlap(amr_region_lo_all(:,k), &
                & amr_region_hi_all(:,k), amr_region_lo_all(:,j), amr_region_hi_all(:,j))) then
                p = j
                return
            end if
        end do

    end function f_amr_parent_block

    !> Pack (dir=+1) / unpack (dir=-1) the fine cells of slot's q_cons over [dlo:dhi] in dim d, full transverse, all sys_size, in a
    !! fixed (i, d-index, transverse) order so a packer and unpacker with matching extents align cell-for-cell. GPU: only this
    !! buff_size-deep near-seam slab is moved device<->host (host<-device before a pack, device<-host after an unpack), interior
    !! transverse (0:fm) only, exactly the cells touched below, so the round-trip equals a full-field update at a tiny fraction of
    !! the volume (the halo runs per stage).
    impure subroutine s_amr_fine_slice(slot, d, dlo, dhi, buf, dir)

        integer, intent(in)                 :: slot, d, dlo, dhi, dir
        real(wp), intent(inout), contiguous :: buf(:)
        integer                             :: i, a, b, c, fm(3), na, nb, nc, loc

        fm(1) = amr_slots(slot)%m; fm(2) = amr_slots(slot)%n; fm(3) = amr_slots(slot)%p
        loc = amr_loc_of(slot)
        nc = dhi - dlo + 1
        ! Pack (dir=1) / unpack (dir=-1) the near-seam slab on the device straight into the contiguous buffer buf, then move only
        ! buf
        ! host<->device. flang miscomputes a strided section (seam dim d < num_dims) of a block's conserved field in a
        ! target-update map clause, corrupting the 2D+ np>1 seam ghosts; the base-grid halo (s_mpi_sendrecv_variables_buffers)
        ! device-packs into a contiguous buffer for the same reason. buf index runs a fastest, then b, then c, then i, so a pack and
        ! an unpack with matching extents align cell-for-cell (na/nb are the transverse fine sizes, nc the slab depth).
        #:for D, TA, TB in [(1, 2, 3), (2, 1, 3), (3, 1, 2)]
            #:set IDX = {1: 'c, a, b', 2: 'a, c, b', 3: 'a, b, c'}[D]
            if (d == ${D}$) then
                na = fm(${TA}$) + 1; nb = fm(${TB}$) + 1  ! scalars; kernel loop bounds must use na-1/nb-1, not fm(..), so no host
                !                     array is referenced in the device region (nvfortran/Cray demand it present)
                if (dir == 1) then  ! host <- device: pack on the device, copyout moves the contiguous buffer to host
                    $:GPU_PARALLEL_LOOP(collapse=4, copyout='[buf]')
                    do i = 1, sys_size
                        do c = dlo, dhi
                            do b = 0, nb - 1
                                do a = 0, na - 1
                                    buf(1 + a + na*(b + nb*(c - dlo + nc*(i - 1)))) = real(amr_cons_st(${IDX}$, i, loc), wp)
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                else  ! device <- host: copyin moves the contiguous buffer to device, then unpack on the device
                    $:GPU_PARALLEL_LOOP(collapse=4, copyin='[buf]')
                    do i = 1, sys_size
                        do c = dlo, dhi
                            do b = 0, nb - 1
                                do a = 0, na - 1
                                    amr_cons_st(${IDX}$, i, loc) = real(buf(1 + a + na*(b + nb*(c - dlo + nc*(i - 1)))), stp)
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
            end if
        #:endfor

    end subroutine s_amr_fine_slice

    !> Convert a physical-space bbox to a global coarse-index bbox padded by pad_cells.
    pure subroutine s_lag_phys_to_cells(pmin, pmax, pad_cells, blo, bhi)

        real(wp), dimension(3), intent(in) :: pmin, pmax
        integer, intent(in)                :: pad_cells
        integer, intent(out)               :: blo(3), bhi(3)

        blo(1) = int((pmin(1) - glb_bounds(1)%beg)/dx(0)) - pad_cells
        bhi(1) = int((pmax(1) - glb_bounds(1)%beg)/dx(0)) + pad_cells
        blo(2) = 0; bhi(2) = 0; blo(3) = 0; bhi(3) = 0
        if (n_glb > 0) then
            blo(2) = int((pmin(2) - glb_bounds(2)%beg)/dy(min(1, n))) - pad_cells
            bhi(2) = int((pmax(2) - glb_bounds(2)%beg)/dy(min(1, n))) + pad_cells
        end if
        if (p_glb > 0) then
            blo(3) = int((pmin(3) - glb_bounds(3)%beg)/dz(0)) - pad_cells
            bhi(3) = int((pmax(3) - glb_bounds(3)%beg)/dz(0)) + pad_cells
        end if

    end subroutine s_lag_phys_to_cells

    !> Rank-local per-stage guard: the local bubbles' padded bbox must stay clear of the current block. Catches an overlapping
    !! initial placement on the first stage and drift that outran the regrid margin afterwards.
    impure subroutine s_amr_check_lag_clear()

        real(wp), dimension(3) :: pmin_loc, pmax_loc
        integer                :: blo(3), bhi(3)
        logical                :: ovl

        if (.not. bubbles_lagrange) return
        call s_lag_cloud_bbox_local(pmin_loc, pmax_loc)
        if (pmin_loc(1) > pmax_loc(1)) return  ! no bubbles on this rank
        call s_lag_phys_to_cells(pmin_loc, pmax_loc, mapCells + 2, blo, bhi)
        ovl = blo(1) <= amr_slots(amr_cur)%region%hi(1) .and. bhi(1) >= amr_slots(amr_cur)%region%lo(1)
        if (n_glb > 0) ovl = ovl .and. blo(2) <= amr_slots(amr_cur)%region%hi(2) .and. bhi(2) >= amr_slots(amr_cur)%region%lo(2)
        if (p_glb > 0) ovl = ovl .and. blo(3) <= amr_slots(amr_cur)%region%hi(3) .and. bhi(3) >= amr_slots(amr_cur)%region%lo(3)
        if (ovl) then
            call s_mpi_abort('amr with Lagrangian bubbles: the bubble cloud (positions + smearing support) ' &
                             & // 'overlaps an active fine block, where two-way coupling would be lost. Keep the initial ' &
                             & // 'block clear of the cloud; under dynamic regrid, reduce amr_regrid_int or increase ' &
                             & // 'amr_buf so the exclusion margin covers the cloud drift between regrids')
        end if

    end subroutine s_amr_check_lag_clear

    !> minmod slope limiter: 0 if a,b differ in sign, else the smaller-magnitude argument.
    pure elemental function minmod(a, b) result(m)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: a, b
        real(wp)             :: m

        if (a*b <= 0._wp) then
            m = 0._wp
        else if (abs(a) < abs(b)) then
            m = a
        else
            m = b
        end if

    end function minmod

    !> Allocate the pooled q_prim/rhs advance scratch (idempotent). The lockstep driver argument-associates the scratch for every
    !! block, owned or not, so it must exist on every rank, including one that never allocates a slot. Called from the one point per
    !! mode where mbuf* are final: end of s_initialize_amr_module when l0_ntile == 0 (pure AMR), and after s_l0_tiles_init's mbuf
    !! union when tiles exist (pure-L0 and coexist; the union can enlarge mbuf* past the fine-only values, so an earlier allocation
    !! would undersize the scratch). rhs mirrors the per-slot igr widening.
    impure subroutine s_amr_scr_init()

        integer                         :: i, slo(3), shi(3)
        type(scalar_field), allocatable :: tmp_p(:), tmp_r(:)

        if (allocated(amr_scr_prim)) return
        ! the batched advance stacks amr_br_batch blocks along the last active dimension (see amr_cons_br)
        slo = [mbuf1_lo, mbuf2_lo, mbuf3_lo]; shi = [mbuf1_hi, mbuf2_hi, mbuf3_hi]
        shi(num_dims) = slo(num_dims) + amr_br_batch*(shi(num_dims) - slo(num_dims) + 1) - 1
        ! CCE OpenMP-offload leaves a bare module-scope derived-type allocatable's descriptor uninitialized, so a direct
        ! allocate here aborts with `lib-4425 INTERNAL ERROR-Unitialized descriptor for ALLOCATE statement argument` on the
        ! host allocate, not the device map. Same defect, same workaround as amr_cg above: allocate a local, which gets a
        ! valid descriptor, then hand it over with move_alloc and map afterwards. GPU_DECLARE alone does not avoid it.
        allocate (tmp_p(1:sys_size)); call move_alloc(tmp_p, amr_scr_prim)
        allocate (tmp_r(1:sys_size)); call move_alloc(tmp_r, amr_scr_rhs)
        $:GPU_ENTER_DATA(create='[amr_scr_prim, amr_scr_rhs]')
        do i = 1, sys_size
            @:ALLOCATE(amr_scr_prim(i)%sf(mbuf1_lo:shi(1), mbuf2_lo:shi(2), mbuf3_lo:shi(3)))
            if (igr) then
                @:ALLOCATE(amr_scr_rhs(i)%sf(mbuf1_lo:shi(1), min(mbuf2_lo, -1):max(shi(2), 1), min(mbuf3_lo, -1):max(shi(3), 1)))
            else
                @:ALLOCATE(amr_scr_rhs(i)%sf(mbuf1_lo:shi(1), mbuf2_lo:shi(2), mbuf3_lo:shi(3)))
            end if
            @:ACC_SETUP_SFs(amr_scr_prim(i))
            @:ACC_SETUP_SFs(amr_scr_rhs(i))
        end do
        if (ib .and. amr_batched_advance) then
            allocate (tmp_p(1:sys_size)); call move_alloc(tmp_p, amr_scr_prim_blk)
            $:GPU_ENTER_DATA(create='[amr_scr_prim_blk]')
            do i = 1, sys_size
                @:ALLOCATE(amr_scr_prim_blk(i)%sf(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi))
                @:ACC_SETUP_SFs(amr_scr_prim_blk(i))
            end do
        end if

    end subroutine s_amr_scr_init

    !> Slot index of regrid-managed fine block k in the shared pool: tiles occupy [1..l0_slot_off], fine blocks [l0_slot_off+1..].
    !! Identity (l0_slot_off=0) until L0 tiles + AMR coexist.
    pure integer function f_l0_slot(k) result(s)

        integer, intent(in) :: k

        s = l0_slot_off + k

    end function f_l0_slot

end module m_amr_state
