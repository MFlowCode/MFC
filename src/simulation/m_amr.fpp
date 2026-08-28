!>
!!@file
!!@brief Contains module m_amr

#:include 'macros.fpp'

!> @brief Block-structured AMR: up to amr_max_blocks refined blocks (2:1 or 4:1 per amr_ref_ratio), optionally nested to
!! amr_max_level, advanced with the shared solver via grid-state swap and conservatively coupled to each block's parent level (ghost
!! prolongation, Berger-Colella flux reflux, restriction); optional dt/2 subcycling and dynamic regrid.
module m_amr

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
    use m_amr_xchg_audit  ! I1a: per-call-site accounting of every AMR p2p transfer (s_xa_rec + XA_* site ids)
    use m_ibm, only: s_ibm_alloc_fine, s_ibm_setup_fine, s_ibm_swap_to_fine, s_ibm_restore_from_fine, s_ibm_correct_state, &
        & s_update_mib, moving_immersed_boundary_flag, num_gps, ib_markers
    use m_hypoelastic, only: s_hypoelastic_update_fd_coeffs
    use m_weno, only: s_compute_weno_coefficients
    use m_active_box, only: ab_active
    use m_bubbles_EL, only: s_lag_cloud_bbox_local
    use m_igr, only: jac, jac_old

    implicit none

    private
    public :: t_level, amr_maxc, amr_maxc_fit, s_initialize_amr_module, s_populate_amr_fine, s_interpolate_coarse_to_fine, &
        & s_restrict_fine_to_coarse, s_finalize_amr_module, s_amr_stage_fill_wave, s_amr_parent_fill_wave, &
        & s_amr_fine_stage_advance, s_amr_fine_fine_halo, s_amr_advance_fine_subcycle_all, s_set_amr_fine_geometry, &
        & s_amr_relax_fine, s_amr_setup_ib, s_amr_p2p_reflux_faces, s_amr_reflux_faces_wave, s_amr_freg_wave, &
        & s_amr_restrict_wave, s_amr_convert_prim_batch, amr_prim_batch, s_amr_reflux_to_parent, s_l0_tiles_init, &
        & s_l0_advance_stage, s_l0_advance_stage_rhs, s_l0_advance_stage_rk, s_l0_copy_coarse_to_tiles, &
        & s_l0_scatter_tiles_to_coarse, s_l0_add_reflux_to_tiles, s_l0_restrict_to_tiles, s_l0_tiles_finalize, s_l0_forced_remap, &
        & s_l0_rebalance, s_l0_fill_tiles_from_coarse
    ! s_amr_swap_to_fine / s_amr_restore_coarse / s_amr_fill_fine_ghosts / amr_dt_fine are internal (no external caller); keeping
    ! them private makes "the swap has exactly these audited call sites" a compiler guarantee, not a convention.
    !> Block/slot state and fine-distribution services consumed by m_amr_regrid and m_amr_restart (the drivers split out of this
    !! module). State stays HERE - only the drivers moved.
    public :: s_amr_br_load_all, s_amr_br_store_all, amr_br_w, amr_br_batch, amr_rg_gather
    public :: s_amr_build_gather_plan, amr_gpl_valid
    public :: s_amr_gather_chunk_post, s_amr_gather_chunk_send, s_amr_gather_consume_box, amr_gath_chunk, s_amr_cov_note, &
        & amr_cad_tot, amr_cad_esc, amr_cad_armed
    public :: amr_slots, amr_cons_st, amr_stor_st, amr_loc_of, amr_seam_pairs_dirty, amr_mesh_epoch, amr_tag_base, &
        & amr_xchg_coarse_ghosts, amr_cpat_mar, s_amr_alloc_slot, s_amr_alloc_slot_stash, s_amr_prereserve_stash, &
        & s_amr_free_slot, s_amr_reconcile_slots, s_amr_reduce_xchg_flag, s_amr_assign_block_owners, s_amr_gather_coarse_patch, &
        & s_amr_gather_send_flush, s_amr_gather_coarse_patch_pbmv, s_amr_prolong_pbmv, s_amr_exchange_coarse_cons_halo, &
        & s_lag_phys_to_cells, s_amr_body_bbox, s_amr_expand_box_over_bodies, s_amr_tile_box, f_amr_seam_dim, &
        & f_amr_boxes_overlap, f_l0_slot

    !> Fine-level time step for subcycling (= 0.5*dt after init; 0 when amr is off).
    real(wp) :: amr_dt_fine = 0._wp

    !> Realizability floor for prolonged Euler-Euler bubble POSITIVE moments (radius nR, non-polytropic partial pressure npb / vapor
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
        !> conserved state lives in the FLAT STORE amr_cons_st, indexed by this slot's dense local index amr_loc_of - not here.
        !> SSP-RK stage storage lives in the FLAT STORE amr_stor_st, indexed by this slot's amr_loc_of.
        type(scalar_field), allocatable :: q_prim(:)  !< primitive stage (fine advance, same bounds)
        type(scalar_field), allocatable :: rhs(:)     !< RHS (fine interior only: 0:m, 0:n, 0:p)
        !> non-polytropic QBMM quadrature side-state on the block (nnode x nb per cell). pb/mv evolve cell-locally (their rhs reads
        !! only the local cell + the block's own moment fluxes), so the fine treatment is prolong -> advance -> restrict with no
        !! reflux; ghosts feed the widened-idwint conversions and are prolonged piecewise-constant (CHyQMOM realizability, like the
        !! moments).
        type(pres_field) :: pb_f, mv_f        !< fine pb/mv (ghost-inclusive)
        type(pres_field) :: pb_stor, mv_stor  !< SSP-RK step-entry backup (also the regrid bounce)
        !> subcycle ghost-lerp sources at coarse t^n / t^{n+1} (ghost shell only): ghost pb feeds the mixture pressure in the
        !! widened conversion, so it needs the same time fidelity as q_cons
        type(pres_field) :: pb_ghost_a, mv_ghost_a
        type(pres_field) :: pb_ghost_b, mv_ghost_b
    end type t_level

    !> Fixed pool of refined-block slots (at init one slot is active; dynamic regrid activates up to amr_max_blocks). The working
    !! slot amr_cur (m_global_parameters) selects which slot every per-block routine operates on.
    type(t_level), allocatable :: amr_slots(:)

    !> DENSE LOCAL INDEX for live slots. `amr_slots` is indexed by GLOBAL block index and allocated lazily, so live slots are SPARSE
    !! across 1:amr_max_blocks (1024 by default). A contiguous per-block field store - the layout AMReX's MultiFab uses and the
    !! prerequisite for batching one kernel over all blocks - must be indexed DENSELY by a local index instead, or it would have to
    !! be sized for the whole global pool. `amr_loc_of(g)` is the local index of global slot g (0 if not live), `amr_loc_n` is the
    !! high-water mark, and freed indices are recycled through `amr_loc_free` so the dense range stays tight under regrid churn.
    !! Pure bookkeeping: nothing reads it yet.
    integer, allocatable :: amr_loc_of(:)    !< global slot -> dense local index, 0 if not live
    integer, allocatable :: amr_loc_free(:)  !< stack of recycled local indices
    integer              :: amr_loc_n = 0    !< high-water mark of local indices handed out
    !> Last high-water REPORTED by the store trip-wire. Module scope (not `save`) so the wire prints one line per NEW high-water
    !! instead of one per slot alloc (~224/regrid/rank would both flood stderr and perturb the run).
    integer :: amr_st_hw = 0
    !> TRACK S instrumentation: bytes each rank ALLOCATES for, and RECEIVES from, global collectives during one regrid. These are
    !! the quantities that must stay O(1) in problem size; wall time at a single problem size cannot see them (rg:clus measured
    !! 0.6%% of wall), so they are counted explicitly. amr_gb_tag counted the level-1 tag ALLGATHERV and is now permanently 0: S3.1
    !! deleted that collective, so a nonzero value means the gather came back. amr_gb_win = the gwin-pair ALLGATHERVs (levels >= 2,
    !! still live -- S3.3), amr_gb_cost = the per-box cost ALLREDUCE.
    integer(8) :: amr_gb_tag = 0, amr_gb_win = 0, amr_gb_cost = 0
    !> S3.0a/S3.0b instrumentation: shape and reduction cost of the Berger-Rigoutsos clustering tree (`s_amr_cluster`). Named
    !! amr_cl_ (clustering), NOT amr_br_ -- that prefix already means the batched BRIDGE in this module. Tree DEPTH decides whether
    !! S3 can fuse one collective per tree LEVEL; BR splits at signature holes, not midpoints, so the tree is not balanced by
    !! construction and depth must be measured. TWO independent maxima are kept: one running max with a tie-break would pair a deep
    !! tiny subtree with the wrong leaf count and could fake either verdict. amr_cl_maxdep / amr_cl_maxdep_leaf = the DEEPEST call
    !! and its leaf count (catches a one-box-at-a-time peel chain). amr_cl_lmax / amr_cl_ldepth = the LARGEST call and its depth
    !! (the tree that actually dominates the cost). amr_cl_nodes = nodes visited = the collectives a per-node distributed recursion
    !! pays. amr_cl_rb = bytes that recursion would ALLREDUCE: one fused reduction per node carrying the 1D signature of every
    !! splittable axis. Sized on the UNTRIMMED box, because trim only shrinks to the contained tags own bbox, so the trimmed
    !! signature is a slice of the untrimmed one and that single reduction serves trim, count AND split.
    integer    :: amr_cl_maxdep = 0, amr_cl_maxdep_leaf = 0, amr_cl_lmax = 0, amr_cl_ldepth = 0
    integer(8) :: amr_cl_nodes = 0, amr_cl_rb = 0, amr_cl_rb_now = 0
    !> S3.2a: does a tree node's box fit inside ONE rank's subdomain? If so the owning rank holds every tag in the subtree and the
    !! whole remaining subtree needs NO communication -- that is the property the S3.2 design rests on, and this measures how much
    !! of the tree has it. shr = nodes spanning >1 rank (must be exchanged), loc = rank-local (free). amr_cl_shr_maxdep is the depth
    !! of the deepest node still shared, i.e. the number of levels S3.2 must communicate. SPLIT BY PATH. The level-2 forest clusters
    !! inside a parent window, so its boxes are small and mostly rank-local BY CONSTRUCTION; folding it in with the level-1 tree
    !! would inflate the local fraction and make S3.2's premise look true even if level-1 -- the only path that reduces today -- is
    !! mostly shared. `_r` = the reducing (level-1) path.
    integer(8) :: amr_cl_shr_nodes = 0, amr_cl_shr_rb = 0, amr_cl_loc_nodes = 0, amr_cl_loc_rb = 0
    integer(8) :: amr_cl_shr_nodes_r = 0, amr_cl_shr_rb_r = 0, amr_cl_loc_nodes_r = 0, amr_cl_loc_rb_r = 0
    integer    :: amr_cl_shr_maxdep = 0, amr_cl_shr_maxdep_r = 0
    !> S3.2a-2: what the shallow phase would cost THIS rank. shr_rb_r counts every shared node and so prices the ALLREDUCE form,
    !! where each node's whole signature lands on every rank. Under S3.2's sparse per-depth exchange a rank touches only the shared
    !! nodes its own subdomain OVERLAPS, so these count that subset -- the quantity that has to be sublinear in P for W4.
    integer(8) :: amr_cl_me_nodes_r = 0, amr_cl_me_rb_r = 0
    !> TRACK T (T0b gate): regrid migration volume. An old block is ISENT to EVERY new-owner rank whose box overlaps it, so the cost
    !! is fan-out x block bytes, not one send per block. amr_mig_blk counts blocks that had to move at all, amr_mig_snd counts the
    !! sends, amr_gb_mig the bytes. fan-out = snd/blk is the reducible quantity: if it is ~1 the volume is inherent and hysteresis
    !! buys nothing.
    integer(8) :: amr_gb_mig = 0, amr_mig_snd = 0, amr_mig_blk = 0
    public :: amr_gb_tag, amr_gb_win, amr_gb_cost
    public :: amr_cl_maxdep, amr_cl_maxdep_leaf, amr_cl_lmax, amr_cl_ldepth, amr_cl_nodes, amr_cl_rb, amr_cl_rb_now
    public :: amr_cl_shr_nodes, amr_cl_shr_rb, amr_cl_loc_nodes, amr_cl_loc_rb, amr_cl_shr_maxdep
    public :: amr_cl_shr_nodes_r, amr_cl_shr_rb_r, amr_cl_loc_nodes_r, amr_cl_loc_rb_r, amr_cl_shr_maxdep_r
    public :: amr_cl_me_nodes_r, amr_cl_me_rb_r
    public :: s_amr_ranks_overlapping  !< exported for the S3.2a scope measurement in m_amr_regrid
    public :: amr_gb_mig, amr_mig_snd, amr_mig_blk
    integer :: amr_loc_nfree = 0  !< depth of the recycle stack

    !> FLAT PER-BLOCK FIELD STORE, indexed (x, y, z, var, LOCAL slot) by the dense index above. One contiguous module array
    !! replacing a per-slot vector of independently allocated scalar_fields: the layout AMReX's MultiFab uses, and the prerequisite
    !! for running ONE kernel over every live block instead of one kernel per block. Every slot's arrays carry the same mbuf
    !! extents, so a single array serves them all. The ghost pair exists only under amr_subcycle. Sized by s_amr_st_reserve.
    real(stp), allocatable, dimension(:,:,:,:,:) :: amr_cons_st, amr_stor_st, amr_gst_a, amr_gst_b
    $:GPU_DECLARE(create='[amr_cons_st, amr_stor_st, amr_gst_a, amr_gst_b]')
    !> local slots the store is sized for; grows and never shrinks, but plateaus at the rebuild-transient high-water because
    !! s_amr_compact_store re-densifies the index space every reconcile
    integer :: amr_st_cap = 0

    !> 2a PRIM LANDING ZONE for the batched cons->prim conversion (s_amr_convert_prim_batch): the COMPUTED prim vars only - the
    !! contiguous eqn_idx%mom%beg..eqn_idx%E range (velocities + pressure), var dim 1..num_vels+1. The aliased prim vars (cont, adv,
    !! c, psi) ride the cons copy-in inside s_compute_rhs as always. Per-STAGE scratch: rewritten by every batch call, so store
    !! growth discards it (no staging round trip). Sized with the store in s_amr_st_reserve, only under the amr_prim_batch gate.
    real(stp), allocatable, dimension(:,:,:,:,:) :: amr_prim_st
    $:GPU_DECLARE(create='[amr_prim_st]')
    !> per-dense-slot batch metadata: participating fine slot + its idwbuff window (host-filled each batch call)
    integer, allocatable :: amr_bt_lo(:,:), amr_bt_hi(:,:)  !< (3, loc)
    logical, allocatable :: amr_bt_on(:)
    $:GPU_DECLARE(create='[amr_bt_lo, amr_bt_hi, amr_bt_on]')
    !> 2a master gate, derived once at init: the batched conversion covers the plain multi-fluid configs (5/6-eq, WENO, with/without
    !! viscous); every feature that adds conversion write-set members or changes its inputs (igr, chemistry, relativity,
    !! hypoelasticity, mhd, cont_damage, ib, Lagrangian bubbles, subcycle) falls back to the per-block conversion path unchanged.
    logical :: amr_prim_batch = .false.

    !> COPY BRIDGE to the shared solver. s_compute_rhs, s_ibm_correct_state, s_pressure_relaxation_procedure and
    !! s_infinite_relaxation_k all take type(scalar_field), dimension(sys_size) and serve the monolithic path too, so the flat store
    !! cannot be handed to them; a pointer view into the store is not attachable on the OpenMP-offload backend (measured - see
    !! docs/documentation/amr_block_batching.md). One block-shaped scalar_field array bridges instead: load it from the store, call,
    !! store it back. All four dummies are intent(inout) - s_compute_rhs writes the buffer region through
    !! s_populate_variables_buffers - so BOTH directions are required at every crossing.
    type(scalar_field), allocatable :: amr_cons_br(:)
    !> Batched-bridge geometry. amr_br_w is one block's buffered k-width; block loc of a batch sits at k-offset (loc-1)*amr_br_w,
    !! carrying its own ghost shell, so consecutive blocks are separated by TWO ghost shells and no block's stencil can reach
    !! another's interior - the property the batched advance's correctness rests on.
    integer :: amr_br_w = 0, amr_br_nblk = 0
    !> P1 pooled advance scratch: the fused per-block fine advance (rhs then rk on ONE block, s_amr_fine_stage_advance) leaves no
    !! cross-block q_prim/rhs lifetime, so every fine block shares this one slot-shaped pair instead of carrying per-slot arrays
    !! (~2x105 MiB per live slot at the S0 point - the np>=8 live-footprint blocker AND the alloc/free churn that fed the
    !! libomptarget retention plateau). Same shared-scratch pattern as amr_rhs_pb_f/amr_cg. L0 tile slots are the exception and keep
    !! per-slot arrays (see s_amr_alloc_slot).
    type(scalar_field), allocatable :: amr_scr_prim(:), amr_scr_rhs(:)
    !> True only while the REGRID path is inside s_amr_gather_coarse_patch, so the WAITALL bracket attributes to rb:wait rather than
    !! mixing in the per-step gather that shares this routine.
    logical :: amr_rg_gather = .false.
    !> Blocks per batched call. BOUNDED on purpose. Sizing the bridge per-block (one slot for every block) OOMed the device on the
    !! 400^3 case: a buffered block is ~110 MB at cap 64 and the slot cap grows geometrically and never shrinks, so the bridge alone
    !! reached multiple GB against a working set already near 43 GiB/GCD. A fixed window keeps bridge memory O(1) in the block count
    !! - which the exascale goal requires anyway - and still collapses launches by this factor.
    !> 1 = the batched path is DORMANT. G0.2 measured PH_RHS at 54-57%% GPU-busy, so batching the fine advance is a ~1.09x item
    !! (plan rule: >50%% -> Track B waits for Track R); a batch >1 would allocate 8x the bridge for nothing.
    integer, parameter :: amr_br_batch = 1
    integer            :: amr_maxc(3)  !< max coarse block cells per dim: (m_glb+1)/2 etc.; 1 for collapsed dims

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
    !! amr_mesh_epoch now lives in m_global_parameters (m_amr_registers keys its participation-map rebuild on it and
    !! cannot use m_amr); it is use-associated here and re-exported, so importers of m_amr are unchanged.
    !! per-family plan message tag bases (families F1..F7, amr_plan_based_exchange.md): amr_max_blocks + 100*f keeps
    !! the plan tag space disjoint from the legacy per-box space (tags in [1..amr_max_blocks]) while families convert;
    !! the epoch is folded in as base + mod(amr_mesh_epoch, 100). The init MPI_TAG_UB assert is the scale tripwire: it caps
    !! GLOBAL blocks near 2**21/1 ~ 2.1e6, i.e. about 28k ranks at ~75 boxes/rank, which is the SECOND scaling wall after W4.
    !! The amr_max_blocks term can only go once NO family uses per-box tags. Verified 2026-08-27: 19 of 41 AMR p2p call sites
    !! still tag per box (F1's unconverted path uses amr_cur, migration uses the column index), and the subcycle sites are an
    !! EXPLICIT deferral to increment I8, not I7 -- I7's own boundary is that any family left per-box keeps its tables.
    integer :: amr_tag_base(7) = 0
    !! cached per-block P2P overlap-rank lists (rebuilt with the seam list - same dirty flag): amr_ovl_gather(:,k) = ranks whose
    !! owned coarse range (s_amr_rank_coarse_range) intersects block k's amr_cpat_mar-padded patch box (gather contributors);
    !! amr_ovl_scatter(:,k) = ranks whose coarse interior (s_amr_rank_interior) intersects block k's region box (restrict-scatter
    !! destinations). Built by O(overlap) inversion (s_amr_ranks_overlapping), rank-ASCENDING and NOT owner-excluded (consumers keep
    !! their owner skip), so iterating a list reproduces the replaced per-call 0..num_procs-1 scan's MPI send/recv order exactly.
    !> (max-overlap, amr_max_blocks); sized in s_amr_build_seam_pairs
    integer, allocatable :: amr_ovl_gather(:,:), amr_ovl_scatter(:,:)
    integer, allocatable :: amr_ovl_gather_n(:), amr_ovl_scatter_n(:)  !< per-block list lengths
    !> Rebuild gather PLAN (gather-batching step 1): the whole rebuild's gather message set, derived up front by
    !! s_amr_build_gather_plan from the replicated caches. Per level-1 slot: contributor count/ranks/message sizes (owner excluded,
    !! list order = amr_ovl_gather order). Per level>=2 slot: the parent-owner source rank (-1 when co-located, no message) and its
    !! message size. The per-box gather ASSERTS its inline derivation against this plan (guarded on amr_rg_gather + amr_gpl_valid,
    !! so per-step gathers never consult it); the chunked exchange of step 2 may trust the plan only because those asserts prove it
    !! reproduces today's message set exactly.
    integer, allocatable :: amr_gpl_nsrc(:), amr_gpl_src(:,:), amr_gpl_sz(:,:), amr_gpl_psrc(:), amr_gpl_psz(:)
    logical              :: amr_gpl_valid = .false.  !< true only between plan build and the end of the rebuild box loop
    !> Chunked rebuild gather (step 2, amr_regrid_gather_batching.md): the rebuild box loop runs in chunks of amr_gath_chunk boxes -
    !! every owned box's recvs (level-1 contributor slices AND split level>=2 parent patches) are pre-posted from the plan into one
    !! flat pool, this rank's sends are issued (level>=2 only when the parent was consumed in an EARLIER chunk - a same-chunk
    !! parent's store is unbuilt until its own consume, so that send stays at the child's consume position), then boxes are consumed
    !! in order with a per-box wait. Requests are appended in box order, so each box's recvs are the contiguous run amr_gcr_r0 :
    !! +amr_gcr_nr-1. The pool/request arrays grow monotonically and are reusable across chunks ONLY because the consume phase waits
    !! every owned box's requests unconditionally inside its own chunk.
    integer, parameter    :: amr_gath_chunk = 32  !< boxes per chunk: staging memory vs message batching
    real(wp), allocatable :: amr_gcr_pool(:)  !< flat recv staging for one chunk
    integer, allocatable  :: amr_gcr_req(:), amr_gcr_off(:)  !< request handle + pool offset per posted recv
    integer               :: amr_gcr_r0(amr_gath_chunk), amr_gcr_nr(amr_gath_chunk)  !< per chunk-local box: first recv, count
    logical               :: amr_gcr_sent(amr_gath_chunk)  !< chunk-local: level>=2 send already issued in the send phase
    integer               :: amr_gcr_n = 0  !< posted recvs in the current chunk
    !> I2a stage-fill WAVE (plan-based exchange, amr_plan_based_exchange.md I2): the non-subcycle level-1 per-stage fill's F1 q_cons
    !! + F3 pb/mv gathers as ONE per-(peer, family) aggregated exchange per RK stage, replacing the per-box owner-WAITALL /
    !! contributor-flush rendezvous chain. Transfer records are SoA flat arrays (no derived types); the wire layout of each peer
    !! message is the ascending-box concatenation of [XA_NH header | slab], which sender and receiver derive independently from the
    !! replicated caches (rank coarse ranges x patch boxes), so no metadata is exchanged. Plans are rebuilt every wave (caching on
    !! amr_mesh_epoch is increment I6). All scratch is high-water and its contents never survive a wave; the rank-indexed build
    !! counters (amr_fw_map/nx/pq/pp) are re-zeroed for touched ranks after each build so they stay all-zero between builds.
    integer               :: amr_fw_snx = 0, amr_fw_rnx = 0, amr_fw_snp = 0, amr_fw_rnp = 0
    integer, allocatable  :: amr_fw_sblk(:), amr_fw_sbl(:,:), amr_fw_sbh(:,:), amr_fw_spi(:), amr_fw_sqo(:), amr_fw_spo(:)
    integer, allocatable  :: amr_fw_rblk(:), amr_fw_rbl(:,:), amr_fw_rbh(:,:), amr_fw_rpi(:), amr_fw_rqo(:), amr_fw_rpo(:)
    integer, allocatable  :: amr_fw_sprank(:), amr_fw_sqsz(:), amr_fw_spsz(:), amr_fw_snxp(:), amr_fw_sqbase(:), amr_fw_spbase(:)
    integer, allocatable  :: amr_fw_rprank(:), amr_fw_rqsz(:), amr_fw_rpsz(:), amr_fw_rnxp(:), amr_fw_rqbase(:), amr_fw_rpbase(:)
    integer, allocatable  :: amr_fw_map(:), amr_fw_nx(:), amr_fw_pq(:), amr_fw_pp(:)  !< rank-indexed build scratch (0:num_procs-1)
    real(wp), allocatable :: amr_fw_sq(:), amr_fw_sp(:), amr_fw_rq(:), amr_fw_rp(:)  !< wire pools (live across the ISENDs)
    integer, allocatable  :: amr_fw_req(:), amr_fw_reqw(:)  !< requests + expected recv word counts (-1 for sends; debug check)
    !> [amr-cov] dead-byte accounting (expert-audit increment, amr_action_plan.md 2026-08-22): partition each gather family's patch
    !! words into live vs provably-dead. (1) step-fill: the ghost-fill kernel reads only floor(f/rr) +- 1, so the patch's interior
    !! core (region shrunk by 1 per face) is never read. (2) rebuild level-1: the same-level carry-forward overwrites prolonged
    !! cells covered by old blocks (shrunk by 1 for the minmod stencil - a conservative under-count). (3) rebuild level>=2: no
    !! carry-forward exists, so the patch is live by construction. Deterministic - identical across reruns; decides ring/coverage
    !! clipping vs T1 (pre-registered: dead > 50% on either family promotes clipping).
    !> 1=step-fill, 2=rebuild L1, 3=rebuild L>=2 [amr-cad] regrid-cadence containment audit: level-1 tags counted at each regrid,
    !! and how many fell OUTSIDE the pre-regrid level-1 coverage (a feature that evolved unrefined since the last regrid - the tag
    !! buffer amr_buf did not cover its drift). Zero escaped validates the (amr_regrid_int, amr_buf) pair for that run; the case
    !! validator only warns (the CFL <= 1 worst case is too strict for low-CFL cases). Incremented by m_amr_regrid, reported by
    !! s_amr_cov_report.
    integer(8) :: amr_cov_tot(3) = 0, amr_cov_dead(3) = 0
    integer(8) :: amr_cad_tot = 0, amr_cad_esc = 0
    logical    :: amr_cad_armed = .false.  !< first regrid (hierarchy population) is skipped - see s_amr_cad_count
    !> (0:num_procs-1) SFC Morton-key upper bound per rank from the cost-weighted split; owner = cut-search (f_amr_owner). The O(P)
    !! computed replacement for the O(global_blocks) amr_block_owner table (validated against it during bring-up).
    integer(kind=8), allocatable :: amr_owner_cut(:)
    !> (0:num_procs-1, 1:amr_max_level) companion cuts for FINE-block (level>=1) owners: ONE INDEPENDENT CUT PER LEVEL. Each level's
    !! boxes are balanced across all ranks on their own weight, so a deep refinement tower no longer pins its whole subtree to one
    !! rank. Per level rather than one mixed cut because same-level boxes are disjoint and so have DISTINCT Morton keys, which the
    !! cut-point binary search requires; a mixed cut would let a child share its parent's region_lo and make the search ambiguous.
    !! In no-tile AMR level 1's cut also mirrors amr_owner_cut, but in coexist amr_owner_cut is overwritten by the TILE cut, so the
    !! fine cuts are kept here for f_amr_owner (fine blocks straddle tiles and cannot be derived from the tile cut).
    integer(kind=8), allocatable :: amr_fine_cut(:,:)

    !> Regrid box size cap per dim (fixed for the run, identical on all ranks; 1 in collapsed dims): a box of at most min-over-ranks
    !! of (local extent + 1)/2 cells intersects EVERY rank in at most (its extent + 1)/2 cells, so the per-rank scratch constraint
    !! 2*(isect cells) - 1 <= local extent holds by construction. Equals amr_maxc at np=1.
    integer :: amr_maxc_fit(3) = 1

    !> SWAP CONTRACT (s_amr_swap_to_fine / s_amr_restore_coarse). The fine advance runs the shared solver on a fine block by
    !! swapping these coarse-grid globals to the block's values and restoring after: m/n/p, idwint/idwbuff, the nine coordinate
    !! arrays (sw_x_cb..sw_dz below), acoustic_source, ab_active; WENO/hypoelastic/IGR spacing coefficients are recomputed for the
    !! active grid rather than saved. RULE for anyone adding grid-dependent state: any module-level variable DERIVED from
    !! m/n/p/idwint/idwbuff/coords that a kernel reads on the fine grid must be swapped here OR refreshed on every fine call at its
    !! use site - and if GPU_DECLARE'd, its DEVICE copy too. A stale device copy of coarse bounds reads out of range on the fine
    !! grid under CCE OpenACC (the ab_int regression, fixed by a per-call GPU_UPDATE in s_compute_rhs; see m_rhs.fpp and
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
    !! each restriction, read IDENTICALLY by the device kernel and host scatter path so np=1 == np>=2 stays element-exact. Allocated
    !! only for cyl_coord. Its DEVICE copy is refreshed (GPU_UPDATE) only in s_restrict_fine_to_coarse; s_amr_restrict_to_parent
    !! reads it WITHOUT refreshing, safe ONLY because m_checker.fpp forbids cyl_coord with amr_max_level > 1 - lifting that gate
    !! makes this an ab_int-class stale-device bug (the parent fold needs a fresh per-block radius table). See the SWAP CONTRACT
    !! note above.
    real(wp), allocatable :: amr_rvw(:)
    $:GPU_DECLARE(create='[amr_rvw]')

    !> Non-polytropic QBMM fine rhs scratch, shared across slots (slots advance sequentially). Module-level raw arrays mirror the
    !! coarse rhs_pb/rhs_mv pattern: derived-type component actuals here tripped nvfortran's component-section data clauses on
    !! device.
    real(wp), allocatable :: amr_rhs_pb_f(:,:,:,:,:), amr_rhs_mv_f(:,:,:,:,:)
    $:GPU_DECLARE(create='[amr_rhs_pb_f, amr_rhs_mv_f]')

    !> Swap nesting depth. Only the OUTERMOST swap saves the coarse state into the sw_* bounce buffers, and only its matching
    !! restore puts it back; an inner swap re-installs and an inner restore is a no-op. Every nested swap site swaps to the same
    !! slot (amr_cur), so the enclosing frame's view is unchanged either way. Replaces the old paired-swap logical, whose assert
    !! forbade nesting outright - the nested sites in the RK pass are what block hoisting the restore across blocks (see
    !! docs/documentation/amr_block_batching.md).
    integer :: amr_swap_depth = 0
    !> True when the coarse grid is nonuniform (stretched grids, or 2D-axisymmetric's half-width axis cell): the spacing-dependent
    !! WENO coefficients are then recomputed for the ACTIVE grid on every block swap (the fine block's grid is itself nonuniform
    !! under stretching) and restored after. False on fully uniform grids - the recompute is skipped, behavior bit-identical.
    logical :: amr_weno_coef_recompute = .false.
    logical :: amr_grid_stretched = .false.  !< stretched coarse spacing (beyond the axisym axis half-cell; set at init)
    !> Persistent GLOBAL coarse cell-boundary arrays (indices -1:X_glb), assembled once at init. The fine-distribution owner
    !! reconstructs whole-block fine coordinates from these (its fine cells cover coarse cells it does not own the coordinate slice
    !! for). Exact on any grid.
    real(wp), allocatable :: amr_gxcb(:), amr_gycb(:), amr_gzcb(:)
    real(wp), allocatable :: sw_x_cb(:), sw_x_cc(:), sw_dx(:)
    real(wp), allocatable :: sw_y_cb(:), sw_y_cc(:), sw_dy(:)
    real(wp), allocatable :: sw_z_cb(:), sw_z_cc(:), sw_dz(:)

    !> Conservation-defect baselines (level-0 interior integrals at init; per-fluid masses + energy)

    !> True (identically on all ranks) iff some rank's fine ghost-fill stencil reads its coarse GHOST cells - the solver populates
    !! only PRIM ghosts, so the CONS ghosts the fill prolongs from must be halo-exchanged first. Never true at np=1 (block faces sit
    !! >= buff_size inside the domain).
    logical :: amr_xchg_coarse_ghosts = .false.
    !> local (un-reduced) accumulator behind amr_xchg_coarse_ghosts. s_set_amr_fine_geometry ORs each block's answer in here and
    !! s_amr_reduce_xchg_flag performs ONE allreduce for the whole scan. Previously the reduction sat inside the routine, so a
    !! regrid over nboxes blocks issued nboxes global collectives - the dominant term in the assignment's cost at scale, and one
    !! measured at 7.4-13 ms per call because it absorbs the spread in the owner-only work that precedes it.
    integer :: amr_xchg_bad = 0

    !> Per-block gathered coarse patch (fine-level distribution). The block owner may not hold the coarse cells its block refines,
    !! so before each prolongation/ghost-fill the coarse patch spanning region_lo-amr_cpat_mar : region_hi+amr_cpat_mar (the full
    !! coarse-cell reach of every prolongation stencil) is gathered here POINT-TO-POINT from the coarse-owners
    !! (s_amr_gather_coarse_patch). Stored in amr_cg as stp scalar_fields (a drop-in for the coarse q_cons in the prolong/ghost-fill
    !! kernels) in a block-LOCAL frame: amr_cg cell 0 is GLOBAL coarse cell amr_cpat_off(d). Messages carry wp, cast to stp
    !! (identity for stp coarse), so at np=1 (owner copies its own coarse) the patch equals the local coarse read bit-for-bit. Sized
    !! to the largest block.
    type(scalar_field), allocatable :: amr_cg(:)
    integer                         :: amr_cpat_mar = 0    !< coarse-cell stencil reach = (buff_size+1)/2 + 1 (matches nmar)
    integer                         :: amr_cpat_hi(3) = 0  !< amr_cg upper local bounds per dim (0 in collapsed dims)

    !> Deferred-send pool for the per-box coarse-patch gather.
    !!
    !! The gather is called once per BOX (794 of them at 400^3) and the non-owner side used a BLOCKING MPI_SEND, so every rank had
    !! to rendezvous with the owner 794 times per rebuild, in lockstep. That serialisation was measured at 45% of regrid and ~25%
    !! of total runtime across this routine's two call sites. Sends are now non-blocking and completed in batches, so a rank that
    !! only contributes data can run ahead instead of blocking on each box.
    !!
    !! The pool owns the buffers because MPI_ISEND requires them to stay live until completion - the old code deallocated sbuf
    !! immediately after the blocking send, which is exactly what must NOT happen here.
    integer, parameter    :: amr_gsnd_max = 64    !< pending sends before a forced drain (bounds pool memory)
    real(wp), allocatable :: amr_gsnd_pool(:,:)
    integer, allocatable  :: amr_gsnd_req(:)
    integer               :: amr_gsnd_n = 0
    integer               :: amr_cpat_off(3) = 0  !< GLOBAL coarse index of amr_cg local cell 0 (region_lo - amr_cpat_mar)
    !> Gathered coarse pb/mv patch for non-polytropic QBMM (analogue of amr_cg): the block's coarse-side pb/mv side-state,
    !! P2P-gathered from the coarse-cell owners into the block owner in the amr_cg patch-local frame (cell 0 == amr_cpat_off). Read
    !! by the pb/mv prolong + ghost-fill so np>=2 couples to the correct coarse rank. Allocated only for non-polytropic QBMM.
    real(stp), allocatable, dimension(:,:,:,:,:) :: amr_cg_pb, amr_cg_mv
    $:GPU_DECLARE(create='[amr_cg_pb, amr_cg_mv]')

    !> L0-AS-BLOCKS SPIKE (l0_ntile > 0, amr off). Feasibility probe for AMReX-style dynamic load balancing: tile the base grid into
    !! l0_ntile**num_dims base-resolution (refinement-ratio-1) blocks and advance each through the SAME swap-based per-block solver
    !! the AMR fine overlay uses, with tile-tile same-level seam halos (s_amr_fine_fine_halo, fmul=1) at interior faces and the
    !! physical BC at domain-edge faces. Correctness bar: l0_ntile>0 must be BYTE-IDENTICAL to l0_ntile=0 (monolithic). Reuses the
    !! full amr_slots/region/owner/seam machinery; the tiles ARE level-1 blocks with amr_ref_ratio 1. Off (l0_ntile=0) => no effect.
    integer :: l0_ntiles_tot = 0  !< total tiles = l0_ntile**num_dims (0 when off)
    integer :: l0_nt(3) = 1       !< tiles per dim (1 in collapsed dims)
    ! Unification pool layout (amr_max_fine, l0_slot_off) lives in m_global_parameters beside amr_num_blocks/amr_max_blocks, so
    ! m_amr_regrid (a separate module) sees it too.
    !> per-dim global periodicity, allreduced from periodic_bc in s_l0_tiles_init (periodic_bc is set on rank 0 only -
    !! s_read_input_file is rank-0-guarded - so it must be made consistent for the wrap-seam decision in f_amr_seam /
    !! s_l0_edge_bc_tile, which every rank must agree on)
    logical :: l0_periodic(3) = .false.
    !> tiles are persistent: L0 seeds them ONCE (first timestep); set at init, cleared after fill
    logical :: l0_tiles_need_fill = .false.
    !> rank owning each tile's L0 STORAGE cells (fixed = init owner); scatter routes the compute-owner's interior back to this rank
    !! when a tile has MIGRATED (owner != l0_owner)
    integer, allocatable :: amr_tile_l0_owner(:)
    !> accumulated MEASURED compute time per owned tile since the last rebalance (GPU-synced wall time); allreduced to a replicated
    !! cost vector that drives s_l0_rebalance, then reset
    real(wp), allocatable :: amr_tile_cost(:)
    !> per-tile exponential moving average of the (replicated) measured cost, smoothed across rebalance windows so GPU
    !! launch-latency timing noise does not drive churn
    real(wp), allocatable :: amr_tile_cost_ema(:)

contains

    !> Build the static refined level-1 block. No-op unless amr. Called after the level-0 grid (x_cb/dx ready) and time-steppers
    !! (sys_size/buff_size set). Per-slot fine arrays allocated lazily (s_amr_reconcile_slots) - only the blocks a rank owns.
    impure subroutine s_initialize_amr_module()

        integer                         :: i, d, islot
        integer                         :: sidx(3), ext(3), maxc_loc(3), bad_loc, bad_glb, fit_d
        integer                         :: blk_lo(3), blk_hi(3)
        type(scalar_field), allocatable :: tmp_cg(:)

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

        amr_dt_fine = 0.5_wp*dt

        ! 2a gate (see amr_prim_batch's declaration). bubbles_euler is checker-prohibited under amr; every listed
        ! feature either extends the conversion write set beyond mom..E or changes its inputs, and falls back to
        ! the per-block conversion.
        ! 2a VERDICT (2026-08-26, ledger 27): OFF by default. The batched conversion kernel is cheap on
        ! device (0.5 s / 5-step probe) and byte-identical, but the per-block PRIM BRIDGE-LOADS that land
        ! its output in the m_rhs scratch cost ~4x wall on the amdflang offload host path (launch/mapping
        ! burden; a flag-off build with the kernel still compiled in prices at baseline, so it is not
        ! codegen poisoning). Bridge loads are exactly what full 2b deletes - partial batching through a
        ! bridge is negative value, so the machinery stays for the 2b store-native consumption experiment
        ! and this gate stays false until then.
        amr_prim_batch = .false.

        ! Fine-block cap = the case amr_max_blocks; the shared pool adds the L0 tile prefix (l0_slot_off, 0 when l0_ntile=0) ahead
        ! of
        ! it, so both AMR fine blocks and any L0 tiles draw from one amr_slots allocation.
        amr_max_fine = amr_max_blocks  ! fine/regrid cap = the case budget
        amr_max_blocks = l0_slot_off + amr_max_fine  ! total shared pool (l0_slot_off=0 when no tiles -> unchanged)

        ! fixed pool of amr_max_blocks slots; init activates exactly one (amr_cur = f_l0_slot(1), the initial fine-block slot);
        ! regrid clusters into up to amr_max_blocks
        allocate (amr_slots(1:amr_max_blocks))
        call s_amr_loc_index_init()
        allocate (amr_region_lo_all(3, amr_max_blocks), amr_region_hi_all(3, amr_max_blocks))
        allocate (amr_isect_lo_all(3, amr_max_blocks), amr_isect_hi_all(3, amr_max_blocks))
        allocate (amr_owns_all(amr_max_blocks))
        allocate (amr_block_owner(amr_max_blocks))
        allocate (amr_owner_cut(0:num_procs - 1)); amr_owner_cut = -1_8
        allocate (amr_fine_cut(0:num_procs - 1,1:max(amr_max_level, 1))); amr_fine_cut = -1_8
        allocate (amr_block_level(amr_max_blocks))
        ! amr_ovl_gather/scatter (the 2D rank lists) are allocated to the computed max overlap in s_amr_build_seam_pairs; only the
        ! per-block counts are sized here.
        allocate (amr_ovl_gather_n(amr_max_blocks), amr_ovl_scatter_n(amr_max_blocks))
        amr_region_lo_all = 0; amr_region_hi_all = 0; amr_isect_lo_all = 0; amr_isect_hi_all = 0; amr_owns_all = .false.
        amr_block_owner = 0
        amr_block_level = 1  ! init default (level-1); regrid re-tags each block's level for nesting
        amr_num_levels = 1
        amr_num_blocks = f_l0_slot(1)
        amr_cur = f_l0_slot(1)

        ! fine-level load balance is capped at min(num blocks, amr_max_blocks) ranks: the SFC map spreads whole blocks, so with
        ! fewer
        ! blocks than ranks some ranks own no fine work. Warn when the pool itself is the limit (raise amr_max_blocks).
        if (proc_rank == 0 .and. num_procs > amr_max_fine) then
            print '(A,I0,A,I0,A)', ' [amr] WARNING: amr_max_blocks (', amr_max_blocks, ') < num_procs (', num_procs, &
                & '): the fine level can occupy at most amr_max_blocks ranks - raise amr_max_blocks for better fine-level balance'
        end if

        ! Lock-step advances every fine block at the coarse dt, but a level-l cell is amr_ref_ratio**l smaller, so its CFL limit is
        ! amr_ref_ratio**amr_max_level tighter than the coarse grid's. The dt (fixed, or the coarse-only cfl_dt estimate) is NOT
        ! scaled for that, so a coarse-CFL dt silently runs the finest block unstable (subcycling instead advances each level at
        ! dt/amr_ref_ratio and is stable by construction). The true CFL is unknown at init, so warn rather than abort - a small
        ! enough dt is valid.
        if (proc_rank == 0 .and. .not. amr_subcycle .and. (amr_ref_ratio > 2 .or. amr_max_level > 1)) then
            print '(A,I0,A)', &
                & ' [amr] WARNING: lock-step (amr_subcycle = F) advances fine blocks at the coarse dt, but the ' &
                & // 'finest cell is amr_ref_ratio**amr_max_level = ', amr_ref_ratio**amr_max_level, &
                & 'x smaller - ensure dt satisfies the FINEST cell CFL (roughly the coarse-stable dt divided by that ' &
                & // 'factor), or enable amr_subcycle, else the fine block may go unstable'
        end if

        ! CONFIGURATION ADVISORIES. Measured 2026-08-02 (3D, np=8, cap-64 optimum); see
        ! @ref amr_per_level_distribution, "Configuration guidance". These are ADVICE, not constraints -
        ! every setting below is legal and sometimes correct, so they warn rather than abort.
        if (proc_rank == 0) then
            ! amr_regrid_int = 0 is STATIC AMR: the block set never changes. That is a legitimate mode
            ! (and the only one supported above amr_max_level = 2), but a user who set `amr = T` expecting
            ! adaptivity gets none, silently.
            if (amr_regrid_int == 0) then
                print '(A)', &
                    & ' [amr] NOTE: amr_regrid_int = 0 - the block set is STATIC and never adapts. ' &
                    & // 'Set amr_regrid_int > 0 (4-8 is a reasonable start) for adaptive refinement.'
            end if
            ! The derived cap is the min-over-ranks local half-extent, so it SHRINKS as ranks are added -
            ! the wrong direction for strong scaling, and it makes the box set (hence the answer, within
            ! tolerance) depend on the rank count. Measured 3.0x slower than a pinned cap of 64 in 3D.
            if (amr_max_grid_size == 0 .and. num_procs > 1) then
                print '(A)', &
                    & ' [amr] NOTE: amr_max_grid_size = 0 derives the block cap from the ' &
                    & // 'decomposition, so it SHRINKS as ranks are added and the box set depends on rank ' &
                    & // 'count. Pinning it (64 measured best in 3D on MI250X, memory-bounded) was 3.0x ' &
                    & // 'faster and makes the box set rank-invariant.'
            end if
            ! Lock-step integrates the COARSE level at the finest-stable dt, i.e. amr_ref_ratio**level
            ! times more often than its own stability requires. Measured 1.55x at amr_max_level = 2.
            ! NB this changes the time integration - it is not a free optimization.
            if (.not. amr_subcycle .and. amr_max_level >= 1 .and. amr_regrid_int > 0) then
                print '(A,I0,A)', ' [amr] NOTE: amr_subcycle = F integrates the coarse level ', amr_ref_ratio**amr_max_level, &
                    & 'x more often than its own CFL requires. ' &
                    & // 'amr_subcycle = T was 1.55x faster per unit physical time (it is a DIFFERENT ' &
                    & // 'time integration, not a drop-in).'
            end if

            ! Frequent regridding is dominated by the per-cell tag sweep, which is flat in box count.
            if (amr_regrid_int > 0 .and. amr_regrid_int < 4) then
                print '(A,I0,A)', ' [amr] NOTE: amr_regrid_int = ', amr_regrid_int, &
                    & ' regrids often; the tag sweep is per-CELL and flat in box count, so interval 8 ' &
                    & // 'measured 1.39x faster. Raise it unless the refined feature moves quickly.'
            end if
        end if

        ! Mirror decomposition: each rank holds the fine cells covering block /\ its own subdomain (np=1: the intersection is the
        ! whole block). buff_size is not available at checker time, so the geometric aborts below must live here.
        sidx = 0; ext = 0
        sidx(1) = start_idx(1); ext(1) = m
        if (n_glb > 0) then; sidx(2) = start_idx(2); ext(2) = n; end if
        if (p_glb > 0) then; sidx(3) = start_idx(3); ext(3) = p; end if
        call s_amr_compute_isect(amr_block_beg, amr_block_end)

        ! the fine ghost shell and reflux outside cells must stay inside the global domain (identical inputs on all ranks; every
        ! rank takes the same branch)
        if (amr_block_beg(1) < buff_size .or. amr_block_end(1) > m_glb - buff_size .or. (n_glb > 0 .and. (amr_block_beg(2) &
            & < buff_size .or. amr_block_end(2) > n_glb - buff_size)) .or. (p_glb > 0 .and. (amr_block_beg(3) < buff_size &
            & .or. amr_block_end(3) > p_glb - buff_size))) then
            call s_mpi_abort('amr block must lie at least buff_size cells inside the domain boundaries')
        end if

        ! Scratch constraint: the fine advance reuses the solver scratch (m_rhs/WENO/Riemann work arrays) and the global coordinate
        ! arrays, all sized to THIS rank's local grid. Fine-level distribution gives a block WHOLE to its owner, so the WHOLE
        ! block's fine extent (2*block-1) must fit every rank's local extent (a big block cannot be whole-owned; it must be split
        ! into <= local-half boxes - the mirror model instead split the block ACROSS ranks). Checked on the replicated block box so
        ! all ranks agree. (np=1: local extent = global, so 2*block-1 <= m_glb always holds.) non-IB: the block is TILED into <=
        ! amr_maxc_fit sub-blocks (each fits every rank's scratch), so no cap is needed. IB keeps a single contiguous block per
        ! body,
        ! so an IB block must itself fit a rank's local half-extent.
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
        ! clamped box satisfies every rank's scratch constraint and can move freely across ranks. That cap SHRINKS as ranks are
        ! added, which tiles a fixed feature into more and more blocks the further you scale - and per-block cost is ~fixed
        ! regardless of block size, so the block count is what costs. It also makes the box set (and so the answer, within
        ! tolerance) depend on the rank count. Setting amr_max_grid_size > 0 pins the cap to an absolute number of coarse cells
        ! instead, exactly like AMReX's max_grid_size: the box set is then IDENTICAL at every rank count.
        amr_maxc_fit = amr_maxc
        do d = 1, num_dims
            call s_mpi_allreduce_integer_min((ext(d) + 1)/amr_ref_ratio, fit_d)
            if (amr_max_grid_size > 0) then
                ! Rank-independent cap, INDEPENDENT of fit_d. The fine advance still borrows this rank's solver scratch, but that
                ! scratch is now sized to the cap rather than to the subdomain (idwbuff_alloc and m/n/p_alloc in
                ! m_global_parameters give it amr_ref_ratio*amr_max_grid_size - 1 fine cells plus the ghost shell), so a block at
                ! the cap fits however small the subdomain becomes. This is what decouples the box set from the rank count:
                ! without it the cap could only ever shrink as ranks grow, which is backwards for strong scaling.
                ! Cost: per-rank scratch is O(cap**num_dims), constant in rank count - bounded by choosing the cap, which is
                ! precisely what this parameter is for.
                amr_maxc_fit(d) = min(amr_maxc(d), amr_max_grid_size)
            else
                amr_maxc_fit(d) = min(amr_maxc(d), fit_d)
            end if
        end do

        ! preallocation cap for MY fine arrays: a block is owned WHOLE, so any rank must hold an entire block. regrid clamps every
        ! box to amr_maxc_fit, so amr_maxc_fit (NOT the global-half amr_maxc) is the true max block a rank can own; sizing to it
        ! right-sizes the fine/coord arrays. At np=1 amr_maxc_fit == amr_maxc, so the sizing (and everything) is unchanged.
        ! NOTE amr_maxc_fit is no longer bounded by the local half-extent when amr_max_grid_size > 0: the solver scratch is sized to
        ! the cap instead (see above), so a rank can own a block LARGER than half its own subdomain. Derived-cap runs
        ! (amr_max_grid_size = 0) still take the min-over-ranks local-half and are unaffected.
        maxc_loc = amr_maxc_fit

        ! max fine extents and buffered bounds for preallocation
        max_f1 = amr_ref_ratio*maxc_loc(1) - 1
        max_f2 = 0; max_f3 = 0
        if (n_glb > 0) max_f2 = amr_ref_ratio*maxc_loc(2) - 1
        if (p_glb > 0) max_f3 = amr_ref_ratio*maxc_loc(3) - 1

        amr_seam_pairs_dirty = .true.; amr_seam_pairs_nblk = -1  ! force a seam-list build on the first fine-fine halo
        amr_mesh_epoch = amr_mesh_epoch + 1
        mbuf1_lo = -buff_size; mbuf1_hi = max_f1 + buff_size
        mbuf2_lo = 0; mbuf2_hi = 0; mbuf3_lo = 0; mbuf3_hi = 0
        if (n_glb > 0) then; mbuf2_lo = -buff_size; mbuf2_hi = max_f2 + buff_size; end if
        if (p_glb > 0) then; mbuf3_lo = -buff_size; mbuf3_hi = max_f3 + buff_size; end if
        ! with tiles, s_l0_tiles_init's mbuf UNION below may still enlarge these - the scratch waits for it (see s_amr_scr_init)
        if (l0_ntile == 0) call s_amr_scr_init()

        ! MEMORY DEMAND, reported not guessed. There is no portable way to ask how much device (or
        ! host) memory is available - hipMemGetInfo / cudaMemGetInfo / nothing-on-CPU across four
        ! compilers and three offload backends - so do NOT try to pick a cap from a memory budget.
        ! What IS exactly known here is the DEMAND: a block costs 2 per-slot field families (q_cons,
        ! q_cons_stor; q_prim/rhs are POOLED - one shared scratch pair, not per block) plus, under
        ! amr_subcycle only, 2 more in the flat store (amr_gst_a/amr_gst_b, sized per LOCAL slot)
        ! x sys_size arrays on the mbuf extents. Print it
        ! and let the reader compare against hardware they know.
        !
        ! Why this matters more than a default: the OPTIMAL cap is set by the largest slot that
        ! fits, and slot volume goes as cap**num_dims - so the best cap measured 1024 in 2D and 64
        ! in 3D, a 16x difference, while the SLOT VOLUMES agreed to 1.7x. One number cannot serve
        ! both dimensions; the volume is the invariant. Exceeding it aborts inside
        ! __tgt_target_data_begin_mapper, which PRESENTS AS A HANG (one rank dies, the rest block in
        ! MPI), so a silent over-large cap is expensive to diagnose.
        if (proc_rank == 0) then
            block
                real(wp) :: slot_gib, cells, nfam
                cells = real(mbuf1_hi - mbuf1_lo + 1, wp)
                if (n_glb > 0) cells = cells*real(mbuf2_hi - mbuf2_lo + 1, wp)
                if (p_glb > 0) cells = cells*real(mbuf3_hi - mbuf3_lo + 1, wp)
                nfam = 2._wp; if (amr_subcycle) nfam = 4._wp
                slot_gib = cells*real(sys_size, wp)*nfam*real(storage_size(1._wp)/8, wp)/1024._wp**3
                print '(A,I0,A,I0,A,ES10.3,A,F8.3,A)', ' [amr] per-block slot: ', nint(cells), ' cells x sys_size x ', &
                    & nint(nfam), ' fields = ', cells*real(sys_size, wp)*nfam, ' words (', slot_gib, ' GiB per owned block)'
                print '(A,F9.2,A,I0,A)', ' [amr]   worst case if one rank owned every block: ', slot_gib*real(amr_max_blocks, &
                    & wp), ' GiB (amr_max_blocks = ', amr_max_blocks, '). Typical is amr_max_blocks/num_procs blocks per rank.'
            end block
        end if

        ! bounce buffers for copy-based coord swap (GPU-safe; same bounds as the base-level global arrays, which are sized on
        ! *_alloc - these are whole-array assigned to/from x_cb etc., so the shapes must agree)
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
        if (cyl_coord .and. n_glb > 0) then
            @:ALLOCATE(amr_rvw(0:max_f2))
        end if

        ! Grid uniformity policy. Both spacing-uniformity consumers are handled exactly: fine-block ghost-shell coordinates extend
        ! by exact parent-cell bisection (reads sw_*_cb), and the spacing-dependent WENO reconstruction coefficients are recomputed
        ! for the active grid on every swap/restore when the grid is nonuniform anywhere (stretched grids, or 2D-axisymmetric's
        ! half-width axis cell dy(0) = dy/2). Tolerance is epsilon-scaled: an absolute 1e-12 would sit below single-precision grid
        ! roundoff and classify every grid as stretched (spuriously tripping the stretched-combo gates). On uniform grids the flag
        ! stays false and behavior is bit-identical to the reuse path. The stretch_* flags are pre_process-only, so the grid itself
        ! is checked (this also catches externally generated grids).
        if (maxval(dx(0:m)) - minval(dx(0:m)) > 1.e3_wp*epsilon(1._wp)*maxval(dx(0:m))) then
            amr_weno_coef_recompute = .true.; amr_grid_stretched = .true.
        end if
        if (n_glb > 0) then
            ! interior nonuniformity is stretching; a lone dy(0) deviation is stretching only when it is NOT the axisymmetric
            ! half-width axis cell
            if (n > 0 .and. maxval(dy(1:n)) - minval(dy(1:n)) > 1.e3_wp*epsilon(1._wp)*maxval(dy(1:n))) then
                amr_weno_coef_recompute = .true.; amr_grid_stretched = .true.
            end if
            if (abs(dy(0) - dy(min(1, n))) > 1.e3_wp*epsilon(1._wp)*dy(min(1, n))) then
                amr_weno_coef_recompute = .true.
                if (.not. cyl_coord) amr_grid_stretched = .true.
            end if
        end if
        if (p_glb > 0) then
            if (maxval(dz(0:p)) - minval(dz(0:p)) > 1.e3_wp*epsilon(1._wp)*maxval(dz(0:p))) then
                amr_weno_coef_recompute = .true.; amr_grid_stretched = .true.
            end if
        end if
        if (weno_order == 1 .or. igr) amr_weno_coef_recompute = .false.  ! order 1 / IGR: no grid-dependent WENO coefficients

        ! persistent global coarse boundaries: the fine-distribution owner rebuilds whole-block fine coordinates from these (needed
        ! once the fine level is decoupled from the coarse decomposition; harmless otherwise)
        call s_amr_build_global_cb()
        ! Fail closed on stretched grid + Lagrangian/IB-dynamic-regrid. TWO independent blockers (both confirmed by experiment):
        !  (1) the position->global-cell-index conversions here use int((x-beg)/dx(0)), inexact on a stretched grid and
        !      rank-inconsistent (dx(0) is rank-local). FIXABLE: assemble the global cell-boundary arrays (allreduce-MAX of owned
        !      boundaries) and bisection-search them - correct on any grid, identical on every rank.
        !  (2) THE HARDER BLOCKER: IB/Lagrangian floor buff_size (10/6), but s_amr_recompute_weno_coefs (armed only on nonuniform
        !      grids) indexes poly_coef_cb* over -buff_size:m+buff_size while m_weno sized those arrays with a smaller buff_size at
        !      init -> OOB write in s_compute_weno_coefficients (gfortran bounds trap at m_weno.fpp weno5 branch). Needs the WENO
        !      coefficient arrays sized to the final buff_size (or the recompute clamped to the module's true bounds) before this
        !      gate can lift. Fix (1) alone is insufficient.
        if (amr_grid_stretched .and. (bubbles_lagrange .or. (ib .and. amr_regrid_int > 0))) then
            call s_mpi_abort('amr on a stretched grid does not support ' &
                             & // 'Lagrangian bubbles or dynamic regrid with immersed bodies: their ' &
                             & // 'position-to-cell-index conversions assume uniform spacing')
        end if

        ! per-slot field arrays are allocated by s_amr_alloc_slot / freed by s_amr_free_slot (sized to the max buffered block). The
        ! lazy owned-only reconcile that keeps a rank's fine memory ~1/num_procs of the pool follows. The QBMM RHS scratch
        ! (amr_rhs_pb_f/mv_f) is single - allocate once.
        allocate (amr_slot_live(amr_max_blocks)); amr_slot_live = .false.
        if (qbmm .and. .not. polytropic) then
            @:ALLOCATE(amr_rhs_pb_f(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi, 1:nnode, 1:nb))
            @:ALLOCATE(amr_rhs_mv_f(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi, 1:nnode, 1:nb))
        end if
        ! per-slot field arrays are allocated lazily by s_amr_reconcile_slots once ownership is known (after the block setup +
        ! s_amr_assign_block_owners below), so a rank holds only its owned blocks' fine arrays - not all amr_max_blocks slots.

        ! fine-level distribution: coarse-patch gather buffer (see decl). Sized to the largest block's coarse footprint (block
        ! coarse
        ! cells + 2*nmar halo, block-local frame). Device-mapped so the runtime ghost-fill reads it on the owner.
        amr_cpat_mar = (buff_size + amr_ref_ratio - 1)/amr_ref_ratio + 1
        amr_cpat_hi = 0
        amr_cpat_hi(1) = maxc_loc(1) - 1 + 2*amr_cpat_mar
        if (n_glb > 0) amr_cpat_hi(2) = maxc_loc(2) - 1 + 2*amr_cpat_mar
        if (p_glb > 0) amr_cpat_hi(3) = maxc_loc(3) - 1 + 2*amr_cpat_mar
        ! CCE OpenMP-offload leaves a bare module-scope derived-type (scalar_field) allocatable's descriptor uninitialized, so a
        ! direct allocate(amr_cg(1:sys_size)) aborts with lib-4425 at program start (verified by an early module-init probe; a LOCAL
        ! scalar_field array and a GPU_DECLARE'd module one like q_prim_vf both allocate fine - only this bare module array does
        ! not). Allocate a local, which gets a valid descriptor, and hand it to the module variable via move_alloc, then map.
        ! OpenACC is unaffected but takes the same path correctly.
        allocate (tmp_cg(1:sys_size))
        call move_alloc(tmp_cg, amr_cg)
        $:GPU_ENTER_DATA(create='[amr_cg]')
        do i = 1, sys_size
            @:ALLOCATE(amr_cg(i)%sf(0:amr_cpat_hi(1), 0:amr_cpat_hi(2), 0:amr_cpat_hi(3)))
            amr_cg(i)%sf = 0._stp  ! padding beyond a block's valid patch extent is never read; keep it finite for the device copy
            @:ACC_SETUP_SFs(amr_cg(i))
        end do

        ! non-polytropic QBMM: gathered coarse pb/mv patch (analogue of amr_cg, same footprint + trailing (nnode, nb) dims). Plain
        ! 5D
        ! arrays (amr_rhs_pb_f idiom): the module GPU_DECLARE + @:ALLOCATE handle device mapping - no @:ACC_SETUP_SFs.
        if (qbmm .and. .not. polytropic) then
            @:ALLOCATE(amr_cg_pb(0:amr_cpat_hi(1), 0:amr_cpat_hi(2), 0:amr_cpat_hi(3), 1:nnode, 1:nb))
            @:ALLOCATE(amr_cg_mv(0:amr_cpat_hi(1), 0:amr_cpat_hi(2), 0:amr_cpat_hi(3), 1:nnode, 1:nb))
            amr_cg_pb = 0._stp; amr_cg_mv = 0._stp
        end if

        ! the coarse decomposition (each rank's coarse start_idx + local m/n/p) is a structured cartesian split, computed O(1) per
        ! rank by s_amr_rank_decomp - no replicated table, no allgather. Validate the formula against this rank's actual values.
        call s_amr_validate_decomp()

        ! per-slot fine-grid IB marker fields (static-body AMR); sized to the same max buffered fine extents as q_cons so the fine
        ! IB pipeline can resolve the body on the block
        if (ib) call s_ibm_alloc_fine(amr_max_blocks, mbuf1_lo, mbuf1_hi, mbuf2_lo, mbuf2_hi, mbuf3_lo, mbuf3_hi)

        ! set geometry (region, m/n/p, idwbuff, coordinates) for the initial block (amr_cur = f_l0_slot(1), the initial fine-block
        ! slot). Under dynamic regrid with bodies
        ! the initial block gets the same body-containment expansion regrid boxes get (the moving-body containment guard requires it
        ! from step 1); for a static block (amr_regrid_int = 0) the user's placement is authoritative. max_grid_size tiling: the
        ! initial block splits into <= amr_maxc_fit sub-blocks (at np=1 amr_maxc_fit == amr_maxc so a normal block stays a single
        ! tile - unchanged), one per slot; IB keeps a single contiguous block.
        blk_lo = amr_block_beg; blk_hi = amr_block_end
        if (ib .and. amr_regrid_int > 0) call s_amr_expand_box_over_bodies(blk_lo, blk_hi)
        block
            type(t_box), allocatable :: tiled(:)
            integer                  :: nt, capt, kk
            allocate (tiled(amr_max_blocks)); nt = 0; capt = 0
            if (ib) then
                nt = 1; tiled(1)%lo = blk_lo; tiled(1)%hi = blk_hi
            else
                call s_amr_tile_box(blk_lo, blk_hi, tiled, nt, amr_max_fine, capt)
            end if
            amr_num_blocks = f_l0_slot(nt)  ! fine blocks occupy [l0_slot_off+1 .. l0_slot_off+nt] in the shared pool
            ! set block regions FIRST so the owner assignment (reads amr_region_*_all) runs BEFORE the owner-dependent geometry -
            ! else s_set_amr_fine_geometry would size the whole-block owner from a stale (default) amr_block_owner
            do kk = 1, nt
                amr_region_lo_all(:,f_l0_slot(kk)) = tiled(kk)%lo; amr_region_hi_all(:,f_l0_slot(kk)) = tiled(kk)%hi
            end do
            call s_amr_assign_block_owners()  ! assign each block's single owner rank (fine-dist map)
            call s_amr_reconcile_slots()  ! allocate this rank's owned initial blocks (owner-guarded geometry writes below)
            do kk = 1, nt
                amr_cur = f_l0_slot(kk)
                call s_set_amr_fine_geometry(tiled(kk)%lo, tiled(kk)%hi)
            end do
            call s_amr_reduce_xchg_flag()
            call s_amr_select_slot(f_l0_slot(1))  ! refresh the per-block mirrors (geometry loop left them on the last tile)
            deallocate (tiled)
        end block

        ! plan-based exchange (I0): per-family tag bases sit above the legacy per-box tag space so both can coexist
        ! while families convert one increment at a time.
        block
            integer :: f
            do f = 1, size(amr_tag_base)
                amr_tag_base(f) = amr_max_blocks + 100*f
            end do
        end block
#ifdef MFC_MPI
        block
            integer(kind=MPI_ADDRESS_KIND) :: tag_ub
            logical                        :: tag_ub_set
            integer                        :: ierr
            call MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, tag_ub, tag_ub_set, ierr)
            @:ASSERT(tag_ub_set, "MPI_TAG_UB attribute unavailable")
            @:ASSERT(amr_tag_base(size(amr_tag_base)) + 100 <= tag_ub, &
                     & "AMR tag space exceeds MPI_TAG_UB: amr_max_blocks is too large for this MPI's tag range")
        end block
#endif

    end subroutine s_initialize_amr_module

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

    !> Fine cell coordinates of block k in dimension d, rebuilt from the GLOBAL coarse boundaries gcb by replaying k's ancestor
    !! chain - touching no other block's slot arrays. A level-l block's grid is l nested midpoint subdivisions of the L0 boundaries,
    !! and every box in the chain is known from REPLICATED metadata (amr_region_*_all + the global amr_ref_ratio), so any rank can
    !! reproduce it. Bit-identical to bisecting the parent's stored coords, because that array is itself the same subdivision.
    !!
    !! This exists because the direct form - bisecting amr_slots(parent)%x_cb - reads the PARENT's slot, which is
    !! ALLOCATED ONLY ON THE PARENT'S OWNER. Under tower co-location the child's owner was always the parent's owner too, so it
    !! worked; under per-level distribution a level>=2 block can be owned by a rank holding no part of its parent, and bisecting an
    !! unallocated array there produced garbage cell widths and NaNs a few steps later.
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

    !> Compute this rank's per-dim intersection of the box lo:hi with its subdomain (GLOBAL indices, mirrored to amr_isect_lo/hi)
    !! and whether it holds fine cells (amr_rank_owns_block: nonempty in all active dims). Must be called with the COARSE grid state
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

    !> Fine-level distribution: assemble the current block's coarse patch on its owner. The patch covers GLOBAL coarse cells
    !! region_lo-amr_cpat_mar : region_hi+amr_cpat_mar (the full reach of every prolongation/ghost-fill stencil) for all sys_size
    !! variables, stored in amr_cg in a block-LOCAL frame (cell 0 == global amr_cpat_off). POINT-TO-POINT: the owner receives the
    !! patch cells it does not hold from exactly the (SFC-local) coarse-owners that hold them - each rank's contribution is the
    !! patch intersected with its contiguous owned coarse range (s_amr_rank_coarse_range, = the f_amr_own_coarse set, computed from
    !! the cartesian decomposition). Non-participants send/recv nothing (no global collective). At np=1 the owner just copies its
    !! own coarse over the patch, bit-for-bit. Runtime (pull_host) packs/unpacks the overlap boxes on the DEVICE (q_coarse
    !! device-current with valid ghosts); init/regrid fills from the host (host-current with valid ghosts). Packed data is wp, cast
    !! to stp into amr_cg (identity for stp coarse), device-current on exit. INVARIANT: "coarse" here means the block's PARENT level
    !! (level l-1), NOT the base grid (level 0). For a level-1 block the parent IS L0, but a level>=2 block folds to/from its parent
    !! block's fine array; the C<->F prolong/restrict/gather routines all operate in the parent-fine frame, not the L0 frame. TWIN
    !! s_amr_gather_coarse_patch_pbmv (q<->pb/mv): same P2P skeleton (rank-range, intersection, pack/send/recv/unpack) and
    !! patch-local frame - keep lockstep.
    !> Make room for one more pending gather send, draining the pool first if it is full. Draining is a WAITALL, so the pool size
    !! sets how far a contributing rank may run ahead of the owners.
    impure subroutine s_amr_gsnd_reserve(slotsz)

        integer, intent(in) :: slotsz

        if (.not. allocated(amr_gsnd_pool)) then
            allocate (amr_gsnd_pool(slotsz, amr_gsnd_max), amr_gsnd_req(amr_gsnd_max))
            amr_gsnd_n = 0
        else if (size(amr_gsnd_pool, 1) < slotsz) then
            call s_amr_gather_send_flush()  ! outstanding sends reference the old buffer - complete them before resizing
            deallocate (amr_gsnd_pool)
            allocate (amr_gsnd_pool(slotsz, amr_gsnd_max))
        end if
        if (amr_gsnd_n >= amr_gsnd_max) call s_amr_gather_send_flush()

    end subroutine s_amr_gsnd_reserve

    !> Complete every pending gather send. MUST be called before the send buffers are reused or the routine returns to a caller that
    !! will free them - an ISEND whose buffer is overwritten in flight silently corrupts the receiver's patch.
    impure subroutine s_amr_gather_send_flush()

        integer :: ierr

        if (amr_gsnd_n == 0) return
#ifdef MFC_MPI
        call MPI_WAITALL(amr_gsnd_n, amr_gsnd_req(1:amr_gsnd_n), MPI_STATUSES_IGNORE, ierr)
#endif
        amr_gsnd_n = 0

    end subroutine s_amr_gather_send_flush

    !> Gather-batching step 1: derive the ENTIRE rebuild gather message set up front - per level-1 box its contributor ranks and
    !! message sizes, per level>=2 box its parent source and size - from the same replicated caches the per-box path reads
    !! (amr_region_*_all, amr_ovl_gather, amr_block_owner, rank coarse ranges, s_amr_parent_foot). Exchange behavior is UNCHANGED by
    !! this step: the per-box gather asserts against the plan, box by box (see amr_gpl_* declarations). Caller
    !! (s_amr_regrid_rebuild_slots) clears amr_gpl_valid when its box loop ends.
    impure subroutine s_amr_build_gather_plan(nboxes)

        integer, intent(in) :: nboxes
        integer             :: k, ks, idx, r, nsrc, mo, pblk
        integer             :: v1hi, v2hi, v3hi, plo(3), phi(3), crlo(3), crhi(3), bl(3), bh(3), w(3)

        ! the per-box path lazily rebuilds the overlap lists inside the first gather; force the SAME rebuild here so the plan
        ! and the boxes read identical lists

        if (amr_seam_pairs_dirty .or. amr_seam_pairs_nblk /= amr_num_blocks) call s_amr_build_seam_pairs()
        mo = size(amr_ovl_gather, 1)
        if (allocated(amr_gpl_src)) then
            if (size(amr_gpl_src, 1) < mo) deallocate (amr_gpl_src, amr_gpl_sz)
        end if
        if (.not. allocated(amr_gpl_nsrc)) allocate (amr_gpl_nsrc(amr_max_blocks), amr_gpl_psrc(amr_max_blocks), &
            & amr_gpl_psz(amr_max_blocks))
        if (.not. allocated(amr_gpl_src)) allocate (amr_gpl_src(mo, amr_max_blocks), amr_gpl_sz(mo, amr_max_blocks))
        do k = 1, nboxes
            ks = f_l0_slot(k)
            amr_gpl_nsrc(ks) = 0; amr_gpl_psrc(ks) = -1; amr_gpl_psz(ks) = 0
            if (amr_block_level(ks) >= 2) then
                pblk = f_amr_parent_block(ks)
                if (amr_block_owner(pblk) /= amr_block_owner(ks)) then
                    call s_amr_parent_foot(ks, pblk, plo, phi)
                    w = 0
                    w(1) = (phi(1) - plo(1)) + 2*amr_cpat_mar
                    if (n_glb > 0) w(2) = (phi(2) - plo(2)) + 2*amr_cpat_mar
                    if (p_glb > 0) w(3) = (phi(3) - plo(3)) + 2*amr_cpat_mar
                    amr_gpl_psrc(ks) = amr_block_owner(pblk)
                    amr_gpl_psz(ks) = sys_size*(w(1) + 1)*(w(2) + 1)*(w(3) + 1)
                end if
            else
                ! level-1 patch box: same arithmetic as the gather's patch-frame block (collapsed dims stay 0)
                plo = 0
                plo(1) = amr_region_lo_all(1, ks) - amr_cpat_mar
                if (n_glb > 0) plo(2) = amr_region_lo_all(2, ks) - amr_cpat_mar
                if (p_glb > 0) plo(3) = amr_region_lo_all(3, ks) - amr_cpat_mar
                v1hi = (amr_region_hi_all(1, ks) - amr_region_lo_all(1, ks)) + 2*amr_cpat_mar
                v2hi = 0; v3hi = 0
                if (n_glb > 0) v2hi = (amr_region_hi_all(2, ks) - amr_region_lo_all(2, ks)) + 2*amr_cpat_mar
                if (p_glb > 0) v3hi = (amr_region_hi_all(3, ks) - amr_region_lo_all(3, ks)) + 2*amr_cpat_mar
                phi(1) = plo(1) + v1hi; phi(2) = plo(2) + v2hi; phi(3) = plo(3) + v3hi
                nsrc = 0
                do idx = 1, amr_ovl_gather_n(ks)
                    r = amr_ovl_gather(idx, ks)
                    if (r == amr_block_owner(ks)) cycle
                    call s_amr_rank_coarse_range(r, crlo, crhi)
                    call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
                    nsrc = nsrc + 1
                    amr_gpl_src(nsrc, ks) = r
                    amr_gpl_sz(nsrc, ks) = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                end do
                amr_gpl_nsrc(ks) = nsrc
            end if
        end do
        amr_gpl_valid = .true.

    end subroutine s_amr_build_gather_plan

    !> Step-2 phase A: pre-post every recv this rank needs for boxes [c_lo, c_hi] - level-1 contributor slices and split level>=2
    !! parent patches - straight from the plan into the flat chunk pool, tag = slot, appended in box order so box k's requests are
    !! one contiguous run. Ownership from amr_block_owner ONLY (amr_owns_all / amr_rank_owns_block still mirror the previous
    !! generation until the consume phase's geometry call). Contains no MPI waits; reallocating the pool here is safe because every
    !! recv posted for the previous chunk was completed inside that chunk's consume phase.
    impure subroutine s_amr_gather_chunk_post(c_lo, c_hi)

        integer, intent(in) :: c_lo, c_hi
        integer             :: k, ks, cb, idx, need, nreq, off, ierr

        @:ASSERT(amr_gpl_valid, "chunk gather: no plan")
        call s_phase_tic(PH_RBPOST)
        need = 0; nreq = 0
        do k = c_lo, c_hi
            ks = f_l0_slot(k)
            if (amr_block_owner(ks) /= proc_rank) cycle
            ! + XA_NH per message: the I1b identity header rides ahead of each payload (zero in production)
            if (amr_block_level(ks) >= 2) then
                if (amr_gpl_psrc(ks) >= 0) then
                    need = need + amr_gpl_psz(ks) + XA_NH; nreq = nreq + 1
                end if
            else
                do idx = 1, amr_gpl_nsrc(ks)
                    need = need + amr_gpl_sz(idx, ks) + XA_NH
                end do
                nreq = nreq + amr_gpl_nsrc(ks)
            end if
        end do
        if (allocated(amr_gcr_pool)) then
            if (size(amr_gcr_pool) < need) deallocate (amr_gcr_pool)
        end if
        if (need > 0 .and. .not. allocated(amr_gcr_pool)) allocate (amr_gcr_pool(need))
        if (allocated(amr_gcr_req)) then
            if (size(amr_gcr_req) < nreq) deallocate (amr_gcr_req, amr_gcr_off)
        end if
        if (nreq > 0 .and. .not. allocated(amr_gcr_req)) allocate (amr_gcr_req(nreq), amr_gcr_off(nreq))

        amr_gcr_n = 0; off = 0
        amr_gcr_nr(:) = 0; amr_gcr_sent(:) = .false.
        do k = c_lo, c_hi
            ks = f_l0_slot(k)
            cb = k - c_lo + 1
            amr_gcr_r0(cb) = amr_gcr_n + 1
            if (amr_block_owner(ks) /= proc_rank) cycle
#ifdef MFC_MPI
            if (amr_block_level(ks) >= 2) then
                if (amr_gpl_psrc(ks) >= 0) then
                    amr_gcr_n = amr_gcr_n + 1
                    amr_gcr_off(amr_gcr_n) = off
                    call s_xa_rec(XA_F2_RCV, 2, amr_gpl_psz(ks), ks)
                    call MPI_IRECV(amr_gcr_pool(off + 1), amr_gpl_psz(ks) + XA_NH, mpi_p, amr_gpl_psrc(ks), ks, MPI_COMM_WORLD, &
                                   & amr_gcr_req(amr_gcr_n), ierr)
                    off = off + amr_gpl_psz(ks) + XA_NH
                    amr_gcr_nr(cb) = 1
                end if
            else
                do idx = 1, amr_gpl_nsrc(ks)
                    amr_gcr_n = amr_gcr_n + 1
                    amr_gcr_off(amr_gcr_n) = off
                    call s_xa_rec(XA_F1_RCV, 2, amr_gpl_sz(idx, ks), ks)
                    call MPI_IRECV(amr_gcr_pool(off + 1), amr_gpl_sz(idx, ks) + XA_NH, mpi_p, amr_gpl_src(idx, ks), ks, &
                                   & MPI_COMM_WORLD, amr_gcr_req(amr_gcr_n), ierr)
                    off = off + amr_gpl_sz(idx, ks) + XA_NH
                end do
                amr_gcr_nr(cb) = amr_gpl_nsrc(ks)
            end if
#endif
        end do
        call s_phase_toc(PH_RBPOST)

    end subroutine s_amr_gather_chunk_post

    !> Step-2 phase B: issue this rank's sends for boxes [c_lo, c_hi]. Level-1: pack the host slice of q_coarse (host is truth
    !! during rebuild) and ISEND through the deferred pool - the per-box contributor path verbatim, driven by the plan. Level>=2
    !! split pairs: send ONLY when the parent was consumed in an earlier chunk (pblk < f_l0_slot(c_lo), monotone slot map) - a
    !! same-chunk parent's new-generation store is not built until its own consume iteration, so that send stays at the child's
    !! consume position, where parents-first ordering guarantees the parent is complete. Geometry from the replicated caches only:
    !! no s_set_amr_fine_geometry swap, no amr_cur.
    impure subroutine s_amr_gather_chunk_send(q_coarse, c_lo, c_hi)

        type(scalar_field), dimension(sys_size), intent(in) :: q_coarse
        integer, intent(in)                                 :: c_lo, c_hi
        integer                                             :: k, ks, cb, idx, i, g1, g2, g3, o1, o2, o3, boxsz, maxsz, pblk, ierr
        integer                                             :: v1hi, v2hi, v3hi, plo(3), phi(3), crlo(3), crhi(3), bl(3), bh(3)
        logical                                             :: contrib

        o1 = start_idx(1); o2 = 0; o3 = 0
        if (n_glb > 0) o2 = start_idx(2)
        if (p_glb > 0) o3 = start_idx(3)
        do k = c_lo, c_hi
            ks = f_l0_slot(k)
            cb = k - c_lo + 1
            if (amr_block_level(ks) >= 2) then
                if (amr_gpl_psrc(ks) < 0) cycle  ! co-located: no message
                pblk = f_amr_parent_block(ks)
                if (amr_block_owner(pblk) /= proc_rank) cycle  ! not the sender
                if (pblk >= f_l0_slot(c_lo)) cycle  ! same-chunk parent: send at the child's consume position
                call s_phase_tic(PH_PGSEND)
                call s_amr_gather_from_parent_field_cons(ks, pblk, amr_loc_of(pblk), .true.)
                call s_phase_toc(PH_PGSEND)
                amr_gcr_sent(cb) = .true.
            else
                if (amr_block_owner(ks) == proc_rank) cycle  ! the owner receives
                contrib = .false.
                do idx = 1, amr_gpl_nsrc(ks)
                    if (amr_gpl_src(idx, ks) == proc_rank) contrib = .true.
                end do
                if (.not. contrib) cycle
                ! same patch-frame arithmetic as the plan builder and the per-box gather
                plo = 0
                plo(1) = amr_region_lo_all(1, ks) - amr_cpat_mar
                if (n_glb > 0) plo(2) = amr_region_lo_all(2, ks) - amr_cpat_mar
                if (p_glb > 0) plo(3) = amr_region_lo_all(3, ks) - amr_cpat_mar
                v1hi = (amr_region_hi_all(1, ks) - amr_region_lo_all(1, ks)) + 2*amr_cpat_mar
                v2hi = 0; v3hi = 0
                if (n_glb > 0) v2hi = (amr_region_hi_all(2, ks) - amr_region_lo_all(2, ks)) + 2*amr_cpat_mar
                if (p_glb > 0) v3hi = (amr_region_hi_all(3, ks) - amr_region_lo_all(3, ks)) + 2*amr_cpat_mar
                phi(1) = plo(1) + v1hi; phi(2) = plo(2) + v2hi; phi(3) = plo(3) + v3hi
                call s_amr_rank_coarse_range(proc_rank, crlo, crhi)
                call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
                boxsz = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                maxsz = sys_size*(v1hi + 1)*(v2hi + 1)*(v3hi + 1)
                call s_phase_tic(PH_RBRSV)
                call s_amr_gsnd_reserve(maxsz + XA_NH)
                call s_phase_toc(PH_RBRSV)
                amr_gsnd_n = amr_gsnd_n + 1
                call s_phase_tic(PH_RBPACK)
                if (XA_NH > 0) call s_xa_hdr_pack(amr_gsnd_pool(:,amr_gsnd_n), XA_F1_SND, ks, bl, bh)
                idx = XA_NH
                do i = 1, sys_size
                    do g3 = bl(3), bh(3)
                        do g2 = bl(2), bh(2)
                            do g1 = bl(1), bh(1)
                                idx = idx + 1
                                amr_gsnd_pool(idx, amr_gsnd_n) = real(q_coarse(i)%sf(g1 - o1, g2 - o2, g3 - o3), wp)
                            end do
                        end do
                    end do
                end do
                call s_phase_toc(PH_RBPACK)
#ifdef MFC_MPI
                call s_phase_tic(PH_RBSEND)
                call s_xa_rec(XA_F1_SND, 1, boxsz, ks)
                call MPI_ISEND(amr_gsnd_pool(1, amr_gsnd_n), boxsz + XA_NH, mpi_p, amr_block_owner(ks), ks, MPI_COMM_WORLD, &
                               & amr_gsnd_req(amr_gsnd_n), ierr)
                call s_phase_toc(PH_RBSEND)
#endif
            end if
        end do

    end subroutine s_amr_gather_chunk_send

    !> Step-2 phase C: the per-box gather body with the exchange already in flight - fill amr_cg for the current box (amr_cur,
    !! geometry already set by the caller) from the own slice plus the chunk pool's pre-posted recvs. Level-1 owner: own-box host
    !! copy, one WAITALL on this box's contiguous request run, host unpack per contributor (plan order = posting order), device
    !! push. Level>=2: co-located parent = local device copy; split parent = the parent owner packs and sends HERE when the parent
    !! shares this chunk (amr_gcr_sent marks the ones phase B already covered), the child owner waits and device-unpacks its single
    !! pre-posted recv. Every rank passes through for every box (the caller's owner-cycle comes after), which is what makes the
    !! request arrays reusable next chunk.
    impure subroutine s_amr_gather_consume_box(q_coarse, k, c_lo)

        type(scalar_field), dimension(sys_size), intent(in) :: q_coarse
        integer, intent(in) :: k, c_lo
        integer :: cb, idx, i, r, g1, g2, g3, o1, o2, o3, boxsz, pblk, w1, w2, w3, ierr, r0, nr, off
        integer :: v1hi, v2hi, v3hi, plo(3), phi(3), crlo(3), crhi(3), bl(3), bh(3)

        cb = k - c_lo + 1
        r0 = amr_gcr_r0(cb); nr = amr_gcr_nr(cb)

        if (amr_block_level(amr_cur) >= 2) then
            call s_phase_tic(PH_PGALL)
            pblk = f_amr_parent_block(amr_cur)
            ! the deferred same-chunk send below reads the parent's store, valid ONLY because parents-first ordering already
            ! consumed the parent (D1 in amr_regrid_gather_batching.md) - trip immediately if the ordering is ever violated
            @:ASSERT(pblk < f_l0_slot(k), "chunk gather: parent box not before child")
            if (amr_gpl_psrc(amr_cur) < 0) then
                ! co-located: the owner's local device copy (the field routine detects co-location itself)
                if (amr_block_owner(amr_cur) == proc_rank) then
                    call s_phase_tic(PH_PGSEND)
                    call s_amr_gather_from_parent_field_cons(amr_cur, pblk, amr_loc_of(pblk), .true.)
                    call s_phase_toc(PH_PGSEND)
                end if
            else if (amr_block_owner(pblk) == proc_rank) then
                ! split, parent side: a same-chunk parent could not be packed in the send phase (its store was unbuilt);
                ! parents-first ordering means it is complete now
                if (.not. amr_gcr_sent(cb)) then
                    call s_phase_tic(PH_PGSEND)
                    call s_amr_gather_from_parent_field_cons(amr_cur, pblk, amr_loc_of(pblk), .true.)
                    call s_phase_toc(PH_PGSEND)
                end if
            else if (amr_block_owner(amr_cur) == proc_rank) then
                ! split, child side: wait on the pre-posted parent patch and unpack on the device
                call s_amr_parent_foot(amr_cur, pblk, plo, phi)
                amr_cpat_off = 0
                amr_cpat_off(1) = plo(1) - amr_cpat_mar
                if (n_glb > 0) amr_cpat_off(2) = plo(2) - amr_cpat_mar
                if (p_glb > 0) amr_cpat_off(3) = plo(3) - amr_cpat_mar
                w1 = (phi(1) - plo(1)) + 2*amr_cpat_mar
                w2 = 0; w3 = 0
                if (n_glb > 0) w2 = (phi(2) - plo(2)) + 2*amr_cpat_mar
                if (p_glb > 0) w3 = (phi(3) - plo(3)) + 2*amr_cpat_mar
#ifdef MFC_MPI
                call s_phase_tic(PH_PGRECV)
                call MPI_WAITALL(nr, amr_gcr_req(r0:r0 + nr - 1), MPI_STATUSES_IGNORE, ierr)
                call s_phase_toc(PH_PGRECV)
                off = amr_gcr_off(r0)
                boxsz = amr_gpl_psz(amr_cur)
                if (XA_NH > 0) call s_xa_hdr_check(amr_gcr_pool(off + 1:off + XA_NH), XA_F2_SND, amr_cur, plo, phi)
                call s_amr_unpack_parent_patch_device(w1, w2, w3, amr_gcr_pool(off + XA_NH + 1:off + XA_NH + boxsz), .true.)
#endif
            end if
            call s_phase_toc(PH_PGALL)
            return
        end if

        ! level-1: same patch frame and own fill as the per-box gather; the recvs are already posted
        amr_cpat_off = 0
        amr_cpat_off(1) = amr_region_lo_all(1, amr_cur) - amr_cpat_mar
        if (n_glb > 0) amr_cpat_off(2) = amr_region_lo_all(2, amr_cur) - amr_cpat_mar
        if (p_glb > 0) amr_cpat_off(3) = amr_region_lo_all(3, amr_cur) - amr_cpat_mar
        v1hi = (amr_region_hi_all(1, amr_cur) - amr_region_lo_all(1, amr_cur)) + 2*amr_cpat_mar
        v2hi = 0; v3hi = 0
        if (n_glb > 0) v2hi = (amr_region_hi_all(2, amr_cur) - amr_region_lo_all(2, amr_cur)) + 2*amr_cpat_mar
        if (p_glb > 0) v3hi = (amr_region_hi_all(3, amr_cur) - amr_region_lo_all(3, amr_cur)) + 2*amr_cpat_mar
        plo = amr_cpat_off
        phi(1) = amr_cpat_off(1) + v1hi; phi(2) = amr_cpat_off(2) + v2hi; phi(3) = amr_cpat_off(3) + v3hi

        if (amr_block_owner(amr_cur) /= proc_rank) return  ! contributor sends were the send phase's job

        o1 = start_idx(1); o2 = 0; o3 = 0
        if (n_glb > 0) o2 = start_idx(2)
        if (p_glb > 0) o3 = start_idx(3)
        call s_amr_rank_coarse_range(proc_rank, crlo, crhi)
        call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
        call s_phase_tic(PH_RBOWN)
        call s_amr_unpack_patch(q_coarse, bl, bh, o1, o2, o3)
        call s_phase_toc(PH_RBOWN)
#ifdef MFC_MPI
        if (nr > 0) then
            call s_phase_tic(PH_RBWAIT)
            call MPI_WAITALL(nr, amr_gcr_req(r0:r0 + nr - 1), MPI_STATUSES_IGNORE, ierr)
            call s_phase_toc(PH_RBWAIT)
            call s_phase_tic(PH_RBUNPK)
            do idx = 1, nr
                ! plan order = posting order; recompute each contributor's slice box exactly as the plan builder did
                call s_amr_rank_coarse_range(amr_gpl_src(idx, amr_cur), crlo, crhi)
                call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
                off = amr_gcr_off(r0 + idx - 1)
                if (XA_NH > 0) call s_xa_hdr_check(amr_gcr_pool(off + 1:off + XA_NH), XA_F1_SND, amr_cur, bl, bh)
                r = XA_NH
                do i = 1, sys_size
                    do g3 = bl(3), bh(3)
                        do g2 = bl(2), bh(2)
                            do g1 = bl(1), bh(1)
                                r = r + 1
                                amr_cg(i)%sf(g1 - amr_cpat_off(1), g2 - amr_cpat_off(2), &
                                       & g3 - amr_cpat_off(3)) = real(amr_gcr_pool(off + r), stp)
                            end do
                        end do
                    end do
                end do
            end do
            call s_phase_toc(PH_RBUNPK)
        end if
#endif
        call s_phase_tic(PH_RBUPD)
        do i = 1, sys_size
            $:GPU_UPDATE(device='[amr_cg(i)%sf]')
        end do
        call s_phase_toc(PH_RBUPD)

    end subroutine s_amr_gather_consume_box

    !> [amr-cov]: accumulate the step-fill coverage split for the current block (owner calls, once per fill). The whole patch is
    !! shipped/copied; the ghost-fill kernel reads only the margin plus ONE interior cell per face (floor(f/rr) +- 1), so the
    !! interior core - the patch frame's region shrunk by 1 per face - is provably dead. Level>=2 blocks fill from the parent-fine
    !! foot patch; same split in that frame.
    impure subroutine s_amr_cov_note_fill()

        integer :: r1, r2, r3, pblk, plo(3), phi(3)

        if (amr_block_level(amr_cur) >= 2) then
            pblk = f_amr_parent_block(amr_cur)
            call s_amr_parent_foot(amr_cur, pblk, plo, phi)
            r1 = phi(1) - plo(1) + 1
            r2 = 1; r3 = 1
            if (n_glb > 0) r2 = phi(2) - plo(2) + 1
            if (p_glb > 0) r3 = phi(3) - plo(3) + 1
        else
            r1 = amr_region_hi_all(1, amr_cur) - amr_region_lo_all(1, amr_cur) + 1
            r2 = 1; r3 = 1
            if (n_glb > 0) r2 = amr_region_hi_all(2, amr_cur) - amr_region_lo_all(2, amr_cur) + 1
            if (p_glb > 0) r3 = amr_region_hi_all(3, amr_cur) - amr_region_lo_all(3, amr_cur) + 1
        end if
        amr_cov_tot(1) = amr_cov_tot(1) + int(sys_size, 8)*int(r1 + 2*amr_cpat_mar, 8)*int(merge(r2 + 2*amr_cpat_mar, 1, &
                    & n_glb > 0), 8)*int(merge(r3 + 2*amr_cpat_mar, 1, p_glb > 0), 8)
        amr_cov_dead(1) = amr_cov_dead(1) + int(sys_size, 8)*int(max(r1 - 2, 0), 8)*int(merge(max(r2 - 2, 0), 1, n_glb > 0), &
                     & 8)*int(merge(max(r3 - 2, 0), 1, p_glb > 0), 8)

    end subroutine s_amr_cov_note_fill

    !> [amr-cov]: accumulate the rebuild-gather coverage split for the current box (owner calls, once per owned box). Level-1: patch
    !! words vs words prolonged into cells the same-level carry-forward overwrites (old same-level regions are disjoint, so the
    !! intersections sum exactly; each shrunk by 1 per face for the minmod stencil - conservative). Level>=2: patch words only - no
    !! carry-forward exists, the patch is live by construction.
    impure subroutine s_amr_cov_note(old_np, old_ilo, old_ext, old_level)

        integer, intent(in) :: old_np, old_ilo(:,:), old_ext(:,:), old_level(:)
        integer             :: kk, pblk, plo(3), phi(3), olo(3), ohi(3), bl(3), bh(3), r1, r2, r3
        integer(8)          :: words

        if (amr_block_level(amr_cur) >= 2) then
            pblk = f_amr_parent_block(amr_cur)
            call s_amr_parent_foot(amr_cur, pblk, plo, phi)
            r1 = phi(1) - plo(1) + 1
            r2 = 1; r3 = 1
            if (n_glb > 0) r2 = phi(2) - plo(2) + 1
            if (p_glb > 0) r3 = phi(3) - plo(3) + 1
            amr_cov_tot(3) = amr_cov_tot(3) + int(sys_size, 8)*int(r1 + 2*amr_cpat_mar, 8)*int(merge(r2 + 2*amr_cpat_mar, 1, &
                        & n_glb > 0), 8)*int(merge(r3 + 2*amr_cpat_mar, 1, p_glb > 0), 8)
            return
        end if
        r1 = amr_region_hi_all(1, amr_cur) - amr_region_lo_all(1, amr_cur) + 1
        r2 = 1; r3 = 1
        if (n_glb > 0) r2 = amr_region_hi_all(2, amr_cur) - amr_region_lo_all(2, amr_cur) + 1
        if (p_glb > 0) r3 = amr_region_hi_all(3, amr_cur) - amr_region_lo_all(3, amr_cur) + 1
        amr_cov_tot(2) = amr_cov_tot(2) + int(sys_size, 8)*int(r1 + 2*amr_cpat_mar, 8)*int(merge(r2 + 2*amr_cpat_mar, 1, &
                    & n_glb > 0), 8)*int(merge(r3 + 2*amr_cpat_mar, 1, p_glb > 0), 8)
        do kk = 1, old_np
            if (old_level(kk) /= amr_block_level(amr_cur)) cycle
            olo = old_ilo(:,kk)
            ohi = olo
            ohi(1) = olo(1) + (old_ext(1, kk) + 1)/amr_ref_ratio**old_level(kk) - 1
            if (n_glb > 0) ohi(2) = olo(2) + (old_ext(2, kk) + 1)/amr_ref_ratio**old_level(kk) - 1
            if (p_glb > 0) ohi(3) = olo(3) + (old_ext(3, kk) + 1)/amr_ref_ratio**old_level(kk) - 1
            call s_amr_box_isect(amr_region_lo_all(:,amr_cur), amr_region_hi_all(:,amr_cur), olo, ohi, bl, bh)
            if (bl(1) > bh(1) .or. bl(2) > bh(2) .or. bl(3) > bh(3)) cycle
            words = int(max(bh(1) - bl(1) + 1 - 2, 0), 8)
            if (n_glb > 0) words = words*int(max(bh(2) - bl(2) + 1 - 2, 0), 8)
            if (p_glb > 0) words = words*int(max(bh(3) - bl(3) + 1 - 2, 0), 8)
            amr_cov_dead(2) = amr_cov_dead(2) + int(sys_size, 8)*words
        end do

    end subroutine s_amr_cov_note

    !> [amr-cov] report: SUM-allreduce the counters and print once on rank 0. Collective - the caller (s_finalize_amr_module) runs
    !! it before the amr early-return so every rank participates (all-zero when amr is off).
    impure subroutine s_amr_cov_report()

        integer(8)                  :: tot(3), dead(3), cad(2), cadr(2)
        integer                     :: ierr, f
        character(len=8), parameter :: nm(3) = ["stepfill", "rb-L1   ", "rb-L2   "]

        tot = amr_cov_tot; dead = amr_cov_dead
        cad(1) = amr_cad_tot; cad(2) = amr_cad_esc; cadr = cad
#ifdef MFC_MPI
        call MPI_ALLREDUCE(amr_cov_tot, tot, 3, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(amr_cov_dead, dead, 3, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(cad, cadr, 2, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
        if (proc_rank == 0) then
            do f = 1, 3
                if (tot(f) > 0) write (0, '(A,A,A,I0,A,I0,A,F6.3)') ' [amr-cov] ', nm(f), ' words ', tot(f), ' dead ', dead(f), &
                    & ' frac ', real(dead(f))/real(tot(f))
            end do
            ! cadence containment: escaped > 0 means a feature outran amr_buf between regrids (see the decl)
            if (cadr(1) > 0) write (0, '(A,I0,A,I0,A,F6.3)') ' [amr-cad] L1 tags ', cadr(1), ' escaped ', cadr(2), ' frac ', &
                & real(cadr(2))/real(cadr(1))
        end if

    end subroutine s_amr_cov_report

    impure subroutine s_amr_gather_coarse_patch(q_coarse, pull_host)

        type(scalar_field), dimension(sys_size), intent(in) :: q_coarse
        !> runtime callers pass .true. (coarse device-current); init/regrid pass .false. (host is truth)
        logical, intent(in)   :: pull_host
        integer               :: i, g1, g2, g3, o1, o2, o3, owner, r, idx, boxsz, maxsz, nsrc, ierr
        integer               :: v1hi, v2hi, v3hi, plo(3), phi(3), crlo(3), crhi(3), bl(3), bh(3)
        real(wp), allocatable :: rbuf(:,:), sbuf(:)
        integer, allocatable  :: reqs(:), srank(:)

        ! multi-level: a level>=2 block's coarse side is its PARENT block's fine cells, not the L0 base grid q_coarse - gather
        ! amr_cg
        ! from the parent's fine array in the parent-fine frame (isect already parent-fine from s_set_amr_fine_geometry). np=1 is a
        ! local copy; the np>=2 P2P version (parent owner -> block owner, mirroring the L0 path) is future work.

        if (amr_block_level(amr_cur) >= 2) then
            if (amr_rg_gather) call s_phase_tic(PH_PGALL)
            call s_amr_gather_from_parent(pull_host)
            if (amr_rg_gather) call s_phase_toc(PH_PGALL)
            return
        end if

        ! block-local patch frame (cell 0 == global region_lo-nmar; collapsed dims -> 0) + its GLOBAL cell range [plo:phi]
        amr_cpat_off = 0
        amr_cpat_off(1) = amr_region_lo_all(1, amr_cur) - amr_cpat_mar
        if (n_glb > 0) amr_cpat_off(2) = amr_region_lo_all(2, amr_cur) - amr_cpat_mar
        if (p_glb > 0) amr_cpat_off(3) = amr_region_lo_all(3, amr_cur) - amr_cpat_mar
        v1hi = (amr_region_hi_all(1, amr_cur) - amr_region_lo_all(1, amr_cur)) + 2*amr_cpat_mar
        v2hi = 0; v3hi = 0
        if (n_glb > 0) v2hi = (amr_region_hi_all(2, amr_cur) - amr_region_lo_all(2, amr_cur)) + 2*amr_cpat_mar
        if (p_glb > 0) v3hi = (amr_region_hi_all(3, amr_cur) - amr_region_lo_all(3, amr_cur)) + 2*amr_cpat_mar
        plo = amr_cpat_off
        phi(1) = amr_cpat_off(1) + v1hi; phi(2) = amr_cpat_off(2) + v2hi; phi(3) = amr_cpat_off(3) + v3hi

        owner = amr_block_owner(amr_cur)
        o1 = start_idx(1); o2 = 0; o3 = 0
        if (n_glb > 0) o2 = start_idx(2)
        if (p_glb > 0) o3 = start_idx(3)
        maxsz = sys_size*(v1hi + 1)*(v2hi + 1)*(v3hi + 1)

        ! np=1: the sole owner holds every covered coarse cell, so copy q_coarse->amr_cg on-device (same index map as
        ! s_amr_unpack_patch), skipping the device->host->device round-trip. Only for pull_host; init/regrid (.not. pull_host) falls
        ! through to the host path (device copy may be stale).
        if (num_procs == 1 .and. pull_host) then
            call s_amr_rank_coarse_range(owner, crlo, crhi)
            call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
            call s_amr_gather_own_box_device(q_coarse, bl, bh, o1, o2, o3)  ! same kernel the np>1 owner path uses
            return
        end if

        ! np>1 runtime (pull_host): NO full-field host pull - the owner's own-box copy, the non-owner pack, and the received-box
        ! unpacks all run on the DEVICE over only the overlap boxes, so just the contiguous wire buffers cross PCIe (MPI stays on
        ! host buffers). Init/regrid (.not. pull_host): host is truth, so the host pack/unpack paths below read it directly.

        ! block set changed: rebuild the cached overlap-rank lists (same lazy trigger as s_amr_fine_fine_halo; local, replicated)
        if (amr_rg_gather) call s_phase_tic(PH_RBSEAM)
        if (amr_seam_pairs_dirty .or. amr_seam_pairs_nblk /= amr_num_blocks) call s_amr_build_seam_pairs()
        if (amr_rg_gather) call s_phase_toc(PH_RBSEAM)

        if (proc_rank == owner) then
            ! fill the cells this rank holds locally (own box), then receive the rest from the other coarse-owners
            call s_amr_rank_coarse_range(proc_rank, crlo, crhi)
            call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
            if (pull_host) then
                ! runtime: q_coarse is device-current - copy the own box on the device (same index map/assignment as the host path)
                call s_amr_gather_own_box_device(q_coarse, bl, bh, o1, o2, o3)
            else
                if (amr_rg_gather) call s_phase_tic(PH_RBOWN)
                call s_amr_unpack_patch(q_coarse, bl, bh, o1, o2, o3)  ! local read: q_coarse own frame -> amr_cg patch frame
                if (amr_rg_gather) call s_phase_toc(PH_RBOWN)
            end if
            ! count + post recvs from every OTHER rank whose owned range overlaps the patch (cached list; every listed rank
            ! overlaps by construction)
            if (amr_rg_gather) call s_phase_tic(PH_RBPOST)
            nsrc = 0
            do idx = 1, amr_ovl_gather_n(amr_cur)
                if (amr_ovl_gather(idx, amr_cur) /= owner) nsrc = nsrc + 1
            end do
            if (amr_rg_gather .and. amr_gpl_valid) then
                @:ASSERT(nsrc == amr_gpl_nsrc(amr_cur), "gather plan: contributor count mismatch")
            end if
            if (amr_rg_gather) call s_phase_toc(PH_RBPOST)
            if (nsrc > 0) then
                if (amr_rg_gather) call s_phase_tic(PH_RBALLOC)
                allocate (rbuf(maxsz + XA_NH, nsrc), reqs(nsrc), srank(nsrc))
                if (amr_rg_gather) call s_phase_toc(PH_RBALLOC)
                if (amr_rg_gather) call s_phase_tic(PH_RBPOST)
                nsrc = 0
                do idx = 1, amr_ovl_gather_n(amr_cur)
                    r = amr_ovl_gather(idx, amr_cur)
                    if (r == owner) cycle
                    call s_amr_rank_coarse_range(r, crlo, crhi)
                    call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
                    boxsz = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                    nsrc = nsrc + 1; srank(nsrc) = r
                    if (amr_rg_gather .and. amr_gpl_valid) then
                        @:ASSERT(amr_gpl_src(nsrc, amr_cur) == r .and. amr_gpl_sz(nsrc, amr_cur) == boxsz, &
                                 & "gather plan: recv entry mismatch")
                    end if
#ifdef MFC_MPI
                    call s_xa_rec(XA_F1_RCV, 2, boxsz, amr_cur)
                    call MPI_IRECV(rbuf(1, nsrc), boxsz + XA_NH, mpi_p, r, amr_cur, MPI_COMM_WORLD, reqs(nsrc), ierr)
#endif
                end do
                if (amr_rg_gather) call s_phase_toc(PH_RBPOST)
#ifdef MFC_MPI
                if (amr_rg_gather) call s_phase_tic(PH_RBWAIT)
                call MPI_WAITALL(nsrc, reqs, MPI_STATUSES_IGNORE, ierr)
                if (amr_rg_gather) call s_phase_toc(PH_RBWAIT)
#endif
                if (amr_rg_gather) call s_phase_tic(PH_RBUNPK)
                do idx = 1, nsrc
                    call s_amr_rank_coarse_range(srank(idx), crlo, crhi)
                    call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
                    if (XA_NH > 0) call s_xa_hdr_check(rbuf(:,idx), XA_F1_SND, amr_cur, bl, bh)
                    if (pull_host) then
                        ! runtime: unpack ONLY this box's wire buffer on the device (same order/cast as the host unpack below)
                        boxsz = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                        call s_amr_unpack_box_device(bl, bh, rbuf(XA_NH + 1:XA_NH + boxsz,idx))
                        cycle
                    end if
                    ! unpack in the SAME (i, g3, g2, g1) order the sender packed; place at amr_cg patch-local index
                    r = XA_NH
                    do i = 1, sys_size
                        do g3 = bl(3), bh(3)
                            do g2 = bl(2), bh(2)
                                do g1 = bl(1), bh(1)
                                    r = r + 1
                                    amr_cg(i)%sf(g1 - amr_cpat_off(1), g2 - amr_cpat_off(2), g3 - amr_cpat_off(3)) = real(rbuf(r, &
                                           & idx), stp)
                                end do
                            end do
                        end do
                    end do
                end do
                if (amr_rg_gather) call s_phase_toc(PH_RBUNPK)
                if (amr_rg_gather) call s_phase_tic(PH_RBALLOC)
                deallocate (rbuf, reqs, srank)
                if (amr_rg_gather) call s_phase_toc(PH_RBALLOC)
            end if
            ! host path only: the runtime device path wrote amr_cg on the device directly (host amr_cg stays stale, as at np=1 -
            ! runtime consumers read the device copy)
            if (.not. pull_host) then
                if (amr_rg_gather) call s_phase_tic(PH_RBUPD)
                do i = 1, sys_size
                    $:GPU_UPDATE(device='[amr_cg(i)%sf]')
                end do
                if (amr_rg_gather) call s_phase_toc(PH_RBUPD)
            end if
        else
            ! non-owner: if my owned coarse range overlaps the patch, pack my slice (wp) and send it to the owner
            call s_amr_rank_coarse_range(proc_rank, crlo, crhi)
            call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
            if (bl(1) <= bh(1) .and. bl(2) <= bh(2) .and. bl(3) <= bh(3)) then
                boxsz = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                if (amr_rg_gather .and. amr_gpl_valid) then
                    ! this rank must appear in the plan's contributor list for this box, with this exact size
                    nsrc = 0
                    do idx = 1, amr_gpl_nsrc(amr_cur)
                        if (amr_gpl_src(idx, amr_cur) == proc_rank .and. amr_gpl_sz(idx, amr_cur) == boxsz) nsrc = 1
                    end do
                    @:ASSERT(nsrc == 1, "gather plan: send entry mismatch")
                end if
                if (amr_rg_gather) call s_phase_tic(PH_RBRSV)
                call s_amr_gsnd_reserve(maxsz + XA_NH)
                if (amr_rg_gather) call s_phase_toc(PH_RBRSV)
                amr_gsnd_n = amr_gsnd_n + 1
                if (pull_host) then
                    ! runtime: pack the overlap box on the device straight into the pool slot (only the box crosses PCIe);
                    ! the slice leaves the I1b header words ahead of the data (kernel untouched)
                    call s_amr_pack_box_device(q_coarse, bl, bh, o1, o2, o3, amr_gsnd_pool(XA_NH + 1:,amr_gsnd_n))
                else
                    if (amr_rg_gather) call s_phase_tic(PH_RBPACK)
                    idx = XA_NH
                    do i = 1, sys_size
                        do g3 = bl(3), bh(3)
                            do g2 = bl(2), bh(2)
                                do g1 = bl(1), bh(1)
                                    idx = idx + 1
                                    amr_gsnd_pool(idx, amr_gsnd_n) = real(q_coarse(i)%sf(g1 - o1, g2 - o2, g3 - o3), wp)
                                end do
                            end do
                        end do
                    end do
                    if (amr_rg_gather) call s_phase_toc(PH_RBPACK)
                end if
#ifdef MFC_MPI
                ! NON-BLOCKING: the owner's per-box IRECV/WAITALL still orders the data correctly, but this rank no longer
                ! rendezvouses on every box. Completed by s_amr_gather_send_flush (caller) or the drain in s_amr_gsnd_reserve.
                if (amr_rg_gather) call s_phase_tic(PH_RBSEND)
                if (XA_NH > 0) call s_xa_hdr_pack(amr_gsnd_pool(:,amr_gsnd_n), XA_F1_SND, amr_cur, bl, bh)
                call s_xa_rec(XA_F1_SND, 1, boxsz, amr_cur)
                call MPI_ISEND(amr_gsnd_pool(1, amr_gsnd_n), boxsz + XA_NH, mpi_p, owner, amr_cur, MPI_COMM_WORLD, &
                               & amr_gsnd_req(amr_gsnd_n), ierr)
                if (amr_rg_gather) call s_phase_toc(PH_RBSEND)
#endif
            end if
        end if

    end subroutine s_amr_gather_coarse_patch

    !> Non-polytropic QBMM analogue of s_amr_gather_coarse_patch: gather the current block's coarse pb/mv patch into amr_cg_pb/mv
    !! (patch frame, cell 0 == amr_cpat_off), P2P from the coarse-cell owners into the block owner. Per-cell payload = 2*nnode*nb
    !! (pb block then mv block). Single-level only (level>=2 QBMM np>=2 is checker-gated); wire is wp, cast to stp on unpack
    !! (identity for stp coarse), so at np=1 the owner copies its own coarse over the patch bit-for-bit. TWIN
    !! s_amr_gather_coarse_patch (pb/mv<->q): mirrors the q_cons gather's P2P skeleton and patch-local frame - keep lockstep.
    impure subroutine s_amr_gather_coarse_patch_pbmv(pb_coarse, mv_coarse, pull_host)

        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(in) :: pb_coarse, mv_coarse
        !> runtime callers pass .true. (coarse device-current); init/regrid pass .false. (host is truth)
        logical, intent(in)   :: pull_host
        integer               :: q, ib_, g1, g2, g3, o1, o2, o3, owner, r, idx, boxsz, maxsz, nsrc, ierr
        integer               :: v1hi, v2hi, v3hi, plo(3), phi(3), crlo(3), crhi(3), bl(3), bh(3), cellsz
        real(wp), allocatable :: rbuf(:,:), sbuf(:)
        integer, allocatable  :: reqs(:), srank(:)

        ! single-level only: a level>=2 block's coarse side is its parent's fine pb/mv, distributed only at np=1 - the checker gate
        ! keeps multi-level QBMM np>=2 fail-closed, so this must never be reached at level>=2.

        if (amr_block_level(amr_cur) >= 2) return

        cellsz = 2*nnode*nb

        ! block-local patch frame (cell 0 == global region_lo-nmar; collapsed dims -> 0) + its GLOBAL cell range [plo:phi]
        amr_cpat_off = 0
        amr_cpat_off(1) = amr_region_lo_all(1, amr_cur) - amr_cpat_mar
        if (n_glb > 0) amr_cpat_off(2) = amr_region_lo_all(2, amr_cur) - amr_cpat_mar
        if (p_glb > 0) amr_cpat_off(3) = amr_region_lo_all(3, amr_cur) - amr_cpat_mar
        v1hi = (amr_region_hi_all(1, amr_cur) - amr_region_lo_all(1, amr_cur)) + 2*amr_cpat_mar
        v2hi = 0; v3hi = 0
        if (n_glb > 0) v2hi = (amr_region_hi_all(2, amr_cur) - amr_region_lo_all(2, amr_cur)) + 2*amr_cpat_mar
        if (p_glb > 0) v3hi = (amr_region_hi_all(3, amr_cur) - amr_region_lo_all(3, amr_cur)) + 2*amr_cpat_mar
        plo = amr_cpat_off
        phi(1) = amr_cpat_off(1) + v1hi; phi(2) = amr_cpat_off(2) + v2hi; phi(3) = amr_cpat_off(3) + v3hi

        owner = amr_block_owner(amr_cur)
        o1 = start_idx(1); o2 = 0; o3 = 0
        if (n_glb > 0) o2 = start_idx(2)
        if (p_glb > 0) o3 = start_idx(3)
        maxsz = cellsz*(v1hi + 1)*(v2hi + 1)*(v3hi + 1)

        ! np=1: the sole owner holds every covered coarse cell, so copy pb_coarse/mv_coarse->amr_cg_pb/mv on-device over the
        ! in-domain patch. Only for pull_host; init/regrid (.not. pull_host) falls through to the host path (device copy may be
        ! stale).
        if (num_procs == 1 .and. pull_host) then
            call s_amr_rank_coarse_range(owner, crlo, crhi)
            call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
            call s_amr_gather_own_box_pbmv_device(pb_coarse, mv_coarse, bl, bh, o1, o2, o3)  ! same kernel the np>1 owner path uses
            return
        end if

        ! np>1 runtime (pull_host): NO full-field host pull - the owner's own-box copy, the non-owner pack, and the received-box
        ! unpacks all run on the DEVICE over only the overlap boxes (mirror of s_amr_gather_coarse_patch). Init/regrid
        ! (.not. pull_host): host is truth, so the host pack/unpack paths below read it directly.

        ! block set changed: rebuild the cached overlap-rank lists (same lazy trigger as s_amr_fine_fine_halo; local, replicated)
        if (amr_seam_pairs_dirty .or. amr_seam_pairs_nblk /= amr_num_blocks) call s_amr_build_seam_pairs()

        if (proc_rank == owner) then
            ! fill the cells this rank holds locally (own box), then receive the rest from the other coarse-owners
            call s_amr_rank_coarse_range(proc_rank, crlo, crhi)
            call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
            if (pull_host) then
                ! runtime: pb/mv device-current - copy the own box on the device (same index map/assignment as the host path)
                call s_amr_gather_own_box_pbmv_device(pb_coarse, mv_coarse, bl, bh, o1, o2, o3)
            else
                do ib_ = 1, nb
                    do q = 1, nnode
                        do g3 = bl(3), bh(3)
                            do g2 = bl(2), bh(2)
                                do g1 = bl(1), bh(1)
                                    amr_cg_pb(g1 - amr_cpat_off(1), g2 - amr_cpat_off(2), g3 - amr_cpat_off(3), q, &
                                              & ib_) = pb_coarse(g1 - o1, g2 - o2, g3 - o3, q, ib_)
                                    amr_cg_mv(g1 - amr_cpat_off(1), g2 - amr_cpat_off(2), g3 - amr_cpat_off(3), q, &
                                              & ib_) = mv_coarse(g1 - o1, g2 - o2, g3 - o3, q, ib_)
                                end do
                            end do
                        end do
                    end do
                end do
            end if
            ! count + post recvs from every OTHER rank whose owned range overlaps the patch (cached list; every listed rank
            ! overlaps by construction)
            nsrc = 0
            do idx = 1, amr_ovl_gather_n(amr_cur)
                if (amr_ovl_gather(idx, amr_cur) /= owner) nsrc = nsrc + 1
            end do
            if (nsrc > 0) then
                allocate (rbuf(maxsz + XA_NH, nsrc), reqs(nsrc), srank(nsrc))
                nsrc = 0
                do idx = 1, amr_ovl_gather_n(amr_cur)
                    r = amr_ovl_gather(idx, amr_cur)
                    if (r == owner) cycle
                    call s_amr_rank_coarse_range(r, crlo, crhi)
                    call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
                    boxsz = cellsz*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                    nsrc = nsrc + 1; srank(nsrc) = r
#ifdef MFC_MPI
                    call s_xa_rec(XA_F3_RCV, 2, boxsz, amr_cur)
                    call MPI_IRECV(rbuf(1, nsrc), boxsz + XA_NH, mpi_p, r, amr_cur, MPI_COMM_WORLD, reqs(nsrc), ierr)
#endif
                end do
#ifdef MFC_MPI
                call MPI_WAITALL(nsrc, reqs, MPI_STATUSES_IGNORE, ierr)
#endif
                do idx = 1, nsrc
                    call s_amr_rank_coarse_range(srank(idx), crlo, crhi)
                    call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
                    if (XA_NH > 0) call s_xa_hdr_check(rbuf(:,idx), XA_F3_SND, amr_cur, bl, bh)
                    if (pull_host) then
                        ! runtime: unpack ONLY this box's wire buffer on the device (same order/cast as the host unpack below)
                        boxsz = cellsz*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                        call s_amr_unpack_box_pbmv_device(bl, bh, rbuf(XA_NH + 1:XA_NH + boxsz,idx))
                        cycle
                    end if
                    ! unpack in the SAME (ib_, q, g3, g2, g1) order the sender packed - pb block then mv block
                    r = XA_NH
                    do ib_ = 1, nb
                        do q = 1, nnode
                            do g3 = bl(3), bh(3)
                                do g2 = bl(2), bh(2)
                                    do g1 = bl(1), bh(1)
                                        r = r + 1
                                        amr_cg_pb(g1 - amr_cpat_off(1), g2 - amr_cpat_off(2), g3 - amr_cpat_off(3), q, &
                                                  & ib_) = real(rbuf(r, idx), stp)
                                    end do
                                end do
                            end do
                        end do
                    end do
                    do ib_ = 1, nb
                        do q = 1, nnode
                            do g3 = bl(3), bh(3)
                                do g2 = bl(2), bh(2)
                                    do g1 = bl(1), bh(1)
                                        r = r + 1
                                        amr_cg_mv(g1 - amr_cpat_off(1), g2 - amr_cpat_off(2), g3 - amr_cpat_off(3), q, &
                                                  & ib_) = real(rbuf(r, idx), stp)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
                deallocate (rbuf, reqs, srank)
            end if
            ! host path only: the runtime device path wrote amr_cg_pb/mv on the device directly (host copies stay stale, as at np=1
            ! -
            ! runtime consumers read the device copy)
            if (.not. pull_host) then
                $:GPU_UPDATE(device='[amr_cg_pb, amr_cg_mv]')
            end if
        else
            ! non-owner: if my owned coarse range overlaps the patch, pack my slice (wp) and send it to the owner
            call s_amr_rank_coarse_range(proc_rank, crlo, crhi)
            call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
            if (bl(1) <= bh(1) .and. bl(2) <= bh(2) .and. bl(3) <= bh(3)) then
                boxsz = cellsz*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                allocate (sbuf(boxsz + XA_NH))
                if (pull_host) then
                    ! runtime: pack the overlap box on the device straight into sbuf (only the box crosses PCIe);
                    ! the slice leaves the I1b header words ahead of the data (kernel untouched)
                    call s_amr_pack_box_pbmv_device(pb_coarse, mv_coarse, bl, bh, o1, o2, o3, sbuf(XA_NH + 1:))
                else
                    idx = XA_NH
                    do ib_ = 1, nb
                        do q = 1, nnode
                            do g3 = bl(3), bh(3)
                                do g2 = bl(2), bh(2)
                                    do g1 = bl(1), bh(1)
                                        idx = idx + 1; sbuf(idx) = real(pb_coarse(g1 - o1, g2 - o2, g3 - o3, q, ib_), wp)
                                    end do
                                end do
                            end do
                        end do
                    end do
                    do ib_ = 1, nb
                        do q = 1, nnode
                            do g3 = bl(3), bh(3)
                                do g2 = bl(2), bh(2)
                                    do g1 = bl(1), bh(1)
                                        idx = idx + 1; sbuf(idx) = real(mv_coarse(g1 - o1, g2 - o2, g3 - o3, q, ib_), wp)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end if
#ifdef MFC_MPI
                if (XA_NH > 0) call s_xa_hdr_pack(sbuf, XA_F3_SND, amr_cur, bl, bh)
                call s_xa_rec(XA_F3_SND, 1, boxsz, amr_cur)
                call MPI_SEND(sbuf, boxsz + XA_NH, mpi_p, owner, amr_cur, MPI_COMM_WORLD, ierr)
#endif
                deallocate (sbuf)
            end if
        end if

        ! host-consumer callers (init/regrid prolong) need the gathered patch on the host
        if (.not. pull_host) then
            $:GPU_UPDATE(host='[amr_cg_pb, amr_cg_mv]')
        end if

    end subroutine s_amr_gather_coarse_patch_pbmv

    !> Multi-level gather: fill amr_cg (the current level>=2 block's coarse patch) from its PARENT block's fine array, in the
    !! parent-fine cell frame (amr_isect_lo/hi already parent-fine from s_set_amr_fine_geometry). np=1 = a local copy on the owner
    !! (which also owns the parent); the np>=2 P2P version (parent owner -> block owner) is future work.
    impure subroutine s_amr_gather_from_parent(pull_host)

        logical, intent(in) :: pull_host
        integer             :: pblk

        pblk = f_amr_parent_block(amr_cur)
        if (amr_rg_gather .and. amr_gpl_valid) then
            if (amr_block_owner(pblk) == amr_block_owner(amr_cur)) then
                @:ASSERT(amr_gpl_psrc(amr_cur) == -1, "gather plan: expected co-located parent")
            else
                @:ASSERT(amr_gpl_psrc(amr_cur) == amr_block_owner(pblk), "gather plan: parent source mismatch")
            end if
        end if
        ! lock-step fill: gather from the parent's CURRENT fine state. pull_host stays in the signature for the level-1 path.
        ! Owner-guard at the CALL SITE: the parent slot is allocated only on ITS owner, and passing its store slot on any
        ! other rank would dereference an unallocated slot. So both participants enter - the parent's owner to pack and send, the
        ! block's owner to receive - and every other rank stays out. When the two coincide (np=1, or a co-located tower) this is the
        ! old local-copy path unchanged. to_host = .not. pull_host: init/regrid (pull_host=F) feed the host prolong/self-test;
        ! runtime (pull_host=T) reads amr_cg on the device in the C/F ghost-fill, so skip the device->host copy.
        if (amr_block_owner(pblk) == proc_rank) then
            ! parent owner: local device copy when it also owns the block, otherwise pack and send.
            if (amr_rg_gather) call s_phase_tic(PH_PGSEND)
            call s_amr_gather_from_parent_field_cons(amr_cur, pblk, amr_loc_of(pblk), .not. pull_host)
            if (amr_rg_gather) call s_phase_toc(PH_PGSEND)
        else if (amr_rank_owns_block) then
            ! block owner only: receive. Deliberately does NOT take the parent field - amr_slots(pblk) is unallocated here.
            if (amr_rg_gather) call s_phase_tic(PH_PGRECV)
            call s_amr_recv_parent_patch(pblk, .not. pull_host)
            if (amr_rg_gather) call s_phase_toc(PH_PGRECV)
        end if

    end subroutine s_amr_gather_from_parent

    !> Gather amr_cg (the current level>=2 block's coarse patch) from a SPECIFIC parent snapshot field qp, in the parent-fine cell
    !! frame (amr_isect_lo/hi already parent-fine from s_set_amr_fine_geometry). The subcycle recursion calls this twice per parent
    !! substep - qp = the parent slot's q_cons_stor (t^n bracket) then q_cons (t^{n+1} bracket) - to build the child's two
    !! ghost-lerp sources. np=1 = a local copy on the owner (which also owns the parent); np>=2 P2P (parent owner -> block owner) is
    !! future work. Two sources, one body: the parent's conserved state lives in the flat store (`_st`, keyed by its slot), its
    !! SSP-RK stage backup q_cons_stor is still a per-slot scalar_field array (`_sf`).
    #:for GSFX, GARR in [('cons', 'amr_cons_st'), ('stor', 'amr_stor_st')]
        impure subroutine s_amr_gather_from_parent_field_${GSFX}$(cblk, pblk, qp, to_host)

            !> the child block (explicit, NOT amr_cur: the chunked send phase calls this before the consume phase's geometry, when
            !! amr_cur points at another box)
            integer, intent(in) :: cblk
            integer, intent(in) :: pblk
            integer, intent(in) :: qp       !< parent's flat-store slot
            logical, intent(in) :: to_host  !< host copy of amr_cg needed (init/regrid), not runtime
            integer             :: w1, w2, w3, powner, cowner, boxsz, ierr
            integer             :: plo(3), phi(3)

            ! Patch box in the PARENT-FINE frame. Both the child owner and the parent owner must agree on it, so derive it from
            ! REPLICATED metadata (amr_region_*_all + the global amr_ref_ratio) rather than from amr_isect_lo/hi, which is the empty
            ! footprint on a non-owner of this block. On the child owner the two agree by construction (s_set_amr_fine_geometry).

            call s_amr_parent_foot(cblk, pblk, plo, phi)
            amr_cpat_off = 0
            amr_cpat_off(1) = plo(1) - amr_cpat_mar
            if (n_glb > 0) amr_cpat_off(2) = plo(2) - amr_cpat_mar
            if (p_glb > 0) amr_cpat_off(3) = plo(3) - amr_cpat_mar
            w1 = (phi(1) - plo(1)) + 2*amr_cpat_mar
            w2 = 0; w3 = 0
            if (n_glb > 0) w2 = (phi(2) - plo(2)) + 2*amr_cpat_mar
            if (p_glb > 0) w3 = (phi(3) - plo(3)) + 2*amr_cpat_mar

            cowner = amr_block_owner(cblk); powner = amr_block_owner(pblk)
            if (powner == cowner) then
                ! co-located (always true at np=1, and under tower co-location): straight device copy, bit-for-bit as before.
                call s_amr_copy_parent_patch_${GSFX}$(qp, w1, w2, w3, to_host)
                return
            end if

#ifdef MFC_MPI
            ! Split ownership, parent side: exactly one destination (the block's owner) and one box, so a blocking pair suffices -
            ! no
            ! overlap map and no collective, matching the L0<->L1 gather's "non-participants send/recv nothing" property.
            ! NON-BLOCKING, via the same deferred pool the level-1 gather uses (see s_amr_gsnd_reserve). The old code used a
            ! BLOCKING MPI_SEND here, so the parent's owner rendezvoused with the child's owner once per box, in lockstep - the
            ! defect m_amr.fpp:256 records for the level-1 path and fixes there, never applied to this one. Level>=2 is the
            ! MAJORITY of boxes (160 of 224 at cap 64), and it measured 420 ms per send / 8.6% of wall at the production point.
            ! The pool owns the buffer because an ISEND requires it to stay live until completion; the drain is the existing
            ! s_amr_gather_send_flush after the rebuild's box loop.
            boxsz = sys_size*(w1 + 1)*(w2 + 1)*(w3 + 1)
            ! guard on the plan alone, NOT amr_rg_gather (nothing sets it since the chunked rebuild landed): this is the ONE
            ! step-1 assert on the chunked path's live route - a send packed short of the plan-sized recv completes short and
            ! the consume unpacks stale pool bytes, the silent-wrong-answer class. amr_gpl_valid is false outside the rebuild
            ! box loop, so subcycle/per-step calls never consult the plan.
            if (amr_gpl_valid) then
                @:ASSERT(amr_gpl_psz(cblk) == boxsz, "gather plan: parent send size mismatch")
            end if
            call s_amr_gsnd_reserve(boxsz + XA_NH)
            amr_gsnd_n = amr_gsnd_n + 1
            ! header written on the host AFTER the device pack lands (copyout) - data at XA_NH+1 via the slice
            call s_amr_pack_parent_patch_device_${GSFX}$(qp, w1, w2, w3, amr_gsnd_pool(XA_NH + 1:,amr_gsnd_n))
            if (XA_NH > 0) call s_xa_hdr_pack(amr_gsnd_pool(:,amr_gsnd_n), XA_F2_SND, cblk, plo, phi)
            call s_xa_rec(XA_F2_SND, 1, boxsz, cblk)
            call MPI_ISEND(amr_gsnd_pool(1, amr_gsnd_n), boxsz + XA_NH, mpi_p, cowner, cblk, MPI_COMM_WORLD, &
                           & amr_gsnd_req(amr_gsnd_n), ierr)
#endif

        end subroutine s_amr_gather_from_parent_field_${GSFX}$
    #:endfor

    !> Receive side of the split-ownership parent gather: fill amr_cg from the parent's owner. Takes only pblk - the parent slot is
    !! NOT allocated on this rank, so the parent field must not appear in the signature. Recomputes the patch box from the same
    !! replicated metadata the sender uses, so the two agree without a handshake.
    impure subroutine s_amr_recv_parent_patch(pblk, to_host)

        integer, intent(in)   :: pblk
        logical, intent(in)   :: to_host
        integer               :: w1, w2, w3, powner, boxsz, ierr, plo(3), phi(3)
        real(wp), allocatable :: xbuf(:)

        call s_amr_parent_foot(amr_cur, pblk, plo, phi)
        amr_cpat_off = 0
        amr_cpat_off(1) = plo(1) - amr_cpat_mar
        if (n_glb > 0) amr_cpat_off(2) = plo(2) - amr_cpat_mar
        if (p_glb > 0) amr_cpat_off(3) = plo(3) - amr_cpat_mar
        w1 = (phi(1) - plo(1)) + 2*amr_cpat_mar
        w2 = 0; w3 = 0
        if (n_glb > 0) w2 = (phi(2) - plo(2)) + 2*amr_cpat_mar
        if (p_glb > 0) w3 = (phi(3) - plo(3)) + 2*amr_cpat_mar

#ifdef MFC_MPI
        powner = amr_block_owner(pblk)
        boxsz = sys_size*(w1 + 1)*(w2 + 1)*(w3 + 1)
        if (amr_rg_gather .and. amr_gpl_valid) then
            @:ASSERT(amr_gpl_psrc(amr_cur) == powner .and. amr_gpl_psz(amr_cur) == boxsz, "gather plan: parent recv mismatch")
        end if
        allocate (xbuf(boxsz + XA_NH))
        call s_xa_rec(XA_F2_RCV, 2, boxsz, amr_cur)
        call MPI_RECV(xbuf, boxsz + XA_NH, mpi_p, powner, amr_cur, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        if (XA_NH > 0) call s_xa_hdr_check(xbuf, XA_F2_SND, amr_cur, plo, phi)
        call s_amr_unpack_parent_patch_device(w1, w2, w3, xbuf(XA_NH + 1:XA_NH + boxsz), to_host)
        deallocate (xbuf)
#endif

    end subroutine s_amr_recv_parent_patch

    !> DEVICE pack of the parent's fine patch into a flat buffer. Same index map as s_amr_copy_parent_patch, writing the send buffer
    !! instead of amr_cg, so the two sides of the P2P gather cannot drift apart.
    #:for GSFX, GARR in [('cons', 'amr_cons_st'), ('stor', 'amr_stor_st')]
        #:set QP = lambda ix: GARR + '(g1 + o1, g2 + o2, g3 + o3, ' + ix + ', qp)'
        impure subroutine s_amr_pack_parent_patch_device_${GSFX}$(qp, w1, w2, w3, buf)

            integer, intent(in)                 :: qp  !< parent's flat-store slot
            integer, intent(in)                 :: w1, w2, w3
            real(wp), intent(inout), contiguous :: buf(:)
            integer                             :: i, g1, g2, g3, o1, o2, o3, n1, n2, n3

            o1 = amr_cpat_off(1); o2 = amr_cpat_off(2); o3 = amr_cpat_off(3)
            n1 = w1 + 1; n2 = w2 + 1; n3 = w3 + 1
            $:GPU_PARALLEL_LOOP(collapse=4, copyout='[buf]')
            do i = 1, sys_size
                do g3 = 0, w3
                    do g2 = 0, w2
                        do g1 = 0, w1
                            buf(1 + g1 + n1*(g2 + n2*(g3 + n3*(i - 1)))) = real(${QP('i')}$, wp)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

        end subroutine s_amr_pack_parent_patch_device_${GSFX}$
    #:endfor

    !> DEVICE unpack of a received parent patch into amr_cg. Inverse of s_amr_pack_parent_patch_device; to_host mirrors
    !! s_amr_copy_parent_patch (init/regrid host consumers need the host copy, runtime reads amr_cg on the device).
    impure subroutine s_amr_unpack_parent_patch_device(w1, w2, w3, buf, to_host)

        integer, intent(in)              :: w1, w2, w3
        real(wp), intent(in), contiguous :: buf(:)
        logical, intent(in)              :: to_host
        integer                          :: i, g1, g2, g3, n1, n2, n3

        n1 = w1 + 1; n2 = w2 + 1; n3 = w3 + 1
        $:GPU_PARALLEL_LOOP(collapse=4, copyin='[buf]')
        do i = 1, sys_size
            do g3 = 0, w3
                do g2 = 0, w2
                    do g1 = 0, w1
                        amr_cg(i)%sf(g1, g2, g3) = buf(1 + g1 + n1*(g2 + n2*(g3 + n3*(i - 1))))
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        if (to_host) then
            do i = 1, sys_size
                $:GPU_UPDATE(host='[amr_cg(i)%sf]')
            end do
        end if

    end subroutine s_amr_unpack_parent_patch_device

    !> Device kernel for s_amr_gather_from_parent: copy the parent block's fine patch into amr_cg over [amr_cpat_off : + w]. amr_cg
    !! is then synced to host for host consumers (init self-test's restrict-prolong check). Two sources, one body - see
    !! s_amr_gather_from_parent_field_st/_sf.
    #:for GSFX, GARR in [('cons', 'amr_cons_st'), ('stor', 'amr_stor_st')]
        #:set QP = lambda ix: GARR + '(g1 + o1, g2 + o2, g3 + o3, ' + ix + ', qp)'
        impure subroutine s_amr_copy_parent_patch_${GSFX}$(qp, w1, w2, w3, to_host)

            integer, intent(in) :: qp  !< parent's flat-store slot
            integer, intent(in) :: w1, w2, w3
            !> .true. only for the init/regrid HOST consumers (whole-block host prolong + restrict-prolong self-test). The runtime
            !! C/F ghost-fill reads amr_cg on the DEVICE (filled by the kernel below), so no device->host copy is needed.
            logical, intent(in) :: to_host
            integer             :: i, g1, g2, g3, o1, o2, o3

            o1 = amr_cpat_off(1); o2 = amr_cpat_off(2); o3 = amr_cpat_off(3)
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do g3 = 0, w3
                    do g2 = 0, w2
                        do g1 = 0, w1
                            amr_cg(i)%sf(g1, g2, g3) = ${QP('i')}$
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            ! amr_cg is now device-current for the runtime C/F ghost-fill. Sync to host only when a host consumer follows.
            if (to_host) then
                do i = 1, sys_size
                    $:GPU_UPDATE(host='[amr_cg(i)%sf]')
                end do
            end if

        end subroutine s_amr_copy_parent_patch_${GSFX}$
    #:endfor

    !> Sub-box variants of the parent-patch pack/unpack/copy for the ring-clipped parent-fill wave (cons only - the wave ships
    !! q_cons; pb/mv runs keep the full-patch contract). Bounds are PATCH-LOCAL cell ranges; the buffer holds the sub-box in the
    !! same (g1 fastest, sys_size outermost) layout as the full-patch kernels, so both wire sides agree by construction.
    impure subroutine s_amr_pack_parent_box_device_cons(qp, bl, bh, buf)

        integer, intent(in)                 :: qp  !< parent's flat-store slot
        integer, intent(in)                 :: bl(3), bh(3)
        real(wp), intent(inout), contiguous :: buf(:)
        integer                             :: i, g1, g2, g3, o1, o2, o3, n1, n2, n3, l1, l2, l3, u1, u2, u3

        o1 = amr_cpat_off(1); o2 = amr_cpat_off(2); o3 = amr_cpat_off(3)
        l1 = bl(1); l2 = bl(2); l3 = bl(3); u1 = bh(1); u2 = bh(2); u3 = bh(3)
        n1 = u1 - l1 + 1; n2 = u2 - l2 + 1; n3 = u3 - l3 + 1
        $:GPU_PARALLEL_LOOP(collapse=4, copyout='[buf]')
        do i = 1, sys_size
            do g3 = l3, u3
                do g2 = l2, u2
                    do g1 = l1, u1
                        buf(1 + (g1 - l1) + n1*((g2 - l2) + n2*((g3 - l3) + n3*(i - 1)))) = real(amr_cons_st(g1 + o1, g2 + o2, &
                            & g3 + o3, i, qp), wp)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_pack_parent_box_device_cons

    impure subroutine s_amr_unpack_parent_box_device(bl, bh, buf)

        integer, intent(in)              :: bl(3), bh(3)
        real(wp), intent(in), contiguous :: buf(:)
        integer                          :: i, g1, g2, g3, n1, n2, n3, l1, l2, l3, u1, u2, u3

        l1 = bl(1); l2 = bl(2); l3 = bl(3); u1 = bh(1); u2 = bh(2); u3 = bh(3)
        n1 = u1 - l1 + 1; n2 = u2 - l2 + 1; n3 = u3 - l3 + 1
        $:GPU_PARALLEL_LOOP(collapse=4, copyin='[buf]')
        do i = 1, sys_size
            do g3 = l3, u3
                do g2 = l2, u2
                    do g1 = l1, u1
                        amr_cg(i)%sf(g1, g2, g3) = buf(1 + (g1 - l1) + n1*((g2 - l2) + n2*((g3 - l3) + n3*(i - 1))))
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_unpack_parent_box_device

    impure subroutine s_amr_copy_parent_box_cons(qp, bl, bh)

        integer, intent(in) :: qp  !< parent's flat-store slot
        integer, intent(in) :: bl(3), bh(3)
        integer             :: i, g1, g2, g3, o1, o2, o3, l1, l2, l3, u1, u2, u3

        o1 = amr_cpat_off(1); o2 = amr_cpat_off(2); o3 = amr_cpat_off(3)
        l1 = bl(1); l2 = bl(2); l3 = bl(3); u1 = bh(1); u2 = bh(2); u3 = bh(3)
        $:GPU_PARALLEL_LOOP(collapse=4)
        do i = 1, sys_size
            do g3 = l3, u3
                do g2 = l2, u2
                    do g1 = l1, u1
                        amr_cg(i)%sf(g1, g2, g3) = amr_cons_st(g1 + o1, g2 + o2, g3 + o3, i, qp)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_copy_parent_box_cons

    !> Rank r's coarse-grid decomposition (start_idx + local extent m/n/p), computed O(1) from its cartesian coords instead of the
    !! replicated amr_decomp table. The domain cart is MPI_CART_CREATE(reorder=.false., dims=[num_procs_x,num_procs_y,num_procs_z]),
    !! so ranks keep MPI_COMM_WORLD order and r's coords are row-major in r. The per-dim split is the integer split-with-remainder
    !! on CELL INDICES, identical to s_mpi_decompose_computational_domain (independent of grid stretching, which changes coords not
    !! indices).
    pure subroutine s_amr_rank_decomp(r, sidx, ext)

        integer, intent(in)  :: r
        integer, intent(out) :: sidx(3), ext(3)
        integer              :: coords(3), gd(3), pd(3), base, rem, d

        gd(1) = m_glb; gd(2) = n_glb; gd(3) = p_glb
        pd(1) = num_procs_x; pd(2) = num_procs_y; pd(3) = num_procs_z
        ! row-major: r = cx*(py*pz) + cy*pz + cz  (py=pz=1 when that dim is not decomposed)
        coords(1) = r/(pd(2)*pd(3))
        coords(2) = mod(r, pd(2)*pd(3))/pd(3)
        coords(3) = mod(r, pd(3))
        sidx = 0; ext = 0
        do d = 1, 3
            if (d == 2 .and. n_glb == 0) cycle
            if (d == 3 .and. p_glb == 0) cycle
            base = (gd(d) + 1)/pd(d)
            rem = mod(gd(d) + 1, pd(d))
            ext(d) = base - 1 + merge(1, 0, coords(d) < rem)
            sidx(d) = coords(d)*base + min(coords(d), rem)
        end do

    end subroutine s_amr_rank_decomp

    !> Permanent guard: the computed accessor must reproduce THIS rank's actual decomposition. Every rank checks its own entry, so
    !! collectively the formula is proven for all r; fires immediately if it ever drifts from s_mpi_decompose_computational_domain.
    !! O(1), always-on.
    impure subroutine s_amr_validate_decomp()

        integer :: sidx(3), ext(3)
        logical :: ok

        call s_amr_rank_decomp(proc_rank, sidx, ext)
        ! nested guards, NOT a single .and./.or.: Fortran does not short-circuit, so start_idx(2)/start_idx(3) (start_idx is sized
        ! num_dims) would be read out of bounds in 1D/2D even though the guard is false - a GNU-reldebug bounds abort.
        ok = (sidx(1) == start_idx(1) .and. ext(1) == m)
        if (n_glb > 0) then
            if (sidx(2) /= start_idx(2) .or. ext(2) /= n) ok = .false.
        end if
        if (p_glb > 0) then
            if (sidx(3) /= start_idx(3) .or. ext(3) /= p) ok = .false.
        end if
        if (.not. ok) then
            call s_mpi_abort('s_amr_rank_decomp does not reproduce this rank''s decomposition - computed split disagrees with ' &
                             & // 's_mpi_decompose_computational_domain')
        end if

    end subroutine s_amr_validate_decomp

    !> Closed-form inverse of the per-dim split-with-remainder used by s_amr_rank_decomp: the cart coord owning global coarse cell
    !! g. base=(gd+1)/pd, rem=mod(gd+1,pd). Exact for g in [0,gd]; callers clamp to [0,pd-1] for out-of-range g (ghost reach).
    !! base==0 (more ranks than cells) => the g<rem*(base+1) branch always holds and c=g.
    pure integer function f_amr_cell_coord(g, base, rem) result(c)

        integer, intent(in) :: g, base, rem

        if (g < rem*(base + 1)) then
            c = g/(base + 1)
        else
            c = rem + (g - rem*(base + 1))/base
        end if

    end function f_amr_cell_coord

    !> Per-dim contiguous rank-coord range [clo:chi] whose owned coarse slab intersects box [blo:bhi], clamped to [0,pd-1]. The
    !! boundary coords' coarse_range is extended by buff_size exactly at the domain edge (s_amr_rank_coarse_range), so this clamped
    !! interior-frame range reproduces BOTH the interior (scatter) and the coarse_range (gather) intersection sets: a box reaching
    !! the ghost zone clamps to the boundary coord whose extended slab contains it, and there is no rank beyond that coord.
    !! Collapsed dims (n_glb==0 / p_glb==0) contribute coord 0.
    pure subroutine s_amr_coord_range(blo, bhi, clo, chi)

        integer, intent(in)  :: blo(3), bhi(3)
        integer, intent(out) :: clo(3), chi(3)
        integer              :: gd(3), pd(3), base, rem, d

        gd(1) = m_glb; gd(2) = n_glb; gd(3) = p_glb
        pd(1) = num_procs_x; pd(2) = num_procs_y; pd(3) = num_procs_z
        clo = 0; chi = 0
        do d = 1, 3
            if (d == 2 .and. n_glb == 0) cycle
            if (d == 3 .and. p_glb == 0) cycle
            base = (gd(d) + 1)/pd(d)
            rem = mod(gd(d) + 1, pd(d))
            clo(d) = min(max(f_amr_cell_coord(blo(d), base, rem), 0), pd(d) - 1)
            chi(d) = min(max(f_amr_cell_coord(bhi(d), base, rem), 0), pd(d) - 1)
        end do

    end subroutine s_amr_coord_range

    !> Ascending rank list overlapping coarse box [blo:bhi]. Enumerates the coord brick cx->cy->cz so r = cx*(Py*Pz)+cy*Pz+cz is
    !! monotonic => ascending, reproducing the r=0..num_procs-1 scan order. Owner NOT excluded (consumers keep their own owner
    !! skip). Caller sizes ranks(:) >= f_amr_overlap_count(blo,bhi).
    pure subroutine s_amr_ranks_overlapping(blo, bhi, ranks, nr)

        integer, intent(in)  :: blo(3), bhi(3)
        integer, intent(out) :: ranks(:)
        integer, intent(out) :: nr
        integer              :: clo(3), chi(3), cx, cy, cz, pyz

        call s_amr_coord_range(blo, bhi, clo, chi)
        pyz = num_procs_y*num_procs_z
        nr = 0
        do cx = clo(1), chi(1)
            do cy = clo(2), chi(2)
                do cz = clo(3), chi(3)
                    nr = nr + 1
                    ranks(nr) = cx*pyz + cy*num_procs_z + cz
                end do
            end do
        end do

    end subroutine s_amr_ranks_overlapping

    !> Overlap count only (allocation sizing), = product of the per-dim coord-range widths.
    pure integer function f_amr_overlap_count(blo, bhi) result(nr)

        integer, intent(in) :: blo(3), bhi(3)
        integer             :: clo(3), chi(3)

        call s_amr_coord_range(blo, bhi, clo, chi)
        nr = (chi(1) - clo(1) + 1)*(chi(2) - clo(2) + 1)*(chi(3) - clo(3) + 1)

    end function f_amr_overlap_count

    !> Rank r's contiguous owned coarse-cell range per dim from the computed decomposition (s_amr_rank_decomp): interior
    !! [start:start+ext] plus its physical-boundary ghosts (buff_size cells only where the subdomain touches the domain edge). Equal
    !! to the set where f_amr_own_coarse is true, but as one contiguous span so box intersections identify contributors without a
    !! per-cell scan.
    pure subroutine s_amr_rank_coarse_range(r, crlo, crhi)

        integer, intent(in)  :: r
        integer, intent(out) :: crlo(3), crhi(3)
        integer              :: sidx(3), ext(3)

        call s_amr_rank_decomp(r, sidx, ext)
        crlo = 0; crhi = 0
        crlo(1) = sidx(1); if (sidx(1) == 0) crlo(1) = -buff_size
        crhi(1) = sidx(1) + ext(1); if (crhi(1) == m_glb) crhi(1) = crhi(1) + buff_size
        if (n_glb > 0) then
            crlo(2) = sidx(2); if (sidx(2) == 0) crlo(2) = -buff_size
            crhi(2) = sidx(2) + ext(2); if (crhi(2) == n_glb) crhi(2) = crhi(2) + buff_size
        end if
        if (p_glb > 0) then
            crlo(3) = sidx(3); if (sidx(3) == 0) crlo(3) = -buff_size
            crhi(3) = sidx(3) + ext(3); if (crhi(3) == p_glb) crhi(3) = crhi(3) + buff_size
        end if

    end subroutine s_amr_rank_coarse_range

    !> Per-dim intersection of two global boxes [alo:ahi] and [blo:bhi] -> [olo:ohi] (empty when olo > ohi in some dim).
    pure subroutine s_amr_box_isect(alo, ahi, blo, bhi, olo, ohi)

        integer, intent(in)  :: alo(3), ahi(3), blo(3), bhi(3)
        integer, intent(out) :: olo(3), ohi(3)

        olo = max(alo, blo); ohi = min(ahi, bhi)

    end subroutine s_amr_box_isect

    !> Do two coarse-index boxes [alo:ahi] and [blo:bhi] overlap? Collapsed dims (n_glb/p_glb == 0) never disqualify.
    pure logical function f_amr_boxes_overlap(alo, ahi, blo, bhi) result(ov)

        integer, intent(in) :: alo(3), ahi(3), blo(3), bhi(3)

        ov = alo(1) <= bhi(1) .and. ahi(1) >= blo(1)
        if (n_glb > 0) ov = ov .and. alo(2) <= bhi(2) .and. ahi(2) >= blo(2)
        if (p_glb > 0) ov = ov .and. alo(3) <= bhi(3) .and. ahi(3) >= blo(3)

    end function f_amr_boxes_overlap

    !> Multi-level nesting: index of the covering level-(level(k)-1) block that block k refines - its coarse parent - or 0 when
    !! block k is level 1 (parent is the L0 base grid). Regions are in L0 cell indices at every level, so the parent is the
    !! level-below block whose box contains k's; proper nesting guarantees exactly one, and the first overlap is returned.
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

    !> Copy this rank's own coarse cells (box [bl:bh] GLOBAL, read from q_coarse at its own start-idx frame o1/o2/o3) into amr_cg in
    !! the block-local patch frame. stp -> stp, exact.
    impure subroutine s_amr_unpack_patch(q_coarse, bl, bh, o1, o2, o3)

        type(scalar_field), dimension(sys_size), intent(in) :: q_coarse
        integer, intent(in)                                 :: bl(3), bh(3), o1, o2, o3
        integer                                             :: i, g1, g2, g3

        do i = 1, sys_size
            do g3 = bl(3), bh(3)
                do g2 = bl(2), bh(2)
                    do g1 = bl(1), bh(1)
                        amr_cg(i)%sf(g1 - amr_cpat_off(1), g2 - amr_cpat_off(2), g3 - amr_cpat_off(3)) = q_coarse(i)%sf(g1 - o1, &
                               & g2 - o2, g3 - o3)
                    end do
                end do
            end do
        end do

    end subroutine s_amr_unpack_patch

    !> Runtime device analogue of s_amr_unpack_patch: copy the owner's own coarse box [bl:bh] GLOBAL from q_coarse (device) into
    !! amr_cg (device) in the patch-local frame - no host round-trip. Same index map and direct stp assignment as the host path.
    !! TWIN s_amr_gather_own_box_pbmv_device (q<->pb/mv): same own-box index map; keep lockstep.
    impure subroutine s_amr_gather_own_box_device(q_coarse, bl, bh, o1, o2, o3)

        type(scalar_field), dimension(sys_size), intent(in) :: q_coarse
        integer, intent(in)                                 :: bl(3), bh(3), o1, o2, o3
        integer                                             :: i, g1, g2, g3, bl1, bl2, bl3, bh1, bh2, bh3, coff1, coff2, coff3

        ! scalar copies: no host array may be referenced inside the device region (nvfortran/Cray demand it PRESENT)

        bl1 = bl(1); bh1 = bh(1); bl2 = bl(2); bh2 = bh(2); bl3 = bl(3); bh3 = bh(3)
        coff1 = amr_cpat_off(1); coff2 = amr_cpat_off(2); coff3 = amr_cpat_off(3)
        $:GPU_PARALLEL_LOOP(collapse=4)
        do i = 1, sys_size
            do g3 = bl3, bh3
                do g2 = bl2, bh2
                    do g1 = bl1, bh1
                        amr_cg(i)%sf(g1 - coff1, g2 - coff2, g3 - coff3) = q_coarse(i)%sf(g1 - o1, g2 - o2, g3 - o3)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_gather_own_box_device

    !> Runtime device pack of the overlap box [bl:bh] GLOBAL from q_coarse (device) into the contiguous wire buffer buf (host, via
    !! copyout) - only the box crosses PCIe, not the full field. Explicit-loop linear buf indexing (g1 fastest, then g2, g3, i) and
    !! the wp cast match the host pack in s_amr_gather_coarse_patch element-for-element, so the receiver's unpack is layout- and
    !! byte-identical (same discipline as s_amr_fine_slice: no array-section syntax near the device map). TWIN
    !! s_amr_pack_box_pbmv_device (q<->pb/mv): same wire linear order + wp cast; keep lockstep.
    impure subroutine s_amr_pack_box_device(q_coarse, bl, bh, o1, o2, o3, buf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_coarse
        integer, intent(in)                                 :: bl(3), bh(3), o1, o2, o3
        real(wp), intent(inout), contiguous                 :: buf(:)
        integer                                             :: i, g1, g2, g3, bl1, bl2, bl3, bh1, bh2, bh3, n1, n2, n3

        bl1 = bl(1); bh1 = bh(1); bl2 = bl(2); bh2 = bh(2); bl3 = bl(3); bh3 = bh(3)
        n1 = bh1 - bl1 + 1; n2 = bh2 - bl2 + 1; n3 = bh3 - bl3 + 1
        $:GPU_PARALLEL_LOOP(collapse=4, copyout='[buf]')
        do i = 1, sys_size
            do g3 = bl3, bh3
                do g2 = bl2, bh2
                    do g1 = bl1, bh1
                        buf(1 + (g1 - bl1) + n1*((g2 - bl2) + n2*((g3 - bl3) + n3*(i - 1)))) = real(q_coarse(i)%sf(g1 - o1, &
                            & g2 - o2, g3 - o3), wp)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_pack_box_device

    !> Runtime device unpack of a received overlap box [bl:bh] GLOBAL from the contiguous wire buffer buf (host, via copyin) into
    !! amr_cg (device) in the patch-local frame - only the box crosses PCIe. Same linear order and stp cast as the host unpack in
    !! s_amr_gather_coarse_patch. TWIN s_amr_unpack_box_pbmv_device (q<->pb/mv): same wire linear order + stp cast; keep lockstep.
    impure subroutine s_amr_unpack_box_device(bl, bh, buf)

        integer, intent(in)              :: bl(3), bh(3)
        real(wp), intent(in), contiguous :: buf(:)
        integer                          :: i, g1, g2, g3, bl1, bl2, bl3, bh1, bh2, bh3, n1, n2, n3, coff1, coff2, coff3

        bl1 = bl(1); bh1 = bh(1); bl2 = bl(2); bh2 = bh(2); bl3 = bl(3); bh3 = bh(3)
        n1 = bh1 - bl1 + 1; n2 = bh2 - bl2 + 1; n3 = bh3 - bl3 + 1
        coff1 = amr_cpat_off(1); coff2 = amr_cpat_off(2); coff3 = amr_cpat_off(3)
        $:GPU_PARALLEL_LOOP(collapse=4, copyin='[buf]')
        do i = 1, sys_size
            do g3 = bl3, bh3
                do g2 = bl2, bh2
                    do g1 = bl1, bh1
                        amr_cg(i)%sf(g1 - coff1, g2 - coff2, &
                               & g3 - coff3) = real(buf(1 + (g1 - bl1) + n1*((g2 - bl2) + n2*((g3 - bl3) + n3*(i - 1)))), stp)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_unpack_box_device

    !> Runtime device own-box copy for the pbmv gather: pb/mv (device) -> amr_cg_pb/mv (device) over [bl:bh] GLOBAL in the
    !! patch-local frame. Same index map and direct stp assignment as the host path in s_amr_gather_coarse_patch_pbmv. TWIN
    !! s_amr_gather_own_box_device (pb/mv<->q): q_cons sibling of this own-box copy; keep the index map lockstep.
    impure subroutine s_amr_gather_own_box_pbmv_device(pb_coarse, mv_coarse, bl, bh, o1, o2, o3)

        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(in) :: pb_coarse, mv_coarse
        integer, intent(in) :: bl(3), bh(3), o1, o2, o3
        integer :: q, ib_, g1, g2, g3, bl1, bl2, bl3, bh1, bh2, bh3, coff1, coff2, coff3

        bl1 = bl(1); bh1 = bh(1); bl2 = bl(2); bh2 = bh(2); bl3 = bl(3); bh3 = bh(3)
        coff1 = amr_cpat_off(1); coff2 = amr_cpat_off(2); coff3 = amr_cpat_off(3)
        $:GPU_PARALLEL_LOOP(collapse=5)
        do ib_ = 1, nb
            do q = 1, nnode
                do g3 = bl3, bh3
                    do g2 = bl2, bh2
                        do g1 = bl1, bh1
                            amr_cg_pb(g1 - coff1, g2 - coff2, g3 - coff3, q, ib_) = pb_coarse(g1 - o1, g2 - o2, g3 - o3, q, ib_)
                            amr_cg_mv(g1 - coff1, g2 - coff2, g3 - coff3, q, ib_) = mv_coarse(g1 - o1, g2 - o2, g3 - o3, q, ib_)
                        end do
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_gather_own_box_pbmv_device

    !> Runtime device pack for the pbmv gather: pb block then mv block of the overlap box [bl:bh] GLOBAL into the contiguous wire
    !! buffer buf (host, via copyout). Linear order (g1 fastest, then g2, g3, q, ib_; mv offset by half the message) and wp cast
    !! match the host pack in s_amr_gather_coarse_patch_pbmv element-for-element. TWIN s_amr_pack_box_device (pb/mv<->q): q_cons
    !! sibling; keep the wire linear order + wp cast lockstep.
    impure subroutine s_amr_pack_box_pbmv_device(pb_coarse, mv_coarse, bl, bh, o1, o2, o3, buf)

        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(in) :: pb_coarse, mv_coarse
        integer, intent(in) :: bl(3), bh(3), o1, o2, o3
        real(wp), intent(inout), contiguous :: buf(:)
        integer :: q, ib_, g1, g2, g3, bl1, bl2, bl3, bh1, bh2, bh3, n1, n2, n3, half

        bl1 = bl(1); bh1 = bh(1); bl2 = bl(2); bh2 = bh(2); bl3 = bl(3); bh3 = bh(3)
        n1 = bh1 - bl1 + 1; n2 = bh2 - bl2 + 1; n3 = bh3 - bl3 + 1
        half = n1*n2*n3*nnode*nb
        $:GPU_PARALLEL_LOOP(collapse=5, copyout='[buf]')
        do ib_ = 1, nb
            do q = 1, nnode
                do g3 = bl3, bh3
                    do g2 = bl2, bh2
                        do g1 = bl1, bh1
                            buf(1 + (g1 - bl1) + n1*((g2 - bl2) + n2*((g3 - bl3) + n3*((q - 1) + nnode*(ib_ - 1))))) &
                                & = real(pb_coarse(g1 - o1, g2 - o2, g3 - o3, q, ib_), wp)
                            buf(half + 1 + (g1 - bl1) + n1*((g2 - bl2) + n2*((g3 - bl3) + n3*((q - 1) + nnode*(ib_ - 1))))) &
                                & = real(mv_coarse(g1 - o1, g2 - o2, g3 - o3, q, ib_), wp)
                        end do
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_pack_box_pbmv_device

    !> Runtime device unpack for the pbmv gather: received wire buffer buf (host, via copyin) -> amr_cg_pb/mv (device) over [bl:bh]
    !! GLOBAL in the patch-local frame. Same linear order and stp cast as the host unpack in s_amr_gather_coarse_patch_pbmv. TWIN
    !! s_amr_unpack_box_device (pb/mv<->q): q_cons sibling; keep the wire linear order + stp cast lockstep.
    impure subroutine s_amr_unpack_box_pbmv_device(bl, bh, buf)

        integer, intent(in)              :: bl(3), bh(3)
        real(wp), intent(in), contiguous :: buf(:)
        integer                          :: q, ib_, g1, g2, g3, bl1, bl2, bl3, bh1, bh2, bh3, n1, n2, n3, half, coff1, coff2, coff3

        bl1 = bl(1); bh1 = bh(1); bl2 = bl(2); bh2 = bh(2); bl3 = bl(3); bh3 = bh(3)
        n1 = bh1 - bl1 + 1; n2 = bh2 - bl2 + 1; n3 = bh3 - bl3 + 1
        half = n1*n2*n3*nnode*nb
        coff1 = amr_cpat_off(1); coff2 = amr_cpat_off(2); coff3 = amr_cpat_off(3)
        $:GPU_PARALLEL_LOOP(collapse=5, copyin='[buf]')
        do ib_ = 1, nb
            do q = 1, nnode
                do g3 = bl3, bh3
                    do g2 = bl2, bh2
                        do g1 = bl1, bh1
                            amr_cg_pb(g1 - coff1, g2 - coff2, g3 - coff3, q, &
                                      & ib_) = real(buf(1 + (g1 - bl1) + n1*((g2 - bl2) + n2*((g3 - bl3) + n3*((q - 1) &
                                      & + nnode*(ib_ - 1))))), stp)
                            amr_cg_mv(g1 - coff1, g2 - coff2, g3 - coff3, q, &
                                      & ib_) = real(buf(half + 1 + (g1 - bl1) + n1*((g2 - bl2) + n2*((g3 - bl3) + n3*((q - 1) &
                                      & + nnode*(ib_ - 1))))), stp)
                        end do
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_unpack_box_pbmv_device

    !> True iff rank r is a reflux applier for the current block: it owns the coarse cell layer just OUTSIDE some block face AND its
    !! subdomain overlaps the block transversely. Mirrors s_amr_reflux_face_flags, but parameterized by r's subdomain from the
    !! computed decomposition (s_amr_rank_decomp, so the block owner can decide which ranks to send freg to, and each rank agrees on
    !! whether it receives). Uses amr_region_lo/hi (the current block, set on every rank by s_amr_select_slot). NOTE: deliberately
    !! NO f_amr_face_is_seam clip (unlike the flags) - and the participation-map build (s_amr_reg_prepare, m_amr_registers) copies
    !! THIS unclipped formula for its clause (c), because both exchange paths gate their freg receives on it. Keep lockstep.
    pure logical function f_amr_reflux_participates(r) result(part)

        integer, intent(in) :: r
        integer             :: sidx(3), ext(3), d, t
        logical             :: tv(3), tvd

        call s_amr_rank_decomp(r, sidx, ext)
        tv(1) = amr_region_lo(1) <= sidx(1) + ext(1) .and. amr_region_hi(1) >= sidx(1)
        tv(2) = (n_glb == 0) .or. (amr_region_lo(2) <= sidx(2) + ext(2) .and. amr_region_hi(2) >= sidx(2))
        tv(3) = (p_glb == 0) .or. (amr_region_lo(3) <= sidx(3) + ext(3) .and. amr_region_hi(3) >= sidx(3))
        part = .false.
        do d = 1, num_dims
            tvd = .true.
            do t = 1, num_dims
                if (t /= d) tvd = tvd .and. tv(t)
            end do
            if (tvd .and. amr_region_lo(d) - 1 >= sidx(d) .and. amr_region_lo(d) - 1 <= sidx(d) + ext(d)) part = .true.
            if (tvd .and. amr_region_hi(d) + 1 >= sidx(d) .and. amr_region_hi(d) + 1 <= sidx(d) + ext(d)) part = .true.
        end do

    end function f_amr_reflux_participates

    !> Per-face refinement of f_amr_reflux_participates: the faces of the CURRENT block that rank r actually APPLIES - it owns the
    !! outside coarse layer with transverse overlap, minus fine-fine seam faces - mirroring s_amr_reflux_face_flags term for term
    !! (same ownership formula, same f_amr_face_is_seam exclusion). The reflux-faces wave ships exactly these: sender and every
    !! receiver derive the identical set from replicated data, so a face a rank never applies never rides the wire.
    pure subroutine s_amr_reflux_faces_for(r, s_lo, s_hi)

        integer, intent(in)  :: r
        logical, intent(out) :: s_lo(3), s_hi(3)
        integer              :: sidx(3), ext(3), d, t
        logical              :: tv(3), tvd

        call s_amr_rank_decomp(r, sidx, ext)
        tv(1) = amr_region_lo(1) <= sidx(1) + ext(1) .and. amr_region_hi(1) >= sidx(1)
        tv(2) = (n_glb == 0) .or. (amr_region_lo(2) <= sidx(2) + ext(2) .and. amr_region_hi(2) >= sidx(2))
        tv(3) = (p_glb == 0) .or. (amr_region_lo(3) <= sidx(3) + ext(3) .and. amr_region_hi(3) >= sidx(3))
        s_lo = .false.; s_hi = .false.
        do d = 1, num_dims
            tvd = .true.
            do t = 1, num_dims
                if (t /= d) tvd = tvd .and. tv(t)
            end do
            s_lo(d) = tvd .and. amr_region_lo(d) - 1 >= sidx(d) .and. amr_region_lo(d) - 1 <= sidx(d) + ext(d) &
                 & .and. .not. f_amr_face_is_seam(d, -1)
            s_hi(d) = tvd .and. amr_region_hi(d) + 1 >= sidx(d) .and. amr_region_hi(d) + 1 <= sidx(d) + ext(d) &
                 & .and. .not. f_amr_face_is_seam(d, 1)
        end do

    end subroutine s_amr_reflux_faces_for

    !> Fine-level distribution: deliver the current block's fine flux registers freg (captured by the owner during the fine advance)
    !! to exactly the (SFC-local) coarse-outside-owners that apply the reflux - POINT-TO-POINT, replacing the global broadcast. The
    !! owner sends its whole freg slot (block-relative; each applier reads its own transverse slice) to every participant; non-owner
    !! participants receive it. Device-resident: owner stages its slot to host, receivers push it back. No-op without MPI/at np=1.
    impure subroutine s_amr_p2p_reflux_faces()

#ifdef MFC_MPI
        integer              :: owner, r, ierr, nreq, cnt, idx, ncand
        integer              :: cand(num_procs), glo(3), ghi(3)
        integer, allocatable :: reqs(:)

        if (.not. amr) return
        if (num_procs == 1) return
        owner = amr_block_owner(amr_cur)
        if (proc_rank == owner) then
            #:for D in [1, 2, 3]
                if (${D}$ <= num_dims) then
                    $:GPU_UPDATE(host='[freg(' + str(D) + ')%lo(:, :, :, amr_reg_cur), freg(' + str(D) &
                                 & + ')%hi(:, :, :, amr_reg_cur)]')
                end if
            #:endfor
            ! participating ranks by O(overlap) inversion (region grown by 1) filtered by the UNCHANGED predicate, in place of the
            ! O(P) rank scan. Ascending (owner-excluded at use), so the ISENDs match the receivers exactly as the scan did.
            glo = 0; ghi = 0
            glo(1) = amr_region_lo(1) - 1; ghi(1) = amr_region_hi(1) + 1
            if (n_glb > 0) then; glo(2) = amr_region_lo(2) - 1; ghi(2) = amr_region_hi(2) + 1; end if
            if (p_glb > 0) then; glo(3) = amr_region_lo(3) - 1; ghi(3) = amr_region_hi(3) + 1; end if
            call s_amr_ranks_overlapping(glo, ghi, cand, ncand)
            nreq = 0
            do idx = 1, ncand
                r = cand(idx)
                if (r /= owner .and. f_amr_reflux_participates(r)) nreq = nreq + 1
            end do
            if (nreq > 0) then
                allocate (reqs(2*num_dims*nreq))
                nreq = 0
                do idx = 1, ncand
                    r = cand(idx)
                    if (r == owner .or. .not. f_amr_reflux_participates(r)) cycle
                    #:for D in [1, 2, 3]
                        if (${D}$ <= num_dims) then
                            cnt = size(freg(${D}$)%lo(:,:,:,amr_reg_cur))
                            nreq = nreq + 1
                            call s_xa_rec(XA_F5_FACE_SND, 1, cnt, ${2*D}$)
                            call MPI_ISEND(freg(${D}$)%lo(:,:,:,amr_reg_cur), cnt, mpi_p, r, ${2*D}$, MPI_COMM_WORLD, reqs(nreq), &
                                           & ierr)
                            nreq = nreq + 1
                            call s_xa_rec(XA_F5_FACE_SND, 1, cnt, ${2*D + 1}$)
                            call MPI_ISEND(freg(${D}$)%hi(:,:,:,amr_reg_cur), cnt, mpi_p, r, ${2*D + 1}$, MPI_COMM_WORLD, &
                                           & reqs(nreq), ierr)
                        end if
                    #:endfor
                end do
                call s_phase_tic(PH_RFWAIT)
                call MPI_WAITALL(nreq, reqs, MPI_STATUSES_IGNORE, ierr)
                call s_phase_toc(PH_RFWAIT)
                deallocate (reqs)
            end if
        else if (f_amr_reflux_participates(proc_rank)) then
            ! Post all 2*num_dims receives, then ONE wait. The blocking form serialised the six faces
            ! against the owner's send order and cost 6.3%% of wall (rf:recv, 4010 calls). The slice
            ! (:,:,:,amr_reg_cur) is contiguous - the dense register slot is the last dimension - and the owner side
            ! already ISENDs the identical shape, so no temporary-buffer hazard is introduced.
            call s_phase_tic(PH_RFRECV)
            allocate (reqs(2*num_dims))
            nreq = 0
            #:for D in [1, 2, 3]
                if (${D}$ <= num_dims) then
                    cnt = size(freg(${D}$)%lo(:,:,:,amr_reg_cur))
                    nreq = nreq + 1
                    call s_xa_rec(XA_F5_FACE_RCV, 2, cnt, ${2*D}$)
                    call MPI_IRECV(freg(${D}$)%lo(:,:,:,amr_reg_cur), cnt, mpi_p, owner, ${2*D}$, MPI_COMM_WORLD, reqs(nreq), ierr)
                    nreq = nreq + 1
                    call s_xa_rec(XA_F5_FACE_RCV, 2, cnt, ${2*D + 1}$)
                    call MPI_IRECV(freg(${D}$)%hi(:,:,:,amr_reg_cur), cnt, mpi_p, owner, ${2*D + 1}$, MPI_COMM_WORLD, reqs(nreq), &
                                   & ierr)
                end if
            #:endfor
            call MPI_WAITALL(nreq, reqs, MPI_STATUSES_IGNORE, ierr)
            deallocate (reqs)
            ! Device update only AFTER the wait: the buffers hold nothing valid until then.
            #:for D in [1, 2, 3]
                if (${D}$ <= num_dims) then
                    $:GPU_UPDATE(device='[freg(' + str(D) + ')%lo(:, :, :, amr_reg_cur), freg(' + str(D) &
                                 & + ')%hi(:, :, :, amr_reg_cur)]')
                end if
            #:endfor
            call s_phase_toc(PH_RFRECV)
        end if
#endif

    end subroutine s_amr_p2p_reflux_faces

    !> I5-F5a: the per-stage level-1 reflux-face exchange as ONE wave (amr_plan_based_exchange.md). The per-box form posted and
    !! WAITALLed per box on both sides - an O(boxes) rendezvous chain. Here every rank walks the level-1 slots ascending: all
    !! receives post first (ZERO-COPY, directly into the freg host mirrors - each box owns a register slot, so no pool is needed and
    !! the message count is unchanged BY DESIGN), then the owners stage + multicast, then ONE waitall, then the receivers push to
    !! device. All messages between a fixed (owner, participant) pair share tag amr_tag_base(5)+epoch and match by non-overtaking in
    !! the (ascending box, ascending dim, lo-then-hi) order both sides derive from the same predicates the per-box form used
    !! (amr_block_owner + the replicated region mirrors + the pure participation test). Under MFC_DEBUG the identity headers travel
    !! as separate 8-word COMPANION messages, one per (box, peer) group ahead of its payloads (a prefix cannot ride a zero-copy
    !! payload); they are never recorded in [amr-xa], so the family words stay exactly comparable. The register arrays are sized UP
    !! FRONT: the apply can REALLOCATE them, so nothing may post into freg before s_amr_reg_prepare.
    impure subroutine s_amr_reflux_faces_wave()

#ifdef MFC_MPI
        use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        integer  :: k, r, ierr, nreq, cnt, idx, ncand, tq, nhr, nhs, j
        integer  :: cand(num_procs), glo(3), ghi(3)
        logical  :: s_lo(3), s_hi(3), u_lo(3), u_hi(3)
        logical  :: cl(3, num_procs), ch(3, num_procs)
        real(wp) :: nanv

        if (num_procs == 1) return
        call s_amr_reg_prepare()
        tq = amr_tag_base(5) + int(mod(amr_mesh_epoch, 50_8))
        nanv = ieee_value(0._wp, ieee_quiet_nan)
        nreq = 0; nhr = 0; nhs = 0
        do k = 1, amr_num_blocks
            if (amr_block_level(k) /= 1) cycle
            call s_amr_select_slot(k)
            if (amr_block_owner(k) == proc_rank) cycle
            if (.not. f_amr_reflux_participates(proc_rank)) cycle
            ! face-selective multicast: receive exactly the faces THIS rank applies (s_amr_reflux_faces_for mirrors the
            ! apply's own_lo/own_hi + seam gates); the owner derives the same set per participant, so the pairing is
            ! exact with no metadata exchange. Debug arm: unreceived faces are NaN-flooded so any hidden reader aborts.
            call s_amr_reflux_faces_for(proc_rank, s_lo, s_hi)
            if (XA_NH > 0) then
                nhr = nhr + 1
                call s_amr_fw_szr(amr_fw_rq, XA_NH*nhr); call s_amr_fw_szi(amr_fw_rblk, nhr)
                amr_fw_rblk(nhr) = k
                nreq = nreq + 1
                call s_amr_fw_szi(amr_fw_req, nreq); call s_amr_fw_szi(amr_fw_reqw, nreq)
                amr_fw_reqw(nreq) = XA_NH
                call MPI_IRECV(amr_fw_rq(XA_NH*(nhr - 1) + 1), XA_NH, mpi_p, amr_block_owner(k), tq, MPI_COMM_WORLD, &
                               & amr_fw_req(nreq), ierr)
            end if
            #:for D in [1, 2, 3]
                if (${D}$ <= num_dims) then
                    cnt = size(freg(${D}$)%lo(:,:,:,amr_reg_cur))
                    if (s_lo(${D}$)) then
                        nreq = nreq + 1
                        call s_amr_fw_szi(amr_fw_req, nreq); call s_amr_fw_szi(amr_fw_reqw, nreq)
                        amr_fw_reqw(nreq) = cnt
                        call s_xa_rec(XA_F5W_FACE_RCV, 2, cnt, tq)
                        call MPI_IRECV(freg(${D}$)%lo(:,:,:,amr_reg_cur), cnt, mpi_p, amr_block_owner(k), tq, MPI_COMM_WORLD, &
                                       & amr_fw_req(nreq), ierr)
#ifdef MFC_DEBUG
                    else
                        freg(${D}$)%lo(:,:,:,amr_reg_cur) = nanv
                        $:GPU_UPDATE(device='[freg(' + str(D) + ')%lo(:, :, :, amr_reg_cur)]')
#endif
                    end if
                    if (s_hi(${D}$)) then
                        nreq = nreq + 1
                        call s_amr_fw_szi(amr_fw_req, nreq); call s_amr_fw_szi(amr_fw_reqw, nreq)
                        amr_fw_reqw(nreq) = cnt
                        call s_xa_rec(XA_F5W_FACE_RCV, 2, cnt, tq)
                        call MPI_IRECV(freg(${D}$)%hi(:,:,:,amr_reg_cur), cnt, mpi_p, amr_block_owner(k), tq, MPI_COMM_WORLD, &
                                       & amr_fw_req(nreq), ierr)
#ifdef MFC_DEBUG
                    else
                        freg(${D}$)%hi(:,:,:,amr_reg_cur) = nanv
                        $:GPU_UPDATE(device='[freg(' + str(D) + ')%hi(:, :, :, amr_reg_cur)]')
#endif
                    end if
                end if
            #:endfor
        end do
        do k = 1, amr_num_blocks
            if (amr_block_level(k) /= 1) cycle
            call s_amr_select_slot(k)
            if (amr_block_owner(k) /= proc_rank) cycle
            glo = 0; ghi = 0
            glo(1) = amr_region_lo(1) - 1; ghi(1) = amr_region_hi(1) + 1
            if (n_glb > 0) then; glo(2) = amr_region_lo(2) - 1; ghi(2) = amr_region_hi(2) + 1; end if
            if (p_glb > 0) then; glo(3) = amr_region_lo(3) - 1; ghi(3) = amr_region_hi(3) + 1; end if
            call s_amr_ranks_overlapping(glo, ghi, cand, ncand)
            ! face-selective multicast: each participant's ship set is its apply set (s_amr_reflux_faces_for), derived
            ! here per candidate; the device->host pull covers only the UNION of shipped faces (previously every owned
            ! block pulled all six faces every stage, participants or not).
            u_lo = .false.; u_hi = .false.
            do idx = 1, ncand
                r = cand(idx)
                cl(:,idx) = .false.; ch(:,idx) = .false.
                if (r == proc_rank .or. .not. f_amr_reflux_participates(r)) cycle
                call s_amr_reflux_faces_for(r, s_lo, s_hi)
                cl(:,idx) = s_lo; ch(:,idx) = s_hi
                u_lo = u_lo .or. s_lo; u_hi = u_hi .or. s_hi
            end do
            #:for D in [1, 2, 3]
                if (${D}$ <= num_dims) then
                    if (u_lo(${D}$)) then
                        $:GPU_UPDATE(host='[freg(' + str(D) + ')%lo(:, :, :, amr_reg_cur)]')
                    end if
                    if (u_hi(${D}$)) then
                        $:GPU_UPDATE(host='[freg(' + str(D) + ')%hi(:, :, :, amr_reg_cur)]')
                    end if
                end if
            #:endfor
            do idx = 1, ncand
                r = cand(idx)
                if (r == proc_rank .or. .not. f_amr_reflux_participates(r)) cycle
                if (XA_NH > 0) then
                    nhs = nhs + 1
                    call s_amr_fw_szr(amr_fw_sq, XA_NH*nhs)
                    call s_xa_hdr_pack(amr_fw_sq(XA_NH*(nhs - 1) + 1:XA_NH*nhs), XA_F5W_FACE_SND, k, [0, 0, 0], [0, 0, 0])
                    nreq = nreq + 1
                    call s_amr_fw_szi(amr_fw_req, nreq); call s_amr_fw_szi(amr_fw_reqw, nreq)
                    amr_fw_reqw(nreq) = -1
                    call MPI_ISEND(amr_fw_sq(XA_NH*(nhs - 1) + 1), XA_NH, mpi_p, r, tq, MPI_COMM_WORLD, amr_fw_req(nreq), ierr)
                end if
                #:for D in [1, 2, 3]
                    if (${D}$ <= num_dims) then
                        cnt = size(freg(${D}$)%lo(:,:,:,amr_reg_cur))
                        if (cl(${D}$, idx)) then
                            nreq = nreq + 1
                            call s_amr_fw_szi(amr_fw_req, nreq); call s_amr_fw_szi(amr_fw_reqw, nreq)
                            amr_fw_reqw(nreq) = -1
                            call s_xa_rec(XA_F5W_FACE_SND, 1, cnt, tq)
                            call MPI_ISEND(freg(${D}$)%lo(:,:,:,amr_reg_cur), cnt, mpi_p, r, tq, MPI_COMM_WORLD, &
                                           & amr_fw_req(nreq), ierr)
                        end if
                        if (ch(${D}$, idx)) then
                            nreq = nreq + 1
                            call s_amr_fw_szi(amr_fw_req, nreq); call s_amr_fw_szi(amr_fw_reqw, nreq)
                            amr_fw_reqw(nreq) = -1
                            call s_xa_rec(XA_F5W_FACE_SND, 1, cnt, tq)
                            call MPI_ISEND(freg(${D}$)%hi(:,:,:,amr_reg_cur), cnt, mpi_p, r, tq, MPI_COMM_WORLD, &
                                           & amr_fw_req(nreq), ierr)
                        end if
                    end if
                #:endfor
            end do
        end do
        if (nreq > 0) then
#ifdef MFC_DEBUG
            block
                integer :: st(MPI_STATUS_SIZE, nreq), gotw, q
                call s_phase_tic(PH_RFWAIT)
                call MPI_WAITALL(nreq, amr_fw_req, st, ierr)
                call s_phase_toc(PH_RFWAIT)
                do q = 1, nreq
                    if (amr_fw_reqw(q) < 0) cycle
                    call MPI_GET_COUNT(st(:,q), mpi_p, gotw, ierr)
                    @:ASSERT(gotw == amr_fw_reqw(q), "reflux-faces wave: received message length differs from the plan")
                end do
            end block
#else
            call s_phase_tic(PH_RFWAIT)
            call MPI_WAITALL(nreq, amr_fw_req, MPI_STATUSES_IGNORE, ierr)
            call s_phase_toc(PH_RFWAIT)
#endif
        end if
        call s_phase_tic(PH_RFRECV)
        j = 0
        do k = 1, amr_num_blocks
            if (amr_block_level(k) /= 1) cycle
            call s_amr_select_slot(k)
            if (amr_block_owner(k) == proc_rank) cycle
            if (.not. f_amr_reflux_participates(proc_rank)) cycle
            if (XA_NH > 0) then
                j = j + 1
                @:ASSERT(amr_fw_rblk(j) == k, "reflux-faces wave: header slot order broke")
                call s_xa_hdr_check(amr_fw_rq(XA_NH*(j - 1) + 1:XA_NH*j), XA_F5W_FACE_SND, k, [0, 0, 0], [0, 0, 0])
            end if
            ! push only the received faces; an unreceived face keeps its device content (never applied here - and
            ! NaN-poisoned in debug, so any hidden reader aborts)
            call s_amr_reflux_faces_for(proc_rank, s_lo, s_hi)
            #:for D in [1, 2, 3]
                if (${D}$ <= num_dims) then
                    if (s_lo(${D}$)) then
                        $:GPU_UPDATE(device='[freg(' + str(D) + ')%lo(:, :, :, amr_reg_cur)]')
                    end if
                    if (s_hi(${D}$)) then
                        $:GPU_UPDATE(device='[freg(' + str(D) + ')%hi(:, :, :, amr_reg_cur)]')
                    end if
                end if
            #:endfor
        end do
        call s_phase_toc(PH_RFRECV)
#endif

    end subroutine s_amr_reflux_faces_wave

    !> I5-F5b: the split-ownership level>=2 freg exchange as ONE wave, run once before the reflux fold (the registers are final
    !! after the advance, and the applies keep their per-box reverse-order position). Replaces one fully BLOCKING SEND/RECV pair per
    !! dim per split child. Same zero-copy, companion-header, single-tag-by-order design as the faces wave; tag block
    !! amr_tag_base(5)+50 keeps it disjoint from the faces wave's within the family. The subcycle path keeps its per-box exchange
    !! inside s_amr_reflux_to_parent (do_xchg).
    impure subroutine s_amr_freg_wave()

#ifdef MFC_MPI
        use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        integer  :: k, ierr, nreq, cnt, pblk, cowner, powner, tq, nhr, nhs, j
        real(wp) :: w_lo(3), w_hi(3), nanv

        if (num_procs == 1) return
        call s_amr_reg_prepare()
        tq = amr_tag_base(5) + 50 + int(mod(amr_mesh_epoch, 50_8))
        nanv = ieee_value(0._wp, ieee_quiet_nan)
        nreq = 0; nhr = 0; nhs = 0
        do k = 1, amr_num_blocks
            if (amr_block_level(k) < 2) cycle
            call s_amr_select_slot(k)
            pblk = f_amr_parent_block(k)
            cowner = amr_block_owner(k); powner = amr_block_owner(pblk)
            if (cowner == powner .or. powner /= proc_rank) cycle
            ! seam clip: a face weighted 0 by the sibling-seam rule is never consumed by the parent-side reflux apply
            ! (s_amr_reflux_to_parent multiplies it away), so it never ships - both sides derive the identical skip from
            ! s_amr_sibling_face_weights on replicated metadata. Debug arm: skipped-face mirrors are NaN-flooded so any
            ! OTHER consumer of an unshipped face aborts within the step.
            call s_amr_sibling_face_weights(k, pblk, w_lo, w_hi)
            if (XA_NH > 0) then
                nhr = nhr + 1
                call s_amr_fw_szr(amr_fw_rq, XA_NH*nhr); call s_amr_fw_szi(amr_fw_rblk, nhr)
                amr_fw_rblk(nhr) = k
                nreq = nreq + 1
                call s_amr_fw_szi(amr_fw_req, nreq); call s_amr_fw_szi(amr_fw_reqw, nreq)
                amr_fw_reqw(nreq) = XA_NH
                call MPI_IRECV(amr_fw_rq(XA_NH*(nhr - 1) + 1), XA_NH, mpi_p, cowner, tq, MPI_COMM_WORLD, amr_fw_req(nreq), ierr)
            end if
            #:for D in [1, 2, 3]
                if (${D}$ <= num_dims) then
                    cnt = size(freg(${D}$)%lo(:,:,:,amr_reg_cur))
                    if (w_lo(${D}$) > 0._wp) then
                        nreq = nreq + 1
                        call s_amr_fw_szi(amr_fw_req, nreq); call s_amr_fw_szi(amr_fw_reqw, nreq)
                        amr_fw_reqw(nreq) = cnt
                        call s_xa_rec(XA_F5W_FREG_RCV, 2, cnt, tq)
                        call MPI_IRECV(freg(${D}$)%lo(:,:,:,amr_reg_cur), cnt, mpi_p, cowner, tq, MPI_COMM_WORLD, &
                                       & amr_fw_req(nreq), ierr)
#ifdef MFC_DEBUG
                    else
                        freg(${D}$)%lo(:,:,:,amr_reg_cur) = nanv
                        $:GPU_UPDATE(device='[freg(' + str(D) + ')%lo(:, :, :, amr_reg_cur)]')
#endif
                    end if
                    if (w_hi(${D}$) > 0._wp) then
                        nreq = nreq + 1
                        call s_amr_fw_szi(amr_fw_req, nreq); call s_amr_fw_szi(amr_fw_reqw, nreq)
                        amr_fw_reqw(nreq) = cnt
                        call s_xa_rec(XA_F5W_FREG_RCV, 2, cnt, tq)
                        call MPI_IRECV(freg(${D}$)%hi(:,:,:,amr_reg_cur), cnt, mpi_p, cowner, tq, MPI_COMM_WORLD, &
                                       & amr_fw_req(nreq), ierr)
#ifdef MFC_DEBUG
                    else
                        freg(${D}$)%hi(:,:,:,amr_reg_cur) = nanv
                        $:GPU_UPDATE(device='[freg(' + str(D) + ')%hi(:, :, :, amr_reg_cur)]')
#endif
                    end if
                end if
            #:endfor
        end do
        do k = 1, amr_num_blocks
            if (amr_block_level(k) < 2) cycle
            call s_amr_select_slot(k)
            pblk = f_amr_parent_block(k)
            cowner = amr_block_owner(k); powner = amr_block_owner(pblk)
            if (cowner == powner .or. cowner /= proc_rank) cycle
            ! seam clip, send side: the identical weight derivation as the recv walk (replicated metadata), so the
            ! posted sends pair the posted recvs exactly. Skipped faces also skip their device->host pulls.
            call s_amr_sibling_face_weights(k, pblk, w_lo, w_hi)
            #:for D in [1, 2, 3]
                if (${D}$ <= num_dims) then
                    if (w_lo(${D}$) > 0._wp) then
                        $:GPU_UPDATE(host='[freg(' + str(D) + ')%lo(:, :, :, amr_reg_cur)]')
                    end if
                    if (w_hi(${D}$) > 0._wp) then
                        $:GPU_UPDATE(host='[freg(' + str(D) + ')%hi(:, :, :, amr_reg_cur)]')
                    end if
                end if
            #:endfor
            if (XA_NH > 0) then
                nhs = nhs + 1
                call s_amr_fw_szr(amr_fw_sq, XA_NH*nhs)
                call s_xa_hdr_pack(amr_fw_sq(XA_NH*(nhs - 1) + 1:XA_NH*nhs), XA_F5W_FREG_SND, k, [0, 0, 0], [0, 0, 0])
                nreq = nreq + 1
                call s_amr_fw_szi(amr_fw_req, nreq); call s_amr_fw_szi(amr_fw_reqw, nreq)
                amr_fw_reqw(nreq) = -1
                call MPI_ISEND(amr_fw_sq(XA_NH*(nhs - 1) + 1), XA_NH, mpi_p, powner, tq, MPI_COMM_WORLD, amr_fw_req(nreq), ierr)
            end if
            #:for D in [1, 2, 3]
                if (${D}$ <= num_dims) then
                    cnt = size(freg(${D}$)%lo(:,:,:,amr_reg_cur))
                    if (w_lo(${D}$) > 0._wp) then
                        nreq = nreq + 1
                        call s_amr_fw_szi(amr_fw_req, nreq); call s_amr_fw_szi(amr_fw_reqw, nreq)
                        amr_fw_reqw(nreq) = -1
                        call s_xa_rec(XA_F5W_FREG_SND, 1, cnt, tq)
                        call MPI_ISEND(freg(${D}$)%lo(:,:,:,amr_reg_cur), cnt, mpi_p, powner, tq, MPI_COMM_WORLD, &
                                       & amr_fw_req(nreq), ierr)
                    end if
                    if (w_hi(${D}$) > 0._wp) then
                        nreq = nreq + 1
                        call s_amr_fw_szi(amr_fw_req, nreq); call s_amr_fw_szi(amr_fw_reqw, nreq)
                        amr_fw_reqw(nreq) = -1
                        call s_xa_rec(XA_F5W_FREG_SND, 1, cnt, tq)
                        call MPI_ISEND(freg(${D}$)%hi(:,:,:,amr_reg_cur), cnt, mpi_p, powner, tq, MPI_COMM_WORLD, &
                                       & amr_fw_req(nreq), ierr)
                    end if
                end if
            #:endfor
        end do
        if (nreq > 0) then
#ifdef MFC_DEBUG
            block
                integer :: st(MPI_STATUS_SIZE, nreq), gotw, q
                call MPI_WAITALL(nreq, amr_fw_req, st, ierr)
                do q = 1, nreq
                    if (amr_fw_reqw(q) < 0) cycle
                    call MPI_GET_COUNT(st(:,q), mpi_p, gotw, ierr)
                    @:ASSERT(gotw == amr_fw_reqw(q), "freg wave: received message length differs from the plan")
                end do
            end block
#else
            call MPI_WAITALL(nreq, amr_fw_req, MPI_STATUSES_IGNORE, ierr)
#endif
        end if
        j = 0
        do k = 1, amr_num_blocks
            if (amr_block_level(k) < 2) cycle
            call s_amr_select_slot(k)
            pblk = f_amr_parent_block(k)
            cowner = amr_block_owner(k); powner = amr_block_owner(pblk)
            if (cowner == powner .or. powner /= proc_rank) cycle
            if (XA_NH > 0) then
                j = j + 1
                @:ASSERT(amr_fw_rblk(j) == k, "freg wave: header slot order broke")
                call s_xa_hdr_check(amr_fw_rq(XA_NH*(j - 1) + 1:XA_NH*j), XA_F5W_FREG_SND, k, [0, 0, 0], [0, 0, 0])
            end if
            ! push only the faces that shipped; a skipped face keeps its device content (dead under weight 0 - and
            ! NaN-poisoned in debug, so any other reader aborts)
            call s_amr_sibling_face_weights(k, pblk, w_lo, w_hi)
            #:for D in [1, 2, 3]
                if (${D}$ <= num_dims) then
                    if (w_lo(${D}$) > 0._wp) then
                        $:GPU_UPDATE(device='[freg(' + str(D) + ')%lo(:, :, :, amr_reg_cur)]')
                    end if
                    if (w_hi(${D}$) > 0._wp) then
                        $:GPU_UPDATE(device='[freg(' + str(D) + ')%hi(:, :, :, amr_reg_cur)]')
                    end if
                end if
            #:endfor
        end do
#endif

    end subroutine s_amr_freg_wave

    !> True iff this rank is the AUTHORITATIVE holder of global coarse cell g in one dimension (o = interior origin start_idx, ext =
    !! interior extent m/n/p, glb = global last index). A cell is owned by exactly one rank: its interior owner, or - for a
    !! physical-exterior ghost (g < 0 or g > glb) - the boundary-adjacent rank that holds it as a ghost. Inter-rank ghosts are
    !! deliberately NOT claimed (the neighbour's interior owns them), so the sentinel-MAX gather has no double-contribution and
    !! needs no coarse-ghost halo exchange for correctness across rank seams.
    pure logical function f_amr_own_coarse(g, o, ext, glb) result(mine)

        integer, intent(in) :: g, o, ext, glb
        ! interior left physical ghost (leftmost rank)

        mine = (g >= o .and. g <= o + ext) .or. (g < 0 .and. o == 0 .and. g >= -buff_size) .or. (g > glb .and. o + ext == glb &
                & .and. g <= glb + buff_size)  ! right physical ghost (rightmost rank)

    end function f_amr_own_coarse

    !> Per-block measured-cost weight over each block's LEVEL-0 footprint, replicated on every rank: each rank sums the load-weight
    !! cost model (base 1 + K_ib per IB-marked cell + K_pc per phase-change Newton iteration, when that diagnostic array is live)
    !! over its owned coarse cells inside the footprint, then one MPI_ALLREDUCE(SUM) makes the vector identical everywhere. No cost
    !! signals -> cost(k) = footprint cell count exactly (pure-geometry fallback). The Lagrangian cloud is excluded from blocks by
    !! construction, so K_bub never applies here; pc_iter_count is populated only when a load-weight diagnostic writer is on (enable
    !! load_weight_wrt to make the balance phase-change-aware) - guarded by allocated().
    impure subroutine s_amr_block_cost(cost)

        real(wp), intent(out) :: cost(:)
        integer               :: k, j, kk, l, lo(3), hi(3)
        real(wp)              :: c

#ifdef MFC_MPI
        integer :: ierr
#endif

        ! one host refresh per regrid: both signal fields advance on the device between regrids
        if (ib) then
            $:GPU_UPDATE(host='[ib_markers%sf]')
        end if
        if (allocated(pc_iter_count)) then
            $:GPU_UPDATE(host='[pc_iter_count]')
        end if
        do k = 1, amr_num_blocks
            ! block footprint /\ this rank's coarse subdomain, in local interior indices (empty -> no-trip loops)
            lo = 0; hi = 0
            lo(1) = max(amr_region_lo_all(1, k) - start_idx(1), 0); hi(1) = min(amr_region_hi_all(1, k) - start_idx(1), m)
            if (n_glb > 0) then
                lo(2) = max(amr_region_lo_all(2, k) - start_idx(2), 0); hi(2) = min(amr_region_hi_all(2, k) - start_idx(2), n)
            end if
            if (p_glb > 0) then
                lo(3) = max(amr_region_lo_all(3, k) - start_idx(3), 0); hi(3) = min(amr_region_hi_all(3, k) - start_idx(3), p)
            end if
            c = 0._wp
            do l = lo(3), hi(3)
                do kk = lo(2), hi(2)
                    do j = lo(1), hi(1)
                        c = c + 1._wp
                        if (ib) then
                            if (ib_markers%sf(j, kk, l) /= 0) c = c + K_ib
                        end if
                        if (allocated(pc_iter_count)) c = c + K_pc*real(pc_iter_count(j, kk, l), wp)
                    end do
                end do
            end do
            cost(k) = c
        end do
#ifdef MFC_MPI
        amr_gb_cost = amr_gb_cost + int(amr_num_blocks, 8)*8_8
        call MPI_ALLREDUCE(MPI_IN_PLACE, cost, amr_num_blocks, mpi_p, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

    end subroutine s_amr_block_cost

    !> Cost-weighted SFC partition of n items into num_procs contiguous Morton-key ranges. Sorts item indices by ascending Morton
    !! key (insertion sort - n is small), walks them accumulating wt, and advances the owner rank when the running weight crosses
    !! the next even share of the total. Writes owner(1:n) (rank per item) and cut(0:num_procs-1) (running Morton upper bound per
    !! rank; ranks that receive no item inherit the predecessor's bound so cut is non-decreasing, its top the global max key).
    !! Shared by the fine-block anchor split and the L0-tile split so the two owner maps cannot drift.
    subroutine s_amr_sfc_cut(keys, wt, n, cut, owner)

        integer, intent(in)          :: n
        integer(kind=8), intent(in)  :: keys(n)
        real(wp), intent(in)         :: wt(n)
        integer(kind=8), intent(out) :: cut(0:num_procs - 1)
        integer, intent(out)         :: owner(n)
        integer                      :: ord(n), mrg(n), k, r
        integer                      :: width, lo_m, mid_m, hi_m, i_m, j_m, t_m
        real(wp)                     :: total, cum, tgt, tol

        cut = -1_8
        if (n < 1) return

        ! Sort item indices by Morton key. Bottom-up MERGE sort: O(n log n), STABLE (ties keep their original
        ! order), iterative, and a pure function of the input - every rank must produce byte-identical
        ! order or the assignment diverges and s_amr_validate_owner aborts. The previous insertion sort
        ! was O(n^2) with the comment "n small"; that holds at the box counts benchmarked so far (512
        ! boxes = 1.3e5 ops) and fails at the counts this design targets (1e5 boxes = 5e9 ops, seconds
        ! per regrid), since the whole strategy is boxes_per_level >> num_procs.
        do k = 1, n
            ord(k) = k
        end do
        width = 1
        do while (width < n)
            lo_m = 1
            do while (lo_m <= n - width)
                mid_m = lo_m + width - 1
                hi_m = min(lo_m + 2*width - 1, n)
                i_m = lo_m; j_m = mid_m + 1; t_m = lo_m
                do while (i_m <= mid_m .and. j_m <= hi_m)
                    ! <= keeps the left run first on ties: stability
                    if (keys(ord(i_m)) <= keys(ord(j_m))) then
                        mrg(t_m) = ord(i_m); i_m = i_m + 1
                    else
                        mrg(t_m) = ord(j_m); j_m = j_m + 1
                    end if
                    t_m = t_m + 1
                end do
                do while (i_m <= mid_m); mrg(t_m) = ord(i_m); i_m = i_m + 1; t_m = t_m + 1; end do
                do while (j_m <= hi_m); mrg(t_m) = ord(j_m); j_m = j_m + 1; t_m = t_m + 1; end do
                ord(lo_m:hi_m) = mrg(lo_m:hi_m)
                lo_m = lo_m + 2*width
            end do
            width = 2*width
        end do

        total = 0._wp
        do k = 1, n
            total = total + wt(k)
        end do

        ! chains-on-chains over the items in SFC order; advance the owner rank when the cumulative weight crosses the next even
        ! share.
        ! All-real arithmetic on replicated weights in a fixed order, so every rank computes the identical assignment.
        r = 0; cum = 0._wp
        do k = 1, n
            tgt = real(r + 1, wp)*total/real(num_procs, wp)
            ! cum is an n-term ACCUMULATION while tgt is CLOSED FORM over another n-term sum, so at an exact share boundary the
            ! two differ by rounding rather than by intent, and the comparison turns on 1 ULP. This was correct only by luck:
            ! every cost term to date is integer-valued (footprint cells, K_ib, K_pc x integer counts), so the arithmetic was
            ! exact. MEASURED with a fractional cost term: 32 IDENTICAL weights over 8 ranks split 5/3 instead of 4/4, reporting
            ! max/mean 1.250 where the same case with integer weights reports exactly 1.000. Tolerance is the accumulated
            ! rounding bound, O(n) ULP of the target; far from a boundary it is negligible and the greedy is unchanged.
            tol = spacing(tgt)*real(n, wp)
            if (cum >= tgt - tol .and. r < num_procs - 1) r = r + 1
            owner(ord(k)) = r
            cut(r) = keys(ord(k))  ! items visited in ascending Morton key => running upper bound for rank r
            cum = cum + wt(ord(k))
        end do
        ! ranks that received no item (or trail the last-assigned rank) inherit the predecessor's bound so the search never lands on
        ! them: cut is non-decreasing and its top equals the global max key.
        do r = 1, num_procs - 1
            if (cut(r) < cut(r - 1)) cut(r) = cut(r - 1)
        end do

    end subroutine s_amr_sfc_cut

    !> Fine-level distribution map: assigns each active block a single owner rank by chains-on-chains balancing of fine-work weight
    !! in Morton order of the block's low corner - the same SFC idea m_sfc_partition uses for the base grid, at block granularity.
    !! Fine work = measured coarse-footprint cost (s_amr_block_cost) x the level's refinement factor per active dim, so blocks
    !! concentrating IB or phase-change work weigh more than equal-size quiescent ones. The cost vector is allreduced (one
    !! collective; every rank must call this), after which the assignment is deterministic and identical on every rank.
    !! s_set_amr_fine_geometry applies it as amr_rank_owns_block = (amr_block_owner(amr_cur) == proc_rank).
    impure subroutine s_amr_assign_block_owners()

        integer :: k, a, lev, maxlev, na
        ! heap, not stack: these are O(global boxes) and at 1e6 blocks the seven together are ~48 MB,
        ! which overflows a default stack long before the box count itself becomes a problem
        integer, allocatable         :: aidx(:), aown(:)
        integer(kind=8), allocatable :: key(:), akey(:)
        real(wp), allocatable        :: wt(:), cost(:), awt(:)

        if (amr_num_blocks < 1) return

        allocate (aidx(amr_num_blocks), aown(amr_num_blocks), key(amr_num_blocks), akey(amr_num_blocks), wt(amr_num_blocks), &
                  & cost(amr_num_blocks), awt(amr_num_blocks))

        call s_amr_block_cost(cost)

        ! per-block own fine-work weight = footprint cost x amr_ref_ratio**(level*active dims). A level-l block is amr_ref_ratio**l
        ! finer than L0 per dim, so its work is the footprint cost x rr**(l*d). The level factor now only scales blocks WITHIN a
        ! level relative to each other (each level is cut separately), but it stays because a level's boxes can differ in footprint.
        ! With no cost signals this reduces to the fine cell count (geometry only).
        do k = 1, amr_num_blocks
            wt(k) = cost(k)*real(amr_ref_ratio, wp)**amr_block_level(k)
            if (n_glb > 0) wt(k) = wt(k)*real(amr_ref_ratio, wp)**amr_block_level(k)
            if (p_glb > 0) wt(k) = wt(k)*real(amr_ref_ratio, wp)**amr_block_level(k)
            key(k) = f_morton(amr_region_lo_all(1, k), amr_region_lo_all(2, k), amr_region_lo_all(3, k))
        end do

        ! PER-LEVEL DISTRIBUTION: balance EVERY level independently, each block on its OWN weight. Tower co-location is gone - a
        ! level-1 block and its descendants are assigned separately, so a deep tower no longer pins its whole subtree (weight
        ! cost*rr**(l*d)) to one rank, which is what capped granularity at depth. The parent<->child gather/restrict/reflux paths
        ! are P2P, so a split tower costs messages rather than correctness.
        !
        ! One cut PER LEVEL, not one mixed cut over all fine blocks: same-level boxes are disjoint and so have DISTINCT Morton
        ! keys, which the cut-point binary search in f_amr_owner needs. Mixed, a level-2 block sharing its parent's region_lo would
        ! collide with it and the search could not tell them apart.
        !
        ! Each level's cut goes into amr_fine_cut(:, lev): fine blocks straddle tiles, so their owner is not tile-cut-derivable and
        ! f_amr_owner reads amr_fine_cut for them. amr_owner_cut mirrors LEVEL 1 only without tiles, where the two are the same
        ! authority. Under coexist amr_owner_cut holds the TILE cut that s_l0_tiles_init built; overwriting it here is harmless at
        ! init (the assigner runs first) but at REGRID time would clobber the tile cut.
        ! Fine blocks occupy slots (l0_slot_off, amr_num_blocks]; slots [1, l0_slot_off] are the L0 TILE
        ! prefix. At init the assigner runs BEFORE s_l0_tiles_init, so those prefix slots are still
        ! uninitialized - level reads 1 and region_lo is all zeros, i.e. Morton key 0. Including them
        ! fed 8 phantom key-0 "level-1 blocks" into the level-1 cut, which the cut then split across
        ! ranks. A key-0 block can only ever resolve to rank 0 (cut is non-decreasing and the search
        ! returns the first r with key <= cut(r)), so any phantom placed on a higher rank is
        ! unrecoverable and s_amr_validate_owner aborts. This was latent: with the old weights all the
        ! phantoms happened to land on rank 0 and the validator agreed by luck.
        maxlev = maxval(amr_block_level(l0_slot_off + 1:amr_num_blocks))
        do lev = 1, maxlev
            na = 0
            do k = l0_slot_off + 1, amr_num_blocks
                if (amr_block_level(k) /= lev) cycle
                na = na + 1
                akey(na) = key(k); awt(na) = wt(k); aidx(na) = k
            end do
            if (na < 1) cycle
            call s_amr_sfc_cut(akey, awt, na, amr_fine_cut(:,lev), aown)
            do a = 1, na
                amr_block_owner(aidx(a)) = aown(a)
            end do
            if (lev == 1 .and. l0_slot_off == 0) amr_owner_cut = amr_fine_cut(:,1)
        end do

        call s_amr_validate_owner()
        call s_amr_report_balance(wt, maxlev)

        deallocate (aidx, aown, key, akey, wt, cost, awt)

    end subroutine s_amr_assign_block_owners

    !> Per-level and total load-balance report: max/mean assigned block weight over ranks, the metric the balancer is actually
    !! trying to minimise. Without it a distribution change can only be judged by end-to-end s/step, which cannot separate
    !! "balanced" from "uniformly slow" - the step-4 A/B measured flat and was uninterpretable for exactly that reason.
    !!
    !! Needs no MPI: wt, amr_block_level and amr_block_owner are replicated and identical on every rank (the cost vector is
    !! allreduced in s_amr_block_cost), so every rank computes the same numbers and rank 0 prints. ratio == 1 is perfect balance;
    !! ratio == num_procs means one rank holds everything at that level. no_blocks_ranks counts ranks holding no block AT THIS
    !! LEVEL, which is the granularity floor showing up directly: a level with fewer boxes than ranks CANNOT balance, however good
    !! the cut is. It is NOT an idleness measure - those ranks still own level-0 work (level 0 covers every rank) and may own
    !! blocks at other levels. Only m_rank_timing measures idleness; do not read this counter as one.
    impure subroutine s_amr_report_balance(wt, maxlev)

        real(wp), intent(in)  :: wt(:)
        integer, intent(in)   :: maxlev
        real(wp), allocatable :: rw(:), tw(:), rc(:)
        real(wp)              :: mx, mean, cmx, cmean
        integer               :: k, lev, nb, empty

        if (.not. load_weight_wrt) return
        if (proc_rank /= 0) return

        ! heap, not automatic: these are num_procs long and this is the routine that runs AT SCALE - two automatic wp arrays would
        ! put O(num_procs) on the stack, which is where the module already puts amr_owner_cut / amr_fine_cut on the heap instead
        allocate (rw(0:num_procs - 1), tw(0:num_procs - 1), rc(0:num_procs - 1))
        tw = 0._wp
        do lev = 1, maxlev
            rw = 0._wp; rc = 0._wp
            nb = 0
            do k = 1, amr_num_blocks
                if (amr_block_level(k) /= lev) cycle
                rw(amr_block_owner(k)) = rw(amr_block_owner(k)) + wt(k)
                rc(amr_block_owner(k)) = rc(amr_block_owner(k)) + 1._wp
                nb = nb + 1
            end do
            if (nb == 0) cycle
            tw = tw + rw
            mx = maxval(rw); mean = sum(rw)/real(num_procs, wp)
            empty = count(rw <= 0._wp)
            ! BOX COUNT imbalance beside weight imbalance. cost(k) is a footprint CELL count, but the measured per-block advance
            ! costs ~1x a full monolithic step almost regardless of block size and does not amortize - so if per-block cost is
            ! largely FIXED, a rank's true load tracks how many boxes it holds, not how many cells. Equal cells with unequal box
            ! counts then reads as perfectly balanced and runs badly skewed, which is exactly the model/measured gap observed
            ! (model 1.03 vs rank_time 1.26 at np=8). Printing both makes that hypothesis falsifiable in one run.
            cmx = maxval(rc); cmean = sum(rc)/real(num_procs, wp)
            ! NOT merge(): merge is a function, so BOTH arms are evaluated and the mean == 0 arm would still divide by zero - the
            ! guard would not guard. no_blocks_ranks counts ranks holding no block AT THIS LEVEL; they are not idle (they still
            ! own level-0 work and possibly other levels), they just take no share of this level's.
            ! cmean cannot be zero here: nb >= 1, so sum(rc) = nb >= 1. No guard needed, and emphatically not merge().
            if (mean > 0._wp) print '(A,I0,A,I0,A,I0,A,F8.3,A,F8.3,A,I0,A,I0)', ' [amr-balance] level ', lev, ': boxes ', nb, &
                & '/ranks ', num_procs, ' max/mean ', mx/mean, ' boxes_max/mean ', cmx/cmean, ' no_blocks_ranks ', empty, ' of ', &
                & num_procs
        end do
        ! per-rank WEIGHT (fine cells) and BOX COUNT, so rhs time can be regressed against actual load
        if (proc_rank == 0) then
            write (*, '(A)', advance='no') ' [amr-balance] per-rank fine_work :'
            do k = 0, num_procs - 1; write (*, '(I12)', advance='no') nint(tw(k), kind=8); end do
            write (*, '(A)') ''
        end if
        mean = sum(tw)/real(num_procs, wp)
        ! fine_work = sum of the assigned weights = the FINE CELLS advanced per step (with no cost signals wt is exactly that).
        ! Without it an AMR-vs-uniform wall-clock ratio is uninterpretable: 5x the time is the expected outcome of advancing 5x
        ! the cells and the failure mode of 5x per-cell overhead, and the two are indistinguishable from the ratio alone.
        if (mean > 0._wp) print '(A,F8.3,A,I0,A,I0)', ' [amr-balance] TOTAL   : max/mean ', maxval(tw)/mean, &
            & ' ranks_with_no_fine_block ', count(tw <= 0._wp), ' fine_work ', nint(sum(tw), kind=8)
        deallocate (rw, tw, rc)

    end subroutine s_amr_report_balance

    !> Owner rank of block k from the O(num_procs) SFC cut-points: binary-search k's OWN Morton key in the owning authority's cut. A
    !! level-0 TILE resolves against amr_owner_cut (the tile cut in tiled modes); a FINE block (level>=1) resolves against its own
    !! level's cut amr_fine_cut(:, level) (level 1's == amr_owner_cut in no-tile AMR). No tower-anchor resolution: under per-level
    !! distribution a block's owner depends on its own key and level alone. Reproduces s_amr_assign_block_owners' / the tile split's
    !! cost-weighted SFC assignment exactly.
    pure integer function f_amr_owner(k) result(r)

        integer, intent(in) :: k
        integer(kind=8)     :: mk, cut(0:num_procs - 1)
        integer             :: a, lo, hi, mid

        a = k
        if (amr_block_level(k) == 0) then
            cut = amr_owner_cut  ! tile: own Morton key vs the tile cut
        else
            cut = amr_fine_cut(:,amr_block_level(k))
        end if
        mk = f_morton(amr_region_lo_all(1, a), amr_region_lo_all(2, a), amr_region_lo_all(3, a))
        lo = 0; hi = num_procs - 1
        do while (lo < hi)
            mid = (lo + hi)/2
            if (mk <= cut(mid)) then
                hi = mid
            else
                lo = mid + 1
            end if
        end do
        r = lo

    end function f_amr_owner

    !> TRANSITIONAL: the SFC cut-point accessor must reproduce the stored owner table exactly (removed once the table is deleted).
    impure subroutine s_amr_validate_owner()

        integer :: k

        do k = 1, amr_num_blocks
            ! Skip the L0 tile prefix when its cut has not been built yet (amr_owner_cut still -1): those
            ! slots are uninitialized at assigner time and carry a stale level with Morton key 0. The
            ! tile-init call site populates both cuts and validates them there.
            if (k <= l0_slot_off .and. amr_owner_cut(num_procs - 1) < 0_8) cycle
            ! every block resolves: tiles (level 0) via amr_owner_cut (tile cut), fine blocks (level>=1) via amr_fine_cut. The
            ! caller
            ! guarantees the relevant cut is populated for the blocks present at each call site (assigner: fine cut; tile init:
            ! both).
            if (f_amr_owner(k) /= amr_block_owner(k)) &
                & call s_mpi_abort('SFC cut-point owner disagrees with amr_block_owner - cut capture or search is wrong')
        end do

    end subroutine s_amr_validate_owner

    impure subroutine s_amr_compute_isect(lo, hi)

        integer, intent(in) :: lo(3), hi(3)
        integer             :: sidx(3), ext(3), d

        sidx = 0; ext = 0
        sidx(1) = start_idx(1); ext(1) = m
        if (n_glb > 0) then; sidx(2) = start_idx(2); ext(2) = n; end if
        if (p_glb > 0) then; sidx(3) = start_idx(3); ext(3) = p; end if
        do d = 1, 3
            amr_isect_lo(d) = max(lo(d), sidx(d))
            amr_isect_hi(d) = min(hi(d), sidx(d) + ext(d))
        end do
        amr_rank_owns_block = amr_isect_lo(1) <= amr_isect_hi(1)
        if (n_glb > 0) amr_rank_owns_block = amr_rank_owns_block .and. amr_isect_lo(2) <= amr_isect_hi(2)
        if (p_glb > 0) amr_rank_owns_block = amr_rank_owns_block .and. amr_isect_lo(3) <= amr_isect_hi(3)

    end subroutine s_amr_compute_isect

    !> Set the fine level's geometry (region, intersection, extents, bounds, coordinates) for the box lo:hi. Arrays are preallocated
    !! at max size; this only updates metadata and refills coords. Collective: ALL ranks must call together (init and regrid do) -
    !! it also refreshes the allreduced amr_xchg_coarse_ghosts flag for the new box. INVARIANT: a level-l block's fine extent is
    !! amr_ref_ratio**l * (coarse-region width) - 1, NOT amr_ref_ratio*width. (amr_ref_ratio*width holds only for the level-1
    !! initial block; nested boxes compound by amr_ref_ratio per level.) Every fine-extent computation - here, the restart-reader
    !! check, the load-weight, the fmul - uses amr_ref_ratio**level; assuming amr_ref_ratio*width rejects level>=2 blocks as corrupt
    !! (the exact bug that bit the multi-level restart reader).
    impure subroutine s_set_amr_fine_geometry(lo, hi)

        integer, intent(in) :: lo(3), hi(3)
        integer             :: sidx(3), ext(3), nmar, bad_loc, pblk, d, rr

        amr_slots(amr_cur)%region%lo = lo; amr_slots(amr_cur)%region%hi = hi
        amr_region_lo = lo; amr_region_hi = hi  ! global mirror for m_amr_registers (no use-cycle)
        amr_region_lo_all(:,amr_cur) = lo; amr_region_hi_all(:,amr_cur) = hi

        ! FINE-LEVEL DISTRIBUTION: a block is owned WHOLE by amr_block_owner(k). The owner holds fine cells for the ENTIRE block;
        ! every other rank holds none. amr_isect_lo/hi records the block's coarse footprint (= the whole block on the owner) - it
        ! drives the coarse<->fine gather/scatter (which coarse cells the owner pulls in / pushes back). At np=1 the owner is rank 0
        ! and the footprint is the whole domain-resident block, so this reduces exactly to the old mirror (block \cap subdomain ==
        ! whole block).
        amr_rank_owns_block = (amr_block_owner(amr_cur) == proc_rank)
        pblk = 0
        if (amr_rank_owns_block) then
            amr_isect_lo = lo; amr_isect_hi = hi
            if (amr_block_level(amr_cur) >= 2) then
                ! multi-level: express the coarse footprint in the PARENT block's fine-cell frame (a level-l block's coarse side is
                ! level l-1). parent-fine index of L0 cell c is rr*(c - R1.lo) where rr is the parent's amr_ref_ratio; the block
                ! spans rr fine cells per parent-covered L0 cell. m below then gets amr_ref_ratio*(footprint) cells, as for a
                ! level-1
                ! block over L0. amr_cg / the prolong read this frame, so no other coupling code changes for the local (np=1) path.
                ! rr is the GLOBAL amr_ref_ratio, not amr_slots(pblk)%amr_ref_ratio: that field is written by s_amr_alloc_slot,
                ! which a rank owning this block but NOT its parent never calls for pblk, so it would read undefined. The two agree
                ! wherever both are defined - only an L0 tile carries a per-slot ratio of 1, and a level>=2 block's parent is never
                ! an L0 tile. This is the same footprint s_amr_parent_foot derives from replicated metadata.
                pblk = f_amr_parent_block(amr_cur)
                call s_amr_parent_foot(amr_cur, pblk, amr_isect_lo, amr_isect_hi)
            end if
        else
            amr_isect_lo = 1; amr_isect_hi = 0  ! empty footprint
            if (n_glb > 0) then; amr_isect_lo(2) = 1; amr_isect_hi(2) = 0; end if
            if (p_glb > 0) then; amr_isect_lo(3) = 1; amr_isect_hi(3) = 0; end if
        end if
        amr_isect_lo_all(:,amr_cur) = amr_isect_lo; amr_isect_hi_all(:,amr_cur) = amr_isect_hi
        amr_owns_all(amr_cur) = amr_rank_owns_block
        ! fine extents cover the WHOLE block on the owner; -1 (empty) on non-owners
        amr_slots(amr_cur)%m = amr_ref_ratio*max(amr_isect_hi(1) - amr_isect_lo(1) + 1, 0) - 1
        amr_slots(amr_cur)%n = 0; amr_slots(amr_cur)%p = 0
        if (n_glb > 0) amr_slots(amr_cur)%n = amr_ref_ratio*max(amr_isect_hi(2) - amr_isect_lo(2) + 1, 0) - 1
        if (p_glb > 0) amr_slots(amr_cur)%p = amr_ref_ratio*max(amr_isect_hi(3) - amr_isect_lo(3) + 1, 0) - 1
        amr_slots(amr_cur)%idwbuff(1)%beg = -buff_size; amr_slots(amr_cur)%idwbuff(1)%end = amr_slots(amr_cur)%m + buff_size
        amr_slots(amr_cur)%idwbuff(2)%beg = 0; amr_slots(amr_cur)%idwbuff(2)%end = 0
        amr_slots(amr_cur)%idwbuff(3)%beg = 0; amr_slots(amr_cur)%idwbuff(3)%end = 0
        if (n_glb > 0) then
            amr_slots(amr_cur)%idwbuff(2)%beg = -buff_size; amr_slots(amr_cur)%idwbuff(2)%end = amr_slots(amr_cur)%n + buff_size
        end if
        if (p_glb > 0) then
            amr_slots(amr_cur)%idwbuff(3)%beg = -buff_size; amr_slots(amr_cur)%idwbuff(3)%end = amr_slots(amr_cur)%p + buff_size
        end if
        ! coord building only on ranks with fine cells (others never read their coord arrays)
        if (amr_rank_owns_block) then
            ! Every level builds the same way: replay the ancestor chain from the GLOBAL L0 boundaries. The owner may hold no part
            ! of the coarse slice it refines, and (level>=2) may not own the parent at all, so neither the local coarse coords nor
            ! the parent's slot can be read here. At level 1 the chain is one step and this is the previous global-boundary form.
            call s_amr_build_block_coords(amr_cur, amr_gxcb, amr_slots(amr_cur)%x_cb, amr_slots(amr_cur)%x_cc, &
                                          & amr_slots(amr_cur)%dx, 1)
            if (n_glb > 0) call s_amr_build_block_coords(amr_cur, amr_gycb, amr_slots(amr_cur)%y_cb, amr_slots(amr_cur)%y_cc, &
                & amr_slots(amr_cur)%dy, 2)
            if (p_glb > 0) call s_amr_build_block_coords(amr_cur, amr_gzcb, amr_slots(amr_cur)%z_cb, amr_slots(amr_cur)%z_cc, &
                & amr_slots(amr_cur)%dz, 3)
        end if

        ! Fine ghost prolongation reads up to nmar coarse cells past each face of the intersection; if that stencil leaves ANY
        ! rank's
        ! interior (block near/at/across a rank boundary), the coarse CONS ghosts it reads must be halo-exchanged before every fill
        ! (the solver populates only PRIM ghosts). All ranks agree on the flag, so the pairwise exchanges are called consistently.
        sidx = 0; ext = 0
        sidx(1) = start_idx(1); ext(1) = m
        if (n_glb > 0) then; sidx(2) = start_idx(2); ext(2) = n; end if
        if (p_glb > 0) then; sidx(3) = start_idx(3); ext(3) = p; end if
        nmar = (buff_size + amr_ref_ratio - 1)/amr_ref_ratio + 1
        bad_loc = 0
        if (amr_rank_owns_block) then
            if (amr_isect_lo(1) - sidx(1) < nmar .or. sidx(1) + ext(1) - amr_isect_hi(1) < nmar) bad_loc = 1
            if (n_glb > 0 .and. (amr_isect_lo(2) - sidx(2) < nmar .or. sidx(2) + ext(2) - amr_isect_hi(2) < nmar)) bad_loc = 1
            if (p_glb > 0 .and. (amr_isect_lo(3) - sidx(3) < nmar .or. sidx(3) + ext(3) - amr_isect_hi(3) < nmar)) bad_loc = 1
        end if
        ! ACCUMULATE, do not reduce: the caller closes the scan with s_amr_reduce_xchg_flag. Every caller loops over blocks and
        ! wants "does ANY block need the exchange", so a per-block collective was both O(nboxes) collectives and, at two call
        ! sites, WRONG - those loops kept the LAST block's answer rather than the OR, so an earlier block needing the exchange
        ! could be masked by a later one that did not.
        amr_xchg_bad = max(amr_xchg_bad, bad_loc)

    end subroutine s_set_amr_fine_geometry

    !> Close a geometry scan: ONE allreduce of the accumulated flag, then reset so the next scan starts clean. Must be called after
    !! every s_set_amr_fine_geometry loop (or single call) - the flag it sets is read by the fine advance
    !! (s_amr_exchange_coarse_cons_halo).
    impure subroutine s_amr_reduce_xchg_flag()

        integer :: bad_glb

        call s_mpi_allreduce_integer_max(amr_xchg_bad, bad_glb)
        amr_xchg_coarse_ghosts = bad_glb == 1
        amr_xchg_bad = 0

    end subroutine s_amr_reduce_xchg_flag

    !> Conservative-linear prolongation for a single variable pair. Reads coarse interior/ghost from qc; writes fine interior to qf.
    !! Minmod-limited slopes.
    impure subroutine s_prolong_one_var(qc, loc, ivar, pos, inject)

        type(scalar_field), intent(in) :: qc
        integer, intent(in)            :: loc, ivar  !< flat-store slot and variable of the fine target
        logical, optional, intent(in)  :: pos        !< floor the child at bub_pos_frac*u0 (bubble radius-moment realizability)
        logical, optional, intent(in)  :: inject     !< piecewise-constant (child = u0): QBMM moment realizability preservation
        integer                        :: fi, fj, fk, ci, cj, ck, ox, oy, oz, rrat, mm, nn, pp, il1, il2, il3
        real(wp)                       :: u0, sx, sy, sz, xix, xiy, xiz, child, bpf
        logical                        :: floor_pos, pw_const, d2, d3

        floor_pos = .false.; if (present(pos)) floor_pos = pos
        pw_const = .false.; if (present(inject)) pw_const = inject

        ! coarse source qc is the gathered block-local patch amr_cg (fine-level distribution): amr_isect_lo is GLOBAL and equals
        ! region_lo on the owner, so amr_isect_lo + f/rr - amr_cpat_off = nmar + f/rr is the patch-local coarse index.
        ! DEVICE kernel (device-side rebuild): reads the patch's device mirror (pushed once per prolong dispatch by
        ! s_interpolate_coarse_to_fine), writes the fine slot in place - the per-slot full push at every caller is deleted.
        ! CPU builds compile this to the identical plain loop, so CPU results are unchanged.
        ox = amr_cpat_off(1); oy = amr_cpat_off(2); oz = amr_cpat_off(3)
        rrat = amr_slots(amr_cur)%amr_ref_ratio
        mm = amr_slots(amr_cur)%m; nn = amr_slots(amr_cur)%n; pp = amr_slots(amr_cur)%p
        il1 = amr_isect_lo(1); il2 = amr_isect_lo(2); il3 = amr_isect_lo(3)
        d2 = n_glb > 0; d3 = p_glb > 0
        bpf = bub_pos_frac
        $:GPU_PARALLEL_LOOP(collapse=3, private='[ci, cj, ck, xix, xiy, xiz, u0, sx, sy, sz, child]')
        do fk = 0, pp
            do fj = 0, nn
                do fi = 0, mm
                    ck = il3 + fk/rrat - oz; if (.not. d3) ck = 0
                    xiz = 0._wp; if (d3) xiz = (real(mod(fk, rrat), wp) - real(rrat - 1, wp)*0.5_wp)/real(rrat, wp)
                    cj = il2 + fj/rrat - oy; if (.not. d2) cj = 0
                    xiy = 0._wp; if (d2) xiy = (real(mod(fj, rrat), wp) - real(rrat - 1, wp)*0.5_wp)/real(rrat, wp)
                    ci = il1 + fi/rrat - ox
                    xix = (real(mod(fi, rrat), wp) - real(rrat - 1, wp)*0.5_wp)/real(rrat, wp)
                    u0 = real(qc%sf(ci, cj, ck), wp)
                    sx = minmod(real(qc%sf(ci + 1, cj, ck), wp) - u0, u0 - real(qc%sf(ci - 1, cj, ck), wp))
                    sy = 0._wp
                    if (d2) sy = minmod(real(qc%sf(ci, cj + 1, ck), wp) - u0, u0 - real(qc%sf(ci, cj - 1, ck), wp))
                    sz = 0._wp
                    if (d3) sz = minmod(real(qc%sf(ci, cj, ck + 1), wp) - u0, u0 - real(qc%sf(ci, cj, ck - 1), wp))
                    if (pw_const) then
                        sx = 0._wp; sy = 0._wp; sz = 0._wp
                    end if
                    child = u0 + sx*xix + sy*xiy + sz*xiz
                    if (floor_pos) child = max(child, bpf*u0)
                    amr_cons_st(fi, fj, fk, ivar, loc) = child
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_prolong_one_var

    !> Conservative-linear prolongation: fill amr_fine interior from coarse (level-0), minmod-limited. Symmetric child offsets
    !! (+/-1/4 of a coarse cell) => the amr_ref_ratio^d children average to the coarse value. Multi-fluid volume fractions take the
    !! sum-preserving closure path instead (single-fluid runs never branch, so their prolongation is untouched). TWIN
    !! s_amr_prolong_pbmv (q<->pb/mv): pb/mv sibling of this prolongation (piecewise-constant there); keep the child-offset frame
    !! and volume-fraction closure lockstep.
    impure subroutine s_interpolate_coarse_to_fine()

        integer :: i, bstride

        ! the prolong kernels read the gathered patch's DEVICE mirror; the level-1 patch is host-filled by the gather
        ! unpack, so push it once per dispatch (patch-sized - 8x smaller than the full-slot pushes this replaces; for a
        ! level>=2 block the patch was device-produced and this re-push of the pulled bytes is redundant but harmless)

        do i = 1, sys_size
            $:GPU_UPDATE(device='[amr_cg(i)%sf]')
        end do
        bstride = 1
        if (bubbles_euler) bstride = (eqn_idx%bub%end - eqn_idx%bub%beg + 1)/nb
        do i = 1, sys_size
            ! Lagrangian bubbles: alphas sum to the LOCAL liquid fraction beta (not 1), so the sum-to-one closure would corrupt the
            ! EL state; each alpha prolongs plainly instead
            if (num_fluids > 1 .and. (.not. bubbles_lagrange) .and. i >= eqn_idx%adv%beg .and. i <= eqn_idx%adv%end) cycle
            if (chemistry .and. i >= eqn_idx%species%beg .and. i <= eqn_idx%species%end) cycle  ! sum/positivity closure below
            ! QBMM carries a bivariate 6-moment set per R0 bin whose CHyQMOM inversion requires realizability (variance c20 =
            ! m20/m00 - (m10/m00)^2 > 0); per-component minmod prolongation can break that joint constraint, so the whole bub block
            ! is injected piecewise-constant (each child inherits the coarse cell's realizable moment set exactly). Non-QBMM
            ! Euler-Euler bubbles instead floor their POSITIVE moments (radius nR, non-polytropic partial pressure npb / vapor mass
            ! nmv); the signed velocity moment nV (offset 1 in each bin's stride) prolongs freely.
            call s_prolong_one_var(amr_cg(i), amr_loc_of(amr_cur), i, &
                                   & pos=bubbles_euler .and. .not. qbmm .and. i >= eqn_idx%bub%beg .and. i <= eqn_idx%bub%end &
                                   & .and. mod(i - eqn_idx%bub%beg, bstride) /= 1, &
                                   & inject=qbmm .and. i >= eqn_idx%bub%beg .and. i <= eqn_idx%bub%end)
        end do
        if (num_fluids > 1 .and. (.not. bubbles_lagrange)) call s_prolong_alphas_closure(amr_cg, amr_loc_of(amr_cur))
        if (chemistry) call s_prolong_species_closure(amr_cg, amr_loc_of(amr_cur))

    end subroutine s_interpolate_coarse_to_fine

    !> Sum-preserving volume-fraction prolongation (num_fluids > 1): fluids adv%beg..adv%end-1 are interpolated with minmod slopes
    !! under a SHARED per-cell limiter switch (a sign change for ANY fluid in a dim zeroes that dim's slope for ALL fluids, so the
    !! closure fluid's effective slope is limited consistently) and clamped to [0,1]; the last fluid closes alpha_n = 1 -
    !! sum(others), so sum(alpha) = 1 on the fine level by construction. For two fluids the closure is also in [0,1]; for >2 fluids
    !! any residual closure undershoot is handled by mpp_lim (required by the checker). Same fine/coarse index mapping as
    !! s_prolong_one_var.
    impure subroutine s_prolong_alphas_closure(qc, loc)

        type(scalar_field), dimension(sys_size), intent(in) :: qc
        integer, intent(in)                                 :: loc
        integer                                             :: fi, fj, fk, ci, cj, ck, ox, oy, oz, i
        integer                                             :: rrat, mm, nn, pp, il1, il2, il3, advb, adve
        real(wp)                                            :: xix, xiy, xiz, u0, sx, sy, sz, av, asum
        logical                                             :: shx, shy, shz, d2, d3

        ! coarse source qc is the gathered block-local patch amr_cg (fine-level distribution): patch-frame offset.
        ! DEVICE kernel (device-side rebuild); the shared limiter switch (s_alpha_shared_switch) is inlined verbatim.

        ox = amr_cpat_off(1); oy = amr_cpat_off(2); oz = amr_cpat_off(3)
        rrat = amr_slots(amr_cur)%amr_ref_ratio
        mm = amr_slots(amr_cur)%m; nn = amr_slots(amr_cur)%n; pp = amr_slots(amr_cur)%p
        il1 = amr_isect_lo(1); il2 = amr_isect_lo(2); il3 = amr_isect_lo(3)
        d2 = n_glb > 0; d3 = p_glb > 0
        advb = eqn_idx%adv%beg; adve = eqn_idx%adv%end
        $:GPU_PARALLEL_LOOP(collapse=3, private='[ci, cj, ck, xix, xiy, xiz, u0, sx, sy, sz, av, asum, shx, shy, shz, i]')
        do fk = 0, pp
            do fj = 0, nn
                do fi = 0, mm
                    ck = il3 + fk/rrat - oz; if (.not. d3) ck = 0
                    xiz = 0._wp; if (d3) xiz = (real(mod(fk, rrat), wp) - real(rrat - 1, wp)*0.5_wp)/real(rrat, wp)
                    cj = il2 + fj/rrat - oy; if (.not. d2) cj = 0
                    xiy = 0._wp; if (d2) xiy = (real(mod(fj, rrat), wp) - real(rrat - 1, wp)*0.5_wp)/real(rrat, wp)
                    ci = il1 + fi/rrat - ox
                    xix = (real(mod(fi, rrat), wp) - real(rrat - 1, wp)*0.5_wp)/real(rrat, wp)
                    ! shared per-cell limiter switch (inlined s_alpha_shared_switch): per dim, slopes stay on only if NO
                    ! fluid's centered differences change sign there (symmetric in the fluids, incl. the closure fluid)
                    shx = .true.; shy = d2; shz = d3
                    do i = advb, adve
                        u0 = real(qc(i)%sf(ci, cj, ck), wp)
                        if ((real(qc(i)%sf(ci + 1, cj, ck), wp) - u0)*(u0 - real(qc(i)%sf(ci - 1, cj, ck), &
                            & wp)) <= 0._wp) shx = .false.
                        if (d2) then
                            if ((real(qc(i)%sf(ci, cj + 1, ck), wp) - u0)*(u0 - real(qc(i)%sf(ci, cj - 1, ck), &
                                & wp)) <= 0._wp) shy = .false.
                        end if
                        if (d3) then
                            if ((real(qc(i)%sf(ci, cj, ck + 1), wp) - u0)*(u0 - real(qc(i)%sf(ci, cj, ck - 1), &
                                & wp)) <= 0._wp) shz = .false.
                        end if
                    end do
                    asum = 0._wp
                    do i = advb, adve - 1
                        u0 = real(qc(i)%sf(ci, cj, ck), wp)
                        sx = 0._wp
                        if (shx) sx = minmod(real(qc(i)%sf(ci + 1, cj, ck), wp) - u0, u0 - real(qc(i)%sf(ci - 1, cj, ck), wp))
                        sy = 0._wp
                        if (d2 .and. shy) sy = minmod(real(qc(i)%sf(ci, cj + 1, ck), wp) - u0, u0 - real(qc(i)%sf(ci, cj - 1, &
                            & ck), wp))
                        sz = 0._wp
                        if (d3 .and. shz) sz = minmod(real(qc(i)%sf(ci, cj, ck + 1), wp) - u0, u0 - real(qc(i)%sf(ci, cj, &
                            & ck - 1), wp))
                        av = min(max(u0 + sx*xix + sy*xiy + sz*xiz, 0._wp), 1._wp)
                        amr_cons_st(fi, fj, fk, i, loc) = av
                        asum = asum + av
                    end do
                    amr_cons_st(fi, fj, fk, adve, loc) = 1._wp - asum
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_prolong_alphas_closure

    !> Species mass-fraction prolongation closure (chemistry): each partial density rho*Y_k is minmod-prolonged and clamped
    !! non-negative, then all species are rescaled so sum_k(rho*Y_k) equals the (already prolonged) continuity density at the fine
    !! cell. This keeps the fine species realizable (Y_k >= 0, and sum(Y_k) = 1 exactly under the cons->prim recovery rho = sum
    !! rho*Y_k) and consistent with the continuity variable the reaction source reads. Same index mapping as s_prolong_one_var; cont
    !! is prolonged in the main loop before this runs.
    impure subroutine s_prolong_species_closure(qc, loc)

        type(scalar_field), dimension(sys_size), intent(in) :: qc
        integer, intent(in)                                 :: loc
        integer                                             :: fi, fj, fk, ci, cj, ck, ox, oy, oz, i
        integer                                             :: rrat, mm, nn, pp, il1, il2, il3, spb, spe, cte
        real(wp)                                            :: xix, xiy, xiz, u0, sx, sy, sz, av, rsum, rscale
        logical                                             :: d2, d3

        ! coarse source qc is the gathered block-local patch amr_cg (fine-level distribution): patch-frame offset.
        ! DEVICE kernel (device-side rebuild); the rescale re-reads only this thread's own cell, so the loop nest is safe.

        ox = amr_cpat_off(1); oy = amr_cpat_off(2); oz = amr_cpat_off(3)
        rrat = amr_slots(amr_cur)%amr_ref_ratio
        mm = amr_slots(amr_cur)%m; nn = amr_slots(amr_cur)%n; pp = amr_slots(amr_cur)%p
        il1 = amr_isect_lo(1); il2 = amr_isect_lo(2); il3 = amr_isect_lo(3)
        d2 = n_glb > 0; d3 = p_glb > 0
        spb = eqn_idx%species%beg; spe = eqn_idx%species%end; cte = eqn_idx%cont%end
        $:GPU_PARALLEL_LOOP(collapse=3, private='[ci, cj, ck, xix, xiy, xiz, u0, sx, sy, sz, av, rsum, rscale, i]')
        do fk = 0, pp
            do fj = 0, nn
                do fi = 0, mm
                    ck = il3 + fk/rrat - oz; if (.not. d3) ck = 0
                    xiz = 0._wp; if (d3) xiz = (real(mod(fk, rrat), wp) - real(rrat - 1, wp)*0.5_wp)/real(rrat, wp)
                    cj = il2 + fj/rrat - oy; if (.not. d2) cj = 0
                    xiy = 0._wp; if (d2) xiy = (real(mod(fj, rrat), wp) - real(rrat - 1, wp)*0.5_wp)/real(rrat, wp)
                    ci = il1 + fi/rrat - ox
                    xix = (real(mod(fi, rrat), wp) - real(rrat - 1, wp)*0.5_wp)/real(rrat, wp)
                    rsum = 0._wp
                    do i = spb, spe
                        u0 = real(qc(i)%sf(ci, cj, ck), wp)
                        sx = minmod(real(qc(i)%sf(ci + 1, cj, ck), wp) - u0, u0 - real(qc(i)%sf(ci - 1, cj, ck), wp))
                        sy = 0._wp
                        if (d2) sy = minmod(real(qc(i)%sf(ci, cj + 1, ck), wp) - u0, u0 - real(qc(i)%sf(ci, cj - 1, ck), wp))
                        sz = 0._wp
                        if (d3) sz = minmod(real(qc(i)%sf(ci, cj, ck + 1), wp) - u0, u0 - real(qc(i)%sf(ci, cj, ck - 1), wp))
                        av = max(u0 + sx*xix + sy*xiy + sz*xiz, 0._wp)
                        amr_cons_st(fi, fj, fk, i, loc) = av
                        rsum = rsum + av
                    end do
                    rscale = real(amr_cons_st(fi, fj, fk, cte, loc), wp)/max(rsum, 1.e-30_wp)
                    do i = spb, spe
                        amr_cons_st(fi, fj, fk, i, loc) = real(amr_cons_st(fi, fj, fk, i, loc), wp)*rscale
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_prolong_species_closure

    !> Disblock prolongation. Guard: no-op unless amr.
    impure subroutine s_populate_amr_fine(q_cons_base)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_base
        integer                                                :: i, islot

        if (.not. amr) return
        ! Prolong EVERY block (max_grid_size tiling can make several) from its gathered coarse patch. The P2P gather pulls each
        ! patch's inter-rank coarse cells from neighbour interiors, so no coarse-ghost halo exchange is needed; host q_cons_base
        ! holds
        ! the ICs here (this runs before s_initialize_gpu_vars). ALL ranks call the gather (P2P); only owners prolong.
        do islot = f_l0_slot(1), amr_num_blocks
            call s_amr_select_slot(islot)
            call s_amr_gather_coarse_patch(q_cons_base, .false.)
            call s_amr_gather_send_flush()  ! keep this site's original blocking semantics
            ! non-polytropic QBMM: gather the coarse pb/mv patch too (ALL ranks - P2P; owners prolong from it below)
            if (qbmm .and. .not. polytropic) call s_amr_gather_coarse_patch_pbmv(pb_ts(1)%sf, mv_ts(1)%sf, .false.)
            if (amr_rank_owns_block) then
                ! the prolong is a device kernel now (writes the slot in place); no push - a host->device push here would
                ! clobber the device result with the stale host mirror
                call s_interpolate_coarse_to_fine()
                ! non-polytropic QBMM: seed the block's quadrature side-state from the coarse fields
                if (qbmm .and. .not. polytropic) call s_amr_prolong_pbmv()
            end if
        end do
        if (amr_max_level >= 2) call s_amr_build_static_multilevel(q_cons_base)
        call s_amr_select_slot(f_l0_slot(1))

    end subroutine s_populate_amr_fine

    !> Build the STATIC multi-level hierarchy (amr_regrid_int = 0): nest exactly one level-2 block inside level-1 block 1 by a fixed
    !! geometric inset (a regrid would place it by sensor-on-fine instead), prolong the parent state into it, and keep it persistent
    !! so the advance driver steps it every timestep. The restrict/reflux identity it relies on is protected by the static
    !! multi-level goldens (75AD6885 et al.) and the runtime conservation-defect probe.
    impure subroutine s_amr_build_static_multilevel(q_cons_base)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_base
        integer                                                :: L2, n1, i, par, inset(3)

        if (amr_max_level < 2) return
        n1 = amr_num_blocks
        if (n1 < 1) return
        ! the static hierarchy nests exactly one level-2 block; without pool room it would SILENTLY refine only to level 1 (an
        ! under-resolved but "successful" run). n1 (the level-1 tile count) is only known here, not at checker time, so abort at the
        ! point of failure. Replicated inputs -> every rank takes the same branch (collective-safe).
        if (n1 + 1 > l0_slot_off + amr_max_fine) call s_mpi_abort('amr static multi-level (amr_max_level > 1, ' &
            & // 'amr_regrid_int = 0): amr_max_blocks is too small to nest the level-2 block (need >= level-1 block count + 1); ' &
            & // 'increase amr_max_blocks')
        L2 = n1 + 1
        ! PARENT is the first FINE block, f_l0_slot(1) - NOT slot 1, which under coexist is the first level-0 TILE. Insetting a
        ! tile instead put the level-2 box in the wrong place and sized it off the tile: with one tile (the whole base grid) that
        ! tripped the amr_maxc_fit cap below, and with two it produced a plausible-looking box that silently corrupted the run.
        par = f_l0_slot(1)
        inset = 0
        inset(1) = max((amr_region_hi_all(1, par) - amr_region_lo_all(1, par) + 1)/4, amr_cpat_mar)
        if (n_glb > 0) inset(2) = max((amr_region_hi_all(2, par) - amr_region_lo_all(2, par) + 1)/4, amr_cpat_mar)
        if (p_glb > 0) inset(3) = max((amr_region_hi_all(3, par) - amr_region_lo_all(3, par) + 1)/4, amr_cpat_mar)
        amr_region_lo_all(:,L2) = amr_region_lo_all(:,par) + inset
        amr_region_hi_all(:,L2) = amr_region_hi_all(:,par) - inset
        ! Guard the fixed-inset box against configs this single-block static builder cannot represent - the dynamic regrid path has
        ! the analogous checks (proper-nesting skip + amr_maxc_fit/2 clamp), but the static path bypasses them. Replicated inputs ->
        ! every rank takes the same branch (collective-safe). (a) a level-1 block smaller than 2*inset inverts the box; (b) a
        ! level-2
        ! L0-extent > amr_maxc_fit/2 makes its parent-fine transverse extent (2*L0) overrun the creg register (allocated
        ! 0:amr_maxc_fit-1), a silent out-of-bounds device write in the L2->L1 reflux capture.
        if (amr_region_lo_all(1, L2) > amr_region_hi_all(1, L2) .or. (n_glb > 0 .and. amr_region_lo_all(2, &
            & L2) > amr_region_hi_all(2, L2)) .or. (p_glb > 0 .and. amr_region_lo_all(3, L2) > amr_region_hi_all(3, &
            & L2))) call s_mpi_abort('amr static multi-level: level-1 block 1 is too small to nest a level-2 block (the fixed ' &
            & // 'inset inverts the box); enlarge the base amr block or reduce amr_cpat_mar')
        if (amr_ref_ratio*(amr_region_hi_all(1, L2) - amr_region_lo_all(1, &
            & L2) + 1) > amr_maxc_fit(1) .or. (n_glb > 0 .and. amr_ref_ratio*(amr_region_hi_all(2, L2) - amr_region_lo_all(2, &
            & L2) + 1) > amr_maxc_fit(2)) .or. (p_glb > 0 .and. amr_ref_ratio*(amr_region_hi_all(3, L2) - amr_region_lo_all(3, &
            & L2) + 1) > amr_maxc_fit(3))) &
            & call s_mpi_abort('amr static multi-level: the nested level-2 block exceeds the per-rank scratch cap ' &
            & // '(2*L0-extent > amr_maxc_fit); static multi-level does not tile the level-2 block - use a smaller base amr ' &
            & // 'block or the dynamic regrid path (amr_regrid_int > 0)')
        amr_block_level(L2) = 2
        amr_block_owner(L2) = amr_block_owner(par)
        amr_num_blocks = L2; amr_num_levels = 2
        call s_amr_reconcile_slots()
        amr_cur = L2
        call s_set_amr_fine_geometry(amr_region_lo_all(:,L2), amr_region_hi_all(:,L2))
        call s_amr_reduce_xchg_flag()
        call s_amr_gather_coarse_patch(q_cons_base, .false.)  ! q_coarse ignored for level>=2 (reads the parent block); pass the
        call s_amr_gather_send_flush()  ! keep this site's original blocking semantics
        ! always-allocated base field, not amr_slots(1) (the parent slot is unallocated on a non-owner rank at np>1)
        if (amr_rank_owns_block) then
            ! the prolong is a device kernel now: the persistent L2 block's device q_cons is valued in place (the historical
            ! host-loop NaN hazard this site used to push against is gone; a push here would clobber the device result)
            call s_interpolate_coarse_to_fine()
        end if
        ! persistent L2 block: KEEP the level-2 block in the active set (amr_num_blocks = L2, amr_num_levels = 2) so the advance
        ! driver steps it across timesteps; no free/revert.
        ! restore amr_cg + the patch frame (amr_cpat_off) to the first FINE block: the L2 gather above overwrote them with the
        ! parent-fine frame, and the normal single-block conservation check that follows reads that block's frame. f_l0_slot(1),
        ! not slot 1 - under coexist slot 1 is a level-0 TILE, and selecting it here left the grid globals describing tile
        ! geometry for the rest of init, so s_initialize_weno_module (m_start_up, called after this) sized its device-mapped
        ! coefficient tables off the wrong bounds and faulted in __tgt_target_data_begin_mapper.
        call s_amr_select_slot(f_l0_slot(1))
        call s_amr_gather_coarse_patch(q_cons_base, .false.)
        call s_amr_gather_send_flush()  ! keep this site's original blocking semantics

    end subroutine s_amr_build_static_multilevel

    !> Volume-weighted restriction: each covered coarse cell = volume-weighted average of its amr_ref_ratio^d fine children (equal
    !! weight on Cartesian grids where children share a volume; radius-weighted by fine y_cc on cyl_coord, where cell volume ~
    !! radius - amr_rvw, single-sourced so device and host paths agree bit-for-bit). Writes the caller's coarse target - in
    !! production the level-0 state q_cons_ts(1)%vf (the deliberate fold-back of fine data each step, plus coarse pb/mv for
    !! non-polytropic QBMM); init-time diagnostics pass a scratch buffer instead. Device kernel.
    impure subroutine s_restrict_fine_to_coarse(coarse_tgt)

        type(scalar_field), dimension(sys_size), intent(inout) :: coarse_tgt
        integer :: nchild, rr, dj_hi, dk_hi, o1, o2, o3, owner, r, idx, boxsz, maxsz, nsrc, ierr
        integer :: rlo(3), rhi(3), ilo(3), ihi(3), bl(3), bh(3)
        real(wp), allocatable :: sbuf(:,:), rbuf(:)
        integer, allocatable :: reqs(:), drank(:)

        if (rank_time_wrt .and. amr_rank_owns_block) call s_rank_time_tic()

        ! multi-level: a level>=2 block folds back into its PARENT block's fine array (the coarse side of level l is level l-1), not
        ! the L0 coarse_tgt. Same restriction kernel, targeted at the parent in the parent-fine frame. When child and parent sit on
        ! different ranks the fold is a P2P pair, so BOTH participants must enter or the receiver never posts.
        if (amr_block_level(amr_cur) >= 2) then
            if (amr_rank_owns_block .or. amr_block_owner(f_amr_parent_block(amr_cur)) == proc_rank) then
                call s_amr_restrict_to_parent()
            end if
            if (rank_time_wrt .and. amr_rank_owns_block) call s_rank_time_toc()
            return
        end if

        ! whole-block-per-rank fold-back: the block owner restricts its fine block to coarse averages over the covered cells
        ! [region_lo:region_hi] and SCATTERS them POINT-TO-POINT to the coarse-cell owners - the owner overwrites the covered cells
        ! it holds locally and SENDS each other coarse-owner exactly its covered slice (all sys_size in one message). Covered cells
        ! are in-domain (no ghosts), so each is owned by exactly one interior owner. At np=1 the owner owns every covered cell,
        ! sends
        ! nothing, and overwrites locally with the same child-sum -> bit-identical.
        rr = amr_slots(amr_cur)%amr_ref_ratio
        nchild = rr; if (n_glb > 0) nchild = nchild*rr; if (p_glb > 0) nchild = nchild*rr
        dj_hi = merge(rr - 1, 0, n_glb > 0); dk_hi = merge(rr - 1, 0, p_glb > 0)
        rlo = 0; rhi = 0
        rlo(1) = amr_region_lo_all(1, amr_cur); rhi(1) = amr_region_hi_all(1, amr_cur)
        if (n_glb > 0) then; rlo(2) = amr_region_lo_all(2, amr_cur); rhi(2) = amr_region_hi_all(2, amr_cur); end if
        if (p_glb > 0) then; rlo(3) = amr_region_lo_all(3, amr_cur); rhi(3) = amr_region_hi_all(3, amr_cur); end if
        owner = amr_block_owner(amr_cur)
        o1 = start_idx(1); o2 = 0; o3 = 0
        if (n_glb > 0) o2 = start_idx(2)
        if (p_glb > 0) o3 = start_idx(3)
        maxsz = sys_size*(rhi(1) - rlo(1) + 1)*(rhi(2) - rlo(2) + 1)*(rhi(3) - rlo(3) + 1)

        ! cyl_coord: fine radial volume weights = this block's fine cell-center radii, pushed to device for the restriction kernels.
        ! Only the owner restricts (device overwrite + device scatter pack), so only the owner needs them; single-sourced from y_cc
        ! so the owner-local and scattered child-averages are bit-identical.
        if (cyl_coord .and. proc_rank == owner) then
            amr_rvw(0:amr_slots(amr_cur)%n) = amr_slots(amr_cur)%y_cc(0:amr_slots(amr_cur)%n)
            $:GPU_UPDATE(device='[amr_rvw]')
        end if

        ! block set changed: rebuild the cached overlap-rank lists (same lazy trigger as s_amr_fine_fine_halo; local, replicated)
        if (amr_seam_pairs_dirty .or. amr_seam_pairs_nblk /= amr_num_blocks) call s_amr_build_seam_pairs()

        if (proc_rank == owner) then
            ! overwrite the covered cells this rank owns, then send each other coarse-owner its covered slice
            call s_amr_rank_interior(proc_rank, ilo, ihi)
            call s_amr_box_isect(rlo, rhi, ilo, ihi, bl, bh)
            if (num_procs == 1) then
                ! np=1 device-native fold-back: restrict the fine block (device) into the coarse (device) over the COVERED cells
                ! only - no host round-trip. The old path pulled fine to host, restricted on host, then pushed the WHOLE coarse
                ! array back to the device (GPU_UPDATE device coarse_tgt), clobbering the device-advanced NON-covered coarse cells
                ! with the stale host copy - a GPU-only divergence (invisible on CPU where host==device) that IGR/MHD/acoustic
                ! amplify. The owner holds every covered cell at np=1.
                if (bl(1) <= bh(1) .and. bl(2) <= bh(2) .and. bl(3) <= bh(3)) call s_amr_restrict_overwrite_device_sf(coarse_tgt, &
                    & amr_loc_of(amr_cur), bl, bh, o1, o2, o3, rlo, rr, dj_hi, dk_hi, nchild)
                if (qbmm .and. .not. polytropic .and. amr_rank_owns_block) call s_restrict_pbmv(pb_ts(1)%sf, mv_ts(1)%sf, &
                    & amr_slots(amr_cur)%pb_f%sf, amr_slots(amr_cur)%mv_f%sf)
                if (rank_time_wrt .and. amr_rank_owns_block) call s_rank_time_toc()
                return
            end if
            ! owner-local covered cells: restrict fine(device) -> coarse(device) touching ONLY those cells (no whole-coarse device
            ! push, which clobbered the device-advanced non-covered coarse cells - the same GPU-only bug fixed at np=1)
            if (bl(1) <= bh(1) .and. bl(2) <= bh(2) .and. bl(3) <= bh(3)) call s_amr_restrict_overwrite_device_sf(coarse_tgt, &
                & amr_loc_of(amr_cur), bl, bh, o1, o2, o3, rlo, rr, dj_hi, dk_hi, nchild)
            ! cached destination list (every listed rank's interior overlaps the region by construction)
            nsrc = 0
            do idx = 1, amr_ovl_scatter_n(amr_cur)
                if (amr_ovl_scatter(idx, amr_cur) /= owner) nsrc = nsrc + 1
            end do
            if (nsrc > 0) then
                allocate (sbuf(maxsz, nsrc), reqs(nsrc), drank(nsrc))
                nsrc = 0
                do idx = 1, amr_ovl_scatter_n(amr_cur)
                    r = amr_ovl_scatter(idx, amr_cur)
                    if (r == owner) cycle
                    call s_amr_rank_interior(r, ilo, ihi)
                    call s_amr_box_isect(rlo, rhi, ilo, ihi, bl, bh)
                    nsrc = nsrc + 1; drank(nsrc) = r
                    boxsz = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                    ! pack this destination's covered slice on the DEVICE: restrict averages straight into the wire buffer (same
                    ! child-sum order and wp values as the device overwrite above) - no full-field host pull
                    call s_amr_restrict_pack_device(amr_loc_of(amr_cur), bl, bh, rlo, rr, dj_hi, dk_hi, nchild, sbuf(1:boxsz,nsrc))
#ifdef MFC_MPI
                    call s_xa_rec(XA_F7A_SND, 1, boxsz, amr_cur)
                    call MPI_ISEND(sbuf(1, nsrc), boxsz, mpi_p, r, amr_cur, MPI_COMM_WORLD, reqs(nsrc), ierr)
#endif
                end do
#ifdef MFC_MPI
                call MPI_WAITALL(nsrc, reqs, MPI_STATUSES_IGNORE, ierr)
#endif
                deallocate (sbuf, reqs, drank)
            end if
        else
            ! coarse-owner: if I hold covered cells, receive my slice from the owner and overwrite my local coarse
            call s_amr_rank_interior(proc_rank, ilo, ihi)
            call s_amr_box_isect(rlo, rhi, ilo, ihi, bl, bh)
            if (bl(1) <= bh(1) .and. bl(2) <= bh(2) .and. bl(3) <= bh(3)) then
                boxsz = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                allocate (rbuf(boxsz))
#ifdef MFC_MPI
                call s_xa_rec(XA_F7A_RCV, 2, boxsz, amr_cur)
                call MPI_RECV(rbuf, boxsz, mpi_p, owner, amr_cur, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
#endif
                ! DEVICE unpack of the covered box, writing only those cells (a whole-array push would clobber the device-advanced
                ! non-covered coarse cells with this rank's stale host copy - the GPU-only bug fixed at np=1). This must NOT be a
                ! host unpack followed by a strided GPU_UPDATE(device=) of the box: a non-contiguous 3-D array section in an OpenMP
                ! target update is copied by AMD flang as size(section) CONTIGUOUS elements, so only the first row lands on the
                ! cells it names and the remainder overwrites neighbouring cells with stale host data - silently, and only at
                ! np >= 2 with a block whose owner holds none of its covered cells. The wire layout (ci fastest, then cj, ck, i)
                ! is exactly s_l0_pack_unpack_block's, so it unpacks s_amr_restrict_pack_device's buffer as-is.
                call s_l0_pack_unpack_block_sf(coarse_tgt, bl(1) - o1, bl(2) - o2, bl(3) - o3, bh(1) - bl(1), bh(2) - bl(2), &
                                               & bh(3) - bl(3), rbuf, .false.)
                deallocate (rbuf)
            end if
        end if

        ! non-polytropic QBMM: distributed pb/mv fold-back - the owner restricts the covered cells it holds and scatters each other
        ! coarse-owner its slice (mirror of the q_cons scatter above). Called on ALL ranks so the P2P send/recv pair up (np=1
        ! handled
        ! by the direct s_restrict_pbmv in the num_procs==1 branch above, which returns before reaching here)
        if (qbmm .and. .not. polytropic) call s_amr_scatter_pbmv(amr_slots(amr_cur)%pb_f%sf, amr_slots(amr_cur)%mv_f%sf)
        if (rank_time_wrt .and. amr_rank_owns_block) call s_rank_time_toc()

    end subroutine s_restrict_fine_to_coarse

    !> Multi-level restriction: fold the current level>=2 block's fine averages back into its PARENT block's fine array over the
    !! covered cells. Same child-sum kernel as the L0 fold-back, targeted at the parent in the parent-fine frame (amr_isect already
    !! parent-fine; offset 0 = the parent's local fine indexing). np=1 local; the np>=2 P2P scatter is future work.
    impure subroutine s_amr_restrict_to_parent()

        integer               :: pblk, rr, nchild, dj_hi, dk_hi, cowner, powner, boxsz, ierr
        integer               :: plo(3), phi(3)
        real(wp), allocatable :: xbuf(:)

        ! NOTE: reads amr_rvw's device copy without a GPU_UPDATE - safe only while cyl_coord + amr_max_level > 1 is checker-gated
        ! (this path never runs under cyl_coord). If that gate lifts, refresh amr_rvw here first (see its declaration).

        pblk = f_amr_parent_block(amr_cur)
        cowner = amr_block_owner(amr_cur); powner = amr_block_owner(pblk)
        if (proc_rank /= cowner .and. proc_rank /= powner) return  ! not a participant

        ! Same replicated-metadata box as the gather, so the folding child and the receiving parent agree without a handshake.
        call s_amr_parent_foot(amr_cur, pblk, plo, phi)
        if (plo(1) > phi(1) .or. plo(2) > phi(2) .or. plo(3) > phi(3)) return  ! empty footprint

        rr = amr_ref_ratio
        nchild = rr; if (n_glb > 0) nchild = nchild*rr; if (p_glb > 0) nchild = nchild*rr
        dj_hi = merge(rr - 1, 0, n_glb > 0); dk_hi = merge(rr - 1, 0, p_glb > 0)

        if (powner == cowner) then
            ! co-located (np=1, or a co-located tower): fold straight into the parent, bit-for-bit as before.
            call s_amr_restrict_overwrite_device_st(amr_loc_of(pblk), amr_loc_of(amr_cur), plo, phi, 0, 0, 0, plo, rr, dj_hi, &
                                                    & dk_hi, nchild)
            return
        end if

#ifdef MFC_MPI
        ! Split ownership: the CHILD restricts locally and ships COARSE cells - rr**num_dims fewer values than shipping its fine
        ! block - which restriction being an overwrite (not an accumulate) makes correct. Reuses the L0<->L1 scatter's pack and
        ! unpack verbatim; their wire layout (ci fastest, then cj, ck, i) is already documented as compatible.
        boxsz = sys_size*(phi(1) - plo(1) + 1)*(phi(2) - plo(2) + 1)*(phi(3) - plo(3) + 1)
        allocate (xbuf(boxsz))
        if (proc_rank == cowner) then
            call s_amr_restrict_pack_device(amr_loc_of(amr_cur), plo, phi, plo, rr, dj_hi, dk_hi, nchild, xbuf)
            call s_xa_rec(XA_F7B_SND, 1, boxsz, amr_cur)
            call MPI_SEND(xbuf, boxsz, mpi_p, powner, amr_cur, MPI_COMM_WORLD, ierr)
        else
            call s_xa_rec(XA_F7B_RCV, 2, boxsz, amr_cur)
            call MPI_RECV(xbuf, boxsz, mpi_p, cowner, amr_cur, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
            ! DEVICE unpack of just the covered box - never a host unpack plus a strided GPU_UPDATE (see the L0 scatter's note: AMD
            ! flang copies a non-contiguous 3-D section as contiguous elements and silently corrupts neighbouring cells).
            call s_l0_pack_unpack_block_st(amr_loc_of(pblk), plo(1), plo(2), plo(3), phi(1) - plo(1), phi(2) - plo(2), &
                                           & phi(3) - plo(3), xbuf, .false.)
        end if
        deallocate (xbuf)
#endif

    end subroutine s_amr_restrict_to_parent

    !> I5b: the lock-step restrict fold as plan-based waves (amr_plan_based_exchange.md). Replaces the per-box reverse fold loop on
    !! the lock-step np>1 path: per level (finest first) the split-ownership child->parent folds become ONE aggregated message per
    !! (child-owner, parent-owner) peer (s_amr_restrict_parent_wave), followed by the per-box reflux-to-parent applies (their freg
    !! was already exchanged by s_amr_freg_wave); then the level-1 -> L0 covered-cell scatter becomes one aggregated message per
    !! (owner, coarse-owner) peer (s_amr_restrict_l1_wave). Both sides derive identical transfer lists from replicated metadata (the
    !! region-box x rank-interior intersections the per-box path already used), so the wire layout needs no handshake and F7 family
    !! words stay EXACT; what is removed is the per-box ISEND+WAITALL / blocking-RECV chain, whose length scales with the GLOBAL
    !! block count. Level order (finest first) preserves child-before-parent folding; within a level the covered targets are
    !! disjoint, and the restrict/reflux interleave is order-free because sibling-shared faces carry weight 0
    !! (s_amr_sibling_face_weights). Subcycle and np=1 keep the per-box loop (the per-box routines remain).
    impure subroutine s_amr_restrict_wave(coarse_tgt, dt_reflux)

        type(scalar_field), dimension(sys_size), intent(inout) :: coarse_tgt
        real(wp), intent(in)                                   :: dt_reflux
        integer                                                :: lev, k

        do lev = amr_max_level, 2, -1
            if (relax) then
                do k = amr_num_blocks, 1, -1
                    if (amr_block_level(k) /= lev) cycle
                    call s_amr_select_slot(k)
                    call s_amr_relax_fine()
                end do
            end if
            call s_phase_tic(PH_RESTR); call s_phase_tic(PH_RSREST)
            call s_amr_restrict_parent_wave(lev)
            call s_phase_toc(PH_RSREST); call s_phase_toc(PH_RESTR)
            do k = amr_num_blocks, 1, -1
                if (amr_block_level(k) /= lev) cycle
                call s_amr_select_slot(k)
                call s_phase_tic(PH_RESTR); call s_phase_tic(PH_RSRFP)
                call s_amr_reflux_to_parent(dt_reflux, .false.)
                call s_phase_toc(PH_RSRFP); call s_phase_toc(PH_RESTR)
            end do
        end do
        if (relax) then
            do k = amr_num_blocks, 1, -1
                if (amr_block_level(k) /= 1) cycle
                call s_amr_select_slot(k)
                call s_amr_relax_fine()
            end do
        end if
        call s_phase_tic(PH_RESTR); call s_phase_tic(PH_RSREST)
        call s_amr_restrict_l1_wave(coarse_tgt)
        call s_phase_toc(PH_RSREST); call s_phase_toc(PH_RESTR)
        ! non-polytropic QBMM pb/mv fold-back keeps its per-box pairing (every rank walks the same reverse order, so the
        ! blocking pairs match exactly as they did inside the per-box fold)
        if (qbmm .and. .not. polytropic) then
            do k = amr_num_blocks, 1, -1
                if (amr_block_level(k) /= 1) cycle
                call s_amr_select_slot(k)
                call s_amr_scatter_pbmv(amr_slots(amr_cur)%pb_f%sf, amr_slots(amr_cur)%mv_f%sf)
            end do
        end if
        call s_amr_select_slot(1)

    end subroutine s_amr_restrict_wave

    !> One per-level wave of the split-ownership level>=2 child->parent restrict folds (the F7B pairs). Co-located folds run inline
    !! on the child-owner (bit-for-bit the per-box kernel); every cross-rank fold ships its parent-frame covered box in one
    !! aggregated message per (child-owner, parent-owner) peer. Both sides walk the same replicated block list ascending with
    !! per-peer running offsets, so the wire layout agrees with no metadata exchange.
    impure subroutine s_amr_restrict_parent_wave(lev)

        integer, intent(in) :: lev

#ifdef MFC_MPI
        integer :: k, pblk, cowner, powner, rr, nchild, dj_hi, dk_hi, ierr, ip, idx, r, cnt, boff, qbase, nreq, tq
        integer :: plo(3), phi(3)

        rr = amr_ref_ratio
        nchild = rr; if (n_glb > 0) nchild = nchild*rr; if (p_glb > 0) nchild = nchild*rr
        dj_hi = merge(rr - 1, 0, n_glb > 0); dk_hi = merge(rr - 1, 0, p_glb > 0)
        tq = amr_tag_base(7) + int(mod(amr_mesh_epoch, 50_8))
        if (.not. allocated(amr_fw_map)) then
            allocate (amr_fw_map(0:num_procs - 1), amr_fw_nx(0:num_procs - 1), amr_fw_pq(0:num_procs - 1), &
                      & amr_fw_pp(0:num_procs - 1))
            amr_fw_map = 0; amr_fw_nx = 0; amr_fw_pq = 0; amr_fw_pp = 0
        end if
        ! send plan + co-located folds (child-owner side)
        amr_fw_snx = 0; amr_fw_snp = 0
        do k = 1, amr_num_blocks
            if (amr_block_level(k) /= lev) cycle
            cowner = amr_block_owner(k)
            if (proc_rank /= cowner) cycle
            pblk = f_amr_parent_block(k)
            powner = amr_block_owner(pblk)
            call s_amr_parent_foot(k, pblk, plo, phi)
            if (plo(1) > phi(1) .or. plo(2) > phi(2) .or. plo(3) > phi(3)) cycle
            if (cowner == powner) then
                call s_amr_restrict_overwrite_device_st(amr_loc_of(pblk), amr_loc_of(k), plo, phi, 0, 0, 0, plo, rr, dj_hi, &
                                                        & dk_hi, nchild)
                cycle
            end if
            if (amr_fw_map(powner) == 0) then
                amr_fw_snp = amr_fw_snp + 1
                call s_amr_fw_szi(amr_fw_sprank, amr_fw_snp); call s_amr_fw_szi(amr_fw_sqsz, amr_fw_snp)
                call s_amr_fw_szi(amr_fw_snxp, amr_fw_snp); call s_amr_fw_szi(amr_fw_sqbase, amr_fw_snp)
                amr_fw_map(powner) = amr_fw_snp
                amr_fw_sprank(amr_fw_snp) = powner
            end if
            cnt = sys_size*(phi(1) - plo(1) + 1)*(phi(2) - plo(2) + 1)*(phi(3) - plo(3) + 1)
            amr_fw_snx = amr_fw_snx + 1
            call s_amr_fw_szi(amr_fw_sblk, amr_fw_snx); call s_amr_fw_szi3(amr_fw_sbl, amr_fw_snx)
            call s_amr_fw_szi3(amr_fw_sbh, amr_fw_snx); call s_amr_fw_szi(amr_fw_spi, amr_fw_snx)
            call s_amr_fw_szi(amr_fw_sqo, amr_fw_snx); call s_amr_fw_szi(amr_fw_spo, amr_fw_snx)
            amr_fw_sblk(amr_fw_snx) = k
            amr_fw_sbl(:,amr_fw_snx) = plo; amr_fw_sbh(:,amr_fw_snx) = phi
            amr_fw_spo(amr_fw_snx) = cnt
            amr_fw_spi(amr_fw_snx) = amr_fw_map(powner)
            amr_fw_sqo(amr_fw_snx) = amr_fw_pq(powner) + amr_fw_nx(powner)*XA_NH
            amr_fw_pq(powner) = amr_fw_pq(powner) + cnt
            amr_fw_nx(powner) = amr_fw_nx(powner) + 1
        end do
        qbase = 0
        do ip = 1, amr_fw_snp
            r = amr_fw_sprank(ip)
            amr_fw_snxp(ip) = amr_fw_nx(r)
            amr_fw_sqsz(ip) = amr_fw_pq(r) + amr_fw_nx(r)*XA_NH
            amr_fw_sqbase(ip) = qbase; qbase = qbase + amr_fw_sqsz(ip)
            amr_fw_map(r) = 0; amr_fw_nx(r) = 0; amr_fw_pq(r) = 0
        end do
        call s_amr_fw_szr(amr_fw_sq, qbase)
        ! recv plan (parent-owner side): the same replicated walk, so per-peer transfer order matches the sender's
        amr_fw_rnx = 0; amr_fw_rnp = 0
        do k = 1, amr_num_blocks
            if (amr_block_level(k) /= lev) cycle
            pblk = f_amr_parent_block(k)
            cowner = amr_block_owner(k); powner = amr_block_owner(pblk)
            if (cowner == powner .or. proc_rank /= powner) cycle
            call s_amr_parent_foot(k, pblk, plo, phi)
            if (plo(1) > phi(1) .or. plo(2) > phi(2) .or. plo(3) > phi(3)) cycle
            if (amr_fw_map(cowner) == 0) then
                amr_fw_rnp = amr_fw_rnp + 1
                call s_amr_fw_szi(amr_fw_rprank, amr_fw_rnp); call s_amr_fw_szi(amr_fw_rqsz, amr_fw_rnp)
                call s_amr_fw_szi(amr_fw_rnxp, amr_fw_rnp); call s_amr_fw_szi(amr_fw_rqbase, amr_fw_rnp)
                amr_fw_map(cowner) = amr_fw_rnp
                amr_fw_rprank(amr_fw_rnp) = cowner
            end if
            cnt = sys_size*(phi(1) - plo(1) + 1)*(phi(2) - plo(2) + 1)*(phi(3) - plo(3) + 1)
            amr_fw_rnx = amr_fw_rnx + 1
            call s_amr_fw_szi(amr_fw_rblk, amr_fw_rnx); call s_amr_fw_szi3(amr_fw_rbl, amr_fw_rnx)
            call s_amr_fw_szi3(amr_fw_rbh, amr_fw_rnx); call s_amr_fw_szi(amr_fw_rpi, amr_fw_rnx)
            call s_amr_fw_szi(amr_fw_rqo, amr_fw_rnx); call s_amr_fw_szi(amr_fw_rpo, amr_fw_rnx)
            amr_fw_rblk(amr_fw_rnx) = k
            amr_fw_rbl(:,amr_fw_rnx) = plo; amr_fw_rbh(:,amr_fw_rnx) = phi
            amr_fw_rpo(amr_fw_rnx) = cnt
            amr_fw_rpi(amr_fw_rnx) = amr_fw_map(cowner)
            amr_fw_rqo(amr_fw_rnx) = amr_fw_pq(cowner) + amr_fw_nx(cowner)*XA_NH
            amr_fw_pq(cowner) = amr_fw_pq(cowner) + cnt
            amr_fw_nx(cowner) = amr_fw_nx(cowner) + 1
        end do
        qbase = 0
        do ip = 1, amr_fw_rnp
            r = amr_fw_rprank(ip)
            amr_fw_rnxp(ip) = amr_fw_nx(r)
            amr_fw_rqsz(ip) = amr_fw_pq(r) + amr_fw_nx(r)*XA_NH
            amr_fw_rqbase(ip) = qbase; qbase = qbase + amr_fw_rqsz(ip)
            amr_fw_map(r) = 0; amr_fw_nx(r) = 0; amr_fw_pq(r) = 0
        end do
        call s_amr_fw_szr(amr_fw_rq, qbase)
        nreq = amr_fw_snp + amr_fw_rnp
        if (nreq == 0) return
        call s_amr_fw_szi(amr_fw_req, nreq); call s_amr_fw_szi(amr_fw_reqw, nreq)
        nreq = 0
        do ip = 1, amr_fw_rnp
            call s_xa_rec(XA_F7BW_RCV, 2, amr_fw_rqsz(ip) - amr_fw_rnxp(ip)*XA_NH, tq)
            nreq = nreq + 1; amr_fw_reqw(nreq) = amr_fw_rqsz(ip)
            call MPI_IRECV(amr_fw_rq(amr_fw_rqbase(ip) + 1), amr_fw_rqsz(ip), mpi_p, amr_fw_rprank(ip), tq, MPI_COMM_WORLD, &
                           & amr_fw_req(nreq), ierr)
        end do
        do idx = 1, amr_fw_snx
            cnt = amr_fw_spo(idx)
            boff = amr_fw_sqbase(amr_fw_spi(idx)) + amr_fw_sqo(idx)
            call s_amr_restrict_pack_device(amr_loc_of(amr_fw_sblk(idx)), amr_fw_sbl(:,idx), amr_fw_sbh(:,idx), amr_fw_sbl(:, &
                                            & idx), rr, dj_hi, dk_hi, nchild, amr_fw_sq(boff + XA_NH + 1:boff + XA_NH + cnt))
            if (XA_NH > 0) call s_xa_hdr_pack(amr_fw_sq(boff + 1:boff + XA_NH), XA_F7BW_SND, amr_fw_sblk(idx), amr_fw_sbl(:,idx), &
                & amr_fw_sbh(:,idx))
        end do
        do ip = 1, amr_fw_snp
            call s_xa_rec(XA_F7BW_SND, 1, amr_fw_sqsz(ip) - amr_fw_snxp(ip)*XA_NH, tq)
            nreq = nreq + 1; amr_fw_reqw(nreq) = -1
            call MPI_ISEND(amr_fw_sq(amr_fw_sqbase(ip) + 1), amr_fw_sqsz(ip), mpi_p, amr_fw_sprank(ip), tq, MPI_COMM_WORLD, &
                           & amr_fw_req(nreq), ierr)
        end do
#ifdef MFC_DEBUG
        block
            integer :: st(MPI_STATUS_SIZE, nreq), gotw, q
            call MPI_WAITALL(nreq, amr_fw_req, st, ierr)
            do q = 1, nreq
                if (amr_fw_reqw(q) < 0) cycle
                call MPI_GET_COUNT(st(:,q), mpi_p, gotw, ierr)
                @:ASSERT(gotw == amr_fw_reqw(q), "restrict parent wave: received message length differs from the plan")
            end do
        end block
#else
        call MPI_WAITALL(nreq, amr_fw_req, MPI_STATUSES_IGNORE, ierr)
#endif
        do idx = 1, amr_fw_rnx
            cnt = amr_fw_rpo(idx)
            boff = amr_fw_rqbase(amr_fw_rpi(idx)) + amr_fw_rqo(idx)
            if (XA_NH > 0) call s_xa_hdr_check(amr_fw_rq(boff + 1:boff + XA_NH), XA_F7BW_SND, amr_fw_rblk(idx), amr_fw_rbl(:, &
                & idx), amr_fw_rbh(:,idx))
            ! DEVICE unpack of the covered box only (the strided-update flang trap - see s_restrict_fine_to_coarse)
            call s_l0_pack_unpack_block_st(amr_loc_of(f_amr_parent_block(amr_fw_rblk(idx))), amr_fw_rbl(1, idx), amr_fw_rbl(2, &
                                           & idx), amr_fw_rbl(3, idx), amr_fw_rbh(1, idx) - amr_fw_rbl(1, idx), amr_fw_rbh(2, &
                                           & idx) - amr_fw_rbl(2, idx), amr_fw_rbh(3, idx) - amr_fw_rbl(3, idx), &
                                           & amr_fw_rq(boff + XA_NH + 1:boff + XA_NH + cnt), .false.)
        end do
#endif

    end subroutine s_amr_restrict_parent_wave

    !> The level-1 -> L0 covered-cell scatter (F7A) as one wave: every owned level-1 block's covered slabs for every listed
    !! coarse-owner ship in one aggregated message per peer; the owner-local covered overwrite and the (cyl_coord) amr_rvw push stay
    !! grouped per block during the pack walk. The receiver plan is my-interior x region(k) over the level-1 blocks I do not own -
    !! by construction (s_amr_ranks_overlapping) exactly the sender's list membership.
    impure subroutine s_amr_restrict_l1_wave(coarse_tgt)

        type(scalar_field), dimension(sys_size), intent(inout) :: coarse_tgt

#ifdef MFC_MPI
        integer :: k, owner, rr, nchild, dj_hi, dk_hi, ierr, ip, idx, r, cnt, boff, qbase, nreq, tq, o1, o2, o3, cur
        integer :: rlo(3), rhi(3), ilo(3), ihi(3), milo(3), mihi(3), bl(3), bh(3)

        tq = amr_tag_base(7) + 50 + int(mod(amr_mesh_epoch, 50_8))
        o1 = start_idx(1); o2 = 0; o3 = 0
        if (n_glb > 0) o2 = start_idx(2)
        if (p_glb > 0) o3 = start_idx(3)
        ! block set changed: rebuild the cached overlap-rank lists (same lazy trigger as the per-box path)
        if (amr_seam_pairs_dirty .or. amr_seam_pairs_nblk /= amr_num_blocks) call s_amr_build_seam_pairs()
        if (.not. allocated(amr_fw_map)) then
            allocate (amr_fw_map(0:num_procs - 1), amr_fw_nx(0:num_procs - 1), amr_fw_pq(0:num_procs - 1), &
                      & amr_fw_pp(0:num_procs - 1))
            amr_fw_map = 0; amr_fw_nx = 0; amr_fw_pq = 0; amr_fw_pp = 0
        end if
        call s_amr_rank_interior(proc_rank, milo, mihi)
        ! send plan (block-owner side): the same (interior x region) covered slabs the per-box path sent, k-grouped
        amr_fw_snx = 0; amr_fw_snp = 0
        do k = 1, amr_num_blocks
            if (amr_block_level(k) /= 1) cycle
            if (amr_block_owner(k) /= proc_rank) cycle
            rlo = 0; rhi = 0
            rlo(1) = amr_region_lo_all(1, k); rhi(1) = amr_region_hi_all(1, k)
            if (n_glb > 0) then; rlo(2) = amr_region_lo_all(2, k); rhi(2) = amr_region_hi_all(2, k); end if
            if (p_glb > 0) then; rlo(3) = amr_region_lo_all(3, k); rhi(3) = amr_region_hi_all(3, k); end if
            do idx = 1, amr_ovl_scatter_n(k)
                r = amr_ovl_scatter(idx, k)
                if (r == proc_rank) cycle
                call s_amr_rank_interior(r, ilo, ihi)
                call s_amr_box_isect(rlo, rhi, ilo, ihi, bl, bh)
                if (bl(1) > bh(1) .or. bl(2) > bh(2) .or. bl(3) > bh(3)) cycle
                if (amr_fw_map(r) == 0) then
                    amr_fw_snp = amr_fw_snp + 1
                    call s_amr_fw_szi(amr_fw_sprank, amr_fw_snp); call s_amr_fw_szi(amr_fw_sqsz, amr_fw_snp)
                    call s_amr_fw_szi(amr_fw_snxp, amr_fw_snp); call s_amr_fw_szi(amr_fw_sqbase, amr_fw_snp)
                    amr_fw_map(r) = amr_fw_snp
                    amr_fw_sprank(amr_fw_snp) = r
                end if
                cnt = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                amr_fw_snx = amr_fw_snx + 1
                call s_amr_fw_szi(amr_fw_sblk, amr_fw_snx); call s_amr_fw_szi3(amr_fw_sbl, amr_fw_snx)
                call s_amr_fw_szi3(amr_fw_sbh, amr_fw_snx); call s_amr_fw_szi(amr_fw_spi, amr_fw_snx)
                call s_amr_fw_szi(amr_fw_sqo, amr_fw_snx); call s_amr_fw_szi(amr_fw_spo, amr_fw_snx)
                amr_fw_sblk(amr_fw_snx) = k
                amr_fw_sbl(:,amr_fw_snx) = bl; amr_fw_sbh(:,amr_fw_snx) = bh
                amr_fw_spo(amr_fw_snx) = cnt
                amr_fw_spi(amr_fw_snx) = amr_fw_map(r)
                amr_fw_sqo(amr_fw_snx) = amr_fw_pq(r) + amr_fw_nx(r)*XA_NH
                amr_fw_pq(r) = amr_fw_pq(r) + cnt
                amr_fw_nx(r) = amr_fw_nx(r) + 1
            end do
        end do
        qbase = 0
        do ip = 1, amr_fw_snp
            r = amr_fw_sprank(ip)
            amr_fw_snxp(ip) = amr_fw_nx(r)
            amr_fw_sqsz(ip) = amr_fw_pq(r) + amr_fw_nx(r)*XA_NH
            amr_fw_sqbase(ip) = qbase; qbase = qbase + amr_fw_sqsz(ip)
            amr_fw_map(r) = 0; amr_fw_nx(r) = 0; amr_fw_pq(r) = 0
        end do
        call s_amr_fw_szr(amr_fw_sq, qbase)
        ! recv plan (coarse-owner side): my interior x region(k) over level-1 blocks I do not own
        amr_fw_rnx = 0; amr_fw_rnp = 0
        do k = 1, amr_num_blocks
            if (amr_block_level(k) /= 1) cycle
            owner = amr_block_owner(k)
            if (owner == proc_rank) cycle
            rlo = 0; rhi = 0
            rlo(1) = amr_region_lo_all(1, k); rhi(1) = amr_region_hi_all(1, k)
            if (n_glb > 0) then; rlo(2) = amr_region_lo_all(2, k); rhi(2) = amr_region_hi_all(2, k); end if
            if (p_glb > 0) then; rlo(3) = amr_region_lo_all(3, k); rhi(3) = amr_region_hi_all(3, k); end if
            call s_amr_box_isect(rlo, rhi, milo, mihi, bl, bh)
            if (bl(1) > bh(1) .or. bl(2) > bh(2) .or. bl(3) > bh(3)) cycle
            if (amr_fw_map(owner) == 0) then
                amr_fw_rnp = amr_fw_rnp + 1
                call s_amr_fw_szi(amr_fw_rprank, amr_fw_rnp); call s_amr_fw_szi(amr_fw_rqsz, amr_fw_rnp)
                call s_amr_fw_szi(amr_fw_rnxp, amr_fw_rnp); call s_amr_fw_szi(amr_fw_rqbase, amr_fw_rnp)
                amr_fw_map(owner) = amr_fw_rnp
                amr_fw_rprank(amr_fw_rnp) = owner
            end if
            cnt = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
            amr_fw_rnx = amr_fw_rnx + 1
            call s_amr_fw_szi(amr_fw_rblk, amr_fw_rnx); call s_amr_fw_szi3(amr_fw_rbl, amr_fw_rnx)
            call s_amr_fw_szi3(amr_fw_rbh, amr_fw_rnx); call s_amr_fw_szi(amr_fw_rpi, amr_fw_rnx)
            call s_amr_fw_szi(amr_fw_rqo, amr_fw_rnx); call s_amr_fw_szi(amr_fw_rpo, amr_fw_rnx)
            amr_fw_rblk(amr_fw_rnx) = k
            amr_fw_rbl(:,amr_fw_rnx) = bl; amr_fw_rbh(:,amr_fw_rnx) = bh
            amr_fw_rpo(amr_fw_rnx) = cnt
            amr_fw_rpi(amr_fw_rnx) = amr_fw_map(owner)
            amr_fw_rqo(amr_fw_rnx) = amr_fw_pq(owner) + amr_fw_nx(owner)*XA_NH
            amr_fw_pq(owner) = amr_fw_pq(owner) + cnt
            amr_fw_nx(owner) = amr_fw_nx(owner) + 1
        end do
        qbase = 0
        do ip = 1, amr_fw_rnp
            r = amr_fw_rprank(ip)
            amr_fw_rnxp(ip) = amr_fw_nx(r)
            amr_fw_rqsz(ip) = amr_fw_pq(r) + amr_fw_nx(r)*XA_NH
            amr_fw_rqbase(ip) = qbase; qbase = qbase + amr_fw_rqsz(ip)
            amr_fw_map(r) = 0; amr_fw_nx(r) = 0; amr_fw_pq(r) = 0
        end do
        call s_amr_fw_szr(amr_fw_rq, qbase)
        nreq = amr_fw_snp + amr_fw_rnp
        call s_amr_fw_szi(amr_fw_req, max(nreq, 1)); call s_amr_fw_szi(amr_fw_reqw, max(nreq, 1))
        nreq = 0
        do ip = 1, amr_fw_rnp
            call s_xa_rec(XA_F7W_RCV, 2, amr_fw_rqsz(ip) - amr_fw_rnxp(ip)*XA_NH, tq)
            nreq = nreq + 1; amr_fw_reqw(nreq) = amr_fw_rqsz(ip)
            call MPI_IRECV(amr_fw_rq(amr_fw_rqbase(ip) + 1), amr_fw_rqsz(ip), mpi_p, amr_fw_rprank(ip), tq, MPI_COMM_WORLD, &
                           & amr_fw_req(nreq), ierr)
        end do
        ! owner-local covered overwrites + device packs, grouped per owned block: amr_rvw is a single device mirror, so a
        ! block's (cyl_coord) radii push must immediately precede that block's overwrite/pack kernels; the transfer list is
        ! k-grouped by construction, so a monotone cursor drains each block's sends inside its group
        cur = 1
        do k = 1, amr_num_blocks
            if (amr_block_level(k) /= 1) cycle
            if (amr_block_owner(k) /= proc_rank) cycle
            rr = amr_slots(k)%amr_ref_ratio
            nchild = rr; if (n_glb > 0) nchild = nchild*rr; if (p_glb > 0) nchild = nchild*rr
            dj_hi = merge(rr - 1, 0, n_glb > 0); dk_hi = merge(rr - 1, 0, p_glb > 0)
            if (cyl_coord) then
                amr_rvw(0:amr_slots(k)%n) = amr_slots(k)%y_cc(0:amr_slots(k)%n)
                $:GPU_UPDATE(device='[amr_rvw]')
            end if
            rlo = 0; rhi = 0
            rlo(1) = amr_region_lo_all(1, k); rhi(1) = amr_region_hi_all(1, k)
            if (n_glb > 0) then; rlo(2) = amr_region_lo_all(2, k); rhi(2) = amr_region_hi_all(2, k); end if
            if (p_glb > 0) then; rlo(3) = amr_region_lo_all(3, k); rhi(3) = amr_region_hi_all(3, k); end if
            call s_amr_box_isect(rlo, rhi, milo, mihi, bl, bh)
            if (bl(1) <= bh(1) .and. bl(2) <= bh(2) .and. bl(3) <= bh(3)) call s_amr_restrict_overwrite_device_sf(coarse_tgt, &
                & amr_loc_of(k), bl, bh, o1, o2, o3, rlo, rr, dj_hi, dk_hi, nchild)
            do while (cur <= amr_fw_snx)
                if (amr_fw_sblk(cur) /= k) exit
                cnt = amr_fw_spo(cur)
                boff = amr_fw_sqbase(amr_fw_spi(cur)) + amr_fw_sqo(cur)
                call s_amr_restrict_pack_device(amr_loc_of(k), amr_fw_sbl(:,cur), amr_fw_sbh(:,cur), rlo, rr, dj_hi, dk_hi, &
                                                & nchild, amr_fw_sq(boff + XA_NH + 1:boff + XA_NH + cnt))
                if (XA_NH > 0) call s_xa_hdr_pack(amr_fw_sq(boff + 1:boff + XA_NH), XA_F7W_SND, k, amr_fw_sbl(:,cur), &
                    & amr_fw_sbh(:,cur))
                cur = cur + 1
            end do
        end do
        do ip = 1, amr_fw_snp
            call s_xa_rec(XA_F7W_SND, 1, amr_fw_sqsz(ip) - amr_fw_snxp(ip)*XA_NH, tq)
            nreq = nreq + 1; amr_fw_reqw(nreq) = -1
            call MPI_ISEND(amr_fw_sq(amr_fw_sqbase(ip) + 1), amr_fw_sqsz(ip), mpi_p, amr_fw_sprank(ip), tq, MPI_COMM_WORLD, &
                           & amr_fw_req(nreq), ierr)
        end do
        if (nreq > 0) then
#ifdef MFC_DEBUG
            block
                integer :: st(MPI_STATUS_SIZE, nreq), gotw, q
                call MPI_WAITALL(nreq, amr_fw_req, st, ierr)
                do q = 1, nreq
                    if (amr_fw_reqw(q) < 0) cycle
                    call MPI_GET_COUNT(st(:,q), mpi_p, gotw, ierr)
                    @:ASSERT(gotw == amr_fw_reqw(q), "restrict L1 wave: received message length differs from the plan")
                end do
            end block
#else
            call MPI_WAITALL(nreq, amr_fw_req, MPI_STATUSES_IGNORE, ierr)
#endif
        end if
        do idx = 1, amr_fw_rnx
            cnt = amr_fw_rpo(idx)
            boff = amr_fw_rqbase(amr_fw_rpi(idx)) + amr_fw_rqo(idx)
            if (XA_NH > 0) call s_xa_hdr_check(amr_fw_rq(boff + 1:boff + XA_NH), XA_F7W_SND, amr_fw_rblk(idx), amr_fw_rbl(:,idx), &
                & amr_fw_rbh(:,idx))
            ! DEVICE unpack of the covered box only (the strided-update flang trap - see s_restrict_fine_to_coarse)
            call s_l0_pack_unpack_block_sf(coarse_tgt, amr_fw_rbl(1, idx) - o1, amr_fw_rbl(2, idx) - o2, amr_fw_rbl(3, idx) - o3, &
                                           & amr_fw_rbh(1, idx) - amr_fw_rbl(1, idx), amr_fw_rbh(2, idx) - amr_fw_rbl(2, idx), &
                                           & amr_fw_rbh(3, idx) - amr_fw_rbl(3, idx), &
                                           & amr_fw_rq(boff + XA_NH + 1:boff + XA_NH + cnt), .false.)
        end do
#endif

    end subroutine s_amr_restrict_l1_wave

    !> Deliver the current level>=2 block's fine flux registers to its PARENT block's owner, which holds the matching creg and
    !! applies the correction. One blocking send/recv pair per dimension, mirroring s_amr_p2p_reflux_faces:
    !! freg(d)%lo/hi(:,:,:,slot) is CONTIGUOUS (trailing slot index fixed, leading dims full), so it goes on the wire with no pack
    !! and its GPU_UPDATE is a contiguous transfer. Several remote children of one parent reuse these tags, which is safe because
    !! MPI does not overtake between a fixed (source, tag, comm) triple and both owners walk the sibling loop in the same replicated
    !! block order. Tag base is disjoint from s_amr_p2p_reflux_faces so an L0/L1 delivery can never be mistaken for a parent
    !! delivery.
    impure subroutine s_amr_p2p_freg_to_parent(pblk)

        integer, intent(in) :: pblk

#ifdef MFC_MPI
        integer :: cowner, powner, cnt, ierr

        cowner = amr_block_owner(amr_cur)
        powner = amr_block_owner(pblk)
        if (proc_rank == cowner) then
            #:for D in [1, 2, 3]
                if (${D}$ <= num_dims) then
                    cnt = size(freg(${D}$)%lo(:,:,:,amr_reg_cur))
                    $:GPU_UPDATE(host='[freg(' + str(D) + ')%lo(:, :, :, amr_reg_cur), freg(' + str(D) &
                                 & + ')%hi(:, :, :, amr_reg_cur)]')
                    call s_xa_rec(XA_F5_FREG_SND, 1, cnt, ${40 + 2*D}$)
                    call MPI_SEND(freg(${D}$)%lo(:,:,:,amr_reg_cur), cnt, mpi_p, powner, ${40 + 2*D}$, MPI_COMM_WORLD, ierr)
                    call s_xa_rec(XA_F5_FREG_SND, 1, cnt, ${41 + 2*D}$)
                    call MPI_SEND(freg(${D}$)%hi(:,:,:,amr_reg_cur), cnt, mpi_p, powner, ${41 + 2*D}$, MPI_COMM_WORLD, ierr)
                end if
            #:endfor
        else
            #:for D in [1, 2, 3]
                if (${D}$ <= num_dims) then
                    cnt = size(freg(${D}$)%lo(:,:,:,amr_reg_cur))
                    call s_xa_rec(XA_F5_FREG_RCV, 2, cnt, ${40 + 2*D}$)
                    call MPI_RECV(freg(${D}$)%lo(:,:,:,amr_reg_cur), cnt, mpi_p, cowner, ${40 + 2*D}$, MPI_COMM_WORLD, &
                                  & MPI_STATUS_IGNORE, ierr)
                    call s_xa_rec(XA_F5_FREG_RCV, 2, cnt, ${41 + 2*D}$)
                    call MPI_RECV(freg(${D}$)%hi(:,:,:,amr_reg_cur), cnt, mpi_p, cowner, ${41 + 2*D}$, MPI_COMM_WORLD, &
                                  & MPI_STATUS_IGNORE, ierr)
                    $:GPU_UPDATE(device='[freg(' + str(D) + ')%lo(:, :, :, amr_reg_cur), freg(' + str(D) &
                                 & + ')%hi(:, :, :, amr_reg_cur)]')
                end if
            #:endfor
        end if
#endif

    end subroutine s_amr_p2p_freg_to_parent

    !> Sibling-seam face weights for level>=2 block kb under parent pblk: 0 on a face shared with a same-parent sibling tile
    !! (fine-fine, not c/f - refluxing there double-writes and leaks; the outside parent cell is covered by the sibling's restrict),
    !! 1 otherwise. REPLICATED metadata only (f_amr_parent_block + f_amr_seam read amr_region_*_all), so every rank derives the same
    !! weights - the reflux apply AND the freg wave's wire-skip both call this, so what ships and what is consumed cannot drift
    !! apart.
    impure subroutine s_amr_sibling_face_weights(kb, pblk, w_lo, w_hi)

        integer, intent(in)   :: kb, pblk
        real(wp), intent(out) :: w_lo(3), w_hi(3)
        integer               :: y, d

        w_lo = 1._wp; w_hi = 1._wp
        do y = 1, amr_num_blocks  ! block outer, dim inner: f_amr_parent_block (a linear scan) is evaluated once per sibling
            if (y == kb) cycle
            if (f_amr_parent_block(y) /= pblk) cycle  ! same-parent sibling tile only (guarantees same level)
            do d = 1, num_dims
                if (f_amr_seam(kb, y, d)) w_hi(d) = 0._wp  ! sibling just above -> shared high face
                if (f_amr_seam(y, kb, d)) w_lo(d) = 0._wp  ! sibling just below -> shared low face
            end do
        end do

    end subroutine s_amr_sibling_face_weights

    !> Multi-level reflux: apply the Berger-Colella C/F flux correction from the current level>=2 block into its PARENT block's
    !! cells just OUTSIDE the block footprint, in the parent-fine frame (mirror of the L0 s_amr_apply_reflux targeted at the parent
    !! - "the coarse" is level l-1). State form: q_parent(outside) += dt*(F_coarse - Fbar_fine)/dxf on the low face and +=
    !! dt*(Fbar_fine - F_coarse)/dxf on the high face, where Fbar_fine is the child-averaged fine register. creg/freg key off this
    !! block's slot. Per-face parent-fine dx (stretched-grid safe).
    !!
    !! The PARENT's owner applies: it holds the parent field and the parent-side creg (captured over its own advance). Only freg
    !! crosses the wire, and only when the two owners differ - under tower co-location they never do, so this is byte-identical to
    !! the previous owner-local form. BOTH participants must reach this routine or the P2P pair deadlocks (cf. the restrict).
    impure subroutine s_amr_reflux_to_parent(dt_reflux, do_xchg)

        real(wp), intent(in) :: dt_reflux
        !> exchange the split-ownership freg here (the subcycle per-box path); the lock-step driver runs s_amr_freg_wave first
        logical, intent(in) :: do_xchg
        integer             :: pblk, d, olo(3), ohi(3), glo(3), ghi(3), woff(3), plo(3), phi(3)
        real(wp)            :: w_lo(3), w_hi(3), mlo(3), mhi(3)
        logical             :: own_child, own_parent

        pblk = f_amr_parent_block(amr_cur)
        own_child = amr_rank_owns_block
        own_parent = (amr_block_owner(pblk) == proc_rank)
        if (.not. (own_child .or. own_parent)) return
        if (do_xchg .and. (own_child .neqv. own_parent)) call s_amr_p2p_freg_to_parent(pblk)
        if (.not. own_parent) return
        ! max_grid_size tiling of a level>=2 feature: a face shared with an ADJACENT sibling tile (same parent) is fine-fine, not a
        ! c/f boundary - its "outside" parent cell is covered by the sibling's restrict, so refluxing there double-writes and leaks.
        ! Skip those faces (weight 0); the fine-fine halo already matched the shared seam flux. No siblings -> all weights 1
        ! (no-op).
        call s_amr_sibling_face_weights(amr_cur, pblk, w_lo, w_hi)
        ! parent-fine frame for the shared reflux kernel: outside cell = isect boundary +/-1; creg-local loop range 0:extent;
        ! transverse write at the isect origin. Per-face parent-fine cell widths - dx at the low/high OUTSIDE cell (olo/ohi),
        ! mirroring the L0/L1 s_amr_apply_reflux_state so a stretched parent grid corrects each C/F face with its own width (on a
        ! uniform grid dx is constant, so this is byte-identical to the previous single-dxf form).
        ! Footprint from REPLICATED metadata (s_amr_parent_foot), not amr_isect_lo/hi: on the parent's owner the child's own isect
        ! is
        ! the empty non-owner sentinel whenever the two differ. Identical box while co-located. rr likewise comes from the global
        ! amr_ref_ratio rather than amr_slots(amr_cur), whose slot need not be allocated on this rank.
        call s_amr_parent_foot(amr_cur, pblk, plo, phi)
        olo = 0; ohi = 0; glo = 0; ghi = 0; woff = 0; mlo = 1._wp; mhi = 1._wp
        do d = 1, num_dims
            olo(d) = plo(d) - 1; ohi(d) = phi(d) + 1
            ghi(d) = phi(d) - plo(d)
            woff(d) = plo(d)
        end do
        mlo(1) = amr_slots(pblk)%dx(olo(1)); mhi(1) = amr_slots(pblk)%dx(ohi(1))
        if (n_glb > 0) then; mlo(2) = amr_slots(pblk)%dy(olo(2)); mhi(2) = amr_slots(pblk)%dy(ohi(2)); end if
        if (p_glb > 0) then; mlo(3) = amr_slots(pblk)%dz(olo(3)); mhi(3) = amr_slots(pblk)%dz(ohi(3)); end if
        call s_amr_br_load(amr_loc_of(pblk))
        call s_amr_reflux_apply_faces(amr_cons_br, amr_reg_cur, amr_ref_ratio, dt_reflux, olo, ohi, glo, ghi, woff, w_lo, w_hi, &
                                      & mlo, mhi)
        call s_amr_br_store(amr_loc_of(pblk))

    end subroutine s_amr_reflux_to_parent

    !> Rank r's coarse INTERIOR box (global) from the computed decomposition (s_amr_rank_decomp, no ghosts). Covered coarse cells
    !! are in-domain, so restriction targets are identified by interior overlap alone.
    pure subroutine s_amr_rank_interior(r, ilo, ihi)

        integer, intent(in)  :: r
        integer, intent(out) :: ilo(3), ihi(3)
        integer              :: sidx(3), ext(3)

        call s_amr_rank_decomp(r, sidx, ext)
        ilo = 0; ihi = 0
        ilo(1) = sidx(1); ihi(1) = sidx(1) + ext(1)
        if (n_glb > 0) then; ilo(2) = sidx(2); ihi(2) = sidx(2) + ext(2); end if
        if (p_glb > 0) then; ilo(3) = sidx(3); ihi(3) = sidx(3) + ext(3); end if

    end subroutine s_amr_rank_interior

    !> Device-native restriction overwrite: restrict the fine block q_fine (DEVICE) to coarse averages over the covered coarse cells
    !! [bl:bh] GLOBAL and write coarse_tgt (DEVICE) directly - no host round-trip, only the covered cells touched (the old
    !! whole-coarse device push clobbered non-covered cells). Same child-sum order (ddk, ddj, then fi0 and fi0+1; /nchild; stp cast)
    !! as the old host restrict path, so bit-identical to it on CPU and matches the coarse restriction. q_fine (== the flat store)
    !! and coarse_tgt are device-resident. TWIN: s_amr_restrict_pack_device runs this same child-sum into a wire buffer - any change
    !! to the loop order, arithmetic, or casts here must be mirrored there byte-identically (owner-local and scattered coarse cells
    !! must match bit-for-bit). TWIN(q<->pb/mv) s_amr_restrict_pbmv_box_device runs this same child-sum on pb/mv - keep the stencil
    !! lockstep. Two targets, one body: the coarse destination is the level-0 monolithic field (`_sf`) or a parent BLOCK in the flat
    !! store (`_st`); the fine source is always a block, so it is always the store.
    #:for SFX, CT in [('sf', ''), ('st', 'amr_cons_st')]
        #:set CW = (lambda ix: CT + '(ci - o1, cj - o2, ck - o3, ' + ix + ', ctloc)') if CT else (lambda ix: 'coarse_tgt(' + ix &
                    & + ')%sf(ci - o1, cj - o2, ck - o3)')
        impure subroutine s_amr_restrict_overwrite_device_${SFX}$(${'ctloc' if CT else 'coarse_tgt'}$, loc, bl, bh, o1, o2, o3, &
            & rlo, rr, dj_hi, dk_hi, nchild)

            #:if CT
                integer, intent(in) :: ctloc
            #:else
                type(scalar_field), dimension(sys_size), intent(inout) :: coarse_tgt
            #:endif
            integer, intent(in) :: loc
            integer, intent(in) :: bl(3), bh(3), o1, o2, o3, rlo(3), rr, dj_hi, dk_hi, nchild
            integer             :: i, ci, cj, ck, fi0, fj0, fk0, ddi, ddj, ddk, bl1, bl2, bl3, bh1, bh2, bh3, rl1, rl2, rl3
            real(wp)            :: acc, wacc, w

            bl1 = bl(1); bl2 = bl(2); bl3 = bl(3); bh1 = bh(1); bh2 = bh(2); bh3 = bh(3)
            rl1 = rlo(1); rl2 = rlo(2); rl3 = rlo(3)
            if (cyl_coord) then
                ! axisymmetric volume-weighted fold-back: weight each fine child by its cell-center radius (amr_rvw = fine y_cc, on
                ! device). Same child order as the Cartesian path and the scatter pack, so CPU==GPU and np=1==np>=2.
                $:GPU_PARALLEL_LOOP(collapse=4, private='[fi0, fj0, fk0, ddi, ddj, ddk, acc, wacc, w]')
                do i = 1, sys_size
                    do ck = bl3, bh3
                        do cj = bl2, bh2
                            do ci = bl1, bh1
                                fi0 = (ci - rl1)*rr; fj0 = (cj - rl2)*rr; fk0 = (ck - rl3)*rr
                                acc = 0._wp; wacc = 0._wp
                                do ddk = 0, dk_hi
                                    do ddj = 0, dj_hi
                                        w = amr_rvw(fj0 + ddj)
                                        do ddi = 0, rr - 1
                                            acc = acc + real(amr_cons_st(fi0 + ddi, fj0 + ddj, fk0 + ddk, i, loc), wp)*w
                                            wacc = wacc + w
                                        end do
                                    end do
                                end do
                                ${CW('i')}$ = real(acc/wacc, stp)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
                return
            end if
            $:GPU_PARALLEL_LOOP(collapse=4, private='[fi0, fj0, fk0, ddi, ddj, ddk, acc]')
            do i = 1, sys_size
                do ck = bl3, bh3
                    do cj = bl2, bh2
                        do ci = bl1, bh1
                            fi0 = (ci - rl1)*rr; fj0 = (cj - rl2)*rr; fk0 = (ck - rl3)*rr
                            acc = 0._wp
                            do ddk = 0, dk_hi
                                do ddj = 0, dj_hi
                                    do ddi = 0, rr - 1
                                        acc = acc + real(amr_cons_st(fi0 + ddi, fj0 + ddj, fk0 + ddk, i, loc), wp)
                                    end do
                                end do
                            end do
                            ${CW('i')}$ = real(acc/real(nchild, wp), stp)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

        end subroutine s_amr_restrict_overwrite_device_${SFX}$
    #:endfor

    !> Device pack of one destination's covered restrict slice (np>=2 scatter): restrict the fine block q_fine (DEVICE) over the
    !! covered coarse box [bl:bh] GLOBAL straight into the contiguous wire buffer buf (host, via copyout) - only the slice crosses
    !! PCIe, not the full fine field. Same child-sum order and wp values as s_amr_restrict_overwrite_device (no stp cast: the wire
    !! carries wp and the receiver casts), packed with ci fastest, then cj, ck, i, matching the receiver's sequential unpack. TWIN:
    !! s_amr_restrict_overwrite_device runs this same child-sum in place - any change to the loop order, arithmetic, or casts here
    !! must be mirrored there byte-identically (owner-local and scattered coarse cells must match bit-for-bit). TWIN(q<->pb/mv)
    !! s_amr_restrict_pbmv_pack_device runs this same child-sum into a wire buffer - keep lockstep.
    impure subroutine s_amr_restrict_pack_device(loc, bl, bh, rlo, rr, dj_hi, dk_hi, nchild, buf)

        integer, intent(in) :: loc
        integer, intent(in) :: bl(3), bh(3), rlo(3), rr, dj_hi, dk_hi, nchild
        real(wp), intent(inout), contiguous :: buf(:)
        integer :: i, ci, cj, ck, fi0, fj0, fk0, ddi, ddj, ddk, bl1, bl2, bl3, bh1, bh2, bh3, rl1, rl2, rl3, n1, n2, n3
        real(wp) :: acc, wacc, w

        bl1 = bl(1); bl2 = bl(2); bl3 = bl(3); bh1 = bh(1); bh2 = bh(2); bh3 = bh(3)
        rl1 = rlo(1); rl2 = rlo(2); rl3 = rlo(3)
        n1 = bh1 - bl1 + 1; n2 = bh2 - bl2 + 1; n3 = bh3 - bl3 + 1
        if (cyl_coord) then
            ! axisymmetric volume-weighted pack (amr_rvw = fine y_cc, on device); same child order as the overwrite kernel
            $:GPU_PARALLEL_LOOP(collapse=4, private='[fi0, fj0, fk0, ddi, ddj, ddk, acc, wacc, w]', copyout='[buf]')
            do i = 1, sys_size
                do ck = bl3, bh3
                    do cj = bl2, bh2
                        do ci = bl1, bh1
                            fi0 = (ci - rl1)*rr; fj0 = (cj - rl2)*rr; fk0 = (ck - rl3)*rr
                            acc = 0._wp; wacc = 0._wp
                            do ddk = 0, dk_hi
                                do ddj = 0, dj_hi
                                    w = amr_rvw(fj0 + ddj)
                                    do ddi = 0, rr - 1
                                        acc = acc + real(amr_cons_st(fi0 + ddi, fj0 + ddj, fk0 + ddk, i, loc), wp)*w
                                        wacc = wacc + w
                                    end do
                                end do
                            end do
                            buf(1 + (ci - bl1) + n1*((cj - bl2) + n2*((ck - bl3) + n3*(i - 1)))) = acc/wacc
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            return
        end if
        $:GPU_PARALLEL_LOOP(collapse=4, private='[fi0, fj0, fk0, ddi, ddj, ddk, acc]', copyout='[buf]')
        do i = 1, sys_size
            do ck = bl3, bh3
                do cj = bl2, bh2
                    do ci = bl1, bh1
                        fi0 = (ci - rl1)*rr; fj0 = (cj - rl2)*rr; fk0 = (ck - rl3)*rr
                        acc = 0._wp
                        do ddk = 0, dk_hi
                            do ddj = 0, dj_hi
                                do ddi = 0, rr - 1
                                    acc = acc + real(amr_cons_st(fi0 + ddi, fj0 + ddj, fk0 + ddk, i, loc), wp)
                                end do
                            end do
                        end do
                        buf(1 + (ci - bl1) + n1*((cj - bl2) + n2*((ck - bl3) + n3*(i - 1)))) = acc/real(nchild, wp)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_restrict_pack_device

    !> Apply phase-change relaxation (relax) to the current fine block's interior, BEFORE restriction. Relaxation is a cell-local,
    !! mass/energy-conserving equilibration (no stencil, no ghosts), so it needs no coarse/fine coupling: it just runs over the fine
    !! interior. Swaps m/n/p to the fine extents so s_infinite_relaxation_k's 0:m,0:n,0:p loop covers this block. Matches the coarse
    !! timing (once per full step; the coarse relax runs once after s_tvd_rk on q_cons_ts(1)) but on the fine solution so the fine
    !! dynamics equilibrate at fine resolution, not only the restricted coarse average.
    impure subroutine s_amr_relax_fine()

        if (.not. amr_rank_owns_block) return
        call s_amr_swap_to_fine()
        call s_amr_br_load(amr_loc_of(amr_cur))
        call s_infinite_relaxation_k(amr_cons_br)
        call s_amr_br_store(amr_loc_of(amr_cur))
        call s_amr_restore_coarse()

    end subroutine s_amr_relax_fine

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
        ! A moving body must stay inside its block (a body overlapping the block edge gets silently clipped forcing - abort
        ! instead).
        ! Under dynamic regrid the expansion contained it with margin max(amr_buf,4) and we require body + image-point stencil reach
        ! (2 coarse cells) to remain contained between regrids; on a STATIC block the user's placement is authoritative (validated
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

    !> Non-polytropic QBMM: piecewise-constant prolongation of the block's pb/mv interior from the gathered coarse side-state
    !! amr_cg_pb/mv (patch-local frame; the callers run s_amr_gather_coarse_patch_pbmv on ALL ranks first, so np>=2 couples to the
    !! correct coarse rank). HOST loops + device push; the gather is host-current (.not. pull_host). TWIN
    !! s_interpolate_coarse_to_fine (pb/mv<->q): pb/mv is piecewise-constant where q_cons is minmod-limited, but the child-offset
    !! frame and realizability/closure policy track it; keep those lockstep.
    impure subroutine s_amr_prolong_pbmv()

        integer :: fi, fj, fk, q, ib_, ci, cj, ck, rr, lo1, lo2, lo3, ox, oy, oz

        ! HOST prolongation (both call paths make the coarse pb/mv host mirrors current first: init writes them on the host, regrid
        ! refreshes them from the device); the device copy of the fine side-state is pushed at the end

        ! coarse pb/mv are read from the gathered block-local patch amr_cg_pb/mv (fine-level distribution): patch-local frame, cell
        ! 0
        ! == amr_cpat_off (matching s_prolong_one_var). The gather is a host loop (.not. pull_host) done by the callers.

        ox = amr_cpat_off(1); oy = amr_cpat_off(2); oz = amr_cpat_off(3)
        rr = amr_slots(amr_cur)%amr_ref_ratio
        lo1 = amr_isect_lo(1); lo2 = amr_isect_lo(2); lo3 = amr_isect_lo(3)
        do ib_ = 1, nb
            do q = 1, nnode
                do fk = 0, amr_slots(amr_cur)%p
                    ck = 0
                    if (p_glb > 0) ck = lo3 + fk/rr - oz
                    do fj = 0, amr_slots(amr_cur)%n
                        cj = 0
                        if (n_glb > 0) cj = lo2 + fj/rr - oy
                        do fi = 0, amr_slots(amr_cur)%m
                            ci = lo1 + fi/rr - ox
                            amr_slots(amr_cur)%pb_f%sf(fi, fj, fk, q, ib_) = amr_cg_pb(ci, cj, ck, q, ib_)
                            amr_slots(amr_cur)%mv_f%sf(fi, fj, fk, q, ib_) = amr_cg_mv(ci, cj, ck, q, ib_)
                        end do
                    end do
                end do
            end do
        end do
        $:GPU_UPDATE(device='[amr_slots(amr_cur)%pb_f%sf, amr_slots(amr_cur)%mv_f%sf]')

    end subroutine s_amr_prolong_pbmv

    !> Non-polytropic QBMM: piecewise-constant prolongation of the fine pb/mv GHOST shell from the coarse side-state (device kernel;
    !! interior untouched). The ghosts feed the widened-idwint conversions and the qbmm rhs over the shell, mirroring the q_cons
    !! ghost fill. All four arrays are assumed-shape dummies with %sf pointer-member actuals (the device-proven pb_ts pattern); only
    !! raw derived-type 5D members as actuals tripped nvfortran's component-section data clauses on device. TWIN
    !! s_amr_fill_fine_ghosts (pb/mv<->q): q_cons sibling; keep the ghost-fill mapping lockstep.
    impure subroutine s_amr_fill_fine_ghosts_pbmv(pb_c, mv_c, pb_t, mv_t)

        !> coarse pb/mv read from the gathered block-local patch amr_cg_pb/mv (0-based patch frame, cell 0 == amr_cpat_off): the
        !! callers run s_amr_gather_coarse_patch_pbmv on ALL ranks first, so np>=2 reads the correct coarse rank's side-state
        real(stp), dimension(0:,0:,0:,1:,1:), intent(in) :: pb_c, mv_c

        real(stp), dimension(amr_slots(amr_cur)%idwbuff(1)%beg:,amr_slots(amr_cur)%idwbuff(2)%beg:, &
             & amr_slots(amr_cur)%idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_t, mv_t
        integer               :: fi, fj, fk, q, ib_, ci, cj, ck, rr, lo1, lo2, lo3, ox, oy, oz
        integer               :: s, ns, l1, u1, l2, u2, l3, u3, ss, g, r, n1, n2, stot
        integer, dimension(6) :: sb1, se1, sb2, se2, sb3, se3, soff, scnt
        logical               :: d2, d3

        ox = amr_cpat_off(1); oy = amr_cpat_off(2); oz = amr_cpat_off(3)
        d2 = n_glb > 0; d3 = p_glb > 0
        rr = amr_slots(amr_cur)%amr_ref_ratio
        lo1 = amr_isect_lo(1); lo2 = amr_isect_lo(2); lo3 = amr_isect_lo(3)
        call s_amr_build_ghost_slabs(ns, sb1, se1, sb2, se2, sb3, se3)
        ! flat index over the concatenated DISJOINT slabs - one kernel instead of ns; see s_amr_fill_fine_ghosts
        soff(1) = 0
        do s = 1, ns
            scnt(s) = (se1(s) - sb1(s) + 1)*(se2(s) - sb2(s) + 1)*(se3(s) - sb3(s) + 1)
            if (s < ns) soff(s + 1) = soff(s) + scnt(s)
        end do
        stot = soff(ns) + scnt(ns)
        $:GPU_PARALLEL_LOOP(collapse=3, copyin='[sb1, se1, sb2, se2, sb3, se3, soff, scnt]', private='[s, ss, r, n1, n2, fi, fj, &
                            & fk, ci, cj, ck]')
        do ib_ = 1, nb
            do q = 1, nnode
                do g = 0, stot - 1
                    s = 1
                    do ss = 2, ns
                        if (g >= soff(ss)) s = ss
                    end do
                    r = g - soff(s)
                    n1 = se1(s) - sb1(s) + 1; n2 = se2(s) - sb2(s) + 1
                    fi = sb1(s) + mod(r, n1)
                    fj = sb2(s) + mod(r/n1, n2)
                    fk = sb3(s) + r/(n1*n2)
                    ck = 0
                    if (d3) ck = lo3 + floor(real(fk, wp)/real(rr, wp)) - oz
                    cj = 0
                    if (d2) cj = lo2 + floor(real(fj, wp)/real(rr, wp)) - oy
                    ci = lo1 + floor(real(fi, wp)/real(rr, wp)) - ox
                    pb_t(fi, fj, fk, q, ib_) = pb_c(ci, cj, ck, q, ib_)
                    mv_t(fi, fj, fk, q, ib_) = mv_c(ci, cj, ck, q, ib_)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_fill_fine_ghosts_pbmv

    !> Non-polytropic QBMM: device copy of the block's pb/mv into the step-entry backup (SSP-RK). TWIN s_amr_copy_fine_fields
    !! (pb/mv<->q): q_cons sibling of this step-entry backup; keep lockstep.
    impure subroutine s_amr_backup_pbmv(pb_s, mv_s, pb_d, mv_d)

        real(stp), dimension(amr_slots(amr_cur)%idwbuff(1)%beg:,amr_slots(amr_cur)%idwbuff(2)%beg:, &
             & amr_slots(amr_cur)%idwbuff(3)%beg:,1:,1:), intent(in) :: pb_s, mv_s
        real(stp), dimension(amr_slots(amr_cur)%idwbuff(1)%beg:,amr_slots(amr_cur)%idwbuff(2)%beg:, &
             & amr_slots(amr_cur)%idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_d, mv_d
        integer :: fi, fj, fk, q, ib_, b1, e1, b2, e2, b3, e3

        b1 = amr_slots(amr_cur)%idwbuff(1)%beg; e1 = amr_slots(amr_cur)%idwbuff(1)%end
        b2 = amr_slots(amr_cur)%idwbuff(2)%beg; e2 = amr_slots(amr_cur)%idwbuff(2)%end
        b3 = amr_slots(amr_cur)%idwbuff(3)%beg; e3 = amr_slots(amr_cur)%idwbuff(3)%end
        $:GPU_PARALLEL_LOOP(collapse=5)
        do ib_ = 1, nb
            do q = 1, nnode
                do fk = b3, e3
                    do fj = b2, e2
                        do fi = b1, e1
                            pb_d(fi, fj, fk, q, ib_) = pb_s(fi, fj, fk, q, ib_)
                            mv_d(fi, fj, fk, q, ib_) = mv_s(fi, fj, fk, q, ib_)
                        end do
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_backup_pbmv

    !> Non-polytropic QBMM: SSP-RK stage update of the block's pb/mv (device kernel, interior only; mirror of the coarse pb_ts/mv_ts
    !! stage combination in m_time_steppers). TWIN s_amr_fine_rk_update + s_tvd_rk (m_time_steppers): same SSP-RK stage combination
    !! on pb/mv; keep all three lockstep.
    impure subroutine s_amr_fine_rk_update_pbmv(pb_u, mv_u, pb_s, mv_s, rpb, rmv, c1, c2, c3, c4, dtl)

        real(stp), dimension(amr_slots(amr_cur)%idwbuff(1)%beg:,amr_slots(amr_cur)%idwbuff(2)%beg:, &
             & amr_slots(amr_cur)%idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_u, mv_u
        real(stp), dimension(amr_slots(amr_cur)%idwbuff(1)%beg:,amr_slots(amr_cur)%idwbuff(2)%beg:, &
             & amr_slots(amr_cur)%idwbuff(3)%beg:,1:,1:), intent(in) :: pb_s, mv_s
        real(wp), dimension(amr_slots(amr_cur)%idwbuff(1)%beg:,amr_slots(amr_cur)%idwbuff(2)%beg:, &
             & amr_slots(amr_cur)%idwbuff(3)%beg:,1:,1:), intent(in) :: rpb, rmv
        real(wp), intent(in) :: c1, c2, c3, c4, dtl
        integer              :: fi, fj, fk, q, ib_, fm, fn, fp

        fm = amr_slots(amr_cur)%m; fn = amr_slots(amr_cur)%n; fp = amr_slots(amr_cur)%p
        $:GPU_PARALLEL_LOOP(collapse=5)
        do ib_ = 1, nb
            do q = 1, nnode
                do fk = 0, fp
                    do fj = 0, fn
                        do fi = 0, fm
                            pb_u(fi, fj, fk, q, ib_) = (c1*pb_u(fi, fj, fk, q, ib_) + c2*pb_s(fi, fj, fk, q, &
                                 & ib_) + c3*dtl*rpb(fi, fj, fk, q, ib_))/c4
                            mv_u(fi, fj, fk, q, ib_) = (c1*mv_u(fi, fj, fk, q, ib_) + c2*mv_s(fi, fj, fk, q, &
                                 & ib_) + c3*dtl*rmv(fi, fj, fk, q, ib_))/c4
                        end do
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_fine_rk_update_pbmv

    !> Non-polytropic QBMM: volume-weighted restriction of the block's pb/mv onto the coarse side-state under the block (device
    !! kernel; same equal-weight child average as the q_cons restrict).
    impure subroutine s_restrict_pbmv(pb_c, mv_c, pb_fin, mv_fin)

        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_c, mv_c

        real(stp), dimension(amr_slots(amr_cur)%idwbuff(1)%beg:,amr_slots(amr_cur)%idwbuff(2)%beg:, &
             & amr_slots(amr_cur)%idwbuff(3)%beg:,1:,1:), intent(in) :: pb_fin, mv_fin
        integer  :: ci, cj, ck, q, ib_, fi0, fj0, fk0, ddi, ddj, ddk, nchild, ox, oy, oz, rr
        integer  :: c1lo, c1hi, c2lo, c2hi, c3lo, c3hi, dj_hi, dk_hi
        real(wp) :: accp, accm

        ox = start_idx(1); oy = 0; oz = 0
        if (n_glb > 0) oy = start_idx(2)
        if (p_glb > 0) oz = start_idx(3)
        rr = amr_slots(amr_cur)%amr_ref_ratio
        nchild = rr
        if (n_glb > 0) nchild = nchild*rr
        if (p_glb > 0) nchild = nchild*rr
        c1lo = amr_isect_lo(1); c1hi = amr_isect_hi(1)
        c2lo = amr_isect_lo(2); c2hi = merge(amr_isect_hi(2), amr_isect_lo(2), n_glb > 0)
        c3lo = amr_isect_lo(3); c3hi = merge(amr_isect_hi(3), amr_isect_lo(3), p_glb > 0)
        dj_hi = merge(rr - 1, 0, n_glb > 0); dk_hi = merge(rr - 1, 0, p_glb > 0)
        $:GPU_PARALLEL_LOOP(collapse=5, private='[fi0, fj0, fk0, ddi, accp, accm, ddj, ddk]')
        do ib_ = 1, nb
            do q = 1, nnode
                do ck = c3lo, c3hi
                    do cj = c2lo, c2hi
                        do ci = c1lo, c1hi
                            fi0 = (ci - c1lo)*rr; fj0 = (cj - c2lo)*rr; fk0 = (ck - c3lo)*rr
                            accp = 0._wp; accm = 0._wp
                            do ddk = 0, dk_hi
                                do ddj = 0, dj_hi
                                    do ddi = 0, rr - 1
                                        accp = accp + real(pb_fin(fi0 + ddi, fj0 + ddj, fk0 + ddk, q, ib_), wp)
                                        accm = accm + real(mv_fin(fi0 + ddi, fj0 + ddj, fk0 + ddk, q, ib_), wp)
                                    end do
                                end do
                            end do
                            pb_c(ci - ox, cj - oy, ck - oz, q, ib_) = accp/real(nchild, wp)
                            mv_c(ci - ox, cj - oy, ck - oz, q, ib_) = accm/real(nchild, wp)
                        end do
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_restrict_pbmv

    !> Non-polytropic QBMM np>=2: DEVICE restriction of pb/mv over the covered box [bl:bh] GLOBAL into pb_c/mv_c (device), fine
    !! origin (ci-rlo)*rr, LOCAL coarse index cell - o. Only the covered cells the owner holds are touched (no whole-array push -
    !! same GPU-only clobber the q_cons device overwrite avoids). Same child-sum as s_restrict_pbmv. TWIN:
    !! s_amr_restrict_pbmv_pack_device runs this same child-sum into a wire buffer - any change to the loop order, arithmetic, or
    !! casts here must be mirrored there byte-identically (local and scattered coarse pb/mv must match bit-for-bit). TWIN(pb/mv<->q)
    !! s_amr_restrict_overwrite_device is the q_cons sibling of this child-sum - keep the stencil lockstep.
    impure subroutine s_amr_restrict_pbmv_box_device(pb_c, mv_c, pb_fin, mv_fin, bl, bh, o1, o2, o3, rlo, rr, dj_hi, dk_hi, nchild)

        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_c, mv_c

        real(stp), dimension(amr_slots(amr_cur)%idwbuff(1)%beg:,amr_slots(amr_cur)%idwbuff(2)%beg:, &
             & amr_slots(amr_cur)%idwbuff(3)%beg:,1:,1:), intent(in) :: pb_fin, mv_fin
        integer, intent(in) :: bl(3), bh(3), o1, o2, o3, rlo(3), rr, dj_hi, dk_hi, nchild
        integer             :: ci, cj, ck, q, ib_, fi0, fj0, fk0, ddi, ddj, ddk, bl1, bl2, bl3, bh1, bh2, bh3, rl1, rl2, rl3
        real(wp)            :: accp, accm

        bl1 = bl(1); bl2 = bl(2); bl3 = bl(3); bh1 = bh(1); bh2 = bh(2); bh3 = bh(3)
        rl1 = rlo(1); rl2 = rlo(2); rl3 = rlo(3)
        $:GPU_PARALLEL_LOOP(collapse=5, private='[fi0, fj0, fk0, ddi, accp, accm, ddj, ddk]')
        do ib_ = 1, nb
            do q = 1, nnode
                do ck = bl3, bh3
                    do cj = bl2, bh2
                        do ci = bl1, bh1
                            fi0 = (ci - rl1)*rr; fj0 = (cj - rl2)*rr; fk0 = (ck - rl3)*rr
                            accp = 0._wp; accm = 0._wp
                            do ddk = 0, dk_hi
                                do ddj = 0, dj_hi
                                    do ddi = 0, rr - 1
                                        accp = accp + real(pb_fin(fi0 + ddi, fj0 + ddj, fk0 + ddk, q, ib_), wp)
                                        accm = accm + real(mv_fin(fi0 + ddi, fj0 + ddj, fk0 + ddk, q, ib_), wp)
                                    end do
                                end do
                            end do
                            pb_c(ci - o1, cj - o2, ck - o3, q, ib_) = real(accp/real(nchild, wp), stp)
                            mv_c(ci - o1, cj - o2, ck - o3, q, ib_) = real(accm/real(nchild, wp), stp)
                        end do
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_restrict_pbmv_box_device

    !> Non-polytropic QBMM np>=2 scatter pack: DEVICE restriction of pb/mv over the covered box [bl:bh] GLOBAL straight into the
    !! contiguous wire buffer buf (host, via copyout; pb block then mv block, ci fastest) - only the slice crosses PCIe, not the
    !! full fine side-state. Same child-sum as s_amr_restrict_pbmv_box_device; the wire carries wp and the receiver casts to stp.
    !! TWIN: s_amr_restrict_pbmv_box_device runs this same child-sum in place - any change to the loop order, arithmetic, or casts
    !! here must be mirrored there byte-identically (local and scattered coarse pb/mv must match bit-for-bit). TWIN(pb/mv<->q)
    !! s_amr_restrict_pack_device is the q_cons sibling of this packed child-sum - keep lockstep.
    impure subroutine s_amr_restrict_pbmv_pack_device(pb_fin, mv_fin, bl, bh, rlo, rr, dj_hi, dk_hi, nchild, buf)

        real(stp), dimension(amr_slots(amr_cur)%idwbuff(1)%beg:,amr_slots(amr_cur)%idwbuff(2)%beg:, &
             & amr_slots(amr_cur)%idwbuff(3)%beg:,1:,1:), intent(in) :: pb_fin, mv_fin
        integer, intent(in) :: bl(3), bh(3), rlo(3), rr, dj_hi, dk_hi, nchild
        real(wp), intent(inout), contiguous :: buf(:)
        integer :: ci, cj, ck, q, ib_, fi0, fj0, fk0, ddi, ddj, ddk, bl1, bl2, bl3, bh1, bh2, bh3, rl1, rl2, rl3
        integer :: n1, n2, n3, half
        real(wp) :: accp, accm

        bl1 = bl(1); bl2 = bl(2); bl3 = bl(3); bh1 = bh(1); bh2 = bh(2); bh3 = bh(3)
        rl1 = rlo(1); rl2 = rlo(2); rl3 = rlo(3)
        n1 = bh1 - bl1 + 1; n2 = bh2 - bl2 + 1; n3 = bh3 - bl3 + 1
        half = n1*n2*n3*nnode*nb
        $:GPU_PARALLEL_LOOP(collapse=5, private='[fi0, fj0, fk0, ddi, ddj, ddk, accp, accm]', copyout='[buf]')
        do ib_ = 1, nb
            do q = 1, nnode
                do ck = bl3, bh3
                    do cj = bl2, bh2
                        do ci = bl1, bh1
                            fi0 = (ci - rl1)*rr; fj0 = (cj - rl2)*rr; fk0 = (ck - rl3)*rr
                            accp = 0._wp; accm = 0._wp
                            do ddk = 0, dk_hi
                                do ddj = 0, dj_hi
                                    do ddi = 0, rr - 1
                                        accp = accp + real(pb_fin(fi0 + ddi, fj0 + ddj, fk0 + ddk, q, ib_), wp)
                                        accm = accm + real(mv_fin(fi0 + ddi, fj0 + ddj, fk0 + ddk, q, ib_), wp)
                                    end do
                                end do
                            end do
                            buf(1 + (ci - bl1) + n1*((cj - bl2) + n2*((ck - bl3) + n3*((q - 1) + nnode*(ib_ - 1))))) &
                                & = accp/real(nchild, wp)
                            buf(half + 1 + (ci - bl1) + n1*((cj - bl2) + n2*((ck - bl3) + n3*((q - 1) + nnode*(ib_ - 1))))) &
                                & = accm/real(nchild, wp)
                        end do
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_restrict_pbmv_pack_device

    !> Device unpack of a received pb/mv covered slice into the coarse fields - exact inverse of s_amr_restrict_pbmv_pack_device's
    !! wire layout (ci fastest, then cj, ck, q, ib_, all of pb followed by all of mv). Unpacking on the DEVICE is required, not a
    !! convenience: a host unpack plus a strided GPU_UPDATE(device=) of the covered box is miscopied as a contiguous run - see the
    !! note at the q_cons unpack in s_restrict_fine_to_coarse, of which this is the pb/mv twin.
    impure subroutine s_amr_restrict_pbmv_unpack_device(pb_c, mv_c, bl, bh, o1, o2, o3, buf)

        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_c, mv_c
        integer, intent(in) :: bl(3), bh(3), o1, o2, o3
        real(wp), intent(inout), contiguous :: buf(:)
        integer :: ci, cj, ck, q, ib_, bl1, bl2, bl3, bh1, bh2, bh3, n1, n2, n3, half

        bl1 = bl(1); bl2 = bl(2); bl3 = bl(3); bh1 = bh(1); bh2 = bh(2); bh3 = bh(3)
        n1 = bh1 - bl1 + 1; n2 = bh2 - bl2 + 1; n3 = bh3 - bl3 + 1
        half = n1*n2*n3*nnode*nb
        $:GPU_PARALLEL_LOOP(collapse=5, copyin='[buf]')
        do ib_ = 1, nb
            do q = 1, nnode
                do ck = bl3, bh3
                    do cj = bl2, bh2
                        do ci = bl1, bh1
                            pb_c(ci - o1, cj - o2, ck - o3, q, &
                                 & ib_) = real(buf(1 + (ci - bl1) + n1*((cj - bl2) + n2*((ck - bl3) + n3*((q - 1) + nnode*(ib_ &
                                 & - 1))))), stp)
                            mv_c(ci - o1, cj - o2, ck - o3, q, &
                                 & ib_) = real(buf(half + 1 + (ci - bl1) + n1*((cj - bl2) + n2*((ck - bl3) + n3*((q - 1) &
                                 & + nnode*(ib_ - 1))))), stp)
                        end do
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_restrict_pbmv_unpack_device

    !> Non-polytropic QBMM: distributed fine->coarse fold-back of the block's pb/mv onto the coarse side-state pb_ts/mv_ts,
    !! mirroring the q_cons scatter in s_restrict_fine_to_coarse. The owner restricts the covered cells it holds (device, owned box)
    !! and SENDS each other coarse-owner its covered pb/mv slice (device-packed averages, pb block then mv block); that owner
    !! overwrites its local coarse. Called on ALL ranks at np>=2 (owner/non-owner split inside) so the P2P send/recv pair up; np=1
    !! is handled locally by the direct s_restrict_pbmv call (this routine is reached only when num_procs > 1). TWIN
    !! s_restrict_fine_to_coarse (scatter half, pb/mv<->q): same scatter skeleton (owner sends covered coarse slices, each
    !! coarse-owner overwrites its local cells) - keep lockstep.
    impure subroutine s_amr_scatter_pbmv(pb_fin, mv_fin)

        real(stp), dimension(amr_slots(amr_cur)%idwbuff(1)%beg:,amr_slots(amr_cur)%idwbuff(2)%beg:, &
             & amr_slots(amr_cur)%idwbuff(3)%beg:,1:,1:), intent(in) :: pb_fin, mv_fin
        integer               :: nchild, rr, dj_hi, dk_hi, o1, o2, o3, owner, r, idx, boxsz, maxsz, nsrc, ierr, cellsz
        integer               :: rlo(3), rhi(3), ilo(3), ihi(3), bl(3), bh(3)
        real(wp), allocatable :: sbuf(:,:), rbuf(:)
        integer, allocatable  :: reqs(:), drank(:)

        rr = amr_slots(amr_cur)%amr_ref_ratio
        nchild = rr; if (n_glb > 0) nchild = nchild*rr; if (p_glb > 0) nchild = nchild*rr
        dj_hi = merge(rr - 1, 0, n_glb > 0); dk_hi = merge(rr - 1, 0, p_glb > 0)
        cellsz = 2*nnode*nb
        rlo = 0; rhi = 0
        rlo(1) = amr_region_lo_all(1, amr_cur); rhi(1) = amr_region_hi_all(1, amr_cur)
        if (n_glb > 0) then; rlo(2) = amr_region_lo_all(2, amr_cur); rhi(2) = amr_region_hi_all(2, amr_cur); end if
        if (p_glb > 0) then; rlo(3) = amr_region_lo_all(3, amr_cur); rhi(3) = amr_region_hi_all(3, amr_cur); end if
        owner = amr_block_owner(amr_cur)
        o1 = start_idx(1); o2 = 0; o3 = 0
        if (n_glb > 0) o2 = start_idx(2)
        if (p_glb > 0) o3 = start_idx(3)
        maxsz = cellsz*(rhi(1) - rlo(1) + 1)*(rhi(2) - rlo(2) + 1)*(rhi(3) - rlo(3) + 1)

        ! block set changed: rebuild the cached overlap-rank lists (same lazy trigger as s_amr_fine_fine_halo; local, replicated)
        if (amr_seam_pairs_dirty .or. amr_seam_pairs_nblk /= amr_num_blocks) call s_amr_build_seam_pairs()

        if (proc_rank == owner) then
            ! overwrite the covered cells this rank owns (device, owned box), then send each other coarse-owner its covered slice
            call s_amr_rank_interior(proc_rank, ilo, ihi)
            call s_amr_box_isect(rlo, rhi, ilo, ihi, bl, bh)
            if (bl(1) <= bh(1) .and. bl(2) <= bh(2) .and. bl(3) <= bh(3)) call s_amr_restrict_pbmv_box_device(pb_ts(1)%sf, &
                & mv_ts(1)%sf, pb_fin, mv_fin, bl, bh, o1, o2, o3, rlo, rr, dj_hi, dk_hi, nchild)
            ! cached destination list (every listed rank's interior overlaps the region by construction)
            nsrc = 0
            do idx = 1, amr_ovl_scatter_n(amr_cur)
                if (amr_ovl_scatter(idx, amr_cur) /= owner) nsrc = nsrc + 1
            end do
            if (nsrc > 0) then
                allocate (sbuf(maxsz, nsrc), reqs(nsrc), drank(nsrc))
                nsrc = 0
                do idx = 1, amr_ovl_scatter_n(amr_cur)
                    r = amr_ovl_scatter(idx, amr_cur)
                    if (r == owner) cycle
                    call s_amr_rank_interior(r, ilo, ihi)
                    call s_amr_box_isect(rlo, rhi, ilo, ihi, bl, bh)
                    nsrc = nsrc + 1; drank(nsrc) = r
                    boxsz = cellsz*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                    ! pack this destination's covered pb/mv slice on the DEVICE (restrict averages straight into the wire
                    ! buffer, same child-sum as the device overwrite above) - no full-field fine host pull
                    call s_amr_restrict_pbmv_pack_device(pb_fin, mv_fin, bl, bh, rlo, rr, dj_hi, dk_hi, nchild, sbuf(1:boxsz,nsrc))
#ifdef MFC_MPI
                    call s_xa_rec(XA_F7C_SND, 1, boxsz, amr_cur)
                    call MPI_ISEND(sbuf(1, nsrc), boxsz, mpi_p, r, amr_cur, MPI_COMM_WORLD, reqs(nsrc), ierr)
#endif
                end do
#ifdef MFC_MPI
                call MPI_WAITALL(nsrc, reqs, MPI_STATUSES_IGNORE, ierr)
#endif
                deallocate (sbuf, reqs, drank)
            end if
        else
            ! coarse-owner: if I hold covered cells, receive my pb/mv slice from the owner and overwrite my local coarse
            call s_amr_rank_interior(proc_rank, ilo, ihi)
            call s_amr_box_isect(rlo, rhi, ilo, ihi, bl, bh)
            if (bl(1) <= bh(1) .and. bl(2) <= bh(2) .and. bl(3) <= bh(3)) then
                boxsz = cellsz*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                allocate (rbuf(boxsz))
#ifdef MFC_MPI
                call s_xa_rec(XA_F7C_RCV, 2, boxsz, amr_cur)
                call MPI_RECV(rbuf, boxsz, mpi_p, owner, amr_cur, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
#endif
                ! DEVICE unpack, writing only the covered cells (a whole-array push would clobber device-advanced non-covered
                ! coarse cells with this rank's stale host copy)
                call s_amr_restrict_pbmv_unpack_device(pb_ts(1)%sf, mv_ts(1)%sf, bl, bh, o1, o2, o3, rbuf)
                deallocate (rbuf)
            end if
        end if

    end subroutine s_amr_scatter_pbmv

    !> Swap the global grid state to the fine block. MUST be paired with s_amr_restore_coarse.
    impure subroutine s_amr_swap_to_fine()

        ! Saving on a NESTED swap would overwrite the sw_* bounce buffers with FINE state, and the eventual restore would install
        ! fine extents as the coarse grid - silent corruption of everything after. Hence every save below is depth-guarded; the
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
        ! active-box windows are COARSE cell indices: applying them on the swapped fine grid would window the wrong cells. Blocks
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
        ! Extend the fine grid into the ghost shell (s_build_level_coords only fills the interior 0:m). Ghost cells use the EXACT
        ! parent-cell bisection - the same formula as the interior, with floor division for negative indices. Fine-level
        ! distribution: the owner may not hold the block's coarse coordinate slice locally, so ghost parent boundaries come from the
        ! GLOBAL boundaries amr_g?cb (cl is a GLOBAL coarse index, region_lo + floor(jg/rr)), matching the interior build. Blocks
        ! stay buff_size inside the domain, so every ghost parent is an in-domain coarse cell with exact coords.
        block
            integer               :: jg, cl, pblk2, k, rr, pnf
            real(wp), allocatable :: cxb(:), cyb(:), czb(:), tcc(:), tdx(:)
            rr = amr_slots(amr_cur)%amr_ref_ratio
            ! ghost parent boundaries: a level>=2 block's coarse side is its PARENT's fine grid (indexed in the parent-fine
            ! amr_isect
            ! frame, matching the interior s_build_level_coords), NOT the L0 global boundaries. amr_isect_lo is a parent-fine index,
            ! so indexing amr_g?cb (sized for L0) reads OUT OF BOUNDS -> garbage on host, NaN on the device copy. Source the
            ! parent's
            ! fine coords for level>=2, the global L0 boundaries for level 1.
            if (amr_block_level(amr_cur) >= 2) then
                ! REBUILD the parent's fine boundaries from replicated metadata - do NOT read amr_slots(pblk2)%x_cb. That array is
                ! allocated only on the PARENT's owner, and under per-level distribution this block's owner need not be it; taking
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
        ! sync the swapped extents/bounds/coordinates to the device: RHS kernels read the device copies of these GPU_DECLARE'd
        ! globals (stale coarse bounds = OOB kernels)
        call s_amr_sync_grid_state_to_device()
        ! hypoelastic stress sources use grid-spacing-dependent FD coefficients: recompute them from the (now fine) grid, else every
        ! fine velocity gradient is halved
        if (hypoelasticity) call s_hypoelastic_update_fd_coeffs()
        ! nonuniform coarse grid (stretched, or the axisymmetric axis half-cell): the per-cell WENO coefficients must be rebuilt for
        ! the block's own grid (no-op flag on uniform grids)
        if (amr_weno_coef_recompute) call s_amr_recompute_weno_coefs()

        ! IGR: save the coarse sigma state and seed the fine solve. jac holds THIS stage's converged coarse sigma (the coarse RHS
        ! ran
        ! first), so its parent values are both the best initial guess and the frozen Dirichlet ghost data for the block-local
        ! Jacobi
        ! solve (the per-iteration BC/halo populate is skipped under amr_in_fine_advance). Piecewise-constant parent injection over
        ! the full buffered fine range.
        if (igr) call s_amr_igr_swap_sigma()

    end subroutine s_amr_swap_to_fine

    !> Restore the global grid state saved by s_amr_swap_to_fine.
    impure subroutine s_amr_restore_coarse()

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
        call s_amr_sync_grid_state_to_device()
        if (hypoelasticity) call s_hypoelastic_update_fd_coeffs()
        if (amr_weno_coef_recompute) call s_amr_recompute_weno_coefs()
        if (igr) call s_amr_igr_restore_sigma()

    end subroutine s_amr_restore_coarse

    !> Save the coarse jac/jac_old and seed the (already swapped-in) fine block's sigma state by piecewise-constant parent injection
    !! from the saved coarse sigma: interior = initial guess, ghost shell = frozen Dirichlet coupling data for the block-local
    !! solve.
    impure subroutine s_amr_igr_swap_sigma()

        integer :: j, k, l, ci, cj, ck
        integer :: cb1, ce1, cb2, ce2, cb3, ce3, fb1, fe1, fb2, fe2, fb3, fe3
        integer :: lo1, lo2, lo3, ox, oy, oz

        ! bounds/offsets hoisted to scalars: sw_idwbuff (and friends) are host-only module state - referencing them inside the
        ! kernels makes OpenACC's present lookup fail (OpenMP's implicit map(to) tolerates it, which is why only acc lanes crashed)

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
        ! SAVE the coarse sigma - outermost swap only. A nested swap must not re-save, or sw_jac would take FINE state and both the
        ! seed below and s_amr_igr_restore_sigma would work from it. The seed that follows is NOT guarded: it reads sw_jac, which
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

    !> Recompute the WENO reconstruction coefficient arrays from the CURRENT grid globals (the fine block's after a swap, the coarse
    !! grid's after a restore). s_compute_weno_coefficients reads the live cell-boundary arrays, refreshes uniform_grid, and pushes
    !! its own device updates; the coefficient arrays were sized to the coarse local ranges at init, which the fine ranges never
    !! exceed (the scratch-fit abort guarantees it).
    impure subroutine s_amr_recompute_weno_coefs()

        type(int_bounds_info) :: is1, is2, is3

        is1%beg = -buff_size; is1%end = m + buff_size
        call s_compute_weno_coefficients(1, is1)
        if (n_glb > 0) then
            is2%beg = -buff_size; is2%end = n + buff_size
            call s_compute_weno_coefficients(2, is2)
        end if
        if (p_glb > 0) then
            is3%beg = -buff_size; is3%end = p + buff_size
            call s_compute_weno_coefficients(3, is3)
        end if

    end subroutine s_amr_recompute_weno_coefs

    !> Push the (host-side) global grid state to its device copies after a swap/restore. m/n/p, idwint/idwbuff, and the coordinate
    !! arrays are GPU_DECLARE'd; kernels read the device copies. No-op on CPU.
    impure subroutine s_amr_sync_grid_state_to_device()

        $:GPU_UPDATE(device='[m, n, p, idwint, idwbuff]')
        $:GPU_UPDATE(device='[x_cb, x_cc, dx]')
        if (n_glb > 0) then
            $:GPU_UPDATE(device='[y_cb, y_cc, dy]')
        end if
        if (p_glb > 0) then
            $:GPU_UPDATE(device='[z_cb, z_cc, dz]')
        end if

    end subroutine s_amr_sync_grid_state_to_device

    !> Decompose the current fine block's ghost shell (buffered extent minus interior) into ns disjoint face slabs whose union is
    !! exactly the non-interior cells, so the ghost-fill kernels do O(surface) work instead of masking the full buffered volume. x
    !! slabs span the full transverse extent; y slabs restrict x to the interior; z slabs restrict x and y. Collapsed dims
    !! (n_glb/p_glb == 0) contribute no slabs.
    pure subroutine s_amr_build_ghost_slabs(ns, sb1, se1, sb2, se2, sb3, se3)

        integer, intent(out)               :: ns
        integer, dimension(6), intent(out) :: sb1, se1, sb2, se2, sb3, se3
        integer                            :: fm, fn, fp, b1, e1, b2, e2, b3, e3

        fm = amr_slots(amr_cur)%m; fn = amr_slots(amr_cur)%n; fp = amr_slots(amr_cur)%p
        b1 = amr_slots(amr_cur)%idwbuff(1)%beg; e1 = amr_slots(amr_cur)%idwbuff(1)%end
        b2 = amr_slots(amr_cur)%idwbuff(2)%beg; e2 = amr_slots(amr_cur)%idwbuff(2)%end
        b3 = amr_slots(amr_cur)%idwbuff(3)%beg; e3 = amr_slots(amr_cur)%idwbuff(3)%end
        ns = 2
        sb1(1) = b1; se1(1) = -1; sb1(2) = fm + 1; se1(2) = e1
        sb2(1:2) = b2; se2(1:2) = e2; sb3(1:2) = b3; se3(1:2) = e3
        if (n_glb > 0) then
            ns = 4
            sb2(3) = b2; se2(3) = -1; sb2(4) = fn + 1; se2(4) = e2
            sb1(3:4) = 0; se1(3:4) = fm; sb3(3:4) = b3; se3(3:4) = e3
        end if
        if (p_glb > 0) then
            ns = 6
            sb3(5) = b3; se3(5) = -1; sb3(6) = fp + 1; se3(6) = e3
            sb1(5:6) = 0; se1(5:6) = fm; sb2(5:6) = 0; se2(5:6) = fn
        end if

    end subroutine s_amr_build_ghost_slabs

    !> Fill the fine ghost shell of q_fine by conservative-linear prolongation from q_coarse - the gathered block-local coarse patch
    !! amr_cg (fine-level distribution; the caller gathers the source first). Device kernel: reads the patch and writes the fine
    !! target in device memory. floor/modulo mapping is valid for negative fine indices (ghosts). Interior untouched. Multi-fluid
    !! volume fractions get the same sum-preserving closure as the interior prolongation (second kernel). TWIN
    !! s_amr_fill_fine_ghosts_pbmv (q<->pb/mv): pb/mv sibling; keep the mapping lockstep.
    !!
    !! ONE body, several targets. The prolongation is identical whatever it writes into, and the target differs only in the write
    !! EXPRESSION, so the variants are generated from a single source body with a Fypp accessor lambda (the idiom
    !! m_riemann_solver_hlld already uses for its per-direction stencil variants). Branching on the target inside one region is NOT
    !! an option: a dummy referenced in ANY branch of a target region is still mapped, which is the per-region tax this whole
    !! exercise exists to remove. `_sf` writes a scalar_field vector (q_cons, until it migrates); `_gsta`/`_gstb` write the flat
    !! per-block store at dense local index `loc`.
    #:for SFX, TGT in [('cons', 'amr_cons_st'), ('gsta', 'amr_gst_a'), ('gstb', 'amr_gst_b')]
        #:set QF = lambda ix: TGT + '(fi, fj, fk, ' + ix + ', loc)'
        impure subroutine s_amr_fill_fine_ghosts_${SFX}$(q_coarse, loc)

            type(scalar_field), dimension(sys_size), intent(in) :: q_coarse
            integer, intent(in)                                 :: loc
            integer                                             :: i, fi, fj, fk, ci, cj, ck, ox, oy, oz
            integer                                             :: rr, lo1, lo2, lo3
            integer                                             :: advb, adve, bbeg, bend, bstride
            integer                                             :: s, ns, l1, u1, l2, u2, l3, u3
            integer                                             :: ss, g, r, n1, n2, stot
            integer, dimension(6)                               :: sb1, se1, sb2, se2, sb3, se3, soff, scnt
            logical                                             :: d2, d3, multi, shx, shy, shz, bubEE
            real(wp)                                            :: u0, sx, sy, sz, xix, xiy, xiz, av, asum

            ! q_coarse is the gathered block-local patch amr_cg (fine-level distribution); amr_isect_lo (GLOBAL, == region_lo on the
            ! owner) + f/rr - amr_cpat_off is the patch-local coarse index. Fine indices are LOCAL to this block.

            ox = amr_cpat_off(1); oy = amr_cpat_off(2); oz = amr_cpat_off(3)
            d2 = n_glb > 0; d3 = p_glb > 0
            rr = amr_slots(amr_cur)%amr_ref_ratio
            lo1 = amr_isect_lo(1); lo2 = amr_isect_lo(2); lo3 = amr_isect_lo(3)
            multi = num_fluids > 1 .and. (.not. bubbles_lagrange)  ! EL alphas sum to beta, not 1: no sum-to-one closure
            advb = eqn_idx%adv%beg; adve = eqn_idx%adv%end
            bubEE = bubbles_euler; bbeg = eqn_idx%bub%beg; bend = eqn_idx%bub%end
            bstride = 1; if (bubEE) bstride = (bend - bbeg + 1)/nb
            call s_amr_build_ghost_slabs(ns, sb1, se1, sb2, se2, sb3, se3)
            ! ONE kernel over the concatenation of the ns face slabs instead of one kernel each. The slabs are disjoint and their
            ! union
            ! is exactly the ghost shell (s_amr_build_ghost_slabs), so every ghost cell is still written exactly once and the result
            ! is
            ! byte-identical however the flat index is ordered. NOT the padded-hull form of s_amr_capture_creg_dense_batch: the x
            ! slabs
            ! span the full transverse extent, so a hull over all slabs is the whole buffered volume and masking it would throw away
            ! the
            ! O(surface) decomposition this routine exists to get.
            soff(1) = 0
            do s = 1, ns
                scnt(s) = (se1(s) - sb1(s) + 1)*(se2(s) - sb2(s) + 1)*(se3(s) - sb3(s) + 1)
                if (s < ns) soff(s + 1) = soff(s) + scnt(s)
            end do
            stot = soff(ns) + scnt(ns)
            $:GPU_PARALLEL_LOOP(collapse=2, copyin='[sb1, se1, sb2, se2, sb3, se3, soff, scnt]', private='[s, ss, r, n1, n2, fi, &
                                & fj, fk, ci, cj, ck, xix, xiy, xiz, u0, sx, sy, sz]')
            do i = 1, sys_size
                do g = 0, stot - 1
                    s = 1  ! decode the flat index: ns <= 6, so a scan beats storing a per-cell slab map
                    do ss = 2, ns
                        if (g >= soff(ss)) s = ss
                    end do
                    r = g - soff(s)
                    n1 = se1(s) - sb1(s) + 1; n2 = se2(s) - sb2(s) + 1
                    fi = sb1(s) + mod(r, n1)
                    fj = sb2(s) + mod(r/n1, n2)
                    fk = sb3(s) + r/(n1*n2)
                    ! the slabs cover exactly the ghost shell; multi-fluid, skip the volume fractions (closure kernel below)
                    if (.not. (multi .and. i >= advb .and. i <= adve)) then
                        ck = 0; xiz = 0._wp
                        if (d3) then
                            ck = lo3 + floor(real(fk, wp)/real(rr, wp)) - oz
                            xiz = (real(modulo(fk, rr), wp) - real(rr - 1, wp)*0.5_wp)/real(rr, wp)
                        end if
                        cj = 0; xiy = 0._wp
                        if (d2) then
                            cj = lo2 + floor(real(fj, wp)/real(rr, wp)) - oy
                            xiy = (real(modulo(fj, rr), wp) - real(rr - 1, wp)*0.5_wp)/real(rr, wp)
                        end if
                        ci = lo1 + floor(real(fi, wp)/real(rr, wp)) - ox
                        xix = (real(modulo(fi, rr), wp) - real(rr - 1, wp)*0.5_wp)/real(rr, wp)
                        u0 = real(q_coarse(i)%sf(ci, cj, ck), wp)
                        sx = minmod(real(q_coarse(i)%sf(ci + 1, cj, ck), wp) - u0, u0 - real(q_coarse(i)%sf(ci - 1, cj, ck), wp))
                        sy = 0._wp
                        if (d2) sy = minmod(real(q_coarse(i)%sf(ci, cj + 1, ck), wp) - u0, u0 - real(q_coarse(i)%sf(ci, cj - 1, &
                            & ck), wp))
                        sz = 0._wp
                        if (d3) sz = minmod(real(q_coarse(i)%sf(ci, cj, ck + 1), wp) - u0, u0 - real(q_coarse(i)%sf(ci, cj, &
                            & ck - 1), wp))
                        ! QBMM: inject the bub block piecewise-constant (child = u0) so the ghost inherits the coarse cell's
                        ! realizable 6-moment set (CHyQMOM needs variance c20 > 0; per-component minmod slopes would break
                        ! that joint constraint). Non-QBMM Euler-Euler bubbles instead floor their positive moments (nR /
                        ! npb / nmv); the signed velocity moment nV (offset 1) is skipped.
                        if (qbmm .and. i >= bbeg .and. i <= bend) then
                            sx = 0._wp; sy = 0._wp; sz = 0._wp
                        end if
                        ${QF('i')}$ = u0 + sx*xix + sy*xiy + sz*xiz
                        if (bubEE .and. .not. qbmm .and. i >= bbeg .and. i <= bend) then
                            if (mod(i - bbeg, bstride) /= 1) ${QF('i')}$ = max(real(${QF('i')}$, wp), bub_pos_frac*u0)
                        end if
                    end if
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            ! multi-fluid volume-fraction ghosts: per-cell closure mirroring s_prolong_alphas_closure (shared limiter switch over
            ! all
            ! fluids; interpolate + clamp fluids advb..adve-1; alpha_n = 1 - sum)
            if (multi) then
                ! same flat-index fusion as the prolongation loop above, over the same disjoint slabs
                $:GPU_PARALLEL_LOOP(copyin='[sb1, se1, sb2, se2, sb3, se3, soff, scnt]', private='[s, ss, r, n1, n2, fi, fj, fk, &
                                    & i, ci, cj, ck, xix, xiy, xiz, u0, sx, sy, sz, av, asum, shx, shy, shz]')
                do g = 0, stot - 1
                    s = 1
                    do ss = 2, ns
                        if (g >= soff(ss)) s = ss
                    end do
                    r = g - soff(s)
                    n1 = se1(s) - sb1(s) + 1; n2 = se2(s) - sb2(s) + 1
                    fi = sb1(s) + mod(r, n1)
                    fj = sb2(s) + mod(r/n1, n2)
                    fk = sb3(s) + r/(n1*n2)
                    ck = 0; xiz = 0._wp
                    if (d3) then
                        ck = lo3 + floor(real(fk, wp)/real(rr, wp)) - oz
                        xiz = (real(modulo(fk, rr), wp) - real(rr - 1, wp)*0.5_wp)/real(rr, wp)
                    end if
                    cj = 0; xiy = 0._wp
                    if (d2) then
                        cj = lo2 + floor(real(fj, wp)/real(rr, wp)) - oy
                        xiy = (real(modulo(fj, rr), wp) - real(rr - 1, wp)*0.5_wp)/real(rr, wp)
                    end if
                    ci = lo1 + floor(real(fi, wp)/real(rr, wp)) - ox
                    xix = (real(modulo(fi, rr), wp) - real(rr - 1, wp)*0.5_wp)/real(rr, wp)
                    shx = .true.; shy = d2; shz = d3
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = advb, adve
                        u0 = real(q_coarse(i)%sf(ci, cj, ck), wp)
                        if ((real(q_coarse(i)%sf(ci + 1, cj, ck), wp) - u0)*(u0 - real(q_coarse(i)%sf(ci - 1, cj, ck), &
                            & wp)) <= 0._wp) shx = .false.
                        if (d2) then
                            if ((real(q_coarse(i)%sf(ci, cj + 1, ck), wp) - u0)*(u0 - real(q_coarse(i)%sf(ci, cj - 1, ck), &
                                & wp)) <= 0._wp) shy = .false.
                        end if
                        if (d3) then
                            if ((real(q_coarse(i)%sf(ci, cj, ck + 1), wp) - u0)*(u0 - real(q_coarse(i)%sf(ci, cj, ck - 1), &
                                & wp)) <= 0._wp) shz = .false.
                        end if
                    end do
                    asum = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = advb, adve - 1
                        u0 = real(q_coarse(i)%sf(ci, cj, ck), wp)
                        sx = 0._wp
                        if (shx) sx = minmod(real(q_coarse(i)%sf(ci + 1, cj, ck), wp) - u0, u0 - real(q_coarse(i)%sf(ci - 1, cj, &
                            & ck), wp))
                        sy = 0._wp
                        if (shy) sy = minmod(real(q_coarse(i)%sf(ci, cj + 1, ck), wp) - u0, u0 - real(q_coarse(i)%sf(ci, cj - 1, &
                            & ck), wp))
                        sz = 0._wp
                        if (shz) sz = minmod(real(q_coarse(i)%sf(ci, cj, ck + 1), wp) - u0, u0 - real(q_coarse(i)%sf(ci, cj, &
                            & ck - 1), wp))
                        av = min(max(u0 + sx*xix + sy*xiy + sz*xiz, 0._wp), 1._wp)
                        ${QF('i')}$ = av
                        asum = asum + av
                    end do
                    ${QF('adve')}$ = 1._wp - asum
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

        end subroutine s_amr_fill_fine_ghosts_${SFX}$
    #:endfor

    !> Lerp the fine ghost shell of q_tgt between the coarse t^n and t^{n+1} ghost sources - block loc's slices of the flat store
    !! amr_gst_a/amr_gst_b - at time fraction th (device kernel). Interior untouched. TWIN s_amr_lerp_fine_ghosts_pbmv (q<->pb/mv):
    !! pb/mv sibling of this ghost lerp; keep lockstep.
    impure subroutine s_amr_lerp_fine_ghosts(loc, th)

        integer, intent(in)   :: loc
        real(wp), intent(in)  :: th
        integer               :: i, fi, fj, fk, s, ns, l1, u1, l2, u2, l3, u3
        integer               :: ss, g, r, n1, n2, stot
        integer, dimension(6) :: soff, scnt
        integer, dimension(6) :: sb1, se1, sb2, se2, sb3, se3

        call s_amr_build_ghost_slabs(ns, sb1, se1, sb2, se2, sb3, se3)
        ! flat index over the concatenated DISJOINT slabs - one kernel instead of ns; see s_amr_fill_fine_ghosts
        soff(1) = 0
        do s = 1, ns
            scnt(s) = (se1(s) - sb1(s) + 1)*(se2(s) - sb2(s) + 1)*(se3(s) - sb3(s) + 1)
            if (s < ns) soff(s + 1) = soff(s) + scnt(s)
        end do
        stot = soff(ns) + scnt(ns)
        $:GPU_PARALLEL_LOOP(collapse=2, copyin='[sb1, se1, sb2, se2, sb3, se3, soff, scnt]', &
                            & private='[s, ss, r, n1, n2, fi, fj, fk]')
        do i = 1, sys_size
            do g = 0, stot - 1
                s = 1
                do ss = 2, ns
                    if (g >= soff(ss)) s = ss
                end do
                r = g - soff(s)
                n1 = se1(s) - sb1(s) + 1; n2 = se2(s) - sb2(s) + 1
                fi = sb1(s) + mod(r, n1)
                fj = sb2(s) + mod(r/n1, n2)
                fk = sb3(s) + r/(n1*n2)
                amr_cons_st(fi, fj, fk, i, loc) = (1._wp - th)*real(amr_gst_a(fi, fj, fk, i, loc), wp) + th*real(amr_gst_b(fi, &
                            & fj, fk, i, loc), wp)
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_lerp_fine_ghosts

    !> Non-polytropic QBMM twin of s_amr_lerp_fine_ghosts: lerp the block's pb/mv ghost shell between the coarse t^n and t^{n+1}
    !! sources at the substage time (device kernel; interior untouched). Ghost pb feeds the mixture pressure in the widened
    !! conversion, so it gets the same time treatment as the conservative ghosts. TWIN s_amr_lerp_fine_ghosts (pb/mv<->q): q_cons
    !! sibling; keep lockstep.
    impure subroutine s_amr_lerp_fine_ghosts_pbmv(pb_t, mv_t, pga, mga, pgb, mgb, th)

        real(stp), dimension(amr_slots(amr_cur)%idwbuff(1)%beg:,amr_slots(amr_cur)%idwbuff(2)%beg:, &
             & amr_slots(amr_cur)%idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_t, mv_t
        real(stp), dimension(amr_slots(amr_cur)%idwbuff(1)%beg:,amr_slots(amr_cur)%idwbuff(2)%beg:, &
             & amr_slots(amr_cur)%idwbuff(3)%beg:,1:,1:), intent(in) :: pga, mga, pgb, mgb
        real(wp), intent(in)  :: th
        integer               :: fi, fj, fk, q, ib_, s, ns, l1, u1, l2, u2, l3, u3, ss, g, r, n1, n2, stot
        integer, dimension(6) :: sb1, se1, sb2, se2, sb3, se3, soff, scnt

        call s_amr_build_ghost_slabs(ns, sb1, se1, sb2, se2, sb3, se3)
        ! flat index over the concatenated DISJOINT slabs - one kernel instead of ns; see s_amr_fill_fine_ghosts
        soff(1) = 0
        do s = 1, ns
            scnt(s) = (se1(s) - sb1(s) + 1)*(se2(s) - sb2(s) + 1)*(se3(s) - sb3(s) + 1)
            if (s < ns) soff(s + 1) = soff(s) + scnt(s)
        end do
        stot = soff(ns) + scnt(ns)
        $:GPU_PARALLEL_LOOP(collapse=3, copyin='[sb1, se1, sb2, se2, sb3, se3, soff, scnt]', &
                            & private='[s, ss, r, n1, n2, fi, fj, fk]')
        do ib_ = 1, nb
            do q = 1, nnode
                do g = 0, stot - 1
                    s = 1
                    do ss = 2, ns
                        if (g >= soff(ss)) s = ss
                    end do
                    r = g - soff(s)
                    n1 = se1(s) - sb1(s) + 1; n2 = se2(s) - sb2(s) + 1
                    fi = sb1(s) + mod(r, n1)
                    fj = sb2(s) + mod(r/n1, n2)
                    fk = sb3(s) + r/(n1*n2)
                    pb_t(fi, fj, fk, q, ib_) = (1._wp - th)*real(pga(fi, fj, fk, q, ib_), wp) + th*real(pgb(fi, fj, fk, q, ib_), wp)
                    mv_t(fi, fj, fk, q, ib_) = (1._wp - th)*real(mga(fi, fj, fk, q, ib_), wp) + th*real(mgb(fi, fj, fk, q, ib_), wp)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_lerp_fine_ghosts_pbmv

    !> Device copy q_src -> q_dst over [b1:e1, b2:e2, b3:e3] for all sys_size fields (RK step-entry backup). TWIN s_amr_backup_pbmv
    !! (q<->pb/mv): pb/mv sibling of this step-entry backup; keep lockstep.
    impure subroutine s_amr_copy_fine_fields(loc, b1, e1, b2, e2, b3, e3)

        integer, intent(in) :: loc  !< flat-store slot: source (amr_cons_st) and destination (amr_stor_st) are the same block
        integer, intent(in) :: b1, e1, b2, e2, b3, e3
        integer             :: i, fi, fj, fk

        $:GPU_PARALLEL_LOOP(collapse=4)
        do i = 1, sys_size
            do fk = b3, e3
                do fj = b2, e2
                    do fi = b1, e1
                        amr_stor_st(fi, fj, fk, i, loc) = amr_cons_st(fi, fj, fk, i, loc)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_copy_fine_fields

    !> Device RK stage update over the fine interior: q = (c1*q + c2*q_stor + c3*dt_in*rhs)/c4 (compute in wp, store stp). Mirrors
    !! the coarse non-IGR rk_coef form in s_tvd_rk. TWIN s_amr_fine_rk_update_pbmv + s_tvd_rk (m_time_steppers): same SSP-RK stage
    !! combination; keep all three lockstep.
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

    !> Exchange the coarse conservative ghost layers at internal rank boundaries (physical-boundary ghosts untouched; per direction
    !! beg then end, mirroring s_populate_variables_buffers' disblock). The solver never fills CONS ghosts (only prim), so ranks
    !! whose fine ghost-fill or prolongation stencil leaves their interior need this first. ALL ranks must call together (pairwise
    !! exchange per internal neighbor).
    impure subroutine s_amr_exchange_coarse_cons_halo(q_cons)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons

        if (bc_x%beg >= 0) call s_mpi_sendrecv_variables_buffers(q_cons, 1, -1, sys_size)
        if (bc_x%end >= 0) call s_mpi_sendrecv_variables_buffers(q_cons, 1, 1, sys_size)
        if (n_glb > 0) then
            if (bc_y%beg >= 0) call s_mpi_sendrecv_variables_buffers(q_cons, 2, -1, sys_size)
            if (bc_y%end >= 0) call s_mpi_sendrecv_variables_buffers(q_cons, 2, 1, sys_size)
        end if
        if (p_glb > 0) then
            if (bc_z%beg >= 0) call s_mpi_sendrecv_variables_buffers(q_cons, 3, -1, sys_size)
            if (bc_z%end >= 0) call s_mpi_sendrecv_variables_buffers(q_cons, 3, 1, sys_size)
        end if

    end subroutine s_amr_exchange_coarse_cons_halo

    !> True iff dim d is globally periodic. Uses l0_periodic (periodic_bc allreduced to all ranks in s_l0_tiles_init); periodic_bc
    !! itself is captured from the ORIGINAL bc before MFC folds periodicity into the MPI-cart topology, but only on rank 0, so the
    !! allreduced copy is the one that is consistent on every rank - required since f_amr_seam builds a replicated seam list.
    pure logical function f_l0_dim_periodic(d) result(per)

        integer, intent(in) :: d

        per = l0_periodic(d)

    end function f_l0_dim_periodic

    !> True iff sub-block yb sits immediately above sub-block xb on xb's HIGH face in dim d: yb's low coarse face is xb's high face
    !! + 1, and (tiling produces a regular grid) they share the transverse coarse extents exactly. Each fine-fine seam is exactly
    !! one such ordered (xb, yb) pair (the lower block is xb). For an L0-tile periodic dim there is ALSO a wrap-seam: xb at the
    !! domain HIGH face (region_hi == gcell) and yb at the domain LOW face (region_lo == 0) with matching transverse - the fine-fine
    !! halo then fills xb's high ghost from yb's low interior and vice versa, exactly the periodic connection. Gated on l0_ntile>0
    !! so the AMR fine-block path (blocks never touch the domain edge) is byte-unchanged.
    pure logical function f_amr_seam(xb, yb, d) result(seam)

        integer, intent(in) :: xb, yb, d
        integer             :: t, gc
        logical             :: adj

        adj = amr_region_lo_all(d, yb) == amr_region_hi_all(d, xb) + 1
        if (l0_ntile > 0 .and. f_l0_dim_periodic(d)) then
            gc = merge(m_glb, merge(n_glb, p_glb, d == 2), d == 1)
            adj = adj .or. (amr_region_hi_all(d, xb) == gc .and. amr_region_lo_all(d, yb) == 0)
        end if
        seam = adj
        do t = 1, 3
            if (t /= d) seam = seam .and. amr_region_lo_all(t, xb) == amr_region_lo_all(t, yb) .and. amr_region_hi_all(t, &
                & xb) == amr_region_hi_all(t, yb)
        end do

    end function f_amr_seam

    !> Pack (dir=+1) / unpack (dir=-1) the fine cells of slot's q_cons over [dlo:dhi] in dim d, full transverse, all sys_size, in a
    !! fixed (i, d-index, transverse) order so a packer and unpacker with matching extents align cell-for-cell. GPU: only this
    !! buff_size-deep near-seam slab is moved device<->host (host<-device before a pack, device<-host after an unpack), interior
    !! transverse (0:fm) only - exactly the cells touched below, so the round-trip is byte-identical to a full-field update at a
    !! tiny fraction of the volume (the halo runs per stage, 6x per fine step).
    impure subroutine s_amr_fine_slice(slot, d, dlo, dhi, buf, dir)

        integer, intent(in)                 :: slot, d, dlo, dhi, dir
        real(wp), intent(inout), contiguous :: buf(:)
        integer                             :: i, a, b, c, fm(3), na, nb, nc, loc

        fm(1) = amr_slots(slot)%m; fm(2) = amr_slots(slot)%n; fm(3) = amr_slots(slot)%p
        loc = amr_loc_of(slot)
        nc = dhi - dlo + 1
        ! Pack (dir=1) / unpack (dir=-1) the near-seam slab ON THE DEVICE straight into the contiguous buffer buf, then move ONLY
        ! buf
        ! host<->device. flang miscomputes a STRIDED section (seam dim d < num_dims) of a block's conserved field in a
        ! target-update map clause, corrupting the 2D+ np>1 seam ghosts; the base-grid halo (s_mpi_sendrecv_variables_buffers)
        ! device-packs into a contiguous buffer for the same reason. buf index runs a fastest, then b, then c, then i, so a pack and
        ! an unpack with matching extents align cell-for-cell (na/nb are the transverse fine sizes, nc the slab depth).
        #:for D, TA, TB in [(1, 2, 3), (2, 1, 3), (3, 1, 2)]
            #:set IDX = {1: 'c, a, b', 2: 'a, c, b', 3: 'a, b, c'}[D]
            if (d == ${D}$) then
                na = fm(${TA}$) + 1; nb = fm(${TB}$) + 1  ! scalars; kernel loop bounds MUST use na-1/nb-1, not fm(..), so no host
                !                     array is referenced in the device region (nvfortran/Cray demand it PRESENT)
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

    !> Same-rank seam exchange for ONE pair, both directions in ONE device kernel: xb's high near-seam interior -> yb's low seam
    !! ghost, and yb's low interior -> xb's high ghost. Replaces the four s_amr_fine_slice calls the local branch used to make -
    !! that path routed a purely LOCAL copy through the (since-removed) shared seam buffers, costing four kernels and four blocking
    !! device<->host round trips per pair per stage. Byte-identical to it: the same cells in the same order, and its wp buffer only
    !! widened and re-narrowed the stp values. Fusing the two directions is safe because all four slabs are disjoint - a block's
    !! seam GHOST lies outside its own interior, and the two reads are from different slots than the two writes - so neither
    !! direction can observe the other's store. Both blocks are now addressed by their flat-store slot, so batching across PAIRS is
    !! no longer structurally blocked (a runtime slot index into the module store is a plain subscript, not the null deref the
    !! per-slot scalar_field layout forced). Index order matches s_amr_fine_slice. BATCHED over pairs: one kernel for every
    !! same-rank seam pair on this rank, instead of one per pair (measured 546 of the 7134 AMR-local launches per run). Both blocks
    !! of a pair are addressed by their flat-store slot, so a runtime pair index is now a plain subscript - the per-slot
    !! scalar_field layout is what used to make this a null deref. The per-pair seam dim varies, so the index decode is a runtime
    !! select rather than the Fypp per-dim unroll the single-pair version used. Threads are launched over the MAX pair extent and
    !! masked, rather than prefix-summed: pair sizes are equal under uniform tiling, so the waste is ~0 and it avoids an O(npair)
    !! per-thread offset search.
    impure subroutine s_amr_fine_seam_exchange(npair, plx, ply, pd, pxhi, pfm, ndep)

        integer, intent(in) :: npair
        integer, intent(in) :: plx(:), ply(:), pd(:), pxhi(:), pfm(:,:)  !< per-pair: slots, seam dim, high extent, fine extents
        integer, intent(in) :: ndep
        integer             :: i, a, b, t, na, nb, pr, g, gmax, lx, ly, d, xhi, cnt
        integer             :: ix1, ix2, ix3, iy1, iy2, iy3

        ! max thread extent over the pairs on this rank (transverse product x seam depth)

        gmax = 0
        do pr = 1, npair
            select case (pd(pr))
            case (1); na = pfm(2, pr) + 1; nb = pfm(3, pr) + 1
            case (2); na = pfm(1, pr) + 1; nb = pfm(3, pr) + 1
            case default; na = pfm(1, pr) + 1; nb = pfm(2, pr) + 1
            end select
            gmax = max(gmax, na*nb*ndep)
        end do

        $:GPU_PARALLEL_LOOP(collapse=3, copyin='[plx, ply, pd, pxhi, pfm]', private='[lx, ly, d, xhi, na, nb, cnt, a, b, t, ix1, &
                            & ix2, ix3, iy1, iy2, iy3]')
        do pr = 1, npair
            do i = 1, sys_size
                do g = 0, gmax - 1
                    d = pd(pr)
                    if (d == 1) then
                        na = pfm(2, pr) + 1; nb = pfm(3, pr) + 1
                    else if (d == 2) then
                        na = pfm(1, pr) + 1; nb = pfm(3, pr) + 1
                    else
                        na = pfm(1, pr) + 1; nb = pfm(2, pr) + 1
                    end if
                    cnt = na*nb*ndep
                    if (g < cnt) then
                        lx = plx(pr); ly = ply(pr); xhi = pxhi(pr)
                        a = mod(g, na); b = mod(g/na, nb); t = g/(na*nb)
                        ! (a, b) are the transverse indices, t the depth into the seam; place them per the pair's seam dim
                        if (d == 1) then
                            ix1 = xhi - ndep + 1 + t; ix2 = a; ix3 = b
                            iy1 = -ndep + t; iy2 = a; iy3 = b
                        else if (d == 2) then
                            ix1 = a; ix2 = xhi - ndep + 1 + t; ix3 = b
                            iy1 = a; iy2 = -ndep + t; iy3 = b
                        else
                            ix1 = a; ix2 = b; ix3 = xhi - ndep + 1 + t
                            iy1 = a; iy2 = b; iy3 = -ndep + t
                        end if
                        ! xb high interior -> yb low ghost, then yb low interior -> xb high ghost. All four slabs are disjoint (a
                        ! block's seam ghost lies outside its own interior, and the reads are from different slots than the
                        ! writes), so fusing the two directions cannot let one observe the other's store.
                        amr_cons_st(iy1, iy2, iy3, i, ly) = amr_cons_st(ix1, ix2, ix3, i, lx)
                        if (d == 1) then
                            ix1 = xhi + 1 + t; iy1 = t
                        else if (d == 2) then
                            ix2 = xhi + 1 + t; iy2 = t
                        else
                            ix3 = xhi + 1 + t; iy3 = t
                        end if
                        amr_cons_st(ix1, ix2, ix3, i, lx) = amr_cons_st(iy1, iy2, iy3, i, ly)
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_fine_seam_exchange

    !> Seam dimension of the ordered pair (xb, yb): the dim d in which yb is the immediate high-face neighbour of xb at matched
    !! resolution (same level), or 0 if not a same-level fine-fine seam. Face adjacency requires transverse overlap, so a pair is
    !! adjacent in at most one dim; the last-true assignment reproduces the original inline scan in s_amr_fine_fine_halo.
    pure integer function f_amr_seam_dim(xb, yb) result(d)

        integer, intent(in) :: xb, yb

        d = 0
        if (xb == yb) return
        if (amr_block_level(xb) /= amr_block_level(yb)) return
        if (f_amr_seam(xb, yb, 1)) d = 1
        if (n_glb > 0) then; if (f_amr_seam(xb, yb, 2)) d = 2; end if
        if (p_glb > 0) then; if (f_amr_seam(xb, yb, 3)) d = 3; end if

    end function f_amr_seam_dim

    !> Rebuild the cached same-level seam-pair list (amr_seam_pairs): one O(nblocks^2) scan per regrid/restart in place of the same
    !! scan every RK stage (6x per fine step). Same (xb, yb) nested-loop order on all ranks (replicated region metadata) so the
    !! paired MPI_SENDRECVs in s_amr_fine_fine_halo stay matched. Count then fill for an exact-size list (no cap, no overflow). Also
    !! rebuilds the per-block gather/scatter overlap-rank lists (amr_ovl_gather/scatter) by O(overlap) inversion of the computed
    !! decomposition (s_amr_ranks_overlapping), sized to the max overlap - no O(num_procs) scan or table.
    !> Binary-search the Morton-sorted block lo corners (ord/mkey) for the block whose region lo equals clo, verify the full
    !! same-level seam predicate against xb, and record it in (mb, md, nm). Blocks are disjoint, so at most one block carries a
    !! given lo corner at a level; f_amr_seam_dim takes the LAST true dim, so a block already recorded is raised to the higher d
    !! rather than duplicated.
    pure subroutine s_amr_seam_probe(xb, d, clo, mkey, ord, nblk, mb, md, nm)

        integer, intent(in)         :: xb, d, clo(3), nblk
        integer(kind=8), intent(in) :: mkey(:)
        integer, intent(in)         :: ord(:)
        integer, intent(inout)      :: mb(3), md(3), nm
        integer                     :: yb, t, i, lo_s, hi_s, mid_s
        integer(kind=8)             :: ck
        logical                     :: ok

        ck = f_morton(clo(1), clo(2), clo(3))
        lo_s = 1; hi_s = nblk
        do while (lo_s <= hi_s)
            mid_s = (lo_s + hi_s)/2
            if (mkey(ord(mid_s)) < ck) then
                lo_s = mid_s + 1
            else
                hi_s = mid_s - 1
            end if
        end do
        ! lo_s is the first index with key >= ck; scan the equal-key run. f_morton keeps 21 bits/dim, so
        ! above 2**21 cells/dim a run can hold distinct corners - the explicit lo check rejects those.
        do i = lo_s, nblk
            if (mkey(ord(i)) /= ck) exit
            yb = ord(i)
            if (yb == xb) cycle
            if (amr_block_level(yb) /= amr_block_level(xb)) cycle
            ok = .true.
            do t = 1, 3
                if (amr_region_lo_all(t, yb) /= clo(t)) ok = .false.
                if (t /= d) then
                    if (amr_region_hi_all(t, yb) /= amr_region_hi_all(t, xb)) ok = .false.
                end if
            end do
            if (.not. ok) cycle
            do t = 1, nm
                if (mb(t) == yb) then
                    md(t) = max(md(t), d)
                    return
                end if
            end do
            nm = nm + 1
            mb(nm) = yb; md(nm) = d
            return
        end do

    end subroutine s_amr_seam_probe

    impure subroutine s_amr_build_seam_pairs()

        integer                      :: xb, d, np, k, mx, pass, nm, im, jm, tb, td, nb
        integer                      :: plo(3), phi(3), rlo(3), rhi(3)
        integer                      :: mb(3), md(3), clo(3), gc
        integer                      :: width, lo_m, mid_m, hi_m, i_m, j_m, t_m
        integer(kind=8), allocatable :: mkey(:)
        integer, allocatable         :: ord(:), mrg(:)

        if (allocated(amr_seam_pairs)) deallocate (amr_seam_pairs)

        ! Blocks are disjoint, so (level, region lo) names one uniquely - and the seam predicate FIXES the
        ! neighbour's lo corner: transverse lo equal to xb's, and lo(d) = hi(d, xb) + 1. The all-pairs
        ! O(nblocks^2) scan is therefore a lookup. Morton-sort the lo corners once, then binary-search the
        ! single candidate per (block, dim) and verify the predicate on it. Emission order is unchanged (xb
        ! ascending, yb ascending within xb), which the paired MPI_SENDRECVs in s_amr_fine_fine_halo depend
        ! on - a reordered list mismatches sends to receives and deadlocks.
        nb = max(amr_num_blocks, 1)
        allocate (mkey(nb), ord(nb), mrg(nb))
        do k = 1, amr_num_blocks
            mkey(k) = f_morton(amr_region_lo_all(1, k), amr_region_lo_all(2, k), amr_region_lo_all(3, k))
            ord(k) = k
        end do

        ! Bottom-up STABLE merge sort by Morton key (same form as s_amr_sfc_cut): a pure function of the
        ! replicated region metadata, so every rank builds the identical order.
        width = 1
        do while (width < amr_num_blocks)
            lo_m = 1
            do while (lo_m <= amr_num_blocks - width)
                mid_m = lo_m + width - 1
                hi_m = min(lo_m + 2*width - 1, amr_num_blocks)
                i_m = lo_m; j_m = mid_m + 1; t_m = lo_m
                do while (i_m <= mid_m .and. j_m <= hi_m)
                    if (mkey(ord(i_m)) <= mkey(ord(j_m))) then
                        mrg(t_m) = ord(i_m); i_m = i_m + 1
                    else
                        mrg(t_m) = ord(j_m); j_m = j_m + 1
                    end if
                    t_m = t_m + 1
                end do
                do while (i_m <= mid_m); mrg(t_m) = ord(i_m); i_m = i_m + 1; t_m = t_m + 1; end do
                do while (j_m <= hi_m); mrg(t_m) = ord(j_m); j_m = j_m + 1; t_m = t_m + 1; end do
                ord(lo_m:hi_m) = mrg(lo_m:hi_m)
                lo_m = lo_m + 2*width
            end do
            width = 2*width
        end do

        ! pass 1 counts, pass 2 fills: keeps amr_seam_pairs exactly sized, as the nested scan did
        do pass = 1, 2
            np = 0
            do xb = 1, amr_num_blocks
                nm = 0
                do d = 1, 3
                    if (d == 2 .and. n_glb <= 0) cycle
                    if (d == 3 .and. p_glb <= 0) cycle
                    clo = amr_region_lo_all(:,xb)
                    clo(d) = amr_region_hi_all(d, xb) + 1
                    call s_amr_seam_probe(xb, d, clo, mkey, ord, amr_num_blocks, mb, md, nm)
                    ! periodic wrap (l0 tiling only): xb on the domain high face pairs with lo(d) = 0
                    if (l0_ntile > 0) then
                        if (f_l0_dim_periodic(d)) then
                            gc = merge(m_glb, merge(n_glb, p_glb, d == 2), d == 1)
                            if (amr_region_hi_all(d, xb) == gc) then
                                clo(d) = 0
                                call s_amr_seam_probe(xb, d, clo, mkey, ord, amr_num_blocks, mb, md, nm)
                            end if
                        end if
                    end if
                end do
                ! ascending yb within xb, matching the old inner loop's emission order (nm <= 3)
                do im = 2, nm
                    tb = mb(im); td = md(im); jm = im - 1
                    do while (jm >= 1)
                        if (mb(jm) <= tb) exit
                        mb(jm + 1) = mb(jm); md(jm + 1) = md(jm); jm = jm - 1
                    end do
                    mb(jm + 1) = tb; md(jm + 1) = td
                end do
                do im = 1, nm
                    np = np + 1
                    if (pass == 2) then
                        amr_seam_pairs(1, np) = xb; amr_seam_pairs(2, np) = mb(im); amr_seam_pairs(3, np) = md(im)
                    end if
                end do
            end do
            if (pass == 1) then
                amr_num_seam_pairs = np
                allocate (amr_seam_pairs(3, max(np, 1)))
            end if
        end do
        deallocate (mkey, ord, mrg)
        ! per-block P2P overlap-rank lists by O(overlap) inversion (gather: rank coarse range vs the amr_cpat_mar-padded patch box;
        ! scatter: rank interior vs the region box), rank-ascending so iterating a list preserves the replaced 0..num_procs-1 scans'
        ! MPI send/recv order. The clamped interior-frame coord range reproduces both frames (see s_amr_coord_range). Bounded first
        ! dim = max overlap over all blocks (dealloc-realloc each build, like amr_seam_pairs above), NOT num_procs - a block spans
        ! O(1) ranks, so this collapses the old O(num_procs*amr_max_blocks) table. Every consumer runs behind a build_seam_pairs
        ! guard, so the arrays are always sized before they are read.
        mx = 1
        do k = 1, amr_num_blocks
            plo = 0; phi = 0; rlo = 0; rhi = 0
            plo(1) = amr_region_lo_all(1, k) - amr_cpat_mar; phi(1) = amr_region_hi_all(1, k) + amr_cpat_mar
            rlo(1) = amr_region_lo_all(1, k); rhi(1) = amr_region_hi_all(1, k)
            if (n_glb > 0) then
                plo(2) = amr_region_lo_all(2, k) - amr_cpat_mar; phi(2) = amr_region_hi_all(2, k) + amr_cpat_mar
                rlo(2) = amr_region_lo_all(2, k); rhi(2) = amr_region_hi_all(2, k)
            end if
            if (p_glb > 0) then
                plo(3) = amr_region_lo_all(3, k) - amr_cpat_mar; phi(3) = amr_region_hi_all(3, k) + amr_cpat_mar
                rlo(3) = amr_region_lo_all(3, k); rhi(3) = amr_region_hi_all(3, k)
            end if
            mx = max(mx, f_amr_overlap_count(plo, phi), f_amr_overlap_count(rlo, rhi))
        end do
        if (allocated(amr_ovl_gather)) deallocate (amr_ovl_gather)
        if (allocated(amr_ovl_scatter)) deallocate (amr_ovl_scatter)
        allocate (amr_ovl_gather(mx, amr_max_blocks), amr_ovl_scatter(mx, amr_max_blocks))
        do k = 1, amr_num_blocks
            plo = 0; phi = 0; rlo = 0; rhi = 0
            plo(1) = amr_region_lo_all(1, k) - amr_cpat_mar; phi(1) = amr_region_hi_all(1, k) + amr_cpat_mar
            rlo(1) = amr_region_lo_all(1, k); rhi(1) = amr_region_hi_all(1, k)
            if (n_glb > 0) then
                plo(2) = amr_region_lo_all(2, k) - amr_cpat_mar; phi(2) = amr_region_hi_all(2, k) + amr_cpat_mar
                rlo(2) = amr_region_lo_all(2, k); rhi(2) = amr_region_hi_all(2, k)
            end if
            if (p_glb > 0) then
                plo(3) = amr_region_lo_all(3, k) - amr_cpat_mar; phi(3) = amr_region_hi_all(3, k) + amr_cpat_mar
                rlo(3) = amr_region_lo_all(3, k); rhi(3) = amr_region_hi_all(3, k)
            end if
            call s_amr_ranks_overlapping(plo, phi, amr_ovl_gather(:,k), amr_ovl_gather_n(k))
            call s_amr_ranks_overlapping(rlo, rhi, amr_ovl_scatter(:,k), amr_ovl_scatter_n(k))
        end do
        amr_seam_pairs_nblk = amr_num_blocks
        amr_seam_pairs_dirty = .false.

    end subroutine s_amr_build_seam_pairs

    !> Block-to-block fine-fine halo (max_grid_size tiling): overwrite each sub-block's seam ghost cells (faces shared with an
    !! ADJACENT sub-block) with the neighbour's stage-entry fine interior, so the shared fine flux matches on both sides
    !! (coarse-prolonged seam ghosts would be non-conservative). For each seam pair (xb below, yb above, dim d) the two owners
    !! exchange the buff_size-deep near-seam interior (MPI_Sendrecv, or a local copy when one rank owns both). Buffer is wp, cast to
    !! stp on unpack (identity for stp fields). No-op with a single block / no adjacent pairs (any untiled case, any np).
    impure subroutine s_amr_fine_fine_halo(lev_only)

        !> level to exchange, or 0 for ALL levels. The subcycled level-2 child advance needs to reconcile ONLY its own level's
        !! seams: it runs inside one of the parent's substeps, when the level-1 blocks are mid-substep and must not be touched.
        integer, intent(in) :: lev_only
        integer             :: xb, yb, d, rX, rY, cnt, xm(3), tsz, ierr, fmul, idx
        integer             :: ip, boff, tq, nreq, qbase, r, sblk, sdlo, sdhi, ublk, udlo, udhi, eblk, edlo, edhi
        !> same-rank pairs are collected here and exchanged in ONE kernel; cross-rank pairs stay serial (each is an MPI_SENDRECV)
        integer              :: nsame
        integer, allocatable :: plx(:), ply(:), pd(:), pxhi(:), pfm(:,:)

        if (.not. amr .and. l0_ntile == 0) return
        if (amr_num_blocks < 2) return

        ! iterate the cached same-level seam list (rebuilt only when the topology changes) instead of the old O(nblocks^2)
        ! f_amr_seam
        ! rescan every stage; the list preserves the original (xb, yb) order so paired MPI_SENDRECVs still match.
        if (amr_seam_pairs_dirty .or. amr_seam_pairs_nblk /= amr_num_blocks) call s_amr_build_seam_pairs()
        ! device<->host of the fine state is done per-seam inside s_amr_fine_slice, moving only the buff_size-deep near-seam slab
        ! each
        ! pack/unpack touches (not the whole block) - a large PCIe saving since this runs per stage (6x per fine step)
        allocate (plx(amr_num_seam_pairs), ply(amr_num_seam_pairs), pd(amr_num_seam_pairs), pxhi(amr_num_seam_pairs), pfm(3, &
                  & amr_num_seam_pairs))
        ! I5-F6 WAVE (plan-based exchange, amr_plan_based_exchange.md): the cross-rank pairs - previously one blocking
        ! MPI_SENDRECV each, in pair-list order, through the two shared seam buffers - become ONE aggregated message per
        ! (peer, direction): all recvs posted, all packs into per-transfer pool slices, all sends, one WAITALL, then all
        ! unpacks. Every pair contributes one send and one recv transfer on each of its two owners; both ranks walk the
        ! SAME replicated pair list ascending with per-peer running offsets, so the wire layout agrees with no metadata
        ! exchange (the property the paired SENDRECVs relied on, made explicit). Wire is wp with the same stp cast on
        ! unpack; under MFC_DEBUG each slab carries the identity header [site, sending slot, (d, dlo, dhi), (cnt, 0, 0)].
        ! Reuses the fill waves' scratch/pools (the waves never overlap in time). Same-rank pairs keep the batched kernel.
        if (.not. allocated(amr_fw_map)) then
            allocate (amr_fw_map(0:num_procs - 1), amr_fw_nx(0:num_procs - 1), amr_fw_pq(0:num_procs - 1), &
                      & amr_fw_pp(0:num_procs - 1))
            amr_fw_map = 0; amr_fw_nx = 0; amr_fw_pq = 0; amr_fw_pp = 0
        end if
        tq = amr_tag_base(6) + int(mod(amr_mesh_epoch, 100_8))
        nsame = 0
        amr_fw_snx = 0; amr_fw_snp = 0
        do idx = 1, amr_num_seam_pairs
            xb = amr_seam_pairs(1, idx); yb = amr_seam_pairs(2, idx); d = amr_seam_pairs(3, idx)
            if (lev_only > 0 .and. amr_block_level(xb) /= lev_only) cycle  ! pairs are same-level, so xb's level is the pair's
            rX = amr_block_owner(xb); rY = amr_block_owner(yb)
            if (proc_rank /= rX .and. proc_rank /= rY) cycle
            ! fine extents from the REPLICATED region metadata (not amr_slots%m/n/p: at np>1 this rank may own only one of the pair,
            ! and the transverse size (used for the buffer count) must be valid for both). A level-L block's region is in L0-coarse
            ! cells but its own grid is 2**L finer (each level halves dx), so fine = 2**L*(coarse extent)-1; xb, yb share the level
            ! (same-level seam). 2**1 keeps L1 byte-identical; L2 tiles need 2**2 (an L1-frame 2* mislocates the seam slice to half
            ! the block, filling the seam ghost from the wrong cells - the source of the L2-L2 leak).
            fmul = amr_ref_ratio**amr_block_level(xb)
            xm(1) = fmul*(amr_region_hi_all(1, xb) - amr_region_lo_all(1, xb) + 1) - 1
            xm(2) = merge(fmul*(amr_region_hi_all(2, xb) - amr_region_lo_all(2, xb) + 1) - 1, 0, n_glb > 0)
            xm(3) = merge(fmul*(amr_region_hi_all(3, xb) - amr_region_lo_all(3, xb) + 1) - 1, 0, p_glb > 0)
            ! transverse fine size (dims /= d); xb and yb share it (exact-match seam)
            tsz = 1
            if (d /= 1) tsz = tsz*(xm(1) + 1)
            if (d /= 2 .and. n_glb > 0) tsz = tsz*(xm(2) + 1)
            if (d /= 3 .and. p_glb > 0) tsz = tsz*(xm(3) + 1)
            cnt = sys_size*buff_size*tsz
            if (rX == rY) then  ! same rank owns both: defer to the ONE batched kernel below (no host buffer, no per-pair launch)
                nsame = nsame + 1
                plx(nsame) = amr_loc_of(xb); ply(nsame) = amr_loc_of(yb)
                pd(nsame) = d; pxhi(nsame) = xm(d); pfm(:,nsame) = xm
                cycle
            end if
            ! cross-rank: append this side's SEND transfer (the matching recv is appended in the second pair walk below)
            if (proc_rank == rX) then
                r = rY; sblk = xb; sdlo = xm(d) - buff_size + 1; sdhi = xm(d)
            else
                r = rX; sblk = yb; sdlo = 0; sdhi = buff_size - 1
            end if
            if (amr_fw_map(r) == 0) then
                amr_fw_snp = amr_fw_snp + 1
                call s_amr_fw_szi(amr_fw_sprank, amr_fw_snp); call s_amr_fw_szi(amr_fw_sqsz, amr_fw_snp)
                call s_amr_fw_szi(amr_fw_snxp, amr_fw_snp); call s_amr_fw_szi(amr_fw_sqbase, amr_fw_snp)
                amr_fw_map(r) = amr_fw_snp
                amr_fw_sprank(amr_fw_snp) = r
            end if
            amr_fw_snx = amr_fw_snx + 1
            call s_amr_fw_szi(amr_fw_sblk, amr_fw_snx); call s_amr_fw_szi3(amr_fw_sbl, amr_fw_snx)
            call s_amr_fw_szi(amr_fw_spi, amr_fw_snx); call s_amr_fw_szi(amr_fw_sqo, amr_fw_snx)
            call s_amr_fw_szi(amr_fw_spo, amr_fw_snx)
            amr_fw_sblk(amr_fw_snx) = sblk
            amr_fw_sbl(1, amr_fw_snx) = d; amr_fw_sbl(2, amr_fw_snx) = sdlo; amr_fw_sbl(3, amr_fw_snx) = sdhi
            amr_fw_spo(amr_fw_snx) = cnt
            amr_fw_spi(amr_fw_snx) = amr_fw_map(r)
            amr_fw_sqo(amr_fw_snx) = amr_fw_pq(r) + amr_fw_nx(r)*XA_NH
            amr_fw_pq(r) = amr_fw_pq(r) + cnt
            amr_fw_nx(r) = amr_fw_nx(r) + 1
        end do
        qbase = 0
        do ip = 1, amr_fw_snp
            r = amr_fw_sprank(ip)
            amr_fw_snxp(ip) = amr_fw_nx(r)
            amr_fw_sqsz(ip) = amr_fw_pq(r) + amr_fw_nx(r)*XA_NH
            amr_fw_sqbase(ip) = qbase; qbase = qbase + amr_fw_sqsz(ip)
            amr_fw_map(r) = 0; amr_fw_nx(r) = 0; amr_fw_pq(r) = 0
        end do
        call s_amr_fw_szr(amr_fw_sq, qbase)
        ! RECV transfers, second walk over the same pairs: unpack destination is MY block's ghost slab; the expected header
        ! is the PEER's pack (its block + its interior slab bounds), derived from the same replicated metadata
        amr_fw_rnx = 0; amr_fw_rnp = 0
        do idx = 1, amr_num_seam_pairs
            xb = amr_seam_pairs(1, idx); yb = amr_seam_pairs(2, idx); d = amr_seam_pairs(3, idx)
            if (lev_only > 0 .and. amr_block_level(xb) /= lev_only) cycle
            rX = amr_block_owner(xb); rY = amr_block_owner(yb)
            if (rX == rY) cycle
            if (proc_rank /= rX .and. proc_rank /= rY) cycle
            fmul = amr_ref_ratio**amr_block_level(xb)
            xm(1) = fmul*(amr_region_hi_all(1, xb) - amr_region_lo_all(1, xb) + 1) - 1
            xm(2) = merge(fmul*(amr_region_hi_all(2, xb) - amr_region_lo_all(2, xb) + 1) - 1, 0, n_glb > 0)
            xm(3) = merge(fmul*(amr_region_hi_all(3, xb) - amr_region_lo_all(3, xb) + 1) - 1, 0, p_glb > 0)
            tsz = 1
            if (d /= 1) tsz = tsz*(xm(1) + 1)
            if (d /= 2 .and. n_glb > 0) tsz = tsz*(xm(2) + 1)
            if (d /= 3 .and. p_glb > 0) tsz = tsz*(xm(3) + 1)
            cnt = sys_size*buff_size*tsz
            if (proc_rank == rX) then
                ! yb's low interior arrives -> xb's high ghost
                r = rY; ublk = xb; udlo = xm(d) + 1; udhi = xm(d) + buff_size
                eblk = yb; edlo = 0; edhi = buff_size - 1
            else
                ! xb's high interior arrives -> yb's low ghost
                r = rX; ublk = yb; udlo = -buff_size; udhi = -1
                eblk = xb; edlo = xm(d) - buff_size + 1; edhi = xm(d)
            end if
            if (amr_fw_map(r) == 0) then
                amr_fw_rnp = amr_fw_rnp + 1
                call s_amr_fw_szi(amr_fw_rprank, amr_fw_rnp); call s_amr_fw_szi(amr_fw_rqsz, amr_fw_rnp)
                call s_amr_fw_szi(amr_fw_rnxp, amr_fw_rnp); call s_amr_fw_szi(amr_fw_rqbase, amr_fw_rnp)
                amr_fw_map(r) = amr_fw_rnp
                amr_fw_rprank(amr_fw_rnp) = r
            end if
            amr_fw_rnx = amr_fw_rnx + 1
            call s_amr_fw_szi(amr_fw_rblk, amr_fw_rnx); call s_amr_fw_szi3(amr_fw_rbl, amr_fw_rnx)
            call s_amr_fw_szi3(amr_fw_rbh, amr_fw_rnx); call s_amr_fw_szi(amr_fw_rpi, amr_fw_rnx)
            call s_amr_fw_szi(amr_fw_rqo, amr_fw_rnx); call s_amr_fw_szi(amr_fw_rpo, amr_fw_rnx)
            amr_fw_rblk(amr_fw_rnx) = ublk
            amr_fw_rbl(1, amr_fw_rnx) = d; amr_fw_rbl(2, amr_fw_rnx) = udlo; amr_fw_rbl(3, amr_fw_rnx) = udhi
            amr_fw_rbh(1, amr_fw_rnx) = eblk; amr_fw_rbh(2, amr_fw_rnx) = edlo; amr_fw_rbh(3, amr_fw_rnx) = edhi
            amr_fw_rpo(amr_fw_rnx) = cnt
            amr_fw_rpi(amr_fw_rnx) = amr_fw_map(r)
            amr_fw_rqo(amr_fw_rnx) = amr_fw_pq(r) + amr_fw_nx(r)*XA_NH
            amr_fw_pq(r) = amr_fw_pq(r) + cnt
            amr_fw_nx(r) = amr_fw_nx(r) + 1
        end do
        qbase = 0
        do ip = 1, amr_fw_rnp
            r = amr_fw_rprank(ip)
            amr_fw_rnxp(ip) = amr_fw_nx(r)
            amr_fw_rqsz(ip) = amr_fw_pq(r) + amr_fw_nx(r)*XA_NH
            amr_fw_rqbase(ip) = qbase; qbase = qbase + amr_fw_rqsz(ip)
            amr_fw_map(r) = 0; amr_fw_nx(r) = 0; amr_fw_pq(r) = 0
        end do
        call s_amr_fw_szr(amr_fw_rq, qbase)
        nreq = amr_fw_snp + amr_fw_rnp
        call s_amr_fw_szi(amr_fw_req, nreq); call s_amr_fw_szi(amr_fw_reqw, nreq)

        nreq = 0
#ifdef MFC_MPI
        do ip = 1, amr_fw_rnp
            call s_xa_rec(XA_F6W_RCV, 2, amr_fw_rqsz(ip) - amr_fw_rnxp(ip)*XA_NH, tq)
            nreq = nreq + 1; amr_fw_reqw(nreq) = amr_fw_rqsz(ip)
            call MPI_IRECV(amr_fw_rq(amr_fw_rqbase(ip) + 1), amr_fw_rqsz(ip), mpi_p, amr_fw_rprank(ip), tq, MPI_COMM_WORLD, &
                           & amr_fw_req(nreq), ierr)
        end do
#endif
        do idx = 1, amr_fw_snx
            cnt = amr_fw_spo(idx)
            boff = amr_fw_sqbase(amr_fw_spi(idx)) + amr_fw_sqo(idx)
            call s_amr_fine_slice(amr_fw_sblk(idx), amr_fw_sbl(1, idx), amr_fw_sbl(2, idx), amr_fw_sbl(3, idx), &
                                  & amr_fw_sq(boff + XA_NH + 1:boff + XA_NH + cnt), 1)
            if (XA_NH > 0) call s_xa_hdr_pack(amr_fw_sq(boff + 1:boff + XA_NH), XA_F6W_SND, amr_fw_sblk(idx), amr_fw_sbl(:,idx), &
                & [cnt, 0, 0])
        end do
#ifdef MFC_MPI
        do ip = 1, amr_fw_snp
            call s_xa_rec(XA_F6W_SND, 1, amr_fw_sqsz(ip) - amr_fw_snxp(ip)*XA_NH, tq)
            nreq = nreq + 1; amr_fw_reqw(nreq) = -1
            call MPI_ISEND(amr_fw_sq(amr_fw_sqbase(ip) + 1), amr_fw_sqsz(ip), mpi_p, amr_fw_sprank(ip), tq, MPI_COMM_WORLD, &
                           & amr_fw_req(nreq), ierr)
        end do
        if (nreq > 0) then
#ifdef MFC_DEBUG
            block
                integer :: st(MPI_STATUS_SIZE, nreq), gotw, q
                call MPI_WAITALL(nreq, amr_fw_req, st, ierr)
                do q = 1, nreq
                    if (amr_fw_reqw(q) < 0) cycle
                    call MPI_GET_COUNT(st(:,q), mpi_p, gotw, ierr)
                    @:ASSERT(gotw == amr_fw_reqw(q), "seam wave: received message length differs from the plan")
                end do
            end block
#else
            call MPI_WAITALL(nreq, amr_fw_req, MPI_STATUSES_IGNORE, ierr)
#endif
        end if
#endif
        do idx = 1, amr_fw_rnx
            cnt = amr_fw_rpo(idx)
            boff = amr_fw_rqbase(amr_fw_rpi(idx)) + amr_fw_rqo(idx)
            if (XA_NH > 0) call s_xa_hdr_check(amr_fw_rq(boff + 1:boff + XA_NH), XA_F6W_SND, amr_fw_rbh(1, idx), [amr_fw_rbl(1, &
                & idx), amr_fw_rbh(2, idx), amr_fw_rbh(3, idx)], [cnt, 0, 0])
            call s_amr_fine_slice(amr_fw_rblk(idx), amr_fw_rbl(1, idx), amr_fw_rbl(2, idx), amr_fw_rbl(3, idx), &
                                  & amr_fw_rq(boff + XA_NH + 1:boff + XA_NH + cnt), -1)
        end do
        ! Every same-rank pair in ONE launch. Two things make this safe. Fusing ACROSS pairs: the four slabs of a pair are disjoint
        ! and no pair writes another's source. DEFERRING past the MPI pairs above (a reordering the per-pair loop did not do):
        ! every seam operation reads only INTERIOR cells and writes only GHOST cells - the packs read [xhi-buff+1:xhi] / [0:buff-1]
        ! and the unpacks write [xhi+1:xhi+buff] / [-buff:-1] - so no seam operation can observe another's write, in either path.
        if (nsame > 0) call s_amr_fine_seam_exchange(nsame, plx(1:nsame), ply(1:nsame), pd(1:nsame), pxhi(1:nsame), pfm(:, &
            & 1:nsame), buff_size)
        deallocate (plx, ply, pd, pxhi, pfm)
        call s_amr_select_slot(1)

    end subroutine s_amr_fine_fine_halo

    !> Ring clip: decompose the hollow shell of patch [plo:phi] minus its OPEN core [clo:chi] (same integer frame; collapsed dims
    !! pass plo=phi=clo=chi=0 so they never cut the core) into at most 6 DISJOINT slabs in a FIXED order: x-low/x-high spanning the
    !! full transverse extent, then y-low/high restricted to the core's x-interval, then z-low/high restricted in x and y. Width<=2
    !! cores are empty (shell = whole patch, legal); the width-1 double-cover is resolved by clamping the high slab past the low
    !! one. Both sides of every clipped exchange derive the same list from replicated metadata, so the wire layout needs no
    !! handshake (docs/documentation/amr_stepfill_ring_clip.md).
    impure subroutine s_amr_shell_slabs(plo, phi, clo, chi, ns, sb1, se1, sb2, se2, sb3, se3, cells)

        integer, intent(in)  :: plo(3), phi(3), clo(3), chi(3)
        integer, intent(out) :: ns, sb1(6), se1(6), sb2(6), se2(6), sb3(6), se3(6), cells
        integer              :: cb(3, 6), ce(3, 6), s, ss
        integer(8)           :: words, patchw, corew

        cb(:,1) = plo; ce(:,1) = [clo(1) - 1, phi(2), phi(3)]
        cb(:,2) = [max(chi(1) + 1, clo(1)), plo(2), plo(3)]; ce(:,2) = phi
        cb(:,3) = [clo(1), plo(2), plo(3)]; ce(:,3) = [chi(1), clo(2) - 1, phi(3)]
        cb(:,4) = [clo(1), max(chi(2) + 1, clo(2)), plo(3)]; ce(:,4) = [chi(1), phi(2), phi(3)]
        cb(:,5) = [clo(1), clo(2), plo(3)]; ce(:,5) = [chi(1), chi(2), clo(3) - 1]
        cb(:,6) = [clo(1), clo(2), max(chi(3) + 1, clo(3))]; ce(:,6) = [chi(1), chi(2), phi(3)]
        ns = 0; words = 0
        do s = 1, 6
            if (cb(1, s) > ce(1, s) .or. cb(2, s) > ce(2, s) .or. cb(3, s) > ce(3, s)) cycle
            ns = ns + 1
            sb1(ns) = cb(1, s); se1(ns) = ce(1, s)
            sb2(ns) = cb(2, s); se2(ns) = ce(2, s)
            sb3(ns) = cb(3, s); se3(ns) = ce(3, s)
            words = words + int(se1(ns) - sb1(ns) + 1, 8)*int(se2(ns) - sb2(ns) + 1, 8)*int(se3(ns) - sb3(ns) + 1, 8)
        end do
        ! the slabs must tile the shell EXACTLY: pairwise disjoint, cells summing to patch - core
        do s = 1, ns - 1
            do ss = s + 1, ns
                @:ASSERT(max(sb1(s), sb1(ss)) > min(se1(s), se1(ss)) .or. max(sb2(s), sb2(ss)) > min(se2(s), &
                         & se2(ss)) .or. max(sb3(s), sb3(ss)) > min(se3(s), se3(ss)), "shell slabs: overlap")
            end do
        end do
        patchw = int(phi(1) - plo(1) + 1, 8)*int(phi(2) - plo(2) + 1, 8)*int(phi(3) - plo(3) + 1, 8)
        corew = int(max(chi(1) - clo(1) + 1, 0), 8)*int(max(chi(2) - clo(2) + 1, 0), 8)*int(max(chi(3) - clo(3) + 1, 0), 8)
        @:ASSERT(words == patchw - corew, "shell slabs: coverage mismatch")
        cells = int(words)

    end subroutine s_amr_shell_slabs

    !> Intersect the shell-slab list with box [bl:bh]: the surviving clipped slabs in the SAME fixed order (each exchange side
    !! derives an identical list from replicated data, so empties drop symmetrically) plus their total cell count.
    impure subroutine s_amr_shell_clip(ns, sb1, se1, sb2, se2, sb3, se3, bl, bh, ms, tb1, te1, tb2, te2, tb3, te3, cells)

        integer, intent(in)  :: ns, sb1(6), se1(6), sb2(6), se2(6), sb3(6), se3(6), bl(3), bh(3)
        integer, intent(out) :: ms, tb1(6), te1(6), tb2(6), te2(6), tb3(6), te3(6), cells
        integer              :: s, l1, u1, l2, u2, l3, u3

        ms = 0; cells = 0
        do s = 1, ns
            l1 = max(sb1(s), bl(1)); u1 = min(se1(s), bh(1))
            l2 = max(sb2(s), bl(2)); u2 = min(se2(s), bh(2))
            l3 = max(sb3(s), bl(3)); u3 = min(se3(s), bh(3))
            if (l1 > u1 .or. l2 > u2 .or. l3 > u3) cycle
            ms = ms + 1
            tb1(ms) = l1; te1(ms) = u1; tb2(ms) = l2; te2(ms) = u2; tb3(ms) = l3; te3(ms) = u3
            cells = cells + (u1 - l1 + 1)*(u2 - l2 + 1)*(u3 - l3 + 1)
        end do

    end subroutine s_amr_shell_clip

    !> Validation arm (debug builds only): flood the current patch extent of amr_cg with quiet NaN BEFORE a clipped gather writes
    !! its shell, so a consumer read of any unshipped cell - core OR a missed shell slab - NaNs the ghost fill within a step
    !! (mandated by docs/documentation/amr_stepfill_ring_clip.md).
    impure subroutine s_amr_poison_patch_device(w1, w2, w3)

        use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        integer, intent(in) :: w1, w2, w3
        integer             :: i, g1, g2, g3
        real(stp)           :: nanv

        nanv = real(ieee_value(0._wp, ieee_quiet_nan), stp)
        $:GPU_PARALLEL_LOOP(collapse=4)
        do i = 1, sys_size
            do g3 = 0, w3
                do g2 = 0, w2
                    do g1 = 0, w1
                        amr_cg(i)%sf(g1, g2, g3) = nanv
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_poison_patch_device

    !> Ring-clipped runtime own-box copy (device): the owner's shell-slab / own-box intersections [tb:te] GLOBAL from q_coarse into
    !! amr_cg in the patch-local frame - no host round-trip, ONE fused kernel over the slab concatenation (the ghost-fill kernel's
    !! flat-index idiom). Same index map and direct stp assignment as the host path in s_amr_gather_coarse_patch.
    impure subroutine s_amr_gather_own_shell_device(q_coarse, ms, tb1, te1, tb2, te2, tb3, te3, o1, o2, o3)

        type(scalar_field), dimension(sys_size), intent(in) :: q_coarse
        integer, intent(in)                                 :: ms, tb1(6), te1(6), tb2(6), te2(6), tb3(6), te3(6), o1, o2, o3
        integer                                             :: lb1(6), le1(6), lb2(6), le2(6), lb3(6), le3(6), soff(6), scnt(6)
        integer                                             :: i, s, ss, g, r, n1, n2, g1, g2, g3, stot, coff1, coff2, coff3

        ! scalar/local copies: no host array may be referenced inside the device region (nvfortran/Cray demand it PRESENT)

        coff1 = amr_cpat_off(1); coff2 = amr_cpat_off(2); coff3 = amr_cpat_off(3)
        soff(1) = 0
        do s = 1, ms
            lb1(s) = tb1(s); le1(s) = te1(s); lb2(s) = tb2(s); le2(s) = te2(s); lb3(s) = tb3(s); le3(s) = te3(s)
            scnt(s) = (te1(s) - tb1(s) + 1)*(te2(s) - tb2(s) + 1)*(te3(s) - tb3(s) + 1)
            if (s < ms) soff(s + 1) = soff(s) + scnt(s)
        end do
        stot = soff(ms) + scnt(ms)
        $:GPU_PARALLEL_LOOP(collapse=2, copyin='[lb1, le1, lb2, le2, lb3, le3, soff, scnt]', &
                            & private='[s, ss, r, n1, n2, g1, g2, g3]')
        do i = 1, sys_size
            do g = 0, stot - 1
                s = 1  ! decode the flat index: ms <= 6, so a scan beats storing a per-cell slab map
                do ss = 2, ms
                    if (g >= soff(ss)) s = ss
                end do
                r = g - soff(s)
                n1 = le1(s) - lb1(s) + 1; n2 = le2(s) - lb2(s) + 1
                g1 = lb1(s) + mod(r, n1)
                g2 = lb2(s) + mod(r/n1, n2)
                g3 = lb3(s) + r/(n1*n2)
                amr_cg(i)%sf(g1 - coff1, g2 - coff2, g3 - coff3) = q_coarse(i)%sf(g1 - o1, g2 - o2, g3 - o3)
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_gather_own_shell_device

    !> Non-subcycle per-stage LEVEL-1 fill as one exchange WAVE (I2a, amr_plan_based_exchange.md): derive this stage's full (box,
    !! contributor) transfer set from the replicated caches, exchange one aggregated message per (peer, family) - F1 q_cons and,
    !! under non-polytropic QBMM, the F3 pb/mv twin - with all recvs posted first, then packs, then sends, then ONE waitall, and
    !! finally consume owned boxes in ascending slot order through the single amr_cg patch (own-box device copy + per-slab device
    !! unpack + ghost fill). The per-box path's operations, re-ordered: no per-box owner WAITALL, no per-box contributor flush, no
    !! F3 blocking sends. Level>=2 blocks keep the per-box F2 path (increment I3); the subcycle sites keep theirs (I8). Under
    !! MFC_DEBUG every slab carries the I1b identity header, verified at consume, and each received message length is checked
    !! against the plan.
    impure subroutine s_amr_stage_fill_wave(q_cons_coarse, pb_in, mv_in)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_coarse
        real(stp), dimension(:,:,:,:,:), intent(inout) :: pb_in, mv_in
        logical :: do_pbmv
        integer :: k, r, idx, ix, ip, owner, o1, o2, o3, qsz, psz, cellsz, tq, tp, nreq, qbase, pbase, ierr
        integer :: v1hi, v2hi, v3hi, plo(3), phi(3), crlo(3), crhi(3), bl(3), bh(3), boff
        integer :: clo(3), chi(3), nsh, msl, isl, scells
        integer :: shb1(6), she1(6), shb2(6), she2(6), shb3(6), she3(6), tb1(6), te1(6), tb2(6), te2(6), tb3(6), te3(6)

        if (amr_num_blocks <= 0) return
        @:ASSERT(.not. amr_subcycle, "stage-fill wave: the subcycle path keeps its per-box sites")
        @:ASSERT(amr_gsnd_n == 0, "stage-fill wave: the deferred gather-send pool must be drained")

        do_pbmv = qbmm .and. .not. polytropic
        cellsz = 0
        if (do_pbmv) cellsz = 2*nnode*nb
        o1 = start_idx(1); o2 = 0; o3 = 0
        if (n_glb > 0) o2 = start_idx(2)
        if (p_glb > 0) o3 = start_idx(3)
        tq = amr_tag_base(1) + int(mod(amr_mesh_epoch, 100_8))
        tp = amr_tag_base(3) + int(mod(amr_mesh_epoch, 100_8))

        call s_phase_tic(PH_GATHER)
        call s_phase_tic(PH_GWPLAN)
        ! block set changed: rebuild the cached overlap-rank lists BEFORE reading them (same lazy trigger as the per-box path)
        if (amr_seam_pairs_dirty .or. amr_seam_pairs_nblk /= amr_num_blocks) call s_amr_build_seam_pairs()
        if (.not. allocated(amr_fw_map)) then
            allocate (amr_fw_map(0:num_procs - 1), amr_fw_nx(0:num_procs - 1), amr_fw_pq(0:num_procs - 1), &
                      & amr_fw_pp(0:num_procs - 1))
            amr_fw_map = 0; amr_fw_nx = 0; amr_fw_pq = 0; amr_fw_pp = 0
        end if

        ! SEND side: for every level-1 box someone else owns, my coarse-range slice of its padded patch box. The per-box lag
        ! guard and slot selection run here so every rank still visits every level-1 slot once per stage, as before.
        amr_fw_snx = 0; amr_fw_snp = 0
        do k = 1, amr_num_blocks
            if (amr_block_level(k) /= 1) cycle
            call s_amr_select_slot(k)
            call s_amr_check_lag_clear()
            owner = amr_block_owner(k)
            if (owner == proc_rank) cycle
            plo(1) = amr_region_lo_all(1, k) - amr_cpat_mar; plo(2) = 0; plo(3) = 0
            if (n_glb > 0) plo(2) = amr_region_lo_all(2, k) - amr_cpat_mar
            if (p_glb > 0) plo(3) = amr_region_lo_all(3, k) - amr_cpat_mar
            v1hi = (amr_region_hi_all(1, k) - amr_region_lo_all(1, k)) + 2*amr_cpat_mar
            v2hi = 0; v3hi = 0
            if (n_glb > 0) v2hi = (amr_region_hi_all(2, k) - amr_region_lo_all(2, k)) + 2*amr_cpat_mar
            if (p_glb > 0) v3hi = (amr_region_hi_all(3, k) - amr_region_lo_all(3, k)) + 2*amr_cpat_mar
            phi(1) = plo(1) + v1hi; phi(2) = plo(2) + v2hi; phi(3) = plo(3) + v3hi
            call s_amr_rank_coarse_range(proc_rank, crlo, crhi)
            call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
            if (bl(1) > bh(1) .or. bl(2) > bh(2) .or. bl(3) > bh(3)) cycle
            ! ring clip (runtime q-only path): consumers of amr_cg read only the patch's hollow shell
            ! (docs/documentation/amr_stepfill_ring_clip.md), so ship only the shell's intersection with this rank's
            ! slice - as up to 6 sub-slab transfers, derived identically on both sides from replicated metadata. The
            ! pbmv gather keeps its full-box wire contract, so qbmm+non-polytropic runs stay unclipped (full slab).
            if (do_pbmv) then
                msl = 1
                tb1(1) = bl(1); te1(1) = bh(1); tb2(1) = bl(2); te2(1) = bh(2); tb3(1) = bl(3); te3(1) = bh(3)
            else
                clo = 0; chi = 0
                clo(1) = amr_region_lo_all(1, k) + 1; chi(1) = amr_region_hi_all(1, k) - 1
                if (n_glb > 0) then; clo(2) = amr_region_lo_all(2, k) + 1; chi(2) = amr_region_hi_all(2, k) - 1; end if
                if (p_glb > 0) then; clo(3) = amr_region_lo_all(3, k) + 1; chi(3) = amr_region_hi_all(3, k) - 1; end if
                call s_amr_shell_slabs(plo, phi, clo, chi, nsh, shb1, she1, shb2, she2, shb3, she3, scells)
                call s_amr_shell_clip(nsh, shb1, she1, shb2, she2, shb3, she3, bl, bh, msl, tb1, te1, tb2, te2, tb3, te3, scells)
                if (msl == 0) cycle
            end if
            do isl = 1, msl
                bl = [tb1(isl), tb2(isl), tb3(isl)]; bh = [te1(isl), te2(isl), te3(isl)]
                qsz = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                psz = cellsz*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                if (amr_fw_map(owner) == 0) then
                    amr_fw_snp = amr_fw_snp + 1
                    call s_amr_fw_szi(amr_fw_sprank, amr_fw_snp); call s_amr_fw_szi(amr_fw_sqsz, amr_fw_snp)
                    call s_amr_fw_szi(amr_fw_spsz, amr_fw_snp); call s_amr_fw_szi(amr_fw_snxp, amr_fw_snp)
                    call s_amr_fw_szi(amr_fw_sqbase, amr_fw_snp); call s_amr_fw_szi(amr_fw_spbase, amr_fw_snp)
                    amr_fw_map(owner) = amr_fw_snp
                    amr_fw_sprank(amr_fw_snp) = owner
                end if
                amr_fw_snx = amr_fw_snx + 1
                call s_amr_fw_szi(amr_fw_sblk, amr_fw_snx); call s_amr_fw_szi3(amr_fw_sbl, amr_fw_snx)
                call s_amr_fw_szi3(amr_fw_sbh, amr_fw_snx); call s_amr_fw_szi(amr_fw_spi, amr_fw_snx)
                call s_amr_fw_szi(amr_fw_sqo, amr_fw_snx); call s_amr_fw_szi(amr_fw_spo, amr_fw_snx)
                amr_fw_sblk(amr_fw_snx) = k; amr_fw_sbl(:,amr_fw_snx) = bl; amr_fw_sbh(:,amr_fw_snx) = bh
                amr_fw_spi(amr_fw_snx) = amr_fw_map(owner)
                amr_fw_sqo(amr_fw_snx) = amr_fw_pq(owner) + amr_fw_nx(owner)*XA_NH
                amr_fw_spo(amr_fw_snx) = amr_fw_pp(owner) + amr_fw_nx(owner)*XA_NH
                amr_fw_pq(owner) = amr_fw_pq(owner) + qsz
                amr_fw_pp(owner) = amr_fw_pp(owner) + psz
                amr_fw_nx(owner) = amr_fw_nx(owner) + 1
            end do
        end do
        qbase = 0; pbase = 0
        do ip = 1, amr_fw_snp
            r = amr_fw_sprank(ip)
            amr_fw_snxp(ip) = amr_fw_nx(r)
            amr_fw_sqsz(ip) = amr_fw_pq(r) + amr_fw_nx(r)*XA_NH
            amr_fw_spsz(ip) = amr_fw_pp(r) + amr_fw_nx(r)*XA_NH
            amr_fw_sqbase(ip) = qbase; qbase = qbase + amr_fw_sqsz(ip)
            amr_fw_spbase(ip) = pbase; pbase = pbase + amr_fw_spsz(ip)
            amr_fw_map(r) = 0; amr_fw_nx(r) = 0; amr_fw_pq(r) = 0; amr_fw_pp(r) = 0
        end do
        call s_amr_fw_szr(amr_fw_sq, qbase)
        if (do_pbmv) call s_amr_fw_szr(amr_fw_sp, pbase)

        ! RECV side: for every level-1 box I own, each listed contributor's slice (owner excluded - the own box is a device
        ! copy at consume). Both sides enumerate boxes ascending with per-rank running offsets, so the offsets agree.
        amr_fw_rnx = 0; amr_fw_rnp = 0
        do k = 1, amr_num_blocks
            if (amr_block_level(k) /= 1) cycle
            if (amr_block_owner(k) /= proc_rank) cycle
            plo(1) = amr_region_lo_all(1, k) - amr_cpat_mar; plo(2) = 0; plo(3) = 0
            if (n_glb > 0) plo(2) = amr_region_lo_all(2, k) - amr_cpat_mar
            if (p_glb > 0) plo(3) = amr_region_lo_all(3, k) - amr_cpat_mar
            v1hi = (amr_region_hi_all(1, k) - amr_region_lo_all(1, k)) + 2*amr_cpat_mar
            v2hi = 0; v3hi = 0
            if (n_glb > 0) v2hi = (amr_region_hi_all(2, k) - amr_region_lo_all(2, k)) + 2*amr_cpat_mar
            if (p_glb > 0) v3hi = (amr_region_hi_all(3, k) - amr_region_lo_all(3, k)) + 2*amr_cpat_mar
            phi(1) = plo(1) + v1hi; phi(2) = plo(2) + v2hi; phi(3) = plo(3) + v3hi
            ! ring clip: the shell is a per-box property; clip each contributor's slice against it (mirror of the
            ! send walk, so both sides derive the identical sub-slab list)
            if (.not. do_pbmv) then
                clo = 0; chi = 0
                clo(1) = amr_region_lo_all(1, k) + 1; chi(1) = amr_region_hi_all(1, k) - 1
                if (n_glb > 0) then; clo(2) = amr_region_lo_all(2, k) + 1; chi(2) = amr_region_hi_all(2, k) - 1; end if
                if (p_glb > 0) then; clo(3) = amr_region_lo_all(3, k) + 1; chi(3) = amr_region_hi_all(3, k) - 1; end if
                call s_amr_shell_slabs(plo, phi, clo, chi, nsh, shb1, she1, shb2, she2, shb3, she3, scells)
            end if
            do idx = 1, amr_ovl_gather_n(k)
                r = amr_ovl_gather(idx, k)
                if (r == proc_rank) cycle
                call s_amr_rank_coarse_range(r, crlo, crhi)
                call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
                if (do_pbmv) then
                    msl = 1
                    tb1(1) = bl(1); te1(1) = bh(1); tb2(1) = bl(2); te2(1) = bh(2); tb3(1) = bl(3); te3(1) = bh(3)
                else
                    call s_amr_shell_clip(nsh, shb1, she1, shb2, she2, shb3, she3, bl, bh, msl, tb1, te1, tb2, te2, tb3, te3, &
                                          & scells)
                    if (msl == 0) cycle
                end if
                do isl = 1, msl
                    bl = [tb1(isl), tb2(isl), tb3(isl)]; bh = [te1(isl), te2(isl), te3(isl)]
                    qsz = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                    psz = cellsz*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                    if (amr_fw_map(r) == 0) then
                        amr_fw_rnp = amr_fw_rnp + 1
                        call s_amr_fw_szi(amr_fw_rprank, amr_fw_rnp); call s_amr_fw_szi(amr_fw_rqsz, amr_fw_rnp)
                        call s_amr_fw_szi(amr_fw_rpsz, amr_fw_rnp); call s_amr_fw_szi(amr_fw_rnxp, amr_fw_rnp)
                        call s_amr_fw_szi(amr_fw_rqbase, amr_fw_rnp); call s_amr_fw_szi(amr_fw_rpbase, amr_fw_rnp)
                        amr_fw_map(r) = amr_fw_rnp
                        amr_fw_rprank(amr_fw_rnp) = r
                    end if
                    amr_fw_rnx = amr_fw_rnx + 1
                    call s_amr_fw_szi(amr_fw_rblk, amr_fw_rnx); call s_amr_fw_szi3(amr_fw_rbl, amr_fw_rnx)
                    call s_amr_fw_szi3(amr_fw_rbh, amr_fw_rnx); call s_amr_fw_szi(amr_fw_rpi, amr_fw_rnx)
                    call s_amr_fw_szi(amr_fw_rqo, amr_fw_rnx); call s_amr_fw_szi(amr_fw_rpo, amr_fw_rnx)
                    amr_fw_rblk(amr_fw_rnx) = k; amr_fw_rbl(:,amr_fw_rnx) = bl; amr_fw_rbh(:,amr_fw_rnx) = bh
                    amr_fw_rpi(amr_fw_rnx) = amr_fw_map(r)
                    amr_fw_rqo(amr_fw_rnx) = amr_fw_pq(r) + amr_fw_nx(r)*XA_NH
                    amr_fw_rpo(amr_fw_rnx) = amr_fw_pp(r) + amr_fw_nx(r)*XA_NH
                    amr_fw_pq(r) = amr_fw_pq(r) + qsz
                    amr_fw_pp(r) = amr_fw_pp(r) + psz
                    amr_fw_nx(r) = amr_fw_nx(r) + 1
                end do
            end do
        end do
        qbase = 0; pbase = 0
        do ip = 1, amr_fw_rnp
            r = amr_fw_rprank(ip)
            amr_fw_rnxp(ip) = amr_fw_nx(r)
            amr_fw_rqsz(ip) = amr_fw_pq(r) + amr_fw_nx(r)*XA_NH
            amr_fw_rpsz(ip) = amr_fw_pp(r) + amr_fw_nx(r)*XA_NH
            amr_fw_rqbase(ip) = qbase; qbase = qbase + amr_fw_rqsz(ip)
            amr_fw_rpbase(ip) = pbase; pbase = pbase + amr_fw_rpsz(ip)
            amr_fw_map(r) = 0; amr_fw_nx(r) = 0; amr_fw_pq(r) = 0; amr_fw_pp(r) = 0
        end do
        call s_amr_fw_szr(amr_fw_rq, qbase)
        if (do_pbmv) call s_amr_fw_szr(amr_fw_rp, pbase)
        nreq = (amr_fw_snp + amr_fw_rnp)*(1 + merge(1, 0, do_pbmv))
        call s_amr_fw_szi(amr_fw_req, nreq); call s_amr_fw_szi(amr_fw_reqw, nreq)
        call s_phase_toc(PH_GWPLAN)

        ! post ALL recvs, then pack ALL sends (device kernels into contiguous host pool slices), then post ALL sends, then one
        ! waitall (amr_plan_based_exchange.md order-of-operations rule). [amr-xa] records payload words only, so the family
        ! totals stay comparable to the per-box baseline; the message counts drop to the peer-pair count by design.
        nreq = 0
#ifdef MFC_MPI
        do ip = 1, amr_fw_rnp
            call s_xa_rec(XA_F1W_RCV, 2, amr_fw_rqsz(ip) - amr_fw_rnxp(ip)*XA_NH, tq)
            nreq = nreq + 1; amr_fw_reqw(nreq) = amr_fw_rqsz(ip)
            call MPI_IRECV(amr_fw_rq(amr_fw_rqbase(ip) + 1), amr_fw_rqsz(ip), mpi_p, amr_fw_rprank(ip), tq, MPI_COMM_WORLD, &
                           & amr_fw_req(nreq), ierr)
            if (do_pbmv) then
                call s_xa_rec(XA_F3W_RCV, 2, amr_fw_rpsz(ip) - amr_fw_rnxp(ip)*XA_NH, tp)
                nreq = nreq + 1; amr_fw_reqw(nreq) = amr_fw_rpsz(ip)
                call MPI_IRECV(amr_fw_rp(amr_fw_rpbase(ip) + 1), amr_fw_rpsz(ip), mpi_p, amr_fw_rprank(ip), tp, MPI_COMM_WORLD, &
                               & amr_fw_req(nreq), ierr)
            end if
        end do
#endif
        call s_phase_tic(PH_GWPACK)
        do ix = 1, amr_fw_snx
            bl = amr_fw_sbl(:,ix); bh = amr_fw_sbh(:,ix)
            qsz = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
            boff = amr_fw_sqbase(amr_fw_spi(ix)) + amr_fw_sqo(ix)
            if (XA_NH > 0) call s_xa_hdr_pack(amr_fw_sq(boff + 1:boff + XA_NH), XA_F1W_SND, amr_fw_sblk(ix), bl, bh)
            call s_amr_pack_box_device(q_cons_coarse, bl, bh, o1, o2, o3, amr_fw_sq(boff + XA_NH + 1:boff + XA_NH + qsz))
            if (do_pbmv) then
                psz = cellsz*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                boff = amr_fw_spbase(amr_fw_spi(ix)) + amr_fw_spo(ix)
                if (XA_NH > 0) call s_xa_hdr_pack(amr_fw_sp(boff + 1:boff + XA_NH), XA_F3W_SND, amr_fw_sblk(ix), bl, bh)
                call s_amr_pack_box_pbmv_device(pb_in, mv_in, bl, bh, o1, o2, o3, amr_fw_sp(boff + XA_NH + 1:boff + XA_NH + psz))
            end if
        end do
        call s_phase_toc(PH_GWPACK)
#ifdef MFC_MPI
        do ip = 1, amr_fw_snp
            call s_xa_rec(XA_F1W_SND, 1, amr_fw_sqsz(ip) - amr_fw_snxp(ip)*XA_NH, tq)
            nreq = nreq + 1; amr_fw_reqw(nreq) = -1
            call MPI_ISEND(amr_fw_sq(amr_fw_sqbase(ip) + 1), amr_fw_sqsz(ip), mpi_p, amr_fw_sprank(ip), tq, MPI_COMM_WORLD, &
                           & amr_fw_req(nreq), ierr)
            if (do_pbmv) then
                call s_xa_rec(XA_F3W_SND, 1, amr_fw_spsz(ip) - amr_fw_snxp(ip)*XA_NH, tp)
                nreq = nreq + 1; amr_fw_reqw(nreq) = -1
                call MPI_ISEND(amr_fw_sp(amr_fw_spbase(ip) + 1), amr_fw_spsz(ip), mpi_p, amr_fw_sprank(ip), tp, MPI_COMM_WORLD, &
                               & amr_fw_req(nreq), ierr)
            end if
        end do
        call s_phase_tic(PH_GWWAIT)
        if (nreq > 0) then
#ifdef MFC_DEBUG
            block
                integer :: st(MPI_STATUS_SIZE, nreq), gotw, q
                call MPI_WAITALL(nreq, amr_fw_req, st, ierr)
                do q = 1, nreq
                    if (amr_fw_reqw(q) < 0) cycle
                    call MPI_GET_COUNT(st(:,q), mpi_p, gotw, ierr)
                    @:ASSERT(gotw == amr_fw_reqw(q), "stage-fill wave: received message length differs from the plan")
                end do
            end block
#else
            call MPI_WAITALL(nreq, amr_fw_req, MPI_STATUSES_IGNORE, ierr)
#endif
        end if
        call s_phase_toc(PH_GWWAIT)
#endif
        call s_phase_toc(PH_GATHER)

        ! CONSUME, ascending slot order: per owned box, patch frame + own-box device copy + per-slab device unpack (recv
        ! transfers were appended box-major, so each box's slabs are the next contiguous run), then the ghost fills - the
        ! same operations the per-box path ran, minus its rendezvous.
        ix = 1
        do k = 1, amr_num_blocks
            if (amr_block_level(k) /= 1) cycle
            call s_amr_select_slot(k)
            if (.not. amr_rank_owns_block) cycle
            call s_phase_tic(PH_GATHER)
            amr_cpat_off = 0
            amr_cpat_off(1) = amr_region_lo_all(1, k) - amr_cpat_mar
            if (n_glb > 0) amr_cpat_off(2) = amr_region_lo_all(2, k) - amr_cpat_mar
            if (p_glb > 0) amr_cpat_off(3) = amr_region_lo_all(3, k) - amr_cpat_mar
            v1hi = (amr_region_hi_all(1, k) - amr_region_lo_all(1, k)) + 2*amr_cpat_mar
            v2hi = 0; v3hi = 0
            if (n_glb > 0) v2hi = (amr_region_hi_all(2, k) - amr_region_lo_all(2, k)) + 2*amr_cpat_mar
            if (p_glb > 0) v3hi = (amr_region_hi_all(3, k) - amr_region_lo_all(3, k)) + 2*amr_cpat_mar
            plo = amr_cpat_off
            phi(1) = plo(1) + v1hi; phi(2) = plo(2) + v2hi; phi(3) = plo(3) + v3hi
            call s_amr_rank_coarse_range(proc_rank, crlo, crhi)
            call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
            if (do_pbmv) then
                call s_amr_gather_own_box_device(q_cons_coarse, bl, bh, o1, o2, o3)
                call s_amr_gather_own_box_pbmv_device(pb_in, mv_in, bl, bh, o1, o2, o3)
            else
#ifdef MFC_DEBUG
                ! validation arm: flood the patch with NaN BEFORE the clipped writes, so a consumer read of any
                ! unshipped cell - core OR a missed shell slab - NaNs the ghost fill within a step
                call s_amr_poison_patch_device(v1hi, v2hi, v3hi)
#endif
                clo = 0; chi = 0
                clo(1) = amr_region_lo_all(1, k) + 1; chi(1) = amr_region_hi_all(1, k) - 1
                if (n_glb > 0) then; clo(2) = amr_region_lo_all(2, k) + 1; chi(2) = amr_region_hi_all(2, k) - 1; end if
                if (p_glb > 0) then; clo(3) = amr_region_lo_all(3, k) + 1; chi(3) = amr_region_hi_all(3, k) - 1; end if
                call s_amr_shell_slabs(plo, phi, clo, chi, nsh, shb1, she1, shb2, she2, shb3, she3, scells)
                call s_amr_shell_clip(nsh, shb1, she1, shb2, she2, shb3, she3, bl, bh, msl, tb1, te1, tb2, te2, tb3, te3, scells)
                if (msl > 0) call s_amr_gather_own_shell_device(q_cons_coarse, msl, tb1, te1, tb2, te2, tb3, te3, o1, o2, o3)
            end if
            do while (ix <= amr_fw_rnx)
                if (amr_fw_rblk(ix) /= k) exit
                bl = amr_fw_rbl(:,ix); bh = amr_fw_rbh(:,ix)
                qsz = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                boff = amr_fw_rqbase(amr_fw_rpi(ix)) + amr_fw_rqo(ix)
                if (XA_NH > 0) call s_xa_hdr_check(amr_fw_rq(boff + 1:boff + XA_NH), XA_F1W_SND, k, bl, bh)
                call s_amr_unpack_box_device(bl, bh, amr_fw_rq(boff + XA_NH + 1:boff + XA_NH + qsz))
                if (do_pbmv) then
                    psz = cellsz*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                    boff = amr_fw_rpbase(amr_fw_rpi(ix)) + amr_fw_rpo(ix)
                    if (XA_NH > 0) call s_xa_hdr_check(amr_fw_rp(boff + 1:boff + XA_NH), XA_F3W_SND, k, bl, bh)
                    call s_amr_unpack_box_pbmv_device(bl, bh, amr_fw_rp(boff + XA_NH + 1:boff + XA_NH + psz))
                end if
                ix = ix + 1
            end do
            call s_phase_toc(PH_GATHER)
            if (rank_time_wrt) call s_rank_time_tic()
            call s_amr_cov_note_fill()
            call s_phase_tic(PH_GFILL)
            call s_amr_fill_fine_ghosts_cons(amr_cg, amr_loc_of(amr_cur))
            call s_phase_toc(PH_GFILL)
            if (do_pbmv) call s_amr_fill_fine_ghosts_pbmv(amr_cg_pb, amr_cg_mv, amr_slots(amr_cur)%pb_f%sf, &
                & amr_slots(amr_cur)%mv_f%sf)
            if (rank_time_wrt) call s_rank_time_toc()
        end do
        @:ASSERT(ix == amr_fw_rnx + 1, "stage-fill wave: unconsumed recv transfers")

    end subroutine s_amr_stage_fill_wave

    !> The parent-fill wave's per-box transfer list in the PATCH-LOCAL frame: the padded patch's hollow-shell slabs (the same
    !! dead-byte proof as the stepfill clip - the runtime consumer is the same amr_cg ghost fill, so the open interior of the parent
    !! footprint [mar+1, w-mar-1] never ships), or the single full patch when non-polytropic QBMM keeps the full-box contract. Send
    !! walk, recv walk, and consume all derive the list HERE, so the wire layout cannot drift between sides.
    impure subroutine s_amr_parent_shell(w1, w2, w3, full, msl, tb1, te1, tb2, te2, tb3, te3)

        integer, intent(in)  :: w1, w2, w3
        logical, intent(in)  :: full
        integer, intent(out) :: msl, tb1(6), te1(6), tb2(6), te2(6), tb3(6), te3(6)
        integer              :: clo(3), chi(3), scells

        if (full) then
            msl = 1
            tb1(1) = 0; te1(1) = w1; tb2(1) = 0; te2(1) = w2; tb3(1) = 0; te3(1) = w3
        else
            clo = 0; chi = 0
            clo(1) = amr_cpat_mar + 1; chi(1) = w1 - amr_cpat_mar - 1
            if (n_glb > 0) then; clo(2) = amr_cpat_mar + 1; chi(2) = w2 - amr_cpat_mar - 1; end if
            if (p_glb > 0) then; clo(3) = amr_cpat_mar + 1; chi(3) = w3 - amr_cpat_mar - 1; end if
            call s_amr_shell_slabs([0, 0, 0], [w1, w2, w3], clo, chi, msl, tb1, te1, tb2, te2, tb3, te3, scells)
        end if

    end subroutine s_amr_parent_shell

    !> Per-step LEVEL-lev fill as one exchange WAVE (I3, amr_plan_based_exchange.md): the F2 parent gather for every level-lev block
    !! in one aggregated exchange - each split child is its s_amr_parent_shell transfer list (ring-clipped shell slabs, or one full
    !! patch under the pbmv contract) from its parent's owner to its own owner, so the plan is a pair list, not an overlap map. Same
    !! skeleton as s_amr_stage_fill_wave (whose scratch arrays it reuses - the two never overlap in time): plans from replicated
    !! metadata (f_amr_parent_block + s_amr_parent_foot + amr_block_owner ONLY - the per-owner mirrors lag and are empty on
    !! non-owners), recvs-packs-sends-one-WAITALL, box-major consume through the single amr_cg. Called per level ASCENDING, so a
    !! level-(lev-1) parent's own ghost fill is complete before this wave reads its interior (the same parent-before-child guarantee
    !! slot order gave the per-box loop). Co-located parent-child stays a consume-phase device copy with no wire transfer. The
    !! regrid keeps the chunked F2 path; the subcycle keeps its per-box sites (I8); init/static keep the per-box
    !! s_amr_gather_from_parent.
    impure subroutine s_amr_parent_fill_wave(lev)

        integer, intent(in) :: lev
        integer             :: k, r, ix, ip, pblk, powner, cowner, boxsz, tq, nreq, qbase, ierr
        integer             :: w1, w2, w3, plo(3), phi(3), boff, bl(3), bh(3)
        integer             :: msl, isl
        integer             :: tb1(6), te1(6), tb2(6), te2(6), tb3(6), te3(6)
        logical             :: do_pbmv

        if (amr_num_blocks <= 0) return
        @:ASSERT(.not. amr_subcycle, "parent-fill wave: the subcycle path keeps its per-box sites")
        @:ASSERT(amr_gsnd_n == 0, "parent-fill wave: the deferred gather-send pool must be drained")
        do_pbmv = qbmm .and. .not. polytropic

        tq = amr_tag_base(2) + int(mod(amr_mesh_epoch, 100_8))

        call s_phase_tic(PH_GATHER)
        ! SEND side: every level-lev block whose parent I own but whose child-owner is another rank. The per-box lag guard and
        ! slot selection run here so every rank still visits every level-lev slot once per stage, as before.
        amr_fw_snx = 0; amr_fw_snp = 0
        if (.not. allocated(amr_fw_map)) then
            allocate (amr_fw_map(0:num_procs - 1), amr_fw_nx(0:num_procs - 1), amr_fw_pq(0:num_procs - 1), &
                      & amr_fw_pp(0:num_procs - 1))
            amr_fw_map = 0; amr_fw_nx = 0; amr_fw_pq = 0; amr_fw_pp = 0
        end if
        do k = 1, amr_num_blocks
            if (amr_block_level(k) /= lev) cycle
            call s_amr_select_slot(k)
            call s_amr_check_lag_clear()
            pblk = f_amr_parent_block(k)
            powner = amr_block_owner(pblk); cowner = amr_block_owner(k)
            if (powner == cowner .or. powner /= proc_rank) cycle
            call s_amr_parent_foot(k, pblk, plo, phi)
            w1 = (phi(1) - plo(1)) + 2*amr_cpat_mar
            w2 = 0; w3 = 0
            if (n_glb > 0) w2 = (phi(2) - plo(2)) + 2*amr_cpat_mar
            if (p_glb > 0) w3 = (phi(3) - plo(3)) + 2*amr_cpat_mar
            call s_amr_parent_shell(w1, w2, w3, do_pbmv, msl, tb1, te1, tb2, te2, tb3, te3)
            do isl = 1, msl
                bl = [tb1(isl), tb2(isl), tb3(isl)]; bh = [te1(isl), te2(isl), te3(isl)]
                boxsz = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                if (amr_fw_map(cowner) == 0) then
                    amr_fw_snp = amr_fw_snp + 1
                    call s_amr_fw_szi(amr_fw_sprank, amr_fw_snp); call s_amr_fw_szi(amr_fw_sqsz, amr_fw_snp)
                    call s_amr_fw_szi(amr_fw_snxp, amr_fw_snp); call s_amr_fw_szi(amr_fw_sqbase, amr_fw_snp)
                    amr_fw_map(cowner) = amr_fw_snp
                    amr_fw_sprank(amr_fw_snp) = cowner
                end if
                amr_fw_snx = amr_fw_snx + 1
                call s_amr_fw_szi(amr_fw_sblk, amr_fw_snx); call s_amr_fw_szi3(amr_fw_sbl, amr_fw_snx)
                call s_amr_fw_szi3(amr_fw_sbh, amr_fw_snx); call s_amr_fw_szi(amr_fw_spi, amr_fw_snx)
                call s_amr_fw_szi(amr_fw_sqo, amr_fw_snx)
                amr_fw_sblk(amr_fw_snx) = k; amr_fw_sbl(:,amr_fw_snx) = bl; amr_fw_sbh(:,amr_fw_snx) = bh
                amr_fw_spi(amr_fw_snx) = amr_fw_map(cowner)
                amr_fw_sqo(amr_fw_snx) = amr_fw_pq(cowner) + amr_fw_nx(cowner)*XA_NH
                amr_fw_pq(cowner) = amr_fw_pq(cowner) + boxsz
                amr_fw_nx(cowner) = amr_fw_nx(cowner) + 1
            end do
        end do
        qbase = 0
        do ip = 1, amr_fw_snp
            r = amr_fw_sprank(ip)
            amr_fw_snxp(ip) = amr_fw_nx(r)
            amr_fw_sqsz(ip) = amr_fw_pq(r) + amr_fw_nx(r)*XA_NH
            amr_fw_sqbase(ip) = qbase; qbase = qbase + amr_fw_sqsz(ip)
            amr_fw_map(r) = 0; amr_fw_nx(r) = 0; amr_fw_pq(r) = 0
        end do
        call s_amr_fw_szr(amr_fw_sq, qbase)

        ! RECV side: every level-lev block I own whose parent lives on another rank - the box's shell-slab transfers (or its
        ! one full-patch transfer under the pbmv contract). Both sides enumerate boxes ascending, slabs in the fixed
        ! s_amr_parent_shell order, with per-rank running offsets, so the wire layout agrees with no metadata exchange.
        amr_fw_rnx = 0; amr_fw_rnp = 0
        do k = 1, amr_num_blocks
            if (amr_block_level(k) /= lev) cycle
            if (amr_block_owner(k) /= proc_rank) cycle
            pblk = f_amr_parent_block(k)
            powner = amr_block_owner(pblk)
            if (powner == proc_rank) cycle
            call s_amr_parent_foot(k, pblk, plo, phi)
            w1 = (phi(1) - plo(1)) + 2*amr_cpat_mar
            w2 = 0; w3 = 0
            if (n_glb > 0) w2 = (phi(2) - plo(2)) + 2*amr_cpat_mar
            if (p_glb > 0) w3 = (phi(3) - plo(3)) + 2*amr_cpat_mar
            call s_amr_parent_shell(w1, w2, w3, do_pbmv, msl, tb1, te1, tb2, te2, tb3, te3)
            do isl = 1, msl
                bl = [tb1(isl), tb2(isl), tb3(isl)]; bh = [te1(isl), te2(isl), te3(isl)]
                boxsz = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                if (amr_fw_map(powner) == 0) then
                    amr_fw_rnp = amr_fw_rnp + 1
                    call s_amr_fw_szi(amr_fw_rprank, amr_fw_rnp); call s_amr_fw_szi(amr_fw_rqsz, amr_fw_rnp)
                    call s_amr_fw_szi(amr_fw_rnxp, amr_fw_rnp); call s_amr_fw_szi(amr_fw_rqbase, amr_fw_rnp)
                    amr_fw_map(powner) = amr_fw_rnp
                    amr_fw_rprank(amr_fw_rnp) = powner
                end if
                amr_fw_rnx = amr_fw_rnx + 1
                call s_amr_fw_szi(amr_fw_rblk, amr_fw_rnx); call s_amr_fw_szi3(amr_fw_rbl, amr_fw_rnx)
                call s_amr_fw_szi3(amr_fw_rbh, amr_fw_rnx); call s_amr_fw_szi(amr_fw_rpi, amr_fw_rnx)
                call s_amr_fw_szi(amr_fw_rqo, amr_fw_rnx)
                amr_fw_rblk(amr_fw_rnx) = k; amr_fw_rbl(:,amr_fw_rnx) = bl; amr_fw_rbh(:,amr_fw_rnx) = bh
                amr_fw_rpi(amr_fw_rnx) = amr_fw_map(powner)
                amr_fw_rqo(amr_fw_rnx) = amr_fw_pq(powner) + amr_fw_nx(powner)*XA_NH
                amr_fw_pq(powner) = amr_fw_pq(powner) + boxsz
                amr_fw_nx(powner) = amr_fw_nx(powner) + 1
            end do
        end do
        qbase = 0
        do ip = 1, amr_fw_rnp
            r = amr_fw_rprank(ip)
            amr_fw_rnxp(ip) = amr_fw_nx(r)
            amr_fw_rqsz(ip) = amr_fw_pq(r) + amr_fw_nx(r)*XA_NH
            amr_fw_rqbase(ip) = qbase; qbase = qbase + amr_fw_rqsz(ip)
            amr_fw_map(r) = 0; amr_fw_nx(r) = 0; amr_fw_pq(r) = 0
        end do
        call s_amr_fw_szr(amr_fw_rq, qbase)
        nreq = amr_fw_snp + amr_fw_rnp
        call s_amr_fw_szi(amr_fw_req, nreq); call s_amr_fw_szi(amr_fw_reqw, nreq)

        nreq = 0
#ifdef MFC_MPI
        do ip = 1, amr_fw_rnp
            call s_xa_rec(XA_F2W_RCV, 2, amr_fw_rqsz(ip) - amr_fw_rnxp(ip)*XA_NH, tq)
            nreq = nreq + 1; amr_fw_reqw(nreq) = amr_fw_rqsz(ip)
            call MPI_IRECV(amr_fw_rq(amr_fw_rqbase(ip) + 1), amr_fw_rqsz(ip), mpi_p, amr_fw_rprank(ip), tq, MPI_COMM_WORLD, &
                           & amr_fw_req(nreq), ierr)
        end do
#endif
        ! pack: the parent-patch pack kernel reads amr_cpat_off from module scope, so set the CHILD's frame per transfer
        ! (the consume loop recomputes it per box; phases are sequential, so the global is single-writer at any time).
        ! sbl/sbh hold the transfer's PATCH-LOCAL slab bounds; the frame comes from the box's parent foot.
        do ix = 1, amr_fw_snx
            k = amr_fw_sblk(ix)
            call s_amr_parent_foot(k, f_amr_parent_block(k), plo, phi)
            amr_cpat_off = 0
            amr_cpat_off(1) = plo(1) - amr_cpat_mar
            if (n_glb > 0) amr_cpat_off(2) = plo(2) - amr_cpat_mar
            if (p_glb > 0) amr_cpat_off(3) = plo(3) - amr_cpat_mar
            bl = amr_fw_sbl(:,ix); bh = amr_fw_sbh(:,ix)
            boxsz = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
            boff = amr_fw_sqbase(amr_fw_spi(ix)) + amr_fw_sqo(ix)
            call s_amr_pack_parent_box_device_cons(amr_loc_of(f_amr_parent_block(k)), bl, bh, &
                                                   & amr_fw_sq(boff + XA_NH + 1:boff + XA_NH + boxsz))
            if (XA_NH > 0) call s_xa_hdr_pack(amr_fw_sq(boff + 1:boff + XA_NH), XA_F2W_SND, k, bl, bh)
        end do
#ifdef MFC_MPI
        do ip = 1, amr_fw_snp
            call s_xa_rec(XA_F2W_SND, 1, amr_fw_sqsz(ip) - amr_fw_snxp(ip)*XA_NH, tq)
            nreq = nreq + 1; amr_fw_reqw(nreq) = -1
            call MPI_ISEND(amr_fw_sq(amr_fw_sqbase(ip) + 1), amr_fw_sqsz(ip), mpi_p, amr_fw_sprank(ip), tq, MPI_COMM_WORLD, &
                           & amr_fw_req(nreq), ierr)
        end do
        if (nreq > 0) then
#ifdef MFC_DEBUG
            block
                integer :: st(MPI_STATUS_SIZE, nreq), gotw, q
                call MPI_WAITALL(nreq, amr_fw_req, st, ierr)
                do q = 1, nreq
                    if (amr_fw_reqw(q) < 0) cycle
                    call MPI_GET_COUNT(st(:,q), mpi_p, gotw, ierr)
                    @:ASSERT(gotw == amr_fw_reqw(q), "parent-fill wave: received message length differs from the plan")
                end do
            end block
#else
            call MPI_WAITALL(nreq, amr_fw_req, MPI_STATUSES_IGNORE, ierr)
#endif
        end if
#endif
        call s_phase_toc(PH_GATHER)

        ! CONSUME, ascending slot order: per owned level-lev box, patch frame + the shell-slab parent copies (co-located) or
        ! the box's received transfers, then the ghost fills - the per-box path's operations, minus its rendezvous.
        ix = 1
        do k = 1, amr_num_blocks
            if (amr_block_level(k) /= lev) cycle
            call s_amr_select_slot(k)
            if (.not. amr_rank_owns_block) cycle
            call s_phase_tic(PH_GATHER)
            pblk = f_amr_parent_block(k)
            call s_amr_parent_foot(k, pblk, plo, phi)
            amr_cpat_off = 0
            amr_cpat_off(1) = plo(1) - amr_cpat_mar
            if (n_glb > 0) amr_cpat_off(2) = plo(2) - amr_cpat_mar
            if (p_glb > 0) amr_cpat_off(3) = plo(3) - amr_cpat_mar
            w1 = (phi(1) - plo(1)) + 2*amr_cpat_mar
            w2 = 0; w3 = 0
            if (n_glb > 0) w2 = (phi(2) - plo(2)) + 2*amr_cpat_mar
            if (p_glb > 0) w3 = (phi(3) - plo(3)) + 2*amr_cpat_mar
#ifdef MFC_DEBUG
            ! validation arm (mirror of the stepfill clip): NaN-flood the patch before the shell writes land, so a consumer
            ! read of any unshipped cell - the clipped core or a missed slab - NaNs the ghost fill within a step
            if (.not. do_pbmv) call s_amr_poison_patch_device(w1, w2, w3)
#endif
            if (amr_block_owner(pblk) == proc_rank) then
                call s_amr_parent_shell(w1, w2, w3, do_pbmv, msl, tb1, te1, tb2, te2, tb3, te3)
                do isl = 1, msl
                    call s_amr_copy_parent_box_cons(amr_loc_of(pblk), [tb1(isl), tb2(isl), tb3(isl)], [te1(isl), te2(isl), &
                                                    & te3(isl)])
                end do
            else
                @:ASSERT(ix <= amr_fw_rnx .and. amr_fw_rblk(ix) == k, "parent-fill wave: missing recv transfer")
                do while (ix <= amr_fw_rnx)
                    if (amr_fw_rblk(ix) /= k) exit
                    bl = amr_fw_rbl(:,ix); bh = amr_fw_rbh(:,ix)
                    boxsz = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                    boff = amr_fw_rqbase(amr_fw_rpi(ix)) + amr_fw_rqo(ix)
                    if (XA_NH > 0) call s_xa_hdr_check(amr_fw_rq(boff + 1:boff + XA_NH), XA_F2W_SND, k, bl, bh)
                    call s_amr_unpack_parent_box_device(bl, bh, amr_fw_rq(boff + XA_NH + 1:boff + XA_NH + boxsz))
                    ix = ix + 1
                end do
            end if
            call s_phase_toc(PH_GATHER)
            if (rank_time_wrt) call s_rank_time_tic()
            call s_amr_cov_note_fill()
            call s_phase_tic(PH_GFILL)
            call s_amr_fill_fine_ghosts_cons(amr_cg, amr_loc_of(amr_cur))
            call s_phase_toc(PH_GFILL)
            if (qbmm .and. .not. polytropic) call s_amr_fill_fine_ghosts_pbmv(amr_cg_pb, amr_cg_mv, amr_slots(amr_cur)%pb_f%sf, &
                & amr_slots(amr_cur)%mv_f%sf)
            if (rank_time_wrt) call s_rank_time_toc()
        end do
        @:ASSERT(ix == amr_fw_rnx + 1, "parent-fill wave: unconsumed recv transfers")

    end subroutine s_amr_parent_fill_wave

    !> High-water sizing for the wave's plan scratch. Callers size-then-write at APPEND time, so a grow must preserve the entries
    !! already appended this wave.
    impure subroutine s_amr_fw_szi(a, n)

        integer, allocatable, intent(inout) :: a(:)
        integer, intent(in)                 :: n
        integer, allocatable                :: tmp(:)

        if (.not. allocated(a)) then
            allocate (a(max(n, 64)))
            return
        end if
        if (size(a) >= n) return
        call move_alloc(a, tmp)
        allocate (a(max(n, 2*size(tmp))))
        a(1:size(tmp)) = tmp

    end subroutine s_amr_fw_szi

    impure subroutine s_amr_fw_szi3(a, n)

        integer, allocatable, intent(inout) :: a(:,:)
        integer, intent(in)                 :: n
        integer, allocatable                :: tmp(:,:)

        if (.not. allocated(a)) then
            allocate (a(3, max(n, 64)))
            return
        end if
        if (size(a, 2) >= n) return
        call move_alloc(a, tmp)
        allocate (a(3, max(n, 2*size(tmp, 2))))
        a(:,1:size(tmp, 2)) = tmp

    end subroutine s_amr_fw_szi3

    !> Wire pools: preserving on grow (the F5 waves append debug header slots incrementally; the other waves size once).
    impure subroutine s_amr_fw_szr(a, n)

        real(wp), allocatable, intent(inout) :: a(:)
        integer, intent(in)                  :: n
        real(wp), allocatable                :: tmp(:)

        if (.not. allocated(a)) then
            allocate (a(max(n, 64)))
            return
        end if
        if (size(a) >= n) return
        call move_alloc(a, tmp)
        allocate (a(max(n, 2*size(tmp))))
        a(1:size(tmp)) = tmp

    end subroutine s_amr_fw_szr

    !> ADVANCE phase of a fine RK stage: fine RHS + RK update (+ QBMM/6eq/IB) for the current block. Owner-only. Reads the block's
    !! ghost shell (coarse prolong + fine-fine halo already applied by the fill + halo phases).
    !> One fine-block RK stage = RHS pass then RK pass. Fused wrapper: the AMR fine blocks (m_time_steppers) call this so their
    !! rhs+rk stay back-to-back (byte-identical). The coexist tile path (s_l0_advance_stage) instead calls the two passes directly
    !! with the reflux-delta copy-back interposed between them, so the corrected coarse rhs reaches the tile before its RK update.
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
        ! 2a: this block's computed prim vars (mom, E) were already produced by the stage-top batched conversion;
        ! land them and let s_compute_rhs skip its per-block conversion. L0 tile slots (level 0) are not in the
        ! batch and keep the per-block conversion.
        if (amr_prim_batch .and. amr_block_level(amr_cur) >= 1) then
            call s_amr_prim_load(q_prim_qp%vf, amr_loc_of(amr_cur))
            amr_prim_preloaded = .true.
        end if
        if (qbmm .and. .not. polytropic) then
            ! the block's OWN side-state and rhs scratch: the coarse pb_in/rhs_pb must not be touched at fine indices (the coarse
            ! stage consumes them after this fine stage)
            call s_compute_rhs(amr_cons_br, q_T_sf, q_prim_b, bc_type, rhs_b, amr_slots(amr_cur)%pb_f%sf, amr_rhs_pb_f, &
                               & amr_slots(amr_cur)%mv_f%sf, amr_rhs_mv_f, t_step, s)
        else
            call s_compute_rhs(amr_cons_br, q_T_sf, q_prim_b, bc_type, rhs_b, pb_in, rhs_pb, mv_in, rhs_mv, t_step, s)
        end if
        amr_prim_preloaded = .false.
        call s_amr_br_store(amr_loc_of(amr_cur))
        call s_phase_toc(PH_RHS)
        call s_phase_tic(PH_SWAP)  ! the other half of the swap pair - keep the bracket symmetric
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
        ! RK stage update (device kernel; mirror of the coarse form - under IGR the rhs already embeds dt, matching the coarse igr
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

    !> Per-block SETUP for the transposed subcycle advance (amr_subcycle): exchange valid coarse ghosts, gather+prolong the selected
    !! block's two time-lerp ghost sources (parent t^n in q_ghost_a, t^{n+1} in q_ghost_b), and zero its flux registers. The
    !! collective exchanges/gathers run on ALL ranks; the owner-only fills and register-zero are guarded. Called once per level-1
    !! block before the transposed stage loop (which reuses the prepared ghost sources every substep).
    impure subroutine s_amr_subcycle_setup_block(q_old, q_new, pb_old, mv_old, pb_in, mv_in)

        type(scalar_field), dimension(sys_size), intent(inout)                                  :: q_old, q_new
        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(in) :: pb_old, mv_old
        real(stp), dimension(:,:,:,:,:), intent(inout)                                          :: pb_in, mv_in

        ! valid coarse CONS ghosts on both lerp sources (ALL ranks call: pairwise halo); the exchanged t^n / t^{n+1} ghost layers
        ! make the prolonged block-boundary ghosts correct even at rank boundaries

        ! the two lerp sources' coarse halos are LOOP-INVARIANT over the setup loop (both read, neither written), so the
        ! exchanges are hoisted to s_amr_advance_fine_subcycle_all - same defect as the lock-step path, doubled for two fields

        call s_amr_check_lag_clear()  ! EVERY rank: non-owner bubbles can reach the block across a seam

        ! fine-level distribution: gather each lerp source's coarse patch (collective - ALL ranks) then prolong its ghost shell on
        ! the owner. Interleaved so the single amr_cg buffer is consumed by the fill before the next gather overwrites it.
        ! non-polytropic QBMM: the pb/mv ghost shell gets the same two-source time-lerp treatment, gathered + filled INTERLEAVED
        ! with
        ! q_cons so the single amr_cg_pb/mv buffer is consumed by each fill before the next gather overwrites it. Gathers collective
        ! (ALL ranks - P2P); fills owner-only.
        call s_amr_gather_coarse_patch(q_old, .true.)
        call s_amr_gather_send_flush()  ! keep this site's original blocking semantics
        if (qbmm .and. .not. polytropic) call s_amr_gather_coarse_patch_pbmv(pb_old, mv_old, .true.)
        if (amr_rank_owns_block) call s_amr_fill_fine_ghosts_gsta(amr_cg, amr_loc_of(amr_cur))
        if (amr_rank_owns_block .and. qbmm .and. .not. polytropic) call s_amr_fill_fine_ghosts_pbmv(amr_cg_pb, amr_cg_mv, &
            & amr_slots(amr_cur)%pb_ghost_a%sf, amr_slots(amr_cur)%mv_ghost_a%sf)
        call s_amr_gather_coarse_patch(q_new, .true.)
        call s_amr_gather_send_flush()  ! keep this site's original blocking semantics
        if (qbmm .and. .not. polytropic) call s_amr_gather_coarse_patch_pbmv(pb_in, mv_in, .true.)
        if (amr_rank_owns_block) call s_amr_fill_fine_ghosts_gstb(amr_cg, amr_loc_of(amr_cur))
        if (amr_rank_owns_block .and. qbmm .and. .not. polytropic) call s_amr_fill_fine_ghosts_pbmv(amr_cg_pb, amr_cg_mv, &
            & amr_slots(amr_cur)%pb_ghost_b%sf, amr_slots(amr_cur)%mv_ghost_b%sf)
        if (.not. amr_rank_owns_block) return

        ! registers accumulate over all six stages of the transposed loop, so zero them once at setup (the stage-1 overwrite trick
        ! cannot span two substeps)
        call s_amr_zero_fine_registers()

    end subroutine s_amr_subcycle_setup_block

    !> Subcycled fine advance (amr_subcycle) over ALL level-1 blocks, TRANSPOSED: instead of each block running its full 2x3-stage
    !! subcycle in turn, every same-level block advances stage-by-stage in LOCKSTEP with the block-to-block fine-fine seam halo
    !! (s_amr_fine_fine_halo) interposed between the ghost lerp and the RHS at each stage. That makes max_grid_size-tiled ADJACENT
    !! sub-blocks (which appear at np>1 when a feature exceeds a rank's slot) compute a MATCHING shared-face flux, so the subcycle
    !! conserves at the seam - the per-block order did not run the halo and leaked there. Two dt/2 SSP-RK3 substeps AFTER the coarse
    !! step: q_old/q_new are the coarse t^n / t^{n+1} states; each stage's ghosts are the linear time interpolation at stage time
    !! theta = (substep-1 + c_s)/2 with SSP-RK3 abscissae c = [0, 1, 1/2]. Level-1 blocks drive their level-2 children per substep
    !! (s_amr_advance_children), which applies this same transposed shape at every deeper level, so L2-L2 seams are reconciled by
    !! the level-filtered halo too. A single owned level-1 block is byte-identical to the old per-block subcycle (the halo is a
    !! no-op with < 2 adjacent same-level blocks).
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

        ! valid coarse CONS ghosts on BOTH lerp sources, ONCE for the whole setup loop (ALL ranks call: pairwise halo). Neither
        ! source is written below, so this is the loop-invariant hoist described in s_amr_subcycle_setup_block.
        if (amr_xchg_coarse_ghosts) then
            call s_amr_exchange_coarse_cons_halo(q_old)
            call s_amr_exchange_coarse_cons_halo(q_new)
        end if

        ! SETUP: each level-1 block prepares its two time-lerp ghost sources and zeros its registers (collective; ALL ranks call)
        do islot = 1, amr_num_blocks
            if (amr_block_level(islot) /= 1) cycle
            call s_amr_select_slot(islot)
            call s_amr_subcycle_setup_block(q_old, q_new, pb_old, mv_old, pb_in, mv_in)
        end do

        do sub = 1, 2
            do s = 1, 3
                th = (real(sub - 1, wp) + c_abs(s))*0.5_wp
                ! lerp every block's ghost shell to the stage time (+ substep-entry backup) BEFORE the seam halo reads interiors
                do islot = 1, amr_num_blocks
                    if (amr_block_level(islot) /= 1) cycle
                    call s_amr_select_slot(islot)
                    if (.not. amr_rank_owns_block) cycle
                    call s_amr_subtree_stage_lerp(s, th)
                end do
                ! reconcile shared seam ghosts among ADJACENT same-level blocks so both sides compute a matching flux. Tiling can
                ! split a wide feature into adjacent sub-blocks at ANY rank count (amr_maxc_fit caps a box at half the global extent
                ! even at np=1), so the halo runs unconditionally - it self-no-ops when there are no seam pairs, keeping every
                ! untiled
                ! case byte-identical.
                call s_amr_fine_fine_halo(0)
                ! RHS + RK update every block from the reconciled ghost shell
                do islot = 1, amr_num_blocks
                    if (amr_block_level(islot) /= 1) cycle
                    call s_amr_select_slot(islot)
                    if (.not. amr_rank_owns_block) cycle
                    call s_amr_subtree_stage_advance(amr_dt_fine, coefs, bc_type, q_T_sf, pb_in, rhs_pb, mv_in, rhs_mv, t_step, &
                                                     & s, th)
                end do
            end do
            ! after this substep EVERY level-1 block is at t_b (q_cons) with t_a in q_cons_stor: level 2 subcycles within [t_a, t_b]
            ! then folds back (restrict + Berger-Colella reflux). ONE level-wide call, not one per parent - the level-2 seam halo
            ! inside it spans all parents, so every owner must arrive at it together. No-op for single-level.
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
        if (qbmm .and. .not. polytropic) then
            ! the block's OWN side-state and rhs scratch (the coarse arrays stay untouched)
            call s_compute_rhs(amr_cons_br, q_T_sf, amr_scr_prim, bc_type, amr_scr_rhs, amr_slots(amr_cur)%pb_f%sf, amr_rhs_pb_f, &
                               & amr_slots(amr_cur)%mv_f%sf, amr_rhs_mv_f, t_step, s)
        else
            call s_compute_rhs(amr_cons_br, q_T_sf, amr_scr_prim, bc_type, amr_scr_rhs, pb_in, rhs_pb, mv_in, rhs_mv, t_step, s)
        end if
        call s_amr_br_store(amr_loc_of(amr_cur))
        call s_amr_restore_coarse()
        amr_in_fine_advance = .false.

        ! RK stage update at the FINE time step (device kernel)
        call s_amr_fine_rk_update(amr_loc_of(amr_cur), amr_scr_rhs, coefs(s, 1), coefs(s, 2), coefs(s, 3), coefs(s, 4), dt_sub)
        if (qbmm .and. .not. polytropic) then
            call s_amr_fine_rk_update_pbmv(amr_slots(amr_cur)%pb_f%sf, amr_slots(amr_cur)%mv_f%sf, amr_slots(amr_cur)%pb_stor%sf, &
                                           & amr_slots(amr_cur)%mv_stor%sf, amr_rhs_pb_f, amr_rhs_mv_f, coefs(s, 1), coefs(s, 2), &
                                           & coefs(s, 3), coefs(s, 4), dt_sub)
        end if
        ! 6-equation model: per-substage pressure relaxation (instantaneous equilibration - per stage at fine dt is the same
        ! infinite-rate limit the coarse applies per stage)
        if (model_eqns == model_eqns_6eq .and. (.not. relax)) call s_amr_pressure_relax_fine()
        ! moving body: rebuild the fine-block IB state at the body's fine sub-time position (th matches the fluid-ghost lerp)
        if (moving_immersed_boundary_flag) call s_amr_update_mib_fine(th)
        ! IB state correction on the fine block after each substep RK update (no-op unless ib)
        call s_amr_ib_correct_fine(amr_scr_prim)
        if (rank_time_wrt) call s_rank_time_toc()

    end subroutine s_amr_subtree_stage_advance

    !> Recursively subcycle EVERY block at level plev+1 - across ALL parents at once - within one of the parents' substeps [t_a,
    !! t_b] (duration dt_sub). Every level-plev block has just finished that substep: q_cons = parent @ t_b, q_cons_stor = parent @
    !! t_a. Per child: gather its two ghost-lerp sources from its OWN parent's two snapshots (parent-fine frame), recurse into level
    !! plev+2 at dt_sub/2 (a child takes amr_ref_ratio substeps covering [t_a, t_b]), then fold back into its parent - restrict the
    !! covered cells and apply the Berger-Colella C/F flux correction (s_amr_reflux_to_parent over dt_sub, consuming the child's
    !! freg + the parent-side creg captured during THIS substep). The registers already carry the matching per-substep time weights
    !! (freg 1/r*rk3_w, creg rk3_w), so conservation closes with no register changes.
    !!
    !! LEVEL-WIDE, NOT PER-PARENT. Driving one parent's whole subtree to completion before the next parent's put the interposed
    !! s_amr_fine_fine_halo(clev) out of lockstep: the halo exchanges EVERY level-clev seam pair, so a pair whose two blocks sit
    !! under different parents had only one side present. Co-location hid that (both ends of such a pair landed on one rank, making
    !! the exchange a local device copy), which is why multi-level subcycling was fail-closed at np>1. Walking the whole level
    !! together puts every owner at the same halo. At np=1 this only re-orders independent per-parent work - each child reads solely
    !! its own parent's finished snapshots and its same-level neighbours - so results are unchanged.
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
        ! SETUP each child: its two ghost-lerp sources from ITS OWN parent's substep endpoints (parent-fine frame) + zeroed
        ! registers
        do kc = 1, amr_num_blocks
            if (amr_block_level(kc) /= clev) cycle
            call s_amr_select_slot(kc)  ! amr_cur = kc; mirrors (isect already parent-fine)
            pblk = f_amr_parent_block(kc)
            if (.not. (amr_rank_owns_block .or. amr_block_owner(pblk) == proc_rank)) cycle
            ! Two P2P pairs (parent @ t_a, then @ t_b), so BOTH owners must arrive or the receiver never posts. The parent owner
            ! packs and sends from its own slot; the child owner receives WITHOUT naming the parent field - amr_slots(pblk) is
            ! unallocated there. Co-located (np=1, or parent and child on one rank) takes the local device-copy path unchanged.
            ! Both sends carry tag amr_cur; MPI non-overtaking on a fixed (source, tag, comm) keeps t_a ahead of t_b.
            if (amr_block_owner(pblk) == proc_rank) then
                call s_amr_gather_from_parent_field_stor(amr_cur, pblk, amr_loc_of(pblk), .false.)  ! parent @ t_a (device C/F fill)
                call s_amr_gather_send_flush()  ! keep this site's original blocking semantics - NO drain follows this loop
            else
                call s_amr_recv_parent_patch(pblk, .false.)
            end if
            if (amr_rank_owns_block) call s_amr_fill_fine_ghosts_gsta(amr_cg, amr_loc_of(kc))
            if (amr_block_owner(pblk) == proc_rank) then
                call s_amr_gather_from_parent_field_cons(amr_cur, pblk, amr_loc_of(pblk), .false.)  ! parent @ t_b (device C/F fill)
                call s_amr_gather_send_flush()  ! keep this site's original blocking semantics - NO drain follows this loop
            else
                call s_amr_recv_parent_patch(pblk, .false.)
            end if
            if (.not. amr_rank_owns_block) cycle
            call s_amr_fill_fine_ghosts_gstb(amr_cg, amr_loc_of(kc))
            call s_amr_zero_fine_registers()
        end do
        ! ADVANCE the level TRANSPOSED - every level-clev block through each substep together, with the level-clev seam halo
        ! interposed - exactly as s_amr_advance_fine_subcycle_all does at level 1. Advancing each child's whole subtree in turn (the
        ! previous shape) left adjacent blocks unable to see each other, so their shared face carried mismatched fluxes; that is why
        ! the regrid clamped subcycle to ONE capped child per box. The halo is level-filtered because this runs INSIDE one of the
        ! parents' substeps, when level plev is mid-substep and must not be touched.
        do sub = 1, 2
            do s = 1, 3
                th = (real(sub - 1, wp) + c_abs(s))*0.5_wp
                do kc = 1, amr_num_blocks
                    if (amr_block_level(kc) /= clev) cycle
                    call s_amr_select_slot(kc)
                    if (.not. amr_rank_owns_block) cycle
                    call s_amr_subtree_stage_lerp(s, th)
                end do
                call s_amr_fine_fine_halo(clev)
                do kc = 1, amr_num_blocks
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
        ! FOLD each child back into its parent (relax the fine phase first, matching the driver's relax -> restrict order)
        do kc = 1, amr_num_blocks
            if (amr_block_level(kc) /= clev) cycle
            call s_amr_select_slot(kc)
            ! The restrict and the reflux are each a P2P pair when child and parent are on different ranks, so BOTH participants
            ! must reach them or the receiver never posts and the pair deadlocks. Owner-only work (relax) stays behind the guard.
            if (amr_rank_owns_block) then
                if (relax) call s_amr_relax_fine()
            end if
            if (amr_rank_owns_block .or. amr_block_owner(f_amr_parent_block(kc)) == proc_rank) then
                call s_amr_restrict_to_parent()
                call s_amr_reflux_to_parent(dt_sub, .true.)
            end if
        end do

    end subroutine s_amr_advance_children

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

    !> Expand a candidate regrid box (global indices) to fully contain every immersed body it overlaps, with a buff_size margin (the
    !! IB image-point stencils need resolved surroundings). Expansion is re-clamped to the domain interior by the caller's own
    !! guards; a body too large for the per-rank block cap aborts with a named message. The bbox reads the live centroid, so a
    !! moving body's box tracks its current position; between regrids s_amr_update_mib_fine guards containment.

    !> Margin-padded global coarse-index bounding box of immersed body i (supported analytic geometries only; aborts on others).
    !! Reads the body's LIVE centroid, so a moving body's box tracks its current position.
    impure subroutine s_amr_body_bbox(i, mrg, blo, bhi)

        integer, intent(in)  :: i, mrg
        integer, intent(out) :: blo(3), bhi(3)
        real(wp)             :: c(3), half(3)

        c = [patch_ib(i)%x_centroid, patch_ib(i)%y_centroid, patch_ib(i)%z_centroid]
        select case (patch_ib(i)%geometry)
        case (2, 8, 10)  ! circle, sphere, cylinder: radius-bounded (cylinder length adds below)
            half = patch_ib(i)%radius
            if (patch_ib(i)%geometry == 10) then
                half(1) = max(half(1), 0.5_wp*patch_ib(i)%length_x)
                half(2) = max(half(2), 0.5_wp*patch_ib(i)%length_y)
                half(3) = max(half(3), 0.5_wp*patch_ib(i)%length_z)
            end if
        case (3, 9)  ! rectangle, box
            half = 0.5_wp*[patch_ib(i)%length_x, patch_ib(i)%length_y, patch_ib(i)%length_z]
        case default
            call s_mpi_abort('amr dynamic regrid with ib: unsupported body geometry for the ' &
                             & // 'containment bounding box (supported: circle/rectangle/sphere/box/cylinder)')
        end select
        ! physical bbox -> global coarse indices: uniform spacing only (stretched grids with ib-dynamic-regrid/Lagrangian are
        ! aborted
        ! at init; the axisymmetric half axis cell only shrinks dy(0), so the floor is still conservative)
        blo(1) = int((c(1) - half(1) - glb_bounds(1)%beg)/dx(0)) - mrg
        bhi(1) = int((c(1) + half(1) - glb_bounds(1)%beg)/dx(0)) + mrg
        blo(2) = 0; bhi(2) = 0; blo(3) = 0; bhi(3) = 0
        if (n_glb > 0) then
            blo(2) = int((c(2) - half(2) - glb_bounds(2)%beg)/dy(min(1, n))) - mrg
            bhi(2) = int((c(2) + half(2) - glb_bounds(2)%beg)/dy(min(1, n))) + mrg
        end if
        if (p_glb > 0) then
            blo(3) = int((c(3) - half(3) - glb_bounds(3)%beg)/dz(0)) - mrg
            bhi(3) = int((c(3) + half(3) - glb_bounds(3)%beg)/dz(0)) + mrg
        end if

    end subroutine s_amr_body_bbox

    impure subroutine s_amr_expand_box_over_bodies(lo, hi)

        integer, intent(inout) :: lo(3), hi(3)
        integer                :: i, d, blo(3), bhi(3), mrg
        logical                :: ovl

        ! containment margin: the IB image-point stencil reaches a few cells beyond the surface (the validated static-block goldens
        ! keep ~5); buff_size (floored to 10 by ib) would exceed the per-rank block cap for ordinary bodies. For amr_max_level > 1
        ! the
        ! body must survive every child nesting inset (amr_cpat_mar per level down to amr_max_level), so the parent block clears the
        ! body by that many extra cells - keeping the finest C/F boundary a full image-point stencil off the surface (refining the
        ! surface, not the interior).

        mrg = max(amr_buf, 4) + max(0, amr_max_level - 1)*amr_cpat_mar

        do i = 1, num_ibs
            call s_amr_body_bbox(i, mrg, blo, bhi)
            ! blocks must stay buff_size inside the domain: a body whose margin-padded bbox does not fit cannot be contained - fail
            ! with a named message instead of a clipped body
            if (blo(1) < buff_size .or. bhi(1) > m_glb - buff_size .or. (n_glb > 0 .and. (blo(2) < buff_size .or. bhi(2) > n_glb &
                & - buff_size)) .or. (p_glb > 0 .and. (blo(3) < buff_size .or. bhi(3) > p_glb - buff_size))) then
                call s_mpi_abort('amr dynamic regrid with ib: the immersed body plus its containment ' &
                                 & // 'margin does not fit inside the refinable domain interior (blocks stay buff_size off the edges)')
            end if
            ovl = lo(1) <= bhi(1) .and. hi(1) >= blo(1)
            if (n_glb > 0) ovl = ovl .and. lo(2) <= bhi(2) .and. hi(2) >= blo(2)
            if (p_glb > 0) ovl = ovl .and. lo(3) <= bhi(3) .and. hi(3) >= blo(3)
            if (.not. ovl) cycle
            do d = 1, num_dims
                lo(d) = min(lo(d), blo(d))
                hi(d) = max(hi(d), bhi(d))
            end do
            if (hi(1) - lo(1) + 1 > amr_maxc_fit(1) .or. (n_glb > 0 .and. hi(2) - lo(2) + 1 > amr_maxc_fit(2)) .or. (p_glb > 0 &
                & .and. hi(3) - lo(3) + 1 > amr_maxc_fit(3))) then
                call s_mpi_abort('amr dynamic regrid with ib: containing the immersed body plus margin ' &
                                 & // 'exceeds the per-rank block size cap; use fewer ranks or a larger amr_maxc_fit')
            end if
        end do

    end subroutine s_amr_expand_box_over_bodies

    !> max_grid_size tiling: split box [lo:hi] into a grid of contiguous sub-boxes each <= amr_maxc_fit per dim (the max a rank can
    !! whole-own), appending them to out(nt+1:). Tiles are ADJACENT (share fine seams) - the block-to-block fine-fine halo makes
    !! those seams conservative. Even split: ntl = ceil(ext/amr_maxc_fit) tiles, each of size ceil(ext/ntl) <= amr_maxc_fit. Sets
    !! capped=1 and stops if the amr_max_blocks cap is hit. Collapsed dims stay [0:0].
    pure subroutine s_amr_tile_box(lo, hi, out, nt, cap, capped, tsz)

        integer, intent(in)        :: lo(3), hi(3), cap
        type(t_box), intent(inout) :: out(:)
        integer, intent(inout)     :: nt, capped
        !> per-dim tile size (default amr_maxc_fit; a level-lev caller passes amr_maxc_fit/amr_ref_ratio**(lev-1) - the slot holds
        !! amr_ref_ratio*amr_maxc_fit fine cells and a level-lev block spans amr_ref_ratio**lev per coarse cell)
        integer, intent(in), optional :: tsz(3)
        integer                       :: ntl(3), s(3), t1, t2, t3, qlo(3), qhi(3), tc(3)

        tc = amr_maxc_fit; if (present(tsz)) tc = tsz
        tc = max(tc, 1)  ! a level>=2 caller passes amr_maxc_fit/2, which is 0 when a rank's fine half-extent is 1 (small subdomain
        !                  at high np) - a 0 tile size would divide-by-zero below; a 1-cell tile is the valid floor
        ntl = 1; s = 1
        ntl(1) = (hi(1) - lo(1) + tc(1))/tc(1); s(1) = (hi(1) - lo(1) + ntl(1))/ntl(1)
        if (n_glb > 0) then
            ntl(2) = (hi(2) - lo(2) + tc(2))/tc(2); s(2) = (hi(2) - lo(2) + ntl(2))/ntl(2)
        end if
        if (p_glb > 0) then
            ntl(3) = (hi(3) - lo(3) + tc(3))/tc(3); s(3) = (hi(3) - lo(3) + ntl(3))/ntl(3)
        end if
        do t3 = 0, ntl(3) - 1
            qlo(3) = 0; qhi(3) = 0
            if (p_glb > 0) then; qlo(3) = lo(3) + t3*s(3); qhi(3) = min(lo(3) + (t3 + 1)*s(3) - 1, hi(3)); end if
            do t2 = 0, ntl(2) - 1
                qlo(2) = 0; qhi(2) = 0
                if (n_glb > 0) then; qlo(2) = lo(2) + t2*s(2); qhi(2) = min(lo(2) + (t2 + 1)*s(2) - 1, hi(2)); end if
                do t1 = 0, ntl(1) - 1
                    if (nt >= cap) then; capped = 1; return; end if
                    qlo(1) = lo(1) + t1*s(1); qhi(1) = min(lo(1) + (t1 + 1)*s(1) - 1, hi(1))
                    nt = nt + 1; out(nt)%lo = qlo; out(nt)%hi = qhi
                end do
            end do
        end do

    end subroutine s_amr_tile_box

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

    !> Allocate slot islot's per-block field arrays (coords + the 6 device-resident field vectors + non-poly QBMM side-state), sized
    !! to the max buffered block. Idempotent (no-op if already live). The single QBMM RHS scratch amr_rhs_pb_f/mv_f and the global
    !! amr_cg are NOT per-slot and stay in init/finalize.
    !> Allocate/reset the dense local-index maps. Called from BOTH pool-allocation sites - s_initialize_amr_module and
    !! s_l0_tiles_init - because pure-L0 mode (amr = F) returns early from the former yet still calls s_amr_alloc_slot. Idempotent
    !! so either order is safe.
    impure subroutine s_amr_loc_index_init()

        if (.not. allocated(amr_loc_of)) allocate (amr_loc_of(1:amr_max_blocks))
        if (.not. allocated(amr_loc_free)) allocate (amr_loc_free(1:amr_max_blocks))
        amr_loc_of = 0; amr_loc_free = 0; amr_loc_n = 0; amr_loc_nfree = 0

    end subroutine s_amr_loc_index_init

    !> Move one local slot's store data src -> dst on the DEVICE, in place within each live store array (no staging copy: holding a
    !! second store-sized array on device is exactly the transient that OOMs at the S0 operating point). s_amr_compact_store's
    !! ascending-source ordering guarantees dst's previous contents are already consumed. The HOST copy goes stale, which is the
    !! store's normal state between rebuilds (device-authoritative; host readers pull per slot).
    impure subroutine s_amr_st_move_slot(src, dst)

        integer, intent(in) :: src, dst
        integer             :: i, j, k, l
        logical             :: want(4)

        want = [.true., .true., amr_subcycle, amr_subcycle]
        #:for ST, IDX in [('amr_cons_st', 1), ('amr_stor_st', 2), ('amr_gst_a', 3), ('amr_gst_b', 4)]
            if (want(${IDX}$)) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = 1, sys_size
                    do l = mbuf3_lo, mbuf3_hi
                        do k = mbuf2_lo, mbuf2_hi
                            do j = mbuf1_lo, mbuf1_hi
                                ${ST}$(j, k, l, i, dst) = ${ST}$(j, k, l, i, src)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

    end subroutine s_amr_st_move_slot

    !> Re-densify the flat store's local INDEX SPACE after every reconcile: the W8 weak-scaling fix. Slot indices are GLOBAL, so at
    !! np>=2 a rank's owned set is a shifting SFC window plus received migration slots; without re-densification `amr_loc_n`
    !! ratchets run-long (~210 MiB/slot of capacity) until the device OOMs at FIXED per-rank work (measured: S0 np=4). Renumbering
    !! every reconcile pins loc_n to the live count, so the capacity high-water plateaus at the rebuild transient (old + new block
    !! generations coexist mid-rebuild, ~2x live) instead of growing without bound. The allocation itself is NOT shrunk - the moves
    !! are device-side and in place, processed in ascending old-index order, which makes in-place safe: each destination (the rank
    !! of its source among the live indices) is <= its source, and every pending source lies above the current destination.
    !!
    !! Called at the END of s_amr_reconcile_slots, where every reader is finished: the rebuild's overlap carry-forward
    !! (which reads amr_stor_st at the OLD blocks' local indices) has completed, and the next rebuild has not started.
    impure subroutine s_amr_compact_store()

        integer              :: k, v, newloc
        integer, allocatable :: inv(:)

        if (amr_st_cap <= 0) return
        ! invert the map: inv(local index) = global slot, 0 if free. The walk covers the L0 tile prefix too - those slots
        ! hold store data even though the reconcile walk skips them.
        allocate (inv(amr_st_cap)); inv = 0
        do k = 1, amr_max_blocks
            if (amr_loc_of(k) > 0) inv(amr_loc_of(k)) = k
        end do
        newloc = 0
        do v = 1, amr_st_cap
            if (inv(v) == 0) cycle
            newloc = newloc + 1
            amr_loc_of(inv(v)) = newloc
            if (newloc /= v) call s_amr_st_move_slot(v, newloc)
        end do
        deallocate (inv)
#ifdef MFC_DEBUG
        if (newloc < amr_loc_n) write (0, '(A,I0,A,I0,A,I0,A,I0)') '[amr-compact] rank ', proc_rank, ' loc_n ', amr_loc_n, &
            & ' -> ', newloc, '  cap ', amr_st_cap
#endif
        amr_loc_n = newloc
        amr_loc_nfree = 0  ! every recycled index is invalid after renumbering
        amr_st_hw = newloc  ! let the trip-wire report the post-compaction trajectory

    end subroutine s_amr_compact_store

    !> Size the flat store for at least nloc local slots, growing 1.25x and never shrinking, so a run pays few reallocations no
    !! matter how the block count churns (with the index space re-densified every reconcile, growth fires only on a new
    !! rebuild-transient high-water). Growth must PRESERVE the live blocks' fields, and those live on the device, so each array
    !! stages through a DEVICE-side temporary (two on-device copies, no PCIe round trip) - the device-native remake the
    !! store-vs-AMReX analysis called for. The old device->host->device trip existed only to carry the migration stash's host writes
    !! across a growth; the stash chain is device-side now, so no reader depends on the host mirror through a grow.
    impure subroutine s_amr_st_reserve(nloc)

        integer, intent(in) :: nloc
        integer             :: oldcap, newcap, i
        integer             :: c5, i4, k3, j2, i1
        !> device-native staging doubles the array's device footprint transiently; above this many old columns, fall back to the
        !! host round trip (device peak max(old, new)). 32 covers the startup/early-regrid growth events whose round trips dominated
        !! short runs, while near-limit late growth keeps the OOM-safe path.
        integer, parameter     :: amr_grow_dev_cap = 32
        logical                :: want(4)
        real(stp), allocatable :: tmp(:,:,:,:,:), hstage(:,:,:,:,:)

        ! CONTRACT: the store is DEVICE-authoritative at every call; growth preserves the DEVICE contents only, and the host
        ! mirror comes out of a growth UNDEFINED (host readers pull per slot before reading - the store's normal state between
        ! rebuilds anyway, cf. s_amr_compact_store). A caller that has just written the store on the HOST (restart read) must
        ! still push its slot to the device before the next s_amr_alloc_slot, or that host data is lost.

        ! TRIP-WIRE on the store trajectory (kept from the W8 ratchet investigation). STDERR because stdout is buffered and
        ! lost on abort - the OOM run's stdout showed nothing at all.

        if (nloc > amr_st_hw) then
            amr_st_hw = nloc
#ifdef MFC_DEBUG
            write (0, '(A,I0,A,I0,A,I0,A,I0)') '[amr-store] rank ', proc_rank, ' NEW high-water nloc ', nloc, ' cap ', &
                   & amr_st_cap, ' recycle-depth ', amr_loc_nfree
#endif
        end if
        if (nloc <= amr_st_cap) return
        oldcap = amr_st_cap
        ! grow 1.25x with the increment CAPPED at 16 slots: a proportional increment is itself store-scaled, and at a large cap
        ! the +25% transient is what tips a near-limit device over (measured: S0 np=4 plateaued flat at 57.3 GiB, then one
        ! late-run growth event's +67-slot overshoot OOMed it). The +8 floor keeps early growth cheap when oldcap is tiny.
        newcap = max(oldcap + max(min(oldcap/4, 16), 8), nloc)
        want = [.true., .true., amr_subcycle, amr_subcycle]

        #:for ST, IDX in [('amr_cons_st', 1), ('amr_stor_st', 2), ('amr_gst_a', 3), ('amr_gst_b', 4)]
            if (want(${IDX}$)) then
                if (oldcap > amr_grow_dev_cap) then
                    ! near-limit fallback: the device-native staging below transiently holds old + tmp = 2*oldcap columns
                    ! on the device, and growth fires exactly at the memory high-water mark - a measured OOM class (a
                    ! +25%-increment transient alone tipped a 57.3 GiB np=4 run; see the cap-sweep history). Above the
                    ! threshold, keep the host round trip: slow (full PCIe both ways) but its device peak is max(old, new).
                    $:GPU_UPDATE(host='[' + ST + ']')
                    allocate (hstage(mbuf1_lo:mbuf1_hi,mbuf2_lo:mbuf2_hi,mbuf3_lo:mbuf3_hi,1:sys_size,1:oldcap))
                    hstage = ${ST}$(:,:,:,:,1:oldcap)
                    @:DEALLOCATE(${ST}$)
                    @:ALLOCATE(${ST}$(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi, 1:sys_size, 1:newcap))
                    ${ST}$(:,:,:,:,1:oldcap) = hstage
                    ${ST}$(:,:,:,:,oldcap + 1:newcap) = 0._stp
                    deallocate (hstage)
                    $:GPU_UPDATE(device='[' + ST + ']')
                else
                    if (oldcap > 0) then
                        ! stage the live columns on the DEVICE (tmp is device-mapped by @:ALLOCATE); no PCIe traffic
                        @:ALLOCATE(tmp(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi, 1:sys_size, 1:oldcap))
                        $:GPU_PARALLEL_LOOP(collapse=4)
                        do c5 = 1, oldcap
                            do i4 = 1, sys_size
                                do k3 = mbuf3_lo, mbuf3_hi
                                    do j2 = mbuf2_lo, mbuf2_hi
                                        do i1 = mbuf1_lo, mbuf1_hi
                                            tmp(i1, j2, k3, i4, c5) = ${ST}$(i1, j2, k3, i4, c5)
                                        end do
                                    end do
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()
                        @:DEALLOCATE(${ST}$)
                    end if
                    @:ALLOCATE(${ST}$(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi, 1:sys_size, 1:newcap))
                    ! restore the preserved columns and zero the rest, both on the device; the host mirror stays undefined
                    ! (see the contract above - every host reader pulls its slot first). Two kernels so the zero-only path
                    ! (oldcap == 0) never references the unallocated tmp.
                    if (oldcap > 0) then
                        $:GPU_PARALLEL_LOOP(collapse=4)
                        do c5 = 1, oldcap
                            do i4 = 1, sys_size
                                do k3 = mbuf3_lo, mbuf3_hi
                                    do j2 = mbuf2_lo, mbuf2_hi
                                        do i1 = mbuf1_lo, mbuf1_hi
                                            ${ST}$(i1, j2, k3, i4, c5) = tmp(i1, j2, k3, i4, c5)
                                        end do
                                    end do
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()
                        @:DEALLOCATE(tmp)
                    end if
                    $:GPU_PARALLEL_LOOP(collapse=4)
                    do c5 = oldcap + 1, newcap
                        do i4 = 1, sys_size
                            do k3 = mbuf3_lo, mbuf3_hi
                                do j2 = mbuf2_lo, mbuf2_hi
                                    do i1 = mbuf1_lo, mbuf1_hi
                                        ${ST}$(i1, j2, k3, i4, c5) = 0._stp
                                    end do
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
            end if
        #:endfor

        amr_st_cap = newcap

        ! the bridge spans a BOUNDED batch of blocks along k (see amr_br_batch), so one s_compute_rhs call can advance a
        ! whole batch instead of one block; it rides the same pool lifetime and is rebuilt only if the window changes
        amr_br_w = mbuf3_hi - mbuf3_lo + 1
        if (.not. allocated(amr_cons_br)) then
            @:ALLOCATE(amr_cons_br(1:sys_size))
            do i = 1, sys_size
                @:ALLOCATE(amr_cons_br(i)%sf(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_lo + amr_br_batch*amr_br_w - 1))
                @:ACC_SETUP_SFs(amr_cons_br(i))
            end do
            amr_br_nblk = amr_br_batch
        end if

        ! 2a prim landing zone + batch metadata: per-STAGE scratch (rewritten by every s_amr_convert_prim_batch
        ! call), so growth discards contents - no device staging round trip, unlike the stores above.
        if (amr_prim_batch) then
            if (allocated(amr_prim_st)) then
                @:DEALLOCATE(amr_prim_st)
                @:DEALLOCATE(amr_bt_lo)
                @:DEALLOCATE(amr_bt_hi)
                @:DEALLOCATE(amr_bt_on)
            end if
            allocate (amr_prim_st(mbuf1_lo:mbuf1_hi,mbuf2_lo:mbuf2_hi,mbuf3_lo:mbuf3_hi,1:num_vels + 1,1:newcap))
            allocate (amr_bt_lo(3, newcap), amr_bt_hi(3, newcap), amr_bt_on(newcap))
        end if

    end subroutine s_amr_st_reserve

    !> 2a: ONE batched cons->prim conversion over every owned fine block (all levels), straight from the flat cons store into the
    !! flat prim landing zone. Runs once per RK stage after the fill + seam phases; each block's store bytes there are identical to
    !! what its per-block conversion point would read (advances write only their own slots), and the kernel is per-cell with no
    !! reductions, so the result is bit-identical to the per-block path it replaces. The cell body below is PINNED to
    !! s_convert_conservative_to_primitive_variables (m_variables_conversion.fpp) restricted to the amr_prim_batch gate's configs:
    !! species fractions (s_compute_species_fraction inlined against the store; igr/bubbles_euler excluded by the gate), mixture
    !! properties, velocity + dynamic pressure, and pressure. Change the conversion and this must follow.
    impure subroutine s_amr_convert_prim_batch()

        integer :: g, loc, i, j, k, l
        integer :: nl, nv, b1l, b1h, b2l, b2h, b3l, b3h

        #:if USING_AMD and not MFC_CASE_OPTIMIZATION
            real(wp), dimension(3) :: alpha_K, alpha_rho_K
            real(wp)               :: rhoYks_b(1:10)
        #:else
            real(wp), dimension(num_fluids) :: alpha_K, alpha_rho_K
            real(wp)                        :: rhoYks_b(1:num_species)
        #:endif
        real(wp) :: Re_K(2)
        real(wp) :: rho_K, gamma_K, pi_inf_K, qv_K, dyn_pres_K, alpha_K_sum, pres, T, pmag

        if (amr_loc_n == 0) return
        call s_phase_tic(PH_CVTB)
        amr_bt_on(1:amr_loc_n) = .false.
        do g = 1, amr_num_blocks
            if (amr_block_level(g) < 1) cycle
            if (amr_block_owner(g) /= proc_rank) cycle
            loc = amr_loc_of(g)
            if (loc <= 0) cycle
            amr_bt_on(loc) = .true.
            do i = 1, 3
                amr_bt_lo(i, loc) = amr_slots(g)%idwbuff(i)%beg
                amr_bt_hi(i, loc) = amr_slots(g)%idwbuff(i)%end
            end do
        end do
        $:GPU_UPDATE(device='[amr_bt_on, amr_bt_lo, amr_bt_hi]')
        ! bounds through local scalars, never GPU_DECLARE'd module state (the CCE-acc stale-device-bounds class)
        nl = amr_loc_n; nv = num_vels
        b1l = mbuf1_lo; b1h = mbuf1_hi; b2l = mbuf2_lo; b2h = mbuf2_hi; b3l = mbuf3_lo; b3h = mbuf3_hi
        $:GPU_PARALLEL_LOOP(collapse=4, private='[alpha_K, alpha_rho_K, Re_K, rhoYks_b, rho_K, gamma_K, pi_inf_K, qv_K, &
                            & dyn_pres_K, alpha_K_sum, pres, T, pmag]', copyin='[nl, nv, b1l, b1h, b2l, b2h, b3l, b3h]')
        do loc = 1, nl
            do l = b3l, b3h
                do k = b2l, b2h
                    do j = b1l, b1h
                        if (.not. amr_bt_on(loc)) cycle
                        if (j < amr_bt_lo(1, loc) .or. j > amr_bt_hi(1, loc) .or. k < amr_bt_lo(2, loc) .or. k > amr_bt_hi(2, &
                            & loc) .or. l < amr_bt_lo(3, loc) .or. l > amr_bt_hi(3, loc)) cycle
                        if (num_fluids == 1) then
                            alpha_rho_K(1) = amr_cons_st(j, k, l, eqn_idx%cont%beg, loc)
                            alpha_K(1) = amr_cons_st(j, k, l, eqn_idx%adv%beg, loc)
                        else
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                alpha_rho_K(i) = amr_cons_st(j, k, l, i, loc)
                                alpha_K(i) = amr_cons_st(j, k, l, eqn_idx%adv%beg + i - 1, loc)
                            end do
                        end if
                        if (mpp_lim) then
                            alpha_K_sum = 0._wp
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                alpha_rho_K(i) = max(0._wp, alpha_rho_K(i))
                                alpha_K(i) = min(max(0._wp, alpha_K(i)), 1._wp)
                                alpha_K_sum = alpha_K_sum + alpha_K(i)
                            end do
                            ! explicit loop, not array syntax: an inline whole-array expression in a target
                            ! region is a per-thread temporary on amdflang
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                alpha_K(i) = alpha_K(i)/max(alpha_K_sum, 1.e-16_wp)
                            end do
                        end if
                        call s_convert_species_to_mixture_variables_kernel(rho_K, gamma_K, pi_inf_K, qv_K, alpha_K, alpha_rho_K, &
                            & Re_K)
                        if (enforce_density_floor_vc) rho_K = max(rho_K, sgm_eps)
                        dyn_pres_K = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, nv
                            amr_prim_st(j, k, l, i, loc) = amr_cons_st(j, k, l, eqn_idx%mom%beg + i - 1, loc)/rho_K
                            dyn_pres_K = dyn_pres_K + 5.e-1_wp*amr_cons_st(j, k, l, eqn_idx%mom%beg + i - 1, loc)*amr_prim_st(j, &
                                & k, l, i, loc)
                        end do
                        pmag = 0._wp
                        call s_compute_pressure(amr_cons_st(j, k, l, eqn_idx%E, loc), amr_cons_st(j, k, l, eqn_idx%alf, loc), &
                                                & dyn_pres_K, pi_inf_K, gamma_K, rho_K, qv_K, rhoYks_b, pres, T, pres_mag=pmag)
                        amr_prim_st(j, k, l, nv + 1, loc) = pres
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        call s_phase_toc(PH_CVTB)

    end subroutine s_amr_convert_prim_batch

    !> 2a: land the current block's batch-computed prim vars (the contiguous mom%beg..E range) from the prim store into the m_rhs
    !! conversion scratch, replacing that block's per-block conversion.
    impure subroutine s_amr_prim_load(q_prim_b, loc)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_b
        integer, intent(in)                                    :: loc
        integer                                                :: i, lb(3)

        ! One launch per var through a PLAIN contiguous array dummy: a scalar_field-array dummy makes the
        ! target region map the derived-type descriptors per launch (the measured per-region attach cost this
        ! per-block path must not pay). The dummy rebases to 1, so offset from the actual's own lower bounds.

        do i = eqn_idx%mom%beg, eqn_idx%E
            lb = lbound(q_prim_b(i)%sf)
            call s_amr_prim_load_one(q_prim_b(i)%sf, i - eqn_idx%mom%beg + 1, loc, amr_slots(amr_cur)%idwbuff(1)%beg, &
                                     & amr_slots(amr_cur)%idwbuff(1)%end, amr_slots(amr_cur)%idwbuff(2)%beg, &
                                     & amr_slots(amr_cur)%idwbuff(2)%end, amr_slots(amr_cur)%idwbuff(3)%beg, &
                                     & amr_slots(amr_cur)%idwbuff(3)%end, 1 - lb(1), 1 - lb(2), 1 - lb(3))
        end do

    end subroutine s_amr_prim_load

    impure subroutine s_amr_prim_load_one(dst, pv, loc, j1l, j1h, j2l, j2h, j3l, j3h, o1, o2, o3)

        real(stp), dimension(:,:,:), contiguous, intent(inout) :: dst
        integer, intent(in)                                    :: pv, loc, j1l, j1h, j2l, j2h, j3l, j3h, o1, o2, o3
        integer                                                :: j, k, l

        $:GPU_PARALLEL_LOOP(collapse=3, copyin='[pv, loc, j1l, j1h, j2l, j2h, j3l, j3h, o1, o2, o3]')
        do l = j3l, j3h
            do k = j2l, j2h
                do j = j1l, j1h
                    dst(j + o1, k + o2, l + o3) = amr_prim_st(j, k, l, pv, loc)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_prim_load_one

    !> Allocate the pooled q_prim/rhs advance scratch (idempotent). The lockstep driver argument-associates the scratch for every
    !! block, owned or not, so it must exist on EVERY rank - including one that never allocates a slot. Called from the ONE point
    !! per mode where mbuf* are final: end of s_initialize_amr_module when l0_ntile == 0 (pure AMR), and after s_l0_tiles_init's
    !! mbuf UNION when tiles exist (pure-L0 AND coexist - the union can enlarge mbuf* past the fine-only values, so an earlier
    !! allocation would undersize the scratch). rhs mirrors the per-slot igr widening.
    impure subroutine s_amr_scr_init()

        integer :: i

        if (allocated(amr_scr_prim)) return
        @:ALLOCATE(amr_scr_prim(1:sys_size))
        @:ALLOCATE(amr_scr_rhs(1:sys_size))
        do i = 1, sys_size
            @:ALLOCATE(amr_scr_prim(i)%sf(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi))
            if (igr) then
                @:ALLOCATE(amr_scr_rhs(i)%sf(mbuf1_lo:mbuf1_hi, min(mbuf2_lo, -1):max(mbuf2_hi, 1), min(mbuf3_lo, &
                           & -1):max(mbuf3_hi, 1)))
            else
                @:ALLOCATE(amr_scr_rhs(i)%sf(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi))
            end if
            @:ACC_SETUP_SFs(amr_scr_prim(i))
            @:ACC_SETUP_SFs(amr_scr_rhs(i))
        end do

    end subroutine s_amr_scr_init

    !> Move block loc's conserved state between the flat store and the bridge. One kernel each way over the whole buffered box:
    !! every slot's arrays carry the same mbuf extents, so this is exactly the box the per-slot q_cons used to own.
    #:set BR = 'amr_cons_br(i)%sf(j, k, l)'
    #:set ST = 'amr_cons_st(j, k, l, i, loc)'
    #:for DIR in ['load', 'store']
        #:set LHS = BR if DIR == 'load' else ST
        #:set RHS = ST if DIR == 'load' else BR
        impure subroutine s_amr_br_${DIR}$(loc)

            integer, intent(in) :: loc
            integer             :: i, j, k, l

            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do l = mbuf3_lo, mbuf3_hi
                    do k = mbuf2_lo, mbuf2_hi
                        do j = mbuf1_lo, mbuf1_hi
                            ${LHS}$ = ${RHS}$
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

        end subroutine s_amr_br_${DIR}$
    #:endfor

    !> BATCHED bridge move: a whole batch of blocks in ONE kernel instead of one kernel per block. This is what the flat store was
    !! built for - its own comment calls a single contiguous array "the prerequisite for running ONE kernel over every live block
    !! instead of one kernel per block". Block `loc` of the batch lands at k-offset (loc-1)*amr_br_w with its ghost shell intact, so
    !! the bridge holds a concatenation of independent blocks two ghost shells apart. The per-block s_amr_br_load/store above are
    !! KEPT unchanged: reflux and the relaxation paths still work one block at a time and read the bridge at offset 0.
    #:set BRA = 'amr_cons_br(i)%sf(j, k, l + (loc - 1)*amr_br_w)'
    #:set STA = 'amr_cons_st(j, k, l, i, base + loc)'
    #:for DIR in ['load', 'store']
        #:set LHS = BRA if DIR == 'load' else STA
        #:set RHS = STA if DIR == 'load' else BRA
        impure subroutine s_amr_br_${DIR}$_all(base, nblk)

            integer, intent(in) :: base  !< dense local slot preceding this batch
            integer, intent(in) :: nblk  !< blocks in this batch, <= amr_br_batch
            integer             :: i, j, k, l, loc

            $:GPU_PARALLEL_LOOP(collapse=5)
            do loc = 1, nblk
                do i = 1, sys_size
                    do l = mbuf3_lo, mbuf3_hi
                        do k = mbuf2_lo, mbuf2_hi
                            do j = mbuf1_lo, mbuf1_hi
                                ${LHS}$ = ${RHS}$
                            end do
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

        end subroutine s_amr_br_${DIR}$_all
    #:endfor

    !> Free the flat store and the dense-index maps. Mirrors s_amr_loc_index_init: called from BOTH finalize paths, because either
    !! pool-allocation site can have created them. Idempotent.
    impure subroutine s_amr_st_finalize()

        integer :: i

        if (allocated(amr_loc_of)) deallocate (amr_loc_of, amr_loc_free)
        #:for ST in ['amr_cons_st', 'amr_stor_st', 'amr_gst_a', 'amr_gst_b']
            if (allocated(${ST}$)) then
                @:DEALLOCATE(${ST}$)
            end if
        #:endfor
        if (allocated(amr_prim_st)) then
            @:DEALLOCATE(amr_prim_st)
            @:DEALLOCATE(amr_bt_lo)
            @:DEALLOCATE(amr_bt_hi)
            @:DEALLOCATE(amr_bt_on)
        end if
        if (allocated(amr_cons_br)) then
            do i = 1, sys_size
                @:ACC_TEARDOWN_SFs(amr_cons_br(i))
                @:DEALLOCATE(amr_cons_br(i)%sf)
            end do
            @:DEALLOCATE(amr_cons_br)
        end if
        if (allocated(amr_scr_prim)) then
            do i = 1, sys_size
                @:ACC_TEARDOWN_SFs(amr_scr_prim(i))
                @:DEALLOCATE(amr_scr_prim(i)%sf)
                @:ACC_TEARDOWN_SFs(amr_scr_rhs(i))
                @:DEALLOCATE(amr_scr_rhs(i)%sf)
            end do
            @:DEALLOCATE(amr_scr_prim)
            @:DEALLOCATE(amr_scr_rhs)
        end if
        amr_st_cap = 0

    end subroutine s_amr_st_finalize

    !> One-shot pre-reserve before a batch of s_amr_alloc_slot_stash calls: grows the store AT MOST ONCE, to the batch's exact final
    !! size, instead of once per 8-16 incrementally allocated slots. Every growth restages the WHOLE store (both arrays) through the
    !! host (s_amr_st_reserve), so a migration wave allocating tens of replica slots pays that restage repeatedly without this.
    !! getk(k) marks dense fine-block index k (slot f_l0_slot(k)) for a coming stash alloc.
    impure subroutine s_amr_prereserve_stash(getk, nblk)

        integer, intent(in) :: nblk
        logical, intent(in) :: getk(nblk)
        integer             :: i, nneed

        nneed = 0
        do i = 1, nblk
            if (getk(i) .and. .not. amr_slot_live(f_l0_slot(i))) nneed = nneed + 1
        end do
        ! mirrors the alloc loop exactly: the first amr_loc_nfree allocs pop the recycle stack and do
        ! not raise amr_loc_n; only the remainder grow
        call s_amr_st_reserve(amr_loc_n + max(nneed - amr_loc_nfree, 0))

    end subroutine s_amr_prereserve_stash

    !> Assign slot islot a dense store index WITHOUT the per-block field arrays: a migration REPLICA only ever has its amr_stor_st
    !! slot written (receive-unpack) and read (the rebuild's overlap carry-forward), so q_prim/rhs would roughly double its device
    !! cost across the np-scaled replica set of a migration-heavy regrid (the W8 gate's np=4 arm OOMed on exactly that storm).
    !! s_amr_alloc_slot upgrades a stash-only slot in place when the same global slot becomes an owned block; s_amr_free_slot
    !! handles both flavors.
    impure subroutine s_amr_alloc_slot_stash(islot)

        integer, intent(in) :: islot

        if (amr_slot_live(islot)) return
        if (amr_loc_nfree > 0) then
            amr_loc_of(islot) = amr_loc_free(amr_loc_nfree)
            amr_loc_nfree = amr_loc_nfree - 1
        else
            amr_loc_n = amr_loc_n + 1
            amr_loc_of(islot) = amr_loc_n
        end if
        amr_slot_live(islot) = .true.
        call s_amr_st_reserve(amr_loc_n)

    end subroutine s_amr_alloc_slot_stash

    impure subroutine s_amr_alloc_slot(islot)

        integer, intent(in) :: islot
        integer             :: i

        ! full-vs-stash discriminator is the grid arrays: every full slot (tile or fine) has x_cb; a stash-only slot has none
        ! (fine slots no longer carry q_prim - see the pooled scratch amr_scr_prim/amr_scr_rhs)

        if (amr_slot_live(islot) .and. allocated(amr_slots(islot)%x_cb)) return
        ! recycle a freed local index if one is available, else extend the dense range; a live STASH-ONLY slot upgrading to a
        ! full one keeps the index (and so the stor data) it already holds
        if (.not. amr_slot_live(islot)) then
            if (amr_loc_nfree > 0) then
                amr_loc_of(islot) = amr_loc_free(amr_loc_nfree)
                amr_loc_nfree = amr_loc_nfree - 1
            else
                amr_loc_n = amr_loc_n + 1
                amr_loc_of(islot) = amr_loc_n
            end if
        end if
        amr_slots(islot)%amr_ref_ratio = amr_ref_ratio
        amr_slots(islot)%buff_size = buff_size
        allocate (amr_slots(islot)%x_cb(-1:max_f1), amr_slots(islot)%x_cc(0:max_f1), amr_slots(islot)%dx(0:max_f1))
        if (n_glb > 0) allocate (amr_slots(islot)%y_cb(-1:max_f2), amr_slots(islot)%y_cc(0:max_f2), amr_slots(islot)%dy(0:max_f2))
        if (p_glb > 0) allocate (amr_slots(islot)%z_cb(-1:max_f3), amr_slots(islot)%z_cc(0:max_f3), amr_slots(islot)%dz(0:max_f3))
        ! P1 pooling: fine blocks advance through the shared scratch (amr_scr_prim/amr_scr_rhs) - the fused per-block advance
        ! leaves no cross-block q_prim/rhs lifetime. L0 tile slots are the exception: all owned tiles' rhs coexist across the
        ! MPI-synchronized reflux point (s_l0_add_reflux_to_tiles between the whole-set RHS and RK passes), and a tile's q_prim
        ! written by the RHS pass is read in the later RK pass (IB correction), so tiles keep per-slot rhs always and per-slot
        ! q_prim exactly when s_compute_rhs's copy-out gate writes it (m_rhs.fpp end-of-rhs gate).
        if (islot <= l0_slot_off) then
            @:ALLOCATE(amr_slots(islot)%rhs(1:sys_size))
            if (run_time_info .or. probe_wrt .or. ib .or. bubbles_lagrange) then
                @:ALLOCATE(amr_slots(islot)%q_prim(1:sys_size))
            end if
            do i = 1, sys_size
                ! rhs is ghost-inclusive (mbuf); igr widens to -1:+1 per dim including collapsed ones (coarse rhs_vf is -1:m+1 etc.)
                if (igr) then
                    @:ALLOCATE(amr_slots(islot)%rhs(i)%sf(mbuf1_lo:mbuf1_hi, min(mbuf2_lo, -1):max(mbuf2_hi, 1), min(mbuf3_lo, &
                               & -1):max(mbuf3_hi, 1)))
                else
                    @:ALLOCATE(amr_slots(islot)%rhs(i)%sf(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi))
                end if
                @:ACC_SETUP_SFs(amr_slots(islot)%rhs(i))
                if (allocated(amr_slots(islot)%q_prim)) then
                    @:ALLOCATE(amr_slots(islot)%q_prim(i)%sf(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi))
                    @:ACC_SETUP_SFs(amr_slots(islot)%q_prim(i))
                end if
            end do
        end if
        if (qbmm .and. .not. polytropic) then
            #:for PF in ['pb_f', 'mv_f', 'pb_stor', 'mv_stor']
                @:ALLOCATE(amr_slots(islot)%${PF}$%sf(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi, 1:nnode, 1:nb))
                @:ACC_SETUP_SFs(amr_slots(islot)%${PF}$)
            #:endfor
            if (amr_subcycle) then
                #:for PF in ['pb_ghost_a', 'mv_ghost_a', 'pb_ghost_b', 'mv_ghost_b']
                    @:ALLOCATE(amr_slots(islot)%${PF}$%sf(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi, 1:nnode, 1:nb))
                    @:ACC_SETUP_SFs(amr_slots(islot)%${PF}$)
                #:endfor
            end if
        end if
        amr_slot_live(islot) = .true.
        call s_amr_st_reserve(amr_loc_n)

    end subroutine s_amr_alloc_slot

    !> Free slot islot's per-block field arrays (inverse of s_amr_alloc_slot). Idempotent (no-op if not live).
    impure subroutine s_amr_free_slot(islot)

        integer, intent(in) :: islot
        integer             :: i

        if (.not. amr_slot_live(islot)) return
        if (amr_loc_of(islot) > 0) then
            amr_loc_nfree = amr_loc_nfree + 1
            amr_loc_free(amr_loc_nfree) = amr_loc_of(islot)
            amr_loc_of(islot) = 0
        end if
        ! Undo each field's ACC_SETUP_SFs (Cray descriptor + %sf copyin) BEFORE the @:DEALLOCATE - Cray 'exit data delete'
        ! decrements
        ! the ref count, so the lone @:DEALLOCATE would leave the descriptor and the ACC_SETUP %sf ref dangling; the leaked host
        ! address is later reused (e.g. by Gs_rs at restart), tripping a Cray "Error placing / already present" present-table crash
        ! (gpu-acc). A STASH-ONLY slot (s_amr_alloc_slot_stash) has none of these arrays - only the index bookkeeping above.
        if (allocated(amr_slots(islot)%q_prim)) then
            do i = 1, sys_size
                @:ACC_TEARDOWN_SFs(amr_slots(islot)%q_prim(i))
                @:DEALLOCATE(amr_slots(islot)%q_prim(i)%sf)
            end do
            @:DEALLOCATE(amr_slots(islot)%q_prim)
        end if
        if (allocated(amr_slots(islot)%rhs)) then
            do i = 1, sys_size
                @:ACC_TEARDOWN_SFs(amr_slots(islot)%rhs(i))
                @:DEALLOCATE(amr_slots(islot)%rhs(i)%sf)
            end do
            @:DEALLOCATE(amr_slots(islot)%rhs)
        end if
        if (qbmm .and. .not. polytropic .and. associated(amr_slots(islot)%pb_f%sf)) then
            #:for PF in ['pb_f', 'mv_f', 'pb_stor', 'mv_stor']
                @:ACC_TEARDOWN_SFs(amr_slots(islot)%${PF}$)
                @:DEALLOCATE(amr_slots(islot)%${PF}$%sf)
            #:endfor
            if (amr_subcycle) then
                #:for PF in ['pb_ghost_a', 'mv_ghost_a', 'pb_ghost_b', 'mv_ghost_b']
                    @:ACC_TEARDOWN_SFs(amr_slots(islot)%${PF}$)
                    @:DEALLOCATE(amr_slots(islot)%${PF}$%sf)
                #:endfor
            end if
        end if
        if (allocated(amr_slots(islot)%x_cb)) deallocate (amr_slots(islot)%x_cb, amr_slots(islot)%x_cc, amr_slots(islot)%dx)
        if (allocated(amr_slots(islot)%y_cb)) deallocate (amr_slots(islot)%y_cb, amr_slots(islot)%y_cc, amr_slots(islot)%dy)
        if (allocated(amr_slots(islot)%z_cb)) deallocate (amr_slots(islot)%z_cb, amr_slots(islot)%z_cc, amr_slots(islot)%dz)
        amr_slot_live(islot) = .false.

    end subroutine s_amr_free_slot

    !> Reconcile the allocated per-slot field arrays to the CURRENT ownership: allocate every active block this rank owns, free
    !! everything else. Call after ownership is set (init/regrid/restart). A rank ends holding only its owned blocks' fine arrays
    !! (~amr_num_blocks/num_procs of the pool), not all amr_max_blocks. Regrid must alloc its transient (received/old) slots BEFORE
    !! this call, since it frees anything not currently owned.
    impure subroutine s_amr_reconcile_slots()

        integer :: k, nliv, nfr, nal, nfree_in
        logical :: needed

        nliv = 0; nfr = 0; nal = 0; nfree_in = amr_loc_nfree
        do k = 1, amr_max_blocks
            ! Skip the L0 tile prefix [1..l0_slot_off]: those slots are owned + sized by s_l0_tiles_init (rr=1 tile geometry), not
            ! by the AMR fine-block reconcile. Without this, at coexist init the tile-prefix owner defaults to 0, so RANK 0 would
            ! alloc these slots here with the FINE mbuf* sizing; s_amr_alloc_slot is idempotent, so s_l0_build_tile_slot could not
            ! then resize them - a latent tile-undersizing landmine (benign only while fine mbuf* >= tile). l0_slot_off = 0 with no
            ! tiles, so this is a no-op for pure AMR. (Audit IMP-3.)
            if (k <= l0_slot_off) cycle
            needed = k <= amr_num_blocks
            if (needed) needed = amr_block_owner(k) == proc_rank
            if (needed) then
                if (.not. amr_slot_live(k)) nal = nal + 1
                call s_amr_alloc_slot(k)
                nliv = nliv + 1
            else
                if (amr_slot_live(k)) nfr = nfr + 1
                call s_amr_free_slot(k)
            end if
        end do
        ! stderr (survives an abort). live = what the run actually needs; loc_n = indices ever handed out;
        ! the gap between them IS the leak. stack_in/out shows whether frees accumulate for the next
        ! rebuild to recycle, and 'new' counts allocs that had to EXTEND rather than recycle.
#ifdef MFC_DEBUG
        write (0, '(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') '[amr-recon] rank ', proc_rank, ' live ', nliv, ' loc_n ', amr_loc_n, &
               & ' freed ', nfr, ' newalloc ', nal, ' stack_in ', nfree_in, ' stack_out ', amr_loc_nfree
#endif
        ! I1a invariant: no stash-only slot survives a reconcile - every migration replica was freed (early-free or the walk
        ! above) or upgraded to a full slot by the owned-path alloc. A survivor would reach the solver with no geometry arrays.
        do k = 1, amr_max_blocks
            if (amr_slot_live(k)) then
                @:ASSERT(allocated(amr_slots(k)%x_cb), "a stash-only replica slot survived reconcile")
            end if
        end do
        ! every reader of a stale slot has finished by here (the rebuild's overlap carry-forward is done and the next rebuild
        ! has not started), which is what makes the renumbering safe - see s_amr_compact_store.
        call s_amr_compact_store()
        ! TRACK S: per-rank store capacity is the W8 invariant (device memory = f(live local boxes)); wall time cannot see it,
        ! so report it like [amr-scale]. live == loc_n after the compaction above; cap - live is the rebuild-transient envelope.
        if (rank_time_wrt) write (0, '(A,I0,A,I0,A,I0)') '[amr-cap] rank ', proc_rank, ' live ', amr_loc_n, ' cap ', amr_st_cap
        amr_mesh_epoch = amr_mesh_epoch + 1  ! local slot indices may have been renumbered: plans that baked them are stale

    end subroutine s_amr_reconcile_slots

    ! L0-AS-BLOCKS SPIKE (l0_ntile > 0):
    ! Base grid tiled into l0_ntiles_tot refinement-ratio-1 blocks advanced through the shared swap-based per-block solver; the
    ! bit-identity oracle is l0_ntile=0 (monolithic). See the l0_ntiles_tot declaration above for the design.

    !> Low global-cell index of tile `it` (0-based) when `ncell` cells split into `nt` balanced tiles: the first mod(ncell,nt) tiles
    !! get one extra cell. Tile `it` spans [f_l0_lo(it) : f_l0_lo(it+1)-1]; the widest tile is ceil(ncell/nt) cells.
    pure integer function f_l0_lo(ncell, nt, it)

        integer, intent(in) :: ncell, nt, it

        f_l0_lo = it*(ncell/nt) + min(it, mod(ncell, nt))

    end function f_l0_lo

    !> True if a domain-face BC code is a PHYSICAL boundary the spike does not support. Supported: extrapolation (BC_GHOST_EXTRAP),
    !! reflective (BC_REFLECTIVE), and periodic (BC_PERIODIC) - each has a cons-space tile fill matching the monolithic prim-space
    !! BC. The characteristic / slip-wall / dirichlet family (bc < BC_GHOST_EXTRAP) has no tile fill yet and is rejected. MPI
    !! processor boundaries (bc >= 0, incl. a periodic wrap-neighbour rank at np>1) are interior seams handled by the fine-fine
    !! halo, never a physical face here.
    pure logical function f_l0_bc_unsupported(bc)

        integer, intent(in) :: bc

        f_l0_bc_unsupported = (bc < BC_GHOST_EXTRAP)

    end function f_l0_bc_unsupported

    !> Slot index of regrid-managed fine block k in the shared pool: tiles occupy [1..l0_slot_off], fine blocks [l0_slot_off+1..].
    !! Identity (l0_slot_off=0) until L0 tiles + AMR coexist.
    pure integer function f_l0_slot(k) result(s)

        integer, intent(in) :: k

        s = l0_slot_off + k

    end function f_l0_slot

    !> Build the base-grid tiling: allocate the slot/region/seam machinery for l0_ntiles_tot rr=1 tiles covering L0, set each tile's
    !! geometry (extent, L0-slice coords, whole-tile footprint, owner) and copy the current L0 state in. Standalone (does not call
    !! s_initialize_amr_module) so the spike lifts out cleanly. np=1 spike: rank 0 owns every tile.
    impure subroutine s_l0_tiles_init()

        integer :: nt(3), ix, iy, iz, k, j, r, e
        integer :: tlo(3), thi(3)
        integer :: rsidx(3), rext(3)
        integer :: ierr

        if (l0_ntile <= 0) return

        ! periodic_bc is set on rank 0 only (s_read_input_file is rank-0-guarded), so make it globally consistent: every rank must
        ! build
        ! the SAME wrap-seam list in f_amr_seam / apply the same periodic edge fill, else an unmatched MPI_SENDRECV deadlocks.
        ! MPI_LOR:
        ! rank 0's .true. wins on all ranks (others hold the .false. default). No-op at np=1.
        l0_periodic = periodic_bc
#ifdef MFC_MPI
        call MPI_ALLREDUCE(MPI_IN_PLACE, l0_periodic, 3, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
#endif

        ! Supported physical faces (any np): extrapolation (BC_GHOST_EXTRAP), reflective (BC_REFLECTIVE), periodic (BC_PERIODIC) -
        ! each
        ! has a cons-space tile fill that commutes with the cons->prim convert so a tile matches the monolithic prim-space BC
        ! bit-for-
        ! bit. The characteristic/slip/dirichlet family (bc < BC_GHOST_EXTRAP) is not yet handled. Validate once here (host), not
        ! per
        ! stage.
        if (f_l0_bc_unsupported(bc_x%beg) .or. f_l0_bc_unsupported(bc_x%end) .or. (n_glb > 0 .and. (f_l0_bc_unsupported(bc_y%beg) &
            & .or. f_l0_bc_unsupported(bc_y%end))) .or. (p_glb > 0 .and. (f_l0_bc_unsupported(bc_z%beg) &
            & .or. f_l0_bc_unsupported(bc_z%end)))) then
            call s_mpi_abort('l0_ntile spike: unsupported physical BC (only extrapolation, reflective, periodic are handled)')
        end if
        ! the monolithic L0 RHS is skipped for l0_ntile>0, so the global q_prim_vf it populated is stale: run-time-info and probes
        ! (which read it at stage 1) are not supported in the spike.
        if (run_time_info .or. probe_wrt) then
            call s_mpi_abort('l0_ntile spike does not support run_time_info or probe_wrt (monolithic q_prim_vf is not maintained)')
        end if

        ! rr=1 makes the swap ghost-coord bisection and the fine-fine-halo fmul (=amr_ref_ratio**level) both identity; slots inherit
        ! it via s_amr_alloc_slot. Only clobber the GLOBAL in pure-L0 (amr off): under coexist the global must stay the real 2/4 so
        ! fine blocks size correctly (s_set_amr_fine_geometry etc.), and tiles get rr=1 via the per-slot override in
        ! s_l0_build_tile_slot (Option 2: level-0 tiles are rr=1 regardless of the global).
        if (.not. amr) amr_ref_ratio = 1

        ! Tiles are PER-RANK: each rank's local chunk is split into nt(:) pieces; the global tile table (region + owner) is the
        ! union
        ! over ranks, REPLICATED on every rank. Total = num_procs * nt(1)*nt(2)*nt(3). At np=1 this is identical to the old global
        ! tiling. Each rank allocates slot DATA (fields + coords) only for its OWN tiles; the seam-pair scan and fine-fine halo see
        ! the
        ! full table and exchange cross-rank seams over MPI.
        ! l0_nt/l0_ntiles_tot/l0_slot_off are computed once by s_initialize_amr_module (which always runs first, per
        ! m_start_up.fpp) so both inits agree on the shared-pool layout; just read them here.
        nt = l0_nt

        amr_num_levels = 1
        amr_cur = 1

        if (.not. amr) then
            ! l0-only mode: this routine owns the pool exactly as before (fine budget = 0)
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
            ! L0 tiles are the BASE level (Option 2: fine blocks become level>=1 on a tile)
            allocate (amr_block_level(amr_max_blocks)); amr_block_level = 0
            ! 2D rank lists sized to the computed max overlap in s_amr_build_seam_pairs; only the per-block counts are sized here.
            allocate (amr_ovl_gather_n(amr_max_blocks), amr_ovl_scatter_n(amr_max_blocks))
            allocate (amr_slot_live(amr_max_blocks)); amr_slot_live = .false.
            amr_region_lo_all = 0; amr_region_hi_all = 0; amr_isect_lo_all = 0; amr_isect_hi_all = 0; amr_owns_all = .false.
        else
            ! coexist mode: s_initialize_amr_module already allocated the shared pool sized l0_slot_off+amr_max_fine.
            ! Only allocate the TILE-specific side tables here, and do NOT touch amr_slots / amr_region_* / amr_owns_all /
            ! amr_block_owner / amr_ovl_* - those are shared with AMR and already sized/allocated.
            allocate (amr_tile_l0_owner(amr_max_blocks)); amr_tile_l0_owner = 0
            allocate (amr_tile_cost(amr_max_blocks)); amr_tile_cost = 0._wp
            allocate (amr_tile_cost_ema(amr_max_blocks)); amr_tile_cost_ema = 0._wp
            ! tiles are level 0 in slots [1..l0_ntiles_tot]; set that band without disturbing the fine slots
            amr_block_level(1:l0_ntiles_tot) = 0
        end if

        ! the per-rank coarse decomposition (global origin + local extent) that the tile geometry and max-tile-extent sizing below
        ! read for every rank is computed O(1) by s_amr_rank_decomp - no table, no allgather. In l0-only mode
        ! s_initialize_amr_module
        ! did not run, so validate the formula against this rank's actual decomposition here (coexist validates it in that routine).
        if (.not. amr) call s_amr_validate_decomp()

        ! max tile extent per dim over ALL ranks (= widest per-rank chunk split by nt); slots + seam buffers are sized to this
        ! global
        ! max so every rank's buffers match. Chunk r has s_amr_rank_decomp ext(d)+1 cells in dim d.
        ! Coexist: the FINE sizing s_initialize_amr_module just computed must SURVIVE - these are module-level and are what
        ! s_amr_alloc_slot reads, so zeroing them here leaves the shared pool sized to the tile extent. The static block's slot
        ! is allocated before this routine and so escapes, but every slot a REGRID allocates afterwards is a fine block cut to
        ! tile size, which then writes past its own bounds (an out-of-range device write, not a host abort). Accumulate the max
        ! of both instead. l0-only (.not. amr): s_initialize_amr_module returned early, so no fine sizing exists - start at 0.
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
        call s_amr_scr_init()  ! mbuf* now FINAL (fine/tile union under coexist); scratch must exist on every rank

        amr_seam_pairs_dirty = .true.; amr_seam_pairs_nblk = -1
        amr_mesh_epoch = amr_mesh_epoch + 1

        ! swap bounce buffers (same bounds as the L0 global coord arrays). SHARED with s_initialize_amr_module (identical m/n/p
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
        call s_l0_build_extended_global_cb()  ! global L0 boundaries EXTENDED into the domain ghost shell (edge tiles need it)
        amr_cpat_mar = (buff_size + amr_ref_ratio - 1)/amr_ref_ratio + 1
        amr_xchg_coarse_ghosts = .false.  ! tiles never prolong from a coarser level

        ! per-tile geometry: level-1 rr=1 blocks. Loop is rank-major then z,y,x (x fastest) so intra-rank neighbours are contiguous;
        ! the global region scan finds cross-rank seams. A rank writes region+owner for EVERY tile but allocates slot data only for
        ! its
        ! own (r == proc_rank). Domain-edge detection uses the global region indices (region_lo == 0 / region_hi == m_glb).
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
                        amr_block_owner(k) = r
                        amr_tile_l0_owner(k) = r  ! L0 storage owner = init owner; stays fixed under migration
                        amr_owns_all(k) = (r == proc_rank)
                        amr_region_lo_all(:,k) = tlo; amr_region_hi_all(:,k) = thi
                        amr_isect_lo_all(:,k) = tlo; amr_isect_hi_all(:,k) = thi  ! footprint = whole tile on the owner (rr=1)
                        ! slot data is NOT allocated here: it must follow the SFC COMPUTE owner assigned below, which need not be
                        ! this cartesian owner r. Every rank writes region+owner metadata for every tile.
                    end do
                end do
            end do
        end do

        ! COMPUTE owner = SFC cost-split over the tiles (fills the O(num_procs) amr_owner_cut), superseding the cartesian owner set
        ! in
        ! the loop above; amr_tile_l0_owner keeps the cartesian STORAGE assignment (where the L0 field data physically lives). Init
        ! cost is geometric (whole-tile cell count), uniform on a uniform grid. The SFC compute owner may differ from the cartesian
        ! storage owner - whether it does depends on how the cartesian split direction lines up with Morton order, so it is a
        ! property of the grid SHAPE (a 2:1 grid at np=2 splits in y and agrees; a square grid tie-breaks to x and diverges).
        ! s_l0_copy_coarse_to_tiles routes the initial fill storage-owner -> compute-owner when they differ, so no precondition is
        ! asserted here; the rebalancer's routed migration already handles post-init divergence.
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
                amr_block_owner(kk) = sfco(kk)
                amr_owns_all(kk) = (sfco(kk) == proc_rank)
            end do
            ! allocate slot data for the tiles this rank COMPUTES (deferred from the cartesian loop above). s_l0_build_tile_slot
            ! reads only replicated region metadata and the global amr_g?cb, so it is valid for any tile on any rank. When the SFC
            ! cut agrees with the cartesian order this allocates exactly the same set as before, just later.
            do kk = 1, l0_ntiles_tot
                if (amr_block_owner(kk) == proc_rank) call s_l0_build_tile_slot(kk)
            end do
        end block
        ! validate the full picture now that tiles exist: tiles vs amr_owner_cut (tile cut, just built), fine blocks vs amr_fine_cut
        ! (the level-1 cut the assigner saved at init). Runs in coexist too (amr_fine_cut is populated by the earlier assigner
        ! call).
        call s_amr_validate_owner()

        call s_amr_select_slot(1)
        ! tiles are PERSISTENT: L0 seeds them once at the first timestep (s_l0_copy_coarse_to_tiles self-gates on this flag), then
        ! they carry their own state across stages/timesteps. No init-time device copy (q_cons device state is not live at module
        ! init).
        l0_tiles_need_fill = .true.

    end subroutine s_l0_tiles_init

    !> Global L0 cell boundaries EXTENDED into the domain ghost shell (-1-buff_size : G+buff_size), unlike s_amr_build_global_cb
    !! (-1:G). The swap rebuilds a block's ghost-shell coordinates from these, and a tile that touches the domain boundary reaches
    !! indices beyond G (AMR fine blocks never do, being buff_size inside). Sourced from the monolithic x_cb, whose ghost cells
    !! already hold the domain's ghost coordinates - so a tile's ghost coords match the monolithic grid's bit-for-bit.
    impure subroutine s_l0_build_extended_global_cb()

        integer             :: j
        real(wp), parameter :: sentinel = -huge(1._wp)

        ! Under coexist, s_amr_build_global_cb (called by s_initialize_amr_module) already allocated amr_g?cb at the NON-extended
        ! bounds (-1:G) for the fine-block geometry. The tiles need the EXTENDED bounds (-1-buff:G+buff) since an edge tile reaches
        ! into the domain ghost shell; the extended array is a value-consistent superset (same x_cb source over the overlap, and the
        ! fine geometry already copied its coords into the slots), so replace it. Without this the second allocate is a fatal error
        ! on gfortran and a silent double-allocate (leaked non-extended buffer) on flang.

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

    !> Build tile k's slot on THIS rank from its (already-set, replicated) region metadata: allocate the field/coord arrays and set
    !! the local extents, idwbuff, and rr=1 cell coordinates sliced from the global amr_g?cb. Shared by s_l0_tiles_init (initial
    !! owned tiles) and s_l0_migrate_tile (a tile arriving on its new owner). Requires amr_gxcb/gycb/gzcb + mbuf*/max_f* already
    !! set.
    impure subroutine s_l0_build_tile_slot(k)

        integer, intent(in) :: k
        integer             :: j, tlo(3), thi(3)

        tlo = amr_region_lo_all(:,k); thi = amr_region_hi_all(:,k)
        call s_amr_alloc_slot(k)  ! sizes to mbuf*, sets slot%amr_ref_ratio = amr_ref_ratio
        ! a base-level tile is rr=1 regardless of the global refinement ratio (Option 2: global may be 2/4 for fine blocks)
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
        ! rr=1: tile cell j (right boundary) IS the global L0 boundary amr_g?cb(tlo + j); interior coords only (the swap extends the
        ! ghost shell from amr_g?cb identically).
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
    !! compute owner is also its L0-storage owner is seeded by a local device copy (the common case, and the ENTIRE path when the
    !! SFC cut agrees with the cartesian order, so byte-identical to before). When the SFC compute owner differs from the cartesian
    !! storage owner the seed is ROUTED: the L0-storage owner device-packs its chunk and sends it to the compute owner, which
    !! unpacks into its tile slot. Exact reverse of s_l0_scatter_tiles_to_coarse, and sound because a tile is built by subdividing
    !! ONE rank's cartesian chunk (s_l0_tiles_init), so it never spans two L0-storage ranks.
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
    !! s_l0_copy_coarse_to_tiles because the coexist SUBCYCLE path round-trips through L0 every step (tiles -> L0, fine fold writes
    !! L0, L0 -> tiles), so it needs this unconditionally, while the seed must still happen exactly once.
    impure subroutine s_l0_fill_tiles_from_coarse(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        integer                                                :: k, o1, o2, o3, fm1, fm2, fm3, bown, lown, cnt, ierr
        real(wp), allocatable                                  :: buf(:)

        do k = 1, l0_ntiles_tot
            bown = amr_block_owner(k); lown = amr_tile_l0_owner(k)
            if (bown == lown) then  ! compute owner holds the L0 cells: local device copy (unchanged path)
                if (bown /= proc_rank) cycle
                call s_l0_tile_l0_offsets(k, o1, o2, o3)
                fm1 = amr_slots(k)%m; fm2 = amr_slots(k)%n; fm3 = amr_slots(k)%p
                call s_l0_copy_block(amr_loc_of(k), q_cons_vf, o1, o2, o3, fm1, fm2, fm3, .true.)
                cycle
            end if
            ! routed seed: extents come from the REPLICATED region (the L0 owner has no slot for this tile)
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
    !! also its L0-storage owner writes locally (device kernel - the common case, and the ENTIRE no-migration path, so
    !! byte-identical to before). A MIGRATED tile (owner != l0_owner) has its interior sent by the compute owner to the L0-storage
    !! owner over MPI, which writes it into L0 - keeping the fixed L0 decomposition (hence output/restart) correct after migration.
    !! Ghosts are not scattered (the tile path never reads L0 ghosts). GPU-correct: the MPI branch device-packs/unpacks via
    !! s_l0_pack_unpack_block, so the receiver writes L0 ON THE DEVICE - it survives the GPU_UPDATE(host) that s_save_data does
    !! before writing.
    impure subroutine s_l0_scatter_tiles_to_coarse(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        integer                                                :: k, o1, o2, o3, fm1, fm2, fm3, bown, lown, cnt, ierr
        real(wp), allocatable                                  :: buf(:)

        ! Precondition: tiles are the authoritative store. Before the first seed (s_l0_copy_coarse_to_tiles) the tile slots hold
        ! uninitialized (zero) state and L0 still holds the initial condition, so there is nothing to refresh - scattering here
        ! would overwrite the IC with zeros (zero density -> NaN once the coexist L0 coarse RHS consumes it). Skip until seeded.

        if (l0_tiles_need_fill) return

        do k = 1, l0_ntiles_tot
            bown = amr_block_owner(k); lown = amr_tile_l0_owner(k)
            fm1 = amr_region_hi_all(1, k) - amr_region_lo_all(1, k)
            fm2 = 0; if (n_glb > 0) fm2 = amr_region_hi_all(2, k) - amr_region_lo_all(2, k)
            fm3 = 0; if (p_glb > 0) fm3 = amr_region_hi_all(3, k) - amr_region_lo_all(3, k)
            if (bown == lown) then  ! not migrated: local device copy (unchanged path)
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

    !> Device ADD of an L0 block [o+0:o+fm] into a tile rhs interior [0:fm] (q_rhs += q_l0). Additive twin of s_l0_copy_block's
    !! to_tile branch, in wp (the rhs is computed in wp, stored stp). Slot rhs is a dummy so the kernel reads a valid mapped
    !! descriptor (indexing module amr_slots%rhs in a kernel is a null deref - see s_l0_copy_block / s_amr_fine_slice).
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

    !> Device ADD of the contiguous MPI buffer buf into a tile rhs interior [0:fm] (q_rhs += buf). Additive twin of
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

    !> Coexist reflux copy-back: ADD the fixed-L0-frame Berger-Colella reflux delta into each tile's per-slot rhs on its (possibly
    !! migrated) compute owner, so the tile RK update sees the same c/f-face correction the monolithic coarse update would. The
    !! delta lives in rhs_delta (the L0 rhs, which s_tvd_rk zeroed before the fine loop so s_amr_apply_reflux filled it with the
    !! PURE delta). Reverse of s_l0_scatter_tiles_to_coarse: source is the fixed L0-storage owner (amr_tile_l0_owner), dest is the
    !! compute owner (amr_block_owner); local when they coincide, else P2P (L0-owner packs the tile's L0 region, compute-owner adds
    !! it). Whole tile interior is routed - the delta is zero outside the c/f reflux shell, so the add is identity elsewhere.
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
            else if (proc_rank == bown) then  ! compute owner: recv, device-ADD the delta into the tile rhs
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
    !! OVERWRITE each covering tile's matching cells with those restricted values on the tile's (possibly migrated) compute owner -
    !! the coexist twin of the monolithic level-0 covered-cell overwrite. Only the covered footprint moves (non-covered tile cells
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
    !! receiver, free it on the sender, and update the replicated owner map + seam topology. ALL ranks call with the same (k,
    !! new_owner). This is the load-balance MIGRATION primitive; the DECISION of which tile moves where is made by the caller (a
    !! forced remap in the spike; a cost-driven trigger later). Ghosts are not moved (refilled by edge-BC + fine-fine halo before
    !! the next stage). GPU-correct: interior is device-packed/unpacked via s_l0_pack_unpack_block (wp buffer, cast to/from stp).
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

        ! replicated ownership update on EVERY rank; mark the seam topology dirty so the next halo rebuilds pair/overlap lists.
        ! The epoch bump matters most HERE: ownership changed with NO regrid, which the consumed boolean cannot express to a
        ! cached exchange plan.
        amr_block_owner(k) = new_owner
        amr_owns_all(k) = (new_owner == proc_rank)
        amr_seam_pairs_dirty = .true.
        amr_mesh_epoch = amr_mesh_epoch + 1

    end subroutine s_l0_migrate_tile

    !> Spike test hook: at t_step == l0_migrate_step, force-migrate the LAST tile (initially owned by rank num_procs-1) to rank 0,
    !! exercising the migration primitive + seam-topology rebuild. Output must stay byte-identical to the no-migration run. No-op at
    !! np=1 (the last tile already lives on rank 0). ALL ranks call with identical arguments.
    impure subroutine s_l0_forced_remap()

        integer :: k

        k = l0_ntiles_tot
        if (amr_block_owner(k) /= 0) call s_l0_migrate_tile(k, 0)

    end subroutine s_l0_forced_remap

    !> Closed-loop rebalancer driven by MEASURED per-tile compute time. Each rank accumulated amr_tile_cost for its OWN tiles since
    !! the last rebalance; an allreduce(SUM) makes the full cost vector REPLICATED and bit-identical on every rank (each tile has
    !! exactly one nonzero contributor, so the sum is exact) -> every rank runs the identical greedy and issues MATCHING P2P
    !! migrations. Greedy: move the tile on the heaviest rank that most reduces the max-min load gap onto the lightest, bounded,
    !! with a small relative deadband so timing noise does not cause churn. Because migration is bit-preserving and the decision
    !! touches no field data, OUTPUT is byte-identical regardless of the (run-to-run nondeterministic) measured schedule - the
    !! decomposition-invariance proven for the tile path is exactly what keeps a nondeterministic cost signal golden-safe. Costs
    !! reset after each rebalance. No-op at np=1.
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

        ! SFC weighted re-cut (replaces the greedy min-max mover): Morton-sort the tiles, then cumulative-split the smoothed cost
        ! into num_procs contiguous SFC ranges - identical partition logic to s_amr_assign_block_owners' cut. Ownership stays
        ! SFC-contiguous and O(num_procs) cut-derivable (the precondition for retiring the amr_block_owner table), and locality is
        ! preserved. Deadband kept: skip the re-cut while the load gap is already within tol (no churn on sub-5% imbalance).
        if (gap0 > tol) then
            block
                integer(kind=8) :: tkey(l0_ntiles_tot), cut_try(0:num_procs - 1)
                integer         :: newo_try(l0_ntiles_tot)
                real(wp)        :: load_try(0:num_procs - 1)
                do k = 1, l0_ntiles_tot
                    tkey(k) = f_morton(amr_region_lo_all(1, k), amr_region_lo_all(2, k), amr_region_lo_all(3, k))
                end do
                ! shared SFC cut: cumulative-split the smoothed cost into num_procs contiguous Morton ranges. Because the cut is
                ! restricted to CONTIGUOUS ranges, the finest correction it can make is one whole tile; when tiles-per-rank is small
                ! that quantum exceeds the gap the deadband admits, and the re-cut can return a partition WORSE than the current one
                ! (measured 5.738E-02 -> 3.898E-01 at 8 tiles / 2 ranks, and it never recovered). Evaluate into a temporary and
                ! reject only a STRICT WORSENING. Not "accept only a strict improvement": with near-uniform tile costs the Morton
                ! partition and the cartesian one have EQUAL gaps, and rejecting those would silently kill the migration that
                ! golden "L0 tiles -> 2D -> SFC re-cut rebalance np=2" exists to cover (it would still PASS, because migration is
                ! bit-neutral - the coverage would just be gone). amr_owner_cut must move WITH newo: f_amr_owner resolves tile
                ! ownership against it, so refreshing it without migrating would leave it disagreeing with amr_block_owner.
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
        ! tol is printed because without it "deadband skipped the re-cut" (gap0 <= tol) and "the re-cut ran and was REJECTED for not
        ! improving the gap" (gap0 > tol, guard above) emit byte-identical output - gap0 -> gap0, 0 migrations. Those are different
        ! behaviours and one of them is the guard doing its job; a verification run cannot tell them apart otherwise.
        if (proc_rank == 0) print '(A,I0,A,ES10.3,A,ES10.3,A,ES10.3,A,I0,A)', ' [l0 rebalance] t_step=', t_step, ' load-gap ', &
            & gap0, ' -> ', gap1, ' (deadband ', tol, ', ', nmig, ' migrations)'

        amr_tile_cost = 0._wp  ! reset the measurement window

    end subroutine s_l0_rebalance

    !> Device copy between a tile interior [0:fm] and the L0 field [o+0:o+fm]. to_tile=T copies L0->tile, F copies tile->L0. The
    !! slot is addressed through the flat store, which is a plain GPU_DECLARE'd module array (the old per-slot layout made a null
    !! deref - see s_amr_fine_slice).
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

    !> Device pack (to_buf=T) / unpack (F) of a field's interior block [o+0:o+fm] <-> the contiguous MPI buffer buf, for the P2P
    !! migration + migrated-tile scatter. Follows s_amr_fine_slice: the pack/unpack runs ON THE DEVICE with copyout/copyin moving
    !! only buf host<->device (no strided %sf section in a map clause - flang miscomputes those). buf index runs j fastest then
    !! k,l,i so a matching pack/unpack aligns cell-for-cell. wp buffer, cast to/from stp (identity at double) - matches the
    !! fine-fine halo. Two targets, one body: `_st` packs a block out of the flat store (the migration/scatter paths), `_sf` packs
    !! the level-0 monolithic q_cons_vf, which is a real scalar_field array and not in the store.
    #:for SFX, TGT in [('st', 'amr_cons_st'), ('sf', '')]
        #:set QB = (lambda ix: TGT + '(o1 + j, o2 + k, o3 + l, ' + ix + ', loc)') if TGT else (lambda ix: 'q(' + ix &
                    & + ')%sf(o1 + j, o2 + k, o3 + l)')
        impure subroutine s_l0_pack_unpack_block_${SFX}$(${'loc' if TGT else 'q'}$, o1, o2, o3, fm1, fm2, fm3, buf, to_buf)

            #:if TGT
                integer, intent(in) :: loc
            #:else
                type(scalar_field), dimension(sys_size), intent(inout) :: q
            #:endif
            integer, intent(in)                 :: o1, o2, o3, fm1, fm2, fm3
            real(wp), intent(inout), contiguous :: buf(:)
            logical, intent(in)                 :: to_buf
            integer                             :: i, j, k, l

            if (to_buf) then
                $:GPU_PARALLEL_LOOP(collapse=4, copyout='[buf]')
                do i = 1, sys_size
                    do l = 0, fm3
                        do k = 0, fm2
                            do j = 0, fm1
                                buf(1 + j + (fm1 + 1)*(k + (fm2 + 1)*(l + (fm3 + 1)*(i - 1)))) = real(${QB('i')}$, wp)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else
                $:GPU_PARALLEL_LOOP(collapse=4, copyin='[buf]')
                do i = 1, sys_size
                    do l = 0, fm3
                        do k = 0, fm2
                            do j = 0, fm1
                                ${QB('i')}$ = real(buf(1 + j + (fm1 + 1)*(k + (fm2 + 1)*(l + (fm3 + 1)*(i - 1)))), stp)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

        end subroutine s_l0_pack_unpack_block_${SFX}$
    #:endfor

    !> Fill each tile's DOMAIN-EDGE face ghosts with the physical BC (interior-seam faces are overwritten by s_amr_fine_fine_halo
    !! afterward). Milestone 1 supports extrapolation (bc <= BC_GHOST_EXTRAP); other codes abort. A face (d,side) is a domain edge
    !! iff the tile touches the global boundary there.
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
    !! the prim ghost the monolithic path produces. The transverse loop spans the FACE interior only (dimension-split scheme reads
    !! no corner ghost).
    impure subroutine s_l0_edge_bc_tile(loc, rlo, rhi, gcell, fm, d, bcbeg, bcend)

        integer, intent(in) :: loc
        integer, intent(in) :: rlo, rhi, gcell, fm(3), d, bcbeg, bcend

        ! BC support is validated once at init (s_l0_tiles_init); here we only apply it at domain-edge faces. Periodicity is read
        ! from
        ! the global periodic_bc(d) (not bcbeg, which becomes a wrap-neighbour RANK at a decomposed periodic boundary): a periodic
        ! dim
        ! wraps - a tile that SPANS it (rlo==0 .and. rhi==gcell) self-wraps here, a partial tile's periodic faces are cross-tile
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

    !> Reflective (symmetry) tile face ghosts in dim d, side (-1 low / +1 high): ghost cell 1..buff_size MIRRORS the near-edge
    !! interior (ghost -jg <- interior jg-1 low; md+jg <- md-(jg-1) high) with the NORMAL-direction momentum (eqn_idx%mom%beg + d -
    !! 1) negated, all other conserved variables copied. Done on q_cons; negating conserved normal momentum commutes with the
    !! cons->prim convert (velocity flips, rho and mom**2 - hence pressure - are unchanged), so this reproduces the monolithic
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

    !> Periodic self-wrap for a tile that SPANS dim d (its low AND high faces are both the domain boundary, i.e. l0_ntile==1 in d):
    !! fill both ghost shells from the opposite-end interior of the SAME tile - low ghost -jg <- interior md-(jg-1), high ghost
    !! md+jg <- interior jg-1 (a pure copy, matching the monolithic prim-space s_periodic; copy commutes with cons->prim convert).
    !! Partial tiles never reach here (their periodic faces are cross-tile wrap-seams handled by s_amr_fine_fine_halo). Face
    !! interior only.
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

    !> Fill a tile's MULTI-DIM ghost cells (>= 2 dims in the ghost region: 2D diagonal corners, 3D ghost edges + corners) from the
    !! nearest INTERIOR cell (each index clamped to [0:m]/[0:n]/[0:p]). The dimension-split face fills only set single-ghost face
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
    !! fine-block phase structure (m_time_steppers) with prolong/reflux dropped - a base-res tile has no coarser level.
    !> Advance every owned tile one RK stage. Fused wrapper (pure-L0): RHS pass then RK pass back-to-back = byte-identical to the
    !! single-pass form. Under coexist (amr) s_tvd_rk instead calls the two passes directly with s_l0_add_reflux_to_tiles between
    !! them, so each tile's coarse rhs is Berger-Colella corrected (from the fixed-L0-frame reflux) before its RK update.
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
        ! stages, reset at each rebalance. Timing is a pure side-channel - it never touches field data, so output stays
        ! bit-identical.

        measure = (l0_rebalance_interval > 0)

        call s_l0_fill_edge_bc()
        call s_amr_fine_fine_halo(0)
        ! Fill the multi-dim ghost cells (2D diagonal corners; 3D also the ghost edges) that the dimension-split face fills
        ! (s_l0_fill_edge_bc + s_amr_fine_fine_halo) deliberately leave unset. The RHS never reads them, but the cons->prim
        ! convert processes the whole buffered range, and an unset ghost (all void fractions 0 -> gamma 0) is a 0/0 that traps
        ! under -ffpe-trap and is a stray NaN otherwise. Each is copied from the nearest INTERIOR cell (valid state; the face
        ! ghosts the RHS reads are single-ghost and untouched, so output is unchanged).
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
            ! tiles fill their OWN rhs (it must survive the whole-set RHS pass through the reflux point to the RK pass); the
            ! per-slot q_prim exists exactly when the copy-out gate writes it - otherwise the pooled scratch takes the (unread,
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
        ! Only free the shared slot pool in pure-L0 mode: under coexist (amr) s_finalize_amr_module runs FIRST (m_start_up call
        ! order) and already freed every slot 1..amr_max_blocks AND deallocated amr_slot_live, so re-running s_amr_free_slot here
        ! would read the deallocated amr_slot_live (use-after-free).
        if (.not. amr) then
            do islot = 1, amr_max_blocks
                call s_amr_free_slot(islot)
            end do
        end if
        if (allocated(amr_seam_pairs)) deallocate (amr_seam_pairs)
        ! amr_slots, amr_region_*, amr_isect_*, amr_owns_all, amr_block_owner, amr_block_level, amr_ovl_*, and
        ! amr_slot_live are SHARED with s_initialize_amr_module/s_finalize_amr_module: when amr, that pair owns them, so only
        ! free them here in l0-only mode to avoid a coexist double-free. amr_tile_l0_owner/amr_tile_cost/amr_tile_cost_ema are
        ! TILE-ONLY and always freed here.
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

    impure subroutine s_finalize_amr_module()

        integer :: i, islot

        ! BEFORE the amr early-return: the report's conservation allreduce is collective, and the L0 tile
        ! families can fire with amr = F. All ranks take the same path either way.

        call s_xa_report()
        call s_amr_cov_report()
        if (.not. amr) return
        do islot = 1, amr_max_blocks
            call s_amr_free_slot(islot)
        end do
        if (qbmm .and. .not. polytropic) then
            @:DEALLOCATE(amr_rhs_pb_f)
            @:DEALLOCATE(amr_rhs_mv_f)
            @:DEALLOCATE(amr_cg_pb)
            @:DEALLOCATE(amr_cg_mv)
        end if
        deallocate (amr_slot_live)
        call s_amr_st_finalize()
        if (allocated(amr_seam_pairs)) deallocate (amr_seam_pairs)
        if (allocated(amr_ovl_gather)) deallocate (amr_ovl_gather)
        if (allocated(amr_ovl_scatter)) deallocate (amr_ovl_scatter)
        deallocate (amr_ovl_gather_n, amr_ovl_scatter_n)
        if (allocated(amr_gpl_nsrc)) deallocate (amr_gpl_nsrc, amr_gpl_src, amr_gpl_sz, amr_gpl_psrc, amr_gpl_psz)
        if (allocated(amr_gcr_pool)) deallocate (amr_gcr_pool)
        if (allocated(amr_gcr_req)) deallocate (amr_gcr_req, amr_gcr_off)
        ! per-array guards, NOT grouped on a lead member: the wave-scratch arrays of one group allocate
        ! independently (spsz/rpsz are sized only by the qbmm pb/mv wave branch), so a non-qbmm np>1 run
        ! reaches here with a group partially allocated. gfortran/ifx abort on deallocating an unallocated
        ! array (amdflang silently tolerates it) - the CI probe's "Restart roundtrip run failed" crash class.
        #:for A in ['amr_fw_sblk', 'amr_fw_sbl', 'amr_fw_sbh', 'amr_fw_spi', 'amr_fw_sqo', 'amr_fw_spo', &
            'amr_fw_rblk', 'amr_fw_rbl', 'amr_fw_rbh', 'amr_fw_rpi', 'amr_fw_rqo', 'amr_fw_rpo', &
            'amr_fw_sprank', 'amr_fw_sqsz', 'amr_fw_spsz', 'amr_fw_snxp', 'amr_fw_sqbase', 'amr_fw_spbase', &
            'amr_fw_rprank', 'amr_fw_rqsz', 'amr_fw_rpsz', 'amr_fw_rnxp', 'amr_fw_rqbase', 'amr_fw_rpbase', &
            'amr_fw_map', 'amr_fw_nx', 'amr_fw_pq', 'amr_fw_pp']
            if (allocated(${A}$)) deallocate (${A}$)
        #:endfor
        if (allocated(amr_fw_sq)) deallocate (amr_fw_sq)
        if (allocated(amr_fw_sp)) deallocate (amr_fw_sp)
        if (allocated(amr_fw_rq)) deallocate (amr_fw_rq)
        if (allocated(amr_fw_rp)) deallocate (amr_fw_rp)
        if (allocated(amr_fw_req)) deallocate (amr_fw_req, amr_fw_reqw)
        do i = 1, sys_size
            @:DEALLOCATE(amr_cg(i)%sf)
        end do
        @:DEALLOCATE(amr_cg)
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
        if (cyl_coord .and. n_glb > 0) then
            @:DEALLOCATE(amr_rvw)
        end if

    end subroutine s_finalize_amr_module

end module m_amr
