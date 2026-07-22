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
    use m_constants, only: num_fluids_max, model_eqns_6eq, mapCells, K_ib, K_pc
    use m_pressure_relaxation, only: s_pressure_relaxation_procedure
    use m_mpi_proxy, only: s_mpi_abort
    use m_mpi_common, only: s_mpi_allreduce_integer_min, s_mpi_allreduce_integer_max, s_mpi_allreduce_sum, s_mpi_allreduce_min, &
        & s_mpi_allreduce_max, s_mpi_allreduce_integer_sum, s_mpi_sendrecv_variables_buffers, s_mpi_allreduce_array_max
    use m_rhs, only: s_compute_rhs
    use m_phase_change, only: s_infinite_relaxation_k, pc_iter_count
    use m_amr_registers, only: s_amr_zero_fine_registers, s_amr_reflux_apply_faces, freg, creg
    use m_rank_timing, only: s_rank_time_tic, s_rank_time_toc
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
        & s_restrict_fine_to_coarse, s_finalize_amr_module, s_amr_fine_stage_fill, s_amr_fine_stage_advance, &
        & s_amr_fine_fine_halo, s_amr_advance_fine_subcycle_all, s_set_amr_fine_geometry, s_amr_relax_fine, s_amr_setup_ib, &
        & s_amr_p2p_reflux_faces, s_amr_reflux_to_parent
    ! s_amr_swap_to_fine / s_amr_restore_coarse / s_amr_fill_fine_ghosts / amr_dt_fine are internal (no external caller); keeping
    ! them private makes "the swap has exactly these audited call sites" a compiler guarantee, not a convention.
    !> Block/slot state and fine-distribution services consumed by m_amr_regrid and m_amr_restart (the drivers split out of this
    !! module). State stays HERE - only the drivers moved.
    public :: amr_slots, amr_seam_pairs_dirty, amr_xchg_coarse_ghosts, amr_cpat_mar, s_amr_alloc_slot, s_amr_reconcile_slots, &
        & s_amr_assign_block_owners, s_amr_gather_coarse_patch, s_amr_gather_coarse_patch_pbmv, s_amr_prolong_pbmv, &
        & s_amr_exchange_coarse_cons_halo, s_lag_phys_to_cells, s_amr_body_bbox, s_amr_expand_box_over_bodies, s_amr_tile_box, &
        & f_amr_seam_dim, f_amr_boxes_overlap

    !> Fine-level time step for subcycling (= 0.5*dt after init; 0 when amr is off).
    real(wp) :: amr_dt_fine = 0._wp

    !> Realizability floor for prolonged Euler-Euler bubble POSITIVE moments (radius nR, non-polytropic partial pressure npb / vapor
    !! mass nmv): a positive fraction of the coarse parent so derived R = nR/n, pb, mv stay >= 0. Minmod keeps a positive field
    !! positive, so this fires only under floating-point edge cases (conservation defect ~0 otherwise).
    real(wp), parameter :: bub_pos_frac = 1.0e-10_wp

    !> One refined level: its own grid + conservative fields. Field arrays are device-resident (@:ALLOCATE); coords/metadata
    !! host-only.
    type t_level
        integer                         :: amr_ref_ratio
        type(t_box)                     :: region          !< block extent in parent (level-0) cell indices
        integer                         :: m, n, p         !< this level's interior extents
        integer                         :: buff_size
        type(int_bounds_info)           :: idwbuff(3)
        real(wp), allocatable           :: x_cb(:), x_cc(:), dx(:)
        real(wp), allocatable           :: y_cb(:), y_cc(:), dy(:)
        real(wp), allocatable           :: z_cb(:), z_cc(:), dz(:)
        type(scalar_field), allocatable :: q_cons(:)
        type(scalar_field), allocatable :: q_cons_stor(:)  !< stage storage (fine advance, same bounds as q_cons)
        type(scalar_field), allocatable :: q_prim(:)       !< primitive stage (fine advance, same bounds as q_cons)
        type(scalar_field), allocatable :: rhs(:)          !< RHS (fine interior only: 0:m, 0:n, 0:p)
        type(scalar_field), allocatable :: q_ghost_a(:)    !< subcycle ghost-lerp source at coarse t^n (ghost shell only)
        type(scalar_field), allocatable :: q_ghost_b(:)    !< subcycle ghost-lerp source at coarse t^{n+1} (ghost shell only)
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
    integer                    :: amr_maxc(3)  !< max coarse block cells per dim: (m_glb+1)/2 etc.; 1 for collapsed dims

    !> Per-slot field-array sizing (module-scope, used by s_amr_alloc_slot/s_amr_free_slot): max fine cells per dim (2*maxc_loc-1)
    !! and the buffered array bounds. amr_slot_live(k) tracks whether slot k's field arrays are allocated - lazy owned-only sizing
    !! keeps a rank's fine memory ~1/num_procs of the pool.
    integer              :: max_f1, max_f2, max_f3
    integer              :: mbuf1_lo, mbuf1_hi, mbuf2_lo, mbuf2_hi, mbuf3_lo, mbuf3_hi
    logical, allocatable :: amr_slot_live(:)
    !! fine-fine seam pack buffers, hoisted out of the per-seam s_amr_fine_fine_halo loop (allocated once at max seam extent)
    real(wp), allocatable :: amr_seambuf_x(:), amr_seambuf_y(:)
    !! cached same-level adjacent-seam list (3, npairs) = (xb, yb, seam-dim), so s_amr_fine_fine_halo iterates O(#seams) instead of
    !! rescanning all O(nblocks^2) pairs every RK stage. Block topology changes only at regrid/restart, so the list is rebuilt only
    !! when amr_seam_pairs_dirty is set (or the block count changes - a tripwire).
    integer, allocatable :: amr_seam_pairs(:,:)
    integer              :: amr_num_seam_pairs, amr_seam_pairs_nblk
    logical              :: amr_seam_pairs_dirty
    !! cached per-block P2P overlap-rank lists (rebuilt with the seam list - same dirty flag): amr_ovl_gather(:,k) = ranks whose
    !! owned coarse range (s_amr_rank_coarse_range) intersects block k's amr_cpat_mar-padded patch box (gather contributors);
    !! amr_ovl_scatter(:,k) = ranks whose coarse interior (s_amr_rank_interior) intersects block k's region box (restrict-scatter
    !! destinations). Rank-ASCENDING and NOT owner-excluded (consumers keep their owner skip), so iterating a list reproduces the
    !! replaced per-call 0..num_procs-1 scan's MPI send/recv order exactly.
    integer, allocatable :: amr_ovl_gather(:,:), amr_ovl_scatter(:,:)  !< (num_procs, amr_max_blocks)
    integer, allocatable :: amr_ovl_gather_n(:), amr_ovl_scatter_n(:)  !< per-block list lengths

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
    !! .claude/rules/common-pitfalls.md). amr_swapped guards against a nested/unpaired swap.
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

    logical :: amr_swapped = .false.  !< paired-swap guard: a nested swap would corrupt the bounce buffers silently
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

    !> Per-block gathered coarse patch (fine-level distribution). The block owner may not hold the coarse cells its block refines,
    !! so before each prolongation/ghost-fill the coarse patch spanning region_lo-amr_cpat_mar : region_hi+amr_cpat_mar (the full
    !! coarse-cell reach of every prolongation stencil) is gathered here POINT-TO-POINT from the coarse-owners
    !! (s_amr_gather_coarse_patch). Stored in amr_cg as stp scalar_fields (a drop-in for the coarse q_cons in the prolong/ghost-fill
    !! kernels) in a block-LOCAL frame: amr_cg cell 0 is GLOBAL coarse cell amr_cpat_off(d). Messages carry wp, cast to stp
    !! (identity for stp coarse), so at np=1 (owner copies its own coarse) the patch equals the local coarse read bit-for-bit. Sized
    !! to the largest block.
    type(scalar_field), allocatable :: amr_cg(:)
    integer                         :: amr_cpat_mar = 0     !< coarse-cell stencil reach = (buff_size+1)/2 + 1 (matches nmar)
    integer                         :: amr_cpat_hi(3) = 0   !< amr_cg upper local bounds per dim (0 in collapsed dims)
    integer                         :: amr_cpat_off(3) = 0  !< GLOBAL coarse index of amr_cg local cell 0 (region_lo - amr_cpat_mar)
    !> Gathered coarse pb/mv patch for non-polytropic QBMM (analogue of amr_cg): the block's coarse-side pb/mv side-state,
    !! P2P-gathered from the coarse-cell owners into the block owner in the amr_cg patch-local frame (cell 0 == amr_cpat_off). Read
    !! by the pb/mv prolong + ghost-fill so np>=2 couples to the correct coarse rank. Allocated only for non-polytropic QBMM.
    real(stp), allocatable, dimension(:,:,:,:,:) :: amr_cg_pb, amr_cg_mv
    $:GPU_DECLARE(create='[amr_cg_pb, amr_cg_mv]')

    !> Replicated coarse-decomposition table (fine-level distribution, point-to-point coupling): amr_decomp(1:3, r) = rank r's
    !! coarse start_idx per dim, amr_decomp(4:6, r) = its coarse extent (m/n/p). Allgathered once at init (coarse decomposition
    !! fixed for the run). Lets any rank compute which ranks own a given global coarse-cell range, so gather/scatter/reflux coupling
    !! is owner <-> only the (SFC-local) coarse-owners - no global collective per block per stage.
    integer, allocatable :: amr_decomp(:,:)  !< (1:6, 0:num_procs-1)

contains

    !> Build the static refined level-1 block. No-op unless amr. Called after the level-0 grid (x_cb/dx ready) and time-steppers
    !! (sys_size/buff_size set). Per-slot fine arrays allocated lazily (s_amr_reconcile_slots) - only the blocks a rank owns.
    impure subroutine s_initialize_amr_module()

        integer                         :: i, d, islot, tcap
        integer                         :: sidx(3), ext(3), maxc_loc(3), bad_loc, bad_glb, fit_d
        integer                         :: blk_lo(3), blk_hi(3)
        type(scalar_field), allocatable :: tmp_cg(:)

        if (.not. amr) return

        amr_dt_fine = 0.5_wp*dt

        ! fixed pool of amr_max_blocks slots; init activates exactly one (amr_cur = 1); regrid clusters into up to amr_max_blocks
        allocate (amr_slots(1:amr_max_blocks))
        allocate (amr_region_lo_all(3, amr_max_blocks), amr_region_hi_all(3, amr_max_blocks))
        allocate (amr_isect_lo_all(3, amr_max_blocks), amr_isect_hi_all(3, amr_max_blocks))
        allocate (amr_owns_all(amr_max_blocks))
        allocate (amr_block_owner(amr_max_blocks))
        allocate (amr_block_level(amr_max_blocks))
        allocate (amr_ovl_gather(num_procs, amr_max_blocks), amr_ovl_gather_n(amr_max_blocks))
        allocate (amr_ovl_scatter(num_procs, amr_max_blocks), amr_ovl_scatter_n(amr_max_blocks))
        amr_region_lo_all = 0; amr_region_hi_all = 0; amr_isect_lo_all = 0; amr_isect_hi_all = 0; amr_owns_all = .false.
        amr_block_owner = 0
        amr_block_level = 1  ! init default (level-1); regrid re-tags each block's level for nesting
        amr_num_levels = 1
        amr_num_blocks = 1
        amr_cur = 1

        ! fine-level load balance is capped at min(num blocks, amr_max_blocks) ranks: the SFC map spreads whole blocks, so with
        ! fewer
        ! blocks than ranks some ranks own no fine work. Warn when the pool itself is the limit (raise amr_max_blocks).
        if (proc_rank == 0 .and. num_procs > amr_max_blocks) then
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

        ! regrid size cap: min over ranks of the local half-extent (see the declaration; = amr_maxc at np=1), so any clamped box
        ! satisfies every rank's scratch constraint and can move freely across ranks
        amr_maxc_fit = amr_maxc
        do d = 1, num_dims
            call s_mpi_allreduce_integer_min((ext(d) + 1)/amr_ref_ratio, fit_d)
            amr_maxc_fit(d) = min(amr_maxc(d), fit_d)
        end do

        ! preallocation cap for MY fine arrays: a block is owned WHOLE, so any rank must hold an entire block - but the
        ! scratch-constraint abort above (allreduced over every rank's local half-extent) guarantees no block exceeds amr_maxc_fit =
        ! min-over-ranks local-half, and regrid clamps boxes to it. So amr_maxc_fit (NOT the global-half amr_maxc) is the true max
        ! block a rank can own; sizing to it right-sizes the fine/coord arrays (~1/num_procs the memory of global-half at scale). At
        ! np=1 amr_maxc_fit == amr_maxc, so the sizing (and everything) is unchanged.
        maxc_loc = amr_maxc_fit

        ! max fine extents and buffered bounds for preallocation
        max_f1 = amr_ref_ratio*maxc_loc(1) - 1
        max_f2 = 0; max_f3 = 0
        if (n_glb > 0) max_f2 = amr_ref_ratio*maxc_loc(2) - 1
        if (p_glb > 0) max_f3 = amr_ref_ratio*maxc_loc(3) - 1

        ! fine-fine seam pack buffers: sized ONCE to the largest possible seam (sys_size*buff_size * max transverse fine face), so
        ! s_amr_fine_fine_halo reuses them instead of allocating per seam per stage. Transverse cap = the largest fine face over the
        ! choice of seam dimension: 3D -> max pairwise product, 2D -> the larger single dim, 1D -> 1.
        tcap = 1
        if (n_glb > 0 .and. p_glb > 0) then
            tcap = max((max_f1 + 1)*(max_f2 + 1), (max_f2 + 1)*(max_f3 + 1), (max_f1 + 1)*(max_f3 + 1))
        else if (n_glb > 0) then
            tcap = max(max_f1, max_f2) + 1
        end if
        allocate (amr_seambuf_x(sys_size*buff_size*tcap), amr_seambuf_y(sys_size*buff_size*tcap))
        amr_seam_pairs_dirty = .true.; amr_seam_pairs_nblk = -1  ! force a seam-list build on the first fine-fine halo
        mbuf1_lo = -buff_size; mbuf1_hi = max_f1 + buff_size
        mbuf2_lo = 0; mbuf2_hi = 0; mbuf3_lo = 0; mbuf3_hi = 0
        if (n_glb > 0) then; mbuf2_lo = -buff_size; mbuf2_hi = max_f2 + buff_size; end if
        if (p_glb > 0) then; mbuf3_lo = -buff_size; mbuf3_hi = max_f3 + buff_size; end if

        ! bounce buffers for copy-based coord swap (GPU-safe; same bounds as the base-level global arrays)
        allocate (sw_x_cb(-1 - buff_size:m + buff_size))
        allocate (sw_x_cc(-buff_size:m + buff_size))
        allocate (sw_dx(-buff_size:m + buff_size))
        if (n_glb > 0) then
            allocate (sw_y_cb(-1 - buff_size:n + buff_size))
            allocate (sw_y_cc(-buff_size:n + buff_size))
            allocate (sw_dy(-buff_size:n + buff_size))
        end if
        if (p_glb > 0) then
            allocate (sw_z_cb(-1 - buff_size:p + buff_size))
            allocate (sw_z_cc(-buff_size:p + buff_size))
            allocate (sw_dz(-buff_size:p + buff_size))
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

        ! replicated coarse-decomposition table for point-to-point coupling (see decl). Coarse decomposition is fixed for the run,
        ! so allgathered once. Row r = [start_idx(1:3), m, n, p] for rank r (collapsed dims pinned to 0).
        allocate (amr_decomp(6,0:num_procs - 1))
        block
            integer :: myrow(6), ierr
            myrow = 0
            myrow(1) = start_idx(1); myrow(4) = m
            if (n_glb > 0) then; myrow(2) = start_idx(2); myrow(5) = n; end if
            if (p_glb > 0) then; myrow(3) = start_idx(3); myrow(6) = p; end if
#ifdef MFC_MPI
            call MPI_ALLGATHER(myrow, 6, MPI_INTEGER, amr_decomp, 6, MPI_INTEGER, MPI_COMM_WORLD, ierr)
#else
            amr_decomp(:,0) = myrow
#endif
        end block

        ! per-slot fine-grid IB marker fields (static-body AMR); sized to the same max buffered fine extents as q_cons so the fine
        ! IB pipeline can resolve the body on the block
        if (ib) call s_ibm_alloc_fine(amr_max_blocks, mbuf1_lo, mbuf1_hi, mbuf2_lo, mbuf2_hi, mbuf3_lo, mbuf3_hi)

        ! set geometry (region, m/n/p, idwbuff, coordinates) for the initial block (amr_cur = 1). Under dynamic regrid with bodies
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
                call s_amr_tile_box(blk_lo, blk_hi, tiled, nt, amr_max_blocks, capt)
            end if
            amr_num_blocks = nt
            ! set block regions FIRST so the owner assignment (reads amr_region_*_all) runs BEFORE the owner-dependent geometry -
            ! else s_set_amr_fine_geometry would size the whole-block owner from a stale (default) amr_block_owner
            do kk = 1, nt
                amr_region_lo_all(:,kk) = tiled(kk)%lo; amr_region_hi_all(:,kk) = tiled(kk)%hi
            end do
            call s_amr_assign_block_owners()  ! assign each block's single owner rank (fine-dist map)
            call s_amr_reconcile_slots()  ! allocate this rank's owned initial blocks (owner-guarded geometry writes below)
            do kk = 1, nt
                amr_cur = kk
                call s_set_amr_fine_geometry(tiled(kk)%lo, tiled(kk)%hi)
            end do
            call s_amr_select_slot(1)  ! refresh the per-block mirrors for slot 1 (geometry loop left them on the last tile)
            deallocate (tiled)
        end block

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
    !! patch intersected with its contiguous owned coarse range (s_amr_rank_coarse_range, = the f_amr_own_coarse set), from the
    !! replicated amr_decomp table. Non-participants send/recv nothing (no global collective). At np=1 the owner just copies its own
    !! coarse over the patch, bit-for-bit. Runtime (pull_host) packs/unpacks the overlap boxes on the DEVICE (q_coarse
    !! device-current with valid ghosts); init/regrid fills from the host (host-current with valid ghosts). Packed data is wp, cast
    !! to stp into amr_cg (identity for stp coarse), device-current on exit. INVARIANT: "coarse" here means the block's PARENT level
    !! (level l-1), NOT the base grid (level 0). For a level-1 block the parent IS L0, but a level>=2 block folds to/from its parent
    !! block's fine array; the C<->F prolong/restrict/gather routines all operate in the parent-fine frame, not the L0 frame. TWIN
    !! s_amr_gather_coarse_patch_pbmv (q<->pb/mv): same P2P skeleton (rank-range, intersection, pack/send/recv/unpack) and
    !! patch-local frame - keep lockstep.
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
            call s_amr_gather_from_parent(pull_host)
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
        if (amr_seam_pairs_dirty .or. amr_seam_pairs_nblk /= amr_num_blocks) call s_amr_build_seam_pairs()

        if (proc_rank == owner) then
            ! fill the cells this rank holds locally (own box), then receive the rest from the other coarse-owners
            call s_amr_rank_coarse_range(proc_rank, crlo, crhi)
            call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
            if (pull_host) then
                ! runtime: q_coarse is device-current - copy the own box on the device (same index map/assignment as the host path)
                call s_amr_gather_own_box_device(q_coarse, bl, bh, o1, o2, o3)
            else
                call s_amr_unpack_patch(q_coarse, bl, bh, o1, o2, o3)  ! local read: q_coarse own frame -> amr_cg patch frame
            end if
            ! count + post recvs from every OTHER rank whose owned range overlaps the patch (cached list; every listed rank
            ! overlaps by construction)
            nsrc = 0
            do idx = 1, amr_ovl_gather_n(amr_cur)
                if (amr_ovl_gather(idx, amr_cur) /= owner) nsrc = nsrc + 1
            end do
            if (nsrc > 0) then
                allocate (rbuf(maxsz, nsrc), reqs(nsrc), srank(nsrc))
                nsrc = 0
                do idx = 1, amr_ovl_gather_n(amr_cur)
                    r = amr_ovl_gather(idx, amr_cur)
                    if (r == owner) cycle
                    call s_amr_rank_coarse_range(r, crlo, crhi)
                    call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
                    boxsz = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                    nsrc = nsrc + 1; srank(nsrc) = r
#ifdef MFC_MPI
                    call MPI_IRECV(rbuf(1, nsrc), boxsz, mpi_p, r, amr_cur, MPI_COMM_WORLD, reqs(nsrc), ierr)
#endif
                end do
#ifdef MFC_MPI
                call MPI_WAITALL(nsrc, reqs, MPI_STATUSES_IGNORE, ierr)
#endif
                do idx = 1, nsrc
                    call s_amr_rank_coarse_range(srank(idx), crlo, crhi)
                    call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
                    if (pull_host) then
                        ! runtime: unpack ONLY this box's wire buffer on the device (same order/cast as the host unpack below)
                        boxsz = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                        call s_amr_unpack_box_device(bl, bh, rbuf(1:boxsz,idx))
                        cycle
                    end if
                    ! unpack in the SAME (i, g3, g2, g1) order the sender packed; place at amr_cg patch-local index
                    r = 0
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
                deallocate (rbuf, reqs, srank)
            end if
            ! host path only: the runtime device path wrote amr_cg on the device directly (host amr_cg stays stale, as at np=1 -
            ! runtime consumers read the device copy)
            if (.not. pull_host) then
                do i = 1, sys_size
                    $:GPU_UPDATE(device='[amr_cg(i)%sf]')
                end do
            end if
        else
            ! non-owner: if my owned coarse range overlaps the patch, pack my slice (wp) and send it to the owner
            call s_amr_rank_coarse_range(proc_rank, crlo, crhi)
            call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
            if (bl(1) <= bh(1) .and. bl(2) <= bh(2) .and. bl(3) <= bh(3)) then
                boxsz = sys_size*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                allocate (sbuf(boxsz))
                if (pull_host) then
                    ! runtime: pack the overlap box on the device straight into sbuf (only the box crosses PCIe)
                    call s_amr_pack_box_device(q_coarse, bl, bh, o1, o2, o3, sbuf)
                else
                    idx = 0
                    do i = 1, sys_size
                        do g3 = bl(3), bh(3)
                            do g2 = bl(2), bh(2)
                                do g1 = bl(1), bh(1)
                                    idx = idx + 1; sbuf(idx) = real(q_coarse(i)%sf(g1 - o1, g2 - o2, g3 - o3), wp)
                                end do
                            end do
                        end do
                    end do
                end if
#ifdef MFC_MPI
                call MPI_SEND(sbuf, boxsz, mpi_p, owner, amr_cur, MPI_COMM_WORLD, ierr)
#endif
                deallocate (sbuf)
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
                allocate (rbuf(maxsz, nsrc), reqs(nsrc), srank(nsrc))
                nsrc = 0
                do idx = 1, amr_ovl_gather_n(amr_cur)
                    r = amr_ovl_gather(idx, amr_cur)
                    if (r == owner) cycle
                    call s_amr_rank_coarse_range(r, crlo, crhi)
                    call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
                    boxsz = cellsz*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                    nsrc = nsrc + 1; srank(nsrc) = r
#ifdef MFC_MPI
                    call MPI_IRECV(rbuf(1, nsrc), boxsz, mpi_p, r, amr_cur, MPI_COMM_WORLD, reqs(nsrc), ierr)
#endif
                end do
#ifdef MFC_MPI
                call MPI_WAITALL(nsrc, reqs, MPI_STATUSES_IGNORE, ierr)
#endif
                do idx = 1, nsrc
                    call s_amr_rank_coarse_range(srank(idx), crlo, crhi)
                    call s_amr_box_isect(plo, phi, crlo, crhi, bl, bh)
                    if (pull_host) then
                        ! runtime: unpack ONLY this box's wire buffer on the device (same order/cast as the host unpack below)
                        boxsz = cellsz*(bh(1) - bl(1) + 1)*(bh(2) - bl(2) + 1)*(bh(3) - bl(3) + 1)
                        call s_amr_unpack_box_pbmv_device(bl, bh, rbuf(1:boxsz,idx))
                        cycle
                    end if
                    ! unpack in the SAME (ib_, q, g3, g2, g1) order the sender packed - pb block then mv block
                    r = 0
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
                allocate (sbuf(boxsz))
                if (pull_host) then
                    ! runtime: pack the overlap box on the device straight into sbuf (only the box crosses PCIe)
                    call s_amr_pack_box_pbmv_device(pb_coarse, mv_coarse, bl, bh, o1, o2, o3, sbuf)
                else
                    idx = 0
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
                call MPI_SEND(sbuf, boxsz, mpi_p, owner, amr_cur, MPI_COMM_WORLD, ierr)
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
        ! lock-step fill: gather from the parent's CURRENT fine state. pull_host stays in the signature for the level-1 path.
        ! Owner-guard at the CALL SITE: on a non-owner rank the parent slot is unallocated (co-located tower - owner holds both
        ! block
        ! and parent), and passing amr_slots(pblk)%q_cons would dereference it before the callee's internal early-return. to_host =
        ! .not. pull_host: init/regrid (pull_host=F) feed the host prolong/self-test; runtime (pull_host=T) reads amr_cg on the
        ! device in the C/F ghost-fill, so skip the device->host copy.
        if (amr_rank_owns_block) call s_amr_gather_from_parent_field(pblk, amr_slots(pblk)%q_cons, .not. pull_host)

    end subroutine s_amr_gather_from_parent

    !> Gather amr_cg (the current level>=2 block's coarse patch) from a SPECIFIC parent snapshot field qp, in the parent-fine cell
    !! frame (amr_isect_lo/hi already parent-fine from s_set_amr_fine_geometry). The subcycle recursion calls this twice per parent
    !! substep - qp = the parent slot's q_cons_stor (t^n bracket) then q_cons (t^{n+1} bracket) - to build the child's two
    !! ghost-lerp sources. np=1 = a local copy on the owner (which also owns the parent); np>=2 P2P (parent owner -> block owner) is
    !! future work.
    impure subroutine s_amr_gather_from_parent_field(pblk, qp, to_host)

        integer, intent(in)                                 :: pblk
        type(scalar_field), dimension(sys_size), intent(in) :: qp
        logical, intent(in)                                 :: to_host  !< host copy of amr_cg needed (init/regrid), not runtime
        integer                                             :: w1, w2, w3

        amr_cpat_off = 0
        amr_cpat_off(1) = amr_isect_lo(1) - amr_cpat_mar
        if (n_glb > 0) amr_cpat_off(2) = amr_isect_lo(2) - amr_cpat_mar
        if (p_glb > 0) amr_cpat_off(3) = amr_isect_lo(3) - amr_cpat_mar
        w1 = (amr_isect_hi(1) - amr_isect_lo(1)) + 2*amr_cpat_mar
        w2 = 0; w3 = 0
        if (n_glb > 0) w2 = (amr_isect_hi(2) - amr_isect_lo(2)) + 2*amr_cpat_mar
        if (p_glb > 0) w3 = (amr_isect_hi(3) - amr_isect_lo(3)) + 2*amr_cpat_mar
        if (.not. amr_rank_owns_block) return  ! np=1: the owner holds both this block and its parent
        ! copy the parent's fine patch into amr_cg with a DEVICE kernel (qp passed as an argument, present-table safe like
        ! s_amr_restrict_overwrite_device). np>=2 P2P is future work.
        call s_amr_copy_parent_patch(qp, w1, w2, w3, to_host)

    end subroutine s_amr_gather_from_parent_field

    !> Device kernel for s_amr_gather_from_parent: copy the parent block's fine patch into amr_cg over [amr_cpat_off : + w]. The
    !! parent q_cons is passed as the qp ARGUMENT (not indexed as amr_slots(pblk) inside the kernel) so its deep %sf attach resolves
    !! present-table safe, like s_amr_restrict_overwrite_device. amr_cg is then synced to host for host consumers (init self-test's
    !! restrict-prolong check).
    impure subroutine s_amr_copy_parent_patch(qp, w1, w2, w3, to_host)

        type(scalar_field), dimension(sys_size), intent(in) :: qp
        integer, intent(in)                                 :: w1, w2, w3
        !> .true. only for the init/regrid HOST consumers (whole-block host prolong + restrict-prolong self-test). The runtime C/F
        !! ghost-fill reads amr_cg on the DEVICE (filled by the kernel below), so no device->host copy is needed.
        logical, intent(in) :: to_host
        integer             :: i, g1, g2, g3, o1, o2, o3

        o1 = amr_cpat_off(1); o2 = amr_cpat_off(2); o3 = amr_cpat_off(3)
        $:GPU_PARALLEL_LOOP(collapse=4)
        do i = 1, sys_size
            do g3 = 0, w3
                do g2 = 0, w2
                    do g1 = 0, w1
                        amr_cg(i)%sf(g1, g2, g3) = qp(i)%sf(g1 + o1, g2 + o2, g3 + o3)
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

    end subroutine s_amr_copy_parent_patch

    !> Rank r's contiguous owned coarse-cell range per dim from the replicated amr_decomp table: interior [start:start+ext] plus its
    !! physical-boundary ghosts (buff_size cells only where the subdomain touches the domain edge). Equal to the set where
    !! f_amr_own_coarse is true, but as one contiguous span so box intersections identify contributors without a per-cell scan.
    pure subroutine s_amr_rank_coarse_range(r, crlo, crhi)

        integer, intent(in)  :: r
        integer, intent(out) :: crlo(3), crhi(3)

        crlo = 0; crhi = 0
        crlo(1) = amr_decomp(1, r); if (amr_decomp(1, r) == 0) crlo(1) = -buff_size
        crhi(1) = amr_decomp(1, r) + amr_decomp(4, r); if (crhi(1) == m_glb) crhi(1) = crhi(1) + buff_size
        if (n_glb > 0) then
            crlo(2) = amr_decomp(2, r); if (amr_decomp(2, r) == 0) crlo(2) = -buff_size
            crhi(2) = amr_decomp(2, r) + amr_decomp(5, r); if (crhi(2) == n_glb) crhi(2) = crhi(2) + buff_size
        end if
        if (p_glb > 0) then
            crlo(3) = amr_decomp(3, r); if (amr_decomp(3, r) == 0) crlo(3) = -buff_size
            crhi(3) = amr_decomp(3, r) + amr_decomp(6, r); if (crhi(3) == p_glb) crhi(3) = crhi(3) + buff_size
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
    !! replicated amr_decomp table (so the block owner can decide which ranks to send freg to, and each rank agrees on whether it
    !! receives). Uses amr_region_lo/hi (the current block, set on every rank by s_amr_select_slot).
    pure logical function f_amr_reflux_participates(r) result(part)

        integer, intent(in) :: r
        integer             :: sidx(3), ext(3), d, t
        logical             :: tv(3), tvd

        sidx = 0; ext = 0
        sidx(1) = amr_decomp(1, r); ext(1) = amr_decomp(4, r)
        if (n_glb > 0) then; sidx(2) = amr_decomp(2, r); ext(2) = amr_decomp(5, r); end if
        if (p_glb > 0) then; sidx(3) = amr_decomp(3, r); ext(3) = amr_decomp(6, r); end if
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

    !> Fine-level distribution: deliver the current block's fine flux registers freg (captured by the owner during the fine advance)
    !! to exactly the (SFC-local) coarse-outside-owners that apply the reflux - POINT-TO-POINT, replacing the global broadcast. The
    !! owner sends its whole freg slot (block-relative; each applier reads its own transverse slice) to every participant; non-owner
    !! participants receive it. Device-resident: owner stages its slot to host, receivers push it back. No-op without MPI/at np=1.
    impure subroutine s_amr_p2p_reflux_faces()

#ifdef MFC_MPI
        integer              :: owner, r, ierr, nreq, cnt
        integer, allocatable :: reqs(:)

        if (.not. amr) return
        if (num_procs == 1) return
        owner = amr_block_owner(amr_cur)
        if (proc_rank == owner) then
            #:for D in [1, 2, 3]
                if (${D}$ <= num_dims) then
                    $:GPU_UPDATE(host='[freg(' + str(D) + ')%lo(:, :, :, amr_cur), freg(' + str(D) + ')%hi(:, :, :, amr_cur)]')
                end if
            #:endfor
            nreq = 0
            do r = 0, num_procs - 1
                if (r /= owner .and. f_amr_reflux_participates(r)) nreq = nreq + 1
            end do
            if (nreq > 0) then
                allocate (reqs(2*num_dims*nreq))
                nreq = 0
                do r = 0, num_procs - 1
                    if (r == owner .or. .not. f_amr_reflux_participates(r)) cycle
                    #:for D in [1, 2, 3]
                        if (${D}$ <= num_dims) then
                            cnt = size(freg(${D}$)%lo(:,:,:,amr_cur))
                            nreq = nreq + 1
                            call MPI_ISEND(freg(${D}$)%lo(:,:,:,amr_cur), cnt, mpi_p, r, ${2*D}$, MPI_COMM_WORLD, reqs(nreq), ierr)
                            nreq = nreq + 1
                            call MPI_ISEND(freg(${D}$)%hi(:,:,:,amr_cur), cnt, mpi_p, r, ${2*D + 1}$, MPI_COMM_WORLD, reqs(nreq), &
                                           & ierr)
                        end if
                    #:endfor
                end do
                call MPI_WAITALL(nreq, reqs, MPI_STATUSES_IGNORE, ierr)
                deallocate (reqs)
            end if
        else if (f_amr_reflux_participates(proc_rank)) then
            #:for D in [1, 2, 3]
                if (${D}$ <= num_dims) then
                    cnt = size(freg(${D}$)%lo(:,:,:,amr_cur))
                    call MPI_RECV(freg(${D}$)%lo(:,:,:,amr_cur), cnt, mpi_p, owner, ${2*D}$, MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                                  & ierr)
                    call MPI_RECV(freg(${D}$)%hi(:,:,:,amr_cur), cnt, mpi_p, owner, ${2*D + 1}$, MPI_COMM_WORLD, &
                                  & MPI_STATUS_IGNORE, ierr)
                    $:GPU_UPDATE(device='[freg(' + str(D) + ')%lo(:, :, :, amr_cur), freg(' + str(D) + ')%hi(:, :, :, amr_cur)]')
                end if
            #:endfor
        end if
#endif

    end subroutine s_amr_p2p_reflux_faces

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
        call MPI_ALLREDUCE(MPI_IN_PLACE, cost, amr_num_blocks, mpi_p, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

    end subroutine s_amr_block_cost

    !> Fine-level distribution map: assigns each active block a single owner rank by chains-on-chains balancing of fine-work weight
    !! in Morton order of the block's low corner - the same SFC idea m_sfc_partition uses for the base grid, at block granularity.
    !! Fine work = measured coarse-footprint cost (s_amr_block_cost) x the level's refinement factor per active dim, so blocks
    !! concentrating IB or phase-change work weigh more than equal-size quiescent ones. The cost vector is allreduced (one
    !! collective; every rank must call this), after which the assignment is deterministic and identical on every rank.
    !! s_set_amr_fine_geometry applies it as amr_rank_owns_block = (amr_block_owner(amr_cur) == proc_rank).
    impure subroutine s_amr_assign_block_owners()

        integer         :: k, kk, r, a, lev, maxlev, ord(amr_num_blocks)
        integer(kind=8) :: key(amr_num_blocks), tmpk
        real(wp)        :: wt(amr_num_blocks), twt(amr_num_blocks), cost(amr_num_blocks), cum, tgt, total
        integer         :: tmpo

        if (amr_num_blocks < 1) return

        call s_amr_block_cost(cost)

        ! per-block own fine-work weight = footprint cost x amr_ref_ratio**(level*active dims). A level-l block is amr_ref_ratio**l
        ! finer than L0 per dim, so its work is the footprint cost x rr**(l*d); the level factor keeps the co-located-tower load
        ! balance honest. With no cost signals this reduces to the fine cell count (geometry only).
        do k = 1, amr_num_blocks
            wt(k) = cost(k)*real(amr_ref_ratio, wp)**amr_block_level(k)
            if (n_glb > 0) wt(k) = wt(k)*real(amr_ref_ratio, wp)**amr_block_level(k)
            if (p_glb > 0) wt(k) = wt(k)*real(amr_ref_ratio, wp)**amr_block_level(k)
            key(k) = f_morton(amr_region_lo_all(1, k), amr_region_lo_all(2, k), amr_region_lo_all(3, k))
            ord(k) = k
        end do

        ! CO-LOCATE refinement towers: a level-1 block and all its nested descendants are owned WHOLE by one rank so every
        ! parent<->child gather/restrict/reflux stays LOCAL (no new MPI - only L0<->L1 crosses ranks). Roll each block's own work up
        ! onto its top-level (level-1) ancestor, SFC-balance only the level-1 anchors, then let descendants inherit. Single-level
        ! (all blocks level 1): twt == wt and the inherit pass is empty, so this is byte-identical to before.
        twt = 0._wp
        do k = 1, amr_num_blocks
            a = k
            do while (amr_block_level(a) > 1)
                if (f_amr_parent_block(a) < 1) exit  ! proper-nesting invariant broken; stop before indexing amr_block_level(0)
                a = f_amr_parent_block(a)
            end do
            twt(a) = twt(a) + wt(k)
        end do

        ! sort block indices by Morton key (insertion sort - amr_num_blocks is small, <= amr_max_blocks)
        do k = 2, amr_num_blocks
            tmpk = key(ord(k)); tmpo = ord(k); kk = k - 1
            do while (kk >= 1)
                if (key(ord(kk)) <= tmpk) exit
                ord(kk + 1) = ord(kk); kk = kk - 1
            end do
            ord(kk + 1) = tmpo
        end do

        ! chains-on-chains over the level-1 anchors in SFC order, weighted by whole-tower work; advance the owner rank when the
        ! cumulative tower weight crosses the next even share. All-real arithmetic on the replicated (allreduced) weights, fixed
        ! loop
        ! order, so every rank computes the identical assignment.
        total = 0._wp
        do k = 1, amr_num_blocks
            if (amr_block_level(k) == 1) total = total + twt(k)
        end do
        r = 0; cum = 0._wp
        do k = 1, amr_num_blocks
            if (amr_block_level(ord(k)) /= 1) cycle
            tgt = real(r + 1, wp)*total/real(num_procs, wp)
            if (cum >= tgt .and. r < num_procs - 1) r = r + 1
            amr_block_owner(ord(k)) = r
            cum = cum + twt(ord(k))
        end do

        ! descendants inherit their parent's owner (top-down, level by level - parents already assigned when their children run)
        maxlev = maxval(amr_block_level(1:amr_num_blocks))
        do lev = 2, maxlev
            do k = 1, amr_num_blocks
                if (amr_block_level(k) == lev) amr_block_owner(k) = amr_block_owner(f_amr_parent_block(k))
            end do
        end do

    end subroutine s_amr_assign_block_owners

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
        integer             :: sidx(3), ext(3), nmar, bad_loc, bad_glb, pblk, d, rr

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
                pblk = f_amr_parent_block(amr_cur); rr = amr_slots(pblk)%amr_ref_ratio
                do d = 1, 3
                    amr_isect_lo(d) = rr*(lo(d) - amr_region_lo_all(d, pblk))
                    amr_isect_hi(d) = rr*(hi(d) - amr_region_lo_all(d, pblk)) + (rr - 1)
                end do
                if (n_glb == 0) then; amr_isect_lo(2) = 0; amr_isect_hi(2) = 0; end if
                if (p_glb == 0) then; amr_isect_lo(3) = 0; amr_isect_hi(3) = 0; end if
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
        ! coord building only on ranks with fine cells (others never read their coord arrays); the parent origin is this rank's
        ! INTERSECTION start, converted to LOCAL indexing so the bisection reads its x_cb slice
        if (amr_rank_owns_block .and. amr_block_level(amr_cur) >= 2) then
            ! level >= 2: bisect the PARENT block's fine coords (its x_cb, lbound -1), parent origin = the parent-fine footprint
            call s_build_level_coords(amr_slots(pblk)%x_cb, -1, amr_isect_lo(1), amr_slots(amr_cur)%m, amr_slots(amr_cur)%x_cb, &
                                      & amr_slots(amr_cur)%x_cc, amr_slots(amr_cur)%dx)
            if (n_glb > 0) call s_build_level_coords(amr_slots(pblk)%y_cb, -1, amr_isect_lo(2), amr_slots(amr_cur)%n, &
                & amr_slots(amr_cur)%y_cb, amr_slots(amr_cur)%y_cc, amr_slots(amr_cur)%dy)
            if (p_glb > 0) call s_build_level_coords(amr_slots(pblk)%z_cb, -1, amr_isect_lo(3), amr_slots(amr_cur)%p, &
                & amr_slots(amr_cur)%z_cb, amr_slots(amr_cur)%z_cc, amr_slots(amr_cur)%dz)
        else if (amr_rank_owns_block) then
            ! level 1: whole-block fine coords from the GLOBAL L0 boundaries (owner may not hold the coarse coordinate slice for
            ! cells
            ! it now refines). amr_gxcb has lbound -1; the parent origin is the block's GLOBAL low corner.
            call s_build_level_coords(amr_gxcb, -1, amr_isect_lo(1), amr_slots(amr_cur)%m, amr_slots(amr_cur)%x_cb, &
                                      & amr_slots(amr_cur)%x_cc, amr_slots(amr_cur)%dx)
            if (n_glb > 0) call s_build_level_coords(amr_gycb, -1, amr_isect_lo(2), amr_slots(amr_cur)%n, &
                & amr_slots(amr_cur)%y_cb, amr_slots(amr_cur)%y_cc, amr_slots(amr_cur)%dy)
            if (p_glb > 0) call s_build_level_coords(amr_gzcb, -1, amr_isect_lo(3), amr_slots(amr_cur)%p, &
                & amr_slots(amr_cur)%z_cb, amr_slots(amr_cur)%z_cc, amr_slots(amr_cur)%dz)
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
        call s_mpi_allreduce_integer_max(bad_loc, bad_glb)
        amr_xchg_coarse_ghosts = bad_glb == 1

    end subroutine s_set_amr_fine_geometry

    !> Conservative-linear prolongation for a single variable pair. Reads coarse interior/ghost from qc; writes fine interior to qf.
    !! Minmod-limited slopes.
    impure subroutine s_prolong_one_var(qc, qf, pos, inject)

        type(scalar_field), intent(in)    :: qc
        type(scalar_field), intent(inout) :: qf
        logical, optional, intent(in)     :: pos     !< floor the child at bub_pos_frac*u0 (bubble radius-moment realizability)
        logical, optional, intent(in)     :: inject  !< piecewise-constant (child = u0): QBMM moment realizability preservation
        integer                           :: fi, fj, fk, ci, cj, ck, ox, oy, oz
        real(wp)                          :: u0, sx, sy, sz, xix, xiy, xiz, child
        logical                           :: floor_pos, pw_const

        floor_pos = .false.; if (present(pos)) floor_pos = pos
        pw_const = .false.; if (present(inject)) pw_const = inject

        ! coarse source qc is the gathered block-local patch amr_cg (fine-level distribution): amr_isect_lo is GLOBAL and equals
        ! region_lo on the owner, so amr_isect_lo + f/rr - amr_cpat_off = nmar + f/rr is the patch-local coarse index
        ox = amr_cpat_off(1); oy = amr_cpat_off(2); oz = amr_cpat_off(3)
        do fk = 0, amr_slots(amr_cur)%p
            ck = amr_isect_lo(3) + fk/amr_slots(amr_cur)%amr_ref_ratio - oz; if (p_glb == 0) ck = 0
            xiz = 0._wp; if (p_glb > 0) xiz = (real(mod(fk, amr_slots(amr_cur)%amr_ref_ratio), &
                             & wp) - real(amr_slots(amr_cur)%amr_ref_ratio - 1, &
                             & wp)*0.5_wp)/real(amr_slots(amr_cur)%amr_ref_ratio, wp)
            do fj = 0, amr_slots(amr_cur)%n
                cj = amr_isect_lo(2) + fj/amr_slots(amr_cur)%amr_ref_ratio - oy; if (n_glb == 0) cj = 0
                xiy = 0._wp; if (n_glb > 0) xiy = (real(mod(fj, amr_slots(amr_cur)%amr_ref_ratio), &
                                 & wp) - real(amr_slots(amr_cur)%amr_ref_ratio - 1, &
                                 & wp)*0.5_wp)/real(amr_slots(amr_cur)%amr_ref_ratio, wp)
                do fi = 0, amr_slots(amr_cur)%m
                    ci = amr_isect_lo(1) + fi/amr_slots(amr_cur)%amr_ref_ratio - ox
                    xix = (real(mod(fi, amr_slots(amr_cur)%amr_ref_ratio), wp) - real(amr_slots(amr_cur)%amr_ref_ratio - 1, &
                           & wp)*0.5_wp)/real(amr_slots(amr_cur)%amr_ref_ratio, wp)
                    u0 = real(qc%sf(ci, cj, ck), wp)
                    sx = minmod(real(qc%sf(ci + 1, cj, ck), wp) - u0, u0 - real(qc%sf(ci - 1, cj, ck), wp))
                    sy = 0._wp
                    if (n_glb > 0) sy = minmod(real(qc%sf(ci, cj + 1, ck), wp) - u0, u0 - real(qc%sf(ci, cj - 1, ck), wp))
                    sz = 0._wp
                    if (p_glb > 0) sz = minmod(real(qc%sf(ci, cj, ck + 1), wp) - u0, u0 - real(qc%sf(ci, cj, ck - 1), wp))
                    if (pw_const) then
                        sx = 0._wp; sy = 0._wp; sz = 0._wp
                    end if
                    child = u0 + sx*xix + sy*xiy + sz*xiz
                    if (floor_pos) child = max(child, bub_pos_frac*u0)
                    qf%sf(fi, fj, fk) = child
                end do
            end do
        end do

    end subroutine s_prolong_one_var

    !> Conservative-linear prolongation: fill amr_fine interior from coarse (level-0), minmod-limited. Symmetric child offsets
    !! (+/-1/4 of a coarse cell) => the amr_ref_ratio^d children average to the coarse value. Multi-fluid volume fractions take the
    !! sum-preserving closure path instead (single-fluid runs never branch, so their prolongation is untouched). TWIN
    !! s_amr_prolong_pbmv (q<->pb/mv): pb/mv sibling of this prolongation (piecewise-constant there); keep the child-offset frame
    !! and volume-fraction closure lockstep.
    impure subroutine s_interpolate_coarse_to_fine()

        integer :: i, bstride

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
            call s_prolong_one_var(amr_cg(i), amr_slots(amr_cur)%q_cons(i), &
                                   & pos=bubbles_euler .and. .not. qbmm .and. i >= eqn_idx%bub%beg .and. i <= eqn_idx%bub%end &
                                   & .and. mod(i - eqn_idx%bub%beg, bstride) /= 1, &
                                   & inject=qbmm .and. i >= eqn_idx%bub%beg .and. i <= eqn_idx%bub%end)
        end do
        if (num_fluids > 1 .and. (.not. bubbles_lagrange)) call s_prolong_alphas_closure(amr_cg, amr_slots(amr_cur)%q_cons)
        if (chemistry) call s_prolong_species_closure(amr_cg, amr_slots(amr_cur)%q_cons)

    end subroutine s_interpolate_coarse_to_fine

    !> Sum-preserving volume-fraction prolongation (num_fluids > 1): fluids adv%beg..adv%end-1 are interpolated with minmod slopes
    !! under a SHARED per-cell limiter switch (a sign change for ANY fluid in a dim zeroes that dim's slope for ALL fluids, so the
    !! closure fluid's effective slope is limited consistently) and clamped to [0,1]; the last fluid closes alpha_n = 1 -
    !! sum(others), so sum(alpha) = 1 on the fine level by construction. For two fluids the closure is also in [0,1]; for >2 fluids
    !! any residual closure undershoot is handled by mpp_lim (required by the checker). Same fine/coarse index mapping as
    !! s_prolong_one_var.
    impure subroutine s_prolong_alphas_closure(qc, qf)

        type(scalar_field), dimension(sys_size), intent(in)    :: qc
        type(scalar_field), dimension(sys_size), intent(inout) :: qf
        integer                                                :: fi, fj, fk, ci, cj, ck, ox, oy, oz, i
        real(wp)                                               :: xix, xiy, xiz, u0, sx, sy, sz, av, asum
        logical                                                :: shx, shy, shz

        ! coarse source qc is the gathered block-local patch amr_cg (fine-level distribution): patch-frame offset

        ox = amr_cpat_off(1); oy = amr_cpat_off(2); oz = amr_cpat_off(3)
        do fk = 0, amr_slots(amr_cur)%p
            ck = amr_isect_lo(3) + fk/amr_slots(amr_cur)%amr_ref_ratio - oz; if (p_glb == 0) ck = 0
            xiz = 0._wp; if (p_glb > 0) xiz = (real(mod(fk, amr_slots(amr_cur)%amr_ref_ratio), &
                             & wp) - real(amr_slots(amr_cur)%amr_ref_ratio - 1, &
                             & wp)*0.5_wp)/real(amr_slots(amr_cur)%amr_ref_ratio, wp)
            do fj = 0, amr_slots(amr_cur)%n
                cj = amr_isect_lo(2) + fj/amr_slots(amr_cur)%amr_ref_ratio - oy; if (n_glb == 0) cj = 0
                xiy = 0._wp; if (n_glb > 0) xiy = (real(mod(fj, amr_slots(amr_cur)%amr_ref_ratio), &
                                 & wp) - real(amr_slots(amr_cur)%amr_ref_ratio - 1, &
                                 & wp)*0.5_wp)/real(amr_slots(amr_cur)%amr_ref_ratio, wp)
                do fi = 0, amr_slots(amr_cur)%m
                    ci = amr_isect_lo(1) + fi/amr_slots(amr_cur)%amr_ref_ratio - ox
                    xix = (real(mod(fi, amr_slots(amr_cur)%amr_ref_ratio), wp) - real(amr_slots(amr_cur)%amr_ref_ratio - 1, &
                           & wp)*0.5_wp)/real(amr_slots(amr_cur)%amr_ref_ratio, wp)
                    call s_alpha_shared_switch(qc, ci, cj, ck, shx, shy, shz)
                    asum = 0._wp
                    do i = eqn_idx%adv%beg, eqn_idx%adv%end - 1
                        u0 = real(qc(i)%sf(ci, cj, ck), wp)
                        sx = 0._wp
                        if (shx) sx = minmod(real(qc(i)%sf(ci + 1, cj, ck), wp) - u0, u0 - real(qc(i)%sf(ci - 1, cj, ck), wp))
                        sy = 0._wp
                        if (n_glb > 0 .and. shy) sy = minmod(real(qc(i)%sf(ci, cj + 1, ck), wp) - u0, u0 - real(qc(i)%sf(ci, &
                            & cj - 1, ck), wp))
                        sz = 0._wp
                        if (p_glb > 0 .and. shz) sz = minmod(real(qc(i)%sf(ci, cj, ck + 1), wp) - u0, u0 - real(qc(i)%sf(ci, cj, &
                            & ck - 1), wp))
                        av = min(max(u0 + sx*xix + sy*xiy + sz*xiz, 0._wp), 1._wp)
                        qf(i)%sf(fi, fj, fk) = av
                        asum = asum + av
                    end do
                    qf(eqn_idx%adv%end)%sf(fi, fj, fk) = 1._wp - asum
                end do
            end do
        end do

    end subroutine s_prolong_alphas_closure

    !> Species mass-fraction prolongation closure (chemistry): each partial density rho*Y_k is minmod-prolonged and clamped
    !! non-negative, then all species are rescaled so sum_k(rho*Y_k) equals the (already prolonged) continuity density at the fine
    !! cell. This keeps the fine species realizable (Y_k >= 0, and sum(Y_k) = 1 exactly under the cons->prim recovery rho = sum
    !! rho*Y_k) and consistent with the continuity variable the reaction source reads. Same index mapping as s_prolong_one_var; cont
    !! is prolonged in the main loop before this runs.
    impure subroutine s_prolong_species_closure(qc, qf)

        type(scalar_field), dimension(sys_size), intent(in)    :: qc
        type(scalar_field), dimension(sys_size), intent(inout) :: qf
        integer                                                :: fi, fj, fk, ci, cj, ck, ox, oy, oz, i
        real(wp)                                               :: xix, xiy, xiz, u0, sx, sy, sz, av, rsum, rscale

        ! coarse source qc is the gathered block-local patch amr_cg (fine-level distribution): patch-frame offset

        ox = amr_cpat_off(1); oy = amr_cpat_off(2); oz = amr_cpat_off(3)
        do fk = 0, amr_slots(amr_cur)%p
            ck = amr_isect_lo(3) + fk/amr_slots(amr_cur)%amr_ref_ratio - oz; if (p_glb == 0) ck = 0
            xiz = 0._wp; if (p_glb > 0) xiz = (real(mod(fk, amr_slots(amr_cur)%amr_ref_ratio), &
                             & wp) - real(amr_slots(amr_cur)%amr_ref_ratio - 1, &
                             & wp)*0.5_wp)/real(amr_slots(amr_cur)%amr_ref_ratio, wp)
            do fj = 0, amr_slots(amr_cur)%n
                cj = amr_isect_lo(2) + fj/amr_slots(amr_cur)%amr_ref_ratio - oy; if (n_glb == 0) cj = 0
                xiy = 0._wp; if (n_glb > 0) xiy = (real(mod(fj, amr_slots(amr_cur)%amr_ref_ratio), &
                                 & wp) - real(amr_slots(amr_cur)%amr_ref_ratio - 1, &
                                 & wp)*0.5_wp)/real(amr_slots(amr_cur)%amr_ref_ratio, wp)
                do fi = 0, amr_slots(amr_cur)%m
                    ci = amr_isect_lo(1) + fi/amr_slots(amr_cur)%amr_ref_ratio - ox
                    xix = (real(mod(fi, amr_slots(amr_cur)%amr_ref_ratio), wp) - real(amr_slots(amr_cur)%amr_ref_ratio - 1, &
                           & wp)*0.5_wp)/real(amr_slots(amr_cur)%amr_ref_ratio, wp)
                    rsum = 0._wp
                    do i = eqn_idx%species%beg, eqn_idx%species%end
                        u0 = real(qc(i)%sf(ci, cj, ck), wp)
                        sx = minmod(real(qc(i)%sf(ci + 1, cj, ck), wp) - u0, u0 - real(qc(i)%sf(ci - 1, cj, ck), wp))
                        sy = 0._wp
                        if (n_glb > 0) sy = minmod(real(qc(i)%sf(ci, cj + 1, ck), wp) - u0, u0 - real(qc(i)%sf(ci, cj - 1, ck), wp))
                        sz = 0._wp
                        if (p_glb > 0) sz = minmod(real(qc(i)%sf(ci, cj, ck + 1), wp) - u0, u0 - real(qc(i)%sf(ci, cj, ck - 1), wp))
                        av = max(u0 + sx*xix + sy*xiy + sz*xiz, 0._wp)
                        qf(i)%sf(fi, fj, fk) = av
                        rsum = rsum + av
                    end do
                    rscale = real(qf(eqn_idx%cont%end)%sf(fi, fj, fk), wp)/max(rsum, 1.e-30_wp)
                    do i = eqn_idx%species%beg, eqn_idx%species%end
                        qf(i)%sf(fi, fj, fk) = real(qf(i)%sf(fi, fj, fk), wp)*rscale
                    end do
                end do
            end do
        end do

    end subroutine s_prolong_species_closure

    !> Shared per-cell limiter switch for the volume-fraction closure prolongation: per dim, slopes stay on only if NO fluid's
    !! centered differences change sign there (symmetric in the fluids, incl. the closure fluid).
    pure subroutine s_alpha_shared_switch(qc, ci, cj, ck, shx, shy, shz)

        type(scalar_field), dimension(sys_size), intent(in) :: qc
        integer, intent(in)                                 :: ci, cj, ck
        logical, intent(out)                                :: shx, shy, shz
        integer                                             :: i
        real(wp)                                            :: u0

        shx = .true.; shy = n_glb > 0; shz = p_glb > 0
        do i = eqn_idx%adv%beg, eqn_idx%adv%end
            u0 = real(qc(i)%sf(ci, cj, ck), wp)
            if ((real(qc(i)%sf(ci + 1, cj, ck), wp) - u0)*(u0 - real(qc(i)%sf(ci - 1, cj, ck), wp)) <= 0._wp) shx = .false.
            if (n_glb > 0) then
                if ((real(qc(i)%sf(ci, cj + 1, ck), wp) - u0)*(u0 - real(qc(i)%sf(ci, cj - 1, ck), wp)) <= 0._wp) shy = .false.
            end if
            if (p_glb > 0) then
                if ((real(qc(i)%sf(ci, cj, ck + 1), wp) - u0)*(u0 - real(qc(i)%sf(ci, cj, ck - 1), wp)) <= 0._wp) shz = .false.
            end if
        end do

    end subroutine s_alpha_shared_switch

    !> Disblock prolongation. Guard: no-op unless amr.
    impure subroutine s_populate_amr_fine(q_cons_base)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_base
        integer                                                :: i, islot

        if (.not. amr) return
        ! Prolong EVERY block (max_grid_size tiling can make several) from its gathered coarse patch. The P2P gather pulls each
        ! patch's inter-rank coarse cells from neighbour interiors, so no coarse-ghost halo exchange is needed; host q_cons_base
        ! holds
        ! the ICs here (this runs before s_initialize_gpu_vars). ALL ranks call the gather (P2P); only owners prolong.
        do islot = 1, amr_num_blocks
            call s_amr_select_slot(islot)
            call s_amr_gather_coarse_patch(q_cons_base, .false.)
            ! non-polytropic QBMM: gather the coarse pb/mv patch too (ALL ranks - P2P; owners prolong from it below)
            if (qbmm .and. .not. polytropic) call s_amr_gather_coarse_patch_pbmv(pb_ts(1)%sf, mv_ts(1)%sf, .false.)
            if (amr_rank_owns_block) then
                call s_interpolate_coarse_to_fine()
                do i = 1, sys_size
                    $:GPU_UPDATE(device='[amr_slots(amr_cur)%q_cons(i)%sf]')
                end do
                ! non-polytropic QBMM: seed the block's quadrature side-state from the coarse fields
                if (qbmm .and. .not. polytropic) call s_amr_prolong_pbmv()
            end if
        end do
        if (amr_max_level >= 2) call s_amr_build_static_multilevel(q_cons_base)
        call s_amr_select_slot(1)

    end subroutine s_populate_amr_fine

    !> Build the STATIC multi-level hierarchy (amr_regrid_int = 0): nest exactly one level-2 block inside level-1 block 1 by a fixed
    !! geometric inset (a regrid would place it by sensor-on-fine instead), prolong the parent state into it, and keep it persistent
    !! so the advance driver steps it every timestep. The restrict/reflux identity it relies on is protected by the static
    !! multi-level goldens (75AD6885 et al.) and the runtime conservation-defect probe.
    impure subroutine s_amr_build_static_multilevel(q_cons_base)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_base
        integer                                                :: L2, n1, i, inset(3)

        if (amr_max_level < 2) return  ! np>=2: the L2 is co-located with block 1
        n1 = amr_num_blocks
        if (n1 < 1) return
        ! the static hierarchy nests exactly one level-2 block; without pool room it would SILENTLY refine only to level 1 (an
        ! under-resolved but "successful" run). n1 (the level-1 tile count) is only known here, not at checker time, so abort at the
        ! point of failure. Replicated inputs -> every rank takes the same branch (collective-safe).
        if (n1 + 1 > amr_max_blocks) call s_mpi_abort('amr static multi-level (amr_max_level > 1, amr_regrid_int = 0): ' &
            & // 'amr_max_blocks is too small to nest the level-2 block (need >= level-1 block count + 1); increase amr_max_blocks')
        L2 = n1 + 1
        inset = 0
        inset(1) = max((amr_region_hi_all(1, 1) - amr_region_lo_all(1, 1) + 1)/4, amr_cpat_mar)
        if (n_glb > 0) inset(2) = max((amr_region_hi_all(2, 1) - amr_region_lo_all(2, 1) + 1)/4, amr_cpat_mar)
        if (p_glb > 0) inset(3) = max((amr_region_hi_all(3, 1) - amr_region_lo_all(3, 1) + 1)/4, amr_cpat_mar)
        amr_region_lo_all(:,L2) = amr_region_lo_all(:,1) + inset
        amr_region_hi_all(:,L2) = amr_region_hi_all(:,1) - inset
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
        amr_block_owner(L2) = amr_block_owner(1)
        amr_num_blocks = L2; amr_num_levels = 2
        call s_amr_reconcile_slots()
        amr_cur = L2
        call s_set_amr_fine_geometry(amr_region_lo_all(:,L2), amr_region_hi_all(:,L2))
        call s_amr_gather_coarse_patch(q_cons_base, .false.)  ! q_coarse ignored for level>=2 (reads the parent block); pass the
        ! always-allocated base field, not amr_slots(1) (the parent slot is unallocated on a non-owner rank at np>1)
        if (amr_rank_owns_block) then
            call s_interpolate_coarse_to_fine()
            ! push the host-side prolong to the device (mirror s_populate_amr_fine): s_prolong_one_var is a host loop, so without
            ! this
            ! the persistent L2 block's device q_cons is never valued (NaN) - a GPU-only failure invisible on CPU (host==device)
            do i = 1, sys_size
                $:GPU_UPDATE(device='[amr_slots(amr_cur)%q_cons(i)%sf]')
            end do
        end if
        ! persistent L2 block: KEEP the level-2 block in the active set (amr_num_blocks = L2, amr_num_levels = 2) so the advance
        ! driver steps it across timesteps; no free/revert.
        ! restore amr_cg + the patch frame (amr_cpat_off) to block 1: the L2 gather above overwrote them with the parent-fine frame,
        ! and the normal single-block conservation check that follows reads block 1's frame.
        call s_amr_select_slot(1)
        call s_amr_gather_coarse_patch(q_cons_base, .false.)

    end subroutine s_amr_build_static_multilevel

    !> Volume-weighted restriction: each covered coarse cell = volume-weighted average of its amr_ref_ratio^d fine children (equal
    !! weight on Cartesian grids where children share a volume; radius-weighted by fine y_cc on cyl_coord, where cell volume ~
    !! radius - amr_rvw, single-sourced so device and host paths agree bit-for-bit). Writes the caller's coarse target - in
    !! production the level-0 state q_cons_ts(1)%vf (the deliberate fold-back of fine data each step, plus coarse pb/mv for
    !! non-polytropic QBMM); init-time diagnostics pass a scratch buffer instead. Device kernel.
    impure subroutine s_restrict_fine_to_coarse(coarse_tgt)

        type(scalar_field), dimension(sys_size), intent(inout) :: coarse_tgt
        integer :: i, ci, cj, ck, nchild, rr, dj_hi, dk_hi, o1, o2, o3, owner, r, idx, boxsz, maxsz, nsrc, ierr
        integer :: rlo(3), rhi(3), ilo(3), ihi(3), bl(3), bh(3)
        real(wp), allocatable :: sbuf(:,:), rbuf(:)
        integer, allocatable :: reqs(:), drank(:)

        if (rank_time_wrt .and. amr_rank_owns_block) call s_rank_time_tic()

        ! multi-level: a level>=2 block folds back into its PARENT block's fine array (the coarse side of level l is level l-1), not
        ! the L0 coarse_tgt. Same restriction kernel, targeted at the parent in the parent-fine frame. Co-located towers keep a
        ! level>=2 block and its parent on the same rank, so the restrict is always local (no cross-rank P2P needed).
        if (amr_block_level(amr_cur) >= 2) then
            if (amr_rank_owns_block) call s_amr_restrict_to_parent()
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
                if (bl(1) <= bh(1) .and. bl(2) <= bh(2) .and. bl(3) <= bh(3)) call s_amr_restrict_overwrite_device(coarse_tgt, &
                    & amr_slots(amr_cur)%q_cons, bl, bh, o1, o2, o3, rlo, rr, dj_hi, dk_hi, nchild)
                if (qbmm .and. .not. polytropic .and. amr_rank_owns_block) call s_restrict_pbmv(pb_ts(1)%sf, mv_ts(1)%sf, &
                    & amr_slots(amr_cur)%pb_f%sf, amr_slots(amr_cur)%mv_f%sf)
                if (rank_time_wrt .and. amr_rank_owns_block) call s_rank_time_toc()
                return
            end if
            ! owner-local covered cells: restrict fine(device) -> coarse(device) touching ONLY those cells (no whole-coarse device
            ! push, which clobbered the device-advanced non-covered coarse cells - the same GPU-only bug fixed at np=1)
            if (bl(1) <= bh(1) .and. bl(2) <= bh(2) .and. bl(3) <= bh(3)) call s_amr_restrict_overwrite_device(coarse_tgt, &
                & amr_slots(amr_cur)%q_cons, bl, bh, o1, o2, o3, rlo, rr, dj_hi, dk_hi, nchild)
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
                    call s_amr_restrict_pack_device(amr_slots(amr_cur)%q_cons, bl, bh, rlo, rr, dj_hi, dk_hi, nchild, &
                                                    & sbuf(1:boxsz,nsrc))
#ifdef MFC_MPI
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
                call MPI_RECV(rbuf, boxsz, mpi_p, owner, amr_cur, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
#endif
                idx = 0
                do i = 1, sys_size
                    do ck = bl(3), bh(3)
                        do cj = bl(2), bh(2)
                            do ci = bl(1), bh(1)
                                idx = idx + 1
                                coarse_tgt(i)%sf(ci - o1, cj - o2, ck - o3) = real(rbuf(idx), stp)
                            end do
                        end do
                    end do
                end do
                deallocate (rbuf)
                ! push ONLY the covered cells to the device - a whole-array update would clobber the device-advanced non-covered
                ! coarse cells with this rank's stale host copy (the same GPU-only bug fixed at np=1)
                do i = 1, sys_size
                    $:GPU_UPDATE(device='[coarse_tgt(i)%sf(bl(1) - o1:bh(1) - o1, bl(2) - o2:bh(2) - o2, bl(3) - o3:bh(3) - o3)]')
                end do
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

        integer :: pblk, rr, nchild, dj_hi, dk_hi

        if (.not. amr_rank_owns_block) return  ! np>=2: child+parent co-located; only the owner folds locally
        ! NOTE: reads amr_rvw's device copy without a GPU_UPDATE - safe only while cyl_coord + amr_max_level > 1 is checker-gated
        ! (this path never runs under cyl_coord). If that gate lifts, refresh amr_rvw here first (see its declaration).
        pblk = f_amr_parent_block(amr_cur)
        rr = amr_slots(amr_cur)%amr_ref_ratio
        nchild = rr; if (n_glb > 0) nchild = nchild*rr; if (p_glb > 0) nchild = nchild*rr
        dj_hi = merge(rr - 1, 0, n_glb > 0); dk_hi = merge(rr - 1, 0, p_glb > 0)
        if (amr_isect_lo(1) <= amr_isect_hi(1) .and. amr_isect_lo(2) <= amr_isect_hi(2) .and. amr_isect_lo(3) <= amr_isect_hi(3)) &
            & call s_amr_restrict_overwrite_device(amr_slots(pblk)%q_cons, amr_slots(amr_cur)%q_cons, amr_isect_lo, amr_isect_hi, &
            & 0, 0, 0, amr_isect_lo, rr, dj_hi, dk_hi, nchild)

    end subroutine s_amr_restrict_to_parent

    !> Multi-level reflux: apply the Berger-Colella C/F flux correction from the current level>=2 block into its PARENT block's
    !! cells just OUTSIDE the block footprint, in the parent-fine frame (mirror of the L0 s_amr_apply_reflux targeted at the parent
    !! - "the coarse" is level l-1). State form: q_parent(outside) += dt*(F_coarse - Fbar_fine)/dxf on the low face and +=
    !! dt*(Fbar_fine - F_coarse)/dxf on the high face, where Fbar_fine is the child-averaged fine register. creg/freg key off this
    !! block's slot. np=1 local; the np>=2 P2P freg delivery is future work. Per-face parent-fine dx (stretched-grid safe).
    impure subroutine s_amr_reflux_to_parent(dt_reflux)

        real(wp), intent(in) :: dt_reflux
        integer              :: pblk, d, y, olo(3), ohi(3), glo(3), ghi(3), woff(3)
        real(wp)             :: w_lo(3), w_hi(3), mlo(3), mhi(3)

        if (.not. amr_rank_owns_block) return  ! np>=2: child+parent co-located; only the owner refluxes locally
        pblk = f_amr_parent_block(amr_cur)
        ! max_grid_size tiling of a level>=2 feature: a face shared with an ADJACENT sibling tile (same parent) is fine-fine, not a
        ! c/f boundary - its "outside" parent cell is covered by the sibling's restrict, so refluxing there double-writes and leaks.
        ! Skip those faces (weight 0); the fine-fine halo already matched the shared seam flux. No siblings -> all weights 1
        ! (no-op).
        w_lo = 1._wp; w_hi = 1._wp
        do y = 1, amr_num_blocks  ! block outer, dim inner: f_amr_parent_block (a linear scan) is evaluated once per sibling
            if (y == amr_cur) cycle
            if (f_amr_parent_block(y) /= pblk) cycle  ! same-parent sibling tile only (guarantees same level)
            do d = 1, num_dims
                if (f_amr_seam(amr_cur, y, d)) w_hi(d) = 0._wp  ! sibling just above -> shared high face
                if (f_amr_seam(y, amr_cur, d)) w_lo(d) = 0._wp  ! sibling just below -> shared low face
            end do
        end do
        ! parent-fine frame for the shared reflux kernel: outside cell = isect boundary +/-1; creg-local loop range 0:extent;
        ! transverse write at the isect origin. Per-face parent-fine cell widths - dx at the low/high OUTSIDE cell (olo/ohi),
        ! mirroring the L0/L1 s_amr_apply_reflux_state so a stretched parent grid corrects each C/F face with its own width (on a
        ! uniform grid dx is constant, so this is byte-identical to the previous single-dxf form).
        olo = 0; ohi = 0; glo = 0; ghi = 0; woff = 0; mlo = 1._wp; mhi = 1._wp
        do d = 1, num_dims
            olo(d) = amr_isect_lo(d) - 1; ohi(d) = amr_isect_hi(d) + 1
            ghi(d) = amr_isect_hi(d) - amr_isect_lo(d)
            woff(d) = amr_isect_lo(d)
        end do
        mlo(1) = amr_slots(pblk)%dx(olo(1)); mhi(1) = amr_slots(pblk)%dx(ohi(1))
        if (n_glb > 0) then; mlo(2) = amr_slots(pblk)%dy(olo(2)); mhi(2) = amr_slots(pblk)%dy(ohi(2)); end if
        if (p_glb > 0) then; mlo(3) = amr_slots(pblk)%dz(olo(3)); mhi(3) = amr_slots(pblk)%dz(ohi(3)); end if
        call s_amr_reflux_apply_faces(amr_slots(pblk)%q_cons, amr_cur, amr_slots(amr_cur)%amr_ref_ratio, dt_reflux, olo, ohi, &
                                      & glo, ghi, woff, w_lo, w_hi, mlo, mhi)

    end subroutine s_amr_reflux_to_parent

    !> Rank r's coarse INTERIOR box (global) from the replicated amr_decomp table (no ghosts). Covered coarse cells are in-domain,
    !! so restriction targets are identified by interior overlap alone.
    pure subroutine s_amr_rank_interior(r, ilo, ihi)

        integer, intent(in)  :: r
        integer, intent(out) :: ilo(3), ihi(3)

        ilo = 0; ihi = 0
        ilo(1) = amr_decomp(1, r); ihi(1) = amr_decomp(1, r) + amr_decomp(4, r)
        if (n_glb > 0) then; ilo(2) = amr_decomp(2, r); ihi(2) = amr_decomp(2, r) + amr_decomp(5, r); end if
        if (p_glb > 0) then; ilo(3) = amr_decomp(3, r); ihi(3) = amr_decomp(3, r) + amr_decomp(6, r); end if

    end subroutine s_amr_rank_interior

    !> Device-native restriction overwrite: restrict the fine block q_fine (DEVICE) to coarse averages over the covered coarse cells
    !! [bl:bh] GLOBAL and write coarse_tgt (DEVICE) directly - no host round-trip, only the covered cells touched (the old
    !! whole-coarse device push clobbered non-covered cells). Same child-sum order (ddk, ddj, then fi0 and fi0+1; /nchild; stp cast)
    !! as the old host restrict path, so bit-identical to it on CPU and matches the coarse restriction. q_fine (==
    !! amr_slots(amr_cur)%q_cons) and coarse_tgt are device-resident (ACC_SETUP_SFs). TWIN: s_amr_restrict_pack_device runs this
    !! same child-sum into a wire buffer - any change to the loop order, arithmetic, or casts here must be mirrored there
    !! byte-identically (owner-local and scattered coarse cells must match bit-for-bit). TWIN(q<->pb/mv)
    !! s_amr_restrict_pbmv_box_device runs this same child-sum on pb/mv - keep the stencil lockstep.
    impure subroutine s_amr_restrict_overwrite_device(coarse_tgt, q_fine, bl, bh, o1, o2, o3, rlo, rr, dj_hi, dk_hi, nchild)

        type(scalar_field), dimension(sys_size), intent(inout) :: coarse_tgt
        type(scalar_field), dimension(sys_size), intent(in) :: q_fine
        integer, intent(in) :: bl(3), bh(3), o1, o2, o3, rlo(3), rr, dj_hi, dk_hi, nchild
        integer :: i, ci, cj, ck, fi0, fj0, fk0, ddi, ddj, ddk, bl1, bl2, bl3, bh1, bh2, bh3, rl1, rl2, rl3
        real(wp) :: acc, wacc, w

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
                                        acc = acc + real(q_fine(i)%sf(fi0 + ddi, fj0 + ddj, fk0 + ddk), wp)*w
                                        wacc = wacc + w
                                    end do
                                end do
                            end do
                            coarse_tgt(i)%sf(ci - o1, cj - o2, ck - o3) = real(acc/wacc, stp)
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
                                    acc = acc + real(q_fine(i)%sf(fi0 + ddi, fj0 + ddj, fk0 + ddk), wp)
                                end do
                            end do
                        end do
                        coarse_tgt(i)%sf(ci - o1, cj - o2, ck - o3) = real(acc/real(nchild, wp), stp)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_restrict_overwrite_device

    !> Device pack of one destination's covered restrict slice (np>=2 scatter): restrict the fine block q_fine (DEVICE) over the
    !! covered coarse box [bl:bh] GLOBAL straight into the contiguous wire buffer buf (host, via copyout) - only the slice crosses
    !! PCIe, not the full fine field. Same child-sum order and wp values as s_amr_restrict_overwrite_device (no stp cast: the wire
    !! carries wp and the receiver casts), packed with ci fastest, then cj, ck, i, matching the receiver's sequential unpack. TWIN:
    !! s_amr_restrict_overwrite_device runs this same child-sum in place - any change to the loop order, arithmetic, or casts here
    !! must be mirrored there byte-identically (owner-local and scattered coarse cells must match bit-for-bit). TWIN(q<->pb/mv)
    !! s_amr_restrict_pbmv_pack_device runs this same child-sum into a wire buffer - keep lockstep.
    impure subroutine s_amr_restrict_pack_device(q_fine, bl, bh, rlo, rr, dj_hi, dk_hi, nchild, buf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_fine
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
                                        acc = acc + real(q_fine(i)%sf(fi0 + ddi, fj0 + ddj, fk0 + ddk), wp)*w
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
                                    acc = acc + real(q_fine(i)%sf(fi0 + ddi, fj0 + ddj, fk0 + ddk), wp)
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
        call s_infinite_relaxation_k(amr_slots(amr_cur)%q_cons)
        call s_amr_restore_coarse()

    end subroutine s_amr_relax_fine

    !> 6-equation model: apply the per-stage pressure relaxation to the fine block's interior (cell-local equilibration, no
    !! stencil), mirroring the coarse per-stage call. Swaps the grid so the routine's 0:m,0:n,0:p loop covers this block.
    impure subroutine s_amr_pressure_relax_fine()

        if (.not. amr_rank_owns_block) return
        call s_amr_swap_to_fine()
        call s_pressure_relaxation_procedure(amr_slots(amr_cur)%q_cons)
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
    impure subroutine s_amr_ib_correct_fine()

        if (.not. ib) return
        if (.not. amr_rank_owns_block) return
        call s_amr_swap_to_fine()
        call s_ibm_swap_to_fine(amr_cur, gps_on_device=.true.)
        if (qbmm .and. .not. polytropic) then
            ! mirror the coarse correct-state: non-polytropic QBMM also corrects the block's own pb/mv side-state at the body ghost
            ! points (bounds match the swapped fine idwbuff)
            call s_ibm_correct_state(amr_slots(amr_cur)%q_cons, amr_slots(amr_cur)%q_prim, amr_slots(amr_cur)%pb_f%sf, &
                                     & amr_slots(amr_cur)%mv_f%sf)
        else
            call s_ibm_correct_state(amr_slots(amr_cur)%q_cons, amr_slots(amr_cur)%q_prim)
        end if
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
        integer               :: s, ns, l1, u1, l2, u2, l3, u3
        integer, dimension(6) :: sb1, se1, sb2, se2, sb3, se3
        logical               :: d2, d3

        ox = amr_cpat_off(1); oy = amr_cpat_off(2); oz = amr_cpat_off(3)
        d2 = n_glb > 0; d3 = p_glb > 0
        rr = amr_slots(amr_cur)%amr_ref_ratio
        lo1 = amr_isect_lo(1); lo2 = amr_isect_lo(2); lo3 = amr_isect_lo(3)
        call s_amr_build_ghost_slabs(ns, sb1, se1, sb2, se2, sb3, se3)
        do s = 1, ns
            l1 = sb1(s); u1 = se1(s); l2 = sb2(s); u2 = se2(s); l3 = sb3(s); u3 = se3(s)
            $:GPU_PARALLEL_LOOP(collapse=5, private='[ci, cj, ck]')
            do ib_ = 1, nb
                do q = 1, nnode
                    do fk = l3, u3
                        do fj = l2, u2
                            do fi = l1, u1
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
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do

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
        integer :: ci, cj, ck, q, ib_, nchild, rr, dj_hi, dk_hi, o1, o2, o3, owner, r, idx, boxsz, maxsz, nsrc, ierr, cellsz
        integer :: rlo(3), rhi(3), ilo(3), ihi(3), bl(3), bh(3)
        real(wp), allocatable :: sbuf(:,:), rbuf(:)
        integer, allocatable :: reqs(:), drank(:)

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
                call MPI_RECV(rbuf, boxsz, mpi_p, owner, amr_cur, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
#endif
                idx = 0
                do ib_ = 1, nb
                    do q = 1, nnode
                        do ck = bl(3), bh(3)
                            do cj = bl(2), bh(2)
                                do ci = bl(1), bh(1)
                                    idx = idx + 1
                                    pb_ts(1)%sf(ci - o1, cj - o2, ck - o3, q, ib_) = real(rbuf(idx), stp)
                                end do
                            end do
                        end do
                    end do
                end do
                do ib_ = 1, nb
                    do q = 1, nnode
                        do ck = bl(3), bh(3)
                            do cj = bl(2), bh(2)
                                do ci = bl(1), bh(1)
                                    idx = idx + 1
                                    mv_ts(1)%sf(ci - o1, cj - o2, ck - o3, q, ib_) = real(rbuf(idx), stp)
                                end do
                            end do
                        end do
                    end do
                end do
                deallocate (rbuf)
                ! push ONLY the covered cells to the device (a whole-array update would clobber device-advanced non-covered coarse)
                $:GPU_UPDATE(device='[pb_ts(1)%sf(bl(1) - o1:bh(1) - o1, bl(2) - o2:bh(2) - o2, bl(3) - o3:bh(3) - o3, :, :)]')
                $:GPU_UPDATE(device='[mv_ts(1)%sf(bl(1) - o1:bh(1) - o1, bl(2) - o2:bh(2) - o2, bl(3) - o3:bh(3) - o3, :, :)]')
            end if
        end if

    end subroutine s_amr_scatter_pbmv

    !> Swap the global grid state to the fine block. MUST be paired with s_amr_restore_coarse.
    impure subroutine s_amr_swap_to_fine()

        ! a nested swap would overwrite the sw_* bounce buffers with FINE state, and the eventual restore would install fine extents
        ! as the coarse grid - silent corruption of everything after
        @:ASSERT(.not. amr_swapped, "nested s_amr_swap_to_fine (swap/restore must pair)")
        amr_swapped = .true.
        sw_m = m; sw_n = n; sw_p = p
        sw_idwint = idwint; sw_idwbuff = idwbuff
        ! the acoustic source's precomputed spatials are coarse-grid cell indices: applying them on the fine block would inject at
        ! wrong cells (or out of bounds). The support is guaranteed not to overlap the block (checked at startup), so the fine RHS
        ! correctly skips the source.
        sw_acoustic_source = acoustic_source; acoustic_source = .false.
        ! active-box windows are COARSE cell indices: applying them on the swapped fine grid would window the wrong cells. Blocks
        ! are
        ! contained in the active window (init check + regrid clamp), so the fine advance legitimately treats its whole block as
        ! active.
        sw_ab_active = ab_active; ab_active = .false.
        $:GPU_UPDATE(device='[ab_active]')
        m = amr_slots(amr_cur)%m; n = amr_slots(amr_cur)%n; p = amr_slots(amr_cur)%p
        idwint(1)%beg = 0; idwint(1)%end = m
        idwint(2)%beg = 0; idwint(2)%end = n
        idwint(3)%beg = 0; idwint(3)%end = p
        idwbuff = amr_slots(amr_cur)%idwbuff
        ! save coarse coords to bounce buffers, then copy fine coords into global arrays
        sw_x_cb = x_cb; sw_x_cc = x_cc; sw_dx = dx
        if (n_glb > 0) then; sw_y_cb = y_cb; sw_y_cc = y_cc; sw_dy = dy; end if
        if (p_glb > 0) then; sw_z_cb = z_cb; sw_z_cc = z_cc; sw_dz = dz; end if
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
            integer               :: jg, cl, pblk2, k, rr
            real(wp), allocatable :: cxb(:), cyb(:), czb(:)
            rr = amr_slots(amr_cur)%amr_ref_ratio
            ! ghost parent boundaries: a level>=2 block's coarse side is its PARENT's fine grid (indexed in the parent-fine
            ! amr_isect
            ! frame, matching the interior s_build_level_coords), NOT the L0 global boundaries. amr_isect_lo is a parent-fine index,
            ! so indexing amr_g?cb (sized for L0) reads OUT OF BOUNDS -> garbage on host, NaN on the device copy. Source the
            ! parent's
            ! fine coords for level>=2, the global L0 boundaries for level 1.
            if (amr_block_level(amr_cur) >= 2) then
                pblk2 = f_amr_parent_block(amr_cur)
                allocate (cxb(lbound(amr_slots(pblk2)%x_cb, 1):ubound(amr_slots(pblk2)%x_cb, 1))); cxb = amr_slots(pblk2)%x_cb
                if (n_glb > 0) then
                    allocate (cyb(lbound(amr_slots(pblk2)%y_cb, 1):ubound(amr_slots(pblk2)%y_cb, 1))); cyb = amr_slots(pblk2)%y_cb
                end if
                if (p_glb > 0) then
                    allocate (czb(lbound(amr_slots(pblk2)%z_cb, 1):ubound(amr_slots(pblk2)%z_cb, 1))); czb = amr_slots(pblk2)%z_cb
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

        @:ASSERT(amr_swapped, "s_amr_restore_coarse without a matching s_amr_swap_to_fine")
        amr_swapped = .false.
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
    impure subroutine s_amr_fill_fine_ghosts(q_coarse, q_fine)

        type(scalar_field), dimension(sys_size), intent(in)    :: q_coarse
        type(scalar_field), dimension(sys_size), intent(inout) :: q_fine
        integer                                                :: i, fi, fj, fk, ci, cj, ck, ox, oy, oz
        integer                                                :: rr, lo1, lo2, lo3
        integer                                                :: advb, adve, bbeg, bend, bstride
        integer                                                :: s, ns, l1, u1, l2, u2, l3, u3
        integer, dimension(6)                                  :: sb1, se1, sb2, se2, sb3, se3
        logical                                                :: d2, d3, multi, shx, shy, shz, bubEE
        real(wp)                                               :: u0, sx, sy, sz, xix, xiy, xiz, av, asum

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
        do s = 1, ns
            l1 = sb1(s); u1 = se1(s); l2 = sb2(s); u2 = se2(s); l3 = sb3(s); u3 = se3(s)
            $:GPU_PARALLEL_LOOP(collapse=4, private='[ci, cj, ck, xix, xiy, xiz, u0, sx, sy, sz]')
            do i = 1, sys_size
                do fk = l3, u3
                    do fj = l2, u2
                        do fi = l1, u1
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
                                sx = minmod(real(q_coarse(i)%sf(ci + 1, cj, ck), wp) - u0, u0 - real(q_coarse(i)%sf(ci - 1, cj, &
                                            & ck), wp))
                                sy = 0._wp
                                if (d2) sy = minmod(real(q_coarse(i)%sf(ci, cj + 1, ck), wp) - u0, u0 - real(q_coarse(i)%sf(ci, &
                                    & cj - 1, ck), wp))
                                sz = 0._wp
                                if (d3) sz = minmod(real(q_coarse(i)%sf(ci, cj, ck + 1), wp) - u0, u0 - real(q_coarse(i)%sf(ci, &
                                    & cj, ck - 1), wp))
                                ! QBMM: inject the bub block piecewise-constant (child = u0) so the ghost inherits the coarse cell's
                                ! realizable 6-moment set (CHyQMOM needs variance c20 > 0; per-component minmod slopes would break
                                ! that joint constraint). Non-QBMM Euler-Euler bubbles instead floor their positive moments (nR /
                                ! npb / nmv); the signed velocity moment nV (offset 1) is skipped.
                                if (qbmm .and. i >= bbeg .and. i <= bend) then
                                    sx = 0._wp; sy = 0._wp; sz = 0._wp
                                end if
                                q_fine(i)%sf(fi, fj, fk) = u0 + sx*xix + sy*xiy + sz*xiz
                                if (bubEE .and. .not. qbmm .and. i >= bbeg .and. i <= bend) then
                                    if (mod(i - bbeg, bstride) /= 1) q_fine(i)%sf(fi, fj, fk) = max(real(q_fine(i)%sf(fi, fj, &
                                        & fk), wp), bub_pos_frac*u0)
                                end if
                            end if
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do

        ! multi-fluid volume-fraction ghosts: per-cell closure mirroring s_prolong_alphas_closure (shared limiter switch over all
        ! fluids; interpolate + clamp fluids advb..adve-1; alpha_n = 1 - sum)
        if (multi) then
            do s = 1, ns
                l1 = sb1(s); u1 = se1(s); l2 = sb2(s); u2 = se2(s); l3 = sb3(s); u3 = se3(s)
                $:GPU_PARALLEL_LOOP(collapse=3, private='[i, ci, cj, ck, xix, xiy, xiz, u0, sx, sy, sz, av, asum, shx, shy, shz]')
                do fk = l3, u3
                    do fj = l2, u2
                        do fi = l1, u1
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
                                    if ((real(q_coarse(i)%sf(ci, cj + 1, ck), wp) - u0)*(u0 - real(q_coarse(i)%sf(ci, cj - 1, &
                                        & ck), wp)) <= 0._wp) shy = .false.
                                end if
                                if (d3) then
                                    if ((real(q_coarse(i)%sf(ci, cj, ck + 1), wp) - u0)*(u0 - real(q_coarse(i)%sf(ci, cj, &
                                        & ck - 1), wp)) <= 0._wp) shz = .false.
                                end if
                            end do
                            asum = 0._wp
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = advb, adve - 1
                                u0 = real(q_coarse(i)%sf(ci, cj, ck), wp)
                                sx = 0._wp
                                if (shx) sx = minmod(real(q_coarse(i)%sf(ci + 1, cj, ck), wp) - u0, &
                                    & u0 - real(q_coarse(i)%sf(ci - 1, cj, ck), wp))
                                sy = 0._wp
                                if (shy) sy = minmod(real(q_coarse(i)%sf(ci, cj + 1, ck), wp) - u0, u0 - real(q_coarse(i)%sf(ci, &
                                    & cj - 1, ck), wp))
                                sz = 0._wp
                                if (shz) sz = minmod(real(q_coarse(i)%sf(ci, cj, ck + 1), wp) - u0, u0 - real(q_coarse(i)%sf(ci, &
                                    & cj, ck - 1), wp))
                                av = min(max(u0 + sx*xix + sy*xiy + sz*xiz, 0._wp), 1._wp)
                                q_fine(i)%sf(fi, fj, fk) = av
                                asum = asum + av
                            end do
                            q_fine(adve)%sf(fi, fj, fk) = 1._wp - asum
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end do
        end if

    end subroutine s_amr_fill_fine_ghosts

    !> Lerp the fine ghost shell of q_tgt between q_a (coarse t^n) and q_b (coarse t^{n+1}) at time fraction th (device kernel).
    !! Interior untouched. TWIN s_amr_lerp_fine_ghosts_pbmv (q<->pb/mv): pb/mv sibling of this ghost lerp; keep lockstep.
    impure subroutine s_amr_lerp_fine_ghosts(q_a, q_b, q_tgt, th)

        type(scalar_field), dimension(sys_size), intent(in)    :: q_a, q_b
        type(scalar_field), dimension(sys_size), intent(inout) :: q_tgt
        real(wp), intent(in)                                   :: th
        integer                                                :: i, fi, fj, fk, s, ns, l1, u1, l2, u2, l3, u3
        integer, dimension(6)                                  :: sb1, se1, sb2, se2, sb3, se3

        call s_amr_build_ghost_slabs(ns, sb1, se1, sb2, se2, sb3, se3)
        do s = 1, ns
            l1 = sb1(s); u1 = se1(s); l2 = sb2(s); u2 = se2(s); l3 = sb3(s); u3 = se3(s)
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do fk = l3, u3
                    do fj = l2, u2
                        do fi = l1, u1
                            q_tgt(i)%sf(fi, fj, fk) = (1._wp - th)*real(q_a(i)%sf(fi, fj, fk), wp) + th*real(q_b(i)%sf(fi, fj, &
                                  & fk), wp)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do

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
        integer               :: fi, fj, fk, q, ib_, s, ns, l1, u1, l2, u2, l3, u3
        integer, dimension(6) :: sb1, se1, sb2, se2, sb3, se3

        call s_amr_build_ghost_slabs(ns, sb1, se1, sb2, se2, sb3, se3)
        do s = 1, ns
            l1 = sb1(s); u1 = se1(s); l2 = sb2(s); u2 = se2(s); l3 = sb3(s); u3 = se3(s)
            $:GPU_PARALLEL_LOOP(collapse=5)
            do ib_ = 1, nb
                do q = 1, nnode
                    do fk = l3, u3
                        do fj = l2, u2
                            do fi = l1, u1
                                pb_t(fi, fj, fk, q, ib_) = (1._wp - th)*real(pga(fi, fj, fk, q, ib_), wp) + th*real(pgb(fi, fj, &
                                     & fk, q, ib_), wp)
                                mv_t(fi, fj, fk, q, ib_) = (1._wp - th)*real(mga(fi, fj, fk, q, ib_), wp) + th*real(mgb(fi, fj, &
                                     & fk, q, ib_), wp)
                            end do
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do

    end subroutine s_amr_lerp_fine_ghosts_pbmv

    !> Device copy q_src -> q_dst over [b1:e1, b2:e2, b3:e3] for all sys_size fields (RK step-entry backup). TWIN s_amr_backup_pbmv
    !! (q<->pb/mv): pb/mv sibling of this step-entry backup; keep lockstep.
    impure subroutine s_amr_copy_fine_fields(q_src, q_dst, b1, e1, b2, e2, b3, e3)

        type(scalar_field), dimension(sys_size), intent(in)    :: q_src
        type(scalar_field), dimension(sys_size), intent(inout) :: q_dst
        integer, intent(in)                                    :: b1, e1, b2, e2, b3, e3
        integer                                                :: i, fi, fj, fk

        $:GPU_PARALLEL_LOOP(collapse=4)
        do i = 1, sys_size
            do fk = b3, e3
                do fj = b2, e2
                    do fi = b1, e1
                        q_dst(i)%sf(fi, fj, fk) = q_src(i)%sf(fi, fj, fk)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_copy_fine_fields

    !> Device RK stage update over the fine interior: q = (c1*q + c2*q_stor + c3*dt_in*rhs)/c4 (compute in wp, store stp). Mirrors
    !! the coarse non-IGR rk_coef form in s_tvd_rk. TWIN s_amr_fine_rk_update_pbmv + s_tvd_rk (m_time_steppers): same SSP-RK stage
    !! combination; keep all three lockstep.
    impure subroutine s_amr_fine_rk_update(q_upd, q_stor, q_rhs, c1, c2, c3, c4, dt_in)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_upd
        type(scalar_field), dimension(sys_size), intent(in)    :: q_stor, q_rhs
        real(wp), intent(in)                                   :: c1, c2, c3, c4, dt_in
        integer                                                :: i, fi, fj, fk, fm, fn, fp

        fm = amr_slots(amr_cur)%m; fn = amr_slots(amr_cur)%n; fp = amr_slots(amr_cur)%p
        $:GPU_PARALLEL_LOOP(collapse=4)
        do i = 1, sys_size
            do fk = 0, fp
                do fj = 0, fn
                    do fi = 0, fm
                        q_upd(i)%sf(fi, fj, fk) = (c1*real(q_upd(i)%sf(fi, fj, fk), wp) + c2*real(q_stor(i)%sf(fi, fj, fk), &
                              & wp) + c3*dt_in*real(q_rhs(i)%sf(fi, fj, fk), wp))/c4
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

    !> True iff sub-block yb sits immediately above sub-block xb on xb's HIGH face in dim d: yb's low coarse face is xb's high face
    !! + 1, and (tiling produces a regular grid) they share the transverse coarse extents exactly. Each fine-fine seam is exactly
    !! one such ordered (xb, yb) pair (the lower block is xb).
    pure logical function f_amr_seam(xb, yb, d) result(seam)

        integer, intent(in) :: xb, yb, d
        integer             :: t

        seam = amr_region_lo_all(d, yb) == amr_region_hi_all(d, xb) + 1
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
    impure subroutine s_amr_fine_slice(slot, q_cons, d, dlo, dhi, buf, dir)

        integer, intent(in)                                    :: slot, d, dlo, dhi, dir
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons  ! amr_slots(slot)%q_cons: pass the
        !                                    scalar_field array in so the device kernel reads the mapped %sf via a dummy - amr_slots
        !                                    itself is not GPU_DECLARE'd, so indexing the module array amr_slots(slot)%q_cons inside
        !                                    a kernel dereferences a null descriptor
        real(wp), intent(inout), contiguous :: buf(:)
        integer                             :: i, a, b, c, fm(3), na, nb, nc

        fm(1) = amr_slots(slot)%m; fm(2) = amr_slots(slot)%n; fm(3) = amr_slots(slot)%p
        nc = dhi - dlo + 1
        ! Pack (dir=1) / unpack (dir=-1) the near-seam slab ON THE DEVICE straight into the contiguous buffer buf, then move ONLY
        ! buf
        ! host<->device. flang miscomputes a STRIDED section (seam dim d < num_dims) of the doubly-nested amr_slots%q_cons%sf in a
        ! target-update map clause, corrupting the 2D+ np>1 seam ghosts; the base-grid halo (s_mpi_sendrecv_variables_buffers)
        ! device-packs into a contiguous buffer for the same reason. buf index runs a fastest, then b, then c, then i, so a pack and
        ! an unpack with matching extents align cell-for-cell (na/nb are the transverse fine sizes, nc the slab depth).
        #:for D, TA, TB in [(1, 2, 3), (2, 1, 3), (3, 1, 2)]
            #:set IDX = {1: '(c, a, b)', 2: '(a, c, b)', 3: '(a, b, c)'}[D]
            if (d == ${D}$) then
                na = fm(${TA}$) + 1; nb = fm(${TB}$) + 1  ! scalars; kernel loop bounds MUST use na-1/nb-1, not fm(..), so no host
                !                     array is referenced in the device region (nvfortran/Cray demand it PRESENT)
                if (dir == 1) then  ! host <- device: pack on the device, copyout moves the contiguous buffer to host
                    $:GPU_PARALLEL_LOOP(collapse=4, copyout='[buf]')
                    do i = 1, sys_size
                        do c = dlo, dhi
                            do b = 0, nb - 1
                                do a = 0, na - 1
                                    buf(1 + a + na*(b + nb*(c - dlo + nc*(i - 1)))) = real(q_cons(i)%sf${IDX}$, wp)
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
                                    q_cons(i)%sf${IDX}$ = real(buf(1 + a + na*(b + nb*(c - dlo + nc*(i - 1)))), stp)
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
            end if
        #:endfor

    end subroutine s_amr_fine_slice

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
    !! rebuilds the per-block gather/scatter overlap-rank lists (amr_ovl_gather/scatter): one O(nblocks*num_procs) scan of the
    !! static amr_decomp table here in place of the same scan per block per RK stage in the P2P gather/scatter coupling.
    impure subroutine s_amr_build_seam_pairs()

        integer :: xb, yb, d, np, k, r
        integer :: plo(3), phi(3), rlo(3), rhi(3), clo(3), chi(3), bl(3), bh(3)

        if (allocated(amr_seam_pairs)) deallocate (amr_seam_pairs)
        np = 0
        do xb = 1, amr_num_blocks
            do yb = 1, amr_num_blocks
                if (f_amr_seam_dim(xb, yb) > 0) np = np + 1
            end do
        end do
        amr_num_seam_pairs = np
        allocate (amr_seam_pairs(3, max(np, 1)))
        np = 0
        do xb = 1, amr_num_blocks
            do yb = 1, amr_num_blocks
                d = f_amr_seam_dim(xb, yb)
                if (d > 0) then
                    np = np + 1
                    amr_seam_pairs(1, np) = xb; amr_seam_pairs(2, np) = yb; amr_seam_pairs(3, np) = d
                end if
            end do
        end do
        ! per-block P2P overlap-rank lists (gather: rank coarse range vs the amr_cpat_mar-padded patch box; scatter: rank coarse
        ! interior vs the region box - the exact tests the consumers ran per call), rank-ascending so iterating a list preserves the
        ! replaced 0..num_procs-1 scans' MPI send/recv order.
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
            amr_ovl_gather_n(k) = 0; amr_ovl_scatter_n(k) = 0
            do r = 0, num_procs - 1
                call s_amr_rank_coarse_range(r, clo, chi)
                call s_amr_box_isect(plo, phi, clo, chi, bl, bh)
                if (bl(1) <= bh(1) .and. bl(2) <= bh(2) .and. bl(3) <= bh(3)) then
                    amr_ovl_gather_n(k) = amr_ovl_gather_n(k) + 1
                    amr_ovl_gather(amr_ovl_gather_n(k), k) = r
                end if
                call s_amr_rank_interior(r, clo, chi)
                call s_amr_box_isect(rlo, rhi, clo, chi, bl, bh)
                if (bl(1) <= bh(1) .and. bl(2) <= bh(2) .and. bl(3) <= bh(3)) then
                    amr_ovl_scatter_n(k) = amr_ovl_scatter_n(k) + 1
                    amr_ovl_scatter(amr_ovl_scatter_n(k), k) = r
                end if
            end do
        end do
        amr_seam_pairs_nblk = amr_num_blocks
        amr_seam_pairs_dirty = .false.

    end subroutine s_amr_build_seam_pairs

    !> Block-to-block fine-fine halo (max_grid_size tiling): overwrite each sub-block's seam ghost cells (faces shared with an
    !! ADJACENT sub-block) with the neighbour's stage-entry fine interior, so the shared fine flux matches on both sides
    !! (coarse-prolonged seam ghosts would be non-conservative). For each seam pair (xb below, yb above, dim d) the two owners
    !! exchange the buff_size-deep near-seam interior (MPI_Sendrecv, or a local copy when one rank owns both). Buffer is wp, cast to
    !! stp on unpack (identity for stp fields). No-op with a single block / no adjacent pairs (any untiled case, any np).
    impure subroutine s_amr_fine_fine_halo()

        integer :: xb, yb, d, rX, rY, cnt, xm(3), ym(3), tsz, ierr, fmul, idx

        if (.not. amr) return
        if (amr_num_blocks < 2) return

        ! iterate the cached same-level seam list (rebuilt only when the topology changes) instead of the old O(nblocks^2)
        ! f_amr_seam
        ! rescan every stage; the list preserves the original (xb, yb) order so paired MPI_SENDRECVs still match.
        if (amr_seam_pairs_dirty .or. amr_seam_pairs_nblk /= amr_num_blocks) call s_amr_build_seam_pairs()
        ! device<->host of the fine state is done per-seam inside s_amr_fine_slice, moving only the buff_size-deep near-seam slab
        ! each
        ! pack/unpack touches (not the whole block) - a large PCIe saving since this runs per stage (6x per fine step)
        do idx = 1, amr_num_seam_pairs
            xb = amr_seam_pairs(1, idx); yb = amr_seam_pairs(2, idx); d = amr_seam_pairs(3, idx)
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
            ym(1) = fmul*(amr_region_hi_all(1, yb) - amr_region_lo_all(1, yb) + 1) - 1
            ym(2) = merge(fmul*(amr_region_hi_all(2, yb) - amr_region_lo_all(2, yb) + 1) - 1, 0, n_glb > 0)
            ym(3) = merge(fmul*(amr_region_hi_all(3, yb) - amr_region_lo_all(3, yb) + 1) - 1, 0, p_glb > 0)
            ! transverse fine size (dims /= d); xb and yb share it (exact-match seam)
            tsz = 1
            if (d /= 1) tsz = tsz*(xm(1) + 1)
            if (d /= 2 .and. n_glb > 0) tsz = tsz*(xm(2) + 1)
            if (d /= 3 .and. p_glb > 0) tsz = tsz*(xm(3) + 1)
            cnt = sys_size*buff_size*tsz
            if (rX == rY) then  ! same rank owns both: pack each near-seam interior, unpack into the other's seam ghost
                ! xb high interior
                call s_amr_fine_slice(xb, amr_slots(xb)%q_cons, d, xm(d) - buff_size + 1, xm(d), amr_seambuf_x(1:cnt), 1)
                call s_amr_fine_slice(yb, amr_slots(yb)%q_cons, d, 0, buff_size - 1, amr_seambuf_y(1:cnt), 1)  ! yb low interior
                call s_amr_fine_slice(yb, amr_slots(yb)%q_cons, d, -buff_size, -1, amr_seambuf_x(1:cnt), -1)  ! -> yb low ghost
                ! -> xb high ghost
                call s_amr_fine_slice(xb, amr_slots(xb)%q_cons, d, xm(d) + 1, xm(d) + buff_size, amr_seambuf_y(1:cnt), -1)
            else if (proc_rank == rX) then
                ! send xb high interior
                call s_amr_fine_slice(xb, amr_slots(xb)%q_cons, d, xm(d) - buff_size + 1, xm(d), amr_seambuf_x(1:cnt), 1)
#ifdef MFC_MPI
                call MPI_SENDRECV(amr_seambuf_x(1:cnt), cnt, mpi_p, rY, 4200, amr_seambuf_y(1:cnt), cnt, mpi_p, rY, 4201, &
                                  & MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
#endif
                ! recv yb low interior -> xb high ghost
                call s_amr_fine_slice(xb, amr_slots(xb)%q_cons, d, xm(d) + 1, xm(d) + buff_size, amr_seambuf_y(1:cnt), -1)
            else  ! proc_rank == rY
                ! send yb low interior
                call s_amr_fine_slice(yb, amr_slots(yb)%q_cons, d, 0, buff_size - 1, amr_seambuf_y(1:cnt), 1)
#ifdef MFC_MPI
                call MPI_SENDRECV(amr_seambuf_y(1:cnt), cnt, mpi_p, rX, 4201, amr_seambuf_x(1:cnt), cnt, mpi_p, rX, 4200, &
                                  & MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
#endif
                ! recv xb high interior -> yb low ghost
                call s_amr_fine_slice(yb, amr_slots(yb)%q_cons, d, -buff_size, -1, amr_seambuf_x(1:cnt), -1)
            end if
        end do
        call s_amr_select_slot(1)

    end subroutine s_amr_fine_fine_halo

    !> FILL phase of a fine RK stage (same dt as level-0, no subcycling; called BETWEEN the coarse RHS and the coarse RK update, so
    !! q_cons_coarse is the coarse STAGE-ENTRY state): valid coarse ghosts + gather the block's coarse patch + prolong the block's
    !! fine ghost shell from it. ALL ranks call (gather is collective P2P); only the owner prolongs. Leaves the block's INTERIOR at
    !! its stage-entry value (the advance phase updates it), so a later fine-fine halo reads stage-entry neighbour data.
    impure subroutine s_amr_fine_stage_fill(q_cons_coarse, pb_in, mv_in)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_coarse
        real(stp), dimension(:,:,:,:,:), intent(inout)         :: pb_in, mv_in

        if (.not. amr) return
        ! valid coarse CONS ghosts for the ghost prolongation below (ALL ranks call: pairwise halo)
        if (amr_xchg_coarse_ghosts) call s_amr_exchange_coarse_cons_halo(q_cons_coarse)
        ! the Lagrangian-cloud guard runs on EVERY rank (rank-local bubbles vs the block's GLOBAL box): a non-owner rank's bubbles
        ! can reach into a block across a rank seam
        call s_amr_check_lag_clear()
        ! fine-level distribution: gather the block's coarse patch onto its owner (collective - ALL ranks, before the owner-only
        ! return). The device ghost-fill below then reads the gathered patch amr_cg instead of local coarse.
        call s_amr_gather_coarse_patch(q_cons_coarse, .true.)
        ! non-polytropic QBMM: gather the coarse pb/mv patch too (collective - ALL ranks, before the owner return; owners fill
        ! below)
        if (qbmm .and. .not. polytropic) call s_amr_gather_coarse_patch_pbmv(pb_in, mv_in, .true.)
        if (.not. amr_rank_owns_block) return

        ! ghost prolongation from the coarse stage-entry conservative state (device kernel reads the gathered patch amr_cg);
        ! rank_time brackets cover the fine-advance compute segments and pause across the MPI exchanges (the inner s_compute_rhs
        ! pair
        ! nests to a no-op)
        if (rank_time_wrt) call s_rank_time_tic()
        call s_amr_fill_fine_ghosts(amr_cg, amr_slots(amr_cur)%q_cons)
        ! the coarse-prolonged ghost shell is the final ghost state EXCEPT at faces shared with an adjacent sub-block (tiling),
        ! which
        ! the block-to-block fine-fine halo overwrites with the neighbour's fine interior after this fill.

        ! non-polytropic QBMM: prolong the block's pb/mv ghost shell from the gathered coarse patch amr_cg_pb/mv (its rhs is
        ! cell-local, so ghosts only feed the widened-idwint conversions)
        if (qbmm .and. .not. polytropic) call s_amr_fill_fine_ghosts_pbmv(amr_cg_pb, amr_cg_mv, amr_slots(amr_cur)%pb_f%sf, &
            & amr_slots(amr_cur)%mv_f%sf)
        if (rank_time_wrt) call s_rank_time_toc()

    end subroutine s_amr_fine_stage_fill

    !> ADVANCE phase of a fine RK stage: fine RHS + RK update (+ QBMM/6eq/IB) for the current block. Owner-only. Reads the block's
    !! ghost shell (coarse prolong + fine-fine halo already applied by the fill + halo phases).
    impure subroutine s_amr_fine_stage_advance(s, coefs, bc_type, q_T_sf, pb_in, rhs_pb, mv_in, rhs_mv, t_step)

        integer, intent(in)                                        :: s, t_step
        real(wp), intent(in)                                       :: coefs(4)
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        type(scalar_field), intent(inout)                          :: q_T_sf
        real(stp), dimension(:,:,:,:,:), intent(inout)             :: pb_in, mv_in
        real(wp), dimension(:,:,:,:,:), intent(inout)              :: rhs_pb, rhs_mv

        if (.not. amr) return
        if (.not. amr_rank_owns_block) return
        if (rank_time_wrt) call s_rank_time_tic()

        ! step-entry backup for the SSP-RK combination (device copy over the current buffered extents)
        if (s == 1) then
            call s_amr_copy_fine_fields(amr_slots(amr_cur)%q_cons, amr_slots(amr_cur)%q_cons_stor, &
                                        & amr_slots(amr_cur)%idwbuff(1)%beg, amr_slots(amr_cur)%idwbuff(1)%end, &
                                        & amr_slots(amr_cur)%idwbuff(2)%beg, amr_slots(amr_cur)%idwbuff(2)%end, &
                                        & amr_slots(amr_cur)%idwbuff(3)%beg, amr_slots(amr_cur)%idwbuff(3)%end)
            if (qbmm .and. .not. polytropic) call s_amr_backup_pbmv(amr_slots(amr_cur)%pb_f%sf, amr_slots(amr_cur)%mv_f%sf, &
                & amr_slots(amr_cur)%pb_stor%sf, amr_slots(amr_cur)%mv_stor%sf)
        end if

        amr_in_fine_advance = .true.
        call s_amr_swap_to_fine()
        idwint = amr_slots(amr_cur)%idwbuff  ! widen the conversion range to the ghost shell (restored by s_amr_restore_coarse)
        $:GPU_UPDATE(device='[idwint]')
        if (qbmm .and. .not. polytropic) then
            ! the block's OWN side-state and rhs scratch: the coarse pb_in/rhs_pb must not be touched at fine indices (the coarse
            ! stage consumes them after this fine stage)
            call s_compute_rhs(amr_slots(amr_cur)%q_cons, q_T_sf, amr_slots(amr_cur)%q_prim, bc_type, amr_slots(amr_cur)%rhs, &
                               & amr_slots(amr_cur)%pb_f%sf, amr_rhs_pb_f, amr_slots(amr_cur)%mv_f%sf, amr_rhs_mv_f, t_step, s)
        else
            call s_compute_rhs(amr_slots(amr_cur)%q_cons, q_T_sf, amr_slots(amr_cur)%q_prim, bc_type, amr_slots(amr_cur)%rhs, &
                               & pb_in, rhs_pb, mv_in, rhs_mv, t_step, s)
        end if
        call s_amr_restore_coarse()
        amr_in_fine_advance = .false.

        ! RK stage update (device kernel; mirror of the coarse form - under IGR the rhs already embeds dt, matching the coarse igr
        ! update, so the dt factor is 1)
        call s_amr_fine_rk_update(amr_slots(amr_cur)%q_cons, amr_slots(amr_cur)%q_cons_stor, amr_slots(amr_cur)%rhs, coefs(1), &
                                  & coefs(2), coefs(3), coefs(4), merge(1._wp, dt, igr))
        if (qbmm .and. .not. polytropic) call s_amr_fine_rk_update_pbmv(amr_slots(amr_cur)%pb_f%sf, amr_slots(amr_cur)%mv_f%sf, &
            & amr_slots(amr_cur)%pb_stor%sf, amr_slots(amr_cur)%mv_stor%sf, amr_rhs_pb_f, amr_rhs_mv_f, coefs(1), coefs(2), &
            & coefs(3), coefs(4), dt)
        ! 6-equation model: per-stage pressure relaxation on the block (before IB correct, coarse order)
        if (model_eqns == model_eqns_6eq .and. (.not. relax)) call s_amr_pressure_relax_fine()
        ! moving body: rebuild the fine-block IB state at the current (lockstep-stage) body position before the correct-state
        if (moving_immersed_boundary_flag) call s_amr_update_mib_fine(-1._wp)
        ! IB state correction on the fine block (mirrors the coarse per-stage correct-state; no-op unless ib)
        call s_amr_ib_correct_fine()
        if (rank_time_wrt) call s_rank_time_toc()

    end subroutine s_amr_fine_stage_advance

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

        if (amr_xchg_coarse_ghosts) then
            call s_amr_exchange_coarse_cons_halo(q_old)
            call s_amr_exchange_coarse_cons_halo(q_new)
        end if
        call s_amr_check_lag_clear()  ! EVERY rank: non-owner bubbles can reach the block across a seam

        ! fine-level distribution: gather each lerp source's coarse patch (collective - ALL ranks) then prolong its ghost shell on
        ! the owner. Interleaved so the single amr_cg buffer is consumed by the fill before the next gather overwrites it.
        ! non-polytropic QBMM: the pb/mv ghost shell gets the same two-source time-lerp treatment, gathered + filled INTERLEAVED
        ! with
        ! q_cons so the single amr_cg_pb/mv buffer is consumed by each fill before the next gather overwrites it. Gathers collective
        ! (ALL ranks - P2P); fills owner-only.
        call s_amr_gather_coarse_patch(q_old, .true.)
        if (qbmm .and. .not. polytropic) call s_amr_gather_coarse_patch_pbmv(pb_old, mv_old, .true.)
        if (amr_rank_owns_block) call s_amr_fill_fine_ghosts(amr_cg, amr_slots(amr_cur)%q_ghost_a)
        if (amr_rank_owns_block .and. qbmm .and. .not. polytropic) call s_amr_fill_fine_ghosts_pbmv(amr_cg_pb, amr_cg_mv, &
            & amr_slots(amr_cur)%pb_ghost_a%sf, amr_slots(amr_cur)%mv_ghost_a%sf)
        call s_amr_gather_coarse_patch(q_new, .true.)
        if (qbmm .and. .not. polytropic) call s_amr_gather_coarse_patch_pbmv(pb_in, mv_in, .true.)
        if (amr_rank_owns_block) call s_amr_fill_fine_ghosts(amr_cg, amr_slots(amr_cur)%q_ghost_b)
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
    !! (s_amr_advance_children); the L2-L2 seam halo is future work (s_amr_check_seam_topology aborts if an L2+ seam pair is reached
    !! under subcycle, e.g. via a restart mode-switch from a lockstep-produced layout). A single owned level-1 block is
    !! byte-identical to the old per-block subcycle (the halo is a no-op with < 2 adjacent same-level blocks).
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
                call s_amr_fine_fine_halo()
                ! RHS + RK update every block from the reconciled ghost shell
                do islot = 1, amr_num_blocks
                    if (amr_block_level(islot) /= 1) cycle
                    call s_amr_select_slot(islot)
                    if (.not. amr_rank_owns_block) cycle
                    call s_amr_subtree_stage_advance(amr_dt_fine, coefs, bc_type, q_T_sf, pb_in, rhs_pb, mv_in, rhs_mv, t_step, &
                                                     & s, th)
                end do
            end do
            ! after this substep each level-1 block is at t_b (q_cons) with t_a in q_cons_stor: its level-2 children subcycle within
            ! [t_a, t_b] then fold back (restrict + Berger-Colella reflux). No-op for single-level.
            if (amr_max_level >= 2) then
                do islot = 1, amr_num_blocks
                    if (amr_block_level(islot) /= 1) cycle
                    call s_amr_select_slot(islot)
                    if (.not. amr_rank_owns_block) cycle
                    call s_amr_advance_children(islot, amr_dt_fine, coefs, bc_type, q_T_sf, pb_in, rhs_pb, mv_in, rhs_mv, t_step)
                end do
            end if
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
        call s_amr_lerp_fine_ghosts(amr_slots(amr_cur)%q_ghost_a, amr_slots(amr_cur)%q_ghost_b, amr_slots(amr_cur)%q_cons, th)
        if (qbmm .and. .not. polytropic) call s_amr_lerp_fine_ghosts_pbmv(amr_slots(amr_cur)%pb_f%sf, amr_slots(amr_cur)%mv_f%sf, &
            & amr_slots(amr_cur)%pb_ghost_a%sf, amr_slots(amr_cur)%mv_ghost_a%sf, amr_slots(amr_cur)%pb_ghost_b%sf, &
            & amr_slots(amr_cur)%mv_ghost_b%sf, th)

        ! substep-entry backup for the SSP-RK combination (device copy, interior only)
        if (s == 1) then
            call s_amr_copy_fine_fields(amr_slots(amr_cur)%q_cons, amr_slots(amr_cur)%q_cons_stor, 0, amr_slots(amr_cur)%m, 0, &
                                        & amr_slots(amr_cur)%n, 0, amr_slots(amr_cur)%p)
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
        if (qbmm .and. .not. polytropic) then
            ! the block's OWN side-state and rhs scratch (the coarse arrays stay untouched)
            call s_compute_rhs(amr_slots(amr_cur)%q_cons, q_T_sf, amr_slots(amr_cur)%q_prim, bc_type, amr_slots(amr_cur)%rhs, &
                               & amr_slots(amr_cur)%pb_f%sf, amr_rhs_pb_f, amr_slots(amr_cur)%mv_f%sf, amr_rhs_mv_f, t_step, s)
        else
            call s_compute_rhs(amr_slots(amr_cur)%q_cons, q_T_sf, amr_slots(amr_cur)%q_prim, bc_type, amr_slots(amr_cur)%rhs, &
                               & pb_in, rhs_pb, mv_in, rhs_mv, t_step, s)
        end if
        call s_amr_restore_coarse()
        amr_in_fine_advance = .false.

        ! RK stage update at the FINE time step (device kernel)
        call s_amr_fine_rk_update(amr_slots(amr_cur)%q_cons, amr_slots(amr_cur)%q_cons_stor, amr_slots(amr_cur)%rhs, coefs(s, 1), &
                                  & coefs(s, 2), coefs(s, 3), coefs(s, 4), dt_sub)
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
        call s_amr_ib_correct_fine()
        if (rank_time_wrt) call s_rank_time_toc()

    end subroutine s_amr_subtree_stage_advance

    !> Advance the currently-selected fine block (amr_cur) through its subcycle: two amr_dt_fine SSP-RK3 substeps whose ghost shell
    !! is the two-source time-lerp of the block's q_ghost_a (parent t^n) / q_ghost_b (parent t^{n+1}) sources, filled by the caller
    !! before this call. Fine flux registers are zeroed here and accumulate over all six stages so the end-of-step reflux sees the
    !! time-averaged effective fine flux. RECURSIVE: after each of this block's substeps, its level+1 children subcycle within that
    !! substep (s_amr_advance_children) at dt_sub/amr_ref_ratio, so per-level dt falls out of the recursion. With no children this
    !! is exactly the single-level L0<->L1 advance. Used for level>=2 children (per-block); level-1 blocks run the transposed,
    !! seam-halo'd s_amr_advance_fine_subcycle_all instead.
    recursive subroutine s_amr_advance_subtree(dt_sub, coefs, bc_type, q_T_sf, pb_in, rhs_pb, mv_in, rhs_mv, t_step)

        real(wp), intent(in) :: dt_sub                 !< this block's substep dt (parent step / amr_ref_ratio)
        real(wp), dimension(:,:), intent(in) :: coefs  !< rk_coef(1:3, 1:4)
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        type(scalar_field), intent(inout) :: q_T_sf
        real(stp), dimension(:,:,:,:,:), intent(inout) :: pb_in, mv_in
        real(wp), dimension(:,:,:,:,:), intent(inout) :: rhs_pb, rhs_mv
        integer, intent(in) :: t_step
        real(wp), parameter :: c_abs(3) = [0._wp, 1._wp, 0.5_wp]
        integer :: sub, s
        real(wp) :: th

        call s_amr_zero_fine_registers()
        do sub = 1, 2
            do s = 1, 3
                th = (real(sub - 1, wp) + c_abs(s))*0.5_wp
                call s_amr_subtree_stage_lerp(s, th)
                call s_amr_subtree_stage_advance(dt_sub, coefs, bc_type, q_T_sf, pb_in, rhs_pb, mv_in, rhs_mv, t_step, s, th)
            end do
            ! multi-level: this block is now at t_b (q_cons) with t_a preserved in q_cons_stor. Each level+1 child subcycles WITHIN
            ! this substep - gathering its two ghost-lerp sources from those two snapshots - then folds back (restrict +
            ! Berger-Colella
            ! reflux) into this block over dt_sub. No-op when this block has no children (single-level / finest).
            if (amr_max_level >= 2) call s_amr_advance_children(amr_cur, dt_sub, coefs, bc_type, q_T_sf, pb_in, rhs_pb, mv_in, &
                & rhs_mv, t_step)
        end do

    end subroutine s_amr_advance_subtree

    !> Recursively subcycle every level+1 child of parent slot pslot within ONE of the parent's substeps [t_a, t_b] (duration
    !! dt_sub). The parent has just finished that substep: q_cons = parent @ t_b, q_cons_stor = parent @ t_a. For each child: gather
    !! its two ghost-lerp sources from those two parent snapshots (parent-fine frame), recurse s_amr_advance_subtree at dt_sub/2
    !! (the child takes amr_ref_ratio substeps covering [t_a, t_b]), then fold the child back into the parent - restrict the covered
    !! cells and apply the Berger-Colella C/F flux correction (s_amr_reflux_to_parent over dt_sub, consuming the child's freg + the
    !! parent-side creg captured during THIS substep). The registers already carry the matching per-substep time weights (freg
    !! 1/r*rk3_w, creg rk3_w), so conservation closes with no register changes. np=1 (child co-owned with parent); np>=2 P2P is
    !! future work (#27).
    recursive subroutine s_amr_advance_children(pslot, dt_sub, coefs, bc_type, q_T_sf, pb_in, rhs_pb, mv_in, rhs_mv, t_step)

        integer, intent(in)                                        :: pslot
        real(wp), intent(in)                                       :: dt_sub
        real(wp), dimension(:,:), intent(in)                       :: coefs
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        type(scalar_field), intent(inout)                          :: q_T_sf
        real(stp), dimension(:,:,:,:,:), intent(inout)             :: pb_in, mv_in
        real(wp), dimension(:,:,:,:,:), intent(inout)              :: rhs_pb, rhs_mv
        integer, intent(in)                                        :: t_step
        integer                                                    :: kc, pslot_l

        ! pslot is argument-associated with the module variable amr_cur (the caller passes amr_cur as pslot); the
        ! s_amr_select_slot(kc) below reassigns amr_cur and would silently corrupt pslot. Copy it to a local read before any slot
        ! switch.

        pslot_l = pslot
        do kc = 1, amr_num_blocks
            if (amr_block_level(kc) /= amr_block_level(pslot_l) + 1) cycle
            if (f_amr_parent_block(kc) /= pslot_l) cycle
            call s_amr_select_slot(kc)  ! amr_cur = kc; mirrors (isect already parent-fine)
            if (.not. amr_rank_owns_block) cycle  ! np>=2: child on another rank - future work (#27)
            ! two ghost-lerp sources from the parent's substep endpoints (parent-fine frame)
            call s_amr_gather_from_parent_field(pslot_l, amr_slots(pslot_l)%q_cons_stor, .false.)  ! parent @ t_a (device C/F fill)
            call s_amr_fill_fine_ghosts(amr_cg, amr_slots(kc)%q_ghost_a)
            call s_amr_gather_from_parent_field(pslot_l, amr_slots(pslot_l)%q_cons, .false.)  ! parent @ t_b (device C/F fill)
            call s_amr_fill_fine_ghosts(amr_cg, amr_slots(kc)%q_ghost_b)
            ! recurse: the child subcycles its amr_ref_ratio substeps (and its own children) over [t_a, t_b]
            call s_amr_advance_subtree(dt_sub*0.5_wp, coefs, bc_type, q_T_sf, pb_in, rhs_pb, mv_in, rhs_mv, t_step)
            ! fold the child back into the parent (relax the fine phase first, matching the driver's relax -> restrict order)
            if (relax) call s_amr_relax_fine()
            call s_amr_restrict_to_parent()
            call s_amr_reflux_to_parent(dt_sub)
        end do
        call s_amr_select_slot(pslot_l)  ! restore the parent's mirrors for its next substep

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

        integer, intent(in)           :: lo(3), hi(3), cap
        type(t_box), intent(inout)    :: out(:)
        integer, intent(inout)        :: nt, capped
        integer, intent(in), optional :: tsz(3)  !< per-dim tile size (default amr_maxc_fit; level>=2 passes amr_maxc_fit/2)
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
    !! amr_cg/amr_decomp are NOT per-slot and stay in init/finalize.
    impure subroutine s_amr_alloc_slot(islot)

        integer, intent(in) :: islot
        integer             :: i

        if (amr_slot_live(islot)) return
        amr_slots(islot)%amr_ref_ratio = amr_ref_ratio
        amr_slots(islot)%buff_size = buff_size
        allocate (amr_slots(islot)%x_cb(-1:max_f1), amr_slots(islot)%x_cc(0:max_f1), amr_slots(islot)%dx(0:max_f1))
        if (n_glb > 0) allocate (amr_slots(islot)%y_cb(-1:max_f2), amr_slots(islot)%y_cc(0:max_f2), amr_slots(islot)%dy(0:max_f2))
        if (p_glb > 0) allocate (amr_slots(islot)%z_cb(-1:max_f3), amr_slots(islot)%z_cc(0:max_f3), amr_slots(islot)%dz(0:max_f3))
        @:ALLOCATE(amr_slots(islot)%q_cons(1:sys_size))
        @:ALLOCATE(amr_slots(islot)%q_cons_stor(1:sys_size))
        @:ALLOCATE(amr_slots(islot)%q_prim(1:sys_size))
        @:ALLOCATE(amr_slots(islot)%rhs(1:sys_size))
        @:ALLOCATE(amr_slots(islot)%q_ghost_a(1:sys_size))
        @:ALLOCATE(amr_slots(islot)%q_ghost_b(1:sys_size))
        do i = 1, sys_size
            @:ALLOCATE(amr_slots(islot)%q_cons(i)%sf(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi))
            @:ALLOCATE(amr_slots(islot)%q_cons_stor(i)%sf(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi))
            @:ALLOCATE(amr_slots(islot)%q_prim(i)%sf(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi))
            ! rhs is ghost-inclusive (mbuf); igr widens to -1:+1 per dim including collapsed ones (coarse rhs_vf is -1:m+1 etc.)
            if (igr) then
                @:ALLOCATE(amr_slots(islot)%rhs(i)%sf(mbuf1_lo:mbuf1_hi, min(mbuf2_lo, -1):max(mbuf2_hi, 1), min(mbuf3_lo, &
                           & -1):max(mbuf3_hi, 1)))
            else
                @:ALLOCATE(amr_slots(islot)%rhs(i)%sf(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi))
            end if
            @:ALLOCATE(amr_slots(islot)%q_ghost_a(i)%sf(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi))
            @:ALLOCATE(amr_slots(islot)%q_ghost_b(i)%sf(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi))
            @:ACC_SETUP_SFs(amr_slots(islot)%q_cons(i))
            @:ACC_SETUP_SFs(amr_slots(islot)%q_prim(i))
            @:ACC_SETUP_SFs(amr_slots(islot)%rhs(i))
            @:ACC_SETUP_SFs(amr_slots(islot)%q_cons_stor(i))
            @:ACC_SETUP_SFs(amr_slots(islot)%q_ghost_a(i))
            @:ACC_SETUP_SFs(amr_slots(islot)%q_ghost_b(i))
        end do
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

    end subroutine s_amr_alloc_slot

    !> Free slot islot's per-block field arrays (inverse of s_amr_alloc_slot). Idempotent (no-op if not live).
    impure subroutine s_amr_free_slot(islot)

        integer, intent(in) :: islot
        integer             :: i

        if (.not. amr_slot_live(islot)) return
        ! Undo each field's ACC_SETUP_SFs (Cray descriptor + %sf copyin) BEFORE the @:DEALLOCATE - Cray 'exit data delete'
        ! decrements
        ! the ref count, so the lone @:DEALLOCATE would leave the descriptor and the ACC_SETUP %sf ref dangling; the leaked host
        ! address is later reused (e.g. by Gs_rs at restart), tripping a Cray "Error placing / already present" present-table crash
        ! (gpu-acc).
        do i = 1, sys_size
            @:ACC_TEARDOWN_SFs(amr_slots(islot)%q_cons(i))
            @:DEALLOCATE(amr_slots(islot)%q_cons(i)%sf)
            @:ACC_TEARDOWN_SFs(amr_slots(islot)%q_cons_stor(i))
            @:DEALLOCATE(amr_slots(islot)%q_cons_stor(i)%sf)
            @:ACC_TEARDOWN_SFs(amr_slots(islot)%q_prim(i))
            @:DEALLOCATE(amr_slots(islot)%q_prim(i)%sf)
            @:ACC_TEARDOWN_SFs(amr_slots(islot)%rhs(i))
            @:DEALLOCATE(amr_slots(islot)%rhs(i)%sf)
            @:ACC_TEARDOWN_SFs(amr_slots(islot)%q_ghost_a(i))
            @:DEALLOCATE(amr_slots(islot)%q_ghost_a(i)%sf)
            @:ACC_TEARDOWN_SFs(amr_slots(islot)%q_ghost_b(i))
            @:DEALLOCATE(amr_slots(islot)%q_ghost_b(i)%sf)
        end do
        @:DEALLOCATE(amr_slots(islot)%q_cons)
        @:DEALLOCATE(amr_slots(islot)%q_cons_stor)
        @:DEALLOCATE(amr_slots(islot)%q_prim)
        @:DEALLOCATE(amr_slots(islot)%rhs)
        @:DEALLOCATE(amr_slots(islot)%q_ghost_a)
        @:DEALLOCATE(amr_slots(islot)%q_ghost_b)
        if (qbmm .and. .not. polytropic) then
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

        integer :: k
        logical :: needed

        do k = 1, amr_max_blocks
            needed = k <= amr_num_blocks
            if (needed) needed = amr_block_owner(k) == proc_rank
            if (needed) then
                call s_amr_alloc_slot(k)
            else
                call s_amr_free_slot(k)
            end if
        end do

    end subroutine s_amr_reconcile_slots

    impure subroutine s_finalize_amr_module()

        integer :: i, islot

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
        if (allocated(amr_seambuf_x)) deallocate (amr_seambuf_x, amr_seambuf_y)
        if (allocated(amr_seam_pairs)) deallocate (amr_seam_pairs)
        deallocate (amr_ovl_gather, amr_ovl_gather_n, amr_ovl_scatter, amr_ovl_scatter_n)
        do i = 1, sys_size
            @:DEALLOCATE(amr_cg(i)%sf)
        end do
        @:DEALLOCATE(amr_cg)
        deallocate (amr_decomp)
        deallocate (amr_slots)
        deallocate (amr_region_lo_all, amr_region_hi_all, amr_isect_lo_all, amr_isect_hi_all, amr_owns_all)
        if (allocated(sw_x_cb)) deallocate (sw_x_cb, sw_x_cc, sw_dx)
        if (allocated(sw_y_cb)) deallocate (sw_y_cb, sw_y_cc, sw_dy)
        if (allocated(sw_z_cb)) deallocate (sw_z_cb, sw_z_cc, sw_dz)
        if (allocated(amr_block_owner)) deallocate (amr_block_owner)
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
