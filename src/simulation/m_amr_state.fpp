!>
!!@file
!!@brief Contains module m_amr_state

#! AMD OpenMP lane: assert allocatables present on every kernel here (see OMP_DEFAULT_STR). A kernel naming an unallocated
#! module array aborts, so conditionally allocated ones (sw_jac/jac: igr) are named only under their allocation condition.
#:set MFC_OMP_PRESENT_ALLOCATABLE = True
#:include 'macros.fpp'

!> @brief Module-scope AMR state (block pool, flat field store, wave buffers) and the leaf helpers every AMR module needs.
module m_amr_state

    use m_derived_types  ! scalar_field, t_box, int_bounds_info
    use m_global_parameters
    use m_amr_restart_io, only: s_amr_subdivide_cb
    use m_constants, only: mapCells
    use m_mpi_proxy, only: s_mpi_abort
    use m_mpi_common, only: s_mpi_allreduce_array_max
    use m_amr_registers, only: s_amr_parent_foot
    use m_bubbles_EL, only: s_lag_cloud_bbox_local

    implicit none

    !> One refined level: its own grid and metadata (host-only). Fields live in the flat stores amr_cons_st/amr_stor_st, indexed by
    !! amr_loc_of.
    type t_level
        integer               :: amr_ref_ratio
        type(t_box)           :: region   !< block extent in parent (level-0) cell indices
        integer               :: m, n, p  !< this level's interior extents
        integer               :: buff_size
        type(int_bounds_info) :: idwbuff(3)
        real(wp), allocatable :: x_cb(:), x_cc(:), dx(:)
        real(wp), allocatable :: y_cb(:), y_cc(:), dy(:)
        real(wp), allocatable :: z_cb(:), z_cc(:), dz(:)
    end type t_level

    !> Fixed pool of refined-block slots; the working slot amr_cur selects which one per-block routines operate on.
    type(t_level), allocatable :: amr_slots(:)

    !> Dense local index for live slots: amr_slots is sparse over 1:amr_max_blocks, the field store is dense. Freed indices are
    !! recycled through amr_loc_free.
    integer, allocatable :: amr_loc_of(:)      !< global slot -> dense local index, 0 if not live
    integer, allocatable :: amr_loc_free(:)    !< stack of recycled local indices
    integer              :: amr_loc_n = 0      !< high-water mark of local indices handed out
    integer              :: amr_loc_nfree = 0  !< depth of the recycle stack

    !> Flat per-block field store (x, y, z, var, local slot) so one kernel can run over every live block; all slots share the mbuf
    !! extents. Sized by s_amr_st_reserve.
    real(stp), allocatable, dimension(:,:,:,:,:) :: amr_cons_st, amr_stor_st
    $:GPU_DECLARE(create='[amr_cons_st, amr_stor_st]')
    !> local slots the store is sized for; grows only (s_amr_compact_store keeps the index space dense)
    integer :: amr_st_cap = 0

    !> Copy bridge to the shared solver, whose routines take type(scalar_field) and serve the monolithic path too (a pointer view
    !! into the flat store is not attachable on OpenMP offload): load from the store, call, store back, in both directions since the
    !! dummies are intent(inout). Not GPU_DECLARE'd: a declare binds the present-table entry to the descriptor that move_alloc then
    !! swaps out (CCE find_in_present_table failed); allocate via move_alloc + GPU_ENTER_DATA instead (CCE lib-4425 otherwise).
    type(scalar_field), allocatable :: amr_cons_br(:)
    !> Batched advance: the bridge spans amr_br_batch blocks along the last active dimension, member i at offset (i-1)*amr_bat_w
    !! with its own ghost shell, so no block's stencil reaches another's interior. amr_bat_loc = members' flat-store columns
    !! (device-resident).
    integer :: amr_bat_loc(amr_bat_max) = 0
    $:GPU_DECLARE(create='[amr_bat_loc]')
    !> Pooled advance scratch shared by every fine block (the fused advance leaves no cross-block lifetime); same pattern as amr_cg.
    type(scalar_field), allocatable :: amr_scr_prim(:), amr_scr_rhs(:)
    !> block-frame prim scratch for the batched per-member IB correction (ib only); the slab prim stacks members along amr_bat_sd
    type(scalar_field), allocatable :: amr_scr_prim_blk(:)
    !> Blocks per batched s_compute_rhs call (amr_bat_max once initialized); bounded to avoid exhausting device memory.
    integer :: amr_br_batch = 1
    integer :: amr_maxc(3)  !< max coarse block cells per dim: (m_glb+1)/2 etc.; 1 for collapsed dims

    !> Per-slot sizing: max fine cells per dim and the buffered array bounds. amr_slot_live(k): slot k's field arrays are allocated
    !! (lazy, owned-only).
    integer              :: max_f(3), mbuf_lo(3), mbuf_hi(3)
    logical, allocatable :: amr_slot_live(:)
    !> cached same-level seam list (xb, yb, seam-dim) so s_amr_fine_fine_post iterates O(#seams); rebuilt when amr_seam_pairs_dirty
    !! is set or the block count changes
    integer, allocatable :: amr_seam_pairs(:,:)
    integer              :: amr_num_seam_pairs, amr_seam_pairs_nblk
    logical              :: amr_seam_pairs_dirty
    ! amr_mesh_epoch lives in m_global_parameters (m_amr_registers keys on it and cannot use m_amr)
    !> Tag base of the clusterer's narrow-node reductions, above the per-box tag space [1..amr_max_blocks]; checked against
    !! MPI_TAG_UB at init.
    integer :: amr_tag_narrow = 0
    !> Keyed wave tags: amr_tag_base + band*65536 + gen*4096 + seq, checked against MPI_TAG_UB. gen (mod 16) bumps per wave entry on
    !! every rank; seq is the message's index in the pair's canonically ordered transfer list (block id, dim, lo before hi), so
    !! matching is independent of posting order.
    integer, parameter :: amr_band_reflux = 0, amr_band_freg = 1, amr_band_parent_fill = 2, amr_band_stage_fill = 3, &
        & amr_band_migrate = 4, amr_band_seam = 5, amr_band_restrict_l1 = 6, amr_band_restrict_parent = 7
    integer, parameter   :: amr_n_bands = 8
    integer              :: amr_tag_base = 0
    integer              :: amr_tag_gen(0:amr_n_bands - 1) = 0
    integer, allocatable :: amr_tsq(:,:)      !< (0:np-1, dir) in-wave per-peer seq counters; touched-reset
    integer, allocatable :: amr_tsq_tch(:,:)  !< touched peers per dir, so the reset is O(active peers), not O(P)
    integer              :: amr_n_tsq(2) = 0
    !> cached per-block P2P overlap-rank lists, rebuilt with the seam list: amr_ovl_gather(:,k) = ranks whose owned coarse range
    !! meets block k's padded patch box; amr_ovl_scatter(:,k) = ranks whose coarse interior meets k's region box. Rank-ascending and
    !! not owner-excluded, so iterating gives the same MPI order as a 0..num_procs-1 scan.
    !> (max-overlap, amr_max_blocks); sized in s_amr_build_seam_pairs
    integer, allocatable :: amr_ovl_gather(:,:), amr_ovl_scatter(:,:)
    !> Block indices this rank owns, ascending (the paired MPI_SENDRECVs depend on the order); built in s_amr_assign_block_owners.
    integer, allocatable :: amr_my_blk(:)
    integer              :: amr_n_my = 0
    !> Level-1 blocks this rank receives restriction for (not its own, overlapping its interior), rebuilt once per mesh epoch. Keyed
    !! on amr_mesh_epoch, not amr_myblk_dirty (owner writes only). Ownership here is amr_block_owner(k) == proc_rank, not the
    !! multi-owner amr_rank_owns_block.
    integer, allocatable :: amr_l1r_blk(:)
    integer, allocatable :: amr_l1p_blk(:)  !< padded variant (region +/- amr_cpat_mar vs my coarse range): stage-fill gather
    integer              :: amr_n_l1p = 0
    integer              :: amr_n_l1r = 0
    integer, allocatable :: amr_fch_blk(:)  !< foreign children of my owned parents (level >= 2, powner == me, cowner /= me)
    integer              :: amr_n_fch = 0
    integer, allocatable :: amr_own_blk(:)  !< blocks this rank intersects (amr_owns_all, the multi-owner notion, all levels)
    integer              :: amr_n_own = 0
    !> cached f_amr_parent_block (0 for level <= 1), avoiding an O(global blocks) scan per block per stage
    integer, allocatable :: amr_parent_blk(:)
    integer, allocatable :: amr_child_ptr(:)   !< children of p = amr_child_idx(amr_child_ptr(p-1)+1 : amr_child_ptr(p)), ascending
    integer, allocatable :: amr_child_idx(:)   !< (amr_child_ptr(0) = 0)
    integer(8)           :: amr_l1r_epoch = -1_8
    integer              :: amr_l1r_nblk = -1  !< second key half, matching s_amr_reg_prepare's (epoch, num_blocks) pair
    !> Set wherever amr_block_owner is written (assigner, inheritance, restart/migration), so lists built from it never go stale.
    logical              :: amr_myblk_dirty = .true.
    integer, allocatable :: amr_ovl_gather_n(:), amr_ovl_scatter_n(:)  !< per-block list lengths
    !> Wire pools of the pooled waves (m_amr_wave lays them out; contents never survive a wave). The fill and restrict waves share
    !! amr_fw_sq/rq; the seam wave is posted before the fills and drained after, so it has its own pair. The zero-copy reflux waves
    !! list their blocks in amr_fw_rblk.
    real(wp), allocatable :: amr_fw_sq(:), amr_fw_rq(:), amr_sw_sq(:), amr_sw_rq(:)
    real(wp), allocatable :: amr_rf_sq(:), amr_rf_rq(:)  !< device-resident pools of the packed reflux-faces wave
    integer, allocatable  :: amr_fw_rblk(:)
    !> Device-resident wire pools: with rdma_mpi MPI sends/receives them by device address, and the step-path kernels find them
    !! present, so slices do not cross PCIe per wave.
    logical :: amr_fw_dev = .false.
    !> seam pairs both of whose blocks this rank owns: exchanged by one batched kernel at drain, no wire
    integer              :: amr_sw_nsame = 0
    integer, allocatable :: amr_sw_plx(:), amr_sw_ply(:), amr_sw_pd(:), amr_sw_pxhi(:), amr_sw_pfm(:,:)
    !> Fused exchange packs (amr_device_pack): one row per wave transfer (slab corner (1:3), extents (4:6), payload offset in the
    !! wire pool (7); the child pack also has source store slot (8) and the child's patch frame (9:11)) plus the exclusive element
    !! prefix, so one kernel walks a family's transfers. Rebuilt per wave.
    integer, allocatable :: amr_fx_pl(:,:), amr_fx_pre(:)
    !> (0:num_procs-1, 1:amr_max_level) SFC Morton-key upper bound per rank, one independent cost-weighted cut per level (owner =
    !! cut search, f_amr_owner); O(P) storage. Per level because same-level boxes have distinct Morton keys, which the binary search
    !! needs.
    integer(kind=8), allocatable :: amr_fine_cut(:,:)

    !> Regrid box size cap per dim (fixed for the run, identical on all ranks; 1 in collapsed dims): a box of at most min-over-ranks
    !! of (local extent + 1)/2 cells intersects every rank in at most (its extent + 1)/2 cells, so the per-rank scratch constraint
    !! 2*(isect cells) - 1 <= local extent holds by construction. Equals amr_maxc at np=1.
    integer :: amr_maxc_fit(3) = 1

    !> Swap contract (s_amr_swap_to_fine / s_amr_restore_coarse): the fine advance swaps these coarse globals to the block's values
    !! and restores after: m/n/p, idwint/idwbuff, the nine coordinate arrays (sw_x_cb..sw_dz), acoustic_source; WENO/hypoelastic/IGR
    !! spacing coefficients are recomputed instead. Any new module-level state derived from m/n/p/idwint/idwbuff/coords that a
    !! kernel reads on the fine grid must be swapped here or refreshed per fine call, including its device copy (a stale device copy
    !! reads out of range under CCE OpenACC; cv_int is refreshed by a per-call GPU_UPDATE in s_compute_rhs). amr_swap_depth makes
    !! the swap re-entrant and guards an unpaired restore.
    !> Saved coarse-level global state for swap/restore
    integer               :: sw_m, sw_n, sw_p
    type(int_bounds_info) :: sw_idwint(3), sw_idwbuff(3)
    logical               :: sw_acoustic_source

    !> IGR sigma-state bounce: the fine solve reuses jac/jac_old at fine indices, so the coarse contents (jac_old is the Jacobi warm
    !! start) are saved here across the fine advance.
    real(wp), allocatable :: sw_jac(:,:,:), sw_jac_old(:,:,:)
    $:GPU_DECLARE(create='[sw_jac, sw_jac_old]')

    !> Swap nesting depth: only the outermost swap saves the coarse state and only its restore puts it back; inner swaps re-install,
    !! inner restores are no-ops.
    integer :: amr_swap_depth = 0
    !> Global coarse cell-boundary arrays (-1:X_glb), assembled once at init; the fine-distribution owner rebuilds whole-block fine
    !! coordinates from them.
    real(wp), allocatable :: amr_gxcb(:), amr_gycb(:), amr_gzcb(:)
    real(wp), allocatable :: sw_x_cb(:), sw_x_cc(:), sw_dx(:)
    real(wp), allocatable :: sw_y_cb(:), sw_y_cc(:), sw_dy(:)
    real(wp), allocatable :: sw_z_cb(:), sw_z_cc(:), sw_dz(:)

    !> True (identically on all ranks) iff some rank's fine ghost fill reads its coarse ghost cells: the solver populates only prim
    !! ghosts, so the cons ghosts must be halo-exchanged first. Never true at np=1.
    logical :: amr_xchg_coarse_ghosts = .false.
    !> local accumulator behind amr_xchg_coarse_ghosts; s_amr_reduce_xchg_flag does one allreduce per scan
    integer :: amr_xchg_bad = 0

    !> Per-block gathered coarse patch: the owner may not hold the coarse cells its block refines, so before each prolongation/ghost
    !! fill the patch region_lo-amr_cpat_mar : region_hi+amr_cpat_mar is gathered point-to-point into amr_cg (stp scalar_fields,
    !! block-local frame: cell 0 = global coarse cell amr_cpat_off(d)). Messages carry wp cast to stp, so at np=1 the patch equals
    !! the local coarse read bitwise. Sized to the largest block. amr_slab_tab: device-resident slab table for the shell/ghost
    !! kernels, rows lo(1:3),hi(1:3),soff,scnt over <= 6 slabs, refreshed by one GPU_UPDATE per launch.
    integer, allocatable            :: amr_slab_tab(:,:)
    integer                         :: amr_slab_n = 0  !< rows queued in amr_slab_tab (s_amr_slab_add)
    integer                         :: amr_cg_np = 1  !< patches amr_cg holds, stacked along dim 3 (s_amr_cg_reserve)
    integer                         :: amr_cg_nuse = 0  !< patches the level-1 pool gather filled (s_amr_l1_pool_gather)
    integer, allocatable            :: amr_cg_pz(:)  !< dim-3 offset of each owned level-1 box's patch in the pool
    type(scalar_field), allocatable :: amr_cg(:)
    integer                         :: amr_cpat_mar = 0  !< coarse-cell stencil reach (f_amr_cpat_mar)
    integer                         :: amr_cpat_hi(3) = 0  !< amr_cg upper local bounds per dim (0 in collapsed dims)
    integer                         :: amr_cpat_off(3) = 0  !< global coarse index of amr_cg local cell 0 (region_lo - amr_cpat_mar)
    integer                         :: amr_wcur = 1  !< consume cursor into the open fill wave's receive transfers (box-major)

    !> Subcycling patch cache: the pooled level-1 patches (amr_cg, amr_cg_nuse patches) at the step start (time level 0) then the
    !! step end (1), each flat with g1 fastest, then g2, g3, variable.
    real(stp), allocatable :: amr_pc(:)

contains

    !> Fill fcb/fcc/fdx by subdividing parent cells (s_amr_subdivide_cb); pcb_lb is lbound(pcb, 1). Only 0..nfine is filled.
    subroutine s_build_level_coords(pcb, pcb_lb, lo, nfine, fcb, fcc, fdx)

        real(wp), intent(in)                 :: pcb(:)
        integer, intent(in)                  :: pcb_lb, lo, nfine
        real(wp), allocatable, intent(inout) :: fcb(:), fcc(:), fdx(:)
        integer                              :: fi

        call s_amr_subdivide_cb(pcb, pcb_lb, lo, nfine, amr_slots(amr_cur)%amr_ref_ratio, fcb)
        do fi = 0, nfine
            fdx(fi) = fcb(fi) - fcb(fi - 1)
            fcc(fi) = 0.5_wp*(fcb(fi - 1) + fcb(fi))
        end do

    end subroutine s_build_level_coords

    !> Block k's region box in global coarse cells; collapsed dimensions read 0:0.
    pure subroutine s_amr_region_box(k, rlo, rhi)

        integer, intent(in)  :: k
        integer, intent(out) :: rlo(3), rhi(3)

        rlo = merge(amr_region_lo_all(:,k), 0, amr_dim); rhi = merge(amr_region_hi_all(:,k), 0, amr_dim)

    end subroutine s_amr_region_box

    !> Block k's padded coarse patch box [plo, phi] (region +/- amr_cpat_mar; collapsed dims 0:0).
    pure subroutine s_amr_patch_box(k, plo, phi)

        integer, intent(in)  :: k
        integer, intent(out) :: plo(3), phi(3)

        plo = merge(amr_region_lo_all(:,k) - amr_cpat_mar, 0, amr_dim); phi = merge(amr_region_hi_all(:,k) + amr_cpat_mar, 0, &
                    & amr_dim)

    end subroutine s_amr_patch_box

    !> How many parent cells inside a child's footprint the ghost prolongation reads: 1 for limited linear, 2 for the quartic.
    pure integer function f_amr_prolong_reach()

        f_amr_prolong_reach = merge(2, 1, amr_prolong == 2)

    end function f_amr_prolong_reach

    !> Coarse-patch margin amr_cpat_mar: the coarse cells a buff_size fine ghost layer spans plus the prolongation reach.
    pure integer function f_amr_cpat_mar()

        f_amr_cpat_mar = (buff_size + amr_ref_ratio - 1)/amr_ref_ratio + f_amr_prolong_reach()

    end function f_amr_cpat_mar

    !> Core [clo, chi] of the patch box that the ghost fill never reads, so the fill waves ship only the shell around it.
    pure subroutine s_amr_patch_core(k, clo, chi)

        integer, intent(in)  :: k
        integer, intent(out) :: clo(3), chi(3)

        clo = merge(amr_region_lo_all(:,k) + f_amr_prolong_reach(), 0, amr_dim)
        chi = merge(amr_region_hi_all(:,k) - f_amr_prolong_reach(), 0, amr_dim)

    end subroutine s_amr_patch_core

    !> Fine cell coordinates of block k in dimension d, rebuilt from the global coarse boundaries gcb by replaying k's ancestor
    !! chain from replicated metadata, so any rank can reproduce them. Do not bisect amr_slots(parent)%x_cb instead: the parent's
    !! slot exists only on its owner, which a level>=2 block's owner need not be.
    impure subroutine s_amr_build_block_coords(k, gcb, glb, fcb, fcc, fdx, d)

        integer, intent(in)                  :: k, d
        real(wp), intent(in)                 :: gcb(:)  !< global L0 cell boundaries
        integer, intent(in)                  :: glb     !< index of gcb's first element (-1, or -1 - buff_size for a padded build)
        real(wp), allocatable, intent(inout) :: fcb(:), fcc(:), fdx(:)
        integer                              :: chain(0:amr_max_level), lev, j, a, lo, nf, span, rr, plo(3), phi(3), clb
        real(wp), allocatable                :: cur(:), scb(:), scc(:), sdx(:)

        rr = amr_ref_ratio
        lev = amr_block_level(k)
        a = k
        do j = lev, 1, -1  ! chain(j) = k's ancestor at level j; chain(lev) = k
            chain(j) = a
            if (j > 1) a = f_amr_parent_block(a)
        end do

        ! level 0: the allocatable assignment takes the dummy's lbound (1), so cur's true first index travels in clb; levels below
        ! are built into (-1:nf) arrays
        cur = gcb; clb = glb
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
                call s_build_level_coords(cur, clb, lo, nf, fcb, fcc, fdx)
            else
                if (allocated(scb)) deallocate (scb, scc, sdx)
                allocate (scb(-1:nf), scc(0:nf), sdx(0:nf))
                call s_build_level_coords(cur, clb, lo, nf, scb, scc, sdx)
                cur = scb; clb = -1
            end if
        end do
        if (allocated(scb)) deallocate (scb, scc, sdx)

    end subroutine s_amr_build_block_coords

    !> Assemble the global coarse cell-boundary arrays over (-1-pad : G+pad): each rank writes its own cells' boundaries into a
    !! sentinel-filled array and an elementwise MAX allreduce recovers the global one. A larger-pad rebuild gives a consistent
    !! superset.
    impure subroutine s_amr_build_global_cb(pad)

        integer, intent(in) :: pad
        integer             :: j
        real(wp), parameter :: sentinel = -huge(1._wp)

        #:for D, X, E in [(1, 'x', 'm'), (2, 'y', 'n'), (3, 'z', 'p')]
            if (amr_dim(${D}$)) then
                if (allocated(amr_g${X}$cb)) deallocate (amr_g${X}$cb)
                allocate (amr_g${X}$cb(-1 - pad:${E}$_glb + pad)); amr_g${X}$cb = sentinel
                do j = -1 - pad, ${E}$ + pad
                    amr_g${X}$cb(start_idx(${D}$) + j) = ${X}$_cb(j)
                end do
                call s_mpi_allreduce_array_max(amr_g${X}$cb, ${E}$_glb + 2 + 2*pad)
            end if
        #:endfor

    end subroutine s_amr_build_global_cb

    !> Do two coarse-index boxes [alo:ahi] and [blo:bhi] overlap? Collapsed dims never disqualify.
    pure logical function f_amr_boxes_overlap(alo, ahi, blo, bhi) result(ov)

        integer, intent(in) :: alo(3), ahi(3), blo(3), bhi(3)

        ov = all((alo <= bhi .and. ahi >= blo) .or. .not. amr_dim)

    end function f_amr_boxes_overlap

    !> Index of the level-(level(k)-1) block that block k refines, or 0 for level 1 (parent is the L0 grid). Regions are in L0 cell
    !! indices at every level; proper nesting guarantees one overlapping parent.
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

    !> Pack (dir=+1) / unpack (dir=-1) the fine cells of slot's q_cons over [dlo:dhi] in dim d (interior transverse, all sys_size)
    !! in a fixed order, so a packer and unpacker with matching extents align cell-for-cell. On GPU only this near-seam slab moves
    !! device<->host.
    impure subroutine s_amr_fine_slice(slot, d, dlo, dhi, buf, dir)

        integer, intent(in)                 :: slot, d, dlo, dhi, dir
        real(wp), intent(inout), contiguous :: buf(:)
        integer                             :: i, a, b, c, fm(3), na, nb, nc, loc

        fm(1) = amr_slots(slot)%m; fm(2) = amr_slots(slot)%n; fm(3) = amr_slots(slot)%p
        loc = amr_loc_of(slot)
        nc = dhi - dlo + 1
        ! Pack on the device into the contiguous buf and move only buf: flang miscomputes a strided section of a block's field in a
        ! target-update map clause (corrupts 2D+ np>1 seam ghosts). buf runs a fastest, then b, c, i; na/nb = transverse sizes.
        #:for D, TA, TB, IDX in [(1, 2, 3, 'c, a, b'), (2, 1, 3, 'a, c, b'), (3, 1, 2, 'a, b, c')]
            if (d == ${D}$) then
                na = fm(${TA}$) + 1; nb = fm(${TB}$) + 1  ! scalars: kernel bounds use na-1/nb-1, not host fm(..) (nvfortran/Cray)
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

    !> Convert a physical-space bbox to a global coarse-index bbox padded by pad_cells (uniform spacing; the axisymmetric half axis
    !! cell only shrinks dy(0), so the floor stays conservative).
    pure subroutine s_amr_phys_to_cells(pmin, pmax, pad_cells, blo, bhi)

        real(wp), dimension(3), intent(in) :: pmin, pmax
        integer, intent(in)                :: pad_cells
        integer, intent(out)               :: blo(3), bhi(3)

        blo = 0; bhi = 0
        #:for D, X in [(1, 'x'), (2, 'y'), (3, 'z')]
            if (amr_dim(${D}$)) then
                blo(${D}$) = int((pmin(${D}$) - glb_bounds(${D}$)%beg)/d${X}$(${'min(1, n)' if D == 2 else 0}$)) - pad_cells
                bhi(${D}$) = int((pmax(${D}$) - glb_bounds(${D}$)%beg)/d${X}$(${'min(1, n)' if D == 2 else 0}$)) + pad_cells
            end if
        #:endfor

    end subroutine s_amr_phys_to_cells

    !> Rank-local per-stage guard: the local bubbles' padded bbox must stay clear of the current block (initial overlap, or drift
    !! beyond the regrid margin).
    impure subroutine s_amr_check_lag_clear()

        real(wp), dimension(3) :: pmin_loc, pmax_loc
        integer                :: blo(3), bhi(3)

        if (.not. bubbles_lagrange) return
        call s_lag_cloud_bbox_local(pmin_loc, pmax_loc)
        if (pmin_loc(1) > pmax_loc(1)) return  ! no bubbles on this rank
        call s_amr_phys_to_cells(pmin_loc, pmax_loc, mapCells + 2, blo, bhi)
        if (f_amr_boxes_overlap(blo, bhi, amr_slots(amr_cur)%region%lo, amr_slots(amr_cur)%region%hi)) then
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

    !> Minmod-limited linear value at fine offset (xix, xiy, xiz) in coarse cell u0, from its face neighbours (m, p) per direction;
    !! pass u0 for both neighbours of a direction without a slope.
    pure real(wp) function f_amr_lin(u0, xm, xp, ym, yp, zm, zp, xix, xiy, xiz) result(v)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: u0, xm, xp, ym, yp, zm, zp, xix, xiy, xiz

        v = u0 + minmod(xp - u0, u0 - xm)*xix + minmod(yp - u0, u0 - ym)*xiy + minmod(zp - u0, u0 - zm)*xiz

    end function f_amr_lin

    !> Admissible fine cell: positive density and (except model 4) positive internal energy; a NaN state is inadmissible.
    pure logical function f_amr_cons_ok(fi, fj, fk, loc) result(ok)

        $:GPU_ROUTINE(parallelism='[seq]')
        integer, intent(in) :: fi, fj, fk, loc
        real(wp)            :: rho, ke
        integer             :: v

        rho = 0._wp; ke = 0._wp
        do v = eqn_idx%cont%beg, eqn_idx%cont%end
            rho = rho + real(amr_cons_st(fi, fj, fk, v, loc), wp)
        end do
        ok = rho > 0._wp
        if (ok .and. model_eqns /= 4) then
            do v = eqn_idx%mom%beg, eqn_idx%mom%end
                ke = ke + real(amr_cons_st(fi, fj, fk, v, loc), wp)**2
            end do
            ok = real(amr_cons_st(fi, fj, fk, eqn_idx%E, loc), wp) - 0.5_wp*ke/rho > 0._wp
        end if

    end function f_amr_cons_ok

    !> Average over one half (xi < 0 left, > 0 right; |xi| = 1/4) of the middle of five cell averages q: the conservative quartic,
    !! or MC-limited linear where either half leaves the range of the three nearest averages widened by their smallest second
    !! difference. Both halves are tested so the two children share one reconstruction (their mean stays q(0)).
    pure real(wp) function f_amr_quartic_half(q, xi) result(v)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: q(-2:2), xi
        real(wp)             :: c, a, b, d, lo, hi

        d = (3._wp*(q(-2) - q(2)) - 22._wp*(q(-1) - q(1)))/128._wp
        c = min(abs(q(-2) - 2._wp*q(-1) + q(0)), abs(q(-1) - 2._wp*q(0) + q(1)), abs(q(0) - 2._wp*q(1) + q(2)))
        lo = min(q(-1), q(0), q(1)) - c; hi = max(q(-1), q(0), q(1)) + c
        if (min(q(0) - abs(d), q(0) + abs(d)) >= lo .and. q(0) + abs(d) <= hi) then
            v = q(0) + sign(1._wp, xi)*d
        else
            a = q(1) - q(0); b = q(0) - q(-1)
            v = q(0)
            if (a*b > 0._wp) v = q(0) + sign(min(0.5_wp*abs(a + b), 2._wp*min(abs(a), abs(b))), a)*xi
        end if

    end function f_amr_quartic_half

    !> Allocate the pooled q_prim/rhs advance scratch (idempotent). The lockstep driver associates it for every block, so it must
    !! exist on every rank, including one that never allocates a slot. rhs gets the igr widening.
    impure subroutine s_amr_scr_init()

        integer                         :: i, slo(3), shi(3)
        type(scalar_field), allocatable :: tmp_p(:), tmp_r(:)

        if (allocated(amr_scr_prim)) return
        ! blocks are stacked along the last active dimension (see amr_cons_br)
        slo = [mbuf_lo(1), mbuf_lo(2), mbuf_lo(3)]; shi = [mbuf_hi(1), mbuf_hi(2), mbuf_hi(3)]
        shi(num_dims) = slo(num_dims) + amr_br_batch*(shi(num_dims) - slo(num_dims) + 1) - 1
        ! CCE OpenMP-offload: a bare module-scope scalar_field array aborts on allocate (lib-4425); allocate a local, move_alloc.
        allocate (tmp_p(1:sys_size)); call move_alloc(tmp_p, amr_scr_prim)
        allocate (tmp_r(1:sys_size)); call move_alloc(tmp_r, amr_scr_rhs)
        $:GPU_ENTER_DATA(create='[amr_scr_prim, amr_scr_rhs]')
        do i = 1, sys_size
            @:ALLOCATE(amr_scr_prim(i)%sf(mbuf_lo(1):shi(1), mbuf_lo(2):shi(2), mbuf_lo(3):shi(3)))
            if (igr) then
                @:ALLOCATE(amr_scr_rhs(i)%sf(mbuf_lo(1):shi(1), min(mbuf_lo(2), -1):max(shi(2), 1), min(mbuf_lo(3), &
                           & -1):max(shi(3), 1)))
            else
                @:ALLOCATE(amr_scr_rhs(i)%sf(mbuf_lo(1):shi(1), mbuf_lo(2):shi(2), mbuf_lo(3):shi(3)))
            end if
            @:ACC_SETUP_SFs(amr_scr_prim(i))
            @:ACC_SETUP_SFs(amr_scr_rhs(i))
        end do
        if (ib) then
            allocate (tmp_p(1:sys_size)); call move_alloc(tmp_p, amr_scr_prim_blk)
            $:GPU_ENTER_DATA(create='[amr_scr_prim_blk]')
            do i = 1, sys_size
                @:ALLOCATE(amr_scr_prim_blk(i)%sf(mbuf_lo(1):mbuf_hi(1), mbuf_lo(2):mbuf_hi(2), mbuf_lo(3):mbuf_hi(3)))
                @:ACC_SETUP_SFs(amr_scr_prim_blk(i))
            end do
        end if

    end subroutine s_amr_scr_init

    !> Grow-on-demand sizing for append-only tables and pools: a grow preserves the entries already appended.
    impure subroutine s_amr_size_int(a, n)

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

    end subroutine s_amr_size_int

    impure subroutine s_amr_size_int3(a, n)

        integer, allocatable, intent(inout) :: a(:,:)
        integer, intent(in)                 :: n

        call s_amr_size_rows(a, 3, n)

    end subroutine s_amr_size_int3

    !> s_amr_size_int for a table of n columns of r ints.
    impure subroutine s_amr_size_rows(a, r, n)

        integer, allocatable, intent(inout) :: a(:,:)
        integer, intent(in)                 :: r, n
        integer, allocatable                :: tmp(:,:)

        if (.not. allocated(a)) then
            allocate (a(r, max(n, 64)))
            return
        end if
        if (size(a, 2) >= n) return
        call move_alloc(a, tmp)
        allocate (a(r, max(n, 2*size(tmp, 2))))
        a(:,1:size(tmp, 2)) = tmp

    end subroutine s_amr_size_rows

    !> Wire pools, preserving on grow. dev keeps the pool device-resident across (re)allocation: the old image is deleted from the
    !! device before it is freed and the new one created after; contents never survive a wave, so nothing is copied.
    impure subroutine s_amr_size_real(a, n, dev)

        real(wp), allocatable, intent(inout) :: a(:)
        integer, intent(in)                  :: n
        logical, intent(in)                  :: dev
        real(wp), allocatable                :: tmp(:)

        if (.not. allocated(a)) then
            allocate (a(max(n, 64)))
            if (dev) then
                $:GPU_ENTER_DATA(create='[a]')
            end if
            return
        end if
        if (size(a) >= n) return
        if (dev) then
            $:GPU_EXIT_DATA(delete='[a]')
        end if
        call move_alloc(a, tmp)
        allocate (a(max(n, 2*size(tmp))))
        a(1:size(tmp)) = tmp
        if (dev) then
            $:GPU_ENTER_DATA(create='[a]')
        end if

    end subroutine s_amr_size_real

end module m_amr_state
