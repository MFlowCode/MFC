# AMR fine-level distribution (design note)

**Status:** in progress on `up/mega` (part of the AMR effort). Phase 1 landed (map computed, not applied); Phase 2 is the apply + gather/scatter checkpoint.

## Problem

MFC's AMR uses a **mirror decomposition**: the fine level reuses the coarse rank
decomposition, so each rank owns exactly the fine cells over its coarse subdomain's
intersection with the blocks (`s_amr_compute_isect` → `amr_rank_owns_block`). This makes
coarse↔fine coupling entirely rank-local (zero communication) but leaves fine work
**unbalanced**: if refinement concentrates in a sub-region, only the ranks owning that
region do the fine advance while the rest idle. At scale this caps AMR speedup.

## Target

Give the fine level its **own distribution**, decoupled from the coarse decomposition, so
fine blocks spread across all ranks by measured work — the AMReX per-level
`DistributionMapping` idea. MFC already enforces an invariant that makes this far simpler
than the general case: **blocks are kept ≥ `buff_size` apart** (regrid merge), so with
**whole-block-per-rank** assignment there is *no fine–fine halo at all*
(`s_mpi_sendrecv_amr_fine_halo` disappears). The only coupling that becomes communication
is coarse↔fine, block-granular.

## Coupling changes (whole-block-per-rank)

| Step | Today (mirror, local) | Decoupled |
|---|---|---|
| Ownership | intersection, multi-owner (`amr_rank_owns_block`) | single owner (`amr_block_owner(k) == proc_rank`) |
| Ghost fill / prolong (coarse→fine) | read local coarse `q_cons` | **gather** coarse patch (block + `buff_size` halo) from coarse owner(s) to block owner |
| Reflux + restriction (fine→coarse) | write local coarse cells | **scatter** corrections from block owner back to coarse owner(s) |
| Fine–fine halo | coarse Cartesian neighbors | none (blocks separated) |
| Regrid | cluster (global) → mirror | cluster → **assign** (`amr_block_owner`) → **migrate** fine state |

The gather/scatter is the only new communication surface. SFC-ordered assignment keeps a
block's coarse patch on an SFC-nearby rank, bounding the cost. Conservation exactness is
preserved because the scatter applies the *same* reflux add / restrict overwrite the
mirror model applied locally — a pure data-movement change, not a numerics change.

## Reuse

`m_load_weight` (per-cell cost field) and `m_sfc_partition` (Morton order + chains-on-chains
balanced partition) already exist as an init-time diagnostic. They become the fine-level
distributor: feed the block list + per-block fine-work weight, return `amr_block_owner(:)`.

## Phasing (each independently mergeable)

1. **Distribution map, computed but not applied** (behavior-preserving). Add
   `amr_block_owner(:)` + SFC assignment + an imbalance diagnostic; mirror ownership and all
   coupling unchanged. Goldens bit-identical. *(this note's first increment)*
2. **Apply the map**: switch `amr_rank_owns_block` to single-owner, add the coarse↔fine
   gather/scatter, drop the fine–fine halo. Validate conservation + assignment-independence.
3. **Dynamic rebalance each regrid**; optionally drop the fixed max-size slot pool for
   right-sized boxes.

Phase 2 carries the correctness risk (gather/scatter exactness) and is the checkpoint gate.

## Phase 2 implementation status (branch `amr-fine-dist-wip`)

Done and np=1 bit-identical (the gather reduces to the local read when the owner is the sole rank):

- **Gather** (`s_amr_gather_coarse_patch`): per-block coarse patch `[region_lo-nmar : region_hi+nmar]`
  assembled by a sentinel-MAX allreduce into `amr_cg` (a drop-in `scalar_field`). Contribution rule
  `f_amr_own_coarse` claims interior cells + physical-boundary ghosts only (never inter-rank ghosts), so
  there is one authoritative contributor per cell and no coarse-ghost halo dependence. `pull_host` flag
  stages device-resident coarse to host for the runtime callers.
- **Read side rerouted**: init/regrid host prolongation (`s_interpolate_coarse_to_fine`) and the runtime
  device ghost-fill (`s_amr_fill_fine_ghosts`, stage + subcycle) now prolong from `amr_cg`.
- **Restriction scatter** (`s_restrict_fine_to_coarse`): owner restricts on host into the `amr_cg_wp`
  scratch (sole contributor), sentinel-MAX assembles, each rank overwrites the covered coarse cells it
  owns. `s_restrict_all_vars` (old local device kernel) deleted.
- **Fine–fine halo dropped**: `s_mpi_sendrecv_amr_fine_halo` calls removed (whole-block owner has no
  continuation faces; blocks ≥ buff_size apart). Routine in `m_mpi_proxy` now dead — remove with the
  reflux rework.

Remaining (the conservation crux — reflux is NOT yet decoupling-correct):

- **Reflux redesign.** `s_amr_reflux_face_flags` derives transverse participation `tv` from `amr_isect`,
  which is owner-only in the decoupled model, so `own_lo/own_hi` (gating BOTH `creg` capture and the
  reflux apply) is true only on the block owner — which does not own the outside coarse cells. Fix:
  (1) recompute `tv` from the replicated block range (`amr_region_lo/hi(t)` overlap with `sidx:sidx+ext`),
  so any coarse-outside-owner participates; (2) the owner captures `freg` (already correct) and BROADCASTS
  it (allreduce, replacing the cart-neighbor `s_mpi_sendrecv_amr_reflux_faces`); (3) coarse-outside-owners
  capture `creg` locally and apply `(freg-creg)` to their owned outside cells — verify the apply's
  block-relative transverse indexing maps onto each receiver's local coarse slice. Conservation-critical.
- **QBMM pb/mv** gather/scatter for np≥2 (currently local; np≥2 QBMM+AMR ungated, untested — gate or scatter).
- np=2 conservation + assignment-independence validation once reflux lands (acceptance goldens:
  BD21A5C0, 5EFB3277, 79B334C7 — they pass on `up/mega` via the mirror).

## Phase 2 status: DONE and validated (reflux + scalability)

The reflux redesign above landed (block-relative frame, participation from the replicated block range), and the
whole coupling is now **point-to-point + right-sized** — the fine level distributes with no global collective per
stage and per-rank memory that scales with the decomposition.

- **Reflux**: `s_amr_reflux_face_flags` derives `tv`/`tlo`/`thi` from the replicated block range; `creg` capture +
  both applies index block-relative; the owner's `freg` is delivered to the applying coarse-outside-owners.
- **Point-to-point coupling** (replaced the correctness-first global collectives): a replicated coarse-decomposition
  table `amr_decomp` (allgathered once at init) lets any rank compute which ranks hold a coarse-cell range.
  - gather: owner `Irecv`s patch slices from the coarse-owners (`f_amr_rank_coarse_range`).
  - restriction: owner `Isend`s each coarse-owner its covered slice (`f_amr_rank_interior`).
  - reflux: owner `Isend`s `freg` only to participants (`f_amr_reflux_participates` = `own_lo/hi` parameterized by
    rank via `amr_decomp`); `s_amr_p2p_reflux_faces` in `m_amr`.
  Non-participants exchange nothing.
- **Right-sized memory**: fine/coord/register arrays sized to `amr_maxc_fit` (min-over-ranks local half-extent = the
  max block a rank can own, enforced by the scratch-constraint abort), not the global half-domain `amr_maxc` — about
  `1/num_procs` (and `1/num_procs^(d-1)` for the face registers) the memory at scale.
- **Coordinate fix**: `s_amr_swap_to_fine` extends fine ghost coords from the global boundaries `amr_g?cb` (the owner
  need not hold the block's local coarse coord slice).

Validated: np=1 1D/2D/3D goldens bit-identical; np=2 conservation exact in 1D (energy 1.9e-16) and 2D with
transverse-face appliers (mass 1.4e-14), bounds-clean. `s_mpi_allreduce_array_max` remains only for the
once-at-init coordinate assembly.

## Phase 3

1. **Block-splitting / tiling** (the coverage gap) — DONE.
   - **Tiling**: `s_amr_tile_box` splits any box > `amr_maxc_fit` into contiguous ≤`amr_maxc_fit` sub-blocks, wired into
     the initial-block setup and the regrid pipeline (non-IB); `s_populate_amr_fine` loops over all blocks refreshing
     per-block mirrors via `s_amr_select_slot`.
   - **Block-to-block fine-fine halo**: the fine advance is three driver phases — **fill all** (gather + coarse
     ghost-fill), **`s_amr_fine_fine_halo`** (overwrite each seam ghost with the neighbour's *stage-entry* interior,
     buff_size-deep, `MPI_Sendrecv` between owners or a local copy when one rank owns both; adjacency from the
     replicated `amr_region_*_all` via `f_amr_seam`), **advance all** (RHS + RK). `s_advance_amr_fine_stage` split into
     `s_amr_fine_stage_fill` / `s_amr_fine_stage_advance`.
   - **Reflux seam-exclusion**: `f_amr_face_is_seam` drops a sub-block face shared with another fine sub-block from
     `s_amr_reflux_face_flags` `own_lo/own_hi` (fine-fine, not c/f).
   - **Owner-ordering fix**: `s_amr_assign_block_owners` runs BEFORE the owner-dependent `s_set_amr_fine_geometry` in
     init, regrid, AND restart — the stale default (rank-0) owner map otherwise sized multi-block owners wrong (only
     surfaces once tiling makes several blocks). Restart (`s_read_amr_restart`) is a two-pass read (regions → assign →
     place data) for both the parallel_io and non-parallel_io paths.
2. **Regrid cross-rank fine-state migration** — DONE (the stretched-np≥2 correctness fix). The regrid overlap-copy
   preserves a covering old block's fine detail by reading `amr_slots(kk)%q_cons_stor`, but that was local-only
   (`if (.not. old_owns(kk)) cycle`). Under fine-level distribution an old block can be owned by a rank *other* than the
   one now owning a covering new block, so the copy was silently skipped and the new block kept only its
   coarse-prolonged values — **exact on uniform grids** (prolongation reproduces the fine field) but **~O(1e-4) wrong on
   stretched grids** (prolongation is not exact). Fix: broadcast every old block's stashed fine state from its owner
   into the replicated slot's `q_cons_stor` on all ranks (`num_procs>1`; no-op at np=1), copy the overlap from *every*
   covering old block regardless of ownership, and take `old_ilo`/`old_ext` from the global replicated region (not the
   owner-only isect). Correctness-first collective; a per-block P2P version mirroring `s_amr_gather_coarse_patch` is
   future work. Validated (`-b mpirun`): 79B334C7 (1D stretched dynamic-regrid, 2 ranks) passes; F0DDE1B4 (np=1 twin),
   5EFB3277 + BD21A5C0 (uniform np=2 tiled + restart), and the np=1 AMR batch all pass.

## Known open

- **B7704247** (AMR 2D stretched-y dynamic regrid, np=1) fails at abs 4.4e-11 (> 1e-12 tol) on `cons.1.00`. Isolated:
  **not** tiling (persists with y-tiling disabled — the block stays whole), **not** the migration (np=1; migration is
  `num_procs>1`-guarded), and it fails at the committed state too. So it is a **whole-block 2D-stretched-y** regression
  from the earlier fine-distribution commits — the whole-block `amr_gycb` y-coordinates or the reflux redesign's y-face
  metric on a stretched grid. Small but real; a distinct, narrower investigation, and it may also affect np≥2
  2D-stretched.
- **QBMM pb/mv** gather/scatter for np≥2 (still local; np≥2 QBMM+AMR ungated/untested — gate or scatter). The regrid
  QBMM overlap-copy is still local-only (same bug class as the q_cons migration above, deferred with the rest of QBMM).
- Minor: the P2P routines scan all `num_procs` to find participants (integer-only, O(P) per block per stage); the regrid
  migration broadcasts (correctness-first, → P2P); patch-only device transfers instead of whole-field host round-trips
  (GPU only); own-only slot allocation.
