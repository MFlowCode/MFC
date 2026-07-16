@page amr_multilevel Multi-level AMR nesting

# Multi-level AMR nesting — design and implementation plan

> **Design record / implementation note.** This documents the internal design and development of
> multi-level AMR nesting. For user-facing behavior and parameters, see @ref amr.

Status: **multi-level nesting implemented.** The block-structured AMR core supports arbitrary
refinement depth (L0, L1, …, L`amr_max_level`, 2:1 per level): static AMR (`amr_regrid_int = 0`)
nests one level-2 block, dynamic regrid (`amr_regrid_int > 0`) nests deeper and per-level. This
document is the design record. It was implemented in behavior-preserving increments: at
`amr_max_level = 1` the code is bit-identical to the single-level core.

## The one assumption to generalize

Everywhere in `m_amr`, "coarse" means the **L0 base grid**:

- `t_level%%region` is a box in **L0 cell indices**.
- Every coupling routine reads/writes the L0 fields (`q_cons_base` / `q_cons_ts(1)`):
  gather (`s_amr_gather_coarse_patch`), prolong (`s_interpolate_coarse_to_fine`),
  restriction (`s_restrict_fine_to_coarse`), reflux (fine flux corrects L0).
- The advance driver (`m_time_steppers`) loops the single fine-block pool once per L0 step.

Multi-level replaces "coarse = L0" with "**the coarse side of level `l+1` is level `l`**".

## Target architecture

- Levels `0 .. amr_max_level`. A level-`l` block (`l ≥ 1`) refines a covering level-`(l-1)`
  region by `ref_ratio` (2 today). Level `l+1` must be **properly nested** inside level `l`
  (surrounded by level-`l` cells; never adjacent to level `l-1`).
- **Flat pool + per-block level tag** (chosen over a per-level pool array): `amr_slots` stays
  one pool; each block gains `level` and `parent` (the covering coarser block, or 0 for an
  L1 block whose parent is L0). Rationale: the distribution machinery (SFC owners, P2P
  gather/restriction/reflux/migration, lazy owned-only slot allocation, repartition-on-restart)
  all operate on the flat pool and are level-agnostic, so they carry over with near-zero change.
  A per-level pool array would force rewriting every `amr_slots(k)` reference (high churn,
  silent-bug-prone — the same reason the #4 lazy-alloc avoided a global→local pool remap).
- **Recursive subcycling** (chosen over lock-step): `advance(l)` advances level `l` by `dt_l`,
  then for `s = 1..ref_ratio` recursively advances level `l+1` by `dt_{l+1} = dt_l/ref_ratio`
  with time-interpolated C/F boundary data from level `l` (extends today's 2-level subcycle,
  `q_ghost_a/b`, recursively), then refluxes `l+1 → l`. Most accurate and efficient; the
  standard Berger–Colella time integration.

## Increments (each validated np=1 bit-identical + np≥2 conservation + GPU)

1. **Foundation (bit-identical).** `amr_max_level` namelist param (default 1), gated
   `amr_max_level > 1` fail-closed in `m_checker` until the recursion lands. *(this increment)*
2. **Recursive coupling.** Add the per-block `%%level`/`%%parent` tags (done: `amr_block_level`
   in 2a; `f_amr_parent_block(k)` finds the covering coarser block by region overlap), then
   parameterize gather/prolong/restriction/reflux by a coarse-level source: the level-`l` block
   data instead of only `q_cons_base`. The L1↔L0 path is the `l = 0` case and stays byte-identical.

   **The one hard detail — the coarse frame.** The prolong reads `amr_cg` via
   (`amr_isect_lo`, `amr_cpat_off`, `ref_ratio`) and is reused unchanged as long as `amr_cg` and
   that frame are in **parent-level cell indices**. Level-1 block: parent level is L0, frame is
   L0 indices (today). Level-2 block with L0 region `R2`, parent L1 block `p` with L0 region `R1`:
   - parent-fine index of L0 cell `c` is `2*(c - R1.lo)` (ref_ratio per level);
   - coarse footprint in the parent's fine frame: `isect = 2*(R2 - R1.lo)`;
   - fine extent `m = ref_ratio*(parent-fine footprint) - 1 = ref_ratio^2 * |R2| - 1`;
   - `amr_cg` holds the parent's fine cells `[isect.lo - nmar : isect.hi + nmar]`, gathered from
     `amr_slots(p)%%q_cons` (np=1: local copy; np≥2: P2P via `amr_block_owner`);
   - "coarse coords" are the parent's fine coords `amr_slots(p)%%x_cb`, not `amr_gxcb` — needed
     only for the advance/stretched grids, so a **uniform-grid np=1 operator self-check**
     (`prolong → restrict` conserves) is the first testable milestone and can skip coords.
   `s_set_amr_fine_geometry` and `s_amr_gather_coarse_patch` choose this frame; each gains a
   `level == 1 ? L0 : parent-fine` branch.
3. **Advance (make level-2 evolve).** 2b built the coupling *into* a level-2 block
   (gather-from-parent, prolong); this builds the coupling *out* of it plus the driver:
   - **restrict-to-parent** — level-aware `s_restrict_fine_to_coarse` folds a level-2 block's
     fine averages back into its parent L1 block's fine array (mirror of gather-from-parent,
     same parent-fine frame);
   - **reflux-to-parent** — the Berger-Colella C/F flux correction from L2 into L1's cells
     (the reflux registers key off "the coarse", which becomes level l-1);
   - a **persistent static L2 block** (replacing the non-intrusive self-test), and a
     **level-loop driver** in `m_time_steppers`: for level 1..maxlevel, fill+advance from the
     parent then restrict+reflux to it.
   Both restrict-to-parent and reflux-to-parent are required for conservation.
   **Lock-step first** (all levels advance at L0's `dt`, interleaved per RK stage — extends
   today's non-subcycle mode; no time-interpolation/recursion): first milestone is a np=1
   static 2-level case that runs several steps and conserves (~1e-13). **Then** add recursive
   subcycling (level l+1 takes `ref_ratio` substeps with time-interpolated ghosts) on top,
   generalizing `s_advance_amr_fine_substeps` (already Berger-Colella for L0↔L1).
4. **Per-level regrid + proper nesting.** Tag each level for the next-finer one; build level
   `l+1` boxes clustered + tiled + nested inside level `l`; distribute (reusing the SFC map).
5. **Restart / distribution / GPU per level.** Extend the fine-restart record with a per-block
   level; validate repartition-on-restart and the device present-table across levels.

## Design decisions on record

- Flat pool + per-block `level`/`parent` tag (not a per-level pool array).
- Recursive subcycling with time-interpolated C/F BCs (not lock-step).
- `ref_ratio` stays 2 for now (ratio-4 is separate, banked work); nesting/coupling are written
  ratio-generic so ratio-4 drops in later.
- Load balance: the existing SFC map distributes blocks across **all** levels from the flat
  pool; the `amr_max_blocks < num_procs` warning applies per the total block count.
