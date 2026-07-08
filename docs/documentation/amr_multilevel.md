# Multi-level AMR nesting — design and implementation plan

Status: **in progress.** The block-structured AMR core today supports a single refined
level (L1, 2:1 over the L0 base grid). This document is the roadmap for generalizing it to
arbitrary refinement depth (L0, L1, …, L`amr_max_level`). It is implemented in
behavior-preserving increments: at `amr_max_level = 1` the code is bit-identical to the
single-level core.

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
2. **Recursive coupling.** Add the per-block `%%level`/`%%parent` tags, then parameterize
   gather/prolong/restriction/reflux by a coarse-level source: the level-`l` block data instead
   of only `q_cons_base`. The L1↔L0 path is the `l = 0` special case and stays byte-identical.
3. **Recursive subcycle driver.** Generalize the `m_time_steppers` AMR advance into a
   recursive `advance(level)` with per-level `dt` and time-interpolated coarse ghosts.
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
