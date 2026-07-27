@page amr_block_batching AMR per-rank block batching

# AMR per-rank block batching (design note)

> **Design record / implementation note.** This documents the measured cost of the swap-based
> per-block advance and the design it implies. For user-facing behavior and parameters, see
> @ref amr. For the fine-level distribution across ranks, see @ref amr_fine_distribution.

Status: **not implemented.** This note records the measurement, the state inventory, and the
increment plan, so the work can be picked up without re-deriving any of it.

## Problem

Within a rank, owned blocks advance **sequentially**: `s_amr_swap_to_fine` overwrites the
solver's global grid state with one block's, the block advances, `s_amr_restore_coarse` puts
the coarse state back, and the next block repeats. Per-rank wall time therefore scales with
the *sum* of its blocks' work, so strong scaling saturates exactly when the fine-level
distribution succeeds at giving ranks many small blocks.

## Measurement

`l0_ntile` tiles the base grid into `l0_ntile**d` blocks advanced through the *same*
swap-based per-block solver, with a correctness bar of byte-identical to `l0_ntile = 0`.
Grid, flops, and answer are fixed; only the number of swap/advance cycles varies. So
`cost(ntile) - cost(monolithic)` isolates the per-block overhead exactly.

Measured on MI250X (amdflang AFAR 23.2.0, OpenMP offload), 2D uniform grid, np=1, 30 steps.
Byte-identity verified at every point (restart md5 identical across tile counts).

| `l0_ntile` | blocks | TimeAvg (s) | vs monolithic |
|---|---|---|---|
| 0 | monolithic | 0.04502 | 1.00x |
| 1 | 1 | 0.04527 | 1.01x |
| 2 | 4 | 0.19811 | 4.40x |
| 4 | 16 | 0.75420 | 16.75x |

cost/blocks = 1.01, 1.10, 1.05: **linear in block count, each block advance costing about as
much as advancing the whole grid** while holding 1/16 the cells.

It does not amortize with problem size:

| N | monolithic (s) | 16 blocks (s) | ratio |
|---|---|---|---|
| 128^2 | 0.0427 | 0.7675 | 17.96x |
| 256^2 | 0.0450 | 0.7542 | 16.75x |
| 512^2 | 0.0511 | 0.7843 | 15.36x |

The 16-block time is flat (~0.77 s) across a 16x increase in cells, and so is the monolithic
time. Both are fixed-overhead bound at these sizes; the GPU is far from saturated. The
penalty is still 15x at 512^2.

**Conclusion: per-block cost is independent of block size.** It is fixed overhead per block
per RK stage, not work — dominated by per-block kernel launch count plus the per-swap device
syncs, none of which shrink as blocks shrink.

## What this rules out

@ref amr says "per-slot state instead of the global swap" is the lever. That is necessary but
**not sufficient**: it removes the per-swap state traffic, not the per-block launch count.
Two specific traps:

- Caching the swap-recomputed coefficient tables per slot removes the *host* recompute but
  **not** the device transfers, because the WENO kernels read module-global arrays. Removing
  the transfers requires the device-side arrays to be per-slot and the kernels to index by
  slot.
- Running K blocks concurrently in separate lanes multiplies the O(N^d) working set by K
  without reducing launch count.

The lever is **batching many blocks into single kernel launches over a block list**, not
per-lane duplication.

## State inventory

Per-block *data* is already per-slot (`t_level`: `region`, `m/n/p`, `buff_size`, `idwbuff`,
the nine coordinate arrays, `q_cons`, `q_cons_stor`, `q_prim`, `rhs`, `q_ghost_a/b`, and the
QBMM side-state). What is still global, and must be indexed by block for a batched kernel:

| Category | Items |
|---|---|
| Solver geometry (the SWAP CONTRACT) | `m/n/p`, `idwint`, `idwbuff`, nine coordinate arrays, `acoustic_source`, `ab_active` |
| Derived per-grid tables | WENO coefficients (`poly_coef_*`, `d_cb*`, `beta_coef_*`), hypoelastic FD coefficients, IGR `jac`/`jac_old` (bounced via `sw_jac`) |
| Scratch justified by sequential advance | `amr_cg` + `amr_cpat_off/hi`, `amr_rvw`, `amr_rhs_pb_f`/`amr_rhs_mv_f` (the last is explicitly commented "shared across slots (slots advance sequentially)") |
| RHS working set | 24 module allocatables across `m_weno` (16), `m_rhs` (5), `m_riemann_solvers` (2), `m_viscous` (1), sized to the base subdomain |

The RHS working set being base-subdomain-sized is also why `amr_maxc` caps a block at about
half the subdomain per dimension: the fine advance borrows the rank-local solver scratch.

## Why this gates the other arcs

The same cap is what blocks multi-level coexist: relaxing the `amr_max_level > 1` clause of
the `m_checker.fpp` `l0_ntile > 0 .and. amr` gate produces
`the nested level-2 block exceeds the per-rank scratch cap (2*L0-extent > amr_maxc_fit);
static multi-level does not tile the level-2 block`. So the coexist gate and this note's
redesign are blocked on the same single-working-slot architecture, not on independent work.
Sequencing the batching arc *after* the unification arc therefore risks reworking it.

## Increments

1. **Per-slot derived tables.** Give each slot its own WENO/FD coefficient storage so a swap
   selects rather than rebuilds. Bit-identical; measurable only on stretched/axisymmetric
   grids (`amr_weno_coef_recompute`), where it also drops ~18 device transfers per swap.
2. **Per-slot device grid state.** Make the device copies of `m/n/p`, `idwint/idwbuff`, and
   the coordinate arrays per-slot so `s_amr_sync_grid_state_to_device` becomes an index
   change rather than a transfer.
3. **Batched block kernels.** Replace the per-block launch of each RHS/RK kernel with one
   launch over a block list, with the per-block geometry read from the slot arrays produced
   by increments 1 and 2. This is the increment that actually removes the measured cost, and
   the one that retires the `amr_swapped` paired-swap guard.

Re-measure with the `l0_ntile` sweep above after each increment: it is cheap, byte-identity
checked, and needs no new instrumentation. Two harness traps: the LAST `Time Avg` line in a
run log is a "Saving" line reporting 0.0 (parse the `Time step` lines), and `D/` is only
populated by post_process (checksum `restart_data/` instead).
