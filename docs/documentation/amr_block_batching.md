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

**Confirmed 2026-08-02 with AMR ACTIVE (the sweeps above are uniform-only / one fine block).**
2D 4096x2048, np=2, corner-refined so each rank already advances 16 fine blocks, varying only
`l0_ntile`. The right independent variable is TOTAL block advances per rank, \f$I = 16 + T\f$
with \f$T = \texttt{l0\_ntile}^{2}\f$ — not the tile count alone:

| `l0_ntile` | \f$T\f$ | advances \f$I\f$ | \f$I/I(1)\f$ | measured | agreement |
|---|---|---|---|---|---|
| 0 | 1 | 17 | 1.000 | 1.000 | — |
| 2 | 4 | 20 | 1.176 | 1.362 | +16% |
| 3 | 9 | 25 | 1.471 | 1.610 | +9% |
| 4 | 16 | 32 | **1.882** | **1.932** | **+3%** |

Cost tracks total block advances, converging to 3% at the largest point. The excess at small
\f$T\f$ is expected: an L0 tile starts far larger than a fine block and costs more per advance,
converging as it shrinks toward fine-block size — i.e. per-advance cost is *weakly* size-dependent,
not perfectly flat.

> **Recorded because I got it wrong.** Fitting the same law on \f$T\f$ ALONE — ignoring the 16
> fine-block advances already in the step — predicts 2.77x at \f$T=16\f$ against 1.932x measured, a
> 30% overprediction, and briefly looked like a falsification of the per-block floor itself. The floor
> is real; the fit used the wrong variable. Any future sweep of this kind must count EVERY block advance
> in the step, not just the ones being varied.

> **SCOPE LIMIT added 2026-08-01: the launch-count premise below holds only at LOW block counts.**
> Profiling a 320-block 3D case (`rocprofv3`, kernels + memory copies, np=2) puts kernel execution at
> 8.5% of the span, device-to-device copies at 23.5% (2937729 of them, 20.5 per kernel), and the GPU
> idle ~68%. Mean kernel duration is 69 us, and 143k launches at 10 us each is ~1.4 s against ~82 s of
> idle - launch count is not in the top two cost terms at that scale. The measurements in this section
> were taken on a 16-tile / 1-fine-block case and remain valid there, but **ranking increments by
> "launches removed" does not generalize to the many-block regime this arc exists to fix.** See
> @ref amr_per_level_distribution, queue item 14, including its sizing caveat.

## Direct confirmation: the cost is kernel launch count

Counted with `rocprofv3 --kernel-trace` (2D 128^2, 6 steps, np=1, under `srun`):

| | monolithic | 16 blocks (`l0_ntile=4`) | ratio |
|---|---|---|---|
| GPU kernel launches | 381 | 7619 | **20.0x** |
| wall time (256^2 sweep) | — | — | 16.75x |

Launch count scales 20x with block count while wall time scales 16.75x, so the per-block
cost tracks launches directly. This is measured, not inferred: combined with the refuted
hoist below (removing swap traffic changed nothing), it closes the argument that the fixed
per-block cost IS the launch count. Note each tile advance issues ~450 launches against the
monolithic step's 381 — the per-block path adds its own ghost fills and halo work on top of
the same RHS kernel sequence.

### Where the launches go, and what rides on them

Repeated at 2048x1024 with dynamic regrid (6 steps, np=1), tracing memory copies as well —
note `--kernel-trace` **alone leaves the copy table empty**, so pass `--memory-copy-trace`
too or the transfers below are invisible:

| | uniform | AMR | ratio |
|---|---|---|---|
| kernel launches | 381 | 4619 | 12.1x |
| device-to-device copies | 7,094 | 88,667 | 12.5x |
| kernel busy time | 297 ms | 601 ms | 2.0x |
| wall span | 716 ms | 4264 ms | **10.9x** |

Wall time tracks launches (12.1x vs 10.9x), not kernel work (2.0x): about 86% of the AMR span
is GPU idle between kernels. Splitting AMR's launches by whether the kernel also appears in
the uniform run:

- **shared-solver kernels** — the per-block `s_compute_rhs` re-invocation — **3165 launches,
  565 ms: 69% of launches and 94% of kernel time**;
- AMR-only kernels (ghost fills, seam halo, RK update, flux capture, reflux): 1454 launches
  but only 36 ms, **6%** of kernel time.

So further fusion of the AMR-side kernels is worth much less than its launch share suggests.
The lever is reducing `s_compute_rhs` invocations.

The ~89k copies are **not** a separable AMR inefficiency: normalised per launch they are 18.6
(uniform) and 19.2 (AMR) — kernel-argument marshalling that the uniform solver pays too.
There is nothing to batch there. This reconciles the refuted hoist below rather than
contradicting it: transfers are not independently expensive, they are the mechanism by which
launch count costs time. Cutting `s_compute_rhs` calls cuts launches and their copy traffic
together.

### Regrid is 7-10% of runtime, not half

An earlier revision of this note inferred "regrid costs about half the runtime, roughly one
time step per regrid" from an A/B on `amr_regrid_int` (0.2974 s -> 0.1525 s per step at np=4
going from every-2-steps to every-10). **That inference is wrong and is retracted.** Varying
`amr_regrid_int` does not hold the block set fixed: regridding less often leaves staler,
generally smaller boxes, so the A/B prices the blocks a regrid *produces* together with the
regrid itself. Per-phase `system_clock` instrumentation inside `s_amr_regrid` measures the
cost directly (2048x1024, 20 steps, `amr_regrid_int=2`, so 10 regrid calls):

| | np=1 | np=4 | np=8 |
|---|---|---|---|
| full rebuilds / early-outs | 2 / 8 | 2 / 8 | 2 / 8 |
| total regrid time | 0.955 s | 0.459 s | 0.602 s |
| total run time | 13.0 s | 5.8 s | 6.2 s |
| **regrid share of runtime** | **7.3%** | **7.9%** | **9.7%** |

Most regrid calls are cheap: 8 of 10 find the box set unchanged and return after tagging
(~0.028 s at np=1, ~0.007 s at np=8), so only the two full rebuilds cost anything. Regrid is
therefore **not** comparable to the block advance, which remains the dominant lever rather than
one of two equal halves. (This sentence originally named the *packed super-grid* as that lever.
Packing is DISPROVED below - the lever is the per-invocation cost of the advance itself, which
packing cannot remove.)

Where a full rebuild's time actually goes differs with rank count:

- **np=1** is the host round trip of fine block state: pulling each old block to the host and
  bouncing it through `q_cons_stor` (0.159 s), the host-side `s_interpolate_coarse_to_fine`
  (0.042-0.091 s) and overlap copy (0.066 s), and pushing the new state back (0.043-0.146 s).
  All of it is serial host work on data that is already resident on the device.
- **np>=4** is dominated by *per-box synchronization*, not by any collective's own cost. The
  per-box `MPI_ALLREDUCE` at the end of `s_set_amr_fine_geometry` measures 7.4 ms (np=4) to
  13 ms (np=8) per call, ~700x a real one-integer allreduce — because it is absorbing the
  spread in the owner-only `s_amr_alloc_slot` work that precedes it. Inserting an
  `MPI_BARRIER` at the top of the rebuild loop collapses that phase from 0.089 s to 0.0011 s
  (np=4) and 0.119 s to 0.0001 s (np=8), with the barrier picking up the difference.
  **Batching those allreduces into one buys nothing AT THIS BOX COUNT** (14-21 boxes, np<=8):
  the time is the wait, not the reduction, so hoisting only moves the rendezvous. The rebuild
  loop would have to stop synchronizing per box at all to help here.

  **That conclusion is regime-limited and does not carry to scale.** With `nboxes` collectives the
  cost has a floor of `nboxes x latency` independent of any imbalance - at 10^5 boxes that is ~0.5 s
  of pure rendezvous per regrid, before a single byte of imbalance. The reduction was hoisted anyway
  (@ref amr_per_level_distribution, "What actually blocks exascale") because it removes an O(nboxes)
  term from the assignment; **expect no measurable change at the sizes benchmarked in this note**,
  which is what this paragraph predicts.

### The box set is not rank-invariant

Same case, same steps, varying only rank count: the regrid produces 14 boxes at np=1, 2, and
4, but **21 at np=8**. `amr_maxc_fit` is the min over ranks of the local half-extent, so more
ranks shrink the cap until boxes that fit at np<=4 must be tiled. Total boxes rise 50% while
per-rank boxes fall only 25% (3.5 -> 2.6), and every per-box collective in the rebuild loop
runs over *all* boxes on *all* ranks. This is a concrete, previously unattributed mechanism
for the np=8 turnover, and it ties that turnover to the same base-subdomain scratch cap — not to
MPI collective volume. (Originally "the cap the packed super-grid has to address"; packing is
DISPROVED below, and the cap was instead addressed directly by pinning `amr_max_grid_size`, which
is what the "Resolved" note immediately following records.)

**Resolved (2026-08-01): that was the DERIVED cap, and with `amr_max_grid_size` pinned the AMR
machinery is exactly rank-invariant.** Measured on a 3D sharp-interface case (128^3, two slabs at
a density ratio of 8, cap 32) at np = 1, 2, 4, 8: the level-1 tag sets are **byte-identical** —
31752 tags, zero symmetric difference in the tagged cell coordinates — and that carries through to
16 level-1 boxes and `fine_work` = 691200 at every rank count. The tagger, the `MPI_ALLGATHERV` tag
union, the coarse-CONS halo exchange, and the signature clusterer are all rank-invariant as designed,
including at np=1 where active-box windowing is active.

**But a badly conditioned CASE will still make the box set look rank-dependent, and that is not a
code defect.** The 2D benchmark pair (`hg8_*`, `sc_*`) uses the `hcid 299` multi-mode analytic
perturbation, which leaves a large fraction of cells sitting marginally at `amr_tag_eps` = 0.02.
Any floating-point-level difference then flips them wholesale: at np = 1/2/4/8 the level-1 tag
counts were 26849 / 37674 / 53924 / 61924 and the tag SETS were nearly disjoint (np=2 vs np=8 shared
**zero** cells; np=1 vs np=2 shared 77 of ~27000). Only 1.5–6.7% of tags lay on internal rank-seam
planes, so seam effects were a minor contributor — the bulk was threshold marginality across the
whole domain. Consequence: `fine_work` varied 5.4x with rank count (35.0M at np=1 to 188.8M at np=8),
which makes any AMR strong-scaling number from that case meaningless.

**Benchmark rule that follows.** An AMR scaling case must be verified rank-invariant BEFORE it is
timed — compare `fine_work` and the per-level box counts across the rank counts to be used, and
reject the case if they move. Prefer a sharp discontinuity (unambiguous tagging) over a smooth
perturbation. The literature is explicit that AMR strong scaling is hard to interpret for exactly
this reason (Athena++: work "is highly variable and depends on the refinement criteria"; PLUTO:
strong scaling is "difficult to interpret in the case of AMR computations"), and the standard metric
is zone-cycles/s — cell-updates summed over ALL levels, divided by wall time — which is what
`fine_work` + base cells provides here. A second, complementary option is to freeze the hierarchy
(`amr_regrid_int = 0`, valid up to `amr_max_level = 2`), which holds work constant by construction
and isolates solver and halo scaling from AMR-decision variability.

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
| RHS working set | module allocatables across `m_weno`, `m_rhs` (including its `vector_field` set, which dominates), and `m_riemann_solvers`, all sized to the base subdomain. `m_viscous` holds no volumetric scratch and does not contribute. |

The RHS working set being base-subdomain-sized is also why a block is capped at about half the
subdomain per dimension: the fine advance borrows the rank-local solver scratch. The enforced
cap is `amr_maxc_fit` (min over ranks of the local half-extent), **not** `amr_maxc` — the
latter is the global half and is read nowhere outside its own computation. For non-IB a box
exceeding the cap is **tiled**, not rejected; that tiling is one of the ways the many-blocks
regime this note addresses is produced. Only IB, which must own a body's block whole and
un-tiled, aborts.

## Why this gates the other arcs

The same cap is what blocks multi-level coexist: relaxing the `amr_max_level > 1` clause of
the `m_checker.fpp` `l0_ntile > 0 .and. amr` gate produces
`the nested level-2 block exceeds the per-rank scratch cap (2*L0-extent > amr_maxc_fit);
static multi-level does not tile the level-2 block`. So the coexist gate and this note's
redesign are blocked on the same single-working-slot architecture, not on independent work.
Sequencing the batching arc *after* the unification arc therefore risks reworking it.

## Strong-scaling baseline (the "before" curve)

Track 3's evidence artifact is a strong-scaling curve. The baseline, measured on the same
machine, 2D 256^2, static single-level AMR, 20 steps:

| np | compact block (32x32 coarse) | wide block (128x32 coarse) |
|---|---|---|
| 1 | 0.1009 s (1.00x) | 0.1009 s (1.00x) |
| 2 | 0.0980 s (1.03x) | 0.1190 s (0.85x) |
| 4 | 0.0983 s (1.03x) | 0.1267 s (**0.80x**) |

**AMR does not strong-scale, and wide blocks anti-scale.** Both cases hold exactly one fine
block (a static block cannot exceed `amr_maxc_fit`, so tiling never triggers), and single-owner
distribution puts all of its work on one rank whatever `np` is. The compact case is
therefore flat: extra ranks only split the coarse grid, which is not the bottleneck. The
wide case is worse than flat because the block spans more ranks' coarse subdomains, so each
added rank buys more coarse<->fine P2P gather/scatter with no fine parallelism to offset it.

Caveat: this measures the single-block regime. It does not test multi-block distribution —
reaching several blocks needs dynamic regrid, since a static block is capped at
`amr_maxc_fit`.

### The multi-block curve

Measured with dynamic regrid and many disjoint tag clusters, 2D 2048x1024, 20 steps, median
of the solver's steady `Time/step` (**not** `Time Avg`, which is a running mean still
carrying startup and therefore a function of run length):

| np | AMR | speedup | uniform control | speedup |
|---|---|---|---|---|
| 1 | 0.6634 s | 1.00x | 0.0611 s | 1.00x |
| 2 | 0.3971 s | 1.67x (84%) | 0.0509 s | 1.20x (60%) |
| 4 | 0.2974 s | 2.23x (56%) | 0.0463 s | 1.32x (33%) |
| 8 | 0.3720 s | 1.78x (22%) | 0.0427 s | 1.43x (18%) |

This **supersedes the single-block conclusion above**: with several blocks to distribute, AMR
does strong-scale, and scales better than the uniform control through np=4. It then turns
over at np=8, a 25% regression the uniform arm does not show.

The number batching must move is the np=1 column: **AMR costs 10.9x the uniform solver on the
same grid**, and the ratio grows with problem size (5.7x at 256x128). That is the per-block
serial advance.

Two cautions for anyone repeating this. Run the **uniform control** arm — without it a rising
AMR curve cannot be told apart from a problem too small to scale, and at 256x128 AMR appears
to *anti*-scale (+65% from np=1 to np=8) purely because fixed per-step overhead dominates a
32k-cell workload. And measure at a realistic size: the small case gives a qualitatively
wrong answer.

## Swap topology (and what blocks the obvious optimization)

The hot path has exactly **one `s_amr_swap_to_fine` / `s_amr_restore_coarse` pair per block
per RK stage**, wrapping `s_compute_rhs` in `s_amr_fine_stage_rhs` (and its subcycle twin).
The RK pass (`s_amr_fine_stage_rk`) does not swap — it works on slot arrays at slot extents.
The fill phase does not swap either; it reads the gathered patch `amr_cg` and slot arrays.
So a timestep costs `3 * nblocks` swap round-trips, each carrying
`s_amr_sync_grid_state_to_device` plus, on nonuniform grids, `s_amr_recompute_weno_coefs`.

The obvious cheap win — hoist the restore, letting consecutive blocks swap fine->fine and
restoring once per phase — has **two** blockers, not one.

1. *Reflux in the coarse frame.* `m_time_steppers` called `s_amr_p2p_reflux_faces` and
   `s_amr_apply_reflux` between blocks inside the advance loop, both operating on the coarse
   `rhs_vf` at coarse indices. **Resolved:** the advance and reflux loops are now split
   (phase 3 advances all blocks, phase 4 refluxes), so the reflux runs in the coarse frame
   after the advances rather than interleaved with them.
2. *Nested swap sites inside the RK pass.* Still open, and the larger of the two.
   `s_amr_fine_stage_rk` reads no grid globals itself (slot arrays plus `dt` and feature
   flags), but it calls `s_amr_pressure_relax_fine`, `s_amr_ib_correct_fine`, and
   `s_amr_update_mib_fine`, and **each opens its own `s_amr_swap_to_fine` /
   `s_amr_restore_coarse` pair**. They assume the coarse frame on entry. With the restore
   deferred they would swap while already swapped and trip
   `@:ASSERT(.not. amr_swapped, "nested s_amr_swap_to_fine (swap/restore must pair)")`.

   **Resolved:** the swap is now depth-counted and re-entrant — only the outermost swap
   saves into the `sw_*` bounce buffers, only the outermost restore puts them back, and an
   inner swap re-installs the same slot idempotently. No nesting occurs yet, so depth never
   exceeds 1 and behavior is unchanged.
3. *The IGR sigma bounce is save-and-seed in one routine.* Still open. `s_amr_swap_to_fine`
   ends with `if (igr) call s_amr_igr_swap_sigma()`, which BOTH saves the coarse
   `jac`/`jac_old` into `sw_jac`/`sw_jac_old` AND seeds the block's sigma from the parent —
   and it reads `sw_idwbuff` for the coarse extent. Once nesting is actually used, an inner
   swap would re-run it and overwrite `sw_jac` with fine state.
   **Resolved:** the save loop is depth-guarded, the seed loop is not. The seed reads
   `sw_jac`, which still holds the coarse state, so every nested block seeds from the correct
   parent — guarding only the save is both necessary and sufficient.

   Note this routine already caused an OpenACC-only crash (its own comment: the `sw_*`
   host-only module state "makes OpenACC's present lookup fail (OpenMP's implicit map(to)
   tolerates it, which is why only acc lanes crashed)"). It is a worked example of why an
   OpenMP-offload pass does not validate this area, and it is the reason the Frontier acc
   lane is the real gate for anything that changes how this routine nests.

With all three resolved, hoisting the restore became a small change, and it was **built,
measured, and abandoned**. Record of that, so it is not retried:

A `s_amr_swap_hold` entry point (save the coarse state without installing a slot) was added
and the per-tile RHS loop in `s_l0_advance_stage_rhs` was bracketed with it, so consecutive
tiles went fine->fine with one coarse restore after the loop instead of one per tile. Output
stayed byte-identical. The `l0_ntile` sweep, 2D 256^2, np=1:

| `l0_ntile` | tiles | with hoist | without |
|---|---|---|---|
| 0 | monolithic | 1.00x | 1.00x |
| 1 | 1 | 1.01x | 1.01x |
| 2 | 4 | 4.32x | 4.40x |
| 4 | 16 | **16.87x** | **16.75x** |

**No effect.** The coarse round trip between blocks — its `sw_*` copies, its device sync, its
WENO/FD coefficient rebuild — is not a measurable share of the per-block cost, even at 16
blocks where the total penalty is ~17x. The change was reverted rather than carried as
unused machinery.

This sharpens the conclusion above: the per-block cost really is **kernel launch count**
inside `s_compute_rhs`, not swap traffic. Reducing swap traffic in any form (per-slot
coefficient tables, per-slot device grid state, hoisted restores) should be expected to do
nothing on its own. Only increment 3 — one launch over a block list instead of one launch
set per block — addresses the measured cost.

The three prerequisites remain resolved in the code, so a future batching attempt does not
have to redo them; note they currently have no caller and are unused generality.

Note the SWAP CONTRACT warning in `.claude/rules/common-pitfalls.md`: a stale device copy of
coarse bounds reads out of range on the fine grid under **CCE OpenACC only**. A CPU-only or
NVHPC-acc pass proves nothing for this class, so any swap-contract change needs a Cray GPU
run before it can be trusted.

## Increments

### 0. Same-rank seam exchange — LANDED (`2e7151fa`)

`s_amr_fine_seam_exchange` replaced the four `s_amr_fine_slice` calls in the same-rank branch of
`s_amr_fine_fine_halo` with ONE fused device kernel doing both directions. Byte-identical.
Measured (2D 64x32, `l0_ntile=4` = 16 tiles / 24 seam pairs, MI250X): launches 10561 -> 7969
(`rocprofv3`), Time Avg 0.8179 -> 0.7121 s, n=3 per arm, **ranges disjoint**; run-to-run variance
also collapsed 3.3% -> 0.5%.

**Two rules this increment established, which govern everything below:**

- **Rank increments by LAUNCHES REMOVED. Credit transfer savings at zero.** The intermediate step
  (device-to-device but still 2 kernels/pair) was worth **+0.03%** — removing four blocking
  device<->host round trips per pair per stage bought *nothing*, while going 2 -> 1 kernels bought
  13%. On this machine launch count is the cost and small-slab PCIe traffic is free.
- **Measure with `rocprofv3 --kernel-trace`, not wall time.** Back-to-back runs scatter 0.17%; runs
  separated in time scatter 3%+. A 1.7% "win" was once reported here that was pure noise.

### Why this arc is now the main line

Per-level distribution is done and its cost model is fixed (@ref amr_per_level_distribution, steps
4-7). Measured imbalance is now 1.012 at np=4 — 98.8% of perfect — while parallel efficiency on the
same run is only 62%. **Balance has been eliminated as the limiter by fixing it**, and the residual
is the per-block fixed overhead measured below. This arc is what remains.

### DISPROVED (2026-08-01): the packed super-grid cannot work, because blocks are already slot-sized

The packed super-grid — lay \f$P\f$ same-level blocks contiguously along x in one slot and call the
unmodified `s_compute_rhs` once — is **arithmetically impossible for tiled blocks**, at every level
and every rank count. It is not a scheduling or memory-budget problem; there is no \f$P>1\f$ to reach.

Each packed block keeps its full buffered extent, so \f$P\f$ blocks of fine extent \f$f\f$ need
\f$P(f + 2\,\texttt{buff\_size}) - 2\,\texttt{buff\_size} \le \texttt{max\_f1} + 1\f$. Writing the
block's coarse extent as \f$b_c\f$ and the cap as \f$C\f$ (so \f$\texttt{max\_f1}+1 = 2C\f$):

\f[ P_{\max} = \left\lfloor \frac{C - b_c}{b_c + \texttt{buff\_size}} \right\rfloor + 1 \f]

so \f$P>1\f$ requires \f$b_c \lesssim C/2\f$. But `s_amr_tile_box` splits **evenly**:
\f$n_{tl} = \lceil e/t_c \rceil\f$ then \f$s = \lceil e/n_{tl}\rceil\f$, which for \f$n_{tl}\ge 2\f$
gives \f$s > t_c(1 - 1/n_{tl}) \ge t_c/2\f$. **Every tiled block therefore exceeds half the tile
size, and \f$P_{\max} = 1\f$ identically.** The tiler's purpose is to make blocks as large as the
slot permits, so a slot sized to hold one maximal block can never hold two.

Measured on the 75.5M case (2D, np=8, `amr_max_blocks` = 4096) with a temporary probe printing each
level's block extents against `max_f1` — level-1 blocks come out *exactly* slot-sized at every cap:

| `amr_max_grid_size` | `max_f1` | level-1 `fx_min`/`fx_max` | \f$P\f$ | blocks packable 2-up |
|---|---|---|---|---|
| 128  | 255  | 240 / 256   | 1 | 0 of 4096 |
| 256  | 511  | 496 / 512   | 1 | 0 of 1152 |
| 512  | 1023 | 1008 / 1024 | 1 | 0 of 288 |
| 1024 | 2047 | 2032 / 2048 | 1 | 0 of 72 |

Level 2 leaves 2–3 blocks per run small enough to pack, out of 128–1156 — under 0.3%.

This also retires the "packing is 7x at np=1 but collapses to ~2x at np>=4" table. That table was
never about rank count in the way it read: it came from a case (2047x1023, 7 stripes) whose blocks
were 108 coarse cells against a 1024 cap, i.e. *untiled* clustered boxes far below the cap — a
regime that does not survive a pinned cap at production scale. Separately, the rank-dependence it
worried about is genuinely gone: `amr_max_grid_size > 0` makes both the pack slot
(`amr_maxc_fit(d) = min(amr_maxc(d), amr_max_grid_size)`) and the solver scratch
(`idwbuff_alloc`, `m/n/p_alloc`) independent of the decomposition. \f$P_{\max}=1\f$ regardless.

Reviving packing would require a pack buffer and solver scratch sized \f$P\times\f$ a maximal block —
a genuinely new allocation, which is exactly the cost the design's "no new allocation is needed"
scope cut claimed to avoid. That cut relied on the np=1 coincidence \f$\texttt{max\_f1} = m_{glb}\f$,
which the pinned cap removes by construction.

### The lever that does exist: block size, bounded by device memory

Since per-block cost is near-fixed, the way to spend it on fewer blocks is to make each block bigger —
which is the `amr_max_grid_size` parameter, not new code. Same case and ranks, varying only the cap
(@ref amr_per_level_distribution, "The per-block floor"): 70.3 s/step at cap 128 versus 8.5–9.5 at
cap 1024, a ~20x span, *while the cap-1024 arm advances 2.8x more cells*.

That does not extend indefinitely. Per-rank scratch is \f$O(C^{\,\texttt{num\_dims}})\f$, and at cap
2048 this case aborts inside `__tgt_target_data_begin_mapper` — device out of memory — before
completing a single step. The failure mode is worth knowing: one rank core-dumps and the rest block
until the job is killed, so it presents as a hang. **The efficient regime is the largest cap that
fits device memory**, which for this case at np=8 in 2D lies between 1024 and 2048.

Even at the best measured cap the residual is ~7.3x uniform per cell (11.4 vs 1.55 ns/cell-update),
so per-block overhead remains the target — but it must be attacked by making each block's advance
cheaper (launch fusion, increments 1–3 below), not by combining blocks.

### Remaining increments

1. **Per-slot derived tables.** Give each slot its own WENO/FD coefficient storage so a swap
   selects rather than rebuilds. Bit-identical; measurable only on stretched/axisymmetric
   grids (`amr_weno_coef_recompute`), where it also drops ~18 device transfers per swap.
2. **Per-slot device grid state.** Make the device copies of `m/n/p`, `idwint/idwbuff`, and
   the coordinate arrays per-slot so `s_amr_sync_grid_state_to_device` becomes an index
   change rather than a transfer.
3. **Batched block kernels.** Replace the per-block launch of each RHS/RK kernel with one
   launch over a block list, with the per-block geometry read from the slot arrays produced
   by increments 1 and 2. This is the increment that actually removes the measured cost, and
   the one that retires the `amr_swapped` paired-swap guard. With the packed super-grid
   disproved above, this is again the *only* route to batching — and it still needs the flat
   backing store it always did, since `amr_slots` is not `GPU_DECLARE`'d and a runtime slot
   index inside a kernel is a null dereference. Cost and risk are unchanged from the original
   assessment: `ACC_SETUP_SFs` would perform overlapping partial mappings of one shared array,
   which is undefined for the present table and untestable outside Cray. Settle that with a
   build before committing to the increment.

Re-measure with the `l0_ntile` sweep above after each increment: it is cheap, byte-identity
checked, and needs no new instrumentation. Two harness traps: the LAST `Time Avg` line in a
run log is a "Saving" line reporting 0.0 (parse the `Time step` lines), and `D/` is only
populated by post_process (checksum `restart_data/` instead).
