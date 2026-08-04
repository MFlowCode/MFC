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

### MEASURED 2026-08-02: the per-block tax is DUMMY ARGUMENTS, not launch count or descriptors

Attributing 337,431 device copies to 16,698 dispatches (rocprofv3, temporal attribution, np=1 2D AMR):

| kernel | dispatches | copies/dispatch |
|---|---|---|
| `m_amr::fine_rk_update` | 678 | 54.0 |
| `m_riemann_solver_hllc` | 696 | 44.0 |
| `m_rhs::compute_rhs` | 696 | 37.7 |
| **`m_amr_registers::capture_boundary_flux`** | **1356** | **0.0** |
| **`m_amr_registers::apply_reflux`** | 240 | **0.0** |

Kernels reading MODULE-LEVEL `GPU_DECLARE`'d state (`freg`/`creg`, `m_amr_registers.fpp:56-57`) pay
ZERO; kernels taking fields as DUMMY ARGUMENTS pay 33-54. Copy count matches the deep-member count
exactly (`3*sys_size + q_T_sf + bc_type(2,2)` = 20 at `sys_size=5`, measured 20.2/dispatch). Total copy
time 2.0 s against ~1.5 s of kernel time - a first-order cost.

A 4-translation-unit reproducer (`amr-bench/attach/mtu/`) isolates the penalty at **4.3-4.7x per
region**, including per-advance re-attach, and confirms cross-TU `declare target` residency, runtime
slot indexing on device, and alternating attach targets all work.

**Two hypotheses this KILLS.** Flattening the `scalar_field` interface (3,840 `%%sf` sites across
`common/`) buys nothing - flat-array dummies measure the same as deep-allocatable ones. And batching
attacks dispatch COUNT while the tax is per-dispatch argument mapping.

**Conversion attempted on `s_amr_fine_rk_update` and REVERTED at 74/76.** Four failures, in order:
`GPU_EXIT_DATA(detach=)` emitted `map(always,from:)` (copies back, drops residency -> NaNs);
`GPU_ENTER_DATA(attach=)` emitted `map(always,to:)` (clobbers live device state with stale host values
-> wrong answers); allocating the views behind `.not. amr` (pure-L0 reaches the routine via
`s_l0_advance_stage_rk` -> segfault). The first two were real MFC bugs in clauses with ZERO in-tree
users, now fixed to `map(always,alloc:)` / `map(release:)`.

**THE OPEN BLOCKER, isolated by a control experiment: there is no correct ATTACH primitive for this
backend.** Keeping the attach calls but reverting the kernel to read its dummy arguments fails
*byte-identically* (`icfl Inf` on `00EB793A`, `8.357478479233075E+22` on `EF58E377`), so the kernel
change is innocent and **the attach operation itself corrupts the slot's device data**:

- `map(always,to:)` copies the STALE HOST array over live device state.
- `map(always,alloc:)` forces a fresh, uninitialised device allocation, orphaning the live data.

A device-side probe comparing the view against the dummy inside one kernel reads exactly zero
difference, which is consistent rather than contradictory: both resolve to the same *wrong* buffer.
That probe is worth keeping in mind - "the view resolves correctly" and "the view resolves to the right
memory" are different questions.

Only the multi-level cases fail because they regrid *and* advance most often; the static multi-level
np=2 golden `09E0D257` passes. A standalone 6-phase reproducer
(`amr-bench/attach/mtu/f_lifetime.f90`: attach-only, attach+detach, slot recycling, varying extents,
realloc-at-new-extent) passes every phase, so this is an INTEGRATION property, not a property of the
technique - an MWE de-risks the mechanism, not the integration.

Before retrying: establish an attach primitive that neither copies nor reallocates on already-present
data (candidates: `map(to:)` / `map(alloc:)` WITHOUT `always`, or `omp_target_associate_ptr`), and
validate it against `00EB793A` specifically. Note `GPU_ENTER_DATA(attach=)` has 8 existing users in
`m_rhs.fpp`/`m_igr.fpp`, so its expansion cannot be changed casually - they currently rely on the
`always,to` behaviour.

### HOW THE FIELD SOLVES THIS, and the plan to match it (2026-08-03)

This is a named, solved problem elsewhere, and the consensus answer is the one this note originally
proposed as increment 3.

- **AMReX** fuses across boxes: *"If there are 512 patches of 32^3 cells each, only one GPU kernel is
  launched to work on all 512 patches, which enables it to achieve similar performance as if operating
  on a single patch of 256^3 cells."* They quantify the problem on OUR hardware: *"a simple kernel
  running on an AMD MI250X will only achieve ~10% of its peak memory bandwidth on small boxes of 32^3
  cells."* They apply the same fusion to halo pack/unpack, where *"the dominant cost was kernel launch
  latency"*, reducing it to ONE launch per rank. (arXiv 2403.12179)
- **Parthenon** does the same via `MeshBlockPack`, with the pack size *"hardware and problem dependent,
  and so may be set at runtime."* (arXiv 2202.12309)
- **Castro/AMReX gridding guidance**: *"Best performance is obtained with bigger boxes, so setting
  `amr.max_grid_size = 128` and `amr.blocking_factor = 32` can give good performance"*; *"too small
  max_grid_size may ruin the code performance."* This independently corroborates the cap finding.
- **OpenMP offload specifically** carries higher per-region overhead than CUDA/HIP because of its device
  runtime layer; LLVM offers a reduced "bare-metal" kernel mode. MFC pays more per operation than an
  equivalent HIP code would, which raises the value of reducing operation COUNT.

**MEASURED HERE (`amr-bench/attach/serial/batch.f90`), the same mechanism, same arithmetic, 64 blocks:**

| | span | operations | copies/kernel |
|---|---|---|---|
| one launch PER BLOCK (MFC today) | 433.5 ms | 15370 | 0.0 |
| one launch over ALL BLOCKS (AMReX/Parthenon) | **35.0 ms** | **250** | 0.0 |

**12.4x, with byte-identical results.** Note BOTH rows show 0 copies per kernel because the store is a
PLAIN CONTIGUOUS ARRAY. That resolves a contradiction earlier in this note: the per-region map traffic
is not about module-scope versus dummy argument, it is about **plain array versus derived type with a
pointer component**. Module-scope `scalar_field` measured 32.8 copies/kernel; a plain module array
measures 0. AMReX's MultiFab is precisely the plain-contiguous-across-boxes layout.

**PLAN TO MATCH THEM.** Three pieces, in dependency order, each golden-gated:

1. **Flat contiguous backing store** - `real(stp) :: store(cell, var, block)` for the per-block field
   families, replacing per-slot `scalar_field` allocations. This is the enabler for 2 and it removes
   the descriptor traffic by itself. Sizing must come from max live blocks per rank, NOT
   `amr_max_blocks` (1024 would OOM) - **step 1a landed this as `amr_loc_of`/`amr_loc_n`, byte-identical
   across 76 goldens (`905e9d7c`)**.

   **SCOPING CONSTRAINT, found by attempting it: a field family cannot be migrated independently if it
   shares a consumer routine with a family that is not migrating.** `q_ghost_a`/`q_ghost_b` look like
   the smallest possible first conversion - 9 references, one file, two consumers - but
   `s_amr_fill_fine_ghosts` is not ghost-specific: it is a general prolong-coarse-into-fine routine
   whose THIRD caller targets `q_cons` (`m_amr.fpp:4597`). Branching inside the routine does not help,
   because a dummy referenced in ANY branch of a target region is still mapped, so the migrated path
   would keep paying the tax and the conversion would buy nothing. **The natural conversion unit is
   therefore `{q_cons, q_ghost_a, q_ghost_b}` together.** Counting references to a family NAME
   understates its blast radius; count the consumers and then the consumers' other targets.

   Separately confirmed while scoping: `q_ghost_a/b` are allocated UNCONDITIONALLY while their only
   consumers sit on subcycle-only paths (`s_amr_subcycle_setup_block`, `s_amr_advance_fine_subcycle_all`,
   `s_amr_advance_children`), and their `pb`/`mv` twins already carry an `if (amr_subcycle)` guard.
   Adding the matching guard frees roughly `2 x sys_size x mbuf x live_slots` of device memory on every
   non-subcycle run - memory that competes directly with `amr_max_grid_size`, the ~20x lever.
2. **Batched block kernels** - one launch over the block list instead of per block, for the RHS and RK
   paths. This is the 12.4x above and it is what `amr_max_grid_size` is currently substituting for.
3. **Fused halo pack/unpack** - AMReX reduced this to one launch per rank; MFC's seam/ghost fills are
   the same shape of work.

Target, from the 60%-of-GCD profile: busy is 8.2 s of a 28.4 s span, so eliminating serialization is
worth **~3.5x on that case** before counting the bandwidth gain AMReX reports from larger effective
kernels. Expect more where blocks are smaller, since the penalty scales with block count.

**RETRACTION.** Earlier in this session I wrote that batching "attacks dispatch count, not the
per-dispatch tax" and set it aside. The profile then showed dispatch count IS the governing quantity,
and AMReX and Parthenon both converge on exactly this mechanism. Increment 3 below was right and the
dismissal was wrong. Note the separate "packing is disproved" result concerns MFC's SLOT-based packing
into one working slot - a different mechanism from cross-block fusion, and only the former is dead.

### THE GOVERNING LAW (2026-08-03, profiled): wall time is set by GPU OPERATION COUNT

Gap analysis of a real AMR run (`rocprofv3`, kernel + copy traces, np=1 2D, 6 steps):

| | |
|---|---|
| span | 13516.7 ms |
| kernel busy | 1823.3 ms (13.5%) |
| copy busy | 2000.4 ms (14.8%) |
| union busy | 3823.6 ms (28.3%) |
| **GPU IDLE** | **9693.0 ms (71.7%)** |

The idle is not a few big stalls - it is ~354,000 operations each separated by a **~15 us median gap**:
261k gaps under 20 us (27.6% of span), 89k between 20-100 us (21.1%), and only 143 gaps over 1 ms
(16.7%, worth attributing separately). **95% of those operations are COPIES** (337k of 354k), so the
map traffic costs ~2.0 s of transfer PLUS ~5 s of serialization - closer to half the span than the 15%
a busy-time reading suggests.

A faithful reproducer (`amr-bench/attach/serial/`) matches the signature - 81.6% idle, 18.7 us median
gap against MFC's 71.7% / 15.4 us - and a clean sweep holding work and arrays fixed while varying only
how many regions the work is split into gives:

| regions fused | operations | span |
|---|---|---|
| 1 | 6555 | 202.4 ms |
| 2 | 3315 | 120.2 ms |
| 4 | 1695 | 73.8 ms |
| 8 | 885 | 56.0 ms |
| 24 | 345 | **37.6 ms** |

Net of ~29 ms fixed startup the variable time falls **20x for a 19x reduction in operations**, and the
median gap is ~18 us at every point. So:

> **span = busy + (operations x ~18 us).** The per-operation gap is irreducible; the COUNT moves, and
> so does the WORK each operation carries.

**The count is only half of it - WORK PER KERNEL decides which regime you are in.** Holding the operation
count fixed at 721 kernels / 5834 copies and varying only kernel duration:

| mean kernel | busy % | regime |
|---|---|---|
| 14 us | 20.5 | latency-bound |
| 16 us | 21.1 | latency-bound |
| 37 us | 25.6 | latency-bound |
| **201 us** | **51.2** | **crossover** |
| 763 us | 77.4 | work-bound |

Crossover is ~200 us of kernel work per ~8 copies. Scaled to MFC's ~20 copies per kernel that is
**~500 us**, and MFC's mean kernel is **109 us** - so MFC runs **4-5x below** the point where work
begins to dominate latency. Separately, raising work per region 256x at small sizes (16 -> 4096 cells)
changed the span by 5%: down in that regime the work is free and only the operation count matters.

**CHECKED AT PRODUCTION OCCUPANCY - the regime does NOT change with problem size.** Every measurement
above ran at ~15% of one GCD, so the obvious objection is that short kernels are an artifact of a toy
problem. Re-profiled at 4x the cells (8.4M -> 33.6M, ~60% of a 64 GB GCD):

| | 15% of GCD | 60% of GCD |
|---|---|---|
| span | 13516.7 ms | 28389.8 ms |
| busy | 28.3% | **28.9%** |
| GPU idle | 71.7% | **71.1%** |
| operations | 354,129 | 639,276 |
| median gap | 15.4 us | **15.5 us** |
| copies/kernel | 20.2 | 21.8 |

Busy fraction and per-operation gap are INVARIANT to problem size. The reason matters: under AMR the
work per kernel is set by `amr_max_grid_size`, **not** by the global grid - a larger domain yields MORE
BLOCKS OF THE SAME SIZE, so it moves along the operation-count axis while leaving work-per-kernel
untouched. Scaling up does not walk you out of the latency-bound regime; only raising the cap does.
(The 60% run exited 134 on the same `getTargetPointer returned null` mapping abort seen at 15% after
producing 28 s of trace; the ratios match the smaller run, but that abort is worth chasing on its own.)

So the governing relation is two-dimensional:

> **efficiency = kernel work / (kernel work + operations x ~21 us)**

This is why `amr_max_grid_size` is worth ~20x (70.3 s/step at cap 128 versus 8.5-9.5 at cap 1024) and
it is not a mysterious constant: a larger cap puts more cells under each kernel, moving the solver from
the latency-bound regime toward the work-bound one. It also upgrades kernel fusion, which improves BOTH
terms at once - a fused region carries more work and costs fewer operations.

**`nowait` is REFUTED as a lever.** Enqueuing the same dependent regions with `nowait` + `depend` and one
`taskwait` per advance measured 213 ms against 202 ms baseline, with the median gap unchanged (18.6 vs
18.4 us). Asynchronous enqueue does not collapse the serialization.

**Consequences, in priority order.**
1. **Fewer block advances** - `amr_max_grid_size`. Every advance multiplies the whole operation count.
   Already measured ~20x and already shipped with runtime advisories. Still the largest available win.
2. **Kernel fusion inside `s_compute_rhs`** (16 direct regions). The sweep above is the evidence that
   this pays, and it pays on BOTH factors: a fused region removes its own gap and its maps together.
   Note fusing raises copies-per-kernel slightly (8.1 -> 10.1 measured) because a fused region touches
   more arrays; total operations still fell 19x, so the trade is strongly favourable.
3. **The 143 stalls over 1 ms** (2253 ms, 16.7% of span, one of 543 ms) are host-side and unattributed.
   Cheap to investigate and potentially a large easy win.

What this retires: anything that reduces the COST of an operation rather than the NUMBER of them. See
the retraction below - seven such mechanisms were measured and all failed.

### RETRACTED 2026-08-03: the promotion plan below is DEAD, and so is its premise

**The premise was an attribution artifact.** "Dummy-argument kernels pay 33-54 maps per dispatch while
module-level `GPU_DECLARE`'d kernels pay 0" came from TEMPORAL attribution: correlation ids do not match
between the kernel and copy traces, so each copy was assigned to the next dispatch by timestamp. Copies
CLUSTER - 11.5% of inter-dispatch intervals contain zero copies and one contains 2893 - so a kernel
dispatched in a burst shows "0" whether or not it needs maps. The `m_amr_registers` kernels that
measured 0 are ~8% of dispatches, squarely inside that zero-interval share. The MEDIAN of 20 copies per
interval is real and matches the deep-member count, so the ~20 maps/dispatch AVERAGE stands; the
per-kernel split does not.

**Every mechanism was then measured directly by COPY COUNT** (`amr-bench/attach/clauses/`, one variant
per process under `rocprofv3`, no attribution needed, immune to the launch-time noise floor):

| mechanism | copies/kernel | verdict |
|---|---|---|
| deep-type dummies (MFC today) | 9.4 | baseline |
| flat-array dummies (no derived type) | 9.4 | no gain - retires "flatten the interface" for good |
| pointer dummies | 9.1 | no gain |
| `defaultmap(present:aggregate)` | 9.1 | no effect |
| `map(present,alloc:)` | 9.1 | no effect |
| `!$omp declare mapper` | 39.1 | **4x WORSE** |
| module state - static, allocatable, or the exact `freg` pattern | 32.8 | **3.5x WORSE** |

**Per-region map traffic is irreducible by any interface or clause change on this compiler, and moving
fields to module scope makes it worse.** That kills promotion, flattening, and the view/attach design
together. It also explains why `freg`-style code looked free: it was never measured, only attributed.

**What survives.** The tax itself is real and first-order (2.0 s of copy time against ~1.5 s of kernel
time). Since it cannot be reduced PER REGION, the only remaining levers are fewer regions and fewer
block advances:

1. **Fewer block advances** - `amr_max_grid_size`. Already measured at ~20x (70.3 s/step at cap 128
   versus 8.5-9.5 at cap 1024) and already shipped as a default plus runtime advisories. This remains
   by far the largest available win and needs no code.
2. **Fewer regions per advance** - kernel fusion inside `s_compute_rhs` (16 direct regions). Untested,
   touches numerics-adjacent code, and is the only untried lever left.

Do not re-attempt promotion, flattening, mapper, or attach without first re-measuring the table above;
all seven were tested on 2026-08-03 and all seven failed.

### PLAN (2026-08-03): promote the RHS working set from dummy arguments to module scope

**The insight that removes the blocker.** The attach hunt was solving the wrong problem. The tax is not
"pointers are slow" - it is that *fields arrive as dummy arguments*, so every target region re-maps their
deep members. Nothing needs to be aliased, attached, or flattened: the fields simply need to BE module
state. `m_rhs` already declares ~12 module-scope `GPU_DECLARE`'d derived-type arrays (`q_cons_qp`,
`q_prim_qp`, `flux_n`, `flux_src_n`, `qL_prim`/`qR_prim`, `tau_Re_vf`) and those kernels measure ZERO
maps per dispatch. Only `q_cons_vf`, `q_prim_vf`, `rhs_vf` arrive as dummies - exactly the ones taxed.

**The change, per routine: declarations and call sites only. Loop bodies are untouched.** Promote the
dummy to a module variable OF THE SAME NAME, so `q_cons_vf(i)%%sf(j,k,l)` still resolves - now to module
state. This is what makes the change tractable and what makes it safe for a bit-identical requirement:
the arithmetic is not edited at all. (Contrast the flat-store design, which threads a slot index through
every reference - a genuine 500-loop rewrite.)

**How AMR then feeds it.** The monolithic path pays NOTHING: the module arrays simply are its arrays.
The AMR path copies block *k*'s fields into them before the advance and `rhs` back after -
device-to-device, ~1.3 MB for a 128^2 2D block at `sys_size=5`, about **1.3 us**, against the measured
**~2.9 ms** of argument mapping per `s_compute_rhs` invocation (2.0 s / 696). Roughly a 2000:1 trade.
This is the `s_amr_swap_to_fine` idiom already used for grid state, extended from geometry to fields.

**Re-entrancy: VERIFIED SAFE, and it was the risk that could have killed the design.** Module-scope state
is only sound if no two block advances overlap. `s_amr_advance_children` is `recursive`, but the
recursion sits OUTSIDE the stage loop - every level-`clev` block completes all three stages, and only
then does the routine recurse into level `clev+1`. `s_l0_advance_stage_rhs` and the fine-advance loops
are likewise serial over blocks. No nesting, so no working-set stack is required.

#### MEASURED 2026-08-03, AFTER step 1a landed: where the operations actually are

Step 1a (`905e9d7c`) and the flat store itself (`cacc14ec`) are in. Before converting 78 call sites, the
per-module split was read out of four independent `rocprofv3` kernel-stats profiles
(`amr-bench/scratch/{l0p-hVcToX/p0,l0p-hVcToX/p2,att-dQ6b8U,sat2-UvUxRt}`). The two metrics disagree, and
which one governs is the whole question:

| | kernel TIME | LAUNCH count |
|---|---|---|
| `s_compute_rhs` tree (`m_variables_conversion`, `m_weno`, `m_riemann_solver_hllc`, `m_riemann_state`, `m_rhs`) | 88.7 - 91.1% | 60.6 - 66.7% |
| AMR-local (`m_amr`, `m_amr_registers`) | 5.1 - 10.1% | 32.7 - 38.2% |

**Launch count governs**, because the regime is latency-bound (cost proportional to operations, 71.7%
idle) - so the AMR-local kernels this plan set out to batch are worth about **one third of the operation
budget**, not the ~8% their kernel time suggests. That third is real and worth taking. The other two
thirds are inside the RHS tree and this plan does not reach them.

**THE BLOCKER, and it is an interface problem, not a mechanism problem.** Splitting the 78
`amr_slots(...)%%q_cons` references by what consumes them:

| kind | count | can the interface change? |
|---|---|---|
| element access `%%q_cons(i)%%sf(a,b,c)` | 15 | yes - free |
| whole-array pass to an AMR-local callee | 40 | yes - the callee is ours |
| whole-array pass to a SHARED solver routine | 8 | **no** |

The 8 are `s_compute_rhs` (x4), `s_ibm_correct_state` (x2), `s_pressure_relaxation_procedure`, and
`s_infinite_relaxation_k` - the last living in `src/common/`, shared by all three executables. All take
`type(scalar_field), dimension(sys_size)` and all are used by the MONOLITHIC path too, so their signatures
cannot follow `q_cons` into a flat store without dragging `q_cons_ts(i)%%vf` and the whole of `m_rhs` with
them. The pointer-view bridge that would avoid this (a `scalar_field` whose `%%sf` points into the store)
is the MWE's variant D1/D2, already shown to fail silently on this backend - the device sees a descriptor
holding a host address.

**Consequence for the plan below: sub-steps 4 and 5 as written are not reachable.** `q_cons` cannot
migrate while those 8 sites exist. And because `s_amr_lerp_fine_ghosts` and the 4579 `s_amr_fill_fine_ghosts`
call both terminate in a `q_cons` write, even the ghost-only migration stops short of a batchable path.
The remaining forks are recorded below; do not restart sub-step 4 without picking one.

#### DECISION + STATE (2026-08-03): fork 1 chosen - AMR-local migration behind a copy bridge

Of the three forks above the choice was **fork 1**: migrate the AMR-local sites to the flat store and
copy device-to-device around the 8 shared-solver call sites, leaving `s_compute_rhs` and friends with
their `scalar_field` interfaces untouched.

**Landed and golden-clean (76/76 each):**

| commit | what |
|---|---|
| `905e9d7c` | step 1a - dense local index |
| `cacc14ec` | the flat store + `s_amr_st_reserve` (device-preserving geometric growth) + `s_amr_st_finalize` |
| `bbf8ec8e` | **2a - the subcycle ghost sources migrated**; first change where data really flows through the store |

2a is +10 net lines because `s_amr_fill_fine_ghosts` is now generated for three targets (`_sf`, `_gsta`,
`_gstb`) from ONE source body using a Fypp accessor lambda -
``#:set QF = (lambda ix: ...) if TGT == '' else (lambda ix: ...)`` - the idiom `m_riemann_solver_hlld`
already uses for its per-direction stencil variants. Branching on the target inside a single region is
not an option: a dummy referenced in ANY branch is still mapped.

**2b - `q_cons` + the copy bridge. LANDED (see the results subsection below). It was ATOMIC and bigger than the original plan said.**

- **70 code sites** (not 59), across `m_amr.fpp` 47, `m_amr_regrid.fpp` 17, `m_amr_restart.fpp` 6,
  `m_time_steppers.fpp` 1.
- They cascade into **~17 AMR-local callee signatures**: `f_amr_rho_tot`, `s_l0_pack_unpack_block`,
  `s_amr_fine_slice`, `s_l0_edge_bc_tile`, `s_l0_copy_block`, `s_amr_restrict_pack_device`,
  `s_amr_gather_from_parent_field`, `s_amr_fine_rk_update`, `s_amr_copy_fine_fields`,
  `s_prolong_species_closure`, `s_prolong_one_var`, `s_prolong_alphas_closure`,
  `s_l0_fill_ghost_corners`, `s_amr_restrict_overwrite_device`, `s_amr_reflux_apply_faces`,
  `s_amr_fine_seam_exchange`, `s_amr_fill_fine_ghosts_sf`. Some take TWO block fields (parent+child,
  or two seam neighbours) and become two `loc` arguments - which is better for batching, not worse.
  Some take a MIX (`s_amr_copy_fine_fields(q_cons, q_cons_stor, ...)`); the unmigrated dummy is still
  mapped there, so that call site keeps its per-region tax until `q_cons_stor` follows.
- `m_amr_restart.fpp`'s 6 sites are HOST-side I/O with explicit `GPU_UPDATE` - retarget them at the
  store (whole-store update is wasteful but restart is rare).
- **It cannot be split.** The store is either authoritative for `q_cons` or it is not; a dual-write
  staging period is exactly the silent-divergence trap. Convert in dependency order, build often, run
  the goldens once at the end, and bisect by reverting individual routines if they fail.

**The bridge costs more than first assumed, so the payoff is smaller.** All four shared routines are
`intent(inout)`, and `s_compute_rhs` reaches `s_populate_variables_buffers`, which writes `q_cons`'s
buffer region - so the bridge must copy BOTH directions at all 8 sites. Each crossing is a region taking
a `scalar_field` dummy (~20 argument maps) while the batched kernels it enables cost 0. Per block
advance that turns ~9.3 AMR-local launches into ~2.4, i.e. **~25-30% of total operations, not the ~35%
launch share**. Quote the 25-30%.

**Harness note:** the session scratchpad rotates and ate a UUID list mid-run (the failure looks like a
kill, not a missing file). Keep test lists and logs under `amr-bench/logs/store/`.

#### 2b RESULT (2026-08-03): the store is authoritative for `q_cons`

Converted in dependency order, one build per group, goldens once at the end. Net **+86 lines** across
three files (`m_amr.fpp`, `m_amr_regrid.fpp`, `m_amr_restart.fpp`); `t_level%%q_cons` and its
`@:ALLOCATE`/`ACC_SETUP_SFs`/`@:DEALLOCATE` are gone.

**The site count was 80, not 70.** The earlier inventory grepped `%%q_cons` and filtered out lines
containing `q_cons_stor` - which also dropped every line carrying BOTH, hiding `s_amr_copy_fine_fields`
and `s_amr_fine_rk_update` (4 sites, 2 signatures). When inventorying a rename, filter on the token, not
on the line.

**Three callees turned out to be polymorphic** over "a block" and "the level-0 monolithic field", so they
became two Fypp-generated variants from one body rather than one converted routine:

| routine | `_st` (flat store) | `_sf` (scalar_field) |
|---|---|---|
| `s_l0_pack_unpack_block` | migration / scatter of a block | `q_cons_vf`, `coarse_tgt`, `rhs_delta` |
| `s_amr_restrict_overwrite_device` | fold a child into a parent BLOCK | fold into level-0 |
| `f_amr_rho_tot` | fine sensor on a block | coarse sensor on `q_cons_base` |

`s_amr_gather_from_parent_field` (and its two callees) needed the same split for a different reason: it
is called once with the parent's `q_cons` (store) and once with its `q_cons_stor` (still a
`scalar_field`), which is the subcycle's two-bracket gather.

**`s_amr_fill_fine_ghosts`'s `_sf` variant DISAPPEARED**, folding into the `_gsta`/`_gstb` pattern from
2a - all three targets are now `(q_coarse, loc)`. Same for `s_amr_lerp_fine_ghosts`, which lost its
`q_tgt` argument entirely. That is where most of the deleted lines came from.

**The bridge is 6 crossings, not 8**, because the two `s_compute_rhs` sites and the two
`s_ibm_correct_state` sites are if/else arms - the load/store hoists around the branch, so one round trip
per stage, not per arm. The sixth is `s_amr_reflux_apply_faces`, which is NOT one of the four shared
solver routines: it lives in `m_amr_registers`, which `m_amr` already `use`s, so reaching the store from
it would be a circular dependency. The bridge solved that at zero extra design cost.

**Bridge invariant: crossings must not nest.** `amr_cons_br` is a single shared buffer, so a load inside
a load would silently clobber the outer block's state. None of the six call chains re-enters AMR code;
check that before adding a seventh.

**Two traps hit:**

- ``#:for DIR, LHS, RHS in [('load', A, B), ('store', B, A)]`` fails `lint_source.py` - it reads the
  flattened list and reports `A` and `B` as duplicate entries. Make the direction a flag and derive the
  operands (``#:set LHS = BR if DIR == 'load' else ST``), which is clearer anyway.
- `./mfc.sh build ... | grep error:` reports GREP's exit code, not the build's. Redirect to a log and
  echo `$?` on its own line. (Same class as the `--only` comma trap already recorded.)

**Confirmation that 2b was the right prerequisite:** `s_amr_fine_seam_exchange`'s own comment named this
change as its blocker - *"Batching further, across PAIRS, is NOT possible without a flat per-slot backing
array: each slot's `q_cons(i)%%sf` is an independent allocation and `amr_slots` is not `GPU_DECLARE`'d, so
a runtime slot index inside a kernel is a null deref."* Both blocks are now addressed by a plain integer
subscript into a `GPU_DECLARE`'d module array. Same for `s_l0_pack_unpack_block` and
`s_l0_fill_ghost_corners`, whose dummies existed only to dodge that null deref.

**2b is a toll, not a win.** It adds two whole-block device copies per crossing and, per the attach/map
measurements, saves none of the ~20 argument maps at the crossing itself. All of its value is unlocking
step 2. Measure step 2's operation count with `amr-bench/attach/serial/gaps.py` before believing any of
the projected 25-30%.

#### EXECUTION PLAN for 1b + 2 (written 2026-08-03, SUPERSEDED IN PART - read the section above first)

**1b - migrate `{q_cons, q_ghost_a, q_ghost_b}` to the flat store.** They move together because
`s_amr_fill_fine_ghosts` targets all three; see the scoping constraint above.

```fortran
real(stp), allocatable :: amr_cons_st(:,:,:,:,:)   ! (x, y, z, var, LOCAL slot)
real(stp), allocatable :: amr_gst_a(:,:,:,:,:), amr_gst_b(:,:,:,:,:)
$:GPU_DECLARE(create='[amr_cons_st, amr_gst_a, amr_gst_b]')
integer :: amr_st_cap = 0
```

Order of work, each sub-step compiling and golden-clean before the next:

1. Declare the stores + `s_amr_st_reserve()` (geometric growth on `amr_loc_n`). Call it from
   `s_amr_alloc_slot` after `amr_slot_live = .true.`. Guard the ghost stores on `amr_subcycle`.
   `q_cons` growth MUST preserve device contents (unlike the ghosts) - copy old -> new before freeing.
2. Convert `s_amr_fill_fine_ghosts`: replace the `q_fine` dummy with `(target, loc)` where target
   selects `amr_cons_st`/`amr_gst_a`/`amr_gst_b`. **Duplicate the target region per branch** rather than
   branching inside one region, or the unmigrated dummy is still mapped and the win is lost.
3. Convert `s_amr_lerp_fine_ghosts` (1 write site) and the 5 ghost call sites.
4. Convert the `q_cons` readers/writers. 59 references across 4 files - do it file by file, building
   between each. `amr_slots(k)%%q_cons(i)%%sf(a,b,c)` becomes `amr_cons_st(a,b,c,i,amr_loc_of(k))`.
5. Delete `q_cons`/`q_ghost_a`/`q_ghost_b` from `t_level` and from alloc/free.

**2 - batch the per-block kernels.** Replace `do k = ...; call s_amr_select_slot(k); <one kernel>` with
a single region carrying `loc` as an extra collapsed dimension over `1:amr_loc_n`. Start with
`s_amr_fine_rk_update` (self-contained, 2 regions), then the RHS path. Requires that per-block loop
bounds (`amr_slots(k)%%m/n/p`) also be readable on device - stage them into a small
`GPU_DECLARE`'d `integer :: amr_blk_m(:), amr_blk_n(:), amr_blk_p(:)` indexed by `loc`.

**Verification protocol, non-negotiable.** Every sub-step is a pure data-layout or scheduling change,
so **all 76 AMR/L0 goldens must stay byte-identical**. A diff means the change is wrong, not that the
goldens need regenerating. In addition, after step 2 confirm the mechanism actually engaged by
re-running the gap analysis (`amr-bench/attach/serial/gaps.py` on a `rocprofv3` dir) and checking that
the operation count fell - a passing golden set with an unchanged operation count means the batching
did not take effect.

**Traps already paid for, do not rediscover:**
- `amr_slots` is allocated in TWO places (`s_initialize_amr_module`, behind `if (.not. amr) return`, and
  `s_l0_tiles_init` for pure-L0). Anything allocated beside it must follow the POOL. Fingerprint: only
  the `L0 tiles` goldens fail.
- Print the text being REPLACED on every non-trivial substitution. A `max(..., bub_pos_frac*u0)` clamp
  was nearly rewritten to `max(..., 0._wp)` - a silent physics change that only a Lagrange-bubbles
  golden would have caught, hours later.
- Do NOT pass the store as a dummy argument. A plain-array dummy costs the same per-region map traffic
  as a deep-type one (9.4 vs 9.4); only a plain MODULE array reaches 0.
- The tool timeout is 10 minutes and an MFC rebuild exceeds it - background the build and poll a log
  marker, and never write `pgrep -f "<string in your own command line>"`, which matches the poller.

#### De-risking MWEs, each run BEFORE the matching phase

Build them as one matrix in a single compile (`amr-bench/attach/mtu/`), not one hypothesis per MFC
rebuild - a 6-phase matrix compiled in ~20 s eliminated four hypotheses at once and correctly predicted
an MFC failure, where the same questions cost ~20 minutes each in-tree.

| MWE | de-risks | pass condition |
|---|---|---|
| M1 promotion, multi-TU | module state referenced from a DIFFERENT translation unit than it is declared in | correct result; per-region time at the ~22 us module level, not the ~110 us dummy level |
| M2 copy-in/advance/copy-out | the AMR feed path, with host and device data DELIBERATELY DIVERGED | correct result; copy cost measured and compared against the mapping it replaces |
| M3 serial reuse, varying extents | one shared working set reused by many blocks of DIFFERENT sizes | correct for every block, in any order |
| M4 sub-block extents | working arrays sized at `mbuf` max while a block uses a sub-range | no read/write outside the block's own range |

M2's divergence requirement is not optional: uniform host/device data hid a data-movement defect through
seven passing variants in the session that produced this note.

#### MWE RESULTS (2026-08-03, `amr-bench/attach/promote/`): correctness de-risked, 5/5

| phase | de-risks | result |
|---|---|---|
| S1 | promotion through a 3-deep call chain, cross-TU, over a sub-range | PASS |
| S2 | POISONED halo beyond the block untouched - no out-of-range write | PASS |
| S3 | 25 advances with host/device DELIBERATELY diverged, plus unrelated target regions between them | PASS |
| S4 | 48 blocks of VARYING extents cycled through one shared working set | PASS |
| S5 | working set REALLOCATED mid-run | PASS |

**Mechanism detail that matters for implementation.** A module-scope deep type needs the parent
descriptor mapped AND its members mapped; `map(to:)` on the parent alone leaves `%%sf` null on device
and the first kernel faults at a nil address. With promotion this setup runs **once at init** and the
working set is never re-pointed - which is exactly why promotion sidesteps the attach problem that
killed the view-based design.

**Timing could NOT be measured cleanly and no speedup is claimed here.** The node showed a ~134 us
per-region launch floor, uniform across all 8 GCDs (129-135 us) with the GPUs otherwise idle, against
~20 us for the same construct earlier the same day. Within-program A/B gave 448 vs 501 us/region: 53 us
saved for 2 dummy arrays, ~26 us per array per region. Extrapolating instead from MFC's own profile
(485 copies / 24 regions = 20 maps per region, 2.0 s / 696 invocations = ~6 us per map = ~120 us per
region, against a 69 us mean kernel) implies roughly a **2-3x per-region reduction**. Treat the 4.4x
from the earlier quiet-node reproducer as an UPPER BOUND and re-measure on an idle node before quoting.

#### Phases, each golden-gated and independently revertible

1. **`s_amr_fine_rk_update`** - AMR-only, one file, three arrays. Smallest in-tree proof of the pattern.
   Gate: all AMR goldens byte-identical AND its copies/dispatch measured 54 -> ~0.
2. **`s_compute_rhs`'s three field dummies** - 5 call sites. Benefits the uniform solver too (it pays
   18.6 copies/launch to AMR's 19.2). Gate: FULL suite byte-identical, plus a uniform-case timing.
3. **Walk the call chain** - `m_riemann_solver_hllc` (44 maps/dispatch), `m_riemann_state` (33),
   `m_viscous`, `m_weno`; ~57 call sites pass these three arrays. Gate: full suite, per-stage.

Stop after any phase whose measured gain does not survive its own run-to-run spread, and say so.

#### What would falsify this plan

- M1 shows promotion does not remove the maps once the reference crosses a TU boundary.
- The copy-in/out is not device-to-device (a host round trip would cost far more than the tax).
- Phase 1 changes any golden. This is a storage-location change with identical arithmetic; a diff means
  the model is wrong, not that the goldens need regenerating.

### SETTLED 2026-08-03: concurrency (option C) is DEAD; fusion (option B) is MEASURED

Both settled on one harness, `amr-bench/attach/serial/conc.f90` (wall clock via `omp_get_wtime`;
`NBLK` blocks x `NREG=15` dependent regions x `NADV=30`, deep-allocatable dummies = the MFC interface).

**Option C - concurrent block advances - produces ZERO overlap and is slower.** Each block was given
its OWN `a/b/c/d` arrays AND its own `depend` token, so the `NBLK` chains were genuinely independent:
the best possible case, with nothing to serialize on.

| NBLK | serial | conc | |
|---|---|---|---|
| 1 | 0.059 s | 0.072 s | 1.23x SLOWER |
| 4 | 0.235 | 0.262 | 1.11x slower |
| 16 | 0.948 | 1.031 | 1.09x slower |
| 64 | 3.91 | 4.13 | 1.06x slower |

Both scale perfectly linearly in `NBLK` - no hint of overlap - and a `rocprofv3` trace of the `conc`
variant confirms it directly: **7201 kernels, 0.0% overlapped**. `nowait` plus independent per-block
`depend` tokens produce no device-side concurrency at all on this backend (amdflang/AFAR, gfx90a).

**This is NOT a repeat of the earlier `nowait` result.** `s.f90`'s `k_nowait` put `depend(inout: d)` on
ONE shared token, serializing every region into a single chain; it asked whether a SERIAL chain
pipelines. The question of whether SEPARATE chains overlap was open until now. It is now closed.

Corollary: the ~17 module-scope `GPU_DECLARE`'d families in `m_rhs` that would have to be replicated
per stream are MOOT - concurrency fails before sharing ever becomes the constraint. **Do not re-attempt
concurrency, streams, or async task graphs.** With this, TEN mechanisms for reducing or hiding
per-region cost have been measured and all ten failed.

**Option B - fewer regions - is the whole game, and wall time is LINEAR in region count.** Same
arithmetic in every variant (`result=21.00`), only the region count differs:

| regions/advance | wall | vs 15 regions | s per region |
|---|---|---|---|
| 15 (= serial) | 0.945-0.955 s | 1.00x | 0.063 |
| 5 | 0.313-0.429 | 3.0x | 0.074 |
| 3 | 0.188-0.190 | 5.0x | 0.063 |
| 1 | 0.064 | 14.8x | 0.064 |

`wall ~ 0.064 s x regions`, flat to +-15% across a 15x range (the 5-region point is the one noisy
sample: 0.313 and 0.429 on two runs). `FUSE=1` reproduces `serial` (0.955 vs 0.945), which is the
harness checking itself.

**Re-pricing option B against the real profile.** The harness is 97.5% idle with ~6 us kernels; the real
AMR run is 72% idle with 91 us average kernels, so fusion helps proportionally less there. Using the
real split - the `s_compute_rhs` tree is 55.5% of launches at ~14.7 regions per block-stage - a 3x
region reduction inside it removes ~37% of all operations, i.e. **~26% of wall**. That is roughly 2.5x
the ENTIRE remaining step-2 batching programme (~10% across three 2b-sized field migrations), and it
also benefits uniform runs, which AMR-local batching does not.

**What this means for the plan.** Every surviving lever is the same lever: FEWER REGIONS. Fusion (B),
batching (A), and bigger blocks (`amr_max_grid_size`, shipped, ~20x) are three applications of it.
Nothing that reduces per-region COST works, and nothing that OVERLAPS regions works. Rank by measured
value: B ~26% > A-remainder ~10% > C = 0. The concrete next target is `m_weno` + `m_riemann_state` -
3216 launches each, ~4.2 regions per block-stage apiece, adjacent in the pipeline with a direct
producer/consumer dependency.

Method notes worth keeping: `cpu_time` is WRONG for this measurement (it sums CPU across threads and
made `conc` look merely 2-16% slow when the wall-clock gap was different); use `omp_get_wtime`.
`rocprofv3` needs `--output-format csv` or it writes only a `.db`. And the pre-existing `k_fused` in
`s.f90` is NOT arithmetically equivalent to its `k_base` (`result=12` vs `30`), so its timings cannot
be used to price fusion - `conc.f90`'s variant was written to be equivalent and checks it.

### MEASURED 2026-08-03 AT PRODUCTION SIZE (3D, 400^3): the cap is exhausted, balance is state of
### the art, and AMR captures ~11% of its own potential

Everything below is on `amr-bench/cases/sc3dx_amr.py` - 400^3 = 64M base cells, np=8, ~48% of a 64 GB
GCD per rank at the ~1M points / 2 GB rule. **All prior AMR pricing in this document came from a 2D case
at ~4% of GCD capacity and should not be trusted where it disagrees with this section.**

**`amr_max_grid_size` is EXHAUSTED in 3D.** Clean U-curve with the optimum exactly at the shipped
default, walled off above by device OOM rather than by tuning:

| cap | boxes | boxes/rank | s/step | ns/cell |
|---|---|---|---|---|
| 12 | 1089 | 136 | 22.18 | 258.2 |
| 16 | 2401 | 300 | 27.34 | 299.5 |
| 24 | 1089 | 136 | 16.76 | 171.3 |
| **32** | **625** | **78** | **9.5-9.8** | **93.9-96.1** |
| 48 | 324 | 40 | 16.73 | 158.5 |
| 64, 96, 128 | - | - | **OOM** | - |

cap 32 reproduced across two independent sweeps to 2%. cap 16 is anomalous (2401 boxes, MORE than
cap 12's 1089) - the clusterer tiles pathologically there; an oddity, not a trend.

**This OVERTURNS the recorded 2D conclusion** that bigger blocks are "by far the largest available win"
(~20x, cap 128 -> 1024). That is a 2D result. In 3D the cap costs memory as cap^3 (the validator says so
explicitly: solver scratch is sized to the cap, "growing as the cap raised to the dimension count"), so
the base problem and the cap compete for the same memory and the ceiling arrives immediately.

**CORRECTION: OOM does NOT present as a hang.** It aborts at exit 134 with
`HSA_STATUS_ERROR_OUT_OF_RESOURCES` and an explicit message. The wrapper's `rc=143` is the OUTER
timeout, two layers up - reading it as a timeout was wrong.

**AMR vs uniform, both converged.** Uniform 400^3 needed ~100 steps to converge (Time Avg 0.493 at 21
steps -> **0.7068** converged, a 43% error); AMR needed ~60 (9.54 single-sample -> **~10.24** Time Avg,
still drifting down). **20-step runs are NOT converged - the cap table above therefore holds as a
RELATIVE comparison (identical protocol per arm) but its absolute ns/cell understates AMR by ~7%.**

| | cell-updates/step | s/step | ns/cell |
|---|---|---|---|
| uniform 400^3 | 64.0M | 0.7068 | **11.05** |
| AMR cap 32 | 101.5M (64.0M base + 37.5M fine) | ~10.24 | **~100.9** |

- **AMR per-cell penalty = ~9.1x** - NOT the ~31x on record, which came from the small 2D case.
- Matching AMR's finest resolution uniformly is 1600^3 = 4.096e9 cells: a **40.4x** cell-update saving.
- AMR converts that into **~4.4x** of wall clock, i.e. it captures **~11% of its own geometric
  potential**. The missing ~9x IS the per-region cost.
- AMR is unambiguously WINNING here (an earlier worry that it might be a net loss is retracted), and
  1600^3 uniform would need ~8.2 TB against 512 GB available - it does not merely lose, it does not fit.

**LOAD BALANCING IS NOT THE PROBLEM - it is state of the art.** Measured on this case:
`max/mean 1.016-1.018, ranks_with_no_fine_block 0` -> **efficiency 0.982-0.984** on AMReX's own metric
(mean load / max load). Published AMReX values: knapsack 0.97-1.00, Painter's SFC 0.92-1.00, original
SFC 0.80-0.95 (arXiv 2505.15122). **We match knapsack and beat their original SFC.** That paper also
explicitly declines to claim production wall-clock savings for ANY algorithm - so our own 0.4% recovery
from rebalancing is not a sign of being behind, it is what converting a good epsilon into seconds looks
like. Do not spend effort here: the remaining 2% of balance is noise against 72% idle.

**THE REAL GAP vs state of the art is BLOCK COUNT.** AMReX operates at 4-16 boxes per rank ("at least 4
... more than 16 starts to cause performance issues"). Our OPTIMUM is 78 boxes/rank, and every cheaper
cap is worse (300/rank at cap 16). We are ~5x outside their regime at our best setting and cannot move
toward it, because per-block memory is cap^3 x ~17 module scratch families where their MultiFab is ONE
contiguous array per level. Same root cause as the 19.5 copies/launch: per-block replication instead of
a packed representation.

**What this means for priorities.** Fusion (~26%) and the step-2 batching remainder (~10%) are margin
work against a ~9x gap between AMR and its own potential. Both of those numbers are also 2D-derived and
must be re-priced here before they are quoted. The axis that matters is per-block overhead and block
count - which is what the flat store began.

**Harness rules this cost four self-inflicted failures to learn** (three already documented and hit
anyway): never run a benchmark while another `mfc.sh run` is live (they contend for the same GPUs and
the control silently never starts); never `pgrep -f` a pattern that appears in the poller's own command
line (an `until ! pgrep` loop then never exits - three stalled monitors and counting); pass `mfcrun.sh`
an ABSOLUTE case path (it `cd`s to the tree first, so a relative path resolves there); and never edit a
shell script while an instance of it is executing (bash reads scripts incrementally, so the running copy
then fails to parse, at a line number that no longer means anything). `mfcrun.sh` now enforces the first three.

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
