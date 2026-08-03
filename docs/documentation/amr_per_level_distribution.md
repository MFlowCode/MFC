@page amr_per_level_distribution AMR per-level distribution

# AMR per-level distribution

Design for the next phase of block-structured AMR: distribute each refinement level
independently, so AMR strong-scales instead of degrading as ranks are added. Written
2026-07-29 against `up/mega` @ `7b1e4933`. Companion to @ref amr_block_batching, which covers
the per-block launch cost, and @ref amr_multilevel, which covers the nesting.

**How to read this.** "What actually blocks exascale" is the primary open work item and comes first.
"Status" says what has landed. "The cost model is regime-dependent" and "Superseded measurements"
together record what measurement has and has not established - several earlier conclusions in this
file were retracted, and the retractions are kept rather than deleted. "Sequencing" is the queue.

## Status

All four design steps are landed. The measured outcome did **not** confirm the performance thesis,
and the reason matters more than the steps did — see "The cost model is regime-dependent" and
"Superseded measurements".

| Step | Commit | Result |
|---|---|---|
| 1. Scratch decoupling | `86782249` | landed |
| 2. Rank-independent cap | `a108dd37` | np=4->8 flipped 1.16x slower to 0.84x faster (single level) — **but see the caveat below** |
| 3. P2P parent<->child | `6832d299`, `d53fac46` | landed; 3d `26e0d080` hoists the level advance to lockstep |
| 4. Per-level mapping | `cfdd2847` | landed; correctness proven with a split tower, **no measured speedup**; the A/B that produced that result was run too small to interpret — see "Superseded measurements" |
| 5. Balance metric | `fc53e097`, `3780f30a` | landed; `[amr-balance]` per-level max/mean behind `load_weight_wrt`. Its *conclusion* ("box supply binds") is **retracted** — see below |
| 6. Model vs measured | — | done; **the cost model does not predict measured time**. Redirects the effort to step 7 |

Landed alongside, from auditing the above (2026-07-31):

| Change | Commit | Note |
|---|---|---|
| Multi-level subcycle at np>1 un-gated | `8dca4b65` | its stated blocker had been removed by `3db24df0`; golden `C45DBB52` |
| Level-general child slot cap | `be94db38`, `1e07eb65` | fixed `/2` was level-2-specific; **verified to fail without** — old cap core-dumps at `amr_max_level=3` |
| Brand-new-region box skipped the cap | `cfbacebe` | **heap corruption** under the pinned cap; golden `00EB793A` |

### Caveat on step 2, and the coverage hole behind it

`amr_max_grid_size` had **zero test coverage** until `00EB793A` — not one golden or example case set
it, only the schema, validator and parameter definitions. It is the mechanism this whole step rests
on (pin the cap → many small boxes → `boxes_per_level >> num_procs`), and nothing exercised it.

That is not hypothetical. `s_amr_regrid`'s brand-new-region branch emitted level≥2 boxes without ever
consulting `amr_maxc_fit` — the one box emitter that skipped `s_amr_tile_box` — so a child over the
cap made `s_amr_build_block_coords` write past `x_cb`, corrupting the heap on every regrid and
surfacing much later as `corrupted size vs. prev_size` inside an unrelated `free()`. Under `-O2` it is
silent. It needs all four of: pinned cap ≡ 0 (mod 8), `amr_max_level ≥ 2`, `num_procs ≥ 2`, and a
scattered IC (a Sod-like patch grows regions that already have fine data, so the branch never fires).

**Therefore step 2's measured numbers were collected on a code path that could silently corrupt the
heap, and a corrupting run still produces plausible timings.** The conclusion is probably still right
— the defect is a bounds overrun, not a change to the distribution — but it is no longer a clean
measurement, and it should be re-run on fixed code before being cited as settled.

The general lesson is the one worth keeping: **a parameter with no golden is a parameter whose every
failure mode is invisible**, and the plan leaned on this one for months.
| `no_blocks_ranks` is not idleness | `3ab65ab6` | see "What binds at scale" |

Two findings from the step-4 A/B that redirect this work:

1. **Per-level distribution shows no benefit at np<=8** on a 3-heavy-tower case (flat within
   noise at every rank count). Expected from the granularity floor below, but **not attributable**
   — there is no imbalance metric, so a flat result cannot distinguish "nothing to fix" from
   "balancer did nothing".
2. **Per-box overhead, not balance, is the dominant cost.** At np=1 the same case runs ~40x the
   coarse solve while the refinement it adds accounts for only ~2.5x. That ~16x residual belongs
   to @ref amr_block_batching, not here.
3. **The sweep was run too small above np=2.** The uniform control reaches only 18% parallel
   efficiency at np=8, so both arms are latency-bound exactly where distribution should have
   mattered — finding 1 is therefore inconclusive rather than negative. Check the control's
   efficiency before trusting any A/B run against it.

The consequence for sequencing: **instrumentation precedes further distribution work.** Steps 1-4
were gated on correctness and end-to-end s/step, neither of which measures balance.

## Where this stands (2026-08-02), and what large scale actually demands

**The goal is state-of-the-art AMR with load balancing at large scale.** This section is the current
measured position against that goal; everything below it is the history that produced it, including
several retracted conclusions.

### Measured position

| question | answer | evidence |
|---|---|---|
| Does AMR strong-scale? | **yes, 0.790 at np=8** (was 0.426) | 3D 256^3, rank-invariant box set, n=3 |
| Does AMR weak-scale? | **yes, 0.846 at np=8** | constant 2.10M base cells/rank, work/rank held to 1.00-1.10x |
| Is load balance the constraint? | **no** | imbalance 1.000-1.071 and ZERO idle ranks at every cap, even at 152 blocks/rank |
| Is AMR efficient? | **no - 4.2x uniform per cell** at the best cap | cap sweep, np=8 |
| Biggest single lever | `amr_max_grid_size`: 64 is **3.2x** faster than 32 | cap sweep |
| Beyond one node? | **UNKNOWN** | nothing here exceeds 8 GCDs on one node |

Two cautions on those numbers. The scaling figures are AMR compared to ITSELF at constant work; the
uniform control at these sizes runs at ~4% of device capacity, so comparisons against it are
meaningless (an earlier version of this document claimed "AMR scales better than uniform" - that was
an artifact of a starved baseline and is retracted). And the strong/weak efficiencies come from a
single case at a single cap on a single node.

### What produced the gain, and what that implies about this plan's method

The 0.426 -> 0.790 improvement did NOT come from per-level distribution. It came from removing a
**loop-invariant call**: `s_amr_exchange_coarse_cons_halo` sat inside per-block loops in both the
lock-step and subcycle paths, exchanging a coarse field that every fill reads and none writes -
~960 identical whole-subdomain exchanges per step where 3 suffice. Because the count tracked the
GLOBAL box count while each exchange costs the same, it did not distribute at all.

It was found by top-down phase attribution in one run. It was NOT predicted by any cost model in this
document, and the hypotheses this document did advance - replicated metadata as "THE exascale
blocker", per-box collectives, cross-rank adjacency, kernel launch count - were each measured and
found to be minor or wrong (regrid is 2.4% of the step and SCALES DOWN; launches are 8.5% of GPU
time). **Attribute before hypothesizing.** `amr-bench/audit.sh` encodes the measurement failure modes
that produced the wrong conclusions.

### The tension that governs large scale

Load-balance freedom requires `boxes_per_level >> num_procs`. Per-block cost is approximately FIXED
per block. These pull in opposite directions, and rank count decides which wins:

- At np=8, cap 64 gives 10 blocks/rank - enough for perfect balance (1.000, no idle ranks) and few
  enough that per-block cost is tolerable (4.2x uniform per cell).
- At 10^3-10^4 ranks, holding 10 blocks/rank means 10^4-10^5 blocks GLOBALLY. That simultaneously
  (a) multiplies the fixed per-block cost by the block count, and (b) reactivates every O(global
  boxes) term, including limit 3's replicated metadata (~560 MB/rank at 10^7 boxes).

**So per-block cost is escapable at 8 ranks and unavoidable at 10^4.** Raising `amr_max_grid_size`
is the lever today, but it is bounded by device memory (cap 96 device-OOMs in 3D) and by the block
supply the balancer needs. At large scale the code is forced into the many-blocks-per-rank regime,
which is precisely the regime where per-block cost dominates. **That, not balance, is what stands
between this implementation and large-scale AMR.**

### Configuration guidance: the parameters are worth more than the code

Measured 2026-08-02, 3D, np=8, 256^3 two-slab interface, `fine_work` IDENTICAL in every arm
(17765376) so all comparisons are like-for-like refinement.

| lever | setting | gain | what it costs |
|---|---|---|---|
| `amr_max_grid_size` | 32 -> **64** | **3.2x** | nothing measurable here; bounded by device memory (96 OOMs in 3D) |
| `amr_regrid_int` | 2 -> **8** | **1.39x** | the box set lags a moving feature |
| `amr_subcycle` | F -> **T** | **1.55x** | a DIFFERENT time integration, not a free optimization |
| **combined (MEASURED, not multiplied)** | | **7.012x** | no code changes |

Measured cumulatively at np=8, each row adding one lever, `fine_work` identical across the last three
arms (17765376) so the comparison is like-for-like refinement:

| arm | s per unit PHYSICAL time | vs default | marginal |
|---|---|---|---|
| default (cap 32, lock-step, regrid 2) | 6.716e4 | 1.000x | - |
| + `amr_max_grid_size` 64 | 2.210e4 | 3.039x | 3.04x |
| + `amr_subcycle` | 1.408e4 | 4.771x | 1.57x |
| + `amr_regrid_int` 8 | **9.578e3** | **7.012x** | 1.47x |

**The three levers compose cleanly** - every marginal contribution reproduces its standalone value
(3.2 / 1.55 / 1.39), because they act on independent things: the cap cuts per-block invocation count,
subcycling cuts how often the coarse level is integrated, and the regrid interval cuts tag-sweep
frequency.

**A retraction.** An earlier draft of this section warned that naive multiplication (~6.9x) was
"impossible, because it would put AMR faster per cell than a uniform grid". That objection was wrong:
it compared against a uniform per-cell rate computed at FIXED dt, but subcycling CHANGES dt, so the
comparison rather than the composition was invalid. The measured total is 7.0x.

For comparison, the loop-invariant halo hoist landed the same day is worth 2.4x, and the entire
remaining code agenda - flat backing store (~1.13x), P2P batching (~1.08x), per-phase rebalancing
(<=1.24x) - totals under 1.6x even if all three land. **The parameters are worth more than every code
lever combined, several times over.**

### Per-phase load imbalance: the balancer optimises the wrong quantity

Whole-step balance reads 1.000-1.010 at every rank count tested, including 32 ranks on 4 nodes. But
per-phase imbalance (`[amr-imbal]`, same `rank_time_wrt` gate) is 1.12-1.64, and it worsens off-node:

| phase | np=8 max/mean | np=16 max/mean | wait s/step at np=16 |
|---|---|---|---|
| fill | 1.229 | 1.475 | 0.128 |
| coarse_rhs | - | 1.447 | 0.085 |
| reflux | 1.212 | 1.348 | 0.080 |
| halo | 1.357 | 1.611 | 0.065 |
| restrict | 1.453 | 1.643 | 0.059 |
| coarse_halo | 1.133 | 1.476 | 0.017 |

**19.4% of the np=16 step is ranks waiting at phase syncs.** The cause is structural: the balancer
distributes TOTAL weight, but every phase has its own barrier and each is driven by a DIFFERENT
quantity - `reflux` touches only level-1 blocks, `restrict` is per-block, `coarse_rhs` is per-cell. A
rank can carry exactly the right total weight and still hold the wrong number of level-1 blocks, so it
arrives late at that phase's barrier while every other rank waits. Perfect per-phase balance is
unattainable (blocks are indivisible), so 1.24x is a ceiling, not a target - but it is the largest
measured code-side lever, and it is confined to the balancer rather than the solver.



**`amr_max_grid_size`.** The default of 0 derives the cap from the decomposition, so it SHRINKS as
ranks are added - backwards for strong scaling. AMReX defaults to 32 in 3D but its GPU guidance pushes
to 128; measured here, 64 is 3.2x faster than 32 and 10.8x faster than 16, with imbalance 1.000 and
ZERO idle ranks at the fastest point. Granularity, not balance, is the binding constraint in 3D.

**`amr_regrid_int`.** `regrid` is ~9% of the step and, importantly, FLAT in global box count (1.8x
across a 15x box-count range), so it is dominated by the per-cell tag sweep rather than the box
machinery. Regridding less often recovers more than regrid's own share (1.39x for interval 2 -> 8),
which means it triggers downstream work - slot rebuild, migration - that the phase timers do not
attribute to it. Risk: a fast-moving feature can outrun a stale box set. This case's refinement is
stable so it costs nothing; a strong moving shock would under-refine.

**`amr_subcycle`.** THE METRIC MATTERS. Lock-step advances every level at the FINEST-stable dt, so the
coarse grid is integrated 2**amr_max_level = 4x more often than its own stability requires. Berger-
Oliger subcycling advances L0 at its own dt, which is 4x larger: per unit physical time the work drops
from 34.5M to 22.0M cell-updates, predicting 1.57x - measured 1.551x, agreement to 1%. On a per-STEP
basis the subcycled arm looks 2.6x WORSE (2.8564 vs 1.1079 s/step) because each of its steps covers 4x
the physical time; comparing s/step would reject the better algorithm. **This is not a drop-in:** the
two are different time integrations and the answers differ beyond roundoff.

### Gap against the state of the art

| | AMReX | Parthenon | MFC today |
|---|---|---|---|
| per-block compute | box bounds passed as ARGUMENTS to a lambda over a POD view | `MeshBlockPack` - many blocks in ONE kernel launch | swaps GLOBAL state (`m/n/p`, `idwint/idwbuff`, coords) and calls the monolithic solver PER BLOCK |
| ghost fill | `FillPatch` over a whole MultiFab: one aggregated communication phase for the level | device-resident, packed | per BLOCK, with every rank participating in every block's gather |
| data residency | device-resident | "all data in device memory" | per-slot allocations; `~626 descriptor copies per RHS call` |
| demonstrated scale | production exascale | **92% weak scaling to 73,728 GPUs** | 8 GCDs, one node |

The architectural difference is *granularity of both compute and communication*: both reference codes
operate on a whole level at once, MFC operates per block. Parthenon adopted packing for exactly the
symptom measured here - kernel runtime smaller than the per-invocation overhead. That is the target
to match, and it is the same flat-backing-store restructuring @ref amr_block_batching describes -
but justified by the LARGE-SCALE argument above rather than by the packing rationale, which was
disproved (see "the packed super-grid cannot work").

### Open work to finish the AMR + load-balancing plan

Scored against the three tracks' own "done when" criteria, 2026-08-02.

**Track 2 (multi-level at np>1): COMPLETE.** Towers no longer co-locate (see the per-level
distribution comment in `s_amr_assign_block_owners`), P2P parent<->child landed, and the L2 seam
abort is gone - `s_amr_check_seam_topology` now rejects only genuinely illegal geometry.

**Track 3 (batch per-rank block advances): CLOSE AS SUPERSEDED.** Its criterion was a strong-scaling
curve showing blocks no longer serialize through one slot. Strong scaling reached **0.790 at np=8
(from 0.426) with the single-slot model untouched** - the limiter was a loop-invariant coarse-halo
exchange. Packing is disproved; the flat backing store measures ~1.13x at the cap the code should run
at. The mechanism was never built and no longer needs to be.

**OPEN 1 - ANSWERED, negatively: the coarse-grid balancer does not pay for itself.** Every earlier
measurement in this document ran `l0_ntile = 0`, so **every balance figure here (1.000-1.010, zero idle
ranks, to 32 ranks on 4 nodes) is a FINE-LEVEL result only**; the coarse level was a fixed Cartesian
decomposition. That path has now been benchmarked on a corner-concentrated 2D case (4096x2048, blob IC
so refinement sits entirely in one rank's half, np=2, `l0_ntile = 2` -> 8 tiles / 2 ranks, 40 steps,
3 reps):

| arm | `l0_ntile` | `l0_rebalance_interval` | med s/step | spread | vs A |
|---|---|---|---|---|---|
| A monolithic | 0 | - | 0.8892 | 0.9% | 1.000x |
| B tiled, no rebalance | 2 | 0 | 1.2045 | 0.8% | **0.738x** |
| D tiled + rebalance | 2 | 4 | 1.1991 | **11.3%** | 0.742x |

**Tiling costs 35%; rebalancing recovers 0.4% of it - inside D's own spread.** The rebalancer was live
(9 invocations, confirmed in-log, not a wiring failure). It migrated exactly once in 40 steps, and that
migration made the load gap **6.8x worse** and never recovered:

```
t_step=32  load-gap 5.738E-02 -> 3.898E-01  (1 migrations)
t_step=36  load-gap 3.488E-01 -> 3.488E-01  (0 migrations)
```

The cause is structural, not tuning. The deadband admits gaps above 5% of mean load, but the SFC re-cut
is restricted to CONTIGUOUS Morton ranges, so its finest correction is one whole tile - **25% of a
rank's load at 4 tiles/rank**. A correction quantum 5x larger than the smallest gap worth correcting can
only overshoot. The escape is self-defeating: correcting at the 5% scale needs ~20+ tiles/rank, and
tiles are exactly what cost the 35%. **The granularity required for useful coarse balancing costs more
than the imbalance it corrects.**

Fixed in passing: the re-cut committed its partition unconditionally, never comparing the resulting gap
to the current one. It now evaluates into a temporary and commits only on strict improvement (the cut
array must move with the owner map - `f_amr_owner` resolves tile ownership against `amr_owner_cut`, so
refreshing one without the other splits ownership silently).

**Scope limit - what this does NOT disprove.** One case class: 2D, single refinement level (the blob
tags level 1 only, confirmed in-log), hydro with uniform per-cell coarse cost, 2 ranks. Coarse work is
near-uniform there by construction, which is *why* there was nothing to recover. Cases where per-cell
coarse cost genuinely varies - IB ghost points, chemistry stiffness, Lagrangian bubbles - are untested
and remain the only place this feature could still pay. The burden of proof has moved: it costs 35% and
has no demonstrated benefit.

**OPEN 2 - DEMOTED by OPEN 1's answer.** It was gated on coarse rebalancing being worth having; it is
not, on the evidence above. Recorded for the case that an imbalanced-coarse-cost case revives it.
`amr_tile_l0_owner` is fixed to the
Cartesian init owner and stays fixed under migration, while COMPUTE ownership follows the cut. Every
migrated tile therefore pays a scatter-back each step through four routines that branch on
`bown == lown`: `s_l0_fill_tiles_from_coarse`, `s_l0_scatter_tiles_to_coarse`,
`s_l0_add_reflux_to_tiles`, `s_l0_restrict_to_tiles`. This is the plan's own product claim - the
coarse grid cannot yet rebalance in production. The `l0_ntile > 0 .and. amr` gate IS lifted, so the
feature runs for dynamic regrid, subcycle, and multi-level.

**A sequencing claim in the original goal is falsified.** It recorded as verified fact that Track 1's
remaining gate and Track 3's redesign were "ONE blocker" and that "Track 1/2 cannot finish underneath
the single-working-slot model". Track 2 finished under exactly that model, and strong scaling improved
1.85x without Track 3.

### Queue

1. **Multi-node weak scaling (8/16/32 ranks, 1/2/4 nodes).** The open exascale question. Nothing
   measured so far leaves one node, so inter-node MPI has never been exercised. Harness:
   `amr-bench/multinode.sbatch` - no uniform control for the scaling claim, oversubscription guard,
   and per-run host verification (GPUs here are not Slurm-managed, and ranks silently packing onto one
   node once produced a fake interconnect cliff).
2. **Sibling-defect audit of the remaining per-block loops.** The hoist found ONE loop-invariant call.
   Three per-block communication sites remain - `s_amr_gather_coarse_patch`, the fine-fine halo, and
   reflux - each entered by every rank per block. This is the only vector that has actually produced a
   large win, and it is cheap.
3. **`amr_max_grid_size` guidance.** The default is 0 = the DERIVED cap, which SHRINKS as ranks grow -
   backwards for strong scaling. Measured: 64 is 3.2x faster than 32 and 10.8x faster than 16 in 3D.
   Documenting the curve is nearly free; changing the default moves every AMR golden and needs its own
   commit.
4. **Then per-block cost**, re-derived at the cap the code should actually run at (the descriptor
   figures were taken at cap 32, where blocks/rank is 4x higher), and checked on a second compiler
   first - the ~626 copies per call may be an amdflang OpenMP mapping artifact rather than structural.

**Correction on limit 3.** An earlier draft of this queue demoted it as "2.4% of the step and scales
down". That measurement was taken at **80 global boxes** and is structurally incapable of seeing the
effect. Balance requires roughly constant blocks/rank, so the GLOBAL box count grows linearly with
rank count, and with it every O(global boxes) term:

| ranks | global boxes (10/rank) | replicated metadata |
|---|---|---|
| 8 | 80 | 0.9 MB/rank |
| 32 | 320 | 3.5 MB/rank |
| 10^3 | 10,000 | 110 MB/rank |
| 10^4 | 100,000 | **1.1 GB/rank** |

Per-block COMPUTE per rank stays constant at fixed blocks/rank - that scales. What grows is the
O(global) work. So limit 3 is not refuted, it is unmeasurable at the scales tested here, and absence
of evidence at 8 ranks is not evidence of absence at 10^4.

**It is testable without 10^4 ranks**, by decoupling global box count from rank count: at fixed np=8
the cap sweep already varies global boxes 1217 -> 80. Run it WITH phase attribution and watch which
phases grow with the GLOBAL box count rather than with blocks/rank - `regrid` growing implicates the
tag allgatherv and assignment, `gather_patch` growing implicates the per-block communication, and
`advance` growing implicates only per-block compute, which is not a scaling problem. That is a
cheaper and strictly more targeted probe than multi-node, which at 32 ranks reaches only 320 global
boxes.

Deliberately NOT next: packing (disproved), launch fusion (8.5% of GPU time).

## What actually blocks exascale: the assignment and regrid cost scale with the GLOBAL box set

**Limit 1 is implemented (pending goldens); limits 2 and 3 remain open.** It was previously recorded only as a diagnosis inside "The
problem, measured" and had no design entry, no queue entry, and no owner. Audited against the code
2026-07-31.

The strategy demands `boxes_per_level >> num_procs`. Three costs in the implementation grew with the
*global* box count, so the strategy and the implementation were in direct contradiction. Limits 1
and 2 are now fixed; limit 3 is open and deliberately not next:

| # | limit | code | 512 boxes (today) | 10^5 boxes |
|---|---|---|---|---|
| 1 | ~~one global collective PER BOX per regrid~~ **HOISTED** | was: `s_set_amr_fine_geometry` reduced inside `do k = 1, nboxes`. Now accumulates into `amr_xchg_bad`; `s_amr_reduce_xchg_flag` reduces ONCE per scan | 512 -> **1** | 10^5 -> **1** |
| 2 | ~~O(n^2) sort in the cut~~ **FIXED** | was: `s_amr_sfc_cut` insertion sort, comment "n small". Now a bottom-up stable merge sort | 1.3e5 -> 4.6e3 ops | 5e9 -> 1.7e6 ops |
| 3 | O(global boxes) replicated metadata per rank | `amr_region_lo_all`, `region_hi_all`, `isect_lo_all`, `isect_hi_all`, `block_owner`, `block_level`, all sized `amr_max_blocks` | 5.6 MB/rank | ~560 MB/rank at 10^7 |

**What the hoist will and will not do.** @ref amr_block_batching measured that batching these
allreduces "buys nothing" at 14-21 boxes, because the 7.4-13 ms per call is absorbing load-imbalance
spread and a barrier collapses it - the time is the wait, not the reduction. That is correct in that
regime and the hoist should NOT be expected to speed up current benchmarks. It removes the O(nboxes)
term: with `nboxes` collectives the cost has a floor of `nboxes x latency` regardless of imbalance,
~0.5 s per regrid at 10^5 boxes before any imbalance at all.

Limit 1 is the worst, and it is worse than box count alone suggests: @ref amr_block_batching measured
that allreduce at **7.4 ms (np=4) to 13 ms (np=8) per call, ~700x a real one-integer allreduce**,
because it absorbs the spread in the owner-only work preceding it. Inserting a barrier collapses the
phase from 0.089 s to 0.0011 s. So its cost grows with rank count as well as with box count.

### Limit 1: how it was removed, and a latent bug it was hiding

The reduction answers "is any block too close to a subdomain edge to prolong its ghosts". Every
per-box answer is immediately OR-ed into one accumulator and only the accumulator survives:

    any_xchg = .false.
    do k = 1, nboxes
        call s_set_amr_fine_geometry(...)          ! ends in s_mpi_allreduce_integer_max
        any_xchg = any_xchg .or. amr_xchg_coarse_ghosts
    end do
    ...
    amr_xchg_coarse_ghosts = any_xchg              ! m_amr_regrid.fpp:1456 - ONLY the OR is kept

So `nboxes` collectives became **one**: each call ORs into the module accumulator `amr_xchg_bad`, and
`s_amr_reduce_xchg_flag` performs a single allreduce to close the scan. All five call sites close
their scan explicitly.

**Two of those call sites were also wrong.** The loops in `s_initialize_amr_module` (L0 tiles) and in
both restart paths kept only the LAST block's answer rather than the OR - an earlier block needing the
coarse-ghost exchange could be masked by a later one that did not, silently skipping an exchange the
fine advance depends on. The accumulator fixes that by construction, since every block ORs in. `amr_xchg_coarse_ghosts` is module state also read in the fine advance
(`m_amr.fpp:4426`, `:4558`), so the flag must still end up global - which is exactly what the hoisted
single reduction produces. Roughly ten lines: split the local test from the reduction and hoist.

### The AMReX model, and which half is adopted

The ownership half is done. The communication half is where the gap is.

| AMReX property | MFC | evidence |
|---|---|---|
| whole-box ownership | yes | always had it |
| absolute small box cap | yes | `amr_max_grid_size` (`3a718392`) |
| per-box scratch, not subdomain-sized | yes | `idwbuff_alloc` (`7b1e4933`) |
| per-level box list AND rank mapping | yes | `cfdd2847` |
| mapping computed redundantly, no communication | yes | `s_amr_sfc_cut` is all-real arithmetic on replicated weights in a fixed order |
| data movement P2P, not collective | yes | `s_amr_gather_coarse_patch`: "Non-participants send/recv nothing (no global collective)" |
| **regrid free of per-box collectives** | **NO** | limit 1 above |

Note what is NOT wrong: the gather is already P2P, and the owner mapping is already computed without
communication. The single defect is the reduction, which was never a design decision - it is a
convenience that is invisible at 512 boxes.

### Sequencing for this arc

1. **Hoist the per-box reduction** (limit 1). Small, mechanically safe, removes the dominant term.
   Correctness bar: byte-identical goldens, since the OR is unchanged.
2. ~~**Replace the insertion sort**~~ **DONE.** Bottom-up merge sort: O(n log n), iterative, and
   stable (`<=` keeps the left run first on ties), so the order is a pure function of the input and
   every rank still produces the identical permutation - which is what the assignment depends on and
   what `s_amr_validate_owner` checks. All 72 AMR goldens byte-identical, as expected: same-level
   boxes have distinct Morton keys, so any correct sort yields the same order.
3. **Only then consider the replicated metadata** (limit 3). It is the least urgent: 5.6 MB/rank at
   10^5 boxes is tolerable, and removing it means giving up the redundant-mapping property that makes
   the assignment communication-free. Do not trade that away without a measurement showing it binds.

None of this is visible on a single node: 512 boxes is three orders of magnitude below where limit 1
bites, which is why the measurements in this document could not have found it and a code audit did.

## The per-block floor: AMR costs ~31x uniform per cell, and the overhead is per-BLOCK

The three limits above are about *scalability* - costs that grow with the global box count. This one
is about *efficiency*, it is a constant factor, and it is much larger than anything else here. A
constant factor of 31 is not fixed by scaling better.

**The measurement.** 75.5M cells, 2D, np=8, `run_time_info=F`, 8 steps, identical grid and IC in both
arms (`amr-bench/cases/hg8_amr.py` and `hg8_ctl.py`); the AMR arm additionally advances a refined
level, so the arms are compared per cell-update, not per step:

| arm | s/step | cell-updates/step | ns per cell-update |
|---|---|---|---|
| uniform | 0.117 | 75.5M | **1.55** |
| AMR (`amr_max_grid_size` 256) | 10.2 | 75.5M + 134.1M fine | **48.7** |

That ratio was unreadable until the balance report started printing `fine_work` (the summed assigned
weight, which with no cost signals is exactly the fine cells advanced). Without it, "AMR is 87x
slower per step" cannot be separated into "it advances 2.78x the cells" and "it pays 31x per cell",
and only the second is a defect.

**The cause is the block count, not the block size.** Same case, same ranks, varying only
`amr_max_grid_size` - which changes how finely the same feature is tiled (`amr_max_blocks` raised to
4096 so it never binds):

| `amr_max_grid_size` | boxes | fine cells | s/step | s per box | ns per cell-update |
|---|---|---|---|---|---|
| 128 | 4096 | 268M | 70.3 | 0.0172 | 204 |
| 256 | 934 | 414M | 40.3 (39.8) | 0.0431 | 82 (81) |
| 512 | 458 | 511M | 12.9 (12.6) | 0.0282 | 22 (21) |
| 1024 | 144 | 761M | 8.5 | 0.0588 | 10 |

Parenthesised values are a repetition of those two rows: 1.3% and 2.9% apart, so everything below is
far outside run-to-run spread.

**Read the first three columns together: the arm that advances the MOST cells is the FASTEST.** 144
boxes over 761M fine cells takes 8.5 s/step; 4096 boxes over 268M fine cells takes 70.3 s/step. Time
rises monotonically with box count while the work done falls monotonically. Cost is set by the number
of blocks advanced, not by the amount of data in them - no curve fitting required to see it.

This confirms at production scale what @ref amr_block_batching measured in the small: a per-block
advance costs a substantial fraction of a full monolithic step regardless of block size, and does not
amortize. Note the per-block cost is not a clean constant either - `s per box` is non-monotonic
(0.0431 at 934 boxes vs 0.0282 at 458), reproducibly so, so it is not simply `a + b*cells`. Something
super-linear in box count is present as well - adjacency-driven work such as fine-fine seam exchange
is the obvious suspect and is not yet isolated.

**Consequences.**

- @ref amr_block_batching's premise is CONFIRMED — per-block overhead is what costs, and the table
  sizes the prize at ~20x between the finest and coarsest tiling of the same feature. Its premise was
  previously argued from efficiency figures since found contaminated (starved GPUs, `run_time_info`
  device syncs), so it had been left open. **But confirming the premise did not validate the planned
  fix:** the packed super-grid was subsequently DISPROVED (2026-08-01), because `s_amr_tile_box`
  splits evenly and so leaves every tiled block larger than half the tile size, making blocks exactly
  slot-sized and \f$P_{\max} = 1\f$. See @ref amr_block_batching, "the packed super-grid cannot work".
- The per-level cap `amr_max_grid_size` is a first-order performance knob, not just a correctness or
  portability one. Small caps are ~20x worse per cell. With packing disproved this is the *only*
  presently-available lever on the per-block floor, and it is bounded by device memory: cap 2048 on
  this case runs the GCD out of memory before the first step completes.
- This interacts with balance in the wrong direction. More boxes per level is what gives the balancer
  freedom (`boxes_per_level >> num_procs`), and it is exactly what costs. The two goals are in
  tension until batching removes the per-block floor, and any future balance work that buys evenness
  by splitting boxes further must be measured against this table.

**What this measurement is not.** The arms do not refine identical areas - a block is a rectangle, so
coarser tiling over-covers, which is why fine cells RISE with the cap. So this is not a clean
single-variable A/B on tiling, and the ns-per-cell column mixes the two effects. The argument
deliberately does not rest on that column: it rests on time and work moving in OPPOSITE directions,
which no confound between them can produce. Repeated rows agree to 1-3%.

## The problem, measured

Cost tracks the **total** number of refined boxes, not the number a rank owns. On a fixed 2D
2048x1024 case with seven refined features, varying only rank count:

| np | `amr_maxc_fit` | boxes | per-rank boxes | s/step |
|---|---|---|---|---|
| 1 | 1024 x 512 | 14 | 14 | 0.6552 |
| 4 | 256 x 512 | 14 | 3.5 | 0.2892 |
| 8 | 256 x 256 | 21 | 2.6 | 0.3479 |

Two things go wrong together. **Box count grows with rank count** — `amr_maxc_fit` is the
min-over-ranks local half-extent, so adding ranks shrinks the cap until a fixed feature must be
tiled into more pieces (a 624-cell-tall feature needs 2 tiles at a cap of 512 and 3 at 256; at
np=32 the cap reaches 128 and it needs 5, giving 35 boxes). And **cost follows total box count
near-linearly**: holding the box set fixed at 35 with `amr_max_grid_size` and comparing against
the default arm at matched rank counts gives time ratios of 2.25x, 2.09x and 1.42x against box
ratios of 2.50x, 2.50x and 1.67x — an exponent of 0.68 to 0.93.

That combination is why np=4 -> 8 regresses. Per-rank boxes *fall* 3.5 -> 2.6, so per-rank work
should drop to ~0.217 s/step; measured it *rose* to 0.3479. The ~1.6x gap is work proportional
to the global box set: every per-box collective in the regrid rebuild loop runs over all boxes on
all ranks.

## Why whole-box ownership is not the problem

An earlier reading of this concluded that a rank must hold an entire block, so per-rank AMR
memory cannot strong-scale, and therefore blocks must be split across ranks. **That conclusion
was wrong.** AMReX, Chombo and BoxLib all keep whole-box ownership — a box belongs to exactly one
rank. They scale because of three properties MFC only partly has:

1. an **absolute, small** box size cap, so boxes are many and small;
2. **per-box scratch**, sized to the box, rather than one working set sized to the subdomain;
3. **per-level distribution** — each level has its own box list *and its own rank mapping*.

Under (1) and (2), per-rank memory is `O(boxes_per_rank * cap^d)`, and `boxes_per_rank` falls as
ranks are added, so it strong-scales. The MFC ceiling was never ownership; it was a cap *derived
from the subdomain* combined with solver scratch *sized to the subdomain*. Those are
`amr_max_grid_size` (landed, `3a718392`) and `idwbuff_alloc` (landed, `7b1e4933`).

## Target design

Adopt (3). Each level keeps its own box list and its own owner mapping, chosen without reference
to the coarse decomposition. Inter-level coupling becomes general point-to-point.

The two rejected alternatives, for the record:

- **Domain-aligned ("mirror") ownership** — a block's fine cells live on whichever rank owns the
  underlying coarse region (named at `m_amr.fpp:326-330`). Parent<->child and reflux become
  rank-local and the communication problem largely disappears, but load balance is then dictated
  by the coarse cut, so a locally-refined region overloads a few ranks. That is precisely the
  problem AMR load balancing exists to solve.
- **Block sub-decomposition** — split a large block across a rank subset. Smallest change from
  today, but it adds a second decomposition layer beneath the existing one and keeps
  block-as-atom, the concept causing the trouble.

## What MFC already has

This is less of a leap than it appears; the branch has been converging on it.

- `s_amr_sfc_cut` / `f_amr_owner` already compute an owner mapping independently of the
  Cartesian coarse decomposition.
- P2P coarse<->fine gather and scatter already exist (`s_amr_gather_coarse_patch`,
  `s_amr_scatter_pbmv`), including the device-side unpack.
- Per-level box lists already exist (`amr_block_level`, `box_level`).
- Level-selective seam halo already exists (`1895c0ef`).

## What was missing — Track 2 (now landed)

Gather, restrict and flux-register delivery were all np=1-local. Multi-rank worked only because a
refinement tower co-located on its level-1 anchor, which made a whole tower the balancer's smallest
atom (weight `cost * rr^(l*d)`) and capped granularity at depth. Under per-level distribution
co-location is not a constraint to be relaxed separately — it disappears as a consequence, which is
why Track 2 was the spine of this phase rather than one item on a queue.

All three paths are now P2P (`6832d299`, `d53fac46`), plus a fourth constraint that was **not** on
the original list and blocked step 4:

- **3d, the per-level lockstep advance** (`26e0d080`). `s_amr_advance_children` drove one parent's
  subtree to completion while the `s_amr_fine_fine_halo(clev)` interposed in its stage loop spans
  EVERY parent, so a seam pair straddling two parents had only one side present. Co-location hid
  the MPI half by putting both ends on one rank. The codebase already knew — `m_checker.fpp`
  fail-closes `amr_subcycle .and. amr_regrid_int > 0 .and. num_procs > 1 .and. amr_max_level > 1`
  citing exactly this. Lifting that PROHIBIT is now unblocked and is the only thing that would
  cover the subcycle SETUP gather added in `cfdd2847`, which currently has **no golden coverage**.

The general lesson, since it cost four bugs: **whenever an invariant says two owners always
coincide, every read that depends on it is load-bearing and invisible.** Grep for the invariant,
not for the symptom.

## Load balancing

Distribution decides *where* a box lives; balancing decides whether that placement spreads work
evenly. Per-level distribution is a prerequisite for the second, not a substitute — this section
states the mechanism, the metric, and the limits that bind at scale.

### Mechanism

Three pieces, only the third of which is AMR-specific:

| Piece | Where | Role |
|---|---|---|
| L0 Cartesian rebalance | `m_load_balance.fpp`, `s_load_balance_rebalance` (called once from `m_start_up`) | weighted splits of the base grid; `load_balance` |
| Cost model + imbalance metric | `m_load_weight.fpp`, `m_rank_timing.fpp` | per-cell weights and a measured per-rank compute time; `load_weight_wrt` |
| AMR block assignment | `s_amr_block_cost` -> `s_amr_sfc_cut` (`m_amr.fpp`) | owner map for level>=1 blocks |

The AMR path weighs each block by `cost(k) * rr**(level*d)`, where `cost` sums a per-cell model
over the block's L0 footprint: base 1, plus `K_ib` per IB-marked cell, plus `K_pc` per
phase-change Newton iteration (`m_constants`: 2 and 3; `K_bub` = 50 applies to Lagrangian bubbles,
which are excluded from blocks by construction). One `MPI_ALLREDUCE(SUM)` makes the vector
identical on every rank, after which the assignment is deterministic and rank-independent.
`s_amr_sfc_cut` then does a chains-on-chains split of the blocks in Morton order of their low
corner, one independent cut **per level**.

**Step 4 is the load-balancing enabler, not merely a distribution change.** Under tower
co-location the balancer's smallest atom was a whole refinement tower, weight `cost * rr**(l*d)`.
No assignment can fix an imbalance whose atom is larger than the imbalance itself, so a single
deep tower pinned work that adding ranks could not relieve. Cutting each level independently is
what makes the cost model actionable.

### The metric

Balance is `max_r W(r) / mean_r W(r)` over per-rank assigned weight `W`, reported per level and
for the total. Two properties are required, and only the second is about scale:

1. **Quality** — imbalance below a fixed tolerance at a given rank count.
2. **Scale invariance** — imbalance must not *grow* with rank count. A scheme that is well
   balanced at np=8 and degrades monotonically by np=1024 has not solved the problem; this is the
   property to test, and it is not visible on a single-node sweep.

`m_rank_timing` supplies the measured counterpart (per-rank RHS + relaxation time), which is the
honest check on the *model*: a cost model that predicts balance while measured times diverge is
wrong, and the model is what the assignment actually uses.

### What binds at scale

These are structural, not tuning knobs, and each sets a floor on achievable balance:

- **Granularity floor.** A level cannot be balanced across more ranks than it has boxes. In
  `s_amr_sfc_cut` a level holding one box always lands on rank 0 (the loop assigns `r = 0` and
  advances only once cumulative weight crosses a share boundary), so a shallow hierarchy leaves
  ranks holding *no block at that level* by construction. Scaling therefore requires
  `boxes_per_level >> num_procs`, which is what the absolute cap `amr_max_grid_size` exists to
  deliver — the cap and the balancer are one mechanism, not two.

  **This is not idleness, and the `[amr-balance]` counter must not be read as if it were.** A rank
  holding no level-`L` block still owns level-0 work — level 0 covers every rank always — and may
  own blocks at other levels. Even the TOTAL line's counter means only "owns no *fine* block". The
  counters are named `no_blocks_ranks` / `ranks_with_no_fine_block` for exactly this reason. Idleness
  is measurable only from `m_rank_timing`; a granularity-floor count is an upper bound on the harm,
  never a measurement of it.
- **Indivisible atom.** The cut is contiguous in Morton order and cannot split a box, so
  imbalance is bounded below by the heaviest single block's weight. As ranks grow, mean per-rank
  weight falls while that floor does not, so imbalance rises unless the cap falls with it.
- **Static within a regrid interval.** Cost is sampled at regrid; work that migrates between
  regrids is not tracked. The relevant knob is `amr_regrid_int`, and its cost is itself
  rank-scaling work (see the box-count analysis above).
- **Cost-model blindness by default.** With no live signals `cost(k)` degenerates to the
  footprint cell count — pure geometry. `pc_iter_count` is populated only when `load_weight_wrt`
  is enabled. A run with heterogeneous per-cell cost and `load_weight_wrt` off is balanced on
  geometry alone, which will look correct and be wrong.
- **Assignment cost.** `s_amr_block_cost` is an allreduce over the global block vector every
  regrid - one call, not one per box (limit 1), and the cut is now O(n log n) (limit 2). What
  remains is that the vector itself is sized by the *global* box set, the same term identified in
  "The problem, measured"; that is limit 3.

### Acceptance criteria

1. Per-level and total imbalance reported at np = 1..N, with imbalance flat or falling in `N`.
2. Model vs measured agreement: `m_rank_timing` per-rank times consistent with assigned weight.
3. A cost-heterogeneous case (IB or phase change) with `load_weight_wrt` on, showing the weighted
   assignment beating the pure-geometry fallback. Geometry-only cases cannot demonstrate this.
4. Deep-tower case: imbalance must not depend on refinement depth once towers may split.

## The cost model is regime-dependent, and production is the cell-dominated regime

Measured 2026-07-31 across three problem sizes. Which metric predicts measured `[rank_time]` flips
with how much work a box carries:

| regime | boxes/rank | box size | measured time tracks |
|---|---|---|---|
| 511x255, `mgs=64`, np=8 | ~2 | tiny | **box count** (1.308 vs measured 1.259; cell weight flat at 1.050) |
| 75.5M, `mgs=256`, np=8 | 13-19 | large | **cell weight** (1.058 vs measured 1.058; box count 1.121) |

Per-block overhead is a FIXED cost - a per-block advance costs ~1x a monolithic step regardless of
size (@ref amr_block_batching). When boxes are tiny that fixed term dominates and load tracks box
count. When boxes carry real work the cell term grows until the two metrics nearly coincide: at 75.5M
with clustered refinement, level-2 weight imbalance is 1.058 and box count 1.121, so the two
objectives barely conflict.

**A prediction that failed, recorded because it constrains the model.** Reading "measured tracks
weight, not box count" at 75.5M, the expectation was that a fixed per-box term would make balance
WORSE by steering toward the metric that does not predict runtime. It did not: `K_box = 8` on the
clustered case moved measured imbalance 1.031 -> 1.022 (np=2, 4) and 1.058 -> 1.038 (np=8), pulling
box-count imbalance 1.121 -> 1.047 while weight imbalance stayed ~1.05. So at scale the two objectives
are close enough that improving one does not cost the other. The regime difference is real but it is a
CONVERGENCE of the two metrics, not a reversal.

Consequence for `K_box` (not landed, see `0e3418ef`): the gain at production scale is ~2% imbalance,
comfortably inside the run-to-run variance measured on this machine (np=4 moved 1.084 -> 1.170 between
two runs of an identical configuration). That does not justify a tunable constant that also requires
the ULP cut fix as a prerequisite. It remains worth revisiting if a future regime pushes
boxes-per-rank back down - which `amr_max_grid_size` does directly.

Three measurement flaws had to be removed before any of this was visible, each of which produced a
confident wrong conclusion first:

1. **Starved GPUs.** 511x255 at np=8 is 16k cells/rank; an MI250X GCD needs O(1e6) to be compute
   bound. `mfcrun.sh` now refuses below 1e5 cells/rank.
2. **`run_time_info=T`** forces a device->host sync and a global reduction EVERY step. It made a
   uniform 75.5M control look like 6% parallel efficiency at np=8. Now off in every generated case.
3. **A low-discrepancy IC.** `gen_blobs.py` places blob centres with a Weyl sequence - evenly spread
   BY DESIGN (24 blobs land 3,3,2,4,3,2,4,3 across eight x-slabs). The workload is inherently
   balanced, so the case cannot discriminate between cost models. `gen_clustered.py` concentrates the
   refinement instead.

And a fourth, which is a result rather than a flaw: **spatial clustering does not create imbalance**,
because per-level distribution decouples ownership from position. With all refinement in one quadrant
the balancer still spreads boxes over every rank (`no_blocks_ranks 0 of 8`, imbalance 1.03-1.06).
That is step 4 working as designed.

## Sequencing

Steps 1-4 are LANDED (see "Status"); they are kept here because the ordering constraints between
them are load-bearing and a reader retracing the work needs them. Steps 5-7 are the live queue.

### Landed

1. **Finish the scratch decoupling.** Convert the `m`/`n`/`p`-keyed allocation family that an
   `idwbuff` grep does not find: `m_riemann_solvers` sizes on `-1:m,-1:n,-1:p` and `m_weno` on
   `is*_weno` (`is1_weno%%end = m - is1_weno%%beg`). **Only then** relax the
   `amr_max_grid_size > fit_d` abort in `m_amr.fpp`. Relaxing it first writes out of bounds
   silently rather than failing to compile.
2. **Measure with the cap pinned below the derived cap at every rank count.** This is what shows
   the box set going rank-invariant and the np=8 turnover softening. Run the uniform control arm
   at every point; without it an AMR curve at 65k cells/rank is unreadable.
3. **Generalise parent<->child to P2P**, level by level, with the level mapping still equal to
   the coarse one — behaviour-preserving, so the existing goldens gate it.
4. **Give each level its own mapping** and drop tower co-location.

   Splitting towers made four latent reads reachable, all one shape: a **parent-slot field
   guarded on owning the CHILD**, safe only while the two owners always coincided. Undefined
   `amr_slots(pblk)%%amr_ref_ratio`, a bisected unallocated `%%x_cb`, and `lbound`/`ubound` of that
   same unallocated array. Symptom was garbage cell widths and NaNs at the rank seam, several
   steps later. Anything reading another block's slot must be guarded on owning **that** block, or
   derive from replicated metadata (`s_amr_parent_foot`, `s_amr_build_block_coords`).

### Next

5. **Instrument balance.** LANDED (`fc53e097`, fixed in `3780f30a`). `s_amr_report_balance` prints
   per-level and total `max/mean` assigned weight, box-count imbalance, and the no-block rank count,
   gated behind `load_weight_wrt`. Needs no MPI - the inputs are replicated, so rank 0 prints.
   Its first conclusion ("box supply binds, not the owner mapping") was RETRACTED: it was inferred
   from box counts and model imbalance on a starved case, with no measured counterpart.
   `m_rank_timing` was already implemented and wired the whole time - acceptance criterion 2 was a run
   to perform, not code to write.
6. **Model vs measured (acceptance criterion 2).** DONE. At production scale the model predicts the
   measurement: 1.000 model vs 1.006-1.058 measured across np = 2, 4, 8 at 75.5M cells. On starved
   cases they diverged 8x, which is what motivated step 7.
7. **Reweight `s_amr_block_cost` on a fixed per-box term.** ATTEMPTED, MEASURED, **REJECTED.**
   `K_box` (fixed per-box cost in mean-block-cell units, added after the allreduce, fine blocks only)
   improved balance substantially on small cases and delivers ~2% at production scale - inside this
   machine's run-to-run variance. Not worth a tunable constant that also requires the ULP cut fix as a
   prerequisite. Removed in `0e3418ef`; see "The cost model is regime-dependent" below for the
   measurements and for a prediction of mine that failed.

   **The experiment paid for itself anyway.** Perturbing the weights exposed two latent partitioner
   bugs, both fixed independently and both of which had survived only because the arithmetic happened
   to be benign:
   - `2051aa19` - `s_amr_sfc_cut` compared an n-term accumulation against a closed-form target, so an
     exact share boundary turned on 1 ULP. Correct only because every cost term is integer-valued.
   - `c3364a5b` - uninitialized L0 tile-prefix slots entered the fine-level cut as phantom key-0
     blocks. Correct only because they all happened to land on rank 0.

### The live queue

Two questions are now settled, and they narrow the remaining work considerably.

**Balance is solved at production scale.** At 75.5M cells with refinement clustered in one quadrant,
the balancer holds 1.03-1.06 measured imbalance and puts blocks on every rank. Spatial concentration
does NOT create imbalance, because per-level distribution decouples ownership from position - step 4
working as designed. The cell-based cost model is adequate there; an experimental fixed per-box term
bought ~2%, inside run-to-run variance, and was rejected (`0e3418ef`).

**The scaling ceiling is MFC's, not this branch's and not AMR's.** Same portable uniform case
(`amr-bench/gen_uniform.py`, no AMR in the path), same harness, run on MFC master `bfdc8f5e` and on
this branch, 3 reps alternating between trees:

| np | master median | branch median | ratio | within-tree spread |
|---|---|---|---|---|
| 4 | 0.3285 s | 0.3359 s | **1.022** | 6.9% / 2.0% |
| 8 | 0.2326 s | 0.2380 s | **1.023** | 11.1% / 5.2% |

The branch costs ~2% on a uniform problem - indistinguishable from noise at an 11% within-tree
spread. Uniform efficiency is ~60% at np=8 on both trees. **No AMR work can move that**, and a
single-run comparison earlier suggested 15% at np=4 purely as an artifact: three reps put it at 2%.

Consequence for measurement discipline: at 11% run-to-run spread on this machine, **any A/B claiming
less than ~10% needs repetitions**, taken alternating between arms so drift hits both equally.

8. ~~**Hoist the per-box reduction in the regrid rebuild loop.**~~ **DONE** (`3f754893`), and limit 2
   with it (`83bf4572`).
9. ~~**What does AMR actually cost relative to uniform, at a loaded size?**~~ **ANSWERED: ~31x the
   per-cell cost of uniform**, and the overhead is per-BLOCK, not per-cell. See "The per-block floor"
   below - this is now the largest single number in this document.
10. ~~**Then, and only then, revisit @ref amr_block_batching.**~~ **Its premise is CONFIRMED**, and by
   a measurement rather than by the contaminated efficiency figures it originally rested on. Per-block
   overhead does dominate; item 9 sizes the prize at ~20x. **The planned fix, however, is DISPROVED:**
   the packed super-grid needs blocks smaller than half the cap, and the even-split tiler guarantees
   the opposite, so \f$P_{\max} = 1\f$ on every measured configuration. Batching survives only via the
   flat backing store and its unresolved `ACC_SETUP_SFs` aliasing risk. Confirming that a cost is
   per-block says nothing about whether a given mechanism can remove it.

13. **A proper 3D scaling benchmark.** The 2D 75.5M pair is retired for scaling work: it cannot
    strong-scale (per-rank AMR memory is ~refined-volume/ranks, so it OOMs at np <= 4 at most caps —
    memory pressure rises as ranks are REMOVED) and its refinement is not reproducible across rank
    counts (@ref amr_block_batching, "the box set is not rank-invariant"). The replacement is a 3D
    sharp-interface case, verified rank-invariant first. The current one (128^3) is still small — at
    np=8 both arms sit near the starvation floor — and forms no level-2 blocks, so it exercises a
    single refinement level. Deepen and enlarge it before treating it as the primary benchmark.

    First valid AMR-vs-uniform strong scaling, identical box set at every rank count (16 level-1
    boxes, `fine_work` 691200, imbalance 1.000):

    | np | AMR s/step | uniform s/step | AMR eff | uniform eff | AMR ns/cell-upd | uniform ns/cell-upd | ratio |
    |---|---|---|---|---|---|---|---|
    | 1 | 1.132 | 0.1142 | 1.000 | 1.000 | 405.9 | 54.4 | 7.46 |
    | 2 | 0.712 | 0.0841 | 0.795 | 0.679 | 255.3 | 40.1 | 6.37 |
    | 4 | 0.656 | 0.0714 | 0.431 | 0.399 | 235.3 | 34.1 | 6.91 |
    | 8 | 0.454 | 0.0679 | 0.312 | 0.210 | 162.7 | 32.4 | 5.03 |

    Two results worth keeping. **AMR scales BETTER than uniform at every rank count here** (0.431 vs
    0.399 at np=4), because at 2.1M cells the uniform arm starves first — AMR carries more work per
    rank and hides latency longer. So AMR is not the scaling limiter at these sizes; problem size is.
    And the per-cell overhead of 5.0–7.5x independently reproduces the ~7.3x measured in 2D at the
    best cap, from a different case, dimension, and cap — the first time that figure has been
    confirmed twice.

14. **THE SCALING LIMITER IS CROSS-RANK BLOCK ADJACENCY, not the per-block constant.** Same 3D case
    enlarged to 256^3 (16.8M base), rank-invariant (64 L1 + 256 L2 boxes, `fine_work` 15160320
    identical at every rank count), np=1 device-OOMs so efficiencies are normalized to np=2:

    | np | AMR s/step | uniform s/step | AMR eff | uniform eff | AMR ns/cell-upd | uniform ns/cell-upd | per-cell |
    |---|---|---|---|---|---|---|---|
    | 2 | 13.121 | 0.3146 | 1.000 | 1.000 | 410.8 | 18.75 | 21.9x |
    | 4 | 9.071 | 0.1763 | 0.723 | 0.892 | 284.0 | 10.51 | 27.0x |
    | 8 | 7.544 | 0.1281 | 0.435 | 0.614 | 236.2 | 7.63 | 30.9x |

    Against the 128^3 case (16 blocks) run identically, AMR there scaled BETTER than uniform and its
    per-cell overhead SHRANK with ranks (7.46 -> 5.03). Here, with 320 blocks, AMR scales WORSE than
    uniform and per-cell overhead GROWS (21.9 -> 30.9). Same code, machine, cap and dimension; the
    only variable is block count.

    **CORRECTION - the mechanism is NOT growing communication.** This section first attributed the
    rising overhead to cross-rank block adjacency. Direct attribution with `rank_time_wrt` (which
    accumulates per-rank RHS + relaxation compute behind a device sync) refutes that. Same case, same
    ranks, 20 steps:

    | np | s/step | compute s/step | non-compute s/step | non-compute share | compute efficiency |
    |---|---|---|---|---|---|
    | 2 | 13.048 | 5.472 | 7.576 | 58.1% | 1.000 |
    | 4 | 9.003 | 2.794 | 6.209 | 69.0% | 0.979 |
    | 8 | 7.227 | 1.447 | 5.780 | 80.0% | 0.946 |

    **The solver scales almost perfectly (0.946 at np=8). The AMR-side work barely distributes:**
    7.58 -> 5.78 s/step for 4x the ranks, a 1.31x reduction against an ideal 4x, rising from 58% to
    80% of runtime. So the overhead is not growing with rank count - it is roughly CONSTANT in
    absolute time and fails to parallelize. The per-cell overhead trend (21.9 -> 30.9x) is that fixed
    cost becoming a larger share as compute shrinks. This is Amdahl, not a communication blow-up, and
    it caps AMR strong scaling regardless of how cheap the per-block advance becomes.

    **Mechanism, measured with `rocprofv3` - and it is NOT kernel launch count.** One profiled run in
    the 320-block steady state at np=2, tracing kernels and memory copies together:

    | stream | count | GPU time | share of span | mean |
    |---|---|---|---|---|
    | kernel execution | 143290 | 9.91 s | 8.5% | 69.2 us |
    | device-to-device copies | **2937729** | 28.60 s | 23.5% | 9.7 us |
    | GPU idle (host / MPI) | - | ~82 s | **~68%** | - |

    Mean kernel duration is 69 us, so these are not launch-latency-bound kernels, and 143k launches at
    a generous 10 us each is ~1.4 s - it cannot explain ~82 s of idle. The striking number is **20.5
    device-to-device copies per kernel**, ~306 per block per RK stage, costing 2.9x more GPU time than
    all kernel execution combined. Ranking: host/MPI serialization first, small device-to-device
    copies second, compute a distant third.

    **This undercuts @ref amr_block_batching's core premise.** That arc ranks increments by launches
    removed, on the strength of a 16-tile / 1-fine-block measurement where launch count was the cost.
    At 320 blocks launch count is not in the top two terms, so launch fusion cannot address the
    limiter. Re-derive the increment ranking against copies and host time before spending on it.

15. **THE ACTUAL LIMITER WAS A LOOP-INVARIANT CALL, and the biggest lever is a PARAMETER.** Two
    results supersede much of the framing in item 14, both found by top-down phase attribution rather
    than by reasoning about mechanisms.

    **(a) `s_amr_exchange_coarse_cons_halo` ran once per BLOCK.** It sat inside `s_amr_fine_stage_fill`
    (lock-step) and `s_amr_subcycle_setup_block` (subcycle, twice - two lerp sources), both called per
    block per stage. It exchanges the COARSE field, which every fill READS and none WRITES, so it is
    loop-invariant: ~960 identical whole-subdomain exchanges per step where 3 suffice. Because the
    count tracked the GLOBAL block count while each exchange costs the same, it did not distribute -
    adding ranks made it absolutely SLOWER (phase efficiency 0.155, 58% of the step at np=8).
    Hoisting it to the two call sites gives, on the 256^3 rank-invariant case:

    | np | before | after | strong-scaling efficiency |
    |---|---|---|---|
    | 2 | 13.48 | 10.37 | 1.000 -> 1.000 |
    | 4 | 9.54 | 6.00 | 0.707 -> 0.864 |
    | 8 | 7.91 | 3.28 | **0.426 -> 0.790** |

    `fine_work` is identical at every rank count, and 61 AMR goldens are byte-identical. AMR now
    strong-scales BETTER than the uniform control (0.790 vs 0.432, the control starving first).

    **(b) `amr_max_grid_size` is worth more than the code fix.** Sweeping it at np=8 (2 reps,
    `amr_max_blocks` pinned so it cannot bind):

    | cap | boxes | blocks/rank | idle ranks | imbalance | s/step | ns/cell-update |
    |---|---|---|---|---|---|---|
    | 16 | 1217 | 152.1 | 0 | 1.017 | 11.82 | 428.1 |
    | 32 | 320 | 40.0 | 0 | 1.000 | 3.47 | 108.6 |
    | 48 | 180 | 22.5 | 0 | 1.071 | 2.17 | 65.5 |
    | **64** | **80** | **10.0** | **0** | **1.000** | **1.09** | **31.7** |
    | 96 | - | - | - | - | device OOM | - |

    Monotone to 64, then a MEMORY wall (confirmed `__tgt_target_data_begin_mapper`; slot and solver
    scratch both go as cap^num_dims). **Cap 64 is 3.2x faster than 32 - which is AMReX's 3D DEFAULT -
    and 10.8x faster than 16**, with perfect load balance at the fastest point. So granularity, not
    balance, is the binding constraint in 3D here; AMReX's own GPU guidance (`max_grid_size` 128) says
    the same. Against the accuracy-matched baseline (brute-force uniform 1024^3 at 8.19 s/step) AMR
    captures 3% of the available 33.6x benefit before the fix, 7% after, and **22% at cap 64**.

    **What this demotes.** Item 14's "replicated metadata / Limit 3" hypothesis is WRONG: regrid is
    2.4% of the step and scales DOWN. The per-invocation descriptor traffic is real and is now the
    largest GPU-side term (66% of GPU time), but at cap 64 there are 10 blocks/rank rather than 40, so
    it matters ~4x less than the cap-32 measurements suggested. Re-derive its value at the cap the code
    should actually run at before committing to the flat backing store.

    **ROOT CAUSE - `s_compute_rhs` has a FIXED per-invocation cost, and AMR pays it per BLOCK.**
    Profiling the uniform control on the identical base grid isolates it. Per `s_compute_rhs` call the
    cost is ~15-31 kernel dispatches and ~305-609 device-to-device copies, essentially independent of
    how many cells that call covers. Uniform pays it 3x per step; AMR pays it 963x (320 blocks x 3
    stages + coarse):

    | | uniform | AMR (320 blocks) | ratio |
    |---|---|---|---|
    | RHS invocations / step | 3 | 963 | **321x** |
    | kernel dispatches (10 steps) | 936 | 143290 | 153x |
    | memory copies (10 steps) | 18260 | 2937743 | 161x |
    | **copy GPU time** | **0.26 s** | **28.53 s** | **109x** |
    | kernel GPU time | 2.58 s | 9.86 s | 3.8x |
    | mean kernel duration | 2761 us | 69 us | 0.025x |

    So the per-block floor is not a mystery constant: it is the RHS machinery's fixed setup, paid once
    per block per stage rather than once per stage. The copies dominate - negligible at 0.26 s in the
    uniform run, they become the largest single GPU consumer at 28.5 s, and with their HSA signalling
    (`hsa_amd_memory_async_copy_on_engine` fires 2928329 times) they account for ~88% of all 33.2M host
    GPU-API calls.

    **What the copies ARE - CONFIRMED by varying `sys_size`.** `q_cons`/`q_prim`/`rhs` are arrays of
    `scalar_field` with pointer `%%sf` components, so the suspicion was OpenMP map/descriptor traffic
    rather than bulk data movement. Test: 3D uniform control, `num_fluids` 1 -> 2 (`sys_size` 6 -> 8),
    everything else fixed.

    | `num_fluids` | `sys_size` | kernels / call | copies / call |
    |---|---|---|---|
    | 1 | 6 | 30 | 626 |
    | 2 | 8 | 30 | 719 |

    **Copies scale with `sys_size`; kernel count does not** (the `sys_size` loops live INSIDE the
    kernels). Solving the two points gives **~46.5 copies per field variable plus a ~347 fixed base**.

    That correlation alone does NOT identify descriptors - any per-field operation in a
    `do i = 1, sys_size` loop scales identically. Two further tests settle it. Varying grid size at
    fixed `sys_size` (63^3 / 127^3 / 255^3, a 64x volume change) gives **626 copies per call at every
    size - exactly constant** - so the count is not data-proportional. And the duration distribution at
    255^3 is sharply bimodal:

    | bucket | count | share of copy time |
    |---|---|---|
    | ~5.4 us (p50) | 11243 (99.7%) | 29.4% |
    | >1 ms | 29 (0.3%) | 70.6% |

    The 29 millisecond-scale copies are bulk data (initial load and save, ~5 per step, not per call);
    they are what made the MEAN appear to grow with grid size. The per-call population is 99.7% flat
    ~5.4 us copies, p10-p99 spanning only 5.0-7.0 us across that 64x volume change - latency-bound, so
    bytes not megabytes. Count fixed in volume, duration flat in volume, count linear in `sys_size`:
    that is metadata/descriptor traffic per (field x mapped argument), ~3.4 ms per `s_compute_rhs`
    entry. Uniform pays it 3x/step (~10 ms, invisible); AMR pays it 963x.

    **Why it is addressable:** `sys_size` is fixed once `s_read_input_file` completes, so every array
    shaped by it has known, unchanging extents from initialization onward. Re-establishing these
    mappings on every RHS entry is redundant work, not an inherent cost.

    **Consequence for the fix.** Making these mappings persistent instead of per-call is the lever, and
    it is worth noting this is the same flat-backing-store restructuring
    (`%%sf => amr_qall(:,:,:,i,k)`) that @ref amr_block_batching proposed for an unrelated reason
    (enabling batched kernels). The batching rationale is dead (see "the packed super-grid cannot
    work"); the descriptor-traffic rationale is live and is supported by direct measurement. Its
    `ACC_SETUP_SFs` aliasing risk on Cray is unchanged and still needs a build to settle.

    **SIZING CAVEAT - these magnitudes are provisional.** Measured density on this machine is 19.34 GiB
    for 16.8M uniform 3D points = 0.87M points/GiB, so a 64 GiB GCD holds ~56M points. This case runs
    2.1M points/rank at np=8, about **4% of capacity**, which inflates the non-compute share (compute
    is tiny while AMR overhead is roughly fixed). The DIRECTION - AMR-side work does not distribute -
    is robust, but every magnitude above needs re-measuring at ~30-40M points/rank (roughly 640^3 for
    the uniform control, with headroom for AMR slots) before being quoted. Profiling also costs ~16%
    (15.9 vs 13.7 s/step), inflating the idle fraction somewhat, though not enough to change the
    ranking.

12. **Where does the residual ~7.3x go at the best cap?** Even at cap 1024 - the largest that fits
    device memory here - AMR is 11.4 ns/cell-update against uniform's 1.55. Packing cannot close it
    and raising the cap further OOMs, so this is the standing efficiency question. Candidates, in the
    order they are worth measuring: remaining unfused ghost-slab loops (`s_amr_lerp_fine_ghosts` is
    the valuable one - it runs per SUBSTEP, not per stage), the super-linear-in-box-count term from
    item 9, and per-block occupancy at the RHS kernels themselves.
11. **Multi-node scale invariance.** This is the actual exascale question and single-node np<=8 cannot
    answer it: the granularity floor and the indivisible atom only bind once ranks approach the box
    count per level. Everything measured here tops out at 8 ranks with 13-64 boxes per rank - a
    regime where the balancer is comfortable. The interesting regime is the one where it is not.

## Superseded measurements, and why they are kept out of the body

Three earlier sections were removed in the 2026-07-31 reorganisation because their numbers are known
to be contaminated. They are recorded here so a reader meeting them in git history knows not to trust
them, not because they are still evidence.

- **"Balance is no longer the limiter"** quoted 62% parallel efficiency at np=4 and 38% at np=8, and
  concluded that per-box overhead had to be the remaining cost. Both figures came from runs with
  `run_time_info = T`, which forces a device->host sync and a global reduction every step. The same
  configuration measured 6% efficiency at np=8 on a *uniform* 75.5M-cell control, i.e. the figure was
  measuring a diagnostic barrier.
- **"The cost model weighs the wrong quantity"** reported measured time tracking box count (1.308)
  while cell weight stayed flat (1.050) at np=8. That is real, but only in the starved regime: at
  511x255 a rank holds ~2 boxes of 16k cells. See "The cost model is regime-dependent" for the
  measurement at production size, where the two metrics converge.
- **The step-4 A/B "Measured outcome"** (2047x1023, co-location vs per-level) was run at a size where
  the uniform control reached 18% parallel efficiency at np=8, so both arms were latency-bound
  exactly where the effect should have appeared. It cannot distinguish "no imbalance to fix" from
  "balancer did nothing".

The general lesson, which cost most of a day: **a timing number is only as good as the configuration
that produced it.** Check `run_time_info`, cells per rank, and the uniform control's own efficiency
before interpreting any of them.

## Validation

The bar this branch already uses: np=1 bit-identical, np>=2 conservation-exact, plus the goldens'
tolerance compare. Cross-compiler coverage (CCE, gfortran, nvfortran, ifx) is CI's job and runs on
the PR — the **local** gate is correctness under OMP GPU offload with the AMD AFAR compilers. Do
not hold work waiting on a Frontier run. Expect CCE to be where a regression in this area first
appears, though: every one so far has been CCE-only, so treat a CCE failure as a real bug rather
than flaky infrastructure.

For steps 3 and 4, prefer the **conservation ladder** over field diffs: authoritative mass after
each stage, the pure reflux delta, coarse-equivalent fine mass per block, and covered-cell mass
per block, all allreduced and printed from both arms. That is what localised the np>=2 restrict
bug (`ab87d49e`) when field diffs could not; many-to-many coupling is exactly where that class of
bug hides.

Balance is validated separately from correctness, and neither substitutes for the other. A run can
be conservation-exact and golden-clean while every fine block sits on one rank; nothing in the
golden suite measures distribution quality, because the goldens are single- and two-rank
tolerance compares of field data. Use the acceptance criteria in "Load balancing" above, and treat
end-to-end s/step as a *consequence* of balance rather than evidence of it — a case whose towers
happen to spread evenly will scale well regardless of whether the balancer did anything.

Note that `./mfc.sh test --only AMR` does not match the coexist goldens — their trace token is
"AMR + L0 tiles". Run `1F074C5D 8D466A94 83CC5C6D 33060D84 FD056B71 98AA6EDB D99F85F8 09E0D257
93EFC4F6` by UUID as well, or they silently go untested. (A `-l | grep -i amr` filter does catch
them, since the token still contains "AMR"; the suite is 71 cases as of `1e07eb65`.)

**A new golden must be shown to fail without the change it protects.** `3db24df0` set this bar and
it has caught real self-deception since: two attempted counterfactuals for the level-3 slot cap
produced byte-identical output and proved nothing, because the case was too small and then because
only the grid was scaled and not `amr_buf` — boxes track the *feature*, not the domain. The third
attempt crashed the old code outright. Where a path genuinely cannot be covered, say so in the test
comment with the reason, as `C45DBB52` does for the level-2 seam over MPI; silent non-coverage reads
exactly like coverage a year later.