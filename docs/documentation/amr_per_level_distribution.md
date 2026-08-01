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

8. **Hoist the per-box reduction in the regrid rebuild loop.** THE PRIMARY WORK ITEM - see "What
   actually blocks exascale" above. `nboxes` global collectives per regrid become one; the per-box
   results are already discarded. ~10 lines, byte-identical goldens expected.
9. **What does AMR actually cost relative to uniform, at a loaded size?** STILL UNKNOWN and the
   gating measurement. The AMR np=1 point at 75.5M exceeded the harness's 1800 s cap, so there is no
   baseline and therefore no overhead ratio. Raise the cap or cut the step count for that point.
   Without this number there is no way to size the prize for any AMR-side optimisation.
10. **Then, and only then, revisit @ref amr_block_batching.** Its premise - that per-block overhead
   dominates - was argued from efficiency figures since found contaminated (starved GPUs at 16k
   cells/rank, and `run_time_info` forcing a device sync every step). It may still be right, but
   against a ~60% uniform ceiling the recoverable headroom is smaller than the original framing
   assumed, and step 8 is what sizes it.
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