# AMR per-level distribution

Design for the next phase of block-structured AMR: distribute each refinement level
independently, so AMR strong-scales instead of degrading as ranks are added. Written
2026-07-29 against `up/mega` @ `7b1e4933`. Companion to @ref amr_block_batching, which covers
the per-block launch cost, and @ref amr_multilevel, which covers the nesting.

## Status

All four design steps are landed. The measured outcome did **not** confirm the performance thesis,
and the reason matters more than the steps did — see "Measured outcome".

| Step | Commit | Result |
|---|---|---|
| 1. Scratch decoupling | `86782249` | landed |
| 2. Rank-independent cap | `a108dd37` | np=4->8 flipped 1.16x slower to 0.84x faster (single level) — **but see the caveat below** |
| 3. P2P parent<->child | `6832d299`, `d53fac46` | landed; 3d `26e0d080` hoists the level advance to lockstep |
| 4. Per-level mapping | `cfdd2847` | landed; correctness proven with a split tower, **no measured speedup — now explained**, see "The cost model weighs the wrong quantity" |
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
  regrid, and the cut is O(nblocks^2) in the insertion sort. Both grow with the *global* box set,
  the same term identified in "The problem, measured".

### Acceptance criteria

1. Per-level and total imbalance reported at np = 1..N, with imbalance flat or falling in `N`.
2. Model vs measured agreement: `m_rank_timing` per-rank times consistent with assigned weight.
3. A cost-heterogeneous case (IB or phase change) with `load_weight_wrt` on, showing the weighted
   assignment beating the pure-geometry fallback. Geometry-only cases cannot demonstrate this.
4. Deep-tower case: imbalance must not depend on refinement depth once towers may split.

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

5. **Instrument balance.** LANDED (`fc53e097`, fixed in `3780f30a`): `s_amr_report_balance` prints
   per-level and total `max/mean` assigned weight plus the no-block rank count, gated behind
   `load_weight_wrt`. It needs no MPI — `wt`, `amr_block_level` and `amr_block_owner` are
   replicated, so rank 0 prints.

   First result: at np=2 the metric is already ~1.005-1.05, i.e. step 4 had no headroom to recover;
   at np=8 the benchmark collapses to ~2 boxes per level and both distribution schemes produce
   byte-identical assignments. **Box supply binds, not the owner mapping** — which points at
   `amr_max_grid_size` and regrid box production, not at distribution.

   **Correction (supersedes an earlier revision of this file).** A previous version of this entry
   said "the `m_rank_timing` half is NOT done". That was wrong: `m_rank_timing` is fully implemented
   and wired — `s_rank_time_tic`/`toc`/`s_report_rank_time`, gated on `rank_time_wrt`, with the AMR
   fine advance bracketed at eight call sites in `m_amr.fpp` under `amr_rank_owns_block`, and the
   allreduce already printing `[rank_time] imbalance(max/mean)`. Acceptance criterion 2 is therefore
   not work to be written; it is a run to be performed. Harness: `amr-bench/model_vs_measured.sh`.
6. **Model vs measured (acceptance criterion 2).** DONE — and it overturned the working thesis. See
   "The cost model weighs the wrong quantity" below. Harness: `amr-bench/model_vs_measured.sh`.
7. **Reweight `s_amr_block_cost` on a fixed per-box term.** MEASURED GOOD, **NOT LANDED — blocked by
   a coexist regression.** `K_box = 8` in
   mean-block-cell units, added after the allreduce. MEASURED on `[rank_time]`, same case:

   | np | measured, `K_box=0` | measured, `K_box=8` | excess removed |
   |---|---|---|---|
   | 2 | 1.087 | **1.020** | 77% |
   | 4 | 1.170 | **1.012** | 93% |
   | 8 | 1.259 | **1.144** | 44% |

   `K_box=8`, not 4: a CPU sweep of assignment quality (`amr-bench/sweep_kbox.sh`) shows level-2 box
   imbalance at 1.362 for `K=0`, **1.106 across {1,2,4}**, and **1.021 across {8,16}** — two plateaus,
   saturating at 8. An earlier value of 4 sat in the lower plateau and bought nothing over `K=1`.

   **Requires the ULP cut fix (`2051aa19`) as a prerequisite.** Without it `K_box` makes the weights
   non-representable and level 1 breaks from 1.000 to 1.250 — the reweighting was unshippable until
   the partitioner stopped being tie-fragile.

   Caveat: one run per point, and np=4 measured 1.084 vs 1.170 across two baseline runs, so variance
   is real. The improvement holds at all three rank counts, which is stronger than any single point,
   but the magnitudes are not settled.

   **It exposed a pre-existing bug, now fixed: phantom tile-prefix slots in the fine cut.** Four
   coexist (L0 tiles + AMR) goldens at np=2 — `8D466A94`, `33060D84`, `98AA6EDB`, `09E0D257` — aborted
   with `SFC cut-point owner disagrees with amr_block_owner`. Dumping the block table showed why:

       fine_cut(:,1) =      0   33792        owner_cut =     -1      -1   (tile cut NOT built yet)
       blk 1..8   lev 1  key 0      lo (0,0,0)     <- the 8 tile-prefix slots, UNINITIALIZED
       blk 9      lev 1  key 5120   lo (16,8,0)    <- real fine block
       blk 10     lev 1  key 33792  lo (32,8,0)    <- real fine block

   Fine blocks occupy slots `(l0_slot_off, amr_num_blocks]`; slots `[1, l0_slot_off]` are the L0 tile
   prefix. At init the assigner runs BEFORE `s_l0_tiles_init`, so those prefix slots are still
   uninitialized — level reads 1, `region_lo` is all zeros, i.e. **Morton key 0**. The assigner was
   looping from slot 1, so it fed eight phantom key-0 "level-1 blocks" into the level-1 cut and split
   them across ranks. **A key-0 block can only ever resolve to rank 0** — the cut is non-decreasing and
   the search returns the first `r` with `key <= cut(r)` — so any phantom placed on a higher rank is
   unrecoverable, and the validator correctly aborts.

   **This was latent, not new.** With the old weights all eight phantoms happened to land on rank 0 and
   the validator agreed by luck; `K_box` moved the cut and exposed it. Same shape as the ULP tie
   (`2051aa19`): a fragility that survived only because the arithmetic happened to be benign.

   Fix: the assigner and the validator both restrict to `(l0_slot_off, amr_num_blocks]`, with the
   validator skipping the prefix only while `amr_owner_cut` is unbuilt (the tile-init call site
   populates both cuts and validates there).

   **`s_amr_validate_owner` is load-bearing.** It is marked TRANSITIONAL in-source ("removed once the
   table is deleted"); it is the only thing that caught this, and it caught it immediately. Do not
   delete it while the cut and the owner table are both live.

   Two things this already establishes, independent of the eventual fix:

   - **The benchmark case never exercises coexist.** The np sweep that produced the numbers above runs
     plain AMR, so the defect was invisible to every measurement and visible only to the suite. Run the
     goldens *before* reporting a performance result, not after.
   - **`s_amr_validate_owner` is load-bearing.** It is described in-source as TRANSITIONAL ("removed
     once the table is deleted"); it is the only thing that caught this, and it caught it immediately.
     Do not delete it while the cut and the owner table are both live.

## Balance is no longer the limiter, and that redirects the effort

At np=4 measured imbalance is **1.012** — 98.8% of perfect — while parallel efficiency on the same
run is **62%** (2.46x on 4 ranks; `t_max` 41.4 -> 16.8 s). At np=8, imbalance 1.144 and efficiency
**38%**. Imbalance accounts for almost none of the loss.

So the remaining cost is a non-scaling term that distribution cannot touch, and it is already
measured: a per-block advance costs ~1.05x a full monolithic step **regardless of block size** and
**does not amortize with problem size** (16 tiles = 16.75x monolithic; cost is linear in block
count). That is the same fixed per-box quantity `K_box` had to model in order to balance correctly.

**The main line moves to @ref amr_block_batching.** Steps 8-9 below remain useful but are refinements
to a balancer that is now close to optimal on this class of case; they are not where the machine is
being lost.
8. **Exercise the heterogeneous terms.** A cost-heterogeneous benchmark (IB or phase change) with
   `load_weight_wrt` on. Every AMR benchmark to date is geometry-only, so `K_ib`/`K_pc` have never
   influenced a measured assignment — and they are corrections to a base term that is itself wrong.
9. **Scale-invariance run.** Sweep rank count past the box count per level and confirm imbalance
   does not grow. Single-node np<=8 cannot show this — the granularity floor and the indivisible
   atom only bind once ranks approach box count.

## The cost model weighs the wrong quantity

Measured 2026-07-31 with `load_weight_wrt` and `rank_time_wrt` both on. 511x255, 24 blobs (`hcid`
299), `amr_max_level=2`, `amr_max_grid_size=64` (pinned, so boxes are many), 30 steps.

| np | weight `max/mean` (what the balancer optimizes) | **box count** `max/mean` | measured `[rank_time]` |
|---|---|---|---|
| 2 | 1.015 | **1.130** | 1.087 |
| 4 | 1.054 | **1.176** | 1.170 |
| 8 | 1.050 | **1.308** | 1.259 |

**Box count predicts measured time; cell weight does not.** The cell-based model sits flat near 1.05
while both box count and measured time climb together, and the gap widens with rank count — the
scale-invariance failure this document names as the real requirement. Level 1 is perfect in both
metrics (32 equal boxes); all the divergence is at level 2, where 104 boxes over 8 ranks give weight
1.050 and count 1.308.

The mechanism was already in our own data: a per-block advance costs ~1x a full monolithic step
**regardless of block size and does not amortize** (@ref amr_block_batching). If per-block cost is
essentially fixed, a rank's load is how many boxes it holds, not how many cells. Equal cells with
unequal counts reads as perfectly balanced and runs badly skewed.

**Two consequences.** First, **step 4's flat A/B now has an explanation**: per-level distribution
redistributed *cells*, and cells were never the binding cost — it was balancing the wrong quantity,
not proving distribution irrelevant. Second, **step 5's conclusion that "box supply binds, not the
owner mapping" does not survive**; it was inferred from box counts and model imbalance with no
measured counterpart, and the model it rested on does not track reality.

Caveats, since this redirects the effort: single run per rank count (np=4 measured 1.084 in one run
and 1.170 in another, so variance is real — the *pattern* held in both); `[rank_time]` brackets
compute regions, so some spread may be MPI wait rather than work; one case, one node, np<=8. The
direction and its growth with scale are actionable, the exact ratios are not.

## Measured outcome

hpcfund MI250X, amdflang OMP offload, 1 rank/GCD, 2047x1023, `amr_max_level=2`, `amr_subcycle=F`,
`amr_max_grid_size=128`, 8 stripes of which 3 are sharp (3 deep towers among 8 level-1 blocks),
20 steps, 3 reps, median s/step with the first 2 steps dropped. Arms are two worktrees:
`26e0d080` (co-location) vs `cfdd2847` (per-level). Harness: `amr-bench/sweep_ml.sh`.

| np | co-located | per-level | ratio | uniform control |
|---|---|---|---|---|
| 1 | 2.6403 | 2.6188 | 0.99x | 0.0648 |
| 2 | 1.5315 | 1.5372 | 1.00x | 0.0533 |
| 4 | 1.0393 | 1.0943 | 1.05x | 0.0486 |
| 8 | 0.9186 | 0.9128 | 0.99x | 0.0458 |

**Flat at every rank count** (rep scatter 2-6%). Two readings, and the doc should not pretend to
choose between them without the metric:

- *Consistent with design.* With 3 heavy towers and np<=8, co-location can already place each
  tower on its own rank. The granularity floor says per-level cannot help until ranks approach the
  box count per level, so np<=8 is the wrong regime to see it.
- *Not demonstrated.* Equally consistent with the assignment not changing at all. **A flat
  end-to-end time is not evidence either way** — that is what step 5 exists to resolve, and it is
  why step 5 now precedes further distribution work.

### The overhead finding

The uniform control solves the same base grid with no refinement, so it bounds AMR's *overhead*,
not its efficiency (the efficiency baseline would be a uniform grid at the finest resolution).
Against it, AMR costs ~20x at np=8. The refinement added accounts for roughly

    1 + 0.20*(4-1) + 0.075*(16-4) ~= 2.5x

(level 1 over ~20% of the domain at 4x cells; level 2 over ~7.5% at 16x). Compare at **np=1**,
where both arms are still compute-bound: 2.6403 / 0.0648 = **40.7x** for ~2.5x the nominal work, a
residual of **~16x**. Regrid is only 7-10% of runtime, so it is not regrid. That points at
per-block launch cost — @ref amr_block_batching — as the dominant term for "efficient AMR", which
per-level distribution does not address. Treat the arithmetic as an order-of-magnitude estimate:
the area fractions are nominal and a 16x-cell uniform grid would not cost exactly 16x on a GPU.

Do **not** read this ratio at np=8 (20x, apparent residual 8x). The uniform control degrades faster
than AMR there — see the regime warning below — so the gap narrows for a reason that has nothing to
do with AMR overhead.

### Regime warning: this sweep is too small above np=2

Parallel efficiency, from the same runs:

| arm | np=1 | np=2 | np=4 | np=8 |
|---|---|---|---|---|
| uniform | 100% | 61% | 33% | **18%** |
| co-located | 100% | 86% | 64% | 36% |
| per-level | 100% | 85% | 60% | 36% |

The **uniform** control — no AMR, no blocks, no balancing — reaches only 1.41x on 8 GCDs. At np=8
the case holds 262k cells/rank, far below what saturates an MI250X GCD, so per-step fixed costs
dominate and the run is latency-bound. That is the machine floor for this problem size: AMR cannot
scale better than the uniform grid on the same mesh.

**Consequence: the flat A/B above is not conclusive.** At the rank counts where per-level
distribution should help most, both arms are dominated by costs that distribution cannot affect,
so the measurement could not have detected the effect even if present. Any future distribution A/B
must keep ~1M+ cells/rank at its TOP rank count (np=2 here, at 1.05M cells/rank and 86%
efficiency, is the last trustworthy point) — for np=8 that means a 4096x2048 case. Always run the
uniform control and check its efficiency FIRST: if the control is not scaling, nothing measured
against it means anything.

**Implication for the thesis.** Per-level distribution removes a structural ceiling (proven: a
tower can now split across ranks and stay correct). It has not been shown to lift a ceiling that
was binding on any case measured so far. Both statements should survive into whatever comes next.

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
