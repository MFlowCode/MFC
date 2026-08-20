# AMR performance plan (2026-08-19 — post-diagnosis execution plan)

**Mission: drive the AMR infrastructure tax toward zero.** Physics (`rhs`, `coarse`, `rk`) is
untouchable; everything else is overhead to be removed. This version supersedes the 2026-08-18
rewrite (git history) now that the WHY is established — the findings live in
`amr_slowness_analysis.md` (causal model, five-reviewer panel) and `amr_tax_review.md` (measurement
audit); this document is only the work list, its gates, and its decision rules.

## Where we stand

- **Landed:** R1 (level-2 gather blocking-SEND -> ISEND pool): **-17%/-22% wall**, regrid -39%,
  goldens 76/76. Phase instrumentation: 44 brackets, call counts, per-rank output; zero cost when
  `rank_time_wrt` is off.
- **Measured taxes:** matched point 6.30x -> **5.21x** post-R1. Production point: **not a number** —
  4.02x in the 40-80 step window growing to 12.14x in 80-160 at constant mesh; the growth costs
  ~61% of a 160-step run. Long runs at `regrid_int=2` **OOM by regrid count** (40 regrids dies).
- **Store fix A' applied** (one line, growth policy): the OOM mechanism is diagnosed as a plateau
  overshot by doubling, not a leak - Phase 1. Verification in flight on the case that dies today.
- **The causal model** (`amr_slowness_analysis.md` sec. 3): regrid churn -> rank-local grow-only
  store ratchet -> (a) OOM at the 64->128 doubling, (b) VRAM pressure -> slow per-launch alloc path
  -> rank-local rhs divergence -> **convoy amplification** through the per-box blocking lattice.
  Link (b) is the leading hypothesis (E-H1), under adjudication now.
- **Two laws that gate every fix below:** (1) a per-box rendezvous is also a BARRIER — deferral
  without a downstream sync relocates cost (measured, twice); aggregation removes it. (2) Code-read
  attributions run ~2-for-11 here — every fix is gated on a bracket or counter, never on a reading.

## Phase 0 — finish the adjudication (hours, in flight)

| item | action | decides |
|---|---|---|
| 0.1 | Read RK_160 per-rank data (running) against pre-registered signatures: rhs bimodal at flat fine_work + same slow ranks = E-H1; rhs tracks block/L2 count = composition; rotating slow ranks = thermal | which Phase-1 premise holds |
| 0.2 | `map(alloc:)` revert A/B at 160 steps (2 lines, temporary) | E-H1 (VRAM-pressure alloc path) vs heap fragmentation |
| 0.3 | Land trip-wire (stderr: `amr_loc_n`, `amr_st_cap`, live slots per rank per regrid) + `PH_RESTR` bracket on the post-stage restrict/reflux chain (~15 LOC) | makes the ratchet and the invisible 4-6% observable |
| 0.4 | Commit the verified instrumentation (reflux + per-rank brackets, goldens green) | keeps the tree honest |

## Phase 1 — the store fix: DIAGNOSED, A' APPLIED, verification in flight

### What the counters actually showed (2026-08-19, `logs/recon-0819_1351`)
```
rank 5  live 29  loc_n 48  freed 19  stack_in 0  stack_out 19
rank 5  live 32  loc_n 89  freed 57  stack_in 0  stack_out 57
rank 5  live 29  loc_n 89  freed 60  stack_in 0  stack_out 60
rank 0  live 26  loc_n 31  freed  4  stack_in 1  stack_out 5
```
**RETRACTION: it is not an index leak.** `loc_n = live + stack_out` EXACTLY on every line - nothing is
lost. It is a **PLATEAU**: `stack_in` is 0 every time because the rebuild drains the whole recycle
stack then needs more, since frees happen AFTER allocs. `loc_n` settles at ~3x live (29 live + 60
parked); the floor with this ordering is ~2x, old and new being concurrently live during a rebuild.

**The OOM comes from DOUBLING overshooting the plateau**, not the plateau itself: `loc_n` tops at 89,
`newcap = max(2*oldcap, nloc)` jumps 64 -> **128** = 28.3 GB of store holding 6.4 GB of live data,
total 60.4 GB of 68.7.

**THE ASYMMETRY IS CHURN, NOT WORK:** rank 0 frees 4-5 slots/regrid, rank 5 frees **57-60** - 12x, at
identical `fine_work`. Which subdomain the feature migrates through is a deterministic geometric fact
and **nothing balances churn**. Rank 5 is also the 4x `rhs` straggler in both per-rank runs.

### FIX A' — applied, one line — FAILED ALONE, still in tree (see the attribution caveat above)
`newcap = max(oldcap + max(oldcap/4, 8), nloc)`. Trajectory 8,16,24,32,40,50,62,77,**96** (~21 GB,
total ~53 GB) instead of 16,32,64,**128**. Growth POLICY only - slot lifetime, index assignment and
the device-authoritative contract untouched, so it cannot produce a wrong answer, only a different
allocation trajectory. Checked: `newcap >= nloc` always (never undersized), `newcap > oldcap` always
(cannot stall).

**A' RESULT: FAILED (2026-08-19, `logs/fixA-0819_1456`).** The case still aborts, same
`HSA_STATUS_ERROR_OUT_OF_RESOURCES` on rank 5, never reaching the step loop. A' worked as designed -
cap 128 -> 120, lower ranks 64 -> 50-77 - but that reclaims only 1.8 GB. TWO errors in the projection:
(a) `loc_n` does NOT plateau at 89; that was measured over 20 regrids and over 40 it reaches **103**,
so growth policy alone cannot bound it; (b) the non-store footprint was estimated at 32.1 GB from one
historical figure, while the VRAM sampler measured **63.9 GB peak**, back-solving to ~35.6 GB. At cap
120: 35.6 + 26.5 = **62.1 GB** against an effective ceiling near 64 - it dies by a nose. A' is kept
(strictly better than doubling, safe, one line) but is NOT sufficient. **Fix B is therefore required,
not optional**, and is under test now.

**PASS BAR:** the 80-step / 40-regrid matched case - which aborts today with
`HSA_STATUS_ERROR_OUT_OF_RESOURCES` - must COMPLETE with zero OOM strings; then goldens 76/76.
Live evidence mid-run: ranks at `nloc` 56-63 sit at cap 62 where doubling would have forced 128; only
the two churniest reached 96; **no rank has requested 128**; and `nloc` 89 on rank 5 reproduces the
plateau measured independently in the reconcile run.

**If it passes**, three things unlock together: long runs at realistic regrid intervals, the
differenced matched-tax measurement that was never possible, and a genuine test of whether relieving
capacity pressure touches the rank-5 straggler (E-H1's real prediction, still untested).

### FIX B — APPLIED AND VERIFIED ON THE OOM CASE (2026-08-19)

`s_amr_compact_store(nlive)`: dense renumber of `amr_loc_of` + shrink realloc of all four store
arrays, called at the end of `s_amr_reconcile_slots`. Trigger `amr_st_cap > 3*nlive`, target
`2*nlive`.

**Result: the case that reliably OOMs today COMPLETED, rc=0, 814.570 s step-loop wall, and the
full suite is 76/76.** Compaction fires as designed (`loc_n 96 -> 29, cap 58`).

> **CORRECTION (2026-08-19): the OOM case is NOT the cap-96 case.** It is
> `regrid_int=2` / `amr_max_grid_size=64` / `weno_order=1` / `riemann_solver=5` / 80 steps — a
> deliberate high-CHURN stress case (40 regrids) with a cheap scheme so store churn dominates the
> run. "cap 96 OOMs" is a **separate, independent** finding from a different experiment. An earlier
> revision of this document conflated the two; whether Fix B unblocks cap 96 is **untested**.

> **ATTRIBUTION CAVEAT — the verified configuration is A' + B TOGETHER, not B alone.** Fix A'
> (1.25x growth, `m_amr.fpp` `s_amr_st_reserve`) was left in the tree when B was built and tested.
> What is established: A' *alone* does not prevent the OOM; A'+B does. **B alone is untested.**
> Do not claim "compaction fixes the OOM" without either running B-alone or keeping both. The two
> are plausibly complementary — A' limits overshoot on the way up, B reclaims on the way down —
> but A' also makes reallocs *more frequent*, and each one is a multi-GB host round trip, so A'
> may be a net time cost once B exists. Resolve this by measurement, not reasoning.

Peak local indices still reach `nloc 103 / cap 112` on the hot rank, so the plateau is *reduced,
not removed* — expected, since the trigger is deliberately hysteretic (see next block for why, and
for the fix that removes the hysteresis).

### The hysteresis is an artifact of a host round trip — external review (2026-08-19)

Reviewers on AMReX, Parthenon/Athena++ and Chombo/SAMRAI converge on the same two points, and both
are actionable:

1. **Our compaction and growth both stage the whole store through the host**
   (`GPU_UPDATE(host)` -> `tmp = store` -> realloc -> `GPU_UPDATE(device)`). Multi-GB over PCIe, so
   we cannot afford to run it often — hence the 3x trigger and the 2x target. AMReX's
   `RemakeLevel` fills the new `MultiFab` **device-side** from the still-live old one. Making our
   remap a device-to-device gather lets compaction run unconditionally at `newcap = nlive`, taking
   steady state from 2-3x live to **1x**, with a transient peak of 2x *live* rather than today's
   2x *cap* — cheaper than the current peak even before the steady-state win.
2. **Neither AMReX nor Parthenon allocates a local index at all.** AMReX's `localindex` is a binary
   search over the sorted owned-box list; Parthenon's `lid = n - nbs` is recomputed contiguous every
   regrid. The index space *is* the live set by construction, so it cannot ratchet. This is Option 1
   (stable slot identity), and it deletes `amr_loc_free`/`amr_loc_nfree`/`amr_loc_n` as concepts
   rather than mitigating them.

Do (1) then (2): (1) is the larger memory win and is local to two routines; (2) removes the bug
class. Together they reproduce AMReX's invariant — *storage size == live working set, every regrid,
with no memory of the past* — while keeping the one contiguous device array the offload result
requires.

Two caveats banked from the same review, against over-claiming the batching payoff:

- Parthenon's headline 82x packing win is **launch-latency-bound buffer fill** (kernels copying ~8
  numbers). MFC's blocks are ~2.3M cells; a per-block kernel runs for milliseconds and cannot be
  launch-latency-bound. Our overhead is a *different* mechanism (offload descriptor mapping), so
  that number does not transfer.
- Parthenon-VIBE, on a code that already has full pack fusion, still measures **4.4% GPU-busy** on
  a 3-level AMR case — worse than ours. Kernel fusion removes one term; it does not fix host-side
  mesh bookkeeping.

### FIX B — original deferred plan (superseded by the block above)
Compact post-reconcile so `loc_n` returns to `live` (peak ~2x not ~3x): recovers ~15.5 GB against
A''s ~7.3 GB. Safe in principle (post-reconcile all readers are done) but it MOVES LIVE DATA - the
silent-corruption class. Price only after A' is measured; it may not be needed.

### REFUTED — kept as a warning
"Free stale slots BEFORE the rebuild's alloc loop" (~20 LOC) was this phase's original fix. Its
premise - stale slots have no readers by rebuild time - is FALSE: the overlap carry-forward reads
`amr_stor_st(..., amr_loc_of(kks))` for OLD blocks (m_amr_regrid.fpp ~1422) while the stash is
WRITTEN at those same old local indices (~1221), and `s_amr_free_slot` zeroes `amr_loc_of` and
recycles the index. Silent wrong answers. Patch parked at `amr-bench/phase1_patch.py` (raises on
execution). **Caught by tracing runtime reads, not by the anchor check - the eleventh code-read
attribution refuted, and the first caught before burning a build.**

## THE RUN QUEUE (2026-08-19, one node / 8 GCDs, serial)

Ordering is the scarce resource: one node means one timing arm at a time, and the overlap policy
forbids building during a timing arm. So the queue is ordered by **how much of the remaining plan
each run invalidates**, not by how promising each idea is.

Every entry names the decision it changes. A run that cannot change a decision does not belong here.

| # | run | cost | decides |
|---|---|---|---|
| **R1** | **A5 discriminator** (RUNNING) | build + 80 + 160 steps, ~40 min | **Whether the balancer track exists at all** |
| **R2** | **Cap sweep 64 / 96 / 128** | 3 runs, ~45 min | Whether Fix B unlocked a faster operating point |
| R3a | A7 SFC-cut hysteresis A/B | 40-60 LOC + build + 2 runs | *gated on R1* |
| R3b | A1 per-block cost constant | ~10 LOC + build + 2 runs | *gated on R1* |
| R4 | Fix B alone vs A'+B | build + 1 run, ~30 min | One line (`s_amr_st_reserve` growth factor) |
| R5 | Device-side store remap | days | Independent; steady-state store 2-3x live -> 1x |

### R1 — the gate (running)

Reads per-rank `gather` against per-rank `rhs` at 80 and 160 steps, so the discriminator is read on
the **onset difference** rather than one endpoint (at 80 steps there is no straggler at all).

- `gather` skews with `rhs` -> **L0-sourcing**: the sick rank sources its coarse gathers from a
  *static* L0 subdomain. No fine-block balancer can move that work. **A1/A7 are dead**; go to A6.
- `gather` flat while `rhs` is 4x -> L0-sourcing is out; the churn correlation survives; A7/A1 live.

The binary also prints per-rank `[amr-recon]` churn and `[amr-store]` peaks, which supply the
reviewer's *blocks-arrived/departed* counter for free. The one thing still missing is
*coarse-gather bytes sourced*; add it only if R1 comes back ambiguous.

**Why this must be read carefully:** the apparent refutation of L0-sourcing (`rhs` 2.91 vs `gather`
1.127) compares numbers from **two different cases** — 2.91 is cap-64 production at 160 steps, 1.127
is the cap-32 matched case over a 10-warm-step window. R1 measures both on one case in one run,
which is the only reason its answer is worth anything.

### R1 + R1b RESULT — SETTLED: the straggler is a STORE-CAPACITY artifact

Same node (k004-004), same case, **byte-identical `fine_work` sequence**, the only difference being
Fix A' + Fix B:

| | wall | ns/fine-cell-step | rhs max/mean | rank-5 rhs | max cap |
|---|---|---|---|---|---|
| pre-fix (reverted) | 1693.195 s | 56.80 | **2.694** | 923.6 s | **128** |
| post-fix (A'+B) | **712.032 s** | **23.89** | **1.069** | 221.3 s | **77** |

**2.378x wall at identical work.** The straggler returns when the fix is removed and vanishes when
it is restored, on the same hardware. Mechanism: caps of 128 drove two GCDs to **63.8 / 62.5 GiB of
64 (97-99%)**; ungoverned store growth exhausts VRAM and pushes onto a slow path. It is
state-dependent — the first 80 steps match across nodes to 0.6%, and all divergence is in 80->160.

**Consequences for this plan:**
- **L0-sourcing REFUTED.** `gather` never skewed to rank 5 (argmax was rank 6). The gate that was
  blocking A1/A7 is lifted — but see below, because the *reason* to do them has changed.
- **"Cost grows with sim time" is not intrinsic**: 3.68x pre-fix, **1.17x post-fix** (linear).
- **The balancer track is deprioritised on its merits, not by the gate.** Work was balanced to 1.4%
  throughout and the straggler had nothing to do with partitioning. A1/A7 remain cheap and
  defensible, but they are no longer aimed at a known 2.4x problem.
- **The device-side remap (R5) is promoted.** If capacity is what costs 2.378x, then getting steady
  state from 2-3x live down to 1x live is the highest-value remaining store item.

**Still unexplained, do not paper over:** ranks 4 and 5 both reach cap 128 pre-fix but only rank 5
is slow (294 s vs 924 s). Capacity is necessary, not sufficient. One run with a per-rank VRAM print
would close it.

### R1 — PRE-REGISTERED READ (written before the 160-step arm returned)

Recorded in advance because the failure mode of this campaign has been fitting a story to numbers
after seeing them. The 80-step baseline is in: `rhs` imbalance **1.089** (pathology absent), and
per-rank `gather` **56.86 on rank 5 against a 52.14 mean** — mildly above, but the argmax is rank 6.

**The direction of the gather skew is the discriminator, not its magnitude.** `gather` contains MPI,
so a slow rank makes *other* ranks wait inside it. The two hypotheses therefore predict opposite
signs, which is a far stronger test than "is gather skewed":

| | rank 5's `gather` | other ranks' `gather` |
|---|---|---|
| **L0-sourcing** (rank 5 does the work) | **HIGH** | low |
| **Consequence of `rhs` skew** (others wait on rank 5) | **LOW** | high |

The `reflux` column already shows what a pure wait-sink looks like: rank 5 is the **lowest** (20.87
vs 31.83 mean, 0.656x) precisely because it is the one computing while everyone else waits. So if
`gather` were only absorbing rank 5's lateness, rank 5 would be *lowest* there too. **At 80 steps it
is not — it is above the mean.** That is a genuine, if weak, point against pure wait-artifact.

Decision rule, fixed now:

- **L0-sourcing CONFIRMED** — at 160 steps rank 5's `gather`/mean rises materially from 1.09 *and*
  rank 5 is the argmax. Balancer track (A1/A7) is dead; go to A6.
- **L0-sourcing REFUTED** — `gather` imbalance stays near 1.1-1.2 or its argmax is another rank,
  while `rhs` imbalance reaches ~2.9. Note this outcome does **not** confirm the churn hypothesis;
  the rank-4/rank-5 control pair below already damages that independently.
- **AMBIGUOUS** — `gather` rises to roughly 1.3-1.5, or rank 5 goes *below* mean (wait-artifact
  contaminating the signal). Then the *coarse-gather bytes sourced* counter is required, and no
  balancer work starts until it is run.

**The control pair, which stands regardless of the above.** At 80 steps ranks 4 and 5 carry
near-identical churn (91 vs 88) and store high-water (nloc 65 vs 64) but `rhs` of 0.998x vs 1.089x
mean. **Churn does not determine `rhs` time.** This reproduces, within a single run, the anomaly
previously noted across two ("rank 4 also reaches cap 128 and is NOT slow") and materially weakens
the churn/store-plateau mechanism as a *cause* of the straggler.

**Do not quote the 80-step correlations** (churn-rhs +0.54, gather-rhs +0.47, churn-gather +0.77).
With n=8 and a 1.089x signal they are not distinguishable from noise; they are recorded only as a
baseline for the 80 -> 160 change.

### R4 RESULT — both store fixes are load-bearing; they fix TWO different costs

| config | wall | ns/fine-cell-step | rhs imbalance | max cap |
|---|---|---|---|---|
| neither | 1693.2 s | 56.80 | **2.694** | 128 |
| **B alone** | 941.8 s | 31.59 | **1.106** | **128** |
| A'+B (2 runs) | **729.7 s** | 24.5 | 1.06 | **77** |

**B alone 1.80x | A' on top 1.29x | combined 2.32x.** All four at byte-identical `fine_work`.
Dropping A' costs **29.1%**, far outside the 4.96% noise floor. **Commit both.**

The mechanism was predicted before the run and confirmed: compaction fires at `cap > 3*nlive`
(trigger ~87), and 2x growth doubles 64 -> 128 in one step straight past it — so B alone overshoots
to cap 128/116 and cleans up afterward (5 compaction events), while A' prevents the overshoot.

**The non-obvious result: these are two separate costs.** B alone already removes the straggler
(2.694 -> 1.106) yet remains 29% slow. So (1) the capacity *cliff* causes the straggler and
compaction fixes it, and (2) *carrying* excess capacity costs a further 1.29x with no straggler
involved. "Just add compaction" would have left a third of the win unclaimed.

### R2 RESULT — cap 64 is an interior optimum; stop sweeping upward

| cap | wall | ns/fine-cell-step | work imbalance | outcome |
|---|---|---|---|---|
| 64 | 747.333 s | **25.07** | 1.027 | **best** |
| 96 | 888.174 s | 29.58 (+18%) | 1.140 (max 1.237) | runs |
| 128 | - | - | - | **device OOM** |

Work matched to +0.7% between 64 and 96, so the 18% is overhead, not less useful work. With the
prior cap 32 -> 64 result (2.26x), the curve is worse-below / best-at-64 / worse-above /
infeasible-beyond.

**Three results, one of them a refutation of my own premise:**

1. **Fix B generalizes.** cap 96 completes (rc=0) where it previously OOMed twice. Fix B was
   verified on the `regrid_int=2` churn case; this is independent confirmation.
2. **THE NOISE FLOOR IS 4.96%.** The cap-64 arm re-ran the exact configuration measured at
   712.032 s and returned 747.333 s. This is the first run-to-run spread for this case. It licenses
   the 2.378x store result (far outside it) and is now the yardstick: **anything under ~5% is not a
   result.**
3. **cap 128's OOM is a DIFFERENT failure from the pre-Fix-B ones.** `HSA_STATUS_ERROR_OUT_OF_
   RESOURCES` at **nloc 16, cap 16** - immediately, at a tiny slot count. A cap-128 level-1 block
   spans 256^3 fine cells, so a few slots plus the stash exhaust 64 GiB. This is the raw per-block
   working set, **not** the capacity ratchet, and Fix B cannot help. Do not conflate the two.

**My premise for R2 was wrong and the direction is closed.** I ranked it "most likely to produce an
outright win" on WarpX's measured optimum of ~9 boxes/GPU. It does not transfer: at cap 96 we sit
near their optimum and are measurably worse, because larger blocks degrade partitioning granularity
(1.027 -> 1.140) faster than they save per-block overhead. Different code, different per-block cost
structure.

**Caveat on box counts:** the parser's `est blocks = fine_work/cap^3` gives 711 at cap 64 against
~224 recorded elsewhere. It scales correctly across arms (711 -> 212) so it is internally
consistent, but its absolute value disagrees - `fine_work` may not be plain fine-cell count.
Compare caps on ns/fine-cell-step, which needs no box count.

### R2 — original premise (superseded by the result above)

**Premise, stated correctly.** `amr_max_grid_size = 96` OOMed in an earlier, separate experiment.
Fix B was verified on a *different* case (the `regrid_int=2` churn stress case) — so **whether Fix B
unblocks cap 96 is exactly what R2 tests**, not something it assumes. If cap 96 still OOMs, that is
itself the result, and it says compaction's 3x/2x hysteresis is too loose for that configuration.

External evidence says block count matters more than any balancer change:

- WarpX measured 150 boxes/GPU giving *better* load-balance efficiency but *worse* walltime than 9
  boxes/GPU (1040 s vs 896 s) — "the overhead associated with a greater number of boxes outweighs
  the performance benefit of improved load balance efficiency". Their optimum was **9 boxes/GPU**.
- Three independent sources converge on **boxes per rank >= 3-4** as the parameter that actually
  governs balance quality in Chombo/BoxLib-family codes.
- **We currently run ~224 blocks over 8 ranks = ~28 boxes/GCD** — far above WarpX's optimum and far
  above the 3-4 floor. Every per-block overhead we have measured is multiplied by that count.
- Our own prior A/B: cap 32 -> 64 gave 4.9x fewer boxes, lower memory, and **2.26x less wall**.

So the direction is established and the next step up was blocked only by the OOM that Fix B just
removed. Sweep 64 / 96 / 128 at the production point, matched on everything else, and report
wall + per-rank phases + box count.

**Caveat to apply when reading it:** larger blocks refine more volume than the tags require, so
arithmetic work rises even as overhead falls. Report **wall AND fine-cell count**, and compare
ns/cell, not wall alone — otherwise a cap that wins on wall by doing less useful work will look
like a free win. This is the same denominator error that produced the retracted 1.57x tax.

### Gating logic, stated once

R3a/R3b are **not** queued behind R1 as a formality — if R1 says L0-sourcing, the cost model and the
cut hysteresis both operate on a quantity that is not the bottleneck, and running them would produce
two clean null results that teach nothing. R4 and R5 are independent of R1 and can fill any gap.

## Migration cost is a SECOND quantity we do not model (external review, 2026-08-19)

Reviewers on p4est/Dendro/GAMER/Enzo-E, Zoltan/ParMETIS and Chombo/SAMRAI/Uintah converge on a
result we have no analogue for. Three independent groups, three decades, same order of magnitude:

| study | migration reduction | cut-quality cost |
|---|---|---|
| Schloegel/Karypis/Kumar (JPDC 1997), multilevel diffusion vs from-scratch | >10x total, >3.3x max-rank | within 5% |
| Walshaw (JOSTLE), 9 successive adaptive meshes | 20-100x (0.54-3.76% vs ~94% migrated) | within ~8% |
| Zheng/Bhatele/Meneses/Kale (IJHPCA 2011), Charm++ RefineLB vs GreedyLB | 16-30x, and it GROWS with scale | - |

And the blunt one, Bak et al. CCGrid 2018: **"GreedyLB doesn't work well even compared to execution
runs without load balancing."** Unconditional optimal remapping is a *net loss*.

**The concrete gap: ParMETIS carries TWO per-object arrays — `vwgt` (work, the balance constraint)
and `vsize` (redistribution cost). `s_amr_block_cost` has only the `vwgt` analogue.** Every scheme
in this literature needs both. For MFC the second is cheap and exact: a block's migration cost is
`sys_size x fine cells`, known from its region and level.

Two structural rules we violate, both nearly free to adopt:

- **t8code**: "only send a local tree to a process p if this tree is not already local on p."
  `s_amr_regrid_stash_migrate` has no such guard.
- **p4est/t8code, stated as a named design rule** — *Principle 2.1 (Complementarity principle)*:
  "A collective mesh operation shall either change the local element sizes within the existing
  partition boundary, or change the partition boundary and keep the elements the same, **but not
  both**." Its stated payoff is exactly our problem: it "simplifies the projection and transfer of
  simulation data" — interpolation runs processor-local because the partition has not moved yet,
  and the subsequent partition step transfers data without further changing the refinement pattern.
  **`s_amr_regrid_stash_migrate` does both at once**: it changes the box set, reassigns owners, and
  moves data in a single pass.

### Reconciling "partition is cheap" with "migration is expensive"

The two reviews appear to contradict each other. p4est says partitioning is "usually negligible in
terms of execution time, so we suggest to call it whenever load balance is desirable"; ParMETIS and
Charm++ say migration must be actively penalized or it costs more than it saves.

**Both are right, and the complementarity principle is why.** p4est's cheap operation is the SFC cut
over a forest of *elements* — pure topology, no field data. The *transfer* is a separate, explicitly
invoked step. So "the cut is cheap" and "the transfer is expensive" are statements about two
different operations that p4est keeps apart and MFC runs together. Fusing them means we pay the
expensive one every time we want the cheap one, and we cannot price them separately because no
instrumentation boundary separates them. Splitting the two is therefore a prerequisite for
*measuring* the migration cost, not merely a tidiness argument.

### COMPETING HYPOTHESIS for the sick rank — and the run that settles it

A reviewer on SAMRAI/Chombo/Uintah proposes a mechanism that would **invalidate the entire
load-balancer track**, so it must be tested before any balancer work:

> Cells balanced to 1.4% with one rank 4x slower is more consistent with *intrinsic* churn — the
> feature crosses that rank's **level-0 coarse subdomain**, so it sources all the per-block coarse
> gathers and migration packs — than with partitioner-induced ownership flipping. **The L0
> decomposition is static, so no fine-block ownership balancer can move that work.**

**Status: not yet discriminated, and our existing data does not settle it.** The competing
predictions differ in *which phase* carries the skew:

| hypothesis | predicts skew in |
|---|---|
| Reviewer's L0-sourcing | `gather` (and regrid migration) |
| Our churn/store correlation | `rhs` |

The phase budget shows `rhs` imbalance **2.91** against `gather` **1.127**, which looks like a clean
refutation of the reviewer — **but it is not**, and claiming so would repeat a mistake we have
already made twice. The 1.127 comes from a 10-warm-step window; the rank-5 pathology **develops
between steps 80 and 160** (at 80 steps `rhs` imbalance is 1.06 and there is no straggler). The two
numbers describe different regimes and cannot be compared.

**The settling run is cheap and the instrument already exists.** The per-rank gather reporting added
this session (`PR_ID = [PH_RHS, PH_REFLUX, PH_GATHER, PH_SEAM]` in `m_phase_timing.fpp`) prints
per-rank `gather` alongside per-rank `rhs`. A single 160-step run with the current binary answers it:
if rank 5's `gather` is flat while its `rhs` is 4x, the L0-sourcing mechanism is out; if `gather`
skews with it, the balancer track is dead and A6/B-class fixes are the ones that pay.

Add the two counters the reviewer suggests to make it unambiguous: per rank per regrid,
**blocks-arrived/blocks-departed** (partitioner churn) versus **coarse-gather bytes sourced**
(subdomain churn).

### Earlier hypothesis (SUPERSEDED by the above — keep for the reasoning, not the conclusion)

We have measured, at *identical* `fine_work`: rank 5 spends 976 s in `rhs` vs ~245 s for the other
seven, and separately that rank 0 frees 4-5 slots/regrid while rank 5 frees 57-60 — **12x the
churn**. Rank 5 is also the highest-`nloc` rank. The literature above says churn is a first-class
cost, which makes "the sick rank is the high-churn rank" a candidate mechanism rather than a
coincidence. It is currently only a correlation across two separately-measured quantities; the test
is to plot per-rank churn against per-rank `rhs` time in a single run. See
[[amr-one-sick-rank]] and [[amr-reflux-is-the-sink-not-the-source]].

### Corrections to earlier reading

- **Zoltan's `LB_APPROACH` applies only to `GRAPH`/`HYPERGRAPH`** — there is no migration-vs-quality
  dial for RCB/RIB/HSFC. The geometric equivalent is `RCB_REUSE` (default 0), which has no
  published payoff.
- The ITR dial's real useful range is **~2x**, not orders of magnitude: sweeping it over six
  decades buys 2.33x worse cut for 2.0x less redistribution, with the knee between 0.1 and 1.
- **Enzo-E has no measured load-balancing result anywhere**; its "cello" balancer is three lines of
  equal-count SFC cut. Its architectural claims should not be cited as measurements.
- **`GreedyRefine` semantic trap**: the old `+LBPercentMoves` was a *migration* budget; the current
  TreeLB `tolerance` is an *imbalance* budget (default 1). Cite the old parameter if citing it as
  prior art for a migration budget.

**Caveat to carry**: Walshaw also finds that when load changes *dramatically*, repartitioning from
scratch can win outright. This is a regime question, and our regrids are frequent and incremental —
which is the regime where diffusion wins — but that should be checked, not assumed.

What is genuinely absent from the literature is an end-to-end **wall-clock** attribution on a GPU
AMR code. MFC would be producing that datum, not reproducing it.

## Ranked items from the SAMRAI/Chombo/Uintah review (2026-08-19)

**The reframing:** none of SAMRAI, Chombo or Uintah puts migration cost in the partitioning
objective — all three balancers are pure `(boxes, weights) -> owners` with no `previous_owner`
input. `s_amr_sfc_cut` has the same signature shape. **MFC is not doing anything unusual; it is
missing the compensating mechanisms**, which live elsewhere: metadata-only movement, a gain
threshold, and anti-churn logic in the *regridder* rather than the balancer.

| # | item | LOC | confidence |
|---|---|---|---|
| A7 | **Hysteresis on the SFC cut.** Keep the previous `amr_fine_cut(:,lev)`, adopt the new one only if it improves max-rank load past a threshold. This is Uintah's `gainThreshold` (0.05) and WarpX's ratio threshold (**measured optimum 10%**; 5% churns too much, 15% too little). No correctness risk — any owner map is valid. | 30-50 | high |
| A1 | **Per-block constant term in `s_amr_block_cost`.** Uintah ships `patchCost = 16` cells; SAMRAI added `minimum_patch_load` for exactly the regime "uniformity in patch count matters more than in cell count". Our per-block constant is almost certainly large given the measured per-region overhead. | ~10 | high |
| A10 | **A dedicated migration phase bracket.** Uintah instruments regrid / LB / schedule / migrate separately, and that separation is what let them identify migration as the AMR-ICE scaling limiter. Makes the A5 question answerable in one run. | ~10 | high |
| A2 | **Fit the cost coefficients to measured per-block time** (Uintah `ModelLS`: regress time on `[cells, ib_cells, pc_iters, 1]`, Cholesky on the normal equations, smooth the *coefficients* with `alpha = 2/(min(iter,T)+1)`, T~20). Yields A1's constant automatically and is far more robust than using raw per-block times as weights. | 80-120 | medium |
| A6 | **Batch the per-block regrid data motion.** SAMRAI fills a whole new level with **one `RefineSchedule`**; Chombo with **one `copyTo` + one `Copier`**. We do one `s_amr_gather_coarse_patch` per block with its own allocs and P2P set. *"SAMRAI can afford a churn-blind partitioner precisely because moving load meant moving only metadata, not voluminous mesh data."* Converts churn from O(blocks) to O(1) plans. | weeks | high value |
| B3 | **Free each old slot as its overlap contribution is consumed**, rather than all of them after. Drops the peak toward `max(live_old, live_new)` instead of the sum. A reordering, not a redesign — but see sec. 9.2: the stash is read at old indices throughout the rebuild, so this needs per-block liveness, not a bulk move. | ~30 | medium |

**A4 — the negative result to heed before investing in measured cost.** Uintah's dissertation found
the measured filters *"perform worse than the model immediately following a regrid... and load
balancing will always occur immediately following a regrid."* WarpX's tuned optimum used the
**heuristic**, not the measured GPU clock, and CUPTI-based measurement made their run **2x slower**.
If we regrid every 2 steps, a fitted algorithmic model (A2) is defensible; direct per-block timing as
the weight is the risky one.

**A3 — keep cost keyed to SPACE, not block identity. We already do this and should not stop.**
Uintah deliberately profiles on a fixed lattice of regions, *"so the patch set can change without
needing to migrate forecasting data between the changing patch sets."* Our footprint-integral
formulation is right. Seed new blocks with the **domain-average** weight, as Uintah does.

**A relevant measured counterpoint to more/smaller blocks (WarpX):** 150 boxes/GPU gave *better*
load-balance efficiency but *worse* walltime than 9 boxes/GPU (1040 s vs 896 s), "because the
overhead associated with a greater number of boxes outweighs the performance benefit." Three
independent sources converge on **boxes per rank >= 3-4** as the parameter that actually governs
balance quality.

### Two source cautions

- **Chombo's design document does not describe the code that ships.** It says `LoadBalance` uses
  Kernighan-Lin knapsack; the KL knapsack has **zero callers outside a unit test**. The shipping path
  is single-pass greedy bin-fill, and the Morton sort its comment assumes is a separate function the
  application must call and the flagship examples do not. Do not cite Chombo as knapsack prior art.
- **Chombo's petascale numbers are not evidence about load balancing** — the 196K-core run
  *replaced* the balancer with lexicographic assignment. IPDPS'09 measured that load balance and
  communication volume were *"the lowest priority"* impediments.

## Phase 2 — the convoy pilot (days; ~40-60 LOC)

freg pre-post: post ALL participant IRECVs at stage start with per-block tags (buffers are already
per-block: `freg(D)%%lo(:,:,:,amr_cur)`), one WAITALL before the apply loop; owner hoists its
per-block WAITALL identically. Keep the subcycle twin (`m_time_steppers.fpp:793`) blocking — the R1
lesson: grep every generated/secondary call site, not the wrapper.

Gate on counts + per-call time: reflux ms/call at 160 steps should collapse toward its 80-step
value if the convoy mechanism is real. **This pilot is the go/no-go for the whole aggregation
family (Phase 3). If it is null, STOP and re-diagnose before spending Phase 3's ~500 LOC.**

## Phase 3 — per-level aggregation family (weeks; conditional on Phase 2)

In value/risk order, re-priced against post-Phase-1 budgets before each:

| item | today's share | LOC | risk |
|---|---|---|---|
| 3.1 L1 gather recv batching (`rb:wait`) | ~5% | 60-100 | medium |
| 3.2 post-stage restrict/freg chain (bracket first via 0.3) | ~4-6% + kills the step-to-step skew feedback | 100-150 | med-high |
| 3.3 fill-gather aggregation (per-block `amr_cg` patches, ~0.5 GB; cached overlap lists exist) | gather 13-17% | 150-300 | high |
| 3.4 seam halo aggregation | 2.5-4% | 60-100 | medium |

All four share Phase 2's pattern: post-all, per-block tags, one drain, cached lists, aggregation
WITH a sync point (never bare deferral).

## Phase 4 — regrid modernisation (weeks)

4.1 Device-side rebuild (prolong + stash on device): `rg:mig` 6.7% + `rb:ovl/slot` ~6%; ~120-200
LOC; HIGH risk — the documented history (device kernels silently not emitted in large routines,
the 3-change NaN). Piecewise, probe-first, one change per A/B. 4.2 capture-sweep precompute (~40
LOC, low-med) when rhs's non-busy share is re-measured post-Phase-1.

## Phase 5 — endgame (after the tax is stationary)

5.1 Cost-weighted balancer (AMReX `makeKnapSack` shape; infra exists in `m_load_weight`) — ONLY if
per-rank data still shows skew after Phases 1-3; under convoys it is absorbed, and after the ratchet
fix there may be nothing left to balance. 5.2 The honest AMReX head-to-head: matched cap 64, same
regrid_int, DIFFERENCED protocol both sides, `amr_buf` scaled with interval and matched `fine_work`
as an iterated gate — meaningful only once MFC's tax is a number again. 5.3 Multi-node scaling
(the per-box lattice anti-scales; Phase 3 is what changes that slope).

## Standing rules (the short list; canon in `amr_tax_review.md` sec. 8 and the memory index)

- Bracket before believing; counts before wall; one clock per ratio; never compose across runs.
- Every number carries its operating point AND its time window.
- From-scratch runs; VRAM settle gate between GPU runs; stderr for anything that must survive an
  abort; `gate_phases.py` before every build; grep generated twins for every touched call site.
- np=2 dynamic-regrid is the minimum correctness bar for anything touching slots, stash, or
  exchanges — np=1 and single-rank goldens prove nothing about this code.

## Expected trajectory (estimates, not promises)

Phase 1: production tax becomes a NUMBER again (~4x static at int=20) and long runs work.
Phases 2-3: the ~27% static per-box communication share compresses -> ~2.5-3x. Phase 4: regrid's
residual halves -> ~2-2.5x, at or below AMReX's 3.40x on its own protocol. Floor: the physics is
already at parity (per-cell arithmetic 0.96x uniform at cap 64) — everything above 1.0x is
infrastructure, and every line of it is now attributed.
