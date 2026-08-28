# AMR performance plan (2026-08-19 — post-diagnosis execution plan)

> **2026-08-20 RE-FOUNDING: read `amr_endstate.md` first.** The program is derived from the
> end-state architecture (four pillars, weak-scaling invariants W1-W8), not from phase shares at
> the matched point. This document is the evidence ledger and detailed work list; where its
> sequencing conflicts with the endstate ladder (notably: the batched advance is REINSTATED as
> Phase 2, the "kills batching-the-advance" reading was an operating-point artifact), the endstate
> document wins.
>
> ## CURRENT STATE (2026-08-27) — read this before picking work
>
> **The live increment plan is `amr_plan_based_exchange.md`** (~3100 LOC over I0-I8, priced, with
> per-increment gates). Work FROM it; do not re-derive scope. Landed: I0, I1a, I1b, I2a, I3, I4a,
> I4b, I5. **Outstanding: I2b, I5b (~250), I6 (~200), I7 (~600), I8 (unpriced).**
>
> **W4 (per-cell global collectives) — half done, and RESTATED in ledger 40.** S3.1 deleted the
> level-1 tag ALLGATHERV (`ntag_bytes` 185 MB -> 0 per rank per regrid, box set bit-identical). The
> residual wall is **not volume**: differencing per regrid shows the level-1 reducing path is FLAT in P
> in both per-rank bytes (1.39 M) and collective count (16,383), because that tree is cap-limited. It is
> **~16,383 global `MPI_ALLREDUCE`s per regrid, of which 99.93% are on nodes lying wholly inside one
> rank** and therefore buy nothing but a full-machine sync. The `rbytes` growth of 1.67x/1.83x per
> doubling that earlier motivated "S3.1 is still O(P)" lives entirely in the **level-2 forest**, which
> does zero collectives today. **W4 is NOT met**; the thing to remove is the collective COUNT. Remaining:
> **S3.2's design contract is written below ("S3.2 DESIGN CONTRACT") — read it before writing any
> S3.2 code; it records that depth-fusion ALONE is still O(P) per rank.** Sequencing matters and is
> NOT free: the level-2 forest runs `reduce = .false.` on a replicated list and performs **zero
> collectives today**, so converting it first (S3.3) would ADD ~22-28k collectives per regrid with no
> mechanism to scope them. **S3.2's mechanism must therefore be proven on the level-1 path, where the
> collectives already exist, and S3.3 then moves the forest onto it.**
>
> **W5 is the second wall, at ~28k ranks**, and it is BIGGER than one increment: 19 of 41 AMR p2p
> call sites still tag per box, F1 is only partially converted (I2b unlanded), and the subcycle
> sites are an explicit deferral to I8. See `amr_plan_based_exchange.md`.
>
> **Adopted ordering (ledger 41, with its premise CORRECTED by the Arm B verdict — BOTH the forest and
> level 1 are O(P); S3.3 went first on absolute size and readiness, not on a difference in exponent):** B0 -> **B1
> (canonicalise the merge)** -> **B0b (`amr_blocking_factor` default 1 -> 4, so the box cap stops being
> load-bearing)** -> **S3.3 (forest ownership split — deletes `gwin_bytes`, 1.44 GB/rank at np32 and
> doubling with P; also retires W1's pass-2 loops)** -> S3.2b (level-1 depth batching; a latency
> constant, 16,383 flat collectives per regrid) -> W5 (I2b/I7/I8) ->
> W3 (F5 aggregation, measured live at **71.7%% of all messages** and the smallest at ~0.8 MB/msg) ->
> W6 (chunked device-side store growth) -> W7 (needs a cost model chosen first). **W2 (batched
> advance) is FLAT IN P and is therefore an efficiency/parity item, NOT on the scaling critical
> path** — a scope decision, since it is also the largest and least-bounded item.
>
> **Working practice (measured 2026-08-27):** `amr-bench/fcheck.sh` structurally checks one .fpp in
> ~4 s versus ~18 min for a GPU build and catches what `ffmt`/`lint_source` pass; a gfortran CPU
> build in an isolated worktree is ~13 s incremental; build in one worktree and measure in another.

## 2026-08-28 (41) — np32 closes the sweep, and it REORDERS the program: the forest, not level 1, is the O(P) term

Three points, per rank, weak-scaled (400^3 -> 800x400x400 -> 800x800x400):

| | np8 | np16 | np32 | growth |
|---|---|---|---|---|
| **`gwin_bytes`** (level>=2 window `ALLGATHERV`) | 360 MB | 719 MB | **1,440 MB** | **1.997x, 2.001x = O(P)** |
| level-1 global `ALLREDUCE`s / regrid | 16,383 | 16,383 | 16,383 | flat |
| level-1 shared nodes / regrid | 7 | 11 | 15 | **+4 per doubling = O(log P)** |
| level-1 shared bytes / regrid | 19 KB | 35 KB | 58 KB | 1.865x then 1.645x, i.e. ~O(log P) |
| `shr_maxdep` (exchange rounds) | 3 | 5 | 7 | +2 per doubling = 2 log2 P - 3 |
| forest nodes that are rank-local | 99.79% | 99.79% | 99.77% | **constant** |

**`gwin_bytes` is an O(P) term in W4** (and, when this was written, believed to be the ONLY one --
**CORRECTED by the Arm B verdict below: level 1 is O(P) too; the 'flat' reading was a pinned cap**).
It is It is 1.44 GB per rank per run at np32 (~144 MB per
rank per regrid) and it doubles exactly with P. Extrapolated to 1e5 ranks that is ~450 GB per rank per
regrid. It is also LARGER than the 185 MB `ntag_bytes` gather S3.1 deleted -- S3.1 removed the smaller
of the two per-cell gathers.

**Level 1 is not the exascale wall.** Its collective COUNT is flat (cap-fixed at 16,383) and everything
about its shared set grows logarithmically. It is a latency constant worth removing, not an exponent.

**I made the two-point error again, and the third point caught it again.** Earlier this session I read
shared nodes as "1.57x per doubling ~ P^0.65, the rank-boundary surface" from np8 and np16 alone, and
extrapolated ~3,000 shared nodes per regrid at 1e5 ranks. With np32 the sequence is 7, 11, 15 --
**linear in log2 P**, giving ~63 at 1e5 ranks. Same failure mode as the retracted `P^0.73`: a ratio
between two points is not an exponent. The rule now is that no growth exponent enters this ledger on
fewer than three points.

**`shr_maxdep` came back 7** -- the value ledger 39 pre-registered as "the design needs rethinking".
It does not, and ledger 40 demoted that statistic on mechanism grounds BEFORE this number existed
(written while the np32 job was still running with zero scope lines in its log). The demotion stands
on its own reasoning: 7 = 2 log2 P - 3 gives ~31 exchange ROUNDS at 1e5 ranks, and rounds are not
volume.


### CORRECTION 2026-08-28 — "level-1 collectives are FLAT in P" is a CAP ARTIFACT, not a property

Prompted by "are you sure you're running the right input decks?", the decks were audited. **They are
correct**, and the audit is worth recording because it validates every number above:

- IC `hcid 306` is `1 + 4*exp(-(cos(pi x)^2 + cos(pi y)^2 + cos(pi z)^2)/0.15)` -- **periodic with period 1**,
  so it is a BLOB LATTICE: 8 blobs at 2x2x2 (np8), 16 at 4x2x2 (np16), 32 at 4x4x2 (np32). The feature
  count scales exactly with P.
- Resolution is constant at 200 cells per unit length in all three rungs, and the final box count scales
  with P: **594 / 1206 / 2405 = 2.03x, 1.99x.**
- Restart decks match their rungs exactly: `lustre_0.dat` is 6 vars x cells x 8 B (3.07 / 6.14 / 12.29 GB),
  and `x/y/z_cb` hold m+2 / n+2 / p+2 doubles. `simulation.inp` as-run confirms m/n/p, `t_step_stop = 200`,
  `amr_regrid_int = 20` identical, and `amr_blocking_factor` 1 in all three.

**But `amr_max_blocks = 8192` is held FIXED across the rungs, and the tree SATURATES it in every one:**
`[amr-tree] lmax` is exactly 8192, the `capped` warning fires on 10 of 10 regrids, and the node count is
**16,383 = 2 x 8192 - 1 -- exactly a full binary tree pinned by the cap.**

**So level-1's collective count could not have grown with P in this experiment no matter what the
algorithm does.** "Flat in P" above is measuring the cap, not the clusterer. The final box count DOES
scale with P (594 -> 2405), so the natural uncapped tree is larger at np32 and is being clipped back to
8192 leaves in every rung.

**What this does and does not change.**
- **`gwin_bytes` O(P) STANDS.** It is driven by tagged cells inside parent windows, which scale with the
  blob lattice; the cap bounds BOXES, not tagged cells. Its 1.997x / 2.001x is genuine.
- **The claim that level 1 is asymptotically harmless does NOT stand.** After B0b makes the recursion
  self-terminating, level-1 node count should track the physics and grow with P -- which would make the
  level-1 collectives O(P) too and weaken the argument for putting S3.3 ahead of S3.2b.
- Ledger 41's "+4 per doubling" for SHARED nodes is a ratio between rungs whose trees are the same size
  by construction, so it measures how the FIXED tree splits across more ranks. That is still the right
  quantity for the shared/local split, but it is not evidence about how the tree itself scales.

**FIRST RESULT IN (job 388669, np8, bf=4): B0b's premise CONFIRMED.** Zero `clustering capped`
warnings against 10-of-10 at bf=1, and `lmax` 6872-7438 rather than pinned at 8192 -- **the recursion
self-terminates at bf=4.** Level-1 nodes/regrid 13,825; boxes 576 vs bf=1's 594 (3%). Cross-check that
matters: **`gwin_bytes` per regrid is ~40 MB at BOTH bf=1 and bf=4**, so the forest gather is insensitive
to the blocking factor -- confirming it is driven by tagged cells, not by the tree or the cap, and that
its O(P) is not a cap artifact.

**Note the ordering decision is ROBUST to how this comes out.** Both S3.2b and S3.3 are required for W4;
the experiment decides the JUSTIFICATION for doing S3.3 first, not the action. S3.3 removes a measured
1.44 GB/rank term with a verified np>1 golden gate and a bounded design, so it is the right next code
change whether or not level 1 also turns out to be O(P). Do not block S3.3 on these jobs.

**Caveat on B0b at scale:** bf=4 gives `lmax` 7438 against a cap of 8192 at np8. If the tree scales with
P, a fixed `amr_max_blocks = 8192` will bind again at np16/np32 -- `amr_max_blocks` is a user-set case
parameter and is NOT auto-scaled. B0b makes the recursion self-terminating at the sizes the goldens and
this benchmark use; it does not remove the need to size the cap with the problem.

**ARM A CONFIRMS THE CONFOUND DIRECTLY (jobs 388669/388670).** At `bf=4` with the cap left at 8192:
np8 gives 13,825 level-1 nodes/regrid and **0** capped warnings; np16 gives **16,383 and 5 capped
warnings** -- it RE-SATURATED. So the level-1 tree genuinely grows with P (np8 sat below the cap, np16
wanted more and was clipped back to exactly 2*8192-1), and `bf=4` alone does NOT keep the cap out of the
way at larger P. **This is direct evidence that the earlier "flat in P" was the cap, not the clusterer.**
Arm B (`amr_max_blocks = 65536`, jobs 388698/388699/388700) is therefore the only unconfounded arm and is
the one to read for the growth rate.

**Consequence for B0b:** raising the blocking factor makes the recursion self-terminating at a GIVEN
problem size; it does not stop a fixed `amr_max_blocks` from binding as the problem grows. B0b remains
correct and necessary (it is what makes `force` inert so S3.2b/S3.3 can be bit-identical at the sizes we
test), but it is not a substitute for sizing the cap with the problem.

**Settled by experiment, not argument:** jobs 388669/388670/388671 re-run np8/16/32 with
`amr_blocking_factor = 4` (a RUNTIME namelist parameter -- no rebuild needed) on the pinned B1 binary
`simulation-b1-54da9155`. If level-1 node count then scales with P, S3.2b returns to the critical path
and the S3.3-first ordering must be re-argued on absolute size rather than on exponent.

### The reordering

The plan sequenced S3.2 (level 1) before S3.3 (forest) on the argument that "the forest performs zero
collectives today, so converting it first would ADD ~22-28k collectives with no mechanism to scope
them". That argument modelled S3.3 as *making the forest reduce-driven like level 1*. The measurement
says the forest does not need that:

- The forest is **99.8% rank-local, and the fraction is constant in P** (0.206%, 0.208%, 0.234%
  shared). The `ALLGATHERV` exists only so that EVERY rank can cluster EVERY parent's window --
  replicated global work whose input is gathered globally.
- So the bulk of S3.3 is **"each rank clusters only the windows of the parents it owns"**, which adds
  no collectives at all, deletes the gather for the 99.8% case, and deletes the replicated work with
  it. Only the ~0.2% of windows that straddle ranks need the shallow mechanism.
- Rank-invariance of the resulting box set then comes from **one `ALLGATHERV` of BOXES plus B1's
  canonical order** -- which is why B1 was worth landing first, and is already written.

**Caveat, stated so the reordering is not oversold: `gwin_bytes` is not today's wall-time bottleneck.**
At np32 it is ~144 MB per rank per regrid, roughly 15 ms of wire time -- small against a 42 s regrid.
Its phase (`rg:build`, 9.2% at np32) does not stand out, and regrid wall is not even monotone across
the sweep (400 / 278 / 420 s), which is the MPI-wait variance this ledger has documented before. The
case for S3.3 is **asymptotic**: 144 MB/rank/regrid doubling with every doubling of P reaches ~450 GB
per rank per regrid at 1e5 ranks, and no other measured term does that. Prioritising it is a bet on
the exascale goal, not a claim about the current operating point -- and it is the right bet only
because the goal is the end-state architecture rather than the matched-point ladder.

**Revised order: B1 -> B0b -> S3.3 (forest ownership split) -> S3.2b (level-1 depth batching).**
S3.3 moves first because it carries a MEASURED O(P) term and, at 144 MB per rank per regrid, the largest
absolute number on the board. S3.2b stays, demoted to a latency-constant fix. S3.3 also retires W1's
regrid pass-2 loops (parents x global tagged cells), so it closes two invariants.

## 2026-08-28 (40) — B1 landed, and W4's level-1 problem restated: it is COLLECTIVE COUNT, not volume

**Per-regrid decomposition of the level-1 reducing path** (differencing consecutive `[amr-scope-r]`
prints instead of reading the cumulative totals — the counters turn out to be near-constant regrid to
regrid, 7,7,7,... at np8, so a 60-step run measures this as well as a 200-step one):

| per regrid | np8 | np16 | ratio |
|---|---|---|---|
| global `MPI_ALLREDUCE`s | 16,383 | 16,383 | 1.00x |
| ...on nodes spanning >1 rank | 7 | 11 | 1.57x |
| local bytes | 1.39 M | 1.39 M | 1.00x |
| shared bytes | 18,960 | 35,360 | 1.87x |
| `shr_maxdep` | 3 | 5 | +2 |

**99.93% of the level-1 global collectives are on nodes that lie wholly inside one rank.** That rank
already holds every tag in the subtree, so the reduction changes nothing and costs a full-machine
synchronisation. The count and the per-rank volume are both EXACTLY FLAT in P under weak scaling
(16,376 -> 16,372 local nodes; 1.39 M -> 1.39 M bytes) because the level-1 tree is cap-limited.

**This corrects how W4 was being scored.** The wall is not bytes — per-rank level-1 volume is already
flat. It is ~16,383 global collectives per regrid, each synchronising all P ranks, so the cost grows
as (fixed count) x (latency ~ log P). At 1e5 ranks that is ~1.6 s per regrid of pure synchronisation.
The earlier `[amr-tree] rbytes` growth of 1.67x/1.83x per doubling that motivated "S3.1 is still O(P)"
is real but lives entirely in the **level-2 forest**, which runs `reduce = .false.` and performs zero
collectives today. Level 1 -- the only path that communicates -- is flat.

**Consequence for the S3.2b gate.** The gate recorded below as "per-rank reduced bytes FLAT across
np8/16/32" measures a quantity that is ALREADY flat and therefore cannot discriminate. **The gate is
now: level-1 global collectives per regrid drop from 16,383 to the shared-node count (7/11/...), and
those residual exchanges are sparse (only the overlapping ranks) rather than over MPI_COMM_WORLD.**
Box-set identity and `[amr-tree] nodes` unchanged still stand as correctness gates.

**`shr_maxdep` is DEMOTED, and this is a deliberate revision of a pre-registered rule.** Ledger 39
pre-registered "np32 = 6 -> the log2 model holds and S3.2 proceeds; 7 -> the design needs rethinking".
That rule treated `shr_maxdep` as a proxy for how much of the tree stays shared. It is not one: it is
a MAXIMUM over a geometry-dependent quantity (a box straddling a rank boundary stays shared however
small it gets), and in the sparse per-depth exchange it sets only the number of latency ROUNDS, since
all shared nodes at one depth exchange concurrently. At 1e5 ranks, `log2 P` rounds is ~17 and
`2 log2 P` is ~34; both are negligible against 16,383 global collectives. **So np32 = 7 would not
collapse the design.** The revision is grounded in mechanism read out of the code (`s_amr_find_split_sig`
scans a node's whole signature, so a participating rank needs the full 1D signature -- a box EXTENT,
not a volume -- and depth-batching makes rounds = maxdep), not in the number being inconvenient; and
it makes np32 LESS decisive rather than more. The statistic that does decide is per-rank
participation, which ledger 40 instruments below.

**S3.2a-2 (instrumentation, golden-neutral).** `amr_cl_me_nodes_r` / `amr_cl_me_rb_r` count the shared
nodes whose box actually reaches into THIS rank's subdomain -- what the sparse exchange would cost a
rank, as opposed to `shr_*_r` which prices the allreduce form where every shared node lands on every
rank. Reduced with MAX across ranks (rank 0 owns a domain corner and understates it) and printed as
`[amr-scope-me]`. **Pre-registered reading:** the design holds if `me_rb_max` grows at roughly the
domain-extent rate (~1.26x/doubling = P^(1/3) in 3D) while `shr_rb_all` grows ~1.87x, i.e. the ratio
me/all FALLS with P; it is in trouble if `me_rb_max` tracks `shr_rb_all`, meaning ranks participate in
nearly every shared node and sparsity buys little.

**What the golden suite can and cannot say about B1.** The 69 AMR goldens ran clean through B1 with no
box-set movement. That is evidence B1 BREAKS nothing; it is NOT evidence the property is exercised. The
merge's order-dependence needs ambiguous chains of too-close boxes, which need many boxes: the test cases
carry a handful, the benchmark carries 594-2405. So the suite cannot distinguish B1 from a no-op, and the
evidence that B1 is NEEDED remains the code reading of `nacc`'s three uses, not a test. Do not let a green
suite be reported as validating the canonicalisation.

**B1 (canonicalise the merge) implemented.** `s_amr_cluster` now sorts the accepted boxes by Morton of
`lo` (a stable insertion sort; accepted boxes are disjoint so their `lo` corners are distinct) before
the min-separation pass, and the pass removes a fused box by shifting down rather than swapping with
the last entry. The merge therefore depends only on the box SET, not on the traversal order that
produced it -- the precondition for S3.2 completing local subtrees in parallel. +34/-4 lines. **Moves
the box set; goldens regenerated in the same commit and the diff explained.**

## 2026-08-27 (39) — B0 + S3.2a: GOLDEN GATE CLOSED, and the np32 question that decides S3.2

`320fe106` was committed with the caveat that the 69-test AMR suite had not finished before the node
allocation expired. **It finished: 69 passed, 0 failed.** B0 (at the `bf = 1` default) and the S3.2a
scope instrumentation are therefore gated by the full suite, not only by bit-identity. The caveat
still stands for `bf > 1`, which no test exercises — that remains evidence-backed by the sweep only.

**B0 priced (np8, per regrid):** `bf=1` 38,255 nodes / 594 boxes; `bf=8` 6,475 / 576; `bf=16`
1,466 / 576. **26x fewer tree nodes produce the SAME answer**, and the cap-saturation warning stops
at `bf >= 4`. But B0 is NOT a W4 lever: at np8 it cuts LOCAL (free) nodes 5x and leaves the SHARED
set untouched (7 -> 8.8 nodes, 18,966 -> 21,547 B), because it prunes leaves and leaves are already
rank-local. Its value is latency, plus making `force` inert so B1 and S3.2 can be gated.

**THE OPEN QUESTION, and this machine can answer it.** `shr_maxdep` went 3 (np8) -> 5 (np16), i.e.
**+2 per doubling, not the +1 that log2(P) predicts**. If that rate persists, at 1e5 ranks the shared
phase would be ~31 levels — essentially the whole tree — and S3.2's premise collapses. np32
(job 388355) discriminates: **6 means the log2 model holds and S3.2 proceeds as designed; 7 means +2
is real and the design needs rethinking BEFORE it is written.**

**Do NOT spend big-machine allocation yet.** The need for the scoping technique is already
established here with a mechanism, not an extrapolation: the forest term grows 1.97x/1.99x (exactly
P) while the cap-fixed level-1 term is constant, so S3.1 is asymptotically O(P). Frontier would
refine WHERE it breaks, not WHETHER. The right moment for big-machine time is AFTER S3.2+S3.3, to
demonstrate per-rank reduced bytes actually flat at 1-4k ranks — a short counters-only run. The one
exception: if np32 returns `shr_maxdep = 7`, a small probe (np64-np256, minutes each) is the cheapest
way to distinguish "+2 forever" from "+2 then flattening".

## S3.2 DESIGN CONTRACT (2026-08-27) — scoping the clustering reductions

Written after S3.1 measured asymptotically O(P). This is the next design contract; S3.3 then applies
the same mechanism to the level-2 forest, which is the term that actually scales.

### Depth-fusion alone FAILS — the arithmetic that kills the obvious design

Fusing all nodes at tree depth d into one allreduce gives ~32 collectives per regrid instead of
16k-127k, which fixes LATENCY and nothing else: **an allreduce delivers the whole buffer to every
rank**, so per-rank volume is the TOTAL signature volume at that depth, `~3 * L * N^(2/3)`. At fixed
per-rank work (`L ~ P^(1/3)`, `N ~ P`) that is **O(P) per rank — the slope S3.1 already has.**
Recorded because it is the design the earlier "~17 fused collectives per regrid" note implies, and it
would have failed its own gate only after being built.

### The property that makes scoping possible

A rank can hold tags only inside its own subdomain, so:

1. near the root a node spans many ranks, but there are FEW such nodes (at most 2^d at depth d);
2. **below roughly depth log2(P) a node's box fits inside ONE rank's subdomain**, that rank holds
   every tag in the subtree, and **the whole remaining subtree needs ZERO communication** — it
   recurses locally exactly as the serial code does today.

At np32 that is ~5 shared levels against a measured `ldepth` of 32: **most of the tree should never
touch MPI at all.**

### The mechanism: sparse per-depth exchange, no communicator per node

Within the shallow phase a rank needs the reduced signature only for nodes its subdomain overlaps —
its own path down the tree, O(1) per depth, not all 2^d. So: **one sparse exchange per shared depth**
(~log2 P per regrid), each rank contributing and receiving only its own nodes. Per-rank volume becomes
`O(log P * local extent)` — **sublinear in P, which is what W4 requires.** Implement with
`MPI_Neighbor_alltoallv` on a graph communicator built once per regrid, or point-to-point to the
overlapping ranks; the overlap inversion already exists (`s_amr_ranks_overlapping`, `amr_ovl_gather`).
**`MPI_Comm_split` is the wrong tool and is not needed anywhere** — that open question dissolves.

### Assembling the box list

No rank holds the whole set, so close with ONE `ALLGATHERV` of the resulting BOXES — per-box global
data, which the endstate explicitly permits. Determinism needs a canonical sort key (Morton of `lo`,
then level) so every rank builds an identically ordered list regardless of arrival order; `f_morton`
is injective below 2^21 per dimension.

### What this does NOT solve — stated up front

- **`force` reads GLOBAL `nacc + nwork`** against the block cap, which a rank finishing a local
  subtree cannot see. **B0 is the dependency:** if a minimum box size stops the bisection reaching the
  cap, `force` never fires. If it does not, the cap needs a per-subtree budget, which moves the box
  set and breaks bit-identity.
- **Bit-identity is NOT free.** The serial recursion visits nodes in a deterministic LIFO order and
  accept/force depend on `nacc`; local subtrees completing in parallel change that order. Whether the
  box SET is identical must be PROVEN, and it is the gate.

### PREREQUISITE FOUND 2026-08-27: the merge is ORDER-DEPENDENT, so S3.2 cannot be bit-identical

Reading every use of `nacc` in `s_amr_cluster` — the state that carries traversal order into the
result — finds three, and the third is fatal to the bit-identity plan:

1. `force = (nacc + nwork + 1 >= cap)` — order-dependent, but **B0 makes it inert** by keeping the
   bisection away from the cap.
2. `if (nacc < cap)` — likewise inert once the cap is not approached.
3. **the min-separation merge** — scans the accepted boxes IN INSERTION ORDER, fuses the FIRST
   too-close pair, removes it by swapping with the last element, and restarts. If A-B and B-C are
   both within `thr` the outcome depends on which is found first, and the swap-with-last removal
   permutes the list every iteration.

**So the merged box set is a function of the order boxes were accepted.** S3.2 completes local
subtrees in parallel, changing that order. Bit-identity is therefore not merely unproven — it is
**unachievable** while the merge consumes insertion order.

**B1, its own increment, BEFORE S3.2.** Canonicalise the merge input: sort the accepted boxes by a
deterministic key (Morton of `lo`, then extent) before the min-separation pass, and make removal
order-preserving instead of swap-with-last. The merge then depends only on the box SET.
**This moves the box set once and is NOT bit-identical**, so it lands alone with goldens regenerated
and the diff explained. Folding it into S3.2 would confound "did the scoping break something" with
"did canonicalisation move the boxes", and every S3 gate so far has relied on bit-identity being a
clean signal.

**Revised S3 order: B0 -> B1 -> S3.2a -> S3.2b -> S3.3.** B0 and B1 together remove every order
dependence from the clusterer, which is what makes S3.2's gate mean anything.

### PREREQUISITE FOUND 2026-08-28: at the DEFAULT `amr_blocking_factor = 1` the cap is load-bearing, so S3.2b cannot be bit-identical

Ledger 39 priced B0 and noted the cap-saturation warning "stops at `bf >= 4`", but left the DEFAULT at
`bf = 1`. Re-reading `s_amr_cluster` against the logs shows what that costs S3.2b:

- `force = (nacc + nwork + 1 >= cap)` accepts a node unsplit to protect the `amr_max_blocks` array
  bound. It reads **global** `nacc`, i.e. how many boxes the whole traversal has accepted so far.
- At `bf = 1` the recursion does not converge on its own -- it splits until the cap stops it. The
  np16 run emits `[amr] WARNING: tag clustering capped` on **10 of 10 regrids**, and `[amr-tree] lmax`
  is exactly 8192 = `amr_max_blocks`. So `force` is live at the default, on every regrid.
- S3.2b's deep phase is PRIVATE: a rank finishing its own subtree cannot see global `nacc`. Enforcing
  the cap per rank is a different box set; ignoring it overflows.

**So B0 is not optional and not merely a latency win: making the recursion self-terminating is a
PRECONDITION for S3.2b's bit-identity gate to be meetable at all.** This was scoped wrongly in ledger 39,
which treated B0 as "latency, plus making `force` inert" without noticing that the shipped default
leaves it un-inert.

**B0b, its own increment: change the `amr_blocking_factor` default from 1 to 4.** Ledger 39's sweep is
the evidence -- saturation stops at `bf >= 4`, and `bf = 8` produces 576 boxes against `bf = 1`'s 594
(a 0.7% change in the final box set) from 26x fewer tree nodes. **4 rather than 8 or 16** because the
69 AMR goldens run on small grids (128^2 and below), where a minimum child extent of 8-16 cells would
distort refinement; 4 is the smallest value the measurement supports.

**Moves the box set, so it lands with goldens regenerated -- as its own commit, after B1.** Two
regenerations rather than one bundled commit: B1 (Morton merge order) and B0b (minimum box size) move
the boxes for unrelated reasons, and a combined diff could not be explained.

**Revised S3 order: B0 -> B1 -> B0b -> S3.2a -> S3.2b -> S3.3.**

### S3.3 DESIGN CONTRACT (written 2026-08-28, from reading the level->=2 block at m_amr_regrid.fpp:1118-1300)

**What the code does today.** For each level >= 2, with `plo = 1, phi = nboxes` spanning EVERY box at the
previous level globally:

- **Pass 1** builds each parent's dense window `gwin`, tags it from old level-(lev-1) blocks, but
  contributes only from blocks this rank owns (`if (amr_owns_all(ob))`), and packs those cells as
  `(int8 linear index, parent kb)` pairs.
- **One `MPI_ALLGATHERV` per level** sends every rank every parent's pairs. **This is `gwin_bytes`:
  360 / 719 / 1440 MB per rank at np8/16/32, exactly O(P).**
- **Pass 2 runs on every rank for every parent**, and rebuilds each parent's window with
  `do i = 1, ntot_g; if (gkb(i) /= kb) cycle`. That is a full scan of the GLOBAL pair list once per
  parent: **O(parents x global tags) replicated work** -- ~2,400 parents x ~9 M pairs per regrid at np32.
  This is exactly the "regrid pass-2 loops (parents x global tagged cells)" that W1 names.

**The change: give each parent an OWNER and exchange point-to-point instead of gathering globally.**

1. Assign parent `kb` to a rank deterministically and identically on every rank -- round-robin
   `mod(kb - 1, num_procs)` both balances and needs no communication to agree on.
2. Pass 1 packs as today, but buckets each pair by `owner(kb)` instead of appending to one list.
3. Replace the `ALLGATHERV` with **`MPI_Alltoallv`**. Per-rank send volume becomes O(this rank's own
   tagged cells) and receive volume O(the tags of the parents it was assigned) -- both O(local), which
   is exactly what W4 asks for. The measurement says 99.8% of forest nodes are rank-local at a constant
   fraction in P, so most of this traffic is self-to-self and never hits the wire.
4. Pass 2 loops over **assigned parents only**. The `do i = 1, ntot_g` filter disappears with it: the
   received buffer already contains only this rank's parents, so sorting it once by `kb` replaces the
   per-parent rescan. **That retires W1's pass-2 term in the same change.**
5. Close with **one `ALLGATHERV` of the resulting child BOXES** -- per-box global data, which the
   endstate permits, ~24 B x 2,400 boxes = 58 KB against 144 MB today -- then B1's canonical Morton
   order makes the assembled list identical on every rank regardless of arrival order.

**Why this does not need S3.2b first.** Ledger 39 sequenced S3.2 ahead of S3.3 on the argument that the
forest "performs zero collectives today, so converting it first would ADD ~22-28k collectives". That
modelled S3.3 as making the forest reduce-driven per node. It is not: ownership plus one alltoallv adds
NO per-node collectives, and the residual cross-rank windows are the same ~0.2% the scope instrument
already measures.

**SPLIT INTO TWO GATED INCREMENTS (2026-08-28), because you cannot target the exchange until you know
who needs what:**

- **S3.3a — pass 2 by ownership, keeping the `ALLGATHERV`.** Assign `powner(kb) = mod(kb - plo, num_procs)`
  (a pure function of `kb`, so every rank agrees with no communication); pass 2 processes ONLY its own
  parents. This deletes the per-parent rescan of the whole gathered pair list — `do i = 1, ntot_g` inside
  `do kb = plo, phi`, i.e. **O(parents x global tags) on EVERY rank**, which is W1's pass-2 term — and the
  replicated clustering with it. The global `nboxes + 1 > amr_max_fine` guard cannot be evaluated per rank
  any more, so children are emitted uncapped into a local list and **replayed in (kb, emission) order**
  after ONE `ALLGATHERV` of the child BOXES (7 ints each, ~67 KB against the 144 MB per-cell gather). Each
  parent has exactly one owner emitting in order, so a stable counting sort by `kb` reproduces the serial
  append order exactly — **same box list, same truncation, same `box_level`. Gate: bit-identity.**
- **S3.3b — replace the pass-1 `ALLGATHERV` with `MPI_Alltoallv`.** Now that a rank only needs its own
  parents' tags, bucket the packed pairs by `powner(skb(i))` (stable counting sort), exchange counts with
  `MPI_Alltoall`, then `Alltoallv`. Send volume becomes O(this rank's tagged cells), receive volume
  O(its assigned parents' tags). **`gwin_bytes` stops being O(P).** Pass 2's remaining per-parent scan
  becomes O(parents/P x tags/P) = **O(1) per rank under weak scaling**, closing W1's term completely.






### W5 SCOPED FROM THE CODE 2026-08-28: the wall is the tag BASE, not only the per-box sites

Auditing every AMR point-to-point tag argument (continuations joined, 57 calls in `m_amr.fpp`):

| tag form | calls | verdict |
|---|---|---|
| `tq = amr_tag_base(f) + mod(amr_mesh_epoch, 50 or 100)` | 22 | **already the end-state (family, epoch) form** |
| `tp = amr_tag_base(3) + mod(...)` | 2 | same |
| `amr_cur` (mesh slot index) | 11 | bounded by the slot pool, not the global block count |
| `k`, `4400 + k`, `cblk`, `ks` | ~12 | per-object |

The per-object sites are concentrated, not scattered: the **L0 tile family** --
`s_l0_fill_tiles_from_coarse`, `s_l0_scatter_tiles_to_coarse`, `s_l0_add_reflux_to_tiles`,
`s_l0_restrict_to_tiles` (2 calls each, and all four use BLOCKING `MPI_SEND`/`MPI_RECV`) -- plus the
chunked gather (`s_amr_gather_chunk_post`/`_send`, `s_amr_gather_from_parent_field_*`).

**But the binding constraint is the BASE, and this is the part worth internalising:**

```
amr_tag_base(f) = amr_max_blocks + 100*f          ! m_amr.fpp:839
@:ASSERT(amr_tag_base(size(amr_tag_base)) + 100 <= tag_ub, ...)   ! :849, the tripwire
```

**Every family's tag base carries an `amr_max_blocks` offset**, reserving `[1..amr_max_blocks]` for the
legacy per-box space while families convert. So even a FULLY converted family has a base that grows with
the problem, and the `MPI_TAG_UB` assert caps global blocks near `2**21` ~ 2.1e6 -- **~28k ranks at the
measured ~75 boxes/rank**, which is exactly the wall the module comment records.

**So W5 is two steps, and only the second removes the wall:**
1. Convert the remaining families off per-object tags -- and note this is NOT a tag rename: per-object
   tags exist to match multiple messages between the same pair, so each family needs its messages
   AGGREGATED per peer first, the same shape as I2a/I3/I5. The L0 tile family additionally needs its
   blocking SEND/RECV pairs converted.
2. **Then drop the `amr_max_blocks` term from `amr_tag_base`.** The module comment states the condition
   exactly: "The amr_max_blocks term can only go once NO family uses per-box tags." After that the tag
   space is O(families x epoch) = O(1) and the 28k-rank wall disappears.

Step 2 is one line and is the whole point; step 1 is the work that unblocks it.

### ARM B VERDICT 2026-08-28 (three points, cap raised to 65536 so it never binds): LEVEL 1 IS O(P)

| per regrid | np8 | np16 | np32 | growth |
|---|---|---|---|---|
| level-1 tree nodes = **global `ALLREDUCE`s** | 13,825 | 27,537 | **54,961** | **1.992x, 1.996x = O(P)** |
| level-1 SHARED nodes | 11 | 23 | **47** | 2.09x, 2.04x = **O(P)** |
| `shr_maxdep` | 3 | 5 | 7 | +2 per doubling |
| `capped` warnings | 0 | 0 | 0 | cap never binds |

**This overturns TWO of my earlier readings, both of which were pure cap artifact:**

1. "Level-1 collective count is FLAT in P (16,383)" -- it is **O(P)**. The 16,383 was `2 x 8192 - 1`, a
   full binary tree pinned by `amr_max_blocks`, identical in every rung by construction.
2. "The level-1 shared set grows as O(log P) (+4 per doubling: 7, 11, 15)" -- it is **O(P)** (11, 23, 47).
   That sequence was measured on cap-pinned trees, so it was reporting how a FIXED tree splits across more
   ranks, not how the tree grows.

**Consequences.**
- **S3.2b is not a latency constant.** It removes an O(P) count of global collectives -- 54,961 per regrid
  at np32 and doubling with P. It is as important as S3.3, not a follow-on.
- **The reordering argument in ledger 41 is dead.** It rested on "the forest is O(P) while level 1 is
  flat". Both are O(P). S3.3 was still the right thing to land first -- it removed a *measured* 1.44 GB
  per rank of O(P) volume with a bounded, bit-identical design -- but the justification is now absolute
  size and readiness, not a difference in exponent.
- **Depth-batching becomes essential rather than an optimisation.** With shared nodes O(P), an
  un-batched sparse exchange would be O(P) collectives; batching by depth gives `shr_maxdep` (~2 log2 P,
  so ~31 rounds at 1e5 ranks) regardless. The remaining question S3.2b must answer is the shared VOLUME
  per rank, which `[amr-scope-me]` now instruments.
- **B0b is more important than recorded.** At the shipped `bf = 1` the cap binds on every regrid, which is
  precisely why this was invisible for so long.

**Method note.** Three "flat in P" readings tonight were a pinned limit, and the tell each time was a
round number (16,383 = 2 x 8192 - 1). Added to the three-point rule: before believing a flat result, ask
whether a cap, pool or table was fixed across the sweep -- and check whether the value is suspiciously
round.

### W1 SIZED PROPERLY 2026-08-28: ~48 full scans of the GLOBAL block list PER STEP

Following the replicated-metadata finding, an audit of every `do ... = 1, amr_num_blocks` in the AMR
sources: **44 sites across 22 routines.** The ones that matter are on the RK stage path, called from
`m_time_steppers.fpp` (`s_amr_stage_fill_wave` :579, `s_amr_parent_fill_wave` :581 per level,
`s_amr_reflux_faces_wave` :633, `s_amr_freg_wave` :819), not in the regrid:

| routine | scans | per |
|---|---|---|
| `s_amr_stage_fill_wave` | 2 | stage |
| `s_amr_parent_fill_wave` | 2 | stage x (levels-1) |
| `s_amr_reflux_faces_wave` | 3 | stage |
| `s_amr_freg_wave` | 3 | stage |
| `s_amr_restrict_l1_wave` | 3 | stage |
| `s_amr_restrict_parent_wave` | 2 | stage |
| `s_amr_convert_prim_batch` | 1 | stage |

~16 scans per stage at `amr_max_level = 2`, so **~48 per step under RK3**. Every one walks the whole
machine's block list and filters to this rank's own ~75:

```
do k = 1, amr_num_blocks
    if (amr_block_level(k) /= 1) cycle
    if (amr_block_owner(k) == proc_rank) cycle     ! or /= , for the send side
```

At np32 that is ~2,500 blocks x 48 = 120 k predicate evaluations per step -- invisible. At 1e5 ranks with
the measured 75 blocks/rank it is 7.5 M blocks x 48 = **~360 M evaluations per rank per step**, roughly
0.4-0.7 s of pure metadata scanning per step, growing linearly in P while the real work stays fixed.
**This is the direct W1 violation** -- "per-rank STEP cost = f(local cells, peers)" -- and it is a bigger
practical problem than the regrid pass-1 O(P^2), because it runs every step rather than every 20th.

**Two fixes, and the cheap one is mechanical.**

- **W1a (cheap, no architecture change): an owned-block index list.** Maintain `amr_my_blocks(:)` /
  `amr_n_my`, rebuilt once per regrid, and convert the ownership-filtered loops to
  `do kk = 1, amr_n_my; k = amr_my_blocks(kk)`. **10 sites are directly convertible** (grep-verified:
  reflux_faces_wave, freg_wave x3, restrict_l1_wave x2, stage_fill_wave, parent_fill_wave x2,
  convert_prim_batch), plus `restrict_parent_wave` whose filter is written `proc_rank /= cowner`. The
  RECEIVE-side loops (`amr_block_owner(k) == proc_rank` -> cycle) need the complementary list: the peer
  blocks this rank participates in, which is exactly the "peers" W1 allows -- build it in the same pass.
- **W1b (architectural): distribute the metadata itself.** See the root-cause section above. W1a takes
  per-step cost to O(local + peers) but leaves the O(P) MEMORY; only W1b removes that.

**Free win found in the same read:** `s_amr_reflux_faces_wave` BUILDS `amr_fw_rblk` inside its first
global scan (`amr_fw_rblk(nhr) = k`, :2744) and then its SECOND global scan (:2870) walks the whole block
list AGAIN purely to re-derive the same index order, asserting agreement
(`@:ASSERT(amr_fw_rblk(j) == k, ...)`). The same shape repeats in `s_amr_freg_wave` (:2930 / :3042).
**Those second scans are removable outright** -- iterate `j = 1, nhr` over `amr_fw_rblk` directly. That is
4 of the ~48 scans deleted with no new data structure and no behaviour change.

### W1's ROOT CAUSE, found 2026-08-28 while scoping S3.3c: the block metadata is REPLICATED GLOBALLY

Scoping S3.3c sent me to `5ec798ed` (the seam all-pairs fix) for reusable machinery. It does not transfer
-- that fix worked because `f_amr_seam`'s predicate FIXES the neighbour's `lo` corner exactly, making it a
lookup rather than a search, whereas pass 1 needs a genuine range-overlap query. But reading it led
somewhere more important:

```
allocate (amr_region_lo_all(3, amr_max_blocks), amr_region_hi_all(3, amr_max_blocks))
allocate (amr_owns_all(amr_max_blocks))
allocate (amr_block_level(amr_max_blocks))
```

**Every rank holds the geometry, ownership and level of EVERY block in the machine.** That is O(P) memory
per rank by construction -- at 1e5 ranks with the measured ~75 blocks/rank, ~7.5 M blocks x 8 ints x 4 B
is **~240 MB per rank of pure metadata**, before any of it is scanned. The O(P^2) pass-1 scan is a
SYMPTOM of this: a loop over all blocks is only writable because all blocks are locally addressable.

**So W1 cannot be fully met by restructuring loops.** Every "iterate `1, amr_num_blocks`" site in the
regrid is bounded below by this array's existence. Closing W1 needs the metadata itself distributed --
each rank holding its own blocks plus a halo of neighbours, with remote queries going through the
existing owner map (`f_amr_owner`'s Morton cut-point search) instead of a local array read. That is an
architectural increment on the scale of the store flattening, not a loop fix, and it should be scoped
separately before any more W1 loop work is done.

**Revised S3.3c (still worth doing, and small).** It takes the pass-1 scan from O(P^2) to O(P) without
touching the metadata layout:
- Compute `mine(kb)` -- does this rank hold data for parent `kb` -- by iterating this rank's OWNED blocks
  only, which is O(local blocks) = O(1) per rank.
- Compute `covered(kb)` the same way, then **one `MPI_ALLREDUCE(MPI_LOR)` over the parents array**
  (~2,400 logicals, per-box global data, which the endstate permits) instead of every rank scanning every
  block. This is the step that removes the O(P^2).
- Guard the dense `gwin` allocation on `mine(kb)`, which removes the O(P) window-zeroing traffic for free.

Residual after S3.3c is O(parents) per rank -- the loop over `plo:phi` itself -- which is exactly the
replicated-metadata bound above and cannot go lower until that is fixed.

### S3.3c FOUND 2026-08-28 (reviewing my own S3.3a/b diff): pass 1 is STILL fully replicated, and it is O(P^2)

S3.3a/b fix the gather and pass 2. **They do not close W1.** Re-reading pass 1 after writing them, every
rank still does this for EVERY parent, regardless of ownership or of whether it holds any relevant data:

```
do kb = plo, phi
    allocate (gwin(mlo(1):mhi(1), mlo(2):mhi(2), mlo(3):mhi(3)))   ! dense window, per parent, per rank
    gwin = .false.
    do ob = 1, amr_num_blocks                                       ! the REPLICATED GLOBAL block list
```

Two residual terms, neither touched by S3.3a/b:

1. **`allocate`+zero of a dense window per parent** -> O(parents x window volume) = **O(P)** memory traffic
   per rank per regrid.
2. **`do ob = 1, amr_num_blocks` inside `do kb = plo, phi`** -> O(parents x global blocks). Parents scale
   with P and the global block list scales with P, so this is **O(P^2) per rank** -- asymptotically the
   WORST term in the regrid, worse than the `gwin_bytes` gather S3.3b removes. It is only cheap today
   because it is a metadata test (~5.8 M iterations per level per regrid at np32).

**Why it cannot simply be skipped:** `covered(kb)` is replicated metadata -- every rank must agree on it --
and it is set inside that same `ob` loop, so the loop cannot be made owner-only as written.

**S3.3c, its own increment.** Two independent fixes:
- **Invert the loop.** Iterate this rank's OWNED old blocks and tag the parents each one overlaps, instead
  of iterating all parents x all blocks. Cost becomes O(owned blocks x parents touched).
- **Compute `covered` once, spatially.** This is the same shape as the seam all-pairs scan that W1 records
  as already FIXED (`5ec798ed`), so the binning machinery to do it exists and the precedent is set.
- Then guard the `gwin` allocation on "this rank actually contributes", which removes term 1 for free.

**Do not report S3.3 as closing W1 until S3.3c lands.** S3.3a removes W1's pass-2 rescan
(O(parents x global TAGS)); the pass-1 terms above are what remain.

**Two implementation details verified by reading the code, because both are easy to get wrong.**
(1) The `Alltoallv` needs the pairs contiguous per destination, but `s_amr_pack_gwin_pairs` appends them
grouped by `kb` in increasing order. A stable counting sort by `owner(skb(i))` fixes that in O(n) and is
the only reordering needed. (2) **That reordering is safe**: pass 2 rebuilds a DENSE `gwin` from the
pairs and re-extracts in `(k,j,i)` order, so `ctags` is already independent of the order pairs arrive in
— the existing comment's "setting .true. once per gathered cell reproduces the old per-parent dedup" is
exactly that property. So S3.3 does not perturb `ctags` and its box set should be bit-identical to the
post-B1/B0b baseline, which makes the gate a clean signal.

**Gate.** `gwin_bytes` falls from 1.44 GB to O(local) and stops doubling with P (the three-point sweep
re-run is the proof); child box set bit-identical to the post-B1/B0b baseline; 69/69 AMR goldens; and the
suite ALREADY carries the right np>1 gate — verified 2026-08-28 by reading `cases.py`: case (m)
**"AMR -> 1D -> multi-level dynamic regrid np=2"** runs `amr_regrid_int=2, amr_max_level=2` at `ppn=2`,
which is precisely the cross-rank dynamic nesting path that carries the `ALLGATHERV`. Cases (l), (l'),
(n) and the `amr_max_level=3` case add static and 2D/3-level coverage at `ppn=2`. So S3.3 does NOT need
a new test — a single-rank golden could not see a broken alltoallv, but these can.

### S3.2b DESIGN, REVISED 2026-08-28 — depth-batched, and much smaller than the version above

The measurement changes the design. With 99.93% of level-1 nodes rank-local and only 7-11 shared per
regrid, the expensive part is not making the shared exchange sparse — it is *not doing the other
16,372*. That splits the clusterer in two, and the shallow half can stay on `MPI_COMM_WORLD`:

- **Shallow phase (`novr > 1`), replicated.** Every rank already agrees on this part of the tree: it
  descends from the global bbox, and `s_amr_ranks_overlapping` is a pure function of the box and the
  decomposition, so `novr` needs no communication to evaluate. Walk it **breadth-first by depth** and
  issue **ONE fused `ALLREDUCE` per depth** carrying the concatenated signatures of every shared node
  at that depth. All ranks hold identical shallow state, so they agree on the node list and hence on
  the concatenation layout. Collectives per regrid become `shr_maxdep` — **3 at np8, 5 at np16** —
  against 16,383 today.
- **Deep phase (`novr == 1`), private.** The moment a node fits inside one rank, that rank holds every
  tag in the subtree: it finishes the subtree alone with **zero communication**, and no other rank
  descends into it at all.
- **Close with one `ALLGATHERV` of the resulting BOXES** — per-box global data, which the endstate
  explicitly permits, and which B1 makes order-independent.

**What this drops from the earlier design.** No graph communicator, no `MPI_Neighbor_alltoallv`, no
per-node sparsity. Those were aimed at making each shared node's exchange cheap; depth-batching instead
makes the NUMBER of exchanges `O(shr_maxdep)` regardless, and the measured volume that has to move is
~35 KB per regrid, which no amount of sparsity would meaningfully improve. Sparsity becomes a later
refinement, not a prerequisite.

**Why this survives growth.** Shared nodes grow ~1.57x per np-doubling (~P^0.65, the rank-boundary
SURFACE, exactly as geometry predicts), so ~3,000 shared nodes per regrid at 1e5 ranks. Un-batched that
would be 3,000 collectives; batched by depth it is `shr_maxdep` ~ 2 log2 P ~ **34**. This is why
`shr_maxdep` still matters and why its growth rate is worth knowing -- but it sets the ROUND count, not
the volume, and 34 rounds is not a wall.

**Gate.** Level-1 global collectives per regrid = `shr_maxdep` (not 16,383); box set identical to the
post-B1 baseline; `[amr-tree] nodes` unchanged; `[amr-scope-me]` confirming the deep phase is genuinely
private. Bit-identity against post-B1 is achievable here because B0 made `force`/the cap inert and B1
made the merge order-independent -- which is the whole reason those two landed first.

### Increments

- **S3.2a (next, cheap):** measure the depth at which nodes become rank-local, and the node/volume
  split above versus below it. The depth instrumentation already exists. **This sizes the whole
  design before any of it is written.**
- **S3.2b:** shallow phase depth-batched, deep phase local — see "S3.2b DESIGN, REVISED 2026-08-28"
  below, which supersedes the neighbour-communicator sketch above and the bytes-flat gate (level-1
  per-rank bytes are ALREADY flat, so that gate cannot discriminate).
- **S3.3:** the same mechanism on the level-2 forest.

## 2026-08-27 (38) — S3.1 LANDED: THE LEVEL-1 W4 COLLECTIVE IS GONE, AND THE SCALING LAW MEASURED

**S3.1.** `s_amr_union_gtag` and its ALLGATHERV are deleted, replaced by `s_amr_local_tags` (each
rank scans only its own disjoint interior) plus ONE fused `MPI_ALLREDUCE` per tree node carrying the
per-axis signatures. Reading all three routines showed one reduction suffices for trim, count AND
split: `s_amr_trim_box` takes the per-axis MIN/MAX of contained tags, which IS the first/last nonzero
of the signature, and its ok test IS sum > 0. Gated four ways: restart output byte-identical, tree
reproduced exactly (`lmax 8192 ldepth 31 nodes 66354` in both arms), **`ntag_bytes` 185,241,240 -> 0**,
and **69/69 AMR goldens** (which cover the 1D/2D, np=1 and multi-level paths the np8 probe cannot).

**The scaling law (np8 -> np16, matched at 10 regrids).** ntag 92.77 -> 185.53 MB/rank/regrid =
2.000x per doubling. rbytes 4.53 -> 7.53 MB = **1.660x = P^0.73**. nboxes 594 -> 1206. So S3 improves
the COEFFICIENT (~20-25x) **and** the EXPONENT — better than the audit feared — but is **not flat, so
W4 is not satisfied**. At 1e5 ranks: ntag would have been ~1.16 TB/rank/regrid, rbytes ~4.4 GB.
Mechanism confirmed twice over: `rbytes ~ nboxes^(2/3)` predicts 1.60x against 1.66x measured.

**THREE-POINT LADDER, and it OVERTURNS the two-point law.** np8/16/32 at matched 10 regrids:
ntag 92.77 -> 185.53 -> 371.04 MB/rank/regrid (2.000x, 2.000x -- clean P^1.0). rbytes 4.53 -> 7.53
-> 13.75 MB = **1.660x then 1.827x**, i.e. the implied exponent RISING 0.73 -> 0.87. Decomposing the
node count settles why: the level-1 tree is CONSTANT at 16,383 nodes (cap-fixed at 8192 leaves) while
the forest goes 28,308 -> 55,872 -> 111,112, i.e. **1.97x, 1.99x -- exactly proportional to P**. Total
= constant + O(P), so the constant dilutes and the exponent climbs toward 1.0.

**So S3.1 is asymptotically O(P): a COEFFICIENT cut (20.5x -> 24.7x -> 27.0x, growing slowly but
converging), NOT an exponent improvement.** The earlier "P^0.73" was a small-np artifact of the
cap-fixed level-1 term and is RETRACTED -- it was published on two points, and the third refuted it.
The standing rule (>=3 points at fixed per-rank work before any asymptotic claim) is what caught it.

**Consequence for S3.2/S3.3.** The forest is the term that grows, so it is where W4 is actually won;
but it performs zero collectives today, so the split-communicator mechanism must still be PROVEN on
the level-1 path first and then applied to the forest by S3.3. Both statements hold at once.

**What `rbytes` actually measures, corrected.** `amr_cl_rb` is incremented OUTSIDE the `if (reduce)`
guard, so it prices every tree node in BOTH paths -- but only the level-1 path reduces today. It is
therefore a **projection of the fully distributed design (S3.1 + S3.3), not today's wire volume**.
That is the right quantity for projecting W4 and the wrong one to call a measurement of current cost.
Today's actual added collectives are the level-1 tree only: ~16,383 per regrid at cap 8192, not the
44,691 total node count.

**Three corrections to my own earlier claims, all from measurement:**
1. **The unfused collective count is 44,691/regrid at np8 and 72,255 at np16** — the design estimate
   was ~1,200, wrong by 50x. That is what makes S3.2 fusion mandatory rather than an optimisation.
2. **S3's value is MEMORY, not wall.** `rg:clus` is 0.6 percent of wall; `regrid` is 53.7 percent but
   `rg:build` is 37.8 percent of it. My guess that cap saturation might explain regrid's share is
   **WRONG and retracted** — saturation drives node count, and node count lives in rg:clus.
3. **The block cap is the WRONG KNOB.** Five caps at fixed domain: 16x cap gives 2.10x nodes and
   2.10x reduction bytes for a **3.7 percent** change in the box set (nboxes 598/598/598/594/576).
   The cap bounds only the level-1 tree; the forest is parent-driven and cap-independent, measured
   FLAT at ~21,800 nodes/regrid across all five caps. Absolute split at cap 8192: level-1 ~16,383
   (43 percent), forest ~21,872 (57 percent) -- at cap 1024 the forest is 91 percent. **An earlier
   "~85 percent" figure was a marginal-increment reading, not the absolute split; corrected.**

**The clusterer saturates EVERY cap** (warning on every regrid at all five settings): with
`amr_cluster_eff = 0.9` and no minimum box size the bisection never converges on its own, it splits
until something stops it and the merge undoes most of it. `rg:clus` grows 13x across the sweep while
nodes grow 2.1x — superlinear, consistent with the O(n^2) merge running at n = cap. **That is B0's
case: not wall share, but waste and the lever on nboxes.**

**AMReX does not solve W4 either** (read from source): `TagBoxArray::collate` Gathers every tagged
cell to the I/O rank, `ClusterList::chop` runs SERIALLY there, and the BoxList is broadcast — with a
hard `amrex::Abort` above INT_MAX tags, marked `xxxxx todo`, whose suggested mitigation is a larger
blocking factor (i.e. S1). There is no distributed clustering in the SOTA to copy, and on this axis
MFC was WORSE: every rank received what only AMReX's root did.

**Deferred deliberately: the IB ownership gate.** It is correct and reduces `amr_gb_win`, but it
lives behind `amr_max_level >= 2` and every case available runs at level 1, so no test could execute
it. Shipping it unverified was the wrong trade; it goes with S3.3, which needs it. The new ppn=2
AMR+IB case (`E4F6CE1E`) landed anyway — it is the only distributed AMR+IB coverage in the suite.

**Process failures and the guards they bought.** (a) Two builds ran concurrently in one tree; the
survivor was probably correct but not provably so. (b) A revert left unbalanced `end do`/`end block`
that `ffmt` and `lint_source` BOTH passed — neither parses Fortran — so the compiler was the first
thing that could say no, and the failing build still reported success to a wrapper. (c) My own
binary-provenance guard passed a STALE artifact, because symbol presence is necessary and not
sufficient; the discriminating question is "is this a different artifact than the last pin". Guards
added: RUNBOOK 3b/3c, `format.sh` now surfaces ffmt stderr and fails on unmatched structure, and
`amr-bench/fcheck.sh` structurally checks one .fpp in ~4 s versus ~18 min.

## 2026-08-27 (37) — AUDIT OF THE S3 PLAN: A GATE THAT ONLY TESTED ONE DIRECTION, AND AN UNVERIFIED DEPTH ASSUMPTION

**The gate on ledger (36)'s migration was one-sided.** 15/15 tests confirmed that VALID cases still
run; nothing confirmed that INVALID cases are still rejected. Having deleted 53 checks, that is only
half a gate: every one of them could have been a no-op and the suite would still be green.
Negative-testing closes it - perturbing `examples/2D_amr_droplet` and running `./mfc.sh validate`:
`time_stepper = 1` -> "amr requires time_stepper = 3" (the B8 fix); `cfl_const_dt = T` alongside
`amr_subcycle` -> "amr_subcycle requires a fixed dt", **which is exactly the case the OLD validator
let through and the solver aborted on** - the derived-`cfl_dt` fix is doing real work. A fourth
probe (`amr_max_level = 3`) correctly reported nothing: the deleted check bounded STATIC AMR
(`amr_regrid_int == 0`) and that example sets `amr_regrid_int = 10`, so three levels is legal.
**STANDING RULE: deleting or relocating a check requires a negative test. A suite that only proves
good cases still pass cannot distinguish a working check from a deleted one.**

**S3's central performance claim is UNVERIFIED, and it is now the first increment.** The design
argues fusion gives "~17 fused collectives per regrid" from tree depth being O(log nboxes). But
Berger-Rigoutsos splits at signature holes and Laplacian inflections, not midpoints, so the tree is
not balanced by construction; its worst case peels one box at a time, giving depth O(nboxes) and
reinstating the latency blow-up fusion exists to prevent. **Nothing in this ledger measures BR tree
depth.** So S3.0a is now instrumentation, not code: carry a depth per stack entry in `s_amr_cluster`,
report max depth and node count per regrid, run it over the S0 arms. Pre-registered rule: depth
within ~`4*log2(nboxes)` -> fusion proceeds as designed; depth growing like `nboxes` -> the design
needs a different shape (cap the depth, batch sibling subtrees, or fall back to S1 as a coefficient
stopgap). **Do not write S3.2 before that number exists.**

**Also marked INFERRED rather than established:** the claim that S3 deletes the remaining W1
quadratic (the regrid pass-2 loops over parents x global tagged cells). That was reasoned from their
structure, not from reading them. Confirm at `m_amr_regrid.fpp:1118` and `:1130` before relying on
it; if those loops need the global list for more than locating each parent's slice, W1 keeps a
separate increment.

**Ordering adopted:** S3.0a (depth measurement, go/no-go) -> S3.0b (emit reduced-buffer bytes) ->
S3.1 (path 1 on reductions, judged on bit-identity ONLY, expected slower) -> S3.2 (fusion, blocked on
S3.0a, watch the `amr_max_blocks` force path whose `nacc + nwork + 1` test means different things
under a LIFO stack and a depth level) -> S3.3 (path 2 as a fused forest) -> S3.4 (delete the
gathers). Then W5 (I7), W3 (F5), W2 (2b), W6/W7. The L0 deletion and the post-process reader fix are
unblocked but not on the exascale path; the upstream `s_check_inputs_weno`/`_muscl` marker issue is
to be reported, not patched here.

**A tooling pattern worth naming, five instances in one day.** Each was a check that looked like
verification while measuring something else: ``pgrep -fc 'a\|b'`` matching nothing (ERE, so `\|` is a
literal pipe) read as "the job died"; two `grep -c` calls whose zero-match exit status read as
"build failed"; a `-x` test performed before the `cd` that would invalidate the path; and a negative
probe whose injected text broke Python syntax, so all four cases "failed validation" for the wrong
reason. None were reasoning errors. **Rule: when a check misleads, harden the check - do not resolve
to be more careful.** Done for `bitcmp_probe.sh` (absolute paths; "did not run" now distinguished
from "differs") and for the scorecard (documents the per-rank and double-count traps, refuses
mismatched arms). The counts themselves are never the evidence: print what matched.

## 2026-08-27 (36) — THE HYGIENE BATCH: CHECKER MIGRATION, STACK AUTOMATICS, AND THE SEAM SCAN DE-QUADRATIFIED

Four commits, all byte-identical or golden-holding, all gated.

**`82dd12f5` - the merge regression fixed.** Upstream `55fb1b14` migrated input-only checks out of
`m_checker` into `case_validator` and added a lint rule to enforce the split; our merge resolved the
conflicted region as "ours", restoring 19 of them, then allowlisted `s_check_inputs` so the new rule
could not fire. The 19 are deleted, with the now-unused `m_helper_basic` import and the already-dead
`muscl_order_first_order` one. One of the 19 was doing real work the validator did not replicate:
the Fortran `pign` check fires on the `dflt_real` SENTINEL while the validator fired only on an
ABSENT key, so a case writing `rburn%%pign` as the sentinel passed `./mfc.sh validate` and then
ignited everywhere from t=0. The validator now rejects any non-positive `pign`, matching the pattern
`rburn%%k`/`rburn%%pref` already use - **that fix landed BEFORE the Fortran check was removed.**
Gates: 10/10 (three reactive-burn + chemistry AMR), precheck.

**`f702fb27` - AMR input constraints moved to Python** (user: "i like the AMR prohibits in the mfc
toolchain (python) when possible"). 53 of the 57 prohibits in `s_check_inputs` deleted; the 4 that
survive all read `num_procs` (MPI_Comm_size, invisible to the validator) and each carries a
`! lint: runtime-check` marker. With only those left, `s_check_inputs` came OFF
`RUNTIME_CHECKER_SUBROUTINES` - **upstream's rule is armed here for the first time**, so a new
input-only check added Fortran-side is now flagged. Nine validator checks were corrected FIRST,
because deleting their counterparts would otherwise have opened silent holes: `cfl_dt` and
`bodyForces` are DERIVED in `s_read_input_file` (from `cfl_adap_dt`/`cfl_const_dt`, and from
`bf_spatial_support`) while the validator read only the literal keys - so `amr_subcycle` with
`cfl_const_dt` passed validation and aborted in the solver; `time_stepper`/`model_eqns` default to
the `dflt_int` sentinel, which Fortran rejects and the validator skipped on an absent key; and
`amr_tag_eps`/`amr_buf` were over-strict the other way, erroring on unset keys whose Fortran
defaults are valid. **Trade recorded:** the validator runs on the `./mfc.sh run` path, so a built
binary invoked directly against a hand-edited `simulation.inp` no longer gets these checks.
**Also removed undefined behaviour:** `s_check_inputs` read `num_dims` - not assigned until
`s_initialize_parallel_io_common` 25 lines later - and sliced `amr_block_beg(1:num_dims)` with it.
Gates: 15/15 (AMR + active_box + coexist + load-balance + diagnostic writers), all examples valid.

**`4e21d0f3` - W1, the crash half.** Twelve O(global boxes) AUTOMATIC arrays moved to the heap:
seven in `s_amr_assign_block_owners` (~48 bytes/block between them, ~48 MB at 1e6 blocks) and five
in `s_amr_regrid`. These OVERFLOW a default stack long before the box count itself becomes the
bottleneck - a crash, not a slowdown. Gates: bitcmp byte-identical, AMR 10/10.

**`5ec798ed` - W1, one of the two quadratic paths.** `s_amr_build_seam_pairs` ran two nested loops
over `amr_num_blocks` per regrid (count then fill) = 2N^2 in the GLOBAL block count. **It was never
a search:** `f_amr_seam` demands exact equality of region lo AND hi on both transverse dims plus
`lo(d,yb) = hi(d,xb)+1`, and blocks are disjoint, so `(level, lo)` names a block uniquely and those
conditions FIX the neighbour's lo corner given `xb` and `d`. Morton-sort the lo corners once (the
stable merge sort `s_amr_sfc_cut` already uses - whose own comment makes exactly this O(n^2)
argument about its predecessor), then binary-search the single candidate per (block, dim) and verify
the predicate: **O(N log N)**. Emission order preserved exactly (xb ascending, each block's <=3
matches insertion-sorted into ascending yb) because the list drives paired MPI_SENDRECVs - a
reordered list DEADLOCKS rather than producing wrong numbers. Equal-key runs are rescanned against
the full lo vector (f_morton keeps 21 bits/dim). Gates: bitcmp byte-identical; 17 AMR tests covering
tiled seams, tiled L2, multi-block, three-level, np=2 multi-level, and the coexist set that
exercises the l0 periodic-wrap branch.

**Remaining in W1:** the regrid nesting loops (`m_amr_regrid.fpp`), whose pass 2 is parents x GLOBAL
TAGGED CELLS - a worse shape than boxes x boxes. **Do not fix it standalone:** it iterates the global
tag list that S3 deletes, so it is subsumed by the W4 work.

**Upstream issue to file (NOT patched here):** `s_check_inputs_weno`/`_muscl` carry
`! lint: runtime-check m/n/p are per-rank extents after MPI decomposition` markers, but they are
called from `s_check_inputs` at `m_start_up.fpp:1070`, BEFORE
`s_mpi_decompose_computational_domain` at `:1097` - so m/n/p are the GLOBAL case-file values there
and those checks are input-only, not runtime. Upstream's checks, upstream's markers: report it,
do not quietly change upstream checker semantics inside a merge branch. That is how the regression
in the first paragraph happened.

## 2026-08-27 (35) — THE LADDER IS NOT THE SCORECARD: W1/W2/W4/W7 UNMET; D-l0 = DELETE; D-phase2 = NOT NEXT

**Retraction.** (34) concluded "declare the performance program done," reasoned from the six-rung
matched-point ladder under a parity framing. The user has confirmed the goal is the end-state
architecture and exascale scaling ("this is just the computer I have for you to use"). The ladder
measures wall degradation per np-doubling on ONE node and does not measure any of W1-W7's defining
quantities except W3/W4 - both of which the same jobs' counters show still growing. The scorecard is
`amr_endstate.md` section 3.

**Invariants re-audited against code (the section 3 "today" column was a week stale).** W1 **UNMET
and worse than documented**: two O(global boxes^2) paths - `m_amr.fpp:6154-6170` (seam-pairs, nested
`do xb`/`do yb` over `amr_num_blocks`, run twice per regrid) and `m_amr_regrid.fpp:1017,1033,1118,1130`
(nesting; pass 2 is parents x global tagged cells) - plus per-stage O(global boxes) plan walks and
O(boxes)/O(P) **automatic STACK** arrays (`m_amr_regrid.fpp:638-642`, `m_amr.fpp:3188-3190,2687-2689`)
which overflow (crash) before anything merely slows. W2 **UNMET, structurally unchanged**: the advance
is still one `s_compute_rhs` per box per stage; 2a is compiled but hard-gated (`amr_prim_batch = .false.`,
`m_amr.fpp:437`). W4 **UNMET**: `s_amr_union_gtag` still ALLGATHERVs every tagged cell globally and
allocates `allidx` at the GLOBAL count (`m_amr_regrid.fpp:438-444`). W5 **PARTIAL**: per-step families
use (family, epoch) bases but the base is still offset by the global cap
(`amr_tag_base(f) = amr_max_blocks + 100*f`, `m_amr.fpp:803-804`), with a loud init assert at ~2.1M
blocks on Cray. W6 **PARTIAL**: store growth reverts to a full HOST round trip above 32 columns
(`m_amr.fpp:7897`) and ranks hold 72-75 boxes, so the host path is the one that runs. W7 **UNMET**:
weight is work-only, no migration term, no hysteresis. W3 **PARTIAL** (per-box WAITALL chain dead;
F5 still per-(box,face,participant)). W8 holds through np32.

**W4 is the top blocker on a STRUCTURAL argument, not a measured coefficient.** Each rank receives
every rank's tagged cells, so at fixed per-rank work the received volume is O(P) per rank - derivable
from the code, and the reason it cannot be tuned away. **NOW MEASURED INDEPENDENTLY AND CONFIRMED**, via the new
`amr-bench/invariant_scorecard.py` on the matched np16/np32 arms of job 386892 (9 regrids each,
cells/rank constant at 8M): ntag **39.29 -> 78.59 MiB per rank per regrid = 2.00x**, gwin
8.60 -> 17.06 = 1.98x. Extrapolating the confirmed law from np32: **~240 GiB per rank per regrid at
1e5 ranks.** Two traps the tool now documents, because both produced wrong readings during this
audit: `amr_gb_tag` accumulates the GLOBAL tag list and is printed from rank 0, but every rank
RECEIVES that list, so it is ALREADY per-rank - do not divide by np; and a `**/*.log` glob already
matches `run.log`, so scanning both double-counts regrids. The tool refuses to compare arms with
mismatched regrid counts rather than silently averaging them.

**THE INVERSION.** Phase 2b - "the largest single investment in the program" - addresses W2, which is
FLAT in P (75 boxes/rank at np8/16/32 alike). It is a ~173x launch constant: a permanent efficiency
tax at every scale, but not what fails at 1e5 ranks. **Corrected order: W4 -> W1 (de-quadratify + heap the stack
automatics) -> W5 (I7) -> W3 (F5) -> W2 (2b) -> W6/W7.**

**CORRECTION to this entry's first draft: S1 is NOT the W4 fix.** W4 requires the tag exchange be
O(local + peers). Block-lattice coarsening shrinks each rank's contribution but leaves the
ALLGATHERV in place - every rank still receives every rank's tags, so received volume stays O(P).
It moves the curve down ~512x (about nine doublings of headroom, reaching ~1e4 ranks from np32);
it does not change the slope. Its own pre-registered gate - "judged on ntag bytes/rank going flat"
(line ~1464) - is therefore UNSATISFIABLE by it: under weak scaling the domain grows with P, so a
coarsened lattice grows with P too. Nobody noticed because S1 was never attempted.
**The genuine W4 fix is architectural, and the code says why.** `s_amr_union_gtag`
(`m_amr_regrid.fpp:408-412`) all-gathers so that "every rank decodes the same gathered index set
into the same list, so the bisection is rank-invariant" - the collective buys RANK-INVARIANCE of
the box decomposition, which the goldens depend on. A neighbor-scoped exchange breaks that
property outright: each rank would cluster from a different tag subset. The fix is therefore
**S3 (local clustering + boundary reconciliation), already scheduled in endstate Phase 3** - pull
it forward, with determinism as its hardest requirement, not as an afterthought.
**Recommendation: SKIP S1 and go straight to S3.** Coarsening the tag lattice changes which cells
are tagged, hence the box set, hence the goldens; S3 changes the box set too. Doing both means two
golden churns for one architectural outcome. S1 remains available as a stopgap ONLY if a deadline
demands headroom before S3 is ready, and must then be labelled a coefficient fix with an explicit
headroom number. D-phase2 resolves: 2b stays on the roadmap
as the parity item, NOT as the next increment. Ordering it later is a claim about *when*, not about
whether it is needed.

**Why 2b is so large (scoped):** the per-box advance SWAPS module-global grid state -
`s_amr_swap_to_fine` (`m_amr.fpp:5270`, 17 call sites) replaces m/n/p, idwint, idwbuff, the whole
x_cb/x_cc/dx family, and feature flags - then calls the ordinary `s_compute_rhs`. So ~270 GPU parallel
regions across ~10,500 lines read grid extents from module globals. The swap is also a CORRECTNESS
contract (it disables coarse-indexed `acoustic_source`/`ab_active` and depth-guards nesting) that any
box-indexed rewrite must reproduce explicitly. Ledger (27)'s 2a bridge-load regression IS that swap
boundary appearing as runtime cost - so any increment batching kernels WITHOUT removing the swap will
regress.

**D-l0 = DELETE.** MFC's level 0 is balanced BY CONSTRUCTION: a Cartesian decomposition hands every
rank one equal-sized chunk. AMReX boxes level 0 only because there boxes ARE the decomposition, which
is also why the literature's ">=4 boxes/rank" floor does not transfer - MFC meets its intent with one
*exactly equal* box per rank. Tiling costs 28-35% wall (three reps) and rebalancing recovers 0.6%,
inside noise. Every surveyed framework ships base-level rebalancing off or slow, and over-decomposition
is expensive on GPUs. **Correction to the study that produced this: its claim that coexist computes the
level-0 RHS twice is WRONG** - `m_amr.fpp:8569` states the monolithic L0 RHS is skipped for
`l0_ntile>0` (which is exactly why the spike aborts on `run_time_info`/`probe_wrt`). The conclusion
survives without that pillar. **RE-ENTRY CONDITION, and it is live rather than hypothetical:** a
measured level-0 work imbalance above ~10% on an IB, chemistry, or Lagrangian-bubble case, with
`[amr-balance]` extended to level 0. MFC's own fine-block cost model already carries `K_ib`/`K_pc`
terms - the code encodes that per-cell cost is NOT uniform once those physics are on; level 0 simply
gets no equivalent treatment. Before deleting, check what the 12 L0/coexist tests cover incidentally.
Consequence: the L0 restart-writer filter is DROPPED - do not fix save/restart for machinery being
removed.

**Landed here: the last stage ifdefs are gone from src/common (`caabd06d`).** `load_weight_wrt` and
`sfc_partition_wrt` are registered for all three targets with their default in
`s_assign_common_defaults`; `amr_in_fine_advance` is declared in `m_global_parameters_common`. `src/`
now has ZERO `MFC_PRE_PROCESS`/`MFC_SIMULATION`/`MFC_POST_PROCESS`, matching master. Net -4 LOC. Also
fixes a latent gap where post_process declared both flags with no initializer and no default. Gates:
all three targets build clean (pre_process is the real gate), 71/71 (AMR-68 + diagnostic writers
`4D5E2869` + load-balance + SFC), precheck clean.

**Folded in, per the standing rule that a pending verification lands with the next change: the np8
wall-neutrality A/B for participation-local registers is RESOLVED and NEUTRAL.** Differenced 240-40
arms on k004-005: reglocal 7.836 s/step vs mergereplay 7.652, +2.4% - inside the ~5% noise floor. The
`[amr-xa]` exchange audit is byte-identical between arms (same msgs, same words, at both step counts),
so the change is comms-neutral as designed. Its payoff was never np8 wall; it was the ~16 GiB that
unblocked np32 (ledger 31).

**A REASONING-ERROR CLASS, recorded because we made it twice.** Ledger (15) "I6 plan caching REFUTED -
the plan walks are free" was correct at 600-2400 boxes and wrong as a SCALING conclusion: those walks
are `do k = 1, amr_num_blocks` per family per stage. (34)'s "performance program complete" is the same
error at larger scale. **Rule: a measurement at the current operating point can retire a COST, never an
ASYMPTOTE.** Asymptotic claims need a law measured over >=3 points at fixed per-rank work, or an
analytical argument with a measured coefficient. Corollary adopted this session: no agent-sourced or
tool-sourced specific enters this ledger, a commit message, or a decision without independent
verification of that exact claim - two of six audit agents today returned confident load-bearing
claims that were false, and both were caught only by spot-checking.

**Verification limit on this machine.** S0 sweeps ranks at fixed per-rank work on ONE node: it exposes
O(boxes) and O(P) terms but cannot exercise inter-node collectives, network topology, or the tag
ceiling. W4's law is measurable here; W5's ceiling is not reachable here at all. Decide per invariant
what "done" means before starting it.

## 2026-08-27 (34) — MERGE AUDIT + STRATEGY SHIFT: THE PERFORMANCE PROGRAM IS COMPLETE; WHAT REMAINS IS CORRECTNESS AND MERGE HYGIENE

> **STRATEGY SECTION RETRACTED by (35).** The merge-audit and coverage findings below stand. The
> strategy did not: it declared the performance program complete on a *benchmark-parity* framing.
> The goal is the END-STATE ARCHITECTURE and exascale scaling, on which the matched-point ladder is
> a waypoint. See (35) for the invariant scorecard and the corrected work order.

**Strategy.** The scaling ladder is finished. Six rungs: 1.594 -> 1.544 -> 1.368 -> 1.343 ->
1.241 -> 1.274. Five sit at or below the AMReX bar, and the third doubling (np16->np32, 1.274 vs
AMReX's own 1.278) degrades slightly LESS than AMReX. Exactly one rung is above a bar: np8->np16
at 1.241 vs 1.192, a 4.1% excess, addressed only by the optional rf:wait overlap increment
(option B, `notes/rfwait_overlap_design.md`). That excess is smaller than the ~5% single-run
noise floor these runs have historically shown; the rungs use a tighter differenced within-job
pairwise design so it is not pure noise, but it is close enough to the floor that chasing it
spends effort without producing evidence. **Recommendation adopted: declare the performance
program done and spend the remaining effort on correctness and merge hygiene.** The claim to
carry into review is "at SOTA on the third doubling, within 4% on the second."

**Merge audit of `5312e834` (base `0c9a1d43`, master `d74cc378`).** Structure is sound: every
upstream commit is a true ancestor of HEAD, all 8 master-deleted files are gone, all master-added
files present. Of 33 files both sides modified, an automated two-way scan (master deletions that
reappear at HEAD / master additions missing at HEAD) flagged 14; three clusters were then audited
individually and cleared. The Riemann 4-way split is intact (all six modules, `hypo_hlld` differs
from master by one line, every hat_R buffer allocated under renamed `_alloc` bounds, no duplicated
solver bodies); `model_eqns_4eq` removal fully honored and the symbol appears nowhere in src/ or
toolchain/; `qbmm_idx`, `m_variables_conversion`, IB neighborhood (`f5d2b4cc` + `b7d78838` both
intact), `module_categories.json` (a 4->2 space reindent, all 83 modules present), and the
HardcodedIC `pRef` casing fix all survive. Scan false-positive modes, for future reuse: line-level
matching against code master RELOCATED, and a whole-file reindent making an unchanged file look
rewritten.

**ONE CONFIRMED REGRESSION.** Master `55fb1b14` (#1717) migrated ~120 `@:PROHIBIT` lines out of
`src/simulation/m_checker.fpp` into `case_validator.py` AND added a lint rule
(`check_checker_input_constraints`) that fails precheck on any `@:PROHIBIT` in `m_checker*.fpp`
unless allowlisted or marked `! lint: runtime-check`. Our merge resolved the conflicted region as
"ours", restoring 21 of them (`m_checker.fpp:39-59, :230-231, :235-258`) — `reactive_burn` count
is 0 at master, 21 at HEAD, plus 16 in the validator — and then added `"s_check_inputs"` to
`RUNTIME_CHECKER_SUBROUTINES`, silencing the guard master built for exactly this. That exemption
sits directly beneath the INHERITED docstring stating such mixed subroutines are "deliberately NOT
listed" because allowlisting one "would have let its replacement back in unnoticed." It is a
disabled alarm, not static debt: it already admitted two new prohibits of our own
(`m_checker.fpp:92` active_box+synthetic turbulence, `:129` AMR+CBC, both from `c5b461f9`, both
input-only and both duplicated into the validator). Two of the duplicated pairs have already
drifted, in both cases with the validator STRICTER than the solver: `reactive_burn` gamma/pi_inf
uses `f_approx_equal` (|a-b|/(|a|+|b|)) in Fortran vs `math.isclose` (rel_tol) in Python, ~2x
different, so a case in the 1e-10..2e-10 relative band is rejected by `validate` yet accepted by
the solver; and `rburn%%pign` fires on the `dflt_real` sentinel in Fortran but only on an ABSENT
key in the validator, so a case writing `rburn%%pign: -1e6` passes `./mfc.sh validate` and then
aborts at runtime. `muscl_order_first_order` is imported at `m_checker.fpp:15` and used nowhere —
the fingerprint of an import hunk resolved as "ours" without rechecking uses.

**Stage guards in src/common (pre-existing, NOT merge damage).** Master carries ZERO
`MFC_PRE_PROCESS`/`MFC_SIMULATION`/`MFC_POST_PROCESS` anywhere under src/; we carry three
(`m_phase_change.fpp:99, :203`; `m_boundary_common.fpp:99`). Provenance traced: the count was 4 at
the pre-merge tip `c5abe1b9`, 4 through the merge, 3 after `c5b461f9` — our own additions, never a
merge conflict. Not a silent-else bug: the build still defines `MFC_<TARGET>`
(`cmake/MFCTargets.cmake:88-92`). They exist because `load_weight_wrt`/`sfc_partition_wrt` are
generated for simulation and post_process but NOT pre_process, which also compiles
`src/common/m_phase_change.fpp`. Fix is a choice: register those two params for pre_process (one
line in `definitions.py`, collapses two of the three) or thread arguments as
`s_initialize_phasechange_module` now does. `amr_in_fine_advance` (simulation-only,
`m_global_parameters.fpp:330`) needs the argument/policy route regardless.

**COVERAGE GAP: the post_process AMR reader has no value-level regression coverage.** Under `-a`,
post_process does run per test with `parallel_io = T`, but its output is only written to
`out_post.txt` and passed through `_process_silo_file` (`test.py:692-706`) — it is NEVER compared
against a golden; goldens hold simulation `D/*.dat` only. So the `-a` lane catches crashes and
non-zero exits and is structurally blind to wrong values. A clean pre-fix baseline over the six
np=2 AMR cases (27DEC5B6, 4644A339, 78314D65, ADA042A2, 2C46C59A, F57C3A5B) is 6/6 GREEN, which
therefore proves only that the ownership bug's spurious-abort half does not fire in those configs
— its silent halves (stream desync, displaced block coordinates) could never have been caught by
that suite at any np. **Consequence for the post-reader fix: do not gate it on the suite.** Gate
it on a direct np=2 A/B of the actual silo overlay values and block coordinates, pre- and
post-fix, and say plainly in the commit message that the reader has no value-level coverage.
Worth an upstream issue independent of this PR.

**Work order adopted.** (1) The checker/lint regression: verify per-check that the validator
covers each of the 21, delete them plus the two now-dead imports, remove `"s_check_inputs"` from
the allowlist and mark the genuinely-runtime AMR prohibits with `! lint: runtime-check`, and move
our own two input-only prohibits to validator-only — the re-armed lint rule is itself the gate,
plus AMR-68 and a bitcmp (startup-path only, expect byte-identical). (2) Audit the newer merge
`9c9e0a75` with the same method while it is fresh. (3) Post-reader ownership fix, gated as above
(patch drafted, `notes/postreader_final.md`). (4) L0 restart-writer filter (`notes/l0filter_patch.md`)
— the one item here with a real discriminator, since coexist tests under `-a` should hit post's
`lvl < 1` abort today, and a crash is precisely what that lane CAN catch. Two open decisions for
the user: whether to upstream the branch-local `damage_energy_cutoff` gate (ours since `91a6f6b4`,
living inside master's own `s_compute_hypoelastic_interface_energy`, which already damage-scales G
— a silent numerical divergence inside an upstream routine is how the NEXT merge gets resolved
wrong), and which stage-guard route to take.

**Method note.** Two `./mfc.sh test` sweeps ran concurrently on this tree for ~20 minutes because
``pgrep -fc 'mfc\|flang\|ld.lld'`` matched nothing — pgrep patterns are EREs, so `\|` is a literal
pipe — and the resulting "no processes" reading was taken as "the job died with the node." Both
sweeps shared the build/staging dir, the `tests/<UUID>/` dirs, and (via `>`) one log path; both
results were discarded and the baseline rerun alone. Same class as the HypoShearContact
cross-contamination in ledger (33). Standing rules added: never conclude a process died from a
count alone (print the matches, `pgrep -fa`), one sweep at a time per tree, and a distinct log
path per run so a second writer is self-evident.

## 2026-08-27 (33) — ELASTIC-GATE PER-SIDE FIX LANDED; ONE REVIEW PROHIBIT RETRACTED AS WRONG

The elastic-gate fix (review finding #10) landed: `s_compute_hypoelastic_interface_energy` now
gates each side's elastic energy on its own state — E_L on `G_L > verysmall` and the left damage
cutoff, E_R on `G_R > verysmall` and the right damage cutoff. The old hybrid joint-gated E_L on
BOTH sides (a damage-collapsed or fluid right state suppressed the left side's stored energy) and
gave E_R no damage cutoff at all (a collapsed right modulus kept contributing energy — the NaN
mechanism the cutoff exists to prevent). One golden moved, exactly the test that exercises the
changed regime: AMR 1D hypoelastic static block with `cont_damage` (D731AB7A), deltas
concentrated in the stress and damage fields near the damaged region (worst rel 2.3e-3 on the
stress trace; densities at O(1e-6)). Regenerated with that delta analysis as justification.

**Retraction:** the batch (32) prohibit on HLLD-hypoelasticity + CBC is REMOVED. Its gate run
blocked an UPSTREAM test (#1414's "CBC uniform preservation", 2274CB4E) whose own comment records
that a pre-fix version of the case drifted O(1e-3) in five steps — proof the CBC path is live and
consequential under HLLD, which the review finding's "corrections silently dropped" reading said
was impossible. Fresh code walk: the claim is only half-true. The conservative rows DO read the
separately-stored `flux_n` (finalized before `s_cbc` corrects `flux_rsx_vf`), but the
advection-source path consumes the s_cbc-corrected `flux_src_rsx_vf`, and #1414 deliberately
supports, fixed, and tests the combination with documented usage limits ("characteristic
boundaries are limited to fluid regions"). Outlawing an upstream-shipped configuration on a
partial static reading was wrong for this PR; the conservative-flux question is an upstream
discussion. The two prohibits covering THIS PR's own findings — AMR + CBC (the interior
fine-block-edge pollution proven by the 9640CE7F golden forensics) and active_box + synthetic
turbulence — stand untouched. LESSON (pairs with (32)'s): a prohibit that trips an existing
green test indicts the golden OR the prohibit — this time it was the prohibit, and the
discriminator was the same both times: read what the test's own history and comments claim,
then verify the mechanism end-to-end, not just the half that supports the finding.

**Suite forensics:** the 131-test gate run (hypo-65 + AMR-68) had two more failures, both false:
the HypoShearContact convergence pair (26660584 HLLD, D0AC8751 HLLC) ran concurrently under
`-j`, and their sweeps — both staged from `examples/2D_hypo_shear_contact/` — cross-contaminated:
each log held the OTHER test's Linf values (one value bit-identical across both logs was the
tell). Rerun solo: HLLD fits order 2.000 (needs 2.0±0.1), HLLC 0.998 (needs 1.0±0.2). Same class
as the known chemistry-test example-cache clobbering; convergence tests sharing an example
directory must not run concurrently — rerun solo before believing a convergence failure.

Gates: S0 np8 bitcmp vs the accepted binary **byte-identical** (pinned
`bins/simulation-elasticgate-22c01bd7`; the retraction is validator-Python only, no binary
change); hypo-65 + AMR-68 all green after triage (127 direct + D731AB7A on the regenerated
golden + 2274CB4E unblocked + the two convergence tests solo); precheck clean.

## 2026-08-27 (32) — FULL-PR REVIEW: NINE CONFIRMED FINDINGS FIXED, GATED, AND LANDED (c5b461f9)

A full adversarially-verified review of PR #1628's diff produced 12 confirmed findings; the 9
pure correctness/hardening ones landed in one batch (the elastic-gate per-side fix, the
post-process reader ownership fix, and the restart-writer level-0 filter have banked designs in
`amr-bench/notes/mrhs_merge_plan.md` and land separately with their own gates). The batch:

- **`_alloc` widening**: four dual-pass/`nc_iface` allocations in `m_riemann_solvers.fpp` and
  three UVM host-pinned fallbacks in `m_igr.fpp` were sized to `m/n/p` — under AMR's fine-grid
  swap those are ONE block's dims at allocation time, an out-of-bounds write for any later,
  larger block. Now sized to `m_alloc/n_alloc/p_alloc` (`idwbuff_alloc` where the buffered range
  is indexed).
- **Phase-change stage independence**: `s_initialize_phasechange_module` takes its allocation
  upper bound as an explicit argument (sim passes the `_alloc` maxima when load-balance/SFC
  output can reshape the grid, `[-1,-1,-1]` otherwise), removing a stage `#ifdef` from
  `src/common/`.
- **Moving-IB `num_gps` device refresh**; **coexist-safe acoustic overlap check** (skip
  `amr_block_level(kb) == 0` tiles); **64-bit regrid clustering arithmetic** (volume accumulator
  and two `jrem` linearized-index sites).
- **Three prohibits** (Fortran checker + `case_validator.py` mirrors): AMR + CBC (`bc = -5..-12`)
  — `s_cbc` is gated only on the GLOBAL bc codes with no `amr_in_fine_advance` guard, so the
  fine-block advance fires it at interior block edges; `active_box` + `synthetic_turbulence`;
  HLLD-hypoelasticity + CBC (the dual pass reads the pre-CBC Riemann fluxes).

The AMR+CBC prohibit caught one of our own tests: `AMR -> 2D -> axisymmetric` (9640CE7F) used
`bc_y%%end = -6` while its header claims machine-zero conservation "on a closed axisymmetric box"
— the nonreflecting buffer contradicted the test's own documented intent. Outer BC switched to
reflective; the regenerated golden CONFIRMS the diagnosis: every worst-case delta clusters at one
location (the static fine block's upper y-edge), where the old golden held density 0.9947 in a
should-be-uniform field, spurious radial momentum -7.1e-3 (now 6.7e-4), and an unphysical volume
fraction 1.00135 > 1 (now 0.99998) — the interior-edge CBC pollution, removed.

Gates: build clean; S0 np8 bitcmp vs `reglocal2-66ee40c7` BYTE-IDENTICAL (the batch touches no
golden-path code); AMR-68 = 68/68 after the axisymmetric BC fix; precheck clean. Pinned binary:
`bins/simulation-reviewfix-67bf6bf3`.

## 2026-08-27 (31) — SIXTH RUNG: np16->np32 = 1.274x vs THE AMReX np32 BAR 1.278x — AT SOTA; np32 UNBLOCKED

Acceptance job 386892 (reglocal binary via its byte-identical pre-format twin, k004-[002-003,
007-008]): ALL FOUR ARMS rc=0. np16 differenced (1385.551-138.608)/100 = **12.47 s/step**;
np32 (1752.075-163.192)/100 = **15.89** -> rung = **1.274x** vs AMReX's own np16->np32
degradation **1.278x** (job 385522) - **excess 0.997x, AT the bar on the third doubling.**
LADDER: 1.594 -> 1.544 -> 1.368 -> 1.343 -> 1.241 (np8->np16, excess 1.041) -> 1.274
(np16->np32, excess 0.997). Caveat: 140-step arms ((140-40)/100), not the 240-step
ledger-grade form; within-job pairwise design as always; cadence tags match the 386617
baseline runs EXACTLY (55,650,766 np16/140; 111,298,856 np32/140; escaped 0) so the physics
trajectory is identical. VRAM PLATEAU CONFIRMED: np32/140 ends at mean 37-44 / max 44.7-47.8
GiB and FLATTENS over its final 20 min (np16-plateau shape) vs the baseline's 63.9/64 crash
(rc=137). ~16 GiB reclaimed - far above the ~3.5 GB static register estimate, which
retro-explains the baseline's discrete +12 GiB jump events: they were REG_GROW device-staged
capacity doublings toward the GLOBAL count, a growth class the participation-local map
deletes outright. The np32 memory blocker is CLOSED; the residual-term worry (ledger 29) is
retired. Cross-node wall observations (np16/140 -16.6%, np32/40 -11.7% vs baseline) are
DIRECTIONAL (different node set); the same-node np8 A/B (reglocal 7.84 s/step, mergereplay
arm in flight on k004-005) is the honest neutrality read. Remaining above a bar: ONLY the
second doubling (1.241 vs 1.192); the rf:wait overlap increment (option B,
notes/rfwait_overlap_design.md) is the one candidate left, and it is optional.

## 2026-08-27 (30) — PARTICIPATION-LOCAL REGISTERS LANDED; THE ONE RED CI TEST WAS A RUNNER SLOT CAP

Register locality landed `0acad7bd` (all gates: S0 np8 bitcmp byte-identical; [amr-xa]
families identical; debug/NaN-poison rc=0; AMR-68 68/68; formatted-tree rebuild twin-bitcmp'd
byte-identical, so the queued acceptance evidence transfers). The 12 creg/freg arrays now
index by a dense participation-local slot (own / parent-owned / reflux-face-participate) with
a mesh-epoch-keyed lazy rebuild (s_amr_reg_prepare) and dense kernel sweeps; the mesh epoch
moved to m_global_parameters (m_amr re-exports). DESIGN CATCH during the site audit: both
freg exchange paths gate receives on f_amr_reflux_participates' UNCLIPPED formula (no seam
clip), so the map's clause (c) copies that formula - the seam-clipped flags predicate would
have under-covered a seam-only participant (an IRECV into dense slot 0). ACCEPTANCE = job
386892 (np32_accept.sbatch, np32/140 must plateau like np16), queued; pre-fix baseline:
np32/140 rc=137 at 63.9/64 GiB (386617, trajectory in notes/reg_locality_design.md).
np8 wall-neutrality pair (reglocal vs mergereplay, k004-005) in flight.

CI: #1628's first cycle reduced to ONE failing test - D127EC91 (AMR churn np=4, the suite's
only np=4 case; upstream has none) - ROOT-CAUSED as an Open MPI slot cap on 2-core GitHub
runners: syscheck's `mpirun -np 4` is refused before MFC runs ("not enough slots"), on every
GNU/Intel/NVHPC-cpu lane, both opt levels. NOT a code bug: the test passes on the local
gfortran CPU lane. Fix `60b931f2`: the default template passes --oversubscribe when the
launcher is Open MPI (MPICH has no slot cap and no such flag), validated locally through the
mpirun path. Teardown fix f81239a3 confirmed clearing the other 4 CPU crashes in CI.
Build & Verify green on #1628 (the lychee exclusion held). Frontier lanes still the site
outage. Upstream moved (#1763 Phoenix CI infra; merged into up/mega and the probe by the
user mid-evening - both my commits rebased on top cleanly).

## 2026-08-27 (29) — MASTER MERGED INTO up/mega: PR #1628 IS CI-LIVE

Merge `5312e834` (upstream master `d74cc378`, unchanged since the probe) reuses the
resolution proven on CI probe PR #1765: flat standard sweep keeps our rs arrays/AMR
capture/active-box; master's dual-pass HLLD dispatches from m_rhs on its own vf fields;
master's 4-way hypoelastic split is the only implementation; s_check_inputs_time_stepping
follows master's deletion. Beyond the probe: the 2a prim-batch m_rhs guards re-applied; the
hllc `extraOmpArgs map(alloc:)` workaround re-grafted onto BOTH dual-emission call sites
(the probe had silently dropped it); D731AB7A golden regenerated for master's #1414
numerics; `s_check_inputs` added to the runtime-checker lint list (validator migration =
standing merge debt). Two root-caused fixes ride along: (a) the probe's --no-verify mystery
— master's check_coverage_map_health.py ran git with the hook-inherited
GIT_DIR/GIT_INDEX_FILE, so under the pre-commit hook its unit tests queried the real repo
instead of their throwaway repos; a GIT_* scrub in the script fixes it (hook-env lint
470/470; candidate upstream PR); (b) the docs Build & Verify lane — all 41 lychee errors
are Doxygen auto-linking backticked `*.md` prose mentions to `*_8md.html` file pages it
never generates; one `.lychee.toml` exclusion (`1f38c02f`). Merge-mechanics lesson: git's
auto-merge of m_riemann_solvers/hypo_hlld silently took master's versions, which do not
compile against our one-arg s_initialize_riemann_solver — the probe's hand-edits to those
NON-conflicted files were required; a probe resolution is more than its conflict set, so
the gate was a whole-tree diff against the probe (result: working tree == probe content +
post-fork commits exactly, per file). Gates: precheck standalone AND hook-env; amdflang GPU
build; S0 np8 bitcmp byte-identical vs cdefer2-591fe7dc — and the lustre_5.dat hash equals
the probe-era bitcmp's, one identical output across cdefer2, the probe merge, and the
replayed merge; AMR-68 68/68. Binary pinned bins/simulation-mergereplay-71cd6f77.
VRAM probe 386617 (same day): np16/140 plateaus (mean ~46-51 GiB/GCD, max 57.7, flat 20
min); np32/140 climbs through that plateau in discrete regrid-event jumps to 63.9/64 GiB —
the O(global boxes) growth signature; analysis appended to notes/reg_locality_design.md,
whose acceptance test (np32/140 plateauing like np16) is unchanged.

## 2026-08-26 (28) — FIFTH RUNG (cdefer): 1.241x — WITHIN 4% OF THE AMReX BAR

Job 386573 (cdefer2-591fe7dc = 44edcc5b's code, node set k004-[001,008]): np8 differenced
(1979.970-119.175)/200 = **9.30 s/step**, np16 (2471.288-161.901)/200 = **11.55** -> rung =
**1.241x** (ladder 1.594 -> 1.544 -> 1.368 -> 1.343 -> 1.241; bar 1.192x; excess **1.041x**).
The coarse deferral did what the (27) design predicted: np16 -8.5% at flat np8; the coarse
phase's growth fell +0.76 -> +0.59 and its imbalance 1.58 -> 1.29 (absolute np16 cost 1.32 ->
1.13 s/step). Residual growth 2.24 s/step and fully DIFFUSE: coarse +0.59, restr +0.43
(rs:rest +0.22, rs:wave +0.20), rf:wait +0.33, regrid +0.32, gather +0.27, halo +0.16 - no
family holds more than 27%. The remaining ~1 s/step of wait/comm growth is the
overlap/skew-class program (options drafted in the bench notes); with it and the
register-locality np32 unblock, bar parity at 2 and 4 nodes is the near-term close-out.

Same day, the CI-exposure lane (probe PR #1765) caught and fixed three latent defects:
allocated() on a pointer member (c37dea4c), a stale AMR-hypo golden vs upstream #1414, and
grouped wave-scratch deallocates crashing every non-qbmm np>1 CPU teardown (f81239a3;
amdflang silently tolerated all three). The master merge is fully resolved, compile-verified,
byte-identical on S0, and exonerated of all CI failures - the replay onto up/mega is next.

## 2026-08-26 (27) — 2a PRICED AND GATED OFF; the coarse call deferred; np32 is MEMORY-BLOCKED

**Phase 2a's mandate was pricing, and it priced decisively.** Landed 3ed9f573 (batched cons->prim,
one launch/stage, all gates green: bitcmp byte-identical, per-rank convert launches 387 -> 30 in
the 5-step probe, AMR-68), then the pricing pair caught a ~4x WALL regression (40-arm 458.9 s vs
~118 baseline) while total DEVICE kernel time IMPROVED (8.79 -> 6.65 s/5-step; batch kernel 0.5 s).
Root-caused by discrimination probe: a flag-off build with the kernel still compiled in prices at
baseline (112.2 s) -> NOT the amdflang Attributor/codegen class; the cost is the RUNTIME host
launch/mapping burden of the per-block PRIM BRIDGE-LOADS (both the scalar_field-array-dummy form
and the per-var plain-array form are catastrophic). Those bridge loads are exactly what
store-native 2b deletes — **partial batching through a bridge is NEGATIVE value on this backend;
the finding is 2b's strongest argument, per the pre-registered decision rule.** Verdict landed
70fea5e6: gate defaults OFF, machinery kept for the 2b experiment.

**Coarse deferral landed 44edcc5b** (from the (26) wait-cluster program): the coexist L0 coarse
RHS moves past the fine advance — its rhs values are discarded and its only external products are
the creg captures phase 4 reads + a self-contained prim ghost fill (m_amr.fpp:338-341: the fine
fill prolongs from CONS ghosts via its OWN exchange), so it now absorbs stage skew where ranks are
best synchronized (rhs imb 1.14) instead of the stage top (coarse imb 1.58 at np16). bitcmp
byte-identical; AMR-68 68/68; np16 price rides the next rung. The OTHER coarse fix class
(shell-restricted L0 RHS) is DEAD BY ARITHMETIC at cap 64: a block's stencil-widened surface
bands are ~84% of its volume, so the sweep cannot shrink meaningfully — the np8 coarse base is
not reducible this way at this operating point.

**np32 is memory-blocked, not harness-blocked**: both the 240-step and 140-step np32 arms die
rc=137 (HSA_STATUS_ERROR_OUT_OF_RESOURCES, device OOM) between step 40 and 140, while per-rank
live blocks match np16 (~75-84) — something in the allocation path grows with GLOBAL box count
past a GCD at np32. The endstate's O(boxes)/O(boxes x ranks) allocation class (Phase-3 I7/S4)
is now a measured 4-node blocker, not just an audit finding. np16 probe arm: 13.7 s/step
(140-40 differenced, probe grade). The np32 ladder point waits on that fix.

Toolchain lesson recorded: mfc.sh prints "Terminated" + exit 143 on BUILD FAILURE — three
"mystery SIGTERMs" were one link error hidden by tail windows; always capture the full log.

## 2026-08-25 (26) — FOURTH RUNG (f7wave): 1.343x; THE CHAIN IS DEAD; residual growth is DIFFUSE wait

Job 385918, node set k004-[005-006], f7wave binary (bc91de59): np8 differenced
(1997.589-118.606)/200 = **9.40 s/step**, np16 (2680.162-156.279)/200 = **12.62** ->
rung = **1.343x** (ladder 1.594 -> 1.544 -> 1.368 -> 1.343; bar 1.192x; excess
**1.127x**). **The F7 wave verdict: restr np8 1.70 -> 0.53 (-69%), np16 3.74 -> 0.98
(-74%); rs:rest growth +1.69 -> +0.18 — the per-global-box chain is dead.** Same-day
k004-001 pair (lpair-k1-f7wave): reps {10.40, 9.48} vs f5asel {10.41, 9.92} — np8
~-4.5% rep-matched there; the rung-to-rung np8 delta (11.67 -> 9.40 on k004-005) mixes
the F7 gain (restr -1.18 + reflux -0.48 by phase) with day/conditions — within-job
ratios remain the honest ladder.

Residual growth (np16-np8 = 3.22 s/step): coarse +0.76 (ratio 2.35, imb 1.58, #1 as
the coarse draft predicted), rf:wait +0.65 (imb 1.51), restr +0.46 (of which rs:wave
+0.27 at imb 2.02 — the freg-wave rendezvous absorbing advance skew; rs:rest +0.18;
rs:rfp flat), gather +0.33, regrid +0.32, seam +0.24, rhs +0.24 (ratio 1.05 — physics
at the bar), halo +0.19. **No single dominator: the wait cluster (coarse-wait +
rf:wait + rs:wave ~ 1.7 of 3.22) is the top pool — the stage-skew/overlap class
re-emerges as the #1 design now that the chains are gone**, with the coarse increment
(amr-bench/notes/coarse_design_draft.md: shell-restricted L0 RHS for the np8 base,
overlap/reorder for the np16 wait) as its first concrete piece. np32 rung 385968
(f7wave, 4-node, mirrors the AMReX (4,4,2) bar 1.278x) queued the same night.

## 2026-08-25 (25) — THIRD RUNG (f5asel): 1.368x; restr growth is the F7 PER-BOX CHAIN, not skew

Job 385698, same pairwise design, node set k004-[005,008] (NOT the (24) set — cross-rung
phase deltas are directional only; the within-job ratio is the honest number): np8
differenced (2452.939-119.330)/200 = **11.67 s/step**, np16 (3355.078-163.324)/200 =
**15.96 s/step** -> rung = **1.368x** (ladder: 1.594x -> 1.544x -> 1.368x; bar 1.192x;
excess **1.147x**, was 1.30x). F5a moved np16 -10.4% while np8 stayed flat (+1.1%,
noise) — the wire-cuts-pay-more-inter-node rule now proven twice.

**The rs:\* verdict (the (24) design gate) REDIRECTS the program.** Differenced growth
np16-np8: restr +2.04 of which **rs:rest +1.69 (ratio 2.20, imbalance 1.13)**,
rs:wave +0.35 (imb 1.82), rs:rfp +0.00; coarse +1.02 (imb 1.76); rf:wait +0.31;
regrid +0.26; gather +0.24; halo +0.17; seam +0.13; rhs +0.12 (ratio 1.03 — physics
scales BELOW the bar). Closure exact (restr = rs:rest+rs:wave+rs:rfp = 3.742).
Low imbalance right after the synchronizing freg wave means rs:rest growth is genuine
per-rank cost, not skew — and per-rank restrict WORK is weak-scaling-constant, so the
cost that doubles is the **fold loop's per-global-box serialized P2P chain**
(s_restrict_fine_to_coarse: owner ISEND+WAITALL per box, receivers blocking RECV per
box, every rank walking ALL global slots in order — the per-box-rendezvous root cause
surviving in the F7 family). The assumed ~4.7 s/step "stage-skew pool" was really
~0.65 of wait (rs:wave + rf:wait) + a chain that the wave discipline already knows how
to kill. **Stage-skew/overlap mega-design DEMOTED; next increment = I5b F7 restrict
wave (design ready: amr-bench/notes/i5b_f7_design.md, its ladder gate now satisfied).**
coarse (+1.02, imb 1.76, wait INSIDE the L0 s_compute_rhs) is #2 and stays skew-class.

Same-day k004-001 anchor (user-held node, lpair-k1floor-f5asel-0825_1808, f5asel
binary): reps 9.92 / 10.90 -> **floor 10.41 s/step**, rep spread 9.9% (2x the historic
noise floor; both reps ran beside the rung job — cross-job lustre/fabric contention.
A/B probes on this node must run arms back-to-back). Node-era spread is large: same
binary measures np8 = 9.92-10.90 (k004-001) vs 11.67 (k004-005) — same-node A/B only.

## 2026-08-25 (24) — SECOND INTER-NODE RUNG (f2clip): 1.544x per doubling; the clip pays MORE inter-node

Job 385612, same self-contained pairwise design as (21), node set k004-[006-007]:
np8 differenced (2428.241-120.665)/200 = **11.54 s/step**, np16 weak-scaled
(3724.576-162.293)/200 = **17.81 s/step** -> rung = **1.544x** (was 1.594x on the
clip16k binary; bar 1.192x; excess 1.30x, was 1.34x). The F2 clip cut np8 -6.1% but
np16 **-9.0%** (19.57 -> 17.81): halved parent-fill wire buys more where inter-node
bandwidth is the scarce resource — every remaining wire cut is worth MORE at scale
than its np8 price suggests. Growth decomposition (np16-np8, s/step, vs the (21)
table): restr +2.41 (was +2.84, still #1), coarse +1.30 (was +1.53), rf:wait +0.95
(unchanged — F5 untouched in this binary), regrid +0.46, gather +0.32 (was +0.44),
seam +0.26, halo +0.25; rhs ratio 1.063 (physics is BELOW the bar). The three
wait/skew families (restr + coarse + rf:wait) are 4.7 of the 6.27 s/step growth ->
the stage-skew/overlap increment stays the #1 major target; the rs:* sub-brackets
riding the f5-complete rung (job 385698, SAME node set, queued behind this one) will
decompose restr into wave/restrict/reflux-to-parent at np16 before that design is
committed.

## 2026-08-25 (23) — F5a FACE-SELECTIVE MULTICAST: F5 halves again (-49.5%), day total F5 -58%

The L1 reflux-faces wave multicast ALL SIX full faces of every rank-boundary block to
every participating rank, but a participant applies only the faces whose outside
coarse layer it owns minus fine-fine seams (own_lo/own_hi + f_amr_face_is_seam in
s_amr_reflux_face_flags) - typically 1-2 of 6. New s_amr_reflux_faces_for(r) mirrors
the apply's gates term for term, parameterized by rank from the replicated
decomposition; the owner ships exactly each participant's apply set and each
participant posts exactly its own (identical derivation both sides, fixed order, so
the shared-tag pairing is exact). The device->host pull now covers only the UNION of
shipped faces (previously EVERY owned L1 block pulled all six faces every stage,
participants or not); the consume pushes only received faces; unreceived faces are
NaN-flooded in debug. GATES: byte-identical (bitcmp), F5 547.4M -> 276.6M words
(-49.5%, msgs 5568 -> 2814), poison probe rc=0, AMR-75. Day total for F5 (seam clip +
face-selectivity): 659.4M -> 276.6M = -58%. Third rung (f5-complete binary) to queue
behind 385612 for the pairwise ladder: clip16k -> f2clip -> f5-complete.

## 2026-08-25 (22) — F5b SEAM-CLIPPED + restr SUB-BRACKETS; k004-004 ERA FLOOR = 10.48

The freg wave (F5b) shipped ALL SIX faces of every split level>=2 child to its
parent's owner, but the parent-side apply (s_amr_reflux_to_parent) multiplies
sibling-seam faces by weight 0 — dead wire. The weight computation is extracted into
s_amr_sibling_face_weights (replicated metadata only) and BOTH wave sides now derive
the identical skip from it: skipped faces post no send/recv, skip their PCIe
pulls/pushes, and are NaN-flooded in debug so any hidden consumer aborts. GATES:
byte-identical (bitcmp), F5 words -17% on the 5-step probe (msgs drop with words —
whole faces), poison probe rc=0, AMR-75 75/75. The modest probe cut says F5a (the L1
face multicast) dominates F5 -> the next comm increment is F5a FACE-SELECTIVE
multicast (scoping: notes/overnight_0825_pricing.md — each participant applies only
the 1-2 of 6 faces whose outside coarse layer it owns; predicates pure+replicated).

Also in this batch: restr sub-brackets rs:wave/rs:rest/rs:rfp (the np16 rung made
restr the largest inter-node growth; instrument-before-optimize). And the k004-004
held-node era FLOOR (aggressive-progress protocol): np8 differenced 10.588/10.376 ->
**10.48 s/step** on the f2clip HEAD — prediction 10.3-10.8 CONFIRMED; trajectory
14.15 -> 11.64 -> 10.48. Profile: rhs 42%, restr 1.61 + rf:wait 1.42 = the top pool,
rg:build 0.12 (the rebuild pipeline is done), gather 0.79 post-F2-clip.

## 2026-08-25 (21) — THE FIRST INTER-NODE RUNG: 1.594x per doubling vs bar 1.192x

Job 385465 (2 nodes, clip16k binary = HEAD-before-F2-clip + pi16k), weak-scaled S0
family, all four arms rc=0, cadence clean: np8 differenced 12.28 s/step -> np16 19.58
= **1.594x per weak doubling across the node boundary**, vs the AMReX bar's measured
inter-node 1.192x (384936). Excess 1.34x — BETTER than the intra-node np4->np8 excess
(1.44x): crossing the fabric did not blow up. Full phase accounting (sums match walls):
the +7.29 s/step decomposes as restr +2.84 (2.56x — the NEW #1: the freg wave/
restriction), coarse +1.53 (2.96x, imb 1.67 — the L0 RHS's internal halo absorbs
inter-node skew; the shell-restricted-L0-RHS candidate is hereby PROMOTED), rf:wait
+0.94, regrid +0.48 (rb:gath 2.4x), gather +0.44, halo +0.31. rhs = 1.09x — the
physics weak-scales nearly perfectly. Wire: every family grows 2.3-2.5x per doubling;
F2 alone is 504 G-words at np16 (~half of ALL wire) — the just-landed F2 clip (-54%)
is not in this measurement; a rung rerun with the f2clip binary is queued (385612).

## 2026-08-25 (20) — F2 RING-CLIPPED: the LARGEST wire family cut -54%, byte-identical

The parent-fill wave (F2, 215 G-words at np8/240 — 3.3x F1 pre-clip) ships the padded
parent patch whose runtime consumer is the SAME amr_cg ghost fill the stepfill clip's
dead-byte proof covers (`amr_stepfill_ring_clip.md` is about the READ side of amr_cg,
independent of which family filled it; the wave's consume calls only
s_amr_fill_fine_ghosts_*, with prolong-feeding init/regrid gathers on their own
unclipped paths and the subcycle asserted away). Each child's one full-patch transfer
becomes its shell-slab list — derived by ONE shared function (s_amr_parent_shell) that
the send walk, recv walk, and consume all call, so the wire layout cannot drift between
sides; the pbmv contract keeps the full patch exactly as F1 does. Three bounded
kernels (pack/unpack/local-copy over a patch-local sub-box) subsume the full-patch
case, so the wave has one uniform path. The co-located parent copy is shell-only too,
under the same debug NaN-poison arm.

GATES: F2 words 1,628,915,568 -> 746,280,528 (-54.2%) with msgs unchanged (524) and
every other family byte-identical; step-5 output BYTE-IDENTICAL (bitcmp); debug poison
probe rc=0 with headers live on every transfer (int=2, so the untouched regrid F2 path
ran too); AMR-75. The seeded-bug arm is deliberately omitted: it probes two-sided
plan-derivation drift, which the single shared derivation function structurally
excludes, and the live debug headers still verify every transfer's bounds on the wire.

## 2026-08-25 (19) — REG_GROW GOES DEVICE-NATIVE (the store-grow twin)

`s_amr_reg_reserve`'s REG_GROW macro was the last growth path that round-tripped device
data over PCIe: every register doubling pulled the WHOLE array (all 12 creg/freg
arrays, ~1.4 MB/slot) to the host, restaged, and pushed back. It now mirrors
s_amr_st_reserve exactly: below a 512-slot transient threshold (byte-equivalent of the
store's 32-column cap) growth stages through a device temporary (copy/restore/zero
kernels, no PCIe); above it the old host round trip remains as the OOM-safe path. The
registers were already device-authoritative — every host consumer pulls its slot with
GPU_UPDATE(host=...) immediately before reading (audited: m_amr.fpp wave/exchange
paths, no whole-array host assignments exist) — so the host mirror coming out of a
device-path growth undefined is the store's existing contract. The 5-step probe
exercises the device branch on every rank (floor 64 -> immediate growth to the global
block count at startup). GATES: bitcmp PASS vs the clip16k binary (byte-identical),
AMR-75 75/75.

## 2026-08-25 (18) — THE OVERNIGHT PAIRS PRICED THE ATTRIBUTOR CLIFF; WORKAROUND ADOPTED IN-TREE

The three k004-008 np8 pairs (differenced 240-40, 2 reps each) priced the 2026-08-24
commits: rbatch (5ac7b6fa) 11.79 s/step — wall-neutral vs the 11.64 post-wave floor
(its win is calls/rank 759->15, which pays at higher np); devmig (62bf3e11) 16.69
(+43%!); clip HEAD (34296def) 16.37. The regression decomposes entirely into phases
the commits NEVER touched — rhs +3.05 s/step (+70%), rk +222%, restr +52% — while
every phase they DID touch improved (rg:mig -0.45, mg:unpk 74->11 ms/call, rg:build
-0.51, rb:gath -47%, pg:recv -70%). Diagnosis: the migration kernels crossed the
amdflang Attributor AAPointerInfo cap and the WHOLE image recompiled to slow ISA (the
known link-time instability, MFlowCode/MFC#1759); none of the overnight binaries
carried the pi16k flag (verified in every staging link.txt). Full table:
`amr-bench/notes/overnight_0825_pricing.md`.

ACTION: the workaround is adopted in-tree (`cmake/MFCTargets.cmake`, LLVMFlang link
line, mirror of PR #1759) so every future build is immune. Fresh-configure rebuild of
HEAD: flag on the link line, step-5 output BYTE-IDENTICAL to the unflagged clip binary
(the flag is optimization-only). Repriced pair + the repaired np16 rung (the 384728
failure was staging — ic/rung_np*/ lacked simulation.inp) submitted on the pinned
`bins/simulation-clip16k-f066397c`. AMReX 2-node bar measured meanwhile: np8->np16
weak-scaling 36.85 -> 43.93 s = 1.192x per doubling ACROSS THE NODE BOUNDARY — the bar
holds inter-node and is what the MFC np16 rung reports against.

## 2026-08-24 (17) — RING CLIP RE-LANDED ON THE WAVES: F1 wire -61%, output byte-identical

The reverted clip (proven correct, killed by the amdflang codegen bug, root-caused and
workaround-verified since) is reimplemented on the wave plan walks: after each pair's
box intersection, the slab is clipped against the patch's hollow shell (open core
[region_lo+1, region_hi-1] is provably dead — `amr_stepfill_ring_clip.md` survived the
revert and carries the proof), yielding up to 6 sub-slab transfers derived identically
on both sides from replicated metadata. The primitives (`s_amr_shell_slabs`,
`s_amr_shell_clip`, the debug NaN-poison arm, the shell-only own-box copy) are lifted
VERBATIM from the reverted implementation; pack/unpack/consume were already generic
over (bl, bh), so sub-slabs are just more transfers. The pbmv gather keeps its
full-box wire contract (qbmm+non-polytropic runs stay unclipped, as the original
deliberately did). Messages stay at the per-peer count (381); only payload drops.

GATES: F1 words 1,071,084,168 -> 416,141,172 (-61%) with msgs unchanged and every
other family byte-identical; step-5 output BYTE-IDENTICAL (the dead-byte proof holds
on the waves); debug probe with the NaN-poison arm LIVE (any consumer read of an
unshipped cell aborts); seeded header arm; AMR-75; wall not in the slow codegen class
(the original blocker); the wall pair queues behind the fix/feature pairs.

## 2026-08-24 (16) — THE ADVERSARIAL REVIEW ROUND: three findings on the day's four commits, all fixed

The four rapid increments got their overdue adversarial review; it returned one
corruption-class, one crash-class, and one contract-gap finding (migration and prolong
came back clean against every attack): (1) the IB path re-breaks the merge invariant
the batched reflux apply relies on — `s_amr_expand_box_over_bodies` runs after
clustering and the follow-up merge fused only OVERLAPPING pairs, so two boxes left
with a 1-cell gap have COINCIDENT outside coarse cells -> an unsynchronized
read-modify-write inside one batched kernel (amr+ib configs only; no gate exercises
them). FIX: the IB merge now fuses pairs closer than a 2-cell gap, restoring the
separation the clusterer guarantees everywhere else. (2) the device-native grow
transiently held old + tmp = 2x the array on device, at exactly the memory high-water
mark — a measured OOM class. FIX: device-native staging only up to 32 old columns
(the startup/early events where the -57% short-run win lives); above that, the old
host round trip (device peak max(old, new)). (3) restart pushes full padded columns
after writing only interiors, and the host pad bytes are undefined since the
device-native grow. FIX: zero the host column before each restart read.

## 2026-08-24 (15) — I6 (plan caching) REFUTED BY MEASUREMENT: the plan walks are free

New gw:plan/gw:pack/gw:wait sub-brackets inside the stage-fill wave (this commit).
Probe verdict: gw:plan = 0.01 ms/call — the replicated-list walks that I6 proposed to
cache cost NOTHING; there is nothing to cache. The wave's real residue is the pack
copyouts (12.4 ms/call) and the WAITALL (15.9 ms/call), plus the other PH_GATHER
paths (the parent-fill wave shares the bracket). I6 is REMOVED from the T1 list. The
evidence now points at RING-CLIP-ON-WAVES as the next comm increment: clipping the
71%-dead stepfill bytes out of the wave slabs shrinks pack, wire, and wait together,
and the compiler blocker that reverted the original clip is workaround-verified.

## 2026-08-24 (14) — STORE GROWTH GOES DEVICE-NATIVE: -57% probe wall; the last host-staged store operation is gone

`s_amr_st_reserve` grew the store by round-tripping the ENTIRE store over PCIe — full
device->host pull, host realloc, full push, for each of up to four arrays per growth
event. The round trip's one remaining justification (carrying the migration stash's
host writes across a growth) died with (12), so growth now stages through a
DEVICE-side temporary: two on-device copies, zero PCIe. The host mirror comes out of a
growth UNDEFINED — within the existing contract (every host reader pulls its slot
first; compaction already leaves the host stale by design). The 5-step np8 probe fell
100.8 -> 43.4 s (-57%!) with BYTE-IDENTICAL output — the growth round trips were a
huge UNTIMED cost in every short run (all earlier probe comparisons stand: both arms
carried it equally). Honest expectation for the operating point: much smaller (int=20
reaches its high-water early; late growth is rare) — the win concentrates in short
runs, tests (suite -1.3 min), regrid-heavy transients, and restart spin-up.
mg:slot 3.72 -> 0.60 s on the probe (-84%). Follow-up: the registers' REG_GROW macro
keeps the same pull/push pattern (smaller arrays; same fix applies).

GATES: step-5 output BYTE-IDENTICAL (np8 probe, sha256); AMR-75 75/75. Its wall pair
joins the queue when a QOS slot frees.

## 2026-08-24 (13) — PROLONG GOES DEVICE-SIDE: the rebuild now builds every slot where the store lives

The last host stage of the rebuild: `s_prolong_one_var` and the two closure prolongs
(alphas sum-to-one, species realizability) were host loops writing the cons host mirror,
forcing a full-slot H2D push per built box at THREE call sites (the rebuild, the startup
populate, the persistent-L2 build). All three are now GPU kernels (minmod was already
GPU_ROUTINE-decorated; the shared alpha limiter switch is inlined and its helper
deleted): the gathered patch's device mirror is pushed once per dispatch (patch-sized,
8x smaller than the slot pushes it replaces; redundant-but-harmless for level>=2 where
the patch is device-produced), the slot is built in place, and every full-slot push is
deleted. With (12), the ENTIRE rebuild data path — stash, migration pack/unpack,
prolong, overlap carry-forward — runs device-side: the "realloc stages through the
HOST" structural gap vs AMReX identified in the store-fix analysis is closed for the
data plane. CPU builds compile the kernels to the identical plain loops, so CPU results
are unchanged; GPU arithmetic moves device-side (not bit-comparable to the host prolong
by construction) and gates on golden tolerance.

GATES: np8 S0 probe rc=0 with the [amr-xa] table byte-exact vs the landed baseline
(slice touches no MPI); AMR-75 75/75; wall 100.1 s (parity-class). Wall pair queued
when a QOS slot frees (the two pairs ahead of it price (11) and (12)).

## 2026-08-24 (12) — MIGRATION GOES DEVICE-SIDE: the stash chain stops staging through the host

The regrid budget's migration half (rg:mig, 0.80 s/step at np8) was host-staged end to
end: full-slot D2H + host copy + full-slot H2D per owned old block (stash creation),
serial host cast loops for the MPI pack/unpack (74 ms/block unpack), a full-slot H2D
push per received replica, and a host overlap carry-forward. All four now run where the
store is authoritative: a device cons->stor stash kernel (both full-slot transfers
deleted), device pack/unpack kernels whose copyin/copyout stage EXACTLY the packed
interior (wire layout byte-identical, so the message set and [amr-xa] F4 totals are
unchanged), and a device overlap kernel — with the per-box prolong push HOISTED ahead
of it (the prolong is still a host loop; same final device state). The stash never
touches the host mirror any more, which also retires the two grow-hazard push sites and
their footguns (s_amr_st_reserve's device->host round trip now preserves the stash by
construction).

GATES: step-5 full-state output BYTE-IDENTICAL to the batched-apply binary (np8 S0
probe); [amr-xa] byte-identical; probe rg:mig 9.66 -> 5.65 s (-41%), mg:unpk 74 -> 11
ms/call; local AMR-75 goldens 75/75 (CPU) + the GPU dynamic-regrid case direct. The
differenced wall pair is queued behind the batched-apply pair on the floor node.

FOUND EN ROUTE (new compiler trap, banked in .claude/rules/common-pitfalls.md):
amdflang SILENTLY DROPS target regions nested inside Fortran BLOCK constructs from the
device image — the host registers them and the first launch dies with
HSA_STATUS_ERROR_INVALID_SYMBOL_NAME naming the missing __omp_offloading_* symbol.
Kernels must live in ordinary (module) subroutines.

## 2026-08-24 (11) — REFLUX APPLY BATCHED: one kernel per face direction, 759 -> 15 calls/rank on the probe

The verdict below ranked the reflux APPLY loop (10.4% of the np8 step — the exchange is
the wave since I5, what remained was the per-box form's up-to-3 tiny face kernels per
block per step) as the #2 overhead. `s_amr_apply_reflux` is now BATCHED: a host
precompute walks the level-1 slots with the SAME select_slot + face-flags logic the
per-box form used, fills per-slot descriptor arrays, and ONE kernel per face direction
corrects the coarse rhs for every block — the exact mirror of the capture-side batching
(`s_amr_capture_creg_dense_batch`) that has served the creg fill since the flat store
landed. The per-box form is DELETED (single call site); the subcycle path's
`s_amr_apply_reflux_state` is untouched. Legality: block corrections are disjoint (the
merge invariant keeps blocks >= buff_size apart) and a block's x/y/z outside layers are
distinct cells, so the per-direction regrouping preserves every read-modify-write —
per-(face, eq, cell) arithmetic and child-sum order are IDENTICAL.

GATES: step-5 full-state output BYTE-IDENTICAL to the pre-batch binary (np8 S0 probe,
parallel-IO restart compared by sha256); [amr-xa] byte-identical (no MPI touched);
reflux calls/rank 759 -> 15 on the 5-step probe, phase mean 2.844 -> 2.486 s
(directional; the operating-point wall pair is queued); local AMR-75 goldens 75/75.

## 2026-08-24 (10) — THE POST-WAVE VERDICT: np8 -17.7%, top rung 1.99x -> 1.65x; rhs is now the largest phase

The differenced 240-40 measurement of the full wave stack (I2a+I3+I5, pinned binary at
this commit's parent, jobs 383882/383883): **np8 on the floor node = 11.174/12.105 ->
mean 11.64 s/step vs the pre-wave 14.15 floor (-17.7%, clearing the 7.7% spread).**
Ladder: np2 5.29 (unchanged, -1.3%), np4 7.04 (-5.6%) — gains concentrate at high rank
counts, the comm-wave signature. Doublings vs the AMReX bar (1.20x/1.15x): np2->np4
1.331x (was 1.392x), np4->np8 1.653x (was 1.993x); top-rung excess over the bar 1.73x
-> 1.44x. Caveat: np8 on k004-008, np2/np4 on k004-001 (the historical ladder's node
shape); the np8 pair is the hard claim.

The np8 240-step phase budget re-ranks the program: **rhs 37.9% — physics is the
largest share for the first time.** Overheads: regrid 17.2% (build 9.1% incl. the
rebuild gather 4.5%; migration 7.2%), reflux 10.4% (the per-box APPLY loop — the
exchange is now the wave), gather 6.5%, coarse halo 5.2%, seam 3.6% (imbalance 1.465 —
residual skew parks there), ~13% residual. Next by size: (1) regrid build+migration,
(2) reflux apply, (3) I6 plan caching for the per-stage plan rebuilds, (4) I5b/F7 is
too small to appear — deprioritized (design note kept in amr-bench/notes/).

## 2026-08-23 (9) — I5-F5 LANDED: reflux faces + split-ownership freg as waves; message count unchanged BY DESIGN

F5 was the last per-box rendezvous chain: the level-1 face exchange posted and WAITALLed
once PER BOX on both sides (`s_amr_p2p_reflux_faces`), and the level>=2 split-ownership
freg handoff did the same inside `s_amr_reflux_to_parent`. Both are now single waves:
`s_amr_reflux_faces_wave` (F5a, per stage-final step) posts every receive ZERO-COPY into
the freg register host mirrors — each box owns a register slot, so no pool or unpack
exists to aggregate into — then the owner D2H + multicast ISENDs, ONE WAITALL, receivers
push H2D. `s_amr_freg_wave` (F5b) does the cowner->powner freg pairs the same way before
the reverse apply fold. `s_amr_reflux_to_parent` takes `do_xchg` (subcycle keeps its
per-box exchange, `.true.` at the dt_sub site; the lock-step driver passes `.false.`
because the wave already ran). The MESSAGE COUNT IS UNCHANGED (6708) by design — what the
wave removes is the O(boxes) chain of separate rendezvous, the same structure I2a/I3
removed for the fills. `s_amr_reg_reserve` hoists ahead of both waves (the apply can
REALLOCATE the registers under it — the map finding).

GATES: [amr-xa] F5 payload words EXACT vs baseline (659,423,232), msgs 6708 unchanged BY
DESIGN; F1/F2/F4/F6/F7 byte-identical (F6 at its landed 378); zero-copy means headers
cannot prefix the payload, so under MFC_DEBUG each transfer gets a COMPANION 8-word
message ([site, blk, ...]; never s_xa_rec'd, keeping family words comparable) plus
per-message MPI_Get_count length asserts, clean rc=0; seeded blk-shift arm aborts at the
companion check (site 39); adversarial review of the F5+F6 pair; local AMR-75 goldens
75/75. The owner-side multicast membership is the conjunction cand(region+-1 overlap)
AND f_amr_reflux_participates — reproduced exactly, header-verified.

## 2026-08-23 (8) — I5-F6 LANDED: the seam halo's cross-rank pairs as one wave; shared seam buffers DELETED

The cross-rank branch of `s_amr_fine_fine_halo` — one blocking MPI_SENDRECV per pair in
pair-list order through two shared buffers, the chain that ABSORBED the skew the fill
waves freed (I2a probe: seam 2.34 -> 3.47 s while gather fell) — is now one aggregated
message per (peer, direction) per call: every pair contributes one send and one recv
transfer on each owner, both ranks walk the SAME replicated pair list ascending, and the
per-peer offsets agree with no metadata exchange (the property the paired SENDRECVs
relied on positionally, made explicit and header-verified). Same-rank pairs keep the
batched device kernel. Converting INSIDE the routine covers all four call sites — the
per-stage driver, both subcycle sites (shape-preserving, so the I8 scope stands), and
the L0-advance coexist site the I5 row required. The two shared seam buffers and their
tile-grow reconciliation block became dead and are DELETED (net +110 LOC). Wire stays
wp with the stp cast on unpack (F6 converts; F5 must NOT — kept separate).

GATES: [amr-xa] F6 payload words EXACT vs baseline (380,849,184), msgs 2646 -> 378
(-86%, the largest collapse of the program); F1/F2/F4/F5/F7 byte-identical; debug live
headers ([site, sending slot, (d, dlo, dhi), (cnt,0,0)]) + per-message length asserts,
clean rc=0; seeded offset-shift arm aborts at the header check (site 37); adversarial
review: all ten invariants clean, three cleanliness fixes applied (orphaned comment,
stale sizing note, the pre-existing dead ym compute); local AMR-75 goldens. F5 (faces +
freg) remains per-box - next: faces need debug-only companion header messages to keep
the zero-copy recv-into-register design, and s_amr_reg_reserve must hoist ahead of any
wave that posts into freg (the apply can REALLOCATE the registers - map finding).

## 2026-08-23 (7) — I3 LANDED: level>=2 parent gathers as per-level F2 waves; per-box stage fill DELETED

`s_amr_parent_fill_wave(lev)` (m_amr.fpp) converts the per-step level>=2 F2 parent
gather — previously one pooled ISEND + one BLOCKING MPI_RECV per box, per stage, on the
majority of boxes — to one aggregated message per (parent-owner -> child-owner) pair per
level per RK stage, levels ascending (preserving the parent-before-child guarantee).
Each split child is exactly one transfer, so the plan is a pair list derived on all
ranks from replicated metadata only (f_amr_parent_block + s_amr_parent_foot +
amr_block_owner — never the lagging per-owner mirrors). Reuses the I2a wave's scratch
arrays and helpers (the two waves never overlap in time); zero new device kernels.
`s_amr_fine_stage_fill` lost its last caller and is DELETED (net +155 LOC for the
increment). Regrid keeps the chunked F2 path; subcycle keeps its per-box sites (I8);
init/static keep per-box. Two defects found and fixed before any gate ran: (1) restart
left `amr_num_levels` at its default 1 until the first regrid — the old per-box loop
keyed on per-slot levels and was immune, the new per-level driver would have silently
skipped every level>=2 fill; fixed in s_read_amr_restart (and exercised by the
multi-level-restart np=2 golden). (2) The architecture doc's timestep diagram cited the
deleted routine (review finding); updated to the wave structure.

GATES: [amr-xa] F2 payload words EXACT vs baseline (1,628,915,568), msgs 1646 -> 524
(the per-step per-box portion collapsed; remainder = regrid chunked + init), F1
unchanged from I2a (381 msgs, words exact), F4-F7 byte-identical; debug live headers on
every F2W transfer + per-message length asserts, clean rc=0; seeded one-word offset
shift aborts at the header check; adversarial review: all ten checked invariants clean,
no correctness findings; local AMR-75 goldens. Directional: the 5-step np=8 probe wall
dropped 130.6 -> 100.9 s vs the I2a-only binary (transient-dominated, n=1 — the honest
number is the differenced steady pair at the operating point).

## 2026-08-23 (6) — I2a LANDED: the level-1 stage-fill WAVE (plan-based exchange, first payoff increment)

`s_amr_stage_fill_wave` (m_amr.fpp) replaces the per-box rendezvous chain of the
non-subcycle level-1 per-stage fill — owner IRECV+WAITALL per box, contributor
pack+ISEND+flush per box, F3 blocking MPI_SEND per box — with ONE wave per RK stage:
transfer plans derived on every rank from the replicated caches (identical enumeration
order on both sides, so offsets agree with no metadata exchange), one aggregated message
per (peer, family) on the runtime tag bases (`amr_tag_base` + epoch fold), all IRECVs
posted first, device packs into contiguous host pool slices via the EXISTING per-box
kernels, all ISENDs, one WAITALL, then box-major consume through the single `amr_cg`.
Level>=2 keeps the per-box F2 path (I3); subcycle keeps its sites (I8). The driver's
parent-before-child order is strictly refined (all level-1 fills before any level>=2
gather — legal per the design doc's inset-children argument).

SCOPE DEVIATION from the I2 row, recorded: per-owned-box patch storage + `amr_cpat_off`
threading is NOT needed for the wave itself — box-major consume after the WAITALL reuses
the single patch buffer. The storage restructure only buys cross-box batched unpack
kernels (launch-count reduction) and is deferred to **I2b, contingent** on the post-I2a
phase budget showing launch/map overhead (not wait) left in the gather share. Zero new
device kernels — the amdflang codegen-lottery rule is not triggered.

GATES: (1) `[amr-xa]` exactness — F1 payload words EXACT vs baseline (1,071,084,168),
msgs 858 -> 381 (the per-step per-box messages collapsed to peer aggregates; the
remainder is the untouched regrid-chunked + init path), F2/F4/F5/F6/F7 byte-identical
in msgs AND words; (2) debug live headers on every wave transfer + per-message
MPI_Get_count length asserts, clean probe rc=0, family totals identical to production;
(3) seeded-bug arm — one-word shift of a consume offset aborted with "header mismatch:
expected site 31" (XA_F1W_SND), seed reverted; (3b) the F3 twin's wave path — never
exercised by the S0 probes (no QBMM) and goldens compile headers out — ran the np=2
QBMM-nonpolytropic case (test DDD79C8B's case) on the DEBUG binary with the [amr-xa]
report enabled: F3 38 msgs / 2,736 words send==recv exact, live headers on every
transfer, rc=0; (4) adversarial review — ONE critical
finding (the high-water grow helpers discarded already-appended transfers on growth;
invisible at S0 scale, corrupting beyond 64 transfers/rank) — fixed with copy-preserving
move_alloc BEFORE any gate ran; every other checked invariant clean; (5) goldens: the local AMR-75 subset 75/75 on k004-004 AND the FULL 708-test suite
708/708, zero failures (job 383666, mi2101x, baseline worktree at e531d354+diff) —
committed as bdb00d5d with 383666 recorded outstanding, now CLOSED; (6) wall:
2x differenced np=8 int=20 pairs with the PINNED COPY binary on k004-008 (job 383667),
priced against the SAME-NODE 3-repeat floor from job 383518: steady s/step
14.817 / 13.730 / 13.892 (mean 14.15, spread 7.7% — WIDER than the historic 5.3%; wall
claims on this case must clear it).

Ops lessons banked to memory: install-dir mtime is NOT binary freshness (a stale Aug-21
binary got probed first — pick by file mtime or explicit path); batch jobs must pin a
COPY of the binary, never a build-tree path (383518 pinned the path I overwrote mid-job;
all reps exec'd before the overwrite, so the floor is clean, but only by minutes).

## 2026-08-23 (5) — I1b-gather LANDED: identity headers live on F1/F2/F3, tripwire proven

Per the binding in amr_plan_based_exchange.md: `XA_NH` (8 under MFC_DEBUG, else 0) +
`s_xa_hdr_pack`/`s_xa_hdr_check` in m_amr_xchg_audit; every gather-trio wire site
prepends [site, blk, bl, bh] to its payload and the receiver verifies before unpacking.
Device pack/unpack kernels are UNTOUCHED (offset via argument slices) — no codegen-lottery
exposure; production arithmetic adds +0 everywhere. `[amr-xa]` records payload words only,
so its baselines stay comparable across debug/production.

GATES, all passed: (1) production build 75/75 AMR goldens; (2) debug live-header probe —
5-step np=8 S0, headers verified on all 858 F1 + 1646 F2 messages, zero mismatches,
[amr-xa] send==recv exact in every family; (3) SEEDED-BUG counterfactual — consume order
reversed locally, headers aborted with the full diagnostic (expected-vs-got slab), seed
reverted. The tiling assert (I1b's other half) needs the plan builder's periodic-wrap
analysis first and rides with I2 prep. Remaining families get headers with their
conversion increments (F5/F6 with I5, F7 with I5b). Ops note: a fresh build VARIANT
configured without FC bakes gfortran into its CMakeCache and every retry reuses it
(rc=143, looks like a SIGTERM race) — purge build/staging/<variant> and verify FC inline.

**I2 (F1+F3 plans, ~600 LOC) is now unblocked and next.**

## 2026-08-23 (4) — T1 RE-PRICED AT int=20: I2 is the program's highest-value increment

Against the int=20 np=8 steady budget, the plan-based-exchange ladder re-prices as:
**I2 (F1+F3 gather plans) targets 23% of wall scaling 3.84x/doubling — first**; I5
(reflux+seam, 20%, worst rung 6.84x) and I5b (restr, 10.8%) behind it; I3 small; I4's
"parallelise the pack" item is REFUTED (see (2) below) and its right-sizing half already
landed as I4a. Already complete from the ladder's prerequisites: I0, I1a, the mandatory
ppn=4 dynamic-regrid case, I4a, I4b-a. **NEXT SESSION: I1b-gather (headers on F1/F2/F3 +
tiling assert + seeded-bug counterfactual — the gate I2's validation requires), then I2.**
Implementation binding with line numbers: amr_plan_based_exchange.md "I1b implementation
binding". Overnight: 3-repeat np=8 int=20 pairs queued (job 383518) for the variance floor
at the new operating point.

## 2026-08-23 (3) — CADENCE MOVE VALIDATED: int=20 cuts 41% of wall; gather is the new front

**The operating point moves: amr_regrid_int 2 → 20 (amr_buf stays 4), certified and adopted**
(user approval on record). S0's feature speed ~1 at fine-CFL ~0.05 means <1 fine cell of
travel per 20-step interval, and `[amr-cad]` certifies containment at BOTH the transient
(40-step arm: 4.6M tags escaped 0) and steady state (240-step arm: 51.0M tags over 12
regrids, escaped 0). Steady box count identical to int=2 (594).

**Measurement protocol note (the transient trap):** a 40-step arm at int=20 regrids only at
steps 20/40, so the mesh is coarse for half the run — its 156.8 s wall is the
hierarchy-population transient, NOT a budget. The honest number is DIFFERENCED from-scratch
arms (240-step minus 40-step = 200 steady steps), the measure-the-step-loop instrument.
Checkpoint-restart is not used (it froze the hierarchy historically).

**np=8 steady result (k004-009, post-I4b-a binary, logs/cad20{,b}-0823_*):
14.86 s/step at int=20 vs 25.2 s/step at int=2 = −41% wall.** Cross-check: the regrid
family was ~44% of int=2 wall; amortizing ~90% of it predicts 15.2 s/step — matches.
Steady top-level budget at int=20: rhs 29.3%, gather 20.2%, reflux 14.2% (rf:p2p ≈ all of
it, rf:wait 8.5%), regrid 12.0%, restr 10.8%, seam 5.7%, coarse 3.0% (sums ~99.5%).

Consequences:
- **I4b-b and rb:slot are DEAD at the operating point** (mg:unpk 1.6%, mg:slot 0.8% of
  steady wall) — the deferral rule fired exactly as designed.
- **The ring clip's value ~doubles**: gather is now the #1 overhead and the parked clip
  removes 64–72% of gather wire words. The AMD compiler-bug report is the gate.
- **restr (10.8%) enters the target list** — the per-step restrict-to-parent chain,
  never optimized, same per-box blocking P2P shape as reflux.
- **AMReX bar protocol**: AMReX stays at ITS operating point (int=2, n_error_buf=1 at
  CFL 0.7 — an int=20 AMReX arm would need ~14-cell buffers, a different code); the
  1.20x/1.15x doubling bar is internal to each code, so it stands unchanged. MFC's
  np=2/np=4 int=20 differenced ladder gives our side (s0cad sweep).

**THE LADDER AT int=20 (s0cad-0823_1225 + cad20{,b}, all differenced 240−40, all
escaped=0, 72 boxes/rank at every np): steady s/step 5.36 (np2) / 7.46 (np4) / 14.86
(np8) → doubling ratios 1.392x (np2→4) and 1.993x (np4→8), vs 1.598x/2.63x at int=2 and
the AMReX 1.20x/1.15x bar.** Per-phase rungs (np4→np8): rhs 1.11x (AT the bar; physics
5.3 s/step and near-flat), gather 3.84x (3.01 s/step at np8 — the largest overhead),
reflux 6.84x, regrid 3.17x, seam 3.12x, restr 2.39x. Attribution complete (phase sums
match walls to 0.3%). **If the five comm families scaled at 1.2x the np4→8 rung would be
1.15x — the remaining program IS the exchange scaling: T1 waves + plan-based exchange +
the parked ring clip.** The cadence win grows with np (−15%/−26%/−41% at np2/4/8),
confirming regrid was the worst-scaling family.

Audit notes on the numbers: (a) per-step cross-check passes for gather (3.49 int=2 vs
3.01 int=20 s/step) but **reflux RISES 1.34 → 2.10 s/step (+57%, rf:wait 1.26 of it) —
frequent regrids were acting as skew synchronizers; at int=20 the skew accumulates across
20 steps and the per-step P2P families absorb it.** Already included in the −41%; it is
T1's mechanism, now visible per-step. (b) The int=2 comparator (25.2 s/step) includes its
cheap early transient, so steady int=2 is slightly higher and −41% is conservative.
(c) Cadence also relieves memory: store caps 77 at int=20 vs 125–189 at int=2 (fewer
migration waves → the replica ratchet barely engages).

## 2026-08-23 (2) — I4b PRICED AND HALF-RETIRED: store growth was the cost, not the pack

The I4a note below left I4b as "pack parallelization, blocked on the host-macro question."
The pricing measurement (new `mg:slot/pack/unpk/push` brackets splitting `rg:move`; np=8 S0
on k004-009, budgets in amr-bench/logs/i4b-price-0823_{1011,1106}) REFUTED that framing:
the pack (27.2 s) and unpack (44.3 s) are secondary — the dominant term was **mg:slot =
57.9 s mean / 90.2 s max (imbalance 1.56)**: each migration wave's replica-slot allocs grew
the shared store in 8–16-slot increments, and every growth event restages the WHOLE store
(both arrays, ~30 GiB at end-of-run caps) through the host in `s_amr_st_reserve`
(D2H + three host passes + H2D), ~12 events on the worst rank.

**I4b-a LANDED (`180cdf17`, +26 LOC, zero new device code):** `s_amr_prereserve_stash`
does one exact-target reserve per wave (at most ONE restage), and `s_amr_st_reserve` zeroes
only the NEW columns. Same-day A/B: mg:slot 57.9→14.2 s (−75%), **mg:wait 53.8→15.6 s
(−71%: the growth events were desynchronizing ranks and charging the delay to MPI waits)**,
rg:mig 206.2→126.0 s (−39%), wall 1048.9→1006.0 s (−4.1%, phase deltas far outside the 4%
noise floor). Gate: step-40 outputs BIT-IDENTICAL (coarse + 30 GB AMR hierarchy),
`[amr-scale]` trajectory identical, `[amr-cad]` tag-for-tag identical (87,993,659 /
escaped 0), 75/75 AMR goldens. Caps drop slightly (122–182 vs 125–189), live identical.

**I4b-b (pack/unpack, 71.5 s combined) is DEFERRED behind the cadence move**: those shares
are int=2 shares, and at int≈20 the whole migration family amortizes ~10x — decide at the
new operating point before building host-threading machinery. The rebuild's `rb:slot`
(25.5 s, same growth mechanism, interleaved frees make the pre-count harder) is likewise
deferred to the post-cadence budget.

## 2026-08-23 — T1/I4a LANDED: migration pools right-sized

`s_amr_regrid_stash_migrate`'s pools were sized for the worst case times every old block:
`spack/rpack(maxcnt, old_np)` (GBs per regrid at production counts, nearly all columns
unused) and `rq(old_np*num_procs)` — the O(boxes x ranks) array the endstate's W-invariants
forbid. Now sized to the blocks actually sent/received (dense column maps from a pre-pass
that reuses the send loop's own destination criterion) and the exact request count.
Message set/sizes/tags/order UNCHANGED — gated on exact `[amr-xa]` F4 equality
(39 msgs / 417,682,080 words both directions on the S0 5-step probe), identical
`[amr-cad]` counts, fast-class rhs (22.1 ms/call), and the 75-case AMR golden subset.
Remaining I4: pack parallelization (I4b) — blocked on the host-threading-macro question
(raw `!$omp` is lint-forbidden; a device pack would need the ISA-stat probe per the
codegen-instability rule) and priced by a pack-time baseline still to be measured.

## 2026-08-22 (final) — RING CLIP PARKED: an amdflang whole-image codegen regression

**VERDICT: the clip is REVERTED (this commit) and PARKED pending a toolchain fix — the
wall regression is NOT the clip's algorithm, it is the compiler.** Root cause, proven via
rocprofv3 per-dispatch ISA records (logs/rcab3-k009): amdflang (AFAR 23.2.1) generates all
device ISA in a whole-image link (`-flto-partitions`, cmake/MFCTargets.cmake), and adding
the clip's 7 target regions deterministically degraded UNTOUCHED kernels' code — weno
scratch 28→140 B (5x spill), the riemann solver's register file flipped into accumulation
AGPRs (VGPR 128→40, AccVGPR 8→136), LDS 2048→2560 image-wide — making every compute
kernel 2.4–4.5x slower (rhs 22.1→41 ms/call; np=4 wall 407→553 s). Reproduced
deterministically in a second tree; partition-count-independent (`-flto-partitions=1`
identical); constant from step one; both nodes. This is the same AGPR-blowup family as
`~/work/software/weno-agpr-repro` (which fixed the SLICE variant in 23.2.1) — a NEW
trigger in the same pipeline. Same-day diagnostics eliminated: node drift, map churn,
mapping aliasing, DPM clocks, GPU sharing, MFC_DEBUG poison, amr_shl.

**What survives the revert:** the clip's correctness is fully proven (output bit-identity
at np=4 AND np=8 incl. hierarchy files, zero transport-assert trips, wire words −64 to
−72%, gather family −33% at np=8) and the two-reviewer design (`amr_stepfill_ring_clip.md`)
is implementation-complete at commits dc6d4129+bd85c792+a7970743 (this revert's parents).
**Return trigger:** a toolchain drop whose 3-minute 5-step probe (dirs in
logs/rcab2-0822 + rcab3-k009; HEAD rhs ≈ 22 ms/call bar) plus per-dispatch ISA stats
(scratch/VGPR/LDS unchanged for untouched kernels) come back clean — then re-land the
three commits and rerun the validation ladder (bit-identity, goldens, np=1 arm, MFC_DEBUG
poison sweep, same-day A/B vs the 1.2x AMReX bar). Measurement rule from this episode,
standing: any increment that adds/removes GPU kernels must include the ISA-stat probe —
a wall A/B alone cannot distinguish algorithm cost from codegen lottery.

## 2026-08-22 (late, superseded by the verdict above) — RING CLIP LANDED, CORRECT, AND WALL-REGRESSED (open hunt)

**The ring clip (dc6d4129 + bd85c792 + a7970743) implements `amr_stepfill_ring_clip.md` in
full and is CORRECT**: output bit-identity at np=4 AND np=8 vs the pre-clip binary
(logs/rcgate-0822_1516, including the 14/31 GB hierarchy files), zero trips of the always-on
transport/coverage/frame asserts, message set unchanged, wire words F1 −72/−68% and
F2 −68/−64% (np4/np8), gather family −33% at np=8 (185.8→124.9 s).

**BUT the clip binary carries a constant ~1.7–1.9x per-launch tax on EVERY compute kernel**
(rhs 22.1→37–42 ms/call, rk ~3x, restrict ~1.8x — same call counts, bit-identical data,
phases the clip never touches), making np=4 wall 553 s vs HEAD-same-day 407 s. Refuted by
same-day A/B (all in logs/rcab*-*): node drift (HEAD reruns at 407/22.1 on BOTH k004-001 and
k004-009), per-launch metadata-map churn (bd85c792 removed all such maps — no change),
exact-width slice aliasing (a7970743 + the tax is constant from step one), the MFC_DEBUG
poison (verified off in flags), amr_shl itself (predated by the 556 s measurement). The DPM
clock theory is REOPENED (the 2-s SCLK sampling that "refuted" it cannot see ms-scale ramp).
**A 3-minute discriminator exists**: the 5-step np=4 probe (HEAD rhs 22.1 ms/call, clip
37–42; dirs under logs/rcab2-0822 and rcab3-k009). Current step: rocprofv3 kernel-trace on
both arms to split the tax into in-kernel device time vs launch-gap time, then a targeted
fix — or, if irreducible on this stack, REVERT the clip (delivered wall is the metric; the
design, proof, and gate evidence stay banked for after the launch-path tax is fixed).
Still queued behind the verdict: goldens, np=1 bit-identity (HEAD binary prebuilt in
mfc-amr-baseline), the MFC_DEBUG poison sweep, np=8 same-day A/B.

## 2026-08-22 — STEP 2 LANDED + THE EXPERT-AUDIT RE-AIM (read this before picking new work)

**Gather-batching step 2 (chunked plan-then-execute, both families) LANDED 01cc4318 + 3de4724e
and is FULLY VALIDATED** — the complete verdict, including the correctness evidence (output
bit-identity at np=4 AND np=8 across binaries), the on-node differenced walls (np=4 -1.9%,
np=8 -0.8%), the `pg:recv`-is-dataflow attribution correction, and the k004-001-GCD6 node
confound, lives in `amr_regrid_gather_batching.md` "STEP-2 VERDICT". The load-bearing
conclusion for THIS ledger: **the parent-gather wait is a dataflow dependency created by the
rebuild itself — no further MPI protocol work on that family can pay.** The rendezvous share
of `rb:gath` was ~20 s at np=8, and step 2 already collected it.

**Expert audit (2 independent reviewers: MPI internals + AMR architecture, both verified
against source) re-aims the increment ladder:**

1. **DONE same day (7b7ae5c9, logs/p1gate-0822_1259) — THE RULE FIRED: promote clipping
   ahead of T1.** Measured at np=8: stepfill 67.65e9 words/run **71.2% dead**, rb-L1
   1.48e9 words **67.9% dead**, rb-L2 4.90e9 words live-by-construction; np=4 agrees
   (70.9%/55.6%). stepfill is 91% of all gather-family words (~590 GB of host-DRAM traffic
   per 40-step np=8 run, 71% dead). Both families clear the pre-registered 50% bar ->
   **next code increment = ring-clip the step fill** (gather only the margin + one interior
   cell per face), then coverage-clip rb-L1; T1 re-ranks after. Same sweep: first
   production run of the hcid dicts (clean), and np=8 wall 1186.4 vs 1235.6 s across
   back-to-back identical runs = 4% run-to-run variance on k004-001 - never quote
   single-run walls here. Original design note follows.
   **`[amr-cov]` dead-byte counter FIRST (~40 LOC, deterministic, piggybacks any S0 run).**
   VERIFIED in code: `s_amr_fine_stage_fill` (m_amr.fpp:4982) gathers the FULL coarse patch
   every RK stage while the consumer prolongs only the ghost shell — the interior of every
   message in the 14.3%-of-np=8-wall gather family is never read. Rebuild side: the
   carry-forward overwrite is level-1-only (m_amr_regrid.fpp:1498), so the rebuild-clip
   claim covers F1 (~5% of wall) — level-2 parent bytes stay live until detail-preserving
   L2 carry-forward exists. PRE-REGISTERED RULE: dead fraction > 50% on either family
   promotes ring/coverage clipping AHEAD of T1 (clipping est. 8-15% of wall for ~300-500
   LOC and it COMPOUNDS with T1; T1's floor is 8-12% for ~3100 LOC).
2. **G-B CLOSED same day (logs/amrexs0-0822_1330, amrex_s0.sh; Advection_AmrCore's IC
   replaced with the SAME periodic-blob field as S0): THE BAR IS 1.20x (np2->4) / 1.15x
   (np4->8) per doubling** at fixed per-rank work — and the comparison is unusually clean:
   AMReX's per-rank advance is EXACTLY flat (76.8M cells/rank/step at every np) and nearly
   equals MFC's S0 per-rank fine work (77.1M/rank). MFC's ratios: 1.598x (np2->4) and
   2.63x (np4->8) -> excess 1.33x and ~2.3x. AMReX does this carrying 27.6k boxes at np=8
   (46x MFC's block count) with regrid_int=2, no subcycle, reflux on, max_grid_size 32
   (its own practice). Caveats: one linear scalar vs 6-var multiphase (ratio is internal,
   physics mostly cancels — and heavier compute HIDES comm, so MFC's worse ratio is more
   damning, not less); single runs (ratios sit far outside the 4% node variance). Axis-2
   "done" is now defined: MFC's np-doubling must fall from 2.63x toward ~1.15-1.2x.
3. **Regrid cadence is a benchmark-definition question.** Production AMR regrids every
   ~8-16 finest steps; our int=2 operating point inflates the exact phase being optimized
   (tax 27.2x at int=2 vs 7.2x at int=20 while AMReX barely moves). Before moving the
   operating point: land the validity coupling (validator rule `amr_buf >= CFL x interval`
   per level + a regrid-time containment assert, ~30 LOC) and re-baseline box counts (TRAP:
   `amr_buf` also feeds merge topology via `thr = buff_size + 2*amr_buf`,
   m_amr_regrid.fpp:596).
4. **T1 migration waves re-ranked after (1) prices clipping.** rg:mig is 213-218 s at np=8
   — still the largest single family — but clipping may be cheaper per saved second.
5. The `amr_gcr_*` chunk machinery is scaffolding for the v2 plan-based exchange (I2/I6
   cached per-peer schedules): absorb and delete it when I2 lands — two parallel exchange
   frameworks is the D2 failure mode as code.

**Matched-tax rung at HEAD (logs/tax-0822_1154, tax_analyze.py protocol-exact,
stationarity 229.6 blocks/step in window): TAX 7.06x -> 6.81x, payoff 2.15x -> 2.23x,
excess over AMReX's 3.13x now 2.18x.** Ladder: 23.92 -> 11.03 -> 7.06 -> 6.81x. The new
rung was measured on k004-001 (the slow-GCD6 node) so it is conservative; the differenced
L2 window fell 332.0 -> 302.8 s (-8.8%) despite the node. Uniform-denominator caveat
carries: honest range ~6.8-7.8x.

Still open and unchanged: the uniform-denominator re-run (13% tax error bar), the P2 idle
re-decomposition (gates the Phase-2 contract), and I1b.

## 2026-08-21 — PHASE 0 MEASUREMENTS (endstate ladder): three verdicts in one night

**M2 mechanism split — THE IDLE IS LOCAL; P2 (batched advance) confirmed as the parity lever.**
Arm A = 200^3 at np=1 (zero MPI, same coarse cells/GCD as the matched 400^3 np=8 arm B).
The discriminator could not be cleaner: **rhs per-call is IDENTICAL with and without MPI —
17.10 ms (np=1) vs 17.48 ms (np=8)**. The advance's cost does not contain MPI; it is the local
launch path, exactly as the swap-overhead and mapped-entity measurements said. Arm A still spends
68.1% of wall in rhs and 16.7% in regrid with no MPI at all. Meanwhile regrid per-call is 5.45x
more expensive at np=8 (1737.6 -> 9472.8 ms) — that excess is the MPI/aggregation share, P3's
target. Each pillar's lever confirmed by the arm that isolates it. (Per the design: no tax is
quoted from arm A — its block structure differs.)

**CMA-off transport control — waits did NOT balloon; the skew/bandwidth mechanism is confirmed.**
Two-copy vader (sender-progress-gated) raises per-call waits only ~15-18% (gather 1.53 -> 1.77 ms,
rb:wait 583.5 -> 670.6 ms, mg:wait 3227 -> 3822 ms) and wall to 877.2 s vs the 782.9-825.5 band.
A sender-progress-dominated wait would have multiplied, not added 15%. The reflux waits moved
NOT AT ALL (rf:recv 12.56 -> 11.71 ms, rf:wait 14.72 -> 13.43 — slightly lower, i.e. noise), so
the transport term lives only in the large-payload gather/migration families and even there is
small. Production (CMA-on) waits are NOT sender-progress: mechanism (a) posting-order skew is
minor, (b) pack/arrival skew + (c) node-bandwidth dominate — the v2 design's corrected mechanism
section stands, and its payoff floor framing (8-12% for T1) is calibrated, not pessimistic.

**S0 weak-scaling harness — first data ever; the W4 violation is now a measured number.**
Fixed 200^3 cells/rank, one blob per unit cell, periodic. boxes/rank exactly flat at 72 (the
construction works). Per-rank per-regrid collective volume: **ntag 0 -> 176.9 MiB and gwin
0 -> 72.4 MiB going just np=1 -> np=2** (the global tag set is replicated to every rank; grows
~np). cost bytes 1288 -> 2568 (tracks global boxes — S2's target). Wall 215.1 -> 232.2 s = weak
efficiency **0.926 at np=2** at fixed per-rank work. S1/S3/S2 now have exact baselines and a pass
bar (flat per-rank bytes). Harness: `amr-bench/s0_sweep.sh` + `s0_report.py`.

**np=4 arm DIAGNOSED + the v2 sweep MEASURED IT (2026-08-21): device OOM from a NEW weak-scaling
violation — per-GCD memory GROWS with np at fixed per-rank work.** Ranks land on distinct GCDs
(binding fine); the OOM presents as srun task Killed/Aborted, not an OOM message. The v2 sweep
(one periodic cos(pi*x)**2 IC — v1's per-arm distinct analytic ICs each forced a FULL solver
rebuild mid-sweep, since analytic ICs compile into the binary):

| np | boxes | boxes/rank | fine_work/rank | ntag MiB/rank/regrid | wall s | peak GiB/GCD |
|---|---|---|---|---|---|---|
| 1 | 72 | 72 | 76.6M | 0 | 227.8 | **49.9** |
| 2 | 144 | 72 | 76.9M (imb 1.004) | 375.4 | 263.0 (eff 0.866) | **56.4** |
| 4 | ~288 | 72 | — | — | **device OOM** | **>=60 (ceiling)** |
| 8 | — | — | — | — | **device OOM** | >=58 (uneven, died mid-spike) |

**POST-W8-FIX (2026-08-21, commit 9bcc9865, node k004-008 — the first complete weak-scaling
pair in the project's history):**

| np | boxes (last regrid) | ntag MiB/rank/regrid | gwin MiB | cost B | wall s (step loop) | peak GiB/GCD |
|---|---|---|---|---|---|---|
| 2 | 144 | 393.6 | 160.7 | 3720 | 255.6 | 55.3 flat |
| 4 | 288 | 787.3 | 320.0 | 7432 | **662.3** | 63.6 flat (hot card) |

Same build, same case family, same 40 steps at fixed 200^3 cells/rank. np=1 was not re-run on
this build; the np=8 arm has NEVER run post-fix and is the next S0 gate. Three readings:
1. **MEMORY weak scaling now holds.** Per-rank live boxes 72 at BOTH np (nboxes doubles exactly
   with np), VRAM plateaus dead flat mid-run at both, and the fix is wall-neutral at np=2
   (255.6 vs the 249-263 s band across all variants). W8: from violated to holding.
2. **TIME weak scaling is now the measured blocker: wall 2.59x per np-doubling at fixed
   per-rank work.** This was invisible before today — the np=4 arm always died. The doubling
   global entity counts are the prime suspects (regrid's O(global boxes) loops, the replicated
   collectives below, the migration storm) but the split is UNMEASURED — the phase brackets
   exist; running the np=2/np=4 pair with a phase-budget diff is the first S0 follow-up.
3. **W4/S2 baselines at np=4 are exact:** ntag and gwin and cost all double per np-doubling
   (787.3 MiB/rank/regrid of replicated tag volume at np=4 — S1's lattice tags target), final
   regrid migration 938.7 MB. Every S-track item now has a two-point scaling curve to beat.

Per-rank fine work is np-INVARIANT (1.004), so the ~6.5 GiB per np-doubling is replicated
per-GLOBAL-entity device memory. **ATTRIBUTED (code audit + VRAM-trace corroboration + exact
arithmetic): it is the STORE-CAPACITY RATCHET in its weak-scaling form.** Slot indices are
GLOBAL (f_l0_slot over the global box list); at np>=2 a rank's owned set is a shifting SFC
WINDOW of that index space, and migration allocates received non-owned slots too — so
amr_loc_n (a high-water that never decrements) grows toward the union of the windows, far
above the live count. Each capacity step costs **210.6 MiB per slot** (amr_cons_st AND
amr_stor_st at 132^3 x sys_size x 8 B each), and Fix B's compaction never fires here (gate
cap > 3*nlive = 216; cap only reaches ~107). The VRAM traces are dispositive: np=1 sits FLAT
at 49.9 GiB from startup (owned window static, ratchet inert), np=2 ratchets mid-run
50.30 -> 55.47 -> 56.33 in steps matching the A' growth trajectory exactly. At np=1 the
ratchet cannot engage at all — which is why a month of single-rank-window measurements never
saw this term. CORRECTION to the previous revision: freg/creg's np-delta is only ~90 MiB
(they were pre-over-provisioned to 128 slots by their own doubling), not "the ~1 GiB class of
the growth". Separate real find, host-side: the migration pack buffers spack/rpack are
allocated at (per-block bytes x GLOBAL old block count) = ~29 GB virtual at np=2 — the I4
right-sizing item, confirmed.

**Consequence — the W8 fix is ALREADY DESIGNED and is hereby PROMOTED: P1's device-side store
remap (kills the host round trip that forces the loose 3x/2x hysteresis, letting compaction
run at every reconcile to cap = nlive) plus local-index derivation (the index space IS the
live set; amr_loc_n stops existing). Those two items are the np>=4 unblocker and the first
S-track increment in practice.**

**W8 FIX LANDED (2026-08-21, same day): the S0 np=4 arm COMPLETES (rc=0, step loop 662.3 s,
nboxes 288 = exactly 2x np=2's 144 — per-rank boxes flat, the weak-scaling construction holds).**
What landed differs instructively from the promoted design. The full device-side remake
(gather to a staging array + realloc) was built first and REFUTED by its own gate: holding
old + staging on device is a store-sized transient, and the np=4 arm OOMed in init on exactly
that; worse, the rebuild's overlap carry-forward host-reads `amr_stor_st` WITHOUT pulling,
relying on the old growth's device->host round trip for coherence — reverting that contract
NaN'd BOTH churn goldens (27DEC5B6/D127EC91 caught it; the np>2 golden mandate paid off the
day after it landed). The fix that survived is four smaller mechanisms, each forced by a
measurement:
1. **In-place index re-densification every reconcile** (`s_amr_st_move_slot` + rewritten
   `s_amr_compact_store`): ascending-source-order device moves, no realloc, no staging. Kills
   the ratchet — np=2 AND np=4 VRAM plateau dead flat, per-rank live 72 at both (the W8
   invariant measured true for the first time). The allocation plateaus at the rebuild-transient
   high-water instead of growing run-long.
2. **Capped growth increment** (`min(oldcap/4, 16)` slots): the proportional +25% growth was
   itself store-scaled — one late-run growth event's +67-slot transient is what tipped a
   59.5 GiB card over 64.
3. **Early-free of consumed old slots** in the rebuild loop (`last_use` per old block, freed
   once the last covering new box is built): replicas held until reconcile peaked the transient
   union; this cut the np=4 plateau 57.3 -> 51.9 GiB. Freed dense indices recycle into the very
   next allocs.
4. **Stash-only replica slots** (`s_amr_alloc_slot_stash` + in-place upgrade in
   `s_amr_alloc_slot`): a received replica only ever touches its `amr_stor_st` half, but the
   full alloc gave it device-resident q_prim/rhs too (~doubling per-slot cost) across the
   np-scaled migration storm of the first dynamic regrid. This was found by the new
   `[amr-cap]` instrument (rank_time_wrt-gated live/cap line at every reconcile), which
   REFUTED the "store is the hog" theory: the store was already flat at 77 slots ~ 16 GiB;
   the other ~35 GiB was solver base + per-slot q_prim/rhs.
Validation: churn goldens 27DEC5B6 (np=2) + D127EC91 (np=4) pass, AMR subset 67/67, S0 gate
np=2 rc=0 (peak 55.3) and np=4 rc=0 (peak 63.6). **Margin note: np=4's hot card peaks at
63.6 of 64 GiB — the gate passes but with ~0.4 GiB headroom. The next memory term is per-slot
q_prim/rhs pooling (P1 proper) and stor-only stashes; do not grow this operating point without
them.** Still OPEN from the promoted design: local-index DERIVATION (deleting
amr_loc_free/amr_loc_nfree as concepts) — the re-densification makes the recycle stack
short-lived rather than gone; it remains a cleanliness increment, no longer a memory one.

**What the W8 arc changes about the ladder (2026-08-21):**
- **The weak-scaling blocker moves from MEMORY to TIME.** With np=4 completing, wall 2.59x per
  np-doubling at fixed per-rank work is now the measured gap (table above). That is the
  constitution's predicted per-global-entity anti-scaling, finally on a curve. First follow-up
  is cheap: re-run the np=2/np=4 pair with the phase-budget diff to split the 2.59x among
  regrid's O(global) loops, the replicated collectives, and migration — BEFORE building
  anything. Do not repeat the phase-share mistake of optimizing the unmeasured.
- **P1's q_prim/rhs pooling is promoted to the next MEMORY lever.** The [amr-cap] instrument
  showed the flat store is only ~16 of the ~52 GiB plateau; per-slot q_prim/rhs is a second
  store-sized term, and the np=4 hot card sits 0.4 GiB from the ceiling. Any operating-point
  growth (np=8, deeper levels, bigger blocks) needs it first. It is the same code motion P2's
  batched advance needs anyway — the two pillars now share their first concrete increment.
- **I1 (exchange validator) remains next in P3** — unchanged by today; the inventory
  corrections stand. The migration path gained two facts I1 should encode: replicas are
  stash-only slots now (an invariant the validator can assert), and the first dynamic regrid
  is a migration STORM (init ownership vs first SFC balance) — any plan-based migration wave
  design must budget for it, not for the steady-state trickle.
- **S0 harness upgrades owed:** np=8 arm post-fix (the next W8-class gate), np=1 re-run on the
  current build (table hygiene), and promote `[amr-cap]` reading into `s0_report.py` so cap
  flatness is asserted by the harness, not eyeballed from VRAM CSVs.
- **Test-capital lesson, now policy:** the churn goldens (27DEC5B6/D127EC91) caught BOTH wrong
  designs same-day — a store-transient OOM and a host-coherence NaN — and the np>2 mandate
  paid off one day after landing. Every future slot/stash/exchange change runs them FIRST
  (they are ~2 min each), before any gate sweep.

## PHASE 1 ATTACK PLAN (2026-08-21) — the operational sequence, pre-registered

Written the day the W8 gate passed, per the constitution's just-in-time rule (one phase ahead,
decision rules pre-registered BEFORE the measurements that trigger them). The architecture is
`amr_endstate.md`; this section is only the order of operations and its gates.

### M3 — attribution of the 2.59x: ALREADY DONE (the gate logs carried full phase budgets)

np=2 -> np=4 step-loop delta = +413.3 s (387.1 accounted; both arms same build 9bcc9865,
identical per-rank work, 40 steps, 20 regrids):

| family | phases | delta | share | reading |
|---|---|---|---|---|
| **wait/skew** | reflux +104.8 (imb 1.81, 0.18 -> 28.8 ms/call = pure wait), coarse +44.4 (imb 2.85), seam +8.2 | **+157** | **38%** | the "reflux is the sink" pattern is BACK at np=4 — skew originates elsewhere and drains here. Do NOT optimize reflux. |
| **rhs per-call** | rhs +110.5 — ms/call 17.29 -> 31.34 (x1.81 in the MEAN, not a straggler), imb 1.009 -> 1.465 | **+110** | **27%** | mechanism OPEN — see M4a. M2 said rhs per-call is np-invariant at the matched 400^3 point, so this is either topology or case-specific. |
| **regrid scaling** | regrid +91.8 (mig +42.7, build +37.8 incl. rb:gath +17.2, clus +5.1), gather +18.9 | **+111** | **27%** | the P3/P4 targets, now with exact per-sub-phase numbers. |

### AUDIT of the attribution (2026-08-21 evening, mid-M4) — what the existing instruments settle

Re-derived with the same-build pair (both arms of gate 1156): rhs 17.87 -> 31.34 ms/call, x1.75 —
the families stand. Deeper reads that CHANGE the interpretation:

1. **The rhs bracket is pure local advance** (`br_load + s_compute_rhs + br_store`, no MPI), and
   the phase brackets are disjoint (per-arm phase sums = 95-97% of wall), so the rhs growth is
   neither nested-wait contamination nor exchange time.
2. **The `[phase-rank]` instrument (it existed all along) decomposes the family: a FIXED straggler
   plus a uniform floor.** np=4 per-rank rhs = [212, 206, **362**, 207] (gate) and
   [188, 196, **353**, 192] (M4 armA) — rank 2 both times. Reflux waiters are exactly ranks 1
   and 3 (~170-190 s each; rank 0 ~0.4 s — no cross-rank reflux pairs by topology): the
   wait/skew family is the straggler's SHADOW plus the np-scaled cross-rank pair count (reflux
   calls/rank double at 2x2). Non-straggler ranks still run ~24.4 ms/call vs np=2's 17.9 — a
   uniform +36% floor remains after the straggler.
3. **Refuted by existing data** (each by a measurement, not an argument): GCD-occupancy contention
   (M2's matched np=8 ran ALL EIGHT GCDs at 17.48 ms/call — verified from budget_full.txt);
   VRAM fullness as a monotone cause (np=2 plateaus HIGHER than np=4, 56.5 vs 51.9 GiB, and is
   fast); host-thread oversubscription (1 host thread/rank, unbound over 128 idle cores,
   observed live); work imbalance (fine_work imb 1.004, [amr-balance]); store capacity
   ([amr-cap] equal, 77 slots every rank).
4. **The rank->GCD map is NOT identity and is load-bearing:** np=2's ranks light rocm-smi cards
   2,3 — a HIP-vs-rocm-smi renumbering. Under the natural PCI/HSA permutation (0->2, 1->3,
   2->0, 3->1), straggler rank 2 sits on card0 — the GCD pinned at 63.3-63.6 GiB (98-99% full)
   for the WHOLE run in both runs — and np=2 is fast because it happens to sit on the two
   comfortable GCDs. Leading hypotheses, discriminated inside M4: (a) a near-ceiling ROCm
   allocator cliff on that one GCD (rank-asymmetric replicas/transients put it there);
   (b) a degraded physical GCD on this node (the node-identity confound). Arm B's vram.csv
   reveals the permutation directly (RVD={0,1,4,5}: whichever cards light up define the map);
   arm D relocates ranks across GCDs (straggler follows the GCD => hardware; follows the rank
   or vanishes with headroom => memory cliff).
5. **Consequence if either holds: the "rhs np-scaling family" is NOT an intrinsic scaling law**
   — it is a per-GCD asymmetry that np=4 exposes — and the honest time-weak-scaling gap
   shrinks toward the wait/skew shadow + the regrid family. P1 pooling gains a third
   justification: headroom removes near-ceiling operation entirely.
6. **Protocol lessons, standing:** benchmark harnesses must run a PINNED binary — `mfc.sh run`
   auto-rebuilds from the live tree, so code edits during a sweep leak into later arms (armA ran
   pre-I1a code, arms B onward compile the I1a counters; negligible here — integer adds — but
   the class is dangerous). And: read the instruments you already have before designing new
   runs — [phase-rank] and [amr-balance] answered M4b for free, and M3 itself cost zero runs.

### M4 — three cheap discriminators BEFORE any building (one allocation session)

- **M4a (rhs doubling — hardware or code?):** np=2 pinned to the two GCDs of ONE physical
  MI250X vs two GCDs on DIFFERENT cards (HIP_VISIBLE_DEVICES pinning, same case). If the
  same-card pair reproduces ~31 ms/call -> the growth is HBM/Infinity-Fabric contention, i.e.
  a hardware property of dense occupancy: the weak-scaling BASELINE must be defined at full-node
  occupancy and the 2.59x shrinks accordingly. If not -> code-side; profile one np=4 rank's rhs
  before touching anything. **Rule: no rhs-side work is scheduled until M4a lands.**
- **M4b (whose skew?):** per-rank fine_work + per-level block counts at np=4 (load_weight_wrt
  already on). If L2 ownership is imbalanced -> the cost-model/balancer item (S2, migration-aware
  LB) moves up. If work is flat but waits are big -> arrival skew, T1's calibrated territory
  (floor 8-12%). **Rule: balancer work is gated on M4b showing work imbalance > 1.15.**
- **M4c (table completion):** np=8 arm post-fix + np=1 on the current build. **Rule: np=8
  device-OOM promotes P1 pooling from "next memory lever" to "immediately before anything
  else"; np=8 completing parks pooling after I1.**

### M4 VERDICTS (2026-08-21 night, sweep m4-0821_1839; armA ran the committed W8 binary,
### arms B/D/np1/np8 the same + uncommitted I1a counters — physics-neutral, marked in PIN)

- **M4a: hardware is INERT.** Three GCD arrangements (0,1,2,3 / 0,1,4,5 spread-NUMA / 0,2,4,6
  one-per-MCM): walls 593/605/581 s, rhs 29.5/30.0/29.1 ms/call, the SAME rank-2 straggler
  (353/357/358 s) every time, and each rank's VRAM fingerprint (56.1 / 62.5 / 63.3-63.6 /
  58.3 GiB) reproduced on whatever GCD it landed on. No sick silicon, no NUMA effect, no MCM
  pairing. The straggler is rank-attached software state.
- **The cliff is sharp and now bounded:** rank 1 at 62.5 GiB is FAST, rank 2 at 63.0-63.6 is
  1.84x slow — the slow path engages between 97.7% and 98.4% device-full. Standing hypothesis
  (probe `m4mem.sh`, pre-registered): cumulative per-slot q_prim/rhs alloc/free churn (rank 2 =
  SFC-middle churns most) ratchets the OpenMP runtime's RETAINED device pool to the ceiling;
  under ~1 GiB free, every target region's transient allocation hits the allocator slow path.
  Zero-code knobs found IN OUR RUNTIME BINARY: LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD=0 and
  OMPX_AMD_MEMORY_MANAGER_THRESHOLD_EXP_2.
- **M4c: np=1 re-run matches history exactly** (228.5 s, rhs 18.46 ms/call, 49.9 GiB): the
  same-code curve is 18.5 -> 17.9 -> 29-31 ms/call at np=1/2/4. **np=8 FAILS (rc=143) with ALL
  EIGHT cards pinned at 62.3-63.7 GiB** — at 2x2x2 every rank is churn-heavy, so the ratchet
  that singles out rank 2 at np=4 hits everyone, exactly as churn-retention predicts. **The
  pre-registered rule fires: P1 q_prim/rhs pooling is promoted to IMMEDIATE** (it removes the
  per-slot alloc/free churn categorically, not just the footprint).
- Sequencing consequence: I1a lands first (in flight), then `m4mem.sh` (5-min mechanism proof
  from a clean tree), then P1 pooling as the next code increment. S-track/exchange fronts wait
  for the post-P1 re-baseline: if the straggler+shadow (~half the 2.59x) falls with P1, the
  remaining honest gap is the regrid family + the unexplained uniform +36% rhs floor
  (np-attached, all np=4 arms, unexplained by every hypothesis tested today — named open
  question, do not hand-wave it).

### The increment sequence (each = one landed permanent piece, one commit, goldens-first)

1. **I1 exchange validator** (P3, ~500 LOC) — starts now, gated on nothing: it is pure code,
   prerequisite to every exchange conversion, and M4 cannot change its scope. v2 contract + the
   inventory corrections (F1/F2 shared pool). New asserts earned this week: replicas are
   stash-only slots; `[amr-cap]` flatness (promote into `s0_report.py`). Verification includes
   the I0 lesson: seed a deliberate staleness bug and confirm the validator TRIPS.
2. **P1 q_prim/rhs pooling — REDESIGNED by the pre-increment audit (consumer-lifetime sweep,
   2026-08-21 night; full report: amr-bench/qprim_rhs_lifetime_audit.md).** The sweep established: fine-block
   `rhs` is single-block scratch (fused per-block compute+update, `m_amr.fpp:4910`; reflux
   never reads slot rhs on that path) but **L0 tiles need per-slot rhs** across the
   MPI-synchronized reflux point (`m_time_steppers.fpp:635-637`); `q_prim` is written only
   under `run_time_info|probe_wrt|ib|bubbles_lagrange` (`m_rhs.fpp:849`) and read only by IB
   correction + the hyper_cleaning term.
   **Before-increment re-audit (2026-08-21, this session) verified the shape in source and
   sharpened it:**
   - Every fine slot's q_prim/rhs is allocated at the SAME uniform `mbuf*` extents
     (`s_amr_alloc_slot`) — a shared scratch is a drop-in, no per-block reshaping. The module
     already carries this exact pattern: `amr_rhs_pb_f`/`amr_rhs_mv_f` and `amr_cg` are shared
     scratch by design (docstring at m_amr.fpp:5484), and `amr_cons_br` (5625) is the
     batch-folded precedent; `amr_br_batch = 1` today, so slot-shaped scratch IS batch-shaped
     until P2 widens it.
   - The tile path calls the SAME stage routines (`s_amr_select_slot` + `s_amr_fine_stage_rhs/rk`
     at 7052/7073). Caller list is COMPLETE and verified: fused fine advance 4910-4911,
     subtree advance 5165/5168/5176, tile loops 7052/7073, m_time_steppers:581. So the
     increment's shape is a thin dispatcher: stage bodies take q_prim/rhs as dummies; tile
     callers pass `amr_slots(k)%...`, fine callers pass the scratch. `s_amr_ib_correct_fine`
     (3325/3327) gets the arrays the same way.
   - Tile slots allocate through the SAME `s_amr_alloc_slot` (tile-geometry mbuf at call time,
     per the reconcile skip-comment at 5851) — the alloc site needs the tile discriminator
     (`islot <= l0_slot_off .or. .not. amr`): tiles keep per-slot rhs; fine slots allocate
     neither.
   - q_prim gate SHARPENED: the only reader of a SLOT q_prim is `s_ibm_correct_state` — so
     per-slot q_prim exists only for TILE slots under `ib` (the RK-pass read spans other
     tiles' RHS work); every other configuration shares the scratch (a gated copy-out into
     scratch nobody reads is harmless). QBMM `pb_f/mv_f/pb_stor/mv_stor` are persistent
     per-block side-state and stay per-slot.
   - The three `allocated(q_prim)` discriminators (alloc-entry upgrade check 5746, free-path
     teardown guard, reconcile stash assert 5879) move to `allocated(x_cb)` — grid arrays are
     allocated for every full slot (tile or fine) and absent on stash-only slots; the free
     path's teardown guards become per-array (q_prim/rhs/QBMM independently).
   - Verify during implementation: `s_amr_copy_fine_fields` (s==1 backup) touches only the
     store, not q_prim/rhs.
   Expected ~15 GiB/rank saved at the S0 point (72 slots x ~2x105 MiB). Same gate as before
   (goldens + subset + S0 np=4 peak >= 10 GiB lower; np=8 after). The carry-forward device
   conversion + device-side growth remain bundled. **Bug candidate RESOLVED for our path: `m_rhs.fpp:772` reads q_prim(psi)
   before the gated copy-out, but amr+hyper_cleaning is transitively unreachable (amr+mhd is
   1D-only per m_checker.fpp:136; hyper_cleaning forbids 1D per case_validator.py:1158). On
   the MONOLITHIC path the read may see stage-stale psi when none of the gate flags is set —
   an upstream accuracy question, filed as a candidate for an upstream issue, not on our
   critical path.**

**Config decision closed (P4 arm): threshold=2^20 measured IDENTICAL to =0** (rhs
[145.0, 140.0, 147.4, 149.7] vs [145.4, 146.5, 152.3, 142.8]) — hoarding of the big freed
blocks was the entire story; descriptor caching adds nothing. **=0 is the STANDING DEFAULT,
exported in amr-bench/env.sh (dated note there); every benchmark from 2026-08-21 on runs with
it, and any runtime upgrade must re-verify it.** The peak-vs-plateau distinction from the P4
VRAM trace (rank 2 still spikes to 62.4 at the migration storm — LIVE transient,
knob-independent — while performance tracks the retention plateau) is the shape memory work
is judged on from now on. **Before quoting any S-track/regrid share: re-baseline the S0
np=2/np=4 pair WITH the knob** — the 2.59x table above is the pre-knob world; the post-knob
gap (~1.59x per doubling) is the one the remaining fronts compete over.

**S0 POST-KNOB RE-BASELINE — DONE (2026-08-21, logs/s0knob-0821_2057 np=2 + m4mem-0821_1959/P1
np=4; same build 9751e479, knob on both, k004-004).** np=2 wall 253.7 s (knob inert below the
cliff: peaks 43.3/43.7 GiB; pre-knob np=2 was 255.6 on the pre-I1a binary). np=4 wall 405.4 s
= **1.598x per doubling, gap +151.7 s**, and the split is decisive:
regrid +86.6 s (57%% of the gap: rg:mig 12.4->54.1 = +41.7; rb:gath 1.7->18.1 = +16.4 at
5.5x/call; rb:tail 1.9->23.6 = +21.7 with imb 2.44; rg:clus +4.5) + exchange-wait families
+42.7 (gather +17.5 at 3.7x/call, reflux +11.0 pure wait, coarse +8.2 imb 2.17, seam +6.0)
+ rhs only +7.8 (per-call 17.64->18.63 ms = +5.6%%: the allocator fix fully retired the rhs
story). rb:slot and rb:ovl FLAT (-1.7/+0.2). **Post-P1 front CONFIRMED: the per-box gather
(rb:gath + its rb:tail shadow) and migration (rg:mig, T1 waves) own the remaining scaling
excess; every S-track share quoted from here on uses THIS pair, not the 2.59x table.**
3. **Third increment — SELECTED by the post-knob re-baseline (2026-08-21, supersedes the
   M4-directed menu below): rb:gath + rb:tail (per-box gather batching), scope SHARPENED by
   the np=8 sub-brackets: rb:gath = pg:recv 99.2 s (level>=2 parent gather's BLOCKING per-box
   MPI_RECV — the unconverted half of R1) + rb:wait 58.9 s (level-1 WAITALL); pack/unpack/
   alloc all ~0. BOTH families go through one chunked plan-then-execute (design:
   amr_regrid_gather_batching.md, updated with the S0 evidence — its old level-1-only scope
   was matched-point-specific).** Judged on rb:gath ms/call flattening across np and the
   rb:tail/reflux/seam wait shadows shrinking with it; increment 1 (plan reproduces today's
   message set, asserted via the I1a XA counters) is the safety net before any batching. **Fourth: T1 migration waves**
   (rg:mig +41.7 s, the largest single item; design v2 reviewed, floor 8-12%%). The old menu
   (S1 lattice tags / I2-I3 waves / S2 balancer) stays written below as the fallback rules if
   either increment's gate fails.
   Original menu: S1 block-lattice tags (if regrid/collective scaling is the
   chosen front — kills the measured ntag doubling; judged on ntag bytes/rank going flat), OR
   I2-I3 exchange waves (if wait/skew traces to per-box exchange arrival), OR the S2/T1
   balancer path (only if M4b's rule fires). One front at a time; the other rules stay written
   so the choice is a lookup, not a debate.

**P1 LANDED (f5f99337) AND GATED (2026-08-21 night, logs/p1gate-0821_2144, pin e964b45f):**
subset 67/67 incl. both churn goldens. **np=4 peaks 30.2-39.5 GiB vs 63.6 pre-P1 (~27 GiB
drop; gate was >=10). np=8 COMPLETED FOR THE FIRST TIME EVER: rc=0, peaks 39.0-51.3 GiB
(predicted ~48), live 72-75 boxes/rank — the MEMORY axis (W8) now holds through np=8 at
fixed 200^3/rank.** First np=8 phase budget (wall 1106.0 s = 2.63x over np=4): regrid 48.1%%
(rg:mig 204.8 s over 11 migration events vs 4 at np=4; rg:build 287.9 with rb:gath 168.5 s =
15.2%% of wall at 6104 calls, 27.6 ms/call) + gather 14.3%% (2.42 ms/call, 3.2x np=4) +
reflux 5.3%% + seam 4.5%%; rhs healthy at 15.7%% (21.3 ms/call). **The np=8 wall IS the
per-box gather + migration — increments 3 (rb:gath) and 4 (T1) confirmed with force.**
Two predictions corrected, honestly: rb:slot did NOT collapse (21.1 s at np=4, unchanged —
that bracket is dominated by host-staged store growth, not the removed per-slot allocs), and
np=4 wall came out 421.2 s vs 405.4 pre-P1 (+3.9%%) — **ADJUDICATED AS VARIANCE by a repeat
arm (410.7 s, rhs 18.88 ms/call vs 18.63 pre-P1 / 19.65 first run; VRAM peaks byte-identical
across both post-P1 runs). P1 is wall-neutral at np=4, as the mapped-entity law predicts
(same dummy count per region, different backing). Post-P1 np=4 wall band: 410.7-421.2 s.**
- **G-A: CLOSED same night (logs/tax-0821_2241, tax_matched_p1.out; the analyzer reproduces
  the historical 11.03x on the old file — protocol-exact, stationarity verified at 229.6
  blocks/step in BOTH windows). MATCHED TAX 11.03x -> 7.06x; PAYOFF vs uniform-at-finest
  1.38x -> 2.15x; excess over AMReX's 3.13x: 3.52x -> 2.26x.** Delta between the pairs =
  knob + P1 only (zero solver-algorithm change): the differenced L2 window fell 489.0 ->
  332.0 s (-32.1%%) at identical fine_work; the uniform window is stable (10.51 -> 11.16 s).
  Carry-forward caveat: the uniform denominator is still the one flagged 13%% high (open
  row: uniform re-run at cc59ad38) — honest tax range ~7.0-8.0x. P2's old 14.91x idle
  factor is now UNMEASURED at the matched point (that decomposition predates both fixes) —
  re-decompose before ranking P2 against increments 3/4.
- **G-B: "no performance tax over SOTA" has NO weak-scaling bar.** We compare AMReX at one
  operating point but judge scaling only against ourselves (1.598x/doubling). Run the
  AMReX S0-equivalent at np=2/4/8 on this node (amrex_tax.sh is most of the harness) to
  define the target our np-doubling ratio must meet. Until it exists, axis-2 "done" is
  undefined.
4. **Phase-2 (P2 batched advance) contract is WRITTEN during increment 2-3** per the
   constitution — seeded by the 2a prototype (batch convert_conservative_to_primitive) and
   M2's verdict; it does not start until its contract passes the audit ritual.

### Re-audit cadence (user directive 2026-08-21: "reaudit things as needed regularly")

An evidence audit is a ZERO-RUN activity — re-read the instruments, re-check same-build/
same-node provenance (the PIN files), re-derive the headline ratios from the raw tables — and
the 2026-08-21 audit overturned two working theories and found the straggler without a single
new run. It is now scheduled, not occasional:

- **Before starting any increment:** does the evidence still support its premise? (The device-
  remake increment failed exactly this check in hindsight — its premise ignored a measured
  transient.)
- **After every sweep:** do the new data change any pre-registered rule's input? Read
  [phase-rank]/[amr-balance]/[amr-cap]/[amr-xa] per rank, not just the means.
- **At phase boundaries:** a full re-derivation of the operating picture against the
  constitution's invariants, with anything refuted banner-marked in the docs AND memory.
- **Confound checklist at each audit:** build identity (PIN), node identity (PIN), bracket
  nesting (phase-sum vs wall), cumulative-vs-stationary phenomena (rule 16), rank->device
  permutation (rule 6).

### Standing constraints for this phase (the do-not list)

- Do not optimize reflux (it is the sink, not the source — twice-confirmed).
- Do not grow the operating point (np=8, deeper levels, bigger blocks) before P1 pooling.
- Do not re-litigate P2's position (constitution D-phase2: full commitment, after Phase 1).
- No multi-node work until single-node np=8 weak scaling is clean (D-node).
- Every slot/stash/exchange change: churn goldens FIRST, then subset, then gates; wall claims
  only from from-scratch same-build pairs; counters over wall wherever a counter exists.
- Watch CI on 9bcc9865/7ca673a0: `s_amr_st_move_slot` + stash-only slots are new device code
  paths for the other three compilers (local gate is amdflang-only by standing rule).

**Mission: drive the AMR infrastructure tax toward zero.** Physics (`rhs`, `coarse`, `rk`) is
untouchable; everything else is overhead to be removed. This version supersedes the 2026-08-18
rewrite (git history) now that the WHY is established — the findings live in
`amr_slowness_analysis.md` (causal model, five-reviewer panel) and `amr_tax_review.md` (measurement
audit); this document is only the work list, its gates, and its decision rules.

## 2026-08-20 EVENING — WHAT THE MEASUREMENTS DID TO THE PLAN

Four items closed or downgraded in one session. **Read this before acting on the ladder below.**

| item | status | evidence |
|---|---|---|
| T0a non-blocking reflux recv | **correct but 0% measured** | -2.0% wall, inside a 5.3% floor; a later same-config run came in HIGHER than baseline |
| M4 t8code guard | **no-op, already exists** | `rr /= proc_rank` at the send loop already excludes the only rank holding the block |
| T0b send aggregation | **dead** | fan-out measured **1.048** — each moved block has one destination |
| T0b SFC hysteresis | **alive but marginal** | `rg:move` is 50% wait / 50% volume; payoff 4.7-8.4% against a 5.3% floor |

**THE NOISE FLOOR IS 5.3%**, from four runs of the identical configuration (782.856 / 799.022 /
818.498 / 825.481). Every remaining T0 item is a 4-8% candidate. **T0 will not produce a defensible
multiplier**, and the earlier 11.03 -> 8.54x ladder projection is withdrawn.

### The one good methodological find

`[amr-mig] ... bytes 10748501376` was **identical to the digit** across two runs. Migration volume is
deterministic, so a hysteresis A/B can be judged on **bytes moved** — exact, zero noise — instead of
wall time. Confirm the volume cut first; spend a timed run only if it lands. Apply this pattern
wherever a fix targets a deterministic quantity: it converts a 5.3%-noise experiment into an exact one.

### Where the effort should go

Infrastructure is 78.2% of wall against a 6.6% parity budget (11.9x reduction needed). Nibbling
single-digit percentages cannot get there. **T1 (per-block -> per-level conversion) is the only item
that touches the 78.2%, and it is the same refactor as S4 (distributed box metadata).** That is the
program; T0 is not.

## TWO PROGRAMS: TAX (T) and SCALING (S)

Written 2026-08-20 after the full budget. These have **different success metrics and must not be
measured against each other** — the scaling fixes score ~0 on the matched benchmark, and the tax
fixes do nothing for O(N) growth.

| | Track T — tax | Track S — scaling |
|---|---|---|
| goal | matched tax 11.03x -> 3.13x (AMReX) | per-rank memory and collective volume O(1) or O(log P) in box count |
| metric | ns per cell-step, matched point | per-rank bytes + bytes-per-collective as f(problem size) |
| harness | exists (`tax_matched.sh`) | **does not exist** - S0 below |
| current | 11.03x, infra 78.2%% of wall | `union_gtag` is **64 GiB/rank at 4096^3** |

### THE KEY STRUCTURAL POINT: the two tracks CONVERGE at the top

T's big item (convert per-box operations to per-level plan-then-execute) and S's big item (stop
replicating global box metadata) **are the same refactor**. A per-level plan only needs *this rank's
boxes plus its neighbours*; building it is what makes the global tables unnecessary. Doing them
separately would mean writing the same code twice and merging conflicts. **Co-design T1 and S4.**

The cheap items on both tracks are genuinely independent and can proceed in any order.

### Track T — tax

| id | item | LOC | expected |
|---|---|---|---|
| T0a | non-blocking reflux receives (`rf:recv` 6.3%%) | ~30 | 11.03 -> 10.5x |
| T0b | migration volume: t8code guard + SFC hysteresis (`rg:move` 12.4%%) | ~50 | -> 9.9x |
| T0c | chunked gather batching (`rb:wait` 9.0%% + the skew sinks it feeds) | ~300 | -> 8.8x |
| T0d | flatten `q_prim`/`rhs` per-slot allocs (`rb:slot` 3.4%%) | ~150 | -> 8.5x |
| **T1** | **per-level plan-then-execute throughout** | months | **-> ~4-5x** |
| T2 | fused advance; attacks the 2.41x per-cell inefficiency | months | -> ~3.5x |

**T0 as a whole is worth ~1.3x and does NOT reach AMReX.** Infrastructure is 78.2%% of wall and the
parity budget is 6.6%% - an 11.9x reduction. T0 delivers ~1.4x of that. Do T0 because it is cheap and
de-risks measurement, but the parity claim rests entirely on T1.

### Track S — scaling

| id | item | LOC | why |
|---|---|---|---|
| S0 | **weak-scaling harness** | ~150 | Track S has no metric today. Report per-rank peak bytes and per-collective volume vs problem size. Without it S1-S4 are unfalsifiable. |
| S1 | coarsen tags to a block lattice before `union_gtag` | ~60 | **512x volume cut at 8^3.** 491 MiB/rank -> 0.96 MiB now; 64 GiB -> 0.12 GiB at 4096^3. AMReX `blocking_factor`, Uintah region lattice. |
| S2 | `MPI_SCAN` prefix weights instead of `ALLREDUCE(cost, nboxes)` | ~40 | Carries **8 B/rank regardless of box count**, O(log P). p4est's approach. |
| S3 | local clustering + boundary reconciliation (SAMRAI) | ~400 | Removes the global tag set entirely rather than shrinking it. |
| **S4** | **distributed box metadata** | months | `amr_block_owner`/`amr_slots`/`amr_region_*_all` are sized `1:amr_max_blocks` on EVERY rank, and the rebuild loops `do k = 1, nboxes` globally. O(N) memory and O(N) trip count per rank. **Co-design with T1.** |

### The three global collectives, measured

| site | carries | current | at 4096^3 |
|---|---|---|---|
| `s_amr_union_gtag` ALLGATHERV | tagged **cell** indices | 491 MiB/rank, 1.53 GiB/regrid | **64 GiB/rank** |
| `s_amr_pack_gwin_pairs` 2x ALLGATHERV | window pairs | same class | same class |
| `s_amr_block_cost` ALLREDUCE | `cost(amr_num_blocks)` | 1.8 KiB/rank | 1.8 MiB/rank |

**Warning for whoever measures S1/S2:** they cost ~0 on the matched benchmark (`rg:clus` 2.2%%,
`rg:tag` 0.5%%). Judging them by wall time at 400^3 will make them look worthless. Judge them by S0's
metrics.

## THE LEDGER — every idea, its status, and what would move it

Audited 2026-08-20. This is the index; the chronological evidence is below. **Read this section
first and do not re-open a CLOSED row without new evidence** — several have been re-proposed more
than once.

### The frame: what the tax actually decomposes into

The matched head-to-head gives **23.92x = 1.60x ARITHMETIC x 14.91x IDLE**. The arithmetic factor is
benign (MFC does real multiphysics against a linear-advection benchmark). **The excess is idle**, and
its shape is explained by launch count: ~224 blocks x 3 RK stages x ~21 kernels ~= 14,100, against a
measured 14,091/step versus AMReX's 81. AMReX launches once per *level* because a MultiFab kernel
spans every box; **we launch once per block**.

Every alive item below attacks one of five idle sources:

| # | idle source | measured share | attacked by |
|---|---|---|---|
| I1 | per-block kernel launches (not fused) | 173x launch ratio | B1, B2, B4 |
| I2 | MPI progress spinning (np=8) | ~51% of host CPU | C1, C2, C3 |
| I3 | per-region GPU map/descriptor cost | 2.00 copies per mapped array | **IRREDUCIBLE — see closed** |
| I4 | per-box gather in regrid | 16.9% of wall | R2, R3 |
| I5 | regrid frequency | 37% of wall at int=2 | R1, R4, R5 |

**Unresolved conflict that gates I1 vs I2**: `PH_RHS` is 54-57% GPU-busy at `int=20` but overall busy
is 5.4% at the matched `int=2` point. Different operating points. Until reconciled, we do not know
whether fusing the advance (I1) or cutting regrid (I5) is the larger lever. **M2 below settles it.**

### CLOSED — landed

| id | item | result |
|---|---|---|
| L1 | Level-2 parent-gather ISEND pool | -17 to -22% wall, regrid -39.3%, 76/76 |
| L2 | Loop-invariant coarse-halo hoist | ~2.4x on its path |
| L3 | **Store: bounded growth + compaction** | **2.32x** (1.80x compact x 1.29x growth) |
| L4 | Flat store is authoritative for q_cons | prerequisite for all batching |
| L5 | max_grid_size heap corruption | fixed + golden |
| L6 | Multi-level + subcycle + np>1 PROHIBIT lift | covered by C45DBB52 |

### CLOSED — refuted, with the evidence that killed each

| id | item | killed by |
|---|---|---|
| X1 | Packed super-grid / block packing | even-split tiler leaves blocks exactly slot-sized |
| X2 | Cost-weighted balancer | work balanced to 1.4%; imbalance 1.015-1.027 all run |
| X3 | Straggler = gather ownership / co-location | gather barrier 29.6 s vs 0.185 s exchange (160:1) |
| X4 | L0-sourcing of coarse patches | gather never skewed to the sick rank (argmax rank 6) |
| X5 | Churn causes the straggler | ranks 4/5: identical churn, 0.998x vs 1.089x rhs |
| X6 | "Cost grows with sim time" is intrinsic | 3.68x pre-fix vs **1.17x post-fix** |
| X7 | Batch-convert remaining blocking calls | hoisting the step gather drain made wall WORSE |
| X8 | Bigger blocks (cap > 64) | cap 96 +18%/cell; cap 128 device-OOM on per-block size |
| X9 | Reducing per-region GPU map cost | **7 mechanisms measured, ALL fail** — treat as a floor |
| X10 | Optimising reflux | 99.6% comms, and it is the SINK not the source |
| X11 | rocprof-sys-causal profiling | rejects MFC's .text FUNCs as ineligible |

### ALIVE — ranked by (value x confidence) / cost

**Tier 1 — cheap, high confidence, do regardless of pending measurements**

| id | item | LOC | why |
|---|---|---|---|
| S1 | **Device-side store remap** | ~80 | Removes the host round trip that forces compaction's 3x/2x hysteresis. Steady state 2-3x live -> **1x**. AMReX `RemakeLevel` does exactly this. Given L3 measured that *carrying* capacity costs 1.29x on its own, this is directly valuable. |
| S2 | Derive the local index, delete the recycle stack | ~40 | AMReX `localindex` = binary search over the owned-box list; Parthenon `lid = n - nbs`. Index space *is* the live set, so the ratchet cannot return. Mostly redundant with S1 but deletes the bug class. |
| A10 | Migration phase bracket | ~10 | We cannot currently price migration at all — regrid and redistribute are fused (see M3). Uintah's separate brackets are what found their scaling limiter. |
| R5 | Grid-identity early-out | exists | Already have the `same` check; verify it covers the level set too. |

**Tier 2 — the main line, gated on M2**

| id | item | LOC | why |
|---|---|---|---|
| B1 | **Device-resident descriptor array + one fused kernel** | weeks | Parthenon's `BndInfo`/SparsePack shape. *Fully compatible with our contiguous store* — each descriptor degenerates to an integer slot index. This is the direct attack on I1. |
| B2 | Rebuild-descriptors-only-on-invalidation | ~30 | Parthenon rebuilds packs only at remesh. Prevents the fused path re-acquiring per-step host cost. |
| R2 | **Batch the per-block regrid data motion** | weeks | SAMRAI fills a whole level with ONE RefineSchedule; Chombo with ONE copyTo. We do one `s_amr_gather_coarse_patch` per block (16.9% of wall). Converts churn from O(blocks) to O(1) plans. |
| R3 | Plan-then-execute for the parent gather | ~200 | Narrower version of R2 targeting the level>=2 path specifically. |

**Tier 3 — load balance; cheap but no longer aimed at a known problem**

| id | item | LOC | why / caveat |
|---|---|---|---|
| A1 | Per-block constant in `s_amr_block_cost` | ~10 | Uintah ships `patchCost = 16` cells; SAMRAI added `minimum_patch_load` for exactly our regime. **Caveat: work is already balanced to 1.4%, so this buys little today.** |
| A7 | SFC-cut hysteresis (gain threshold) | ~40 | Uintah `gainThreshold` 0.05; WarpX measured optimum 10%. Zero correctness risk. |
| A2 | Fit cost coefficients to measured per-block time | ~100 | Uintah `ModelLS`. **Heed A4**: Uintah's own dissertation found measured filters *worse* than the fitted model immediately after a regrid — and LB always follows a regrid. |
| M3 | `vsize` migration-cost companion to block cost | ~30 | Every scheme in the literature needs work AND redistribution cost; we model only work. Cheap and exact for us: `sys_size x fine cells`. |
| M4 | t8code guard: don't send a block to a rank that already has it | ~10 | `s_amr_regrid_stash_migrate` has no such guard. |
| M5 | Decouple adapt from repartition (p4est Principle 2.1) | ~200 | We fuse them, so we pay the expensive half every time we want the cheap half — and cannot price them separately. **Prerequisite for measuring migration at all.** |

**Tier 4 — worth trying, low effort, independent**

| id | item | why |
|---|---|---|
| C1 | Ranks per GCD (MPS-analogue) | Parthenon Table 1: ~1.8-2x independently, near-zero code change |
| C4 | Per-thread-block scratch instead of per-block | Parthenon-VIBE: 8.858 GB -> 0.138 GB (modeled, not end-to-end) |
| T1 | `amr_tag_eps` per-level thresholds | a single threshold populates exactly ONE level on smooth features |
| F1 | Delete the flux families | in progress; its own premise was wrong (delete, don't flatten) |

### BLOCKED — needs a measurement before it can be ranked

| id | question | experiment | status |
|---|---|---|---|
| M1 | What is the tax NOW, post-store-fix? | matched arm, differenced 40->80 | **RUNNING** |
| M2 | Is the excess launch-bound (I1) or MPI-bound (I2)? | **np=1 run of the MATCHED case** | designed, never run — the single highest-value open experiment |
| M6 | Why is rank 5 slow when rank 4 has the same capacity? | per-rank `hipMemGetInfo` at each regrid | ~10 LOC |
| M7 | Does slot size vary enough to justify per-box sizing? | max/mean box-volume ratio | derivable from a log |

### M2 DESIGN — separating local launch cost from MPI progress

The two idle findings (85% survives at np=1 with zero MPI; host 51% busy in MPI progress at np=8)
were measured on different cases and have never been reconciled. The obvious experiment — run the
matched case at np=1 — **does not work, and it is worth recording why** so it is not attempted again:

- The matched case is 400^3. At np=8 we measured 44-64 GiB in use *per GCD*. Halving the rank count
  doubles per-GCD memory, so np=4 needs ~88-128 GiB against a 64 GiB device. **np<8 on this domain
  OOMs immediately**; np=1 is off by 8x.

**The design that does work: hold per-GCD work constant instead of holding the domain constant.**

| arm | domain | ranks | cells/GCD | MPI |
|---|---|---|---|---|
| A | 200^3 | 1 | same as B | **none** |
| B | 400^3 | 8 | same as A | yes |

200^3 at np=1 puts exactly the same number of cells on one GCD as 400^3 at np=8. Compare
**ns per fine-cell-step** between the arms:

- A ~= B -> the idle is **local** (launch path / descriptor mapping). Fusing kernels (B1) is the lever.
- A << B -> the idle is **MPI progress**. The convoy work (C-tier) is the lever, and fusing the
  advance would be aimed at the wrong term.

This needs no rocprof: ns/fine-cell-step is already the reported quantity, and the phase brackets
carry `GPU_WAIT` on both ends so they measure completed device work.

**Known imperfection, stated up front:** a 200^3 domain does not produce the same AMR block
*structure* as 400^3 (fewer blocks, different boundary fraction), so this is a controlled comparison
of the idle *mechanism*, not of the tax. Do not quote a tax from arm A.

### THE STRATEGIC FACT nobody should lose sight of

**AMR payoff is currently < 1 at single-node scale — uniform refinement beats our AMR by 4-5x.** Every
item above is about making AMR *cost less*, not about whether it currently *pays*. The tax must fall
by roughly an order of magnitude before AMR is the right choice on one node; the case for it is
multi-node, where uniform refinement does not fit in memory. **Do not let a local win obscure that
the headline payoff is still negative.**


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
