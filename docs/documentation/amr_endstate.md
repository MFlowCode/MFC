# The AMR end-state: weak-scaling architecture and the re-derived program

Written 2026-08-20, from the ground up, after a month of measurement. This document is the
**constitution** for the AMR performance program: it states the architecture we are building toward,
disposes of every piece of work done so far (keep / carry / delete / never-revisit), and re-derives
the increment ladder so that **every increment is a permanent piece of the end state**. The
day-to-day work list remains `amr_action_plan.md`; the exchange-layer contract remains
`amr_plan_based_exchange.md`. Findings update this document; they do not re-found it.

## 1. The diagnosis — why a month of work has not closed the gap

The campaign produced real wins (store fix 2.32x, ISEND pool -17 to -22%, tax 23.92x -> 11.03x,
payoff now above 1) and a definitive evidence base. It also has a structural failure worth naming,
because it is the answer to "something critical is going on":

**D1 — Sequencing by phase share, not by architecture.** Each cycle attacked the biggest phase at
the current operating point. But the measured decomposition proves no single-phase attack can reach
parity: 23.92x = 1.60x arithmetic x 14.91x idle, and parity needs ~41% GPU busy while *deleting all
AMR infrastructure outright* only reaches 16.4% — the advance itself is ~22% busy from per-block
launches. Two independent overheads, **both required**. Yet "batching the advance" was deferred
because "regrid is 47.4% now" — an operating-point artifact, not a refutation. The result: a month
of pathology removals (correct, kept) while the two structural constants — per-box kernels and
per-box messages — remained untouched.

**D2 — Four efforts that are one architecture, never unified.** The flat store (landed), batched
AMR-local kernels (proven 12.4x on the AMR-local ops), the plan-based exchange (designed, v2), the
fused advance (T2, never started), and the S-track scaling items are the four pillars of a single
architecture — the one AMReX/Parthenon/SAMRAI all implement. Tracked as separate tracks with
separate justifications, each kept being individually deprioritized against the others.

**D3 — The stated goal has no instrument.** The goal is weak-scaling AMR; every measurement is
np=8 on one node. The weak-scaling harness (S0, ~150 LOC) does not exist; the known scaling walls
come from code audit; no increment is currently gated on a scaling metric. A program cannot
converge on a goal it never measures.

The methodology fixes already adopted (5% noise floor, deterministic byte/count gates,
phase-tree mechanization, pre-registered decision rules) stay. The strategic fix is this document.

## 2. Should we keep this implementation at all?

Three options were weighed:

- **(a) Rearchitect in place along the four pillars below — CHOSEN.** MFC already has the
  design AMReX/Chombo/BoxLib share: whole-box ownership, absolute cap, replicated box metadata with
  communication-free assignment, P2P data movement. The earlier "exascale needs a rewrite"
  conclusion was retracted after audit (the real ceiling was subdomain-derived caps, since removed).
  Two pillars are partially landed (flat store; plan-based exchange designed and I0 in
  verification). The gap is *granularity of execution*, not topology of ownership.
- **(b) Adopt AMReX as the mesh layer.** Rejected: a C++/HIP dependency across four CI compilers
  plus AMD flang; and the hardest item — rewriting MFC's solver kernels to operate on batched box
  sets — is *contained in* option (b) anyway, plus an integration layer. All of (b)'s value with
  less risk is available by adopting AMReX's *patterns* (which this document does, explicitly).
- **(c) Continue phase-share optimization.** Rejected — it is diagnosis D1.

## 3. The end state, stated as invariants

> **The "today" column below is a SNAPSHOT and drifts.** It was written 2026-08-20 and was a week
> stale when re-audited on 2026-08-27 (ledger 35: W1 is O(global boxes^2), not O(global boxes); W8
> now holds through np32). Do not trust it - **run `amr-bench/invariant_scorecard.py <arm dirs>`**,
> which computes the invariant-facing quantities from the `[amr-scale]` counters the code already
> emits. The program drifted for five days onto a wall-clock metric that is structurally blind to
> W4; a scorecard you have to remember to consult is one that stops being consulted.

A weak-scaling AMR arm satisfies, at fixed per-rank work as ranks and problem grow together:

| # | invariant | today | end state |
|---|---|---|---|
| W1 | per-rank step cost = f(local cells, peers) | seam all-pairs scan FIXED (5ec798ed). Pass-2's O(parents x global tagged cells) rescan is retired by S3.3a. **ROOT CAUSE FOUND 2026-08-28: the block metadata is REPLICATED GLOBALLY** -- `amr_region_lo_all`/`amr_region_hi_all`/`amr_owns_all`/`amr_block_level` are all allocated to `amr_max_blocks`, so every rank holds every block's geometry, ownership and level (~240 MB/rank of metadata at 1e5 ranks x 75 blocks/rank). Every `do ob = 1, amr_num_blocks` site is only writable BECAUSE of this, so W1 cannot be closed by restructuring loops: pass 1 is still O(parents x global blocks) = **O(P^2) per rank**, and S3.3c takes it to O(P) but no lower. Closing W1 needs the metadata DISTRIBUTED (own blocks + neighbour halo, remote queries through `f_amr_owner`) -- an architectural increment on the scale of the store flattening, to be scoped before further W1 loop work | f(local, peers) |
| W2 | kernel launches/step = O(levels x stages) | 14,091/step (O(boxes x stages x kernels)). **FLAT IN P** (75 boxes/rank at np8/16/32), so this is an efficiency/parity item and NOT on the scaling critical path -- also the largest and least-bounded item | ~10^2 |
| W3 | MPI messages/step = O(peers x families x levels) | measured live 2026-08-27: **F5 is 71.7 percent of ALL messages** (10,944 of 15,258) while carrying 9.8 percent of the words -- its messages average 0.8 MB against 8-65 MB for every other family. A LATENCY problem; F1/F2 aggregation is the template | O(peers) |
| W4 | no per-CELL global collectives | **level-1 gather DELETED by S3.1** (`s_amr_union_gtag` gone; ntag_bytes 185 MB -> 0 per rank per regrid, box set bit-identical). **RESTATED 2026-08-28 by differencing per regrid: the residual wall is COLLECTIVE COUNT, not volume.** The level-1 reducing path is FLAT in P in both per-rank bytes (1.39 MB) and collective count (**16,383 global `MPI_ALLREDUCE`s per regrid**, cap-fixed) -- and **99.93% of those are on nodes lying wholly inside ONE rank**, where the reducing rank already holds every tag in the subtree, so they buy nothing but a full-machine sync. Only 7 (np8) / 11 (np16) span >1 rank. The earlier `rbytes` growth of 1.67x/1.83x per doubling lives entirely in the **level-2 forest**, which does zero collectives today. **2026-08-28, CORRECTED: BOTH paths are O(P).** The level->=2 window `ALLGATHERV`, `gwin_bytes` = 360/719/1440 MB per rank at np8/16/32 — exactly 2.00x per doubling, and larger than the 185 MB gather S3.1 deleted. Level 1 is O(P) too: with the block cap raised so it never binds, level-1 tree nodes (= global ALLREDUCEs per regrid) go 13,825 / 27,537 / 54,961 = 1.992x, 1.996x, and its shared set goes 11 / 23 / 47. The earlier 'flat 16,383' and 'O(log P) shared set' readings were BOTH artifacts of a fixed amr_max_blocks = 8192 saturating in every rung (16,383 = 2*8192-1). So S3.3 (each rank clusters only the parent windows it OWNS — the forest is 99.8% rank-local at a CONSTANT fraction in P) goes FIRST, then S3.2b for level-1's flat 16,383-collective latency constant** -- and S3.3 also retires W1's regrid pass-2 loops. Prerequisites B0 and B1 have landed | tag exchange O(local + peers) |
| W5 | tag space independent of box count | the SECOND wall, at ~28k ranks (2^21 / ~75 boxes per rank). Verified 2026-08-27: **19 of 41 AMR p2p call sites still tag per box**; F1 only partially converted (I2b unlanded), migration tags by column index, and the subcycle sites are an EXPLICIT deferral to I8. Needs I2b + I7 + I8, not one increment | (family, epoch) tag bases |
| W6 | store lifecycle device-resident | a device-native growth path EXISTS; the host round trip is a deliberate fallback above `amr_grow_dev_cap = 32` columns because device-native staging transiently doubles the footprint (a measured OOM sits behind it). Fix is CHUNKED device-side migration, not deleting the host path | device-side remake, index derived |
| W7 | migration priced and bounded | verified: `s_compute_load_weight` is 1.0 inside the active box and 0.0 outside -- a CELL COUNT, not a cost model. 170 LOC, so the code surface is tiny; the blocker is choosing the model (the ledger records that a*cells + b*boxes failed its residual test and greedy remapping was worse than none). Goes LAST | work + migration terms, decoupled, hysteresis |
| W8 | per-rank DEVICE memory = f(live local boxes), never f(global index space) | **HOLDS THROUGH np=8 as of f5f99337 (2026-08-21 night, P1 pooling): np=4 peaks 30.2-39.5 GiB (was 63.6), np=8 COMPLETES for the first time at 39.0-51.3 GiB with live 72-75/rank — >12 GiB headroom on every card** | landed: in-place re-densify + capped growth + early-free + stash-only replicas (9bcc9865) + pooled q_prim/rhs advance scratch (f5f99337) |

Deliberately **kept global**: the replicated box list + owner map (~tens of bytes/box on every
rank). That is AMReX's own design — it buys communication-free assignment and is tolerable to
~10^7 boxes. Do not trade it away (audited conclusion, 2026-08-01). What must go is everything
sized per-cell globally, iterated per-box globally, or allocated O(boxes x ranks).

## 4. The four pillars, with existence proofs

**P1 — Storage: level-major flat store, device-authoritative.**
One contiguous device array per field family, dense slot index, remade device-side at regrid,
local index *derived* (search over the owned-box list), no recycle stack.
*Foreign proof:* AMReX MultiFab + RemakeLevel; Parthenon lid = n - nbs.
*Status:* flat store landed and authoritative for q_cons; bounded growth + compaction landed
(2.32x); **remaining:** device-side remap (kills the host round trip that forces the 3x/2x
hysteresis; steady state 2-3x live -> 1x), then index derivation (deletes the ratchet bug class).
**PROMOTED 2026-08-21: these two are the measured W8 fix** — at np>=2 the ratchet runs over
GLOBAL slot indices (a shifting SFC window plus received migration slots) and device-OOMs the
weak-scaling sweep at np=4; invisible at np=1 where the owned window is static. They are the
first S-track increment in practice and slot in alongside Phase 1 (independent of the exchange
conversion; ~120 LOC per the old R5/S2-tier estimates).
**LANDED same day — S0 np=4 completes; per-rank live boxes AND store capacity measured flat
across np (the W8 invariant holds).** What survived contact differs from the design: the full
device-side remake OOMed on its own staging transient and broke a host-coherence contract the
rebuild's carry-forward depends on (both churn goldens caught it). The landed form is in-place
index re-densification every reconcile + capped growth increments + early-free of consumed old
slots + stash-only replica slots, plus the `[amr-cap]` invariant instrument — full narrative in
`amr_action_plan.md` "W8 FIX LANDED". np=4's hot card peaks 63.6/64 GiB: the next memory term
is pooling per-slot q_prim/rhs (this pillar, proper) before growing the operating point.
Local-index derivation stays open as a cleanliness increment (no longer a memory one).

**P2 — Compute: one kernel per (stage, level) over the local box set.**
Flattened-prefix indexing with binary search over a cell-count prefix (the form already shipping
and mandated in the exchange design), explicit-shape or module-array access — free under the
mapped-entity law — with map(alloc:) for unavoidable private arrays.
*Foreign proof:* AMReX ParallelFor over a FabArray; Parthenon packs rebuilt only on remesh.
*Why it is required, measured three ways:* a per-block advance costs ~1x a full monolithic step
regardless of block size (launch count confirmed by kernel trace, 20x launches at 16 tiles); the
mapped-entity law prices every mapped array at 2.00 copies / ~31 us per launch, irreducible by
clauses (7 mechanisms refuted); batching the AMR-local kernels over the flat store measured
**433.5 -> 35.0 ms (12.4x), operations 15,370 -> 250**. The advance (s_compute_rhs tree: 88-91% of
kernel time, 60-66% of launches) is the one place this form has not been applied.
*Status:* proven on AMR-local ops; **the advance conversion is the largest open item and the
parity item.**
*Caveat carried from Parthenon-VIBE:* full kernel fusion alone still measured 4.4% busy on a
3-level case — fusion removes the launch term, not host-side mesh bookkeeping. P2 without P3/P4 is
not sufficient either. The pillars compose; none substitutes for another.

**P3 — Communication: plan-based, per-(family, level) exchange waves.**
The complete contract is `amr_plan_based_exchange.md` (v2, four-auditor review): SoA plans, epoch
staleness, aggregated per-peer messages, mandated pack-kernel form, persistent wire buffers,
per-family precision, hardened validator, increments I0-I8.
*Foreign proof:* AMReX FillPatch/FillBoundary; Chombo copyTo; SAMRAI RefineSchedule.
*Status:* designed; I0 (epoch, tags, asserts, the migration-stash corruption fix) in verification.
Payoff floor 8-12% of wall at the matched point; **the scale case (W3/W4/W5) is unconditional.**

**P4 — Regrid + load balance: local, scan-based, migration-aware.**
Block-lattice tag coarsening (512x volume cut at 8^3 blocks — AMReX blocking_factor), MPI_SCAN
prefix weights (8 B/rank regardless of box count — p4est), local clustering with boundary
reconciliation (SAMRAI), adapt decoupled from repartition (p4est Principle 2.1), SFC cut with a
migration-cost term and gain hysteresis (ParMETIS vsize; Uintah gainThreshold; Charm++ RefineLB —
three decades of literature agree migration-blind rebalancing is a net loss, 10-30x more movement).
*Status:* not started. **Warning to the measurer:** these score ~0 on the matched benchmark
(rg:tag 0.5%, rg:clus 2.2%); they are judged on S0 scaling metrics and deterministic bytes, never
on 400^3 wall time.

## 5. Why this reaches the goal — the arithmetic, honestly

Matched point today: 799.022 s, tax 11.03x, AMReX 3.13x (AMReX's own arithmetic floor 1.96x).

- Regrid made free (P3+P4 ceiling): 11.03x -> ~5.8x.
- Advance at our own uniform arm's efficiency (P2 target; uniform MFC is 80.7% busy, launches 81
  per step): the decomposition gives 1.60x arithmetic floor — *below* AMReX's 1.96x. The AMR
  algorithm was never the problem.
- Composite: **estimate <=5x after Phase 2 lands, 3-4x band after Phase 3** — estimates, not
  promises; each phase carries its own measured gate.
- Weak scaling: W1-W7 are what AMReX demonstrates at 10^4 GPUs. Nothing in MFC's ownership
  topology forbids them; every violation is an implementation artifact with a named fix above.

Two cheap discriminators calibrate the ordering before the big investments (Phase 0): the M2
mechanism split (is the matched-point idle local-launch-bound or MPI-progress-bound — 200^3 np=1
arm at equal cells/GCD vs the existing np=8 budget) and the CMA-off transport control (splits
posting-skew from pack/bandwidth in the exchange families).

## 6. Disposition of everything written so far

**KEEP — already end-state pieces:** flat store + device-authoritative contract; growth bound A'
+ compaction B; level-2 ISEND pool; loop-invariant coarse-halo hoist; P2P parent/child + tower
de-colocation; per-box allreduce hoist; merge-sort SFC cut; map(alloc:) clauses (hllc, weno);
phase instrumentation + phase_tree.py + [amr-scale]/[amr-mig] counters; the v2 exchange design;
seam topology relaxations; store-fix regression case.

**CARRY until replaced by the ladder:** per-box gather paths (replaced in I2-I5b); per-box advance
(replaced in Phase 2); whole-block migration sends (I4 keeps them for the exact-bytes gate;
clipping is I4b); subcycle per-box call sites (I8).

**DELETE — candidates, each a user decision:**
- The L0 coarse-tiling machinery (~1014 LOC of s_l0_ routines): measured 28-35% cost with perfect
  placement, rebalancing recovers 0.4% (inside noise), and the granularity argument is structural
  (the finest correction is one tile = 25% of a rank's load). Its goldens would go with it.
- The recycle stack + amr_loc_free/amr_loc_n bookkeeping — deleted by P1 index derivation.
- The flux families (deletion already in progress on its own evidence).

**NEVER REVISIT (closed by measurement; the ledger in `amr_action_plan.md` holds the evidence):**
packed super-grid; cap above 64; greedy/unconditional remapping; per-region map-cost reduction by
clauses; optimizing reflux (it is the sink, not the source); batch-converting blocking calls
without a downstream sync; T0 micro-items as a program; L0-sourcing; "cost grows with sim time".

## 7. The ladder, re-derived

Ordering rule: correctness first; then what unblocks measurement; then the scaling walls; then
constant factors. Every increment lands green and reports S0 metrics once S0 exists. **The local
test gate is the AMR subset (--only AMR plus the four coexist UUIDs) plus a small cross-section
when a change touches shared code — never the full 706 locally (user directive, stated three
times); the full matrix and the other compilers are CI's job on push.**

**Phase 0 — instruments and in-flight verification (days).** Results: `amr_action_plan.md`,
"2026-08-21 PHASE 0 MEASUREMENTS".
| id | item | status |
|---|---|---|
| 0.1 | migration-stash fix: verify + commit | **DONE** ca360af2, subset 65/65 |
| 0.2 | CMA-off control | **DONE** — waits +15-18% only: skew/bandwidth mechanism confirmed, sender-progress refuted |
| 0.3 | M2 mechanism split: 200^3 np=1 arm | **DONE** — rhs per-call IDENTICAL with/without MPI (17.10 vs 17.48 ms): the idle is LOCAL; P2 confirmed as the parity lever, regrid's 5.45x np=8 per-call excess is P3's |
| 0.4 | **S0 weak-scaling harness** (`amr-bench/s0_sweep.sh` + `s0_report.py`, gate arms `s0_w8gate.sh`) | **DONE, and the W8 blocker it found is FIXED (9bcc9865)** — boxes/rank and fine_work/rank flat by construction (imb 1.004); post-fix np=2/np=4 both complete with flat VRAM plateaus (55.3 / 63.6 GiB) and live 72/rank at both. **The measured gap is now TIME: wall 2.59x per np-doubling at fixed per-rank work (255.6 -> 662.3 s), with ntag/gwin/cost all doubling per np-doubling (W4). Next S0 runs: phase-budget diff of the np=2/np=4 pair (split the 2.59x BEFORE building anything), np=8 arm, np=1 on the current build.** (The earlier freg/creg ~1 GiB attribution guess was WRONG — their np-delta is ~90 MiB; the ratchet + per-slot q_prim/rhs were the real terms.) |
| 0.5 | uniform-baseline re-run (13% discrepancy) | pending |

**Phase 1 — P3 exchange (in flight; contract = `amr_plan_based_exchange.md`).**
I0 -> I6 as staged there, plus the two independent cheap scaling items folded in early: S1
block-lattice tag coarsening and S2 scan-based weights (~100 LOC combined, judged on S0 metrics).
Exit: messages O(peers); plans cached on epoch; validator green at ppn=4; bitwise goldens.

**Phase 2 — P2 batched advance (the parity item).**
- 2a prototype: batch convert_conservative_to_primitive over each level's boxes (the single
  biggest kernel in BOTH arms, 562 launches/step -> one per stage per level). Prices the full
  conversion before commitment; gate: launch count + matched wall + byte-identity.
- 2b the rhs tree: (stage, level) kernels over the flat store; flattened-prefix box indexing;
  explicit-shape/module access; map(alloc:) private arrays. The largest single investment in the
  program (m_rhs/m_weno/m_riemann, three offload backends). Risk bounded by 2a, the l0_ntile
  byte-identity sweep (the established instrument for exactly this cost), and per-family
  conversion order.
Exit: launches/step O(levels x stages); l0_ntile sweep ratio ~1 (today 16.75x at 16 tiles);
matched tax at or below ~5x.

**Phase 3 — P4 regrid + load balance.**
S3 local clustering + reconciliation; M5 decouple adapt/repartition (prerequisite for pricing
migration at all); M3 migration-cost term + A7 hysteresis (judged on deterministic bytes — the
10,748,501,376-byte reproducibility makes this an exact experiment); I7/S4 shrink the O(boxes)
loops and the O(boxes x ranks) allocations whose consumers are all converted.
Exit: per-rank regrid cost O(local + peers); S0 metrics flat across the size sweep.

**Phase 4 — completion.**
I8 subcycle conversion; execute the deletion list; full suite + four CI compilers + AMD flang;
**multi-node weak-scaling validation** (needs an allocation decision — see below); upstream PR
strategy for #1628.

## 8. Decision points (user)

- **D-phase2: RE-SEQUENCED 2026-08-27 — 2b is required but is NOT next.** The 2026-08-20 decision
  ("yes full commitment we need to get this done") stands on WHETHER; ledger 35 changed WHEN. W2 is
  FLAT in P (75 boxes/rank at np8/16/32 alike): a ~173x launch constant is a permanent efficiency
  tax at every scale, but it is not what fails at 1e5 ranks. W4 and W1 are, and both have far
  smaller fixes. Order: W4 -> W1 -> W5 -> W3 -> W2 (2b) -> W6/W7. 2a stays gated off
  (`amr_prim_batch = .false.`); its bridge-load regression is the argument FOR store-native 2b, not
  against Phase 2.
- **D-node: CONSTRAINED 2026-08-20 — nodes are scarce here.** Weak-scaling validation is designed
  single-node-first: S0 sweeps problem size and rank count (1..8 GCDs) at fixed per-rank work,
  which exposes every O(boxes) and O(P) term without a second node. Multi-node becomes a final
  spot-check if and when an allocation window exists; nothing in the ladder blocks on it.
- **D-l0: DECIDED 2026-08-27 — DELETE.** MFC's level 0 is balanced BY CONSTRUCTION (a Cartesian
  decomposition hands every rank one equal-sized chunk), so L0 tiling load-balances something that
  cannot be imbalanced in the uniform-cost case. AMReX boxes level 0 only because there boxes ARE
  the decomposition, which is also why the literature's ">=4 boxes/rank" floor does not transfer.
  Measured: tiling costs 28-35% wall, rebalancing recovers 0.6% (inside noise). **RE-ENTRY
  CONDITION, and it is live rather than hypothetical:** a measured level-0 work imbalance above
  ~10% on an IB, chemistry, or Lagrangian-bubble case, with `[amr-balance]` extended to level 0 —
  MFC's own fine-block cost model already carries `K_ib`/`K_pc` terms, so the code encodes that
  per-cell cost is NOT uniform once those physics are on; level 0 just gets no equivalent. Before
  deleting, check what the 12 L0/coexist tests cover incidentally. The L0 restart-writer filter is
  DROPPED - do not repair save/restart for machinery being removed.

## 9. Planning discipline: just-in-time contracts, one phase ahead

Only one phase carries a detailed implementation contract at a time (as of 2026-08-27: the W4 fix =
`amr-bench/notes/s3_distributed_clustering_design.md`; Phase 1's `amr_plan_based_exchange.md` is
landed and no longer the active contract). When phase N is roughly 70% landed, phase N+1 gets its own
contract at the same resolution — family/kernel inventory, data-layout contract, portability
constraints, increments with gates — written against the code as it exists *then*, and put
through the independent multi-reviewer audit ritual before its first increment (v1 of the
exchange design was wrong seven ways when audited; the ritual also found a live corruption bug).
Planning further ahead than one phase is recorded here only at the pillar/invariant level.
Concretely: the Phase 2 contract is written during Phase 1's I4-I6 stretch, seeded by the landed
pack-kernel form, the 2a prototype's measured price, and the M2 discriminator's verdict.

## 10. Standing rules (unchanged, restated)

Noise floor ~5%: single-run wall deltas below it are not results; judge deterministic bytes/counts
where possible. Pre-register decision rules before data lands. Mechanize error-prone readings.
An unexplained golden diff is a bug report. Smallest correct change; increments land green,
one commit each. Findings update this document's status columns; re-founding it requires evidence
that an invariant in section 3 is wrong.
