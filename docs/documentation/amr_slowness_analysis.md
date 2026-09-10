# Why MFC is slow with AMR: the unified analysis (2026-08-19)

Synthesis of the full measurement campaign plus a five-reviewer panel (advance-path architecture,
communication structure, memory/regrid, AMReX source comparison, rhs-skew mechanisms), each reviewer
reading actual code — MFC's and AMReX's — under read-only constraints, with the campaign's refuted-
ideas list supplied so nothing dead was re-proposed. Companions: `amr_action_plan.md` (what to do),
`amr_tax_review.md` (measurement audit trail). This document answers WHY.

Reviewer priors were committed to disk before any report was read; section 7 records where they
were wrong.

---

## 1. The one-paragraph answer

MFC's AMR treats **the block** as the unit of communication, scheduling, and memory, walked on **one
global ordered list** that every rank traverses in lockstep with a blocking element at nearly every
position. That design (i) pays fixed per-box costs O(boxes) times per step where AMReX pays them
O(neighbor-ranks) times per level with cached communication graphs; (ii) **amplifies** any rank-local
cost divergence into convoys — head-of-line blocking on rendezvous-class messages makes everyone's
per-call wait grow ~10x while counts and bytes stay flat; and (iii) **ratchets memory** to the
high-water union of every mesh the run has ever had, which both OOMs long runs (terminal, at the
64-to-128 capacity doubling) and — via VRAM pressure on the per-launch allocation path — appears to
be the *source* of the rank-local slowdown that the convoys amplify. AMReX's units are **the level**
(aggregated, cached, nonblocking exchanges with one waitall) and **the arena** (fully reusable
memory), which is why its tax is flat in regrid frequency and simulation time where MFC's grows.

## 2. The measured facts any explanation must fit

| # | fact | status |
|---|---|---|
| F1 | Tax grows with sim time at CONSTANT mesh: 4.02x (steps 40-80) to 12.14x (80-160), identical 224 boxes, fine_work, rhs calls | replicated |
| F2 | The growth is per-call WAIT: reflux ms/call 2.5 to 28 (10x) at 2.33x calls; participation only +17% | measured |
| F3 | Growth is paid at blocking sites with LOW, falling imbalance (reflux 1.32 to 1.16) | measured |
| F4 | The divergence ORIGINATES in rhs: imbalance 1.09 to 2.90; PH_RHS is MPI-free and GPU_WAIT-bracketed, so it is genuine rank-local time | measured + code-verified |
| F5 | Coarse rhs echoes it small: +77% with near-flat imbalance (1.04 to 1.06) | measured |
| F6 | Memory OOMs by REGRID COUNT: 40 regrids dies, 8 over 2x the steps survives | confirmed behind a VRAM gate |
| F7 | Removing one per-box rendezvous (R1, level-2 gather) bought -17/-22% wall AND shrank the untouched allreduce 54-84% | measured, causal prediction held twice |
| F8 | Removing a rendezvous with NO downstream sync (step-1 hoist) improved its phase 14% and WORSENED wall — cost relocated | measured, both operating points |
| F9 | Cells balanced 1.02, boxes 1.15, while rhs seconds diverge 2.9x | measured |
| F10 | Per-regrid cost ~6-26 s vs AMReX ~0.027 s | measured |

## 3. The causal chain (best-supported model)

```
 regrid churn (blob advects)
   └─> rank-local slot high-water creeps up (amr_loc_n never decrements; frees only recycle)
         └─> store capacity doubles rank-locally (s_amr_st_reserve, newcap = max(2*oldcap, nloc), no shrink)
               ├─> TERMINAL: the 64->128 doubling must allocate 28.2 GB into ~14 GB free  ==> OOM (F6)
               └─> VRAM pressure on ratcheted ranks
                     └─> per-launch map(alloc:) allocations in the rhs kernel family hit a slow path
                           └─> rank-local rhs slowdown, ~2100 launches/step exposed (F4), coarse echo at ~75 launches (F5)
                                 └─> the per-box blocking lattice (rendezvous msgs, fixed tags, one global order)
                                     CONVOYS the divergence: everyone's wait grows, imbalance stays low (F1,F2,F3)
```

Two links are proven (ratchet arithmetic; convoy amplification is the only reading consistent with
F2+F3+F7+F8). The middle link — VRAM pressure causing the rhs slowdown via the per-launch allocation
path — is the leading hypothesis (E-H1), currently being adjudicated by a per-rank instrumented run,
and has a 2-line decisive A/B (revert the two `map(alloc:)` clauses, rerun 160 steps).

## 3b. THE PER-RANK RESULT: it is ONE SICK RANK (measured 2026-08-19, after this model was drafted)

Per-rank phase times, production point, 160 steps (`amr-bench/logs/rank-0819_1156`):

| rank | 0 | 1 | 2 | 3 | 4 | **5** | 6 | 7 |
|---|---|---|---|---|---|---|---|---|
| fine_work (M cells) | 23.5 | 23.5 | 23.0 | 23.1 | 23.5 | **23.6** | 23.5 | 22.5 |
| **rhs (s)** | 225 | 240 | 240 | 259 | 253 | **976** | 240 | 247 |
| **reflux (s)** | 853 | 844 | 820 | 804 | 808 | **71** | 823 | 820 |

Rank 5 does IDENTICAL work (+1.4%% on cells) and takes **4x longer in rhs**; the other seven then
wait ~820 s in reflux while rank 5 waits 71 s. Totals are near-equal (1047 vs ~1065 s) - a clean
mirror. `max/mean` = 2.91, which IS the aggregate 2.90. At 80 steps there is no straggler at all
(imbalance 1.06), so the pathology **develops on one rank between steps 80 and 160**.

**This sharpens the model in section 3 rather than replacing it.** The chain still holds, but the
final link is now specific: it is not a diffuse per-rank divergence, it is ONE rank falling off a
cliff, with the convoy converting that into seven idle ranks. Consequences:
- **Composition is dead** as an explanation (work is balanced to 1.4%%).
- **The cost-weighted balancer (section 6.4) is DEPRIORITISED** - the load is already balanced;
  rebalancing cannot help a rank that is slow at equal work.
- **The convoy is confirmed as transmission, not source** - and still worth fixing, because it is
  what turns one sick rank into eight.
- E-H1's pre-registered signature (fine_work flat + rhs bimodal + rank-local) is matched exactly.

**Decisive test in flight:** the store trip-wire asks whether rank 5 carries a HIGHER `amr_st_cap`
than its peers. If yes, the memory ratchet and the rank pathology are ONE mechanism and section 6.1
becomes the whole fix. If no, E-H1 dies and the rank-local cause must be found elsewhere.

## 4. The four structural deltas vs AMReX (from source, both codes)

1. **Exchange granularity and caching.** AMReX `FillPatchTwoLevels`: one coarse-patch MultiFab, one
   `ParallelCopy` in, one out; the communication graph is CACHED per (BoxArray, DistributionMapping)
   and re-used; `PostRcvs` posts ONE Irecv per source RANK with all boxes packed. MFC: one
   gather per BOX per stage (~2,500-3,500 rendezvous-class p2p messages/step at 224 boxes vs
   ~150-250 for the aggregated equivalent), overlap geometry recomputed per call. Explains F7, F10,
   and why MFC's tax responds to regrid_int while AMReX's does not.
2. **Reflux.** AMReX `FluxRegister::Reflux`: six aggregated per-LEVEL `copyTo` calls plus one fused
   kernel. MFC: per BOX, 6 ISENDs per participant + WAITALL on the owner, 6 blocking RECVs per
   participant. This is the lattice where F1-F3 are paid.
3. **Regrid.** AMReX: unchanged levels are not touched at all; `RemakeLevel` fills new MultiFabs with
   one FillPatch, interpolation on device, migration implicit; old fabs return to a reusing `CArena`.
   MFC: per-box collective gather in the rebuild, host-based prolong history, stash/migrate host
   round trips, grow-only uniform-slot store. Explains F6, F10.
4. **Load model.** AMReX's DEFAULT is also cell-weighted SFC — so "MFC weights the wrong quantity" is
   NOT by itself the delta. The delta: aggregated nonblocking comm TOLERATES skew (absorbed at one
   waitall); MFC's per-box blocking chains AMPLIFY it. AMReX additionally ships runtime-cost
   rebalancing (`makeKnapSack(rcost_local, ...)`) for when apps need it; MFC has no such hook.

## 5. Ranked root-cause hypotheses (with adjudication status)

| rank | hypothesis | evidence for | discriminator | status |
|---|---|---|---|---|
| 1 | **Convoy amplification** on the global lockstep block list (4 walks per stage-cycle, blocking element in each) | F2, F3, F7, F8; two reviewers independently | timestamp rank arrivals at 3 list positions; spread grows along the list | strongly supported |
| 2 | **Store ratchet -> VRAM pressure -> slow per-launch alloc path** as the rhs-divergence source | F4, F5, F6; unifies memory + skew | per-rank data: fine_work flat, rhs bimodal, same slow ranks both windows; then the 2-line map(alloc:) revert | leading; being adjudicated |
| 3 | **Device-heap fragmentation** (same ratchet, different mechanism) | same as #2 | map(alloc:) revert helps (#2) vs does not (#3); fix identical either way | fallback for #2 |
| 4 | **Ownership composition drift** (L2 share / block count per rank re-dealt at each regrid) | windows contain 2 vs 4 regrids; box imbal 1.15 | per-rank regression: rhs seconds vs block count and cells | partial at best (1.15 cannot give 2.9x) |
| 5 | Data-dependent kernel cost (blob steepening) | grows with time by construction | per-rank: slow ranks track the interface and rotate | largely ruled out (smooth single-fluid advection; coarse imbalance flat; gfx90a full-rate denormals) |
| 6 | GCD thermal/clock divergence | long runs | slow ranks correlate with clocks, not regrid count; rotate across replicas | disfavored (replication matched to 4%) |

Also surfaced, not yet priced: an **unbracketed per-step per-box blocking chain** after the stage
loop (`restrict_to_parent`, `freg_to_parent`, `restrict_fine_to_coarse`, m_time_steppers ~:783-800)
— same defect shape, invisible to the budget, feeding each step's exit skew into the next step's
entry (a candidate for the super-linearity); and a per-rhs-call **capture sweep over the global
block list** with an ~11-array metadata push on every rank (m_amr_registers ~:586-616) — a fixed
launch-tax term at O(global boxes) per call.

## 6. Fix architecture, in order

1. **Stop the bleed (bounded, this week).** (a) Free stale slots BEFORE allocating new ones in the
   rebuild (or reconcile mid-loop) so the union spike never reaches the store-growth path; add a
   shrink-on-reconcile or a capacity cap with eviction. (b) Per-regrid stderr trip-wire:
   `amr_loc_n`, `amr_st_cap`, live-slot count per rank. (c) Bracket the unbracketed restrict/reflux
   chain (PH_RESTR) so 6% of wall stops being invisible.
2. **The decisive cheap experiments.** (a) The per-rank run in flight (adjudicates #2 vs #4/#5).
   (b) The 2-line `map(alloc:)` revert at 160 steps (adjudicates #2 vs #3). (c) B's freg pre-post
   pilot: IRECVs for all participating blocks posted at stage start with per-block tags, one
   WAITALL before the apply loop — buffers are already per-block indexed, so this is aggregation-
   lite with no refactor, and it doubles as the convoy discriminator.
3. **The structural fix (the campaign's endgame).** Per-LEVEL aggregation of the four per-box
   exchange families (fill gather, seam halo, reflux, restrict) with cached communication graphs —
   the FillPatch/FluxRegister shape. This is what F7/F8 jointly demand: aggregation removes cost
   AND keeps a sync point, where piecewise deferral only relocates cost.
4. **Balancer, conditional.** Only after #2's verdict: if composition matters, runtime-cost
   weighting (AMReX `makeKnapSack` shape). Not before — with convoys active, a better balance would
   be largely absorbed by the lattice anyway.

## 7. What the pre-registered priors got wrong (recorded for honesty)

- "Data-dependent kernel cost, 55%" — too high for THIS case: relax/igr/adap_dt/chemistry are all
  off, the blob is smooth single-fluid, and the coarse-grid echo has flat imbalance. Demoted.
- "Low imbalance at reflux = real work" (earlier session reading) — wrong inference: a convoy
  distributes waiting uniformly, so low imbalance at a sink is what a chain-propagated straggler
  looks like. The conclusion (skew paid at blocking sites) survived; the reasoning did not.
- The balancer-currency delta ("MFC weights cells, cost is physics") — weakened: AMReX's default is
  also cell-weighted. The operative delta is skew TOLERANCE of the communication design.
- Missed entirely until reviewer A: the unbracketed post-stage chain and the windows' unequal
  regrid counts (2 vs 4), which entangle composition drift with time in F1.
