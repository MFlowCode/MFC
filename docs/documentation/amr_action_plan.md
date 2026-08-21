# AMR performance plan (2026-08-19 — post-diagnosis execution plan)

> **2026-08-20 RE-FOUNDING: read `amr_endstate.md` first.** The program is now derived from the
> end-state architecture (four pillars, weak-scaling invariants W1-W7), not from phase shares at
> the matched point. This document remains the evidence ledger and detailed work list; where its
> sequencing conflicts with the endstate ladder (notably: the batched advance is REINSTATED as
> Phase 2, the "kills batching-the-advance" reading was an operating-point artifact), the endstate
> document wins.

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

### The increment sequence (each = one landed permanent piece, one commit, goldens-first)

1. **I1 exchange validator** (P3, ~500 LOC) — starts now, gated on nothing: it is pure code,
   prerequisite to every exchange conversion, and M4 cannot change its scope. v2 contract + the
   inventory corrections (F1/F2 shared pool). New asserts earned this week: replicas are
   stash-only slots; `[amr-cap]` flatness (promote into `s0_report.py`). Verification includes
   the I0 lesson: seed a deliberate staleness bug and confirm the validator TRIPS.
2. **P1 q_prim/rhs pooling** — position set by M4c (see rule above). Scope: pool the two
   block-sized per-slot families into level-major store arrays; convert the rebuild overlap
   carry-forward to a device kernel IN THE SAME increment (that deletes the host-coherence
   contract, which then legalizes device-side growth). Pre-registered gate: churn goldens +
   subset 67/67 + S0 np=4 hot-card peak drops >= 10 GiB (expected ~48 GiB); np=8 attempted
   after. QBMM side-state stays per-slot (out of scope).
3. **M4-directed third increment:** S1 block-lattice tags (if regrid/collective scaling is the
   chosen front — kills the measured ntag doubling; judged on ntag bytes/rank going flat), OR
   I2-I3 exchange waves (if wait/skew traces to per-box exchange arrival), OR the S2/T1
   balancer path (only if M4b's rule fires). One front at a time; the other rules stay written
   so the choice is a lookup, not a debate.
4. **Phase-2 (P2 batched advance) contract is WRITTEN during increment 2-3** per the
   constitution — seeded by the 2a prototype (batch convert_conservative_to_primitive) and
   M2's verdict; it does not start until its contract passes the audit ritual.

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
