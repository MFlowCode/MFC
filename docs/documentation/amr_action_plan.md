# AMR performance action plan (rewritten 2026-08-18)

The working plan for closing MFC's AMR tax. This version supersedes the 2026-08-15 plan (in git
history) after two rounds of external audit; it folds in the audit corrections, the 2026-08-18
phase budgets, and the rb:mem/rb:unpk refutation. Companions: `amr_tax_review.md` (the audit and
its reasoning), `amr_block_batching.md` (chronological research log). Every number below carries
its operating point; a number without one is a bug in this document.

---

## 1. Where we are

**The gap.** At matched block size (cap 64 both codes, 400^3, np=8, 2 levels, volumetric blob,
weno1 + Lax-Friedrichs, regrid_int=2), MFC's AMR tax is **6.89x** vs AMReX's **3.40x** — a **2.03x
excess**. The once-quoted 7.64x excess was an unmatched comparison (MFC cap 32 vs AMReX 64^3) and
is retired. The tax metric structurally rewards heavier physics, so production numerics (WENO5 +
HLLC) will report a LOWER tax with nothing changed; quote both operating points.

**What is already won** (landed, measured): cap 32 -> 64 = 2.24-2.32x wall at lower memory;
`flux_n` deletion = -25%; hllc + weno `map(alloc:)` = -26.4% (3D AMR case, byte-identical);
loop-invariant coarse halo; regrid stash steps A/B; the flat store. The per-cell arithmetic penalty
of AMR is GONE at cap 64 (a = 3.50 ns/cell vs uniform 3.63) — what remains is per-box overhead and
the regrid path.

**The budget everything below allocates against** (2026-08-18, exclusive node, 400^3, np=8,
`regrid_int=2`, weno1+LF, from-scratch, `ph_wall_total`; cap 32 wall 656.975 s / cap 64 293.730 s):

| phase | cap 32 | cap 64 | inside it (cap 64) |
|---|---|---|---|
| **regrid** | 35.2% | **42.2%** | rg:build 26.1% of wall (rb:gath 15.8 = wait 7.4 + UNBRACKETED host ~8.4); rg:mig 9.4%; rest <4% |
| per-box per-step MPI (gather+reflux+seam) | 29.0% | 25.1% | gather 11.7, reflux 9.8, seam 3.6 |
| rhs (fine-block advance) | 23.1% | 18.7% | internal GPU-busy ~45-50% (ESTIMATE — see G0) |
| coarse (monolithic L0 RHS) | 2.6% | 4.9% | ~53-65% busy |
| halo / gfill / rk | 2.4% | 2.5% | |
| unbracketed | 7.6% | 6.7% | contains the per-block swap — never measured |

**Regrid's share RISES with the cap** (35.2 -> 42.2%): raising the cap made regrid relatively more
dominant, not less. And the whole table is at `regrid_int=2` + the cheapest numerics: at int=20
regrid roughly halves and rhs roughly doubles (cap-32 sweep), and production numerics push rhs
further. **The production-point budget does not exist and the priority ordering depends on it** —
which is why the plan starts with measurements, not code.

**Settled negatively, this week:** the regrid host cost is NOT the per-box allocate (rb:mem
0.002-0.010 s) and NOT the unpack (rb:unpk 0.026-0.104 s) — both brackets landed 2026-08-18 at both
caps; the scratch-hoist fix is cancelled. The ~25-64 s host remainder is in the still-unbracketed
send-side pack, per-(box,source) geometry scan, and send-pool drains. Third refuted code-read
attribution of the campaign; only brackets decide.

---

## 1b. GATE G0 — RESULTS (measured 2026-08-18, `amr-bench/logs/g0-0818_1056`)

G0 ran and **routed the plan away from where it was pointing.** Three from-scratch runs, cap 64,
400^3, np=8, one binary, exclusive node. A = matched point (`regrid_int=2`, WENO1+LF, wall 305.881 s);
B = **production point** (`regrid_int=20`, WENO5+HLLC, wall 139.940 s); C = A + kernel trace.
A and B carry **identical box counts** (64 L1 + 160 L2) and `fine_work` within 9%, so the mesh-lag
confound that invalidated earlier interval comparisons is absent.

| phase | A (int=2, WENO1+LF) | **B (PRODUCTION: int=20, WENO5+HLLC)** |
|---|---|---|
| **regrid** | **39.4%** | **37.6%** |
| gather (per-step, per-box) | 12.7% | **16.6%** |
| rhs | 17.9% | 15.0% |
| coarse (monolithic L0) | 5.4% | 13.8% |
| reflux | 10.9% | 6.8% |
| seam | 4.0% | 2.9% |
| **swap** | **0.8%** | **0.5%** |
| unbracketed | 6.1% | 3.7% |

**G0.1 — the per-block swap is 0.8% / 0.5% of wall.** It was never measured before and was assumed
large. **This refutes the central premise of the batched-advance design**: per-block grid-state
reconfiguration is not a material cost, so removing it cannot pay for the batching machinery.

**G0.2 — `PH_RHS` is 54-57%% GPU-busy**, measured on all 8 ranks with numerator and denominator from
ONE run (rocprofv3 kernel trace; coarse and fine `s_compute_rhs` separated by the exact signature
that the coarse call fires 3x per step at the rank's 200^3 subdomain while fine blocks spread over
many smaller grids; `br_load`/`br_store` at ~77/step = 28 blocks x 3 stages corroborates the split).
Tracing inflates the denominator, so this is a LOWER bound. **Batching ceiling = rhs overhead 7.8%%
+ swap 0.8%% = ~1.09x.** The audited estimate of 44-49%% was close and slightly conservative.

**G0.3 — the priority flip that the review predicted does NOT happen.** Regrid is ~38%% at BOTH
operating points and, at the production point, **2.5x larger than rhs**. Production numerics raise
`coarse` (5.4 -> 13.8%%) and `gather` (12.7 -> 16.6%%), not `rhs`.

**G0.4 — all four round-1 gather candidates measured ~ZERO** (`rb:own` 0.010, `rb:upd` 0.017,
`rb:pack` 0.043, `rb:rsv` **0.000** s at the production point). The `s_amr_gsnd_reserve` force-drain
hypothesis is refuted outright. **70%% of `rb:gath` and 12.5%% of wall inside `rg:build` remain
unexplained** - round-2 brackets (`rb:seam`, `rb:post`, `rb:geo`, `rb:slot`, `rb:tail`) are in.

**BONUS RESULT — regrid frequency is a weak lever.** Cutting frequency 10x (int 2 -> 20) cut total
regrid cost only **2.3x** (120.6 -> 52.7 s), because per-regrid cost rose **4.4x** (6.03 -> 26.33 s)
as the mesh drifts further between rebuilds. "Regrid less often" is not a mitigation, and the old
"tax 27.2x -> 7.2x" framing overstates what the interval buys.

### ROUND 2 (same day): rg:build CLOSES; the gather's host half does not

Five more brackets, both arms re-run. **`rg:build` now accounts to 100.0%%** (residuals 0.003 s and
0.001 s) - the 22.5%% that was unexplained is two regions:

| bracket | matched | **production** | imbalance | what it is |
|---|---|---|---|---|
| **`rb:tail`** | 5.1%% | **8.8%%** | **2.6-2.7** | post-loop `send_flush` (WAITALL) + `reduce_xchg_flag` (ALLREDUCE) + reconcile + topology check |
| `rb:slot` | 1.5%% | 4.2%% | 1.21 | `s_amr_alloc_slot` |
| `rb:seam` | 0.002 s | 0.001 s | - | **REFUTED** - the O(N^2) `s_amr_build_seam_pairs` was the top candidate |
| `rb:post` | 0.028 s | 0.018 s | - | REFUTED |
| `rb:geo` | 0.002 s | 0.001 s | - | REFUTED - `s_set_amr_fine_geometry` runs per box on every rank and is free |

**`rb:tail`'s imbalance of 2.6-2.7 is the highest in the entire budget** (mean 11.7 s, max ~30 s at
the production point). A WAITALL plus an ALLREDUCE with that spread is ranks ARRIVING at different
times and blocking, not work - the per-box rendezvous skew resurfacing at the drain point.

**`rb:gath`'s host half is STILL unexplained: 11.9%% of wall at the production point**, the largest
single unknown in the budget. Only 32.1%% of `rb:gath` is accounted there. **That is SEVEN code-read
attributions refuted by their own brackets against two real finds** - a ~22%% hit rate, which is the
strongest argument in this document for bracketing before fixing. Round 3 brackets the only two
things left in that routine (the non-owner `MPI_ISEND` on 1.5-3 MB rendezvous-sized messages, and two
scalar geometry calls) and splits `rb:tail` into its four collectives.

### ROUNDS 3-4: THE ROOT CAUSE, AND THE REGRID PATH NOW CLOSES 100%%

**Round 3** killed the last two candidates (`rb:send` 0.017 s; `rb:tail` is **99.98%% one
`MPI_ALLREDUCE`**, `rb:flush` 0.000 s) and left 11.9%% of wall unexplained with no code left to blame.

**Round 4 found why every candidate measured zero.** `s_amr_gather_coarse_patch` returns at its FIRST
branch for level >= 2 blocks, into `s_amr_gather_from_parent` - so all ten rb:* brackets instrument
only the level-1 path, **64 of 224 boxes. The other 160 (71%%) were never measured.**

Confirmed against a PRE-REGISTERED prediction (`amr-bench/PREREGISTERED_pgather.md`): predicted
`pg:all` = 11.9%% of wall at the production point, **measured 11.9%%** (matched point: predicted
22.9 s, measured 26.6 s, +16%%, inside the 20%% rule). `pg:send + pg:recv = pg:all` exactly, and
`rb:gath`, `rg:build`, `rb:tail` and `pg:all` now all account to **100.0%%**.

**THE DEFECT: the level-2 parent gather uses a BLOCKING `MPI_SEND`/`MPI_RECV` once per box**
(`m_amr.fpp:1305` and `s_amr_recv_parent_patch`). The codebase documents this exact pathology at
`m_amr.fpp:256` - *"the non-owner side used a BLOCKING MPI_SEND, so every rank had to rendezvous with
the owner 794 times per rebuild, in lockstep ... measured at 45%% of regrid and ~25%% of total
runtime"* - and the fix (`MPI_ISEND` + `amr_gsnd_pool` + one drain per rebuild, `m_amr.fpp:1025`)
**was applied to the level-1 path only.** Level 2 is the majority path.

### TRACK R WORK LIST - measured, production point (int=20, WENO5+HLLC, cap 64, np=8)

| target | %% wall | calls/rank | ms/call | imbal | fix |
|---|---|---|---|---|---|
| **`rb:xchg`** one-flag ALLREDUCE | **8.9%%** | 2 | **5857** | 2.62 | SYMPTOM - pure arrival skew; fix upstream, never the allreduce |
| **`pg:send`** blocking MPI_SEND | **8.6%%** | 27 | **420** | 2.02 | ISEND + pool + one drain (in-tree precedent at :1025) |
| `rg:mig` stash + migrate | 6.7%% | - | - | 1.05 | not yet decomposed |
| `rb:wait` per-box WAITALL (L1) | 5.0%% | 15 | 441 | 1.54 | batch the recv side |
| `rb:slot` `s_amr_alloc_slot` | 4.2%% | - | - | 1.21 | |
| **`pg:recv`** blocking MPI_RECV | **3.3%%** | 19 | 228 | 2.95 | IRECV + drain |

**`pg:send` + `pg:recv` + `rb:wait` = 16.9%% of wall is per-box BLOCKING point-to-point, and
`rb:xchg`'s 8.9%% is the arrival skew that serialisation produces. Together ~26%% of wall - two thirds
of the whole regrid phase - is one defect with one known fix.**

Not yet priced, same family: blocking `MPI_SEND` at `m_amr.fpp:1219` (level-1 host path), `:2874`
(reflux to parent), `:2908/2909` (flux registers to parent).

### What G0 decided

- **Track R (regrid) is FIRST**, at both operating points, by a factor of ~2.5 over anything else.
- **Track B (batching) is DEFERRED** by the plan's own >50%%-busy rule. The widened bridge is set
  dormant (`amr_br_batch = 1`, restoring the one-block allocation); the full hunk is preserved in
  `amr-bench/bridge_and_brackets.patch` if it is ever revived.
- **Track M rises**: per-step `gather` is 16.6%% at the production point, now second-largest.

---

## 2. Gate G0 — measurements that route the plan (hours; do before any track)

1. **Bracket the per-block swap.** `s_amr_swap_to_fine()` sits OUTSIDE `PH_RHS` while
   `s_amr_restore_coarse()` sits INSIDE (`m_amr.fpp`, `s_amr_fine_stage_rhs`) — the bracket is
   charged half a swap pair, and the batched design's single largest claimed benefit has never been
   measured. Add a symmetric `PH_SWAP` around both halves; move the restore out of `PH_RHS`.
2. **Kernel time INSIDE `PH_RHS`, one run, one clock.** rocprofv3 kernel trace intersected with the
   bracket on the current binary. The composed estimate (~45-50% busy) went through three audit
   layers and shifted each time; it must be replaced, not refined. This number decides Track B's
   ceiling: busy < 30% -> batching is first-order here; busy > 50% -> batching waits.
3. **The production-point budget.** Same 400^3 case at `regrid_int=20`, WENO5 + HLLC, cap 64,
   from-scratch, phase report on. One run, no comparison arm. This is the budget the paper's claims
   will live at, and every priority below is provisional until it exists.
4. **Bracket the gather's remaining host work**: send-side pack, geometry scan, pool drains (gate
   with the existing `amr_rg_gather` flag). Completes the rb decomposition that E2 started.

---

## 2b. R1 LANDED AND MEASURED (2026-08-18) — one call site, -17.3%% wall

`m_amr.fpp:1305`, level-2 parent gather: blocking `MPI_SEND` -> `MPI_ISEND` into the existing
`amr_gsnd_pool`, drained by the `s_amr_gather_send_flush` already present after the rebuild's box
loop. The level-1 path has used exactly this since the defect note at `m_amr.fpp:256`; R1 applies it
to the majority path.

**Matched point, before -> after** (wall **338.545 -> 280.042 s, -17.3%%**):

| phase | before | after | delta | calls (both) | ms/call |
|---|---|---|---|---|---|
| `pg:send` (the converted site) | 18.368 | 4.663 | **-74.6%%** | 89 | 206.4 -> **52.4** |
| `pg:recv` | 8.263 | 3.956 | -52.1%% | 63 | 131.2 -> **62.8** |
| **`rb:xchg`** (never touched) | 17.552 | 8.020 | **-54.3%%** | 5 | 3510 -> **1604** |
| **`rb:wait`** (level-1 path, never touched) | 23.335 | 23.052 | **-1.2%%** | - | NULL CONTROL |
| `regrid` | 128.379 | 101.150 | -21.2%% | 20 | - |

**Why this is more than one noisy pair.** (a) **Call counts are IDENTICAL** (89/63/5) - the workload
is unchanged in every countable respect, so only per-call cost moved. (b) **A built-in null control**:
`rb:wait` is the sibling level-1 receive path, untouched, and it moved -1.2%% while the converted path
moved -74.6%% - global drift cannot explain that. (c) **A pre-stated prediction held**: if the per-box
blocking sends CREATE the arrival skew, then `rb:xchg` should shrink although never touched. It halved.
That establishes the causal chain (serialisation upstream -> skew downstream), not just correlation.

**Correctness.** 76/76 AMR/store goldens pass, and four of them (`multi-level static np=2` 1D and 2D,
`multi-level dynamic regrid np=2`, `2D dynamic regrid np=2`) genuinely execute the converted line -
at np=1 `powner == cowner` takes the local device-copy branch and the send is never reached, so
single-rank tests would have proved nothing.

**A DEFECT IN R1, FOUND BY CALL-SITE AUDIT AND FIXED SEPARATELY.** `s_amr_subcycle_setup_level`
(`m_amr.fpp:5140/5146`) calls the two Fypp-generated variants DIRECTLY, bypassing
`s_amr_gather_from_parent`, and that loop has NO drain: sends would accumulate to the pool's 64-slot
force-drain, whose WAITALL blocks until every matching receiver has posted - a receiver itself stuck
behind another rank's force-drain closes a deadlock cycle. **No test covers it**: every subcycle
golden is 1D at np=1, i.e. the co-located branch. Fixed by draining immediately at both sites (the
same "keep this site's original blocking semantics" treatment seven other runtime sites already
carry), so the regrid path keeps the benefit and this path keeps its original behaviour.
**Lesson: auditing the wrapper's callers is NOT enough - Fypp-generated routines have call sites the
wrapper does not.**

**Caveat:** one before/after pair against a ~12%% node noise floor. The phase deltas and ms/call
figures are the robust claims; the -17.3%% wall is suggestive and needs interleaved pairs before it is
quoted anywhere external.

---

## 3. Track R — the regrid rebuild (the largest measured item)

Regrid is 42.2% of wall at cap 64/int=2 and ~18% at int=20 (cap-32 measurement); MFC pays ~12 s per
regrid where AMReX pays ~0.027 s for the same logical operation — two to three orders of headroom.
Fixing it also removes the regrid-frequency asymmetry (tax 27.2x -> 7.2x going int 2 -> 20), which
currently makes the reported tax hostage to a tuning knob.

- **R1. Aggregate the rebuild gather across boxes.** `s_amr_gather_coarse_patch` runs once per box
  (~224-1207 per rebuild) with an owner-side IRECV + immediate WAITALL each. The deferred-SEND pool
  with a single post-loop drain already exists and works on this path; the recv side needs the same
  treatment: post all boxes' receives, then one wait-and-unpack pass (needs per-box `amr_cg`
  landing buffers or direct-to-slot unpack). Targets rb:wait (7.4% of wall) plus the serialization
  skew it creates (rb:wait imbalance 1.59).
- **R2. Whatever G0.4 finds in the ~8.4% unbracketed host remainder.** No fix is chosen until the
  bracket lands — the last three pre-chosen fixes on this path were all refuted by their brackets.
- **R3. Device-resident rebuild.** The rebuild is host-based end to end: host prolong
  (`s_prolong_one_var`), host overlap carry-forward, and the stash/migrate round trip (rg:mig,
  9.4% of wall) that moves full slots through the host even at np=1 where migration never fires.
  The conditional-migration design (pull to host only blocks that actually migrate) is recorded and
  believed right but HAS NEVER EXECUTED; steps A/B de-risked the kernel-placement traps (device
  kernels must live in small separate routines; gate on run rc, not `strings`). This is the
  structural fix and the one that approaches AMReX's 0.027 s.

Verification: counts first (gather Waitalls per rebuild must drop ~boxes-fold), wall second, at
BOTH int=2 and int=20; goldens + an np=2 run (migration paths are invisible at np=1).

---

## 4. Track B — the batched fine-block advance

The design is gated and groundwork landed (flat store, widened bridge, BC escape hatch verified,
uniform slot shapes, Cartesian-only coordinate audit). The batched-kernel MWE returns 12.4x — but
that is an UPPER BOUND (plain contiguous array, no bridge staging), and the in-situ ceiling at
cap 64/int=2 is only ~1.11-1.20x if `PH_RHS` is really ~half busy. **G0.2 and G0.3 size this track;
do not start it before they land.**

If it proceeds: goldens on the widened bridge first (never run); re-measure the bridge penalty with
>= 8 interleaved pairs on an exclusive node (the +24.8% at n=3 was a converging transient); sweep
`amr_br_batch` 2/4/8 with the batched path unused; then wire `s_amr_br_load_all` -> tall (m,n,p) ->
one `s_compute_rhs` -> `s_amr_br_store_all`, batch 8, per-block fallback retained. **Gate on
counts** (rhs-family launches/step must fall ~batch-fold, GPU busy must rise); wall at pilot scale
is below the noise floor. If G0 defers this track, REVERT the widened bridge — it allocates 8x
bridge memory (~880 MB at cap 64) that nothing uses, measured slower, goldens never run.

---

## 5. Track M — per-box per-step MPI (25.1% of wall at cap 64)

The per-step fill gathers (11.7%), reflux (9.8%, ~11 Recv per box), and seam halo (3.6%, imbalance
1.9) are per-box blocking chains — the same shape Track R1 fixes at regrid time. The level-
aggregated exchange (post-all per level, one drain — the FillPatch pattern) is the fix. **Sequence
it after R or B frees the host**: the Waitall-hoist experiment removed 53.8% of Waitalls for zero
wall because the host was busy anyway; these waits become critical path only once it is not.
Verify by counts (Waitalls/step), which reproduce to ratio 1.000 across 52%-different walls.

---

## 6. Track C — cheap configuration and comparison items (any order, fill idle node time)

- **C1. Cap 96, one from-scratch run.** "96/128 OOM" is UNVERIFIED — it descends from the sweep
  table that also listed cap 64 as OOM (the checkpoint-restart confound). Expect marginal (~1.48x
  memory headroom, scratch scales as cap^3, ~66 boxes over 8 ranks is poor granularity), but 32->64
  bought 2.24x for zero code. Gate on exit 134 + `HSA_STATUS_ERROR_OUT_OF_RESOURCES`, not a hang.
- **C2. The honest AMReX comparison at the production point.** Cap 64 both sides, same realistic
  regrid_int both sides, landed fixes. Close the mesh-lag confound by scaling `amr_buf` with the
  interval — and because `amr_buf` doubles as the box-MERGE threshold (`thr = buff_size +
  2*amr_buf`), matched `fine_work` must be an iterated GATE with box count reported beside it.
- **C3. Split-finder midpoint fallback.** Real defect: convex blobs return their bounding box
  (~91% over-cover for a sphere) and `amr_cluster_eff` is dead code. No valid break-even model
  exists (`b` is not a constant) — decide by direct A/B on a fixed problem.
- **C4. np=1 of the matched case.** Splits MPI-progress idle from local launch idle at the
  operating point that matters; sets expectations for what R/B can recover at np=8.
- **C5. Cross-backend probe (needs an allocation).** Same case on CCE offload or an NVIDIA lane:
  launches/step, dead time per launch, tax. Bounds the amdflang/ROCm share of the toll.

---

## 7. Decision rules

- ~~G0.3 reallocates everything~~ **FIRED: regrid leads at the production point (37.6% vs rhs 15.0%).
  Track R is first.**
- ~~G0.2: busy < 30% -> Track B first-order~~ **FIRED: 54-57% busy -> Track B waits, bridge dormant.**
- Remaining live rule: no fix on Tracks R/M is chosen until its bracket lands. **Rounds 1-2 refuted
  SEVEN code-read candidates and found two real ones (`rb:tail` 8.8%%, `rb:slot` 4.2%%).**
- **Two gates `precheck` cannot provide, now run before every build** (it passed 9/9 twice on
  non-compiling trees): assert `PH_N` equals the `PH_NAME` entry count, and assert every `PH_*` used
  in a file is imported by that file's `use m_phase_timing, only:` list.
- Any fix on Tracks R/M lands only with: counts moved as predicted, goldens green, an np=2 case
  exercised, and wall measured at both int=2 and int=20 on an exclusive node.
- A code-read attribution is a hypothesis for a BRACKET, never for a fix. Three were refuted by
  their brackets this week alone.

## 8. DO NOT DO (all measured; do not re-litigate without new evidence)

- Load balancing (oracle ceiling 5.1% of wall; box/weight balance 1.008; skew explains 7.6% of MPI).
- `nowait` chains and cross-chain concurrency (both disproved in dedicated MWEs).
- Mapping-clause mechanisms (all seven fail; module-scope derived types are 3.5x WORSE), and
  further descriptor-copy reduction beyond the landed hllc/weno clauses (convert_conservative was
  correct, mechanism-verified, and bought +0.1% — reverted).
- Byte-count reduction as a goal (removing 12% of bytes bought zero wall) and Waitall-count
  reduction alone (removing 53.8% bought zero wall — hence Track M's sequencing).
- Runtime knob tuning (12 configurations; the harmful ones registered 1.8-3.2x, the rest are null).
- Scoping anything from cap-32 phase shares or from the retired 7.64x/23.92x/443x numbers.

## 9. Metric discipline (the canon; violations produced every retraction)

- Perturb and regress; never decompose a residual.
- One clock per ratio, and **never compose a ratio from two runs** — this produced a 2x-wrong
  batching ceiling during the 2026-08-18 audit, and its correction was itself off (wall used as
  kernel) until a third pass. If numerator and denominator cannot come from one run, the honest
  output is an experiment.
- A tax figure carries its arithmetic intensity and its operating point (cap, regrid_int, np,
  fixes, physics).
- Judge cap changes by WALL on a fixed physical problem, not ns/cell.
- Union-of-intervals for kernel busy time, never median x count.
- Never take wall from an instrumented run (rocprofv3 inflates MFC 1.3-2.4x).
- Counts before wall: counts reproduce to ratio 1.000 across runs whose walls differ 52%.
- From-scratch runs for anything that redetermines the hierarchy; a restart pins the state the
  parameter is supposed to set.
- ~12% single-run noise floor on an exclusive node; include a known-null arm in every sweep.
- Zero-time phases are silently dropped from the report; a missing line means a stale binary, not
  missing code.
