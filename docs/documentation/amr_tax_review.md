# AMR tax campaign: an outside review (2026-08-18)

A fresh-eyes review of the whole investigation into why MFC's AMR arm pays a much larger tax than
AMReX, written after re-reading the full evidence base. Companions:
`docs/documentation/amr_action_plan.md` (the action list, authoritative on current numbers; rewritten 2026-08-18 to incorporate this review and its audits) and
`docs/documentation/amr_block_batching.md` (the chronological research log). Where this document and
the action plan disagree on a number, the action plan wins; this document's job is synthesis,
critique, and a prioritized path forward.

**Revised 2026-08-18 (same day) after an audit of this document's own arithmetic.** One claim made
in the first draft — that the `rhs` bracket is 62-74% GPU-busy and batching's ceiling is ~1.13x — was
WRONG: it charged the coarse-grid `s_compute_rhs` (which has its own `PH_COARSE` bracket) to the fine
`PH_RHS` bracket. Corrected numbers are in section 3; the batching prize is about twice what that
draft said. The error was a cross-run composition, which is the campaign's own documented failure
mode, and it is the reason section 7 now opens with two direct measurements instead of an argument.

---

## 1. Verdict

**The diagnosis is done, and it is correct.** The problem was never "one bug to find," and that is
precisely why weeks of increasingly good instrumentation kept producing nulls. What the evidence
shows is a **design property, not a defect**: MFC's AMR arm advances the hierarchy one block at a
time through a host-orchestrated, fully synchronous pipeline, and every per-block operation — kernel
launch, argument mapping, grid-state swap, metadata sync, blocking gather — carries a fixed host
toll of order 15-260 us that a 70-110 us kernel cannot hide. The cost is smeared across thousands of
operations per step of roughly five different kinds. There is no single culprit because **the sum is
the culprit.**

AMReX shows the same shape from the other side, and the QUALITATIVE contrast is what to carry: it
also launches per-box, but its communication is aggregated per level and its per-box host work is a
pointer capture rather than a global grid-state reconfiguration, so its per-launch dead time is
roughly two orders of magnitude below ours.

**Do not quote the specific figures that once expressed this** — "1,237 us dead per launch",
"2,582 vs 14,091 launches/step", "443x more dead time". All are from the cap-32-vs-cap-64 mismatch
described in section 2, and at matched block size the excess is 2.03x, not 7.64x. The mechanism
survived the correction; those numbers did not.

**The campaign made real progress even though it felt like none:**

- cap 32 -> 64: **2.32x wall** on the same physical problem, lower memory (Tier 0.1, landed logic);
- `flux_n`/`flux_gsrc_n` deletion: **-25% wall** (724ef4fd);
- hllc + weno `map(alloc:)`: **-26.4% wall** cumulative on the 3D AMR case, byte-identical;
- the loop-invariant halo hoist (~2.4x on its case), the regrid stash fix, the flat store;
- and, critically, the headline gap itself collapsed under scrutiny: the "7.64x excess vs AMReX"
  was an unmatched comparison (MFC at cap 32 vs AMReX at 64^3 blocks). **At matched cap 64 the
  excess is 2.03x** (MFC 6.89x vs AMReX 3.40x, each vs its own uniform arm).

**Recommendation in one line: stop diagnosing, build the level-batched advance pilot, and gate it
on operation counts rather than wall.** The measurement machine is now better than the questions
being asked of it; the decision-relevant uncertainty has moved to "does batching deliver in situ
what the MWE promises," and only building it answers that.

---

## 2. The current numbers (retire the stale ones)

Quoted numbers in circulation span three operating points and two eras of fixes. The ones to carry
forward, each with its operating point stamped:

| quantity | value | operating point |
|---|---|---|
| MFC tax vs own uniform | **6.37-6.89x** | cap 64, matched blob, regrid_int 2, np=8 |
| AMReX tax vs own uniform | 3.40x | same, matched cap |
| **MFC excess over AMReX** | **2.03x** | matched cap 64 (4.15x at matched cap 32) |
| arithmetic term | **MFC 1.90x vs AMReX 1.92x — identical** | **cap 32**, union-of-intervals estimator |
| per-cell coefficient a | 3.50 ns/cell at **cap 64**, **equal to uniform's 3.63** | 4-cap sweep, clean instrument |
| per-box coefficient b | 6.29 -> 15.24 ms/box, **rises with cap** | same sweep; b is not a constant |
| regrid-frequency asymmetry | MFC tax 27.2x -> 7.2x going int 2 -> 20; AMReX 8.57 -> 8.13 | cap 32, 400^3 |
| per-regrid cost | MFC ~12 s constant; AMReX ~0.027 s | cap 32; MFC side is phase-measured |
| GPU busy, AMR arm | 15.1% at cap 64 (8.8% at cap 32) | uniform arm ~80% |

Consequences of the corrected picture:

- **The gap to AMReX is a factor ~2, not ~8.** That changes the tone of the campaign from "something
  is catastrophically broken" to "one structural factor of two remains, and we know where it lives."
- **The per-cell penalty of AMR is already gone at cap 64** (a = 3.50 vs uniform 3.63). What remains
  is entirely per-box/per-launch overhead plus the regrid path.
- **Those two rows are the SAME physical quantity at different caps — reconcile them, do not read
  them as a contradiction.** MFC's AMR arithmetic penalty is 1.90x at cap 32 and 0.96x at cap 64
  (`a` is kernel time per cell: 0.944 s/step / 3.50 ns = 269.7M cells = 64M base + 205.9M
  `fine_work`, an exact match, so `a` is kernel-based and directly comparable to the uniform arm's
  3.633). The arithmetic penalty is a pure block-size effect and it vanishes at the cap we would
  ship. Any decomposition quoting `tax = 1.90x arithmetic x N idle` is a cap-32 statement.
- **The tax metric structurally rewards heavier physics.** The matched case runs weno1 + LF, the
  cheapest numerics MFC has, which maximizes the reported tax. Production numerics lower it while
  changing nothing. Always quote both.
- The regrid-frequency result suggests the practical gap at realistic operating points may already be
  small — but that comparison is confounded (MFC's mesh lags at int=20, doing 2.78x less fine work),
  so it is a hypothesis, not a claim. Experiment E1 in section 7 is designed to settle it.

---

## 3. The phase budget at the shipping cap, and what it does to the priorities

Measured 2026-08-18, exclusive node, 400^3, np=8, 2 levels, `regrid_int=2`, 40 steps
(`amr-bench/logs/rgbuild-0818_0914`). Cap 32 wall 656.975 s; cap 64 wall 293.730 s.

| phase | cap 32 | cap 64 |
|---|---|---|
| **regrid** | 35.2% | **42.2%** |
| gather + reflux + seam (per-box MPI) | 29.0% | 25.1% |
| **rhs** (fine-block advance) | 23.1% | **18.7%** |
| coarse (monolithic L0 RHS) | 2.6% | 4.9% |
| halo / gfill / rk | 2.4% | 2.5% |
| unbracketed | 7.6% | **6.7%** |

**(a) Regrid's share RISES with the cap.** Wall falls 2.24x going 32 -> 64 but regrid falls only
1.87x, so its share goes 35.2% -> 42.2%. Raising the cap remains the best single change made in this
campaign, and it makes regrid relatively MORE dominant, not less.

**(b) The `rhs` bracket is only about half GPU-busy, which bounds what batching can return.**
Total kernel at cap 64 is 0.944 s/step of 6.238 (15.1%); the `s_compute_rhs` tree is 88.7-91.1% of
kernel time; the coarse call is its own `PH_COARSE` bracket and accounts for ~0.233 s/step (the whole
uniform 400^3 arm). So the fine `compute_rhs` inside `PH_RHS` is ~0.61-0.63 s/step against a
1.374 s/step bracket:

| | value |
|---|---|
| `PH_RHS` internal GPU-busy | **44-46%** |
| removable overhead inside `PH_RHS` | **10.2-10.5% of wall** |
| ceiling, perfect batching of `rhs` | **~1.11x** |
| ceiling, + all unbracketed time | **~1.20x** |
| for comparison: regrid alone | 1.73x |

Sanity check that was not tuned to come out right: the same method gives `PH_COARSE` 65% busy
(0.233 of 0.358 s/step), a sensible figure for a monolithic 400^3 advance.

**CAVEAT, and it is why section 7 opens with a measurement.** This composes kernel time from the
4-cap clean-instrument run (wall 6.238 s/step) with phase seconds from this run (7.343 s/step),
violating the action plan's own "one clock per ratio" rule. The first draft of this document got
the same calculation wrong by 2x by charging the coarse kernel to the wrong bracket — and a second
audit found the CORRECTION itself used the uniform arm's WALL (0.2325 s/step) as the coarse KERNEL,
when kernel measurements span 0.187-0.228 s/step across instruments, so the defensible band is
**44-49% busy**, and the "PH_COARSE 65% busy" sanity figure is really 53-65%. A third weakness: the
88.7-91.1% kernel share was profiled on a DIFFERENT case (the 256^3 flat-store campaign), so this
estimate stacks two cross-run compositions. Three audits, three layers of the same defect. **Treat
~45-50% as an estimate that must be replaced by a direct measurement, not as a result.**

**(c) `PH_RHS` does not contain the thing batching exists to remove.** In `m_amr.fpp`,
`s_amr_swap_to_fine()` is called BEFORE `s_phase_tic(PH_RHS)` while its partner
`s_amr_restore_coarse()` is INSIDE the bracket — so the bracket is charged half a swap pair, and the
per-block grid-state reconfiguration lands in the 6.7% unbracketed remainder. The Phase-3 `do islot`
advance loop in `m_time_steppers.fpp` has no enclosing bracket at all. The time stepper's own
comment at Phase 4 names that swap/restore round trip as what blocks batching. **The single largest
claimed benefit of the batched design has never been measured, and cannot be read off any existing
budget.**

**(d) The operating point cuts both ways — do not over-read this table.** These runs use
`regrid_int=2`, which this document argues elsewhere is unrealistic, and `weno_order=1` +
Lax-Friedrichs, the cheapest numerics MFC has (chosen to match AMReX's linear advection). Both
choices inflate the non-physics shares: at `regrid_int=20` regrid roughly halves and `rhs` roughly
doubles, and production WENO5+HLLC multiplies kernel work per cell, raising `rhs` further. **So this
table justifies neither "regrid first" nor "batching first" on its own** — the ordering genuinely
flips between the AMReX-comparison point and the production point, and no budget exists at the
production point. That is E0 and E1 in section 7.

---

## 4. What is established, with trust grades

**Solid (multi-instrument, counts-based, or arithmetically closed):**

1. Uniform MFC is healthy: ~80% GPU busy, rank-invariant launch counts. The tax is AMR machinery.
2. The idle closes arithmetically: dead-time-per-launch x launches ~= the whole gap. Nothing large
   is hiding.
3. The host is the saturated resource, not the network and not the GPU: perf shows ~80% host CPU
   busy (~51% Open MPI progress spinning + ~20% HSA runtime at np=8; pure HSA spinning at np=1).
4. The mapped-entity law (controlled MWE, R^2 = 1.0000, cross-validated in situ): each private array
   or assumed-shape dummy costs ~2 copies and ~31 us per launch; scalars, module-direct arrays and
   explicit-shape dummies are free. Copy-count predictions transfer to MFC to within 0.4%.
5. Every intervention that moved wall reduced toll x count: the cap change, the flux deletion, the
   two `map(alloc:)` clauses. Every intervention that did not, didn't.
6. The batched-kernel MWE: one kernel over all blocks = 433.5 -> 35.0 ms (**12.4x**), operations
   15,370 -> 250. The work-bound crossover is ~500 us of kernel work per region; MFC's mean AMR
   kernel is ~109 us, so **batch >= 8 blocks pushes the mean kernel past the crossover.** Batch 8 is
   not arbitrary; it is the physics of the toll.

**Established negatives (dead levers — do not re-litigate without new evidence):**

- `nowait` chains (213 vs 202 ms, no pipelining) and cross-chain concurrency (disproved in the
  dedicated MWE).
- All seven mapping-clause mechanisms (flat dummies, pointer dummies, present clauses, declare
  mapper, module-scope scalar_field, ...); module scope for derived types is 3.5x WORSE.
- Twelve runtime knob configurations (the harmful ones registered at 1.8-3.2x, which is what makes
  the nulls on the rest meaningful).
- Byte reduction (removing 12% of transferred bytes bought zero wall) and Waitall reduction alone
  (removing 53.8% of Waitalls bought zero wall).
- Load balancing (oracle ceiling 5.1% of wall; regrid imbalance 1.003).

**Suspect — do not build on these without re-measuring:**

- Anything measured at cap 32 (4.9x more boxes than we would ship) or quoted from shared-node runs.
- Absolute s/step across sessions (10% level shifts between option-hash builds; ~12% single-run
  noise floor on this node).
- **"Caps 96 and 128 OOM."** Repeated in `amr_action_plan.md` as settled, but there is NO
  from-scratch run at either cap anywhere in `amr-bench/logs` — the claim descends from the same
  sweep table (`amr_block_batching.md`) that also lists cap 64 as OOM, and cap 64 demonstrably runs.
  It may still be true for an independent reason: per-block solver scratch scales as cap^3, and
  cap 64 already sits at 43.1 of 64 GiB per GCD, leaving only ~1.48x headroom. Worth ONE
  from-scratch check, tempered by two countervailing signals — `b` (ms/box) rises with the cap, and
  ~66 boxes over 8 ranks is poor granularity.
- Anything derived by composing kernel time from one run with phase seconds from another (section 3b
  is explicitly flagged; the first draft of this document got exactly this wrong by 2x).
- The +24.8% widened-bridge penalty (n=3, arms converging — a settling transient, not a number).
- Any MFC-vs-AMReX comparison at regrid_int=20 (mesh-lag confound, 2.78x work difference).

---

## 5. Why it "should have been straightforward" and was not

Worth recording, because the feeling of no progress is itself diagnostic:

1. **This is close to the hardest profiling regime there is.** Four abstraction layers (Fortran
   runtime -> libomptarget -> HSA -> hardware queue), where every instrument distorts the layer
   below it by 1.3-2.4x; run-to-run node noise up to 2.3x shared / 12% exclusive; and three distinct
   transients (first-touch mapping decay, device-pool growth, mesh lag) that make any short window
   lie. In every timing retraction of the campaign, the paired runs agreed in sign — noise
   masquerading as corroboration.
2. **The implicit model "find THE bottleneck" was wrong for this failure mode.** In a latency-bound
   serial chain with several comparable toll terms, removing any one term produces a null, which
   reads as "wrong hypothesis" when it actually means "right hypothesis, wrong granularity." The
   nulls were evidence, correctly interpreted only in aggregate.
3. **Operating-point drift.** The benchmark shipped regrid_int=2, making regrid look like the story
   (37% of wall); at realistic intervals the advance dominates. Different sessions optimized
   different denominators and their numbers contradicted each other while all being true. Hence the
   standing rule: every phase budget and tax figure carries its operating point (cap, regrid_int,
   np, fixes present).
4. **A large fraction of the fortnight went into building the metrology** (counts-first gates, phase
   brackets at 0.1% overhead, settled-tail protocol, from-scratch discipline, apparatus controls).
   That was the price of admission, and it is paid now.

---

## 6. Hypothesis slate going forward

**H1 (dominant, effectively the consensus of all instruments): host-issue-bound per-block advance.**
Fix = fewer, larger units of GPU work per step: the level-batched advance. Design gates are passed
(BC escape hatch exists and is honored; slots are uniform-shape; coordinates reduce to differentials
on Cartesian grids), the widened bridge is built, and the MWE prices the mechanism at 12.4x. The
batch also amortizes three watch-list costs for free: per-swap GPU_UPDATE metadata syncs, bridge
staging (2 full block copies per block per stage today), and the in-order device queue's
serialization of tiny operations.

**But H1 now carries a bound it did not have.** Section 3 puts the ceiling at ~1.11x from the `rhs`
bracket and ~1.20x including all unbracketed time, at cap 64 / `regrid_int=2` — well short of what
"the main event" implies, and far short of the 12.4x the MWE returns. Two things could raise it and
neither is measured: a realistic `regrid_int` roughly doubles `rhs`'s share, and the unbracketed
per-block swap (section 3c) is unknown in size. **The 12.4x MWE figure is an UPPER BOUND** — its
kernel ran over one plain contiguous array with no bridge staging, whereas MFC's batched path must
still pay `s_amr_br_load_all`/`s_amr_br_store_all`.

**H2: per-box blocking MPI is the second wall, currently hidden.** The Waitall-hoist null says the
comms are not rate-limiting NOW — because the host is busy anyway. Once batching frees the host, the
per-box gather/reflux/seam chains become the critical path. Plan the level-aggregated exchange
(post-all, one drain per level — the FillPatch pattern) as stage 2, not never. Note it is smaller
than it looked: cap 64 already removed 4.9x of the boxes.

**H3: vendor-lane share — genuinely untested.** The per-operation toll is a property of amdflang +
libomptarget + ROCm as much as of MFC. Nobody has run the AMR arm on another lane (Cray CCE offload,
or OpenACC on an NVIDIA node). If the same code shows a several-fold cheaper launch path elsewhere,
part of the residual tax is a vendor-runtime tax — which changes how much restructuring is worth
versus waiting on ROCm, and matters for a code that gates on four compilers.

**H4: clustering over-cover.** The split finder has no midpoint fallback, so a convex blob returns
its bounding box (~91% over-cover for a sphere) and `amr_cluster_eff` is dead code (0.7/0.9/0.98
give byte-identical box counts). The action plan correctly marks the break-even model VOID; the only
way to price this is the direct experiment — add the fallback, measure wall on a fixed problem.

**H5: the regrid path** is the same disease on a different limb (per-box, host-based end to end,
~12 s per regrid vs AMReX's ~0.027 s), but its share is operating-point dependent (37% of wall at
int=2, 18% at int=20). **Update, same day: the rb:mem vs rb:unpk discrimination landed (cap 32) and
refuted BOTH code-derived candidates** — the per-box allocate/free brackets to 0.002 s and the
unpack to 0.104 s, against a host half of ~64 s (rb:gath 109.4 s, rb:wait 45.7 s, 16.7% of wall,
shares reproducing the previous session exactly). The planned scratch-hoist fix is dead before it
was built — the third time in this campaign that a code-read attribution failed its bracket. The
host cost lives in the still-unbracketed remainder of `s_amr_gather_coarse_patch`: the non-owner
SEND-side packing loops (byte-proportional, which fits the observed near-flat cap scaling), the
per-(box,source) geometry scan, and the send-pool reserve/drain. Next bracket goes there. The cap 64 arm confirms it (rb:mem 0.010 s, rb:unpk 0.026 s vs a ~25 s host remainder; rb:gath 46.5 s = 15.8% of wall) - both caps agree, the refutation is settled.

**H4a: raise the cap again.** Cap 32 -> 64 returned 2.24x for zero code, per-cell efficiency
saturates there, and the "96/128 OOM" claim is unverified (section 4). Cheap to test, plausible to
fail on scratch scaling.

**H6 (watch-list, post-batching):** single in-order device queue (kernel/copy overlap is exactly
0.000 s today), residual mapped entities on the remaining hot regions, MPI progress spinning.

---

## 7. Experiments, in order of information per hour

**E0 — two direct measurements that decide the ordering (hours, do these first).** Section 3's
estimate of what batching can return rests on a cross-run composition, and an earlier draft of that
same estimate was wrong by 2x. Both gaps close cheaply on the current binary:

  1. **Kernel time INSIDE `PH_RHS`.** One rocprofv3 kernel trace, intersected with the bracket. If
     `rhs` is already ~45% busy the batching ceiling is ~1.11x at this operating point; if it is
     launch-dominated as the cap-32-era numbers suggested, E3 is correctly prioritized. Nothing else
     separates those two worlds.
  2. **Bracket the per-block swap.** Add a phase around `s_amr_swap_to_fine()` (and move
     `s_amr_restore_coarse()` out of `PH_RHS`, or bracket it separately, so the pair is symmetric —
     today the bracket is charged half a swap). This is the batched design's single largest claimed
     benefit and it has never been measured.

**E1 — the one missing baseline (highest priority among the campaign questions, ~half a day).** Cap 64, realistic regrid
interval, landed fixes, both codes, one table — with the mesh-lag confound closed by scaling the
error buffer (`amr_buf`) with the interval so fine coverage matches. **Make matched `fine_work` an
iterated GATE, not a reported control, and report box count beside it**: `amr_buf` does double duty,
since `m_amr_regrid.fpp` uses `thr = buff_size + 2*amr_buf` as the min-separation MERGE threshold, so
raising it fuses boxes as well as growing coverage. Two arms can match on `fine_work` and still
differ in box count, which is the quantity the tax is most sensitive to. This is the number the whole campaign is nominally about, and it does not
exist yet. It also tests the cheerful hypothesis that MFC is already near parity at realistic
operating points (at cap 32 the int=20 taxes were MFC 7.17x vs AMReX 8.13x — confounded, but
suggestive).

**E2 — DONE for cap 32 (see H5): both candidates refuted; scratch hoist cancelled.** Follow-up
bracket: send-side pack + geometry scan inside the gather. Price any fix at a realistic regrid
interval before celebrating.

**E3 — the batching pilot (the main STRUCTURAL change; size it with E0 first).** Wire `s_amr_br_load_all` -> tall (m,n,p) -> one
`s_compute_rhs` -> `s_amr_br_store_all` for a batch of 8, Cartesian-only, per-block path retained as
fallback. Prerequisites, in order: run goldens on the widened bridge (never done); re-measure the
bridge penalty with >= 8 interleaved pairs on an exclusive node (the +24.8% must not stand at n=3);
sweep `amr_br_batch` 2/4/8 with the batched path still unused to separate a memory-driven penalty
from a fixed mapping cost. **Note the tree currently carries the widened bridge uncommitted** — it
allocates 8x the bridge (~880 MB at cap 64) while nothing calls the batched path, and goldens have
never been run on it. If E0 defers E3, revert it rather than leave a measured-slower, untested hunk
in place. **Gate the pilot on counts, not wall**: launches/step on the rhs family
must fall ~batch-fold and GPU busy must rise; wall at pilot scale sits below the noise floor and a
wall gate would produce a false negative.

**E4 — np=1 of the matched case (cheap, explicitly missing).** Cleanly splits MPI-progress idle from
local launch idle at the operating point that matters, and sets the expectation for what batching
alone can recover at np=8.

**E5 — split-finder midpoint fallback A/B (direct experiment, no model).** Fixes real dead code
either way; on convex features it may buy a real fraction of the box count; on fronts and sheets it
will matter more.

**E5a — cap 96 from scratch (one run).** Settles the unverified "96/128 OOM" claim (section 4).
Cap 32 -> 64 returned 2.24x for zero code; ghost overhead keeps falling ((132/128)^3 = 1.096 ->
(196/192)^3 = 1.063) and cap 64 used LESS memory than cap 32. Expect it to be marginal: only ~1.48x
memory headroom remains, per-block scratch scales as cap^3, and box granularity gets poor at ~66
boxes over 8 ranks. Gate on the documented signature — exit 134 with
`HSA_STATUS_ERROR_OUT_OF_RESOURCES`, not a hang.

**E6 — cross-backend probe (if an allocation exists).** Same case on CCE offload or an NVIDIA lane;
three numbers only: launches/step, dead time per launch, tax. Bounds the vendor share of the toll.

---

## 8. Process rules worth keeping

The action plan's "Metric discipline" section is the canon; four additions from this review:

- **Stamp every number with its operating point** (cap, regrid_int, np, fixes present, physics).
  At least three of the campaign's internal contradictions were two true measurements of different
  operating points arguing with each other.
- **Prefer count gates to wall gates** for any intervention whose predicted wall effect is under
  ~2x the in-session noise floor. Counts reproduced to ratio 1.000 across runs whose walls differed
  by 52%.
- **Never compose a ratio from two runs, including when reviewing.** Section 3b does it knowingly
  and flags it; an earlier draft did it unknowingly and reported a 2x-wrong ceiling that argued
  against the campaign's own main line. If the numerator and denominator cannot come from one run,
  the honest output is an experiment, not an estimate.
- **A null after removing one term of a sum is information, not refutation.** Bank it as "term N is
  not independently rate-limiting" and keep the structural hypothesis alive until the structure
  itself has been changed once. The campaign's one structural change at MWE scale (batching)
  delivered 12.4x; none of the term-removals delivered more than 1.36x.
