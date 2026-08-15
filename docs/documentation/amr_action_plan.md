# AMR performance: what to change, and what not to (2026-08-15)

Companion to `amr_block_batching.md`, which is the chronological research log. This file is the
**action list**: what to change in the code, in priority order, with the measurement behind each and
an explicit list of things NOT to do. It supersedes the three-track plan appended to
`amr_block_batching.md` earlier the same day — that plan was built on a configuration (cap 32) now
known to be 2.74x off the achievable one, and on a headline gap (7.64x) now known to be an
unmatched comparison.

---

## The numbers this rests on

Matched case: 400^3, np=8, 2 levels, ref_ratio 2, volumetric Gaussian blob, `weno_order=1` +
Lax-Friedrichs (chosen to sit as close to AMReX's linear advection as MFC can, which is the
conservative choice for a tax ratio - see "Metric discipline" below).

**Tax at MATCHED block size** (each code vs its OWN uniform baseline):

| cap | MFC | AMReX | MFC excess |
|---|---|---|---|
| 32 | 24.24x | 5.84x | **4.15x** |
| 64 | 6.89x | 3.40x | **2.03x** |

The previously reported "MFC 23.92x vs AMReX 3.13x = 7.64x excess" compared **MFC at cap 32 against
AMReX at cap 64** - a `sed` override in the harness of an inputs file whose committed value was 32.
At matched cap the excess is 2.0-4.2x. AMReX's tax is **not** flat in the cap either (5.84 -> 3.40),
so per-box cost is not unique to MFC; ours is simply larger.

**Cost model `wall = a*cells + b*boxes`: REJECTED, and now REPLACED BY DIRECT MEASUREMENT.**

The fit failed its own residual test (spread 30.8% at 4 caps vs a ~12% noise floor), and the direct
disproof needed no fit: cap 40 -> cap 48 has BYTE-IDENTICAL cells, 11% fewer boxes, 31% less wall.
Both break-even numbers derived from it (113k and 573k cells/box) are VOID.

**MEASURED decomposition (clean protocol - see "How this is measured" below):**

| | cap 32 | cap 64 |
|---|---|---|
| step wall | 14.485 s | **6.238 s (2.32x faster)** |
| kernel (trace union, differenced) | 1.278 s | 0.944 s |
| **a** | **6.09 ns/cell** | **3.50 ns/cell** |
| **b** | **6.29 ms/box** | **15.24 ms/box** |
| residual | 5.613 s (39%) | 1.880 s (30%) |
| boxes / cells | 1207 / 209.7M | 224 / 269.6M |

**`a` depends on block size - 1.74x across ONE cap step**, which is what broke the fit. And cap 64's
a = 3.50 is essentially the UNIFORM arm's 3.63: **at 132^3 blocks the per-cell efficiency penalty of
AMR vanishes.** At 68^3 it is 1.74x. So bigger blocks win on kernel efficiency, not only box count.

**`b` moves the WRONG WAY - 6.29 -> 15.24 ms/box.** Per-box cost is NOT a fixed toll; it scales with
block volume, because packing, prolongation and ghost fills inside those phases are per-CELL work.
Total per-box time still falls (7.594 -> 3.414 s/step) since 5.4x fewer boxes beats 2.4x more each.
But "a fixed ~10 ms per box" was wrong, and the fitted model was wrong in BOTH terms.

**The residual (30-39%)** is mostly phases excluded from the per-box bucket - chiefly `rhs`'s
non-kernel time (25.5% of wall, and rhs is only ~5% kernel). The FULL top-level phase set covers
91.9% of cap-32 wall with an 8.1% residual, so this is launch overhead inside the physics (Tier 2.1),
not an unknown.

**HOW THIS IS MEASURED (do not regress to the old way):**
- wall = `ph_wall_total` ("step-loop wall", needs rank_time_wrt=T). It brackets only
  s_perform_time_step, INCLUDES regrid, and excludes pre_process/init/finalize/checkpoint I/O.
  `Total-time` differencing gave **34% spread on byte-identical work** and must not be used.
- difference two FROM-SCRATCH runs. Do NOT warm-start: a restart restores a PARTIAL hierarchy and the
  first regrid rebuilds it (level-2 boxes 342 -> 760 -> 792), moving preprocessing INTO the window.
- sum only TOP-LEVEL phases; the nested rg*/rb*/gwait brackets are inside their parents and
  double-count (the printed RESIDUAL goes to -55%).
- kernel time by UNION of intervals, differenced between the two step counts so the transient cancels.
- Runner: `amr-bench/expt/measure_clean.sh`. Reader: `amr-bench/fit_cost.py`.

---

## TIER 0 - configuration, no code

### 0.1 Raise `amr_max_grid_size` 32 -> 64
**2.32x wall on the same physical problem** (14.485 -> 6.238 s/step, clean instrument), 5.4x fewer
boxes (1207 -> 224), and **lower** peak memory (50.1 -> 43.1 GiB/GCD). Caps 96 and 128 genuinely OOM,
so **64 is the maximum - which is exactly AMReX's block size.**
The earlier "2.74x" came from Total-time differencing (34% noise) and is superseded.
Bigger blocks win TWICE: 5.4x fewer boxes AND per-cell efficiency recovering to the uniform arm's
(a 6.09 -> 3.50 vs uniform 3.63).

This was blocked for weeks by a recorded claim that the cap was exhausted and larger caps OOM. That
claim was a **checkpoint-restart confound**: restarting each cap from a cap-32 checkpoint freezes the
box count while enlarging every slot, so init allocated cap-64-sized slots for cap-32's 1250 boxes -
a combination that exists in no real run. From scratch, caps 40/48/64 all run.

**Gate before changing the default:** goldens, and settle whether cap 32 is *under-refining* level 2
rather than cap 64 over-refining it. `fine_work` goes 139.6M -> 205.9M (+47%) with the cap, level-1
coverage looks cap-independent (343 = 7^3 and 64 = 4^3 tile the same ~200-coarse-cell extent), and
the nesting margin (`amr_cpat_mar` = 2 coarse cells at buff_size 2) predicts only 1.23x of a >=1.47x
effect - **so the mechanism is NOT yet pinned**. If cap 32 under-refines, the two caps are not
solving the same problem and part of the 2.74x is accuracy, not speed.

---

## TIER 1 - small, local, high confidence

All four are verified in source, none touches the solver, and all attack `b` (per-box overhead).
Sizes below come from cap-32 measurements and **must be re-measured at cap 64** before scoping.

### 1.1 Skip ranks with no overlapping coarse data in the gather
`s_amr_fine_stage_fill` calls `s_amr_gather_coarse_patch` **before** the
`if (.not. amr_rank_owns_block) return`, so every rank enters it for all ~1250 boxes.

**This is not simply "8x redundant"** - non-owners legitimately `MPI_ISEND` coarse data they hold, so
they must participate. The correct early-out tests the **cached `amr_ovl_gather(:, amr_cur)` overlap
list** (already O(1), replicated): a rank with no overlap for this box can skip the geometry, the
allocations and the message entirely. Only ranks that actually contribute should do work.

### 1.2 Hoist the per-call heap allocations; size `rbuf` correctly
`m_amr.fpp:910` allocates `rbuf(maxsz, nsrc), reqs, srank` **per gather call** (~625 calls/step), and
`maxsz` is the **whole patch** while each source contributes only its intersection box (~4x
over-allocated). Hoist to a module-scope pool sized once per regrid; size at `boxsz`.

### 1.3 Use GPU-aware MPI in the AMR gather
`use_device_addr` appears **0 times** in `m_amr.fpp`, while the base halo already wraps its buffers in
`GPU_HOST_DATA(use_device_addr=...)` under `use_rdma_transport` (`m_mpi_common.fpp:860`). Every AMR box
currently goes device -> host (`copyout`) -> MPI -> host -> device (`copyin`). All 8 ranks are on one
node over XGMI; device-to-device is the natural path. Per the mapped-entity law each `copyout`/`copyin`
on a target region is itself a mapped array costing ~2 copies + ~27 HSA calls per launch.

### 1.4 Cache the coarse-restore grid sync
`s_amr_swap_to_fine` / `s_amr_restore_coarse` each end in `s_amr_sync_grid_state_to_device`: four
`GPU_UPDATE(device=)` statements covering 14 entities, run **per block per RK stage** (~470 pairs per
rank-step). The restore re-uploads **byte-identical coarse coordinates** every time. Cache the
last-synced slot id and skip. Better but larger: hoist the restore out of the per-block loop entirely
(swap b1 -> advance -> swap b2 -> ... -> restore once); `m_time_steppers.fpp:592` already identifies
the swap/restore interleaving as what blocks batching.

---

## TIER 2 - larger, and one is not an AMR problem at all

### 2.1 `convert_conservative_to_primitive` runs at 1.4% of memory bandwidth
34.03 ms median for ~8M cells/rank - **11x slower than structurally similar kernels in the same code**
(riemann/weno/tvd_rk all land near 15% of peak). It is **54.6% of all uniform-arm GPU time** and 33.2%
of the AMR arm's.

The ratio between MFC's own kernels needs no bandwidth assumption: same case, same cells, same rank.
Suspected register spill - the kernel privatizes **16 variables including 5 arrays** per thread. Get
`-Rpass-analysis=kernel-resource-usage` (or the ROCm equivalent) and confirm spilling before changing
anything.

**This lives in `src/common/`, runs in all three executables, and is independent of AMR.** It is
probably the single largest non-AMR win available, and it partially undermines the claim that the
27.6x per-cell gap vs AMReX's uniform arm "is physics, not a defect" - some of it is inefficiency.

### 2.2 Per-level batching of the gather (the AMReX pattern)
AMReX's `FillPatchTwoLevels` is called **2.23x/step** and services a whole level in one
`ParallelCopy`; MFC's gather is called **~625x/step**, one box at a time, each with its own
`MPI_WAITALL`. Call-count ratio **281x**.

This is the structural fix, but it is now **much smaller than it looked**: at cap 64 there are 4.9x
fewer boxes, so most of what batching would remove is already gone. Tiers 0+1 and this one overlap
heavily and **must not be multiplied**.

---

## DO NOT DO

- **Load balancing.** 0.98, matching AMReX's knapsack. `[amr-balance] max/mean` 1.003-1.009 and regrid
  imbalance 1.003. Ruled out by measurement.
- **[UNDECIDED - the number justifying this is VOID] Tighter clustering / adding the missing
  midpoint fallback to `s_amr_find_split`.**
  The defect is real - the split finder tries a zero-signature hole, then a Laplacian inflection, and
  gives up (`ok=.false.`) with **no midpoint fallback**, so a convex blob returns its bounding box
  (~91% over-cover for a sphere) and `amr_cluster_eff` is unreachable dead code (verified: 0.7 / 0.9 /
  0.98 give byte-identical box counts). **The break-even argument that justified NOT fixing this is UNSOUND as of 2026-08-15.** `b/a` came
  from one arbitrary pair of caps; a second pair from the same sweep disagrees badly:

  | fitted from | a (ns/cell) | b (ms/box) | break-even |
  |---|---|---|---|
  | cap 32 + 64 | 16.65 | 9.54 | **573k cells/box** |
  | cap 32 + 40 | 42.34 | 4.78 | **113k cells/box** |

  `a` differs 2.5x, `b` 2x. Both fits are exact by construction (2 points, 2 unknowns, zero degrees
  of freedom) so neither is testable - and **they STRADDLE the decision boundary**: a cap-32 block is
  262k cells, BELOW 573k (do not split) but ABOVE 113k (splitting could pay). **Do not act on this
  entry in either direction until caps 48 and 64 land.** If `fit_cost.py` then reports residuals
  monotone in cap, the two-term model is missing a term and NEITHER break-even is meaningful.

  Log the defect regardless - it makes `amr_cluster_eff` dead code, and it will matter on non-convex
  features (shock fronts, sheets) where holes exist and the finder does work.
- **Descriptor-copy reduction.** 19.54 copies/launch, but the initiating HSA call is 0.9% of the gap.
- **Reducing kernel/region count for its own sake.** AMReX launches 7.9x more kernels than MFC on one
  case and 5.5x fewer on the matched one; per-launch cost dominates, not count.
- **Anything scoped from cap-32 phase shares.** Regrid 38.0% / gather 12.4% / reflux 11.4% / seam 3.7%
  were all measured with 4.9x more boxes than we would ship.

---

## Metric discipline (adopt these; they are why the numbers above are trustworthy)

- **Perturb and regress; never decompose a residual.** Every retracted claim in this campaign was a
  residual (`idle = wall - kernel`, `dead/launch = (wall - kernel)/launches`) that absorbed whatever
  it did not understand and could never be falsified. `a` and `b` come from independent perturbations.
- **A tax figure must carry its arithmetic intensity.** `tax = 1 + Overhead/(a*cells)`, so the metric
  *structurally rewards doing more work per cell*. Switching to weno5+HLLC would lower our tax while
  changing nothing about the overhead. Quote weno1+LF for the AMReX comparison (conservative), and
  weno5+HLLC separately for production relevance.
- **Judge cap changes by WALL on a fixed physical problem, not ns/cell.** A cap that over-covers
  inflates its own denominator: cap 64 looks 3.5x better on ns/cell but is 2.74x better on wall.
- **One clock per ratio.** GPU-busy computed as traced-kernel-time over untraced wall gave rank 4
  **101.2%** - proof the terms came from different runs.
- **Union-of-intervals for kernel busy time, never median x count.** MFC never overlaps kernels
  (sum == union exactly); AMReX does (sum/union = 1.270). Under a consistent estimator the arithmetic
  terms are **MFC 1.90x vs AMReX 1.92x - identical**, not the "1.60x vs 1.96x, our algorithm is
  better" that median x count produced.
- **Never take wall from an instrumented run.** rocprofv3 inflates MFC 1.3-2.4x; TinyProfiler inflates
  AMReX 37%.
- **Single runs on this node have a ~12% noise floor** (three byte-identical-work runs: 645/676/725 s).
  Include a known-null setting in every sweep to measure it in-session.
- **Every experiment needs an apparatus control.** A restart-based A/B silently pins the state the
  parameter is supposed to redetermine; that one mistake cost weeks.
