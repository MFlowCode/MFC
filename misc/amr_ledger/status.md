
# AMR performance status (2026-09-15)

This page is the current truth for the AMR performance program and replaces the ledgers as the place to
look first. `action_plan.md` is the lab notebook (165 dated entries, append-only, with retractions);
`endstate.md` is the architecture (four pillars, invariants W1-W8). Nothing here that disagrees with
them is a typo: this page wins, and it is revised in the same change that lands or retracts a result.

## The three statements

| # | statement | status | evidence |
|---|---|---|---|
| 1 | weak scaling: wall per rank-doubling at fixed work per rank, 8 GPUs/node | **1.28x per doubling np8 -> np32** (1.23x, then 1.32x) against AMReX's 1.23x on matched cross-node rungs; one rep per rung, np32 biased by its own straggler; the growth is the base-grid advance and base halo, not the regrid family | ledger 128, `amr-bench/notes/l128_curve_so_far.txt` |
| 2 | per-GPU AMR overhead: excess of the AMR step over the code's own uniform step, same node and hour, vs AMReX on the same problem | **PARITY. MFC 0.42 s/step vs AMReX 0.38-0.39, 1.08x** on a verified-healthy node (k004-005, A-B-A-A, 3 reps/arm, all arms tenancy- and stall-clean); the parent commit reads 0.42-0.45; the control-to-control spread of 0.03 is the instrument's floor | ledger 165, `amr-bench/logs/twocode-u5-*-419849-*.log` |
| 3 | correctness at scale: no silent NaNs, no configuration accepted but not honoured | complete: keyed tags on every wave family, np=2 oracle and seed controls, 71 AMR goldens | ledgers 56-60 |

Per cell-update on the statement-2 deck (400^3 base, two levels, ratio 2, regrid every 20 steps, RK3, no
subcycling, one ideal gas, advected density blob):

| | uniform ns per cell-update | AMR ns per cell-update | AMR / own uniform |
|---|---|---|---|
| MFC (WENO5 mapped + mp_weno, HLLC, 5-equation model) | 2.16 | 3.84 | 1.78x |
| AMReX-CNS (second-order, 5 variables) | 1.28 | 2.38 | 1.86x |

The AMR machinery inflates MFC's step by less than AMReX's inflates its own. MFC's AMR step is slower in
absolute terms (0.96 vs 0.82 s) because its base solver costs 1.7x per cell, which is a scheme choice.

## What the remaining 0.42 s/step is (mean rank, differenced, ledger 165 run B)

| term | s/step | ceiling if removed | why it stays |
|---|---|---|---|
| fine-solver inflation (ghost planes per block, batch padding, dispatch floor) | ~0.09 | ~0.05 | launch count is spent (ledger 163); private-array tax taken (160); recon window clipped (164) |
| MPI wait, all families | ~0.19 | ~0.08 | ~14 rendezvous/step; the skew is per-segment work on two ranks (161); rebalancing adds launches (159); reordering does nothing (137) |
| AMR work, non-wait (regrid 0.06, restrict/gather/reflux packing 0.05, seam/fill/migration 0.03) | ~0.14 | ~0.05 | wire pools already device-resident (162); regrid is 1.2 s per event, per-event AMReX time never measured |

Realistic total yield 0.05-0.10 s/step against a 0.03 floor. The measurement program for statement 2 is
closed.

## What was tried and rejected, one line each

Store growth ratchet (fixed, 2.32x), per-box rendezvous (waves), blocking level-2 sends (ISEND pool),
level-1 gather collective (deleted), host-staged wire pools (device-resident), private work arrays
(block-local), reconstruction over ghost planes (clipped), batch-leader ordering (neutral), cost-weighted
balancer (worse than none), greedy remapping (worse), reflux post/drain reorder (sign unresolved),
asynchronous kernels via nowait (runs inline on this stack), single/mixed precision for speed (forbidden:
the comparison is at double), kernel-as-device-routine (device-side loss), deferred sends (only where a
downstream sync absorbs the drift), equal-tile clustering and load-balance feedback (instruments, default
off), subcycling (parity at matched fidelity).

## Why the numbers looked unclear for a week

The excess is a difference of two noisy walls (floor +-0.03, 7 % of the answer); the ratio carries the
node (the same code reads 1.08x on k004-005 and 1.5-1.9x on k004-004/006, because MFC's host-bound
overhead is node-sensitive and AMReX's device-bound step is not); the target moved four times across
GOAL versions; and levers were priced on mean rows rather than the heaviest rank's own chain until
ledger 164. Rules that now hold: node-matched control and treatment in one allocation, parent commit as
control, pre-registration before submission, stall detection on every arm, no post-hoc statistic changes.

## What blocks the PR (as of 2026-09-15)

1. The tip did not compile on the GNU CI lanes: three `#ifdef` lines indented by the ledger-160 commit
   (fixed in the same change that adds this page).
2. The four-compiler test suite has not run to completion on the branch since 2026-08-02; every run since
   was cancelled by the next push. The NVHPC and Frontier lanes have to be read before anything else.
3. The branch bundles the AMR with kernel restructurings, a hypoelastic HLLD solver, reactive-burn
   substeps, IB collisions, probes, and load-balance/SFC/active-box modules: 1048 commits, 856 files. It
   needs splitting before review.
4. `m_amr.fpp` carries every experiment alive. Removed 2026-09-15 (commit after ledger 165): the
   `amr_equal_tiles` and `amr_lb_beta` instruments and their goldens, the `amr_batched_gather` pooled path,
   the `[amr-cov]` dead-byte counters, and the `amr_bat_pad` knob (now the measured constant 0.10). Still
   in the tree and each a product decision: the per-block fine advance (it is the AMR path for every physics
   the batched advance excludes: subcycling, stretched or cylindrical grids, MHD, IGR, bubbles,
   hypoelasticity, chemistry, relaxation, surface tension, the 6-equation model; 32 of the 70 goldens run
   it), subcycling (parity with lock-step at matched fidelity, 18 goldens), the level-0 tiling subsystem
   (`l0_ntile`, 26 routines, no production deck uses it, 9 goldens), and the load-balance / SFC /
   active-box modules from before the AMR.

## Finish line (user-gated, in order)

Frontier CCE suite; the constant-density 8-GPU/node ladder with AMReX at matched rungs; the upstream
landing in the constitution's order.
