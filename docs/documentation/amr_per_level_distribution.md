# AMR per-level distribution

Design for the next phase of block-structured AMR: distribute each refinement level
independently, so AMR strong-scales instead of degrading as ranks are added. Written
2026-07-29 against `up/mega` @ `7b1e4933`. Companion to @ref amr_block_batching, which covers
the per-block launch cost, and @ref amr_multilevel, which covers the nesting.

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

## What is missing

**P2P parent<->child coupling — Track 2.** Gather, restrict and flux-register delivery are
np=1-local (future work #27; `m_amr.fpp:629, 967, 988, 1006, 2350, 2372, 4354`). Multi-rank works
today only because a refinement tower co-locates on its level-1 anchor, which makes a whole tower
the balancer's smallest atom (weight `cost * rr^(l*d)`, `m_amr.fpp:1673-1679`) and caps
granularity at depth.

Under per-level distribution, tower co-location is not a constraint to be relaxed separately — it
disappears as a consequence. Track 2 therefore stops being one item on a queue and becomes the
spine of this phase.

## Sequencing

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

## Validation

The bar this branch already uses: np=1 bit-identical, np>=2 conservation-exact, plus the goldens'
tolerance compare, and Frontier CCE as the cross-compiler gate — this work touches device
residency of grid state, where every regression so far has been CCE-only.

For steps 3 and 4, prefer the **conservation ladder** over field diffs: authoritative mass after
each stage, the pure reflux delta, coarse-equivalent fine mass per block, and covered-cell mass
per block, all allreduced and printed from both arms. That is what localised the np>=2 restrict
bug (`ab87d49e`) when field diffs could not; many-to-many coupling is exactly where that class of
bug hides.

Note that `./mfc.sh test --only AMR` does not match the coexist goldens — their trace token is
"AMR + L0 tiles". Run `1F074C5D 8D466A94 83CC5C6D 33060D84 FD056B71 98AA6EDB D99F85F8 09E0D257
93EFC4F6` by UUID as well, or they silently go untested.
