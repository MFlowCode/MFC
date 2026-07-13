# Banked AMR increments — design + plan (for fresh-session execution)

Three independent, conservation-critical AMR follow-ups, deferred from the multi-level-IB
session (2026-07-13) to be executed **fresh, one at a time** (silent-wrong-answer risk +
context fatigue). Each is its own brainstorm-refine → SDD cycle. **Execution order: #29 →
#30a → #30b** (self-contained first; ref_ratio=4 broadest last). Branch `amr-multilevel`
(PR #6). Each currently **fail-closed** in `m_checker.fpp`.

Ground rules (all three): golden-file validation on CPU **and** GPU (V100, targeted
`--only <UUIDs> -- -b mpirun`, NOT the full-suite sbatch which flakes on SIGILL); existing
goldens byte-identical where the feature is inactive; `./mfc.sh format` → build → precheck-on-commit.

---

## #29 — Distributed pb/mv coupling for non-polytropic QBMM at np≥2

**Problem.** Non-polytropic QBMM carries a per-cell quadrature side-state `pb` (bubble
pressure) and `mv` (vapor mass), `nnode × nb` per cell. It **evolves cell-locally** (no
face flux), but its AMR coarse↔fine coupling is done LOCALLY, exact only at np=1. At np≥2
the SFC block owner needn't hold the block's coarse cells, so pb/mv couple to the wrong
rank's coarse side-state = silent wrong answer. Gated: `m_checker.fpp:106`
`@:PROHIBIT(qbmm .and. .not. polytropic .and. num_procs > 1, ...)`.

**Current local sites (`m_amr.fpp`):**
- Prolong (coarse→fine): `s_amr_prolong_pbmv` (~1884), called ~1347. Reads coarse pb/mv locally.
- Restrict (fine→coarse fold-back): `s_restrict_pbmv` (~1884 def; called 1509/1590). Local; comment 1589 "np>=2 QBMM fold-back is not yet distributed."

**Approach: mirror the `q_cons` P2P distribution for pb/mv (no reflux — pb/mv aren't fluxed).**
The pattern to copy: `s_amr_gather_coarse_patch` (the coarse-patch P2P gather:
`f_amr_rank_coarse_range` + `MPI_IRECV`/`MPI_ISEND`, owner assembles `amr_cg`) and the
restrict scatter in `s_restrict_fine_to_coarse` (owner restricts, scatters covered coarse
slices to coarse-cell owners). Two differences from q_cons: (a) per-cell payload is
`nnode*nb` (both pb and mv), so buffer sizes scale by that; (b) NO Berger-Colella reflux
(pb/mv have no C/F flux correction — prolong sets the fine ghost/interior, restrict folds
back; that's the whole coupling).

**Task plan (SDD):**
1. **Distributed pb/mv gather** — assemble the block's coarse pb/mv patch on the owner via
   P2P (mirror `s_amr_gather_coarse_patch`; a parallel `amr_cg`-like pb/mv buffer, or extend
   the gather to carry pb/mv alongside q_cons). Prolong reads that instead of local coarse.
   Gate: np=1 byte-identical (single owner → local copy path unchanged); build.
2. **Distributed pb/mv restrict scatter** — owner restricts fine pb/mv → coarse averages,
   scatters covered coarse slices to coarse-cell owners (mirror the q_cons restrict P2P).
   Gate: np=1 byte-identical.
3. **Lift the gate + golden** — drop the `num_procs > 1` term from `m_checker.fpp:106`; add a
   np=2 non-polytropic-QBMM + AMR golden. Validate: per-node pb/mv moments conserved
   (machine-zero for the conservative moments; pb/mv are a side-state so define the gate as
   "np=2 trajectory == np=1 trajectory" on a bit-uniform grid, mirroring the q_cons np-cross
   check). CPU + GPU.

**Open questions for the fresh brainstorm:** whether to fold pb/mv into the existing q_cons
gather/scatter buffers (one MPI round) or a separate exchange (simpler, more messages);
device-buffer handling for the larger payload (pb/mv are `pres_field`, GPU_DECLARE'd).
**Validation oracle:** np=2 == np=1 trajectory on a bit-uniform grid (the WENO-table-ulp
finding means non-bit-uniform grids diverge at ulp).

---

## #30a — Lazy owned-only IB marker sizing

**Problem.** The multi-level-IB increment sized the declare-target `ib_markers` (and the park
slots) to `2**amr_max_level * base_block_extent` — the **global** deepest extent, right at
np=1. At np≥2 (once multi-level IB is un-gated there — a separate future item) that
over-allocates the never-realloc'd device field on every rank to the global deepest, even
ranks owning only shallow/no fine blocks. Task #30's "lazy owned-only sizing" = size the
marker field to the deepest fine block **a given rank actually owns**, allocated lazily.

**Current state:** `s_ibm_marker_bounds` (m_ibm.fpp, added 2026-07-13) computes the deepest
bound from `amr_max_level` + `amr_block_beg/end` (global). Multi-level IB is currently np=1
only, so this is not yet a live cost — **#30a is a memory optimization that pairs with
un-gating multi-level IB at np>1** (itself deferred). Low priority until np>1 IB lands.

**Approach.** Replace the global `2**amr_max_level` bound with the per-rank owned deepest
level. Because `ib_markers` is a device declare-target that must NOT be reallocated after
mapping, "lazy" means: at init, size to the deepest level the rank's INITIAL decomposition
could own; if a later regrid would need deeper, that's the same never-realloc constraint —
so either (a) size to the rank's static owned-region deepest possible, or (b) accept the
global bound at np=1 (status quo) and only optimize once np>1 IB + repartition-on-restart
lands. **Recommend deferring #30a until np>1 multi-level IB exists** — optimizing an
allocation for a configuration the checker doesn't yet admit is speculative.

**Task plan:** small — parameterize `s_ibm_marker_bounds` by a per-rank owned-level input;
validate memory footprint drops at np>1 with no golden change. Blocked on np>1 IB.

---

## #30b — ref_ratio = 4

**Problem.** Refinement ratio is hard-coded 2:1 per level. One 4:1 level reaches the
resolution of two 2:1 levels with one fewer coupling layer. The `2*`/`2**level` extent
factors are threaded through every coupling kernel.

**Approach.** Introduce a runtime `ref_ratio ∈ {2,4}` (param, default 2 = byte-identical).
De-hardcode the fine-extent, prolong, restrict, halo, and reflux stencils from `2` to
`ref_ratio` (and `2**level` to `ref_ratio**level`). The prolong stencil widens (4:1 injection
+ the multi-fluid/species closure over a 4-wide child block); restrict averages `ref_ratio^d`
children; the fine-fine seam halo and reflux child-face counts scale by `ref_ratio`.

**Highest-risk feature of the three** — it touches every coupling kernel's stencil, and a
wrong factor conserves-to-machine-zero while being physically wrong.

**Task plan (SDD):**
1. **Param + gate** — add `ref_ratio` (definitions.py + descriptions.py + checker default 2);
   gate `ref_ratio ∉ {2,4}` and any unsupported combos fail-closed. Byte-identical at
   ref_ratio=2. `amr_slots(:)%ref_ratio` already exists (per-block) — audit that it's
   populated from the param, not a `2` literal.
2. **De-hardcode extents** — sweep `2*(...)`/`2**level` → `ref_ratio*(...)`/`ref_ratio**level`
   in fine-geometry sizing, marker sizing, `old_ext`, tower weight, `s_amr_fine_fine_halo`
   (the level-aware `2**level` from #35), reflux child-face counts. Gate: ref_ratio=2
   byte-identical (every site reduces to `2`).
3. **Prolong/restrict stencils** — generalize `s_prolong_one_var` (injection/interp for a
   `ref_ratio`-wide child), `s_restrict_one_var` (`ref_ratio^d` child sum), and the alpha/
   species closures. Gate: ref_ratio=2 byte-identical.
4. **ref_ratio=4 golden** — a single-level ref_ratio=4 case; validate conservation
   machine-zero (inviscid) + a resolution check vs two 2:1 levels. CPU + GPU.

**Open questions:** buff_size / stencil reach for a 4-wide prolong (does the coarse patch
margin `amr_cpat_mar` suffice?); interaction with multi-level (ref_ratio=4 AND amr_max_level>1
= 16:1 two levels — likely gate to one at a time first); whether `ref_ratio` is global or
per-level.

---

## Cross-cutting notes
- All three keep existing goldens byte-identical when inactive (param defaults / gated).
- GPU validate on V100 via targeted `-- -b mpirun` (the full-suite sbatch SIGILL-flakes on
  arch-mismatched nodes; see the `phoenix-srun-mpirun` memory + the hpcx-bin-on-PATH note).
- Conservation gates: q_cons machine-zero; pb/mv and any non-conservative side-state use the
  "np=2 == np=1 trajectory on a bit-uniform grid" oracle.
