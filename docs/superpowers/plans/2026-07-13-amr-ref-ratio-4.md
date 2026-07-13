# Runtime ref_ratio ∈ {2,4} for single-level AMR — Implementation Plan (#30b)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax. This is the HIGHEST-RISK banked AMR increment — a wrong refinement factor conserves-to-machine-zero while being physically misplaced. Every task gates on "existing goldens byte-identical at ref_ratio=2 (the default)". Start FRESH with context headroom.

**Goal:** Make the AMR refinement ratio a runtime parameter `ref_ratio ∈ {2,4}` (default 2), so one fine level can reach 4:1 resolution instead of 2:1 — de-hardcoding every 2:1 assumption threaded through the coupling kernels.

**Architecture:** The per-block field `amr_slots%ref_ratio` already exists and is used in ~49 sites, but it is pinned `= 2` at one assignment (`m_amr.fpp:5418`) and many kernels still bake in the literal 2. Drive the field from a new param, then fix every hardcoded-2 site. Restrict to **single-level, global, lock-step** for v1 (see Global Constraints); multi-level 4:1 (16:1 towers) and subcycled 4:1 are deferred follow-ups.

**Tech Stack:** Fortran 2008 + Fypp, OpenACC/OpenMP-offload GPU, `src/simulation/m_amr.fpp` + `m_amr_registers.fpp`, param system (`toolchain/mfc/params/definitions.py`).

## Global Constraints
- **Scope = single-level, global, lock-step, {2,4}.** ref_ratio applies to every level (global); v1 gates `ref_ratio ≠ 2 ⇒ amr_max_level = 1 .and. .not. amr_subcycle`. Multi-level 4:1 and subcycled 4:1 are separate follow-ups.
- **Default ref_ratio = 2 ⇒ every changed site reduces to the current literal.** Existing goldens MUST stay byte-identical after EVERY task. This is the primary per-task gate and the safety net for this high-risk work.
- **Prolong offset is conservation-visible.** Unlike most factors (which conserve-but-misplace when wrong), the prolong child-offset must keep `mean over children = 0` or the fine children don't average to the coarse value → conservation breaks. So the ref_ratio=4 golden's conservation check catches an offset bug. Use this.
- MFC discipline: GPU via `GPU_*` Fypp macros only; precision `wp`/`stp`; abort via `s_mpi_abort`/`@:PROHIBIT`; `@:ALLOCATE`↔`@:DEALLOCATE`; new param via `definitions.py` `_r()`+`_nv()` (+ `case_validator.py` if physics-constrained). Build: `./mfc.sh format` → build → precheck-on-commit.
- **Ops (this env):** builds need the Bash sandbox OFF (else SIGTERM/143); CPU mpirun via `source ./mfc.sh load -c p -m c`; GPU via `load -c p -m g` + `clean --gpu acc` before an incremental GPU build (nvlink "newer than"); do NOT run interactive builds while a batch GPU build shares `build/`. Golden validate `-- -b mpirun`.

---

## The hardcoded-2 inventory (exhaustive sweep, 2026-07-13)
The change point is `amr_slots(islot)%ref_ratio = 2` (`m_amr.fpp:5418`) → `= ref_ratio`. Every site below assumes 2:1 and MUST be generalized. `rr` = `amr_slots(amr_cur)%ref_ratio` (or the local `rr`/`ref_ratio` already in scope).

**Cat 1 — fine-extent sizing (`2*E-1` → `rr*E-1`):** `m_amr.fpp` 265, 267, 268 (`max_f{1,2,3}`), 1289, 1291, 1292 (`%m/%n/%p`), 1599–1601 (L2-fits-parent guard — multi-level, gated, but generalize for correctness). `m_amr_registers.fpp` 147, 149, 150 (freg sizing — MUST match m_amr's max_f).

**Cat 2 — `2**amr_block_level` → `rr**amr_block_level`:** `m_amr.fpp` 1146–1148 (tower weight), 3210 + 3211–3216 (`fmul` fine-fine halo extent), 4580 + 4581–4583 (`old_ext` regrid stash).

**Cat 3 — fine→coarse index `/2` and `maxc /2` → `/rr`:** `m_amr.fpp` 244, 246, 247 (`amr_maxc = (glb+1)/2`), 253 (fit cap allreduce), 459 (`c = lo + fi/2` in `s_build_level_coords`), 2609, 2620, 2632, 2643, 2656, 2667 (`floor(jg/2.)` ghost-coord maps in `s_amr_swap_to_fine`), 2757, 2759, 2760 (`floor(j/2.)` in `s_amr_igr_swap_sigma` — IGR only).

**Cat 4 — prolong child-offset `(mod(f,rr)-0.5)*0.5` → general (see Task 2 for exact form):** `m_amr.fpp` 1363, 1366, 1369 (`s_prolong_one_var`), 1437, 1440, 1443 + 1484, 1487, 1490 (`s_interpolate_coarse_to_fine` alpha + bubble paths), 2871, 2876, 2879 + 2920, 2925, 2928 (`s_amr_fill_fine_ghosts` GPU kernels). **15 sites.** CONSERVATION-CRITICAL.

**Cat 7 — coarse-patch margin:** `m_amr.fpp` 362 (`amr_cpat_mar = (buff_size+1)/2+1`), 1331 (`nmar`, same). → `(buff_size + rr - 1)/rr + 1` (byte-identical at rr=2). The `2*amr_cpat_mar` two-sided additions (364–366, 542–545, 692–695, 896–899) are rr-INDEPENDENT — leave them.

**Cat 8 — fine-grid COORDINATE generation (NOT mechanical — redesign):** `m_amr.fpp` 463 (`mod(fi,2)==0` left/right-half branch in `s_build_level_coords`), 2610, 2621, 2633, 2644, 2657, 2668 (`mod(jg,2)` coarse-midpoint-vs-edge branch in `s_amr_swap_to_fine`). For rr=2 a fine boundary is either the coarse midpoint (even) or coarse edge (odd); for rr general there are `rr` uniform sub-boundaries per coarse cell at fractions `k/rr`. The even/odd branch must become a general `k = mod(idx, rr)` fractional placement `xcb(c-1) + (k/rr)*(xcb(c)-xcb(c-1))`.

**Cat 5/6 — coarse-RHS reflux (`s_amr_apply_reflux`, `m_amr_registers.fpp`):** 549–552, 560–561 (x-faces `nch`/`dd*_hi`/`f*0 = 2*c`), 582–584, 592–593, 596 (y-faces), 614, 622–623, 625–626 (z-faces, `nch=4`, `do dd=0,1`). This routine has its OWN hardcoded-2 loops SEPARATE from the already-general shared kernel `s_amr_reflux_apply_faces` (which takes `rr`). Preferred fix: **route `s_amr_apply_reflux` through the shared kernel** (completes the #39 dedup); fallback: thread `rr` into its loops.

**Subcycle (DEFERRED — v1 gates lock-step):** 182 (`amr_dt_fine = 0.5*dt`), 3426, 3566 (`do sub = 1, 2`), 3428, 3568 (`th = ...*0.5`), 3620 (`dt_sub*0.5`). Out of v1 scope (gated `.not. amr_subcycle`); a follow-up increment generalizes these to `ref_ratio` substeps.

**Multi-level inset (OUT of scope — gated):** 1585–1587, 4539–4541 (`/4` L2 placement heuristic). amr_max_level=1 in v1.

**Confirmed SAFE (do not touch):** the shared reflux kernel `s_amr_reflux_apply_faces` + capture routines (already `rr`-general); `nchild`/child-loops in the restrict + pbmv + tagger paths (already `rr`); Morton `iand` (SFC); the `0.5*(xcb+xcb)` cell-center formulas in coord gen (mesh-spacing, not ratio); the `/(2*r0)` gradient sensor (finite-difference, not ratio); the `2*amr_cpat_mar` two-sided margin adds.

---

## Task 1: Param + gate + mechanical de-hardcode (Cat 1, 2, 3, 7)

**Files:** `toolchain/mfc/params/definitions.py`, `src/simulation/m_amr.fpp`, `src/simulation/m_amr_registers.fpp`, `src/common/m_derived_types.fpp` (if a scalar decl is needed — likely not, ref_ratio is a scalar namelist var), `src/simulation/m_checker.fpp`, `src/common/m_global_parameters_common.fpp` (default assignment).

- [ ] **Add the param.** `definitions.py`: `_r('ref_ratio')` + `_nv()` NAMELIST_VARS registration; default 2 (in `s_assign_default_values_to_user_inputs`). Re-run cmake so the Fortran decl + namelist binding are generated.
- [ ] **Checker gate** (`m_checker.fpp`, AMR block): `@:PROHIBIT(ref_ratio /= 2 .and. ref_ratio /= 4, "ref_ratio must be 2 or 4")` and `@:PROHIBIT(ref_ratio /= 2 .and. (amr_max_level > 1 .or. amr_subcycle), "ref_ratio /= 2 is only supported at amr_max_level = 1 without subcycling (v1)")`.
- [ ] **Drive the per-block field:** `m_amr.fpp:5418` `amr_slots(islot)%ref_ratio = 2` → `= ref_ratio`.
- [ ] **Cat 1/2/3/7 swaps** at every line listed above: `2*` → `ref_ratio*`, `2**amr_block_level` → `ref_ratio**amr_block_level`, `/2`/`/2._wp` → `/ref_ratio`/`/real(ref_ratio,wp)`, `(m_glb+1)/2` → `(m_glb+1)/ref_ratio`, margin `(buff_size+1)/2+1` → `(buff_size+ref_ratio-1)/ref_ratio+1`. In `m_amr_registers.fpp` use the module `ref_ratio` (add `use` if needed) or the level-1 value. **Do NOT touch Cat 8 (Task 3), Cat 5/6 reflux (Task 4), or subcycle sites.**
- [ ] **Build (sandbox off) + gate:** existing AMR goldens byte-identical at ref_ratio=2 (default). Run a representative subset (`BCBA6E74`, a core single-level AMR golden, a multi-fluid AMR golden) `-- -b mpirun`. Commit `amr(ref_ratio): param + gate + mechanical de-hardcode of the 2:1 extent/index/margin factors (ref_ratio=2 byte-identical)`.

## Task 2: Prolong child-offset generalization (Cat 4 — 15 sites)

**Files:** `src/simulation/m_amr.fpp`.

Current (rr=2): `xi = (real(mod(f, rr), wp) - 0.5_wp)*0.5_wp` → children at ±0.25.
General: `xi = (real(mod(f, rr), wp) - real(rr - 1, wp)*0.5_wp) / real(rr, wp)`.
Verify byte-identity at rr=2: `(mod-0.5)/2` = `(0-0.5)/2=-0.25`, `(1-0.5)/2=0.25` — bit-identical to the old form (0.5, 0.25 exact in binary). For rr=4: `(0,1,2,3 - 1.5)/4` = `-0.375,-0.125,0.125,0.375` — symmetric (mean 0 ⇒ conservative), correct.

- [ ] Replace the offset at all 15 sites (1363, 1366, 1369, 1437, 1440, 1443, 1484, 1487, 1490, 2871, 2876, 2879, 2920, 2925, 2928) with the general form, using the `rr`/`ref_ratio` already in scope at each site. The QBMM `inject`/`pw_const` path zeroes the slopes so it is offset-independent — but fix it anyway for the non-inject variables that share the loop.
- [ ] **Build + gate:** existing goldens byte-identical (the algebra is exact at rr=2 — verify with a golden diff, not just PASS). Commit `amr(ref_ratio): generalize prolong child-offset stencil (rr=2 bit-identical, symmetric+conservative at rr=4)`.

## Task 3: Fine-grid coordinate generation (Cat 8 — redesign)

**Files:** `src/simulation/m_amr.fpp` — `s_build_level_coords` (~459–472), `s_amr_swap_to_fine` (~2609–2668).

- [ ] **`s_build_level_coords`:** the `c = lo + fi/2` map is Task 1 (`/ref_ratio`); the `mod(fi,2)` left/right-half branch (463) needs the general sub-boundary placement. A fine cell-boundary index `fi` sits at coarse cell `c = lo + fi/rr`, sub-index `k = mod(fi, rr)`, fractional position `k/rr` within `[xcb(c-1), xcb(c)]`: `x_cb(fi) = xcb(c-1) + (real(k,wp)/real(rr,wp))*(xcb(c) - xcb(c-1))`. Read the full routine and replace the even/odd branch with this. rr=2: k∈{0,1} → fractions 0, 0.5 → coarse edge, coarse midpoint — matches the old branch (verify the index convention).
- [ ] **`s_amr_swap_to_fine` ghost-coord maps (2609–2668):** same generalization for the 6 ghost-extension branches (hi/lo × x/y/z) — replace `mod(jg,2)` midpoint-vs-edge with the `k=mod(jg,rr)` fractional placement.
- [ ] **Build + gate:** existing goldens byte-identical (rr=2 reduces to bisection). This is the trickiest task — diff the `x_cb`/grid output of a golden to confirm bit-identity, not just physics PASS. Commit `amr(ref_ratio): generalize fine-grid coordinate subdivision to rr-way (rr=2 bit-identical)`.

## Task 4: Coarse-RHS reflux generalization (Cat 5/6)

**Files:** `src/simulation/m_amr_registers.fpp` — `s_amr_apply_reflux` (~519–639).

- [ ] Read `s_amr_apply_reflux` AND the shared `s_amr_reflux_apply_faces` (~721–818, already `rr`-general). **Preferred:** refactor `s_amr_apply_reflux` to call the shared kernel (completing the #39 dedup) so there is ONE reflux averaging path. **Fallback if the call shapes don't line up:** thread `rr` into `s_amr_apply_reflux`'s own loops — `nch` products (549–552, 582–584, 614) use `rr`; `f*0 = 2*c` (560–561, 592–593, 622–623) → `rr*c`; `do dd = 0, 1` (596, 625–626) → `do dd = 0, rr-1`.
- [ ] **Build + gate:** existing goldens byte-identical + per-fluid mass conservation ~machine-zero (reflux is conservation-critical). Commit `amr(ref_ratio): generalize coarse-RHS reflux child-averaging to rr (dedup into shared kernel)`.

## Task 5: ref_ratio=4 golden + validation

**Files:** `toolchain/mfc/test/cases.py`, new `tests/<UUID>/`.

- [ ] Add a **single-level, lock-step, ref_ratio=4** golden: a 1D (or 2D) inviscid case with `amr=T, amr_max_level=1, amr_subcycle=F, ref_ratio=4`, a static block over a smooth feature, bit-uniform grid. Mirror an existing single-level AMR golden's params + set `ref_ratio=4`.
- [ ] **Generate + validate CPU:** (a) conservation machine-zero (inviscid, no IB) — the primary gate; a prolong-offset or reflux bug shows here. (b) A **resolution check**: the 4:1 block's fine solution should match a two-level 2:1 tower (amr_max_level=2, ref_ratio=2) to discretization order on the same feature — confirms the 4:1 refinement is physically placed, not just conservative. (c) Optionally an amr_max_level=1 ref_ratio=2 vs a coarser ref_ratio=4-of-half-the-blocks equivalence.
- [ ] **Validate GPU (V100):** `clean --gpu acc` build, run the new golden + a ref_ratio=2 regression at np=1 `-- -b mpirun`. Commit `amr(ref_ratio): single-level ref_ratio=4 golden + conservation/resolution validation`.

---

## Cross-cutting notes
- After EACH task, a golden DIFF (not just PASS) is the byte-identity proof for the rr=2 default — this high-risk work relies on that invariant holding at every step.
- The conservation-visible prolong offset (Task 2) + per-fluid mass check (Task 4) are the two places a silent ratio bug becomes loud — lean on them.
- Deferred follow-ups (own increments): subcycled ref_ratio=4 (subcycle site generalization + 4-substep ghost-lerp/reflux conservation); multi-level ref_ratio=4 (16:1 towers — the `/4` inset heuristics, `rr**level` cumulative factors, fine-fine seam at mixed depth). Both gated closed by Task 1's checker.
