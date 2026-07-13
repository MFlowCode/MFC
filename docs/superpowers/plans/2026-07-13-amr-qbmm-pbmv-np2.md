# Distributed pb/mv coupling for non-polytropic QBMM at np≥2 — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the non-polytropic QBMM pb/mv quadrature side-state couple correctly across MPI ranks under AMR, so `qbmm .and. .not. polytropic` runs on more than one rank (currently fail-closed in the checker).

**Architecture:** pb/mv (bubble pressure, vapor mass; `nnode × nb` per cell) evolve cell-locally — no face flux, therefore **no Berger-Colella reflux**. Their only coarse↔fine coupling is prolong (coarse→fine ghosts+interior) and restrict (fine→coarse fold-back). Today both read/write the block **owner's LOCAL** coarse pb/mv, which is exact only at np=1. This plan mirrors the proven q_cons P2P distribution for pb/mv as a **separate exchange** (q_cons buffers untouched → q_cons stays byte-identical by construction): a gathered coarse pb/mv patch (`amr_cg_pb`/`amr_cg_mv`, the analogue of `amr_cg`) that the prolong/ghost-fill read, and a restrict scatter that folds fine averages back to the coarse-cell owners.

**Tech Stack:** Fortran 2008 + Fypp, OpenACC/OpenMP-offload GPU, MPI P2P (`MPI_ISEND`/`MPI_IRECV`), MFC AMR module `src/simulation/m_amr.fpp`.

## Global Constraints

- **Scope = single-level np≥2 only.** Mirror `s_amr_gather_coarse_patch` (the L0↔L1 gather), NOT the parent gather. Multi-level QBMM np≥2 (`amr_max_level>1`) stays fail-closed.
- **np=1 gate goldens (Tasks 1-3):** `BCBA6E74` (AMR → 1D → bubbles QBMM → nonpolytropic) is the primary gate; `B0A5D230` (same + regrid subcycle) covers the subcycle ghost-fill path. Run both at np=1 after each task.
- **np=1 byte-identical.** Every task's np=1 path must reduce to the current local read/write (single owner holds every covered/coarse cell). Existing goldens unchanged. This is the primary per-task gate.
- **No reflux for pb/mv.** pb/mv have no flux; do NOT add any reflux/flux-correction. Prolong sets fine, restrict folds back — that is the entire coupling.
- **GPU macros only** (`GPU_*` Fypp), precision via `wp`/`stp`, abort via `s_mpi_abort`/`@:PROHIBIT`, `@:ALLOCATE`↔`@:DEALLOCATE` paired. MPI type `mpi_p` ↔ `wp` on the wire (messages carry `wp`, cast to `stp` on unpack — mirror q_cons at `m_amr.fpp:615`, `:1575`).
- **Validation oracle for pb/mv** (a non-conservative side-state): **np=2 trajectory == np=1 trajectory on a bit-uniform grid** — the same np-cross check q_cons uses. Non-bit-uniform grids diverge at ulp (WENO-table finding), so the golden grid must be bit-uniform.
- Reference frames (memorize both):
  - **Local coarse frame:** `ci = amr_isect_lo(d) + f/rr - start_idx(d)` (what pb/mv reads TODAY — owner's own coarse array).
  - **Patch-local frame:** `ci = amr_isect_lo(d) + f/rr - amr_cpat_off(d)` (what q_cons prolong reads from `amr_cg` — the gathered patch; see `s_prolong_one_var` `m_amr.fpp:1149-1159`).
  The whole change is: read/write the gathered patch in the patch-local frame instead of the local coarse array.
- Workflow: `./mfc.sh format` → build (`-t simulation --gpu acc` for GPU tasks) → precheck-on-commit (never `--no-verify`) → targeted golden. GPU validate on V100 via `--only <UUID> -- -b mpirun` (NOT full-suite sbatch — it SIGILL-flakes on arch-mismatched nodes; see `phoenix-srun-mpirun` memory, prepend hpcx bin to PATH after `source ./mfc.sh load -c p -m g`).

---

## Task 1: Gathered coarse pb/mv patch + distributed prolong (init/regrid path)

**Files:**
- Modify: `src/simulation/m_amr.fpp` — module state (near `amr_cg` decl ~line 153), allocation (~368-375), finalize dealloc, new gather routine (after `s_amr_gather_coarse_patch` ~651), `s_amr_prolong_pbmv` (~1887), its callers (1340 init, 4316/4355/4799 regrid).
- Test: existing QBMM+AMR golden (find via `./mfc.sh test -l | grep -i qbmm` — the non-polytropic AMR case) run at np=1.

**Interfaces:**
- Produces: `s_amr_gather_coarse_patch_pbmv(pb_coarse, mv_coarse, pull_host)` — gathers the current block's coarse pb/mv patch into module arrays `amr_cg_pb`/`amr_cg_mv` (block-local/patch frame, cell 0 = global `amr_cpat_off(d)`). `pull_host=.false.` ⇒ host copy current (init/regrid host prolong); `.true.` ⇒ device current (runtime). Signature/semantics mirror `s_amr_gather_coarse_patch(q_coarse, pull_host)`.
- Consumes: `f_amr_rank_coarse_range`, `amr_cpat_off`, `amr_cpat_hi`, `amr_isect_lo`, `amr_slots(amr_cur)`, `nnode`, `nb` (all existing).

- [ ] **Step 1: Add module state.** After the `amr_cg` declaration block (~line 153-156), add:
```fortran
    !! Gathered coarse pb/mv patch for non-polytropic QBMM (analogue of amr_cg): the block's coarse-side pb/mv side-state,
    !! P2P-gathered from the coarse-cell owners into the block owner in the amr_cg patch-local frame (cell 0 == amr_cpat_off).
    !! Read by the pb/mv prolong + ghost-fill so np>=2 couples to the correct coarse rank. Allocated only for non-polytropic QBMM.
    real(stp), allocatable, dimension(:,:,:,:,:) :: amr_cg_pb, amr_cg_mv
    $:GPU_DECLARE(create='[amr_cg_pb, amr_cg_mv]')
```

- [ ] **Step 2: Allocate/map** next to `amr_cg` (after line 375), gated on non-polytropic QBMM, sized to the same patch footprint as `amr_cg` plus the pb/mv `(nnode, nb)` trailing dims:
```fortran
        if (qbmm .and. .not. polytropic) then
            @:ALLOCATE(amr_cg_pb(0:amr_cpat_hi(1), 0:amr_cpat_hi(2), 0:amr_cpat_hi(3), 1:nnode, 1:nb))
            @:ALLOCATE(amr_cg_mv(0:amr_cpat_hi(1), 0:amr_cpat_hi(2), 0:amr_cpat_hi(3), 1:nnode, 1:nb))
            amr_cg_pb = 0._stp; amr_cg_mv = 0._stp
        end if
```
NOTE for implementer: `amr_cg_pb`/`amr_cg_mv` are plain 5D arrays, NOT `scalar_field`s. Follow the **`amr_rhs_pb_f` idiom EXACTLY** (`m_amr.fpp:115-116` module `GPU_DECLARE(create=...)` + `m_amr.fpp:349` `@:ALLOCATE`): the module `GPU_DECLARE` (Step 1) + `@:ALLOCATE` handle device mapping. Do NOT use `@:ACC_SETUP_SFs` (that is only for `scalar_field` `%sf` pointers, e.g. `amr_cg(i)`). Do NOT invent a scalar_field wrapper. Note the payload is `stp` (matches `pb_f%sf`/`pb_ts`), whereas `amr_rhs_pb_f` is `wp` — use `real(stp)` for `amr_cg_pb`/`amr_cg_mv`.

- [ ] **Step 3: Finalize dealloc.** Where `amr_rhs_pb_f`/`amr_cg` are deallocated in `s_finalize_amr_module`, add the paired `@:DEALLOCATE(amr_cg_pb)` / `@:DEALLOCATE(amr_cg_mv)` under the same `qbmm .and. .not. polytropic` guard (grep the finalize routine for `amr_rhs_pb_f` to co-locate).

- [ ] **Step 4: Write the gather routine.** Add `s_amr_gather_coarse_patch_pbmv(pb_coarse, mv_coarse, pull_host)` after `s_amr_gather_coarse_patch` (~651). Structure it as a **verbatim structural mirror** of `s_amr_gather_coarse_patch` (read `m_amr.fpp:505-651` first), with these payload differences:
  - Dummies: `real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(in) :: pb_coarse, mv_coarse` (the coarse-level pb/mv, e.g. `pb_ts(1)%sf`).
  - Assert single-level: this routine is single-level only. Early-guard `if (amr_block_level(amr_cur) >= 2) return` — multi-level QBMM np≥2 is gated in the checker; the routine must never be reached at level≥2 (Task 4 gate enforces it).
  - Per-cell payload = `2*nnode*nb` (pb then mv). Buffer sizes: replace `sys_size` with `2*nnode*nb` in `boxsz`/`maxsz`; pack loop order `do ib_=1,nb; do q=1,nnode` for pb, then the same for mv (or interleave — keep sender/receiver order identical, as q_cons does over `i=1,sys_size`).
  - np=1 device-kernel local-copy path: mirror lines 543-566 (`num_procs == 1` branch) — a `GPU_PARALLEL_LOOP` copying `pb_coarse`→`amr_cg_pb`, `mv_coarse`→`amr_cg_mv` over the in-domain patch, same index map `amr_cg_pb(g1-coff1, g2-coff2, g3-coff3, q, ib_) = pb_coarse(g1-o1, g2-o2, g3-o3, q, ib_)`. Hoist `coff*`/`o*` like the q_cons kernel. Device-current; sync host only when `.not. pull_host` (mirror the `to_host` gate at 625/726).
  - np≥2 P2P path: mirror 577-651 — own-slice local unpack + `MPI_IRECV`/`MPI_ISEND` over `f_amr_rank_coarse_range`, `wp` on the wire, cast to `stp` into `amr_cg_pb`/`amr_cg_mv`. `GPU_UPDATE(device=...)` after receive; `GPU_UPDATE(host=...)` before the host-consumer path when `.not. pull_host`.

- [ ] **Step 5: Rewrite `s_amr_prolong_pbmv`** (1887) to read the gathered patch. Two edits:
  - Add a gather call at the top: `call s_amr_gather_coarse_patch_pbmv(pb_ts(1)%sf, mv_ts(1)%sf, .false.)` (host-current: this is a host loop).
  - Change the coarse read from local `pb_ts(1)%sf(ci,cj,ck,q,ib_)` to `amr_cg_pb(ci,cj,ck,q,ib_)` and the index origin from `ox = start_idx(1)` (etc.) to `ox = amr_cpat_off(1); oy = amr_cpat_off(2); oz = amr_cpat_off(3)` (patch-local frame, matching `s_prolong_one_var:1151`). Keep `ci = lo1 + fi/rr - ox` shape.

- [ ] **Step 6: Build** (`./mfc.sh build -t simulation --gpu acc -j 8`). Expected: compiles clean.

- [ ] **Step 7: np=1 golden byte-identical.** Run the existing non-polytropic-QBMM+AMR golden at np=1 (`./mfc.sh test --only <UUID> -- -b mpirun`). Expected: PASS unchanged (single owner ⇒ gather's np=1 branch copies local coarse→patch, prolong reads the same values it read before, just via the patch frame).

- [ ] **Step 8: Commit** `amr(qbmm): gather coarse pb/mv patch + distributed prolong (init/regrid, np=1 byte-identical)`.

---

## Task 2: Distributed pb/mv ghost-fill (runtime advance path)

**Files:**
- Modify: `src/simulation/m_amr.fpp` — `s_amr_fill_fine_ghosts_pbmv` (~1925) index frame, and its callers (2864, 2959, 2960, 3053) to gather-then-read.
- Test: same np=1 QBMM+AMR golden.

**Interfaces:**
- Consumes: `s_amr_gather_coarse_patch_pbmv` (Task 1), `amr_cg_pb`/`amr_cg_mv`, `amr_cpat_off`.

- [ ] **Step 1: Change `s_amr_fill_fine_ghosts_pbmv` index origin** (1925). Its dummies `pb_c`/`mv_c` will now receive `amr_cg_pb`/`amr_cg_mv`. Change the dummy declaration bounds to the patch array bounds `real(stp), dimension(0:,0:,0:,1:,1:), intent(in) :: pb_c, mv_c` (the gathered patch is 0-based; DO NOT keep the `idwbuff(1)%beg:` bounds — the patch has no ghost origin). Change `ox = start_idx(1)` (etc.) to `ox = amr_cpat_off(1); oy = amr_cpat_off(2); oz = amr_cpat_off(3)`. Keep the floor-based `ci = lo1 + floor(fi/rr) - ox` map.

- [ ] **Step 2: Gather before each runtime ghost-fill.** At each caller, insert a device-current pb/mv gather immediately before the fill, reading the SAME coarse snapshot the paired q_cons gather uses, and pass the patch:
  - Line 2864: after the q_cons `s_amr_fill_fine_ghosts(amr_cg, ...)` at 2858, replace `call s_amr_fill_fine_ghosts_pbmv(pb_in, mv_in, ...)` with `call s_amr_gather_coarse_patch_pbmv(pb_in, mv_in, .true.)` then `call s_amr_fill_fine_ghosts_pbmv(amr_cg_pb, amr_cg_mv, ...)`.
  - Lines 2959/2960: pair with the q_cons gathers at 2951 (`q_old`) / 2953 (`q_new`) — gather `(pb_old,mv_old,.true.)` before the 2959 fill, `(pb_in,mv_in,.true.)` before the 2960 fill.
  - Line 3053 (`s_amr_lerp_fine_ghosts_pbmv`): **CONFIRMED reads pre-filled ghost shells `pb_ghost_a/b`, NOT a coarse array** (verified `m_amr.fpp:3051-3055`). NO new gather here — leave it unchanged.

- [ ] **Step 3: Build.** Expected: clean.

- [ ] **Step 4: np=1 golden byte-identical.** Same UUID, `-- -b mpirun`. Expected: PASS unchanged.

- [ ] **Step 5: Commit** `amr(qbmm): distributed pb/mv ghost-fill on the runtime advance path (np=1 byte-identical)`.

---

## Task 3: Distributed pb/mv restrict scatter (fine→coarse fold-back)

**Files:**
- Modify: `src/simulation/m_amr.fpp` — replace the two local `s_restrict_pbmv` calls (1509 np=1 branch, 1590 np≥2 branch) with a distributed scatter; add a `s_amr_scatter_pbmv` routine (or extend `s_restrict_pbmv` to scatter).
- Test: same np=1 QBMM+AMR golden.

**Interfaces:**
- Consumes: `f_amr_rank_interior`, `s_amr_box_isect`, `amr_region_lo_all`/`amr_region_hi_all`, `amr_block_owner`, the restrict child-average arithmetic in `s_restrict_pbmv`.

- [ ] **Step 1: Read the q_cons scatter** (`s_restrict_fine_to_coarse:1479-1587`) — the owner restricts covered cells locally + sends each other coarse-owner its covered slice; the coarse-owner receives and overwrites local coarse. This is the exact pattern to mirror for pb/mv.

- [ ] **Step 2: Add `s_amr_scatter_pbmv(pb_coarse, mv_coarse, pb_fin, mv_fin)`** mirroring the q_cons scatter structure but with the pb/mv child-average (from `s_restrict_pbmv:2059-2070`) as the per-cell restrict value and `2*nnode*nb` payload:
  - Owner branch: over covered cells `[region_lo:region_hi]`, for cells this rank owns (`f_amr_rank_interior(owner) ∩ region`), overwrite `pb_coarse`/`mv_coarse` locally with the child average (device kernel, mirror `s_amr_restrict_overwrite_device`; or reuse the existing `s_restrict_pbmv` for the owner-local covered box). For every OTHER coarse-owner, pack its covered slice (child averages, `wp`) and `MPI_ISEND`.
  - Coarse-owner branch: if it holds covered cells, `MPI_RECV` its slice, cast `wp`→`stp` into local `pb_coarse`/`mv_coarse`, `GPU_UPDATE(device=...)` only the covered slice (mirror 1583-1585 — never whole-array, which clobbers device-advanced non-covered coarse).
  - np=1: owner owns every covered cell, sends nothing, overwrites locally ⇒ **identical child-sum to the current `s_restrict_pbmv`** (bit-identical).

- [ ] **Step 3: Wire it in.** Replace line 1509 (`call s_restrict_pbmv(...)`, np=1 device branch) and line 1590 (np≥2 branch) with `call s_amr_scatter_pbmv(pb_ts(1)%sf, mv_ts(1)%sf, amr_slots(amr_cur)%pb_f%sf, amr_slots(amr_cur)%mv_f%sf)`. Remove the `amr_rank_owns_block` guard from the call site — the scatter must run on ALL ranks (owner sends, coarse-owners receive); do the owner/non-owner split INSIDE the routine (like `s_restrict_fine_to_coarse`, which is called unconditionally and branches on `proc_rank == owner`). Keep the level≥2 early-return upstream (1473-1477) unchanged — multi-level fold-back stays gated.

- [ ] **Step 4: Build.** Expected: clean.

- [ ] **Step 5: np=1 golden byte-identical.** Same UUID, `-- -b mpirun`. Expected: PASS unchanged.

- [ ] **Step 6: Commit** `amr(qbmm): distributed pb/mv restrict scatter (np=1 byte-identical)`.

---

## Task 4: Lift the single-level gate + np=2 golden

**Files:**
- Modify: `src/simulation/m_checker.fpp` (~106) — narrow the gate to multi-level only.
- Modify: `toolchain/mfc/test/cases.py` — add a np=2 non-polytropic-QBMM + single-level-AMR golden on a bit-uniform grid.
- Test: new golden UUID (CPU + GPU).

- [ ] **Step 1: Narrow the checker gate** (`m_checker.fpp:106`). Change:
```fortran
@:PROHIBIT(qbmm .and. (.not. polytropic) .and. num_procs > 1, "...")
```
to gate only the still-unsupported MULTI-level case:
```fortran
@:PROHIBIT(qbmm .and. (.not. polytropic) .and. amr_max_level > 1 .and. num_procs > 1, &
           & "amr with non-polytropic QBMM on more than one MPI rank is only supported at amr_max_level = 1: the multi-level (parent-side) pb/mv coarse/fine coupling is not yet distributed. Run multi-level non-polytropic QBMM on a single rank.")
```
Update the explanatory comment above it (101-105) to say single-level np≥2 is now distributed; multi-level remains TODO.

- [ ] **Step 2: Add the np=2 golden** in `cases.py`, mirroring the `BCBA6E74` case (trace `"AMR -> 1D -> bubbles nonpolytropic"` region at `cases.py:3440`; find the QBMM-nonpolytropic sibling near it). Requirements: `amr_max_level=1`, run at **np=2** with a domain (raise `m`) that decomposes so the AMR block **straddles the rank boundary** — a coarse cell under the block owned by a DIFFERENT rank than the block owner. This is what exercises the gather/scatter; if the block sits entirely on one rank the P2P path is dead and the test proves nothing. **Bit-uniform grid** (uniform `dx`, no stretching) so the np-cross oracle holds. Distinct trace noting np=2.

- [ ] **Step 3: Generate + validate on CPU.** Build CPU (`./mfc.sh build -t simulation -j 8`), `./mfc.sh test --generate --only <new UUID> -- -b mpirun`. Then **prove the oracle by hand**: run the same case at np=1 and np=2, confirm the pb/mv moments (and q_cons) match to machine-zero (bit-uniform ⇒ exact). If np=2 ≠ np=1, the scatter/gather is wrong — do NOT accept the golden; return to Task 1-3.

- [ ] **Step 4: Validate on GPU (V100).** `source ./mfc.sh load -c p -m g`, prepend hpcx bin to PATH, `./mfc.sh build -t simulation --gpu acc -j 8`, `./mfc.sh test --only <new UUID> -- -b mpirun` at np=2. Expected: PASS (GPU np=2 == CPU golden). Also re-run the pre-existing QBMM+AMR np=1 golden on GPU to confirm no regression.

- [ ] **Step 5: Commit** `amr(qbmm): admit single-level non-polytropic QBMM at np>=2 + np=2 golden`.

---

## Cross-cutting notes
- The three coupling edits (gather+prolong, ghost-fill, scatter) are each independently np=1-byte-identical, so a regression in any one is caught by the existing QBMM+AMR golden before Task 4 ever runs at np≥2.
- If the existing QBMM+AMR golden's block does not straddle a rank boundary at np=2, Task 4's new case MUST arrange straddle — otherwise the np≥2 paths are dead-code-covered and the "distributed" claim is unproven.
- Do NOT touch `s_amr_restrict_to_parent` / `s_amr_gather_from_parent` (multi-level) — out of scope, kept gated.
