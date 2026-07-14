# AMR on AMD flang (OpenMP target offload) — known failure and root cause

**Status:** the block-structured AMR GPU path does **not** yet run correctly on the AMD flang
OpenMP-target-offload backend (`gpu-omp` on Frontier). It works on NVHPC (OpenACC), NVHPC
(OpenMP offload), and **Cray CCE (OpenMP offload)**. This note records the root-cause analysis so the
fix can be made deliberately; until then the AMR tests are expected to fail only on the AMD `gpu-omp` lane.

## Symptom (CI: Frontier AMD `gpu-omp`)
AMR tests at **np ≥ 2** fail two ways, while the same tests pass on every other backend (including
Cray CCE `gpu-omp`):
- **Self-abort** — `244B1E42` (2D single-level subcycle np=2), `ADA042A2` (2D multi-level static np=2):
  the simulation calls `MPI_Abort` (`s_mpi_abort`) because a runtime guard trips on an **unphysical
  state** (NaN / negative density or pressure).
- **Tolerance mismatch** — the 2D stretched-grid **dynamic-regrid** cases (`660FFBFE`, `B7704247`):
  complete but the field diverges beyond golden tolerance.

np = 1 AMR passes on AMD; the failure is in the **np ≥ 2-only** machinery (seam-halo /
gather-scatter kernels that dereference nested derived-type slot members
`amr_slots(slot)%q_cons(i)%sf`). All P2P is host-staged (`GPU_UPDATE(host)` → MPI → `GPU_UPDATE(device)`);
there is no GPU-aware-MPI / `use_device_ptr` path, so this is a **device data-mapping** problem, not a GTL one.

## Root cause (primary)
`src/common/include/omp_macros.fpp` gives each backend a different default data-mapping clause on
`GPU_PARALLEL_LOOP` (`OMP_DEFAULT_STR` / the start-directive builder):
- **NVHPC:** `defaultmap(firstprivate:scalar)` + `defaultmap(tofrom:allocatable/pointer)`
- **Cray CCE:** `defaultmap(firstprivate:scalar)` + `defaultmap(present:allocatable) defaultmap(present:pointer)`
- **AMD flang:** *bare directive — NO `defaultmap` clause at all.*

With no `defaultmap`, flang applies the OpenMP **implicit `tofrom`** rule to every unlisted
variable in the kernel, including the AMR kernels' **assumed-shape / derived-type slot dummies**.
Implicitly mapping such a descriptor re-copies a **host** descriptor (whose embedded data pointer is a
host address) over the device mapping that `@:ALLOCATE` established, so the kernel reads garbage →
NaN → abort. CCE avoids this by asserting `present`; NVHPC by OpenACC-style auto-attach.

Two contributing factors compound it:
- **`ACC_SETUP_SFs/VFs`** (the manual device-descriptor attach for nested slot members) is
  `#ifdef _CRAYFTN` **only** — flang gets no descriptor attach, only the leaf `%sf` data is mapped.
- **`ACC_TEARDOWN_SFs`** (delete the mapping before a slot is freed) is likewise Cray-only, so on
  **dynamic regrid** flang reuses a stale present-table entry for a freed-then-reallocated `%sf` —
  which is exactly why the abort cases are static/subcycle while the divergence cases are dynamic-regrid.

## Fix direction
1. In `src/common/include/omp_macros.fpp`, give **AMD** the same defaultmap string as CCE:
   `defaultmap(firstprivate:scalar) defaultmap(present:allocatable) defaultmap(present:pointer)`
   (or add explicit `present=[...]` to the offending AMR kernels).
2. Add an AMD/flang branch to `ACC_SETUP_SFs/VFs` and `ACC_TEARDOWN_SFs` (`src/common/include/macros.fpp`)
   emitting the equivalent `!$omp target enter/exit data map(...)` for the descriptor + `%sf`.

## Confirming experiment (Frontier, no code edit)
Run one failing np=2 AMR case on the AMD `gpu-omp` build with mapping diagnostics:
```
export LIBOMPTARGET_INFO=-1
export OMP_TARGET_OFFLOAD=MANDATORY
./mfc.sh test --only 244B1E42 -- -b mpirun -n 2
```
If the abort is immediately preceded by an **implicit map** / "present table entry not found" on an
`amr_slots(...)%q_cons(...)%sf` descriptor, the empty-defaultmap is confirmed. A/B: temporarily set the
AMD branch of the start directive to the CCE string and re-run — if it passes (or the abort moves), done.

## Full ranked top-10 (summary)
1. **flang `GPU_PARALLEL_LOOP` emits no `defaultmap`** → implicit `tofrom` re-maps derived-type/assumed-shape slot dummies. (`omp_macros.fpp:183-184`, `:30-31`) — **primary.**
2. `ACC_SETUP_SFs` descriptor attach is `#ifdef _CRAYFTN` only → flang gets no attach for nested slot members. (`macros.fpp:111-124`, used in `s_amr_alloc_slot`)
3. `GPU_UPDATE(host/device)` on a partial array-section of a doubly-nested derived-type member in the seam halo → flang miscomputes section offset/stride. (`m_amr.fpp:3165/3196`) — best fit for the tolerance-mismatch pair.
4. Missing `defaultmap(firstprivate:scalar)` on AMD → module-scope loop-bound scalars implicitly `tofrom`-mapped/stale. (`omp_macros.fpp:184`)
5. `GPU_UPDATE(host=)` on an assumed-shape dummy (`pb_coarse`/`mv_coarse`, the documented #29 ⚠) — QBMM/bubble subset. (`m_amr.fpp:738`)
6. `declare target` for deferred-alloc module arrays (`amr_cg_pb/mv`, `amr_rhs_*`) — flang links descriptor weakly. (`m_amr.fpp:110/116/161`)
7. `GPU_LOOP` empty on both CCE and AMD (shared trap; only compounds AMD when combined with #1/#4). Low.
8. `GPU_UPDATE` has no barrier/`nowait` — possible flang `target update`↔host-MPI ordering race (would be run-to-run flaky). (`omp_macros.fpp:305-312`)
9. Derived-type-expression lower bounds on mapped assumed-shape dummies (`…%idwbuff(1)%beg:`) → flang descriptor-bound mishandling → off-by-buff indexing. (`m_amr.fpp:3031-3034`)
10. `map(always,alloc)` vs `present,alloc` on freed-then-realloc'd slots across regrid, with the Cray-only teardown missing → stale mapping. (`omp_macros.fpp:60/78`) — ties to the dynamic-regrid divergence.

**Single most likely:** #1 (empty defaultmap), with #2 the same disease from the allocation side; fix together.
