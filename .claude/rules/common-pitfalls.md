# Common Pitfalls

Traps that live in no other doc, collected because their failure modes are silent —
wrong answers and wrong indices, not error messages. General development pitfalls are
covered in `docs/documentation/contributing.md`.

## Indexing and Ghost Cells

- Grid dimensions `m`, `n`, `p` (cells in x, y, z); 1D: n=p=0, 2D: p=0. Interior `0:m`;
  ghost region `-buff_size:m+buff_size`; bounds structs `idwint(1:3)` (interior) and
  `idwbuff(1:3)` (with ghosts); cell boundaries `x_cb(-1-buff_size:m+buff_size)`.
- `buff_size` is **not** a single formula: it's set per reconstruction scheme in
  `s_configure_coordinate_bounds` (`m_helper_basic.fpp`) and floored higher for Lagrange
  bubbles and IB. Read that routine for the current value rather than assuming one.
- Riemann solvers: left state at `j`, right state at `j+1`.
- All equation indices live in the `eqn_idx` struct (`eqn_idx_info` in
  `m_derived_types.fpp`, populated in `m_global_parameters.fpp`): `%cont`, `%mom`, `%E`,
  `%adv`, plus optional ranges (`%bub`, `%stress`, `%species`, `%B`, ...). The old
  `contxb`/`momxb` shorthands are gone. Index positions depend on `model_eqns` and
  enabled features — changing either moves ALL indices; never hard-code one.

## GPU

- WARNING: do NOT wrap `GPU_LOOP` in `GPU_PARALLEL` for spatial loops — `GPU_LOOP` emits
  empty directives on Cray and AMD, causing silent serial execution. Spatial loops always
  use `GPU_PARALLEL_LOOP`/`END_GPU_PARALLEL_LOOP`. Macro API:
  `docs/documentation/gpuParallelization.md`; signatures:
  `src/common/include/parallel_macros.fpp`. Never call the `ACC_*`/`OMP_*`
  implementation layers directly.
- Only `src/simulation/` is GPU-accelerated. Backends: OpenACC (nvfortran primary, Cray)
  and OpenMP offload (Cray primary, AMD flang, nvfortran). The CPU-only build must always
  work — every `#ifdef` needs a path for all configurations (CPU, ACC, OMP, with/without
  MPI). Gates: `MFC_GPU`, `MFC_OpenACC`, `MFC_OpenMP`, `MFC_MPI`, `MFC_DEBUG`,
  `MFC_SINGLE_PRECISION`/`MFC_MIXED_PRECISION`, `MFC_PRE_PROCESS`/`MFC_SIMULATION`/
  `MFC_POST_PROCESS`, and compiler macros (`_CRAYFTN`, `__PGI`, ...).
- `@:ACC_SETUP_VFs(...)`/`@:ACC_SETUP_SFs(...)` GPU pointer setup compiles only under
  Cray. Around MPI: `GPU_UPDATE(host=...)` before send, `GPU_UPDATE(device=...)` after
  receive.

## Parameters

- Adding one: `_r()` definition + `_nv()` `NAMELIST_VARS` registration in
  `toolchain/mfc/params/definitions.py`; `case_validator.py` only if physics-constrained
  (with a `PHYSICS_DOCS` entry). Fortran declarations and namelist bindings are
  auto-generated at CMake configure time — re-run cmake (or `./mfc.sh build`) after editing.
- Still manual: array variables and derived-type members (declare in
  `src/*/m_global_parameters.fpp` / the relevant type), and case-optimization parameters
  (`CASE_OPT_PARAMS` + the `#:else` block in `src/simulation/m_global_parameters.fpp`).
  Gotcha: under `--case-optimization` those are baked into the binary and dropped from
  the namelist, so changing one needs a *rebuild*, not a case edit.
- Runtime checks (`@:PROHIBIT`) go where they run: shared →
  `src/common/m_checker_common.fpp`; simulation-only → `src/simulation/m_checker.fpp`;
  pre/post-only → `src/{pre,post}_process/m_checker.fpp` (their `s_check_inputs` are
  currently empty — that IS the right place, not m_checker_common).
- Analytic ICs are compiled into the binary (syntax error = build failure, not runtime).
  Each IC variable maps to an `eqn_idx%…` expression in `QPVF_IDX_VARS`
  (`toolchain/mfc/case.py`); a new patch-settable conserved variable means updating that
  map AND the Fortran `eqn_idx` builder to agree — a mismatch is a silent wrong index.
  Variables available in expressions: `docs/documentation/case.md`.

## Tests

- Tests are generated programmatically in `toolchain/mfc/test/cases.py` (parameter
  modifications on `BASE_CFG` via the `CaseGeneratorStack` push/pop pattern); test UUID =
  CRC32 of the trace string; `./mfc.sh test -l` lists all.
- Golden files are tolerance-compared. Regenerate only the affected tests
  (`./mfc.sh test --generate --only <tests>`) — an unexplained golden-file diff is a bug
  report, not noise to be regenerated away.
