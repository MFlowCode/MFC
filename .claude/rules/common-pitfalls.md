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
  `m_derived_types.fpp`, populated by `s_initialize_eqn_idx` in `m_global_parameters_common.fpp`): `%cont`, `%mom`, `%E`,
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
- An array whose bound is a device global (`dimension(num_fluids)`, `dimension(num_species)`) may be
  passed to a device routine **from a parallel-loop body, but not from inside another
  `GPU_ROUTINE(parallelism='[seq]')`**. CCE OpenACC rejects the second form with
  `ftn-7066 ... Global in accelerator routine without declare -- num_fluids`, and reports it at
  whatever line it gave up on: remove one trigger and the message *walks forward* to the next call,
  so the reported line is not the cause. Only the plain lanes fail - under `--case-optimization`
  those bounds are `parameter`s, so a green Case Opt lane beside a failing plain one is the
  signature. Every accepted call site in the tree already obeys this (`m_cbc`, `m_ibm`,
  `m_bubbles_EL`, `s_compute_cell_state`): form such a call in the loop body and pass scalars
  deeper. Neither `cray_inline` nor a `num_fluids_max` bound nor dropping optional dummies helps -
  all three were measured.
- nvfortran 23.11/24.1 segfault (`fort2 TERMINATED by signal 11`) on a caller that passes a
  `parameter` array from `m_thermochem` (e.g. `molecular_weights`) into a declare-target routine.
  Read such arrays directly in the kernel, or pass a plain local computed from them.
- The same "call it from the loop body" rule covers `m_thermochem`: calling `get_species_*` from
  inside a `GPU_ROUTINE` rather than from the kernel gave CCE OpenMP a runtime
  `Memory access fault by GPU node-N ... Reason: Unknown` on the first step (exit 134), while every
  other backend ran. Evaluate them at the call site and pass the arrays in. Note this one only shows
  at runtime, and only on a case that reaches the path - the build is clean.

## Parameters

- Adding one: `_r()` definition + `_nv()` `NAMELIST_VARS` registration in
  `toolchain/mfc/params/definitions.py`; `case_validator.py` only if physics-constrained
  (with a `PHYSICS_DOCS` entry). Fortran declarations and namelist bindings are
  auto-generated at build time (ninja-tracked custom command) — re-run cmake (or `./mfc.sh build`) after editing.
- Still manual: derived-type `TYPE` member definitions in `src/common/m_derived_types.fpp`;
  default-value assignments in `s_assign_default_values_to_user_inputs`; the
  `CASE_OPT_EXTRA_LINES` literal in `toolchain/mfc/params/generators/fortran_gen.py` (covers `num_dims`,
  `num_vels`, `weno_polyn`, `muscl_polyn`, `weno_num_stencils`, `wenojs`);
  multi-variable declaration lines (`bc_x/y/z`, `x/y/z_domain`, `x/y/z_output`, post's
  `G`); and the MPI broadcast residue in `m_mpi_proxy` (computed variables that are not
  namelist-bound: `m_glb`/`n_glb`/`p_glb`, `cfl_dt`, `bc_io`, and complex struct-member
  array loops — these cannot be auto-generated and stay hand-listed). Everything else — scalar declarations, plain arrays (`FORTRAN_ARRAY_DIMS`
  table in `definitions.py`), derived-type namelist declarations including `GPU_DECLARE`
  lines and Doxygen descs (`TYPED_DECLS` table in `definitions.py`), the simulation
  case-optimization declaration block, and the per-target MPI broadcast lists for all
  namelist-registry scalars (`generated_bcast.fpp`) — is regenerated at build time by a
  ninja-tracked custom command (editing `params/*.py` triggers regeneration automatically).
  Gotcha: ADDING a new file under `toolchain/mfc/params/` needs one reconfigure
  (the custom command's DEPENDS list is globbed at configure time). Under `--case-optimization` the baked-in constants are dropped from the
  namelist, so changing one needs a *rebuild*, not a case edit.
- Derived-type params (`chem_params`, `lag_params`, `rburn`) are NOT auto-broadcast:
  `generated_bcast.fpp` covers namelist *scalars* only. Each type needs a hand-written
  `_emit_<name>` in `toolchain/mfc/params/generators/fortran_gen.py` plus its call site in
  the `target == "sim"` block, and — if it is read on device — an explicit
  `$:GPU_UPDATE(device='[name]')` in BOTH `m_global_parameters.fpp` and `m_start_up.fpp`
  (`GPU_DECLARE` alone does not make it device-resident). Regrouping scalars into a derived
  type silently drops their broadcast, so every non-root rank keeps the `dflt_real`
  sentinel; single-rank goldens cannot see this, so pair such a change with a `ppn=2` test
  and confirm it fails without the emitter.
- A `patch_ib` member that any `m_ibm` ghost-point code reads must ALSO be set in
  `s_add_cloud_particle` (`src/simulation/m_particle_cloud.fpp`): `particle_cloud_ibs` is
  `allocate`d without default initialization, and `s_reduce_ib_patch_array` copies the whole
  struct into `patch_ib`, overwriting the defaults from
  `s_assign_default_values_to_user_inputs`. Anything left unset reaches the solver as
  uninitialized memory, and only where the allocation is not already zero-filled — a
  garbage `v_blow` failed Frontier AMD with `ICFL is NaN` while every NVIDIA lane and all
  local CPU/GPU runs passed. A platform-only NaN is the signature of this class.
- Shared-state pattern: namelist declarations (`#:include 'generated_decls.fpp'`), the
  `eqn_idx`/`sys_size` state variables, and the common defaults
  core all live in `src/common/m_global_parameters_common.fpp`. Each per-target
  `m_global_parameters.fpp` does `use m_global_parameters_common` (default-public), so
  `use m_global_parameters` continues to work for all downstream modules without change.
  `src/common/` carries no `MFC_PRE_PROCESS`/`MFC_SIMULATION`/`MFC_POST_PROCESS` guards:
  stage-varying behavior is passed in as an explicit argument or initialization policy, and
  device residency for generated simulation scalars is emitted from `SIM_GPU_DECL_VARS`
  (`toolchain/mfc/params/generators/fortran_gen.py`). Generated includes
  (`generated_decls.fpp`, `generated_bcast.fpp`, `generated_case_opt_decls.fpp`) must exist for every target — pre/post
  get the common computed scalars (`num_dims`, `num_vels`, `weno_polyn`, `muscl_polyn`), so
  a common file that includes one will compile for pre/post too.
- Runtime checks (`@:PROHIBIT`) go where they run: shared →
  `src/common/m_checker_common.fpp`; simulation-only → `src/simulation/m_checker.fpp`;
  pre/post-only → `src/{pre,post}_process/m_checker.fpp` (their `s_check_inputs` are
  currently empty — that IS the right place, not m_checker_common).
- Analytic ICs are compiled into the binary. Expressions are AST-validated at case load
  (syntax errors and unknown variables are immediate, named errors; bare `e` is not a
  variable — write `exp(1.0)`).
  Each IC variable maps to an `eqn_idx%…` expression in `QPVF_IDX_VARS`
  (`toolchain/mfc/case.py`); a new patch-settable conserved variable means updating that
  map AND the Fortran `eqn_idx` builder to agree — a mismatch is a silent wrong index.
  Variables available in expressions: `docs/documentation/case.md`.

## Tests

- Tests are generated programmatically in `toolchain/mfc/test/cases.py` (parameter
  modifications on `BASE_CFG` via the `CaseGeneratorStack` push/pop pattern); test UUID =
  CRC32 of the trace string; `./mfc.sh test -l` lists all.
- `--only` matches whole trace *elements*, not substrings, and `_filter_only`
  (`toolchain/mfc/test/test.py`) **ANDs labels while ORing UUIDs**. So `--only bubbles` matches
  nothing (the element is `Bubbles`), and `--only low_Mach=1 low_Mach=2` asks for cases carrying
  both and also matches nothing. It then exits **143**, which reads like an external kill rather
  than an empty filter. Pass UUIDs whenever you want the union of several groups.
- Sibling `define_case_d` calls off the same stack level are never *combined*. Two switches that
  only matter together (`avg_state=1` needs `wave_speeds=2` to be read at all) therefore get zero
  effective coverage unless something pushes one and defines the other beneath it. Check
  reachability before trusting that a flag is tested.
- `--no-build` silently runs whatever binary is on disk for a configuration it did not build.
  Chemistry has its own config (`gpu-mp-chem-*`) that a plain `./mfc.sh build` never produces, so
  a `--no-build` run reports failures from stale binaries and hides real compile breaks. Run
  chemistry-touching sets without it.
- Pick the newest binary by the *binary's* mtime (`ls -t build/install/*/bin/simulation`), not the
  install directory's - a stale config's directory can be newer than a fresh build's.
- The pre-commit hook lives in the main repo's `.git/hooks/` and git exports `GIT_DIR` there
  during a commit, so from a worktree the toolchain lint enumerates the *other* checkout and
  fails. Reproduce with `GIT_DIR=<main>/.git ./mfc.sh precheck`. Run precheck by hand and commit
  with `--no-verify`.
- `/tmp` is node-local: scratch does not survive a compute-node change, and its absence is
  silence, not an error. Keep patches and resource baselines on a shared filesystem.
- Golden files are tolerance-compared. Regenerate only the affected tests
  (`./mfc.sh test --generate --only <tests>`) — an unexplained golden-file diff is a bug
  report, not noise to be regenerated away.
