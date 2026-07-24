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

## AMR levels (silent-index traps)

- **A level-`l` block's fine extent is `amr_ref_ratio**l * (coarse-region width) - 1`, NOT
  `amr_ref_ratio*width`.** The `amr_ref_ratio*width` form is correct only for the level-1 initial
  block; nested boxes compound by `amr_ref_ratio` per level (`amr_ref_ratio**level`). Every
  fine-extent computation uses `amr_ref_ratio**amr_block_level` — geometry
  (`s_set_amr_fine_geometry`), the restart-reader extent check, load-weight, `fmul`.
  Assuming `amr_ref_ratio*width` rejects level≥2 blocks as corrupt (the exact bug that bit the
  multi-level restart reader).
- **"coarse" in the AMR coupling routines means the block's PARENT level (`l-1`), not the
  base grid (level 0).** For a level-1 block the parent IS L0; for level≥2 the block folds
  to/from its parent block's fine array. `s_amr_gather_coarse_patch`,
  `s_interpolate_coarse_to_fine`, and the restrict/reflux path all operate in the
  parent-fine frame — assuming L0 silently corrupts level≥2 coupling.
- **The fine advance SWAPS the coarse grid globals (`m/n/p`, `idwint/idwbuff`, coords,
  `acoustic_source`, `ab_active`) to a fine block and restores them after — see the SWAP
  CONTRACT block at the `sw_*` declarations in `m_amr.fpp`.** Any module-level variable
  DERIVED from the grid that a kernel reads during the fine advance must be swapped there or
  refreshed per fine call at its use site; if it is `GPU_DECLARE`'d, its DEVICE copy must be
  refreshed too. A stale device copy of coarse bounds reads out of range on the fine grid
  under **CCE OpenACC only** (NVHPC/CCE-omp evaluate bounds host-side) — this was the `ab_int`
  regression, fixed by an unconditional `GPU_UPDATE` in `s_compute_rhs`. `amr_rvw` (cyl_coord
  radius weights) is the next candidate, currently safe only via a `m_checker.fpp` gate.
  A CPU-only or NVHPC-acc pass proves NOTHING here; this class is CCE-acc-specific.

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
- Shared-state pattern: namelist declarations (`#:include 'generated_decls.fpp'`), the
  `eqn_idx`/`sys_size`/`b_size`/`tensor_size` state variables, and the common defaults
  core all live in `src/common/m_global_parameters_common.fpp`. Each per-target
  `m_global_parameters.fpp` does `use m_global_parameters_common` (default-public), so
  `use m_global_parameters` continues to work for all downstream modules without change.
  Sim-only declarations (GPU_DECLARE, Re_idx allocation) stay in
  `m_global_parameters_common` behind `#ifdef MFC_SIMULATION`. Generated includes
  (`generated_decls.fpp`, `generated_bcast.fpp`, `generated_case_opt_decls.fpp`) must exist for every target — the build
  emits stubs where the content is sim-only, so a common file that includes one will
  compile for pre/post too.
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
- Golden files are tolerance-compared. Regenerate only the affected tests
  (`./mfc.sh test --generate --only <tests>`) — an unexplained golden-file diff is a bug
  report, not noise to be regenerated away.
