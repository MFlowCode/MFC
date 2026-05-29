# Common Pitfalls

## Array Bounds & Ghost Cells
- Grid dimensions: `m`, `n`, `p` (cells in x, y, z). 1D: n=p=0. 2D: p=0.
- Interior domain: `0:m`, `0:n`, `0:p`
- Buffer/ghost region: `-buff_size:m+buff_size` (similar for n, p in multi-D)
- `buff_size` is set per reconstruction scheme in `s_configure_coordinate_bounds` (`m_helper_basic.fpp`):
  WENO `weno_polyn + 2` (`2*weno_polyn + 2` if `viscous`), MUSCL `muscl_polyn + 2`, IGR `(igr_order-1)/2 + 2`;
  then floored to `max(.,6)` with Lagrange bubbles and `max(.,10)` with IB. Don't assume `2*weno_polyn + 2`.
- Domain bounds: `idwint(1:3)` (interior `0:m`), `idwbuff(1:3)` (with ghost cells)
- Cell-center coords: `x_cc(-buff_size:m+buff_size)`, `y_cc(...)`, `z_cc(...)`
- Cell-boundary coords: `x_cb(-1-buff_size:m+buff_size)`
- Riemann solver indexing: left state at `j`, right state at `j+1`
- Off-by-one errors in ghost cell regions are a common source of bugs

## Field Variable Indexing
- Conserved variables: `q_cons_vf(1:sys_size)`. Primitive: `q_prim_vf(1:sys_size)`.
- All equation indices live in the unified `eqn_idx` struct (`eqn_idx_info` type in `m_derived_types.fpp`).
  Index ranges depend on `model_eqns` and enabled features (set in `m_global_parameters.fpp`):
  - `eqn_idx%cont` — continuity range (partial densities, one per fluid)
  - `eqn_idx%mom` — momentum range
  - `eqn_idx%E` — total energy (scalar)
  - `eqn_idx%adv` — volume fractions (advection equations)
  - `eqn_idx%bub`, `eqn_idx%stress`, `eqn_idx%xi`, `eqn_idx%species`, `eqn_idx%B` — optional
  - `eqn_idx%gamma`, `eqn_idx%pi_inf`, `eqn_idx%alf`, `eqn_idx%int_en` — additional scalars/ranges
- Use `eqn_idx%cont%beg`/`eqn_idx%cont%end`, `eqn_idx%mom%beg`/`eqn_idx%mom%end`, etc. (old `contxb`/`contxe`, `momxb`/`momxe` shorthands are gone)
- `sys_size` = total number of conserved variables (computed at startup)
- Changing `model_eqns` or enabling features changes ALL index positions

## Blast Radius
- `src/common/` is shared by ALL three executables (pre_process, simulation, post_process)
- Any change to common/ requires testing all three targets
- Public subroutine signature changes affect all callers across all targets
- Parameter default changes affect all existing case files

## Physics Consistency
- Pressure formula MUST match `model_eqns` setting
- Model-specific conservative ↔ primitive conversion paths exist
- Volume fractions must sum to 1.0
- Boundary condition symmetry requirements must be maintained

## Compiler-Specific Issues
- See the compiler-backend matrix in `.claude/rules/gpu-and-mpi.md` for which compilers
  are CI-gated and which backends each supports.
- Each compiler has different strictness levels and warning behavior
- Fypp macros must expand correctly for both GPU and CPU builds

## Test System
- Tests are generated **programmatically** in `toolchain/mfc/test/cases.py`, not standalone files
- Each test is a parameter modification on top of `BASE_CFG` defaults
- Test UUID = CRC32 hash of the test's trace string; `./mfc.sh test -l` lists all
- To add a test: modify `cases.py` using `CaseGeneratorStack` push/pop pattern
- Golden files: `tests/<UUID>/golden.txt` — tolerance-based comparison, not exact match
- If your change intentionally modifies output, regenerate golden files:
  `./mfc.sh test --generate --only <affected_tests> -j 8`
- Do not regenerate ALL golden files unless you understand every output change

## PR Checklist
The base loop (format → precheck → build → test → one logical commit) is the
Development Workflow Contract in `CLAUDE.md`. Beyond it, watch for:
- [ ] If adding parameters: definitions.py (_r + _nv) updated; cmake reconfigured; case_validator.py if constraints
- [ ] If modifying `src/common/`: all three targets tested
- [ ] If changing output: golden files regenerated for affected tests (`./mfc.sh test --generate --only <tests>`)
