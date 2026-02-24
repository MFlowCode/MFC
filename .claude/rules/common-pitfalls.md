# Common Pitfalls

## Array Bounds & Ghost Cells
- Grid dimensions: `m`, `n`, `p` (cells in x, y, z). 1D: n=p=0. 2D: p=0.
- Interior domain: `0:m`, `0:n`, `0:p`
- Buffer/ghost region: `-buff_size:m+buff_size` (similar for n, p in multi-D)
- `buff_size` depends on WENO order and features (typically `2*weno_polyn + 2`)
- Domain bounds: `idwint(1:3)` (interior `0:m`), `idwbuff(1:3)` (with ghost cells)
- Cell-center coords: `x_cc(-buff_size:m+buff_size)`, `y_cc(...)`, `z_cc(...)`
- Cell-boundary coords: `x_cb(-1-buff_size:m+buff_size)`
- Riemann solver indexing: left state at `j`, right state at `j+1`
- Off-by-one errors in ghost cell regions are a common source of bugs

## Field Variable Indexing
- Conserved variables: `q_cons_vf(1:sys_size)`. Primitive: `q_prim_vf(1:sys_size)`.
- Index ranges depend on `model_eqns` and enabled features (set in `m_global_parameters.fpp`):
  - `cont_idx` — continuity (partial densities, one per fluid)
  - `mom_idx` — momentum components
  - `E_idx` — total energy (scalar)
  - `adv_idx` — volume fractions (advection equations)
  - `bub_idx`, `stress_idx`, `xi_idx`, `species_idx`, `B_idx`, `c_idx` — optional
- Shorthand scalars: `momxb`/`momxe`, `contxb`/`contxe`, `advxb`/`advxe`, etc.
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
- Code must compile on gfortran, nvfortran, Cray ftn, and Intel ifx
- Each compiler has different strictness levels and warning behavior
- Fypp macros must expand correctly for both GPU and CPU builds
- GPU builds only work with nvfortran, Cray ftn, and AMD flang

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
Before submitting a PR:
- [ ] `./mfc.sh format -j 8` (auto-format)
- [ ] `./mfc.sh precheck -j 8` (5 CI lint checks)
- [ ] `./mfc.sh build -j 8` (compiles)
- [ ] `./mfc.sh test --only <relevant> -j 8` (tests pass)
- [ ] If adding parameters: all 4 locations updated
- [ ] If modifying `src/common/`: all three targets tested
- [ ] If changing output: golden files regenerated for affected tests
- [ ] One logical change per commit
