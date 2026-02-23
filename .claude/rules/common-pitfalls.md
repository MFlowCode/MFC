# Common Pitfalls

## Array Bounds
- Arrays use non-unity lower bounds with ghost cells
- Riemann solver indexing: left state at `j`, right state at `j+1`
- Off-by-one errors in ghost cell regions are a common source of bugs

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

## Test Golden Files
- Tests compare output against golden files in `tests/<hash>/golden.txt`
- If your change intentionally modifies output, regenerate golden files:
  `./mfc.sh test --generate --only <affected_tests> -j 8`
- Do not regenerate ALL golden files unless you understand every output change
- Golden file diffs are compared with tolerance, not exact match

## PR Checklist
Before submitting a PR:
- [ ] `./mfc.sh format -j 8` (auto-format)
- [ ] `./mfc.sh precheck -j 8` (5 CI lint checks)
- [ ] `./mfc.sh build -j 8` (compiles)
- [ ] `./mfc.sh test --only <relevant> -j 8` (tests pass)
- [ ] If adding parameters: all 3 locations updated
- [ ] If modifying `src/common/`: all three targets tested
- [ ] If changing output: golden files regenerated for affected tests
- [ ] One logical change per commit
