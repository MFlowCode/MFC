# Parameter System

## Overview
MFC has ~3,400 simulation parameters defined in Python and read by Fortran via namelist files.

## Parameter Flow: Python → Fortran

1. **Definition**: `toolchain/mfc/params/definitions.py` — source of truth
   - Parameters are indexed families: `patch_icpp(i)%attr`, `fluid_pp(i)%attr`, etc.
   - Each has type, default, constraints, and tags

2. **Validation** (two layers):
   - `toolchain/mfc/case.py` / `toolchain/mfc/params/registry.py` — JSON schema validation
     via fastjsonschema (type checking, defaults)
   - `toolchain/mfc/case_validator.py` — Physics constraint checking
     (e.g., volume fractions sum to 1, dependency validation)

3. **Input Generation**: `toolchain/mfc/run/input.py`
   - Python case dict → Fortran namelist `.inp` file
   - Format: `&user_inputs` ... `&end/`

4. **Fortran Reading**: `src/*/m_start_up.fpp`
   - Reads `&user_inputs` namelist
   - Each parameter must be declared in the namelist statement

## Adding a New Parameter (4-location checklist)

YOU MUST update the first 3 locations. Missing any causes silent failures or compile errors.
Location 4 is required only if the parameter has physics constraints.

1. **`toolchain/mfc/params/definitions.py`**: Add parameter with type, default, constraints
2. **`src/*/m_global_parameters.fpp`**: Declare the Fortran variable in the relevant
   target(s). If the param is used by simulation only, add it there. If shared, add to
   all three targets' m_global_parameters.fpp.
3. **`src/*/m_start_up.fpp`**: Add to the Fortran `namelist` declaration in the relevant
   target(s).
4. **`toolchain/mfc/case_validator.py`**: Add validation rules if the parameter has
   physics constraints. Include `PHYSICS_DOCS` entry with title, category, explanation.

## Case Files
- Case files are Python scripts (`.py`) that define a dict of parameters
- Validated with `./mfc.sh validate case.py`
- Examples in `examples/` directory
- Create new cases with `./mfc.sh new <name>`
- Search parameters with `./mfc.sh params <query>`

## Fortran-Side Runtime Validation
Each target has `m_checker*.fpp` files (e.g., `src/simulation/m_checker.fpp`,
`src/common/m_checker_common.fpp`) containing runtime parameter validation using
`@:PROHIBIT(condition, message)`. When adding parameters with physics constraints,
add Fortran-side checks here in addition to `case_validator.py`.

## Analytical Initial Conditions
String expressions in parameters become Fortran code via `case.py.__get_analytic_ic_fpp()`.
These are compiled into the binary, so syntax errors cause build failures, not runtime errors.

Available variables in analytical IC expressions:
- `x`, `y`, `z` — cell-center coordinates (mapped to `x_cc(i)`, `y_cc(j)`, `z_cc(k)`)
- `xc`, `yc`, `zc` — patch centroid coordinates
- `lx`, `ly`, `lz` — patch lengths
- `r` — patch radius; `eps`, `beta` — vortex parameters
- `e` — Euler's number (2.71828...)
- Standard Fortran math intrinsics available: `sin`, `cos`, `exp`, `sqrt`, `abs`, etc.
- For moving immersed boundaries: `t` (simulation time) is also available

Example: `'patch_icpp(1)%vel(2)': '(x - xc) * exp(-((x-xc)**2 + (y-yc)**2))'`
