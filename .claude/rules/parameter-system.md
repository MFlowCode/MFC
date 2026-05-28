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

## Adding a New Parameter (2-location checklist)

Fortran declarations and namelist bindings are now auto-generated from definitions.py
at CMake configure time — no manual Fortran edits needed for simple scalar parameters.

1. **`toolchain/mfc/params/definitions.py`**: Add parameter with `_r()` (type, default,
   constraints) AND add it to `NAMELIST_VARS` via `_nv()` for the relevant target(s).
   After editing, re-run cmake (or `./mfc.sh build`) to regenerate the Fortran includes.
2. **`toolchain/mfc/case_validator.py`**: Add validation rules if the parameter has
   physics constraints. Include `PHYSICS_DOCS` entry with title, category, explanation.

**Exceptions — still require manual Fortran edits:**
- Array variables (e.g. `logical, dimension(num_fluids_max)`) → declare in `src/*/m_global_parameters.fpp`
- Derived-type members (`fluid_pp%attr`, `patch_icpp(i)%attr`) → declare in the relevant derived type
- Case-optimization parameters → add to `CASE_OPT_PARAMS` and the `#:else` block in `src/simulation/m_global_parameters.fpp`

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
