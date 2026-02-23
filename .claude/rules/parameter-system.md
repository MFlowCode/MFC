# Parameter System

## Overview
MFC has ~3,300 simulation parameters defined in Python and read by Fortran via namelist files.

## Parameter Flow: Python → Fortran

1. **Definition**: `toolchain/mfc/params/definitions.py` — source of truth
   - Parameters are indexed families: `patch_icpp(i)%attr`, `fluid_pp(i)%attr`, etc.
   - Each has type, default, constraints, and tags

2. **Validation**: `toolchain/mfc/case_validator.py`
   - JSON schema validation via fastjsonschema
   - Physics constraint checking (e.g., volume fractions sum to 1)
   - Dependency validation (required/recommended params)

3. **Input Generation**: `toolchain/mfc/run/input.py`
   - Python case dict → Fortran namelist `.inp` file
   - Format: `&user_inputs / ... / &end`

4. **Fortran Reading**: `src/*/m_start_up.fpp`
   - Reads `&user_inputs` namelist
   - Each parameter must be declared in the namelist statement

## Adding a New Parameter (3-location checklist)

YOU MUST update all 3 locations. Missing any causes silent failures.

1. **`toolchain/mfc/params/definitions.py`**: Add parameter with type, default, constraints
2. **`src/*/m_start_up.fpp`**: Add to the Fortran `namelist` declaration in the relevant
   target(s). If the param is used by simulation only, add it there. If shared, add to
   all three targets' m_start_up.fpp.
3. **`toolchain/mfc/case_validator.py`**: Add validation rules if the parameter has
   physics constraints. Include `PHYSICS_DOCS` entry with title, category, explanation.

## Case Files
- Case files are Python scripts (`.py`) that define a dict of parameters
- Validated with `./mfc.sh validate case.py`
- Examples in `examples/` directory
- Create new cases with `./mfc.sh new <name>`
- Search parameters with `./mfc.sh params <query>`

## Analytical Initial Conditions
String expressions in parameters become Fortran code via `case.py.__get_analytic_ic_fpp()`.
These are compiled into the binary, so syntax errors cause build failures, not runtime errors.
