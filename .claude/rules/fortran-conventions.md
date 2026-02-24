# Fortran Conventions

## File Format
- Source files use `.fpp` extension (Fortran + Fypp preprocessor macros)
- Fypp preprocesses `.fpp` → `.f90` at build time via CMake
- Fypp supports conditional compilation, code generation, and regex macros

## Module Structure
Every Fortran module follows this pattern:
- File: `m_<feature>.fpp`
- Module: `module m_<feature>`
- `implicit none` required
- Explicit `intent(in)`, `intent(out)`, or `intent(inout)` on ALL subroutine/function arguments
- Initialization subroutine: `s_initialize_<feature>_module`
- Finalization subroutine: `s_finalize_<feature>_module`

## Naming
- Modules: `m_<feature>`
- Public subroutines: `s_<verb>_<noun>`
- Public functions: `f_<verb>_<noun>`
- Private/local variables: no prefix required
- Constants: descriptive names, not ALL_CAPS

## Forbidden Patterns

Caught by `./mfc.sh precheck` (source lint step 4/5):
- `dsqrt`, `dexp`, `dlog`, `dble`, `dabs` → use `sqrt`, `exp`, `log`, `real(..., wp)`
- `real(8)`, `real(4)` → use `wp` or `stp` kind parameters
- Raw `!$acc` or `!$omp` directives → use Fypp GPU_* macros from `parallel_macros.fpp`

Enforced by convention/code review (not automated):
- `goto`, `COMMON` blocks, global `save` variables
- `stop`, `error stop` → use `call s_mpi_abort()`

## Precision Types
- `wp` (working precision): used for computation. Double by default.
- `stp` (storage precision): used for I/O. Double by default.
- In single-precision mode (`--single`): both become single.
- In mixed-precision mode (`--mixed`): wp=double, stp=half.
- MPI type matching: `mpi_p` must match `wp`, `mpi_io_p` must match `stp`.
- Always use generic intrinsics: `sqrt` not `dsqrt`, `abs` not `dabs`.
- Cast with `real(..., wp)` or `real(..., stp)`, never `dble(...)`.

## Size Guidelines (soft)
- Subroutine: ≤500 lines
- Helper routine: ≤150 lines
- Function: ≤100 lines
- File: ≤1000 lines
- Arguments: ≤6 preferred
