# MFC — Multi-component Flow Code

MFC is an exascale multi-physics CFD solver written in modern Fortran 2008+ with Fypp
preprocessing. It has three executables (pre_process, simulation, post_process), a Python
toolchain for building/running/testing, and supports GPU acceleration via OpenACC and
OpenMP target offload. It must compile with gfortran, nvfortran, Cray ftn, and Intel ifx (CI-gated).
AMD flang is additionally supported for OpenMP target offload GPU builds.

## Commands

Prefer using `./mfc.sh` as the entry point for building, running, testing, formatting,
and linting. It handles virtual environments, module loading, dependency bootstrapping,
and build configuration. Avoid invoking CMake, Python toolchain scripts, or Fortran
compilers directly unless you have a specific reason.

All commands run from the repo root via `./mfc.sh`.

```bash
# Building
./mfc.sh build -j 8                        # Build all 3 targets (pre_process, simulation, post_process)
./mfc.sh build -t simulation -j 8          # Build only simulation
./mfc.sh build --gpu acc -j 8              # Build with OpenACC GPU support
./mfc.sh build --gpu mp -j 8              # Build with OpenMP target offload GPU support
./mfc.sh build --debug -j 8                # Debug build
./mfc.sh build -i case.py --case-optimization -j 8  # Case-optimized build (10x speedup)

# Running
./mfc.sh run case.py -n 4                  # Run case with 4 MPI ranks
./mfc.sh run case.py --no-build            # Run without rebuilding
./mfc.sh run case.py -e batch -N 2 -n 4 -c phoenix -a ACCOUNT  # Batch submit on Phoenix

# Testing
./mfc.sh test -j 8                         # Run full test suite (560+ tests)
./mfc.sh test --only 1D -j 8              # Only 1D tests
./mfc.sh test --only 2D Bubbles -j 8      # Only 2D bubble tests
./mfc.sh test --only <UUID> -j 8          # Run one specific test by UUID
./mfc.sh test -l                           # List all tests with UUIDs and traces
./mfc.sh test -% 10 -j 8                  # Run 10% random sample
./mfc.sh test --generate --only <feature>  # Regenerate golden files after intentional output change

# Verification (pre-commit CI checks)
./mfc.sh precheck -j 8                     # Run all 6 lint checks (same as CI gate)
./mfc.sh format -j 8                       # Auto-format Fortran (.fpp/.f90) + Python
./mfc.sh lint                              # Pylint + Python unit tests
./mfc.sh spelling                          # Spell check

# Module loading (HPC clusters only — must use `source`)
source ./mfc.sh load -c p -m g             # Load Phoenix GPU modules
source ./mfc.sh load -c f -m g             # Load Frontier GPU modules
source ./mfc.sh load -c p -m c             # Load Phoenix CPU modules

# Other
./mfc.sh validate case.py                  # Validate case file without running
./mfc.sh params <query>                    # Search 3,400 case parameters
./mfc.sh clean                             # Remove build artifacts
./mfc.sh new <name>                        # Create new case from template
```

## System Identification and Module Loading

MFC targets HPC clusters. Before building on a cluster, load the correct modules
via `source ./mfc.sh load -c <slug> -m <mode>`.

To identify the current system, check multiple signals — hostname alone is not always
sufficient (compute nodes may differ from login nodes):

```bash
hostname                    # e.g., login-phoenix-gnr-2.pace.gatech.edu
echo $LMOD_SYSHOST          # e.g., "phoenix" (most reliable when set)
echo $CRAY_LD_LIBRARY_PATH  # Non-empty → Cray system (Frontier, Carpenter Cray)
echo $MODULESHOME           # Confirms module system is available
```

Supported systems and their slugs (full list in `toolchain/modules`):

| Slug | System | GPU Backend | Example |
|------|--------|-------------|---------|
| `p` | GT Phoenix | OpenACC (nvfortran) | `source ./mfc.sh load -c p -m g` |
| `f` | OLCF Frontier | OpenACC/OpenMP (Cray ftn) | `source ./mfc.sh load -c f -m g` |
| `tuo` | LLNL Tuolumne | OpenMP (Cray ftn) | `source ./mfc.sh load -c tuo -m g` |
| `d` | NCSA Delta | OpenACC (nvfortran) | `source ./mfc.sh load -c d -m g` |
| `b` | PSC Bridges2 | OpenACC (nvfortran) | `source ./mfc.sh load -c b -m g` |
| `cc` | DoD Carpenter (Cray) | CPU only | `source ./mfc.sh load -c cc -m c` |
| `c` | DoD Carpenter (GNU) | CPU only | `source ./mfc.sh load -c c -m c` |
| `o` | Brown Oscar | OpenACC (nvfortran) | `source ./mfc.sh load -c o -m g` |
| `h` | UF HiPerGator | OpenACC (nvfortran) | `source ./mfc.sh load -c h -m g` |

The `-m` flag selects mode: `g`/`gpu` for GPU builds, `c`/`cpu` for CPU-only.
Batch job templates for `./mfc.sh run -e batch -c <system>` are in `toolchain/templates/`.

IMPORTANT: `source` (or `.`) is required for `load` — it sets environment variables
in the current shell. Using `./mfc.sh load` without `source` will error.

## Development Workflow Contract

IMPORTANT: Follow this loop for ALL code changes. Do not skip steps.

1. **Read first** — Read and understand relevant code before modifying it.
2. **Plan** — For multi-file changes, outline your approach before implementing.
3. **Implement** — Make small, focused changes. One logical change per commit.
4. **Format** — Run `./mfc.sh format -j 8` to auto-format.
5. **Verify** — Run `./mfc.sh precheck -j 8` (same 6 checks as CI lint gate).
6. **Build** — Run `./mfc.sh build -j 8` to verify compilation.
7. **Test** — Run relevant tests: `./mfc.sh test --only <feature> -j 8`.
   For changes to `src/common/`, test ALL three targets: `./mfc.sh test -j 8`.
8. **Commit** — Only after steps 4-7 pass. Do not commit untested code.

YOU MUST run `./mfc.sh precheck` before any commit. This is enforced by pre-commit hooks.
YOU MUST run tests relevant to your changes before claiming work is done.
NEVER commit code that does not compile or fails tests.

## Architecture

```
src/
  common/         # Shared code (used by ALL three executables — wide blast radius)
  pre_process/    # Grid generation and initial conditions
  simulation/     # CFD solver (GPU-accelerated via OpenACC / OpenMP target offload)
  post_process/   # Data output and visualization
toolchain/        # Python CLI, build system, testing, parameter management
  mfc/params/definitions.py   # ~3,400 parameter definitions (source of truth)
  mfc/case_validator.py       # Physics constraint validation
  mfc/test/                   # Test runner and case generation
examples/         # Example simulation cases (case.py files)
tests/            # 560+ regression test golden files
```

Source files are `.fpp` (Fortran + Fypp macros), preprocessed to `.f90` by CMake.

## Critical Rules

NEVER use raw OpenACC/OpenMP pragmas (`!$acc`, `!$omp`). Use `GPU_*` Fypp macros instead.
  Raw `#ifdef`/`#ifndef` preprocessor guards for feature/compiler/library gating ARE normal.
NEVER use double-precision intrinsics: `dsqrt`, `dexp`, `dlog`, `dble`, `dabs`, `real(8)`, `real(4)`.
  Use generic intrinsics (`sqrt`, `exp`, `log`) and precision types (`wp`, `stp`).
NEVER use `d` exponent literals (`1.0d0`). Use `1.0_wp` instead.
NEVER use `stop` or `error stop`. Use `call s_mpi_abort()` or `@:PROHIBIT()`/`@:ASSERT()`.
NEVER use `goto`, `COMMON` blocks, or global `save` variables.

Every `@:ALLOCATE(...)` MUST have a matching `@:DEALLOCATE(...)`.
Every new parameter MUST be added in at least 3 places (4 if it has constraints):
  1. `toolchain/mfc/params/definitions.py` (parameter definition)
  2. Fortran variable declaration in `src/*/m_global_parameters.fpp`
  3. Fortran namelist in `src/*/m_start_up.fpp` (namelist binding)
  4. `toolchain/mfc/case_validator.py` (only if parameter has physics constraints)

Changes to `src/common/` affect ALL three executables. Test comprehensively.

## Naming Conventions

- Modules: `m_<feature>` (e.g., `m_bubbles`)
- Public subroutines: `s_<verb>_<noun>` (e.g., `s_compute_pressure`)
- Public functions: `f_<verb>_<noun>`
- 2-space indentation, lowercase keywords, explicit `intent` on all arguments

## Precision System

- `wp` = working precision (computation). `stp` = storage precision (field data arrays and I/O).
- Default: both double. Single mode: both single. Mixed: wp=double, stp=half.
- MPI types must match: `mpi_p` ↔ `wp`, `mpi_io_p` ↔ `stp`.

## Code Review Priorities

When reviewing PRs, prioritize in this order:
1. Correctness (logic bugs, numerical issues, array bounds)
2. Precision discipline (stp vs wp mixing)
3. Memory management (@:ALLOCATE/@:DEALLOCATE pairing, GPU pointer setup)
4. MPI correctness (halo exchange, buffer sizing, GPU_UPDATE calls)
5. GPU code (GPU_* Fypp macros only, no raw pragmas)
6. Physics consistency (pressure formula matches model_eqns)
7. Compiler portability (4 CI-gated compilers + AMD flang for GPU)
