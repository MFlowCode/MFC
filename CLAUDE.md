# MFC — Multi-component Flow Code

MFC is an exascale multi-physics CFD solver written in modern Fortran 2008+ with Fypp
preprocessing. It has three executables (pre_process, simulation, post_process), a Python
toolchain for building/running/testing, and supports GPU acceleration via OpenACC and
OpenMP target offload. It must compile with gfortran, nvfortran, Cray ftn, and Intel ifx (CI-gated).
AMD flang is additionally supported for OpenMP target offload GPU builds.

## Working Style

Make surgical changes: every changed line should trace to the request. Don't refactor,
reformat, or "improve" adjacent code — in a four-compiler, golden-file-gated codebase,
incidental edits are how regressions slip in. For general behavioral guidance (simplicity,
surfacing assumptions, verifiable success criteria), invoke the `karpathy-guidelines` skill.

## Commands

Prefer using `./mfc.sh` as the entry point for building, running, testing, formatting,
and linting. It handles virtual environments, module loading, dependency bootstrapping,
and build configuration. Avoid invoking CMake, Python toolchain scripts, or Fortran
compilers directly unless you have a specific reason.

All commands run from the repo root via `./mfc.sh`.

Run `./mfc.sh <command> --help` for the full flag set; the most-used invocations:

```bash
# Build / run / test  (-j N = parallel jobs)
./mfc.sh build -j 8                 # all 3 targets; flags: -t <target>, --gpu acc|mp, --debug,
                                    #   -i case.py --case-optimization (10x speedup)
./mfc.sh run case.py -n 4           # run with 4 MPI ranks; --no-build; -e batch (toolchain/templates/)
./mfc.sh test -j 8                  # full suite (560+); --only <1D|Bubbles|UUID>, -l, -% N (sample),
                                    #   --generate (regenerate golden files after an intended output change)

# Verify before committing
./mfc.sh precheck -j 8              # all 7 CI lint checks
./mfc.sh format -j 8               # auto-format Fortran (.fpp/.f90) + Python
./mfc.sh lint                       # ruff lint + Python unit tests   (spelling: ./mfc.sh spelling)

# Case files
./mfc.sh validate case.py           # validate without running
./mfc.sh params <query>             # search 3,400 case parameters
./mfc.sh new <name>                 # new case from template       (clean: ./mfc.sh clean)
```

Module loading (`source ./mfc.sh load -c <slug> -m <mode>`) is covered under System Identification below.

## System Identification and Module Loading

On an HPC cluster, load modules before building: `source ./mfc.sh load -c <slug> -m <mode>`
(`-m g`/`gpu` or `c`/`cpu`). The `source` is required — plain `./mfc.sh load` errors, since
the command sets environment variables in the current shell.

Slugs live in `toolchain/modules` (e.g. `p` Phoenix, `f` Frontier, `tuo` Tuolumne, `d` Delta,
`b` Bridges2, `c`/`cc` Carpenter GNU/Cray, `o` Oscar, `h` HiPerGator; GPU backend per system
is defined there). To identify the current system, check `$LMOD_SYSHOST` (most reliable),
then a non-empty `$CRAY_LD_LIBRARY_PATH` (→ Cray: Frontier / Carpenter-Cray), then `hostname`
— login and compute nodes may differ. Batch templates for `./mfc.sh run -e batch -c <system>`
are in `toolchain/templates/`.

## Development Workflow Contract

IMPORTANT: Follow this loop for ALL code changes. Do not skip steps.

1. **Read first** — Read and understand relevant code before modifying it.
2. **Plan** — For multi-file changes, outline your approach before implementing.
3. **Implement** — Make small, focused changes. One logical change per commit.
4. **Format** — Run `./mfc.sh format -j 8` to auto-format.
5. **Verify** — Run `./mfc.sh precheck -j 8` (same 7 checks as CI lint gate).
6. **Build** — Run `./mfc.sh build -j 8` to verify compilation.
7. **Test** — Run relevant tests: `./mfc.sh test --only <feature> -j 8`.
   For changes to `src/common/`, test ALL three targets: `./mfc.sh test -j 8`.
8. **Commit** — Only after steps 4-7 pass. Do not commit untested code.

YOU MUST run `./mfc.sh precheck` before any commit. This is enforced by pre-commit hooks.
YOU MUST run tests relevant to your changes before claiming work is done.
NEVER commit code that does not compile or fails tests.
NEVER use heredocs for git commit messages. Use simple `git commit -m "message"` instead.

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
Every new parameter MUST be added in at least 2 places (3 if it has constraints):
  1. `toolchain/mfc/params/definitions.py` (parameter definition + NAMELIST_VARS target set)
  2. `toolchain/mfc/case_validator.py` (only if parameter has physics constraints)
  Note: Fortran declarations and namelist bindings are auto-generated from definitions.py
  at CMake configure time. Simple scalars need no manual Fortran edits. Array/derived-type
  variables still require a manual declaration in `src/*/m_global_parameters.fpp`.

Changes to `src/common/` affect ALL three executables. Test comprehensively.

## Naming Conventions

- Modules: `m_<feature>` (e.g., `m_bubbles`)
- Public subroutines: `s_<verb>_<noun>` (e.g., `s_compute_pressure`)
- Public functions: `f_<verb>_<noun>`
- Private/local variables: no prefix required. Constants: descriptive names, not ALL_CAPS.
- 2-space indentation, lowercase keywords, explicit `intent` on all arguments

## Precision System

- `wp` = working precision (computation). `stp` = storage precision (field data arrays and I/O).
- Both double by default. See `.claude/rules/fortran-conventions.md` for single/mixed
  modes, casting rules, and MPI type matching (`mpi_p` ↔ `wp`, `mpi_io_p` ↔ `stp`).

## Code Review Priorities

When reviewing PRs, prioritize in this order:
1. Correctness (logic bugs, numerical issues, array bounds)
2. Precision discipline (stp vs wp mixing)
3. Memory management (@:ALLOCATE/@:DEALLOCATE pairing, GPU pointer setup)
4. MPI correctness (halo exchange, buffer sizing, GPU_UPDATE calls)
5. GPU code (GPU_* Fypp macros only, no raw pragmas)
6. Physics consistency (pressure formula matches model_eqns)
7. Compiler portability (4 CI-gated compilers + AMD flang for GPU)
