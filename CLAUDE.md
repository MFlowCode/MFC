# MFC — Multi-component Flow Code

MFC is an exascale multi-physics CFD solver written in modern Fortran 2008+ with Fypp
preprocessing. It has three executables (pre_process, simulation, post_process), a Python
toolchain for building/running/testing, and supports GPU acceleration via OpenACC and
OpenMP target offload. It must compile with gfortran, nvfortran, Cray ftn, and Intel ifx (CI-gated).
AMD flang is additionally supported for OpenMP target offload GPU builds.

## How to Work Here

You are editing this repository as a conservative maintainer. Three facts shape everything:

1. **Four compilers, one truth.** Every line must compile and behave identically under four
   CI-gated compilers and three GPU configurations. Code that is merely clever on one is
   broken on another. Prefer the established idiom over the elegant one.
2. **Failures here are silent.** The worst bugs in a CFD code are not crashes — they are
   answers that are subtly wrong: a wrong index, a serial loop that should be parallel, a
   precision mix. Golden-file regression tests are the safety net, and incidental edits are
   how regressions slip past them.
3. **Trust the toolchain.** `./mfc.sh` already solves environments, module loading,
   dependencies, and build configuration; the lint (`./mfc.sh precheck`) encodes the
   project's forbidden patterns. When unsure of a flag or rule, ask the tooling
   (`--help`, precheck) rather than guessing.

**Primary rule: make the smallest correct change. Prefer deleting code to adding code.**

Diff discipline:
- Every changed line should trace to the request. Keep changes localized.
- Do not rewrite, reformat, or rename unrelated code.
- Do not add new files, dependencies, or net-positive LOC without clear justification.

Do not add bloat:
- abstractions for one call site, or defensive code for hypothetical future callers
- validation layers around trusted internal data — constraint checks belong in
  `case_validator.py` and the `m_checker*.fpp` files, not scattered through the solver
- optional parameters or config knobs for behavior with one correct value
- broad try/except blocks in toolchain Python
- comments explaining obvious code; compatibility shims unless explicitly requested

Tests: add one only when it protects real behavior — it would fail before your change and
covers behavior a real case depends on. Prefer one targeted case over many broad ones.

Before editing, state: the smallest viable fix, what changes, what deliberately does not
change, and whether a test is needed. After editing, report: files changed, net LOC, and
anything that could still be deleted or simplified. If a patch exceeds roughly 100 net new
lines, stop and justify before continuing.

For general behavioral guidance (simplicity, surfacing assumptions, verifiable success
criteria), invoke the `karpathy-guidelines` skill.

## Commands

All commands run from the repo root via `./mfc.sh`; avoid invoking CMake, Python toolchain
scripts, or Fortran compilers directly. Run `./mfc.sh <command> --help` for flags
(`-j N` = parallel jobs).

```bash
./mfc.sh build -j 8                 # all 3 targets; -t <target>, --gpu acc|mp, --debug
./mfc.sh run case.py -n 4           # run with 4 MPI ranks; -e batch for clusters
./mfc.sh test -j 8                  # full suite; --only <filter>, -l to list,
                                    #   --generate to refresh golden files after an intended output change
./mfc.sh format -j 8                # auto-format Fortran + Python
./mfc.sh precheck -j 8              # all CI lint checks — run before every commit
./mfc.sh validate case.py           # validate a case without running
./mfc.sh params <query>             # search case parameters
```

On HPC clusters, load modules before building: `source ./mfc.sh load -c <slug> -m <g|c>`
(the `source` is required). Slugs and per-system GPU backends: `toolchain/modules`; batch
templates: `toolchain/templates/`. Identify the system via `$LMOD_SYSHOST`, then a
non-empty `$CRAY_LD_LIBRARY_PATH` (→ Cray), then `hostname`.

## Development Workflow Contract

For ALL code changes: read the relevant code first, plan multi-file changes before
implementing, then **format → precheck → build → test → commit**, in that order, with one
logical change per commit.

- YOU MUST run `./mfc.sh precheck` before any commit (pre-commit hooks enforce this).
- YOU MUST run the tests relevant to your change before claiming work is done.
  Changes to `src/common/` are shared by all three executables — test all three.
- NEVER commit code that does not compile or fails tests.
- NEVER use heredocs for git commit messages. Use simple `git commit -m "message"`.

## Architecture

```
src/common/       # shared by ALL three executables — wide blast radius
src/pre_process/  # grid generation and initial conditions
src/simulation/   # CFD solver (the only GPU-accelerated target)
src/post_process/ # data output and visualization
toolchain/        # Python CLI; params/definitions.py is the parameter source of truth
examples/         # example cases (case.py); tests/ holds regression golden files
```

Source files are `.fpp` (Fortran + Fypp macros), preprocessed to `.f90` by CMake.

## Critical Rules

The exhaustive forbidden-pattern list is `toolchain/mfc/lint_source.py` (enforced by
precheck and CI). The rules to internalize — they exist because MFC must behave
identically across compilers, precisions, and GPU backends:

- GPU directives only via `GPU_*` Fypp macros — NEVER raw `!$acc`/`!$omp` pragmas.
  (Raw `#ifdef`/`#ifndef` guards for feature/compiler/library gating ARE normal.)
- Precision only via `wp`/`stp` kinds and generic intrinsics — NEVER `dsqrt`/`dble`/
  `real(8)` or `d` exponent literals. Write `sqrt`, `1.0_wp`, `real(..., wp)`.
- Abort only via `call s_mpi_abort()` or `@:PROHIBIT()`/`@:ASSERT()` — NEVER `stop`/`error stop`.
- NEVER `goto`, `COMMON` blocks, or global `save` variables.
- Every `@:ALLOCATE(...)` MUST have a matching `@:DEALLOCATE(...)`.
- New parameters are defined in `toolchain/mfc/params/definitions.py` (plus
  `case_validator.py` if physics-constrained); Fortran declarations and namelist bindings
  are auto-generated. Exceptions: `.claude/rules/common-pitfalls.md`.

## Naming and Style

Modules `m_<feature>` with `s_initialize_/s_finalize_<feature>_module` pairs; public
subroutines `s_<verb>_<noun>`, functions `f_<verb>_<noun>`; 2-space indent, lowercase
keywords, explicit `intent` on all arguments; constants get descriptive names, not
ALL_CAPS. Full hard/soft rule tables: `docs/documentation/contributing.md`.

## Precision

`wp` = working precision (computation); `stp` = storage precision (field arrays and I/O).
Both double by default; `--single` → both single; `--mixed` → wp=double, stp=half — so
wp/stp mixing is a silent precision bug, not a style issue. `scalar_field%sf` and all new
field arrays use `stp`. MPI types must match: `mpi_p` ↔ `wp`, `mpi_io_p` ↔ `stp`.

## Where Things Are Documented

`docs/documentation/` is freshness-checked by precheck (`lint_docs.py`) — prefer pointing
there over restating it. Most relevant: `contributing.md` (standards, architecture,
general pitfalls), `gpuParallelization.md` (GPU macro API), `testing.md` (test system),
`case.md` (case parameters, analytic ICs). MFC-specific traps with silent failure modes
live in `.claude/rules/common-pitfalls.md` — read it before touching indexing, GPU loops,
parameters, or tests.

## Code Review Priorities

When reviewing PRs, prioritize in this order:
1. Correctness (logic bugs, numerical issues, array bounds)
2. Precision discipline (stp vs wp mixing)
3. Memory management (@:ALLOCATE/@:DEALLOCATE pairing, GPU pointer setup)
4. MPI correctness (halo exchange, buffer sizing, GPU_UPDATE calls)
5. GPU code (GPU_* Fypp macros only, no raw pragmas)
6. Physics consistency (pressure formula matches model_eqns)
7. Compiler portability (4 CI-gated compilers + AMD flang for GPU)
