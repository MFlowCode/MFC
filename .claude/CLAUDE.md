This file is the single entry point for contributor and review instructions (`AGENTS.md` and `.github/copilot-instructions.md` point here).
It holds policy and pointers only; facts live in the linked docs, and `./mfc.sh precheck` fails if a path or command below goes stale.

* New branches cannot be made on MFlowCode/MFC, they are made on forks
* PRs:
    * are made from those MFC forks
    * that change CFD results need verification that the PR is correct
    * follow template (`.github/pull_request_template.md`)
    * that break a feature but promise a followup PR to fix it are rejected
* Commands:
    * MFC should almost always build and run using ./mfc.sh; `./mfc.sh <command> --help` lists options
    * Load a compiler and MPI before building. On a supported cluster, `source ./mfc.sh load -c <computer> -m <g|c>` loads a known-working module set. Pass both flags; omitting either waits for interactive input. Computer slugs are the first column of `toolchain/modules`; elsewhere, `module load` the equivalents
    * Agent shell tools run each command in a fresh shell, so a one-off `module load` or `source ./mfc.sh load` is gone by the next command. Chain it (`source ./mfc.sh load -c p -m g && ./mfc.sh build -j 8`) or put the steps in a script
    * Running mfc.sh commands can create a lock file in build/ that is sticky, be careful (`toolchain/mfc/lock.py`)
    * Run `./mfc.sh precheck` before pushing; it runs the CI lint gate except the slow rendering and Fortran-compile unit tests
    * Iterate small: `./mfc.sh build -t simulation` builds one target; `./mfc.sh test --only <UUID>` runs one test; `./mfc.sh test --select-enforce` runs only the tests whose coverage touches your changed files
* Programming and Design:
    * New code should follow the DRY principle and also make side-effect code DRY as well
    * Comments should be as short as possible without sacrificing value
    * GPU macros should follow the source's existing GPU macro principles and patterns
    * Functions/subroutines/modules shorter is better while being correct, fast, and separating concerns
    * When a fix fails twice, stop and re-read in full the code you are changing, its callers and callees, and where each global it uses is set, before the next attempt. Repeated failures usually come from a detail remembered wrong

## Where to look

| Topic | Read |
|---|---|
| Layout, coding standards, how-to guides | `docs/documentation/contributing.md` |
| Review checklist | Common Pitfalls in `docs/documentation/contributing.md` |
| GPU macros, compiler traps, NVHPC/Cray/AMD issues | `docs/documentation/gpuParallelization.md` (macros: `src/common/include/parallel_macros.fpp`) |
| GPU code that passes on CPU but fails on GPU | Silent-Failure Traps in `docs/documentation/gpuParallelization.md`; read before any GPU change |
| Tests and golden files | `docs/documentation/testing.md` (cases: `toolchain/mfc/test/cases.py`) |
| Adding a parameter | How-to in `docs/documentation/contributing.md`; `toolchain/mfc/params/definitions.py`, `toolchain/mfc/case_validator.py` |
| Build, run and test failures | `docs/documentation/troubleshooting.md` |
| CI jobs | `.github/workflows/test.yml` |
