This file is the single entry point for agents and AI reviewers (`AGENTS.md` and `.github/copilot-instructions.md` point here).
It holds policy and pointers only; facts live in the linked docs, and `./mfc.sh precheck` fails if a path or command below goes stale.

* New branches cannot be made on MFlowCode/MFC, they are made on forks
* PRs:
    * made using AI tools like Claude Code and Codex should say so.
    * are made from those MFC forks
    * that change CFD results need verification that the PR is correct
    * follow template (`.github/pull_request_template.md`)
    * that break a feature but promise a followup PR to fix it are rejected
* Commands:
    * MFC should almost always build and run using ./mfc.sh; `./mfc.sh <command> --help` lists options
    * Running mfc.sh commands can create a lock file in build/ that is sticky, be careful (`toolchain/mfc/lock.py`)
    * Run `./mfc.sh precheck` before pushing; it runs the same lint gate as CI
* Programming and Design:
    * New code should follow the DRY principle and also make side-effect code DRY as well
    * Comments should be as short as possible without sacrificing value
    * GPU macros should follow the source's existing GPU macro principles and patterns
    * Functions/subroutines/modules shorter is better while being correct, fast, and separating concerns

## Where to look

| Topic | Read |
|---|---|
| Layout, coding standards, how-to guides | `docs/documentation/contributing.md` |
| Review checklist (humans and AI reviewers) | Common Pitfalls in `docs/documentation/contributing.md` |
| GPU macros, compiler traps, NVHPC/Cray/AMD issues | `docs/documentation/gpuParallelization.md` (macros: `src/common/include/parallel_macros.fpp`) |
| Tests and golden files | `docs/documentation/testing.md` (cases: `toolchain/mfc/test/cases.py`) |
| Adding a parameter | How-to in `docs/documentation/contributing.md`; `toolchain/mfc/params/definitions.py`, `toolchain/mfc/case_validator.py` |
| Build, run and test failures | `docs/documentation/troubleshooting.md` |
| CI jobs | `.github/workflows/test.yml` |
