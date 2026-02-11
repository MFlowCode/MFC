# AI Code Review Instructions for MFC

These instructions guide AI code reviewers (GitHub Copilot, CodeRabbit, Qodo, Cubic, CodeAnt, etc.) when evaluating pull requests in this repository.

**Coding standards, common pitfalls, and contribution guidelines:** `docs/documentation/contributing.md`
**GPU macro API reference:** `docs/documentation/gpuParallelization.md`

Formatting and linting are enforced by pre-commit hooks. Focus review effort on correctness, not style.

---

## 1  Project Context

* **Project:** MFC (exascale many-physics solver) written in **modern Fortran 2008+**, preprocessed with **Fypp**.
* **Directory layout:**
  * Sources in `src/`, tests in `tests/`, examples in `examples/`, Python toolchain in `toolchain/`.
  * Most source files are `.fpp` (Fypp templates); CMake transpiles them to `.f90`.
* **Fypp macros** are in `src/<subprogram>/include/`, where `<subprogram>` is `simulation`, `common`, `pre_process`, or `post_process`.
* Only `simulation` (plus its `common` dependencies) is GPU-accelerated via **OpenACC**.
* Code must compile with **GNU gfortran**, **NVIDIA nvfortran**, **Cray ftn**, and **Intel ifx**.
* Precision modes: double (default), single, and mixed (`wp` = working precision, `stp` = storage precision).
* **Python toolchain** requires **Python 3.10+** — do not suggest `from __future__` imports or other backwards-compatibility shims.

---

## 2  What to Review

See the **Common Pitfalls** section of `docs/documentation/contributing.md` for the full reference. Key review priorities:

1. **Correctness over style** — logic bugs, numerical issues, array bounds (non-unity lower bounds with ghost cells)
2. **Precision mixing** — `stp` vs `wp` conversions, no double-precision intrinsics (`dsqrt`, `dble`, `real(8)`, etc.)
3. **Memory** — `@:ALLOCATE`/`@:DEALLOCATE` pairing, GPU pointer setup (`@:ACC_SETUP_VFs`/`@:ACC_SETUP_SFs`)
4. **MPI** — halo exchange pack/unpack offsets, `GPU_UPDATE` around MPI calls, buffer sizing
5. **GPU** — no raw OpenACC/OpenMP pragmas (use Fypp GPU macros), `private(...)` on loop-local variables, no `stop`/`error stop` in device code
6. **Physics** — pressure formula must match `model_eqns`, `case_validator.py` constraints for new parameters
7. **Compiler portability** — all four compilers, Fypp macros for GPU and CPU builds
8. **`src/common/` blast radius** — changes affect all three executables

---

## 3  PR-Pattern Triggers

Flag these patterns when reviewing a pull request:

* PR adds a parameter in `toolchain/mfc/params/definitions.py` but does not update the `namelist /user_inputs/` in `src/*/m_start_up.fpp` or declare it in `src/*/m_global_parameters.fpp`
* PR adds a parameter with cross-parameter constraints but does not add validation in `toolchain/mfc/case_validator.py`
* PR modifies files in `src/common/` but does not mention testing all three targets (pre_process, simulation, post_process)
* PR adds `@:ALLOCATE` calls without matching `@:DEALLOCATE` in the corresponding finalization subroutine
* PR renames or moves files referenced in `docs/documentation/contributing.md` or `docs/documentation/gpuParallelization.md` without updating those docs
* PR adds a new `.fpp` module missing `implicit none` or missing `private`/`public` declarations
* PR modifies MPI pack/unpack logic in one sweep direction without updating all directions
* PR adds a new physics feature or model without a corresponding regression test
