# GitHub Copilot Instructions – Pull-Request Reviews for MFC

These instructions guide **GitHub Copilot Code Review** and **Copilot Chat** when they evaluate pull requests in this repository.

---

## 1  Project Context (always include)

* **Project:** MFC (exascale many-physics solver) written in **modern Fortran 2008+**, generated with **Fypp**.
* **Directory layout:**  
  * Sources in `src/`, tests in `tests/`, examples in `examples/`.  
  * Most source files are templated `.fpp`; CMake transpiles them to `.f90`.
* **Fypp macros** are in `src/<subprogram>/include/`, where `<subprogram>` is `simulation`, `common`, `pre_process`, or `post_process`. Review these first.
* Only `simulation` (plus its `common` dependencies) is GPU-accelerated with **OpenACC**.

> **Copilot, when reviewing:**  
> * Treat the codebase as free-form Fortran 2008+ with `implicit none`, explicit `intent`, and standard intrinsics.  
> * Prefer `module … contains … subroutine foo()` over legacy constructs; flag uses of `COMMON`, file-level `include`, or other global state.

---

## 2  Incremental-Change Workflow

Copilot, when reviewing:
* Encourage small, buildable commits

---

## 3  Style & Naming Conventions (*.fpp / *.f90)

| Element | Rule |
|---------|------|
| Indentation | 2 spaces; continuation lines align beneath &. |
| Case | Lower-case keywords & intrinsics (do, end subroutine, …). |
| Modules | m_<feature> (e.g. m_transport). |
| Public subroutines | s_<verb>_<noun> (s_compute_flux). |
| Public functions | f_<verb>_<noun>. |
| Routine size | subroutine ≤ 500 lines, helper ≤ 150, function ≤ 100, file ≤ 1000. |
| Arguments | ≤ 6; else use a derived-type params struct. |
| Forbidden | goto (except legacy), COMMON, save globals. |
| Variables | Every arg has explicit intent; use dimension/allocatable/pointer as appropriate. |
| Errors | Call s_mpi_abort(<msg>), never stop or error stop. |

Copilot, when reviewing:
* Flag violations of any cell above.
* Suggest refactors when size or argument limits are exceeded.
* Ensure private helpers stay inside their defining module and avoid nested procedures.

---

## 4  OpenACC Guidelines (for GPU kernels)

Wrap tight loops:

```fortran
!$acc parallel loop gang vector default(present) reduction(...)
```

* Add collapse(n) when safe.
* Declare loop-local variables with private(...).
* Allocate large arrays with managed or move them into a persistent !$acc enter data region at start-up.
* Avoid stop/error stop inside device code.
* Code must compile with Cray ftn, NVIDIA nvfortran, GNU gfortran, and Intel ifx/ifort.

---

## 5  Review Checklist (what Copilot should verify)

1.	Buildability: PR includes build instructions or CI passes the staged build.
2.	Tests: Focused tests are added/updated.
3.	Style compliance: All rules in §3 are satisfied.
4.	Module hygiene: No new global state; proper namespacing.
5.	Performance: GPU code follows §4; no large host/device transfers in hot loops.
6.	Documentation: Updated in-code comments and, when needed, README or docs site.
7.	Regressions: No changes to outputs of golden tests without clear justification.
