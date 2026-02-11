@page contributing Contributing

# Contributing to MFC

We welcome contributions of all kinds: bug fixes, new features, documentation, tests, and issue triage.
This guide covers everything you need to get started and get your changes merged.

## Getting Set Up

1. **Fork and clone**
   ```bash
   git clone https://github.com/<your-user>/MFC.git
   cd MFC
   git remote add upstream https://github.com/MFlowCode/MFC.git
   ```
2. **Build MFC** (see @ref getting-started for full details):
   ```bash
   ./mfc.sh build -j $(nproc)
   ```
3. **Run the test suite** to verify your environment:
   ```bash
   ./mfc.sh test -j $(nproc)
   ```

## Development Workflow

| Step | Command / Action |
|------|-----------------|
| Sync your fork | `git checkout master && git pull upstream master` |
| Create a branch on your fork | `git checkout -b feature/<short-name>` |
| Code, test, document | Follow the standards below |
| Run tests | `./mfc.sh test` |
| Commit | Clear, atomic commits (see below) |
| Push to your fork | `git push origin feature/<short-name>` |
| Open a PR | From your fork to `MFlowCode/MFC:master`. Every push triggers CI -- bundle changes to avoid flooding the queue |

### Commit Messages

- Start with a concise (50 chars or fewer) summary in imperative mood: `Fix out-of-bounds in EOS module`
- Add a blank line, then a detailed explanation if needed
- Reference related issues: `Fixes #123` or `Part of #456`

## Coding Standards

MFC is written in modern Fortran 2008+ with [Fypp](https://github.com/aradi/fypp) metaprogramming.
The standards below are split into **hard rules** (enforced in CI and review) and **soft guidelines** (goals for new code).

### Hard Rules

These are enforced. CI and reviewers will flag violations.

| Element | Rule |
|---------|------|
| Formatting | Enforced automatically by pre-commit hook (`./mfc.sh format` and `./mfc.sh lint`) |
| Indentation | 2 spaces; continuation lines align beneath `&` |
| Case | Lowercase keywords and intrinsics (`do`, `end subroutine`, ...) |
| Module naming | `m_<feature>` (e.g. `m_riemann_solvers`) |
| Public subroutines | `s_<verb>_<noun>` (e.g. `s_compute_flux`) |
| Public functions | `f_<verb>_<noun>` (e.g. `f_create_bbox`) |
| Variables | Every argument has explicit `intent`; use `implicit none` |
| Forbidden | `goto`, `COMMON` blocks, global `save` variables |
| Error handling | Call `s_mpi_abort(<msg>)` -- never `stop` or `error stop` |
| Compiler support | Code must compile with GNU gfortran, NVIDIA nvfortran, Cray ftn, and Intel ifx |

### Soft Guidelines

Aim for these in new and modified code. Existing code may not meet all of them.

| Element | Guideline |
|---------|-----------|
| File length | Prefer 1000 lines or fewer per file |
| Subroutine length | Prefer 500 lines or fewer (helpers: 150) |
| Function length | Prefer 100 lines or fewer |
| Arguments | Prefer 6 or fewer; consider a derived-type struct for more |
| DRY | Avoid duplicating logic; factor shared code into helpers |

## Fypp and GPU

MFC uses Fypp macros (in `src/*/include/`) to generate accelerator-specific Fortran for OpenACC and OpenMP backends.
Key points:

- Use the project's GPU macros (`@:GPU_PARALLEL_LOOP`, `@:ALLOCATE`, etc.) -- raw OpenACC/OpenMP pragmas are not allowed
- Keep macros simple and readable
- Only `simulation` (plus its `common` dependencies) is GPU-accelerated
- See @ref gpuParallelization for the full macro API reference

## Testing

MFC has 500+ regression tests. See @ref testing for the full guide.

- **Add tests** for any new feature or bug fix
- Use `./mfc.sh test --generate` to create golden files for new cases
- Keep tests fast: use small grids and short runtimes
- Test with `--test-all` to include post-processing validation

## Documentation

- Add or update **Doxygen docstrings** in source files for new public routines
- Update **markdown docs** under `docs/` if user-facing behavior changes
- Provide a minimal **example case** in `examples/` for new features when practical

## Submitting a Pull Request

1. **PRs come from your fork.** Do not create branches on `MFlowCode/MFC` directly. Push to your fork and open a PR from there against `MFlowCode/MFC:master`.
2. **One PR = one logical change.** Split large changes into focused PRs.
3. **Fill out the PR template.** Remove checklist items that don't apply.
4. **Link issues** with `Fixes #<id>` or `Part of #<id>`.
5. **Ensure CI passes** before requesting review. Run `./mfc.sh test` locally first. Formatting and linting are handled automatically by the pre-commit hook.
6. **Describe your testing**: what you ran, which compilers/platforms you used.

If your change touches GPU code (`src/simulation/`), see the GPU checklist in the PR template.

## Code Review and Merge

- Respond to reviewer comments promptly
- Push focused updates; each push re-runs CI, so batch your fixes
- A maintainer will merge your PR once all reviews are approved and CI is green

If your PR is large or architectural, consider opening an issue first to discuss the approach.
