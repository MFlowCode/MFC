# Contributing to the MFC Codebase (Multi‑Component Flow Code)

**Multi‑Component Flow Code (MFC)** is an open‑source, high‑performance code for simulating compressible multi‑component, multi‑phase flows.
We welcome contributions of all kinds—bug fixes, new features, documentation, tests, and issue triage—from both newcomers and experienced developers.
This guide explains how to set up your environment, follow MFC's coding standards, and navigate the pull-request (PR) process so your work can be merged smoothly.

---

## 1. Setting Up Your Development Environment

1. **Fork and clone**
   ```bash
   git clone https://github.com/<your‑user>/MFC.git
   cd MFC
   git remote add upstream https://github.com/MFlowCode/MFC.git
   ```
2. **Build MFC** – follow the [documentation](https://mflowcode.github.io/documentation/md_getting-started.html). For example:
   ```bash
   ./mfc.sh build -j 8   # parallel build with 8 threads
   ```
3. **Run the test suite** to verify your environment:
   ```bash
   ./mfc.sh test -j 8
   ```

---

## 2. Development Workflow

| Step | Action | Notes |
|------|--------|-------|
| 1 | **Sync your fork**: `git checkout master && git pull upstream master` | Stay up‑to‑date to avoid merge conflicts. |
| 2 | **Create a branch**: `git checkout -b feature/<short‑name>` | Keep each branch focused on one logical change. |
| 3 | **Code, test, document** | Follow the guidelines in §3. |
| 4 | **Format & lint**: `./mfc.sh format` | CI will re‑check; make it pass locally first. |
| 5 | **Run tests**: `./mfc.sh test` | All existing and new tests must pass. |
| 6 | **Commit** (see *Commit Messages* below) | Write clear, atomic commits. |
| 7 | **Push** & open a **PR** | Be mindful: *every push triggers CI*. Bundle fixes together to avoid dozens of CI runs. |

### Commit Messages

- Start with a concise (≤50 chars) summary in imperative mood: `Fix out‑of‑bounds in EOS module`.
- Add a blank line, then a detailed explanation.
- Reference related issues or PRs, e.g., `Fixes #123`.

### Managing CI Runs

Each push to a branch with an open PR runs the full CI matrix (which can take hours).
Plan your pushes—run tests locally and group changes—so the CI queue is not flooded.

---

## 3. Coding Guidelines and Best Practices

### 3.1 Style, Formatting & Linting
MFC enforces a project‑wide Fortran style:
- **Formatter**: `./mfc.sh format` auto‑formats your changes.
- **Linter**: CI runs several linter checks that spot common Fortran-gotchas (implicit typing, shadowed variables, unused locals, etc.). Fix issues before pushing or the linter will often catch them.

### 3.2 Fypp Metaprogramming

MFC uses [**Fypp**](https://github.com/aradi/fypp), a lightweight Python-based preprocessor, to generate repetitive or accelerator-specific Fortran.
Key points:
- Fypp macros live in `include/` directories nested within `src/`.
- Run `./mfc.sh format` to format the example case files and the source code.
- When editing `.fpp`, maintain readability, prefer simple macros over deeply nested constructs.

### 3.3 Documentation

- Add or update Doxygen comments in source files.
- Update Markdown docs under `docs/` if user‑facing behavior changes.
- Provide a minimal example in `examples/` for each new feature when practical.

### 3.4 Testing

- Add regression tests that fail before your change and pass after.
- Use `./mfc.sh test --generate` to create golden files for new cases.
- Keep tests fast; favor small grids and short runtimes.

### 3.5 GPU & Performance

- Ensure code compiles for CPU *and* GPU targets (NVHPC for NVIDIA, Cray for AMD).
- Profile critical kernels; avoid introducing bottlenecks.

---

## 4. Preparing Your Pull Request

1. **One PR = One logical change**. If you plan a follow‑up change, open an issue describing it and assign yourself for visibility.
2. **Fill out the PR template**. Remove checkboxes that do **not** apply.
3. **Describe testing** – list commands, compilers, and any profiling.
4. **Link issues** – `Fixes #<id>` or `Part of #<id>`.
5. **Ensure CI passes** before requesting review.

> **Tip** If your change is large, consider splitting it into smaller PRs. Document the intent in an issue so reviewers understand the overall roadmap.

---

## 5. Code Review & Merge

- Respond promptly to reviewer comments.
- Push focused updates; each push re‑runs CI.
- When all reviews are approved and CI is green, a maintainer will merge your PR.

---

## 6. Issue Triage

If you prefer helping with issue management:
- Comment to clarify reproduction steps.
- Label issues when you have triage rights.
- Close fixed issues and reference the PR.

