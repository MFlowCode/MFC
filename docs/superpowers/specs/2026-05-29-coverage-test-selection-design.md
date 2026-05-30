# Coverage-based test selection, done right

**Status:** Design / awaiting review
**Date:** 2026-05-29
**Supersedes:** the removed gcov coverage-cache (`--only-changes` / `--build-coverage-cache`), deleted in MFlowCode/MFC#1460.

## Problem

CI turnaround is dominated by running the full ~610-test suite on every PR, across a slow self-hosted compiler matrix. We want to run only the tests a change can actually affect — but the previous attempt was unreliable and was removed. The post-mortem (see #1460) found the *idea* (execution-based selection) was sound; the *implementation* rotted:

- keyed on `UUID = CRC32(trace)`, so any `cases.py` edit silently orphaned entries;
- the committed cache shipped ~3 months stale and nothing ever auto-refreshed it;
- every failure path silently fell back to the full suite, so breakage was invisible;
- `git merge-base` failed on shallow self-hosted clones, defeating selection exactly on the slow jobs.

We also established that **only execution-based selection is sound**: static feature/label maps miss transitive side effects (a change to a shared module silently alters behavior of a load-bearing module like `m_bubbles.fpp`), and a static `use`/`#:include` dependency graph over-approximates to "all tests" because linking ≠ execution. Therefore selection must use real per-test execution coverage (gcov), and the engineering challenge is making the freshness and reliability *trustworthy*.

## Invariants (the north stars)

Every design choice serves these two rules:

1. **Soundness — the selector may run MORE tests than necessary, never FEWER.** Under-inclusion is structurally impossible; every uncertainty resolves to "run it." Only over-inclusion (wasted time) is permitted.
2. **Anti-rot — staleness is always LOUD, never silent.** No path silently degrades to the full suite without surfacing it. A stale or under-covering map is impossible to *not* notice.

## Non-goals

- Changing the golden-file identity model (`UUID = CRC32(trace)`); coverage keying is fully decoupled from it.
- Line-level coverage. File-level (`.fpp`) is sufficient and far cheaper to collect and store.
- Selection on GPU/Cray-specific behavior — the map is CPU/gfortran-only; GPU paths are handled by the conservative ladder + the master backstop, not by selection.
- Eliminating the cost of refreshing the map: refreshing inherently runs the full suite under instrumentation. We amortize it; we do not avoid it.

## Architecture

### Components

| Component | Responsibility |
|---|---|
| **Coverage map** (`tests/coverage_map.json.gz`, committed) | The single authoritative artifact: `param_hash → covered .fpp files`. Its git history *is* the freshness record. |
| **Map builder** (`./mfc.sh test --build-coverage-map`) | Resurrected gcov collection from the deleted `coverage.py`, re-keyed by `param_hash`. |
| **Selector** (`./mfc.sh test --only-changes`) | Applies the conservative ladder to choose which tests run; emits a visible summary. |
| **Refresh workflow** (master) | Rebuilds the map on a cadence + relevant triggers and commits it back to master via a bot. |
| **Map-health workflow** (master, scheduled) | Fails LOUD (red + notify) if the map is stale or under-covers. The anti-rot enforcement, deliberately kept OFF the PR path so it can never wedge a PR. |

### Data model

The map is gzipped JSON:

```json
{
  "<param_hash>": ["src/simulation/m_rhs.fpp", "src/common/m_helper.fpp", "..."],
  "_meta": { "built_at": "<ISO-8601 UTC>", "git_sha": "<sha>", "gfortran_version": "<str>", "n_tests": 610 }
}
```

**`param_hash`** = SHA-256 over the test's *canonical* defining param dict — sorted keys, normalized float repr — as produced by `cases.py`, **excluding** coverage-run overrides (`t_step_stop=1`, etc.). Properties this gives us:

- Cosmetic trace edits (renaming/reordering a `CaseGeneratorStack` segment) leave params unchanged → **same key** (robust to the thing that rotted the old system).
- A real param change (the thing that changes behavior) → **new key** → not in map → run it (correct *and* conservative).
- Coverage keying is independent of the golden-file UUID.

### The conservative ladder (soundness, made concrete)

A test runs if **any** rung holds; rungs are checked in order, each resolving uncertainty toward "run":

1. **Changed-file list unavailable** → run **ALL** (loud annotation).
2. A changed file is in **`ALWAYS_RUN_ALL`**: the GPU/Fypp macro & include files (`src/common/include/**`), `CMakeLists.txt`, `toolchain/cmake/**`, the parameter codegen (`toolchain/mfc/params/**`), and build/run scripts → run **ALL**. (Note: ordinary `src/common/**.fpp` *modules* are NOT here — they are coverage-tracked and selectable across all three targets, which is the whole point; only the macro/codegen/build inputs that the CPU map can't attribute line-for-line are forced to run all.)
3. A changed **`.f90`/`.f`** under `src/` (the map tracks `.fpp` only) → run **ALL**.
4. A changed **`.fpp` covered by *zero* tests** in the map → run **ALL** (loud). This single rung closes the GPU blind spot: a file compiled only under `#ifdef MFC_OpenACC` is never covered by the CPU map, so a change to it forces a full run instead of silently selecting nothing.
5. A test's **`param_hash` is absent** from the map (new/changed test) → run **it**.
6. A test's covered files **∩ changed `.fpp` ≠ ∅** → run **it**.
7. Otherwise → skip.

The only residual gap is one shared by *all* coverage-based test-impact-analysis: a change that newly makes test *T* reach file *Y* when *T* did not previously execute the edited file. This is caught by the **master full-suite backstop** below.

### Gate model

- **PR required check = the selected subset** (fast).
- **Master push runs the FULL suite** — the backstop for (a) the reachability-change gap above and (b) GPU/Cray paths the CPU map cannot see.
- A rare escape is caught on master and forward-fixed. (Chosen explicitly over a full-suite merge gate, for maximum PR-time savings.)

### Changed-file detection

The selector consumes the changed-file list from CI's existing `dorny/paths-filter` step (already computed reliably) rather than `git merge-base` on a shallow runner — eliminating the old shallow-clone failure. For **local** runs it uses `git diff` against the merge-base with a self-healing `git fetch --deepen` retry; if that still fails, it runs ALL (rung 1), loudly.

## Data flows

**Refresh (on master):**
```
trigger (weekly floor | push touching cases.py or use/#:include graph)
  → build --gcov  →  test --build-coverage-map  (full suite under gcov, per-test GCOV_PREFIX isolation)
  → write tests/coverage_map.json.gz  →  bot commits to master (guarded from re-triggering the test matrix)
```

**PR selection:**
```
paths-filter → changed files
  → load committed map (loud if missing/corrupt → run ALL, annotated)
  → conservative ladder → selected subset
  → emit summary (ran N/M; map age; rung reasons)
```

## Observability

Every PR test job prints a one-line summary to the job log/summary, e.g.:
```
Coverage selection: ran 47/610 tests · map age 2d · always-run rungs: none · skipped 563 (no overlap)
```
This makes per-PR selection behavior visible at a glance. Breakage shows up as "ran 610/610 (rung 1: changed-file list unavailable)" rather than a silent full run.

## Anti-rot enforcement

The **map-health workflow** runs on a schedule on master and **fails red + notifies** if either:
- `_meta.built_at` is older than the **staleness threshold**, or
- the map covers fewer than the **minimum-coverage fraction** of current `param_hash`es (too many conservative-includes ⇒ selection is degrading).

It is intentionally NOT a PR-required check, so a refresh outage can never block unrelated PRs — PRs degrade to run-all *visibly* instead. Rot surfaces as a red master workflow + notification.

## CI integration summary

- `test.yml` PR jobs: pass the paths-filter changed-file list to `./mfc.sh test --only-changes`; selected subset is the required check.
- `test.yml` master push: full suite (unchanged behavior = backstop).
- New `coverage-refresh.yml` (master): rebuild + bot-commit the map. Guard the map-only commit from re-triggering the full matrix (`paths-ignore` on the map file and/or a `[skip ci]` convention).
- New `coverage-health.yml` (scheduled): the loud anti-rot monitor.

## Build sequence

1. **Resurrect + re-key the collector** by `param_hash`. Unit test: cosmetic trace edits keep the key; a param change alters it; collection is deterministic.
2. **Selector + observability**, landed **disabled (shadow mode)**: it computes and prints what it *would* select but the suite still runs in full. This yields real evidence that selection never under-selects, at zero risk.
3. **Map builder + first committed map**, built on master.
4. **Refresh workflow** (bot commits) + **map-health workflow** (loud monitor).
5. **Flip the selector on** as the PR gate once shadow data shows it never under-selects.

Shadow mode (step 2) is the safety valve: selection proves itself against reality before it gates anything.

## Tunable parameters (defaults; confirm during review)

- **Refresh cadence:** weekly floor (Mon 06:00 UTC) **plus** master pushes touching `cases.py` or the `use`/`#:include` dependency graph.
- **Staleness threshold:** 10 days (3 days' grace beyond the weekly refresh).
- **Minimum-coverage fraction:** 80% of current `param_hash`es present in the map.
- **Map storage:** committed at `tests/coverage_map.json.gz`, refreshed by bot push to master (matches the existing homebrew-release push pattern). Alternative if master-history churn is unwanted: a dedicated `coverage-map` branch or a release asset.

## Risks and mitigations

- **Refresh cost (full gcov run, ~hours on Phoenix).** Inherent to coverage; mitigated by amortizing (weekly + targeted triggers, never per-PR) and made visible by the health workflow.
- **CPU-only map misses GPU-only files.** Mitigated by ladder rung 4 (zero-coverage changed file → run all) + the master backstop.
- **Reachability-change gap.** Mitigated by the master full-suite backstop; inherent to all coverage TIA.
- **Bot push to master.** Uses the existing CI-push pattern; guarded against CI re-trigger loops.
