# Coverage-Based Test Selection Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Re-introduce sound, execution-based test selection (`./mfc.sh test --only-changes`) whose freshness is git-visible and whose failures are loud, replacing the removed gcov coverage-cache.

**Architecture:** A committed coverage map (`param_hash → covered .fpp files`) is built on master by a refresh workflow and consumed by a conservative selector on PRs. The selector can only over-include, never under-include. A scheduled health workflow fails loudly if the map goes stale. The selector ships in "shadow mode" first (prints what it would select; runs full suite) and is only flipped to gate PRs after shadow data proves it never under-selects.

**Tech Stack:** Python 3.9+ (toolchain), gcov/gfortran (collection), GitHub Actions (refresh/health/CI), `dorny/paths-filter` (changed-file detection).

**Spec:** `docs/superpowers/specs/2026-05-29-coverage-test-selection-design.md`

**Recovering the old collector:** the deleted gcov-collection code is at `git show 574e53d4~1:toolchain/mfc/test/coverage.py`. Tasks 8–9 reuse its mechanics verbatim except for the re-keying.

---

## File Structure

- Create `toolchain/mfc/test/coverage.py` — PR-path logic: `param_hash`, map load/save, `ALWAYS_RUN_ALL`, `get_changed_files`, `select_tests` (the ladder), `format_summary`. Pure/unit-tested.
- Create `toolchain/mfc/test/coverage_build.py` — master-path gcov collection: `build_coverage_map`. Integration; recovered from git + re-keyed.
- Create `toolchain/mfc/test/test_coverage_unit.py` — unit tests for `coverage.py`.
- Create `.github/scripts/check_coverage_map_health.py` — health check used by the health workflow. Unit-tested.
- Modify `toolchain/mfc/cli/commands.py` — re-add `test` flags: `--only-changes`, `--build-coverage-map`, `--changed-files`, `--changes-branch`.
- Modify `toolchain/mfc/test/test.py` — wire selector (shadow) + map builder dispatch.
- Create `.github/workflows/coverage-refresh.yml` — rebuild + bot-commit map on master.
- Create `.github/workflows/coverage-health.yml` — scheduled loud staleness monitor.
- Modify `.github/workflows/test.yml`, `.github/file-filter.yml`, `.github/workflows/common/test.sh` — pass changed files; enable selection (final task only).

---

## Task 1: `param_hash` — stable, trace-independent test key

**Files:**
- Create: `toolchain/mfc/test/coverage.py`
- Test: `toolchain/mfc/test/test_coverage_unit.py`

- [ ] **Step 1: Write the failing test**

```python
# toolchain/mfc/test/test_coverage_unit.py
from mfc.test.coverage import param_hash


def test_param_hash_is_order_independent():
    a = param_hash({"m": 100, "weno_order": 5, "bubbles_euler": "T"})
    b = param_hash({"weno_order": 5, "bubbles_euler": "T", "m": 100})
    assert a == b


def test_param_hash_changes_with_value():
    a = param_hash({"weno_order": 5})
    b = param_hash({"weno_order": 3})
    assert a != b


def test_param_hash_is_hex_and_short():
    h = param_hash({"m": 1})
    assert len(h) == 16 and all(c in "0123456789abcdef" for c in h)


def test_param_hash_nested_order_independent():
    a = param_hash({"patch": {"x": 1, "y": 2}})
    b = param_hash({"patch": {"y": 2, "x": 1}})
    assert a == b
```

- [ ] **Step 2: Run test to verify it fails**

Run: `build/venv/bin/python3 -m pytest toolchain/mfc/test/test_coverage_unit.py -q`
Expected: FAIL — `ModuleNotFoundError`/`ImportError: cannot import name 'param_hash'`.

- [ ] **Step 3: Write minimal implementation**

```python
# toolchain/mfc/test/coverage.py
"""Execution-coverage-based test selection (PR path).

Selection is sound (only over-includes) and its failures are loud. See
docs/superpowers/specs/2026-05-29-coverage-test-selection-design.md.
"""

import hashlib
import json


def param_hash(params: dict) -> str:
    """Stable 16-hex key identifying a test by its defining params.

    Independent of dict ordering and of the human-readable trace, so cosmetic
    cases.py edits don't change the key; a real param change does.
    """
    canonical = json.dumps(params, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()[:16]
```

- [ ] **Step 4: Run test to verify it passes**

Run: `build/venv/bin/python3 -m pytest toolchain/mfc/test/test_coverage_unit.py -q`
Expected: PASS (4 passed).

- [ ] **Step 5: Commit**

```bash
git add toolchain/mfc/test/coverage.py toolchain/mfc/test/test_coverage_unit.py
git commit -m "feat(test): add stable param_hash key for coverage selection"
```

---

## Task 2: Map load/save with freshness metadata

**Files:**
- Modify: `toolchain/mfc/test/coverage.py`
- Test: `toolchain/mfc/test/test_coverage_unit.py`

- [ ] **Step 1: Write the failing test**

```python
# append to test_coverage_unit.py
import gzip, json, tempfile
from pathlib import Path
from mfc.test.coverage import save_map, load_map


def test_save_then_load_roundtrip():
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "m.json.gz"
        save_map(p, {"abc": ["src/simulation/m_rhs.fpp"]}, n_tests=1, git_sha="deadbee", gfortran_version="13")
        entries, meta = load_map(p)
        assert entries == {"abc": ["src/simulation/m_rhs.fpp"]}
        assert meta["n_tests"] == 1 and meta["git_sha"] == "deadbee"
        assert "built_at" in meta


def test_load_missing_returns_none():
    assert load_map(Path("/nonexistent/m.json.gz")) == (None, None)


def test_load_corrupt_returns_none():
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "m.json.gz"
        p.write_bytes(b"not gzip")
        assert load_map(p) == (None, None)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `build/venv/bin/python3 -m pytest toolchain/mfc/test/test_coverage_unit.py -q -k save or load`
Expected: FAIL — `ImportError: cannot import name 'save_map'`.

- [ ] **Step 3: Write minimal implementation**

```python
# add to coverage.py
import datetime
import gzip
from pathlib import Path
from typing import Optional, Tuple

COVERAGE_MAP_PATH = Path("tests/coverage_map.json.gz")


def save_map(path: Path, entries: dict, *, n_tests: int, git_sha: str, gfortran_version: str) -> None:
    payload = dict(entries)
    payload["_meta"] = {
        "built_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "git_sha": git_sha,
        "gfortran_version": gfortran_version,
        "n_tests": n_tests,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, sort_keys=True)


def load_map(path: Path) -> Tuple[Optional[dict], Optional[dict]]:
    """Return (entries_without_meta, meta), or (None, None) if missing/corrupt."""
    if not Path(path).exists():
        return None, None
    try:
        with gzip.open(path, "rt", encoding="utf-8") as f:
            data = json.load(f)
    except (OSError, gzip.BadGzipFile, json.JSONDecodeError, UnicodeDecodeError):
        return None, None
    if not isinstance(data, dict) or "_meta" not in data:
        return None, None
    meta = data.pop("_meta")
    return data, meta
```

- [ ] **Step 4: Run test to verify it passes**

Run: `build/venv/bin/python3 -m pytest toolchain/mfc/test/test_coverage_unit.py -q`
Expected: PASS (7 passed).

- [ ] **Step 5: Commit**

```bash
git add toolchain/mfc/test/coverage.py toolchain/mfc/test/test_coverage_unit.py
git commit -m "feat(test): coverage map load/save with freshness metadata"
```

---

## Task 3: `ALWAYS_RUN_ALL` classification

**Files:**
- Modify: `toolchain/mfc/test/coverage.py`
- Test: `toolchain/mfc/test/test_coverage_unit.py`

- [ ] **Step 1: Write the failing test**

```python
# append to test_coverage_unit.py
from mfc.test.coverage import is_always_run_all


def test_macro_file_forces_all():
    assert is_always_run_all({"src/common/include/parallel_macros.fpp"})


def test_cmake_forces_all():
    assert is_always_run_all({"CMakeLists.txt"})
    assert is_always_run_all({"toolchain/cmake/foo.cmake"})


def test_param_codegen_forces_all():
    assert is_always_run_all({"toolchain/mfc/params/definitions.py"})


def test_ordinary_common_module_does_not_force_all():
    assert not is_always_run_all({"src/common/m_helper.fpp"})


def test_ordinary_sim_module_does_not_force_all():
    assert not is_always_run_all({"src/simulation/m_rhs.fpp"})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `build/venv/bin/python3 -m pytest toolchain/mfc/test/test_coverage_unit.py -q -k always_run or force_all`
Expected: FAIL — `ImportError: cannot import name 'is_always_run_all'`.

- [ ] **Step 3: Write minimal implementation**

```python
# add to coverage.py
ALWAYS_RUN_ALL_EXACT = frozenset(["CMakeLists.txt"])
ALWAYS_RUN_ALL_PREFIXES = (
    "src/common/include/",     # GPU/Fypp macro & include files (CPU map can't line-attribute)
    "toolchain/cmake/",        # build system
    "toolchain/mfc/params/",   # parameter codegen → emits Fortran broadly
    "toolchain/bootstrap/",    # build/run scripts
)


def is_always_run_all(changed_files: set) -> bool:
    """True if any changed file forces the full suite (un-attributable by the CPU map)."""
    if changed_files & ALWAYS_RUN_ALL_EXACT:
        return True
    return any(f.startswith(ALWAYS_RUN_ALL_PREFIXES) for f in changed_files)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `build/venv/bin/python3 -m pytest toolchain/mfc/test/test_coverage_unit.py -q`
Expected: PASS (12 passed).

- [ ] **Step 5: Commit**

```bash
git add toolchain/mfc/test/coverage.py toolchain/mfc/test/test_coverage_unit.py
git commit -m "feat(test): ALWAYS_RUN_ALL classification for coverage selection"
```

---

## Task 4: The conservative ladder — `select_tests`

**Files:**
- Modify: `toolchain/mfc/test/coverage.py`
- Test: `toolchain/mfc/test/test_coverage_unit.py`

- [ ] **Step 1: Write the failing test**

```python
# append to test_coverage_unit.py
from mfc.test.coverage import select_tests


class _Case:
    def __init__(self, ph, params=None):
        self._ph = ph
        self.params = params or {}

    def coverage_key(self):
        return self._ph


def _cases(*phs):
    return [_Case(p) for p in phs]


def test_rung1_no_changed_files_runs_all():
    cases = _cases("a", "b")
    run, skip, reason = select_tests(cases, {"a": ["src/x.fpp"]}, None)
    assert len(run) == 2 and skip == [] and reason.startswith("rung1")


def test_rung2_always_run_all():
    cases = _cases("a", "b")
    run, skip, reason = select_tests(cases, {"a": [], "b": []}, {"CMakeLists.txt"})
    assert len(run) == 2 and reason.startswith("rung2")


def test_rung3_f90_change_runs_all():
    cases = _cases("a")
    run, skip, reason = select_tests(cases, {"a": []}, {"src/common/m_precision_select.f90"})
    assert len(run) == 1 and reason.startswith("rung3")


def test_rung4_changed_fpp_with_zero_coverage_runs_all():
    cases = _cases("a")
    # m_gpu_only.fpp is covered by no test in the map
    run, skip, reason = select_tests(cases, {"a": ["src/simulation/m_rhs.fpp"]}, {"src/simulation/m_gpu_only.fpp"})
    assert len(run) == 1 and reason.startswith("rung4")


def test_rung5_unmapped_test_is_included():
    cases = _cases("a", "new")  # 'new' not in map
    run, skip, _ = select_tests(cases, {"a": ["src/simulation/m_rhs.fpp"]}, {"src/simulation/m_rhs.fpp"})
    assert {c.coverage_key() for c in run} == {"a", "new"}


def test_rung6_and_7_overlap_selects_subset():
    cases = _cases("hit", "miss")
    cov = {"hit": ["src/simulation/m_bubbles_EE.fpp"], "miss": ["src/simulation/m_rhs.fpp"]}
    run, skip, _ = select_tests(cases, cov, {"src/simulation/m_bubbles_EE.fpp"})
    assert [c.coverage_key() for c in run] == ["hit"]
    assert [c.coverage_key() for c in skip] == ["miss"]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `build/venv/bin/python3 -m pytest toolchain/mfc/test/test_coverage_unit.py -q -k select or rung`
Expected: FAIL — `ImportError: cannot import name 'select_tests'`.

- [ ] **Step 3: Write minimal implementation**

```python
# add to coverage.py
def _covered_fpp(coverage_map: dict) -> set:
    files = set()
    for cov in coverage_map.values():
        files.update(cov)
    return files


def select_tests(cases, coverage_map, changed_files):
    """Return (to_run, skipped, reason). Conservative ladder — only over-includes.

    `cases` items expose `.coverage_key()`. `changed_files` is a set of repo-relative
    paths, or None if detection failed.
    """
    # Rung 1: no changed-file info -> run all.
    if changed_files is None:
        return list(cases), [], "rung1: changed-file list unavailable"

    # Rung 2: macro/codegen/build inputs -> run all.
    if is_always_run_all(changed_files):
        return list(cases), [], "rung2: macro/codegen/build input changed"

    # Rung 3: changed .f90/.f under src/ (map tracks .fpp only) -> run all.
    if any(f.startswith("src/") and f.endswith((".f90", ".f")) for f in changed_files):
        return list(cases), [], "rung3: hand-written .f90/.f changed"

    changed_fpp = {f for f in changed_files if f.endswith(".fpp")}
    if not changed_fpp:
        return [], list(cases), "rung7: no Fortran source changed"

    # Rung 4: a changed .fpp that no test covers -> run all (GPU-only blind spot).
    covered = _covered_fpp(coverage_map)
    if changed_fpp - covered:
        return list(cases), [], "rung4: changed .fpp not covered by any test"

    # Rungs 5-7: per-test.
    to_run, skipped = [], []
    for case in cases:
        key = case.coverage_key()
        cov = coverage_map.get(key)
        if cov is None:  # rung 5: unmapped/new test
            to_run.append(case)
        elif set(cov) & changed_fpp:  # rung 6: overlap
            to_run.append(case)
        else:  # rung 7: skip
            skipped.append(case)
    return to_run, skipped, f"selected {len(to_run)}/{len(cases)} by coverage overlap"
```

- [ ] **Step 4: Run test to verify it passes**

Run: `build/venv/bin/python3 -m pytest toolchain/mfc/test/test_coverage_unit.py -q`
Expected: PASS (18 passed).

- [ ] **Step 5: Commit**

```bash
git add toolchain/mfc/test/coverage.py toolchain/mfc/test/test_coverage_unit.py
git commit -m "feat(test): conservative coverage selection ladder"
```

---

## Task 5: `coverage_key()` on the test case object

**Files:**
- Modify: `toolchain/mfc/test/case.py` (the `TestCase` class — find it via `grep -n "def get_uuid" toolchain/mfc/test/case.py`)
- Test: `toolchain/mfc/test/test_coverage_unit.py`

- [ ] **Step 1: Write the failing test**

```python
# append to test_coverage_unit.py
def test_case_coverage_key_matches_param_hash():
    from mfc.test.case import TestCase
    from mfc.test.coverage import param_hash
    tc = TestCase("1D -> Foo", {"m": 100, "weno_order": 5})
    assert tc.coverage_key() == param_hash({"m": 100, "weno_order": 5})


def test_case_coverage_key_ignores_trace():
    from mfc.test.case import TestCase
    a = TestCase("1D -> Foo", {"m": 100})
    b = TestCase("totally -> different -> trace", {"m": 100})
    assert a.coverage_key() == b.coverage_key()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `build/venv/bin/python3 -m pytest toolchain/mfc/test/test_coverage_unit.py -q -k coverage_key`
Expected: FAIL — `AttributeError: 'TestCase' object has no attribute 'coverage_key'`.

- [ ] **Step 3: Write minimal implementation**

Add this method to the `TestCase` class in `toolchain/mfc/test/case.py` (next to `get_uuid`):

```python
    def coverage_key(self) -> str:
        from .coverage import param_hash
        return param_hash(self.params)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `build/venv/bin/python3 -m pytest toolchain/mfc/test/test_coverage_unit.py -q -k coverage_key`
Expected: PASS (2 passed).

- [ ] **Step 5: Commit**

```bash
git add toolchain/mfc/test/case.py toolchain/mfc/test/test_coverage_unit.py
git commit -m "feat(test): TestCase.coverage_key() = param_hash(params)"
```

---

## Task 6: Changed-file resolution (CI list + self-healing git)

**Files:**
- Modify: `toolchain/mfc/test/coverage.py`
- Test: `toolchain/mfc/test/test_coverage_unit.py`

- [ ] **Step 1: Write the failing test**

```python
# append to test_coverage_unit.py
from unittest.mock import patch
import types as _types
from mfc.test.coverage import get_changed_files


def test_changed_files_prefers_explicit_list():
    files = get_changed_files("/repo", "master", explicit="src/a.fpp\nsrc/b.fpp\n")
    assert files == {"src/a.fpp", "src/b.fpp"}


def test_changed_files_deepens_then_recovers():
    state = {"deepened": False}

    def fake_run(cmd, **kw):
        sub = cmd[1] if len(cmd) > 1 else ""
        if sub == "fetch":
            state["deepened"] = True
            return _types.SimpleNamespace(returncode=0, stdout="", stderr="")
        if sub == "merge-base":
            return _types.SimpleNamespace(returncode=0 if state["deepened"] else 1, stdout="base\n", stderr="")
        if sub == "diff":
            return _types.SimpleNamespace(returncode=0, stdout="src/x.fpp\n", stderr="")
        return _types.SimpleNamespace(returncode=0, stdout="", stderr="")

    with patch("subprocess.run", fake_run):
        assert get_changed_files("/repo", "master") == {"src/x.fpp"}


def test_changed_files_returns_none_when_unrecoverable():
    def fake_run(cmd, **kw):
        rc = 1 if (len(cmd) > 1 and cmd[1] == "merge-base") else 0
        return _types.SimpleNamespace(returncode=rc, stdout="", stderr="boom")

    with patch("subprocess.run", fake_run):
        assert get_changed_files("/repo", "master") is None
```

- [ ] **Step 2: Run test to verify it fails**

Run: `build/venv/bin/python3 -m pytest toolchain/mfc/test/test_coverage_unit.py -q -k changed_files`
Expected: FAIL — `ImportError: cannot import name 'get_changed_files'`.

- [ ] **Step 3: Write minimal implementation**

```python
# add to coverage.py
import subprocess
from typing import Optional


def _git(args, cwd, timeout=60):
    return subprocess.run(["git", *args], capture_output=True, text=True, cwd=cwd, timeout=timeout, check=False)


def _merge_base(cwd, branch):
    for ref in (branch, f"origin/{branch}"):
        r = _git(["merge-base", ref, "HEAD"], cwd)
        if r.returncode == 0 and r.stdout.strip():
            return r.stdout.strip()
    return None


def get_changed_files(root_dir, compare_branch="master", explicit: Optional[str] = None):
    """Set of changed repo-relative paths, or None if undeterminable (-> run all).

    `explicit` is a newline-separated list from CI (paths-filter); preferred when given.
    Otherwise use git merge-base, self-healing a shallow clone with a deepen+retry.
    """
    if explicit is not None:
        return {f for f in explicit.splitlines() if f.strip()}
    try:
        base = _merge_base(root_dir, compare_branch)
        if base is None:
            _git(["fetch", "origin", f"{compare_branch}:{compare_branch}", "--depth=1"], root_dir, 120)
            _git(["fetch", "--deepen=200"], root_dir, 120)
            base = _merge_base(root_dir, compare_branch)
        if base is None:
            return None
        diff = _git(["diff", base, "HEAD", "--name-only", "--no-color"], root_dir)
        if diff.returncode != 0:
            return None
        return {f for f in diff.stdout.splitlines() if f.strip()}
    except (subprocess.TimeoutExpired, OSError):
        return None
```

- [ ] **Step 4: Run test to verify it passes**

Run: `build/venv/bin/python3 -m pytest toolchain/mfc/test/test_coverage_unit.py -q`
Expected: PASS (all).

- [ ] **Step 5: Commit**

```bash
git add toolchain/mfc/test/coverage.py toolchain/mfc/test/test_coverage_unit.py
git commit -m "feat(test): robust changed-file detection (CI list + self-healing git)"
```

---

## Task 7: Selection summary line (observability)

**Files:**
- Modify: `toolchain/mfc/test/coverage.py`
- Test: `toolchain/mfc/test/test_coverage_unit.py`

- [ ] **Step 1: Write the failing test**

```python
# append to test_coverage_unit.py
from mfc.test.coverage import format_summary


def test_summary_mentions_counts_age_reason():
    s = format_summary(ran=47, total=610, reason="selected 47/610 by coverage overlap",
                       meta={"built_at": "2026-05-20T00:00:00+00:00"}, now="2026-05-29T00:00:00+00:00")
    assert "47/610" in s and "9d" in s and "coverage overlap" in s


def test_summary_handles_missing_meta():
    s = format_summary(ran=610, total=610, reason="rung1: changed-file list unavailable", meta=None, now="2026-05-29T00:00:00+00:00")
    assert "610/610" in s and "map age unknown" in s
```

- [ ] **Step 2: Run test to verify it fails**

Run: `build/venv/bin/python3 -m pytest toolchain/mfc/test/test_coverage_unit.py -q -k summary`
Expected: FAIL — `ImportError: cannot import name 'format_summary'`.

- [ ] **Step 3: Write minimal implementation**

```python
# add to coverage.py
def format_summary(*, ran, total, reason, meta, now) -> str:
    if meta and meta.get("built_at"):
        built = datetime.datetime.fromisoformat(meta["built_at"])
        age_days = (datetime.datetime.fromisoformat(now) - built).days
        age = f"map age {age_days}d"
    else:
        age = "map age unknown"
    return f"Coverage selection: ran {ran}/{total} tests · {age} · {reason}"
```

- [ ] **Step 4: Run test to verify it passes**

Run: `build/venv/bin/python3 -m pytest toolchain/mfc/test/test_coverage_unit.py -q`
Expected: PASS (all).

- [ ] **Step 5: Commit**

```bash
git add toolchain/mfc/test/coverage.py toolchain/mfc/test/test_coverage_unit.py
git commit -m "feat(test): coverage selection summary line"
```

---

## Task 8: Recover the gcov collector, re-keyed by param_hash

**Files:**
- Create: `toolchain/mfc/test/coverage_build.py`

This task reuses the deleted collector's gcov mechanics. It is integration-level (needs a `--gcov` build) and is validated in Task 11, not unit-tested.

- [ ] **Step 1: Recover the old collector's mechanics**

Run: `git show 574e53d4~1:toolchain/mfc/test/coverage.py > /tmp/old_coverage.py`
This file contains `find_gcov_binary`, `find_gcno_files`, `_parse_gcov_json_output`, `_compute_gcov_prefix_strip`, `_collect_single_test_coverage`, `_run_single_test_direct`, `_prepare_test`, and `build_coverage_cache`. Copy those functions into `toolchain/mfc/test/coverage_build.py` unchanged EXCEPT `build_coverage_cache` (next step).

- [ ] **Step 2: Re-key and rename the entry point**

Replace the old `build_coverage_cache`'s cache-writing tail. Where it previously did `cache[uuid] = coverage` keyed by `case.get_uuid()`, key by `param_hash` instead, and write via `save_map`. Rename the function to `build_coverage_map`:

```python
# in coverage_build.py — the collection loop stores by param_hash:
#   key = case.coverage_key()
#   cache[key] = sorted(coverage)
# and the tail becomes:
from .coverage import COVERAGE_MAP_PATH, save_map

# ... after collection, replacing the old gzip-write block:
    n_tests = sum(1 for v in cache.values() if v is not None)
    if n_tests == 0:
        raise MFCException("Coverage build produced zero coverage. Check the --gcov build and gcov binary.")
    import subprocess
    git_sha = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, cwd=root_dir).stdout.strip()
    save_map(COVERAGE_MAP_PATH, cache, n_tests=n_tests, git_sha=git_sha, gfortran_version=_get_gcov_version(gcov_bin))
    cons.print(f"[bold green]Coverage map written to {COVERAGE_MAP_PATH}[/bold green] ({n_tests} tests)")
```

- [ ] **Step 3: Verify it imports**

Run: `build/venv/bin/python3 -c "import sys; sys.path.insert(0,'toolchain'); from mfc.test import coverage_build"`
Expected: no output, exit 0 (imports cleanly).

- [ ] **Step 4: Commit**

```bash
git add toolchain/mfc/test/coverage_build.py
git commit -m "feat(test): gcov coverage-map collector, keyed by param_hash"
```

---

## Task 9: Wire `--only-changes` (shadow) and `--build-coverage-map` into the CLI

**Files:**
- Modify: `toolchain/mfc/cli/commands.py` (the `TEST_COMMAND` arguments list — find via `grep -n "name=\"shard\"" toolchain/mfc/cli/commands.py`)
- Modify: `toolchain/mfc/test/test.py`

- [ ] **Step 1: Re-add the flags**

In `commands.py`, add to the `TEST_COMMAND` `arguments=[...]` list (after the `shard` Argument):

```python
        Argument(name="build-coverage-map", dest="build_coverage_map", action=ArgAction.STORE_TRUE, default=False,
                 help="Build the gcov coverage map (requires a prior --gcov build). Master-side only."),
        Argument(name="only-changes", dest="only_changes", action=ArgAction.STORE_TRUE, default=False,
                 help="Select only tests whose covered files overlap changed files (shadow mode unless --select-enforce)."),
        Argument(name="select-enforce", dest="select_enforce", action=ArgAction.STORE_TRUE, default=False,
                 help="With --only-changes, actually skip unselected tests (otherwise shadow: print selection, run all)."),
        Argument(name="changed-files", dest="changed_files", type=str, default=None,
                 help="Newline- or comma-separated changed-file list (from CI paths-filter). Overrides git detection."),
        Argument(name="changes-branch", dest="changes_branch", type=str, default="master",
                 help="Branch to diff against for --only-changes."),
```

- [ ] **Step 2: Add the build-map dispatch in test.py**

In `toolchain/mfc/test/test.py`, in `test()`, before the normal `cases, skipped_cases = __filter(cases)` line, add:

```python
    if ARG("build_coverage_map"):
        from .coverage_build import build_coverage_map
        import itertools
        from ..build import PRE_PROCESS, SIMULATION, POST_PROCESS, build
        all_cases = [b.to_case() for b in cases]
        unique = set()
        for case, code in itertools.product(all_cases, [PRE_PROCESS, SIMULATION, POST_PROCESS]):
            slug = code.get_slug(case.to_input_file())
            if slug not in unique:
                build(code, case.to_input_file())
                unique.add(slug)
        build_coverage_map(common.MFC_ROOT_DIR, all_cases, n_jobs=int(ARG("jobs")))
        return
```

- [ ] **Step 3: Add the selector (shadow-capable) in `__filter`**

In `__filter` in `test.py`, after the `--only` handling and before returning, add:

```python
    if ARG("only_changes"):
        import datetime
        from .. import common
        from .coverage import COVERAGE_MAP_PATH, load_map, get_changed_files, select_tests, format_summary
        entries, meta = load_map(COVERAGE_MAP_PATH)
        if entries is None:
            cons.print("[yellow]Coverage selection: map missing/corrupt — running full suite.[/yellow]")
        else:
            changed = get_changed_files(common.MFC_ROOT_DIR, ARG("changes_branch"), explicit=ARG("changed_files"))
            to_run, to_skip, reason = select_tests(cases, entries, changed)
            cons.print(format_summary(ran=len(to_run), total=len(cases), reason=reason, meta=meta,
                                      now=datetime.datetime.now(datetime.timezone.utc).isoformat()))
            if ARG("select_enforce"):
                skipped_cases += to_skip
                cases = to_run
            else:
                cons.print("[dim](shadow mode: running full suite; pass --select-enforce to actually skip)[/dim]")
```

- [ ] **Step 4: Verify the CLI parses and shadow mode runs**

Run: `./mfc.sh test --only-changes --no-build -l 2>&1 | grep -i "coverage selection\|shadow" | head`
Expected: prints a "Coverage selection: ..." line (or the map-missing message if no map yet) and the "(shadow mode ...)" note. The full list still prints.

- [ ] **Step 5: Run the unit + CLI test suites**

Run: `build/venv/bin/python3 -m pytest toolchain/mfc -q`
Expected: PASS (existing + new tests).

- [ ] **Step 6: Commit**

```bash
git add toolchain/mfc/cli/commands.py toolchain/mfc/test/test.py
git commit -m "feat(test): wire --only-changes (shadow) and --build-coverage-map"
```

---

## Task 10: Map-health check script

**Files:**
- Create: `.github/scripts/check_coverage_map_health.py`
- Test: `toolchain/mfc/test/test_coverage_unit.py`

- [ ] **Step 1: Write the failing test**

```python
# append to test_coverage_unit.py
from mfc.test.coverage import map_health


def test_health_ok():
    ok, msg = map_health(meta={"built_at": "2026-05-28T00:00:00+00:00", "n_tests": 600},
                         current_keys=set(str(i) for i in range(600)),
                         mapped_keys=set(str(i) for i in range(580)),
                         now="2026-05-29T00:00:00+00:00", max_age_days=10, min_fraction=0.8)
    assert ok, msg


def test_health_stale_fails():
    ok, msg = map_health(meta={"built_at": "2026-05-01T00:00:00+00:00", "n_tests": 600},
                         current_keys=set(["a"]), mapped_keys=set(["a"]),
                         now="2026-05-29T00:00:00+00:00", max_age_days=10, min_fraction=0.8)
    assert not ok and "stale" in msg.lower()


def test_health_undercoverage_fails():
    ok, msg = map_health(meta={"built_at": "2026-05-28T00:00:00+00:00", "n_tests": 10},
                         current_keys=set(str(i) for i in range(100)),
                         mapped_keys=set(str(i) for i in range(50)),
                         now="2026-05-29T00:00:00+00:00", max_age_days=10, min_fraction=0.8)
    assert not ok and "coverage" in msg.lower()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `build/venv/bin/python3 -m pytest toolchain/mfc/test/test_coverage_unit.py -q -k health`
Expected: FAIL — `ImportError: cannot import name 'map_health'`.

- [ ] **Step 3: Add `map_health` to coverage.py**

```python
# add to coverage.py
def map_health(*, meta, current_keys, mapped_keys, now, max_age_days, min_fraction):
    """Return (ok, message). Loud anti-rot check used by the health workflow."""
    if not meta or not meta.get("built_at"):
        return False, "Coverage map has no build metadata."
    age = (datetime.datetime.fromisoformat(now) - datetime.datetime.fromisoformat(meta["built_at"])).days
    if age > max_age_days:
        return False, f"Coverage map is STALE: {age}d old (max {max_age_days}d). Refresh workflow may be broken."
    if current_keys:
        frac = len(current_keys & mapped_keys) / len(current_keys)
        if frac < min_fraction:
            return False, f"Coverage map under-covers: {frac:.0%} of current tests mapped (min {min_fraction:.0%})."
    return True, f"Coverage map healthy: {age}d old."
```

- [ ] **Step 4: Write the CLI wrapper script**

```python
# .github/scripts/check_coverage_map_health.py
"""Fail loudly if the committed coverage map is stale or under-covers. Used by coverage-health.yml."""
import datetime
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "toolchain"))
from mfc.test.coverage import COVERAGE_MAP_PATH, load_map, map_health  # noqa: E402
from mfc.test.cases import generate_cases  # noqa: E402  (returns the current test list)

MAX_AGE_DAYS = 10
MIN_FRACTION = 0.80

entries, meta = load_map(COVERAGE_MAP_PATH)
if entries is None:
    sys.exit("Coverage map missing or corrupt.")
current_keys = {c.coverage_key() for c in (b.to_case() for b in generate_cases())}
ok, msg = map_health(meta=meta, current_keys=current_keys, mapped_keys=set(entries),
                     now=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                     max_age_days=MAX_AGE_DAYS, min_fraction=MIN_FRACTION)
print(msg)
sys.exit(0 if ok else 1)
```

Note: confirm the test-list generator name with `grep -nE "def generate_cases|def list_cases|return cases" toolchain/mfc/test/cases.py` and adjust the import if it differs.

- [ ] **Step 5: Run tests**

Run: `build/venv/bin/python3 -m pytest toolchain/mfc/test/test_coverage_unit.py -q -k health`
Expected: PASS (3 passed).

- [ ] **Step 6: Commit**

```bash
git add toolchain/mfc/test/coverage.py .github/scripts/check_coverage_map_health.py toolchain/mfc/test/test_coverage_unit.py
git commit -m "feat(test): coverage map health check (loud anti-rot)"
```

---

## Task 11: First real map build (manual, validates Task 8)

**Files:** none (produces `tests/coverage_map.json.gz`)

- [ ] **Step 1: Build instrumented + collect**

Run:
```bash
./mfc.sh build --gcov -j 8
./mfc.sh test --build-coverage-map --gcov -j 8
```
Expected: ends with "Coverage map written to tests/coverage_map.json.gz (N tests)".

- [ ] **Step 2: Sanity-check the artifact**

Run:
```bash
build/venv/bin/python3 -c "
from pathlib import Path; import sys; sys.path.insert(0,'toolchain')
from mfc.test.coverage import load_map, COVERAGE_MAP_PATH
e,m = load_map(COVERAGE_MAP_PATH)
print('tests:', len(e), 'meta:', m)
assert len(e) > 100 and all(isinstance(v, list) for v in e.values())
print('OK')
"
```
Expected: prints a plausible test count and `OK`.

- [ ] **Step 3: Shadow-verify selection against a known change**

Run:
```bash
./mfc.sh test --only-changes --changed-files "src/simulation/m_bubbles_EE.fpp" --no-build -l 2>&1 | grep -i "coverage selection"
```
Expected: a summary showing a *subset* selected (bubble tests), not all and not zero.

- [ ] **Step 4: Commit the map**

```bash
git add tests/coverage_map.json.gz
git commit -m "test: initial coverage map"
```

---

## Task 12: Refresh + health workflows

**Files:**
- Create: `.github/workflows/coverage-refresh.yml`
- Create: `.github/workflows/coverage-health.yml`

- [ ] **Step 1: Refresh workflow**

```yaml
# .github/workflows/coverage-refresh.yml
name: 'Coverage Map Refresh'
on:
  schedule:
    - cron: '0 6 * * 1'      # weekly floor
  push:
    branches: [master]
    paths:
      - 'toolchain/mfc/test/cases.py'
      - 'src/**/*.fpp'
  workflow_dispatch:
concurrency:
  group: coverage-refresh
  cancel-in-progress: true
jobs:
  refresh:
    if: github.repository == 'MFlowCode/MFC'
    timeout-minutes: 240
    runs-on:
      group:  phoenix
      labels: gt
    steps:
      - uses: actions/checkout@v4
        with: { clean: false }
      - name: Build + collect coverage map (SLURM)
        run: bash .github/scripts/submit-slurm-job.sh .github/workflows/common/coverage-refresh.sh cpu none phoenix
      - name: Commit refreshed map
        run: |
          if ! git diff --quiet tests/coverage_map.json.gz; then
            git config user.name  "mfc-bot"
            git config user.email "mfc-bot@users.noreply.github.com"
            git add tests/coverage_map.json.gz
            git commit -m "test: refresh coverage map [skip ci]"
            git push origin HEAD:master
          else
            echo "Coverage map unchanged."
          fi
```

Also create `.github/workflows/common/coverage-refresh.sh`:

```bash
#!/bin/bash
set -e
NJOBS="${SLURM_CPUS_ON_NODE:-24}"; [ "$NJOBS" -gt 64 ] && NJOBS=64
./mfc.sh clean
source .github/scripts/retry-build.sh
retry_build ./mfc.sh build --gcov -j 8
./mfc.sh test --build-coverage-map --gcov -j "$NJOBS"
```

- [ ] **Step 2: Health workflow**

```yaml
# .github/workflows/coverage-health.yml
name: 'Coverage Map Health'
on:
  schedule:
    - cron: '0 7 * * *'      # daily; loud if the refresh stopped working
  workflow_dispatch:
jobs:
  health:
    if: github.repository == 'MFlowCode/MFC'
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with: { python-version: '3.12' }
      - name: Initialize MFC
        run: ./mfc.sh init
      - name: Check coverage map freshness
        run: build/venv/bin/python3 .github/scripts/check_coverage_map_health.py
```

- [ ] **Step 3: Validate YAML**

Run:
```bash
build/venv/bin/python3 -c "import yaml; [yaml.safe_load(open(f)) for f in ['.github/workflows/coverage-refresh.yml','.github/workflows/coverage-health.yml']]; print('YAML OK')"
bash -n .github/workflows/common/coverage-refresh.sh && echo "sh OK"
```
Expected: `YAML OK` and `sh OK`.

- [ ] **Step 4: Commit**

```bash
git add .github/workflows/coverage-refresh.yml .github/workflows/coverage-health.yml .github/workflows/common/coverage-refresh.sh
git commit -m "ci: coverage map refresh + health workflows"
```

---

## Task 13: Pass changed files into PR test jobs (shadow on PRs)

**Files:**
- Modify: `.github/file-filter.yml`
- Modify: `.github/workflows/test.yml`
- Modify: `.github/workflows/common/test.sh`

- [ ] **Step 1: Emit the changed-file list from paths-filter**

In `.github/file-filter.yml`, add a filter that lists Fortran/source files (paths-filter exposes matched files via `*_files` outputs when `list-files: shell` is set). In `test.yml`'s `file-changes` job, set `list-files: shell` on the `dorny/paths-filter` step and add an output `changed_files: ${{ steps.changes.outputs.checkall_files }}`.

- [ ] **Step 2: Pass it to the test command (shadow)**

In `test.yml` github job and `common/test.sh`, add `--only-changes --changed-files "<list>"` to the `./mfc.sh test` invocation on `pull_request` events (NO `--select-enforce` yet — shadow). Plumb the list via an env var, e.g. `CHANGED_FILES: ${{ needs.file-changes.outputs.changed_files }}`.

- [ ] **Step 3: Validate YAML**

Run: `build/venv/bin/python3 -c "import yaml; yaml.safe_load(open('.github/workflows/test.yml')); print('OK')"`
Expected: `OK`.

- [ ] **Step 4: Commit**

```bash
git add .github/file-filter.yml .github/workflows/test.yml .github/workflows/common/test.sh
git commit -m "ci: run coverage selection in shadow mode on PRs"
```

- [ ] **Step 5: Open the PR (shadow). Do NOT enforce yet.**

Push the branch and open a PR. Let several PRs run with shadow output. Confirm in the logs that the printed "would-select" subsets never *exclude* a test that the full run shows is actually affected. Only after that evidence accrues, do Task 14.

---

## Task 14: Flip selection on as the PR gate (separate PR, after shadow validation)

**Files:**
- Modify: `.github/workflows/test.yml`, `.github/workflows/common/test.sh`

- [ ] **Step 1: Add `--select-enforce` to PR test invocations**

Add `--select-enforce` alongside the existing `--only-changes --changed-files ...` on `pull_request` events only. Master pushes keep running the full suite (the backstop).

- [ ] **Step 2: Validate + commit + PR**

```bash
build/venv/bin/python3 -c "import yaml; yaml.safe_load(open('.github/workflows/test.yml')); print('OK')"
git add .github/workflows/test.yml .github/workflows/common/test.sh
git commit -m "ci: enforce coverage-based test selection on PRs"
```

Open this as its own PR so the enforcement flip is an isolated, revertable change.

---

## Self-Review

**Spec coverage:** param_hash keying (T1, T5) ✓; map load/save + meta (T2) ✓; conservative ladder rungs 1-7 (T4) ✓; ALWAYS_RUN_ALL incl. macro/codegen, excl. ordinary common modules (T3) ✓; changed-file detection robust + CI list (T6) ✓; observability summary (T7) ✓; gcov collector re-keyed (T8) ✓; CLI wiring + shadow mode (T9) ✓; health check loud (T10) ✓; refresh-on-master + bot-commit + guard (T12) ✓; health workflow off PR path (T12) ✓; shadow-then-enforce build sequence (T13, T14) ✓; gate model (master full backstop) preserved (T13/T14 leave master pushes full) ✓.

**Placeholder scan:** Task 8 (recover-from-git) and Task 13 (paths-filter output wiring) reference existing code rather than inlining it — both give the exact recovery command / output name and the specific change; acceptable since the source is verbatim-recoverable. No TBDs.

**Type consistency:** `coverage_key()` (T5) is used by `select_tests` (T4) and the health script (T10); `param_hash` (T1) underlies both; `load_map`/`save_map`/`COVERAGE_MAP_PATH` (T2) used by T8/T9/T10; `select_tests` return shape `(to_run, skipped, reason)` consistent across T4/T9. Consistent.

**Tunables:** `MAX_AGE_DAYS=10`, `MIN_FRACTION=0.80` live in one place (`check_coverage_map_health.py`); refresh cadence in `coverage-refresh.yml`. Matches the spec defaults.
