"""Execution-coverage-based test selection (PR path).

Selection is sound (only over-includes) and its failures are loud. See
docs/superpowers/specs/2026-05-29-coverage-test-selection-design.md.
"""

import datetime
import gzip
import hashlib
import json
import subprocess
from pathlib import Path
from typing import Optional, Tuple


def param_hash(params: dict) -> str:
    """Stable 16-hex key identifying a test by its defining params.

    Independent of dict ordering and of the human-readable trace, so cosmetic
    cases.py edits don't change the key; a real param change does.
    """
    canonical = json.dumps(params, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()[:16]


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


# Test-definition file: changing it adds/modifies tests, but only the tests it touches
# (their param_hash changes -> not in map -> rung 5 runs them). NOT in ALWAYS_RUN_ALL so a
# test addition doesn't blanket-run the whole suite.
CASES_PY = "toolchain/mfc/test/cases.py"

ALWAYS_RUN_ALL_EXACT = frozenset(
    [
        "CMakeLists.txt",
    ]
)
ALWAYS_RUN_ALL_PREFIXES = (
    "src/common/include/",  # GPU/Fypp macro & include files (CPU map can't line-attribute)
    "toolchain/cmake/",  # build system
    "toolchain/bootstrap/",  # build/run scripts
)


def is_always_run_all(changed_files: set) -> bool:
    """True if any changed file forces the full suite (un-attributable by the CPU map)."""
    if changed_files & ALWAYS_RUN_ALL_EXACT:
        return True
    if any(f.startswith(ALWAYS_RUN_ALL_PREFIXES) for f in changed_files):
        return True
    # Any toolchain/mfc/*.py change (params/, run/, test infra, case.py, common.py,
    # build.py, state.py, sched.py, …) affects EVERY test's generation or execution and
    # cannot be attributed to individual tests by the coverage map.  Treat the entire
    # toolchain/mfc/ subtree as run-all EXCEPT cases.py, which is handled precisely by
    # rung 5 (new/modified tests have a fresh param_hash absent from the map and run
    # individually; unchanged tests have no .fpp overlap and are skipped).
    if any(f.startswith("toolchain/mfc/") and f.endswith(".py") and f != CASES_PY for f in changed_files):
        return True
    # gcov rolls #:include'd .fpp into the parent compilation unit, so include files
    # (inline_*.fpp, HardcodedIC, macros) are not reliably attributed in the map. Force a
    # full run for ANY src/**/include/ change so this attribution gap can never cause
    # under-inclusion — by rule, not by relying on the file being absent from the map.
    return any(f.startswith("src/") and "/include/" in f and f.endswith(".fpp") for f in changed_files)


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


def _covered_fpp(coverage_map: dict) -> set:
    files = set()
    for cov in coverage_map.values():
        files.update(cov)
    return files


def select_tests(cases, coverage_map, changed_files):
    """Return (to_run, skipped, reason). Conservative ladder -- only over-includes.

    `cases` items expose `.coverage_key()`. `changed_files` is a set of repo-relative
    paths, or None if detection failed.
    """
    # Rung 1: no changed-file info -> run all.
    if changed_files is None:
        return list(cases), [], "rung1: changed-file list unavailable"

    # Rung 2: macro/codegen/build inputs -> run all.
    if is_always_run_all(changed_files):
        return list(cases), [], "rung2: macro/codegen/build input changed"

    # Rung 3: changed hand-written Fortran source under src/ (map tracks .fpp only) -> run
    # all.  Match case-insensitively to catch .F90, .F95, .FOR, etc.
    _FORTRAN_EXTS = (".f90", ".f", ".f95", ".f03", ".f08", ".for")
    if any(f.startswith("src/") and f.lower().endswith(_FORTRAN_EXTS) for f in changed_files):
        return list(cases), [], "rung3: hand-written .f90/.f changed"

    changed_fpp = {f for f in changed_files if f.endswith(".fpp")}
    # Skip-all only when nothing test-relevant changed. If cases.py changed (no .fpp), fall
    # through to per-test: new/modified tests have a fresh param_hash absent from the map and
    # run via rung 5; unchanged tests have no .fpp overlap and are skipped.
    if not changed_fpp and CASES_PY not in changed_files:
        return [], list(cases), "rung7: no Fortran or test-definition change"

    # Rung 4: a changed .fpp that no test covers -> run all (GPU-only blind spot).
    covered = _covered_fpp(coverage_map)
    if changed_fpp - covered:
        return list(cases), [], "rung4: changed .fpp not covered by any test"

    # Rungs 5-7: per-test.
    to_run, skipped = [], []
    for case in cases:
        key = case.coverage_key()
        cov = coverage_map.get(key)
        if not cov:  # rung 5: unmapped/new test, or empty (uncertain) coverage -> run
            to_run.append(case)
        elif set(cov) & changed_fpp:  # rung 6: overlap
            to_run.append(case)
        else:  # rung 7: skip
            skipped.append(case)
    return to_run, skipped, f"selected {len(to_run)}/{len(cases)} by coverage overlap"


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

    `explicit` is a changed-file list from CI (paths-filter); preferred when given. It may
    be separated by newlines, spaces, or commas (paths-filter's shell output is space-sep).
    Otherwise use git merge-base, self-healing a shallow clone with a deepen+retry.
    """
    if explicit is not None:
        files = {f for f in explicit.replace(",", " ").split() if f.strip()}
        if files:
            return files
        # explicit was given but empty/whitespace -> ambiguous (a paths-filter/env failure vs
        # genuinely nothing). Per the soundness invariant, uncertainty must run, not skip:
        # fall through to git detection, ultimately None -> rung 1 (run all). Never return
        # an empty set here, which would be read as "nothing changed -> skip all".
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


def format_summary(*, ran, total, reason, meta, now) -> str:
    if meta and meta.get("built_at"):
        built = datetime.datetime.fromisoformat(meta["built_at"])
        age_days = (datetime.datetime.fromisoformat(now) - built).days
        age = f"map age {age_days}d"
    else:
        age = "map age unknown"
    return f"Coverage selection: ran {ran}/{total} tests · {age} · {reason}"


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
