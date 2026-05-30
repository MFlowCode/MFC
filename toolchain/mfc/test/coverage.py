"""Execution-coverage-based test selection (PR path).

Selection is sound (only over-includes) and its failures are loud. See
docs/superpowers/specs/2026-05-29-coverage-test-selection-design.md.
"""

import datetime
import gzip
import hashlib
import json
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


ALWAYS_RUN_ALL_EXACT = frozenset(["CMakeLists.txt"])
ALWAYS_RUN_ALL_PREFIXES = (
    "src/common/include/",  # GPU/Fypp macro & include files (CPU map can't line-attribute)
    "toolchain/cmake/",  # build system
    "toolchain/mfc/params/",  # parameter codegen -> emits Fortran broadly
    "toolchain/bootstrap/",  # build/run scripts
)


def is_always_run_all(changed_files: set) -> bool:
    """True if any changed file forces the full suite (un-attributable by the CPU map)."""
    if changed_files & ALWAYS_RUN_ALL_EXACT:
        return True
    return any(f.startswith(ALWAYS_RUN_ALL_PREFIXES) for f in changed_files)


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
