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
