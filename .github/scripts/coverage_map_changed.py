"""Decide whether a rebuilt coverage map is worth committing. Used by coverage-refresh.yml.

The refresh job rebuilds the map on every master push touching src/**/*.fpp or cases.py,
but the rebuild is usually identical in substance -- only _meta (built_at, git_sha) moves.
A file-level `git diff` therefore always reports a change, so the bot pushed a no-op commit
on nearly every run. Compare the coverage entries instead.

Exit codes: 0 = entries changed, commit. UNCHANGED (10) = skip. Anything else = error.

`unchanged` deliberately is NOT 1: an uncaught exception exits 1, and the caller must not
be able to mistake a crashed comparison for "nothing to push" -- that would turn a broken
refresh into a silently green job. Every code outside {0, 10} fails the workflow.
"""

import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "toolchain"))
from mfc.test.coverage import COVERAGE_MAP_PATH, entries_equal, load_map  # noqa: E402

UNCHANGED = 10

new_entries, _ = load_map(COVERAGE_MAP_PATH)
if new_entries is None:
    print(f"Rebuilt {COVERAGE_MAP_PATH} is missing or corrupt.", file=sys.stderr)
    sys.exit(2)

committed = subprocess.run(["git", "show", f"HEAD:{COVERAGE_MAP_PATH}"], capture_output=True, check=False)
if committed.returncode != 0:
    print("No coverage map committed at HEAD -> commit the rebuilt one.")
    sys.exit(0)

with tempfile.NamedTemporaryFile(suffix=".json.gz") as tmp:
    tmp.write(committed.stdout)
    tmp.flush()
    old_entries, _ = load_map(Path(tmp.name))

if old_entries is None:
    print("Committed coverage map is corrupt -> commit the rebuilt one.")
    sys.exit(0)

if entries_equal(old_entries, new_entries):
    print(f"Coverage entries unchanged ({len(new_entries)} tests) -> nothing to push.")
    sys.exit(UNCHANGED)

added = len(set(new_entries) - set(old_entries))
removed = len(set(old_entries) - set(new_entries))
recovered = sum(1 for k in set(old_entries) & set(new_entries) if sorted(old_entries[k]) != sorted(new_entries[k]))
print(f"Coverage entries changed: +{added} -{removed}, {recovered} with different coverage -> committing.")
sys.exit(0)
