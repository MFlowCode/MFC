"""Fail loudly if the committed coverage map is stale or under-covers. Used by coverage-health.yml."""
import datetime
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "toolchain"))
from mfc.test.coverage import COVERAGE_MAP_PATH, load_map, map_health  # noqa: E402
from mfc.test.cases import list_cases  # noqa: E402  (returns the current test list)

MAX_AGE_DAYS = 10
MIN_FRACTION = 0.80

# Must mirror the `paths:` trigger of coverage-refresh.yml. A superset would raise false
# alarms (a change that never triggers a refresh would look like a missed refresh).
COVERAGE_RELEVANT_PATHS = [":(glob)src/**/*.fpp", "toolchain/mfc/test/cases.py"]


def built_after_last_change(git_sha):
    """Was the map built at or after the last coverage-relevant commit? None if unknowable.

    This is the direct question the wall-clock age check only approximated: it stays quiet
    over a genuinely quiet repo (no refresh needed, so no heartbeat commit needed either)
    and fires on the next relevant push after a refresh breaks, rather than 10 days later.
    """
    if not git_sha:
        return None
    last = subprocess.run(["git", "log", "-1", "--format=%H", "--", *COVERAGE_RELEVANT_PATHS], capture_output=True, text=True, check=False)
    if last.returncode != 0 or not last.stdout.strip():
        return None  # shallow clone or no such commit -> fall back to the age rule
    ancestor = subprocess.run(["git", "merge-base", "--is-ancestor", last.stdout.strip(), git_sha], capture_output=True, check=False)
    return {0: True, 1: False}.get(ancestor.returncode)  # anything else -> None (unknown sha, shallow history)

entries, meta = load_map(COVERAGE_MAP_PATH)
if entries is None:
    sys.exit("Coverage map missing or corrupt.")
# Compute each current test's coverage key. Loading a case executes its case
# file; some (e.g. chemistry examples) import optional deps like cantera that are
# not installed in this lightweight job. Skip any case that fails to load instead
# of crashing — map_health measures the fraction of *loadable* current tests that
# are mapped, so a smaller current_keys cannot produce a false "stale" result.
current_keys = set()
unloadable = []
for b in list_cases():
    try:
        current_keys.add(b.to_case().coverage_key())
    except Exception as exc:  # noqa: BLE001 — a case file that won't import must not crash the health check
        last_line = (str(exc).strip().splitlines() or ["(no message)"])[-1][:140]
        unloadable.append((getattr(b, "trace", repr(b)), last_line))
if unloadable:
    print(f"Note: {len(unloadable)} case(s) could not be loaded in this lightweight job (excluded from the check):")
    for trace, err in unloadable[:15]:
        print(f"  - {trace}: {err}")
ok, msg = map_health(
    meta=meta,
    current_keys=current_keys,
    mapped_keys=set(entries),
    now=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    max_age_days=MAX_AGE_DAYS,
    min_fraction=MIN_FRACTION,
    built_after_last_change=built_after_last_change(meta.get("git_sha")),
)
print(msg)
sys.exit(0 if ok else 1)
