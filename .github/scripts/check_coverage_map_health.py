"""Fail loudly if the committed coverage map is stale or under-covers. Used by coverage-health.yml."""
import datetime
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "toolchain"))
from mfc.test.coverage import COVERAGE_MAP_PATH, load_map, map_health  # noqa: E402
from mfc.test.cases import list_cases  # noqa: E402  (returns the current test list)

MAX_AGE_DAYS = 10
MIN_FRACTION = 0.80

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
        unloadable.append((getattr(b, "trace", repr(b)), str(exc).strip().splitlines()[-1][:140]))
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
)
print(msg)
sys.exit(0 if ok else 1)
