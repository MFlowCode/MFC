#!/usr/bin/env python3
"""Say whether a scaling run is a regression.

    scaling-check.py <new line> [scaling.jsonl]

The new line is scaling-record.py's output; scaling.jsonl is the history before it. Exits 1 and
prints one line per problem when, against the median of the last four accepted runs on the lane:

    efficiency at some size has fallen by more than 5 points,
    grind time at the smallest size has risen by more than 5%, or
    a size those runs measured has no result,

or, with or without history, weak-scaling efficiency at some size is under 80%. Any other failure
exits 2. A run that regressed is not accepted, so it does not become the usual: the problem is
reported until a run is back within tolerance, or a run recorded with "accept": true (a manual
run of scaling.yml with that input) starts the baseline again from itself.
docs/dashboard/ci-pages.js (scalingStatus) shows the same tolerances on the pages.
"""

import json
import statistics
import sys

POINTS, GRIND, WEAK_FLOOR, RUNS = 5.0, 5.0, 80.0, 4


def measures(record: dict) -> dict:
    """(kind, nodes) -> efficiency in percent, and (kind, "grind") -> grind time at the smallest size."""
    out = {}
    for kind in ("weak", "strong"):
        runs = sorted((r for r in record["runs"] if r["kind"] == kind), key=lambda r: r["nodes"])
        if not runs:
            continue
        base = runs[0]
        out[(kind, "grind")] = base["grind"]
        for r in runs[1:]:
            out[(kind, r["nodes"])] = base["grind"] * base["ranks"] / (r["grind"] * r["ranks"]) * 100
    return out


def against(now: dict, past: list) -> list:
    """What `now` has lost against the measures of the runs in `past`."""
    out = []
    for key in sorted({k for m in past for k in m} - now.keys(), key=str):
        out.append(f"{key[0]} scaling: no result {'at the smallest size' if key[1] == 'grind' else f'at {key[1]} nodes'}")
    for key, value in now.items():
        kind, size = key
        usual = [m[key] for m in past if key in m]
        if not usual:
            continue
        usual = statistics.median(usual)
        if size == "grind":
            if (value / usual - 1) * 100 > GRIND:
                out.append(f"{kind} scaling: grind time at the smallest size is {value:.3g} ns, {(value / usual - 1) * 100:.1f}% above the usual {usual:.3g} ns")
        elif usual - value > POINTS:
            out.append(f"{kind} scaling at {size} nodes: efficiency {value:.0f}%, down from the usual {usual:.0f}%")
    return out


def problems(new: dict, history: list) -> list:
    accepted = []
    for r in history:
        if (r.get("cluster"), r.get("interface")) != (new.get("cluster"), new.get("interface")) or "runs" not in r:
            continue
        m = measures(r)
        if r.get("accept"):
            accepted = [m]
        elif not against(m, accepted[-RUNS:]):
            accepted.append(m)
    now = measures(new)
    out = [] if new.get("accept") else against(now, accepted[-RUNS:])
    said = {p.split(":")[0] for p in out}
    for (kind, size), value in now.items():
        if kind == "weak" and size != "grind" and value < WEAK_FLOOR and f"weak scaling at {size} nodes" not in said:
            out.append(f"weak scaling at {size} nodes: efficiency {value:.0f}%, under {WEAK_FLOOR:.0f}%")
    return out


if __name__ == "__main__":
    if len(sys.argv) not in (2, 3):
        sys.exit(__doc__)
    try:
        with open(sys.argv[1]) as f:
            new = json.loads(f.read())
        history = []
        if len(sys.argv) == 3:
            with open(sys.argv[2]) as f:
                history = [json.loads(line) for line in f if line.strip()]
        found = problems(new, [r for r in history if r.get("sha") != new.get("sha") or r.get("date") != new.get("date")])
    except (OSError, ValueError, KeyError, TypeError, ZeroDivisionError) as e:
        print(f"scaling-check.py could not read its input: {e!r}", file=sys.stderr)
        sys.exit(2)
    print("\n".join(f"- {p}" for p in found))
    sys.exit(1 if found else 0)
