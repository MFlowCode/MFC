#!/usr/bin/env python3
"""Say whether a scaling run is a regression.

    scaling-check.py <new line> [scaling.jsonl]

The new line is scaling-record.py's output; scaling.jsonl is the history before it. Exits 1 and
prints one line per problem when, against the median of the last four runs on the same lane:

    efficiency at some size has fallen by more than 5 points, or
    grind time at the smallest size has risen by more than 5%,

or, with or without history, weak-scaling efficiency at some size is under 80%.
docs/dashboard/ci-pages.js (scalingStatus) applies the same rules to the pages.
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


def problems(new: dict, history: list) -> list:
    now = measures(new)
    past = [measures(r) for r in history if (r.get("cluster"), r.get("interface")) == (new.get("cluster"), new.get("interface"))][-RUNS:]
    out = []
    for key, value in now.items():
        kind, size = key
        usual = [m[key] for m in past if key in m]
        usual = statistics.median(usual) if usual else None
        if size == "grind":
            if usual and (value / usual - 1) * 100 > GRIND:
                out.append(f"{kind} scaling: grind time at the smallest size is {value:.3g} ns, {(value / usual - 1) * 100:.1f}% above the usual {usual:.3g} ns")
            continue
        if usual is not None and usual - value > POINTS:
            out.append(f"{kind} scaling at {size} nodes: efficiency {value:.0f}%, down from the usual {usual:.0f}%")
        elif kind == "weak" and value < WEAK_FLOOR:
            out.append(f"weak scaling at {size} nodes: efficiency {value:.0f}%, under {WEAK_FLOOR:.0f}%")
    return out


if __name__ == "__main__":
    if len(sys.argv) not in (2, 3):
        sys.exit(__doc__)
    with open(sys.argv[1]) as f:
        new = json.loads(f.read())
    history = []
    if len(sys.argv) == 3:
        with open(sys.argv[2]) as f:
            history = [json.loads(line) for line in f if line.strip()]
    found = problems(new, [r for r in history if r.get("sha") != new.get("sha") or r.get("date") != new.get("date")])
    print("\n".join(f"- {p}" for p in found))
    sys.exit(1 if found else 0)
