#!/usr/bin/env python3
"""Say whether a scaling run is a regression.

    scaling-check.py [--partial] <new line> [scaling.jsonl]

The new line is scaling-record.py's output; scaling.jsonl is the history before it. Exits 1 and
prints one line per problem when, against the median of the last four accepted runs on the lane:

    efficiency at some size has fallen by more than 5 points,
    grind time at the smallest size has risen by more than 5%, or
    a size those runs measured has no result,

or, with or without history, weak-scaling efficiency at some size is under 80%. --partial is for a
run at other sizes: a size without a result is not a problem. Any other failure exits 2. A run
that regressed is not accepted, so it does not become the usual: the problem is
reported until a run is back within tolerance, or a run recorded with "accept": true (a manual
run of scaling.yml with that input) starts the baseline again from itself.
docs/dashboard/ci-pages.js (scalingStatus) applies the same rules to the pages.
"""

import json
import statistics
import sys

POINTS, GRIND, WEAK_FLOOR, RUNS = 5.0, 5.0, 80.0, 4


def measures(record: dict, bases: dict = None) -> dict:
    """(kind, nodes) -> efficiency in percent against the base size, and (kind, "grind") -> the base's grind time.

    The base is the size in `bases` (kind -> nodes), or the smallest the record has. A kind whose
    base has no result gives (kind, "no base") -> the base's nodes and nothing else: an efficiency
    against another size would not compare with the history's.
    """
    out = {}
    for kind in ("weak", "strong"):
        runs = {r["nodes"]: r for r in record["runs"] if r["kind"] == kind}
        nodes = (bases or {}).get(kind) or min(runs, default=None)
        if nodes is None:
            continue
        if nodes not in runs:
            out[(kind, "no base")] = nodes
            continue
        base = runs[nodes]
        out[(kind, "base")] = nodes
        out[(kind, "grind")] = base["grind"]
        for n, r in runs.items():
            if n != nodes:
                out[(kind, n)] = base["grind"] * base["ranks"] / (r["grind"] * r["ranks"]) * 100
    return out


def bases(measured: list) -> dict:
    """kind -> the base size of the latest of these runs."""
    return {kind: nodes for m in measured[-1:] for (kind, what), nodes in m.items() if what == "base"}


def against(now: dict, past: list, partial: bool = False) -> list:
    """What `now` has lost against the measures of the runs in `past`. `partial`: a size without a result is not a loss."""
    out = []
    for (kind, what), nodes in sorted(now.items(), key=str):
        if what == "no base" and not partial:
            out.append(f"{kind} scaling: no result at {nodes} nodes, the size the others are measured against")
    lost = {k for m in past for k in m if k[1] not in ("base", "grind", "no base")} - now.keys()
    for kind, nodes in sorted(lost, key=str):
        if not partial and (kind, "no base") not in now:
            out.append(f"{kind} scaling: no result at {nodes} nodes")
    for key, value in now.items():
        kind, size = key
        usual = [m[key] for m in past if key in m]
        if size in ("base", "no base") or not usual:
            continue
        usual = statistics.median(usual)
        if size == "grind":
            if (value / usual - 1) * 100 > GRIND:
                out.append(f"{kind} scaling: grind time at the smallest size is {value:.3g} ns, {(value / usual - 1) * 100:.1f}% above the usual {usual:.3g} ns")
        elif usual - value > POINTS:
            out.append(f"{kind} scaling at {size} nodes: efficiency {value:.0f}%, down from the usual {usual:.0f}%")
    return out


def problems(new: dict, history: list, partial: bool = False) -> list:
    accepted = []
    for r in history:
        if (r.get("cluster"), r.get("interface")) != (new.get("cluster"), new.get("interface")) or "runs" not in r:
            continue
        if r.get("accept"):
            accepted = [measures(r)]
        else:
            m = measures(r, bases(accepted))
            if not against(m, accepted[-RUNS:]):
                accepted.append(m)
    now = measures(new, bases(accepted))
    out = [] if new.get("accept") else against(now, accepted[-RUNS:], partial)
    said = {p.split(":")[0] for p in out}
    for (kind, size), value in now.items():
        if kind == "weak" and isinstance(size, int) and value < WEAK_FLOOR and f"weak scaling at {size} nodes" not in said:
            out.append(f"weak scaling at {size} nodes: efficiency {value:.0f}%, under {WEAK_FLOOR:.0f}%")
    return out


if __name__ == "__main__":
    partial = "--partial" in sys.argv
    sys.argv = [a for a in sys.argv if a != "--partial"]
    if len(sys.argv) not in (2, 3):
        sys.exit(__doc__)
    try:
        with open(sys.argv[1]) as f:
            new = json.loads(f.read())
        history = []
        if len(sys.argv) == 3:
            with open(sys.argv[2]) as f:
                history = [json.loads(line) for line in f if line.strip()]
        found = problems(new, [r for r in history if r.get("sha") != new.get("sha") or r.get("date") != new.get("date")], partial)
    except (OSError, ValueError, KeyError, TypeError, ZeroDivisionError) as e:
        print(f"scaling-check.py could not read its input: {e!r}", file=sys.stderr)
        sys.exit(2)
    print("\n".join(f"- {p}" for p in found))
    sys.exit(1 if found else 0)
