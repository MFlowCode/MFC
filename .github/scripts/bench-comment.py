#!/usr/bin/env python3
"""Write a pull request's benchmark comment: its grind times against master's.

    bench-comment.py <artifacts dir> [benchmarks.jsonl]

The directory holds one folder per lane, bench-compare-<cluster>-<device>-<interface>, each with
pr/bench-*.yaml and master/bench-*.yaml from the same benchmark job (minutes apart; on Phoenix the same node).
benchmarks.jsonl is master's history (MFlowCode/MFC-ci-data): it says how much a lane's grind
time moves between runs with no change at all, and a difference is only called a change when it
is at least 3% and three times that. Prints Markdown, or nothing when no lane has both results.
The files come from a pull request's own run: they are read as numbers and names, nothing more.
"""

import json
import math
import os
import re
import statistics
import sys

import yaml

MARK = "<!-- mfc-bench-comment -->"
LANE = re.compile(r"^bench-compare-(frontier_amd|frontier|phoenix)-(gpu|cpu)-(acc|omp|none)$")
NAMES = {"frontier": "Frontier CCE", "frontier_amd": "Frontier AMD", "phoenix": "NVIDIA"}
OFFLOAD = {"acc": "OpenACC", "omp": "OpenMP", "none": "CPU"}
FLOOR, SCATTERS = 3.0, 3.0


def grinds(path: str) -> dict:
    """Case -> grind time, for the cases of a bench.yaml that have one."""
    try:
        with open(path) as f:
            cases = (yaml.safe_load(f) or {}).get("cases") or {}
    except (OSError, yaml.YAMLError):
        return {}
    out = {}
    for slug, case in cases.items():
        grind = (((case or {}).get("output_summary") or {}).get("simulation") or {}).get("grind")
        if isinstance(grind, (int, float)) and math.isfinite(grind) and grind > 0 and re.fullmatch(r"[\w.-]{1,80}", str(slug)):
            out[str(slug)] = grind
    return out


def scatter(history: str) -> dict:
    """(cluster, device, interface, case) -> typical run-to-run change of master's grind time, in percent.

    Only between runs alike in everything but the commit: the same GPU model (Phoenix has several, a
    few times apart in speed), the same case definition and the same modules.
    """
    series = {}
    try:
        with open(history) as f:
            for line in f:
                r = json.loads(line)
                if r.get("approx"):
                    continue
                for slug, case in r["cases"].items():
                    lane = (r["cluster"], r["device"], r["interface"], slug)
                    series.setdefault(lane, {}).setdefault((r.get("hardware"), case.get("def"), r.get("env")), []).append(case["grind"])
    except (OSError, ValueError, KeyError):
        return {}
    out = {}
    for key, alike in series.items():
        steps = [abs(b / a - 1) * 100 for values in alike.values() for a, b in zip(values, values[1:]) if a > 0][-40:]
        if len(steps) >= 10:
            out[key] = statistics.median(steps)
    return out


def comment(root: str, history: str = "") -> str:
    noise, rows, flagged = scatter(history) if history else {}, [], 0
    for folder in sorted(os.listdir(root)):
        m = LANE.match(folder)
        if not m:
            continue
        name = f"bench-{m.group(2)}-{m.group(3)}.yaml"
        pr, master = grinds(os.path.join(root, folder, "pr", name)), grinds(os.path.join(root, folder, "master", name))
        for slug in sorted(pr.keys() & master.keys()):
            pct, typical = (pr[slug] / master[slug] - 1) * 100, noise.get((*m.groups(), slug))
            real = abs(pct) >= max(FLOOR, SCATTERS * (typical or 0))
            flagged += real and pct > 0
            verdict = ("**slower**" if pct > 0 else "faster") if real else "within noise"
            rows.append(f"| {NAMES[m.group(1)]} · {OFFLOAD[m.group(3)]} | `{slug}` | {master[slug]:.3g} | {pr[slug]:.3g} | {pct:+.1f}% | {f'±{typical:.1f}%' if typical is not None else '–'} | {verdict} |")
    if not rows:
        return ""
    head = f"**{flagged} benchmark {'case is' if flagged == 1 else 'cases are'} slower than master beyond noise.**" if flagged else "No benchmark case is slower than master beyond noise."
    return "\n".join(
        [
            MARK,
            "### Benchmarks: this pull request against master",
            "",
            head,
            "",
            "| Lane | Case | master | PR | Change | Usual run-to-run | |",
            "|---|---|---|---|---|---|---|",
            *rows,
            "",
            f"Grind time, ns per grid point, equation and RK stage; lower is faster. Both sides ran in the same benchmark job, minutes apart (on Phoenix on the same node). A change counts when it is at least {FLOOR:.0f}% and {SCATTERS:.0f} times the lane's usual run-to-run change on master. [History](https://mflowcode.github.io/dashboard/benchmarks.html)",
        ]
    )


if __name__ == "__main__":
    if len(sys.argv) not in (2, 3):
        sys.exit(__doc__)
    print(comment(*sys.argv[1:]))
