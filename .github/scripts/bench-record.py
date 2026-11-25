#!/usr/bin/env python3
"""Turn one benchmark result file into one line of the benchmark history.

    bench-record.py <bench.yaml> <cluster> <device> <interface> <sha> <date> <subject> [hardware] [compare.json]

Run from the root of the benchmarked checkout. Prints a JSON object: the commit (hash, date,
subject, pull request), the lane (with the GPU model where a cluster has several), and for each case its grind time (ns per grid point, equation and RK
stage) and `def`, a hash of what the case is: its input file, its arguments and the
benchmark invocation. A grind time is only comparable with another of the same `def`, and of
the same `env`: a hash of the cluster's own lines of toolchain/modules at that commit.
With compare.json (the API's `compare/<previously benchmarked commit>...<sha>`), the line also
lists `covers`: every commit since the previous benchmark, as [sha, subject, pull request]. Any
of them can be behind a change in the numbers.
The history (one such line per lane per scheduled run) feeds docs/dashboard/benchmarks.html.
"""

import hashlib
import json
import math
import re
import sys

import yaml


def definition(description: dict, invocation: list, read=None) -> str:
    """Hash of what a case is. `read(path)` returns the input file's bytes (default: from disk)."""
    # Left out because they vary with the runner or node, not the case: the output path
    # (-o <file>), the GPU ids (-g <ids...>) and the checkout the input file sits in.
    def stable(args):
        out, skip = [], None
        for a in args or []:
            if skip == "-o" or (skip == "-g" and a.isdigit()):
                skip = None if skip == "-o" else skip
                continue
            skip = a if a in ("-o", "--output", "-g", "--gpus") else None
            if skip in ("--output", "--gpus"):
                skip = skip[1:3]
            if skip is None:
                out.append(a)
        return out

    path = description.get("path") or ""
    path = path[path.rfind("benchmarks/") :] if "benchmarks/" in path else path
    h = hashlib.sha1(json.dumps([stable(description.get("args")), stable(invocation)]).encode())
    try:
        h.update(read(path) if read else open(path, "rb").read())
    except (OSError, ValueError):
        pass
    return h.hexdigest()[:10]


PR = re.compile(r"\(#(\d+)\)\s*$")  # a squash merge's subject ends with its pull request number


def covers(compare_path: str) -> list:
    try:
        with open(compare_path) as f:
            commits = json.load(f).get("commits", [])
    except (OSError, ValueError):
        return []
    out = []
    for c in commits:
        subject = (c.get("commit", {}).get("message") or "").split("\n")[0]
        pr = PR.search(subject)
        out.append([c.get("sha"), subject, int(pr.group(1)) if pr else None])
    return out


SLUGS = {"frontier": "f", "frontier_amd": "famd", "phoenix": "p"}  # each cluster's entries in toolchain/modules


def environment(cluster: str, read=None) -> str:
    """Hash of the cluster's own lines of toolchain/modules: its compiler and module set as declared here."""
    try:
        text = (read("toolchain/modules") if read else open("toolchain/modules", "rb").read()).decode(errors="replace")
    except (OSError, ValueError):
        return ""
    mine = [line.strip() for line in text.splitlines() if re.match(rf"^{SLUGS.get(cluster, cluster)}(-|\s)", line)]
    return hashlib.sha1("\n".join(mine).encode()).hexdigest()[:10]


def record(path: str, cluster: str, device: str, interface: str, sha: str, date: str, subject: str, hardware: str = "", compare: str = "", read=None) -> dict:
    with open(path) as f:
        data = yaml.safe_load(f)

    invocation = (data.get("metadata") or {}).get("invocation") or []
    cases = {}
    for slug, case in (data.get("cases") or {}).items():
        grind = ((case.get("output_summary") or {}).get("simulation") or {}).get("grind")
        if isinstance(grind, (int, float)) and math.isfinite(grind):
            cases[slug] = {"grind": grind, "def": definition(case.get("description") or {}, invocation, read)}

    pr = PR.search(subject)
    return {"sha": sha, "date": date, "subject": subject, "pr": int(pr.group(1)) if pr else None, "cluster": cluster, "device": device, "interface": interface, "hardware": hardware, "env": environment(cluster, read), "covers": covers(compare) if compare else [], "cases": cases}


if __name__ == "__main__":
    if len(sys.argv) not in (8, 9, 10):
        sys.exit(__doc__)
    rec = record(*sys.argv[1:])
    if not rec["cases"]:
        sys.exit(f"{sys.argv[1]}: no case with a grind time")
    print(json.dumps(rec, separators=(",", ":"), sort_keys=True))
