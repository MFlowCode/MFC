#!/usr/bin/env python3
"""Turn one scaling run's logs into one line of the scaling history.

    scaling-record.py <log dir> <sha> <date> <subject>

The directory holds <weak|strong>-<nodes>.out, the simulation logs scaling.sh leaves. Prints one
JSON object, or exits non-zero when no log has a result:

    sha, date, subject    the master commit
    cluster, interface    where and how it ran
    runs                  per log: kind (weak|strong), nodes, ranks, cells (in all), grind (ns per
                          grid point, equation and RK stage) and step (seconds per time step)

Grind time is wall time per grid point over all ranks, so with perfect scaling it falls as
1/ranks in both cases (weak holds the cells per rank fixed, strong the cells in all), and the
efficiency against the smallest run is grind(smallest)·ranks(smallest)/(grind(N)·ranks(N)).
docs/dashboard/benchmarks.html draws it.
"""

import glob
import json
import os
import re
import sys

LOG = re.compile(r"(weak|strong)-(\d+)\.out$")
NUMBER = r"\s*([0-9.]+(?:[eE][+-]?\d+)?)"


def last(pattern: str, text: str):
    found = re.findall(pattern, text)
    return float(found[-1]) if found else None


def runs(root: str) -> list:
    out = []
    for path in sorted(glob.glob(os.path.join(root, "*.out"))):
        m = LOG.search(path)
        if not m:
            continue
        with open(path, errors="replace") as f:
            text = f.read()
        grind, step = last(r"Performance:" + NUMBER, text), last(r"\(avg" + NUMBER + r"s\)", text)
        size = re.search(r"Simulating a [\w-]+ (\d+)x(\d+)x(\d+) case on (\d+) rank", text)
        if grind is None or step is None or not size:
            continue
        cells = int(size.group(1)) * int(size.group(2)) * int(size.group(3))
        out.append({"kind": m.group(1), "nodes": int(m.group(2)), "ranks": int(size.group(4)), "cells": cells, "grind": grind, "step": step})
    return sorted(out, key=lambda r: (r["kind"], r["nodes"]))


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    found = runs(sys.argv[1])
    if not found:
        sys.exit(f"{sys.argv[1]}: no log with a result")
    print(json.dumps({"sha": sys.argv[2], "date": sys.argv[3], "subject": sys.argv[4], "cluster": "frontier_amd", "interface": "omp", "runs": found}, separators=(",", ":")))
