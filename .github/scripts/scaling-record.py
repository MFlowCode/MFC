#!/usr/bin/env python3
"""Turn one scaling run's logs into one line of the scaling history.

    scaling-record.py <log dir> <sha> <date> <subject>

The directory holds <weak|strong>-<nodes>.out, the simulation logs scaling.sh leaves. Prints one
JSON object, or exits non-zero when no log has a result:

    sha, date, subject    the master commit
    cluster, interface    where and how it ran
    runs                  per log: kind (weak|strong), nodes, ranks, cells (in all), grind (ns per
                          grid point, equation and RK stage) and step (seconds per time step)

Weak scaling holds the cells per rank fixed, so its efficiency against the smallest run is
step(smallest)/step(N); strong scaling holds the cells fixed, so its efficiency is
step(smallest)·nodes(smallest)/(step(N)·N). docs/benchmarks.html draws both.
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
        grind, step = last(r"Performance:" + NUMBER, text), last(r"Time Avg =" + NUMBER, text)
        if grind is None or step is None:
            continue
        cells = [int(v) + 1 for k in "mnp" for v in re.findall(rf"^\s*{k}\s*=?\s*(\d+)\s*$", text, re.M)[:1]]
        size = cells[0] * cells[1] * cells[2] if len(cells) == 3 else None
        out.append({"kind": m.group(1), "nodes": int(m.group(2)), "ranks": int(m.group(2)) * 8, "cells": size, "grind": grind, "step": step})
    return sorted(out, key=lambda r: (r["kind"], r["nodes"]))


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    found = runs(sys.argv[1])
    if not found:
        sys.exit(f"{sys.argv[1]}: no log with a result")
    print(json.dumps({"sha": sys.argv[2], "date": sys.argv[3], "subject": sys.argv[4], "cluster": "frontier", "interface": "acc", "runs": found}, separators=(",", ":")))
