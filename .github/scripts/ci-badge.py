#!/usr/bin/env python3
"""Write the README's status badge from the CI history.

    ci-badge.py <MFC-ci-data checkout>

Prints a shields.io endpoint object: how many cluster lanes pass on master, taking each
lane's latest result that is not a cancellation, over the last 30 days. A lane is a
machine, compiler and offload model, counted once for its tests and once for its
case-optimization job (as docs/status.html counts them).
"""

import glob
import json
import os
import sys
from datetime import datetime, timedelta

RANK = {"success": 0, "cancelled": 1, "timed_out": 2, "failure": 3}


def lanes(root: str) -> dict:
    """Lane -> (when, passed) of its latest master result."""
    since = (datetime.utcnow() - timedelta(days=30)).strftime("%Y-%m-%dT%H:%M:%SZ")
    latest = {}
    for path in sorted(glob.glob(os.path.join(root, "ci", "*.jsonl")))[-2:]:
        with open(path) as f:
            for line in f:
                run = json.loads(line)
                if not run.get("master") or run.get("workflow") != "test" or run["date"] < since:
                    continue
                attempt = max(j["a"] for j in run["jobs"])
                cells = {}
                for j in run["jobs"]:
                    if j["a"] == attempt and j["k"] in ("test", "case-opt"):
                        key = f"{j['l']}|{j['k']}"
                        cells[key] = max(cells.get(key, "success"), j["c"], key=lambda c: RANK.get(c, 1))
                for key, worst in cells.items():
                    if worst != "cancelled" and (key not in latest or run["date"] >= latest[key][0]):
                        latest[key] = (run["date"], worst == "success")
    return latest


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    state = lanes(sys.argv[1])
    ok, n = sum(passed for _, passed in state.values()), len(state)
    color = "lightgrey" if not n else "brightgreen" if ok == n else "orange" if ok * 2 >= n else "red"
    print(json.dumps({"schemaVersion": 1, "label": "master", "message": f"{ok}/{n} lanes passing" if n else "no runs", "color": color}))
