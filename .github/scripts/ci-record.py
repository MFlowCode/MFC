#!/usr/bin/env python3
"""Turn one workflow run into one line of the CI history.

    ci-record.py <run.json> <jobs.json> [--usage usage.txt] [--compare compare.json] [--logs DIR]

The inputs are the GitHub API's answers for the run and for its jobs (`actions/runs/<id>` and
`actions/runs/<id>/jobs?filter=all`). Prints nothing for a run with no cluster job. Otherwise
one JSON object: the run (event, commit, title, pull request for a master commit) and, per
cluster job, a compact entry:

    l  lane, e.g. frontier_amd-gpu-omp      k  kind: test, case-opt or bench
    s  shard ("1/2") or null                m  "smoke" or "full"
    a  run attempt                          c  conclusion
    w  seconds waiting for a runner         d  seconds from start to finish
    id GitHub job id                        u  SLURM jobs: [name, nodes, elapsed s, queued s, state]

and, for a test job whose log is in DIR as <job id>.txt:

    n    tests that passed
    f    tests that failed: [uuid, name, error class]
    r    tests that failed and then passed on a retry: [uuid, name, error class]

and, for any failed job whose log is there:

    why  what failed: tests, build, cluster (down), queue (never started), timeout (SLURM), node, checkout (git),
         upload (artifact), runner, or other. All but tests, build and other are the machine
         or the plumbing, not the change: ci-record.yml re-runs a run whose failures are all those.

--usage: lines "<GitHub job id> <slurm-usage notice>", the notices run_monitored_slurm_job.sh
leaves on each job. --compare: the API's `compare/<previous run's commit>...<this run's commit>`;
a master run then lists the `merges` it is the first to test, as [sha, title, pull request]: a
scheduled run covers every merge since the scheduled run before it.
The history feeds docs/ci.html and docs/tests.html.
"""

import argparse
import json
import os
import re
import sys
from datetime import datetime

LANE = re.compile(r"(Frontier|Phoenix) \((CCE|AMD|NVHPC)\) \((gpu|cpu)-(acc|omp|none)(?: \[(\d+/\d+)\])?\)")
CLUSTERS = {("Frontier", "CCE"): "frontier", ("Frontier", "AMD"): "frontier_amd", ("Phoenix", "NVHPC"): "phoenix"}


def seconds(start, end):
    if not start or not end:
        return None
    t = [datetime.strptime(x, "%Y-%m-%dT%H:%M:%SZ") for x in (start, end)]
    return max(0, int((t[1] - t[0]).total_seconds()))


def usage(lines) -> dict:
    """GitHub job id -> its SLURM jobs, from 'id|name|nodes|elapsed|submit|start|state' notices."""
    out = {}
    for line in lines:
        gh_job, _, notice = line.strip().partition(" ")
        f = notice.split("|")
        if len(f) != 7 or not f[2].isdigit() or not f[3].isdigit():
            continue
        queued = int(f[5]) - int(f[4]) if f[4].isdigit() and f[5].isdigit() else None
        out.setdefault(gh_job, []).append([f[1], int(f[2]), int(f[3]), queued, f[6].split()[0]])
    return out


PR = re.compile(r"\(#(\d+)\)\s*$")  # a squash merge's title ends with its pull request number


def merges(compare: dict) -> list:
    out = []
    for c in (compare or {}).get("commits", []):
        title = (c.get("commit", {}).get("message") or "").split("\n")[0]
        pr = PR.search(title)
        out.append([c.get("sha"), title, int(pr.group(1)) if pr else None])
    return out


STAMP = re.compile(r"^\d{4}-\d\d-\d\dT[\d:.]+Z ")
ANSI = re.compile(r"\x1b\[[0-9;]*m")
PASSED = re.compile(r"^\s*\(\s*\d+/\d+\)\s+(.*?)\s+\d+\.\d+\s+([0-9A-F]{8})\s*$")
RECOVERED = re.compile(r"^\s*recovered on attempt \d+ \(([^)]*)\): (.*)$")


def error_class(text: str) -> str:
    """The test runner's own buckets (toolchain/mfc/test: classify_error), from an error message."""
    t = text.lower()
    for needle, name in (("tolerance", "tolerance mismatch"), ("golden", "tolerance mismatch"), ("timeout", "timeout"), ("nan", "NaN detected"), ("memory access fault by gpu", "GPU memory fault"), ("failed", "execution failed")):
        if needle in t:
            return name
    return "other"


def tests(log: str, failed_job: bool) -> dict:
    """What a test job's log says about its tests, and for a failed job what failed."""
    lines = [ANSI.sub("", STAMP.sub("", line)) for line in log.splitlines()]
    passed, failed, recovered = {}, {}, {}
    for i, line in enumerate(lines):
        m = PASSED.match(line)
        if m:
            passed[m.group(2)] = m.group(1)
            continue
        m = RECOVERED.match(line)
        if m:
            recovered[m.group(2).strip()] = m.group(1)
            continue
        if "\u2717 FAILED:" in line:
            name, block = line.split("FAILED:", 1)[1].strip(), "\n".join(lines[i + 1 : i + 5])
            uuid, err = re.search(r"UUID: ([0-9A-F]{8})", block), re.search(r"Error: (.*)", block)
            if uuid:
                failed[uuid.group(1)] = [uuid.group(1), name, error_class(err.group(1) if err else "")]
    # A retry line names the test in full; the pass line that follows has its id and a name cut to fit.
    retried = []
    for name, cls in recovered.items():
        uuid = next((u for u, shown in passed.items() if name.startswith(shown.rstrip(".").rstrip())), "")
        retried.append([uuid, name, cls])

    out = {"n": len(passed), "f": list(failed.values()), "r": retried}
    if failed_job:
        text = "\n".join(lines)
        out["why"] = (
            "tests" if failed or "CRITICAL:" in text
            else "cluster" if " unavailable::" in text
            else "queue" if "QUEUE STARVATION" in text
            else "build" if "Build failed after" in text
            else "timeout" if "state=TIMEOUT" in text
            else "node" if "Preflight failed" in text or "failed preflight" in text
            else "checkout" if "The process '/usr/bin/git' failed" in text
            else "upload" if "Failed to CreateArtifact" in text or "Failed to FinalizeArtifact" in text
            else "runner" if re.search(r"runner has received a shutdown signal|operation was canceled|has timed out after|transport endpoint", text)
            else "other"
        )
    return out


def record(run: dict, jobs: list, slurm: dict = None, compare: dict = None, logs: str = None):
    entries = []
    for job in jobs:
        name = job.get("name") or ""
        m = LANE.search(name)
        if not m or job.get("conclusion") in (None, "skipped"):
            continue
        entry = {}
        log = os.path.join(logs, f"{job.get('id')}.txt") if logs else ""
        if os.path.isfile(log):
            with open(log, errors="replace") as f:
                entry = tests(f.read(), job["conclusion"] == "failure")
            if name.startswith(("Bench", "Case Opt")):  # not a test job: only why it failed
                entry = {k: v for k, v in entry.items() if k == "why"}
        entries.append(
            {
                "l": f"{CLUSTERS[(m.group(1), m.group(2))]}-{m.group(3)}-{m.group(4)}",
                "k": "bench" if name.startswith("Bench") else "case-opt" if name.startswith("Case Opt") else "test",
                "s": m.group(5),
                "m": "smoke" if "[smoke]" in name else "full",
                "a": job.get("run_attempt", 1),
                "c": job["conclusion"],
                "w": seconds(job.get("created_at"), job.get("started_at")),
                "d": seconds(job.get("started_at"), job.get("completed_at")),
                "id": job.get("id"),
                "u": (slurm or {}).get(str(job.get("id")), []),
                **entry,
            }
        )
    if not entries:
        return None

    master = run.get("head_branch") == "master" and run.get("event") in ("push", "schedule", "workflow_dispatch")
    pr = PR.search(run.get("display_title") or "") if master else None
    return {
        "run": run["id"],
        "workflow": "bench" if run.get("path", "").endswith("bench.yml") else "test",
        "event": run.get("event"),
        "master": master,
        "sha": run.get("head_sha"),
        "title": run.get("display_title"),
        "pr": int(pr.group(1)) if pr else None,
        "date": run.get("created_at"),
        "merges": merges(compare) if master else [],
        "jobs": entries,
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("jobs")
    ap.add_argument("--usage")
    ap.add_argument("--compare")
    ap.add_argument("--logs")
    args = ap.parse_args()

    def load(path, default):
        if not path:
            return default
        with open(path) as f:
            return json.load(f)

    slurm = {}
    if args.usage:
        with open(args.usage) as f:
            slurm = usage(f)
    rec = record(load(args.run, {}), load(args.jobs, {})["jobs"], slurm, load(args.compare, None), args.logs)
    if rec:
        print(json.dumps(rec, separators=(",", ":"), sort_keys=True))
