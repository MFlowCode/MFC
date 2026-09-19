#!/usr/bin/env python3
"""Assert the shared-target invariant that prebuild-case-optimization.sh relies on.

That script builds syscheck, pre_process and post_process once, from benchmarks[0], behind a
done-marker so concurrent shards skip them. That serialization only covers the benchmarks[0]
slug. It is correct while every benchmark case hashes those three targets to that same slug:
if two cases in *different* shards shared some other slug, both shards would build that
staging directory at the same time, with no marker to serialize them.

A target's slug covers Case.get_fpp, so it moves when the case changes any compile-time axis
baked into it -- _prepend() bakes chemistry and eos_state_dependent (a Mie-Gruneisen, JWL or
Vinet fluid flips the latter), and the pre_process half bakes analytic initial conditions.
None of the current benchmarks trip any of those, which is precisely the property this check
pins down so a future benchmark cannot quietly remove it.

Usage: check_caseopt_shared_slugs.py <acc|mp|no> <case.py> [<case.py> ...]
"""

import os
import sys

SHARED_TARGETS = ("syscheck", "pre_process", "post_process")


def main(argv):
    if len(argv) < 3:
        print(f"usage: {os.path.basename(argv[0])} <acc|mp|no> <case.py> [<case.py> ...]", file=sys.stderr)
        return 2

    interface, cases = argv[1], argv[2:]

    sys.path.insert(0, os.path.join(os.getcwd(), "toolchain"))
    from mfc import state
    from mfc.state import MFCConfig

    state.gCFG = MFCConfig(gpu=interface)
    # Only the keys the loader and the case files actually read: load() forwards ARGS() to each
    # case.py as its DICT, and the benchmarks read gpu/nodes/tasks_per_node from it.
    state.gARG = {
        "gpu": interface != "no",
        "nodes": 1,
        "tasks_per_node": 1,
        "case_optimization": True,
        "input": None,
        "--": [],
        "rdma_mpi": False,
    }

    from mfc.build import get_target
    from mfc.run import input as mfc_input

    slugs = {}
    for case in cases:
        state.gARG["input"] = case
        ifile = mfc_input.load(case, [], {}, do_print=False)
        slugs[case] = {name: get_target(name).get_slug(ifile) for name in SHARED_TARGETS}

    reference = cases[0]
    mismatches = []
    for case in cases[1:]:
        for name in SHARED_TARGETS:
            if slugs[case][name] != slugs[reference][name]:
                mismatches.append((name, case, slugs[case][name], slugs[reference][name]))

    if mismatches:
        print("ERROR: case-optimization benchmarks no longer share one slug per shared target.", file=sys.stderr)
        print(f"       Reference case (benchmarks[0]): {reference}", file=sys.stderr)
        for name, case, got, want in mismatches:
            print(f"         {name}: {case}", file=sys.stderr)
            print(f"           got  {got}", file=sys.stderr)
            print(f"           want {want}", file=sys.stderr)
        print(
            "       prebuild-case-optimization.sh builds these three targets once, from the reference\n"
            "       case, and lets other shards skip them -- which serializes only the reference slug.\n"
            "       Give every shard its own shared-target build, or drop the offending case, before\n"
            "       relying on the marker.",
            file=sys.stderr,
        )
        return 1

    print(f"Shared-target slugs agree across {len(cases)} benchmark case(s):")
    for name in SHARED_TARGETS:
        print(f"  {name}: {slugs[reference][name]}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
