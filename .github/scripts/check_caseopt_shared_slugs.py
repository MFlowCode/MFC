#!/usr/bin/env python3
"""Assert every benchmark case shares benchmarks[0]'s slug for syscheck, pre_process and post_process.

prebuild-case-optimization.sh builds those once from benchmarks[0] behind a done-marker; a case with a different
slug (changed chemistry, eos_state_dependent or analytic ICs) would build its staging dir concurrently, unserialized.

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
    # Only keys load() forwards to case.py as DICT (benchmarks read gpu/nodes/tasks_per_node).
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
