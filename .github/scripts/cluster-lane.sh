#!/bin/bash
# One cluster lane's SLURM work, shared by test.yml and cluster-bisect.sh.
#   cluster-lane.sh build <device> <interface> <cluster>
#   cluster-lane.sh test  <device> <interface> <cluster> [shard ...]
# build: one job, except Frontier AMD, whose variants (MFC_BUILD_VARIANTS, default "base chem-h2o2 chem-carbon eos")
#        are concurrent jobs so the amdflang LTO links do not serialize past the 2 h walltime, followed by a
#        short job that builds any chemistry mechanism those did not name.
# test:  one job per shard, concurrently, from this workspace and the binaries build left in it
#        (test.sh passes --no-build; each shard writes only its own tests/<UUID> directories).
#        No shard: one job. MFC_TEST_EXTRA_OPTS is appended to the test command.
set -euo pipefail

step=$1 device=$2 interface=$3 cluster=$4
shift 4
submit=(bash .github/scripts/submit-slurm-job.sh)

# Phoenix tests run on the partition their build ran on (see below).
if [ "$step" = test ] && [ "$cluster" = phoenix ]; then
    MFC_GPU_PARTITION="$(cat .ci-build-partition)"
    export MFC_GPU_PARTITION
fi

# A closed partition is the machine's failure, not the change's: say so before queueing.
case "$cluster" in
    phoenix) partitions="${MFC_GPU_PARTITION:-gpu-h200,gpu-h100,gpu-a100,gpu-v100}" ;;
    *)       partitions="g1" ;;
esac
if avail=$(sinfo -h -p "$partitions" -o %a 2>/dev/null) && [ -n "$avail" ] && ! grep -q up <<<"$avail"; then
    # ci-record.py classes the job as a cluster failure by "cluster unavailable" in the log, where
    # an annotation's title does not appear.
    echo "::error title=$cluster unavailable::cluster unavailable: SLURM partition $partitions on $cluster is not up. This is the machine, not this change."
    exit 78
fi

wait_all() {
    local rc=0 p
    for p in "$@"; do wait "$p" || rc=1; done
    return $rc
}

case "$step" in
    build)
        if [ "$cluster" = "frontier_amd" ]; then
            pids=()
            for v in ${MFC_BUILD_VARIANTS:-base chem-h2o2 chem-carbon eos}; do
                "${submit[@]}" .github/workflows/common/build.sh "$device" "$interface" "$cluster" "" "$v" &
                pids+=($!)
            done
            wait_all "${pids[@]}"
            # Sweep: `-o Chemistry` finds the named mechanisms built and builds any other, so a test
            # with a new mechanism needs no change here (lint_test_suite.py counts on that).
            "${submit[@]}" .github/workflows/common/build.sh "$device" "$interface" "$cluster" "" chem
        else
            # build.sh records the partition it ran on: a binary built on one Phoenix node type can
            # SIGILL on another, so the tests are pinned to it.
            rm -f .ci-build-partition
            "${submit[@]}" .github/workflows/common/build.sh "$device" "$interface" "$cluster"
        fi
        ;;
    test)
        rm -f tests/failed_uuids*.txt
        [ $# -gt 0 ] || set -- ""
        pids=()
        for s in "$@"; do
            "${submit[@]}" .github/workflows/common/test.sh "$device" "$interface" "$cluster" "$s" &
            pids+=($!)
        done
        rc=0
        wait_all "${pids[@]}" || rc=1
        # Read from the job log by phoenix-master-alert.sh to bisect a red master.
        if failed=$(cat tests/failed_uuids*.txt 2>/dev/null) && [ -n "$failed" ]; then
            echo "FAILED_TESTS: $(echo $failed)"
        fi
        exit $rc
        ;;
    *)
        echo "Usage: $0 build|test <device> <interface> <cluster> [shard ...]"
        exit 2
        ;;
esac
