#!/bin/bash
# Prove this node can run MFC before spending an allocation on it.
#
# Runs syscheck -- the same binary the suite already builds -- as the first
# thing inside the allocation that will execute the tests. syscheck initialises
# MPI, creates a device context, launches a kernel and reads the result back, so
# a node with a dead GPU, a broken MPI layer, or a binary built for a different
# microarchitecture fails here in seconds instead of after the build.
#
# Why at the top of the *test* allocation: over 2026-08-18..31, 48 of 58
# measurable jobs built on one node and tested on another, and in every ECC
# failure the build node was healthy while the tests landed on a bad one. A
# probe placed after the build checks the wrong machine. (Phoenix's combined
# build-and-test allocation avoids the split; Frontier still submits the two
# separately and lands elsewhere ~90% of the time.)
#
# Usage: preflight.sh <cluster> <device>
#
# Exit codes:
#   0   node looks healthy, carry on
#   77  node-local fault -- caller should exclude this node and resubmit

set -uo pipefail

cluster="${1:-}"
device="${2:-}"

if [ -z "$cluster" ] || [ -z "$device" ]; then
    echo "Usage: $0 <cluster> <device>" >&2
    exit 2
fi

EXIT_HEALTHY=0
EXIT_NODE_FAULT=77

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
node="${SLURMD_NODENAME:-$(hostname -s 2>/dev/null || hostname)}"

# Only judge a node from inside its own allocation. `mfc.sh load` is also used
# for building on login nodes -- bench.yml and frontier/build.sh both load the
# GPU module set there -- and a login node has no GPU to probe. Reporting a node
# fault in that context would have the wrapper exclude a login node and requeue,
# which is both wrong and hard to diagnose. Nothing calls the probe from a login
# node today; this makes that a property of the probe rather than a convention
# every future caller has to remember.
if [ -z "${SLURM_JOB_ID:-}" ]; then
    echo "Preflight: not inside a SLURM allocation; skipping the node probe."
    exit $EXIT_HEALTHY
fi

# --- Node health ---
# Pick the *newest* install matching this job's device (build/install is named
# e.g. gpu-acc-<hash>, gpu-mp-<hash>). Both halves matter: the device filter
# avoids probing another variant's binary, and newest-wins avoids probing a
# leftover from an earlier job. Not every caller nukes build/ first -- bench.sh
# only does so on Phoenix -- and a stale binary compiled for a different
# microarchitecture dies with SIGILL, which would be reported as a bad node and
# get a perfectly healthy one excluded.
newest_syscheck() {
    find "$@" -name syscheck -type f -printf '%T@ %p\n' 2>/dev/null \
        | sort -rn | head -1 | cut -d' ' -f2-
}

syscheck_bin=$(newest_syscheck build/install -path "*${device}*")
if [ -z "$syscheck_bin" ]; then
    syscheck_bin=$(newest_syscheck build/install)
fi

if [ -z "$syscheck_bin" ]; then
    # Nothing to probe with. A missing binary is a build problem, not a bad
    # node: requeuing would land somewhere healthy and fail the same way, so
    # let the build or test step report it instead.
    echo "Preflight: no syscheck binary under build/install; skipping node probe on $node."
    exit $EXIT_HEALTHY
fi

# A binary built for a device this allocation did not ask for tells us nothing
# about the node. The case-optimization pre-build is submitted as cpu (it is a
# --dry-run that only builds) while producing GPU binaries, so its syscheck
# asserts a device exists and fails on a GPU-less node by design. Read as a node
# fault, that condemned three healthy Phoenix nodes and excluded two of them.
# The mismatch is a property of the job, never of the machine.
case "$syscheck_bin" in
    *gpu-*) binary_device="gpu" ;;
    *)      binary_device="cpu" ;;
esac
if [ "$device" = "cpu" ] && [ "$binary_device" = "gpu" ]; then
    echo "Preflight: $device allocation but the available syscheck is a GPU build"
    echo "  ($syscheck_bin); it cannot pass here, so skipping rather than judging $node."
    exit $EXIT_HEALTHY
fi

echo "Preflight: probing $node with $syscheck_bin"

# Launch the probe the way this cluster launches everything else. Phoenix uses
# mpirun -- its openmpi predates the PMIx that shipped with its Slurm upgrade,
# so a bare MPI binary misreads the environment and aborts in MPI_Init. Frontier
# and frontier_amd use srun and Cray MPICH ships no mpirun at all, so running
# one there fails 127 no matter how healthy the node is. See
# toolchain/templates/{phoenix,frontier,frontier_amd}.mako.
case "$cluster" in
    # --bind-to none: a single-rank health probe has nothing to bind against,
    # and Open MPI's default binding fails outright on some Phoenix nodes
    # ("hwloc_set_cpubind returned Error for bitmap 0"), killing the process
    # before the binary is even launched. That is a launcher problem, not a
    # node problem -- but it condemned three healthy nodes before being caught.
    phoenix)               launcher=(mpirun --bind-to none -np 1) ;;
    frontier|frontier_amd) launcher=(srun -n1) ;;
    *)                     launcher=() ;;
esac

# A launcher missing from PATH says nothing about the node. Probing bare is a
# weaker test, but calling a healthy node bad is far worse: it costs three
# allocations and blacklists three good nodes before giving up.
if [ "${#launcher[@]}" -gt 0 ] && ! command -v "${launcher[0]}" >/dev/null 2>&1; then
    echo "Preflight: ${launcher[0]} is not on PATH; probing without a launcher."
    launcher=()
fi

# Output goes to the log verbatim. Only the exit status decides the verdict:
# PMIX_ERR_NO_PERMISSIONS and friends from dstore_base.c are benign and appear
# in more passing jobs than failing ones, so matching on log text would fail
# healthy nodes.
# Captured to a variable, not a temp file: this runs before any module set is
# guaranteed and mktemp is not always on PATH here.
run_probe() {
    probe_rc=0
    if [ "$#" -eq 0 ]; then
        probe_out=$("$syscheck_bin" 2>&1) || probe_rc=$?
    else
        probe_out=$("$@" "$syscheck_bin" 2>&1) || probe_rc=$?
    fi
}

run_probe "${launcher[@]}"

# If this launcher does not take the flags we added, drop them and probe again
# rather than reporting a verdict about the node. Otherwise a launcher that
# rejects an option would fail every probe, and -- because a failed launch is
# treated as inconclusive below -- would silently switch the preflight off
# instead of failing loudly.
case "$probe_out" in
    *"unrecognized option"*|*"unrecognized argument"*|*"Unknown option"*|*"invalid option"*)
        if [ "${#launcher[@]}" -gt 1 ]; then
            echo "Preflight: ${launcher[0]} rejected the probe's options; retrying with none of them."
            run_probe "${launcher[0]}"
        fi
        ;;
esac

printf '%s\n' "$probe_out"

if [ "$probe_rc" -eq 0 ]; then
    echo "Preflight: $node passed."
    exit $EXIT_HEALTHY
fi

# Only a binary that RAN and failed says anything about this node. When the
# launcher never got as far as starting it, the verdict is about mpirun or the
# allocation, and excluding the node is both wrong and expensive -- three
# healthy Phoenix nodes were excluded this way, two jobs deep, before the run
# gave up. Judge nothing on a launch that never happened.
case "$probe_out" in
    *"The specified application failed to start"*|*"unable to start the specified application"*|*"was killed without launching the target application"*)
        echo "Preflight: the launcher could not start $syscheck_bin on $node;"
        echo "  that is a launcher or allocation problem, not evidence about the node. Continuing."
        exit $EXIT_HEALTHY
        ;;
esac

echo "::error::Preflight failed on $node: syscheck could not run MFC here."
echo "This is an INFRASTRUCTURE fault, not a code or test failure."
echo "MFC_FAULT_NODE=$node"
exit $EXIT_NODE_FAULT
