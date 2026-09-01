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
#   78  cluster-wide outage already recorded -- caller should skip, not requeue

set -uo pipefail

cluster="${1:-}"
device="${2:-}"

if [ -z "$cluster" ] || [ -z "$device" ]; then
    echo "Usage: $0 <cluster> <device>" >&2
    exit 2
fi

EXIT_HEALTHY=0
EXIT_NODE_FAULT=77
EXIT_OUTAGE=78

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

# --- Cluster-wide outage: requeuing cannot help, so skip rather than retry ---
outage_rc=0
bash "$SCRIPT_DIR/ci-outage.sh" check "$cluster" || outage_rc=$?
if [ "$outage_rc" -eq 1 ]; then
    echo "Preflight: skipping on $node because $cluster is known to be down."
    exit $EXIT_OUTAGE
elif [ "$outage_rc" -ne 0 ]; then
    # Only exit 1 means "tripped". Anything else means the breaker could not be
    # read at all (missing script, unreadable state dir), which says nothing
    # about the cluster -- treating it as an outage would halt CI on a bug here.
    echo "Preflight: could not read the outage breaker (exit $outage_rc); continuing."
fi

# --- Node health ---
# Prefer an install matching this job's device (build/install is named e.g.
# gpu-acc-<hash>, gpu-mp-<hash>), so a leftover install from another variant is
# not probed instead. Fall back to any syscheck if none matches.
syscheck_bin=$(find build/install -path "*${device}*" -name syscheck -type f 2>/dev/null | head -1)
if [ -z "$syscheck_bin" ]; then
    syscheck_bin=$(find build/install -name syscheck -type f 2>/dev/null | head -1)
fi

if [ -z "$syscheck_bin" ]; then
    # Nothing to probe with. A missing binary is a build problem, not a bad
    # node: requeuing would land somewhere healthy and fail the same way, so
    # let the build or test step report it instead.
    echo "Preflight: no syscheck binary under build/install; skipping node probe on $node."
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
    phoenix)               launcher=(mpirun -np 1) ;;
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
probe_rc=0
if [ "${#launcher[@]}" -eq 0 ]; then
    "$syscheck_bin" 2>&1 || probe_rc=$?
else
    "${launcher[@]}" "$syscheck_bin" 2>&1 || probe_rc=$?
fi

if [ "$probe_rc" -eq 0 ]; then
    echo "Preflight: $node passed."
    exit $EXIT_HEALTHY
fi

echo "::error::Preflight failed on $node: syscheck could not run MFC here."
echo "This is an INFRASTRUCTURE fault, not a code or test failure."
echo "MFC_FAULT_NODE=$node"
exit $EXIT_NODE_FAULT
