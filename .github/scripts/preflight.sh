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

# --- Cluster-wide outage: requeuing cannot help, so skip rather than retry ---
if ! bash "$SCRIPT_DIR/ci-outage.sh" check "$cluster"; then
    echo "Preflight: skipping on $node because $cluster is known to be down."
    exit $EXIT_OUTAGE
fi

# --- Node health ---
syscheck_bin=$(find build/install -name syscheck -type f 2>/dev/null | head -1)

if [ -z "$syscheck_bin" ]; then
    # Nothing to probe with. A missing binary is a build problem, not a bad
    # node: requeuing would land somewhere healthy and fail the same way, so
    # let the build or test step report it instead.
    echo "Preflight: no syscheck binary under build/install; skipping node probe on $node."
    exit $EXIT_HEALTHY
fi

echo "Preflight: probing $node with $syscheck_bin"

# Launch under mpirun, not bare. Phoenix's openmpi/4.1.5 predates the PMIx that
# shipped with its Slurm upgrade, and an MPI binary started bare inside an
# allocation misreads the PMIx environment and aborts in MPI_Init.
#
# Output goes to the log verbatim. Only the exit status decides the verdict:
# PMIX_ERR_NO_PERMISSIONS and friends from dstore_base.c are benign and appear
# in more passing jobs than failing ones, so matching on log text would fail
# healthy nodes.
if mpirun -np 1 "$syscheck_bin" 2>&1; then
    echo "Preflight: $node passed."
    exit $EXIT_HEALTHY
fi

echo "::error::Preflight failed on $node: syscheck could not run MFC here."
echo "This is an INFRASTRUCTURE fault, not a code or test failure."
echo "MFC_FAULT_NODE=$node"
exit $EXIT_NODE_FAULT
