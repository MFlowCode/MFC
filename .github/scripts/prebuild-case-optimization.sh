#!/bin/bash

# Pre-builds all benchmark cases with --case-optimization on the login node.
# Called by the case-optimization CI job before SLURM submission.
# Usage: bash .github/scripts/prebuild-case-optimization.sh <cluster> <device> <interface>

set -e

cluster=$1
job_device=$2
job_interface=$3

# Derive module flag from cluster name
case "$cluster" in
    phoenix)      flag="p" ;;
    frontier)     flag="f" ;;
    frontier_amd) flag="famd" ;;
    *) echo "ERROR: Unknown cluster '$cluster'"; exit 1 ;;
esac

. ./mfc.sh load -c "$flag" -m g
source .github/scripts/gpu-opts.sh

# Case-optimized GPU builds are memory-intensive (nvfortran/CCE + target offload).
# Login nodes have per-user cgroup memory limits (e.g., 4GB on Phoenix) that
# cause OOM kills at higher parallelism.
for case in benchmarks/*/case.py; do
    echo "=== Pre-building: $case ==="
    ./mfc.sh build -i "$case" --case-optimization $gpu_opts -j 2
done
