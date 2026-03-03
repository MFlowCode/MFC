#!/bin/bash

# Pre-builds all benchmark cases with --case-optimization.
# Can run in two modes:
#   1. Direct (Frontier login nodes): pass cluster/device/interface as args
#   2. Inside SLURM (Phoenix): uses $job_device/$job_interface from submit.sh
# Usage: bash prebuild-case-optimization.sh [<cluster> <device> <interface>]

set -e

# Support both positional args (direct invocation) and env vars (SLURM via submit.sh)
cluster="${1:-${job_cluster:-phoenix}}"
job_device="${2:-$job_device}"
job_interface="${3:-$job_interface}"

# Derive module flag from cluster name
case "$cluster" in
    phoenix)      flag="p" ;;
    frontier)     flag="f" ;;
    frontier_amd) flag="famd" ;;
    *) echo "ERROR: Unknown cluster '$cluster'"; exit 1 ;;
esac

. ./mfc.sh load -c "$flag" -m g
source .github/scripts/gpu-opts.sh

for case in benchmarks/*/case.py; do
    echo "=== Pre-building: $case ==="
    ./mfc.sh build -i "$case" --case-optimization $gpu_opts -j 8
done
