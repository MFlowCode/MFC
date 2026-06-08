#!/bin/bash

# Pre-builds all benchmark cases with --case-optimization using --dry-run so
# binaries are cached before the GPU run job. No simulation is executed.
# Can run in two modes:
#   1. Direct (Frontier login nodes): pass cluster/device/interface as args
#   2. Inside SLURM (Phoenix/frontier_amd): uses $job_device/$job_interface
# Usage: bash prebuild-case-optimization.sh [<cluster> <device> <interface>]

set -e

# Support both positional args (direct invocation) and env vars (SLURM)
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

# Phoenix starts fresh (no prior dep build); other clusters pre-build deps via
# build.sh first, so we must preserve them and only clean MFC target staging.
if [ "$cluster" = "phoenix" ]; then
    source .github/scripts/clean-build.sh
    clean_build
else
    find build/staging -maxdepth 1 -regex '.*/[0-9a-f]+' -type d -exec rm -rf {} + 2>/dev/null || true
    find build/install -maxdepth 1 -regex '.*/[0-9a-f]+' -type d -exec rm -rf {} + 2>/dev/null || true
fi

. ./mfc.sh load -c "$flag" -m g

case "$job_interface" in
    acc) gpu_opts="--gpu acc" ;;
    omp) gpu_opts="--gpu mp" ;;
    *)   echo "ERROR: prebuild requires gpu interface (acc or omp)"; exit 1 ;;
esac

for case in benchmarks/*/case.py; do
    echo "=== Pre-building: $case ==="
    ./mfc.sh run "$case" --case-optimization $gpu_opts -j 8 --dry-run
done
