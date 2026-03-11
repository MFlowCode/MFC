#!/bin/bash

# Pre-builds all benchmark cases with --case-optimization.
# No GPU hardware needed — compilation only.
# Can run in two modes:
#   1. Direct (Frontier login nodes): pass cluster/device/interface as args
#   2. Inside SLURM (Phoenix): uses $job_device/$job_interface from submit-slurm-job.sh
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

rm -rf build

. ./mfc.sh load -c "$flag" -m g

# Set GPU build flags from interface — this is always a GPU build.
# Don't use gpu-opts.sh since $job_device may be "cpu" when submitted
# to a CPU SLURM partition (no GPU hardware needed for compilation).
case "$job_interface" in
    acc) gpu_opts="--gpu acc" ;;
    omp) gpu_opts="--gpu mp" ;;
    *)   echo "ERROR: prebuild requires gpu interface (acc or omp)"; exit 1 ;;
esac

for case in benchmarks/*/case.py; do
    echo "=== Pre-building: $case ==="
    ./mfc.sh build -i "$case" --case-optimization $gpu_opts -j 8
done
