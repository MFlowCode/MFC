#!/bin/bash
# Unified test script for all clusters.
# Runs inside a SLURM job via submit-slurm-job.sh.
# Expects env vars: $job_device, $job_interface, $job_shard, $job_cluster

set -euo pipefail

source .github/scripts/gpu-opts.sh
build_opts="$gpu_opts"

# --- Phoenix TMPDIR setup ---
# Phoenix compute nodes have a small /tmp. With 8 parallel test threads each
# spawning MPI processes, it fills up and ORTE session dir creation fails.
# Redirect TMPDIR to project storage, same as bench.sh.
if [ "$job_cluster" = "phoenix" ]; then
    tmpbuild=/storage/project/r-sbryngelson3-0/sbryngelson3/mytmp_build
    currentdir=$tmpbuild/run-$(( RANDOM % 9000 ))
    mkdir -p $tmpbuild
    mkdir -p $currentdir
    export TMPDIR=$currentdir
    trap 'rm -rf "$currentdir" || true' EXIT
fi

# --- Build (if not pre-built on login node) ---
# Phoenix builds inside SLURM; Frontier pre-builds via build.sh on the login node.
# Phoenix builds inside SLURM on heterogeneous compute nodes — always start fresh
# to avoid SIGILL from stale binaries compiled on a different microarchitecture.
if [ "$job_cluster" = "phoenix" ]; then
    source .github/scripts/clean-build.sh
    clean_build
fi

if [ ! -d "build" ]; then
    source .github/scripts/retry-build.sh

    # Phoenix: smoke-test the syscheck binary to catch architecture mismatches
    # (SIGILL from binaries compiled on a different compute node).
    validate_cmd=""
    if [ "$job_cluster" = "phoenix" ]; then
        validate_cmd='syscheck_bin=$(find build/install -name syscheck -type f 2>/dev/null | head -1); [ -z "$syscheck_bin" ] || "$syscheck_bin" > /dev/null 2>&1'
    fi

    RETRY_VALIDATE_CMD="$validate_cmd" \
        retry_build ./mfc.sh test -v --dry-run -j 8 $build_opts || exit 1
fi

# --- GPU detection and thread count ---
device_opts=""
rdma_opts=""
shard_opts=""

case "$job_cluster" in
    phoenix)      n_test_threads=8 ;;
    *)            n_test_threads=32 ;;
esac

if [ "$job_device" = "gpu" ]; then
    source .github/scripts/detect-gpus.sh

    case "$job_cluster" in
        phoenix)
            device_opts="-g $gpu_ids"
            n_test_threads=$((ngpus * 2))
            ;;
        *)
            # Frontier: --gpu flag is already in $build_opts; no extra device opts needed
            device_opts=""
            n_test_threads=$ngpus
            ;;
    esac

    # RDMA for Frontier CCE (not frontier_amd)
    if [ "$job_cluster" = "frontier" ]; then
        rdma_opts="--rdma-mpi"
    fi
else
    device_opts="--no-gpu"
fi

# --- Sharding (Frontier only) ---
if [ -n "${job_shard:-}" ]; then
    shard_opts="--shard $job_shard"
fi

./mfc.sh test -v --max-attempts 3 -a -j $n_test_threads $rdma_opts $device_opts $build_opts $shard_opts -- -c $job_cluster
