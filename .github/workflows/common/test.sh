#!/bin/bash
# Test-only script for all clusters.
# Runs inside a SLURM job via submit-slurm-job.sh.
# Assumes MFC is already built (by a prior build.sh SLURM job).
# Expects env vars: $job_device, $job_interface, $job_shard, $job_cluster

set -euo pipefail

source .github/scripts/gpu-opts.sh
build_opts="$gpu_opts"

# --- Phoenix TMPDIR setup ---
if [ "$job_cluster" = "phoenix" ]; then
    tmpbuild=/storage/project/r-sbryngelson3-0/sbryngelson3/mytmp_build
    currentdir=$tmpbuild/run-$(( RANDOM % 9000 ))
    mkdir -p $tmpbuild
    mkdir -p $currentdir
    export TMPDIR=$currentdir
    trap 'rm -rf "$currentdir" || true' EXIT
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

# Only prune tests on PRs; master pushes must run the full suite.
prune_flag=""
if [ "${GITHUB_EVENT_NAME:-}" = "pull_request" ]; then
    prune_flag="--only-changes"
fi

./mfc.sh test -v --max-attempts 3 --no-build $prune_flag -a -j $n_test_threads $rdma_opts $device_opts $build_opts $shard_opts -- -c $job_cluster
