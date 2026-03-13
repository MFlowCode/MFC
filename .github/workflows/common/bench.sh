#!/bin/bash
# Unified benchmark script for all clusters.
# Runs inside a SLURM job via submit-slurm-job.sh.
# Expects env vars: $job_device, $job_interface, $job_slug, $job_cluster

set -euo pipefail

source .github/scripts/bench-preamble.sh

# Cap parallel jobs at 64 to avoid overwhelming MPI daemons on large nodes
# (GNR nodes have 192 cores but nproc is too aggressive for build).
n_jobs=$(( $(nproc) > 64 ? 64 : $(nproc) ))

# --- Phoenix TMPDIR setup ---
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
# Phoenix: always nuke stale builds (heterogeneous compute nodes → ISA mismatch risk).
if [ "$job_cluster" = "phoenix" ]; then
    source .github/scripts/clean-build.sh
    clean_build
fi

if [ ! -d "build" ]; then
    source .github/scripts/retry-build.sh
    retry_build ./mfc.sh build -j $n_jobs $build_opts || exit 1
fi

# --- Bench cluster flag ---
if [ "$job_cluster" = "phoenix" ]; then
    bench_cluster="phoenix-bench"
else
    bench_cluster="$job_cluster"
fi

# --- Run benchmark ---
if [ "$job_device" = "gpu" ]; then
    ./mfc.sh bench --mem 4 -o "$job_slug.yaml" -- -c $bench_cluster $device_opts -n $n_ranks
else
    ./mfc.sh bench --mem 1 -o "$job_slug.yaml" -- -c $bench_cluster $device_opts -n $n_ranks
fi

# --- Phoenix cleanup (trap EXIT handles rm -rf "$currentdir") ---
if [ "$job_cluster" = "phoenix" ]; then
    sleep 10
    unset TMPDIR
fi
