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

# Optional sharding (format "i/N", e.g. "1/2"), set by submit-slurm-job.sh's
# [shard] argument via $job_shard: shard i builds every Nth case of the sorted
# case list. Unset = build all cases in one job (default; other clusters).
shard="${job_shard:-}"
if [ -n "$shard" ]; then
    shard_idx="${shard%%/*}"
    shard_count="${shard##*/}"
    case "${shard_idx}${shard_count}" in
        ''|*[!0-9]*) echo "ERROR: bad shard '$shard' (expected i/N)"; exit 1 ;;
    esac
    if [ "$shard" != "$shard_idx/$shard_count" ] || [ "$shard_idx" -lt 1 ] || [ "$shard_idx" -gt "$shard_count" ]; then
        echo "ERROR: bad shard '$shard' (expected i/N with 1 <= i <= N)"; exit 1
    fi
fi

# Phoenix starts fresh (no prior dep build); other clusters pre-build deps via
# build.sh first, so we must preserve them and only clean MFC target staging.
# Sharded jobs share one workspace and run concurrently, so the workflow
# cleans once before submitting them — cleaning here would wipe a sibling
# shard's in-progress build.
if [ "$cluster" = "phoenix" ]; then
    source .github/scripts/clean-build.sh
    clean_build
elif [ -z "$shard" ]; then
    find build/staging -maxdepth 1 -regex '.*/[0-9a-f]+' -type d -exec rm -rf {} + 2>/dev/null || true
    find build/install -maxdepth 1 -regex '.*/[0-9a-f]+' -type d -exec rm -rf {} + 2>/dev/null || true
fi

. ./mfc.sh load -c "$flag" -m g

case "$job_interface" in
    acc) gpu_opts="--gpu acc" ;;
    omp) gpu_opts="--gpu mp" ;;
    *)   echo "ERROR: prebuild requires gpu interface (acc or omp)"; exit 1 ;;
esac

# Case-optimized simulation builds land in per-case hash-named staging dirs,
# but syscheck/pre_process/post_process hash identically across these cases.
# Concurrent shards must not build those shared staging dirs simultaneously:
# shard 1 builds them first; the other shards wait for its marker, after
# which their builds no-op in the shared dirs.
if [ -n "$shard" ] && [ "$shard_count" -gt 1 ]; then
    shared_marker="build/.prebuild-shared-targets-done"
    set -- benchmarks/*/case.py
    first_case="$1"
    if [ "$shard_idx" -eq 1 ]; then
        echo "=== Shard 1/$shard_count: building shared targets ==="
        ./mfc.sh build -i "$first_case" -t syscheck pre_process post_process --case-optimization $gpu_opts -j 8
        touch "$shared_marker"
    else
        echo "=== Shard $shard_idx/$shard_count: waiting for shard 1 to build shared targets ==="
        waited=0
        until [ -f "$shared_marker" ]; do
            if [ "$waited" -ge 5400 ]; then
                echo "ERROR: timed out waiting for $shared_marker"; exit 1
            fi
            sleep 30
            waited=$((waited + 30))
        done
    fi
fi

idx=0
for case in benchmarks/*/case.py; do
    idx=$((idx + 1))
    if [ -n "$shard" ] && [ $(((idx - 1) % shard_count)) -ne $((shard_idx - 1)) ]; then
        continue
    fi
    echo "=== Pre-building: $case ==="
    ./mfc.sh run "$case" --case-optimization $gpu_opts -j 8 --dry-run
done
