#!/bin/bash
# Build-only script for all clusters.
# Runs inside a SLURM job via submit-slurm-job.sh.
# Builds MFC without running tests (--dry-run).
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

# --- Build ---
# Phoenix builds everything inside SLURM (no login-node build step).
# Frontier/Frontier AMD: deps already fetched on login node via --deps-only;
# source code is built here on the compute node.
# Phoenix: always start fresh to avoid SIGILL from stale binaries compiled
# on a different microarchitecture.
if [ "$job_cluster" = "phoenix" ]; then
    source .github/scripts/clean-build.sh
    clean_build
fi

source .github/scripts/retry-build.sh

# Phoenix: smoke-test the syscheck binary to catch architecture mismatches
# (SIGILL from binaries compiled on a different compute node).
validate_cmd=""
if [ "$job_cluster" = "phoenix" ]; then
    validate_cmd='syscheck_bin=$(find build/install -name syscheck -type f 2>/dev/null | head -1); [ -z "$syscheck_bin" ] || "$syscheck_bin" > /dev/null 2>&1'
fi

RETRY_VALIDATE_CMD="$validate_cmd" \
    retry_build ./mfc.sh test -v --dry-run -a -j 8 $build_opts || exit 1
