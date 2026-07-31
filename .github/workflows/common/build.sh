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

# --- Variant selection ---
# The suite needs two simulation binaries: the base one and a chemistry one (the
# mechanism is compiled in). On amdflang each device link is ~1 hour, so building
# both serially exceeds the 2 h walltime. $job_variant lets the caller split them
# into concurrent jobs that write to disjoint staging directories:
#   base -> plain targets, what every non-chemistry test runs
#   chem -> the chemistry variant only
# Unset builds everything in this job (default for every other cluster).
case "${job_variant:-}" in
    base) build_cmd=(./mfc.sh build -j 8 $build_opts) ;;
    chem) build_cmd=(./mfc.sh test -v --dry-run -a -j 8 -o Chemistry $build_opts) ;;
    "")   build_cmd=(./mfc.sh test -v --dry-run -a -j 8 $build_opts) ;;
    *)    echo "ERROR: unknown job_variant '$job_variant'"; exit 1 ;;
esac

RETRY_VALIDATE_CMD="$validate_cmd" \
    retry_build "${build_cmd[@]}" || exit 1
