#!/bin/bash

# Ignore SIGHUP to survive login node session drops
trap '' HUP

# Determine compiler flag from directory name
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cluster_name="$(basename "$SCRIPT_DIR")"
case "$cluster_name" in
    frontier)     compiler_flag="f" ;;
    frontier_amd) compiler_flag="famd" ;;
    *) echo "ERROR: Unknown cluster '$cluster_name'"; exit 1 ;;
esac

job_device=$1
job_interface=$2
source .github/scripts/gpu-opts.sh
build_opts="$gpu_opts"

. ./mfc.sh load -c $compiler_flag -m $([ "$job_device" = "gpu" ] && echo "g" || echo "c")

source .github/scripts/clean-build.sh
clean_build

source .github/scripts/retry-build.sh

# This login-node step is where Frontier's dependency install actually happens,
# and so where a PyPI outage actually lands -- 17 jobs spent ~33 minutes each
# rediscovering one on 2026-08-28. Tee the output and classify a failure so the
# first job to hit it records it and the rest of the matrix can skip.
# No set -e in this script, so capture the status rather than toggling it.
deps_log="deps-${cluster_name}-${job_device}-${job_interface}.log"
retry_build ./mfc.sh build --deps-only -j 8 $build_opts 2>&1 | tee "$deps_log"
deps_rc=${PIPESTATUS[0]}

if [ "$deps_rc" -ne 0 ]; then
    cls=0
    bash .github/scripts/classify-build-failure.sh "$deps_log" "$cluster_name" || cls=$?
    if [ "$cls" -ne 0 ]; then
        exit "$cls"
    fi
    exit 1
fi
