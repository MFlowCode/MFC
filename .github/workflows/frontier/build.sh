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
run_bench=$3
source .github/scripts/gpu-opts.sh
build_opts="$gpu_opts"

. ./mfc.sh load -c $compiler_flag -m $([ "$job_device" = "gpu" ] && echo "g" || echo "c")

source .github/scripts/clean-build.sh
clean_build

source .github/scripts/retry-build.sh
if [ "$run_bench" == "bench" ]; then
    retry_build ./mfc.sh build -j 8 $build_opts || exit 1
else
    retry_build ./mfc.sh test -v -a --dry-run $([ "$cluster_name" = "frontier" ] && echo "--rdma-mpi") -j 8 $build_opts || exit 1
fi
