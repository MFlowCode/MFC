#!/bin/bash

source .github/scripts/gpu-opts.sh
build_opts="$gpu_opts"

# Set up persistent build cache
source .github/scripts/setup-build-cache.sh phoenix "$job_device" "$job_interface"

# Build with retry; smoke-test cached binaries to catch architecture mismatches
# (SIGILL from binaries compiled on a different compute node).
source .github/scripts/retry-build.sh
RETRY_VALIDATE_CMD='syscheck_bin=$(find build/install -name syscheck -type f 2>/dev/null | head -1); [ -z "$syscheck_bin" ] || "$syscheck_bin" > /dev/null 2>&1' \
    retry_build ./mfc.sh test -v --dry-run -j 8 $build_opts || exit 1

n_test_threads=8

if [ "$job_device" = "gpu" ]; then
    source .github/scripts/detect-gpus.sh
    device_opts="-g $gpu_ids"
    n_test_threads=$((ngpus * 2))
fi

./mfc.sh test -v --max-attempts 3 -a -j $n_test_threads $device_opts -- -c phoenix
