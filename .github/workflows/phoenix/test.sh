#!/bin/bash

build_opts=""
if [ "$job_device" = "gpu" ]; then
    build_opts="--gpu"
    if [ "$job_interface" = "omp" ]; then
      build_opts+=" mp"
    elif [ "$job_interface" = "acc" ]; then
      build_opts+=" acc"
    fi
fi

# Set up persistent build cache
source .github/scripts/setup-build-cache.sh phoenix "$job_device" "$job_interface"

max_attempts=3
attempt=1
while [ $attempt -le $max_attempts ]; do
    echo "Build attempt $attempt of $max_attempts..."
    if ./mfc.sh test -v --dry-run -j 8 $build_opts; then
        echo "Build succeeded on attempt $attempt."
        break
    fi

    if [ $attempt -lt $max_attempts ]; then
        echo "Build failed on attempt $attempt. Clearing staging/install and retrying in 30s..."
        rm -rf build/staging build/install build/lock.yaml
        sleep 30
    else
        echo "Build failed after $max_attempts attempts."
        exit 1
    fi
    attempt=$((attempt + 1))
done

# Release the cache lock before running tests. Tests only read installed
# binaries and can take hours — no need to block other builds.
if [ "${_cache_locked:-false}" = true ]; then
    flock --unlock 9
    exec 9>&-
    echo "Released build cache lock before tests"
fi

n_test_threads=8

if [ "$job_device" = "gpu" ]; then
    gpu_count=$(nvidia-smi -L | wc -l)        # number of GPUs on node
    gpu_ids=$(seq -s ' ' 0 $(($gpu_count-1))) # 0,1,2,...,gpu_count-1
    device_opts="-g $gpu_ids"
    n_test_threads=`expr $gpu_count \* 2`
fi

./mfc.sh test -v --max-attempts 3 -a -j $n_test_threads $device_opts -- -c phoenix

