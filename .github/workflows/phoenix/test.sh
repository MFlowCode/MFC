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

        # Smoke-test the cached binaries to catch architecture mismatches
        # (SIGILL from binaries compiled on a different compute node).
        syscheck_bin=$(find build/install -name syscheck -type f 2>/dev/null | head -1)
        if [ -n "$syscheck_bin" ] && ! "$syscheck_bin" > /dev/null 2>&1; then
            echo "WARNING: syscheck binary crashed — cached install is stale."
            if [ $attempt -lt $max_attempts ]; then
                echo "Clearing cache and rebuilding..."
                rm -rf build/staging build/install build/lock.yaml
                sleep 5
                attempt=$((attempt + 1))
                continue
            else
                echo "ERROR: syscheck still failing after $max_attempts attempts."
                exit 1
            fi
        fi

        break
    fi

    if [ $attempt -lt $max_attempts ]; then
        echo "Build failed on attempt $attempt. Clearing cache and retrying in 30s..."
        rm -rf build/staging build/install build/lock.yaml
        sleep 30
    else
        echo "Build failed after $max_attempts attempts."
        exit 1
    fi
    attempt=$((attempt + 1))
done

# Use up to 64 parallel test threads on CPU (GNR nodes have 192 cores).
# Cap at 64 to avoid overwhelming MPI's ORTE daemons with concurrent launches.
n_test_threads=$(( SLURM_CPUS_ON_NODE > 64 ? 64 : ${SLURM_CPUS_ON_NODE:-8} ))

if [ "$job_device" = "gpu" ]; then
    gpu_count=$(nvidia-smi -L | wc -l)        # number of GPUs on node
    gpu_ids=$(seq -s ' ' 0 $(($gpu_count-1))) # 0,1,2,...,gpu_count-1
    device_opts="-g $gpu_ids"
    n_test_threads=`expr $gpu_count \* 2`
fi

./mfc.sh test -v --max-attempts 3 -a -j $n_test_threads $device_opts -- -c phoenix
