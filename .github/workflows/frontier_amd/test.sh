#!/bin/bash

gpus=`rocm-smi --showid | awk '{print $1}' | grep -Eo '[0-9]+' | uniq | tr '\n' ' '`
ngpus=`echo "$gpus" | tr -d '[:space:]' | wc -c`

device_opts=""
if [ "$job_device" = "gpu" ]; then
    device_opts+="--gpu"
    if [ "$job_interface" = "acc" ]; then
        device_opts+=" acc"
    elif [ "$job_interface" = "omp" ]; then
        device_opts+=" mp"
    fi
else
    device_opts+=" --no-gpu"
fi

# Build source code on compute node (deps already fetched on login node)
max_attempts=3
attempt=1
while [ $attempt -le $max_attempts ]; do
    echo "Build attempt $attempt of $max_attempts..."
    if ./mfc.sh test -v -a --dry-run -j 8 $device_opts; then
        echo "Build succeeded on attempt $attempt."
        break
    fi

    if [ $attempt -lt $max_attempts ]; then
        echo "Build failed on attempt $attempt. Cleaning and retrying in 30s..."
        ./mfc.sh clean
        sleep 30
    else
        echo "Build failed after $max_attempts attempts."
        exit 1
    fi
    attempt=$((attempt + 1))
done

# Run tests
if [ "$job_device" = "gpu" ]; then
    ./mfc.sh test -v -a --max-attempts 3 -j $ngpus $device_opts -- -c frontier_amd
else
    ./mfc.sh test -v -a --max-attempts 3 -j 32 --no-gpu -- -c frontier_amd
fi
