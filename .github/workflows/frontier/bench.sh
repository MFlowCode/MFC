#!/bin/bash

n_ranks=12
build_opts=""
device_opts=""
if [ "$job_device" = "gpu" ]; then
    gpus=$(rocm-smi --showid | awk '{print $1}' | grep -Eo '[0-9]+' | uniq | tr '\n' ' ')
    n_ranks=$(echo "$gpus" | wc -w)         # number of GPUs on node
    gpu_ids=$(echo "$gpus" | tr ' ' '\n' | tr '\n' ' ' | sed 's/ $//')  # GPU IDs from rocm-smi
    build_opts+="--gpu"
    if [ "$job_interface" = "acc" ]; then
        build_opts+=" acc"
    elif [ "$job_interface" = "omp" ]; then
        build_opts+=" mp"
    fi
    device_opts="$build_opts -g $gpu_ids"
fi

# Build case-optimized binaries on compute node (deps already fetched on login node)
max_attempts=3
attempt=1
while [ $attempt -le $max_attempts ]; do
    echo "Build attempt $attempt of $max_attempts..."
    build_cmd_ok=true
    for dir in benchmarks/*/; do
        if ! ./mfc.sh run -v "$dir/case.py" --case-optimization -j 8 --dry-run $build_opts; then
            build_cmd_ok=false
            break
        fi
    done

    if [ "$build_cmd_ok" = true ]; then
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

if [ "$job_device" = "gpu" ]; then
    ./mfc.sh bench --mem 12 -j $n_ranks -o "$job_slug.yaml" -- -c frontier $device_opts -n $n_ranks
else
    ./mfc.sh bench --mem 1 -j $(nproc) -o "$job_slug.yaml" -- -c frontier $device_opts -n $n_ranks
fi
