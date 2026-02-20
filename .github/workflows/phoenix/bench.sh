#!/bin/bash

n_ranks=12

echo "My interface is:" $job_interface
build_opts=""
device_opts=""
if [ "$job_device" = "gpu" ]; then
    n_ranks=$(nvidia-smi -L | wc -l)        # number of GPUs on node
    gpu_ids=$(seq -s ' ' 0 $(($n_ranks-1))) # 0,1,2,...,gpu_count-1
    build_opts+="--gpu"
    if [ "$job_interface" = "acc" ]; then
        build_opts+=" acc"
    elif [ "$job_interface" = "omp" ]; then
        build_opts+=" mp"
    fi
    device_opts="$build_opts -g $gpu_ids"
fi

tmpbuild=/storage/scratch1/6/sbryngelson3/mytmp_build
currentdir=$tmpbuild/run-$(( RANDOM % 900 ))
mkdir -p $tmpbuild
mkdir -p $currentdir

export TMPDIR=$currentdir

if [ "$job_device" = "gpu" ]; then
    bench_opts="--mem 12"
else
    bench_opts="--mem 1"
fi

max_attempts=3
attempt=1
while [ $attempt -le $max_attempts ]; do
    echo "Build attempt $attempt of $max_attempts..."
    build_cmd_ok=true
    for dir in benchmarks/*/; do
        if ! ./mfc.sh run -v "$dir/case.py" --case-optimization -j $(nproc) --dry-run $build_opts; then
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

./mfc.sh bench $bench_opts -j $(nproc) -o "$job_slug.yaml" -- -c phoenix-bench $device_opts -n $n_ranks

sleep 10
rm -rf "$currentdir" || true

unset TMPDIR
