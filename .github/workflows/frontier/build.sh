#!/bin/bash

# Ignore SIGHUP to survive login node session drops
trap '' HUP

job_device=$1
job_interface=$2
run_bench=$3
build_opts=""
if [ "$job_device" = "gpu" ]; then
  build_opts+="--gpu"
  if [ "$job_interface" = "acc" ]; then
      build_opts+=" acc"
  elif [ "$job_interface" = "omp" ]; then
      build_opts+=" mp"
  fi
fi

. ./mfc.sh load -c f -m g

max_attempts=3
attempt=1
while [ $attempt -le $max_attempts ]; do
    echo "Build attempt $attempt of $max_attempts..."
    if [ "$run_bench" == "bench" ]; then
        build_cmd_ok=true
        for dir in benchmarks/*/; do
            dirname=$(basename "$dir")
            if ! ./mfc.sh run -v "$dir/case.py" --case-optimization -j 8 --dry-run $build_opts; then
                build_cmd_ok=false
                break
            fi
        done
    else
        if ./mfc.sh test -v -a --dry-run --rdma-mpi -j 8 $build_opts; then
            build_cmd_ok=true
        else
            build_cmd_ok=false
        fi
    fi

    if [ "$build_cmd_ok" = true ]; then
        echo "Build succeeded on attempt $attempt."
        exit 0
    fi

    if [ $attempt -lt $max_attempts ]; then
        echo "Build failed on attempt $attempt. Cleaning and retrying in 30s..."
        ./mfc.sh clean
        sleep 30
    fi
    attempt=$((attempt + 1))
done

echo "Build failed after $max_attempts attempts."
exit 1
