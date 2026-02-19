#!/bin/bash

set -e

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

. ./mfc.sh load -c famd -m g

# Clean stale build artifacts from previous CI runs
./mfc.sh clean

if [ "$run_bench" == "bench" ]; then
    for dir in benchmarks/*/; do
        ./mfc.sh run -v "$dir/case.py" --case-optimization -j 4 --dry-run $build_opts
    done
else
    ./mfc.sh test -v -a --dry-run -j 4 $build_opts
fi
