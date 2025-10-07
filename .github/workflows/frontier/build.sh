#!/bin/bash

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

if [ "$run_bench" == "bench" ]; then
    for dir in benchmarks/*/; do
        dirname=$(basename "$dir")
        ./mfc.sh run "$dir/case.py" --case-optimization -j 8 --dry-run $build_opts
    done
else
    ./mfc.sh test -a --dry-run --rdma-mpi -j 8 $build_opts
fi

