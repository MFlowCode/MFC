#!/bin/bash

source .github/scripts/detect-gpus.sh
source .github/scripts/gpu-opts.sh
device_opts="$gpu_opts"

shard_opts=""
if [ -n "$job_shard" ]; then
    shard_opts="--shard $job_shard"
fi

if [ "$job_device" = "gpu" ]; then
    rdma_opts=""
    if [ "$job_cluster" = "frontier" ]; then
        rdma_opts="--rdma-mpi"
    fi
    ./mfc.sh test -v -a $rdma_opts --max-attempts 3 -j $ngpus $device_opts $shard_opts -- -c $job_cluster
else
    ./mfc.sh test -v -a --max-attempts 3 -j 32 --no-gpu $shard_opts -- -c $job_cluster
fi
