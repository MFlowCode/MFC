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
fi

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
