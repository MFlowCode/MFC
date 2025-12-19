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

if [ "$job_device" = "gpu" ]; then
    ./mfc.sh test -a --rdma-mpi --max-attempts 3 -j $ngpus $device_opts -- -c frontier
else
    ./mfc.sh test -a --max-attempts 3 -j 32 --no-gpu -- -c frontier
fi
