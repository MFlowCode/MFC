#!/bin/bash

n_ranks=12

if [ "$job_device" = "gpu" ]; then
    gpus=$(rocm-smi --showid | awk '{print $1}' | grep -Eo '[0-9]+' | uniq | tr '\n' ' ')
    n_ranks=$(echo "$gpus" | wc -w)         # number of GPUs on node
    gpu_ids=$(echo "$gpus" | tr ' ' '\n' | tr '\n' ' ' | sed 's/ $//')  # GPU IDs from rocm-smi
    device_opts="--gpu -g $gpu_ids"
fi

if [ "$job_device" = "gpu" ]; then
    ./mfc.sh bench --mem 12 -j $n_ranks -o "$job_slug.yaml" -- -c frontier $device_opts -n $n_ranks
else
    ./mfc.sh bench --mem 1 -j $(nproc) -o "$job_slug.yaml" -- -c frontier $device_opts -n $n_ranks
fi
