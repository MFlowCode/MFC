#!/bin/bash

source .github/scripts/detect-gpus.sh
source .github/scripts/gpu-opts.sh

n_ranks=12
device_opts=""
if [ "$job_device" = "gpu" ]; then
    n_ranks=$ngpus
    device_opts="$gpu_opts -g $gpu_ids"
fi

if [ "$job_device" = "gpu" ]; then
    ./mfc.sh bench --mem 4 -j $n_ranks -o "$job_slug.yaml" -- -c $job_cluster $device_opts -n $n_ranks
else
    ./mfc.sh bench --mem 1 -j $(nproc) -o "$job_slug.yaml" -- -c $job_cluster $device_opts -n $n_ranks
fi
