#!/bin/bash


if [ -z "$job_device" ] || [ -z "$job_precision" ]; then
    echo "Usage: $0 [cpu|gpu] [single|double]"
    exit 1
fi

n_ranks=12
precision_flag=""

if [ "$job_precision" == "single" ]; then
    precision_flag="--single"
fi

if [ "$job_device" == "gpu" ]; then
    n_ranks=$(nvidia-smi -L | wc -l)
    gpu_ids=$(seq -s ' ' 0 $(($n_ranks-1)))
    device_opts="--gpu -g $gpu_ids"
else
    device_opts=""
fi

mem_value=1
if [ "$job_device" == "gpu" ]; then
    mem_value=12
fi

./mfc.sh bench --mem $mem_value -j $(nproc) -o "bench-${job_device}-${job_precision}.yaml" -- $precision_flag -c phoenix $device_opts -n $n_ranks

