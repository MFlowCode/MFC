#!/bin/bash

n_ranks=1

if [ "$job_device" == "gpu" ]; then
    n_ranks=$(nvidia-smi -L | wc -l)        # number of GPUs on node
    gpu_ids=$(seq -s ' ' 0 $(($n_ranks-1))) # 0,1,2,...,gpu_count-1
    device_opts="--gpu -g $gpu_ids"
fi

# ./mfc.sh bench -j $(nproc) -o "$job_slug.yaml" -- -c phoenix $device_opts -n $n_ranks

gpu_count=$(nvidia-smi -L | wc -l)        # number of GPUs on node
gpu_ids=$(seq -s ' ' 0 $(($gpu_count-1))) # 0,1,2,...,gpu_count-1

# ./mfc.sh test -j 1 --gpu -g $gpu_ids -- -b mpirun

./mfc.sh bench -j $(nproc) -o "$job_slug.yaml" --gpu -- -b mpirun $device_opts -n $n_ranks
