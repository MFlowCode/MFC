#!/bin/bash

build_opts=""
if [ "$job_device" = "gpu" ]; then
    build_opts="--gpu"
fi

./mfc.sh test --dry-run -j 8 $build_opts

n_test_threads=8

if [ "$job_device" = "gpu" ]; then
    gpu_count=$(nvidia-smi -L | wc -l)        # number of GPUs on node
    gpu_ids=$(seq -s ' ' 0 $(($gpu_count-1))) # 0,1,2,...,gpu_count-1
    device_opts="-g $gpu_ids"
    n_test_threads=`expr $gpu_count \* 2`
fi

./mfc.sh test --max-attempts 3 -a -j $n_test_threads $device_opts -- -c phoenix

