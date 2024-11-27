#!/bin/bash

build_opts=""
if [ "$job_device" == "gpu" ]; then
    build_opts="--gpu"
fi

./mfc.sh test --dry-run -j 20 $build_opts

n_test_threads=8

if [ "$job_device" == "gpu" ]; then
    gpu_count=$(nvidia-smi -L | wc -l)        # number of GPUs on node
    gpu_ids=$(seq -s ' ' 0 $(($gpu_count-1))) # 0,1,2,...,gpu_count-1
    device_opts="-g $gpu_ids"
    n_test_threads=`expr $gpu_count \* 2`
fi

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/sw/spack/deltas11-2023-03/apps/linux-rhel8-zen3/nvhpc-24.1/openmpi-4.1.5-zkiklxi/lib/
./mfc.sh test --max-attempts 3 -a -j $n_test_threads $device_opts -- -c delta

