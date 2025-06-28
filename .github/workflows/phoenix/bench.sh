#!/bin/bash

n_ranks=12

if [ "$job_device" = "gpu" ]; then
    n_ranks=$(nvidia-smi -L | wc -l)        # number of GPUs on node
    gpu_ids=$(seq -s ' ' 0 $(($n_ranks-1))) # 0,1,2,...,gpu_count-1
    device_opts="--gpu -g $gpu_ids"
fi

tmpbuild=/storage/scratch1/6/sbryngelson3/mytmp_build
currentdir=$tmpbuild/run-$(( RANDOM % 900 ))
mkdir -p $tmpbuild
mkdir -p $currentdir

export TMPDIR=$currentdir

if [ "$job_device" = "gpu" ]; then
    ./mfc.sh bench --mem 12 -j $(nproc) -o "$job_slug.yaml" -- -c phoenix-bench $device_opts -n $n_ranks
else
    ./mfc.sh bench --mem 1 -j $(nproc) -o "$job_slug.yaml" -- -c phoenix-bench $device_opts -n $n_ranks
fi

sleep 10
rm -rf "$currentdir" || true

unset TMPDIR
