#!/bin/bash

source .github/scripts/gpu-opts.sh
source .github/scripts/detect-gpus.sh

n_ranks=12
echo "My interface is:" $job_interface
build_opts="$gpu_opts"
device_opts=""
if [ "$job_device" = "gpu" ]; then
    n_ranks=$ngpus
    device_opts="$gpu_opts -g $gpu_ids"
fi

tmpbuild=/storage/scratch1/6/sbryngelson3/mytmp_build
currentdir=$tmpbuild/run-$(( RANDOM % 900 ))
mkdir -p $tmpbuild
mkdir -p $currentdir

export TMPDIR=$currentdir

if [ "$job_device" = "gpu" ]; then
    bench_opts="--mem 4"
else
    bench_opts="--mem 1"
fi

source .github/scripts/retry-build.sh
RETRY_CLEAN_CMD="./mfc.sh clean" retry_build ./mfc.sh build -j $(nproc) $build_opts || exit 1

./mfc.sh bench $bench_opts -j $(nproc) -o "$job_slug.yaml" -- -c phoenix-bench $device_opts -n $n_ranks

sleep 10
rm -rf "$currentdir" || true

unset TMPDIR
