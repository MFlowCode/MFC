#!/bin/bash
# Shared preamble for benchmark scripts: detects GPUs, sets build/device opts.
# Sets: $gpu_opts, $build_opts, $device_opts, $n_ranks, $ngpus, $gpu_ids
# Usage: source .github/scripts/bench-preamble.sh

source .github/scripts/detect-gpus.sh
source .github/scripts/gpu-opts.sh

n_ranks=12
build_opts="$gpu_opts"
device_opts=""
if [ "$job_device" = "gpu" ]; then
    n_ranks=$ngpus
    device_opts="$gpu_opts -g $gpu_ids"
fi
