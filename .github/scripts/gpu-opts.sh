#!/bin/bash
# Sets $gpu_opts from $job_device and $job_interface.
# Usage: source .github/scripts/gpu-opts.sh

gpu_opts=""
if [ "$job_device" = "gpu" ]; then
    gpu_opts="--gpu"
    if [ "$job_interface" = "omp" ]; then
        gpu_opts+=" mp"
    elif [ "$job_interface" = "acc" ]; then
        gpu_opts+=" acc"
    fi
fi
