#!/bin/bash
# Run a single benchmark on a Frontier compute node (build already done on login node).
# Usage: frontier_bench_config.sh <cluster> <device> <interface>
# Runs inside a SLURM allocation on an ssh'd compute node.

set -e
set -x

cluster=$1; device=$2; interface=$3

flag="f"; [ "$cluster" = "frontier_amd" ] && flag="famd"

. ./mfc.sh load -c "$flag" -m g

# Benchmark
job_slug="bench-${device}-${interface}"
n_ranks=12
device_opts=""
if [ "$device" = "gpu" ]; then
    gpus=$(rocm-smi --showid | awk '{print $1}' | grep -Eo '[0-9]+' | uniq | tr '\n' ' ')
    n_ranks=$(echo "$gpus" | wc -w)
    gpu_ids=$(echo "$gpus" | tr ' ' '\n' | tr '\n' ' ' | sed 's/ $//')
    device_opts="--gpu"
    [ "$interface" = "acc" ] && device_opts+=" acc"
    [ "$interface" = "omp" ] && device_opts+=" mp"
    device_opts+=" -g $gpu_ids"
fi

if [ "$device" = "gpu" ]; then
    ./mfc.sh bench --mem 12 -j $n_ranks -o "$job_slug.yaml" -- -c "$cluster" $device_opts -n $n_ranks
else
    ./mfc.sh bench --mem 1 -j $(nproc) -o "$job_slug.yaml" -- -c "$cluster" $device_opts -n $n_ranks
fi
