#!/bin/bash
# Run a single test on a Frontier compute node (build already done on login node).
# Usage: frontier_test_config.sh <cluster> <device> <interface>
# Runs inside a SLURM allocation on an ssh'd compute node.

set -e
set -x

cluster=$1; device=$2; interface=$3

flag="f"; [ "$cluster" = "frontier_amd" ] && flag="famd"
mode="g"; [ "$device" = "cpu" ] && mode="c"

. ./mfc.sh load -c "$flag" -m "$mode"

# Device options
device_opts=""
if [ "$device" = "gpu" ]; then
    device_opts="--gpu"
    [ "$interface" = "acc" ] && device_opts+=" acc"
    [ "$interface" = "omp" ] && device_opts+=" mp"
fi

rdma=""
[ "$cluster" = "frontier" ] && [ "$device" = "gpu" ] && rdma="--rdma-mpi"

# Test
gpus=$(rocm-smi --showid | awk '{print $1}' | grep -Eo '[0-9]+' | uniq | tr '\n' ' ')
ngpus=$(echo "$gpus" | wc -w)

if [ "$device" = "gpu" ]; then
    ./mfc.sh test -v -a $rdma --max-attempts 3 -j $ngpus $device_opts -- -c "$cluster"
else
    ./mfc.sh test -v -a --max-attempts 3 -j 32 --no-gpu -- -c "$cluster"
fi
