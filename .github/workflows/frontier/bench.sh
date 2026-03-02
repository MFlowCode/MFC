#!/bin/bash

source .github/scripts/bench-preamble.sh

if [ "$job_device" = "gpu" ]; then
    ./mfc.sh bench --mem 4 -j $n_ranks -o "$job_slug.yaml" -- -c $job_cluster $device_opts -n $n_ranks
else
    ./mfc.sh bench --mem 1 -j $(nproc) -o "$job_slug.yaml" -- -c $job_cluster $device_opts -n $n_ranks
fi
