#!/bin/bash

source .github/scripts/bench-preamble.sh

# Cap parallel jobs at 64 to avoid overwhelming MPI daemons on large nodes.
n_jobs=$(( $(nproc) > 64 ? 64 : $(nproc) ))

if [ "$job_device" = "gpu" ]; then
    ./mfc.sh bench --mem 4 -j $n_ranks -o "$job_slug.yaml" -- -c $job_cluster $device_opts -n $n_ranks
else
    ./mfc.sh bench --mem 1 -j $n_jobs -o "$job_slug.yaml" -- -c $job_cluster $device_opts -n $n_ranks
fi
