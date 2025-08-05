#!/bin/bash

gpus=`rocm-smi --showid | awk '{print $1}' | grep -Eo '[0-9]+' | uniq | tr '\n' ' '`
ngpus=`echo "$gpus" | tr -d '[:space:]' | wc -c`

if [ "$job_device" = "gpu" ]; then
    ./mfc.sh test -a --rdma-mpi --max-attempts 3 -j $ngpus -- -c frontier
else
    ./mfc.sh test -a --rdma-mpi --max-attempts 3 -j 32 -- -c frontier
fi
