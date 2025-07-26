#!/bin/bash

gpus=`rocm-smi --showid | awk '{print $1}' | grep -Eo '[0-9]+' | uniq | tr '\n' ' '`
ngpus=`echo "$gpus" | tr -d '[:space:]' | wc -c`

if [ "$job_device" = "gpu" ]; then
    ./mfc.sh test --max-attempts 3 -j $ngpus --sched-debug -- -c frontier
else
    ./mfc.sh test --max-attempts 3 -j 32 --sched-debug -- -c frontier
fi
