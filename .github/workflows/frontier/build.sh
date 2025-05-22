#!/bin/bash

build_opts=""
if [ "$job_device" == "gpu" ]; then
    build_opts="--gpu"
fi

. ./mfc.sh load -c f -m g
./mfc.sh test --dry-run -j 8 $build_opts
