#!/bin/bash

build_opts=""
if [ "$1" == "gpu" ]; then
    build_opts="--gpu"
fi

. ./mfc.sh load -c f -m g
./mfc.sh test --dry-run -j 8 $build_opts
