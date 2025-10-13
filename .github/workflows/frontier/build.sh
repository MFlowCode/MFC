#!/bin/bash

build_opts=""
if [ "$1" = "gpu" ]; then
    build_opts="--gpu"
fi

. ./mfc.sh load -c f -m g

if [ "$2" == "bench" ]; then
    for dir in benchmarks/*/; do
        dirname=$(basename "$dir")
        ./mfc.sh run "$dir/case.py" --case-optimization -j 8 --dry-run $build_opts
    done
else
    ./mfc.sh test -a --dry-run --rdma-mpi -j 8 $build_opts
fi

