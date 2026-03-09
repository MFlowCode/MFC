#!/bin/bash

source .github/scripts/bench-preamble.sh

# Cap parallel jobs at 64 to avoid overwhelming MPI daemons on large nodes
# (GNR nodes have 192 cores but nproc is too aggressive for build/bench).
n_jobs=$(( $(nproc) > 64 ? 64 : $(nproc) ))

tmpbuild=/storage/project/r-sbryngelson3-0/sbryngelson3/mytmp_build
currentdir=$tmpbuild/run-$(( RANDOM % 900 ))
mkdir -p $tmpbuild
mkdir -p $currentdir

export TMPDIR=$currentdir

if [ "$job_device" = "gpu" ]; then
    bench_opts="--mem 4"
else
    bench_opts="--mem 1"
fi

rm -rf build

source .github/scripts/retry-build.sh
retry_build ./mfc.sh build -j $n_jobs $build_opts || exit 1

./mfc.sh bench $bench_opts -j $n_jobs -o "$job_slug.yaml" -- -c phoenix-bench $device_opts -n $n_ranks

sleep 10
rm -rf "$currentdir" || true

unset TMPDIR
