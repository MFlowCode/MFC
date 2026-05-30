#!/bin/bash
set -e
NJOBS="${SLURM_CPUS_ON_NODE:-24}"; [ "$NJOBS" -gt 64 ] && NJOBS=64
./mfc.sh clean
source .github/scripts/retry-build.sh
retry_build ./mfc.sh build --gcov -j 8
./mfc.sh test --build-coverage-map --gcov -j "$NJOBS"
