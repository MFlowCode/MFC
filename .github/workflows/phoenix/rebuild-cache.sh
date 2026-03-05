#!/bin/bash
set -e

# Number of parallel jobs: use SLURM allocation or default to 24.
# Cap at 64 to avoid overwhelming OpenMPI daemons and OS process limits with concurrent launches.
NJOBS="${SLURM_CPUS_ON_NODE:-24}"
if [ "$NJOBS" -gt 64 ]; then NJOBS=64; fi

# Clean stale build artifacts: the self-hosted runner may have a cached
# GPU build (e.g. --gpu mp) whose CMake flags are incompatible with gcov.
./mfc.sh clean

# Source retry_build() for NFS stale file handle resilience (3 attempts).
source .github/scripts/retry-build.sh

# Build MFC with gcov coverage instrumentation (CPU-only, gfortran).
retry_build ./mfc.sh build --gcov -j 8

# Run all tests in parallel, collecting per-test coverage data.
# Each test gets an isolated GCOV_PREFIX directory so .gcda files
# don't collide. Coverage is collected per-test after all tests finish.
# --gcov is required so the internal build step preserves instrumentation.
./mfc.sh test --build-coverage-cache --gcov -j "$NJOBS"
