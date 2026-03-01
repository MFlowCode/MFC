#!/bin/bash

# Number of parallel jobs: use SLURM allocation or default to 24.
# Cap at 64 to avoid overwhelming MPI's ORTE daemons with concurrent launches.
NJOBS="${SLURM_CPUS_ON_NODE:-24}"
if [ "$NJOBS" -gt 64 ]; then NJOBS=64; fi

# Build MFC with gcov coverage instrumentation (CPU-only, gfortran).
# -j 8 for compilation (memory-heavy, more cores doesn't help much).
./mfc.sh build --gcov -j 8

# Run all tests in parallel, collecting per-test coverage data.
# Each test gets an isolated GCOV_PREFIX directory so .gcda files
# don't collide. Coverage is collected per-test after all tests finish.
# --gcov is required so the internal build step preserves instrumentation.
./mfc.sh test --build-coverage-cache --gcov -j "$NJOBS"
