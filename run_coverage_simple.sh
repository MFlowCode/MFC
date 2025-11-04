#!/usr/bin/env bash
# Simplified Full Coverage Comparison
# This script runs the full test suite with coverage and reports the results

set -euo pipefail

LOGDIR="build/coverage_comparison"
mkdir -p "$LOGDIR"

# Write immediately to file
exec > >(tee -a "$LOGDIR/run.log") 2>&1

echo "=========================================="
echo "MFC Full Coverage Comparison - STARTED"
echo "Time: $(date)"
echo "=========================================="

# Phase 1: Clean and build
echo ""
echo "[1/4] Cleaning previous builds..."
./mfc.sh clean

echo ""
echo "[2/4] Building with coverage instrumentation..."
echo "This may take 10-15 minutes..."
./mfc.sh build --gcov --no-gpu --debug -t pre_process simulation post_process -j $(sysctl -n hw.ncpu)

# Phase 2: Run tests
echo ""
echo "[3/4] Running FULL test suite (100% of all tests)..."
echo "This will take 2-4 hours..."
echo "Test run started at: $(date)"

GCOV_EXEC=$(which gcov-15 2>/dev/null || which gcov-14 2>/dev/null || which gcov)
export GCOV_PREFIX=${PWD}/build/staging
export GCOV_PREFIX_STRIP=0

time ./mfc.sh test --no-examples -j $(sysctl -n hw.ncpu) || {
    echo "WARNING: Some tests failed, continuing with coverage analysis..."
}

echo "Test run completed at: $(date)"

# Phase 3: Generate report
echo ""
echo "[4/4] Generating coverage report..."
mkdir -p "$LOGDIR/results"

gcovr build/staging --root . \
    --gcov-executable "${GCOV_EXEC}" \
    --filter 'src/.*' \
    -j 1 \
    --gcov-ignore-parse-errors=suspicious_hits.warn \
    --print-summary | tee "$LOGDIR/results/summary.txt"

echo ""
echo "=========================================="
echo "COVERAGE RUN COMPLETE"
echo "Time: $(date)"
echo "=========================================="
echo ""
echo "Summary saved to: $LOGDIR/results/summary.txt"
echo ""
echo "Coverage results:"
cat "$LOGDIR/results/summary.txt"





