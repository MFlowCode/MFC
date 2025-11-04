#!/usr/bin/env bash
# Run baseline coverage with ORIGINAL test suite (before expansion)

set -euo pipefail

LOGDIR="build/coverage_baseline"
mkdir -p "$LOGDIR"

echo "=========================================="
echo "MFC Baseline Coverage (Original 790 Tests)"
echo "=========================================="
echo "Started: $(date)"
echo ""

# Step 1: Restore original cases.py
echo "[1/5] Restoring original test suite..."
git checkout HEAD -- toolchain/mfc/test/cases.py
echo "Original test suite restored"
echo ""

# Step 2: Clean
echo "[2/5] Cleaning previous builds..."
./mfc.sh clean

# Step 3: Build with coverage
echo "[3/5] Building with coverage instrumentation..."
START_BUILD=$(date +%s)
./mfc.sh build --gcov --no-gpu --debug -t pre_process simulation post_process -j $(sysctl -n hw.ncpu)
END_BUILD=$(date +%s)
BUILD_TIME=$((END_BUILD - START_BUILD))
echo "Build time: ${BUILD_TIME} seconds"
echo ""

# Step 4: Run ALL tests (100%) WITH post-processing
echo "[4/5] Running 100% of ORIGINAL test suite (WITH post-processing)..."
echo "Started at: $(date)"
START_TEST=$(date +%s)

# Run without --no-examples to enable post-processing
./mfc.sh test -j $(sysctl -n hw.ncpu) || {
    echo "WARNING: Some tests failed, continuing..."
}

END_TEST=$(date +%s)
TEST_TIME=$((END_TEST - START_TEST))
echo "Test time: ${TEST_TIME} seconds ($((TEST_TIME/60)) minutes)"
echo ""

# Step 5: Generate coverage reports
echo "[5/5] Generating coverage reports..."
mkdir -p "$LOGDIR"

GCOV_EXEC=$(which gcov-15 2>/dev/null || which gcov-14 2>/dev/null || which gcov)
echo "Using: ${GCOV_EXEC}"

# Find build directories
BUILD_DIRS=$(find build/staging -type d -name "CMakeFiles" 2>/dev/null | sed 's|/CMakeFiles||' | head -20)

# Generate reports for each build directory
for BUILD_DIR in $BUILD_DIRS; do
    echo "Processing: $BUILD_DIR"
    
    gcovr "$BUILD_DIR" \
        --root . \
        --gcov-executable "${GCOV_EXEC}" \
        --filter 'src/.*' \
        -j 1 \
        --gcov-ignore-parse-errors=suspicious_hits.warn \
        --print-summary 2>&1 | tee -a "$LOGDIR/summary.txt" || true
done

echo ""
echo "=========================================="
echo "BASELINE COVERAGE COMPLETE"
echo "=========================================="
echo "Completed: $(date)"
echo "Build time: ${BUILD_TIME} seconds"
echo "Test time: ${TEST_TIME} seconds ($((TEST_TIME/60)) minutes)"
echo "Total time: $((BUILD_TIME + TEST_TIME)) seconds ($(((BUILD_TIME + TEST_TIME)/60)) minutes)"
echo ""
echo "Summary saved to: $LOGDIR/summary.txt"
echo ""
cat "$LOGDIR/summary.txt" | grep -A 5 "TOTAL" || true

