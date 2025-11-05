#!/usr/bin/env bash
# MFC Full Coverage Run - 100% of all tests
# Based on the working coverage_fixed.sh script

set -euo pipefail

PERCENT=100
JOBS=${JOBS:-$(sysctl -n hw.ncpu 2>/dev/null || echo 4)}

GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

LOGFILE="build/full_coverage_run.log"
mkdir -p build

# Redirect all output to both console and log file
exec > >(tee "$LOGFILE") 2>&1

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}MFC FULL Coverage Run (100% of tests)${NC}"
echo -e "${GREEN}========================================${NC}"
echo "Started: $(date)"
echo "Jobs: ${JOBS}"
echo ""

# Step 1: Clean
echo -e "${YELLOW}[1/4] Cleaning previous builds...${NC}"
./mfc.sh clean
echo "Clean complete at: $(date)"
echo ""

# Step 2: Build with coverage
echo -e "${YELLOW}[2/4] Building with coverage instrumentation...${NC}"
echo "This will take 10-15 minutes..."
START_BUILD=$(date +%s)
./mfc.sh build --gcov --no-gpu --debug -t pre_process simulation post_process -j ${JOBS}
END_BUILD=$(date +%s)
BUILD_TIME=$((END_BUILD - START_BUILD))
echo "Build complete at: $(date)"
echo "Build time: ${BUILD_TIME} seconds"
echo ""

# Step 3: Run ALL tests
echo -e "${YELLOW}[3/4] Running 100% of test suite...${NC}"
echo "This will take 2-4 hours..."
echo "Test run started at: $(date)"
START_TEST=$(date +%s)

./mfc.sh test --no-examples --no-build -j ${JOBS} || {
    echo -e "${YELLOW}WARNING: Some tests failed, continuing with coverage analysis...${NC}"
}

END_TEST=$(date +%s)
TEST_TIME=$((END_TEST - START_TEST))
echo "Test run complete at: $(date)"
echo "Test time: ${TEST_TIME} seconds ($((TEST_TIME/60)) minutes)"
echo ""

# Step 4: Generate coverage reports
echo -e "${YELLOW}[4/4] Generating coverage reports...${NC}"
mkdir -p build/coverage_full

GCOV_EXEC=$(which gcov-15 2>/dev/null || which gcov-14 2>/dev/null || which gcov)
echo "Using gcov: ${GCOV_EXEC}"
echo ""

# Find .gcda files
GCDA_COUNT=$(find build/staging -name "*.gcda" 2>/dev/null | wc -l)
GCNO_COUNT=$(find build/staging -name "*.gcno" 2>/dev/null | wc -l)

echo "Coverage data files found:"
echo "  .gcda files: ${GCDA_COUNT}"
echo "  .gcno files: ${GCNO_COUNT}"
echo ""

if [ "${GCDA_COUNT}" -eq 0 ]; then
    echo -e "${RED}ERROR: No .gcda files found!${NC}"
    echo "Coverage data was not collected."
    exit 1
fi

# Search all build directories for coverage data
echo "Searching for build directories with coverage data..."
BUILD_DIRS=$(find build/staging -type d -name "CMakeFiles" 2>/dev/null | sed 's|/CMakeFiles||' | head -20)

echo "Found $(echo "$BUILD_DIRS" | wc -l) build directories"
echo ""

# Generate report for each build directory
for BUILD_DIR in $BUILD_DIRS; do
    echo "Processing: $BUILD_DIR"
    
    gcovr "$BUILD_DIR" \
        --root . \
        --gcov-executable "${GCOV_EXEC}" \
        --filter 'src/.*' \
        -j 1 \
        --gcov-ignore-parse-errors=suspicious_hits.warn \
        --print-summary 2>&1 | tee -a build/coverage_full/summary.txt || {
        echo "  (had issues, continuing...)"
    }
    echo ""
done

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}FULL COVERAGE RUN COMPLETE${NC}"
echo -e "${GREEN}========================================${NC}"
echo "Completed: $(date)"
echo "Total time: $((BUILD_TIME + TEST_TIME)) seconds ($((( BUILD_TIME + TEST_TIME)/60)) minutes)"
echo ""
echo "Results saved to:"
echo "  build/coverage_full/summary.txt"
echo "  ${LOGFILE}"
echo ""
echo "=== COVERAGE SUMMARY ==="
cat build/coverage_full/summary.txt 2>/dev/null || echo "No summary generated"
echo ""





