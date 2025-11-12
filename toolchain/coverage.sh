#!/usr/bin/env bash
# MFC Code Coverage Assessment Script
# This script builds MFC with coverage instrumentation, runs tests, and generates reports

set -euo pipefail

# Configuration
PERCENT=${PERCENT:-25}
MIN_LINES=${MIN_LINES:-65}
MIN_BRANCHES=${MIN_BRANCHES:-50}
JOBS=${JOBS:-$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4)}

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}===================================${NC}"
echo -e "${GREEN}MFC Coverage Assessment${NC}"
echo -e "${GREEN}===================================${NC}"
echo ""

# Step 1: Clean previous builds
echo -e "${YELLOW}[1/6] Cleaning previous builds...${NC}"
./mfc.sh clean

# Step 2: Build with coverage instrumentation
echo -e "${YELLOW}[2/6] Building MFC with coverage (--gcov --no-gpu --debug)...${NC}"
./mfc.sh build --gcov --no-gpu --debug -t pre_process simulation post_process -j ${JOBS}

# Step 3: Set GCOV_PREFIX to ensure coverage data is written to build directory
# This is critical - without this, installed binaries won't write .gcda files back
export GCOV_PREFIX=${PWD}/build/staging
export GCOV_PREFIX_STRIP=0

echo -e "${YELLOW}[3/6] Running ${PERCENT}% of test suite (no --no-build flag)...${NC}"
echo "Note: Tests will run slower because each test rebuilds, but this ensures coverage data is collected."

# Run tests without --no-build to ensure coverage data is collected
# We need to temporarily modify how tests run to collect coverage properly
./mfc.sh test --no-examples -% ${PERCENT} -j ${JOBS} || {
    echo -e "${RED}Warning: Some tests may have failed, but continuing with coverage analysis${NC}"
}

# Step 4: Generate coverage reports
echo -e "${YELLOW}[4/6] Generating coverage reports...${NC}"
mkdir -p build/coverage

# Find the correct gcov executable that matches gfortran
GCOV_EXEC=$(which gcov-15 2>/dev/null || which gcov-14 2>/dev/null || which gcov)
echo "Using gcov executable: ${GCOV_EXEC}"

# Generate coverage reports
echo "Attempting coverage analysis..."
gcovr build/staging --root . \
    --gcov-executable "${GCOV_EXEC}" \
    --filter 'src/.*' \
    --html --html-details -o build/coverage/index.html \
    --xml-pretty -o build/coverage/coverage.xml \
    --txt -o build/coverage/summary.txt \
    --print-summary || {
    echo -e "${RED}Coverage generation failed. This may be due to:${NC}"
    echo "  1. No tests were executed successfully"
    echo "  2. Coverage data files (.gcda) were not written"
    echo "  3. Mismatch between .gcno and .gcda files"
    echo ""
    echo "Checking for coverage data files..."
    echo "Number of .gcda files: $(find build/staging -name '*.gcda' | wc -l)"
    echo "Number of .gcno files: $(find build/staging -name '*.gcno' | wc -l)"
    
    # Try to find at least one .gcda file and process it directly
    SAMPLE_GCDA=$(find build/staging -name '*.gcda' | head -1)
    if [ -n "$SAMPLE_GCDA" ]; then
        echo ""
        echo "Sample .gcda file: $SAMPLE_GCDA"
        echo "Attempting direct gcov on sample file..."
        GCDA_DIR=$(dirname "$SAMPLE_GCDA")
        (cd "$GCDA_DIR" && ${GCOV_EXEC} -o . *.gcda 2>&1 | head -20)
    fi
    
    exit 1
}

# Step 5: Display summary
echo ""
echo -e "${YELLOW}[5/6] Coverage Summary:${NC}"
cat build/coverage/summary.txt || echo "Summary file not generated"

# Step 6: Check thresholds
echo ""
echo -e "${YELLOW}[6/6] Checking coverage thresholds...${NC}"
echo "Minimum lines: ${MIN_LINES}%"
echo "Minimum branches: ${MIN_BRANCHES}%"

# Try to apply thresholds
gcovr build/staging --root . \
    --gcov-executable "${GCOV_EXEC}" \
    --filter 'src/.*' \
    --fail-under-line ${MIN_LINES} \
    --fail-under-branch ${MIN_BRANCHES} 2>/dev/null || {
    echo -e "${RED}Coverage below thresholds!${NC}"
    echo "To improve coverage:"
    echo "  1. Add unit tests for untested functions"
    echo "  2. Expand regression tests to cover more code paths"
    echo "  3. Review build/coverage/index.html for details"
}

echo ""
echo -e "${GREEN}===================================${NC}"
echo -e "${GREEN}Coverage report generated!${NC}"
echo -e "${GREEN}===================================${NC}"
echo ""
echo "View detailed HTML report:"
echo "  open build/coverage/index.html"
echo ""
echo "Or view text summary:"
echo "  cat build/coverage/summary.txt"
echo ""

