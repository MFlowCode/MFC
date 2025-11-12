#!/usr/bin/env bash
# Simplified MFC Code Coverage Script
# More robust version with better error handling

set -euo pipefail

# Configuration
PERCENT=${PERCENT:-10}
JOBS=${JOBS:-$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4)}

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo -e "${GREEN}=== MFC Coverage (Simple Mode) ===${NC}"
echo "Running ${PERCENT}% of tests with ${JOBS} jobs"
echo ""

# Step 1: Clean
echo -e "${YELLOW}[1/5] Cleaning...${NC}"
./mfc.sh clean

# Step 2: Build with coverage
echo -e "${YELLOW}[2/5] Building with coverage...${NC}"
./mfc.sh build --gcov --no-gpu --debug -t pre_process simulation post_process -j ${JOBS}

# Step 3: Run tests with GCOV_PREFIX
echo -e "${YELLOW}[3/5] Running tests...${NC}"
export GCOV_PREFIX=${PWD}/build/staging
export GCOV_PREFIX_STRIP=0

./mfc.sh test --no-examples -% ${PERCENT} -j ${JOBS} || {
    echo -e "${RED}Some tests failed, but continuing with coverage analysis${NC}"
}

# Step 4: Find .gcda files
echo -e "${YELLOW}[4/5] Locating coverage data...${NC}"
GCDA_COUNT=$(find build -name "*.gcda" 2>/dev/null | wc -l | tr -d ' ')
GCNO_COUNT=$(find build -name "*.gcno" 2>/dev/null | wc -l | tr -d ' ')

echo "Found ${GCDA_COUNT} .gcda files"
echo "Found ${GCNO_COUNT} .gcno files"

if [ "$GCDA_COUNT" -eq 0 ]; then
    echo -e "${RED}Error: No .gcda files found. Tests may not have run with coverage.${NC}"
    exit 1
fi

# Step 5: Generate reports with simpler options
echo -e "${YELLOW}[5/5] Generating coverage reports (simplified)...${NC}"
mkdir -p build/coverage

# Find correct gcov
GCOV_EXEC=$(which gcov-15 2>/dev/null || which gcov-14 2>/dev/null || which gcov)
echo "Using: ${GCOV_EXEC}"

# Generate text summary first (most robust)
echo "Generating text summary..."
gcovr build/staging --root . \
    --gcov-executable "${GCOV_EXEC}" \
    --filter 'src/.*' \
    --print-summary 2>&1 | tee build/coverage/summary.txt || {
    echo -e "${YELLOW}Warning: Text summary generation had issues${NC}"
}

# Try HTML report (single threaded for stability)
echo "Generating HTML report..."
gcovr build/staging --root . \
    --gcov-executable "${GCOV_EXEC}" \
    --filter 'src/.*' \
    --html -o build/coverage/index.html \
    -j 1 || {
    echo -e "${YELLOW}Warning: HTML report generation failed${NC}"
}

# Try XML report
echo "Generating XML report..."
gcovr build/staging --root . \
    --gcov-executable "${GCOV_EXEC}" \
    --filter 'src/.*' \
    --xml -o build/coverage/coverage.xml \
    -j 1 || {
    echo -e "${YELLOW}Warning: XML report generation failed${NC}"
}

echo ""
echo -e "${GREEN}=== Coverage Complete ===${NC}"
echo ""
echo "View reports:"
echo "  open build/coverage/index.html"
echo "  cat build/coverage/summary.txt"
echo ""

# Show summary if available
if [ -f build/coverage/summary.txt ]; then
    echo "Summary:"
    cat build/coverage/summary.txt
fi






