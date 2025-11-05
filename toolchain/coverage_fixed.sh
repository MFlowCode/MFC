#!/usr/bin/env bash
# MFC Coverage - Fixed for gcovr path issues
# Don't use GCOV_PREFIX - let gcovr find files naturally

set -euo pipefail

PERCENT=${PERCENT:-10}
JOBS=${JOBS:-$(sysctl -n hw.ncpu 2>/dev/null || echo 4)}

GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo -e "${GREEN}=== MFC Coverage (Fixed) ===${NC}"
echo "Running ${PERCENT}% of tests"
echo ""

# Step 1: Clean
echo -e "${YELLOW}[1/4] Cleaning...${NC}"
./mfc.sh clean

# Step 2: Build with coverage
echo -e "${YELLOW}[2/4] Building with coverage...${NC}"
./mfc.sh build --gcov --no-gpu --debug -t pre_process simulation post_process -j ${JOBS}

# Step 3: Run tests WITHOUT GCOV_PREFIX (let them write to build dirs naturally)
echo -e "${YELLOW}[3/4] Running ${PERCENT}% of tests...${NC}"
./mfc.sh test --no-examples --no-build -% ${PERCENT} -j ${JOBS} || {
    echo "Some tests failed, continuing..."
}

# Step 4: Generate reports - point gcovr at the build directories
echo -e "${YELLOW}[4/4] Generating reports...${NC}"
mkdir -p build/coverage

GCOV_EXEC=$(which gcov-15 2>/dev/null || which gcov-14 2>/dev/null || which gcov)
echo "Using: ${GCOV_EXEC}"

# Search all build directories for coverage data
BUILD_DIRS=$(find build/staging -type d -name "CMakeFiles" 2>/dev/null | sed 's|/CMakeFiles||' | head -10)

echo "Found build directories:"
echo "$BUILD_DIRS"
echo ""

# Try gcovr on each build directory
for BUILD_DIR in $BUILD_DIRS; do
    echo "Processing: $BUILD_DIR"
    
    gcovr "$BUILD_DIR" \
        --root . \
        --gcov-executable "${GCOV_EXEC}" \
        --filter 'src/.*' \
        -j 1 \
        --print-summary 2>&1 | tee -a build/coverage/summary.txt || true
done

echo ""
echo -e "${GREEN}=== Coverage Complete ===${NC}"
echo ""
echo "Summary saved to: build/coverage/summary.txt"
cat build/coverage/summary.txt







