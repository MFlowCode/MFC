#!/usr/bin/env bash
# Run coverage with post-process validation (-a flag)
# This will test the post_process binary and collect its coverage

set -euo pipefail

RESULTS_DIR="coverage_results_postprocess"
mkdir -p "$RESULTS_DIR"

JOBS=$(sysctl -n hw.ncpu 2>/dev/null || echo 4)
GCOV_EXEC=$(which gcov-15 2>/dev/null || which gcov-14 2>/dev/null || which gcov)

echo "========================================" | tee "$RESULTS_DIR/progress.log"
echo "MFC Post-Process Coverage Run" | tee -a "$RESULTS_DIR/progress.log"
echo "Started: $(date)" | tee -a "$RESULTS_DIR/progress.log"
echo "========================================" | tee -a "$RESULTS_DIR/progress.log"
echo "" | tee -a "$RESULTS_DIR/progress.log"

# Clean
echo "[1/5] Cleaning..." | tee -a "$RESULTS_DIR/progress.log"
./mfc.sh clean >> "$RESULTS_DIR/progress.log" 2>&1

# Build
echo "[2/5] Building with coverage instrumentation..." | tee -a "$RESULTS_DIR/progress.log"
BUILD_START=$(date +%s)
./mfc.sh build --gcov --no-gpu --debug -t pre_process simulation post_process -j ${JOBS} >> "$RESULTS_DIR/build.log" 2>&1
BUILD_END=$(date +%s)
BUILD_TIME=$((BUILD_END - BUILD_START))
echo "  Build time: ${BUILD_TIME}s" | tee -a "$RESULTS_DIR/progress.log"

# Count tests
TEST_COUNT=$(./mfc.sh test --list 2>&1 | grep -E "^ *[A-F0-9]{8} " | wc -l | tr -d ' ')
echo "  Test count: ${TEST_COUNT}" | tee -a "$RESULTS_DIR/progress.log"

# Run tests WITH -a flag for post-processing validation
echo "[3/5] Running tests with POST-PROCESSING validation (-a flag)..." | tee -a "$RESULTS_DIR/progress.log"
echo "  This tests the post_process binary on all test outputs" | tee -a "$RESULTS_DIR/progress.log"
TEST_START=$(date +%s)
./mfc.sh test -a --no-build -j ${JOBS} >> "$RESULTS_DIR/tests.log" 2>&1 || {
    echo "  Some tests failed, continuing..." | tee -a "$RESULTS_DIR/progress.log"
}
TEST_END=$(date +%s)
TEST_TIME=$((TEST_END - TEST_START))
echo "  Test time: ${TEST_TIME}s ($((TEST_TIME / 60))m)" | tee -a "$RESULTS_DIR/progress.log"

# Generate coverage
echo "[4/5] Generating coverage report..." | tee -a "$RESULTS_DIR/progress.log"
gcovr build/staging \
    --root . \
    --gcov-executable "${GCOV_EXEC}" \
    --filter 'src/.*' \
    --gcov-ignore-parse-errors=suspicious_hits.warn \
    --print-summary \
    --txt -o "$RESULTS_DIR/coverage.txt" \
    --html --html-details -o "$RESULTS_DIR/index.html" \
    -j 1 2>&1 | tee -a "$RESULTS_DIR/progress.log"

# Extract summary
echo "" | tee -a "$RESULTS_DIR/progress.log"
echo "[5/5] Summary" | tee -a "$RESULTS_DIR/progress.log"
echo "========================================" | tee -a "$RESULTS_DIR/progress.log"
tail -20 "$RESULTS_DIR/coverage.txt" | tee -a "$RESULTS_DIR/progress.log"
echo "========================================" | tee -a "$RESULTS_DIR/progress.log"
echo "" | tee -a "$RESULTS_DIR/progress.log"
echo "Completed: $(date)" | tee -a "$RESULTS_DIR/progress.log"
echo "" | tee -a "$RESULTS_DIR/progress.log"
echo "Results in: $RESULTS_DIR/" | tee -a "$RESULTS_DIR/progress.log"
echo "  - coverage.txt (text report)" | tee -a "$RESULTS_DIR/progress.log"
echo "  - index.html (HTML report)" | tee -a "$RESULTS_DIR/progress.log"
echo "  - tests.log (test output)" | tee -a "$RESULTS_DIR/progress.log"
echo "  - build.log (build output)" | tee -a "$RESULTS_DIR/progress.log"

# Create comparison with baseline
if [ -f coverage_results/baseline_coverage.txt ]; then
    echo "" | tee -a "$RESULTS_DIR/progress.log"
    echo "Comparing with baseline..." | tee -a "$RESULTS_DIR/progress.log"
    
    BASELINE_LINES=$(grep "^lines:" coverage_results/baseline_coverage.txt | awk '{print $2}')
    POSTPROC_LINES=$(grep "^lines:" "$RESULTS_DIR/coverage.txt" | awk '{print $2}')
    
    cat > "$RESULTS_DIR/COMPARISON.md" << EOF
# Coverage Comparison: Baseline vs. Post-Processing

## Baseline Run (no -a flag)
- Tests: 528
- Coverage from: baseline_coverage.txt
- Line coverage: ${BASELINE_LINES}

## Post-Processing Run (with -a flag)
- Tests: ${TEST_COUNT}
- Coverage from: coverage.txt
- Line coverage: ${POSTPROC_LINES}

## Post-Process Module Coverage

### Baseline:
\`\`\`
$(grep -A 5 "src/post_process" coverage_results/baseline_coverage.txt | head -10)
\`\`\`

### With -a flag:
\`\`\`
$(grep -A 5 "src/post_process" "$RESULTS_DIR/coverage.txt" | head -10)
\`\`\`

## Full Comparison

View detailed reports:
- Baseline: coverage_results/baseline_coverage.txt
- Post-process: $RESULTS_DIR/coverage.txt

Open HTML reports:
- open coverage_results/baseline_coverage.html (if exists)
- open $RESULTS_DIR/index.html
EOF
    
    cat "$RESULTS_DIR/COMPARISON.md" | tee -a "$RESULTS_DIR/progress.log"
fi

echo "" | tee -a "$RESULTS_DIR/progress.log"
echo "========================================" | tee -a "$RESULTS_DIR/progress.log"
echo "POST-PROCESS COVERAGE RUN COMPLETE" | tee -a "$RESULTS_DIR/progress.log"
echo "========================================" | tee -a "$RESULTS_DIR/progress.log"





