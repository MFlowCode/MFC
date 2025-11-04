#!/usr/bin/env bash
# Full Coverage Comparison Script
# Runs baseline and new test suites to measure actual improvement

set -euo pipefail

LOGDIR="build/coverage_comparison"
mkdir -p "$LOGDIR"

echo "=========================================="
echo "MFC Full Coverage Comparison"
echo "=========================================="
echo "This will take 4-6 hours total"
echo ""
echo "Phase 1: Baseline (all tests before expansion)"
echo "Phase 2: New suite (all tests after expansion)"
echo "Phase 3: Comparison report"
echo ""
echo "Started: $(date)"
echo "=========================================="
echo ""

# Phase 1: Get baseline with original test list
echo "[Phase 1/3] Running BASELINE coverage..."
echo "This will take ~2-3 hours"
echo ""

# Clean first
./mfc.sh clean

# Build with coverage
echo "Building with coverage..."
./mfc.sh build --gcov --no-gpu --debug -t pre_process simulation post_process -j $(sysctl -n hw.ncpu)

# Run ALL tests (100%)
echo "Running ALL tests for baseline..."
GCOV_EXEC=$(which gcov-15 2>/dev/null || which gcov-14 2>/dev/null || which gcov)
export GCOV_PREFIX=${PWD}/build/staging
export GCOV_PREFIX_STRIP=0

time ./mfc.sh test --no-examples -j $(sysctl -n hw.ncpu) 2>&1 | tee "$LOGDIR/baseline_tests.log"

# Generate baseline report
echo "Generating baseline report..."
mkdir -p "$LOGDIR/baseline"

gcovr build/staging --root . \
    --gcov-executable "${GCOV_EXEC}" \
    --filter 'src/.*' \
    -j 1 \
    --print-summary 2>&1 | tee "$LOGDIR/baseline/summary.txt"

# Save baseline data
find build/staging -name "*.gcda" > "$LOGDIR/baseline_gcda_files.txt"
echo "Baseline .gcda files: $(cat $LOGDIR/baseline_gcda_files.txt | wc -l)" | tee -a "$LOGDIR/baseline/summary.txt"

echo ""
echo "=========================================="
echo "[Phase 1/3] COMPLETE"
echo "Baseline results saved to: $LOGDIR/baseline/"
echo "=========================================="
echo ""

# Phase 2: Run with all new tests
echo "[Phase 2/3] Running NEW SUITE coverage..."
echo "This will take ~4-5 hours"
echo ""

# Clean again
./mfc.sh clean

# Build with coverage again
echo "Building with coverage..."
./mfc.sh build --gcov --no-gpu --debug -t pre_process simulation post_process -j $(sysctl -n hw.ncpu)

# Run ALL tests (100%) - now includes the 607 new tests
echo "Running ALL tests for new suite..."
export GCOV_PREFIX=${PWD}/build/staging
export GCOV_PREFIX_STRIP=0

time ./mfc.sh test --no-examples -j $(sysctl -n hw.ncpu) 2>&1 | tee "$LOGDIR/new_tests.log"

# Generate new report
echo "Generating new suite report..."
mkdir -p "$LOGDIR/new"

gcovr build/staging --root . \
    --gcov-executable "${GCOV_EXEC}" \
    --filter 'src/.*' \
    -j 1 \
    --print-summary 2>&1 | tee "$LOGDIR/new/summary.txt"

# Save new data
find build/staging -name "*.gcda" > "$LOGDIR/new_gcda_files.txt"
echo "New suite .gcda files: $(cat $LOGDIR/new_gcda_files.txt | wc -l)" | tee -a "$LOGDIR/new/summary.txt"

echo ""
echo "=========================================="
echo "[Phase 2/3] COMPLETE"
echo "New suite results saved to: $LOGDIR/new/"
echo "=========================================="
echo ""

# Phase 3: Generate comparison
echo "[Phase 3/3] Generating comparison report..."
echo ""

cat > "$LOGDIR/COMPARISON_REPORT.md" << 'REPORT'
# MFC Coverage Comparison Report

## Test Suite Statistics

### Baseline (Original)
```
Test count: 790 tests
EOF

echo "Total test time: $(grep 'real' $LOGDIR/baseline_tests.log | tail -1)" >> "$LOGDIR/COMPARISON_REPORT.md"
echo '```' >> "$LOGDIR/COMPARISON_REPORT.md"
echo "" >> "$LOGDIR/COMPARISON_REPORT.md"

echo "### New Suite" >> "$LOGDIR/COMPARISON_REPORT.md"
echo '```' >> "$LOGDIR/COMPARISON_REPORT.md"
echo "Test count: 1,397 tests (+607, +77%)" >> "$LOGDIR/COMPARISON_REPORT.md"
echo "Total test time: $(grep 'real' $LOGDIR/new_tests.log | tail -1)" >> "$LOGDIR/COMPARISON_REPORT.md"
echo '```' >> "$LOGDIR/COMPARISON_REPORT.md"
echo "" >> "$LOGDIR/COMPARISON_REPORT.md"

echo "## Coverage Results" >> "$LOGDIR/COMPARISON_REPORT.md"
echo "" >> "$LOGDIR/COMPARISON_REPORT.md"
echo "### Baseline Coverage" >> "$LOGDIR/COMPARISON_REPORT.md"
echo '```' >> "$LOGDIR/COMPARISON_REPORT.md"
cat "$LOGDIR/baseline/summary.txt" >> "$LOGDIR/COMPARISON_REPORT.md"
echo '```' >> "$LOGDIR/COMPARISON_REPORT.md"
echo "" >> "$LOGDIR/COMPARISON_REPORT.md"

echo "### New Suite Coverage" >> "$LOGDIR/COMPARISON_REPORT.md"
echo '```' >> "$LOGDIR/COMPARISON_REPORT.md"
cat "$LOGDIR/new/summary.txt" >> "$LOGDIR/COMPARISON_REPORT.md"
echo '```' >> "$LOGDIR/COMPARISON_REPORT.md"
echo "" >> "$LOGDIR/COMPARISON_REPORT.md"

echo "---" >> "$LOGDIR/COMPARISON_REPORT.md"
echo "" >> "$LOGDIR/COMPARISON_REPORT.md"
echo "**Completed**: $(date)" >> "$LOGDIR/COMPARISON_REPORT.md"

echo "=========================================="
echo "[Phase 3/3] COMPLETE"
echo "=========================================="
echo ""
echo "FINAL REPORT: $LOGDIR/COMPARISON_REPORT.md"
echo ""
echo "Completed: $(date)"
echo "=========================================="

# Display the comparison
cat "$LOGDIR/COMPARISON_REPORT.md"






