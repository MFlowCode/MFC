#!/usr/bin/env bash
# Direct coverage comparison - no buffering issues
# Simplified approach that writes directly to files

set -euo pipefail

RESULTS_DIR="coverage_results"
mkdir -p "$RESULTS_DIR"
touch "$RESULTS_DIR/progress.log"

echo "========================================" > "$RESULTS_DIR/progress.log"
echo "MFC Coverage Comparison - Started $(date)" >> "$RESULTS_DIR/progress.log"
echo "========================================" >> "$RESULTS_DIR/progress.log"

JOBS=$(sysctl -n hw.ncpu 2>/dev/null || echo 4)
GCOV_EXEC=$(which gcov-15 2>/dev/null || which gcov-14 2>/dev/null || which gcov)

echo "Jobs: ${JOBS}" >> "$RESULTS_DIR/progress.log"
echo "Gcov: ${GCOV_EXEC}" >> "$RESULTS_DIR/progress.log"
echo "" >> "$RESULTS_DIR/progress.log"

# ============================================================================
# PHASE 1: Baseline Coverage
# ============================================================================

echo "PHASE 1: Baseline Coverage - $(date)" >> "$RESULTS_DIR/progress.log"

# Clean
echo "  Cleaning..." >> "$RESULTS_DIR/progress.log"
./mfc.sh clean >> "$RESULTS_DIR/progress.log" 2>&1

# Build
echo "  Building with coverage - $(date)" >> "$RESULTS_DIR/progress.log"
BASELINE_BUILD_START=$(date +%s)
./mfc.sh build --gcov --no-gpu --debug -t pre_process simulation post_process -j ${JOBS} >> "$RESULTS_DIR/baseline_build.log" 2>&1
BASELINE_BUILD_END=$(date +%s)
BASELINE_BUILD_TIME=$((BASELINE_BUILD_END - BASELINE_BUILD_START))
echo "  Build complete: ${BASELINE_BUILD_TIME}s" >> "$RESULTS_DIR/progress.log"

# Count tests
BASELINE_COUNT=$(./mfc.sh test --list 2>&1 | grep -E "^ *[A-F0-9]{8} " | wc -l | tr -d ' ')
echo "  Baseline tests: ${BASELINE_COUNT}" >> "$RESULTS_DIR/progress.log"

# Run tests
echo "  Running tests - $(date)" >> "$RESULTS_DIR/progress.log"
BASELINE_TEST_START=$(date +%s)
./mfc.sh test --no-build -j ${JOBS} >> "$RESULTS_DIR/baseline_tests.log" 2>&1 || echo "Some tests failed" >> "$RESULTS_DIR/progress.log"
BASELINE_TEST_END=$(date +%s)
BASELINE_TEST_TIME=$((BASELINE_TEST_END - BASELINE_TEST_START))
echo "  Tests complete: ${BASELINE_TEST_TIME}s" >> "$RESULTS_DIR/progress.log"

# Generate coverage
echo "  Generating coverage report - $(date)" >> "$RESULTS_DIR/progress.log"
gcovr build/staging \
    --root . \
    --gcov-executable "${GCOV_EXEC}" \
    --filter 'src/.*' \
    --gcov-ignore-parse-errors=suspicious_hits.warn \
    --print-summary \
    -j 1 > "$RESULTS_DIR/baseline_coverage.txt" 2>&1 || echo "Coverage had issues" >> "$RESULTS_DIR/progress.log"

echo "Phase 1 complete - $(date)" >> "$RESULTS_DIR/progress.log"
echo "" >> "$RESULTS_DIR/progress.log"

# ============================================================================
# PHASE 2: Add Tests
# ============================================================================

echo "PHASE 2: Adding safe new tests - $(date)" >> "$RESULTS_DIR/progress.log"

# Backup original
cp toolchain/mfc/test/cases.py toolchain/mfc/test/cases.py.original

# Add new test functions - directly modify the file
python3 << 'EOF'
import re

with open('toolchain/mfc/test/cases.py', 'r') as f:
    content = f.read()

# Find alter_riemann_solvers and modify it to include solvers 3 and 4
# Original: for riemann_solver in [1, 5, 2]:
# New: for riemann_solver in [1, 5, 2, 3]:
content = content.replace(
    'for riemann_solver in [1, 5, 2]:',
    'for riemann_solver in [1, 5, 2, 3]:'
)

# Add new test functions before list_cases()
new_functions = '''
def alter_time_integrators():
    """Test different time integrators (safe configurations only)"""
    # RK2, RK4, RK5 (skip Euler and TVD-RK3 which may have issues)
    for ts in [2, 4, 5]:
        cases.append(define_case_d(stack, f"time_stepper={ts}",
            {'time_stepper': ts, 't_step_stop': 5}))

def alter_cfl_adaptive():
    """Test adaptive CFL"""
    cases.append(define_case_d(stack, "cfl_adap_dt=T",
        {'cfl_adap_dt': 'T', 'cfl_target': 0.5, 't_step_stop': 10}))

'''

insert_pos = content.find('def list_cases()')
if insert_pos > 0:
    content = content[:insert_pos] + new_functions + content[insert_pos:]

# Add calls in foreach_dimension
# Find the alter_muscl() call and add our new tests after it
pattern = r'(alter_muscl\(\))'
replacement = r'\1\n        alter_time_integrators()\n        alter_cfl_adaptive()'
content = re.sub(pattern, replacement, content)

with open('toolchain/mfc/test/cases.py', 'w') as f:
    f.write(content)

print("Successfully added tests")
EOF

# Count expanded tests
EXPANDED_COUNT=$(./mfc.sh test --list 2>&1 | grep -E "^ *[A-F0-9]{8} " | wc -l | tr -d ' ')
NEW_COUNT=$((EXPANDED_COUNT - BASELINE_COUNT))
echo "  Expanded tests: ${EXPANDED_COUNT} (+${NEW_COUNT})" >> "$RESULTS_DIR/progress.log"
echo "" >> "$RESULTS_DIR/progress.log"

# ============================================================================
# PHASE 3: Expanded Coverage
# ============================================================================

echo "PHASE 3: Expanded Coverage - $(date)" >> "$RESULTS_DIR/progress.log"

# Clean
echo "  Cleaning..." >> "$RESULTS_DIR/progress.log"
./mfc.sh clean >> "$RESULTS_DIR/progress.log" 2>&1

# Build
echo "  Building with coverage - $(date)" >> "$RESULTS_DIR/progress.log"
EXPANDED_BUILD_START=$(date +%s)
./mfc.sh build --gcov --no-gpu --debug -t pre_process simulation post_process -j ${JOBS} >> "$RESULTS_DIR/expanded_build.log" 2>&1
EXPANDED_BUILD_END=$(date +%s)
EXPANDED_BUILD_TIME=$((EXPANDED_BUILD_END - EXPANDED_BUILD_START))
echo "  Build complete: ${EXPANDED_BUILD_TIME}s" >> "$RESULTS_DIR/progress.log"

# Run tests
echo "  Running tests - $(date)" >> "$RESULTS_DIR/progress.log"
EXPANDED_TEST_START=$(date +%s)
./mfc.sh test --no-build -j ${JOBS} >> "$RESULTS_DIR/expanded_tests.log" 2>&1 || echo "Some tests failed" >> "$RESULTS_DIR/progress.log"
EXPANDED_TEST_END=$(date +%s)
EXPANDED_TEST_TIME=$((EXPANDED_TEST_END - EXPANDED_TEST_START))
echo "  Tests complete: ${EXPANDED_TEST_TIME}s" >> "$RESULTS_DIR/progress.log"

# Generate coverage
echo "  Generating coverage report - $(date)" >> "$RESULTS_DIR/progress.log"
gcovr build/staging \
    --root . \
    --gcov-executable "${GCOV_EXEC}" \
    --filter 'src/.*' \
    --gcov-ignore-parse-errors=suspicious_hits.warn \
    --print-summary \
    -j 1 > "$RESULTS_DIR/expanded_coverage.txt" 2>&1 || echo "Coverage had issues" >> "$RESULTS_DIR/progress.log"

echo "Phase 3 complete - $(date)" >> "$RESULTS_DIR/progress.log"
echo "" >> "$RESULTS_DIR/progress.log"

# ============================================================================
# PHASE 4: Comparison
# ============================================================================

echo "PHASE 4: Generating comparison report - $(date)" >> "$RESULTS_DIR/progress.log"

# Restore original cases.py
mv toolchain/mfc/test/cases.py.original toolchain/mfc/test/cases.py

# Generate report
cat > "$RESULTS_DIR/FINAL_REPORT.md" << EOF_REPORT
# MFC Coverage Comparison - Final Report

**Generated**: $(date)

## Test Suite Comparison

| Metric | Baseline | Expanded | Change |
|--------|----------|----------|--------|
| **Test Count** | ${BASELINE_COUNT} | ${EXPANDED_COUNT} | +${NEW_COUNT} (+$(awk "BEGIN {printf \"%.1f\", ($NEW_COUNT*100/$BASELINE_COUNT)}")%) |
| **Build Time** | ${BASELINE_BUILD_TIME}s | ${EXPANDED_BUILD_TIME}s | $((EXPANDED_BUILD_TIME - BASELINE_BUILD_TIME))s |
| **Test Time** | ${BASELINE_TEST_TIME}s | ${EXPANDED_TEST_TIME}s | $((EXPANDED_TEST_TIME - BASELINE_TEST_TIME))s |
| **Total Time** | $((BASELINE_BUILD_TIME + BASELINE_TEST_TIME))s | $((EXPANDED_BUILD_TIME + EXPANDED_TEST_TIME))s | $((EXPANDED_BUILD_TIME + EXPANDED_TEST_TIME - BASELINE_BUILD_TIME - BASELINE_TEST_TIME))s |

## Coverage Results

### Baseline Coverage
\`\`\`
$(cat "$RESULTS_DIR/baseline_coverage.txt")
\`\`\`

### Expanded Coverage
\`\`\`
$(cat "$RESULTS_DIR/expanded_coverage.txt")
\`\`\`

## Tests Added

1. **Time Integrators**: RK2, RK4, RK5 (3 variants × dimensions)
2. **Adaptive CFL**: cfl_adap_dt=T (1 variant × dimensions)
3. **Riemann Solver 3**: Exact Riemann solver (added to existing tests)

**Total new tests**: ${NEW_COUNT}

## Files Generated

- \`progress.log\` - Execution timeline
- \`baseline_build.log\` - Baseline build output
- \`baseline_tests.log\` - Baseline test output
- \`baseline_coverage.txt\` - Baseline coverage report
- \`expanded_build.log\` - Expanded build output
- \`expanded_tests.log\` - Expanded test output
- \`expanded_coverage.txt\` - Expanded coverage report

---

**Completed**: $(date)
EOF_REPORT

echo "========================================" >> "$RESULTS_DIR/progress.log"
echo "COMPLETE - $(date)" >> "$RESULTS_DIR/progress.log"
echo "========================================" >> "$RESULTS_DIR/progress.log"
echo "" >> "$RESULTS_DIR/progress.log"
echo "Results in: $RESULTS_DIR/FINAL_REPORT.md" >> "$RESULTS_DIR/progress.log"

# Print summary to console
echo ""
echo "========================================="
echo "MFC Coverage Comparison COMPLETE"
echo "========================================="
echo ""
cat "$RESULTS_DIR/FINAL_REPORT.md"

