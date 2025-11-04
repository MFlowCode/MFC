#!/usr/bin/env bash
# Comprehensive Coverage Comparison Script
# This script will:
# 1. Run baseline coverage with original test suite (with post-processing)
# 2. Add new tests (excluding broken ones)
# 3. Run expanded coverage with new test suite (with post-processing)
# 4. Compare results and runtimes

set -euo pipefail

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

RESULTS_DIR="build/coverage_comparison_full"
mkdir -p "$RESULTS_DIR"

# Redirect all output
exec > >(tee "$RESULTS_DIR/full_run.log") 2>&1

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}MFC Comprehensive Coverage Comparison${NC}"
echo -e "${GREEN}========================================${NC}"
echo "Started: $(date)"
echo ""

# Configuration
JOBS=${JOBS:-$(sysctl -n hw.ncpu 2>/dev/null || echo 4)}
GCOV_EXEC=$(which gcov-15 2>/dev/null || which gcov-14 2>/dev/null || which gcov)

echo "Configuration:"
echo "  Jobs: ${JOBS}"
echo "  Gcov: ${GCOV_EXEC}"
echo ""

# ============================================================================
# PHASE 1: Baseline Coverage (Original Test Suite with Post-Processing)
# ============================================================================

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}PHASE 1: Baseline Coverage${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

# Back up current cases.py
cp toolchain/mfc/test/cases.py toolchain/mfc/test/cases.py.backup

echo -e "${YELLOW}[1.1] Cleaning build directory...${NC}"
./mfc.sh clean

echo -e "${YELLOW}[1.2] Building with coverage instrumentation...${NC}"
BASELINE_BUILD_START=$(date +%s)
./mfc.sh build --gcov --no-gpu --debug -t pre_process simulation post_process -j ${JOBS}
BASELINE_BUILD_END=$(date +%s)
BASELINE_BUILD_TIME=$((BASELINE_BUILD_END - BASELINE_BUILD_START))
echo "Baseline build time: ${BASELINE_BUILD_TIME}s ($((BASELINE_BUILD_TIME / 60))m)"

echo -e "${YELLOW}[1.3] Counting baseline tests...${NC}"
BASELINE_TEST_COUNT=$(./mfc.sh test --list 2>&1 | grep -c "^tests/" || echo "unknown")
echo "Baseline test count: ${BASELINE_TEST_COUNT}"

echo -e "${YELLOW}[1.4] Running baseline tests (WITH post-processing, NO --no-examples)...${NC}"
BASELINE_TEST_START=$(date +%s)
./mfc.sh test --no-build -j ${JOBS} > "$RESULTS_DIR/baseline_test_output.txt" 2>&1 || {
    echo -e "${YELLOW}Some baseline tests failed, continuing...${NC}"
}
BASELINE_TEST_END=$(date +%s)
BASELINE_TEST_TIME=$((BASELINE_TEST_END - BASELINE_TEST_START))
echo "Baseline test time: ${BASELINE_TEST_TIME}s ($((BASELINE_TEST_TIME / 60))m)"

echo -e "${YELLOW}[1.5] Generating baseline coverage report...${NC}"
mkdir -p "$RESULTS_DIR/baseline"

# Find all build directories
BUILD_DIRS=$(find build/staging -type d -name "CMakeFiles" 2>/dev/null | sed 's|/CMakeFiles||' | head -20)

echo "Processing $(echo "$BUILD_DIRS" | wc -l) build directories for baseline..."

for BUILD_DIR in $BUILD_DIRS; do
    gcovr "$BUILD_DIR" \
        --root . \
        --gcov-executable "${GCOV_EXEC}" \
        --filter 'src/.*' \
        -j 1 \
        --gcov-ignore-parse-errors=suspicious_hits.warn \
        --print-summary 2>&1 | tee -a "$RESULTS_DIR/baseline/summary.txt" || {
        echo "  (had issues with $BUILD_DIR, continuing...)"
    }
done

echo ""
echo -e "${GREEN}Phase 1 Complete!${NC}"
echo "  Build time: ${BASELINE_BUILD_TIME}s"
echo "  Test time: ${BASELINE_TEST_TIME}s"
echo "  Total time: $((BASELINE_BUILD_TIME + BASELINE_TEST_TIME))s"
echo ""

# ============================================================================
# PHASE 2: Add New Tests (Excluding Broken Ones)
# ============================================================================

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}PHASE 2: Adding New Tests${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

echo -e "${YELLOW}[2.1] Modifying cases.py to add working new tests...${NC}"

# Create modified cases.py with new tests (excluding broken ones)
cat > /tmp/test_additions.py << 'PYTHON_CODE'
# New test functions to add (excluding broken configurations)

def alter_time_integrators():
    """Test different Runge-Kutta time integrators"""
    # time_stepper: 1=Euler, 2=RK2, 3=RK3 (default), 4=RK4, 5=RK5, 23=TVD-RK3
    # NOTE: Excluding configurations that might conflict
    for time_stepper in [1, 2, 4, 5, 23]:
        cases.append(define_case_d(stack, f"time_stepper={time_stepper}",
            {'time_stepper': time_stepper, 't_step_stop': 5}))

def alter_riemann_solvers_extended(num_fluids):
    """Extended Riemann solver testing"""
    # Include solvers 3 and 4, but ONLY with compatible configurations
    # Solver 3 (Exact) - don't use with wave_speeds parameter
    # Solver 4 (HLLD) - skip (requires MHD which we're not testing)
    
    for riemann_solver in [1, 5, 2]:  # Keep original working solvers
        stack.push(f"riemann_solver={riemann_solver}", {'riemann_solver': riemann_solver})
        
        if num_fluids <= 2:
            cases.append(define_case_d(stack, f"mixture_err", {'mixture_err': 'T'}))
        
        if riemann_solver != 2:
            cases.append(define_case_d(stack, f"avg_state=1", {'avg_state': 1}))
        
        if riemann_solver != 3:  # Solver 3 doesn't support wave_speeds
            cases.append(define_case_d(stack, f"wave_speeds=2", {'wave_speeds': 2}))
        
        if riemann_solver == 1:
            if num_fluids == 2:
                cases.append(define_case_d(stack, f"mpp_lim", {'mpp_lim': 'T'}))
        
        if riemann_solver == 2:
            if num_fluids <= 2:
                cases.append(define_case_d(stack, f"avg_state=1", {'avg_state': 1}))
            
            cases.append(define_case_d(stack, f"model_eqns=3", {'model_eqns': 3}))
            
            if num_fluids == 2:
                cases.append(define_case_d(stack, f"alt_soundspeed", {'alt_soundspeed': 'T'}))
                cases.append(define_case_d(stack, f"mpp_lim", {'mpp_lim': 'T'}))
        
        cases.append(define_case_d(stack, f"low_Mach=1", {'low_Mach': 1}))
        
        if riemann_solver == 2:
            cases.append(define_case_d(stack, f"low_Mach=2", {'low_Mach': 2}))
        
        stack.pop()

def alter_cfl_modes_safe():
    """Test CFL adaptation modes with safe parameters"""
    # Adaptive CFL with proper t_stop
    cases.append(define_case_d(stack, "cfl_adap_dt=T",
        {'cfl_adap_dt': 'T', 'cfl_target': 0.5, 't_step_stop': 10}))
    
    # Note: cfl_const_dt tests were causing t_stop <= 0 errors, so we'll skip them

def alter_model_equations_safe():
    """Test different model equation formulations (safe versions)"""
    # model_eqns=1 doesn't support num_fluids, so skip
    # Test only models 2 and 3 which are more flexible
    for model_eqns in [2, 3]:
        cases.append(define_case_d(stack, f"model_eqns={model_eqns}",
            {'model_eqns': model_eqns}))

def alter_grid_stretching_safe(dimInfo):
    """Test grid stretching (only in supported dimensions)"""
    # x_stretch seems to have issues in 3D, so only test in 1D and 2D
    if len(dimInfo[0]) <= 2:
        cases.append(define_case_d(stack, "x_stretch=T",
            {'x_stretch': 'T', 'a_x': 1.5, 'x_a': -1.0, 'x_b': 1.0}))
PYTHON_CODE

# Now insert these functions into cases.py
python3 << 'PYTHON_SCRIPT'
import re

# Read original file
with open('toolchain/mfc/test/cases.py', 'r') as f:
    content = f.read()

# Read new functions
with open('/tmp/test_additions.py', 'r') as f:
    new_functions = f.read()

# Find the foreach_dimension function and add calls to new test functions
# Look for the pattern in foreach_dimension where we call alter functions

pattern = r'(def foreach_dimension\(\):.*?)(alter_lag_bubbles\(dimInfo\))'
replacement = r'\1\2\n        alter_time_integrators()\n        alter_cfl_modes_safe()\n        alter_model_equations_safe()\n        alter_grid_stretching_safe(dimInfo)'

content_modified = re.sub(pattern, replacement, content, flags=re.DOTALL)

# Also update alter_riemann_solvers to call the extended version
# But actually, let's keep it simple and not modify the existing function

# Insert new functions before list_cases
insert_position = content_modified.find('def list_cases()')
if insert_position > 0:
    content_final = content_modified[:insert_position] + new_functions + '\n\n' + content_modified[insert_position:]
else:
    print("ERROR: Could not find list_cases function")
    content_final = content_modified

# Write modified file
with open('toolchain/mfc/test/cases.py', 'w') as f:
    f.write(content_final)

print("Successfully modified cases.py")
PYTHON_SCRIPT

echo "Modified cases.py with new test functions"

echo -e "${YELLOW}[2.2] Verifying expanded test count...${NC}"
EXPANDED_TEST_COUNT=$(./mfc.sh test --list 2>&1 | grep -c "^tests/" || echo "unknown")
echo "Expanded test count: ${EXPANDED_TEST_COUNT}"
NEW_TESTS_ADDED=$((EXPANDED_TEST_COUNT - BASELINE_TEST_COUNT))
echo "New tests added: ${NEW_TESTS_ADDED}"

# ============================================================================
# PHASE 3: Expanded Coverage (New Test Suite with Post-Processing)
# ============================================================================

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}PHASE 3: Expanded Coverage${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

echo -e "${YELLOW}[3.1] Cleaning build directory...${NC}"
./mfc.sh clean

echo -e "${YELLOW}[3.2] Building with coverage instrumentation...${NC}"
EXPANDED_BUILD_START=$(date +%s)
./mfc.sh build --gcov --no-gpu --debug -t pre_process simulation post_process -j ${JOBS}
EXPANDED_BUILD_END=$(date +%s)
EXPANDED_BUILD_TIME=$((EXPANDED_BUILD_END - EXPANDED_BUILD_START))
echo "Expanded build time: ${EXPANDED_BUILD_TIME}s ($((EXPANDED_BUILD_TIME / 60))m)"

echo -e "${YELLOW}[3.3] Running expanded tests (WITH post-processing, NO --no-examples)...${NC}"
EXPANDED_TEST_START=$(date +%s)
./mfc.sh test --no-build -j ${JOBS} > "$RESULTS_DIR/expanded_test_output.txt" 2>&1 || {
    echo -e "${YELLOW}Some expanded tests failed, continuing...${NC}"
}
EXPANDED_TEST_END=$(date +%s)
EXPANDED_TEST_TIME=$((EXPANDED_TEST_END - EXPANDED_TEST_START))
echo "Expanded test time: ${EXPANDED_TEST_TIME}s ($((EXPANDED_TEST_TIME / 60))m)"

echo -e "${YELLOW}[3.4] Generating expanded coverage report...${NC}"
mkdir -p "$RESULTS_DIR/expanded"

# Find all build directories
BUILD_DIRS=$(find build/staging -type d -name "CMakeFiles" 2>/dev/null | sed 's|/CMakeFiles||' | head -20)

echo "Processing $(echo "$BUILD_DIRS" | wc -l) build directories for expanded..."

for BUILD_DIR in $BUILD_DIRS; do
    gcovr "$BUILD_DIR" \
        --root . \
        --gcov-executable "${GCOV_EXEC}" \
        --filter 'src/.*' \
        -j 1 \
        --gcov-ignore-parse-errors=suspicious_hits.warn \
        --print-summary 2>&1 | tee -a "$RESULTS_DIR/expanded/summary.txt" || {
        echo "  (had issues with $BUILD_DIR, continuing...)"
    }
done

echo ""
echo -e "${GREEN}Phase 3 Complete!${NC}"
echo "  Build time: ${EXPANDED_BUILD_TIME}s"
echo "  Test time: ${EXPANDED_TEST_TIME}s"
echo "  Total time: $((EXPANDED_BUILD_TIME + EXPANDED_TEST_TIME))s"
echo ""

# ============================================================================
# PHASE 4: Comparison and Analysis
# ============================================================================

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}PHASE 4: Comparison and Analysis${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

# Generate comprehensive comparison report
cat > "$RESULTS_DIR/COMPARISON_REPORT.md" << 'REPORT_EOF'
# MFC Coverage Comparison Report

**Generated**: $(date)

## Test Suite Statistics

### Baseline (Original Test Suite)
- **Test Count**: BASELINE_TEST_COUNT_PLACEHOLDER
- **Build Time**: BASELINE_BUILD_TIME_PLACEHOLDERs (BASELINE_BUILD_MIN_PLACEHOLDERm)
- **Test Execution Time**: BASELINE_TEST_TIME_PLACEHOLDERs (BASELINE_TEST_MIN_PLACEHOLDERm)
- **Total Time**: BASELINE_TOTAL_TIME_PLACEHOLDERs (BASELINE_TOTAL_MIN_PLACEHOLDERm)

### Expanded (With New Tests)
- **Test Count**: EXPANDED_TEST_COUNT_PLACEHOLDER (+NEW_TESTS_ADDED_PLACEHOLDER, +PERCENT_INCREASE_PLACEHOLDER%)
- **Build Time**: EXPANDED_BUILD_TIME_PLACEHOLDERs (EXPANDED_BUILD_MIN_PLACEHOLDERm)
- **Test Execution Time**: EXPANDED_TEST_TIME_PLACEHOLDERs (EXPANDED_TEST_MIN_PLACEHOLDERm)
- **Total Time**: EXPANDED_TOTAL_TIME_PLACEHOLDERs (EXPANDED_TOTAL_MIN_PLACEHOLDERm)

### Time Comparison
- **Build Time Change**: BUILD_TIME_DIFF_PLACEHOLDERs (BUILD_TIME_PERCENT_PLACEHOLDER%)
- **Test Time Change**: TEST_TIME_DIFF_PLACEHOLDERs (TEST_TIME_PERCENT_PLACEHOLDER%)
- **Total Time Change**: TOTAL_TIME_DIFF_PLACEHOLDERs (TOTAL_TIME_PERCENT_PLACEHOLDER%)

## Coverage Results

### Baseline Coverage
```
BASELINE_COVERAGE_PLACEHOLDER
```

### Expanded Coverage
```
EXPANDED_COVERAGE_PLACEHOLDER
```

## New Tests Added

The following test categories were added (excluding broken configurations):

1. **Time Integrators** (5 variants):
   - Euler (time_stepper=1)
   - RK2 (time_stepper=2)
   - RK4 (time_stepper=4)
   - RK5 (time_stepper=5)
   - TVD-RK3 (time_stepper=23)

2. **CFL Modes** (1 variant):
   - Adaptive CFL (cfl_adap_dt=T)
   - Note: Constant CFL tests excluded due to parameter conflicts

3. **Model Equations** (2 variants):
   - Pi-gamma model (model_eqns=2)
   - 5-equation model (model_eqns=3)
   - Note: Gamma model (model_eqns=1) excluded due to num_fluids conflict

4. **Grid Stretching** (1D and 2D only):
   - Non-uniform grid (x_stretch=T)
   - Note: 3D stretching tests excluded due to incompatibility

## Excluded Tests

The following test configurations were identified as broken and excluded:

1. **model_eqns=1** - Conflicts with num_fluids parameter
2. **riemann_solver=3** (Exact Riemann) - Conflicts with wave_speeds parameter
3. **riemann_solver=4** (HLLD) - Requires MHD mode (not tested)
4. **cfl_const_dt=T** - Causes t_stop <= 0 errors
5. **x_stretch=T in 3D** - Property not allowed
6. **loops_x=2** - Various configuration errors

## Analysis

### Coverage Improvement
- **Baseline**: See detailed numbers above
- **Expanded**: See detailed numbers above
- **Net Change**: Calculate from above

### Runtime Impact
- **Per-test overhead**: Calculate average test time
- **Scalability**: Assess if test time scales linearly with test count

## Recommendations

1. **Keep expanded test suite** if coverage improvement is significant (>5%)
2. **Investigate excluded tests** to understand parameter constraints
3. **Document parameter incompatibilities** to prevent future issues
4. **Consider post-processing tests** - now included in both runs

---

**Files Generated**:
- `baseline/summary.txt` - Baseline coverage data
- `expanded/summary.txt` - Expanded coverage data
- `baseline_test_output.txt` - Baseline test output
- `expanded_test_output.txt` - Expanded test output
- `full_run.log` - Complete execution log
REPORT_EOF

# Fill in the placeholders
BASELINE_TOTAL=$((BASELINE_BUILD_TIME + BASELINE_TEST_TIME))
EXPANDED_TOTAL=$((EXPANDED_BUILD_TIME + EXPANDED_TEST_TIME))
BUILD_DIFF=$((EXPANDED_BUILD_TIME - BASELINE_BUILD_TIME))
TEST_DIFF=$((EXPANDED_TEST_TIME - BASELINE_TEST_TIME))
TOTAL_DIFF=$((EXPANDED_TOTAL - BASELINE_TOTAL))

PERCENT_INCREASE=$(awk "BEGIN {printf \"%.1f\", ($NEW_TESTS_ADDED / $BASELINE_TEST_COUNT) * 100}")
BUILD_PERCENT=$(awk "BEGIN {printf \"%.1f\", ($BUILD_DIFF / $BASELINE_BUILD_TIME) * 100}")
TEST_PERCENT=$(awk "BEGIN {printf \"%.1f\", ($TEST_DIFF / $BASELINE_TEST_TIME) * 100}")
TOTAL_PERCENT=$(awk "BEGIN {printf \"%.1f\", ($TOTAL_DIFF / $BASELINE_TOTAL) * 100}")

sed -i.bak \
    -e "s/BASELINE_TEST_COUNT_PLACEHOLDER/$BASELINE_TEST_COUNT/g" \
    -e "s/BASELINE_BUILD_TIME_PLACEHOLDER/$BASELINE_BUILD_TIME/g" \
    -e "s/BASELINE_BUILD_MIN_PLACEHOLDER/$((BASELINE_BUILD_TIME / 60))/g" \
    -e "s/BASELINE_TEST_TIME_PLACEHOLDER/$BASELINE_TEST_TIME/g" \
    -e "s/BASELINE_TEST_MIN_PLACEHOLDER/$((BASELINE_TEST_TIME / 60))/g" \
    -e "s/BASELINE_TOTAL_TIME_PLACEHOLDER/$BASELINE_TOTAL/g" \
    -e "s/BASELINE_TOTAL_MIN_PLACEHOLDER/$((BASELINE_TOTAL / 60))/g" \
    -e "s/EXPANDED_TEST_COUNT_PLACEHOLDER/$EXPANDED_TEST_COUNT/g" \
    -e "s/NEW_TESTS_ADDED_PLACEHOLDER/$NEW_TESTS_ADDED/g" \
    -e "s/PERCENT_INCREASE_PLACEHOLDER/$PERCENT_INCREASE/g" \
    -e "s/EXPANDED_BUILD_TIME_PLACEHOLDER/$EXPANDED_BUILD_TIME/g" \
    -e "s/EXPANDED_BUILD_MIN_PLACEHOLDER/$((EXPANDED_BUILD_TIME / 60))/g" \
    -e "s/EXPANDED_TEST_TIME_PLACEHOLDER/$EXPANDED_TEST_TIME/g" \
    -e "s/EXPANDED_TEST_MIN_PLACEHOLDER/$((EXPANDED_TEST_TIME / 60))/g" \
    -e "s/EXPANDED_TOTAL_TIME_PLACEHOLDER/$EXPANDED_TOTAL/g" \
    -e "s/EXPANDED_TOTAL_MIN_PLACEHOLDER/$((EXPANDED_TOTAL / 60))/g" \
    -e "s/BUILD_TIME_DIFF_PLACEHOLDER/$BUILD_DIFF/g" \
    -e "s/BUILD_TIME_PERCENT_PLACEHOLDER/$BUILD_PERCENT/g" \
    -e "s/TEST_TIME_DIFF_PLACEHOLDER/$TEST_DIFF/g" \
    -e "s/TEST_TIME_PERCENT_PLACEHOLDER/$TEST_PERCENT/g" \
    -e "s/TOTAL_TIME_DIFF_PLACEHOLDER/$TOTAL_DIFF/g" \
    -e "s/TOTAL_TIME_PERCENT_PLACEHOLDER/$TOTAL_PERCENT/g" \
    "$RESULTS_DIR/COMPARISON_REPORT.md"

# Insert coverage data
BASELINE_COV=$(tail -20 "$RESULTS_DIR/baseline/summary.txt" | grep -A 10 "TOTAL" | head -15 || echo "No coverage data")
EXPANDED_COV=$(tail -20 "$RESULTS_DIR/expanded/summary.txt" | grep -A 10 "TOTAL" | head -15 || echo "No coverage data")

# This is a bit tricky with sed, so let's use a different approach
python3 << PYTHON_EOF
import re

with open('$RESULTS_DIR/COMPARISON_REPORT.md', 'r') as f:
    content = f.read()

baseline_cov = '''$BASELINE_COV'''
expanded_cov = '''$EXPANDED_COV'''

content = content.replace('BASELINE_COVERAGE_PLACEHOLDER', baseline_cov)
content = content.replace('EXPANDED_COVERAGE_PLACEHOLDER', expanded_cov)

with open('$RESULTS_DIR/COMPARISON_REPORT.md', 'w') as f:
    f.write(content)
PYTHON_EOF

# Restore original cases.py
mv toolchain/mfc/test/cases.py.backup toolchain/mfc/test/cases.py

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}ALL PHASES COMPLETE!${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo "Completed: $(date)"
echo ""
echo "Summary:"
echo "--------"
echo "Baseline:"
echo "  Tests: ${BASELINE_TEST_COUNT}"
echo "  Total Time: ${BASELINE_TOTAL}s ($((BASELINE_TOTAL / 60))m)"
echo ""
echo "Expanded:"
echo "  Tests: ${EXPANDED_TEST_COUNT} (+${NEW_TESTS_ADDED}, +${PERCENT_INCREASE}%)"
echo "  Total Time: ${EXPANDED_TOTAL}s ($((EXPANDED_TOTAL / 60))m)"
echo ""
echo "Time Difference: ${TOTAL_DIFF}s (+${TOTAL_PERCENT}%)"
echo ""
echo "Results saved to: $RESULTS_DIR/"
echo "  - COMPARISON_REPORT.md"
echo "  - baseline/summary.txt"
echo "  - expanded/summary.txt"
echo ""
echo -e "${GREEN}View the full report:${NC}"
echo "  cat $RESULTS_DIR/COMPARISON_REPORT.md"





