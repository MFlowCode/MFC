# Coverage Troubleshooting Guide

## Issue: gcovr Worker Thread Exception

### Symptoms
```
RuntimeError: Worker thread raised exception, workers canceled.
```

Tests ran successfully (`.gcda` files created), but report generation failed.

---

## Root Cause

gcovr's parallel processing can fail when:
1. Processing too many files simultaneously
2. Memory constraints
3. File path complexity
4. gcov subprocess issues

---

## Solutions (In Order of Preference)

### Solution 1: Use Simplified Script ✅ **RECOMMENDED**

I created a more robust version:

```bash
cd /Users/spencer/Downloads/MFC
PERCENT=10 ./toolchain/coverage_simple.sh
```

This version:
- Uses single-threaded gcovr (`-j 1`)
- Better error handling  
- Simpler report generation
- More stable

**Status**: Currently running in background

---

### Solution 2: Reduce Parallelism in Original Script

Edit `toolchain/coverage.sh` line 56-62:

```bash
# BEFORE:
gcovr build/staging --root . \
    --gcov-executable "${GCOV_EXEC}" \
    --filter 'src/.*' \
    --html --html-details -o build/coverage/index.html \
    ...

# AFTER: Add -j 1
gcovr build/staging --root . \
    --gcov-executable "${GCOV_EXEC}" \
    --filter 'src/.*' \
    -j 1 \
    --html --html-details -o build/coverage/index.html \
    ...
```

---

### Solution 3: Generate Reports Separately

If you have `.gcda` files from a previous run:

```bash
cd /Users/spencer/Downloads/MFC

# Find correct gcov
GCOV_EXEC=$(which gcov-15 2>/dev/null || which gcov-14 2>/dev/null || which gcov)
echo "Using: ${GCOV_EXEC}"

# Create output directory
mkdir -p build/coverage

# Step 1: Text summary (most stable)
echo "Generating text summary..."
gcovr build/staging --root . \
    --gcov-executable "${GCOV_EXEC}" \
    --filter 'src/.*' \
    --print-summary | tee build/coverage/summary.txt

# Step 2: Simple HTML (no details, single thread)
echo "Generating HTML..."
gcovr build/staging --root . \
    --gcov-executable "${GCOV_EXEC}" \
    --filter 'src/.*' \
    -j 1 \
    --html -o build/coverage/index.html

# Step 3: XML for CI
echo "Generating XML..."
gcovr build/staging --root . \
    --gcov-executable "${GCOV_EXEC}" \
    --filter 'src/.*' \
    -j 1 \
    --xml -o build/coverage/coverage.xml
```

---

### Solution 4: Process Subset of Files

If reports are still failing, process files incrementally:

```bash
cd /Users/spencer/Downloads/MFC
GCOV_EXEC=$(which gcov-15 2>/dev/null || which gcov)

# Just simulation files
gcovr build/staging --root . \
    --gcov-executable "${GCOV_EXEC}" \
    --filter 'src/simulation/.*' \
    -j 1 \
    --html -o build/coverage/simulation.html

# Just pre_process files  
gcovr build/staging --root . \
    --gcov-executable "${GCOV_EXEC}" \
    --filter 'src/pre_process/.*' \
    -j 1 \
    --html -o build/coverage/pre_process.html
```

---

### Solution 5: Use gcov Directly (Most Basic)

Process one file at a time:

```bash
cd /Users/spencer/Downloads/MFC/build/staging/<some_build_hash>

# Find a .gcda file
ls CMakeFiles/simulation.dir/fypp/simulation/*.gcda

# Run gcov on it
gcov-15 -o . CMakeFiles/simulation.dir/fypp/simulation/m_time_steppers.fpp.f90.gcda

# View the .gcov file
ls *.gcov
head -50 m_time_steppers.fpp.f90.gcov
```

---

## Prevention: Best Practices

### 1. Start Small
```bash
# Use fewer tests initially
PERCENT=5 ./toolchain/coverage.sh

# Then scale up
PERCENT=10 ./toolchain/coverage.sh
PERCENT=25 ./toolchain/coverage.sh
```

### 2. Use Single-Threaded gcovr
Always add `-j 1` to gcovr commands for stability:
```bash
gcovr ... -j 1 ...
```

### 3. Generate Reports Separately
Don't combine `--html`, `--html-details`, `--xml` in one command:
```bash
# INSTEAD OF:
gcovr ... --html --html-details --xml ...

# DO:
gcovr ... --html -o index.html
gcovr ... --xml -o coverage.xml
gcovr ... --txt -o summary.txt
```

### 4. Clean Between Runs
```bash
./mfc.sh clean
rm -rf build/coverage
```

---

## Verification

### Check if .gcda Files Exist
```bash
find build/staging -name "*.gcda" | wc -l
```
**Expected**: > 0 (e.g., 162)

### Check if .gcno Files Exist
```bash
find build/staging -name "*.gcno" | wc -l  
```
**Expected**: > .gcda count (e.g., 195)

### Test gcov Manually
```bash
GCOV_EXEC=$(which gcov-15 2>/dev/null || which gcov)
$GCOV_EXEC --version
```
**Expected**: Version matching gfortran (e.g., 15.x.x)

---

## Current Status

### What Worked ✅
1. Build with coverage (`--gcov`) ✅
2. Tests ran successfully ✅  
3. `.gcda` files created (162 files) ✅
4. `.gcno` files present (195 files) ✅

### What Failed ❌
- gcovr parallel report generation ❌

### Current Solution 🔄
- Running simplified script (`coverage_simple.sh`)
- Single-threaded gcovr
- 10% of tests for faster results
- **Status**: In progress

---

## Quick Commands Reference

```bash
# Check if coverage run is active
ps aux | grep coverage

# View simplified script log
tail -f /tmp/coverage_simple.log

# Run manually (10% tests, ~10-15 min)
cd /Users/spencer/Downloads/MFC
PERCENT=10 ./toolchain/coverage_simple.sh

# Generate reports from existing data
GCOV_EXEC=$(which gcov-15 || which gcov)
gcovr build/staging --root . \
    --gcov-executable "${GCOV_EXEC}" \
    --filter 'src/.*' \
    -j 1 \
    --html -o build/coverage/index.html
```

---

## Expected Results (When Working)

### Text Summary
```
------------------------------------------------------------------------------
GCC Code Coverage Report
Directory: .
------------------------------------------------------------------------------
File                                       Lines    Exec  Cover   Missing
------------------------------------------------------------------------------
src/simulation/m_time_steppers.fpp          234     198   84%   ...
src/simulation/m_riemann_solvers.fpp       512     410   80%   ...
...
------------------------------------------------------------------------------
TOTAL                                     15234   10156   67%
------------------------------------------------------------------------------
Branches                                   4523    2034   45%
------------------------------------------------------------------------------
```

### HTML Report
- Color-coded source files
- Line-by-line coverage
- Per-file and per-function breakdown
- Interactive navigation

---

## If All Else Fails

The test suite expansion is still valuable even without exact coverage numbers:

✅ **607 new tests added** (790 → 1,397)  
✅ **Critical gaps addressed**:
  - Time integrators (was 0% tested)
  - Riemann solvers (expanded from 3 to 5)
  - Grid stretching (was 0% tested)
  - CFL modes
  - Model equations

The infrastructure is in place. You can:
1. Run tests normally: `./mfc.sh test`
2. Generate coverage when needed
3. Add more tests following `REGRESSION_TEST_EXPANSION.md`

---

**Last Updated**: Nov 1, 2025, 3:00 PM  
**Issue**: gcovr worker thread exception  
**Status**: Workaround implemented (`coverage_simple.sh`)  
**Next**: Wait for simplified run to complete or run manually






