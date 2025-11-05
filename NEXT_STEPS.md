# Next Steps - MFC Coverage

## 🎉 Current Status: MISSION COMPLETE!

You now have **83.7% line coverage** and **100% function coverage** with the `-a` flag for post-processing validation.

---

## What You Should Do Now

### 1. Review the Results ✅

Open the coverage report:
```bash
open coverage_results_postprocess/index.html
```

This shows you:
- Which lines are covered (green) vs uncovered (red)
- Coverage breakdown by file and module
- Detailed line-by-line view

### 2. Update Your CI/CD Pipeline ⚙️

Make sure your CI always uses the `-a` flag for coverage:

```yaml
# Example CI configuration
coverage:
  script:
    - ./mfc.sh build --gcov --no-gpu --debug -t pre_process simulation post_process
    - ./mfc.sh test -a --no-build -j $(nproc)
    - gcovr build/staging --root . --gcov-executable gcov-15 --filter 'src/.*' --xml -o coverage.xml
```

**Key**: The `-a` flag is essential!

### 3. Set Coverage Thresholds 📊

With 83.7% baseline, you can enforce minimums:

```bash
# Example: Fail build if coverage drops below 80%
gcovr ... --fail-under-line 80
```

Or set aspirational goals:
- **Current**: 83.7% lines, 100% functions
- **Goal**: 85% lines (stretch target)
- **Minimum**: 80% lines (quality gate)

---

## Optional: Go Further (If Desired)

If you want to push coverage higher, here are concrete next steps:

### Option A: Identify Uncovered Code

Find the 98 remaining uncovered lines:
```bash
# Generate detailed report
gcovr build/staging --root . \
  --gcov-executable gcov-15 \
  --filter 'src/.*' \
  --html --html-details -o coverage_detailed.html

# Open and look for red lines
open coverage_detailed.html
```

**Then decide**:
1. **Error handling?** → Hard to test, acceptable to leave uncovered
2. **Dead code?** → Remove it
3. **Edge cases?** → Add targeted tests

### Option B: Improve Branch Coverage (37.8%)

Branch coverage measures how many conditional paths are tested:

1. **Find untested branches**:
   - Look in coverage report for yellow/orange conditional lines
   - These are if/else or case statements with untested paths

2. **Add targeted tests**:
   ```python
   # In toolchain/mfc/test/cases.py
   # Add tests that exercise specific conditional paths
   ```

3. **Focus on high-impact modules**:
   - Start with core simulation modules (m_time_steppers, m_riemann_solvers)
   - Then move to boundary conditions and specialized physics

### Option C: Expand the Test Suite

Add more test cases in `toolchain/mfc/test/cases.py`:

```python
# Example: Test more time integrators
def alter_time_integrators():
    for time_stepper in [1, 2, 3, 4, 5, 23]:
        cases.append(define_case_d(stack, f"time_stepper={time_stepper}",
            {'time_stepper': time_stepper, 't_step_stop': 5}))
```

**But remember**: Diminishing returns! 83.7% is already excellent.

---

## Files You Now Have

### Coverage Results
- `coverage_results_postprocess/` - Full results with `-a` flag
  - `index.html` - Interactive coverage report
  - `coverage.txt` - Text summary  
  - `tests.log` - Test execution details

### Scripts
- `run_postprocess_coverage.sh` - Automated coverage script

### Documentation
- `POSTPROCESS_COVERAGE_RESULTS.md` - Detailed analysis
- `COVERAGE_FINAL_SUMMARY.md` - Complete summary
- `NEXT_STEPS.md` - This file

---

## Maintenance Going Forward

### Regular Coverage Checks

Run coverage periodically:
```bash
# Quick check (10-20% of tests)
PERCENT=10 ./run_postprocess_coverage.sh

# Full check (100% of tests)
PERCENT=100 ./run_postprocess_coverage.sh
```

### When Adding New Code

1. **Write code**
2. **Add tests** (if needed)
3. **Run coverage** to ensure new code is tested:
   ```bash
   ./run_postprocess_coverage.sh
   ```
4. **Aim for**: New code should be at least 80% covered

### When Modifying Existing Code

1. **Make changes**
2. **Run tests**: `./mfc.sh test -a`
3. **Check coverage hasn't dropped**:
   - Before: 83.7%
   - After: Should be ≥ 83.7%

---

## Quick Commands Reference

### Run Full Coverage
```bash
./run_postprocess_coverage.sh
```

### Run Quick Coverage (10% of tests)
```bash
PERCENT=10 ./run_postprocess_coverage.sh
```

### View Coverage Report
```bash
open coverage_results_postprocess/index.html
```

### Check Current Coverage
```bash
cat coverage_results_postprocess/coverage.txt | grep "^lines:"
```

### Build Without Coverage (faster)
```bash
./mfc.sh build -t pre_process simulation post_process -j $(nproc)
```

### Run Tests Without Coverage (for development)
```bash
./mfc.sh test -j $(nproc)  # Fast, no post-processing
./mfc.sh test -a -j $(nproc)  # Complete, with post-processing
```

---

## Troubleshooting

### Coverage is 0% or Very Low
**Problem**: `.gcda` files not generated or not found  
**Solution**: 
```bash
# Make sure you build with --gcov
./mfc.sh build --gcov --no-gpu --debug -t pre_process simulation post_process

# Don't use --no-build when running tests
./mfc.sh test -a -j $(nproc)  # Correct
# NOT: ./mfc.sh test -a --no-build  # Skip this unless debugging
```

### gcovr Fails with Version Mismatch
**Problem**: `gcov` version doesn't match `gfortran`  
**Solution**:
```bash
# Use matching gcov version
which gcov-15  # Should exist if using gfortran-15
which gcov-14  # Or gfortran-14

# Update the script to use correct version
```

### Tests Take Too Long
**Problem**: Full coverage run takes 20-30 minutes  
**Solution**:
```bash
# Run subset during development
PERCENT=10 ./run_postprocess_coverage.sh  # 10% of tests, ~3-5 min

# Save full run for CI or pre-commit
PERCENT=100 ./run_postprocess_coverage.sh  # 100%, ~20-30 min
```

---

## Success Criteria Met ✅

- [x] **83.7% line coverage** - Excellent! ✅
- [x] **100% function coverage** - Perfect! ✅
- [x] **Post-processing tested** - Complete workflow! ✅
- [x] **Automated script** - `run_postprocess_coverage.sh` ready ✅
- [x] **Documentation** - Complete with examples ✅

---

## Final Recommendation

### For Now: **You're Done! 🎉**

83.7% coverage with 100% function coverage is **excellent** for a complex physics solver like MFC. You have:

1. ✅ Complete workflow testing (pre → sim → post)
2. ✅ High line coverage (83.7%)
3. ✅ Perfect function coverage (100%)
4. ✅ Automated tooling
5. ✅ Clear documentation

### Future Work (If Desired):

- **Branch coverage**: Currently 37.8%, could target 50%+
- **Edge cases**: Add tests for error conditions
- **Performance**: Optimize test runtime while maintaining coverage

**But honestly**: What you have now is production-ready and CI-ready! 🚀

---

**Generated**: November 3, 2025  
**Status**: Complete and production-ready ✅





