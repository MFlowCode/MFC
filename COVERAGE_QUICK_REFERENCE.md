# MFC Coverage - Quick Reference Card

## 📊 Current Status

- **Line Coverage**: 83.7% (504/602 lines)
- **Function Coverage**: 100% (15/15 functions)
- **Branch Coverage**: 37.8% (1943/5146 branches)
- **Test Count**: 528 tests
- **Flag Required**: `-a` (for post-processing validation)

---

## 🚀 Quick Commands

### Run Full Coverage (Recommended for CI)
```bash
./run_postprocess_coverage.sh
```
*Takes ~20-30 minutes, tests complete workflow*

### Run Quick Coverage (10% of tests)
```bash
PERCENT=10 ./run_postprocess_coverage.sh
```
*Takes ~3-5 minutes, good for quick checks*

### View Results
```bash
open coverage_results_postprocess/index.html
```

---

## 📝 Manual Coverage Run

```bash
# 1. Clean and build with coverage
./mfc.sh clean
./mfc.sh build --gcov --no-gpu --debug \
  -t pre_process simulation post_process \
  -j $(sysctl -n hw.ncpu)

# 2. Run tests WITH -a flag (essential!)
./mfc.sh test -a -j $(sysctl -n hw.ncpu)

# 3. Generate report
gcovr build/staging --root . \
  --gcov-executable gcov-15 \
  --filter 'src/.*' \
  --html --html-details -o coverage.html \
  --print-summary
```

---

## ⚠️ Critical Points

### ALWAYS Use `-a` Flag
```bash
./mfc.sh test -a  # ✅ Correct - 83.7% coverage
./mfc.sh test     # ❌ Wrong - Only 62.1% coverage
```

The `-a` flag enables post-processing validation, which is **essential** for complete coverage.

### Use Matching gcov Version
```bash
# If using gfortran-15:
--gcov-executable gcov-15  # ✅ Correct

# Check which you have:
which gfortran    # Shows version
which gcov-15     # Should exist
```

---

## 📂 Generated Files

After running coverage, you'll find:

```
coverage_results_postprocess/
├── index.html        # ← Open this for visual report
├── coverage.txt      # Text summary
├── tests.log         # Test execution details
├── build.log         # Build log
└── progress.log      # Timeline
```

---

## 🎯 Coverage Targets

| Metric | Current | Target | Status |
|--------|---------|--------|--------|
| Lines | 83.7% | ≥80% | ✅ PASS |
| Functions | 100% | ≥90% | ✅ PASS |
| Branches | 37.8% | ≥40% | ⚠️ Close |

---

## 🔧 Common Issues

### Issue: Coverage shows 0%
```bash
# Solution: Don't use --no-build flag
./mfc.sh test -a -j $(nproc)  # ✅ Correct
```

### Issue: "Version mismatch" error
```bash
# Solution: Use matching gcov version
--gcov-executable gcov-15
```

### Issue: Tests fail
```bash
# Check logs:
cat coverage_results_postprocess/tests.log
```

---

## 📚 Documentation Files

- `COVERAGE_FINAL_SUMMARY.md` - Complete analysis
- `POSTPROCESS_COVERAGE_RESULTS.md` - Detailed results
- `NEXT_STEPS.md` - What to do next
- `COVERAGE_QUICK_REFERENCE.md` - This file

---

## ✅ Success Checklist

For a proper coverage run:

- [ ] Built with `--gcov --no-gpu --debug`
- [ ] All targets included: `pre_process simulation post_process`
- [ ] Used `-a` flag when running tests
- [ ] Used `gcov-15` (matching compiler version)
- [ ] Generated HTML report
- [ ] Reviewed uncovered lines

---

## 💡 Pro Tips

### Fast Iteration During Development
```bash
# Skip coverage for quick tests:
./mfc.sh build -t pre_process simulation -j $(nproc)
./mfc.sh test -j $(nproc)
```

### Full Validation Before Commit
```bash
# Run complete coverage:
./run_postprocess_coverage.sh
# Check coverage hasn't dropped
```

### CI Configuration
```bash
# In your CI pipeline:
./mfc.sh build --gcov --no-gpu --debug -t pre_process simulation post_process
./mfc.sh test -a -j $(nproc)
gcovr build/staging --root . --gcov-executable gcov-15 --filter 'src/.*' --xml -o coverage.xml
```

---

## 📞 Quick Answers

**Q: What's the minimum acceptable coverage?**  
A: 80% lines, 90% functions (you're at 83.7% and 100%!)

**Q: Why is the -a flag important?**  
A: It runs post-processing validation. Without it, you only test 2/3 of the workflow.

**Q: How long does it take?**  
A: ~20-30 minutes for 100% of tests, ~3-5 minutes for 10%

**Q: Is 83.7% good enough?**  
A: YES! This is excellent for a complex physics solver.

**Q: Should I aim for 100% coverage?**  
A: No. Diminishing returns. 80-90% is the sweet spot.

---

**Last Updated**: November 3, 2025  
**Status**: Production Ready ✅





