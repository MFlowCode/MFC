# Coverage Comparison Run - ACTIVE

**Status**: ✅ Running  
**Started**: Monday, November 3, 2025, 12:34 AM EST  
**Process ID**: 25547  
**Current Phase**: Phase 1 - Building baseline with coverage  
**Estimated Total Time**: 4-6 hours

---

## What's Running

A comprehensive coverage comparison that will:

1. **Phase 1**: Run baseline test suite (528 tests) with post-processing and generate coverage report
2. **Phase 2**: Add new safe tests (time integrators, CFL adaptive, Riemann solver 3)
3. **Phase 3**: Run expanded test suite with post-processing and generate coverage report
4. **Phase 4**: Generate comparison report with runtime and coverage analysis

---

## Monitoring

### Progress Log
```bash
tail -f coverage_results/progress.log
```

### Full Output
```bash
tail -f /tmp/coverage_run.log
```

### Check if Running
```bash
ps -p 25547
```

---

## Expected Timeline

| Phase | Task | Duration | ETA |
|-------|------|----------|-----|
| 1 | Baseline build | ~3 min | 12:38 AM |
| 1 | Baseline tests (528 tests) | ~1-2 hours | 1:30-2:30 AM |
| 1 | Generate baseline coverage | ~5 min | 1:35-2:35 AM |
| 2 | Modify test suite | ~1 min | 1:36-2:36 AM |
| 3 | Expanded build | ~3 min | 1:39-2:39 AM |
| 3 | Expanded tests (~550-600 tests) | ~1-2 hours | 2:40-4:40 AM |
| 3 | Generate expanded coverage | ~5 min | 2:45-4:45 AM |
| 4 | Create comparison report | ~1 min | 2:46-4:46 AM |
| **TOTAL** | **All phases** | **~2.5-4.5 hours** | **~3-5 AM** |

---

## Output Files

When complete, you'll find in `coverage_results/`:

- **`FINAL_REPORT.md`** ⭐ Main results with comparison
- `progress.log` - Timeline of execution
- `baseline_build.log` - Baseline build output
- `baseline_tests.log` - Baseline test execution log
- `baseline_coverage.txt` - Baseline coverage report
- `expanded_build.log` - Expanded build output
- `expanded_tests.log` - Expanded test execution log
- `expanded_coverage.txt` - Expanded coverage report

---

## Current Status Details

### What Just Happened
- ✅ Script launched successfully
- ✅ Cleaning old build
- 🔄 Building MFC with GCC coverage instrumentation

### What's Next
- Build will complete in ~3 minutes
- Then run 528 baseline tests (includes post-processing!)
- Each test runs: syscheck → pre_process → simulation → post_process
- Coverage data (.gcda files) collected throughout

---

## Key Improvements in This Run

1. **Post-processing included** - Previous runs used `--no-examples`, this doesn't
2. **Safe test additions** - Only adding tests known to work
3. **Direct logging** - No buffering issues
4. **Outside build/** - Logs won't be deleted by `./mfc.sh clean`
5. **Proper test counting** - Fixed UUID pattern matching

---

## Tests Being Added (Phase 2)

- **Time integrators**: RK2, RK4, RK5 (avoiding problematic ones)
- **Adaptive CFL**: cfl_adap_dt=T
- **Riemann solver 3**: Exact Riemann solver

Expected to add ~20-70 new tests depending on dimensions.

---

**Last Updated**: 12:35 AM EST  
**Next Check**: 12:45 AM (see if build complete and tests starting)





