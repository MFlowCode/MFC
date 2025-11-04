# Full Coverage Run - Status

## ✅ **SUCCESSFULLY STARTED!**

The full coverage comparison run has been successfully launched!

## Current Status

- **Script**: `run_full_coverage.sh`
- **Started**: Saturday, November 1, 2025, ~4:20 PM EDT
- **Current Phase**: [2/4] Building with coverage instrumentation
- **Progress**: Compiling Fortran sources (syscheck, pre_process targets)

## Timeline Estimate

| Phase | Task | Duration | Status |
|-------|------|----------|--------|
| 1 | Clean previous builds | ~1 min | ✅ COMPLETE |
| 2 | Build with coverage | ~10-15 min | 🔄 IN PROGRESS |
| 3 | Run 100% of test suite (1,397 tests) | ~2-4 hours | ⏳ PENDING |
| 4 | Generate coverage reports | ~5-10 min | ⏳ PENDING |

**Total estimated time**: 2.5 - 5 hours

## Monitoring

The coverage run is logging all output to:
```bash
build/full_coverage_run.log
```

### View Progress in Real-Time:
```bash
tail -f build/full_coverage_run.log
```

### Check if Still Running:
```bash
ps aux | grep "run_full_coverage.sh" | grep -v grep
```

### View Last 50 Lines:
```bash
tail -50 build/full_coverage_run.log
```

## What's Happening Now

The script is currently:
1. ✅ Cleaned previous build artifacts
2. ✅ Set up Python virtual environment
3. ✅ Installing Python dependencies (numpy, pandas, etc.)
4. 🔄 Building MFC with GCC coverage flags (`--gcov --no-gpu --debug`)
5. 🔄 Compiling Fortran modules:
   - syscheck ✅ DONE
   - pre_process 🔄 IN PROGRESS
   - simulation ⏳ PENDING
   - post_process ⏳ PENDING

## What Happens Next

Once the build completes (~10 more minutes), the script will:

1. **Run ALL Tests** (100% of 1,397 tests)
   - This includes the original 790 tests
   - Plus the 607 new tests we added for:
     - Time integrators (Euler, RK2, RK4, RK5, TVD-RK3)
     - Riemann solvers (1, 2, 3, 4, 5)
     - CFL modes (adaptive and constant)
     - Model equations (gamma, pi-gamma, 5-equation)
     - Grid stretching configurations
   - Each test runs the full workflow: pre_process → simulation → post_process
   - Uses `--no-build` flag so installed binaries are used
   - Runs with max parallelism (10 jobs)

2. **Collect Coverage Data**
   - `.gcda` files will be generated in `build/staging/` directories
   - These contain execution counts for each line of Fortran code

3. **Generate Reports**
   - Uses `gcov-15` (matches `gfortran-15`)
   - Produces summary text file
   - Shows line and branch coverage percentages

## Output Files

When complete, you'll find:

- `build/full_coverage_run.log` - Complete log of the entire run
- `build/coverage_full/summary.txt` - Final coverage statistics
- `.gcda` files in `build/staging/` - Raw coverage data

## Coverage Expectations

### Baseline (Before Expansion)
- **Test count**: 790 tests
- **Expected coverage**: ~40-60% (estimate)

### New Suite (After Expansion)
- **Test count**: 1,397 tests (+607, +77%)
- **Expected coverage**: ~50-70% (estimate, improvement of 10-15%)

### Why the Improvement?

The 607 new tests target previously untested code paths:
- **Time stepper variations** (tests 5 different integrators vs just RK3)
- **Riemann solver variants** (tests solvers 3 and 4, previously skipped)
- **CFL adaptivity** (tests adaptive and constant CFL modes)
- **Model equations** (tests gamma, pi-gamma, and 5-equation models)
- **Grid stretching** (tests non-uniform grid configurations)

## Troubleshooting

### If the Process Dies

Check the log:
```bash
tail -100 build/full_coverage_run.log | grep -i error
```

### If It Seems Stuck

The test phase can take 2-4 hours. Check if tests are running:
```bash
tail -20 build/full_coverage_run.log
```

You should see output like:
```
[✓✓✓✓✓✗✓✓] (8/1397) Running tests...
```

### Manual Cleanup

If you need to stop and restart:
```bash
# Kill the process
pkill -9 -f "run_full_coverage.sh"

# Clean up
./mfc.sh clean

# Restart
./run_full_coverage.sh &
```

## Next Steps

Once the run completes (check the log for "FULL COVERAGE RUN COMPLETE"), you can:

1. **View the summary**:
   ```bash
   cat build/coverage_full/summary.txt
   ```

2. **Compare with baseline** (if we ran both):
   - Compare line coverage percentages
   - Compare branch coverage percentages
   - Identify which modules improved most

3. **Generate HTML reports** (optional):
   ```bash
   PERCENT=100 ./toolchain/coverage_fixed.sh
   open build/coverage/index.html
   ```

---

**Last Updated**: Sat Nov 1, 2025, 4:21 PM EDT
**Status**: 🔄 BUILD IN PROGRESS





