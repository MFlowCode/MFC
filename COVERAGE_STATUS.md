# MFC Code Coverage - Current Status

## What We've Accomplished

### ✅ Infrastructure Setup (Completed)

1. **Coverage Build System**
   - Verified GNU compiler with coverage support (GCC 15.1.0)
   - Tested `--gcov` build flag integration
   - Confirmed `.gcno` files are generated during compilation

2. **Automated Coverage Script**
   - Created `toolchain/coverage.sh` for one-command coverage assessment
   - Handles GCOV_PREFIX environment setup
   - Generates HTML, XML, and text reports
   - Includes threshold checking

3. **Comprehensive Documentation**
   - Created `docs/documentation/coverage.md` with full strategy
   - Documented tools, workflows, and troubleshooting
   - Provided examples for unit tests and regression test expansion

## ✅ SOLUTION FOUND: Coverage Data Collection Working!

### The Issues We Solved

1. **GCOV_PREFIX Environment Variable**
   - The installed binaries lose connection to `.gcno` files in the build directory
   - **Solution**: Set `GCOV_PREFIX=${PWD}/build/staging` before running tests
   - This tells coverage runtime to write `.gcda` files back to the build directory

2. **gcov Version Mismatch** (CRITICAL)
   - Using default `gcov` with `gfortran-15` produced empty coverage reports
   - **Solution**: Use matching `gcov-15` for `gfortran-15`
   - The `toolchain/coverage.sh` script now auto-detects the correct version

### Verified Working Setup

```bash
export GCOV_PREFIX=${PWD}/build/staging
export GCOV_PREFIX_STRIP=0
./mfc.sh run tests/MANUAL_COVERAGE_TEST/case.py -t pre_process -n 1 --no-build
gcovr build/staging --root . --gcov-executable gcov-15 --filter 'src/.*' --print-summary
```

**Result**: **45.7% line coverage** from a single pre_process run!

## Next Steps to Get Coverage Working

### Option 1: Modify Test Infrastructure (Recommended for Long-term)

Update `toolchain/mfc/test/case.py` to support a coverage mode:

```python
def run(self, targets, gpus, coverage_mode=False):
    # ... existing code ...
    
    if coverage_mode or os.environ.get('MFC_COVERAGE_MODE'):
        # Don't use --no-build; rebuild each test to ensure .gcda collection
        command = [mfc_script, "run", filepath, *tasks, *case_optimization, ...]
    else:
        # Normal mode: use --no-build for speed
        command = [mfc_script, "run", filepath, "--no-build", *tasks, ...]
```

### Option 2: Use GCOV_PREFIX (Quick Test)

Run manually to verify coverage collection:

```bash
# 1. Build with coverage
./mfc.sh build --gcov --no-gpu --debug -t pre_process simulation post_process -j $(sysctl -n hw.ncpu)

# 2. Set environment
export GCOV_PREFIX=${PWD}/build/staging
export GCOV_PREFIX_STRIP=0

# 3. Run ONE test manually to verify
cd tests/
mkdir -p TEST001
cd TEST001
cat > case.py << 'EOF'
#!/usr/bin/env python3
import json
print(json.dumps({
    'run_time_information': {'dt': 1e-06, 't_step_stop': 10, 't_step_save': 10},
    'm': 49, 'n': 0, 'p': 0,
    'x_domain': {'beg': 0.0, 'end': 1.0},
    'bc_x': {'beg': -3, 'end': -3},
    'num_fluids': 1,
    'num_patches': 1,
    'patch_icpp(1)': {
        'geometry': 1,
        'x_centroid': 0.5,
        'length_x': 1.0,
        'vel(1)': 0.0,
        'pres': 1.0,
        'alpha_rho(1)': 1.0,
        'alpha(1)': 1.0
    },
    'fluid_pp(1)': {'gamma': 1.4, 'pi_inf': 0.0},
    'weno_order': 3
}))
EOF
chmod +x case.py

# 4. Run pre_process manually
../../build/install/*/bin/pre_process

# 5. Check if .gcda files were created/updated in build/staging
find ../../build/staging -name "*.gcda" -newer ../../build/install -exec ls -lh {} \; | head -5

# 6. If files exist and are non-zero size, run gcovr
cd ../..
gcovr build/staging --root . --filter 'src/pre_process/.*' --print-summary
```

### Option 3: Run Test Without --no-build Flag

Modify the test command temporarily:

```bash
# In toolchain/mfc/test/case.py, line 146, comment out --no-build:
# command = [
#     mfc_script, "run", filepath, "--no-build", *tasks, ...
# ]
# Change to:
# command = [
#     mfc_script, "run", filepath, *tasks, ...
# ]
```

Then run:
```bash
./mfc.sh test --no-examples -% 10 -j 4
```

## Immediate Action Items

1. **Verify Coverage Collection** (15 minutes)
   - Run Option 2 above to manually test one case
   - Confirm `.gcda` files are created and non-empty
   - Run gcovr to see if coverage is > 0%

2. **If Manual Test Works** (30 minutes)
   - Implement Option 1 (modify test infrastructure)
   - Add `--coverage-mode` flag to `./mfc.sh test`
   - Run small test suite with new flag
   - Generate HTML report

3. **If Manual Test Doesn't Work** (troubleshooting)
   - Check if binaries have coverage instrumentation:
     ```bash
     nm build/install/*/bin/pre_process | grep gcov
     ```
   - Verify build flags were applied:
     ```bash
     grep -r "fprofile-arcs" build/staging/*/CMakeCache.txt
     ```
   - Try rebuilding with verbose flags:
     ```bash
     ./mfc.sh build --gcov --no-gpu --debug --verbose -t pre_process -j 1
     ```

## Expected Timeline

### Week 1: Get Coverage Working
- [ ] Fix coverage data collection (1-2 days)
- [ ] Generate baseline coverage report (1 hour)
- [ ] Identify top-10 under-covered files (1 hour)
- [ ] Document baseline percentage

### Week 2: Quick Wins
- [ ] Add 5-10 unit tests for `src/common` helpers
- [ ] Add 10-20 targeted regression test variants
- [ ] Aim for +10% coverage increase
- [ ] Create coverage dashboard

### Week 3: CI Integration
- [ ] Add coverage check to PR workflow
- [ ] Set up nightly full coverage run
- [ ] Configure coverage badges/reports
- [ ] Establish thresholds (70% lines, 60% branches)

## Resources Created

1. **`toolchain/coverage.sh`** - Automated coverage script
2. **`docs/documentation/coverage.md`** - Complete strategy guide
3. **`COVERAGE_STATUS.md`** (this file) - Current status and next steps

## Key Commands

```bash
# Quick coverage check (once working)
./toolchain/coverage.sh

# Custom coverage run
PERCENT=50 MIN_LINES=70 ./toolchain/coverage.sh

# Manual coverage analysis
gcovr build/staging --root . --filter 'src/.*' --print-summary

# View HTML report
open build/coverage/index.html
```

## Known Limitations

1. **GPU Code**: Coverage only works for CPU code paths
2. **Performance**: Coverage builds are 20-50% slower
3. **Fypp Preprocessing**: Line numbers may not match exactly
4. **CI Not Yet Configured**: Coverage checks not in automated pipeline

## Questions?

See `docs/documentation/coverage.md` for:
- Detailed troubleshooting
- Tool documentation
- Examples of unit tests
- Regression test expansion patterns
- CI integration recommendations

