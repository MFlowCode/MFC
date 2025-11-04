# Code Coverage Strategy for MFC

## Overview

This document outlines the strategy for assessing and improving code coverage in MFC using a combination of unit tests and regression tests.

## Current Setup

### Coverage Infrastructure

MFC has GCC coverage support built into the CMake configuration:
- **Build flag**: `--gcov` enables coverage instrumentation
- **Compiler flags**: `-fprofile-arcs -ftest-coverage -O1` 
- **Link flags**: `-lgcov --coverage`
- **Supported compiler**: GNU gfortran (tested with GCC 15.1.0)

### Existing Tests

- **Regression tests**: `./mfc.sh test` runs parameterized test cases
- **Post-process tests**: `./mfc.sh test -a` includes HDF5/SILO validation
- **Test generator**: `toolchain/mfc/test/cases.py` creates test variants programmatically

## Quick Start: Assess Current Coverage

### Method 1: Use the Coverage Script (Recommended)

```bash
# Run with default settings (25% of tests, 65% line threshold, 50% branch threshold)
./toolchain/coverage.sh

# Run with custom settings
PERCENT=50 MIN_LINES=70 MIN_BRANCHES=60 ./toolchain/coverage.sh

# Run full test suite
PERCENT=100 ./toolchain/coverage.sh
```

### Method 2: Manual Assessment

```bash
# 1. Clean and build with coverage
./mfc.sh clean
./mfc.sh build --gcov --no-gpu --debug -t pre_process simulation post_process -j $(sysctl -n hw.ncpu)

# 2. Set environment for coverage collection
export GCOV_PREFIX=${PWD}/build/staging
export GCOV_PREFIX_STRIP=0

# 3. Run tests (without --no-build to collect coverage)
./mfc.sh test --no-examples -% 25 -j $(sysctl -n hw.ncpu)

# 4. Generate reports
mkdir -p build/coverage
gcovr build/staging --root . \
    --filter 'src/.*' \
    --html --html-details -o build/coverage/index.html \
    --xml-pretty -o build/coverage/coverage.xml \
    --print-summary

# 5. View results
open build/coverage/index.html  # macOS
# or
xdg-open build/coverage/index.html  # Linux
```

## Important Notes on Coverage Collection

### Why GCOV_PREFIX is Critical

When MFC binaries are installed to `build/install/`, they lose the direct path to the `.gcno` files in `build/staging/`. Setting `GCOV_PREFIX` ensures that:
1. Coverage data (`.gcda` files) are written back to the build directory
2. `.gcda` files are co-located with `.gcno` files for proper analysis

### Test Execution Strategy

**Current issue**: Running tests with `--no-build` flag (the default) means:
- Tests execute installed binaries from `build/install/`
- These binaries don't write `.gcda` files to the correct location
- Coverage analysis fails or shows 0% coverage

**Solution**: Either:
1. Use `GCOV_PREFIX` environment variable (implemented in `coverage.sh`)
2. Modify test infrastructure to support coverage mode
3. Run tests without `--no-build` (slower but ensures coverage collection)

## Identifying Under-Tested Code

### Using HTML Reports

1. Open `build/coverage/index.html`
2. Sort by "Lines Uncovered" column
3. Focus on files with:
   - High total line count
   - Low coverage percentage
   - Critical functionality

### Using Text Reports

```bash
# Generate sorted summary by uncovered lines
gcovr build/staging --root . --filter 'src/.*' --txt --sort-uncovered > build/coverage/sorted.txt
cat build/coverage/sorted.txt
```

### Per-Directory Analysis

```bash
# Check coverage by subsystem
gcovr build/staging --root . --filter 'src/common/.*' --print-summary
gcovr build/staging --root . --filter 'src/simulation/.*' --print-summary
gcovr build/staging --root . --filter 'src/pre_process/.*' --print-summary
gcovr build/staging --root . --filter 'src/post_process/.*' --print-summary
```

## Strategies to Improve Coverage

### 1. Unit Tests with pFUnit

**Target**: Pure/elemental functions and isolated logic in `src/common/`

**Setup** (to be implemented):
```cmake
# Add to CMakeLists.txt
option(MFC_UNIT_TESTS "Build unit tests" OFF)

if (MFC_UNIT_TESTS)
    find_package(PFUNIT REQUIRED)
    enable_testing()
    add_subdirectory(tests/unit)
endif()
```

**High-value targets**:
- `m_finite_differences.fpp`: Pure mathematical routines
- `m_variables_conversion.fpp`: Conversion functions
- `m_helper_basic.fpp`: Basic utilities
- `m_boundary_common.fpp`: Boundary indexing logic

**Example pFUnit test structure**:
```fortran
module test_m_finite_differences
    use pFUnit_mod
    use m_finite_differences
    use m_precision_select, only: wp
    implicit none

contains

    @test
    subroutine test_first_derivative_1d()
        real(wp), dimension(5) :: f, df_dx
        real(wp) :: dx
        integer :: i
        
        ! Setup: f(x) = x^2
        dx = 0.1_wp
        do i = 1, 5
            f(i) = real(i-1, wp)**2 * dx**2
        end do
        
        ! Execute
        call s_compute_fd_gradient_x(f, df_dx, dx, size(f))
        
        ! Assert: df/dx ≈ 2x
        @assertEqual(2.0_wp*dx, df_dx(2), tolerance=1.0e-6_wp)
    end subroutine

end module
```

### 2. Expand Regression Tests

**Target**: Untested code paths in physics modules

**Process**:
1. Review `build/coverage/index.html` to identify gaps
2. Edit `toolchain/mfc/test/cases.py`
3. Add targeted test variants
4. Generate golden files for new tests only

**High-impact additions**:

```python
# In toolchain/mfc/test/cases.py

# Add more time-stepping variants
def alter_time_steppers():
    for time_stepper in [1, 2, 3]:
        for cfl_mode in ['cfl_adap_dt', 'cfl_const_dt']:
            stack.push(f"time_stepper={time_stepper}, {cfl_mode}", {
                'time_stepper': time_stepper,
                cfl_mode: 'T',
                'cfl_target': 0.5,
                't_step_stop': 10
            })
            cases.append(define_case_d(stack, '', {}))
            stack.pop()

# Add rare boundary condition combinations
def alter_rare_bcs(dimInfo):
    for bc_combo in [(-13, -14), (-18, -19)]:  # Less-tested BC types
        stack.push(f"bc_combo={bc_combo}", {
            'bc_x%beg': bc_combo[0],
            'bc_x%end': bc_combo[1]
        })
        cases.append(define_case_d(stack, '', {}))
        stack.pop()
```

**Generate golden files**:
```bash
# After adding new cases
./mfc.sh test --generate -o "time_stepper" -j $(sysctl -n hw.ncpu)
```

### 3. Exclude Unreachable Code

For error-handling branches or boilerplate that shouldn't count against coverage:

```fortran
! Single line exclusion
if (error_condition) then
    call s_mpi_abort("Error message")  ! GCOVR_EXCL_LINE
end if

! Block exclusion
! GCOVR_EXCL_START
if (impossible_condition) then
    ! Defensive programming that should never execute
    call s_mpi_abort("This should never happen")
end if
! GCOVR_EXCL_STOP
```

## CI Integration Recommendations

### PR Coverage Check (Fast)

```bash
# In CI: Run on pull requests
./mfc.sh clean
./mfc.sh build --gcov --no-gpu --debug -t pre_process simulation post_process -j 4
export GCOV_PREFIX=${PWD}/build/staging
export GCOV_PREFIX_STRIP=0
./mfc.sh test --no-examples -% 25 -j 4
gcovr build/staging --root . --filter 'src/.*' \
    --fail-under-line 70 \
    --fail-under-branch 60 \
    --xml-pretty -o build/coverage/coverage.xml
```

### Nightly Full Coverage (Comprehensive)

```bash
# In CI: Run nightly or weekly
PERCENT=100 MIN_LINES=75 MIN_BRANCHES=65 ./toolchain/coverage.sh
# Upload build/coverage/index.html as artifact
# Upload build/coverage/coverage.xml to Codecov/Coveralls
```

### Diff Coverage (PR-specific)

```bash
# Requires diff-cover package
pip install diff-cover
gcovr build/staging --root . --filter 'src/.*' --xml-pretty -o coverage.xml
git fetch origin main
diff-cover coverage.xml --compare-branch origin/main --fail-under 80
```

## Tools and Dependencies

### Required
- **gfortran**: GCC 12+ (tested with GCC 15.1.0)
- **gcov**: Bundled with GCC (use matching version: `gcov-15` for `gfortran-15`)
- **gcovr**: `pip install gcovr` (tested with 8.3)

**Important**: You must use the gcov version that matches your gfortran compiler. For example:
```bash
# If using gfortran-15:
gcovr --gcov-executable gcov-15 ...

# Find the correct version:
which gcov-15 || which gcov-14 || which gcov
```

### Optional
- **pFUnit**: For unit tests (https://github.com/Goddard-Fortran-Ecosystem/pFUnit)
- **diff-cover**: For PR diff coverage (`pip install diff-cover`)
- **lcov**: Alternative to gcovr (`brew install lcov` or apt-get)

## Limitations and Considerations

### GPU Code

Coverage analysis only works for CPU code paths. GPU kernels (OpenACC/OpenMP) are not instrumented by gcov.

**Strategy**:
- Build CPU-only (`--no-gpu`) for coverage
- Rely on regression test numerics to validate GPU paths
- Unit test the host logic that drives GPU kernels

### Fypp-Generated Code

MFC uses Fypp preprocessing, which can complicate coverage analysis:
- Line numbers in `.f90` files may not match source `.fpp` files
- CMake configuration includes `--line-numbering` to help
- gcovr may need `--gcov-ignore-parse-errors` for complex macros

### Performance Impact

Coverage instrumentation adds overhead:
- **Build time**: ~10-20% slower with `-fprofile-arcs -ftest-coverage`
- **Runtime**: ~20-50% slower due to instrumentation
- **File I/O**: `.gcda` files written on every test run

**Mitigation**:
- Use `-O1` optimization (already configured) instead of `-O0`
- Run coverage checks on a subset of tests for fast feedback
- Reserve full coverage runs for nightly CI

## Next Steps

1. **Immediate** (Week 1):
   - Run `./toolchain/coverage.sh` to establish baseline
   - Review `build/coverage/index.html` and identify top-10 under-covered files
   - Document current coverage percentage

2. **Short-term** (Weeks 2-3):
   - Set up pFUnit infrastructure in `tests/unit/`
   - Add 5-10 unit tests for `src/common` pure functions
   - Add 10-20 targeted regression test variants for under-covered modules
   - Aim for +10-15% coverage increase

3. **Medium-term** (Month 2):
   - Integrate coverage checks into CI (PR and nightly)
   - Set and enforce coverage thresholds (start at 70% lines, 60% branches)
   - Create coverage dashboard/badge

4. **Long-term** (Ongoing):
   - Maintain coverage as new features are added
   - Increase thresholds gradually (target 80% lines, 70% branches)
   - Refactor large functions to improve testability

## Troubleshooting

### Problem: Coverage shows 0% despite tests running

**Cause**: Either `.gcda` files not written, or gcov version mismatch

**Solution**:
```bash
# 1. Check if .gcda files exist and have been updated
find build/staging -name "*.gcda" -exec ls -lh {} \; | head -10

# 2. Verify you're using the matching gcov version
gfortran --version  # Note the version number
which gcov-15       # Use gcov-XX matching gfortran-XX

# 3. Ensure GCOV_PREFIX is set when running tests
export GCOV_PREFIX=${PWD}/build/staging
export GCOV_PREFIX_STRIP=0

# 4. Run one test manually to verify
cd tests/MANUAL_COVERAGE_TEST
${PWD}/build/install/*/bin/pre_process
# Check if .gcda files are created/updated in build/staging

# 5. Try gcovr with the correct gcov executable
gcovr build/staging --root . --gcov-executable gcov-15 --filter 'src/.*' --print-summary
```

### Problem: gcovr fails with "no coverage data found"

**Cause**: Mismatch between `.gcno` and `.gcda` files or incorrect paths

**Solution**:
```bash
# Verify .gcno and .gcda are in the same directory
ls build/staging/*/CMakeFiles/*/dir/*/*.{gcno,gcda}

# Try running gcovr from the build staging directory
cd build/staging/XXXXXXXX  # Your build hash
gcovr . --root ${PWD}/../../.. --print-summary
```

### Problem: Line numbers don't match source files

**Cause**: Fypp preprocessing can alter line numbers

**Solution**:
- Fypp already configured with `--line-numbering` in CMake
- Use `--gcov-ignore-parse-errors` with gcovr if needed
- Focus on function-level coverage rather than line-by-line analysis

## References

- **GCC Coverage**: https://gcc.gnu.org/onlinedocs/gcc/Gcov.html
- **gcovr**: https://gcovr.com/en/stable/
- **pFUnit**: https://github.com/Goddard-Fortran-Ecosystem/pFUnit
- **MFC Testing Docs**: `docs/documentation/testing.md`
- **MFC Build System**: `CMakeLists.txt`, `toolchain/mfc/build.py`

## Contact

For questions or issues with coverage analysis:
1. Review this document and `docs/documentation/testing.md`
2. Check existing CI coverage reports (if available)
3. Open an issue with coverage report and logs

