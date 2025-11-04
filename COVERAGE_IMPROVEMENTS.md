# MFC Code Coverage Improvement Plan - In Progress

## Current Status (Nov 1, 2025)

### Infrastructure Completed ✓
1. **Coverage tooling setup**
   - `toolchain/coverage.sh` - automated coverage collection and reporting
   - Auto-detects correct `gcov` version (gcov-15 for gfortran-15)
   - Proper `GCOV_PREFIX` configuration for `.gcda` file collection
   - Integration with `gcovr` for HTML, XML, and text reports

2. **Documentation created**
   - `docs/documentation/coverage.md` - comprehensive coverage guide
   - `README_COVERAGE.md` - quick start guide
   - `REGRESSION_TEST_EXPANSION.md` - detailed test expansion plan

3. **Build system updates**
   - CMake option `MFC_UNIT_TESTS` added
   - Support for pFUnit unit testing framework (in progress)
   - Coverage flags properly configured via `MFC_GCov`

### Running Now
- Baseline coverage assessment with 5% of test suite
- Command: `PERCENT=5 ./toolchain/coverage.sh`
- Output being logged to `build/coverage_full.log`

## Next Steps to Increase Coverage

### 1. Immediate: Add More Regression Test Variants

Based on analysis of `toolchain/mfc/test/cases.py`, we can add:

#### Time Integrators (High Impact)
```python
# Currently only testing time_stepper=3
# Add time_stepper=1,2,4,5,23
for ts in [1, 2, 4, 5, 23]:
    cases.append(define_case_d(stack, f"time_stepper={ts}", 
        {'time_stepper': ts, 't_step_stop': 5}))
```

#### WENO Schemes (Medium-High Impact)
```python
# Currently limited WENO testing
# Add weno_order=3,5,7 with various weno_eps values
for order in [3, 5, 7]:
    for eps in [1.0E-6, 1.0E-16, 1.0E-40]:
        cases.append(define_case_d(stack, f"weno_order={order}_eps={eps}", 
            {'weno_order': order, 'weno_eps': eps}))
```

#### Riemann Solvers (Medium Impact)
```python
# Test all riemann_solver options: 1-4
for rs in [1, 2, 3, 4]:
    cases.append(define_case_d(stack, f"riemann_solver={rs}", 
        {'riemann_solver': rs}))
```

#### Physics Combinations (High Impact)
```python
# Viscous flow variations
for visc_type in [1, 2, 3]:  # molecular, shear, bulk
    cases.append(define_case_d(stack, f"fluid_pp(1)%Re(1)=100_visc={visc_type}",
        {'fluid_pp(1)%Re(1)': 100.0, 'viscous': 'T', ...}))

# Surface tension
for st_model in [1, 2, 3]:
    cases.append(define_case_d(stack, f"surface_tension_model={st_model}",
        {'surface_tension': 'T', ...}))

# Phase change combinations
for pc_model in [1, 2, 3, 4, 5]:
    cases.append(define_case_d(stack, f"phase_change_model={pc_model}",
        {'phase_change': 'T', ...}))
```

#### Boundary Condition Combinations
```python
# Mix different BCs on different boundaries
bc_combos = [
    {'bc_x%beg': -1, 'bc_x%end': -2},  # periodic + reflective
    {'bc_x%beg': -3, 'bc_x%end': -1},  # ghost + periodic
    # ... more combinations
]
```

### 2. Unit Tests for Untested Modules

#### m_helper_basic.fpp
```fortran
! Test approximate equality functions
! Test default value checkers
! Test array utilities
```

#### m_precision_select.fpp
```fortran
! Verify precision settings
! Test integer/real/char default values
```

#### m_constants.fpp
```fortran
! Verify mathematical constants
! Test physical constants
```

### 3. Edge Cases in Existing Tests

#### Dimensional extremes
- Very small domains (m=5, n=5, p=5)
- Large aspect ratios (m=500, n=10, p=10)
- Single-cell dimensions (m=1, n=10, p=10)

#### Timestep variations
- Very small dt (t_step_stop=1000 with small dt)
- Adaptive timestep boundary cases
- CFL limit testing

### 4. Post-Process Test Expansion

Currently minimal post-process coverage. Add:
```python
# Parallel I/O tests
cases.append(define_case_d(stack, "parallel_io=T", 
    {'parallel_io': 'T', 'format': 1}))

# Different output formats
for fmt in [1, 2, 3]:  # Binary, ASCII, HDF5
    cases.append(define_case_d(stack, f"format={fmt}",
        {'format': fmt}))

# Slice output in all directions
for dim in ['x', 'y', 'z']:
    cases.append(define_case_d(stack, f"slice_{dim}",
        {f'{dim}_cb%beg': 10, f'{dim}_cb%end': 10}))
```

## Coverage Targets

### Phase 1 (Current - Baseline)
- **Goal**: Establish baseline coverage percentage
- **Method**: Run 5-10% of existing tests with coverage
- **Expected**: 40-60% line coverage

### Phase 2 (Regression Expansion)
- **Goal**: Reach 75% line coverage
- **Method**: Add 200-300 new regression test variants
- **Focus**: Time integrators, WENO, Riemann solvers

### Phase 3 (Unit Tests)
- **Goal**: Reach 85% line coverage
- **Method**: Add pFUnit unit tests for helper modules
- **Focus**: m_helper_basic, m_precision_select, m_constants

### Phase 4 (Edge Cases)
- **Goal**: Reach 90%+ line coverage
- **Method**: Add edge cases and error condition tests
- **Focus**: Boundary conditions, extreme parameters

## Metrics to Track

1. **Line Coverage** - % of lines executed
2. **Branch Coverage** - % of conditional branches taken
3. **Function Coverage** - % of functions called
4. **Module Coverage** - Per-module breakdown

## CI Integration

Coverage should:
- Run on every PR
- Generate diff coverage (only new/changed lines)
- Block merge if coverage drops >2%
- Upload reports to coverage service (Codecov/Coveralls)

## Files Modified

1. `/Users/spencer/Downloads/MFC/toolchain/coverage.sh` - Created
2. `/Users/spencer/Downloads/MFC/docs/documentation/coverage.md` - Created
3. `/Users/spencer/Downloads/MFC/CMakeLists.txt` - Modified (added MFC_UNIT_TESTS option)
4. `/Users/spencer/Downloads/MFC/tests/unit/` - Directory created (unit tests pending)

## Commands Reference

```bash
# Quick baseline (5% of tests)
PERCENT=5 MIN_LINES=10 ./toolchain/coverage.sh

# Full coverage (50% of tests)
PERCENT=50 MIN_LINES=65 MIN_BRANCHES=50 ./toolchain/coverage.sh

# View results
open build/coverage/index.html
cat build/coverage/summary.txt
```

---

**Last Updated**: Nov 1, 2025
**Status**: Infrastructure complete, baseline assessment in progress
**Next**: Add regression test variants based on REGRESSION_TEST_EXPANSION.md







