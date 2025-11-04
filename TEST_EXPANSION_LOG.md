# MFC Test Suite Expansion Log

## Overview
This document tracks all test additions made to improve code coverage.

## Session: November 1, 2025

### Starting Point
- **Original test count**: 790 tests
- **Original coverage**: ~50-60% (estimated)

---

## Test Additions

### Round 1: Time Integrators (NEW FEATURE)
**Lines Modified**: `toolchain/mfc/test/cases.py:174-179`

**Function Added**:
```python
def alter_time_integrators():
    # Test different Runge-Kutta time integrators
    # time_stepper: 1=Euler, 2=RK2, 3=RK3 (default), 4=RK4, 5=RK5, 23=TVD-RK3
    for time_stepper in [1, 2, 4, 5, 23]:
        cases.append(define_case_d(stack, f"time_stepper={time_stepper}", 
            {'time_stepper': time_stepper, 't_step_stop': 5}))
```

**Impact**:
- Tests added: 15 (5 schemes × 3 dimensions)
- Code coverage target: `src/simulation/m_time_steppers.fpp`
- **Problem solved**: Previously had ZERO tests for non-default time steppers
- **Coverage gain**: Estimated +5-8%

---

### Round 2: Riemann Solver Expansion
**Lines Modified**: `toolchain/mfc/test/cases.py:198`

**Change**:
```python
# BEFORE:
for riemann_solver in [1, 5, 2]:

# AFTER:
for riemann_solver in [1, 5, 2, 3, 4]:
```

**Impact**:
- Tests added: ~571 (2 new solvers × all test combinations)
- Code coverage target: `src/simulation/m_riemann_solvers.fpp`
- **Problem solved**: Missing tests for solvers 3 and 4 (HLLD for MHD)
- **Coverage gain**: Estimated +3-5%
- **Test count after**: 1,376 tests

---

### Round 3: CFL Adaptation Modes (NEW FEATURE)
**Lines Modified**: `toolchain/mfc/test/cases.py:181-187`

**Function Added**:
```python
def alter_cfl_modes():
    # Test CFL adaptation and constant CFL modes
    cases.append(define_case_d(stack, "cfl_adap_dt=T", 
        {'cfl_adap_dt': 'T', 'cfl_target': 0.5, 't_step_stop': 10}))
    cases.append(define_case_d(stack, "cfl_const_dt=T", 
        {'cfl_const_dt': 'T', 'cfl_target': 0.3, 't_step_stop': 10}))
```

**Impact**:
- Tests added: 6 (2 modes × 3 dimensions)
- Code coverage target: `src/simulation/m_time_steppers.fpp` (CFL computation)
- **Problem solved**: Limited testing of adaptive/constant CFL modes
- **Coverage gain**: Estimated +2-3%

---

### Round 4: Model Equations (NEW FEATURE)
**Lines Modified**: `toolchain/mfc/test/cases.py:189-194`

**Function Added**:
```python
def alter_model_equations():
    # Test different model equation formulations
    # 1=gamma model, 2=pi-gamma model, 3=5-equation model
    for model_eqns in [1, 2, 3]:
        cases.append(define_case_d(stack, f"model_eqns={model_eqns}",
            {'model_eqns': model_eqns}))
```

**Impact**:
- Tests added: 9 (3 models × 3 dimensions)
- Code coverage target: Multiple files (equation handling throughout)
- **Problem solved**: Sparse testing of different equation models
- **Coverage gain**: Estimated +3-4%

---

### Round 5: Grid Stretching (NEW FEATURE)
**Lines Modified**: `toolchain/mfc/test/cases.py:196-201`

**Function Added**:
```python
def alter_grid_stretching():
    # Test grid stretching options (for non-uniform grids)
    cases.append(define_case_d(stack, "x_stretch=T",
        {'x_stretch': 'T', 'a_x': 1.5, 'x_a': -1.0, 'x_b': 1.0}))
    cases.append(define_case_d(stack, "loops_x=2",
        {'loops_x': 2}))
```

**Impact**:
- Tests added: 6 (2 grid options × 3 dimensions)
- Code coverage target: `src/pre_process/m_grid.fpp`, `src/simulation/m_start_up.fpp`
- **Problem solved**: Grid stretching was COMPLETELY UNTESTED
- **Coverage gain**: Estimated +2-3%

---

## Summary Statistics

### Test Count Growth
| Stage | Tests | Change |
|-------|-------|--------|
| Original | 790 | baseline |
| + Time integrators | 805 | +15 |
| + Riemann solvers | 1,376 | +571 |
| + CFL modes | 1,382 | +6 |
| + Model equations | 1,391 | +9 |
| + Grid stretching | 1,397 | +6 |
| **TOTAL** | **~1,397** | **+607 (+77%)** |

### Coverage Impact Summary
| Feature | Tests Added | Est. Coverage Gain | Priority |
|---------|-------------|-------------------|----------|
| Time integrators | 15 | +5-8% | HIGH ✓ |
| Riemann solvers | 571 | +3-5% | HIGH ✓ |
| CFL modes | 6 | +2-3% | MEDIUM ✓ |
| Model equations | 9 | +3-4% | MEDIUM ✓ |
| Grid stretching | 6 | +2-3% | MEDIUM ✓ |
| **TOTAL** | **607** | **+15-23%** | |

### Expected Coverage
- **Before**: ~50-60% line coverage
- **After**: ~65-75% line coverage
- **Improvement**: **+15-23 percentage points**

---

## Code Areas Now Covered

### 1. Time Integration (`m_time_steppers.fpp`)
- ✅ Euler (RK1) 
- ✅ RK2
- ✅ RK3 (was already tested as default)
- ✅ RK4
- ✅ RK5
- ✅ TVD-RK3

### 2. Riemann Solvers (`m_riemann_solvers.fpp`)
- ✅ Solver 1 (HLL)
- ✅ Solver 2 (HLLC)
- ✅ Solver 3 (NEW)
- ✅ Solver 4 (HLLD for MHD) (NEW)
- ✅ Solver 5 (Viscous)

### 3. CFL Computation (`m_time_steppers.fpp`)
- ✅ Fixed dt (was already tested)
- ✅ Adaptive CFL (NEW)
- ✅ Constant CFL (NEW)

### 4. Model Equations (Multiple files)
- ✅ Gamma model (model_eqns=1) (NEW)
- ✅ Pi-gamma model (model_eqns=2) (NEW)
- ✅ 5-equation model (model_eqns=3) (NEW)

### 5. Grid Generation (`m_grid.fpp`)
- ✅ Uniform grids (was already tested)
- ✅ Stretched grids (NEW)
- ✅ Multiple grid loops (NEW)

---

## Still Available for Future Expansion

### High Priority (Not Yet Done)
1. **Post-process output variations** (+8-12% coverage)
   - Parallel I/O options
   - Different file formats (Binary, ASCII, HDF5, Silo)
   - Slice outputs in all directions
   - Estimated: 20-40 new tests

2. **Physics combinations** (+10-15% coverage)
   - Viscous + bubbles
   - Surface tension variations
   - Phase change models
   - Hypoelasticity options
   - Estimated: 100-200 new tests

3. **Boundary condition combinations** (+5-8% coverage)
   - Mixed BCs on different boundaries
   - Ghost cell + wall combinations
   - Periodic + non-periodic
   - Estimated: 50-100 new tests

### Medium Priority (Not Yet Done)
4. **Dimensional edge cases** (+3-5% coverage)
   - Very small domains (m=5, n=5, p=5)
   - Large aspect ratios
   - Single-cell dimensions
   - Estimated: 30-50 new tests

5. **Monopole source variations** (+2-3% coverage)
   - Different pulse types
   - Multiple sources
   - Various support values
   - Estimated: 20-30 new tests

---

## Integration with Main Test Loop

All new test functions are called in `foreach_dimension()`:

```python
def foreach_dimension():
    for dimInfo, dimParams in get_dimensions():
        stack.push(f"{len(dimInfo[0])}D", dimParams)
        alter_bcs(dimInfo)
        alter_grcbc(dimInfo)
        alter_weno(dimInfo)
        alter_time_integrators()      # NEW
        alter_cfl_modes()              # NEW
        alter_model_equations()        # NEW
        alter_grid_stretching()        # NEW
        alter_muscl()
        alter_num_fluids(dimInfo)
        if len(dimInfo[0]) == 2:
            alter_2d()
        if len(dimInfo[0]) == 3:
            alter_3d()
        alter_lag_bubbles(dimInfo)
        # ... rest of loop
```

Each function is called for 1D, 2D, and 3D test variants, multiplying the impact.

---

## Verification Commands

### Count Total Tests
```bash
./mfc.sh test -l | wc -l
# Expected: ~1,397 (up from 790)
```

### List New Time Integrator Tests
```bash
./mfc.sh test -l | grep time_stepper
# Should show: time_stepper=1, 2, 4, 5, 23
```

### List New Riemann Solver Tests
```bash
./mfc.sh test -l | grep "riemann_solver=3\|riemann_solver=4"
# Should show many results
```

### List New CFL Tests
```bash
./mfc.sh test -l | grep "cfl_adap_dt\|cfl_const_dt"
# Should show adaptive and constant CFL tests
```

### List New Model Equation Tests
```bash
./mfc.sh test -l | grep model_eqns
# Should show model_eqns=1, 2, 3
```

### List New Grid Stretching Tests
```bash
./mfc.sh test -l | grep "x_stretch\|loops_x"
# Should show grid stretching variations
```

---

## Coverage Analysis Status

### Current Run
- **Started**: Nov 1, 2025, 2:24 PM
- **Command**: `PERCENT=50 MIN_LINES=50 MIN_BRANCHES=30 ./toolchain/coverage.sh`
- **Test percentage**: 50% of 1,397 tests (~698 tests)
- **Expected duration**: 30-60 minutes
- **Output**: `build/coverage/index.html` (HTML), `build/coverage/summary.txt` (text)

### What to Expect
The coverage report will show:
- **Line coverage**: Should be 65-75% (up from ~50-60%)
- **Branch coverage**: Should be 45-55%
- **Function coverage**: Should be 70-80%

Significant improvements expected in:
- `src/simulation/m_time_steppers.fpp` (time integration)
- `src/simulation/m_riemann_solvers.fpp` (Riemann solvers)
- `src/pre_process/m_grid.fpp` (grid stretching)
- `src/simulation/m_start_up.fpp` (initialization with new options)

---

## Next Steps After This Run

1. **Analyze the HTML report**: `open build/coverage/index.html`
   - Identify files with lowest coverage
   - Look for red (uncovered) lines in key files

2. **Add post-process tests** (high impact, minimal current coverage)
   - Follow `REGRESSION_TEST_EXPANSION.md` lines 350-410

3. **Add physics combination tests** (medium-high impact)
   - Follow `REGRESSION_TEST_EXPANSION.md` lines 180-280

4. **Target 90% coverage**:
   - Run full coverage again: `PERCENT=100 ./toolchain/coverage.sh`
   - Add unit tests for helper modules
   - Cover remaining edge cases

---

**Log Updated**: Nov 1, 2025, 2:35 PM
**Status**: 607 new tests added (+77%), coverage run in progress
**Files Modified**: `toolchain/mfc/test/cases.py` (5 new test functions)






