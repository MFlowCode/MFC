# Regression Test Expansion Plan

Based on coverage analysis, here are targeted regression test additions to improve coverage of under-tested modules.

## Priority 1: Time Stepping & CFL Modes

**Target modules**: `m_time_steppers.fpp`, `m_rhs.fpp`

**Current gap**: Only `time_stepper=3` is tested by default. CFL adaptive and constant modes are rarely exercised.

**Add to `toolchain/mfc/test/cases.py`**:

```python
def alter_time_steppers(dimInfo):
    """Test different time stepping schemes and CFL modes"""
    # Test different time steppers
    for time_stepper in [1, 2, 3]:
        stack.push(f"time_stepper={time_stepper}", {
            'time_stepper': time_stepper,
            't_step_stop': 20,
            't_step_save': 20,
            'dt': 1e-6
        })
        cases.append(define_case_d(stack, '', {}))
        stack.pop()
    
    # Test CFL adaptive time stepping
    stack.push('CFL Adaptive', {
        'cfl_adap_dt': 'T',
        'cfl_target': 0.5,
        'cfl_max': 0.8,
        'cfl_min': 0.01,
        't_step_stop': 20,
        't_step_save': 20,
        'n_start': 0
    })
    cases.append(define_case_d(stack, '', {}))
    stack.pop()
    
    # Test CFL constant mode
    stack.push('CFL Constant', {
        'cfl_const_dt': 'T',
        'cfl_target': 0.3,
        't_step_stop': 20,
        't_step_save': 20,
        'n_start': 0
    })
    cases.append(define_case_d(stack, '', {}))
    stack.pop()
```

**Integration point**: Add call in `foreach_dimension()` after `alter_num_fluids()`.

---

## Priority 2: Rarely-Tested Boundary Conditions

**Target modules**: `m_cbc.fpp`, `m_compute_cbc.fpp`, `m_boundary_common.fpp`

**Current gap**: Only common BCs (-1 through -12) are well-tested. Rare BCs need coverage.

**Add to `toolchain/mfc/test/cases.py`**:

```python
def alter_rare_bcs(dimInfo):
    """Test less-common boundary conditions"""
    # BC types that are under-tested
    rare_bcs = [
        (-13, "Riemann extrapolation"),
        (-14, "Characteristic decomposition"),
        (-18, "Far-field"),
        (-19, "Custom BC"),
    ]
    
    for bc_type, bc_name in rare_bcs:
        # Only test in appropriate dimensions
        if len(dimInfo[0]) >= 1:
            stack.push(f"bc={bc_type} ({bc_name})", get_bc_mods(bc_type, dimInfo))
            cases.append(define_case_d(stack, '', {}))
            stack.pop()
```

**Integration point**: Add call in `foreach_dimension()` alongside existing `alter_bcs()`.

---

## Priority 3: Viscous Flux Computation Variants

**Target modules**: `m_viscous.fpp`

**Current gap**: Limited testing of `weno_avg` and different viscous flux computations.

**Add to `toolchain/mfc/test/cases.py`**:

```python
def alter_viscous_variants():
    """Test viscous flux computation variants"""
    # Already tested in some cases, but add focused variants
    stack.push("Viscous Minimal", {
        'viscous': 'T',
        'fluid_pp(1)%Re(1)': 100.0,
        'dt': 1e-9,
        't_step_stop': 10,
        't_step_save': 10
    })
    
    # Test weno_avg without weno_Re_flux
    stack.push("weno_avg only", {'weno_avg': 'T', 'weno_Re_flux': 'F'})
    cases.append(define_case_d(stack, '', {}))
    stack.pop()
    
    # Test different finite difference orders
    for fd_order in [1, 2, 4]:
        stack.push(f"fd_order={fd_order}", {'fd_order': fd_order})
        cases.append(define_case_d(stack, '', {}))
        stack.pop()
    
    stack.pop()
```

**Integration point**: Add in `foreach_dimension()` within 1-fluid section.

---

## Priority 4: Riemann Solver Edge Cases

**Target modules**: `m_riemann_solvers.fpp`

**Current gap**: Not all solver + state combinations tested.

**Add to `toolchain/mfc/test/cases.py`**:

```python
def alter_riemann_variants():
    """Test Riemann solver edge cases"""
    # Riemann solver + wave speed combinations
    for solver in [1, 2, 5]:
        for wave_speeds in [1, 2]:
            if solver in [1, 2]:  # Only applicable to HLL/HLLC
                stack.push(f"riemann={solver}, wave_speeds={wave_speeds}", {
                    'riemann_solver': solver,
                    'wave_speeds': wave_speeds,
                    't_step_stop': 15,
                    't_step_save': 15
                })
                cases.append(define_case_d(stack, '', {}))
                stack.pop()
    
    # Test average state variants
    for avg_state in [1, 2]:
        stack.push(f"avg_state={avg_state}", {
            'riemann_solver': 1,
            'avg_state': avg_state,
            't_step_stop': 15,
            't_step_save': 15
        })
        cases.append(define_case_d(stack, '', {}))
        stack.pop()
```

**Integration point**: Add in `foreach_dimension()` as standalone function.

---

## Priority 5: Data Output Variants

**Target modules**: `m_data_output.fpp`

**Current gap**: Limited testing of different output formats and precision modes.

**Add to `toolchain/mfc/test/cases.py`**:

```python
def alter_output_formats():
    """Test different output formats and precision"""
    # Test precision modes
    for precision in [1, 2]:
        stack.push(f"precision={precision}", {
            'precision': precision,
            't_step_stop': 10,
            't_step_save': 10,
            'parallel_io': 'T'
        })
        cases.append(define_case_d(stack, '', {}))
        stack.pop()
    
    # Test format modes
    for format_mode in [1, 2]:
        stack.push(f"format={format_mode}", {
            'format': format_mode,
            't_step_stop': 10,
            't_step_save': 10
        })
        cases.append(define_case_d(stack, '', {}))
        stack.pop()
```

**Integration point**: Add in `foreach_dimension()` near the end.

---

## Priority 6: MUSCL Reconstruction Variants

**Target modules**: `m_muscl.fpp`

**Current gap**: Some MUSCL limiter combinations not fully exercised.

**Current implementation**: Already exists in `alter_muscl()` but could be expanded.

**Enhancement**:

```python
def alter_muscl():
    # ... existing code ...
    
    # Add test with muscl_order=0 (disable MUSCL)
    stack.push("muscl_order=0", {
        'muscl_order': 0,
        'recon_type': 2,
        'weno_order': 0
    })
    cases.append(define_case_d(stack, '', {}))
    stack.pop()
```

---

## Priority 7: Surface Tension Edge Cases

**Target modules**: `m_surface_tension.fpp`

**Current gap**: Minimal testing of surface tension implementation.

**Add to `toolchain/mfc/test/cases.py`**:

```python
def alter_surface_tension_variants():
    """Test surface tension with different configurations"""
    stack.push("Surface Tension", {
        'surface_tension': 'T',
        'sigma': 1.0,
        'model_eqns': 3,
        'num_fluids': 2,
        'fluid_pp(2)%gamma': 2.5,
        'fluid_pp(2)%pi_inf': 0.0,
        'patch_icpp(1)%cf_val': 1,
        'patch_icpp(2)%cf_val': 0,
        'patch_icpp(3)%cf_val': 1,
        't_step_stop': 10,
        't_step_save': 10
    })
    
    # Test with different surface tension coefficients
    for sigma in [0.1, 1.0, 10.0]:
        stack.push(f"sigma={sigma}", {'sigma': sigma})
        cases.append(define_case_d(stack, '', {}))
        stack.pop()
    
    stack.pop()
```

**Integration point**: Add in `foreach_dimension()` for 2D and 3D cases only.

---

## Implementation Guide

### Step 1: Add Functions to `cases.py`

Add the new test generator functions to `toolchain/mfc/test/cases.py`:

```python
# Add after line 365 (after alter_ib definition):

def alter_time_steppers(dimInfo):
    # ... code from Priority 1 ...

def alter_rare_bcs(dimInfo):
    # ... code from Priority 2 ...

# etc.
```

### Step 2: Integrate into `foreach_dimension()`

Update `foreach_dimension()` function (around line 972):

```python
def foreach_dimension():
    for dimInfo, dimParams in get_dimensions():
        stack.push(f"{len(dimInfo[0])}D", dimParams)
        
        # Existing tests
        alter_bcs(dimInfo)
        alter_grcbc(dimInfo)
        alter_weno(dimInfo)
        alter_muscl()
        alter_num_fluids(dimInfo)
        
        # NEW: Add targeted coverage tests
        alter_time_steppers(dimInfo)      # Priority 1
        alter_rare_bcs(dimInfo)            # Priority 2
        alter_riemann_variants()           # Priority 4
        alter_output_formats()             # Priority 5
        
        if len(dimInfo[0]) == 2:
            alter_2d()
            alter_surface_tension_variants()  # Priority 7 (2D only)
        
        if len(dimInfo[0]) == 3:
            alter_3d()
        
        # ... rest of existing code ...
        
        stack.pop()
```

### Step 3: Generate Golden Files

After adding new tests:

```bash
# List new tests to verify they were added
./mfc.sh test -l | grep -E "(time_stepper|CFL|bc=-1[3489]|weno_avg|sigma)"

# Generate golden files for ONLY the new tests
./mfc.sh test --generate -o "time_stepper" -j 8
./mfc.sh test --generate -o "CFL" -j 8
./mfc.sh test --generate -o "bc=-13" -j 8
./mfc.sh test --generate -o "bc=-14" -j 8
# ... etc for each new test pattern
```

### Step 4: Verify New Tests

```bash
# Run the new tests
./mfc.sh test -o "time_stepper" -j 8
./mfc.sh test -o "CFL Adaptive" -j 8

# Check coverage improvement
PERCENT=100 ./toolchain/coverage.sh
```

---

## Expected Coverage Improvements

### Before (Current)
- `m_time_steppers.fpp`: ~30-40%
- `m_cbc.fpp`: ~20-30%
- `m_viscous.fpp`: ~40-50%
- `m_riemann_solvers.fpp`: ~60-70%
- **Overall**: ~45-55%

### After (With All Priorities)
- `m_time_steppers.fpp`: ~70-80% (+30-40%)
- `m_cbc.fpp`: ~50-60% (+30%)
- `m_viscous.fpp`: ~65-75% (+20%)
- `m_riemann_solvers.fpp`: ~80-90% (+20%)
- **Overall**: ~65-75% (+20%)

---

## Minimal Quick Wins (Priority 1-2 Only)

If time is limited, implement **only Priority 1 and 2**:
- Add ~15-20 new test cases
- Takes ~2-3 hours to implement and generate golden files
- Should increase coverage by ~10-15%
- Focuses on most under-covered critical modules

---

## Notes

1. **Keep grids small**: Use default dimensions or smaller (m=49 for 1D, m=n=24 for 2D, etc.)
2. **Short runs**: `t_step_stop=10-20` is sufficient for code coverage
3. **Tolerance**: New tests may need custom tolerances; adjust in `TestCase.compute_tolerance()`
4. **Test time**: Each test should run in < 5 seconds; total addition ~2-5 minutes to test suite

---

## Validation Checklist

After adding new tests:

- [ ] All new tests generate golden files successfully
- [ ] All new tests pass validation
- [ ] No duplicate UUIDs (checked automatically by `list_cases()`)
- [ ] Coverage report shows improvement in target modules
- [ ] Test suite total runtime increased by < 20%
- [ ] Documentation updated (if new features exposed)

---

## Future Expansions

Additional areas for future test expansion (lower priority):

- MHD solver variants (different magnetic field configurations)
- Bubbles Euler/Lagrange edge cases
- Phase change with different fluids
- Hypo/hyperelastic material parameters
- Acoustic source combinations
- IGR solver variants

These can be added incrementally as coverage targets increase.







