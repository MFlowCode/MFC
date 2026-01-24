# Phase 1 Validation: Multiphase Chemistry Coupling

## Overview

This document describes the validation tests for Phase 1 of the multiphase 
chemistry coupling implementation. Phase 1 enables:

1. Chemistry reactions only in gas-phase regions
2. Automatic transfer of evaporated mass to fuel species
3. Input validation for coupling parameters

## Test Cases

### Test 1: Chemistry Skipping in Liquid Cells

**Objective:** Verify that chemistry reactions are not computed in liquid-dominated cells.

**Setup:**
- 1D domain with liquid droplet (alpha_liquid = 0.99) on left, gas on right
- Chemistry enabled with reactions
- `chem_params%multiphase = T`
- `chem_params%gas_phase_threshold = 0.01`

**Expected Outcomes:**

| Region | alpha_liquid | alpha_gas | Chemistry Active? | Expected Behavior |
|--------|--------------|-----------|-------------------|-------------------|
| Liquid core | 0.99 | 0.01 | NO | No species production rates |
| Interface | 0.50 | 0.50 | YES | Normal chemistry |
| Gas phase | 0.01 | 0.99 | YES | Normal chemistry |

**Verification:**
- In liquid cells: `omega_k = 0` for all species
- In gas cells: `omega_k != 0` when reactants present
- No NaN or Inf values anywhere

---

### Test 2: Evaporation Mass Transfer to Species

**Objective:** Verify that evaporated liquid mass is added to the fuel species.

**Setup:**
- 1D domain with evaporating interface
- Phase change enabled (`relax = T`, `relax_model = 6`)
- Chemistry enabled with `chem_params%multiphase = T`
- `chem_params%fuel_species_idx = 1` (fuel is first species)

**Expected Outcomes:**

| Quantity | Before Evaporation | After Evaporation | Conservation Check |
|----------|-------------------|-------------------|-------------------|
| Liquid mass (m1) | m1_initial | m1_final < m1_initial | - |
| Fuel species (rho*Y_fuel) | Y_fuel_initial | Y_fuel_final > Y_fuel_initial | - |
| Mass evaporated | - | dm = m1_initial - m1_final | dm > 0 |
| Species mass gained | - | dY = Y_fuel_final - Y_fuel_initial | dY ≈ dm |

**Verification:**
- `dm_evap = m1_old - m1_new` is positive when evaporation occurs
- Fuel species mass increases by approximately `dm_evap`
- Total mass is conserved: `sum(alpha_i * rho_i) + sum(rho * Y_k)` constant

---

### Test 3: Input Validation

**Objective:** Verify that invalid configurations are caught by the checker.

**Test Matrix:**

| Test | Configuration | Expected Result |
|------|--------------|-----------------|
| 3a | `multiphase=T`, `relax=F` | ERROR: "requires relax = T" |
| 3b | `multiphase=T`, `num_fluids=1` | ERROR: "requires num_fluids >= 2" |
| 3c | `liquid_phase_idx=0` | ERROR: "must be in range [1, num_fluids]" |
| 3d | `liquid_phase_idx=5`, `num_fluids=3` | ERROR: "must be in range [1, num_fluids]" |
| 3e | `fuel_species_idx=0` | ERROR: "must be in range [1, num_species]" |
| 3f | `gas_phase_threshold=-0.1` | ERROR: "must be in range [0, 1]" |
| 3g | `gas_phase_threshold=1.5` | ERROR: "must be in range [0, 1]" |
| 3h | Valid configuration | SUCCESS: No errors |

---

### Test 4: Boundary Conditions

**Objective:** Verify behavior at domain boundaries with mixed phases.

**Setup:**
- Reflective or periodic boundaries
- Droplet near boundary

**Expected Outcomes:**

| Boundary Type | Expected Behavior |
|---------------|-------------------|
| Reflective | No chemistry in ghost cells if liquid-dominated |
| Periodic | Consistent behavior across periodic boundary |
| Outflow | Species can exit domain normally |

---

### Test 5: Conservation Properties

**Objective:** Verify mass and energy conservation.

**Quantities to Track:**

| Quantity | Formula | Should Be |
|----------|---------|-----------|
| Total mass | `sum(alpha_i * rho_i)` over domain | Constant (closed system) |
| Species mass | `sum(rho * Y_k)` for each species | Conserved per reactions |
| Total energy | `sum(E)` over domain | Constant (adiabatic) |
| Element mass | e.g., H atoms, O atoms | Strictly conserved |

---

### Test 6: Threshold Sensitivity

**Objective:** Verify that `gas_phase_threshold` behaves correctly.

**Test Matrix:**

| Threshold | Cell alpha_gas | Chemistry? | Notes |
|-----------|----------------|------------|-------|
| 0.01 | 0.005 | NO | Below threshold |
| 0.01 | 0.015 | YES | Above threshold |
| 0.01 | 0.01 | NO | At threshold (< not <=) |
| 0.10 | 0.05 | NO | Higher threshold |
| 0.10 | 0.15 | YES | Above higher threshold |
| 0.00 | 0.001 | YES | Zero threshold = always active in any gas |

---

## Validation Metrics

### Quantitative Checks

1. **Mass Conservation Error:**
   ```
   error_mass = |M(t) - M(0)| / M(0)
   ```
   Expected: < 1e-10 (machine precision)

2. **Species Balance:**
   ```
   error_species = |sum(rho*Y_k) - expected| / expected
   ```
   Expected: < 1e-8

3. **Energy Conservation (adiabatic):**
   ```
   error_energy = |E(t) - E(0)| / E(0)
   ```
   Expected: < 1e-8

### Qualitative Checks

1. No NaN or Inf values in any field
2. All volume fractions remain in [0, 1]
3. All mass fractions remain in [0, 1]
4. Pressure and temperature remain positive
5. Chemistry only active where expected

---

## Test Case Parameters

### Minimal 1D Evaporating Interface

```python
# Domain
m = 199
n = 0
p = 0
x_domain_beg = 0.0
x_domain_end = 1.0e-3  # 1 mm

# Fluids (3-fluid model)
num_fluids = 3
# Fluid 1: Liquid fuel
# Fluid 2: Fuel vapor  
# Fluid 3: Oxidizer (air)

# Phase change
relax = True
relax_model = 6

# Chemistry
chemistry = True
cantera_file = "h2o2.yaml"
num_species = 10

# Multiphase coupling (Phase 1)
chem_params_multiphase = True
chem_params_liquid_phase_idx = 1
chem_params_fuel_species_idx = 1  # H2
chem_params_gas_phase_threshold = 0.01

# Initial conditions
# Patch 1: Liquid droplet (left half)
#   alpha = [0.99, 0.005, 0.005]
#   Y = [0, 0, ..., 0]  # No species in liquid
# Patch 2: Gas (right half)  
#   alpha = [0.0, 0.0, 1.0]
#   Y = [0, 0, 0, 0.233, 0, 0, 0, 0, 0, 0.767]  # Air (O2 + N2)

# Time stepping
t_step_stop = 1000
dt = 1e-9
```

---

## Success Criteria

Phase 1 validation is **PASSED** if:

1. All Test 1-6 pass without errors
2. Mass conservation error < 1e-10
3. No numerical instabilities (NaN, Inf)
4. Chemistry correctly skips liquid cells
5. Evaporated mass appears in fuel species
6. Invalid configurations are rejected

Phase 1 validation is **FAILED** if:

1. Any test crashes or produces NaN/Inf
2. Mass conservation violated
3. Chemistry runs in liquid cells
4. Evaporated mass not transferred to species
5. Invalid configurations accepted

---

## Running the Tests

```bash
# Navigate to test directory
cd examples/2D_burning_droplet

# Run validation test
./mfc.sh run ./test_phase1_validation.py -t pre_process simulation -j $(nproc)

# Check results
python3 validate_results.py
```

---

## Known Limitations (Phase 1)

1. **Single fuel species:** Only one species receives evaporated mass
2. **No heat release coupling:** Chemistry heat doesn't affect evaporation rate
3. **Simple threshold:** Binary on/off based on volume fraction
4. **No diffusion modification:** Diffusion flux not adjusted at interfaces

These limitations will be addressed in Phase 2.
