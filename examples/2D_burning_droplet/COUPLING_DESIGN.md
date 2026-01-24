# Design Document: Coupling Phase Change with Chemistry in MFC

## Overview

This document describes the code modifications required to enable simulation of a **burning liquid droplet** where:
1. Liquid fuel evaporates (phase change)
2. Fuel vapor mixes with oxidizer (diffusion)
3. Mixture combusts (chemistry)

## Current Architecture

### Phase Change Module (`m_phase_change.fpp`)
- Uses `num_fluids > 1` (typically 3: liquid, vapor, gas)
- Tracks **volume fractions** (α_i) for each fluid
- Uses pT or pTg equilibrium relaxation
- Called after RK time stepping: `s_infinite_relaxation_k()`

### Chemistry Module (`m_chemistry.fpp`)
- Uses `num_fluids = 1` with species mass fractions
- Tracks **species mass fractions** (Y_k) globally
- Adds species equations from `chemxb` to `chemxe`
- Computes reaction source terms via Cantera

### Why They Don't Mix
| Aspect | Phase Change | Chemistry |
|--------|-------------|-----------|
| Fluid tracking | Volume fractions α_i | Single fluid |
| Species tracking | None | Mass fractions Y_k |
| Density | α_i * ρ_i per fluid | ρ * Y_k per species |
| Model equations | 5-eqn or 6-eqn | 5-eqn with species |

## Proposed Coupling Approach

### Conceptual Model

For a burning droplet with 3 fluids (liquid fuel, fuel vapor, oxidizer/air):

```
Fluid 1: Liquid Fuel    - Pure fuel, no species tracking needed
Fluid 2: Fuel Vapor     - Contains fuel species (e.g., H2, CH4)  
Fluid 3: Oxidizer/Air   - Contains O2, N2/AR, products (H2O, CO2)
```

The key insight: **species exist within the gas phase only** (fluids 2 and 3).

### Equation System Extension

Current 6-equation model variables:
```
q_cons = [α₁ρ₁, α₂ρ₂, α₃ρ₃,    # Partial densities (3 fluids)
          ρu, ρv, ρw,           # Momentum (3 components)
          E,                     # Total energy
          α₁, α₂, α₃,           # Volume fractions (3 fluids)
          α₁ρ₁e₁, α₂ρ₂e₂, α₃ρ₃e₃]  # Internal energies (6-eqn only)
```

Proposed extension for chemistry:
```
q_cons_extended = [... existing variables ...,
                   α_g ρ_g Y₁,  # Gas-phase species 1
                   α_g ρ_g Y₂,  # Gas-phase species 2
                   ...
                   α_g ρ_g Yₙ]  # Gas-phase species n
```

Where `α_g = α₂ + α₃` (vapor + air volume fractions).

## Required Code Modifications

### 1. Global Parameters (`m_global_parameters.fpp`)

```fortran
! Add new flags
logical :: multiphase_chemistry  ! Enable coupled phase change + chemistry

! Extend species indices to track gas-phase species
! chemxb:chemxe now represent gas-phase species: α_g * ρ_g * Y_k
```

**Lines to modify:** ~1158-1162 (species index setup)

### 2. Pre-processor (`m_initial_condition.fpp`, `m_assign_variables.fpp`)

Add initialization of species within gas phases:
```fortran
! For each cell
! Species are initialized based on which fluid dominates
if (alpha_liquid > 0.5) then
    ! Liquid cell: no species (or trace amounts)
    Y(:) = 0
    Y(fuel_species) = eps  ! Trace fuel
else
    ! Gas cell: full species composition
    Y(:) = patch_icpp%Y(:)
endif
```

**Files to modify:**
- `src/pre_process/m_initial_condition.fpp`
- `src/pre_process/m_assign_variables.fpp`

### 3. Chemistry Module (`m_chemistry.fpp`)

#### 3.1 Modify diffusion flux computation

```fortran
subroutine s_compute_chemistry_diffusion_flux(...)
    ! Add check for gas-phase region
    alpha_gas = q_prim_qp(advxb+1)%sf(x,y,z) + q_prim_qp(advxb+2)%sf(x,y,z)
    
    if (alpha_gas > gas_threshold) then
        ! Compute diffusion only in gas phase
        ! Scale fluxes by alpha_gas
        Mass_Diffu_Flux = alpha_gas * rho_gas * D * dY/dx
    endif
end subroutine
```

**Lines to modify:** ~166-407

#### 3.2 Modify reaction flux computation

```fortran
subroutine s_compute_chemistry_reaction_flux(...)
    ! Only compute reactions in gas phase
    alpha_gas = sum(q_cons_qp(advxb+1:advxe)%sf(x,y,z)) - alpha_liquid
    
    if (alpha_gas > gas_threshold) then
        ! Get gas-phase temperature and species
        call get_gas_phase_properties(...)
        
        ! Compute reaction rates
        call get_net_production_rates(rho_gas, T_gas, Ys_gas, omega)
        
        ! Scale by gas volume fraction
        rhs_vf(eqn)%sf(x,y,z) = rhs_vf(eqn)%sf(x,y,z) + alpha_gas * omega_m
    endif
end subroutine
```

**Lines to modify:** ~118-163

### 4. Phase Change Module (`m_phase_change.fpp`)

#### 4.1 Add species source terms for evaporation

```fortran
subroutine s_infinite_relaxation_k(q_cons_vf)
    ! Existing phase change logic...
    
    ! NEW: When liquid evaporates, add fuel species to gas phase
    if (multiphase_chemistry) then
        dm_evap = m1_new - m1_old  ! Mass transferred from liquid to vapor
        
        ! Add evaporated mass to fuel vapor species
        q_cons_vf(chemxb + fuel_idx - 1)%sf(j,k,l) = &
            q_cons_vf(chemxb + fuel_idx - 1)%sf(j,k,l) + dm_evap
    endif
end subroutine
```

**Lines to modify:** ~83-276 (in `s_infinite_relaxation_k`)

### 5. Variables Conversion (`m_variables_conversion.fpp`)

#### 5.1 Add gas-phase species conversion

```fortran
subroutine s_convert_gas_species_to_mass_fractions(q_vf, i, j, k, Ys_gas)
    ! Compute gas phase volume fraction
    alpha_gas = 0._wp
    rho_gas = 0._wp
    do fl = 2, num_fluids  ! Skip liquid (fluid 1)
        alpha_gas = alpha_gas + q_vf(fl + advxb - 1)%sf(i,j,k)
        rho_gas = rho_gas + q_vf(fl + contxb - 1)%sf(i,j,k)
    end do
    
    ! Convert species conservative to mass fractions
    do sp = chemxb, chemxe
        Ys_gas(sp - chemxb + 1) = q_vf(sp)%sf(i,j,k) / max(rho_gas, eps)
    end do
end subroutine
```

**New subroutine to add**

### 6. Riemann Solvers (`m_riemann_solvers.fpp`)

Modify species flux computation to account for gas-phase only:

```fortran
if (chemistry .and. multiphase_chemistry) then
    ! Species flux is weighted by gas volume fraction
    alpha_gas_L = sum(qL_prim(advxb+1:advxe)) - qL_prim(advxb)
    alpha_gas_R = sum(qR_prim(advxb+1:advxe)) - qR_prim(advxb)
    
    do i = chemxb, chemxe
        flux_rs(i) = 0.5*(alpha_gas_L + alpha_gas_R) * species_flux
    end do
endif
```

**Multiple locations in HLL/HLLC solvers**

### 7. Time Stepper (`m_time_steppers.fpp`)

Ensure proper ordering:
1. RK substeps with chemistry RHS
2. Phase change relaxation (moves mass between phases)
3. Species update for evaporated mass

```fortran
! In s_perform_time_step
call s_tvd_rk(t_step, time_avg, time_stepper)  ! Includes chemistry RHS

if (relax) then
    call s_infinite_relaxation_k(q_cons_ts(1)%vf)  ! Phase change
    
    if (multiphase_chemistry) then
        call s_update_species_from_evaporation(q_cons_ts(1)%vf)  ! NEW
    endif
endif
```

### 8. Input/Output

#### 8.1 Case Validator (`toolchain/mfc/case_validator.py`)

```python
def check_multiphase_chemistry(self):
    if self.get('chemistry') == 'T' and self.get('relax') == 'T':
        self.require(self.get('multiphase_chemistry') == 'T',
                    "Combining phase change with chemistry requires multiphase_chemistry = T")
        self.require(num_fluids >= 2,
                    "Multiphase chemistry requires at least 2 fluids")
```

#### 8.2 Case File Parameters

New parameters needed:
```python
"multiphase_chemistry": "T",           # Enable coupled mode
"liquid_fluid_idx": 1,                 # Which fluid is liquid
"fuel_species_idx": 1,                 # Which species is fuel
"gas_phase_threshold": 0.01,           # Min gas fraction for chemistry
```

## Implementation Priority

### Phase 1: Minimal Viable Product
1. Add `multiphase_chemistry` flag
2. Modify chemistry to skip liquid-dominated cells
3. Add evaporation source term to species equations
4. Test with simple 2-fluid (liquid + gas) case

### Phase 2: Full Coupling
1. Proper gas-phase averaging for thermodynamic properties
2. Multi-species diffusion with phase boundaries
3. Heat release coupling back to phase change
4. Validation against experimental data

### Phase 3: Optimization
1. GPU acceleration for coupled terms
2. Adaptive mesh refinement near flame
3. Subgrid evaporation models

## Estimated Effort

| Component | Files | Complexity | Estimated Lines |
|-----------|-------|------------|-----------------|
| Global params | 1 | Low | ~50 |
| Pre-processor | 2 | Medium | ~100 |
| Chemistry | 1 | High | ~200 |
| Phase change | 1 | High | ~150 |
| Variables conv. | 1 | Medium | ~100 |
| Riemann solvers | 1 | High | ~200 |
| Time stepper | 1 | Low | ~30 |
| Validation | 1 | Low | ~50 |
| **Total** | **9** | | **~880** |

## Testing Strategy

1. **Unit tests**: Each modified subroutine
2. **Regression tests**: Existing phase change and chemistry cases still work
3. **Validation case**: 1D evaporating interface with reaction
4. **Full case**: 2D burning droplet comparison with d² law

## References

1. Turns, S.R. "An Introduction to Combustion" - Droplet burning theory
2. Law, C.K. "Combustion Physics" - Evaporation models
3. MFC Documentation - Existing module interfaces
4. Cantera - Chemistry integration
