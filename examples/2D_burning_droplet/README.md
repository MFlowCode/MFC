# 2D Burning Droplet Simulation

This example demonstrates how to simulate a burning droplet using MFC's chemistry capabilities. The simulation models a fuel vapor "droplet" surrounded by oxidizer, where species diffusion and chemical reactions lead to combustion.

## Physics Overview

### Droplet Combustion
In a burning droplet, several coupled processes occur:
1. **Vaporization**: Liquid fuel evaporates from the droplet surface
2. **Diffusion**: Fuel vapor diffuses outward while oxidizer diffuses inward
3. **Mixing**: At the flame front, fuel and oxidizer mix at stoichiometric ratios
4. **Reaction**: Chemical reactions release heat, sustaining combustion
5. **Heat Transfer**: Heat from the flame preheats incoming reactants

### Simplified Gas-Phase Model
This example uses a **gas-phase approximation** where:
- The "droplet" is a region of pre-vaporized fuel vapor with smooth transitions
- A tanh profile creates realistic fuel-oxidizer mixing layers
- Chemistry handles species transport and reactions
- No explicit liquid-vapor phase change (simplified model)

This approach captures the essential mixing and combustion physics while avoiding the complexity of multiphase coupling.

### Current Limitation: Phase Change vs Chemistry

MFC has two physics modules that use different approaches:

| Feature | Phase Change | Chemistry |
|---------|-------------|-----------|
| Tracking | Volume fractions (α) | Mass fractions (Y) |
| Fluids | `num_fluids > 1` (liquid, vapor, gas) | `num_fluids = 1` |
| Model | 5/6-equation + relaxation | 5-equation + species |
| Physics | Evaporation/condensation | Reactions/diffusion |

**The modules are not yet coupled.** This means:
- Phase change can vaporize a droplet but won't react the vapor
- Chemistry can burn fuel vapor but needs pre-vaporized initial conditions

### Workarounds for True Burning Droplet

**Option 1: Gas-Phase Approximation** (this example, `case.py`)
- Assume droplet is already vaporized
- Use smooth fuel profile to mimic diffusion from droplet surface
- Chemistry handles mixing and combustion
- Best for studying flame dynamics and combustion chemistry

**Option 2: Two-Stage Simulation**
1. Run `case_liquid_droplet.py` to get vapor distribution over time
2. Extract vapor field at specific time
3. Use as initial condition for chemistry case
4. Repeat to capture quasi-steady behavior

**Option 3: Future MFC Development**
Coupling phase change with chemistry would require:
- Tracking species mass fractions within each fluid phase
- Modifying the equation system to include both α and Y
- This is an area of active research/development

## Running the Simulation

### Prerequisites
- MFC compiled with chemistry support
- Cantera installed with the `h2o2.yaml` mechanism

### Quick Test Run
```bash
# Fast mode for testing (coarser grid, shorter time)
./mfc.sh run examples/2D_burning_droplet/case.py -t pre_process simulation -j $(nproc) -- --fast
```

### Full Simulation
```bash
# Pre-process and run full simulation
./mfc.sh run examples/2D_burning_droplet/case.py -t pre_process simulation -j $(nproc)
```

### Command-Line Options
```bash
# Scale up resolution (2x finer grid)
./mfc.sh run examples/2D_burning_droplet/case.py -- --scale 2

# Adjust initial droplet temperature (controls ignition)
./mfc.sh run examples/2D_burning_droplet/case.py -- --T_droplet 2000

# Disable chemistry (for debugging/pure mixing)
./mfc.sh run examples/2D_burning_droplet/case.py -- --no-chem

# Disable diffusion (for debugging)
./mfc.sh run examples/2D_burning_droplet/case.py -- --no-diffusion

# Fast test mode (coarse grid, short time)
./mfc.sh run examples/2D_burning_droplet/case.py -- --fast
```

## Configuration Details

### Chemistry Mechanism
The simulation uses the `h2o2.yaml` mechanism from Cantera, which includes:
- **Species**: H2, H, O, O2, OH, H2O, HO2, H2O2, AR, N2
- **Reactions**: Detailed hydrogen-oxygen combustion kinetics

For hydrocarbon fuel droplets, consider using:
- `gri30.yaml` - Natural gas (methane) combustion
- Custom mechanisms for heavier fuels (heptane, dodecane)

### Initial Conditions
| Region | Composition | Temperature |
|--------|-------------|-------------|
| Droplet | Pure H2 (fuel) | 1200 K (ignited) |
| Ambient | O2 + AR (oxidizer) | 300 K |

### Key Parameters
| Parameter | Value | Description |
|-----------|-------|-------------|
| `droplet_radius` | 1 mm | Initial droplet size |
| `domain_size` | 10 mm | Computational domain |
| `P_ambient` | 1 atm | Ambient pressure |
| `transport_model` | 2 | Unity-Lewis number |

## Expected Results

### Flame Structure
- **Flame front**: Thin reaction zone at the stoichiometric radius
- **Temperature peak**: ~2500-3000 K in the reaction zone
- **Products**: H2O accumulates in the flame region

### Observable Phenomena
1. Initial mixing layer development
2. Ignition and flame establishment
3. Quasi-steady flame propagation
4. Species distribution evolution

## Extending the Example

### Different Fuels
To simulate different fuels:
1. Change `ctfile` to appropriate mechanism
2. Update species mass fractions in `Y_droplet` and `Y_ambient`
3. Adjust molecular weights and species indices

### 3D Simulations
For true spherical droplet combustion:
```python
# Modify domain to 3D
"p": Nz,
"z_domain%beg": z_beg,
"z_domain%end": z_end,
# Use spherical patch
"patch_icpp(2)%geometry": 8,  # Sphere
```

### Phase Change + Combustion (Future/Advanced)
For liquid-vapor phase change combined with combustion, there are several approaches:

**Approach 1: Sequential Simulation**
1. Run phase change simulation to get vapor distribution over time
2. Use vapor profile as initial condition for chemistry simulation
3. Repeat to capture quasi-steady behavior

**Approach 2: Source Term Modeling (Advanced)**
Add evaporation source terms to the chemistry simulation using the d² law:
```
dm/dt = -π * d * D * ρ * B * ln(1 + B) / (1 + B^0.7)
```
where B is the transfer number.

**Approach 3: Future MFC Development**
Coupling the phase change module (num_fluids > 1) with chemistry (multispecies) 
would require extending MFC to handle both volume fractions AND species mass 
fractions simultaneously. This is an area of active research.

See `2D_phasechange_bubble` for phase change physics examples.

## Theory: Classical Droplet Burning

The Burke-Schumann flame sheet model predicts:
- **Flame radius**: r_f = r_s * sqrt(1 + B)
- **Burning rate**: proportional to ln(1 + B)

Where B is the transfer number combining thermodynamic and kinetic effects.

## Files in This Example

| File | Description |
|------|-------------|
| `case.py` | **Recommended**: Gas-phase H2-O2 combustion with smooth fuel profile |
| `case_liquid_droplet.py` | Phase change droplet (vaporization only, no chemistry) |
| `case_hydrocarbon.py` | Template for CH4/GRI30 mechanism |
| `viz.py` | Visualization script |
| `README.md` | This documentation |

## Which Case Should I Use?

**For combustion physics (flame, reactions):** Use `case.py`
- Pre-vaporized fuel with smooth transition profile
- Full chemistry with reactions and diffusion
- Captures flame dynamics and species evolution

**For vaporization physics (liquid-vapor):** Use `case_liquid_droplet.py`  
- True liquid droplet with phase change
- Three fluids: liquid, vapor, air
- No chemistry reactions (phase change only)

## Available Chemistry Mechanisms

MFC uses Cantera for chemistry. Common built-in mechanisms:

| Mechanism | Fuel | Species | Description |
|-----------|------|---------|-------------|
| `h2o2.yaml` | H2 | 10 | Hydrogen-oxygen combustion |
| `gri30.yaml` | CH4 | 53 | Natural gas (methane) combustion |

For heavier hydrocarbons (heptane, dodecane, kerosene), you need external mechanisms:
- UC San Diego mechanism: https://web.eng.ucsd.edu/mae/groups/combustion/mechanism.html
- Lawrence Livermore mechanisms: https://combustion.llnl.gov/mechanisms

## Combining Vaporizing Droplet with Combustion

If you have a working vaporizing droplet case and want to add combustion:

### Current Limitation
MFC's phase change (num_fluids > 1) and chemistry (num_fluids = 1) use different 
formulations that aren't directly coupled. Here are workarounds:

### Workaround 1: Gas-Phase Approximation (This Example)
Use the chemistry module with a pre-vaporized fuel profile:
- Set up smooth fuel-oxidizer distribution (tanh profile)
- Use high initial temperature for ignition
- Chemistry handles diffusion and reactions

### Workaround 2: Sequential Simulation
1. Run phase change simulation (your vaporizing droplet case)
2. Extract vapor volume fraction at desired time
3. Convert to species mass fraction for chemistry simulation
4. Run chemistry simulation with this initial condition

### Workaround 3: External Coupling
Use a script to:
1. Run phase change for small time step
2. Extract evaporated fuel mass
3. Add to chemistry domain as source term
4. Repeat

## Troubleshooting

### Chemistry doesn't ignite
- Increase `--T_droplet` (try 1500-2500 K)
- Ensure fuel and oxidizer overlap (check transition sharpness)
- Verify species mass fractions sum to 1.0

### Simulation crashes
- Reduce CFL (edit `cfl` in case.py)
- Use `--fast` mode first
- Check boundary conditions for acoustic reflections

### Wrong species
- Verify mechanism file exists and is correct
- Check species indices match mechanism ordering
- Use Cantera Python to verify: `ct.Solution('h2o2.yaml').species_names`

## References

1. Williams, F.A. (1985) "Combustion Theory"
2. Law, C.K. (2006) "Combustion Physics" 
3. Turns, S.R. (2011) "An Introduction to Combustion"
4. MFC Documentation: https://mflowcode.github.io
5. Cantera Documentation: https://cantera.org
