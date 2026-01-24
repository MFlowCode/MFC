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
- The "droplet" is a region of pre-vaporized fuel vapor
- Chemistry handles species transport and reactions
- No explicit liquid-vapor phase change (simplified model)

This approach captures the essential mixing and combustion physics while avoiding the complexity of multiphase coupling.

## Running the Simulation

### Prerequisites
- MFC compiled with chemistry support
- Cantera installed with the `h2o2.yaml` mechanism

### Basic Run
```bash
# Pre-process and run simulation
./mfc.sh run examples/2D_burning_droplet/case.py -t pre_process simulation -j $(nproc)
```

### Command-Line Options
```bash
# Scale up resolution (2x)
./mfc.sh run examples/2D_burning_droplet/case.py -t pre_process simulation -- --scale 2

# Disable chemistry (for debugging)
./mfc.sh run examples/2D_burning_droplet/case.py -t pre_process simulation -- --no-chem

# Disable diffusion (for debugging)
./mfc.sh run examples/2D_burning_droplet/case.py -t pre_process simulation -- --no-diffusion
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

### Phase Change (Advanced)
For liquid-vapor phase change combined with combustion:
- Requires coupling MFC's phase change module with chemistry
- Currently experimental - see `2D_phasechange_bubble` for phase change physics

## Theory: Classical Droplet Burning

The Burke-Schumann flame sheet model predicts:
- **Flame radius**: r_f = r_s * sqrt(1 + B)
- **Burning rate**: proportional to ln(1 + B)

Where B is the transfer number combining thermodynamic and kinetic effects.

## References

1. Williams, F.A. (1985) "Combustion Theory"
2. Law, C.K. (2006) "Combustion Physics" 
3. Turns, S.R. (2011) "An Introduction to Combustion"
4. MFC Documentation: https://mflowcode.github.io
