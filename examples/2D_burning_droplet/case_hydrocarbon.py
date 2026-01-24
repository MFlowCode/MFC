#!/usr/bin/env python3
"""
2D Hydrocarbon Burning Droplet Example for MFC

This example demonstrates a configuration more suitable for hydrocarbon fuel
droplet combustion (e.g., methane, which is the primary component of natural gas).

Uses the GRI30 mechanism (gri30.yaml) which includes:
- Methane (CH4) oxidation chemistry
- CO/CO2 formation pathways
- H2O formation
- NOx chemistry (if N2 is present)

The setup mimics a fuel droplet surrounded by air where:
- Center: Fuel-rich (CH4)
- Ambient: Oxidizer (O2 + N2)
- Smooth transition profile for realistic mixing

Note: For heavier hydrocarbons (heptane, dodecane, etc.), you would need
appropriate Cantera mechanism files which may need to be sourced separately.

Author: MFC Development Team
"""

import json
import argparse
import math
import sys

parser = argparse.ArgumentParser(
    prog="2D_burning_droplet_hydrocarbon",
    description="Hydrocarbon fuel droplet combustion simulation",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT",
                    help="MFC's toolchain's internal state.")
parser.add_argument("--scale", type=float, default=1.0, 
                    help="Grid resolution scale factor.")
parser.add_argument("--T_droplet", type=float, default=1800.0,
                    help="Initial droplet temperature (K) for ignition.")
parser.add_argument("--fast", action="store_true",
                    help="Run faster test simulation.")

args = parser.parse_args()

# =============================================================================
# MECHANISM FILE
# =============================================================================
# GRI30 mechanism for natural gas (methane) combustion
# Contains 53 species and 325 reactions
ctfile = "gri30.yaml"

# For methane-air combustion, key species are:
# CH4, O2, N2, H2O, CO, CO2, H2, OH, H, O, etc.
# Check gri30.yaml for full species list and their order

# =============================================================================
# PHYSICAL PARAMETERS
# =============================================================================
Ru = 8.314462  # J/(mol·K)

# Droplet parameters
droplet_radius = 0.5e-3   # 0.5 mm droplet core
domain_size = 3.0e-3      # 3 mm domain (6 droplet radii)

# Thermodynamic conditions
T_droplet = args.T_droplet  # K (hot for ignition)
T_ambient = 300.0           # K
P_ambient = 101325.0        # Pa (1 atm)

# Transition sharpness
transition_k = 3000.0  # 1/m

# =============================================================================
# SIMPLIFIED SPECIES CONFIGURATION
# =============================================================================
# For GRI30, there are 53 species. Key ones for CH4-air combustion:
# We'll focus on the main reactants and products.

# Note: When using GRI30, you need to set up all 53 species
# Here we provide a template - adjust indices based on your mechanism

# Approximate molecular weights for key species (kg/mol)
MW_CH4 = 16.04e-3
MW_O2 = 32.0e-3
MW_N2 = 28.0e-3
MW_H2O = 18.02e-3
MW_CO2 = 44.01e-3
MW_AR = 39.95e-3

# For simplicity with GRI30, we define the initial compositions
# The species indices depend on the specific mechanism file

# =============================================================================
# NUMERICAL PARAMETERS
# =============================================================================
if args.fast:
    Nx = Ny = 75
    t_end = 5e-6
else:
    Nx = int(150 * args.scale)
    Ny = int(150 * args.scale)
    t_end = 1e-4

L = domain_size
x_beg, x_end = -L, L
y_beg, y_end = 0.0, L

dx = (x_end - x_beg) / Nx
dy = (y_end - y_beg) / Ny

# Time stepping
gamma_gas = 1.3  # Slightly lower for hydrocarbon mixtures
c_sound = math.sqrt(gamma_gas * Ru / MW_N2 * T_ambient)
cfl = 0.1
dt = cfl * min(dx, dy) / c_sound

Nt = int(t_end / dt)
Ns = max(1, Nt // (20 if args.fast else 50))

print(f"# Hydrocarbon Droplet Simulation", file=sys.stderr)
print(f"# Mechanism: {ctfile}", file=sys.stderr)
print(f"# Grid: {Nx} x {Ny}", file=sys.stderr)
print(f"# Time step: dt = {dt:.2e} s, Nt = {Nt}", file=sys.stderr)

# =============================================================================
# SPATIAL PROFILES
# =============================================================================
r_expr = "sqrt(x*x + y*y)"
r0 = droplet_radius
phi_expr = f"0.5*(1.0 + tanh({transition_k}*({r_expr} - {r0})))"

# Density estimates
# Fuel-rich (CH4 + small amount of products): use ideal gas law
rho_fuel = P_ambient / (Ru / MW_CH4 * T_droplet)
rho_air = P_ambient / (Ru / MW_N2 * T_ambient)

rho_expr = f"({rho_fuel}*(1.0 - {phi_expr}) + {rho_air}*{phi_expr})"

# =============================================================================
# CASE DICTIONARY
# =============================================================================
# Note: This is a TEMPLATE. The actual species configuration depends on
# how gri30.yaml orders its species. You may need to adjust species indices.

case = {
    "run_time_info": "T",
    
    # Domain
    "x_domain%beg": x_beg,
    "x_domain%end": x_end,
    "y_domain%beg": y_beg,
    "y_domain%end": y_end,
    "m": Nx,
    "n": Ny,
    "p": 0,
    
    # Time
    "dt": float(dt),
    "t_step_start": 0,
    "t_step_stop": Nt,
    "t_step_save": Ns,
    "t_step_print": max(1, Ns // 5),
    
    # Algorithm
    "model_eqns": 2,
    "num_fluids": 1,
    "num_patches": 1,
    "mpp_lim": "F",
    "mixture_err": "F",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1e-16,
    "weno_avg": "T",
    "mapped_weno": "T",
    "mp_weno": "T",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    
    # Boundary conditions
    "bc_x%beg": -6,
    "bc_x%end": -6,
    "bc_y%beg": -2,
    "bc_y%end": -6,
    
    # Chemistry
    "chemistry": "T",
    "chem_params%diffusion": "T",
    "chem_params%reactions": "T",
    "chem_params%transport_model": 2,
    "cantera_file": ctfile,
    
    # Output
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "parallel_io": "T" if args.mfc.get("mpi", True) else "F",
    "chem_wrt_T": "T",
    
    # Patch
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": y_end / 2,
    "patch_icpp(1)%length_x": 2 * L,
    "patch_icpp(1)%length_y": L,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": P_ambient,
    "patch_icpp(1)%alpha(1)": 1.0,
    "patch_icpp(1)%alpha_rho(1)": rho_expr,
    
    # Fluid
    "fluid_pp(1)%gamma": 1.0 / (gamma_gas - 1.0),
    "fluid_pp(1)%pi_inf": 0.0,
}

# =============================================================================
# SPECIES SETUP (TEMPLATE)
# =============================================================================
# IMPORTANT: For GRI30, you need to configure all 53 species!
# This is a template showing the structure. In practice, you should:
# 1. Import cantera and load the mechanism
# 2. Get species list and indices
# 3. Set up proper mass fractions

# Example structure (species indices must match your mechanism):
# For gri30.yaml, typical species indices are:
# CH4 ~ index 13, O2 ~ index 3, N2 ~ index 47, etc.
# Check your specific mechanism file for correct indices.

print("""
WARNING: This is a template for hydrocarbon combustion.
To run this case, you need to:
1. Verify species indices match gri30.yaml ordering
2. Set mass fractions for all 53 species
3. Adjust initial conditions as needed

For quick testing, use the H2-O2 case (case.py) which uses
the simpler h2o2.yaml mechanism with only 10 species.
""", file=sys.stderr)

# Placeholder species configuration
# Replace with actual gri30 species setup
num_species_gri30 = 53  # GRI30 has 53 species

# For now, print the H2O2 case structure as the primary example
print(json.dumps({
    "NOTE": "This is a template. See case.py for working H2-O2 example.",
    "mechanism": ctfile,
    "num_species": num_species_gri30
}))
