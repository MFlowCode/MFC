#!/usr/bin/env python3
"""
2D Burning Droplet Simulation for MFC

This example simulates a burning "droplet" of fuel vapor in an oxidizer environment.
The simulation models:
1. A spherical (cylindrical in 2D) region of fuel vapor with smooth transition
2. Surrounding ambient oxidizer (air: O2 + AR)
3. Species diffusion allowing fuel-oxidizer mixing
4. Chemical reactions for combustion

This is a gas-phase combustion model where:
- The "droplet" represents vaporized fuel (e.g., from a d² evaporation law)
- A smooth fuel concentration profile mimics diffusion from evaporation
- Chemistry handles species transport and reactions at the flame front

For true liquid-vapor phase change, see the phase change examples.
This example focuses on the gas-phase combustion physics.

References:
+ Williams, F.A. "Combustion Theory" - Classical droplet combustion theory
+ Law, C.K. "Combustion Physics" - Fundamentals of droplet burning
+ MFC Chemistry Documentation

Author: MFC Development Team
"""

import json
import argparse
import math

parser = argparse.ArgumentParser(
    prog="2D_burning_droplet",
    description="2D simulation of a burning fuel droplet in oxidizer environment",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

parser.add_argument(
    "--mfc",
    type=json.loads,
    default="{}",
    metavar="DICT",
    help="MFC's toolchain's internal state.",
)
parser.add_argument(
    "--no-chem",
    dest="chemistry",
    default=True,
    action="store_false",
    help="Disable chemistry (for debugging without reactions).",
)
parser.add_argument(
    "--no-diffusion",
    dest="diffusion",
    default=True,
    action="store_false",
    help="Disable diffusion (for debugging without mixing).",
)
parser.add_argument(
    "--scale",
    type=float,
    default=1.0,
    help="Grid resolution scale factor.",
)
parser.add_argument(
    "--T_droplet",
    type=float,
    default=1500.0,
    help="Initial droplet temperature (K). Set high (>1000K) for ignition.",
)
parser.add_argument(
    "--fast",
    action="store_true",
    help="Run a faster, coarser simulation for testing.",
)

args = parser.parse_args()

# =============================================================================
# MECHANISM FILE
# =============================================================================
# Using H2-O2 mechanism (h2o2.yaml) from Cantera's built-in mechanisms
# For hydrocarbon fuels, use gri30.yaml (methane) or other mechanisms
ctfile = "h2o2.yaml"

# =============================================================================
# PHYSICAL PARAMETERS
# =============================================================================
# Universal gas constant
Ru = 8.314462  # J/(mol·K)

# Droplet parameters
droplet_radius = 0.5e-3    # 0.5 mm radius (core)
flame_radius = 1.5e-3      # 1.5 mm flame standoff (approximate)
domain_size = 5.0e-3       # 5 mm domain

# Thermodynamic conditions
T_droplet = args.T_droplet  # K (hot for ignition, ~1500K)
T_ambient = 300.0           # K (room temperature)
P_ambient = 101325.0        # Pa (1 atm)

# Smoothing parameter for fuel-oxidizer transition
# Larger value = sharper transition, smaller = more diffuse
transition_sharpness = 2000.0  # 1/m (controls mixing layer thickness)

# =============================================================================
# SPECIES CONFIGURATION
# =============================================================================
# H2-O2 mechanism species (from h2o2.yaml):
# Species order: H2, H, O, O2, OH, H2O, HO2, H2O2, AR, N2
# Indices:        1,  2, 3,  4,  5,   6,   7,    8,  9, 10
num_species = 10

# Molecular weights (kg/mol)
MW_list = [2.016e-3, 1.008e-3, 15.999e-3, 31.998e-3, 17.007e-3,
           18.015e-3, 33.006e-3, 34.014e-3, 39.948e-3, 28.014e-3]

# =============================================================================
# INITIAL CONDITIONS WITH SMOOTH TRANSITION
# =============================================================================
# Use analytical expressions for smooth fuel-oxidizer profile
# This mimics a diffusion flame with pre-mixing from evaporation

# Core (droplet center): Rich in fuel (H2), diluted with some products
# The smooth transition uses a tanh profile based on radius

# For the initial condition, we use a hardcoded patch (hcid) approach
# to define spatially varying species mass fractions

# Alternatively, use analytical expressions in patch definitions
# MFC supports expressions like "0.5*(1 - tanh(2000*(sqrt(x^2+y^2) - 0.0005)))"

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================
def compute_mixture_MW(Y):
    """Compute mixture molecular weight from mass fractions."""
    sum_Y_over_M = sum(Y[i] / MW_list[i] for i in range(len(Y)) if Y[i] > 0 and MW_list[i] > 0)
    return 1.0 / sum_Y_over_M if sum_Y_over_M > 0 else MW_list[0]

def compute_density(P, T, Y):
    """Compute density from ideal gas equation of state."""
    M_mix = compute_mixture_MW(Y)
    R_specific = Ru / M_mix
    return P / (R_specific * T)

# =============================================================================
# REFERENCE COMPOSITIONS
# =============================================================================
# Droplet (fuel-rich) composition at center
Y_fuel_rich = [0.0] * num_species
Y_fuel_rich[0] = 0.8   # H2 (hydrogen fuel) - not pure to avoid numerical issues
Y_fuel_rich[8] = 0.2   # AR (argon as diluent)

# Ambient (oxidizer) composition at infinity
Y_oxidizer = [0.0] * num_species
Y_oxidizer[3] = 0.233   # O2 (oxygen)
Y_oxidizer[8] = 0.767   # AR (argon)

# Compute reference densities
rho_fuel = compute_density(P_ambient, T_droplet, Y_fuel_rich)
rho_oxidizer = compute_density(P_ambient, T_ambient, Y_oxidizer)

# =============================================================================
# NUMERICAL PARAMETERS
# =============================================================================
# Grid resolution
if args.fast:
    Nx_base = 100
    Ny_base = 100
    t_end = 1e-5  # Very short for testing
else:
    Nx_base = 200
    Ny_base = 200
    t_end = 2e-4  # 0.2 ms for flame development

Nx = int(Nx_base * args.scale)
Ny = int(Ny_base * args.scale)

# Domain
L = domain_size
x_beg = -L
x_end = L
y_beg = 0.0  # Symmetric about y=0, simulate upper half
y_end = L

# Grid spacing
dx = (x_end - x_beg) / Nx
dy = (y_end - y_beg) / Ny

# Time stepping
gamma_gas = 1.4
R_specific = Ru / compute_mixture_MW(Y_oxidizer)
c_sound = math.sqrt(gamma_gas * R_specific * max(T_ambient, T_droplet))

# CFL-limited time step
cfl = 0.15
dt = cfl * min(dx, dy) / c_sound

Nt = int(t_end / dt)
save_count = 50 if args.fast else 100
Ns = max(1, Nt // save_count)

# Diagnostic output
import sys
print(f"# Grid: {Nx} x {Ny}", file=sys.stderr)
print(f"# Domain: [{x_beg*1000:.1f}, {x_end*1000:.1f}] x [{y_beg*1000:.1f}, {y_end*1000:.1f}] mm", file=sys.stderr)
print(f"# Time step: dt = {dt:.2e} s", file=sys.stderr)
print(f"# Total steps: {Nt}, saves every {Ns} steps", file=sys.stderr)
print(f"# Simulation time: {t_end*1e6:.1f} μs", file=sys.stderr)
print(f"# Fuel density: {rho_fuel:.4f} kg/m³ at T={T_droplet:.0f} K", file=sys.stderr)
print(f"# Oxidizer density: {rho_oxidizer:.4f} kg/m³ at T={T_ambient:.0f} K", file=sys.stderr)

# =============================================================================
# DEFINE SPATIAL PROFILES USING ANALYTICAL EXPRESSIONS
# =============================================================================
# MFC supports analytical expressions in patch definitions
# We use a smooth tanh profile to transition from fuel to oxidizer

# Define radius expression (distance from origin in 2D)
r_expr = "sqrt(x*x + y*y)"

# Smooth step function: 0 at center, 1 at infinity
# phi = 0.5 * (1 + tanh(k*(r - r0))) where k controls sharpness
r0 = droplet_radius
k = transition_sharpness
phi_expr = f"0.5*(1.0 + tanh({k}*({r_expr} - {r0})))"

# Species mass fractions (interpolated between fuel and oxidizer)
# Y(r) = Y_fuel * (1 - phi) + Y_oxidizer * phi
def species_expr(Y_fuel_val, Y_ox_val):
    """Generate analytical expression for species mass fraction."""
    if abs(Y_fuel_val - Y_ox_val) < 1e-10:
        return Y_fuel_val  # Constant
    return f"({Y_fuel_val}*(1.0 - {phi_expr}) + {Y_ox_val}*{phi_expr})"

# Temperature profile (smooth transition)
T_expr = f"({T_droplet}*(1.0 - {phi_expr}) + {T_ambient}*{phi_expr})"

# Density is computed from pressure and temperature
# rho = P / (R_mix * T), but R_mix depends on composition
# For simplicity, use a weighted average
rho_expr = f"({rho_fuel}*(1.0 - {phi_expr}) + {rho_oxidizer}*{phi_expr})"

# =============================================================================
# CASE DICTIONARY
# =============================================================================
case = {
    # -------------------------------------------------------------------------
    # Logistics
    # -------------------------------------------------------------------------
    "run_time_info": "T",
    
    # -------------------------------------------------------------------------
    # Computational Domain
    # -------------------------------------------------------------------------
    "x_domain%beg": x_beg,
    "x_domain%end": x_end,
    "y_domain%beg": y_beg,
    "y_domain%end": y_end,
    "m": Nx,
    "n": Ny,
    "p": 0,  # 2D simulation
    
    # -------------------------------------------------------------------------
    # Time Stepping
    # -------------------------------------------------------------------------
    "dt": float(dt),
    "t_step_start": 0,
    "t_step_stop": Nt,
    "t_step_save": Ns,
    "t_step_print": max(1, Ns // 10),
    
    # -------------------------------------------------------------------------
    # Simulation Algorithm
    # -------------------------------------------------------------------------
    "model_eqns": 2,         # 5-equation model
    "num_fluids": 1,         # Single fluid with multiple species (chemistry)
    "num_patches": 1,        # Single patch with analytical profiles
    "mpp_lim": "F",
    "mixture_err": "F",
    "time_stepper": 3,       # 3rd order TVD Runge-Kutta
    "weno_order": 5,         # 5th order WENO
    "weno_eps": 1e-16,
    "weno_avg": "T",         # Average for stability
    "mapped_weno": "T",
    "mp_weno": "T",
    "riemann_solver": 2,     # HLLC
    "wave_speeds": 1,
    "avg_state": 2,
    
    # -------------------------------------------------------------------------
    # Boundary Conditions
    # -------------------------------------------------------------------------
    "bc_x%beg": -6,  # Non-reflecting subsonic buffer
    "bc_x%end": -6,  # Non-reflecting subsonic buffer
    "bc_y%beg": -2,  # Reflective (symmetry plane)
    "bc_y%end": -6,  # Non-reflecting subsonic buffer
    
    # -------------------------------------------------------------------------
    # Chemistry
    # -------------------------------------------------------------------------
    "chemistry": "T" if args.chemistry else "F",
    "chem_params%diffusion": "T" if args.diffusion else "F",
    "chem_params%reactions": "T" if args.chemistry else "F",
    "chem_params%transport_model": 2,  # Unity-Lewis (simpler, more stable)
    "cantera_file": ctfile,
    
    # -------------------------------------------------------------------------
    # Output
    # -------------------------------------------------------------------------
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "parallel_io": "T" if args.mfc.get("mpi", True) else "F",
    "chem_wrt_T": "T",  # Write temperature field
    
    # -------------------------------------------------------------------------
    # Patch 1: Entire Domain with Smooth Profiles
    # -------------------------------------------------------------------------
    "patch_icpp(1)%geometry": 3,  # Rectangle
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": y_end / 2,
    "patch_icpp(1)%length_x": 2 * L,
    "patch_icpp(1)%length_y": L,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": P_ambient,
    "patch_icpp(1)%alpha(1)": 1.0,
    "patch_icpp(1)%alpha_rho(1)": rho_expr,
    
    # -------------------------------------------------------------------------
    # Fluid Properties
    # -------------------------------------------------------------------------
    "fluid_pp(1)%gamma": 1.0 / (gamma_gas - 1.0),
    "fluid_pp(1)%pi_inf": 0.0,
}

# =============================================================================
# SPECIES MASS FRACTIONS (with analytical expressions)
# =============================================================================
if args.chemistry:
    for i in range(num_species):
        case[f"chem_wrt_Y({i + 1})"] = "T"  # Write all species
        Y_expr = species_expr(Y_fuel_rich[i], Y_oxidizer[i])
        case[f"patch_icpp(1)%Y({i + 1})"] = Y_expr

# =============================================================================
# OUTPUT
# =============================================================================
if __name__ == "__main__":
    print(json.dumps(case))
