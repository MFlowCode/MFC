#!/usr/bin/env python3
"""
2D Burning Droplet Simulation for MFC

This example simulates a burning "droplet" of fuel vapor in an oxidizer environment.
The simulation models:
1. A spherical (cylindrical in 2D) region of fuel vapor
2. Surrounding ambient oxidizer (air: O2 + N2 or AR)
3. Species diffusion allowing fuel-oxidizer mixing
4. Chemical reactions for combustion

This is a gas-phase combustion model where:
- The "droplet" is pre-vaporized fuel (H2 or other fuel species)
- Chemistry enables reactions when fuel and oxidizer mix at proper stoichiometry
- Diffusion transports species and heat across the domain

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
    "--ignite",
    default=True,
    action="store_true",
    help="Enable ignition (hot droplet to initiate combustion).",
)

args = parser.parse_args()

# =============================================================================
# MECHANISM FILE
# =============================================================================
# Using H2-O2 mechanism (h2o2.yaml) from Cantera's built-in mechanisms
# This is a detailed mechanism for hydrogen-oxygen combustion
# For hydrocarbon fuels, use mechanisms like: gri30.yaml, n-heptane, etc.
ctfile = "h2o2.yaml"

# =============================================================================
# PHYSICAL PARAMETERS
# =============================================================================
# Universal gas constant
Ru = 8.314462  # J/(mol·K)

# Droplet parameters
droplet_radius = 0.001  # 1 mm radius droplet
domain_size = 0.01      # 10 mm domain (10x droplet radius)

# Thermodynamic conditions
# Droplet (fuel) region - hydrogen-rich, hot to initiate combustion
T_droplet = 1200.0 if args.ignite else 300.0  # K (hot for ignition)
P_ambient = 101325.0   # Pa (1 atm)

# Ambient (oxidizer) region - air-like mixture at room temperature
T_ambient = 300.0      # K

# =============================================================================
# SPECIES MASS FRACTIONS
# =============================================================================
# H2-O2 mechanism species (from h2o2.yaml):
# Species order: H2, H, O, O2, OH, H2O, HO2, H2O2, AR, N2
# Indices:        1,  2, 3,  4,  5,   6,   7,    8,  9, 10

# Note: The number of species depends on the mechanism file
# For h2o2.yaml, there are typically 10 species
num_species = 10

# Droplet composition: Pure hydrogen fuel
# Y_H2 = 1.0 (fuel only)
Y_droplet = [0.0] * num_species
Y_droplet[0] = 1.0  # H2 (hydrogen fuel)

# Ambient composition: Stoichiometric H2-Air or H2-O2-AR mixture
# For H2 + 0.5 O2 -> H2O, stoichiometric mass ratio
# Y_O2 ≈ 0.233 (typical air), Y_N2 or AR ≈ 0.767
Y_ambient = [0.0] * num_species
Y_ambient[3] = 0.233   # O2 (oxygen)
Y_ambient[8] = 0.767   # AR (argon) - using AR instead of N2 for simpler kinetics

# =============================================================================
# THERMODYNAMIC STATE CALCULATIONS
# =============================================================================
# For ideal gas: P = rho * R_specific * T
# where R_specific = Ru / M_mixture

# Molecular weights (g/mol -> kg/mol for SI)
MW = {
    'H2': 2.016e-3,
    'H': 1.008e-3,
    'O': 15.999e-3,
    'O2': 31.998e-3,
    'OH': 17.007e-3,
    'H2O': 18.015e-3,
    'HO2': 33.006e-3,
    'H2O2': 34.014e-3,
    'AR': 39.948e-3,
    'N2': 28.014e-3
}

MW_list = [MW['H2'], MW['H'], MW['O'], MW['O2'], MW['OH'], 
           MW['H2O'], MW['HO2'], MW['H2O2'], MW['AR'], MW['N2']]

def compute_mixture_MW(Y):
    """Compute mixture molecular weight from mass fractions."""
    sum_Y_over_M = sum(Y[i] / MW_list[i] for i in range(len(Y)) if Y[i] > 0)
    return 1.0 / sum_Y_over_M if sum_Y_over_M > 0 else MW_list[0]

def compute_density(P, T, Y):
    """Compute density from equation of state."""
    M_mix = compute_mixture_MW(Y)
    R_specific = Ru / M_mix
    return P / (R_specific * T)

# Compute densities
rho_droplet = compute_density(P_ambient, T_droplet, Y_droplet)
rho_ambient = compute_density(P_ambient, T_ambient, Y_ambient)

# =============================================================================
# NUMERICAL PARAMETERS
# =============================================================================
# Grid resolution
Nx_base = 200   # Base resolution in x
Ny_base = 200   # Base resolution in y
Nx = int(Nx_base * args.scale)
Ny = int(Ny_base * args.scale)

# Domain
L = domain_size
x_beg = -L / 2
x_end = L / 2
y_beg = 0.0  # Symmetric about y=0, simulate upper half
y_end = L / 2

# Grid spacing
dx = (x_end - x_beg) / Nx
dy = (y_end - y_beg) / Ny

# Time stepping
# CFL-based time step estimation
# Characteristic velocity ~ sound speed ~ sqrt(gamma * R * T)
gamma_gas = 1.4
R_specific = Ru / compute_mixture_MW(Y_ambient)
c_sound = math.sqrt(gamma_gas * R_specific * T_ambient)

# For chemistry, need smaller time step due to stiff reactions
cfl = 0.2
dt = cfl * min(dx, dy) / c_sound

# Simulation time
# Characteristic diffusion time ~ R^2 / D, where D ~ 1e-4 m^2/s for gases
D_approx = 1e-4  # m^2/s (approximate diffusivity)
t_diffusion = droplet_radius**2 / D_approx

# Run for multiple diffusion times to observe combustion
t_end = 0.001  # 1 ms (adjustable based on physics)

Nt = int(t_end / dt)
save_count = 100
Ns = max(1, Nt // save_count)

print(f"# Grid: {Nx} x {Ny}", file=__import__('sys').stderr)
print(f"# Time step: dt = {dt:.2e} s", file=__import__('sys').stderr)
print(f"# Total steps: {Nt}", file=__import__('sys').stderr)
print(f"# Droplet density: {rho_droplet:.4f} kg/m³", file=__import__('sys').stderr)
print(f"# Ambient density: {rho_ambient:.4f} kg/m³", file=__import__('sys').stderr)

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
    "num_patches": 2,        # Background + droplet
    "mpp_lim": "F",
    "mixture_err": "F",
    "time_stepper": 3,       # 3rd order TVD Runge-Kutta
    "weno_order": 5,         # 5th order WENO
    "weno_eps": 1e-16,
    "weno_avg": "F",
    "mapped_weno": "T",
    "mp_weno": "T",
    "riemann_solver": 2,     # HLLC
    "wave_speeds": 1,
    "avg_state": 2,
    
    # -------------------------------------------------------------------------
    # Boundary Conditions
    # -------------------------------------------------------------------------
    # Non-reflecting subsonic buffers to allow outflow
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
    # Patch 1: Background (Oxidizer/Ambient)
    # -------------------------------------------------------------------------
    "patch_icpp(1)%geometry": 3,  # Rectangle
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": y_end / 2,
    "patch_icpp(1)%length_x": L,
    "patch_icpp(1)%length_y": L / 2,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": P_ambient,
    "patch_icpp(1)%alpha(1)": 1.0,
    "patch_icpp(1)%alpha_rho(1)": rho_ambient,
    
    # -------------------------------------------------------------------------
    # Patch 2: Droplet (Fuel Vapor)
    # -------------------------------------------------------------------------
    "patch_icpp(2)%geometry": 2,  # Circle
    "patch_icpp(2)%x_centroid": 0.0,
    "patch_icpp(2)%y_centroid": 0.0,
    "patch_icpp(2)%radius": droplet_radius,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": P_ambient,
    "patch_icpp(2)%alpha(1)": 1.0,
    "patch_icpp(2)%alpha_rho(1)": rho_droplet,
    
    # -------------------------------------------------------------------------
    # Fluid Properties
    # -------------------------------------------------------------------------
    # For chemistry, gamma and pi_inf are computed from species properties
    # These are approximate values for initialization
    "fluid_pp(1)%gamma": 1.0 / (gamma_gas - 1.0),
    "fluid_pp(1)%pi_inf": 0.0,
}

# =============================================================================
# SPECIES MASS FRACTIONS
# =============================================================================
if args.chemistry:
    for i in range(num_species):
        case[f"chem_wrt_Y({i + 1})"] = "T"  # Write all species
        case[f"patch_icpp(1)%Y({i + 1})"] = Y_ambient[i]
        case[f"patch_icpp(2)%Y({i + 1})"] = Y_droplet[i]

# =============================================================================
# OUTPUT
# =============================================================================
if __name__ == "__main__":
    print(json.dumps(case))
