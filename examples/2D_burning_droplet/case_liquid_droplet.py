#!/usr/bin/env python3
"""
2D Liquid Fuel Droplet Combustion - Experimental Case

This is an EXPERIMENTAL case attempting to combine:
1. Multiphase flow with phase change (liquid -> vapor)
2. Chemistry for combustion (vapor + O2 -> products)

IMPORTANT: This combination may not be fully supported in MFC.
The phase change module uses volume fractions (num_fluids > 1)
while chemistry uses species mass fractions (typically num_fluids = 1).

This case sets up:
- Fluid 1: Liquid fuel (high pi_inf for stiff liquid EOS)
- Fluid 2: Fuel vapor 
- Fluid 3: Oxidizer (O2 + diluent)

With chemistry attempting to react the vapor with oxidizer.

If this doesn't work, see the README for alternative approaches:
1. Gas-phase approximation (case.py) 
2. Two-stage simulation (phase change then chemistry)

Author: MFC Development Team
"""

import json
import argparse
import math
import sys

parser = argparse.ArgumentParser(
    prog="2D_liquid_burning_droplet",
    description="Experimental liquid droplet combustion case",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT")
parser.add_argument("--scale", type=float, default=1.0)
parser.add_argument("--test-only", action="store_true", 
                    help="Generate minimal case for validation testing")

args = parser.parse_args()

# =============================================================================
# WARNING
# =============================================================================
print("""
================================================================================
EXPERIMENTAL CASE: Liquid Droplet + Chemistry

This case attempts to combine phase change (relax_model) with chemistry.
This combination may not be fully supported in MFC.

If validation fails, use the gas-phase case (case.py) instead.
================================================================================
""", file=sys.stderr)

# =============================================================================
# PHYSICAL PARAMETERS
# =============================================================================
# Universal constants
Ru = 8.3144598  # J/(mol·K)

# Using water/steam as a model system since phase change params are known
# For fuel droplet, you would use appropriate fuel properties

# Pressure and Temperature
P0 = 101325.0   # Pa
T0 = 373.15     # K (boiling point of water at 1 atm)

# Domain
droplet_radius = 0.5e-3  # 0.5 mm
domain_size = 2.0e-3     # 2 mm

# =============================================================================
# FLUID PROPERTIES (3-fluid model for phase change)
# =============================================================================
# Based on 2D_phasechange_bubble example parameters

# Fluid 1: Liquid (water-like, with stiff EOS)
piwl = 1.0e9      # pi_inf for liquid
cvwl = 1816       # cv
cpwl = 4267       # cp
gamwl = cpwl / cvwl
qvwl = -1167000   # energy reference
qvpwl = 0.0

# Fluid 2: Vapor (water vapor-like)
piwv = 0.0        # pi_inf for vapor
gamwv = 1.4
Rv = Ru / (18.01528e-3)  # gas constant for water vapor
cpwv = Rv * gamwv / (gamwv - 1)
cvwv = cpwv / gamwv
qvwv = 2030000
qvpwv = -23400

# Fluid 3: Air/Oxidizer
pia = 0.0
gama = 1.4
Ra = Ru / (28.966e-3)  # gas constant for air
cpa = Ra * gama / (gama - 1)
cva = cpa / gama
qva = 0.0
qvpa = 0.0

# Compute densities
rho_liquid = (P0 + piwl) / ((gamwl - 1) * cvwl * T0)
rho_vapor = (P0 + piwv) / ((gamwv - 1) * cvwv * T0)
rho_air = (P0 + pia) / ((gama - 1) * cva * T0)

# Volume fractions
eps = 1e-8  # Small value for numerical stability

# Background (ambient air)
alpha_wl_bg = eps
alpha_wv_bg = eps
alpha_a_bg = 1.0 - alpha_wl_bg - alpha_wv_bg

# Droplet (mostly liquid with small vapor at interface)
alpha_wl_drop = 1.0 - 2*eps
alpha_wv_drop = eps
alpha_a_drop = eps

# =============================================================================
# NUMERICAL PARAMETERS
# =============================================================================
if args.test_only:
    Nx = Ny = 50
    t_end = 1e-7
else:
    Nx = int(100 * args.scale)
    Ny = int(100 * args.scale)
    t_end = 1e-5

L = domain_size
x_beg, x_end = -L, L
y_beg, y_end = 0.0, L

dx = (x_end - x_beg) / Nx
dy = (y_end - y_beg) / Ny

# Time step
c_sound = math.sqrt(gamwl * (P0 + piwl) / rho_liquid)
cfl = 0.3
dt = cfl * min(dx, dy) / c_sound

Nt = int(t_end / dt)
Ns = max(1, Nt // 20)

print(f"# Grid: {Nx} x {Ny}", file=sys.stderr)
print(f"# dt = {dt:.2e} s, Nt = {Nt}", file=sys.stderr)
print(f"# Liquid density: {rho_liquid:.2f} kg/m³", file=sys.stderr)
print(f"# Vapor density: {rho_vapor:.4f} kg/m³", file=sys.stderr)
print(f"# Air density: {rho_air:.4f} kg/m³", file=sys.stderr)

# =============================================================================
# CASE DICTIONARY
# =============================================================================
case = {
    # Logistics
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
    
    # Algorithm - 3-fluid model with phase change
    "model_eqns": 3,          # 6-equation model for phase change
    "num_fluids": 3,          # liquid, vapor, air
    "num_patches": 2,
    "mpp_lim": "T",
    "mixture_err": "T",
    "time_stepper": 3,
    "weno_order": 3,
    "weno_eps": 1e-16,
    "mapped_weno": "T",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    
    # Phase change relaxation
    "relax": "T",
    "relax_model": 6,         # pTg-equilibrium
    "palpha_eps": 1e-2,
    "ptgalpha_eps": 1e-2,
    
    # Boundary conditions
    "bc_x%beg": -6,
    "bc_x%end": -6,
    "bc_y%beg": -2,
    "bc_y%end": -6,
    
    # Output
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "parallel_io": "T" if args.mfc.get("mpi", True) else "F",
    
    # Patch 1: Background (ambient oxidizer)
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": y_end / 2,
    "patch_icpp(1)%length_x": 2 * L,
    "patch_icpp(1)%length_y": L,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": P0,
    "patch_icpp(1)%alpha_rho(1)": alpha_wl_bg * rho_liquid,
    "patch_icpp(1)%alpha_rho(2)": alpha_wv_bg * rho_vapor,
    "patch_icpp(1)%alpha_rho(3)": alpha_a_bg * rho_air,
    "patch_icpp(1)%alpha(1)": alpha_wl_bg,
    "patch_icpp(1)%alpha(2)": alpha_wv_bg,
    "patch_icpp(1)%alpha(3)": alpha_a_bg,
    
    # Patch 2: Droplet (liquid fuel)
    "patch_icpp(2)%geometry": 2,  # Circle
    "patch_icpp(2)%x_centroid": 0.0,
    "patch_icpp(2)%y_centroid": 0.0,
    "patch_icpp(2)%radius": droplet_radius,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": P0,
    "patch_icpp(2)%alpha_rho(1)": alpha_wl_drop * rho_liquid,
    "patch_icpp(2)%alpha_rho(2)": alpha_wv_drop * rho_vapor,
    "patch_icpp(2)%alpha_rho(3)": alpha_a_drop * rho_air,
    "patch_icpp(2)%alpha(1)": alpha_wl_drop,
    "patch_icpp(2)%alpha(2)": alpha_wv_drop,
    "patch_icpp(2)%alpha(3)": alpha_a_drop,
    
    # Fluid properties
    # Fluid 1: Liquid
    "fluid_pp(1)%gamma": 1.0 / (gamwl - 1),
    "fluid_pp(1)%pi_inf": gamwl * piwl / (gamwl - 1),
    "fluid_pp(1)%cv": cvwl,
    "fluid_pp(1)%qv": qvwl,
    "fluid_pp(1)%qvp": qvpwl,
    
    # Fluid 2: Vapor
    "fluid_pp(2)%gamma": 1.0 / (gamwv - 1),
    "fluid_pp(2)%pi_inf": gamwv * piwv / (gamwv - 1),
    "fluid_pp(2)%cv": cvwv,
    "fluid_pp(2)%qv": qvwv,
    "fluid_pp(2)%qvp": qvpwv,
    
    # Fluid 3: Air/Oxidizer
    "fluid_pp(3)%gamma": 1.0 / (gama - 1),
    "fluid_pp(3)%pi_inf": gama * pia / (gama - 1),
    "fluid_pp(3)%cv": cva,
    "fluid_pp(3)%qv": qva,
    "fluid_pp(3)%qvp": qvpa,
}

# =============================================================================
# CHEMISTRY CONFIGURATION (EXPERIMENTAL)
# =============================================================================
# Uncomment below to attempt adding chemistry to the phase change case
# This may not work and could cause errors

# ctfile = "h2o2.yaml"
# num_species = 10
# 
# case["chemistry"] = "T"
# case["chem_params%diffusion"] = "T"
# case["chem_params%reactions"] = "T"
# case["chem_params%transport_model"] = 2
# case["cantera_file"] = ctfile
# case["chem_wrt_T"] = "T"
# 
# # Species mass fractions would need to be defined here
# # This is where the coupling becomes complex

# =============================================================================
# OUTPUT
# =============================================================================
if __name__ == "__main__":
    print(json.dumps(case))
