#!/usr/bin/env python3
"""
Simple gas-only chemistry test to verify chemistry module works.

This test uses a single fluid (num_fluids=1) with H2-O2 chemistry.
No phase change, no multiphase coupling - just pure gas-phase combustion.

Setup:
- 1D domain with H2 on left, O2+N2 (air) on right
- Temperature high enough for ignition (~1200 K)
- Should see reaction at the interface
"""

import json
import argparse

parser = argparse.ArgumentParser(
    prog="test_chemistry_only",
    description="Simple gas-phase chemistry test",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

parser.add_argument("--T", type=float, default=1200.0,
                    help="Initial temperature (K)")
parser.add_argument("--no-reactions", action="store_true",
                    help="Disable reactions (diffusion only)")

args, _ = parser.parse_known_args()

# =============================================================================
# DOMAIN PARAMETERS
# =============================================================================
Nx = 199
Lx = 1.0e-2  # 1 cm domain (larger for diffusion)

# =============================================================================
# THERMODYNAMIC STATE
# =============================================================================
T0 = args.T  # Temperature (K)
p0 = 1.01325e5  # Pressure (Pa) - 1 atm

# Gas constant
R_universal = 8314.46  # J/(kmol·K)

# Approximate properties for H2-air mixture
gamma = 1.4
W_mix = 20.0  # Approximate molecular weight (kg/kmol)
R_gas = R_universal / W_mix  # J/(kg·K)

# Density from ideal gas law
rho0 = p0 / (R_gas * T0)

# =============================================================================
# CHEMISTRY CONFIGURATION
# =============================================================================
ctfile = "h2o2.yaml"
num_species = 10

# Species indices in h2o2.yaml:
# 1: H2, 2: H, 3: O, 4: O2, 5: OH, 6: H2O, 7: HO2, 8: H2O2, 9: AR, 10: N2
idx_H2 = 1
idx_O2 = 4
idx_N2 = 10

# Mass fractions
# Left side: pure H2
Y_H2_left = 1.0

# Right side: air (O2 + N2)
Y_O2_right = 0.233
Y_N2_right = 0.767

# =============================================================================
# TIME STEPPING
# =============================================================================
dt = 1.0e-8  # 10 ns time step
t_stop = 100  # 100 steps
t_save = 10   # Save every 10 steps

# =============================================================================
# CASE DICTIONARY
# =============================================================================
case = {
    # -------------------------------------------------------------------------
    # Logistics
    # -------------------------------------------------------------------------
    "run_time_info": "T",
    
    # -------------------------------------------------------------------------
    # Domain
    # -------------------------------------------------------------------------
    "m": Nx,
    "n": 0,
    "p": 0,
    
    "x_domain%beg": 0.0,
    "x_domain%end": Lx,
    
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": t_stop,
    "t_step_save": t_save,
    
    # -------------------------------------------------------------------------
    # Model - Single fluid for chemistry
    # -------------------------------------------------------------------------
    "model_eqns": 2,  # 5-equation model (standard for chemistry)
    "num_fluids": 1,
    "num_patches": 2,
    "mpp_lim": "F",  # Not supported with num_fluids=1
    "mixture_err": "T",
    "time_stepper": 3,  # 3rd order TVD RK
    
    # -------------------------------------------------------------------------
    # Numerics
    # -------------------------------------------------------------------------
    "weno_order": 5,
    "weno_eps": 1.0e-16,
    "mapped_weno": "T",
    "mp_weno": "T",
    "riemann_solver": 2,  # HLLC
    "wave_speeds": 1,
    "avg_state": 2,
    
    # -------------------------------------------------------------------------
    # Boundary Conditions
    # -------------------------------------------------------------------------
    "bc_x%beg": -3,  # Reflective
    "bc_x%end": -3,  # Reflective
    
    # -------------------------------------------------------------------------
    # Patch 1: Background - Air (O2 + N2)
    # -------------------------------------------------------------------------
    "patch_icpp(1)%geometry": 1,  # Line
    "patch_icpp(1)%x_centroid": Lx / 2,
    "patch_icpp(1)%length_x": Lx,
    
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%pres": p0,
    
    # Single fluid - volume fraction = 1
    "patch_icpp(1)%alpha(1)": 1.0,
    "patch_icpp(1)%alpha_rho(1)": rho0,
    
    # -------------------------------------------------------------------------
    # Patch 2: Fuel region (left 30%)
    # -------------------------------------------------------------------------
    "patch_icpp(2)%geometry": 1,  # Line
    "patch_icpp(2)%x_centroid": 0.15 * Lx,
    "patch_icpp(2)%length_x": 0.3 * Lx,
    "patch_icpp(2)%alter_patch(1)": "T",
    
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%pres": p0,
    
    "patch_icpp(2)%alpha(1)": 1.0,
    "patch_icpp(2)%alpha_rho(1)": rho0,
    
    # -------------------------------------------------------------------------
    # Fluid Properties (ideal gas)
    # -------------------------------------------------------------------------
    "fluid_pp(1)%gamma": 1.0 / (gamma - 1),
    "fluid_pp(1)%pi_inf": 0.0,
    
    # -------------------------------------------------------------------------
    # Chemistry
    # -------------------------------------------------------------------------
    "chemistry": "T",
    "chem_params%diffusion": "T",
    "chem_params%reactions": "F" if args.no_reactions else "T",
    "chem_params%gamma_method": 2,  # Use cp/cv for gamma
    "chem_params%transport_model": 2,  # Unity Lewis number
    "cantera_file": ctfile,
    
    # -------------------------------------------------------------------------
    # Output
    # -------------------------------------------------------------------------
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "parallel_io": "T",
    "cons_vars_wrt": "T",
    "chem_wrt_T": "T",  # Write temperature field
}

# =============================================================================
# SPECIES MASS FRACTIONS
# =============================================================================
# Patch 1: Air (O2 + N2)
for i in range(1, num_species + 1):
    case[f"patch_icpp(1)%Y({i})"] = 0.0
case[f"patch_icpp(1)%Y({idx_O2})"] = Y_O2_right
case[f"patch_icpp(1)%Y({idx_N2})"] = Y_N2_right

# Patch 2: Fuel (H2)
for i in range(1, num_species + 1):
    case[f"patch_icpp(2)%Y({i})"] = 0.0
case[f"patch_icpp(2)%Y({idx_H2})"] = Y_H2_left

# =============================================================================
# OUTPUT
# =============================================================================
if __name__ == "__main__":
    print(json.dumps(case))
