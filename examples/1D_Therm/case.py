#!/usr/bin/env python3
import json
import math
import cantera as ct

Lx = 0.05

ctfile = "h2o2.yaml"
sol_L = ct.Solution(ctfile)
# Initialize domain with pure H2 (setting temp to match the left wall initially)
sol_L.TPY = 600, ct.one_atm, "H2:1.0"

# Configuring case dictionary
case = {    
    "run_time_info": "T",
    "x_domain%beg": 0.0,
    "x_domain%end": Lx,
    "m": 199,
    "n": 0,     # Set n and p to 0 for 1D
    "p": 0,
    "dt": 8.0e-08,
    "t_step_start": 0,
    "t_step_stop": 20000,
    "t_step_save": 1000,
    "num_patches": 1,
    "model_eqns": 2,
    "alt_soundspeed": "F",
    "num_fluids": 1,
    "mpp_lim": "F",
    "mixture_err": "T",
    "time_stepper": 3,
    "mp_weno": "F",
    "weno_order": 5,
    "weno_eps": 1e-16,
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    
    # --- 1D Boundary Conditions ---
    # Left Wall: 600 K
    "bc_x%beg": -16, 
    "bc_x%isothermal_in": "T",
    "bc_x%Twall_in" : 600.0,
    
    # Right Wall: 900 K
    "bc_x%end": -16,
    "bc_x%isothermal_out": "T",
    "bc_x%Twall_out" : 900.0,
    
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "parallel_io": "F",
    
    # --- Physics & Chemistry ---
    "chemistry": "T" ,
    "chem_params%diffusion": "T",
    "chem_params%reactions": "F",
    "cantera_file": ctfile,
    "chem_wrt_T": "T",
    "viscous": "T",
    "fluid_pp(1)%Re(1)": 100000,
    "fluid_pp(1)%gamma": 1.0e00 / (1.4e00 - 1.0e00),
    "fluid_pp(1)%pi_inf": 0.0e00,
    
    # --- 1D Initial Patch ---
    "patch_icpp(1)%geometry": 1,
    "patch_icpp(1)%hcid": 191,
    "patch_icpp(1)%x_centroid": Lx/2,
    "patch_icpp(1)%length_x": Lx,
    "patch_icpp(1)%vel(1)": 0,
    "patch_icpp(1)%pres": 101325,
    "patch_icpp(1)%alpha_rho(1)": 1.00,
    "patch_icpp(1)%alpha(1)": 1,
}

# Automatically append the mass fractions for all species in the Cantera file
for i in range(len(sol_L.Y)):
    case[f"chem_wrt_Y({i + 1})"] = "T"
    case[f"patch_icpp(1)%Y({i+1})"] = sol_L.Y[i]
   
if __name__ == "__main__":
    print(json.dumps(case, indent=4))
