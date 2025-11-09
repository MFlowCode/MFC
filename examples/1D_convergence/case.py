#!/usr/bin/env python3
import math
import json
import argparse

# Parsing command line arguments
parser = argparse.ArgumentParser(description="Generate JSON case configuration for two-fluid convergence simulation.")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC's toolchain's internal state.")
parser.add_argument("--order", type=int, default=5, help="WENO order (default: 5)")
parser.add_argument("--meqns", type=int, default=2, help="Model equations (default: 2 (five-equation model))")
parser.add_argument("--rs", type=int, default=2, help="Riemann solver (default: 2 (HLLC))")
parser.add_argument("-N", type=int, default=1024, help="Number of grid points (default: 1024)")

args = parser.parse_args()

# Numerical setup
Nx = args.N - 1
dx = 1.0 / (1.0 * (Nx + 1))

Tend, Nt = 5.0, 200000
mydt = Tend / (1.0 * Nt)

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0e00,
            "x_domain%end": 1.0e00,
            "m": Nx,
            "n": 0,
            "p": 0,
            "dt": mydt,
            "t_step_start": 0,
            "t_step_stop": int(Nt),
            "t_step_save": int(math.ceil(Nt / 10.0)),
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": args.meqns,
            "alt_soundspeed": "F",
            "num_fluids": 2,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 3,
            "weno_order": args.order,
            "weno_eps": dx**2,
            "weno_Re_flux": "F",
            "weno_avg": "F",
            "mapped_weno": "F",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": args.rs,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -1,
            "bc_x%end": -1,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "F",
            # Patch 1 L
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%vel(1)": 1.0,
            "patch_icpp(1)%pres": 1.0,
            "patch_icpp(1)%alpha_rho(1)": f"0.5 - 0.5*sin(2*pi*x)",
            "patch_icpp(1)%alpha(1)": f"0.5 - 0.5*sin(2*pi*x)",
            "patch_icpp(1)%alpha_rho(2)": f"0.5 + 0.5*sin(2*pi*x)",
            "patch_icpp(1)%alpha(2)": f"0.5 + 0.5*sin(2*pi*x)",
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (1.4 - 1.0e00),
            "fluid_pp(1)%pi_inf": 0.0,
            "fluid_pp(2)%gamma": 1.0e00 / (1.4 - 1.0e00),
            "fluid_pp(2)%pi_inf": 0.0,
        }
    )
)
