#!/usr/bin/env python3
import json
import argparse
import math

import cantera as ct

parser = argparse.ArgumentParser(prog="1D_Premixed_Flame", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC's toolchain's internal state.")
parser.add_argument("--no-chem", dest="chemistry", default=True, action="store_false", help="Disable chemistry.")

args = parser.parse_args()

ctfile = "h2o2.yaml"
sol_L = ct.Solution(ctfile)
sol_L.TPX = 300, 8000, "O2:2,N2:2,H2O:5"

L = 0.08
Nx = 2000
dx = L / Nx
dt = 0.25e-7
Tend = 0.60e-2

NT = int(Tend / dt)
SAVE_COUNT = 2000
NS = 2000
case = {
    "run_time_info": "T",
    "x_domain%beg": 0,
    "x_domain%end": +L,
    "m": Nx,
    "n": 0,
    "p": 0,
    "dt": float(dt),
    "t_step_start": 0,
    "t_step_stop": NT,
    "t_step_save": NS,
    "t_step_print": NS,
    "parallel_io": "F",
    "model_eqns": 2,
    "num_fluids": 1,
    "num_patches": 1,
    "mpp_lim": "F",
    "mixture_err": "F",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1e-16,
    "weno_avg": "F",
    "mapped_weno": "T",
    "mp_weno": "T",
    "riemann_solver": 2,
    "wave_speeds": 2,
    "avg_state": 1,
    "bc_x%beg": -8,
    "bc_x%end": -8,
    "viscous": "T",
    "chemistry": "T" if not args.chemistry else "T",
    "chem_params%diffusion": "T",
    "chem_params%reactions": "T",
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "patch_icpp(1)%geometry": 1,
    "patch_icpp(1)%x_centroid": L / 2,
    "patch_icpp(1)%length_x": L,
    "patch_icpp(1)%vel(1)": "0",
    "patch_icpp(1)%pres": 1.01325e5,
    "patch_icpp(1)%alpha(1)": 1,
    "patch_icpp(1)%Y(1)": "(0.02851162527578208 + (1.882823170077901*10**(-5) - 0.02851162527578208) * 0.5*(1-tanh(1250*(x-0.04))))",
    "patch_icpp(1)%Y(2)": "(0.0 + (3.3821733749391013*10**(-8) - 0.0) * 0.5*(1-tanh(1250*(x-0.04))))",
    "patch_icpp(1)%Y(4)": "(0.22626859761271706 + (0.00013139577035943594 - 0.22626859761271706) * 0.5*(1-tanh(1250*(x-0.04))))",
    "patch_icpp(1)%Y(3)": "(0.0 + (9.488275739032257*10**(-8) - 0.0) * 0.5*(1-tanh(1250*(x-0.04))))",
    "patch_icpp(1)%Y(5)": "(0.0 + (3.868433894115635*10**(-5)- 0.0) * 0.5*(1-tanh(1250*(x-0.04))))",
    "patch_icpp(1)%Y(7)": "(0.0 + (3.087047007163217*10**(-9) - 0.0) * 0.5*(1-tanh(1250*(x-0.04))))",
    "patch_icpp(1)%Y(8)": "(0.0 + (1.7396153620760586*10**(-9) - 0.0) * 0.5*(1-tanh(1250*(x-0.04))))",
    "patch_icpp(1)%Y(6)": "(0.0 + (0.2545911810163394 - 0.0) * 0.5*(1-tanh(1250*(x-0.04))))",
    "patch_icpp(1)%Y(10)": "0.0",
    "patch_icpp(1)%Y(10)": "(0.7452197771115008 + (0.7452197771115058 - 0.7452197771115008) * 0.5*(1-tanh(1250*(x-0.04))))",
    "patch_icpp(1)%alpha_rho(1)": "1.01325d0*10.0d0**(5.0d0)/(8.3144626d0*1000.0d0*(300 + (1600 - 300) * 0.5 * (1 - tanh(1250 * (x - 0.04))))* "
    + "((0.02851162527578208 + (1.882823170077901*10**(-5) - 0.02851162527578208) * 0.5*(1-tanh(1250*(x-0.04))))/2.016  + ( 3.3821733749391013*10**(-8) * 0.5*(1-tanh(1250*(x-0.04))))/1.008 +   "
    + "(0.22626859761271706 + (0.00013139577035943594 - 0.22626859761271706) * 0.5*(1-tanh(1250*(x-0.04))))/31.998 +  (9.488275739032257*10**(-8)  * 0.5*(1-tanh(1250*(x-0.04))))/16.00 +  "
    + "(0.0 + (3.868433894115635*10**(-5)- 0.0) * 0.5*(1-tanh(1250*(x-0.04))))/17.007 +  (0.0 + (3.087047007163217*10**(-9) - 0.0) * 0.5*(1-tanh(1250*(x-0.04))))/33.006  +  "
    + " (0.0 + (1.7396153620760586*10**(-9) - 0.0) * 0.5*(1-tanh(1250*(x-0.04))))/34.014+ (0.0 + (0.2545911810163394 - 0.0) * 0.5*(1-tanh(1250*(x-0.04))))/18.015 + "
    + "(0.7452197771115008 + (0.7452197771115058 - 0.7452197771115008) * 0.5*(1-tanh(1250*(x-0.04))))/28.014))   ",
    "fluid_pp(1)%gamma": 1.0e00 / (1.9326e00 - 1.0e00),
    "fluid_pp(1)%pi_inf": 0,
    "fluid_pp(1)%Re(1)": 20000000,
    "cantera_file": ctfile,
}

if __name__ == "__main__":
    print(json.dumps(case))
