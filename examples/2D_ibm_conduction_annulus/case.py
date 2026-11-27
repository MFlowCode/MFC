#!/usr/bin/env python3
"""
Steady conduction across a 2D annulus between two isothermal immersed boundaries: image-point interpolation of the wall temperature.

The inner body is a circle (radius R = 1, patch_ib%geometry = 2, Twall = 2), the outer body is a flat ring read from ring.stl whose
inner edge is the circle r = 2 R (geometry 5, Twall = 3). The gas between them is quiescent with constant conductivity, so the steady
temperature is exactly T(r) = T1 + (T2 - T1) ln(r / R) / ln(r2 / R). The initial temperature is the exact profile on 10 shells and the run is advanced for 8
slowest conduction times. --interp selects the image-point weights of the ghost-cell method: 1 inverse distance, 2 triangular-shaped cloud
(ib_ip_interp = 2, the solver default).

    ./mfc.sh run examples/2D_ibm_conduction_annulus/case.py -n 1 --gpu                    # inverse distance
    ./mfc.sh run examples/2D_ibm_conduction_annulus/case.py -n 1 --gpu -- --interp 2      # TSC
    python3 examples/2D_ibm_conduction_annulus/check_annulus.py <case directory>

Defaults are a CI-sized run: 12 cells per radius (64 x 64 cells). Convergence runs use --cells_per_R 10, 20, 40 (README).
Nondimensional units: gamma = 1.4, c_v = 1, rho = 1, thermal diffusivity 0.1. For a physical annulus (R = 50 um, 600 K and 300 K walls,
air at 1 atm) scale lengths by R and temperatures by T2 - T1 and use a temperature-dependent conductivity (needs chemistry = T).
"""

import argparse
import json
import math

parser = argparse.ArgumentParser(prog="2D_ibm_conduction_annulus")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC toolchain state.")
parser.add_argument("--cells_per_R", type=int, default=12, help="Cells per inner radius.")
parser.add_argument("--interp", type=int, default=1, choices=(1, 2), help="Image-point weights: 1 inverse distance, 2 TSC.")
parser.add_argument("--steps", type=int, default=0, help="Number of time steps (0 = 8 slowest conduction times).")
args = parser.parse_args()

gam, cv = 1.4, 1.0
T1, T2 = 2.0, 3.0  # inner and outer wall temperature
R, r2 = 1.0, 2.0  # ring.stl has its inner edge at r = 2
alpha = 0.1  # thermal diffusivity k / (rho cv)
rho = 1.0  # density at the mean wall temperature; the pressure follows from it
p0 = (gam - 1.0) * cv * rho * 0.5 * (T1 + T2)
c0 = (gam * p0 / rho) ** 0.5
NSHELL = 10  # concentric shells of the initial profile (patch_icpp holds at most 10 patches)


def T_exact(r):
    return T1 + (T2 - T1) * math.log(r / R) / math.log(r2 / R)


dx = R / args.cells_per_R
half = math.ceil((r2 + max(0.3 * R, 8 * dx)) / dx)  # cells from the centre to the domain edge
L = half * dx
tau = ((r2 - R) / math.pi) ** 2 / alpha  # slowest conduction mode of the gap
dt = min(0.25 * dx / c0, 0.15 * dx**2 / alpha)
steps = args.steps or math.ceil(8.0 * tau / dt)

case = {
    "run_time_info": "F",
    "x_domain%beg": -L,
    "x_domain%end": L,
    "y_domain%beg": -L,
    "y_domain%end": L,
    "m": 2 * half - 1,
    "n": 2 * half - 1,
    "p": 0,
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": steps,
    "t_step_save": steps,
    "num_patches": NSHELL,
    "model_eqns": 2,
    "alt_soundspeed": "F",
    "num_fluids": 1,
    "mpp_lim": "F",
    "mixture_err": "F",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1.0e-16,
    "mapped_weno": "F",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "bc_y%beg": -3,
    "bc_y%end": -3,
    "fd_order": 2,
    "viscous": "F",
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "F",
    "parallel_io": "T",
    "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(1)%eos": "ideal_gas",
    "fluid_pp(1)%cv": cv,
    "fluid_pp(1)%k_therm": alpha * rho * cv,
    # Isothermal no-slip walls: the inner circle and the outer ring
    "ib": "T",
    "num_ibs": 2,
    "patch_ib(1)%geometry": 2,
    "patch_ib(1)%x_centroid": 0.0,
    "patch_ib(1)%y_centroid": 0.0,
    "patch_ib(1)%radius": R,
    "patch_ib(1)%slip": "F",
    "patch_ib(1)%thermal_bc": 1,
    "patch_ib(1)%Twall": T1,
    "num_stl_models": 1,
    "stl_models(1)%model_filepath": "ring.stl",
    "stl_models(1)%model_threshold": 0.5,
    "patch_ib(2)%geometry": 5,
    "patch_ib(2)%model_id": 1,
    "patch_ib(2)%x_centroid": 0.0,
    "patch_ib(2)%y_centroid": 0.0,
    "patch_ib(2)%slip": "F",
    "patch_ib(2)%thermal_bc": 1,
    "patch_ib(2)%Twall": T2,
}
# Initial temperature: the exact profile sampled at the middle of NSHELL concentric shells of the gap (constant patches, outermost shell
# first; each smaller circle alters the earlier ones), at uniform pressure. Starting from the exact profile instead of a uniform
# temperature leaves the gas almost no acoustic ringing to damp, which an inviscid run would otherwise keep as a 1 % temperature error.
for j in range(1, NSHELL + 1):
    shell = NSHELL + 1 - j  # 1 = innermost
    T_shell = T_exact(R + (r2 - R) * (shell - 0.5) / NSHELL)
    pre = f"patch_icpp({j})%"
    case.update(
        {
            pre + "geometry": 3 if j == 1 else 2,
            pre + "x_centroid": 0.0,
            pre + "y_centroid": 0.0,
            pre + "vel(1)": 0.0,
            pre + "vel(2)": 0.0,
            pre + "pres": p0,
            pre + "alpha_rho(1)": p0 / ((gam - 1.0) * cv * T_shell),
            pre + "alpha(1)": 1.0,
        }
    )
    if j == 1:
        case.update({pre + "length_x": 2 * L, pre + "length_y": 2 * L})
    else:
        case[pre + "radius"] = R + (r2 - R) * shell / NSHELL
        case.update({pre + f"alter_patch({i})": "T" for i in range(1, j)})
case["ib_ip_interp"] = args.interp  # always written: the solver default is 2

if __name__ == "__main__":
    print(json.dumps(case))
