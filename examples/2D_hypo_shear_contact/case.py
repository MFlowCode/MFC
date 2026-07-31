#!/usr/bin/env python3
"""Hypoelastic shear Riemann problem with a swept tau_xy jump.

Normalized density contact (rho 1|2 at x = 0.5) with uniform pressure p0 = 1,
zero normal velocity, uniform tangential velocity v0, tau_xx = 0, and a
tangential-stress jump tau_xy = 0|A, which launches shear waves. The
amplitude-order regression (`Convergence -> HypoShearContact` in ./mfc.sh
test) sweeps A at fixed grid: HLLD's momentum and energy fluxes carry the same
tangential double-star state, cancelling the velocity-linear mismatch and
leaving an O(A^2) pressure response; HLLC's mismatched weights leave an O(A)
response. The uniform v0 /= 0 is essential — with v0 = 0 both solvers respond
at O(A^2).

Both fluids are physically identical (gamma = 1.4, pi_inf = 0, G = 1) because
hypoelastic HLLD validation currently requires exactly two components.
"""

import argparse
import json
import math

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--mfc", type=json.loads, default={})
parser.add_argument("--solver", choices=("hllc", "hlld"), default="hlld")
parser.add_argument("--amp", type=float, default=1.0e-3, help="Right-state tau_xy")
parser.add_argument("-N", type=int, default=64)
parser.add_argument("--ny", type=int, default=8)
parser.add_argument("--steps", type=int, default=1)
parser.add_argument("--v0", type=float, default=1.0)
parser.add_argument("--cfl", type=float, default=0.5)
args = parser.parse_args()

gamma = 1.4
p0 = 1.0
G = 1.0
rho_left = 1.0
rho_right = 2.0
eps_alpha = 1.0e-12

dx = 1.0 / args.N
c_max = math.sqrt((gamma * p0 + 4.0 * G / 3.0) / min(rho_left, rho_right))
dt = args.cfl * dx / c_max

print(
    json.dumps(
        {
            "run_time_info": "F",
            "x_domain%beg": 0.0,
            "x_domain%end": 1.0,
            "y_domain%beg": 0.0,
            "y_domain%end": 1.0,
            "m": args.N - 1,
            "n": args.ny - 1,
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": args.steps,
            "t_step_save": args.steps,
            "num_patches": 2,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 2,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 1,
            "recon_type": 2,
            "muscl_order": 1,
            "muscl_lim": 1,
            "riemann_solver": {"hllc": 2, "hlld": 4}[args.solver],
            "wave_speeds": 1,
            "avg_state": 2,
            "riemann_hypo_ADC": "F",
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -1,
            "bc_y%end": -1,
            "format": 2,
            "precision": 2,
            "parallel_io": "F",
            "hypoelasticity": "T",
            # Required with hypoelasticity: the stress sources are finite-differenced
            # with fd_order-sized coefficient arrays (m_hypoelastic.fpp).
            "fd_order": 4,
            "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
            "fluid_pp(1)%pi_inf": 0.0,
            "fluid_pp(1)%G": G,
            "fluid_pp(2)%gamma": 1.0 / (gamma - 1.0),
            "fluid_pp(2)%pi_inf": 0.0,
            "fluid_pp(2)%G": G,
            # Left/background state.
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%y_centroid": 0.5,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%length_y": 1.0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": args.v0,
            "patch_icpp(1)%pres": p0,
            "patch_icpp(1)%tau_e(1)": 0.0,
            "patch_icpp(1)%tau_e(2)": 0.0,
            "patch_icpp(1)%tau_e(3)": 0.0,
            "patch_icpp(1)%alpha_rho(1)": rho_left * (1.0 - eps_alpha),
            "patch_icpp(1)%alpha(1)": 1.0 - eps_alpha,
            "patch_icpp(1)%alpha_rho(2)": rho_left * eps_alpha,
            "patch_icpp(1)%alpha(2)": eps_alpha,
            # Right state: density and tau_xy jump, same material.
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%x_centroid": 0.75,
            "patch_icpp(2)%y_centroid": 0.5,
            "patch_icpp(2)%length_x": 0.5,
            "patch_icpp(2)%length_y": 1.0,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": args.v0,
            "patch_icpp(2)%pres": p0,
            "patch_icpp(2)%tau_e(1)": 0.0,
            "patch_icpp(2)%tau_e(2)": args.amp,
            "patch_icpp(2)%tau_e(3)": 0.0,
            "patch_icpp(2)%alpha_rho(1)": rho_right * eps_alpha,
            "patch_icpp(2)%alpha(1)": eps_alpha,
            "patch_icpp(2)%alpha_rho(2)": rho_right * (1.0 - eps_alpha),
            "patch_icpp(2)%alpha(2)": 1.0 - eps_alpha,
        }
    )
)
