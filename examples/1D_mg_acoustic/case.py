"""
Right-moving acoustic pulse in a single Mie-Gruneisen fluid at its reference state.
The pulse is a simple wave (drho, dp = c^2 drho, du = c drho/rho0) so only the right-going
characteristic carries it; the harness measures its speed against the general analytic c.
"""

import argparse
import json
import math

parser = argparse.ArgumentParser(description="1D Mie-Gruneisen acoustic pulse")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT")
parser.add_argument("-N", type=int, default=200)
parser.add_argument("--cfl", type=float, default=0.4)
args = parser.parse_args()

rho0, p0, c0, s, gruneisen = 1.0, 1.0, 1.0, 1.5, 0.4
c = math.sqrt(c0**2 + (1.0 + gruneisen) * p0 / rho0)  # the frozen speed at the reference state
amp, x0, width = 1.0e-4, 0.3, 0.05
N, L, T_end = args.N, 1.0, 0.4
dt = args.cfl * (L / N) / c
Nt = math.ceil(T_end / dt)
dt = T_end / Nt
pulse = f"{amp}*exp(-((x - {x0})/{width})**2)"

print(
    json.dumps(
        {
            "run_time_info": "F",
            "x_domain%beg": 0.0,
            "x_domain%end": L,
            "m": N - 1,
            "n": 0,
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": Nt,
            "t_step_save": Nt,
            "num_patches": 1,
            "model_eqns": 2,
            "num_fluids": 1,
            "time_stepper": 3,
            "recon_type": "weno",
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "mapped_weno": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "F",
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%length_x": L,
            "patch_icpp(1)%alpha_rho(1)": f"{rho0} + {pulse}",
            "patch_icpp(1)%alpha(1)": 1.0,
            "patch_icpp(1)%vel(1)": f"{c}/{rho0}*{pulse}",
            "patch_icpp(1)%pres": f"{p0} + {c}**2*{pulse}",
            "fluid_pp(1)%eos": "mie_gruneisen",
            "fluid_pp(1)%mg_rho0": rho0,
            "fluid_pp(1)%mg_c0": c0,
            "fluid_pp(1)%mg_s": s,
            "fluid_pp(1)%mg_gruneisen": gruneisen,
        }
    )
)
