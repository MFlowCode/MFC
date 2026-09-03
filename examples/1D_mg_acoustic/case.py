"""
Right-moving acoustic pulse in a single Mie-Gruneisen fluid at its reference state.
A rectangle patch carries a simple-wave perturbation (drho, dp = c^2 drho, du = c drho/rho0), so
only the right-going characteristic is excited; the harness tracks the centroid of drho against
the general analytic c. Patches rather than an analytic IC: an IC expression is compiled in, and
every distinct one costs the test suite a full rebuild.
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
amp, x0, width = 1.0e-4, 0.25, 0.1
N, L, T_end = args.N, 1.0, 0.3
dt = args.cfl * (L / N) / c
Nt = math.ceil(T_end / dt)
dt = T_end / Nt

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
            "num_patches": 2,
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%length_x": L,
            "patch_icpp(1)%alpha_rho(1)": rho0,
            "patch_icpp(1)%alpha(1)": 1.0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%pres": p0,
            "patch_icpp(2)%geometry": 1,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": x0,
            "patch_icpp(2)%length_x": width,
            "patch_icpp(2)%alpha_rho(1)": rho0 + amp,
            "patch_icpp(2)%alpha(1)": 1.0,
            "patch_icpp(2)%vel(1)": c / rho0 * amp,
            "patch_icpp(2)%pres": p0 + c**2 * amp,
            "fluid_pp(1)%eos": "mie_gruneisen",
            "fluid_pp(1)%mg_rho0": rho0,
            "fluid_pp(1)%mg_c0": c0,
            "fluid_pp(1)%mg_s": s,
            "fluid_pp(1)%mg_gruneisen": gruneisen,
        }
    )
)
