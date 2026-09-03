"""
Symmetric impact of two Mie-Gruneisen slabs approaching at relative speed U.
Each slab is brought to rest by a shock with particle-velocity jump U/2, so the shock state
lies exactly on the Hugoniot u_s = c0 + s u_p; the harness checks the shock speed and the
plateau density against that relation.
"""

import argparse
import json
import math

parser = argparse.ArgumentParser(description="1D Mie-Gruneisen symmetric impact")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT")
parser.add_argument("-N", type=int, default=800)
parser.add_argument("--U", type=float, default=1.0, help="closing speed of the two slabs")
parser.add_argument("--cfl", type=float, default=0.4)
args = parser.parse_args()

rho0, p0, c0, s, gruneisen = 1.0, 1.0e-3, 1.0, 1.5, 0.4
N, L, T_end = args.N, 1.0, 0.2
dt = args.cfl * (L / N) / (c0 + (s + 1.0) * args.U)
Nt = math.ceil(T_end / dt)
dt = T_end / Nt

case = {
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
    "num_patches": 2,
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
    "fluid_pp(1)%eos": "mie_gruneisen",
    "fluid_pp(1)%mg_rho0": rho0,
    "fluid_pp(1)%mg_c0": c0,
    "fluid_pp(1)%mg_s": s,
    "fluid_pp(1)%mg_gruneisen": gruneisen,
}
for pid, (x_c, vel) in enumerate([(0.25, 0.5 * args.U), (0.75, -0.5 * args.U)], start=1):
    case.update(
        {
            f"patch_icpp({pid})%geometry": 1,
            f"patch_icpp({pid})%x_centroid": x_c,
            f"patch_icpp({pid})%length_x": 0.5 * L,
            f"patch_icpp({pid})%alpha_rho(1)": rho0,
            f"patch_icpp({pid})%alpha(1)": 1.0,
            f"patch_icpp({pid})%vel(1)": vel,
            f"patch_icpp({pid})%pres": p0,
        }
    )
print(json.dumps(case))
