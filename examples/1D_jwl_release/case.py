"""
Isentropic release of a JWL fluid: a Riemann problem between two states of the same fluid.
Everything left of the contact keeps the left state's entropy, so the rarefaction fan and the
left star state must lie on the JWL isentrope through (rho0, p0), which is closed form.
"""

import argparse
import json
import math

parser = argparse.ArgumentParser(description="1D JWL isentropic release")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT")
parser.add_argument("-N", type=int, default=400)
parser.add_argument("--cfl", type=float, default=0.4)
args = parser.parse_args()

rho0, p0, rho_r, p_r = 1.0, 1.0, 0.3, 0.1
jwl = {"a": 6.0, "b": 0.15, "r1": 4.0, "r2": 1.0, "omega": 0.3, "rho0": rho0}
N, L, T_end = args.N, 1.0, 0.15
c_max = math.sqrt(((1.0 + jwl["omega"]) * p0 + jwl["a"] + jwl["b"]) / rho0)  # generous bound on c + |u|
dt = args.cfl * (L / N) / c_max
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
    "fluid_pp(1)%eos": "jwl",
    **{f"fluid_pp(1)%jwl_{k}": v for k, v in jwl.items()},
}
for pid, (x_c, rho, pres) in enumerate([(0.25, rho0, p0), (0.75, rho_r, p_r)], start=1):
    case.update(
        {
            f"patch_icpp({pid})%geometry": 1,
            f"patch_icpp({pid})%x_centroid": x_c,
            f"patch_icpp({pid})%length_x": 0.5 * L,
            f"patch_icpp({pid})%alpha_rho(1)": rho,
            f"patch_icpp({pid})%alpha(1)": 1.0,
            f"patch_icpp({pid})%vel(1)": 0.0,
            f"patch_icpp({pid})%pres": pres,
        }
    )
print(json.dumps(case))
