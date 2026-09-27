"""
Shock-initiated detonation of a stiffened-gas reactant into JWL products loaded from a material file.
The products' Q sets both phase qv values, so the probes at 2, 3 and 4 cm should time the front at the
CJ speed of the fit, 7404 m/s.
"""

import argparse
import json

parser = argparse.ArgumentParser(description="1D JWL detonation")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT")
parser.add_argument("-N", type=int, default=1500)
parser.add_argument("--tend", type=float, default=6.0e-6)
args = parser.parse_args()

rho0, p0, p_init, x_init = 1600.0, 1.0e5, 5.0e9, 0.004
Gamma, Pi = 3.0, 6.0e8
# A trace of expanded product keeps every phase inside its equation of state.
eps, rho_trace = 1.0e-6, 1.0
N, L = args.N, 0.06
dt = 0.2 * (L / N) / 12000.0
Nt = int(args.tend / dt)


def patch(i, x_centroid, length, pres):
    return {
        f"patch_icpp({i})%geometry": 1,
        f"patch_icpp({i})%x_centroid": x_centroid,
        f"patch_icpp({i})%length_x": length,
        f"patch_icpp({i})%vel(1)": 0.0,
        f"patch_icpp({i})%pres": pres,
        f"patch_icpp({i})%alpha_rho(1)": rho0 * (1.0 - eps),
        f"patch_icpp({i})%alpha_rho(2)": rho_trace * eps,
        f"patch_icpp({i})%alpha(1)": 1.0 - eps,
        f"patch_icpp({i})%alpha(2)": eps,
    }


case = {
    "run_time_info": "F",
    "x_domain%beg": 0.0,
    "x_domain%end": L,
    "m": N,
    "n": 0,
    "p": 0,
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": Nt,
    "t_step_save": Nt,
    "parallel_io": "F",
    "model_eqns": 2,
    "num_fluids": 2,
    "num_patches": 2,
    "mpp_lim": "T",
    "mixture_err": "T",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1.0e-16,
    "mapped_weno": "T",
    "mp_weno": "T",
    "riemann_solver": "hllc",
    "wave_speeds": "direct",
    "avg_state": "arithmetic",
    "bc_x%beg": -2,
    "bc_x%end": -3,
    "reactive_burn": "T",
    "rburn%k": 5.0e6,
    "rburn%pign": 5.0e8,
    "rburn%pref": 1.0e9,
    "rburn%n": 1.0,
    "format": "silo",
    "precision": "double",
    "prim_vars_wrt": "T",
    "probe_wrt": "T",
    "fd_order": 1,
    "num_probes": 3,
    **{f"probe({i + 1})%x": x for i, x in enumerate([0.02, 0.03, 0.04])},
    **patch(1, L / 2, L, p0),
    **patch(2, x_init / 2, x_init, p_init),
    "patch_icpp(2)%alter_patch(1)": "T",
    "fluid_pp(1)%eos": "stiffened_gas",
    "fluid_pp(1)%gamma": 1.0 / (Gamma - 1.0),
    "fluid_pp(1)%pi_inf": Gamma * Pi / (Gamma - 1.0),
    "fluid_pp(2)%material_file": "jwl_products.yaml",
}

if __name__ == "__main__":
    print(json.dumps(case))
