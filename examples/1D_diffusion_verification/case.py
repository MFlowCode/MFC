#!/usr/bin/env python3
"""
1D verification of mixture-averaged species diffusion and heat conduction (h2o2.yaml, reactions off, viscosity off).

--mode erf_species  trace H2 in N2, a step of mass fraction 1e-6: exact solution (A/2) erfc((x - x0)/sqrt(4 D t))
--mode erf_heat     pure N2, a step of 0.01 K at 300 K: the same erfc with the thermal diffusivity
--mode mixing       H2/N2 at 300 K meets O2/H2O/N2 at 1000 K at constant pressure (the interface starts as a step)
--mode stab         periodic, a one-cell temperature spike; dt is set by chem_params%transport_dt. Run with --cfl above the
                    limit and the grid-scale mode grows (the default is stable)

    ./mfc.sh run examples/1D_diffusion_verification/case.py -t pre_process simulation post_process -n 1 -- --mode erf_species --n 256
    python3 examples/1D_diffusion_verification/check_diffusion.py examples/1D_diffusion_verification --case-args '--mode erf_species --n 256'

The erf runs end when the front width sqrt(4 D t) is --width of the domain (6.4 cells at --n 128), so a grid sweep is a
convergence study. Defaults are a CI-sized run; the README gives the grids of its tables.
"""

import argparse
import json
import math

import cantera as ct

parser = argparse.ArgumentParser(prog="1D_diffusion_verification", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC toolchain state.")
parser.add_argument("--mode", default="erf_species", choices=("erf_species", "erf_heat", "mixing", "stab"))
parser.add_argument("--n", type=int, default=128, help="Cells.")
parser.add_argument("--width", type=float, default=0.05, help="erf modes: front width sqrt(4 D t_end) at the end, as a fraction of the domain.")
parser.add_argument("--tm", type=int, default=1, choices=(1, 2), help="chem_params%%transport_model: 1 mixture-averaged, 2 unity Lewis.")
parser.add_argument("--cfl", type=float, default=0.5, help="stab mode: cfl_target (the transport limit is at about 0.6).")
parser.add_argument("--mixture", default="n2", choices=("n2", "mix"), help="stab mode: N2 at 300 K or H2/O2/H2O/N2 at 1050 K.")
parser.add_argument("--bc_out", type=int, default=-3, help="mixing mode: right boundary (-3 extrapolation, -8 characteristic outflow at p0).")
args = parser.parse_args()

MECH = "h2o2.yaml"
P0 = 101325.0
gas = ct.Solution(MECH)
gas.transport_model = "mixture-averaged"


def state(T, X):
    gas.TPX = T, P0, X
    return gas.density, gas.Y.copy()


# Two or three uniform states side by side: (x_end of the region, T, X); the first region starts at 0.
if args.mode == "erf_species":
    A = 1.0e-6
    L = 1.0e-3
    T0 = 300.0
    gas.TPX = T0, P0, "N2:1"
    gas.TPY = T0, P0, {"H2": A, "N2": 1.0 - A}
    D = gas.mix_diff_coeffs[gas.species_index("H2")]
    regions = [(L / 2, T0, {"H2": A, "N2": 1.0 - A}), (L, T0, {"N2": 1.0})]
elif args.mode == "erf_heat":
    A = 0.01
    L = 1.0e-3
    T0 = 300.0
    gas.TPX = T0, P0, "N2:1"
    D = gas.thermal_conductivity / (gas.density * gas.cp_mass)
    regions = [(L / 2, T0 + A, {"N2": 1.0}), (L, T0, {"N2": 1.0})]
elif args.mode == "mixing":
    L = 3.0e-3
    regions = [(1.0e-3, 300.0, "H2:0.3,N2:0.7"), (L, 1000.0, "O2:0.2,H2O:0.2,N2:0.6")]
else:
    X = "N2:1" if args.mixture == "n2" else "H2:0.15,O2:0.1,H2O:0.1,N2:0.65"
    T0 = 300.0 if args.mixture == "n2" else 1050.0
    gas.TPX = T0, P0, X
    alpha_v = gas.thermal_conductivity / (gas.density * gas.cv_mass)  # conduction diffuses T at lambda / (rho c_v)
    present = gas.Y > 0.0  # a species at Y = 0 has no flux, so only the species present limit the step
    Dmax = alpha_v if args.tm == 2 else max(gas.mix_diff_coeffs[present].max(), alpha_v)
    # dx = Dmax / (10 c): the diffusion limit is ten times tighter than the acoustic one
    L = args.n * Dmax / (10.0 * gas.sound_speed)
    regions = None

dx = L / args.n
c = 650.0 if args.mode == "mixing" else gas.sound_speed
case = {
    "run_time_info": "T",
    "x_domain%beg": 0.0,
    "x_domain%end": L,
    "m": args.n - 1,
    "n": 0,
    "p": 0,
    "parallel_io": "F",
    "model_eqns": 2,
    "num_fluids": 1,
    "mpp_lim": "F",
    "mixture_err": "F",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1e-16,
    "mapped_weno": "T",
    "mp_weno": "T",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    "viscous": "F",
    "chemistry": "T",
    "chem_params%diffusion": "T",
    "chem_params%reactions": "F",
    "chem_params%transport_model": args.tm,
    "cantera_file": MECH,
    "fluid_pp(1)%gamma": 1.0 / (1.4 - 1.0),
    "fluid_pp(1)%eos": "stiffened_gas",
    "fluid_pp(1)%pi_inf": 0.0,
    "format": "silo",
    "precision": "double",
    "prim_vars_wrt": "T",
    "chem_wrt_T": "T",
}
for i in range(gas.n_species):
    case[f"chem_wrt_Y({i + 1})"] = "T"

if args.mode == "stab":
    # A spike of 0.1 K in the middle cell of 64 on a periodic domain; 400 steps at the transport-limited dt.
    dt_tr = args.cfl * dx**2 / Dmax
    t_stop = 400 * dt_tr
    regions = [((args.n // 2) * dx, T0, X), ((args.n // 2 + 1) * dx, T0 + 0.1, X), (L, T0, X)]
    case.update(
        {
            "bc_x%beg": -1,
            "bc_x%end": -1,
            "cfl_adap_dt": "T",
            "cfl_target": args.cfl,
            "n_start": 0,
            "t_stop": t_stop,
            "t_save": t_stop,
            "chem_params%transport_dt": "T",
        }
    )
else:
    t_end = (args.width * L) ** 2 / (4.0 * D) if args.mode != "mixing" else 2.0e-5
    steps = math.ceil(t_end / (0.4 * dx / c))
    save = steps // 2 if args.mode == "mixing" else steps
    steps = save * (2 if args.mode == "mixing" else 1)
    case.update(
        {
            "bc_x%beg": -2 if args.mode == "mixing" else -3,
            "bc_x%end": args.bc_out if args.mode == "mixing" else -3,
            "dt": t_end / steps,
            "t_step_start": 0,
            "t_step_stop": steps,
            "t_step_save": save,
            "t_step_print": save,
        }
    )
    if args.mode == "mixing" and args.bc_out == -8:
        case.update({"bc_x%grcbc_out": "T", "bc_x%pres_out": P0})

case["num_patches"] = len(regions)
beg = 0.0
for k, (end, T, X) in enumerate(regions, start=1):
    rho, Y = state(T, X)
    case.update(
        {
            f"patch_icpp({k})%geometry": 1,
            f"patch_icpp({k})%x_centroid": 0.5 * (beg + end),
            f"patch_icpp({k})%length_x": end - beg,
            f"patch_icpp({k})%vel(1)": 0.0,
            f"patch_icpp({k})%pres": P0,
            f"patch_icpp({k})%alpha(1)": 1.0,
            f"patch_icpp({k})%alpha_rho(1)": rho,
        }
    )
    for i, y in enumerate(Y):
        case[f"patch_icpp({k})%Y({i + 1})"] = float(y)
    beg = end

if __name__ == "__main__":
    print(json.dumps(case))
