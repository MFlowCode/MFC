#!/usr/bin/env python3
# Freely propagating stoichiometric H2/air flame (h2o2.yaml, Tu = 300 K, 1 atm, mixture-averaged transport, no Soret) started from the
# Cantera FreeFlame profile and compared with the Cantera flame speed by compare_flame_speed.py.
#
# The profile is computed with Cantera at the first run and written to IC/ (read through hcid 170); the flame starts at x = 6.5 mm and
# moves in -x into gas at rest. Reflective wall at x = 0, GRCBC subsonic outflow (bc -8) with the outlet pressure held at 1 atm.
# Without the pressure target the closed domain pressurizes (+13 % unburnt density by 1 ms); extrapolation at both ends lets the whole
# domain accelerate.
#
# Production grids: dx = 40, 20, 10, 5 um (--dx_um; 8 to 66 cells per thermal thickness of 329 um), run to 1.5 ms (--tend).
# Reference: Cantera FreeFlame, S_L = 2.3317 m/s (width 2 cm, refine slope = curve = 0.02).
import argparse
import json
import os

import cantera as ct
import numpy as np

parser = argparse.ArgumentParser(prog="1D_h2_air_flame_speed", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC's toolchain's internal state.")
parser.add_argument("--dx_um", type=float, default=40.0, help="Cell size [um].")
parser.add_argument("--tend", type=float, default=1.5e-3, help="End time [s].")
parser.add_argument("--cfl", type=float, default=0.25, help="dt = cfl dx / 1100 m/s (1100 m/s bounds |u| + c in the burnt gas).")
parser.add_argument("--substeps", type=int, default=2, help="chem_params%%reaction_substeps.")
args = parser.parse_args()

here = os.path.dirname(os.path.abspath(__file__))
ctfile = "h2o2.yaml"
L = 8.0e-3  # domain length [m]
X_FLAME = 6.5e-3  # initial front position [m]
dx = args.dx_um * 1e-6
Nx = int(round(L / dx))
dt = args.cfl * dx / 1100.0
NT = int(round(args.tend / dt))
NSAVE = 150


def write_ic():
    """Cantera FreeFlame on the cell centres, shifted to X_FLAME and to the frame of the unburnt gas (skipped if IC/ is current)."""
    key = os.path.join(here, "IC", ".cache_key.json")
    if os.path.exists(key) and json.load(open(key)) == {"dx_um": args.dx_um}:
        return
    gas = ct.Solution(ctfile)
    gas.TPX = 300.0, ct.one_atm, "H2:2, O2:1, N2:3.76"
    f = ct.FreeFlame(gas, width=0.02)
    f.transport_model = "mixture-averaged"
    f.soret_enabled = False
    f.set_refine_criteria(ratio=2.0, slope=0.02, curve=0.02, prune=0.0)
    f.solve(loglevel=0, auto=True)
    z, T = f.grid, f.T
    zf = np.interp(0.5 * (T[0] + T[-1]), T, z)
    xc = (np.arange(Nx) + 0.5) * dx
    zq = xc - X_FLAME + zf
    cols = [np.interp(zq, z, f.density), np.interp(zq, z, f.velocity) - f.velocity[0], np.full_like(xc, f.P), np.ones_like(xc)]
    cols += [np.interp(zq, z, y) for y in f.Y]
    os.makedirs(os.path.join(here, "IC"), exist_ok=True)
    for k, v in enumerate(cols, start=1):
        np.savetxt(os.path.join(here, "IC", f"prim.{k}.00.000000.dat"), np.column_stack([xc, v]), fmt="%.17e")
    # reference values for compare_flame_speed.py: S_L [m/s], thermal thickness [m], unburnt density, Y_H2 unburnt and burnt
    iH2 = gas.species_index("H2")
    np.savetxt(os.path.join(here, "IC", "reference.dat"), [f.velocity[0], (T[-1] - T[0]) / np.gradient(T, z).max(), f.density[0], f.Y[iH2, 0], f.Y[iH2, -1]])
    json.dump({"dx_um": args.dx_um}, open(key, "w"))


if __name__ == "__main__":
    write_ic()

case = {
    # Logistics
    "run_time_info": "F",
    # Computational Domain Parameters
    "x_domain%beg": 0.0,
    "x_domain%end": L,
    "m": Nx - 1,
    "n": 0,
    "p": 0,
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": NT,
    "t_step_save": max(1, NT // NSAVE),
    "t_step_print": max(1, NT // NSAVE),
    "parallel_io": "F",
    # Simulation Algorithm Parameters
    "model_eqns": "5eq",
    "num_fluids": 1,
    "num_patches": 1,
    "mpp_lim": "F",
    "mixture_err": "F",
    "time_stepper": "rk3",
    "weno_order": 5,
    "weno_eps": 1e-16,
    "mapped_weno": "T",
    "riemann_solver": "hllc",
    "wave_speeds": "direct",
    "avg_state": "arithmetic",
    "bc_x%beg": -2,
    "bc_x%end": -8,
    "bc_x%grcbc_out": "T",
    "bc_x%pres_out": 101325.0,
    "viscous": "F",
    # Formatted Database Files Structure Parameters
    "format": "silo",
    "precision": "double",
    "prim_vars_wrt": "T",
    "chem_wrt_T": "T",
    "files_dir": os.path.join(here, "IC"),
    "file_extension": "000000",
    # Patch 1: profile from IC/
    "patch_icpp(1)%geometry": 1,
    "patch_icpp(1)%hcid": 170,
    "patch_icpp(1)%x_centroid": 0.5 * L,
    "patch_icpp(1)%length_x": L,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%pres": 101325.0,
    "patch_icpp(1)%alpha(1)": 1,
    "patch_icpp(1)%alpha_rho(1)": 1.0,
    # Fluids Physical Parameters
    "fluid_pp(1)%gamma": 1.0 / (1.4 - 1.0),
    "fluid_pp(1)%eos": "stiffened_gas",
    "fluid_pp(1)%pi_inf": 0,
    # Chemistry
    "chemistry": "T",
    "chem_params%diffusion": "T",
    "chem_params%reactions": "T",
    "chem_params%transport_model": 1,
    "chem_params%reaction_substeps": args.substeps,
    "cantera_file": ctfile,
}
for i in range(ct.Solution(ctfile).n_species):
    case[f"chem_wrt_Y({i + 1})"] = "T"

if __name__ == "__main__":
    print(json.dumps(case))
