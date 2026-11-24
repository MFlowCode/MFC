#!/usr/bin/env python3
"""
3D premixed H2/air flame (phi = 1, 1 atm, 300 K) hit by a vortex ring, no immersed boundary: the 3D companion of 2D_flame_vortex_pair.

    ./mfc.sh run examples/3D_flame_vortex_ring/case.py -n 1 -t pre_process -- --N 32
    python3 examples/3D_flame_vortex_ring/ic3d.py examples/3D_flame_vortex_ring
    ./mfc.sh run examples/3D_flame_vortex_ring/case.py -n 1 -t simulation post_process -- --N 32
    python3 examples/3D_flame_vortex_ring/analyze_fv3d.py examples/3D_flame_vortex_ring

Cube [0, L]^3, L = 3.2 mm: slip wall at x = 0, characteristic subsonic outflow at x = L (bc_x%end = -8, GRCBC to 1 atm), periodic in y
and z. pre_process writes a uniform state and ic3d.py overwrites restart_data/lustre_0.dat with the Cantera flame profile plus the ring
(radius 0.4 mm, Gaussian core 0.15 mm, circulation 0.022 m^2/s, axis along x, centre at x = 0.9 mm), so the case is not auto-registered
as a test. Fixed dt = 5e-9 s on every grid. Production values of the validation study: N = 64 to 256, --tend 1e-4 (--nyz 32 and --gamma 0
in ic3d.py give the flat-flame control).
"""

import argparse
import importlib.util
import json
import os

here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("fv_common", os.path.join(here, "..", "2D_flame_vortex_pair", "fv_common.py"))
fv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fv)

parser = argparse.ArgumentParser(prog="3D_flame_vortex_ring")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC toolchain state.")
parser.add_argument("--N", type=int, default=32, help="Cells per direction (x extent L).")
parser.add_argument("--nyz", type=int, default=0, help="Cells in y and z (0: N, a cube); the y-z extent keeps dx.")
parser.add_argument("--tend", type=float, default=1.0e-4, help="End time [s].")
parser.add_argument("--save-us", type=float, default=10.0, help="Save interval [us].")
parser.add_argument("--dt", type=float, default=5.0e-9, help="Fixed time step [s].")
parser.add_argument("--set", action="append", default=[], metavar="KEY=VAL", help="Override a case parameter.")
args, _ = parser.parse_known_args()

N = args.N
NYZ = args.nyz or N
nsave = round(args.save_us * 1e-6 / args.dt)
g = fv.gas_unburnt()
case = fv.numerics(g, args.mfc)
case.update(
    {
        "m": N - 1,
        "n": NYZ - 1,
        "p": NYZ - 1,
        "dt": args.dt,
        "t_step_start": 0,
        "t_step_stop": nsave * round(args.tend / (nsave * args.dt)),
        "t_step_save": nsave,
        "t_step_print": nsave,
        "bc_y%beg": -1,
        "bc_y%end": -1,
        "bc_z%beg": -1,
        "bc_z%end": -1,
        "num_patches": 1,
        "patch_icpp(1)%geometry": 9,
        "patch_icpp(1)%vel(1)": 0.0,
        "patch_icpp(1)%vel(2)": 0.0,
        "patch_icpp(1)%vel(3)": 0.0,
        "patch_icpp(1)%pres": fv.P0,
        "patch_icpp(1)%alpha_rho(1)": g.density,
        "patch_icpp(1)%alpha(1)": 1.0,
    }
)
for ax, ext in zip("xyz", (fv.L, fv.L * NYZ / N, fv.L * NYZ / N)):
    case.update({f"{ax}_domain%beg": 0.0, f"{ax}_domain%end": ext, f"patch_icpp(1)%{ax}_centroid": 0.5 * ext, f"patch_icpp(1)%length_{ax}": ext})
for i, y in enumerate(g.Y):
    case[f"patch_icpp(1)%Y({i + 1})"] = float(y)
    case[f"chem_wrt_Y({i + 1})"] = "T"
fv.apply_set(case, args.set)

if __name__ == "__main__":
    print(json.dumps(case))
