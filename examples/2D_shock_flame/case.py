#!/usr/bin/env python3
# 2D reactive shock-flame interaction (Richtmyer-Meshkov). A planar shock is driven
# through hot burned gas and strikes the flame -- the density interface between the
# hot, light products and the cold, dense fresh 2H2+O2+7Ar. The impulsive
# acceleration of that interface drives the Richtmyer-Meshkov instability: the
# initially-flat front (given a small transverse seed) wrinkles and rolls up, its
# surface area grows, and the burn accelerates -- the mechanism behind flame
# acceleration and deflagration-to-detonation transition. Chemistry stays stable
# through the shocked interface via operator-split reaction sub-stepping.
import argparse
import json
import sys

import cantera as ct

parser = argparse.ArgumentParser(prog="2D_shock_flame")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC toolchain state.")
parser.add_argument("--scale", type=float, default=1.0, help="Grid multiplier.")
parser.add_argument("--shockp", type=float, default=10.0, help="Shock driver over-pressure factor.")
parser.add_argument("--shockvel", type=float, default=900.0, help="Driver (post-shock) velocity [m/s], pushes the shock toward the flame.")
parser.add_argument("--amp", type=float, default=0.1, help="Flame-interface wrinkle amplitude (fraction of Ly) for the RM instability.")
parser.add_argument("--kmode", type=int, default=4, help="Transverse perturbation wavelengths.")
parser.add_argument("--tend", type=float, default=1.0e-4, help="Physical end time [s].")
args = parser.parse_args()

ctfile = "h2o2.yaml"
X = "H2:2,O2:1,AR:7"
T0, P0 = 300.0, 6670.0

fresh = ct.Solution(ctfile)
fresh.TPX = T0, P0, X
burned = ct.Solution(ctfile)  # hot products behind the flame
burned.TPX = T0, P0, X
burned.equilibrate("HP")
driver = ct.Solution(ctfile)  # over-pressured burned gas that launches the shock
driver.TPX = T0, P0, X
driver.equilibrate("HP")
driver.SP = driver.entropy_mass, args.shockp * driver.P
print(f"fresh T={fresh.T:.0f} rho={fresh.density:.4f} | burned T={burned.T:.0f} rho={burned.density:.4f} | driver P={driver.P:.2e}", file=sys.stderr)

Ly = 0.03
Lx = 4.0 * Ly
Ny = int(160 * args.scale)
Nx = int(4 * Ny)
dx = Lx / Nx
x_drv = 0.10 * Lx  # shock driver at the far left
x_flame = 0.40 * Lx  # mean flame-interface position (burned | fresh)
A_flame = args.amp * Ly  # wrinkle amplitude of the flame interface

dt = 0.05 * dx / (1600.0 + driver.sound_speed)
NT = int(args.tend / dt)
NS = max(1, NT // 90)

case = {
    "run_time_info": "T",
    "x_domain%beg": 0.0,
    "x_domain%end": Lx,
    "y_domain%beg": 0.0,
    "y_domain%end": Ly,
    "m": Nx,
    "n": Ny,
    "p": 0,
    "dt": float(dt),
    "t_step_start": 0,
    "t_step_stop": NT,
    "t_step_save": NS,
    "t_step_print": NS,
    "parallel_io": "T",
    "model_eqns": "5eq",
    "num_fluids": 1,
    "num_patches": 3,
    "mpp_lim": "F",
    "mixture_err": "T",
    "weno_avg": "F",
    "time_stepper": "rk3",
    "weno_order": 5,
    "weno_eps": 1e-16,
    "mapped_weno": "T",
    "mp_weno": "T",
    "riemann_solver": "hllc",
    "wave_speeds": "direct",
    "avg_state": "arithmetic",
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "bc_y%beg": -1,
    "bc_y%end": -1,
    "chemistry": "T",
    "chem_params%diffusion": "F",
    "chem_params%reactions": "T",
    "chem_params%reaction_substeps": 10,
    "cantera_file": ctfile,
    "format": "silo",
    "precision": "double",
    "prim_vars_wrt": "T",
    "chem_wrt_T": "T",
    # Patch 1: fresh reactants (cold, dense), whole domain
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": Lx / 2,
    "patch_icpp(1)%y_centroid": Ly / 2,
    "patch_icpp(1)%length_x": Lx,
    "patch_icpp(1)%length_y": Ly,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": fresh.P,
    "patch_icpp(1)%alpha(1)": 1.0,
    "patch_icpp(1)%alpha_rho(1)": fresh.density,
    # Patch 2: hot burned products (light) left of the flame; hcid 280 carves a
    # sinusoidal (burned | fresh) interface so the shock has a wrinkle to amplify.
    "patch_icpp(2)%geometry": 3,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%x_centroid": (x_flame + A_flame) / 2,
    "patch_icpp(2)%y_centroid": Ly / 2,
    "patch_icpp(2)%length_x": x_flame + A_flame,
    "patch_icpp(2)%length_y": Ly,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": burned.P,
    "patch_icpp(2)%alpha(1)": 1.0,
    "patch_icpp(2)%alpha_rho(1)": burned.density,
    "patch_icpp(2)%hcid": 275,
    "patch_icpp(2)%a(2)": x_flame,
    "patch_icpp(2)%a(3)": A_flame,
    "patch_icpp(2)%a(4)": float(args.kmode),
    # Patch 3: over-pressured driver at the far left -> launches a shock toward the flame.
    # It sits inside the burned gas, so it must overwrite both patch 1 and patch 2 cells.
    "patch_icpp(3)%geometry": 3,
    "patch_icpp(3)%alter_patch(1)": "T",
    "patch_icpp(3)%alter_patch(2)": "T",
    "patch_icpp(3)%x_centroid": x_drv / 2,
    "patch_icpp(3)%y_centroid": Ly / 2,
    "patch_icpp(3)%length_x": x_drv,
    "patch_icpp(3)%length_y": Ly,
    "patch_icpp(3)%vel(1)": args.shockvel,
    "patch_icpp(3)%vel(2)": 0.0,
    "patch_icpp(3)%pres": driver.P,
    "patch_icpp(3)%alpha(1)": 1.0,
    "patch_icpp(3)%alpha_rho(1)": driver.density,
    "fluid_pp(1)%gamma": 1.0 / (1.4 - 1.0),
    "fluid_pp(1)%pi_inf": 0.0,
}

for i in range(len(fresh.Y)):
    case[f"chem_wrt_Y({i + 1})"] = "T"
    case[f"patch_icpp(1)%Y({i + 1})"] = float(fresh.Y[i])
    case[f"patch_icpp(2)%Y({i + 1})"] = float(burned.Y[i])
    case[f"patch_icpp(3)%Y({i + 1})"] = float(driver.Y[i])

if __name__ == "__main__":
    print(json.dumps(case))
