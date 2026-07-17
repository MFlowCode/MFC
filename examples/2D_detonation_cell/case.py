#!/usr/bin/env python3
# 2D cellular detonation in argon-diluted H2/O2 (2H2 + O2 + 7Ar), driven by the
# automated stiff-chemistry pipeline:
#
#   * cfl_const_dt  -- the flow timestep is computed once from the hydrodynamic CFL
#                      condition (dt = cfl_target * min dx/(|u|+c)); no hand-tuned dt.
#                      Because the reaction source is operator-split out of the flow
#                      RHS (alpha-QSS, below), the CFL sees a clean sound speed and
#                      takes the full *hydro*-stable step.
#   * adap_substeps -- the operator-split alpha-QSS reaction integrator chooses one
#                      sub-step count per flow step, per rank, sized from the stiffest
#                      cell that rank owns: a rank of only burned/unburned gas sits at
#                      the floor (reaction_substeps), and any rank spanning the
#                      detonation front ramps to the ceiling (reaction_substeps_max).
#
# Integrating the reaction source explicitly in the RK3 RHS (reaction_substeps=0)
# caps this detonation near CFL ~0.06; alpha-QSS runs the *same* physics at CFL ~0.5
# -- an ~8x larger timestep -- because stiffness then sets a small chemistry *sub-step*
# rather than a small *flow step*. Heavy Ar dilution -> regular cells. Diffusion OFF
# (shock-driven). Reactions ON.
import argparse
import json
import sys

import cantera as ct

parser = argparse.ArgumentParser(prog="2D_detonation_cell")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC toolchain state.")
parser.add_argument("--scale", type=float, default=1.0, help="Grid multiplier (calibration knob).")
parser.add_argument("--ndim", type=int, default=2, choices=(2, 3), help="Spatial dimensions (2 or 3).")
parser.add_argument("--tend", type=float, default=200e-6, help="Physical end time [s].")
parser.add_argument("--cfl-target", type=float, default=0.5, help="Target hydrodynamic CFL; sets the (constant) flow timestep.")
parser.add_argument("--substeps", type=int, default=2, help="alpha-QSS sub-step floor -- the count a rank uses when none of its cells are stiff.")
parser.add_argument("--substeps-max", type=int, default=16, help="alpha-QSS sub-step ceiling (used across the detonation front).")
parser.add_argument(
    "--overdrive",
    type=float,
    default=1.5,
    help="Isentropic compression of the UV-equilibrium driver (1 = plain; >1 = overdriven, thermo-consistent). Kept mild so the driver/fresh contact stays stable.",
)
parser.add_argument(
    "--drivervel",
    type=float,
    default=1200.0,
    help="Driver forward (piston) velocity [m/s]; a velocity-driven shock heats fresh gas above the ~1100 K H2/O2 crossover so it ignites and couples into a CJ detonation.",
)
parser.add_argument("--seed", dest="seed", default=True, action="store_true", help="Add a coherent transverse-velocity perturbation to seed regular detonation cells.")
parser.add_argument("--no-seed", dest="seed", action="store_false")
parser.add_argument("--riemann", default="hllc", choices=("hllc", "hll"), help="Riemann solver (hll is more dissipative / carbuncle-free).")
parser.add_argument("--seedamp", type=float, default=40.0, help="Cell-seed transverse-velocity amplitude [m/s] (lower it in 3D, where seeding y and z together is stronger).")
args = parser.parse_args()

ctfile = "h2o2.yaml"
X = "H2:2,O2:1,AR:7"  # Ar-diluted -> regular cells
T0, P0 = 298.0, 6670.0  # low p lengthens the reaction zone

# Fresh premixed reactants (fills the domain).
fresh = ct.Solution(ctfile)
fresh.TPX = T0, P0, X

# Detonation driver: constant-volume explosion products -> a thermodynamically
# self-consistent hot, high-pressure burned state (T, P, rho all from Cantera).
# --overdrive (>1) isentropically compresses the products for extra drive, scaling
# P and rho together so T stays in the thermo range.
driver = ct.Solution(ctfile)
driver.TPX = T0, P0, X
driver.equilibrate("UV")
if args.overdrive > 1.0:
    driver.SP = driver.entropy_mass, args.overdrive * driver.P
P_drive = driver.P
print(f"driver init: T={driver.T:.1f} K  P={driver.P:.1f} Pa  rho={driver.density:.5f} kg/m^3", file=sys.stderr)

# --- Domain: x = propagation, y (and z in 3D) = transverse (periodic). ---
is_3d = args.ndim == 3
Ly = 0.03  # channel height [m]
Lx = 4.0 * Ly  # run-up length
Lz = Ly
Ny = int(200 * args.scale)
Nx = int(4 * Ny)
Nz = Ny if is_3d else 0
geom = 9 if is_3d else 3  # 3D box vs 2D rectangle

# Hot high-pressure driver slug fills the left 25%; its expansion (plus the piston-velocity IC)
# launches a finite-energy blast that locks into a self-sustained CJ detonation.
x_driver = 0.25 * Lx

case = {
    "run_time_info": "T",
    # Domain
    "x_domain%beg": 0.0,
    "x_domain%end": Lx,
    "y_domain%beg": 0.0,
    "y_domain%end": Ly,
    "m": Nx,
    "n": Ny,
    "p": Nz,
    # Automated flow timestep: computed once from the CFL target (no fixed dt)
    "cfl_const_dt": "T",
    "cfl_target": args.cfl_target,
    "n_start": 0,
    "t_save": args.tend / 100.0,
    "t_stop": args.tend,
    "parallel_io": "F",
    # Algorithm
    "model_eqns": "5eq",
    "num_fluids": 1,
    "num_patches": 2,
    "mpp_lim": "F",
    "mixture_err": "F",  # WENO5 + mp_weno bound the reconstruction, so unphysical p/T never form
    "weno_avg": "F",
    "time_stepper": "rk3",
    "weno_order": 5,
    "weno_eps": 1e-16,
    "mapped_weno": "T",
    "mp_weno": "T",
    "riemann_solver": args.riemann,
    "wave_speeds": "direct",
    "avg_state": "arithmetic",
    # BCs: x extrapolation (open), y periodic (clean regular cells)
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "bc_y%beg": -1,
    "bc_y%end": -1,
    # Chemistry: operator-split alpha-QSS reaction integrator with per-rank adaptive sub-stepping
    "chemistry": "T",
    "chem_params%diffusion": "F",
    "chem_params%reactions": "T",
    "chem_params%reaction_substeps": args.substeps,
    "chem_params%reaction_substeps_max": args.substeps_max,
    "chem_params%adap_substeps": "T",
    "cantera_file": ctfile,
    # Output
    "format": "silo",
    "precision": "double",
    "prim_vars_wrt": "T",
    "chem_wrt_T": "T",
    # Patch 1: fresh reactants, whole domain
    "patch_icpp(1)%geometry": geom,
    "patch_icpp(1)%x_centroid": Lx / 2,
    "patch_icpp(1)%y_centroid": Ly / 2,
    "patch_icpp(1)%length_x": Lx,
    "patch_icpp(1)%length_y": Ly,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": fresh.P,
    "patch_icpp(1)%alpha(1)": 1.0,
    "patch_icpp(1)%alpha_rho(1)": fresh.density,
    # Patch 2: overdriven driver (hot products), left region, full cross-section
    "patch_icpp(2)%geometry": geom,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%x_centroid": x_driver / 2,
    "patch_icpp(2)%y_centroid": Ly / 2,
    "patch_icpp(2)%length_x": x_driver,
    "patch_icpp(2)%length_y": Ly,
    "patch_icpp(2)%vel(1)": args.drivervel,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": P_drive,
    "patch_icpp(2)%alpha(1)": 1.0,
    "patch_icpp(2)%alpha_rho(1)": driver.density,
    # Fluid EOS (ideal-gas closure is bypassed by chemistry, but gamma must be set)
    "fluid_pp(1)%gamma": 1.0 / (1.4 - 1.0),
    "fluid_pp(1)%pi_inf": 0.0,
}

# 3D: add the periodic z direction and give both patches full z-extent.
if is_3d:
    case.update(
        {
            "z_domain%beg": 0.0,
            "z_domain%end": Lz,
            "bc_z%beg": -1,
            "bc_z%end": -1,
            "patch_icpp(1)%z_centroid": Lz / 2,
            "patch_icpp(1)%length_z": Lz,
            "patch_icpp(1)%vel(3)": 0.0,
            "patch_icpp(2)%z_centroid": Lz / 2,
            "patch_icpp(2)%length_z": Lz,
            "patch_icpp(2)%vel(3)": 0.0,
        }
    )

# Seed transverse cells with a COHERENT sinusoidal transverse-velocity perturbation on the
# driver: this imposes one clean cellular wavelength (a per-cell random perturbation instead
# injects grid-scale white noise that advects into y-striations, not physical cells).
kmode = 3  # transverse wavelengths across the channel
amp = args.seedamp  # perturbation amplitude [m/s], a few % of the piston velocity
if args.seed:
    case["patch_icpp(2)%vel(2)"] = f"{amp}*sin(2*pi*{kmode}*y/{Ly})"
    if is_3d:
        case["patch_icpp(2)%vel(3)"] = f"{amp}*sin(2*pi*{kmode}*z/{Lz})"

# Species mass fractions per patch + per-species output
for i in range(len(fresh.Y)):
    case[f"chem_wrt_Y({i + 1})"] = "T"
    case[f"patch_icpp(1)%Y({i + 1})"] = float(fresh.Y[i])
    case[f"patch_icpp(2)%Y({i + 1})"] = float(driver.Y[i])

if __name__ == "__main__":
    print(json.dumps(case))
