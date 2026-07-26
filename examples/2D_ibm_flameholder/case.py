#!/usr/bin/env python3
# 2D bluff-body flame holder: a solid cylinder (immersed boundary) injects H2 off
# its surface into a hot O2/Ar crossflow. The freestream carries no fuel, so it
# cannot auto-ignite; the flame anchors only where the injected H2 meets the hot
# oxidizer -- in the recirculating wake behind the cylinder. This is the classic
# bluff-body-stabilized (afterburner/ramjet) flame, and it exercises the IBM
# surface-injection knobs (patch_ib%v_blow, %inj_species) together with the
# IBM+chemistry ghost-state fix.
import argparse
import json
import sys

import cantera as ct

parser = argparse.ArgumentParser(prog="2D_ibm_flameholder")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC toolchain state.")
parser.add_argument("--tend", type=float, default=3.0e-4, help="Physical end time [s].")
parser.add_argument("--res", type=float, default=1.0, help="Resolution multiplier: m=300*res, n=112*res.")
parser.add_argument("--frames", type=int, default=80, help="Number of saved output frames.")
args = parser.parse_args()

ctfile = "h2o2.yaml"

# Hot oxidizer crossflow (fills the domain, enters at bc_x%beg). No fuel -> the
# freestream cannot burn on its own; T above the ~1100 K H2/O2 crossover so the
# injected fuel ignites promptly on contact.
ox = ct.Solution(ctfile)
ox.TPX = 1200.0, 101325.0, "O2:1,AR:3"
u_ox = 60.0  # crossflow velocity [m/s]
c0 = ox.sound_speed
mu0 = ox.viscosity

# Domain: 8 cm x 3 cm channel, cylinder near the inlet.
Lx, Ly = 0.08, 0.03
m = int(300 * args.res)
n = int(112 * args.res)
dx = Lx / m

r_cyl = 0.004
x_cyl = 0.02
y_cyl = Ly / 2.0
v_blow = 12.0  # wall-normal H2 injection speed off the cylinder [m/s] (gentle: strong
#   injection drives an unstable, super-heated ignition transient that NaNs)

print(f"oxidizer: T={ox.T:.1f} K rho={ox.density:.4f} c={c0:.1f} m/s", file=sys.stderr)

# Acoustic CFL; small for explicit-chemistry stability (reacting hot spots need headroom).
cfl = 0.05
dt = cfl * dx / (u_ox + c0)
NT = int(args.tend / dt)
NS = max(1, NT // args.frames)

case = {
    "run_time_info": "T",
    "x_domain%beg": 0.0,
    "x_domain%end": Lx,
    "y_domain%beg": 0.0,
    "y_domain%end": Ly,
    "m": m,
    "n": n,
    "p": 0,
    "cyl_coord": "F",
    "dt": float(dt),
    "t_step_start": 0,
    "t_step_stop": NT,
    "t_step_save": NS,
    "t_step_print": max(1, NT // 20),
    "parallel_io": "T",
    # Algorithm
    "model_eqns": 2,
    "num_fluids": 1,
    "num_patches": 1,
    "mpp_lim": "F",
    "mixture_err": "T",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1e-16,
    "mapped_weno": "T",
    "mp_weno": "T",
    "weno_avg": "T",
    "weno_Re_flux": "T",
    "null_weights": "F",
    "riemann_solver": "hllc",
    "wave_speeds": "direct",
    "avg_state": "arithmetic",
    "fd_order": 2,
    "viscous": "T",
    # BCs: characteristic inflow, extrapolated outflow, slip channel walls
    "bc_x%beg": -7,
    "bc_x%end": -3,
    "bc_y%beg": -15,
    "bc_y%end": -15,
    # Chemistry + diffusion + reactions ON (fuel/oxidizer mixing and combustion)
    "chemistry": "T",
    "chem_params%diffusion": "T",
    "chem_params%reactions": "T",
    "cantera_file": ctfile,
    "chem_wrt_T": "T",
    # Immersed boundary: solid cylinder that injects H2 off its surface
    "ib": "T",
    "num_ibs": 1,
    "patch_ib(1)%geometry": 2,
    "patch_ib(1)%x_centroid": x_cyl,
    "patch_ib(1)%y_centroid": y_cyl,
    "patch_ib(1)%radius": r_cyl,
    "patch_ib(1)%slip": "F",
    "patch_ib(1)%v_blow": v_blow,
    "patch_ib(1)%inj_species": 1,  # inject pure H2
    # Output
    "format": "silo",
    "precision": "double",
    "prim_vars_wrt": "T",
    "ib_state_wrt": "T",
    # Patch: hot O2/Ar oxidizer fills the domain
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": Lx / 2,
    "patch_icpp(1)%y_centroid": Ly / 2,
    "patch_icpp(1)%length_x": Lx,
    "patch_icpp(1)%length_y": Ly,
    "patch_icpp(1)%vel(1)": u_ox,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": ox.P,
    "patch_icpp(1)%alpha_rho(1)": ox.density,
    "patch_icpp(1)%alpha(1)": 1.0,
    # Fluid EOS (calorically perfect closure; chemistry supplies the real thermo)
    "fluid_pp(1)%gamma": 1.0 / (1.4 - 1.0),
    "fluid_pp(1)%pi_inf": 0.0,
    "fluid_pp(1)%Re(1)": 1.0 / mu0,
}

for i in range(len(ox.Y)):
    case[f"patch_icpp(1)%Y({i + 1})"] = float(ox.Y[i])
    case[f"chem_wrt_Y({i + 1})"] = "T"

if __name__ == "__main__":
    print(json.dumps(case))
