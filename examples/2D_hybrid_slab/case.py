#!/usr/bin/env python3
# 2D hybrid-rocket fuel slab: a boundary-layer diffusion flame sustained by
# prescribed surface fuel injection. Hot O2/Ar crossflow meets H2 blown in
# through a Dirichlet strip on the bottom wall, approximating steady solid-fuel
# regression without modeling the solid itself. Chemistry ON with diffusion
# ON (unlike the shock-driven detonation case) since this is a diffusion flame.
import argparse
import json
import sys

import cantera as ct

parser = argparse.ArgumentParser(prog="2D_hybrid_slab")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC toolchain state.")
parser.add_argument("--tend", type=float, default=1.2e-3, help="Physical end time [s] (~1 crossflow flow-through Lx/u_ox so the wall-injected flame anchors, not just the startup transient).")
parser.add_argument("--res", type=float, default=1.0, help="Resolution multiplier: m=200*res, n=50*res.")
parser.add_argument("--frames", type=int, default=5, help="Number of saved output frames.")
args = parser.parse_args()

ctfile = "h2o2.yaml"

# --- Domain: x = crossflow direction, y = wall-normal. ---
Lx, Ly = 0.08, 0.02  # 8 cm x 2 cm
Nx = int(200 * args.res)
Ny = int(50 * args.res)
dx, dy = Lx / Nx, Ly / Ny

# Oxidizer crossflow (fills domain, enters at bc_x%beg): hot enough to ignite
# H2 on contact (H2/O2 chain-branching crossover is ~1100 K; 1500 K gives a
# short ignition delay so the flame anchors close to the injection strip).
oxidizer = ct.Solution(ctfile)
oxidizer.TPX = 1500.0, 101325.0, "O2:1,AR:3"
u_ox = 75.0  # crossflow velocity [m/s]

# Fuel blown normal to the wall through a strip on bc_y%beg (models steady
# regression of a solid fuel slab without simulating the solid).
fuel = ct.Solution(ctfile)
fuel.TPX = 600.0, 101325.0, "H2:1"
v_blow = 5.0  # wall-normal blowing velocity [m/s]

print(f"oxidizer: T={oxidizer.T:.1f} K rho={oxidizer.density:.4f} c={oxidizer.sound_speed:.1f} m/s", file=sys.stderr)
print(f"fuel:     T={fuel.T:.1f} K rho={fuel.density:.4f} c={fuel.sound_speed:.1f} m/s", file=sys.stderr)

# Fuel injection strip: x = [2 cm, 6 cm] on the bottom wall.
x_fuel_beg, x_fuel_end = 0.02, 0.06

# Explicit-solver CFL from the fastest signal (sound speed dominates; the
# blowing/crossflow velocities are subsonic).
c_max = max(oxidizer.sound_speed, fuel.sound_speed)
u_max = max(u_ox, v_blow)
cfl = 0.3
dt = cfl * dx / (u_max + c_max)

NT = int(args.tend / dt)
NS = max(1, NT // args.frames)
NP = max(1, NT // 20)

case = {
    "run_time_info": "T",
    # Domain
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
    "t_step_print": NP,
    "parallel_io": "T",
    # Algorithm
    "model_eqns": "5eq",
    "num_fluids": 1,
    "num_patches": 2,
    "mpp_lim": "F",
    "mixture_err": "F",
    "weno_avg": "F",
    "time_stepper": "rk3",
    "weno_order": 5,
    "weno_eps": 1e-16,
    "mapped_weno": "T",
    "mp_weno": "T",
    "riemann_solver": "hllc",
    "wave_speeds": "direct",
    "avg_state": "arithmetic",
    "viscous": "T",
    # BCs: hot oxidizer Dirichlet inflow, extrapolated outlet, no-slip wall with
    # a Dirichlet fuel-injection strip, extrapolated (open) top.
    "bc_x%beg": -17,
    "bc_x%end": -3,
    "bc_y%beg": -16,
    "bc_y%end": -3,
    "num_bc_patches": 1,
    "patch_bc(1)%geometry": 1,
    "patch_bc(1)%type": -17,
    "patch_bc(1)%dir": 2,
    "patch_bc(1)%loc": -1,
    "patch_bc(1)%centroid(1)": (x_fuel_beg + x_fuel_end) / 2,
    "patch_bc(1)%length(1)": x_fuel_end - x_fuel_beg,
    # Chemistry
    "chemistry": "T",
    "chem_params%diffusion": "T",
    "chem_params%reactions": "T",
    "cantera_file": ctfile,
    # Output
    "format": "silo",
    "precision": "double",
    "prim_vars_wrt": "T",
    "chem_wrt_T": "T",
    # Patch 1: hot oxidizer crossflow, whole domain
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": Lx / 2,
    "patch_icpp(1)%y_centroid": Ly / 2,
    "patch_icpp(1)%length_x": Lx,
    "patch_icpp(1)%length_y": Ly,
    "patch_icpp(1)%vel(1)": u_ox,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": oxidizer.P,
    "patch_icpp(1)%alpha(1)": 1.0,
    "patch_icpp(1)%alpha_rho(1)": oxidizer.density,
    # Patch 2: fuel injection strip, first wall-adjacent cell row only (this
    # row is what gets frozen into the Dirichlet ghost-cell buffer).
    "patch_icpp(2)%geometry": 3,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%x_centroid": (x_fuel_beg + x_fuel_end) / 2,
    "patch_icpp(2)%y_centroid": dy / 2,
    "patch_icpp(2)%length_x": x_fuel_end - x_fuel_beg,
    "patch_icpp(2)%length_y": dy,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": v_blow,
    "patch_icpp(2)%pres": fuel.P,
    "patch_icpp(2)%alpha(1)": 1.0,
    "patch_icpp(2)%alpha_rho(1)": fuel.density,
    # Fluid EOS (ideal-gas closure is bypassed by chemistry, but gamma/Re must be set)
    "fluid_pp(1)%gamma": 1.0 / (1.4 - 1.0),
    "fluid_pp(1)%pi_inf": 0.0,
    "fluid_pp(1)%Re(1)": 1.0 / oxidizer.viscosity,  # mu = 1/Re; use the physical O2/Ar viscosity at T0
}

# Species mass fractions per patch + per-species output
for i in range(len(oxidizer.Y)):
    case[f"chem_wrt_Y({i + 1})"] = "T"
    case[f"patch_icpp(1)%Y({i + 1})"] = float(oxidizer.Y[i])
    case[f"patch_icpp(2)%Y({i + 1})"] = float(fuel.Y[i])

if __name__ == "__main__":
    print(json.dumps(case))
