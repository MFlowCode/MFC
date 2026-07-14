#!/usr/bin/env python3
# 2D solid-propellant "burning grain": a solid cylinder (immersed boundary) that
# INJECTS fuel (pure H2) off its surface into a hot O2/Ar oxidizer chamber. The
# injected H2 auto-ignites against the hot oxidizer and anchors a diffusion flame
# on the injecting surface -- the essential solid-rocket-motor picture (surface
# fuel injection + combustion). Enabled by two new IB knobs:
#   patch_ib(i)%v_blow      -- wall-normal surface blowing speed [m/s]
#   patch_ib(i)%inj_species -- injected species index (1 = H2 in this mechanism)
# Both default off, so non-burning IBM is unchanged. Builds on the IBM+chemistry
# ghost-state fix (m_ibm.fpp): the injecting surface sets a thermodynamically
# consistent reacting ghost state.
#
# --burn_exp n makes the injection pressure-coupled (Vieille's law: v_blow scales
# with the local surface pressure as (p/p0)^n). In this closed chamber that gives
# internal-ballistics feedback -- combustion raises the pressure, which raises the
# burn rate, which raises the pressure -- an accelerating self-pressurization.
import argparse
import json

import cantera as ct

parser = argparse.ArgumentParser(prog="2D_ibm_burning_grain")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC toolchain state.")
parser.add_argument("--burn_exp", type=float, default=0.0, help="Pressure exponent n in Vieille's law v_blow*(p/p0)^n; 0 = constant injection.")
parser.add_argument("--tend", type=float, default=1.5e-4, help="Physical end time [s].")
args = parser.parse_args()

ctfile = "h2o2.yaml"
# Hot oxidizer chamber (no fuel): injected H2 meets hot O2 and ignites.
X = "O2:1,AR:3"
T0, P0 = 1200.0, 101325.0  # above the H2/O2 crossover -> prompt surface ignition

ox = ct.Solution(ctfile)
ox.TPX = T0, P0, X
rho0 = ox.density
mu0 = ox.viscosity
c0 = ox.sound_speed

# Chamber: 4 cm x 3 cm, solid fuel cylinder near the center.
Lx, Ly = 0.04, 0.03
m, n = 160, 120
dx = Lx / m

r_cyl = 0.004
x_cyl = Lx / 2.0
y_cyl = Ly / 2.0

v_blow = 20.0  # surface fuel-injection (blowing) speed [m/s]

# dt from the acoustic CFL; small for explicit chemistry stability.
cfl = 0.05
dt = cfl * dx / (v_blow + c0)
tend = args.tend
NT = int(tend / dt)
NS = max(1, NT // 60)

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
    "t_step_print": NS,
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
    # Closed chamber: reflective/slip walls all around so it pressurizes as it burns.
    "bc_x%beg": -2,
    "bc_x%end": -2,
    "bc_y%beg": -2,
    "bc_y%end": -2,
    # Chemistry + diffusion + reactions ON (fuel/oxidizer mixing and combustion)
    "chemistry": "T",
    "chem_params%diffusion": "T",
    "chem_params%reactions": "T",
    "cantera_file": ctfile,
    "chem_wrt_T": "T",
    # Immersed boundary: solid fuel cylinder that injects H2 off its surface
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
    # Patch: hot O2/Ar oxidizer fills the chamber
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": Lx / 2,
    "patch_icpp(1)%y_centroid": Ly / 2,
    "patch_icpp(1)%length_x": Lx,
    "patch_icpp(1)%length_y": Ly,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": P0,
    "patch_icpp(1)%alpha_rho(1)": rho0,
    "patch_icpp(1)%alpha(1)": 1.0,
    # Fluid EOS (calorically perfect closure; chemistry supplies the real thermo)
    "fluid_pp(1)%gamma": 1.0 / (1.4 - 1.0),
    "fluid_pp(1)%pi_inf": 0.0,
    "fluid_pp(1)%Re(1)": 1.0 / mu0,
}

# Pressure-coupled burn rate (internal-ballistics demo): v_blow -> v_blow*(p/P0)^n.
if args.burn_exp > 0.0:
    case["patch_ib(1)%burn_rate_pref"] = P0
    case["patch_ib(1)%burn_rate_exp"] = args.burn_exp

for i in range(len(ox.Y)):
    case[f"patch_icpp(1)%Y({i + 1})"] = float(ox.Y[i])
    case[f"chem_wrt_Y({i + 1})"] = "T"

if __name__ == "__main__":
    print(json.dumps(case))
