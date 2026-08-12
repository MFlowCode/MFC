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
parser.add_argument("--res", type=float, default=1.0, help="Resolution multiplier: 2D m=160*res, n=120*res; 3D m=120*res, n=p=90*res.")
parser.add_argument("--ndim", type=int, default=2, choices=(2, 3), help="Spatial dimensions: 2 = cylindrical grain (circle IB), 3 = spherical grain (sphere IB).")
args = parser.parse_args()
is_3d = args.ndim == 3

ctfile = "h2o2.yaml"
# Hot oxidizer chamber (no fuel): injected H2 meets hot O2 and ignites.
X = "O2:1,AR:3"
T0, P0 = 1200.0, 101325.0  # above the H2/O2 crossover -> prompt surface ignition

ox = ct.Solution(ctfile)
ox.TPX = T0, P0, X
rho0 = ox.density
mu0 = ox.viscosity
c0 = ox.sound_speed

# Chamber: 4 cm x 3 cm (x 3 cm in 3D), solid fuel grain at the center.
Lx, Ly, Lz = 0.04, 0.03, 0.03
r_cyl = 0.004
if is_3d:
    # Octant symmetry: a centered sphere in a slip-walled cubic chamber is symmetric about
    # all three midplanes, so simulate ONE eighth -- the sphere sits at the inner corner and
    # the three inner faces are slip (= symmetry) walls, identical to the closed-chamber walls.
    # This puts 2x the resolution per dimension on the sphere/flame for the same cost; the
    # renderer mirrors the octant back to the full sphere.
    Dx, Dy, Dz = Lx / 2, Ly / 2, Lz / 2
    m, n, p = int(160 * args.res), int(120 * args.res), int(120 * args.res)
    x_cyl, y_cyl, z_cyl = 0.0, 0.0, 0.0  # sphere centroid on the symmetry corner
else:
    Dx, Dy, Dz = Lx, Ly, Lz
    m, n, p = int(160 * args.res), int(120 * args.res), 0
    x_cyl, y_cyl, z_cyl = Lx / 2, Ly / 2, Lz / 2
dx = Dx / m
ib_geom = 8 if is_3d else 2  # sphere (3D) / circle (2D)
patch_geom = 9 if is_3d else 3  # box (3D) / rectangle (2D)

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
    "x_domain%end": Dx,
    "y_domain%beg": 0.0,
    "y_domain%end": Dy,
    "m": m,
    "n": n,
    "p": p,
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
    "patch_ib(1)%geometry": ib_geom,
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
    # Patch: hot O2/Ar oxidizer fills the (octant) chamber
    "patch_icpp(1)%geometry": patch_geom,
    "patch_icpp(1)%x_centroid": Dx / 2,
    "patch_icpp(1)%y_centroid": Dy / 2,
    "patch_icpp(1)%length_x": Dx,
    "patch_icpp(1)%length_y": Dy,
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

# 3D: add the z direction, close the z walls, and give the sphere / oxidizer box full z-extent.
if is_3d:
    case.update(
        {
            "z_domain%beg": 0.0,
            "z_domain%end": Dz,
            "bc_z%beg": -2,
            "bc_z%end": -2,
            "patch_ib(1)%z_centroid": z_cyl,
            "patch_icpp(1)%z_centroid": Dz / 2,
            "patch_icpp(1)%length_z": Dz,
            "patch_icpp(1)%vel(3)": 0.0,
        }
    )

# Pressure-coupled burn rate (internal-ballistics demo): v_blow -> v_blow*(p/P0)^n.
if args.burn_exp > 0.0:
    case["patch_ib(1)%burn_rate_pref"] = P0
    case["patch_ib(1)%burn_rate_exp"] = args.burn_exp

for i in range(len(ox.Y)):
    case[f"patch_icpp(1)%Y({i + 1})"] = float(ox.Y[i])
    case[f"chem_wrt_Y({i + 1})"] = "T"

if __name__ == "__main__":
    print(json.dumps(case))
