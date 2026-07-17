#!/usr/bin/env python3
"""3D spherical PETN ZND detonation, octant symmetry, grid-convergence study.

One octant of a spherical PETN charge (Kuhl et al., UCRL-PROC-225822 (2006)
Appendix JWL fit; rho0 = 1.0 g/cc, Kuhl's own pressed-charge reference density
-- see toolchain/mfc/jwl_products.py) sits unreacted at ambient pressure in
free air. Reflective symmetry planes at x/y/z = 0 (bc_*%beg = -2) recover the
full spherical geometry at 1/8 the cost -- the same trick already used in
examples/2D_sedov_taylor_blast and benchmarks/3D_jwl_spherical_tnt_free_air_validation,
here with MFC's genuine 3D machinery rather than an unimplemented native
spherical (1/r^2) source term.

Unlike the instantaneous-burst TNT benchmark, this case uses a resolved
self-propagating JWL++ detonation (jwl_reactive, Souers 2000) with the Garno
et al. (2020) Eq. 17 reactant energy offset (fluid_pp(1)%jwl_delta_e), so the
detonation carries genuine ZND structure (von Neumann spike -> reaction zone
-> CJ plane -> Taylor rarefaction) rather than a kinematic or instantaneous
release. A pressure-only hot spot was found (in 1D calibration trials) to be
sub-critical for this material's CJ pressure -- it fizzles rather than
transitioning to a self-sustained wave, consistent with the documented
spherical-divergence quenching risk (senior-cfd-physics.md P4/P6). The
booster patch is therefore seeded already-reacted (rxn_val = 1) for robust
direct initiation.

The self-sustained CJ state of this material (Hugoniot/Rayleigh tangency,
EOS only) is D_CJ ~ 5384 m/s, P_CJ ~ 7.35 GPa; Kuhl reports only T_CJ ~ 4600
K, not P_CJ/D_CJ. The 1D calibration run is overdriven by its 30 GPa booster
plus closed rear boundary (it propagates at ~6245 m/s / ~16.7 GPa, exactly
the overdriven lambda=1 Hugoniot state at that speed), so jwl_G here sets the
reaction-zone width but the front speed is piston-set, not CJ. See
examples/3D_jwl_znd_spherical_octant/VALIDATION.md.

Grid convergence (P8 protocol: 3 grids x factor 2 in dx):
    ./mfc.sh run examples/3D_jwl_znd_spherical_octant/case.py -n 8 -- --grid 39   # dx = 1.00 mm
    ./mfc.sh run examples/3D_jwl_znd_spherical_octant/case.py -n 8 -- --grid 79   # dx = 0.50 mm
    ./mfc.sh run examples/3D_jwl_znd_spherical_octant/case.py -n 8 -- --grid 159  # dx = 0.25 mm
Each --grid value writes to its own examples/3D_jwl_znd_spherical_octant/grid_<N>/
subdirectory (see run_convergence.sh); analyze_convergence.py reduces the probe
histories to a front-arrival-time fit at each resolution and reports the
observed order of convergence on the front speed.
"""

import argparse
import json

from mfc.jwl_products import AIR, PETN_KUHL, ambient_fluid, jwl_fluid, znd_delta_e

parser = argparse.ArgumentParser(description="3D JWL PETN spherical ZND octant case")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC toolchain state")
parser.add_argument("--grid", type=int, default=79, help="MFC m/n/p index (dx = L/(grid+1))")
parser.add_argument("--t-stop", type=float, default=6.0e-6, help="final simulation time in seconds")
args = parser.parse_args()

eps = 1.0e-8

# --- Geometry: octant of a spherical PETN charge in free air ---
L = 0.04  # octant edge length [m]; full domain if mirrored: [-L, L]^3
charge_radius = 0.015
booster_radius = 0.0025
p_amb = 101325.0

dx = L / (args.grid + 1)
probe_yz = 0.5 * dx


def nearest_cell_center(x):
    i = int(round(x / dx - 0.5))
    i = max(0, min(args.grid, i))
    return (i + 0.5) * dx


# Radial probes spanning charge interior, charge edge, and breakout/air region.
target_probe_r = (0.3 * charge_radius, 0.6 * charge_radius, charge_radius, 1.5 * charge_radius, 2.0 * charge_radius)
probe_r = tuple(nearest_cell_center(r) for r in target_probe_r)

# --- PETN JWL products (Kuhl et al. 2006 fit) + reactant energy offset ---
jwl_delta_e = znd_delta_e(PETN_KUHL, p_amb)

# --- JWL++ rate: dl/dt = jwl_G * p^jwl_b_exp * (1 - lambda) ---
# Calibrated in 1D (see module docstring); this G gives a ~1.25 mm reaction
# zone, i.e. O(1) cells at the coarsest grid level here and increasingly
# resolved at finer levels -- the grid-convergence study below is exactly
# what demonstrates the front-speed/structure trend as that zone becomes
# resolved.
jwl_G = 1.0e-13
jwl_b_exp = 2.0

params = {
    "run_time_info": "T",
    "x_domain%beg": 0.0,
    "x_domain%end": L,
    "y_domain%beg": 0.0,
    "y_domain%end": L,
    "z_domain%beg": 0.0,
    "z_domain%end": L,
    "m": args.grid,
    "n": args.grid,
    "p": args.grid,
    "t_step_start": 0,
    "t_step_stop": 200000,
    "t_step_save": 200000,
    "cfl_adap_dt": "T",
    "cfl_target": 0.3,
    "n_start": 0,
    "t_stop": args.t_stop,
    "t_save": args.t_stop / 20.0,
    "num_patches": 3,
    "model_eqns": 2,
    "num_fluids": 2,
    "mpp_lim": "T",
    "mixture_err": "T",
    "time_stepper": 3,
    "weno_order": 3,
    "weno_eps": 1.0e-16,
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "F",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    "bc_x%beg": -2,
    "bc_x%end": -3,
    "bc_y%beg": -2,
    "bc_y%end": -3,
    "bc_z%beg": -2,
    "bc_z%end": -3,
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "parallel_io": "F",
    "probe_wrt": "T",
    "fd_order": 1,
    "num_probes": len(probe_r),
    # Self-propagating JWL++ reactive burn (lambda = 0 everywhere initially).
    "jwl_reactive": "T",
    "jwl_G": jwl_G,
    "jwl_b_exp": jwl_b_exp,
    # Patch 1: ambient air filling the octant.
    "patch_icpp(1)%geometry": 9,
    "patch_icpp(1)%x_centroid": L / 2.0,
    "patch_icpp(1)%y_centroid": L / 2.0,
    "patch_icpp(1)%z_centroid": L / 2.0,
    "patch_icpp(1)%length_x": L,
    "patch_icpp(1)%length_y": L,
    "patch_icpp(1)%length_z": L,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%vel(3)": 0.0,
    "patch_icpp(1)%pres": p_amb,
    "patch_icpp(1)%alpha_rho(1)": eps * PETN_KUHL["rho0"],
    "patch_icpp(1)%alpha_rho(2)": (1.0 - eps) * AIR["rho0"],
    "patch_icpp(1)%alpha(1)": eps,
    "patch_icpp(1)%alpha(2)": 1.0 - eps,
    # Patch 2: PETN charge sphere at the symmetry corner, unreacted, ambient pressure.
    "patch_icpp(2)%geometry": 8,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%x_centroid": 0.0,
    "patch_icpp(2)%y_centroid": 0.0,
    "patch_icpp(2)%z_centroid": 0.0,
    "patch_icpp(2)%radius": charge_radius,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%vel(3)": 0.0,
    "patch_icpp(2)%pres": p_amb,
    "patch_icpp(2)%alpha_rho(1)": (1.0 - eps) * PETN_KUHL["rho0"],
    "patch_icpp(2)%alpha_rho(2)": eps * AIR["rho0"],
    "patch_icpp(2)%alpha(1)": 1.0 - eps,
    "patch_icpp(2)%alpha(2)": eps,
    # Patch 3: central booster, seeded already-reacted (rxn_val = 1) for robust
    # direct initiation -- a pressure-only hot spot is sub-critical here (see
    # module docstring).
    "patch_icpp(3)%geometry": 8,
    "patch_icpp(3)%alter_patch(1)": "T",
    "patch_icpp(3)%alter_patch(2)": "T",
    "patch_icpp(3)%x_centroid": 0.0,
    "patch_icpp(3)%y_centroid": 0.0,
    "patch_icpp(3)%z_centroid": 0.0,
    "patch_icpp(3)%radius": booster_radius,
    "patch_icpp(3)%vel(1)": 0.0,
    "patch_icpp(3)%vel(2)": 0.0,
    "patch_icpp(3)%vel(3)": 0.0,
    "patch_icpp(3)%pres": 25.0e9,
    "patch_icpp(3)%alpha_rho(1)": (1.0 - eps) * PETN_KUHL["rho0"],
    "patch_icpp(3)%alpha_rho(2)": eps * AIR["rho0"],
    "patch_icpp(3)%alpha(1)": 1.0 - eps,
    "patch_icpp(3)%alpha(2)": eps,
    "patch_icpp(3)%rxn_val": 1.0,
    # Fluid 1: PETN JWL products (+ ambient references); Fluid 2: air.
    **jwl_fluid(1, PETN_KUHL, AIR, delta_e=jwl_delta_e),
    **ambient_fluid(2, AIR),
}

for i, r in enumerate(probe_r, start=1):
    params[f"probe({i})%x"] = r
    params[f"probe({i})%y"] = probe_yz
    params[f"probe({i})%z"] = probe_yz

print(json.dumps(params))
