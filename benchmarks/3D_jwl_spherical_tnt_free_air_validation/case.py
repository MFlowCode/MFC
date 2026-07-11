#!/usr/bin/env python3
"""3D JWL spherical TNT free-air benchmark.

This case initializes one octant of a spherical TNT detonation-products region
in ambient air. The symmetry planes at x/y/z = 0 recover the full spherical
free-air blast while using a doubled default grid for sharper local visualization.

Two initiation modes (--mode):
  burst (default): the charge starts already-detonated, at the JWL reference
    (V=1) pressure -- the original baseline this benchmark validates.
  prog_burn: the charge starts unreacted at ambient pressure and is ignited by
    a kinematic prog_burn front at the TNT CJ speed, with jwl_afterburn mixing
    layered on top -- the same source combination validated in
    examples/2D_jwl_prog_burn_afterburn_mixing/case.py. This tests whether a
    real detonation front (rather than an instantaneous burst) changes the
    near-field match against Kinney-Graham, which the burst-IC baseline's own
    near-field error (see validate_kingery.py) has been suspected to reflect.

Run:
    ./mfc.sh run benchmarks/3D_jwl_spherical_tnt_free_air_validation/case.py -n 4
    ./mfc.sh run benchmarks/3D_jwl_spherical_tnt_free_air_validation/case.py -n 4 -- --mode prog_burn
"""

import argparse
import json
import math

parser = argparse.ArgumentParser(description="3D JWL spherical TNT free-air benchmark")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC toolchain state")
parser.add_argument(
    "--grid",
    type=int,
    default=63,
    help="MFC m/n/p index. The default gives 64^3 cells in the octant.",
)
parser.add_argument("--t-stop", type=float, default=2.0e-4, help="final simulation time in seconds")
parser.add_argument(
    "--mode",
    choices=["burst", "prog_burn"],
    default="burst",
    help="burst: direct-pressure detonation-products IC (default, existing baseline). "
    "prog_burn: ambient-pressure charge ignited by a kinematic prog_burn front plus "
    "jwl_afterburn mixing (see module docstring).",
)
args = parser.parse_args()

eps = 1.0e-6

# Geometry and ambient air.
L = 0.5
charge_radius = 0.05
p_amb = 101325.0
rho_air = 1.225
gamma_air = 1.4

# TNT/JWL constants used by the existing MFC JWL examples and benchmarks.
rho_tnt = 1630.0
A = 3.712e11
B = 3.231e9
R1 = 4.15
R2 = 0.95
omega = 0.30
E0 = 1.0089e10

# Initialize products with the JWL reference energy at V = rho0/rho = 1.
V = 1.0
p_cold = A * (1.0 - omega / (R1 * V)) * math.exp(-R1 * V)
p_cold += B * (1.0 - omega / (R2 * V)) * math.exp(-R2 * V)
p_charge = p_cold + omega * E0

dx = L / (args.grid + 1)
probe_yz = 0.5 * dx
target_probe_x = (0.15, 0.25, 0.35, 0.45)

# prog_burn reaction-zone band width: ~3 cells, matching the ratio validated in
# examples/2D_jwl_prog_burn_afterburn_mixing/case.py. This case uses cfl_adap_dt,
# so m_checker.fpp's pb_D_cj*dt <= pb_width startup check does not apply (it only
# fires for fixed-dt runs) -- sizing the band a few cells wide keeps the front from
# outrunning its own deposition band under the adaptive dt in practice.
pb_D_cj = 6930.0  # TNT CJ velocity (m/s)
pb_width = 3.0 * dx


def nearest_cell_center(x):
    i = int(round(x / dx - 0.5))
    i = max(0, min(args.grid, i))
    return (i + 0.5) * dx


probe_x = tuple(nearest_cell_center(x) for x in target_probe_x)

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
    "t_step_stop": 100000,
    "t_step_save": 100000,
    "cfl_adap_dt": "T",
    "cfl_target": 0.3,
    "n_start": 0,
    "t_stop": args.t_stop,
    "t_save": 1.0e-5,
    "num_patches": 2,
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
    "rho_wrt": "T",
    "pres_wrt": "T",
    "parallel_io": "F",
    "probe_wrt": "T",
    "fd_order": 1,
    "num_probes": len(probe_x),
    # Ambient air filling the full octant.
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
    "patch_icpp(1)%alpha_rho(1)": eps * rho_tnt,
    "patch_icpp(1)%alpha_rho(2)": (1.0 - eps) * rho_air,
    "patch_icpp(1)%alpha(1)": eps,
    "patch_icpp(1)%alpha(2)": 1.0 - eps,
    # TNT detonation-products sphere at the symmetry corner.
    "patch_icpp(2)%geometry": 8,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%x_centroid": 0.0,
    "patch_icpp(2)%y_centroid": 0.0,
    "patch_icpp(2)%z_centroid": 0.0,
    "patch_icpp(2)%radius": charge_radius,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%vel(3)": 0.0,
    "patch_icpp(2)%pres": p_amb if args.mode == "prog_burn" else p_charge,
    "patch_icpp(2)%alpha_rho(1)": (1.0 - eps) * rho_tnt,
    "patch_icpp(2)%alpha_rho(2)": eps * rho_air,
    "patch_icpp(2)%alpha(1)": 1.0 - eps,
    "patch_icpp(2)%alpha(2)": eps,
    # Fluid 1: TNT detonation products with JWL EOS.
    "fluid_pp(1)%eos": 2,
    "fluid_pp(1)%gamma": 1.0 / omega,
    "fluid_pp(1)%pi_inf": 0.0,
    "fluid_pp(1)%cv": 613.5,
    "fluid_pp(1)%jwl_A": A,
    "fluid_pp(1)%jwl_B": B,
    "fluid_pp(1)%jwl_R1": R1,
    "fluid_pp(1)%jwl_R2": R2,
    "fluid_pp(1)%jwl_omega": omega,
    "fluid_pp(1)%jwl_rho0": rho_tnt,
    "fluid_pp(1)%jwl_E0": E0,
    "fluid_pp(1)%jwl_air_e0": 2.5575e5,
    "fluid_pp(1)%jwl_air_rho0": rho_air,
    # Fluid 2: ambient air, ideal gas in MFC's stiffened-gas convention.
    "fluid_pp(2)%eos": 1,
    "fluid_pp(2)%gamma": 1.0 / (gamma_air - 1.0),
    "fluid_pp(2)%pi_inf": 0.0,
    "fluid_pp(2)%cv": 717.5,
}

for i, x in enumerate(probe_x, start=1):
    params[f"probe({i})%x"] = x
    params[f"probe({i})%y"] = probe_yz
    params[f"probe({i})%z"] = probe_yz

if args.mode == "prog_burn":
    params.update(
        {
            # Kinematic front, ignited at the symmetry corner (the charge center).
            "prog_burn": "T",
            "pb_D_cj": pb_D_cj,
            "pb_width": pb_width,
            "pb_x_det": 0.0,
            "pb_y_det": 0.0,
            "pb_z_det": 0.0,
            # Mixing-rate afterburn as products contact the ambient air, same model
            # and rate constants validated in 2D_jwl_prog_burn_afterburn_mixing.
            "jwl_afterburn": "T",
            "jwl_ab_model": 1,
            "jwl_q_ab": 3.0e6,
            "jwl_ab_tau": 3.0e-6,
        }
    )

print(json.dumps(params))
