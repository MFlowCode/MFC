#!/usr/bin/env python3
"""Automatic ib_neighborhood_radius on a decomposition whose ranks are not cubes.

A thin plate held in a long, narrow channel. The domain is 20 chords long and 8 by 6 across, and the cell
counts (400 x 50 x 50) are chosen so MFC's topology search settles on 16 x 2 x 2 at 64 ranks. That makes the
x ranks 1.25 chords wide while the y and z ranks are 4.0 and 3.0 -- an ordinary situation for a wake, a jet or
a channel, where the flow direction is resolved far more finely than the cross-stream ones.

`ib_neighborhood_radius` is deliberately left unset, so MFC chooses it at start-up and prints

    Automatic choice of ib_neighborhood_radius selected:  N

The radius counts rank hops, and the body must be reachable within that many hops in *every* direction, so
the hop that matters is the one crossing the thinnest rank. The plate's half-extent is 1.301 chords and the
thinnest rank is 1.250 wide, so it needs 2 hops, not 1.

See README.md for the measured before/after and the argument.
"""

import json
import math
import os

Re, Ma, gamma, U, rho = 1000.0, 0.1, 1.4, 1.0, 1.0
cs = U / Ma
P = rho * cs**2 / gamma
c = 1.0  # chord
SPAN = float(os.environ.get("SPAN", 2.4))  # spanwise length of the plate
THICK = 0.1 * c

M = int(os.environ.get("M", 400))  # cells in x; 400/16 ranks = 25, the stencil minimum
N = int(os.environ.get("N", 50))  # cells in y; 50/2 ranks = 25
PZ = int(os.environ.get("PZ", 50))  # cells in z; 50/2 ranks = 25
TSTOP = float(os.environ.get("TSTOP", 0.5))

x0, x1 = -10.0 * c, 10.0 * c
y0, y1 = -4.0 * c, 4.0 * c
z0, z1 = -3.0 * c, 3.0 * c

case = {
    "run_time_info": "T",
    "parallel_io": "T",
    "prim_vars_wrt": "T",
    "ib_state_wrt": "T",
    "format": "silo",
    "precision": "double",
    "x_domain%beg": x0,
    "x_domain%end": x1,
    "y_domain%beg": y0,
    "y_domain%end": y1,
    "z_domain%beg": z0,
    "z_domain%end": z1,
    "m": M - 1,
    "n": N - 1,
    "p": PZ - 1,
    "cyl_coord": "F",
    "cfl_adap_dt": "T",
    "cfl_target": 0.4,
    "n_start": 0,
    "t_save": TSTOP,
    "t_stop": TSTOP,
    "num_patches": 1,
    "num_fluids": 1,
    "model_eqns": "5eq",
    "alt_soundspeed": "F",
    "mpp_lim": "F",
    "mixture_err": "T",
    "time_stepper": "rk3",
    "weno_order": 5,
    "weno_eps": 1.0e-10,
    "weno_Re_flux": "T",
    "weno_avg": "T",
    "avg_state": "arithmetic",
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "F",
    "riemann_solver": "hllc",
    "low_Mach": 2,
    "wave_speeds": "direct",
    "viscous": "T",
    "fd_order": 4,
    "patch_icpp(1)%geometry": 9,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%z_centroid": 0.0,
    "patch_icpp(1)%length_x": x1 - x0,
    "patch_icpp(1)%length_y": y1 - y0,
    "patch_icpp(1)%length_z": z1 - z0,
    "patch_icpp(1)%vel(1)": U,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%vel(3)": 0.0,
    "patch_icpp(1)%pres": P,
    "patch_icpp(1)%alpha_rho(1)": rho,
    "patch_icpp(1)%alpha(1)": 1.0,
    "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
    "fluid_pp(1)%eos": "ideal_gas",
    "fluid_pp(1)%Re(1)": Re / (c * U),
    "bc_x%beg": -7,
    "bc_x%grcbc_in": "T",
    "bc_x%vel_in(1)": U,
    "bc_x%vel_in(2)": 0.0,
    "bc_x%vel_in(3)": 0.0,
    "bc_x%pres_in": P,
    "bc_x%alpha_rho_in(1)": rho,
    "bc_x%alpha_in(1)": 1.0,
    "bc_x%end": -8,
    "bc_y%beg": -8,
    "bc_y%end": -8,
    "bc_z%beg": -8,
    "bc_z%end": -8,
    # ib_neighborhood_radius deliberately NOT set: this case exists to exercise the automatic choice
    "ib": "T",
    "num_ibs": 1,
    "patch_ib(1)%geometry": 9,  # cuboid
    "patch_ib(1)%x_centroid": 0.0,
    "patch_ib(1)%y_centroid": 0.0,
    "patch_ib(1)%z_centroid": 0.0,
    "patch_ib(1)%length_x": c,
    "patch_ib(1)%length_y": SPAN,
    "patch_ib(1)%length_z": THICK,
    "patch_ib(1)%slip": "F",
    "patch_ib(1)%moving_ibm": 0,
}

if __name__ == "__main__":
    if os.environ.get("SUMMARY"):
        bound = 0.5 * math.sqrt(c**2 + SPAN**2 + THICK**2)
        print(f"plate half-extent (s_get_ib_bound, geometry 9): {bound:.4f}")
        print(f"rank extents at 64 ranks (16 x 2 x 2): x {(x1 - x0) / 16:.3f}, " f"y {(y1 - y0) / 2:.3f}, z {(z1 - z0) / 2:.3f}")
        print(f"hops needed across the thinnest rank: ceil(1.1 * {bound:.4f} / {(x1 - x0) / 16:.3f}) = " f"{max(1, math.ceil(1.1 * bound / ((x1 - x0) / 16)))}")
    else:
        print(json.dumps(case, indent=4))
