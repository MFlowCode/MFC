#!/usr/bin/env python3
# 2D high-fidelity JWL++ reactive PETN detonation launching a DENSE 4x4 pack of
# immersed free-flight particles (16 cylinders, near-touching at t=0).
#
#   - PETN charge (full density, unreacted, at ambient p) on the left axis.
#   - A small central hot spot ignites a self-propagating JWL++ reactive burn
#     (jwl_reactive, Souers 2000): dl/dt = G p^b (1-l), sweeping outward at the
#     PETN CJ speed (~8300 m/s). The Garno delta_e puts the unreacted explosive on
#     a stiffer Hugoniot; on this dx=1 mm mesh the reaction zone is sub-cell (no
#     resolved VN spike, expected) but the front self-sustains at ~D_CJ because
#     jwl_G is fast.
#   - The products vent into air and drive a blast onto a dense 4x4 grid of rigid
#     IB cylinders (moving_ibm = 2, two-way coupled), spaced only 5% apart at
#     t=0. Unlike the sparse 4-particle pack, here the front row shields the
#     rear rows from the blast, and soft-sphere IB-IB collisions (collision_model
#     = 1) fire continuously as the compressed cluster is driven downstream and
#     jostles against itself, this is the real stress test for the collision
#     model: dozens of simultaneous near-contacts instead of one clean pair.
#
# Mirror-symmetric about y = 0 (charge, hot spot, and the 4x4 grid all
# symmetric with an even row count), so a correct run keeps net y-force and
# torque at round-off and the cluster drifts in +x only (net y-scatter from
# particle-particle collisions is expected once the cluster starts tumbling).
# Non-reactive air ambient; extrapolation BCs (-3, CBC is prohibited with
# eos_jwl). Sized to run in a few minutes on 8 CPU cores.
import json

from mfc.jwl_products import AIR, PETN, ambient_fluid, jwl_fluid, znd_delta_e

rho0 = PETN["rho0"]
air_rho = AIR["rho0"]
p_amb = 101325.0

# Garno Eq. 17 reactant offset for a stiff-Hugoniot unreacted state.
jwl_delta_e = znd_delta_e(PETN, p_amb)

# JWL++ rate: dl/dt = G p^b (1 - lambda). Fast G keeps the reaction zone sub-cell
# so the detonation self-sustains at the CJ speed on this coarse mesh.
jwl_G = 5.0e-14
jwl_b_exp = 2.0

# --- geometry (metres) ---
CH_X, CH_R = -0.03, 0.015  # charge centre and radius
HS_R = 0.005  # hot-spot radius (still PETN, high pressure)
PR = 0.006  # IB particle radius (smaller than the sparse case
# so a 4x4 grid packs into a similar footprint)
GAP = 1.05  # 5% initial gap between neighbours (not touching,
# so pre_process does not reject overlapping IBs)
SPACING = 2.0 * PR * GAP
NROWS, NCOLS = 4, 4
CLUSTER_X = 0.045  # cluster centroid; ~30 mm standoff from the charge edge
# areal density matched to the sparse 4-particle case (mass / (pi r^2) = 9.945
# kg/m^2 there) so per-particle inertia is physically comparable at this radius.
AREAL_RHO = 2.0e-3 / (3.14159265358979 * 0.008**2)
PMASS = AREAL_RHO * 3.14159265358979 * PR**2

PACK = []
for j in range(NROWS):
    y = (j - 0.5 * (NROWS - 1)) * SPACING
    for i in range(NCOLS):
        x = CLUSTER_X + (i - 0.5 * (NCOLS - 1)) * SPACING
        PACK.append((x, y))

# domain bounds must contain the whole 4x4 cluster's footprint at t=0
PACK_X0 = min(x for x, _ in PACK) - PR
PACK_X1 = max(x for x, _ in PACK) + PR
PACK_Y0 = min(y for _, y in PACK) - PR
PACK_Y1 = max(y for _, y in PACK) + PR

# Sized so the launched cluster stays INSIDE the domain for the whole run: MFC
# does not collide immersed bodies with the domain walls (the boundary is a
# fluid BC, not a rigid stop), so a body given real momentum leaves any finite
# window unless the window is large enough. A dense cluster shares the blast
# impulse across 16x the mass of one sparse particle and the front rows shield
# the rear, so it accelerates less and disperses less than the sparse pack;
# still sized generously since the collision cascade can eject individual
# particles at speeds well above the cluster mean.
XB, XE = -0.05, 0.32
YB, YE = -0.14, 0.14


def ib(i, x, y):
    p = f"patch_ib({i})%"
    return {
        p + "geometry": 2,
        p + "x_centroid": x,
        p + "y_centroid": y,
        p + "radius": PR,
        p + "slip": "F",
        p + "moving_ibm": 2,
        p + "vel(1)": 0.0,
        p + "vel(2)": 0.0,
        p + "angles(1)": 0.0,
        p + "angles(2)": 0.0,
        p + "angles(3)": 0.0,
        p + "angular_vel(1)": 0.0,
        p + "angular_vel(2)": 0.0,
        p + "angular_vel(3)": 0.0,
        p + "mass": PMASS,
    }


ibs = {}
for k, (x, y) in enumerate(PACK, start=1):
    ibs.update(ib(k, x, y))

case = {
    "run_time_info": "T",
    "x_domain%beg": XB,
    "x_domain%end": XE,
    "y_domain%beg": YB,
    "y_domain%end": YE,
    "m": 600,
    "n": 500,
    "p": 0,  # dx = dy = 1 mm
    # Adaptive CFL time stepping: the PETN products over-expand to near-vacuum
    # behind the fast-launched cluster, spiking |u|+c; a fixed dt blows the CFL
    # there, and the collision cascade adds its own stiff transients.
    "cfl_adap_dt": "T",
    "cfl_target": 0.2,
    "n_start": 0,
    "t_stop": 8.0e-5,  # 80 us
    "t_save": 2.0e-6,  # 40 frames
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
    "riemann_solver": 2,  # HLLC required for jwl_reactive
    "wave_speeds": 1,
    "avg_state": 2,
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "bc_y%beg": -3,
    "bc_y%end": -3,
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "parallel_io": "T",
    "jwl_wrt": "T",  # write T, Y_products, lambda arrays
    "ib_state_wrt": "T",
    # self-propagating JWL++ reactive burn
    "jwl_reactive": "T",
    "jwl_G": jwl_G,
    "jwl_b_exp": jwl_b_exp,
    # dense immersed particle pack
    "ib": "T",
    "num_ibs": len(PACK),
    "fd_order": 2,
    # Soft-sphere IB-IB collisions: at 5% initial gap the blast slams the front
    # row into the rear rows almost immediately, so this cluster runs through
    # many simultaneous contacts per step (unlike the sparse case's single
    # pairwise collision). Keep the same soft, damped contact tuning that held
    # the adaptive CFL there; a stiffer spring is even more likely to blow up
    # with this many concurrent contacts.
    "collision_model": 1,  # soft-sphere (the only model MFC ships)
    "coefficient_of_restitution": 0.98,  # damped rebound, dissipates impact energy
    "collision_time": 3.0e-6,  # soft enough to hold the adaptive dt
    "ib_coefficient_of_friction": 0.1,
    # Patch 1: ambient air filling the domain.
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.5 * (XB + XE),
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": XE - XB,
    "patch_icpp(1)%length_y": YE - YB,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": p_amb,
    "patch_icpp(1)%alpha_rho(1)": rho0 * 1.0e-8,
    "patch_icpp(1)%alpha_rho(2)": air_rho * (1.0 - 1.0e-8),
    "patch_icpp(1)%alpha(1)": 1.0e-8,
    "patch_icpp(1)%alpha(2)": 0.99999999,
    # Patch 2: PETN charge disc, unreacted, at ambient pressure.
    "patch_icpp(2)%geometry": 2,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%x_centroid": CH_X,
    "patch_icpp(2)%y_centroid": 0.0,
    "patch_icpp(2)%radius": CH_R,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": p_amb,
    "patch_icpp(2)%alpha_rho(1)": rho0 * (1.0 - 1.0e-8),
    "patch_icpp(2)%alpha_rho(2)": air_rho * 1.0e-8,
    "patch_icpp(2)%alpha(1)": 0.99999999,
    "patch_icpp(2)%alpha(2)": 1.0e-8,
    # Patch 3: central hot spot (still PETN) igniting the burn.
    "patch_icpp(3)%geometry": 2,
    "patch_icpp(3)%alter_patch(1)": "T",
    "patch_icpp(3)%alter_patch(2)": "T",
    "patch_icpp(3)%x_centroid": CH_X,
    "patch_icpp(3)%y_centroid": 0.0,
    "patch_icpp(3)%radius": HS_R,
    "patch_icpp(3)%vel(1)": 0.0,
    "patch_icpp(3)%vel(2)": 0.0,
    "patch_icpp(3)%pres": 40.0e9,
    "patch_icpp(3)%alpha_rho(1)": rho0 * (1.0 - 1.0e-8),
    "patch_icpp(3)%alpha_rho(2)": air_rho * 1.0e-8,
    "patch_icpp(3)%alpha(1)": 0.99999999,
    "patch_icpp(3)%alpha(2)": 1.0e-8,
    # Fluid 1: PETN JWL products (+ air references); Fluid 2: air.
    **jwl_fluid(1, PETN, AIR, delta_e=jwl_delta_e),
    **ambient_fluid(2, AIR),
}
case.update(ibs)
print(json.dumps(case))
