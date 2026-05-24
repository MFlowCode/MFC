#!/usr/bin/env python3
# 3D tri-periodic advection of a cubic lattice of Lagrangian bubbles
# with log-normally distributed radii. Background flow is uniform along
# the (1,1,1) diagonal.
import json
import math
import os

import numpy as np

# Reference values for nondimensionalization
L0 = 1e-3  # length - m
rho0 = 950  # density - kg/m3
c0 = 1449.0  # speed of sound - m/s
p0 = rho0 * c0 * c0  # pressure - Pa
T0 = 298  # temperature - K

# Domain: cube [-1, 1] mm, nondimensionalized to [-1, 1]
xb, xe = -1.0e-3 / L0, 1.0e-3 / L0
yb, ye = -1.0e-3 / L0, 1.0e-3 / L0
zb, ze = -1.0e-3 / L0, 1.0e-3 / L0

# Host (water)
gamma_host = 6.12
pi_inf_host = 3.43e8
mu_host = 0.001
rho_host = 950

# Lagrangian bubbles' properties (matched to the reference 3D lag case)
R_uni = 8314
MW_g = 28.0
MW_v = 18.0
gamma_g = 1.4
gamma_v = 1.333
pv = 2350
cp_g = 1.0e3
cp_v = 2.1e3
k_g = 0.025
k_v = 0.02
diffVapor = 2.5e-5
sigBubble = 0.020
mu_g = 1.48e-5
rho_g = 1

patm = 1e5
rho_host_nd = rho_host / rho0
rho_g_nd = rho_g / rho0
pres = patm / p0

# Diagonal advection: equal (nondim) velocity components. ~Mach 0.01.
u_adv = 0.05
v_adv = 0.05
w_adv = 0.05

# Time stepping. With m=n=p=49 (dx=0.04) and unit nondim sound speed,
# CFL=0.5 at dt=0.02. One full traversal of the domain along x at
# u_adv=0.01 takes (xe-xb)/u_adv = 200 nondim time units.
Ncells = 50
m = n = p = Ncells - 1
dx = (xe - xb) / Ncells
dt = 0.02
t_traverse = (xe - xb) / u_adv
t_step_stop = int(round(t_traverse / dt))
t_step_save = t_step_stop // 100

# Lattice and log-normal radius distribution
N_side = 2
nBubs = N_side**3
R_median_nd = 0.02  # 20 micron median radius
sigma_ln = 0.3

data = {
    "run_time_info": "T",
    "x_domain%beg": xb,
    "x_domain%end": xe,
    "y_domain%beg": yb,
    "y_domain%end": ye,
    "z_domain%beg": zb,
    "z_domain%end": ze,
    "m": m,
    "n": n,
    "p": p,
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": t_step_stop,
    "t_step_save": t_step_save,
    "model_eqns": 2,
    "alt_soundspeed": "F",
    "mixture_err": "F",
    "mpp_lim": "T",
    "time_stepper": 3,
    "weno_order": 5,
    "mapped_weno": "T",
    "mp_weno": "F",
    "avg_state": 2,
    "weno_eps": 1e-16,
    "riemann_solver": 2,
    "wave_speeds": 1,
    # Tri-periodic boundaries
    "bc_x%beg": -1,
    "bc_x%end": -1,
    "bc_y%beg": -1,
    "bc_y%end": -1,
    "bc_z%beg": -1,
    "bc_z%end": -1,
    "num_patches": 1,
    "num_fluids": 2,
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "parallel_io": "T",
    "lag_db_wrt": "T",
    "lag_rad_wrt": "T",
    # Fluid 1: water host
    "fluid_pp(1)%gamma": 1.0 / (gamma_host - 1.0),
    "fluid_pp(1)%pi_inf": gamma_host * (pi_inf_host / p0) / (gamma_host - 1.0),
    # Fluid 2: gas (inside bubbles)
    "fluid_pp(2)%gamma": 1.0 / (gamma_g - 1.0),
    "fluid_pp(2)%pi_inf": 0.0,
    # Bubble physical parameters
    "bub_pp%R0ref": 1.0,
    "bub_pp%p0ref": 1.0,
    "bub_pp%rho0ref": 1.0,
    "bub_pp%T0ref": 1.0,
    "bub_pp%ss": sigBubble / (rho0 * L0 * c0 * c0),
    "bub_pp%pv": pv / p0,
    "bub_pp%vd": diffVapor / (L0 * c0),
    "bub_pp%mu_l": mu_host / (rho0 * L0 * c0),
    "bub_pp%gam_v": gamma_v,
    "bub_pp%gam_g": gamma_g,
    "bub_pp%M_v": MW_v,
    "bub_pp%M_g": MW_g,
    "bub_pp%k_v": k_v * (T0 / (L0 * rho0 * c0 * c0 * c0)),
    "bub_pp%k_g": k_g * (T0 / (L0 * rho0 * c0 * c0 * c0)),
    "bub_pp%cp_v": cp_v * (T0 / (c0 * c0)),
    "bub_pp%cp_g": cp_g * (T0 / (c0 * c0)),
    "bub_pp%R_v": (R_uni / MW_v) * (T0 / (c0 * c0)),
    "bub_pp%R_g": (R_uni / MW_g) * (T0 / (c0 * c0)),
    # Viscosity
    "viscous": "T",
    "fluid_pp(1)%Re(1)": 1.0 / (mu_host / (rho0 * c0 * L0)),
    "fluid_pp(2)%Re(1)": 1.0 / (mu_g / (rho0 * c0 * L0)),
    # Single patch: uniform background flow advecting diagonally
    "patch_icpp(1)%geometry": 9,
    "patch_icpp(1)%x_centroid": (xb + xe) / 2,
    "patch_icpp(1)%y_centroid": (yb + ye) / 2,
    "patch_icpp(1)%z_centroid": (zb + ze) / 2,
    "patch_icpp(1)%length_x": (xe - xb),
    "patch_icpp(1)%length_y": (ye - yb),
    "patch_icpp(1)%length_z": (ze - zb),
    "patch_icpp(1)%vel(1)": u_adv,
    "patch_icpp(1)%vel(2)": v_adv,
    "patch_icpp(1)%vel(3)": w_adv,
    "patch_icpp(1)%pres": pres,
    "patch_icpp(1)%alpha_rho(1)": rho_host_nd,
    "patch_icpp(1)%alpha_rho(2)": 0,
    "patch_icpp(1)%alpha(1)": 1,
    "patch_icpp(1)%alpha(2)": 0,
    # Lagrangian bubbles
    "bubbles_lagrange": "T",
    "fd_order": 4,
    "bubble_model": 3,
    "thermal": 3,
    "polytropic": "F",
    "adap_dt": "T",
    "lag_params%nBubs_glb": nBubs,
    "lag_params%vel_model": 1,
    "lag_params%drag_model": 1,
    "lag_params%solver_approach": 1,
    "lag_params%cluster_type": 2,
    "lag_params%pressure_corrector": "F",
    "lag_params%smooth_type": 1,
    "lag_params%epsilonb": 1.0,
    "lag_params%valmaxvoid": 0.9,
    "lag_params%write_bubbles": "T",
    "lag_params%write_bubbles_stats": "T",
    "lag_params%write_void_evol": "T",
}


def write_lag_bubbles_file():
    # Cubic lattice spanning the periodic domain. Lattice spacing equals
    # the domain length divided by N_side; positions are shifted by
    # spacing/2 so that no bubble sits exactly on a periodic face.
    Lx = xe - xb
    Ly = ye - yb
    Lz = ze - zb
    sx = Lx / N_side
    sy = Ly / N_side
    sz = Lz / N_side

    xs = xb + 0.5 * sx + np.arange(N_side) * sx
    ys = yb + 0.5 * sy + np.arange(N_side) * sy
    zs = zb + 0.5 * sz + np.arange(N_side) * sz
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    x_pos = X.flatten()
    y_pos = Y.flatten()
    z_pos = Z.flatten()

    # Log-normal radii (nondim by L0). Clip the tails so a freak draw
    # cannot overlap with neighbors on the lattice.
    rng = np.random.default_rng(seed=42)
    R = rng.lognormal(mean=math.log(R_median_nd), sigma=sigma_ln, size=nBubs)
    R_max_allowed = 0.45 * min(sx, sy, sz)
    R = np.minimum(R, R_max_allowed)

    vx = 0.0 * np.full(nBubs, u_adv)
    vy = 0.0 * np.full(nBubs, v_adv)
    vz = 0.0 * np.full(nBubs, w_adv)
    # Small initial radial velocity to kick off bubble dynamics.
    # For R~0.02 at atmospheric pres, Minnaert ω_n ≈ 0.73 (nondim),
    # so rdot=1e-3 → ε/R₀ ≈ rdot/(R₀·ω_n) ≈ 7% amplitude oscillations.
    rdot = 0.0 * np.full(nBubs, 1.0e-3)

    input_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "input")
    os.makedirs(input_dir, exist_ok=True)
    with open(os.path.join(input_dir, "lag_bubbles.dat"), "w") as f:
        for i in range(nBubs):
            f.write(f"{x_pos[i]:12.6f}\t{y_pos[i]:12.6f}\t{z_pos[i]:12.6f}\t{vx[i]:12.6f}\t{vy[i]:12.6f}\t{vz[i]:12.6f}\t{R[i]:12.6f}\t{rdot[i]:12.6f}\n")


write_lag_bubbles_file()

print(json.dumps(data, indent=4))
