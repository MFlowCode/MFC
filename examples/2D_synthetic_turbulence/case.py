#!/usr/bin/env python3
"""
2D synthetic turbulence forcing demonstration.

A uniform air stream flows in the +x direction through a rectangular channel.
A Gaussian synthetic-turbulence source is placed near the inlet and injects
broadband velocity fluctuations that advect downstream at the mean flow speed.

Physical setup
--------------
  Domain      : x in [0, 4] m,  y in [-1, 1] m
  Background  : air at rest (p=101325 Pa, rho=1.225 kg/m^3), U=10 m/s in x
  Source zone : Gaussian centered at (x=1, y=0), half-widths Lx=Ly=0.6 m
  Turbulence  : 3 energy shells, wave vectors drawn uniformly on the unit circle

Energy shell layout
-------------------
  Shell 1  k = 2*pi/0.6 ~ 10.5  rad/m   4 waves   A = 0.5 m/s
  Shell 2  k = 2*pi/0.3 ~ 21.0  rad/m   6 waves   A = 0.3 m/s
  Shell 3  k = 2*pi/0.15 ~ 42.0 rad/m  10 waves   A = 0.1 m/s

All wave directions are random (seed=42); the field advects in +x at U_inf.
"""

import json
import math

# --------------------------------------------------------------------------
# Physical constants
# --------------------------------------------------------------------------
gam   = 1.4
p0    = 101325.0    # Pa
rho0  = 1.225       # kg/m^3
c0    = math.sqrt(gam * p0 / rho0)   # ~340 m/s
U_inf = 1.0        # mean +x velocity (M ~ 0.03, well subsonic)

# --------------------------------------------------------------------------
# Domain
# --------------------------------------------------------------------------
Lx, Ly = 4.0, 2.0
Nx, Ny = 400, 200   # dx = dy = 0.02 m

dx = Lx / Nx
cfl = 0.8
dt  = cfl * dx / c0

t_end  = 2.0 * Lx / U_inf  # ~0.2 s  (half flow-through)
t_save = t_end / 100.0

# --------------------------------------------------------------------------
# Synthetic turbulence shell parameters
# --------------------------------------------------------------------------
n_shells = 3

k_shells   = [2 * math.pi / 0.60,   # shell 1 — large eddies
              2 * math.pi / 0.30,   # shell 2 — medium eddies
              2 * math.pi / 0.15]   # shell 3 — small eddies

amp_shells = [0.2,   # m/s  ~5 % of U_inf
              0.1,   # m/s  ~3 %
              0.05]   # m/s  ~1 %

waves_per_shell = [3, 5, 6]   # random directions per shell

# Gaussian source zone
src_x, src_y = 0.4, 0.0   # center
src_Lx, src_Ly = 0.3, 1.5  # full extents fed to synth_L

# --------------------------------------------------------------------------
# Build the case dict
# --------------------------------------------------------------------------
case = {
    # -- Logistics --
    "run_time_info": "T",

    # -- Domain --
    "x_domain%beg": 0.0,
    "x_domain%end": Lx,
    "y_domain%beg": -Ly / 2,
    "y_domain%end":  Ly / 2,
    "m": Nx,
    "n": Ny,
    "p": 0,

    # -- Time stepping --
    "dt": dt,
    "cfl_adap_dt": "T",
    "cfl_target": cfl,
    "t_stop":  t_end,
    "t_save":  t_save,
    "n_start": 0,

    # -- Numerics --
    "num_patches":   1,
    "model_eqns":    2,
    "num_fluids":    1,
    "mpp_lim":       "F",
    "mixture_err":   "T",
    "time_stepper":  3,
    "weno_order":    5,
    "weno_eps":      1.0e-16,
    "weno_Re_flux": "T",
    "mapped_weno":   "T",
    "mp_weno":       "F",
    "weno_avg":      "F",
    "null_weights":  "F",
    "riemann_solver": 2,
    "wave_speeds":   1,
    "avg_state":     2,
    "viscous": "T",

    # -- Boundary conditions --
    # Extrapolation on x (open channel), periodic in y
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "bc_y%beg": -1,
    "bc_y%end": -1,

    # -- Output --
    "format":        1,
    "precision":     2,
    "prim_vars_wrt": "T",
    "parallel_io":   "T",

    # -- Initial condition: uniform stream --
    "patch_icpp(1)%geometry":   3,
    "patch_icpp(1)%x_centroid": Lx / 2,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x":   Lx,
    "patch_icpp(1)%length_y":   Ly,
    "patch_icpp(1)%vel(1)":     U_inf,
    "patch_icpp(1)%vel(2)":     0.0,
    "patch_icpp(1)%pres":       p0,
    "patch_icpp(1)%alpha_rho(1)": rho0,
    "patch_icpp(1)%alpha(1)":   1.0,

    # -- Fluid --
    "fluid_pp(1)%gamma":  1.0 / (gam - 1.0),
    "fluid_pp(1)%pi_inf": 0.0,

    # -- Synthetic turbulence --
    "synthetic_turbulence":  "T",
    "synth_seed":            42,
    "synth_n_shells":        n_shells,
    "synth_U_inf":           U_inf,
    "num_turbulent_sources": 1,

    # Energy shells (wave-number magnitudes, amplitudes, wave counts)
    "synth_k_shell(1)":          k_shells[0],
    "synth_k_shell(2)":          k_shells[1],
    "synth_k_shell(3)":          k_shells[2],
    "synth_amp_shell(1)":        amp_shells[0],
    "synth_amp_shell(2)":        amp_shells[1],
    "synth_amp_shell(3)":        amp_shells[2],
    "synth_n_waves_per_shell(1)": waves_per_shell[0],
    "synth_n_waves_per_shell(2)": waves_per_shell[1],
    "synth_n_waves_per_shell(3)": waves_per_shell[2],

    # Gaussian forcing zone (source 1)
    "turb_pos(1,1)": src_x,
    "turb_pos(1,2)": src_y,
    "turb_pos(1,3)": 0.0,   # z not used in 2-D
    "synth_L(1,1)":  src_Lx,
    "synth_L(1,2)":  src_Ly,
    "synth_L(1,3)":  1.0,   # z extent irrelevant in 2-D; set nonzero to avoid div-by-zero

    # Patch: Airfoil Immersed Boundary
    "patch_ib(1)%geometry": 4,
    "patch_ib(1)%x_centroid": 2.0,
    "patch_ib(1)%y_centroid": 0.0,
    "patch_ib(1)%length_z": 0.75,
    "patch_ib(1)%c": 0.5,
    "patch_ib(1)%t": 0.15,
    "patch_ib(1)%p": 0.4,
    "patch_ib(1)%m": 0.02,
    "patch_ib(1)%angular_vel(3)": "-10.0 * 1.0 * pi * sin(1.0 * pi * t) * pi / 180.",
    "patch_ib(1)%angles(3)": -5.0 * math.pi / 180.,
    "patch_ib(1)%moving_ibm": 1,
    "patch_ib(1)%slip": "F",

    "fluid_pp(1)%gamma": 1.0e00 / (gam - 1.0e00),  # 2.50(Not 1.40)
    "fluid_pp(1)%pi_inf": 0,
    "fluid_pp(1)%Re(1)": 100000,
}

print(json.dumps(case))
