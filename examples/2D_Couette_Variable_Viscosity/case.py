#!/usr/bin/env python3
"""Variable-viscosity planar Couette flow -- MFC input deck.

Validates the temperature-dependent (Arrhenius) viscosity feature
fluid_pp%visc_model = 1, mu(T) = exp(C + D/T), against the exact steady solution
in reference.py. See README.md for the full description.

Setup (all constants from couette_config.py):
  - 2D channel, periodic in x, walls top and bottom in y.
  - Bottom wall (y=0): no-slip, fixed, isothermal at T0 (cold -> viscous).
  - Top wall   (y=H): no-slip, sliding at U, isothermal at T1 (hot -> thin).
  - Single fluid; model_eqns = 3 + riemann_solver = 2 (HLLC).
  - viscous on, and Fourier conduction (fluid_pp(1)%k_therm > 0) pins the wall
    temperatures and sustains the temperature gradient that drives mu(T).
  - Initial condition: quiescent (u=0), uniform pressure, with a y-varying
    density that encodes the linear conduction temperature profile via the EOS.
    The flow then develops from rest toward the curved steady Couette profile.

The wall-normal grid n is set by the COUETTE_N environment variable so the
validation harness can run a spatial-refinement sweep on a single deck.
"""

import json
import os

import couette_config as cfg

# Wall-normal resolution (refined by the convergence sweep); streamwise grid is
# fixed and small because there is no x-gradient.
# WENO5 needs >= 25 cells per dimension, so the streamwise (x) grid is pinned at
# that minimum (there is no x-gradient) and the sweep keeps n >= 32.
n_cells = int(os.environ.get("COUETTE_N", "64"))
m_cells = int(os.environ.get("COUETTE_M", "24"))

# Square cells: size the (periodic, gradient-free) x-extent so dx == dy on every
# grid. This keeps the acoustic CFL set by dy alone and the time step ~ 1/n.
dy = cfg.H / (n_cells + 1)
Lx = dy * (m_cells + 1)
dx = Lx / (m_cells + 1)

# Time step: minimum of acoustic, viscous, and thermal explicit-stability limits.
c_max = (cfg.gam * (cfg.gam - 1.0) * cfg.cv * cfg.T1) ** 0.5  # sound speed, hot wall
dt_acoustic = 0.4 * min(dx, dy) / c_max
dt_viscous = 0.4 * dy**2 / (2.0 * cfg.nu_max)
alpha_max = cfg.k_therm / (cfg.rho_of_T(cfg.T0) * cfg.cp)  # thermal diffusivity, cold wall
dt_thermal = 0.4 * dy**2 / (2.0 * alpha_max)
dt = min(dt_acoustic, dt_viscous, dt_thermal)

t_step_stop = int(round(cfg.t_end / dt))
if "COUETTE_NSTEPS" in os.environ:  # short override for timing probes
    t_step_stop = int(os.environ["COUETTE_NSTEPS"])
t_step_save = max(1, t_step_stop // 40)  # ~40 snapshots to confirm steadiness

# Linear conduction temperature profile T(y) = T0 + (T1-T0)*y/H imposed through
# the density proxy rho(y) = (p0+p_inf)/((gam-1)*cv*T(y)). Written with plain
# decimals (no 'e' literals -- MFC would expand a bare 'e' to Euler's number).
gm1_cv = (cfg.gam - 1.0) * cfg.cv
slope = (cfg.T1 - cfg.T0) / cfg.H
rho_expr = f"({cfg.p0 + cfg.p_inf})/(({gm1_cv})*({cfg.T0} + ({slope})*y))"

case = {
    # -- Run info --
    "run_time_info": "F",
    # -- Domain (periodic x, walls in y) --
    "x_domain%beg": 0.0,
    "x_domain%end": Lx,
    "y_domain%beg": 0.0,
    "y_domain%end": cfg.H,
    "m": m_cells,
    "n": n_cells,
    "p": 0,
    "cyl_coord": "F",
    # -- Time stepping --
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": t_step_stop,
    "t_step_save": t_step_save,
    # -- Model: 6-equation + HLLC (visc_model = 1 also runs with model_eqns = 2 or HLL) --
    "model_eqns": 3,
    "num_fluids": 1,
    "alt_soundspeed": "F",
    "mpp_lim": "F",
    "mixture_err": "T",
    "time_stepper": 3,
    # -- Spatial scheme --
    "weno_order": 5,
    "weno_eps": 1e-16,
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "F",
    # Compact (non-WENO) viscous gradient. weno_Re_flux is meant for discontinuous
    # viscosity at material interfaces; for a single fluid with smooth mu(T) it
    # leaves the 2*dy (odd-even) mode of the wall-normal velocity undamped -- in
    # this advection-free flow there is no convective upwinding to remove it -- so
    # u(y) develops a grid-scale checkerboard. The compact central gradient damps
    # it and the profile matches the exact solution. (Both must be set together:
    # weno_avg=F with weno_Re_flux=T disables the viscous coupling entirely.)
    "weno_Re_flux": os.environ.get("COUETTE_WENO_RE", "F"),
    "weno_avg": os.environ.get("COUETTE_WENO_AVG", "F"),
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    # -- Boundary conditions: periodic x; no-slip isothermal walls in y --
    "bc_x%beg": -1,
    "bc_x%end": -1,
    "bc_y%beg": -16,  # no-slip wall, bottom (fixed)
    "bc_y%end": -16,  # no-slip wall, top (sliding)
    "bc_y%ve1": cfg.U,  # top-wall x-velocity (the Couette drive)
    "bc_y%isothermal_in": "T",
    "bc_y%Twall_in": cfg.T0,
    "bc_y%isothermal_out": "T",
    "bc_y%Twall_out": cfg.T1,
    # -- Output --
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "parallel_io": "T",
    # -- Physics --
    "viscous": "T",
    "chemistry": "F",
    # -- Fluid: ideal gas with Arrhenius mu(T) = exp(C + D/T) --
    "fluid_pp(1)%gamma": 1.0 / (cfg.gam - 1.0),
    "fluid_pp(1)%pi_inf": cfg.gam * cfg.p_inf / (cfg.gam - 1.0),
    "fluid_pp(1)%cv": cfg.cv,
    "fluid_pp(1)%Re(1)": 1.0 / cfg.mu_ref,  # nominal Re; visc_model=1 overrides with mu(T)
    "fluid_pp(1)%k_therm": cfg.k_therm,
    "fluid_pp(1)%visc_model": 1,
    "fluid_pp(1)%visc_c": cfg.C,
    "fluid_pp(1)%visc_d": cfg.D,
    # -- Initial condition: quiescent, uniform p, linear-T density profile --
    "num_patches": 1,
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": Lx / 2,
    "patch_icpp(1)%y_centroid": cfg.H / 2,
    "patch_icpp(1)%length_x": Lx,
    "patch_icpp(1)%length_y": cfg.H,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": cfg.p0,
    "patch_icpp(1)%alpha_rho(1)": rho_expr,
    "patch_icpp(1)%alpha(1)": 1.0,
}

if __name__ == "__main__":
    print(json.dumps(case))
