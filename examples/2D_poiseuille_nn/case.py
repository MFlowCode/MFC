#!/usr/bin/env python3
"""
2D Poiseuille Flow with Herschel-Bulkley Non-Newtonian Fluid

Pressure-driven channel flow between two no-slip walls, driven by a constant
body force in the streamwise (x) direction. Periodic BCs in x, no-slip in y.

HB model: mu = (tau0/gdot)*(1 - exp(-m*gdot)) + K * gdot^(n-1)
  - tau0: yield stress
  - K:    consistency index
  - n:    flow behavior index (< 1 shear-thinning, > 1 shear-thickening)
  - m:    Papanastasiou regularization parameter

For Newtonian Poiseuille validation, set tau0=0, nn=1, K=mu.
The analytical solution is: u(y) = (G/(2*mu)) * y * (H - y)
where G = rho * g_x is the effective pressure gradient.
"""
import json
import math

# === Channel geometry (square domain) ===
L = 1.0         # Channel length (streamwise, x)
H = 1.0         # Channel height (wall-normal, y)

# === Grid resolution ===
Nx = 24         # Cells in x (streamwise, minimal — periodic)
Ny = 81         # Cells in y (wall-normal)

# === Fluid properties ===
rho = 1.0       # Density
p0 = 1e5        # Reference pressure (high for low Mach)
gamma = 1.4     # Ratio of specific heats

# Sound speed and CFL
c = math.sqrt(gamma * p0 / rho)
dx = L / (Nx + 1)
cfl = 0.3
dt = cfl * dx / c

# === Body force (pressure gradient substitute) ===
# G = rho * g_x acts as dp/dx driving force
g_x = 0.5

# === HB non-Newtonian model parameters ===
tau0 = 0.0          # Yield stress (set 0 for power-law)
K = 0.1            # Consistency index
nn = 2.0          # Flow behavior index (< 1 = shear-thinning)
hb_m = 1000.0       # Papanastasiou regularization parameter
mu_min = 1e-4       # Minimum viscosity bound
mu_max = 10.0       # Maximum viscosity bound
mu_bulk = 0.0       # Bulk viscosity

# Reference Re based on consistency index (used as baseline)
Re_ref = 1.0 / K    # = 100

# === Time control ===
t_end = 10.0        # End time (allow flow to reach steady state)
t_save = 5.0        # Save interval

eps = 1e-6

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0,
            "x_domain%end": L,
            "y_domain%beg": 0.0,
            "y_domain%end": H,
            "m": Nx,
            "n": Ny,
            "p": 0,
            "cfl_adap_dt": "T",
            "cfl_target": cfl,
            "n_start": 0,
            "t_stop": t_end,
            "t_save": t_save,
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 2,
            "mpp_lim": "F",
            "mixture_err": "T",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1e-16,
            "mapped_weno": "T",
            "weno_Re_flux": "T",
            "mp_weno": "T",
            "weno_avg": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            # Boundary Conditions
            # Periodic in x (streamwise), no-slip walls in y
            "bc_x%beg": -1,
            "bc_x%end": -1,
            "bc_y%beg": -16,
            "bc_y%end": -16,
            # Viscous
            "viscous": "T",
            # Body Force (drives the flow like a pressure gradient)
            "bf_x": "T",
            "g_x": g_x,
            "k_x": 0.0,
            "w_x": 0.0,
            "p_x": 0.0,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "omega_wrt(3)": "T",
            "fd_order": 4,
            "parallel_io": "T",
            # Patch 1: Entire channel domain (initially at rest)
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": L / 2.0,
            "patch_icpp(1)%y_centroid": H / 2.0,
            "patch_icpp(1)%length_x": L,
            "patch_icpp(1)%length_y": H,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": p0,
            "patch_icpp(1)%alpha_rho(1)": rho * 0.5,
            "patch_icpp(1)%alpha(1)": 0.5,
            "patch_icpp(1)%alpha_rho(2)": rho * 0.5,
            "patch_icpp(1)%alpha(2)": 0.5,
            # Fluid 1: HB non-Newtonian fluid
            "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
            "fluid_pp(1)%pi_inf": 0.0,
            "fluid_pp(1)%Re(1)": Re_ref,
            "fluid_pp(1)%non_newtonian": "T",
            "fluid_pp(1)%tau0": tau0,
            "fluid_pp(1)%K": K,
            "fluid_pp(1)%nn": nn,
            "fluid_pp(1)%hb_m": hb_m,
            "fluid_pp(1)%mu_min": mu_min,
            "fluid_pp(1)%mu_max": mu_max,
            "fluid_pp(1)%mu_bulk": mu_bulk,
            # Fluid 2: same properties (single-phase effectively)
            "fluid_pp(2)%gamma": 1.0 / (gamma - 1.0),
            "fluid_pp(2)%pi_inf": 0.0,
            "fluid_pp(2)%Re(1)": Re_ref,
            "fluid_pp(2)%non_newtonian": "T",
            "fluid_pp(2)%tau0": tau0,
            "fluid_pp(2)%K": K,
            "fluid_pp(2)%nn": nn,
            "fluid_pp(2)%hb_m": hb_m,
            "fluid_pp(2)%mu_min": mu_min,
            "fluid_pp(2)%mu_max": mu_max,
            "fluid_pp(2)%mu_bulk": mu_bulk,
        }
    )
)
