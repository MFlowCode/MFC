#!/usr/bin/env python3
"""
2-D lumped-wall example: a hot, inert cylinder cooling in a cold gas stream.

The cylinder is an immersed boundary with thermal_bc = 3: one solid body at a
single temperature, advanced by its heat balance with the gas. Its heat capacity
is scaled down so the temperature relaxes within a few flow-through times; see
README.md. ib_heat_wrt records the heat the body exchanges every step to
D/ib_heat.dat, which heat_balance.py reads.

Required mechanism file (shared with 2D_ibm_reacting_surface):
    carbon_gasphase_reduced_gri11.yaml
"""

import json
import os
import sys

import cantera as ct

#
# USER PARAMETERS
#

current_dir = os.path.dirname(os.path.abspath(__file__))
ctfile = os.path.normpath(os.path.join(current_dir, "..", "2D_ibm_reacting_surface", "carbon_gasphase_reduced_gri11.yaml"))

# Freestream
p_inf = 101325.0  # Pa
T_inf = 300.0  # K
u_inf = 20.0  # m/s
Y_inf = {"N2": 0.760557, "O2": 0.233000, "H2O": 0.005800, "CO2": 0.000643}

# Lumped cylinder: initial temperature, and the fraction of the run its thermal time takes
Twall0 = 900.0  # K
tau_ratio = 0.25  # thermal time / run time (estimated, see below)
rho_solid = 1000.0  # kg/m^3

# Geometry and grid
D = 0.002  # cylinder diameter (m)
Lx = 7.0 * D
Ly = 5.0 * D
cells_per_D = 24
cfl = 0.8
thermal_layer = 0.25 * D  # initial conduction layer around the hot cylinder

# Run length in flow-through times D/u_inf
Tstar = 20.0
frames = 10

#
# FREESTREAM AND THERMAL TIME ESTIMATE
#

gas = ct.Solution(ctfile)
gas.TPY = T_inf, p_inf, Y_inf
Y_inf_vec = gas.Y.copy()
a_inf = gas.sound_speed
rho_inf = gas.density
mu_inf = gas.viscosity

# Film-temperature properties and the Churchill-Bernstein Nusselt number for a cylinder in crossflow
gas.TPY = 0.5 * (T_inf + Twall0), p_inf, Y_inf
Re_D = gas.density * u_inf * D / gas.viscosity
Pr = gas.cp_mass * gas.viscosity / gas.thermal_conductivity
Nu = 0.3 + 0.62 * Re_D**0.5 * Pr ** (1.0 / 3.0) / (1.0 + (0.4 / Pr) ** (2.0 / 3.0)) ** 0.25 * (1.0 + (Re_D / 282000.0) ** (5.0 / 8.0)) ** 0.8
h_est = Nu * gas.thermal_conductivity / D

Tend = Tstar * D / u_inf
t_save = Tend / frames

# m c_s dT/dt = -h A (T - T_gas) with m/A = rho_solid D/4 per unit depth: choose c_s so the thermal time is tau_ratio * Tend
tau = tau_ratio * Tend
cp_solid = tau * h_est / (rho_solid * D / 4.0)

dx = D / cells_per_D
Nx = int(round(Lx / dx))
Ny = int(round(Ly / dx))

print(
    f"""Lumped cylinder in crossflow
----------------------------
Re_D = {Re_D:.1f}   Pr = {Pr:.3f}   Nu (Churchill-Bernstein) = {Nu:.2f}
h (estimate)          = {h_est:.1f} W/m^2/K
rho_solid, cp_solid   = {rho_solid:.1f} kg/m^3, {cp_solid:.3f} J/kg/K  (scaled so the thermal time is {tau:.3e} s)
Run time              = {Tend:.3e} s   Grid = {Nx} x {Ny}
""",
    file=sys.stderr,
)

#
# MFC CASE
#

case = {
    "run_time_info": "T",
    # Domain
    "x_domain%beg": 0.0,
    "x_domain%end": Lx,
    "y_domain%beg": 0.0,
    "y_domain%end": Ly,
    "m": Nx - 1,
    "n": Ny - 1,
    "p": 0,
    "cyl_coord": "F",
    # Time integration
    "cfl_adap_dt": "T",
    "cfl_target": cfl,
    "n_start": 0,
    "t_save": t_save,
    "t_stop": Tend,
    # Numerics
    "num_patches": 1,
    "num_fluids": 1,
    "model_eqns": 2,
    "alt_soundspeed": "F",
    "mpp_lim": "F",
    "mixture_err": "T",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1.0e-10,
    "weno_Re_flux": "T",
    "weno_avg": "T",
    "avg_state": 2,
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "F",
    "riemann_solver": 2,
    "low_Mach": 0,
    "wave_speeds": 1,
    "viscous": "T",
    "fd_order": 4,
    # Characteristic subsonic inflow at x-beg, outflow at x-end, periodic in y
    "bc_x%beg": -7,
    "bc_x%end": -3,
    "bc_y%beg": -1,
    "bc_y%end": -1,
    "bc_x%grcbc_in": "T",
    "bc_x%vel_in(1)": u_inf,
    "bc_x%vel_in(2)": 0.0,
    "bc_x%pres_in": p_inf,
    "bc_x%alpha_rho_in(1)": rho_inf,
    "bc_x%alpha_in(1)": 1.0,
    # Chemistry: gas-phase transport carries the heat; no reactions at these temperatures
    "chemistry": "T",
    "chem_params%diffusion": "T",
    "chem_params%reactions": "F",
    "chem_params%transport_model": 1,
    "cantera_file": ctfile,
    "chem_wrt_T": "T",
    # Initial freestream
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.5 * Lx,
    "patch_icpp(1)%y_centroid": 0.5 * Ly,
    "patch_icpp(1)%length_x": Lx,
    "patch_icpp(1)%length_y": Ly,
    "patch_icpp(1)%vel(1)": u_inf,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": p_inf,
    "patch_icpp(1)%alpha_rho(1)": rho_inf,
    "patch_icpp(1)%alpha(1)": 1.0,
    # Lumped cylinder, offset slightly from the symmetry line to seed the wake instability
    "ib": "T",
    "num_ibs": 1,
    "patch_ib(1)%geometry": 2,
    "patch_ib(1)%x_centroid": 2.0 * D,
    "patch_ib(1)%y_centroid": 0.5 * Ly + dx / 3.0,
    "patch_ib(1)%radius": 0.5 * D,
    "patch_ib(1)%slip": "F",
    "patch_ib(1)%thermal_bc": 3,
    "patch_ib(1)%Twall": Twall0,
    "patch_ib(1)%thermal_layer": thermal_layer,
    "patch_ib(1)%rho_solid": rho_solid,
    "patch_ib(1)%cp_solid": cp_solid,
    # Fluid properties
    "fluid_pp(1)%gamma": 1.0 / (1.4 - 1.0),
    "fluid_pp(1)%eos": "stiffened_gas",
    "fluid_pp(1)%pi_inf": 0.0,
    "fluid_pp(1)%Re(1)": 1.0 / mu_inf,
    # Output
    "parallel_io": "T",
    "format": "silo",
    "precision": "single",
    "prim_vars_wrt": "T",
    "ib_heat_wrt": "T",
}

for i in range(len(Y_inf_vec)):
    case[f"patch_icpp(1)%Y({i + 1})"] = float(Y_inf_vec[i])

if __name__ == "__main__":
    print(json.dumps(case))
