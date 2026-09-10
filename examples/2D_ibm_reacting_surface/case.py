#!/usr/bin/env python3
"""
2-D heterogeneous reacting-surface example.

An oxidizing gas mixture flows past a stationary circular carbon cylinder
represented by an immersed boundary. Heterogeneous surface reactions consume
gas-phase reactants and produce carbon-containing species at the cylinder
surface.

This example exercises:
    - immersed-boundary surface chemistry
    - species diffusion
    - homogeneous gas-phase chemistry
    - heterogeneous surface species-flux boundary conditions
    - prescribed surface temperature

Gas-phase chemistry:
    reduced 11-species GRI-based mechanism

Surface chemistry:
    heterogeneous carbon oxidation mechanism

Required mechanism files:
    carbon_gasphase_reduced_gri11.yaml
    carbon_surface_bradley_11species.yaml
"""

import json
import os
import sys

import cantera as ct

#
# USER PARAMETERS
#

# Chemistry mechanisms
current_dir = os.path.dirname(os.path.abspath(__file__))
ctfile = os.path.join(current_dir, "carbon_gasphase_reduced_gri11.yaml")
surface_ctfile = os.path.join(current_dir, "carbon_surface_bradley_11species.yaml")

# Reacting immersed-boundary surface
# thermal_bc: 0 = zero normal temperature gradient
#             1 = prescribed surface temperature
#             2 = coupled surface energy balance
thermal_bc = 1
surface_reaction = 1  # 0 = disabled; 1 = enabled

# Freestream conditions
p_inf = 101325.0  # Pa
T_inf = 298.0  # K
u_inf = 18.9  # m/s

# Freestream mass fractions
Y_inf = {
    "N2": 0.760557,
    "O2": 0.233000,
    "H2O": 0.005800,
    "CO2": 0.000643,
}

# Prescribed cylinder surface temperature for thermal_bc = 1
Twall = 1200.0  # K

# Cylinder and domain geometry
D = 0.002  # cylinder diameter (m)
Lx = 7.0 * D
Ly = 6.0 * D
x_cyl = 2.0 * D
y_cyl = Ly / 2.0
r_cyl = 0.5 * D

# Grid and simulation duration
cells_per_D = 40
cfl = 0.8

# Convective nondimensional time: Tstar = t * u_inf / D
Tstar = 12.0
frames = 10

# Chemistry and transport
chemistry = "T"
diffusion = "T"
reactions = "T"
transport_model = 1
substeps = 2
adap = "F"
substeps_max = 8

# Numerics
low_Mach = 0
fd_order = 4

#
# FREESTREAM AND FILM PROPERTIES
#

gas = ct.Solution(ctfile)

# Freestream properties
gas.TPY = T_inf, p_inf, Y_inf
Y_inf_vec = gas.Y.copy()
a_inf = gas.sound_speed
rho_inf = gas.density
mu_inf = gas.viscosity
Mach_inf = u_inf / a_inf
Re_D_inf = rho_inf * u_inf * D / mu_inf

# Film properties evaluated at T_film = (T_inf + Twall)/2
T_film = 0.5 * (T_inf + Twall)
gas.TPY = T_film, p_inf, Y_inf
rho_film = gas.density
mu_film = gas.viscosity
Re_D_film = rho_film * u_inf * D / mu_film

#
# GRID
#

dx = D / cells_per_D
dy = dx
Nx = int(round(Lx / dx))
Ny = int(round(Ly / dy))
m = Nx - 1
n = Ny - 1

#
# TIME INTEGRATION
#

# dt and NT provide diagnostic estimates only; MFC uses CFL-based
# adaptive time stepping for the actual simulation.
dt_est = cfl * dx / (u_inf + a_inf)
Tend = Tstar * D / u_inf
NT_est = int(Tend / dt_est)
t_save = Tend / frames

#
# DIAGNOSTICS
#

print(
    f"""Reacting carbon-cylinder example
--------------------------------
Gas mechanism       = {os.path.basename(ctfile)}
Surface mechanism   = {os.path.basename(surface_ctfile)}
Freestream           = {Y_inf}
Pressure             = {p_inf:.1f} Pa
Freestream T         = {T_inf:.1f} K
Film T               = {T_film:.1f} K
Mach number          = {Mach_inf:.4f}
Freestream velocity  = {u_inf:.4f} m/s
Re_D (freestream)    = {Re_D_inf:.1f}
Re_D (film)          = {Re_D_film:.1f}
Surface T            = {Twall:.1f} K
Cylinder diameter    = {D:.4e} m
Cells across D       = {cells_per_D}
Grid                 = {Nx} x {Ny}
Domain in x          = [0, {Lx*1e3:.2f}] mm
Domain in y          = [0, {Ly*1e3:.2f}] mm
Total grid pts       = {Nx*Ny}
Tstar                = {Tstar:.2f}
Estimated time steps = {NT_est}
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
    "m": m,
    "n": n,
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
    "low_Mach": low_Mach,
    "wave_speeds": 1,
    "viscous": "T",
    "fd_order": fd_order,
    # Boundary conditions:
    # characteristic subsonic inflow at x-beg,
    # outflow at x-end, and periodic boundaries in y
    "bc_x%beg": -7,
    "bc_x%end": -3,
    "bc_y%beg": -1,
    "bc_y%end": -1,
    # Prescribed state for the characteristic subsonic inflow
    "bc_x%grcbc_in": "T",
    "bc_x%vel_in(1)": u_inf,
    "bc_x%vel_in(2)": 0.0,
    "bc_x%pres_in": p_inf,
    "bc_x%alpha_rho_in(1)": rho_inf,
    "bc_x%alpha_in(1)": 1.0,
    # Chemistry and transport
    "chemistry": chemistry,
    "chem_params%diffusion": diffusion,
    "chem_params%reactions": reactions,
    "chem_params%transport_model": transport_model,
    "chem_params%reaction_substeps": substeps,
    "chem_params%adap_substeps": adap,
    "chem_params%reaction_substeps_max": substeps_max,
    "cantera_file": ctfile,
    "surface_cantera_file": surface_ctfile,
    "surface_phase": "carbon_surface",
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
    # Reacting carbon cylinder
    #   Offset the cylinder slightly from the grid symmetry line
    #   to seed the wake instability.
    "ib": "T",
    "num_ibs": 1,
    "patch_ib(1)%geometry": 2,
    "patch_ib(1)%x_centroid": x_cyl,
    "patch_ib(1)%y_centroid": y_cyl + dy / 3.0,
    "patch_ib(1)%radius": r_cyl,
    "patch_ib(1)%slip": "F",
    "patch_ib(1)%thermal_bc": thermal_bc,
    "patch_ib(1)%Twall": Twall,
    "patch_ib(1)%surface_reaction": surface_reaction,
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
    "ib_state_wrt": "T",
}


# Initialize all gas species in Cantera mechanism order and request
# each species mass fraction in the simulation output.
for i in range(len(Y_inf_vec)):
    case[f"patch_icpp(1)%Y({i + 1})"] = float(Y_inf_vec[i])
    case[f"chem_wrt_Y({i + 1})"] = "T"


if __name__ == "__main__":
    print(json.dumps(case))
