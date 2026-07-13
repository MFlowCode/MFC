#!/usr/bin/env python3
# 2D immersed-boundary flame-holder geometry: a solid cylinder in a hot,
# near-ignition H2/O2/Ar premix. IBM + chemistry with temperature-dependent
# transport used to hang/diverge; that is now fixed in m_ibm.fpp (species and a
# consistent energy are set at immersed-boundary points). Chemistry + diffusion
# are on here (species transport around the cylinder + wake). Reactions are off
# to avoid a separate explicit-chemistry ignition-stiffness issue; enable
# chem_params%reactions with a smaller dt to attempt a wake-anchored flame.
import json

import cantera as ct

ctfile = "h2o2.yaml"
X = "H2:2,O2:1,AR:5"
T0, P0 = 1000.0, 101325.0  # near-ignition inflow temperature, 1 atm
u0 = 40.0  # inflow velocity [m/s]

gas = ct.Solution(ctfile)
gas.TPX = T0, P0, X
rho0 = gas.density
mu0 = gas.viscosity
c0 = gas.sound_speed

# Domain: 8 cm x 3 cm channel, cylinder near the inlet.
Lx, Ly = 0.08, 0.03
m, n = 200, 75
dx = Lx / m

r_cyl = 0.004
x_cyl = 0.016
y_cyl = Ly / 2.0

# dt from the acoustic CFL (H2-rich mixture -> high sound speed, low Mach).
cfl = 0.05
dt = cfl * dx / (u0 + c0)
tend = 4.0e-4  # a few cylinder flow-through times (D/u0 ~ 0.2 ms)
NT = int(tend / dt)
NS = max(1, NT // 4)

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
    "dt": float(dt),
    "t_step_start": 0,
    "t_step_stop": NT,
    "t_step_save": NS,
    "t_step_print": NS,
    "parallel_io": "T",
    # Algorithm
    "model_eqns": 2,
    "alt_soundspeed": "F",
    "num_fluids": 1,
    "num_patches": 1,
    "mpp_lim": "F",
    "mixture_err": "T",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1e-16,
    "mapped_weno": "T",
    "mp_weno": "T",
    "weno_avg": "T",
    "weno_Re_flux": "T",
    "null_weights": "F",
    "riemann_solver": "hllc",
    "wave_speeds": "direct",
    "avg_state": "arithmetic",
    "fd_order": 2,
    "viscous": "T",
    # BCs: characteristic inflow, extrapolation outflow, slip channel walls
    "bc_x%beg": -7,
    "bc_x%end": -3,
    "bc_y%beg": -15,
    "bc_y%end": -15,
    # Chemistry on (species transport); reactions off (see header note)
    "chemistry": "T",
    "chem_params%diffusion": "T",
    "chem_params%reactions": "F",
    "cantera_file": ctfile,
    "chem_wrt_T": "T",
    # Immersed boundary: inert solid cylinder (flame-holder)
    "ib": "T",
    "num_ibs": 1,
    "patch_ib(1)%geometry": 2,
    "patch_ib(1)%x_centroid": x_cyl,
    "patch_ib(1)%y_centroid": y_cyl,
    "patch_ib(1)%radius": r_cyl,
    "patch_ib(1)%slip": "F",
    # Output
    "format": "silo",
    "precision": "double",
    "prim_vars_wrt": "T",
    "ib_state_wrt": "T",
    # Patch: hot H2/O2/Ar premix fills the domain
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": Lx / 2,
    "patch_icpp(1)%y_centroid": Ly / 2,
    "patch_icpp(1)%length_x": Lx,
    "patch_icpp(1)%length_y": Ly,
    "patch_icpp(1)%vel(1)": u0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": P0,
    "patch_icpp(1)%alpha_rho(1)": rho0,
    "patch_icpp(1)%alpha(1)": 1.0,
    # Fluid EOS (calorically perfect gas matched to the mixture at T0, P0)
    "fluid_pp(1)%gamma": 1.0 / (1.4 - 1.0),
    "fluid_pp(1)%pi_inf": 0.0,
    "fluid_pp(1)%Re(1)": 1.0 / mu0,
}

for i in range(len(gas.Y)):
    case[f"patch_icpp(1)%Y({i + 1})"] = float(gas.Y[i])

if __name__ == "__main__":
    print(json.dumps(case))
