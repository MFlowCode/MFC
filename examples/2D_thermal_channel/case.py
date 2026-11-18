#!/usr/bin/env python3
"""
2D plane channel between isothermal no-slip walls (bc_y = -16) at Tw, periodic in x, driven at a bulk mass flux
(mean_mflux). The thermostat with mean_temp_wall = Tw holds the gas in the thermally fully developed state
T = Tw + psi(y) exp(-beta x), with the mixing-cup temperature at T0. Constant-property reference (Shah & London 1978):
Nu = q_w D_h / (k (Tw - T_m)) = 7.541, D_h = 2H, T_m the mixing-cup temperature; axial conduction changes it by
under 0.02 % for Pe_Dh > 100. q_w = -mean_qdot V_f / (2 L_x) at steady state (mean_temp.dat). A hot/cold wall pair
(--tw above and below --t0) cancels the viscous heating to first order.
"""

import argparse
import json
import math

import cantera as ct

parser = argparse.ArgumentParser()
parser.add_argument("--ny", type=int, default=32, help="cells across the gap H")
parser.add_argument("--nx", type=int, default=25, help="cells along x (WENO5 needs at least 25 per direction)")
parser.add_argument("--tw", type=float, default=330.0, help="wall temperature [K]")
parser.add_argument("--t0", type=float, default=300.0, help="initial gas and target mixing-cup temperature [K]")
parser.add_argument("--mach", type=float, default=0.1, help="bulk Mach number")
parser.add_argument("--re", type=float, default=150.0, help="Re_Dh = rho U_b 2H / mu")
parser.add_argument("--tend", type=float, default=0.5, help="run time [H^2 / alpha]; developed in about 0.1")
parser.add_argument("--saves", type=int, default=10)
args, _ = parser.parse_known_args()

gas = ct.Solution("h2o2.yaml")
gas.TPY = args.t0, ct.one_atm, "O2:0.233,N2:0.767"
rho, mu, k = gas.density, gas.viscosity, gas.thermal_conductivity
a = math.sqrt(gas.cp_mass / gas.cv_mass * ct.one_atm / rho)
U = args.mach * a
H = args.re * mu / (2.0 * rho * U)
dy = H / args.ny
Lx = args.nx * dy
dt = 0.4 * dy / (U + a)
nt = args.saves * math.ceil(args.tend * H**2 * rho * gas.cp_mass / k / dt / args.saves)

case = {
    "run_time_info": "F",
    "x_domain%beg": 0.0,
    "x_domain%end": Lx,
    "y_domain%beg": 0.0,
    "y_domain%end": H,
    "m": args.nx - 1,
    "n": args.ny - 1,
    "p": 0,
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": nt,
    "t_step_save": nt // args.saves,
    "num_patches": 1,
    "model_eqns": 2,
    "num_fluids": 1,
    "mpp_lim": "F",
    "mixture_err": "T",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1e-16,
    "mapped_weno": "T",
    "weno_Re_flux": "T",
    "mp_weno": "F",
    "weno_avg": "T",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    "fd_order": 4,
    "bc_x%beg": -1,
    "bc_x%end": -1,
    "bc_y%beg": -16,
    "bc_y%end": -16,
    "bc_y%isothermal_in": "T",
    "bc_y%isothermal_out": "T",
    "bc_y%Twall_in": args.tw,
    "bc_y%Twall_out": args.tw,
    "viscous": "T",
    "chemistry": "T",
    "chem_params%diffusion": "T",
    "chem_params%reactions": "F",
    "chem_params%transport_model": 1,
    "cantera_file": "h2o2.yaml",
    "chem_wrt_T": "T",
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "parallel_io": "F",
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.5 * Lx,
    "patch_icpp(1)%y_centroid": 0.5 * H,
    "patch_icpp(1)%length_x": Lx,
    "patch_icpp(1)%length_y": H,
    "patch_icpp(1)%vel(1)": U,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": ct.one_atm,
    "patch_icpp(1)%alpha_rho(1)": rho,
    "patch_icpp(1)%alpha(1)": 1.0,
    "fluid_pp(1)%gamma": 1.0 / (1.4 - 1.0),
    "fluid_pp(1)%eos": "stiffened_gas",
    "fluid_pp(1)%pi_inf": 0.0,
    "fluid_pp(1)%Re(1)": 1.0 / mu,
    # Poiseuille drive and the developed wall loss (Nu 7.541) as initial controller states
    "mean_dpdx(1)": -12.0 * mu * U / H**2,
    "mean_mflux(1)": rho * U,
    "mean_dpdx_gain": 0.05,
    "mean_qdot": -7.541 * k * (args.tw - args.t0) / H**2,
    "mean_temp": args.t0,
    "mean_temp_gain": 0.05,
    "mean_temp_wall": args.tw,
}
for i, y in enumerate(gas.Y):
    case[f"patch_icpp(1)%Y({i + 1})"] = y
    case[f"chem_wrt_Y({i + 1})"] = "T"

print(json.dumps(case))
