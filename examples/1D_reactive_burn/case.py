#!/usr/bin/env python3
# 1D condensed-phase detonation via programmed pressure burn (reactive_burn).
# A dense stiffened-gas "reactant" is converted to a "product" fluid that shares
# the same mechanical EOS (gamma, pi_inf) but has a lower reference energy qv --
# so reactant -> product releases (qv_reactant - qv_product) per unit mass through
# the mixture EOS, with no explicit energy source. A high-pressure initiator at the
# left end launches a shock; where the shock raises the pressure above rburn_pign
# the reactant burns, and the energy release sustains a self-propagating detonation
# (von Neumann spike + Taylor rarefaction, near-CJ speed). Multi-fluid model,
# chemistry OFF -- this deliberately bypasses the Cantera num_fluids=1 lock.
import argparse
import json

parser = argparse.ArgumentParser(prog="1D_reactive_burn")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC toolchain state.")
parser.add_argument("--res", type=float, default=1.0, help="Resolution multiplier: m=3000*res.")
parser.add_argument("--tend", type=float, default=8.0e-6, help="Physical end time [s].")
parser.add_argument("--frames", type=int, default=80, help="Number of saved output frames.")
args = parser.parse_args()

# --- Stiffened-gas EOS shared by reactant and product (mechanically identical) ---
Gamma = 3.0  # physical adiabatic exponent
Pi = 6.0e8  # stiffening pressure [Pa]
gamma_p = 1.0 / (Gamma - 1.0)  # MFC gamma parameter
pi_inf_p = Gamma * Pi / (Gamma - 1.0)  # MFC pi_inf parameter
Q = 4.0e6  # heat of reaction [J/kg] (reactant qv; product qv = 0)

rho0 = 1600.0  # reactant density [kg/m^3]
p0 = 1.0e5  # ambient pressure [Pa]
c0 = (Gamma * (p0 + Pi) / rho0) ** 0.5  # reactant sound speed

# Domain
L = 0.06
m = int(3000 * args.res)
dx = L / m
x_init = 0.004  # initiator occupies the left 4 mm
p_init = 5.0e9  # initiator pressure [Pa] -> launches an igniting shock

# CJ-scale detonation speed sets dt.
D_guess = 8000.0
cfl = 0.3
dt = cfl * dx / (D_guess + c0)
NT = int(args.tend / dt)
NS = max(1, NT // args.frames)

case = {
    "run_time_info": "T",
    "x_domain%beg": 0.0,
    "x_domain%end": L,
    "m": m,
    "n": 0,
    "p": 0,
    "dt": float(dt),
    "t_step_start": 0,
    "t_step_stop": NT,
    "t_step_save": NS,
    "t_step_print": max(1, NT // 20),
    "parallel_io": "T",
    # Algorithm
    "model_eqns": 2,
    "num_fluids": 2,
    "num_patches": 2,
    "mpp_lim": "T",
    "mixture_err": "T",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1e-16,
    "mapped_weno": "T",
    "mp_weno": "T",
    "weno_avg": "F",
    "riemann_solver": "hllc",
    "wave_speeds": "direct",
    "avg_state": "arithmetic",
    # Reactant reflects at the initiator wall; open at the far end.
    "bc_x%beg": -2,
    "bc_x%end": -3,
    # Condensed-phase reactive burn
    "reactive_burn": "T",
    "rburn_k": 5.0e6,  # rate coefficient [1/s]
    "rburn_pign": 5.0e8,  # ignition pressure threshold [Pa]
    "rburn_pref": 1.0e9,  # reference pressure for the pressure drive [Pa]
    "rburn_n": 1.0,  # pressure exponent
    # Output
    "format": "silo",
    "precision": "double",
    "prim_vars_wrt": "T",
    # Patch 1: unreacted reactant fills the domain (fluid 1)
    "patch_icpp(1)%geometry": 1,
    "patch_icpp(1)%x_centroid": L / 2,
    "patch_icpp(1)%length_x": L,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%pres": p0,
    "patch_icpp(1)%alpha_rho(1)": rho0,
    "patch_icpp(1)%alpha_rho(2)": 0.0,
    "patch_icpp(1)%alpha(1)": 1.0,
    "patch_icpp(1)%alpha(2)": 0.0,
    # Patch 2: high-pressure initiator (reactant), left end
    "patch_icpp(2)%geometry": 1,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%x_centroid": x_init / 2,
    "patch_icpp(2)%length_x": x_init,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%pres": p_init,
    "patch_icpp(2)%alpha_rho(1)": rho0,
    "patch_icpp(2)%alpha_rho(2)": 0.0,
    "patch_icpp(2)%alpha(1)": 1.0,
    "patch_icpp(2)%alpha(2)": 0.0,
    # Fluid EOS: reactant (1) and product (2) share gamma/pi_inf, differ only in qv.
    "fluid_pp(1)%gamma": gamma_p,
    "fluid_pp(1)%pi_inf": pi_inf_p,
    "fluid_pp(1)%qv": Q,
    "fluid_pp(2)%gamma": gamma_p,
    "fluid_pp(2)%pi_inf": pi_inf_p,
    "fluid_pp(2)%qv": 0.0,
}

if __name__ == "__main__":
    print(json.dumps(case))
