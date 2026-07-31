#!/usr/bin/env python3
"""One-step regression for the hypoelastic HLLD inner-fan fallback.

The stationary x-interface separates a normalized solid (rho=1, G=1) from a
light gas (rho=0.01, G=0) at equal pressure and zero normal velocity.  In the
solid-anchored pass, S_M=0 and C_NC~=1, so the raw right inner speed is about
sqrt(C_NC)/rho_R*=100 while the outer speed is about sqrt(1.4/0.01)=11.8.
The resulting five-wave fan is invalid and must use the HLL fallback.

The tangential velocity jump makes the HLL fallback observably different from
the former speed-only clamp.  The nearly pure gas anchor has C_NC below the
existing degeneracy floor and independently exercises the unchanged inner-wave
collapse in the other pass.
"""

import argparse
import json
import math

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--mfc", type=json.loads, default={})
parser.parse_args()

gamma = 1.4
p0 = 1.0
G_solid = 1.0
rho_solid = 1.0
rho_gas = 1.0e-2
eps_alpha = 1.0e-12

nx = 32
ny = 8
dx = 1.0 / nx
c_outer = math.sqrt(gamma * p0 / rho_gas)
dt = 0.2 * dx / c_outer

print(
    json.dumps(
        {
            "run_time_info": "F",
            "x_domain%beg": 0.0,
            "x_domain%end": 1.0,
            "y_domain%beg": 0.0,
            "y_domain%end": 1.0,
            "m": nx - 1,
            "n": ny - 1,
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": 1,
            "t_step_save": 1,
            "num_patches": 2,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 2,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 1,
            "recon_type": 1,
            "weno_order": 1,
            "weno_eps": 1.0e-16,
            "riemann_solver": 4,
            "wave_speeds": 1,
            "avg_state": 2,
            "riemann_hypo_ADC": "F",
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -1,
            "bc_y%end": -1,
            "format": 1,
            "precision": 2,
            "parallel_io": "F",
            "hypoelasticity": "T",
            "fd_order": 4,
            "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
            "fluid_pp(1)%pi_inf": 0.0,
            "fluid_pp(1)%G": G_solid,
            "fluid_pp(2)%gamma": 1.0 / (gamma - 1.0),
            "fluid_pp(2)%pi_inf": 0.0,
            "fluid_pp(2)%G": 0.0,
            # Solid background/left state.
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%y_centroid": 0.5,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%length_y": 1.0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": p0,
            "patch_icpp(1)%tau_e(1)": 0.0,
            "patch_icpp(1)%tau_e(2)": 0.0,
            "patch_icpp(1)%tau_e(3)": 0.0,
            "patch_icpp(1)%alpha_rho(1)": rho_solid * (1.0 - eps_alpha),
            "patch_icpp(1)%alpha(1)": 1.0 - eps_alpha,
            "patch_icpp(1)%alpha_rho(2)": rho_solid * eps_alpha,
            "patch_icpp(1)%alpha(2)": eps_alpha,
            # Light-gas right state with tangential slip.
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%x_centroid": 0.75,
            "patch_icpp(2)%y_centroid": 0.5,
            "patch_icpp(2)%length_x": 0.5,
            "patch_icpp(2)%length_y": 1.0,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 1.0,
            "patch_icpp(2)%pres": p0,
            "patch_icpp(2)%tau_e(1)": 0.0,
            "patch_icpp(2)%tau_e(2)": 0.0,
            "patch_icpp(2)%tau_e(3)": 0.0,
            "patch_icpp(2)%alpha_rho(1)": rho_gas * eps_alpha,
            "patch_icpp(2)%alpha(1)": eps_alpha,
            "patch_icpp(2)%alpha_rho(2)": rho_gas * (1.0 - eps_alpha),
            "patch_icpp(2)%alpha(2)": 1.0 - eps_alpha,
        }
    )
)
