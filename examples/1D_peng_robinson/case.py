#!/usr/bin/env python3
"""Small periodic acoustic wave in single-phase supercritical nitrogen, SI units."""

import json
import math

R = 296.803052
TC = 126.192
PC = 3.3958e6
OMEGA = 0.0372
CV = 2.5 * R
RHO = 150.0
TEMP = 200.0
LENGTH = 1.0
AMPLITUDE = 1.0e-5
CELLS = 64
STEPS = 400


def reference_state():
    """PR pressure and isentropic acoustic speed at the uniform base state."""
    b = 0.07780 * R * TC / PC
    ac = 0.45724 * (R * TC) ** 2 / PC
    kappa = 0.37464 + 1.54226 * OMEGA - 0.26992 * OMEGA**2
    factor = 1 + kappa * (1 - math.sqrt(TEMP / TC))
    a = ac * factor**2
    da = -ac * factor * kappa / math.sqrt(TEMP * TC)
    d2a = ac * kappa * (1 + kappa) / (2 * math.sqrt(TC) * TEMP**1.5)
    denominator = 1 + 2 * b * RHO - (b * RHO) ** 2
    integral = math.log((1 + (1 + math.sqrt(2)) * b * RHO) / (1 + (1 - math.sqrt(2)) * b * RHO)) / (2 * math.sqrt(2) * b)
    pressure = RHO * R * TEMP / (1 - b * RHO) - a * RHO**2 / denominator
    cv = CV + TEMP * d2a * integral
    dpdt = RHO * R / (1 - b * RHO) - da * RHO**2 / denominator
    dpdrho = R * TEMP / (1 - b * RHO) ** 2 - 2 * a * RHO * (1 + b * RHO) / denominator**2
    sound = math.sqrt(dpdrho + TEMP * dpdt**2 / (RHO**2 * cv))
    return pressure, sound


def case():
    pressure, sound = reference_state()
    sine = "sin(2*pi*x)"
    return {
        "run_time_info": "F",
        "x_domain%beg": 0.0,
        "x_domain%end": LENGTH,
        "m": CELLS - 1,
        "n": 0,
        "p": 0,
        "dt": LENGTH / sound / STEPS,
        "t_step_start": 0,
        "t_step_stop": STEPS,
        "t_step_save": STEPS // 4,
        "num_patches": 1,
        "model_eqns": "5eq",
        "num_fluids": 1,
        "time_stepper": "rk3",
        "weno_order": 5,
        "weno_eps": 1e-16,
        "mapped_weno": "T",
        "riemann_solver": "hllc",
        "wave_speeds": "direct",
        "avg_state": "arithmetic",
        "bc_x%beg": -1,
        "bc_x%end": -1,
        "format": "silo",
        "precision": "double",
        "parallel_io": "F",
        "prim_vars_wrt": "T",
        "cons_vars_wrt": "T",
        "T_wrt": "T",
        "patch_icpp(1)%geometry": 1,
        "patch_icpp(1)%x_centroid": LENGTH / 2,
        "patch_icpp(1)%length_x": LENGTH,
        "patch_icpp(1)%alpha(1)": 1.0,
        "patch_icpp(1)%alpha_rho(1)": f"{RHO}*(1+{AMPLITUDE}*{sine})",
        "patch_icpp(1)%vel(1)": f"{sound * AMPLITUDE}*{sine}",
        "patch_icpp(1)%pres": f"{pressure}+{sound**2 * RHO * AMPLITUDE}*{sine}",
        "fluid_pp(1)%eos": "peng_robinson",
        "fluid_pp(1)%pr_tc": TC,
        "fluid_pp(1)%pr_pc": PC,
        "fluid_pp(1)%pr_omega": OMEGA,
        "fluid_pp(1)%pr_r": R,
        "fluid_pp(1)%cv": CV,
    }


if __name__ == "__main__":
    print(json.dumps(case()))
