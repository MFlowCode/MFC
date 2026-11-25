#!/usr/bin/env python3
"""Periodic advection of an isobaric xenon density interface in SI units."""

import json
import math

from scipy.optimize import brentq

CELLS = 128
WIDTH = 0.005  # Interface width divided by domain length
COLD_T = 330.0
PRESSURE = 10e6
SPEED = 20.0
LENGTH = 0.01
R = 8.31446261815324 / 0.131293
TC, PC, OMEGA = 289.7325683307143, 5841913.640364157, 0.00363
CV = 1.5 * R
STEPS = 100 * math.ceil(64 * CELLS / 100)  # One transit with 100 equal save intervals


def light_stream_density():
    b = 0.07780 * R * TC / PC
    kappa = 0.37464 + 1.54226 * OMEGA - 0.26992 * OMEGA**2
    a = 0.45724 * (R * TC) ** 2 / PC * (1 + kappa * (1 - math.sqrt(COLD_T / TC))) ** 2

    def residual(rho):
        return rho * R * COLD_T / (1 - b * rho) - a * rho**2 / (1 + 2 * b * rho - (b * rho) ** 2) - PRESSURE

    return 0.5 * brentq(residual, 1e-6, 0.95 / b)


def case():
    light_density = light_stream_density()
    return {
        "run_time_info": "F",
        "x_domain%beg": 0.0,
        "x_domain%end": LENGTH,
        "m": CELLS - 1,
        "n": 0,
        "p": 0,
        "dt": LENGTH / SPEED / STEPS,
        "t_step_start": 0,
        "t_step_stop": STEPS,
        "t_step_save": STEPS // 100,
        "num_patches": 1,
        "model_eqns": "5eq",
        "num_fluids": 1,
        "time_stepper": "rk3",
        "weno_order": 5,
        "weno_eps": 1e-10,
        "mapped_weno": "F",
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
        # Runtime parameters for the analytic density profile.
        "patch_icpp(1)%beta": light_density,
        "patch_icpp(1)%epsilon": WIDTH,
        "patch_icpp(1)%alpha_rho(1)": "beta*(1.5-0.5*tanh(cos(2*pi*x/lx)/(2*pi*eps)))",
        "patch_icpp(1)%vel(1)": SPEED,
        "patch_icpp(1)%pres": PRESSURE,
        "fluid_pp(1)%eos": "peng_robinson",
        "fluid_pp(1)%pr_tc": TC,
        "fluid_pp(1)%pr_pc": PC,
        "fluid_pp(1)%pr_omega": OMEGA,
        "fluid_pp(1)%pr_r": R,
        "fluid_pp(1)%cv": CV,
    }


if __name__ == "__main__":
    print(json.dumps(case()))
