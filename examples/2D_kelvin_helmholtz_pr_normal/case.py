#!/usr/bin/env python3
"""Periodic xenon shear-layer roll-up; SI inputs, matched nondimensional flow."""

import json
import math
import os

import numpy as np
from scipy.optimize import brentq

# Edit grid and duration here; environment overrides are useful for smoke tests.
NX = int(os.getenv("MFC_KH_NX", "128"))
NY = int(os.getenv("MFC_KH_NY", str(NX)))
FLOW_TIMES = float(os.getenv("MFC_KH_FLOW_TIMES", "2.0"))
SAVES = 100
PARALLEL_IO = os.getenv("MFC_KH_PARALLEL_IO", "T")
LENGTH = 0.01  # m
SHEAR_WIDTH = 0.025  # fraction of LENGTH
PERTURBATION = 0.1  # fraction of reference velocity
R = 8.31446261815324 / 0.131293
TC, PC, OMEGA = 289.7325683307143, 5841913.640364157, 0.00363
CV = 1.5 * R  # monatomic ideal-gas limit, gamma = 5/3
IDEAL = False
SUPERCRITICAL = False
PRESSURE = 101325.0
COLD_TEMPERATURE = 300.0  # dense-stream T: ideal at 1 atm, PR at supercritical pressure
B = 0.07780 * R * TC / PC
AC = 0.45724 * (R * TC) ** 2 / PC
KAPPA = 0.37464 + 1.54226 * OMEGA - 0.26992 * OMEGA**2
C, D = 1 + KAPPA, KAPPA / math.sqrt(TC)


def departure(rho):
    return np.log((1 + (1 + np.sqrt(2)) * B * rho) / (1 + (1 - np.sqrt(2)) * B * rho)) / (2 * np.sqrt(2) * B)


def properties(rho, temperature, ideal=False):
    """Pressure, energy, cv, cp and acoustic speed from the selected closure."""
    if ideal:
        return rho * R * temperature, CV * temperature, CV, CV + R, np.sqrt((1 + R / CV) * R * temperature)
    factor = C - D * np.sqrt(temperature)
    a = AC * factor**2
    da = -AC * factor * D / np.sqrt(temperature)
    d2a = AC * C * D / (2 * temperature**1.5)
    den = 1 + 2 * B * rho - (B * rho) ** 2
    integral = departure(rho)
    pressure = rho * R * temperature / (1 - B * rho) - a * rho**2 / den
    energy = CV * temperature + (temperature * da - a) * integral
    cv = CV + temperature * d2a * integral
    pt = rho * R / (1 - B * rho) - da * rho**2 / den
    pr = R * temperature / (1 - B * rho) ** 2 - 2 * a * rho * (1 + B * rho) / den**2
    cp = cv + temperature * pt**2 / (rho**2 * pr)
    sound = np.sqrt(pr + temperature * pt**2 / (rho**2 * cv))
    return pressure, energy, cv, cp, sound


def temperature_from_energy(rho, energy, ideal=False):
    if ideal:
        return energy / CV
    integral = departure(rho)
    linear = AC * C * D * integral
    constant = energy + AC * C**2 * integral
    root = 2 * constant / (linear + np.sqrt(linear**2 + 4 * CV * constant))
    return root**2


def temperature_from_pressure(rho, pressure):
    if IDEAL:
        return pressure / (rho * R)
    den = 1 + 2 * B * rho - (B * rho) ** 2
    qa = rho * R / (1 - B * rho) - AC * D**2 * rho**2 / den
    qb = 2 * AC * C * D * rho**2 / den
    qc = -pressure - AC * C**2 * rho**2 / den
    return (-2 * qc / (qb + np.sqrt(qb**2 - 4 * qa * qc))) ** 2


def reference():
    # The dilute pair shares rho, p and u exactly, not temperature or Mach.
    rho = PRESSURE / (2 * R * COLD_TEMPERATURE)
    if SUPERCRITICAL:
        rho = 0.5 * brentq(lambda r: properties(r, COLD_TEMPERATURE)[0] - PRESSURE, 1e-6, 0.9 / B)
    velocity = math.sqrt(PRESSURE / (2.5 * rho))
    return dict(rho=rho, velocity=velocity, length=LENGTH, time=LENGTH / velocity)


def case():
    s = reference()
    steps = SAVES * math.ceil(FLOW_TIMES * max(NX, NY) / 0.025 / SAVES)
    layer = f"(0.5*(1-tanh(cos(2*pi*y/ly)/(2*pi*{SHEAR_WIDTH}))))"
    gaussian = "(exp(-((y/ly-0.25)/0.05)**2)+exp(-((y/ly-0.75)/0.05)**2))"
    result = {
        "run_time_info": "F",
        "x_domain%beg": 0.0,
        "x_domain%end": LENGTH,
        "y_domain%beg": 0.0,
        "y_domain%end": LENGTH,
        "m": NX*3 - 1 +2,
        "n": NY*3 - 1 +2,
        "p": 0,
        # Identical nondimensional save times, conservative acoustic time step.
        "dt": FLOW_TIMES * s["time"] / steps,
        "t_step_start": 0,
        "t_step_stop": steps,
        "t_step_save": steps // SAVES,
        "num_patches": 1,
        "model_eqns": "5eq",
        "num_fluids": 1,
        "time_stepper": "rk3",
        "weno_order": 5,
        "weno_eps": 1e-10,
        "null_weights": "F",
        "mp_weno": "F",
        "riemann_solver": "hllc",
        "wave_speeds": "direct",
        "avg_state": "arithmetic",
        "bc_x%beg": -1,
        "bc_x%end": -1,
        "bc_y%beg": -1,
        "bc_y%end": -1,
        "format": "silo",
        "precision": "double",
        "parallel_io": PARALLEL_IO,
        "prim_vars_wrt": "T",
        "T_wrt": "T",
        "patch_icpp(1)%geometry": 3,
        "patch_icpp(1)%x_centroid": LENGTH / 2,
        "patch_icpp(1)%y_centroid": LENGTH / 2,
        "patch_icpp(1)%length_x": LENGTH,
        "patch_icpp(1)%length_y": LENGTH,
        "patch_icpp(1)%alpha(1)": 1.0,
        "patch_icpp(1)%alpha_rho(1)": f"{s['rho']}*(1+{layer})",
        "patch_icpp(1)%vel(1)": f"{s['velocity']}*({layer}-0.5)",
        "patch_icpp(1)%vel(2)": f"{PERTURBATION*s['velocity']}*sin(4*pi*x/lx)*{gaussian}",
        "patch_icpp(1)%pres": PRESSURE,
        "fluid_pp(1)%cv": CV,
    }
    if IDEAL:
        result.update({"fluid_pp(1)%eos": "ideal_gas", "fluid_pp(1)%gamma": CV / R})
    else:
        result.update({"fluid_pp(1)%eos": "peng_robinson", "fluid_pp(1)%pr_tc": TC, "fluid_pp(1)%pr_pc": PC, "fluid_pp(1)%pr_r": R, "fluid_pp(1)%pr_omega": OMEGA})
    return result


if __name__ == "__main__":
    print(json.dumps(case()))
