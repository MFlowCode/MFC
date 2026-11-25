#!/usr/bin/env python3
"""Matched-Reynolds/Mach isothermal nitrogen flat plates (SI units)."""

import json
import math

import numpy as np
from scipy.optimize import brentq

R = 296.803052
TC, PC, OMEGA = 126.192, 3.3958e6, 0.0372
CV = 2.5 * R
MU, PRANDTL, MACH, REYNOLDS = 1.8e-5, 0.72, 0.2, 10000.0
NX, NY = 249, 169
# Edit these settings independently in this case file.
SAVES = 100  # Set to e.g. 80 for 81 saved animation frames, including the initial state.
FLOW_TIMES = 2.0
PARALLEL_IO = "T"  # Use --mpi; set to "F" for a --no-mpi run.
IDEAL = False
TEMPERATURE = 300.0  # K
PRESSURE = 7500000.0  # Pa
WALL_TEMPERATURE = 270.0  # K
STRETCH, CORE = 4.0, 0.25
B = 0.07780 * R * TC / PC
AC = 0.45724 * (R * TC) ** 2 / PC
KAPPA = 0.37464 + 1.54226 * OMEGA - 0.26992 * OMEGA**2
C, D = 1 + KAPPA, KAPPA / math.sqrt(TC)


def departure(rho):
    return np.log((1 + (1 + np.sqrt(2)) * B * rho) / (1 + (1 - np.sqrt(2)) * B * rho)) / (2 * np.sqrt(2) * B)


def properties(rho, temperature, ideal=False):
    """Pressure, energy, cv, cp and acoustic speed from the selected closure."""
    if ideal:
        return rho * R * temperature, CV * temperature, CV, CV + R, np.sqrt(1.4 * R * temperature)
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


def reference():
    ideal = IDEAL
    temperature, pressure = TEMPERATURE, PRESSURE
    rho = pressure / (R * temperature) if ideal else brentq(lambda r: properties(r, temperature)[0] - pressure, 1e-6, 0.9 / B)
    _, _, _, cp, sound = properties(rho, temperature, ideal)
    velocity = MACH * sound
    length = REYNOLDS * MU / (rho * velocity)
    grid_height = 0.2 * length
    height = grid_height * (1 + (math.log(math.cosh(STRETCH * (1 + CORE))) + math.log(math.cosh(STRETCH * (1 - CORE))) - 2 * math.log(math.cosh(STRETCH * CORE))) / STRETCH)
    spacing = min(length / NX, grid_height / NY)
    conductivity = MU * cp / PRANDTL
    # Sum of directional acoustic rates and conservative ideal-cv diffusion bound.
    dt = min(0.4 / ((velocity + sound) / (length / NX) + sound / (grid_height / NY)), 0.2 * spacing**2 * rho / max(MU, conductivity / CV))
    steps = SAVES * math.ceil(FLOW_TIMES * length / velocity / dt / SAVES)
    return dict(
        ideal=ideal,
        temperature=temperature,
        pressure=pressure,
        rho=rho,
        cp=cp,
        sound=sound,
        velocity=velocity,
        length=length,
        height=height,
        grid_height=grid_height,
        conductivity=conductivity,
        dt=dt,
        steps=steps,
        save_steps=steps // SAVES,
        wall_temperature=WALL_TEMPERATURE,
    )


def case():
    s = reference()
    result = {
        "run_time_info": "F",
        "x_domain%beg": 0.0,
        "x_domain%end": s["length"],
        "y_domain%beg": 0.0,
        "y_domain%end": s["grid_height"],
        "stretch_y": "T",
        "a_y": STRETCH,
        "y_a": -CORE * s["grid_height"],
        "y_b": CORE * s["grid_height"],
        "m": NX - 1,
        "n": NY - 1,
        "p": 0,
        "dt": s["dt"],
        "t_step_start": 0,
        "t_step_stop": s["steps"],
        "t_step_save": s["save_steps"],
        "num_patches": 1,
        "model_eqns": "5eq",
        "num_fluids": 1,
        "time_stepper": "rk3",
        "weno_order": 5,
        "weno_eps": 1e-16,
        "mapped_weno": "T",
        "weno_Re_flux": "T",
        "weno_avg": "T",
        "riemann_solver": "hllc",
        "wave_speeds": "direct",
        "avg_state": "arithmetic",
        "fd_order": 4,
        "bc_x%beg": -17,
        "bc_x%end": -3,
        "bc_y%beg": -16,
        "bc_y%end": -17,
        "bc_y%isothermal_in": "T",
        "bc_y%Twall_in": s["wall_temperature"],
        "viscous": "T",
        "format": "silo",
        "precision": "double",
        "parallel_io": PARALLEL_IO,
        "prim_vars_wrt": "T",
        "T_wrt": "T",
        "patch_icpp(1)%geometry": 3,
        "patch_icpp(1)%x_centroid": s["length"] / 2,
        "patch_icpp(1)%y_centroid": s["height"] / 2,
        "patch_icpp(1)%length_x": s["length"],
        "patch_icpp(1)%length_y": s["height"],
        "patch_icpp(1)%alpha(1)": 1.0,
        "patch_icpp(1)%alpha_rho(1)": s["rho"],
        "patch_icpp(1)%vel(1)": s["velocity"],
        "patch_icpp(1)%vel(2)": 0.0,
        "patch_icpp(1)%pres": s["pressure"],
        "fluid_pp(1)%cv": CV,
        "fluid_pp(1)%Re(1)": 1 / MU,
        "fluid_pp(1)%k_therm": s["conductivity"],
    }
    result.update({"fluid_pp(1)%eos": "peng_robinson", "fluid_pp(1)%pr_tc": TC, "fluid_pp(1)%pr_pc": PC, "fluid_pp(1)%pr_r": R, "fluid_pp(1)%pr_omega": OMEGA})
    return result


if __name__ == "__main__":
    print(json.dumps(case()))
