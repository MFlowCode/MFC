#!/usr/bin/env python3
"""Matched-Reynolds/Mach isothermal nitrogen flat plates (SI units)."""

import json
import math

import numpy as np

R = 296.803052
TC, PC, OMEGA = 126.192, 3.3958e6, 0.0372
CV = 2.5 * R
MU, PRANDTL, MACH, REYNOLDS = 1.8e-5, 0.72, 0.2, 10000.0
NX, NY = 249, 169
# Edit these settings independently in this case file.
SAVES = 100  # Set to e.g. 80 for 81 saved animation frames, including the initial state.
FLOW_TIMES = 2.0
PARALLEL_IO = "T"  # Use --mpi; set to "F" for a --no-mpi run.
IDEAL = True
TEMPERATURE = 300.0  # K
PRESSURE = 101325.0  # Pa
WALL_TEMPERATURE = 270.0  # K
STRETCH, CORE = 4.0, 0.25


def properties(rho, temperature, ideal=True):
    """Calorically perfect nitrogen: pressure, energy, cv, cp, sound speed."""
    return rho * R * temperature, CV * temperature, CV, CV + R, np.sqrt(1.4 * R * temperature)


def temperature_from_energy(rho, energy, ideal=True):
    return energy / CV


def reference():
    ideal = IDEAL
    temperature, pressure = TEMPERATURE, PRESSURE
    rho = pressure / (R * temperature)
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
    result.update({"fluid_pp(1)%eos": "ideal_gas", "fluid_pp(1)%gamma": CV / R})
    return result


if __name__ == "__main__":
    print(json.dumps(case()))
