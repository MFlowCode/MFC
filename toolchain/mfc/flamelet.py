"""Cantera-only initial profiles shared by MFC's reacting mixing-layer examples.

Hot profiles use a unity-Lewis counterflow diffusion flame tabulated against
Bilger mixture fraction, then mapped onto the prescribed tanh mixing layer.
This replaces the former self-consistent scalar-dissipation flamelet solve;
the counterflow's nominal strain rate is an explicit initialization parameter.
"""

import json
import os
from dataclasses import dataclass

import cantera as ct
import numpy as np

from .thermochem.fingerprint import mechanism_fingerprint  # noqa: F401

INITIALIZER_VERSION = "cantera-counterflow-v1"


@dataclass
class SimulationFields:
    mixture_fraction: np.ndarray
    temperature: np.ndarray
    pressure: np.ndarray
    velocity: np.ndarray
    mass_fractions: np.ndarray


def stoichiometric_mixture_fraction(sol, y_ox, y_fu):
    return 1.0 / (1.0 + sol.stoich_air_fuel_ratio(y_fu, y_ox, basis="mass"))


def density(sol, pressure, temperature, mass_fractions):
    inv_weight = np.sum(mass_fractions / sol.molecular_weights[:, None], axis=0)
    return pressure / (ct.gas_constant * temperature * inv_weight)


def _hot_mass_fractions(sol, pres, temp_ox, temp_fu, y_ox, y_fu, width, strain_rate, mixture_fraction):
    if not np.isfinite(strain_rate) or strain_rate <= 0:
        raise ValueError("Flame initialization requires a positive, finite strain rate")
    original_state = sol.TPY
    original_transport = sol.transport_model
    try:
        flame = ct.CounterflowDiffusionFlame(sol, width=width)
        flame.P = pres
        flame.fuel_inlet.T = temp_fu
        flame.fuel_inlet.Y = y_fu
        flame.oxidizer_inlet.T = temp_ox
        flame.oxidizer_inlet.Y = y_ox
        inlet_speed = strain_rate * width / 2
        sol.TPY = temp_fu, pres, y_fu
        flame.fuel_inlet.mdot = sol.density * inlet_speed
        sol.TPY = temp_ox, pres, y_ox
        flame.oxidizer_inlet.mdot = sol.density * inlet_speed
        flame.transport_model = "unity-Lewis-number"
        flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.2)
        flame.solve(loglevel=0, auto=True)
        if flame.extinct():
            raise ValueError("Hot initialization converged to an extinguished flame; reduce flame_strain_rate")
        # Normalize against the inlet streams, not the flame-domain extrema:
        # diffusion can make the boundary composition differ from the inlet.
        z = np.empty(len(flame.grid))
        flame_y = flame.Y
        for i, y in enumerate(flame_y.T):
            sol.Y = y
            z[i] = sol.mixture_fraction(y_fu, y_ox, basis="mass")
        if np.any(np.diff(z) > 1e-7):
            raise ValueError("Counterflow mixture fraction is not monotone; cannot map the flame profile")
        # Supply exact stream endpoints and discard roundoff outside their range.
        interior = (z > 1e-12) & (z < 1 - 1e-12)
        z = np.concatenate(([0.0], z[interior], [1.0]))
        y = np.column_stack((y_ox, flame_y[:, interior], y_fu))
        z, indices = np.unique(z, return_index=True)
        y = y[:, indices]
        mapped = np.array([np.interp(mixture_fraction, z, species) for species in y])
        if np.min(mapped) < -1e-8:
            raise ValueError("Counterflow solve produced negative mass fractions")
        mapped = np.maximum(mapped, 0)
        return mapped / mapped.sum(axis=0)
    finally:
        sol.transport_model = original_transport
        sol.TPY = original_state


def create_simulation_fields(sol, pres, temp_ox, temp_fu, cross_coord, vort_thickness, stream_ox, stream_fu, strain_rate, cold):
    """Map cold enthalpy mixing or a hot counterflow flame onto the MFC grid."""
    if not np.isfinite(vort_thickness) or vort_thickness <= 0:
        raise ValueError("Vorticity thickness must be positive and finite")
    z = 0.5 * (1 - np.tanh(2 * cross_coord / vort_thickness))
    h_ox, y_ox, _, _, u_ox, _ = stream_ox
    h_fu, y_fu, _, _, u_fu, _ = stream_fu
    y = y_ox[:, None] + (y_fu - y_ox)[:, None] * z
    if not cold:
        y = _hot_mass_fractions(sol, pres, temp_ox, temp_fu, y_ox, y_fu, 10 * vort_thickness, strain_rate, z)
    h = h_ox + (h_fu - h_ox) * z
    temperature = np.empty_like(z)
    for i in range(len(z)):
        sol.HPY = h[i], pres, y[:, i]
        temperature[i] = sol.T
    if not np.all(np.isfinite(temperature)) or np.any(temperature <= 0):
        raise ValueError("Invalid temperature in mixing-layer initialization")
    return SimulationFields(z, temperature, np.full_like(z, pres), u_ox + (u_fu - u_ox) * z, y)


def streams(sol, fuel, pres, temp_ox, temp_fu, molefrac_ox, molefrac_fu, vort_thickness, mach_c, mach_fu=None):
    """Thermodynamic state and velocities for the oxidizer and fuel streams
    (temporal evolution: symmetric convective frame, u_ox = -u_fu)."""
    molefrac_di = 1 - molefrac_ox
    sol.TPX = temp_ox, pres, f"O2:{molefrac_ox}, N2:{molefrac_di}"
    y_ox = sol.Y
    h_ox = sol.enthalpy_mass
    c_ox = np.sqrt((sol.cp_mass / sol.cv_mass) * sol.P / sol.density)
    nu_ox = sol.viscosity / sol.density
    rho_ox = sol.density

    molefrac_di = 1 - molefrac_fu
    sol.TPX = temp_fu, pres, f"{fuel}:{molefrac_fu}, N2:{molefrac_di}"
    y_fu = sol.Y
    h_fu = sol.enthalpy_mass
    c_fu = np.sqrt((sol.cp_mass / sol.cv_mass) * sol.P / sol.density)
    nu_fu = sol.viscosity / sol.density
    rho_fu = sol.density

    z_st = stoichiometric_mixture_fraction(sol, y_ox, y_fu)
    print(f"Stoichiometric mixture fraction: Z_st = {z_st:.3f}")

    u_ox = 0.5 * mach_c * (c_ox + c_fu) if mach_fu is None else mach_c * c_ox
    u_fu = -u_ox if mach_fu is None else mach_fu * c_fu
    delta_u = u_ox - u_fu

    print(f"Convective Mach: Ma = {delta_u / (c_ox + c_fu)}")
    print(f"Reynolds number: Re = {0.5 * delta_u * vort_thickness / nu_ox}")
    return (
        (h_ox, y_ox, c_ox, rho_ox, u_ox, nu_ox),
        (h_fu, y_fu, c_fu, rho_fu, u_fu, nu_fu),
        z_st,
    )


def reference_fluid_properties(sol, temperature_ox, pressure, mole_fraction_ox):
    """Cheap (no equilibration/solve) reference gamma and viscosity for fluid_pp(1)."""
    sol.TPX = temperature_ox, pressure, f"O2:{mole_fraction_ox}, N2:{1 - mole_fraction_ox}"
    return {"gamma": float(sol.cp_mass / sol.cv_mass), "viscosity": float(sol.viscosity)}


def ic_cache_valid(ic_dir, file_extension, expected_lines, cache_key=None):
    """True only if IC/ has prim.1.<ext>.dat with exactly the current grid's expected
    line count AND (if `cache_key` is given) a matching .cache_key.json. A bare "IC/ is
    non-empty" check isn't enough: this same case.py is invoked with different --scale
    values by different toolchain paths (e.g. `./mfc.sh validate` during precheck uses no
    args/default scale, while a registered test passes its own --scale) that can share
    this directory, so a cache populated by one grid size must not be silently reused by a
    run expecting a different one -- reading it would desync the Fortran reader (hcid=273
    expects exactly len(cross_coord) lines) from the actual grid. The `cache_key` further
    guards against silently reusing an IC generated with a different mode (--hot vs cold)
    or different physical parameters that leave the line count unchanged."""
    path = os.path.join(ic_dir, f"prim.1.00.{file_extension}.dat")
    if not os.path.isfile(path):
        return False
    with open(path) as fh:
        if sum(1 for _ in fh) != expected_lines:
            return False
    if cache_key is not None:
        key_path = os.path.join(ic_dir, ".cache_key.json")
        if not os.path.isfile(key_path):
            return False
        try:
            with open(key_path) as fh:
                stored = json.load(fh)
        except (OSError, ValueError):
            return False
        if stored != cache_key:
            return False
    return True


def write_cache_key(ic_dir, cache_key):
    """Record the parameters an IC/ was generated with, so ic_cache_valid can detect a
    stale cache (different --hot/cold mode or physical parameters at the same grid size)."""
    with open(os.path.join(ic_dir, ".cache_key.json"), "w") as fh:
        json.dump(cache_key, fh, sort_keys=True)
