"""Flamelet-based initial condition generator for a temporal reacting mixing layer.

Adapted from an external Cantera + Pyrometheus(flamelets) + JAX tool that used to write
MFC's internal binary restart format directly. Here the same 1-D flamelet solve instead
writes plain-text `prim.<n>.00.<ext>.dat` files read by hcid=273
(`src/common/include/2dHardcodedIC.fpp`), so `pre_process` — not this script — owns the
binary layout.

AXIS CONVENTION (do not "fix" this to look more natural -- it is required by hcid=273):
    MFC x  ->  cross-stream / flame coordinate (the flamelet profile varies here)
    MFC y  ->  streamwise (periodic, statistically homogeneous; extruded by hcid=273)
This is the reverse of the more intuitive x=streamwise labeling, because
`ExtrusionHardcodedIC.fpp` mechanically requires the profile to vary with MFC's x and be
replicated across MFC's y. See docs/documentation/case.md's hcid=273 entry.

The streamwise mean velocity (which must vary along MFC's x, i.e. the profile axis) is
written into the mom%beg (index 2) file slot -- the slot that would otherwise hold the
(always-zero, in this unperturbed base state) cross-stream velocity. `case(273)` in
2dHardcodedIC.fpp swaps it into its correct place (mom%end) at load time.
"""

import contextlib
import json
import os
import sys
import time
from dataclasses import dataclass

import jax.numpy as jnp
import numpy as np
from pyrometheus.flamelets.domain import Domain, DomainConfig
from pyrometheus.flamelets.solver import FlameletSolver
from pyrometheus.flamelets.state import FlameletState
from pyrometheus.flamelets.utils import bell_profile, stoichiometric_mixture_fraction
from scipy.interpolate import interp1d


@dataclass
class SimulationFields:
    """1-D flamelet profiles on the flame/cross-stream (MFC x) coordinate."""

    mixture_fraction: jnp.ndarray
    temperature: jnp.ndarray
    pressure: jnp.ndarray
    velocity: jnp.ndarray  # streamwise (MFC y) velocity
    mass_fractions: jnp.ndarray  # shape (Ns, nx)


def diffusivity(pyro_gas, pressure, temperature, mass_fractions):
    """Unity-Lewis diffusivity (scalar or array)."""
    k = pyro_gas.get_mixture_thermal_conductivity_mixavg(temperature, mass_fractions)
    rho = pyro_gas.get_density(pressure, temperature, mass_fractions)
    cp = pyro_gas.get_mixture_specific_heat_cp_mass(temperature, mass_fractions)
    return k / (rho * cp)


def new_dissipation_profile(z, z_st, val_st, bval):
    """Bell-shaped scalar dissipation profile with value val_st at z_st, bval at ends."""
    profile = val_st * bell_profile(z) / bell_profile(z_st)
    for ib in (0, -1):
        profile = profile.at[ib].set(bval)
    return profile


def configure_flamelet_solver():
    return {
        "max_attempts": 10,
        "bdf": {
            "maxsteps": 10,
            "time_step": 1e-5,
            "newton": {"maxiter": 20, "tol": 1e-9},
        },
        "newton": {"maxiter": 10, "tol": 1e-8},
        "eos": {
            "maxiter": 40,
            "tol": 1e-8,
            "update_size": 0.1,
            "update_method": "gauss_newton",
        },
    }


def streams(sol, fuel, pres, temp_ox, temp_fu, molefrac_ox, molefrac_fu, vort_thickness, mach_c):
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

    u_ox = 0.5 * mach_c * (c_ox + c_fu)
    u_fu = -0.5 * mach_c * (c_ox + c_fu)
    delta_u = 2 * u_ox

    print(f"Convective Mach: Ma = {delta_u / (c_ox + c_fu)}")
    print(f"Reynolds number: Re = {0.5 * delta_u * vort_thickness / nu_ox}")
    return (
        (h_ox, y_ox, c_ox, rho_ox, u_ox, nu_ox),
        (h_fu, y_fu, c_fu, rho_fu, u_fu, nu_fu),
        z_st,
    )


def equilibrium_profiles(sol, pres, h_ox, h_fu, y_ox, y_fu, mixture_fraction):
    """Equilibrium flamelet state at every mixture-fraction point."""
    num_z = len(mixture_fraction)
    linear_enthalpy = h_ox + (h_fu - h_ox) * mixture_fraction
    linear_mass_frac = (y_ox + (y_fu - y_ox) * mixture_fraction[:, None]).T

    state_equil = FlameletState(
        enthalpy=linear_enthalpy,
        mass_fractions=linear_mass_frac,
    )
    temp_equil = jnp.zeros(num_z, dtype=jnp.float64)
    for i in range(num_z):
        sol.HPY = (state_equil.enthalpy[i], pres, state_equil.mass_fractions[:, i])
        sol.equilibrate("HP")
        temp_equil = temp_equil.at[i].set(sol.T)
        state_equil.mass_fractions = state_equil.mass_fractions.at[:, i].set(sol.Y)

    print(f"Equilibrium profiles: Max temperature T_max = {temp_equil.max():.3f} [K]")
    return state_equil, temp_equil


def baseline_s_curve_state(z_st, pres, h_ox, h_fu, state_guess, temp_guess, mixture_fraction, solver_options, solver):
    """Solve the flamelet at two high dissipation rates (also JIT-warms the solver)."""

    def _diss(val_st):
        return new_dissipation_profile(mixture_fraction, z_st, val_st, 2 * solver.domain.jac[0] ** 2)

    def _visc():
        return new_dissipation_profile(mixture_fraction, z_st, 0, 0)

    t0 = time.time()
    wmp = solver.warmup(
        "flamelet_newton_step",
        state_guess,
        _diss(100),
        _visc(),
        temp_guess,
        pres,
        h_ox,
        h_fu,
    )
    wmp.enthalpy.block_until_ready()
    print(f"Warm up time: {(time.time() - t0):.4e} s")

    t0 = time.time()
    wmp, _, _ = solver.warmup(
        "flamelet_time_step",
        state_guess,
        solver_options["bdf"]["newton"]["maxiter"],
        solver_options["bdf"]["newton"]["tol"],
        state_guess,
        solver_options["bdf"]["time_step"],
        _diss(100),
        _visc(),
        temp_guess,
        pres,
        h_ox,
        h_fu,
    )
    wmp.enthalpy.block_until_ready()
    print(f"Warm up time: {(time.time() - t0):.4e} s")

    def _solve(val_st, try_newton, t_in, s_in):
        return solver.solve(
            True,
            solver_options["newton"]["maxiter"],
            solver_options["newton"]["tol"],
            solver_options["bdf"]["newton"]["maxiter"],
            solver_options["bdf"]["newton"]["tol"],
            solver_options["bdf"]["time_step"],
            solver_options["bdf"]["maxsteps"],
            try_newton,
            solver_options["max_attempts"],
            _diss(val_st),
            _visc(),
            t_in,
            pres,
            h_ox,
            h_fu,
            s_in,
        )

    t0 = time.time()
    state_base, temp_base = _solve(100, False, temp_guess, state_guess)
    state_base.enthalpy.block_until_ready()
    print(f"Run time: {(time.time() - t0):.4e} s")

    t0 = time.time()
    state_base, temp_base = _solve(1000, True, temp_base, state_base)
    state_base.enthalpy.block_until_ready()
    print(f"Run time: {(time.time() - t0):.4e} s")

    return state_base, temp_base


def find_flame_dissipation_rate(pyro_gas, z_st, mixture_fraction, cross_coord, num_iter, sim_fields_cold, state_guess, temp_guess, pres, h_ox, h_fu, solver_options, solver):
    """Iterate the scalar dissipation rate to match the cold mixture-fraction field."""
    cold_diff = diffusivity(
        pyro_gas,
        sim_fields_cold.pressure,
        sim_fields_cold.temperature,
        sim_fields_cold.mass_fractions,
    )
    dz_dy = np.gradient(sim_fields_cold.mixture_fraction, cross_coord)
    diss_rate = 2 * cold_diff * (dz_dy**2)

    sim_fields = SimulationFields(
        mixture_fraction=sim_fields_cold.mixture_fraction,
        velocity=sim_fields_cold.velocity,
        pressure=sim_fields_cold.pressure,
        temperature=sim_fields_cold.temperature,
        mass_fractions=sim_fields_cold.mass_fractions,
    )
    state_it = state_guess
    temp_it = temp_guess
    try_newton = False

    for j_chi in range(num_iter):
        _diss = interp1d(sim_fields.mixture_fraction, diss_rate, fill_value="extrapolate")(mixture_fraction)
        _diss[0] = 2 * solver.domain.jac[0] ** 2
        _diss[-1] = 2 * solver.domain.jac[0] ** 2

        state_it, temp_it = solver.solve(
            True,
            solver_options["newton"]["maxiter"],
            solver_options["newton"]["tol"],
            solver_options["bdf"]["newton"]["maxiter"],
            solver_options["bdf"]["newton"]["tol"],
            solver_options["bdf"]["time_step"],
            solver_options["bdf"]["maxsteps"],
            try_newton,
            solver_options["max_attempts"],
            _diss,
            jnp.zeros_like(_diss),
            temp_it,
            pres,
            h_ox,
            h_fu,
            state_it,
        )
        if not try_newton:
            try_newton = True

        sim_fields.temperature = interp1d(mixture_fraction, temp_it)(sim_fields.mixture_fraction)
        sim_fields.mass_fractions = interp1d(mixture_fraction, state_it.mass_fractions)(sim_fields.mixture_fraction)

        diff = diffusivity(pyro_gas, pres, sim_fields.temperature, sim_fields.mass_fractions)
        diss_rate = 2 * diff * (dz_dy**2)
        chi_st = interp1d(sim_fields.mixture_fraction, diss_rate)(z_st)
        print(f"---> j = {j_chi}, chi_max = {diss_rate.max():.4f}, " f"chi_st = {chi_st:.4f}, T_max = {temp_it.max():.4f}")

    return sim_fields


def compute_grid(vort_thickness, cross_min, cross_max, points_per_cross, stream_min, stream_max, num_y):
    """Pure grid arithmetic -- no Cantera/JAX. Always cheap to call, including on an
    IC/ cache hit, so case.py never needs to run the flamelet solve just to learn its
    own domain size.

    Returns
    -------
    cross_coord : ndarray (nx,), cross-stream cell centres (MFC x)
    grid : dict, m/x_domain/n/y_domain for case.py's case dict
    """
    cross_lo = cross_min * vort_thickness
    cross_hi = cross_max * vort_thickness
    dx = vort_thickness / points_per_cross
    # weno_order=5 needs m+1 >= num_stcls_min*weno_order (25); floor with margin so a
    # small --scale for cheap test runs can't shrink the grid below that.
    num_x = max(int(round((cross_hi - cross_lo) / dx)), 32)

    # Cell centers, matching MFC's own x_cc(i) = x_domain%beg + (i+0.5)*dx convention
    # exactly, so pre_process's grid lines up with this array cell-for-cell.
    cross_coord = cross_lo + (np.arange(num_x) + 0.5) * dx

    grid = {
        "m": num_x - 1,
        "x_domain_beg": float(cross_lo),
        "x_domain_end": float(cross_lo + num_x * dx),
        "n": num_y - 1,
        "y_domain_beg": float(stream_min * vort_thickness),
        "y_domain_end": float(stream_max * vort_thickness),
    }
    return cross_coord, grid


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


def create_simulation_fields(pyro_gas, sol, pres, temp_ox, temp_fu, cross_coord, vort_thickness, stream_ox, stream_fu, z_st, num_iter, cold):
    """1-D flamelet profiles on the cross-stream (MFC x) coordinate `cross_coord`."""
    mollifier = 0.5 * (1 - np.tanh(2 * cross_coord / vort_thickness))

    h_ox, y_ox, _, _, u_ox, nu_ox = stream_ox
    h_fu, y_fu, _, _, u_fu, _ = stream_fu

    sim_fields_cold = SimulationFields(
        mixture_fraction=mollifier,
        temperature=pyro_gas.get_temperature_from_enthalpy(
            h_ox + (h_fu - h_ox) * mollifier,
            (y_ox + (y_fu - y_ox) * mollifier[:, None]).T,
            temp_ox + (temp_fu - temp_ox) * mollifier,
        ),
        pressure=pres * jnp.ones_like(mollifier),
        velocity=(u_ox + (u_fu - u_ox) * mollifier),
        mass_fractions=(y_ox + (y_fu - y_ox) * mollifier[:, None]).T,
    )

    if cold:
        return sim_fields_cold

    domain_config = DomainConfig(num_x=101, x_l=0, x_r=1)
    domain = Domain(domain_config)
    solver = FlameletSolver(domain, pyro_gas, (y_ox, y_fu))
    mixture_fraction = jnp.array(domain.x)

    state_equil, temp_equil = equilibrium_profiles(sol, pres, h_ox, h_fu, y_ox, y_fu, mixture_fraction)
    solver_options = configure_flamelet_solver()
    # baseline_s_curve_state's return value is unused -- its purpose is JIT warmup
    # for the JAX solver calls that follow (matches the upstream tool's structure).
    baseline_s_curve_state(
        z_st,
        pres,
        h_ox,
        h_fu,
        state_equil,
        temp_equil,
        mixture_fraction,
        solver_options,
        solver,
    )

    sim_fields = find_flame_dissipation_rate(
        pyro_gas,
        z_st,
        mixture_fraction,
        cross_coord,
        num_iter,
        sim_fields_cold,
        state_equil,
        temp_equil,
        pres,
        h_ox,
        h_fu,
        solver_options,
        solver,
    )
    return sim_fields


def write_hcid_ic(output_dir, cross_coord, density, streamwise_velocity, pressure, mass_fractions, file_extension="000000"):
    """Write hcid=273 IC text files: prim.<n>.00.<ext>.dat, one `x value` pair per line.

    File order (matches eqn_idx for model_eqns=2, num_fluids=1, chemistry=T, skipping the
    mom%end slot that @:HardcodedReadValues() always zeros):
        1: density (alpha_rho(1))
        2: mom%beg slot, REPURPOSED to carry the streamwise-velocity profile
        3: pressure
        4: alpha(1) = 1
        5..4+Ns: species mass fractions, Cantera species order
    """
    os.makedirs(output_dir, exist_ok=True)
    num_species = mass_fractions.shape[0]

    columns = [density, streamwise_velocity, pressure, np.ones_like(density)]
    columns += [mass_fractions[k] for k in range(num_species)]

    for n, values in enumerate(columns, start=1):
        path = os.path.join(output_dir, f"prim.{n}.00.{file_extension}.dat")
        with open(path, "w") as fh:
            for x, v in zip(cross_coord, values):
                fh.write(f"{float(x)!r} {float(v)!r}\n")

    return len(columns)


def generate_ic_files(
    *, output_dir, sol, pyro_gas, cross_coord, pressure, temperature_ox, temperature_fu, fuel, mole_fraction_ox, mole_fraction_fu, vort_thickness, mach_c, num_iter, cold, file_extension="000000"
):
    """Run the flamelet solve (expensive when cold=False) and write hcid=273 IC
    files on `cross_coord` (from `compute_grid`, so the file spacing exactly matches
    the grid case.py declares).

    All Cantera/JAX/Pyrometheus stdout diagnostics are redirected to stderr: case.py's
    contract requires its entire stdout to be exactly one JSON line.
    """
    with contextlib.redirect_stdout(sys.stderr):
        stream_ox, stream_fu, z_st = streams(
            sol,
            fuel,
            pressure,
            temperature_ox,
            temperature_fu,
            mole_fraction_ox,
            mole_fraction_fu,
            vort_thickness,
            mach_c,
        )

        sim_fields = create_simulation_fields(
            pyro_gas,
            sol,
            pressure,
            temperature_ox,
            temperature_fu,
            cross_coord,
            vort_thickness,
            stream_ox,
            stream_fu,
            z_st,
            num_iter,
            cold,
        )

        temperature_1d = np.array(sim_fields.temperature)
        pressure_1d = np.array(sim_fields.pressure)
        velocity_1d = np.array(sim_fields.velocity)
        mass_fractions_1d = np.array(sim_fields.mass_fractions)
        density_1d = np.array(pyro_gas.get_density(pressure_1d, temperature_1d, mass_fractions_1d))

        # Fail at generation time rather than writing a non-finite IC that would only
        # surface downstream as a cryptic VCFL=Inf crash (e.g. a diverged --hot solve).
        if not all(np.all(np.isfinite(a)) for a in (temperature_1d, pressure_1d, velocity_1d, mass_fractions_1d, density_1d)):
            raise ValueError("flamelet IC solve produced non-finite values; refusing to write IC")

        write_hcid_ic(output_dir, cross_coord, density_1d, velocity_1d, pressure_1d, mass_fractions_1d, file_extension=file_extension)

        print(f"[flamelet_ic] Wrote IC to {output_dir}: " f"nx={len(cross_coord)}, T_max={temperature_1d.max():.1f} K")
