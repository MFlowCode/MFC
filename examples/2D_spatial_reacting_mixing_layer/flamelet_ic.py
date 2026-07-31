"""Flamelet-based initial/inflow condition generator for a spatially-evolving reacting
mixing layer.

The 1-D flamelet solver core here is identical to (and kept in sync by hand with)
examples/2D_reacting_mixing_layer/flamelet_ic.py -- duplicated rather than imported
across example directories, since each MFC example is a self-contained, copyable unit
(no example currently imports from a sibling directory). Only the writer and grid
functions differ, because this case needs a genuinely full 2D field, not an extrusion.

AXIS CONVENTION (the reverse of the temporal example's x/y roles, and deliberately so):
    MFC x  ->  streamwise (matches spatial_bf's hardcoded advecting-direction convention
               in src/simulation/m_body_forces.fpp: theta_x includes a conv_vel*t term,
               theta_y does not)
    MFC y  ->  cross-stream (the flamelet profile varies here)
Uses hcid=274 (src/common/include/2dHardcodedIC.fpp): a full (x,y) field read with no
extrusion assumption, unlike hcid=273 (which mechanically requires cross-stream=x and
would conflict with spatial_bf's fixed axis convention -- see case.py's module
docstring).
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
    """1-D flamelet profiles on the flame/cross-stream (MFC y) coordinate."""

    mixture_fraction: jnp.ndarray
    temperature: jnp.ndarray
    pressure: jnp.ndarray
    velocity: jnp.ndarray  # streamwise (MFC x) velocity
    mass_fractions: jnp.ndarray  # shape (Ns, ny)


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


def streams(sol, fuel, pres, temp_ox, temp_fu, molefrac_ox, molefrac_fu, vort_thickness, mach_ox, mach_fu):
    """Thermodynamic state and velocities for the oxidizer and fuel streams, spatial
    (lab-frame) evolution: u_ox = mach_ox*c_ox, u_fu = mach_fu*c_fu -- both streams
    co-flow at their own absolute speed, giving a nonzero mean/convective velocity
    (unlike the temporal case's symmetric convective frame, where u_ox = -u_fu by
    construction and the mean is always zero -- unusable for spatial_bf%conv_vel,
    which needs to be nonzero to drive the forcing)."""
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

    u_ox = mach_ox * c_ox
    u_fu = mach_fu * c_fu
    delta_u = u_ox - u_fu

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
    run expecting a different one -- reading it would desync the Fortran reader (hcid=274
    expects exactly (m_glb+1)*(n_glb+1) lines) from the actual grid. The `cache_key`
    further guards against silently reusing an IC generated with a different mode (--hot vs
    cold) or different physical parameters that leave the line count unchanged."""
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
    """1-D flamelet profiles on the cross-stream (MFC y) coordinate `cross_coord`."""
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


def compute_grid_spatial(vort_thickness, cross_min, cross_max, points_per_cross, stream_min, stream_max, points_per_stream):
    """Pure grid arithmetic -- no Cantera/JAX. Always cheap to call, including on an
    IC/ cache hit, so case.py never needs to run the flamelet solve just to learn its
    own domain size.

    Returns
    -------
    stream_coord : ndarray (nx,), streamwise cell centres (MFC x)
    cross_coord  : ndarray (ny,), cross-stream cell centres (MFC y)
    grid : dict, m/x_domain/n/y_domain for case.py's case dict
    """

    def _axis(lo_mult, hi_mult, points_per):
        lo = lo_mult * vort_thickness
        hi = hi_mult * vort_thickness
        dx = vort_thickness / points_per
        # weno_order=5 needs m+1 (or n+1) >= num_stcls_min*weno_order (25).
        num = max(int(round((hi - lo) / dx)), 32)
        coord = lo + (np.arange(num) + 0.5) * dx
        return coord, lo, dx, num

    stream_coord, stream_lo, stream_dx, num_x = _axis(stream_min, stream_max, points_per_stream)
    cross_coord, cross_lo, cross_dx, num_y = _axis(cross_min, cross_max, points_per_cross)

    grid = {
        "m": num_x - 1,
        "x_domain_beg": float(stream_lo),
        "x_domain_end": float(stream_lo + num_x * stream_dx),
        "n": num_y - 1,
        "y_domain_beg": float(cross_lo),
        "y_domain_end": float(cross_lo + num_y * cross_dx),
    }
    return stream_coord, cross_coord, grid


def write_hcid274_ic(output_dir, stream_coord, cross_coord, density, streamwise_velocity, pressure, mass_fractions, file_extension="000000"):
    """Write hcid=274 IC text files: prim.<n>.00.<ext>.dat, one `x y value` triple per
    line, x-major order (outer loop streamwise/x, inner loop cross-stream/y). This is a
    genuinely full 2D field -- no extrusion, no zeroed/repurposed component, all
    sys_size variables written directly in eqn_idx order:
        1: density (alpha_rho(1))
        2: streamwise velocity (mom%beg, MFC x-velocity)
        3: cross-stream velocity (mom%end, MFC y-velocity) -- 0 everywhere at t=0
        4: pressure
        5: alpha(1) = 1
        6..5+Ns: species mass fractions, Cantera species order

    The cross-stream profile is uniform along the streamwise axis at t=0 (the flow only
    develops streamwise variation once the simulation -- inflow BC + spatial_bf forcing --
    starts evolving it).
    """
    os.makedirs(output_dir, exist_ok=True)
    num_species = mass_fractions.shape[0]

    columns = [density, streamwise_velocity, np.zeros_like(density), pressure, np.ones_like(density)]
    columns += [mass_fractions[k] for k in range(num_species)]

    for n, values in enumerate(columns, start=1):
        path = os.path.join(output_dir, f"prim.{n}.00.{file_extension}.dat")
        with open(path, "w") as fh:
            for x in stream_coord:
                for y, v in zip(cross_coord, values):
                    fh.write(f"{float(x)!r} {float(y)!r} {float(v)!r}\n")

    return len(columns)


def generate_ic_files_spatial(
    *,
    output_dir,
    sol,
    pyro_gas,
    stream_coord,
    cross_coord,
    pressure,
    temperature_ox,
    temperature_fu,
    fuel,
    mole_fraction_ox,
    mole_fraction_fu,
    vort_thickness,
    mach_ox,
    mach_fu,
    num_iter,
    cold,
    file_extension="000000",
):
    """Run the flamelet solve (expensive when cold=False) and write hcid=274 IC files.
    The t=0 field is still purely a cross-stream profile, uniform along the streamwise
    axis -- it only develops streamwise variation once simulation starts evolving it
    via the inflow BC (bc_x%beg=-17) and spatial_bf forcing.

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
            mach_ox,
            mach_fu,
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

        write_hcid274_ic(output_dir, stream_coord, cross_coord, density_1d, velocity_1d, pressure_1d, mass_fractions_1d, file_extension=file_extension)

        print(f"[flamelet_ic] Wrote spatial IC to {output_dir}: " f"nx={len(stream_coord)}, ny={len(cross_coord)}, T_max={temperature_1d.max():.1f} K")
