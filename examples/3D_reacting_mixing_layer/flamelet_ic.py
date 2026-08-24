"""Flamelet-based initial condition generator for a 3D temporal reacting mixing layer.

Shares its flamelet solve (Cantera + Pyrometheus(flamelets) + JAX) with
`examples/2D_reacting_mixing_layer/flamelet_ic.py`; only the grid and IC-writing
functions differ, to extrude the 1-D profile across two directions instead of one and to
match hcid=370's file format instead of hcid=273's.

AXIS CONVENTION:
    MFC x  ->  streamwise (periodic, statistically homogeneous)
    MFC y  ->  cross-stream / flame coordinate (the flamelet profile varies here)
    MFC z  ->  spanwise (periodic, statistically homogeneous)
This is the "intuitive" x=streamwise labeling -- unlike the 2D case, hcid=370 (3D
extrusion of a 2D profile from external data, `src/common/include/3dHardcodedIC.fpp`)
reads genuine (x, y) coordinates from file and only mechanically requires the *z* axis to
be the uniform/extruded one, so x and y are free to take their natural roles here.

In-plane (x,y) velocity perturbations are baked directly into this file (`perturb_xy`,
ported from ../tmp-mixlyr/perturb.py's perturb_temporally_evolving_layer -- computes
wavenumbers relative to domain_length/vort_thickness, so it's dimensionally correct
regardless of the case's absolute length-scale choice, unlike MFC's own
`mixlayer_perturb` which hardcodes an absolute wavenumber range unsuited to this
millimeter-scale domain). The spanwise (z) component isn't -- and can't be, since hcid=370
extrudes this (x,y) file uniformly across z -- so hcid=371 (`src/common/include/
3dHardcodedIC.fpp`) additionally modulates the already-written (x,y) perturbation by a
closed-form z-dependent factor at IC-assignment time, using only the z domain's own extent
(no new case parameters needed).
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
    # Was 1000 -- for this case's stream temperatures (T_ox=500K, T_fu=300K) and
    # sandiego.yaml mechanism, chi=1000 1/s stalls the Newton solve (residual oscillates,
    # never converges): it's past extinction for this fuel/oxidizer pair, not a numerics
    # bug. This value is pure JIT-warmup scaffolding (state_base/temp_base's return here
    # is discarded -- see this function's docstring), not a physical target, so any
    # comfortably-converging value works; 150 is safely below both this case's actual
    # chi_st (~286 1/s at vort_thickness=1e-3, from find_flame_dissipation_rate) and
    # extinction.
    state_base, temp_base = _solve(150, True, temp_base, state_base)
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


def compute_grid_3d(vort_thickness, cross_min, cross_max, points_per_cross, stream_min, stream_max, num_x, span_min, span_max, num_z):
    """Pure grid arithmetic -- no Cantera/JAX. Always cheap to call, including on an
    IC/ cache hit, so case.py never needs to run the flamelet solve just to learn its
    own domain size.

    Returns
    -------
    cross_coord : ndarray (ny,), cross-stream cell centres (MFC y) -- the flamelet
        profile axis.
    x_coord : ndarray (nx,), streamwise cell centres (MFC x), uniform.
    grid : dict, m/n/p and x/y/z_domain for case.py's case dict.
    """
    cross_lo = cross_min * vort_thickness
    cross_hi = cross_max * vort_thickness
    dy = vort_thickness / points_per_cross
    # weno_order=5 needs m/n/p+1 >= num_stcls_min*weno_order (25) in every direction with
    # more than one cell; floor with margin so a small --scale can't shrink below that.
    num_y = max(int(round((cross_hi - cross_lo) / dy)), 32)
    cross_coord = cross_lo + (np.arange(num_y) + 0.5) * dy

    stream_lo = stream_min * vort_thickness
    stream_hi = stream_max * vort_thickness
    num_x = max(num_x, 32)
    dx = (stream_hi - stream_lo) / num_x
    x_coord = stream_lo + (np.arange(num_x) + 0.5) * dx

    span_lo = span_min * vort_thickness
    span_hi = span_max * vort_thickness
    num_z = max(num_z, 32)

    grid = {
        "m": num_x - 1,
        "x_domain_beg": float(stream_lo),
        "x_domain_end": float(stream_hi),
        "n": num_y - 1,
        "y_domain_beg": float(cross_lo),
        "y_domain_end": float(cross_lo + num_y * dy),
        "p": num_z - 1,
        "z_domain_beg": float(span_lo),
        "z_domain_end": float(span_hi),
    }
    return cross_coord, x_coord, grid


def reference_fluid_properties(sol, temperature_ox, pressure, mole_fraction_ox):
    """Cheap (no equilibration/solve) reference gamma and viscosity for fluid_pp(1)."""
    sol.TPX = temperature_ox, pressure, f"O2:{mole_fraction_ox}, N2:{1 - mole_fraction_ox}"
    return {"gamma": float(sol.cp_mass / sol.cv_mass), "viscosity": float(sol.viscosity)}


def ic_cache_valid(ic_dir, file_extension, expected_lines, cache_key=None):
    """True only if IC/ has prim.1.<ext>.dat with exactly the current grid's expected
    line count AND (if `cache_key` is given) a matching .cache_key.json. See
    examples/2D_reacting_mixing_layer/flamelet_ic.py's docstring for why both checks
    are needed."""
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


def perturb_xy(x_coord, cross_coord, vort_thickness, delta_u, num_modes=10, num_blocks=5, seed=None):
    """Solenoidal (x,y) velocity perturbation, ported from
    ../tmp-mixlyr/perturb.py's perturb_temporally_evolving_layer. Wavenumbers are
    computed relative to domain_length (x) and vort_thickness (y), not a hardcoded
    absolute range, so the result is dimensionally sensible regardless of the case's
    length-scale choice.

    Returns
    -------
    u_p, v_p : ndarray (nx, ny) -- streamwise/cross-stream velocity perturbations.
    """
    rng = np.random.default_rng(seed)
    domain_length = x_coord[-1] - x_coord[0]
    domain_height = cross_coord[-1] - cross_coord[0]

    stream_coord = x_coord[:, None]  # (nx, 1)
    flame_coord = cross_coord[None, :]  # (1, ny)
    dx = x_coord[1] - x_coord[0]
    dy = cross_coord[1] - cross_coord[0]

    phases = np.pi * (2 * rng.uniform(0, 1, size=num_modes) - 1)
    phase_jitter = rng.uniform(0, 1, size=num_blocks)

    x_int = np.clip(np.array(stream_coord // (domain_length / num_blocks), dtype=int), 0, num_blocks - 1)
    conditions = [phase_jitter[x_int] < 0.33, (phase_jitter[x_int] > 0.33) * (phase_jitter[x_int] < 0.66), phase_jitter[x_int] > 0.66]
    f = np.select(conditions, [1, -1, 0])

    k_x = 2 / domain_length
    k_y = 0.1 / vort_thickness

    modes = np.stack([np.cos(2 * np.pi * i * k_x * stream_coord + 2 * np.pi * k_y * flame_coord + (phases[i - 1] + np.pi * f)) for i in range(1, num_modes + 1)])
    potential = np.sum(modes, axis=0)
    dp_dx, dp_dy = np.gradient(potential, dx, dy)
    u_p = -dp_dy
    v_p = dp_dx

    # Was 0.1 (10% of delta_u, 20% of u_ox) -- negative species mass fractions and a
    # species mass fraction slightly >1 appeared by t_step=5000 at the flame sheet, likely
    # from this large a seed imposing violent initial gradients on an already-thin
    # reacting interface. 0.03 still comfortably seeds the instability (it grows
    # exponentially regardless of seed size) while reducing the initial strain.
    fac = 0.03
    midplane = np.abs(flame_coord.ravel()) <= vort_thickness
    weight = np.sqrt(np.mean(u_p[:, midplane] ** 2 + v_p[:, midplane] ** 2))
    u_pp = (fac * delta_u / weight) * u_p
    v_pp = (fac * delta_u / weight) * v_p

    alpha = 1000
    beta = 8000
    mollifier = 0.5 * (
        np.tanh(alpha * (flame_coord / domain_height) + beta * (vort_thickness / domain_height)) - np.tanh(alpha * (flame_coord / domain_height) - beta * (vort_thickness / domain_height))
    )
    return mollifier * u_pp, mollifier * v_pp


def write_hcid370_ic(output_dir, x_coord, cross_coord, density, streamwise_velocity, pressure, mass_fractions, u_perturb, v_perturb, file_extension="000000"):
    """Write hcid=370 IC text files: prim.<n>.00.<ext>.dat, one `x y value` triple per
    line, x-major/y-minor order -- the format `HardcodedReadValues()`'s num_dims==3
    branch (`src/common/include/ExtrusionHardcodedIC.fpp`) reads and extrudes uniformly
    across z. The 1-D profiles (functions of `cross_coord` only) are broadcast uniformly
    across `x_coord`, then `u_perturb`/`v_perturb` (nx, ny), from `perturb_xy`, are added
    to the velocity columns -- so unlike the base (unperturbed) hcid=370 usage, the
    written field genuinely varies with x, not just y (still uniform in z; see hcid=371
    in 3dHardcodedIC.fpp for the z-modulation that adds z-structure at IC-assignment time).

    File order (matches eqn_idx for model_eqns=2, num_fluids=1, chemistry=T, skipping the
    mom%end slot that @:HardcodedReadValues() always zeros -- here the spanwise/w
    component, since MFC x is streamwise/mom%beg and MFC y is cross-stream/mom%beg+1, so
    unlike hcid=273 no axis swap is needed):
        1: density (alpha_rho(1))
        2: mom%beg (streamwise velocity, base profile varies with y, plus u_perturb(x,y))
        3: mom%beg+1 (cross-stream velocity, base=0, plus v_perturb(x,y))
        4: pressure
        5: alpha(1) = 1
        6..5+Ns: species mass fractions, Cantera species order
    """
    os.makedirs(output_dir, exist_ok=True)
    num_species = mass_fractions.shape[0]
    nx = len(x_coord)
    ny = len(cross_coord)

    ones_2d = np.ones((nx, ny))
    u_2d = streamwise_velocity[None, :] + u_perturb
    v_2d = v_perturb
    columns_2d = [density[None, :] * ones_2d, u_2d, v_2d, pressure[None, :] * ones_2d, ones_2d]
    columns_2d += [mass_fractions[k][None, :] * ones_2d for k in range(num_species)]

    x_grid, y_grid = np.meshgrid(x_coord, cross_coord, indexing="ij")
    for n, values_2d in enumerate(columns_2d, start=1):
        path = os.path.join(output_dir, f"prim.{n}.00.{file_extension}.dat")
        with open(path, "w") as fh:
            for row_x, row_y, row_v in zip(x_grid, y_grid, values_2d):
                for x, y, v in zip(row_x, row_y, row_v):
                    fh.write(f"{float(x)!r} {float(y)!r} {float(v)!r}\n")

    return len(columns_2d), ny


def generate_ic_files(
    *,
    output_dir,
    sol,
    pyro_gas,
    cross_coord,
    x_coord,
    pressure,
    temperature_ox,
    temperature_fu,
    fuel,
    mole_fraction_ox,
    mole_fraction_fu,
    vort_thickness,
    mach_c,
    num_iter,
    cold,
    file_extension="000000",
):
    """Run the flamelet solve (expensive when cold=False) and write hcid=370 IC
    files on `cross_coord`/`x_coord` (from `compute_grid_3d`, so the file spacing exactly
    matches the grid case.py declares).

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

        delta_u = float(velocity_1d.max() - velocity_1d.min())
        u_perturb, v_perturb = perturb_xy(x_coord, np.asarray(cross_coord), vort_thickness, delta_u)
        if not (np.all(np.isfinite(u_perturb)) and np.all(np.isfinite(v_perturb))):
            raise ValueError("perturb_xy produced non-finite values; refusing to write IC")

        write_hcid370_ic(output_dir, x_coord, cross_coord, density_1d, velocity_1d, pressure_1d, mass_fractions_1d, u_perturb, v_perturb, file_extension=file_extension)

        print(f"[flamelet_ic] Wrote IC to {output_dir}: " f"nx={len(x_coord)}, ny={len(cross_coord)}, T_max={temperature_1d.max():.1f} K")
