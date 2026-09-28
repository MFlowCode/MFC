"""Grid and MFC file layout for Cantera mixing-layer initial profiles.

The common thermodynamic and flame solve lives in toolchain/mfc/flamelet.py.
"""

import contextlib
import os
import sys
from pathlib import Path

import numpy as np

# Also support direct `python case.py` from a source checkout.
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "toolchain"))
from mfc.flamelet import (
    INITIALIZER_VERSION,
    create_simulation_fields,
    density,
    ic_cache_valid,
    mechanism_fingerprint,
    reference_fluid_properties,
    streams,
    write_cache_key,
)


def compute_grid_3d(vort_thickness, cross_min, cross_max, points_per_cross, stream_min, stream_max, num_x, span_min, span_max, num_z):
    """Pure grid arithmetic -- no Cantera. Always cheap to call, including on an
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


def perturb_xy(x_coord, cross_coord, vort_thickness, delta_u, seed, num_modes=10, num_blocks=5):
    """Solenoidal (x,y) velocity perturbation. Wavenumbers are computed relative to
    domain_length (x) and vort_thickness (y), not a hardcoded absolute range, so the
    result is dimensionally sensible regardless of the case's length-scale choice.

    `seed` is mandatory: IC/ is gitignored and regenerated on every checkout, so a
    default would reseed the mode phases from OS entropy and make the example
    irreproducible.

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

    # Bin relative to x_coord[0]: x starts negative here, and binning the raw coordinate
    # clips every x < 0 into block 0, leaving two of five blocks unused.
    x_int = np.clip(np.array((stream_coord - x_coord[0]) // (domain_length / num_blocks), dtype=int), 0, num_blocks - 1)
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

    # savetxt, not a per-element write: the 210x560 grid is ~1.6M rows across 14 files.
    # "%.17g" round-trips float64, so the recovered values match repr()'s exactly.
    x_grid, y_grid = np.meshgrid(x_coord, cross_coord, indexing="ij")
    xy_flat = np.column_stack((x_grid.ravel(), y_grid.ravel()))
    for n, values_2d in enumerate(columns_2d, start=1):
        path = os.path.join(output_dir, f"prim.{n}.00.{file_extension}.dat")
        rows = np.column_stack((xy_flat, np.asarray(values_2d, dtype=float).ravel()))
        np.savetxt(path, rows, fmt="%.17g")

    return len(columns_2d), ny


def generate_ic_files(
    *,
    output_dir,
    sol,
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
    strain_rate,
    cold,
    perturb_seed,
    file_extension="000000",
):
    """Run the flamelet solve (expensive when cold=False) and write hcid=370 IC
    files on `cross_coord`/`x_coord` (from `compute_grid_3d`, so the file spacing exactly
    matches the grid case.py declares).

    All Cantera stdout diagnostics are redirected to stderr: case.py's
    contract requires its entire stdout to be exactly one JSON line.
    """
    with contextlib.redirect_stdout(sys.stderr):
        stream_ox, stream_fu, _ = streams(
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
            sol,
            pressure,
            temperature_ox,
            temperature_fu,
            cross_coord,
            vort_thickness,
            stream_ox,
            stream_fu,
            strain_rate,
            cold,
        )

        temperature_1d = np.array(sim_fields.temperature)
        pressure_1d = np.array(sim_fields.pressure)
        velocity_1d = np.array(sim_fields.velocity)
        mass_fractions_1d = np.array(sim_fields.mass_fractions)
        density_1d = np.array(density(sol, pressure_1d, temperature_1d, mass_fractions_1d))

        # Fail at generation time rather than writing a non-finite IC that would only
        # surface downstream as a cryptic VCFL=Inf crash (e.g. a diverged --hot solve).
        if not all(np.all(np.isfinite(a)) for a in (temperature_1d, pressure_1d, velocity_1d, mass_fractions_1d, density_1d)):
            raise ValueError("flamelet IC solve produced non-finite values; refusing to write IC")

        delta_u = float(velocity_1d.max() - velocity_1d.min())
        u_perturb, v_perturb = perturb_xy(x_coord, np.asarray(cross_coord), vort_thickness, delta_u, perturb_seed)
        if not (np.all(np.isfinite(u_perturb)) and np.all(np.isfinite(v_perturb))):
            raise ValueError("perturb_xy produced non-finite values; refusing to write IC")

        write_hcid370_ic(output_dir, x_coord, cross_coord, density_1d, velocity_1d, pressure_1d, mass_fractions_1d, u_perturb, v_perturb, file_extension=file_extension)

        print(f"[flamelet_ic] Wrote IC to {output_dir}: " f"nx={len(x_coord)}, ny={len(cross_coord)}, T_max={temperature_1d.max():.1f} K")
