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


def compute_grid(vort_thickness, cross_min, cross_max, points_per_cross, stream_min, stream_max, num_y):
    """Pure grid arithmetic -- no Cantera. Always cheap to call, including on an
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
    *, output_dir, sol, cross_coord, pressure, temperature_ox, temperature_fu, fuel, mole_fraction_ox, mole_fraction_fu, vort_thickness, mach_c, strain_rate, cold, file_extension="000000"
):
    """Run the flamelet solve (expensive when cold=False) and write hcid=273 IC
    files on `cross_coord` (from `compute_grid`, so the file spacing exactly matches
    the grid case.py declares).

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

        write_hcid_ic(output_dir, cross_coord, density_1d, velocity_1d, pressure_1d, mass_fractions_1d, file_extension=file_extension)

        print(f"[flamelet_ic] Wrote IC to {output_dir}: " f"nx={len(cross_coord)}, T_max={temperature_1d.max():.1f} K")
