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


def compute_grid_spatial(vort_thickness, cross_min, cross_max, points_per_cross, stream_min, stream_max, points_per_stream):
    """Pure grid arithmetic -- no Cantera. Always cheap to call, including on an
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
    strain_rate,
    cold,
    file_extension="000000",
):
    """Run the flamelet solve (expensive when cold=False) and write hcid=274 IC files.
    The t=0 field is still purely a cross-stream profile, uniform along the streamwise
    axis -- it only develops streamwise variation once simulation starts evolving it
    via the inflow BC (bc_x%beg=-17) and spatial_bf forcing.

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
            mach_ox,
            mach_fu,
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

        write_hcid274_ic(output_dir, stream_coord, cross_coord, density_1d, velocity_1d, pressure_1d, mass_fractions_1d, file_extension=file_extension)

        print(f"[flamelet_ic] Wrote spatial IC to {output_dir}: " f"nx={len(stream_coord)}, ny={len(cross_coord)}, T_max={temperature_1d.max():.1f} K")
