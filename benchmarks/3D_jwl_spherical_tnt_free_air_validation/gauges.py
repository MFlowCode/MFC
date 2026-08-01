#!/usr/bin/env python3
"""Reduce MFC probe histories to free-air blast metrics."""

import argparse
import csv
import math
from pathlib import Path

P_AMB = 101325.0
CHARGE_MASS_KG = 0.8534660042252272
DOMAIN_LENGTH = 0.5
TARGET_PROBE_X = (0.15, 0.25, 0.35, 0.45)
TIME_COL = 0
PRES_COL = 5
ARRIVAL_FRAC = 0.05


def read_probe(path):
    rows = []
    with path.open() as stream:
        for line in stream:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            try:
                rows.append([float(value) for value in stripped.split()])
            except ValueError:
                continue
    if not rows:
        raise ValueError(f"no numeric probe rows found in {path}")
    return rows


def positive_impulse(times, pressures):
    impulse = 0.0
    for t0, t1, p0, p1 in zip(times, times[1:], pressures, pressures[1:]):
        dp0 = max(p0 - P_AMB, 0.0)
        dp1 = max(p1 - P_AMB, 0.0)
        impulse += 0.5 * (dp0 + dp1) * (t1 - t0)
    return impulse


def nearest_cell_center(x, grid):
    dx = DOMAIN_LENGTH / (grid + 1)
    i = int(round(x / dx - 0.5))
    i = max(0, min(grid, i))
    return (i + 0.5) * dx


def probe_locations(grid):
    dx = DOMAIN_LENGTH / (grid + 1)
    yz = 0.5 * dx
    probes = []
    for probe_id, target_x in enumerate(TARGET_PROBE_X, start=1):
        probes.append({"probe_id": probe_id, "x": nearest_cell_center(target_x, grid), "y": yz, "z": yz})
    return probes


def reduce_probe(path, arrival_frac):
    rows = read_probe(path)
    times = [row[TIME_COL] for row in rows if len(row) > PRES_COL]
    pressures = [row[PRES_COL] for row in rows if len(row) > PRES_COL]
    if not pressures:
        raise ValueError(f"probe rows in {path} do not contain pressure column {PRES_COL}")

    baseline = pressures[0]
    overpressures = [pressure - baseline for pressure in pressures]
    peak_overpressure = max(overpressures)
    threshold = arrival_frac * P_AMB
    arrival = math.nan
    for time, overpressure in zip(times, overpressures):
        if overpressure > threshold:
            arrival = time
            break

    corrected_pressures = [P_AMB + overpressure for overpressure in overpressures]
    return arrival, peak_overpressure, positive_impulse(times, corrected_pressures)


def fmt(value):
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "nan"
    return f"{value:.6e}"


def main():
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dir", type=Path, default=script_dir / "D", help="directory containing probe<i>_prim.dat")
    parser.add_argument(
        "--output",
        type=Path,
        default=script_dir / "gauge_results.csv",
        help="CSV output path",
    )
    parser.add_argument(
        "--arrival-frac",
        type=float,
        default=ARRIVAL_FRAC,
        help="arrival threshold as a fraction above ambient pressure",
    )
    parser.add_argument("--grid", type=int, default=63, help="grid index used by case.py")
    args = parser.parse_args()

    rows = []
    w13 = CHARGE_MASS_KG ** (1.0 / 3.0)
    for probe in probe_locations(args.grid):
        radius = math.sqrt(probe["x"] ** 2 + probe["y"] ** 2 + probe["z"] ** 2)
        row = {
            "probe_id": probe["probe_id"],
            "x": probe["x"],
            "y": probe["y"],
            "z": probe["z"],
            "r": radius,
            "scaled_distance": radius / w13,
            "arrival_time_s": math.nan,
            "peak_overpressure_pa": math.nan,
            "impulse_pa_s": math.nan,
            "status": "missing",
        }
        probe_path = args.dir / f"probe{probe['probe_id']}_prim.dat"
        if probe_path.is_file():
            arrival, peak, impulse = reduce_probe(probe_path, args.arrival_frac)
            row.update(
                {
                    "arrival_time_s": arrival,
                    "peak_overpressure_pa": peak,
                    "impulse_pa_s": impulse,
                    "status": "ok",
                }
            )
        rows.append(row)

    fieldnames = [
        "probe_id",
        "x",
        "y",
        "z",
        "r",
        "scaled_distance",
        "arrival_time_s",
        "peak_overpressure_pa",
        "impulse_pa_s",
        "status",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print("| Gauge | r (m) | Z (m/kg^(1/3)) | MFC arrival (s) | MFC peak dp (Pa) | MFC impulse (Pa s) | Status |")
    print("|---:|---:|---:|---:|---:|---:|---|")
    for row in rows:
        print(
            f"| {row['probe_id']} | {row['r']:.6f} | {row['scaled_distance']:.6f} | "
            f"{fmt(row['arrival_time_s'])} | {fmt(row['peak_overpressure_pa'])} | "
            f"{fmt(row['impulse_pa_s'])} | {row['status']} |"
        )
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
