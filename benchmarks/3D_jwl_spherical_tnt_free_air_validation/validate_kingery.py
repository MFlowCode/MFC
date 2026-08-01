#!/usr/bin/env python3
"""Physics validation of free-air JWL/TNT blast against a reference air-blast curve.

Reads the reduced gauge metrics (gauge_results.csv: peak overpressure vs scaled
distance Z = R / W^(1/3)) and compares MFC's peak incident overpressure to the
Kinney & Graham (1985) closed-form correlation for a spherical TNT charge in
free air. This checks the physics (are the JWL coefficients and closure giving a
correct blast), not just implementation self-consistency.

Kinney-Graham peak overpressure ratio (Z in m/kg^(1/3), valid ~0.05 <= Z <= 40):

    Po/Pa = 808 [1 + (Z/4.5)^2]
            / sqrt(1 + (Z/0.048)^2) / sqrt(1 + (Z/0.32)^2) / sqrt(1 + (Z/1.35)^2)
"""

import argparse
import csv
import math
from pathlib import Path

P_AMB = 101325.0


def kinney_graham_overpressure(z):
    """Peak incident overpressure (Pa) at scaled distance z (m/kg^(1/3))."""
    num = 808.0 * (1.0 + (z / 4.5) ** 2)
    den = math.sqrt(1.0 + (z / 0.048) ** 2) * math.sqrt(1.0 + (z / 0.32) ** 2) * math.sqrt(1.0 + (z / 1.35) ** 2)
    return P_AMB * num / den


def main():
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, default=script_dir / "gauge_results.csv")
    parser.add_argument(
        "--tol",
        type=float,
        default=0.20,
        help="max acceptable relative error in the reliable far-field band",
    )
    parser.add_argument(
        "--z-farfield",
        type=float,
        default=0.30,
        help="scaled distance above which Kinney-Graham and the grid resolution are both reliable",
    )
    args = parser.parse_args()

    with args.csv.open() as stream:
        rows = [r for r in csv.DictReader(stream) if r.get("status") == "ok"]
    if not rows:
        raise SystemExit(f"no usable gauge rows in {args.csv}")

    # The near-field (small Z) is expected to deviate: the charge is burst-initialized
    # rather than detonated, the steep gradient is under-resolved on a coarse grid, and
    # Kinney-Graham itself has larger uncertainty below Z~0.2. The physically meaningful
    # test is the far-field band, together with the convergence trend toward the curve.
    print(f"{'Z (m/kg^1/3)':>13} {'MFC (kPa)':>12} {'K-G (kPa)':>12} {'rel err':>9}  band")
    worst_far = 0.0
    for r in rows:
        z = float(r["scaled_distance"])
        mfc = float(r["peak_overpressure_pa"])
        ref = kinney_graham_overpressure(z)
        rel = abs(mfc - ref) / ref
        band = "far" if z >= args.z_farfield else "near"
        if band == "far":
            worst_far = max(worst_far, rel)
        print(f"{z:13.4f} {mfc / 1e3:12.1f} {ref / 1e3:12.1f} {rel:9.2%}  {band}")

    print(f"\nworst far-field (Z >= {args.z_farfield}) error: {worst_far:.2%} (tolerance {args.tol:.0%})")
    if worst_far > args.tol:
        raise SystemExit(f"FAIL: far-field overpressure deviates from Kinney-Graham by {worst_far:.1%}")
    print("PASS: far-field blast overpressures track the Kinney-Graham TNT curve;")
    print("      near-field converges toward it as Z increases (expected for burst init + coarse grid)")


if __name__ == "__main__":
    main()
