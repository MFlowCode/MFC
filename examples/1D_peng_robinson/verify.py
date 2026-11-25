#!/usr/bin/env python3
"""Check the saved acoustic period and conservation after running case.py."""

from pathlib import Path

import numpy as np
from case import AMPLITUDE, CELLS, RHO, STEPS, reference_state


def main():
    directory = Path(__file__).resolve().parent / "D"
    pressure, sound = reference_state()
    scales = [RHO * AMPLITUDE, sound * AMPLITUDE, sound**2 * RHO * AMPLITUDE]
    for variable, name, scale in zip((1, 2, 3), ("density", "velocity", "pressure"), scales):
        initial = np.loadtxt(directory / f"prim.{variable}.00.000000.dat")[:, 1]
        errors = []
        for step in range(STEPS // 4, STEPS + 1, STEPS // 4):
            final = np.loadtxt(directory / f"prim.{variable}.00.{step:06d}.dat")[:, 1]
            assert np.isfinite(final).all(), f"Nonfinite {name}"
            expected = np.roll(initial, CELLS * step // STEPS)
            errors.append(np.sqrt(np.mean((final - expected) ** 2)) / scale)
        relative_error = max(errors)
        assert relative_error < 0.005, f"{name} wave error {relative_error:.3e} exceeds 0.5% of amplitude"
        print(f"{name} RMS error / wave amplitude = {relative_error:.3e}")
    for variable, name in ((1, "mass"), (2, "momentum"), (3, "energy")):
        initial = np.loadtxt(directory / f"cons.{variable}.00.000000.dat")[:, 1]
        final = np.loadtxt(directory / f"cons.{variable}.00.{STEPS:06d}.dat")[:, 1]
        relative_drift = abs(np.sum(final) - np.sum(initial)) / np.sum(np.abs(initial))
        assert relative_drift < 1e-9, f"{name} conservation error {relative_drift:.3e}"
        print(f"{name} conservation error = {relative_drift:.3e}")
    print(f"PASS at base pressure {pressure:.6g} Pa and PR sound speed {sound:.6g} m/s")


if __name__ == "__main__":
    main()
