"""Report pressure-equilibrium errors over every saved frame (serial CPU run)."""

import argparse
from pathlib import Path

import numpy as np
from case import CELLS, LENGTH, PRESSURE, SPEED, STEPS, WIDTH, light_stream_density

root = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--tolerance", type=float, default=1e-8, help="Relative peak pressure-error acceptance criterion")
args = parser.parse_args()
rows = []
for step in range(0, STEPS + 1, STEPS // 100):
    x, p = np.loadtxt(root / "D" / f"prim.3.00.{step:06d}.dat", unpack=True)
    if x.size != CELLS or not np.isfinite(p).all() or np.any(p <= 0):
        raise ValueError(f"Invalid pressure field or grid at step {step}")
    err = (p - PRESSURE) / PRESSURE
    if step == 0:
        np.testing.assert_allclose(p, PRESSURE, rtol=1e-10, atol=0)
        velocity = np.loadtxt(root / "D/prim.2.00.000000.dat")[:, 1]
        np.testing.assert_allclose(velocity, SPEED, rtol=1e-10)
        light = light_stream_density()
        density = np.loadtxt(root / "D/prim.1.00.000000.dat")[:, 1]
        # MFC pairs each cell average with its right face in the 1-D D files.
        np.testing.assert_allclose(x, LENGTH * np.arange(1, CELLS + 1) / CELLS, atol=1e-13)
        centers = x - 0.5 * LENGTH / CELLS
        expected = light * (1.5 - 0.5 * np.tanh(np.cos(2 * np.pi * centers / LENGTH) / (2 * np.pi * WIDTH)))
        np.testing.assert_allclose(density, expected, rtol=2e-6)
    rows.append([step, step / STEPS, np.max(np.abs(err)), np.sqrt(np.mean(err**2))])
rows = np.asarray(rows)
np.savetxt(root / "pressure_history.csv", rows, delimiter=",", header="step,transits,relative_Linf,relative_RMS", comments="")
i = np.argmax(rows[:, 2])
print(f"Peak relative Linf = {rows[i,2]:.6e} at {rows[i,1]:.2f} transits; peak relative RMS = {max(rows[:,3]):.6e}")
if rows[i, 2] > args.tolerance:
    raise SystemExit(f"FAIL: pressure equilibrium exceeded tolerance {args.tolerance:g}")
print("PASS: pressure equilibrium within tolerance")
