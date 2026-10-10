#!/usr/bin/env python3
# Writes <outdir>/input/particles.dat for the particle half-ring (default outdir: this example's directory).
# case.py calls it when pre_process runs; the test suite calls it for the Example test.
import math
import os
import random
import sys

# Magnesium particles, radius 50 um, in a half-ring r in [r_in, r_out] above the wall y = 0, no overlaps
rp = 50e-6
r_in, r_out = 0.01105, 0.01265
n_particles = 1612


def write_particles(outdir):
    """Random sequential (seeded) placement with a minimum spacing of one particle diameter, one line per particle:
    x, y, z, u, v, w, radius."""
    random.seed(1)
    cell = 2 * rp  # bucket size = exclusion distance
    buckets = {}
    pts = []
    while len(pts) < n_particles:
        r = math.sqrt(random.uniform(r_in**2, r_out**2))
        th = random.uniform(0.0, math.pi)
        x, y = r * math.cos(th), r * math.sin(th)
        if y < rp:
            continue
        i, j = int(x // cell), int(y // cell)
        near = (p for di in (-1, 0, 1) for dj in (-1, 0, 1) for p in buckets.get((i + di, j + dj), ()))
        if any((x - px) ** 2 + (y - py) ** 2 < cell**2 for px, py in near):
            continue
        buckets.setdefault((i, j), []).append((x, y))
        pts.append((x, y))

    os.makedirs(os.path.join(outdir, "input"), exist_ok=True)
    with open(os.path.join(outdir, "input", "particles.dat"), "w") as f:
        for x, y in pts:
            f.write(f"{x:.16e} {y:.16e} 0.0 0.0 0.0 0.0 {rp:.16e}\n")


if __name__ == "__main__":
    write_particles(sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__)))
