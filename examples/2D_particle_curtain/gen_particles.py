#!/usr/bin/env python3
# Writes <outdir>/input/particles.dat for the particle curtain (default outdir: this example's directory).
# case.py calls it when pre_process runs; the test suite calls it for the Example test.
import math
import os
import random
import sys

# Glass particles, radius 57.5 um, volume fraction 0.21 in the band x in [0, 2 mm] across the channel y in [0, 15 mm].
# In 2D each particle stands for a slab of depth charwidth, so projected particles may overlap (collisions are not modeled).
rp = 57.5e-6
x_curtain = (0.0, 2.0e-3)
y_channel = (0.0, 0.015)
vf = 0.21
charwidth = 2.5e-4  # the grid spacing
n_particles = round(vf * (x_curtain[1] - x_curtain[0]) * (y_channel[1] - y_channel[0]) * charwidth / (4.0 / 3.0 * math.pi * rp**3))


def write_particles(outdir):
    """Uniform random (seeded) positions in the curtain band, one line per particle: x, y, z, u, v, w, radius."""
    random.seed(1)
    os.makedirs(os.path.join(outdir, "input"), exist_ok=True)
    with open(os.path.join(outdir, "input", "particles.dat"), "w") as f:
        for _ in range(n_particles):
            x = random.uniform(*x_curtain)
            y = random.uniform(y_channel[0] + rp, y_channel[1] - rp)
            f.write(f"{x:.16e} {y:.16e} 0.0 0.0 0.0 0.0 {rp:.16e}\n")


if __name__ == "__main__":
    write_particles(sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__)))
