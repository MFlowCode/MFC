#!/usr/bin/env python3
"""2D TNT detonation blast: pressure field and azimuthally-averaged radial front.

Reads the MFC pressure output (prim.5 for a 2-fluid 2D case) and shows (left) the
circular blast pressure field at a late frame and (right) the azimuthally-averaged
radial pressure profile at several times -- the outward-propagating detonation
front, with the thin azimuthal spread confirming radial symmetry.

Run from this directory:  python3 plot_blast.py  ->  figures/blast.png
"""

import glob
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "D")
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)
P_AMB = 101325.0


def load_field(step):
    d = np.loadtxt(os.path.join(D, f"prim.5.00.{step:06d}.dat"))
    xs = np.unique(d[:, 0])
    ys = np.unique(d[:, 1])
    P = d[:, 2].reshape(len(xs), len(ys))
    return xs, ys, P


steps = sorted(int(os.path.basename(f).split(".")[3]) for f in glob.glob(os.path.join(D, "prim.5.00.*.dat")))
dt, t_save = 2.5e-8, 1000  # from case.py

fig, ax = plt.subplots(1, 2, figsize=(13, 5.5))

# (left) pressure field at a late frame (overpressure in MPa, capped for contrast)
sL = steps[len(steps) * 3 // 4]
xs, ys, P = load_field(sL)
op = (P - P_AMB) / 1e6
im = ax[0].pcolormesh(xs, ys, op.T, cmap="inferno", shading="auto", vmin=0, vmax=np.percentile(op, 99.5))
ax[0].set_aspect("equal")
ax[0].set_title(f"(A) Blast overpressure field  t = {sL * dt * 1e6:.1f} $\\mu$s")
ax[0].set_xlabel("x [m]")
ax[0].set_ylabel("y [m]")
fig.colorbar(im, ax=ax[0], label="overpressure [MPa]")

# (right) azimuthally-averaged radial pressure at several times
X, Y = np.meshgrid(xs, ys, indexing="ij")
R = np.sqrt(X**2 + Y**2)
bins = np.linspace(0, 1.0, 120)
rc = 0.5 * (bins[:-1] + bins[1:])
for s in steps[1::2]:
    _, _, P = load_field(s)
    idx = np.digitize(R.ravel(), bins) - 1
    prof = np.array([P.ravel()[idx == i].mean() if np.any(idx == i) else np.nan for i in range(len(rc))])
    ax[1].plot(rc, (prof - P_AMB) / 1e6, lw=1.6, label=f"{s * dt * 1e6:.0f} $\\mu$s")
ax[1].set_title("(B) Azimuthally-averaged radial overpressure\n(outward-propagating detonation front)")
ax[1].set_xlabel("radius [m]")
ax[1].set_ylabel("overpressure [MPa]")
ax[1].legend(title="time", fontsize=9)
ax[1].grid(alpha=0.3)

fig.suptitle("2D TNT detonation (program burn) in air -- JWL products expansion", fontsize=13)
fig.tight_layout(rect=(0, 0, 1, 0.96))
out = os.path.join(FIG, "blast.png")
fig.savefig(out, dpi=130)
print(f"wrote {out}")

# front position vs time (peak overpressure radius) -> front speed sanity check
print("  front (peak-overpressure radius) vs time:")
for s in steps[1:]:
    _, _, P = load_field(s)
    idx = np.digitize(R.ravel(), bins) - 1
    prof = np.array([P.ravel()[idx == i].mean() if np.any(idx == i) else 0.0 for i in range(len(rc))])
    rf = rc[np.nanargmax(prof)]
    print(f"    t = {s * dt * 1e6:6.1f} us   front R = {rf:.3f} m   peak dP = {(np.nanmax(prof) - P_AMB) / 1e6:6.2f} MPa")
