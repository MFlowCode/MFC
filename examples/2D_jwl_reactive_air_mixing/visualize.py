#!/usr/bin/env python3
"""Visualize the 2D JWL reactive-burn products-air mixing run.

Produces two figures from the D/ output:
  contours_final.png     -- pressure, products mass fraction Y, reaction
                            progress lambda, and a synthetic schlieren, at
                            the final save (the developed products-air mixing).
  contours_evolution.png -- pressure at four times, shared log scale, showing
                            initiation -> outward detonation -> products-air
                            expansion.

Field layout (num_fluids=2, model_eqns=2, jwl_reactive):
  prim.5 = pressure, prim.3/4 = velocity, prim.8 = reaction progress lambda,
  cons.1/2 = alpha_rho of TNT products / air, so Y = cons.1 / (cons.1+cons.2).

Usage: ./visualize.py [case_dir]
"""

import glob
import os
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle

HERE = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "D")
DT = 1.5e-8
T_SAVE = 250
CHARGE_R = 0.05


def load(pre, i, t):
    a = np.loadtxt(f"{D}/{pre}.{i}.00.{t:06d}.dat")
    xs, ys = np.unique(a[:, 0]), np.unique(a[:, 1])
    return xs, ys, a[:, 2].reshape(len(xs), len(ys))


def fields(t):
    xs, ys, p = load("prim", 5, t)
    _, _, u = load("prim", 3, t)
    _, _, v = load("prim", 4, t)
    _, _, lam = load("prim", 8, t)
    _, _, c1 = load("cons", 1, t)
    _, _, c2 = load("cons", 2, t)
    rho = c1 + c2
    return xs, ys, p, np.hypot(u, v), c1 / rho, lam, rho


def schlieren(rho, dx):
    gx, gy = np.gradient(rho, dx)
    g = np.hypot(gx, gy)
    return np.exp(-25.0 * g / max(g.max(), 1e-30))


GRAY, MUTED = "#444444", "#777777"
plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.edgecolor": MUTED, "axes.linewidth": 0.6, "xtick.color": GRAY, "ytick.color": GRAY, "text.color": GRAY, "axes.labelcolor": GRAY})

steps = sorted(int(f.rsplit(".", 2)[1]) for f in glob.glob(f"{D}/prim.5.00.*.dat"))
if not steps:
    sys.exit(f"no prim output found in {D} -- run the case first")

# figure 1: four fields at the final save
t = steps[-1]
xs, ys, p, umag, Y, lam, rho = fields(t)
dx = xs[1] - xs[0]
ext = [xs[0], xs[-1], ys[0], ys[-1]]

fig, ax = plt.subplots(2, 2, figsize=(11.5, 10.4), constrained_layout=True)
fig.suptitle(f"2D TNT reactive detonation into air  •  400×400  •  t = {t*DT*1e6:.0f} µs", fontsize=13, color="#222222")

im = ax[0, 0].imshow(p.T, origin="lower", extent=ext, cmap="magma", norm=LogNorm(vmin=1.0e5, vmax=max(p.max(), 2.0e6)))
ax[0, 0].set_title("pressure  p [Pa]  (log)")
fig.colorbar(im, ax=ax[0, 0], shrink=0.85)

im = ax[0, 1].imshow(Y.T, origin="lower", extent=ext, cmap="viridis", vmin=0.0, vmax=1.0)
ax[0, 1].set_title("products mass fraction  Y  (0 = air, 1 = products)")
fig.colorbar(im, ax=ax[0, 1], shrink=0.85)

im = ax[1, 0].imshow(lam.T, origin="lower", extent=ext, cmap="inferno", vmin=0.0, vmax=1.0)
ax[1, 0].set_title("reaction progress  λ  (0 = unreacted, 1 = burnt)")
fig.colorbar(im, ax=ax[1, 0], shrink=0.85)

im = ax[1, 1].imshow(schlieren(rho, dx).T, origin="lower", extent=ext, cmap="gray")
ax[1, 1].set_title("synthetic schlieren  (|∇ρ|)")
fig.colorbar(im, ax=ax[1, 1], shrink=0.85)

for a in ax.flat:
    a.add_patch(Circle((0, 0), CHARGE_R, fill=False, ec="cyan", lw=0.7, ls="--"))
    a.set_xlabel("x [m]")
    a.set_ylabel("y [m]")
out1 = os.path.join(HERE, "contours_final.png")
fig.savefig(out1, dpi=140)
print(f"wrote {out1}")

# figure 2: pressure evolution
pick = [steps[k] for k in np.linspace(0, len(steps) - 1, 4).astype(int)]
fig, ax = plt.subplots(1, 4, figsize=(15.5, 4.4), constrained_layout=True)
fig.suptitle("Pressure: initiation → outward detonation → products-air expansion", fontsize=12, color="#222222")
for a, tt in zip(ax, pick):
    xs, ys, p, *_ = fields(tt)
    im = a.imshow(p.T, origin="lower", extent=ext, cmap="magma", norm=LogNorm(vmin=1.0e5, vmax=2.0e10))
    a.add_patch(Circle((0, 0), CHARGE_R, fill=False, ec="cyan", lw=0.7, ls="--"))
    a.set_title(f"t = {tt*DT*1e6:.0f} µs")
    a.set_xlabel("x [m]")
ax[0].set_ylabel("y [m]")
fig.colorbar(im, ax=ax, shrink=0.85, label="p [Pa]")
out2 = os.path.join(HERE, "contours_evolution.png")
fig.savefig(out2, dpi=140)
print(f"wrote {out2}")
