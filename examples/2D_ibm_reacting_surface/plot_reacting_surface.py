#!/usr/bin/env python3
"""Plot temperature, streamwise velocity, and CO mass fraction."""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import PowerNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

case_dir = Path(__file__).resolve().parent
mfc_root = case_dir.parents[1]
sys.path.insert(0, str(mfc_root / "toolchain"))

from mfc.viz.silo_reader import assemble_silo

p0_dir = case_dir / "silo_hdf5" / "p0"
steps = sorted(int(f.stem) for f in p0_dir.glob("*.silo") if f.stem.isdigit())
if not steps:
    raise RuntimeError(f"No Silo output found in {p0_dir}")

step = steps[-1]
data = assemble_silo(str(case_dir), step)

x = np.asarray(data.x_cc).squeeze()
y = np.asarray(data.y_cc).squeeze()
V = data.variables


def orient(q):
    q = np.asarray(q).squeeze()
    return q.T if q.shape == (len(x), len(y)) else q


def first_field(*names):
    for name in names:
        if name in V:
            return name
    raise KeyError(f"Could not find any of: {', '.join(names)}")


T = orient(V[first_field("T", "temp", "temperature", "chem_T")])
u = orient(V["vel1"])
Y_CO = orient(V[first_field("CO", "Y_CO", "YCO", "chem_Y_CO", "chem_YCO", "Y1", "Y_1", "chem_Y1", "chem_Y_1")])

co_positive = Y_CO[Y_CO > 1.0e-12]
co_max = np.nanpercentile(co_positive, 99.0) if co_positive.size else np.nanmax(Y_CO)


def levels(q, vmin=None, vmax=None):
    qmin = np.nanmin(q) if vmin is None else vmin
    qmax = np.nanmax(q) if vmax is None else vmax
    return np.linspace(qmin, qmax, 31)


fields = [
    ("Temperature (K)", T, "turbo", levels(T), None),
    ("Streamwise velocity (m/s)", u, "RdBu_r", levels(u), None),
    ("CO mass fraction", Y_CO, "inferno", levels(Y_CO, 0.0, co_max), PowerNorm(gamma=0.5, vmin=0.0, vmax=co_max)),
]

fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.2))

for ax, (title, q, cmap, lev, norm) in zip(axes, fields):
    im = ax.contourf(x * 1e3, y * 1e3, q, levels=lev, cmap=cmap, norm=norm)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title(title)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="4%", pad=0.08)
    fig.colorbar(im, cax=cax)

fig.suptitle("Heterogeneous reacting flow past a carbon cylinder", fontsize=14)
plt.tight_layout(rect=[0, 0, 1, 0.95])

outfile = case_dir / f"reacting_surface_T_u_CO_{step:06d}.png"
fig.savefig(outfile, bbox_inches="tight")
print(f"Saved: {outfile}")
plt.show()
