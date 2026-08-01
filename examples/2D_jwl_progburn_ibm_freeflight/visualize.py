#!/usr/bin/env python3
# Post-process figure for 2D_jwl_progburn_ibm_freeflight.
#
# Renders a four-panel pressure time series (log-scaled, perceptually uniform
# colormap) with the moving cylinder outline overlaid at its actual free-flight
# position and the charge origin marked. Run the case first, then:
#
#   python3 visualize.py
#
# Reads this case's own silo_hdf5/ output (any rank count; MPI-decomposed ranks
# are stitched by their mesh coordinates) and restart_data/ib_state_*.dat for the
# body position. Writes jwl_progburn_freeflight.png next to this script.
import glob
import os
import re

import h5py
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle

HERE = os.path.dirname(os.path.abspath(__file__))
SILO = os.path.join(HERE, "silo_hdf5")
IBD = os.path.join(HERE, "restart_data")
DT = 2.5e-8  # matches case.py dt
# Domain and grid from case.py (x,y in [-0.1, 0.1], m = n = 149 -> 150 cells).
X0, X1, NX = -0.1, 0.1, 150


def _rank_dirs():
    return sorted(glob.glob(os.path.join(SILO, "p*")))


def steps_available():
    d0 = _rank_dirs()[0]
    s = [int(re.findall(r"(\d+)\.silo", f)[0]) for f in glob.glob(os.path.join(d0, "*.silo"))]
    return sorted(s)


def _read_rank(fname):
    """Return (pressure_field[y,x], x_centers, y_centers) for one rank file."""
    with h5py.File(fname, "r") as h:
        g = h["/.silo"]
        fields, coords = [], []
        for k in sorted(g.keys(), key=lambda s: int(s.strip("#")) if s.strip("#").isdigit() else -1):
            d = g[k]
            if not isinstance(d, h5py.Dataset):
                continue
            if d.ndim == 2:
                fields.append(np.array(d).astype(float))
            elif d.ndim == 1 and 2 < d.size < 4000:
                coords.append(np.array(d).astype(float))
        # Pressure is the only 2D field whose max exceeds 5e4 Pa (density max ~5,
        # velocities ~1e3, alpha in [0,1]). Key on max, not min: rarefaction
        # behind the blast drives pressure below 1e4 Pa in some ranks, so a min
        # floor would drop those ranks and leave white gaps.
        pres = None
        for a in fields:
            if a.max() > 5e4 and a.max() < 1e10 and a.min() > -1e3:
                if pres is None or a.max() > pres.max():
                    pres = a
        # Two longest 1D arrays are the x and y cell boundaries; cell centers between them.
        cb = sorted(coords, key=len)[-2:]
        cb = sorted(cb, key=lambda a: a.size)  # smaller = y (n), larger = x (m) if non-square; else order-safe
    return pres, cb


def read_pressure(step):
    """Stitch all ranks into one global [y,x] pressure field.

    Each rank's silo carries a few ghost cells that overlap its neighbours'
    interior. Placing rank fields by index alone lets a rank's stale ghost cells
    overwrite a neighbour's real data, which corrupts the high-gradient blast
    region (visible as a broken, non-mirror-symmetric field). Resolve overlaps by
    keeping, for each global cell, the contribution that is most interior to its
    own rank (largest distance to that rank's array edge) -- ghost cells always
    lose to a neighbour's interior. Works for any rank count and decomposition.
    """
    dx = (X1 - X0) / NX
    glob_p = np.full((NX, NX), np.nan)
    glob_d = np.full((NX, NX), -1.0)  # interiorness of the current owner
    for rd in _rank_dirs():
        f = os.path.join(rd, f"{step}.silo")
        if not os.path.exists(f):
            continue
        p, cb = _read_rank(f)
        if p is None:
            continue
        by, bx = cb[0], cb[-1]
        xc = 0.5 * (bx[:-1] + bx[1:])
        yc = 0.5 * (by[:-1] + by[1:])
        p = p[: yc.size, : xc.size]  # p is [y, x]
        ny, nx = p.shape
        jj = np.round((yc[:ny] - X0) / dx - 0.5).astype(int)
        ii = np.round((xc[:nx] - X0) / dx - 0.5).astype(int)
        di = np.minimum(np.arange(nx), nx - 1 - np.arange(nx))
        dj = np.minimum(np.arange(ny), ny - 1 - np.arange(ny))
        for a in range(ny):
            j = jj[a]
            if j < 0 or j >= NX:
                continue
            for b in range(nx):
                i = ii[b]
                if i < 0 or i >= NX:
                    continue
                d = dj[a] if dj[a] < di[b] else di[b]
                if d > glob_d[j, i]:
                    glob_d[j, i] = d
                    glob_p[j, i] = p[a, b]
    return glob_p


def cylinder_x(step):
    f = os.path.join(IBD, f"ib_state_{step}.dat")
    if not os.path.exists(f):
        return 0.0
    return float(np.fromfile(f, dtype="<f8")[16])  # x-displacement


def main():
    steps = steps_available()
    # Four frames spread across the run: lit charge -> impact -> pushed -> downstream.
    pick = [steps[int(round(f * (len(steps) - 1)))] for f in (0.10, 0.30, 0.60, 1.00)]
    ext = [X0 * 1e3, X1 * 1e3, X0 * 1e3, X1 * 1e3]  # mm
    R = 15.0  # cylinder radius, mm
    vmin, vmax = 2e5, 2e8

    fig, axes = plt.subplots(1, 4, figsize=(15, 4.4), constrained_layout=True)
    im = None
    for ax, st in zip(axes, pick):
        P = np.clip(read_pressure(st), vmin, vmax)
        im = ax.imshow(P, origin="lower", extent=ext, cmap="inferno",
                       norm=LogNorm(vmin=vmin, vmax=vmax), aspect="equal", interpolation="bilinear")
        xc = cylinder_x(st) * 1e3
        ax.add_patch(Circle((xc, 0), R, fill=True, fc="0.15", alpha=0.55, zorder=2))
        ax.add_patch(Circle((xc, 0), R, fill=False, ec="white", lw=1.6, zorder=3))
        ax.plot(-55, 0, marker="*", ms=9, mfc="cyan", mec="k", mew=0.5, zorder=4)
        ax.set_title(rf"t = {st * DT * 1e6:.0f} $\mu$s", fontsize=12)
        ax.set_xlabel("x (mm)", fontsize=11)
        ax.set_xlim(-100, 100)
        ax.set_ylim(-100, 100)
        ax.tick_params(labelsize=9)
    axes[0].set_ylabel("y (mm)", fontsize=11)
    for ax in axes[1:]:
        ax.set_yticklabels([])
    cb = fig.colorbar(im, ax=axes, shrink=0.82, pad=0.01)
    cb.set_label("Pressure (Pa)", fontsize=11)
    fig.suptitle("2D JWL program-burn detonation driving a free-flight cylinder  "
                 r"(reduced-density TNT charge $\star$, mirror-symmetric about y=0)", fontsize=13)
    out = os.path.join(HERE, "jwl_progburn_freeflight.png")
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    print("wrote", out)


if __name__ == "__main__":
    main()
