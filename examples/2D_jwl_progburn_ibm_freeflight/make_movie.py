#!/usr/bin/env python3
# Render every saved step of the free-flight case to a frame and encode an mp4.
# Reuses the field-stitching logic from the example's visualize.py.
import glob
import os
import re
import subprocess

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
FRAMES = os.path.join(HERE, "frames")
DT = 2.5e-8
X0, X1, NX = -0.1, 0.1, 150
R = 15.0
VMIN, VMAX = 2e5, 2e8


def _rank_dirs():
    return sorted(glob.glob(os.path.join(SILO, "p*")))


def steps_available():
    d0 = _rank_dirs()[0]
    s = [int(re.findall(r"(\d+)\.silo", f)[0]) for f in glob.glob(os.path.join(d0, "*.silo"))]
    return sorted(s)


def _read_rank(fname):
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
        cb = sorted(coords, key=len)[-2:]
        cb = sorted(cb, key=lambda a: a.size)
    return pres, cb


def read_pressure(step):
    # Resolve rank overlaps by keeping the most-interior contribution per global
    # cell, so ghost cells never overwrite a neighbour's real data (see the fuller
    # note in the example's visualize.py).
    dx = (X1 - X0) / NX
    glob_p = np.full((NX, NX), np.nan)
    glob_d = np.full((NX, NX), -1.0)
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
        p = p[: yc.size, : xc.size]
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
    return float(np.fromfile(f, dtype="<f8")[16])


def main():
    os.makedirs(FRAMES, exist_ok=True)
    steps = steps_available()
    ext = [X0 * 1e3, X1 * 1e3, X0 * 1e3, X1 * 1e3]
    for i, st in enumerate(steps):
        P = np.clip(read_pressure(st), VMIN, VMAX)
        fig, ax = plt.subplots(figsize=(5.6, 5.2), constrained_layout=True)
        im = ax.imshow(P, origin="lower", extent=ext, cmap="inferno",
                       norm=LogNorm(vmin=VMIN, vmax=VMAX), aspect="equal", interpolation="bilinear")
        xc = cylinder_x(st) * 1e3
        ax.add_patch(Circle((xc, 0), R, fill=True, fc="0.12", alpha=0.6, zorder=2))
        ax.add_patch(Circle((xc, 0), R, fill=False, ec="white", lw=1.8, zorder=3))
        ax.plot(-55, 0, marker="*", ms=10, mfc="cyan", mec="k", mew=0.5, zorder=4)
        ax.set_xlim(-100, 100)
        ax.set_ylim(-100, 100)
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        ax.set_title(rf"2D JWL program-burn + free-flight cylinder    t = {st * DT * 1e6:5.1f} $\mu$s")
        cb = fig.colorbar(im, ax=ax, shrink=0.9, pad=0.02)
        cb.set_label("Pressure (Pa)")
        fig.savefig(os.path.join(FRAMES, f"frame_{i:04d}.png"), dpi=140, facecolor="white")
        plt.close(fig)
    print("rendered", len(steps), "frames")
    out = os.path.join(HERE, "jwl_progburn_freeflight.mp4")
    subprocess.run([
        "ffmpeg", "-y", "-framerate", "20", "-i", os.path.join(FRAMES, "frame_%04d.png"),
        "-c:v", "libx264", "-pix_fmt", "yuv420p", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", out,
    ], check=True)
    print("wrote", out)


if __name__ == "__main__":
    main()
