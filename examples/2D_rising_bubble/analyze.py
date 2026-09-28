#!/usr/bin/env python3
"""
Benchmark quantities of the rising bubble against the reference data (reference/README.md): the
centroid height y_c and rise velocity v_c (bubble-volume-fraction weighted), and the circularity,
the perimeter of the circle of equal area over the length of the alpha = 0.5 contour. Reads the
restart files of a run of case.py.

usage: analyze.py CASE_DIR [--case 1|2] [--plot FILE]
"""

import argparse
import glob
import math
import os
import re

import numpy as np

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("dir", help="case directory holding simulation.inp and restart_data/")
parser.add_argument("--case", type=int, default=1, choices=[1, 2], help="benchmark test case run (default: %(default)s)")
parser.add_argument("--plot", help="write the comparison figure to this file")
args = parser.parse_args()

here = os.path.dirname(os.path.abspath(__file__))
groups = {"TP2D": "g1l7" if args.case == 1 else "g1l8", "FreeLIFE": "g2l3", "MooNMD": "g3l4"}
ref = {name: np.loadtxt(os.path.join(here, "reference", f"c{args.case}{tag}.txt")) for name, tag in groups.items()}
ref_shape = {name: np.loadtxt(os.path.join(here, "reference", f"c{args.case}{tag}s.txt")) for name, tag in groups.items()}

with open(os.path.join(args.dir, "simulation.inp")) as f:
    inp = f.read()
t_save = float(re.search(r"^\s*t_save\s*=\s*(\S+)", inp, re.M).group(1))
nx, ny = (int(re.search(rf"^\s*{v}\s*=\s*(\S+)", inp, re.M).group(1)) + 1 for v in "mn")
dx, dy = 1.0 / nx, 2.0 / ny
x, y = (np.arange(nx) + 0.5) * dx, (np.arange(ny) + 0.5) * dy


def contour_length(field, level=0.5):
    """Length of the level contour by marching squares over cell centers (linear interpolation on each cell edge)."""
    f = field - level
    length = 0.0
    for j in range(ny - 1):
        for i in range(nx - 1):
            c = [f[j, i], f[j, i + 1], f[j + 1, i + 1], f[j + 1, i]]
            p = [(x[i], y[j]), (x[i + 1], y[j]), (x[i + 1], y[j + 1]), (x[i], y[j + 1])]
            pts = []
            for e in range(4):
                a, b = c[e], c[(e + 1) % 4]
                if (a < 0) != (b < 0):
                    s = a / (a - b)
                    pts.append((p[e][0] + s * (p[(e + 1) % 4][0] - p[e][0]), p[e][1] + s * (p[(e + 1) % 4][1] - p[e][1])))
            for k in range(0, len(pts) - 1, 2):
                length += math.hypot(pts[k + 1][0] - pts[k][0], pts[k + 1][1] - pts[k][1])
    return length


steps = sorted(int(f.split("_")[-1][:-4]) for f in glob.glob(os.path.join(args.dir, "restart_data", "lustre_[0-9]*.dat")))
T, YC, VC, CIRC = [], [], [], []
for s in steps:
    # alpha_rho(1:2), mom(1:2), E, alpha(1:2), color function; x fastest
    q = np.fromfile(os.path.join(args.dir, "restart_data", f"lustre_{s}.dat"), dtype=np.float64).reshape(-1, ny, nx)
    a2 = q[6]
    area = a2.sum() * dx * dy
    T.append(s * t_save)
    YC.append((a2 * y[:, None]).sum() * dx * dy / area)
    VC.append((a2 * q[3] / (q[0] + q[1])).sum() * dx * dy / area)
    CIRC.append(2.0 * math.sqrt(math.pi * area) / contour_length(a2))
T, YC, VC, CIRC = map(np.array, (T, YC, VC, CIRC))


def at(series_t, series_v, t):
    return float(np.interp(t, series_t, series_v))


tp = ref["TP2D"]
print(f"{'quantity':<22} {'MFC':>9} {'TP2D':>9}")
for name, v, col in (("min circularity", CIRC, 2), ("max rise velocity", VC, 4)):
    k = int(np.argmin(v) if col == 2 else np.argmax(v))
    kr = int(np.argmin(tp[:, col]) if col == 2 else np.argmax(tp[:, col]))
    print(f"{name:<22} {v[k]:9.4f} {tp[kr, col]:9.4f}")
    print(f"{'  at t':<22} {T[k]:9.4f} {tp[kr, 0]:9.4f}")
print(f"{'y_c(t = 3)':<22} {at(T, YC, 3.0):9.4f} {at(tp[:, 0], tp[:, 3], 3.0):9.4f}")

if args.plot:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(2, 2, figsize=(10, 8))
    ax = axs.ravel()
    for axi, v, col, label in ((ax[0], YC, 3, "centroid height $y_c$"), (ax[1], VC, 4, "rise velocity $v_c$"), (ax[2], CIRC, 2, "circularity")):
        for (name, r), style in zip(ref.items(), ("-", "--", ":")):
            axi.plot(r[:, 0], r[:, col], "k" + style, lw=1, label=name)
        axi.plot(T, v, "C3", lw=1.5, label="MFC, all-Mach projection")
        axi.set_xlabel("t")
        axi.set_ylabel(label)
    ax[0].legend(fontsize=8)
    # Final shape, zoomed to the bubble with equal axes
    q = np.fromfile(os.path.join(args.dir, "restart_data", f"lustre_{steps[-1]}.dat"), dtype=np.float64).reshape(-1, ny, nx)
    for (name, r), mark in zip(ref_shape.items(), ("k.", "C0.", "C2.")):
        ax[3].plot(r[:, 0], r[:, 1], mark, ms=1, label=name)
    cs = ax[3].contour(x, y, q[6], levels=[0.5], colors="C3", linewidths=1.5)
    pts = np.vstack([r for r in ref_shape.values()] + [seg for seg in cs.allsegs[0] if len(seg)])
    lo, hi = pts.min(axis=0), pts.max(axis=0)
    pad = 0.1 * (hi - lo).max()
    ax[3].set_xlim(lo[0] - pad, hi[0] + pad)
    ax[3].set_ylim(lo[1] - pad, hi[1] + pad)
    ax[3].set_aspect("equal")
    ax[3].set_xlabel("x")
    ax[3].set_ylabel("y")
    ax[3].set_title(f"bubble shape, t = {T[-1]:.2f}")
    ax[3].plot([], [], "C3", lw=1.5, label="MFC")
    ax[3].legend(fontsize=8, markerscale=6)
    fig.tight_layout()
    fig.savefig(args.plot, dpi=150)
