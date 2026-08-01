#!/usr/bin/env python3
"""Render a mid-plane (z = 0) contour movie of the 3D JWL++ detonation burst.

Reads the per-rank ASCII output (`D/prim.<var>.<rank>.<step>.dat`, columns
x y z value), extracts the z = 0 slice of each saved frame, and animates two
panels side by side:

  * left  -- pressure [GPa], log-scaled `inferno`: the detonation front and the
             blast vented into the air;
  * right -- TNT/products mass fraction, `magma`, with the lambda = 0.5 reaction
             front drawn as a contour: the charge burning and mixing into air.

Output: detonation_burst.mp4 (falls back to .gif if ffmpeg is unavailable).

Usage: ./make_movie.py [case_dir] [--fps N] [--subframes K]
  --subframes K  temporally interpolate K-1 frames between saves for smoothness.
"""

import glob
import os
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.colors import LogNorm

# --- CLI and configuration ---
HERE = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("-") else os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "D")
FPS = int(sys.argv[sys.argv.index("--fps") + 1]) if "--fps" in sys.argv else 8
SUB = int(sys.argv[sys.argv.index("--subframes") + 1]) if "--subframes" in sys.argv else 3

DT = 3.0e-8  # timestep, matches case.py (for the frame time labels)

# prim field layout (num_fluids=2, 3D): alpha_rho1, alpha_rho2, u, v, w, p, a1, a2, lambda
P_IDX, ARHO_TNT, ARHO_AIR, LAM_IDX = 6, 1, 2, 9


# --- per-rank slice reader ---
def read_var_slice(var, step, z_target=0.0):
    """Assemble the z-nearest-0 plane of `var` at `step` from all rank files."""
    chunks = [np.loadtxt(f) for f in sorted(glob.glob(f"{D}/prim.{var}.??.{step:06d}.dat"))]

    # pick the z grid-line closest to the target plane (ranks may straddle it)
    zbest = None
    for a in chunks:
        znear = np.unique(a[:, 2])[np.argmin(np.abs(np.unique(a[:, 2]) - z_target))]
        if zbest is None or abs(znear - z_target) < abs(zbest - z_target):
            zbest = znear

    # gather the (x, y, value) triples on that plane and pivot onto the 2D grid
    on_plane = [a[np.abs(a[:, 2] - zbest) < 1e-9] for a in chunks]
    x = np.concatenate([a[:, 0] for a in on_plane])
    y = np.concatenate([a[:, 1] for a in on_plane])
    v = np.concatenate([a[:, 3] for a in on_plane])
    xu, yu = np.unique(x), np.unique(y)
    grid = np.full((len(xu), len(yu)), np.nan)
    grid[np.searchsorted(xu, x), np.searchsorted(yu, y)] = v
    return xu, yu, grid


# --- load every saved frame as (pressure[GPa], mass-fraction, lambda) ---
steps = sorted(int(f.rsplit(".", 2)[1]) for f in glob.glob(f"{D}/prim.{P_IDX}.00.*.dat"))
if not steps:
    sys.exit(f"no output in {D} -- run the case first (./mfc.sh run ... --mpi -n 8)")
print(f"{len(steps)} frames: steps {steps[0]}..{steps[-1]}")

frames = []
for s in steps:
    xu, yu, p = read_var_slice(P_IDX, s)
    _, _, art = read_var_slice(ARHO_TNT, s)
    _, _, ara = read_var_slice(ARHO_AIR, s)
    _, _, lam = read_var_slice(LAM_IDX, s)
    Y = art / np.maximum(art + ara, 1e-30)  # TNT-material (products + reactant) mass fraction
    frames.append((p / 1e9, Y, np.clip(lam, 0, 1)))
    print(f"  loaded step {s}: p_max = {p.max()/1e9:.2f} GPa, lambda_max = {lam.max():.3f}")

xu, yu = xu * 100.0, yu * 100.0  # m -> cm
extent = [yu[0], yu[-1], xu[0], xu[-1]]
pmax = max(f[0].max() for f in frames)
pfloor = 1.013e-4  # GPa, ambient ~1 atm (LogNorm floor)


# --- expand saves into a smooth playback track by temporal interpolation ---
track = []
for i in range(len(frames) - 1):
    for k in range(SUB):
        a, b, w = frames[i], frames[i + 1], k / SUB
        t = steps[i] * DT + w * (steps[i + 1] - steps[i]) * DT
        track.append(tuple((1 - w) * a[j] + w * b[j] for j in range(3)) + (t,))
track.append(frames[-1] + (steps[-1] * DT,))


# --- figure (dark theme) ---
fig, ax = plt.subplots(1, 2, figsize=(12.4, 6.0), constrained_layout=True)
fig.patch.set_facecolor("#0d0d0f")
for a in ax:
    a.set_facecolor("#0d0d0f")
    a.set_xlabel("y [cm]", color="#cfcfcf")
    a.tick_params(colors="#9a9a9a")
    for sp in a.spines.values():
        sp.set_color("#3a3a3a")
ax[0].set_ylabel("x [cm]", color="#cfcfcf")

# left: pressure (log)
im0 = ax[0].imshow(track[0][0].T, origin="lower", extent=extent, cmap="inferno", norm=LogNorm(vmin=pfloor, vmax=pmax), interpolation="bilinear", aspect="equal")
cb0 = fig.colorbar(im0, ax=ax[0], fraction=0.046, pad=0.02)
cb0.set_label("pressure [GPa]", color="#cfcfcf")
cb0.ax.tick_params(colors="#9a9a9a")

# right: products/air mass fraction, with the lambda = 0.5 burn front overlaid in update()
im1 = ax[1].imshow(track[0][1].T, origin="lower", extent=extent, cmap="magma", vmin=0, vmax=1, interpolation="bilinear", aspect="equal")
cb1 = fig.colorbar(im1, ax=ax[1], fraction=0.046, pad=0.02)
cb1.set_label("TNT / products mass fraction", color="#cfcfcf")
cb1.ax.tick_params(colors="#9a9a9a")
lam_contour = [None]

title = fig.suptitle("", color="#f0f0f0", fontsize=14)


# --- animate and save ---
def update(n):
    p, Y, lam, t = track[n]
    im0.set_data(p.T)
    im1.set_data(Y.T)
    if lam_contour[0] is not None:  # clear the previous frame's burn-front contour
        lam_contour[0].remove()
        lam_contour[0] = None
    if lam.max() > 0.5 > lam.min():  # only draw when a lambda = 0.5 crossing exists
        lam_contour[0] = ax[1].contour(yu, xu, lam, levels=[0.5], colors="#39d0ff", linewidths=1.0, alpha=0.8)
    title.set_text(f"3D TNT JWL++ detonation  •  z = 0 mid-plane  •  t = {t*1e6:5.1f} µs   (p$_{{max}}$ = {p.max():.1f} GPa)")


anim = animation.FuncAnimation(fig, update, frames=len(track), blit=False)
out = os.path.join(HERE, "detonation_burst.mp4")
try:
    anim.save(out, writer=animation.FFMpegWriter(fps=FPS, bitrate=6000), dpi=140)
except Exception as e:
    out = out.replace(".mp4", ".gif")
    print(f"ffmpeg failed ({e}); writing gif instead")
    anim.save(out, writer=animation.PillowWriter(fps=FPS), dpi=100)
print(f"\nwrote {out}  ({len(track)} played frames @ {FPS} fps)")
