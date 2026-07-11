#!/usr/bin/env python3
"""Per-save-step snapshots for the prog_burn + jwl_afterburn products-air mixing case.

This case does NOT use jwl_reactive (JWL++) -- it uses the kinematic prog_burn
front plus the mixing-rate jwl_afterburn source. There is therefore no reaction
progress lambda field; the afterburn progress variable b (prim.8) plays the
analogous "how far has the local energy release gotten" role, but it is gated
by (1-Y) so it only activates where products have started mixing with air, not
uniformly behind the front the way lambda does in the jwl_reactive case
(see examples/2D_jwl_reactive_air_mixing for that variant).

One 3-panel PNG per saved step (pressure, products mass fraction Y, afterburn
progress b), with a title stage-labeled from the known front kinematics:
  r_front(t) = pb_D_cj * t   (exact -- prog_burn is a prescribed kinematic front)
  charge is fully swept once r_front >= CHARGE_R -> "program burn" phase ends
  and everything after is products/air mixing with the afterburn source riding
  along the interface.

Usage: ./visualize_stages.py [case_dir]
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
OUT = os.path.join(HERE, "frames")
os.makedirs(OUT, exist_ok=True)

DT = 5.0e-8
CHARGE_R = 0.05
PB_D_CJ = 6930.0


def load_merged(pre, i, t):
    files = sorted(glob.glob(f"{D}/{pre}.{i}.*.{t:06d}.dat"))
    if not files:
        return None
    rows = np.vstack([np.loadtxt(f) for f in files])
    xs, ys = np.unique(rows[:, 0]), np.unique(rows[:, 1])
    grid = np.full((len(xs), len(ys)), np.nan)
    xi = np.searchsorted(xs, rows[:, 0])
    yi = np.searchsorted(ys, rows[:, 1])
    grid[xi, yi] = rows[:, 2]
    return xs, ys, grid


def stage_label(t):
    r_front = PB_D_CJ * t
    if r_front < CHARGE_R:
        return f"program burn  (front r={r_front*100:.1f} cm < charge {CHARGE_R*100:.0f} cm)"
    return f"products/air mixing + afterburn  (front r={r_front*100:.1f} cm, past charge edge)"


GRAY, MUTED = "#444444", "#777777"
plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.edgecolor": MUTED, "axes.linewidth": 0.6, "xtick.color": GRAY, "ytick.color": GRAY, "text.color": GRAY, "axes.labelcolor": GRAY})

steps = sorted({int(f.rsplit(".", 2)[1]) for f in glob.glob(f"{D}/cons.1.*.*.dat")})
if not steps:
    sys.exit(f"no cons output found in {D} -- run the case first")

written = []
for t in steps:
    xs, ys, c1 = load_merged("cons", 1, t)
    _, _, c2 = load_merged("cons", 2, t)
    Y = c1 / (c1 + c2)
    ext = [xs[0], xs[-1], ys[0], ys[-1]]

    prim5 = load_merged("prim", 5, t)
    prim8 = load_merged("prim", 8, t)

    # Y*(1-Y): the same gate m_jwl_sources.fpp uses to activate the afterburn
    # rate (Y*(1-Y)*(1-b) > sgm_eps) -- zero in pure air/pure products, peaked
    # exactly where products and air are actually mixed. This is a sharper
    # "where is afterburn active" diagnostic than b itself, because a trace
    # background product seeding (alpha_rho ~ 1.6e-5 kg/m^3 in the nominal air
    # patch, standard practice to avoid an exact-zero field) keeps the gate
    # barely open everywhere, so b saturates to ~1 domain-wide within a few
    # jwl_ab_tau once the run has advanced -- b answers "has this cell's
    # afterburn clock run out", not "is this cell at the interface".
    mixedness = Y * (1.0 - Y)

    ncols = 1 + (prim5 is not None) + 1 + (prim8 is not None)
    fig, ax = plt.subplots(1, ncols, figsize=(4.4 * ncols, 4.2), constrained_layout=True)
    if ncols == 1:
        ax = [ax]
    fig.suptitle(f"t = {t*DT*1e6:.2f} µs  (step {t})  —  {stage_label(t*DT)}", fontsize=11, color="#222222")

    col = 0
    if prim5 is not None:
        _, _, p = prim5
        im = ax[col].imshow(p.T, origin="lower", extent=ext, cmap="magma", norm=LogNorm(vmin=1.0e5, vmax=2.0e10))
        ax[col].set_title("pressure  p [Pa]  (log)")
        fig.colorbar(im, ax=ax[col], shrink=0.85)
        col += 1

    im = ax[col].imshow(Y.T, origin="lower", extent=ext, cmap="viridis", vmin=0.0, vmax=1.0)
    ax[col].set_title("products mass fraction  Y")
    fig.colorbar(im, ax=ax[col], shrink=0.85)
    col += 1

    im = ax[col].imshow(mixedness.T, origin="lower", extent=ext, cmap="cividis", vmin=0.0, vmax=0.25)
    ax[col].set_title("Y(1-Y)  —  afterburn-active (mixed) zone")
    fig.colorbar(im, ax=ax[col], shrink=0.85)
    col += 1

    if prim8 is not None:
        _, _, b = prim8
        im = ax[col].imshow(b.T, origin="lower", extent=ext, cmap="inferno", vmin=0.0, vmax=1.0)
        ax[col].set_title("afterburn progress  b  (saturates fast, see caption)")
        fig.colorbar(im, ax=ax[col], shrink=0.85)
        col += 1

    for a in ax:
        a.add_patch(Circle((0, 0), CHARGE_R, fill=False, ec="cyan", lw=0.7, ls="--"))
        a.set_xlabel("x [m]")
        a.set_ylabel("y [m]")

    out = os.path.join(OUT, f"step_{t:06d}.png")
    fig.savefig(out, dpi=130)
    plt.close(fig)
    written.append(out)
    print(f"wrote {out}")

print(f"\n{len(written)} frames written to {OUT}/")
