#!/usr/bin/env python3
"""Measure a thermocapillary drop's rise from MFC restart data, for both validation cases.

One tool, two reporting modes (auto-detected from the domain, or forced with a 2nd argument):

  fig5  (case_Ma_0p001.py, Samareh Sec 4.1.1 / Fig 5) -- v/v_YGB vs t/t_r. The drop rises in +y and is
        tracked by its color function. In the slip-wall box (Samareh's geometry, default) the box is
        the rest frame, so the lab-frame rise velocity is the measure; in an open box we report the
        drift-corrected U = u_drop - u_far. Reports the terminal velocity (post-overshoot minimum
        of the smoothed curve), the overshoot peak, the endpoint, and the late-time drift slope.
        Compare against Samareh 2D ~ 0.80.
  fig7  (case_Ma_20.py, Sec 4.1.2 / Fig 7, Nas & Tryggvason) -- U*=U/U_r vs t*=t/t_r, with
        U_r=|sigma_T gradT| r0/mu_b and t_r=mu_b/|sigma_T gradT|. Reports the overshoot peak and the
        terminal (final-t_r-window) value. Compare the peak against ~0.13.

Mode auto-detection uses the rise-axis extent Ly: fig5 is Ly=7.5D, fig7 is Ly=4D. Pass
fig5|fig7 as the 2nd argument to override.

All run-dependent constants come from simulation.inp / pre_process.inp so this can't silently
disagree with the data. The conserved-variable layout is (model_eqns=3, num_fluids=2,
surface_tension): index 0,1 = partial densities, 3 = y-momentum, and the color function is the last
conserved variable.

Usage:  python3 measure.py [case_dir] [fig5|fig7]
Writes: <case_dir>/viz/<mode>*.png  and prints a JSON summary line (tag: RESULT_JSON).
"""

import glob
import json
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

case_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
mode_arg = sys.argv[2].lower() if len(sys.argv) > 2 else None

R = 0.5  # droplet radius (D = 1)


def read_namelist(path):
    """Parse a Fortran namelist file's plain "name = value" lines into a dict (lowercase keys)."""
    out = {}
    with open(path) as fh:
        for line in fh:
            if "=" in line:
                name, value = line.split("=", 1)
                out[name.strip().lower()] = value.strip().rstrip(",")
    return out


params = read_namelist(os.path.join(case_dir, "simulation.inp"))
patches = read_namelist(os.path.join(case_dir, "pre_process.inp"))


def param(name, src=params):
    """Look up a parameter by name and return it as a float."""
    return float(src[name.lower()])


nx = int(param("m")) + 1  # box-width cells (x)
ny = int(param("n")) + 1  # rise-axis cells (y)
nz = int(param("p")) + 1  # short-axis cells (z); 1 in 2D
dt = param("dt")
Wx = param("x_domain%end", patches) - param("x_domain%beg", patches)  # box width
Ly = param("y_domain%end", patches) - param("y_domain%beg", patches)  # rise-axis extent
dim = 3 if int(param("p")) > 0 else 2

# Mode: explicit 2nd arg, else auto-detect from the rise-axis extent.
if mode_arg in ("fig5", "fig7"):
    mode = mode_arg
elif Ly > 6.0:
    mode = "fig5"  # 7.5D box
else:
    mode = "fig7"  # 4D box

restart_dir = os.path.join(case_dir, "restart_data")
# Cell-center y positions from the boundary file (last ny+1 boundaries are the interior).
yb = np.fromfile(os.path.join(restart_dir, "lustre_y_cb.dat"), np.float64)[-(ny + 1) :]
y = 0.5 * (yb[:-1] + yb[1:])  # shape (ny,)

steps = sorted(int(m.group(1)) for f in glob.glob(os.path.join(restart_dir, "lustre_*.dat")) if (m := re.search(r"lustre_(\d+)\.dat$", f)))
if not steps:
    sys.exit(f"no lustre_<step>.dat snapshots in {restart_dir!r} -- run the case first")

cells = nx * ny * nz
nvars = np.fromfile(os.path.join(restart_dir, f"lustre_{steps[0]}.dat"), np.float64).size // cells
c_idx = nvars - 1  # color function is the last conserved variable


def field(snapshot, i):
    """Variable i of a snapshot, reshaped to (nz, ny, nx) -- x is MFC's fast (innermost) index."""
    return snapshot[i * cells : (i + 1) * cells].reshape(nz, ny, nx)


# Shared core: walk the snapshots, building color-weighted time series. Identical math in 2D/3D.
yb3 = y[None, :, None]  # broadcast y over (nz, ny, nx)
times, y_centroid, u_lab, u_far, rho_drop_t = [], [], [], [], []
for s in steps:
    snap = np.fromfile(os.path.join(restart_dir, f"lustre_{s}.dat"), np.float64)
    rho = field(snap, 0) + field(snap, 1)  # total density (sum of partial densities)
    vy = field(snap, 3) / rho  # rise (y) velocity
    c = np.clip(field(snap, c_idx), 0.0, 1.0)  # color function (1 inside the drop)
    csum = c.sum()
    times.append(s * dt)
    y_centroid.append((c * yb3).sum() / csum)  # color-weighted y-centroid
    u_lab.append((c * vy).sum() / csum)  # color-weighted drop rise velocity
    rho_drop_t.append((c * rho).sum() / csum)  # color-weighted drop density
times, y_centroid, u_lab, rho_drop_t = map(np.array, (times, y_centroid, u_lab, rho_drop_t))

viz_dir = os.path.join(case_dir, "viz")
os.makedirs(viz_dir, exist_ok=True)


def reference_scales():
    """Marangoni velocity/time scales U_r, t_r and the YGB speed, from the run's constants."""
    mu = 1.0 / param("fluid_pp(1)%re(1)")  # MFC takes Re = 1/mu
    dsigma_dT = param("sigma_dtdt")  # sigma(T) slope
    # gradT: prefer the imposed wall temperatures (conduction cases); else the IC construction 1/Ly.
    if "bc_y%twall_out" in params and "bc_y%twall_in" in params:
        gradT = abs(param("bc_y%twall_out") - param("bc_y%twall_in")) / Ly
    else:
        gradT = 2.0 / 15.0 if mode == "fig5" else 1.0 / Ly  # fig5 uses 1/(7.5D)=2/15 exactly
    G = abs(dsigma_dT * gradT)
    return mu, dsigma_dT, gradT, G, G * R / mu, mu / G  # mu, sigma_T, gradT, G, U_r, t_r


if mode == "fig5":
    mu, dsigma_dT, gradT, _, _, t_r = reference_scales()
    wall = int(param("bc_y%beg")) in (-2, -15, -16)  # walled box (Samareh) vs open box
    far_y = 0.75 * (Ly / 2.0)
    is_far = np.abs(y) > far_y  # a cell is "far field" once |y| exceeds this (open-box drift)
    # far-field y-velocity per snapshot (open-box return-flow drift correction)
    u_far = []
    for s in steps:
        snap = np.fromfile(os.path.join(restart_dir, f"lustre_{s}.dat"), np.float64)
        u_far.append((field(snap, 3) / (field(snap, 0) + field(snap, 1)))[:, is_far, :].mean())
    u_far = np.array(u_far)
    U = u_lab if wall else u_lab - u_far
    v_YGB = (2.0 / 15.0) * (-dsigma_dT) * gradT * R / mu  # YGB terminal speed (mu* = k* = 1)
    rho_drop = float(rho_drop_t[0])
    tau = rho_drop * R**2 / mu  # viscous time
    samareh_ratio = 0.80

    # Terminal velocity = post-overshoot minimum of the curve smoothed over ~one viscous time
    # (the comparison to Samareh's invariant-T terminal velocity; the endpoint drifts as T distorts).
    ratio_t = U / v_YGB
    dt_snap = times[1] - times[0] if len(times) > 1 else tau
    smooth_w = max(1, int(round(tau / dt_snap)))
    if smooth_w > 1:  # edge-corrected boxcar so the ends aren't depressed by zero-padding
        k = np.ones(smooth_w)
        sm = np.convolve(ratio_t, k, mode="same") / np.convolve(np.ones_like(ratio_t), k, mode="same")
    else:
        sm = ratio_t
    peak_i = int(np.argmax(np.where(times < 0.6 * t_r, sm, -np.inf)))
    plateau_i = peak_i + int(np.argmin(sm[peak_i:]))
    ratio_plateau, t_plateau_tr, overshoot, ratio_final = float(sm[plateau_i]), float(times[plateau_i] / t_r), float(sm[peak_i]), float(ratio_t[-1])
    in_tail = times >= times[-1] - tau
    slope_per_tr = float(np.polyfit(times[in_tail] / t_r, ratio_t[in_tail], 1)[0]) if in_tail.sum() > 1 else 0.0

    mode_label = "slip-wall box (lab-frame U)" if wall else "open box (drift-corrected U)"
    print(f"[fig5] dim={dim}D  nx={nx} ny={ny} nz={nz}  cells={cells}  nvars={nvars}  {mode_label}")
    print(f"rho_drop={rho_drop:.4f}  tau={tau:.3f}  t_r={t_r:.3f}  v_YGB={v_YGB:.6f}  run length={times[-1] / t_r:.2f} t_r\n")
    print(f"{'step':>7} {'t':>8} {'t/t_r':>6} {'y_drop':>9} {'u_lab':>9} {'u_far':>9} {'U':>9} {'U/vYGB':>8}")
    for i, t in enumerate(times):
        print(f"{steps[i]:>7} {t:>8.4f} {t / t_r:>6.2f} {y_centroid[i]:>9.5f} {u_lab[i]:>9.5f} {u_far[i]:>9.5f} {U[i]:>9.5f} {U[i] / v_YGB:>+8.2f}")
    print(f"\nrising (+y, toward hot top): {y_centroid[-1] > y_centroid[0]}")
    print(f"terminal velocity  v_t/v_YGB = {ratio_plateau:+.3f}  (at t/t_r = {t_plateau_tr:.2f})   [Samareh {dim}D ~ {samareh_ratio:.2f}]")
    print(f"overshoot peak = {overshoot:+.3f}   endpoint = {ratio_final:+.3f}   late-time drift = {slope_per_tr:+.3f}/t_r")
    if abs(slope_per_tr) > 0.05:
        print(f"  NOTE: endpoint is still drifting ({slope_per_tr:+.3f}/t_r) -- the plateau, not the endpoint, is the Samareh comparison")

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.plot(times / t_r, ratio_t, "o-", color="C0", ms=3.0, lw=1.0, label="MFC " + ("(lab frame)" if wall else "(drift-corrected)"))
    ax.axhline(1.0, ls=":", color="0.45", lw=1.3)
    ax.text(0.02, 1.0, r" $v_{\mathrm{YGB}}$ (zero-Ma Stokes, sphere)", va="bottom", ha="left", color="0.4", fontsize=9)
    ax.axhline(samareh_ratio, ls="--", color="C3", lw=1.3, label=rf"Samareh {dim}D $\approx$ {samareh_ratio:.2f}")
    ax.plot(t_plateau_tr, ratio_plateau, "*", color="k", ms=14, zorder=5, label=rf"terminal velocity = {ratio_plateau:.2f}")
    ax.set_xlabel(r"$t / t_r$  ($t_r = \mu / |\sigma_T \nabla T|$, Samareh time scale)")
    ax.set_ylabel(r"normalized rise velocity  $v / v_{\mathrm{YGB}}$")
    ax.set_xlim(left=0.0)
    ax.set_ylim(0.0, max(1.1, 1.05 * ratio_t.max()))
    ax.set_title(f"{dim}D thermocapillary rise ({'Samareh slip-wall box' if wall else 'open box'})")
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right", fontsize=9)
    fig.tight_layout()
    out_png = os.path.join(viz_dir, "rise_velocity.png")
    fig.savefig(out_png, dpi=150)
    print(f"saved figure -> {out_png}")
    summary = {
        "mode": "fig5",
        "dim": dim,
        "nx_width": nx,
        "cells_per_D": nx / 5.0,
        "cells": cells,
        "mu": mu,
        "rho_drop": rho_drop,
        "dsigma_dT": dsigma_dT,
        "gradT": gradT,
        "r": R,
        "v_YGB": v_YGB,
        "tau": tau,
        "t_r": t_r,
        "wall": wall,
        "drift_corrected": not wall,
        "samareh_ratio": samareh_ratio,
        "ratio_plateau": ratio_plateau,
        "t_plateau_tr": t_plateau_tr,
        "overshoot": overshoot,
        "ratio_final": ratio_final,
        "slope_per_tr": slope_per_tr,
        "u_lab_final": float(u_lab[-1]),
        "u_far_final": float(u_far[-1]),
        "rises": bool(y_centroid[-1] > y_centroid[0]),
        "t_end": float(times[-1]),
        "t_end_tr": float(times[-1] / t_r),
    }

else:  # fig7
    mu_b, sigma_T, gradT, _, U_r, t_r = reference_scales()
    NAS_TRYGGVASON_PEAK = 0.13
    t_star, U_star = times / t_r, u_lab / U_r
    peak_i = int(np.argmax(U_star))
    peak, t_peak = float(U_star[peak_i]), float(t_star[peak_i])
    in_tail = times >= times[-1] - t_r
    terminal = float(U_star[in_tail].mean())
    print(f"[fig7] nx={nx} ny={ny} nz={nz}  cells={cells}  cells/D={nx / (Wx / 1.0):.1f}  nvars={nvars}  color idx={c_idx}")
    print(f"mu_b={mu_b:.4f}  sigma_T={sigma_T:.5f}  gradT={gradT:.5f}  U_r={U_r:.5f}  t_r={t_r:.4f}  run={t_star[-1]:.2f} t_r\n")
    print(f"{'step':>8} {'t*':>7} {'y_drop':>9} {'U':>10} {'U*':>9}")
    for i in range(len(times)):
        print(f"{steps[i]:>8} {t_star[i]:>7.2f} {y_centroid[i]:>9.5f} {u_lab[i]:>10.6f} {U_star[i]:>9.4f}")
    print(f"\nrising (+y, toward hot top): {y_centroid[-1] > y_centroid[0]}")
    print(f"overshoot peak  U* = {peak:.4f}  at t* = {t_peak:.2f}   [Nas & Tryggvason / Samareh Fig 7 peak ~ {NAS_TRYGGVASON_PEAK:.2f}]")
    print(f"terminal (final t_r window)  U* = {terminal:.4f}")
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.plot(t_star, U_star, "o-", color="C0", ms=3.0, lw=1.2, label=f"MFC ({nx}/box width, {nx / (Wx / 1.0):.0f}/$D$)")
    ax.axhline(NAS_TRYGGVASON_PEAK, ls="--", color="C3", lw=1.3, label=rf"Nas & Tryggvason peak $\approx$ {NAS_TRYGGVASON_PEAK:.2f}")
    ax.plot(t_peak, peak, "*", color="k", ms=14, zorder=5, label=rf"MFC peak = {peak:.3f}")
    ax.set_xlabel(r"$t^* = t / t_r$  ($t_r = \mu_b / |\sigma_T \nabla T|$)")
    ax.set_ylabel(r"$U^* = U / U_r$  ($U_r = |\sigma_T \nabla T|\, r_0 / \mu_b$)")
    ax.set_title("Fig 7: finite-Ma thermocapillary migration (Re=5, Ma=20, Ca=0.01666)")
    ax.set_xlim(left=0.0)
    ax.set_ylim(bottom=min(0.0, 1.05 * U_star.min()))
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right", fontsize=9)
    fig.tight_layout()
    out_png = os.path.join(viz_dir, "fig7_migration.png")
    fig.savefig(out_png, dpi=150)
    print(f"saved figure -> {out_png}")
    summary = {
        "mode": "fig7",
        "nx_width": nx,
        "cells_per_D": nx / (Wx / 1.0),
        "cells": cells,
        "mu_b": mu_b,
        "sigma_T": sigma_T,
        "gradT": gradT,
        "U_r": U_r,
        "t_r": t_r,
        "peak": peak,
        "t_peak_tr": t_peak,
        "terminal": terminal,
        "nas_tryggvason_peak": NAS_TRYGGVASON_PEAK,
        "rises": bool(y_centroid[-1] > y_centroid[0]),
        "t_end_tr": float(t_star[-1]),
    }

print("RESULT_JSON " + json.dumps(summary))
