#!/usr/bin/env python3
"""
Mean gas velocity through the shock-curtain run of 2D_ibm_shock_curtain against a published profile.

    check_curtain.py CASE_DIR [--png profiles.png] [--ref sridhar2025_ms166_phi02.txt]
    check_curtain.py CASE_DIR --drag --mach 1.22        # a --single run: unsteady drag of one body

CASE_DIR is a run with --silo and post_process (an adaptive-step run, --steps 0, so that save k is at k t_save). For every save near a
reference panel t/tau_L = 0.5, 1, 2 (tau_L = L/W_s with t = 0 at the shock's arrival at the curtain) it prints the Favre-averaged
velocity u~(x) = <I rho u>/<I rho> over the cross-stream slabs (I = 1 for cell centres outside every body), in units of the
post-shock velocity, and the mean absolute difference from the reference profile over -0.35 < (x - x0)/L < 1, excluding +-0.05 L around
the shocks. The reference is the viscous particle-resolved DNS of Sridhar, Fox and Capecelatro (2025, Case 5, Ms 1.66,
phi 0.2); its values were read from the published figure and carry about 0.01 u_ps of reading error. In 2D the comparison is
qualitative (cylinders against spheres at equal solid fraction).

--drag reads D/ib_forces.dat of a --single run and prints Cd(t) = F_x/(1/2 rho2 u2^2 A) (A = pi D^2/4 for a sphere), its peak and the
time of the peak in t/tau (tau = D/W_s, t = 0 at the shock's arrival at the nose), and the mean over the last tau. The inviscid
single-sphere peak of Mehta et al. (2016, Fig. 6) is tabulated in MEHTA_PEAK.
"""

# Peak Cd and its time t/tau of the inviscid single sphere, Mehta et al. (2016), Fig. 6, read from the figure
MEHTA_PEAK = {1.22: (7.49, 0.36), 1.5: (3.72, 0.35), 2.0: (2.34, 0.34), 3.0: (1.71, 0.38), 6.0: (1.44, 0.34)}

import argparse
import csv
import glob
import math
import os
import re

import numpy as np


def read_inp(run):
    out = {}
    for name in ("pre_process.inp", "simulation.inp"):
        for line in open(os.path.join(run, name)):
            m = re.match(r"\s*([^=\s]+)\s*=\s*(\S+)", line)
            if m:
                out[m.group(1)] = m.group(2).strip("'\"")
    return out


def gas_mask(inp, xs):
    """I = 1 in gas: cell centres outside every body (periodic across the flow)."""
    f = lambda k: float(inp[k].replace("d", "e").replace("D", "e"))  # noqa: E731
    dim = 3 if int(inp["p"]) > 0 else 2
    n, R, W = int(inp["num_ibs"]), f("patch_ib(1)%radius"), f("y_domain%end")
    grid = np.meshgrid(*xs, indexing="ij")
    inside = np.zeros(grid[0].shape, bool)
    for i in range(n):
        c = [f(f"patch_ib({i + 1})%{k}_centroid") for k in "xyz"[:dim]]
        r2 = (grid[0] - c[0]) ** 2
        for g, ci in zip(grid[1:], c[1:]):
            r2 = r2 + np.minimum(np.abs(g - ci), W - np.abs(g - ci)) ** 2
        inside |= r2 < R * R
    return ~inside


def drag(a, inp):
    f = lambda k: float(inp[k].replace("d", "e").replace("D", "e"))  # noqa: E731
    D, g = 2 * f("patch_ib(1)%radius"), 1.4
    rho1, p1 = f("patch_icpp(2)%alpha_rho(1)"), f("patch_icpp(2)%pres")
    rho2, u2 = f("patch_icpp(1)%alpha_rho(1)"), f("patch_icpp(1)%vel(1)")
    Ws = a.mach * math.sqrt(g * p1 / rho1)
    area = math.pi * D**2 / 4 if int(inp["p"]) > 0 else D
    d = np.loadtxt(os.path.join(a.run, "D", "ib_forces.dat"), ndmin=2)
    d = d[np.argsort(d[:, 1])]
    t, cd = d[:, 1] * Ws / D, d[:, 2] / (0.5 * rho2 * u2**2 * area)
    i = int(np.argmax(cd))
    late = t >= t[-1] - 1.0
    print(f"Ms {a.mach:g}  {D / (f('x_domain%end') - f('x_domain%beg')) * (int(inp['m']) + 1):.1f} cells/D  t_end {t[-1]:.2f} tau")
    print(f"  Cd peak {cd[i]:.3f} at t/tau {t[i]:.3f}   mean Cd over the last tau {cd[late].mean():.3f}")
    if a.mach in MEHTA_PEAK:
        ref, tr = MEHTA_PEAK[a.mach]
        print(f"  Mehta et al. inviscid: peak {ref} at {tr}: {100 * (cd[i] / ref - 1):+.1f} %")


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("--png", default=None, help="Write the profiles at the save nearest each reference panel to this PNG.")
    ap.add_argument("--drag", action="store_true", help="Unsteady drag of a --single run.")
    ap.add_argument("--mach", type=float, default=1.66, help="Shock Mach number.")
    ap.add_argument("--h", type=float, default=None, help="Curtain thickness [D] (default: from the body positions).")
    ap.add_argument("--ref", default=os.path.join(here, "sridhar2025_ms166_phi02.txt"))
    a = ap.parse_args()
    inp = read_inp(a.run)
    if a.drag:
        return drag(a, inp)
    from mfc.viz.silo_reader import assemble_silo

    f = lambda k: float(inp[k].replace("d", "e").replace("D", "e"))  # noqa: E731
    ref = {}
    for row in csv.reader(l for l in open(a.ref) if not l.startswith("#")):
        if row[0] != "t_over_tauL":
            ref.setdefault(float(row[0]), []).append((float(row[1]), float(row[2])))
    ref = {t: np.array(v) for t, v in ref.items()}
    D, g = 2 * f("patch_ib(1)%radius"), 1.4
    rho1, p1, u_ps = f("patch_icpp(2)%alpha_rho(1)"), f("patch_icpp(2)%pres"), f("patch_icpp(1)%vel(1)")
    Ws = a.mach * math.sqrt(g * p1 / rho1)
    # Curtain thickness: the bodies fill [D/2, h - D/2] unless --h (in D) is given
    h = a.h * D if a.h else max(f(f"patch_ib({i + 1})%x_centroid") for i in range(int(inp["num_ibs"]))) + 0.5 * D
    tau, shift, t_save = h / Ws, 0.5 * D / Ws, f("t_save")
    steps = sorted(int(os.path.basename(p)[:-5]) for p in glob.glob(os.path.join(a.run, "silo_hdf5", "p0", "*.silo")) if os.path.basename(p)[:-5].isdigit())
    print(f"curtain {h / D:.1f} D thick, {inp['num_ibs']} bodies, u_ps {u_ps:.1f} m/s, tau_L {tau:.3e} s")
    best = {}
    print("| save | t/tau_L | ref panel | in-curtain mean u~/u_ps | mean abs diff from reference |")
    print("| ---: | ---: | ---: | ---: | ---: |")
    for k in steps:
        if k == 0:
            continue
        t = (k * t_save - shift) / tau
        tr = min(ref, key=lambda s: abs(s - t))
        if abs(tr - t) > 0.25:
            continue
        d = assemble_silo(a.run, k)
        xs = [np.asarray(c).squeeze() for c in (d.x_cc, d.y_cc, d.z_cc)[: d.ndim]]
        rho, rhou = (np.asarray(d.variables[v]).squeeze() for v in ("rho", "vel1"))
        mask = gas_mask(inp, xs)
        ax = tuple(range(1, d.ndim))
        u = (rho * rhou * mask).sum(ax) / (rho * mask).sum(ax) / u_ps
        xn = xs[0] / h
        rx, ru = ref[tr][:, 0], ref[tr][:, 1]
        ui = np.interp(rx, xn, u)
        keep = (rx > -0.35) & (rx < 1.0)
        for xj in rx[np.where(np.abs(np.diff(ru)) > 0.15)[0]]:
            keep &= np.abs(rx - xj) > 0.05
        inc = (xn > 0.1) & (xn < 0.9)
        print(f"| {k} | {t:.3f} | {tr:g} | {u[inc].mean():.3f} | {np.mean(np.abs(ui[keep] - ru[keep])):.3f} |")
        if tr not in best or abs(tr - t) < best[tr][0]:
            best[tr] = (abs(tr - t), t, xn, u)
    if a.png:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axs = plt.subplots(1, len(ref), figsize=(9, 3), sharey=True, constrained_layout=True)
        for ax, tr in zip(np.atleast_1d(axs), sorted(ref)):
            ax.axvspan(0, 1, color="0.92")
            ax.plot(ref[tr][:, 0], ref[tr][:, 1], "o", mfc="none", mec="0.4", ms=3, label="Sridhar et al. (2025), DNS")
            if tr in best:
                ax.plot(best[tr][2], best[tr][3], "-", color="C0", lw=1.5, label=f"MFC, t/$\\tau_L$ = {best[tr][1]:.2f}")
            ax.set(xlim=(-0.4, 1.35), ylim=(-0.05, 1.5), xlabel="(x - x0) / L", title=f"t/$\\tau_L$ = {tr:g}")
            ax.legend(fontsize=7, frameon=False, loc="upper left")
        np.atleast_1d(axs)[0].set_ylabel("Favre-averaged $u$ / $u_{ps}$")
        fig.savefig(a.png, dpi=110)


if __name__ == "__main__":
    main()
