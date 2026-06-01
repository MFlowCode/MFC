#!/usr/bin/env python3
"""Benchmark MFC's single-material JWL shock tube against the exact Riemann solution."""
import os
import sys
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from jwl_exact_riemann import JWL, solve

CASE = os.path.join(os.path.dirname(__file__), "..",
                    "examples", "1D_jwl_single_material_shocktube")
D = os.path.join(CASE, "D")
OUT = os.path.join(os.path.dirname(__file__), "..", "jwl_verification_out")
os.makedirs(OUT, exist_ok=True)

# --- case parameters (must mirror case.py) ---
rho0 = 1770.0
jwl = JWL(A=3.712e11, B=3.231e9, R1=4.15, R2=0.95, omega=0.30, rho0=rho0)
x0 = 0.5
dt = 2.0e-8
# prim var indices: 1=rho, 2=u, 3=p
VARS = {"rho": 1, "u": 2, "p": 3}


def load(step):
    out = {}
    for name, idx in VARS.items():
        f = os.path.join(D, f"prim.{idx}.00.{step:06d}.dat")
        data = np.loadtxt(f)
        out["x"] = data[:, 0]
        out[name] = data[:, 1]
    return out


def steps_available():
    files = glob.glob(os.path.join(D, "prim.1.00.*.dat"))
    return sorted(int(os.path.basename(f).split(".")[3]) for f in files)


def benchmark_step(step):
    sim = load(step)
    t = step * dt
    x = sim["x"]
    WL = (sim["rho"][0], sim["u"][0], sim["p"][0])
    WR = (sim["rho"][-1], sim["u"][-1], sim["p"][-1])
    ex = solve(jwl, jwl, WL, WR, x0, t, x)

    # error norms (L1, L2) relative
    errs = {}
    for k in ("rho", "u", "p"):
        s, e = sim[k], ex[k]
        scale = np.max(np.abs(e)) + 1e-30
        l1 = np.mean(np.abs(s - e)) / scale
        l2 = np.sqrt(np.mean((s - e) ** 2)) / scale
        errs[k] = (l1, l2)
    return sim, ex, t, errs


def plot_step(step, sim, ex, t, errs, zoom=None):
    fig, axes = plt.subplots(3, 1, figsize=(9, 11), sharex=True)
    panels = [("rho", r"Density $\rho$ [kg/m$^3$]"),
              ("u", r"Velocity $u$ [m/s]"),
              ("p", r"Pressure $p$ [Pa]")]
    for ax, (k, lab) in zip(axes, panels):
        ax.plot(sim["x"], ex[k], "-", color="k", lw=2, label="Exact Riemann (same JWL EOS)")
        ax.plot(sim["x"], sim[k], "o", ms=3.5, mfc="none", color="crimson",
                label=f"MFC (MUSCL, HLLC, {len(sim['x'])} cells)")
        ax.set_ylabel(lab)
        ax.grid(alpha=0.3)
        l1, l2 = errs[k]
        ax.set_title(f"{k}: rel L1={l1:.2e}, rel L2={l2:.2e}", fontsize=9, loc="right")
        if zoom:
            ax.set_xlim(*zoom)
    axes[0].legend(loc="best", fontsize=9)
    axes[-1].set_xlabel("x [m]")
    wL = "shock" if ex["left_shock"] else "rarefaction"
    wR = "shock" if ex["right_shock"] else "rarefaction"
    fig.suptitle(f"Single-material JWL shock tube @ t={t*1e6:.3f} us  "
                 f"(left {wL} | contact | right {wR})\n"
                 f"p*={ex['pstar']/1e9:.3f} GPa, u*={ex['ustar']:.1f} m/s",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    tag = "zoom" if zoom else "full"
    path = os.path.join(OUT, f"jwl_step{step:06d}_{tag}.png")
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def main():
    steps = steps_available()
    print("steps:", steps)
    summary = []
    for step in steps:
        if step == 0:
            continue
        sim, ex, t, errs = benchmark_step(step)
        p1 = plot_step(step, sim, ex, t, errs)
        # determine wave extent for zoom
        span = 0.06 + 4000.0 * t  # heuristic widen with time
        p2 = plot_step(step, sim, ex, t, errs, zoom=(x0 - span, x0 + span))
        summary.append((step, t, errs))
        print(f"step {step:4d} t={t*1e6:6.3f}us  "
              f"rho L1={errs['rho'][0]:.2e}  u L1={errs['u'][0]:.2e}  "
              f"p L1={errs['p'][0]:.2e}  p*={ex['pstar']/1e9:.3f}GPa")

    # convergence-style summary plot of L1 vs time
    fig, ax = plt.subplots(figsize=(8, 5))
    ts = [s[1] * 1e6 for s in summary]
    for k, c in zip(("rho", "u", "p"), ("navy", "green", "crimson")):
        ax.plot(ts, [s[2][k][0] for s in summary], "o-", color=c, label=f"{k} rel L1")
    ax.set_xlabel("t [us]")
    ax.set_ylabel("relative L1 error")
    ax.set_yscale("log")
    ax.grid(alpha=0.3, which="both")
    ax.legend()
    ax.set_title("MFC vs exact JWL Riemann: error history")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "error_history.png"), dpi=130)
    plt.close(fig)
    print("\nWrote plots to", os.path.abspath(OUT))


if __name__ == "__main__":
    main()
