#!/usr/bin/env python3
"""Compare MFC's variable-viscosity Couette runs to the exact solution.

This script does NOT run MFC. Run the cases first with run_suite.sh (which leaves
each grid in runs/n<N>/), then run this to analyze them:

    ./run_suite.sh 32 64 96
    python3 validate.py            # analyzes runs/n32, runs/n64, runs/n96

For each grid it reads back the steady velocity and temperature profiles MFC
saved, compares them to the exact coupled-BVP solution from reference.py, reports
the error, and -- across grids -- the observed spatial order of accuracy. It
writes figures/ and summary.json.

The headline check is the velocity profile: mu(T) = exp(C + D/T) makes the hot
(thin) fluid near the moving wall shear faster than the cold (viscous) fluid near
the fixed wall, so u(y) is curved. A constant-viscosity solver would give a
straight line; the curvature is the mu(T) signal, and it must match the exact
profile to the discretization error.
"""

import glob
import json
import os
import re
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import couette_config as cfg
import matplotlib.pyplot as plt
import reference

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "runs")
FIG = os.path.join(HERE, "figures")
SUMMARY = os.path.join(HERE, "summary.json")

# model_eqns = 3, num_fluids = 1, 2D conserved-field layout (0-based):
#   0: rho   1: rho*u   2: rho*v   3: E   4: alpha   5: rho*e (internal energy)
I_RHO, I_MOMX, I_RHOE = 0, 1, 5


def read_inp(rundir):
    inp = {}
    for line in open(os.path.join(rundir, "simulation.inp")):
        if "=" in line:
            k, v = line.split("=", 1)
            inp[k.strip().lower()] = v.strip().rstrip(",")
    return inp


def steps_of(rundir):
    rd = os.path.join(rundir, "restart_data")
    steps = []
    for f in glob.glob(os.path.join(rd, "lustre_*.dat")):
        match = re.search(r"lustre_(\d+)\.dat$", os.path.basename(f))
        if match:
            steps.append(int(match.group(1)))
    return sorted(steps)


def load_step(rundir, m, n, step):
    a = np.fromfile(os.path.join(rundir, "restart_data", f"lustre_{step}.dat"), np.float64)
    return a.reshape(a.size // (m * n), 1, n, m)  # (nvars, p=1, n, m), x innermost


def profiles(fields):
    """x-averaged wall-normal profiles u(y) and T(y) from the conserved fields."""
    rho = fields[I_RHO, 0]  # (n, m)
    u = fields[I_MOMX, 0] / rho
    T = fields[I_RHOE, 0] / (rho * cfg.cv)  # ideal gas (p_inf = 0): T = e/cv
    return u.mean(axis=1), T.mean(axis=1)  # average over x -> (n,)


def analyze_grid(N):
    rundir = os.path.join(RUNS, f"n{N}")
    if not os.path.isdir(os.path.join(rundir, "restart_data")):
        raise SystemExit(f"missing {rundir}/restart_data -- run ./run_suite.sh {N} first")
    inp = read_inp(rundir)
    m, n = int(inp["m"]) + 1, int(inp["n"]) + 1
    yc = (np.arange(n) + 0.5) * cfg.H / n
    steps = steps_of(rundir)

    u_last, T_last = profiles(load_step(rundir, m, n, steps[-1]))
    u_prev, _ = profiles(load_step(rundir, m, n, steps[-2]))
    # Steadiness: how much the velocity profile still moves between the last two
    # saved snapshots, relative to the wall speed.
    unsteadiness = float(np.max(np.abs(u_last - u_prev)) / cfg.U)

    u_ex = reference.u_at(yc)
    T_ex = reference.T_at(yc)
    err_u = float(np.sqrt(np.mean((u_last - u_ex) ** 2)) / cfg.U)
    err_u_inf = float(np.max(np.abs(u_last - u_ex)) / cfg.U)
    err_T = float(np.sqrt(np.mean((T_last - T_ex) ** 2)) / (cfg.T1 - cfg.T0))
    return {
        "N": N,
        "ny": n,
        "dy": cfg.H / n,
        "nsteps": steps[-1],
        "unsteadiness": unsteadiness,
        "err_u_L2": err_u,
        "err_u_Linf": err_u_inf,
        "err_T_L2": err_T,
        "yc": yc,
        "u": u_last,
        "T": T_last,
        "u_ex": u_ex,
        "T_ex": T_ex,
    }


def main(grids):
    os.makedirs(FIG, exist_ok=True)
    results = [analyze_grid(N) for N in grids]
    results.sort(key=lambda r: r["ny"])

    nys = np.array([r["ny"] for r in results])
    errs_u = np.array([r["err_u_L2"] for r in results])
    errs_T = np.array([r["err_T_L2"] for r in results])
    # error ~ dy^order; dy = H/Ny, so fit log(err) vs log(dy) -> positive order
    dys = cfg.H / nys
    order_u = float(np.polyfit(np.log(dys), np.log(errs_u), 1)[0])

    # --- console report ---
    print(f"\nVariable-viscosity Couette validation  (Re={cfg.Re:.0f}, Ma={cfg.Ma:.2f}, Br={cfg.Br:.3f}, mu-contrast={cfg.mu_of_T(cfg.T0) / cfg.mu_of_T(cfg.T1):.2f})")
    print(f"{'Ny':>5} {'dy':>10} {'unsteady':>11} {'L2(u)/U':>11} {'Linf(u)/U':>11} {'L2(T)/dT':>11}")
    for r in results:
        print(f"{r['ny']:>5} {r['dy']:>10.3e} {r['unsteadiness']:>11.2e} {r['err_u_L2']:>11.3e} {r['err_u_Linf']:>11.3e} {r['err_T_L2']:>11.3e}")
    print(f"observed spatial order (L2 velocity) = {order_u:.2f}")

    # --- figure 1: profile overlay on the finest grid ---
    rf = results[-1]
    fig, (axu, axt) = plt.subplots(1, 2, figsize=(11, 4.6))
    axu.plot(rf["u_ex"], rf["yc"], "-", color="k", lw=2, label="exact (coupled BVP)")
    axu.plot(rf["u"], rf["yc"], "o", ms=4, color="C0", label=f"MFC (Ny={rf['ny']})")
    axu.plot(cfg.U * rf["yc"] / cfg.H, rf["yc"], "--", color="C3", lw=1, label=r"constant-$\mu$ (straight)")
    axu.set_xlabel("u / U"), axu.set_ylabel("y / H")
    axu.set_title("velocity: curvature is the mu(T) signal"), axu.legend(fontsize=8)
    axt.plot(rf["T_ex"], rf["yc"], "-", color="k", lw=2, label="exact")
    axt.plot(rf["T"], rf["yc"], "s", ms=4, color="C1", label="MFC")
    axt.set_xlabel("T [K]"), axt.set_ylabel("y / H")
    axt.set_title("temperature profile"), axt.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "couette_profiles.png"), dpi=130)

    # --- figure 2: convergence ---
    fig2, ax = plt.subplots(figsize=(6, 5))
    ax.loglog(nys, errs_u, "o-", label=f"L2 velocity (order {order_u:.2f})")
    ax.loglog(nys, errs_T, "s-", label="L2 temperature")
    ref2 = errs_u[0] * (nys[0] / nys) ** 2
    ax.loglog(nys, ref2, "k--", lw=1, label="2nd-order reference")
    ax.set_xlabel("Ny (wall-normal cells)"), ax.set_ylabel("relative L2 error")
    ax.set_title("grid convergence, mu(T) Couette"), ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    fig2.tight_layout()
    fig2.savefig(os.path.join(FIG, "couette_convergence.png"), dpi=130)

    # --- figure 3: spatial relative error ---
    # Left: signed relative velocity error along the gap, one curve per grid --
    # shows where the error lives and that it shrinks under refinement. Right: the
    # same error scaled by (Ny/Ny_coarsest)^2; the curves collapse onto one shape,
    # so the 2nd-order convergence holds pointwise, not just in the L2 norm.
    ny0 = nys[0]
    fig3, (axe, axs) = plt.subplots(1, 2, figsize=(11, 4.6))
    for r in results:
        rel = 1e4 * (r["u"] - r["u_ex"]) / cfg.U  # in units of 1e-4 for readable ticks
        axe.plot(rel, r["yc"], "-o", ms=3, label=f"Ny={r['ny']}")
        axs.plot(rel * (r["ny"] / ny0) ** 2, r["yc"], "-o", ms=3, label=f"Ny={r['ny']}")
    for a in (axe, axs):
        a.axvline(0.0, color="k", lw=0.6)
        a.set_ylabel("y / H")
        a.legend(fontsize=8)
        a.xaxis.set_major_locator(plt.MaxNLocator(5))
    axe.set_xlabel(r"$(u_{\rm MFC} - u_{\rm exact})\,/\,U\quad[\times 10^{-4}]$")
    axe.set_title("relative velocity error")
    axs.set_xlabel(r"$(u_{\rm MFC} - u_{\rm exact})/U \times (N_y/N_{y,0})^2\quad[\times 10^{-4}]$")
    axs.set_title("grid-scaled error (collapse = 2nd order in space)")
    fig3.tight_layout()
    fig3.savefig(os.path.join(FIG, "couette_error.png"), dpi=130)

    # --- summary.json ---
    summary = {
        "config": {
            "Re": cfg.Re,
            "Ma": cfg.Ma,
            "Br": cfg.Br,
            "Pr": cfg.Pr,
            "mu_contrast": cfg.mu_of_T(cfg.T0) / cfg.mu_of_T(cfg.T1),
            "T0": cfg.T0,
            "T1": cfg.T1,
            "U": cfg.U,
            "u_mid_exact": float(reference.u_at(np.array([cfg.H / 2]))[0]),
            "u_mid_constant_mu": cfg.U / 2,
        },
        "observed_order_u_L2": order_u,
        "grids": [{k: r[k] for k in ("N", "ny", "dy", "nsteps", "unsteadiness", "err_u_L2", "err_u_Linf", "err_T_L2")} for r in results],
    }
    json.dump(summary, open(SUMMARY, "w"), indent=2)
    print(f"\nwrote {SUMMARY}")
    print(f"wrote {FIG}/couette_profiles.png, couette_convergence.png, couette_error.png")


if __name__ == "__main__":
    grids = [int(a) for a in sys.argv[1:]] or [32, 64, 96]
    main(grids)
