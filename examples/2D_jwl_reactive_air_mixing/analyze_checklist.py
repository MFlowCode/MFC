#!/usr/bin/env python3
"""Reactive-JWL validation checklist for the 2D PETN products-air mixing run.

Runs four quantitative checks directly from the D/ field output, by binning the
(nearly axisymmetric) 2D field onto radius r = hypot(x, y):

  checklist_1_propagation.png -- radial profiles p, rho, u_r, lambda, Y at
                                 successive times (initiation and propagation).
  checklist_2_energy.png      -- volume-integrated explosive mass (conserved),
                                 unreacted mass rho Y (1-lambda) (must fall to
                                 zero), and total energy (energy-release check).
  checklist_3_breakout.png    -- shock trajectory R_s(t) with breakout, the
                                 peak-pressure envelope p_shock(r) with a fitted
                                 decay slope, and specific impulse I(r).
  checklist_4_closure.png     -- the mixing/closure region 1e-6 < Y < 1-1e-6
                                 and the (p, rho, e) thermodynamic state there.

House style: Agg backend, perceptually-uniform colormaps, sim as markers where a
reference line exists, explicit units, constrained_layout.

Usage: ./analyze_checklist.py [case_dir]
"""

import glob
import os
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

HERE = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "D")

DT_STEP = 1.5e-8  # s per solver step
SAVE_EVERY = 250  # steps between dumps
P_AMB = 101325.0  # Pa
RHO0 = 1000.0  # kg/m^3, PETN_KUHL reference density
CHARGE_R = 0.05  # m
DX = 0.4 / 400.0  # m, uniform cell size
CELL_A = DX * DX  # m^2 (per unit depth)
# Self-sustained CJ state of the standard-form JWL with Kuhl's PETN constants,
# from the Hugoniot/Rayleigh tangency (rate-independent, EOS only). The same
# construction reproduces literature CJ for TNT/PETN/Comp-B to a few percent.
# Kuhl reports only T_CJ ~ 4600 K, no P_CJ/D_CJ; these are computed, not his.
D_CJ = 5384.0  # m/s
P_CJ = 7.35e9  # Pa


def load(pre, i, t):
    """Concatenate every rank's chunk (parallel_io=F) -> (x, y, value) columns."""
    files = sorted(glob.glob(f"{D}/{pre}.{i}.*.{t:06d}.dat"))
    a = np.vstack([np.loadtxt(f) for f in files])
    return a[:, 0], a[:, 1], a[:, 2]


def frame(t):
    """All fields at dump t as flat cell arrays, plus radius r."""
    x, y, p = load("prim", 5, t)
    _, _, u = load("prim", 3, t)
    _, _, v = load("prim", 4, t)
    _, _, lam = load("prim", 8, t)
    _, _, ar1 = load("cons", 1, t)  # alpha_rho products
    _, _, ar2 = load("cons", 2, t)  # alpha_rho air
    _, _, Etot = load("cons", 5, t)  # total energy density rho*E
    rho = ar1 + ar2
    r = np.hypot(x, y)
    ur = np.where(r > 0, (u * x + v * y) / np.maximum(r, 1e-30), 0.0)
    Y = ar1 / np.maximum(rho, 1e-30)
    e = (Etot - 0.5 * rho * (u * u + v * v)) / np.maximum(rho, 1e-30)
    return dict(r=r, p=p, rho=rho, ur=ur, lam=lam, Y=Y, ar1=ar1, Etot=Etot, e=e, t=t)


def radial_mean(r, q, edges):
    idx = np.digitize(r, edges)
    nb = len(edges) - 1
    out = np.full(nb, np.nan)
    for b in range(1, nb + 1):
        m = idx == b
        if m.any():
            out[b - 1] = q[m].mean()
    return out


steps = sorted(int(f.rsplit(".", 2)[1]) for f in glob.glob(f"{D}/prim.5.00.*.dat"))
if not steps:
    sys.exit(f"no prim output found in {D} -- run the case first")
times = np.array(steps) * DT_STEP

GRAY, MUTED = "#444444", "#777777"
plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.edgecolor": MUTED, "axes.linewidth": 0.6, "xtick.color": GRAY, "ytick.color": GRAY, "text.color": GRAY, "axes.labelcolor": GRAY})

edges = np.linspace(0.0, 0.28, 141)  # ~2 mm radial bins out to the domain corner
rc = 0.5 * (edges[:-1] + edges[1:])

# Load every frame once; radial profiles + integrals reused across figures.
frames = [frame(t) for t in steps]
sel = np.linspace(0, len(frames) - 1, 6).astype(int)  # 6 times for profile panels
cmap = matplotlib.colormaps["viridis"]

# Shock front R_s(t) and the time the blast first touches the (open) domain edge.
# Beyond that time, volume integrals are not closed -- mass and energy leave
# through the non-reflecting boundary, so conservation checks apply only before it.
EDGE_R = 0.19  # m, just inside the +/-0.2 m domain edge
Rs = np.array([f["r"][f["p"] > 2 * P_AMB].max() if (f["p"] > 2 * P_AMB).any() else 0.0 for f in frames])
i_edge = int(np.argmax(Rs > EDGE_R)) if (Rs > EDGE_R).any() else len(frames) - 1
t_edge = times[i_edge]


def checklist_1():
    fig, ax = plt.subplots(5, 1, figsize=(7.2, 12.0), sharex=True, constrained_layout=True)
    fig.suptitle("Radial profiles of the diverging PETN detonation (2D, 400×400, dx = 1 mm)", fontsize=11, color="#222222")
    keys = [("p", "p [Pa]  (log)", True), ("rho", "ρ [kg/m³]  (log)", True), ("ur", "u_r [m/s]", False), ("lam", "λ  (reaction progress)", False), ("Y", "Y  (products mass fraction)", False)]
    for row, (k, ylabel, logy) in enumerate(keys):
        for j in sel:
            f = frames[j]
            prof = radial_mean(f["r"], f[k], edges)
            c = cmap(j / (len(frames) - 1))
            ax[row].plot(rc * 1000, prof, color=c, lw=1.3, label=f"t = {f['t'] * DT_STEP * 1e6:.0f} µs" if row == 0 else None)
        ax[row].axvline(CHARGE_R * 1000, color=MUTED, ls=":", lw=0.8)
        ax[row].set_ylabel(ylabel)
        if logy:
            ax[row].set_yscale("log")
    ax[0].axhline(P_CJ, color="crimson", ls="--", lw=0.8, label="EOS CJ pressure ≈ 7.4 GPa")
    ax[0].text(CHARGE_R * 1000 + 2, P_CJ * 1.3, "charge edge", color=MUTED, fontsize=7)
    ax[0].legend(fontsize=7, ncol=2, loc="upper right")
    ax[-1].set_xlabel("radius  r [mm]")
    out = os.path.join(HERE, "checklist_1_propagation.png")
    fig.savefig(out, dpi=140)
    print(f"wrote {out}")


def checklist_2():
    m_expl = np.array([f["ar1"].sum() * CELL_A for f in frames])  # conserved
    m_unrx = np.array([(f["ar1"] * (1 - f["lam"])).sum() * CELL_A for f in frames])
    e_tot = np.array([f["Etot"].sum() * CELL_A for f in frames])  # total energy
    # conservation applies only before the blast reaches the open boundary
    pre = slice(0, i_edge + 1)
    m_dev = 100 * (m_expl[pre].max() - m_expl[pre].min()) / m_expl[0]
    fig, ax = plt.subplots(1, 2, figsize=(11.0, 4.2), constrained_layout=True)
    fig.suptitle("Explosive consumption and energy release", fontsize=11, color="#222222")
    for a in ax:
        a.axvspan(t_edge * 1e6, times[-1] * 1e6, color=MUTED, alpha=0.12)
        a.axvline(t_edge * 1e6, color=MUTED, ls="--", lw=0.8)
    ax[0].text(t_edge * 1e6 + 0.5, 0.45, "blast reaches open\nboundary (efflux)", color=MUTED, fontsize=7)
    ax[0].plot(times * 1e6, m_expl / m_expl[0], "o-", color="#2E6B69", ms=4, label="explosive mass  ∫ρ_expl dV  (flat until efflux)")
    ax[0].plot(times * 1e6, m_unrx / m_unrx[0], "s-", color="#B8420F", ms=4, label="unreacted mass  ∫ρY(1-λ) dV  (falls to zero)")
    ax[0].set_xlabel("time  t [µs]")
    ax[0].set_ylabel("fraction of initial")
    ax[0].axhline(0, color=MUTED, ls=":", lw=0.6)
    ax[0].legend(fontsize=7)
    ax[0].set_title(f"explosive mass conserved to {m_dev:.3f}% before efflux; " f"unreacted fraction {100 * m_unrx[i_edge] / m_unrx[0]:.1f}%")
    ax[1].plot(times * 1e6, e_tot / e_tot[0], "^-", color="#154c79", ms=4, label="total energy  ∫ρE dV")
    ax[1].set_xlabel("time  t [µs]")
    ax[1].set_ylabel("∫ρE dV / initial")
    ax[1].legend(fontsize=7)
    e_plateau = 100 * (e_tot[3 : i_edge + 1].max() - e_tot[3 : i_edge + 1].min()) / e_tot[3]
    ax[1].set_title(f"chemical release ×{e_tot[i_edge] / e_tot[0]:.2f} embedded in EOS, " f"then flat to {e_plateau:.2f}%")
    out = os.path.join(HERE, "checklist_2_energy.png")
    fig.savefig(out, dpi=140)
    print(f"wrote {out}")
    return m_expl, m_unrx, e_tot


def checklist_3():
    # shock front R_s(t) is the module-level Rs (p > 2 x ambient front radius)
    # peak-pressure envelope over radius: max p in each radial bin over all frames
    penv = np.full(len(rc), np.nan)
    for f in frames:
        pr = radial_mean(f["r"], f["p"], edges)
        penv = np.fmax(penv, pr)
    # specific impulse I(r) = integral of (p - p_amb) dt at each radius
    P_of_rt = np.array([radial_mean(f["r"], f["p"], edges) for f in frames])  # (nt, nr)
    over = np.clip(P_of_rt - P_AMB, 0, None)
    I_r = np.trapz(over, times, axis=0) if hasattr(np, "trapz") else np.trapezoid(over, times, axis=0)

    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.3), constrained_layout=True)
    fig.suptitle("Blast breakout and near-field decay", fontsize=11, color="#222222")

    ax[0].plot(times * 1e6, Rs * 1000, "o", mfc="none", color="#154c79", ms=6, label="MFC shock front (p > 2 atm)")
    tt = np.linspace(times.min(), times.max(), 100)
    # reference slope: charge-edge breakout ray at the self-sustained planar CJ
    # speed; a diverging front runs below it (curvature velocity deficit)
    ib = np.argmax(Rs > CHARGE_R)
    t_bo = times[ib]
    ax[0].axhline(CHARGE_R * 1000, color=MUTED, ls=":", lw=0.8)
    ax[0].axvline(t_bo * 1e6, color="crimson", ls="--", lw=0.8, label=f"breakout t ≈ {t_bo * 1e6:.0f} µs")
    ax[0].plot(tt * 1e6, (CHARGE_R + D_CJ * (tt - t_bo)) * 1000, color="crimson", lw=1.0, alpha=0.6, label="planar CJ ray (D ≈ 5.4 km/s)")
    ax[0].set_xlabel("time  t [µs]")
    ax[0].set_ylabel("shock radius  R_s [mm]")
    ax[0].legend(fontsize=7, loc="lower right")
    ax[0].set_title("front trajectory and breakout")

    # fit window: beyond the charge edge, inside the boundary-contaminated far
    # field (skill plotting-validation §1.2 -- fit only where the law holds)
    R_FIT_HI = 0.18
    show = np.isfinite(penv) & (penv > 2 * P_AMB) & (rc > CHARGE_R)
    fitw = show & (rc < R_FIT_HI)
    ax[1].loglog(rc[show] * 1000, penv[show], "o", mfc="none", color="#B8420F", ms=5, label="MFC peak-pressure envelope")
    if fitw.sum() > 3:
        b, a = np.polyfit(np.log(rc[fitw]), np.log(penv[fitw]), 1)
        rrf = rc[fitw]
        ax[1].loglog(rrf * 1000, np.exp(a) * rrf**b, color="#B8420F", lw=1.4, label=f"fit  p ∝ r^{b:.2f}  (fit window)")
    ax[1].axvspan(CHARGE_R * 1000, R_FIT_HI * 1000, color="#B8420F", alpha=0.07)
    ax[1].set_xlabel("radius  r [mm]")
    ax[1].set_ylabel("peak pressure  p_shock [Pa]")
    ax[1].legend(fontsize=7)
    ax[1].set_title("peak-pressure decay (fit charge edge to 180 mm)")

    ax[2].plot(rc * 1000, I_r, "-", color="#2E6B69", lw=1.4)
    ax[2].axvline(CHARGE_R * 1000, color=MUTED, ls=":", lw=0.8)
    ax[2].set_xlabel("radius  r [mm]")
    ax[2].set_ylabel("specific impulse  I = ∫(p−p_amb) dt  [Pa·s]")
    ax[2].set_title("radial impulse distribution")
    out = os.path.join(HERE, "checklist_3_breakout.png")
    fig.savefig(out, dpi=140)
    print(f"wrote {out}")
    return Rs, t_bo


def checklist_4():
    # pick a mid/late frame with a developed mixing layer
    f = frames[len(frames) * 3 // 4]
    mix = (f["Y"] > 1e-6) & (f["Y"] < 1 - 1e-6)
    fig, ax = plt.subplots(1, 3, figsize=(15.0, 4.4), constrained_layout=True)
    fig.suptitle(f"Products-air mixing region: cells with 1e-6 < Y < 1-1e-6  " f"(t = {f['t'] * DT_STEP * 1e6:.0f} µs, {mix.sum()} cells)", fontsize=11, color="#222222")
    # spatial location of the closure cells, colored by Y
    x, y, _ = load("prim", 5, f["t"])
    sc = ax[0].scatter(x[mix], y[mix], c=f["Y"][mix], cmap="viridis", s=2, vmin=0, vmax=1)
    ax[0].set_aspect("equal")
    ax[0].set_xlabel("x [m]")
    ax[0].set_ylabel("y [m]")
    ax[0].set_title("where the mixture closure is active")
    fig.colorbar(sc, ax=ax[0], shrink=0.85, label="Y")
    # thermodynamic state in the closure cells: p vs rho colored by Y
    sc2 = ax[1].scatter(f["rho"][mix], f["p"][mix], c=f["Y"][mix], cmap="viridis", s=3, vmin=0, vmax=1)
    ax[1].set_yscale("log")
    ax[1].set_xscale("log")
    ax[1].set_xlabel("ρ [kg/m³]")
    ax[1].set_ylabel("p [Pa]")
    ax[1].set_title("(ρ, p) state in the closure region")
    fig.colorbar(sc2, ax=ax[1], shrink=0.85, label="Y")
    # specific internal energy vs Y in the closure cells
    ax[2].scatter(f["Y"][mix], f["e"][mix], s=3, color="#154c79", alpha=0.5)
    ax[2].set_xlabel("Y  (products mass fraction)")
    ax[2].set_ylabel("specific internal energy  e [J/kg]")
    ax[2].set_title("e across the mixing layer")
    out = os.path.join(HERE, "checklist_4_closure.png")
    fig.savefig(out, dpi=140)
    print(f"wrote {out}")
    return mix.sum()


checklist_1()
m_expl, m_unrx, e_tot = checklist_2()
Rs, t_bo = checklist_3()
nmix = checklist_4()

print("\n=== checklist summary (integrals valid before boundary efflux) ===")
print(f"blast reaches open boundary at t ~ {t_edge * 1e6:.0f} us (frame {i_edge})")
print(f"explosive mass conserved to " f"{100 * (m_expl[:i_edge + 1].max() - m_expl[:i_edge + 1].min()) / m_expl[0]:.3f}% " f"before efflux")
print(f"unreacted mass  {m_unrx[0]:.3e} -> {m_unrx[i_edge]:.3e} kg  " f"({100 * m_unrx[i_edge] / m_unrx[0]:.1f}% remains at efflux time)")
print(f"total energy rises x{e_tot[i_edge] / e_tot[0]:.2f} (chemical release, EOS-embedded)")
print(f"breakout of charge at t ~ {t_bo * 1e6:.1f} us; " f"front at {Rs[i_edge] * 1000:.0f} mm when it reaches the boundary")
print(f"closure cells active in shown frame: {nmix}")
