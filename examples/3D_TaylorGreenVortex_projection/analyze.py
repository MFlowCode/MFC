#!/usr/bin/env python3
"""
Post-processing for run_sweep.sh.

    analyze.py timing CSV... [--plot FILE] cost per step and wall time per convective time, with the
                                           speedup over explicit at the same grid, GPU count and Mach
    analyze.py ke RUN_DIR... [--plot FILE] kinetic energy E_k/E_k(0), dissipation rate -dE_k/dt
                                           (in U0^3/L) and enstrophy over t/tC, read from restart_data; the first
                                           run is the reference the others are compared against (and
                                           divided into in the plot's panel (b))
    analyze.py weak CSV [--plot FILE]      weak-scaling efficiency, s/step on one GPU over s/step on G at
                                           the same cells per GPU (run_sweep.sh weak)
"""

import argparse
import csv
import glob
import os
import re

import numpy as np

U0 = 0.1 * (1.4 * 101325.0) ** 0.5  # the velocity scale of hardcoded IC 380; L = rho0 = 1

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
sub = parser.add_subparsers(dest="mode", required=True)
timing_p = sub.add_parser("timing")
timing_p.add_argument("csv", nargs="+")
timing_p.add_argument("--plot", help="write wall time, speedup and scaling plots to this file")
weak_p = sub.add_parser("weak")
weak_p.add_argument("csv")
weak_p.add_argument("--plot", help="write the efficiency plot to this file")
ke_p = sub.add_parser("ke")
ke_p.add_argument("dirs", nargs="+")
ke_p.add_argument("--plot", help="write E_k and dissipation-rate curves to this file")
args = parser.parse_args()


def pyplot():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def timing(paths, plot):
    rows = []
    for path in paths:
        with open(path) as f:
            rows += [r for r in csv.DictReader(f) if r["s_per_step"]]
    for r in rows:
        r.update(N=int(r["N"]), ngpu=int(r["ngpu"]), mach=float(r["mach"]), sps=float(r["s_per_step"]), spt=float(r["steps_per_tC"]))
    # Explicit cost per step does not depend on Mach; its steps per tC are the projection's times 1 + 1/M
    explicit = {}
    for r in rows:
        if r["solver"] == "explicit":
            explicit.setdefault((r["N"], r["ngpu"]), []).append(r["sps"])
    explicit = {k: np.mean(v) for k, v in explicit.items()}
    proj = [r for r in rows if r["solver"] == "projection"]
    for r in proj:
        r["expl_tc"] = explicit.get((r["N"], r["ngpu"]), np.nan) * r["spt"] * (1 + 1 / r["mach"])
    print(f"{'solver':<10} {'N':>4} {'GPUs':>4} {'Mach':>7} {'s/step':>9} {'steps/tC':>9} {'wall/tC [s]':>11} {'expl/tC [s]':>11} {'speedup':>8}")
    for r in sorted(rows, key=lambda r: (r["N"], r["ngpu"], r["solver"], -r["mach"])):
        per_tc = r["sps"] * r["spt"]
        cols = f"{r['expl_tc']:11.3e} {r['expl_tc'] / per_tc:8.1f}" if "expl_tc" in r else f"{'-':>11} {'-':>8}"
        print(f"{r['solver']:<10} {r['N']:>4} {r['ngpu']:>4} {r['mach']:7.3g} {r['sps']:9.3e} {r['spt']:9.1f} {per_tc:11.3e} {cols}")
    if plot:
        timing_plot(proj, explicit, plot)


def timing_plot(proj, explicit, plot):
    """(a) wall time per tC and (b) speedup against Mach on 1 GPU, (c) strong scaling of the cost per step."""
    plt = pyplot()
    Ns = sorted({r["N"] for r in proj})
    gpus = sorted({r["ngpu"] for r in proj})
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    M = np.logspace(-3, -1, 50)
    for i, N in enumerate(Ns):
        c = f"C{i}"
        rs = sorted((r for r in proj if r["N"] == N and r["ngpu"] == 1), key=lambda r: r["mach"])
        if not rs:
            continue
        m = np.array([r["mach"] for r in rs])
        ax[0].loglog(m, [r["sps"] * r["spt"] for r in rs], "--o", color=c, mfc="none", label=f"projection, ${N}^3$")
        ax[0].loglog(M, explicit[(N, 1)] * rs[0]["spt"] * (1 + 1 / M), "-", color=c, label=f"explicit, ${N}^3$ (cost per step x steps)")
        sp = np.array([r["expl_tc"] / (r["sps"] * r["spt"]) for r in rs])
        ax[1].loglog(m, sp, "--o", color=c, mfc="none", label=f"${N}^3$, 1 GPU")
        for k, (x, y) in enumerate(zip(m, sp)):
            ax[1].annotate(f"{y:.0f}x", (x, y), textcoords="offset points", xytext=(0, 7 if i else -13), ha="center", fontsize=8, color=c)
        for g, a in zip((g for g in gpus if g > 1), (0.45, 0.25)):
            rg = sorted((r for r in proj if r["N"] == N and r["ngpu"] == g), key=lambda r: r["mach"])
            ax[1].loglog([r["mach"] for r in rg], [r["expl_tc"] / (r["sps"] * r["spt"]) for r in rg], ":", color=c, alpha=a + 0.3, label=f"${N}^3$, {g} GPUs")
        # Strong scaling at the middle Mach: cost per step relative to 1 GPU
        mid = rs[len(rs) // 2]["mach"]
        for solver, cost, ls in (
            ("projection", {r["ngpu"]: r["sps"] for r in proj if r["N"] == N and r["mach"] == mid}, "--o"),
            ("explicit", {g: explicit[(N, g)] for g in gpus if (N, g) in explicit}, "-s"),
        ):
            gs = sorted(cost)
            ax[2].plot(gs, [cost[1] / cost[g] for g in gs], ls, color=c, mfc="none", label=f"{solver}, ${N}^3$" + (f" (M = {mid:g})" if solver == "projection" else ""))
    ax[1].axhline(1, color="k", lw=0.8)
    ax[2].plot(gpus, gpus, "k:", lw=1, label="ideal")
    ax[0].set_title("(a) wall time per convective time, 1 GPU", loc="left")
    ax[0].set_ylabel("wall time per $t_C$ [s]")
    ax[1].set_title("(b) projection speedup (dotted: 2 and 4 GPUs)", loc="left")
    ax[1].set_ylim(top=ax[1].get_ylim()[1] * 2)
    ax[1].set_ylabel("explicit wall time / projection wall time")
    ax[2].set_title("(c) strong scaling of the cost per step", loc="left")
    ax[2].set_ylabel("speedup over 1 GPU")
    ax[2].set_xlabel("GPUs")
    ax[2].set_xticks(gpus)
    for a in ax[:2]:
        a.set_xlabel("Mach number")
        a.invert_xaxis()
    for a in ax:
        a.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(plot, dpi=150)


def weak(path, plot):
    with open(path) as f:
        rows = [r for r in csv.DictReader(f) if r["s_per_step"]]
    runs = {}
    for r in rows:
        runs.setdefault((r["solver"], int(r["N"]), float(r["mach"])), {})[int(r["ngpu"])] = float(r["s_per_step"])
    print(f"{'solver':<10} {'N/GPU':>5} {'Mach':>7} {'GPUs':>4} {'s/step':>9} {'efficiency':>10}")
    for (solver, N, mach), t in sorted(runs.items()):
        for g in sorted(t):
            print(f"{solver:<10} {N:>5} {mach:7.3g} {g:>4} {t[g]:9.3e} {t[min(t)] / t[g]:10.2f}")
    if plot:
        plt = pyplot()
        fig, ax = plt.subplots(figsize=(6, 4.5))
        Ns = sorted({k[1] for k in runs})
        for (solver, N, mach), t in sorted(runs.items()):
            g = sorted(t)
            ax.plot(
                g,
                [t[min(t)] / t[x] for x in g],
                "--o" if solver == "projection" else "-s",
                mfc="none",
                color=f"C{Ns.index(N)}",
                alpha=1.0 if solver == "explicit" or mach == min(k[2] for k in runs if k[0] == solver) else 0.5,
                label=f"{solver}, ${N}^3$/GPU" + (f", M = {mach:g}" if solver == "projection" else ""),
            )
        ax.axhline(1, color="k", lw=0.8)
        ax.set_xlabel("GPUs")
        ax.set_ylabel("weak-scaling efficiency (s/step on 1 GPU / on G)")
        ax.set_xticks(sorted({g for t in runs.values() for g in t}))
        ax.set_ylim(0, 1.1)
        ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(plot, dpi=150)


def ke_history(d):
    with open(os.path.join(d, "simulation.inp")) as f:
        inp = f.read()
    dt = float(re.search(r"^\s*dt\s*=\s*(\S+)", inp, re.M).group(1))
    nx, ny, nz = (int(re.search(rf"^\s*{v}\s*=\s*(\S+)", inp, re.M).group(1)) + 1 for v in "mnp")
    dv = (2 * np.pi / nx) * (2 * np.pi / ny) * (2 * np.pi / nz)
    steps = sorted(int(f.split("_")[-1][:-4]) for f in glob.glob(os.path.join(d, "restart_data", "lustre_[0-9]*.dat")))
    kx, ky, kz = (1j * np.fft.fftfreq(n, 1.0 / n) for n in (nx, ny, nz))  # 2 pi periodic
    d_ = [
        lambda f: np.fft.ifft(kx[None, None, :] * np.fft.fft(f, axis=2), axis=2).real,
        lambda f: np.fft.ifft(ky[None, :, None] * np.fft.fft(f, axis=1), axis=1).real,
        lambda f: np.fft.ifft(kz[:, None, None] * np.fft.fft(f, axis=0), axis=0).real,
    ]
    t, ek, zeta = [], [], []
    for s in steps:
        # alpha_rho, mom(1:3), E, alpha; x fastest
        q = np.fromfile(os.path.join(d, "restart_data", f"lustre_{s}.dat"), dtype=np.float64).reshape(-1, nz, ny, nx)
        u = q[1:4] / q[0]
        w2 = sum((d_[b](u[c]) - d_[c](u[b])) ** 2 for b, c in ((1, 2), (2, 0), (0, 1)))  # |curl u|^2, spectral
        t.append(s * dt * U0)
        ek.append(0.5 * (q[0] * (u**2).sum(0)).sum() * dv / (2 * np.pi) ** 3)
        zeta.append(0.5 * w2.mean())
    t, ek, zeta = np.array(t), np.array(ek), np.array(zeta)
    # Enstrophy in U0^2/L^2
    return t, ek / ek[0], -np.gradient(ek, t) / U0**2, zeta / U0**2


def ke(dirs, plot):
    runs = {os.path.basename(os.path.normpath(d)): ke_history(d) for d in dirs}
    (ref_name, (tr, er, dr, _)), *rest = runs.items()
    print(f"reference: {ref_name}, E_k/E_k(0) = {er[-1]:.5f} and peak dissipation {dr.max():.4e} at t/tC = {tr[dr.argmax()]:.2f}")
    print(f"{'run':<34} {'max |dE_k|/E_k(0)':>18} {'peak eps':>10} {'at t/tC':>8}")
    for name, (t, e, dis, _) in rest:
        n = min(len(t), len(tr))
        err = np.abs(np.interp(t[:n], tr, er) - e[:n]).max()
        print(f"{name:<34} {err:18.3e} {dis.max():10.4e} {t[dis.argmax()]:8.2f}")
    if plot:
        plt = pyplot()
        # run_sweep.sh names runs <solver>_N<N>_M<mach>_g<ngpu>[_<variant>]: color by Mach, solid explicit (dash-dot
        # for a variant), dashed projection with staggered markers so overlapping curves stay visible
        tags = {name: re.match(r"(explicit|projection)_N\d+_M([\d.e-]+)_g\d+(?:_(\w+))?$", name) for name in runs}
        machs = sorted({float(m[2]) for m in tags.values() if m}, reverse=True)
        fig, ((ax_e, ax_d), (ax_r, ax_z)) = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
        n_proj = 0
        for name, (t, e, dis, z) in runs.items():
            m = tags[name]
            style = {"label": name}
            if m:
                style = {"color": f"C{machs.index(float(m[2]))}", "ls": "--", "label": f"{m[1]}, M = {m[2]}" + (f" ({m[3]})" if m[3] else "")}
                if m[1] == "explicit":
                    style["ls"] = "-." if m[3] else "-"
                else:
                    style.update(marker="osd^v"[n_proj % 5], ms=5, markevery=(n_proj + 1, 3), mfc="none")
                    n_proj += 1
            ax_e.plot(t, e, **style)
            ax_d.plot(t, dis, **style)
            ax_z.plot(t, z, **style)
            if name != ref_name:
                ratio = np.abs(e / np.interp(t, tr, er) - 1)
                ax_r.semilogy(t[1:], ratio[1:], **style)  # both start at 1
        ref = f"{tags[ref_name][1]}, M = {tags[ref_name][2]}" if tags[ref_name] else ref_name
        ax_e.set_title("(a) kinetic energy", loc="left")
        ax_r.set_title(f"(b) kinetic energy ratio to the reference ({ref})", loc="left", fontsize=10)
        ax_d.set_title("(c) kinetic energy dissipation rate", loc="left")
        ax_z.set_title("(d) enstrophy", loc="left")
        ax_e.set_ylabel("$E_k/E_k(0)$")
        ax_r.set_ylabel(r"$|E_k/E_k^{ref} - 1|$")
        ax_d.set_ylabel(r"$-dE_k/dt$ $[U_0^3/L]$")
        ax_z.set_ylabel(r"$\zeta = \frac{1}{2}\langle|\omega|^2\rangle$ $[U_0^2/L^2]$")
        for a in (ax_r, ax_z):
            a.set_xlabel("$t/t_C$")
        ax_e.legend(fontsize=7, loc="lower left")
        fig.tight_layout()
        fig.savefig(plot, dpi=150)


if args.mode == "timing":
    timing(args.csv, args.plot)
elif args.mode == "weak":
    weak(args.csv, args.plot)
else:
    ke(args.dirs, args.plot)
