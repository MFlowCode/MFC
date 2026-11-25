#!/usr/bin/env python3
"""
Shock radius and peak density of the Sedov-Taylor blast against the self-similar solution.

    python3 check_sedov.py <run dir> [<run dir> ...]

Per save: the blast energy from the data (the total energy above the ambient; the blast is at the origin), the shock radius R (outermost
radius where the shell-averaged density, on the finest cells, is halfway between rho0 and its peak), the exact radius
xi0 (E t^2 / rho0)^(1/5), the signed error, and the peak density over rho0 (strong-shock limit (gamma+1)/(gamma-1) = 6).
xi0 comes from the energy integral of the similarity solution, integrated here.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "2D_amr_isentropic_vortex"))
from amr_composite import _num, leaves, read_inp  # noqa: E402

GAM, RHO0, P0 = 1.4, 1.0, 1.0e-3


def xi0(gam=GAM, n=20000):
    """Sedov constant for a spherical blast: u = R' U, rho = rho0 G, p = rho0 R'^2 P of xi = r/R, integrated inward from the
    strong-shock state at xi = 1, then E = 4 pi rho0 (4/25) xi0^5 E I with I = int (G U^2/2 + P/(gamma-1)) xi^2 dxi."""

    def rhs(x, s):
        U, G, P = s
        w = U - x
        # continuity, momentum, entropy: a 3x3 linear system in (U', G', P')
        A = np.array([[G, w, 0.0], [w, 0.0, 1.0 / G], [0.0, -gam * P / G * w, w]])
        b = np.array([-2.0 * G * U / x, 1.5 * U, 3.0 * P])
        return np.linalg.solve(A, b)

    s = np.array([2.0 / (gam + 1.0), (gam + 1.0) / (gam - 1.0), 2.0 / (gam + 1.0)])
    xs = np.linspace(1.0, 1.0e-4, n)
    integ = 0.0
    f = lambda x, s: (0.5 * s[1] * s[0] ** 2 + s[2] / (gam - 1.0)) * x**2  # noqa: E731
    for x, x2 in zip(xs[:-1], xs[1:]):
        h = x2 - x
        k1 = rhs(x, s)
        k2 = rhs(x + h / 2, s + h / 2 * k1)
        k3 = rhs(x + h / 2, s + h / 2 * k2)
        k4 = rhs(x2, s + h * k3)
        s2 = s + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        integ -= h / 2 * (f(x, s) + f(x2, s2))
        s = s2
    return (25.0 / (16.0 * np.pi * integ)) ** 0.2


def save_times(run_dir, t_save):
    """{save index: time} from the run_time.inf table (step, dt, time at the step's start, ICFL): a save is written after the
    first step whose end time passes a multiple of t_save. The table has 4 significant digits."""
    times, last = {}, 0
    with open(os.path.join(run_dir, "run_time.inf")) as f:
        for line in f:
            cols = line.split()
            if len(cols) < 3 or not cols[0].isdigit():
                continue
            t = float(cols[2]) + float(cols[1])
            k = int(t / t_save + 1e-9)
            if k > last:
                times[k], last = t, k
    return times


def measure(run_dir, k, t, x0):
    parts, _ = leaves(run_dir, k, var=0)
    eparts, _ = leaves(run_dir, k, var=4)
    r, rho, vol, e = [], [], [], 0.0
    for (c, dx, q, _), (_, _, qe, _) in zip(parts, eparts):
        rr = np.sqrt(sum(a**2 for a in c))
        r.append(rr)
        rho.append(q)
        vol.append(np.full(q.size, np.prod(dx)))
        e += float(((qe - P0 / (GAM - 1.0)) * np.prod(dx)).sum())
    r, rho, vol = map(np.concatenate, (r, rho, vol))
    dr = min(p[1].min() for p in parts if p[2].size)
    edges = np.arange(0.0, r.max() + dr, dr)
    idx = np.digitize(r, edges) - 1
    prof = np.bincount(idx, rho * vol, edges.size) / np.maximum(np.bincount(idx, vol, edges.size), 1e-300)
    ok = np.bincount(idx, vol, edges.size) > 0
    centers = edges + dr / 2
    peak = int(np.argmax(np.where(ok, prof, 0.0)))
    half = RHO0 + 0.5 * (prof[peak] - RHO0)
    j = peak
    while j + 1 < edges.size and (not ok[j + 1] or prof[j + 1] > half):
        j += 1
    # linear interpolation of the half-jump crossing between bins j and j+1
    R = centers[j] + dr * (prof[j] - half) / max(prof[j] - prof[j + 1], 1e-300)
    return e, R, x0 * (e * t**2 / RHO0) ** 0.2, prof[peak] / RHO0


def main(dirs):
    x0 = xi0()
    print(f"xi0 (gamma = {GAM}) = {x0:.4f}")
    print(f"{'run':40s} {'t':>7s} {'E':>8s} {'R':>7s} {'R exact':>8s} {'error':>8s} {'rho max':>8s}")
    for d in dirs:
        inp = read_inp(d)
        label = os.path.basename(os.path.normpath(d))
        times = save_times(d, _num(inp["t_save"]))
        # the last save is at t_stop: the run shortens its last step to land there
        times.setdefault(round(_num(inp["t_stop"]) / _num(inp["t_save"])), _num(inp["t_stop"]))
        for k, t in sorted(times.items()):
            if not os.path.isfile(os.path.join(d, "restart_data", f"lustre_{k}.dat")):
                continue
            E, R, Rx, pk = measure(d, k, t, x0)
            print(f"{label:40s} {t:7.4f} {E:8.4f} {R:7.4f} {Rx:8.4f} {100 * (R - Rx) / Rx:+7.2f}% {pk:8.3f}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1:])
