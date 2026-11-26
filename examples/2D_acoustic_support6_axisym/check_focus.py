#!/usr/bin/env python3
"""
Focal pressure of the axisymmetric focused transducer (2D_acoustic_support6_axisym) against the on-axis solution of a
spherical cap.

    check_focus.py RUN_DIR [RUN_DIR ...] [--case-args ARGS] [--tol TOL]

Each RUN_DIR is a finished simulation run with --probe (D/probe1_prim.dat holds the pressure in the probe cell every
step, row k at t = k dt: it is written from the state at the start of step k); --case-args are its case.py arguments (';' separates one entry per run).

The source emits a converging spherical wave from a cap of radius of curvature F (acoustic(1)%foc_length) and
aperture 2a, with pressure p0(t) = mag g(t - delay) at the cap. On the axis, at distance z from the cap apex, the
exact linear solution for a cap with uniform normal velocity (O'Neil 1949) is

    p(z, t) = F / (F - z) [p0(t - z/c) - p0(t - R_e(z)/c)],   R_e = sqrt(a^2 + (z - h)^2),  h = F - sqrt(F^2 - a^2),

which at the focus becomes p(F, t) = (h / c) dp0/dt (t - F/c): two peaks of +-mag h / (c sigma sqrt(e)) for a
Gaussian g. The script evaluates it at the probe cell centre and prints the positive and negative peak values and
times (parabola through the extremum), their errors, the L2 error of the signal relative to its L2 norm, and the
least-squares time shift tau of p_MFC = p_exact(t + tau) with the L2 error left after it (tau < 0: MFC lags).
Exit 1 if the finest run's relative L2 error exceeds --tol.
"""

import argparse
import json
import math
import os
import subprocess
import sys

import numpy as np


def g(case, t, c):
    sig = case.get("acoustic(1)%gauss_sigma_time", None) or case["acoustic(1)%gauss_sigma_dist"] / c
    return case["acoustic(1)%mag"] * np.exp(-0.5 * ((t - case["acoustic(1)%delay"]) / sig) ** 2), sig


def extremum(t, y, sign):
    i = int(np.clip(np.argmax(sign * y), 1, y.size - 2))
    y0, y1, y2 = y[i - 1 : i + 2]
    den = y0 - 2 * y1 + y2
    off = 0.5 * (y0 - y2) / den if den != 0 else 0.0
    return t[i] + off * (t[1] - t[0]), y1 - 0.25 * (y0 - y2) * off


def check(d, case_args):
    case = json.loads(subprocess.run([sys.executable, "case.py", *case_args.split()], cwd=d, capture_output=True, text=True, check=True).stdout)
    G, piinf, rho0, p_amb = case["fluid_pp(1)%gamma"], case.get("fluid_pp(1)%pi_inf", 0.0), case["patch_icpp(1)%alpha_rho(1)"], case["patch_icpp(1)%pres"]
    gam = 1 + 1 / G
    c = math.sqrt(gam * (p_amb + piinf / (gam * G)) / rho0)
    F, a = case["acoustic(1)%foc_length"], 0.5 * case["acoustic(1)%aperture"]
    h = F - math.sqrt(F * F - a * a)
    dx = (case["x_domain%end"] - case["x_domain%beg"]) / (case["m"] + 1)
    xp = case["x_domain%beg"] + (math.floor((case["probe(1)%x"] - case["x_domain%beg"]) / dx) + 0.5) * dx
    z = xp - case["acoustic(1)%loc(1)"]
    data = np.loadtxt(os.path.join(d, "D", "probe1_prim.dat"))
    t = np.arange(data.shape[0]) * case["dt"]
    pm = data[:, 4] - p_amb
    if abs(F - z) > 1e-9 * F:
        re = math.sqrt(a * a + (z - h) ** 2)
        pe = F / (F - z) * (g(case, t - z / c, c)[0] - g(case, t - re / c, c)[0])
    else:
        p0, sig = g(case, t - F / c, c)
        pe = h / c * p0 * -(t - F / c - case["acoustic(1)%delay"]) / sig**2
    sig = g(case, t, c)[1]
    print(f"{d}: dx = {dx:.3e} m ({c * sig / dx:.1f} cells per c sigma), c = {c:.2f} m/s, h = {h * 1e3:.4f} mm, probe z - F = {(z - F) * 1e6:+.2f} um")
    print(f"  focal gain h / (c sigma sqrt(e)) = {h / (c * sig * math.sqrt(math.e)):.4f}")
    for name, s in (("positive", 1), ("negative", -1)):
        tm, am = extremum(t, pm, s)
        te, ae = extremum(t, pe, s)
        print(f"  {name} peak: MFC {am:+.4f} Pa at {tm * 1e9:.2f} ns, exact {ae:+.4f} Pa at {te * 1e9:.2f} ns: ratio - 1 = {am / ae - 1:+.3e}, time error {(tm - te) * 1e9:+.2f} ns")
    l2 = math.sqrt(np.sum((pm - pe) ** 2) / np.sum(pe**2))
    # Least-squares shift tau of pm = pe(t + tau), from the linearization in tau; then the L2 error left after it
    dpe = np.gradient(pe, t)
    tau = float(np.dot(pm - pe, dpe) / np.dot(dpe, dpe))
    l2s = math.sqrt(np.sum((pm - np.interp(t + tau, t, pe)) ** 2) / np.sum(pe**2))
    print(f"  rel. L2 error = {l2:.3e}; best shift tau = {tau * 1e9:+.3f} ns = {tau / case['dt']:+.2f} dt, L2 after it = {l2s:.3e}")
    return dx, l2


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--case-args", default="--probe", help="case.py arguments of the runs; ';' separates one entry per run")
    ap.add_argument("--tol", type=float, default=0.3, help="Max relative L2 error of the finest run")
    args = ap.parse_args()
    cargs = args.case_args.split(";")
    if len(cargs) == 1:
        cargs *= len(args.runs)
    res = sorted((check(os.path.abspath(d), ca) for d, ca in zip(args.runs, cargs)), reverse=True)
    for (h1, e1), (h2, e2) in zip(res, res[1:]):
        print(f"observed order dx {h1:.2e} -> {h2:.2e}: {math.log(e1 / e2) / math.log(h1 / h2):.2f}")
    if res[-1][1] > args.tol:
        print(f"FAIL: rel. L2 error {res[-1][1]:.3e} > {args.tol}")
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    main()
