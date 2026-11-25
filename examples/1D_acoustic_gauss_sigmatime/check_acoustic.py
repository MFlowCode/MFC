#!/usr/bin/env python3
"""
1D planar acoustic sources (acoustic(1)%support = 1) against linear acoustics.

    check_acoustic.py RUN_DIR [RUN_DIR ...] [--case-args ARGS] [--step N] [--tol TOL]

Each RUN_DIR is a finished run (pre_process, simulation, post_process) of one of the 1D acoustic examples
(1D_acoustic_gauss_sigmatime, _gauss_sigmadist, _sine_frequency, _sine_wavelength, _dipole); --case-args are
the arguments its case.py was run with (one string, applied to every RUN_DIR; give one per run with ';').

The source adds S = mag g(t) f(x) to the momentum, S/c to the mass and S c / (gamma - 1) to the energy, so
p_t = c S and the Riemann invariant p + rho c u gains 2 c S while p - rho c u gains nothing: one pulse goes in
the direction sign(dir). With dipole = T the mass and energy terms are doubled and the momentum term dropped,
so both invariants gain 2 c S: two pulses of the same sign. Away from the support, each pulse is

    p'(x, t) = mag sum_j f_j dx_j g(t - delay - |x - x_j| / c)

with f_j the discrete Gaussian weights of the support (standard deviation 1.25 dx_j). The script evaluates
this at the saved step (default: the last one with the whole signal inside the domain), and prints the
relative L2 and max errors of p' outside +-5 support widths, the ratio of the max |p'|, and the amplitude A and
time shift tau of a least-squares fit p'_MFC = A p'_lin(t + tau) (tau < 0: MFC lags). With several runs (a resolution sweep) it also prints the observed order of the L2 error.
Exit 1 if the relative L2 error of the finest run exceeds --tol.
"""

import argparse
import json
import math
import os
import subprocess
import sys

import numpy as np

from mfc.viz import assemble, assemble_silo, discover_format, discover_timesteps


def load(d, step):
    fmt = discover_format(d)
    a = assemble_silo(d, step) if fmt == "silo" else assemble(d, step, fmt)
    return a.x_cc, a.variables["pres"]


def pulse(case, t, c):
    """Temporal waveform g(t) of s_source_temporal (pulse 1 = sine, 2 = Gaussian)."""
    td = t - case.get("acoustic(1)%delay", 0.0)
    if case["acoustic(1)%pulse"] == 1:
        f = case.get("acoustic(1)%frequency", None) or c / case["acoustic(1)%wavelength"]
        return np.where((td >= 0) & (td * f <= case["acoustic(1)%npulse"]), np.sin(2 * math.pi * f * td), 0.0)
    sig = case.get("acoustic(1)%gauss_sigma_time", None) or case["acoustic(1)%gauss_sigma_dist"] / c
    return np.exp(-0.5 * (td / sig) ** 2)


def predict(case, x, xcb, t, c):
    """p' of linear acoustics from the discrete source weights; also the support half-width used to mask the source."""
    x0, mag = case["acoustic(1)%loc(1)"], case["acoustic(1)%mag"]
    dx = np.diff(xcb)
    sig = 1.25 * dx
    w = np.exp(-0.5 * ((x - x0) / sig) ** 2) / (math.sqrt(2 * math.pi) * sig) * dx
    keep = w > 1e-14
    xj, wj = x[keep], w[keep]
    dirs = (1.0, -1.0) if case.get("acoustic(1)%dipole", "F") == "T" else (math.copysign(1.0, case.get("acoustic(1)%dir", 1.0)),)
    p = np.zeros_like(x)
    for s in dirs:
        side = s * (x - x0) > 0
        tau = t - s * (x[side, None] - xj[None, :]) / c
        p[side] += mag * (pulse(case, tau, c) * wj[None, :]).sum(axis=1)
    return p, 5 * sig[np.argmin(np.abs(x - x0))]


def check(d, case_args, step_arg):
    case = json.loads(subprocess.run([sys.executable, "case.py", *case_args.split()], cwd=d, capture_output=True, text=True, check=True).stdout)
    G, piinf, rho0, p0 = case["fluid_pp(1)%gamma"], case.get("fluid_pp(1)%pi_inf", 0.0), case["patch_icpp(1)%alpha_rho(1)"], case["patch_icpp(1)%pres"]
    gam = 1 + 1 / G
    c = math.sqrt(gam * (p0 + piinf / (gam * G)) / rho0)
    steps = discover_timesteps(d, discover_format(d))
    x, _ = load(d, steps[0])
    xcb = np.concatenate([[case["x_domain%beg"]], 0.5 * (x[1:] + x[:-1]), [0.0]])
    xcb[-1] = 2 * x[-1] - xcb[-2]
    if step_arg is None:  # last save with the predicted signal inside the domain
        for s in reversed(steps):
            pp, _ = predict(case, x, xcb, s * case["dt"], c)
            if max(abs(pp[0]), abs(pp[-1])) < 1e-6 * np.abs(pp).max():
                step_arg = s
                break
    t = step_arg * case["dt"]
    x, pres = load(d, step_arg)
    pe, mask_w = predict(case, x, xcb, t, c)
    pm = pres - p0
    far = np.abs(x - case["acoustic(1)%loc(1)"]) > mask_w
    err = pm[far] - pe[far]
    l2 = math.sqrt(np.sum(err**2) / np.sum(pe[far] ** 2))
    linf = np.abs(err).max() / np.abs(pe).max()
    am, ae = np.abs(pm[far]).max(), np.abs(pe[far]).max()
    # Least-squares fit pm = A pe(t + tau): tau from the linearization in tau, then A
    eps = 1e-3 * case["dt"]
    dpe = (predict(case, x, xcb, t + eps, c)[0] - predict(case, x, xcb, t - eps, c)[0])[far] / (2 * eps)
    tau = float(np.dot(err, dpe) / np.dot(dpe, dpe))
    pt = predict(case, x, xcb, t + tau, c)[0][far]
    amp = float(np.dot(pm[far], pt) / np.dot(pt, pt))
    dxs = np.diff(xcb)
    print(f"{d}: N = {x.size}, dx = {dxs.min():.3e}..{dxs.max():.3e} m, c = {c:.2f} m/s, dt = {case['dt']:.3e} s, step {step_arg}, t = {t:.4e} s")
    print(f"  max |p'|: MFC {am:.5f} Pa, linear {ae:.5f} Pa (mag {case['acoustic(1)%mag']}): ratio - 1 = {am / ae - 1:+.3e}")
    print(f"  fit p'_MFC = A p'_lin(t + tau): A - 1 = {amp - 1:+.3e}, tau = {tau * 1e9:+.4f} ns = {tau / case['dt']:+.3f} dt")
    l2fit = math.sqrt(np.sum((pm[far] - amp * pt) ** 2) / np.sum(pt**2))
    print(f"  rel. L2 error = {l2:.3e} ({l2fit:.3e} after the fit), max error / max p' = {linf:.3e}")
    return x.size, l2


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--case-args", default="", help="case.py arguments of the runs; ';' separates one entry per run")
    ap.add_argument("--step", type=int, default=None, help="Saved step to compare (default: last with the signal inside)")
    ap.add_argument("--tol", type=float, default=0.05, help="Max relative L2 error of the finest run")
    args = ap.parse_args()
    cargs = args.case_args.split(";")
    if len(cargs) == 1:
        cargs *= len(args.runs)
    res = [check(os.path.abspath(d), a, args.step) for d, a in zip(args.runs, cargs)]
    for (n1, e1), (n2, e2) in ((a, b) for a, b in zip(res, res[1:]) if a[0] != b[0]):
        print(f"observed order {n1} -> {n2}: {math.log(e1 / e2) / math.log(n2 / n1):.2f}")
    fine = min(res, key=lambda r: r[1])[1] if len(res) > 1 else res[0][1]
    if fine > args.tol:
        print(f"FAIL: rel. L2 error {fine:.3e} > {args.tol}")
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    main()
