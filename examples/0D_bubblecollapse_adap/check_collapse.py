#!/usr/bin/env python3
"""
Bubble collapse of 0D_bubblecollapse_adap against the Keller-Miksis ODE and Rayleigh's collapse time.

    check_collapse.py CASE_DIR [--case-args ARGS] [--tol TOL]

CASE_DIR is a finished run (pre_process, simulation, post_process). The domain is uniform and periodic, so every cell
holds the same bubble; R = nR001/n is read from the Silo saves (cell 0) at t = step * dt. The reference integrates the
Keller-Miksis equation (polytropic gas, surface tension, liquid viscosity, the liquid sound speed of the initial state)
in the case's units with a constant far-field pressure, the initial pressure of the patch. Prints the time and radius
of the first minimum for both (the radius at the same save times), the Rayleigh time 0.915 R0 sqrt(rho_l/(p_inf - p_v)), and max |R_MFC - R_KM|/R0 over the
saves before the first minimum; the same with the ODE driven by MFC's cell pressure, which the collapsing
void fraction lowers (two-way coupling). Exits 1 if that error exceeds --tol.
"""

import argparse
import json
import math
import subprocess
import sys

import numpy as np
from scipy.integrate import solve_ivp

from mfc.viz import assemble_silo, discover_timesteps


def km_rhs(case, p_of_t=None):
    """Keller-Miksis right-hand side as in m_bubbles.fpp; far-field pressure constant, or p_of_t(t)."""
    R0, gam = case["bub_pp%R0ref"], case["bub_pp%gam_g"]
    pv, We_inv, mu = case["bub_pp%pv"], case["bub_pp%ss"], case["bub_pp%mu_l"]
    Ca = case["bub_pp%p0ref"] - pv
    p_inf = case["patch_icpp(1)%pres"]
    rho = case["patch_icpp(1)%alpha_rho(1)"] / (1.0 - case["patch_icpp(1)%alpha(1)"])
    n_tait = 1.0 + 1.0 / case["fluid_pp(1)%gamma"]
    B = case["fluid_pp(1)%pi_inf"] * (n_tait - 1.0) / n_tait
    c = math.sqrt(n_tait * (p_inf + B) / rho)
    pg0 = Ca + 2.0 * We_inv / R0

    def rhs(t, y):
        R, V = y
        p_far = p_inf if p_of_t is None else p_of_t(t)
        pg = pg0 * (R0 / R) ** (3.0 * gam)
        pbw = pg + pv - 2.0 * We_inv / R - 4.0 * mu * V / R
        dpbw = -3.0 * gam * pg * V / R + 2.0 * We_inv * V / R**2 + 4.0 * mu * (V / R) ** 2
        num = 1.5 * V**2 * (V / (3.0 * c) - 1.0) + (1.0 + V / c) * (pbw - p_far) / rho + dpbw * R / (rho * c)
        return [V, num / (R * (1.0 - V / c) + 4.0 * mu / (rho * c))]

    return rhs, rho, p_inf


def first_min(t, R):
    """Time (parabola through the three saves around it) and sampled radius of the first minimum, and its index."""
    i = int(np.argmax((R[1:-1] <= R[:-2]) & (R[1:-1] <= R[2:]))) + 1
    y0, y1, y2 = R[i - 1 : i + 2]
    den = y0 - 2.0 * y1 + y2
    off = 0.5 * (y0 - y2) / den if den != 0 else 0.0
    return t[i] + off * (t[i + 1] - t[i]), y1, i


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case_dir")
    ap.add_argument("--case-args", default="", help="Arguments given to case.py for the run.")
    ap.add_argument("--tol", type=float, default=0.0, help="Fail if max |R_MFC - R_KM|/R0 before the first minimum exceeds this (0: no check).")
    args = ap.parse_args()

    d = args.case_dir
    case = json.loads(subprocess.run([sys.executable, "case.py", *args.case_args.split()], cwd=d, capture_output=True, text=True, check=True).stdout)
    steps = discover_timesteps(d, "silo")
    t = np.array(steps) * case["dt"]
    R, P = [], []
    for s in steps:
        v = {k: float(np.asarray(a).ravel()[0]) for k, a in assemble_silo(d, s).variables.items()}
        R.append(v["nR001"] / v["n"])
        P.append(v["pres"])
    R, P = np.array(R), np.array(P)

    R0 = case["bub_pp%R0ref"]
    rhs, rho, p_inf = km_rhs(case)
    sol = solve_ivp(rhs, (0.0, t[-1]), [case["patch_icpp(1)%r0"], case["patch_icpp(1)%v0"]], method="LSODA", rtol=1e-11, atol=1e-13, dense_output=True)
    R_km = sol.sol(t)[0]
    # the same ODE driven by MFC's own cell pressure, which falls as the void fraction collapses (two-way coupling)
    rhs_p, _, _ = km_rhs(case, lambda s: np.interp(s, t, P))
    R_kmp = solve_ivp(rhs_p, (0.0, t[-1]), sol.y[:, 0], method="LSODA", rtol=1e-11, atol=1e-13, t_eval=t, max_step=t[1] - t[0]).y[0]
    t_c = 0.915 * R0 * math.sqrt(rho / (p_inf - case["bub_pp%pv"]))

    tm, Rm, i = first_min(t, R)
    tk, Rk, _ = first_min(t, R_km)
    err = float(np.max(np.abs(R[: i + 1] - R_km[: i + 1]))) / R0
    tf = np.linspace(t[i - 1], t[i + 1], 4001)
    R_min_km = float(sol.sol(tf)[0].min())
    print(f"saves {len(t)}, dt {case['dt']:.3e}, Rayleigh time t_c = {t_c:.5e}")
    print(f"first minimum: MFC t = {tm:.5e} ({tm / t_c:.4f} t_c), R = {Rm:.4e} R0;  Keller-Miksis t = {tk:.5e} ({tk / t_c:.4f} t_c), R = {Rk:.4e} R0")
    print(f"collapse time MFC vs KM {tm / tk - 1:+.3e}; sampled minimum radius MFC vs KM {Rm / Rk - 1:+.3e} (KM between saves: {R_min_km:.4e} R0)")
    print(f"max |R_MFC - R_KM| / R0 up to the first minimum: {err:.3e}; over the whole run: {np.max(np.abs(R - R_km)) / R0:.3e}")
    print(
        f"cell pressure {P.min() / P[0] - 1:+.2e} at most from its initial value; KM driven by it: max |R_MFC - R_KM(p_MFC)| / R0 {np.max(np.abs(R[: i + 1] - R_kmp[: i + 1])) / R0:.3e}, first minimum at {first_min(t, R_kmp)[0]:.5e}"
    )
    ok = args.tol <= 0 or err <= args.tol
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
