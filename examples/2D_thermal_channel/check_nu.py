#!/usr/bin/env python3
"""
Nusselt number of the thermally developed channel (2D_thermal_channel) against the exact eigenvalue solution.

    check_nu.py HOT_DIR [COLD_DIR] [--tol TOL]

Each directory is a finished run (it holds pre_process.inp, simulation.inp, mean_temp.dat and the last save in p_all/).
Per run, from the mean_temp.dat rows of the last --win fraction (T_m = mixing-cup temperature, qdot = source, V_f = fluid volume):

    Nu_src  = q_w D_h / (k (T_s - T_m)),  q_w = -qdot V_f / (2 L_x)   the source total is the wall heat at steady state
    Nu_grad = the same with q_w = -k dT/dn at both walls from the last save (second-order one-sided, independent of the source)
    hold    = max |T_m / T_target - 1| over the window

k is the conductivity at (T_s + T_m)/2. With two runs (wall temperatures above and below the target) it also prints the pair Nu
from half-differences, which cancels viscous heating and the first-order property asymmetry; exit 1 if the pair Nu_src (or the
single-run Nu_src) is farther than --tol from the exact value. D_h = 2H. The exact value comes from
psi'' + (6 Pe_H eta (1 - eta) b + b^2) psi = 0, psi(0) = psi(1) = 0, b = beta H (axial conduction included; 7.5407 as Pe -> inf).
"""

import argparse
import glob
import os
import re
import sys

import cantera as ct
import numpy as np

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("runs", nargs="+", help="One run directory, or the hot and the cold run directory.")
parser.add_argument("--win", type=float, default=0.1, help="Averaging window, fraction of the mean_temp.dat rows (default 0.1).")
parser.add_argument("--tol", type=float, default=0.01, help="Exit 1 if |Nu_src / Nu_exact - 1| exceeds this (default 1 %%).")
parser.add_argument("--mech", default="h2o2.yaml", help="Gas mechanism of the run (cantera_file in case.py).")
parser.add_argument("--dtype", default="float64", help="Storage type of p_all (float32 for single-precision builds).")
args = parser.parse_args()


def read_inp(path):
    out = {}
    for line in open(path):
        match = re.match(r"\s*([^=\s]+)\s*=\s*(\S+)", line)
        if match:
            out[match.group(1)] = match.group(2).strip("'\"")
    return out


def nu_exact(pe_h, n=400):
    """Exact Nu and beta*H of plane Poiseuille flow between isothermal plates at Peclet number Pe_H = rho cp U_b H / k"""
    eta = np.linspace(0.0, 1.0, n + 1)
    h = eta[1] - eta[0]
    u = 6.0 * eta[1:-1] * (1.0 - eta[1:-1])
    m = n - 1
    d2 = (np.diag(-2.0 * np.ones(m)) + np.diag(np.ones(m - 1), 1) + np.diag(np.ones(m - 1), -1)) / h**2
    # quadratic eigenproblem -D2 psi = (b Pe u + b^2) psi, linearized with phi = b psi
    w, v = np.linalg.eig(np.block([[np.zeros((m, m)), np.eye(m)], [-d2, -pe_h * np.diag(u)]]))
    ok = (np.abs(w.imag) < 1e-9 * np.abs(w)) & (w.real > 0)
    i = np.argmin(np.where(ok, w.real, np.inf))
    psi = np.concatenate([[0.0], v[:m, i].real, [0.0]])
    psi /= psi[n // 2]
    uu = 6.0 * eta * (1.0 - eta)
    psi_m = np.trapezoid(uu * psi, eta) / np.trapezoid(uu, eta)
    return 2.0 * (-3.0 * psi[0] + 4.0 * psi[1] - psi[2]) / (2.0 * h) / psi_m, w[i].real


def analyze(d):
    inp = {**read_inp(os.path.join(d, "pre_process.inp")), **read_inp(os.path.join(d, "simulation.inp"))}
    nx, ny = int(inp["m"]) + 1, int(inp["n"]) + 1
    height, length, t_wall, t_target = (float(inp[k]) for k in ("y_domain%end", "x_domain%end", "bc_y%Twall_in", "mean_temp"))
    gas = ct.Solution(args.mech)
    mass_fr = [float(inp[f"patch_icpp(1)%Y({i + 1})"]) for i in range(gas.n_species)]

    def conductivity(temp):
        gas.TPY = temp, float(inp["patch_icpp(1)%pres"]), mass_fr
        return gas.thermal_conductivity

    rows = np.loadtxt(os.path.join(d, "mean_temp.dat"))
    w = max(1, int(args.win * len(rows)))
    t_m, qdot, vol = (rows[-w:, c].mean() for c in (3, 4, 5))
    q_src = -qdot * vol / (2.0 * length)
    k_f = conductivity(0.5 * (t_wall + t_m))

    # wall gradient from the last save: T from the conservative state through the mechanism
    last = max(glob.glob(os.path.join(d, "p_all", "p0", "*")), key=lambda s: int(os.path.basename(s)))
    q = [np.fromfile(os.path.join(last, f"q_cons_vf{i}.dat"), dtype=args.dtype, offset=4, count=nx * ny).reshape(ny, nx) for i in range(1, 6 + gas.n_species)]
    rho, e_tot = q[0], q[3]
    kin = 0.5 * (q[1] ** 2 + q[2] ** 2) / rho**2
    temp = np.empty_like(rho)
    for j in range(ny):
        for i in range(nx):
            gas.UVY = e_tot[j, i] / rho[j, i] - kin[j, i], 1.0 / rho[j, i], [q[5 + s][j, i] / rho[j, i] for s in range(gas.n_species)]
            temp[j, i] = gas.T
    psi, dy = temp - t_wall, height / ny
    grad = ((9.0 * psi[0] - psi[1]) + (9.0 * psi[-1] - psi[-2])).mean() / (3.0 * dy) / 2.0  # mean dT/dn, n the distance from each wall
    q_grad = -conductivity(t_wall) * grad

    # Pe_H from the held mass flux and the target state
    gas.TPY = t_target, float(inp["patch_icpp(1)%pres"]), mass_fr
    pe_h = float(inp["mean_mflux(1)"]) * gas.cp_mass * height / gas.thermal_conductivity
    exact, beta_h = nu_exact(pe_h)
    nu = lambda qw, tm: qw * 2.0 * height / (k_f * (t_wall - tm))  # noqa: E731
    return dict(
        name=os.path.basename(os.path.abspath(d)),
        t_wall=t_wall, t_m=t_m, q_src=q_src, q_grad=q_grad, k_f=k_f, height=height, pe_h=pe_h, exact=exact, beta_h_exact=beta_h,
        nu_src=nu(q_src, t_m), nu_grad=nu(q_grad, t_m), beta_h=rows[-w:, 6].mean() * height,
        hold=np.abs(rows[-w:, 3] / t_target - 1.0).max(), steps=len(rows),
    )  # fmt: skip


res = [analyze(d) for d in args.runs]
for r in res:
    print(
        f"{r['name']}: T_s {r['t_wall']:.0f} K, T_m {r['t_m']:.3f} K (hold {r['hold']:.1e}), beta H {r['beta_h']:.5f} (exact {r['beta_h_exact']:.5f}), "
        f"Pe_H {r['pe_h']:.1f}, Nu_src {r['nu_src']:.4f} ({r['nu_src'] / r['exact'] - 1:+.2%}), Nu_grad {r['nu_grad']:.4f} ({r['nu_grad'] / r['exact'] - 1:+.2%}), "
        f"exact {r['exact']:.4f}, {r['steps']} steps"
    )
err = abs(res[0]["nu_src"] / res[0]["exact"] - 1.0)
if len(res) == 2:
    h, c = res
    d_t = 0.5 * ((h["t_wall"] - h["t_m"]) - (c["t_wall"] - c["t_m"]))
    k_f = 0.5 * (h["k_f"] + c["k_f"])
    pair = {lab: 0.5 * (h[f"q_{lab}"] - c[f"q_{lab}"]) * 2.0 * h["height"] / (k_f * d_t) for lab in ("src", "grad")}
    for lab, val in pair.items():
        print(f"pair half-difference Nu_{lab} {val:.4f} ({val / h['exact'] - 1:+.2%})")
    err = abs(pair["src"] / h["exact"] - 1.0)
print(f"{'PASS' if err <= args.tol else 'FAIL'}: Nu_src error {err:.2e} (tolerance {args.tol:g})")
sys.exit(0 if err <= args.tol else 1)
