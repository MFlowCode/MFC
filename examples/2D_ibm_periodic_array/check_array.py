#!/usr/bin/env python3
"""
Drag, force statistics and Nusselt number of a periodic fixed array (2D_ibm_periodic_array, 3D_ibm_periodic_array).

    check_array.py CASE_DIR [--avg 0.4] [--nu]

Reads pre_process.inp, simulation.inp, D/ib_forces.dat, mean_dpdx.dat (and D/ib_heat.dat, mean_temp.dat with --nu); no
post_process output is needed. Averages are over the last --avg of the run. It prints

  Q / Q*        fluid-averaged mass flux held by the controller over its target
  balance       sum of the body forces over (-dpdx V_f): the exact global momentum balance of a periodic box, which holds at
                steady state (V_f = fluid cells of the staircase body mask, the controller's definition)
  F (3D)        total force per body over F_St = 3 pi mu D (1 - phi) W, with the total force the surface force over (1 - phi)
                (Tenneti et al. 2011, Eqs. 33 and 38), and Tenneti Eq. 43 at the run's phi and Re_m; 2D has no reference
  Cd (2D)       total force per unit length over (1/2) rho W_s^2 D, with W_s = (1 - phi) W the superficial velocity
  spread        standard deviation over the mean of the time-averaged streamwise force of the bodies
  Nu            --nu: heat the bodies give the gas over N A_p k_f (T_s - T_m), A_p = pi D^2 (3D) or pi D (2D), k_f at the film
                temperature; Gunn (1978) and Sun et al. (2015) correlations in 3D, T_m = mean_temp
"""

import argparse
import math
import os
import re

import numpy as np


def tenneti(phi, re):
    """Tenneti, Garg and Subramaniam (2011), Eq. 43: total force per body over F_St (0.1 <= phi <= 0.5, Re_m <= 300)."""
    f_isol = 1.0 + 0.15 * re**0.687
    f_phi = 5.81 * phi / (1 - phi) ** 3 + 0.48 * phi ** (1 / 3) / (1 - phi) ** 4
    f_phi_re = phi**3 * re * (0.95 + 0.61 * phi**3 / (1 - phi) ** 2)
    return f_isol / (1 - phi) ** 3 + f_phi + f_phi_re


def gunn(phi, re, pr):
    """Gunn (1978) particle Nusselt number of a fixed bed, with epsilon = 1 - phi and Re = Re_m."""
    e = 1 - phi
    return (7 - 10 * e + 5 * e**2) * (1 + 0.7 * re**0.2 * pr ** (1 / 3)) + (1.33 - 2.4 * e + 1.2 * e**2) * re**0.7 * pr ** (1 / 3)


def sun2015(phi, re, pr):
    """Sun, Tenneti and Subramaniam (2015) mixing-cup Nusselt number per cross-section perimeter (4/pi times the surface-area Nu)."""
    e = 1 - phi
    return (-0.46 + 1.77 * e + 0.69 * e**2) / e**3 + (1.37 - 2.4 * e + 1.2 * e**2) * re**0.7 * pr ** (1 / 3)


def read_inp(run):
    out = {}
    for name in ("pre_process.inp", "simulation.inp"):
        for line in open(os.path.join(run, name)):
            m = re.match(r"\s*([^=\s]+)\s*=\s*(\S+)", line)
            if m:
                out[m.group(1)] = m.group(2).strip("'\"")
    return out


def setup(run):
    inp = read_inp(run)
    f = lambda k: float(inp[k].replace("d", "e").replace("D", "e"))  # noqa: E731
    dim = 3 if int(inp["p"]) > 0 else 2
    n = int(inp["num_ibs"])
    R, L, cells = f("patch_ib(1)%radius"), f("x_domain%end"), int(inp["m"]) + 1
    centers = np.array([[f(f"patch_ib({i + 1})%{c}_centroid") for c in "xyz"[:dim]] for i in range(n)])
    # Fluid volume of the staircase mask: cell centres outside every body, with periodic images
    x = (np.arange(cells) + 0.5) * L / cells
    grid = np.meshgrid(*([x] * dim), indexing="ij")
    inside = np.zeros(grid[0].shape, bool)
    for c in centers:
        r2 = sum((np.minimum(np.abs(g - ci), L - np.abs(g - ci))) ** 2 for g, ci in zip(grid, c))
        inside |= r2 < R * R
    Vf = float((~inside).sum()) * (L / cells) ** dim
    return dict(
        dim=dim,
        cells=cells,
        n=n,
        D=2 * R,
        L=L,
        rho=f("patch_icpp(1)%alpha_rho(1)"),
        mu=1.0 / f("fluid_pp(1)%Re(1)"),
        Q=f("mean_mflux(1)"),
        Vf=Vf,
        phi=1 - Vf / L**dim,
        T=f("patch_ib(1)%Twall"),
        Tm=float(inp.get("mean_temp", "nan")),
        p=f("patch_icpp(1)%pres"),
    )


def table(path):
    """A numeric text file; the unfinished last line of a run that stopped early is dropped."""
    lines = [ln for ln in open(path) if ln.strip()]
    n = len(lines[0].split())
    return np.array([[float(x) for x in ln.split()] for ln in lines if len(ln.split()) == n]).reshape(-1, n)


def window(t, y, avg):
    m = t >= t[-1] - avg * (t[-1] - t[0])
    return float(np.average(y[m], weights=np.gradient(t)[m])) if m.sum() > 1 else float(y[-1])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("--avg", type=float, default=0.4)
    ap.add_argument("--nu", action="store_true")
    a = ap.parse_args()
    s = setup(a.run)
    d = table(os.path.join(a.run, "D", "ib_forces.dat"))
    ids = d[:, 0].astype(int)
    times = np.unique(d[:, 1])
    fx = np.array([d[(ids == i + 1), 2] for i in range(s["n"])])[:, : len(times)]  # bodies x records
    mp = table(os.path.join(a.run, "mean_dpdx.dat"))
    t_end = times[-1]
    q, dpdx = window(mp[:, 1], mp[:, 2], a.avg), window(mp[:, 1], mp[:, 5], a.avg)  # columns: step, t, Q(1:3), dpdx(1:3)
    wm = window(times, fx.sum(axis=0), a.avg)
    dim, n, D, phi = s["dim"], s["n"], s["D"], s["phi"]
    W = s["Q"] / s["rho"]
    re_m = s["rho"] * (1 - phi) * W * D / s["mu"]
    f_tot = wm / n / (1 - phi)
    per_body = np.array([window(times, fx[i], a.avg) for i in range(n)])
    print(f"{dim}D  N {n}  phi {phi:.4f}  Re_m {re_m:.2f}  cells/D {s['cells'] * D / s['L']:.1f}  t_end {t_end:.4g} s  (averages over the last {a.avg:.0%})")
    print(f"  Q / Q*       {q / s['Q']:.6f}")
    print(f"  balance      sum(F_x) / (-dpdx V_f) = {wm / (-dpdx * s['Vf']):.4f}   (dpdx {dpdx:.5g} Pa/m)")
    if dim == 3:
        f_st = 3 * math.pi * s["mu"] * D * (1 - phi) * W
        ref = tenneti(phi, re_m)
        print(f"  F/F_St       {f_tot / f_st:.3f}   Tenneti Eq. 43 {ref:.3f}   ({100 * (f_tot / f_st / ref - 1):+.1f} %)")
    else:
        print(f"  Cd           {f_tot / (0.5 * s['rho'] * ((1 - phi) * W) ** 2 * D):.3f}   (no 2D reference)")
    print(f"  spread       std/mean of per-body <F_x> = {per_body.std() / per_body.mean():.3f}   (min {per_body.min() / per_body.mean():.2f}, max {per_body.max() / per_body.mean():.2f})")
    if a.nu:
        import cantera as ct

        h = table(os.path.join(a.run, "D", "ib_heat.dat"))
        ts = np.unique(h[:, 1])
        qs = np.array([h[h[:, 0].astype(int) == i + 1, 4] for i in range(n)])[:, : len(ts)].sum(axis=0)
        q_gas = -window(ts, qs, a.avg)  # heat the bodies give the gas
        g = ct.Solution("h2o2.yaml")
        g.TPX = 0.5 * (s["T"] + s["Tm"]), s["p"], "O2:0.21,N2:0.79"
        k_f, pr = g.thermal_conductivity, g.viscosity * g.cp_mass / g.thermal_conductivity
        area = math.pi * D**2 if dim == 3 else math.pi * D
        nu = q_gas * D / (n * area * k_f * (s["T"] - s["Tm"]))
        print(f"  Nu           {nu:.3f}   (Pr_f {pr:.3f}, T_s {s['T']:.1f} K, T_m {s['Tm']:.1f} K)")
        if dim == 3:
            print(f"               Gunn {gunn(phi, re_m, pr):.3f} ({100 * (nu / gunn(phi, re_m, pr) - 1):+.1f} %)   Sun et al. (surface-area basis) {sun2015(phi, re_m, pr) * math.pi / 4:.3f}")


if __name__ == "__main__":
    main()
