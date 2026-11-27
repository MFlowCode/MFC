#!/usr/bin/env python3
"""
Gas-mass budget, conducted heat and force of the runs of 2D_ibm_moving_mass_budget.

    check_budget.py mass RUN [RUN ...]
    check_budget.py twin --moving RUN --fixed RUN [--avg 0.5]

mass: the gas mass M = sum(rho dV) over the cells whose centre is outside the body (at its prescribed position, with the periodic
      wrap) at every save, as the change from the first save in percent of M, against the distance the body has moved relative
      to the gas. The closed box holds a constant mass, so the change is created by the body. Prints the change per save, and the
      steady rate per diameter travelled after the first diameter (the impulsive start is excluded). The mass oscillates with the
      stair-step geometry, so read the rate over whole cell-crossing periods: the table uses all saves after the first diameter.
twin: a viscous run with --Twall. Conducted heat Q = k_therm R sum(dT/dn) dtheta into the wall, from a cubic fit of T along 720
      rays (the same estimator for both runs, which share no body-frame bookkeeping); the force F_x from D/ib_forces.dat. Prints
      the moving-to-fixed ratios of the means over the last --avg of the run, and the RMS of F_moving(t) - F_fixed(t).

Fields are read from restart_data/lustre_<step>.dat (parallel_io = T, one rank) or, without parallel_io, from
p_all/p0/<step>/q_cons_vf<k>.dat (Fortran unformatted records).
"""

import argparse
import glob
import os
import re

import numpy as np


def read_inp(path):
    out = {}
    for line in open(path):
        m = re.match(r"\s*([^=\s]+)\s*=\s*(\S+)", line)
        if m:
            out[m.group(1)] = m.group(2).strip("'\"")
    return out


class Run:
    def __init__(self, path):
        self.path = path
        self.inp = {**read_inp(os.path.join(path, "pre_process.inp")), **read_inp(os.path.join(path, "simulation.inp"))}
        self.m, self.n = int(self.inp["m"]) + 1, int(self.inp["n"]) + 1
        self.x0, self.y0 = self.f("x_domain%beg"), self.f("y_domain%beg")
        self.Lx, self.Ly = self.f("x_domain%end") - self.x0, self.f("y_domain%end") - self.y0
        self.dx, self.dy = self.Lx / self.m, self.Ly / self.n
        self.dt, self.R = self.f("dt"), self.f("patch_ib(1)%radius")
        self.V = self.f("patch_ib(1)%vel(1)", 0.0) if self.f("patch_ib(1)%moving_ibm", 0.0) else 0.0
        self.saves = self.find_saves()

    def f(self, key, default=None):
        return float(self.inp[key].replace("d", "e")) if key in self.inp else default

    def find_saves(self):
        steps = {int(m.group(1)) for p in glob.glob(os.path.join(self.path, "restart_data", "lustre_*.dat")) if (m := re.search(r"lustre_(\d+)\.dat$", p))}
        steps |= {int(m.group(1)) for p in glob.glob(os.path.join(self.path, "p_all", "p0", "*", "q_cons_vf1.dat")) if (m := re.search(r"/(\d+)/q_cons_vf1\.dat$", p))}
        return sorted(steps)

    def centroid(self, step):
        return self.f("patch_ib(1)%x_centroid") + self.V * step * self.dt, self.f("patch_ib(1)%y_centroid")

    def cells(self):
        x = self.x0 + (np.arange(self.m) + 0.5) * self.dx
        y = self.y0 + (np.arange(self.n) + 0.5) * self.dy
        return np.meshgrid(x, y, indexing="xy")  # arrays [y, x]

    def distance(self, step):
        X, Y = self.cells()
        cx, cy = self.centroid(step)
        ddx = (X - cx + 0.5 * self.Lx) % self.Lx - 0.5 * self.Lx
        ddy = (Y - cy + 0.5 * self.Ly) % self.Ly - 0.5 * self.Ly
        return np.hypot(ddx, ddy)

    def fields(self, step):
        """Conservative variables, array [variable, y, x]."""
        lustre = os.path.join(self.path, "restart_data", f"lustre_{step}.dat")
        if os.path.exists(lustre):
            return np.fromfile(lustre, dtype="<f8").reshape(-1, self.n, self.m)
        out = []
        k = 1
        while os.path.exists(p := os.path.join(self.path, "p_all", "p0", str(step), f"q_cons_vf{k}.dat")):
            raw = np.fromfile(p, dtype=np.uint8)
            out.append(np.frombuffer(raw[4:-4].tobytes(), dtype="<f8").reshape(self.n, self.m))  # 4-byte record markers
            k += 1
        return np.array(out)

    def history(self, name="ib_forces.dat"):
        d = np.loadtxt(os.path.join(self.path, "D", name), ndmin=2)
        return d[np.argsort(d[:, 1])]


def mass_budget(path):
    run = Run(path)
    rows = []
    for s in run.saves:
        gas = run.distance(s) > run.R
        rows.append((abs(run.V) * s * run.dt / (2 * run.R) if run.V else s * run.dt, run.fields(s)[0][gas].sum() * run.dx * run.dy))
    d = np.array(rows)
    change = 100.0 * (d[:, 1] / d[0, 1] - 1.0)
    return d[:, 0], change, run


def conducted_heat(run, step):
    """Conducted heat into the isothermal wall [per unit depth]: k_therm * integral of dT/dn."""
    Tw, k = run.f("patch_ib(1)%Twall"), run.f("fluid_pp(1)%k_therm")
    q = run.fields(step)
    T = (q[3] - 0.5 * (q[1] ** 2 + q[2] ** 2) / q[0]) / q[0] / run.f("fluid_pp(1)%cv", 1.0)
    th = np.linspace(0, 2 * np.pi, 720, endpoint=False)
    ds = np.array([1.5, 2.0, 2.5, 3.0]) * run.dx
    pinv = np.linalg.pinv(np.stack([ds, ds**2, ds**3], axis=1))
    cx, cy = run.centroid(step)
    samples = []
    for d in ds:
        fi = ((cx + (run.R + d) * np.cos(th) - run.x0) % run.Lx) / run.dx - 0.5
        fj = ((cy + (run.R + d) * np.sin(th) - run.y0) % run.Ly) / run.dy - 0.5
        i0, j0 = np.floor(fi).astype(int), np.floor(fj).astype(int)
        ax, ay = fi - i0, fj - j0
        i0, j0 = i0 % run.m, j0 % run.n
        i1, j1 = (i0 + 1) % run.m, (j0 + 1) % run.n
        samples.append((1 - ax) * (1 - ay) * T[j0, i0] + ax * (1 - ay) * T[j0, i1] + (1 - ax) * ay * T[j1, i0] + ax * ay * T[j1, i1])
    slope = (pinv @ (np.array(samples) - Tw))[0]
    return k * run.R * slope.sum() * (th[1] - th[0])


def cmd_mass(args):
    for path in args.runs:
        travel, change, run = mass_budget(path)
        print(f"{path}: {run.m} x {run.n} cells, ib_mass_relax = {run.inp.get('ib_mass_relax', 'default')}, {len(travel)} saves")
        print(f"  {'D travelled':>12} {'dM/M [%]':>12}")
        for t, c in zip(travel, change):
            print(f"  {t:12.3f} {c:+12.5f}")
        late = travel >= 1.0
        if late.sum() >= 2:
            rate = np.polyfit(travel[late], change[late], 1)[0]
            print(f"  steady rate after 1 D (least squares): {rate:+.5f} % of the gas per D travelled")


def cmd_twin(args):
    mv, fx = Run(args.moving), Run(args.fixed)
    sel = lambda r: r.saves[int(len(r.saves) * (1 - args.avg)) :]  # noqa: E731
    q = [np.mean([conducted_heat(r, s) for s in sel(r)]) for r in (mv, fx)]
    print(f"conducted heat, moving / fixed: {q[0]:.5e} / {q[1]:.5e}  ratio {q[0] / q[1]:.4f}")
    fm, ff = mv.history(), fx.history()
    t0 = (1 - args.avg) * fm[-1, 1]
    sm = fm[:, 1] >= t0
    t = fm[sm, 1]
    d = fm[sm, 2] - np.interp(t, ff[:, 1], ff[:, 2])
    print(f"F_x, moving / fixed: {fm[sm, 2].mean():.5e} / {ff[ff[:, 1] >= t0, 2].mean():.5e}  ratio {fm[sm, 2].mean() / ff[ff[:, 1] >= t0, 2].mean():.4f}")
    print(
        f"RMS of F_moving - F_fixed: {np.sqrt(np.mean((d - d.mean()) ** 2)) / abs(ff[ff[:, 1] >= t0, 2].mean()) * 100:.3f} % of the fixed mean (mean difference {d.mean() / abs(ff[ff[:, 1] >= t0, 2].mean()) * 100:+.3f} %)"
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("mass")
    s.add_argument("runs", nargs="+")
    s.set_defaults(fn=cmd_mass)
    s = sub.add_parser("twin")
    s.add_argument("--moving", required=True)
    s.add_argument("--fixed", required=True)
    s.add_argument("--avg", type=float, default=0.5, help="Fraction of the run averaged (the last part).")
    s.set_defaults(fn=cmd_twin)
    args = ap.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
