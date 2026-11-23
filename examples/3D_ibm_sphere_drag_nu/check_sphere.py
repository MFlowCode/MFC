#!/usr/bin/env python3
"""
Drag and Nusselt number of 3D_ibm_sphere_drag_nu against correlations.

    check_sphere.py CASE_DIR [--hot HOT_DIR] [--avg 0.2]

CASE_DIR is a run with --tw_ratio 1 (the drag run); HOT_DIR a run of the same case with --tw_ratio 1.1. It reads
pre_process.inp, simulation.inp, D/ib_forces.dat and D/ib_heat.dat (plain text, no post_process needed) and prints

  Cd_IB    <F_x>/(rho U^2 pi D^2 / 8) from the face-flux force of ib_force_wrt, averaged over the last --avg of the run,
           against Schiller-Naumann Cd = 24/Re (1 + 0.15 Re^0.687)
  Nu       -(Q_hot - Q_cold) / (pi D^2 (Tw - T_inf)) D / k_f with Q the heat the body takes (D/ib_heat.dat column 5); the cold
           twin removes the viscous-heating flux of the wall. k_f is the film conductivity from Cantera. Against Ranz-Marshall
           and Whitaker at the free-stream Re and the film Prandtl number.

If silo_hdf5/ is present (run with --silo and post_process) it also prints, from the saves in the averaging window,
  Cd_CV    the same drag from a steady momentum balance on cubes of half-width 1 and 2 D, which uses only the fields
  L/D      the first u = 0 crossing on the wake axis behind the rear stagnation point, last save (Johnson and Patel: 0.88 at Re 100)
"""

import argparse
import glob
import os
import re

import numpy as np

# Plane-by-plane trapezoid rule needs only 1D integrals, which np.trapz / np.trapezoid provide under either name
trapezoid = getattr(np, "trapezoid", None) or np.trapz


def cd_schiller_naumann(re):
    return 24.0 / re * (1.0 + 0.15 * re**0.687)


def nu_ranz_marshall(re, pr):
    return 2.0 + 0.6 * re**0.5 * pr ** (1.0 / 3.0)


def nu_whitaker(re, pr):
    return 2.0 + (0.4 * re**0.5 + 0.06 * re ** (2.0 / 3.0)) * pr**0.4


def read_inp(path):
    out = {}
    for line in open(path):
        m = re.match(r"\s*([^=\s]+)\s*=\s*(\S+)", line)
        if m:
            out[m.group(1)] = m.group(2).strip("'\"")
    return out


def case_info(run):
    inp = {**read_inp(os.path.join(run, "pre_process.inp")), **read_inp(os.path.join(run, "simulation.inp"))}
    f = lambda k: float(inp[k].replace("d", "e").replace("D", "e"))  # noqa: E731
    info = dict(
        rho=f("patch_icpp(1)%alpha_rho(1)"),
        U=f("patch_icpp(1)%vel(1)"),
        p=f("patch_icpp(1)%pres"),
        mu=1.0 / f("fluid_pp(1)%Re(1)"),
        R=f("patch_ib(1)%radius"),
        center=[f(f"patch_ib(1)%{c}_centroid") for c in "xyz"],
        Tw=f("patch_ib(1)%Twall"),
        dx=(f("x_domain%end") - f("x_domain%beg")) / (int(inp["m"]) + 1),
    )
    info["D"] = 2 * info["R"]
    info["Re"] = info["rho"] * info["U"] * info["D"] / info["mu"]
    info["cells_per_D"] = info["D"] / info["dx"]
    info["q"] = 0.5 * info["rho"] * info["U"] ** 2 * np.pi * info["R"] ** 2
    return info


def history(run, name):
    d = np.loadtxt(os.path.join(run, "D", name), ndmin=2)
    return d[np.argsort(d[:, 1])]


def window_mean(t, y, avg):
    m = t >= t[-1] - avg * (t[-1] - t[0])
    return float(np.average(y[m], weights=np.gradient(t)[m])) if m.sum() > 1 else float(y[-1])


def film(info, Tw, T_inf):
    import cantera as ct

    g = ct.Solution("h2o2.yaml")
    g.TPX = 0.5 * (Tw + T_inf), info["p"], "O2:0.21,N2:0.79"
    return g.thermal_conductivity, g.viscosity * g.cp_mass / g.thermal_conductivity


def nusselt(cold, hot, info_c, info_h, avg):
    T_inf = info_c["Tw"]  # the cold twin's wall is at the free-stream temperature
    k_f, pr = film(info_h, info_h["Tw"], T_inf)
    q = [window_mean(h[:, 1], h[:, 4], avg) for h in (history(hot, "ib_heat.dat"), history(cold, "ib_heat.dat"))]
    nu = -(q[0] - q[1]) / (np.pi * info_h["D"] ** 2) * info_h["D"] / (k_f * (info_h["Tw"] - T_inf))
    re = info_h["Re"]
    return {"Nu": nu, "Ranz-Marshall": nu_ranz_marshall(re, pr), "Whitaker": nu_whitaker(re, pr), "Pr_f": pr}


def fields(run, info, avg):
    from mfc.viz.silo_reader import assemble_silo

    steps = sorted(int(os.path.basename(f)[:-5]) for f in glob.glob(os.path.join(run, "silo_hdf5", "p0", "*.silo")) if os.path.basename(f)[:-5].isdigit())
    if not steps:
        raise SystemExit(f"{run}/silo_hdf5/p0 holds no *.silo saves (run with --silo and post_process)")
    use = [s for s in steps if s >= (1 - avg) * steps[-1]] or steps[-1:]
    cvs = []
    for s in use:
        d = assemble_silo(run, s)
        xs = [np.asarray(c).squeeze() for c in (d.x_cc, d.y_cc, d.z_cc)]
        V = {k: np.asarray(v).squeeze() for k, v in d.variables.items()}
        cvs.append([cv_drag(xs, V, info, h) for h in (1.0, 2.0)])
    return {"Cd_CV(h=1,2)": np.mean(cvs, axis=0), "L/D": recirculation(xs, V["vel1"], info)}


def cv_drag(xs, V, info, half):
    """Cd from a steady momentum balance on the cube |x_i - c_i| <= half D (frozen unsteady and acoustic terms)."""
    idx = [np.searchsorted(x, c + np.array([-half, half]) * info["D"]) for x, c in zip(xs, info["center"])]
    box = tuple(slice(i0, i1 + 1) for i0, i1 in idx)
    xb = [x[s] for x, s in zip(xs, box)]
    vel = [V[f"vel{i + 1}"][box] for i in range(3)]
    rho, p, mu = V["rho"][box], V["pres"][box] - info["p"], info["mu"]
    g = [np.gradient(v, *xb) for v in vel]
    div = g[0][0] + g[1][1] + g[2][2]
    total = 0.0
    for j in range(3):  # face normal
        flux = rho * vel[0] * vel[j] - mu * (g[0][j] + g[j][0]) + (p + 2 / 3 * mu * div) * (j == 0)
        other = [a for a in range(3) if a != j]
        for end, sign in ((-1, 1.0), (0, -1.0)):
            face = np.take(flux, end, axis=j)
            total += sign * trapezoid(trapezoid(face, xb[other[1]], axis=1), xb[other[0]])
    return -total / info["q"]


def recirculation(xs, u, info):
    iy, iz = (int(np.argmin(np.abs(xs[a] - info["center"][a]))) for a in (1, 2))
    xa, ua = xs[0] - info["center"][0], u[:, iy, iz]
    m = xa > info["R"] + 0.5 * info["dx"]
    xa, ua = xa[m], ua[m]
    if ua[0] >= 0:
        return 0.0
    i = int(np.argmax(ua > 0))
    xr = xa[i - 1] - ua[i - 1] * (xa[i] - xa[i - 1]) / (ua[i] - ua[i - 1])
    return float((xr - info["R"]) / info["D"])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case_dir")
    ap.add_argument("--hot", help="Run directory of the heated twin (--tw_ratio 1.1).")
    ap.add_argument("--avg", type=float, default=0.2, help="Averaging window, fraction of the run.")
    args = ap.parse_args()

    info = case_info(args.case_dir)
    print(f"Re_D {info['Re']:.2f}  {info['cells_per_D']:.1f} cells/D")
    a = history(args.case_dir, "ib_forces.dat")
    cd = window_mean(a[:, 1], a[:, 2] / info["q"], args.avg)
    sn = cd_schiller_naumann(info["Re"])
    print(f"  t U/D {a[-1, 1] * info['U'] / info['D']:.1f}  Cd_IB {cd:.4f}  Schiller-Naumann {sn:.4f} ({100 * (cd / sn - 1):+.1f} %)")
    if os.path.isdir(os.path.join(args.case_dir, "silo_hdf5")):
        f = fields(args.case_dir, info, args.avg)
        print("  Cd_CV (half-width 1, 2 D): " + ", ".join(f"{c:.4f}" for c in f["Cd_CV(h=1,2)"]) + f"   L/D {f['L/D']:.3f}")
    if args.hot:
        nu = nusselt(args.case_dir, args.hot, info, case_info(args.hot), args.avg)
        print(
            f"  Nu {nu['Nu']:.3f}  Ranz-Marshall {nu['Ranz-Marshall']:.3f} ({100 * (nu['Nu'] / nu['Ranz-Marshall'] - 1):+.1f} %)  "
            f"Whitaker {nu['Whitaker']:.3f} ({100 * (nu['Nu'] / nu['Whitaker'] - 1):+.1f} %)  Pr_f {nu['Pr_f']:.3f}"
        )


if __name__ == "__main__":
    main()
