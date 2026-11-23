#!/usr/bin/env python3
"""
Drag, Strouhal number, wake length and Nusselt number of 2D_ibm_cylinder_drag_nu against references.

    check_cylinder.py CASE_DIR [--hot HOT_DIR] [--avg 0.2]

CASE_DIR is a run with --tw_ratio 1 (the drag run); HOT_DIR a run of the same case with --tw_ratio 1.1. It reads
pre_process.inp, simulation.inp, D/ib_forces.dat and D/ib_heat.dat (plain text, no post_process needed) and prints

  Cd_IB    <F_x>/(rho U^2 D / 2) from the face-flux force of ib_force_wrt, averaged over the last --avg of the run
  Cl_amp, St   lift amplitude and Strouhal number of the lift history (second half of the run)
  Nu       -(Q_hot - Q_cold) D / (pi D k_f (Tw - T_inf)) with Q the heat the body takes (D/ib_heat.dat column 5); the
           cold twin removes the viscous-heating flux of the wall. k_f is the film conductivity from Cantera.

If silo_hdf5/ is present (run with --silo and post_process) it also prints, from the saves in the averaging window,
  Cd_CV    the same drag from a steady momentum balance on squares of half-width 1, 2 and 4 D, which uses only the fields
  L/D      the first u = 0 crossing on the wake centreline behind the rear stagnation point, last save
"""

import argparse
import glob
import os
import re

import numpy as np

trapezoid = getattr(np, "trapezoid", None) or np.trapz  # np.trapezoid needs NumPy >= 2

# Reference values as tabulated in the papers (see the validation page for the citations)
CD_REF = {20: {"Dennis & Chang": 2.05, "Fornberg": 2.00}, 40: {"Dennis & Chang": 1.52, "Fornberg": 1.50}, 100: {"Henderson": 1.35, "Park et al.": 1.33}}
LD_REF = {20: 0.94, 40: 2.35}  # Dennis & Chang


def st_williamson(re):
    return -3.3265 / re + 0.1816 + 1.6e-4 * re


def nu_lange(re):
    return 0.082 * re**0.5 + 0.734 * re ** (0.05 + 0.226 * re**0.085)


def nu_churchill_bernstein(re, pr):
    return 0.3 + 0.62 * re**0.5 * pr ** (1 / 3) / (1 + (0.4 / pr) ** (2 / 3)) ** 0.25 * (1 + (re / 282000.0) ** 0.625) ** 0.8


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
        center=[f("patch_ib(1)%x_centroid"), f("patch_ib(1)%y_centroid")],
        Tw=f("patch_ib(1)%Twall"),
        dx=(f("x_domain%end") - f("x_domain%beg")) / (int(inp["m"]) + 1),
    )
    info["D"] = 2 * info["R"]
    info["Re"] = info["rho"] * info["U"] * info["D"] / info["mu"]
    info["cells_per_D"] = info["D"] / info["dx"]
    info["q"] = 0.5 * info["rho"] * info["U"] ** 2 * info["D"]
    return info


def history(run, name):
    d = np.loadtxt(os.path.join(run, "D", name), ndmin=2)
    return d[np.argsort(d[:, 1])]


def window_mean(t, y, avg):
    m = t >= t[-1] - avg * (t[-1] - t[0])
    return float(np.average(y[m], weights=np.gradient(t)[m])) if m.sum() > 1 else float(y[-1])


def forces(run, info, avg):
    a = history(run, "ib_forces.dat")
    t, cd, cl = a[:, 1], a[:, 2] / info["q"], a[:, 3] / info["q"]
    late = t >= t[-1] - avg * (t[-1] - t[0])
    out = {"t*": t[-1] * info["U"] / info["D"], "Cd_IB": window_mean(t, cd, avg), "Cl_amp": 0.5 * float(cl[late].max() - cl[late].min())}
    h = t >= 0.5 * t[-1]
    if h.sum() > 16 and out["Cl_amp"] > 1.0e-2:  # a steady wake has no shedding frequency
        tu = np.linspace(t[h][0], t[h][-1], h.sum())
        y = np.interp(tu, t[h], cl[h] - cl[h].mean())
        spec = np.abs(np.fft.rfft(y * np.hanning(len(y))))
        k = int(np.argmax(spec[1:])) + 1
        if k < len(spec) - 1:
            y0, y1, y2 = np.log(spec[k - 1 : k + 2] + 1e-300)
            k += 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)
        out["St"] = k / (len(tu) * (tu[1] - tu[0])) * info["D"] / info["U"]
    return out


def film(info, Tw, T_inf):
    import cantera as ct

    g = ct.Solution("h2o2.yaml")
    g.TPX = 0.5 * (Tw + T_inf), info["p"], "O2:0.21,N2:0.79"
    return g.thermal_conductivity, g.viscosity * g.cp_mass / g.thermal_conductivity


def nusselt(cold, hot, info_c, info_h, avg):
    T_inf = info_c["Tw"]  # the cold twin's wall is at the free-stream temperature
    k_f, pr = film(info_h, info_h["Tw"], T_inf)
    q = [window_mean(h[:, 1], h[:, 4], avg) for h in (history(hot, "ib_heat.dat"), history(cold, "ib_heat.dat"))]
    nu = -(q[0] - q[1]) / (np.pi * info_h["D"]) * info_h["D"] / (k_f * (info_h["Tw"] - T_inf))
    re = info_h["Re"]
    return {"Nu": nu, "Lange": nu_lange(re), "Churchill-Bernstein": nu_churchill_bernstein(re, pr), "Pr_f": pr}


def fields(run, info, avg):
    from mfc.viz.silo_reader import assemble_silo

    steps = sorted(int(os.path.basename(f)[:-5]) for f in glob.glob(os.path.join(run, "silo_hdf5", "p0", "*.silo")) if os.path.basename(f)[:-5].isdigit())
    if not steps:
        raise SystemExit(f"{run}/silo_hdf5/p0 holds no *.silo saves (run with --silo and post_process)")
    use = [s for s in steps if s >= (1 - avg) * steps[-1]] or steps[-1:]
    cvs = []
    for s in use:
        d = assemble_silo(run, s)
        xs = [np.asarray(c).squeeze() for c in (d.x_cc, d.y_cc)]
        V = {k: np.asarray(v).squeeze() for k, v in d.variables.items()}
        cvs.append([cv_drag(xs, V, info, h) for h in (1.0, 2.0, 4.0)])
    ld = recirculation(xs, V["vel1"], info)
    return {"Cd_CV(h=1,2,4)": np.mean(cvs, axis=0), "L/D": ld}


def cv_drag(xs, V, info, half):
    """Cd from a steady momentum balance on the square |x_i - c_i| <= half D (frozen unsteady and acoustic terms)."""
    idx = [np.searchsorted(x, c + np.array([-half, half]) * info["D"]) for x, c in zip(xs, info["center"])]
    box = tuple(slice(i0, i1 + 1) for i0, i1 in idx)
    xb = [x[s] for x, s in zip(xs, box)]
    vel = [V[f"vel{i + 1}"][box] for i in range(2)]
    rho, p, mu = V["rho"][box], V["pres"][box] - info["p"], info["mu"]
    g = [np.gradient(v, *xb) for v in vel]
    div = g[0][0] + g[1][1]
    total = 0.0
    for j in range(2):  # face normal
        flux = rho * vel[0] * vel[j] - mu * (g[0][j] + g[j][0]) + (p + 2 / 3 * mu * div) * (j == 0)
        for end, sign in ((-1, 1.0), (0, -1.0)):
            total += sign * trapezoid(np.take(flux, end, axis=j), xb[1 - j])
    return -total / info["q"]


def recirculation(xs, u, info):
    iy = int(np.argmin(np.abs(xs[1] - info["center"][1])))
    xa, ua = xs[0] - info["center"][0], u[:, iy]
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
    re_key = int(round(info["Re"]))
    print(f"Re_D {info['Re']:.2f}  {info['cells_per_D']:.1f} cells/D")
    out = forces(args.case_dir, info, args.avg)
    for k, v in out.items():
        print(f"  {k:<10s} {v:.4g}")
    for name, v in CD_REF.get(re_key, {}).items():
        print(f"  Cd_IB vs {name} {v}: {100 * (out['Cd_IB'] / v - 1):+.1f} %")
    if "St" in out and re_key == 100:
        print(f"  St vs Williamson {st_williamson(info['Re']):.4f}: {100 * (out['St'] / st_williamson(info['Re']) - 1):+.1f} %")
    if os.path.isdir(os.path.join(args.case_dir, "silo_hdf5")):
        f = fields(args.case_dir, info, args.avg)
        print("  Cd_CV (half-width 1, 2, 4 D): " + ", ".join(f"{c:.4f}" for c in f["Cd_CV(h=1,2,4)"]))
        if re_key in LD_REF:  # the wake of a shedding cylinder has no steady length
            print(f"  L/D {f['L/D']:.3f} (Dennis & Chang {LD_REF[re_key]})")
    if args.hot:
        nu = nusselt(args.case_dir, args.hot, info, case_info(args.hot), args.avg)
        print(
            f"  Nu {nu['Nu']:.3f}  Lange {nu['Lange']:.3f} ({100 * (nu['Nu'] / nu['Lange'] - 1):+.1f} %)  Churchill-Bernstein {nu['Churchill-Bernstein']:.3f} ({100 * (nu['Nu'] / nu['Churchill-Bernstein'] - 1):+.1f} %)  Pr_f {nu['Pr_f']:.3f}"
        )


if __name__ == "__main__":
    main()
