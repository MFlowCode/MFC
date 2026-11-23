#!/usr/bin/env python3
"""
Lumped-temperature cooling of the cylinder (--thermal_bc 3) against the decay time its fixed-wall Nusselt number implies.

    check_lumped_cooling.py LUMPED_DIR --cold COLD_DIR --hot HOT_DIR [--tmin 20]

LUMPED_DIR is a run with --thermal_bc 3 --tw_ratio 1.1 (and --rho_cs); COLD_DIR and HOT_DIR are fixed-wall twins (--tw_ratio 1 and 1.1,
same Re and grid) as for check_cylinder.py. A body of heat capacity rho_s c_s in a stream whose heat-transfer coefficient is
Nu k_f / D relaxes to the adiabatic-wall temperature T_aw with

    tau = rho_s c_s D^2 / (4 Nu k_f)       (2D cylinder: mass per area rho_s c_s D / 4)

Prints the model tau and T_aw from the twins, the fit T_a' + A exp(-t/tau') of Twall(t) from D/ib_heat.dat over t U / D > --tmin with
the implied Nu_fit and its difference from the fixed-wall Nu, and the ratio of the lumped body's heat rate to the one the twins'
linear law Q(Tw) gives at the same wall temperature (the fit-free check).
"""

import argparse
import os

import numpy as np
from check_cylinder import case_info, film, history, nusselt, window_mean


def fit_exponential(t, T):
    """Least-squares T_a + A exp(-t/tau): scan tau, solve (T_a, A) linearly."""
    best = None
    for tau in np.geomspace(0.05 * (t[-1] - t[0]), 20 * (t[-1] - t[0]), 4000):
        M = np.stack([np.ones_like(t), np.exp(-(t - t[0]) / tau)], axis=1)
        c, res, *_ = np.linalg.lstsq(M, T, rcond=None)
        r = float(np.sum((M @ c - T) ** 2))
        if best is None or r < best[0]:
            best = (r, tau, c)
    return best[1], best[2]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("lumped_dir")
    ap.add_argument("--cold", required=True)
    ap.add_argument("--hot", required=True)
    ap.add_argument("--tmin", type=float, default=20.0, help="Fit window start [D/U].")
    ap.add_argument("--avg", type=float, default=0.2, help="Averaging window of the twins, fraction of the run.")
    args = ap.parse_args()

    il, ic, ih = (case_info(d) for d in (args.lumped_dir, args.cold, args.hot))
    T_inf = ic["Tw"]
    nu = nusselt(args.cold, args.hot, ic, ih, args.avg)["Nu"]
    k_f, _ = film(il, il["Tw"], T_inf)
    inp = open(os.path.join(args.lumped_dir, "simulation.inp")).read()
    rho_cs = float(next(l.split("=")[1] for l in inp.splitlines() if "patch_ib(1)%rho_solid" in l).replace("d", "e")) * float(
        next(l.split("=")[1] for l in inp.splitlines() if "patch_ib(1)%cp_solid" in l).replace("d", "e")
    )
    tau = rho_cs * il["D"] ** 2 / (4 * nu * k_f)
    # adiabatic-wall temperature: heat into the body is linear in Tw and vanishes at T_aw
    qc, qh = (window_mean(h[:, 1], h[:, 4], args.avg) for h in (history(args.cold, "ib_heat.dat"), history(args.hot, "ib_heat.dat")))
    Taw = T_inf + (qc / (qc - qh)) * (ih["Tw"] - T_inf)
    h = history(args.lumped_dir, "ib_heat.dat")
    t, Tw = h[:, 1] * il["U"] / il["D"], h[:, 5]
    m = t >= args.tmin
    tau_fit, (Ta, A) = fit_exponential(t[m], Tw[m])
    tau_model = tau * il["U"] / il["D"]
    nu_fit = nu * tau_model / tau_fit
    print(f"Re {il['Re']:.1f}, {il['cells_per_D']:.0f} cells/D, rho_s c_s = {rho_cs:.3e} J/m^3/K, t U/D to {t[-1]:.1f}")
    print(f"fixed-wall twins: Nu = {nu:.3f}, T_aw = {Taw:.2f} K")
    print(f"model: tau U/D = {tau_model:.1f}")
    print(f"fit  : tau' U/D = {tau_fit:.1f}, T_a' = {Ta:.2f} K, Nu_fit = {nu_fit:.3f} ({nu_fit / nu - 1:+.1%} against the twins)")
    # Heat into the body is linear in Tw: the twins give Q(Tw); compare the lumped run's heat rate with it at its own Tw(t)
    q_pred = qc + (qh - qc) * (Tw - T_inf) / (ih["Tw"] - T_inf)
    ratio = h[:, 4][m] / q_pred[m]
    print(f"heat rate / twins' law at the same Tw over the window: mean {np.average(ratio, weights=np.gradient(t[m])):.3f}, first {ratio[0]:.3f}, last {ratio[-1]:.3f}")
    print(f"Twall: {Tw[0]:.2f} K -> {Tw[-1]:.2f} K; model exp(-t/tau) from T_aw: {Taw + (Tw[0] - Taw) * np.exp(-t[-1] * il['D'] / il['U'] / tau):.2f} K")


if __name__ == "__main__":
    main()
