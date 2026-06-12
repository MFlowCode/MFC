#!/usr/bin/env python3
"""
Validation A (n = 1 Newtonian equivalence) for examples/2D_ibm_poiseuille_nn.

A power-law fluid with nn = 1 and tau0 = 0 is analytically a Newtonian fluid
with mu = K, and for a single fluid the non-Newtonian mixture arithmetic reduces
to the same 1/K — so a run with IBM_NN_MODE=newtonian (mu = 0.02, NN code path
OFF) and a run with IBM_NN_MODE=nn1 (K = 0.02, nn = 1, NN code path ON,
including the per-stencil-sample IBM viscosity) must produce the same fields to
near round-off. Both modes use the SAME fixed dt, so the comparison is at
matched time step with identical dt histories.

Usage:
    python check_equivalence.py <newtonian_run_dir> <nn1_run_dir>

Each directory must contain restart_data/lustre_*.dat from a run of case.py
with the corresponding IBM_NN_MODE. Reports the max abs and relative L2
difference of the u and v fields at the last common save; PASS if rel L2 < 1e-8.
"""

import glob
import os
import sys

import numpy as np

# Must match case.py
NX = 25
NY = 96
NVAR = 5  # alpha_rho(1), mom_x, mom_y, E, alpha(1)


def last_save(run_dir):
    dats = sorted(
        glob.glob(os.path.join(run_dir, "restart_data", "lustre_[0-9]*.dat")),
        key=lambda p: int(os.path.basename(p).split("_")[1].split(".")[0]),
    )
    if not dats:
        raise SystemExit(f"No restart files in {run_dir}/restart_data")
    return dats[-1]


def read_uv(path):
    raw = np.fromfile(path, dtype=np.float64).reshape((NVAR, NY, NX))
    return raw[1] / raw[0], raw[2] / raw[0]  # u, v


def main():
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    f_newt, f_nn1 = last_save(sys.argv[1]), last_save(sys.argv[2])
    if os.path.basename(f_newt) != os.path.basename(f_nn1):
        raise SystemExit(f"Save mismatch: {f_newt} vs {f_nn1}")
    print(f"Comparing {f_newt}\n versus   {f_nn1}")

    u_a, v_a = read_uv(f_newt)
    u_b, v_b = read_uv(f_nn1)

    max_abs = max(np.max(np.abs(u_a - u_b)), np.max(np.abs(v_a - v_b)))
    rel_l2 = np.sqrt((np.sum((u_a - u_b) ** 2) + np.sum((v_a - v_b) ** 2)) / (np.sum(u_a**2) + np.sum(v_a**2)))
    print(f"u_max (newtonian run)              : {u_a.max():.6e}")
    print(f"max abs velocity difference        : {max_abs:.3e}")
    print(f"relative L2 velocity difference    : {rel_l2:.3e}")
    print("PASS" if rel_l2 < 1e-8 else "FAIL (> 1e-8: possible IBM mu_eff bug)")


if __name__ == "__main__":
    main()
