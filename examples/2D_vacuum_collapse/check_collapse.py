#!/usr/bin/env python3
"""Check a finished 2D_vacuum_collapse run: positive density and pressure everywhere at every saved step, and total mass and
energy kept to 1e-4 (the limiter is conservative; the drift is the weak outer rarefaction leaving through the boundaries).

Usage: python3 check_collapse.py <case dir> [--fluids 1|2] [--model 5eq|6eq]
Reads the parallel-I/O restart files restart_data/lustre_<step>.dat (conservative variables, Fortran order).
"""

import argparse
import glob
import os
import re
import sys

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("case_dir")
ap.add_argument("--fluids", type=int, default=1)
ap.add_argument("--model", default="5eq")
a = ap.parse_args()

nf = a.fluids
nvar = 2 * nf + 3 + (nf if a.model == "6eq" else 0)  # alpha_rho, 2 momenta, E, alpha (+ internal energies)
gammas = [1 / 0.4, 1 / (2 / 3)][:nf]  # 1/(gamma - 1)
files = sorted(glob.glob(os.path.join(a.case_dir, "restart_data", "lustre_[0-9]*.dat")), key=lambda f: int(re.findall(r"(\d+)\.dat", f)[0]))
if not files:
    sys.exit("no restart_data/lustre_<step>.dat files")
size = os.path.getsize(files[0])
dtype = np.float64 if np.sqrt(size / 8 / nvar) % 1 == 0 else np.float32  # double- or single-precision build
N = int(round(np.sqrt(size / np.dtype(dtype).itemsize / nvar)))
dA = (1.0 / N) ** 2
ok, base = True, None
print(f"{'step':>6} {'min rho':>11} {'min p':>11} {'mass drift':>11} {'energy drift':>12}")
for f in files:
    q = np.fromfile(f, dtype=dtype).astype(np.float64).reshape(nvar, N, N)
    rho = q[:nf].sum(0)
    E = q[nf + 2]
    alpha = q[nf + 3 : 2 * nf + 3]
    G = sum(alpha[i] * gammas[i] for i in range(nf))
    p = (E - 0.5 * (q[nf] ** 2 + q[nf + 1] ** 2) / rho) / G
    tot = np.array([rho.sum() * dA, E.sum() * dA])
    base = tot if base is None else base
    drift = (tot - base) / base
    step = re.findall(r"(\d+)\.dat", f)[0]
    print(f"{step:>6} {rho.min():11.3e} {p.min():11.3e} {drift[0]:11.1e} {drift[1]:12.1e}")
    ok &= bool(rho.min() > 0 and p.min() > 0 and np.all(np.abs(drift) < 1e-4))
print("PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
