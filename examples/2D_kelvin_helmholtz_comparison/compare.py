#!/usr/bin/env python3
"""Compare validated KH runs at identical grids and nondimensional times."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "viz"
os.environ.setdefault("MPLCONFIGDIR", str(OUT / ".matplotlib"))


def main():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    names = ("ideal_gas", "pr_normal", "pr_supercritical")
    labels = ("Ideal gas, 1 atm", "PR, 1 atm", "PR, 10 MPa")
    runs = []
    # Run each standalone verifier in its own process to avoid ambiguous case imports.
    import subprocess
    import sys

    for name in names:
        folder = ROOT.parent / f"2D_kelvin_helmholtz_{name}"
        subprocess.run([sys.executable, str(folder / "verify.py")], check=True)
        runs.append((json.loads((folder / "viz/metrics.json").read_text()), np.load(folder / "viz/fields.npz")))
    base, data0 = runs[0]
    differences = {}
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), constrained_layout=True)
    growth, ga = plt.subplots(figsize=(6, 4), constrained_layout=True)
    for i, (name, label, (metrics, data)) in enumerate(zip(names, labels, runs)):
        assert (metrics["nx"], metrics["ny"]) == (base["nx"], base["ny"]), "Use matching grids"
        np.testing.assert_allclose(data["history"][:, 0], data0["history"][:, 0], rtol=1e-11, atol=1e-12)
        np.testing.assert_allclose(data["initial"][:4], data0["initial"][:4], rtol=1e-5, atol=1e-6)
        f = data["final"]
        diff = f[:4] - data0["final"][:4]
        differences[name] = {field + "_normalized_rms_difference": float(np.sqrt(np.mean(d * d))) for field, d in zip(("density", "u", "v", "pressure"), diff)}
        differences[name]["transverse_ke_relative_difference"] = float(data["history"][-1, 1] / data0["history"][-1, 1] - 1)
        for row, values, title in ((0, f[0], label), (1, diff[0], "Density difference from ideal")):
            im = axes[row, i].imshow(values, origin="lower", extent=(0, 1, 0, 1), cmap="turbo" if row == 0 else "RdBu_r", vmin=0.8 if row == 0 else -0.3, vmax=2.2 if row == 0 else 0.3)
            axes[row, i].set(title=title, xlabel="x/L", ylabel="y/L")
            fig.colorbar(im, ax=axes[row, i], label="rho / rho_ref" if row == 0 else "Delta(rho / rho_ref)")
        ga.plot(data["history"][:, 0], data["history"][:, 1], label=label, linestyle=("-", "--", "-.")[i])
    fig.suptitle(f"Matched Kelvin–Helmholtz cases: {base['nx']} × {base['ny']}, t U/L = {base['flow_times']:.2f}")
    ga.set(xlabel="t U/L", ylabel="Normalized transverse kinetic energy")
    ga.legend()
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / "density_comparison.png", dpi=180)
    growth.savefig(OUT / "growth.png", dpi=180)
    (OUT / "comparison.json").write_text(json.dumps(dict(grid=[base["nx"], base["ny"]], flow_times=base["flow_times"], differences=differences), indent=2) + "\n")
    print(json.dumps(differences, indent=2))
    print("EOS sensitivity comparison, not a grid-convergence or experimental validation.")


if __name__ == "__main__":
    main()
