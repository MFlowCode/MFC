#!/usr/bin/env python3
"""Validate periodic KH output and save normalized fields and growth histories."""

import argparse
import json
import os
from pathlib import Path

import case
import numpy as np
from scipy.io import FortranFile

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "viz"
os.environ.setdefault("MPLCONFIGDIR", str(OUT / ".matplotlib"))


def record(path):
    with FortranFile(path, "r") as stream:
        return stream.read_reals(np.float64)


def analyze():
    inp = dict(line.strip().split(" = ", 1) for line in (ROOT / "simulation.inp").read_text().splitlines() if " = " in line)
    nx, ny = int(inp["m"]) + 1, int(inp["n"]) + 1
    parallel = inp["parallel_io"] == "T"
    assert int(inp["fluid_pp(1)%eos"]) == (2 if case.IDEAL else 6), "EOS differs from case.py"
    for key in ("fluid_pp(1)%gamma", "fluid_pp(1)%cv", "fluid_pp(1)%pr_tc", "fluid_pp(1)%pr_pc", "fluid_pp(1)%pr_r", "fluid_pp(1)%pr_omega"):
        if key in case.case():
            assert np.isclose(float(inp[key]), case.case()[key], rtol=1e-12, atol=0), f"Changed {key}; rerun first"
    s = case.reference()
    for axis, size in (("x", nx), ("y", ny)):
        edges = np.fromfile(ROOT / "restart_data" / f"lustre_{axis}_cb.dat") if parallel else record(ROOT / "p_all/p0/0" / f"{axis}_cb.dat")
        assert len(edges) == size + 1 and np.allclose(edges, np.linspace(0, s["length"], size + 1), atol=1e-14), "Grid mismatch"
    steps = list(range(0, int(inp["t_step_stop"]) + 1, int(inp["t_step_save"])))
    history, initial, totals0, scale = [], None, None, None
    min_t, min_p = np.inf, np.inf
    max_conservation = 0.0
    snapshots = []
    for step in steps:
        if parallel:
            raw = np.fromfile(ROOT / "restart_data" / f"lustre_{step}.dat")
            assert raw.size == 5 * nx * ny, "Expected combined MPI I/O (file_per_process=F)"
            q = raw.reshape(5, ny, nx)
        else:
            q = np.array([record(ROOT / "p_all/p0" / str(step) / f"q_cons_vf{i}.dat").reshape(ny, nx) for i in range(1, 6)])
        rho, mx, my, energy, alpha = q
        assert np.isfinite(q).all() and rho.min() > 0 and np.max(abs(alpha - 1)) < 1e-8
        u, v = mx / rho, my / rho
        t = case.temperature_from_energy(rho, energy / rho - 0.5 * (u * u + v * v), case.IDEAL)
        p, _, cv, cp, sound = case.properties(rho, t, case.IDEAL)
        assert all(np.isfinite(a).all() and np.min(a) > 0 for a in (t, p, np.asarray(cv), np.asarray(cp), sound)), "Invalid thermodynamic state"
        min_t, min_p = min(min_t, float(t.min())), min(min_p, float(p.min()))
        fields = np.array([rho / s["rho"], u / s["velocity"], v / s["velocity"], p / case.PRESSURE, t])
        totals = q[:4].mean(axis=(1, 2))
        if step == 0:
            y = (np.arange(ny) + 0.5) / ny
            x = (np.arange(nx) + 0.5) / nx
            layer = 0.5 * (1 - np.tanh(np.cos(2 * np.pi * y) / (2 * np.pi * case.SHEAR_WIDTH)))
            np.testing.assert_allclose(rho, np.broadcast_to(s["rho"] * (1 + layer[:, None]), rho.shape), rtol=2e-6)
            np.testing.assert_allclose(u, np.broadcast_to(s["velocity"] * (layer[:, None] - 0.5), u.shape), atol=1e-5 * s["velocity"])
            expected_v = case.PERTURBATION * s["velocity"] * np.sin(4 * np.pi * x)[None, :] * (np.exp(-(((y - 0.25) / 0.05) ** 2)) + np.exp(-(((y - 0.75) / 0.05) ** 2)))[:, None]
            np.testing.assert_allclose(v, expected_v, atol=1e-6 * s["velocity"])
            np.testing.assert_allclose(p, case.PRESSURE, rtol=2e-6)
            initial, totals0 = fields, totals
            scale = np.array([totals[0], totals[0] * s["velocity"], totals[0] * s["velocity"], max(abs(totals[3]), case.PRESSURE)])
        error = np.max(abs(totals - totals0) / scale)
        max_conservation = max(max_conservation, float(error))
        assert error < 1e-8, f"Periodic conservation failed at {step}: {error}"
        time = step * float(inp["dt"]) / s["time"]
        transverse_ke = np.mean(0.5 * fields[0] * fields[2] ** 2)
        mode = abs(np.mean(fields[2] * np.exp(-4j * np.pi * (np.arange(nx) + 0.5) / nx)[None, :]))
        history.append([time, transverse_ke, mode, t.min(), p.min(), error])
        snapshots.append(fields[0].astype(np.float32))
    if case.SUPERCRITICAL:
        assert min_t > case.TC and min_p > case.PC, "A saved state left the supercritical region"
    report = dict(
        nx=nx,
        ny=ny,
        flow_times=time,
        saves=len(steps),
        min_temperature_K=min_t,
        min_pressure_Pa=min_p,
        max_conservation_error=max_conservation,
        initial_temperature_range_K=[float(initial[4].min()), float(initial[4].max())],
        final_transverse_ke=history[-1][1],
        transverse_ke_ratio=history[-1][1] / history[0][1],
        reference=s,
        ideal=case.IDEAL,
        supercritical=case.SUPERCRITICAL,
    )
    OUT.mkdir(exist_ok=True)
    np.savetxt(OUT / "history.csv", history, delimiter=",", header="t_star,transverse_ke,mode2,min_T_K,min_p_Pa,conservation_error", comments="")
    np.savez_compressed(OUT / "fields.npz", initial=initial, final=fields, density_frames=snapshots, history=history)
    (OUT / "metrics.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return fields, np.asarray(history), snapshots


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--animate", action="store_true", help="Also write a density GIF")
    args = parser.parse_args()
    fields, history, frames = analyze()
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6), constrained_layout=True)
    for ax, values, label in zip(axes[:2], (fields[0], fields[4]), ("Density / reference density", "Temperature [K]")):
        im = ax.imshow(values, origin="lower", extent=(0, 1, 0, 1), cmap="turbo")
        fig.colorbar(im, ax=ax)
        ax.set(xlabel="x/L", ylabel="y/L", title=label)
    axes[2].plot(history[:, 0], history[:, 1])
    axes[2].set(xlabel="t U/L", ylabel="Transverse kinetic energy (normalized)")
    fig.savefig(OUT / "summary.png", dpi=160)
    plt.close(fig)
    if args.animate:
        from matplotlib.animation import FuncAnimation, PillowWriter

        fig, ax = plt.subplots(figsize=(5, 5))
        im = ax.imshow(frames[0], origin="lower", extent=(0, 1, 0, 1), cmap="turbo", vmin=0.8, vmax=2.2)
        fig.colorbar(im, ax=ax, label="Density / reference density")

        def update(i):
            im.set_data(frames[i])
            ax.set_title(f"t U/L = {history[i, 0]:.2f}")
            return (im,)

        FuncAnimation(fig, update, frames=len(frames)).save(OUT / "density.gif", writer=PillowWriter(fps=15))
    print("PASS: finite stable states and periodic conservation; see viz/")


if __name__ == "__main__":
    main()
