#!/usr/bin/env python3
"""Check this standalone boundary layer and plot its saved profiles."""

import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "viz" / ".matplotlib"))

import matplotlib  # noqa: E402
import numpy as np
from case import MU, NX, NY, PARALLEL_IO, PC, TC, properties, reference, temperature_from_energy
from scipy.integrate import solve_ivp
from scipy.io import FortranFile

VARIANTS = ("normal_pr",)

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def read_record(path):
    with FortranFile(path, "r") as stream:
        return stream.read_reals(np.float64)


def grid(axis):
    if PARALLEL_IO == "T":
        return np.fromfile(ROOT / "restart_data" / f"lustre_{axis}_cb.dat", dtype=np.float64)
    return read_record(ROOT / "p_all" / "p0" / "0" / f"{axis}_cb.dat")


def fields(variant, step):
    if PARALLEL_IO == "T":
        path = ROOT / "restart_data" / f"lustre_{step}.dat"
        data = np.fromfile(path, dtype=np.float64)
        if data.size != 5 * NY * NX:
            raise ValueError(f"{path}: expected five double-precision fields on a {NX} x {NY} grid; use combined MPI I/O")
        q = data.reshape(5, NY, NX)
    else:
        directory = ROOT / "p_all" / "p0" / str(step)
        q = [read_record(directory / f"q_cons_vf{i}.dat").reshape(NY, NX) for i in range(1, 6)]
    rho, mx, my, energy, alpha = q
    assert all(np.isfinite(a).all() for a in q), f"{variant}: nonfinite conservative state"
    assert np.min(rho) > 0, f"{variant}: nonpositive density"
    assert np.max(abs(alpha - 1)) < 1e-8, f"{variant}: single-fluid volume fraction drift"
    u, v = mx / rho, my / rho
    internal = energy / rho - 0.5 * (u**2 + v**2)
    temperature = temperature_from_energy(rho, internal, variant == "normal_ideal")
    pressure, _, cv, cp, sound = properties(rho, temperature, variant == "normal_ideal")
    assert all(np.isfinite(a).all() and np.min(a) > 0 for a in (temperature, pressure, np.asarray(cv), np.asarray(cp), sound))
    return u, v, temperature, pressure


def blasius(eta):
    solution = solve_ivp(lambda _, f: (f[1], f[2], -0.5 * f[0] * f[2]), (0, max(12.0, float(np.max(eta)))), (0, 0, 0.332057336215), rtol=1e-10, atol=1e-12, dense_output=True)
    return solution.sol(eta)[1]


def main(variants=VARIANTS):
    destination = ROOT / "viz"
    destination.mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    field_fig, field_axes = plt.subplots(2, len(variants), figsize=(4 * len(variants), 6), squeeze=False, constrained_layout=True)
    reports, profiles = {}, {}
    for column_index, variant in enumerate(variants):
        s = reference()
        if PARALLEL_IO == "F":
            inlet = read_record(ROOT / "p_all" / "p0" / "0" / "bc_buffers.dat").reshape(NY, 6)
            assert np.allclose(inlet[:, 5], s["temperature"], rtol=1e-12), "Incorrect saved inlet temperature"
        else:
            _, _, initial_temperature, _ = fields(variant, 0)
            assert np.allclose(initial_temperature, s["temperature"], rtol=1e-10), "Incorrect initial temperature"
        assert grid("x").size == NX + 1 and grid("y").size == NY + 1, "Grid differs from case.py"
        u, v, temperature, pressure = fields(variant, s["steps"])
        old_u, _, old_t, _ = fields(variant, s["steps"] - s["save_steps"])
        delta_t = s["temperature"] - s["wall_temperature"]
        steady_u = np.max(abs(u - old_u)) / s["velocity"]
        steady_t = np.max(abs(temperature - old_t)) / delta_t
        assert steady_u < 0.01 and steady_t < 0.02, f"{variant}: not steady: velocity={steady_u:.3e}, temperature={steady_t:.3e}"
        assert abs(np.mean(u[-4:]) / s["velocity"] - 1) < 0.05, f"{variant}: outer velocity disturbed"
        assert abs(np.mean(temperature[-4:]) - s["temperature"]) / delta_t < 0.05, f"{variant}: outer temperature disturbed"
        assert np.min(u) > -0.02 * s["velocity"], f"{variant}: reversed flow"
        if variant == "supercritical_pr":
            assert np.min(temperature) > TC and np.min(pressure) > PC, "State left the supercritical region"
        ix = int(0.75 * NX)
        x = (ix + 0.5) * s["length"] / NX
        y_edges = grid("y")
        y = 0.5 * (y_edges[1:] + y_edges[:-1])
        eta = y * np.sqrt(s["rho"] * s["velocity"] / (MU * x))
        velocity = u[:, ix] / s["velocity"]
        theta = (temperature[:, ix] - s["wall_temperature"]) / delta_t
        # Extrapolate the first two cell centers to the physical wall.
        wall_u = (y[1] * velocity[0] - y[0] * velocity[1]) / (y[1] - y[0])
        wall_t = (y[1] * theta[0] - y[0] * theta[1]) / (y[1] - y[0])
        assert abs(wall_u) < 0.025 and abs(wall_t) < 0.05, f"{variant}: wall extrapolation mismatch {wall_u=:.3e}, {wall_t=:.3e}"
        assert 0 < theta[0] < 0.3 and theta[-1] > 0.95, f"{variant}: missing thermal boundary layer"
        assert velocity[0] < 0.2 and velocity[-1] > 0.95, f"{variant}: missing velocity boundary layer"
        # First crossing of 99% freestream speed avoids integrating far-field disturbances.
        thickness = []
        for column in (NX // 4, ix):
            profile = u[:, column] / s["velocity"]
            crossings = np.flatnonzero(profile >= 0.99)
            assert len(crossings) and crossings[0] > 0, f"{variant}: unresolved boundary-layer edge"
            edge = crossings[0]
            thickness.append(np.interp(0.99, profile[edge - 1 : edge + 1], y[edge - 1 : edge + 1]))
        assert thickness[1] > 1.2 * thickness[0] > 0, f"{variant}: boundary layer does not grow"
        cf = 2 * MU * u[0, ix] / y[0] / (s["rho"] * s["velocity"] ** 2)
        cf_blasius = 0.664 / np.sqrt(s["rho"] * s["velocity"] * x / MU)
        profile_error = np.sqrt(np.mean((velocity[eta < 8] - blasius(eta[eta < 8])) ** 2))
        # Blasius is only an approximate reference for these cooled compressible cases.
        if variant != "supercritical_pr":
            assert profile_error < 0.08 and 0.7 < cf / cf_blasius < 1.4, f"{variant}: large deviation from dilute-gas Blasius reference"
        report = dict(
            velocity_change=steady_u,
            temperature_change=steady_t,
            wall_velocity_extrapolation=wall_u,
            wall_temperature_extrapolation=wall_t,
            blasius_profile_rms=profile_error,
            thickness_99_over_length=thickness[1] / s["length"],
            skin_friction=cf,
            skin_friction_over_blasius=cf / cf_blasius,
            heat_flux_into_fluid=-s["conductivity"] * (temperature[0, ix] - s["wall_temperature"]) / y[0],
        )
        reports[variant] = {k: float(value) for k, value in report.items()}
        profiles[variant] = (velocity, theta)
        print(f"{variant}: {json.dumps(reports[variant], indent=2)}")
        style = {"normal_ideal": "-", "normal_pr": "--", "supercritical_pr": "-."}[variant]
        axes[0].plot(velocity, eta, label=variant, linestyle=style)
        axes[1].plot(theta, eta, label=variant, linestyle=style)
        x_edges = grid("x")
        for row, values in enumerate((u / s["velocity"], (temperature - s["wall_temperature"]) / delta_t)):
            ax = field_axes[row, column_index]
            mesh = ax.pcolormesh(x_edges / s["length"], y_edges / s["length"], values, vmin=0, vmax=1.1 if row == 0 else 1, cmap="viridis", shading="flat")
            ax.set(xlabel="x / L", ylabel="y / L", ylim=(0, 0.12), title=variant.replace("_", " "))
            field_fig.colorbar(mesh, ax=ax, label="u / U∞" if row == 0 else "(T − Tw) / (T∞ − Tw)")
        np.savetxt(destination / f"{variant}_profile.csv", np.column_stack((eta, velocity, theta)), delimiter=",", header="eta,u_over_U,theta")
    eta = np.linspace(0, 10, 200)
    axes[0].plot(blasius(eta), eta, "k--", label="Blasius (approximate reference)")
    for ax, label in zip(axes, ("u / U∞", "(T − Tw) / (T∞ − Tw)")):
        ax.set(xlabel=label, ylabel="η = y √(ρ∞ U∞ / μx)", ylim=(0, 10))
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
    fig.savefig(destination / "profiles.png", dpi=180)
    field_fig.savefig(destination / "fields.png", dpi=180)
    (destination / "verification.json").write_text(json.dumps(reports, indent=2) + "\n")
    print(f"PASS: {len(variants)} boundary-layer smoke test(s)")
    print(f"Plot: {destination / 'profiles.png'}")


if __name__ == "__main__":
    main()
