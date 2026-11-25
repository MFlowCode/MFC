#!/usr/bin/env python3
"""Compare the three standalone runs without modifying their cases or results."""

import csv
import importlib.util
import json
import os
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "viz"
os.environ.setdefault("MPLCONFIGDIR", str(OUT / ".matplotlib"))

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

FOLDERS = ("2D_boundary_layer_ideal_gas", "2D_boundary_layer_pr_normal", "2D_boundary_layer_pr_supercritical")
LABELS = ("Ideal gas, 1 atm", "PR, 1 atm", "PR, 7.5 MPa")
STYLES = ("-", "--", "-.")


def load_run(folder):
    directory = ROOT.parent / folder
    spec = importlib.util.spec_from_file_location(folder, directory / "case.py")
    case = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(case)
    s = case.reference()
    inp = {}
    for line in (directory / "simulation.inp").read_text().splitlines():
        if "=" in line:
            k, value = line.split("=", 1)
            inp[k.strip()] = value.strip()
    # Refuse to interpret old outputs using edited grid, timing or EOS settings.
    for key, value in case.case().items():
        if key in inp and isinstance(value, (int, float)):
            assert np.isclose(float(inp[key]), value, rtol=1e-12, atol=0), f"{folder}: case.py disagrees with simulation.inp for {key}"
    assert inp["parallel_io"] == "T" and inp["precision"] == "2"
    assert int(inp["fluid_pp(1)%eos"]) == (2 if case.IDEAL else 6)
    data_dir = directory / "restart_data"
    steps = sorted(int(p.stem.split("_")[1]) for p in data_dir.glob("lustre_*.dat") if p.stem.split("_")[1].isdigit())
    assert steps == list(range(0, int(inp["t_step_stop"]) + 1, int(inp["t_step_save"]))), f"Missing or stale saves: {folder}"
    xcb = np.fromfile(data_dir / "lustre_x_cb.dat", dtype=np.float64)
    ycb = np.fromfile(data_dir / "lustre_y_cb.dat", dtype=np.float64)
    assert xcb.size == case.NX + 1 and ycb.size == case.NY + 1
    assert np.isclose(xcb[-1] - xcb[0], s["length"], rtol=1e-12, atol=0)
    assert np.isclose(ycb[-1] - ycb[0], s["height"], rtol=1e-12, atol=0)
    x, y = (xcb[1:] + xcb[:-1]) / 2, (ycb[1:] + ycb[:-1]) / 2

    def read(step):
        data = np.fromfile(data_dir / f"lustre_{step}.dat", dtype=np.float64)
        assert data.size == 5 * case.NX * case.NY, f"Unexpected MPI layout: {folder}, {step}"
        rho, mx, my, energy, alpha = data.reshape(5, case.NY, case.NX)
        assert np.isfinite(data).all() and np.min(rho) > 0 and np.max(abs(alpha - 1)) < 1e-8
        u, v = mx / rho, my / rho
        t = case.temperature_from_energy(rho, energy / rho - 0.5 * (u * u + v * v), case.IDEAL)
        p, _, cv, cp, sound = case.properties(rho, t, case.IDEAL)
        assert all(np.isfinite(a).all() and np.min(a) > 0 for a in (t, p, np.asarray(cv), np.asarray(cp), sound))
        return dict(rho=rho, u=u, v=v, t=t, p=p)

    initial = read(0)
    for key, expected in (("rho", s["rho"]), ("u", s["velocity"]), ("t", s["temperature"]), ("p", s["pressure"])):
        assert np.allclose(initial[key], expected, rtol=1e-9), f"{folder}: initial {key} disagrees with case.py"
    final = read(steps[-1])
    if folder.endswith("supercritical"):
        assert final["t"].min() > case.TC and final["p"].min() > case.PC, "Final field left the supercritical region"
    times = np.asarray(steps) * float(inp["dt"]) * s["velocity"] / s["length"]
    old_index = np.argmin(abs(times - (times[-1] - 0.5)))
    old = read(steps[old_index])
    ix = int(np.argmin(abs(x / s["length"] - 0.75)))
    delta_t = s["temperature"] - s["wall_temperature"]

    def wall(f):
        cf = 2 * case.MU * f["u"][0] / y[0] / (s["rho"] * s["velocity"] ** 2)
        qw = s["conductivity"] * (f["t"][0] - s["wall_temperature"]) / y[0]
        nu = qw * x / (s["conductivity"] * delta_t)
        return cf, qw, nu

    cf, qw, nu = wall(final)
    history = []
    for step, time in zip(steps, times):
        f = read(step)
        c, q, n = wall(f)
        history.append((step, time, c[ix], q[ix], n[ix]))
    history = np.asarray(history)
    metrics = dict(
        freestream_temperature_K=s["temperature"],
        wall_temperature_K=s["wall_temperature"],
        pressure_Pa=s["pressure"],
        density_kg_m3=s["rho"],
        velocity_m_s=s["velocity"],
        sound_speed_m_s=s["sound"],
        plate_length_m=s["length"],
        conductivity_W_m_K=s["conductivity"],
        cp_J_kg_K=s["cp"],
        compressibility_Z=s["pressure"] / (s["rho"] * case.R * s["temperature"]),
        mach=case.MACH,
        reynolds=case.REYNOLDS,
        prandtl=case.PRANDTL,
        saves=len(steps),
        nx=case.NX,
        ny=case.NY,
        last_step=steps[-1],
        flow_times=times[-1],
        physical_time_s=steps[-1] * float(inp["dt"]),
        x_over_L=x[ix] / s["length"],
        skin_friction=cf[ix],
        heat_flux_to_wall_W_m2=qw[ix],
        local_Nusselt=nu[ix],
        velocity_max_change_over_U=float(np.max(abs(final["u"] - old["u"])) / s["velocity"]),
        temperature_max_change_over_deltaT=float(np.max(abs(final["t"] - old["t"])) / delta_t),
        cf_relative_change=float(abs(cf[ix] - history[old_index, 2]) / abs(cf[ix])),
        Nu_relative_change=float(abs(nu[ix] - history[old_index, 4]) / abs(nu[ix])),
        comparison_interval_flow_times=times[-1] - times[old_index],
        min_temperature_K=float(final["t"].min()),
        min_pressure_Pa=float(final["p"].min()),
    )
    return dict(folder=folder, case=case, s=s, x=x, y=y, xcb=xcb, ycb=ycb, ix=ix, final=final, cf=cf, qw=qw, nu=nu, history=history, metrics=metrics)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    runs = [load_run(f) for f in FOLDERS]
    for r in runs[1:]:
        for key in ("freestream_temperature_K", "wall_temperature_K", "mach", "reynolds", "prandtl", "nx", "ny", "flow_times"):
            assert np.isclose(r["metrics"][key], runs[0]["metrics"][key], rtol=1e-12), f"Cases differ in {key}"
    metrics = {r["folder"]: r["metrics"] for r in runs}
    differences = {}
    for r in runs[1:]:
        a, b = runs[0], r
        differences[r["folder"]] = {}
        for field, scale in (("u", "velocity"), ("t", None), ("rho", "rho")):

            def normalized(run):
                f, s = run["final"][field], run["s"]
                return f / s[scale] if scale else (f - s["wall_temperature"]) / (s["temperature"] - s["wall_temperature"])

            diff = normalized(b) - normalized(a)
            np.testing.assert_allclose(a["y"] / a["s"]["length"], b["y"] / b["s"]["length"], rtol=1e-12)
            differences[r["folder"]][field + "_profile_max_difference"] = float(np.max(abs(diff[:, a["ix"]])))
            differences[r["folder"]][field + "_field_max_difference"] = float(np.max(abs(diff)))
    (OUT / "metrics.json").write_text(json.dumps(dict(cases=metrics, normalized_differences_from_ideal=differences), indent=2) + "\n")
    with (OUT / "metrics.csv").open("w") as f:
        writer = csv.DictWriter(f, fieldnames=["case"] + list(runs[0]["metrics"]))
        writer.writeheader()
        for r in runs:
            writer.writerow(dict(case=r["folder"], **r["metrics"]))
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5), constrained_layout=True)
    fig2, ax2 = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    fig3, ax3 = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    fig4, ax4 = plt.subplots(2, 3, figsize=(13, 7), constrained_layout=True)
    for col, (r, label, style) in enumerate(zip(runs, LABELS, STYLES)):
        s, f, ix = r["s"], r["final"], r["ix"]
        eta = r["y"] * np.sqrt(s["rho"] * s["velocity"] / (r["case"].MU * r["x"][ix]))
        profiles = (f["u"][:, ix] / s["velocity"], (f["t"][:, ix] - s["wall_temperature"]) / (s["temperature"] - s["wall_temperature"]), f["rho"][:, ix] / s["rho"])
        for ax, profile in zip(axes, profiles):
            ax.plot(profile, eta, style, label=label)
        mask = (r["x"] / s["length"] >= 0.1) & (r["x"] / s["length"] <= 0.9)
        for ax, values in zip(ax2, (r["cf"], r["nu"], r["qw"] / 1000)):
            ax.plot(r["x"][mask] / s["length"], values[mask], style, label=label)
        hist = r["history"]
        keep = hist[:, 1] >= 0.1
        for ax, index in zip(ax3, (2, 4)):
            ax.plot(hist[keep, 1], hist[keep, index], style, label=label)
        for row, values in enumerate((f["u"] / s["velocity"], (f["t"] - s["wall_temperature"]) / (s["temperature"] - s["wall_temperature"]))):
            mesh = ax4[row, col].pcolormesh(r["xcb"] / s["length"], r["ycb"] / s["length"], values, vmin=0, vmax=1.1 if row == 0 else 1, cmap="viridis")
            ax4[row, col].set(xlabel="x/L", ylabel="y/L", ylim=(0, 0.12), title=label)
            fig4.colorbar(mesh, ax=ax4[row, col], label="u/U∞" if row == 0 else "(T−Tw)/(T∞−Tw)")
        np.savetxt(OUT / (r["folder"] + "_profiles.csv"), np.column_stack((eta, *profiles)), delimiter=",", header="eta,u_over_U,theta,rho_over_rhoinf")
        np.savetxt(OUT / (r["folder"] + "_history.csv"), hist, delimiter=",", header="step,tU_over_L,Cf,qwall_W_m2,Nu_x")
    for ax, label in zip(axes, ("u/U∞", "(T−Tw)/(T∞−Tw)", "ρ/ρ∞")):
        ax.set(xlabel=label, ylabel="η = y √(ρ∞ U∞ / μx)", ylim=(0, 10))
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.suptitle(f"Final profiles at x/L = {runs[0]['metrics']['x_over_L']:.3f}; tU∞/L = {runs[0]['metrics']['flow_times']:.3f}")
    for ax, label in zip(ax2, ("Estimated Cf", "Estimated Nu_x", "Heat flux to wall [kW/m²]")):
        ax.set(xlabel="x/L", ylabel=label)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    for ax, label in zip(ax3, ("Estimated Cf at x/L ≈ 0.75", "Estimated Nu_x at x/L ≈ 0.75")):
        ax.set(xlabel="Flow-through times, tU∞/L", ylabel=label)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    for fig, name in ((fig, "profiles"), (fig2, "wall_transfer"), (fig3, "settling"), (fig4, "fields")):
        fig.savefig(OUT / (name + ".png"), dpi=180)
    write_report(runs, differences)
    print(json.dumps(dict(cases=metrics, normalized_differences_from_ideal=differences), indent=2))


def write_report(runs, differences):
    """Write a standalone report alongside the figures and machine-readable data."""

    def table(headers, rows):
        return (
            "<table><thead><tr>"
            + "".join(f"<th>{h}</th>" for h in headers)
            + "</tr></thead><tbody>"
            + "".join("<tr>" + "".join(f"<td>{v}</td>" for v in row) + "</tr>" for row in rows)
            + "</tbody></table>"
        )

    rows = []
    changes = []
    for r, label in zip(runs, LABELS):
        m = r["metrics"]
        rows.append(
            [
                label,
                f"{m['density_kg_m3']:.4f}",
                f"{m['velocity_m_s']:.3f}",
                f"{1000*m['plate_length_m']:.5f}",
                f"{m['compressibility_Z']:.6f}",
                f"{m['conductivity_W_m_K']:.6f}",
                f"{m['skin_friction']:.6f}",
                f"{m['local_Nusselt']:.3f}",
                f"{m['heat_flux_to_wall_W_m2']/1000:.3f}",
            ]
        )
        changes.append(
            [label, f"{100*m['velocity_max_change_over_U']:.2f}%", f"{100*m['temperature_max_change_over_deltaT']:.2f}%", f"{100*m['cf_relative_change']:.2f}%", f"{100*m['Nu_relative_change']:.2f}%"]
        )
    diff_rows = []
    for folder, label in zip(FOLDERS[1:], LABELS[1:]):
        d = differences[folder]
        diff_rows.append([label, *[f"{100*d[k+'_profile_max_difference']:.4f}%" for k in ("u", "t", "rho")]])
    a, b = runs[0]["metrics"], runs[2]["metrics"]
    ratio = b["heat_flux_to_wall_W_m2"] / a["heat_flux_to_wall_W_m2"]
    length_ratio = a["plate_length_m"] / b["plate_length_m"]
    nu_change = 100 * (b["local_Nusselt"] / a["local_Nusselt"] - 1)
    m = runs[0]["metrics"]
    date = datetime.now(ZoneInfo("America/New_York")).strftime("%d %B %Y")
    html = f"""<!doctype html><html lang="en"><meta charset="utf-8"><title>MFC boundary-layer comparison</title>
<style>body{{font:16px/1.55 system-ui,sans-serif;max-width:1200px;margin:40px auto;padding:0 24px;color:#17232e}}
h1,h2{{line-height:1.2}}table{{border-collapse:collapse;width:100%;font-size:14px;margin:18px 0}}td,th{{padding:9px;border-bottom:1px solid #ccd5df;text-align:right}}td:first-child,th:first-child{{text-align:left}}th{{background:#eef3f8}}img{{width:100%;height:auto}}.note{{padding:18px;background:#fff2ce;border-left:5px solid #b57713}}code{{background:#edf1f5;padding:2px 4px}}a{{color:#075f9a}}</style>
<h1>Comparison of the three boundary-layer runs</h1><p>Analysis of existing MPI-I/O outputs, {date}. No simulation settings or saved run data were modified.</p>
<p>All cases: {m['nx']} × {m['ny']} cells; T∞ = {m['freestream_temperature_K']:.0f} K; Tw = {m['wall_temperature_K']:.0f} K; Mach {m["mach"]:g}; Re<sub>L</sub> = {m["reynolds"]:g}; freestream Pr = {m["prandtl"]:g}.
The final save is step {m['last_step']}, tU∞/L = {m['flow_times']:.6f}. Profiles and local wall values are sampled at x/L = {m['x_over_L']:.6f}.</p>
<div class="note"><strong>These are evolving, matched-time results.</strong> The maximum velocity changes by about 4.5% of U∞ over the final half flow-through time. Local Nu changes by about 8%.
The runs have not met the earlier settling criteria of 1% velocity change and 2% temperature change over this interval. Small frame-to-frame changes with 100 saves do not establish steady state.</div>
<h2>Physical conditions and final wall estimates</h2>
{table(['Case','ρ∞ [kg/m³]','U∞ [m/s]','L [mm]','Z∞','k [W/(m K)]','Cf estimate','Nu_x estimate','q to wall [kW/m²]'],rows)}
<p>Z = p/(ρRT). At the 7.5 MPa freestream Z = {b['compressibility_Z']:.4f}, so the pressure-density relation at this temperature is still fairly close to the ideal-gas value of one.
The most visible normalized profile difference is in the density increase near the cooled wall. The high-pressure final field remains above the specified nitrogen critical temperature and pressure.</p>
<h2>Profile differences from ideal gas</h2>
{table(['Case','max |Δ(u/U∞)|','max |Δ((T−Tw)/(T∞−Tw))|','max |Δ(ρ/ρ∞)|'],diff_rows)}
<p>Percentages here are differences expressed as percent of each normalization scale, not pointwise relative errors near the wall. All comparisons use corresponding dimensionless positions.</p><img src="profiles.png" alt="Velocity, temperature and density profiles">
<h2>Interpreting the heat flux</h2>
<p>The dimensional high-pressure heat flux is {ratio:.1f} times the ideal-gas result, but its plate is {length_ratio:.1f} times shorter.
Its conductivity also differs because k = μ cp∞/Pr. At corresponding x/L, q = Nu_x k (T∞−Tw)/x, whereas Nu_x changes by only {nu_change:.2f}%.
The large dimensional flux ratio therefore cannot be described as an intrinsic supercritical heat-transfer enhancement at fixed geometry.</p>
<p>The three cases prescribe μ = 1.8×10⁻⁵ Pa s and constant k within each run. They are controlled EOS comparisons at matched Reynolds/Mach/Prandtl numbers, with different physical sizes and durations, rather than measurements using a state-dependent nitrogen transport model.</p>
<img src="wall_transfer.png" alt="Wall shear and heat transfer along the plate">
<h2>Settling over the last half flow-through time</h2>
{table(['Case','max Δu / U∞','max ΔT / (T∞−Tw)','relative Cf change','relative Nu_x change'],changes)}
<p>Cf and Nu relative changes are measured at x/L ≈ 0.75 and normalized by their final values. Velocity and temperature changes are maxima over the whole domain.</p>
<img src="settling.png" alt="Wall quantities versus flow-through time">
<p>Extend the runs to at least four flow-through times and recheck the same half-flow-time interval; continue longer if needed. A separate grid-convergence check is needed for quantitative wall shear and heat flux. The near-wall velocity shows small spatial oscillations, and the shear estimate amplifies differences in the first-cell velocity.</p>
<h2>Two-dimensional layers</h2><img src="fields.png" alt="Near-wall velocity and temperature fields for all three cases">
<h2>Method and reproducibility</h2>
<p>The script reads five double-precision conservative fields from the combined MPI-I/O restart files, reconstructs temperature and pressure using each case's EOS, and checks grid dimensions, saved-step coverage, initial states, and agreement between case.py and simulation.inp.
All stored states used in the time histories have finite, positive density, temperature and pressure. This checks these outputs; it is not a serial-versus-MPI equivalence test.</p>
<p>Wall quantities are first-cell estimates: τw ≈ μ u₁/y₁; Cf = 2τw/(ρ∞U∞²); q to wall ≈ k(T₁−Tw)/y₁; Nu_x = q x/[k(T∞−Tw)]. They are not directly exported numerical face fluxes.</p>
<p>Reproduce from the repository root: <code>build/venv/bin/python examples/2D_boundary_layer_comparison/compare.py</code>.</p>
<p><a href="metrics.csv">Metrics CSV</a> · <a href="metrics.json">Metrics JSON</a>. Profile and time-history CSV files are also saved beside this report.</p></html>"""
    (OUT / "comparison.html").write_text(html)


if __name__ == "__main__":
    main()
