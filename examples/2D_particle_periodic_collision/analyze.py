#!/usr/bin/env python3
# Compare the contact force of every pair in the 1- and 2-rank runs against the middle pair of the 1-rank run, and render a
# video. Usage: analyze.py <1-rank case dir> <2-rank case dir> <output dir>
import glob
import importlib.util
import math
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation
from matplotlib.patches import Circle

here = os.path.dirname(os.path.abspath(__file__))
sys.argv, case_argv = [sys.argv[0]], sys.argv[1:]  # case.py parses its own arguments
spec = importlib.util.spec_from_file_location("case", os.path.join(here, "case.py"))
case = importlib.util.module_from_spec(spec)
with open(os.devnull, "w") as devnull:
    stdout, sys.stdout = sys.stdout, devnull
    spec.loader.exec_module(case)
    sys.stdout = stdout

pair_names = ["middle (x)", "x boundary", "y boundary", "corner"]
mass = 4.0 / 3.0 * math.pi * case.radius**3 * case.particle_density
restitution = 0.8


def load_run(case_dir):
    """Per particle ID: time, position (n, 2), velocity (n, 2), and the rank whose file holds each row."""
    rows = []
    for path in sorted(glob.glob(os.path.join(case_dir, "D", "lag_particle_evol_*.dat"))):
        rank = int(path.rsplit("_", 1)[1].split(".")[0])
        for line in open(path):
            fields = line.split()
            try:
                values = [float(v) for v in fields]
            except ValueError:
                continue  # header
            rows.append([rank] + values)
    data = np.array(rows)
    particles = {}
    for pid in np.unique(data[:, 2]).astype(int):
        sel = data[data[:, 2] == pid]
        sel = sel[np.argsort(sel[:, 1], kind="stable")]
        _, first = np.unique(sel[:, 1], return_index=True)
        sel = sel[first]
        particles[pid] = {"t": sel[:, 1], "pos": sel[:, 3:5], "vel": sel[:, 6:8], "rank": sel[:, 0].astype(int)}
    return particles


def pair_history(particles, pair):
    """Normal velocities of both particles and the contact force on particle 1 along its normal (m du/dt)."""
    normal = np.array(case.pairs[pair][1])
    first, second = particles[2 * pair + 1], particles[2 * pair + 2]
    u1, u2 = first["vel"] @ normal, second["vel"] @ normal
    t = first["t"]
    return t, u1, u2, mass * np.gradient(u1, t), mass * np.gradient(u2, t)


def analytic_force(t):
    """Linear spring-dashpot: overlap d'' = -k d - eta d' with relative speed 2v at contact; force on particle 1 = m_eff d''."""
    log_e = math.log(restitution)
    eta = -2.0 * log_e / case.collision_time
    omega = math.pi / case.collision_time
    t_contact = case.gap / (2.0 * case.speed)
    tau = t - t_contact
    v_rel = 2.0 * case.speed
    overlap = v_rel / omega * np.exp(-0.5 * eta * tau) * np.sin(omega * tau)
    rate = v_rel * np.exp(-0.5 * eta * tau) * (np.cos(omega * tau) - 0.5 * eta / omega * np.sin(omega * tau))
    force = -0.5 * mass * ((math.pi**2 + log_e**2) / case.collision_time**2 * overlap + eta * rate)
    return np.where((tau > 0) & (tau < case.collision_time), force, 0.0)


def main(one_rank_dir, two_rank_dir, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    runs = {"1 rank": load_run(one_rank_dir), "2 ranks": load_run(two_rank_dir)}
    t_ref, u_ref, _, f_ref, _ = pair_history(runs["1 rank"], 0)
    f_scale = np.abs(f_ref).max()

    print(f"{'run':8s} {'pair':12s} {'max|dF|/Fmax':>13s} {'max|du|/v':>10s} {'F1+F2 /Fmax':>12s} {'rebound/v':>10s}")
    histories = {}
    for run, particles in runs.items():
        for pair, name in enumerate(pair_names):
            t, u1, u2, f1, f2 = pair_history(particles, pair)
            histories[(run, pair)] = (t, f1)
            n = min(len(t), len(t_ref))
            print(
                f"{run:8s} {name:12s} {np.abs(f1[:n] - f_ref[:n]).max() / f_scale:13.3e} "
                f"{np.abs(u1[:n] - u_ref[:n]).max() / case.speed:10.3e} {np.abs(f1 + f2).max() / f_scale:12.3e} "
                f"{-u1[-1] / case.speed:10.4f}"
            )

    t_fine = np.linspace(t_ref[0], t_ref[-1], 4000)
    figure, axis = plt.subplots(figsize=(7, 4.5))
    axis.plot(t_fine, analytic_force(t_fine), color="0.6", lw=4, label="analytic spring-dashpot")
    for (run, pair), (t, f1) in histories.items():
        axis.plot(t, f1, ls="-" if run == "1 rank" else "--", lw=1.2, label=f"{run}: {pair_names[pair]}")
    axis.set_xlim(0.045, 0.06)
    axis.set_xlabel("time")
    axis.set_ylabel("contact force on particle 1 along its normal")
    axis.legend(fontsize=7, ncol=2)
    figure.tight_layout()
    figure.savefig(os.path.join(out_dir, "contact_force.png"), dpi=150)

    render_video(runs["2 ranks"], histories, os.path.join(out_dir, "collisions.mp4"))


def render_video(particles, histories, path):
    """Whole box with the rank owning each particle, a zoom on every contact (periodically unwrapped), and the forces."""
    figure = plt.figure(figsize=(15, 6.5))
    grid = figure.add_gridspec(2, 5)
    box = figure.add_subplot(grid[:, :2])
    zooms = [figure.add_subplot(grid[i // 2, 2 + i % 2]) for i in range(4)]
    force_axis = figure.add_subplot(grid[:, 4])
    for (run, pair), (t, f1) in histories.items():
        force_axis.plot(t, f1, ls="-" if run == "1 rank" else "--", lw=1, label=f"{run}: {pair_names[pair]}")
    force_axis.set_xlim(0.049, 0.057)
    force_axis.set_xlabel("time")
    force_axis.set_title("contact force along each pair's normal", fontsize=9)
    force_axis.legend(fontsize=6)
    cursor = force_axis.axvline(0.0, color="k", lw=0.8)
    rank_colors = ["tab:blue", "tab:orange"]
    times = particles[1]["t"]

    def draw(frame):
        box.clear()
        box.set_xlim(0, 1)
        box.set_ylim(0, 1)
        box.set_aspect("equal")
        box.axvline(0.5, color="0.5", ls="--", lw=1)
        box.set_title(f"2 ranks, t = {times[frame]:.4f}  (color = owning rank, dashed = rank boundary)", fontsize=9)
        for pid, history in particles.items():
            x, y = history["pos"][frame]
            box.plot(x, y, "o", ms=7, color=rank_colors[history["rank"][frame]])
        for pair, ((cx, cy), _) in enumerate(case.pairs):
            label_x = (cx + 0.02) % 1.0 if cx < 0.5 else cx - 0.02
            label_y = (cy + 0.04) % 1.0 if cy < 0.5 else cy - 0.04
            box.annotate(pair_names[pair], (label_x, label_y), fontsize=8, ha="center")

        for pair, axis in enumerate(zooms):
            axis.clear()
            (cx, cy), _ = case.pairs[pair]
            half = 0.02
            axis.set_xlim(cx - half, cx + half)
            axis.set_ylim(cy - half, cy + half)
            axis.set_aspect("equal")
            axis.set_title(pair_names[pair], fontsize=9)
            axis.tick_params(labelsize=6)
            axis.set_xticks([cx - half, cx, cx + half])
            axis.set_yticks([cy - half, cy, cy + half])
            if cx == 0.5:
                axis.axvline(cx, color="0.5", ls="--", lw=1)
            if cx in (0.0, 1.0):
                axis.axvline(cx, color="k", lw=0.8)
            if cy in (0.0, 1.0):
                axis.axhline(cy, color="k", lw=0.8)
            for pid in (2 * pair + 1, 2 * pair + 2):
                x, y = particles[pid]["pos"][frame]
                x += round(cx - x)  # the periodic image nearest the contact
                y += round(cy - y)
                axis.add_patch(Circle((x, y), case.radius, color=rank_colors[particles[pid]["rank"][frame]], alpha=0.8))
        cursor.set_xdata([times[frame]])
        return []

    frames = range(0, len(times), 4)
    writer = animation.FFMpegWriter(fps=20, bitrate=3000)
    animation.FuncAnimation(figure, draw, frames=frames, blit=False).save(path, writer=writer, dpi=110)


if __name__ == "__main__":
    if len(case_argv) != 3:
        sys.exit("usage: analyze.py <1-rank case dir> <2-rank case dir> <output dir>")
    main(*case_argv)
