#!/usr/bin/env python3
"""Plot the max-over-ranks wall time per step of each phase against rank count.

Reads one or more phase_time_data.dat files written by simulation with phase_timing_wrt = T
(a directory is read as <dir>/phase_time_data.dat). Rows from all files are pooled; when
several runs share a rank count, the last one read wins.
"""

import argparse
import os
from collections import defaultdict

import matplotlib.pyplot as plt


def read_rows(paths):
    data = defaultdict(dict)  # phase -> {ranks: (mean, max)}
    for path in paths:
        if os.path.isdir(path):
            path = os.path.join(path, "phase_time_data.dat")
        with open(path) as f:
            next(f)
            for line in f:
                ranks, _, _, _, t_mean, t_max, phase = line.split()
                data[phase][int(ranks)] = (float(t_mean), float(t_max))
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("paths", nargs="+", help="phase_time_data.dat files or the case directories holding them")
    parser.add_argument("--phases", nargs="+", help="phases to plot (default: --top largest)")
    parser.add_argument("--top", type=int, default=10, help="plot the N phases with the largest max time at the largest rank count")
    parser.add_argument("--mean", action="store_true", help="also plot the mean over ranks (dashed)")
    parser.add_argument("--relative", action="store_true", help="divide each phase by its time at the smallest rank count")
    parser.add_argument("-o", "--output", default="phase_timings.png", help="output image")
    args = parser.parse_args()

    data = read_rows(args.paths)
    phases = args.phases or sorted(data, key=lambda ph: data[ph][max(data[ph])][1], reverse=True)[: args.top]

    fig, ax = plt.subplots(figsize=(9, 6))
    for phase in phases:
        ranks = sorted(data[phase])
        t_mean, t_max = zip(*(data[phase][r] for r in ranks))
        scale = t_max[0] if args.relative else 1.0
        (line,) = ax.plot(ranks, [t / scale for t in t_max], "o-", label=phase)
        if args.mean:
            ax.plot(ranks, [t / scale for t in t_mean], "--", color=line.get_color())

    ax.set_xscale("log", base=2)
    ax.set_xlabel("Ranks")
    ax.set_ylabel("max time / time at fewest ranks" if args.relative else "max over ranks [s/step]")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize="small")
    fig.tight_layout()
    fig.savefig(args.output, dpi=150)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
