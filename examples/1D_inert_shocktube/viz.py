import mfc.viz
import os

import subprocess
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm

from case import sol_L as sol

case = mfc.viz.Case(".")

os.makedirs("viz", exist_ok=True)

# sns.set_theme(style=mfc.viz.generate_cpg_style())

Y_VARS = ["H2", "O2", "H2O", "N2"]

variables = [
    ("rho", "prim.1"),
    ("u_x", "prim.2"),
    ("p", "prim.3"),
    ("E", "cons.3"),
    *[(f"Y_{name}", f"prim.{5 + sol.species_index(name)}") for name in Y_VARS],
    ("T", "prim.15"),
]

for variable in tqdm(variables, desc="Loading Variables"):
    case.load_variable(*variable)

for step in tqdm(case.get_timesteps(), desc="Rendering Frames"):
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))

    def pad_ylim(ylim, pad=0.1):
        return (
            ylim[0] - pad * (ylim[1] - ylim[0]),
            ylim[1] + pad * (ylim[1] - ylim[0]),
        )

    case.plot_step(step, "rho", ax=axes[0, 0])
    axes[0, 0].set_ylim(*pad_ylim(case.get_minmax_time("rho")))
    axes[0, 0].set_ylabel("$\\rho$")
    case.plot_step(step, "u_x", ax=axes[0, 1])
    axes[0, 1].set_ylim(*pad_ylim(case.get_minmax_time("u_x")))
    axes[0, 1].set_ylabel("$u_x$")
    case.plot_step(step, "p", ax=axes[1, 0])
    axes[1, 0].set_ylim(*pad_ylim(case.get_minmax_time("p")))
    axes[1, 0].set_ylabel("$p$")
    for y in Y_VARS:
        case.plot_step(step, f"Y_{y}", ax=axes[1, 1], label=y)
    axes[1, 1].set_ylim(0, 1.1 * max(case.get_minmax_time(f"Y_{y}")[1] for y in Y_VARS))
    axes[1, 1].set_ylabel("$Y_k$")
    case.plot_step(step, "T", ax=axes[1, 2])
    axes[1, 2].set_ylim(*pad_ylim(case.get_minmax_time("T")))
    axes[1, 2].set_ylabel("$T$")
    case.plot_step(step, "E", ax=axes[0, 2])
    axes[0, 2].set_ylim(*pad_ylim(case.get_minmax_time("E")))
    axes[0, 2].set_ylabel("$E$")

    plt.tight_layout()
    plt.savefig(f"viz/{step:06d}.png")
    plt.close()

subprocess.run(
    [
        "ffmpeg",
        "-y",
        "-framerate",
        "60",
        "-pattern_type",
        "glob",
        "-i",
        "viz/*.png",
        "-c:v",
        "libx264",
        "-pix_fmt",
        "yuv420p",
        "viz.mp4",
    ]
)
