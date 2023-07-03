#!/usr/bin/env python3

import os
import glob
import math
import argparse
import seaborn as sns
import pandas  as pd
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Visualize the results of a case")
parser.add_argument('case_dir', type=str, help="Path to the case directory")
parser.add_argument('fps',      type=int, help="Frames per second for the video")

ARGS     = vars(parser.parse_args())
DIRPATH  = os.path.abspath(ARGS['case_dir'])
OUTDIR   = os.path.join(DIRPATH, 'viz')
CASENAME = os.path.basename(DIRPATH)

TARGET_ASPECT_RATIO = 16 / 9

PLOT_PADDING = 0.04
PLOT_DIMS    = (15, 10)

if not os.path.isdir(OUTDIR):
    os.mkdir(OUTDIR)

VARS = {
    'Density':  ['prim.1', ['Density'] , '#ed4337'],
    'Velocity': ['prim.2', ['Velocity'], '#3561f2'],
    'Pressure': ['prim.3', ['Pressure'], '#28ed87'],
}

DATA = {}

print(f"Importing data from {os.path.relpath(DIRPATH, os.getcwd())}...")

for name, cfg in VARS.items():
    print(f"* {name} ({cfg[0]}):")
    print(f"  * Reading data files...")

    dfs = {}
    for f in glob.glob(os.path.join(DIRPATH, 'D', f'{cfg[0]}.*.*.dat')):     
        proc, t_step = int(f.split('.')[-3]), int(f.split('.')[-2])
        
        if t_step != dfs:
            dfs[t_step] = []
                
        dfs[t_step].append(pd.read_csv(f, sep='\s+', header=None, names=['x'] + cfg[1]))
    
    print(f"  * Concatenating processors...")
    for t_step in dfs.keys():
        dfs[t_step] = pd.concat(dfs[t_step])
    
    print(f"  * Merging across timesteps...")
    for t_step, df in dfs.items():    
        if t_step not in DATA:
            DATA[t_step] = df
        else:
            DATA[t_step] = pd.merge(DATA[t_step], df)

PLOTLIMITS = {}

print(f"* Calculating min/max...")
for name in list(VARS.keys()) + ['x']:
    print(f" * {name}...")
    
    minmax = (math.inf, -math.inf)
    
    for df in DATA.values():
        minmax = (
            min(minmax[0], df[name].min()),
            max(minmax[1], df[name].max())
        )
    
    PLOTLIMITS[name] = (
        minmax[0] - PLOT_PADDING * (minmax[1] - minmax[0]),
        minmax[1] + PLOT_PADDING * (minmax[1] - minmax[0])
    )

def plot(name: str, t_step: int, ax=None):
    sns.lineplot(
        data=DATA[t_step], x='x', y=name,
        color=VARS[name][2], linewidth=1, alpha=0.8,
        ax=ax)
    
    ax.set_xlim(PLOTLIMITS['x'][0],  PLOTLIMITS['x'][1])
    ax.set_ylim(PLOTLIMITS[name][0], PLOTLIMITS[name][1])

def save_filename(name: str, t_step: int):
    return f"{name}-{t_step:08d}.png"

def save(name: str, t_step: int):
    plt.savefig(os.path.join(OUTDIR, save_filename(name, t_step)))

def calculate_layout(n):
    total_ratio  = n * PLOT_DIMS[0] / PLOT_DIMS[1]
    
    rows    = int((total_ratio / TARGET_ASPECT_RATIO) ** 0.5)
    columns = int(rows * TARGET_ASPECT_RATIO)
    
    while rows * columns < n:
        if rows > columns:
            columns += 1
        else:
            rows += 1
    
    return rows, columns

rows, cols = calculate_layout(len(VARS))

for t_step, df in sorted(DATA.items(), key=lambda x: x[0]):
    print(f"Visualizing t_step = {t_step}:")

    fig, axes = plt.subplots(rows, cols, figsize=PLOT_DIMS)

    for i, name in enumerate(df.keys()[1:]):
        print(f" * Plotting {name}...")
        plot(name, t_step, axes[i % rows, i // rows])

    for i in range(len(DATA), rows * cols):
        print(f" * Skipping empty plot...")
        axes[i % rows, i // rows].axis('off')

    print(f" * Saving plot...")
    plt.suptitle(f"{CASENAME} | t_step = {t_step:08d}")
    plt.tight_layout()
    save("t", t_step)

    print(f" * Cleaning up...")
    plt.close('all')

print(f"\nDone: You can find your plots in the 'viz' directory ({os.path.relpath(OUTDIR, os.getcwd())}).")

os.system(f"ffmpeg -y -framerate {ARGS['fps']} -pattern_type glob -i '{os.path.join(OUTDIR, 't-*.png')}' -c:v libx264 -pix_fmt yuv420p '{CASENAME}.mp4'")
