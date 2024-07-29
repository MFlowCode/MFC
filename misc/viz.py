#!/usr/bin/env python3

import os
import glob
import math
import argparse
import multiprocessing as mp
import seaborn as sns
import pandas  as pd
import matplotlib.pyplot as plt

PLOTS = {
    'Velocity': [('$u_x$',                 'prim.2', '#00FFFF')],
    'Pressure': [('p',                     'prim.3', '#FFFF00')],
    'Density':  [('$\\rho$',               'prim.1', '#FF0000')],
    'Momentum density': [('$\\rho u_x$',                 'cons.2', '#00FF00')],
    'Alpha': [('$\\alpha$',                 'prim.4', '#00FF00')],
}


BG_COLOR = '#1a1a1a'
TX_COLOR = '#FFFFFF'

sns.set_style('dark', {
    'axes.facecolor':    '#121212',
    'axes.edgecolor':    BG_COLOR,
    'axes.labelcolor':   TX_COLOR,
    'text.color':        TX_COLOR,
    'xtick.color':       TX_COLOR,
    'ytick.color':       TX_COLOR,
    'grid.color':        BG_COLOR,
    'figure.facecolor':  BG_COLOR,
    'figure.edgecolor':  BG_COLOR,
    'savefig.facecolor': BG_COLOR,
    'savefig.edgecolor': BG_COLOR,
})

parser = argparse.ArgumentParser(description="Visualize the results of a case")
parser.add_argument('case_dir',                       type=str, help="Path to the case directory")
parser.add_argument('-f', '--fps',        default=24, type=int, help="Frames per second for the video")
parser.add_argument('-n', '--num-frames',             type=int, help="Number of frames to process")
parser.add_argument('-s', '--start-frame',            type=int, help="Starting frame to process")

ARGS        = vars(parser.parse_args())
DIRPATH     = os.path.abspath(ARGS['case_dir'])
OUTDIR      = os.path.join(DIRPATH, 'viz')
CASENAME    = os.path.basename(DIRPATH)
DESCRIPTION = "0D React - Chemistry = T, Advection = T, Diffusion = F, Reactions = T"

TARGET_ASPECT_RATIO = 16 / 9

PLOT_PADDING = 0.04
PLOT_DIMS    = (16, 9)

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

rows, cols = calculate_layout(len(PLOTS.keys()))
rows, cols = 3, 3


DATA = {}
MFC_VARS = set()
VARLIMITS = {}


def plot(name: str, t_step: int, ax=None):
    ax.set_xlabel('$x$', labelpad=10)
    ax.set_ylabel(name,  labelpad=10)

    for p in PLOTS[name]:
        axisname, varname, color = p

        ys = None
        if callable(varname):
            ys = [ varname(x, t_step) for x in DATA[t_step]['x'] ]
        else:
            ys = varname

        sns.lineplot(
            data=DATA[t_step], x='x', y=ys,
            color=color, linewidth=1.5, alpha=1,
            ax=ax, dashes=False, label=axisname)

    try:
        ax.set_xlim(*VARLIMITS['x'])
    except:
        print(f" * Warning: Could not set x-axis limits for {name}.")

    ymin = min([VARLIMITS[p[1]][0] for p in PLOTS[name]])
    ymax = max([VARLIMITS[p[1]][1] for p in PLOTS[name]])

    # 10% less and more
    delta = ymax - ymin
    ymin  = ymin - 0.1 * delta
    ymax  = ymax + 0.1 * delta

    if ymin != ymax:
        try:
            ax.set_ylim(ymin, ymax)
        except:
            print(f" * Warning: Could not set y-axis limits for {name}.")

def save_filename(t_step: int):
    return f"t-{t_step:08d}.png"

def save(t_step: int):
    plt.savefig(os.path.join(OUTDIR, save_filename(t_step)), dpi=300)
    plt.close('all')
    plt.cla()


print(f"Importing data from {os.path.relpath(DIRPATH, os.getcwd())}...")

def load_var_data(var):
    print(f"* {var}:")

    dfs = {}
    if callable(var):
        print(f"  * Using custom function for {var}.")
        for t_step in DATA.keys():
            xs = DATA[t_step]['x']
            ys = [ var(x, t_step) for x in xs ]
            dfs[t_step] = pd.DataFrame({
                'x': xs,
                var.__name__: ys
            })
    else:
        print(f"  * Reading data files...")
        for f in glob.glob(os.path.join(DIRPATH, 'D', f'{var}.*.*.dat')):
            proc, t_step = int(f.split('.')[-3]), int(f.split('.')[-2])

            if t_step not in dfs:
                dfs[t_step] = []

            dfs[t_step].append(pd.read_csv(f, sep=r'\s+', header=None, names=['x', var]))

        if len(dfs) == 0:
            print(f" * Error: No data found for {var}.")
            exit(1)

        print(f"  * Concatenating processors...")
        for t_step in dfs.keys():
            dfs[t_step] = pd.concat(dfs[t_step])

    print(f"  * Merging with existing data & steps...")
    for t_step, df in dfs.items():  
        if t_step not in DATA:
            DATA[t_step] = df
        else:
            DATA[t_step] = pd.merge(DATA[t_step], df)


for name, plots in PLOTS.items():
    for p in plots:
        MFC_VARS.add(p[1])

bFound = False
for var in MFC_VARS:
    if not callable(var):
        load_var_data(var)            
        bFound = True
        break

if not bFound:
    print(f" * Error: No non-analytical variables found.")
    exit(1)

for var in set(MFC_VARS) - set(var):
    load_var_data(var)

# a bunch of renaming :)
MFC_VARS = { _.__name__ if callable(_) else _ for _ in MFC_VARS }
for i, name in enumerate(PLOTS.keys()):
    for j, p in enumerate(PLOTS[name]):
        if callable(p[1]):
            PLOTS[name][j] = (p[0], p[1].__name__, *p[2:])


print(f"* Calculating min/max...")
for var in list(MFC_VARS) + ['x']:
    print(f" * {var}...")

    minmax = (math.inf, -math.inf)

    for step, df in DATA.items():
        minmax = (
            min(minmax[0], df[var].min()),
            max(minmax[1], df[var].max())
        )

    VARLIMITS[var] = (
        minmax[0] - PLOT_PADDING * (minmax[1] - minmax[0]),
        minmax[1] + PLOT_PADDING * (minmax[1] - minmax[0])
    )


def _worker(t_steps: list):
    for i, t_step in enumerate(t_steps):
        print(f"\nProcessing frame {t_step}...")
        _, axes = plt.subplots(rows, cols, figsize=(13,13))

        for i, name in enumerate(PLOTS.keys()):
            plot(name, t_step, axes[i % rows, i // rows])

        for i in range(len(PLOTS.keys()), rows * cols):
            axes[i % rows, i // rows].axis('off')

        plt.suptitle(DESCRIPTION, fontsize=24)
        plt.tight_layout(h_pad=2, w_pad=2)
        save(t_step)

        plt.close('all')


if __name__ == '__main__':
    if not os.path.isdir(DIRPATH):
        print(f"Error: {DIRPATH} is not a case directory.")
        exit(1)

    if not os.path.isdir(OUTDIR):
        os.mkdir(OUTDIR)

    with mp.Pool(mp.cpu_count()) as pool:
        steps       = sorted(DATA.keys())
        step_chunks = [steps[i::mp.cpu_count()] for i in range(mp.cpu_count())]

        pool.map(func=_worker, iterable=step_chunks)

    print(f"\nDone: You can find your plots in the 'viz' directory ({os.path.relpath(OUTDIR, os.getcwd())}).")

    command = f"ffmpeg -y -framerate {ARGS['fps']} -pattern_type glob -i '{os.path.join(OUTDIR, 't-*.png')}' -c:v libx264 -pix_fmt yuv420p '{DIRPATH}/{CASENAME}.mp4'"
    print(command)
    os.system(command)
