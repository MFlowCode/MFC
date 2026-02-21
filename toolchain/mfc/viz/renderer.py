"""
Image and video rendering for MFC visualization.

Produces PNG images (1D line plots, 2D colormaps) and MP4 videos
from assembled MFC data. Uses matplotlib with the Agg backend
for headless rendering.
"""

import os
import subprocess

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # pylint: disable=wrong-import-position
from matplotlib.colors import LogNorm  # pylint: disable=wrong-import-position


def render_1d(x_cc, data, varname, step, output, **opts):  # pylint: disable=too-many-arguments,too-many-positional-arguments
    """Render a 1D line plot and save as PNG."""
    fig, ax = plt.subplots(figsize=opts.get('figsize', (10, 6)))
    ax.plot(x_cc, data, linewidth=1.5)
    ax.set_xlabel('x')
    ax.set_ylabel(varname)
    ax.set_title(f'{varname} (step {step})')

    vmin = opts.get('vmin')
    vmax = opts.get('vmax')
    if vmin is not None or vmax is not None:
        ax.set_ylim(vmin, vmax)

    fig.tight_layout()
    fig.savefig(output, dpi=opts.get('dpi', 150))
    plt.close(fig)


def render_2d(x_cc, y_cc, data, varname, step, output, **opts):  # pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals
    """Render a 2D colormap via pcolormesh and save as PNG."""
    fig, ax = plt.subplots(figsize=opts.get('figsize', (10, 8)))

    cmap = opts.get('cmap', 'viridis')
    vmin = opts.get('vmin')
    vmax = opts.get('vmax')
    log_scale = opts.get('log_scale', False)

    norm = None
    if log_scale:
        lo = vmin if vmin is not None else np.nanmin(data[data > 0]) if np.any(data > 0) else 1e-10
        hi = vmax if vmax is not None else np.nanmax(data)
        norm = LogNorm(vmin=lo, vmax=hi)
        vmin = None
        vmax = None

    # data shape is (nx, ny), pcolormesh expects (ny, nx) when using x_cc, y_cc
    pcm = ax.pcolormesh(x_cc, y_cc, data.T, cmap=cmap, vmin=vmin, vmax=vmax,
                        norm=norm, shading='auto')
    fig.colorbar(pcm, ax=ax, label=varname)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(f'{varname} (step {step})')
    ax.set_aspect('equal', adjustable='box')

    fig.tight_layout()
    fig.savefig(output, dpi=opts.get('dpi', 150))
    plt.close(fig)


def render_3d_slice(assembled, varname, step, output, slice_axis='z',  # pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals,too-many-statements,too-many-branches
                    slice_index=None, slice_value=None, **opts):
    """Extract a 2D slice from 3D data and render as a colormap."""
    data_3d = assembled.variables[varname]

    axis_map = {'x': 0, 'y': 1, 'z': 2}
    axis_idx = axis_map[slice_axis]

    coords = [assembled.x_cc, assembled.y_cc, assembled.z_cc]
    coord_along = coords[axis_idx]

    if slice_index is not None:
        idx = slice_index
    elif slice_value is not None:
        idx = int(np.argmin(np.abs(coord_along - slice_value)))
    else:
        idx = len(coord_along) // 2

    idx = max(0, min(idx, len(coord_along) - 1))

    if axis_idx == 0:
        sliced = data_3d[idx, :, :]
        x_plot, y_plot = assembled.y_cc, assembled.z_cc
        xlabel, ylabel = 'y', 'z'
    elif axis_idx == 1:
        sliced = data_3d[:, idx, :]
        x_plot, y_plot = assembled.x_cc, assembled.z_cc
        xlabel, ylabel = 'x', 'z'
    else:
        sliced = data_3d[:, :, idx]
        x_plot, y_plot = assembled.x_cc, assembled.y_cc
        xlabel, ylabel = 'x', 'y'

    fig, ax = plt.subplots(figsize=opts.get('figsize', (10, 8)))

    cmap = opts.get('cmap', 'viridis')
    vmin = opts.get('vmin')
    vmax = opts.get('vmax')
    log_scale = opts.get('log_scale', False)

    norm = None
    if log_scale:
        lo = vmin if vmin is not None else np.nanmin(sliced[sliced > 0]) if np.any(sliced > 0) else 1e-10
        hi = vmax if vmax is not None else np.nanmax(sliced)
        norm = LogNorm(vmin=lo, vmax=hi)
        vmin = None
        vmax = None

    # sliced shape depends on axis: need to transpose appropriately
    pcm = ax.pcolormesh(x_plot, y_plot, sliced.T, cmap=cmap, vmin=vmin,
                        vmax=vmax, norm=norm, shading='auto')
    fig.colorbar(pcm, ax=ax, label=varname)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    slice_coord = coord_along[idx]
    ax.set_title(f'{varname} (step {step}, {slice_axis}={slice_coord:.4g})')
    ax.set_aspect('equal', adjustable='box')

    fig.tight_layout()
    fig.savefig(output, dpi=opts.get('dpi', 150))
    plt.close(fig)


def render_mp4(case_dir, varname, steps, output, fps=10,  # pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals,too-many-statements
               read_func=None, **opts):
    """
    Generate an MP4 video by iterating over timesteps.

    Args:
        case_dir: Path to the case directory.
        varname: Variable name to plot.
        steps: List of timestep integers.
        output: Output MP4 file path.
        fps: Frames per second.
        read_func: Callable(step) -> AssembledData for loading each frame.
        **opts: Rendering options (cmap, vmin, vmax, dpi, log_scale, figsize,
                slice_axis, slice_index, slice_value).
    """
    if read_func is None:
        raise ValueError("read_func must be provided for MP4 rendering")

    if not steps:
        raise ValueError("No timesteps provided for MP4 generation")

    # Pre-compute vmin/vmax from first and last frames if not provided
    auto_vmin = opts.get('vmin')
    auto_vmax = opts.get('vmax')

    if auto_vmin is None or auto_vmax is None:
        sample_steps = [steps[0]]
        if len(steps) > 1:
            sample_steps.append(steps[-1])
        if len(steps) > 2:
            sample_steps.append(steps[len(steps) // 2])

        all_mins, all_maxs = [], []
        for s in sample_steps:
            ad = read_func(s)
            d = ad.variables.get(varname)
            if d is not None:
                all_mins.append(np.nanmin(d))
                all_maxs.append(np.nanmax(d))

        if auto_vmin is None and all_mins:
            opts['vmin'] = min(all_mins)
        if auto_vmax is None and all_maxs:
            opts['vmax'] = max(all_maxs)

    # Write frames as images to a temp directory
    viz_dir = os.path.join(case_dir, 'viz', '_frames')
    os.makedirs(viz_dir, exist_ok=True)

    try:
        from tqdm import tqdm  # pylint: disable=import-outside-toplevel
        step_iter = tqdm(steps, desc='Rendering frames')
    except ImportError:
        step_iter = steps

    for i, step in enumerate(step_iter):
        assembled = read_func(step)
        frame_path = os.path.join(viz_dir, f'{i:06d}.png')

        if assembled.ndim == 1:
            render_1d(assembled.x_cc, assembled.variables[varname],
                      varname, step, frame_path, **opts)
        elif assembled.ndim == 2:
            render_2d(assembled.x_cc, assembled.y_cc,
                      assembled.variables[varname],
                      varname, step, frame_path, **opts)
        elif assembled.ndim == 3:
            render_3d_slice(assembled, varname, step, frame_path, **opts)

    # Combine frames into MP4 using ffmpeg
    frame_pattern = os.path.join(viz_dir, '%06d.png')
    ffmpeg_cmd = [
        'ffmpeg', '-y',
        '-framerate', str(fps),
        '-i', frame_pattern,
        '-c:v', 'libx264',
        '-pix_fmt', 'yuv420p',
        '-vf', 'pad=ceil(iw/2)*2:ceil(ih/2)*2',
        output,
    ]

    try:
        subprocess.run(ffmpeg_cmd, check=True, capture_output=True)
    except FileNotFoundError:
        print(f"ffmpeg not found. Frames saved to {viz_dir}/")
        print(f"To create video manually: ffmpeg -framerate {fps} -i {frame_pattern} -c:v libx264 -pix_fmt yuv420p {output}")
        return
    except subprocess.CalledProcessError as e:
        print(f"ffmpeg failed: {e.stderr.decode()}")
        print(f"Frames saved to {viz_dir}/")
        return

    # Clean up frames
    for fname in os.listdir(viz_dir):
        os.remove(os.path.join(viz_dir, fname))
    os.rmdir(viz_dir)
