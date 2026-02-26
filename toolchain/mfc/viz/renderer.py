"""
Image and video rendering for MFC visualization.

Produces PNG images (1D line plots, 2D colormaps) and MP4 videos
from assembled MFC data. Uses matplotlib with the Agg backend
for headless rendering.
"""

import math
import os
import tempfile

import numpy as np

import imageio.v2 as imageio

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # pylint: disable=wrong-import-position
from matplotlib.colors import LogNorm  # pylint: disable=wrong-import-position

matplotlib.rcParams.update({
    'mathtext.fontset': 'cm',
    'font.family': 'serif',
})


def render_1d(x_cc, data, varname, step, output, **opts):  # pylint: disable=too-many-arguments,too-many-positional-arguments
    """Render a 1D line plot and save as PNG."""
    fig, ax = plt.subplots(figsize=opts.get('figsize', (10, 6)))
    ax.plot(x_cc, data, linewidth=1.5)
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(varname)
    ax.set_title(f'{varname} (step {step})')
    ax.grid(True, alpha=0.3)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-3, 4), useMathText=True)

    vmin = opts.get('vmin')
    vmax = opts.get('vmax')
    if vmin is not None or vmax is not None:
        ax.set_ylim(vmin, vmax)

    fig.tight_layout()
    fig.savefig(output, dpi=opts.get('dpi', 150))
    plt.close(fig)


def render_1d_tiled(x_cc, variables, step, output, **opts):  # pylint: disable=too-many-locals
    """Render all 1D variables in a tiled subplot grid and save as PNG."""
    varnames = sorted(variables.keys())
    n = len(varnames)
    if n == 0:
        return
    if n == 1:
        render_1d(x_cc, variables[varnames[0]], varnames[0], step, output, **opts)
        return

    ncols = 2 if n <= 8 else 3
    nrows = math.ceil(n / ncols)
    fig_w = 5 * ncols
    fig_h = 2.8 * nrows
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=opts.get('figsize', (fig_w, fig_h)),
                             sharex=True, squeeze=False)

    for idx, vn in enumerate(varnames):
        row, col = divmod(idx, ncols)
        ax = axes[row][col]
        ax.plot(x_cc, variables[vn], linewidth=1.2)
        ax.set_ylabel(vn, fontsize=9)
        ax.tick_params(labelsize=8)
        ax.grid(True, alpha=0.3)

    # Hide unused subplots
    for idx in range(n, nrows * ncols):
        row, col = divmod(idx, ncols)
        axes[row][col].set_visible(False)

    # X-label only on bottom row
    for col in range(ncols):
        bottom_row = min(nrows - 1, (n - 1) // ncols) if col < (n % ncols or ncols) else nrows - 2
        axes[bottom_row][col].set_xlabel(r'$x$', fontsize=9)

    fig.suptitle(f'step {step}', fontsize=11, y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(output, dpi=opts.get('dpi', 150))
    plt.close(fig)


def _figsize_for_domain(x_cc, y_cc, base=10):
    """Compute figure size that matches the physical domain aspect ratio."""
    dx = float(x_cc[-1] - x_cc[0]) if len(x_cc) > 1 else 1.0
    dy = float(y_cc[-1] - y_cc[0]) if len(y_cc) > 1 else 1.0
    aspect = dy / dx if dx > 0 else 1.0
    # Clamp to avoid extremely tall/wide figures
    aspect = max(0.2, min(aspect, 5.0))
    # Extra width for colorbar
    fig_w = base + 1.5
    fig_h = max(base * aspect, 3.0)
    return (fig_w, fig_h)


def render_2d(x_cc, y_cc, data, varname, step, output, **opts):  # pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals
    """Render a 2D colormap via pcolormesh and save as PNG."""
    default_size = _figsize_for_domain(x_cc, y_cc)
    fig, ax = plt.subplots(figsize=opts.get('figsize', default_size))

    cmap = opts.get('cmap', 'viridis')
    vmin = opts.get('vmin')
    vmax = opts.get('vmax')
    log_scale = opts.get('log_scale', False)

    norm = None
    if log_scale:
        lo = vmin if vmin is not None else np.nanmin(data[data > 0]) if np.any(data > 0) else 1e-10
        hi = vmax if vmax is not None else np.nanmax(data)
        if not np.isfinite(hi) or hi <= 0:
            hi = 1.0
        if not np.isfinite(lo) or lo <= 0 or lo >= hi:
            lo = hi * 1e-10
        norm = LogNorm(vmin=lo, vmax=hi)
        vmin = None
        vmax = None

    # data shape is (nx, ny), pcolormesh expects (ny, nx) when using x_cc, y_cc
    pcm = ax.pcolormesh(x_cc, y_cc, data.T, cmap=cmap, vmin=vmin, vmax=vmax,
                        norm=norm, shading='auto')
    fig.colorbar(pcm, ax=ax, label=varname)
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')
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
    if slice_axis not in axis_map:
        raise ValueError(
            f"Invalid slice_axis '{slice_axis}'. Must be one of: 'x', 'y', 'z'."
        )
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
        xlabel, ylabel = r'$y$', r'$z$'
    elif axis_idx == 1:
        sliced = data_3d[:, idx, :]
        x_plot, y_plot = assembled.x_cc, assembled.z_cc
        xlabel, ylabel = r'$x$', r'$z$'
    else:
        sliced = data_3d[:, :, idx]
        x_plot, y_plot = assembled.x_cc, assembled.y_cc
        xlabel, ylabel = r'$x$', r'$y$'

    default_size = _figsize_for_domain(x_plot, y_plot)
    fig, ax = plt.subplots(figsize=opts.get('figsize', default_size))

    cmap = opts.get('cmap', 'viridis')
    vmin = opts.get('vmin')
    vmax = opts.get('vmax')
    log_scale = opts.get('log_scale', False)

    norm = None
    if log_scale:
        pos = sliced[sliced > 0]
        lo = vmin if vmin is not None else np.nanmin(pos) if pos.size > 0 else 1e-10
        hi = vmax if vmax is not None else np.nanmax(sliced)
        if not np.isfinite(hi) or hi <= 0:
            hi = 1.0
        if not np.isfinite(lo) or lo <= 0 or lo >= hi:
            lo = hi * 1e-10
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


def render_mp4(varname, steps, output, fps=10,  # pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals,too-many-statements,too-many-branches
               read_func=None, tiled=False, **opts):
    """
    Generate an MP4 video by iterating over timesteps.

    Args:
        varname: Variable name to plot (ignored when tiled=True).
        steps: List of timestep integers.
        output: Output MP4 file path.
        fps: Frames per second.
        read_func: Callable(step) -> AssembledData for loading each frame.
        tiled: If True, render all 1D variables in a tiled layout per frame.
        **opts: Rendering options (cmap, vmin, vmax, dpi, log_scale, figsize,
                slice_axis, slice_index, slice_value).

    Returns:
        True if the MP4 was successfully written, False on failure
        (e.g., missing imageio dependency or encoding error).
    """
    if read_func is None:
        raise ValueError("read_func must be provided for MP4 rendering")

    if not steps:
        raise ValueError("No timesteps provided for MP4 generation")

    opts = dict(opts)  # avoid mutating the caller's dict

    # Pre-compute vmin/vmax from first, middle, and last frames if not provided
    # (not needed for tiled mode — each subplot auto-scales independently)
    auto_vmin = opts.get('vmin')
    auto_vmax = opts.get('vmax')

    if not tiled and (auto_vmin is None or auto_vmax is None):
        sample_steps = [steps[0]]
        if len(steps) > 1:
            sample_steps.append(steps[-1])
        if len(steps) > 2:
            sample_steps.append(steps[len(steps) // 2])

        all_mins, all_maxs = [], []
        log_scale = opts.get('log_scale', False)
        for s in sample_steps:
            ad = read_func(s)
            d = ad.variables.get(varname)
            if d is not None:
                if log_scale:
                    pos = d[d > 0]
                    if pos.size > 0:
                        all_mins.append(np.nanmin(pos))
                        all_maxs.append(np.nanmax(pos))
                else:
                    all_mins.append(np.nanmin(d))
                    all_maxs.append(np.nanmax(d))

        if auto_vmin is None and all_mins:
            opts['vmin'] = min(all_mins)
        if auto_vmax is None and all_maxs:
            opts['vmax'] = max(all_maxs)

    # Write frames to a unique temp directory to avoid concurrent-run conflicts
    output_dir = os.path.dirname(os.path.abspath(output))
    os.makedirs(output_dir, exist_ok=True)
    viz_dir = tempfile.mkdtemp(dir=output_dir, prefix='_frames_')

    try:
        from tqdm import tqdm  # pylint: disable=import-outside-toplevel
        step_iter = tqdm(steps, desc='Rendering frames')
    except ImportError:
        step_iter = steps

    for i, step in enumerate(step_iter):
        assembled = read_func(step)
        frame_path = os.path.join(viz_dir, f'{i:06d}.png')

        if tiled and assembled.ndim == 1:
            render_1d_tiled(assembled.x_cc, assembled.variables,
                            step, frame_path, **opts)
        elif assembled.ndim == 1:
            var_data = assembled.variables.get(varname)
            if var_data is None:
                continue
            render_1d(assembled.x_cc, var_data,
                      varname, step, frame_path, **opts)
        elif assembled.ndim == 2:
            var_data = assembled.variables.get(varname)
            if var_data is None:
                continue
            render_2d(assembled.x_cc, assembled.y_cc,
                      var_data,
                      varname, step, frame_path, **opts)
        elif assembled.ndim == 3:
            var_data = assembled.variables.get(varname)
            if var_data is None:
                continue
            render_3d_slice(assembled, varname, step, frame_path, **opts)
        else:
            raise ValueError(
                f"Unsupported dimensionality ndim={assembled.ndim} for step {step}. "
                "Expected 1, 2, or 3."
            )

    # Combine frames into MP4 using imageio + imageio-ffmpeg (bundled ffmpeg)
    frame_files = sorted(f for f in os.listdir(viz_dir) if f.endswith('.png'))

    success = False
    try:
        imageio.mimwrite(
            output,
            [imageio.imread(os.path.join(viz_dir, fname)) for fname in frame_files],
            fps=fps, codec='libx264', pixelformat='yuv420p', macro_block_size=2,
            ffmpeg_log_level='error',
        )
        success = True
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"imageio MP4 write failed: {exc}")
    finally:
        # Always clean up temporary frame files
        for fname in frame_files:
            try:
                os.remove(os.path.join(viz_dir, fname))
            except OSError:
                pass
        try:
            os.rmdir(viz_dir)
        except OSError:
            pass
    return success
