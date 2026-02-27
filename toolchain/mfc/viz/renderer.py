"""
Image and video rendering for MFC visualization.

Produces PNG images (1D line plots, 2D colormaps) and MP4 videos
from assembled MFC data. Uses matplotlib with the Agg backend
for headless rendering.
"""

import math
import os
import re
import tempfile

import numpy as np

import imageio

import matplotlib
try:
    matplotlib.use('Agg')
except Exception:  # pylint: disable=broad-except
    pass
import matplotlib.pyplot as plt  # pylint: disable=wrong-import-position
from matplotlib.colors import LogNorm  # pylint: disable=wrong-import-position

matplotlib.rcParams.update({
    'mathtext.fontset': 'cm',
    'font.family': 'serif',
})

# LaTeX-style labels for known MFC variable names
_LABEL_MAP = {
    'pres':     r'$p$',
    'rho':      r'$\rho$',
    'E':        r'$E$',
    'T':        r'$T$',
    'D':        r'$D$',
    'c':        r'$c$',
    'gamma':    r'$\gamma$',
    'pi_inf':   r'$\pi_\infty$',
    'pres_inf': r'$p_\infty$',
    'heat_ratio': r'$\gamma$',
    'schlieren':  r'$|\nabla \rho|$',
    'psi':        r'$\psi$',
    'n':          r'$n$',
    'qm':         r'$q_m$',
    'Bx': r'$B_x$', 'By': r'$B_y$', 'Bz': r'$B_z$',
    'voidFraction': r'void fraction',
    'liutex_mag':   r'$|\lambda|$',
    'damage_state': r'damage',
}

_INDEXED_PATTERNS = [
    (r'^vel(\d+)$',        lambda m: [r'$u$', r'$v$', r'$w$'][int(m.group(1)) - 1]
     if int(m.group(1)) <= 3 else rf'$v_{m.group(1)}$'),
    (r'^mom(\d+)$',        lambda m: rf'$\rho {["u", "v", "w"][int(m.group(1)) - 1]}$'
     if int(m.group(1)) <= 3 else rf'$m_{m.group(1)}$'),
    (r'^alpha(\d+)$',      lambda m: rf'$\alpha_{m.group(1)}$'),
    (r'^alpha_rho(\d+)$',  lambda m: rf'$\alpha_{m.group(1)}\rho_{m.group(1)}$'),
    (r'^alpha_rho_e(\d+)$', lambda m: rf'$\alpha_{m.group(1)}\rho_{m.group(1)}e_{m.group(1)}$'),
    (r'^omega(\d+)$',      lambda m: rf'$\omega_{m.group(1)}$'),
    (r'^tau(\d+)$',        lambda m: rf'$\tau_{m.group(1)}$'),
    (r'^xi(\d+)$',         lambda m: rf'$\xi_{m.group(1)}$'),
    (r'^flux(\d+)$',       lambda m: rf'$F_{m.group(1)}$'),
    (r'^liutex_axis(\d+)$', lambda m: rf'$\lambda_{m.group(1)}$'),
    (r'^rho(\d+)$',        lambda m: rf'$\rho_{m.group(1)}$'),
    (r'^Y_(.+)$',          lambda m: rf'$Y_{{\mathrm{{{m.group(1)}}}}}$'),
    (r'^nR(\d+)$',         lambda m: rf'$nR_{{{m.group(1)}}}$'),
    (r'^nV(\d+)$',         lambda m: rf'$nV_{{{m.group(1)}}}$'),
    (r'^nP(\d+)$',         lambda m: rf'$nP_{{{m.group(1)}}}$'),
    (r'^nM(\d+)$',         lambda m: rf'$nM_{{{m.group(1)}}}$'),
    (r'^color_function(\d+)$', lambda m: rf'color $f_{m.group(1)}$'),
]


def pretty_label(varname):
    """Map an MFC variable name to a LaTeX-style label for plots."""
    if varname in _LABEL_MAP:
        return _LABEL_MAP[varname]
    for pattern, formatter in _INDEXED_PATTERNS:
        m = re.match(pattern, varname)
        if m:
            return formatter(m)
    return varname


def render_1d(x_cc, data, varname, step, output, **opts):  # pylint: disable=too-many-arguments,too-many-positional-arguments
    """Render a 1D line plot and save as PNG."""
    fig, ax = plt.subplots(figsize=opts.get('figsize', (10, 6)))
    label = pretty_label(varname)
    ax.plot(x_cc, data, linewidth=1.5)
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(label)
    ax.set_title(f'{label} (step {step})')
    ax.grid(True, alpha=0.3)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-3, 4), useMathText=True)

    log_scale = opts.get('log_scale', False)
    if log_scale:
        ax.set_yscale('log')
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

    log_scale = opts.get('log_scale', False)
    for idx, vn in enumerate(varnames):
        row, col = divmod(idx, ncols)
        ax = axes[row][col]
        ax.plot(x_cc, variables[vn], linewidth=1.2)
        ax.set_ylabel(pretty_label(vn), fontsize=9)
        ax.tick_params(labelsize=8)
        ax.grid(True, alpha=0.3)
        if log_scale:
            ax.set_yscale('log')

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
    label = pretty_label(varname)
    fig.colorbar(pcm, ax=ax, label=label)
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')
    ax.set_title(f'{label} (step {step})')
    ax.set_aspect('equal', adjustable='box')

    fig.tight_layout()
    fig.savefig(output, dpi=opts.get('dpi', 150))
    plt.close(fig)


def render_2d_tiled(assembled, step, output, **opts):  # pylint: disable=too-many-locals
    """Render all 2D variables in a tiled subplot grid and save as PNG."""
    varnames = sorted(assembled.variables.keys())
    n = len(varnames)
    if n == 0:
        return
    if n == 1:
        render_2d(assembled.x_cc, assembled.y_cc,
                  assembled.variables[varnames[0]], varnames[0], step, output, **opts)
        return

    ncols = min(n, 3)
    nrows = math.ceil(n / ncols)
    cell_w, cell_h = _figsize_for_domain(assembled.x_cc, assembled.y_cc, base=4)
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=opts.get('figsize', (cell_w * ncols, cell_h * nrows)),
                             squeeze=False)

    cmap = opts.get('cmap', 'viridis')
    log_scale = opts.get('log_scale', False)
    for idx, vn in enumerate(varnames):
        row, col = divmod(idx, ncols)
        ax = axes[row][col]
        data = assembled.variables[vn]
        norm = None
        vmin = opts.get('vmin')
        vmax = opts.get('vmax')
        if log_scale and np.any(data > 0):
            lo = float(np.nanmin(data[data > 0]))
            hi = float(np.nanmax(data))
            if np.isfinite(lo) and np.isfinite(hi) and lo < hi:
                norm = LogNorm(vmin=lo, vmax=hi)
                vmin = None
                vmax = None
        pcm = ax.pcolormesh(assembled.x_cc, assembled.y_cc, data.T,
                            cmap=cmap, vmin=vmin, vmax=vmax,
                            norm=norm, shading='auto')
        label = pretty_label(vn)
        fig.colorbar(pcm, ax=ax, label=label)
        ax.set_title(label, fontsize=9)
        ax.set_aspect('equal', adjustable='box')
        ax.tick_params(labelsize=7)

    for idx in range(n, nrows * ncols):
        row, col = divmod(idx, ncols)
        axes[row][col].set_visible(False)

    fig.suptitle(f'step {step}', fontsize=11, y=1.01)
    fig.tight_layout()
    fig.savefig(output, dpi=opts.get('dpi', 150), bbox_inches='tight')
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
    label = pretty_label(varname)
    fig.colorbar(pcm, ax=ax, label=label)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    slice_coord = coord_along[idx]
    ax.set_title(f'{label} (step {step}, {slice_axis}={slice_coord:.4g})')
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

    def _cleanup():
        for fname in sorted(f for f in os.listdir(viz_dir) if f.endswith('.png')):
            try:
                os.remove(os.path.join(viz_dir, fname))
            except OSError:
                pass
        try:
            os.rmdir(viz_dir)
        except OSError:
            pass

    try:
        from tqdm import tqdm  # pylint: disable=import-outside-toplevel
        step_iter = tqdm(steps, desc='Rendering frames')
    except ImportError:
        step_iter = steps

    try:
        for i, step in enumerate(step_iter):
            assembled = read_func(step)
            frame_path = os.path.join(viz_dir, f'{i:06d}.png')

            if tiled and assembled.ndim == 1:
                render_1d_tiled(assembled.x_cc, assembled.variables,
                                step, frame_path, **opts)
            elif tiled and assembled.ndim == 2:
                render_2d_tiled(assembled, step, frame_path, **opts)
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
    except BaseException:
        _cleanup()
        raise

    # Combine frames into MP4 using imageio + imageio-ffmpeg (bundled ffmpeg)
    frame_files = sorted(f for f in os.listdir(viz_dir) if f.endswith('.png'))

    success = False
    try:
        if not frame_files:
            return False

        def _to_rgb(arr):
            """Normalise an image array to uint8 RGB (3-channel).

            imageio may return RGBA (4-ch) or even grayscale depending on the
            PNG source.  libx264/yuv420p requires consistent 3-channel input.
            """
            if arr.ndim == 2:                        # grayscale → RGB
                arr = np.stack([arr, arr, arr], axis=-1)
            elif arr.shape[2] == 4:                  # RGBA → RGB (drop alpha)
                arr = arr[:, :, :3]
            return arr.astype(np.uint8)

        # First pass: find the maximum frame dimensions across all frames.
        # Round up to the nearest even pixel: yuv420p requires even width and height.
        max_h, max_w = 0, 0
        for fname in frame_files:
            arr = imageio.imread(os.path.join(viz_dir, fname))
            max_h = max(max_h, arr.shape[0])
            max_w = max(max_w, arr.shape[1])
        ref_h = max_h + max_h % 2
        ref_w = max_w + max_w % 2

        def _uniform_frame(arr):
            """Convert to RGB and pad with white to (ref_h, ref_w)."""
            arr = _to_rgb(arr)
            h, w = arr.shape[:2]
            if h == ref_h and w == ref_w:
                return arr
            out = np.full((ref_h, ref_w, 3), 255, dtype=np.uint8)
            out[:h, :w] = arr
            return out

        # Second pass: encode.  macro_block_size=1 disables imageio's own resize
        # since we already ensured even dimensions above.
        with imageio.get_writer(
            output, fps=fps, codec='libx264', pixelformat='yuv420p',
            macro_block_size=1, ffmpeg_log_level='error',
        ) as writer:
            for fname in frame_files:
                writer.append_data(_uniform_frame(
                    imageio.imread(os.path.join(viz_dir, fname))
                ))
        success = True
    except Exception as exc:  # pylint: disable=broad-except
        import warnings  # pylint: disable=import-outside-toplevel
        warnings.warn(f"imageio MP4 write failed: {exc}", stacklevel=2)
    finally:
        _cleanup()
    return success
