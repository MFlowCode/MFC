"""
Interactive Dash-based visualization for MFC post-processed data.

Launched via ``./mfc.sh viz <dir> --var <var> --step all --interactive``.
Opens a dark-themed web UI in your browser (or via SSH tunnel) with live
controls for slice position, isosurface thresholds, volume opacity,
colormap, log scale, vmin/vmax, and timestep playback.
"""
# pylint: disable=use-dict-literal,too-many-lines

import atexit
import base64
import concurrent.futures
import logging
import math
import os
import sys
import threading
import time
from typing import List, Callable, Optional

import numpy as np
import plotly.graph_objects as go
import pyvista as pv
from dash import Dash, Patch, dcc, html, Input, Output, State, callback_context, no_update
from skimage.measure import marching_cubes as _marching_cubes  # type: ignore[import]  # pylint: disable=no-name-in-module
from skimage.measure import find_contours as _find_contours  # type: ignore[import]  # pylint: disable=no-name-in-module

from mfc.printer import cons
from . import _step_cache
from ._step_cache import prefetch_one as _prefetch_one

pv.OFF_SCREEN = True

# Suppress VTK's internal error/warning output (e.g. EGL failures) so it
# doesn't flood stderr.  We handle failures at the Python level instead.
try:
    import vtk as _vtk  # type: ignore[import]
    _vtk.vtkOutputWindow.GetInstance().SetGlobalWarningDisplay(0)  # pylint: disable=no-member
except Exception:  # pylint: disable=broad-except
    pass

logger = logging.getLogger(__name__)

# PyTurboJPEG wraps libjpeg-turbo via ctypes. TurboJPEG() opens the native
# shared library at instantiation time; on HPC clusters where libjpeg-turbo
# is not in LD_LIBRARY_PATH it raises OSError even if the pip package is
# installed. Fall back to Pillow in that case.
try:
    from turbojpeg import TurboJPEG as _TurboJPEG  # type: ignore[import]
    _tj = _TurboJPEG()
except (ImportError, OSError):
    _tj = None

# ---------------------------------------------------------------------------
# Fast PNG generation via 256-entry colormap LUT
# ---------------------------------------------------------------------------

# Both caches grow to at most one entry per named colormap (~10 entries max).
# Concurrent writes are benign: both threads produce identical data for the
# same colormap name, so no lock is needed.
_lut_cache: dict = {}    # cmap_name → (256, 3) uint8 LUT
_cscale_cache: dict = {} # cmap_name → plotly colorscale list

# Downsampled-3D cache: (step, var, max_total) → (raw_ds, x_ds, y_ds, z_ds)
_ds3_cache: dict = {}
_DS3_CACHE_MAX = 10
_ds3_lock = threading.Lock()

# JPEG pre-encode cache: (step, var, cmap, vmin_in, vmax_in, log_bool) → src str
_jpeg_cache: dict = {}
_JPEG_CACHE_MAX = 30
_jpeg_lock = threading.Lock()
_jpeg_pool: Optional[concurrent.futures.ThreadPoolExecutor] = None
_jpeg_pool_lock = threading.Lock()

# 3D mesh prefetch cache: key → (vx, vy, vz, fi, fj, fk, intens) or volume data
_mesh3_cache: dict = {}
_MESH3_CACHE_MAX = 20
_mesh3_lock = threading.Lock()
_mesh3_in_flight: set = set()
_mesh3_pool: Optional[concurrent.futures.ThreadPoolExecutor] = None
_mesh3_pool_lock = threading.Lock()


def _get_jpeg_pool() -> concurrent.futures.ThreadPoolExecutor:
    """Return the JPEG prefetch pool, creating it lazily on first use."""
    global _jpeg_pool  # pylint: disable=global-statement
    with _jpeg_pool_lock:
        if _jpeg_pool is None:
            _jpeg_pool = concurrent.futures.ThreadPoolExecutor(
                max_workers=1, thread_name_prefix='mfc_jpeg')
            atexit.register(_jpeg_pool.shutdown, wait=False)
        return _jpeg_pool


def _get_mesh3_pool() -> concurrent.futures.ThreadPoolExecutor:
    """Return the 3D mesh prefetch pool, creating it lazily on first use."""
    global _mesh3_pool  # pylint: disable=global-statement
    with _mesh3_pool_lock:
        if _mesh3_pool is None:
            _mesh3_pool = concurrent.futures.ThreadPoolExecutor(
                max_workers=2, thread_name_prefix='mfc_mesh3')
            atexit.register(_mesh3_pool.shutdown, wait=False)
        return _mesh3_pool


def _get_lut(cmap_name: str) -> np.ndarray:
    """Return a (256, 3) uint8 LUT for the named matplotlib colormap."""
    if cmap_name not in _lut_cache:
        try:
            import matplotlib as mpl  # pylint: disable=import-outside-toplevel
            cm = mpl.colormaps.get_cmap(cmap_name)
        except (ImportError, KeyError):
            import matplotlib.cm as mcm  # pylint: disable=import-outside-toplevel
            cm = mcm.get_cmap(cmap_name)
        t = np.linspace(0.0, 1.0, 256)
        _lut_cache[cmap_name] = (cm(t)[:, :3] * 255 + 0.5).astype(np.uint8)
    return _lut_cache[cmap_name]


def _lut_to_plotly_colorscale(cmap_name: str) -> list:
    """Convert a matplotlib LUT to a Plotly colorscale list.

    Samples 32 evenly-spaced entries from the LUT — enough for a smooth
    colorbar while keeping the JSON payload small.  Result is cached per
    colormap name since it is called on every render callback.
    """
    if cmap_name in _cscale_cache:
        return _cscale_cache[cmap_name]
    lut = _get_lut(cmap_name)
    indices = np.linspace(0, 255, 32, dtype=int)
    result = [
        [float(i) / 31.0, f'rgb({lut[idx,0]},{lut[idx,1]},{lut[idx,2]})']
        for i, idx in enumerate(indices)
    ]
    _cscale_cache[cmap_name] = result
    return result


def _encode_jpeg(rgb: np.ndarray) -> bytes:
    """Encode (h, w, 3) uint8 RGB → JPEG bytes.

    Uses libjpeg-turbo when available; falls back to Pillow otherwise.
    """
    if _tj is not None:
        return _tj.encode(rgb, quality=90)
    from PIL import Image as _PIL  # pylint: disable=import-outside-toplevel
    import io as _io              # pylint: disable=import-outside-toplevel
    buf = _io.BytesIO()
    _PIL.fromarray(rgb, 'RGB').save(buf, format='jpeg', quality=90, optimize=False)
    return buf.getvalue()


def _make_png_source(arr_yx: np.ndarray, cmap_name: str,
                     vmin: float, vmax: float) -> str:
    """Encode a (ny, nx) float array as a colorized base64 JPEG data URI.

    Uses a 256-entry LUT for fast colormap application and libjpeg-turbo
    (when available) for JPEG encoding.

    arr_yx: shape (ny, nx), row 0 = smallest y (physical bottom).
    No vertical flip is applied — go.Image with y0=y_min, dy>0 renders
    row 0 at the bottom of the y-axis, matching physics convention.
    """
    lut = _get_lut(cmap_name)
    scale = max(float(vmax - vmin), 1e-30)
    normed = np.clip(
        (arr_yx - vmin) / scale * 255.0 + 0.5, 0, 255
    ).astype(np.uint8)
    rgb = lut[normed]                      # (ny, nx, 3) uint8
    b64 = base64.b64encode(_encode_jpeg(rgb)).decode()
    return f"data:image/jpeg;base64,{b64}"


def _compute_isomesh(raw_ds: np.ndarray, x_ds: np.ndarray, y_ds: np.ndarray,  # pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals
                     z_ds: np.ndarray, log_fn, ilo: float, ihi: float,
                     iso_n: int):
    """Server-side marching cubes for *iso_n* levels between *ilo* and *ihi*.

    Returns flat arrays (vx, vy, vz, fi, fj, fk, intensity) ready for
    go.Mesh3d.  All faces from all levels are concatenated; intensity holds
    the isosurface level value per vertex so the colorscale matches.

    Much faster than go.Isosurface which runs marching cubes in JavaScript.
    """
    vol = log_fn(raw_ds).astype(np.float64)
    # Replace NaN/inf (e.g. from log of zero/negative values) with a fill
    # value below any isosurface level so those cells produce no triangles.
    finite_mask = np.isfinite(vol)
    if not np.all(finite_mask):
        fill = float(np.nanmin(vol[finite_mask])) - 1.0 if np.any(finite_mask) else -1.0
        vol = np.where(finite_mask, vol, fill)

    dx = float(x_ds[-1] - x_ds[0]) / max(len(x_ds) - 1, 1) if len(x_ds) > 1 else 1.0
    dy = float(y_ds[-1] - y_ds[0]) / max(len(y_ds) - 1, 1) if len(y_ds) > 1 else 1.0
    dz = float(z_ds[-1] - z_ds[0]) / max(len(z_ds) - 1, 1) if len(z_ds) > 1 else 1.0
    spacing = (dx, dy, dz)

    levels = np.linspace(ilo, ihi, max(int(iso_n), 1))
    xs, ys, zs, ii, jj, kk, intens = [], [], [], [], [], [], []
    offset = 0

    for level in levels:
        try:
            verts, faces, _, _ = _marching_cubes(
                vol, level=float(level), spacing=spacing,
                allow_degenerate=False,
            )
        except (ValueError, RuntimeError):
            continue
        if len(faces) == 0:
            continue
        # Shift vertices to physical origin
        verts[:, 0] += float(x_ds[0])
        verts[:, 1] += float(y_ds[0])
        verts[:, 2] += float(z_ds[0])
        xs.append(verts[:, 0]); ys.append(verts[:, 1]); zs.append(verts[:, 2])
        ii.append(faces[:, 0] + offset)
        jj.append(faces[:, 1] + offset)
        kk.append(faces[:, 2] + offset)
        intens.append(np.full(len(verts), float(level), dtype=np.float32))
        offset += len(verts)

    if not xs:
        # No surface found — tiny invisible dummy mesh
        dummy = np.zeros(3, dtype=np.float32)
        return dummy, dummy, dummy, np.array([0]), np.array([1]), np.array([2]), dummy

    # float32 halves the JSON payload vs float64 (~50% smaller .tolist())
    return (np.concatenate(xs).astype(np.float32),
            np.concatenate(ys).astype(np.float32),
            np.concatenate(zs).astype(np.float32),
            np.concatenate(ii).astype(np.int32),
            np.concatenate(jj).astype(np.int32),
            np.concatenate(kk).astype(np.int32),
            np.concatenate(intens).astype(np.float32))


def _downsample_3d(raw: np.ndarray, x_cc: np.ndarray, y_cc: np.ndarray,
                   z_cc: np.ndarray, max_total: int = 150_000):
    """Stride a (nx, ny, nz) array to stay within a total cell budget.

    Uses a **uniform** stride s = ceil((nx*ny*nz / max_total)^(1/3)) so that
    all axes are sampled equally.  For an anisotropic grid like 901×201×201,
    this gives stride=7 → 129×29×29 = 108K cells instead of the per-axis
    strategy which would give stride=18 in x (only 50 pts), causing jagged
    isosurfaces along the long axis.
    """
    nx, ny, nz = raw.shape
    total = nx * ny * nz
    if total <= max_total:
        return raw, x_cc, y_cc, z_cc
    s = max(1, math.ceil((total / max_total) ** (1.0 / 3.0)))
    return raw[::s, ::s, ::s], x_cc[::s], y_cc[::s], z_cc[::s]


def _get_ds3(step, var, raw, x_cc, y_cc, z_cc, max_total):  # pylint: disable=too-many-arguments,too-many-positional-arguments
    """Downsampled 3D array with bounded LRU caching.

    Avoids re-striding the same large array on every iso threshold / volume
    opacity slider move.  Key is (step, var, max_total); value is the tuple
    (raw_ds, x_ds, y_ds, z_ds) returned by _downsample_3d.
    """
    key = (step, var, max_total)
    with _ds3_lock:
        if key in _ds3_cache:
            return _ds3_cache[key]
    result = _downsample_3d(raw, x_cc, y_cc, z_cc, max_total)
    with _ds3_lock:
        if key not in _ds3_cache:
            if len(_ds3_cache) >= _DS3_CACHE_MAX:
                # FIFO eviction: next(iter(dict)) yields the oldest entry by
                # insertion order, guaranteed in CPython 3.7+ and the language
                # spec from Python 3.7.
                _ds3_cache.pop(next(iter(_ds3_cache)))
            _ds3_cache[key] = result
    return result


def _prefetch_jpeg(step, var, get_ad_fn, cmap, vmin_in, vmax_in, log_bool,  # pylint: disable=too-many-arguments,too-many-positional-arguments
                   max_nx=1200, max_ny=600):
    """Pre-encode JPEG for *step/var* in a background thread.

    Called immediately after data prefetch so the JPEG is ready by the time
    the user navigates to the next step.  No-ops if already cached.

    get_ad_fn: callable(step) → AssembledData (wraps _step_cache.load).
    """
    key = (step, var, cmap, vmin_in, vmax_in, log_bool)
    with _jpeg_lock:
        if key in _jpeg_cache:
            return

    def _bg():  # pylint: disable=too-many-locals
        try:
            ad = get_ad_fn(step)
            if var not in ad.variables:
                return
            raw = ad.variables[var]
            # Subsample for range estimation (same logic as main callback)
            if raw.size > 200_000:
                _sr = max(1, math.ceil((raw.size / 200_000) ** (1.0 / raw.ndim)))
                _slc = tuple(slice(None, None, _sr) for _ in range(raw.ndim))
                _rr = raw[_slc]
            else:
                _rr = raw
            if vmin_in is not None:
                vmin = float(vmin_in)
            else:
                _safe = _rr[_rr > 0] if log_bool and np.any(_rr > 0) else _rr
                vmin = float(np.nanmin(_safe))
            vmax = float(vmax_in) if vmax_in is not None else float(np.nanmax(_rr))
            if vmax <= vmin:
                vmax = vmin + 1e-10
            if log_bool:
                cmin = float(np.log10(max(vmin, 1e-300)))
                cmax = float(np.log10(max(vmax, 1e-300)))
                display = np.where(raw > 0, np.log10(np.maximum(raw, 1e-300)), np.nan)
            else:
                cmin, cmax = vmin, vmax
                display = raw.astype(float)
            nx2, ny2 = display.shape
            _sx = max(1, nx2 // max_nx)
            _sy = max(1, ny2 // max_ny)
            if _sx > 1 or _sy > 1:
                display = display[::_sx, ::_sy]
            src = _make_png_source(display.T, cmap, cmin, cmax)
            with _jpeg_lock:
                if key not in _jpeg_cache:
                    if len(_jpeg_cache) >= _JPEG_CACHE_MAX:
                        _jpeg_cache.pop(next(iter(_jpeg_cache)))
                    _jpeg_cache[key] = src
        except Exception:  # pylint: disable=broad-except
            logger.debug("JPEG prefetch failed for step %s var %s", step, var, exc_info=True)

    _get_jpeg_pool().submit(_bg)


def _prefetch_3d_mesh(step, var, get_ad_fn, mode, log_bool,  # pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals,too-many-statements
                      vmin_in, vmax_in,
                      iso_min_frac, iso_max_frac, iso_n,
                      vol_min_frac, vol_max_frac):
    """Pre-compute 3D isomesh or volume data for *step* in a background thread.

    Caches the result so the next playback frame can skip marching cubes.
    No-ops if already cached or in-flight.
    """
    key = (step, var, mode, log_bool, vmin_in, vmax_in,
           iso_min_frac, iso_max_frac, iso_n,
           vol_min_frac, vol_max_frac)
    with _mesh3_lock:
        if key in _mesh3_cache or key in _mesh3_in_flight:
            return
        _mesh3_in_flight.add(key)

    def _bg():  # pylint: disable=too-many-locals
        try:
            ad = get_ad_fn(step)
            if var not in ad.variables:
                return
            raw = ad.variables[var]
            # Compute range
            if raw.size > 200_000:
                _sr = max(1, math.ceil((raw.size / 200_000) ** (1.0 / raw.ndim)))
                _slc = tuple(slice(None, None, _sr) for _ in range(raw.ndim))
                _rr = raw[_slc]
            else:
                _rr = raw
            if vmin_in is not None:
                vmin = float(vmin_in)
            else:
                _safe = _rr[_rr > 0] if log_bool and np.any(_rr > 0) else _rr
                vmin = float(np.nanmin(_safe))
            vmax = float(vmax_in) if vmax_in is not None else float(np.nanmax(_rr))
            if vmax <= vmin:
                vmax = vmin + 1e-10
            if log_bool:
                def _tf(arr):
                    return np.where(arr > 0, np.log10(np.maximum(arr, 1e-300)), np.nan)
                cmin = float(np.log10(max(vmin, 1e-300)))
                cmax = float(np.log10(max(vmax, 1e-300)))
            else:
                def _tf(arr): return arr
                cmin, cmax = vmin, vmax
            rng = cmax - cmin if cmax > cmin else 1.0

            if mode == 'isosurface':
                raw_ds, x_ds, y_ds, z_ds = _get_ds3(
                    step, var, raw, ad.x_cc, ad.y_cc, ad.z_cc, 500_000)
                ilo = cmin + rng * float(iso_min_frac)
                ihi = cmin + rng * max(float(iso_max_frac), float(iso_min_frac) + 0.01)
                result = _compute_isomesh(
                    raw_ds, x_ds, y_ds, z_ds, _tf, ilo, ihi, int(iso_n))
            else:  # volume
                raw_ds, _, _, _ = _get_ds3(
                    step, var, raw, ad.x_cc, ad.y_cc, ad.z_cc, 150_000)
                vf = _tf(raw_ds).ravel().astype(np.float32)
                vlo = cmin + rng * float(vol_min_frac)
                vhi = cmin + rng * max(float(vol_max_frac), float(vol_min_frac) + 0.01)
                result = (vf, vlo, vhi, cmin, cmax)

            with _mesh3_lock:
                if key not in _mesh3_cache:
                    if len(_mesh3_cache) >= _MESH3_CACHE_MAX:
                        _mesh3_cache.pop(next(iter(_mesh3_cache)))
                    _mesh3_cache[key] = result
        except Exception:  # pylint: disable=broad-except
            logger.debug("3D mesh prefetch failed for step %s var %s", step, var, exc_info=True)
        finally:
            with _mesh3_lock:
                _mesh3_in_flight.discard(key)

    _get_mesh3_pool().submit(_bg)


def _get_cached_3d_mesh(step, var, mode, log_bool, vmin_in, vmax_in,  # pylint: disable=too-many-arguments,too-many-positional-arguments
                        iso_min_frac, iso_max_frac, iso_n,
                        vol_min_frac, vol_max_frac):
    """Return cached 3D mesh result or None if not yet computed."""
    key = (step, var, mode, log_bool, vmin_in, vmax_in,
           iso_min_frac, iso_max_frac, iso_n,
           vol_min_frac, vol_max_frac)
    with _mesh3_lock:
        return _mesh3_cache.get(key)


# ---------------------------------------------------------------------------
# PyVista server-side 3D renderer
# ---------------------------------------------------------------------------

# VTK's Cocoa backend (macOS) requires ALL OpenGL operations on the main
# thread.  Dash callbacks run in worker threads, so PyVista server-side
# rendering is only available on Linux/other platforms where VTK uses
# EGL backends.  All PyVista calls are funnelled through a single
# dedicated thread (_pv_thread) so the EGL context is never accessed
# from multiple threads (EGL contexts are thread-bound).
_PV_AVAILABLE = sys.platform != 'darwin'

# Persistent plotter — owned exclusively by _pv_thread.
_pv_plotter: Optional[pv.Plotter] = None
_pv_step_key: Optional[tuple] = None  # (step, var, mode, ...) currently loaded
# Module-level kill switch: once set, _pv_render raises immediately without
# submitting to the thread pool.  This prevents in-flight callbacks from
# queueing behind a hung render thread.
_pv_disabled = threading.Event()

# Dedicated render thread: a single-worker pool that owns the EGL context.
_pv_thread: Optional[concurrent.futures.ThreadPoolExecutor] = None
_pv_thread_lock = threading.Lock()


def _get_pv_thread() -> concurrent.futures.ThreadPoolExecutor:
    """Return the dedicated PyVista render thread, creating it on first use."""
    global _pv_thread  # pylint: disable=global-statement
    with _pv_thread_lock:
        if _pv_thread is None:
            _pv_thread = concurrent.futures.ThreadPoolExecutor(
                max_workers=1, thread_name_prefix='mfc_pv')
            atexit.register(_pv_thread.shutdown, wait=False)
        return _pv_thread


def _pv_ensure_plotter() -> pv.Plotter:
    """Return the persistent plotter, creating it lazily on first use.

    MUST be called on the dedicated _pv_thread — EGL/GLX contexts are
    thread-bound and cannot be shared across threads.

    Temporarily removes DISPLAY so VTK does not try to connect to an
    X server (common on HPC nodes with SSH X-forwarding).
    """
    global _pv_plotter  # pylint: disable=global-statement
    if _pv_plotter is not None:
        return _pv_plotter
    saved = os.environ.pop('DISPLAY', None)
    try:
        pl = pv.Plotter(off_screen=True, window_size=(800, 600),
                        lighting='three lights')
        pl.set_background('#181825')
        pl.render()  # warm up — ensures the render window is fully initialized
        # Verify the render actually produced pixels (catches broken EGL)
        img = pl.screenshot(return_img=True)
        if img is None or img.size == 0:
            raise RuntimeError('PyVista screenshot returned empty image')
        _pv_plotter = pl
    finally:
        if saved is not None:
            os.environ['DISPLAY'] = saved
    return _pv_plotter


def _pv_probe_job():
    """Run on the dedicated thread: init plotter + render a tiny test mesh."""
    pl = _pv_ensure_plotter()
    # Render a small cube to verify the full pipeline (not just init)
    grid = pv.ImageData()
    grid.dimensions = (3, 3, 3)
    grid.point_data['t'] = np.arange(27, dtype=np.float64)
    pl.clear()
    pl.set_background('#181825')
    pl.add_mesh(grid, scalars='t', cmap='viridis', show_scalar_bar=False)
    pl.render()
    img = pl.screenshot(return_img=True)
    pl.clear()
    if img is None or img.size == 0:
        raise RuntimeError('PyVista probe render returned empty image')
    return True


def _pv_probe() -> bool:
    """Test whether PyVista off-screen rendering works on this system.

    Submits a full render job (init + mesh + screenshot) to the dedicated
    render thread.  Returns True if a valid image is produced, False otherwise.
    """
    if not _PV_AVAILABLE:
        return False
    try:
        pool = _get_pv_thread()
        fut = pool.submit(_pv_probe_job)
        fut.result(timeout=15)
        return True
    except Exception as exc:  # pylint: disable=broad-except
        logger.info('PyVista probe failed: %s', exc)
        return False


def _pv_render_on_thread(raw, x_cc, y_cc, z_cc, mode, cmap, log_fn,  # pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals,too-many-branches,too-many-statements
                         cmin, cmax,
                         iso_min_frac, iso_max_frac, iso_n,
                         vol_min_frac, vol_max_frac,
                         camera, step, varname,
                         overlay_raw, overlay_x, overlay_y,
                         overlay_z, overlay_mode,
                         overlay_iso_min, overlay_iso_max,
                         overlay_nlevels, overlay_color,
                         overlay_opacity,
                         overlay_vol_min, overlay_vol_max,
                         overlay_vol_opacity, overlay_vol_nsurf):  # pylint: disable=unused-argument
    """Inner render — runs exclusively on the dedicated _pv_thread."""
    global _pv_step_key  # pylint: disable=global-statement

    rng = cmax - cmin if cmax > cmin else 1.0
    mesh_key = (step, varname, mode, id(raw),
                iso_min_frac, iso_max_frac, iso_n,
                vol_min_frac, vol_max_frac,
                overlay_raw is not None)

    pl = _pv_ensure_plotter()
    need_rebuild = _pv_step_key != mesh_key

    if need_rebuild:
        pl.clear()
        pl.set_background('#181825')

        ds_raw, ds_x, ds_y, ds_z = _downsample_3d(
            raw, x_cc, y_cc, z_cc, 500_000)

        nx, ny, nz = ds_raw.shape
        grid = pv.ImageData()
        grid.dimensions = (nx, ny, nz)
        grid.origin = (float(ds_x[0]), float(ds_y[0]), float(ds_z[0]))
        if nx > 1:
            grid.spacing = (
                float(ds_x[-1] - ds_x[0]) / max(nx - 1, 1),
                float(ds_y[-1] - ds_y[0]) / max(ny - 1, 1),
                float(ds_z[-1] - ds_z[0]) / max(nz - 1, 1),
            )
        scalar = log_fn(ds_raw).astype(np.float64)
        scalar = np.where(np.isfinite(scalar), scalar, np.nan)
        grid.point_data['scalar'] = scalar.flatten(order='F')

        if mode == 'isosurface':
            ilo = cmin + rng * iso_min_frac
            ihi = cmin + rng * max(iso_max_frac, iso_min_frac + 0.01)
            levels = np.linspace(ilo, ihi, max(int(iso_n), 1)).tolist()
            try:
                iso = grid.contour(levels, scalars='scalar')
                if iso.n_points > 0:
                    pl.add_mesh(iso, scalars='scalar', cmap=cmap,
                                clim=[ilo, ihi],
                                smooth_shading=True, opacity=1.0,
                                show_scalar_bar=True,
                                scalar_bar_args=dict(
                                    title=varname, color='#cdd6f4',
                                    label_font_size=10, title_font_size=12))
            except Exception:  # pylint: disable=broad-except
                pass

        elif mode == 'volume':
            vlo = cmin + rng * vol_min_frac
            vhi = cmin + rng * max(vol_max_frac, vol_min_frac + 0.01)
            grid_vol = grid.threshold([vlo, vhi], scalars='scalar')
            if grid_vol.n_points > 0:
                pl.add_volume(grid_vol, scalars='scalar', cmap=cmap,
                              clim=[cmin, cmax],
                              opacity='sigmoid',
                              show_scalar_bar=True,
                              scalar_bar_args=dict(
                                  title=varname, color='#cdd6f4',
                                  label_font_size=10, title_font_size=12))

        if overlay_raw is not None:
            ov_ds, ov_x, ov_y, ov_z = _downsample_3d(
                overlay_raw, overlay_x, overlay_y, overlay_z, 500_000)
            o_nx, o_ny, o_nz = ov_ds.shape
            ov_grid = pv.ImageData()
            ov_grid.dimensions = (o_nx, o_ny, o_nz)
            ov_grid.origin = (float(ov_x[0]), float(ov_y[0]), float(ov_z[0]))
            if o_nx > 1:
                ov_grid.spacing = (
                    float(ov_x[-1] - ov_x[0]) / max(o_nx - 1, 1),
                    float(ov_y[-1] - ov_y[0]) / max(o_ny - 1, 1),
                    float(ov_z[-1] - ov_z[0]) / max(o_nz - 1, 1),
                )
            ov_scalar = log_fn(ov_ds).astype(np.float64)
            ov_scalar = np.where(np.isfinite(ov_scalar), ov_scalar, np.nan)
            ov_grid.point_data['scalar'] = ov_scalar.flatten(order='F')
            ov_vmin = float(np.nanmin(ov_scalar))
            ov_vmax = float(np.nanmax(ov_scalar))
            ov_rng = ov_vmax - ov_vmin if ov_vmax > ov_vmin else 1.0

            if overlay_mode == 'isosurface':
                ov_lo = ov_vmin + ov_rng * overlay_iso_min
                ov_hi = ov_vmin + ov_rng * max(overlay_iso_max,
                                                overlay_iso_min + 0.01)
                ov_lvls = np.linspace(ov_lo, ov_hi,
                                      max(int(overlay_nlevels), 1)).tolist()
                try:
                    ov_iso = ov_grid.contour(ov_lvls, scalars='scalar')
                    if ov_iso.n_points > 0:
                        pl.add_mesh(ov_iso, color=overlay_color,
                                    opacity=float(overlay_opacity),
                                    smooth_shading=True,
                                    show_scalar_bar=False)
                except Exception:  # pylint: disable=broad-except
                    pass
            else:
                ov_vlo = ov_vmin + ov_rng * overlay_vol_min
                ov_vhi = ov_vmin + ov_rng * max(overlay_vol_max,
                                                 overlay_vol_min + 0.01)
                ov_vol = ov_grid.threshold([ov_vlo, ov_vhi],
                                           scalars='scalar')
                if ov_vol.n_points > 0:
                    pl.add_volume(ov_vol, scalars='scalar', cmap=cmap,
                                  clim=[ov_vmin, ov_vmax],
                                  opacity='sigmoid',
                                  show_scalar_bar=False)

        _pv_step_key = mesh_key

    if camera and 'position' in camera:
        pl.camera.position = camera['position']
        pl.camera.focal_point = camera['focal_point']
        pl.camera.up = camera['view_up']
    elif need_rebuild:
        pl.camera_position = 'iso'

    pl.render()
    img = pl.screenshot(return_img=True)

    b64 = base64.b64encode(_encode_jpeg(img)).decode()
    return f"data:image/jpeg;base64,{b64}"


def _pv_render(raw, x_cc, y_cc, z_cc, mode, cmap, log_fn,  # pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals
               cmin, cmax,
               iso_min_frac, iso_max_frac, iso_n,
               vol_opacity, vol_nsurf, vol_min_frac, vol_max_frac,  # pylint: disable=unused-argument
               camera, step, varname,
               overlay_raw=None, overlay_x=None, overlay_y=None,
               overlay_z=None, overlay_mode='isosurface',
               overlay_iso_min=0.2, overlay_iso_max=0.8,
               overlay_nlevels=3, overlay_color='white',
               overlay_opacity=0.6,
               overlay_vol_min=0.0, overlay_vol_max=1.0,
               overlay_vol_opacity=0.1, overlay_vol_nsurf=15,  # pylint: disable=unused-argument
               width=800, height=600):  # pylint: disable=unused-argument
    """Render a 3D scene with PyVista and return a JPEG base64 data URI.

    Submits the render job to the dedicated _pv_thread so that the EGL/GLX
    context is always accessed from the same thread.  Raises on failure so
    callers can fall back to the Plotly path.
    """
    if _pv_disabled.is_set():
        raise RuntimeError('PyVista disabled')
    pool = _get_pv_thread()
    fut = pool.submit(
        _pv_render_on_thread,
        raw, x_cc, y_cc, z_cc, mode, cmap, log_fn,
        cmin, cmax,
        iso_min_frac, iso_max_frac, iso_n,
        vol_min_frac, vol_max_frac,
        camera, step, varname,
        overlay_raw, overlay_x, overlay_y,
        overlay_z, overlay_mode,
        overlay_iso_min, overlay_iso_max,
        overlay_nlevels, overlay_color,
        overlay_opacity,
        overlay_vol_min, overlay_vol_max,
        overlay_vol_opacity, overlay_vol_nsurf,
    )
    return fut.result(timeout=10)


# ---------------------------------------------------------------------------
# Colormaps available in the picker
# ---------------------------------------------------------------------------
_CMAPS = [
    "viridis", "plasma", "inferno", "magma", "cividis",
    "turbo", "jet", "rainbow", "nipy_spectral",
    "RdBu", "RdYlBu", "RdYlGn", "coolwarm", "bwr", "seismic", "Spectral",
    "hot", "afmhot", "gist_heat", "copper",
    "bone", "gray", "spring", "summer", "autumn", "winter", "cool", "pink",
    "Blues", "Greens", "Oranges", "Reds", "Purples", "Greys",
    "twilight", "twilight_shifted", "hsv",
    "tab10", "tab20", "terrain", "ocean", "gist_earth",
    "gnuplot", "gnuplot2", "CMRmap", "cubehelix", "Wistia",
]

# ---------------------------------------------------------------------------
# Catppuccin Mocha palette
# ---------------------------------------------------------------------------
_BG     = '#181825'
_SURF   = '#1e1e2e'
_OVER   = '#313244'
_BORDER   = '#45475a'
_TEXT   = '#cdd6f4'
_SUB    = '#a6adc8'
_MUTED  = '#6c7086'
_ACCENT = '#cba6f7'
_GREEN  = '#a6e3a1'
_RED    = '#f38ba8'
_BLUE   = '#89b4fa'
_TEAL   = '#94e2d5'
_YELLOW = '#f9e2af'

# ---------------------------------------------------------------------------
# Server-side data cache  {step -> AssembledData}  (bounded to avoid OOM)
# ---------------------------------------------------------------------------
_load = _step_cache.load
_CACHE_MAX = _step_cache.CACHE_MAX


# ---------------------------------------------------------------------------
# Layout helpers
# ---------------------------------------------------------------------------

def _section(title, *children):
    return html.Div([
        html.Div(title, style={
            'fontSize': '10px', 'fontWeight': 'bold',
            'textTransform': 'uppercase', 'letterSpacing': '0.08em',
            'color': _MUTED, 'borderBottom': f'1px solid {_OVER}',
            'paddingBottom': '4px', 'marginTop': '16px', 'marginBottom': '6px',
        }),
        *children,
    ])


def _lbl(text):
    return html.Div(text, style={
        'fontSize': '11px', 'color': _SUB,
        'marginBottom': '2px', 'marginTop': '6px',
    })


def _slider(sid, lo, hi, step, val, marks=None):  # pylint: disable=too-many-arguments,too-many-positional-arguments
    return dcc.Slider(
        id=sid, min=lo, max=hi, step=step, value=val,
        marks=marks or {}, updatemode='mouseup',
    )


def _btn(bid, label, color=_TEXT):
    return html.Button(label, id=bid, n_clicks=0, style={
        'flex': '1', 'padding': '5px 8px', 'fontSize': '12px',
        'backgroundColor': _OVER, 'color': color,
        'border': f'1px solid {_BORDER}', 'borderRadius': '4px',
        'cursor': 'pointer', 'fontFamily': 'monospace',
    })


def _num(sid, placeholder='auto'):
    return dcc.Input(
        id=sid, type='number', placeholder=placeholder, debounce=True,
        style={
            'width': '100%', 'backgroundColor': _OVER, 'color': _TEXT,
            'border': f'1px solid {_BORDER}', 'borderRadius': '4px',
            'padding': '4px 6px', 'fontSize': '12px', 'fontFamily': 'monospace',
            'boxSizing': 'border-box', 'colorScheme': 'dark',
        },
    )


# ---------------------------------------------------------------------------
# 3D figure builder
# ---------------------------------------------------------------------------

def _make_cbar(title_text: str, cmin: float, cmax: float, n: int = 6) -> dict:
    """Colorbar dict with Python-formatted tick labels.

    Uses Python's f'{v:.2e}' formatter which zero-pads exponents to two digits
    (e.g. 1.23e+04) unlike Plotly's d3-format which gives variable-width
    exponents (e.g. 1.23e+4).  Custom tickvals/ticktext override d3.
    """
    tickvals = np.linspace(cmin, cmax, n).tolist()
    ticktext = [f'{v:.2e}' for v in tickvals]
    return dict(
        title=dict(text=title_text, font=dict(color=_TEXT)),
        tickfont=dict(color=_TEXT),
        tickmode='array',
        tickvals=tickvals,
        ticktext=ticktext,
        len=0.75, lenmode='fraction',
        yanchor='middle', y=0.5,
        thickness=15,
    )


def _slice_3d(raw, log_fn, x_cc, y_cc, z_cc, slice_axis, slice_pos,  # pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals
              max_pts=(600, 300)):
    """Extract and downsample a 2-D slice from a 3-D array.

    Returns (sliced_ds, coord1_ds, coord2_ds, actual, const_axis_value)
    where sliced_ds is the downsampled surfacecolor array and coord*_ds are
    the downsampled coordinate vectors for the two varying axes.
    """
    axis_coords = {'x': x_cc, 'y': y_cc, 'z': z_cc}
    coords = axis_coords[slice_axis]
    coord_val = coords[0] + (coords[-1] - coords[0]) * slice_pos
    idx = int(np.clip(np.argmin(np.abs(coords - coord_val)), 0, len(coords) - 1))
    actual = float(coords[idx])

    if slice_axis == 'x':
        sliced = log_fn(raw[idx, :, :])   # (ny, nz)
        c1, c2 = y_cc, z_cc
    elif slice_axis == 'y':
        sliced = log_fn(raw[:, idx, :])   # (nx, nz)
        c1, c2 = x_cc, z_cc
    else:
        sliced = log_fn(raw[:, :, idx])   # (nx, ny)
        c1, c2 = x_cc, y_cc

    s1 = max(1, sliced.shape[0] // max_pts[0])
    s2 = max(1, sliced.shape[1] // max_pts[1])
    return sliced[::s1, ::s2], c1[::s1], c2[::s2], actual


def _build_3d(ad, raw, varname, step, mode, cmap,  # pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals,too-many-branches,too-many-statements,unused-argument
              log_fn, cmin, cmax, cbar_title,
              slice_axis, slice_pos,
              iso_min_frac, iso_max_frac, iso_n, _iso_caps,
              vol_opacity, vol_nsurf, vol_min_frac, vol_max_frac,
              iso_solid_color=None, iso_opacity=1.0,
              cached_mesh=None,
              max_total_3d: int = 150_000):
    """Return (trace, title) for a 3D assembled dataset.

    Downsamples iso/volume data via a uniform stride targeting *max_total_3d*
    total cells, and slice data to ≤600×300, keeping JSON payloads small.
    Uses the matplotlib LUT-derived colorscale so colors match exactly.
    """
    cscale = _lut_to_plotly_colorscale(cmap)
    cbar = _make_cbar(cbar_title, cmin, cmax)
    rng = cmax - cmin if cmax > cmin else 1.0

    if mode == 'slice':
        sliced, c1, c2, actual = _slice_3d(
            raw, log_fn, ad.x_cc, ad.y_cc, ad.z_cc,
            slice_axis, slice_pos,
        )
        # Use float32 for all arrays — halves JSON payload vs float64.
        # For z-slices (the default), also pass 1D x/y coordinate vectors:
        # Plotly Surface accepts 1D x (n cols) and 1D y (m rows) with z shape
        # (m, n), saving the two full meshgrid arrays (~164× less coord data).
        # Convention: z[j,i] → (x[i], y[j]), so surfacecolor must be sliced.T.
        # x and y slices keep 2D coords but use float32.
        n1, n2 = len(c1), len(c2)
        if slice_axis == 'x':
            YY, ZZ = np.meshgrid(c1.astype(np.float32), c2.astype(np.float32),
                                 indexing='ij')
            trace = go.Surface(
                x=np.full((n1, n2), actual, dtype=np.float32),
                y=YY, z=ZZ,
                surfacecolor=sliced.astype(np.float32),
                cmin=cmin, cmax=cmax,
                colorscale=cscale, colorbar=cbar, showscale=True,
            )
        elif slice_axis == 'y':
            XX, ZZ = np.meshgrid(c1.astype(np.float32), c2.astype(np.float32),
                                 indexing='ij')
            trace = go.Surface(
                x=XX, y=np.full((n1, n2), actual, dtype=np.float32), z=ZZ,
                surfacecolor=sliced.astype(np.float32),
                cmin=cmin, cmax=cmax,
                colorscale=cscale, colorbar=cbar, showscale=True,
            )
        else:
            # z-slice: x and y both vary — use 1D vectors.
            # Plotly rows=y, cols=x, so z and surfacecolor are shape (n2, n1).
            trace = go.Surface(
                x=c1.astype(np.float32),
                y=c2.astype(np.float32),
                z=np.full((n2, n1), actual, dtype=np.float32),
                surfacecolor=sliced.T.astype(np.float32),
                cmin=cmin, cmax=cmax,
                colorscale=cscale, colorbar=cbar, showscale=True,
            )
        title = f'{varname}  ·  {slice_axis} = {actual:.4g}  ·  step {step}'

    elif mode == 'isosurface':
        # Use a higher-resolution grid for server-side marching cubes.
        # skimage.marching_cubes runs in compiled C — 10–100× faster than
        # Plotly's go.Isosurface which runs marching cubes in JavaScript.
        # Targeting 500K cells gives stride=5 on a 901×201×201 grid (181×41×41),
        # which is 5× finer than the 150K volume budget.
        ilo = cmin + rng * iso_min_frac
        ihi = cmin + rng * max(iso_max_frac, iso_min_frac + 0.01)
        if cached_mesh is not None:
            vx, vy, vz, fi, fj, fk, intens = cached_mesh
        else:
            raw_ds, x_ds, y_ds, z_ds = _get_ds3(step, varname, raw, ad.x_cc, ad.y_cc, ad.z_cc, 500_000)
            vx, vy, vz, fi, fj, fk, intens = _compute_isomesh(
                raw_ds, x_ds, y_ds, z_ds, log_fn, ilo, ihi, iso_n,
            )
        _iso_op = float(iso_opacity) if iso_opacity is not None else 1.0
        if iso_solid_color:
            # Use a uniform-intensity + single-color colorscale so that the
            # mesh rendering path is the same as the variable-colored one
            # (Plotly's Mesh3d `color` prop is unreliable for overriding the
            # default colorscale behavior).
            _solid_cs = [[0, iso_solid_color], [1, iso_solid_color]]
            trace = go.Mesh3d(
                x=vx, y=vy, z=vz, i=fi, j=fj, k=fk,
                intensity=np.zeros(len(vx), dtype=np.float32),
                intensitymode='vertex',
                colorscale=_solid_cs, cmin=0, cmax=1,
                showscale=False, opacity=_iso_op,
                lighting=dict(ambient=0.7, diffuse=0.9, specular=0.3,
                              roughness=0.5, fresnel=0.2),
                lightposition=dict(x=1000, y=500, z=500),
                flatshading=False,
            )
        else:
            trace = go.Mesh3d(
                x=vx, y=vy, z=vz, i=fi, j=fj, k=fk,
                intensity=intens, intensitymode='vertex',
                colorscale=cscale, cmin=ilo, cmax=ihi,
                colorbar=_make_cbar(cbar_title, ilo, ihi), showscale=True,
                opacity=_iso_op,
                lighting=dict(ambient=0.7, diffuse=0.9, specular=0.3,
                              roughness=0.5, fresnel=0.2),
                lightposition=dict(x=1000, y=500, z=500),
                flatshading=False,
            )
        title = f'{varname}  ·  {int(iso_n)} isosurfaces  ·  step {step}'

    else:                                                # volume
        raw_ds, x_ds, y_ds, z_ds = _get_ds3(step, varname, raw, ad.x_cc, ad.y_cc, ad.z_cc, max_total_3d)
        X3, Y3, Z3 = np.meshgrid(x_ds, y_ds, z_ds, indexing='ij')
        vf = log_fn(raw_ds.ravel()).astype(np.float32)
        vlo = cmin + rng * vol_min_frac
        vhi = cmin + rng * max(vol_max_frac, vol_min_frac + 0.01)
        trace = go.Volume(
            x=X3.ravel().astype(np.float32),
            y=Y3.ravel().astype(np.float32),
            z=Z3.ravel().astype(np.float32),
            value=vf,
            isomin=vlo, isomax=vhi,
            opacity=float(vol_opacity), surface_count=int(vol_nsurf),
            colorscale=cscale, cmin=cmin, cmax=cmax, colorbar=cbar,
        )
        title = f'{varname}  ·  volume  ·  step {step}'

    return trace, title


# ---------------------------------------------------------------------------
# Contour overlay helpers
# ---------------------------------------------------------------------------

def _interp_indices(indices, coords):
    """Map fractional array indices to physical coordinates via linear interp."""
    cl = np.clip(indices, 0, len(coords) - 1)
    fl = np.floor(cl).astype(int)
    frac = cl - fl
    ce = np.minimum(fl + 1, len(coords) - 1)
    return coords[fl] * (1 - frac) + coords[ce] * frac


def _compute_contour_traces(data_2d, x_cc, y_cc, nlevels, color, lw):  # pylint: disable=too-many-arguments,too-many-positional-arguments
    """Compute isocontour lines on a 2D array via skimage.find_contours.

    Returns a list of go.Scatter traces (one per contour segment) that can be
    added on top of a heatmap figure.  Coordinates are mapped from array
    indices to physical space using *x_cc* and *y_cc*.
    """
    vmin_c, vmax_c = float(np.nanmin(data_2d)), float(np.nanmax(data_2d))
    if vmax_c <= vmin_c:
        return []
    levels = np.linspace(vmin_c, vmax_c, nlevels + 2)[1:-1]  # exclude endpoints
    traces = []
    for level in levels:
        contours = _find_contours(data_2d, level)
        for contour in contours:
            # contour is (N, 2) in row/col index space
            px = _interp_indices(contour[:, 0], x_cc)
            py = _interp_indices(contour[:, 1], y_cc)
            traces.append(go.Scatter(
                x=px, y=py, mode='lines',
                line=dict(color=color, width=lw),
                showlegend=False, hoverinfo='skip',
            ))
    return traces


def _compute_contour_traces_3d(data_3d, x_cc, y_cc, z_cc,  # pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals
                                slice_axis, slice_pos,
                                nlevels, color, lw):
    """Compute isocontour lines on a 3D slice, returned as Scatter3d traces."""
    axis_coords = {'x': x_cc, 'y': y_cc, 'z': z_cc}
    coords = axis_coords[slice_axis]
    coord_val = coords[0] + (coords[-1] - coords[0]) * slice_pos
    idx = int(np.clip(np.argmin(np.abs(coords - coord_val)), 0, len(coords) - 1))
    actual = float(coords[idx])

    if slice_axis == 'x':
        sliced = data_3d[idx, :, :]   # (ny, nz)
        c1, c2 = y_cc, z_cc
    elif slice_axis == 'y':
        sliced = data_3d[:, idx, :]   # (nx, nz)
        c1, c2 = x_cc, z_cc
    else:
        sliced = data_3d[:, :, idx]   # (nx, ny)
        c1, c2 = x_cc, y_cc

    vmin_c, vmax_c = float(np.nanmin(sliced)), float(np.nanmax(sliced))
    if vmax_c <= vmin_c:
        return []
    levels = np.linspace(vmin_c, vmax_c, nlevels + 2)[1:-1]
    traces = []
    for level in levels:
        contours = _find_contours(sliced, level)
        for contour in contours:
            p1 = _interp_indices(contour[:, 0], c1)
            p2 = _interp_indices(contour[:, 1], c2)
            n = len(p1)
            const = np.full(n, actual, dtype=np.float32)
            if slice_axis == 'x':
                sx, sy, sz = const, p1, p2
            elif slice_axis == 'y':
                sx, sy, sz = p1, const, p2
            else:
                sx, sy, sz = p1, p2, const
            traces.append(go.Scatter3d(
                x=sx, y=sy, z=sz, mode='lines',
                line=dict(color=color, width=lw * 2),  # thicker in 3D
                showlegend=False, hoverinfo='skip',
            ))
    return traces


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run_interactive(  # pylint: disable=too-many-locals,too-many-statements,too-many-arguments,too-many-positional-arguments
        varname: str,
        steps: List[int],
        read_func: Callable,
        port: int = 8050,
        host: str = '127.0.0.1',
        bubble_func: Optional[Callable] = None,
        read_one_var_func: Optional[Callable] = None,
):
    """Launch the interactive Dash visualization server.

    Args:
        bubble_func: Optional callable ``(step: int) -> np.ndarray | None``.
            Must return a float64 array of shape ``(N, 4)`` with columns
            ``[x, y, z, radius]`` in simulation-normalized units, or ``None``
            when no bubble data is available for that step.
    """
    app = Dash(
        __name__,
        title=f'MFC viz · {varname}',
        suppress_callback_exceptions=True,
    )

    # Override Dash's internal component styles for dark theme.
    # Dash components (Dropdown, Input, Slider) render internal DOM that
    # ignores inline style props.  We inject CSS via app.index_string.
    # Build CSS using %-formatting to avoid f-string brace conflicts.
    _V = {'bg': _OVER, 'tx': _TEXT, 'bd': _BORDER, 'ac': _ACCENT}
    _dark_css = """
* { color-scheme: dark; }
/* Dropdowns — target by known IDs + universal child selectors */
#var-sel *, #cmap-sel *, #overlay-var-sel *, #overlay-color-sel *, #iso-solid-color *, #overlay-mode-sel *,
#var-sel > div, #cmap-sel > div, #overlay-var-sel > div, #overlay-color-sel > div, #iso-solid-color > div, #overlay-mode-sel > div {
    background-color: %(bg)s !important;
    color: %(tx)s !important;
    border-color: %(bd)s !important;
}
#var-sel input, #cmap-sel input, #overlay-var-sel input, #overlay-color-sel input, #iso-solid-color input, #overlay-mode-sel input {
    background-color: %(bg)s !important;
    color: %(tx)s !important;
}
/* Inputs */
input, input[type=number], input[type=text] {
    background-color: %(bg)s !important;
    color: %(tx)s !important;
    border: 1px solid %(bd)s !important;
    border-radius: 4px;
}
/* Number input spinner buttons (browser chrome) */
input[type=number]::-webkit-inner-spin-button,
input[type=number]::-webkit-outer-spin-button {
    -webkit-appearance: none;
    margin: 0;
}
input[type=number] { -moz-appearance: textfield; }
.rc-slider-mark-text { color: %(tx)s !important; }
.rc-slider-rail { background-color: %(bd)s !important; }
.rc-slider-track { background-color: %(ac)s !important; }
.rc-slider-handle {
    border-color: %(ac)s !important;
    background-color: %(ac)s !important;
}
.rc-slider-dot { background-color: %(bd)s !important; border-color: %(bd)s !important; }
/* Radio button labels — browser default is dark text; force light */
input[type=radio] + span, label { color: %(tx)s !important; }
/* Dash 4 uses CSS custom properties for all component styling.
   Override them globally so the dropdown panel (rendered in a Radix portal)
   also picks up the dark theme. */
:root, .dash-dropdown, .dash-dropdown-content,
.dash-dropdown-wrapper, .dash-dropdown-search-container {
    --Dash-Fill-Inverse-Strong: %(bg)s;
    --Dash-Text-Strong:         %(tx)s;
    --Dash-Text-Weak:           #a6adc8;
    --Dash-Text-Disabled:       #6c7086;
    --Dash-Stroke-Strong:       %(bd)s;
    --Dash-Fill-Disabled:       %(bd)s;
    --Dash-Fill-Interactive-Strong: %(ac)s;
    --Dash-Fill-Interactive-Weak:   #313244;
    --Dash-Shading-Strong: rgba(0,0,0,0.6);
    --Dash-Shading-Weak:   rgba(0,0,0,0.3);
}
/* Option hover / active / selected states */
.dash-dropdown-option:hover,
.dash-dropdown-option.active {
    background-color: %(bd)s !important;
}
.dash-dropdown-option.selected {
    background-color: %(ac)s !important;
    color: #11111b !important;
}
""" % _V
    app.index_string = (
        '<!DOCTYPE html>\n<html>\n<head>\n'
        '{%metas%}\n<title>{%title%}</title>\n{%favicon%}\n{%css%}\n'
        '<style>\n' + _dark_css + '\n</style>\n'
        '</head>\n<body>\n'
        '{%app_entry%}\n<footer>\n{%config%}\n{%scripts%}\n{%renderer%}\n'
        '</footer>\n</body>\n</html>'
    )

    # Load first step to know dimensionality and available variables
    init = _load(steps[0], read_func)
    ndim = init.ndim
    all_varnames = sorted(init.variables.keys())
    if varname not in all_varnames:
        varname = all_varnames[0] if all_varnames else varname

    var_opts  = [{'label': v, 'value': v} for v in all_varnames]
    cmap_opts = [{'label': c, 'value': c} for c in _CMAPS]

    if ndim == 3:
        mode_opts = [
            {'label': '  Slice',      'value': 'slice'},
            {'label': '  Isosurface', 'value': 'isosurface'},
            {'label': '  Volume',     'value': 'volume'},
        ]
    elif ndim == 2:
        mode_opts = [{'label': '  Heatmap', 'value': 'heatmap'}]
    else:
        mode_opts = [{'label': '  Line', 'value': 'line'}]
    default_mode = mode_opts[0]['value']

    # ------------------------------------------------------------------
    # Sidebar layout
    # ------------------------------------------------------------------
    sidebar = html.Div([

        # Header
        html.Div('MFC viz', style={
            'fontSize': '16px', 'fontWeight': 'bold', 'color': _ACCENT,
        }),
        html.Div(
            f'{ndim}D  ·  {len(steps)} step{"s" if len(steps) != 1 else ""}',
            style={'fontSize': '11px', 'color': _MUTED},
        ),

        # ── Variable ──────────────────────────────────────────────────
        _section('Variable',
            dcc.Dropdown(
                id='var-sel', options=var_opts, value=varname, clearable=False,
                style={'fontSize': '12px', 'backgroundColor': _OVER,
                       'border': f'1px solid {_BORDER}'},
            ),
        ),

        # ── Timestep ──────────────────────────────────────────────────
        _section('Timestep',
            dcc.Slider(
                id='step-sl',
                min=0, max=len(steps) - 1, step=1, value=0,
                marks={i: {'label': str(s), 'style': {'fontSize': '9px', 'color': _MUTED}}
                       for i, s in enumerate(steps)}
                      if len(steps) <= 10
                      else {0: {'label': str(steps[0]), 'style': {'fontSize': '9px', 'color': _MUTED}},
                            len(steps) - 1: {'label': str(steps[-1]), 'style': {'fontSize': '9px', 'color': _MUTED}}},
                updatemode='mouseup',
            ),
            html.Div(id='step-label', style={
                'fontSize': '11px', 'color': _YELLOW, 'textAlign': 'center',
                'marginTop': '2px',
            }),
            # Hidden store that holds the actual step value (not the index)
            dcc.Store(id='step-sel', data=steps[0]),
            html.Div([
                _btn('play-btn', '▶  Play', _GREEN),
                html.Div(style={'width': '6px'}),
                _btn('stop-btn', '■  Stop', _RED),
            ], style={'display': 'flex', 'marginTop': '6px'}),
            _lbl('Playback speed (fps)'),
            _slider('fps-sl', 0.5, 10, 0.5, 2,
                    marks={0.5: '0.5', 2: '2', 5: '5', 10: '10'}),
            dcc.Checklist(
                id='loop-chk',
                options=[{'label': '  Loop', 'value': 'loop'}], value=['loop'],
                style={'fontSize': '12px', 'color': _SUB, 'marginTop': '4px'},
            ),
        ),

        # ── Viz mode ──────────────────────────────────────────────────
        _section('Viz Mode',
            dcc.RadioItems(
                id='mode-sel', options=mode_opts, value=default_mode,
                style={'fontSize': '12px', 'color': _TEXT},
                inputStyle={'marginRight': '6px', 'cursor': 'pointer'},
                labelStyle={'display': 'block', 'marginBottom': '5px',
                            'cursor': 'pointer', 'color': _TEXT},
            ),
        ),

        # ── Slice ─────────────────────────────────────────────────────
        html.Div(id='ctrl-slice', children=[
            _section('Slice',
                _lbl('Axis'),
                dcc.RadioItems(
                    id='slice-axis', options=['x', 'y', 'z'], value='z',
                    inline=True, style={'fontSize': '12px', 'color': _TEXT},
                    inputStyle={'marginRight': '4px'},
                    labelStyle={'marginRight': '14px', 'color': _TEXT},
                ),
                _lbl('Position (0 = start, 1 = end)'),
                _slider('slice-pos', 0.0, 1.0, 0.01, 0.5,
                        marks={0: '0', 0.5: '½', 1: '1'}),
            ),
        ]),

        # ── Isosurface ────────────────────────────────────────────────
        html.Div(id='ctrl-iso', style={'display': 'none'}, children=[
            _section('Isosurface',
                _lbl('Min threshold (fraction of color range)'),
                _slider('iso-min', 0.0, 1.0, 0.01, 0.2,
                        marks={0: '0', 0.5: '0.5', 1: '1'}),
                _lbl('Max threshold (fraction of color range)'),
                _slider('iso-max', 0.0, 1.0, 0.01, 0.8,
                        marks={0: '0', 0.5: '0.5', 1: '1'}),
                _lbl('Number of isosurfaces'),
                _slider('iso-n', 1, 10, 1, 3,
                        marks={1: '1', 3: '3', 5: '5', 10: '10'}),
                _lbl('Opacity'),
                _slider('iso-opacity', 0.05, 1.0, 0.05, 1.0,
                        marks={0.05: '0', 0.5: '0.5', 1.0: '1'}),
                dcc.Checklist(
                    id='iso-caps',
                    options=[{'label': '  Show end-caps', 'value': 'caps'}], value=[],
                    style={'fontSize': '12px', 'color': _SUB, 'marginTop': '6px'},
                ),
                dcc.Checklist(
                    id='iso-solid-chk',
                    options=[{'label': '  Solid color', 'value': 'solid'}], value=[],
                    style={'fontSize': '12px', 'color': _SUB, 'marginTop': '6px'},
                ),
                html.Div(id='ctrl-iso-solid-color', style={'display': 'none'}, children=[
                    _lbl('Surface color'),
                    dcc.Dropdown(
                        id='iso-solid-color',
                        options=[
                            {'label': 'white',  'value': 'white'},
                            {'label': 'gray',   'value': '#6c7086'},
                            {'label': 'red',    'value': '#f38ba8'},
                            {'label': 'cyan',   'value': '#94e2d5'},
                            {'label': 'yellow', 'value': '#f9e2af'},
                            {'label': 'green',  'value': '#a6e3a1'},
                            {'label': 'blue',   'value': '#89b4fa'},
                            {'label': 'mauve',  'value': '#cba6f7'},
                        ],
                        value='#89b4fa', clearable=False,
                        style={'fontSize': '12px', 'backgroundColor': _OVER,
                               'border': f'1px solid {_BORDER}'},
                    ),
                ]),
            ),
        ]),

        # ── Volume ────────────────────────────────────────────────────
        html.Div(id='ctrl-vol', style={'display': 'none'}, children=[
            _section('Volume',
                _lbl('Opacity per shell'),
                _slider('vol-opacity', 0.01, 0.5, 0.01, 0.1,
                        marks={0.01: '0', 0.25: '.25', 0.5: '.5'}),
                _lbl('Number of shells'),
                _slider('vol-nsurf', 3, 30, 1, 15,
                        marks={3: '3', 15: '15', 30: '30'}),
                _lbl('Isomin (fraction of color range)'),
                _slider('vol-min', 0.0, 1.0, 0.01, 0.0,
                        marks={0: '0', 0.5: '0.5', 1: '1'}),
                _lbl('Isomax (fraction of color range)'),
                _slider('vol-max', 0.0, 1.0, 0.01, 1.0,
                        marks={0: '0', 0.5: '0.5', 1: '1'}),
            ),
        ]),

        # ── Color ─────────────────────────────────────────────────────
        _section('Color',
            _lbl('Colormap'),
            dcc.Dropdown(
                id='cmap-sel', options=cmap_opts, value='viridis', clearable=False,
                style={'fontSize': '12px', 'backgroundColor': _OVER,
                       'border': f'1px solid {_BORDER}'},
            ),
            dcc.Checklist(
                id='log-chk',
                options=[{'label': '  Log scale', 'value': 'log'}], value=[],
                style={'fontSize': '12px', 'color': _SUB, 'marginTop': '6px'},
            ),
            html.Div([
                html.Div([_lbl('vmin'), _num('vmin-inp')],
                         style={'flex': 1, 'marginRight': '6px'}),
                html.Div([_lbl('vmax'), _num('vmax-inp')],
                         style={'flex': 1}),
            ], style={'display': 'flex'}),
            html.Button('↺  Auto range', id='reset-btn', n_clicks=0, style={
                'marginTop': '8px', 'padding': '4px 8px', 'fontSize': '11px',
                'width': '100%', 'backgroundColor': _OVER, 'color': _TEAL,
                'border': f'1px solid {_BORDER}', 'borderRadius': '4px',
                'cursor': 'pointer', 'fontFamily': 'monospace',
            }),
        ),

        # ── Overlay ───────────────────────────────────────────────────
        html.Div(id='ctrl-overlay', children=[
            _section('Overlay',
                _lbl('Variable'),
                dcc.Dropdown(
                    id='overlay-var-sel',
                    options=[{'label': 'None', 'value': '__none__'}] + var_opts,
                    value='__none__', clearable=False,
                    style={'fontSize': '12px', 'backgroundColor': _OVER,
                           'border': f'1px solid {_BORDER}'},
                ),
                _lbl('Levels'),
                _slider('overlay-nlevels', 1, 20, 1, 5,
                        marks={1: '1', 5: '5', 10: '10', 20: '20'}),
                _lbl('Color'),
                dcc.Dropdown(
                    id='overlay-color-sel',
                    options=[
                        {'label': 'white',  'value': 'white'},
                        {'label': 'black',  'value': 'black'},
                        {'label': 'red',    'value': '#f38ba8'},
                        {'label': 'cyan',   'value': '#94e2d5'},
                        {'label': 'yellow', 'value': '#f9e2af'},
                        {'label': 'green',  'value': '#a6e3a1'},
                        {'label': 'blue',   'value': '#89b4fa'},
                        {'label': 'mauve',  'value': '#cba6f7'},
                    ],
                    value='white', clearable=False,
                    style={'fontSize': '12px', 'backgroundColor': _OVER,
                           'border': f'1px solid {_BORDER}'},
                ),
                # Contour-specific: line width (hidden in isosurface/volume mode)
                html.Div(id='ctrl-overlay-contour', children=[
                    _lbl('Line width'),
                    _slider('overlay-lw', 0.5, 3.0, 0.5, 1.0,
                            marks={0.5: '0.5', 1.0: '1', 2.0: '2', 3.0: '3'}),
                ]),
                # Overlay mode selector for 3D isosurface/volume modes
                html.Div(id='ctrl-overlay-mode', style={'display': 'none'}, children=[
                    _lbl('Overlay type'),
                    dcc.Dropdown(
                        id='overlay-mode-sel',
                        options=[
                            {'label': 'Isosurface', 'value': 'isosurface'},
                            {'label': 'Isovolume',  'value': 'volume'},
                        ],
                        value='isosurface', clearable=False,
                        style={'fontSize': '12px', 'backgroundColor': _OVER,
                               'border': f'1px solid {_BORDER}'},
                    ),
                ]),
                # Isosurface overlay controls
                html.Div(id='ctrl-overlay-iso', style={'display': 'none'}, children=[
                    _lbl('Min threshold (fraction of range)'),
                    _slider('overlay-iso-min', 0.0, 1.0, 0.01, 0.2,
                            marks={0: '0', 0.5: '0.5', 1: '1'}),
                    _lbl('Max threshold (fraction of range)'),
                    _slider('overlay-iso-max', 0.0, 1.0, 0.01, 0.8,
                            marks={0: '0', 0.5: '0.5', 1: '1'}),
                    _lbl('Opacity'),
                    _slider('overlay-iso-opacity', 0.05, 1.0, 0.05, 0.6,
                            marks={0.05: '0', 0.5: '0.5', 1.0: '1'}),
                    dcc.Checklist(
                        id='overlay-iso-byval',
                        options=[{'label': '  Color by value', 'value': 'byval'}],
                        value=[],
                        style={'fontSize': '12px', 'color': _SUB, 'marginTop': '6px'},
                    ),
                ]),
                # Isovolume overlay controls
                html.Div(id='ctrl-overlay-vol', style={'display': 'none'}, children=[
                    _lbl('Min threshold (fraction of range)'),
                    _slider('overlay-vol-min', 0.0, 1.0, 0.01, 0.0,
                            marks={0: '0', 0.5: '0.5', 1: '1'}),
                    _lbl('Max threshold (fraction of range)'),
                    _slider('overlay-vol-max', 0.0, 1.0, 0.01, 1.0,
                            marks={0: '0', 0.5: '0.5', 1: '1'}),
                    _lbl('Opacity per shell'),
                    _slider('overlay-vol-opacity', 0.01, 0.5, 0.01, 0.1,
                            marks={0.01: '0', 0.25: '.25', 0.5: '.5'}),
                    _lbl('Number of shells'),
                    _slider('overlay-vol-nsurf', 3, 30, 1, 15,
                            marks={3: '3', 15: '15', 30: '30'}),
                ]),
            ),
        ], style={'display': 'block' if ndim >= 2 else 'none'}),

        # ── Status ────────────────────────────────────────────────────
        html.Div(id='status-bar', style={
            'marginTop': 'auto', 'paddingTop': '12px',
            'fontSize': '11px', 'color': _MUTED,
            'borderTop': f'1px solid {_OVER}', 'lineHeight': '1.7',
        }),

    ], style={
        'width': '265px', 'minWidth': '265px',
        'backgroundColor': _SURF, 'padding': '14px',
        'height': '100vh', 'overflowY': 'auto',
        'display': 'flex', 'flexDirection': 'column',
        'fontFamily': 'monospace', 'color': _TEXT,
        'boxSizing': 'border-box',
    })

    app.layout = html.Div([
        sidebar,
        html.Div([
            dcc.Graph(
                id='viz-graph', style={'height': '100vh'},
                config={
                    'displayModeBar': True, 'scrollZoom': True,
                    'modeBarButtonsToRemove': ['select2d', 'lasso2d'],
                    'toImageButtonOptions': {'format': 'png', 'scale': 2},
                },
            ),
            html.Img(
                id='pv-img',
                style={
                    'display': 'none', 'width': '100%', 'height': '100vh',
                    'objectFit': 'contain', 'backgroundColor': _BG,
                    'cursor': 'grab',
                },
            ),
        ], style={'flex': '1', 'overflow': 'hidden', 'backgroundColor': _BG,
                  'position': 'relative'}),

        dcc.Interval(id='play-iv', interval=500, n_intervals=0, disabled=True),
        dcc.Store(id='playing-st', data=False),
        dcc.Store(id='camera-3d', data=None),
    ], style={
        'display': 'flex', 'height': '100vh', 'overflow': 'hidden',
        'backgroundColor': _BG, 'fontFamily': 'monospace',
    })

    # ------------------------------------------------------------------
    # Callbacks
    # ------------------------------------------------------------------

    # Playback flag — set by _toggle_play, read by _update.
    # Avoids adding State('playing-st') to _update (which triggers Dash 4
    # initial-call serialization issues).
    _is_playing = [False]

    @app.callback(
        Output('play-iv', 'disabled'),
        Output('play-iv', 'interval'),
        Output('playing-st', 'data'),
        Output('play-btn', 'children'),
        Input('play-btn', 'n_clicks'),
        Input('stop-btn', 'n_clicks'),
        Input('fps-sl', 'value'),
        State('playing-st', 'data'),
        prevent_initial_call=True,
    )
    def _toggle_play(_, __, fps, is_playing):  # pylint: disable=unused-argument
        iv = max(int(1000 / max(float(fps or 2), 0.1)), 50)
        trig = (callback_context.triggered or [{}])[0].get('prop_id', '')
        if 'stop-btn' in trig:
            _is_playing[0] = False
            return True, iv, False, '▶  Play'
        if 'play-btn' in trig:
            playing = not is_playing
            _is_playing[0] = playing
            return not playing, iv, playing, ('⏸  Pause' if playing else '▶  Play')
        return not is_playing, iv, is_playing, no_update  # fps-only change

    # Playback advances the slider position.
    # _last_update_t rate-limits advances so the browser has time to render
    # each frame before the next is sent.  Without this, the dcc.Interval
    # fires faster than the browser can process, causing frame skipping.
    # _min_frame_gap adapts to actual render time: measured as the wall-clock
    # time between the start of _update and the next _advance_step call.
    _last_update_t = [0.0]
    _min_frame_gap = [0.3]
    # PyVista availability — probed once at startup, disabled on failure.
    # Mutable list so callbacks can flip it off on runtime errors.
    _pv_ok = [False]
    if _PV_AVAILABLE and ndim == 3:
        cons.print('[dim]Probing PyVista off-screen rendering...[/dim]')
        _pv_ok[0] = _pv_probe()
        if _pv_ok[0]:
            cons.print('[dim][green]PyVista OK[/green] — '
                       'server-side 3D rendering enabled.[/dim]')
        else:
            cons.print('[dim][yellow]PyVista probe failed[/yellow] — '
                       'falling back to Plotly WebGL for 3D playback.[/dim]')

    @app.callback(
        Output('step-sl', 'value'),
        Input('play-iv', 'n_intervals'),
        State('step-sl', 'value'),
        State('loop-chk', 'value'),
        prevent_initial_call=True,
    )
    def _advance_step(_, current_idx, loop_val):
        now = time.perf_counter()
        if now - _last_update_t[0] < _min_frame_gap[0]:
            return no_update
        _last_update_t[0] = now
        idx = int(current_idx or 0)
        nxt = idx + 1
        if nxt >= len(steps):
            return 0 if ('loop' in (loop_val or [])) else no_update
        return nxt

    # Slider index → actual step value + label
    @app.callback(
        Output('step-sel', 'data'),
        Output('step-label', 'children'),
        Input('step-sl', 'value'),
    )
    def _sync_step(sl_idx):
        idx = int(sl_idx) if sl_idx is not None else 0
        idx = max(0, min(idx, len(steps) - 1))
        return steps[idx], f'step {steps[idx]}'

    @app.callback(
        Output('ctrl-slice',           'style'),
        Output('ctrl-iso',             'style'),
        Output('ctrl-vol',             'style'),
        Output('ctrl-overlay',         'style'),
        Output('ctrl-overlay-contour', 'style'),
        Output('ctrl-overlay-mode',    'style'),
        Input('mode-sel', 'value'),
    )
    def _toggle_controls(mode):
        show, hide = {'display': 'block'}, {'display': 'none'}
        overlay_ok = mode in ('heatmap', 'slice', 'isosurface', 'volume') and ndim >= 2
        is_3d_surf = mode in ('isosurface', 'volume')
        return (
            show if mode == 'slice'      else hide,
            show if mode == 'isosurface' else hide,
            show if mode == 'volume'     else hide,
            show if overlay_ok           else hide,
            show if not is_3d_surf       else hide,   # contour line controls
            show if is_3d_surf           else hide,   # overlay type selector
        )

    @app.callback(
        Output('ctrl-overlay-iso', 'style'),
        Output('ctrl-overlay-vol', 'style'),
        Input('mode-sel', 'value'),
        Input('overlay-mode-sel', 'value'),
    )
    def _toggle_overlay_subs(mode, ov_mode):
        show, hide = {'display': 'block'}, {'display': 'none'}
        is_3d_surf = mode in ('isosurface', 'volume')
        if not is_3d_surf:
            return hide, hide
        ov = ov_mode or 'isosurface'
        return (
            show if ov == 'isosurface' else hide,
            show if ov == 'volume'     else hide,
        )

    @app.callback(
        Output('ctrl-iso-solid-color', 'style'),
        Input('iso-solid-chk', 'value'),
    )
    def _toggle_iso_solid(chk):
        show, hide = {'display': 'block'}, {'display': 'none'}
        return show if chk and 'solid' in chk else hide

    @app.callback(
        Output('vmin-inp', 'value'),
        Output('vmax-inp', 'value'),
        Input('reset-btn', 'n_clicks'),
        prevent_initial_call=True,
    )
    def _reset_range(_reset):
        return None, None

    @app.callback(
        Output('viz-graph',  'figure'),
        Output('status-bar', 'children'),
        Output('pv-img',     'src'),
        Output('pv-img',     'style'),
        Output('viz-graph',  'style'),
        Input('var-sel',     'value'),
        Input('step-sel',    'data'),
        Input('mode-sel',    'value'),
        Input('slice-axis',  'value'),
        Input('slice-pos',   'value'),
        Input('iso-min',     'value'),
        Input('iso-max',     'value'),
        Input('iso-n',       'value'),
        Input('iso-opacity', 'value'),
        Input('iso-caps',    'value'),
        Input('vol-opacity', 'value'),
        Input('vol-nsurf',   'value'),
        Input('vol-min',     'value'),
        Input('vol-max',     'value'),
        Input('cmap-sel',    'value'),
        Input('log-chk',     'value'),
        Input('vmin-inp',    'value'),
        Input('vmax-inp',    'value'),
        Input('iso-solid-chk',     'value'),
        Input('iso-solid-color',   'value'),
        Input('overlay-var-sel',   'value'),
        Input('overlay-nlevels',   'value'),
        Input('overlay-color-sel', 'value'),
        Input('overlay-lw',        'value'),
        Input('overlay-mode-sel',      'value'),
        Input('overlay-iso-min',     'value'),
        Input('overlay-iso-max',     'value'),
        Input('overlay-iso-opacity', 'value'),
        Input('overlay-iso-byval',   'value'),
        Input('overlay-vol-min',     'value'),
        Input('overlay-vol-max',     'value'),
        Input('overlay-vol-opacity', 'value'),
        Input('overlay-vol-nsurf',   'value'),
        Input('playing-st',          'data'),
    )
    def _update(var_sel, step, mode,  # pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals,too-many-branches,too-many-statements
                slice_axis, slice_pos,
                iso_min_frac, iso_max_frac, iso_n, iso_opacity, iso_caps,
                vol_opacity, vol_nsurf, vol_min_frac, vol_max_frac,
                cmap, log_chk, vmin_in, vmax_in,
                iso_solid_chk, iso_solid_color,
                overlay_var, overlay_nlevels, overlay_color, overlay_lw,
                overlay_mode_sel,
                overlay_iso_min, overlay_iso_max, overlay_iso_opacity,
                overlay_iso_byval,
                overlay_vol_min, overlay_vol_max, overlay_vol_opacity,
                overlay_vol_nsurf,
                playing_st):

        _t0 = time.perf_counter()
        _GRAPH_SHOW = {'height': '100vh', 'display': 'block'}
        _GRAPH_HIDE = {'height': '100vh', 'display': 'none'}
        _PV_SHOW = {
            'display': 'block', 'width': '100%', 'height': '100vh',
            'objectFit': 'contain', 'backgroundColor': _BG, 'cursor': 'grab',
        }
        _PV_HIDE = {
            'display': 'none', 'width': '100%', 'height': '100vh',
            'objectFit': 'contain', 'backgroundColor': _BG, 'cursor': 'grab',
        }
        selected_var = var_sel or varname

        # When the variable selector is the trigger, ignore any stale manual
        # range values — they belong to the previous variable.
        trig = (callback_context.triggered or [{}])[0].get('prop_id', '')
        if 'var-sel' in trig:
            vmin_in = vmax_in = None

        # Use single-variable loading when available: reads only the selected
        # variable from disk instead of all variables, giving an N_vars speedup.
        # Cache key is (step, var) so switching variables is also a fast read.
        try:
            if read_one_var_func is not None:
                ad = _load((step, selected_var),
                           lambda key: read_one_var_func(key[0], key[1]))
            else:
                ad = _load(step, read_func)
        except (OSError, ValueError, EOFError) as exc:
            return (no_update,
                    [html.Span(f' Error loading step {step}: {exc}',
                               style={'color': _RED})],
                    no_update, no_update, no_update)
        _t_load = time.perf_counter()

        # Resolve log/cmap early — needed by the prefetch JPEG encoder below.
        log = bool(log_chk and 'log' in log_chk)
        cmap = cmap or 'viridis'

        # Eagerly pre-load the next steps in the background so that
        # navigation feels instant after the structure cache is warm.
        # prefetch_one is a no-op if the key is already cached or in-flight.
        # For 2D, also pre-encode the JPEG at the current colormap / range
        # settings so the next step callback can skip encoding entirely.
        # For 3D, pre-compute isomesh or volume data so playback is instant.
        # Pre-load 4 steps ahead for 3D (marching cubes takes longer), 2 for 2D.
        _pf_depth = 4 if (ad.ndim == 3 and mode in ('isosurface', 'volume')) else 2
        if read_one_var_func is not None:
            try:
                _idx = steps.index(step)
                for _ns in steps[_idx + 1: _idx + 1 + _pf_depth]:
                    _nk = (_ns, selected_var)
                    _prefetch_one(_nk, lambda k=_nk: read_one_var_func(k[0], k[1]))
                    if ad.ndim == 2:
                        _sv = selected_var
                        _prefetch_jpeg(
                            _ns, _sv,
                            lambda s, sv=_sv: _step_cache.load(
                                (s, sv), lambda k: read_one_var_func(k[0], k[1])),
                            cmap, vmin_in, vmax_in, log,
                        )
                    elif ad.ndim == 3 and mode in ('isosurface', 'volume'):
                        _sv = selected_var
                        _prefetch_3d_mesh(
                            _ns, _sv,
                            lambda s, sv=_sv: _step_cache.load(
                                (s, sv), lambda k: read_one_var_func(k[0], k[1])),
                            mode, log, vmin_in, vmax_in,
                            float(iso_min_frac or 0.2), float(iso_max_frac or 0.8),
                            int(iso_n or 3),
                            float(vol_min_frac or 0.0), float(vol_max_frac or 1.0),

                        )
            except (ValueError, IndexError):
                pass
        else:
            try:
                _idx = steps.index(step)
                for _ns in steps[_idx + 1: _idx + 1 + _pf_depth]:
                    _prefetch_one(_ns, read_func)
                    if ad.ndim == 2:
                        _prefetch_jpeg(
                            _ns, selected_var,
                            lambda s: _step_cache.load(s, read_func),
                            cmap, vmin_in, vmax_in, log,
                        )
                    elif ad.ndim == 3 and mode in ('isosurface', 'volume'):
                        _prefetch_3d_mesh(
                            _ns, selected_var,
                            lambda s: _step_cache.load(s, read_func),
                            mode, log, vmin_in, vmax_in,
                            float(iso_min_frac or 0.2), float(iso_max_frac or 0.8),
                            int(iso_n or 3),
                            float(vol_min_frac or 0.0), float(vol_max_frac or 1.0),

                        )
            except (ValueError, IndexError):
                pass

        if selected_var not in ad.variables:
            avail = ', '.join(sorted(ad.variables))
            return (no_update,
                    [html.Span(
                        f' Variable {selected_var!r} not in step {step} '
                        f'(available: {avail})', style={'color': _RED})],
                    no_update, no_update, no_update)
        raw = ad.variables[selected_var]

        # Color range — subsample large arrays for speed (nanmin/nanmax on
        # a 36M element 3D array takes 44–183ms; on a 200K subsample ~0.5ms).
        # A regular stride subsample captures extrema reliably for smooth CFD
        # fields, and the full array is still used for all actual computations.
        if raw.size > 200_000:
            _sr = max(1, math.ceil((raw.size / 200_000) ** (1.0 / raw.ndim)))
            _slc = tuple(slice(None, None, _sr) for _ in range(raw.ndim))
            _raw_range = raw[_slc]
        else:
            _raw_range = raw

        if vmin_in is not None:
            vmin = float(vmin_in)
        else:
            _safe = _raw_range[_raw_range > 0] if log and np.any(_raw_range > 0) else _raw_range
            vmin = float(np.nanmin(_safe))
        if vmax_in is not None:
            vmax = float(vmax_in)
        else:
            vmax = float(np.nanmax(_raw_range))
        if vmax <= vmin:
            vmax = vmin + 1e-10

        if log:
            def _tf(arr):
                return np.where(arr > 0, np.log10(np.maximum(arr, 1e-300)), np.nan)
            cmin = float(np.log10(max(vmin, 1e-300)))
            cmax = float(np.log10(max(vmax, 1e-300)))
            cbar_title = f'log\u2081\u2080({selected_var})'
        else:
            def _tf(arr): return arr
            cmin, cmax = vmin, vmax
            cbar_title = selected_var

        _t_prep = time.perf_counter()

        # ----------------------------------------------------------------------
        # PyVista server-side rendering for 3D during playback.
        # Renders to a JPEG displayed via html.Img — bypasses Plotly WebGL
        # entirely, giving ~50-100ms full mesh swap and ~5ms camera re-renders.
        # ----------------------------------------------------------------------
        _use_pv = (
            _pv_ok[0]
            and bool(playing_st)
            and ad.ndim == 3
            and mode in ('isosurface', 'volume')
        )
        if _use_pv:
            _has_ov_pv = overlay_var and overlay_var != '__none__'
            _ov_raw_pv = None
            _ov_x_pv = _ov_y_pv = _ov_z_pv = None
            if _has_ov_pv:
                if overlay_var in ad.variables:
                    _ov_raw_pv = ad.variables[overlay_var]
                    _ov_x_pv, _ov_y_pv, _ov_z_pv = ad.x_cc, ad.y_cc, ad.z_cc
                elif read_one_var_func is not None:
                    try:
                        _ov_ad_pv = _load(
                            (step, overlay_var),
                            lambda k: read_one_var_func(k[0], k[1]))
                        if overlay_var in _ov_ad_pv.variables:
                            _ov_raw_pv = _ov_ad_pv.variables[overlay_var]
                            _ov_x_pv = _ov_ad_pv.x_cc
                            _ov_y_pv = _ov_ad_pv.y_cc
                            _ov_z_pv = _ov_ad_pv.z_cc
                    except (OSError, ValueError, EOFError):
                        pass
                else:
                    try:
                        _ov_ad_pv = _load(step, read_func)
                        if overlay_var in _ov_ad_pv.variables:
                            _ov_raw_pv = _ov_ad_pv.variables[overlay_var]
                            _ov_x_pv = _ov_ad_pv.x_cc
                            _ov_y_pv = _ov_ad_pv.y_cc
                            _ov_z_pv = _ov_ad_pv.z_cc
                    except (OSError, ValueError, EOFError):
                        pass

            try:
                pv_src = _pv_render(
                    raw, ad.x_cc, ad.y_cc, ad.z_cc,
                    mode, cmap, _tf, cmin, cmax,
                    float(iso_min_frac or 0.2), float(iso_max_frac or 0.8),
                    int(iso_n or 3),
                    float(vol_opacity or 0.1), int(vol_nsurf or 15),
                    float(vol_min_frac or 0.0), float(vol_max_frac or 1.0),
                    None,  # camera — use default iso view
                    step, selected_var,
                    overlay_raw=_ov_raw_pv,
                    overlay_x=_ov_x_pv, overlay_y=_ov_y_pv,
                    overlay_z=_ov_z_pv,
                    overlay_mode=overlay_mode_sel or 'isosurface',
                    overlay_iso_min=float(overlay_iso_min or 0.2),
                    overlay_iso_max=float(overlay_iso_max or 0.8),
                    overlay_nlevels=int(overlay_nlevels or 3),
                    overlay_color=overlay_color or 'white',
                    overlay_opacity=float(overlay_iso_opacity or 0.6),
                    overlay_vol_min=float(overlay_vol_min or 0.0),
                    overlay_vol_max=float(overlay_vol_max or 1.0),
                    overlay_vol_opacity=float(overlay_vol_opacity or 0.1),
                    overlay_vol_nsurf=int(overlay_vol_nsurf or 15),
                )
            except Exception as _pv_exc:  # pylint: disable=broad-except
                # PyVista rendering failed — disable for the rest of the
                # session and fall through to the Plotly WebGL path.
                _pv_ok[0] = False
                _pv_disabled.set()  # kill switch for concurrent callbacks
                cons.print(
                    '[dim][yellow]PyVista render failed, falling back to '
                    f'Plotly:[/yellow] {_pv_exc}[/dim]')
                pv_src = None
                _t_prep = time.perf_counter()  # reset so Plotly timing is clean

            if pv_src is not None:
                _t_pv = time.perf_counter()
                dmin_pv = float(np.nanmin(_raw_range))
                dmax_pv = float(np.nanmax(_raw_range))
                cons.print(
                    f'[dim]viz timing  step={step}  shape={raw.shape}'
                    f'  load={_t_load-_t0:.3f}s'
                    f'  prep={_t_prep-_t_load:.3f}s'
                    f'  pyvista={_t_pv-_t_prep:.3f}s'
                    f'  total={_t_pv-_t0:.3f}s [PV-3D][/dim]'
                )
                status_pv = html.Div([
                    html.Span(f'step {step}', style={'color': _YELLOW}),
                    html.Span(f'  ·  shape {raw.shape}', style={'color': _MUTED}),
                    html.Br(),
                    html.Span('min ', style={'color': _MUTED}),
                    html.Span(f'{dmin_pv:.4g}', style={'color': _BLUE}),
                    html.Span('  max ', style={'color': _MUTED}),
                    html.Span(f'{dmax_pv:.4g}', style={'color': _RED}),
                ])
                _now_pv = time.perf_counter()
                _last_update_t[0] = _now_pv
                _server_pv = _now_pv - _t0
                _min_frame_gap[0] = max(0.15, min(1.0, _server_pv + 0.05))
                return no_update, status_pv, pv_src, _PV_SHOW, _GRAPH_HIDE

        fig = go.Figure()
        title = ''

        if ad.ndim == 3:
            # ------------------------------------------------------------------
            # 3D Patch() fast path: on step / vmin / vmax changes only update
            # the data arrays without rebuilding the full figure.  Slice mode
            # also patches surfacecolor; iso/vol patch value + thresholds.
            # Colormap, variable, mode, and axis changes always trigger a full
            # render because they require new coordinate arrays or trace types.
            # ------------------------------------------------------------------
            _has_overlay_3d = overlay_var and overlay_var != '__none__'
            _trig3 = {t.get('prop_id', '') for t in (callback_context.triggered or [])}
            _PT_BASE = {'step-sel.data', 'vmin-inp.value', 'vmax-inp.value'}
            _PT_ISO  = _PT_BASE | {'iso-min.value', 'iso-max.value',
                                   'iso-n.value', 'iso-caps.value',
                                   'iso-opacity.value'}
            _PT_VOL  = _PT_BASE | {'vol-min.value', 'vol-max.value',
                                   'vol-opacity.value', 'vol-nsurf.value'}
            # Slice mode always does a full render — Plotly does not reliably
            # re-render go.Surface when surfacecolor is updated via Patch().
            # Overlay forces full render (adds extra traces).
            _do_patch_3d = (
                _trig3
                and '.' not in _trig3
                and not (_has_overlay_3d and mode in ('isosurface', 'volume'))
                and (
                    (mode == 'isosurface' and _trig3.issubset(_PT_ISO))  or
                    (mode == 'volume'     and _trig3.issubset(_PT_VOL))
                )
            )
            if _do_patch_3d:
                _cscale3 = _lut_to_plotly_colorscale(cmap)
                rng3 = cmax - cmin if cmax > cmin else 1.0
                patch = Patch()
                _cache_hit = False
                if mode == 'isosurface':
                    ilo = cmin + rng3 * float(iso_min_frac or 0.2)
                    ihi = cmin + rng3 * max(float(iso_max_frac or 0.8), ilo + 0.01)
                    # Try pre-computed mesh first, fall back to computing now
                    _cached = _get_cached_3d_mesh(
                        step, selected_var, mode, log,
                        vmin_in, vmax_in,
                        float(iso_min_frac or 0.2), float(iso_max_frac or 0.8),
                        int(iso_n or 3),
                        float(vol_min_frac or 0.0), float(vol_max_frac or 1.0),
                    )
                    if _cached is not None:
                        vx, vy, vz, fi, fj, fk, intens = _cached
                        _cache_hit = True
                    else:
                        raw_ds, x_ds3, y_ds3, z_ds3 = _get_ds3(
                            step, selected_var, raw, ad.x_cc, ad.y_cc, ad.z_cc, 500_000,
                        )
                        vx, vy, vz, fi, fj, fk, intens = _compute_isomesh(
                            raw_ds, x_ds3, y_ds3, z_ds3, _tf, ilo, ihi,
                            int(iso_n or 3),
                        )
                    patch['data'][0]['x'] = vx
                    patch['data'][0]['y'] = vy
                    patch['data'][0]['z'] = vz
                    patch['data'][0]['i'] = fi
                    patch['data'][0]['j'] = fj
                    patch['data'][0]['k'] = fk
                    patch['data'][0]['intensity'] = intens
                    patch['data'][0]['cmin'] = ilo
                    patch['data'][0]['cmax'] = ihi
                    patch['data'][0]['colorscale'] = _cscale3
                    patch['data'][0]['opacity'] = float(iso_opacity or 1.0)
                else:  # volume
                    # Try pre-computed volume data first
                    _cached = _get_cached_3d_mesh(
                        step, selected_var, mode, log,
                        vmin_in, vmax_in,
                        float(iso_min_frac or 0.2), float(iso_max_frac or 0.8),
                        int(iso_n or 3),
                        float(vol_min_frac or 0.0), float(vol_max_frac or 1.0),
                    )
                    if _cached is not None:
                        vf, vlo, vhi, _, _ = _cached
                        _cache_hit = True
                    else:
                        raw_ds, _, _, _ = _get_ds3(
                            step, selected_var, raw, ad.x_cc, ad.y_cc, ad.z_cc, 150_000,
                        )
                        vf = _tf(raw_ds).ravel()
                        vlo = cmin + rng3 * float(vol_min_frac or 0.0)
                        vhi = cmin + rng3 * max(float(vol_max_frac or 1.0), vlo + 0.01)
                    patch['data'][0]['value'] = vf
                    patch['data'][0]['isomin'] = vlo
                    patch['data'][0]['isomax'] = vhi
                    patch['data'][0]['opacity'] = float(vol_opacity or 0.1)
                    patch['data'][0]['surface_count'] = int(vol_nsurf or 15)
                    patch['data'][0]['colorscale'] = _cscale3
                patch['layout']['title']['text'] = (
                    f'{selected_var}  ·  step {step}'
                )
                _t_trace = time.perf_counter()
                dmin3, dmax3 = float(np.nanmin(_raw_range)), float(np.nanmax(_raw_range))
                cons.print(
                    f'[dim]viz timing  step={step}  shape={raw.shape}'
                    f'  load={_t_load-_t0:.3f}s'
                    f'  prep={_t_prep-_t_load:.3f}s'
                    f'  patch={_t_trace-_t_prep:.3f}s [PATCH-3D{"·HIT" if _cache_hit else ""}]'
                    f'  total={_t_trace-_t0:.3f}s[/dim]'
                )
                status = html.Div([
                    html.Span(f'step {step}', style={'color': _YELLOW}),
                    html.Span(f'  ·  shape {raw.shape}', style={'color': _MUTED}),
                    html.Br(),
                    html.Span('min ', style={'color': _MUTED}),
                    html.Span(f'{dmin3:.4g}', style={'color': _BLUE}),
                    html.Span('  max ', style={'color': _MUTED}),
                    html.Span(f'{dmax3:.4g}', style={'color': _RED}),
                ])
                _now = time.perf_counter()
                _last_update_t[0] = _now
                # Patch path: browser overhead is low (no figure rebuild).
                _server_s = _now - _t0
                _min_frame_gap[0] = max(0.3, min(2.0, _server_s + 0.25))
                return patch, status, no_update, _PV_HIDE, _GRAPH_SHOW

            # Full 3D render — check mesh prefetch cache for isosurface mode
            _iso_solid = (iso_solid_color or '#89b4fa') if (
                iso_solid_chk and 'solid' in iso_solid_chk) else None
            _cached_primary = None
            if mode == 'isosurface':
                _cached_primary = _get_cached_3d_mesh(
                    step, selected_var, mode, log,
                    vmin_in, vmax_in,
                    float(iso_min_frac or 0.2), float(iso_max_frac or 0.8),
                    int(iso_n or 3),
                    float(vol_min_frac or 0.0), float(vol_max_frac or 1.0),
                )
            trace, title = _build_3d(
                ad, raw, selected_var, step, mode, cmap, _tf, cmin, cmax, cbar_title,
                slice_axis or 'z', float(slice_pos or 0.5),
                float(iso_min_frac or 0.2), float(iso_max_frac or 0.8),
                int(iso_n or 3), bool(iso_caps and 'caps' in iso_caps),
                float(vol_opacity or 0.1), int(vol_nsurf or 15),
                float(vol_min_frac or 0.0), float(vol_max_frac or 1.0),
                iso_solid_color=_iso_solid,
                iso_opacity=float(iso_opacity or 1.0),
                cached_mesh=_cached_primary,
            )
            fig.add_trace(trace)
            # Bubble overlay for 3D
            if bubble_func is not None:
                try:
                    bubbles = bubble_func(step)
                    if bubbles is not None and len(bubbles) > 0:
                        if mode == 'slice':
                            s_axis = slice_axis or 'z'
                            s_col = {'x': 0, 'y': 1, 'z': 2}[s_axis]
                            ax_coords = {'x': ad.x_cc, 'y': ad.y_cc, 'z': ad.z_cc}[s_axis]
                            s_coord = ax_coords[0] + (ax_coords[-1] - ax_coords[0]) * float(slice_pos or 0.5)
                            near = np.abs(bubbles[:, s_col] - s_coord) <= bubbles[:, 3]
                            vis = bubbles[near] if near.any() else None
                        else:
                            vis = bubbles
                        if vis is not None and len(vis) > 0:
                            fig.add_trace(go.Scatter3d(
                                x=vis[:, 0], y=vis[:, 1], z=vis[:, 2],
                                mode='markers',
                                marker=dict(size=4, color='white', opacity=0.6, symbol='circle'),
                                showlegend=False,
                                hovertemplate='x=%{x:.3g}<br>y=%{y:.3g}<br>z=%{z:.3g}<extra>bubble</extra>',
                            ))
                except (OSError, ValueError):
                    pass  # bubble overlay is best-effort; skip on read errors
            # Overlay for 3D: contour lines (slice), isosurface, or isovolume
            _has_overlay = overlay_var and overlay_var != '__none__'
            if _has_overlay and mode in ('slice', 'isosurface', 'volume'):
                _ov_3d_raw = None
                _ov_3d_ad = ad
                if overlay_var in ad.variables:
                    _ov_3d_raw = ad.variables[overlay_var]
                elif read_one_var_func is not None:
                    try:
                        _ov_3d_ad = _load((step, overlay_var),
                                          lambda k: read_one_var_func(k[0], k[1]))
                        if overlay_var in _ov_3d_ad.variables:
                            _ov_3d_raw = _ov_3d_ad.variables[overlay_var]
                    except (OSError, ValueError, EOFError):
                        pass
                else:
                    try:
                        _ov_3d_ad = _load(step, read_func)
                        if overlay_var in _ov_3d_ad.variables:
                            _ov_3d_raw = _ov_3d_ad.variables[overlay_var]
                    except (OSError, ValueError, EOFError):
                        pass
                if _ov_3d_raw is not None and mode == 'slice':
                    _ov3_traces = _compute_contour_traces_3d(
                        _ov_3d_raw, _ov_3d_ad.x_cc, _ov_3d_ad.y_cc, _ov_3d_ad.z_cc,
                        slice_axis or 'z', float(slice_pos or 0.5),
                        int(overlay_nlevels or 5),
                        overlay_color or 'white',
                        float(overlay_lw or 1.0),
                    )
                    for _ct in _ov3_traces:
                        fig.add_trace(_ct)
                elif _ov_3d_raw is not None and mode in ('isosurface', 'volume'):
                    _ov_type = overlay_mode_sel or 'isosurface'
                    _ov_vmin = float(np.nanmin(_ov_3d_raw))
                    _ov_vmax = float(np.nanmax(_ov_3d_raw))
                    _ov_rng = _ov_vmax - _ov_vmin if _ov_vmax > _ov_vmin else 1.0
                    if _ov_type == 'isosurface':
                        # Overlay isosurfaces of the second variable
                        _ov_min_f = float(overlay_iso_min or 0.2)
                        _ov_max_f = float(overlay_iso_max or 0.8)
                        _ov_ilo = _ov_vmin + _ov_rng * _ov_min_f
                        _ov_ihi = _ov_vmin + _ov_rng * max(_ov_max_f, _ov_min_f + 0.01)
                        _ov_ds, _ox, _oy, _oz = _get_ds3(
                            step, overlay_var, _ov_3d_raw,
                            _ov_3d_ad.x_cc, _ov_3d_ad.y_cc, _ov_3d_ad.z_cc,
                            500_000,
                        )
                        _ov_vx, _ov_vy, _ov_vz, _ov_fi, _ov_fj, _ov_fk, _ov_int = \
                            _compute_isomesh(
                                _ov_ds, _ox, _oy, _oz,
                                _tf, _ov_ilo, _ov_ihi,
                                int(overlay_nlevels or 3),
                            )
                        _ov_byval = overlay_iso_byval and 'byval' in overlay_iso_byval
                        _ov_op = float(overlay_iso_opacity or 0.6)
                        if _ov_byval:
                            _ov_cscale = _lut_to_plotly_colorscale(cmap)
                            fig.add_trace(go.Mesh3d(
                                x=_ov_vx, y=_ov_vy, z=_ov_vz,
                                i=_ov_fi, j=_ov_fj, k=_ov_fk,
                                intensity=_ov_int, intensitymode='vertex',
                                colorscale=_ov_cscale, cmin=_ov_ilo, cmax=_ov_ihi,
                                colorbar=_make_cbar(overlay_var, _ov_ilo, _ov_ihi),
                                showscale=True, opacity=_ov_op,
                                lighting=dict(ambient=0.7, diffuse=0.9, specular=0.3,
                                              roughness=0.5, fresnel=0.2),
                                lightposition=dict(x=1000, y=500, z=500),
                                flatshading=False,
                            ))
                        else:
                            _ov_solid_cs = [[0, overlay_color or 'white'],
                                            [1, overlay_color or 'white']]
                            fig.add_trace(go.Mesh3d(
                                x=_ov_vx, y=_ov_vy, z=_ov_vz,
                                i=_ov_fi, j=_ov_fj, k=_ov_fk,
                                intensity=np.zeros(len(_ov_vx), dtype=np.float32),
                                intensitymode='vertex',
                                colorscale=_ov_solid_cs, cmin=0, cmax=1,
                                showscale=False, opacity=_ov_op,
                                lighting=dict(ambient=0.7, diffuse=0.9, specular=0.3,
                                              roughness=0.5, fresnel=0.2),
                                lightposition=dict(x=1000, y=500, z=500),
                                flatshading=False,
                            ))
                    else:
                        # Overlay isovolume of the second variable
                        _ov_vlo_f = float(overlay_vol_min or 0.0)
                        _ov_vhi_f = float(overlay_vol_max or 1.0)
                        _ov_vlo = _ov_vmin + _ov_rng * _ov_vlo_f
                        _ov_vhi = _ov_vmin + _ov_rng * max(_ov_vhi_f, _ov_vlo_f + 0.01)
                        _ov_ds, _ox, _oy, _oz = _get_ds3(
                            step, overlay_var, _ov_3d_raw,
                            _ov_3d_ad.x_cc, _ov_3d_ad.y_cc, _ov_3d_ad.z_cc,
                            150_000,
                        )
                        _ov_X, _ov_Y, _ov_Z = np.meshgrid(
                            _ox, _oy, _oz, indexing='ij')
                        _ov_vf = _tf(_ov_ds.ravel()).astype(np.float32)
                        _ov_vol_op = float(overlay_vol_opacity or 0.1)
                        _ov_vol_ns = int(overlay_vol_nsurf or 15)
                        _ov_cscale = _lut_to_plotly_colorscale(cmap)
                        fig.add_trace(go.Volume(
                            x=_ov_X.ravel().astype(np.float32),
                            y=_ov_Y.ravel().astype(np.float32),
                            z=_ov_Z.ravel().astype(np.float32),
                            value=_ov_vf,
                            isomin=_ov_vlo, isomax=_ov_vhi,
                            opacity=_ov_vol_op, surface_count=_ov_vol_ns,
                            colorscale=_ov_cscale, cmin=_ov_vmin, cmax=_ov_vmax,
                            colorbar=_make_cbar(overlay_var, _ov_vmin, _ov_vmax),
                        ))
            # Compute aspect ratio from domain extents so slices (which
            # have a constant coordinate on one axis) don't collapse that axis.
            dx = float(ad.x_cc[-1] - ad.x_cc[0]) if len(ad.x_cc) > 1 else 1.0
            dy = float(ad.y_cc[-1] - ad.y_cc[0]) if len(ad.y_cc) > 1 else 1.0
            dz = float(ad.z_cc[-1] - ad.z_cc[0]) if len(ad.z_cc) > 1 else 1.0
            max_d = max(dx, dy, dz, 1e-30)
            _xr = [float(ad.x_cc[0]), float(ad.x_cc[-1])]
            _yr = [float(ad.y_cc[0]), float(ad.y_cc[-1])]
            _zr = [float(ad.z_cc[0]), float(ad.z_cc[-1])]
            fig.update_layout(scene=dict(
                xaxis=dict(title='x', range=_xr, autorange=False,
                           backgroundcolor=_SURF, gridcolor=_OVER, color=_TEXT),
                yaxis=dict(title='y', range=_yr, autorange=False,
                           backgroundcolor=_SURF, gridcolor=_OVER, color=_TEXT),
                zaxis=dict(title='z', range=_zr, autorange=False,
                           backgroundcolor=_SURF, gridcolor=_OVER, color=_TEXT),
                bgcolor=_BG,
                aspectmode='manual',
                aspectratio=dict(x=dx/max_d, y=dy/max_d, z=dz/max_d),
            ))

        elif ad.ndim == 2:
            # Downsample for display.  The array may have more pixels than the
            # screen can show (e.g. 3601×801 on a 1920px monitor).  Striding
            # to ~1200×600 max reduces PNG encoding time by 5–10× with no
            # perceptible quality loss.  Color range is computed on the full
            # array above so min/max accuracy is preserved.
            _display_arr = _tf(raw)
            _nx, _ny = _display_arr.shape
            _sx = max(1, _nx // 1200)
            _sy = max(1, _ny // 600)
            if _sx > 1 or _sy > 1:
                _display_arr = _display_arr[::_sx, ::_sy]
                _x_disp = ad.x_cc[::_sx]
                _y_disp = ad.y_cc[::_sy]
            else:
                _x_disp = ad.x_cc
                _y_disp = ad.y_cc
            title = f'{selected_var}  ·  step {step}'

            # Patch fast path: when only step or color-range changes, skip
            # rebuilding the full Plotly figure.  Dash merges the diff with
            # the existing browser figure and Plotly only re-draws the image,
            # saving ~100–200 ms of browser-side figure reconstruction per step.
            # Disabled when contour overlay is active — contour traces change
            # with step and are hard to patch surgically.
            _has_overlay = overlay_var and overlay_var != '__none__'
            _trig = {t.get('prop_id', '') for t in (callback_context.triggered or [])}
            _PATCH_OK = {'step-sel.data', 'vmin-inp.value', 'vmax-inp.value'}
            _do_patch = (
                _trig                          # not initial render (empty set)
                and '.' not in _trig           # not Dash synthetic init trigger
                and _trig.issubset(_PATCH_OK)  # only step/range changed
                and not _has_overlay           # no contour overlay active
            )
            if _do_patch:
                # Check pre-encode cache first (populated by background prefetch)
                _jkey = (step, selected_var, cmap, vmin_in, vmax_in, log)
                with _jpeg_lock:
                    png_src = _jpeg_cache.get(_jkey)
                _jpeg_hit = png_src is not None
                if not _jpeg_hit:
                    png_src = _make_png_source(_display_arr.T, cmap, cmin, cmax)
                patch = Patch()
                patch['data'][0]['source'] = png_src
                # Update colorbar (data[1] = dummy scatter colorbar)
                _cb = _make_cbar(cbar_title, cmin, cmax)
                patch['data'][1]['marker']['colorscale'] = _lut_to_plotly_colorscale(cmap)
                patch['data'][1]['marker']['cmin'] = cmin
                patch['data'][1]['marker']['cmax'] = cmax
                patch['data'][1]['marker']['color'] = [cmin]
                patch['data'][1]['marker']['colorbar'] = _cb
                patch['layout']['title']['text'] = title
                # Update bubble overlay shapes (they change with step)
                if bubble_func is not None:
                    try:
                        bubbles = bubble_func(step)
                        if bubbles is not None and len(bubbles) > 0:
                            patch['layout']['shapes'] = [
                                dict(
                                    type='circle', xref='x', yref='y',
                                    x0=float(b[0] - b[3]), y0=float(b[1] - b[3]),
                                    x1=float(b[0] + b[3]), y1=float(b[1] + b[3]),
                                    line=dict(color='white', width=0.8),
                                    fillcolor='rgba(0,0,0,0)',
                                )
                                for b in bubbles
                            ]
                        else:
                            patch['layout']['shapes'] = []
                    except (OSError, ValueError):
                        patch['layout']['shapes'] = []
                _t_trace = time.perf_counter()
                dmin, dmax = float(np.nanmin(_raw_range)), float(np.nanmax(_raw_range))
                _tag = 'PATCH-JPEG-HIT' if _jpeg_hit else 'PATCH'
                cons.print(
                    f'[dim]viz timing  step={step}  shape={raw.shape}'
                    f'  load={_t_load-_t0:.3f}s'
                    f'  prep={_t_prep-_t_load:.3f}s'
                    f'  png={_t_trace-_t_prep:.3f}s [{_tag}]'
                    f'  payload={len(png_src)//1024}KB'
                    f'  total={_t_trace-_t0:.3f}s[/dim]'
                )
                status = html.Div([
                    html.Span(f'step {step}', style={'color': _YELLOW}),
                    html.Span(f'  ·  shape {raw.shape}', style={'color': _MUTED}),
                    html.Br(),
                    html.Span('min ', style={'color': _MUTED}),
                    html.Span(f'{dmin:.4g}', style={'color': _BLUE}),
                    html.Span('  max ', style={'color': _MUTED}),
                    html.Span(f'{dmax:.4g}', style={'color': _RED}),
                ])
                _last_update_t[0] = time.perf_counter()
                return patch, status, no_update, _PV_HIDE, _GRAPH_SHOW

            # Full render: initial load or structural change (colormap, variable,
            # log scale, etc.)
            # Use go.Image with LUT-colored RGB PNG (px.imshow with
            # binary_string=True only generates grayscale — no colormap).
            _dx_val = float((_x_disp[-1] - _x_disp[0]) / max(len(_x_disp) - 1, 1))
            _dy_val = float((_y_disp[-1] - _y_disp[0]) / max(len(_y_disp) - 1, 1))
            _jkey2 = (step, selected_var, cmap, vmin_in, vmax_in, log)
            with _jpeg_lock:
                png_src = _jpeg_cache.get(_jkey2)
            if png_src is None:
                png_src = _make_png_source(_display_arr.T, cmap, cmin, cmax)
            fig = go.Figure([
                go.Image(
                    source=png_src,
                    x0=float(_x_disp[0]), dx=_dx_val,
                    y0=float(_y_disp[0]), dy=_dy_val,
                    hoverinfo='skip',
                ),
                # Invisible scatter whose sole purpose is rendering the colorbar.
                # go.Image has no native colorscale support.
                go.Scatter(
                    x=[None], y=[None], mode='markers',
                    showlegend=False, hoverinfo='skip',
                    marker=dict(
                        colorscale=_lut_to_plotly_colorscale(cmap), cmin=cmin, cmax=cmax,
                        color=[cmin], showscale=True,
                        colorbar=_make_cbar(cbar_title, cmin, cmax),
                    ),
                ),
            ])
            fig.update_layout(
                xaxis=dict(title='x', color=_TEXT, gridcolor=_OVER,
                           scaleanchor='y', exponentformat='e'),
                yaxis=dict(title='y', color=_TEXT, gridcolor=_OVER,
                           exponentformat='e'),
                plot_bgcolor=_BG,
            )
            # Bubble overlay for 2D
            if bubble_func is not None:
                try:
                    bubbles = bubble_func(step)
                    if bubbles is not None and len(bubbles) > 0:
                        shapes = [
                            dict(
                                type='circle',
                                xref='x', yref='y',
                                x0=float(b[0] - b[3]), y0=float(b[1] - b[3]),
                                x1=float(b[0] + b[3]), y1=float(b[1] + b[3]),
                                line=dict(color='white', width=0.8),
                                fillcolor='rgba(0,0,0,0)',
                            )
                            for b in bubbles
                        ]
                        fig.update_layout(shapes=shapes)
                except (OSError, ValueError):
                    pass  # bubble overlay is best-effort; skip on read errors

            # Contour overlay for 2D
            if _has_overlay and overlay_var in ad.variables:
                _ov_raw = ad.variables[overlay_var]
                # If loaded via single-var mode, we need to load the overlay
                # variable separately.
                _ov_traces = _compute_contour_traces(
                    _ov_raw, ad.x_cc, ad.y_cc,
                    int(overlay_nlevels or 5),
                    overlay_color or 'white',
                    float(overlay_lw or 1.0),
                )
                for _ct in _ov_traces:
                    fig.add_trace(_ct)
            elif _has_overlay and overlay_var not in ad.variables:
                # Single-var loading mode: overlay var needs a separate load
                try:
                    if read_one_var_func is not None:
                        _ov_ad = _load((step, overlay_var),
                                       lambda k: read_one_var_func(k[0], k[1]))
                    else:
                        _ov_ad = _load(step, read_func)
                    if overlay_var in _ov_ad.variables:
                        _ov_raw = _ov_ad.variables[overlay_var]
                        _ov_traces = _compute_contour_traces(
                            _ov_raw, _ov_ad.x_cc, _ov_ad.y_cc,
                            int(overlay_nlevels or 5),
                            overlay_color or 'white',
                            float(overlay_lw or 1.0),
                        )
                        for _ct in _ov_traces:
                            fig.add_trace(_ct)
                except (OSError, ValueError, EOFError):
                    pass  # contour overlay is best-effort

        else:                            # 1D
            plot_y = _tf(raw) if log else raw
            fig.add_trace(go.Scatter(
                x=ad.x_cc, y=plot_y, mode='lines',
                line=dict(color=_ACCENT, width=2), name=selected_var,
            ))
            fig.update_layout(
                xaxis=dict(title='x', color=_TEXT, gridcolor=_OVER, exponentformat='e'),
                yaxis=dict(title=cbar_title, color=_TEXT, gridcolor=_OVER,
                           tickformat='.2e',
                           range=[cmin, cmax] if (vmin_in is not None or vmax_in is not None) else None),
                plot_bgcolor=_BG,
            )
            title = f'{selected_var}  ·  step {step}'

        _t_trace = time.perf_counter()
        fig.update_layout(
            title=dict(text=title, font=dict(color=_TEXT, size=13, family='monospace')),
            paper_bgcolor=_BG,
            font=dict(color=_TEXT, family='monospace'),
            margin=dict(l=0, r=130, t=36, b=0),
            uirevision=mode,   # preserve camera angle within a mode
        )
        _t_layout = time.perf_counter()

        dmin, dmax = float(np.nanmin(_raw_range)), float(np.nanmax(_raw_range))
        cons.print(
            f'[dim]viz timing  step={step}  shape={raw.shape}'
            f'  load={_t_load-_t0:.3f}s'
            f'  prep={_t_prep-_t_load:.3f}s'
            f'  trace={_t_trace-_t_prep:.3f}s'
            f'  layout={_t_layout-_t_trace:.3f}s'
            f'  total={_t_layout-_t0:.3f}s[/dim]'
        )
        status = html.Div([
            html.Span(f'step {step}', style={'color': _YELLOW}),
            html.Span(f'  ·  shape {raw.shape}', style={'color': _MUTED}),
            html.Br(),
            html.Span('min ', style={'color': _MUTED}),
            html.Span(f'{dmin:.4g}', style={'color': _BLUE}),
            html.Span('  max ', style={'color': _MUTED}),
            html.Span(f'{dmax:.4g}', style={'color': _RED}),
        ])


        _now = time.perf_counter()
        _last_update_t[0] = _now
        # Full render is heavier browser-side than a patch — add extra headroom
        # so the browser has time to rebuild the WebGL scene before the next
        # frame is sent.  For 3D with overlays this can take 1–2s.
        _server_s = _now - _t0
        _min_frame_gap[0] = max(0.8, min(3.0, _server_s + 0.7))
        return fig, status, no_update, _PV_HIDE, _GRAPH_SHOW

    # ------------------------------------------------------------------
    cons.print(f'\n[bold green]Interactive viz server:[/bold green] '
               f'[bold]http://{host}:{port}[/bold]')
    if host in ('127.0.0.1', 'localhost'):
        cons.print(
            f'\n[dim]To view from your laptop/desktop, open a [bold]new terminal on your local machine[/bold] and run:[/dim]\n'
            f'\n  [bold]ssh -L {port}:localhost:{port} user@cluster-hostname[/bold]\n'
            f'\n[dim]  Replace [bold]user@cluster-hostname[/bold] with whatever you normally pass to [bold]ssh[/bold][/dim]\n'
            f'[dim]  to reach this cluster (e.g. [bold]jdoe@login.delta.ncsa.illinois.edu[/bold] or an[/dim]\n'
            f'[dim]  alias from your [bold]~/.ssh/config[/bold]).[/dim]\n'
            f'[dim]  Then open [bold]http://localhost:{port}[/bold] in your local browser.[/dim]\n'
            f'[dim]  If you see [bold]Address already in use[/bold], free the port with:[/dim]\n'
            f'  [bold]lsof -ti :{port} | xargs kill[/bold]'
        )
    if ndim == 3 and not _pv_ok[0] and not _PV_AVAILABLE:
        cons.print(
            '[dim][yellow]Note:[/yellow] PyVista server-side rendering '
            'is not available on macOS (VTK Cocoa threading limitation). '
            '3D playback uses the Plotly WebGL path instead.[/dim]'
        )

    cons.print('[dim]\nCtrl+C to stop.[/dim]\n')
    app.run(debug=False, port=port, host=host)
