"""
Terminal UI (TUI) for MFC visualization using Textual + plotext.

Launched via ``./mfc.sh viz <dir> --tui [--var VAR] [--step STEP]``.
Opens a full-terminal interactive viewer that works over SSH with no
browser or port-forwarding required.

Supports 1D line plots and 2D heatmaps only.

Requires: textual>=0.43, textual-plotext, plotext
"""
from __future__ import annotations

from typing import Callable, List, Optional, Tuple

import numpy as np

from rich.color import Color as RichColor
from rich.console import Group as RichGroup
from rich.style import Style
from rich.text import Text as RichText

from textual import on, work
from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import Horizontal, Vertical
from textual.message import Message
from textual.reactive import reactive
from textual.widgets import (
    Digits, Footer, Header, Label, ListItem, ListView, Static,
)
from textual.worker import get_current_worker

from textual_plotext import PlotextPlot

from mfc.common import MFCException
from mfc.printer import cons
from . import _step_cache

# Colormaps available via [c] cycling
_CMAPS: List[str] = [
    'viridis', 'plasma', 'inferno', 'magma', 'cividis',
    'coolwarm', 'RdBu_r', 'seismic', 'gray',
]

_load = _step_cache.load
_CACHE_MAX = _step_cache.CACHE_MAX

# Terminal character cells are approximately twice as tall as they are wide in
# pixels (e.g. 8 px wide × 16 px tall).  A square physical domain should
# therefore occupy a ~2:1 (col:row) character grid to look correct.
_CELL_RATIO: float = 2.0

# Physical domain aspect ratio is clamped to this range so that very elongated
# domains don't produce unusable slivers.
_ASPECT_MIN: float = 0.2
_ASPECT_MAX: float = 5.0

# Rows of header above the heatmap inside the content area.
_HEADER_ROWS: int = 1


# ---------------------------------------------------------------------------
# Plot widget
# ---------------------------------------------------------------------------

class MFCPlot(PlotextPlot):  # pylint: disable=too-many-instance-attributes,too-few-public-methods
    """Plotext plot widget.  Caller sets ._x_cc / ._y_cc / ._data / ._ndim /
    ._varname / ._step before calling .refresh()."""

    # Disable text-selection so Textual doesn't intercept left-button events
    # before they bubble to our on_mouse_up handler.
    ALLOW_SELECT = False

    class Clicked(Message):
        """Posted when the user clicks a heatmap cell (Feature 5)."""
        def __init__(self, x_val: float, y_val: float, val: float) -> None:
            self.x_val = x_val
            self.y_val = y_val
            self.val = val
            super().__init__()

    DEFAULT_CSS = """
    MFCPlot {
        border: solid $accent;
        width: 1fr;
        height: 1fr;
    }
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._x_cc: Optional[np.ndarray] = None
        self._y_cc: Optional[np.ndarray] = None
        self._data: Optional[np.ndarray] = None
        self._ndim: int = 1
        self._varname: str = ""
        self._step: int = 0
        self._cmap_name: str = _CMAPS[0]
        self._log_scale: bool = False
        self._vmin: Optional[float] = None
        self._vmax: Optional[float] = None
        self._last_vmin: float = 0.0
        self._last_vmax: float = 1.0
        self._bubbles: Optional[np.ndarray] = None  # (N,4) x,y,z,r
        # Zoom state: (x_frac0, x_frac1, y_frac0, y_frac1) all in [0,1].
        self._zoom: Tuple[float, float, float, float] = (0.0, 1.0, 0.0, 1.0)
        # Last-render dimensions — used to map click/scroll coords back to data.
        self._last_w_map: int = 0
        self._last_h_plot: int = 0
        self._last_ix: Optional[np.ndarray] = None
        self._last_iy: Optional[np.ndarray] = None

    # ------------------------------------------------------------------
    # Zoom helpers (Feature 6)
    # ------------------------------------------------------------------

    def reset_zoom(self) -> None:
        """Reset to full view."""
        self._zoom = (0.0, 1.0, 0.0, 1.0)
        self.refresh()

    def _zoom_around(  # pylint: disable=too-many-locals
        self, cx_frac: float, cy_frac: float, factor: float
    ) -> None:
        """Zoom by *factor* centred at *(cx_frac, cy_frac)* in [0,1]² of current view."""
        x0, x1, y0, y1 = self._zoom
        x_span = x1 - x0
        y_span = y1 - y0
        new_x_span = x_span * factor
        new_y_span = y_span * factor
        # Enforce minimum zoom (2 % of full range per axis).
        if new_x_span < 0.02 or new_y_span < 0.02:
            return
        cx = x0 + cx_frac * x_span
        cy = y0 + cy_frac * y_span
        new_x0 = max(0.0, cx - cx_frac * new_x_span)
        new_x1 = min(1.0, cx + (1.0 - cx_frac) * new_x_span)
        new_y0 = max(0.0, cy - cy_frac * new_y_span)
        new_y1 = min(1.0, cy + (1.0 - cy_frac) * new_y_span)
        if new_x1 - new_x0 < 0.02 or new_y1 - new_y0 < 0.02:
            return
        self._zoom = (new_x0, new_x1, new_y0, new_y1)

    # ------------------------------------------------------------------
    # Mouse scroll handlers — zoom (Feature 6)
    # ------------------------------------------------------------------

    def _scroll_zoom(self, event, factor: float) -> None:
        """Zoom by *factor* centred on the scroll cursor position."""
        if self._data is None or self._ndim != 2:
            return
        # get_content_offset_capture never returns None — safe at border too.
        offset = event.get_content_offset_capture(self)
        col = offset.x
        row = offset.y - _HEADER_ROWS
        w = max(self._last_w_map, 1)
        h = max(self._last_h_plot, 1)
        cx_frac = max(0.0, min(1.0, col / (w - 1) if w > 1 else 0.5))
        # Row 0 = top = y_max → cy_frac = 1 in zoom space; row h-1 = 0.
        cy_frac = max(0.0, min(1.0, 1.0 - row / (h - 1) if h > 1 else 0.5))
        self._zoom_around(cx_frac, cy_frac, factor=factor)
        event.stop()
        self.refresh()

    def on_mouse_scroll_up(self, event) -> None:  # type: ignore[override]
        self._scroll_zoom(event, factor=0.75)

    def on_mouse_scroll_down(self, event) -> None:  # type: ignore[override]
        self._scroll_zoom(event, factor=1.0 / 0.75)

    def on_mouse_up(self, event) -> None:  # pylint: disable=too-many-locals
        """Feature 5 — post Clicked message with the data value at the heatmap cell."""
        if event.button != 1:
            return
        if self._data is None or self._ndim != 2:
            return
        if self._last_w_map == 0 or self._last_ix is None or self._last_iy is None:
            return
        # Offset relative to this widget's content area (inside the border).
        offset = event.get_content_offset_capture(self)
        col = offset.x
        row = offset.y - _HEADER_ROWS
        col = max(0, min(col, self._last_w_map - 1))
        row = max(0, min(row, self._last_h_plot - 1))
        n_ix = len(self._last_ix)
        n_iy = len(self._last_iy)
        ix_pos = int(np.round(col * (n_ix - 1) / max(self._last_w_map - 1, 1)))
        # Display is y-flipped: row 0 = top = y_max.
        iy_pos = n_iy - 1 - int(np.round(row * (n_iy - 1) / max(self._last_h_plot - 1, 1)))
        xi = int(self._last_ix[np.clip(ix_pos, 0, n_ix - 1)])  # pylint: disable=unsubscriptable-object
        yi = int(self._last_iy[np.clip(iy_pos, 0, n_iy - 1)])  # pylint: disable=unsubscriptable-object
        x_cc = self._x_cc
        y_cc = self._y_cc if self._y_cc is not None else np.array([0.0, 1.0])
        data = self._data
        x_val = float(x_cc[xi])   # type: ignore[index]  # pylint: disable=unsubscriptable-object
        y_val = float(y_cc[yi])
        val = float(data[xi, yi])  # type: ignore[index]  # pylint: disable=unsubscriptable-object
        self.post_message(MFCPlot.Clicked(x_val, y_val, val))

    def render(self):  # pylint: disable=too-many-branches,too-many-locals,too-many-statements
        data = self._data
        x_cc = self._x_cc
        self.plt.clear_figure()

        # 1D: use normal plotext path — gives proper axes and title for free.
        if data is None or x_cc is None or self._ndim == 1:
            if data is not None and x_cc is not None:
                if self._log_scale:
                    plot_y = np.where(data > 0, np.log10(np.maximum(data, 1e-300)), np.nan)
                    ylabel = f"log\u2081\u2080({self._varname})"
                    title_tag = "  [log]"
                else:
                    plot_y = data
                    ylabel = self._varname
                    title_tag = ""
                if self._vmin is not None or self._vmax is not None:
                    title_tag += "  [frozen]"
                finite = plot_y[np.isfinite(plot_y)]
                self._last_vmin = float(finite.min()) if finite.size else 0.0
                self._last_vmax = float(finite.max()) if finite.size else 1.0
                self.plt.plot(x_cc.tolist(), plot_y.tolist())
                self.plt.xlabel("x")
                self.plt.ylabel(ylabel)
                self.plt.title(f"{self._varname}  (step {self._step}){title_tag}")
                if self._vmin is not None or self._vmax is not None:
                    lo = self._vmin if self._vmin is not None else self._last_vmin
                    hi = self._vmax if self._vmax is not None else self._last_vmax
                    self.plt.ylim(lo, hi)
            else:
                self.plt.title("No data loaded")
            return super().render()

        # 2D: pure-Rich heatmap with vertical colorbar.
        import matplotlib                     # pylint: disable=import-outside-toplevel
        import matplotlib.colors as mcolors  # pylint: disable=import-outside-toplevel

        # Content area = widget size minus 1-char border on each side.
        # Reserve 1 row each for header and footer → h_plot rows for the image.
        w_plot = max(self.size.width - 2, 4)
        h_plot_avail = max(self.size.height - 4, 4)  # -2 border, -2 header+footer

        # Right side: gap + gradient strip + value labels.
        # 11 chars fits negative scientific notation e.g. " -2.26e-05".
        _CB_GAP, _CB_W, _CB_LBL = 1, 2, 11
        w_map_avail = max(w_plot - _CB_GAP - _CB_W - _CB_LBL, 4)

        # Preserve the physical x/y aspect ratio.
        y_cc_2d = self._y_cc if self._y_cc is not None else np.array([0.0, 1.0])
        x_extent = max(abs(float(x_cc[-1]) - float(x_cc[0])), 1e-30)  # pylint: disable=unsubscriptable-object
        y_extent = max(abs(float(y_cc_2d[-1]) - float(y_cc_2d[0])), 1e-30)
        domain_ratio = float(np.clip(x_extent / y_extent, _ASPECT_MIN, _ASPECT_MAX))
        char_ratio = domain_ratio * _CELL_RATIO

        w_ideal = int(round(h_plot_avail * char_ratio))
        if w_ideal <= w_map_avail:
            w_map  = max(w_ideal, 4)
            h_plot = h_plot_avail
        else:
            h_plot = max(int(round(w_map_avail / char_ratio)), 4)
            w_map  = w_map_avail

        # Apply zoom window to data index ranges (Feature 6).
        x0_f, x1_f, y0_f, y1_f = self._zoom
        x0_i = int(x0_f * (data.shape[0] - 1))
        x1_i = max(x0_i + 1, int(x1_f * (data.shape[0] - 1)))
        y0_i = int(y0_f * (data.shape[1] - 1))
        y1_i = max(y0_i + 1, int(y1_f * (data.shape[1] - 1)))
        ix = np.linspace(x0_i, x1_i, w_map, dtype=int)
        iy = np.linspace(y0_i, y1_i, h_plot, dtype=int)

        # Cache for click/scroll coordinate mapping (Features 5 & 6).
        self._last_w_map = w_map
        self._last_h_plot = h_plot
        self._last_ix = ix
        self._last_iy = iy

        ds = data[np.ix_(ix, iy)]  # pylint: disable=unsubscriptable-object

        # Compute which screen cells to stamp with an open-circle glyph.
        bubble_cells: set = set()
        bubbles = self._bubbles
        if bubbles is not None and len(bubbles) > 0:
            x_phys = x_cc[ix]   # type: ignore[index]  # pylint: disable=unsubscriptable-object
            y_phys = y_cc_2d[iy]
            x_min, x_max = float(x_phys[0]),  float(x_phys[-1])
            y_min, y_max = float(y_phys[0]),  float(y_phys[-1])
            x_range = max(abs(x_max - x_min), 1e-30)
            y_range = max(abs(y_max - y_min), 1e-30)
            for b in bubbles:  # pylint: disable=not-an-iterable
                bx, by, br = float(b[0]), float(b[1]), float(b[3])
                if bx < x_min - br or bx > x_max + br:
                    continue
                if by < y_min - br or by > y_max + br:
                    continue
                col_c = (bx - x_min) / x_range * (w_map - 1)
                row_c = (y_max - by) / y_range * (h_plot - 1)
                col_r = br / x_range * (w_map - 1)
                row_r = br / y_range * (h_plot - 1)
                if col_r < 0.5 and row_r < 0.5:
                    c, r = int(round(col_c)), int(round(row_c))
                    if 0 <= r < h_plot and 0 <= c < w_map:
                        bubble_cells.add((r, c))
                else:
                    n_pts = min(max(12, int(2 * np.pi * max(col_r, row_r))), 72)
                    for ti in range(n_pts):
                        angle = 2 * np.pi * ti / n_pts
                        c = int(round(col_c + col_r * np.cos(angle)))
                        r = int(round(row_c + row_r * np.sin(angle)))
                        if 0 <= r < h_plot and 0 <= c < w_map:
                            bubble_cells.add((r, c))

        vmin = self._vmin if self._vmin is not None else float(ds.min())
        vmax = self._vmax if self._vmax is not None else float(ds.max())
        if vmax <= vmin:
            vmax = vmin + 1e-10
        cmap = matplotlib.colormaps[self._cmap_name]
        log_active = False
        if self._log_scale:
            pos = ds[ds > 0]
            lo = float(np.nanmin(pos)) if pos.size > 0 else None
            if lo is not None and lo < vmax:
                norm = mcolors.LogNorm(vmin=lo, vmax=vmax)
                vmin = lo
                log_active = True
            else:
                norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
        else:
            norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
        self._last_vmin = vmin
        self._last_vmax = vmax
        # Transpose + flip so y=0 appears at the bottom of the display.
        rgba = cmap(norm(ds.T[::-1]))  # (h_plot, w_map, 4)

        lines = []
        for row in range(h_plot):
            line = RichText()
            for col in range(w_map):
                r = int(rgba[row, col, 0] * 255)
                g = int(rgba[row, col, 1] * 255)
                b = int(rgba[row, col, 2] * 255)
                bg = RichColor.from_rgb(r, g, b)
                if (row, col) in bubble_cells:
                    line.append("○", style=Style(bgcolor=bg, color="white", bold=True))
                else:
                    line.append(" ", style=Style(bgcolor=bg))
            # Gap
            line.append(" " * _CB_GAP)
            # Colorbar gradient strip (t=1 at top = vmax, t=0 at bottom = vmin)
            t = 1.0 - row / max(h_plot - 1, 1)
            cb = cmap(t)
            cr, cg, cbb = int(cb[0] * 255), int(cb[1] * 255), int(cb[2] * 255)
            for _ in range(_CB_W):
                line.append(" ", style=Style(bgcolor=RichColor.from_rgb(cr, cg, cbb)))
            # Value labels at top, middle, bottom
            if row == 0:
                lbl = f" {vmax:.3g}"
            elif row == h_plot - 1:
                lbl = f" {vmin:.3g}"
            elif row == h_plot // 2:
                mid = np.sqrt(vmin * vmax) if (log_active and vmin > 0) else (vmin + vmax) / 2
                lbl = f" {mid:.3g}"
            else:
                lbl = ""
            line.append(lbl.ljust(_CB_LBL)[:_CB_LBL])
            lines.append(line)

        log_tag = "  [log]" if log_active else ("  [log n/a]" if self._log_scale else "")
        frozen_tag = "  [frozen]" if self._vmin is not None else ""
        zoomed_tag = "  [zoom]" if self._zoom != (0.0, 1.0, 0.0, 1.0) else ""
        header = RichText(
            f" {self._varname}  (step {self._step})"
            f"   [{vmin:.3g}, {vmax:.3g}]{log_tag}{frozen_tag}{zoomed_tag}",
            style="bold"
        )
        # Show the visible coordinate range (reflects zoom when active).
        x_lo = float(x_cc[ix[0]])   # type: ignore[index]  # pylint: disable=unsubscriptable-object
        x_hi = float(x_cc[ix[-1]])  # type: ignore[index]  # pylint: disable=unsubscriptable-object
        y_vis = y_cc_2d[iy]
        footer = RichText(
            f" x: [{x_lo:.3f} \u2026 {x_hi:.3f}]"
            f"   y: [{float(y_vis[0]):.3f} \u2026 {float(y_vis[-1]):.3f}]",
            style="dim"
        )
        return RichGroup(header, *lines, footer)


# ---------------------------------------------------------------------------
# Main TUI app
# ---------------------------------------------------------------------------

class MFCTuiApp(App):  # pylint: disable=too-many-instance-attributes
    """Textual TUI for MFC post-processed data."""

    CSS = """
    Screen {
        layers: base;
    }

    #content {
        height: 1fr;
        layout: horizontal;
    }

    #sidebar {
        width: 22;
        border-right: solid $accent;
        padding: 0 1;
    }

    #step-counter {
        width: 1fr;
        height: 4;
        color: $accent;
    }

    #var-title {
        text-style: bold;
        color: $accent;
        padding: 0 0 1 0;
    }

    #var-list {
        height: 1fr;
    }

    #status {
        height: 1;
        background: $panel;
        color: $text-muted;
        padding: 0 1;
    }
    """

    BINDINGS = [
        Binding("q", "quit", "Quit"),
        Binding("comma", "prev_step", "◀ step"),
        Binding("period", "next_step", "step ▶"),
        Binding("left", "prev_step", "◀ step", show=False),
        Binding("right", "next_step", "step ▶", show=False),
        Binding("space", "toggle_play", "▶/⏸"),
        Binding("c", "cycle_cmap", "cmap"),
        Binding("l", "toggle_log", "log"),
        Binding("f", "toggle_freeze", "freeze"),
        Binding("r", "reset_zoom", "zoom↺", show=False),
    ]

    step_idx: reactive[int] = reactive(0, always_update=True)
    var_name: reactive[str] = reactive("", always_update=True)
    cmap_name: reactive[str] = reactive(_CMAPS[0], always_update=True)
    log_scale: reactive[bool] = reactive(False, always_update=True)
    playing: reactive[bool] = reactive(False, always_update=True)

    def __init__(  # pylint: disable=too-many-arguments,too-many-positional-arguments
        self,
        steps: List[int],
        varnames: List[str],
        read_func: Callable,
        ndim: int,
        init_var: Optional[str] = None,
        bubble_func: Optional[Callable] = None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self._steps = steps
        self._varnames = varnames
        self._read_func = read_func
        self._ndim = ndim
        self._bubble_func = bubble_func
        self._init_var = init_var or (varnames[0] if varnames else "")
        self._frozen_range: Optional[Tuple[float, float]] = None
        self._play_timer = None
        self._click_info: str = ""  # last clicked cell; included in status bar

    def compose(self) -> ComposeResult:
        yield Header(show_clock=False)
        with Horizontal(id="content"):
            with Vertical(id="sidebar"):
                yield Digits("0", id="step-counter")
                yield Label("Variables", id="var-title")
                yield ListView(
                    *[ListItem(Label(v), id=f"var-{v}") for v in self._varnames],
                    id="var-list",
                )
            yield MFCPlot(id="plot")
        yield Static(self._status_text(), id="status")
        yield Footer()

    def on_mount(self) -> None:
        # DOM is ready — now safe to set the reactive (fires watcher → _push_data)
        self.var_name = self._init_var
        lv = self.query_one("#var-list", ListView)
        for i, v in enumerate(self._varnames):
            if v == self.var_name:
                lv.index = i
                break

    # ------------------------------------------------------------------
    # Reactive watchers
    # ------------------------------------------------------------------

    def watch_step_idx(self, _old: int, _new: int) -> None:
        self._click_info = ""
        self._push_data()

    def watch_var_name(self, _old: str, _new: str) -> None:
        self._push_data()

    def watch_cmap_name(self, _old: str, _new: str) -> None:
        self._push_data()

    def watch_log_scale(self, _old: bool, _new: bool) -> None:
        self._push_data()

    def watch_playing(self, _old: bool, new: bool) -> None:
        if new:
            self._play_timer = self.set_interval(0.5, self._auto_advance)
        else:
            if self._play_timer is not None:
                self._play_timer.stop()
                self._play_timer = None

    # ------------------------------------------------------------------
    # MFCPlot.Clicked handler — update status bar (Feature 5)
    # ------------------------------------------------------------------

    def on_mfcplot_clicked(self, event: MFCPlot.Clicked) -> None:
        """Receive the heatmap click message and update the status bar."""
        self._click_info = (
            f"  │  x={event.x_val:.4f}  y={event.y_val:.4f}  val={event.val:.6g}"
        )
        self.query_one("#status", Static).update(self._status_text())

    # ------------------------------------------------------------------
    # Background data loading (Feature 4)
    # ------------------------------------------------------------------

    @work(exclusive=True, thread=True)
    def _push_data(self) -> None:
        """Load the current step/var in a background thread and push to the plot."""
        if not self._steps or not self.var_name:
            return
        # Snapshot all reactive state before entering the thread.
        step_idx = min(self.step_idx, len(self._steps) - 1)
        step = self._steps[step_idx]
        var = self.var_name
        cmap = self.cmap_name
        log = self.log_scale
        frozen = self._frozen_range

        try:
            assembled = _load(step, self._read_func)
        except (OSError, ValueError, EOFError) as exc:
            self.call_from_thread(
                self.query_one("#status", Static).update,
                f" [red]Error loading step {step}: {exc}[/red]",
            )
            return

        worker = get_current_worker()
        if worker.is_cancelled:
            return

        data = assembled.variables.get(var)
        bubbles = None
        if self._bubble_func is not None and self._ndim == 2:
            try:
                bubbles = self._bubble_func(step)
            except (OSError, ValueError):
                pass  # bubble overlay is best-effort; skip on read errors

        self.call_from_thread(
            self._apply_data, assembled, data, step, var, cmap, log, frozen, bubbles,
        )

    def _apply_data(  # pylint: disable=too-many-arguments,too-many-positional-arguments
        self,
        assembled,
        data: Optional[np.ndarray],
        step: int,
        var: str,
        cmap: str,
        log: bool,
        frozen: Optional[Tuple[float, float]],
        bubbles: Optional[np.ndarray],
    ) -> None:
        """Apply loaded data to the plot widget.  Runs on the main thread."""
        plot = self.query_one("#plot", MFCPlot)
        plot._x_cc = assembled.x_cc        # pylint: disable=protected-access
        plot._y_cc = assembled.y_cc        # pylint: disable=protected-access
        plot._data = data                  # pylint: disable=protected-access
        plot._ndim = self._ndim            # pylint: disable=protected-access
        plot._varname = var                # pylint: disable=protected-access
        plot._step = step                  # pylint: disable=protected-access
        plot._cmap_name = cmap             # pylint: disable=protected-access
        plot._log_scale = log              # pylint: disable=protected-access
        plot._bubbles = bubbles            # pylint: disable=protected-access
        if frozen is not None:
            plot._vmin, plot._vmax = frozen  # pylint: disable=protected-access
        else:
            plot._vmin = None              # pylint: disable=protected-access
            plot._vmax = None              # pylint: disable=protected-access
        plot.refresh()

        # Update step counter (Feature 2).
        self.query_one("#step-counter", Digits).update(str(step))

        self.query_one("#status", Static).update(self._status_text())

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _status_text(self) -> str:
        step = self._steps[self.step_idx] if self._steps else 0
        total = len(self._steps)
        flags = []
        if self.log_scale:
            flags.append("log")
        if self._frozen_range is not None:
            flags.append("frozen")
        if self.playing:
            flags.append("▶")
        flag_str = ("  " + "  ".join(flags)) if flags else ""
        return (
            f" step {step}  [{self.step_idx + 1}/{total}]"
            f"  var: {self.var_name}"
            f"  cmap: {self.cmap_name}"
            f"{flag_str}"
            f"{self._click_info}"
        )

    # ------------------------------------------------------------------
    # Actions
    # ------------------------------------------------------------------

    @on(ListView.Selected, "#var-list")
    def on_var_selected(self, event: ListView.Selected) -> None:
        item_id = event.item.id or ""
        if item_id.startswith("var-"):
            self.var_name = item_id[4:]

    def action_prev_step(self) -> None:
        if self.step_idx > 0:
            self.step_idx -= 1

    def action_next_step(self) -> None:
        if self.step_idx < len(self._steps) - 1:
            self.step_idx += 1

    def action_cycle_cmap(self) -> None:
        idx = (_CMAPS.index(self.cmap_name) + 1) % len(_CMAPS)
        self.cmap_name = _CMAPS[idx]

    def action_toggle_log(self) -> None:
        self.log_scale = not self.log_scale

    def action_toggle_freeze(self) -> None:
        if self._frozen_range is not None:
            self._frozen_range = None
        else:
            plot = self.query_one("#plot", MFCPlot)
            self._frozen_range = (plot._last_vmin, plot._last_vmax)  # pylint: disable=protected-access
        self._push_data()

    def action_toggle_play(self) -> None:
        self.playing = not self.playing

    def action_reset_zoom(self) -> None:
        """Feature 6 — reset 2D zoom to full view."""
        self.query_one("#plot", MFCPlot).reset_zoom()

    def _auto_advance(self) -> None:
        self.step_idx = (self.step_idx + 1) % len(self._steps)


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def run_tui(
    init_var: Optional[str],
    steps: List[int],
    read_func: Callable,
    ndim: int,
    bubble_func: Optional[Callable] = None,
) -> None:
    """Launch the Textual TUI for MFC visualization (1D/2D only)."""
    if ndim not in (1, 2):
        raise MFCException(
            f"--tui only supports 1D and 2D data (got ndim={ndim}). "
            "Use --interactive for 3D data."
        )

    # Preload first step to discover variables
    first = _load(steps[0], read_func)
    varnames = sorted(first.variables.keys())
    if not varnames:
        raise MFCException("No variables found in data")
    if init_var not in varnames:
        init_var = varnames[0]

    cons.print(
        f"[bold]Launching TUI[/bold] — {len(steps)} step(s), "
        f"{len(varnames)} variable(s)"
    )
    cons.print(
        "[dim]  ,/. or ←/→  step  •  space play  •  l log  •  f freeze"
        "  •  c cmap  •  ↑↓ var  •  scroll zoom  •  r reset zoom  •  click value  •  q quit[/dim]"
    )

    _step_cache.seed(steps[0], first)

    app = MFCTuiApp(
        steps=steps,
        varnames=varnames,
        read_func=read_func,
        ndim=ndim,
        init_var=init_var,
        bubble_func=bubble_func,
    )
    app.run()
