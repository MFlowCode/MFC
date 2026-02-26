"""
Terminal UI (TUI) for MFC visualization using Textual + plotext.

Launched via ``./mfc.sh viz <dir> --tui [--var VAR] [--step STEP]``.
Opens a full-terminal interactive viewer that works over SSH with no
browser or port-forwarding required.

Supports 1D line plots and 2D heatmaps only.

Requires: textual, textual-plotext, plotext
"""
from __future__ import annotations

from typing import Callable, Dict, List, Optional

import numpy as np

from rich.color import Color as RichColor
from rich.console import Group as RichGroup
from rich.style import Style
from rich.text import Text as RichText

from textual import on
from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import Horizontal, Vertical
from textual.reactive import reactive
from textual.widgets import Footer, Header, Label, ListItem, ListView, Static

from textual_plotext import PlotextPlot

from mfc.printer import cons

# Colormaps available via [c] cycling
_CMAPS: List[str] = [
    'viridis', 'plasma', 'inferno', 'magma', 'cividis',
    'coolwarm', 'RdBu_r', 'seismic', 'gray',
]

# ---------------------------------------------------------------------------
# Step cache  {step -> AssembledData}
# ---------------------------------------------------------------------------
_cache: Dict[int, object] = {}


def _load(step: int, read_func: Callable) -> object:
    if step not in _cache:
        _cache[step] = read_func(step)
    return _cache[step]


# ---------------------------------------------------------------------------
# Plot widget
# ---------------------------------------------------------------------------

class MFCPlot(PlotextPlot):  # pylint: disable=too-many-instance-attributes
    """Plotext plot widget.  Caller sets ._x_cc / ._y_cc / ._data / ._ndim /
    ._varname / ._step before calling .refresh()."""

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
        h_plot = max(self.size.height - 4, 4)  # -2 border, -2 header+footer

        # Right side: gap + gradient strip + value labels
        _CB_GAP, _CB_W, _CB_LBL = 1, 2, 9
        w_map = max(w_plot - _CB_GAP - _CB_W - _CB_LBL, 4)

        ix = np.linspace(0, data.shape[0] - 1, w_map, dtype=int)
        iy = np.linspace(0, data.shape[1] - 1, h_plot, dtype=int)
        ds = data[np.ix_(ix, iy)]  # pylint: disable=unsubscriptable-object

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
            # Heatmap cells — one terminal character per data point.
            for col in range(w_map):
                r = int(rgba[row, col, 0] * 255)
                g = int(rgba[row, col, 1] * 255)
                b = int(rgba[row, col, 2] * 255)
                line.append(" ", style=Style(bgcolor=RichColor.from_rgb(r, g, b)))
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
                mid = np.sqrt(vmin * vmax) if (self._log_scale and vmin > 0) else (vmin + vmax) / 2
                lbl = f" {mid:.3g}"
            else:
                lbl = ""
            line.append(lbl.ljust(_CB_LBL)[:_CB_LBL])
            lines.append(line)

        y_cc = self._y_cc if self._y_cc is not None else np.array([0.0, 1.0])
        log_tag = "  [log]" if log_active else ("  [log n/a]" if self._log_scale else "")
        frozen_tag = "  [frozen]" if self._vmin is not None else ""
        header = RichText(
            f" {self._varname}  (step {self._step})"
            f"   [{vmin:.3g}, {vmax:.3g}]{log_tag}{frozen_tag}",
            style="bold"
        )
        footer = RichText(
            f" x: [{x_cc[0]:.3f} \u2026 {x_cc[-1]:.3f}]"  # pylint: disable=unsubscriptable-object
            f"   y: [{y_cc[0]:.3f} \u2026 {y_cc[-1]:.3f}]",
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

    #var-title {
        text-style: bold;
        color: $accent;
        padding: 0 0 1 0;
    }

    #var-list {
        height: 1fr;
    }

    #status {
        dock: bottom;
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
        **kwargs,
    ):
        super().__init__(**kwargs)
        self._steps = steps
        self._varnames = varnames
        self._read_func = read_func
        self._ndim = ndim
        # Store init_var but don't set the reactive yet — the DOM doesn't exist
        # until on_mount, and the watcher calls query_one which needs the DOM.
        self._init_var = init_var or (varnames[0] if varnames else "")
        self._frozen_range: Optional[tuple] = None
        self._play_timer = None

    def compose(self) -> ComposeResult:
        yield Header(show_clock=False)
        with Horizontal(id="content"):
            with Vertical(id="sidebar"):
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
        # Highlight the initial variable in the sidebar list
        lv = self.query_one("#var-list", ListView)
        for i, v in enumerate(self._varnames):
            if v == self.var_name:
                lv.index = i
                break

    # ------------------------------------------------------------------
    # Reactive watchers
    # ------------------------------------------------------------------

    def watch_step_idx(self, _old: int, _new: int) -> None:
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
        )

    def _push_data(self) -> None:
        """Load the current step/var and push data into the plot widget."""
        if not self._steps or not self.var_name:
            return
        step = self._steps[self.step_idx]
        try:
            assembled = _load(step, self._read_func)
        except Exception as exc:  # pylint: disable=broad-except
            self.query_one("#status", Static).update(
                f" [red]Error loading step {step}: {exc}[/red]"
            )
            return

        data = assembled.variables.get(self.var_name)
        plot = self.query_one("#plot", MFCPlot)
        plot._x_cc = assembled.x_cc        # pylint: disable=protected-access
        plot._y_cc = assembled.y_cc        # pylint: disable=protected-access
        plot._data = data                  # pylint: disable=protected-access
        plot._ndim = self._ndim            # pylint: disable=protected-access
        plot._varname = self.var_name      # pylint: disable=protected-access
        plot._step = step                  # pylint: disable=protected-access
        plot._cmap_name = self.cmap_name   # pylint: disable=protected-access
        plot._log_scale = self.log_scale   # pylint: disable=protected-access
        if self._frozen_range is not None:
            plot._vmin, plot._vmax = self._frozen_range  # pylint: disable=protected-access
        else:
            plot._vmin = None              # pylint: disable=protected-access
            plot._vmax = None              # pylint: disable=protected-access
        plot.refresh()

        self.query_one("#status", Static).update(self._status_text())

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
) -> None:
    """Launch the Textual TUI for MFC visualization (1D/2D only)."""
    if ndim not in (1, 2):
        raise ValueError(
            f"--tui only supports 1D and 2D data (got ndim={ndim}). "
            "Use --interactive for 3D data."
        )

    # Preload first step to discover variables
    first = _load(steps[0], read_func)
    varnames = sorted(first.variables.keys())
    if not varnames:
        raise ValueError("No variables found in data")
    if init_var not in varnames:
        init_var = varnames[0]

    cons.print(
        f"[bold]Launching TUI[/bold] — {len(steps)} step(s), "
        f"{len(varnames)} variable(s)"
    )
    cons.print("[dim]  ,/. or ←/→  prev/next step  •  space play  •  l log  •  f freeze  •  ↑↓ variable  •  q quit[/dim]")

    _cache.clear()
    _cache[steps[0]] = first

    app = MFCTuiApp(
        steps=steps,
        varnames=varnames,
        read_func=read_func,
        ndim=ndim,
        init_var=init_var,
    )
    app.run()
