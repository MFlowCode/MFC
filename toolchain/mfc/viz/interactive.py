"""
Interactive Dash-based visualization for MFC post-processed data.

Launched via ``./mfc.sh viz <dir> --var <var> --step all --interactive``.
Opens a dark-themed web UI in your browser (or via SSH tunnel) with live
controls for slice position, isosurface thresholds, volume opacity,
colormap, log scale, vmin/vmax, and timestep playback.
"""
# pylint: disable=use-dict-literal

from typing import List, Callable, Optional

import numpy as np
import plotly.graph_objects as go
from dash import Dash, dcc, html, Input, Output, State, callback_context, no_update

from mfc.printer import cons
from . import _step_cache

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
        tooltip={'placement': 'bottom', 'always_visible': True},
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

def _build_3d(ad, raw, varname, step, mode, cmap,  # pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals,too-many-branches,too-many-statements
              log_fn, cmin, cmax, cbar_title,
              slice_axis, slice_pos,
              iso_min_frac, iso_max_frac, iso_n, iso_caps,
              vol_opacity, vol_nsurf, vol_min_frac, vol_max_frac):
    """Return (trace, title) for a 3D assembled dataset."""
    cbar = dict(
        title=dict(text=cbar_title, font=dict(color=_TEXT)),
        tickfont=dict(color=_TEXT),
        tickformat='.2e',
    )
    rng = cmax - cmin if cmax > cmin else 1.0

    if mode == 'slice':
        axis_coords = {'x': ad.x_cc, 'y': ad.y_cc, 'z': ad.z_cc}
        coords = axis_coords[slice_axis]
        coord_val = coords[0] + (coords[-1] - coords[0]) * slice_pos
        idx = int(np.clip(np.argmin(np.abs(coords - coord_val)), 0, len(coords) - 1))
        actual = float(coords[idx])

        if slice_axis == 'x':
            sliced = log_fn(raw[idx, :, :])            # (ny, nz)
            YY, ZZ = np.meshgrid(ad.y_cc, ad.z_cc, indexing='ij')
            trace = go.Surface(
                x=np.full_like(YY, actual), y=YY, z=ZZ,
                surfacecolor=sliced, cmin=cmin, cmax=cmax,
                colorscale=cmap, colorbar=cbar, showscale=True,
            )
        elif slice_axis == 'y':
            sliced = log_fn(raw[:, idx, :])            # (nx, nz)
            XX, ZZ = np.meshgrid(ad.x_cc, ad.z_cc, indexing='ij')
            trace = go.Surface(
                x=XX, y=np.full_like(XX, actual), z=ZZ,
                surfacecolor=sliced, cmin=cmin, cmax=cmax,
                colorscale=cmap, colorbar=cbar, showscale=True,
            )
        else:                                            # z
            sliced = log_fn(raw[:, :, idx])             # (nx, ny)
            XX, YY = np.meshgrid(ad.x_cc, ad.y_cc, indexing='ij')
            trace = go.Surface(
                x=XX, y=YY, z=np.full_like(XX, actual),
                surfacecolor=sliced, cmin=cmin, cmax=cmax,
                colorscale=cmap, colorbar=cbar, showscale=True,
            )
        title = f'{varname}  ·  {slice_axis} = {actual:.4g}  ·  step {step}'

    elif mode == 'isosurface':
        X3, Y3, Z3 = np.meshgrid(ad.x_cc, ad.y_cc, ad.z_cc, indexing='ij')
        vf = log_fn(raw.ravel())
        ilo = cmin + rng * iso_min_frac
        ihi = cmin + rng * max(iso_max_frac, iso_min_frac + 0.01)
        caps = dict(x_show=iso_caps, y_show=iso_caps, z_show=iso_caps)
        trace = go.Isosurface(
            x=X3.ravel(), y=Y3.ravel(), z=Z3.ravel(), value=vf,
            isomin=ilo, isomax=ihi, surface_count=int(iso_n),
            colorscale=cmap, cmin=cmin, cmax=cmax,
            caps=caps, colorbar=cbar,
        )
        title = f'{varname}  ·  {int(iso_n)} isosurfaces  ·  step {step}'

    else:                                                # volume
        X3, Y3, Z3 = np.meshgrid(ad.x_cc, ad.y_cc, ad.z_cc, indexing='ij')
        vf = log_fn(raw.ravel())
        vlo = cmin + rng * vol_min_frac
        vhi = cmin + rng * max(vol_max_frac, vol_min_frac + 0.01)
        trace = go.Volume(
            x=X3.ravel(), y=Y3.ravel(), z=Z3.ravel(), value=vf,
            isomin=vlo, isomax=vhi,
            opacity=float(vol_opacity), surface_count=int(vol_nsurf),
            colorscale=cmap, colorbar=cbar,
        )
        title = f'{varname}  ·  volume  ·  step {step}'

    return trace, title


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
#var-sel *, #step-sel *, #cmap-sel *,
#var-sel > div, #step-sel > div, #cmap-sel > div {
    background-color: %(bg)s !important;
    color: %(tx)s !important;
    border-color: %(bd)s !important;
}
#var-sel input, #step-sel input, #cmap-sel input {
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
/* Slider tooltip bubble */
.rc-slider-tooltip-inner {
    background-color: %(tx)s !important;
    color: #11111b !important;
    border: none !important;
}
.rc-slider-mark-text { color: %(tx)s !important; }
.rc-slider-rail { background-color: %(bd)s !important; }
.rc-slider-track { background-color: %(ac)s !important; }
.rc-slider-handle {
    border-color: %(ac)s !important;
    background-color: %(ac)s !important;
}
.rc-slider-dot { background-color: %(bd)s !important; border-color: %(bd)s !important; }
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

    step_opts = [{'label': str(s), 'value': s} for s in steps]
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
            dcc.Dropdown(
                id='step-sel', options=step_opts, value=steps[0], clearable=False,
                style={'fontSize': '12px', 'backgroundColor': _OVER,
                       'border': f'1px solid {_BORDER}'},
            ),
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
                labelStyle={'display': 'block', 'marginBottom': '5px', 'cursor': 'pointer'},
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
                    labelStyle={'marginRight': '14px'},
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
                dcc.Checklist(
                    id='iso-caps',
                    options=[{'label': '  Show end-caps', 'value': 'caps'}], value=[],
                    style={'fontSize': '12px', 'color': _SUB, 'marginTop': '6px'},
                ),
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
        ], style={'flex': '1', 'overflow': 'hidden', 'backgroundColor': _BG}),

        dcc.Interval(id='play-iv', interval=500, n_intervals=0, disabled=True),
        dcc.Store(id='playing-st', data=False),
    ], style={
        'display': 'flex', 'height': '100vh', 'overflow': 'hidden',
        'backgroundColor': _BG, 'fontFamily': 'monospace',
    })

    # ------------------------------------------------------------------
    # Callbacks
    # ------------------------------------------------------------------

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
            return True, iv, False, '▶  Play'
        if 'play-btn' in trig:
            playing = not is_playing
            return not playing, iv, playing, ('⏸  Pause' if playing else '▶  Play')
        return not is_playing, iv, is_playing, no_update  # fps-only change

    @app.callback(
        Output('step-sel', 'value'),
        Input('play-iv', 'n_intervals'),
        State('step-sel', 'value'),
        State('loop-chk', 'value'),
        prevent_initial_call=True,
    )
    def _advance_step(_, current, loop_val):
        try:
            idx = steps.index(current)
        except ValueError:
            idx = 0
        nxt = idx + 1
        if nxt >= len(steps):
            return steps[0] if ('loop' in (loop_val or [])) else no_update
        return steps[nxt]

    @app.callback(
        Output('ctrl-slice', 'style'),
        Output('ctrl-iso',   'style'),
        Output('ctrl-vol',   'style'),
        Input('mode-sel', 'value'),
    )
    def _toggle_controls(mode):
        show, hide = {'display': 'block'}, {'display': 'none'}
        return (
            show if mode == 'slice'      else hide,
            show if mode == 'isosurface' else hide,
            show if mode == 'volume'     else hide,
        )

    @app.callback(
        Output('vmin-inp', 'value'),
        Output('vmax-inp', 'value'),
        Input('reset-btn', 'n_clicks'),
        Input('var-sel',   'value'),
        prevent_initial_call=True,
    )
    def _reset_range(_reset, _var):
        return None, None

    @app.callback(
        Output('viz-graph',  'figure'),
        Output('status-bar', 'children'),
        Input('var-sel',     'value'),
        Input('step-sel',    'value'),
        Input('mode-sel',    'value'),
        Input('slice-axis',  'value'),
        Input('slice-pos',   'value'),
        Input('iso-min',     'value'),
        Input('iso-max',     'value'),
        Input('iso-n',       'value'),
        Input('iso-caps',    'value'),
        Input('vol-opacity', 'value'),
        Input('vol-nsurf',   'value'),
        Input('vol-min',     'value'),
        Input('vol-max',     'value'),
        Input('cmap-sel',    'value'),
        Input('log-chk',     'value'),
        Input('vmin-inp',    'value'),
        Input('vmax-inp',    'value'),
    )
    def _update(var_sel, step, mode,  # pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals,too-many-branches,too-many-statements
                slice_axis, slice_pos,
                iso_min_frac, iso_max_frac, iso_n, iso_caps,
                vol_opacity, vol_nsurf, vol_min_frac, vol_max_frac,
                cmap, log_chk, vmin_in, vmax_in):

        selected_var = var_sel or varname
        try:
            ad = _load(step, read_func)
        except (OSError, ValueError, EOFError) as exc:
            return no_update, [html.Span(f' Error loading step {step}: {exc}',
                                         style={'color': _RED})]
        if selected_var not in ad.variables:
            avail = ', '.join(sorted(ad.variables))
            return no_update, [html.Span(
                f' Variable {selected_var!r} not in step {step} '
                f'(available: {avail})', style={'color': _RED})]
        raw = ad.variables[selected_var]
        log = bool(log_chk and 'log' in log_chk)
        cmap = cmap or 'viridis'

        # Color range
        if vmin_in is not None:
            vmin = float(vmin_in)
        else:
            safe = raw[raw > 0] if log and np.any(raw > 0) else raw
            vmin = float(np.nanmin(safe))
        if vmax_in is not None:
            vmax = float(vmax_in)
        else:
            vmax = float(np.nanmax(raw))
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

        fig = go.Figure()
        title = ''

        if ad.ndim == 3:
            trace, title = _build_3d(
                ad, raw, selected_var, step, mode, cmap, _tf, cmin, cmax, cbar_title,
                slice_axis or 'z', float(slice_pos or 0.5),
                float(iso_min_frac or 0.2), float(iso_max_frac or 0.8),
                int(iso_n or 3), bool(iso_caps and 'caps' in iso_caps),
                float(vol_opacity or 0.1), int(vol_nsurf or 15),
                float(vol_min_frac or 0.0), float(vol_max_frac or 1.0),
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
            # Compute aspect ratio from domain extents so slices (which
            # have a constant coordinate on one axis) don't collapse that axis.
            dx = float(ad.x_cc[-1] - ad.x_cc[0]) if len(ad.x_cc) > 1 else 1.0
            dy = float(ad.y_cc[-1] - ad.y_cc[0]) if len(ad.y_cc) > 1 else 1.0
            dz = float(ad.z_cc[-1] - ad.z_cc[0]) if len(ad.z_cc) > 1 else 1.0
            max_d = max(dx, dy, dz, 1e-30)
            fig.update_layout(scene=dict(
                xaxis=dict(title='x', backgroundcolor=_SURF,
                           gridcolor=_OVER, color=_TEXT),
                yaxis=dict(title='y', backgroundcolor=_SURF,
                           gridcolor=_OVER, color=_TEXT),
                zaxis=dict(title='z', backgroundcolor=_SURF,
                           gridcolor=_OVER, color=_TEXT),
                bgcolor=_BG,
                aspectmode='manual',
                aspectratio=dict(x=dx/max_d, y=dy/max_d, z=dz/max_d),
            ))

        elif ad.ndim == 2:
            cbar = dict(title=dict(text=cbar_title, font=dict(color=_TEXT)),
                        tickfont=dict(color=_TEXT), tickformat='.2e')
            fig.add_trace(go.Heatmap(
                x=ad.x_cc, y=ad.y_cc, z=_tf(raw).T,
                zmin=cmin, zmax=cmax, colorscale=cmap, colorbar=cbar,
            ))
            fig.update_layout(
                xaxis=dict(title='x', color=_TEXT, gridcolor=_OVER, scaleanchor='y', exponentformat='e'),
                yaxis=dict(title='y', color=_TEXT, gridcolor=_OVER, exponentformat='e'),
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
            title = f'{selected_var}  ·  step {step}'

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

        fig.update_layout(
            title=dict(text=title, font=dict(color=_TEXT, size=13, family='monospace')),
            paper_bgcolor=_BG,
            font=dict(color=_TEXT, family='monospace'),
            margin=dict(l=0, r=80, t=36, b=0),
            uirevision=mode,   # preserve camera angle within a mode
        )

        dmin, dmax = float(np.nanmin(raw)), float(np.nanmax(raw))
        status = html.Div([
            html.Span(f'step {step}', style={'color': _YELLOW}),
            html.Span(f'  ·  shape {raw.shape}', style={'color': _MUTED}),
            html.Br(),
            html.Span('min ', style={'color': _MUTED}),
            html.Span(f'{dmin:.4g}', style={'color': _BLUE}),
            html.Span('  max ', style={'color': _MUTED}),
            html.Span(f'{dmax:.4g}', style={'color': _RED}),
        ])
        return fig, status

    # ------------------------------------------------------------------
    cons.print(f'\n[bold green]Interactive viz server:[/bold green] '
               f'[bold]http://{host}:{port}[/bold]')
    if host in ('127.0.0.1', 'localhost'):
        cons.print(
            f'\n[dim]To view from your laptop/desktop:[/dim]\n'
            f'[dim]  1. Open a [bold]new terminal[/bold] on your [bold]local[/bold] machine (not the cluster)[/dim]\n'
            f'[dim]  2. Run:[/dim]  [bold]ssh -L {port}:localhost:{port} <cluster-hostname>[/bold]\n'
            f'[dim]     (replace <cluster-hostname> with the host you SSH into, e.g. login-phoenix.pace.gatech.edu)[/dim]\n'
            f'[dim]  3. Open [bold]http://localhost:{port}[/bold] in your local browser[/dim]\n'
            f'[dim]  If your cluster requires 2FA, add to your ~/.ssh/config:[/dim]\n'
            f'[dim]     Host <cluster-hostname>[/dim]\n'
            f'[dim]       ControlMaster auto[/dim]\n'
            f'[dim]       ControlPath ~/.ssh/sockets/%r@%h-%p[/dim]\n'
            f'[dim]       ControlPersist 600[/dim]\n'
            f'[dim]  Then your first SSH session handles 2FA; the tunnel reuses it without re-authenticating.[/dim]\n'
            f'[dim]  (Run [bold]mkdir -p ~/.ssh/sockets[/bold] once if the directory does not exist.)[/dim]'
        )
    cons.print('[dim]\nCtrl+C to stop.[/dim]\n')
    app.run(debug=False, port=port, host=host)
