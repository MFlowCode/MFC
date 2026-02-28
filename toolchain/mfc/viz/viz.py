"""
Main entry point for the ``./mfc.sh viz`` command.

Dispatches to reader + renderer based on CLI arguments.
"""

import os
import importlib
import importlib.util
import shutil
import subprocess
import sys

from mfc.state import ARG
from mfc.common import MFC_ROOT_DIR, MFCException
from mfc.printer import cons


def _ensure_viz_deps() -> None:
    """Install the [viz] optional extras on first use.

    Checks one sentinel per feature group so that a user who has matplotlib
    pre-installed (e.g. via Anaconda) but lacks imageio, textual, or h5py
    still triggers the install.
    """
    _SENTINELS = ("matplotlib", "imageio", "h5py", "textual", "dash", "plotext", "plotly")
    if all(importlib.util.find_spec(p) is not None for p in _SENTINELS):
        return  # all present

    toolchain_path = os.path.join(MFC_ROOT_DIR, "toolchain")
    cons.print("[bold]Installing viz dependencies[/bold] "
               "(first run — this may take a minute)...")

    env = {**os.environ, "UV_LINK_MODE": "copy"}
    if shutil.which("uv"):
        cmd = ["uv", "pip", "install", f"{toolchain_path}[viz]"]
    else:
        cmd = [sys.executable, "-m", "pip", "install", f"{toolchain_path}[viz]"]

    try:
        result = subprocess.run(cmd, env=env, check=False, timeout=300)
    except subprocess.TimeoutExpired as exc:
        raise MFCException(
            "Timed out installing viz dependencies (network may be restricted). "
            f"Try manually: pip install '{toolchain_path}[viz]'"
        ) from exc
    if result.returncode != 0:
        raise MFCException(
            "Failed to install viz dependencies. "
            f"Try manually: pip install '{toolchain_path}[viz]'"
        )

    importlib.invalidate_caches()
    cons.print("[bold green]Viz dependencies installed.[/bold green]\n")


_CMAP_POPULAR = (
    'viridis, plasma, inferno, magma, turbo, '
    'coolwarm, RdBu_r, bwr, hot, jet, gray, seismic'
)


def _validate_cmap(name):
    """Raise a helpful MFCException if *name* is not a known matplotlib colormap."""
    import matplotlib  # pylint: disable=import-outside-toplevel
    if name in matplotlib.colormaps:
        return
    try:
        from rapidfuzz import process  # pylint: disable=import-outside-toplevel
        matches = process.extract(name, list(matplotlib.colormaps), limit=3)
        suggestions = ', '.join(m[0] for m in matches)
        hint = f" Did you mean: {suggestions}?"
    except ImportError:
        hint = f" Popular choices: {_CMAP_POPULAR}."
    raise MFCException(f"Unknown colormap '{name}'.{hint}")


def _steps_hint(steps, n=8):
    """Short inline preview of available steps for error messages."""
    if not steps:
        return "none found"
    if len(steps) <= n:
        return ', '.join(str(s) for s in steps)
    half = n // 2
    head = ', '.join(str(s) for s in steps[:half])
    tail = ', '.join(str(s) for s in steps[-half:])
    return f"{head}, ... [{len(steps)} total] ..., {tail}"


def _parse_step_list(s, available_steps):
    """
    Parse a comma-separated step list, with optional '...' ellipsis expansion.

    Examples:
      "0,100,200,1000"        -> [0, 100, 200, 1000] (intersection with available)
      "0,100,200,...,1000"    -> range(0, 1001, 100) (infers stride=100 from last pair)
    """
    parts = [p.strip() for p in s.split(',')]
    avail_set = set(available_steps)

    if '...' in parts:
        idx = parts.index('...')
        if idx < 2:
            raise MFCException(
                f"Invalid --step value '{s}'. "
                "Ellipsis '...' requires at least two values before it, "
                "e.g. '0,100,...,1000'."
            )
        if idx != len(parts) - 2:
            raise MFCException(
                f"Invalid --step value '{s}'. "
                "Ellipsis '...' must be the second-to-last item, "
                "e.g. '0,100,...,1000'."
            )
        try:
            prefix = [int(p) for p in parts[:idx]]
            end = int(parts[idx + 1])
        except ValueError as exc:
            raise MFCException(
                f"Invalid --step value '{s}': all values must be integers."
            ) from exc

        stride = prefix[-1] - prefix[-2]
        if stride <= 0:
            raise MFCException(
                f"Invalid --step value '{s}': "
                f"ellipsis stride must be positive (got {stride})."
            )
        requested = list(range(prefix[0], end + 1, stride))
    else:
        try:
            requested = [int(p) for p in parts]
        except ValueError as exc:
            raise MFCException(
                f"Invalid --step value '{s}': all values must be integers."
            ) from exc

    return [t for t in requested if t in avail_set]


def _parse_steps(step_arg, available_steps):
    """
    Parse the --step argument into a list of timestep integers.

    Formats:
      - Single int:       "1000"
      - Range:            "0:10000:500"  (start:end:stride)
      - Comma list:       "0,100,200,1000"
      - Ellipsis list:    "0,100,200,...,1000"  (stride inferred from last pair)
      - "last":           last available timestep
      - "all":            all available timesteps
    """
    if step_arg is None or step_arg == 'all':
        return available_steps

    if step_arg == 'last':
        return [available_steps[-1]] if available_steps else []

    s = str(step_arg)

    if ',' in s:
        return _parse_step_list(s, available_steps)

    try:
        if ':' in s:
            parts = s.split(':')
            start = int(parts[0])
            end = int(parts[1])
            stride = int(parts[2]) if len(parts) > 2 else 1
            requested = list(range(start, end + 1, stride))
            return [t for t in requested if t in set(available_steps)]

        single = int(s)
    except ValueError as exc:
        raise MFCException(
            f"Invalid --step value '{step_arg}'. "
            "Expected an integer, a range (start:end:stride), "
            "a comma list (0,100,200), an ellipsis list (0,100,...,1000), "
            "'last', or 'all'."
        ) from exc

    if available_steps and single not in set(available_steps):
        return []
    return [single]


def viz():  # pylint: disable=too-many-locals,too-many-statements,too-many-branches
    """Main viz command dispatcher."""
    _ensure_viz_deps()

    from .reader import discover_format, discover_timesteps, assemble, has_lag_bubble_evol, read_lag_bubbles_at_step  # pylint: disable=import-outside-toplevel
    from .renderer import render_1d, render_1d_tiled, render_2d, render_2d_tiled, render_3d_slice, render_mp4  # pylint: disable=import-outside-toplevel

    case_dir = ARG('input')
    if case_dir is None:
        raise MFCException("Please specify a case directory.")

    # Resolve case directory
    if not os.path.isdir(case_dir):
        raise MFCException(f"Directory not found: {case_dir}")

    # Auto-detect or use specified format
    fmt_arg = ARG('format')
    if fmt_arg:
        if fmt_arg not in ('binary', 'silo'):
            raise MFCException(f"Unknown format '{fmt_arg}'. "
                               "Supported formats: 'binary', 'silo'.")
        fmt = fmt_arg
    else:
        try:
            fmt = discover_format(case_dir)
        except FileNotFoundError as exc:
            msg = str(exc)
            if os.path.isfile(os.path.join(case_dir, 'case.py')):
                msg += (" This looks like an MFC case directory. "
                        "Did you forget to run post_process first?")
            raise MFCException(msg) from exc

    cons.print(f"[bold]Format:[/bold] {fmt}")

    # Handle --list-steps
    if ARG('list_steps'):
        steps = discover_timesteps(case_dir, fmt)
        if not steps:
            cons.print("[yellow]No timesteps found.[/yellow]")
        else:
            cons.print(f"[bold]Available timesteps ({len(steps)}):[/bold]")
            # Print in columns
            line = ""
            for i, s in enumerate(steps):
                line += f"{s:>8}"
                if (i + 1) % 10 == 0:
                    cons.print(line)
                    line = ""
            if line:
                cons.print(line)
        return

    # Handle --list-vars (requires --step)
    if ARG('list_vars'):
        step_arg = ARG('step')
        steps = discover_timesteps(case_dir, fmt)
        if not steps:
            raise MFCException(
                f"No timesteps found in '{case_dir}' ({fmt} format). "
                "Ensure post_process has been run and produced output files.")

        if step_arg is None or step_arg == 'all':
            step = steps[0]
            cons.print(f"[dim]Using first available timestep: {step}[/dim]")
        elif step_arg == 'last':
            step = steps[-1]
            cons.print(f"[dim]Using last available timestep: {step}[/dim]")
        else:
            try:
                step = int(step_arg)
            except ValueError as exc:
                raise MFCException(f"Invalid --step value '{step_arg}'. "
                                   "Expected an integer, 'last', or 'all'.") from exc
            if step not in steps:
                raise MFCException(
                    f"Timestep {step} not found. "
                    f"Available steps: {_steps_hint(steps)}")

        if fmt == 'silo':
            from .silo_reader import assemble_silo  # pylint: disable=import-outside-toplevel
            assembled = assemble_silo(case_dir, step)
        else:
            assembled = assemble(case_dir, step, fmt)

        varnames = sorted(assembled.variables.keys())
        cons.print(f"[bold]Available variables ({len(varnames)}):[/bold]")
        for vn in varnames:
            data = assembled.variables[vn]
            cons.print(f"  {vn:<20s}  min={data.min():.6g}  max={data.max():.6g}")
        return

    # For rendering, --step is required; --var is optional for 1D/2D (shows all in tiled layout)
    varname = ARG('var')
    step_arg = ARG('step')
    tiled = varname is None or varname == 'all'

    if ARG('interactive') or ARG('tui') or ARG('mp4'):
        # Load all steps by default; honour an explicit --step so users can
        # reduce the set for large 3D cases before hitting the step limit.
        if step_arg == 'last':
            step_arg = 'all'

    steps = discover_timesteps(case_dir, fmt)
    if not steps:
        raise MFCException(
            f"No timesteps found in '{case_dir}' ({fmt} format). "
            "Ensure post_process has been run and produced output files.")

    requested_steps = _parse_steps(step_arg, steps)
    if not requested_steps:
        raise MFCException(
            f"No matching timesteps for --step {step_arg!r}. "
            f"Available steps: {_steps_hint(steps)}")

    # Collect rendering options
    render_opts = {
        'cmap':  ARG('cmap'),
        'dpi':   ARG('dpi'),
        'slice_axis': ARG('slice_axis'),
    }
    if ARG('vmin') is not None:
        render_opts['vmin'] = float(ARG('vmin'))
    if ARG('vmax') is not None:
        render_opts['vmax'] = float(ARG('vmax'))
    if ARG('log_scale'):
        render_opts['log_scale'] = True
    if ARG('slice_index') is not None and ARG('slice_value') is not None:
        raise MFCException("--slice-index and --slice-value are mutually exclusive.")
    if ARG('slice_index') is not None:
        render_opts['slice_index'] = int(ARG('slice_index'))
    if ARG('slice_value') is not None:
        render_opts['slice_value'] = float(ARG('slice_value'))

    interactive = ARG('interactive')

    # Lagrange bubble overlay: auto-detect D/lag_bubble_evol_*.dat files
    bubble_func = None
    if has_lag_bubble_evol(case_dir):
        bubble_func = lambda s: read_lag_bubbles_at_step(case_dir, s)  # noqa: E731
        cons.print("[dim]Lagrange bubble positions detected — overlaying on plots.[/dim]")

    # Load all variables when tiled, interactive, or TUI; filter otherwise.
    # TUI needs all vars loaded so the sidebar can switch between them.
    load_all = tiled or interactive or ARG('tui')

    def read_step(step):
        if fmt == 'silo':
            from .silo_reader import assemble_silo  # pylint: disable=import-outside-toplevel
            return assemble_silo(case_dir, step, var=None if load_all else varname)
        return assemble(case_dir, step, fmt, var=None if load_all else varname)

    # Validate variable name / discover available variables
    test_assembled = read_step(requested_steps[0])
    avail = sorted(test_assembled.variables.keys())

    # Guard against loading too many 3D timesteps (memory).
    # Interactive mode caches all steps simultaneously, so use a tighter limit.
    _3d_limit = 50 if interactive else 500
    if test_assembled.ndim == 3 and len(requested_steps) > _3d_limit:
        raise MFCException(
            f"Refusing to load {len(requested_steps)} timesteps for 3D data "
            f"(limit is {_3d_limit}). Use --step with a range or stride to reduce.")

    # Tiled mode works for 1D and 2D.  For 3D, auto-select the first variable.
    if tiled and not interactive and not ARG('tui'):
        if test_assembled.ndim == 3:
            varname = avail[0] if avail else None
            if varname is None:
                raise MFCException("No variables found in output.")
            tiled = False
            cons.print(f"[dim]Auto-selected variable: [bold]{varname}[/bold]"
                       " (use --var to specify)[/dim]")

    if not tiled and not interactive and not ARG('tui') and varname not in test_assembled.variables:
        # test_assembled was loaded with var_filter=varname so its variables dict
        # may be empty. Re-read without filter (errors only, so extra I/O is fine)
        # to build a useful "available variables" list for the error message.
        if not avail:
            if fmt == 'silo':
                from .silo_reader import assemble_silo  # pylint: disable=import-outside-toplevel
                _full = assemble_silo(case_dir, requested_steps[0])
            else:
                _full = assemble(case_dir, requested_steps[0], fmt)
            avail = sorted(_full.variables.keys())
        avail_str = ', '.join(avail) if avail else '(none — check post_process output)'
        raise MFCException(
            f"Variable '{varname}' not found. "
            f"Available: {avail_str}. "
            f"Use --list-vars to see variables at a given step."
        )

    # Validate colormap early so all modes get a clean error for bad --cmap
    cmap_name = ARG('cmap')
    if cmap_name:
        _validate_cmap(cmap_name)

    # TUI mode — launch Textual terminal UI (1D/2D only)
    if ARG('tui'):
        if test_assembled.ndim == 3:
            raise MFCException(
                "--tui only supports 1D and 2D data. "
                "Use --interactive for 3D data."
            )
        from .tui import run_tui  # pylint: disable=import-outside-toplevel
        init_var = varname if varname in avail else (avail[0] if avail else None)
        run_tui(init_var, requested_steps, read_step, ndim=test_assembled.ndim,
                bubble_func=bubble_func)
        return

    # Interactive mode — launch Dash web server
    if interactive:
        from .interactive import run_interactive  # pylint: disable=import-outside-toplevel
        port = ARG('port')
        host = ARG('host')
        # Default to first available variable if --var was not specified
        init_var = varname if varname in avail else (avail[0] if avail else None)
        run_interactive(init_var, requested_steps, read_step,
                        port=int(port), host=str(host),
                        bubble_func=bubble_func)
        return

    # Create output directory
    output_base = ARG('output')
    if output_base is None:
        output_base = os.path.join(case_dir, 'viz')
    os.makedirs(output_base, exist_ok=True)

    # MP4 mode
    if ARG('mp4'):
        fps = float(ARG('fps'))
        _FPS_DEFAULT = 10.0
        _MIN_DURATION = 5.0  # seconds
        n_frames = len(requested_steps)
        if fps == _FPS_DEFAULT and n_frames / fps < _MIN_DURATION:
            fps = max(1.0, n_frames / _MIN_DURATION)
            cons.print(
                f"[dim]Auto-adjusted fps to {fps:.2g} "
                f"so video is at least {_MIN_DURATION:.0f}s "
                f"(use --fps to override).[/dim]"
            )
        label = 'all' if tiled else varname
        mp4_path = os.path.join(output_base, f'{label}.mp4')
        cons.print(f"[bold]Generating MP4:[/bold] {mp4_path} ({n_frames} frames @ {fps:.2g} fps = {n_frames/fps:.1f}s)")
        success = render_mp4(varname, requested_steps, mp4_path,
                             fps=fps, read_func=read_step,
                             tiled=tiled, bubble_func=bubble_func, **render_opts)
        if success:
            cons.print(f"[bold green]Done:[/bold green] {mp4_path}")
        else:
            raise MFCException(f"Failed to generate {mp4_path}. "
                               "Ensure imageio and imageio-ffmpeg are installed.")
        return

    # Single or multiple PNG frames
    try:
        from tqdm import tqdm  # pylint: disable=import-outside-toplevel
        step_iter = tqdm(requested_steps, desc='Rendering')
    except ImportError:
        step_iter = requested_steps

    failures = []
    label = 'all' if tiled else varname
    for step in step_iter:
        try:
            # Reuse the already-loaded probe data for the first step
            assembled = test_assembled if step == requested_steps[0] else read_step(step)
        except (FileNotFoundError, EOFError, ValueError) as exc:
            cons.print(f"[yellow]Warning:[/yellow] Skipping step {step}: {exc}")
            failures.append(step)
            continue

        output_path = os.path.join(output_base, f'{label}_{step}.png')

        # Inject bubble positions for this step
        step_opts = render_opts
        if bubble_func is not None:
            try:
                step_opts = dict(render_opts, bubbles=bubble_func(step))
            except Exception as exc:  # pylint: disable=broad-except
                cons.print(f"[yellow]Warning:[/yellow] Skipping bubble overlay for step {step}: {exc}")

        if tiled and assembled.ndim == 1:
            render_1d_tiled(assembled.x_cc, assembled.variables,
                            step, output_path, **step_opts)
        elif tiled and assembled.ndim == 2:
            render_2d_tiled(assembled, step, output_path, **step_opts)
        elif assembled.ndim == 1:
            render_1d(assembled.x_cc, assembled.variables[varname],
                      varname, step, output_path, **step_opts)
        elif assembled.ndim == 2:
            render_2d(assembled.x_cc, assembled.y_cc,
                      assembled.variables[varname],
                      varname, step, output_path, **step_opts)
        elif assembled.ndim == 3:
            render_3d_slice(assembled, varname, step, output_path, **step_opts)
        else:
            cons.print(f"[yellow]Warning:[/yellow] Unsupported ndim={assembled.ndim} "
                       f"for step {step}, skipping.")
            failures.append(step)
            continue

        if len(requested_steps) == 1:
            cons.print(f"[bold green]Saved:[/bold green] {output_path}")

    rendered = len(requested_steps) - len(failures)
    if failures:
        cons.print(f"[yellow]Warning:[/yellow] {len(failures)} step(s) failed: "
                   f"{failures[:10]}{'...' if len(failures) > 10 else ''}")
    if rendered > 1:
        cons.print(f"[bold green]Saved {rendered} frames to:[/bold green] {output_base}/")
