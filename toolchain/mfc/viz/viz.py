"""
Main entry point for the ``./mfc.sh viz`` command.

Dispatches to reader + renderer based on CLI arguments.
"""

import os

from mfc.state import ARG
from mfc.common import MFCException
from mfc.printer import cons


def _parse_steps(step_arg, available_steps):
    """
    Parse the --step argument into a list of timestep integers.

    Formats:
      - Single int: "1000"
      - Range: "0:10000:500"  (start:end:stride)
      - "all": all available timesteps
    """
    if step_arg is None or step_arg == 'all':
        return available_steps

    if step_arg == 'last':
        return [available_steps[-1]] if available_steps else []

    try:
        if ':' in str(step_arg):
            parts = str(step_arg).split(':')
            start = int(parts[0])
            end = int(parts[1])
            stride = int(parts[2]) if len(parts) > 2 else 1
            requested = list(range(start, end + 1, stride))
            return [s for s in requested if s in set(available_steps)]

        single = int(step_arg)
    except ValueError as exc:
        raise MFCException(f"Invalid --step value '{step_arg}'. "
                           "Expected an integer, a range (start:end:stride), 'last', or 'all'.") from exc

    if available_steps and single not in set(available_steps):
        return []
    return [single]


def viz():  # pylint: disable=too-many-locals,too-many-statements,too-many-branches
    """Main viz command dispatcher."""
    from .reader import discover_format, discover_timesteps, assemble  # pylint: disable=import-outside-toplevel
    from .renderer import render_1d, render_1d_tiled, render_2d, render_3d_slice, render_mp4  # pylint: disable=import-outside-toplevel

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
            raise MFCException(str(exc)) from exc

    cons.print(f"[bold]Format:[/bold] {fmt}")

    # Quick guide when no action is specified
    no_action = (not ARG('list_steps') and not ARG('list_vars')
                 and ARG('var') is None and ARG('step') is None
                 and not ARG('interactive') and not ARG('tui'))
    if no_action:
        cons.print()
        d = case_dir
        cons.print("[bold]Quick start:[/bold]")
        cons.print(f"  [green]./mfc.sh viz {d} --list-steps[/green]"
                   "                   [dim]see available timesteps[/dim]")
        cons.print(f"  [green]./mfc.sh viz {d} --list-vars --step 0[/green]"
                   "            [dim]see available variables[/dim]")
        cons.print(f"  [green]./mfc.sh viz {d} --var pres --step last[/green]"
                   "            [dim]render a PNG[/dim]")
        cons.print(f"  [green]./mfc.sh viz {d} --var pres --step all --mp4[/green]"
                   "     [dim]render an MP4[/dim]")
        cons.print(f"  [green]./mfc.sh viz {d} --var pres --step 0 --slice-axis z[/green]"
                   " [dim]3D midplane slice[/dim]")
        cons.print()
        cons.print("[dim]Run [bold]./mfc.sh viz --help[/bold] for all options.[/dim]")
        return

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
            raise MFCException("No timesteps found.")

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
                    f"Timestep {step} not found. Available range: "
                    f"{steps[0]} to {steps[-1]} ({len(steps)} timesteps)")

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

    # For rendering, --step is required; --var is optional for 1D (shows all)
    varname = ARG('var')
    step_arg = ARG('step')
    tiled = varname is None or varname == 'all'

    if step_arg is None:
        if ARG('interactive') or ARG('tui'):
            step_arg = 'all'   # default to all steps in interactive/TUI mode
        else:
            raise MFCException("--step is required for rendering. "
                               "Use --list-steps to see available timesteps, or pass --step all.")

    steps = discover_timesteps(case_dir, fmt)
    if not steps:
        raise MFCException("No timesteps found.")

    requested_steps = _parse_steps(step_arg, steps)
    if not requested_steps:
        msg = f"No matching timesteps for --step {step_arg}"
        if steps:
            msg += f". Available range: {steps[0]} to {steps[-1]} ({len(steps)} timesteps)"
        raise MFCException(msg)

    # Collect rendering options
    render_opts = {}
    cmap = ARG('cmap')
    if cmap:
        render_opts['cmap'] = cmap
    vmin = ARG('vmin')
    if vmin is not None:
        render_opts['vmin'] = float(vmin)
    vmax = ARG('vmax')
    if vmax is not None:
        render_opts['vmax'] = float(vmax)
    dpi = ARG('dpi')
    if dpi is not None:
        render_opts['dpi'] = int(dpi)
    if ARG('log_scale'):
        render_opts['log_scale'] = True

    slice_axis = ARG('slice_axis')
    slice_index = ARG('slice_index')
    slice_value = ARG('slice_value')
    if slice_axis:
        render_opts['slice_axis'] = slice_axis
    if slice_index is not None:
        render_opts['slice_index'] = int(slice_index)
    if slice_value is not None:
        render_opts['slice_value'] = float(slice_value)

    interactive = ARG('interactive')

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

    # Guard against loading too many 3D timesteps (memory)
    if test_assembled.ndim == 3 and len(requested_steps) > 500:
        raise MFCException(
            f"Refusing to load {len(requested_steps)} timesteps for 3D data "
            "(limit is 500). Use --step with a range or stride to reduce.")

    # Tiled mode for non-TUI, non-interactive rendering only works for 1D
    if tiled and not interactive and not ARG('tui'):
        if test_assembled.ndim != 1:
            raise MFCException("--var is required for 2D/3D rendering. "
                               "Use --list-vars to see available variables.")

    if not tiled and not interactive and not ARG('tui') and varname not in test_assembled.variables:
        raise MFCException(f"Variable '{varname}' not found. "
                           f"Available variables: {', '.join(avail)}")

    # TUI mode — launch Textual terminal UI (1D/2D only)
    if ARG('tui'):
        if test_assembled.ndim == 3:
            raise MFCException(
                "--tui only supports 1D and 2D data. "
                "Use --interactive for 3D data."
            )
        from .tui import run_tui  # pylint: disable=import-outside-toplevel
        init_var = varname if varname in avail else (avail[0] if avail else None)
        run_tui(init_var, requested_steps, read_step, ndim=test_assembled.ndim)
        return

    # Interactive mode — launch Dash web server
    if interactive:
        from .interactive import run_interactive  # pylint: disable=import-outside-toplevel
        port = ARG('port') or 8050
        # Default to first available variable if --var was not specified
        init_var = varname if varname in avail else (avail[0] if avail else None)
        run_interactive(init_var, requested_steps, read_step, port=int(port))
        return

    # Create output directory
    output_base = ARG('output')
    if output_base is None:
        output_base = os.path.join(case_dir, 'viz')
    os.makedirs(output_base, exist_ok=True)

    # MP4 mode
    if ARG('mp4'):
        fps = ARG('fps') or 10
        label = 'all' if tiled else varname
        mp4_path = os.path.join(output_base, f'{label}.mp4')
        cons.print(f"[bold]Generating MP4:[/bold] {mp4_path} ({len(requested_steps)} frames)")
        success = render_mp4(varname, requested_steps, mp4_path,
                             fps=int(fps), read_func=read_step,
                             tiled=tiled, **render_opts)
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
            assembled = read_step(step)
        except (FileNotFoundError, EOFError, ValueError) as exc:
            cons.print(f"[yellow]Warning:[/yellow] Skipping step {step}: {exc}")
            failures.append(step)
            continue

        output_path = os.path.join(output_base, f'{label}_{step}.png')

        if tiled and assembled.ndim == 1:
            render_1d_tiled(assembled.x_cc, assembled.variables,
                            step, output_path, **render_opts)
        elif assembled.ndim == 1:
            render_1d(assembled.x_cc, assembled.variables[varname],
                      varname, step, output_path, **render_opts)
        elif assembled.ndim == 2:
            render_2d(assembled.x_cc, assembled.y_cc,
                      assembled.variables[varname],
                      varname, step, output_path, **render_opts)
        elif assembled.ndim == 3:
            render_3d_slice(assembled, varname, step, output_path, **render_opts)
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
