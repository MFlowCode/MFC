"""
Main entry point for the ``./mfc.sh viz`` command.

Dispatches to reader + renderer based on CLI arguments.
"""

import os
import sys

from mfc.state import ARG
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

    try:
        if ':' in str(step_arg):
            parts = str(step_arg).split(':')
            start = int(parts[0])
            end = int(parts[1])
            stride = int(parts[2]) if len(parts) > 2 else 1
            requested = list(range(start, end + 1, stride))
            return [s for s in requested if s in set(available_steps)]

        single = int(step_arg)
    except ValueError:
        cons.print(f"[bold red]Error:[/bold red] Invalid --step value '{step_arg}'. "
                   "Expected an integer, a range (start:end:stride), or 'all'.")
        sys.exit(1)

    if available_steps and single not in set(available_steps):
        return []
    return [single]


def viz():  # pylint: disable=too-many-locals,too-many-statements,too-many-branches
    """Main viz command dispatcher."""
    from .reader import discover_format, discover_timesteps, assemble  # pylint: disable=import-outside-toplevel
    from .renderer import render_1d, render_2d, render_3d_slice, render_mp4  # pylint: disable=import-outside-toplevel

    case_dir = ARG('input')
    if case_dir is None:
        cons.print("[bold red]Error:[/bold red] Please specify a case directory.")
        sys.exit(1)

    # Resolve case directory
    if not os.path.isdir(case_dir):
        cons.print(f"[bold red]Error:[/bold red] Directory not found: {case_dir}")
        sys.exit(1)

    # Auto-detect or use specified format
    fmt_arg = ARG('format')
    if fmt_arg:
        fmt = fmt_arg
    else:
        fmt = discover_format(case_dir)

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
            cons.print("[bold red]Error:[/bold red] No timesteps found.")
            sys.exit(1)

        if step_arg is None or step_arg == 'all':
            step = steps[0]
            cons.print(f"[dim]Using first available timestep: {step}[/dim]")
        else:
            try:
                step = int(step_arg)
            except ValueError:
                cons.print(f"[bold red]Error:[/bold red] Invalid --step value '{step_arg}'. "
                           "Expected an integer or 'all'.")
                sys.exit(1)

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

    # For rendering, --var and --step are required
    varname = ARG('var')
    step_arg = ARG('step')

    if varname is None:
        cons.print("[bold red]Error:[/bold red] --var is required for rendering. "
                   "Use --list-vars to see available variables.")
        sys.exit(1)

    steps = discover_timesteps(case_dir, fmt)
    if not steps:
        cons.print("[bold red]Error:[/bold red] No timesteps found.")
        sys.exit(1)

    requested_steps = _parse_steps(step_arg, steps)
    if not requested_steps:
        cons.print(f"[bold red]Error:[/bold red] No matching timesteps for --step {step_arg}")
        sys.exit(1)

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

    # Choose read function based on format
    def read_step(step):
        if fmt == 'silo':
            from .silo_reader import assemble_silo  # pylint: disable=import-outside-toplevel
            return assemble_silo(case_dir, step, var=varname)
        return assemble(case_dir, step, fmt, var=varname)

    # Validate variable name by reading the first timestep (without var filter)
    def read_step_all_vars(step):
        if fmt == 'silo':
            from .silo_reader import assemble_silo  # pylint: disable=import-outside-toplevel
            return assemble_silo(case_dir, step)
        return assemble(case_dir, step, fmt)

    test_assembled = read_step_all_vars(requested_steps[0])
    if varname not in test_assembled.variables:
        avail = sorted(test_assembled.variables.keys())
        cons.print(f"[bold red]Error:[/bold red] Variable '{varname}' not found.")
        cons.print(f"[bold]Available variables:[/bold] {', '.join(avail)}")
        sys.exit(1)

    # Create output directory
    output_base = ARG('output')
    if output_base is None:
        output_base = os.path.join(case_dir, 'viz')
    os.makedirs(output_base, exist_ok=True)

    # MP4 mode
    if ARG('mp4'):
        fps = ARG('fps') or 10
        mp4_path = os.path.join(output_base, f'{varname}.mp4')
        cons.print(f"[bold]Generating MP4:[/bold] {mp4_path} ({len(requested_steps)} frames)")
        success = render_mp4(varname, requested_steps, mp4_path,
                             fps=int(fps), read_func=read_step, **render_opts)
        if success:
            cons.print(f"[bold green]Done:[/bold green] {mp4_path}")
        else:
            cons.print(f"[bold red]Error:[/bold red] Failed to generate {mp4_path}. "
                       "Ensure imageio and imageio-ffmpeg are installed.")
        return

    # Single or multiple PNG frames
    try:
        from tqdm import tqdm  # pylint: disable=import-outside-toplevel
        step_iter = tqdm(requested_steps, desc='Rendering')
    except ImportError:
        step_iter = requested_steps

    for step in step_iter:
        assembled = read_step(step)
        output_path = os.path.join(output_base, f'{varname}_{step}.png')

        if assembled.ndim == 1:
            render_1d(assembled.x_cc, assembled.variables[varname],
                      varname, step, output_path, **render_opts)
        elif assembled.ndim == 2:
            render_2d(assembled.x_cc, assembled.y_cc,
                      assembled.variables[varname],
                      varname, step, output_path, **render_opts)
        elif assembled.ndim == 3:
            render_3d_slice(assembled, varname, step, output_path, **render_opts)

        if len(requested_steps) == 1:
            cons.print(f"[bold green]Saved:[/bold green] {output_path}")

    if len(requested_steps) > 1:
        cons.print(f"[bold green]Saved {len(requested_steps)} frames to:[/bold green] {output_base}/")
