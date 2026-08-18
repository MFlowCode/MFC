"""
Help output and onboarding for the MFC toolchain.

This module provides:
- Enhanced help output with Rich formatting
- Onboarding for new users
"""

import os

from rich import box
from rich.panel import Panel

# Import command definitions from CLI schema (SINGLE SOURCE OF TRUTH)
from .cli.commands import COMMANDS
from .common import MFC_ROOT_DIR
from .printer import cons

# ENHANCED HELP OUTPUT


def _truncate_desc(desc: str, max_len: int = 50) -> str:
    """Truncate description to fit compact display."""
    if len(desc) <= max_len:
        return desc
    return desc[: max_len - 3] + "..."


def print_help():
    """Print compact, colorized help overview."""

    # Header (no box)
    cons.print()
    cons.print("[bold cyan]MFC[/bold cyan] - Multi-component Flow Code")
    cons.print("[dim]Exascale CFD solver for compressible multi-phase flows[/dim]")
    cons.print()

    # Commands section - compact format (using COMMANDS as source of truth)
    cons.print("[bold]Commands:[/bold]")

    # Primary commands (shown prominently with aliases)
    primary = ["build", "run", "test", "viz", "validate", "new", "clean"]
    for cmd in primary:
        if cmd not in COMMANDS:
            continue
        info = COMMANDS[cmd]
        alias = info.get("alias") or ""
        alias_str = f" ({alias})" if alias else "    "
        desc = _truncate_desc(info["description"])
        cons.print(f"  [green]{cmd:9}[/green][dim]{alias_str:4}[/dim] {desc}")

    # Secondary commands (dimmed)
    secondary = ["params", "load"]
    for cmd in secondary:
        if cmd not in COMMANDS:
            continue
        desc = _truncate_desc(COMMANDS[cmd]["description"])
        cons.print(f"  [dim]{cmd:13} {desc}[/dim]")

    cons.print()

    # Quick start - single line
    cons.print("[bold]Quick start:[/bold] [cyan]./mfc.sh new my_case[/cyan] → edit case.py → [cyan]./mfc.sh build[/cyan] → [cyan]./mfc.sh run[/cyan]")

    # Footer
    cons.print("[dim]Run ./mfc.sh <command> --help for options[/dim]")
    cons.print()


def print_command_help(command: str, show_argparse: bool = True):
    """Print enhanced help for a specific command."""
    if command not in COMMANDS:
        cons.print(f"[red]Unknown command: {command}[/red]")
        return False

    cmd = COMMANDS[command]
    alias = cmd.get("alias", "")
    alias_str = f" [dim](alias: {alias})[/dim]" if alias else ""

    # Header panel
    cons.print()
    cons.raw.print(Panel(f"[bold cyan]{command}[/bold cyan]{alias_str}\n[dim]{cmd['description']}[/dim]", box=box.ROUNDED, padding=(0, 2)))
    cons.print()

    # Examples
    if cmd.get("examples"):
        cons.print("[bold]Examples:[/bold]")
        for example, desc in cmd["examples"]:
            cons.print(f"  [green]{example}[/green]")
            cons.print(f"      [dim]{desc}[/dim]")
        cons.print()

    # Key options
    if cmd.get("key_options"):
        cons.print("[bold]Key Options:[/bold]")
        for opt, desc in cmd["key_options"]:
            if opt.startswith("-- ") and opt.endswith(" --"):
                cons.print(f"  [bold yellow]{opt}[/bold yellow]")
            else:
                cons.print(f"  [cyan]{opt:24}[/cyan] {desc}")
        cons.print()
        if show_argparse:
            cons.print("[dim]Run with --help for full option list[/dim]")
            cons.print()

    return True


# ONBOARDING FOR NEW USERS


def is_first_time_user() -> bool:
    """Check if this is a first-time user (no build directory)."""
    build_dir = os.path.join(MFC_ROOT_DIR, "build")
    return not os.path.exists(build_dir)


def print_welcome():
    """Print welcome message for new users."""
    cons.print()
    cons.raw.print(
        Panel(
            "[bold cyan]Welcome to MFC![/bold cyan]\n\n"
            "It looks like this is your first time using MFC. Here's how to get started:\n\n"
            "  [green]1.[/green] [bold]Load environment[/bold] (HPC clusters):\n"
            "     [cyan]source ./mfc.sh load -c <cluster> -m <mode>[/cyan]\n"
            "     Example: [dim]source ./mfc.sh load -c p -m g[/dim] (Phoenix, GPU)\n\n"
            "  [green]2.[/green] [bold]Create a new case[/bold]:\n"
            "     [cyan]./mfc.sh new my_first_case[/cyan]\n\n"
            "  [green]3.[/green] [bold]Build MFC[/bold]:\n"
            "     [cyan]./mfc.sh build -j $(nproc)[/cyan]\n\n"
            "  [green]4.[/green] [bold]Run your simulation[/bold]:\n"
            "     [cyan]./mfc.sh run my_first_case/case.py[/cyan]\n\n"
            "[bold yellow]Optional:[/bold yellow] Enable tab completion for your shell:\n"
            "     [cyan]./mfc.sh completion install[/cyan]\n\n"
            "[dim]Run [cyan]./mfc.sh --help[/cyan] for all available commands[/dim]",
            title="[bold]Getting Started[/bold]",
            box=box.DOUBLE,
            border_style="cyan",
            padding=(1, 2),
        )
    )
    cons.print()
