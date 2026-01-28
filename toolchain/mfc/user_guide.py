"""
User guide, help, tips, and onboarding for MFC toolchain.

This module provides:
- Enhanced help output with Rich formatting
- Contextual tips after errors/failures
- Interactive mode with menu
- Onboarding for new users
"""

import os

from rich.panel import Panel
from rich.table import Table
from rich.prompt import Prompt
from rich import box

from .printer import cons
from .common import MFC_ROOT_DIR


# =============================================================================
# COMMAND DEFINITIONS
# =============================================================================

COMMANDS = {
    "build": {
        "description": "Build MFC and its dependencies",
        "examples": [
            ("./mfc.sh build", "Build all targets"),
            ("./mfc.sh build -t simulation", "Build only simulation"),
            ("./mfc.sh build -j 8", "Build with 8 parallel jobs"),
            ("./mfc.sh build --gpu", "Build with GPU (OpenACC) support"),
        ]
    },
    "run": {
        "description": "Run a simulation case",
        "examples": [
            ("./mfc.sh run case.py", "Run a case file"),
            ("./mfc.sh run case.py -n 4", "Run with 4 MPI ranks"),
            ("./mfc.sh run case.py -e batch", "Submit as batch job"),
        ]
    },
    "test": {
        "description": "Run the MFC test suite",
        "examples": [
            ("./mfc.sh test", "Run all tests"),
            ("./mfc.sh test -j 4", "Run tests with 4 parallel jobs"),
            ("./mfc.sh test --only <UUID>", "Run a specific test"),
        ]
    },
    "clean": {
        "description": "Remove build artifacts and cache",
        "examples": [
            ("./mfc.sh clean", "Clean all build files"),
        ]
    },
    "validate": {
        "description": "Check a case file for errors",
        "examples": [
            ("./mfc.sh validate case.py", "Validate case file syntax and constraints"),
        ]
    },
    "init": {
        "description": "Create a new case from a template",
        "examples": [
            ("./mfc.sh init my_case", "Create case with default 1D template"),
            ("./mfc.sh init my_case -t 2D_minimal", "Create case with 2D template"),
            ("./mfc.sh init --list", "List available templates"),
        ]
    },
    "count": {
        "description": "Count lines of code in MFC",
        "examples": [
            ("./mfc.sh count", "Show LOC statistics"),
        ]
    },
    "packer": {
        "description": "Pack/unpack/compare simulation data",
        "examples": [
            ("./mfc.sh packer pack case.py", "Pack case output"),
            ("./mfc.sh packer compare a.pack b.pack", "Compare two packed files"),
        ]
    },
    "load": {
        "description": "Load MFC environment (use with source)",
        "examples": [
            ("source ./mfc.sh load -c p -m g", "Load Phoenix GPU modules"),
            ("source ./mfc.sh load -c f -m c", "Load Frontier CPU modules"),
        ]
    },
}


# =============================================================================
# ENHANCED HELP OUTPUT
# =============================================================================

def print_help():
    """Print enhanced, colorized help overview."""

    # Header
    cons.print()
    cons.raw.print(Panel(
        "[bold cyan]MFC[/bold cyan] - [dim]Multi-component Flow Code[/dim]\n"
        "[dim]Exascale CFD solver for compressible multi-phase flows[/dim]",
        box=box.ROUNDED,
        padding=(0, 2)
    ))
    cons.print()

    # Commands table
    table = Table(
        title="[bold]Commands[/bold]",
        box=box.SIMPLE,
        show_header=True,
        header_style="bold cyan",
        title_justify="left",
        padding=(0, 2)
    )
    table.add_column("Command", style="green", no_wrap=True)
    table.add_column("Description", style="white")

    # Primary commands
    for cmd in ["build", "run", "test", "validate", "init", "clean"]:
        table.add_row(cmd, COMMANDS[cmd]["description"])

    table.add_row("", "")  # Spacer

    # Secondary commands
    for cmd in ["count", "packer", "load"]:
        table.add_row(f"[dim]{cmd}[/dim]", f"[dim]{COMMANDS[cmd]['description']}[/dim]")

    cons.raw.print(table)
    cons.print()

    # Quick start
    cons.raw.print(Panel(
        "[bold]Quick Start[/bold]\n\n"
        "  [green]1.[/green] [cyan]./mfc.sh init my_case[/cyan]           Create a new case\n"
        "  [green]2.[/green] [cyan]vim my_case/case.py[/cyan]             Edit parameters\n"
        "  [green]3.[/green] [cyan]./mfc.sh validate my_case/case.py[/cyan]  Check for errors\n"
        "  [green]4.[/green] [cyan]./mfc.sh build -j $(nproc)[/cyan]      Build MFC\n"
        "  [green]5.[/green] [cyan]./mfc.sh run my_case/case.py[/cyan]    Run simulation",
        box=box.ROUNDED,
        border_style="green",
        padding=(1, 2)
    ))
    cons.print()

    # Footer
    cons.print("[dim]Run [cyan]./mfc.sh <command> --help[/cyan] for detailed options[/dim]")
    cons.print("[dim]Run [cyan]./mfc.sh interactive[/cyan] for guided menu[/dim]")
    cons.print()


def print_command_help(command: str):
    """Print detailed help for a specific command."""
    if command not in COMMANDS:
        cons.print(f"[red]Unknown command: {command}[/red]")
        return

    cmd = COMMANDS[command]

    cons.print()
    cons.print(f"[bold cyan]{command}[/bold cyan] - {cmd['description']}")
    cons.print()

    if cmd.get("examples"):
        cons.print("[bold]Examples:[/bold]")
        for example, desc in cmd["examples"]:
            cons.print(f"  [green]{example}[/green]")
            cons.print(f"    [dim]{desc}[/dim]")
        cons.print()


# =============================================================================
# CONTEXTUAL TIPS
# =============================================================================

class Tips:
    """Contextual tips shown after various events."""

    @staticmethod
    def after_build_failure():
        """Show tips after a build failure."""
        cons.print()
        cons.raw.print(Panel(
            "[bold yellow]Troubleshooting Tips[/bold yellow]\n\n"
            "  [cyan]1.[/cyan] Run with [green]--debug-log[/green] to see detailed output\n"
            "  [cyan]2.[/cyan] Check [green]docs/documentation/troubleshooting.md[/green]\n"
            "  [cyan]3.[/cyan] Ensure required modules are loaded: [green]source ./mfc.sh load -c <cluster> -m <mode>[/green]\n"
            "  [cyan]4.[/cyan] Try [green]./mfc.sh clean[/green] and rebuild",
            box=box.ROUNDED,
            border_style="yellow",
            padding=(0, 2)
        ))

    @staticmethod
    def after_case_error(case_path: str = None):
        """Show tips after a case file error."""
        msg = "[bold yellow]Tip[/bold yellow]\n\n"
        if case_path:
            msg += f"  Run [green]./mfc.sh validate {case_path}[/green] to check your case file for errors"
        else:
            msg += "  Run [green]./mfc.sh validate <case.py>[/green] to check your case file for errors"

        cons.print()
        cons.raw.print(Panel(msg, box=box.ROUNDED, border_style="yellow", padding=(0, 2)))

    @staticmethod
    def after_test_failure(failed_uuids: list = None):
        """Show tips after test failures."""
        lines = [
            "[bold yellow]Next Steps[/bold yellow]\n",
            "  [cyan]1.[/cyan] Check individual test output in [green]tests/<UUID>/[/green]",
            "  [cyan]2.[/cyan] Run specific test: [green]./mfc.sh test --only <UUID>[/green]",
            "  [cyan]3.[/cyan] Update golden files (if changes are intentional): [green]./mfc.sh test --generate[/green]",
        ]

        if failed_uuids and len(failed_uuids) <= 3:
            lines.append("")
            lines.append("  [bold]Failed tests:[/bold]")
            for uuid in failed_uuids:
                lines.append(f"    [red]•[/red] {uuid}")

        cons.print()
        cons.raw.print(Panel("\n".join(lines), box=box.ROUNDED, border_style="yellow", padding=(0, 2)))

    @staticmethod
    def after_run_failure():
        """Show tips after a run failure."""
        cons.print()
        cons.raw.print(Panel(
            "[bold yellow]Troubleshooting Tips[/bold yellow]\n\n"
            "  [cyan]1.[/cyan] Validate your case: [green]./mfc.sh validate case.py[/green]\n"
            "  [cyan]2.[/cyan] Check the output in [green]<case_dir>/[/green]\n"
            "  [cyan]3.[/cyan] Run with [green]--debug-log[/green] for more details\n"
            "  [cyan]4.[/cyan] Check MFC documentation: [green]docs/[/green]",
            box=box.ROUNDED,
            border_style="yellow",
            padding=(0, 2)
        ))

    @staticmethod
    def suggest_validate():
        """Generic suggestion to use validate."""
        cons.print()
        cons.print("[dim]Tip: Run [cyan]./mfc.sh validate case.py[/cyan] to check for errors before running[/dim]")


# =============================================================================
# ONBOARDING FOR NEW USERS
# =============================================================================

def is_first_time_user() -> bool:
    """Check if this is a first-time user (no build directory)."""
    build_dir = os.path.join(MFC_ROOT_DIR, "build")
    return not os.path.exists(build_dir)


def print_welcome():
    """Print welcome message for new users."""
    cons.print()
    cons.raw.print(Panel(
        "[bold cyan]Welcome to MFC![/bold cyan]\n\n"
        "It looks like this is your first time using MFC. Here's how to get started:\n\n"
        "  [green]1.[/green] [bold]Load environment[/bold] (HPC clusters):\n"
        "     [cyan]source ./mfc.sh load -c <cluster> -m <mode>[/cyan]\n"
        "     Example: [dim]source ./mfc.sh load -c p -m g[/dim] (Phoenix, GPU)\n\n"
        "  [green]2.[/green] [bold]Create a new case[/bold]:\n"
        "     [cyan]./mfc.sh init my_first_case[/cyan]\n\n"
        "  [green]3.[/green] [bold]Build MFC[/bold]:\n"
        "     [cyan]./mfc.sh build -j $(nproc)[/cyan]\n\n"
        "  [green]4.[/green] [bold]Run your simulation[/bold]:\n"
        "     [cyan]./mfc.sh run my_first_case/case.py[/cyan]\n\n"
        "[bold yellow]Optional:[/bold yellow] Enable tab completion for your shell:\n"
        "     [cyan]./mfc.sh completion install[/cyan]\n\n"
        "[dim]Run [cyan]./mfc.sh --help[/cyan] for all available commands[/dim]\n"
        "[dim]Run [cyan]./mfc.sh interactive[/cyan] for a guided menu[/dim]",
        title="[bold]Getting Started[/bold]",
        box=box.DOUBLE,
        border_style="cyan",
        padding=(1, 2)
    ))
    cons.print()


# =============================================================================
# INTERACTIVE MODE
# =============================================================================

def interactive_mode():
    """Run interactive menu-driven interface."""

    while True:
        cons.print()
        cons.raw.print(Panel(
            "[bold cyan]MFC Interactive Mode[/bold cyan]",
            box=box.ROUNDED,
            padding=(0, 2)
        ))
        cons.print()

        # Menu options
        options = [
            ("1", "Create a new case", "init"),
            ("2", "Validate a case file", "validate"),
            ("3", "Build MFC", "build"),
            ("4", "Run a simulation", "run"),
            ("5", "Run tests", "test"),
            ("6", "Clean build files", "clean"),
            ("7", "Show help", "help"),
            ("q", "Quit", None),
        ]

        for key, label, _ in options:
            if key == "q":
                cons.print(f"  [red]{key}[/red]) {label}")
            else:
                cons.print(f"  [green]{key}[/green]) {label}")

        cons.print()
        choice = Prompt.ask("[bold]Select an option[/bold]", choices=[o[0] for o in options], default="q")

        if choice == "q":
            cons.print("[dim]Goodbye![/dim]")
            break

        if choice == "7":
            print_help()
            continue

        # Get the command for the selected option
        cmd = next((o[2] for o in options if o[0] == choice), None)
        if cmd is None:
            continue

        cons.print()

        # Dispatch to handler
        handlers = {
            "init": _interactive_init,
            "validate": _interactive_validate,
            "build": _interactive_build,
            "run": _interactive_run,
            "test": _interactive_test,
            "clean": _interactive_clean,
        }
        if cmd in handlers:
            handlers[cmd]()


def _interactive_init():
    """Interactive case creation."""
    cons.print("[bold]Create a New Case[/bold]")
    cons.print()

    # Show templates
    cons.print("Available templates: [cyan]1D_minimal[/cyan], [cyan]2D_minimal[/cyan], [cyan]3D_minimal[/cyan]")
    cons.print("[dim]Or use 'example:<name>' to copy from examples[/dim]")
    cons.print()

    name = Prompt.ask("Case name", default="my_case")
    template = Prompt.ask("Template", default="1D_minimal")

    cmd = f"./mfc.sh init {name} -t {template}"
    cons.print()
    cons.print(f"[dim]Running: {cmd}[/dim]")
    cons.print()
    os.system(cmd)


def _interactive_validate():
    """Interactive case validation."""
    cons.print("[bold]Validate a Case File[/bold]")
    cons.print()

    path = Prompt.ask("Path to case.py")

    cmd = f"./mfc.sh validate {path}"
    cons.print()
    cons.print(f"[dim]Running: {cmd}[/dim]")
    cons.print()
    os.system(cmd)


def _interactive_build():
    """Interactive build."""
    cons.print("[bold]Build MFC[/bold]")
    cons.print()

    jobs = Prompt.ask("Number of parallel jobs", default="4")
    gpu = Prompt.ask("Enable GPU support?", choices=["y", "n"], default="n")

    cmd = f"./mfc.sh build -j {jobs}"
    if gpu == "y":
        cmd += " --gpu"

    cons.print()
    cons.print(f"[dim]Running: {cmd}[/dim]")
    cons.print()
    os.system(cmd)


def _interactive_run():
    """Interactive run."""
    cons.print("[bold]Run a Simulation[/bold]")
    cons.print()

    path = Prompt.ask("Path to case.py")
    ranks = Prompt.ask("Number of MPI ranks", default="1")

    cmd = f"./mfc.sh run {path} -n {ranks}"

    cons.print()
    cons.print(f"[dim]Running: {cmd}[/dim]")
    cons.print()
    os.system(cmd)


def _interactive_test():
    """Interactive test."""
    cons.print("[bold]Run Tests[/bold]")
    cons.print()

    jobs = Prompt.ask("Number of parallel jobs", default="4")

    cmd = f"./mfc.sh test -j {jobs}"

    cons.print()
    cons.print(f"[dim]Running: {cmd}[/dim]")
    cons.print()
    os.system(cmd)


def _interactive_clean():
    """Interactive clean."""
    cons.print("[bold]Clean Build Files[/bold]")
    cons.print()

    confirm = Prompt.ask("Are you sure you want to clean all build files?", choices=["y", "n"], default="n")

    if confirm == "y":
        cmd = "./mfc.sh clean"
        cons.print()
        cons.print(f"[dim]Running: {cmd}[/dim]")
        cons.print()
        os.system(cmd)
    else:
        cons.print("[dim]Cancelled[/dim]")
