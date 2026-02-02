"""
User guide, help, tips, and onboarding for MFC toolchain.

This module provides:
- Enhanced help output with Rich formatting
- Contextual tips after errors/failures
- Interactive mode with menu
- Onboarding for new users
- Topic-based help system
"""

import os

from rich.panel import Panel
from rich.table import Table
from rich.prompt import Prompt
from rich import box

from .printer import cons
from .common import MFC_ROOT_DIR

# Import command definitions from CLI schema (SINGLE SOURCE OF TRUTH)
from .cli.commands import COMMANDS


# =============================================================================
# HELP TOPICS
# =============================================================================

HELP_TOPICS = {
    "gpu": {
        "title": "GPU Configuration",
        "content": """\
[bold cyan]GPU Support in MFC[/bold cyan]

MFC supports GPU acceleration via OpenACC or OpenMP offloading.

[bold]Building with GPU support:[/bold]
  [green]./mfc.sh build --gpu[/green]           Enable OpenACC (default)
  [green]./mfc.sh build --gpu acc[/green]       Explicitly use OpenACC
  [green]./mfc.sh build --gpu mp[/green]        Use OpenMP offloading

[bold]Running on GPUs:[/bold]
  [green]./mfc.sh run case.py --gpu[/green]     Run with GPU support
  [green]./mfc.sh run case.py -g 0 1[/green]    Use specific GPU IDs

[bold]GPU Profiling:[/bold]
  [green]./mfc.sh run case.py --ncu[/green]     NVIDIA Nsight Compute
  [green]./mfc.sh run case.py --nsys[/green]    NVIDIA Nsight Systems
  [green]./mfc.sh run case.py --rcu[/green]     AMD ROCm rocprof-compute
  [green]./mfc.sh run case.py --rsys[/green]    AMD ROCm rocprof-systems

[bold]Environment Setup:[/bold]
  Most HPC systems require loading GPU modules first:
  [cyan]source ./mfc.sh load -c <cluster> -m g[/cyan]

[bold]Supported Compilers:[/bold]
  [yellow]NVIDIA GPUs:[/yellow]
    • NVHPC (nvfortran) - OpenACC or OpenMP
  [yellow]AMD GPUs:[/yellow]
    • Cray (ftn) - OpenACC or OpenMP
    • AMD (amdflang) - OpenMP only

[bold]Troubleshooting:[/bold]
  • Ensure CUDA/ROCm toolkit is in PATH
  • Check GPU visibility: [cyan]nvidia-smi[/cyan] or [cyan]rocm-smi[/cyan]
  • Use [cyan]--debug-log[/cyan] for detailed build output"""
    },
    "clusters": {
        "title": "Cluster Configuration",
        "content": """\
[bold cyan]Supported HPC Clusters[/bold cyan]

MFC includes pre-configured module sets for many clusters.

[bold]Loading Cluster Modules:[/bold]
  [green]source ./mfc.sh load -c <cluster> -m <mode>[/green]

[bold]Available Clusters:[/bold]
  [yellow]ORNL:[/yellow]      [cyan]a[/cyan]=Ascent  [cyan]f[/cyan]=Frontier  [cyan]famd[/cyan]=Frontier AMD  [cyan]s[/cyan]=Summit  [cyan]w[/cyan]=Wombat
  [yellow]LLNL:[/yellow]      [cyan]tuo[/cyan]=Tuolumne
  [yellow]ACCESS:[/yellow]    [cyan]b[/cyan]=Bridges2  [cyan]e[/cyan]=Expanse  [cyan]d[/cyan]=Delta  [cyan]dai[/cyan]=DeltaAI
  [yellow]Georgia Tech:[/yellow] [cyan]p[/cyan]=Phoenix
  [yellow]Caltech:[/yellow]   [cyan]r[/cyan]=Richardson
  [yellow]Brown:[/yellow]     [cyan]o[/cyan]=Oscar
  [yellow]DoD:[/yellow]       [cyan]cc[/cyan]=Carpenter Cray  [cyan]c[/cyan]=Carpenter GNU  [cyan]n[/cyan]=Nautilus
  [yellow]Florida:[/yellow]   [cyan]h[/cyan]=HiPerGator

[bold]Modes:[/bold]
  [cyan]c[/cyan] or [cyan]cpu[/cyan] - CPU only
  [cyan]g[/cyan] or [cyan]gpu[/cyan] - GPU enabled

[bold]Examples:[/bold]
  [green]source ./mfc.sh load -c p -m g[/green]     Phoenix with GPU
  [green]source ./mfc.sh load -c f -m g[/green]     Frontier with GPU (AMD MI250X)
  [green]source ./mfc.sh load -c s -m g[/green]     Summit with GPU (NVIDIA V100)
  [green]source ./mfc.sh load -c d -m c[/green]     Delta CPU-only

[bold]Custom Clusters:[/bold]
  For unlisted clusters, manually load:
  • Fortran compiler (gfortran, nvfortran, amdflang, etc.)
  • MPI implementation (OpenMPI, MPICH, Cray-MPICH)
  • CMake 3.18+, Python 3.11+"""
    },
    "batch": {
        "title": "Batch Job Submission",
        "content": """\
[bold cyan]Submitting Batch Jobs[/bold cyan]

Use [green]-e batch[/green] to submit jobs to a scheduler instead of running interactively.

[bold]Basic Batch Submission:[/bold]
  [green]./mfc.sh run case.py -e batch[/green]

[bold]Common Options:[/bold]
  [cyan]-N, --nodes[/cyan]           Number of nodes
  [cyan]-n, --tasks-per-node[/cyan]  MPI ranks per node
  [cyan]-w, --walltime[/cyan]        Time limit (HH:MM:SS)
  [cyan]-a, --account[/cyan]         Account/allocation to charge
  [cyan]-p, --partition[/cyan]       Queue/partition name
  [cyan]-q, --qos[/cyan]             Quality of service
  [cyan]-@, --email[/cyan]           Email for job notifications
  [cyan]-#, --name[/cyan]            Job name (default: MFC)
  [cyan]-c, --computer[/cyan]        Submission template (e.g., phoenix, frontier)

[bold]Examples:[/bold]
  [green]./mfc.sh run case.py -e batch -N 4 -n 8 -w 02:00:00[/green]
    4 nodes, 8 ranks/node, 2 hour limit

  [green]./mfc.sh run case.py -e batch -a myproject -p gpu -@ user@email.com[/green]
    Submit to 'gpu' partition with email notifications

  [green]./mfc.sh run case.py -e batch -c frontier[/green]
    Use Frontier-specific submission template

[bold]Dry Run:[/bold]
  [green]./mfc.sh run case.py -e batch --dry-run[/green]
  Shows the generated batch script without submitting.

[bold]Wait for Completion:[/bold]
  [green]./mfc.sh run case.py -e batch --wait[/green]
  Blocks until the job finishes."""
    },
    "debugging": {
        "title": "Debugging & Troubleshooting",
        "content": """\
[bold cyan]Debugging MFC[/bold cyan]

[bold]Verbosity Levels:[/bold]
  [green]-v[/green]       Basic verbose output
  [green]-vv[/green]      More detailed output
  [green]-vvv[/green]     Maximum verbosity
  [green]-d[/green]       Debug log (writes to file)

[bold]Debug Builds:[/bold]
  [green]./mfc.sh build --debug[/green]     Debug symbols, reduced optimization
  [green]./mfc.sh build --gcov[/green]      Code coverage instrumentation

[bold]Validating Cases:[/bold]
  [green]./mfc.sh validate case.py[/green]
  Checks case file for syntax errors and constraint violations.

[bold]Common Issues:[/bold]

  [yellow]Build fails with missing MPI:[/yellow]
    → Load cluster modules: [cyan]source ./mfc.sh load -c <cluster> -m <mode>[/cyan]

  [yellow]Build fails with compiler errors:[/yellow]
    → Try [cyan]./mfc.sh clean[/cyan] then rebuild
    → Check compiler compatibility in docs

  [yellow]Tests fail with tolerance errors:[/yellow]
    → May be expected for different compilers/platforms
    → Use [cyan]--generate[/cyan] to update golden files if intentional

  [yellow]Simulation crashes or gives wrong results:[/yellow]
    → Validate case: [cyan]./mfc.sh validate case.py[/cyan]
    → Build with debug: [cyan]./mfc.sh build --debug[/cyan]
    → Check CFL condition and grid resolution

[bold]Getting Help:[/bold]
  • Docs: [cyan]docs/documentation/[/cyan]
  • Issues: [cyan]https://github.com/MFlowCode/MFC/issues[/cyan]"""
    },
}


def print_topic_help(topic: str):
    """Print help for a specific topic."""
    if topic not in HELP_TOPICS:
        cons.print(f"[red]Unknown topic: {topic}[/red]")
        cons.print()
        cons.print("[bold]Available topics:[/bold]")
        for t, info in HELP_TOPICS.items():
            cons.print(f"  [green]{t:12}[/green] {info['title']}")
        cons.print()
        cons.print("[dim]Usage: ./mfc.sh help <topic>[/dim]")
        return

    topic_info = HELP_TOPICS[topic]
    cons.print()
    cons.raw.print(Panel(
        topic_info["content"],
        title=f"[bold]{topic_info['title']}[/bold]",
        box=box.ROUNDED,
        padding=(1, 2)
    ))
    cons.print()


def print_help_topics():
    """Print list of available help topics."""
    cons.print()
    cons.raw.print(Panel(
        "[bold cyan]MFC Help System[/bold cyan]",
        box=box.ROUNDED,
        padding=(0, 2)
    ))
    cons.print()

    table = Table(box=box.SIMPLE, show_header=False, padding=(0, 2))
    table.add_column("Topic", style="green")
    table.add_column("Description")

    for topic, info in HELP_TOPICS.items():
        table.add_row(topic, info["title"])

    cons.raw.print(table)
    cons.print()
    cons.print("[dim]Usage: [cyan]./mfc.sh help <topic>[/cyan][/dim]")
    cons.print("[dim]Example: [cyan]./mfc.sh help gpu[/cyan][/dim]")
    cons.print()


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
    table.add_column("Alias", style="dim", no_wrap=True)
    table.add_column("Description", style="white")

    # Primary commands with aliases
    for cmd in ["build", "run", "test", "validate", "new", "clean"]:
        alias = COMMANDS[cmd].get("alias", "")
        alias_str = alias if alias else ""
        table.add_row(cmd, alias_str, COMMANDS[cmd]["description"])

    table.add_row("", "", "")  # Spacer

    # Secondary commands
    for cmd in ["params", "count", "packer", "load"]:
        table.add_row(f"[dim]{cmd}[/dim]", "", f"[dim]{COMMANDS[cmd]['description']}[/dim]")

    table.add_row("", "", "")  # Spacer
    table.add_row("[dim]help[/dim]", "", "[dim]Show help on a topic (gpu, clusters, batch, debugging)[/dim]")

    cons.raw.print(table)
    cons.print()

    # Quick start
    cons.raw.print(Panel(
        "[bold]Quick Start[/bold]\n\n"
        "  [green]1.[/green] [cyan]./mfc.sh new my_case[/cyan]            Create a new case\n"
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
    cons.print("[dim]Run [cyan]./mfc.sh help <topic>[/cyan] for topic help (gpu, clusters, batch, debugging)[/dim]")
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
    cons.raw.print(Panel(
        f"[bold cyan]{command}[/bold cyan]{alias_str}\n"
        f"[dim]{cmd['description']}[/dim]",
        box=box.ROUNDED,
        padding=(0, 2)
    ))
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
            cons.print(f"  [cyan]{opt:24}[/cyan] {desc}")
        cons.print()
        if show_argparse:
            cons.print("[dim]Run with --help for full option list[/dim]")
            cons.print()

    return True


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
        "     [cyan]./mfc.sh new my_first_case[/cyan]\n\n"
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
            ("1", "Create a new case", "new"),
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
            "new": _interactive_new,
            "validate": _interactive_validate,
            "build": _interactive_build,
            "run": _interactive_run,
            "test": _interactive_test,
            "clean": _interactive_clean,
        }
        if cmd in handlers:
            handlers[cmd]()


def _interactive_new():
    """Interactive case creation."""
    cons.print("[bold]Create a New Case[/bold]")
    cons.print()

    # Show templates
    cons.print("Available templates: [cyan]1D_minimal[/cyan], [cyan]2D_minimal[/cyan], [cyan]3D_minimal[/cyan]")
    cons.print("[dim]Or use 'example:<name>' to copy from examples[/dim]")
    cons.print()

    name = Prompt.ask("Case name", default="my_case")
    template = Prompt.ask("Template", default="1D_minimal")

    cmd = f"./mfc.sh new {name} -t {template}"
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
