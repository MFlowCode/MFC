"""
Help output and onboarding for the MFC toolchain.

This module provides:
- Enhanced help output with Rich formatting
- The cluster table behind ``./mfc.sh help``, derived from toolchain/modules
- Onboarding for new users
"""

import os
import re

from rich import box
from rich.panel import Panel

# Import command definitions from CLI schema (SINGLE SOURCE OF TRUTH)
from .cli.commands import COMMANDS
from .common import MFC_ROOT_DIR
from .printer import cons

# DYNAMIC CLUSTER HELP GENERATION

# Organization mapping based on system name prefixes and known clusters
CLUSTER_ORGS = {
    "OLCF": "ORNL",
    "LLNL": "LLNL",
    "PSC": "ACCESS",
    "SDSC": "ACCESS",
    "NCSA": "ACCESS",
    "GT": "Georgia Tech",
    "Brown": "Brown",
    "DoD": "DoD",
    "Richardson": "Caltech",
    "hipergator": "Florida",
    "CSCS": "CSCS",
}

# Explicit slug-to-org overrides (for cases where modules file naming is inconsistent)
SLUG_ORG_OVERRIDE = {
    "tuo": "LLNL",  # Tuolumne is at LLNL, not ORNL (modules file says "OLCF" incorrectly)
}

# Display name overrides for clusters
SLUG_NAME_OVERRIDE = {
    "h": "HiPerGator",  # Proper capitalization
}

# Display order and colors for organizations
ORG_ORDER = ["ORNL", "LLNL", "ACCESS", "Georgia Tech", "Caltech", "Brown", "DoD", "Florida", "CSCS"]
ORG_COLORS = {
    "ORNL": "yellow",
    "LLNL": "yellow",
    "ACCESS": "yellow",
    "Georgia Tech": "yellow",
    "Caltech": "yellow",
    "Brown": "yellow",
    "DoD": "yellow",
    "Florida": "yellow",
}


def _parse_modules_file():
    """Parse the modules file to extract cluster information.

    Returns a dict: {slug: {"name": full_name, "org": organization}}
    """
    modules_path = os.path.join(MFC_ROOT_DIR, "toolchain", "modules")
    clusters = {}

    try:
        with open(modules_path, "r", encoding="utf-8") as f:
            for raw_line in f:
                line = raw_line.strip()
                # Skip comments and empty lines
                if not line or line.startswith("#"):
                    continue
                # Skip lines with -all, -cpu, -gpu (module definitions)
                if "-all" in line or "-cpu" in line or "-gpu" in line:
                    continue

                # Parse cluster definition lines: "slug     System Name"
                match = re.match(r"^(\S+)\s+(.+)$", line)
                if match:
                    slug = match.group(1)
                    full_name = match.group(2).strip()

                    # Check for explicit org override first
                    if slug in SLUG_ORG_OVERRIDE:
                        org = SLUG_ORG_OVERRIDE[slug]
                    else:
                        # Determine organization from name
                        org = "Other"
                        for prefix, org_name in CLUSTER_ORGS.items():
                            if prefix in full_name or full_name.lower() == prefix.lower():
                                org = org_name
                                break

                    clusters[slug] = {"name": full_name, "org": org}
    except FileNotFoundError:
        # Fallback if modules file not found
        pass

    return clusters


def _get_cluster_short_name(slug, full_name):
    """Get display name for a cluster, with overrides and prefix stripping."""
    if slug in SLUG_NAME_OVERRIDE:
        return SLUG_NAME_OVERRIDE[slug]
    # Strip org prefix if present
    for prefix in CLUSTER_ORGS:
        if full_name.startswith(prefix + " "):
            return full_name[len(prefix) + 1 :]
    return full_name


def _generate_clusters_content():
    """Generate the clusters help content dynamically from modules file."""
    clusters = _parse_modules_file()

    # Group clusters by organization
    org_clusters = {org: [] for org in ORG_ORDER}
    org_clusters["Other"] = []

    for slug, info in clusters.items():
        org = info["org"]
        if org not in org_clusters:
            org_clusters["Other"].append((slug, info["name"]))
        else:
            org_clusters[org].append((slug, info["name"]))

    # Build the cluster list section
    cluster_lines = []
    for org in ORG_ORDER:
        if not org_clusters.get(org):
            continue
        # Format: "  [yellow]ORG:[/yellow]  [cyan]slug[/cyan]=Name  [cyan]slug2[/cyan]=Name2"
        entries = [f"[cyan]{slug}[/cyan]={_get_cluster_short_name(slug, name)}" for slug, name in org_clusters[org]]
        color = ORG_COLORS.get(org, "yellow")
        cluster_lines.append(f"  [{color}]{org}:[/{color}]    " + "  ".join(entries))

    # Handle "Other" if any
    if org_clusters.get("Other"):
        entries = [f"[cyan]{slug}[/cyan]={name}" for slug, name in org_clusters["Other"]]
        cluster_lines.append("  [yellow]Other:[/yellow]    " + "  ".join(entries))

    cluster_list = "\n".join(cluster_lines) if cluster_lines else "  [dim]No clusters found in modules file[/dim]"

    # Return full help content with dynamic cluster list
    return f"""\
[bold cyan]Supported HPC Clusters[/bold cyan]

MFC includes pre-configured module sets for many clusters.

[bold]Loading Cluster Modules:[/bold]
  [green]source ./mfc.sh load -c <cluster> -m <mode>[/green]

[bold]Available Clusters:[/bold]
{cluster_list}

[bold]Modes:[/bold]
  [cyan]c[/cyan] or [cyan]cpu[/cyan] - CPU only
  [cyan]g[/cyan] or [cyan]gpu[/cyan] - GPU enabled

[bold]Examples:[/bold]
  [green]source ./mfc.sh load -c p -m g[/green]     Phoenix with GPU
  [green]source ./mfc.sh load -c f -m g[/green]     Frontier with GPU (AMD MI250X)
  [green]source ./mfc.sh load -c d -m c[/green]     Delta CPU-only

[bold]Custom Clusters:[/bold]
  For unlisted clusters, manually load:
  • Fortran compiler (gfortran, nvfortran, amdflang, etc.)
  • MPI implementation (OpenMPI, MPICH, Cray-MPICH)
  • CMake 3.18+, Python 3.11+"""


def print_clusters_help():
    """Print the cluster configuration table behind ``./mfc.sh help``."""
    cons.print()
    cons.raw.print(Panel(_generate_clusters_content(), title="[bold]Cluster Configuration[/bold]", box=box.ROUNDED, padding=(1, 2)))
    cons.print()


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
    secondary = ["params", "load", "help"]
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
            "[dim]Run [cyan]./mfc.sh --help[/cyan] for all available commands[/dim]\n"
            "[dim]Run [cyan]./mfc.sh help[/cyan] for the list of supported clusters[/dim]",
            title="[bold]Getting Started[/bold]",
            box=box.DOUBLE,
            border_style="cyan",
            padding=(1, 2),
        )
    )
    cons.print()
