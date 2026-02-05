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
import subprocess
import re

from rich.panel import Panel
from rich.table import Table
from rich.prompt import Prompt
from rich.markdown import Markdown
from rich import box

from .printer import cons
from .common import MFC_ROOT_DIR

# Import command definitions from CLI schema (SINGLE SOURCE OF TRUTH)
from .cli.commands import COMMANDS


# =============================================================================
# DYNAMIC CLUSTER HELP GENERATION
# =============================================================================

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
            for line in f:
                line = line.strip()
                # Skip comments and empty lines
                if not line or line.startswith("#"):
                    continue
                # Skip lines with -all, -cpu, -gpu (module definitions)
                if "-all" in line or "-cpu" in line or "-gpu" in line:
                    continue

                # Parse cluster definition lines: "slug     System Name"
                match = re.match(r'^(\S+)\s+(.+)$', line)
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
            return full_name[len(prefix) + 1:]
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
        entries = [
            f"[cyan]{slug}[/cyan]={_get_cluster_short_name(slug, name)}"
            for slug, name in org_clusters[org]
        ]
        color = ORG_COLORS.get(org, 'yellow')
        cluster_lines.append(f"  [{color}]{org}:[/{color}]    " + "  ".join(entries))

    # Handle "Other" if any
    if org_clusters.get("Other"):
        entries = [f"[cyan]{slug}[/cyan]={name}" for slug, name in org_clusters["Other"]]
        cluster_lines.append(f"  [yellow]Other:[/yellow]    " + "  ".join(entries))

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


# =============================================================================
# MARKDOWN-BASED HELP (Single source of truth from docs/)
# =============================================================================

# Mapping of help topics to their source markdown files and optional section
# Format: {"topic": ("file_path", "section_heading" or None for full file)}
MARKDOWN_HELP_FILES = {
    "debugging": ("docs/documentation/troubleshooting.md", None),  # Full file
    "gpu": ("docs/documentation/running.md", "Running on GPUs"),  # Section only
    "batch": ("docs/documentation/running.md", "Batch Execution"),  # Section only
    "performance": ("docs/documentation/expectedPerformance.md", "Achieving Maximum Performance"),
}


def _extract_markdown_section(content: str, section_heading: str) -> str:
    """Extract a specific section from markdown content.

    Extracts from the given heading until the next heading of same or higher level,
    or until a horizontal rule (---).
    """
    # Find the section heading (## or ###)
    # Note: In f-strings, literal braces must be doubled: {{1,3}} -> {1,3}
    pattern = rf'^(#{{1,3}})\s+{re.escape(section_heading)}\s*$'
    match = re.search(pattern, content, re.MULTILINE)
    if not match:
        return None

    start_pos = match.end()

    # Find the end: horizontal rule (---) which separates major sections
    # Note: We use --- instead of heading detection because shell comments
    # inside code blocks (# comment) look like markdown headings to regex
    end_pattern = r'^---'
    end_match = re.search(end_pattern, content[start_pos:], re.MULTILINE)

    if end_match:
        section = content[start_pos:start_pos + end_match.start()]
    else:
        section = content[start_pos:]

    return section.strip()


def _load_markdown_help(topic: str) -> str:
    """Load help content from a markdown file.

    Can load full file or extract a specific section.
    Strips Doxygen-specific syntax and returns clean markdown.
    """
    if topic not in MARKDOWN_HELP_FILES:
        return None

    file_path, section = MARKDOWN_HELP_FILES[topic]
    filepath = os.path.join(MFC_ROOT_DIR, file_path)

    try:
        with open(filepath, "r", encoding="utf-8") as f:
            content = f.read()
    except FileNotFoundError:
        return None

    # Extract section if specified
    if section:
        content = _extract_markdown_section(content, section)
        if content is None:
            return None

    # Strip Doxygen-specific syntax
    # Remove @page directives
    content = re.sub(r'^@page\s+\S+\s+.*$', '', content, flags=re.MULTILINE)
    # Remove @ref, @see directives (but keep the text after them readable)
    content = re.sub(r'@(ref|see)\s+"([^"]+)"', r'\2', content)  # @ref "Text" -> Text
    content = re.sub(r'@(ref|see)\s+(\S+)', '', content)  # @ref name -> (remove)
    # Clean up any resulting empty lines at the start
    content = content.lstrip('\n')

    return content


def _generate_markdown_help(topic: str):
    """Generate a function that loads markdown help for a topic."""
    def loader():
        return _load_markdown_help(topic)
    return loader


# =============================================================================
# HELP TOPICS
# =============================================================================

HELP_TOPICS = {
    "gpu": {
        "title": "Running on GPUs",
        # Content loaded from docs/documentation/running.md "Running on GPUs" section
        "content": _generate_markdown_help("gpu"),
        "markdown": True,
    },
    "clusters": {
        "title": "Cluster Configuration",
        # Content is generated dynamically from toolchain/modules file
        "content": _generate_clusters_content,
    },
    "batch": {
        "title": "Batch Job Submission",
        # Content loaded from docs/documentation/running.md "Batch Execution" section
        "content": _generate_markdown_help("batch"),
        "markdown": True,
    },
    "debugging": {
        "title": "Debugging & Troubleshooting",
        # Content loaded from docs/documentation/troubleshooting.md
        "content": _generate_markdown_help("debugging"),
        "markdown": True,
    },
    "performance": {
        "title": "Performance Optimization",
        # Content loaded from docs/documentation/expectedPerformance.md "Achieving Maximum Performance" section
        "content": _generate_markdown_help("performance"),
        "markdown": True,
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
    # Support callable content for dynamic generation
    content = topic_info["content"]
    if callable(content):
        content = content()

    if content is None:
        cons.print(f"[red]Could not load help for topic: {topic}[/red]")
        return

    cons.print()

    # Check if content should be rendered as markdown
    if topic_info.get("markdown", False):
        # Render markdown content directly (no panel - markdown has its own formatting)
        cons.print(f"[bold cyan]{topic_info['title']}[/bold cyan]")
        cons.print()
        cons.raw.print(Markdown(content))
    else:
        # Render as Rich markup in a panel
        cons.raw.print(Panel(
            content,
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

def _truncate_desc(desc: str, max_len: int = 50) -> str:
    """Truncate description to fit compact display."""
    if len(desc) <= max_len:
        return desc
    return desc[:max_len-3] + "..."


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
    primary = ["build", "run", "test", "validate", "new", "clean"]
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
            "  [cyan]1.[/cyan] Rebuild with [green]--debug[/green] for debug compiler flags and verbose output\n"
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
            "  [cyan]3.[/cyan] Rebuild with [green]--debug[/green] for debug compiler flags\n"
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


def _run_mfc_command(args: list):
    """Run an MFC command safely using subprocess."""
    cmd_str = " ".join(args)
    cons.print()
    cons.print(f"[dim]Running: {cmd_str}[/dim]")
    cons.print()
    try:
        subprocess.run(args, check=False)
    except FileNotFoundError:
        cons.print(f"[red]Command not found: {args[0]}[/red]")


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

    _run_mfc_command(["./mfc.sh", "new", name, "-t", template])


def _interactive_validate():
    """Interactive case validation."""
    cons.print("[bold]Validate a Case File[/bold]")
    cons.print()

    path = Prompt.ask("Path to case.py")

    _run_mfc_command(["./mfc.sh", "validate", path])


def _interactive_build():
    """Interactive build."""
    cons.print("[bold]Build MFC[/bold]")
    cons.print()

    jobs = Prompt.ask("Number of parallel jobs", default="4")
    gpu = Prompt.ask("Enable GPU support?", choices=["y", "n"], default="n")

    args = ["./mfc.sh", "build", "-j", jobs]
    if gpu == "y":
        args.append("--gpu")

    _run_mfc_command(args)


def _interactive_run():
    """Interactive run."""
    cons.print("[bold]Run a Simulation[/bold]")
    cons.print()

    path = Prompt.ask("Path to case.py")
    ranks = Prompt.ask("Number of MPI ranks", default="1")

    _run_mfc_command(["./mfc.sh", "run", path, "-n", ranks])


def _interactive_test():
    """Interactive test."""
    cons.print("[bold]Run Tests[/bold]")
    cons.print()

    jobs = Prompt.ask("Number of parallel jobs", default="4")

    _run_mfc_command(["./mfc.sh", "test", "-j", jobs])


def _interactive_clean():
    """Interactive clean."""
    cons.print("[bold]Clean Build Files[/bold]")
    cons.print()

    confirm = Prompt.ask("Are you sure you want to clean all build files?", choices=["y", "n"], default="n")

    if confirm == "y":
        _run_mfc_command(["./mfc.sh", "clean"])
    else:
        cons.print("[dim]Cancelled[/dim]")
