"""
Generate completion scripts and documentation from CLI schema.

This module regenerates all derived files from the single source of truth
in cli/commands.py. Run `./mfc.sh generate` after modifying commands.
"""

from pathlib import Path

from .printer import cons
from .common import MFC_ROOT_DIR
from .state import ARG
from .cli.commands import MFC_CLI_SCHEMA
from .cli.completion_gen import generate_bash_completion, generate_zsh_completion


def generate():
    """Regenerate completion scripts from CLI schema."""
    check_mode = ARG("check")

    completions_dir = Path(MFC_ROOT_DIR) / "toolchain" / "completions"

    # Generate bash completion
    bash_content = generate_bash_completion(MFC_CLI_SCHEMA)
    bash_path = completions_dir / "mfc.bash"

    if check_mode:
        if not bash_path.exists():
            cons.print(f"[red]ERROR:[/red] {bash_path} does not exist")
            exit(1)
        existing = bash_path.read_text()
        if existing != bash_content:
            cons.print(f"[red]ERROR:[/red] {bash_path} is out of date")
            cons.print("[yellow]Run ./mfc.sh generate to update[/yellow]")
            exit(1)
        cons.print(f"[green]OK[/green] {bash_path.name} is up to date")
    else:
        bash_path.write_text(bash_content)
        cons.print(f"[green]Generated[/green] {bash_path}")

    # Generate zsh completion
    zsh_content = generate_zsh_completion(MFC_CLI_SCHEMA)
    zsh_path = completions_dir / "_mfc"

    if check_mode:
        if not zsh_path.exists():
            cons.print(f"[red]ERROR:[/red] {zsh_path} does not exist")
            exit(1)
        existing = zsh_path.read_text()
        if existing != zsh_content:
            cons.print(f"[red]ERROR:[/red] {zsh_path} is out of date")
            cons.print("[yellow]Run ./mfc.sh generate to update[/yellow]")
            exit(1)
        cons.print(f"[green]OK[/green] {zsh_path.name} is up to date")
    else:
        zsh_path.write_text(zsh_content)
        cons.print(f"[green]Generated[/green] {zsh_path}")

    if not check_mode:
        cons.print()
        cons.print("[bold]Completion scripts regenerated from cli/commands.py[/bold]")
        cons.print("[dim]Commit these files to keep completions in sync[/dim]")
