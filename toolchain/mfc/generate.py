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
from .cli.docs_gen import generate_cli_reference


def _check_or_write(path: Path, content: str, check_mode: bool) -> bool:
    """Check if file is up to date or write new content. Returns True on success."""
    if check_mode:
        if not path.exists():
            cons.print(f"[red]ERROR:[/red] {path} does not exist")
            return False
        if path.read_text() != content:
            cons.print(f"[red]ERROR:[/red] {path} is out of date")
            cons.print("[yellow]Run ./mfc.sh generate to update[/yellow]")
            return False
        cons.print(f"[green]OK[/green] {path.name} is up to date")
    else:
        path.write_text(content)
        cons.print(f"[green]Generated[/green] {path}")
    return True


def generate():
    """Regenerate completion scripts and optionally JSON schema."""
    check_mode = ARG("check")
    json_schema_mode = ARG("json_schema")

    # If only generating JSON schema, do that and return
    if json_schema_mode:
        _generate_json_schema()
        return

    completions_dir = Path(MFC_ROOT_DIR) / "toolchain" / "completions"
    docs_dir = Path(MFC_ROOT_DIR) / "docs"
    docs_dir.mkdir(exist_ok=True)

    # Generate and check/write all files
    files = [
        (completions_dir / "mfc.bash", generate_bash_completion(MFC_CLI_SCHEMA)),
        (completions_dir / "_mfc", generate_zsh_completion(MFC_CLI_SCHEMA)),
        (docs_dir / "cli-reference.md", generate_cli_reference(MFC_CLI_SCHEMA)),
    ]

    for path, content in files:
        if not _check_or_write(path, content, check_mode):
            exit(1)

    if not check_mode:
        cons.print()
        cons.print("[bold]Files regenerated from cli/commands.py[/bold]")
        cons.print("[dim]Commit these files to keep them in sync[/dim]")


def _generate_json_schema():
    """Generate JSON Schema for IDE auto-completion."""
    import json
    from .params.generators.json_schema_gen import generate_json_schema, get_schema_stats
    from .ide import update_vscode_settings

    schema = generate_json_schema(include_descriptions=True)
    output_path = Path(MFC_ROOT_DIR) / "toolchain" / "mfc-case-schema.json"

    with open(output_path, 'w') as f:
        json.dump(schema, f, indent=2)

    # Update VS Code settings
    update_vscode_settings()

    stats = get_schema_stats()

    cons.print(f"[green]Generated[/green] {output_path}")
    cons.print()
    cons.print(f"[bold]JSON Schema Statistics:[/bold]")
    cons.print(f"  Total parameters: {stats['total_params']}")
    cons.print(f"  With constraints: {stats['with_constraints']}")
    cons.print(f"  With descriptions: {stats['with_descriptions']}")
    cons.print()
    cons.print("[bold]Usage:[/bold]")
    cons.print("  VS Code: .vscode/settings.json auto-configured on startup")
    cons.print("  Auto-completion works for case.json and case.yaml files")
    cons.print("  For Python case files, use ./mfc.sh params <name> for parameter info")
