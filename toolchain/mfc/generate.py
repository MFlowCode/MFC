"""
Generate completion scripts and documentation from CLI schema.

This module regenerates all derived files from the single source of truth
in cli/commands.py. Run `./mfc.sh generate` after modifying commands.
"""
# pylint: disable=import-outside-toplevel

import json
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


def _constraint_docs(docs_dir: Path) -> list:
    """Generate constraint documentation files."""
    from .gen_case_constraints_docs import main as gen_case_constraints
    from .gen_physics_docs import render as render_physics
    from .params.ast_analyzer import analyze_case_validator

    validator_path = Path(MFC_ROOT_DIR) / "toolchain" / "mfc" / "case_validator.py"
    rules = analyze_case_validator(validator_path)["rules"]

    return [
        (docs_dir / "case_constraints.md", gen_case_constraints(as_string=True)),
        (docs_dir / "physics_constraints.md", render_physics(rules)),
    ]


def generate():
    """Regenerate completion scripts and optionally JSON schema."""
    from .params.generators.json_schema_gen import generate_json_schema
    from .params.generators.docs_gen import generate_parameter_docs

    check_mode = ARG("check")
    json_schema_mode = ARG("json_schema")

    # If only generating JSON schema, do that and return
    if json_schema_mode:
        _generate_json_schema()
        return

    completions_dir = Path(MFC_ROOT_DIR) / "toolchain" / "completions"
    docs_dir = Path(MFC_ROOT_DIR) / "docs" / "documentation"
    completions_dir.mkdir(exist_ok=True)
    docs_dir.mkdir(exist_ok=True)

    # Generate all derived files
    files = [
        (completions_dir / "mfc.bash", generate_bash_completion(MFC_CLI_SCHEMA)),
        (completions_dir / "_mfc", generate_zsh_completion(MFC_CLI_SCHEMA)),
        (docs_dir / "cli-reference.md", generate_cli_reference(MFC_CLI_SCHEMA)),
        (Path(MFC_ROOT_DIR) / "toolchain" / "mfc-case-schema.json",
         json.dumps(generate_json_schema(include_descriptions=True), indent=2)),
        (docs_dir / "parameters.md", generate_parameter_docs()),
    ] + _constraint_docs(docs_dir)

    all_ok = True
    for path, content in files:
        if not _check_or_write(path, content, check_mode):
            all_ok = False

    if not all_ok:
        exit(1)

    if not check_mode:
        cons.print()
        cons.print("[bold]Files regenerated from cli/commands.py, params/definitions.py, and case_validator.py[/bold]")


def _generate_json_schema():
    """Generate JSON Schema and parameter documentation (standalone mode)."""
    from .params.generators.json_schema_gen import generate_json_schema, get_schema_stats
    from .params.generators.docs_gen import generate_parameter_docs
    from .ide import update_vscode_settings

    # Generate JSON Schema
    schema = generate_json_schema(include_descriptions=True)
    schema_path = Path(MFC_ROOT_DIR) / "toolchain" / "mfc-case-schema.json"
    with open(schema_path, 'w') as f:
        json.dump(schema, f, indent=2)

    # Generate parameter documentation
    docs_path = Path(MFC_ROOT_DIR) / "docs" / "documentation" / "parameters.md"
    docs_path.write_text(generate_parameter_docs())

    # Update VS Code settings
    update_vscode_settings()

    stats = get_schema_stats()

    cons.print(f"[green]Generated[/green] {schema_path}")
    cons.print(f"[green]Generated[/green] {docs_path}")
    cons.print()
    cons.print(f"[bold]Parameter Statistics:[/bold]")
    cons.print(f"  Total parameters: {stats['total_params']}")
    cons.print(f"  With constraints: {stats['with_constraints']}")
    cons.print(f"  With descriptions: {stats['with_descriptions']}")
    cons.print()
    cons.print("[bold]Parameter Lookup:[/bold]")
    cons.print("  CLI: [cyan]./mfc.sh params <query>[/cyan]")
    cons.print("  Docs: [cyan]docs/documentation/parameters.md[/cyan]")
