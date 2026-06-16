"""
MFC Validate Command - Validate a case file without building or running.
"""

import os
import sys

from .case_validator import CaseConstraintError, CaseValidator
from .common import MFCException
from .printer import cons
from .run import input as run_input
from .state import ARG


def validate():
    """Validate a case file without building or running."""
    input_file = ARG("input")

    if not os.path.isfile(input_file):
        cons.print(f"[bold red]Error:[/bold red] File not found: {input_file}")
        sys.exit(1)

    cons.print(f"Validating [bold magenta]{input_file}[/bold magenta]...\n")

    try:
        # Step 1: Load and parse case file (checks syntax)
        case = run_input.load(input_file, do_print=False)
        cons.print("[bold green]✓[/bold green] Syntax valid - case file parsed successfully")
        cons.print(f"  [dim]Loaded {len(case.params)} parameters[/dim]")

        # Step 2: Run constraint validation for each stage
        stages = ["pre_process", "simulation", "post_process"]
        all_passed = True

        for stage in stages:
            try:
                validator = CaseValidator(case.params)
                warnings = validator.validate(stage)
                if warnings:
                    cons.print(f"[bold green]✓[/bold green] {stage} constraints passed (with warnings)")
                    for warning in warnings:
                        cons.print(f"    [yellow]⚠ {warning}[/yellow]")
                else:
                    cons.print(f"[bold green]✓[/bold green] {stage} constraints passed")
            except CaseConstraintError as e:
                all_passed = False
                cons.print(f"[bold yellow]![/bold yellow] {stage} constraints: issues found")
                # Show the constraint violations indented
                for line in str(e).split("\n"):
                    if line.strip():
                        cons.print(f"    [dim]{line}[/dim]")

        # Step 3: Show summary
        cons.print()
        if all_passed:
            cons.print("[bold green]Case validation complete - all checks passed![/bold green]")
        else:
            cons.print("[bold yellow]Case validation complete with warnings.[/bold yellow]")
            cons.print("[dim]Note: Some constraint violations may be OK if you're not using that stage.[/dim]")

        if ARG("migrate") and not all_passed:
            cons.print("[yellow]Skipping --migrate: fix the validation issues above first.[/yellow]")
        elif ARG("migrate"):
            from .params.migrate import migrate_text

            with open(input_file, "r") as f:
                original = f.read()
            new_text, n = migrate_text(original)
            if n == 0:
                cons.print("No integer codes to migrate.")
                return
            with open(input_file, "w") as f:
                f.write(new_text)
            # Safety: the migrated file must load to the same normalized parameters.
            try:
                re_case = run_input.load(input_file, do_print=False)
            except MFCException:
                with open(input_file, "w") as f:
                    f.write(original)
                cons.print("[bold red]Error:[/bold red] Migrated file failed to load; original restored.")
                sys.exit(1)
            if re_case.params != case.params:
                with open(input_file, "w") as f:
                    f.write(original)
                cons.print("[bold red]Error:[/bold red] Migration changed case semantics; file restored.")
                sys.exit(1)
            cons.print(f"[bold green]✓[/bold green] Migrated {n} value(s) to named syntax")

    except MFCException as e:
        cons.print("\n[bold red]✗ Validation failed:[/bold red]")
        cons.print(f"{e}")
        sys.exit(1)
