"""
MFC Validate Command - Validate a case file without building or running.
"""

import os

from .state import ARG
from .printer import cons
from .run import input as run_input
from .case_validator import CaseValidator, CaseConstraintError
from .common import MFCException


def validate():
    """Validate a case file without building or running."""
    input_file = ARG("input")

    if not os.path.isfile(input_file):
        cons.print(f"[bold red]Error:[/bold red] File not found: {input_file}")
        exit(1)

    cons.print(f"Validating [bold magenta]{input_file}[/bold magenta]...\n")

    try:
        # Step 1: Load and parse case file (checks syntax)
        case = run_input.load(input_file, do_print=False)
        cons.print("[bold green]✓[/bold green] Syntax valid - case file parsed successfully")
        cons.print(f"  [dim]Loaded {len(case.params)} parameters[/dim]")

        # Step 2: Run constraint validation for each stage
        stages = ['pre_process', 'simulation', 'post_process']
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
                for line in str(e).split('\n'):
                    if line.strip():
                        cons.print(f"    [dim]{line}[/dim]")

        # Step 3: Show summary
        cons.print()
        if all_passed:
            cons.print("[bold green]Case validation complete - all checks passed![/bold green]")
        else:
            cons.print("[bold yellow]Case validation complete with warnings.[/bold yellow]")
            cons.print("[dim]Note: Some constraint violations may be OK if you're not using that stage.[/dim]")

    except MFCException as e:
        cons.print(f"\n[bold red]✗ Validation failed:[/bold red]")
        cons.print(f"{e}")
        exit(1)
