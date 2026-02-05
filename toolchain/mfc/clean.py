"""
MFC Clean Command - Remove build artifacts.
"""

import os
import shutil

from .printer import cons
from .common import MFC_BUILD_DIR


def clean():
    """Remove the build directory and all build artifacts."""
    if os.path.isdir(MFC_BUILD_DIR):
        cons.print(f"Removing [bold magenta]{MFC_BUILD_DIR}[/bold magenta]...")
        try:
            shutil.rmtree(MFC_BUILD_DIR)
            cons.print("[bold green]Build directory cleaned successfully.[/bold green]")
        except OSError as e:
            cons.print(f"[bold red]Error cleaning build directory:[/bold red] {e}")
    else:
        cons.print("[yellow]Build directory does not exist, nothing to clean.[/yellow]")
