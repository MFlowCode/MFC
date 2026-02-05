"""
Shell completion installation for MFC.

Installs completion scripts to ~/.local/share/mfc/completions/ and
configures the user's shell to source them automatically.
"""

import os
import shutil
from pathlib import Path

from .printer import cons
from .common import MFC_ROOT_DIR


# Installation directory (user-local, independent of MFC clone location)
COMPLETION_INSTALL_DIR = Path.home() / ".local" / "share" / "mfc" / "completions"

# Source files
BASH_COMPLETION_SRC = Path(MFC_ROOT_DIR) / "toolchain" / "completions" / "mfc.bash"
ZSH_COMPLETION_SRC = Path(MFC_ROOT_DIR) / "toolchain" / "completions" / "_mfc"

# Shell RC files
BASHRC = Path.home() / ".bashrc"
ZSHRC = Path.home() / ".zshrc"

# Lines to add to RC files
BASH_SOURCE_LINE = f'[ -f "{COMPLETION_INSTALL_DIR}/mfc.bash" ] && source "{COMPLETION_INSTALL_DIR}/mfc.bash"'
ZSH_FPATH_LINE = f'fpath=("{COMPLETION_INSTALL_DIR}" $fpath)'


def _line_in_file(filepath: Path, line: str) -> bool:
    """Check if a line (or similar) already exists in a file."""
    if not filepath.exists():
        return False

    content = filepath.read_text()
    # Check for the exact line or the key part of it
    return line in content or str(COMPLETION_INSTALL_DIR) in content


def _append_to_file(filepath: Path, line: str, comment: str = None):
    """Append a line to a file if it doesn't already exist."""
    if _line_in_file(filepath, line):
        cons.print(f"  [dim]Already configured in {filepath}[/dim]")
        return False

    with open(filepath, "a") as f:
        f.write(f"\n# {comment}\n" if comment else "\n")
        f.write(f"{line}\n")

    return True


def install_bash():
    """Install bash completion."""
    cons.print("[bold]Installing Bash completion...[/bold]")

    # Create installation directory
    COMPLETION_INSTALL_DIR.mkdir(parents=True, exist_ok=True)

    # Copy completion script
    dest = COMPLETION_INSTALL_DIR / "mfc.bash"
    if BASH_COMPLETION_SRC.exists():
        shutil.copy2(BASH_COMPLETION_SRC, dest)
        cons.print(f"  [green]✓[/green] Copied completion script to [cyan]{dest}[/cyan]")
    else:
        cons.print(f"  [red]✗[/red] Source file not found: {BASH_COMPLETION_SRC}")
        return False

    # Add source line to .bashrc (creates if needed)
    if _append_to_file(BASHRC, BASH_SOURCE_LINE, "MFC shell completion"):
        cons.print(f"  [green]✓[/green] Added source line to [cyan]{BASHRC}[/cyan]")

    cons.print()
    cons.print("[green]Bash completion installed![/green]")
    cons.print()
    cons.print("To activate now (without restarting your shell):")
    cons.print(f"  [cyan]source {dest}[/cyan]")
    cons.print()
    cons.print("Or start a new terminal session.")

    return True


def install_zsh():
    """Install zsh completion."""
    cons.print("[bold]Installing Zsh completion...[/bold]")

    # Create installation directory
    COMPLETION_INSTALL_DIR.mkdir(parents=True, exist_ok=True)

    # Copy completion script
    dest = COMPLETION_INSTALL_DIR / "_mfc"
    if ZSH_COMPLETION_SRC.exists():
        shutil.copy2(ZSH_COMPLETION_SRC, dest)
        cons.print(f"  [green]✓[/green] Copied completion script to [cyan]{dest}[/cyan]")
    else:
        cons.print(f"  [red]✗[/red] Source file not found: {ZSH_COMPLETION_SRC}")
        return False

    # Add fpath line to .zshrc
    if ZSHRC.exists():
        if _append_to_file(ZSHRC, ZSH_FPATH_LINE, "MFC shell completion"):
            cons.print(f"  [green]✓[/green] Added fpath to [cyan]{ZSHRC}[/cyan]")

        # Also need to ensure compinit is called
        compinit_line = "autoload -Uz compinit && compinit"
        if not _line_in_file(ZSHRC, "compinit"):
            _append_to_file(ZSHRC, compinit_line)
            cons.print(f"  [green]✓[/green] Added compinit to [cyan]{ZSHRC}[/cyan]")
    else:
        cons.print(f"  [yellow]![/yellow] {ZSHRC} not found - you may need to configure manually")
        cons.print(f"    Add to your shell config: {ZSH_FPATH_LINE}")

    cons.print()
    cons.print("[green]Zsh completion installed![/green]")
    cons.print()
    cons.print("Start a new terminal session to activate.")

    return True


def install_auto():
    """Auto-detect shell and install appropriate completion."""
    shell = os.environ.get("SHELL", "")

    if "zsh" in shell:
        return install_zsh()

    # Default to bash
    return install_bash()


def uninstall():
    """Remove installed completion files and RC modifications."""
    cons.print("[bold]Uninstalling MFC shell completion...[/bold]")

    # Remove completion directory
    if COMPLETION_INSTALL_DIR.exists():
        shutil.rmtree(COMPLETION_INSTALL_DIR)
        cons.print(f"  [green]✓[/green] Removed [cyan]{COMPLETION_INSTALL_DIR}[/cyan]")

    # Note: We don't automatically remove lines from .bashrc/.zshrc
    # as that's more risky. Just inform the user.
    cons.print()
    cons.print("[yellow]Note:[/yellow] You may want to manually remove the MFC completion lines from:")
    cons.print(f"  [cyan]{BASHRC}[/cyan]")
    cons.print(f"  [cyan]{ZSHRC}[/cyan]")

    return True


def show_status():
    """Show current completion installation status."""
    cons.print("[bold]MFC Shell Completion Status[/bold]")
    cons.print()

    # Check installation directory
    if COMPLETION_INSTALL_DIR.exists():
        cons.print(f"  [green]✓[/green] Install directory: [cyan]{COMPLETION_INSTALL_DIR}[/cyan]")

        bash_installed = (COMPLETION_INSTALL_DIR / "mfc.bash").exists()
        zsh_installed = (COMPLETION_INSTALL_DIR / "_mfc").exists()

        if bash_installed:
            cons.print("  [green]✓[/green] Bash completion installed")
        else:
            cons.print("  [dim]✗ Bash completion not installed[/dim]")

        if zsh_installed:
            cons.print("  [green]✓[/green] Zsh completion installed")
        else:
            cons.print("  [dim]✗ Zsh completion not installed[/dim]")
    else:
        cons.print(f"  [dim]✗ Not installed[/dim]")

    cons.print()

    # Check RC files
    if BASHRC.exists() and _line_in_file(BASHRC, str(COMPLETION_INSTALL_DIR)):
        cons.print(f"  [green]✓[/green] Configured in [cyan]{BASHRC}[/cyan]")
    else:
        cons.print(f"  [dim]✗ Not configured in {BASHRC}[/dim]")

    if ZSHRC.exists() and _line_in_file(ZSHRC, str(COMPLETION_INSTALL_DIR)):
        cons.print(f"  [green]✓[/green] Configured in [cyan]{ZSHRC}[/cyan]")
    else:
        cons.print(f"  [dim]✗ Not configured in {ZSHRC}[/dim]")


def completion():
    """Main entry point for completion command."""
    # pylint: disable=import-outside-toplevel
    from .state import ARG

    action = ARG("completion_action")

    if action == "install":
        shell = ARG("completion_shell")
        if shell == "bash":
            install_bash()
        elif shell == "zsh":
            install_zsh()
        else:
            install_auto()
    elif action == "uninstall":
        uninstall()
    elif action == "status":
        show_status()
    else:
        # Default: show status and usage
        show_status()
        cons.print()
        cons.print("[bold]Usage:[/bold]")
        cons.print("  [cyan]./mfc.sh completion install[/cyan]        Auto-detect shell and install")
        cons.print("  [cyan]./mfc.sh completion install bash[/cyan]   Install bash completion")
        cons.print("  [cyan]./mfc.sh completion install zsh[/cyan]    Install zsh completion")
        cons.print("  [cyan]./mfc.sh completion uninstall[/cyan]      Remove completion files")
        cons.print("  [cyan]./mfc.sh completion status[/cyan]         Show installation status")
