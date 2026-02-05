"""
CLI Schema and Auto-Generation Module.

This package provides the single source of truth for MFC's CLI definition,
with generators for argparse, shell completions, and documentation.

Usage:
    from mfc.cli.commands import MFC_CLI_SCHEMA
    from mfc.cli.argparse_gen import generate_parser
    from mfc.cli.completion_gen import generate_bash_completion, generate_zsh_completion
"""

from .schema import (
    CLISchema,
    Command,
    Argument,
    Positional,
    Example,
    CommonArgumentSet,
    MutuallyExclusiveGroup,
    ArgAction,
    CompletionType,
    Completion,
)

__all__ = [
    "CLISchema",
    "Command",
    "Argument",
    "Positional",
    "Example",
    "CommonArgumentSet",
    "MutuallyExclusiveGroup",
    "ArgAction",
    "CompletionType",
    "Completion",
]
