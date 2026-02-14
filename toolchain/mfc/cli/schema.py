"""
CLI Schema Dataclass Definitions.

This module defines the dataclasses used to describe CLI commands, arguments,
and their properties. These definitions serve as the single source of truth
for generating argparse parsers, shell completions, and documentation.
"""

from dataclasses import dataclass, field
from enum import Enum, auto
from typing import List, Optional, Any, Union


class ArgAction(Enum):
    """Supported argparse actions."""
    STORE = "store"
    STORE_TRUE = "store_true"
    STORE_FALSE = "store_false"
    STORE_CONST = "store_const"
    COUNT = "count"
    APPEND = "append"


class CompletionType(Enum):
    """Types of shell completion behavior."""
    NONE = auto()           # No completion
    FILES = auto()          # All file completion
    FILES_PY = auto()       # Python files only (*.py)
    FILES_PACK = auto()     # Pack files only (*.pack)
    FILES_YAML = auto()     # YAML files only (*.yaml, *.yml)
    DIRECTORIES = auto()    # Directory completion
    CHOICES = auto()        # Static choices from choices list


@dataclass
class Completion:
    """Completion configuration for an argument."""
    type: CompletionType = CompletionType.NONE
    choices: Optional[List[str]] = None


@dataclass
class Argument:  # pylint: disable=too-many-instance-attributes
    """
    Definition of a single CLI option argument (--flag).

    This represents one add_argument() call for a flag-style argument.
    """
    # Identity
    name: str                           # Long form without dashes (e.g., "targets")
    short: Optional[str] = None         # Short form without dash (e.g., "t")

    # Argparse configuration
    help: str = ""
    action: ArgAction = ArgAction.STORE
    type: Optional[type] = None         # str, int, float, etc.
    default: Any = None
    choices: Optional[List[Any]] = None
    nargs: Optional[Union[str, int]] = None  # "+", "*", "?", int, or "..." for REMAINDER
    metavar: Optional[str] = None
    required: bool = False
    dest: Optional[str] = None          # Override destination name
    const: Any = None                   # For store_const action

    # Completion
    completion: Completion = field(default_factory=Completion)

    def get_flags(self) -> List[str]:
        """Return the flag strings for argparse."""
        flags = []
        if self.short:
            flags.append(f"-{self.short}")
        flags.append(f"--{self.name}")
        return flags

    def get_dest(self) -> str:
        """Return the destination name for argparse."""
        if self.dest:
            return self.dest
        return self.name.replace("-", "_")


@dataclass
class Positional:
    """Definition of a positional argument."""
    name: str                           # Metavar and destination
    help: str = ""
    type: type = str
    nargs: Optional[Union[str, int]] = None
    default: Any = None
    choices: Optional[List[str]] = None

    # Completion
    completion: Completion = field(default_factory=Completion)


@dataclass
class Example:
    """A usage example for documentation."""
    command: str
    description: str


@dataclass
class MutuallyExclusiveGroup:
    """A group where only one argument can be specified."""
    arguments: List[Argument] = field(default_factory=list)
    required: bool = False


@dataclass
class Command:  # pylint: disable=too-many-instance-attributes
    """
    Definition of a CLI command/subcommand.

    This is the main building block for the CLI structure.
    """
    # Identity
    name: str
    help: str
    aliases: List[str] = field(default_factory=list)

    # Arguments
    positionals: List[Positional] = field(default_factory=list)
    arguments: List[Argument] = field(default_factory=list)
    mutually_exclusive: List[MutuallyExclusiveGroup] = field(default_factory=list)

    # Inherit common argument sets
    include_common: List[str] = field(default_factory=list)  # e.g., ["targets", "jobs"]

    # Subcommands (for nested commands like "packer pack")
    subcommands: List["Command"] = field(default_factory=list)

    # Documentation
    description: Optional[str] = None   # Long description for docs
    examples: List[Example] = field(default_factory=list)
    key_options: List[tuple] = field(default_factory=list)  # (option, description) pairs

    # Handler module path (for dispatch)
    handler: Optional[str] = None       # Module.function path


@dataclass
class CommonArgumentSet:
    """
    A reusable set of arguments that can be included in multiple commands.

    Replaces the add_common_arguments() function pattern.
    """
    name: str                           # Identifier for include_common
    arguments: List[Argument] = field(default_factory=list)
    # For MFCConfig flags that need --X and --no-X pairs
    mfc_config_flags: bool = False


@dataclass
class CLISchema:
    """
    The complete CLI schema - single source of truth.

    This dataclass contains all information needed to generate:
    - argparse parsers
    - Bash completions
    - Zsh completions
    - User guide help content
    - CLI reference documentation
    """
    prog: str = "./mfc.sh"
    description: str = ""

    # Top-level arguments (like --help)
    arguments: List[Argument] = field(default_factory=list)

    # Commands
    commands: List[Command] = field(default_factory=list)

    # Reusable argument sets
    common_sets: List[CommonArgumentSet] = field(default_factory=list)

    # Help topics (separate from commands)
    help_topics: dict = field(default_factory=dict)

    def get_command(self, name: str) -> Optional[Command]:
        """Get a command by name or alias."""
        for cmd in self.commands:
            if cmd.name == name or name in cmd.aliases:
                return cmd
        return None

    def get_all_command_names(self) -> List[str]:
        """Get all command names and aliases."""
        names = []
        for cmd in self.commands:
            names.append(cmd.name)
            names.extend(cmd.aliases)
        return names

    def get_common_set(self, name: str) -> Optional[CommonArgumentSet]:
        """Get a common argument set by name."""
        for cs in self.common_sets:
            if cs.name == name:
                return cs
        return None
