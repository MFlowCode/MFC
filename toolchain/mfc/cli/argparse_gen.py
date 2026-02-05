"""
Generate argparse parsers from CLI schema.

This module converts the declarative CLI schema into argparse ArgumentParsers.
"""

import argparse
import dataclasses
from typing import Dict, Tuple

from .schema import (
    CLISchema, Command, Argument, Positional,
    ArgAction, CommonArgumentSet
)


def _action_to_argparse(action: ArgAction) -> str:
    """Convert schema action to argparse action string."""
    return action.value


def _add_argument(parser: argparse.ArgumentParser, arg: Argument):
    """Add a single argument to parser."""
    kwargs: Dict = {
        "help": arg.help,
        "default": arg.default,
    }

    if arg.action != ArgAction.STORE:
        kwargs["action"] = _action_to_argparse(arg.action)
    else:
        if arg.type is not None:
            kwargs["type"] = arg.type
        if arg.choices is not None:
            kwargs["choices"] = arg.choices
        if arg.nargs is not None:
            # Handle "..." as REMAINDER
            if arg.nargs == "...":
                kwargs["nargs"] = argparse.REMAINDER
            else:
                kwargs["nargs"] = arg.nargs
        if arg.metavar is not None:
            kwargs["metavar"] = arg.metavar
        if arg.required:
            kwargs["required"] = arg.required

    if arg.dest is not None:
        kwargs["dest"] = arg.dest

    if arg.const is not None:
        kwargs["const"] = arg.const

    flags = arg.get_flags()
    parser.add_argument(*flags, **kwargs)


def _add_positional(parser: argparse.ArgumentParser, pos: Positional):
    """Add a positional argument to parser."""
    kwargs: Dict = {
        "help": pos.help,
        "metavar": pos.name.upper(),
    }

    if pos.type != str:
        kwargs["type"] = pos.type
    if pos.nargs is not None:
        kwargs["nargs"] = pos.nargs
    if pos.default is not None:
        kwargs["default"] = pos.default
    if pos.choices is not None:
        kwargs["choices"] = pos.choices

    parser.add_argument(pos.name, **kwargs)


def _add_mfc_config_arguments(parser: argparse.ArgumentParser, config):
    """
    Add MFCConfig boolean pair arguments dynamically.

    This handles --mpi/--no-mpi, --gpu/--no-gpu, etc. from MFCConfig dataclass.
    """
    # Import here to avoid circular dependency
    from ..state import gpuConfigOptions  # pylint: disable=import-outside-toplevel

    for f in dataclasses.fields(config):
        if f.name == 'gpu':
            parser.add_argument(
                f"--{f.name}",
                action="store",
                nargs='?',
                const=gpuConfigOptions.ACC.value,
                default=gpuConfigOptions.NONE.value,
                dest=f.name,
                choices=[e.value for e in gpuConfigOptions],
                help=f"Turn the {f.name} option to OpenACC or OpenMP."
            )
            parser.add_argument(
                f"--no-{f.name}",
                action="store_const",
                const=gpuConfigOptions.NONE.value,
                dest=f.name,
                help=f"Turn the {f.name} option OFF."
            )
        else:
            parser.add_argument(
                f"--{f.name}",
                action="store_true",
                help=f"Turn the {f.name} option ON."
            )
            parser.add_argument(
                f"--no-{f.name}",
                action="store_false",
                dest=f.name,
                help=f"Turn the {f.name} option OFF."
            )

    parser.set_defaults(**{
        f.name: getattr(config, f.name)
        for f in dataclasses.fields(config)
    })


def _add_common_arguments(
    parser: argparse.ArgumentParser,
    command: Command,
    common_sets: Dict[str, CommonArgumentSet],
    config=None
):
    """Add common arguments to a command parser."""
    for set_name in command.include_common:
        common_set = common_sets.get(set_name)
        if common_set is None:
            continue

        # Handle MFC config flags specially
        if common_set.mfc_config_flags and config is not None:
            _add_mfc_config_arguments(parser, config)
        else:
            for arg in common_set.arguments:
                _add_argument(parser, arg)


def _add_command_subparser(
    subparsers,
    cmd: Command,
    common_sets: Dict[str, CommonArgumentSet],
    config
) -> argparse.ArgumentParser:
    """Add a single command's subparser and return it."""
    subparser = subparsers.add_parser(
        name=cmd.name,
        aliases=cmd.aliases or [],
        help=cmd.help,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Add positional arguments first
    for pos in cmd.positionals:
        _add_positional(subparser, pos)

    # Add common arguments
    _add_common_arguments(subparser, cmd, common_sets, config)

    # Add command-specific arguments
    for arg in cmd.arguments:
        _add_argument(subparser, arg)

    # Add mutually exclusive groups
    for meg in cmd.mutually_exclusive:
        group = subparser.add_mutually_exclusive_group(required=meg.required)
        for arg in meg.arguments:
            _add_argument(group, arg)

    # Handle subcommands (e.g., packer pack, packer compare)
    if cmd.subcommands:
        sub_subparsers = subparser.add_subparsers(dest=cmd.name)
        for subcmd in cmd.subcommands:
            sub_sub = sub_subparsers.add_parser(
                name=subcmd.name,
                help=subcmd.help,
                formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            )
            for pos in subcmd.positionals:
                _add_positional(sub_sub, pos)
            for arg in subcmd.arguments:
                _add_argument(sub_sub, arg)

    return subparser


def generate_parser(
    schema: CLISchema,
    config=None  # MFCConfig instance
) -> Tuple[argparse.ArgumentParser, Dict[str, argparse.ArgumentParser]]:
    """
    Generate complete argparse parser from schema.

    Args:
        schema: The CLI schema definition
        config: Optional MFCConfig for dynamic boolean pairs

    Returns:
        Tuple of (main parser, dict mapping command names to subparsers)
    """
    parser = argparse.ArgumentParser(
        prog=schema.prog,
        description=schema.description,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False,  # We handle --help ourselves
    )

    # Add top-level arguments
    for arg in schema.arguments:
        _add_argument(parser, arg)

    # Build common sets lookup
    common_sets = {cs.name: cs for cs in schema.common_sets}

    # Add subparsers for commands
    subparsers = parser.add_subparsers(dest="command")
    subparser_map: Dict[str, argparse.ArgumentParser] = {}

    for cmd in schema.commands:
        subparser_map[cmd.name] = _add_command_subparser(subparsers, cmd, common_sets, config)

    return parser, subparser_map
