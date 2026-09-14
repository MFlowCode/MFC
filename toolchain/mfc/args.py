"""
MFC CLI argument parsing.

This module provides argument parsing using auto-generated parsers
from the central CLI schema in cli/commands.py.
"""

import os.path
import re
import sys

from .cli.argparse_gen import generate_parser
from .cli.commands import COMMAND_ALIASES, MFC_CLI_SCHEMA
from .common import MFCException
from .printer import cons
from .state import MFCConfig
from .user_guide import (
    is_first_time_user,
    print_command_help,
    print_help,
    print_welcome,
)


def _get_command_from_args(args_list):
    """Extract command name from args list, resolving aliases.

    Scans for the first non-option token to support any top-level options
    that may appear before the command name.
    """
    # Skip the program name and any leading options (starting with '-')
    for token in args_list[1:]:
        if not token.startswith("-"):
            return COMMAND_ALIASES.get(token, token)
    return None


def _handle_enhanced_help(args_list):
    """Handle --help with enhanced output for known commands."""
    if len(args_list) >= 2 and args_list[1] in ("-h", "--help"):
        # ./mfc.sh --help -> show enhanced help
        print_help()
        sys.exit(0)

    if len(args_list) >= 3 and args_list[2] in ("-h", "--help"):
        # ./mfc.sh <command> --help -> show enhanced command help
        command = args_list[1]
        # Resolve alias
        command = COMMAND_ALIASES.get(command, command)
        # Print enhanced help, then let argparse show its help too
        print_command_help(command, show_argparse=True)
        # Return command so argparse can show its help
        return command
    return None


def _multi_value_flags(command_name: str):
    """Yield (flag_tokens, dest) for every list-valued option of a command.

    A list-valued option is one declared with nargs="+" or "*": argparse stores
    these with the default "store" action, so a second occurrence of the flag
    REPLACES whatever the first one collected instead of appending to it.
    Both the command's own arguments and the common sets it includes count.
    """
    command = MFC_CLI_SCHEMA.get_command(command_name)
    if command is None:
        return

    arguments = list(command.arguments)
    for set_name in command.include_common:
        common_set = MFC_CLI_SCHEMA.get_common_set(set_name)
        if common_set is not None:
            arguments.extend(common_set.arguments)

    for argument in arguments:
        if argument.nargs in ("+", "*"):
            yield argument.get_flags(), argument.get_dest()


def _warn_on_repeated_multi_value_flags(command_name: str, cli_argv):
    """Warn when a list-valued flag was passed more than once.

    `-t` takes a space-separated list, so `-t pre_process -t simulation` does not
    append -- argparse keeps only the last occurrence and the earlier targets are
    dropped without a word. That is quiet and lands far from its cause: the
    generated batch script simply has one fewer step than the user expected.
    Repeating a flag appends in many other CLIs, and `-t simulation` on its own is
    a perfectly legitimate invocation (restarting from existing data), so nothing
    downstream can tell the two apart. Warn here, where we can still see that the
    flag was written twice.
    """
    for flags, dest in _multi_value_flags(command_name):
        occurrences = sum(1 for tok in cli_argv if tok in flags or any(tok.startswith(f"{flag}=") for flag in flags))
        if occurrences > 1:
            joined = " / ".join(flags)
            cons.print(
                f"[yellow]{joined} was given {occurrences} times, but it takes a space-separated list "
                f"and only the last occurrence is kept. Earlier values were discarded; "
                f"pass them together as e.g. --{dest.replace('_', '-')} A B.[/yellow]"
            )


def parse(config: MFCConfig):
    """Parse command line arguments using the CLI schema."""
    # Handle enhanced help before argparse
    help_command = _handle_enhanced_help(sys.argv)

    # Generate parser from schema
    parser, subparser_map = generate_parser(MFC_CLI_SCHEMA, config)

    # If enhanced help was printed, also show argparse help and exit
    if help_command and help_command in subparser_map:
        subparser_map[help_command].print_help()
        sys.exit(0)

    try:
        extra_index = sys.argv.index("--")
    except ValueError:
        extra_index = len(sys.argv)

    # Patch subparser error methods to show full help on error
    attempted_command = _get_command_from_args(sys.argv)
    if attempted_command and attempted_command in subparser_map:
        subparser = subparser_map[attempted_command]

        def custom_error(message):
            # Show enhanced help + full argparse help (like -h would)
            print_command_help(attempted_command, show_argparse=False)
            subparser.print_help()
            sys.stdout.flush()  # Ensure help prints before error
            sys.stderr.write(f"\n{subparser.prog}: error: {message}\n")
            sys.exit(2)

        subparser.error = custom_error

    cli_argv = sys.argv[1:extra_index]
    args: dict = vars(parser.parse_args(cli_argv))
    args["--"] = sys.argv[extra_index + 1 :]
    args["targets_explicit"] = any(tok in ("-t", "--targets") for tok in cli_argv)

    if attempted_command:
        _warn_on_repeated_multi_value_flags(attempted_command, cli_argv)

    # Handle --help at top level
    if args.get("help") and args["command"] is None:
        print_help()
        sys.exit(0)

    # Resolve command aliases
    if args["command"] in COMMAND_ALIASES:
        args["command"] = COMMAND_ALIASES[args["command"]]

    # Add default arguments of other subparsers
    # This ensures all argument keys exist even for commands that don't define them
    # Only process subparsers that have common arguments we need
    relevant_subparsers = ["run", "test", "build", "clean", "count", "validate", "viz"]
    for name in relevant_subparsers:
        if args["command"] == name:
            continue
        if name not in subparser_map:
            continue

        subparser = subparser_map[name]
        # Parse with dummy input to get defaults (suppress errors for required positionals)
        try:
            # Commands with required positional input need a dummy value
            if name in ["run", "validate", "viz"]:
                vals, _ = subparser.parse_known_args(["dummy_dir/"])
            elif name == "build":
                vals, _ = subparser.parse_known_args([])
            else:
                vals, _ = subparser.parse_known_args([])
        except SystemExit:
            continue  # Skip if parsing fails

        for key, val in vars(vals).items():
            if key == "input":
                args[key] = args.get(key)
            elif key not in args:
                args[key] = args.get(key, val)

    if args["command"] is None:
        # Show welcome for first-time users, otherwise show enhanced help
        if is_first_time_user():
            print_welcome()
        else:
            print_help()
        sys.exit(0)

    # "Slugify" the name of the job (only for batch jobs, not for new command)
    if args.get("name") is not None and isinstance(args["name"], str) and args["command"] != "new":
        args["name"] = re.sub(r"[\W_]+", "-", args["name"])

    # We need to check for some invalid combinations of arguments because of
    # the limitations of argparse.
    if args["command"] == "build":
        if (args["input"] is not None) ^ args["case_optimization"]:
            raise MFCException("./mfc.sh build's --case-optimization and --input must be used together.")
    if args["command"] == "run":
        if args["binary"] is not None and args["engine"] != "interactive":
            raise MFCException("./mfc.sh run's --binary can only be used with --engine=interactive.")

    # Resolve test case defaults (deferred to avoid slow startup for non-test commands)
    if args["command"] == "test":
        from .test.cases import list_cases

        test_cases = list_cases()
        if args.get("from") is None:
            args["from"] = test_cases[0].get_uuid()
        if args.get("to") is None:
            args["to"] = test_cases[-1].get_uuid()

    # Input files to absolute paths
    for e in ["input", "input1", "input2"]:
        if e not in args:
            continue

        if args.get(e) is not None:
            if isinstance(args[e], list):  # validate accepts multiple case files (nargs="+")
                args[e] = [os.path.abspath(p) for p in args[e]]
            else:
                args[e] = os.path.abspath(args[e])

    return args
