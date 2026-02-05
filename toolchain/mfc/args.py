"""
MFC CLI argument parsing.

This module provides argument parsing using auto-generated parsers
from the central CLI schema in cli/commands.py.
"""

import re
import sys
import os.path

from .common import MFCException
from .state import MFCConfig
from .cli.commands import MFC_CLI_SCHEMA, COMMAND_ALIASES
from .cli.argparse_gen import generate_parser
from .user_guide import (
    print_help, is_first_time_user, print_welcome,
    print_command_help, print_topic_help, print_help_topics,
)


def _get_command_from_args(args_list):
    """Extract command name from args list, resolving aliases.

    Scans for the first non-option token to support any top-level options
    that may appear before the command name.
    """
    # Skip the program name and any leading options (starting with '-')
    for token in args_list[1:]:
        if not token.startswith('-'):
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


# pylint: disable=too-many-locals, too-many-branches, too-many-statements
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
        extra_index = sys.argv.index('--')
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
            sys.stderr.write(f'\n{subparser.prog}: error: {message}\n')
            sys.exit(2)

        subparser.error = custom_error

    args: dict = vars(parser.parse_args(sys.argv[1:extra_index]))
    args["--"] = sys.argv[extra_index + 1:]

    # Handle --help at top level
    if args.get("help") and args["command"] is None:
        print_help()
        sys.exit(0)

    # Handle 'help' command
    if args["command"] == "help":
        topic = args.get("topic")
        if topic:
            print_topic_help(topic)
        else:
            print_help_topics()
        sys.exit(0)

    # Resolve command aliases
    if args["command"] in COMMAND_ALIASES:
        args["command"] = COMMAND_ALIASES[args["command"]]

    # Add default arguments of other subparsers
    # This ensures all argument keys exist even for commands that don't define them
    # Only process subparsers that have common arguments we need
    relevant_subparsers = ["run", "test", "build", "clean", "count", "count_diff", "validate"]
    for name in relevant_subparsers:
        if args["command"] == name:
            continue
        if name not in subparser_map:
            continue

        subparser = subparser_map[name]
        # Parse with dummy input to get defaults (suppress errors for required positionals)
        try:
            # Commands with required positional input need a dummy value
            if name in ["run", "validate"]:
                vals, _ = subparser.parse_known_args(["dummy_input.py"])
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
        args["name"] = re.sub(r'[\W_]+', '-', args["name"])

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
        from .test.cases import list_cases  # pylint: disable=import-outside-toplevel
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
            args[e] = os.path.abspath(args[e])

    return args
