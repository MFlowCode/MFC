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

    if args["command"] in ("build", "run", "test"):
        acc = args.get("gpu") == "acc"
        case_opt_acc = args.get("case_optimization") and acc
        if args.get("fused_rhs"):
            if not args.get("case_optimization"):
                raise MFCException("--fused-rhs requires --case-optimization.")
            if not acc:
                raise MFCException("--fused-rhs requires an OpenACC GPU build.")
        if args.get("active_box"):
            if not args.get("fused_rhs"):
                raise MFCException("--active-box requires --fused-rhs (it windows the fused kernel).")
            if not args.get("case_optimization"):
                raise MFCException("--active-box requires --case-optimization.")
            if not acc:
                raise MFCException("--active-box requires an OpenACC GPU build.")
        if args.get("nvshmem"):
            if not args.get("mpi"):
                raise MFCException("--nvshmem requires --mpi.")
            if not acc:
                raise MFCException("--nvshmem requires an OpenACC GPU build.")
            if args.get("unified"):
                raise MFCException("--nvshmem requires --no-unified.")
        if args.get("async_queue"):
            if not acc:
                raise MFCException("--async-queue requires an OpenACC GPU build.")
            if args.get("unified"):
                raise MFCException("--async-queue requires --no-unified.")
        if args.get("halo_single_round") and not args.get("mpi"):
            raise MFCException("--halo-single-round requires --mpi.")
        if args.get("cuda_graphs"):
            if not acc:
                raise MFCException("--cuda-graphs requires an OpenACC GPU build.")
            if args.get("unified"):
                raise MFCException("--cuda-graphs requires --no-unified.")
            if not args.get("async_queue"):
                raise MFCException("--cuda-graphs requires --async-queue.")
        if args.get("cbc_no_host_sync") and not acc:
            raise MFCException("--cbc-no-host-sync requires an OpenACC GPU build.")
        if args.get("cbc_parallel_scatter") and not args.get("cbc_no_host_sync"):
            raise MFCException("--cbc-parallel-scatter requires --cbc-no-host-sync.")
        if args.get("halo_index_cache") and not args.get("halo_single_round"):
            raise MFCException("--halo-index-cache requires --halo-single-round.")
        if args.get("hllc_register_tuning") and not case_opt_acc:
            raise MFCException("--hllc-register-tuning requires a case-optimized OpenACC GPU build.")
        if args.get("weno_direct_loads"):
            if not case_opt_acc:
                raise MFCException("--weno-direct-loads requires a case-optimized OpenACC GPU build.")
            if not (args.get("fused_rhs") or args.get("hllc_register_tuning")):
                raise MFCException("--weno-direct-loads requires --fused-rhs or --hllc-register-tuning (the restructured WENO module).")
        if args.get("direct_riemann_flux") and not case_opt_acc:
            raise MFCException("--direct-riemann-flux requires a case-optimized OpenACC GPU build.")
        if args.get("persistent_bc") and not acc:
            raise MFCException("--persistent-bc requires an OpenACC GPU build.")
        if args.get("fused_bubble_source") and not acc:
            raise MFCException("--fused-bubble-source requires an OpenACC GPU build.")
        if args.get("qbmm_device_inversion"):
            if not acc:
                raise MFCException("--qbmm-device-inversion requires an OpenACC GPU build.")
            if not args.get("case_optimization"):
                raise MFCException("--qbmm-device-inversion requires --case-optimization.")
        if args.get("accurate_grind_time") and not acc:
            raise MFCException("--accurate-grind-time requires an OpenACC GPU build (NVHPC or CCE).")
        if args.get("async_queue") and not args.get("accurate_grind_time"):
            raise MFCException("--async-queue requires --accurate-grind-time. The host runs ahead of the device, so the printed grind time needs device events.")

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
