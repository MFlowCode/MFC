import re, sys, os.path, argparse, dataclasses

from .common       import MFCException, format_list_to_string
from .state        import gpuConfigOptions, MFCConfig
from .user_guide   import (print_help, is_first_time_user, print_welcome,
                           print_command_help, print_topic_help, print_help_topics,
                           COMMAND_ALIASES, HELP_TOPICS)

# Hardcoded values to avoid heavy imports at startup
# These must be kept in sync with build.py and run/run.py
_TARGET_NAMES = ['fftw', 'hdf5', 'silo', 'lapack', 'hipfort',
                 'pre_process', 'simulation', 'post_process',
                 'syscheck', 'documentation']
_DEFAULT_TARGET_NAMES = ['pre_process', 'simulation', 'post_process']
_TEMPLATE_NAMES = ['bridges2', 'carpenter', 'carpenter-cray', 'default',
                   'delta', 'deltaai', 'frontier', 'hipergator', 'nautilus',
                   'oscar', 'phoenix', 'phoenix-bench', 'santis', 'tuo']


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
    # Handle enhanced help before argparse
    help_command = _handle_enhanced_help(sys.argv)

    parser = argparse.ArgumentParser(
        prog="./mfc.sh",
        description="""\
Welcome to the MFC master script. This tool automates and manages building, testing, \
running, and cleaning of MFC in various configurations on all supported platforms. \
The README documents this tool and its various commands in more detail. To get \
started, run ./mfc.sh build -h.""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False,  # We handle --help ourselves
    )
    parser.add_argument("-h", "--help", action="store_true", help="Show help message")

    # Here are all of the parser arguments that call functions in other python files
    parsers = parser.add_subparsers(dest="command")
    run        = parsers.add_parser(name="run",        help="Run a case with MFC.",                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    test       = parsers.add_parser(name="test",       help="Run MFC's test suite.",                  formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    build      = parsers.add_parser(name="build",      help="Build MFC and its dependencies.",        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    clean      = parsers.add_parser(name="clean",      help="Clean build artifacts.",                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    bench      = parsers.add_parser(name="bench",      help="Benchmark MFC (for CI).",                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    bench_diff = parsers.add_parser(name="bench_diff", help="Compare MFC Benchmarks (for CI).",       formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    count      = parsers.add_parser(name="count",      help="Count LOC in MFC.",                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    count_diff = parsers.add_parser(name="count_diff", help="Count LOC in MFC.",                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    packer     = parsers.add_parser(name="packer",     help="Packer utility (pack/unpack/compare).",  formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    validate   = parsers.add_parser(name="validate",   help="Validate a case file without running.",  formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    init       = parsers.add_parser(name="init",       help="Create a new case from a template.",     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parsers.add_parser(name="interactive", help="Launch interactive menu-driven interface.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    comp       = parsers.add_parser(name="completion",  help="Install shell tab-completion.",            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Help command with topics
    help_parser = parsers.add_parser(name="help", help="Show help on a topic.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    help_parser.add_argument("topic", metavar="TOPIC", type=str, nargs="?", default=None,
                             choices=list(HELP_TOPICS.keys()),
                             help=f"Help topic: {', '.join(HELP_TOPICS.keys())}")

    # Command aliases (b=build, r=run, t=test, v=validate, c=clean)
    parsers.add_parser(name="b", help="Alias for 'build'.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parsers.add_parser(name="r", help="Alias for 'run'.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parsers.add_parser(name="t", help="Alias for 'test'.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parsers.add_parser(name="v", help="Alias for 'validate'.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parsers.add_parser(name="c", help="Alias for 'clean'.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Completion subcommands
    comp.add_argument("completion_action", metavar="ACTION", type=str, nargs="?", default=None,
                      choices=["install", "uninstall", "status"],
                      help="Action: install, uninstall, or status")
    comp.add_argument("completion_shell", metavar="SHELL", type=str, nargs="?", default=None,
                      choices=["bash", "zsh"],
                      help="Shell type: bash or zsh (auto-detected if not specified)")

    # These parser arguments all call BASH scripts, and they only exist so that they show up in the help message
    parsers.add_parser(name="load",       help="Loads the MFC environment with source.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parsers.add_parser(name="lint",       help="Lints all code after editing.",          formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parsers.add_parser(name="format",     help="Formats all code after editing.",        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parsers.add_parser(name="spelling",   help="Runs the spell checker after editing.",  formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    packers = packer.add_subparsers(dest="packer")
    pack = packers.add_parser(name="pack", help="Pack a case into a single file.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    pack.add_argument("input", metavar="INPUT", type=str, default="", help="Input file of case to pack.")
    pack.add_argument("-o", "--output", metavar="OUTPUT", type=str, default=None, help="Base name of output file.")

    compare = packers.add_parser(name="compare", help="Compare two cases.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    compare.add_argument("input1", metavar="INPUT1", type=str, default=None, help="First pack file.")
    compare.add_argument("input2", metavar="INPUT2", type=str, default=None, help="Second pack file.")
    compare.add_argument("-rel", "--reltol", metavar="RELTOL", type=float, default=1e-12, help="Relative tolerance.")
    compare.add_argument("-abs", "--abstol", metavar="ABSTOL", type=float, default=1e-12, help="Absolute tolerance.")

    def add_common_arguments(p: argparse.ArgumentParser, mask = None):
        if mask is None:
            mask = ""

        if "t" not in mask:
            p.add_argument("-t", "--targets", metavar="TARGET", nargs="+", type=str.lower, choices=_TARGET_NAMES,
                           default=_DEFAULT_TARGET_NAMES,
                           help=f"Space separated list of targets to act upon. Allowed values are: {format_list_to_string(_TARGET_NAMES)}.")

        if "m" not in mask:
            for f in dataclasses.fields(config):
                if f.name == 'gpu':
                    p.add_argument(f"--{f.name}", action="store", nargs='?', const= gpuConfigOptions.ACC.value,default=gpuConfigOptions.NONE.value, dest=f.name, choices=[e.value for e in gpuConfigOptions], help=f"Turn the {f.name} option to OpenACC or OpenMP.")
                    p.add_argument(f"--no-{f.name}", action="store_const", const = gpuConfigOptions.NONE.value, dest=f.name, help=f"Turn the {f.name} option OFF.")
                    continue
                p.add_argument(   f"--{f.name}", action="store_true",               help=f"Turn the {f.name} option ON.")
                p.add_argument(f"--no-{f.name}", action="store_false", dest=f.name, help=f"Turn the {f.name} option OFF.")

            p.set_defaults(**{ f.name: getattr(config, f.name) for f in dataclasses.fields(config) })

        if "j" not in mask:
            p.add_argument("-j", "--jobs", metavar="JOBS", type=int, default=1, help="Allows for JOBS concurrent jobs.")

        if "v" not in mask:
            p.add_argument("-v", "--verbose", action="store_true", help="Enables verbose compiler & linker output.")

        if "d" not in mask:
            p.add_argument("-d", "--debug-log", action="store_true", dest="debug_log", help="Enable debug logging for troubleshooting.")

        if "g" not in mask:
            p.add_argument("-g", "--gpus", nargs="+", type=int, default=None, help="(Optional GPU override) List of GPU #s to use (environment default if unspecified).")

    # BUILD
    add_common_arguments(build, "g")
    build.add_argument("-i", "--input", type=str, default=None, help="(GPU Optimization) Build a version of MFC optimized for a case.")
    build.add_argument("--case-optimization", action="store_true", default=False, help="(GPU Optimization) Compile MFC targets with some case parameters hard-coded (requires --input).")

    # TEST
    # Note: list_cases() is deferred until actually needed (when running 'test' command)
    # to avoid slow startup for other commands like './mfc.sh' or './mfc.sh build'

    add_common_arguments(test, "t")
    test.add_argument("-l", "--list",         action="store_true", help="List all available tests.")
    test.add_argument("-f", "--from",         default=None, type=str, help="First test UUID to run.")
    test.add_argument("-t", "--to",           default=None, type=str, help="Last test UUID to run.")
    test.add_argument("-o", "--only",         nargs="+", type=str,     default=[], metavar="L", help="Only run tests with specified properties.")
    test.add_argument("-a", "--test-all",     action="store_true",     default=False,     help="Run the Post Process Tests too.")
    test.add_argument("-%", "--percent",      type=int,                default=100,       help="Percentage of tests to run.")
    test.add_argument("-m", "--max-attempts", type=int,                default=1,         help="Maximum number of attempts to run a test.")
    test.add_argument(      "--rdma-mpi",     action="store_true",     default=False,     help="Run tests with RDMA MPI enabled")
    test.add_argument(      "--no-build",     action="store_true",     default=False,     help="(Testing) Do not rebuild MFC.")
    test.add_argument(      "--no-examples",  action="store_true",     default=False,     help="Do not test example cases." )
    test.add_argument("--case-optimization",  action="store_true",     default=False,     help="(GPU Optimization) Compile MFC targets with some case parameters hard-coded.")
    test.add_argument(      "--dry-run",      action="store_true",     default=False,     help="Build and generate case files but do not run tests.")

    test_meg = test.add_mutually_exclusive_group()
    test_meg.add_argument("--generate",          action="store_true", default=False, help="(Test Generation) Generate golden files.")
    test_meg.add_argument("--add-new-variables", action="store_true", default=False, help="(Test Generation) If new variables are found in D/ when running tests, add them to the golden files.")
    test_meg.add_argument("--remove-old-tests",  action="store_true", default=False, help="(Test Generation) Delete tests directories that are no longer.")

    # RUN
    add_common_arguments(run)
    run.add_argument("input",                      metavar="INPUT",              type=str,                     help="Input file to run.")
    run.add_argument("-e", "--engine",             choices=["interactive", "batch"],              type=str, default="interactive", help="Job execution/submission engine choice.")
    run.add_argument("-p", "--partition",          metavar="PARTITION",          type=str, default="",         help="(Batch) Partition for job submission.")
    run.add_argument("-q", "--quality_of_service", metavar="QOS",          type=str, default="",         help="(Batch) Quality of Service for job submission.")
    run.add_argument("-N", "--nodes",              metavar="NODES",              type=int, default=1,          help="(Batch) Number of nodes.")
    run.add_argument("-n", "--tasks-per-node",     metavar="TASKS",              type=int, default=1,          help="Number of tasks per node.")
    run.add_argument("-w", "--walltime",           metavar="WALLTIME",           type=str, default="01:00:00", help="(Batch) Walltime.")
    run.add_argument("-a", "--account",            metavar="ACCOUNT",            type=str, default="",         help="(Batch) Account to charge.")
    run.add_argument("-@", "--email",              metavar="EMAIL",              type=str, default="",         help="(Batch) Email for job notification.")
    run.add_argument("-#", "--name",               metavar="NAME",               type=str, default="MFC",      help="(Batch) Job name.")
    run.add_argument("-s", "--scratch",            action="store_true",                    default=False,      help="Build from scratch.")
    run.add_argument("-b", "--binary",             choices=["mpirun", "jsrun", "srun", "mpiexec"],             type=str, default=None,       help="(Interactive) Override MPI execution binary")
    run.add_argument(      "--dry-run",            action="store_true",                    default=False,      help="(Batch) Run without submitting batch file.")
    run.add_argument("--case-optimization",        action="store_true",                    default=False,      help="(GPU Optimization) Compile MFC targets with some case parameters hard-coded.")
    run.add_argument(      "--no-build",           action="store_true",                    default=False,      help="(Testing) Do not rebuild MFC.")
    run.add_argument("--wait",                     action="store_true",                    default=False,      help="(Batch) Wait for the job to finish.")
    run.add_argument("-c", "--computer",           metavar="COMPUTER",           type=str, default="default",  help=f"(Batch) Path to a custom submission file template or one of {format_list_to_string(_TEMPLATE_NAMES)}.")
    run.add_argument("-o", "--output-summary",     metavar="OUTPUT",             type=str, default=None,       help="Output file (YAML) for summary.")
    run.add_argument("--clean",                    action="store_true",                    default=False,      help="Clean the case before running.")
    run.add_argument("--ncu",                      nargs=argparse.REMAINDER,     type=str,                     help="Profile with NVIDIA Nsight Compute.")
    run.add_argument("--nsys",                     nargs=argparse.REMAINDER,     type=str,                     help="Profile with NVIDIA Nsight Systems.")
    run.add_argument("--rcu",                      nargs=argparse.REMAINDER,     type=str,                     help="Profile with ROCM rocprof-compute.")
    run.add_argument("--rsys",                     nargs=argparse.REMAINDER,     type=str,                     help="Profile with ROCM rocprof-systems.")

    # BENCH
    add_common_arguments(bench)
    bench.add_argument("-o", "--output", metavar="OUTPUT", default=None, type=str, required="True", help="Path to the YAML output file to write the results to.")
    bench.add_argument("-m", "--mem", metavar="MEM", default=1, type=int, help="Memory per task for benchmarking cases")

    # BENCH_DIFF
    add_common_arguments(bench_diff, "t")
    bench_diff.add_argument("lhs", metavar="LHS", type=str, help="Path to a benchmark result YAML file.")
    bench_diff.add_argument("rhs", metavar="RHS", type=str, help="Path to a benchmark result YAML file.")

    # COUNT
    add_common_arguments(count, "g")

    # COUNT
    add_common_arguments(count_diff, "g")

    # VALIDATE
    add_common_arguments(validate, "tjmgv")  # Only add debug-log flag
    validate.add_argument("input", metavar="INPUT", type=str, help="Path to case file to validate.")

    # INIT
    init.add_argument("name", metavar="NAME", type=str, nargs="?", default=None, help="Name/path for the new case directory.")
    init.add_argument("-t", "--template", type=str, default="1D_minimal", help="Template to use (e.g., 1D_minimal, 2D_minimal, 3D_minimal, or example:<name>).")
    init.add_argument("-l", "--list", action="store_true", help="List available templates.")

    # If enhanced help was printed, also show argparse help and exit
    if help_command:
        subparser_map = {
            "build": build, "run": run, "test": test, "clean": clean,
            "bench": bench, "count": count, "validate": validate, "init": init,
            "packer": packer, "completion": comp
        }
        if help_command in subparser_map:
            subparser_map[help_command].print_help()
        sys.exit(0)

    try:
        extra_index = sys.argv.index('--')
    except ValueError:
        extra_index = len(sys.argv)

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
    for name, subparser in [("run",    run),   ("test",   test), ("build", build),
                            ("clean",  clean), ("count", count), ("count_diff", count_diff),
                            ("validate", validate)]:
        if args["command"] == name:
            continue

        vals, _ = subparser.parse_known_args(["-i", "None"])
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

    # "Slugify" the name of the job (only for batch jobs, not for init command)
    if args.get("name") is not None and isinstance(args["name"], str) and args["command"] != "init":
        args["name"] = re.sub(r'[\W_]+', '-', args["name"])

    # We need to check for some invalid combinations of arguments because of
    # the limitations of argparse.
    if args["command"] == "build":
        if (args["input"] is not None) ^ args["case_optimization"] :
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
