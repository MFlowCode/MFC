import re, sys, os.path, argparse, dataclasses

from .run.run      import get_baked_templates
from .build        import TARGETS, DEFAULT_TARGETS
from .common       import MFCException, format_list_to_string
from .test.cases   import list_cases
from .state        import gpuConfigOptions, MFCConfig

# pylint: disable=too-many-locals, too-many-branches, too-many-statements
def parse(config: MFCConfig):
    parser = argparse.ArgumentParser(
        prog="./mfc.sh",
        description="""\
Welcome to the MFC master script. This tool automates and manages building, testing, \
running, and cleaning of MFC in various configurations on all supported platforms. \
The README documents this tool and its various commands in more detail. To get \
started, run ./mfc.sh build -h.""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

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
            p.add_argument("-t", "--targets", metavar="TARGET", nargs="+", type=str.lower, choices=[ _.name for _ in TARGETS ],
                           default=[ _.name for _ in sorted(DEFAULT_TARGETS, key=lambda t: t.runOrder) ],
                           help=f"Space separated list of targets to act upon. Allowed values are: {format_list_to_string([ _.name for _ in TARGETS ])}.")

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

        if "g" not in mask:
            p.add_argument("-g", "--gpus", nargs="+", type=int, default=None, help="(Optional GPU override) List of GPU #s to use (environment default if unspecified).")

    # BUILD
    add_common_arguments(build, "g")
    build.add_argument("-i", "--input", type=str, default=None, help="(GPU Optimization) Build a version of MFC optimized for a case.")
    build.add_argument("--case-optimization", action="store_true", default=False, help="(GPU Optimization) Compile MFC targets with some case parameters hard-coded (requires --input).")

    # TEST
    test_cases = list_cases()

    add_common_arguments(test, "t")
    test.add_argument("-l", "--list",         action="store_true", help="List all available tests.")
    test.add_argument("-f", "--from",         default=test_cases[0].get_uuid(), type=str, help="First test UUID to run.")
    test.add_argument("-t", "--to",           default=test_cases[-1].get_uuid(), type=str, help="Last test UUID to run.")
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
    run.add_argument("-c", "--computer",           metavar="COMPUTER",           type=str, default="default",  help=f"(Batch) Path to a custom submission file template or one of {format_list_to_string(list(get_baked_templates().keys()))}.")
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

    try:
        extra_index = sys.argv.index('--')
    except ValueError:
        extra_index = len(sys.argv)

    args: dict = vars(parser.parse_args(sys.argv[1:extra_index]))
    args["--"] = sys.argv[extra_index + 1:]

    # Add default arguments of other subparsers
    for name, parser in [("run",    run),   ("test",   test), ("build", build),
                         ("clean",  clean), ("count", count), ("count_diff", count_diff)]:
        if args["command"] == name:
            continue

        vals, _ = parser.parse_known_args(["-i", "None"])
        for key, val in vars(vals).items():
            if key == "input":
                args[key] = args.get(key)
            elif not key in args:
                args[key] = args.get(key, val)

    if args["command"] is None:
        parser.print_help()
        exit(-1)

    # "Slugify" the name of the job
    args["name"] = re.sub(r'[\W_]+', '-', args["name"])

    # We need to check for some invalid combinations of arguments because of
    # the limitations of argparse.
    if args["command"] == "build":
        if (args["input"] is not None) ^ args["case_optimization"] :
            raise MFCException("./mfc.sh build's --case-optimization and --input must be used together.")
    if args["command"] == "run":
        if args["binary"] is not None and args["engine"] != "interactive":
            raise MFCException("./mfc.sh run's --binary can only be used with --engine=interactive.")

    # Input files to absolute paths
    for e in ["input", "input1", "input2"]:
        if e not in args:
            continue

        if args.get(e) is not None:
            args[e] = os.path.abspath(args[e])

    return args
