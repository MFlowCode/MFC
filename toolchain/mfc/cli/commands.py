"""
MFC CLI Command Definitions - SINGLE SOURCE OF TRUTH

All command definitions live here. This file is used to generate:
- argparse parsers
- bash/zsh completions
- user guide help content
- CLI reference documentation

When adding a new command or option, ONLY modify this file.
Then run `./mfc.sh generate` to update completions.
"""
# pylint: disable=too-many-lines

from .schema import (
    CLISchema, Command, Argument, Positional, Example,
    ArgAction, Completion, CompletionType,
    CommonArgumentSet, MutuallyExclusiveGroup
)


# =============================================================================
# CONSTANTS (shared with other modules)
# =============================================================================

TARGET_NAMES = [
    'fftw', 'hdf5', 'silo', 'lapack', 'hipfort',
    'pre_process', 'simulation', 'post_process',
    'syscheck', 'documentation'
]

DEFAULT_TARGET_NAMES = ['pre_process', 'simulation', 'post_process']

TEMPLATE_NAMES = [
    'bridges2', 'carpenter', 'carpenter-cray', 'default',
    'delta', 'deltaai', 'frontier', 'hipergator', 'nautilus',
    'oscar', 'phoenix', 'phoenix-bench', 'santis', 'tuo'
]

GPU_OPTIONS = ['acc', 'mp']

ENGINE_OPTIONS = ['interactive', 'batch']

MPI_BINARIES = ['mpirun', 'jsrun', 'srun', 'mpiexec']


# =============================================================================
# COMMON ARGUMENT SETS
# =============================================================================

COMMON_TARGETS = CommonArgumentSet(
    name="targets",
    arguments=[
        Argument(
            name="targets",
            short="t",
            help="Space-separated list of targets to act upon.",
            nargs="+",
            type=str,
            default=DEFAULT_TARGET_NAMES,
            choices=TARGET_NAMES,
            metavar="TARGET",
            completion=Completion(type=CompletionType.CHOICES, choices=TARGET_NAMES),
        ),
    ]
)

COMMON_JOBS = CommonArgumentSet(
    name="jobs",
    arguments=[
        Argument(
            name="jobs",
            short="j",
            help="Allows for JOBS concurrent jobs.",
            type=int,
            default=1,
            metavar="JOBS",
        ),
    ]
)

COMMON_VERBOSE = CommonArgumentSet(
    name="verbose",
    arguments=[
        Argument(
            name="verbose",
            short="v",
            help="Increase verbosity level. Use -v, -vv, or -vvv for more detail.",
            action=ArgAction.COUNT,
            default=0,
        ),
    ]
)

COMMON_DEBUG_LOG = CommonArgumentSet(
    name="debug_log",
    arguments=[
        Argument(
            name="debug-log",
            short="d",
            help="Enable debug logging for troubleshooting.",
            action=ArgAction.STORE_TRUE,
            dest="debug_log",
        ),
    ]
)

COMMON_GPUS = CommonArgumentSet(
    name="gpus",
    arguments=[
        Argument(
            name="gpus",
            short="g",
            help="(Optional GPU override) List of GPU #s to use (environment default if unspecified).",
            nargs="+",
            type=int,
            default=None,
        ),
    ]
)

# MFCConfig flags are handled specially in argparse_gen.py
# This marker tells the generator to add --mpi/--no-mpi, --gpu/--no-gpu, etc.
COMMON_MFC_CONFIG = CommonArgumentSet(
    name="mfc_config",
    mfc_config_flags=True,
    arguments=[],  # Generated dynamically
)


# =============================================================================
# COMMAND DEFINITIONS
# =============================================================================

BUILD_COMMAND = Command(
    name="build",
    aliases=["b"],
    help="Build MFC and its dependencies.",
    description="Build MFC targets with optional GPU support and case optimization.",
    include_common=["targets", "mfc_config", "jobs", "verbose", "debug_log"],
    arguments=[
        Argument(
            name="input",
            short="i",
            help="(GPU Optimization) Build a version of MFC optimized for a case.",
            type=str,
            default=None,
            completion=Completion(type=CompletionType.FILES_PY),
        ),
        Argument(
            name="case-optimization",
            help="(GPU Optimization) Compile MFC targets with some case parameters hard-coded (requires --input).",
            action=ArgAction.STORE_TRUE,
            default=False,
            dest="case_optimization",
        ),
    ],
    examples=[
        Example("./mfc.sh build", "Build all default targets (CPU)"),
        Example("./mfc.sh build -j 8", "Build with 8 parallel jobs"),
        Example("./mfc.sh build --gpu", "Build with GPU (OpenACC) support"),
        Example("./mfc.sh build -i case.py --case-optimization -j 8", "Case optimization (10x faster!)"),
    ],
    key_options=[
        ("-j, --jobs N", "Number of parallel build jobs"),
        ("-t, --targets", "Targets: pre_process, simulation, post_process"),
        ("--gpu [acc|mp]", "Enable GPU support (OpenACC or OpenMP)"),
        ("--case-optimization", "Hard-code case params for 10x speedup"),
        ("--debug", "Build in debug mode"),
    ],
)

RUN_COMMAND = Command(
    name="run",
    aliases=["r"],
    help="Run a case with MFC.",
    description="Run an MFC simulation case interactively or submit as a batch job.",
    include_common=["targets", "mfc_config", "jobs", "verbose", "debug_log", "gpus"],
    positionals=[
        Positional(
            name="input",
            help="Input file to run.",
            completion=Completion(type=CompletionType.FILES_PY),
        ),
    ],
    arguments=[
        Argument(
            name="engine",
            short="e",
            help="Job execution/submission engine choice.",
            choices=ENGINE_OPTIONS,
            default="interactive",
            completion=Completion(type=CompletionType.CHOICES, choices=ENGINE_OPTIONS),
        ),
        Argument(
            name="partition",
            short="p",
            help="(Batch) Partition for job submission.",
            default="",
            metavar="PARTITION",
        ),
        Argument(
            name="quality_of_service",
            short="q",
            help="(Batch) Quality of Service for job submission.",
            default="",
            metavar="QOS",
        ),
        Argument(
            name="nodes",
            short="N",
            help="(Batch) Number of nodes.",
            type=int,
            default=1,
            metavar="NODES",
        ),
        Argument(
            name="tasks-per-node",
            short="n",
            help="Number of tasks per node.",
            type=int,
            default=1,
            metavar="TASKS",
            dest="tasks_per_node",
        ),
        Argument(
            name="walltime",
            short="w",
            help="(Batch) Walltime.",
            default="01:00:00",
            metavar="WALLTIME",
        ),
        Argument(
            name="account",
            short="a",
            help="(Batch) Account to charge.",
            default="",
            metavar="ACCOUNT",
        ),
        Argument(
            name="email",
            short="@",
            help="(Batch) Email for job notification.",
            default="",
            metavar="EMAIL",
        ),
        Argument(
            name="name",
            short="#",
            help="(Batch) Job name.",
            default="MFC",
            metavar="NAME",
        ),
        Argument(
            name="scratch",
            short="s",
            help="Build from scratch.",
            action=ArgAction.STORE_TRUE,
            default=False,
        ),
        Argument(
            name="binary",
            short="b",
            help="(Interactive) Override MPI execution binary",
            choices=MPI_BINARIES,
            default=None,
            completion=Completion(type=CompletionType.CHOICES, choices=MPI_BINARIES),
        ),
        Argument(
            name="dry-run",
            help="(Batch) Run without submitting batch file.",
            action=ArgAction.STORE_TRUE,
            default=False,
            dest="dry_run",
        ),
        Argument(
            name="case-optimization",
            help="(GPU Optimization) Compile MFC targets with some case parameters hard-coded.",
            action=ArgAction.STORE_TRUE,
            default=False,
            dest="case_optimization",
        ),
        Argument(
            name="no-build",
            help="(Testing) Do not rebuild MFC.",
            action=ArgAction.STORE_TRUE,
            default=False,
            dest="no_build",
        ),
        Argument(
            name="wait",
            help="(Batch) Wait for the job to finish.",
            action=ArgAction.STORE_TRUE,
            default=False,
        ),
        Argument(
            name="computer",
            short="c",
            help="(Batch) Path to a custom submission file template or one of the built-in templates.",
            default="default",
            metavar="COMPUTER",
            completion=Completion(type=CompletionType.CHOICES, choices=TEMPLATE_NAMES),
        ),
        Argument(
            name="output-summary",
            short="o",
            help="Output file (YAML) for summary.",
            default=None,
            metavar="OUTPUT",
            dest="output_summary",
        ),
        Argument(
            name="clean",
            help="Clean the case before running.",
            action=ArgAction.STORE_TRUE,
            default=False,
        ),
        # Profiler arguments with REMAINDER
        Argument(
            name="ncu",
            help="Profile with NVIDIA Nsight Compute.",
            nargs="...",  # REMAINDER
            type=str,
        ),
        Argument(
            name="nsys",
            help="Profile with NVIDIA Nsight Systems.",
            nargs="...",  # REMAINDER
            type=str,
        ),
        Argument(
            name="rcu",
            help="Profile with ROCM rocprof-compute.",
            nargs="...",  # REMAINDER
            type=str,
        ),
        Argument(
            name="rsys",
            help="Profile with ROCM rocprof-systems.",
            nargs="...",  # REMAINDER
            type=str,
        ),
    ],
    examples=[
        Example("./mfc.sh run case.py", "Run interactively with 1 rank"),
        Example("./mfc.sh run case.py -n 4", "Run with 4 MPI ranks"),
        Example("./mfc.sh run case.py --case-optimization -j 8", "10x faster with case optimization!"),
        Example("./mfc.sh run case.py -e batch -N 2 -n 4", "Submit batch job: 2 nodes, 4 ranks/node"),
    ],
    key_options=[
        ("--case-optimization", "Hard-code params for 10x speedup!"),
        ("-n, --tasks-per-node", "MPI ranks per node"),
        ("-N, --nodes", "Number of nodes (batch)"),
        ("-e, --engine", "interactive or batch"),
        ("-a, --account", "Account to charge (batch)"),
        ("-w, --walltime", "Wall time limit (batch)"),
    ],
)

TEST_COMMAND = Command(
    name="test",
    aliases=["t"],
    help="Run MFC's test suite.",
    description="Run MFC's test suite with various filtering and generation options.",
    include_common=["mfc_config", "jobs", "verbose", "debug_log", "gpus"],
    # Note: does NOT include "targets" - test uses different target handling
    arguments=[
        Argument(
            name="list",
            short="l",
            help="List all available tests.",
            action=ArgAction.STORE_TRUE,
        ),
        Argument(
            name="from",
            short="f",
            help="First test UUID to run.",
            default=None,
            type=str,
        ),
        Argument(
            name="to",
            short="t",
            help="Last test UUID to run.",
            default=None,
            type=str,
        ),
        Argument(
            name="only",
            short="o",
            help="Only run tests with specified properties.",
            nargs="+",
            type=str,
            default=[],
            metavar="L",
        ),
        Argument(
            name="test-all",
            short="a",
            help="Run the Post Process Tests too.",
            action=ArgAction.STORE_TRUE,
            default=False,
            dest="test_all",
        ),
        Argument(
            name="percent",
            short="%",
            help="Percentage of tests to run.",
            type=int,
            default=100,
        ),
        Argument(
            name="max-attempts",
            short="m",
            help="Maximum number of attempts to run a test.",
            type=int,
            default=1,
            dest="max_attempts",
        ),
        Argument(
            name="rdma-mpi",
            help="Run tests with RDMA MPI enabled",
            action=ArgAction.STORE_TRUE,
            default=False,
            dest="rdma_mpi",
        ),
        Argument(
            name="no-build",
            help="(Testing) Do not rebuild MFC.",
            action=ArgAction.STORE_TRUE,
            default=False,
            dest="no_build",
        ),
        Argument(
            name="no-examples",
            help="Do not test example cases.",
            action=ArgAction.STORE_TRUE,
            default=False,
            dest="no_examples",
        ),
        Argument(
            name="case-optimization",
            help="(GPU Optimization) Compile MFC targets with some case parameters hard-coded.",
            action=ArgAction.STORE_TRUE,
            default=False,
            dest="case_optimization",
        ),
        Argument(
            name="dry-run",
            help="Build and generate case files but do not run tests.",
            action=ArgAction.STORE_TRUE,
            default=False,
            dest="dry_run",
        ),
    ],
    mutually_exclusive=[
        MutuallyExclusiveGroup(arguments=[
            Argument(
                name="generate",
                help="(Test Generation) Generate golden files.",
                action=ArgAction.STORE_TRUE,
                default=False,
            ),
            Argument(
                name="add-new-variables",
                help="(Test Generation) If new variables are found in D/ when running tests, add them to the golden files.",
                action=ArgAction.STORE_TRUE,
                default=False,
                dest="add_new_variables",
            ),
            Argument(
                name="remove-old-tests",
                help="(Test Generation) Delete test directories that are no longer needed.",
                action=ArgAction.STORE_TRUE,
                default=False,
                dest="remove_old_tests",
            ),
        ]),
    ],
    examples=[
        Example("./mfc.sh test", "Run all tests"),
        Example("./mfc.sh test -j 4", "Run with 4 parallel jobs"),
        Example("./mfc.sh test --only 3D", "Run only 3D tests"),
        Example("./mfc.sh test --generate", "Regenerate golden files"),
    ],
    key_options=[
        ("-j, --jobs N", "Number of parallel test jobs"),
        ("-o, --only PROP", "Run tests matching property"),
        ("-f, --from UUID", "Start from specific test"),
        ("--generate", "Generate/update golden files"),
        ("--no-build", "Skip rebuilding MFC"),
    ],
)

CLEAN_COMMAND = Command(
    name="clean",
    aliases=["c"],
    help="Clean build artifacts.",
    description="Remove build artifacts and cache files.",
    include_common=["targets", "mfc_config", "jobs", "verbose", "debug_log"],
    examples=[
        Example("./mfc.sh clean", "Clean all build files"),
    ],
    key_options=[],
)

VALIDATE_COMMAND = Command(
    name="validate",
    aliases=["v"],
    help="Validate a case file without running.",
    description="Check a case file for syntax errors and constraint violations.",
    include_common=["debug_log"],  # Only debug-log from common
    positionals=[
        Positional(
            name="input",
            help="Path to case file to validate.",
            completion=Completion(type=CompletionType.FILES_PY),
        ),
    ],
    examples=[
        Example("./mfc.sh validate case.py", "Check syntax and constraints"),
        Example("./mfc.sh validate case.py -d", "Validate with debug output"),
    ],
    key_options=[
        ("-d, --debug-log", "Enable debug logging"),
    ],
)

NEW_COMMAND = Command(
    name="new",
    help="Create a new case from a template.",
    description="Create a new simulation case directory from a built-in or example template.",
    positionals=[
        Positional(
            name="name",
            help="Name/path for the new case directory.",
            nargs="?",
            default=None,
            completion=Completion(type=CompletionType.DIRECTORIES),
        ),
    ],
    arguments=[
        Argument(
            name="template",
            short="t",
            help="Template to use (e.g., 1D_minimal, 2D_minimal, 3D_minimal, or example:<name>).",
            default="1D_minimal",
            completion=Completion(
                type=CompletionType.CHOICES,
                choices=["1D_minimal", "2D_minimal", "3D_minimal"],
            ),
        ),
        Argument(
            name="list",
            short="l",
            help="List available templates.",
            action=ArgAction.STORE_TRUE,
        ),
    ],
    examples=[
        Example("./mfc.sh new my_case", "Create with 1D_minimal template"),
        Example("./mfc.sh new my_case -t 2D_minimal", "Create with 2D template"),
        Example("./mfc.sh new my_case -t example:3D_sphbubcollapse", "Copy from example"),
        Example("./mfc.sh new --list", "List available templates"),
    ],
    key_options=[
        ("-t, --template NAME", "Template: 1D_minimal, 2D_minimal, 3D_minimal"),
        ("-l, --list", "List all available templates"),
    ],
)

PACKER_COMMAND = Command(
    name="packer",
    help="Packer utility (pack/unpack/compare).",
    description="Pack simulation output into a single file or compare packed files.",
    subcommands=[
        Command(
            name="pack",
            help="Pack a case into a single file.",
            positionals=[
                Positional(
                    name="input",
                    help="Input file of case to pack.",
                    completion=Completion(type=CompletionType.FILES_PY),
                ),
            ],
            arguments=[
                Argument(
                    name="output",
                    short="o",
                    help="Base name of output file.",
                    default=None,
                    metavar="OUTPUT",
                ),
            ],
        ),
        Command(
            name="compare",
            help="Compare two cases.",
            positionals=[
                Positional(
                    name="input1",
                    help="First pack file.",
                    completion=Completion(type=CompletionType.FILES_PACK),
                ),
                Positional(
                    name="input2",
                    help="Second pack file.",
                    completion=Completion(type=CompletionType.FILES_PACK),
                ),
            ],
            arguments=[
                Argument(
                    name="reltol",
                    short="rel",
                    help="Relative tolerance.",
                    type=float,
                    default=1e-12,
                    metavar="RELTOL",
                ),
                Argument(
                    name="abstol",
                    short="abs",
                    help="Absolute tolerance.",
                    type=float,
                    default=1e-12,
                    metavar="ABSTOL",
                ),
            ],
        ),
    ],
    examples=[
        Example("./mfc.sh packer pack case.py", "Pack case output"),
        Example("./mfc.sh packer compare a.pack b.pack", "Compare two packed files"),
    ],
    key_options=[],
)

COMPLETION_COMMAND = Command(
    name="completion",
    help="Install shell tab-completion.",
    description="Install or manage shell tab-completion for MFC commands.",
    positionals=[
        Positional(
            name="completion_action",
            help="Action: install, uninstall, or status",
            nargs="?",
            default=None,
            choices=["install", "uninstall", "status"],
            completion=Completion(
                type=CompletionType.CHOICES,
                choices=["install", "uninstall", "status"],
            ),
        ),
        Positional(
            name="completion_shell",
            help="Shell type: bash or zsh (auto-detected if not specified)",
            nargs="?",
            default=None,
            choices=["bash", "zsh"],
            completion=Completion(
                type=CompletionType.CHOICES,
                choices=["bash", "zsh"],
            ),
        ),
    ],
    examples=[
        Example("./mfc.sh completion install", "Install tab completion for current shell"),
        Example("./mfc.sh completion install bash", "Install for bash specifically"),
        Example("./mfc.sh completion status", "Check completion installation status"),
    ],
    key_options=[
        ("install", "Install completion scripts"),
        ("uninstall", "Remove completion scripts"),
        ("status", "Check installation status"),
    ],
)

HELP_COMMAND = Command(
    name="help",
    help="Show help on a topic.",
    positionals=[
        Positional(
            name="topic",
            help="Help topic: gpu, clusters, batch, debugging, performance",
            nargs="?",
            default=None,
            choices=["gpu", "clusters", "batch", "debugging", "performance"],
            completion=Completion(
                type=CompletionType.CHOICES,
                choices=["gpu", "clusters", "batch", "debugging", "performance"],
            ),
        ),
    ],
)

# Simple commands (shell scripts, minimal arguments)
LOAD_COMMAND = Command(
    name="load",
    help="Loads the MFC environment with source.",
    description="Load MFC environment modules (use with source).",
    examples=[
        Example("source ./mfc.sh load -c p -m g", "Load Phoenix GPU modules"),
        Example("source ./mfc.sh load -c f -m c", "Load Frontier CPU modules"),
    ],
    key_options=[
        ("-c CLUSTER", "Cluster: p(hoenix), f(rontier), a(ndes), etc."),
        ("-m MODE", "Mode: c(pu), g(pu)"),
    ],
)

LINT_COMMAND = Command(
    name="lint",
    help="Lints and tests all toolchain code.",
    description="Run pylint and unit tests on MFC's toolchain Python code.",
    arguments=[
        Argument(
            name="no-test",
            help="Skip running unit tests (only run pylint).",
            action=ArgAction.STORE_TRUE,
        ),
    ],
    examples=[
        Example("./mfc.sh lint", "Run pylint and unit tests"),
        Example("./mfc.sh lint --no-test", "Run only pylint (skip unit tests)"),
    ],
    key_options=[
        ("--no-test", "Skip unit tests, only run pylint"),
    ],
)

FORMAT_COMMAND = Command(
    name="format",
    help="Formats all code after editing.",
    description="Format all code in the repository.",
    examples=[
        Example("./mfc.sh format", "Format all code"),
        Example("./mfc.sh format -j 8", "Format with 8 parallel jobs"),
    ],
    key_options=[
        ("-j, --jobs N", "Number of parallel formatting jobs"),
    ],
)

SPELLING_COMMAND = Command(
    name="spelling",
    help="Runs the spell checker after editing.",
    description="Check spelling in the codebase.",
    examples=[
        Example("./mfc.sh spelling", "Run spell checker"),
    ],
)

INTERACTIVE_COMMAND = Command(
    name="interactive",
    help="Launch interactive menu-driven interface.",
    description="Launch an interactive menu for MFC operations.",
)

GENERATE_COMMAND = Command(
    name="generate",
    help="Regenerate completion scripts from CLI schema.",
    description="Regenerate shell completion scripts, documentation, and JSON schema from the CLI schema.",
    arguments=[
        Argument(
            name="check",
            help="Check if generated files are up to date (exit 1 if not).",
            action=ArgAction.STORE_TRUE,
            default=False,
        ),
        Argument(
            name="json-schema",
            help="Generate JSON Schema for IDE auto-completion of case files.",
            action=ArgAction.STORE_TRUE,
            default=False,
            dest="json_schema",
        ),
    ],
    examples=[
        Example("./mfc.sh generate", "Regenerate completion scripts"),
        Example("./mfc.sh generate --check", "Check if completions are up to date"),
        Example("./mfc.sh generate --json-schema", "Generate JSON Schema for IDE support"),
    ],
)

# Bench commands (CI-specific)
BENCH_COMMAND = Command(
    name="bench",
    help="Benchmark MFC (for CI).",
    include_common=["targets", "mfc_config", "jobs", "verbose", "debug_log", "gpus"],
    arguments=[
        Argument(
            name="output",
            short="o",
            help="Path to the YAML output file to write the results to.",
            required=True,
            metavar="OUTPUT",
        ),
        Argument(
            name="mem",
            short="m",
            help="Memory per task for benchmarking cases",
            type=int,
            default=1,
            metavar="MEM",
        ),
    ],
    examples=[
        Example("./mfc.sh bench -o results.yaml", "Run benchmarks and save results"),
    ],
    key_options=[
        ("-o, --output FILE", "Output file for benchmark results (required)"),
        ("-m, --mem SIZE", "Memory limit for benchmarks"),
    ],
)

BENCH_DIFF_COMMAND = Command(
    name="bench_diff",
    help="Compare MFC Benchmarks (for CI).",
    include_common=["mfc_config", "jobs", "verbose", "debug_log"],
    # Note: does NOT include "targets"
    positionals=[
        Positional(name="lhs", help="Path to a benchmark result YAML file."),
        Positional(name="rhs", help="Path to a benchmark result YAML file."),
    ],
)

COUNT_COMMAND = Command(
    name="count",
    help="Count LOC in MFC.",
    include_common=["targets", "mfc_config", "jobs", "verbose", "debug_log"],
    examples=[
        Example("./mfc.sh count", "Show LOC statistics"),
    ],
)

COUNT_DIFF_COMMAND = Command(
    name="count_diff",
    help="Compare LOC between branches.",
    include_common=["targets", "mfc_config", "jobs", "verbose", "debug_log"],
)

PARAMS_COMMAND = Command(
    name="params",
    help="Search and explore MFC case parameters.",
    description="Search, list, and get information about MFC's ~3,300 case parameters.",
    positionals=[
        Positional(
            name="query",
            help="Search query (parameter name or pattern to search for).",
            nargs="?",
            default=None,
        ),
    ],
    arguments=[
        Argument(
            name="type",
            short="t",
            help="Filter by type: int, real, log, str.",
            choices=["int", "real", "log", "str"],
            default=None,
            dest="param_type",
        ),
        Argument(
            name="families",
            short="f",
            help="List parameter families (grouped by prefix).",
            action=ArgAction.STORE_TRUE,
            default=False,
        ),
        Argument(
            name="features",
            short="F",
            help="List feature groups (mhd, bubbles, weno, etc.).",
            action=ArgAction.STORE_TRUE,
            default=False,
        ),
        Argument(
            name="feature",
            help="Show all parameters for a feature group (e.g., --feature mhd).",
            type=str,
            default=None,
            metavar="NAME",
        ),
        Argument(
            name="names-only",
            help="Only search parameter names (not descriptions).",
            action=ArgAction.STORE_TRUE,
            default=False,
            dest="names_only",
        ),
        Argument(
            name="count",
            short="c",
            help="Show count statistics only.",
            action=ArgAction.STORE_TRUE,
            default=False,
        ),
        Argument(
            name="limit",
            short="n",
            help="Maximum number of results to show (default: unlimited).",
            type=int,
            default=10000,
        ),
        Argument(
            name="describe",
            short="d",
            help="Show parameter descriptions.",
            action=ArgAction.STORE_TRUE,
            default=False,
        ),
    ],
    examples=[
        Example("./mfc.sh params weno", "Search for parameters (names + descriptions)"),
        Example("./mfc.sh params magnetic", "Find params mentioning 'magnetic'"),
        Example("./mfc.sh params --feature mhd", "Show all MHD-related parameters"),
        Example("./mfc.sh params -F", "List all feature groups"),
        Example("./mfc.sh params -f", "List parameter families (by prefix)"),
        Example("./mfc.sh params -c", "Show parameter count statistics"),
    ],
    key_options=[
        ("--feature NAME", "Show params for a feature (mhd, bubbles, weno, etc.)"),
        ("-F, --features", "List all feature groups"),
        ("-f, --families", "List parameter families (by prefix)"),
        ("-t, --type", "Filter by type (int, real, log, str)"),
        ("--names-only", "Only search names (skip descriptions)"),
    ],
)


# =============================================================================
# HELP TOPICS
# =============================================================================

HELP_TOPICS = {
    "gpu": {
        "title": "GPU Configuration",
        "description": "How to configure GPU builds and runs",
    },
    "clusters": {
        "title": "Cluster Configuration",
        "description": "How to configure MFC for different HPC clusters",
    },
    "batch": {
        "title": "Batch Job Submission",
        "description": "How to submit batch jobs with MFC",
    },
    "debugging": {
        "title": "Debugging & Troubleshooting",
        "description": "Tips for debugging MFC issues",
    },
}


# =============================================================================
# COMPLETE CLI SCHEMA
# =============================================================================

MFC_CLI_SCHEMA = CLISchema(
    prog="./mfc.sh",
    description="""\
Welcome to the MFC master script. This tool automates and manages building, testing, \
running, and cleaning of MFC in various configurations on all supported platforms. \
The README documents this tool and its various commands in more detail. To get \
started, run `./mfc.sh build -h`.""",

    arguments=[
        Argument(
            name="help",
            short="h",
            help="Show help message",
            action=ArgAction.STORE_TRUE,
        ),
    ],

    commands=[
        BUILD_COMMAND,
        RUN_COMMAND,
        TEST_COMMAND,
        CLEAN_COMMAND,
        VALIDATE_COMMAND,
        NEW_COMMAND,
        PARAMS_COMMAND,
        PACKER_COMMAND,
        COMPLETION_COMMAND,
        HELP_COMMAND,
        GENERATE_COMMAND,
        LOAD_COMMAND,
        LINT_COMMAND,
        FORMAT_COMMAND,
        SPELLING_COMMAND,
        INTERACTIVE_COMMAND,
        BENCH_COMMAND,
        BENCH_DIFF_COMMAND,
        COUNT_COMMAND,
        COUNT_DIFF_COMMAND,
    ],

    common_sets=[
        COMMON_TARGETS,
        COMMON_JOBS,
        COMMON_VERBOSE,
        COMMON_DEBUG_LOG,
        COMMON_GPUS,
        COMMON_MFC_CONFIG,
    ],

    help_topics=HELP_TOPICS,
)


# =============================================================================
# DERIVED DATA (for use by other modules)
# =============================================================================

# Command aliases mapping (replaces COMMAND_ALIASES in user_guide.py)
COMMAND_ALIASES = {}
for cmd in MFC_CLI_SCHEMA.commands:
    for alias in cmd.aliases:
        COMMAND_ALIASES[alias] = cmd.name

# Commands dict (replaces COMMANDS in user_guide.py)
def get_commands_dict():
    """Generate COMMANDS dict from schema for user_guide.py compatibility."""
    return {
        cmd.name: {
            "description": cmd.description or cmd.help,
            "alias": cmd.aliases[0] if cmd.aliases else None,
            "examples": [(e.command, e.description) for e in cmd.examples],
            "key_options": list(cmd.key_options),
        }
        for cmd in MFC_CLI_SCHEMA.commands
    }

COMMANDS = get_commands_dict()
