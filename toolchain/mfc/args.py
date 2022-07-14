import argparse


from build import get_mfc_target_names
from build import get_target_names
from build import get_dependencies_names


def parse(mfc):
    from main         import MFCState
    from run.engines  import ENGINES
    from run.mpi_bins import BINARIES

    mfc: MFCState

    parser = argparse.ArgumentParser(
        prog="./mfc.sh",
        description="Wecome to the MFC master script.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    mode_names = [ e.name for e in mfc.user.modes ]

    parsers = parser.add_subparsers(dest="command")

    run   = parsers.add_parser(name="run",   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    test  = parsers.add_parser(name="test",  formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    clean = parsers.add_parser(name="clean", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    def add_common_arguments(p, mask=""):
        if "t" not in mask:
            p.add_argument("-t", "--targets", nargs="+", type=str.lower, choices=get_target_names(),
                           default=get_mfc_target_names(), help="")

        if "m" not in mask:
            p.add_argument("-m", "--mode", type=str.lower, choices=mode_names, default=mfc.lock.mode,
                           help="Mode used to compile, run, and test MFC.")

        if "j" not in mask:
            p.add_argument("-j", "--jobs", metavar="N", type=int, default=int(mfc.user.build.threads),
                           help="Allows for N concurrent jobs.")

        for name in get_dependencies_names():
            p.add_argument(f"--no-{name}", action="store_true", help=f"Do not build the {name} dependency. Use the system version instead.")


    # === CLEAN ===
    add_common_arguments(clean, "j")

    binaries = [ b.bin for b in BINARIES ]

    # === TEST ===
    add_common_arguments(test, "t")
    test.add_argument("-g", "--generate",   action="store_true", help="Generate golden files.")
    test.add_argument("-l", "--list",       action="store_true", help="List all available tests.")
    test.add_argument("-f", "--from",       default=mfc.test.cases[0].get_uuid(), type=str, help="First test UUID to run.")
    test.add_argument("-t", "--to",         default=mfc.test.cases[-1].get_uuid(), type=str, help="Last test UUID to run.")
    test.add_argument("-o", "--only",       nargs="+", type=str, default=[], metavar="L", help="Only run tests with UUIDs or hashes L.")
    test.add_argument("-b", "--binary",     choices=binaries, type=str, default=None, help="(Serial) Override MPI execution binary")
    test.add_argument("-r", "--relentless", action="store_true", default=False, help="Run all tests, even if multiple fail.")

    # === RUN ===
    engines  = [ e.slug for e in ENGINES ]

    add_common_arguments(run)
    run.add_argument("input",                 metavar="INPUT",                 type=str,                                      help="Input file for run.")
    run.add_argument("-e", "--engine",        choices=engines,                 type=str, default=engines[0],                  help="Job execution/submission engine choice.")
    run.add_argument("-p", "--partition",     metavar="PARTITION",             type=str, default=mfc.user.run.partition,      help="(Batch) Partition for job submission.")
    run.add_argument("-N", "--nodes",         metavar="NODES",                 type=int, default=mfc.user.run.nodes,          help="(Batch) Number of nodes.")
    run.add_argument("-n", "--cpus-per-node", metavar="CPUS",                  type=int, default=mfc.user.run.cpus_per_node,  help="           Number of tasks per node.")
    run.add_argument("-g", "--gpus-per-node", metavar="GPUS",                  type=int, default=mfc.user.run.gpus_per_node,  help="(Batch) Number of GPUs  per node.")
    run.add_argument("-w", "--walltime",      metavar="WALLTIME",              type=str, default=mfc.user.run.walltime,       help="(Batch) Walltime.")
    run.add_argument("-a", "--account",       metavar="ACCOUNT",               type=str, default=mfc.user.run.account,        help="(Batch) Account to charge.")
    run.add_argument("-@", "--email",         metavar="EMAIL",                 type=str, default=mfc.user.run.email,          help="(Batch) Email for job notification.")
    run.add_argument("-#", "--name",          metavar="NAME",                  type=str, default=mfc.user.run.name,           help="(Batch) Job name.")
    run.add_argument("-f", "--flags",         metavar="FLAGS",     nargs="+",  type=str, default=mfc.user.run.flags,          help="(Batch) Additional batch options.")
    run.add_argument("-b", "--binary",        choices=binaries, type=str, default=None, help="(Interactive) Override MPI execution binary")
    run.add_argument("-s", "--scratch",       action="store_true", default=False, help="Build from scratch.")
    run.add_argument(      "--dry-run",       action="store_true", default=False, help="(Batch) Run without submitting batch file.")

    args: dict = vars(parser.parse_args())

    # Add default arguments of other subparsers
    def append_defaults_to_data(name: str, parser):
        if args["command"] != name:
            vals, errs = parser.parse_known_args(["-i None"])
            for key,val in vars(vals).items():
                if not key in args:
                    args[key] = val

    for a, b in [("run", run), ("test", test), ("clean", clean)]:
        append_defaults_to_data(a, b)

    if args["command"] is None:
        parser.print_help()
        exit(-1)

    return args
