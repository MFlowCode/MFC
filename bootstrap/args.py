import argparse


def parse(mfc):
    from mfc import MFCState

    mfc: MFCState

    parser = argparse.ArgumentParser(description="Wecome to the MFC master script.", prog="./mfc.sh",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    compiler_configuration_names = [e.name for e in mfc.user.configurations]

    parsers = parser.add_subparsers(dest="command")

    build = parsers.add_parser(name="build", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    run   = parsers.add_parser(name="run",   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    test  = parsers.add_parser(name="test",  formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    def add_common_arguments(p):
        compiler_target_names = [e.name for e in mfc.conf.targets]

        p.add_argument("-t", "--targets", nargs="+", type=str.lower,
                            choices=compiler_target_names, default=["mfc"], help="")

        p.add_argument("-cc", "--compiler-configuration", type=str.lower,
                    choices=compiler_configuration_names, default=mfc.user.general.configuration, help="")
        
        p.add_argument("-j", "--jobs", metavar="N", type=int,
                    help="Allows for N concurrent jobs.", default=int(mfc.user.build.threads))

    # === BUILD ===
    add_common_arguments(build)
    
    build.add_argument("-s", "--scratch", action="store_true",
                        help="Build all targets from scratch.")

    # === TEST ===
    add_common_arguments(test)
    test.add_argument("-g", "--generate", action="store_true", help="Generate golden files.")
    test.add_argument("-o", "--only",     nargs="+", type=str, default=[], metavar="L", help="Only run tests with ids or hashes L.")

    # === RUN ===
    add_common_arguments(run)
    run.add_argument("-e", "--engine", choices=["serial", "parallel"], default="serial",
                        help="Job execution/submission engine choice.")
    run.add_argument("-i", "--input",                               type=str, required=True,                       help="Input file for run.")
    run.add_argument("-p", "--partition",      metavar="PARTITION", type=str, default=mfc.user.run.partition,      help="(Parallel) Partition for job submission.")
    run.add_argument("-N", "--nodes",          metavar="NODES",     type=int, default=mfc.user.run.nodes,          help="(Parallel) Number of nodes.")
    run.add_argument("-n", "--tasks-per-node", metavar="TASKS",     type=int, default=mfc.user.run.tasks_per_node, help="           Number of tasks per node.")
    run.add_argument("-g", "--gpus-per-node",  metavar="GPUS",      type=int, default=mfc.user.run.gpus_per_node,  help="(Parallel) Number of GPUs  per node.")
    run.add_argument("-w", "--walltime",       metavar="WALLTIME",  type=str, default=mfc.user.run.walltime,       help="(Parallel) Walltime.")
    run.add_argument("-a", "--account",        metavar="ACCOUNT",   type=str, default=mfc.user.run.account,        help="(Parallel) Account to charge.")

    args: dict = vars(parser.parse_args())

    # Add default arguments of other subparsers
    def append_defaults_to_data(name: str, parser):
        if args["command"] != name:
            vals, errs = parser.parse_known_args(["-i None"])
            for key,val in vars(vals).items():
                if not key in args:
                    args[key] = val

    for a, b in [("run", run), ("test", test), ("build", build)]:
        append_defaults_to_data(a, b)

    if args["command"] is None:
        parser.print_help()
        exit(-1)

    return args
