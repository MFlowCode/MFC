import argparse
from internal.common import MFCException

import internal.user       as user
import internal.conf       as conf
import internal.objecttree as objecttree

class MFCArgs(objecttree.ObjectTree):
    def get_command(self):
        return self.tree_get("command")

    def __init__(self, conf: conf.MFCConf, user: user.MFCUser):
        parser = argparse.ArgumentParser(description="Wecome to the MFC master script.", prog="./mfc.sh",
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)

        compiler_configuration_names = [e.name for e in user.configurations]

        parsers = parser.add_subparsers(dest="command")

        build = parsers.add_parser(name="build", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        run   = parsers.add_parser(name="run",   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        test  = parsers.add_parser(name="test",  formatter_class=argparse.ArgumentDefaultsHelpFormatter)

        def add_common_arguments(p):
            compiler_target_names = [e.name for e in conf.targets]

            p.add_argument("-t", nargs="+", metavar="TARGET", type=str.lower,
                             choices=compiler_target_names, default=["mfc"], help=f"Targets [{', '.join(compiler_target_names)}]")

            p.add_argument("-c", type=str.lower, metavar="CONFIG",
                     choices=compiler_configuration_names, default=user.general.configuration, help=f"Compiler configuration ({', '.join(compiler_configuration_names)})")
            
            p.add_argument("-j", metavar="N", type=int,
                     help="Allows for N concurrent jobs.", default=int(user.build.threads))

        # === BUILD ===
        add_common_arguments(build)
        
        build.add_argument("-s", action="store_true", help="Build all targets from scratch.")

        # === TEST ===
        add_common_arguments(test)
        test.add_argument("-g", action="store_true", help="Generate golden files.")
        test.add_argument("-o", nargs="+", type=str, default=[], metavar="L", help="Only run tests with ids or hashes L.")

        # === RUN ===
        add_common_arguments(run)
        run.add_argument("-e", metavar="ENGINE", choices=["serial", "parallel"], default="serial",
                           help="Job execution/submission engine choice. [serial, parallel]")
        run.add_argument("-i", metavar="FILEPATH",  type=str, required=True,                   help="            Input file for run.")
        run.add_argument("-p", metavar="PARTITION", type=str, default=user.run.partition,      help="(Parrallel) Partition for job submission.")
        run.add_argument("-N", metavar="NODES",     type=int, default=user.run.nodes,          help="(Parrallel) Number of nodes.")
        run.add_argument("-n", metavar="TASKS",     type=int, default=user.run.tasks_per_node, help="            Number of tasks per node.")
        run.add_argument("-g", metavar="GPUS",      type=int, default=user.run.gpus_per_node,  help="(Parrallel) Number of GPUs  per node.")
        run.add_argument("-w", metavar="WALLTIME",  type=str, default=user.run.walltime,       help="(Parrallel) Walltime.")

        super().__init__(vars(parser.parse_args()))

        # Add default arguments of other subparsers
        def append_defaults_to_data(name: str, parser):
            if self.data["command"] != name:
                vals, errs = parser.parse_known_args(["-i None"])
                for key,val in vars(vals).items():
                    if not self.exists(key):
                        self.data[key] = val

        for a, b in [("run", run), ("test", test), ("build", build)]:
            append_defaults_to_data(a, b)

        if self.tree_get("command") is None:
            parser.print_help()
            exit(-1)
