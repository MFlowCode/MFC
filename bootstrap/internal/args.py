import argparse
from internal.common import MFCException

import internal.user       as user
import internal.conf       as conf
import internal.objecttree as objecttree

class MFCArgs(objecttree.ObjectTree):
    def __init__(self, conf: conf.MFCConf, user: user.MFCUser):
        parser = argparse.ArgumentParser(description="Wecome to the MFC master script.", )

        compiler_configuration_names = [e.name for e in user.configurations]

        action  = parser.add_argument_group(title="Action")
        general = parser.add_argument_group(title="General")
        build   = parser.add_argument_group(title="Build")
        test    = parser.add_argument_group(title="Test")

        test.add_argument("-g", "--generate", action="store_true", help="Generate golden files.")

        grp = action.add_mutually_exclusive_group(required=True)

        grp.add_argument("--build", action="store_true", help="Build targets.")
        grp.add_argument("--test",  action="store_true", help="Test targets.")
        grp.add_argument("--clean", action="store_true", help="Clean MFC targets.")
        grp.add_argument("--set-current", type=str.lower, choices=compiler_configuration_names,
                            help="Select a compiler configuration to use when running MFC.")

        compiler_target_names = [e.name for e in conf.targets]
        action.add_argument("-t", "--targets", nargs="+", type=str.lower,
                            choices=compiler_target_names, default=["mfc"],
                            help="The space-separated targets you wish to have built.")

        general.add_argument("-cc", "--compiler-configuration", type=str.lower,
                            choices=compiler_configuration_names, default=user.defaults.configuration)

        build.add_argument("-j", "--jobs", metavar="N", type=int,
                            help="Allows for N concurrent jobs.", default=int(user.defaults.threads))

        build.add_argument("-s", "--scratch", action="store_true",
                            help="Build all targets from scratch.")

        super().__init__(vars(parser.parse_args()))

        if not self.tree_get("build") and not self.tree_get("test") and not self.tree_get("clean") and not self.tree_get("set_current"):
            parser.print_help()
            exit(-1)

