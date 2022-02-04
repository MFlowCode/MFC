import argparse

import internal.user       as user
import internal.conf       as conf
import internal.objecttree as objecttree

class MFCArgs(objecttree.ObjectTree):
    def __init__(self, conf: conf.MFCConf, user: user.MFCUser):
        parser = argparse.ArgumentParser(description="Wecome to the MFC master script.", )

        compiler_configuration_names = [e.name for e in user.configurations]

        grp_func = parser.add_argument_group(title="Action")
        grp_func = grp_func.add_mutually_exclusive_group()

        grp_func.add_argument("--build", action="store_true", help="Build targets.")
        grp_func.add_argument("--test",  action="store_true", help="Test targets.")
        grp_func.add_argument("--clean", action="store_true", help="Clean MFC targets.")
        grp_func.add_argument("--set-current", type=str, choices=compiler_configuration_names,
                            help="Select a compiler configuration to use when running MFC.")

        compiler_target_names = [e.name for e in conf.targets]
        parser.add_argument("-t", "--targets", nargs="+", type=str,
                            choices=compiler_target_names, default=["MFC"],
                            help="The space-separated targets you wish to have built.")

        parser.add_argument("-cc", "--compiler-configuration", type=str,
                            choices=compiler_configuration_names, default=user.defaults.configuration)

        parser.add_argument("-j", "--jobs", metavar="N", type=int,
                            help="Allows for N concurrent jobs.", default=int(user.defaults.threads))

        parser.add_argument("-s", "--scratch", action="store_true",
                            help="Build all targets from scratch.")

        super().__init__(vars(parser.parse_args()))

        if not self.tree_get("build") and not self.tree_get("test") and not self.tree_get("clean"):
            self.tree_set("build", True)
