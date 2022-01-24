import argparse

import internal.conf   as conf
import internal.common as common

class MFCArgs:
    def __init__(self, conf: conf.MFCConf):
        parser = argparse.ArgumentParser(description="Wecome to the MFC master script.", )

        compiler_configuration_names = [e["name"] for e in conf["configurations"]]

        parser.add_argument("--build", action="store_true", help="Build targets.")
        parser.add_argument("--test",  action="store_true", help="Test targets.")
        parser.add_argument("--clean", action="store_true", help="Clean the targets.")
        parser.add_argument("-sc", "--set-current", type=str, choices=compiler_configuration_names,
                            help="Select a compiler configuration to use when running MFC.")

        compiler_target_names = [e["name"] for e in conf["targets"]]
        parser.add_argument("-t", "--targets", nargs="+", type=str,
                            choices=compiler_target_names, default=["MFC"],
                            help="The space-separated targets you wish to have built.")

        parser.add_argument("-cc", "--compiler-configuration", type=str,
                            choices=compiler_configuration_names, default="release-cpu")

        parser.add_argument("-j", "--jobs", metavar="N", type=int,
                            help="Allows for N concurrent jobs.", default=1)

        parser.add_argument("-s", "--scratch", action="store_true",
                            help="Build all targets from scratch.")

        self.data = vars(parser.parse_args())

    def __getitem__(self, key: str, default=None):
        if key not in self.data:
            if default==None:
                raise common.MFCException(f'MFCArgs: Key "{key}" doesn\'t exist.')
            else:
                return default

        return self.data[key]
