#!/usr/bin/env python3

import rich
#import traceback

import user, conf, lock, args, test, clean, common, run

class MFCState:
    def __init__(self) -> None:
        from build import MFCBuild

        self.conf  = conf.MFCConf(self)
        self.user  = user.MFCUser()
        self.setup_directories()
        self.lock  = lock.MFCLock()
        self.args  = args.parse(self)
        self.clean = clean.MFCClean(self)
        self.build = MFCBuild(self)
        self.test  = test.MFCTest(self)
        self.run   = run.MFCRun(self)

        rich.print(common.MFC_HEADER)

        if self.args["command"] == "test":
            rich.print("[bold][u]Test:[/u][/bold]")
            self.test.test()

        if self.args["command"] == "run":
            self.run.run()

        if self.args["command"] == "clean":
            rich.print("[bold][u]Clean:[/u][/bold]")
            self.clean.run()

        for target_name in [ x.name for x in self.conf.targets ]:
            if target_name in self.args["targets"]:
                if self.args["command"] == "build":
                    rich.print("[bold][u]Build:[/u][/bold]")
                    self.build.build_target(target_name)

        self.lock.save()

    def setup_directories(self):
        common.create_directory(common.MFC_SUBDIR)

        for d in ["src", "build", "log", "temp"]:
            for mode in [ mode.name for mode in self.user.modes ] + ["common"]:
                common.create_directory(f"{common.MFC_SUBDIR}/{mode}/{d}")
                if d == "build":
                    for build_subdir in ["bin", "include", "lib", "share"]:
                        common.create_directory(f"{common.MFC_SUBDIR}/{mode}/{d}/{build_subdir}")

def main():
    mfc = MFCState()

if __name__ == "__main__":
    try:
        main()
    except common.MFCException as exc:
#        traceback.print_exc()
        rich.print(f"[red]> {str(exc)}[/red]")
        exit(1)
    except KeyboardInterrupt as exc:
        exit(1)
