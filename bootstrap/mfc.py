#!/usr/bin/env python3

import rich
import rich.console

import user
import conf
import lock
import args
import clean
import common
import signal

from run   import run
from tests import tests

class MFCState:
    def __init__(self) -> None:
        from build import MFCBuild

        rich.print(common.MFC_HEADER)

        self.conf  = conf.MFCConf(self)
        self.user  = user.MFCUser()
        self.setup_directories()
        self.lock  = lock.MFCLock(self)
        self.args  = args.parse(self)
        self.clean = clean.MFCClean(self)
        self.build = MFCBuild(self)
        self.test  = tests.MFCTest(self)
        self.run   = run.MFCRun(self)

        self.check_mode()

        rich.print(f"\n[yellow]You are currently in the [bold green]{self.lock.mode}[/bold green] mode.[/yellow]\n")

        if self.args["command"] == "test":
            self.test.test()
        elif self.args["command"] == "run":
            self.run.run()
        elif self.args["command"] == "clean":
            self.clean.run()

        if self.args["command"] == "build":
            for target_name in [ x.name for x in self.conf.targets ]:
                if target_name in self.args["targets"]:
                    rich.print("[bold][u]Build:[/u][/bold]")
                    self.build.build_target(target_name)

    def check_mode(self):
        def update_mode():
            if update_mode.triggered:
                return

            update_mode.triggered = True

            rich.print(f'[yellow]Switching to [bold green]{self.args["mode"]}[/bold green] from [bold magenta]{self.lock.mode}[/bold magenta]. Purging references to other modes...[/yellow]')

            # Update mode in mfc.user.yaml
            self.lock.mode = self.args["mode"]
            self.lock.save()

            for mode in self.user.modes:
                mode: user.Mode

                if mode.name == self.lock.mode:
                    return

                # Delete the build directory of other modes
                common.delete_directory_recursive(self.build.get_mode_base_path(mode.name))

        update_mode.triggered = False

        # User requested a new mode using -m as a command-line argument
        if self.args["mode"] != self.lock.mode:
            update_mode()

        for idx, entry in enumerate(self.lock.targets):
            entry: lock.LockTargetHolder

            # There exists a (built) target, which is not a common one, that has different mode
            if entry.target.common_mode == None and entry.metadata.mode != self.args["mode"]:
                update_mode()

                # Remove it
                del self.lock.targets[idx]

                self.lock.save()

    def setup_directories(self):
        common.create_directory(common.MFC_SUBDIR)

        for d in ["src", "build", "log", "temp"]:
            for mode in [ mode.name for mode in self.user.modes ] + ["common"]:
                common.create_directory(f"{common.MFC_SUBDIR}/{mode}/{d}")
                if d == "build":
                    for build_subdir in ["bin", "include", "lib", "share"]:
                        common.create_directory(f"{common.MFC_SUBDIR}/{mode}/{d}/{build_subdir}")

FILE_ISSUE_MSG = f"""\
We apologize for the inconvenience. If you believe this is an issue with MFC, \
please visit https://github.com/MFlowCode/MFC-develop/issues to file an issue.\
"""

if __name__ == "__main__":
    try:
        MFCState()
    except common.MFCException as exc:
        rich.print(f"""
[bold red]FATAL ERROR[/bold red]> {str(exc)}
{FILE_ISSUE_MSG}
""")
        common.quit(signal.SIGTERM)
    except KeyboardInterrupt as exc:
        common.quit(signal.SIGTERM)
    except Exception as exc:
        rich.console.Console().print_exception()
        rich.print(f"""
[bold red]FATAL ERROR[/bold red]> An unexpected exception occurred.
{FILE_ISSUE_MSG}
""")
        common.quit(signal.SIGTERM)
